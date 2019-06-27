/*
 *  This file is part of the MOLECFIT software package.
 *  Copyright (C) 2009-2013 European Southern Observatory
 *
 *  This programme is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This programme is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this programme. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \ingroup molecfit
 */

/**@{*/

/*!
 * \file mf_molecfit.c
 *
 * Top-level routines for MOLECFIT
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  05 Nov 2010
 * \date   03 Nov 2014
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <mf_molecfit.h>

#if __GNUC__ >= 4 || __clang__
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
#endif


/*****************************************************************************
 *                                GLOBALS                                    *
 ****************************************************************************/

/* Definition of global variables */

/* Reference fit parameters */
//cpl_array *reffitpar;
/* Number of MPFIT fitting function calls */
int nfev = 0;
/* Number of LBLRTM calls (wavenumber-restricted subspectra are not
   counted individually) */
int n_code = 0;
/* Total LBLRTM run time */
double t_code = 0.;

extern cpl_array *reffitpar;


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

static cpl_error_code
mf_molecfit_(mfdrv * drvpar, const char * parfile, const char mode,
             cpl_boolean read_kernel, cpl_table** prof_out, cpl_table** res_out, cpl_table** spec_out);

#define MF_OVERRIDE_PAR(NAME, SET_PAR, VALUE) \
    do { \
        cpl_parameter * par__ = cpl_parameterlist_find(drvpar->parlist, NAME); \
        SET_PAR(par__, VALUE); \
        cpl_error_ensure(cpl_parameterlist_find(parlist, NAME) == NULL, \
                         CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_ILLEGAL_INPUT, \
                         NAME " must not be set in input parameters"); \
    } while (0)

/* opaque state structure for multiple calctrans runs
 * current used calctrans code will check if reuse is appropriate already
 * (=wavelength ranges equal to last run otherwise recompute),
 * otherwise one could store relevant parameters here and check they are the
 * same in the next call */
struct DLL_PUBLIC mf_calctrans_state_ {
    char * tmpdir;
};

//void DLL_PUBLIC test() {}

mf_calctrans_state DLL_PUBLIC
*mf_init_calctrans_state(void)
{
    /*!
     * \callgraph
     *
     * Initialize calctrans sate, must be deleted with mf_cleanup
     *
     */
    return cpl_calloc(1, sizeof(mf_calctrans_state));
}

static const char * mf_get_basedir(void) {
    /*!
     * \callgraph
     *
     * Gets molecfit installation directory containing binaries and data.
     * Uses Env variable MOLECFITDIR or hardcoded installation prefix.
     *
     * \b OUTPUT:
     * \param basedir base directory
     *
     */
    const char * basedir = getenv("MOLECFITDIR");
    if (basedir) {
        return basedir;
    } else {
        return MOLECFIT_BASEDIR;
    }
}

const char * mf_get_datadir(void) {
  /*!
   * \callgraph
   *
   * Gets molecfit installation directory containing binaries and data.
   * Uses Env variable MOLECFITDIR_DATA or hardcoded installation prefix.
   *
   * \b OUTPUT:
   * \param datadir data directory
   *
   */
    const char * datadir = getenv("MOLECFITDIR_DATA");
    if (datadir) {
        return datadir;
    } else {
        return MOLECFIT_DATADIR;
    }
}

cpl_error_code DLL_PUBLIC
mf_cleanup(mf_calctrans_state * state)
{
    /*!
     * \callgraph
     *
     * Cleanup calctrans state and remove temporary data
     *
     * \b INPUT:
     * \param state state structure allocated by mf_init_calctrans_state
     *
     */
    if (state == NULL) {
        return CPL_ERROR_NONE;
    }

    if (state->tmpdir) {
        char * cmd = cpl_sprintf("rm -rf \"%s\"", state->tmpdir);
        if (system(cmd) != 0) {
            cpl_msg_error(cpl_func, "Removing temporary directory %s failed.",
                          state->tmpdir);
        }
        cpl_free(cmd);
        cpl_free(state->tmpdir);
    }
    cpl_free(state);

    // Clean global variables
    if (reffitpar) {
        cpl_array_delete(reffitpar);
        reffitpar = NULL;
    }
    nfev = 0;
    n_code = 0;
    t_code = 0.;

    return cpl_error_get_code();
}

cpl_error_code
mf_cleanup_standalone(char* tmpdir)
{
    /*!
     * \callgraph
     *
     * Remove temporary data
     *
     * \b INPUT:
     * \param tmpdir string containing the path of the tmp directory
     *
     */

    if (tmpdir) {
        char * cmd = cpl_sprintf("rm -rf \"%s\"", tmpdir);
        if (system(cmd) != 0) {
            cpl_msg_error(cpl_func, "Removing temporary directory %s failed.",
                          tmpdir);
        }
        cpl_free(cmd);
        cpl_free(tmpdir);
    }

    return cpl_error_get_code();
}

cpl_error_code
mf_override_config(mfdrv * drvpar,
                   cpl_parameterlist * parlist, cpl_propertylist * plist,
                   cpl_table ** inspec, cpl_size nspec,
                   cpl_table * molectab,
                   cpl_table * wlinclude,
                   cpl_table * wlexclude,
                   cpl_table * pixexclude,
                   cpl_matrix * kernel,
                   mf_calctrans_state * state
                  )
{
    /*!
     * \callgraph
     *
     * Override default initialized configuration of ::mfdrv based on pipeline
     * provided data
     *
     * \b INPUT:
     * \param drvpar default initialized ::mfdrv structure
     * \param parlist parameter list of configuration values to apply to drvpar
     * \param plist property list containing values molecfit reads (ESO TEL/INS)
     * \param inspec array of input spectra, each entry is considered a chip
     * \param nspec number input spectra
     * \param molectab table defining molecules to fit
     * \param wlinclude wavelength inclusion table
     * \param wlexclude wavelength exclusion table
     * \param pixexclude pixel exclusion table
     * \param kernel kernel, one row per total input spectra pixel
     * \param state computation state variable
     *
     */
    cpl_parameter * par;
    char * tmpdir = NULL;

    cpl_ensure_code(nspec > 0, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(inspec, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(drvpar, CPL_ERROR_NULL_INPUT);

    if (state && state->tmpdir) {
        tmpdir = cpl_strdup(state->tmpdir);
    }
    else {
        if (getenv("TMPDIR")) {
            tmpdir = cpl_sprintf("%s/molecfit_tmp_XXXXXX", getenv("TMPDIR"));
        }
        else {
            tmpdir = cpl_strdup("/tmp/molecfit_tmp_XXXXXX");
        }

        if (!mf_get_tempdir(tmpdir)) {
            cpl_free(tmpdir);
            return cpl_error_get_code();
        }
    }

    /* TODO check first and lastcall influence in mf_modsim */

    /* setup output directories and make sure we do not write into the base
     * installation */

    MF_OVERRIDE_PAR("basedir", cpl_parameter_set_string, mf_get_basedir());

    /* written to (gdas and lnfl) so it needs to be copied */
    char * datadir = cpl_sprintf("%s/data/", tmpdir);
    MF_OVERRIDE_PAR("datadir", cpl_parameter_set_string, datadir);
    cpl_free(datadir);

    char * gdasdir = cpl_sprintf("%s/data/profiles/grib", tmpdir);
    MF_OVERRIDE_PAR("gdas_dir", cpl_parameter_set_string, gdasdir);
    cpl_free(gdasdir);

    /* molecfit writes into the data directory
     * gdas updates and reusable LNFL files
     * Copy the data directory as the installation one may be read-only
     * TODO the large hitran and aer folders could be symlinked */
    if (state == NULL || (state && state->tmpdir == NULL)) {
        cpl_msg_info(cpl_func, "Creating workspace in %s", tmpdir);
        char * cmd = cpl_sprintf("cp -r %s %s", mf_get_datadir(), tmpdir);
        int status = (system(cmd));
        cpl_free(cmd);
        if (status != 0) {
            cpl_error_code code = cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                                        "Could not copy %s/data to %s",
                                                        mf_get_basedir(), tmpdir);
            cpl_free(tmpdir);
            return code;
        }
    }
    else {
        cpl_msg_info(cpl_func, "Reusing workspace in %s", tmpdir);
    }

    MF_OVERRIDE_PAR("output_dir", cpl_parameter_set_string, tmpdir);

    MF_OVERRIDE_PAR("output_name", cpl_parameter_set_string, "output_prefix");

   /* store inclusions and exclusions as fits files to be used loaded later
     * when code reads the input spectra
     * TODO split spectra and range loading so input can be used directly */

    /* TODO ref_atm support or dissallow changing it for now */
    if (wlinclude) {
        char * fn_wlinc = cpl_sprintf("%s/wrange_include.fits", tmpdir);
        cpl_table_save(wlinclude, NULL, NULL, fn_wlinc, CPL_IO_CREATE);
        MF_OVERRIDE_PAR("wrange_include", cpl_parameter_set_string, fn_wlinc);
        cpl_free(fn_wlinc);
    }
    if (wlexclude) {
        char * fn_wlexc = cpl_sprintf("%s/wrange_exclude.fits", tmpdir);
        cpl_table_save(wlexclude, NULL, NULL, fn_wlexc, CPL_IO_CREATE);
        MF_OVERRIDE_PAR("wrange_exclude", cpl_parameter_set_string, fn_wlexc);
        cpl_free(fn_wlexc);
    }
    if (pixexclude) {
        char * fn_pixexc = cpl_sprintf("%s/prange_exclude.fits", tmpdir);
        cpl_table_save(pixexclude, NULL, NULL, fn_pixexc, CPL_IO_CREATE);
        MF_OVERRIDE_PAR("prange_exclude", cpl_parameter_set_string, fn_pixexc);
        cpl_free(fn_pixexc);
    }

    /* store input spectra to be loaded later again, also store header as by
     * default many parameters are defined by ESO header keys */
    char * fn_spec = cpl_sprintf("%s/inputspec.fits", tmpdir);
    cpl_table_save(inspec[0], plist, NULL, fn_spec, CPL_IO_CREATE);
    for (cpl_size i = 1; i < nspec; i++) {
        cpl_table_save(inspec[i], NULL, NULL, fn_spec, CPL_IO_EXTEND);
    }

    MF_OVERRIDE_PAR("filename", cpl_parameter_set_string, fn_spec);
    cpl_free(fn_spec);

    /* update the parameters that are tables frmo the input */
    /* TODO might be mandatory */
    if (molectab) {
        cpl_table_delete(drvpar->molectab);
        drvpar->molectab = cpl_table_duplicate(molectab);
        mf_par_finalize_lbl_molecs(drvpar);
        cpl_error_ensure(cpl_parameterlist_find(parlist, "list_molec") == NULL &&
                         cpl_parameterlist_find(parlist, "fit_molec") == NULL &&
                         cpl_parameterlist_find(parlist, "relcol") == NULL,
                         CPL_ERROR_ILLEGAL_INPUT, return CPL_ERROR_ILLEGAL_INPUT,
                         "molecule configuration must be provided via "
                         "molectab table");
    }

    /* TODO
     * possibly handle rangetab and chiptab, they are a mix of input parameters
     * from expert mode and values filled after reading the spectra
     * for now expert mode should be useable by adding the correct parameter
     * list keys, but if needed a better API could be made
     */

    /* normalize and use user provided kernel*/
    if (kernel) {
        cpl_size npix = 0;
        for (cpl_size i = 0; i < nspec; i++) {
            npix += cpl_table_get_nrow(inspec[i]);
        }
        /* broadcasting 1 row kernel could be added as convenience */
        if (npix != cpl_matrix_get_nrow(kernel)) {
            cpl_free(tmpdir);
            return cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                                         "Input kernel must have as may rows "
                                         "as total input spectrum has "
                                         "pixels (%ld)", (long)npix);
        }

        cpl_matrix_delete(drvpar->kernel);
        drvpar->kernel = cpl_matrix_duplicate(kernel);
        /* normalize kernel */
        for (cpl_size r = 0; r < cpl_matrix_get_nrow(drvpar->kernel); r++) {
            double sum = 0.;
            for (cpl_size c = 0; c < cpl_matrix_get_ncol(drvpar->kernel); c++) {
                sum += cpl_matrix_get(drvpar->kernel, r, c);
            }
            for (cpl_size c = 0; c < cpl_matrix_get_ncol(drvpar->kernel); c++) {
                cpl_matrix_set(drvpar->kernel, r, c,
                               cpl_matrix_get(drvpar->kernel, r, c) / sum);
            }
        }

        /* needs to be set to something it is actually used in mf_modsim */
        MF_OVERRIDE_PAR("kernel_file", cpl_parameter_set_string, "notnone");
    }

    /* transfer all other parameters from input to the default initialized */
    for (par = cpl_parameterlist_get_first(parlist); par != NULL;
         par = cpl_parameterlist_get_next(parlist)) {
        const char * pname = cpl_parameter_get_name(par);
        cpl_parameter * lpar = cpl_parameterlist_find(drvpar->parlist, pname);
        if (lpar == NULL) {
            cpl_msg_warning(cpl_func, "Unknown parameter %s", pname);
            continue;
        }
        cpl_type tpar = cpl_parameter_get_type(par);
        if (tpar != cpl_parameter_get_type(lpar)) {
            cpl_msg_warning(cpl_func, "Parameter type mismatch %s", pname);
            continue;
        }
        switch (tpar) {
            case CPL_TYPE_BOOL:
                cpl_parameter_set_bool(lpar, cpl_parameter_get_bool(par));
                break;
            case CPL_TYPE_INT:
                cpl_parameter_set_int(lpar, cpl_parameter_get_int(par));
                break;
            case CPL_TYPE_DOUBLE:
                cpl_parameter_set_double(lpar, cpl_parameter_get_double(par));
                break;
            case CPL_TYPE_STRING:
                cpl_parameter_set_string(lpar, cpl_parameter_get_string(par));
                break;
            default:
                cpl_msg_warning(cpl_func, "Unsupported parameter type %s", pname);
                continue;
        }
    }

    if (state) {
        // Make sure to delete memory before we point to somewhere else.
        if (state->tmpdir) {
            cpl_free(state->tmpdir);
        }
        state->tmpdir = cpl_strdup(tmpdir);
    }

    cpl_free(tmpdir);
    return CPL_ERROR_NONE;
}

cpl_error_code DLL_PUBLIC mf_run_molecfit(
    cpl_parameterlist  *parlist,
    cpl_propertylist   *plist,
    cpl_table          **inspec,
    cpl_size           nspec,
    cpl_table          *molectab,
    cpl_table          *wlinclude,
    cpl_table          *wlexclude,
    cpl_table          *pixexclude,
    cpl_matrix         *kernel,
    cpl_table          **prof_out,
    cpl_table          **res_out,
    cpl_table          **spec_out,
    mf_calctrans_state **state)
{
    /*!
     * \callgraph
     *
     * Run molecfit in temporary folder
     *
     * \b INPUT:
     * \param drvpar     default initialized ::mfdrv structure
     * \param parlist    parameter list of configuration values to apply to drvpar
     * \param plist      property list containing values molecfit reads (ESO TEL/INS)
     * \param inspec     array of input spectra, each entry is considered a chip
     * \param nspec      number input spectra
     * \param molectab   table defining molecules to fit
     * \param wlinclude  wavelength inclusion table
     * \param wlexclude  wavelength exclusion table
     * \param pixexclude pixel exclusion table
     * \param kernel     kernel, one row per total input spectra pixel
     * \param prof_out   output atmospheric profile of mf_run_molecfit(...)
     * \param res_out    output results table of mf_run_molecfit(...)
     * \param spec_out   output spec table
     * \param state      computation state variable
     *
     */

    mfdrv drvpar;
    cpl_ensure_code(state, CPL_ERROR_NULL_INPUT);

    mf_par_initall(&drvpar);

    /* molecfit changes cwd during run, save and restore */
    char * orig_cwd = mf_get_cwd();

    mf_override_config(&drvpar,
                       parlist, plist, inspec, nspec,
                       molectab,
                       wlinclude,
                       wlexclude,
                       pixexclude,
                       kernel,
                       *state);

    if (kernel) {
        // Debug
        cpl_msg_info(cpl_func, "Using user-defined Convolution Kernel");
        //cpl_matrix_dump(kernel, NULL);
    }

    cpl_error_code code = mf_molecfit_(&drvpar, "", 'm', kernel == NULL,
        prof_out, res_out, spec_out);

    if (chdir(orig_cwd)) {}
    cpl_free(orig_cwd);

    return code;
}


static cpl_error_code
mf_molecfit_(mfdrv * drvpar, const char * parfile, const char mode,
             cpl_boolean read_kernel, cpl_table** prof_out, cpl_table** res_out,
             cpl_table** spec_out)
{
    /*!
     * \callgraph
     *
     * See mf_molecfit
     *
     */
    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_table *prof;
    mp_result result;
    int plotrange = 0, nrange = 0, i = 0;
    double ts = 0., te = 0., fittime = 0.;

    /* Remove file types which depend on the number of fit ranges */
    mf_molecfit_clean(drvpar);

    /* Read spectral and header data from data file */
    *spec_out = cpl_table_new(0);
    if ((status = mf_readspec(*spec_out, drvpar)) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(*spec_out);
        return status;
    }

    /* Read fixed kernel if provided */
    if (read_kernel &&
        (status = mf_par_readkernel(drvpar)) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(*spec_out);
        return status;
    }

    /* Create atmospheric profiles from reference profiles and GDAS data */
    prof = cpl_table_new(0);
    if ((status = mf_atm_createatm(prof, drvpar)) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(*spec_out);
        cpl_table_delete(prof);
        return status;
    }

    /* Adapt atmospheric profiles to ESO meteo monitor data */
    mf_atm_adaptenv(prof, drvpar);

    /* Scale H2O profile to given PWV value if requested */
    mf_atm_scaletopwv(prof, drvpar);

    /* Run LNFL if required */
    mf_lnfl(drvpar);

    /* Exit programme if errors occurred */
    if ((status = cpl_error_get_code()) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(*spec_out);
        cpl_table_delete(prof);
        return status;
    }

    /* Perform fitting procedure */
    ts = cpl_test_get_walltime();
    if (mode == 'm') {
        /* Use fitting plan */
        status = mf_molecfit_batch(&result, *spec_out, prof, prof_out, res_out, drvpar);
    } else if (mode == 's') {
        /* Single call of CMPFIT */
        status = mf_mpfit(&result, *spec_out, prof, prof_out, res_out, drvpar);
    }

    /* Get CMPFIT run time in min */
    te = cpl_test_get_walltime();
    fittime = (te - ts) / 60;

    /* Print fit results */
    cpl_msg_info(cpl_func, "FIT RESULTS:");
    if (status != CPL_ERROR_NONE && result.status >= 0) {
        result.status = -99;
    }
    cpl_msg_info(cpl_func, "status: %d", result.status);
    if (status == CPL_ERROR_NONE) {
        cpl_msg_info(cpl_func, "npar: %d", result.npar);
        cpl_msg_info(cpl_func, "npix: %d", result.nfunc);
        cpl_msg_info(cpl_func, "niter: %d", result.niter);
        cpl_msg_info(cpl_func, "nfev: %d", result.nfev);
        cpl_msg_info(cpl_func, "fittime: %.2f min", fittime);
        cpl_msg_info(cpl_func, "t_code/n_code: %.2lf s",
                     t_code / ((double)n_code));
        cpl_msg_info(cpl_func, "n_code: %i", n_code);
        cpl_msg_info(cpl_func, "orignorm: %.3e", result.orignorm);
        p = cpl_parameterlist_find(drvpar->parlist, "chi2");
        cpl_msg_info(cpl_func, "bestnorm: %.3e",
                     cpl_parameter_get_double(p));
    }

    /* Write a summary of the fit results into an ASCII file in output
       folder */
    mf_mpfit_writeresults(drvpar, *spec_out, &result, fittime, res_out);

    /* Calculate function for telluric absorption correction if a transmission
       spectrum was fitted */
    mf_molecfit_calctacfunc(*spec_out, drvpar);

    /* Set model flux to 0 for non-fitted pixels if desired */
    mf_molecfit_cleanmodelflux(*spec_out, drvpar);

    /* Write content of CPL table "spec" to ASCII and FITS file in output
       folder */
    mf_molecfit_writefile(*spec_out, drvpar);

    /* Write best-fit parameters to driver file in output folder */
    mf_par_writefile(drvpar, parfile);

    /* Remove file types which depend on the number of fit ranges */
    mf_molecfit_clean(drvpar);

    /* Plot observed and best-fit spectrum (full range and individual fit
       ranges if desired */
    if (status == CPL_ERROR_NONE) {
        mf_molecfit_plot(*spec_out, drvpar, 0);
        p = cpl_parameterlist_find(drvpar->parlist, "plot_range");
        plotrange = cpl_parameter_get_int(p);
        if (plotrange == 1) {
            p = cpl_parameterlist_find(drvpar->parlist, "nrange");
            nrange = cpl_parameter_get_int(p);
            for (i = 0; i < nrange; i++) {
                mf_molecfit_plot(*spec_out, drvpar, i+1);
            }
        }
    }

    /* Free allocated memory */
    cpl_table_delete(prof);
    mf_par_deleteall(drvpar);
    mf_mpfit_freememresult(&result);

    /* Return error code of last error */
    return cpl_error_get_code();
}

cpl_error_code fix_directories(mfdrv* drvpar, char** tmpdir) {

    MF_OVERRIDE_PAR_STANDALONE("basedir", cpl_parameter_set_string, mf_get_basedir());

	// Get a path for the temporary directory
    if (getenv("TMPDIR")) {
        *tmpdir = cpl_sprintf("%s/molecfit_tmp_XXXXXX", getenv("TMPDIR"));
    }
    else {
        *tmpdir = cpl_strdup("/tmp/molecfit_tmp_XXXXXX");
    }

    if (!mf_get_tempdir(*tmpdir)) {
        cpl_free(*tmpdir);
        return cpl_error_get_code();
    }

	// Create the directory
	cpl_msg_info(cpl_func, "Creating workspace in %s", *tmpdir);
    char * cmd = cpl_sprintf("cp -r %s %s", mf_get_datadir(), *tmpdir);
    int status = (system(cmd));
    cpl_free(cmd);
    if (status != 0) {
        cpl_error_code code = cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                                    "Could not copy %s/data to %s",
                                                    mf_get_basedir(), *tmpdir);
        cpl_free(*tmpdir);
        return code;
    }

    char * datadir = cpl_sprintf("%s/data/", *tmpdir);
    MF_OVERRIDE_PAR_STANDALONE("datadir", cpl_parameter_set_string, datadir);
    cpl_free(datadir);

    char * gdasdir = cpl_sprintf("%s/data/profiles/grib", *tmpdir);
    MF_OVERRIDE_PAR_STANDALONE("gdas_dir", cpl_parameter_set_string, gdasdir);
    cpl_free(gdasdir);

    return CPL_ERROR_NONE;
}

cpl_error_code mf_molecfit(const char *parfile_in, const char mode)
{
    /*!
     * \callgraph
     *
     * This is the top-level routine of MOLECFIT that runs the CMPFIT-based
     * fitting procedure for transmission or thermal-IR radiance spectra.
     * It needs an input driver parameter file and requires a CMPFIT run mode.
     * Either a single run or a sequence of runs with different but fixed
     * fit flag settings can be selected. As result of the fitting procedure
     * a series of output files is written: an updated driver parameter file,
     * an ASCII file which summarises the fit results (including the ppmv of
     * the different molecules and especially the water-related PWV in mm), an
     * ASCII file with the best-fit atmospheric profiles, and an ASCII and a
     * FITS file with the observed and best-fit model spectrum (plus weights
     * and deviations). All files have the same name but different suffixes.
     *
     * \b INPUT:
     * \param parfile_in  name of parameter file
     * \param mode        single run ('s') or multiple runs ('m') of CMPFIT
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    mfdrv drvpar;
    char errtxt[MF_MAXLEN], parfile[MF_MAXLEN];

    /* Check mode */
    if (mode != 's' && mode != 'm') {
        sprintf(errtxt, "%s: mode ('s' or 'm' only)", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Read MOLECFIT driver file */
    mf_par_initall(&drvpar);

    mf_basic_absfile(parfile, parfile_in);
    if ((status = mf_par_readfile(&drvpar, parfile, 0)) != CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        return status;
    }

    // Modification to stop using system directories in case of a
    // systenwide installation
    char* tmpdir = NULL;
    if ((status = fix_directories(&drvpar, &tmpdir)) != CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        return status;
    }

    /* Make output directory if required */
    mf_molecfit_makeoutputdir(&drvpar);

    /* Copy MOLECFIT driver file to output directory and rename it */
    if ((status = mf_par_copyfile(&drvpar, parfile)) != CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        return status;
    }

    cpl_table* prof_out = NULL;
    cpl_table* res_out = NULL;
    cpl_table* spec_out = NULL;
    status = mf_molecfit_(&drvpar, parfile, mode, CPL_TRUE, &prof_out, &res_out, &spec_out);
    if (prof_out) cpl_table_delete(prof_out);
    if (res_out) cpl_table_delete(res_out);
    if (spec_out) cpl_table_delete(spec_out);

    mf_cleanup_standalone(tmpdir);

    return status; 
}

cpl_error_code mf_molecfit_makeoutputdir(const mfdrv *drvpar)
{
    /*!
     * Makes output directory as given by the ::mfdrv parameter structure if
     * required.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * -
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    char outdir[MF_MAXLEN] = "";
    char path[MF_MAXLEN] = "";

    /* Get output path */
    // If the path is relative, use the current directory as a root.
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    strcpy(outdir, cpl_parameter_get_string(p));
    // mf_basic_abspath(path, outdir, basedir);
    mf_basic_abspath(path, outdir, curdir);

    /* Create output directory if not present */
    if (access(path, W_OK) != 0) {
        cpl_msg_info(cpl_func, "Make output directory %s", path);
        if (mkdir(path, 0777)) {};
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_molecfit_clean(const mfdrv *drvpar)
{
    /*!
     * Removes file types which depend on the number of fit ranges. In fact,
     * this concerns the output spectra of the selected radiative transfer
     * code and plot files for the different fit ranges. Only those files
     * files are deleted that match the name space given in the ::mfdrv
     * parameter structure.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * -
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    char outdir[MF_MAXLEN], spectype = 'T';
    char outspec[MF_MAXLEN], sys[MF_MAXLEN], outplot[MF_MAXLEN];
    int trans = 0;

    /* Get output folder and name space */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);

    /* Get names of spectra (radiative transfer code) and remove files if
       present */
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);
    if (trans == 0) {
        spectype = 'R';
    } else {
        spectype = 'T';
    }
    // p = cpl_parameterlist_find(drvpar->parlist, "spec_out");
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    sprintf(outspec, "%s%s_*_%c.fits", outdir, cpl_parameter_get_string(p),
            spectype);
    sprintf(sys, "rm -rf %s", outspec);
    if (system(sys)) {};

    /* Get names of plots and remove files if present */
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    sprintf(outplot, "%s%s_fit*.ps", outdir, cpl_parameter_get_string(p));
    sprintf(sys, "rm -rf %s", outplot);
    if (system(sys)) {};

    return CPL_ERROR_NONE;
}


cpl_error_code mf_molecfit_batch(mp_result *result, cpl_table *spec,
                                 cpl_table *prof, cpl_table** prof_out,
                                 cpl_table** res_out,
                                 mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * This routine runs the fitting procedure several times. For each fit
     * the fit flags are changed.
     *
     * \b INPUT:
     * \param spec     CPL table with observed spectrum
     * \param prof     CPL table with atmospheric profiles
     * \param drvpar   ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec     CPL table with observed and best-fit model spectrum
     * \param result   CMPFIT structure for fit results
     *
     * \b ERRORS:
     * - No data
     * - Insufficient memory
     * - see ::mf_mpfit
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_array *fit_molec_0, *fit_molec_1, *fit_cont_0, *fit_cont_1;
    cpl_array *fit_wlc_0, *fit_wlc_1, *fit_res_0, *fit_res_1;
    char errtxt[MF_MAXLEN];
    int nmolec = 0, nrange = 0, trans = 0, ncont = 0, nchip = 0, nwlc = 0;
    int nres = 0, i = 0, runcheck[5] = {1, 1, 1, 1, 1}, niter = 0;
    double orignorm = 0.;

    /* Get size of fit flag arrays */

    nmolec = cpl_table_get_nrow(drvpar->molectab);
    if (nmolec < 1) {
        sprintf(errtxt, "%s: cpl_table drvpar->molectab", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }
    nrange = cpl_table_get_nrow(drvpar->rangetab);
    if (nrange < 1) {
        sprintf(errtxt, "%s: cpl_table drvpar->rangetab", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);
    if (trans == 0) {
        ncont = nrange+1;
    } else {
        ncont = nrange;
    }
    nchip = cpl_table_get_nrow(drvpar->chiptab);
    if (nchip < 1) {
        sprintf(errtxt, "%s: cpl_table drvpar->chiptab", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }
    nwlc = nchip;
    nres = 3;

    /* Initialise fit flag arrays */

    fit_molec_0 = cpl_array_new(nmolec, CPL_TYPE_INT);
    cpl_array_fill_window_int(fit_molec_0, 0, nmolec, 0);
    fit_molec_1 = cpl_array_new(nmolec, CPL_TYPE_INT);
    cpl_array_fill_window_int(fit_molec_1, 0, nmolec, 0);
    fit_cont_0 = cpl_array_new(ncont, CPL_TYPE_INT);
    cpl_array_fill_window_int(fit_cont_0, 0, ncont, 0);
    fit_cont_1 = cpl_array_new(ncont, CPL_TYPE_INT);
    cpl_array_fill_window_int(fit_cont_1, 0, ncont, 0);
    fit_wlc_0 = cpl_array_new(nwlc, CPL_TYPE_INT);
    cpl_array_fill_window_int(fit_wlc_0, 0, nwlc, 0);
    fit_wlc_1 = cpl_array_new(nwlc, CPL_TYPE_INT);
    cpl_array_fill_window_int(fit_wlc_1, 0, nwlc, 0);
    fit_res_0 = cpl_array_new(nres, CPL_TYPE_INT);
    cpl_array_fill_window_int(fit_res_0, 0, nres, 0);
    fit_res_1 = cpl_array_new(nres, CPL_TYPE_INT);
    cpl_array_fill_window_int(fit_res_1, 0, nres, 0);

    /* Fill non-zero fit flag arrays */

    /* Set fit_molec_1 fit flags */
    for (i = 0; i < nmolec; i++) {
        cpl_array_set_int(fit_molec_1, i,
                          cpl_table_get(drvpar->molectab, "fit_molec", i,
                                        NULL));
    }

    /* Set fit_cont_1 fit flags */
    for (i = 0; i < nrange; i++) {
        cpl_array_set_int(fit_cont_1, i,
                          cpl_table_get(drvpar->rangetab, "fit_range", i,
                                        NULL));
    }
    if (trans == 0) {
        /* Fit telescope background if relevant */
        p = cpl_parameterlist_find(drvpar->parlist, "fit_back");
        cpl_array_set_int(fit_cont_1, nrange, cpl_parameter_get_int(p));
    }

    /* Set fit_wlc_1 fit flag */
    for (i = 0; i < nchip; i++) {
        cpl_array_set_int(fit_wlc_1, i,
                          cpl_table_get(drvpar->chiptab, "fit_chip", i,
                                        NULL));
    }

    /* Set fit_res_1 fit flags */
    p = cpl_parameterlist_find(drvpar->parlist, "fit_res_box");
    cpl_array_set_int(fit_res_1, 0, cpl_parameter_get_int(p));
    p = cpl_parameterlist_find(drvpar->parlist, "fit_res_gauss");
    cpl_array_set_int(fit_res_1, 1, cpl_parameter_get_int(p));
    p = cpl_parameterlist_find(drvpar->parlist, "fit_res_lorentz");
    cpl_array_set_int(fit_res_1, 2, cpl_parameter_get_int(p));

    /* Find "non-zero" fit flag arrays with all elements = 0 */
    if (cpl_array_get_max(fit_molec_1) == 0) {
        runcheck[0] = 0;
    }
    if (cpl_array_get_max(fit_cont_1) == 0) {
        runcheck[1] = 0;
    }
    if (cpl_array_get_max(fit_wlc_1) == 0) {
        runcheck[2] = 0;
    }
    if (cpl_array_get_max(fit_res_1) == 0) {
        runcheck[3] = 0;
    }

    /* Run CMPFIT */

    /* Deactivate flag for last modsim call */
    p = cpl_parameterlist_find(drvpar->parlist, "lastcall");
    cpl_parameter_set_int(p, -1);

    /* 1st run: get the continuum first approximation */
    cpl_msg_info(cpl_func, "1st run: get the continuum first "
                 "approximation");
    if (runcheck[1] == 0) {
        cpl_msg_info(cpl_func, "-> skipped!");
    } else {
        mf_par_setfitflags(drvpar, fit_molec_0, fit_cont_1, fit_wlc_0,
                           fit_res_0);
        status = mf_mpfit(result, spec, prof, prof_out, res_out, drvpar);
        niter += result->niter;
        orignorm = result->orignorm;
        mf_mpfit_freememresult(result);
    }

    /* 2nd run: get the wavelength calibration and resolution first
       approximation */
    cpl_msg_info(cpl_func, "2nd run: get the wavelength calibration and "
                 "resolution first approximation");
    if ((runcheck[2] == 0 && runcheck[3] == 0) || status != CPL_ERROR_NONE) {
        cpl_msg_info(cpl_func, "-> skipped!");
    } else {
        mf_par_setfitflags(drvpar, fit_molec_0, fit_cont_0, fit_wlc_1,
                           fit_res_1);
        status = mf_mpfit(result, spec, prof, prof_out, res_out, drvpar);
        niter += result->niter;
        if (runcheck[1] == 0) {
            orignorm = result->orignorm;
        }
        mf_mpfit_freememresult(result);
    }

    /* 3rd run: get the continuum second approximation */
    cpl_msg_info(cpl_func, "3rd run: get the continuum second "
                 "approximation");
    if (runcheck[1] == 0) {
        cpl_msg_info(cpl_func, "-> skipped!");
    } else {
        mf_par_setfitflags(drvpar, fit_molec_0, fit_cont_1, fit_wlc_0,
                           fit_res_0);
        status = mf_mpfit(result, spec, prof, prof_out, res_out, drvpar);
        niter += result->niter;
        if (runcheck[1] == 0) {
            orignorm = result->orignorm;
        }
        mf_mpfit_freememresult(result);
    }

    /* Activate flag for last modsim call */
    p = cpl_parameterlist_find(drvpar->parlist, "lastcall");
    cpl_parameter_set_int(p, 0);

    /* 4th run: get the column densities first approximation */
    cpl_msg_info(cpl_func, "4th run: get the column densities first "
                 "approximation");
    if (runcheck[0] == 0 || status != CPL_ERROR_NONE) {
        cpl_msg_info(cpl_func, "-> skipped!");
    } else {
        mf_par_setfitflags(drvpar, fit_molec_1, fit_cont_0, fit_wlc_0,
                           fit_res_0);
        status = mf_mpfit(result, spec, prof, prof_out, res_out, drvpar);
        niter += result->niter;
        if (runcheck[1] == 0 && runcheck[2] == 0 && runcheck[3] == 0) {
            orignorm = result->orignorm;
        }
        mf_mpfit_freememresult(result);
    }

    /* Deactivate flag for last modsim call */
    p = cpl_parameterlist_find(drvpar->parlist, "lastcall");
    cpl_parameter_set_int(p, -1);

    /* 5th run: redo continuum, wavelength, resolution */
    cpl_msg_info(cpl_func, "5th run: redo continuum, wavelength, "
                 "resolution");
    if ((runcheck[1] == 0 && runcheck[2] == 0 && runcheck[3] == 0) ||
        status != CPL_ERROR_NONE) {
        cpl_msg_info(cpl_func, "-> skipped!");
    } else {
        mf_par_setfitflags(drvpar, fit_molec_0, fit_cont_1, fit_wlc_1,
                           fit_res_1);
        status = mf_mpfit(result, spec, prof, prof_out, res_out, drvpar);
        niter += result->niter;
        if (runcheck[0] == 1) {
            mf_mpfit_freememresult(result);
        }
    }

    /* Activate flag for last modsim call */
    p = cpl_parameterlist_find(drvpar->parlist, "lastcall");
    cpl_parameter_set_int(p, 0);

    /* 6th run: final adjustment */
    cpl_msg_info(cpl_func, "6th run: final adjustment");
    if (runcheck[0] == 0 || status != CPL_ERROR_NONE) {
        cpl_msg_info(cpl_func, "-> skipped!");
    } else {
        mf_par_setfitflags(drvpar, fit_molec_1, fit_cont_1, fit_wlc_1,
                           fit_res_1);
        status = mf_mpfit(result, spec, prof, prof_out, res_out, drvpar);
        niter += result->niter;
    }

    /* Allocate memory for results structure if all fit flags = 0 */
    if (runcheck[0] == 0 && runcheck[1] == 0 && runcheck[2] == 0 &&
        runcheck[3] == 0) {
        mf_mpfit_allocmemresult(result, 1, 1);
        result->status = -19;
        status = -1;
    }

    /* Correct run-dependent parameters of results structure */
    result->niter = niter;
    result->nfev = nfev;
    result->orignorm = orignorm;

    /* Free allocated memory */
    cpl_array_delete(fit_molec_0);
    cpl_array_delete(fit_molec_1);
    cpl_array_delete(fit_cont_0);
    cpl_array_delete(fit_cont_1);
    cpl_array_delete(fit_wlc_0);
    cpl_array_delete(fit_wlc_1);
    cpl_array_delete(fit_res_0);
    cpl_array_delete(fit_res_1);

    return status;
}


cpl_error_code mf_molecfit_calctacfunc(cpl_table *spec, const mfdrv *drvpar)
{
    /*!
     * Adds the column "mtrans", which includes a transmission curve that can
     * be used for telluric absorption correction. It differs from the
     * best-fit model by the neglection of the continuum fit, i.e. the
     * unabsorbed continuum is characterised by a value of 1. No column is
     * created if the fitted spectrum is an telluric emission spectrum.
     *
     * \b INPUT:
     * \param spec    CPL table with observed and modelled spectrum
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    table with additional column for model transmission
     *                curve (only for telluric features in absorption)
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    int trans = 0;

    /* Return if plain sky spectrum was fitted */
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);
    if (trans == 0) {
        return CPL_ERROR_NONE;
    }

    /* Create new column and set values to model flux */
    cpl_table_duplicate_column(spec, "mtrans", spec, "mflux");

    /* Divide model spectrum by continuum scaling function */
    cpl_table_divide_columns(spec, "mtrans", "mscal");

    return CPL_ERROR_NONE;
}


cpl_error_code mf_molecfit_cleanmodelflux(cpl_table *spec,
                                          const mfdrv *drvpar)
{
    /*!
     * Sets the model flux to zero for pixels that were excluded from the
     * fitting procedure. Plots created afterwards will then only show data
     * for the fitted pixels. To switch on this option, the parameter
     * \e clean_mflux has to be set to 1.
     *
     * \b INPUT:
     * \param spec    CPL table with observed and modelled spectrum
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    table with cleaned mflux column if requested
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    int nrow = 0, i = 0;
    double *mflux, *weight;

    /* Return if cleaning is not requested */
    p = cpl_parameterlist_find(drvpar->parlist, "clean_mflux");
    if (cpl_parameter_get_int(p) == 0) {
        return CPL_ERROR_NONE;
    }

    /* Get pointers to CPL table columns mflux and weight */
    mflux = cpl_table_get_data_double(spec, "mflux");
    weight = cpl_table_get_data_double(spec, "weight");

    /* Get number of pixels */
    nrow = cpl_table_get_nrow(spec);

    /* Set mflux to 0 for pixels with weight = 0 */
    for (i = 0; i < nrow; i++) {
        if (weight[i] == 0.) {
            mflux[i] = 0.;
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_molecfit_writefile(cpl_table *spec, const mfdrv *drvpar)
{
    /*!
     * Writes CPL table with observed (column "flux") and best-fit spectrum
     * (column "mflux") into an ASCII and a FITS file. Further output columns
     * are chip ("chip"), wavelength ("lambda"), weight of observed data
     * ("weight"), model weight ("mweight"), and weighted difference between
     * model and observed spectrum ("dev"). For telluric features in
     * absorption, the model transmission curve ("mtrans") is also written
     * out. Writing the ASCII file can be switched off by setting the
     * parameter \e writeascii to 0.
     *
     * \b INPUT:
     * \param spec     CPL table with observed and modelled spectrum
     * \param drvpar   ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * -
     *
     * \b ERRORS:
     * - none
     */

    FILE *stream;
    cpl_parameter *p;
    char outdir[MF_MAXLEN], outname[MF_MAXLEN];
    char outdat[MF_MAXLEN] = "", outfits[MF_MAXLEN];
    char curdir[MF_MAXLEN];

    /* Get output folder and name space */
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);

    /* Write info message */
    cpl_msg_info(cpl_func, "Write fit results into output folder %s", outdir);

    /* Write ASCII file if desired */
    p = cpl_parameterlist_find(drvpar->parlist, "writeascii");
    if (cpl_parameter_get_int(p) == 1) {
        sprintf(outdat, "%s%s_fit.asc", outdir, outname);
        stream = fopen(outdat, "w+");
        // printf("DEBUG outdat=%s %p\n", outdat, (void*) stream);
        cpl_table_dump(spec, 0, cpl_table_get_nrow(spec), stream);
        fclose(stream);
    }

    /* Write FITS file */
    sprintf(outfits, "%s%s_fit.fits", outdir, outname);
    cpl_table_save(spec, NULL, NULL, outfits, CPL_IO_CREATE);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_molecfit_plot(cpl_table *spec, const mfdrv *drvpar,
                                const int range)
{
    /*!
     * Compares a spectrum corrected for telluric absorption with the
     * uncorrected spectrum. GNUPLOT is used for plotting.
     *
     * \b INPUT:
     * \param spec    CPL table with observed and corrected spectrum
     * \param drvpar  ::mfdrv parameter structure
     * \param range   range number
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    FILE *specfile1, *specfile2, *specfile3, *gnufile;
    cpl_parameter *p;
    cpl_table *rangespec;
    char errtxt[MF_MAXLEN], tmpdir[MF_MAXLEN], tmpfile1[MF_MAXLEN];
    char tmpfile2[MF_MAXLEN], tmpfile3[MF_MAXLEN];
    char plottype[MF_MAXLEN], outdir[MF_MAXLEN];
    char outname[MF_MAXLEN], psfile[MF_MAXLEN], tmpfile4[MF_MAXLEN];
    char systemcall[MF_MAXLEN];
    int nrow = 0, i = 0, dummy = 0, j = 0, plotopt = 2;
    double xmin = 0., xmax = 0., ymin = 0., ymax = 0., dy = 0.;
    double lam = 0., flux1 = 0., weight1 = 0., flux2 = 0., weight2 = 0.;
    double dflux = 0.;

    /* Labels of required table columns */
    char collam[] = "lambda";
    char colflux1[] = "flux";
    char colweight1[] = "weight";
    char colflux2[] = "mflux";
    char colweight2[] = "mweight";

    /* Extension of y-axis in per cent */
    double del = 0.05;

    /* Temporary filenames */
    char filename1[] = "obsspec.dat";
    char filename2[] = "modspec.dat";
    char filename3[] = "diffspec.dat";
    char gnuname[] = "plot.gnu";

    /* Plot labels */
    char xlabel[] = "Wavelength [micron]";
    char ylabel[] = "Radiance";
    char title[] = "Comparison of observed and best-fit model spectrum";

    /* Check existence of required columns */
    if (cpl_table_has_column(spec, collam) != 1 ||
        cpl_table_has_column(spec, colflux1) != 1 ||
        cpl_table_has_column(spec, colweight1) != 1 ||
        cpl_table_has_column(spec, colflux2) != 1 ||
        cpl_table_has_column(spec, colweight2) != 1) {
        sprintf(errtxt, "%s: cpl_table *spec (required columns not found)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Skip empty ranges */
    if (range != 0 &&
        cpl_table_get(drvpar->rangetab, "wn_end", range-1, NULL) == 0.) {
        return CPL_ERROR_NONE;
    }

    /* Extract plot range from input spectrum */
    if (range == 0) {
        nrow = cpl_table_get_nrow(spec);
        rangespec = cpl_table_duplicate(spec);
    } else {
        cpl_table_unselect_all(spec);
        nrow = cpl_table_or_selected_int(spec, "mrange", CPL_EQUAL_TO, range);
        rangespec = cpl_table_extract_selected(spec);
        cpl_table_select_all(spec);
    }

    /* Exit if plot range is zero */
    if (nrow == 0) {
        return CPL_ERROR_NONE;
    }

    /* Get wavelength range for x-axis */
    xmin = cpl_table_get(rangespec, collam, 0, NULL);
    xmax = cpl_table_get(rangespec, collam, nrow-1, NULL);

    /* Get flux range for y-axis */
    for (ymin = 1000., ymax = 0., i = 0; i < nrow; i++) {
        flux1 = cpl_table_get(rangespec, colflux1, i, NULL);
        weight1 = cpl_table_get(rangespec, colweight1, i, NULL);
        if (flux1 < ymin && weight1 > 0.) {
            ymin = flux1;
        }
        if (flux1 > ymax && weight1 > 0.) {
            ymax = flux1;
        }
    }

    /* Handle missing model spectrum */
    if (ymax < ymin) {
        ymin = cpl_table_get_column_min(rangespec, colflux1);
        ymax = cpl_table_get_column_max(rangespec, colflux1);
    }

    /* Extend y-axis range */
    dy = ymax - ymin;
    ymin -= del * dy;
    ymax += del * dy;

    /* Create temporary directory */
    sprintf(tmpdir, "__tmpDIRtmp__");
    if (access(tmpdir, W_OK) == 0) {
        cpl_msg_warning(cpl_func, "Directory %s already exists!", tmpdir);
    } else {
        if ((dummy = mkdir(tmpdir, 0777))) {};
    }

    /* Write ASCII files containing observed, fitted, and difference
       spectra */
    sprintf(tmpfile1, "%s/%s", tmpdir, filename1);
    specfile1 = fopen(tmpfile1, "w");
    sprintf(tmpfile2, "%s/%s", tmpdir, filename2);
    specfile2 = fopen(tmpfile2, "w");
    sprintf(tmpfile3, "%s/%s", tmpdir, filename3);
    specfile3 = fopen(tmpfile3, "w");
    for (i = 0; i < nrow; i++) {
        lam = cpl_table_get(rangespec, collam, i, NULL);
        flux1 = cpl_table_get(rangespec, colflux1, i, NULL);
        weight1 = cpl_table_get(rangespec, colweight1, i, NULL);
        flux2 = cpl_table_get(rangespec, colflux2, i, NULL);
        weight2 = cpl_table_get(rangespec, colweight2, i, NULL);
        dflux = flux1 - flux2;
        if (weight1 == 0. || weight2 == 0.) {
            dflux = 0.;
        } else {
            dflux = flux1 - flux2;
        }
        fprintf(specfile1, "%5.6g\t%5.6g\n", lam, flux1);
        fprintf(specfile2, "%5.6g\t%5.6g\n", lam, flux2);
        fprintf(specfile3, "%5.6g\t%5.6g\n", lam, dflux);
    }
    fclose(specfile1);
    fclose(specfile2);
    fclose(specfile3);

    /* Delete temporary CPL table */
    cpl_table_delete(rangespec);

    /* Check plot options */
    p = cpl_parameterlist_find(drvpar->parlist, "plot_creation");
    sprintf(plottype, "%s", cpl_parameter_get_string(p));

    /* Get path and name of POSTSCRIPT file if required */
    if (strchr(plottype, 'P') != NULL) {
        char curdir[MF_MAXLEN];
        p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
        mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
        p = cpl_parameterlist_find(drvpar->parlist, "output_name");
        strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);
        if (range == 0) {
            sprintf(psfile, "%s%s_fit.ps", outdir, outname);
        } else {
            sprintf(psfile, "%s%s_fit_%d.ps", outdir, outname, range);
        }
    }

    /* Create GNUPLOT driver files and run them */

    for (j = 0; j < 3; j++) {

        /* Plot type requested? */
        if (j == 0 && strchr(plottype, 'W') != NULL) {
            plotopt = 1;
        } else if (j == 1 && strchr(plottype, 'X') != NULL) {
            plotopt = 2;
        } else if (j == 2 && strchr(plottype, 'P') != NULL) {
            plotopt = 3;
        } else {
            continue;
        }

        /* Open temporary GNUPLOT file */
        mf_basic_initstring(tmpfile4, MF_MAXLEN);
        sprintf(tmpfile4, "%s/%s", tmpdir, gnuname);
        gnufile = fopen(tmpfile4, "w");

        /* Write lines dependent on plot type */
        if (plotopt == 1) {
            /* Plot on WXT terminal */
            fprintf(gnufile, "set term wxt\n");
            fprintf(gnufile, "set termoption enhanced\n");
        } else if (plotopt == 2) {
            /* Plot on X11 terminal */
            fprintf(gnufile, "set term x11\n");
            fprintf(gnufile, "set termoption enhanced\n");
        } else if (plotopt == 3) {
            /* POSTSCRIPT file */
            fprintf(gnufile, "set term postscript enhanced color\n");
            fprintf(gnufile, "set output \"%s\"\n", psfile);
            fprintf(gnufile, "set termoption font \"Times,12\"\n");
        }

        /* Write lines independent of plot type */
        fprintf(gnufile, "# Plotting\n");
        fprintf(gnufile, "set style line  1 lt 1 lc rgb '#ff0000'\n"); /* red */
        fprintf(gnufile, "set style line  2 lt 1 lc rgb '#0000ff'\n"); /* blue */
        fprintf(gnufile, "set key at screen 0.74, 0.57 autotitle column box samplen 1 left\n");
        fprintf(gnufile, "set tmargin 0\n");
        fprintf(gnufile, "set bmargin 5\n");
        fprintf(gnufile, "set lmargin 12\n");
        fprintf(gnufile, "set rmargin 3\n");
        fprintf(gnufile, "set xrange [%g:%g]\n", xmin, xmax);
        fprintf(gnufile, "set yrange [%g:%g]\n", ymin, ymax);
        fprintf(gnufile, "unset title\n");
        fprintf(gnufile, "set multiplot layout 2,1 title \"%s\" offset 0,-0.03\n", title);
        fprintf(gnufile, "set xlabel \"%s\"\n", xlabel);
        fprintf(gnufile, "set ylabel \"%s\" offset 1,0\n", ylabel);
        fprintf(gnufile, "set style data boxes\n");
        fprintf(gnufile, "plot '%s' using 1:2 title \"Obs. spectrum\"  with lines ls 1, "
        		              "'%s' using 1:2 title \"Model spectrum\" with lines ls 2\n",
                         tmpfile1, tmpfile2);
        fprintf(gnufile, "set nokey\n");
        fprintf(gnufile, "set yrange [*:*]\n");
        fprintf(gnufile, "set ylabel \"Residual (Obs.-Model)\" offset 1,0\n");
        fprintf(gnufile, "set xlabel \"%s\"\n", xlabel);
        fprintf(gnufile, "plot '%s' using 1:2 title \"Residual\" with lines ls 1\n",
        		         tmpfile3);
        fprintf(gnufile, "unset multiplot\n");

        /* Close temoprary GNUPLOT file */
        fclose(gnufile);

        /* Call GNUPLOT */
        sprintf(systemcall, "gnuplot -persist %s", tmpfile4);
        dummy = system(systemcall);

        /* Remove temporary GNUPLOT file */
        dummy = remove(tmpfile4);

    }

    /* Remove ASCII files with spectra and delete temporary directory */
    dummy = remove(tmpfile1);
    dummy = remove(tmpfile2);
    dummy = remove(tmpfile3);
    dummy = rmdir(tmpdir);

    return CPL_ERROR_NONE;
}

/**@}*/
