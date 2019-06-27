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
 * \file mf_trans.c
 *
 * Routines for deriving a transmission curve for telluric absorption
 * correction
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  10 Jul 2012
 * \date   03 Nov 2014
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/
#define _GNU_SOURCE
#include "mf_trans.h"
#include "mf_molecfit.h"

#if __GNUC__ >= 4 || __clang__
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
#endif

/*****************************************************************************
 *                                GLOBALS                                    *
 ****************************************************************************/

/* Definition of global variables */

/* Reference fit parameters */
cpl_array *reffitpar = NULL;
/* Number of fitting function calls */
extern int nfev;
/* Number of LBLRTM calls (wavenumber-restricted subspectra are not
   counted individually) */
extern int n_code;
/* Total LBLRTM run time */
extern double t_code;


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

static cpl_error_code mf_trans_(mfdrv * drvpar, double wn_start, double wn_end,
                                cpl_boolean read_kernel,
                                cpl_table ** output_spec,
                                cpl_table* prof, cpl_table* res_table);
static cpl_error_code mf_trans_conv(mfdrv * drvpar, double wn_start, double wn_end,
        cpl_boolean read_kernel,
        cpl_table ** output_spec,
        cpl_table* prof, cpl_table* res_table,
        cpl_error_code* arr_err_code);
static cpl_error_code mf_trans_lblrtm(mfdrv * drvpar, double wn_start, double wn_end,
                                      cpl_table* prof, cpl_table* res_table,
                                      cpl_error_code** arr_err_code);
static cpl_error_code save_arr_code_stat(mfdrv* drvpar, cpl_error_code* arr_code_stat);
static cpl_error_code load_arr_code_stat(mfdrv* drvpar, cpl_error_code** arr_code_stat);

cpl_error_code DLL_PUBLIC mf_run_calctrans(
    cpl_parameterlist  *parlist,
    cpl_propertylist   *plist,
    cpl_table          **inspec,
    cpl_size           nspec,
    cpl_table          *molectab,
    cpl_matrix         *kernel,
    double             wl_start,
    double             wl_end,
    cpl_table          *atmprof,
    cpl_table          *res_table,
    mf_calctrans_state **state,
    cpl_table          **result)
{
    /*!
     * \callgraph
     *
     * Run calctrans in temporary folder
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
     * \param wl_start   lower wavelength override for LNFL. -1 to use values from inspec pipeline units, converted via wlgtomicron parameter
     * \param wl_end     lower wavelength override for LNFL. -1 to use values from inspec pipeline units, converted via wlgtomicron parameter
     * \param prof_out   input atmospheric profile  of mf_run_molecfit(...)
     * \param res_out    input results table of mf_run_molecfit(...)
     * \param state      computation state variable
     * \param results    output results table of mf_run_calctrans_calctrans (...)
     */

    mfdrv drvpar;
    cpl_ensure_code(result, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(state, CPL_ERROR_NULL_INPUT);
    *result = NULL;

    /* molecfit changes cwd during run, save and restore */
    char * orig_cwd = mf_get_cwd();

    mf_par_initall(&drvpar);

    mf_override_config(&drvpar,
                       parlist, plist, inspec, nspec,
                       molectab,
                       NULL,
                       NULL,
                       NULL,
                       kernel, *state);

    /* update total wavelength range to allow reusing existing lnfl result for
     * multiple input spectra with different kernels */
    cpl_parameter * par = cpl_parameterlist_find(parlist, "wlgtomicron");
    /* lnfl wave number cm^-1 */
    double wn_start = 1e4 / (wl_end * cpl_parameter_get_double(par));
    double wn_end = 1e4 / (wl_start * cpl_parameter_get_double(par));

    cpl_error_code code = mf_trans_(&drvpar, wn_start, wn_end, kernel == NULL, result, atmprof, res_table);

    if (chdir(orig_cwd)) {}
    cpl_free(orig_cwd);

    return code;

}

cpl_error_code DLL_PUBLIC
mf_run_calctrans_lblrtm(cpl_parameterlist * parlist,
                 cpl_propertylist * plist, cpl_table ** inspec, cpl_size nspec,
                 cpl_table * molectab,
                 double wl_start,
                 double wl_end,
                 cpl_table* atmprof,
                 cpl_table* res_table,
                 mf_calctrans_state ** state)
{
    /*!
     * \callgraph
     *
     * Run calctrans in temporary folder
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
     * \param wl_start lower wavelength override for LNFL
     *                 -1 to use values from inspec
     *                 pipeline units, converted via wlgtomicron parameter
     * \param wl_end lower wavelength override for LNFL
                      -1 to use values from inspec
                      pipeline units, converted via wlgtomicron parameter
     * \param state computation state variable
     *
     */
    mfdrv drvpar;
    cpl_ensure_code(state, CPL_ERROR_NULL_INPUT);

    /* molecfit changes cwd during run, save and restore */
    char * orig_cwd = mf_get_cwd();

    mf_par_initall(&drvpar);

    mf_override_config(&drvpar,
                       parlist, plist, inspec, nspec,
                       molectab,
                       NULL,
                       NULL,
                       NULL,
                       NULL, *state);

    /* update total wavelength range to allow reusing existing lnfl result for
     * multiple input spectra with different kernels */
    cpl_parameter * par = cpl_parameterlist_find(parlist, "wlgtomicron");
    /* lnfl wave number cm^-1 */
    double wn_start = 1e4 / (wl_end * cpl_parameter_get_double(par));
    double wn_end = 1e4 / (wl_start * cpl_parameter_get_double(par));

    cpl_error_code* arr_err_code = NULL;
    cpl_error_code code = mf_trans_lblrtm(&drvpar, wn_start, wn_end, atmprof, res_table, &arr_err_code);
    save_arr_code_stat(&drvpar, arr_err_code);
    if (arr_err_code) {
        free(arr_err_code);
    }

    if (chdir(orig_cwd)) {}
    cpl_free(orig_cwd);

    mf_par_deleteall(&drvpar);

    return code;

}

cpl_error_code DLL_PUBLIC
mf_run_calctrans_convolution(cpl_parameterlist * parlist,
                 cpl_propertylist * plist, cpl_table ** inspec, cpl_size nspec,
                 cpl_matrix * kernel,
                 double wl_start,
                 double wl_end,
                 cpl_table* atmprof,
                 cpl_table* res_table,
                 mf_calctrans_state ** state,
                 cpl_table ** result)
{
    /*!
     * \callgraph
     *
     * Run calctrans in temporary folder
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
     * \param wl_start lower wavelength override for LNFL
     *                 -1 to use values from inspec
     *                 pipeline units, converted via wlgtomicron parameter
     * \param wl_end lower wavelength override for LNFL
                      -1 to use values from inspec
                      pipeline units, converted via wlgtomicron parameter
     * \param state computation state variable
     *
     */
    mfdrv drvpar;
    cpl_ensure_code(result, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(state, CPL_ERROR_NULL_INPUT);
    *result = NULL;

    /* molecfit changes cwd during run, save and restore */
    char * orig_cwd = mf_get_cwd();

    mf_par_initall(&drvpar);

    mf_override_config(&drvpar,
                       parlist, plist, inspec, nspec,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       kernel, *state);

    if (kernel) {
        // Debug
        cpl_msg_info(cpl_func, "Using user-defined Convolution Kernel");
        //cpl_matrix_dump(kernel, NULL);
    }

    /* update total wavelength range to allow reusing existing lnfl result for
     * multiple input spectra with different kernels */
    cpl_parameter * par = cpl_parameterlist_find(parlist, "wlgtomicron");
    /* lnfl wave number cm^-1 */
    double wn_start = 1e4 / (wl_end * cpl_parameter_get_double(par));
    double wn_end = 1e4 / (wl_start * cpl_parameter_get_double(par));

    cpl_error_code* arr_err_code = NULL;
    load_arr_code_stat(&drvpar, &arr_err_code);
    cpl_error_code code = mf_trans_conv(&drvpar, wn_start, wn_end, kernel == NULL, result, atmprof, res_table, arr_err_code);
    if (arr_err_code) {
        free(arr_err_code);
    }

    if (chdir(orig_cwd)) {}
    cpl_free(orig_cwd);

    mf_par_deleteall(&drvpar);

    return code;

}

static cpl_error_code mf_trans_lblrtm(mfdrv * drvpar, double wn_start, double wn_end,
                                      cpl_table* prof, cpl_table* res_table,
                                      cpl_error_code** arr_err_code)
{

    nfev = 1;

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_table *spec;  // , *prof;
    cpl_array *fitpar;

    /* Read spectral and header data from data file */
    spec = cpl_table_new(0);
    if ((status = mf_trans_readspec(spec, drvpar)) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(spec);
        return status;
    }

    /* adapt start and end of rangetab to user provided values
     * used to execute lnfl only once for not fully overlapping spectra */
    if (wn_start > 0 && wn_end > wn_start) {
        for (cpl_size i = 0; i < cpl_table_get_nrow(drvpar->rangetab); i++) {
            double o_wnstart = cpl_table_get(drvpar->rangetab, "wn_start", i, NULL);
            double o_wnend = cpl_table_get(drvpar->rangetab, "wn_end", i, NULL);
            if (o_wnstart <= 0 || o_wnend <= 0) {
                continue;
            }
            cpl_table_set(drvpar->rangetab, "wn_start", i, wn_start);
            cpl_table_set(drvpar->rangetab, "wn_end", i, wn_end);

            /* wn_step not actually used in code TODO remove */
        }
    }

    /* Read best-fit atmospheric profile */
    // prof = cpl_table_new(0);
    // Load only if not provided
    short provided = 1;
    if (! prof) {
        provided = 0;
        if ((status = mf_trans_readatm(&prof, drvpar)) != CPL_ERROR_NONE) {
            mf_par_deleteall(drvpar);
            cpl_table_delete(spec);
            cpl_table_delete(prof);
            return status;
        }
    }

    /* Read the summary of the fit results (ASCII file) and write the fit
       parameters in a CPL array */
    fitpar = cpl_array_new(0, CPL_TYPE_DOUBLE);
    // if ((status = mf_trans_readresults(fitpar, drvpar)) != CPL_ERROR_NONE ||
    if ((status = mf_trans_readresults_fromFits(fitpar, drvpar, res_table)) != CPL_ERROR_NONE ||
        cpl_array_get_size(fitpar) == 0) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(spec);
        cpl_table_delete(prof);
        cpl_array_delete(fitpar);
        return status;
    }

    /* Run LNFL if required */
    if ((status = mf_lnfl(drvpar)) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(spec);
        cpl_table_delete(prof);
        cpl_array_delete(fitpar);
        return status;
    }

    /* Set global reference parameter vector */
    reffitpar = cpl_array_duplicate(fitpar);

    /* Set flag for last modsim call to 1 */
    p = cpl_parameterlist_find(drvpar->parlist, "lastcall");
    cpl_parameter_set_int(p, 1);

    /* Calculate model transmission curve */
    if ((status = mf_modsim_lblrtm(prof, drvpar, fitpar, arr_err_code)) != CPL_ERROR_NONE) {
            mf_trans_clean(drvpar, 0);
            mf_par_deleteall(drvpar);
            cpl_table_delete(spec);
            cpl_table_delete(prof);
            cpl_array_delete(fitpar);
            cpl_array_delete(reffitpar);
            return status;
    }

    /* Free allocated
       memory */
    if (! provided) {
        cpl_table_delete(prof);
    }
    cpl_array_delete(fitpar);
    cpl_table_delete(spec);
    cpl_array_delete(reffitpar);
    reffitpar = NULL;

    /* Return error code of last error */
    return cpl_error_get_code();
}

static cpl_error_code mf_trans_conv(mfdrv * drvpar, double wn_start, double wn_end,
                                    cpl_boolean read_kernel,
                                    cpl_table ** output_spec,
                                    cpl_table* prof, cpl_table* res_table,
                                    cpl_error_code* arr_err_code)
{

    nfev = 1;

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_table *spec;  // , *prof;
    cpl_array *fitpar;

    if (output_spec) {
        *output_spec = NULL;
    }

    /* Read spectral and header data from data file */
    spec = cpl_table_new(0);
    if ((status = mf_trans_readspec(spec, drvpar)) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(spec);
        return status;
    }

    /* adapt start and end of rangetab to user provided values
     * used to execute lnfl only once for not fully overlapping spectra */
    if (wn_start > 0 && wn_end > wn_start) {
        for (cpl_size i = 0; i < cpl_table_get_nrow(drvpar->rangetab); i++) {
            double o_wnstart = cpl_table_get(drvpar->rangetab, "wn_start", i, NULL);
            double o_wnend = cpl_table_get(drvpar->rangetab, "wn_end", i, NULL);
            if (o_wnstart <= 0 || o_wnend <= 0) {
                continue;
            }
            cpl_table_set(drvpar->rangetab, "wn_start", i, wn_start);
            cpl_table_set(drvpar->rangetab, "wn_end", i, wn_end);

            /* wn_step not actually used in code TODO remove */
        }
    }

    /* Read fixed kernel if provided */
    if (read_kernel &&
        (status = mf_par_readkernel(drvpar)) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(spec);
        return status;
    }

    /* Read best-fit atmospheric profile */
    // prof = cpl_table_new(0);
    // Load only if not provided
    short provided = 1;
    if (! prof) {
        provided = 0;
        if ((status = mf_trans_readatm(&prof, drvpar)) != CPL_ERROR_NONE) {
            mf_par_deleteall(drvpar);
            cpl_table_delete(spec);
            cpl_table_delete(prof);
            return status;
        }
    }

    /* Read the summary of the fit results (ASCII file) and write the fit
       parameters in a CPL array */
    fitpar = cpl_array_new(0, CPL_TYPE_DOUBLE);
    // if ((status = mf_trans_readresults(fitpar, drvpar)) != CPL_ERROR_NONE ||
    if ((status = mf_trans_readresults_fromFits(fitpar, drvpar, res_table)) != CPL_ERROR_NONE ||
        cpl_array_get_size(fitpar) == 0) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(spec);
        cpl_table_delete(prof);
        cpl_array_delete(fitpar);
        return status;
    }

    /* Set flag for last modsim call to 1 */
    p = cpl_parameterlist_find(drvpar->parlist, "lastcall");
    cpl_parameter_set_int(p, 1);

    if ((status = mf_modsim_conv(spec, drvpar, fitpar, arr_err_code)) != CPL_ERROR_NONE) {
            mf_trans_clean(drvpar, 0);
            mf_par_deleteall(drvpar);
            cpl_table_delete(spec);
            cpl_table_delete(prof);
            cpl_array_delete(fitpar);
            cpl_array_delete(reffitpar);
            reffitpar = NULL;
            return status;
    }

    /* Free allocated
       memory */
    if (! provided) {
       cpl_table_delete(prof);
    }
    cpl_array_delete(fitpar);

    /* Calculate function for telluric absorption correction and perform this
       correction for the input spectrum */
    mf_trans_calctac(spec, drvpar);

    /* Write correction function for telluric absorption to FITS file in
       output folder */
    mf_trans_writefile(spec, drvpar);

    /* Plot input spectrum and resulting spectrum from telluric absoprtion
       correction */
    if (status == CPL_ERROR_NONE) {
        mf_trans_plot(spec, drvpar);
    }

    /* Write results into file of same format as input data file */
    if (status == CPL_ERROR_NONE) {
        mf_trans_writeresults(drvpar);
    }

    if (output_spec) {
        *output_spec = cpl_table_duplicate(spec);
    }

    /* Free allocated memory */
    cpl_table_delete(spec);

    /* Return error code of last error */
    return cpl_error_get_code();
}

static cpl_error_code mf_trans_(mfdrv * drvpar, double wn_start, double wn_end,
                                cpl_boolean read_kernel,
                                cpl_table ** output_spec,
                                cpl_table* prof, cpl_table* res_table)
{

    nfev = 1;

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_table *spec;  // , *prof;
    cpl_array *fitpar;

    if (output_spec) {
        *output_spec = NULL;
    }

    /* Read spectral and header data from data file */
    spec = cpl_table_new(0);
    if ((status = mf_trans_readspec(spec, drvpar)) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(spec);
        return status;
    }

    /* adapt start and end of rangetab to user provided values
     * used to execute lnfl only once for not fully overlapping spectra */
    if (wn_start > 0 && wn_end > wn_start) {
        for (cpl_size i = 0; i < cpl_table_get_nrow(drvpar->rangetab); i++) {
            double o_wnstart = cpl_table_get(drvpar->rangetab, "wn_start", i, NULL);
            double o_wnend = cpl_table_get(drvpar->rangetab, "wn_end", i, NULL);
            if (o_wnstart <= 0 || o_wnend <= 0) {
                continue;
            }
            cpl_table_set(drvpar->rangetab, "wn_start", i, wn_start);
            cpl_table_set(drvpar->rangetab, "wn_end", i, wn_end);

            /* wn_step not actually used in code TODO remove */
        }
    }

    /* Read fixed kernel if provided */
    if (read_kernel &&
        (status = mf_par_readkernel(drvpar)) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(spec);
        return status;
    }

    /* Read best-fit atmospheric profile */
    // prof = cpl_table_new(0);
    // Load only if not provided
    short provided = 1;
    if (! prof) {
        provided = 0;
        if ((status = mf_trans_readatm(&prof, drvpar)) != CPL_ERROR_NONE) {
            mf_par_deleteall(drvpar);
            cpl_table_delete(spec);
            cpl_table_delete(prof);
            return status;
        }
    }

    /* Read the summary of the fit results (ASCII file) and write the fit
       parameters in a CPL array */
    fitpar = cpl_array_new(0, CPL_TYPE_DOUBLE);
    // if ((status = mf_trans_readresults(fitpar, drvpar)) != CPL_ERROR_NONE ||
    if ((status = mf_trans_readresults_fromFits(fitpar, drvpar, res_table)) != CPL_ERROR_NONE ||
        cpl_array_get_size(fitpar) == 0) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(spec);
        cpl_table_delete(prof);
        cpl_array_delete(fitpar);
        return status;
    }

    /* Run LNFL if required */
    if ((status = mf_lnfl(drvpar)) != CPL_ERROR_NONE) {
        mf_par_deleteall(drvpar);
        cpl_table_delete(spec);
        cpl_table_delete(prof);
        cpl_array_delete(fitpar);
        return status;
    }

    /* Set global reference parameter vector */
    reffitpar = cpl_array_duplicate(fitpar);

    /* Set flag for last modsim call to 1 */
    p = cpl_parameterlist_find(drvpar->parlist, "lastcall");
    cpl_parameter_set_int(p, 1);

    /* Calculate model transmission curve */
    cpl_table* prof_out = NULL;
    if ((status = mf_modsim(spec, &prof_out, prof, drvpar, fitpar)) != CPL_ERROR_NONE) {
            mf_trans_clean(drvpar, 1);
            mf_par_deleteall(drvpar);
            cpl_table_delete(spec);
            cpl_table_delete(prof);
            cpl_array_delete(fitpar);
            cpl_array_delete(reffitpar);
            reffitpar = NULL;
            return status;
    }
    if (prof_out) cpl_table_delete(prof_out);

    /* Remove files related to radiative transfer code and free allocated
       memory */
    mf_trans_clean(drvpar, 1);
    if (! provided) {
        cpl_table_delete(prof);
    }
    cpl_array_delete(fitpar);
    cpl_array_delete(reffitpar);
    reffitpar = NULL;

    /* Calculate function for telluric absorption correction and perform this
       correction for the input spectrum */
    mf_trans_calctac(spec, drvpar);

    /* Write correction function for telluric absorption to FITS file in
       output folder */
    if ((status =mf_trans_writefile(spec, drvpar)) != CPL_ERROR_NONE) {
        return status;
    }

    /* Plot input spectrum and resulting spectrum from telluric absoprtion
       correction */
    if (status == CPL_ERROR_NONE) {
        mf_trans_plot(spec, drvpar);
    }

    /* Write results into file of same format as input data file */
    if (status == CPL_ERROR_NONE) {
        mf_trans_writeresults(drvpar);
    }

    if (output_spec) {
        *output_spec = cpl_table_duplicate(spec);
    }

    /* Free allocated memory */
    mf_par_deleteall(drvpar);
    cpl_table_delete(spec);

    /* Return error code of last error */
    return cpl_error_get_code();
}

static cpl_error_code save_arr_code_stat(mfdrv* drvpar, cpl_error_code* arr_code_stat) {

    cpl_parameter *p;
    // char basedir[MF_MAXLEN];
    char output_dir[MF_MAXLEN];
    // const char* output_name = "arr_code_stat";
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    char output_name[MF_MAXLEN];
    strncpy(output_name, cpl_parameter_get_string(p), MF_MAXLEN);

    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(output_dir, cpl_parameter_get_string(p), curdir);
    char outfile[MF_MAXLEN];
    sprintf(outfile, "%s%s_code_stat.flags", output_dir, output_name);

    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    int nrange = cpl_parameter_get_int(p);

    FILE* file = fopen(outfile, "w");
    for (int index = 0; index < nrange; ++index) {
        fprintf(file, "%d\n", arr_code_stat[index]);
    }
    fclose(file);

    return CPL_ERROR_NONE;
}

static cpl_error_code load_arr_code_stat(mfdrv* drvpar, cpl_error_code** arr_code_stat) {

    cpl_parameter *p;
    // char basedir[MF_MAXLEN];
    char output_dir[MF_MAXLEN];
    // const char* output_name = "arr_code_stat";
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    char output_name[MF_MAXLEN];
    strncpy(output_name, cpl_parameter_get_string(p), MF_MAXLEN);

    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(output_dir, cpl_parameter_get_string(p), curdir);
    char outfile[MF_MAXLEN];
    sprintf(outfile, "%s%s_code_stat.flags", output_dir, output_name);

    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    int nrange = cpl_parameter_get_int(p);

    *arr_code_stat = malloc(nrange * sizeof(cpl_error_code));

    FILE* file = fopen(outfile, "r");

    int index = 0;
    char* line = NULL;
    size_t len = 0;
    while (1) {
        ssize_t read = getline(&line, &len, file);
        if (read == -1) {
            break;
        }

        cpl_free(line);
        (*arr_code_stat)[index++] = CPL_ERROR_NONE;
        if (index >= nrange) {
            break;
        }
    }

    fclose(file);

    return CPL_ERROR_NONE;
}

cpl_error_code DLL_PUBLIC clean_arr_code_stat(cpl_parameterlist * parlist) {

    cpl_parameter *p;
    // char basedir[MF_MAXLEN];
    char output_dir[MF_MAXLEN];
    // const char* output_name = "arr_code_stat";
    p = cpl_parameterlist_find(parlist, "output_name");
    char output_name[MF_MAXLEN];
    strncpy(output_name, cpl_parameter_get_string(p), MF_MAXLEN);

    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(parlist, "output_dir");
    mf_basic_abspath(output_dir, cpl_parameter_get_string(p), curdir);
    char outfile[MF_MAXLEN];
    sprintf(outfile, "%s%s_code_stat.flags", output_dir, output_name);

    char cmd[MF_MAXLEN];
    sprintf(cmd, "rm %s", outfile);

    int d;
    if ((d = system(cmd))){};

    return CPL_ERROR_NONE;
}

cpl_error_code mf_trans_lblrtm_(const char *parfile_in) {

    cpl_error_code status = CPL_ERROR_NONE;
    mfdrv drvpar;
    char parfile[MF_MAXLEN];

    /* Read MOLECFIT driver file */
    mf_par_initall(&drvpar);

    mf_basic_absfile(parfile, parfile_in);
    if ((status = mf_par_readfile(&drvpar, parfile, 1)) != CPL_ERROR_NONE) {
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

    cpl_error_code* arr_code_stat = NULL;
    status = mf_trans_lblrtm(&drvpar, -1., -1., NULL, NULL, &arr_code_stat);

    // Save array of error codes
    status = save_arr_code_stat(&drvpar, arr_code_stat);
    if (arr_code_stat) {
        free(arr_code_stat);
    }

    // Remove LBLRTM-related files
    mf_cleanup_standalone(tmpdir);

    mf_par_deleteall(&drvpar);

    return status;
}

cpl_error_code mf_trans_conv_(const char *parfile_in) {

    cpl_error_code status = CPL_ERROR_NONE;
    mfdrv drvpar;
    char parfile[MF_MAXLEN];

    /* Read MOLECFIT driver file */
    mf_par_initall(&drvpar);
    mf_basic_absfile(parfile, parfile_in);
    if ((status = mf_par_readfile(&drvpar, parfile, 1)) != CPL_ERROR_NONE) {
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

    cpl_error_code* arr_code_stat = NULL;
    status = load_arr_code_stat(&drvpar, &arr_code_stat);

    status = mf_trans_conv(&drvpar, -1., -1., CPL_TRUE, NULL, NULL, NULL, arr_code_stat);

    if (arr_code_stat) {
        free(arr_code_stat);
    }

    // mf_trans_clean remove files produced by lblrtm
    mf_cleanup_standalone(tmpdir);
    mf_trans_clean(&drvpar, 0);
    mf_par_deleteall(&drvpar);

    return status;
}

cpl_error_code mf_trans(const char *parfile_in)
{
    /*!
     * \callgraph
     *
     * Derives transmission curve for entire wavelength range of input data by
     * using the fit results of MOLECFIT. Moreover, the input spectrum is
     * corrected for telluric absorption by using this transmission curve. The
     * fit data is identified by the input parameter file. As output, a
     * standard FITS table and a file of the same format as the input data
     * file is written. In the case of FITS images, the transmission curve and
     * the telluric absorption corrected spectrum are written into separate
     * files. Finally, a diagnostic plot is optionally created.
     *
     * \b INPUT:
     * \param parfile_in  name of parameter file
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    mfdrv drvpar;
    char parfile[MF_MAXLEN];

    /* Read MOLECFIT driver file */
    mf_par_initall(&drvpar);

    mf_basic_absfile(parfile, parfile_in);
    if ((status = mf_par_readfile(&drvpar, parfile, 1)) != CPL_ERROR_NONE) {
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

    status = mf_trans_(&drvpar, -1., -1., CPL_TRUE, NULL, NULL, NULL);

    mf_cleanup_standalone(tmpdir);

    return status;
}


cpl_error_code mf_trans_readspec(cpl_table *spec, mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Reads a spectrum in table format and header keyworks if provided.
     * This routine is similar to ::mf_readspec, except for neglecting
     * wavelength ranges provided by additional files. CALCTRANS is forced to
     * compute a transmission curve for the full wavelength range covered by
     * the input spectrum.
     *
     * \b INPUT:
     * \param spec    empty CPL table
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    CPL table with wavelengths, fluxes, and weights
     * \param drvpar  ::mfdrv parameter structure with FITS header values
     *                if the observing instrument is given
     *
     * \b ERRORS:
     * - see ::mf_readspec
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    char wrangefile[MF_MAXLEN], prangefile[MF_MAXLEN];

    /* Change parameters to calculate transmission curve for full spectrum */
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    if (cpl_parameter_get_int(p) == 0) {
        cpl_parameter_set_int(p, 2);
    }
    p = cpl_parameterlist_find(drvpar->parlist, "wrange_include");
    cpl_parameter_set_string(p, "none");
    p = cpl_parameterlist_find(drvpar->parlist, "wrange_exclude");
    strncpy(wrangefile, cpl_parameter_get_string(p), MF_MAXLEN);
    cpl_parameter_set_string(p, "none");
    p = cpl_parameterlist_find(drvpar->parlist, "prange_exclude");
    strncpy(prangefile, cpl_parameter_get_string(p), MF_MAXLEN);
    cpl_parameter_set_string(p, "none");

    /* Read input spectrum and adapt mfdrv parameter structure */
    status = mf_readspec(spec, drvpar);

    /* Set weights to 0 for wavelength ranges provided by a file */
    p = cpl_parameterlist_find(drvpar->parlist, "wrange_exclude");
    cpl_parameter_set_string(p, wrangefile);
    mf_readspec_excluderanges(spec, drvpar, 'w');

    /* Set weights to 0 for pixel ranges provided by a file */
    p = cpl_parameterlist_find(drvpar->parlist, "prange_exclude");
    cpl_parameter_set_string(p, prangefile);
    mf_readspec_excluderanges(spec, drvpar, 'p');

    return status;
}


cpl_error_code mf_trans_readatm(cpl_table **prof, const mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Reads best-fit atmospheric profile and writes it into a CPL table. Path
     * and file name are provided by the input ::mfdrv parameter structure.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param prof    CPL table with best-fit profile.
     *
     * \b ERRORS:
     * - see ::mf_atm_readatm
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    char outdir[MF_MAXLEN], outname[MF_MAXLEN];
    char atmfile[MF_MAXLEN];

    /* Get output folder and name space */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);

    /* Compose full path and name of best-fit atmospheric profile */
    // sprintf(atmfile, "%s%s_fit.atm", outdir, outname);
    sprintf(atmfile, "%s%s_fit.atm.fits", outdir, outname);

    /* Write info message */
    cpl_msg_info(cpl_func, "Read atmospheric profile file %s", atmfile);

    /* Read atmospheric profile */
    // status = mf_atm_readatm(prof, atmfile);
    status = mf_atm_readatm_fromFits(prof, atmfile);
    if (status != CPL_ERROR_NONE) {
        char errtxt[MF_MAXLEN];
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, atmfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    return status;
}

cpl_error_code mf_trans_readresults_fromFits(cpl_array *fitpar, const mfdrv *drvpar, cpl_table* res_table)
{

    cpl_parameter *p;
    char outdir[MF_MAXLEN], outname[MF_MAXLEN];
    char fitsfilename[MF_MAXLEN];

    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);

    /* Get number of fit parameters */
    p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    int nmolec = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "nchip");
    int nchip = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "wlc_n");
    int nwlc = cpl_parameter_get_int(p) + 1;
    p = cpl_parameterlist_find(drvpar->parlist, "cont_n");
    int ncont = cpl_parameter_get_int(p) + 1;
    int npar = nmolec + (nwlc + ncont) * nchip + 3;

    /* Set size of output CPL array */
    cpl_array_set_size(fitpar, npar);

    cpl_array_fill_window_double(fitpar, 0, nmolec, 1.);

    /* Set constant of continuum correction function to 1 and all other
       coefficients to 0 to obtain unscaled transmission curve */
    int contparmin = nmolec + nchip * nwlc;
    int contparmax = contparmin + nchip * ncont - 1;
    for (int j = -1, i = contparmin; i <= contparmax; i++) {
        j++;
        if (j == ncont) {
            j = 0;
        }
        if (j == 0) {
            cpl_array_set(fitpar, i, 1);
        } else {
            cpl_array_set(fitpar, i, 0);
        }
    }

    /* Compose full path and name of file with fit results */
    sprintf(fitsfilename, "%s%s_fit.res.fits", outdir, outname);
    // cpl_table* res_table;

    // Load only if not provided
    short provided = 1;
    if (! res_table) {
        provided = 0;
        cpl_msg_info(cpl_func, "Read fit results file %s", fitsfilename);
        if(!(res_table = cpl_table_load(fitsfilename, 1, 0))) {
            cpl_error_code status = cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                         "Could not open fits file: %s",
                                         fitsfilename);
            cpl_table_delete(res_table);
            return status;
        }
    }
    int wlcparmin = nmolec;
    int wlcparmax = wlcparmin + nchip * nwlc - 1;
    int fitCoefIndex = wlcparmin;

    int status = -9999;
    double boxfwhm = -1;
    double gaussfwhm = -1;
    double lorentzfwhm = -1;
    printf("[ DEBUG ] wlcparmin = %d \n", wlcparmin);
    printf("[ DEBUG ] wlcparmax = %d \n", wlcparmax);
    printf("[ DEBUG ] nwlc=%d\n", nwlc);
    for (unsigned index = 0; index < cpl_table_get_nrow(res_table); ++index) {
        char const* paramName = cpl_table_get_string(res_table, "parameter", index);
        printf("[ DEBUG ] index=%d %s\n", index, paramName);
        if (!strncmp("status", paramName, 6)) {
            status = (int) cpl_table_get_double(res_table, "value", index, NULL);
            if (status < 1) {
                cpl_array_set_size(fitpar, 0);
                cpl_msg_info(cpl_func, "No fit results available -> Return");
                return CPL_ERROR_NONE;
            }
        }
        else if (!strncmp("boxfwhm", paramName, 7)) {
            boxfwhm = cpl_table_get_double(res_table, "value", index, NULL);
            cpl_array_set(fitpar, npar-3, boxfwhm);
        }
        else if (!strncmp("gaussfwhm", paramName, 9)) {
            gaussfwhm = cpl_table_get_double(res_table, "value", index, NULL);
            cpl_array_set(fitpar, npar-2, gaussfwhm);
        }
        else if (!strncmp("lorentzfwhm", paramName, 11)) {
            lorentzfwhm = cpl_table_get_double(res_table, "value", index, NULL);
            cpl_array_set(fitpar, npar-1, lorentzfwhm);
        }
        else if (!strncmp("Chip", paramName, 4)) {
            double value = cpl_table_get_double(res_table, "value", index, NULL);
            printf("[ DEBUG ] param=%s, index=%d, fitCoefIndex=%d\n", paramName, index, fitCoefIndex);
            cpl_array_set(fitpar, fitCoefIndex++, value);
            if(fitCoefIndex > wlcparmax + 1) {
                char errtxt[MF_MAXLEN];
                sprintf(errtxt, "%s: %s Too many coef. for wavelength solution",
                    MF_ERROR_UFS_TXT, fitsfilename);
                return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
            }
        }
    }

    if (status == -9999) {
        char errtxt[MF_MAXLEN];
        sprintf(errtxt, "%s: %s ('Status' not found)",
                    MF_ERROR_UFS_TXT, fitsfilename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
    }
    if (boxfwhm < 0) {
        char errtxt[MF_MAXLEN];
        sprintf(errtxt, "%s: %s ('FWHM of Gaussian' not found)",
                MF_ERROR_UFS_TXT, fitsfilename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }
    if (lorentzfwhm < 0) {
        char errtxt[MF_MAXLEN];
        sprintf(errtxt, "%s: %s ('FWHM of Lorentzian' not found)",
                MF_ERROR_UFS_TXT, fitsfilename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }
    if (fitCoefIndex == wlcparmin) {
        char errtxt[MF_MAXLEN];
        sprintf(errtxt, "%s: %s (coef. for wavelength solution not found)",
                MF_ERROR_UFS_TXT, fitsfilename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    if (! provided) {
        cpl_table_delete(res_table);
    }

    return CPL_ERROR_NONE;

}

cpl_error_code mf_trans_readresults(cpl_array *fitpar, const mfdrv *drvpar)
{
    /*!
     * Reads best-fit parameter values from ASCII file for fit results and
     * writes them into a CPL array. For continuum-related parameters, the
     * default values are taken. Path and file name are provided by the
     * input ::mfdrv parameter structure.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param fitpar  CPL array with best-fit parameter values.
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    cpl_parameter *p;
    mfpar par[MF_MAXPAR];
    char outdir[MF_MAXLEN], outname[MF_MAXLEN];
    char resfile[MF_MAXLEN], errtxt[MF_MAXLEN];
    const char *word1 = "Status:", *word2 = "SPECTRAL";
    int nmolec = 0, nchip = 0, nwlc = 0, ncont = 0, npar = 0, contparmin = 0;
    int contparmax = 0, j = 0, i = 0, nel = 0, wlcparmin = 0, wlcparmax = 0;
    size_t wordlen;

    /* Get number of fit parameters */
    p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "nchip");
    nchip = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "wlc_n");
    nwlc = cpl_parameter_get_int(p) + 1;
    p = cpl_parameterlist_find(drvpar->parlist, "cont_n");
    ncont = cpl_parameter_get_int(p) + 1;
    npar = nmolec + (nwlc + ncont) * nchip + 3;

    /* Set size of output CPL array */
    cpl_array_set_size(fitpar, npar);

    /* Set molecular parameters to 1 (no change of best-fit atmospheric
       profile) */
    cpl_array_fill_window_double(fitpar, 0, nmolec, 1.);

    /* Set constant of continuum correction function to 1 and all other
       coefficients to 0 to obtain unscaled transmission curve */
    contparmin = nmolec + nchip * nwlc;
    contparmax = contparmin + nchip * ncont - 1;
    for (j = -1, i = contparmin; i <= contparmax; i++) {
        j++;
        if (j == ncont) {
            j = 0;
        }
        if (j == 0) {
            cpl_array_set(fitpar, i, 1);
        } else {
            cpl_array_set(fitpar, i, 0);
        }
    }

    /* Get output folder and name space */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);

    /* Compose full path and name of file with fit results */
    sprintf(resfile, "%s%s_fit.res", outdir, outname);

    /* Write info message */
    cpl_msg_info(cpl_func, "Read results file %s", resfile);

    /* Open results file if it exists */
    if ((stream = fopen(resfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, resfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Get length of search word */
    wordlen = strlen(word1);

    /* Search "Status:" in file */
    do {
        mf_basic_readline(stream, par, &nel);
        if (feof(stream) != 0) {
            fclose(stream);
            strcpy(par[0].c, "#");
            sprintf(errtxt, "%s: %s ('Status' not found)",
                    MF_ERROR_UFS_TXT, resfile);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }
    } while (par[0].c[0] == '#' || par[0].c[0] == '\n' ||
             strncmp(par[0].c, word1, wordlen) != 0);

    /* Return in the case of an erroneous fit */
    if (par[1].d <= 0) {
        cpl_array_set_size(fitpar, 0);
        cpl_msg_info(cpl_func, "No fit results available -> Return");
        return CPL_ERROR_NONE;
    }

    /* Get length of search word */
    wordlen = strlen(word2);

    /* Search line "SPECTRAL RESOLUTION:" in file */
    do {
        mf_basic_readline(stream, par, &nel);
        if (feof(stream) != 0) {
            fclose(stream);
            strcpy(par[0].c, "#");
            sprintf(errtxt, "%s: %s ('SPECTRAL RESOLUTION' not found)",
                    MF_ERROR_UFS_TXT, resfile);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }
    } while (par[0].c[0] == '#' || par[0].c[0] == '\n' ||
             strncmp(par[0].c, word2, wordlen) != 0);

    /* Set relative FWHM of boxcar */
    mf_basic_readline(stream, par, &nel);
    if (nel != 9 && nel != 11) {
        fclose(stream);
        sprintf(errtxt, "%s: %s ('Rel. FWHM of boxcar' not found)",
                MF_ERROR_UFS_TXT, resfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }
    cpl_array_set(fitpar, npar-3, par[8].d);

    /* Skip one line */
    mf_basic_readline(stream, par, &nel);

    /* Set FWHM of Gaussian */
    mf_basic_readline(stream, par, &nel);
    if (nel != 6 && nel != 8) {
        fclose(stream);
        sprintf(errtxt, "%s: %s ('FWHM of Gaussian' not found)",
                MF_ERROR_UFS_TXT, resfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }
    cpl_array_set(fitpar, npar-2, par[5].d);

    /* Set FWHM of Lorentzian */
    mf_basic_readline(stream, par, &nel);
    if (nel != 6 && nel != 8) {
        fclose(stream);
        sprintf(errtxt, "%s: %s ('FWHM of Lorentzian' not found)",
                MF_ERROR_UFS_TXT, resfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }
    cpl_array_set(fitpar, npar-1, par[5].d);

    /* Skip two lines (one is expected to be empty) */
    mf_basic_readline(stream, par, &nel);

    /* Set coefficients for correction of wavelength grid */
    wlcparmin = nmolec;
    wlcparmax = wlcparmin + nchip * nwlc - 1;
    for (i = wlcparmin; i <= wlcparmax; i++) {
        mf_basic_readline(stream, par, &nel);
        if (nel != 5 && nel != 7) {
            fclose(stream);
            sprintf(errtxt,
                    "%s: %s (coef. for wavelength solution not found)",
                    MF_ERROR_UFS_TXT, resfile);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }
        cpl_array_set(fitpar, i, par[4].d);
    }

    /* Close results file */
    fclose(stream);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_trans_clean(const mfdrv *drvpar, const char cleanSpec)
{
    /*!
     * Removes the files produced by the selected radiative transfer code.
     * Only transmission spectra are expected and hence deleted.
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
    char outdir[MF_MAXLEN], outspec[MF_MAXLEN];
    char sys[MF_MAXLEN];

    /* Get output folder and name space */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);

    /* Get names of transmission spectra (radiative transfer code) and remove
       files if present */
    if (cleanSpec) {
        // p = cpl_parameterlist_find(drvpar->parlist, "spec_out");
        p = cpl_parameterlist_find(drvpar->parlist, "output_name");
        sprintf(outspec, "%s%s_*_T.fits", outdir, cpl_parameter_get_string(p));
        sprintf(sys, "rm -rf %s", outspec);
        if (system(sys)) {};
    }
    return CPL_ERROR_NONE;
}


cpl_error_code mf_trans_calctac(cpl_table *spec, const mfdrv *drvpar)
{
    /*!
     * Provides model transmission curve in "mtrans" column of output
     * spectrum. Moreover, this function is used to correct the observed flux
     * for telluric absorption. The corrected fluxes are given in column
     * "cflux". The column "qual" indicates whether the correction could be
     * performed for a pixel (= 1) or whether it failed (= 0). No flux
     * correction is carried out in the case of a sky emission spectrum. In
     * this case, the column "cflux" contains the uncorrected input flux.
     *
     * \b INPUT:
     * \param spec    CPL table with observed and modelled spectrum
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    table with columns for model transmission and observed
     *                flux corrected for telluric absorption
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    cpl_parameter *p;
    char errtxt[MF_MAXLEN];
    int trans = 0;
    double minscal = 0., maxscal = 0.;

    /* Check whether a transmission was calculated */
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);
    if (trans == 0) {
        sprintf(errtxt, "%s: trans of mfdrv *drvpar "
                "(not transmission as required)", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Check whether there was no continuum correction */
    minscal = cpl_table_get_column_min(spec, "mscal");
    maxscal = cpl_table_get_column_max(spec, "mscal");
    if (maxscal - minscal > MF_TOL) {
        sprintf(errtxt, "%s: column mscal of cpl_table *spec != 1 "
                "(unwanted continuum correction)", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Write info message */
    if (trans == 2) {
        /* Sky emission spectrum */
        cpl_msg_info(cpl_func, "Sky emission spectrum: no flux correction");
    } else if (trans == 1) {
        /* Transmission spectrum */
        cpl_msg_info(cpl_func, "Perform telluric absorption correction");
    }

    /* Rename flux column */
    cpl_table_name_column(spec, "mflux", "mtrans");

    /* Remove unnecessary columns in spectrum table */
    cpl_table_erase_column(spec, "mrange");
    cpl_table_erase_column(spec, "mscal");
    cpl_table_erase_column(spec, "dev");

    /* Corrects flux for telluric absorption and writes it in new column
       "cflux". The quality of the correction is indicated by the column
       "qual" (0 = bad, 1 = acceptable) */
    mf_corr_calctac(spec, trans);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_trans_writefile(cpl_table *spec, const mfdrv *drvpar)
{
    /*!
     * Writes CPL table with observed flux (column "flux"), best-fit
     * atmospheric transmission (column "mtrans"), and telluric absorption
     * corrected flux (colum "cflux") into an ASCII and a FITS file marked by
     * the suffix "_tac". Further output columns are chip ("chip"),
     * wavelength ("lambda"), weight of observed data ("weight"), model weight
     * ("mweight"), and quality flag for telluric absorption correction
     * ("qual"). Writing the ASCII file can be switched off by setting the
     * parameter \e writeascii to 0.
     *
     * \b INPUT:
     * \param spec     CPL table for telluric absorption correction
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
    cpl_error_code status = CPL_ERROR_NONE;

    /* Get output folder and name space */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);

    /* Write info message */
    cpl_msg_info(cpl_func, "Write results into output folder %s", outdir);

    /* Write ASCII file if desired */
    p = cpl_parameterlist_find(drvpar->parlist, "writeascii");
    if (cpl_parameter_get_int(p) == 1) {
        sprintf(outdat, "%s%s_tac.asc", outdir, outname);
        stream = fopen(outdat, "w+");
        cpl_table_dump(spec, 0, cpl_table_get_nrow(spec), stream);
        fclose(stream);
    }

    /* Write FITS file */
    sprintf(outfits, "%s%s_tac.fits", outdir, outname);
    if ((status = cpl_table_save(spec, NULL, NULL, outfits, CPL_IO_CREATE)) != CPL_ERROR_NONE) {
        char errtxt[MF_MAXLEN];
        sprintf(errtxt, "Could not write telluric absorption FITS file %s", outfits);
        return cpl_error_set_message(cpl_func, MF_ERROR_IO, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_trans_plot(const cpl_table *spec, const mfdrv *drvpar)
{
    /*!
     * Compares a spectrum corrected for telluric absorption with the
     * uncorrected spectrum. GNUPLOT is used for plotting.
     *
     * \b INPUT:
     * \param spec    CPL table with observed and corrected spectrum
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    FILE *specfile1, *specfile2, *gnufile;
    cpl_parameter *p;
    char errtxt[MF_MAXLEN], tmpdir[MF_MAXLEN], tmpfile1[MF_MAXLEN];
    char tmpfile2[MF_MAXLEN], plottype[MF_MAXLEN];
    char outdir[MF_MAXLEN], outname[MF_MAXLEN];
    char psfile[MF_MAXLEN], tmpfile3[MF_MAXLEN];
    char systemcall[MF_MAXLEN];
    int nrow = 0, trans = 0, i = 0, dummy = 0, j = 0;
    int plotopt = 2;
    double xmin = 0., xmax = 0., ymin = 0., ymax = 0., dy = 0.;
    double lam = 0., flux1 = 0., weight1 = 0., flux2 = 0.;

    /* Labels of required table columns */
    char collam[] = "lambda";
    char colflux1[] = "flux";
    char colweight1[] = "weight";
    char colflux2[] = "cflux";
    char coltrans[] = "mtrans";

    /* Extension of y-axis in per cent */
    double del = 0.05;

    /* Temporary filenames */
    char filename1[] = "obsspec.dat";
    char filename2[] = "tacspec.dat";
    char gnuname[] = "plot.gnu";

    /* Plot labels */
    char xlabel[] = "Wavelength [micron]";
    char ylabel[] = "Radiance";
    char title[] = "Quality of telluric absorption correction";

    /* Check existence of required columns */
    if (cpl_table_has_column(spec, collam) != 1 ||
        cpl_table_has_column(spec, colflux1) != 1 ||
        cpl_table_has_column(spec, colweight1) != 1 ||
        cpl_table_has_column(spec, colflux2) != 1 ||
        cpl_table_has_column(spec, coltrans) != 1) {
        sprintf(errtxt, "%s: cpl_table *spec (required columns not found)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get wavelength range for x-axis */
    nrow = cpl_table_get_nrow(spec);
    xmin = cpl_table_get(spec, collam, 0, NULL);
    xmax = cpl_table_get(spec, collam, nrow-1, NULL);

    /* Sky emission: no flux correction */
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);

    /* Get flux range for y-axis */
    ymin = 0.;
    for (ymax = 0., i = 0; i < nrow; i++) {
        flux1 = cpl_table_get(spec, colflux1, i, NULL);
        weight1 = cpl_table_get(spec, colweight1, i, NULL);
        if (flux1 > ymax && weight1 > 0.) {
            ymax = flux1;
        }
    }
    dy = ymax - ymin;
    ymin -= del * dy;
    if (trans == 1) {
        ymax /= cpl_table_get_column_median(spec, coltrans);
    }
    ymax += del * dy;

    /* Create temporary directory */
    sprintf(tmpdir, "__tmpDIRtmp__");
    if (access(tmpdir, W_OK) == 0) {
        cpl_msg_warning(cpl_func, "Directory %s already exists!", tmpdir);
    } else {
        if ((dummy = mkdir(tmpdir, 0777))) {};
    }

    /* Write ASCII files containing observed and corrected spectra */
    sprintf(tmpfile1, "%s/%s", tmpdir, filename1);
    specfile1 = fopen(tmpfile1, "w");
    sprintf(tmpfile2, "%s/%s", tmpdir, filename2);
    specfile2 = fopen(tmpfile2, "w");
    for (i = 0; i < nrow; i++) {
        lam = cpl_table_get(spec, collam, i, NULL);
        flux1 = cpl_table_get(spec, colflux1, i, NULL);
        flux2 = cpl_table_get(spec, colflux2, i, NULL);
        fprintf(specfile1, "%5.6g\t%5.6g\n", lam, flux1);
        fprintf(specfile2, "%5.6g\t%5.6g\n", lam, flux2);
    }
    fclose(specfile1);
    fclose(specfile2);

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
        sprintf(psfile, "%s%s_tac.ps", outdir, outname);
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
        mf_basic_initstring(tmpfile3, MF_MAXLEN);
        sprintf(tmpfile3, "%s/%s", tmpdir, gnuname);
        gnufile = fopen(tmpfile3, "w");

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
        }

        /* Write lines independent of plot type */
        fprintf(gnufile, "set nokey\n");
        fprintf(gnufile, "set tmargin 2\n");
        fprintf(gnufile, "set bmargin 5\n");
        fprintf(gnufile, "set lmargin 13\n");
        fprintf(gnufile, "set rmargin 3\n");
        fprintf(gnufile, "set xrange [ %g : %g ]\n", xmin, xmax);
        fprintf(gnufile, "set yrange [ %g : %g ]\n", ymin, ymax);
        fprintf(gnufile, "set title \"%s\"\n", title);
        fprintf(gnufile, "set xlabel \"%s\"\n", xlabel);
        fprintf(gnufile, "set ylabel \"%s\"\n", ylabel);
        fprintf(gnufile, "plot '%s' using 1:2 title \"Obs. spectrum\" with "
                         "lines lt -1, '%s' using 1:2 title "
                         "\"Corr. spectrum\" with lines lt 8\n",
                         tmpfile1, tmpfile2);
        fprintf(gnufile, "\n");

        /* Close temoprary GNUPLOT file */
        fclose(gnufile);

        /* Call GNUPLOT */
        sprintf(systemcall, "gnuplot -persist %s", tmpfile3);
        dummy = system(systemcall);

        /* Remove temporary GNUPLOT file */
        dummy = remove(tmpfile3);

    }

    /* Remove ASCII files with spectra and delete temporary directory */
    dummy = remove(tmpfile1);
    dummy = remove(tmpfile2);
    dummy = rmdir(tmpdir);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_trans_writeresults(mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Writes results of MOLECFIT and CALCTRANS into a file of the same format
     * as the input data file. In the case of FITS images, the transmission
     * curve and the telluric absorption corrected spectrum are written into
     * separate files.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_table *results = NULL;

    /* Read file with CALCTRANS results and convert it into CPL table */
    results = cpl_table_new(0);
    if ((status = mf_conv_readresults(results, drvpar)) !=
        CPL_ERROR_NONE) {
        cpl_table_delete(results);
        return status;
    }

    /* Write CPL table and CPL property list to file with format of input data
       file (FITS image: second file with transmission curve) */
    mf_conv_writefile(results, drvpar);

    /* Free allocated memory */
    cpl_table_delete(results);

    /* Return error code of last error */
    return cpl_error_get_code();
}

/**@}*/
