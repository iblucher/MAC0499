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
 * \file mf_mpfit.c
 *
 * Routines for handling CMPFIT
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  25 Jun 2010
 * \date   26 Jan 2014
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <mf_mpfit.h>


/*****************************************************************************
 *                                GLOBALS                                    *
 ****************************************************************************/

/* Definition of global variables */

/* Reference fit parameters */
extern cpl_array *reffitpar;

/* Declaration of global variables */

/* Number of MPFIT fitting function calls */
extern int nfev;
/* Number of LBLRTM calls (wavenumber-restricted subspectra are not
   counted individually) */
extern int n_code;
/* Total LBLRTM run time */
extern double t_code;


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code mf_mpfit(mp_result *result, cpl_table *spec, cpl_table *prof,
    cpl_table** prof_out, cpl_table** res_out, mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Handles the fitting routine CMPFIT (see mpfit.c for details on the
     * fitting algorithm). Needs the observed spectrum plus weights,
     * atmospheric profiles for pressure, temperature, and relevant molecules,
     * and fit-related parameters from the ::mfdrv driver parameter structure.
     * The input data are used to compute a model spectrum which is compared
     * to the observed spectrum by deriving a vector of weighted deviations.
     * CMPFIT optimises the model spectrum by manipulating a vector of the
     * fit parameters which is used as input for the model calculation.
     * Information on the fit quality is written into a special results
     * structure. The CPL table with the input spectrum is supplemented by the
     * best-fit model and the weighted deviations taken for the \f${\chi^2}\f$
     * computation. Moreover, an updated driver file, an ASCII file which
     * summarises the fit results (including the ppmv of the different
     * molecules and especially the water-related PWV in mm), and an ASCII
     * file with the best-fit atmospheric profiles are written.
     *
     * \b INPUT:
     * \param spec    CPL table with observed spectrum
     * \param prof    CPL table with atmospheric profiles
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    CPL table with observed and best-fit model spectrum
     * \param result  CMPFIT structure for fit results
     *
     * \b ERRORS:
     * - No data
     * - Insufficient memory
     * - Error in subroutine
     */

    cpl_error_code errstat = CPL_ERROR_NONE;
    cpl_parameter *pp;
    cpl_array *fitpars;
    mp_config config = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    mfpars fitpar;
    mfvars v;
    char errtxt[MF_MAXLEN];
    int m = 0, status = 0;
    double ts = 0., te = 0.;
    double fittime = 0.;

    /* Get number of data points */
    m = cpl_table_get_nrow(spec);

    /* Set parameter vector and parameter constraints structure */
    mf_mpfit_setpar(&fitpar, drvpar);
    if (fitpar.p == NULL || fitpar.pars == NULL) {
        result->status = -99;
        sprintf(errtxt, "%s: mfpars fitpar", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Set global reference parameter vector */
    reffitpar = cpl_array_wrap_double(fitpar.p, fitpar.n);

    /* Set fit precision */
    pp = cpl_parameterlist_find(drvpar->parlist, "ftol");
    config.ftol = cpl_parameter_get_double(pp);
    pp = cpl_parameterlist_find(drvpar->parlist, "xtol");
    config.xtol = cpl_parameter_get_double(pp);

    /* Set maximum number of iterations and function evaluations */
    config.maxiter = MF_MAXITER;
    config.maxfev = config.maxiter * fitpar.n;

    /* Pack all data and parameters into temporary structure in order to
       create a void pointer required by the CMPFIT user function */
    v.spec = spec;
    v.prof = prof;
    v.drvpar = drvpar;

    /* Set return of fit residuals and parameter errors */
    if ((int) mf_mpfit_allocmemresult(result, m, fitpar.n) ==
            (int) MF_ERROR_ISM) {
        return MF_ERROR_ISM;
    }

    /* Info messages */
    cpl_msg_info(cpl_func, "Fitting function calls and fit parameter "
                 "changes:");
    cpl_msg_info(cpl_func, "call par newval    reldev");

    /* Call fitting function for m data points and n parameters */
    ts = cpl_test_get_walltime();
    status = mpfit(mf_mpfit_calcdev, m, fitpar.n, fitpar.p, fitpar.pars,
                   &config, (void *) &v, result, prof_out);

    /* Get CMPFIT run time in min */
    te = cpl_test_get_walltime();
    fittime = (te - ts) / 60;

    /* Write best fit parameters to the MOLECFIT parameter structure */
    mf_mpfit_getpar(drvpar, &fitpar);

    /* Put fit parameters in CPL array */
    fitpars = cpl_array_wrap_double(fitpar.p, fitpar.n);

    /* Set flag for last modsim call to 1 if flag is active */
    pp = cpl_parameterlist_find(drvpar->parlist, "lastcall");
    if (cpl_parameter_get_int(pp) == 0) {
        cpl_parameter_set_int(pp, 1);
    }

    /* Fill CPL table spec with best-fit model spectrum and deviations and
       write best atmospheric profiles to ASCII file */
    errstat = mf_modsim(spec, prof_out, prof, drvpar, fitpars);

    /* Write a summary of the CMPFIT results into an ASCII file */
    mf_mpfit_writeresults(drvpar, spec, result, fittime, res_out);

    /* Free memory */
    mf_mpfit_freemempar(&fitpar);
    cpl_array_unwrap(reffitpar);
    reffitpar = NULL;
    cpl_array_unwrap(fitpars);

    if (status <= 0) {
        sprintf(errtxt, "%s: mpfit", MF_ERROR_EIS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_EIS, "%s", errtxt);
    } else if (errstat != CPL_ERROR_NONE) {
        return errstat;
    }

    return CPL_ERROR_NONE;
}


int mf_mpfit_calcdev(int m, int n, double *p, double *dy, double **dvec,
                     void *vars, cpl_table** prof_out)
{
    /*!
     * \callgraph
     *
     * User function for CMPFIT. Returns weighted deviations between the model
     * (calculated by ::mf_modsim) and the observed spectrum. The syntax of
     * the function is predefined by CMPFIT.
     *
     * \b INPUT:
     * \param m     number of data points
     * \param n     number of parameters
     * \param p     array of fit parameters
     * \param dvec  derivatives (not used)
     * \param vars  private data -> observed spectrum and driver file
     *                              parameters
     *
     * \b OUTPUT:
     * \param dy    array of residuals ([model - obs. spectrum] * weight)
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     */

    cpl_error_code errstat = CPL_ERROR_NONE;
    mfvars *v = (mfvars *) vars;
    mfdrv *drvpar;
    cpl_table *spec, *prof;
    cpl_array *fitpar;
    char errtxt[MF_MAXLEN];
    int j = 0, i = 0;
    double **junk;

    if ((junk = dvec)) {};

    /* Update number of fitting function call */
    nfev++;

    /* Unpack observed spectral data and input parameters */
    spec = v->spec;
    prof = v->prof;
    drvpar = v->drvpar;

    /* Check for parameters with "nan" as value */
    for (j = 0; j < n; j++) {
        if (isnan(p[j]) != 0) {
            sprintf(errtxt, "%s: double *p has nan values", MF_ERROR_IIP_TXT);
            cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
            break;
        }
    }

    /* Put fit parameters in CPL array */
    fitpar = cpl_array_wrap_double(p, n);

    /* Take reference parameter values in the case of errors */
    if (isnan(p[j-1]) != 0) {
        for (j = 0; j < n; j++) {
            cpl_array_set(fitpar, j, cpl_array_get(reffitpar, j, NULL));
        }
    }

    /* Model calculation if parameter values are valid */
    errstat = mf_modsim(spec, prof_out, prof, drvpar, fitpar);

    /* Fill array of residuals */
    for (i = 0; i < m; i++) {
        if (errstat != CPL_ERROR_NONE) {
            dy[i] = 0.;
        } else {
            dy[i] = cpl_table_get(spec, "dev", i, NULL);
        }
    }

    /* Free memory */
    cpl_array_unwrap(fitpar);

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return errstat;
}


cpl_error_code mf_mpfit_setpar(mfpars *fitpar, const mfdrv *drvpar)
{
    /*!
     * Provides a vector of parameters which are variables of the CMPFIT
     * fitting process. Moreover, constraints for these parameters are
     * delivered by the CMPFIT structure mp_par. Both objects are returned by
     * a container structure.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param fitpar  structure containing the fit parameters
     *
     * \b ERRORS:
     * - see ::mf_mpfit_allocmempar
     */

    cpl_error_code err = CPL_ERROR_NONE;
    cpl_parameter *pp;
    cpl_array *wlc_coef = NULL, *cont_coef = NULL;
    char name[MF_LENLINE+2];
    int npar = 0, nmolec = 0, nchip = 0, nwlc = 0, nrange = 0, ncont = 0;
    int trans = 1, i = 0, idx = -1, fitwlc1 = 0, j = 0, fitwlc = 0;
    int fitcont = 0;

    /* Get number of parameters */
    pp = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(pp);
    pp = cpl_parameterlist_find(drvpar->parlist, "nchip");
    nchip = cpl_parameter_get_int(pp);
    pp = cpl_parameterlist_find(drvpar->parlist, "wlc_n");
    nwlc = cpl_parameter_get_int(pp) + 1;
    pp = cpl_parameterlist_find(drvpar->parlist, "nrange");
    nrange = cpl_parameter_get_int(pp);
    pp = cpl_parameterlist_find(drvpar->parlist, "cont_n");
    ncont = cpl_parameter_get_int(pp) + 1;
    npar = nmolec + nchip * nwlc + nrange * ncont + 3;

    /* Transmission or emission? */
    pp = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(pp);
    if (trans == 0) {
        npar++;
    }

    /* Allocate memory for mfpars structure */
    if ((err = mf_mpfit_allocmempar(fitpar, npar)) != CPL_ERROR_NONE) {
        return err;
    }

    /* Set parameter values and constraints for fitting */

    /* Amount of molecular gas */
    for (i = 0; i < nmolec; i++) {
        idx++;
        fitpar->p[idx] = cpl_table_get(drvpar->molectab, "relcol", i, NULL);
        fitpar->pars[idx].fixed = 1 - cpl_table_get(drvpar->molectab,
                                                    "fit_molec", i, NULL);
        fitpar->pars[idx].limited[0] = 1;
        fitpar->pars[idx].limited[1] = 1;
        fitpar->pars[idx].limits[0] = 1e-5;
        fitpar->pars[idx].limits[1] = 1e2;
        fitpar->pars[idx].relstep = 0.01;
        strcpy(name, cpl_table_get_string(drvpar->molectab, "list_molec", i));
        strcpy(fitpar->pars[idx].parname, name);
    }

    /* Wavelength solution */
    pp = cpl_parameterlist_find(drvpar->parlist, "fit_wlc_lin");
    fitwlc1 = cpl_parameter_get_int(pp);
    for (j = 0; j < nchip; j++) {
        wlc_coef = cpl_array_duplicate(cpl_table_get_array(drvpar->chiptab,
                                                           "wlc_coef", j));
        fitwlc = cpl_table_get(drvpar->chiptab, "fit_chip", j, NULL);
        for (i = 0; i < nwlc; i++) {
            idx++;
            fitpar->p[idx] = cpl_array_get_double(wlc_coef, i, NULL);
            if (i == 1 && fitwlc == 1) {
                fitpar->pars[idx].fixed = 1 - fitwlc1;
            } else {
                fitpar->pars[idx].fixed = 1 - fitwlc;
            }
            fitpar->pars[idx].limited[0] = 0;
            fitpar->pars[idx].limited[1] = 0;
            fitpar->pars[idx].limits[0] = 0.;
            fitpar->pars[idx].limits[1] = 0.;
            if (i == 1) {
                fitpar->pars[idx].relstep = 1e-4;
            } else {
                fitpar->pars[idx].relstep = 0.1;
            }
            sprintf(name, "c%d_wlc_coef_%d", j+1, i);
            strcpy(fitpar->pars[idx].parname, name);
        }
        cpl_array_delete(wlc_coef);
    }

    /* Continuum */
    for (j = 0; j < nrange; j++) {
        cont_coef = cpl_array_duplicate(cpl_table_get_array(drvpar->rangetab,
                                                            "cont_coef", j));
        fitcont = cpl_table_get(drvpar->rangetab, "fit_range", j, NULL);
        for (i = 0; i < ncont; i++) {
            idx++;
            fitpar->p[idx] = cpl_array_get_double(cont_coef, i, NULL);
            fitpar->pars[idx].fixed = 1 - fitcont;
            if (i == 0) {
                fitpar->pars[idx].limited[0] = 1;
            } else {
                fitpar->pars[idx].limited[0] = 0;
            }
            fitpar->pars[idx].limited[1] = 0;
            fitpar->pars[idx].limits[0] = MF_TOL * fitpar->p[idx];
            fitpar->pars[idx].limits[1] = 0.;
            fitpar->pars[idx].relstep = 0.01;
            sprintf(name, "c%d_cont_coef_%d", j+1, i);
            strcpy(fitpar->pars[idx].parname, name);
        }
        cpl_array_delete(cont_coef);
    }

    /* Resolution */

    /* Boxcar */
    idx++;
    pp = cpl_parameterlist_find(drvpar->parlist, "relres_box");
    fitpar->p[idx] = cpl_parameter_get_double(pp);
    pp = cpl_parameterlist_find(drvpar->parlist, "fit_res_box");
    fitpar->pars[idx].fixed = 1 - cpl_parameter_get_int(pp);
    fitpar->pars[idx].limited[0] = 1;
    fitpar->pars[idx].limited[1] = 1;
    fitpar->pars[idx].limits[0] = 0.;
    fitpar->pars[idx].limits[1] = 2.;
    fitpar->pars[idx].relstep = 0.1;
    strcpy(name, "relres_box");
    strcpy(fitpar->pars[idx].parname, name);

    /* Gaussian */
    idx++;
    pp = cpl_parameterlist_find(drvpar->parlist, "res_gauss");
    fitpar->p[idx] = cpl_parameter_get_double(pp);
    pp = cpl_parameterlist_find(drvpar->parlist, "fit_res_gauss");
    fitpar->pars[idx].fixed = 1 - cpl_parameter_get_int(pp);
    fitpar->pars[idx].limited[0] = 1;
    fitpar->pars[idx].limited[1] = 1;
    fitpar->pars[idx].limits[0] = 0.;
    fitpar->pars[idx].limits[1] = 100.;
    fitpar->pars[idx].relstep = 0.1;
    strcpy(name, "res_gauss");
    strcpy(fitpar->pars[idx].parname, name);

    /* Lorentzian */
    idx++;
    pp = cpl_parameterlist_find(drvpar->parlist, "res_lorentz");
    fitpar->p[idx] = cpl_parameter_get_double(pp);
    pp = cpl_parameterlist_find(drvpar->parlist, "fit_res_lorentz");
    fitpar->pars[idx].fixed = 1 - cpl_parameter_get_int(pp);
    fitpar->pars[idx].limited[0] = 1;
    fitpar->pars[idx].limited[1] = 1;
    fitpar->pars[idx].limits[0] = 0.;
    fitpar->pars[idx].limits[1] = 100.;
    fitpar->pars[idx].relstep = 0.1;
    strcpy(name, "res_lorentz");
    strcpy(fitpar->pars[idx].parname, name);

    /* Telescope background scaling */
    if (trans == 0) {
        idx++;
        pp = cpl_parameterlist_find(drvpar->parlist, "telback");
        fitpar->p[idx] = cpl_parameter_get_double(pp);
        pp = cpl_parameterlist_find(drvpar->parlist, "fit_back");
        fitpar->pars[idx].fixed = 1 - cpl_parameter_get_int(pp);
        fitpar->pars[idx].limited[0] = 1;
        fitpar->pars[idx].limited[1] = 1;
        fitpar->pars[idx].limits[0] = 0.;
        fitpar->pars[idx].limits[1] = 1. - MF_TOL;
        fitpar->pars[idx].relstep = 0.01;
        strcpy(name, "telback");
        strcpy(fitpar->pars[idx].parname, name);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_mpfit_getpar(mfdrv *drvpar, const mfpars *fitpar)
{
    /*!
     * Writes fit parameter values as provided by an ::mfpars structure into
     * the ::mfdrv parameter structure.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     * \param fitpar  structure containing the fit parameters
     *
     * \b OUTPUT:
     * \param drvpar  updated ::mfdrv parameter structure
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    cpl_parameter *pp;
    cpl_array *wlc_coef = NULL, *cont_coef = NULL;
    char errtxt[MF_MAXLEN];
    int npar = 0, nmolec = 0, nchip = 0, nwlc = 0, nrange = 0, ncont = 0;
    int trans = 1, i = 0, idx = -1, j = 0;

    /* Get number of parameters */
    pp = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(pp);
    pp = cpl_parameterlist_find(drvpar->parlist, "nchip");
    nchip = cpl_parameter_get_int(pp);
    pp = cpl_parameterlist_find(drvpar->parlist, "wlc_n");
    nwlc = cpl_parameter_get_int(pp) + 1;
    pp = cpl_parameterlist_find(drvpar->parlist, "nrange");
    nrange = cpl_parameter_get_int(pp);
    pp = cpl_parameterlist_find(drvpar->parlist, "cont_n");
    ncont = cpl_parameter_get_int(pp) + 1;
    npar = nmolec + nchip * nwlc + nrange * ncont + 3;

    /* Transmission or emission? */
    pp = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(pp);
    if (trans == 0) {
        npar++;
    }

    /* Correct number of parameters? */
    if (fitpar->n != npar) {
        sprintf(errtxt, "%s: fitpar->n != npar derived from mfdrv *drvpar",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Write fit parameter values into mfdrv structure */

    /* Amount of molecular gas */
    for (i = 0; i < nmolec; i++) {
        idx++;
        cpl_table_set_double(drvpar->molectab, "relcol", i, fitpar->p[idx]);
    }

    /* Wavelength solution */
    wlc_coef = cpl_array_new(nwlc, CPL_TYPE_DOUBLE);
    for (j = 0; j < nchip; j++) {
        for (i = 0; i < nwlc; i++) {
            idx++;
            cpl_array_set_double(wlc_coef, i, fitpar->p[idx]);
        }
        cpl_table_set_array(drvpar->chiptab, "wlc_coef", j, wlc_coef);
    }
    cpl_array_delete(wlc_coef);

    /* Continuum */
    cont_coef = cpl_array_new(ncont, CPL_TYPE_DOUBLE);
    for (j = 0; j < nrange; j++) {
        for (i = 0; i < ncont; i++) {
            idx++;
            cpl_array_set_double(cont_coef, i, fitpar->p[idx]);
        }
        cpl_table_set_array(drvpar->rangetab, "cont_coef", j, cont_coef);
    }
    cpl_array_delete(cont_coef);

    /* Resolution */
    idx++;
    pp = cpl_parameterlist_find(drvpar->parlist, "relres_box");
    cpl_parameter_set_double(pp, fitpar->p[idx]);
    idx++;
    pp = cpl_parameterlist_find(drvpar->parlist, "res_gauss");
    cpl_parameter_set_double(pp, fitpar->p[idx]);
    idx++;
    pp = cpl_parameterlist_find(drvpar->parlist, "res_lorentz");
    cpl_parameter_set_double(pp, fitpar->p[idx]);

    /* Telescope background scaling */
    if (trans == 0) {
        idx++;
        pp = cpl_parameterlist_find(drvpar->parlist, "telback");
        cpl_parameter_set_double(pp, fitpar->p[idx]);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_mpfit_allocmempar(mfpars *fitpar, const int npar)
{
    /*!
     * Allocates memory for an ::mfpars structure.
     *
     * \b INPUT:
     * \param npar    number of fit parameters
     *
     * \b OUTPUT:
     * \param fitpar  ::mfpars structure containing npar fit parameters
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     * - Insufficient memory
     */

    cpl_boolean fl_mem = CPL_TRUE;
    char errtxt[MF_MAXLEN];
    int it = 0, i = 0, nchar = MF_LENLINE+2;

    /* Check number of parameters */
    fitpar->n = npar;
    if (fitpar->n < 0) {
        fitpar->n = 0;
        sprintf(errtxt, "%s: npar < 0", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    } else if (fitpar->n == 0) {
        return CPL_ERROR_NONE;
    }

    /* Allocate memory for parameter vector */
    fitpar->p = (double *) calloc(fitpar->n, sizeof(double));
    if (fitpar->p == NULL) {
        fitpar->n = 0;
        sprintf(errtxt, "%s: mfpars *fitpar", MF_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_ISM, "%s", errtxt);
    }

    /* Allocate memory for parameter constraints structure */
    fitpar->pars = (mp_par *) calloc(fitpar->n, sizeof(mp_par));
    if (fitpar->pars == NULL) {
        fitpar->n = 0;
        free(fitpar->p);
        fitpar->p = NULL;
        sprintf(errtxt, "%s: mfpars *fitpar", MF_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_ISM, "%s", errtxt);
    }

    /* Allocate memory for parameter names in parameter constraints
       structure */
    for (it = 0; it < 2; it++) {
        for (i = 0; i < fitpar->n; i++) {
            fitpar->pars[i].parname = (char *) calloc(nchar, sizeof(char));
            if (it == 0 && fitpar->pars[i].parname == NULL) {
                nchar = 0;
                fl_mem = CPL_FALSE;
                continue;
            }
        }
        if (it == 0 && fl_mem == CPL_TRUE) {
            break;
        } else if (it == 1 && fl_mem == CPL_FALSE) {
            fitpar->n = 0;
            free(fitpar->p);
            fitpar->p = NULL;
            free(fitpar->pars);
            fitpar->pars = NULL;
            sprintf(errtxt, "%s: mfpars *fitpar", MF_ERROR_ISM_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_ISM, "%s",
                                         errtxt);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_mpfit_freemempar(mfpars *fitpar)
{
    /*!
     * Frees memory occupied by an ::mfpars structure.
     *
     * \b INPUT:
     * \param fitpar  ::mfpars structure containing npar fit parameters
     *
     * \b OUTPUT:
     * \param fitpar  ::mfpars structure without allocated memory
     *
     * \b ERRORS:
     * - none
     */

    int i;

    /* Free memory occupied by parameter vector */
    if (fitpar->p != NULL) {
        free(fitpar->p);
        fitpar->p = NULL;
    }

    /* Free memory occupied by parameter constraints structure */
    if (fitpar->pars != NULL) {
        for (i = 0; i < fitpar->n; i++) {
            free(fitpar->pars[i].parname);
            fitpar->pars[i].parname = NULL;
        }
        free(fitpar->pars);
        fitpar->pars = NULL;
    }

    /* Set number of parameters to 0 */
    fitpar->n = 0;

    return CPL_ERROR_NONE;
}


cpl_error_code mf_mpfit_allocmemresult(mp_result *result, const int m,
                                       const int n)
{
    /*!
     * Allocates memory for fit residuals and parameter errors in the CMPFIT
     * results structure.
     *
     * \b INPUT:
     * \param m       number of data points
     * \param n       number of parameters
     *
     * \b OUTPUT:
     * \param result  CMPFIT structure for fit results with allocated memory
     *                for residuals and parameter errors
     *
     * \b ERRORS:
     * - Insufficient memory
     */

    char errtxt[MF_MAXLEN];

    /* No consideration of the covariance matrix */
    result->covar = NULL;

    /* Memory allocation for fit residuals */
    result->resid = (double *) calloc(m, sizeof(double));
    if (result->resid == NULL) {
        result->status = -99;
        result->xerror = NULL;
        sprintf(errtxt, "%s: mp_result *result", MF_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_ISM, "%s", errtxt);
    }

    /* Memory allocation for parameter errors */
    result->xerror = (double *) calloc(n, sizeof(double));
    if (result->xerror == NULL) {
        result->status = -99;
        free(result->resid);
        result->resid = NULL;
        sprintf(errtxt, "%s: mp_result *result", MF_ERROR_ISM_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_ISM, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_mpfit_freememresult(mp_result *result)
{
    /*!
     * Frees memory occupied by a CMPFIT structure for fit residuals and
     * parameter errors. Memory allocated for carrying the covariance matrix
     * is not freed.
     *
     * \b INPUT:
     * \param result  CMPFIT structure for fit results with allocated memory
     *                for residuals and parameter errors
     *
     * \b OUTPUT:
     * \param result  CMPFIT structure for fit results without memory for
     *                residuals and parameter errors
     *
     * \b ERRORS:
     * - none
     */

    /* Free memory for residuals */
    if (result->resid != NULL) {
        free(result->resid);
        result->resid = NULL;
    }

    /* Free memory for parameter errors */
    if (result->xerror != NULL) {
        free(result->xerror);
        result->xerror = NULL;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_mpfit_writeresults(mfdrv *drvpar, const cpl_table *spec,
                                     const mp_result *result,
                                     const double fittime,
                                     cpl_table** outtable)
{
    /*!
     * Writes a summary of the CMPFIT results into an ASCII file.
     *
     * \b INPUT:
     * \param drvpar   ::mfdrv parameter structure
     * \param spec     CPL table with observed and best-fit spectrum
     * \param result   CMPFIT structure for fit results
     * \param fittime  fit run time in minutes
     *
     * \b OUTPUT:
     * \param drvpar   ::mfdrv parameter structure with \f${\chi^2}\f$
     *
     *
     * \b ERRORS:
     * - File opening failed
     * - Insufficient data points
     * - Invalid object value(s)
     */

    FILE *stream;
    cpl_error_code err = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_array *wlc_coef = NULL, *cont_coef = NULL;
    char output_dir[MF_MAXLEN], basedir[MF_MAXLEN];
    char output_name[MF_MAXLEN], outfile[MF_MAXLEN];
    char errtxt[MF_MAXLEN];
    char **mol;
    int l = 0, nmr = 0, nmw = 0, nw = 0, nmolec = 0, nchip = 0, nwlc = 0;
    int nrange = 0, ncont = 0, npar = 0, trans = 0, ipar = 0, j = 0, fit = 0;
    int i = 0, range = 0, chip = 0, fith2o = 0;
    const int *mrange;
    int *fitmol;
    double chi2 = 0., chi2red = 0., wsum = 0., wfsum = 0., w2sum = 0.;
    double rms = 0., slitw = 0., pixsc = 0., h2ocol = 0., reldel = 0.;
    double reldelh2o = 0.;
    const double *weight, *mweight, *dev, *flux;
    double *relcol, *ppmv;

    /* Get output file name */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(output_dir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strcpy(output_name, cpl_parameter_get_string(p));
    sprintf(outfile, "%s%s_fit.res", output_dir, output_name);

    p = cpl_parameterlist_find(drvpar->parlist, "nchip");
    nchip = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "wlc_n");
    nwlc = cpl_parameter_get_int(p) + 1;

    /* Get number of parameters */
    p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    nrange = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "cont_n");
    ncont = cpl_parameter_get_int(p) + 1;
    npar = nmolec + nchip * nwlc + nrange * ncont + 3;

    // For the table //
    char outfitsfile[MF_MAXLEN];
    sprintf(outfitsfile, "%s%s_fit.res.fits", output_dir, output_name);
    unsigned nRows = 1 +               // H2O MM
                     nmolec +          // RELATIVE MOLECULAR GAS COLUMNS
                     nmolec +          // MOLECULAR GAS COLUMNS IN PPMV
                     1 +               // telback
                     nrange * ncont +  // CONTINUUM CORRECTION
                     nchip * nwlc +    // WAVELENGTH SOLUTION
                     1 +               // FWHM Lorentzian
                     1 +               // FWHM Gaussian
                     2 +               // FWHM boxcar + in pixels
                     15;               // MPFIT
    
    // If not the first call, rewrite it.
    // Otherwise, we end up with a memory leak.
    if (*outtable) {
        cpl_table_delete(*outtable);
    }
    *outtable = cpl_table_new(nRows);
    cpl_table_new_column(*outtable, "parameter", CPL_TYPE_STRING);
    cpl_table_new_column(*outtable, "value", CPL_TYPE_DOUBLE);
    cpl_table_new_column(*outtable, "uncertainty", CPL_TYPE_DOUBLE);
    unsigned outIndex = 0;

    /* Open output file */
    if ((stream = fopen(outfile, "w+")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, outfile);
        cpl_table_delete(*outtable);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Write spectrum name */
    fprintf(stream, "DATA FILE:\n");
    p = cpl_parameterlist_find(drvpar->parlist, "filename");
    fprintf(stream, "%s\n\n", cpl_parameter_get_string(p));

    /* Write status ID and return in the case of errors */
    fprintf(stream, "MPFIT RESULTS:\n");
    if (result->status <= 0) {
        fprintf(stream, "Status: %d -> ERROR!\n\n", result->status);
        fclose(stream);
        return CPL_ERROR_NONE;
    } else {
        fprintf(stream, "Status:                    %d\n", result->status);
        cpl_table_set_string(*outtable, "parameter", outIndex, "status");
        cpl_table_set_double(*outtable, "value", outIndex, result->status);
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
        ++outIndex;
    }

    /* Get number of data points with valid model fluxes and non-zero weight
       (= DOF+1) and calculate chi^2 */
    weight = cpl_table_get_data_double_const(spec, "weight");
    mrange = cpl_table_get_data_int_const(spec, "mrange");
    mweight = cpl_table_get_data_double_const(spec, "mweight");
    dev = cpl_table_get_data_double_const(spec, "dev");

    cpl_size nrows = cpl_table_get_nrow(spec);
    for (chi2 = 0, l = 0; l < result->nfunc && l < nrows; l++) {
        if (mrange[l] > 0) {
            nmr++;
        }
        if (mweight[l] > 0) {
            nmw++;
            if (weight[l] > 0) {
                nw++;
                chi2 += dev[l] * dev[l];
            }
        }
    }
    if (nw == 1) {
        sprintf(errtxt, "%s: cpl_table *spec "
                "(only 1 data point with weight > 0)", MF_ERROR_ISD_TXT);
        err = cpl_error_set_message(cpl_func, MF_ERROR_ISD, "%s", errtxt);
    }

    /* Write chi^2 into parameter list */
    p = cpl_parameterlist_find(drvpar->parlist, "chi2");
    cpl_parameter_set_double(p, chi2);

    /* Write CMPFIT results */
    fprintf(stream, "Fit parameters:            %d\n", result->npar);
    fprintf(stream, "Data points:               %d\n", result->nfunc);
    fprintf(stream, "Weight > 0:                %d\n", nw);
    fprintf(stream, "Frac. of valid model pix.: %.2f\n",
            ((double) nmw) / ((double) nmr));
    fprintf(stream, "Iterations:                %d\n", result->niter);
    fprintf(stream, "Function evaluations:      %d\n", result->nfev);
    fprintf(stream, "Fit run time in min:       %.2f\n", fittime);
    fprintf(stream, "Avg. LBLRTM run time in s: %.2lf\n",
            t_code / ((double) n_code));
    fprintf(stream, "LBLRTM calls:              %i\n", n_code);
    fprintf(stream, "Initial chi2:              %.3e\n", result->orignorm);
    fprintf(stream, "Best chi2:                 %.3e\n", chi2);
    if (nw <= 1) {
        fprintf(stream, "Reduced chi2:              UNDEF\n");
        fprintf(stream, "RMS rel. to error:         UNDEF\n");
    } else {
        chi2red = chi2 / (nw - 1);
        fprintf(stream, "Reduced chi2:              %.3e\n", chi2red);
        fprintf(stream, "RMS rel. to error:         %.3e\n", sqrt(chi2red));
    }

    cpl_table_set_string(*outtable, "parameter", outIndex, "fit_params");
    cpl_table_set_double(*outtable, "value", outIndex, result->npar);
    cpl_table_set_double(*outtable, "uncertainty", outIndex++, -1);

    cpl_table_set_string(*outtable, "parameter", outIndex, "data_points");
    cpl_table_set_double(*outtable, "value", outIndex, result->nfunc);
    cpl_table_set_double(*outtable, "uncertainty", outIndex++, -1);

    cpl_table_set_string(*outtable, "parameter", outIndex, "positive_weights");
    cpl_table_set_double(*outtable, "value", outIndex, nw);
    cpl_table_set_double(*outtable, "uncertainty", outIndex++, -1);

    cpl_table_set_string(*outtable, "parameter", outIndex, "valid_pix_frac");
    cpl_table_set_double(*outtable, "value", outIndex, ((double) nmw) / ((double) nmr));
    cpl_table_set_double(*outtable, "uncertainty", outIndex++, -1);

    cpl_table_set_string(*outtable, "parameter", outIndex, "iterations");
    cpl_table_set_double(*outtable, "value", outIndex, result->niter);
    cpl_table_set_double(*outtable, "uncertainty", outIndex++, -1);

    cpl_table_set_string(*outtable, "parameter", outIndex, "func_eval");
    cpl_table_set_double(*outtable, "value", outIndex, result->nfev);
    cpl_table_set_double(*outtable, "uncertainty", outIndex++, -1);

    cpl_table_set_string(*outtable, "parameter", outIndex, "lblrtm_calls");
    cpl_table_set_double(*outtable, "value", outIndex, n_code);
    cpl_table_set_double(*outtable, "uncertainty", outIndex++, -1);

    cpl_table_set_string(*outtable, "parameter", outIndex, "initial_chi2");
    cpl_table_set_double(*outtable, "value", outIndex, result->orignorm);
    cpl_table_set_double(*outtable, "uncertainty", outIndex++, -1);

    cpl_table_set_string(*outtable, "parameter", outIndex, "best_chi2");
    cpl_table_set_double(*outtable, "value", outIndex, chi2);
    cpl_table_set_double(*outtable, "uncertainty", outIndex++, -1);

    cpl_table_set_string(*outtable, "parameter", outIndex, "reduced_chi2");
    cpl_table_set_string(*outtable, "parameter", outIndex + 1, "rms_rel_to_err");

    if (nw <= 1) {

        cpl_table_set_double(*outtable, "value", outIndex, -1);
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);

        cpl_table_set_double(*outtable, "value", outIndex + 1, -1);
        cpl_table_set_double(*outtable, "uncertainty", outIndex + 1, -1);

    } else {

        cpl_table_set_double(*outtable, "value", outIndex, chi2red);
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);

        cpl_table_set_double(*outtable, "value", outIndex + 1, sqrt(chi2red));
        cpl_table_set_double(*outtable, "uncertainty", outIndex + 1, -1);
    }

    outIndex += 2;

    /* Compute RMS relative to weighted mean and write it to file */
    flux = cpl_table_get_data_double_const(spec, "flux");
    for (l = 0; l < result->nfunc && l < nrows; l++) {
        if (mweight[l] > 0) {
            wsum += weight[l];
            wfsum += weight[l] * flux[l];
            w2sum += weight[l] * weight[l];
        }
    }
    if (wsum == 0) {
        fprintf(stream, "RMS rel. to mean:          UNDEF\n\n");
        sprintf(errtxt, "%s: cpl_table *spec (all weights = 0)",
                MF_ERROR_IOV_TXT);
        err = cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    } else if (wfsum == 0) {
        fprintf(stream, "RMS rel. to mean:          UNDEF\n\n");
        sprintf(errtxt, "%s: cpl_table *spec (all fluxes = 0)",
                MF_ERROR_IOV_TXT);
        err = cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    } else {
        rms = sqrt(chi2 / w2sum) * (wsum / wfsum);
        fprintf(stream, "RMS rel. to mean:          %.3e\n\n", rms);
    }

    cpl_table_set_string(*outtable, "parameter", outIndex, "rms_rel_to_mean");
    if (wsum == 0) {
        cpl_table_set_double(*outtable, "value", outIndex, -1);
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
    }
    else if (wfsum == 0) {
        cpl_table_set_double(*outtable, "value", outIndex, -1);
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
    }
    else {
        cpl_table_set_double(*outtable, "value", outIndex, rms);
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
    }
    ++outIndex;

    /* Molecular spectrum in emission? */
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);

    /* Write best-fit parameters */
    fprintf(stream, "BEST-FIT PARAMETERS:\n\n");

    /* Spectral resolution */
    fprintf(stream, "SPECTRAL RESOLUTION:\n");
    p = cpl_parameterlist_find(drvpar->parlist, "fit_res_box");
    fit = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "slitw");
    slitw = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find(drvpar->parlist, "pixsc");
    pixsc = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find(drvpar->parlist, "relres_box");
    if (fit == 1) {
        fprintf(stream, "Rel. FWHM of boxcar (slit width = 1): %.3f "
                "+- %.3f\n", cpl_parameter_get_double(p),
                result->xerror[npar-3]);
        fprintf(stream, "FWHM of boxcar in pixels:             %.3f "
                "+- %.3f\n", cpl_parameter_get_double(p) * slitw / pixsc,
                result->xerror[npar-3] * slitw / pixsc);
    } else {
        fprintf(stream, "Rel. FWHM of boxcar (slit width = 1): %.3f\n",
                cpl_parameter_get_double(p));
        fprintf(stream, "FWHM of boxcar in pixels:             %.3f\n",
                cpl_parameter_get_double(p) * slitw / pixsc);
    }
    cpl_table_set_string(*outtable, "parameter", outIndex, "boxfwhm");
    cpl_table_set_string(*outtable, "parameter", outIndex + 1, "boxfwhm_pix");
    cpl_table_set_double(*outtable, "value", outIndex, cpl_parameter_get_double(p));
    cpl_table_set_double(*outtable, "value", outIndex + 1, cpl_parameter_get_double(p) * slitw / pixsc);
    if (fit == 1) {
        cpl_table_set_double(*outtable, "uncertainty", outIndex, result->xerror[npar-3]);
        cpl_table_set_double(*outtable, "uncertainty", outIndex + 1, result->xerror[npar-3] * slitw / pixsc);
    }
    else {
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
        cpl_table_set_double(*outtable, "uncertainty", outIndex + 1, -1);
    }
    outIndex += 2;

    p = cpl_parameterlist_find(drvpar->parlist, "fit_res_gauss");
    fit = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "res_gauss");
    if (fit == 1) {
        fprintf(stream, "FWHM of Gaussian in pixels:           %.3f "
                "+- %.3f\n", cpl_parameter_get_double(p),
                result->xerror[npar-2]);
    } else {
        fprintf(stream, "FWHM of Gaussian in pixels:           %.3f\n",
                cpl_parameter_get_double(p));
    }
    // For the table //
    cpl_table_set_string(*outtable, "parameter", outIndex, "gaussfwhm");
    cpl_table_set_double(*outtable, "value", outIndex, cpl_parameter_get_double(p));
    if (fit == 1) {
         cpl_table_set_double(*outtable, "uncertainty", outIndex, result->xerror[npar-2]);
    }
    else {
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
    }
    ++outIndex;

    p = cpl_parameterlist_find(drvpar->parlist, "fit_res_lorentz");
    fit = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "res_lorentz");
    if (fit == 1) {
        fprintf(stream, "FWHM of Lorentzian in pixels:         %.3f "
                "+- %.3f\n\n", cpl_parameter_get_double(p),
                result->xerror[npar-1]);
    } else {
        fprintf(stream, "FWHM of Lorentzian in pixels:         %.3f\n\n",
                cpl_parameter_get_double(p));
    }
    cpl_table_set_string(*outtable, "parameter", outIndex, "lorentzfwhm");
    cpl_table_set_double(*outtable, "value", outIndex, cpl_parameter_get_double(p));
    if (fit == 1) {
         cpl_table_set_double(*outtable, "uncertainty", outIndex, result->xerror[npar-1]);
    }
    else {
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
    }
    ++outIndex;

    /* Wavelength solution */
    fprintf(stream, "WAVELENGTH SOLUTION:\n");
    ipar = nmolec;
    for (j = 0; j < nchip; j++) {
        wlc_coef = cpl_array_duplicate(cpl_table_get_array(drvpar->chiptab,
                                                           "wlc_coef", j));
        fit = cpl_table_get(drvpar->chiptab, "fit_chip", j, NULL);
        for (i = 0; i < nwlc; i++) {
            if (fit == 1) {
                fprintf(stream, "Chip %d, coef %d: %10.3e +- %.3e\n", j+1, i,
                        cpl_array_get_double(wlc_coef, i, NULL),
                        result->xerror[ipar]);
            } else {
                fprintf(stream, "Chip %d, coef %d: %10.3e\n", j+1, i,
                        cpl_array_get_double(wlc_coef, i, NULL));
            }
            char pname[MF_MAXLEN];
            sprintf(pname, "Chip %d, coef %d", j+1, i);
            cpl_table_set_string(*outtable, "parameter", outIndex, pname);
            cpl_table_set_double(*outtable, "value", outIndex, cpl_array_get_double(wlc_coef, i, NULL));
            if (fit == 1) {
                 cpl_table_set_double(*outtable, "uncertainty", outIndex, result->xerror[ipar]);
            }
            else {
                cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
            }
            ++outIndex;
            ipar++;
        }
        cpl_array_delete(wlc_coef);
    }
    fprintf(stream, "\n");

    /* Continuum correction */
    fprintf(stream, "CONTINUUM CORRECTION:\n");
    for (j = 0; j < nrange; j++) {
        range = j + 1;
        chip = cpl_table_get(drvpar->rangetab, "chip", j, NULL);
        cont_coef = cpl_array_duplicate(cpl_table_get_array(drvpar->rangetab,
                                                            "cont_coef", j));
        fit = cpl_table_get(drvpar->rangetab, "fit_range", j, NULL);
        for (i = 0; i < ncont; i++) {
            if (fit == 1) {
                fprintf(stream, "Range %d, chip %d, coef %d: %10.3e "
                        "+- %.3e\n", range, chip, i,
                        cpl_array_get_double(cont_coef, i, NULL),
                        result->xerror[ipar]);
            } else {
                fprintf(stream, "Range %d, chip %d, coef %d: %10.3e\n",
                        range, chip, i,
                        cpl_array_get_double(cont_coef, i, NULL));
            }
            char pname[MF_MAXLEN];
            sprintf(pname, "Range %d, chip %d, coef %d", range, chip, i);
            cpl_table_set_string(*outtable, "parameter", outIndex, pname);
            cpl_table_set_double(*outtable, "value", outIndex, cpl_array_get_double(cont_coef, i, NULL));
            if (fit == 1) {
                 cpl_table_set_double(*outtable, "uncertainty", outIndex, result->xerror[ipar]);
            }
            else {
                cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
            }
            ++outIndex;
            ipar++;
        }
        cpl_array_delete(cont_coef);
    }
    fprintf(stream, "\n");

    /* Emissivity (for radiance spectrum only) */
    cpl_table_set_string(*outtable, "parameter", outIndex, "telback");
    if (trans == 0) {
        npar++;
        p = cpl_parameterlist_find(drvpar->parlist, "fit_back");
        fit = cpl_parameter_get_int(p);
        p = cpl_parameterlist_find(drvpar->parlist, "telback");
        if (fit == 1) {
            fprintf(stream, "EMISSIVITY: %f +- %f\n\n",
                    cpl_parameter_get_double(p), result->xerror[npar-1]);
        } else {
            fprintf(stream, "EMISSIVITY: %f\n\n",
                    cpl_parameter_get_double(p));
        }
        cpl_table_set_double(*outtable, "value", outIndex, cpl_parameter_get_double(p));
        if (fit == 1) {
            cpl_table_set_double(*outtable, "uncertainty", outIndex, result->xerror[npar-1]);
        }
        else {
            cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
        }
    }
    else {
        cpl_table_set_double(*outtable, "value", outIndex, -1);
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
    }
    ++outIndex;

    /* Relative molecular gas columns */
    fprintf(stream, "RELATIVE MOLECULAR GAS COLUMNS:\n");
    mol = cpl_table_get_data_string(drvpar->molectab, "list_molec");
    fitmol = cpl_table_get_data_int(drvpar->molectab, "fit_molec");
    relcol = cpl_table_get_data_double(drvpar->molectab, "relcol");
    for (i = 0; i < nmolec; i++) {
        char pname[MF_MAXLEN];
        sprintf(pname, "rel_mol_col_%s", mol[i]);
        cpl_table_set_string(*outtable, "parameter", outIndex, pname);
        cpl_table_set_double(*outtable, "value", outIndex, relcol[i]);
        if (fitmol[i] == 1) {
            fprintf(stream, "%3s: %.3f +- %.3f\n", mol[i], relcol[i],
                    result->xerror[i]);
        } else {
            fprintf(stream, "%3s: %.3f\n", mol[i], relcol[i]);
        }
        if (fitmol[i] == 1) {
            cpl_table_set_double(*outtable, "uncertainty", outIndex, result->xerror[i]);
        }
        else {
            cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
        }
        ++outIndex;
    }
    fprintf(stream, "\n");

    /* Molecular gas columns in ppmv */
    fprintf(stream, "MOLECULAR GAS COLUMNS IN PPMV:\n");
    mf_atm_calcol(&h2ocol, drvpar);
    ppmv = cpl_table_get_data_double(drvpar->molectab, "ppmv");
    cpl_size nrows_molectab = cpl_table_get_nrow(drvpar->molectab);
    for (i = 0; i < nmolec && i < nrows_molectab && ppmv; i++) {
        reldel = result->xerror[i] / relcol[i];
        char pname[MF_MAXLEN];
        sprintf(pname, "rel_mol_col_ppmv_%s", mol[i]);
        cpl_table_set_string(*outtable, "parameter", outIndex, pname);
        cpl_table_set_double(*outtable, "value", outIndex, ppmv[i]);
        if (strncmp(mol[i], "H2O", 3) == 0) {
            fith2o = fitmol[i];
            reldelh2o = reldel;
        }
        if (fitmol[i] == 1) {
            fprintf(stream, "%3s: %.3e +- %.3e\n", mol[i], ppmv[i],
                    reldel * ppmv[i]);
        } else {
            fprintf(stream, "%3s: %.3e\n", mol[i], ppmv[i]);
        }
        if (fitmol[i] == 1) {
            cpl_table_set_double(*outtable, "uncertainty", outIndex, reldel * ppmv[i]);
        }
        else {
            cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
        }
        ++outIndex;
    }
    fprintf(stream, "\n");

    /* H2O column in mm */
    if (fith2o == 1) {
        fprintf(stream, "H2O COLUMN IN MM: %.3f +- %.3f\n\n", h2ocol,
                reldelh2o * h2ocol);
    } else {
        fprintf(stream, "H2O COLUMN IN MM: %.3f\n\n", h2ocol);
    }
    cpl_table_set_string(*outtable, "parameter", outIndex, "h2o_col_mm");
    cpl_table_set_double(*outtable, "value", outIndex, h2ocol);
    if (fith2o == 1) {
        cpl_table_set_double(*outtable, "uncertainty", outIndex, reldelh2o * h2ocol);
    }
    else {
        cpl_table_set_double(*outtable, "uncertainty", outIndex, -1);
    }

    fclose(stream);

    cpl_table_save(*outtable, NULL, NULL, outfitsfile, CPL_IO_CREATE);

    return err;
}

/**@}*/
