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
 * \file mf_modsim.c
 *
 * Routines for preparation of model spectra required for the mf_mpfit
 * fitting procedure
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  06 Jul 2010
 * \date   26 Mar 2015
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <mf_mpfit.h>
#include <mf_modsim.h>


/*****************************************************************************
 *                                GLOBALS                                    *
 ****************************************************************************/

/* Declaration of global variables */

/* Reference fit parameters */
extern cpl_array *reffitpar;
/* Number of MPFIT fitting function calls */
extern int nfev;


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code mf_modsim_lblrtm(const cpl_table *prof,
                                const mfdrv *drvpar, const cpl_array *fitpar,
                                cpl_error_code** arr_code_stat) {

    // First Part of modsim: the call to LBLRTM

    cpl_error_code codestat = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_table *modprof = NULL;
    int lastcall = 0, ismolcol = 1, singlespec = 0, nrange = 0, j = 0;
    int nmod = 0;

    /* Get flag for last call */
    p = cpl_parameterlist_find(drvpar->parlist, "lastcall");
    lastcall = cpl_parameter_get_int(p);

    /* Check whether a relative molecular column density was changed
       (only for MOLECFIT not CALCTRANS) */
    if (nfev != 1 || lastcall != 1) {
        mf_modsim_compfitpar(&ismolcol, drvpar, fitpar);
    }

    /* Run LBLRTM in any case if first or last function call */
    if (nfev == 1 || lastcall == 1) {
        ismolcol = 1;
    }

    /* Modify profiles of atmospheric molecules if necessary */
    if (ismolcol == 1) {
        modprof = cpl_table_duplicate(prof);
        cpl_table* modprof_out = NULL;
        mf_modsim_modatm(modprof, &modprof_out, drvpar, fitpar);
        if (modprof_out) cpl_table_delete(modprof_out);
    }

    /* Calculate one molecular spectrum for full/selected wavelength range? */
    p = cpl_parameterlist_find(drvpar->parlist, "singlespec");
    singlespec = cpl_parameter_get_int(p);

    /* Get number of fit ranges */
    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    nrange = cpl_parameter_get_int(p);

    /* Adapt model spectrum for each part (chip) of the observed spectrum */
    *arr_code_stat = malloc(nrange * sizeof(cpl_error_code));
    for (j = 0; j < nrange; j++) {
        (*arr_code_stat)[j] = CPL_ERROR_NONE;
    }

    for (j = 0; j < nrange; j++) {

        /* Skip empty ranges */
        if (cpl_table_get(drvpar->rangetab, "wn_end", j, NULL) == 0.) {
            continue;
        }

        /* Get number of LBLRTM spectra */
        if (singlespec == 1) {
            nmod = 1;
        } else {
            nmod = j + 1;
        }

        /* Compute LBLRTM spectrum if necessary */
        if (ismolcol == 1 &&
            (singlespec == 0 || (singlespec == 1 && j == 0))) {
            codestat = mf_lblrtm_call(drvpar, modprof, nmod);
            (*arr_code_stat)[j] = codestat;
        }

        /* LBLRTM errors -> output: model spectrum = observed spectrum */
        if (codestat != CPL_ERROR_NONE) {
            cpl_table_delete(modprof);
            return codestat;
        }

    }

    /* Free memory */
    if (ismolcol == 1) {
        cpl_table_delete(modprof);
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}

cpl_error_code mf_modsim_conv(cpl_table *spec,
                              const mfdrv *drvpar, const cpl_array *fitpar,
                              cpl_error_code* arr_code_stat){

    // Second Part: Convolution

    cpl_error_code codestat = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_array *origselrows = NULL, *selrows = NULL;
    cpl_table *modspec = NULL, *rangespec = NULL;
    cpl_table *extmodspec = NULL;
    char kernelname[MF_MAXLEN] = "";
    int lastcall = 0, m = 0, singlespec = 0, nrange = 0, j = 0;
    int nmod = 0, nsel = 0, chip = 0, h = 0, i = 0, range = 0;
    double limlam[2] = {0., 0.}, mflux = 0., dev = 0., chi2 = 0.;

    /* Set model columns to default values */
    m = cpl_table_get_nrow(spec);
    //cpl_table_fill_column_window(spec, "iscont", 0, m, -1);
    cpl_table_fill_column_window(spec, "mlambda", 0, m, 0.);
    cpl_table_fill_column_window(spec, "mscal", 0, m, 1.);
    cpl_table_fill_column_window(spec, "mflux", 0, m, 0.);
    cpl_table_fill_column_window(spec, "mweight", 0, m, 0.);
    cpl_table_fill_column_window(spec, "dev", 0, m, 0.);

    /* Calculate one molecular spectrum for full/selected wavelength range? */
    p = cpl_parameterlist_find(drvpar->parlist, "singlespec");
    singlespec = cpl_parameter_get_int(p);

    /* Fixed kernel or synthetic kernel? */
    p = cpl_parameterlist_find(drvpar->parlist, "kernel_file");
    strncpy(kernelname, cpl_parameter_get_string(p), MF_MAXLEN);

    /* Get number of fit ranges */
    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    nrange = cpl_parameter_get_int(p);

    /* Adapt model spectrum for each part (chip) of the observed spectrum */

    for (j = 0; j < nrange; j++) {

        /* Skip empty ranges */
        if (cpl_table_get(drvpar->rangetab, "wn_end", j, NULL) == 0.) {
            continue;
        }

        /* Get number of LBLRTM spectra */
        if (singlespec == 1) {
            nmod = 1;
        } else {
            nmod = j + 1;
        }

        /* Read LBLRTM spectrum and convert wavenumbers to wavelengths */
        if (singlespec == 0 || (singlespec == 1 && j == 0)) {
            modspec = cpl_table_new(0);
            codestat = arr_code_stat[j];
            if (codestat == CPL_ERROR_NONE) {
                codestat = mf_modsim_readlblspec(modspec, drvpar, nmod);
            }
        }

        /* LBLRTM errors -> output: model spectrum = observed spectrum */
        if (codestat != CPL_ERROR_NONE) {
            cpl_table_delete(modspec);
            return codestat;
        }

        /* Extract wavelength grid of observed spectrum for each range */
        cpl_table_unselect_all(spec);
        nsel = cpl_table_or_selected_int(spec, "mrange", CPL_EQUAL_TO, j+1);
        rangespec = cpl_table_extract_selected(spec);
        origselrows = cpl_table_where_selected(spec);
        selrows = cpl_array_cast(origselrows, CPL_TYPE_INT);
        cpl_array_delete(origselrows);
        cpl_table_select_all(spec);
        cpl_table_erase_column(rangespec, "chip");
        cpl_table_name_column(rangespec, "flux", "oflux");
        cpl_table_erase_column(rangespec, "weight");
        cpl_table_name_column(rangespec, "mscal", "scal");
        cpl_table_erase_column(rangespec, "mflux");
        cpl_table_erase_column(rangespec, "mweight");
        cpl_table_erase_column(rangespec, "dev");

        /* Extract chip-related wavelength range from model spectrum */
        limlam[0] = cpl_table_get(rangespec, "lambda", 0, NULL);
        limlam[1] = cpl_table_get(rangespec, "lambda", nsel-1, NULL);
        limlam[0] = 1e4 / (1e4 / limlam[0] + MF_EXTRACOVER);
        limlam[1] = 1e4 / (1e4 / limlam[1] - MF_EXTRACOVER);
        cpl_table_or_selected_double(modspec, "lambda", CPL_GREATER_THAN,
                                     limlam[0]);
        cpl_table_and_selected_double(modspec, "lambda", CPL_NOT_GREATER_THAN,
                                     limlam[1]);
        extmodspec = cpl_table_extract_selected(modspec);

        /* Delete model spectrum if no more required */
        if (singlespec == 0 || (singlespec == 1 && j == nrange-1)) {
            cpl_table_delete(modspec);
        }

        /* Get chip number for selected fit range */
        chip = cpl_table_get(drvpar->rangetab, "chip", j, NULL);

        /* Modify wavelength grid of model spectrum by means of Chebyshev
           polynomials */
        mf_modsim_modwavegrid(extmodspec, drvpar, fitpar, chip);

        /* Rebin model spectrum to wavelength grid of observed spectrum */
        mf_basic_rebin(rangespec, "lambda", "flux", extmodspec, "lambda",
                       "flux");
        mf_basic_rebin(rangespec, "lambda", "mlambda", extmodspec, "lambda",
                       "lambda0");
        cpl_table_delete(extmodspec);

        /* Convolve model spectrum with kernel(s) */
        if (strcmp(kernelname, "none") == 0) {
            /* Wavelength-dependent boxcar, Gaussian, and/or Lorentzian
               kernels */
            printf("[ DEBUG ] Convolving with a synthetic kernel.\n");
            mf_modsim_convolvesynthkernel(rangespec, drvpar, fitpar);
        } else {
            /* Fixed kernel from file */
            printf("[ DEBUG ] Convolving with a user-defined kernel.\n");
            mf_modsim_convolvereadkernel(rangespec, selrows, drvpar);
        }

        /* Add grey body in the case of a radiance spectrum */
        mf_modsim_telback(rangespec, drvpar, fitpar);

        /* Adapt flux units in the case of a radiance spectrum */
        mf_modsim_fluxunits(rangespec, drvpar);

        /* Modify continuum of model spectrum by means of polynomials */
        mf_modsim_modcont(rangespec, drvpar, fitpar, j+1);

        /* Correct model by interpolated continuum of observed spectrum
           (skipped in the case of a radiance spectrum) */
        //mf_modsim_corrobsspec(rangespec, drvpar);
        //cpl_table_multiply_columns(rangespec, "flux", "scal");

        /* Write resulting spectrum in output CPL table */
        for (h = 0; h < nsel; h++) {
            i = cpl_array_get(selrows, h, NULL);
            //cpl_table_set(spec, "iscont", i,
            //              cpl_table_get(rangespec, "iscont", h, NULL));
            cpl_table_set(spec, "mlambda", i,
                          cpl_table_get(rangespec, "mlambda", h, NULL));
            cpl_table_set(spec, "mscal", i,
                          cpl_table_get(rangespec, "scal", h, NULL));
            cpl_table_set(spec, "mflux", i,
                          cpl_table_get(rangespec, "flux", h, NULL));
        }

        /* Free memory */
        cpl_array_delete(selrows);
        cpl_table_delete(rangespec);

    }

    for (i = 0; i < m; i++) {
        range = cpl_table_get(spec, "mrange", i, NULL);
        mflux = cpl_table_get(spec, "mflux", i, NULL);
        if (range == 0 || mflux < 0 || isnan(mflux) != 0) {
            cpl_table_set(spec, "mflux", i, 0.);
            cpl_table_set(spec, "mweight", i, 0.);
        } else {
            cpl_table_set(spec, "mweight", i, 1.);
        }
    }

    /* Calculate weighted deviations between modelled and observed spectrum */
    cpl_table_add_columns(spec, "dev", "mflux");
    cpl_table_subtract_columns(spec, "dev", "flux");
    cpl_table_multiply_columns(spec, "dev", "weight");
    cpl_table_multiply_columns(spec, "dev", "mweight");

    /* Print chi^2 (only for MOLECFIT not CALCTRANS) */
    if (nfev != 1 || lastcall != 1) {
        for (chi2 = 0, i = 0; i < m; i++) {
            dev = cpl_table_get(spec, "dev", i, NULL);
            chi2 += dev * dev;
        }
        cpl_msg_info(cpl_func, "chi2 = %g", chi2);
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim(cpl_table *spec,
    cpl_table **prof_out,
    const cpl_table *prof,
    const mfdrv *drvpar,
    const cpl_array *fitpar)
{
    /*!
     * \callgraph
     *
     * Calculates sky model spectrum for given atmospheric profiles and adapts
     * it to wavelength grid, resolution, and continuum of the observed
     * spectrum. Finally, weighted deviations between the model and the
     * observed spectrum are computed.
     *
     * \b INPUT:
     * \param spec    CPL table with observed spectrum
     * \param prof    CPL table with atmospheric profiles
     * \param drvpar  ::mfdrv parameter structure
     * \param fitpar  CPL array with fit parameters
     *
     * \b OUTPUT:
     * \param spec    CPL table with observed and modelled spectrum
     *
     * \b ERRORS:
     * - see ::mf_lblrtm_call
     */

    cpl_error_code codestat = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_array *origselrows = NULL, *selrows = NULL;
    cpl_table *modprof = NULL, *modspec = NULL, *rangespec = NULL;
    cpl_table *extmodspec = NULL;
    char kernelname[MF_MAXLEN] = "";
    int lastcall = 0, ismolcol = 1, m = 0, singlespec = 0, nrange = 0, j = 0;
    int nmod = 0, nsel = 0, chip = 0, h = 0, i = 0, range = 0;
    double limlam[2] = {0., 0.}, mflux = 0., dev = 0., chi2 = 0.;

    /* Get flag for last call */
    p = cpl_parameterlist_find(drvpar->parlist, "lastcall");
    lastcall = cpl_parameter_get_int(p);

    /* Check whether a relative molecular column density was changed
       (only for MOLECFIT not CALCTRANS) */
    if (nfev != 1 || lastcall != 1) {
        mf_modsim_compfitpar(&ismolcol, drvpar, fitpar);
    }

    /* Run LBLRTM in any case if first or last function call */
    if (nfev == 1 || lastcall == 1) {
        ismolcol = 1;
    }

    /* Modify profiles of atmospheric molecules if necessary */
    if (ismolcol == 1) {
        modprof = cpl_table_duplicate(prof);
        mf_modsim_modatm(modprof, prof_out, drvpar, fitpar);
    }

    /* Set model columns to default values */
    m = cpl_table_get_nrow(spec);
    //cpl_table_fill_column_window(spec, "iscont", 0, m, -1);
    cpl_table_fill_column_window(spec, "mlambda", 0, m, 0.);
    cpl_table_fill_column_window(spec, "mscal", 0, m, 1.);
    cpl_table_fill_column_window(spec, "mflux", 0, m, 0.);
    cpl_table_fill_column_window(spec, "mweight", 0, m, 0.);
    cpl_table_fill_column_window(spec, "dev", 0, m, 0.);

    /* Calculate one molecular spectrum for full/selected wavelength range? */
    p = cpl_parameterlist_find(drvpar->parlist, "singlespec");
    singlespec = cpl_parameter_get_int(p);

    /* Fixed kernel or synthetic kernel? */
    p = cpl_parameterlist_find(drvpar->parlist, "kernel_file");
    strncpy(kernelname, cpl_parameter_get_string(p), MF_MAXLEN);

    /* Get number of fit ranges */
    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    nrange = cpl_parameter_get_int(p);

    /* Adapt model spectrum for each part (chip) of the observed spectrum */

    for (j = 0; j < nrange; j++) {

        /* Skip empty ranges */
        if (cpl_table_get(drvpar->rangetab, "wn_end", j, NULL) == 0.) {
            continue;
        }

        /* Get number of LBLRTM spectra */
        if (singlespec == 1) {
            nmod = 1;
        } else {
            nmod = j + 1;
        }

        /* Compute LBLRTM spectrum if necessary */
        if (ismolcol == 1 &&
            (singlespec == 0 || (singlespec == 1 && j == 0))) {
            codestat = mf_lblrtm_call(drvpar, modprof, nmod);
        }

        /* Read LBLRTM spectrum and convert wavenumbers to wavelengths */
        if (singlespec == 0 || (singlespec == 1 && j == 0)) {
            modspec = cpl_table_new(0);
            if (codestat == CPL_ERROR_NONE) {
                codestat = mf_modsim_readlblspec(modspec, drvpar, nmod);
            }
        }

        /* LBLRTM errors -> output: model spectrum = observed spectrum */
        if (codestat != CPL_ERROR_NONE) {
            cpl_table_delete(modprof);
            cpl_table_delete(modspec);
            return codestat;
        }

        /* Extract wavelength grid of observed spectrum for each range */
        cpl_table_unselect_all(spec);
        nsel = cpl_table_or_selected_int(spec, "mrange", CPL_EQUAL_TO, j+1);
        rangespec = cpl_table_extract_selected(spec);
        origselrows = cpl_table_where_selected(spec);
        selrows = cpl_array_cast(origselrows, CPL_TYPE_INT);
        cpl_array_delete(origselrows);
        cpl_table_select_all(spec);
        cpl_table_erase_column(rangespec, "chip");
        cpl_table_name_column(rangespec, "flux", "oflux");
        cpl_table_erase_column(rangespec, "weight");
        cpl_table_name_column(rangespec, "mscal", "scal");
        cpl_table_erase_column(rangespec, "mflux");
        cpl_table_erase_column(rangespec, "mweight");
        cpl_table_erase_column(rangespec, "dev");

        /* Extract chip-related wavelength range from model spectrum */
        limlam[0] = cpl_table_get(rangespec, "lambda", 0, NULL);
        limlam[1] = cpl_table_get(rangespec, "lambda", nsel-1, NULL);
        limlam[0] = 1e4 / (1e4 / limlam[0] + MF_EXTRACOVER);
        limlam[1] = 1e4 / (1e4 / limlam[1] - MF_EXTRACOVER);
        cpl_table_or_selected_double(modspec, "lambda", CPL_GREATER_THAN,
                                     limlam[0]);
        cpl_table_and_selected_double(modspec, "lambda", CPL_NOT_GREATER_THAN,
                                     limlam[1]);
        extmodspec = cpl_table_extract_selected(modspec);

        /* Delete model spectrum if no more required */
        if (singlespec == 0 || (singlespec == 1 && j == nrange-1)) {
            cpl_table_delete(modspec);
        }

        /* Get chip number for selected fit range */
        chip = cpl_table_get(drvpar->rangetab, "chip", j, NULL);

        /* Modify wavelength grid of model spectrum by means of Chebyshev
           polynomials */
        mf_modsim_modwavegrid(extmodspec, drvpar, fitpar, chip);

        /* Rebin model spectrum to wavelength grid of observed spectrum */
        mf_basic_rebin(rangespec, "lambda", "flux", extmodspec, "lambda",
                       "flux");
        mf_basic_rebin(rangespec, "lambda", "mlambda", extmodspec, "lambda",
                       "lambda0");
        cpl_table_delete(extmodspec);

        /* Convolve model spectrum with kernel(s) */
        if (strcmp(kernelname, "none") == 0) {
            /* Wavelength-dependent boxcar, Gaussian, and/or Lorentzian
               kernels */
            printf("[ DEBUG ] Convolving with a synthetic kernel.\n");
            mf_modsim_convolvesynthkernel(rangespec, drvpar, fitpar);
        } else {
            /* Fixed kernel from file */
            printf("[ DEBUG ] Convolving with a user-defined kernel.\n");
            mf_modsim_convolvereadkernel(rangespec, selrows, drvpar);
        }

        /* Add grey body in the case of a radiance spectrum */
        mf_modsim_telback(rangespec, drvpar, fitpar);

        /* Adapt flux units in the case of a radiance spectrum */
        mf_modsim_fluxunits(rangespec, drvpar);

        /* Modify continuum of model spectrum by means of polynomials */
        mf_modsim_modcont(rangespec, drvpar, fitpar, j+1);

        /* Correct model by interpolated continuum of observed spectrum
           (skipped in the case of a radiance spectrum) */
        //mf_modsim_corrobsspec(rangespec, drvpar);
        //cpl_table_multiply_columns(rangespec, "flux", "scal");

        /* Write resulting spectrum in output CPL table */
        for (h = 0; h < nsel; h++) {
            i = cpl_array_get(selrows, h, NULL);
            //cpl_table_set(spec, "iscont", i,
            //              cpl_table_get(rangespec, "iscont", h, NULL));
            cpl_table_set(spec, "mlambda", i,
                          cpl_table_get(rangespec, "mlambda", h, NULL));
            cpl_table_set(spec, "mscal", i,
                          cpl_table_get(rangespec, "scal", h, NULL));
            cpl_table_set(spec, "mflux", i,
                          cpl_table_get(rangespec, "flux", h, NULL));
        }

        /* Free memory */
        cpl_array_delete(selrows);
        cpl_table_delete(rangespec);

    }

    /* Free memory */
    if (ismolcol == 1) {
        cpl_table_delete(modprof);
    }

    /* Set mflux = 0 and mweight = 0 for pixels outside ranges and with
       negative or nan (not a number) model flux */
    for (i = 0; i < m; i++) {
        range = cpl_table_get(spec, "mrange", i, NULL);
        mflux = cpl_table_get(spec, "mflux", i, NULL);
        if (range == 0 || mflux < 0 || isnan(mflux) != 0) {
            cpl_table_set(spec, "mflux", i, 0.);
            cpl_table_set(spec, "mweight", i, 0.);
        } else {
            cpl_table_set(spec, "mweight", i, 1.);
        }
    }

    /* Calculate weighted deviations between modelled and observed spectrum */
    cpl_table_add_columns(spec, "dev", "mflux");
    cpl_table_subtract_columns(spec, "dev", "flux");
    cpl_table_multiply_columns(spec, "dev", "weight");
    cpl_table_multiply_columns(spec, "dev", "mweight");

    /* Print chi^2 (only for MOLECFIT not CALCTRANS) */
    if (nfev != 1 || lastcall != 1) {
        for (chi2 = 0, i = 0; i < m; i++) {
            dev = cpl_table_get(spec, "dev", i, NULL);
            chi2 += dev * dev;
        }
        cpl_msg_info(cpl_func, "chi2 = %g", chi2);
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_compfitpar(int *ismolcol, const mfdrv *drvpar,
                                    const cpl_array *fitpar)
{
    /*!
     * Checks the modification of a molecular column density fit parameter
     * and returns a flag (yes: 1, no: 0, error: -1). Moreover, the current
     * fit parameters are saved in the reference fit parameter array
     * ::reffitpar.
     *
     * \b INPUT:
     * \param drvpar    ::mfdrv parameter structure
     * \param fitpar    CPL array with fit parameters
     *
     * \b OUTPUT:
     * \param ismolcol  flag indicating the modification of a molecular
     *                  parameter (yes: 1, no: 0, error: -1)
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    cpl_parameter *p;
    char errtxt[MF_MAXLEN];
    int n = 0, nref = 0, nmolec = 0, i = 0;
    double par = 0., refpar = 0., reldev = 0.;

    /* Check number of array values */
    n = cpl_array_get_size(fitpar);
    nref = cpl_array_get_size(reffitpar);
    if (n == 0 && n == nref) {
        *ismolcol = 1;
        return CPL_ERROR_NONE;
    } else if (n != nref) {
        *ismolcol = -1;
        sprintf(errtxt, "%s: cpl_array *fitpar != cpl_array *reffitpar "
                "(size)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get number of molecules in driver parameter structure */
    p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(p);

    /* Check modification of molecular column density parameter and copy
       fitpar to reffitpar */
    *ismolcol = 0;
    for (i = 0; i < n; i++) {
        par = cpl_array_get(fitpar, i, NULL);
        refpar = cpl_array_get(reffitpar, i, NULL);
        reldev = fabs(par - refpar);
        if (refpar == 0.) {
            if (par == 0.) {
                reldev = 0.;
            } else {
                reldev /= par;
            }
        } else {
            reldev /= refpar;
        }
        if (reldev > MF_TOL || nfev == 1) {
            cpl_msg_info("mf_modsim", "%4d %3d %9.5f %8.5f",
                         nfev, i, par, reldev);
        }
        if (i < nmolec && reldev > MF_TOL) {
            *ismolcol = 1;
        }
        cpl_array_set(reffitpar, i, par);
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_modatm(cpl_table *modprof, cpl_table** modprof_out,
    const mfdrv *drvpar,
    const cpl_array *fitpar)
{
    /*!
     * \callgraph
     *
     * Multiplies profiles of molecules listed in the ::mfdrv structure by
     * correction factors given by the CPL array of fit parameters.
     * The resulting CPL table of atmospheric profiles is written into an
     * ASCII output file (parameter "output_name" in ::mfdrv structure +
     * "_fit.atm") which is required as input file for the LBLRTM code.
     *
     * \b INPUT:
     * \param modprof  CPL table with atmospheric profiles
     * \param drvpar   ::mfdrv parameter structure
     * \param fitpar   CPL array with fit parameters
     *
     * \b OUTPUT:
     * \param modprof  atmospheric profiles with corrected molecular
     *                 frequencies
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    cpl_error_code err = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_array *colarray;
    cpl_boolean ext = CPL_FALSE;
    char **colname, **mol, errtxt[MF_MAXLEN], basedir[MF_MAXLEN];
    char output_dir[MF_MAXLEN], output_name[MF_MAXLEN];
    char atmfile[MF_MAXLEN];
    int nmolcol = 0, nmolec = 0, i = 0, j = 0;
    const double *par;

    /* Get number of columns with molecular profiles */
    nmolcol = cpl_table_get_ncol(modprof) - 3;

    /* Get pointer to column names */
    colarray = cpl_table_get_column_names(modprof);
    colname = cpl_array_get_data_string(colarray);

    /* Get number of molecules in driver parameter structure */
    p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(p);

    /* Get pointer to names of molecules in driver parameter structure */
    mol = cpl_table_get_data_string(drvpar->molectab, "list_molec");

    /* Get pointer to CPL array with fit parameters */
    par = cpl_array_get_data_double_const(fitpar);

    /* Modify column densities of molecules */

    for (i = 0; i < nmolec; i++) {
        ext = CPL_FALSE;
        for (j = 0; j < nmolcol; j++) {
           if (strcmp(colname[j+3], mol[i]) == 0) {
                ext = CPL_TRUE;
                cpl_table_multiply_scalar(modprof, colname[j+3], par[i]);
            }
        }
        if (ext == CPL_FALSE) {
            /* Molecule of the driver file is not present in list of CPL
               table column names */
            sprintf(errtxt, "%s: cpl_table *modprof (column %s not found)",
                    MF_ERROR_IOS_TXT, mol[i]);
            err = cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
        }
    }

    /* Write modified CPL table into file -> input for LBLRTM */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(output_dir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strcpy(output_name, cpl_parameter_get_string(p));
    sprintf(atmfile, "%s%s_fit.atm", output_dir, output_name);
    mf_atm_writeatm(modprof, modprof_out, atmfile);

    /* Free memory */
    cpl_array_delete(colarray);

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return err;
}


cpl_error_code mf_modsim_readlblspec(cpl_table *modspec, const mfdrv *drvpar,
                                     const int range)
{
    /*!
     * \callgraph
     *
     * Reads LBLRTM-type spectrum for given fit range from FITS table into CPL
     * table. The fluxes are converted from
     * \f${\rm W\,cm}^{-2}\,{\rm cm\,sr}^{-1}\f$ to
     * \f${\rm phot\,s}^{-1}{\rm m}^{-2}\mu{\rm m}^{-1}{\rm arcsec}^{-1}\f$.
     *
     * \b INPUT:
     * \param drvpar   ::mfdrv parameter structure
     * \param range    range number
     *
     * \b OUTPUT:
     * \param modspec  CPL table with read LBLRTM spectrum in wavelength units
     *
     * \b ERRORS:
     * - Error in subroutine
     * - see ::mf_basic_access
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_table *spec;
    cpl_boolean isnanum = CPL_FALSE;
    char basedir[MF_MAXLEN], outdir[MF_MAXLEN], filename[MF_MAXLEN];
    char errtxt[MF_MAXLEN];
    const char *outfile;
    int trans = 0, i = 0;
    const double conv = 1e-4 /
                 (MF_LAM_UNIT * CPL_PHYS_C * CPL_PHYS_H * MF_SR_IN_ARCSEC2);

    /* Get base directory */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));

    /* Get output directory */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);

    /* Get output filename */
    // p = cpl_parameterlist_find(drvpar->parlist, "spec_out");
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    outfile = cpl_parameter_get_string(p);

    /* Get requested spectrum type */
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);

    /* Compose filename */
    if (trans == 0) {
      sprintf(filename, "%s%s_%d_R.fits", outdir, outfile, range);
    } else {
      sprintf(filename, "%s%s_%d_T.fits", outdir, outfile, range);
    }

    /* Check file existence */
    if ((status = mf_basic_access(filename, F_OK)) != CPL_ERROR_NONE) {
        return status;
    }

    /* Read spectrum from FITS file into CPL table */
    if ((spec = cpl_table_load(filename, 1, 0)) == NULL) {
        cpl_table_delete(spec);
        sprintf(errtxt, "%s: cpl_table_load(%s, 1, 0)", MF_ERROR_EIS_TXT,
                filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_EIS, "%s", errtxt);
    }

    cpl_table_set_size(modspec, cpl_table_get_nrow(spec));
    cpl_table_new_column(modspec, "lambda", CPL_TYPE_DOUBLE);
    cpl_table_new_column(modspec, "flux", CPL_TYPE_DOUBLE);
    cpl_table_copy_data_double(modspec, "lambda",
                               cpl_table_get_data_double(spec, "lambda"));
    cpl_table_copy_data_double(modspec, "flux",
                               cpl_table_get_data_double(spec, "flux"));
    cpl_table_delete(spec);

    /* Convert fluxes if radiance spectrum */
    if (trans == 0) {
        cpl_table_divide_columns(modspec, "flux", "lambda");
        cpl_table_multiply_scalar(modspec, "flux", conv);
    }

    /* Avoid "nan" (not a number) in input table */
    for (i = 0; i < cpl_table_get_nrow(modspec); i++) {
        if (isnan(cpl_table_get(modspec, "flux", i, NULL)) != 0) {
            cpl_table_set(modspec, "flux", i, -9e99);
            isnanum = CPL_TRUE;
        }
    }

    /* Print warning message in the case of "nan" values */
    if (isnanum == CPL_TRUE) {
        cpl_msg_warning(cpl_func, "LBLRTM produces nan values -> weight = 0");
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_modwavegrid(cpl_table *spec, const mfdrv *drvpar,
                                     const cpl_array *fitpar, const int chip)
{
     /*!
     * Modifies the wavelength scale of a model spectrum by means of
     * Chebyshev polynomials. The coefficients are provided by the fit
     * parameter vector (CPL array). The adaption of the wavelength grid is
     * carried out for the given chip only. The input wavelength grid is
     * saved and written into the output column "lambda0".
     *
     * \b INPUT:
     * \param spec    CPL table with model spectrum
     * \param drvpar  ::mfdrv parameter structure
     * \param fitpar  CPL array with fit parameters
     * \param chip    chip number
     *
     * \b OUTPUT:
     * \param spec    model spectrum with modified wavelength grid
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    cpl_parameter *p;
    cpl_array *cheby;
    char errtxt[MF_MAXLEN];
    int nmolec = 0, nwlc = 0, n0 = 0, nlam = 0, i = 0, j = 0;
    double limlam[2] = {0., 0.}, dlam = 0., wmean = 0.;
    const double *par;
    double *lam, *slam, *t;

    /* Find position of wavelength calibration parameters in fit parameter
       CPL array */
    p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "wlc_n");
    nwlc = cpl_parameter_get_int(p) + 1;
    n0 = nmolec + nwlc * (chip - 1);
    if (nwlc-1 < 1) {
        sprintf(errtxt, "%s: wlc_n of mfdrv *drvpar < 1", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Get pointer to CPL array with fit parameters */
    par = cpl_array_get_data_double_const(fitpar);

    /* Save input wavelength grid */
    cpl_table_duplicate_column(spec, "lambda0", spec, "lambda");

    /* Scale wavelengths in the way that the interval [-1,1] is covered */
    limlam[0] = cpl_table_get(drvpar->chiptab, "wl_min", chip-1, NULL);
    limlam[1] = cpl_table_get(drvpar->chiptab, "wl_max", chip-1, NULL);
    dlam = limlam[1] - limlam[0];
    wmean = (limlam[0] + limlam[1]) / 2;
    cpl_table_duplicate_column(spec, "slambda", spec, "lambda");
    cpl_table_subtract_scalar(spec, "slambda", wmean);
    cpl_table_divide_scalar(spec, "slambda", dlam / 2.);

    /* Initialise wavelength column (all values = 0) */
    nlam = cpl_table_get_nrow(spec);
    cpl_table_fill_column_window(spec, "lambda", 0, nlam, 0.);

    /* Get pointers to CPL table columns */
    lam = cpl_table_get_data_double(spec, "lambda");
    slam = cpl_table_get_data_double(spec, "slambda");

    /* Create array for Chebyshev polynomials */
    cheby = cpl_array_new(nwlc, CPL_TYPE_DOUBLE);
    t = cpl_array_get_data_double(cheby);
    t[0] = 1;

    /* Compute Chebyshev polynomials */
    for (i = 0; i < nlam; i++) {
        for (j = 0; j < nwlc; j++) {
            if (j == 1) {
                t[j] = slam[i];
            } else if (j > 1) {
                t[j] = 2 * slam[i] * t[j-1] - t[j-2];
            }
            lam[i] += par[n0 + j] * t[j];
        }
    }

    /* Rescale wavelengths */
    cpl_table_multiply_scalar(spec, "lambda", dlam / 2.);
    cpl_table_add_scalar(spec, "lambda", wmean);

    /* Free memory */
    cpl_array_delete(cheby);
    cpl_table_erase_column(spec, "slambda");

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_convolvesynthkernel(cpl_table *spec,
                                             const mfdrv *drvpar,
                                             const cpl_array *fitpar)
{
    /*!
     * \callgraph
     *
     * Applies synthetic kernels to a model spectrum as provided by a CPL
     * table with columns "lambda" and "flux". Depending on the parameter
     * \e varkern, either a constant or a variable kernel is applied in the
     * case of a synthesis based on a boxcar, Gaussian, and Lorentzian
     * component (see ::mf_modsim_calckerneltable for more details). The
     * variable kernel increases linearly with wavelength, as it is expected
     * for constant resolution. For reasonable results, the input wavelength
     * grid should have a constant step size.
     *
     * \b INPUT:
     * \param spec    input spectrum (CPL table)
     * \param drvpar  ::mfdrv parameter structure
     * \param fitpar  CPL array with fit parameters
     *
     * \b OUTPUT:
     * \param spec    convolved spectrum
     *
     * \b ERRORS:
     * - none
     */

    cpl_table *kerntab = NULL;

    /* Create kernel table */
    kerntab = cpl_table_new(0);

    /* Calculate kernel table */
    mf_modsim_calckerneltable(kerntab, spec, drvpar, fitpar);

    /* Invert kernel table */
    mf_modsim_invertkerneltable(kerntab);

    /* Convolve flux column of spectrum with inverted kernel(s) */
    mf_modsim_convolvekerneltable(spec, kerntab);

    /* Free memory */
    cpl_table_delete(kerntab);

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_calckerneltable(cpl_table *kerntab,
                                         const cpl_table *spec,
                                         const mfdrv *drvpar,
                                         const cpl_array *fitpar)
{
    /*!
     * \callgraph
     *
     * Calculates convolution kernels for different pixel ranges of an input
     * spectrum and stores them in an output table. Each row of this table
     * gives a pixel range and the correponding kernel as array.
     *
     * The kernels are calculated as a combination of a boxcar, Gaussian, and
     * Lorentzian kernel. Depending on the parameter \e kernmode, the latter
     * two components can be computed either independently or via a Voigt
     * profile approximation. The width of the boxcar is computed from the
     * width of the spectrograph slit in pixels times a correction factor
     * provided by the fitpar CPL array. The width of the Gaussian and the
     * Lorentzian are directly taken from the fitpar CPL array. The number of
     * pixels of the latter and the Voigt profile depends on the parameter
     * \e kernfac times FWHM.
     *
     * If the \e varkern parameter is set to 1, the FWHM of all kernel
     * components linearly increase with wavelength, which results in a
     * constant resolution. In this case, the FWHM from the fitpar array are
     * related to the centre of the full wavelength range (considering the
     * data of all chips). The number of kernels provided by the output table
     * depends on ::MF_LIMRELLAMVAR, which provides the relative wavelength
     * change that causes a recalculation. If \e varkern is set to 0, only a
     * single kernel is calculated, i.e. the output table has only one row.
     *
     * For a correct application of the routine, it is required that the
     * wavelength difference of adjacent pixels is a constant for the entire
     * spectrum.
     *
     * \b INPUT:
     * \param kerntab   empty CPL table
     * \param spec      input spectrum (CPL table)
     * \param drvpar    ::mfdrv parameter structure
     * \param fitpar    CPL array with fit parameters
     *
     * \b OUTPUT:
     * \param kerntab   table with kernels for different pixel ranges
     *
     * \b ERRORS:
     * - No data
     * - Invalid object value(s)
     * - Invalid object structure
     */

    cpl_parameter *p;
    cpl_array *kernel0 = NULL, *kernel = NULL, **kernarr = NULL;
    cpl_array *tabkernel = NULL;
    char errtxt[MF_MAXLEN];
    int m = 0, varkern = 0, nchip = 0, nres = 0, trans = 0, kernmode = 0;
    int np0 = 0, p0c = 0, n1r = 0, n2r = 0, npmax = 0, nkern = 0, j = 0;
    int pminj = 0, pmaxj = 0, jo = 0, npj = 0, k = 0;
    int *pmin, *pmax, *p0, *np;
    double wspec[3] = {0., 0., 0.}, dlam = 0., limlam[2] = {0., 0.};
    double reflam = 0., slitw = 0., pixsc = 0., wbox = 0., kernfac = 0.;
    double wgauss = 0., wlorentz = 0., x0 = 0., dx0 = 0., r0 = 0., c = 0.;
    double xmin = 0., xmax = 0., dxmaxr = 0., clam = 0., llam = 0., ulam = 0.;
    double cwbox = 0., cwgauss = 0., cwlorentz = 0.;
    double *kern0, *kern, *tabkern;

    /* Check number of data points in spectrum */
    m = cpl_table_get_nrow(spec);
    if (m <= 0) {
        sprintf(errtxt, "%s: cpl_table *spec", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Get wavelength range, central wavelength, and wavelength step of input
       spectrum */
    wspec[0] = cpl_table_get(spec, "lambda", 0, NULL);
    wspec[1] = cpl_table_get(spec, "lambda", m-1, NULL);
    wspec[2] = 0.5 * (wspec[0] + wspec[1]);
    dlam = (wspec[1] - wspec[0]) / (m - 1);

    /* Kernel: linearly increasing with wavelength or constant? */
    p = cpl_parameterlist_find(drvpar->parlist, "varkern");
    varkern = cpl_parameter_get_int(p);

    /* Get central wavelength */
    nchip = cpl_table_get_nrow(drvpar->chiptab);
    limlam[0] = cpl_table_get(drvpar->chiptab, "wl_min", 0, NULL);
    limlam[1] = cpl_table_get(drvpar->chiptab, "wl_max", nchip-1, NULL);
    reflam = 0.5 * (limlam[0] + limlam[1]);

    /* Boxcar kernel */

    /* Get slit width in pixels (width of boxcar) */
    p = cpl_parameterlist_find(drvpar->parlist, "slitw");
    slitw = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find(drvpar->parlist, "pixsc");
    pixsc = cpl_parameter_get_double(p);
    if (pixsc == 0.) {
        sprintf(errtxt, "%s: pixsc of mfdrv *drvpar = 0", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }
    wbox = slitw / pixsc;

    /* Get relative width of boxcar from CPL array of fit parameters
       and correct slit width */
    nres = cpl_array_get_size(fitpar) - 3;
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);
    if (trans == 0) {
        nres--;
    }
    wbox *= cpl_array_get(fitpar, nres, NULL);

    /* Gaussian, Lorentzian, and Voigt profile kernels */

    /* Get kernel size in FWHM */
    p = cpl_parameterlist_find(drvpar->parlist, "kernfac");
    kernfac = cpl_parameter_get_double(p);

    /* Get FWHM of Gaussian in pixels from CPL array of fit parameters */
    nres++;
    wgauss = cpl_array_get(fitpar, nres, NULL);

    /* Get FWHM of Lorentzian in pixels from CPL array of fit parameters */
    nres++;
    wlorentz = cpl_array_get(fitpar, nres, NULL);

    /* Independent Gaussian and Lorentzian kernels or Voigt profile
       approximation? */
    p = cpl_parameterlist_find(drvpar->parlist, "kernmode");
    kernmode = cpl_parameter_get_int(p);

    /* Create CPL array for central kernel */
    kernel0 = cpl_array_new(0, CPL_TYPE_DOUBLE);

    /* Calculate synthetic kernel for central wavelength */
    mf_modsim_calcsynthkernel(kernel0, wbox, wgauss, wlorentz, kernfac,
                              kernmode);

    /* Get number of kernel pixels for central wavelength */
    np0 = cpl_array_get_size(kernel0);

    /* Get central pixel of input spectrum */
    p0c = floor(0.5 * m);

    /* Derive maximum kernel size and number of kernels */
    if (varkern == 1) {
        /* Analytical solution for a kernel that linearly increases with
           wavelength and that is only calculated for positions differing by
           a constant factor to an integer power; a constant step size of the
           wavelength grid is required */
        x0 = wspec[2];
        dx0 = 0.5 * (wspec[2] / reflam) * np0 * dlam;
        r0 = wspec[2] - wspec[0];
        c = 1 + MF_LIMRELLAMVAR;
        xmin = (x0 - r0) * x0 / (x0 + dx0);
        n1r = ceil(log(x0 / xmin) / log(c) - 0.5);
        xmax = (x0 + r0) * x0 / (x0 - dx0);
        n2r = floor(log(xmax / x0) / log(c) + 0.5);
        dxmaxr = dx0 * pow(c, n2r);
        npmax = ceil(2 * dxmaxr / dlam) + 4;
        nkern = 1 + n1r + n2r;
    } else {
        npmax = np0;
        nkern = 1;
    }

    /* Create CPL array for kernel of maximum size */
    kernel = cpl_array_new(npmax, CPL_TYPE_DOUBLE);
    cpl_array_fill_window_double(kernel, 0, npmax, 0.);
    tabkernel = cpl_array_new(npmax, CPL_TYPE_DOUBLE);
    cpl_array_fill_window_double(tabkernel, 0, npmax, 0.);

    /* Get pointers to kernel arrays */
    kern0 = cpl_array_get_data_double(kernel0);
    tabkern = cpl_array_get_data_double(tabkernel);

    /* Create output table columns if required */
    if (cpl_table_has_column(kerntab, "pixmin") != 1) {
        cpl_table_new_column(kerntab, "pixmin", CPL_TYPE_INT);
    }
    if (cpl_table_has_column(kerntab, "pixmax") != 1) {
        cpl_table_new_column(kerntab, "pixmax", CPL_TYPE_INT);
    }
    if (cpl_table_has_column(kerntab, "pix0") != 1) {
        cpl_table_new_column(kerntab, "pix0", CPL_TYPE_INT);
    }
    if (cpl_table_has_column(kerntab, "npix") != 1) {
        cpl_table_new_column(kerntab, "npix", CPL_TYPE_INT);
    }
    if (cpl_table_has_column(kerntab, "kernel") != 1) {
        cpl_table_new_column_array(kerntab, "kernel", CPL_TYPE_DOUBLE, npmax);
    }

    /* Initialise output table */
    cpl_table_set_size(kerntab, nkern);
    cpl_table_fill_column_window_int(kerntab, "pixmin", 0, nkern, 0);
    cpl_table_fill_column_window_int(kerntab, "pixmax", 0, nkern, 0);
    cpl_table_fill_column_window_int(kerntab, "pix0", 0, nkern, 0);
    cpl_table_fill_column_window_int(kerntab, "npix", 0, nkern, 0);
    cpl_table_fill_column_window_array(kerntab, "kernel", 0, nkern, kernel);

    /* Get pointers to table columns */
    pmin = cpl_table_get_data_int(kerntab, "pixmin");
    pmax = cpl_table_get_data_int(kerntab, "pixmax");
    p0 = cpl_table_get_data_int(kerntab, "pix0");
    np = cpl_table_get_data_int(kerntab, "npix");
    kernarr = cpl_table_get_data_array(kerntab, "kernel");

    /* Fill output table for constant kernel and return */
    if (varkern != 1) {
        pmin[0] = - floor(0.5 * (npmax - 1));
        pmax[0] = m - 1 + floor(0.5 * npmax);
        p0[0] = p0c;
        np[0] = npmax;
        cpl_array_copy_data_double(kernarr[0], kern0);
        cpl_array_delete(kernel0);
        cpl_array_delete(kernel);
        cpl_array_delete(tabkernel);
        return CPL_ERROR_NONE;
    }

    /* Calculate kernels for different pixel ranges */

    for (j = 0; j < nkern; j++) {

        /* Get central, lower, and upper wavelength of range */
        clam = x0 * pow(c, j-n1r);
        llam = clam / sqrt(c);
        ulam = clam * sqrt(c);

        /* Convert wavelengths into pixels by assuming a constant step size
           of the wavelength grid */
        pminj = p0c + floor((llam - x0) / dlam + 0.5);
        pmaxj = p0c + floor((ulam - x0) / dlam + 0.5);

        /* Skip kernel if it is valid for the same pixels as the previous
           loop */
        if (j > 0) {
            if (pminj == pmax[jo]) pminj++;
            if (pminj > pmaxj) continue;
            if (pmaxj == pmax[jo]) continue;
        }

        /* Write pixel ranges into output table */
        pmin[j] = pminj;
        pmax[j] = pmaxj;
        p0[j] = floor(0.5 * (pminj + pmaxj));
        jo = j;

        /* Scale FWHM of kernel components */
        cwbox = wbox * clam / reflam;
        cwgauss = wgauss * clam / reflam;
        cwlorentz = wlorentz * clam / reflam;

        /* Calculate synthetic kernel */
        mf_modsim_calcsynthkernel(kernel, cwbox, cwgauss, cwlorentz, kernfac,
                                  kernmode);

        /* Get kernel size and check whether it is not too large for the array
           column in the output table */
        npj = cpl_array_get_size(kernel);
        if (npj > npmax) {
            cpl_array_delete(kernel0);
            cpl_array_delete(kernel);
            cpl_array_delete(tabkernel);
            sprintf(errtxt, "%s: cpl_array *kernel (size too large for "
                    "column 'kernel' in cpl_table *kerntab)",
                    MF_ERROR_IOS_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }

        /* Write kernel and its size into output table */
        np[j] = npj;
        kern = cpl_array_get_data_double(kernel);
        for (k = 0; k < npj; k++) {
            tabkern[k] = kern[k];
        }
        cpl_array_copy_data_double(kernarr[j], tabkern);

    }

    /* Remove possible empty rows in output table */
    cpl_table_select_all(kerntab);
    cpl_table_and_selected_int(kerntab, "npix", CPL_EQUAL_TO, 0);
    cpl_table_erase_selected(kerntab);

    /* Free memory */
    cpl_array_delete(kernel0);
    cpl_array_delete(kernel);
    cpl_array_delete(tabkernel);

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_calcsynthkernel(cpl_array *kernel,
                                         const double wbox,
                                         const double wgauss,
                                         const double wlorentz,
                                         const double kernfac,
                                         const int kernmode)
{
    /*!
     * Calculates a synthetic kernel consisting of a boxcar, a Gaussian, and a
     * Lorentzian of given FWHM. The convolution of the latter two components
     * can also be substituted by a Voigt profile approximation (parameter
     * \e kernmode). The size of the Gaussian, Lorentzian, or Voigtian kernels
     * are determined by \e kernfac, which gives the size in units of FWHM.
     * The sum of the output kernel values is normalised to 1.
     *
     * \b INPUT:
     * \param wbox      FWHM of boxcar in pixels
     * \param wgauss    FWHM of Gaussian in pixels
     * \param wlorentz  FWHM of Lorentzian in pixels
     * \param kernfac   size of Gaussian/Lorentzian kernel in FWHM
     * \param kernmode  kernel mode (approx. for Voigtian profile = 1,
     *                  separate Gaussian and Lorentzian kernels = 0)
     *
     * \b OUTPUT:
     * \param kernel    CPL array containing the kernel elements
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     */

    cpl_array *bkernel = NULL, *gkernel = NULL, *lkernel = NULL;
    cpl_array *vkernel = NULL;
    char errtxt[MF_MAXLEN];

    /* Invalid FWHM */
    if (wbox < 0. || wgauss < 0. || wlorentz < 0.) {
        kernel = NULL;
        sprintf(errtxt, "%s: wbox < 0 || wgauss < 0 || wlorentz < 0",
                MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Invalid kernel size */
    if (kernfac < 0.) {
        kernel = NULL;
        sprintf(errtxt, "%s: kernfac < 0", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Invalid kernel mode */
    if (kernmode < 0 || kernmode > 1) {
        kernel = NULL;
        sprintf(errtxt, "%s: kernmode != 0 or 1", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Create temporary kernel arrays */
    bkernel = cpl_array_new(0, CPL_TYPE_DOUBLE);
    gkernel = cpl_array_new(0, CPL_TYPE_DOUBLE);
    lkernel = cpl_array_new(0, CPL_TYPE_DOUBLE);
    vkernel = cpl_array_new(0, CPL_TYPE_DOUBLE);

    /* Calculate boxcar kernel */
    mf_modsim_calcboxkernel(bkernel, wbox);

    /* Independent Gaussian and Lorentzian kernels or Voigt profile
       approximation? */

    if (kernmode == 1) {

        /* Calculate Voigt profile kernel */
        mf_modsim_calcvoigtkernel(vkernel, wgauss, wlorentz, kernfac);

    } else {

        /* Calculate Gaussian kernel */
        mf_modsim_calcgausskernel(gkernel, wgauss, kernfac);

        /* Calculate Lorentzian kernel */
        mf_modsim_calclorentzkernel(lkernel, wlorentz, kernfac);

        /* Convolve Gaussian with Lorentzian kernel to get Voigtian kernel */
        mf_basic_convolvekernels(vkernel, gkernel, lkernel);

    }

    /* Convolve boxcar with Voigtian kernel to get output kernel */
    mf_basic_convolvekernels(kernel, bkernel, vkernel);

    /* Free memory */
    cpl_array_delete(bkernel);
    cpl_array_delete(gkernel);
    cpl_array_delete(lkernel);
    cpl_array_delete(vkernel);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_calcboxkernel(cpl_array *kernel, const double fwhm)
{
    /*!
     * Calculates boxcar kernel.
     * The sum of the kernel values is normalised to 1.
     *
     * \b INPUT:
     * \param fwhm    width of boxcar in pixels
     *
     * \b OUTPUT:
     * \param kernel  CPL array containing the kernel elements
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     */

    char errtxt[MF_MAXLEN];
    int nkpix = 0;
    double term = 0., mval = 0., fval = 0., k = 0.;

    /* Invalid FWHM */
    if (fwhm < 0.) {
        kernel = NULL;
        sprintf(errtxt, "%s: fwhm < 0", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Number of kernel pixels */
    term = (fwhm - 1.) / 2.;
    nkpix = 2 * ceil(term) + 1;
    cpl_array_set_size(kernel, nkpix);

    /* Kernel with one pixel only */
    if (nkpix == 1) {
        cpl_array_set(kernel, 0, 1.);
        return CPL_ERROR_NONE;
    }

    /* Set kernel values */

    mval = fmod(term, 1.) / fwhm;
    fval = 1. / fwhm;

    for (k = 0; k < nkpix; k++) {
        if (mval > 0 && (k == 0 || k == nkpix-1)) {
            cpl_array_set(kernel, k, mval);
        } else {
            cpl_array_set(kernel, k, fval);
        }
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_calcvoigtkernel(cpl_array *kernel,
                                         const double wgauss,
                                         const double wlorentz,
                                         const double kernfac)
{
    /*!
     * Calculates kernel based on a Voigt profile approximation.
     * The sum of the kernel values is normalised to 1.
     * The number of kernel pixels is determined by the product of the Voigt
     * FWHM derived from the input FWHM of Gaussian and Lorentzian and \e
     * kernfac (size of kernel in units of FWHM). Then the upper odd integer
     * number is taken.
     *
     * \b INPUT:
     * \param wgauss    FWHM of Gaussian in pixels
     * \param wlorentz  FWHM of Lorentzian in pixels
     * \param kernfac   size of kernel in FWHM
     *
     * \b OUTPUT:
     * \param kernel    CPL array containing the kernel elements
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     */

    char errtxt[MF_MAXLEN];
    int nkpix = 0, k = 0, nx = 0, refpix = 0, i = 0;
    double *kern;
    double gam = 0, wvoigt = 0., wlwv = 0., xmax = 0, xmin = 0, dx = 0;
    double x = 0, xv = 0., xv2 = 0., xv225 = 0., sum = 0;

    /* Invalid FWHM */
    if (wgauss < 0. || wlorentz < 0.) {
        kernel = NULL;
        sprintf(errtxt, "%s: wgauss < 0 || wlorentz < 0", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Invalid kernel size */
    if (kernfac < 0.) {
        kernel = NULL;
        sprintf(errtxt, "%s: kernfac < 0", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* gamma of Lorentzian */
    gam = wlorentz / 2;

    /* FWHM of Voigt profile in pixels */
    wvoigt = gam + sqrt(gam * gam + wgauss * wgauss);

    /* Ratio of Lorentzian and Voigt profile FWHM */
    wlwv = wlorentz / wvoigt;

    /* Number of kernel pixels */
    nkpix = 2 * ceil(wvoigt * kernfac / 2 - 0.5) + 1;
    cpl_array_set_size(kernel, nkpix);

    /* Kernel with one pixel only */
    if (nkpix == 1) {
        cpl_array_set(kernel, 0, 1.);
        return CPL_ERROR_NONE;
    }

    /* Get pointer to CPL array */
    cpl_array_fill_window_double(kernel, 0, nkpix, 0.);
    kern = cpl_array_get_data_double(kernel);

    /* Integration limits for upper bin */
    xmax = 0.5 * nkpix;
    xmin = xmax - 1;

    /* Number of points per pixel for integration of Voigt profile */
    nx = ceil(MF_BINS_PER_FWHM / wvoigt);

    /* Step in pixels for integration of Voigt profile */
    dx = 1. / nx;

    /* Reference pixel for mirroring of kernel values */
    refpix = floor(nkpix / 2.);

    /* Calculate kernel up to reference pixel */

    for (k = nkpix - 1; k >= refpix; k--) {

        if (xmax <= 0.) {

            /* Skip integration */
            kern[k] = 1.;

        } else {

            /* First point in pixels for integration */
            x = xmin + dx / 2;

            /* Perform integration */

            kern[k] = 0.;

            for (i = 0; i < nx; i++) {

                /* Get variables of approximation formula */
                xv = x / wvoigt;
                xv2 = xv * xv;
                xv225 = pow(fabs(xv), 2.25);

                /* Calculate Voigt approximation for integration point */
                kern[k] += (1. - wlwv) * exp(-2.772 * xv2)
                         + wlwv / (1. + 4. * xv2)
                         + 0.016 * (1. - wlwv) * wlwv
                         * (exp(-0.4 * xv225) - 10. / (10. + xv225));

                /* Get next integration point */
                x += dx;

            }

            kern[k] /= (double) nx;

        }

        /* Shift integration limits for next bin */
        xmax = xmin;
        xmin = xmax - 1;

    }

    /* Mirror right wing of kernel */
    for (k = refpix - 1; k >= 0; k--) {
        kern[k] = kern[nkpix - k - 1];
    }

    /* Add all kernel values */
    for (k = 0; k < nkpix; k++) {
        if (kern[k] < 0.) {
            kern[k] = 0.;
        }
        sum += kern[k];
    }

    /* Normalise kernel values */
    for (k = 0; k < nkpix; k++) {
        kern[k] /= sum;
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_calcgausskernel(cpl_array *kernel, const double fwhm,
                                         const double kernfac)
{
    /*!
     * Calculates Gaussian kernel.
     * The sum of the kernel values is normalised to 1.
     * The number of kernel pixels is determined by the product of the input
     * parameters FWHM and \e kernfac (size of kernel in units of FWHM).
     * Then the upper odd integer number is taken.
     *
     * \b INPUT:
     * \param fwhm     FWHM of Gaussian in pixels
     * \param kernfac  size of kernel in FWHM
     *
     * \b OUTPUT:
     * \param kernel   CPL array containing the kernel elements
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     */

    char errtxt[MF_MAXLEN];
    int nkpix = 0, k = 0, nx = 0, refpix = 0, i = 0;
    double *kern;
    double sigma = 0, xmax = 0, xmin = 0, dx = 0, x = 0, sum = 0;

    /* Invalid FWHM */
    if (fwhm < 0.) {
        kernel = NULL;
        sprintf(errtxt, "%s: fwhm < 0", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Invalid kernel size */
    if (kernfac < 0.) {
        kernel = NULL;
        sprintf(errtxt, "%s: kernfac < 0", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* sigma of Gaussian */
    sigma = fwhm / CPL_MATH_FWHM_SIG;

    /* Number of kernel pixels */
    nkpix = 2 * ceil(fwhm * kernfac / 2 - 0.5) + 1;
    cpl_array_set_size(kernel, nkpix);

    /* Kernel with one pixel only */
    if (nkpix == 1) {
        cpl_array_set(kernel, 0, 1.);
        return CPL_ERROR_NONE;
    }

    /* Get pointer to CPL array */
    cpl_array_fill_window_double(kernel, 0, nkpix, 0.);
    kern = cpl_array_get_data_double(kernel);

    /* Integration limits for upper bin */
    xmax = 0.5 * nkpix;
    xmin = xmax - 1;

    /* Number of points per pixel for integration of Gaussian */
    nx = ceil(MF_BINS_PER_FWHM / fwhm);

    /* Step in pixels for integration of Gaussian */
    dx = 1. / nx;

    /* Reference pixel for mirroring of kernel values */
    refpix = floor(nkpix / 2.);

    /* Calculate kernel up to reference pixel */

    for (k = nkpix - 1; k >= refpix; k--) {

        if (xmax <= 0.) {

            /* Skip integration */
            kern[k] = 1.;

        } else {

            /* First point in pixels for integration */
            x = xmin + dx / 2;

            /* Perform integration */

            kern[k] = 0.;

            for (i = 0; i < nx; i++) {
                kern[k] += exp(-0.5 * pow(x / sigma, 2));
                x += dx;
            }

            kern[k] /= (double) nx;

        }

        /* Shift integration limits for next bin */
        xmax = xmin;
        xmin = xmax - 1;

    }

    /* Mirror right wing of kernel */
    for (k = refpix - 1; k >= 0; k--) {
        kern[k] = kern[nkpix - k - 1];
    }


    /* Add all kernel values */
    for (k = 0; k < nkpix; k++) {
        if (kern[k] < 0.) {
            kern[k] = 0.;
        }
        sum += kern[k];
    }

    /* Normalise kernel values */
    for (k = 0; k < nkpix; k++) {
        kern[k] /= sum;
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_calclorentzkernel(cpl_array *kernel,
                                           const double fwhm,
                                           const double kernfac)
{
    /*!
     * Calculates Lorentzian kernel.
     * The sum of the kernel values is normalised to 1.
     * The number of kernel pixels is determined by the product of the input
     * parameters FWHM and \e kernfac (size of kernel in units of FWHM).
     * Then the upper odd integer number is taken.
     *
     * \b INPUT:
     * \param fwhm     FWHM of Lorentzian in pixels
     * \param kernfac  size of kernel in FWHM
     *
     * \b OUTPUT:
     * \param kernel   CPL array containing the kernel elements
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     */

    char errtxt[MF_MAXLEN];
    int nkpix = 0, k = 0, nx = 0, refpix = 0, i = 0;
    double *kern;
    double gam = 0, xmax = 0, xmin = 0, dx = 0, x = 0, sum = 0;

    /* Invalid FWHM */
    if (fwhm < 0.) {
        kernel = NULL;
        sprintf(errtxt, "%s: fwhm < 0", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Invalid kernel size */
    if (kernfac < 0.) {
        kernel = NULL;
        sprintf(errtxt, "%s: kernfac < 0", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* gamma of Lorentzian */
    gam = fwhm / 2;

    /* Number of kernel pixels */
    nkpix = 2 * ceil(fwhm * kernfac / 2 - 0.5) + 1;
    cpl_array_set_size(kernel, nkpix);

    /* Kernel with one pixel only */
    if (nkpix == 1) {
        cpl_array_set(kernel, 0, 1.);
        return CPL_ERROR_NONE;
    }

    /* Get pointer to CPL array */
    cpl_array_fill_window_double(kernel, 0, nkpix, 0.);
    kern = cpl_array_get_data_double(kernel);

    /* Integration limits for upper bin */
    xmax = 0.5 * nkpix;
    xmin = xmax - 1;

    /* Number of points per pixel for integration of Lorentzian */
    nx = ceil(MF_BINS_PER_FWHM / fwhm);

    /* Step in pixels for integration of Lorentzian */
    dx = 1. / nx;

    /* Reference pixel for mirroring of kernel values */
    refpix = floor(nkpix / 2.);

    /* Calculate kernel up to reference pixel */

    for (k = nkpix - 1; k >= refpix; k--) {

        if (xmax <= 0.) {

            /* Skip integration */
            kern[k] = 1.;

        } else {

            /* First point in pixels for integration */
            x = xmin + dx / 2;

            /* Perform integration */

            kern[k] = 0.;

            for (i = 0; i < nx; i++) {
                kern[k] += gam / (x * x + gam * gam);
                x += dx;
            }

            kern[k] /= (double) nx;

        }

        /* Shift integration limits for next bin */
        xmax = xmin;
        xmin = xmax - 1;

    }

    /* Mirror right wing of kernel */
    for (k = refpix - 1; k >= 0; k--) {
        kern[k] = kern[nkpix - k - 1];
    }


    /* Add all kernel values */
    for (k = 0; k < nkpix; k++) {
        if (kern[k] < 0.) {
            kern[k] = 0.;
        }
        sum += kern[k];
    }

    /* Normalise kernel values */
    for (k = 0; k < nkpix; k++) {
        kern[k] /= sum;
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_invertkerneltable(cpl_table *kerntab)
{
    /*!
     * Inverts kernels in a table that are valid for different pixel ranges.
     * The resulting kernels provide the contributions of the different
     * pixels to each output pixel. The applied approximation derives the
     * inverted kernels only for the central pixel of each range and uses
     * these kernels for the full range. This approach is consistent with the
     * kernel creation in ::mf_modsim_calckerneltable, which also calculates
     * the kernels only for the central range pixels. The resulting kernels
     * are composed of contributions from different input kernels if their
     * size is larger than the corresponding pixel ranges. Only in the latter
     * case, the output kernels differ from the input kernels and become
     * (slightly) asymmetric. The output kernels are normalised, i.e. the sum
     * of the kernel values equals 1.
     *
     * \b INPUT:
     * \param kerntab   table with kernels for different pixel ranges
     *
     * \b OUTPUT:
     * \param kerntab   table with inverted kernels
     *
     * \b ERRORS:
     * - Invalid object structure
     * - No data
     */

    cpl_array **kernarr = NULL, **tkernarr = NULL;
    char errtxt[MF_MAXLEN];
    int nkern = 0, npmax = 0, i = 0, kimax = 0, j = 0, kjcen = 0, kjmin = 0;
    int kjmax = 0, nk = 0, kimin = 0, k = 0;
    int *pmin, *pmax, *p0, *np, *tnp, *kcen;
    double sum = 0.;
    double *outkern, *inkern;

    /* Check existence of required columns in kernel table */
    if (cpl_table_has_column(kerntab, "pixmin") != 1 ||
        cpl_table_has_column(kerntab, "pixmax") != 1 ||
        cpl_table_has_column(kerntab, "pix0") != 1 ||
        cpl_table_has_column(kerntab, "npix") != 1 ||
        cpl_table_has_column(kerntab, "kernel") != 1) {
        sprintf(errtxt, "%s: cpl_table *kerntab (missing column(s))",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get number of kernels */
    nkern = cpl_table_get_nrow(kerntab);
    if (nkern == 0) {
        sprintf(errtxt, "%s: cpl_table *kerntab", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Get maximum number of kernel elements */
    npmax = cpl_table_get_column_depth(kerntab, "kernel");
    if (npmax == 0) {
        sprintf(errtxt, "%s: cpl_table *kerntab (empty kernel arrays)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Duplicate kernel size column */
    if (cpl_table_has_column(kerntab, "tnpix") == 1) {
        cpl_table_erase_column(kerntab, "tnpix");
    }
    cpl_table_duplicate_column(kerntab, "tnpix", kerntab, "npix");

    /* Duplicate kernel array column */
    if (cpl_table_has_column(kerntab, "tkernel") == 1) {
        cpl_table_erase_column(kerntab, "tkernel");
    }
    cpl_table_duplicate_column(kerntab, "tkernel", kerntab, "kernel");

    /* Create new column for central kernel pixel if required */
    if (cpl_table_has_column(kerntab, "kerncen") != 1) {
        cpl_table_new_column(kerntab, "kerncen", CPL_TYPE_INT);
    }
    cpl_table_fill_column_window_int(kerntab, "kerncen", 0, nkern, 0);

    /* Get pointers to kernel table columns */
    pmin = cpl_table_get_data_int(kerntab, "pixmin");
    pmax = cpl_table_get_data_int(kerntab, "pixmax");
    p0 = cpl_table_get_data_int(kerntab, "pix0");
    np = cpl_table_get_data_int(kerntab, "npix");
    tnp = cpl_table_get_data_int(kerntab, "tnpix");
    kernarr = cpl_table_get_data_array(kerntab, "kernel");
    tkernarr = cpl_table_get_data_array(kerntab, "tkernel");
    kcen = cpl_table_get_data_int(kerntab, "kerncen");

    /* Invert kernel(s) */

    for (i = 0; i < nkern; i++) {

        /* Get pointer to kernel in table column to be modified */
        outkern = cpl_array_get_data_double(kernarr[i]);

        for (kimax = -1, sum = 0., j = nkern-1; j >= 0; j--) {

            /* Derive contributing range of kernel elements */
            kjcen = floor(0.5 * tnp[j]);
            kjmin = p0[i] - pmax[j] + kjcen;
            kjmax = p0[i] - pmin[j] + kjcen;

            /* Skip non-contributing kernels */
            if (kjmax < 0 || kjmin >= tnp[j]) continue;

            /* Avoid invalid kernel elements */
            if (kjmin < 0) kjmin = 0;
            if (kjmax >= tnp[j]) kjmax = tnp[j] - 1;

            /* Make sure that output kernels are not cut due to too a small
               range of pixels valid for the first or last kernel in the
               table */
            if (j == 0 && kjmax < tnp[j]-1) kjmax = tnp[j] - 1;
            if (j == nkern-1 && kjmin > 0) kjmin = 0;

            /* Elements of output kernel to be filled */
            nk = kjmax - kjmin + 1;
            kimin = kimax + 1;
            kimax = kimin + nk - 1;

            /* Identify central element of output kernel */
            if (j == i) kcen[i] = kimin + kjcen - kjmin;

            /* Make sure that maximum kernel size is not exceeded */
            if (kimax >= npmax) {
                kimax = npmax - 1;
                nk = kimax - kimin + 1;
                if (nk < 1) break;
            }

            /* Get pointer to kernel in temporary table column */
            inkern = cpl_array_get_data_double(tkernarr[j]);

            /* Copy selected kernel elements from input to output kernel and
               sum up kernel values*/
            for (k = 0; k < nk; k++) {
                outkern[kimin + k] = inkern[kjmin + k];
                sum += outkern[kimin + k];
            }

        }

        /* Set size of inverted kernel */
        np[i] = kimax + 1;

        /* Normalise kernel values */
        for (k = 0; k < np[i]; k++) {
            outkern[k] /= sum;
        }

    }

    /* Erase temporary kernel array column */
    cpl_table_erase_column(kerntab, "tnpix");
    cpl_table_erase_column(kerntab, "tkernel");

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_convolvekerneltable(cpl_table *spec,
                                             const cpl_table *kerntab)
{
    /*!
     * Convolves the flux column of a spectrum with kernels from a table. Each
     * input kernel has to be valid for a fixed pixel range, which is also
     * provided by the kernel table. The convolution is performed in an
     * inverted way, i.e. the kernel gives the contributions of the different
     * pixels to each output pixel.
     *
     * \b INPUT:
     * \param spec      input spectrum (CPL table)
     * \param kerntab   table with kernels for different pixel ranges
     *
     * \b OUTPUT:
     * \param spec      input spectrum with convolved flux
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     */

    const cpl_array **kernarr = NULL;
    cpl_array *kernel = NULL, *flux = NULL, *convflux = NULL;
    char errtxt[MF_MAXLEN];
    int m = 0, nkern = 0, npmax = 0, j = 0, range[2] = {0, 0}, k = 0;
    const int *pmin, *pmax, *np, *kcen;
    const double *tabkern;
    double *vec, *kern;

    /* Check number of data points in spectrum */
    m = cpl_table_get_nrow(spec);
    if (m <= 0) {
        sprintf(errtxt, "%s: cpl_table *spec", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check existence of 'flux' column in spectrum table */
    if (cpl_table_has_column(spec, "flux") != 1) {
        sprintf(errtxt, "%s: cpl_table *spec (no 'flux' column)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of required columns in kernel table */
    if (cpl_table_has_column(kerntab, "pixmin") != 1 ||
        cpl_table_has_column(kerntab, "pixmax") != 1 ||
        cpl_table_has_column(kerntab, "npix") != 1 ||
        cpl_table_has_column(kerntab, "kernel") != 1) {
        sprintf(errtxt, "%s: cpl_table *kerntab (missing column(s))",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get number of kernels */
    nkern = cpl_table_get_nrow(kerntab);
    if (nkern == 0) {
        sprintf(errtxt, "%s: cpl_table *kerntab", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Get maximum number of kernel elements */
    npmax = cpl_table_get_column_depth(kerntab, "kernel");
    if (npmax == 0) {
        sprintf(errtxt, "%s: cpl_table *kerntab (empty kernel arrays)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get pointers to kernel table columns */
    pmin = cpl_table_get_data_int_const(kerntab, "pixmin");
    pmax = cpl_table_get_data_int_const(kerntab, "pixmax");
    np = cpl_table_get_data_int_const(kerntab, "npix");
    kernarr = cpl_table_get_data_array_const(kerntab, "kernel");
    kcen = cpl_table_get_data_int_const(kerntab, "kerncen");

    /* Initialise CPL array for kernel */
    kernel = cpl_array_new(npmax, CPL_TYPE_DOUBLE);
    cpl_array_fill_window_double(kernel, 0, npmax, 0.);

    /* Copy flux column into CPL array */
    vec = cpl_table_get_data_double(spec, "flux");
    flux = cpl_array_new(m, CPL_TYPE_DOUBLE);
    cpl_array_copy_data_double(flux, vec);

    /* Create CPL array for convolved flux */
    convflux = cpl_array_new(m, CPL_TYPE_DOUBLE);
    cpl_array_fill_window_double(convflux, 0, m, 0.);

    /* Convolve spectrum with kernel(s) */

    for (j = nkern - 1; j >= 0; j--) {

        /* Set pixel range */
        range[0] = pmin[j];
        range[1] = pmax[j];

        /* Make sure that the pixel ranges are valid */
        if (range[1] < 0 || range[0] >= m) continue;
        if (range[0] < 0) range[0] = 0;
        if (range[1] >= m) range[1] = m - 1;

        /* Adapt size of kernel array and get pointer */
        cpl_array_set_size(kernel, np[j]);
        kern = cpl_array_get_data_double(kernel);

        /* Get pointer to kernel in table */
        tabkern = cpl_array_get_data_double_const(kernarr[j]);

        /* Copy kernel elements into array of correct length */
        for (k = 0; k < np[j]; k++) {
            kern[k] = tabkern[k];
        }

        /* Convolve spectrum with inverted kernel */
        mf_basic_convolvewindow_inv(convflux, flux, range, kernel, kcen[j]);

    }

    /* Copy resulting flux array into "flux" column of input CPL table */
    vec = cpl_array_get_data_double(convflux);
    cpl_table_copy_data_double(spec, "flux", vec);

    /* Free memory */
    cpl_array_delete(kernel);
    cpl_array_delete(flux);
    cpl_array_delete(convflux);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_convolvereadkernel(cpl_table *spec,
                                            const cpl_array *specrows,
                                            const mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Applies kernels read from a file to a model spectrum as provided by a
     * CPL table with columns "lambda" and "flux". For each spectrum row, a
     * kernel is expected. They have to be provided in matrix form via the
     * ::mfdrv parameter structure. All kernels have to have the same size and
     * the kernel centre should always be at the expected median position.
     *
     * \b INPUT:
     * \param spec      input spectrum (CPL table)
     * \param specrows  CPL array with spectrum rows to be convolved
     * \param drvpar    ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec      convolved spectrum
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     */

    cpl_array *kernel = NULL, *flux = NULL, *convflux = NULL;
    char errtxt[MF_MAXLEN];
    int m = 0, nsel = 0, npmax = 0, nrow = 0, i = 0, range[2] = {0, 0}, k = 0;
    int j = 0;
    const int *rows;
    double *matkern, *kern, *vec;

    /* Check number of data points in spectrum */
    m = cpl_table_get_nrow(spec);
    if (m <= 0) {
        sprintf(errtxt, "%s: cpl_table *spec", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check existence of 'flux' column in spectrum table */
    if (cpl_table_has_column(spec, "flux") != 1) {
        sprintf(errtxt, "%s: cpl_table *spec (no 'flux' column)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check number of selected spectrum pixels */
    nsel = cpl_array_get_size(specrows);
    if (nsel != m) {
        sprintf(errtxt, "%s: cpl_array *specrows (inconsistent with size of "
                "cpl_table *spec)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get pointer to array of row numbers */
    rows = cpl_array_get_data_int_const(specrows);

    /* Check validity of row numbers */
    if (rows[0] < 0 || rows[nsel-1] < rows[0] ||
        rows[nsel-1] - rows[0] + 1 != nsel) {
        sprintf(errtxt, "%s: cpl_array *specrows (no list of increasing row "
                "numbers)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get pointer to matrix of kernel elements */
    matkern = cpl_matrix_get_data(drvpar->kernel);

    /* Get maximum kernel size */
    npmax = cpl_matrix_get_ncol(drvpar->kernel);

    /* Check for one element size of kernel matrix -> no convolution */
    nrow = cpl_matrix_get_nrow(drvpar->kernel);
    if (nrow == 1 && npmax == 1 && matkern[0] == 1.) {
        return CPL_ERROR_NONE;
    }

    /* Check number of rows in kernel matrix */
    if (nrow < nsel || nrow - 1 < rows[nsel-1]) {
        sprintf(errtxt, "%s: cpl_matrix drvpar->kernel (inconsistent with "
                "size of cpl_array *specrows)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Initialise CPL array for kernel and get pointer */
    kernel = cpl_array_new(npmax, CPL_TYPE_DOUBLE);
    cpl_array_fill_window_double(kernel, 0, npmax, 0.);
    kern = cpl_array_get_data_double(kernel);

    /* Copy flux column into CPL array */
    vec = cpl_table_get_data_double(spec, "flux");
    flux = cpl_array_new(m, CPL_TYPE_DOUBLE);
    cpl_array_copy_data_double(flux, vec);

    /* Create CPL array for convolved flux */
    convflux = cpl_array_new(m, CPL_TYPE_DOUBLE);
    cpl_array_fill_window_double(convflux, 0, m, 0.);

    /* Convolve spectrum with kernel(s) (pixel by pixel) */

    for (i = 0; i < m; i++) {

        /* Set pixel range */
        range[0] = i;
        range[1] = i;

        /* Copy kernel elements into array of correct length */
        for (k = 0; k < npmax; k++) {
            j = k + npmax * rows[i];
            kern[k] = matkern[j];
        }

        /* Convolve pixel with kernel */
        mf_basic_convolvewindow(convflux, flux, range, kernel);

    }

    /* Copy resulting flux array into "flux" column of input CPL table */
    vec = cpl_array_get_data_double(convflux);
    cpl_table_copy_data_double(spec, "flux", vec);

    /* Free memory */
    cpl_array_delete(kernel);
    cpl_array_delete(flux);
    cpl_array_delete(convflux);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_telback(cpl_table *spec, const mfdrv *drvpar,
                                 const cpl_array *fitpar)
{
    /*!
     * Computes thermal emission by telesope/instrument and adds it to the
     * input spectrum. The routine assumes a grey body depending on emissivity
     * (fit parameter telback) and ambient temperature (from ::mfdrv
     * structure).
     *
     * \note The grey body emission is divided by 1 - emissivity in order to
     * have the same flux level as the sky emission. This approach is correct
     * for flux-calibrated spectra, where the reflection of light by the
     * telescope mirror has been corrected, which is wrong for the telescope
     * emission. For this reason, the apparent telescope emission becomes
     * higher than the true one.
     *
     * \b INPUT:
     * \param spec    input spectrum (CPL table)
     * \param drvpar  ::mfdrv parameter structure
     * \param fitpar  CPL array with fit parameters
     *
     * \b OUTPUT:
     * \param spec    input spectrum + telescope/instrument emission
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    cpl_parameter *p;
    char errtxt[MF_MAXLEN];
    int trans = 0, ntel = 0, nlam = 0., i = 0;
    double c1 = 0., c2 = 0., c3 = 0., c4 = 0., c5 = 0.;
    double temp = 0., eps = 0.;
    double *lam, *flux;

    /* Return immediately in the case of a transmission curve */
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);
    if (trans != 0) {
        return CPL_ERROR_NONE;
    }

    /* Constants and unit conversions */
    c1 = 2 * CPL_PHYS_H * pow(CPL_PHYS_C, 2);
    c2 = CPL_PHYS_H * CPL_PHYS_C / CPL_PHYS_K;
    c3 = MF_LAM_UNIT / (CPL_PHYS_H * CPL_PHYS_C * MF_SR_IN_ARCSEC2);
    c4 = c1 * c3 / pow(MF_LAM_UNIT, 4);
    c5 = c2 / MF_LAM_UNIT;

    /* Get primary mirror temperature */
    p = cpl_parameterlist_find(drvpar->parlist, "m1temp");
    temp = cpl_parameter_get_double(p) + 273.15;
    if (temp < 0.) {
        sprintf(errtxt, "%s: m1temp of mfdrv *drvpar < 0 K",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Get emissivity estimate from CPL array of fit parameters */
    ntel = cpl_array_get_size(fitpar) - 1;
    eps = cpl_array_get(fitpar, ntel, NULL);
    if (eps < 0 || eps >= 1.) {
        sprintf(errtxt, "%s: emissivity of cpl_array *fitpar < 0 or >= 1",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Get pointers to CPL table columns */
    lam = cpl_table_get_data_double(spec, "lambda");
    flux = cpl_table_get_data_double(spec, "flux");

    /* Compute and add grey body emission to spectrum */
    nlam = cpl_table_get_nrow(spec);
    for (i = 0; i < nlam; i++) {
        flux[i] += (eps / (1 - eps)) * c4 /
                   (pow(lam[i], 4) * (expm1(c5 / (lam[i] * temp))));
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_fluxunits(cpl_table *spec, const mfdrv *drvpar)
{
    /*!
     * Conversion of fluxes from phot/(s*m2*mum*as2) (emission spectrum only)
     * to flux unit of observed spectrum.
     * Meaning of the flag flux_unit of the ::mfdrv parameter structure:
     * 0: phot/(s*m^2*mum*as^2) [no conversion]
     * 1: W/(m^2*mum*as^2)
     * 2: erg/(s*cm^2*A*as^2)
     * 3: mJy/as^2
     * For other units, the conversion factor has to be considered as constant
     * term of the continuum fit parameters.
     *
     * \b INPUT:
     * \param spec    CPL table with model spectrum
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    model spectrum with adapted flux units
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    cpl_parameter *p;
    char errtxt[MF_MAXLEN];
    int trans = 0, unit = 0;

    /* Return immediately in the case of a transmission curve */
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);
    if (trans != 0) {
        return CPL_ERROR_NONE;
    }

    /* Get flux unit flag */
    p = cpl_parameterlist_find(drvpar->parlist, "flux_unit");
    unit = cpl_parameter_get_int(p);
    if (unit < 0 || unit > 3) {
        sprintf(errtxt, "%s: flux_unit of mfdrv *drvpar < 0 or > 3",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Convert fluxes */

    if (unit == 1) {
        /* phot/sm2mum -> W/m2mum */
        cpl_table_divide_columns(spec, "flux", "lambda");
        cpl_table_multiply_scalar(spec, "flux",
                                  CPL_PHYS_C * CPL_PHYS_H / MF_LAM_UNIT);
    } else if (unit == 2) {
        /* phot/sm2mum -> erg/scm2A */
        cpl_table_divide_columns(spec, "flux", "lambda");
        cpl_table_multiply_scalar(spec, "flux",
                                 0.1 * CPL_PHYS_C * CPL_PHYS_H / MF_LAM_UNIT);
    } else if (unit == 3) {
        /* phot/sm2mum -> mJy */
        cpl_table_multiply_columns(spec, "flux", "lambda");
        cpl_table_multiply_scalar(spec, "flux",
                                  1e35 * MF_LAM_UNIT * CPL_PHYS_H);
    }

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_modcont(cpl_table *spec, const mfdrv *drvpar,
                                 const cpl_array *fitpar, const int range)
{
    /*!
     * Modifies continuum of a model spectrum by means of polynomials.
     * The coefficients are provided by the fit parameter vector (CPL array).
     * The continuum correction is carried out for the given fit range only.
     *
     * \b INPUT:
     * \param spec    CPL table with model spectrum
     * \param drvpar  ::mfdrv parameter structure
     * \param fitpar  CPL array with fit parameters
     * \param range   range number
     *
     * \b OUTPUT:
     * \param spec    model spectrum with modified continuum
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    cpl_parameter *p;
    char errtxt[MF_MAXLEN];
    int nmolec = 0, nchip = 0, nwlc = 0, ncont = 0, n0 = 0, nlam = 0;
    int i = 0, j = 0;
    double limlam[2] = {0., 0.}, wmean = 0., fac = 0.;
    const double *par;
    double *slam, *scal;

    /* Find position of wavelength calibration parameters in fit parameter
       CPL array */
    p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "nchip");
    nchip = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "wlc_n");
    nwlc = cpl_parameter_get_int(p) + 1;
    p = cpl_parameterlist_find(drvpar->parlist, "cont_n");
    ncont = cpl_parameter_get_int(p) + 1;
    n0 = nmolec + nwlc * nchip + ncont * (range - 1);
    if (ncont-1 < 0) {
        sprintf(errtxt, "%s: cont_n of mfdrv *drvpar < 0", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Get pointer to CPL array with fit parameters */
    par = cpl_array_get_data_double_const(fitpar);

    /* Shift zero point of wavelength scale to centre of wavelength range */
    nlam = cpl_table_get_nrow(spec);
    limlam[0] = cpl_table_get(spec, "lambda", 0, NULL);
    limlam[1] = cpl_table_get(spec, "lambda", nlam-1, NULL);
    wmean = (limlam[0] + limlam[1]) / 2;
    cpl_table_duplicate_column(spec, "slambda", spec, "lambda");
    cpl_table_subtract_scalar(spec, "slambda", wmean);

    /* Get pointers to CPL table columns */
    slam = cpl_table_get_data_double(spec, "slambda");
    scal = cpl_table_get_data_double(spec, "scal");

    /* Compute continuum correction polynomials */
    for (i = 0; i < nlam; i++) {
        for (fac = 0, j = 0; j < ncont; j++) {
            fac += par[n0 + j] * pow(slam[i], j);
        }
        if (fac <= 0) {
            scal[i] = 1.;
        } else {
            scal[i] = fac;
        }
    }

    /* Remove temporary table column */
    cpl_table_erase_column(spec, "slambda");

    /* Scale continuum */
    cpl_table_multiply_columns(spec, "flux", "scal");

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_modsim_corrobsspec(cpl_table *spec, const mfdrv *drvpar)
{
    /*!
     * Scales a model spectrum by means of a linear interpolation of the
     * wavelength windows of the observed spectrum which are not affected by
     * significant molecular absorption. No scaling is performed in the case
     * of molecular emission spectra.
     *
     * \b INPUT:
     * \param spec    CPL table with model spectrum
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    model spectrum with scaled continuum
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    cpl_parameter *p;
    cpl_table *contspec;
    //char errtxt[MF_MAXLEN];
    int trans = 0, nlam = 0, count = 0, i = 0;
    int *iscont;
    double *flux;
    double transmin = 0.93;
    double pixfrac = 0.2;

    /* Return immediately in the case of an atmospheric emission spectrum */
    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);
    if (trans == 0) {
        return CPL_ERROR_NONE;
    }

    /* Create column for interpolated flux */
    nlam = cpl_table_get_nrow(spec);
    cpl_table_duplicate_column(spec, "ipoflux", spec, "oflux");

    /* Get pointers to CPL table columns */
    flux = cpl_table_get_data_double(spec, "flux");
    iscont = cpl_table_get_data_int(spec, "iscont");

    if (iscont[0] == -1) {

        /* Find continuum pixels */
        for (count = 0, i = 0; i < nlam; i++) {
            if (flux[i] < transmin) {
                iscont[i] = 0;
            } else {
                iscont[i] = 1;
                count++;
            }
        }
        printf("count = %d\n", count);

        /* Require minimum number of continuum pixels */
        if ((double) count / ((double) nlam) < pixfrac) {
            cpl_msg_warning(cpl_func, "Not enough pixels for continuum "
                            "interpolation in observed spectrum");
            return CPL_ERROR_NONE;
        }

    }

    /* Extract data for continuum pixels */
    cpl_table_unselect_all(spec);
    cpl_table_or_selected_int(spec, "iscont", CPL_EQUAL_TO, 1);
    contspec = cpl_table_extract_selected(spec);
    cpl_table_select_all(spec);

    /* Interpolate object spectrum between the identified continuum pixels */
    mf_basic_interpolcolumn(contspec, "lambda", "oflux", spec, "lambda",
                            "ipoflux", 0, nlam);

    /* Delete temporary table */
    cpl_table_delete(contspec);

    /* Modify scaling function */
    cpl_table_multiply_columns(spec, "scal", "ipoflux");

    /* Correct model spectrum by interpolated observed spectrum */
    //cpl_table_multiply_columns(spec, "flux", "ipoflux");

    /* Remove temporary table columns */
    cpl_table_erase_column(spec, "ipoflux");

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    return CPL_ERROR_NONE;
}

/**@}*/
