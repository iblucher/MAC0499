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
 * \file mf_par.c
 *
 * Routines for reading and writing the MOLECFIT driver file
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  27 Apr 2010
 * \date   03 Nov 2014
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#define _GNU_SOURCE
#include <regex.h>
#include <unistd.h>
#include <mf_par.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

/*!
 * Utility macro that calls ::mf_par_search.
 * Jumps to parse_end on parse failure.
 */

#define MF_PAR_SEARCH_CHECK(par, npar, parname, parfile) \
    do { \
        mf_par_search(par, npar, parname, parfile); \
        if (cpl_error_get_code() != CPL_ERROR_NONE) { \
            goto parse_end; \
        } \
    } while (0)


/*!
 * Utility macro that calls ::mf_par_search for npar == 2 parameters and fills
 * the result into the parameter with the name \e parname of type TYPE.
 * \e code needs to be the name of the member of the ::mfpar structure.
 * Jumps to parse_end on parse failure.
 */

#define MF_PAR_SEARCH_TYPE(TYPE, code, parname, parfile, parlist) \
    do { \
        int _npar; \
        mfpar _par[MF_MAXPAR]; \
        MF_PAR_SEARCH_CHECK(_par, &_npar, parname, parfile); \
        if (_par[0].c[0] != '#' && _npar == 2) { \
            cpl_parameter * _p = \
                cpl_parameterlist_find(parlist, parname); \
            cpl_parameter_set_##TYPE(_p, _par[1].code); \
        } \
        else if (_npar > 2) { /* zero means use default */ \
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, \
                  "Found %d parameters when parsing for %s, expected 1", \
                  _npar - 1, parname); \
                goto parse_end; \
        } \
    } while (0)


/*!
 * Macro to search and set one string parameter from \e parfile into
 * \e parlist.
 * Jumps to parse_end on parse failure.
 */

#define MF_PAR_SEARCH_S(parname, parfile, parlist) \
    MF_PAR_SEARCH_TYPE(string, c, parname, parfile, parlist)


/*!
 * Macro to search and set one double parameter from \e parfile into
 * \e parlist.
 * Jumps to parse_end on parse failure.
 */

#define MF_PAR_SEARCH_D(parname, parfile, parlist) \
    MF_PAR_SEARCH_TYPE(double, d, parname, parfile, parlist)


/*!
 * Macro to search and set one integer parameter from \e parfile into
 * \e parlist.
 * Jumps to parse_end on parse failure.
 */

#define MF_PAR_SEARCH_I(parname, parfile, parlist) \
    MF_PAR_SEARCH_TYPE(int, i, parname, parfile, parlist)

cpl_error_code mf_par_regex_number(int *npar, const char *parname,
                             const char *parfile);

cpl_error_code mf_par_initall(mfdrv *drvpar)
{
    /*!
     * Initialises the structure for the parameters of the MOLECFIT driver
     * file and fills this structure with default values.
     *
     * \b INPUT:
     * \param drvpar   empty ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param drvpar   ::mfdrv parameter structure filled with default values
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    char mol[4] = "H2O";

    /* Initialise parameter list */

    drvpar->parlist = cpl_parameterlist_new();

    /* Directory parameters */

    /* get current working directory */
    char cwd[MF_MAXLEN] = "", *d = NULL;
    if ((d = getcwd(cwd, sizeof(cwd))) == NULL) {
        return cpl_error_set_message(cpl_func, MF_ERROR_GETCWD,
                                     "%s: %s", MF_ERROR_GETCWD_TXT, cwd);
    }

    p = cpl_parameter_new_value("curdir", CPL_TYPE_STRING, "", "", cwd);
    cpl_parameterlist_append(drvpar->parlist, p);

    p = cpl_parameter_new_value("relbasedir", CPL_TYPE_STRING, "", "", ".");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("basedir", CPL_TYPE_STRING, "", "", ".");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("bindir", CPL_TYPE_STRING, "", "", "bin");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("configdir", CPL_TYPE_STRING, "", "",
                                "config");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("datadir", CPL_TYPE_STRING, "", "", "data");
    cpl_parameterlist_append(drvpar->parlist, p);
    /* only used by GUI, needed to write it out again in expert mode */
    p = cpl_parameter_new_value("user_workdir", CPL_TYPE_STRING, "", "", ".");
    cpl_parameterlist_append(drvpar->parlist, p);

    /* Input file parameters */

    p = cpl_parameter_new_value("filename", CPL_TYPE_STRING, "", "", "none");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("listname", CPL_TYPE_STRING, "", "", "none");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("trans", CPL_TYPE_INT, "", "", 1, 0, 2);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("col_lam", CPL_TYPE_STRING, "", "", "undef");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("col_flux", CPL_TYPE_STRING, "", "", "undef");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("col_dflux", CPL_TYPE_STRING, "", "",
                                "undef");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("col_mask", CPL_TYPE_STRING, "", "", "undef");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("default_error", CPL_TYPE_DOUBLE, "", "",
                                0.01);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("wlgtomicron", CPL_TYPE_DOUBLE, "", "", 1.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("vac_air", CPL_TYPE_STRING, "", "", "vac");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("wrange_include", CPL_TYPE_STRING, "", "",
                                "none");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("wrange_exclude", CPL_TYPE_STRING, "", "",
                                "none");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("prange_exclude", CPL_TYPE_STRING, "", "",
                                "none");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("npix", CPL_TYPE_INT, "", "", 0);
    cpl_parameterlist_append(drvpar->parlist, p);

    /* Output file parameters */

    p = cpl_parameter_new_value("output_dir", CPL_TYPE_STRING, "", "",
                                "output");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("output_name", CPL_TYPE_STRING, "", "",
                                "none");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("plot_creation", CPL_TYPE_STRING, "", "",
                                "none");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("plot_range", CPL_TYPE_INT, "", "", 0, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("writeascii", CPL_TYPE_INT, "", "", 1, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("clean_mflux", CPL_TYPE_INT, "", "", 0, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("writerespar", CPL_TYPE_INT, "", "", 0, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);

    /* Parameters for fitting (molecules, continuum, wavelength, resol.) */

    p = cpl_parameter_new_value("ftol", CPL_TYPE_DOUBLE, "", "", 1e-10);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("xtol", CPL_TYPE_DOUBLE, "", "", 1e-10);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("nmolec", CPL_TYPE_INT, "", "", 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("lbl_molecs", CPL_TYPE_STRING, "", "",
                                "00000000000000000000000000000000000000000000000");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("lnfl", CPL_TYPE_STRING, "", "", "");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("flux_unit", CPL_TYPE_INT, "", "", 0);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("fit_back", CPL_TYPE_INT, "", "", 1, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("telback", CPL_TYPE_DOUBLE, "", "", 0.1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("nrange", CPL_TYPE_INT, "", "", 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    // p = cpl_parameter_new_value("nrange_max", CPL_TYPE_INT, "", "", 8);
    // cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("fit_cont", CPL_TYPE_INT, "", "", 1, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("cont_n", CPL_TYPE_INT, "", "", 0);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("cont_const", CPL_TYPE_DOUBLE, "", "", 1.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("nchip", CPL_TYPE_INT, "", "", 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    // p = cpl_parameter_new_value("nchip_max", CPL_TYPE_INT, "", "", 4);
    // cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("fit_wlc", CPL_TYPE_INT, "", "", 1, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("fit_wlc_lin", CPL_TYPE_INT, "", "", 1, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("wlc_n", CPL_TYPE_INT, "", "", 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("wlc_const", CPL_TYPE_DOUBLE, "", "", 0.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("fit_res_box", CPL_TYPE_INT, "", "", 1, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("relres_box", CPL_TYPE_DOUBLE, "", "", 1.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("kernmode", CPL_TYPE_INT, "", "", 0, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("fit_res_gauss", CPL_TYPE_INT, "", "", 1, 0,
                                1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("res_gauss", CPL_TYPE_DOUBLE, "", "", 1.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("fit_res_lorentz", CPL_TYPE_INT, "", "", 1, 0,
                                1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("res_lorentz", CPL_TYPE_DOUBLE, "", "", 1.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("kernfac", CPL_TYPE_DOUBLE, "", "", 3.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("varkern", CPL_TYPE_INT, "", "", 0, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("kernel_file", CPL_TYPE_STRING, "", "",
                                "none");
    cpl_parameterlist_append(drvpar->parlist, p);

    /* Parameters for fit results */

    p = cpl_parameter_new_range("lastcall", CPL_TYPE_INT, "", "", 0, -1, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("chi2", CPL_TYPE_DOUBLE, "", "", 1e12);
    cpl_parameterlist_append(drvpar->parlist, p);

    /* Parameters for ambient conditions */

    p = cpl_parameter_new_value("obsdate", CPL_TYPE_DOUBLE, "", "", -1.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("obsdate_key", CPL_TYPE_STRING, "", "",
                                "MJD-OBS");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("utc", CPL_TYPE_DOUBLE, "", "", -1.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("utc_key", CPL_TYPE_STRING, "", "", "UTC");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("telalt", CPL_TYPE_DOUBLE, "", "", 90.,
                                -90., 90.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("telalt_key", CPL_TYPE_STRING, "", "",
                                "ESO TEL ALT");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("rhum", CPL_TYPE_DOUBLE, "", "", 15.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("rhum_key", CPL_TYPE_STRING, "", "",
                                "ESO TEL AMBI RHUM");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("pres", CPL_TYPE_DOUBLE, "", "", 750.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("pres_key", CPL_TYPE_STRING, "", "",
                                "ESO TEL AMBI PRES START");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("temp", CPL_TYPE_DOUBLE, "", "", 15.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("temp_key", CPL_TYPE_STRING, "", "",
                                "ESO TEL AMBI TEMP");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("m1temp", CPL_TYPE_DOUBLE, "", "", 15.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("m1temp_key", CPL_TYPE_STRING, "", "",
                                "ESO TEL TH M1 TEMP");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("geoelev", CPL_TYPE_DOUBLE, "", "", 2635.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("geoelev_key", CPL_TYPE_STRING, "", "",
                                "ESO TEL GEOELEV");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("longitude", CPL_TYPE_DOUBLE, "", "",
                                -70.4051);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("longitude_key", CPL_TYPE_STRING, "", "",
                                "ESO TEL GEOLON");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("latitude", CPL_TYPE_DOUBLE, "", "",
                                -24.6276);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("latitude_key", CPL_TYPE_STRING, "", "",
                                "ESO TEL GEOLAT");
    cpl_parameterlist_append(drvpar->parlist, p);

    /* Instrumental parameters */

    p = cpl_parameter_new_value("slitw", CPL_TYPE_DOUBLE, "", "", 0.4);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("slitw_key", CPL_TYPE_STRING, "", "",
                                "ESO INS SLIT1 WID");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("pixsc", CPL_TYPE_DOUBLE, "", "", 0.086);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("pixsc_key", CPL_TYPE_STRING, "", "", "NONE");
    cpl_parameterlist_append(drvpar->parlist, p);

    /* LBLRTM parameters */

    p = cpl_parameter_new_value("spec_out", CPL_TYPE_STRING, "", "", "LBL");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("singlespec", CPL_TYPE_INT, "", "", 0, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);

    /* Atmospheric profiles */

    p = cpl_parameter_new_value("ref_atm", CPL_TYPE_STRING, "", "",
                                "equ.atm");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("gdas_dir", CPL_TYPE_STRING, "", "",
                                "data/profiles/grib");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("gdas_prof", CPL_TYPE_STRING, "", "", "auto");
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_range("layers", CPL_TYPE_INT, "", "", 1, 0, 1);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("emix", CPL_TYPE_DOUBLE, "", "", 5.);
    cpl_parameterlist_append(drvpar->parlist, p);
    p = cpl_parameter_new_value("pwv", CPL_TYPE_DOUBLE, "", "", -1.);
    cpl_parameterlist_append(drvpar->parlist, p);

    /* Create parameter tables */

    /* Table for molecules */

    drvpar->molectab = cpl_table_new(1);
    cpl_table_new_column(drvpar->molectab, "list_molec", CPL_TYPE_STRING);
    cpl_table_fill_column_window_string(drvpar->molectab, "list_molec", 0, 1,
                                        mol);
    cpl_table_new_column(drvpar->molectab, "fit_molec", CPL_TYPE_INT);
    cpl_table_fill_column_window(drvpar->molectab, "fit_molec", 0, 1, 1);
    cpl_table_new_column(drvpar->molectab, "relcol", CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window(drvpar->molectab, "relcol", 0, 1, 1.);

    /* Table for range-related parameters */

    drvpar->rangetab = cpl_table_new(1);
    cpl_table_new_column(drvpar->rangetab, "chip", CPL_TYPE_INT);
    cpl_table_new_column(drvpar->rangetab, "fit_range", CPL_TYPE_INT);
    cpl_table_new_column_array(drvpar->rangetab, "cont_coef", CPL_TYPE_DOUBLE,
                               1);
    cpl_table_new_column(drvpar->rangetab, "pixres", CPL_TYPE_DOUBLE);
    cpl_table_new_column(drvpar->rangetab, "wn_start", CPL_TYPE_DOUBLE);
    cpl_table_new_column(drvpar->rangetab, "wn_end", CPL_TYPE_DOUBLE);
    cpl_table_new_column(drvpar->rangetab, "wn_step", CPL_TYPE_DOUBLE);
    cpl_table_new_column(drvpar->rangetab, "lnfl", CPL_TYPE_STRING);

    /* Table for chip-related parameters */

    drvpar->chiptab = cpl_table_new(1);
    cpl_table_new_column(drvpar->chiptab, "fit_chip", CPL_TYPE_INT);
    cpl_table_new_column_array(drvpar->chiptab, "wlc_coef", CPL_TYPE_DOUBLE,
                               2);
    cpl_table_new_column(drvpar->chiptab, "wl_min", CPL_TYPE_DOUBLE);
    cpl_table_new_column(drvpar->chiptab, "wl_max", CPL_TYPE_DOUBLE);

    /* Matrix for fixed kernel */

    drvpar->kernel = cpl_matrix_new(1, 1);
    cpl_matrix_set(drvpar->kernel, 0, 0, 1.);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_par_readfile(mfdrv *drvpar, const char *parfile, char isCalctrans)
{
    /*!
     * \callgraph
     *
     * Reads the MOLECFIT driver file and puts the parameters in a ::mfdrv
     * structure.
     *
     * \b INPUT:
     * \param drvpar   ::mfdrv parameter structure
     * \param parfile  name of parameter file
     *
     * \b OUTPUT:
     * \param drvpar   ::mfdrv parameter structure with read values
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    cpl_array *cont_coef = NULL, *wlc_coef = NULL;
    mfpar par[MF_MAXPAR];
    cpl_parameter *p, *q;
    char errtxt[MF_MAXLEN];
    char str[MF_MAXLEN] = "", tag[MF_MAXTAGLEN+1] = "";
    int npar = MF_MAXPAR, nmolec = 0, i = 0, j = 0, nrange = 0, nparc = 0;
    int n = 0, writerespar = 0, nchip = 0, nparw = 0;

    /* Fill parameter list */

    /* Directory parameters */

    // MF_PAR_SEARCH_CHECK(par, &npar, "basedir", parfile);
    // if (par[0].c[0] != '#' && npar == 2) {
    //     /* Relative base path */
    //     p = cpl_parameterlist_find(drvpar->parlist, "relbasedir");
    //     cpl_parameter_set_string(p, par[1].c);
    //     /* Absolute base path */
    //     if (getcwd(cwd, sizeof(cwd)) == NULL) {
    //         return cpl_error_set_message(cpl_func, MF_ERROR_GETCWD,
    //                                      MF_ERROR_GETCWD_TXT);
    //     }
    //     mf_basic_abspath(basedir, par[1].c, cwd);
    //     p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    //     cpl_parameter_set_string(p, basedir);
    // }


    MF_PAR_SEARCH_S("bindir", parfile, drvpar->parlist);
    MF_PAR_SEARCH_S("configdir", parfile, drvpar->parlist);
    MF_PAR_SEARCH_S("datadir", parfile, drvpar->parlist);
    /* only used by GUI, needed to write it out again in expert mode */
    MF_PAR_SEARCH_S("user_workdir", parfile, drvpar->parlist);

    /* Input file parameters */

    MF_PAR_SEARCH_S("filename", parfile, drvpar->parlist);
    MF_PAR_SEARCH_S("listname", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("trans", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "columns", parfile);
    if (par[0].c[0] != '#' && npar == 5) {
        p = cpl_parameterlist_find(drvpar->parlist, "col_lam");
        cpl_parameter_set_string(p, par[1].c);
        p = cpl_parameterlist_find(drvpar->parlist, "col_flux");
        cpl_parameter_set_string(p, par[2].c);
        p = cpl_parameterlist_find(drvpar->parlist, "col_dflux");
        cpl_parameter_set_string(p, par[3].c);
        p = cpl_parameterlist_find(drvpar->parlist, "col_mask");
        cpl_parameter_set_string(p, par[4].c);
    }
    MF_PAR_SEARCH_D("default_error", parfile, drvpar->parlist);
    MF_PAR_SEARCH_D("wlgtomicron", parfile, drvpar->parlist);
    MF_PAR_SEARCH_S("vac_air", parfile, drvpar->parlist);
    if (!isCalctrans) {
        MF_PAR_SEARCH_S("wrange_include", parfile, drvpar->parlist);
        MF_PAR_SEARCH_S("wrange_exclude", parfile, drvpar->parlist);
        MF_PAR_SEARCH_S("prange_exclude", parfile, drvpar->parlist);
    }

    /* Output file parameters */

    MF_PAR_SEARCH_S("output_dir", parfile, drvpar->parlist);
    MF_PAR_SEARCH_S("output_name", parfile, drvpar->parlist);
    MF_PAR_SEARCH_S("plot_creation", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("plot_range", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("writeascii", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("clean_mflux", parfile, drvpar->parlist);

    /* Parameters for fitting (molecules, continuum, wavelength, resol.) */

    MF_PAR_SEARCH_D("ftol", parfile, drvpar->parlist);
    MF_PAR_SEARCH_D("xtol", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("flux_unit", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("fit_back", parfile, drvpar->parlist);
    MF_PAR_SEARCH_D("telback", parfile, drvpar->parlist);
    // MF_PAR_SEARCH_I("nrange_max", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("fit_cont", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("cont_n", parfile, drvpar->parlist);
    MF_PAR_SEARCH_D("cont_const", parfile, drvpar->parlist);
    // MF_PAR_SEARCH_I("nchip_max", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("fit_wlc", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "wlc_n", parfile);
    if (par[0].c[0] != '#' && npar == 2) {
        p = cpl_parameterlist_find(drvpar->parlist, "wlc_n");
        q = cpl_parameterlist_find(drvpar->parlist, "fit_wlc_lin");
        if (par[1].i == 0) {
            cpl_parameter_set_int(p, 1);
            cpl_parameter_set_int(q, 0);
        } else {
            cpl_parameter_set_int(p, par[1].i);
            cpl_parameter_set_int(q, 1);
        }
    }
    MF_PAR_SEARCH_D("wlc_const", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("fit_res_box", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "relres_box", parfile);
    if (par[0].c[0] != '#' && npar == 2) {
        if (par[1].d < 0. || par[1].d > 2.) {
            sprintf(errtxt, "%s: Value of relres_box in parameter file (%s)"
                    " out of range [0,2]\n", MF_ERROR_IIP_TXT, parfile);
            return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s",
                                         errtxt);
        }
        p = cpl_parameterlist_find(drvpar->parlist, "relres_box");
        cpl_parameter_set_double(p, par[1].d);
    }
    MF_PAR_SEARCH_I("kernmode", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("fit_res_gauss", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "res_gauss", parfile);
    if (par[0].c[0] != '#' && npar == 2) {
        if (par[1].d < 0. || par[1].d > 100.) {
            sprintf(errtxt, "%s: Value of res_gauss in parameter file (%s)"
                    " out of range [0,100]\n", MF_ERROR_IIP_TXT, parfile);
            return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s",
                                         errtxt);
        }
        p = cpl_parameterlist_find(drvpar->parlist, "res_gauss");
        cpl_parameter_set_double(p, par[1].d);
    }
    MF_PAR_SEARCH_I("fit_res_lorentz", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "res_lorentz", parfile);
    if (par[0].c[0] != '#' && npar == 2) {
        if (par[1].d < 0. || par[1].d > 100.) {
            sprintf(errtxt, "%s: Value of res_lorentz in parameter file (%s)"
                    " out of range [0,100]\n", MF_ERROR_IIP_TXT, parfile);
            return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s",
                                         errtxt);
        }
        p = cpl_parameterlist_find(drvpar->parlist, "res_lorentz");
        cpl_parameter_set_double(p, par[1].d);
    }
    MF_PAR_SEARCH_D("kernfac", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("varkern", parfile, drvpar->parlist);
    MF_PAR_SEARCH_S("kernel_file", parfile, drvpar->parlist);

    /* Parameters for ambient conditions */

    MF_PAR_SEARCH_D("obsdate", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "obsdate_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "obsdate_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }
    MF_PAR_SEARCH_D("utc", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "utc_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "utc_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }
    MF_PAR_SEARCH_D("telalt", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "telalt_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "telalt_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }
    MF_PAR_SEARCH_D("rhum", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "rhum_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "rhum_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }
    MF_PAR_SEARCH_D("pres", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "pres_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "pres_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }
    MF_PAR_SEARCH_D("temp", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "temp_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "temp_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }
    MF_PAR_SEARCH_D("m1temp", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "m1temp_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "m1temp_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }
    MF_PAR_SEARCH_D("geoelev", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "geoelev_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "geoelev_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }
    MF_PAR_SEARCH_D("longitude", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "longitude_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "longitude_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }
    MF_PAR_SEARCH_D("latitude", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "latitude_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "latitude_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }

    /* Instrumental parameters */

    MF_PAR_SEARCH_D("slitw", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "slitw_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "slitw_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }
    MF_PAR_SEARCH_D("pixsc", parfile, drvpar->parlist);
    MF_PAR_SEARCH_CHECK(par, &npar, "pixsc_key", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        p = cpl_parameterlist_find(drvpar->parlist, "pixsc_key");
        mf_basic_joinpar(str, par, 1, npar-1);
        cpl_parameter_set_string(p, str);
    }

    /* Atmospheric profiles */

    MF_PAR_SEARCH_S("ref_atm", parfile, drvpar->parlist);
    MF_PAR_SEARCH_S("gdas_prof", parfile, drvpar->parlist);
    MF_PAR_SEARCH_S("gdas_dir", parfile, drvpar->parlist);
    MF_PAR_SEARCH_S("gdas_prof", parfile, drvpar->parlist);
    MF_PAR_SEARCH_I("layers", parfile, drvpar->parlist);
    MF_PAR_SEARCH_D("emix", parfile, drvpar->parlist);
    MF_PAR_SEARCH_D("pwv", parfile, drvpar->parlist);

    /* Fill table for molecules */

    MF_PAR_SEARCH_CHECK(par, &npar, "list_molec", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        nmolec = npar - 1;
        if (nmolec == 0) {
            sprintf(errtxt, "%s: par. list_molec of %s (no molecules given)",
                    MF_ERROR_UFS_TXT, parfile);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }
        cpl_table_set_size(drvpar->molectab, nmolec);
        for (i = 0; i < nmolec; i++) {
          for (j = 0; j < (int) strlen(par[i+1].c); j++) {
                par[i+1].c[j] = toupper(par[i+1].c[j]);
            }
            cpl_table_set_string(drvpar->molectab, "list_molec", i,
                                 par[i+1].c);
        }
    }
    MF_PAR_SEARCH_CHECK(par, &npar, "fit_molec", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        for (i = 0; i < nmolec; i++) {
            if (i >= npar-1) {
                if (i == npar-1) {
                    cpl_msg_warning(cpl_func, "Par. fit_molec of %s: "
                                    "N < nmolec "
                                    "-> Use last value for missing flags",
                                    parfile);
                }
                cpl_table_set_int(drvpar->molectab, "fit_molec", i,
                                  par[npar-1].i);
            } else {
                cpl_table_set_int(drvpar->molectab, "fit_molec", i,
                                  par[i+1].i);
            }
        }
    }
    MF_PAR_SEARCH_CHECK(par, &npar, "relcol", parfile);
    if (par[0].c[0] != '#' && npar > 1) {
        for (i = 0; i < nmolec; i++) {
            if (i >= npar-1) {
                if (i == npar-1) {
                    cpl_msg_warning(cpl_func, "Par. relcol of %s: N < nmolec "
                                    "-> Use last value for missing factors",
                                    parfile);
                }
                cpl_table_set_double(drvpar->molectab, "relcol", i,
                                     par[npar-1].d);
            } else {
                if (par[i+1].d < 1e-5 || par[i+1].d > 100.) {
                    sprintf(errtxt, "%s: Value #%d of relcol in parameter "
                            "file (%s) out of range [1e-5,1e2]\n",
                            MF_ERROR_IIP_TXT, i+1, parfile);
                    return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s",
                                            errtxt);
                }
                cpl_table_set_double(drvpar->molectab, "relcol", i,
                                     par[i+1].d);
            }
        }
    }
    if (mf_par_finalize_lbl_molecs(drvpar) != CPL_ERROR_NONE) {
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                     errtxt);
    }

    /* Fill table for range-related parameters */

    // Get the number of fit_rangei parameters
    if (mf_par_regex_number(&nrange, "fit_range", parfile)) {
        sprintf(errtxt, "%s: error reading fit_range parameters in %s.",
                MF_ERROR_UFS_TXT, parfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                     errtxt);
    }
    // p = cpl_parameterlist_find(drvpar->parlist, "nrange_max");
    // nrange = cpl_parameter_get_int(p);

    for (nparc = 0, n = 0; n < nrange; n++) {
        sprintf(tag, "fit_range%d", n+1);
        MF_PAR_SEARCH_CHECK(par, &npar, tag, parfile);
        if (par[0].c[0] != '#' && npar == 2) {
            cpl_table_set_size(drvpar->rangetab, n+1);
            cpl_table_set(drvpar->rangetab, "fit_range", n, par[1].i);
        }
        if (par[0].c[0] != '#') {
            writerespar = 1;
        }
        sprintf(tag, "cont_range%d", n+1);
        MF_PAR_SEARCH_CHECK(par, &npar, tag, parfile);
        if (par[0].c[0] != '#' && npar > 1) {
            cpl_table_set_size(drvpar->rangetab, n+1);
            if (npar-1 > nparc) {
                nparc = npar - 1;
                cpl_table_set_column_depth(drvpar->rangetab, "cont_coef",
                                           nparc);
            }
            cont_coef = cpl_array_new(npar-1, CPL_TYPE_DOUBLE);
            for (i = 0; i < npar-1; i++) {
                cpl_array_set(cont_coef, i, par[i+1].d);
            }
            cpl_table_set_array(drvpar->rangetab, "cont_coef", n, cont_coef);
            cpl_array_delete(cont_coef);
        }
        if (par[0].c[0] != '#') {
            writerespar = 1;
        }
    }

    /* Fill table for chip-related parameters */

    // Get the number of fit_chipi parameters
    if (mf_par_regex_number(&nchip, "fit_chip", parfile)) {
        sprintf(errtxt, "%s: error reading fit_chip parameters in %s.",
                MF_ERROR_UFS_TXT, parfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                     errtxt);
    }
    // p = cpl_parameterlist_find(drvpar->parlist, "nchip_max");
    // nchip = cpl_parameter_get_int(p);

    for (nparw = 0, n = 0; n < nchip; n++) {
        sprintf(tag, "fit_chip%d", n+1);
        MF_PAR_SEARCH_CHECK(par, &npar, tag, parfile);
        if (par[0].c[0] != '#' && npar == 2) {
            cpl_table_set_size(drvpar->chiptab, n+1);
            cpl_table_set(drvpar->chiptab, "fit_chip", n, par[1].i);
        }
        if (par[0].c[0] != '#') {
            writerespar = 1;
        }
        sprintf(tag, "wlc_chip%d", n+1);
        MF_PAR_SEARCH_CHECK(par, &npar, tag, parfile);
        if (par[0].c[0] != '#' && npar > 1) {
            cpl_table_set_size(drvpar->chiptab, n+1);
            if (npar-1 > nparw) {
                nparw = npar - 1;
                cpl_table_set_column_depth(drvpar->chiptab, "wlc_coef",
                                           nparw);
            }
            wlc_coef = cpl_array_new(npar-1, CPL_TYPE_DOUBLE);
            for (i = 0; i < npar-1; i++) {
                cpl_array_set(wlc_coef, i, par[i+1].d);
            }
            cpl_table_set_array(drvpar->chiptab, "wlc_coef", n, wlc_coef);
            cpl_array_delete(wlc_coef);
        }
        if (par[0].c[0] != '#') {
            writerespar = 1;
        }
    }

    parse_end:

    /* Exit in case of errors */
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        return cpl_error_get_code();
    }

    /* Write parameter file with fit results? */
    p = cpl_parameterlist_find(drvpar->parlist, "writerespar");
    cpl_parameter_set_int(p, writerespar);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_par_finalize_lbl_molecs(mfdrv * drvpar)
{
    /*!
     * Updates drvpar based on contents of filled molectab.
     * Fills nmolec and lbl_molecs
     *
     * \b INPUT:
     * \param drvpar   ::mfdrv parameter structure
     *
     */
    cpl_size nmolec = cpl_table_get_nrow(drvpar->molectab);
    cpl_parameter * p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    cpl_parameter_set_int(p, nmolec);
    for (cpl_size i = 0; i < nmolec; i++) {
        const char * molec = cpl_table_get_string(drvpar->molectab, "list_molec", i);
        p = cpl_parameterlist_find(drvpar->parlist, "lbl_molecs");
        mf_par_setmolec_all(p, molec);
    }

    return cpl_error_get_code();
}


cpl_error_code mf_par_copyfile(const mfdrv *drvpar, const char *parfile)
{
    /*!
     * Copies MOLECFIT driver file to output directory and renames it.
     *
     * \b INPUT:
     * \param drvpar   ::mfdrv parameter structure
     * \param parfile  name of input MOLECFIT driver file
     *
     * \b ERRORS:
     * - File opening failed
     * - see ::mf_basic_access
     */

    FILE *stream;
    cpl_error_code errcode = CPL_ERROR_NONE;
    cpl_parameter *p;
    char errtxt[MF_MAXLEN];
    char outdir[MF_MAXLEN] = "", outname[MF_MAXLEN] = "";
    char path[MF_MAXLEN] = "", outfile[MF_MAXLEN] = "";
    char sys[MF_MAXLEN] = "";
    int d = 0;

    /* Check existence of input parameter file */
    if ((stream = fopen(parfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, parfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }
    fclose(stream);

    /* Get output path and name */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    strcpy(outdir, cpl_parameter_get_string(p));
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strcpy(outname, cpl_parameter_get_string(p));
    mf_basic_abspath(path, outdir, curdir);
    sprintf(outfile, "%s%s_fit.par", path, outname);

    /* Check existence of output directory */
    if ((errcode = mf_basic_access(path, F_OK)) != CPL_ERROR_NONE) {
        return errcode;
    }

    /* Copy parameter file */
    sprintf(sys, "cp %s %s", parfile, outfile);
    if ((d = system(sys))){};

    return CPL_ERROR_NONE;
}


cpl_error_code mf_par_writefile(const mfdrv *drvpar, const char *parfile)
{
    /*!
     * Writes the ::mfdrv structure that contains the MOLECFIT parameters
     * in a file whose name is a structure parameter itself (bestfit_file).
     * Comments are taken from the input parameter file.
     * Only the values of parameters set in the input MOLECFIT driver file are
     * written into the output file.
     *
     * \b INPUT:
     * \param drvpar   ::mfdrv parameter structure
     * \param parfile  name of input MOLECFIT driver file
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     * - Invalid object structure
     */

    FILE *stream;
    cpl_parameter *p, *q, *r, *s;
    cpl_array *cont_coef = NULL, *wlc_coef = NULL;
    cpl_type type;
    char errtxt[MF_MAXLEN], line[MF_MAXNLINE][MF_MAXLEN], *tagstr = NULL;
    char tag[MF_MAXTAGLEN+1] = "";
    char outdir[MF_MAXLEN] = "", outname[MF_MAXLEN] = "";
    char path[MF_MAXLEN] = "", outfile[MF_MAXLEN] = "", **cp = NULL;
    char str[MF_MAXLEN] = "", strcomp[MF_LENLINE+2] = "";
    char tag0[MF_MAXTAGLEN+1] = "";
    size_t taglen;
    int *ip = NULL;
    int nmolec = 0, nrangemax = 0, nchipmax = 0, i = 0, nline = 0, j = 0;
    int n = 0, nrange = 0, ncont = 0, nchip = 0, nwlc = 0;
    double *dp = NULL;

    /* Check flag for file creation and return if writing is not requested */
    p = cpl_parameterlist_find(drvpar->parlist, "writerespar");
    if (cpl_parameter_get_int(p) != 1) {
        return CPL_ERROR_NONE;
    }

    /* Read MOLECFIT driver file */

    /* Check existence of input parameter file */
    if ((stream = fopen(parfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, parfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Check structure of input parameter file */

    for (i = 0; i < MF_MAXNLINE; i++) {
        if (fgets(line[i], MF_MAXLEN, stream) == NULL) {
            fclose(stream);
            sprintf(errtxt, "%s: %s (reading of line %d failed)",
                    MF_ERROR_UFS_TXT, parfile, i+1);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }
        nline = i + 1;
        if (strncmp(line[i], "end", 3) == 0) {
            fclose(stream);
            break;
        }
    }

    /* Get number of molecules */
    p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(p);

    /* Get maximum number of ranges */
    // p = cpl_parameterlist_find(drvpar->parlist, "nrange_max");
    // nrangemax = cpl_parameter_get_int(p);
    // Get the number of fit_rangei parameters
    if (mf_par_regex_number(&nrangemax, "fit_range", parfile)) {
        sprintf(errtxt, "%s: error reading fit_range parameters in %s.",
                MF_ERROR_UFS_TXT, parfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                     errtxt);
    }

    /* Get maximum number of chips */
    // p = cpl_parameterlist_find(drvpar->parlist, "nchip_max");
    // nchipmax = cpl_parameter_get_int(p);
    // Get the number of fit_chipi parameters
    if (mf_par_regex_number(&nchipmax, "fit_chip", parfile)) {
        sprintf(errtxt, "%s: error reading fit_chip parameters in %s.",
                MF_ERROR_UFS_TXT, parfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                     errtxt);
    }

    /* Get output path and name */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    strcpy(outdir, cpl_parameter_get_string(p));
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strcpy(outname, cpl_parameter_get_string(p));
    mf_basic_abspath(path, outdir, curdir);
    sprintf(outfile, "%s%s_fit.rpar", path, outname);

    /* Write info message */
    cpl_msg_info(cpl_func, "Write driver file %s with best-fit parameters",
                 outfile);

    /* Write output parameter file */

    /* Open output parameter file */
    if ((stream = fopen(outfile, "w+")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, outfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    for (i = 0; i < nline; i++) {

        if (line[i][0] != '#' && strchr(line[i], ':') != NULL) {

            /* Get parameter tag */
            tagstr = strtok(line[i], ":");
            mf_basic_strtrim_inplace(tagstr);
            sprintf(tag, "%s", tagstr);
            taglen = strlen(tag);

            /* Write parameter list to file */

            /* Parameters with special handling */

            // if (strncmp(tag, "basedir", taglen) == 0) {
            //     p = cpl_parameterlist_find(drvpar->parlist, "relbasedir");
            //     fprintf(stream, "%s: %s\n", tag, cpl_parameter_get_string(p));
            // } else 
            if (strncmp(tag, "columns", taglen) == 0) {
                p = cpl_parameterlist_find(drvpar->parlist, "col_lam");
                q = cpl_parameterlist_find(drvpar->parlist, "col_flux");
                r = cpl_parameterlist_find(drvpar->parlist, "col_dflux");
                s = cpl_parameterlist_find(drvpar->parlist, "col_mask");
                fprintf(stream, "%s: %s %s %s %s\n", tag,
                        cpl_parameter_get_string(p),
                        cpl_parameter_get_string(q),
                        cpl_parameter_get_string(r),
                        cpl_parameter_get_string(s));
            } else if (strncmp(tag, "wlc_n", taglen) == 0) {
                p = cpl_parameterlist_find(drvpar->parlist, tag);
                q = cpl_parameterlist_find(drvpar->parlist, "fit_wlc_lin");
                if (cpl_parameter_get_int(q) == 0) {
                    fprintf(stream, "%s: 0\n", tag);
                } else {
                    fprintf(stream, "%s: %d\n", tag,
                            cpl_parameter_get_int(p));
                }
            }

            /* Parameters from table of molecules */

            else if (strncmp(tag, "list_molec", taglen) == 0) {
                cp = cpl_table_get_data_string(drvpar->molectab, tag);
                sprintf(str, "%s: %s", tag, cp[0]);
                if (nmolec > 1) {
                    for (j = 1; j < nmolec; j++) {
                        sprintf(strcomp, " %s", cp[j]);
                        strncat(str, strcomp, strlen(strcomp));
                    }
                }
                strncat(str, "\n", 1);
                fprintf(stream, "%s", str);
            } else if (strncmp(tag, "fit_molec", taglen) == 0) {
                ip = cpl_table_get_data_int(drvpar->molectab, tag);
                sprintf(str, "%s: %d", tag, ip[0]);
                if (nmolec > 1) {
                    for (j = 1; j < nmolec; j++) {
                        sprintf(strcomp, " %d", ip[j]);
                        strncat(str, strcomp, strlen(strcomp));
                    }
                }
                strncat(str, "\n", 1);
                fprintf(stream, "%s", str);
            } else if (strncmp(tag, "relcol", taglen) == 0) {
                dp = cpl_table_get_data_double(drvpar->molectab, tag);
                sprintf(str, "%s: %g", tag, dp[0]);
                if (nmolec > 1) {
                    for (j = 1; j < nmolec; j++) {
                        sprintf(strcomp, " %g", dp[j]);
                        strncat(str, strcomp, strlen(strcomp));
                    }
                }
                strncat(str, "\n", 1);
                fprintf(stream, "%s", str);
            }

            /* Parameters from range-related table */

            else if (strncmp(tag, "fit_range", 9) == 0) {
                nrange = cpl_table_get_nrow(drvpar->rangetab);
                for (n = 0; n < nrangemax; n++) {
                    sprintf(tag0, "fit_range%d", n+1);
                    if (strcmp(tag, tag0) == 0 && n < nrange) {
                        sprintf(str, "%s: %d\n", tag,
                                cpl_table_get_int(drvpar->rangetab,
                                                  "fit_range", n, NULL));
                        fprintf(stream, "%s", str);
                    }
                }
            } else if (strncmp(tag, "cont_range", 10) == 0) {
                nrange = cpl_table_get_nrow(drvpar->rangetab);
                for (n = 0; n < nrangemax; n++) {
                    sprintf(tag0, "cont_range%d", n+1);
                    if (strcmp(tag, tag0) == 0 && n < nrange) {
                        cont_coef = cpl_array_duplicate(cpl_table_get_array(
                                    drvpar->rangetab, "cont_coef", n));
                        dp = cpl_array_get_data_double(cont_coef);
                        sprintf(str, "%s: %g", tag, dp[0]);
                        ncont = cpl_array_get_size(cont_coef);
                        if (ncont > 1) {
                            for (j = 1; j < ncont; j++) {
                                sprintf(strcomp, " %g", dp[j]);
                                strncat(str, strcomp, strlen(strcomp));
                            }
                        }
                        strncat(str, "\n", 1);
                        fprintf(stream, "%s", str);
                        cpl_array_delete(cont_coef);
                    }
                }
            }

            /* Parameters from chip-related table */

            else if (strncmp(tag, "fit_chip", 8) == 0) {
                nchip = cpl_table_get_nrow(drvpar->chiptab);
                for (n = 0; n < nchipmax; n++) {
                    sprintf(tag0, "fit_chip%d", n+1);
                    if (strcmp(tag, tag0) == 0 && n < nchip) {
                        sprintf(str, "%s: %d\n", tag,
                                cpl_table_get_int(drvpar->chiptab,
                                                  "fit_chip", n, NULL));
                        fprintf(stream, "%s", str);
                    }
                }
            } else if (strncmp(tag, "wlc_chip", 8) == 0) {
                nchip = cpl_table_get_nrow(drvpar->chiptab);
                for (n = 0; n < nchipmax; n++) {
                    sprintf(tag0, "wlc_chip%d", n+1);
                    if (strcmp(tag, tag0) == 0 && n < nchip) {
                        wlc_coef = cpl_array_duplicate(cpl_table_get_array(
                                   drvpar->chiptab, "wlc_coef", n));
                        dp = cpl_array_get_data_double(wlc_coef);
                        sprintf(str, "%s: %g", tag, dp[0]);
                        nwlc = cpl_array_get_size(wlc_coef);
                        if (nwlc > 1) {
                            for (j = 1; j < nwlc; j++) {
                                sprintf(strcomp, " %g", dp[j]);
                                strncat(str, strcomp, strlen(strcomp));
                            }
                        }
                        strncat(str, "\n", 1);
                        fprintf(stream, "%s", str);
                        cpl_array_delete(wlc_coef);
                    }
                }
            }

            /* Other parameters */

            else {
                if ((p = cpl_parameterlist_find(drvpar->parlist, tag)) !=
                    NULL) {
                    type = cpl_parameter_get_type(p);
                    if (type == CPL_TYPE_STRING) {
                        fprintf(stream, "%s: %s\n", tag,
                                cpl_parameter_get_string(p));
                    } else if (type == CPL_TYPE_INT) {
                        fprintf(stream, "%s: %d\n", tag,
                                cpl_parameter_get_int(p));
                    } else if (type == CPL_TYPE_DOUBLE) {
                        fprintf(stream, "%s: %g\n", tag,
                                cpl_parameter_get_double(p));
                    } else {
                        fclose(stream);
                        sprintf(errtxt, "%s: mfdrv *drvpar (unexpected "
                                "type of parameter %s)", MF_ERROR_IOS_TXT,
                                tag);
                        return cpl_error_set_message(cpl_func, MF_ERROR_IOS,
                                                     "%s", errtxt);
                    }
                } else {
                    /* Error if no match */
                    fclose(stream);
                    sprintf(errtxt, "%s: %s (unexpected parameter %s)",
                            MF_ERROR_UFS_TXT, outfile, line[i]);
                    return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                                 errtxt);
                }
            }

        } else {

            /* Write comment, empty line, inactive parameter tag, or "end" */
            fprintf(stream, "%s", line[i]);

        }

    }

    /* Close output file */
    fclose(stream);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_par_deleteall(mfdrv *drvpar)
{
    /*!
     * Frees the memory occupied by a ::mfdrv parameter structure.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param drvpar  empty ::mfdrv parameter structure
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameterlist_delete(drvpar->parlist);
    cpl_table_delete(drvpar->rangetab);
    cpl_table_delete(drvpar->chiptab);
    cpl_table_delete(drvpar->molectab);
    cpl_matrix_delete(drvpar->kernel);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_par_search(mfpar par[], int *npar, const char *parname,
                             const char *parfile)
{
    /*!
     * \callgraph
     *
     * Searches for a given parameter in a parameter file.
     * An active parameter has to be marked by the name + a colon without any
     * spaces in between. Then the rest of the line contains the parameter
     * values. They are returned by means of a ::mfpar structure, which
     * provides the parameter name at position 0. Moreover, the read number of
     * parameters is provided.
     * If the desired parameter is not present in the file, the first and
     * only character in the ::mfpar structure will become a # and the routine
     * will return MF_ERROR_UFS.
     *
     * \b INPUT:
     * \param par      array of ::mfpar structures
     *                 (to transfer parameter values)
     * \param parname  parameter ID
     * \param parfile  name of parameter file
     *
     * \b OUTPUT:
     * \param par      ::mfpar array that contains the read value(s) as
     *                 character string ("c"), integer ("i"), or double
     *                 precision floating point number ("d").
     *                 The different data types can be accessed by adding the
     *                 corresponding suffixes c, i, or d to the name of the
     *                 ::mfpar structure.
     * \param npar     found number of parameter values (at one line)
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    char errtxt[MF_MAXLEN], tag[MF_MAXTAGLEN+2];
    size_t taglen;

    /* Set npar to zero */
    *npar = 0;

    /* Default parameter values */
    strcpy(par[0].c, "");
    par[0].i = 0;
    par[0].d = 0.;

    /* Check existence of parameter file */
    if ((stream = fopen(parfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, parfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Search for parameter by reading each line */

    sprintf(tag, "%s:", parname);
    taglen = strlen(tag);

    do {
        cpl_error_code err = mf_basic_readline(stream, par, npar);
        if (err != CPL_ERROR_NONE) {
            fclose(stream);
            strcpy(par[0].c, "#");
            return err;
        }
        if (strncmp(par[0].c, "end", 3) == 0 ||
            strncmp(par[0].c, "eof", 3) == 0) {
            fclose(stream);
            strcpy(par[0].c, "#");
            return MF_ERROR_UFS;
        }
    } while (par[0].c[0] == '#' || par[0].c[0] == '\n' ||
             strncmp(par[0].c, tag, taglen) != 0);

    fclose(stream);

    return CPL_ERROR_NONE;
}

cpl_error_code mf_par_regex_number(int *npar, const char *parname,
                             const char *parfile)
{

    FILE *stream;
    char errtxt[MF_MAXLEN];

    /* Set npar to zero */
    *npar = 0;

    /* Check existence of parameter file */
    if ((stream = fopen(parfile, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, parfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Search for parameter by reading each line */
    regex_t regex;
    char regexstr[MF_MAXLEN];
    sprintf(regexstr, "^%s([1-9]*[0-9]):.*", parname);
    int reti = regcomp(&regex, regexstr, REG_EXTENDED);
    if (reti) {
        return MF_ERROR_UFS;
    }

    while (1) {
        char* buffer = NULL;
        size_t len = 0;
        if (getline(&buffer, &len, stream) == -1) {
            free(buffer);
            break;
        }
        // Regex
        regmatch_t match[2];
        reti = regexec(&regex, buffer, 2, match, 0);
        free(buffer);
        if (!reti) {
            (*npar)++;
            // unsigned strlen = (unsigned)(match[1].rm_eo - match[1].rm_so) + 1;
            // char* tmp = (char*) malloc((strlen + 1) * sizeof(char));
            // strncpy(tmp, buffer + match[1].rm_so, strlen);
            // tmp[strlen - 1] = '\0';
        }

    }
    regfree(&regex);

    fclose(stream);

    return CPL_ERROR_NONE;

}

cpl_error_code mf_par_setfitflags(mfdrv *drvpar, const cpl_array *fit_molec,
                                  const cpl_array *fit_cont,
                                  const cpl_array *fit_wlc,
                                  const cpl_array *fit_res)
{
    /*!
     * Modifies the fit flags (0 or 1) in the MOLECFIT driver parameter
     * structure ::mfdrv. The input CPL arrays have to consist of nmolec
     * (fit_molec), nrange (fit_cont), nchip (fit_wlc), and 3 (fit_res)
     * elements, respectively. In the case of a radiance spectrum a fit flag
     * for the telescope background has to be added to fit_cont.
     *
     * \b INPUT:
     * \param drvpar     ::mfdrv parameter structure
     * \param fit_molec  CPL array for flags for fitting of molecules
     * \param fit_cont   CPL array for flags for continuum fit
     * \param fit_wlc    CPL array for flags for wavelength fit
     * \param fit_res    CPL array for flags for resolution fit
     *
     * \b OUTPUT:
     * \param drvpar     ::mfdrv parameter structure with modified fit flags
     *
     * \b ERRORS:
     * - Invalid object value(s)
     * - Invalid object structure
     */

    cpl_parameter *p;
    char errtxt[MF_MAXLEN];
    int nmolec = 0, nrange = 0, trans = 0, nchip = 0, i = 0, flag = 0;

    /* Get number of molecules and compare with format of parameter table */

    p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(p);
    if (nmolec != cpl_table_get_nrow(drvpar->molectab)) {
        sprintf(errtxt, "%s: nmolec of mfdrv *drvpar (number of molecules)",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Check number of fit_molec array elements */

    if (cpl_array_get_size(fit_molec) != nmolec) {
        sprintf(errtxt, "%s: cpl_array *fit_molec (size != nmolec)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get number of ranges and compare with format of parameter table */

    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    nrange = cpl_parameter_get_int(p);
    if (nrange != cpl_table_get_nrow(drvpar->rangetab)) {
        sprintf(errtxt, "%s: nrange of mfdrv *drvpar (number of ranges)",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Check number of fit_cont array elements */

    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    trans = cpl_parameter_get_int(p);
    if (trans == 0) {
        if (cpl_array_get_size(fit_cont) != nrange+1) {
            sprintf(errtxt, "%s: cpl_array *fit_cont (size != nrange+1)",
                    MF_ERROR_IOS_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }
    } else {
        if (cpl_array_get_size(fit_cont) != nrange) {
            sprintf(errtxt, "%s: cpl_array *fit_cont (size != nrange)",
                    MF_ERROR_IOS_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }
    }

    /* Get number of ranges and compare with format of parameter table */

    p = cpl_parameterlist_find(drvpar->parlist, "nchip");
    nchip = cpl_parameter_get_int(p);
    if (nchip != cpl_table_get_nrow(drvpar->chiptab)) {
        sprintf(errtxt, "%s: nchip of mfdrv *drvpar (number of chips)",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Check number of fit_wlc array elements */

    if (cpl_array_get_size(fit_wlc) != nchip) {
        sprintf(errtxt, "%s: cpl_array *fit_wlc (size != nchip)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check number of fit_res array elements */

    if (cpl_array_get_size(fit_res) != 3) {
        sprintf(errtxt, "%s: cpl_array *fit_res (size != 3)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check validity of fit_molec fit flags */

    for (i = 0; i < cpl_array_get_size(fit_molec); i++) {
        flag = cpl_array_get(fit_molec, i, NULL);
        if (flag < 0 || flag > 1) {
            sprintf(errtxt, "%s: cpl_array *fit_molec (fit flag != 0 or 1)",
                    MF_ERROR_IOV_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s",
                                         errtxt);
        }
    }

    /* Check validity of fit_cont fit flags */

    for (i = 0; i < cpl_array_get_size(fit_cont); i++) {
        flag = cpl_array_get(fit_cont, i, NULL);
        if (flag < 0 || flag > 1) {
            sprintf(errtxt, "%s: cpl_array *fit_cont (fit flag != 0 or 1)",
                    MF_ERROR_IOV_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s",
                                         errtxt);
        }
    }

    /* Check validity of fit_wlc fit flags */

    for (i = 0; i < cpl_array_get_size(fit_wlc); i++) {
        flag = cpl_array_get(fit_wlc, i, NULL);
        if (flag < 0 || flag > 1) {
            sprintf(errtxt, "%s: cpl_array *fit_wlc (fit flag != 0 or 1)",
                    MF_ERROR_IOV_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s",
                                         errtxt);
        }
    }

    /* Check validity of fit_res fit flags */

    for (i = 0; i < cpl_array_get_size(fit_res); i++) {
        flag = cpl_array_get(fit_res, i, NULL);
        if (flag < 0 || flag > 1) {
            sprintf(errtxt, "%s: cpl_array *fit_res (fit flag != 0 or 1)",
                    MF_ERROR_IOV_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s",
                                         errtxt);
        }
    }

    /* Set fit_molec fit flags */

    for (i = 0; i < nmolec; i++) {
        cpl_table_set_int(drvpar->molectab, "fit_molec", i,
                          cpl_array_get(fit_molec, i, NULL));
    }

    /* Set fit_cont fit flags */

    for (i = 0; i < nrange; i++) {
        cpl_table_set_int(drvpar->rangetab, "fit_range", i,
                          cpl_array_get(fit_cont, i, NULL));
    }
    if (trans == 0) {
        /* Fit telescope background if relevant */
        p = cpl_parameterlist_find(drvpar->parlist, "fit_back");
        cpl_parameter_set_int(p, cpl_array_get(fit_cont, nrange, NULL));
    }

    /* Set fit_wlc fit flag */

    for (i = 0; i < nchip; i++) {
        cpl_table_set_int(drvpar->chiptab, "fit_chip", i,
                          cpl_array_get(fit_wlc, i, NULL));
    }

    /* Set fit_res fit flags */

    p = cpl_parameterlist_find(drvpar->parlist, "fit_res_box");
    cpl_parameter_set_int(p, cpl_array_get(fit_res, 0, NULL));
    p = cpl_parameterlist_find(drvpar->parlist, "fit_res_gauss");
    cpl_parameter_set_int(p, cpl_array_get(fit_res, 1, NULL));
    p = cpl_parameterlist_find(drvpar->parlist, "fit_res_lorentz");
    cpl_parameter_set_int(p, cpl_array_get(fit_res, 2, NULL));

    return CPL_ERROR_NONE;
}


void mf_par_setmolec_all(cpl_parameter *p, const char *mol) {
    /*!
     * \brief
     *   Change state of all molecules
     *
     * This function changes the "off"-state of all molecules to an
     * "on"-state, i.e. enables the molecule for LBLRTM/LNFL.
     *
     * \b INPUT:
     * \param p    ::mfdrv parameter lbl_molecs
     * \param mol  current molecule
     *
     * \b OUTPUT:
     * \param p    ::mfdrv parameter lbl_molecs
     */

    /*  New list of Molecules
     
    (M):  AVAILABLE  ( 1)  H2O  ( 2)  CO2  ( 3)    O3 ( 4)   N2O ( 5)    CO ( 6)   CH4 ( 7)    O2
          MOLECULAR  ( 8)   NO  ( 9)  SO2  (10)   NO2 (11)   NH3 (12)  HNO3 (13)    OH (14)    HF
           SPECIES   (15)  HCL  (16)  HBR  (17)    HI (18)   CLO (19)   OCS (20)  H2CO (21)  HOCL
                     (22)   N2  (23)  HCN  (24) CH3CL (25)  H2O2 (26)  C2H2 (27)  C2H6 (28)   PH3
                     (29) COF2  (30)  SF6  (31)   H2S (32) HCOOH (33)   HO2 (34)     O (35)CLONO2
                     (36)  NO+  (37) HOBR  (38)  C2H4 (39) C3HOH (40) CH3Br (41) CH3CN (42)   CF4
                     (43) C4H2  (44) HC3N  (45)    H2 (46)    CS (47)   SO3
     */

    if (strcmp(mol, "H2O") == 0) {
        mf_par_setmolec(p, 0);
    } else if (strcmp(mol, "CO2") == 0) {
        mf_par_setmolec(p, 1);
    } else if (strcmp(mol, "O3") == 0) {
        mf_par_setmolec(p, 2);
    } else if (strcmp(mol, "N2O") == 0) {
        mf_par_setmolec(p, 3);
    } else if (strcmp(mol, "CO") == 0) {
        mf_par_setmolec(p, 4);
    } else if (strcmp(mol, "CH4") == 0) {
        mf_par_setmolec(p, 5);
    } else if (strcmp(mol, "O2") == 0) {
        mf_par_setmolec(p, 6);
    } else if (strcmp(mol, "NO") == 0) {
        mf_par_setmolec(p, 7);
    } else if (strcmp(mol, "SO2") == 0) {
        mf_par_setmolec(p, 8);
    } else if (strcmp(mol, "NO2") == 0) {
        mf_par_setmolec(p, 9);
    } else if (strcmp(mol, "NH3") == 0) {
        mf_par_setmolec(p, 10);
    } else if (strcmp(mol, "HNO3") == 0) {
        mf_par_setmolec(p, 11);
    } else if (strcmp(mol, "OH") == 0) {
        mf_par_setmolec(p, 12);
    } else if (strcmp(mol, "HF") == 0) {
        mf_par_setmolec(p, 13);
    } else if (strcmp(mol, "HCL") == 0) {
        mf_par_setmolec(p, 14);
    } else if (strcmp(mol, "HBR") == 0) {
        mf_par_setmolec(p, 15);
    } else if (strcmp(mol, "HI") == 0) {
        mf_par_setmolec(p, 16);
    } else if (strcmp(mol, "CLO") == 0) {
        mf_par_setmolec(p, 17);
    } else if (strcmp(mol, "OCS") == 0) {
        mf_par_setmolec(p, 18);
    } else if (strcmp(mol, "H2CO") == 0) {
        mf_par_setmolec(p, 19);
    } else if (strcmp(mol, "HOCL") == 0) {
        mf_par_setmolec(p, 20);
    } else if (strcmp(mol, "N2") == 0) {
        mf_par_setmolec(p, 21);
    } else if (strcmp(mol, "HCN") == 0) {
        mf_par_setmolec(p, 22);
    } else if (strcmp(mol, "CH3CL") == 0) {
        mf_par_setmolec(p, 23);
    } else if (strcmp(mol, "H2O2") == 0) {
        mf_par_setmolec(p, 24);
    } else if (strcmp(mol, "C2H2") == 0) {
        mf_par_setmolec(p, 25);
    } else if (strcmp(mol, "C2H6") == 0) {
        mf_par_setmolec(p, 26);
    } else if (strcmp(mol, "PH3") == 0) {
        mf_par_setmolec(p, 27);
    } else if (strcmp(mol, "COF2") == 0) {
        mf_par_setmolec(p, 28);
    } else if (strcmp(mol, "SF6") == 0) {
        mf_par_setmolec(p, 29);
    } else if (strcmp(mol, "H2S") == 0) {
        mf_par_setmolec(p, 30);
    } else if (strcmp(mol, "HCOOH") == 0) {
        mf_par_setmolec(p, 31);
    } else if (strcmp(mol, "HO2") == 0) {
        mf_par_setmolec(p, 32);
    } else if (strcmp(mol, "O") == 0) {
        mf_par_setmolec(p, 33);
    } else if (strcmp(mol, "CLONO2") == 0) {
        mf_par_setmolec(p, 34);
    } else if (strcmp(mol, "NO+") == 0) {
        mf_par_setmolec(p, 35);
    } else if (strcmp(mol, "HOBR") == 0) {
        mf_par_setmolec(p, 36);
    } else if (strcmp(mol, "C2H4") == 0) {
        mf_par_setmolec(p, 37);
    } else if (strcmp(mol, "CH3OH") == 0) {
        mf_par_setmolec(p, 38);
    } else if (strcmp(mol, "CH3Br") == 0) {
        mf_par_setmolec(p, 39);
    } else if (strcmp(mol, "CH3CN") == 0) {
        mf_par_setmolec(p, 40);
    } else if (strcmp(mol, "CF4") == 0) {
        mf_par_setmolec(p, 41);
    } else if (strcmp(mol, "C4H2") == 0) {
        mf_par_setmolec(p, 42);
    } else if (strcmp(mol, "HC3N") == 0) {
        mf_par_setmolec(p, 43);
    } else if (strcmp(mol, "H2") == 0) {
        mf_par_setmolec(p, 44);
    } else if (strcmp(mol, "CS") == 0) {
        mf_par_setmolec(p, 45);
    } else if (strcmp(mol, "SO3") == 0) {
        mf_par_setmolec(p, 46);
    }
}


void mf_par_setmolec(cpl_parameter *p, const int pos) {
    /*!
     * \brief
     *   Change molecule state
     *
     * This function changes the "off"-state of a molecule with number \em pos
     * to an "on"-state, i.e. enables the molecule for LBLRTM/LNFL.
     *
     * \b INPUT:
     * \param p    ::mfdrv parameter lbl_molecs
     * \param pos  number of a molecule
     *
     * \b OUTPUT:
     * \param p    ::mfdrv parameter lbl_molecs
     */

    const char *molecs_ptr = NULL;
    char molecs[48];

    molecs_ptr = cpl_parameter_get_string(p);
    strcpy(molecs, molecs_ptr);
    *(molecs+pos) = '1';
    cpl_parameter_set_string(p, molecs);
}


cpl_error_code mf_par_readkernel(mfdrv *drvpar)
{
    /*!
     * Read fixed kernel(s) from an ASCII file if provided by the user and
     * write it/them into the kernel matrix of the ::mfdrv parameter
     * structure. If no kernel is given ("none"), the default zero size matrix
     * is left untouched. In this case, the kernel the model spectrum will be
     * convolved with a synthetic kernel. A single kernel valid for the entire
     * spectrum can be provided by either all elements on a single line
     * separated by spaces or each element on a separate line. Note that the
     * maximum number of elements on a line is ::MF_MAXPAR. Is is also
     * possible to provide an individual kernel for each pixel of the input
     * spectrum. In this case, each kernel has to be on a separate line of the
     * ASCII file. The kernels have to be sorted by the pixel number. A change
     * of the chip does not restart the counting. It is required that each
     * kernel has the same number of elements. If the number of kernels is
     * lower than the number of pixels, the missing kernels are taken from the
     * last valid kernel. The central kernel element for the convolution is
     * the median element. For an even number of elements, the upper median
     * element is used. If the number of kernel elements is not constant, zero
     * values have to be added to the smaller kernels to fulfil the
     * requirements for a correct interpretation. It is not required that the
     * provided kernels are normalised. The routine will make sure that this
     * will be correct.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * -
     *
     * \b ERRORS:
     * - Invalid object value(s)
     * - File opening failed
     * - Unexpected file structure
     */

    cpl_error_code err = CPL_ERROR_NONE;
    FILE *stream;
    mfpar val[MF_MAXPAR];
    cpl_parameter *p;
    char errtxt[MF_MAXLEN] = "";
    char kernelname[MF_MAXLEN] = "", relkernelname[MF_MAXLEN] = "";
    int npix = 0, ilast = 0, nrow = 0, ncol = 0, i = 0, nval = 0, nkern = 0;
    int k = 0, j = 0, jlast = 0;
    double sum = 0.;
    double *matkern = NULL;

    /* Get name of kernel and return if "none" */
    // p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    // strncpy(basedir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "kernel_file");
    strncpy(kernelname, cpl_parameter_get_string(p), MF_MAXLEN);
    if (strcmp(kernelname, "none") == 0) {
        /* No kernel -> return */
        return CPL_ERROR_NONE;
    } else if (kernelname[0] != '/') {
        // sprintf(relkernelname, "%s/%s", basedir, kernelname);
        char curdir[MF_MAXLEN];
        p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        sprintf(relkernelname, "%s/%s", curdir, kernelname);
        mf_basic_absfile(kernelname, relkernelname);
    }

    /* Get number of pixels in input spectrum */
    p = cpl_parameterlist_find(drvpar->parlist, "npix");
    npix = cpl_parameter_get_int(p);
    if (npix == 0) {
        sprintf(errtxt, "%s: npix of mfdrv *drvpar = 0", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }
    ilast = npix - 1;

    /* Get size of kernel matrix */
    nrow = cpl_matrix_get_nrow(drvpar->kernel);
    ncol = cpl_matrix_get_ncol(drvpar->kernel);

    /* Write info message */
    cpl_msg_info(cpl_func, "Read kernel file %s", kernelname);

    /* Check file existence */
    if ((stream = fopen(kernelname, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, kernelname);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Read kernel elements from file and write them into a CPL matrix */

    for (i = 0; i < npix; i++) {

        /* Read line from file */
        err = mf_basic_readline(stream, val, &nval);
        if (err != CPL_ERROR_NONE) {
            fclose(stream);
            return err;
        }

        /* Leave loop if there are less kernels than pixels */
        if (strncmp(val[0].c, "eof", 3) == 0) {
            ilast = i - 1;
            break;
        }

        /* Set size of kernel matrix and get pointer */
        if (i == 0) {
            nkern = nval;
            cpl_matrix_resize(drvpar->kernel, 0, npix - nrow, 0,
                              nkern - ncol);
            matkern = cpl_matrix_get_data(drvpar->kernel);
        }

        /* Exit in the case of varying kernel size */
        if (nval != nkern) {
            fclose(stream);
            sprintf(errtxt, "%s: %s (varying number of elements per line)",
                    MF_ERROR_UFS_TXT, kernelname);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Are the kernel elements numbers? */
        for (k = 0; k < nkern; k++) {
            if (mf_basic_isnumber(val[k].c) == CPL_FALSE) {
                fclose(stream);
                sprintf(errtxt, "%s: %s (element %d of line %d is not a "
                        "number)", MF_ERROR_UFS_TXT, kernelname, k+1, i+1);
                return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                             errtxt);
            }
        }

        /* Resize kernel matrix if each element is given on a separate line */
        if (nkern == 1 && i > 0) {
            cpl_matrix_resize(drvpar->kernel, 0, 0, 0, 1);
            matkern = cpl_matrix_get_data(drvpar->kernel);
        }

        /* Write kernel into CPL matrix (differences for row- and column-based
           kernels) */
        if (nkern == 1) {
            matkern[i] = val[0].d;
            sum += matkern[i];
        } else {
            for (sum = 0., k = 0; k < nkern; k++) {
                j = k + nkern * i;
                matkern[j] = val[k].d;
                sum += matkern[j];
            }
        }

        /* Check validity of sum */
        if (nkern > 1 && sum <= 0.) {
            fclose(stream);
            sprintf(errtxt, "%s: %s (sum of elements <= 0)",
                    MF_ERROR_UFS_TXT, kernelname);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Normalise kernel if kernel elements are given on a line */
        if (nkern > 1) {
            for (k = 0; k < nkern; k++) {
                j = k + nkern * i;
                matkern[j] /= sum;
            }
        }

    }

    /* Close ASCII file */
    fclose(stream);

    /* Check existence of data */
    if (ilast < 0) {
        sprintf(errtxt, "%s: %s", MF_ERROR_NDA_TXT, kernelname);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check validity of sum */
    if (nkern == 1 && sum <= 0.) {
        sprintf(errtxt, "%s: %s (sum of elements <= 0)",
                MF_ERROR_UFS_TXT, kernelname);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Normalise kernel if each element is given on a separate line */
    if (nkern == 1) {
        for (i = 0; i <= ilast; i++) {
            matkern[i] /= sum;
        }
        /* Convert nkern and ilast */
        nkern = ilast + 1;
        ilast = 0;
    }

    /* Fill empty rows of kernel matrix with data of last valid kernel */
    if (ilast < npix - 1) {
        for (k = 0; k < nkern; k++) {
            jlast = k + nkern * ilast;
            for (i = ilast + 1; i < npix; i++) {
                j = k + nkern * i;
                matkern[j] = matkern[jlast];
            }
        }
    }

    cpl_msg_info(cpl_func, "mf_par_readkernel() success!");
    //cpl_matrix_dump(drvpar->kernel, stdout);

    return CPL_ERROR_NONE;
}


/**@}*/
