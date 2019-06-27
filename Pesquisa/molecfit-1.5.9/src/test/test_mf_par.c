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

/*
 * test_mf_par.c
 *
 * Authors:     Stefan Noll & ESO In-Kind Team Innsbruck
 * Created:     27 Apr 2010
 * Last update: 18 Jul 2012
 *
 * Test programme for mf_par.c
 */

#include "mf_par.h"

#include <cpl_test.h>


static void test_invalid(void)
{
    mfdrv drvpar;
    mf_par_initall(&drvpar);
    mf_par_readfile(&drvpar, "config/long_line.par", 0);
    cpl_test_error(MF_ERROR_UFS);
    mf_par_deleteall(&drvpar);

    mf_par_initall(&drvpar);
    mf_par_readfile(&drvpar, "config/invalid_simple.par", 0);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    mf_par_deleteall(&drvpar);
}

static void test_valid(void)
{
    mfdrv drvpar;

    mf_par_initall(&drvpar);
    mf_par_readfile(&drvpar, "config/valid_simple.par", 0);

    cpl_parameter *p;

    p = cpl_parameterlist_find(drvpar.parlist, "output_dir");
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(cpl_parameter_get_type(p), CPL_TYPE_STRING);
    cpl_test(strcmp(cpl_parameter_get_string(p), "/tmp/") == 0);

    p = cpl_parameterlist_find(drvpar.parlist, "xtol");
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(cpl_parameter_get_type(p), CPL_TYPE_DOUBLE);
    cpl_test_abs(cpl_parameter_get_double(p), 0.53, 1e-12);

    mf_par_deleteall(&drvpar);
}

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    /* Initialize */
    cpl_errorstate errstate = cpl_errorstate_get();

    /* Get time */
    double cs = cpl_test_get_walltime();

    cpl_array *molec_flags;
    cpl_array *cont_flags;
    cpl_array *wlc_flags;
    cpl_array *res_flags;

    char parfile[MF_LENLINE+2] = "config/molecfit_test_crires.par";
    int nmolec = 3;

    /* Read MOLECFIT driver file */
    mfdrv drvpar;

    mf_par_initall(&drvpar);
    mf_par_readfile(&drvpar, parfile, 0);


    /* Set all fit flags to 1 */

    molec_flags = cpl_array_new(nmolec, CPL_TYPE_INT);
    cpl_array_fill_window_int(molec_flags, 0, nmolec, 1);

    cont_flags = cpl_array_new(1, CPL_TYPE_INT);
    cpl_array_fill_window_int(cont_flags, 0, 1, 1);

    wlc_flags = cpl_array_new(1, CPL_TYPE_INT);
    cpl_array_fill_window_int(wlc_flags, 0, 1, 1);

    res_flags = cpl_array_new(3, CPL_TYPE_INT);
    cpl_array_fill_window_int(res_flags, 0, 3, 1);

    mf_par_setfitflags(&drvpar, molec_flags, cont_flags, wlc_flags, res_flags);

    cpl_array_delete(molec_flags);
    cpl_array_delete(cont_flags);
    cpl_array_delete(wlc_flags);
    cpl_array_delete(res_flags);


    /* Write MOLECFIT driver file */

    mf_par_writefile(&drvpar, parfile);

    mf_par_deleteall(&drvpar);

    test_invalid();

    test_valid();

    /* Show time */
    double ce = cpl_test_get_walltime();
    cpl_msg_info(cpl_func, "test_mf_par() -> Run time: %g min\n", (ce - cs) / 60.);

    /* Show errors and return */
    cpl_errorstate_dump(errstate, CPL_TRUE, cpl_errorstate_dump_one);
    return !cpl_errorstate_is_equal(errstate);
}
