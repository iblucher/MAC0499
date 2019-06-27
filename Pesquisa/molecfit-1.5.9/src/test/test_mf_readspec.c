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

/**@{*/

/*
 * test_mf_readspec.c
 *
 * Authors:     Stefan Noll & ESO In-Kind Team Innsbruck
 * Created:     22 Jun 2010
 * Last update: 25 Sep 2012
 *
 * Test programme for mf_readspec.c
 */

#include <mf_readspec.h>

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    /* Initialize */
    cpl_errorstate errstate = cpl_errorstate_get();
    double         cs       = cpl_test_get_cputime();

    /* Read MOLECFIT driver file */
    mfdrv drvpar;
    const char *parfile = MOLECFIT_SOURCEDIR"/config/molecfit_test_crires.par";
    mf_par_initall(&drvpar);
    mf_par_readfile(&drvpar, parfile, 0);

    /* Read spectrum and header data from data file */
    cpl_table *spec = cpl_table_new(0);
    mf_readspec(spec, &drvpar);

    /* Write content of CPL table "spec" to file in "output/" */
    const char *outfile = MOLECFIT_SOURCEDIR"/output/spec.dat";
    cpl_msg_info(cpl_func, "%s", outfile);
    FILE *stream = fopen(outfile, "w+");
    if (stream) {
        cpl_table_dump(spec, 0, cpl_table_get_nrow(spec), stream);
        fclose(stream);
    }

    /* Write updated driver file parameters to file in "output/" */
    mf_par_writefile(&drvpar, parfile);
    cpl_test_error(CPL_ERROR_NONE);

    /* cleanup */
    cpl_table_delete(spec);
    mf_par_deleteall(&drvpar);
    cpl_test_error(CPL_ERROR_NONE);


    /* Show execution time */
    double ce = cpl_test_get_walltime();
    cpl_msg_info(cpl_func, "test_mf_readspec -> Run time: %g s\n", ce - cs);


    /* Show errors and return */
    cpl_errorstate_dump(errstate, CPL_TRUE, cpl_errorstate_dump_one);
    if (cpl_errorstate_is_equal(errstate)) return EXIT_SUCCESS;
    else                                   return EXIT_FAILURE;
}

/**@}*/
