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
 * molecfit.c
 *
 * Authors:     Stefan Noll & ESO In-Kind Team Innsbruck
 * Created:     13 Dec 2010
 * Last update: 28 Aug 2013
 *
 * MOLECFIT main routine
 *
 * The programme has two input parameters (only the 2nd one is optional):
 * (1) Path and name of input parameter file
 * (2) 'm': multiple runs of CMPFIT (fit schedule) [default]
 *     's': single run of CMPFIT
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <mf_molecfit.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

int main(int argc, char *argv[])
{
    cpl_errorstate errstate = cpl_errorstate_get();
    char mode = 'm';
    char parfile[MF_MAXLEN], errtxt[MF_MAXLEN];
    cpl_boolean err = CPL_FALSE;

    cpl_init(CPL_INIT_DEFAULT);

    /* Read input parameters */
    if (argc > 1) {
        strcpy(parfile, (char *) argv[1]);
        /* Print name of parameter file */
        cpl_msg_info(cpl_func, "Driver file: %s", parfile);
        if (argc > 2) {
            mode = (char) argv[2][0];
            if (mode != 's' && mode != 'm') {
                sprintf(errtxt, "%s: mode ('s' or 'm' only)",
                        MF_ERROR_IIP_TXT);
                cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
                err = CPL_TRUE;
            }
        }
    } else {
        sprintf(errtxt, "%s: no parameter file", MF_ERROR_IIP_TXT);
        cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
        err = CPL_TRUE;
    }

    /* Fit data */
    if (err == CPL_FALSE) {
        mf_molecfit(parfile, mode);
    }

    /* Show errors */
    if (cpl_errorstate_is_equal(errstate)) {
        cpl_msg_info(cpl_func, "No errors occurred");
    } else {
        cpl_errorstate_dump(errstate, CPL_TRUE, cpl_errorstate_dump_one);
    }

    if (cpl_errorstate_is_equal(errstate)) {
        cpl_end();
        return EXIT_SUCCESS;
    } else {
        cpl_end();
        return EXIT_FAILURE;
    }
}

/**@}*/
