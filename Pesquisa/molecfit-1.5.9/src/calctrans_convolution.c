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
 * calctrans_convolution.c
 *
 * Authors:     Dr. Gurvan Bazin, ESO
 * Created:     2017-04-04
 * Last update: 2017-04-04
 *
 * Main routine for convolving the spectra.
 *
 * The programme has one input parameter:
 * (1) Path and name of input parameter file
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <mf_trans.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

int main(int argc, char *argv[])
{
    cpl_errorstate errstate = cpl_errorstate_get();
    char parfile[MF_MAXLEN], errtxt[MF_MAXLEN];
    cpl_boolean err = CPL_FALSE;

    cpl_init(CPL_INIT_DEFAULT);

    /* Read input parameters */
    if (argc > 1) {
        strcpy(parfile, (char *) argv[1]);
        /* Print name of parameter file */
        cpl_msg_info(cpl_func, "Driver file: %s", parfile);
    } else {
        sprintf(errtxt, "%s: no parameter file", MF_ERROR_IIP_TXT);
        cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
        err = CPL_TRUE;
    }

    /* Calculate transmission curve and correct input data file */
    if (err == CPL_FALSE) {
        mf_trans_conv_(parfile);
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
