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
 * test_mf_atm_writeatm.c
 *
 *  Created on: Jul 21, 2010
 *      Author: barden
 */

/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <cpl_test.h>

#include "mf_atm.h"

/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    /* Initialize */
    cpl_errorstate errstate = cpl_errorstate_get();

    /* Read ATM profile */
    const char *file_in  = "../../data/profiles/mipas/equ.atm";
    cpl_table *atm_profile = cpl_table_new(1);
    if (mf_atm_readatm(atm_profile, file_in) != CPL_ERROR_NONE) {
        cpl_msg_error(cpl_func, "mf_atm_readatm() failed:");
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
    }
    cpl_test_error(CPL_ERROR_NONE);

    /* Write ATM profile */
    const char *file_out = MOLECFIT_SOURCEDIR"/output/equ_test.atm";
    cpl_table* outprof = NULL;
    if (mf_atm_writeatm(atm_profile, &outprof, file_out) != CPL_ERROR_NONE) {
        cpl_msg_error(cpl_func, "mf_atm_writeatm() failed:");
        cpl_errorstate_dump(errstate, CPL_FALSE, cpl_errorstate_dump_one);
    }
    cpl_test_error(CPL_ERROR_NONE);

    /* Cleanup */
    if (outprof) cpl_table_delete(outprof);
    cpl_table_delete(atm_profile);

    /* Show errors */
    cpl_errorstate_dump(errstate, CPL_TRUE, cpl_errorstate_dump_one);
    if (cpl_errorstate_is_equal(errstate)) return EXIT_SUCCESS;
    else                                   return EXIT_FAILURE;
}

/**@}*/
