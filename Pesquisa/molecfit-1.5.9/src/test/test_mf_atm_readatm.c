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
 * test_mf_atm_readatm.c
 *
 *  Created on: Jul 5, 2010
 *      Author: barden
 */

#include "mf_atm.h"

#include <cpl_test.h>


int main(void) {

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    /* Initialize */
    cpl_errorstate errstate = cpl_errorstate_get();

    /* Get time */
    double cs = cpl_test_get_walltime();

    /* Read standard profile */
    const char *filename = "../../data/profiles/mipas/equ.atm";
    cpl_table  *atm_profile = cpl_table_new(1);
    mf_atm_readatm(atm_profile, filename);
    cpl_test_error(CPL_ERROR_NONE);

    /* Cleanup*/
    cpl_table_delete(atm_profile);
    cpl_test_error(CPL_ERROR_NONE);

    /* Show time */
    double ce = cpl_test_get_walltime();
    cpl_msg_info(cpl_func, "test_mf_atm_readatm() -> Run time: %g min\n", (ce - cs) / 60.);

    /* Show errors and return */
    cpl_errorstate_dump(errstate, CPL_TRUE, cpl_errorstate_dump_one);
    return !cpl_errorstate_is_equal(errstate);
}
