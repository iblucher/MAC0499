/*
 *  This file is part of the MOLECFIT software package.
 *  Copyright (C) 2016 European Southern Observatory
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
 * \file molecfit.h
 *
 * Header for top-level exported routines for MOLECFIT
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *                                INCLUDES                                   *
 ****************************************************************************/

/* MOLECFIT headers */
#include <cpl.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef INCL_MOLECFIT_H
#define INCL_MOLECFIT_H

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

/* opaque state structure for molecfit -> calctrans runs */
typedef struct mf_calctrans_state_ mf_calctrans_state;

mf_calctrans_state* mf_init_calctrans_state(void);

cpl_error_code mf_run_molecfit(
    cpl_parameterlist  *parlist,
    cpl_propertylist   *plist, 
    cpl_table          **inspec, 
    cpl_size           nspec,
    cpl_table          *molectab,
    cpl_table          *wlinclude,
    cpl_table          *wlexclude,
    cpl_table          *pixexclude,
    cpl_matrix         *kernel,
    cpl_table          **prof_out,
    cpl_table          **res_out,
    cpl_table          **spec_out,
    mf_calctrans_state **state);

cpl_error_code mf_run_calctrans(
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
    cpl_table          **result);

cpl_error_code mf_run_calctrans_lblrtm(
    cpl_parameterlist  *parlist,
    cpl_propertylist   *plist,
    cpl_table          **inspec,
    cpl_size           nspec,
    cpl_table          *molectab,
    double             wl_start,
    double             wl_end,
    cpl_table          *atmprof,
    cpl_table          *res_table,
    mf_calctrans_state **state);

cpl_error_code mf_run_calctrans_convolution(
    cpl_parameterlist  *parlist,
    cpl_propertylist   *plist,
    cpl_table          **inspec,
    cpl_size           nspec,
    cpl_matrix         *kernel,
    double             wl_start,
    double             wl_end,
    cpl_table          *atmprof,
    cpl_table          *res_table,
    mf_calctrans_state **state,
    cpl_table          **result);

cpl_error_code clean_arr_code_stat(
    cpl_parameterlist  *parlist);

cpl_error_code mf_cleanup(
    mf_calctrans_state *state);


#endif /* INCL_MOLECFIT_H */

#ifdef __cplusplus
}
#endif

/**@}*/
