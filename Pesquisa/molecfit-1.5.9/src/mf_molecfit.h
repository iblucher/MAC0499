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
 * \file mf_molecfit.h
 *
 * Header for top-level routines for MOLECFIT
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  05 Nov 2010
 * \date   03 Nov 2014
 */

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *                                INCLUDES                                   *
 ****************************************************************************/

/* Config header */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* MOLECFIT headers */

#include "molecfit.h"
#include <mf_basic.h>
#include <mf_par.h>
#include <mf_readspec.h>
#include <mf_atm.h>
#include <mf_mpfit.h>
#include <mf_lnfl.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef MF_MOLECFIT_H
#define MF_MOLECFIT_H

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

#define MF_OVERRIDE_PAR_STANDALONE(NAME, SET_PAR, VALUE) \
    do { \
        cpl_parameter * par__ = cpl_parameterlist_find(drvpar->parlist, NAME); \
        SET_PAR(par__, VALUE); \
    } while (0)

/* Declaration of functions */

cpl_error_code mf_molecfit(const char *parfile_in, const char mode);
cpl_error_code mf_molecfit_makeoutputdir(const mfdrv *drvpar);
cpl_error_code mf_molecfit_clean(const mfdrv *drvpar);
cpl_error_code mf_molecfit_batch(mp_result *result, cpl_table *spec,
                                 cpl_table *prof, cpl_table** prof_out,
                                 cpl_table** res_out,
                                 mfdrv *drvpar);
cpl_error_code mf_molecfit_calctacfunc(cpl_table *spec, const mfdrv *drvpar);
cpl_error_code mf_molecfit_cleanmodelflux(cpl_table *spec,
                                          const mfdrv *drvpar);
cpl_error_code mf_molecfit_writefile(cpl_table *spec, const mfdrv *drvpar);
cpl_error_code mf_molecfit_plot(cpl_table *spec, const mfdrv *drvpar,
                                const int range);

cpl_error_code
mf_override_config(mfdrv * drvpar,
                   cpl_parameterlist * parlist, cpl_propertylist * plist,
                   cpl_table ** inspec, cpl_size nspec,
                   cpl_table * molectab,
                   cpl_table * wlinclude,
                   cpl_table * wlexclude,
                   cpl_table * pixexclude,
                   cpl_matrix * kernel ,
                   mf_calctrans_state * state
                  );

const char * mf_get_datadir(void);

cpl_error_code fix_directories(mfdrv* drvpar, char** tmpdir);
cpl_error_code mf_cleanup_standalone(char* tmpdir);

#endif /* MF_MOLECFIT_H */

#ifdef __cplusplus
}
#endif

/**@}*/
