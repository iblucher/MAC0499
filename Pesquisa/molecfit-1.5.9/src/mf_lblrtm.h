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
 * \file mf_lblrtm.h
 *
 * Header for routines related to the running of the radiative transfer code
 * LBLRTM
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 * \since  20 Jun 2012
 * \date   21 Jul 2013
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

#include <mf_basic.h>
#include <mf_par.h>
#include <mf_lib.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef MF_LBLRTM_H
#define MF_LBLRTM_H

/* Definition of constants */

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_lblrtm_call(const mfdrv *drvpar, const cpl_table *prof,
                              const int range);
cpl_error_code mf_lblrtm_readsetup(cpl_parameterlist *plist,
                                   const mfdrv *drvpar);
cpl_error_code mf_lblrtm_start1(const mfdrv *drvpar,
                                cpl_parameterlist *lblrtm_setup,
                                const cpl_table *prof, double V1, double V2,
                                const int range);
cpl_array *mf_lblrtm_allmolecs(void);
cpl_error_code mf_lblrtm_renameoutput(const char *wdir, const int a,
                                      const char mode);

#endif /* MF_LBLRTM_H */

#ifdef __cplusplus
}
#endif

/**@}*/
