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
 * \file mf_trans.h
 *
 * Header related to routines for deriving a transmission curve for telluric
 * absorption correction
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  10 Jul 2012
 * \date   05 Jun 2013
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
#include <mf_readspec.h>
#include <mf_atm.h>
#include <mf_lnfl.h>
#include <mf_modsim.h>
#include <mf_conv.h>
#include <mf_corr.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef MF_TRANS_H
#define MF_TRANS_H

/* Definition of constants */

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_trans(const char *parfile_in);
cpl_error_code mf_trans_readspec(cpl_table *spec, mfdrv *drvpar);
cpl_error_code mf_trans_readatm(cpl_table **prof, const mfdrv *drvpar);
cpl_error_code mf_trans_readresults_fromFits(cpl_array *fitpar, const mfdrv *drvpar, cpl_table* res_table);
cpl_error_code mf_trans_readresults(cpl_array *fitpar, const mfdrv *drvpar);
cpl_error_code mf_trans_clean(const mfdrv *drvpar, const char cleanSpec);
cpl_error_code mf_trans_calctac(cpl_table *spec, const mfdrv *drvpar);
cpl_error_code mf_trans_writefile(cpl_table *spec, const mfdrv *drvpar);
cpl_error_code mf_trans_plot(const cpl_table *spec, const mfdrv *drvpar);
cpl_error_code mf_trans_writeresults(mfdrv *drvpar);
cpl_error_code mf_trans_lblrtm_(const char *parfile_in);
cpl_error_code mf_trans_conv_(const char *parfile_in);

#endif /* MF_TRANS_H */

#ifdef __cplusplus
}
#endif

/**@}*/
