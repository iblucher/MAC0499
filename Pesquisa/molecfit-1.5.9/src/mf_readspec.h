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
 * \file mf_readspec.h
 *
 * Header for routines related to the reading of observing data
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  09 Jun 2010
 * \date   26 Jan 2014
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
#include <mf_conv.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef MF_READSPEC_H
#define MF_READSPEC_H

/* Definition of constants */

/*! Number of FITS header keywords defined in mf_readspec_header */
#define MF_NKEY 12

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_readspec(cpl_table *spec, mfdrv *drvpar);
cpl_error_code mf_readspec_preptable(mfdrv *drvpar);
cpl_error_code mf_readspec_fits(cpl_table *spec, mfdrv *drvpar);
cpl_error_code mf_readspec_header(mfdrv *drvpar);
cpl_error_code mf_readspec_slitwidth_xshooter(double *slitw,
                                              const cpl_propertylist *header);
cpl_error_code mf_readspec_airtovac(cpl_table *spec, const mfdrv *drvpar);
cpl_error_code mf_readspec_fillpartab(mfdrv *drvpar, cpl_table *spec);
cpl_error_code mf_readspec_cutspec(cpl_table *spec, const mfdrv *drvpar);
cpl_error_code mf_readspec_includeranges(cpl_table *spec, mfdrv *drvpar);
cpl_error_code mf_readspec_replacecoef(mfdrv *drvpar, cpl_table *rangetab,
                                       cpl_table *chiptab);
cpl_error_code mf_readspec_excluderanges(cpl_table *spec,
                                         const mfdrv *drvpar, const char wp);
cpl_error_code mf_readspec_modranges(mfdrv *drvpar, const cpl_table *spec);
cpl_error_code mf_readspec_setweight(cpl_table *spec, const mfdrv *drvpar);
cpl_error_code mf_readspec_calcwaverange(mfdrv *drvpar, cpl_table *spec);

#endif /* MF_READSPEC_H */

#ifdef __cplusplus
}
#endif

/**@}*/
