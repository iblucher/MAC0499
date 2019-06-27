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
 * \file mf_corr.h
 *
 * Header for routines related to telluric absorption correction of a set of
 * files by means of a provided transmission curve
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  29 Aug 2012
 * \date   28 Mar 2013
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

#ifndef MF_CORR_H
#define MF_CORR_H

/* Definition of constants */

/*! Minimum transmission criterion for good quality of telluric absorption
    correction */
#define MF_MINTRANS 1e-2
/*! Value for bad mask pixels if not provided by data */
#define MF_MASKVAL_UNDEF 99.

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_corr_filelist(const char *parfile_in);
cpl_error_code mf_corr_readfilelist(cpl_array *filelist, const mfdrv *drvpar);
cpl_error_code mf_corr_ascii(const cpl_array *filelist,
                             const char *reffilename, const mfdrv *drvpar);
cpl_error_code mf_corr_fitstable(const cpl_array *filelist,
                                 const char *reffilename,
                                 const mfdrv *drvpar);
cpl_error_code mf_corr_fitsimage(const cpl_array *filelist,
                                 const char *reffilename,
                                 const char *transfilename,
                                 const mfdrv *drvpar);
cpl_error_code mf_corr_getcolnames(cpl_array *colnames,
                                   const char *filename);
cpl_error_code mf_corr_performtac_tarr(mftarr *tabdat,
                                       const mftarr *reftabdat,
                                       const char *fluxcolname,
                                       const char *dfluxcolname);
cpl_error_code mf_corr_performtac_varr(mfvarr *vecdat, mfvarr *transdat,
                                       const cpl_table *extnames);
cpl_error_code mf_corr_tarr2tactable(cpl_table *spec, const mftarr *tabdat,
                                     const mftarr *reftabdat,
                                     const char *fluxcolname,
                                     const char *dfluxcolname);
cpl_error_code mf_corr_varr2tactable(cpl_table *spec, const mfvarr *vecdat,
                                     const mfvarr *transdat,
                                     const cpl_table *extnames);
cpl_error_code mf_corr_calctac(cpl_table *spec, const int trans);
cpl_error_code mf_corr_tactable2tarr(mftarr *tabdat, const cpl_table *spec);
cpl_error_code mf_corr_tactable2varr(mfvarr *vecdat, mfvarr *transdat,
                                     const cpl_table *extnames,
                                     const cpl_table *spec);
cpl_error_code mf_corr_tacvarr2iarr(mfiarr *imadat, const mfvarr *vecdat,
                                    const mfvarr *transdat,
                                    const cpl_table *extnames);

#endif /* MF_CORR_H */

#ifdef __cplusplus
}
#endif

/**@}*/
