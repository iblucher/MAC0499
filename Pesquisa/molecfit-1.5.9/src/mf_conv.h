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
 * \file mf_conv.h
 *
 * Header for routines related to file conversion
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  25 Jul 2012
 * \date   24 Jan 2014
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

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef MF_CONV_H
#define MF_CONV_H

/* Definition of constants */

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/*!
 * Structure for data of FITS tables
 *
 * \param next    number of FITS extensions
 * \param **head  array of CPL property lists of size \e next + 1
 * \param **tab   array of CPL tables of size \e next + 1
 */

typedef struct _mftarr_ {
    int next;
    cpl_propertylist **head;
    cpl_table **tab;
} mftarr;

/*!
 * Structure for data of 1D FITS images
 *
 * \param next    number of FITS extensions
 * \param **head  array of CPL property lists of size \e next + 1
 * \param **vec   array of CPL vectors of size \e next + 1
 */

typedef struct _mfvarr_ {
    int next;
    cpl_propertylist **head;
    cpl_vector **vec;
} mfvarr;

/*!
 * Structure for data of 2D FITS images
 *
 * \param next    number of FITS extensions
 * \param **head  array of CPL property lists of size \e next + 1
 * \param **ima   array of CPL images of size \e next + 1
 */

typedef struct _mfiarr_ {
    int next;
    cpl_propertylist **head;
    cpl_image **ima;
} mfiarr;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_conv_preptable(const char *parfile_in);
cpl_error_code mf_conv_readfile(mftarr *tabdat, const mfdrv *drvpar);
cpl_error_code mf_conv_checkfitsformat(int *fitsformat, const char *filename);
cpl_error_code mf_conv_setcolnames(cpl_array *colnames, const mfdrv *drvpar);
cpl_error_code mf_conv_setextnames(cpl_table *extnames, const mfvarr *vecdat,
                                   const mfdrv *drvpar);
cpl_error_code mf_conv_varr2tarr(mftarr *tabdat, const mfvarr *vecdat,
                                 const cpl_table *extnames);
double mf_conv_getwcskey(const cpl_propertylist *plist, const char *key);
cpl_error_code mf_conv_modtable(mftarr *tabdat, const mfdrv *drvpar);
cpl_error_code mf_conv_writetable(const mftarr *tabdat, const mfdrv *drvpar);
cpl_error_code mf_conv_writeresults(const char *parfile_in);
cpl_error_code mf_conv_readresults(cpl_table *results, const mfdrv *drvpar);
cpl_error_code mf_conv_writefile(cpl_table *results, const mfdrv *drvpar);
cpl_error_code mf_conv_readprepfits(mftarr *tabdat, const mfdrv *drvpar);
cpl_error_code mf_conv_results2tarr(mftarr *tabdat, cpl_table *results,
                                    const mfdrv *drvpar);
cpl_error_code mf_conv_erasemaskcol(mftarr *tabdat, const mfdrv *drvpar);
cpl_error_code mf_conv_resultstarr2varr(mfvarr *transdat, mfvarr *vecdat,
                                        const cpl_table *extnames,
                                        const mftarr *tabdat);
cpl_error_code mf_conv_getmaskval(double maskval[2], const mftarr *tabdat,
                                  const cpl_table *extnames);
cpl_error_code mf_conv_iarr2varr(mfvarr *vecdat, const mfiarr *imadat,
                                 const cpl_table *extnames);
cpl_error_code mf_conv_ascii_read(mftarr *tabdat, const char *filename,
                                  const cpl_array *colnames);
cpl_error_code mf_conv_ascii_write(const char *filename,
                                   const mftarr *tabdat);
cpl_error_code mf_conv_tarr_init(mftarr *tabdat, const int next);
cpl_error_code mf_conv_tarr_read(mftarr *tabdat, const char *filename);
cpl_error_code mf_conv_tarr_write(const char *filename, const mftarr *tabdat);
cpl_error_code mf_conv_tarr_delete(mftarr *tabdat);
cpl_error_code mf_conv_varr_init(mfvarr *vecdat, const int next);
cpl_error_code mf_conv_varr_read(mfvarr *vecdat, const char *filename);
cpl_error_code mf_conv_varr_write(const char *filename, const mfvarr *vecdat);
cpl_error_code mf_conv_varr_delete(mfvarr *vecdat);
cpl_error_code mf_conv_iarr_init(mfiarr *imadat, const int next);
cpl_error_code mf_conv_iarr_read(mfiarr *imadat, const char *filename);
cpl_error_code mf_conv_iarr_write(const char *filename, const mfiarr *imadat);
cpl_error_code mf_conv_iarr_delete(mfiarr *imadat);

#endif /* MF_CONV_H */

#ifdef __cplusplus
}
#endif

/**@}*/
