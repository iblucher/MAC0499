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
 * \file mf_lib.h
 *
 * Header for rebinning and flux conversion of wrapper output spectra of the
 * LBLRTM code
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  22 Sep 2009
 * \date   04 Jul 2013
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

#ifndef MF_LIB_H
#define MF_LIB_H

/* Definition of constants */

/*! LBLRTM wrapper output radiance file name */
#define MF_RADNAME "TAPE27"
/*! LBLRTM wrapper output transmission file name */
#define MF_TRANSNAME "TAPE28"

/*! \f${\rm cm}^{-1}\f$ \f$\to\f$ \f$\mu{\rm m}\f$ */
#define MF_CONV_K_LAM 1e4

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_lib_createlibspec(const char *inpath, const char *outpath,
                                    const char *outname, const char *suffix,
                                    const double limlam[2],
                                    const double resol, const mfdrv *drvpar);
cpl_error_code mf_lib_createlamgrid(cpl_table *outspec,
                                    const double limlam[2],
                                    const double resol);
cpl_error_code mf_lib_rebmolspecall(cpl_table *spec, const char *path,
                                    const char *filename);
cpl_error_code mf_lib_findklim(cpl_array *klimall, const char *path,
                               const char *filename);
cpl_error_code mf_lib_rebmolspec(cpl_table *spec, const char *filename,
                                 const double klim[2]);
cpl_error_code mf_lib_interpolspec(cpl_table *spec, const int pixlim[2]);

#endif /* MF_LIB_H */

#ifdef __cplusplus
}
#endif

/**@}*/
