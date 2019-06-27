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
 * \file mf_atm.h
 *
 * Headers for routines related to derivation and handling of atmospheric
 * profiles
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 * \since  20 Jun 2012
 * \date   23 Jan 2014
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
#include <mf_lblrtm.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef MF_ATM_H
#define MF_ATM_H

/* Definition of constants */

/*! Height fraction relative to highest valid GDAS layer used for
    interpolating between GDAS and standard profile */
#define MF_MERGEFRAC 0.2
/*! Gas constant in J/(mol K) */
#define MF_R 8.3144510

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_atm_createatm(cpl_table *out_profile, const mfdrv *drvpar);
cpl_error_code mf_atm_getgdas_auto(cpl_table **gdas_profile1,
                                   cpl_table **gdas_profile2,
                                   cpl_array *timestamps,
                                   const mfdrv *drvpar);
cpl_error_code mf_atm_getgdas_user(cpl_table *gdas_profile,
                                   const mfdrv *drvpar);
void mf_atm_concatgdasname(char *gdasname, const char *gdasdir,
                           const double lon, const double lat,
                           const char *date, const char *hour);
cpl_error_code mf_atm_readatm_fromFits(cpl_table **atm_profile, const char *atm_file);
cpl_error_code mf_atm_readatm(cpl_table *atm_profile, const char *atm_file);
cpl_error_code mf_atm_readgdas(cpl_table *gdas_profile, const char *gdas_file,
                               const char *hgt_units);
cpl_error_code mf_atm_interpolprofile(cpl_table *outprofile,
                                      const cpl_table *profile1,
                                      const cpl_table *profile2,
                                      const cpl_array *timestamps);
cpl_error_code mf_atm_convertgdas(cpl_table *merged_profile,
                                  const cpl_table *atm_profile,
                                  const cpl_table *gdas_profile,
                                  const cpl_array *molecs,
                                  const mfdrv *drvpar);
cpl_error_code mf_atm_convertgdas_fixed(cpl_table *merged_profile,
                                        const cpl_table *atm_profile,
                                        const cpl_table *gdas_profile,
                                        const cpl_array *molecs,
                                        const double geoelev);
cpl_error_code mf_atm_mergelayers(cpl_table *merged_profile,
                                  const cpl_table *atm_profile,
                                  const cpl_table *gdas_profile,
                                  const mfdrv *drvpar);
cpl_error_code mf_atm_mergegdas(cpl_table *merged_profile,
                                const cpl_table *atm_profile,
                                const cpl_table *gdas_profile,
                                const cpl_array *molecs);
cpl_error_code mf_atm_adaptenv(cpl_table *profile, const mfdrv *drvpar);
void mf_atm_adaptenv_basic(cpl_table *profile, const double gelev,
                            const double gpres, const double gtemp,
                            const double ghum, const double emix);
double mf_atm_interpollog(const double *x, cpl_table *profile,
                          const int start, const char *varstr);
cpl_error_code mf_atm_writeatm(cpl_table *atm_profile,
                               cpl_table **atm_profile_out,
                               const char *atm_file);
cpl_error_code mf_atm_scaletopwv(cpl_table *profile, const mfdrv *drvpar);
cpl_error_code mf_atm_calpwv(double *pwv, cpl_table *prof,
                             const double *geoelev);
cpl_error_code mf_atm_calpwv_histo(double *pwv, cpl_table *prof,
                                   const double *geoelev);
cpl_error_code mf_atm_calcol(double *h2ocol, mfdrv *drvpar);
cpl_error_code mf_atm_calcol_histo(double *h2ocol, mfdrv *drvpar);

#endif /* MF_ATM_H */

#ifdef __cplusplus
}
#endif

/**@}*/
