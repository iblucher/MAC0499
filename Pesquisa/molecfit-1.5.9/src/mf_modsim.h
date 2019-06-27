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
 * \file mf_modsim.h
 *
 * Header for routines related to the preparation of model spectra required
 * for the mf_mpfit fitting procedure
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  06 Jul 2010
 * \date   01 Nov 2014
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
#include <mf_atm.h>
#include <mf_lblrtm.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef MF_MODSIM_H
#define MF_MODSIM_H

/* Definition of constants */

/*! Threshold for recalculation of kernel depending on relative change of
    wavelength */
#define MF_LIMRELLAMVAR 0.01
/*! for integral of Gaussian and Lorentzian */
#define MF_BINS_PER_FWHM 200
/*! \f$\mu{\rm m}\f$ */
#define MF_LAM_UNIT 1e-6
/*! steradians in \f${\rm arcsec}^2\f$ */
#define MF_SR_IN_ARCSEC2 4.254517e+10

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_modsim_lblrtm(const cpl_table *prof,
                                const mfdrv *drvpar, const cpl_array *fitpar,
                                cpl_error_code** arr_code_stat);
cpl_error_code mf_modsim_conv(cpl_table *spec,
                              const mfdrv *drvpar, const cpl_array *fitpar,
                              cpl_error_code* arr_code_stat);
cpl_error_code mf_modsim(cpl_table *spec, cpl_table **prof_out,
                         const cpl_table *prof,
                         const mfdrv *drvpar, const cpl_array *fitpar);
cpl_error_code mf_modsim_compfitpar(int *ismolcol, const mfdrv *drvpar,
                                    const cpl_array *fitpar);
cpl_error_code mf_modsim_modatm(cpl_table *modprof, cpl_table** modprof_out,
                                const mfdrv *drvpar,
                                const cpl_array *fitpar);
cpl_error_code mf_modsim_readlblspec(cpl_table *modspec, const mfdrv *drvpar,
                                     const int range);
cpl_error_code mf_modsim_modwavegrid(cpl_table *spec, const mfdrv *drvpar,
                                     const cpl_array *fitpar, const int chip);
cpl_error_code mf_modsim_convolvesynthkernel(cpl_table *spec,
                                             const mfdrv *drvpar,
                                             const cpl_array *fitpar);
cpl_error_code mf_modsim_calckerneltable(cpl_table *kerntab,
                                         const cpl_table *spec,
                                         const mfdrv *drvpar,
                                         const cpl_array *fitpar);
cpl_error_code mf_modsim_calcsynthkernel(cpl_array *kernel,
                                         const double wbox,
                                         const double wgauss,
                                         const double wlorentz,
                                         const double kernfac,
                                         const int kernmode);
cpl_error_code mf_modsim_calcboxkernel(cpl_array *kernel, const double fwhm);
cpl_error_code mf_modsim_calcvoigtkernel(cpl_array *kernel,
                                         const double wgauss,
                                         const double wlorentz,
                                         const double kernfac);
cpl_error_code mf_modsim_calcgausskernel(cpl_array *kernel,
                                         const double fwhm,
                                         const double kernfac);
cpl_error_code mf_modsim_calclorentzkernel(cpl_array *kernel,
                                           const double fwhm,
                                           const double kernfac);
cpl_error_code mf_modsim_invertkerneltable(cpl_table *kerntab);
cpl_error_code mf_modsim_convolvekerneltable(cpl_table *spec,
                                             const cpl_table *kerntab);
cpl_error_code mf_modsim_convolvereadkernel(cpl_table *spec,
                                            const cpl_array *specrows,
                                            const mfdrv *drvpar);
cpl_error_code mf_modsim_telback(cpl_table *spec, const mfdrv *drvpar,
                                 const cpl_array *fitpar);
cpl_error_code mf_modsim_fluxunits(cpl_table *spec, const mfdrv *drvpar);
cpl_error_code mf_modsim_modcont(cpl_table *spec, const mfdrv *drvpar,
                                 const cpl_array *fitpar, const int range);
cpl_error_code mf_modsim_corrobsspec(cpl_table *spec, const mfdrv *drvpar);

#endif /* MF_MODSIM_H */

#ifdef __cplusplus
}
#endif

/**@}*/
