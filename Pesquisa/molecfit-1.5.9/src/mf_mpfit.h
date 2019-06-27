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
 * \file mf_mpfit.h
 *
 * Header for routines related to the handling of CMPFIT
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  25 Jun 2010
 * \date   11 Sep 2013
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

/* Use of CMPFIT */

#include <mpfit.h>

/* MOLECFIT headers */

#include <mf_basic.h>
#include <mf_par.h>
#include <mf_modsim.h>
#include <mf_atm.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef MF_MPFIT_H
#define MF_MPFIT_H

/* Definition of constants */

/*! maximum number of iterations per CMPFIT run */
#define MF_MAXITER 100

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/* Definition of structures */

/*!
 * Structure for holding the spectral data and the driver file parameters
 *
 * \param spec     CPL table with wavelengths, fluxes, and weights
 * \param prof     CPL table with atmospheric profiles
 * \param drvpar   ::mfdrv parameter structure
 */

typedef struct _mfvars_ {
    cpl_table *spec;
    cpl_table *prof;
    mfdrv *drvpar;
} mfvars;

/*!
 * Structure for holding the fit parameters and their constraints
 *
 * \param n     number of parameters
 * \param p     parameter vector for fitting
 * \param pars  structure for parameter constraints
 */

typedef struct _mfpars_ {
    int n;
    double *p;
    mp_par *pars;
} mfpars;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_mpfit(mp_result *result, cpl_table *spec, cpl_table *prof,
                        cpl_table** prof_out, cpl_table** res_out, mfdrv *drvpar);
int mf_mpfit_calcdev(int m, int n, double *p, double *dy, double **dvec,
                     void *vars, cpl_table** prof_out);
cpl_error_code mf_mpfit_setpar(mfpars *fitpar, const mfdrv *drvpar);
cpl_error_code mf_mpfit_getpar(mfdrv *drvpar, const mfpars *fitpar);
cpl_error_code mf_mpfit_allocmempar(mfpars *fitpar, const int npar);
cpl_error_code mf_mpfit_freemempar(mfpars *fitpar);
cpl_error_code mf_mpfit_allocmemresult(mp_result *result, const int m,
                                       const int n);
cpl_error_code mf_mpfit_freememresult(mp_result *result);
cpl_error_code mf_mpfit_writeresults(mfdrv *drvpar, const cpl_table *spec,
                                     const mp_result *result,
                                     const double fittime,
                                     cpl_table** outtable);

#endif /* MF_MPFIT_H */

#ifdef __cplusplus
}
#endif

/**@}*/
