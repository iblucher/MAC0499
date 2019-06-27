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
 * \file mf_par.h
 *
 * Header for MOLECFIT driver file routines
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  27 Apr 2010
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

#include <mf_basic.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef MF_PAR_H
#define MF_PAR_H

/* Definition of constants */

/*! Maximum parameter tag length */
#define MF_MAXTAGLEN 20
/*! Maximum number of lines in MOLECFIT driver file */
#define MF_MAXNLINE 1000

/* Default column names */

/*! Default wavelength column name */
#define MF_DEFLAMCOL "LAMBDA"
/*! Default flux column name */
#define MF_DEFFLUXCOL "FLUX"
/*! Default flux error column name */
#define MF_DEFDFLUXCOL "DFLUX"
/*! Default mask column name */
#define MF_DEFMASKCOL "MASK"

/* Parameters not provided by the driver file */

/*! Minimum expected wavelength in \f$\mu\f$m, warning emitted for smaller wavelength */
#define MF_WAVELENGTH_MIN 0.28
/*! Maximum expected wavelength in \f$\mu\f$m, warning emitted for larger wavelength */
#define MF_WAVELENGTH_MAX 8600.
/*! Extra wavenumber coverage at both sides of the model spectrum in
    \f${\rm cm}^{-1}\f$ */
#define MF_EXTRACOVER 5
/*! Oversampling factor for model spectrum */
#define MF_SAMPFAC 5

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/* Definition of structures */

/*!
 * Structure for storing the parameters of the MOLECFIT driver file
 *
 * \param parlist   CPL parameter list
 * \param chiptab   CPL table for chip-related parameters
 * \param molectab  CPL table for parameters related to molecules
 * \param kernel    CPL matrix for elements of optional kernel for each pixel
 *                  (separate file)
 */

typedef struct _mfdrv_ {
    cpl_parameterlist *parlist;
    cpl_table *rangetab;
    cpl_table *chiptab;
    cpl_table *molectab;
    cpl_matrix *kernel;
} mfdrv;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

/* Basic IO functions for MOLECFIT driver file */

cpl_error_code mf_par_initall(mfdrv *drvpar);
cpl_error_code mf_par_readfile(mfdrv *drvpar, const char *parfile, char isCalctrans);
cpl_error_code mf_par_copyfile(const mfdrv *drvpar, const char *parfile);
cpl_error_code mf_par_writefile(const mfdrv *drvpar, const char *parfile);
cpl_error_code mf_par_deleteall(mfdrv *drvpar);
cpl_error_code mf_par_search(mfpar par[], int *npar, const char *parname,
                             const char *parfile);

/* Special functions for manipulating the MOLECFIT parameter structure */

cpl_error_code mf_par_setfitflags(mfdrv *drvpar, const cpl_array *fit_molec,
                                  const cpl_array *fit_cont,
                                  const cpl_array *fit_wlc,
                                  const cpl_array *fit_res);
cpl_error_code mf_par_finalize_lbl_molecs(mfdrv * drvpar);
void mf_par_setmolec_all(cpl_parameter *p, const char *mol);
void mf_par_setmolec(cpl_parameter *p, const int pos);

/* Functions related to user-provided kernels */

cpl_error_code mf_par_readkernel(mfdrv *drvpar);

#endif /* MF_PAR_H */

#ifdef __cplusplus
}
#endif

/**@}*/
