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
 * \file mf_lnfl.h
 *
 * Header for running LNFL (line file creation)
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 * \since  20 Jun 2012
 * \date   29 Jun 2012
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

#ifndef MF_LNFL_H
#define MF_LNFL_H

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_lnfl(const mfdrv *drvpar);
cpl_error_code mf_lnfl_call(const mfdrv *drvpar, const int range);
cpl_error_code mf_lnfl_files(char *tape3_name, char *tape5_name,
                             const char *dir, const double wn_start,
                             const double wn_end, const char *molecs);

#endif /* MF_LNFL_H */

#ifdef __cplusplus
}
#endif

/**@}*/
