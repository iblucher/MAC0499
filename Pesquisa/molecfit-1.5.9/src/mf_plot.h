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


/*!
 * \file mf_plot.h
 *
 * \brief Header for plotting library
 *
 * \author Wolfgang Kausch & ESO In-Kind Team Innsbruck
 */

/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <mf_basic.h>
#include <mf_par.h>

/*****************************************************************************
 *                                 DEFINES                                   *
 ****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef MF_PLOT_H
#define MF_PLOT_H

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *                                 TYPEDEF                                   *
 ****************************************************************************/

/*****************************************************************************
 *                                 GLOBALS                                   *
 ****************************************************************************/

/*****************************************************************************
 *                                 PROTOTYPES                                *
 ****************************************************************************/

  cpl_error_code mf_plot_single(
      const cpl_table *spec,
      const char      *user_x_label,
      const char      *user_y_label,
      const char      *plot_title,
      const mfdrv     *drvpar);

cpl_error_code mf_plot_double(const cpl_table *spec1, const cpl_table *spec2,
                              const mfdrv *drvpar);

cpl_error_code mf_plot_xy(
    const cpl_table *spec,
    const char      *user_x_label,
    const char      *user_y_label,
    const char      *plot_title,
    const mfdrv     *drvpar);

cpl_error_code mf_plot_hist(const cpl_table *histdat, char *x_col,
                            char *y_col, char *plot_title,
                            const mfdrv *drvpar);

#ifdef __cplusplus
}
#endif

#endif
