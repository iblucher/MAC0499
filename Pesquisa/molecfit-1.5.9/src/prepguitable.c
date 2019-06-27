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

/**@{*/

/*
 * prepguitable.c
 *
 * Authors:     Stefan Noll & ESO In-Kind Team Innsbruck
 * Created:     14 May 2013
 * Last update: 28 Aug 2013
 *
 * Main routine for conversion of input data into a FITS table for the
 * MOLECFIT GUI
 *
 * The programme has one input parameter:
 * (1) Path and name of input parameter file
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <mf_readspec.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

int main(int argc, char *argv[])
{
    cpl_error_code status = CPL_ERROR_NONE;
    cpl_errorstate errstate = cpl_errorstate_get();
    cpl_parameter *p;
    cpl_table *spec;
    mfdrv drvpar;
    mftarr tabdat;
    char parfile_in[MF_MAXLEN], errtxt[MF_MAXLEN], parfile[MF_MAXLEN];
    char basedir[MF_MAXLEN], outdir[MF_MAXLEN], outname[MF_MAXLEN];
    char outfits[MF_MAXLEN];
    int i = 0;
    double weight = 0.;

    cpl_init(CPL_INIT_DEFAULT);

    /* Read input parameters */
    if (argc > 1) {
        strcpy(parfile_in, (char *) argv[1]);
        /* Print name of parameter file */
        cpl_msg_info(cpl_func, "Driver file: %s", parfile_in);
    } else {
        sprintf(errtxt, "%s: no parameter file", MF_ERROR_IIP_TXT);
        cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
        cpl_errorstate_dump(errstate, CPL_TRUE, cpl_errorstate_dump_one);
        cpl_end();
        return EXIT_FAILURE;
    }

    /* Read MOLECFIT driver file */
    mf_par_initall(&drvpar);
    mf_basic_absfile(parfile, parfile_in);
    status = mf_par_readfile(&drvpar, parfile, 0);

    /* Create output table with columns */
    spec = cpl_table_new(0);
    cpl_table_new_column(spec, "chip", CPL_TYPE_INT);
    cpl_table_new_column(spec, "lambda", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "flux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "weight", CPL_TYPE_DOUBLE);

    /* Read input file and convert it to CPL table and CPL property list */
    if (status == CPL_ERROR_NONE) {
        status = mf_conv_readfile(&tabdat, &drvpar);
    }

    /* Write CPL table and CPL property list to FITS table */
    if (status == CPL_ERROR_NONE) {
        mf_conv_writetable(&tabdat, &drvpar);
        mf_conv_tarr_delete(&tabdat);
    }

    /* Read spectroscopic data from previously created FITS table */
    if (status == CPL_ERROR_NONE) {
        status = mf_readspec_fits(spec, &drvpar);
    }

    /* Convert air to vacuum wavelengths if required */
    if (status == CPL_ERROR_NONE) {
        mf_readspec_airtovac(spec, &drvpar);
    }

    /* Convert weight into mask */
    if (status == CPL_ERROR_NONE) {
        cpl_table_new_column(spec, "mask", CPL_TYPE_INT);
        for (i = 0; i < cpl_table_get_nrow(spec); i++) {
            if ((weight = cpl_table_get(spec, "weight", i, NULL)) == 0.) {
                cpl_table_set(spec, "mask", i, 0);
            } else {
                cpl_table_set(spec, "mask", i, 1);
            }
        }
        cpl_table_erase_column(spec, "weight");
    }

    /* Get output file name */
    p = cpl_parameterlist_find(drvpar.parlist, "basedir");
    strncpy(basedir, cpl_parameter_get_string(p), MF_MAXLEN);
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar.parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar.parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar.parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);
    sprintf(outfits, "%s%s_gui.fits", outdir, outname);

    /* Write FITS file */
    cpl_msg_info(cpl_func, "Output FITS table for GUI: %s", outfits);
    cpl_table_save(spec, NULL, NULL, outfits, CPL_IO_CREATE);

    /* Free allocated memory */
    cpl_table_delete(spec);
    mf_par_deleteall(&drvpar);

    /* Show errors */
    if (cpl_errorstate_is_equal(errstate)) {
        cpl_msg_info(cpl_func, "No errors occurred");
    } else {
        cpl_errorstate_dump(errstate, CPL_TRUE, cpl_errorstate_dump_one);
    }

    if (cpl_errorstate_is_equal(errstate)) {
        cpl_end();
        return EXIT_SUCCESS;
    } else {
        cpl_end();
        return EXIT_FAILURE;
    }
}

/**@}*/
