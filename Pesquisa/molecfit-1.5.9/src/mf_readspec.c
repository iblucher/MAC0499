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
 * \file mf_readspec.c
 *
 * Routines for reading observing data
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  09 Jun 2010
 * \date   25 Mar 2015
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <mf_readspec.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code mf_readspec(cpl_table *spec, mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Reads a FITS table that was created by PREPTABLE and puts the read data
     * in a CPL table consisting of the columns "lambda", "flux", and
     * "weight". An input FITS file can have an arbitrary number of extensions
     * with different parts of the spectrum (to be defined in the MOLECFIT
     * driver file). If the hierarchical ESO keywords are given in the header
     * of the provided FITS file, the required parameters describing the
     * telescope site and observing conditions are taken from the header.
     * Their values substitute those given in the MOLECFIT driver file.
     * Moreover, files consisting of wavelength ranges to be included or
     * wavelength or pixel ranges to be excluded are optionally read. These
     * data determine the model range as given by the "mrange" column and the
     * range table of the ::mfdrv parameter structure.
     *
     * \b INPUT:
     * \param spec    empty CPL table
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    CPL table with wavelengths, fluxes, and weights
     * \param drvpar  ::mfdrv parameter structure with FITS header values
     *                if the observing instrument is given
     *
     * \b ERRORS:
     * - No data
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_table *rangetab = NULL, *chiptab = NULL;
    mftarr tabdat;
    char errtxt[MF_MAXLEN];

    /* Create output table columns */
    cpl_table_set_size(spec, 0);
    cpl_table_new_column(spec, "chip", CPL_TYPE_INT);
    cpl_table_new_column(spec, "lambda", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "flux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "weight", CPL_TYPE_DOUBLE);
    //cpl_table_new_column(spec, "iscont", CPL_TYPE_INT);
    cpl_table_new_column(spec, "mrange", CPL_TYPE_INT);
    cpl_table_new_column(spec, "mlambda", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "mscal", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "mflux", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "mweight", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "dev", CPL_TYPE_DOUBLE);

    /* Read input file and convert it to CPL table and CPL property list */
    if ((status = mf_conv_readfile(&tabdat, drvpar)) != CPL_ERROR_NONE) {
        return status;
    }
    
    /* Write CPL table and CPL property list to FITS table */
    mf_conv_writetable(&tabdat, drvpar);

    /* Free allocated memory */
    mf_conv_tarr_delete(&tabdat);

    /* Exit in case of errors */
    if ((status = cpl_error_get_code()) != CPL_ERROR_NONE) {
        return status;
    }

    /* Get spectroscopic data and observing parameters (if available) */
    if (cpl_error_get_code() == CPL_ERROR_NONE) {
        mf_readspec_fits(spec, drvpar);
    }
    if (cpl_error_get_code() == CPL_ERROR_NONE) {
        mf_readspec_header(drvpar);
    }

    /* Exit in case of errors */
    if ((status = cpl_error_get_code()) != CPL_ERROR_NONE) {
        return status;
    }

    /* No data points? */
    if (cpl_table_get_nrow(spec) == 0) {
        sprintf(errtxt, "%s: cpl_table *spec", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Convert air to vacuum wavelengths if required */
    mf_readspec_airtovac(spec, drvpar);

    /* Save read parameter tables before manipulation */
    rangetab = cpl_table_duplicate(drvpar->rangetab);
    chiptab = cpl_table_duplicate(drvpar->chiptab);

    /* Fill the parameter tables dependent on chip number */
    mf_readspec_fillpartab(drvpar, spec);

    /* Optional use of range table to select wavelength ranges for the fit */
    if ((status = mf_readspec_includeranges(spec, drvpar)) !=
        CPL_ERROR_NONE) {
        cpl_table_delete(rangetab);
        cpl_table_delete(chiptab);
        return status;
    }

    /* Replace default parameter values in the range and chip tables by data
       from the parameter file */
    mf_readspec_replacecoef(drvpar, rangetab, chiptab);

    /* Delete temporary tables */
    cpl_table_delete(rangetab);
    cpl_table_delete(chiptab);

    /* Optional use of range table to exclude wavelength ranges from the
       fit */
    if ((status = mf_readspec_excluderanges(spec, drvpar, 'w')) !=
        CPL_ERROR_NONE) {
        return status;
    }

    /* Optional use of range table to exclude pixel ranges from the fit */
    if ((status = mf_readspec_excluderanges(spec, drvpar, 'p')) !=
        CPL_ERROR_NONE) {
        return status;
    }

    /* Modify ranges to calculate only a single model spectrum? */
    mf_readspec_modranges(drvpar, spec);

    /* Set weights for spectra without error column */
    mf_readspec_setweight(spec, drvpar);

    /* All weights = 0? */
    if (cpl_table_get_column_max(spec, "weight") == 0) {
        sprintf(errtxt, "%s: cpl_table *spec (all weights = 0)",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Calculate required wavenumber ranges and store it in mfdrv structure */
    if ((status = mf_readspec_calcwaverange(drvpar, spec) !=
         CPL_ERROR_NONE)) {
        return status;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_preptable(mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Converts input file (ASCII, FITS image, or FITS table) into a FITS
     * table for MOLECFIT. The conversion is skipped if the output file
     * already exists.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b ERRORS:
     * - see subroutines
     */

    FILE *stream;
    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    mfpar par[MF_MAXPAR];
    mftarr tabdat;
    char filename[MF_MAXLEN], outdir[MF_MAXLEN];
    char outname[MF_MAXLEN], outfile[MF_MAXLEN], resfile[MF_MAXLEN];
    int nel = 0;

    /* Get name of original input file from parameter list */
    p = cpl_parameterlist_find(drvpar->parlist, "filename");
    strncpy(filename, cpl_parameter_get_string(p), MF_MAXLEN);

    /* Get name of converted file */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);
    sprintf(outfile, "%s%s.fits", outdir, outname);

    /* Check possible existence of output file */

    status = mf_basic_access(outfile, F_OK);

    if (status == CPL_ERROR_NONE) {

        /* Open results file if it exists and compare data file names */

        sprintf(resfile, "%s%s_fit.res", outdir, outname);

        if ((stream = fopen(resfile, "r")) == NULL) {

            /* Set error state back to no errors */
            cpl_errorstate_set(CPL_ERROR_NONE);

        } else {

            /* Get name of original input file from results file */
            mf_basic_readline(stream, par, &nel);
            mf_basic_readline(stream, par, &nel);

            /* Close results file */
            fclose(stream);

            /* If file names agree, do not perform conversion */
            if (strcmp(par[0].c, filename) == 0 && nel == 1) {
                cpl_msg_info(cpl_func, "Input data file: %s", filename);
                cpl_msg_info(cpl_func, "Converted input file %s already "
                             "exists", outfile);
                return cpl_error_get_code();
            }

        }

    } else if ((int) status == (int) MF_ERROR_NOENT) {

        /* Set error state back to no errors */
        cpl_errorstate_set(CPL_ERROR_NONE);

    } else {

        /* Error -> return */
        return status;

    }

    /* Read file and convert it to CPL table and CPL property list */
    if ((status = mf_conv_readfile(&tabdat, drvpar)) !=
        CPL_ERROR_NONE) {
        return status;
    }

    /* Write CPL table and CPL property list to FITS table */
    mf_conv_writetable(&tabdat, drvpar);

    /* Free allocated memory */
    mf_conv_tarr_delete(&tabdat);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code mf_readspec_fits(cpl_table *spec, mfdrv *drvpar)
{
    /*!
     * Reads a FITS file with tabulated spectroscopic data (wavelength, flux,
     * flux error, mask) that was created by ::mf_readspec_preptable and puts
     * the read data in a CPL table consisting of the columns "lambda",
     * "flux", and "weight". The names of the required FITS table columns are
     * provided by the ::mfdrv parameter structure. The presence of flux error
     * is optional. For skipping such a column the name has to be 'NULL'. The
     * original input file could also miss a mask column (also indicated by
     * 'NULL'). However, ::mf_readspec_preptable makes sure that a suitable
     * mask column indicated by the given column name or ::MF_DEFMASKCOL (in
     * the case of 'NULL') + '_I' is present. The input file can have an
     * arbitrary number of extensions with different parts of the spectrum (to
     * be defined in the MOLECFIT driver file). A weight of zero is taken if a
     * wavelength is excluded by the mask.
     *
     * \b INPUT:
     * \param spec    empty CPL table
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    CPL table with wavelengths, fluxes, and weights
     * \param drvpar  :mfdrv with number of chips
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    cpl_parameter *p;
    cpl_table *intab;
    cpl_array *colnames;
    cpl_boolean exerr = CPL_TRUE;
    char outdir[MF_MAXLEN] = "";
    char outname[MF_MAXLEN] = "", filename[MF_MAXLEN] = "";
    char errtxt[MF_MAXLEN], col_lam[MF_LENLINE+2], col_flux[MF_LENLINE+2];
    char col_dflux[MF_LENLINE+2], col_mask[MF_LENLINE+2];
    char col_imask[MF_LENLINE+2], colname[MF_LENLINE+2];
    int coln[4] = {0, 0, 0, 0};
    int nchip = 0, i = 0, ncolmin = 5, ncol = 0, check = 0;
    int nspec = 0, jmin = 0, nrow = 0, j = 0, k = 0, mask = 0;
    double wlgtomicron = 0., lam = 0., flux = 0., dflux = 0.;
    double min_lam = 1e200, max_lam = 0.;

    /* Get file name */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);
    sprintf(filename, "%s%s.fits", outdir, outname);

    /* Write info message */
    cpl_msg_info(cpl_func, "Read %s", filename);

    /* Check file existence and get number of FITS extensions */
    nchip = cpl_fits_count_extensions(filename);
    if (nchip < 0) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Write number of chips (FITS extensions) into parameter list */
    p = cpl_parameterlist_find(drvpar->parlist, "nchip");
    cpl_parameter_set_int(p, nchip);

    /* Get column labels from parameter list */
    p = cpl_parameterlist_find(drvpar->parlist, "col_lam");
    strncpy(col_lam, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_lam, "NULL") == 0) {
        sprintf(col_lam, "%s", MF_DEFLAMCOL);
    }
    p = cpl_parameterlist_find(drvpar->parlist, "col_flux");
    strncpy(col_flux, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_flux, "NULL") == 0) {
        sprintf(col_flux, "%s", MF_DEFFLUXCOL);
    }
    p = cpl_parameterlist_find(drvpar->parlist, "col_dflux");
    strncpy(col_dflux, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_dflux, "NULL") == 0) {
        exerr = CPL_FALSE;
        ncolmin--;
    }
    p = cpl_parameterlist_find(drvpar->parlist, "col_mask");
    strncpy(col_mask, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_mask, "NULL") == 0) {
        sprintf(col_imask, "%s_I", MF_DEFMASKCOL);
        ncolmin--;
    } else {
        sprintf(col_imask, "%s_I", col_mask);
    }

    /* Get wavelength unit conversion factor */
    p = cpl_parameterlist_find(drvpar->parlist, "wlgtomicron");
    wlgtomicron = cpl_parameter_get_double(p);

    /* Read selected extensions of FITS file, compute weights,
       and fill output CPL table */

    for (i = 0; i < nchip; i++) {

        /* Read FITS extension */
        intab = cpl_table_load(filename, i+1, 0);

        /* Get column labels */
        colnames = cpl_table_get_column_names(intab);
        ncol = cpl_array_get_size(colnames);

        /* Check existence of columns */

        for (check = 0, j = 0; j < ncol; j++) {
            strncpy(colname, cpl_array_get_string(colnames, j),
                    MF_LENLINE+2);
            if (strcmp(colname, col_lam) == 0) {
                coln[0] = j;
                check++;
            } else if (strcmp(colname, col_flux) == 0) {
                coln[1] = j;
                check++;
            } else if (strcmp(colname, col_dflux) == 0) {
                coln[2] = j;
                check++;
            } else if (strcmp(colname, col_mask) == 0) {
                check++;
            } else if (strcmp(colname, col_imask) == 0) {
                coln[3] = j;
                check++;
            }
        }

        if (check < ncolmin) {
            cpl_table_set_size(spec, 0);
            cpl_array_delete(colnames);
            cpl_table_delete(intab);
            sprintf(errtxt, "%s: %s (missing FITS table column(s))",
                    MF_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Resize output table */
        jmin = nspec;
        nrow = cpl_table_get_nrow(intab);
        nspec += nrow;
        cpl_table_set_size(spec, nspec);

        /* Transfer wavelength and flux/transmission and compute weights */

        for (j = 0; j < nrow; j++) {

            k = j + jmin;

            cpl_table_set(spec, "chip", k, i+1);

            lam = cpl_table_get(intab,
                                cpl_array_get_string(colnames, coln[0]), j,
                                NULL) * wlgtomicron;
            min_lam = MF_MIN(lam, min_lam);
            max_lam = MF_MAX(lam, max_lam);
            cpl_table_set(spec, "lambda", k, lam);

            flux = cpl_table_get(intab,
                                 cpl_array_get_string(colnames, coln[1]), j,
                                 NULL);
            cpl_table_set(spec, "flux", k, flux);

            if (exerr == CPL_TRUE) {
                dflux = cpl_table_get(intab,
                                      cpl_array_get_string(colnames, coln[2]),
                                      j, NULL);
            } else {
                dflux = 1.;
            }

            mask = cpl_table_get(intab,
                                 cpl_array_get_string(colnames, coln[3]), j,
                                 NULL);

            if (dflux <= 0. || mask == 0) {
                cpl_table_set(spec, "weight", k, 0.);
            } else {
                cpl_table_set(spec, "weight", k, 1. / dflux);
            }

        }

        /* Delete temporary CPL objects */
        cpl_array_delete(colnames);
        cpl_table_delete(intab);

    }

    /* Write number of pixels (sum of all chips) into parameter list */
    p = cpl_parameterlist_find(drvpar->parlist, "npix");
    cpl_parameter_set_int(p, nspec);

    /* Check wavelength range, warn only (PIPE-5418) */
    if (min_lam < MF_WAVELENGTH_MIN || max_lam > MF_WAVELENGTH_MAX) {
        cpl_msg_warning(cpl_func,
        "Wavelength range is %g to %g micron (with wlgtomicron=%g). "
        "Is this really intended?", min_lam, max_lam, wlgtomicron);
    }
    else if (min_lam <= 0) {
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS,
             "%s: %s (negative wavelength range, wlgtomicron=%g)",
             MF_ERROR_UFS_TXT, filename, wlgtomicron);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_header(mfdrv *drvpar)
{
    /*!
     * Reads keywords from a header of a FITS table that was created by
     * ::mf_readspec_preptable and puts the read values in the ::mfdrv
     * parameter structure. Keywords related to the telescope site and
     * observing conditions are mainly considered.
     *
     * The first required header keyword found defines the FITS extension
     * where all keywords are read from. If an input file does not include any
     * required keywords, all parameters will be read from the ::mfdrv
     * parameter structure. If a keyword is not set, the default value is
     * taken.
     *
     * A parameter is not read from the file if the corresponding keyword name
     * parameter ("_key") is set to "NONE". If an expected keyword (not
     * "NONE") cannot be found in a valid FITS header, an error message will
     * be written to the CPL error structure. An exception is the slit width
     * with the default keyword "ESO INS SLIT1 WID", whose presence is
     * instrument dependent. CRIRES and VISIR data provide this keyword. For
     * X-Shooter frames (which do not contain this keyword), the slit width is
     * obtained in a different way by a special subroutine. For all other
     * cases, the original value of the ::mfdrv parameter structure remains
     * untouched. This will be indicated by an info message.
     *
     * Finally, the \e obsdate parameter is converted into years if the date
     * is provided as MJD.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param drvpar  ::mfdrv parameter structure with FITS header values
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     * - Invalid object value(s)
     */

    cpl_parameter *p;
    cpl_propertylist *header = NULL;
    cpl_property *prop = NULL;
    char outdir[MF_MAXLEN] = "";
    char outname[MF_MAXLEN] = "", filename[MF_MAXLEN] = "";
    char errtxt[MF_MAXLEN], key[MF_MAXLEN] = "", keystr[MF_MAXLEN] = "";
    cpl_boolean isheader = CPL_FALSE;
    int next = 0, ext = 0, i = 0, nerr = 0;
    double val = 0.;

    /* Names of parameters */
    char par[MF_NKEY][MF_LENLINE+2] = {"obsdate", "utc", "telalt", "rhum",
                                       "pres", "temp", "m1temp", "geoelev",
                                       "longitude", "latitude", "slitw",
                                       "pixsc"};

    /* Get file name */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);
    sprintf(filename, "%s%s.fits", outdir, outname);

    /* Find FITS extension with required keywords (first detection decides)
       and write the header of this extension into CPL property list */
    next = cpl_fits_count_extensions(filename);
    for (ext = 0; ext <= next; ext++) {
        header = cpl_propertylist_load(filename, ext);
        if (header == NULL) {
            sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s",
                                         errtxt);
        }
        /* Check existence of keywords */
        for (i = 0; i < MF_NKEY; i++) {
            sprintf(key, "%s_key", par[i]);
            p = cpl_parameterlist_find(drvpar->parlist, key);
            strncpy(keystr, cpl_parameter_get_string(p), MF_MAXLEN);
            prop = cpl_propertylist_get_property(header, keystr);
            if (prop != NULL) {
                break;
            }
        }
        if (prop != NULL) {
            break;
        } else {
            cpl_propertylist_delete(header);
        }
    }

    /* Set header flag and Write info message */
    if (prop != NULL) {
        isheader = CPL_TRUE;
        cpl_msg_info(cpl_func, "Take keywords from FITS extension %d", ext);
    } else {
        isheader = CPL_FALSE;
        cpl_msg_info(cpl_func, "Header keywords are taken from parameter "
                     "file");
    }

    /* Update mfdrv parameter structure with keywords from
       CPL property list */

    for (i = 0; i < MF_NKEY; i++) {

        /* Read property */
        if (isheader == CPL_TRUE) {
            sprintf(key, "%s_key", par[i]);
            p = cpl_parameterlist_find(drvpar->parlist, key);
            strncpy(keystr, cpl_parameter_get_string(p), MF_MAXLEN);
            prop = cpl_propertylist_get_property(header, keystr);
        }

        if (prop == NULL && strcmp(keystr, "NONE") != 0 &&
            isheader == CPL_TRUE) {

            if (i == MF_NKEY - 2) {
                /* Get slit width if XSHOOTER data is provided */
                mf_readspec_slitwidth_xshooter(&val, header);
                if (val > 0) {
                    p = cpl_parameterlist_find(drvpar->parlist, par[i]);
                    cpl_parameter_set_double(p, val);
                }  else {
                    cpl_msg_info(cpl_func, "%s is taken from parameter file",
                                 par[i]);
                }
            } else {
                /* Set error message in the case of missing
                   telescope-independent keyword */
                nerr++;
                sprintf(errtxt, "%s: %s (keyword %s not found)",
                        MF_ERROR_UFS_TXT, filename, keystr);
                cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
            }

        } else if (prop == NULL &&
                   (strcmp(keystr, "NONE") == 0 || isheader == CPL_FALSE)) {

            /* Write info message */
            cpl_msg_info(cpl_func, "%s is taken from parameter file", par[i]);

        } else {

            /* Write info message */
            cpl_msg_info(cpl_func, "Get %s from keyword %s", par[i], keystr);

            /* Set parameter value */
            p = cpl_parameterlist_find(drvpar->parlist, par[i]);
            cpl_parameter_set_double(p, cpl_property_get_double(prop));

        }

    }

    /* Free allocated memory */
    if (isheader == CPL_TRUE) {
        cpl_propertylist_delete(header);
    }

    /* Return MF_ERROR_UFS in the case of keyword mismatch */
    if (nerr > 0) {
        return MF_ERROR_UFS;
    }

    /* Check setting of obsdate, which is -1 by default, and convert MJD in
       years if required */
    p = cpl_parameterlist_find(drvpar->parlist, par[0]);
    val = cpl_parameter_get_double(p);
    if (val < 0) {
        sprintf(errtxt, "%s: obsdate of mfdrv *drvpar < 0 (probably not set)",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    } else if (val >= 10000) {
        /* MJD -> date in years */
        cpl_msg_info(cpl_func, "Convert MJD into date in years");
        val = mf_basic_mjd2fracyear(val);
        cpl_parameter_set_double(p, val);
    }

    /* Check setting of utc, which is -1 by default */
    p = cpl_parameterlist_find(drvpar->parlist, par[1]);
    val = cpl_parameter_get_double(p);
    if (val < 0) {
        sprintf(errtxt, "%s: utc of mfdrv *drvpar < 0 (probably not set)",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_slitwidth_xshooter(double *slitw,
                                              const cpl_propertylist *header)
{
    /*!
     * Obtains slit width of XSHOOTER exposures from FITS header. Returns -1
     * if the header of another instrument is provided.
     *
     * \b INPUT:
     * \param header  CPL property list with FITS header
     *
     * \b OUTPUT:
     * \param slitw   slit width in arcsec
     *
     * \b ERRORS:
     * - none
     */

    const cpl_property *prop;
    char slitkey[MF_LENLINE+2] = "", slitstr[4] = "   ";
    int opti = 0;

    /* Default slit width */
    *slitw = -1;

    /* Check instrument */
    prop = cpl_propertylist_get_property_const(header, "INSTRUME");
    if (strstr(cpl_property_get_string(prop), "SHOOT") == NULL) {
        /* Instrument is not XSHOOTER */
        return CPL_ERROR_NONE;
    }

    /* Check XSHOOTER arm (UVB, VIS, or NIR) */
    prop = cpl_propertylist_get_property_const(header, "ESO SEQ ARM");
    if (strstr(cpl_property_get_string(prop), "UVB") != NULL) {
        opti = 3;
    } else if (strstr(cpl_property_get_string(prop), "VIS") != NULL) {
        opti = 4;
    } else if (strstr(cpl_property_get_string(prop), "NIR") != NULL) {
        opti = 5;
    } else {
        /* Unexpected XSHOOTER arm */
        return CPL_ERROR_NONE;
    }

    /* Get name of fits keyword and write info message */
    sprintf(slitkey, "ESO INS OPTI%d NAME", opti);
    cpl_msg_info(cpl_func, "Get slitw from %s", slitkey);

    /* Get slit width for identified arm */
    prop = cpl_propertylist_get_property_const(header, slitkey);
    strncpy(slitstr, cpl_property_get_string(prop), 3);
    mf_basic_terminatestring(slitstr);
    *slitw = strtod(slitstr, NULL);
    if (*slitw == 0) {
        *slitw = -1;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_airtovac(cpl_table *spec, const mfdrv *drvpar)
{
    /*!
     * Converts air wavelengths to vacuum wavelengths by using the formula of
     * Edlen (1966). No action is performed if the wavelengths of the input
     * spectrum are already for vacuum.
     *
     * \b INPUT:
     * \param spec    spectrum with air wavelengths in \f$\mu{\rm m}\f$
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    spectrum with vacuum wavelengths in \f$\mu{\rm m}\f$
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    char vac_air[4];
    int nrow = 0, i = 0;
    double *lam, sig2 = 0., n = 0.;

    /* Input spectrum: vacuum or air wavelengths? */
    p = cpl_parameterlist_find(drvpar->parlist, "vac_air");
    strcpy(vac_air, cpl_parameter_get_string(p));

    /* No conversion in the case of input vacuum wavelengths */
    if (strncmp(vac_air, "air", 3) != 0) {
        if (strncmp(vac_air, "vac", 3) != 0) {
            cpl_msg_warning(cpl_func, "mfdrv *drvpar: vac_air neither 'vac' "
                            "nor 'air' -> no air to vacuum conversion");
        }
        return CPL_ERROR_NONE;
    }

    /* Get number of rows */
    nrow = cpl_table_get_nrow(spec);

    /* Get pointer to wavelengths */
    lam = cpl_table_get_data_double(spec, "lambda");

    /* Get refractive index and convert wavelengths */
    for (i = 0; i < nrow; i++) {
         sig2 = pow(lam[i], -2);
         n = 8342.13 + 2406030. / (130. - sig2) + 15997. / (38.9 - sig2);
         n = 1. + 1e-8 * n;
         lam[i] *= n;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_fillpartab(mfdrv *drvpar, cpl_table *spec)
{
    /*!
     * Fills the parameter tables of the ::mfdrv parameter structure for the
     * number of chips determined by the file read routine (see
     * ::mf_readspec_fits).
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     * \param spec    CPL table with observed spectrum
     *
     * \b OUTPUT:
     * \param drvpar  ::mfdrv with parameter tables adapted to number of chips
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    cpl_array *cont_coef = NULL, *wlc_coef = NULL;
    cpl_table *chipspec = NULL;
    int nchip = 0, nwlc = 0, ncont = 0, fit_cont = 0, i = 0, fit_wlc = 0;
    int nsel = 0;
    double cont0 = 0., wlc0 = 0.;

    /* Get number of chips and polynomial coefficients */
    p = cpl_parameterlist_find(drvpar->parlist, "nchip");
    nchip = cpl_parameter_get_int(p);
    p = cpl_parameterlist_find(drvpar->parlist, "wlc_n");
    nwlc = cpl_parameter_get_int(p) + 1;
    p = cpl_parameterlist_find(drvpar->parlist, "cont_n");
    ncont = cpl_parameter_get_int(p) + 1;

    /* Set range number to chip number */
    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    cpl_parameter_set_int(p, nchip);

    /* Get fit flag for continuum fit */
    p = cpl_parameterlist_find(drvpar->parlist, "fit_cont");
    fit_cont = cpl_parameter_get_int(p);

    /* Get constant term for continuum fit */
    p = cpl_parameterlist_find(drvpar->parlist, "cont_const");
    cont0 = cpl_parameter_get_double(p);

    /* Resize table for range-related parameters */
    cpl_table_set_size(drvpar->rangetab, nchip);
    cpl_table_set_column_depth(drvpar->rangetab, "cont_coef", ncont);

    /* Fill table for range-related parameters */

    /* Coefficients for continuum correction */
    cont_coef = cpl_array_new(ncont, CPL_TYPE_DOUBLE);
    cpl_array_fill_window(cont_coef, 0, ncont, 0.);
    cpl_array_set(cont_coef, 0, cont0);
    for (i = 0; i < nchip; i++) {
        cpl_table_set(drvpar->rangetab, "chip", i, i + 1);
        cpl_table_set(drvpar->rangetab, "fit_range", i, fit_cont);
        cpl_table_set_array(drvpar->rangetab, "cont_coef", i, cont_coef);
    }
    cpl_array_delete(cont_coef);

    /* Get fit flag for wavelength correction */
    p = cpl_parameterlist_find(drvpar->parlist, "fit_wlc");
    fit_wlc = cpl_parameter_get_int(p);

    /* Get constant term for wavelength correction */
    p = cpl_parameterlist_find(drvpar->parlist, "wlc_const");
    wlc0 = cpl_parameter_get_double(p);

    /* Resize table for chip-related parameters */
    cpl_table_set_size(drvpar->chiptab, nchip);
    cpl_table_set_column_depth(drvpar->chiptab, "wlc_coef", nwlc);

    /* Fill table for chip-related parameters */

    /* Coefficients for wavelength correction */
    wlc_coef = cpl_array_new(nwlc, CPL_TYPE_DOUBLE);
    cpl_array_fill_window(wlc_coef, 0, nwlc, 0.);
    cpl_array_set(wlc_coef, 0, wlc0);
    cpl_array_set(wlc_coef, 1, 1.);
    cpl_table_fill_column_window(drvpar->chiptab, "fit_chip", 0, nchip,
                                 fit_wlc);
    cpl_table_fill_column_window_array(drvpar->chiptab, "wlc_coef", 0, nchip,
                                       wlc_coef);
    cpl_array_delete(wlc_coef);

    /* Minimum and maximum wavelengths */
    for (i = 0; i < nchip; i++) {
        /* Extract wavelength grid of observed spectrum for each chip */
        cpl_table_unselect_all(spec);
        cpl_table_or_selected_int(spec, "chip", CPL_EQUAL_TO, i+1);
        chipspec = cpl_table_extract_selected(spec);
        cpl_table_select_all(spec);
        nsel = cpl_table_get_nrow(chipspec);
        /* Get minimum and maximum wavelengths and write them into table */
        cpl_table_set(drvpar->chiptab, "wl_min", i,
                      cpl_table_get(chipspec, "lambda", 0, NULL));
        cpl_table_set(drvpar->chiptab, "wl_max", i,
                      cpl_table_get(chipspec, "lambda", nsel-1, NULL));
        /* Delete temporary table */
        cpl_table_delete(chipspec);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_cutspec(cpl_table *spec, const mfdrv *drvpar)
{
    /*!
     * Cuts read spectrum at wavelengths in \f$\mu{\rm m}\f$ provided by the
     * ::mfdrv parameter structure.
     *
     * \b INPUT:
     * \param spec    CPL table with observed spectrum
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    CPL table with cut spectrum
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    cpl_parameter *p;
    char errtxt[MF_MAXLEN];
    int m = 0;
    double cut_min = 0., cut_max = 0.;

    /* Read cut limits */
    p = cpl_parameterlist_find(drvpar->parlist, "cut_min");
    cut_min = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find(drvpar->parlist, "cut_max");
    cut_max = cpl_parameter_get_double(p);

    /* Check validity of limits */

    if (cut_max <= cut_min) {
        sprintf(errtxt, "%s: cut_max <= cut_min in mfdrv *drvpar",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }
    if (cut_max <= cpl_table_get(spec, "lambda", 0, NULL)) {
        sprintf(errtxt, "%s: cut_max of mfdrv *drvpar <= lambda_min of "
                "cpl_table *spec", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }
    m = cpl_table_get_nrow(spec);
    if (cut_min >= cpl_table_get(spec, "lambda", m-1, NULL)) {
        sprintf(errtxt, "%s: cut_min of mfdrv *drvpar >= lambda_max of "
                "cpl_table *spec", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Erase undesired wavelength range */
    cpl_table_and_selected_double(spec, "lambda", CPL_LESS_THAN, cut_min);
    cpl_table_or_selected_double(spec, "lambda", CPL_GREATER_THAN, cut_max);
    cpl_table_erase_selected(spec);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_includeranges(cpl_table *spec, mfdrv *drvpar)
{
    /*!
     * Reads wavelength ranges (in \f$\mu{\rm m}\f$) from a two-column ASCII
     * or FITS file (lower and upper limits) and selects these ranges in the
     * CPL table of the observed spectrum for the fitting procedure. Since the
     * radiative transfer code is only run for these ranges, the selection is
     * crucial for the code execution time. If a range table is not available,
     * the parameter 'wrange_include' in the ::mfdrv structure has to be set
     * to 'none'. Apart from pairs of wavelength limits, an ASCII file can
     * also contain empty lines or comment lines starting with '#'.
     *
     * \b INPUT:
     * \param spec    CPL table with observed spectrum
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    spectrum with mrange > 0 for selected ranges
     * \param drvpar  ::mfdrv with modified continuum fit parameters
     *
     * \b ERRORS:
     * - Unexpected file structure
     * - File opening failed
     * - see subroutines
     */

    FILE *stream;
    cpl_parameter *p;
    cpl_table *ranges;
    cpl_array *colnames = NULL, *rangechips, *cont_coef;
    cpl_boolean isvalidnlim = CPL_TRUE, isvalidranges = CPL_TRUE;
    char filename[MF_MAXLEN] = "";
    char relfilename[MF_MAXLEN] = "", errtxt[MF_MAXLEN];
    char str[MF_LENLINE+2];
    int nrow = 0, fitsformat = 0, i = 0, chip = 0, nline = 0, nlim = 0;
    int nrange = 0, j = 0, h = 0, nrangeext = 0, range = 0, chip0 = 0;
    int fit_cont = 0, nw = 0;
    double llim = 0., ulim = 0., minlam = 0., maxlam = 0., lam = 0.;

    /* Default mrange column: full spectrum selected and range number for each
       chip */
    nrow = cpl_table_get_nrow(spec);
    for (i = 0; i < nrow; i++) {
        chip = cpl_table_get(spec, "chip", i, NULL);
        cpl_table_set(spec, "mrange", i, chip);
    }

    /* Get name of range table */
    p = cpl_parameterlist_find(drvpar->parlist, "wrange_include");
    strncpy(filename, cpl_parameter_get_string(p), MF_MAXLEN);
    if (strncmp(filename, "none", 4) == 0) {
        /* Return if name of range table is "none" */
        cpl_msg_info(cpl_func, "Fit full spectrum (no fit ranges)");
        return CPL_ERROR_NONE;
    }
    if (filename[0] != '/') {
        char curdir[MF_MAXLEN];
        p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        sprintf(relfilename, "%s/%s", curdir, filename);
        mf_basic_absfile(filename, relfilename);
    }

    /* Write info message */
    cpl_msg_info(cpl_func, "Read fit ranges from %s", filename);

    /* Get file type of range table */
    mf_conv_checkfitsformat(&fitsformat, filename);
    if (fitsformat == -1) {
        return cpl_error_get_code();
    } else if (fitsformat >= 2) {
        sprintf(errtxt, "%s: %s (unsupported FITS image)", MF_ERROR_UFS_TXT,
                filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Transfer ASCII or FITS data to temporary CPL table */

    if (fitsformat == 0) {

        /* Check file existence */
        if ((stream = fopen(filename, "r")) == NULL) {
            sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s",
                                         errtxt);
        }

        /* Find number of ranges and check number of values per line */
        while (fgets(str, MF_LENLINE+2, stream) != NULL) {
            nline++;
            if (str[0] == '#' || strspn(str, "\n\t ") == strlen(str)) {
                continue;
            }
            nlim = sscanf(str, "%le %le", &llim, &ulim);
            if (nlim < 2) {
                isvalidnlim = CPL_FALSE;
            }
            nrange++;
        }
        rewind(stream);

        /* Return if no ranges are provided */
        if (nrange == 0) {
            fclose(stream);
            return CPL_ERROR_NONE;
        }

        /* Exit in the case of missing range limits */
        if (isvalidnlim == CPL_FALSE) {
            fclose(stream);
            sprintf(errtxt, "%s: %s (not two values per line)",
                    MF_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Initialise CPL table */
        ranges = cpl_table_new(nrange);
        cpl_table_new_column(ranges, "llim", CPL_TYPE_DOUBLE);
        cpl_table_new_column(ranges, "ulim", CPL_TYPE_DOUBLE);

        /* Read range table and write data into CPL table */
        for (j = 0, h = 0; h < nline; h++) {
            if (fgets(str, MF_LENLINE+2, stream)) {};
            if (str[0] == '#' || strspn(str, "\n\t ") == strlen(str)) {
                continue;
            }
            sscanf(str, "%le %le", &llim, &ulim);
            cpl_table_set(ranges, "llim", j, llim);
            cpl_table_set(ranges, "ulim", j, ulim);
            j++;
        }

        fclose(stream);

    } else {

        /* Load FITS file (extension 1!) into CPL table */
        ranges = cpl_table_load(filename, 1, 0);
        if (ranges == NULL) {
            sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s",
                                         errtxt);
        }

        /* Return if no ranges are provided */
        nrange = cpl_table_get_nrow(ranges);
        if (nrange == 0) {
            cpl_table_delete(ranges);
            return CPL_ERROR_NONE;
        }

        /* Exit in the case of a missing range column  */
        nlim = cpl_table_get_ncol(ranges);
        if (nlim != 2) {
            cpl_table_delete(ranges);
            sprintf(errtxt, "%s: %s (number of columns != 2)",
                    MF_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Rename columns if necessary */
        colnames = cpl_table_get_column_names(ranges);
        if (cpl_table_has_column(ranges, "llim") == 0) {
            cpl_table_name_column(ranges, cpl_array_get_string(colnames, 0),
                                  "llim");
        }
        if (cpl_table_has_column(ranges, "ulim") == 0) {
            cpl_table_name_column(ranges, cpl_array_get_string(colnames, 1),
                                  "ulim");
        }

    }

    /* Remove ranges outside valid wavelength range */
    minlam = cpl_table_get_column_min(spec, "lambda");
    maxlam = cpl_table_get_column_max(spec, "lambda");
    cpl_table_unselect_all(ranges);
    cpl_table_or_selected_double(ranges, "llim", CPL_GREATER_THAN, maxlam);
    cpl_table_or_selected_double(ranges, "ulim", CPL_LESS_THAN, minlam);
    cpl_table_erase_selected(ranges);

    /* Check number of remaining ranges */
    nrange = cpl_table_get_nrow(ranges);
    if (nrange == 0) {
        cpl_table_delete(ranges);
        sprintf(errtxt, "%s: %s (no range inside spectrum cpl_table *spec)",
                MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Create array for chip numbers of selected ranges */
    rangechips = cpl_array_new(nrange, CPL_TYPE_INT);

    /* Initialise column for range flag in spectrum table */
    cpl_table_fill_column_window_int(spec, "mrange", 0, nrow, 0);

    /* Set flags for ranges listed in the temporary CPL table */

    for (nrangeext = nrange, range = 0, j = 0; j < nrange; j++) {
        llim = cpl_table_get(ranges, "llim", j, NULL);
        ulim = cpl_table_get(ranges, "ulim", j, NULL);
        if (ulim <= llim) {
            isvalidranges = CPL_FALSE;
            nrangeext--;
            continue;
        }
        range++;
        for (chip0 = -1, i = 0; i < nrow; i++) {
            chip = cpl_table_get(spec, "chip", i, NULL);
            lam = cpl_table_get(spec, "lambda", i, NULL);
            if (lam >= llim && lam <= ulim) {
                if (chip0 == -1) {
                    chip0 = cpl_table_get(spec, "chip", i, NULL);
                }
                if (chip0 != chip) {
                    /* New range number if chip changes */
                    range++;
                    nrangeext++;
                    chip0 = chip;
                    cpl_array_set_size(rangechips, nrangeext);
                }
                cpl_table_set(spec, "mrange", i, range);
                cpl_array_set(rangechips, range-1, chip0);
            }
        }
    }

    /* Set number of ranges in parameter list */
    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    cpl_parameter_set_int(p, nrangeext);

    /* Get fit flag for continuum fit */
    p = cpl_parameterlist_find(drvpar->parlist, "fit_cont");
    fit_cont = cpl_parameter_get_int(p);

    /* Expand table of continuum fit parameters in mfdrv structure */
    cpl_table_set_size(drvpar->rangetab, nrangeext);
    cont_coef = cpl_array_duplicate(cpl_table_get_array(drvpar->rangetab,
                                                        "cont_coef", 0));
    for (j = 0; j < nrangeext; j++) {
        cpl_table_set(drvpar->rangetab, "chip", j,
                      cpl_array_get(rangechips, j, NULL));
        cpl_table_set(drvpar->rangetab, "fit_range", j, fit_cont);
        cpl_table_set_array(drvpar->rangetab, "cont_coef", j, cont_coef);
    }

    /* Free allocated memory */
    if (fitsformat != 0) {
        cpl_array_delete(colnames);
    }
    cpl_array_delete(rangechips);
    cpl_array_delete(cont_coef);
    cpl_table_delete(ranges);

    /* Treat errors */

    if (isvalidranges == CPL_FALSE) {
        sprintf(errtxt, "%s: %s (upper limit <= lower limit)",
                MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    cpl_table_and_selected_double(spec, "weight", CPL_GREATER_THAN, 0.);
    nw = cpl_table_count_selected(spec);
    cpl_table_select_all(spec);
    if (nw == 0) {
        sprintf(errtxt, "%s: cpl_table *spec (all weights = 0)",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    cpl_table_and_selected_int(spec, "mrange", CPL_GREATER_THAN, 0);
    nw = cpl_table_count_selected(spec);
    cpl_table_select_all(spec);
    if (nrange > 0 && nw == 0) {
        sprintf(errtxt, "%s: cpl_table *spec (wavelength inclusion range "
                "selects no values in the input data)", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_replacecoef(mfdrv *drvpar, cpl_table *rangetab,
                                       cpl_table *chiptab)
{
    /*!
     * Replaces entries in the range and chip tables of the ::mfdrv parameter
     * structure by values provided by the parameter file. In the range table,
     * the range-specific fit flags and the coefficients of the continuum fit
     * are updated. In the chip table, the chip-specific fit flags and the
     * coefficients of the wavelength correction are modified. Only those
     * table entries are changed which were explicitly given in the parameter
     * file using the keywords fit_range[range], cont_range[range],
     * fit_chip[chip], and wlc_chip[chip].
     *
     * \b INPUT:
     * \param drvpar    ::mfdrv parameter structure
     * \param rangetab  range table with data from the input parameter file
     * \param chiptab   chip table with data from the input parameter file
     *
     * \b OUTPUT:
     * \param drvpar    ::mfdrv with modified range and chip tables
     *
     * \b ERRORS:
     * - none
     */

    cpl_array *coef_in = NULL, *coef_out = NULL;
    int nrange = 0, ncont = 0, i = 0, fit = 0, j = 0, nchip = 0, nwlc = 0;
    double coef = 0.;

    /* Adapt size of input range table with data to be added */
    nrange = cpl_table_get_nrow(drvpar->rangetab);
    cpl_table_set_size(rangetab, nrange);

    /* Adapt maximum number of continuum coefficients in input range table
       with data to be added */
    ncont = cpl_table_get_column_depth(drvpar->rangetab, "cont_coef");
    cpl_table_set_column_depth(rangetab, "cont_coef", ncont);

    /* Replace range-specific flags for continuum fit if they were set in the
       parameter file */
    for (i = 0; i < nrange; i++) {
        if (cpl_table_is_valid(rangetab, "fit_range", i) == 1) {
            fit = cpl_table_get(rangetab, "fit_range", i, NULL);
            cpl_table_set(drvpar->rangetab, "fit_range", i, fit);
        }
    }

    /* Replace range-specific continuum fit coefficients if they were set in
       the parameter file */
    for (i = 0; i < nrange; i++) {
        if (!cpl_table_is_valid(rangetab, "cont_coef", i)) {
            continue;
        }
        coef_in = cpl_array_duplicate(cpl_table_get_array(rangetab,
                                                          "cont_coef", i));
        coef_out = cpl_array_duplicate(cpl_table_get_array(drvpar->rangetab,
                                                           "cont_coef", i));
        for (j = 0; j < ncont; j++) {
            if (cpl_array_is_valid(coef_in, j) == 1) {
                coef = cpl_array_get(coef_in, j, NULL);
                cpl_array_set(coef_out, j, coef);
            }
        }
        cpl_table_set_array(drvpar->rangetab, "cont_coef", i, coef_out);
        cpl_array_delete(coef_in);
        cpl_array_delete(coef_out);
    }

    /* Adapt size of input chip table with data to be added */
    nchip = cpl_table_get_nrow(drvpar->chiptab);
    cpl_table_set_size(chiptab, nchip);

    /* Adapt maximum number of wavelength coefficients in input chip table
       with data to be added */
    nwlc = cpl_table_get_column_depth(drvpar->chiptab, "wlc_coef");
    cpl_table_set_column_depth(chiptab, "wlc_coef", nwlc);

    /* Replace chip-specific flags for wavelength fit if they were set in the
       parameter file */
    for (i = 0; i < nchip; i++) {
        if (cpl_table_is_valid(chiptab, "fit_chip", i) == 1) {
            fit = cpl_table_get(chiptab, "fit_chip", i, NULL);
            cpl_table_set(drvpar->chiptab, "fit_chip", i, fit);
        }
    }

    /* Replace chip-specific wavelength fit coefficients if they were set in
       the parameter file */
    for (i = 0; i < cpl_table_get_nrow(chiptab); i++) {
        if (!cpl_table_is_valid(chiptab, "wlc_coef", i)) {
            continue;
        }
        coef_in = cpl_array_duplicate(cpl_table_get_array(chiptab,
                                                          "wlc_coef", i));
        coef_out = cpl_array_duplicate(cpl_table_get_array(drvpar->chiptab,
                                                           "wlc_coef", i));
        for (j = 0; j < nwlc; j++) {
            if (cpl_array_is_valid(coef_in, j) == 1) {
                coef = cpl_array_get(coef_in, j, NULL);
                cpl_array_set(coef_out, j, coef);
            }
        }
        cpl_table_set_array(drvpar->chiptab, "wlc_coef", i, coef_out);
        cpl_array_delete(coef_in);
        cpl_array_delete(coef_out);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_excluderanges(cpl_table *spec,
                                         const mfdrv *drvpar, const char wp)
{
    /*!
     * Reads wavelength (in \f$\mu{\rm m}\f$) or pixel ranges (indicated by
     * 'w' or 'p' for input parameter 'wp') from a two-column ASCII or FITS
     * file (lower and upper limits) and sets the weight of these ranges in
     * the CPL table of the observed spectrum to zero. If a range table is not
     * available, the parameter 'wrange_exclude' or 'prange_exclude' in the
     * ::mfdrv structure has to be set to 'none'. Apart from pairs of
     * wavelength or pixel limits, an ASCII file can also contain empty lines
     * or comment lines starting with '#'.
     *
     * \b INPUT:
     * \param spec    CPL table with observed spectrum
     * \param drvpar  ::mfdrv parameter structure
     * \param wp      w[avelength] or p[ixel] ranges
     *
     * \b OUTPUT:
     * \param spec    spectrum with weight = 0 for excluded ranges
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     * - Unexpected file structure
     * - File opening failed
     * - see subroutines
     */

    FILE *stream;
    cpl_parameter *p = NULL;
    cpl_table *ranges;
    cpl_array *colnames = NULL;
    cpl_boolean isvalidnlim = CPL_TRUE, isvalidranges = CPL_TRUE;
    char filename[MF_MAXLEN] = "";
    char relfilename[MF_MAXLEN] = "", errtxt[MF_MAXLEN];
    char str[MF_LENLINE+2];
    int fitsformat = 0, nlim = 0, nline = 0, nrange = 0, j = 0, h = 0;
    int nrow = 0, i = 0, nw = 0;
    double llim = 0, ulim = 0, val = 0;

    /* Check letter for wavelengths or pixels */
    if (wp != 'w' && wp != 'p') {
        sprintf(errtxt, "%s: char wp ('w' or 'p' only)", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Get name of range table (wavelengths or pixels) */
    if (wp == 'w') {
        p = cpl_parameterlist_find(drvpar->parlist, "wrange_exclude");
    } else if (wp == 'p') {
        p = cpl_parameterlist_find(drvpar->parlist, "prange_exclude");
    }
    strncpy(filename, cpl_parameter_get_string(p), MF_MAXLEN);
    if (strncmp(filename, "none", 1) == 0) {
        /* Return if name of range table is "none" */
        return CPL_ERROR_NONE;
    }
    if (filename[0] != '/') {
        char curdir[MF_MAXLEN];
        p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        sprintf(relfilename, "%s/%s", curdir, filename);
        mf_basic_absfile(filename, relfilename);
    }

    /* Write info message */
    if (wp == 'w') {
        cpl_msg_info(cpl_func, "Exclude wavelength ranges provided by %s",
                     filename);
    } else if (wp == 'p') {
        cpl_msg_info(cpl_func, "Exclude pixel ranges provided by %s",
                     filename);
    }

    /* Get file type of range table */
    mf_conv_checkfitsformat(&fitsformat, filename);
    if (fitsformat == -1) {
        return cpl_error_get_code();
    } else if (fitsformat >= 2) {
        sprintf(errtxt, "%s: %s (unsupported FITS image)", MF_ERROR_UFS_TXT,
                filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Transfer ASCII or FITS data to temporary CPL table */

    if (fitsformat == 0) {

        /* Check file existence */
        if ((stream = fopen(filename, "r")) == NULL) {
            sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s",
                                         errtxt);
        }

        /* Find number of ranges and check number of values per line */
        while (fgets(str, MF_LENLINE+2, stream) != NULL) {
            nline++;
            if (str[0] == '#' || strspn(str, "\n\t ") == strlen(str)) {
                continue;
            }
            nlim = sscanf(str, "%lf %lf", &llim, &ulim);
            if (nlim < 2) {
                isvalidnlim = CPL_FALSE;
            }
            nrange++;
        }
        rewind(stream);

        /* Return if no ranges are provided */
        if (nrange == 0) {
            fclose(stream);
            return CPL_ERROR_NONE;
        }

        /* Exit in the case of missing range limits */
        if (isvalidnlim == CPL_FALSE) {
            fclose(stream);
            sprintf(errtxt, "%s: %s (not two values per line)",
                    MF_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Initialise CPL table */
        ranges = cpl_table_new(nrange);
        cpl_table_new_column(ranges, "llim", CPL_TYPE_DOUBLE);
        cpl_table_new_column(ranges, "ulim", CPL_TYPE_DOUBLE);

        /* Read range table and write data into CPL table */
        for (j = 0, h = 0; h < nline; h++) {
            if (fgets(str, MF_LENLINE+2, stream)) {};
            if (str[0] == '#' || strspn(str, "\n\t ") == strlen(str)) {
                continue;
            }
            sscanf(str, "%lf %lf", &llim, &ulim);
            cpl_table_set(ranges, "llim", j, llim);
            cpl_table_set(ranges, "ulim", j, ulim);
            j++;
        }

        fclose(stream);

    } else {

        /* Load FITS file (extension 1!) into CPL table */
        ranges = cpl_table_load(filename, 1, 0);
        if (ranges == NULL) {
            sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s",
                                         errtxt);
        }

        /* Return if no ranges are provided */
        nrange = cpl_table_get_nrow(ranges);
        if (nrange == 0) {
            cpl_table_delete(ranges);
            return CPL_ERROR_NONE;
        }

        /* Exit in the case of a missing range column  */
        nlim = cpl_table_get_ncol(ranges);
        if (nlim != 2) {
            cpl_table_delete(ranges);
            sprintf(errtxt, "%s: %s (number of columns != 2)",
                    MF_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Rename columns if necessary */
        colnames = cpl_table_get_column_names(ranges);
        if (cpl_table_has_column(ranges, "llim") == 0) {
            cpl_table_name_column(ranges, cpl_array_get_string(colnames, 0),
                                  "llim");
        }
        if (cpl_table_has_column(ranges, "ulim") == 0) {
            cpl_table_name_column(ranges, cpl_array_get_string(colnames, 1),
                                  "ulim");
        }

    }

    /* Set weight to zero for ranges listed in the temporary CPL table */
    nrow = cpl_table_get_nrow(spec);
    for (j = 0; j < nrange; j++) {
        llim = cpl_table_get(ranges, "llim", j, NULL);
        ulim = cpl_table_get(ranges, "ulim", j, NULL);
        if (ulim < llim) {
            isvalidranges = CPL_FALSE;
            continue;
        }
        for (i = 0; i < nrow; i++) {
            if (wp == 'w') {
                /* wavelength */
                val = cpl_table_get(spec, "lambda", i, NULL);
            } else if (wp == 'p') {
                /* pixel */
                val = i + 1;
            }
            if (val >= llim && val <= ulim) {
                cpl_table_set(spec, "weight", i, 0.);
            }
        }
    }

    /* Free allocated memory */
    if (fitsformat != 0) {
        cpl_array_delete(colnames);
    }
    cpl_table_delete(ranges);

    /* Treat errors */

    if (isvalidranges == CPL_FALSE) {
        sprintf(errtxt, "%s: %s (upper limit < lower limit)",
                MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    cpl_table_and_selected_double(spec, "weight", CPL_GREATER_THAN, 0.);
    nw = cpl_table_count_selected(spec);
    cpl_table_select_all(spec);
    if (nw == 0) {
        sprintf(errtxt, "%s: cpl_table *spec (exclusion range "
                "rejects all values in the input data)", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_modranges(mfdrv *drvpar, const cpl_table *spec)
{
    /*!
     * Checks whether the selected ranges can be covered by a single model
     * spectrum. If this is the case, the ::mfdrv parameter structure will be
     * modified. As criterion, the largest step in wavenumber between adjacent
     * range pixels is compared to two times ::MF_EXTRACOVER. The purpose of
     * this routine is an optimisation of the code run time.
     *
     * \b INPUT:
     * \param spec    CPL table with observed spectrum
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param drvpar  ::mfdrv with modified range table and singlespec
     *                parameter
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    int nrow = 0, i = 0, range = 0, nrange = 0, j = 0;
    double own = 0., maxstep = 0., nwn = 0., dwn = 0.;

    /* Find largest step in wavenumber for range pixels */
    nrow = cpl_table_get_nrow(spec);
    for (own = 0., maxstep = 0., i = 0; i < nrow; i++) {
        range = cpl_table_get(spec, "mrange", i, NULL);
        if (range > 0) {
            nwn = 1e4 / cpl_table_get(spec, "lambda", i, NULL);
            dwn = own - nwn;
            if (own > 0 && dwn > maxstep) {
                maxstep = dwn;
            }
            own = nwn;
        }
    }

    /* Criterion for single model spectrum */
    p = cpl_parameterlist_find(drvpar->parlist, "singlespec");
    if (maxstep > 2 * MF_EXTRACOVER) {
        cpl_parameter_set_int(p, 0);
        return CPL_ERROR_NONE;
    }
    cpl_parameter_set_int(p, 1);

    /* Get number of ranges */
    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    nrange = cpl_parameter_get_int(p);

    /* Set all obsolete parameters to -1 */
    for (j = 1; j < nrange; j++) {
        cpl_table_set(drvpar->rangetab, "pixres", j, -1);
        cpl_table_set(drvpar->rangetab, "wn_start", j, -1);
        cpl_table_set(drvpar->rangetab, "wn_end", j, -1);
        cpl_table_set(drvpar->rangetab, "wn_step", j, -1);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_readspec_setweight(cpl_table *spec, const mfdrv *drvpar)
{
    /*!
     * Sets the weights for spectra that do not contain an error column.
     * The default relative error provided by the ::mfdrv parameter structure
     * is multiplied by the mean flux of all wavelengths that are used by the
     * fit. The resulting absolute error is taken for all wavelengths with
     * non-zero weight.
     *
     * \b INPUT:
     * \param spec    CPL table with observed spectrum
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param spec    spectrum with optimised weights
     *
     * \b ERRORS:
     * - Invalid object value(s)
     */

    cpl_error_code err = CPL_ERROR_NONE;
    cpl_parameter *p;
    char col_dflux[MF_LENLINE+2], errtxt[MF_MAXLEN];
    int nrow = 0, i = 0, nw = 0;
    double deferr = 0., weight = 0., fsum = 0., weight0 = 0.;

    /* Return if error data exist */
    p = cpl_parameterlist_find(drvpar->parlist, "col_dflux");
    strncpy(col_dflux, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_dflux, "NULL") != 0) {
        return CPL_ERROR_NONE;
    }

    /* Get default relative error */
    p = cpl_parameterlist_find(drvpar->parlist, "default_error");
    deferr = cpl_parameter_get_double(p);

    /* Count data points with non-zero weight and sum up their fluxes */
    nrow = cpl_table_get_nrow(spec);
    for (i = 0; i < nrow; i++) {
        weight = cpl_table_get(spec, "weight", i, NULL);
        if (weight > 0) {
            nw++;
            fsum += cpl_table_get(spec, "flux", i, NULL);
        }
    }

    /* Set errors to default error * mean if error column is missing */
    if (fsum > 0) {
        weight0 = (double) nw / (fsum * deferr);
        for (i = 0; i < nrow; i++) {
            weight = cpl_table_get(spec, "weight", i, NULL);
            if (weight > 0) {
                cpl_table_set(spec, "weight", i, weight0);
            }
        }
    }

    /* Check for errors */

    if (fsum == 0) {
        sprintf(errtxt, "%s: cpl_table *spec (all fluxes = 0)",
                MF_ERROR_IOV_TXT);
        err = cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    if (nw == 0) {
        sprintf(errtxt, "%s: cpl_table *spec (all weights = 0)",
                MF_ERROR_IOV_TXT);
        err = cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    return err;
}


cpl_error_code mf_readspec_calcwaverange(mfdrv *drvpar, cpl_table *spec)
{
    /*!
     * Calculates the model wavenumber limits and step sizes for the different
     * selected wavelength ranges and stores this information in the ::mfdrv
     * parameter structure.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     * \param spec    CPL table with observed spectrum
     *
     * \b OUTPUT:
     * \param drvpar  ::mfdrv parameter structure with wavenumber limits and
     *                and step sizes
     *
     * \b ERRORS:
     * - Invalid object value(s)
     * - No data
     */

    cpl_parameter *p;
    cpl_table *rangespec;
    char errtxt[MF_MAXLEN];
    int m = 0, nrange = 0, singlespec = 0, j = 0, nsel = 0, i = 0;
    double slitw = 0., pixsc = 0., respix = 0., wn_start = 0., wn_end = 0.;
    double hres = 0., olam = 0., nlam = 0., res = 0., wn_step = 0., dwn = 0.;
    double npix = 0.;

    /* Get slit width in pixels
       -> rough estimate/lower limit of resolution in pixels */
    p = cpl_parameterlist_find(drvpar->parlist, "slitw");
    slitw = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find(drvpar->parlist, "pixsc");
    pixsc = cpl_parameter_get_double(p);
    if (pixsc == 0.) {
        sprintf(errtxt, "%s: pixsc of mfdrv *drvpar = 0", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }
    respix = slitw / pixsc;

    /* Require at least 2 pixels per resolution element */
    if (respix < 2.) {
        respix = 2.;
    }

    /* Check existence of data */
    m = cpl_table_get_nrow(spec);
    if (m <= 0) {
        sprintf(errtxt, "%s: cpl_table *spec", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Get number of wavelength ranges */
    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    nrange = cpl_parameter_get_int(p);

    /* Calculate one molecular spectrum for full wavelength range? */
    p = cpl_parameterlist_find(drvpar->parlist, "singlespec");
    singlespec = cpl_parameter_get_int(p);
    if (singlespec == 1) {
        nrange = 1;
    }

    /* Calculate input parameters for radiative transfer code */

    for (j = 0; j < nrange; j++) {

        /* Extract wavelength grid of observed spectrum for each range */
        cpl_table_unselect_all(spec);
        if (singlespec == 1) {
            cpl_table_or_selected_int(spec, "mrange", CPL_NOT_EQUAL_TO, 0);
        } else {
            cpl_table_or_selected_int(spec, "mrange", CPL_EQUAL_TO, j+1);
        }
        rangespec = cpl_table_extract_selected(spec);
        cpl_table_select_all(spec);
        nsel = cpl_table_get_nrow(rangespec);

        /* Case: all weights = 0 */
        if (cpl_table_get_column_max(rangespec, "weight") == 0.) {
            cpl_msg_info(cpl_func, "Range %d: all weights = 0", j+1);
            cpl_table_set(drvpar->rangetab, "pixres", j, 0.);
            cpl_table_set(drvpar->rangetab, "wn_start", j, 0.);
            cpl_table_set(drvpar->rangetab, "wn_end", j, 0.);
            cpl_table_set(drvpar->rangetab, "wn_step", j, 0.);
            cpl_table_delete(rangespec);
            continue;
        }

        /* Derive wavenumber range in cm^-1 */
        wn_start = 1e4 / cpl_table_get(rangespec, "lambda", nsel-1, NULL);
        wn_end = 1e4 / cpl_table_get(rangespec, "lambda", 0, NULL);

        /* Consider additional spectral coverage */
        wn_start -= MF_EXTRACOVER;
        wn_end += MF_EXTRACOVER;

        /* Find highest resolution */
        for (hres = 0., olam = 0., i = 0; i < nsel; i++) {
            nlam = cpl_table_get(rangespec, "lambda", i, NULL);
            res = nlam / ((nlam - olam) * respix);
            if (res > hres) {
              hres = res;
            }
            olam = nlam;
        }

        /* Derive wavenumber step */
        wn_step = wn_end / (hres * respix * MF_SAMPFAC);
        dwn = wn_end - wn_start;
        npix = dwn / wn_step;
        wn_step *= npix / ceil(npix);

        /* Write pixel resolution, wavenumber range and step into mfdrv
           structure */
        cpl_table_set(drvpar->rangetab, "pixres", j, hres * respix);
        cpl_table_set(drvpar->rangetab, "wn_start", j, wn_start);
        cpl_table_set(drvpar->rangetab, "wn_end", j, wn_end);
        cpl_table_set(drvpar->rangetab, "wn_step", j, wn_step);

        /* Delete temporary table */
        cpl_table_delete(rangespec);

    }

   return CPL_ERROR_NONE;
}

/**@}*/
