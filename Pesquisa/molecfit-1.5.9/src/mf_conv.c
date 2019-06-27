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
 * \file mf_conv.c
 *
 * Routines for conversion of files
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  25 Jul 2012
 * \date   29 Jan 2014
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <mf_conv.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code mf_conv_preptable(const char *parfile_in)
{
    /*!
     * \callgraph
     *
     * Converts input file into a FITS table for MOLECFIT.
     *
     * \b INPUT:
     * \param parfile_in  name of parameter file
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    mfdrv drvpar;
    mftarr tabdat;
    char parfile[MF_MAXLEN] = "";

    /* Read MOLECFIT driver file */
    mf_par_initall(&drvpar);
    mf_basic_absfile(parfile, parfile_in);
    if ((status = mf_par_readfile(&drvpar, parfile, 0)) != CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        return status;
    }

    /* Read file and convert it to CPL table and CPL property list */
    if ((status = mf_conv_readfile(&tabdat, &drvpar)) !=
        CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        return status;
    }

    /* Write CPL table and CPL property list to FITS table */
    mf_conv_writetable(&tabdat, &drvpar);

    /* Free allocated memory */
    mf_par_deleteall(&drvpar);
    mf_conv_tarr_delete(&tabdat);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code mf_conv_readfile(mftarr *tabdat, const mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Reads an ASCII file or FITS file (either FITS table or 1D fits image)
     * and write its/their data into a ::mftarr structure, which contains an
     * array of CPL tables (spectra) and CPL property lists (header keywords).
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param tabdat  ::mftarr structure with data of input file
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     * - see subroutines
     */

    cpl_parameter *p;
    cpl_array *colnames = NULL;
    cpl_table *extnames = NULL;
    mfvarr vecdat;
    char filename[MF_MAXLEN] = "";
    char relfilename[MF_MAXLEN] = "", errtxt[MF_MAXLEN] = "";
    int fitsformat = 0;

    /* Get file name and write info message */
    p = cpl_parameterlist_find(drvpar->parlist, "filename");
    strncpy(filename, cpl_parameter_get_string(p), MF_MAXLEN);
    cpl_msg_info(cpl_func, "Input data file: %s", filename);
    if (filename[0] != '/') {
        char curdir[MF_MAXLEN];
        p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        sprintf(relfilename, "%s/%s", curdir, filename);
        mf_basic_absfile(filename, relfilename);
    }

    /* Get file type */
    mf_conv_checkfitsformat(&fitsformat, filename);

    /* Read file dependent on type */
    if (fitsformat == 0) {
        /* Read non-FITS file (ASCII format assumed) */
        colnames = cpl_array_new(0, CPL_TYPE_STRING);
        mf_conv_setcolnames(colnames, drvpar);
        mf_conv_ascii_read(tabdat, filename, colnames);
        cpl_array_delete(colnames);
    } else if (fitsformat == 1) {
        /* Read FITS table */
        mf_conv_tarr_read(tabdat, filename);
    } else if (fitsformat == 2) {
        /* Read 1D FITS image */
        mf_conv_varr_read(&vecdat, filename);
        extnames = cpl_table_new(0);
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            mf_conv_setextnames(extnames, &vecdat, drvpar);
        }
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            mf_conv_varr2tarr(tabdat, &vecdat, extnames);
        }
        mf_conv_varr_delete(&vecdat);
        cpl_table_delete(extnames);
    } else if (fitsformat > 2) {
        sprintf(errtxt, "%s: %s (multidimensional FITS image)",
                MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Adapt data for use by MOLECFIT */
    if (cpl_error_get_code() == CPL_ERROR_NONE) {
        mf_conv_modtable(tabdat, drvpar);
    }

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code mf_conv_checkfitsformat(int *fitsformat, const char *filename)
{
    /*!
     * Identifies the FITS type (table or image) by means of the header
     * keyword XTENSION of the first FITS extension. If an extension is not
     * present, the image format is assumed. The keyword NAXIS is used to
     * distinguish between 1D, 2D, and 3D images. The routine returns 0 for
     * non-FITS format (e.g. ASCII), 1 for FITS table, 2 for 1D FITS image, 3
     * for 2D FITS image, 4 for 3D FITS image, and -1 in the case of errors.
     * It is assumed that all extensions have the same FITS type.
     *
     * \b INPUT:
     * \param filename    path and name of input FITS file
     *
     * \b OUTPUT:
     * \param fitsformat  flag for FITS format
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    cpl_propertylist *header = NULL;
    char errtxt[MF_MAXLEN] = "";
    int next = 0;
    cpl_errorstate prestate = cpl_errorstate_get();

    /* Check file existence */
    if ((stream = fopen(filename, "r")) == NULL) {
        *fitsformat = -1;
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }
    fclose(stream);

    /* Get number of extensions */
    next = cpl_fits_count_extensions(filename);

    /* FITS table or FITS image? */

    if (next == -1) {

        /* Not a FITS file */
        cpl_errorstate_set(prestate);
        *fitsformat = 0;

    } else if (next == 0) {

        /* File header from zeroth extension */
        header = cpl_propertylist_load(filename, 0);
        if (!header) return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "Header wrong, it couldn't be load [file=%s, ext=0]", filename);

        /* FITS tables require at least one extension */
        *fitsformat = 2;

    } else {

        cpl_property *prop;

        /* File header from first extension */
        header = cpl_propertylist_load(filename, 1);
        if (!header) return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "Header wrong, it couldn't be load [file=%s, ext=1]", filename);

        /* Read FITS keyword XTENSION */
        prop = cpl_propertylist_get_property(header, "XTENSION");
        if (prop == NULL) {
            *fitsformat = -1;
            cpl_propertylist_delete(header);
            sprintf(errtxt, "%s: %s (keyword XTENSION not found)",
                        MF_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        } else {
            if (strcmp(cpl_property_get_string(prop), "IMAGE") == 0) {
                /* FITS image */
                *fitsformat = 2;
            } else if (strcmp(cpl_property_get_string(prop), "BINTABLE")
                       == 0) {
                /* FITS table */
                *fitsformat = 1;
            } else {
                *fitsformat = -1;
                cpl_propertylist_delete(header);
                sprintf(errtxt, "%s: %s (keyword XTENSION != IMAGE or "
                        "BINTABLE)", MF_ERROR_UFS_TXT, filename);
                return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                             errtxt);
            }
        }
    }

    /* 1D or 2D FITS image? */

    if (*fitsformat == 2) {
        cpl_property *prop = cpl_propertylist_get_property(header, "NAXIS");
        if (prop == NULL) {
            *fitsformat = -1;
            cpl_propertylist_delete(header);
            sprintf(errtxt, "%s: %s (keyword NAXIS not found)",
                        MF_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        } else {
            if (cpl_property_get_int(prop) == 1) {
                /* 1D FITS image */
                *fitsformat = 2;
            } else if (cpl_property_get_int(prop) == 2) {
                /* 2D FITS image */
                *fitsformat = 3;
            } else if (cpl_property_get_int(prop) == 3) {
                /* 3D FITS image */
                *fitsformat = 4;
            } else {
                *fitsformat = -1;
                cpl_propertylist_delete(header);
                sprintf(errtxt, "%s: %s (keyword NAXIS != 1, 2, or 3)",
                        MF_ERROR_UFS_TXT, filename);
                return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                             errtxt);
            }
        }

    }

    /* Free allocated memory */
    cpl_propertylist_delete(header);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code mf_conv_setcolnames(cpl_array *colnames, const mfdrv *drvpar)
{
    /*!
     * Gets column names for ASCII file from ::mfdrv parameter structure. The
     * mandatory columns are wavelength and flux. Optionally, error and mask
     * can also be provided. If this is not the case, this is indicated by
     * 'NULL' in the parameter structure.
     *
     * \b INPUT:
     * \param colnames  empty CPL array for strings
     * \param drvpar    ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param colnames  CPL array of column names for reading of ASCII file
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    char col_lam[MF_LENLINE+2] = "", col_flux[MF_LENLINE+2] = "";
    char col_dflux[MF_LENLINE+2] = "", col_mask[MF_LENLINE+2] = "";
    int j = 0;

    /* Set default size of array for column names */
    cpl_array_set_size(colnames, 4);

    /* Get names of wavelength and flux column */
    p = cpl_parameterlist_find(drvpar->parlist, "col_lam");
    strncpy(col_lam, cpl_parameter_get_string(p), MF_LENLINE+2);
    cpl_array_set_string(colnames, j, col_lam);
    j++;
    p = cpl_parameterlist_find(drvpar->parlist, "col_flux");
    strncpy(col_flux, cpl_parameter_get_string(p), MF_LENLINE+2);
    cpl_array_set_string(colnames, j, col_flux);
    j++;

    /* Get information on existence and names of error and mask column */
    p = cpl_parameterlist_find(drvpar->parlist, "col_dflux");
    strncpy(col_dflux, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_dflux, "NULL") != 0) {
        cpl_array_set_string(colnames, j, col_dflux);
        j++;
    }
    p = cpl_parameterlist_find(drvpar->parlist, "col_mask");
    strncpy(col_mask, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_mask, "NULL") != 0) {
        cpl_array_set_string(colnames, j, col_mask);
        j++;
    }

    /* Resize array */
    cpl_array_set_size(colnames, j);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_setextnames(cpl_table *extnames, const mfvarr *vecdat,
                                   const mfdrv *drvpar)
{
    /*!
     * Prepares CPL table with FITS extension numbers and names for conversion
     * between ::mfvarr and ::mftarr structures. The column for the extension
     * names is also used for the table column names. Extension/column names
     * are provided by the \e columns parameter of the driver file.
     *
     * \b INPUT:
     * \param extnames  empty CPL table
     * \param vecdat    ::mfvarr structure with data of 1D FITS image
     * \param drvpar    ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b ERRORS:
     * - Invalid object structure
     * - see subroutines
     */

    cpl_parameter *p;
    cpl_boolean none_lam = CPL_FALSE, none_flux = CPL_FALSE;
    char errtxt[MF_MAXLEN] = "", colname[MF_LENLINE+2] = "";
    char extname[MF_LENLINE+2] = "";
    char **col;
    int ncol0 = 4, ncol = 0, next = 0, check = 0, j = 0, h = 0;
    int *extn;

    /* Set default size of output table */
    cpl_table_set_size(extnames, ncol0);

    /* Create columns for extension numbers and names */
    if (cpl_table_has_column(extnames, "col") != 1) {
        cpl_table_new_column(extnames, "col", CPL_TYPE_STRING);
    }
    if (cpl_table_has_column(extnames, "extn") != 1) {
        cpl_table_new_column(extnames, "extn", CPL_TYPE_INT);
    }
    cpl_table_fill_column_window_int(extnames, "extn", 0, ncol0, -1);

    /* Get column/extension names from parameter list */

    ncol = ncol0;

    p = cpl_parameterlist_find(drvpar->parlist, "col_lam");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 0, colname);
    if (strcmp(colname, "NULL") == 0) {
        cpl_table_set_string(extnames, "col", 0, MF_DEFLAMCOL);
        ncol--;
        none_lam = CPL_TRUE;
    }

    p = cpl_parameterlist_find(drvpar->parlist, "col_flux");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 1, colname);
    if (strcmp(colname, "NULL") == 0) {
        cpl_table_set_string(extnames, "col", 1, MF_DEFFLUXCOL);
        none_flux = CPL_TRUE;
    }

    p = cpl_parameterlist_find(drvpar->parlist, "col_dflux");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 2, colname);
    if (strcmp(colname, "NULL") == 0) {
        ncol--;
    }

    p = cpl_parameterlist_find(drvpar->parlist, "col_mask");
    sprintf(colname, "%s", cpl_parameter_get_string(p));
    cpl_table_set_string(extnames, "col", 3, colname);
    if (strcmp(colname, "NULL") == 0) {
        ncol--;
    }

    /* Get pointers to table columns */
    col = cpl_table_get_data_string(extnames, "col");
    extn = cpl_table_get_data_int(extnames, "extn");

    /* Get number of extensions */
    next = vecdat->next;

    /* Check existence of extensions with expected names and set extension
       numbers */

    for (check = 0, j = 0; j <= next; j++) {
        cpl_property *prop =
            cpl_propertylist_get_property(vecdat->head[j], "EXTNAME");
        if (prop == NULL && j > 0) {
            continue;
        } else if (prop == NULL && j == 0) {
            /* Make sure that flux vector is counted independent of existence
               of keyword EXTNAME (requires 'NULL' as extension name)  */
            if (none_flux == CPL_TRUE) {
                extn[1] = 0;
                check++;
            }
        } else {
            strncpy(extname, cpl_property_get_string(prop), MF_LENLINE+2);
            for (h = 0; h < ncol0; h++) {
                if (strncmp(extname, col[h], MF_LENLINE+2) == 0) {
                    extn[h] = j;
                    check++;
                }
            }
        }
    }

    if (check != ncol) {
        if (extn[0] < 0 && none_lam == CPL_FALSE) {
            sprintf(errtxt, "%s: mfvarr *vecdat (wavelength extension '%s' "
                    "not found)", MF_ERROR_IOS_TXT, col[0]);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (extn[1] < 0 && none_flux == CPL_FALSE) {
            sprintf(errtxt, "%s: mfvarr *vecdat (flux extension '%s' not "
                    "found)", MF_ERROR_IOS_TXT, col[1]);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (extn[2] < 0 && strcmp(col[2], "NULL") != 0) {
            sprintf(errtxt, "%s: mfvarr *vecdat (flux error extension '%s' "
                    "not found)", MF_ERROR_IOS_TXT, col[2]);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (extn[3] < 0 && strcmp(col[3], "NULL") != 0) {
            sprintf(errtxt, "%s: mfvarr *vecdat (mask extension '%s' not "
                    "found)", MF_ERROR_IOS_TXT, col[3]);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_varr2tarr(mftarr *tabdat, const mfvarr *vecdat,
                                 const cpl_table *extnames)
{
    /*!
     * Fills CPL tables and CPL property lists of a ::mftarr structure by data
     * read from a 1D FITS image with optional extensions for error, mask, and
     * wavelength. The image data is provided by a ::mfvarr structure. The
     * extension numbers and names are given by a CPL table. Extensions are
     * not read if names are missing or wrong. Images for different chips have
     * to be converted separately.
     *
     * \b INPUT:
     * \param vecdat    ::mfvarr structure with data of 1D FITS image
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b OUTPUT:
     * \param tabdat    ::mftarr structure with data of FITS table
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     * - Invalid object value(s)
     */

    char errtxt[MF_MAXLEN] = "";
    const char **col;
    int ncol = 0, next = 0, h = 0, nrow = 0, j = 0, i = 0;
    const int *extn;
    double crpix = 0., crval = 0., cdelt = 0.;

    /* Initialise CPL tables and CPL property lists for content of mfvarr
       structure (memory allocation) */
    mf_conv_tarr_init(tabdat, 1);
    tabdat->tab[0] = cpl_table_new(0);
    tabdat->tab[1] = cpl_table_new(0);

    /* Transfer FITS header data */
    tabdat->head[0] = cpl_propertylist_duplicate(vecdat->head[0]);
    tabdat->head[1] = cpl_propertylist_new();

    /* Check table for extension numbers and names */
    ncol = cpl_table_get_nrow(extnames);
    if (ncol <= 0) {
        mf_conv_tarr_delete(tabdat);
        sprintf(errtxt, "%s: cpl_table *extnames", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check column names */
    if (cpl_table_has_column(extnames, "col") != 1 ||
        cpl_table_has_column(extnames, "extn") != 1) {
        mf_conv_tarr_delete(tabdat);
        sprintf(errtxt, "%s: cpl_table *extnames (column names)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of data */
    if (cpl_table_get_column_max(extnames, "extn") == -1) {
        mf_conv_tarr_delete(tabdat);
        sprintf(errtxt, "%s: cpl_table *extnames (no extension with valid "
                "data)", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Get number of extensions */
    next = vecdat->next;

    /* Check extension numbers */
    if (cpl_table_get_column_max(extnames, "extn") > next) {
        mf_conv_tarr_delete(tabdat);
        sprintf(errtxt, "%s: cpl_table *extnames (incorrect extension "
                "number)", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Get pointers to table columns */
    col = cpl_table_get_data_string_const(extnames, "col");
    extn = cpl_table_get_data_int_const(extnames, "extn");

    /* Create table columns */
    cpl_table_new_column(tabdat->tab[1], col[0], CPL_TYPE_DOUBLE);
    for (h = 1; h < ncol; h++) {
        if (extn[h] >= 0) {
            cpl_table_new_column(tabdat->tab[1], col[h], CPL_TYPE_DOUBLE);
        }
    }

    /* Set size of data table */
    nrow = cpl_vector_get_size(vecdat->vec[0]);
    if (nrow <= 0) {
        mf_conv_tarr_delete(tabdat);
        sprintf(errtxt, "%s: mfvarr *vecdat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    } else if (nrow > 0 && next > 0) {
        for (j = 1; j <= next; j++) {
            if (cpl_vector_get_size(vecdat->vec[j]) != nrow) {
                mf_conv_tarr_delete(tabdat);
                sprintf(errtxt, "%s: mfvarr *vecdat (vector size differs for "
                        "different extensions", MF_ERROR_IOS_TXT);
                return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                             errtxt);
            }
        }
    }
    cpl_table_set_size(tabdat->tab[1], nrow);

    /* Get wavelength grid from FITS header if a wavelength vector is not
       provided */
    if (extn[0] < 0) {
         crpix = mf_conv_getwcskey(vecdat->head[0], "CRPIX1");
         if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
             mf_conv_tarr_delete(tabdat);
             sprintf(errtxt, "%s: mfvarr *vecdat (keyword CRPIX1 not found)",
                     MF_ERROR_IOS_TXT);
             return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                          errtxt);
         }
         crval = mf_conv_getwcskey(vecdat->head[0], "CRVAL1");
         if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
             mf_conv_tarr_delete(tabdat);
             sprintf(errtxt, "%s: mfvarr *vecdat (keyword CRVAL1 not found)",
                     MF_ERROR_IOS_TXT);
             return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                          errtxt);
         }
         cdelt = mf_conv_getwcskey(vecdat->head[0], "CDELT1");
         if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
             mf_conv_tarr_delete(tabdat);
             sprintf(errtxt, "%s: mfvarr *vecdat (keyword CDELT1 not found)",
                     MF_ERROR_IOS_TXT);
             return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                          errtxt);
         }
    }

    /* Convert vectors to table columns */
    for (h = 0; h < ncol; h++) {
        for (i = 0; i < nrow; i++) {
            if (extn[h] >= 0) {
                cpl_table_set(tabdat->tab[1], col[h], i,
                              cpl_vector_get(vecdat->vec[extn[h]], i));
            } else if (extn[0] < 0 && h == 0) {
                cpl_table_set(tabdat->tab[1], col[h], i,
                              crval + (i - crpix + 1) * cdelt);
            }
        }
    }

    return CPL_ERROR_NONE;
}


double mf_conv_getwcskey(const cpl_propertylist *plist, const char *key)
{
    /*!
     * Retrieves a FITS WCS key from a property list
     *
     * \b INPUT:
     * \param plist  cpl_propertylist
     * \param key    key to retrieve
     *
     * \b RETURN:
     * - FITS WCS key
     *
     * \b ERRORS:
     * - CPL_ERROR_DATA_NOT_FOUND: key does not exist
     */

    double d;
    const cpl_errorstate cleanstate = cpl_errorstate_get();

    d = cpl_propertylist_get_double(plist, key);

    /* WCS keys can be written as ints but should be interpreted as doubles */
    if (cpl_error_get_code() == CPL_ERROR_TYPE_MISMATCH) {
        cpl_errorstate_set(cleanstate);
        d = (double) cpl_propertylist_get_int(plist, key);
    }

    return d;
}


static cpl_table * mf_convert_sdp_table(const cpl_table * inptable,
                                const char * col_lam, const char * col_flux,
                                const char * col_dflux, const char * col_mask)
{
    /*!
     * Convert SDP formatted VOCLASS Spectrum data to internal format.
     * This format consists of one row containing arrays with the data in the
     * columns
     *
     * \b INPUT:
     * \param inptable   input Spectrum table with one row containing arrays
     * \param col_lam    name of wavelength column
     * \param col_flux   name of flux column
     * \param col_dflux  name of flux error column or NULL
     * \param col_mask   name of mask column or NULL
     *
     * \b OUTPUT:
     * new table in with arrays expanded to rows
     *
     * \b ERRORS:
     * - Missing or inconsistent columns
     */
    const cpl_array * alam = cpl_table_get_array(inptable, col_lam, 0);
    if (alam == NULL) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Expected %s "
                              "array column in input", col_lam);;
        return NULL;
    }
    const cpl_array * aflux = cpl_table_get_array(inptable, col_flux, 0);
    if (aflux == NULL) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Expected %s "
                              "array column in input", col_flux);
        return NULL;
    }
    const cpl_array * adflux = NULL, * amask = NULL;
    if (col_dflux) {
        adflux = cpl_table_get_array(inptable, col_dflux, 0);
        if (adflux == NULL) {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Expected %s "
                                  "array column in input", col_dflux);
            return NULL;
        }
    }
    if (col_mask) {
        amask = cpl_table_get_array(inptable, col_mask, 0);
        if (adflux == NULL) {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT, "Expected %s "
                                  "array column in input", col_mask);
            return NULL;
        }
    }
    if (cpl_array_get_size(alam) != cpl_array_get_size(aflux) ||
        (adflux && cpl_array_get_size(alam) != cpl_array_get_size(adflux)) ||
        (amask && cpl_array_get_size(alam) != cpl_array_get_size(amask))) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Table columns do not all have the same size");
        return NULL;
    }

    cpl_table * ntab = cpl_table_new(cpl_array_get_size(alam));

    /* const casts due to PIPE-6075 */
    CPL_DIAG_PRAGMA_PUSH_IGN(-Wcast-qual);
      cpl_array *dalam = cpl_array_cast((cpl_array*)alam, CPL_TYPE_DOUBLE);
    CPL_DIAG_PRAGMA_POP;
    cpl_table_new_column(ntab, col_lam, CPL_TYPE_DOUBLE);
    cpl_table_copy_data_double(ntab, col_lam, cpl_array_get_data_double_const(dalam));
    cpl_array_delete(dalam);

    CPL_DIAG_PRAGMA_PUSH_IGN(-Wcast-qual);
      cpl_array *daflux = cpl_array_cast((cpl_array*)aflux, CPL_TYPE_DOUBLE);
    CPL_DIAG_PRAGMA_POP;
    cpl_table_new_column(ntab, col_flux, CPL_TYPE_DOUBLE);
    cpl_table_copy_data_double(ntab, col_flux, cpl_array_get_data_double_const(daflux));
    cpl_array_delete(daflux);

    if(adflux) {
        CPL_DIAG_PRAGMA_PUSH_IGN(-Wcast-qual);
          cpl_array *dadflux = cpl_array_cast((cpl_array*)adflux, CPL_TYPE_DOUBLE);
        CPL_DIAG_PRAGMA_POP;
        cpl_table_new_column(ntab, col_dflux, CPL_TYPE_DOUBLE);
        cpl_table_copy_data_double(ntab, col_dflux, cpl_array_get_data_double_const(dadflux));
        cpl_array_delete(dadflux);
    }
    if (amask) {
        cpl_table_new_column(ntab, col_mask, CPL_TYPE_INT);
        cpl_table_copy_data_int(ntab, col_mask, cpl_array_get_data_int_const(amask));
    }

    return ntab;
}


cpl_error_code mf_conv_modtable(mftarr *tabdat, const mfdrv *drvpar)
{
    /*!
     * Modifies data in ::mftarr structure in order to be consistent with
     * requirements of MOLECFIT.
     * The routine creates a column with the integer mask values 0 (rejected)
     * and 1 (ok). By default, it is expected that input mask values agree
     * with this definition. If other values are found, 0 is assumed to be ok
     * and all other values cause pixel rejection. Possible nan in the flux
     * and error columns are substituted by zero flux and indicated by
     * mask = 0. The same is performed for negative errors.
     *
     * \b INPUT:
     * \param tabdat  ::mftarr structure with read data
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param tabdat  ::mftarr structure with modified data
     *
     * \b ERRORS:
     * - Invalid object structure
     */

    cpl_parameter *p;
    char col_lam[MF_LENLINE+2] = "", col_flux[MF_LENLINE+2] = "";
    char col_dflux[MF_LENLINE+2] = "", col_mask[MF_LENLINE+2] = "";
    char errtxt[MF_MAXLEN] = "", col_imask[MF_LENLINE+2] = "";
    char col_renamed[MF_LENLINE+2] = "";
    cpl_boolean exerr = CPL_TRUE, exmask = CPL_TRUE, isnanflux = CPL_FALSE;
    cpl_boolean isnanum = CPL_FALSE, isnandflux = CPL_FALSE;
    cpl_boolean isoutrange = CPL_FALSE, isnomask = CPL_FALSE, is0 = CPL_FALSE;
    cpl_boolean ismask0 = CPL_FALSE;
    int nchip = 0, j = 0, nrow = 0, imask = 0, i = 0, nmask0 = 0;
    double flux = 0., dflux = 0., mask = 0.;

    /* Get names of wavelength and flux column */
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

    /* Get information on existence and names of error and mask column */
    p = cpl_parameterlist_find(drvpar->parlist, "col_dflux");
    strncpy(col_dflux, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_dflux, "NULL") == 0) {
        exerr = CPL_FALSE;
    }
    p = cpl_parameterlist_find(drvpar->parlist, "col_mask");
    strncpy(col_mask, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_mask, "NULL") == 0) {
        exmask = CPL_FALSE;
    }

    /* Get number of chips */
    nchip = tabdat->next;

    /* Check data for each chip */

    for (j = 0; j < nchip; j++) {

        /* Check existence of expected columns */
        if (cpl_table_has_column(tabdat->tab[j+1], col_lam) != 1) {
            sprintf(errtxt, "%s: mftarr *tabdat (wavelength column '%s' not "
                    "found in extension %d)", MF_ERROR_IOS_TXT, col_lam, j+1);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (cpl_table_has_column(tabdat->tab[j+1], col_flux) != 1) {
            sprintf(errtxt, "%s: mftarr *tabdat (flux column '%s' not found "
                    "in extension %d)", MF_ERROR_IOS_TXT, col_flux, j+1);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }
        if (exerr == TRUE) {
            if (cpl_table_has_column(tabdat->tab[j+1], col_dflux) != 1) {
                sprintf(errtxt, "%s: mftarr *tabdat (flux error column '%s' "
                        "not found in extension %d)", MF_ERROR_IOS_TXT,
                        col_dflux, j+1);
                return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                             errtxt);
            }
        }
        if (exmask == TRUE) {
            if (cpl_table_has_column(tabdat->tab[j+1], col_mask) != 1) {
                sprintf(errtxt, "%s: mftarr *tabdat (mask column '%s' not "
                        "found in extension %d)", MF_ERROR_IOS_TXT, col_mask,
                        j+1);
                return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                             errtxt);
            }
        }

        /* check for and convert sdp data format */
        if (cpl_propertylist_has(tabdat->head[j + 1], "VOCLASS") &&
            strncmp(cpl_propertylist_get_string(tabdat->head[j + 1],
                                                "VOCLASS"),
                    "SPECTRUM", strlen("SPECTRUM")) == 0) {
            cpl_table * ntab = mf_convert_sdp_table(tabdat->tab[j + 1], col_lam, col_flux,
                                                    exerr ? col_dflux : NULL,
                                                    exmask ? col_mask : NULL);
            if (cpl_error_get_code()) {
                return cpl_error_get_code();
            }
            cpl_table_delete(tabdat->tab[j + 1]);
            tabdat->tab[j + 1] = ntab;
        }

        /* Create integer mask column (and rename already existing column of
           the same name) */
        if (exmask == CPL_TRUE) {
            sprintf(col_imask, "%s_I", col_mask);
        } else {
            sprintf(col_imask, "%s_I", MF_DEFMASKCOL);
        }
        if (cpl_table_has_column(tabdat->tab[j+1], col_imask) == 1) {
            sprintf(col_renamed, "%s_orig", col_imask);
            if (cpl_table_has_column(tabdat->tab[j+1], col_renamed) != 1) {
               cpl_msg_info(cpl_func, "Name of internal integer mask already "
                         "used: Rename %s in %s", col_imask, col_renamed);
               cpl_table_name_column(tabdat->tab[j+1], col_imask,
                                     col_renamed);
            } else {
               cpl_msg_info(cpl_func, "Use of reserved mask column names: "
                            "Erase %s, keep %s", col_imask, col_renamed);
               cpl_table_erase_column(tabdat->tab[j+1], col_imask);
            }
        }
        cpl_table_new_column(tabdat->tab[j+1], col_imask, CPL_TYPE_INT);

        /* Get number of rows */
        nrow = cpl_table_get_nrow(tabdat->tab[j+1]);

        /* Check for nan values and negative errors and correct them */

        for (nmask0 = 0, i = 0; i < nrow; i++) {

            flux = cpl_table_get(tabdat->tab[j+1], col_flux, i, NULL);
            if (isnan(flux) != 0) {
                cpl_table_set(tabdat->tab[j+1], col_flux, i, 0.);
                isnanflux = CPL_TRUE;
                isnanum = CPL_TRUE;
            } else {
                isnanflux = CPL_FALSE;
            }

            if (exerr == CPL_TRUE) {
                dflux = cpl_table_get(tabdat->tab[j+1], col_dflux, i, NULL);
                if (dflux <= 0 || isnan(dflux) != 0) {
                    cpl_table_set(tabdat->tab[j+1], col_dflux, i, 0.);
                    isnandflux = CPL_TRUE;
                    isoutrange = CPL_TRUE;
                } else {
                    isnandflux = CPL_FALSE;
                }
            }

            if (exmask == CPL_TRUE) {
                mask = cpl_table_get(tabdat->tab[j+1], col_mask, i, NULL);
            } else {
                mask = 1;
            }
            if (isnanflux == CPL_TRUE || isnandflux == CPL_TRUE) {
                cpl_table_set(tabdat->tab[j+1], col_imask, i, -1);
                if (mask != 0 && mask != 1) {
                    isnomask = CPL_TRUE;
                } else if (mask == 0) {
                    is0 = CPL_TRUE;
                    nmask0++;
                }
            } else {
                if (mask == 0) {
                    cpl_table_set(tabdat->tab[j+1], col_imask, i, 0);
                    is0 = CPL_TRUE;
                    nmask0++;
                } else {
                    cpl_table_set(tabdat->tab[j+1], col_imask, i, 1);
                    if (mask != 1) {
                        isnomask = CPL_TRUE;
                    }
                }
            }

        }

    }

    /* Return if mask cannot be interpreted */
    if (isnomask == CPL_TRUE && is0 == CPL_FALSE) {
        sprintf(errtxt, "%s: mftarr *tabdat (all mask value(s) != 0 or 1)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Correct integer mask values if required (e.g. reverse definition) */

    for (j = 0; j < nchip; j++) {
        nrow = cpl_table_get_nrow(tabdat->tab[j+1]);
        for (i = 0; i < nrow; i++) {
            imask = cpl_table_get(tabdat->tab[j+1], col_imask, i, NULL);
            if (isnomask == CPL_TRUE) {
                if (imask == 0) {
                    cpl_table_set(tabdat->tab[j+1], col_imask, i, 1);
                } else {
                    cpl_table_set(tabdat->tab[j+1], col_imask, i, 0);
                }
            } else {
                if (imask == -1) {
                    cpl_table_set(tabdat->tab[j+1], col_imask, i, 0);
                } else if (nmask0 == nrow) {
                    cpl_table_set(tabdat->tab[j+1], col_imask, i, 1);
                    ismask0 = CPL_TRUE;
                }
            }
        }
    }

    /* Print info message in the case of bad fluxes, errors, or mask
       values */

    if (isnanum == CPL_TRUE) {
        cpl_msg_info(cpl_func, "Input data: flux(es) = 'nan' "
                     "-> set mask = 0");
    }

    if (isoutrange == CPL_TRUE) {
        cpl_msg_info(cpl_func, "Input data: error(s) <= 0 or 'nan' "
                     "-> set mask = 0");
    }

    if (isnomask == CPL_TRUE) {
        cpl_msg_info(cpl_func, "Input data: mask value(s) != 0 or 1 "
                     "-> reverse definition (0 -> 1; != 0 -> 0)");
    }

    if (ismask0 == CPL_TRUE) {
        cpl_msg_info(cpl_func, "Input data: all mask values = 0 "
                     "-> reverse definition (0 -> 1)");
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_writetable(const mftarr *tabdat, const mfdrv *drvpar)
{
    /*!
     * Writes CPL tables (spectra) and CPL property lists (header keywords)
     * for each chip from an ::mftarr structure into a FITS table. The number
     * of FITS extensions equals the number of chips.
     *
     * \b INPUT:
     * \param tabdat  ::mftarr structure with FITS data
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    char basedir[MF_MAXLEN] = "", outdir[MF_MAXLEN] = "";
    char outname[MF_MAXLEN] = "", outfile[MF_MAXLEN] = "";

    /* Get output file name */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strncpy(basedir, cpl_parameter_get_string(p), MF_MAXLEN);
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);
    sprintf(outfile, "%s%s.fits", outdir, outname);

    /* Write info message */
    cpl_msg_info(cpl_func, "Convert input data file into %s", outfile);

    /* Write new FITS table for first chip and new extension for all chips */
    mf_conv_tarr_write(outfile, tabdat);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_writeresults(const char *parfile_in)
{
    /*!
     * \callgraph
     *
     * Writes results of MOLECFIT and CALCTRANS into a file of the same format
     * as the input data file. In the case of FITS images, the transmission
     * curve and the telluric absorption corrected spectrum are written into
     * separate files.
     *
     * \b INPUT:
     * \param parfile_in  name of parameter file
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_table *results = NULL;
    mfdrv drvpar;
    char parfile[MF_MAXLEN] = "";

    /* Read MOLECFIT driver file */
    mf_par_initall(&drvpar);
    mf_basic_absfile(parfile, parfile_in);
    if ((status = mf_par_readfile(&drvpar, parfile, 1)) != CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        return status;
    }

    /* Read file with CALCTRANS results and convert it into CPL table */
    results = cpl_table_new(0);
    if ((status = mf_conv_readresults(results, &drvpar)) !=
        CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        cpl_table_delete(results);
        return status;
    }

    /* Write CPL table and CPL property list to file with format of input data
       file */
    mf_conv_writefile(results, &drvpar);

    /* Free allocated memory */
    mf_par_deleteall(&drvpar);
    cpl_table_delete(results);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code mf_conv_readresults(cpl_table *results, const mfdrv *drvpar)
{
    /*!
     * Fills CPL table with the results of CALCTRANS saved as a FITS table.
     * Converts micron into initial wavelength units if necessary.
     *
     * \b INPUT:
     * \param results  empty CPL table
     * \param drvpar   ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param results  CPL table with results of CALCTRANS
     *
     * \b ERRORS:
     * - File opening failed
     */

    cpl_parameter *p;
    cpl_table *tmptab = NULL;
    char basedir[MF_MAXLEN] = "", outdir[MF_MAXLEN] = "";
    char outname[MF_MAXLEN] = "", filename[MF_MAXLEN] = "";
    char errtxt[MF_MAXLEN] = "";
    double wlgtomicron = 0.;

    /* Get path and name of results file */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strncpy(basedir, cpl_parameter_get_string(p), MF_MAXLEN);
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);
    sprintf(filename, "%s%s_tac.fits", outdir, outname);

    /* Load FITS table into CPL table */
    tmptab = cpl_table_load(filename, 1, 0);
    if (tmptab == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Copy content of temporary table into results table */
    mf_basic_copytable(results, tmptab);

    /* Get wavelength unit conversion factor */
    p = cpl_parameterlist_find(drvpar->parlist, "wlgtomicron");
    wlgtomicron = cpl_parameter_get_double(p);

    /* Change wavelength units if necessary */
    cpl_table_divide_scalar(results, "lambda", wlgtomicron);

    /* Free allocated memory */
    cpl_table_delete(tmptab);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_writefile(cpl_table *results, const mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Writes results from CALCTRANS into an ASCII file or FITS file (FITS
     * table or 1D FITS image) dependent on the format of the input data file.
     * In the case of FITS images, the telluric absorption corrected spectrum
     * (suffix: 'TAC') and the transmission curve (suffix: 'TRA') are written
     * into separate files.
     *
     * \b INPUT:
     * \param results  CPL table with results of CALCTRANS
     * \param drvpar   ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param results  results of CALCTRANS with intial wavelength grid
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     * - see subroutines
     */

    cpl_parameter *p;
    cpl_table *extnames = NULL;
    mftarr tabdat;
    mfvarr vecdat, transdat;
    char basedir[MF_MAXLEN] = "", filename[MF_MAXLEN] = "";
    char relfilename[MF_MAXLEN] = "", outdir[MF_MAXLEN] = "";
    char indir[MF_MAXLEN] = "", infilename[MF_MAXLEN] = "";
    char suffix[MF_MAXLEN] = "", outfilename[MF_MAXLEN] = "";
    char transfilename[MF_MAXLEN] = "", errtxt[MF_MAXLEN] = "";
    int fitsformat = 0;

    /* Create mftarr structure and read FITS table with data of initial input
       file */
    if (mf_conv_readprepfits(&tabdat, drvpar) != CPL_ERROR_NONE) {
        return cpl_error_get_code();
    }

    /* Write results into mftarr structure */
    if (mf_conv_results2tarr(&tabdat, results, drvpar) != CPL_ERROR_NONE) {
        mf_conv_tarr_delete(&tabdat);
        return cpl_error_get_code();
    }

    /* Get input file name */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strncpy(basedir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "filename");
    strncpy(filename, cpl_parameter_get_string(p), MF_MAXLEN);
    if (filename[0] != '/') {
        char curdir[MF_MAXLEN];
        p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        sprintf(relfilename, "%s/%s", curdir, filename);
        mf_basic_absfile(filename, relfilename);
    }

    /* Get output file names */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    mf_basic_getfilename(indir, infilename, suffix, filename);
    sprintf(outfilename, "%s%s_TAC.%s", outdir, infilename, suffix);
    sprintf(transfilename, "%s%s_TRA.%s", outdir, infilename, suffix);

    /* Write info message */
    cpl_msg_info(cpl_func, "Corrected input data file: %s", outfilename);

    /* Get file type */
    mf_conv_checkfitsformat(&fitsformat, filename);

    /* Write file dependent on input type */
    if (fitsformat == 0) {
        /* Write non-FITS file (ASCII format assumed) */
        mf_conv_erasemaskcol(&tabdat, drvpar);
        mf_conv_ascii_write(outfilename, &tabdat);
    } else if (fitsformat == 1) {
        /* Write FITS table */
        mf_conv_erasemaskcol(&tabdat, drvpar);
        mf_conv_tarr_write(outfilename, &tabdat);
    } else if (fitsformat == 2) {
        /* Write 1D FITS images */
        mf_conv_varr_read(&vecdat, filename);
        extnames = cpl_table_new(0);
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            mf_conv_setextnames(extnames, &vecdat, drvpar);
        }
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            mf_conv_resultstarr2varr(&transdat, &vecdat, extnames, &tabdat);
        }
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            mf_conv_varr_write(outfilename, &vecdat);
            cpl_msg_info(cpl_func, "Transmission frame: %s", transfilename);
            mf_conv_varr_write(transfilename, &transdat);
        }
        mf_conv_varr_delete(&vecdat);
        mf_conv_varr_delete(&transdat);
        cpl_table_delete(extnames);
    } else if (fitsformat > 2) {
        mf_conv_tarr_delete(&tabdat);
        sprintf(errtxt, "%s: %s (multidimensional FITS image)",
                MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Free allocated memory */
    mf_conv_tarr_delete(&tabdat);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code mf_conv_readprepfits(mftarr *tabdat, const mfdrv *drvpar)
{
    /*!
     * Fills CPL tables of a new ::mftarr structure by data read from a FITS
     * table prepared for MOLECFIT.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param tabdat  ::mftarr structure with data of FITS file
     *
     * \b ERRORS:
     * - see ::mf_conv_tarr_read
     */

    cpl_parameter *p;
    char basedir[MF_MAXLEN] = "", outdir[MF_MAXLEN] = "";
    char outname[MF_MAXLEN] = "", filename[MF_MAXLEN] = "";

    /* Get path and name of results file */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strncpy(basedir, cpl_parameter_get_string(p), MF_MAXLEN);
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strncpy(outname, cpl_parameter_get_string(p), MF_MAXLEN);
    sprintf(filename, "%s%s.fits", outdir, outname);

    /* Read FITS table */
    mf_conv_tarr_read(tabdat, filename);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code mf_conv_results2tarr(mftarr *tabdat, cpl_table *results,
                                    const mfdrv *drvpar)
{
    /*!
     * Writes results from CALCTRANS into an ::mftarr structure that already
     * contains data from the input file. The new table columns "mtrans",
     * "tacflux", and "tacqual" are created. They contain the data from the
     * columns "mtrans", "cflux", and "qual" of the results data table. If the
     * the ::mftarr structure has an error column, a new column named
     * "tacdflux" is created, which contains the corrected error related to
     * the "cflux" column, i.e. the input error divided by the transmission in
     * "mtrans". The "chip" column of the results data table is used to
     * transfer the data to the different extensions of the ::mftarr
     * structure.
     *
     * \b INPUT:
     * \param tabdat   ::mftarr structure with data of FITS file prepared for
     *                 MOLECFIT
     * \param results  CPL table with results of CALCTRANS
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param tabdat   input ::mftarr structure extended for CALCTRANS results
     *
     * \b ERRORS:
     * - Invalid object structure
     * - Inconsistent data grids
     */

    cpl_parameter *p;
    cpl_table *chipres = NULL;
    cpl_boolean exerr = CPL_TRUE;
    char errtxt[MF_MAXLEN] = "", col_dflux[MF_LENLINE+2] = "";
    int nchip = 0, i = 0, nrow = 0;

    /* Check existence of required columns in results data table */
    if (cpl_table_has_column(results, "mtrans") != 1 ||
        cpl_table_has_column(results, "cflux") != 1 ||
        cpl_table_has_column(results, "qual") != 1) {
        sprintf(errtxt, "%s: cpl_table *results (required column(s) not "
                "present)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get number of chips */
    nchip = tabdat->next;
    if (nchip <= 0) {
        sprintf(errtxt, "%s: mftarr *tabdat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check non-existence of results columns in input data table
       (assume same column names for all chips) */
    if (cpl_table_has_column(tabdat->tab[1], "mtrans") != 0 ||
        cpl_table_has_column(tabdat->tab[1], "tacflux") != 0 ||
        cpl_table_has_column(tabdat->tab[1], "tacdflux") != 0 ||
        cpl_table_has_column(tabdat->tab[1], "tacqual") != 0) {
        sprintf(errtxt, "%s: mftarr *tabdat (results column(s) already "
                "exist(s))", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get name of error column in input data table */
    p = cpl_parameterlist_find(drvpar->parlist, "col_dflux");
    strncpy(col_dflux, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_dflux, "NULL") == 0) {
        exerr = CPL_FALSE;
    }

    /* Prepare CPL tables for writing */

    for (i = 0; i < nchip; i++) {

        /* Select subtable with data points for a given chip */
        cpl_table_unselect_all(results);
        cpl_table_or_selected_int(results, "chip", CPL_EQUAL_TO, i+1);
        chipres = cpl_table_extract_selected(results);

        /* Get and check number of data points */
        nrow = cpl_table_get_nrow(chipres);
        if (nrow == 0) {
            cpl_table_delete(chipres);
            sprintf(errtxt, "%s: cpl_table *results (chip %d)",
                    MF_ERROR_NDA_TXT, i+1);
            return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s",
                                         errtxt);
        } else if (nrow != cpl_table_get_nrow(tabdat->tab[i+1])) {
            cpl_table_delete(chipres);
            sprintf(errtxt, "%s: cpl_table *results (chip %d) != mftarr "
                    "tabdat->tab[%d]", MF_ERROR_IDG_TXT, i+1, i+1);
            return cpl_error_set_message(cpl_func, MF_ERROR_IDG, "%s",
                                         errtxt);
        }

        /* Copy columns with CALCTRANS results into table with input data */

        /* Copy transmission and flux columns */
        cpl_table_duplicate_column(tabdat->tab[i+1], "mtrans", chipres,
                                   "mtrans");
        cpl_table_duplicate_column(tabdat->tab[i+1], "tacflux", chipres,
                                   "cflux");

        /* Create column for corrected errors if error column exists in input
           data table */
        if (exerr == CPL_TRUE) {
            cpl_table_duplicate_column(tabdat->tab[i+1], "tacdflux",
                                       tabdat->tab[i+1], col_dflux);
            cpl_table_divide_columns(tabdat->tab[i+1], "tacdflux", "mtrans");
        }

        /* Copy quality column */
        cpl_table_duplicate_column(tabdat->tab[i+1], "tacqual", chipres,
                                   "qual");

        /* Delete temporary table */
        cpl_table_delete(chipres);

    }

    /* Select all table rows again */
    cpl_table_select_all(results);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_erasemaskcol(mftarr *tabdat, const mfdrv *drvpar)
{
    /*!
     * Erases the integer mask column for MOLECFIT in an ::mftarr structure.
     * The column was created by ::mf_conv_modtable.
     *
     * \b INPUT:
     * \param tabdat  ::mftarr structure with read data
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param tabdat  ::mftarr structure without integer mask column
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    char col_mask[MF_LENLINE+2] = "", col_imask[MF_LENLINE+2] = "";
    char col_renamed[MF_LENLINE+2] = "";
    int nchip = 0, j = 0;

    /* Get name of integer mask column */
    p = cpl_parameterlist_find(drvpar->parlist, "col_mask");
    strncpy(col_mask, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(col_mask, "NULL") == 0) {
        sprintf(col_imask, "%s_I", MF_DEFMASKCOL);
    } else {
        sprintf(col_imask, "%s_I", col_mask);
    }

    /* Get number of chips */
    nchip = tabdat->next;

    /* Erase integer mask column if it exists */
    for (j = 0; j < nchip; j++) {
        if (cpl_table_has_column(tabdat->tab[j+1], col_imask) == 1) {
            cpl_table_erase_column(tabdat->tab[j+1], col_imask);
        }
    }

    /* Rename a column with the same original name as the internal integer
       mask column (temporary suffix: "_orig") */
    sprintf(col_renamed, "%s_orig", col_imask);
    for (j = 0; j < nchip; j++) {
        if (cpl_table_has_column(tabdat->tab[j+1], col_renamed) == 1) {
            cpl_table_name_column(tabdat->tab[j+1], col_renamed, col_imask);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_resultstarr2varr(mfvarr *transdat, mfvarr *vecdat,
                                        const cpl_table *extnames,
                                        const mftarr *tabdat)
{
    /*!
     * Writes telluric absorption corrected spectrum (and its error function
     * if present) from CALCTRANS into an ::mfvarr structure representing a
     * FITS image of the same format as the input data file. Data points where
     * the correction becomes unreliable are indicated if a mask extension
     * exists. Moreover, an ::mfvarr structure is created that contains the
     * transmission curve required for the correction of other spectra.
     *
     * \b INPUT:
     * \param vecdat    ::mfvarr structure with data of input 1D FITS image
     * \param extnames  CPL table with FITS extension numbers and names
     * \param tabdat    ::mftarr structure with input data and CALCTRANS
     *                  results
     *
     * \b OUTPUT:
     * \param transdat  ::mfvarr structure with model transmission curve
     * \param vecdat    ::mfvarr structure with applied telluric absorption
     *                  correction
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     * - Invalid object value(s)
     * - Inconsistent data grids
     */

    char errtxt[MF_MAXLEN] = "";
    cpl_boolean hascol = CPL_TRUE;
    int next = 0, extn_flux = 0, extn_dflux = 0, extn_mask = 0, nrow = 0;
    int j = 0, i = 0;
    int *qual = NULL;
    double maskval[2] = {0., 1.};
    double *trans = NULL, *flux = NULL, *dflux = NULL, *mask = NULL;
    double *mtrans = NULL, *cflux = NULL, *cdflux = NULL;

    /* Check table for extension numbers and names */
    if (cpl_table_get_nrow(extnames) != 4) {
        sprintf(errtxt, "%s: cpl_table *extnames (number of rows != 4)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check column names */
    if (cpl_table_has_column(extnames, "col") != 1 ||
        cpl_table_has_column(extnames, "extn") != 1) {
        sprintf(errtxt, "%s: cpl_table *extnames (column names)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of data */
    if (cpl_table_get_column_max(extnames, "extn") == -1) {
        sprintf(errtxt, "%s: cpl_table *extnames (no extension with valid "
                "data)", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Get number of extensions */
    next = vecdat->next;

    /* Check extension numbers */
    if (cpl_table_get_column_max(extnames, "extn") > next) {
        sprintf(errtxt, "%s: cpl_table *extnames (incorrect extension "
                "number)", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Get extension numbers */
    extn_flux = cpl_table_get(extnames, "extn", 1, NULL);
    extn_dflux = cpl_table_get(extnames, "extn", 2, NULL);
    extn_mask = cpl_table_get(extnames, "extn", 3, NULL);

    /* Check number of extensions in results table */
    if (tabdat->next != 1) {
        sprintf(errtxt, "%s: mftarr *tabdat (number of extensions != 1)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get number of data points */
    nrow = cpl_table_get_nrow(tabdat->tab[1]);
    if (nrow == 0) {
        sprintf(errtxt, "%s: mftarr *tabdat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check number of data points */
    for (j = 0; j <= next; j++) {
        if (cpl_vector_get_size(vecdat->vec[j]) != nrow) {
                sprintf(errtxt, "%s: mfvarr vecdat->vec[%d] != "
                        "mftarr tabdat->tab[1]", MF_ERROR_IDG_TXT, j);
                return cpl_error_set_message(cpl_func, MF_ERROR_IDG, "%s",
                                             errtxt);
        }
    }

    /* Get the two (extreme) original mask values */
    if (mf_conv_getmaskval(maskval, tabdat, extnames) != CPL_ERROR_NONE) {
        return cpl_error_get_code();
    }

    /* Initialise mfvarr structure for transmission curve */
    mf_conv_varr_init(transdat, 0);

    /* Duplicate header and modify BUNIT keyword for transmission curve */
    transdat->head[0] = cpl_propertylist_duplicate(vecdat->head[0]);
    if (cpl_propertylist_has(transdat->head[0], "BUNIT") == 1) {
        cpl_propertylist_set_string(transdat->head[0], "BUNIT",
                                    "Transmission");
    }

    /* Duplicate flux vector for transmission curve */
    transdat->vec[0] = cpl_vector_duplicate(vecdat->vec[extn_flux]);

    /* Get pointers to image data */
    trans = cpl_vector_get_data(transdat->vec[0]);
    flux = cpl_vector_get_data(vecdat->vec[extn_flux]);
    if (extn_dflux >= 0) {
        dflux = cpl_vector_get_data(vecdat->vec[extn_dflux]);
    }
    if (extn_mask >= 0) {
        mask = cpl_vector_get_data(vecdat->vec[extn_mask]);
    }

    /* Get pointers to columns of results data table */
    mtrans = cpl_table_get_data_double(tabdat->tab[1], "mtrans");
    if (mtrans == NULL) hascol = CPL_FALSE;
    cflux = cpl_table_get_data_double(tabdat->tab[1], "tacflux");
    if (cflux == NULL) hascol = CPL_FALSE;
    if (extn_dflux >= 0) {
        cdflux = cpl_table_get_data_double(tabdat->tab[1], "tacdflux");
        if (cdflux == NULL) hascol = CPL_FALSE;
    }
    qual = cpl_table_get_data_int(tabdat->tab[1], "tacqual");
    if (qual == NULL) hascol = CPL_FALSE;
    if (hascol == CPL_FALSE) {
        sprintf(errtxt, "%s: mftarr *tabdat (required results column(s) not "
                "present)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Write transmission, corrected flux, and mask into image vectors */
    for (i = 0; i < nrow; i++) {
        trans[i] = mtrans[i];
        flux[i] = cflux[i];
        if (extn_dflux >= 0) {
            dflux[i] = cdflux[i];
        }
        if (extn_mask >= 0) {
            /* Correct mask value if necessary (good -> bad) */
            if (mask[i] == maskval[1] && qual[i] == 0) {
                mask[i] = maskval[0];
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_getmaskval(double maskval[2], const mftarr *tabdat,
                                  const cpl_table *extnames)
{
    /*!
     * Gets the bad (maskval[0]) and good (maskval[1]) mask values from an
     * :mftarr structure containing the data of the FITS table prepared for
     * MOLECFIT.
     *
     * \b INPUT:
     * \param tabdat    ::mftarr structure with input data
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b OUTPUT:
     * \param maskval   array for bad and good mask value
     *
     * \b ERRORS:
     * - Invalid object structure
     * - No data
     */

    char errtxt[MF_MAXLEN] = "", col_mask[MF_LENLINE+2] = "";
    char col_imask[MF_LENLINE+2] = "";
    cpl_boolean hascol = CPL_TRUE;
    cpl_size imaskmaxpos = 0;
    int extn_mask = 0;
    double maskmin = 0., maskmax = 0.;

    /* Default mask values */
    maskval[0] = 0;
    maskval[1] = 1;

    /* Check table for extension numbers and names */
    if (cpl_table_get_nrow(extnames) != 4) {
        sprintf(errtxt, "%s: cpl_table *extnames (number of rows != 4)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check column names */
    if (cpl_table_has_column(extnames, "col") != 1 ||
        cpl_table_has_column(extnames, "extn") != 1) {
        sprintf(errtxt, "%s: cpl_table *extnames (column names)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get name and extension of input mask column */
    sprintf(col_mask, "%s", cpl_table_get_string(extnames, "col", 3));
    extn_mask = cpl_table_get_int(extnames, "extn", 3, NULL);

    /* Get name of integer mask column */
    if (extn_mask >= 0) {
        sprintf(col_imask, "%s_I", col_mask);
    } else {
        sprintf(col_imask, "%s_I", MF_DEFMASKCOL);
    }

    /* Check number of extensions in results table */
    if (tabdat->next != 1) {
        sprintf(errtxt, "%s: mftarr *tabdat (number of extensions != 1)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check number of data points */
    if (cpl_table_get_nrow(tabdat->tab[1]) == 0) {
        sprintf(errtxt, "%s: mftarr *tabdat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check existence of mask columns in input data table */
    if (extn_mask >= 0) {
        if (cpl_table_has_column(tabdat->tab[1], col_mask) != 1) {
            hascol = CPL_FALSE;
        }
    }
    if (cpl_table_has_column(tabdat->tab[1], col_imask) != 1) {
        hascol = CPL_FALSE;
    }
    if (hascol == CPL_FALSE) {
        sprintf(errtxt, "%s: mftarr *tabdat (required input column(s) not "
                "present)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get the two (extreme) original mask values */
    if (extn_mask >= 0) {
        maskmin = cpl_table_get_column_min(tabdat->tab[1], col_mask);
        maskmax = cpl_table_get_column_max(tabdat->tab[1], col_mask);
        cpl_table_get_column_maxpos(tabdat->tab[1], col_imask, &imaskmaxpos);
        if (cpl_table_get_int(tabdat->tab[1], col_imask, imaskmaxpos, NULL)
            == 1) {
            maskval[0] = maskmax;
            maskval[1] = maskmin;
        } else {
            maskval[0] = maskmin;
            maskval[1] = maskmax;
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_iarr2varr(mfvarr *vecdat, const mfiarr *imadat,
                                 const cpl_table *extnames)
{
    /*!
     * Fills CPL vectors and CPL property lists of a ::mfvarr structure by
     * data read from a 2D FITS image with optional extensions for error,
     * mask, and wavelength, and added in spatial direction. The image data is
     * provided by a ::mfiarr structure. The extension numbers and names are
     * given by a CPL table. Extensions are not read if names are missing or
     * wrong. Images for different chips have to be converted separately.
     *
     * The flux summation considers only good pixels if a mask extension is
     * provided. If pixels along the spatial direction are not considered, the
     * summed flux is corrected by means of a model of the spatial flux
     * distribution. The latter is derived from a projection of the 2D
     * spectrum. If errors exist, they are added quadratically. Possible
     * correlations of the pixels are not considered.
     *
     * \b INPUT:
     * \param imadat    ::mfiarr structure with data of 2D FITS image
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b OUTPUT:
     * \param vecdat    ::mfvarr structure with summed data of 2D FITS image
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     * - Invalid object value(s)
     * - see subroutines
     */

    cpl_vector *row, *prof, *relsum;
    char errtxt[MF_MAXLEN] = "";
    int next = 0, h = 0, nx = 0, ny = 0, j = 0, i = 0, qual = 0;
    const int *extn;
    double maskval[2] = {0., 0.}, profmin = 0., sum = 0., scale = 0.;
    double mask = 0., val = 0.;

    /* Get number of extensions */
    next = imadat->next;

    /* Initialise CPL vectors and CPL property lists for content of mfiarr
       structure (memory allocation) */
    mf_conv_varr_init(vecdat, next);
    for (h = 0; h <= next; h++) {
        vecdat->vec[h] = cpl_vector_new(1);
    }

    /* Transfer FITS header data */
    for (h = 0; h <= next; h++) {
        vecdat->head[h] = cpl_propertylist_duplicate(imadat->head[h]);
    }

    /* Check table for extension numbers and names */
    if (cpl_table_get_nrow(extnames) != 4) {
        mf_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (number of rows != 4)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check column names */
    if (cpl_table_has_column(extnames, "col") != 1 ||
        cpl_table_has_column(extnames, "extn") != 1) {
        mf_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (column names)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of data */
    if (cpl_table_get_column_max(extnames, "extn") == -1) {
        mf_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (no extension with valid "
                "data)", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Check extension numbers */
    if (cpl_table_get_column_max(extnames, "extn") > next) {
        mf_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: cpl_table *extnames (incorrect extension "
                "number)", MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Get pointer to extension numbers */
    extn = cpl_table_get_data_int_const(extnames, "extn");

    /* Check image size */
    nx = cpl_image_get_size_x(imadat->ima[0]);
    ny = cpl_image_get_size_y(imadat->ima[0]);
    if (nx <= 0 || ny <= 0) {
        mf_conv_varr_delete(vecdat);
        sprintf(errtxt, "%s: mfiarr *imadat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    } else if (nx > 0 && ny > 0 && next > 0) {
        for (h = 1; h <= next; h++) {
            if (cpl_image_get_size_x(imadat->ima[h]) != nx ||
                cpl_image_get_size_y(imadat->ima[h]) != ny) {
                mf_conv_varr_delete(vecdat);
                sprintf(errtxt, "%s: mfiarr *imadat (image size differs for "
                        "different extensions", MF_ERROR_IOS_TXT);
                return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                             errtxt);
            }
        }
    }

    /* Get meaning of mask values (identification of good and bad pixels) */
    if (extn[3] >= 0) {
        if (mf_basic_getmaskval_image(maskval, imadat->ima[extn[3]]) !=
            CPL_ERROR_NONE) {
            return MF_ERROR_IOV;
        }
    }

    /* Create row (x) and profile (y) vector */
    row = cpl_vector_new(nx);
    prof = cpl_vector_new(ny);

    /* Get profile along slit */
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            cpl_vector_set(row, i,
                           cpl_image_get(imadat->ima[extn[1]], i+1, j+1,
                                         &qual));
        }
        cpl_vector_set(prof, j, cpl_vector_get_median(row));
    }

    /* Delete temporary row vector */
    cpl_vector_delete(row);

    /* Avoid negative profile values */
    profmin = cpl_vector_get_min(prof);
    if (profmin < 0.) {
        cpl_vector_subtract_scalar(prof, profmin);
    }

    /* Normalise profile function to get weights */
    for (sum = 0., j = 0; j < ny; j++) {
        sum += cpl_vector_get(prof, j);
    }
    if (sum == 0.) {
        mf_conv_varr_delete(vecdat);
        cpl_vector_delete(prof);
        sprintf(errtxt, "%s: mfiarr *imadat (flux sum = 0)",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }
    cpl_vector_divide_scalar(prof, sum);

    /* Calculate correction function for columns with bad pixels */
    relsum = cpl_vector_new(nx);
    if (extn[3] >= 0) {
        for (i = 0; i < nx; i++) {
            for (sum = 0., j = 0; j < ny; j++) {
                if (cpl_image_get(imadat->ima[extn[3]], i+1, j+1, &qual) <
                    maskval[1] + MF_TOL) {
                    sum += cpl_vector_get(prof, j);
                }
            }
            cpl_vector_set(relsum, i, sum);
        }
    } else {
        cpl_vector_fill(relsum, 1.);
    }

    /* Delete temporary profile vector */
    cpl_vector_delete(prof);

    /* Set size of vectors */
    for (h = 0; h <= next; h++) {
        cpl_vector_set_size(vecdat->vec[h], nx);
    }

    /* Convert mfiarr to mfvarr structure (unidentified extensions are
       handled like flux) */

    for (h = 0; h <= next; h++) {
        for (i = 0; i < nx; i++) {

            /* Get scaling factor for pixel flux */
            scale = cpl_vector_get(relsum, i);

            if (h == extn[3]) {

                /* Set mask values (bad value only if no good pixel in image
                   column) */
                if (scale < 2 * MF_TOL) {
                    cpl_vector_set(vecdat->vec[extn[3]], i, maskval[0]);
                } else {
                    cpl_vector_set(vecdat->vec[extn[3]], i, maskval[1]);
                }

            } else {

                /* Add values of good pixels */

                for (sum = 0., j = 0; j < ny; j++) {

                    /* Get mask value for pixel */
                    if (extn[3] >= 0) {
                        mask = cpl_image_get(imadat->ima[extn[3]], i+1, j+1,
                                             &qual);
                    } else {
                        mask = maskval[1];
                    }

                    /* Get pixel value */
                    val = cpl_image_get(imadat->ima[h], i+1, j+1, &qual);
                    if (val <= 0. && h == extn[0]) {
                        mf_conv_varr_delete(vecdat);
                        cpl_vector_delete(relsum);
                        sprintf(errtxt, "%s: mfiarr *imadat "
                                "(wavelength <= 0)", MF_ERROR_IOV_TXT);
                        return cpl_error_set_message(cpl_func, MF_ERROR_IOV,
                                                     "%s", errtxt);
                    }

                    /* Consider good pixels only if not wavelength */
                    if (mask == maskval[1] || h == extn[0]) {
                        /* Add pixel values */
                        if (h == extn[2]) {
                            /* Squared summation of error pixels */
                            sum += val * val;
                        } else {
                            sum += val;
                        }
                    }

                }

                /* Get resulting value depending on kind of data */
                if (h == extn[0]) {
                    /* Wavelength */
                    sum /= (double) ny;
                } else if (h == extn[2]) {
                    /* Error */
                    if (scale > 0.) {
                        sum = sqrt(sum) / scale;
                    }
                } else {
                    /* Flux */
                    if (scale > 0.) {
                        sum /= scale;
                    }
                }

                /* Write resulting value into output vector */
                cpl_vector_set(vecdat->vec[h], i, sum);
                cpl_vector_get(vecdat->vec[h], i);

            }

        }
    }

    /* Free allocated memory */
    cpl_vector_delete(relsum);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_ascii_read(mftarr *tabdat, const char *filename,
                                  const cpl_array *colnames)
{
    /*!
     * Fills CPL table of a ::mftarr structure by data read from an ASCII
     * file. The number and names of the expected columns are provided by a
     * CPL array. Header lines in the ASCII file are allowed if they are
     * marked by '#'.
     *
     * \b INPUT:
     * \param filename  path and name of input ASCII file
     * \param colnames  CPL array of column names for reading of ASCII file
     *
     * \b OUTPUT:
     * \param tabdat    ::mftarr structure with data of ASCII file
     *
     * \b ERRORS:
     * - No data
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    mfpar val[MF_MAXPAR];
    char errtxt[MF_MAXLEN] = "", str[MF_LENLINE+2] = "";
    const char **col;
    int ncolmin = 0, nrec = 0, j = 0, i = 0;
    int ncol0 = MF_MAXPAR, ncol = MF_MAXPAR;

    /* Get size of column name array */
    ncolmin = cpl_array_get_size(colnames);
    if (ncolmin == 0 || colnames == NULL) {
        sprintf(errtxt, "%s: cpl_array *colnames", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Get pointer for column name array */
    col = cpl_array_get_data_string_const(colnames);

    /* Check file existence and open ASCII file */
    if ((stream = fopen(filename, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Find number of data lines */
    while (fgets(str, MF_LENLINE+2, stream) != NULL) {
        if (str[0] == '#' || str[0] == '\n') {
        } else if (isdigit(str[0]) || isspace(str[0]) || str[0] == '-') {
            nrec++;
        } else {
            fclose(stream);
            sprintf(errtxt, "%s: %s (unexpected first character at line)",
                    MF_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }
    }
    rewind(stream);

    /* No data points */
    if (nrec == 0) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (no data)",
                MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Number of values per line */
    mf_basic_readline(stream, val, &ncol0);
    if (ncol0 < ncolmin) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (too low number of columns)",
                MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }
    rewind(stream);

    /* Initialise CPL table and CPL property list for content of ASCII file
       (memory allocation) */
    mf_conv_tarr_init(tabdat, 1);

    /* Create empty CPL property lists */
    tabdat->head[0] = cpl_propertylist_new();
    tabdat->head[1] = cpl_propertylist_new();

    /* Create CPL tables of required size (put data in first extension of
       output file) */
    tabdat->tab[0] = cpl_table_new(0);
    tabdat->tab[1] = cpl_table_new(nrec);

    /* Create required table columns */
    for (j = 0; j < ncolmin; j++) {
        cpl_table_new_column(tabdat->tab[1], col[j], CPL_TYPE_DOUBLE);
    }

    /* Read spectral data from file and write it to CPL table */

    for (i = 0; i < nrec; i++) {

        mf_basic_readline(stream, val, &ncol);

        if (ncol != ncol0) {
            mf_conv_tarr_delete(tabdat);
            fclose(stream);
            sprintf(errtxt, "%s: %s (unexpected number of values at line)",
                    MF_ERROR_UFS_TXT, filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }

        for (j = 0; j < ncolmin; j++) {
            cpl_table_set(tabdat->tab[1], col[j], i, val[j].d);
        }

    }

    /* Close ASCII file */
    fclose(stream);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_ascii_write(const char *filename, const mftarr *tabdat)
{
    /*!
     * Writes the table data of an ::mftarr structure into an ASCII file. Only
     * one extension is allowed. The column names are written into a header
     * line starting with '#'. The routine supports the column types STRING,
     * INT, FLOAT, and DOUBLE.
     *
     * \b INPUT:
     * \param filename  path and name of output ASCII file
     * \param tabdat    ::mftarr structure with tabulated data
     *
     * \b ERRORS:
     * - Invalid object structure
     * - No data
     * - File opening failed
     * - see subroutines
     */

    FILE *stream;
    cpl_type coltype;
    cpl_array *colnames = NULL;
    char errtxt[MF_MAXLEN] = "", str[MF_MAXLEN] = "";
    char strcomp[MF_LENLINE+2] = "";
    char **col;
    int nrow = 0, ncol = 0, i = 0, j = 0;

    /* Check number of extensions */
    if (tabdat->next != 1) {
        sprintf(errtxt, "%s: mftarr *tabdat (number of extensions != 1)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get and check number of data points */
    nrow = cpl_table_get_nrow(tabdat->tab[1]);
    if (nrow == 0) {
        sprintf(errtxt, "%s: mftarr *tabdat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Get column names in input table */
    colnames = cpl_table_get_column_names(tabdat->tab[1]);
    ncol = cpl_array_get_size(colnames);

    /* Get pointer to array */
    col = cpl_array_get_data_string(colnames);

    /* Open output ASCII file */
    if ((stream = fopen(filename, "w+")) == NULL) {
        cpl_array_delete(colnames);
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Build string of column names */
    sprintf(str, "#");
    for (j = 0; j < ncol; j++) {
        sprintf(strcomp, " %s", col[j]);
        strncat(str, strcomp, strlen(strcomp));
    }

    /* Write column names into ASCII file */
    fprintf(stream, "%s\n", str);

    /* Write data into ASCII file */

    for (i = 0; i < nrow; i++) {

        /* Build string of data values */
        str[0] = '\0';
        for (j = 0; j < ncol; j++) {
            coltype = cpl_table_get_column_type(tabdat->tab[1], col[j]);
            if (coltype == CPL_TYPE_STRING) {
                sprintf(strcomp, "%s", cpl_table_get_string(tabdat->tab[1],
                                                            col[j], i));
            } else if (coltype == CPL_TYPE_INT) {
                sprintf(strcomp, "%d", cpl_table_get_int(tabdat->tab[1],
                                                         col[j], i, NULL));
            } else if (coltype == CPL_TYPE_FLOAT) {
                sprintf(strcomp, "%f", cpl_table_get_float(tabdat->tab[1],
                                                           col[j], i, NULL));
            } else if (coltype == CPL_TYPE_DOUBLE) {
                sprintf(strcomp, "%e", cpl_table_get_double(tabdat->tab[1],
                                                            col[j], i, NULL));
            } else {
                fclose(stream);
                cpl_array_delete(colnames);
                sprintf(errtxt, "%s: mftarr tabdat->tab[1] (unsupported "
                        "column type)", MF_ERROR_IOS_TXT);
                return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                             errtxt);
            }
            strncat(str, strcomp, strlen(strcomp));
            if (j != ncol-1) {
                sprintf(strcomp, " ");
                strncat(str, strcomp, strlen(strcomp));
            }
        }

        /* Write data string into ASCII file */
        fprintf(stream, "%s\n", str);

   }

    /* Close ASCII file */
    fclose(stream);

    /* Free allocated memory */
    cpl_array_delete(colnames);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_tarr_init(mftarr *tabdat, const int next)
{
    /*!
     * Initialises array of CPL tables and CPL property lists as ::mftarr
     * structure depending on the input number of FITS extensions.
     *
     * \b INPUT:
     * \param next    number of FITS extensions
     *
     * \b OUTPUT:
     * \param tabdat  empty ::mftarr structure for FITS table data
     *
     * \b ERRORS:
     * - none
     */

    tabdat->next = next;
    tabdat->head = cpl_calloc(next+1, sizeof(cpl_propertylist *));
    tabdat->tab = cpl_calloc(next+1, sizeof(cpl_table *));

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_tarr_read(mftarr *tabdat, const char *filename)
{
    /*!
     * Reads a FITS table with an arbritrary number of extensions and puts the
     * data into an ::mftarr structure consisting of an array of CPL tables
     * and CPL property lists.
     *
     * \b INPUT:
     * \param filename  path and name of input FITS table
     *
     * \b OUTPUT:
     * \param tabdat    ::mftarr structure with FITS table data
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    char errtxt[MF_MAXLEN] = "";
    int next = 0, fitsformat = 0, i = 0;

    /* Get number of extensions and check file */
    next = cpl_fits_count_extensions(filename);
    if (next < 0) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    } else if (next == 0) {
        sprintf(errtxt, "%s: %s (number of FITS extensions = 0)",
                MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Check FITS format */
    mf_conv_checkfitsformat(&fitsformat, filename);
    if (fitsformat != 1) {
        sprintf(errtxt, "%s: %s (no FITS table)", MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Initialise CPL tables and CPL property lists for content of FITS
       file (memory allocation) */
    mf_conv_tarr_init(tabdat, next);

    /* Read FITS file */
    for (i = 0; i <= next; i++) {
        tabdat->head[i] = cpl_propertylist_load(filename, i);
        if (i == 0) {
            tabdat->tab[i] = cpl_table_new(0);
        } else {
            tabdat->tab[i] = cpl_table_load(filename, i, 0);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_tarr_write(const char *filename, const mftarr *tabdat)
{
    /*!
     * Writes data of an ::mftarr structure into a FITS table. The number of
     * FITS extensions depends on the content of the ::mftarr structure.
     *
     * \b INPUT:
     * \param filename  path and name of output FITS table
     * \param tabdat    ::mftarr structure with FITS table data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    /* Get number of FITS extensions */
    next = tabdat->next;

    /* Write FITS table with next extensions */
    for (i = 0; i < next; i++) {
        if (i == 0) {
            cpl_table_save(tabdat->tab[i+1], tabdat->head[0],
                           tabdat->head[i+1], filename, CPL_IO_CREATE);
        } else {
            cpl_table_save(tabdat->tab[i+1], NULL,
                           tabdat->head[i+1], filename, CPL_IO_EXTEND);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_tarr_delete(mftarr *tabdat)
{
    /*!
     * Deletes an ::mftarr structure, which contains an array of CPL tables
     * and CPL property lists.
     *
     * \b INPUT:
     * \param tabdat  ::mftarr structure with FITS table data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    if (tabdat == NULL) {
        return CPL_ERROR_NONE;
    }

    next = tabdat->next;

    for (i = 0; i <= next; i++) {
        cpl_propertylist_delete(tabdat->head[i]);
        cpl_table_delete(tabdat->tab[i]);
    }

    cpl_free(tabdat->head);
    cpl_free(tabdat->tab);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_varr_init(mfvarr *vecdat, const int next)
{
    /*!
     * Initialises array of CPL vectors and CPL property lists as ::mfvarr
     * structure depending on the input number of FITS extensions.
     *
     * \b INPUT:
     * \param next    number of FITS extensions
     *
     * \b OUTPUT:
     * \param vecdat  empty ::mfvarr structure for 1D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    vecdat->next = next;
    vecdat->head = cpl_calloc(next+1, sizeof(cpl_propertylist *));
    vecdat->vec = cpl_calloc(next+1, sizeof(cpl_vector *));

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_varr_read(mfvarr *vecdat, const char *filename)
{
    /*!
     * Reads a 1D FITS image with an arbritrary number of extensions and puts
     * the data into an ::mfvarr structure consisting of an array of CPL
     * vectors and CPL property lists.
     *
     * \b INPUT:
     * \param filename  path and name of input 1D FITS image
     *
     * \b OUTPUT:
     * \param vecdat    ::mfvarr structure with 1D FITS image data
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    char errtxt[MF_MAXLEN] = "";
    int next = 0, fitsformat = 0, i = 0;

    /* Get number of extensions and check file */
    next = cpl_fits_count_extensions(filename);
    if (next < 0) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Check FITS format */
    mf_conv_checkfitsformat(&fitsformat, filename);
    if (fitsformat != 2) {
        sprintf(errtxt, "%s: %s (no 1D FITS image)", MF_ERROR_UFS_TXT,
                filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Initialise CPL tables and CPL property lists for content of FITS
       file (memory allocation) */
    mf_conv_varr_init(vecdat, next);

    /* Read FITS file */
    for (i = 0; i <= next; i++) {
        vecdat->head[i] = cpl_propertylist_load(filename, i);
        vecdat->vec[i] = cpl_vector_load(filename, i);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_varr_write(const char *filename, const mfvarr *vecdat)
{
    /*!
     * Writes data of an ::mfvarr structure into a 1D FITS image. The number
     * of FITS extensions depends on the content of the ::mfvarr structure.
     *
     * \b INPUT:
     * \param filename  path and name of output 1D FITS image
     * \param vecdat    ::mfvarr structure with 1D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    /* Get number of FITS extensions */
    next = vecdat->next;

    /* Write 1D FITS image with next extensions */
    for (i = 0; i <= next; i++) {
        if (i == 0) {
            cpl_vector_save(vecdat->vec[i], filename, CPL_TYPE_FLOAT,
                            vecdat->head[i], CPL_IO_CREATE);
        } else {
            cpl_vector_save(vecdat->vec[i], filename, CPL_TYPE_FLOAT,
                            vecdat->head[i], CPL_IO_EXTEND);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_varr_delete(mfvarr *vecdat)
{
    /*!
     * Deletes an ::mfvarr structure, which contains an array of CPL vectors
     * and CPL property lists.
     *
     * \b INPUT:
     * \param vecdat  ::mfvarr structure with 1D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    if (vecdat == NULL) {
        return CPL_ERROR_NONE;
    }

    next = vecdat->next;

    for (i = 0; i <= next; i++) {
        cpl_propertylist_delete(vecdat->head[i]);
        cpl_vector_delete(vecdat->vec[i]);
    }

    cpl_free(vecdat->head);
    cpl_free(vecdat->vec);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_iarr_init(mfiarr *imadat, const int next)
{
    /*!
     * Initialises array of CPL images and CPL property lists as ::mfiarr
     * structure depending on the input number of FITS extensions.
     *
     * \b INPUT:
     * \param next    number of FITS extensions
     *
     * \b OUTPUT:
     * \param imadat  empty ::mfiarr structure for 2D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    imadat->next = next;
    imadat->head = cpl_calloc(next+1, sizeof(cpl_propertylist *));
    imadat->ima = cpl_calloc(next+1, sizeof(cpl_image *));

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_iarr_read(mfiarr *imadat, const char *filename)
{
    /*!
     * Reads a 2D FITS image with an arbritrary number of extensions and puts
     * the data into an ::mfiarr structure consisting of an array of CPL
     * images and CPL property lists.
     *
     * \b INPUT:
     * \param filename  path and name of input 2D FITS image
     *
     * \b OUTPUT:
     * \param imadat    ::mfiarr structure with 2D FITS image data
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     */

    char errtxt[MF_MAXLEN] = "";
    int next = 0, fitsformat = 0, i = 0;

    /* Get number of extensions and check file */
    next = cpl_fits_count_extensions(filename);
    if (next < 0) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Check FITS format */
    mf_conv_checkfitsformat(&fitsformat, filename);
    if (fitsformat != 3) {
        sprintf(errtxt, "%s: %s (no 2D FITS image)", MF_ERROR_UFS_TXT,
                filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Initialise CPL tables and CPL property lists for content of FITS
       file (memory allocation) */
    mf_conv_iarr_init(imadat, next);

    /* Read FITS file */
    for (i = 0; i <= next; i++) {
        imadat->head[i] = cpl_propertylist_load(filename, i);
        imadat->ima[i] = cpl_image_load(filename, CPL_TYPE_UNSPECIFIED, 0, i);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_iarr_write(const char *filename, const mfiarr *imadat)
{
    /*!
     * Writes data of an ::mfiarr structure into a 2D FITS image. The number
     * of FITS extensions depends on the content of the ::mfiarr structure.
     *
     * \b INPUT:
     * \param filename  path and name of output 2D FITS image
     * \param imadat    ::mfiarr structure with 2D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    /* Get number of FITS extensions */
    next = imadat->next;

    /* Write 2D FITS image with next extensions */
    for (i = 0; i <= next; i++) {
        if (i == 0) {
            cpl_image_save(imadat->ima[i], filename, CPL_TYPE_UNSPECIFIED,
                           imadat->head[i], CPL_IO_CREATE);
        } else {
            cpl_image_save(imadat->ima[i], filename, CPL_TYPE_UNSPECIFIED,
                           imadat->head[i], CPL_IO_EXTEND);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_conv_iarr_delete(mfiarr *imadat)
{
    /*!
     * Deletes an ::mfiarr structure, which contains an array of CPL images
     * and CPL property lists.
     *
     * \b INPUT:
     * \param imadat  ::mfiarr structure with 2D FITS image data
     *
     * \b ERRORS:
     * - none
     */

    int next = 0, i = 0;

    if (imadat == NULL) {
        return CPL_ERROR_NONE;
    }

    next = imadat->next;

    for (i = 0; i <= next; i++) {
        cpl_propertylist_delete(imadat->head[i]);
        cpl_image_delete(imadat->ima[i]);
    }

    cpl_free(imadat->head);
    cpl_free(imadat->ima);

    return CPL_ERROR_NONE;
}

/**@}*/
