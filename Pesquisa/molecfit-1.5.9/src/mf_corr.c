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
 * \file mf_corr.c
 *
 * Routines for telluric absorption correction of a set of files by means of
 * a provided transmission curve
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  29 Aug 2012
 * \date   09 Apr 2015
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/
#include <mf_molecfit.h>
#include <mf_corr.h>
#include <string.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code mf_corr_filelist(const char *parfile_in)
{
    /*!
     * \callgraph
     *
     * Corrects a list of spectra provided by an ASCII file for telluric
     * absorption by means of a model transmission curve derived from a
     * reference spectrum. The file type, i.e. ASCII file, FITS table, or FITS
     * image, of the input spectra has to be same as for the reference
     * spectrum. The list of FITS images can contain 1D and 2D images. The
     * corresponding transmission curve has to be a 1D FITS image.
     *
     * \b INPUT:
     * \param parfile_in  name of parameter file
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_array *filelist = NULL;
    mfdrv drvpar;
    char parfile[MF_MAXLEN] = "", basedir[MF_MAXLEN] = "";
    char outdir[MF_MAXLEN] = "", filename[MF_MAXLEN] = "";
    char infilename[MF_MAXLEN] = "", suffix[MF_MAXLEN] = "";
    char reffilename[MF_MAXLEN] = "", transfilename[MF_MAXLEN] = "";
    char errtxt[MF_MAXLEN] = "";
    int fitsformat = 0;

    /* Read MOLECFIT driver file */
    mf_par_initall(&drvpar);
    mf_basic_absfile(parfile, parfile_in);
    if ((status = mf_par_readfile(&drvpar, parfile, 0)) != CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        return status;
    }

    // Modification to stop using system directories in case of a
    // systenwide installation
    char* tmpdir = NULL;
    if ((status = fix_directories(&drvpar, &tmpdir)) != CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        return status;
    }

    /* Get file names for reference transmission curve */
    p = cpl_parameterlist_find(drvpar.parlist, "basedir");
    strncpy(basedir, cpl_parameter_get_string(p), MF_MAXLEN);
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar.parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar.parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar.parlist, "filename");
    strncpy(filename, cpl_parameter_get_string(p), MF_MAXLEN);
    mf_basic_getfilename(NULL, infilename, suffix, filename);
    sprintf(reffilename, "%s%s_TAC.%s", outdir, infilename, suffix);
    sprintf(transfilename, "%s%s_TRA.%s", outdir, infilename, suffix);

    /* Get file type */
    if ((status = mf_conv_checkfitsformat(&fitsformat, reffilename)) !=
        CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        return status;
    }

    /* Get array of file names for correction */
    filelist = cpl_array_new(0, CPL_TYPE_STRING);
    if ((status = mf_corr_readfilelist(filelist, &drvpar)) !=
        CPL_ERROR_NONE) {
        mf_par_deleteall(&drvpar);
        return status;
    }

    /* Check existence of file list */
    if (cpl_array_get_size(filelist) == 0) {
        mf_par_deleteall(&drvpar);
        cpl_array_delete(filelist);
        cpl_msg_info(cpl_func, "No file list -> exit");
        return CPL_ERROR_NONE;
    }

    /* Write info message */
    cpl_msg_info(cpl_func, "Read reference file %s", reffilename);

    /* Correct files dependent on type of reference file */
    if (fitsformat == 0) {
        /* Correct non-FITS files (ASCII format assumed) */
        mf_corr_ascii(filelist, reffilename, &drvpar);
    } else if (fitsformat == 1) {
        /* Correct FITS tables */
        mf_corr_fitstable(filelist, reffilename, &drvpar);
    } else if (fitsformat == 2) {
        /* Correct 1D and 2D FITS images */
        cpl_msg_info(cpl_func, "Read transmission file %s", transfilename);
        mf_corr_fitsimage(filelist, reffilename, transfilename, &drvpar);
    } else if (fitsformat > 2) {
        mf_par_deleteall(&drvpar);
        cpl_array_delete(filelist);
        sprintf(errtxt, "%s: %s (multidimensional FITS image)",
                MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Free allocated memory */
    mf_par_deleteall(&drvpar);
    cpl_array_delete(filelist);

    mf_cleanup_standalone(tmpdir);

    /* Return error code of last error */
    return cpl_error_get_code();
}


cpl_error_code mf_corr_readfilelist(cpl_array *filelist, const mfdrv *drvpar)
{
    /*!
     * Reads a list of file names from an ASCII file and writes them into a
     * CPL array. The name of the input ASCII file has to be provided by the
     * MOLECFIT driver file. If the name of the ASCII file is "none", an empty
     * CPL array is returned.
     *
     * \b INPUT:
     * \param filelist  empty CPL string array
     * \param drvpar    ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param filelist  array of file names
     *
     * \b ERRORS:
     * - Invalid object structure
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    cpl_parameter *p;
    char errtxt[MF_MAXLEN] = "", basedir[MF_MAXLEN] = "";
    char listname[MF_MAXLEN] = "", rellistname[MF_MAXLEN] = "";
    char filename[MF_MAXLEN] = "", relfilename[MF_MAXLEN] = "";
    int next = 0, i = 0;

    /* Check correct format of CPL array */
    if (cpl_array_get_type(filelist) != CPL_TYPE_STRING) {
        sprintf(errtxt, "%s: cpl_array *filelist (not a string array)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get name of file list */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strncpy(basedir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "listname");
    strncpy(listname, cpl_parameter_get_string(p), MF_MAXLEN);
    if (strcmp(listname, "none") == 0) {
        /* No file list -> return empty CPL array */
        cpl_array_set_size(filelist, 0);
        return CPL_ERROR_NONE;
    } else if (listname[0] != '/') {
        char curdir[MF_MAXLEN];
        p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        sprintf(rellistname, "%s/%s", curdir, listname);
        mf_basic_absfile(listname, rellistname);
    }

    /* Write info message */
    cpl_msg_info(cpl_func, "Read file list %s", listname);

    /* Check file existence */
    if ((stream = fopen(listname, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, listname);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Check whether a FITS file was given erroneously */
    next = cpl_fits_count_extensions(listname);
    cpl_errorstate_set(CPL_ERROR_NONE);
    if (next != -1) {
        fclose(stream);
        sprintf(errtxt, "%s: %s (FITS file instead of ASCII list of file "
                "names)", MF_ERROR_UFS_TXT, listname);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Read file names from file and write them into a CPL array */

    while (fgets(filename, MF_MAXLEN, stream) != NULL) {
        /* Remove leading and trailing blanks from string */
        mf_basic_strtrim_inplace(filename);
        /* skip empty lines */
        if (strlen(filename) == 0) {
            continue;
        }

        /* Count files */
        i++;

        /* Resize array */
        cpl_array_set_size(filelist, i);

        /* Convert path if not absolute */
        if (filename[0] != '/') {
            char curdir[MF_MAXLEN];
            p = cpl_parameterlist_find(drvpar->parlist, "curdir");
            strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
            sprintf(relfilename, "%s/%s", curdir, filename);
            mf_basic_absfile(filename, relfilename);
        }

        /* Write file name into CPL array */
        cpl_array_set_string(filelist, i-1, filename);

    }

    /* Close ASCII file */
    fclose(stream);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_ascii(const cpl_array *filelist,
                             const char *reffilename, const mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Corrects list of ASCII files provided by a CPL array for telluric
     * absorption. The reference file containing the required transmission
     * curve has to be in ASCII format, too. Moreover, this file has to
     * provide the column names in a header line starting with '#'. The latter
     * is guaranteed by CALCTRANS. The required column is 'mtrans'. The column
     * names of the input files are taken from the parameter file. The columns
     * of these files have to be same as for the input file that was used by
     * MOLECFIT for the calculation of the reference transmission curve.
     *
     * \b INPUT:
     * \param filelist     CPL array of file names
     * \param reffilename  name of reference file
     * \param drvpar       ::mfdrv parameter structure
     *
     * \b ERRORS:
     * - Invalid object structure
     * - No data
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_array *refcolnames, *colnames;
    mftarr reftabdat, tabdat;
    char outdir[MF_MAXLEN] = "", reffluxcol[MF_LENLINE+2] = "";
    char fluxcol[MF_LENLINE+2] = "", dfluxcol[MF_LENLINE+2] = "";
    char errtxt[MF_MAXLEN] = "", infilename[MF_MAXLEN] = "";
    char filename[MF_MAXLEN] = "", suffix[MF_MAXLEN] = "";
    char outfilename[MF_MAXLEN] = "";
    int nfile = 0, i = 0;

    /* Get directory for output files */
    mf_basic_getfilename(outdir, NULL, NULL, reffilename);

    /* Get column names from header line of reference file */
    refcolnames = cpl_array_new(0, CPL_TYPE_STRING);
    if ((status = mf_corr_getcolnames(refcolnames, reffilename)) !=
        CPL_ERROR_NONE) {
        cpl_array_delete(refcolnames);
        return status;
    }

    /* Get name of flux column (second column) of reference file */
    sprintf(reffluxcol, "%s", cpl_array_get_string(refcolnames, 1));

    /* Get column names from parameter file */
    colnames = cpl_array_new(0, CPL_TYPE_STRING);
    mf_conv_setcolnames(colnames, drvpar);

    /* Get name of flux column from parameter file */
    sprintf(fluxcol, "%s", cpl_array_get_string(colnames, 1));

    /* Check agreement of flux column names */
    if (strcmp(fluxcol, reffluxcol) != 0) {
        cpl_array_delete(refcolnames);
        cpl_array_delete(colnames);
        sprintf(errtxt, "%s: cpl_array colnames[1] != "
                "cpl_array refcolnames[1]", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get name of flux error column from parameter file ('NULL' in the case
       of absence) */
    p = cpl_parameterlist_find(drvpar->parlist, "col_dflux");
    strncpy(dfluxcol, cpl_parameter_get_string(p), MF_LENLINE+2);

    /* Read reference file with transmission curve */
    if ((status = mf_conv_ascii_read(&reftabdat, reffilename, refcolnames)) !=
         CPL_ERROR_NONE) {
        cpl_array_delete(refcolnames);
        cpl_array_delete(colnames);
        return status;
    }

    /* Delete array for reference column names */
    cpl_array_delete(refcolnames);

    /* Get number of files for correction */
    nfile = cpl_array_get_size(filelist);
    if (nfile == 0) {
        cpl_array_delete(colnames);
        mf_conv_tarr_delete(&reftabdat);
        sprintf(errtxt, "%s: cpl_array *filelist", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Perform telluric absorption correction for list of ASCII tables */

    for (i = 0; i < nfile; i++) {

        /* Read ASCII file from list */
        sprintf(infilename, "%s", cpl_array_get_string(filelist, i));
        if ((status = mf_conv_ascii_read(&tabdat, infilename, colnames)) !=
            CPL_ERROR_NONE) {
            cpl_array_delete(colnames);
            mf_conv_tarr_delete(&reftabdat);
            return status;
        }

        /* Print message */
        cpl_msg_info(cpl_func, "Correct %s", infilename);

        /* Perform telluric absorption correction */
        if ((status = mf_corr_performtac_tarr(&tabdat, &reftabdat, fluxcol,
                                              dfluxcol)) != CPL_ERROR_NONE) {
            cpl_array_delete(colnames);
            mf_conv_tarr_delete(&reftabdat);
            mf_conv_tarr_delete(&tabdat);
            return status;
        }

        /* Get output file name */
        mf_basic_getfilename(NULL, filename, suffix, infilename);
        sprintf(outfilename, "%s/%s_TAC.%s", outdir, filename, suffix);

        /* Write corrected ASCII file */
        if ((status = mf_conv_ascii_write(outfilename, &tabdat)) !=
            CPL_ERROR_NONE) {
            cpl_array_delete(colnames);
            mf_conv_tarr_delete(&reftabdat);
            mf_conv_tarr_delete(&tabdat);
            return status;
        }

        /* Free allocated memory */
        mf_conv_tarr_delete(&tabdat);

    }

    /* Free allocated memory */
    cpl_array_delete(colnames);
    mf_conv_tarr_delete(&reftabdat);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_fitstable(const cpl_array *filelist,
                                 const char *reffilename, const mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Corrects list of FITS tables provided by a CPL array for telluric
     * absorption. The reference file containing the required transmission
     * curve has to be in FITS table format, too. The name of the flux column
     * (and a possible error column) in the input files is/are taken from the
     * parameter file. The required 'mtrans' column has to be given by the
     * reference file, which has to be produced by MOLECFIT and CALCTRANS from
     * a file with the same input columns as the listed input files that shall
     * be corrected. The input files can have data of multiple chips saved in
     * different FITS extensions.
     *
     * \b INPUT:
     * \param filelist     CPL array of file names
     * \param reffilename  name of reference file
     * \param drvpar       ::mfdrv parameter structure
     *
     * \b ERRORS:
     * - No data
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    mftarr reftabdat, tabdat;
    char outdir[MF_MAXLEN] = "", fluxcolname[MF_LENLINE+2] = "";
    char dfluxcolname[MF_LENLINE+2] = "", errtxt[MF_MAXLEN] = "";
    char infilename[MF_MAXLEN] = "", filename[MF_MAXLEN] = "";
    char suffix[MF_MAXLEN] = "", outfilename[MF_MAXLEN] = "";
    int nfile = 0, i = 0;

    /* Get directory for output files */
    mf_basic_getfilename(outdir, NULL, NULL, reffilename);

    /* Read reference file with transmission curve */
    if ((status = mf_conv_tarr_read(&reftabdat, reffilename)) !=
         CPL_ERROR_NONE) {
        return status;
    }

    /* Get names of flux and error columns (the latter can be 'NULL') */
    p = cpl_parameterlist_find(drvpar->parlist, "col_flux");
    strncpy(fluxcolname, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(fluxcolname, "NULL") == 0) {
        sprintf(fluxcolname, "%s", MF_DEFFLUXCOL);
    }
    p = cpl_parameterlist_find(drvpar->parlist, "col_dflux");
    strncpy(dfluxcolname, cpl_parameter_get_string(p), MF_LENLINE+2);
    if (strcmp(dfluxcolname, "NULL") == 0) {
        sprintf(fluxcolname, "%s", MF_DEFDFLUXCOL);
    }

    /* Get number of files for correction */
    nfile = cpl_array_get_size(filelist);
    if (nfile == 0) {
        mf_conv_tarr_delete(&reftabdat);
        sprintf(errtxt, "%s: cpl_array *filelist", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Perform telluric absorption correction for list of FITS tables */

    for (i = 0; i < nfile; i++) {

        /* Read FITS table from list */
        sprintf(infilename, "%s", cpl_array_get_string(filelist, i));
        if ((status = mf_conv_tarr_read(&tabdat, infilename)) !=
            CPL_ERROR_NONE) {
            mf_conv_tarr_delete(&reftabdat);
            return status;
        }

        /* Print message */
        cpl_msg_info(cpl_func, "Correct %s", infilename);

        /* Perform telluric absorption correction */
        if ((status = mf_corr_performtac_tarr(&tabdat, &reftabdat,
                                              fluxcolname, dfluxcolname)) !=
            CPL_ERROR_NONE) {
            mf_conv_tarr_delete(&reftabdat);
            mf_conv_tarr_delete(&tabdat);
            return status;
        }

        /* Get output file name */
        mf_basic_getfilename(NULL, filename, suffix, infilename);
        sprintf(outfilename, "%s/%s_TAC.%s", outdir, filename, suffix);

        /* Write corrected FITS table */
        if ((status = mf_conv_tarr_write(outfilename, &tabdat)) !=
            CPL_ERROR_NONE) {
            mf_conv_tarr_delete(&reftabdat);
            mf_conv_tarr_delete(&tabdat);
            return status;
        }

        /* Free allocated memory */
        mf_conv_tarr_delete(&tabdat);

    }

    /* Free allocated memory */
    mf_conv_tarr_delete(&reftabdat);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_fitsimage(const cpl_array *filelist,
                                 const char *reffilename,
                                 const char *transfilename,
                                 const mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Corrects list of FITS 1D and/or 2D images provided by a CPL array for
     * telluric absorption. The file containing the required transmission
     * curve has to be a FITS 1D image without an extension as provided by
     * CALCTRANS. The names of the FITS extensions and their numbers in the
     * input files are taken from the parameter and the reference file
     * respectively.
     *
     * \b INPUT:
     * \param filelist       CPL array of file names
     * \param reffilename    name of reference file
     * \param transfilename  name of transmission file
     * \param drvpar         ::mfdrv parameter structure
     *
     * \b ERRORS:
     * - No data
     * - Unexpected file structure
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_table *extnames;
    mfvarr refvecdat, transdat, vecdat;
    mfiarr imadat;
    char outdir[MF_MAXLEN] = "", errtxt[MF_MAXLEN] = "";
    char infilename[MF_MAXLEN] = "", filename[MF_MAXLEN] = "";
    char suffix[MF_MAXLEN] = "", outfilename[MF_MAXLEN] = "";
    int nfile = 0, i = 0, fitsformat = 0;

    /* Get directory for output files */
    mf_basic_getfilename(outdir, NULL, NULL, reffilename);

    /* Read reference file */
    if ((status = mf_conv_varr_read(&refvecdat, reffilename)) !=
         CPL_ERROR_NONE) {
        return status;
    }

    /* Get names and positions of FITS extensions */
    extnames = cpl_table_new(0);
    if ((status = mf_conv_setextnames(extnames, &refvecdat, drvpar)) !=
         CPL_ERROR_NONE) {
        cpl_table_delete(extnames);
        mf_conv_varr_delete(&refvecdat);
        return status;
    }

    /* Delete temporary mfvarr structure */
    mf_conv_varr_delete(&refvecdat);

    /* Read transmission file */
    if ((status = mf_conv_varr_read(&transdat, transfilename)) !=
         CPL_ERROR_NONE) {
        cpl_table_delete(extnames);
        return status;
    }

    /* Get number of files for correction */
    nfile = cpl_array_get_size(filelist);
    if (nfile == 0) {
        cpl_table_delete(extnames);
        mf_conv_varr_delete(&transdat);
        sprintf(errtxt, "%s: cpl_array *filelist", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Perform telluric absorption correction for list of FITS images */

    for (i = 0; i < nfile; i++) {

        /* Get file type */
        sprintf(infilename, "%s", cpl_array_get_string(filelist, i));
        if ((status = mf_conv_checkfitsformat(&fitsformat, infilename)) !=
            CPL_ERROR_NONE) {
            cpl_table_delete(extnames);
            mf_conv_varr_delete(&transdat);
            return status;
        }

        /* Print message */
        cpl_msg_info(cpl_func, "Correct %s", infilename);

        /* Read 1D or 2D FITS image from list */
        if (fitsformat == 2) {
            if ((status = mf_conv_varr_read(&vecdat, infilename)) !=
                CPL_ERROR_NONE) {
                cpl_table_delete(extnames);
                mf_conv_varr_delete(&transdat);
                return status;
            }
        } else if (fitsformat == 3) {
            if ((status = mf_conv_iarr_read(&imadat, infilename)) !=
                CPL_ERROR_NONE) {
                cpl_table_delete(extnames);
                mf_conv_varr_delete(&transdat);
                return status;
            }
        } else {
            cpl_table_delete(extnames);
            mf_conv_varr_delete(&transdat);
            sprintf(errtxt, "%s: %s (no 1D or 2D FITS image)",
                    MF_ERROR_UFS_TXT, infilename);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s",
                                         errtxt);
        }

        /* Convert mfiarr structure (2D) into mfvarr structure (1D) */
        if (fitsformat == 3) {
            if ((status = mf_conv_iarr2varr(&vecdat, &imadat, extnames)) !=
                CPL_ERROR_NONE) {
                cpl_table_delete(extnames);
                mf_conv_varr_delete(&transdat);
                mf_conv_iarr_delete(&imadat);
                return status;
            }
        }

        /* Perform telluric absorption correction */
        if ((status = mf_corr_performtac_varr(&vecdat, &transdat, extnames))
            != CPL_ERROR_NONE) {
            cpl_table_delete(extnames);
            mf_conv_varr_delete(&transdat);
            mf_conv_varr_delete(&vecdat);
            if (fitsformat == 3) {
                mf_conv_iarr_delete(&imadat);
            }
            return status;
        }

        /* Convert mfvarr structure (1D) into mfiarr structure (2D) if
           required */
        if (fitsformat == 3) {
            if ((status = mf_corr_tacvarr2iarr(&imadat, &vecdat, &transdat,
                                               extnames)) != CPL_ERROR_NONE) {
                cpl_table_delete(extnames);
                mf_conv_varr_delete(&transdat);
                mf_conv_varr_delete(&vecdat);
                mf_conv_iarr_delete(&imadat);
                return status;
            }
        }

        /* Get output file name */
        mf_basic_getfilename(NULL, filename, suffix, infilename);
        sprintf(outfilename, "%s/%s_TAC.%s", outdir, filename, suffix);

        /* Write corrected 1D or 2D FITS image */
        if (fitsformat == 2) {
            /* Remove wrong keyword -> bug in XSHOOTER pipeline? */
            if (cpl_propertylist_has(vecdat.head[0], "CTYPE2") == 1) {
                cpl_propertylist_erase(vecdat.head[0], "CTYPE2");
            }
            if ((status = mf_conv_varr_write(outfilename, &vecdat)) !=
                CPL_ERROR_NONE) {
                cpl_table_delete(extnames);
                mf_conv_varr_delete(&transdat);
                mf_conv_varr_delete(&vecdat);
                return status;
            }
        } else if (fitsformat == 3) {
            if ((status = mf_conv_iarr_write(outfilename, &imadat)) !=
                CPL_ERROR_NONE) {
                cpl_table_delete(extnames);
                mf_conv_varr_delete(&transdat);
                mf_conv_varr_delete(&vecdat);
                mf_conv_iarr_delete(&imadat);
                return status;
            }
        }

        /* Free allocated memory */
        mf_conv_varr_delete(&vecdat);
        if (fitsformat == 3) {
            mf_conv_iarr_delete(&imadat);
        }

    }

    /* Free allocated memory */
    cpl_table_delete(extnames);
    mf_conv_varr_delete(&transdat);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_getcolnames(cpl_array *colnames,
                                   const char *filename)
{
    /*!
     * Gets column names of an ASCII file from a header line starting with
     * '#' and put them into a CPL array. The column names have to be
     * separated by spaces.
     *
     * \b INPUT:
     * \param colnames  empty CPL array for strings
     * \param filename  name of ASCII file
     *
     * \b OUTPUT:
     * \param colnames  CPL array of column names for reading of ASCII file
     *
     * \b ERRORS:
     * - Invalid object structure
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    char errtxt[MF_MAXLEN] = "", str[MF_MAXLEN] = "", *colname;
    int j = 0;

    /* Check type of CPL array */
    if (cpl_array_get_type(colnames) != CPL_TYPE_STRING) {
        sprintf(errtxt, "%s: cpl_array *colnames (no string array)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check file existence */
    if ((stream = fopen(filename, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Read header line */
    if (fgets(str, MF_MAXLEN, stream)) {}

    /* Close ASCII file */
    fclose(stream);

    /* Check first character */
    if (str[0] != '#') {
        sprintf(errtxt, "%s: %s (no '#' header line)", MF_ERROR_UFS_TXT,
                filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Split header line into column names and write them into CPL array */
    strtok(str, "\n\t "); // remove hash
    colname = strtok(NULL, "\n\t ");
    while (colname != NULL) {
        j++;
        cpl_array_set_size(colnames, j);
        cpl_array_set_string(colnames, j-1, colname);
        colname = strtok(NULL, "\n\t ");
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_performtac_tarr(mftarr *tabdat,
                                       const mftarr *reftabdat,
                                       const char *fluxcolname,
                                       const char *dfluxcolname)
{
    /*!
     * \callgraph
     *
     * Corrects a spectrum provided as ::mftarr structure for telluric
     * absorption. The required transmission curve has to be provided by
     * another ::mftarr structure with the same flux column name. The column
     * for the transmission curve has to be 'mtrans'. It is copied into the
     * input ::mftarr structure, which is also supplemented by a column for
     * the corrected flux ('tacflux') and for a quality flag ('tacqual'). If
     * an flux error column exists in the input ::mftarr structure, the
     * corrected error is written into a new 'tacdflux' column.
     *
     * \b INPUT:
     * \param tabdat        input ::mftarr structure
     * \param reftabdat     ::mftarr structure containing transmission curve
     * \param fluxcolname   name of flux column in both ::mftarr structures
     * \param dfluxcolname  name of error column in first ::mftarr structure
     *
     * \b OUTPUT:
     * \param tabdat        ::mftarr structure with corrected flux
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_table *spec;

    /* Convert mftarr structure into CPL table for telluric absorption
       correction */
    spec = cpl_table_new(0);
    if ((status = mf_corr_tarr2tactable(spec, tabdat, reftabdat, fluxcolname,
                                        dfluxcolname)) != CPL_ERROR_NONE) {
        cpl_table_delete(spec);
        return status;
    }

    /* Perform telluric absorption correction */
    if ((status = mf_corr_calctac(spec, 1)) != CPL_ERROR_NONE) {
        cpl_table_delete(spec);
        return status;
    }

    /* Convert CPL table for telluric absorption correction into mftarr
       structure */
    if ((status = mf_corr_tactable2tarr(tabdat, spec)) != CPL_ERROR_NONE) {
        cpl_table_delete(spec);
        return status;
    }

    /* Free allocated memory */
    cpl_table_delete(spec);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_performtac_varr(mfvarr *vecdat, mfvarr *transdat,
                                       const cpl_table *extnames)
{
    /*!
     * \callgraph
     *
     * Corrects a spectrum provided as ::mfvarr structure for telluric
     * absorption. The required transmission curve has to be provided by
     * another ::mfvarr structure without extensions. The extension names and
     * numbers have to be provided by a CPL table. Apart from flux correction,
     * the routine adapts a possible quality extension to indicate difficult
     * wavelengths for the flux correction. Finally, the input transmission
     * curve is converted into the resulting correction function.
     *
     * \b INPUT:
     * \param vecdat    input ::mfvarr structure
     * \param transdat  ::mfvarr structure containing transmission curve
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b OUTPUT:
     * \param vecdat    ::mfvarr structure with corrected flux
     * \param transdat  ::mfvarr structure with correction function
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_table *spec;

    /* Convert mfvarr structure into CPL table for telluric absorption
       correction */
    spec = cpl_table_new(0);
    if ((status = mf_corr_varr2tactable(spec, vecdat, transdat, extnames))
        != CPL_ERROR_NONE) {
        cpl_table_delete(spec);
        return status;
    }

    /* Perform telluric absorption correction */
    if ((status = mf_corr_calctac(spec, 1)) != CPL_ERROR_NONE) {
        cpl_table_delete(spec);
        return status;
    }

    /* Convert CPL table for telluric absorption correction into mfvarr
       structure */
    if ((status = mf_corr_tactable2varr(vecdat, transdat, extnames, spec))
        != CPL_ERROR_NONE) {
        cpl_table_delete(spec);
        return status;
    }

    /* Free allocated memory */
    cpl_table_delete(spec);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_tarr2tactable(cpl_table *spec, const mftarr *tabdat,
                                     const mftarr *reftabdat,
                                     const char *fluxcolname,
                                     const char *dfluxcolname)
{
    /*!
     * Prepares a CPL table for the correction of telluric absorption from
     * two ::mftarr structures. The first one has to contain the flux column
     * (and possibly an error column) and the other one has to provide the
     * transmission curve. If there are multiple chips, the individual spectra
     * will be merged.
     *
     * \b INPUT:
     * \param spec          empty CPL table
     * \param tabdat        ::mftarr structure containing spectral flux
     * \param reftabdat     ::mftarr structure containing transmission curve
     * \param fluxcolname   name of flux column in both ::mftarr structures
     * \param dfluxcolname  name of error column in first ::mftarr structure
     *
     * \b OUTPUT:
     * \param spec          CPL table with flux, (error,) and transmission
     *                      columns
     *
     * \b ERRORS:
     * - Invalid object structure
     * - No data
     */

    char errtxt[MF_MAXLEN] = "";
    int next = 0, nrow = 0, ndat = 0, h = 0, j = 0, i = 0;

    /* Check correspondence of extension numbers in input mftarr structures */
    next = tabdat->next;
    if (next <= 0) {
        sprintf(errtxt, "%s: mftarr *tabdat (no extensions)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    } else if (next > 0 && reftabdat->next != next) {
        sprintf(errtxt, "%s: mftarr *reftabdat != mftarr *tabdat (number of "
                "extensions)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check correspondence of data points in input mftarr structures
       (check first extension only) */
    nrow = cpl_table_get_nrow(tabdat->tab[1]);
    if (nrow <= 0) {
        sprintf(errtxt, "%s: mftarr *tabdat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    } else if (nrow > 0 && cpl_table_get_nrow(reftabdat->tab[1]) != nrow) {
        sprintf(errtxt, "%s: mftarr *reftabdat != mftarr *tabdat (number of "
                "data points)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of required flux input column */
    if (cpl_table_has_column(tabdat->tab[1], fluxcolname) != 1) {
        sprintf(errtxt, "%s: mftarr *tabdat (no flux column of expected "
                "name)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check correspondence of flux column names in input mftarr structures */
    if (cpl_table_has_column(reftabdat->tab[1], fluxcolname) != 1) {
        sprintf(errtxt, "%s: mftarr *reftabdat != mftarr *tabdat (name of "
                "flux column)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check existence of required transmission input column */
    if (cpl_table_has_column(reftabdat->tab[1], "mtrans") != 1) {
        sprintf(errtxt, "%s: mftarr *tabdat (no transmission column of "
                "expected name)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Create columns in output CPL table */
    if (cpl_table_has_column(spec, "chip") != 1) {
        cpl_table_new_column(spec, "chip", CPL_TYPE_INT);
    }
    if (cpl_table_has_column(spec, "flux") != 1) {
        cpl_table_new_column(spec, "flux", CPL_TYPE_DOUBLE);
    }
    if (cpl_table_has_column(tabdat->tab[1], dfluxcolname) == 1) {
        /* Create error column only if it is present in input table */
        if (cpl_table_has_column(spec, "dflux") != 1) {
            cpl_table_new_column(spec, "dflux", CPL_TYPE_DOUBLE);
        }
    }
    if (cpl_table_has_column(spec, "mtrans") != 1) {
        cpl_table_new_column(spec, "mtrans", CPL_TYPE_DOUBLE);
    }

    /* Set size of output CPL table */
    ndat = nrow * next;
    cpl_table_set_size(spec, ndat);

    /* Copy data from mftarr structures to output CPL table */
    for (h = 0, j = 0; j < next; j++) {
        for (i = 0; i < nrow; i++) {
            cpl_table_set(spec, "chip", h, j+1);
            cpl_table_set(spec, "flux", h,
                          cpl_table_get(tabdat->tab[j+1], fluxcolname, i,
                                        NULL));
            if (cpl_table_has_column(tabdat->tab[1], dfluxcolname) == 1) {
                cpl_table_set(spec, "dflux", h,
                              cpl_table_get(tabdat->tab[j+1], dfluxcolname, i,
                                            NULL));
            }
            cpl_table_set(spec, "mtrans", h,
                          cpl_table_get(reftabdat->tab[j+1], "mtrans", i,
                                        NULL));
            h++;
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_varr2tactable(cpl_table *spec, const mfvarr *vecdat,
                                     const mfvarr *transdat,
                                     const cpl_table *extnames)
{
    /*!
     * Prepares a CPL table for the correction of telluric absorption from
     * two ::mfvarr structures. The first one has to contain a flux vector and
     * possibly a flux error and/or quality vector. The second one has to
     * provide the transmission curve.
     *
     * \b INPUT:
     * \param spec      empty CPL table
     * \param vecdat    ::mfvarr structure containing spectral flux
     * \param transdat  ::mfvarr structure containing transmission curve
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b OUTPUT:
     * \param spec      CPL table with flux and transmission columns
     *
     * \b ERRORS:
     * - Invalid object structure
     * - No data
     */

    char errtxt[MF_MAXLEN] = "";
    int next = 0, nrow = 0, extn_flux = 0, extn_dflux = 0, i = 0;

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

   /* Check extension numbers in input mfvarr structures */
    next = vecdat->next;

    // printf("DEBUG %d, %d\n", next, cpl_table_get_column_max(extnames, "extn"));
    // if (next != cpl_table_get_column_max(extnames, "extn")) {
    //     sprintf(errtxt, "%s: mfvarr *vecdat (incorrect extension number)",
    //             MF_ERROR_IOS_TXT);
    //     return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    // }
    if (next < cpl_table_get_column_max(extnames, "extn")) {
        sprintf(errtxt, "%s: mfvarr *vecdat (incorrect extension number)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }
    if (transdat->next != 0) {
        sprintf(errtxt, "%s: mfvarr *transdat (extension number != 0)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check correspondence of data points in input mfvarr structures
       (check zeroth extension only) */
    nrow = cpl_vector_get_size(vecdat->vec[0]);
    if (nrow <= 0) {
        sprintf(errtxt, "%s: mfvarr *vecdat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    } else if (nrow > 0 && cpl_vector_get_size(transdat->vec[0]) != nrow) {
        sprintf(errtxt, "%s: mfvarr *transdat != mfvarr *vecdat (number of "
                "data points)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get extension numbers for flux and error */
    extn_flux = cpl_table_get(extnames, "extn", 1, NULL);
    extn_dflux = cpl_table_get(extnames, "extn", 2, NULL);

    /* Create columns in output CPL table */
    if (cpl_table_has_column(spec, "flux") != 1) {
        cpl_table_new_column(spec, "flux", CPL_TYPE_DOUBLE);
    }
    if (extn_dflux >= 0) {
        /* Create error column only if a corresponding FITS extension
           exists */
        if (cpl_table_has_column(spec, "dflux") != 1) {
            cpl_table_new_column(spec, "dflux", CPL_TYPE_DOUBLE);
        }
    }
    if (cpl_table_has_column(spec, "mtrans") != 1) {
        cpl_table_new_column(spec, "mtrans", CPL_TYPE_DOUBLE);
    }

    /* Set size of output CPL table */
    cpl_table_set_size(spec, nrow);

    /* Copy data from mfvarr structures to output CPL table */
    for (i = 0; i < nrow; i++) {
        cpl_table_set(spec, "flux", i,
                      cpl_vector_get(vecdat->vec[extn_flux], i));
        if (extn_dflux >= 0) {
            cpl_table_set(spec, "dflux", i,
                          cpl_vector_get(vecdat->vec[extn_dflux], i));
        }
        cpl_table_set(spec, "mtrans", i,
                      cpl_vector_get(transdat->vec[0], i));
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_calctac(cpl_table *spec, const int trans)
{
    /*!
     * Corrects the observed flux (and error if present) for telluric
     * absorption. The corrected fluxes (and errors) are given in column
     * "cflux" (and "cdflux"). The column "qual" indicates whether the
     * correction is relatively reliable for a pixel (= 1) or whether it
     * probably failed (= 0). The minimum transmission ::MF_MINTRANS is used
     * for the classification. The columns "cflux" and "qual" are created if
     * they do not exist. The existence of the input columns "flux" and
     * "mtrans" is mandatory. The column "cdflux" is only created if it does
     * not exist and an input column labelled "dflux" is present. No flux
     * correction is carried out in the case of a sky emission spectrum
     * (\e trans != 1). In this case, the column(s) "cflux" (and "cdflux")
     * contain(s) the uncorrected input flux (and error respectively).
     *
     * \b INPUT:
     * \param spec   CPL table with observed flux and modelled transmission
     *               curve
     * \param trans  transmission (= 1) or emission (= 0 or 2)
     *
     * \b OUTPUT:
     * \param spec   table with columns for corrected flux, (error,) and
     *               quality flag
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     */

    char errtxt[MF_MAXLEN];
    int nrow = 0, i = 0;
    double minflux = 0., mintrans = 0.;
    double *flux, *mtrans, *cflux, *dflux = NULL, *cdflux = NULL;

    /* Check numer of data points */
    nrow = cpl_table_get_nrow(spec);
    if (nrow <= 0) {
        sprintf(errtxt, "%s: cpl_table *spec", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check existence of required input columns */
    if (cpl_table_has_column(spec, "flux") != 1 ||
        cpl_table_has_column(spec, "mtrans") != 1) {
        sprintf(errtxt, "%s: cpl_table *spec (column names)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Create column for corrected flux and set values to observed flux */
    if (cpl_table_has_column(spec, "cflux") != 1) {
        cpl_table_duplicate_column(spec, "cflux", spec, "flux");
    }

    /* Check existence of error column and create column for corrected error
       if required */
    if (cpl_table_has_column(spec, "dflux") == 1) {
        if (cpl_table_has_column(spec, "cdflux") != 1) {
            cpl_table_duplicate_column(spec, "cdflux", spec, "dflux");
        }
    }

    /* Create column for quality flag of telluric absorption correction */
    if (cpl_table_has_column(spec, "qual") != 1) {
        cpl_table_new_column(spec, "qual", CPL_TYPE_INT);
    }

    /* Get pointers to CPL table columns */
    flux = cpl_table_get_data_double(spec, "flux");
    mtrans = cpl_table_get_data_double(spec, "mtrans");
    cflux = cpl_table_get_data_double(spec, "cflux");
    if (cpl_table_has_column(spec, "dflux") == 1) {
        dflux = cpl_table_get_data_double(spec, "dflux");
        cdflux = cpl_table_get_data_double(spec, "cdflux");
    }

    /* Get criteria for minimum observed flux and transmission */
    minflux = MF_MINTRANS * cpl_table_get_column_median(spec, "flux");
    mintrans = MF_MINTRANS;

    /* Divide model spectrum by transmission curve and set quality flag */
    for (i = 0; i < nrow; i++) {
        if (mtrans[i] < mintrans) {
            cpl_table_set(spec, "qual", i, 0);
        } else {
            cpl_table_set(spec, "qual", i, 1);
        }
        if (flux[i] < minflux && trans == 1) {
            cpl_table_set(spec, "qual", i, 0);
        }
        if (trans != 1) {
            /* No correction in the case of a sky emission spectrum */
            cflux[i] = flux[i];
            if (cpl_table_has_column(spec, "dflux") == 1) {
                cdflux[i] = dflux[i];
            }
        } else {
            if (mtrans[i] < mintrans) {
                 /* Avoid extreme correction factors */
                 cflux[i] = flux[i] / mintrans;
                 if (cpl_table_has_column(spec, "dflux") == 1) {
                     cdflux[i] = dflux[i] / mintrans;
                 }
            } else {
                 cflux[i] = flux[i] / mtrans[i];
                 if (cpl_table_has_column(spec, "dflux") == 1) {
                     cdflux[i] = dflux[i] / mtrans[i];
                 }
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_tactable2tarr(mftarr *tabdat, const cpl_table *spec)
{
    /*!
     * Writes results of telluric absorption correction as provided by a CPL
     * table of fixed format into the input ::mftarr structure. The columns
     * "mtrans" (transmission curve), "tacflux" (corrected flux), and
     * "tacqual" (quality of correction) are added. "tacdflux" is only added
     * if a column for corrected flux errors exists in the results table.
     *
     * \b INPUT:
     * \param tabdat  ::mftarr structure containing spectral flux
     * \param spec    CPL table with results of telluric absorption correction
     *
     * \b OUTPUT:
     * \param tabdat  ::mftarr structure containing corrected flux
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     * - Invalid object value(s)
     */

    char errtxt[MF_MAXLEN] = "";
    int ndat = 0, next = 0, nrow = 0, j = 0, h = 0, i = 0;

    /* Check number of data points in input CPL table */
    ndat = cpl_table_get_nrow(spec);
    if (ndat <= 0) {
        sprintf(errtxt, "%s: cpl_table *spec", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check existence of required columns in input CPL table */
    if (cpl_table_has_column(spec, "chip") != 1 ||
        cpl_table_has_column(spec, "mtrans") != 1 ||
        cpl_table_has_column(spec, "cflux") != 1 ||
        cpl_table_has_column(spec, "qual") != 1) {
        sprintf(errtxt, "%s: cpl_table *spec (missing column(s))",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check number of extensions in mftarr structure */
    next = tabdat->next;
    if (next <= 0) {
        sprintf(errtxt, "%s: mftarr *tabdat (no extensions)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check correspondence of extensions in mftarr structure and chips in
       input CPL table */
    if (cpl_table_get_column_max(spec, "chip") != next) {
        sprintf(errtxt, "%s: mftarr *tabdat != cpl_table *spec (number of "
                "chips)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check number of data points in mftarr structure */
    nrow = cpl_table_get_nrow(tabdat->tab[1]);
    if (nrow <= 0) {
        sprintf(errtxt, "%s: mftarr *tabdat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check correspondence of data point numbers */
    if (nrow * next != ndat) {
        sprintf(errtxt, "%s: mftarr *tabdat != cpl_table *spec (number of "
                "data points)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Create columns in mftarr structure */
    for (j = 0; j < next; j++) {
        if (cpl_table_has_column(tabdat->tab[j+1], "mtrans") != 1) {
            cpl_table_new_column(tabdat->tab[j+1], "mtrans", CPL_TYPE_DOUBLE);
        }
        if (cpl_table_has_column(tabdat->tab[j+1], "tacflux") != 1) {
            cpl_table_new_column(tabdat->tab[j+1], "tacflux",
                                 CPL_TYPE_DOUBLE);
        }
        if (cpl_table_has_column(spec, "cdflux") == 1) {
            /* Create error column if corrected flux errors exist */
            if (cpl_table_has_column(tabdat->tab[j+1], "tacdflux") != 1) {
                cpl_table_new_column(tabdat->tab[j+1], "tacdflux",
                                     CPL_TYPE_DOUBLE);
            }
        }
        if (cpl_table_has_column(tabdat->tab[j+1], "tacqual") != 1) {
            cpl_table_new_column(tabdat->tab[j+1], "tacqual", CPL_TYPE_INT);
        }
    }

    /* Copy data from input CPL table to mftarr structure */
    for (h = 0, j = 0; j < next; j++) {
        for (i = 0; i < nrow; i++) {
            if (cpl_table_get(spec, "chip", h, NULL) != j+1) {
                sprintf(errtxt, "%s: cpl_table *spec (unexpected chip)",
                        MF_ERROR_IOV_TXT);
                return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s",
                                             errtxt);
            }
            cpl_table_set(tabdat->tab[j+1], "mtrans", i,
                          cpl_table_get(spec, "mtrans", h, NULL));
            cpl_table_set(tabdat->tab[j+1], "tacflux", i,
                          cpl_table_get(spec, "cflux", h, NULL));
            if (cpl_table_has_column(spec, "cdflux") == 1) {
                cpl_table_set(tabdat->tab[j+1], "tacdflux", i,
                              cpl_table_get(spec, "cdflux", h, NULL));
            }
            cpl_table_set(tabdat->tab[j+1], "tacqual", i,
                          cpl_table_get(spec, "qual", h, NULL));
            h++;
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_tactable2varr(mfvarr *vecdat, mfvarr *transdat,
                                     const cpl_table *extnames,
                                     const cpl_table *spec)
{
    /*!
     * Writes results of telluric absorption correction as provided by a CPL
     * table of fixed format into the two input ::mfvarr structures. The
     * first one (\e vecdat) will contain the corrected flux, and the corrected
     * error flux and an updated quality flag if present. The second one
     * (\e transdat) will contain the correction function.
     *
     * \b INPUT:
     * \param vecdat    ::mfvarr structure containing spectral flux
     * \param transdat  ::mfvarr structure containing transmission curve
     * \param extnames  CPL table with FITS extension numbers and names
     * \param spec      CPL table with results of telluric absorption
     *                  correction
     *
     * \b OUTPUT:
     * \param vecdat    ::mfvarr structure containing corrected flux
     * \param transdat  ::mfvarr structure containing correction function
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     * - Invalid object value(s)
     */

    char errtxt[MF_MAXLEN] = "";
    int ndat = 0, next = 0, nrow = 0, extn_flux = 0, extn_dflux = 0;
    int extn_mask = 0, i = 0, qual = 0;
    double maskval[2] = {0., 0.}, mask = 0.;

    /* Check number of data points in input CPL table */
    ndat = cpl_table_get_nrow(spec);
    if (ndat <= 0) {
        sprintf(errtxt, "%s: cpl_table *spec", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check existence of required columns in input CPL table */
    if (cpl_table_has_column(spec, "mtrans") != 1 ||
        cpl_table_has_column(spec, "cflux") != 1 ||
        cpl_table_has_column(spec, "qual") != 1) {
        sprintf(errtxt, "%s: cpl_table *spec (missing column(s))",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

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

    /* Check extension numbers in input mfvarr structures */
    next = vecdat->next;

    // printf("DEBUG %d, %d\n", next, cpl_table_get_column_max(extnames, "extn"));
    // if (next != cpl_table_get_column_max(extnames, "extn")) {
    //     sprintf(errtxt, "%s: mfvarr *vecdat (incorrect extension number)",
    //             MF_ERROR_IOS_TXT);
    //     return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    // }
    if (next < cpl_table_get_column_max(extnames, "extn")) {
        sprintf(errtxt, "%s: mfvarr *vecdat (incorrect extension number)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }
    if (transdat->next != 0) {
        sprintf(errtxt, "%s: mfvarr *transdat (extension number != 0)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check correspondence of data points in input mfvarr structures
       (check zeroth extension only) */
    nrow = cpl_vector_get_size(vecdat->vec[0]);
    if (nrow <= 0) {
        sprintf(errtxt, "%s: mfvarr *vecdat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    } else if (nrow > 0 && cpl_vector_get_size(transdat->vec[0]) != nrow) {
        sprintf(errtxt, "%s: mfvarr *transdat != mfvarr *vecdat (number of "
                "data points)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check correspondence of data point numbers in input CPL table and
       mfvarr structures */
    if (nrow != ndat) {
        sprintf(errtxt, "%s: mfvarr *vecdat != cpl_table *spec (number of "
                "data points)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get extension numbers for flux, error, and mask */
    extn_flux = cpl_table_get(extnames, "extn", 1, NULL);
    extn_dflux = cpl_table_get(extnames, "extn", 2, NULL);
    extn_mask = cpl_table_get(extnames, "extn", 3, NULL);

    /* Check existence of flux error column in input CPL table if the mfvarr
       structure contains error data */
    if (extn_dflux >= 0) {
        if (cpl_table_has_column(spec, "cdflux") != 1) {
            sprintf(errtxt, "%s: mfvarr *vecdat != cpl_table *spec (flux "
                    "error data)", MF_ERROR_IOS_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
        }
    }

    /* Get meaning of mask values (maskval[0] -> bad, maskval[1] -> good) */
    if (extn_mask >= 0) {
        if (mf_basic_getmaskval_vector(maskval, vecdat->vec[extn_mask]) !=
            CPL_ERROR_NONE) {
            return MF_ERROR_IOV;
        }
    }

    /* Copy data from input CPL table to mfvarr structures */
    for (i = 0; i < nrow; i++) {
        cpl_vector_set(transdat->vec[0], i,
                       cpl_table_get(spec, "mtrans", i, NULL));
        cpl_vector_set(vecdat->vec[extn_flux], i,
                       cpl_table_get(spec, "cflux", i, NULL));
        if (extn_dflux >= 0) {
            cpl_vector_set(vecdat->vec[extn_dflux], i,
                           cpl_table_get(spec, "cdflux", i, NULL));
        }
        if (extn_mask >= 0) {
            qual = cpl_table_get(spec, "qual", i, NULL);
            mask = cpl_vector_get(vecdat->vec[extn_mask], i);
            if (mask == maskval[1] && qual == 0) {
                if (maskval[0] == maskval[1]) {
                    /* Define bad mask value if unknown */
                    if (maskval[1] == 1.) {
                        cpl_vector_set(vecdat->vec[extn_mask], i, 0.);
                    } else {
                        cpl_vector_set(vecdat->vec[extn_mask], i,
                                       MF_MASKVAL_UNDEF);
                    }
                } else {
                    cpl_vector_set(vecdat->vec[extn_mask], i, maskval[0]);
                }
            }
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_corr_tacvarr2iarr(mfiarr *imadat, const mfvarr *vecdat,
                                    const mfvarr *transdat,
                                    const cpl_table *extnames)
{
    /*!
     * Writes results of telluric absorption correction as provided by the
     * ::mfvarr structures \e vecdat (corrected flux and quality flag,
     * possibly corrected flux error) and \e transdat (correction function)
     * into the input ::mfiarr structure. The 1D correction function is
     * applied to all rows of the 2D spectrum individually. If the quality of
     * the telluric absorption correction is bad, all mask values of a column
     * are set to the corresponding value for non-selection. In the other
     * case, no manipulation of mask pixels is performed.
     *
     * \b INPUT:
     * \param imadat    ::mfiarr structure containing input flux
     * \param vecdat    ::mfvarr structure containing corrected flux
     * \param transdat  ::mfvarr structure containing correction function
     * \param extnames  CPL table with FITS extension numbers and names
     *
     * \b OUTPUT:
     * \param imadat    ::mfiarr structure containing corrected flux
     *
     * \b ERRORS:
     * - No data
     * - Invalid object structure
     * - Invalid object value(s)
     */

    char errtxt[MF_MAXLEN] = "";
    int nx = 0, ny = 0, next = 0, nrow = 0, extn_flux = 0, extn_dflux = 0;
    int extn_mask = 0, i = 0, j = 0, qual = 0;
    double maskval[2] = {0., 0.}, refmaskval[2] = {0., 0.}, mtrans = 0.;
    double refmask = 0., flux = 0., dflux = 0., mask = 0.;

    /* Check number of data points in x and y direction in input mfiarr
       structure */
    nx = cpl_image_get_size_x(imadat->ima[0]);
    ny = cpl_image_get_size_y(imadat->ima[0]);
    if (nx <= 0 || ny <= 0) {
        sprintf(errtxt, "%s: mfiarr *imadat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

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

    /* Check extension numbers in input mfvarr structures */
    next = vecdat->next;
    // if (next != cpl_table_get_column_max(extnames, "extn")) {
    if (next < cpl_table_get_column_max(extnames, "extn")) {
        sprintf(errtxt, "%s: mfvarr *vecdat (incorrect extension number)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }
    if (transdat->next != 0) {
        sprintf(errtxt, "%s: mfvarr *transdat (extension number != 0)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check correspondence of extensions in mfiarr and mfvarr structures */
    if (imadat->next != next) {
        sprintf(errtxt, "%s: mfvarr *vecdat != mfiarr *imadat (number of "
                "extensions)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check correspondence of data points in input mfvarr structures
       (check zeroth extension only) */
    nrow = cpl_vector_get_size(vecdat->vec[0]);
    if (nrow == 0) {
        sprintf(errtxt, "%s: mfvarr *vecdat", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    } else if (nrow > 0 && cpl_vector_get_size(transdat->vec[0]) != nrow) {
        sprintf(errtxt, "%s: mfvarr *transdat != mfvarr *vecdat (number of "
                "data points)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Check correspondence of data point numbers in mfiarr and mfvarr
       structures */
    if (nrow != nx) {
        sprintf(errtxt, "%s: mfvarr *vecdat != mfiarr *imadat (number of "
                "data points)", MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get extension numbers for flux and mask */
    extn_flux = cpl_table_get(extnames, "extn", 1, NULL);
    extn_dflux = cpl_table_get(extnames, "extn", 2, NULL);
    extn_mask = cpl_table_get(extnames, "extn", 3, NULL);

    /* Get meaning of mask values (maskval[0] -> bad, maskval[1] -> good) */
    if (extn_mask >= 0) {
        if (mf_basic_getmaskval_image(maskval, imadat->ima[extn_mask]) !=
            CPL_ERROR_NONE) {
            return MF_ERROR_IOV;
        }
    }

    /* Check correspondence of mask values in mfiarr and mfvarr structures */
    if (extn_mask >= 0) {
        if (mf_basic_getmaskval_vector(refmaskval, vecdat->vec[extn_mask]) !=
            CPL_ERROR_NONE) {
            return MF_ERROR_IOV;
        }
        if ((refmaskval[0] != maskval[0] && refmaskval[0] == 0.) ||
            refmaskval[1] != refmaskval[1]) {
            sprintf(errtxt, "%s: mfvarr *vecdat != mfiarr *imadat (mask "
                    "values)", MF_ERROR_IOV_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s",
                                         errtxt);
        }
    }

    /* Use correction function from mfvarr structure to correct flux in mfiarr
       structure and set bad pixels according to mask in mfvarr structure */
    for (i = 0; i < nx; i++) {
        mtrans = cpl_vector_get(transdat->vec[0], i);
        if (extn_mask >= 0) {
            refmask = cpl_vector_get(vecdat->vec[extn_mask], i);
            if (refmask < 0.) {
                /* For lost bad mask value (temporarily set to -1) */
                refmask = maskval[0];
            }
        }
        for (j = 0; j < ny; j++) {
            flux = cpl_image_get(imadat->ima[extn_flux], i+1, j+1, &qual);
            cpl_image_set(imadat->ima[extn_flux], i+1, j+1, flux / mtrans);
            if (extn_dflux >= 0) {
                dflux = cpl_image_get(imadat->ima[extn_dflux], i+1, j+1,
                                      &qual);
                cpl_image_set(imadat->ima[extn_dflux], i+1, j+1,
                              dflux / mtrans);
            }
            if (extn_mask >= 0) {
                mask = cpl_image_get(imadat->ima[extn_mask], i+1, j+1, &qual);
                if (mask != refmask) {
                    cpl_image_set(imadat->ima[extn_mask], i+1, j+1, refmask);
                }
            }
        }
    }

    return CPL_ERROR_NONE;
}

/**@}*/
