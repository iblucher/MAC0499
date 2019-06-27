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
 * \file mf_lib.c
 *
 * Function library for rebinning of wrapper output spectra of the LBLRTM code
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  16 Sep 2009
 * \date   28 Aug 2013
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include <mf_lib.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code mf_lib_createlibspec(const char *inpath, const char *outpath,
                                    const char *outname, const char *suffix,
                                    const double limlam[2],
                                    const double resol, const mfdrv *drvpar)
{
    /*!
     * \callgraph
     *
     * Merges and rebins a set of LBLRTM radiance ("*R.*") or transmission
     * ("*T.*") spectra. The final spectrum is written to an ASCII ("*.dat" or
     * "*.ascii") or FITS file ("*.fits" or "*.mt").
     *
     * \b INPUT:
     * \param inpath   path name for input LBLRTM data
     *                 (naming scheme examples: TAPE27_8 or TAPE28_16)
     * \param outpath  path name for the output radiance and transmission
     *                 spectra
     * \param outname  name of the output spectra
     *                 (supplemented by "R" and "T", respectively)
     * \param suffix   file name extension
     *                 (defined: "dat", "ascii", "fits", or "mt")
     * \param limlam   wavelength limits for the output spectra
     * \param resol    resolution of output spectra
     *                 \f$\to\f$ bin sizes proportional to wavelength
     * \param drvpar   ::mfdrv parameter structure
     *
     * \b ERRORS:
     * - Invalid file name extension
     * - see subroutines as well
     */

    FILE *stream;
    cpl_error_code status = CPL_ERROR_NONE;
    cpl_table *spec;
    char errtxt[MF_MAXLEN], filename[MF_MAXLEN];
    cpl_parameter *p;

    /* Initialise wavelength grid */
    spec = cpl_table_new(0);
    cpl_table_new_column(spec, "lambda", CPL_TYPE_DOUBLE);
    cpl_table_new_column(spec, "flux", CPL_TYPE_DOUBLE);
    if ((status = mf_lib_createlamgrid(spec, limlam, resol)) !=
        CPL_ERROR_NONE) {
        cpl_table_delete(spec);
        return status;
    }

    p = cpl_parameterlist_find(drvpar->parlist, "trans");
    if (cpl_parameter_get_int(p) != 0) {
        /* Rebin set of transmission spectra */
        mf_lib_rebmolspecall(spec, inpath, MF_TRANSNAME);

        /* Create output file names */
        sprintf(filename, "%s/%s_T.%s", outpath, outname, suffix);
    } else {
        /* Rebin set of radiance spectra */
        mf_lib_rebmolspecall(spec, inpath, MF_RADNAME);

        /* Create output file names */
        sprintf(filename, "%s/%s_R.%s", outpath, outname, suffix);
    }

    /* Check suffix for file class and write to file */
    if (strncmp(suffix, "dat", strlen(suffix)) == 0 ||
        strncmp(suffix, "ascii", strlen(suffix)) == 0 ) {
        /* Write resulting spectrum to ASCII file */
        stream = fopen(filename, "w+");
        cpl_table_dump(spec, 0, cpl_table_get_nrow(spec), stream);
        fclose(stream);
    } else if (strncmp(suffix, "fits", strlen(suffix)) == 0 ||
        strncmp(suffix, "mt", strlen(suffix)) == 0 ) {
        /* Write resulting spectrum to FITS file */
        cpl_table_save(spec, NULL, NULL, filename, CPL_IO_CREATE);
    } else {
        sprintf(errtxt, "%s: %s ('dat', 'ascii', 'fits', or 'mt' only)",
                MF_ERROR_IFE_TXT, suffix);
        status = cpl_error_set_message(cpl_func, MF_ERROR_IFE, "%s", errtxt);
    }

    /* Free allocated memory */
    cpl_table_delete(spec);

    return status;
}


cpl_error_code mf_lib_createlamgrid(cpl_table *outspec,
                                    const double limlam[2],
                                    const double resol)
{
    /*!
     * Initialisation of a wavelength grid that is characterised by constant
     * resolution lambda / delta lambda
     *
     * \b INPUT:
     * \param limlam  lower and upper limit of wavelength range
     * \param resol   resolution lambda / delta lambda
     *
     * \b OUTPUT:
     * \param outspec  CPL table that provides wavelengths and fluxes (= 0.)
     *                 of a spectrum
     *
     * \b ERRORS:
     * - Invalid input parameter(s)
     */

    char errtxt[MF_MAXLEN];
    int nrow = 0, i = 0;
    double lam = 0.;

    /* Test input parameters */
    if (limlam[1] <= limlam[0] || resol <= 0.) {
        /* Return spectrum with zero data points */
        sprintf(errtxt, "%s: limlam[1] <= limlam[0] || resol <= 0. "
                "(wavelength grid)", MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Get number of data points and resize CPL table */
    nrow = ceil(log(limlam[1] / limlam[0]) / log1p(1 / resol));
    cpl_table_set_size(outspec, nrow);

    /* Create wavelength grid */
    lam = limlam[0];
    cpl_table_set(outspec, "lambda", i, lam);
    for (i = 1; i < nrow; i++) {
        lam *= (1 + 1 / resol);
        cpl_table_set(outspec, "lambda", i, lam);
    }

    /* Set flux = 0 */
    cpl_table_fill_column_window(outspec, "flux", 0, nrow, 0.);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_lib_rebmolspecall(cpl_table *spec, const char *path,
                                    const char *filename)
{
    /*!
     * \callgraph
     * Rebins a set of LBLRTM spectra (wrapper output) in wavelength units
     * [\f$\mu{\rm m}\f$] (variable step size possible).
     *
     * \b INPUT:
     * \param spec      CPL table with wavelength grid
     * \param path      path name that contains LBLRTM data
     * \param filename  name of LBLRTM data file (TAPE27 or TAPE28;
     *                  for final file names routine adds "_" and file number)
     *                  consisting of wavelength and radiance/transmission
     *                  columns
     *
     * \b OUTPUT:
     * \param spec      CPL table with input wavelength grid and resulting
     *                  fluxes
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_array *klimall;
    char name[MF_MAXLEN];
    int ifile;
    double klim[2];

    /* Find number of files and wavenumber ranges */
    klimall = cpl_array_new(0, CPL_TYPE_DOUBLE);
    if ((status = mf_lib_findklim(klimall, path, filename)) !=
        CPL_ERROR_NONE) {
        /* Leave programme in case of errors */
        cpl_array_delete(klimall);
        return status;
    }

    /* Rebin spectra */
    for (ifile = 1; ifile <= cpl_array_get_size(klimall) - 1; ifile++) {
        sprintf(name, "%s/%s_%d", path, filename, ifile);
        klim[0] = cpl_array_get(klimall, ifile-1, NULL);
        klim[1] = cpl_array_get(klimall, ifile, NULL);
        cpl_msg_info(cpl_func, "Process LBLRTM output file %s_%d",
                     filename, ifile);
        mf_lib_rebmolspec(spec, name, klim);
    }

    /* Free allocated memory */
    cpl_array_delete(klimall);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_lib_findklim(cpl_array *klimall, const char *path,
                               const char *filename)
{
    /*!
     * Finds the wavenumber ranges of a set of LBLRTM output files
     *
     * \b INPUT:
     * \param path      directory that contains the LBLRTM output files
     * \param filename  standard filename (usually "TAPE27" or "TAPE28";
     *                  the name is supplemented by "_" and the file number)
     *
     * \b OUTPUT:
     * \param klimall   CPL array of wavenumbers
     *
     * \b ERRORS:
     * - File opening failed
     * - Unexpected file structure
     * - Invalid order of data points
     */

    FILE *stream = NULL;
    cpl_boolean v1line = CPL_FALSE, linestruct = CPL_TRUE;
    char name[MF_MAXLEN], errtxt[MF_MAXLEN], line[MF_MAXLEN];
    char *dummy = NULL, *str = NULL;
    int ifile = 0, nfile = 0, j = 0, i = 0;
    double *kmin, *kmax;

    /* Check file existence to get number of files */

    do {

        ifile++;
        sprintf(name, "%s/%s_%d", path, filename, ifile);
        stream = fopen(name, "r");

        if (stream == NULL) {
            if (ifile == 1) {
                cpl_array_set_size(klimall, 0);
                sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
                return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s",
                                             errtxt);
            } else {
                cpl_array_set_size(klimall, ifile);
                nfile = ifile - 1;
                break;
            }
        } else {
            fclose(stream);
        }

    } while (CPL_TRUE);

    /* Allocate memory for wavenumber limits */

    kmin = (double *) calloc(nfile, sizeof(double));
    kmax = (double *) calloc(nfile, sizeof(double));

    /* Read wavenumber limits */

    for (ifile = 1; ifile <= nfile; ifile++) {

        sprintf(name, "%s/%s_%d", path, filename, ifile);
        stream = fopen(name, "r");

        v1line = CPL_FALSE;
        if ((dummy = fgets(line, MF_MAXLEN, stream))){};

        if (line[0] != '1') {

            linestruct = CPL_FALSE;

        } else {

            do {

                if (strstr(line, "V1 =") != NULL) {

                    v1line = CPL_TRUE;

                    for (j = 0; j < 2; j++) {

                        /* Split read line into strings */
                        str = strtok(line, "\n\t =");
                        for (i = 1; i < 3; i++) {
                            str = strtok(NULL, "\n\t =");
                            if (str == NULL) {
                                linestruct = CPL_FALSE;
                            }
                        }

                        /* Convert wavenumber string into double */
                        if (j == 0) {
                            kmin[ifile-1] = strtod(str, NULL);
                            dummy = fgets(line, MF_MAXLEN, stream);
                        } else {
                            kmax[ifile-1] = strtod(str, NULL);
                        }

                    }

                }

                dummy = fgets(line, MF_MAXLEN, stream);

            } while (v1line == CPL_FALSE);

            if (kmin[ifile-1] == 0 || kmax[ifile-1] == 0) {
                linestruct = CPL_FALSE;
            }

        }

        fclose(stream);

    }

    /* Compute mean wavenumber limits (overlaps!) and write them into
       output array "klimall" */

    cpl_array_set(klimall, 0, kmin[0]);
    for (ifile = 1; ifile < nfile; ifile++) {
        cpl_array_set(klimall, ifile, (kmax[ifile-1] + kmin[ifile]) / 2);
    }
    cpl_array_set(klimall, nfile, kmax[nfile-1]);

    /* Free memory occupied by "kmin" and "kmax" */

    free(kmin);
    kmin = NULL;
    free(kmax);
    kmax = NULL;

    /* Handle file structure errors */

    if (linestruct == CPL_FALSE) {
        sprintf(errtxt, "%s: %s (wavenumber limits missing or = 0)",
                MF_ERROR_UFS_TXT, name);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Test correct order of wavenumbers */

    for (ifile = 1; ifile <= nfile; ifile++) {
        if (cpl_array_get(klimall, ifile-1, NULL) >=
            cpl_array_get(klimall, ifile, NULL)) {
            sprintf(errtxt, "%s: cpl_array *klimall (wavenumber limits of "
                    "LBLRTM output files", MF_ERROR_IOD_TXT);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOD, "%s",
                                         errtxt);
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_lib_rebmolspec(cpl_table *spec, const char *filename,
                                 const double klim[2])
{
    /*!
     * \callgraph
     * Rebins LBLRTM spectra (wrapper output) in wavelength units
     * [\f$\mu{\rm m}\f$] (variable step size possible).
     *
     * \b INPUT:
     * \param spec      CPL table with wavelength grid
     * \param filename  path and name of LBLRTM data file
     * \param klim      wavenumber limits for extraction
     *
     * \b OUTPUT:
     * \param spec      CPL table with input wavelength grid and resulting
     *                  fluxes
     *
     * \b ERRORS:
     * - No data
     * - Invalid input parameter(s)
     * - File opening failed
     * - Unexpected file structure
     */

    FILE *stream;
    char errtxt[MF_MAXLEN], str[MF_MAXLEN], *dummy;
    cpl_boolean endhead = CPL_FALSE, israd = CPL_FALSE, empty = CPL_FALSE;
    cpl_boolean usampl = CPL_FALSE;
    int nrow = 0, j = 0, jmin = -1, jmax = 0, num = 0, ncol = 2;
    int jlim[2] = {0, 0};
    double llim[2] = {0., 0.}, lmin0 = 0., lmin = 0., lmax0 = 0., k = 0.;
    double flux = 0., lam = HUGE_VAL;
    double *lamv = NULL, *fluxv = NULL;

    /* Check number of data points in spectrum */
    nrow = cpl_table_get_nrow(spec);
    if (nrow == 0) {
        sprintf(errtxt, "%s: cpl_table *spec", MF_ERROR_NDA_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_NDA, "%s", errtxt);
    }

    /* Check wavenumber range */
    if (klim[1] <= klim[0]) {
        sprintf(errtxt, "%s: klim[1] <= klim[0] (wavenumber limits)",
                MF_ERROR_IIP_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IIP, "%s", errtxt);
    }

    /* Convert wavenumbers to wavelengths */
    llim[0] = MF_CONV_K_LAM / klim[1];
    llim[1] = MF_CONV_K_LAM / klim[0];

    /* Get pointers to CPL table columns */
    lamv = cpl_table_get_data_double(spec, "lambda");
    fluxv = cpl_table_get_data_double(spec, "flux");

    /* Compare wavelength ranges */
    if (llim[1] < lamv[0] || llim[0] > lamv[nrow-1]) {
        /* Desired wavelength range not covered by the input model spectrum */
        return CPL_ERROR_NONE;
    }

    /* Check file existence */
    if ((stream = fopen(filename, "r")) == NULL) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Skip header */
    if ((dummy = fgets(str, MF_MAXLEN, stream))){};
    if (str[0] != '1') {
        fclose(stream);
        sprintf(errtxt, "%s: %s (first character != '1')",
                MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    } else {
        do {
            if (strstr(str, "WAVENUMBER") != NULL) {
                endhead = CPL_TRUE;
                /* Flux conversion? */
                if (strstr(str, "RADIANCE") != NULL) {
                    if ((israd = CPL_TRUE)){};
                }
            }
            dummy = fgets(str, MF_MAXLEN, stream);
        } while (endhead == CPL_FALSE);
    }

    /*
     * Search for limiting pixels and wavelengths of input grid that
     * correspond to the given wavenumber interval
     * pixels : jmin, jmax
     * wavelengths : lmin0, lmax0
     */

    for (j = 0; j < nrow; j++) {
        if (jmin < 0 && lamv[j] > llim[0]) {
            jmin = j;
            if (jmin == 0) {
                lmin0 = 1.5 * lamv[jmin] - 0.5 * lamv[jmin+1];
            } else {
                lmin0 = (lamv[jmin-1] + lamv[jmin]) / 2;
            }
        }

        if (jmin >= 0 && lamv[j] > llim[1]) {
            jmax = j - 1;
            if (jmax == nrow - 1) {
                lmax0 = 1.5 * lamv[jmax] - 0.5 * lamv[jmax-1];
            } else {
                lmax0 = (lamv[jmax] + lamv[jmax+1]) / 2;
            }
            break;
        } else if (j == nrow - 1 && lamv[j] <= llim[1]) {
            jmax = j;
            lmax0 = 1.5 * lamv[jmax] - 0.5 * lamv[jmax-1];
        }
    }

    /*
     * Read wavelengths and fluxes of input file.
     * Average all flux values inside a bin of the output wavelength grid.
     */

    for (j = jmax; j >= jmin; j--) {

        /*
         * Adapt lmin for next bin.
         * (lmin: grid related, variable, lower lambda limit)
         */

        if (j == jmin) {
            lmin = lmin0;
        } else {
            lmin = (lamv[j-1] + lamv[j]) / 2;
        }

        /*
         * Read new data point(s) from file?
         * YES: take saved flux value as first summand
         */

        if (lam > lmin && ncol == 2) {
            empty = CPL_FALSE;
            if (j == jmax) {
                fluxv[j] = 0.;
                num = 0;
            } else {
                fluxv[j] = flux;
                num = 1;
            }
        } else {
            empty = CPL_TRUE;
            num = 0;
        }

        while (empty == CPL_FALSE &&
               (ncol = fscanf(stream, "%le %le", &k, &flux)) == 2) {
            lam = MF_CONV_K_LAM / k;

            if (lam > lmax0) {
                continue;
            }

            if (lam > lmin) {
                /* Sum up fluxes and count data points */
                fluxv[j] += flux;
                num++;
            } else {
                break;
            }
        }

        if (num != 0) {
            /* Average fluxes */
            fluxv[j] /= (double) num;
            /* Avoid negative values */
            if (fluxv[j] < 0.) {
                fluxv[j] = 0.;
            }
        } else {
            /* Mark "empty" bins */
            fluxv[j] = -HUGE_VAL;
            if (usampl == CPL_FALSE) usampl = CPL_TRUE;
        }

    }

    fclose(stream);

    /*
     * Interpolate "empty" bins by means of the closest bins that contain
     * data points (valid flux values)
     */

    jlim[0] = jmin;
    jlim[1] = jmax;

    if (usampl == CPL_TRUE) {
        mf_lib_interpolspec(spec, jlim);
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_lib_interpolspec(cpl_table *spec, const int pixlim[2])
{
    /*!
     * Linear interpolation for "empty" bins in spectra (marked by -HUGE_VAL).
     *
     * \b INPUT:
     * \param spec    CPL table with spectrum (after rebinning)
     * \param pixlim  pixel limits of spectrum
     *
     * \b OUTPUT:
     * \param spec    spectrum corrected by interpolation
     *
     * \b ERRORS:
     * - No data
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_boolean ok = CPL_TRUE;
    char errtxt[MF_MAXLEN];
    int j = 0, jlim[2] = {-1, -1}, i = 0;
    double m = 0.;
    double *lamv = NULL, *fluxv = NULL;

    /* Get pointers to CPL table columns */
    lamv = cpl_table_get_data_double(spec, "lambda");
    fluxv = cpl_table_get_data_double(spec, "flux");

    /* Perform interpolation/extrapolation procedure */

    for (j = pixlim[0]; j <= pixlim[1]; j++) {

        if (fluxv[j] != -HUGE_VAL) {

            if (ok == CPL_TRUE) {
                jlim[0] = j;
            } else {
                jlim[1] = j;
                if (jlim[0] == -1) {
                    /* Constant value for extrapolation of lower margin */
                    for (i = pixlim[0]; i <= jlim[1] - 1; i++) {
                        fluxv[i] = fluxv[jlim[1]];
                    }
                } else {
                    /* Linear interpolation */
                    m = (fluxv[jlim[1]] - fluxv[jlim[0]]) /
                        (lamv[jlim[1]] - lamv[jlim[0]]);
                    for (i = jlim[0] + 1; i <= jlim[1] - 1; i++) {
                        fluxv[i] = m * (lamv[i] - lamv[jlim[0]]) +
                                   fluxv[jlim[0]];
                    }
                }
                jlim[0] = jlim[1];
                ok = CPL_TRUE;
            }

        } else {

            ok = CPL_FALSE;

            if (j == pixlim[1]) {
                if (jlim[0] == -1) {
                    sprintf(errtxt,
                            "%s: cpl_table *spec (no valid flux value)",
                            MF_ERROR_NDA_TXT);
                    status = cpl_error_set_message(cpl_func, MF_ERROR_NDA,
                                                   "%s", errtxt);
                } else {
                    /* Constant value for extrapolation of upper margin */
                    for (i = jlim[0] + 1; i <= pixlim[1]; i++) {
                        fluxv[i] = fluxv[jlim[0]];
                    }
                }
            }

        }

    }

    return status;
}

/**@}*/
