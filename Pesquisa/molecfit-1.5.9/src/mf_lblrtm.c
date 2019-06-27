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
 * \file mf_lblrtm.c
 *
 * Routines for running the radiative transfer code LBLRTM
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 * \since  20 Jun 2012
 * \date   26 Jul 2014
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#include "mf_molecfit.h"
#include <mf_lblrtm.h>


/*****************************************************************************
 *                                GLOBALS                                    *
 ****************************************************************************/

/* Declaration of global variables */

/* Number of LBLRTM calls (wavenumber-restricted subspectra are not
   counted individually) */
extern int n_code;
/* Total LBLRTM run time */
extern double t_code;


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code mf_lblrtm_call(const mfdrv *drvpar, const cpl_table *prof,
                              const int range)
{
    /*!
     * \brief
     *   Create a single LBLRTM sky model spectrum.
     *
     * This routine calls LBLRTM with the parameters specified in the driver
     * parameter structure \em drvpar. Wavelength range and resolution is
     * selected by the input parameter \em range. In addition, an atmospheric
     * profile \em prof is required as input.
     *
     * \b INPUT:
     * \param drvpar  driver parameter structure
     * \param prof    CPL table containing the atmospheric profiles
     * \param range   range number
     *
     * \b OUTPUT:
     * \return        CPL_ERROR_NONE on success, else:
     *                MF_ERROR_IO or MF_ERROR_SUBROUTINE.
     *
     * \b ERRORS:
     * - CPL_ERROR_FILE_NOT_CREATED:  Could not create working directory
     * - MF_ERROR_SUBROUTINE:  Errors in various subroutines (reading LBLRTM
     *                         setup, starting LBLRTM, renaming LBLRTM output)
     */

    char wdir[MF_MAXLEN],
         rmwdir[MF_MAXLEN],
         basedir[MF_MAXLEN],
         bindir[MF_MAXLEN],
         datadir[MF_MAXLEN],
         confdir[MF_MAXLEN],
         outdir[MF_MAXLEN],
         outfile[MF_MAXLEN];

    double res;     /* resolution */

    int j, jmax = 0, d;

    double wn1, wn2, min, max, minc, maxc, delta = 25., limlam[2] = {0, 0}, z;

    cpl_parameterlist *lblrtm_setup;
    cpl_parameter *p;

    cpl_errorstate err_state;

    /* get pixel resolution */
    res = cpl_table_get(drvpar->rangetab, "pixres", range-1, NULL)
        * MF_SAMPFAC;
    if (res > 1e6) {
        res = 1e6;
    }

    /* exit if range is empty */
    if (res == 0.) {
        return CPL_ERROR_NONE;
    }

//    printf("======================Pixel resolution: %lf\n", res);

    /* get base directory */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));

    /* get config directory */
    p = cpl_parameterlist_find(drvpar->parlist, "configdir");
    mf_basic_abspath(confdir, cpl_parameter_get_string(p), basedir);

    /* get data directory */
    p = cpl_parameterlist_find(drvpar->parlist, "datadir");
    // mf_basic_abspath(datadir, cpl_parameter_get_string(p), basedir);
    strcpy(datadir, cpl_parameter_get_string(p));

    /* get bin directory */
    p = cpl_parameterlist_find(drvpar->parlist, "bindir");
    mf_basic_abspath(bindir, cpl_parameter_get_string(p), basedir);

    /* get output directory */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);

    /* get zenith distance */
    p = cpl_parameterlist_find(drvpar->parlist, "telalt");
    z = 90. - cpl_parameter_get_double(p);

    /* get output filename */
    // p = cpl_parameterlist_find(drvpar->parlist, "spec_out");
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    sprintf(outfile, "%s_%d", cpl_parameter_get_string(p), range);

    /* create working directory */
    sprintf(wdir, "%swdir", outdir);
    err_state = cpl_errorstate_get();
    if (mf_basic_access(wdir, F_OK) != CPL_ERROR_NONE) {
        if (mkdir(wdir, 0755)) {
            return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                        "could not create working directory");
        }
    }
    cpl_errorstate_set(err_state);

    /* command for removal of working directory */
    sprintf(rmwdir, "rm -rf %s", wdir);

    /* read default lblrtm_setup */
    lblrtm_setup = cpl_parameterlist_new();
    if (mf_lblrtm_readsetup(lblrtm_setup, drvpar) != CPL_ERROR_NONE) {
        cpl_parameterlist_delete(lblrtm_setup);
        if ((d = chdir(outdir))){};
        if ((d = system(rmwdir))){};
        return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                     "could not read LBLRTM setup");
    }

    /* set lblrtm_setup parameter ANGLE */
    p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.ANGLE");
    cpl_parameter_set_double(p, z);

    /* get wavenumber range */
    wn1 = cpl_table_get(drvpar->rangetab, "wn_start", range-1, NULL);
    wn2 = cpl_table_get(drvpar->rangetab, "wn_end", range-1, NULL);

    min = wn2;
    for (j = 0, max = min; min > wn1; j++) {
        delta = 0.14 * pow(max, 1.1);
        if (delta > 1750.) {
            delta = 1750;
        }
        if (delta < 125.) {
            delta = 125;
        }
        min = max - delta;
        if (min < wn1) {
            min = wn1;
        }
        max = min;
        jmax = j;
    }

    min = wn2;
    for (j = jmax, max = min; min > wn1; j--) {
        delta = 0.14 * pow(max, 1.1);
        if (delta > 1750.) {
            delta = 1750;
        }
        if (delta < 125.) {
            delta = 125;
        }
        min = max - delta;
        if (min < wn1) {
            min = wn1;
        }
        if (j == jmax) {
            maxc = max;
        } else {
            maxc = max + MF_EXTRACOVER;
        }
        if (j == 0) {
            minc = min;
        } else {
            minc = min - MF_EXTRACOVER;
        }
//        printf("%i RANGE: %lf - %lf (%lf)\n",j,min,max,delta);
//        prof=NULL;
//        mf_lblrtm_renameoutput(wdir,0,'t');
        if (mf_lblrtm_start1(drvpar, lblrtm_setup, prof, minc, maxc, range) !=
                CPL_ERROR_NONE) {
            cpl_parameterlist_delete(lblrtm_setup);
            if ((d = chdir(outdir))){};
            if ((d = system(rmwdir))){};
            return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                         "LBLRTM failed for wavelength "
                                         "interval: %5.2f-%5.2f [µm]",
                                         MF_CONV_K_LAM/max,
                                         MF_CONV_K_LAM/min);
        }

        p = cpl_parameterlist_find(drvpar->parlist, "trans");
        if (cpl_parameter_get_int(p) != 0) {
            // Transmission
            if (mf_lblrtm_renameoutput(wdir, j+1, 't') != CPL_ERROR_NONE) {
                cpl_parameterlist_delete(lblrtm_setup);
                if ((d = chdir(outdir))){};
                if ((d = system(rmwdir))){};
                return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                             "Could not rename LBLRTM "
                                             "output");
            }
        } else {
            // Radiation
            if (mf_lblrtm_renameoutput(wdir, j+1, 'r') != CPL_ERROR_NONE) {
                cpl_parameterlist_delete(lblrtm_setup);
                if ((d = chdir(outdir))){};
                if ((d = system(rmwdir))){};
                return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                             "Could not rename LBLRTM "
                                             "output");
            }
        }
        max = min;
    }

    /* Number of mf_lblrtm_call calls */
    n_code++;

    /* clean exit */
    cpl_parameterlist_delete(lblrtm_setup);

    /* combine and convolve spectra */
    limlam[0] = MF_CONV_K_LAM/wn2;
    limlam[1] = MF_CONV_K_LAM/wn1;
    if (mf_lib_createlibspec(wdir, outdir, outfile, "fits", limlam, res,
                             drvpar) != CPL_ERROR_NONE) {
        if ((d = chdir(outdir))){};
        if ((d = system(rmwdir))){};
        return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                     "Error in mf_lib_createlibspec");
    }

//    printf("# run time mf_lib_createlibspec(): %10.5lf s\n",
//           cpl_test_get_cputime() - ts);

    /* remove working directory */
    if ((d = chdir(outdir))){};
    if ((d = system(rmwdir))){};

    return CPL_ERROR_NONE;
}


cpl_error_code mf_lblrtm_readsetup(cpl_parameterlist *plist,
                                   const mfdrv *drvpar)
{
    /*!
     * \brief
     *   Read LBLRTM setup file.
     *
     * This function reads a text file into a \c CPL_PARAMETERLIST containing
     * all parameters required to run LBLRTM.
     *
     * \b INPUT:
     * \param drvpar  driver parameter structure
     *
     * \b OUTPUT:
     * \param plist   extracted variable \c CPL_PARAMETERLIST
     * \return        CPL_ERROR_NONE on success,
     *                MF_ERROR_BADUSERINPUT or CPL_ERROR_FILE_IO else
     */

    FILE *fp;                        /* file pointer */

    char line[MF_MAXLEN];            /* line in setup file */

    char str[MF_MAXLEN],
         errtxt[MF_MAXLEN],
         *cptr;

    double tmp;

    int ICNTNM_success = 0,
        IAERSL_success = 0,
        MPTS_success = 0,
        NPTS_success = 0,
        V1_success = 0,
        V2_success = 0,
        SAMPLE_success = 0,
        ALFAL0_success = 0,
        AVMASS_success = 0,
        DPTMIN_success = 0,
        DPTFAC_success = 0,
        TBOUND_success = 0,
        SREMIS1_success = 0, SREMIS2_success = 0, SREMIS3_success = 0,
        SRREFL1_success = 0, SRREFL2_success = 0, SRREFL3_success = 0,
        MODEL_success = 0,
        ITYPE_success = 0,
        NOZERO_success = 0,
        NOPRNT_success = 0,
        IPUNCH_success = 0,
        RE_success = 0,
        HSPACE_success = 0,
        VBAR_success = 0,
        REF_LAT_success = 0,
        H1_success = 0,
        H2_success = 0,
        ANGLE_success = 0,
        RANGE_success = 0,
        BETA_success = 0,
        LEN_success = 0,
        HOBS_success = 0,
        AVTRAT_success = 0,
        TDIFF1_success = 0,
        TDIFF2_success = 0,
        ALTD1_success = 0,
        ALTD2_success = 0,
        DELV_success = 0;

    /* continua and Rayleigh extinction */
    int ICNTNM = 0, ICNTNM_default = ICNTNM;
    /*
     *  = 0  no continuum calculated
     *  = 1  all continua calculated,
     *       including Rayleigh extinction where applicable
     *  = 2  H2O self not calculated,
     *       all other continua/Rayleigh extinction calculated
     *  = 3  H2O foreign not calculated,
     *       all other continua/Rayleigh extinction calculated
     *  = 4  H2O self and foreign not calculated,
     *       all other continua/Rayleigh extinction calculated
     *  = 5  Rayleigh extinction not calculated, all other continua calculated
     *  = 6  Individual continuum scale factors input (Requires Record 1.2a)
     *
     *  [6 will not be implemented]
     */

    /* aerosols */
    int IAERSL = 0, IAERSL_default = IAERSL;
    /*
     *  = 0  no aerosols used
     *  = 1  internal LOWTRAN aerosol models
     *  = 5  spectral optical depths by layer from file 'in_lblrtm_cld'
     *  = 7  user defined aerosol models
     *  = 9  use precalculated aerosols (TAPE20 from a previous aerosol run)
     *
     *  [5,7,9 will not be implemented]
     */

    /* optical depth values */
    int MPTS = 5, MPTS_default = MPTS;
    /*
     *  number of optical depth values printed for the beginning and
     *  ending of each panel as a result of convolution for current layer
     *  (for MPTS < O, output printing is suppressed)
     */

    /* number of values for each panel */
    int NPTS = 5, NPTS_default = NPTS;
    /*
     *  number of values printed for the beginning and ending of each panel
     *  as result of merge of current layer with previous layers
     *  (optical depth for IEMIT=0; radiance and transmission for IEMIT=1)
     */

    /* beginning wavenumber value for the calculation */
    double V1 = 4167., V1_default = V1; /* -> ~2400nm */

    /* ending wavenumber value for the calculation
     * (V2-V1 must be less than 2020 cm-1) */
    double V2 = 5263., V2_default = V2; /* -> ~1900nm */

    /* number of sample points per mean halfwidth (between 1 and 4)
     * (default = 4) */
    int SAMPLE = 4, SAMPLE_default = SAMPLE;

    /* average collision broadened halfwidth (cm - 1/atm)
     * (default = 0.04) */
    double ALFAL0 = 0.04, ALFAL0_default = ALFAL0;

    /* average molecular mass (amu) for Doppler halfwidth
     * (default = 36) */
    double AVMASS = 36, AVMASS_default = AVMASS;

    /* minimum molecular optical depth below which lines will be rejected
     * (negative value defaults to DPTMIN = 0.0002) */
    double DPTMIN = 0.0002, DPTMIN_default = DPTMIN;

    /* factor multiplying molecular continuum optical depth to
     * determine optical depth below which lines will be rejected
     * (negative value defaults to DPTFAC = 0.001) */
    double DPTFAC = 0.001, DPTFAC_default = DPTFAC;

    /* temperature of boundary [K] */
    double TBOUND = 0, TBOUND_default = TBOUND;

    /* frequency dependent boundary emissivity coefficients
     * EMISSIVITY   = SREMIS1 + SREMIS2*V + SREMIS3*(V**2) */
    double SREMIS1 = 0, SREMIS2 = 0, SREMIS3 = 0, SREMIS1_default = SREMIS1,
        SREMIS2_default = SREMIS2, SREMIS3_default = SREMIS3;

    /* frequency dependent boundary reflectivity coefficients
     * REFLECTIVITY = SRREFL1 + SRREFL2*V + SRREFL3*(V**2) */
    double SRREFL1 = 0, SRREFL2 = 0, SRREFL3 = 0, SRREFL1_default = SRREFL1,
        SRREFL2_default = SRREFL2, SRREFL3_default = SRREFL3;

    /* selects atmospheric profile */
    int MODEL = 0, MODEL_default = MODEL;
    /*
     *  = 0  user supplied atmospheric profile
     *  = 1  tropical model
     *  = 2  midlatitude summer model
     *  = 3  midlatitude winter model
     *  = 4  subarctic summer model
     *  = 5  subarctic winter model
     *  = 6  U.S. standard 1976
     */

    /* selects type of path */
    int ITYPE = 3, ITYPE_default = ITYPE;
    /*
     *  = 1  horizontal path (constant pressure, temperature), use RECORD 3.2H
     *  = 2  slant path from H1 to H2, use RECORD 3.2
     *  = 3  slant path from H1 to space (see HSPACE), use RECORD 3.2
     */

    /* zeroing of small amounts of absorbers */
    int NOZERO = 0, NOZERO_default = NOZERO;
    /*
     *  = 0  zeroes absorber amounts which are less than 0.1 percent of total
     *       (default)
     *  = 1  suppresses zeroing of small amounts
     */

    /* output */
    int NOPRNT = 0, NOPRNT_default = NOPRNT;
    /*
     *  = 0  full printout
     *  = 1  selects short printout
     */

    /* write out layer data */
    int IPUNCH = 0, IPUNCH_default = IPUNCH;
    /*
     *  = 0  layer data not written (default)
     *  = 1  layer data written to unit IPU (TAPE7)
     */

    /* radius of earth [km] */
    double RE = 0, RE_default = RE;
    /*
     *  defaults for RE=0:
     *  a)  MODEL 0,2,3,6    RE = 6371.23 km
     *  b)        1          RE = 6378.39 km
     *  c)        4,5        RE = 6356.91 km
     */

    /* altitude definition for space (default = 100 km) */
    double HSPACE = 120, HSPACE_default = HSPACE;
    /*
     *  internal models defined to 120 km
     */

    /* frequency for refractive geometry calculation */
    double VBAR = 0, VBAR_default = VBAR;
    /*
     *  (default:  VBAR = (V1+V2) / 2 )     (V1,V2 from Record 1.3)
     */

    /* latitude of location of calculation [degrees] */
    double REF_LAT = 0, REF_LAT_default = REF_LAT;
    /*
     *  defaults for REF_LAT = 0:
     *  a) MODEL 0,2,3,6    REF_LAT = 45.0 degrees
     *  b) MODEL 1          REF_LAT = 15.0
     *  c) MODEL 4,5        REF_LAT = 60.0
     */

    /* observer altitude [km] */
    double H1 = 0, H1_default = H1;

    /* upper height limit */
    double H2 = 0, H2_default = H2;
    /*
     *  for ITYPE = 2, H2 is the end point altitude [km]
     *      ITYPE = 3, H2 is the tangent height [km] for H2 .GT. 0.
     *                 if H2 = 0. ANGLE determines tangent height
     */

    /* zenith angle at H1 [degrees] */
    double ANGLE = 0, ANGLE_default = ANGLE;

    /* length of a straight path from H1 to H2 [km] */
    double RANGE = 0, RANGE_default = RANGE;

    /* earth centered angle from H1 to H2 [degrees] */
    double BETA = 0, BETA_default = BETA;

    int LEN = 0, LEN_default = LEN;
    /*
     *  = 0  short path (default)
     *  = 1  long path through a tangent height
     *
     *  LEN is only used for H1 > H2 (ANGLE > 90`)
     *
     *  for ITYPE = 2, only 3 of the first 5 parameters are required to
     *                 specify the path, e.g., H1, H2, ANGLE or H1, H2 and
     *                 RANGElblrtm_setup
     *
     *  for ITYPE = 3, H1 = observer altitude must be specified. Either
     *                 H2 = tangent height or ANGLE must be specified.
     *                 Other parameters are ignored.
     */

    /* Height of observer */
    double HOBS = 0, HOBS_default = HOBS;
    /*
     *  Height of observer, used only for informational purposes in
     *  satellite-type simulations when computing output geometry
     *  above 120 km.
     */

    /* maximum Voigt width ratio across a layer (if zero, default=1.5) */
    double AVTRAT = 1.5, AVTRAT_default = AVTRAT;

    /* maximum layer temperature difference at ALTD1 (if zero, default=5 K) */
    double TDIFF1 = 5, TDIFF1_default = TDIFF1;

    /* maximum layer temperature difference at ALTD2 (if zero, default=8 K) */
    double TDIFF2 = 8, TDIFF2_default = TDIFF2;

    /* altitude of TDIFF1 (if zero, default = 0 Km) */
    double ALTD1 = 0, ALTD1_default = ALTD1;

    /* altitude of TDIFF2 (if zero, default = 100 Km) */
    double ALTD2 = 120, ALTD2_default = ALTD2;

    /* number of wavenumbers [cm-1] per major division. */
    double DELV = 1., DELV_default = DELV;

    cpl_parameter *p;

    char confdir[MF_MAXLEN],
         setup_file[MF_MAXLEN],
         basedir[MF_MAXLEN];

    /* get base directory */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));

    /* get config directory */
    p = cpl_parameterlist_find(drvpar->parlist, "configdir");

    char sharedir[MF_MAXLEN];
    sprintf(sharedir, "%s/..", mf_get_datadir());
    mf_basic_abspath(confdir, cpl_parameter_get_string(p), sharedir);

    /* lblrtm_setup filename */
    sprintf(setup_file, "%slblrtm_setup", confdir);

    /* open file for reading */
    fp = fopen(setup_file, "r");
    if (fp == NULL) {
        return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                     "Could not open lblrtm setup file: %s",
                                     setup_file);
    } else {
        mf_basic_initstring(line, MF_MAXLEN);

        /*
         *  Loop through lines in setup file.
         *  Each line may be empty, contain a comment (starting with "#") or
         *  contain a command. Lines with commands are parsed into the
         *  corresponding variables.
         */
        while (fgets(line, MF_MAXLEN - 1, fp) != NULL) { /* end of file? */
            mf_basic_terminatestring(line);
            mf_basic_strtrim_inplace(line);    /* trim line */
            mf_basic_rmcntrl_inplace(line);    /* remove control characters */

            if (line[0] == '#') {
                continue;                      /* skip over comments */
            }

            if (strlen(line) == 0) {
                continue;                      /* skip over empty lines */
            }

            /* extract variable identifier */
            mf_basic_initstring(str, MF_MAXLEN);
            cptr = strpbrk(line, " =");
            strncpy(str, line, cptr-line);

            /* continua and raleigh extinction */
            mf_basic_processrange(str, "ICNTNM", line, (void *)&ICNTNM,
                                  &ICNTNM_success, t_int, 0, 5);

            /* aerosols */
            mf_basic_processrange(str, "IAERSL", line, (void *)&IAERSL,
                                  &IAERSL_success, t_int, 0, 1);

            /* optical depth values */
            mf_basic_processvariable(str, "MPTS", line, (void *)&MPTS,
                                     &MPTS_success, t_int);

            /* number of values for each panel */
            mf_basic_processvariable(str, "NPTS", line, (void *)&NPTS,
                                     &NPTS_success, t_int);

            /* start wavelength [µm] (will be converted to wavenumber) */
            mf_basic_processvariable(str, "V1", line, (void *)&V1, &V1_success,
                                  t_double);

            /* end wavelength [µm] (will be converted to wavenumber) */
            mf_basic_processvariable(str, "V2", line, (void *)&V2, &V2_success,
                                  t_double);

            /* number of sample points */
            mf_basic_processrange(str, "SAMPLE", line, (void *)&SAMPLE,
                                  &SAMPLE_success, t_int, 1, 4);

            /* average collision broadened halfwidth [cm - 1/atm] */
            mf_basic_processvariable(str, "ALFAL0", line, (void *)&ALFAL0,
                                     &ALFAL0_success, t_double);

            /* average molecular mass [amu] */
            mf_basic_processvariable(str, "AVMASS", line, (void *)&AVMASS,
                                     &AVMASS_success, t_double);

            /* minimum molecular optical depth */
            mf_basic_processvariable(str, "DPTMIN", line, (void *)&DPTMIN,
                                     &DPTMIN_success, t_double);

            /* factor multiplying molecular continuum optical depth to
             * determine optical depth below which lines will be rejected */
            mf_basic_processvariable(str, "DPTFAC", line, (void *)&DPTFAC,
                                     &DPTFAC_success, t_double);

            /* temperature of boundary [K] */
            mf_basic_processvariable(str, "TBOUND", line, (void *)&TBOUND,
                                     &TBOUND_success, t_double);

            /* boundary emissivity coefficient 1 */
            mf_basic_processvariable(str, "SREMIS1", line, (void *)&SREMIS1,
                                     &SREMIS1_success, t_double);

            /* boundary emissivity coefficient 2 */
            mf_basic_processvariable(str, "SREMIS2", line, (void *)&SREMIS2,
                                     &SREMIS2_success, t_double);

            /* boundary emissivity coefficient 3 */
            mf_basic_processvariable(str, "SREMIS3", line, (void *)&SREMIS3,
                                     &SREMIS3_success, t_double);

            /* boundary reflectivity coefficient 1 */
            mf_basic_processvariable(str, "SRREFL1", line, (void *)&SRREFL1,
                                     &SRREFL1_success, t_double);

            /* boundary reflectivity coefficient 2 */
            mf_basic_processvariable(str, "SRREFL2", line, (void *)&SRREFL2,
                                     &SRREFL2_success, t_double);

            /* boundary reflectivity coefficient 3 */
            mf_basic_processvariable(str, "SRREFL3", line, (void *)&SRREFL3,
                                     &SRREFL3_success, t_double);

            /* atmospheric profile */
            mf_basic_processrange(str, "MODEL", line, (void *)&MODEL,
                                  &MODEL_success, t_int, 0, 6);

            /* type of path */
            mf_basic_processrange(str, "ITYPE", line, (void *)&ITYPE,
                                  &ITYPE_success, t_int, 1, 3);

            /* zeroing of small amounts of absorbers */
            mf_basic_processrange(str, "NOZERO", line, (void *)&NOZERO,
                                  &NOZERO_success, t_int, 0, 1);

            /* output */
            mf_basic_processrange(str, "NOPRNT", line, (void *)&NOPRNT,
                                  &NOPRNT_success, t_int, 0, 1);

            /* write out layer data */
            mf_basic_processrange(str, "IPUNCH", line, (void *)&IPUNCH,
                                  &IPUNCH_success, t_int, 0, 1);

            /* radius of earth [km] */
            mf_basic_processvariable(str, "RE", line, (void *)&RE,
                                     &RE_success, t_double);

            /* altitude definition for space */
            mf_basic_processvariable(str, "HSPACE", line, (void *)&HSPACE,
                                     &HSPACE_success, t_double);

            /* frequency for refractive geometry calculation */
            if (strcmp(str, "VBAR") == 0) {
                if (mf_basic_syntaxok(line) != CPL_TRUE &&
                    VBAR_success == 0) {
                    VBAR = (V1 + V2) / 2.;
                } else {
                    VBAR = atof(line);
                    VBAR_success = 1;
                }
            } /* frequency for refractive geometry calculation */

            /* latitude of location of calculation [degrees] */
            mf_basic_processrange(str, "REF_LAT", line, (void *)&REF_LAT,
                                  &REF_LAT_success, t_double, -90, 90);

            /* observer altitude [km] */
            mf_basic_processvariable(str, "H1", line, (void *)&H1,
                                     &H1_success, t_double);

            /* upper height limit */
            mf_basic_processvariable(str, "H2", line, (void *)&H2,
                                     &H2_success, t_double);

            /* zenith angle at H1 [degrees] */
            mf_basic_processvariable(str, "ANGLE", line, (void *)&ANGLE,
                                     &ANGLE_success, t_double);

            /* length of a straight path from H1 to H2 [km] */
            mf_basic_processvariable(str, "RANGE", line, (void *)&RANGE,
                                     &RANGE_success, t_double);

            /* earth centered angle from H14 to H2 [degrees] */
            mf_basic_processvariable(str, "BETA", line, (void *)&BETA,
                                     &BETA_success, t_double);

            /* path length */
            mf_basic_processrange(str, "LEN", line, (void *)&LEN,
                                  &LEN_success, t_int, 0, 1);

            /* height of observer */
            mf_basic_processvariable(str, "HOBS", line, (void *)&HOBS,
                                     &HOBS_success, t_double);

            /* maximum Voigt width ratio across a layer */
            mf_basic_processvariable(str, "AVTRAT", line, (void *)&AVTRAT,
                                     &AVTRAT_success, t_double);

            /* maximum layer temperature difference at ALTD1 */
            mf_basic_processvariable(str, "TDIFF1", line, (void *)&TDIFF1,
                                     &TDIFF1_success, t_double);

            /* maximum layer temperature difference at ALTD2 */
            mf_basic_processvariable(str, "TDIFF2", line, (void *)&TDIFF2,
                                     &TDIFF2_success, t_double);

            /* altitude of TDIFF1 */
            mf_basic_processvariable(str, "ALTD1", line, (void *)&ALTD1,
                                     &ALTD1_success, t_double);

            /* altitude of TDIFF2 */
            mf_basic_processvariable(str, "ALTD2", line, (void *)&ALTD2,
                                     &ALTD2_success, t_double);

            /* number of wavenumbers per major division */
            mf_basic_processvariable(str, "DELV", line, (void *)&DELV,
                                     &DELV_success, t_double);
        }

        fclose(fp);
    }

    /*
     *  ensure that V2-V1 is not too large and that V1/V2 are in the right
     *  order. Also convert them from wavelength to wavenumber
     */
    if (V1 < V2) {
        tmp = V1;
        V1 = V2;
        V2 = tmp;
    }
    V1 = MF_CONV_K_LAM / V1;
    V2 = MF_CONV_K_LAM / V2;
    if (V2 - V1 > 2020 || V1 == V2) {
        V1 = V1_default;
        V2 = V2_default;
        V1_success = 0;
        V2_success = 0;
        return cpl_error_set_message(cpl_func, MF_ERROR_BADUSERINPUT,
                                     "LBLRTM variable delta V: V2-V1 must be "
                                     "less than 2020cm-1");
    }

    /* initialise setup parameter list */

    /* continua and Rayleigh extinction */
    p = cpl_parameter_new_range("lblrtm_config.ICNTNM", CPL_TYPE_INT,
                                "continua and raleigh extinction "
                                "[0,1,2,3,4,5]", "lblrtm_config",
                                ICNTNM_default, 0, 5);
    cpl_parameter_set_int(p, ICNTNM);
    cpl_parameterlist_append(plist, p);
    /*
     *  = 0  no continuum calculated
     *  = 1  all continua calculated,
     *       including Rayleigh extinction where applicable
     *  = 2  H2O self not calculated,
     *       all other continua/Rayleigh extinction calculated
     *  = 3  H2O foreign not calculated,
     *       all other continua/Rayleigh extinction calculated
     *  = 4  H2O self and foreign not calculated,
     *       all other continua/Rayleigh extinction calculated
     *  = 5  Rayleigh extinction not calculated, all other continua calculated
     *  = 6  Individual continuum scale factors input (Requires Record 1.2a)
     *
     *  [6 will not be implemented]
     */

    /* aerosols */
    p = cpl_parameter_new_range("lblrtm_config.IAERSL", CPL_TYPE_INT,
                                "aerosols [0,1]", "lblrtm_config",
                                IAERSL_default, 0, 1);
    cpl_parameter_set_int(p, IAERSL);
    cpl_parameterlist_append(plist, p);
    /*
     *  = 0  no aerosols used
     *  = 1  internal LOWTRAN aerosol models
     *  = 5  spectral optical depths by layer from file 'in_lblrtm_cld'
     *  = 7  user defined aerosol models
     *  = 9  use precalculated aerosols (TAPE20 from a previous aerosol run)
     *
     *  [5,7,9 will not be implemented]
     */

    /* optical depth values */
    p = cpl_parameter_new_value("lblrtm_config.MPTS", CPL_TYPE_INT,
                                "number of optical depth values",
                                "lblrtm_config", MPTS_default);
    cpl_parameter_set_int(p, MPTS);
    cpl_parameterlist_append(plist, p);
    /*
     *  number of optical depth values printed for the beginning and
     *  ending of each panel as a result of convolution for current layer
     *  (for MPTS < O, output printing is suppressed)
     */

    /* number of values for each panel */
    p = cpl_parameter_new_value("lblrtm_config.NPTS", CPL_TYPE_INT,
                                "number of values for each panel",
                                "lblrtm_config", NPTS_default);
    cpl_parameter_set_int(p, NPTS);
    cpl_parameterlist_append(plist, p);
    /*
     *  number of values printed for the beginning and ending of each panel
     *  as result of merge of current layer with previous layers
     *  (optical depth for IEMIT=0; radiance and transmission for IEMIT=1)
     */

    /* beginning wavenumber value for the calculation */
    p = cpl_parameter_new_value("lblrtm_config.V1", CPL_TYPE_DOUBLE,
                                "beginning wavenumber value for the "
                                "calculation", "lblrtm_config", V1_default);
    cpl_parameter_set_double(p, V1);
    cpl_parameterlist_append(plist, p);

    /* ending wavenumber value for the calculation
     * (V2-V1 must be less than 2020 cm-1) */
    p = cpl_parameter_new_value("lblrtm_config.V2", CPL_TYPE_DOUBLE,
                                "ending wavenumber value for the "
                                "calculation", "lblrtm_config", V2_default);
    cpl_parameter_set_double(p, V2);
    cpl_parameterlist_append(plist, p);

    /* number of sample points per mean halfwidth (between 1 and 4)
     * (default = 4) */
    p = cpl_parameter_new_range("lblrtm_config.SAMPLE", CPL_TYPE_INT,
                                "number of sample points per mean halfwidth "
                                "[1.-4.]", "lblrtm_config", SAMPLE_default,
                                1., 4.);
    cpl_parameter_set_int(p, SAMPLE);
    cpl_parameterlist_append(plist, p);

    /* average collision broadened halfwidth [cm - 1/atm]
     * (default = 0.04) */
    p = cpl_parameter_new_value("lblrtm_config.ALFAL0", CPL_TYPE_DOUBLE,
                                "average collision broadened halfwidth "
                                "[cm-1/atm]", "lblrtm_config",
                                ALFAL0_default);
    cpl_parameter_set_double(p, ALFAL0);
    cpl_parameterlist_append(plist, p);

    /* average molecular mass [amu] for Doppler halfwidth
     * (default = 36) */
    p = cpl_parameter_new_value("lblrtm_config.AVMASS", CPL_TYPE_DOUBLE,
                                "average molecular mass [amu] for Doppler "
                                "halfwidth", "lblrtm_config", AVMASS_default);
    cpl_parameter_set_double(p, AVMASS);
    cpl_parameterlist_append(plist, p);

    /* minimum molecular optical depth below which lines will be rejected
     * (negative value defaults to DPTMIN = 0.0002) */
    p = cpl_parameter_new_value("lblrtm_config.DPTMIN", CPL_TYPE_DOUBLE,
                                "minimum molecular optical depth below which "
                                "lines will be rejected", "lblrtm_config",
                                DPTMIN_default);
    cpl_parameter_set_double(p, DPTMIN);
    cpl_parameterlist_append(plist, p);

    /* factor multiplying molecular continuum optical depth to
     * determine optical depth below which lines will be rejected
     * (negative value defaults to DPTFAC = 0.001) */
    p = cpl_parameter_new_value("lblrtm_config.DPTFAC", CPL_TYPE_DOUBLE,
                                "factor multiplying molecular continuum "
                                "optical depth", "lblrtm_config",
                                DPTFAC_default);
    cpl_parameter_set_double(p, DPTFAC);
    cpl_parameterlist_append(plist, p);

    /* temperature of boundary [K] */
    p = cpl_parameter_new_value("lblrtm_config.TBOUND", CPL_TYPE_DOUBLE,
                                "temperature of boundary [K]",
                                "lblrtm_config", TBOUND_default);
    cpl_parameter_set_double(p, TBOUND);
    cpl_parameterlist_append(plist, p);

    /* frequency dependent boundary emissivity coefficients
     * EMISSIVITY   = SREMIS1 + SREMIS2*V + SREMIS3*(V**2) */
    p = cpl_parameter_new_value("lblrtm_config.SREMIS1", CPL_TYPE_DOUBLE,
                                "emissivity coefficient 1", "lblrtm_config",
                                SREMIS1_default);
    cpl_parameter_set_double(p, SREMIS1);
    cpl_parameterlist_append(plist, p);
    p = cpl_parameter_new_value("lblrtm_config.SREMIS2", CPL_TYPE_DOUBLE,
                                "emissivity coefficient 2", "lblrtm_config",
                                SREMIS2_default);
    cpl_parameter_set_double(p, SREMIS2);
    cpl_parameterlist_append(plist, p);
    p = cpl_parameter_new_value("lblrtm_config.SREMIS3", CPL_TYPE_DOUBLE,
                                "emissivity coefficient 3", "lblrtm_config",
                                SREMIS3_default);
    cpl_parameter_set_double(p, SREMIS3);
    cpl_parameterlist_append(plist, p);

    /* frequency dependent boundary reflectivity coefficients
     * REFLECTIVITY = SRREFL1 + SRREFL2*V + SRREFL3*(V**2) */
    p = cpl_parameter_new_value("lblrtm_config.SRREFL1", CPL_TYPE_DOUBLE,
                                "reflectivity coefficient 1",
                                "lblrtm_config", SRREFL1_default);
    cpl_parameter_set_double(p, SRREFL1);
    cpl_parameterlist_append(plist, p);
    p = cpl_parameter_new_value("lblrtm_config.SRREFL2", CPL_TYPE_DOUBLE,
                                "reflectivity coefficient 2",
                                "lblrtm_config", SRREFL2_default);
    cpl_parameter_set_double(p, SRREFL2);
    cpl_parameterlist_append(plist, p);
    p = cpl_parameter_new_value("lblrtm_config.SRREFL3", CPL_TYPE_DOUBLE,
                                "reflectivity coefficient 3",
                                "lblrtm_config", SRREFL3_default);
    cpl_parameter_set_double(p, SRREFL3);
    cpl_parameterlist_append(plist, p);

    /* selects atmospheric profile */
    p = cpl_parameter_new_range("lblrtm_config.MODEL", CPL_TYPE_INT,
                                "atmospheric profile [0,1,2,3,4,5,6]",
                                "lblrtm_config", MODEL_default, 0, 6);
    cpl_parameter_set_int(p, MODEL);
    cpl_parameterlist_append(plist, p);
    /*
     *  = 0  user supplied atmospheric profile
     *  = 1  tropical model
     *  = 2  midlatitude summer model
     *  = 3  midlatitude winter model
     *  = 4  subarctic summer model
     *  = 5  subarctic winter model
     *  = 6  U.S. standard 1976
     */

    /* selects type of path */
    p = cpl_parameter_new_range("lblrtm_config.ITYPE", CPL_TYPE_INT,
                                "type of path [1,2,3]", "lblrtm_config",
                                ITYPE_default, 1, 3);
    cpl_parameter_set_int(p, ITYPE);
    cpl_parameterlist_append(plist, p);
    /*
     *  = 1  horizontal path (constant pressure, temperature), use RECORD 3.2H
     *  = 2  slant path from H1 to H2, use RECORD 3.2
     *  = 3  slant path from H1 to space (see HSPACE), use RECORD 3.2
     */

    /* zeroing of small amounts of absorbers */
    p = cpl_parameter_new_range("lblrtm_config.NOZERO", CPL_TYPE_INT,
                                "zeroing of small amounts of absorbers [0,1]",
                                "lblrtm_config", NOZERO_default, 0, 1);
    cpl_parameter_set_int(p, NOZERO);
    cpl_parameterlist_append(plist, p);
    /*
     *  = 0  zeroes absorber amounts which are less than 0.1 percent of total
     *       (default)
     *  = 1  suppresses zeroing of small amounts
     */

    /* output */
    p = cpl_parameter_new_range("lblrtm_config.NOPRNT", CPL_TYPE_INT,
                                "do not print output? [0,1]", "lblrtm_config",
                                NOPRNT_default, 0, 1);
    cpl_parameter_set_int(p, NOPRNT);
    cpl_parameterlist_append(plist, p);
    /*
     *  = 0  full printout
     *  = 1  selects short printout
     */

    /* write out layer data */
    p = cpl_parameter_new_range("lblrtm_config.IPUNCH", CPL_TYPE_INT,
                                "write out layer data to TAPE7 [0,1]",
                                "lblrtm_config", IPUNCH_default, 0, 1);
    cpl_parameter_set_int(p, IPUNCH);
    cpl_parameterlist_append(plist, p);
    /*
     *  = 0  layer data not written (default)
     *  = 1  layer data written to unit ITAPE7)PU (TAPE7)
     */


    /* radius of earth [km] */
    p = cpl_parameter_new_value("lblrtm_config.RE", CPL_TYPE_DOUBLE,
                                "radius of earth [km]", "lblrtm_config",
                                RE_default);
    cpl_parameter_set_double(p, RE);
    cpl_parameterlist_append(plist, p);
    /*
     *  defaults for RE=0:
     *  a)  MODEL 0,2,3,6    RE = 6371.23 km
     *  b)        1          RE = 6378.39 km
     *  c)        4,5        RE = 6356.91 km
     */

    /* altitude definition for space (default = 100 km) */
    p = cpl_parameter_new_value("lblrtm_config.HSPACE", CPL_TYPE_DOUBLE,
                                "altitude definition for space [km]",
                                "lblrtm_config", HSPACE_default);
    cpl_parameter_set_double(p, HSPACE);
    cpl_parameterlist_append(plist, p);
    /*
     *  internal models defined to 120 km
     */

    /* frequency for refractive geometry calculation */
    p = cpl_parameter_new_value("lblrtm_config.VBAR", CPL_TYPE_DOUBLE,
                                "frequency for refractive geometry "
                                "calculation", "lblrtm_config", VBAR_default);
    cpl_parameter_set_double(p, VBAR);
    cpl_parameterlist_append(plist, p);
    /*
     *  (default:  VBAR = (V1+V2) / 2 )     (V1,V2 from Record 1.3)
     */

    /* latitude of location of calculation [degrees] */
    p = cpl_parameter_new_range("lblrtm_config.REF_LAT", CPL_TYPE_DOUBLE,
                                "latitude of location of calculation "
                                "[degrees] [-90.-90]", "lblrtm_config",
                                REF_LAT_default, -90., 90.);
    cpl_parameter_set_double(p, REF_LAT);
    cpl_parameterlist_append(plist, p);
    /*
     *  defaults for REF_LAT = 0:
     *  a) MODEL 0,2,3,6    REF_LAT = 45.0 degrees
     *  b) MODEL 1          REF_LAT = 15.0
     *  c) MODEL 4,5        REF_LAT = 60.0
     */

    /* observer altitude [km] */
    p = cpl_parameter_new_value("lblrtm_config.H1", CPL_TYPE_DOUBLE,
                                "observer altitude [km]", "lblrtm_config",
                                H1_default);
    cpl_parameter_set_double(p, H1);
    cpl_parameterlist_append(plist, p);

    /* upper height limit */
    p = cpl_parameter_new_value("lblrtm_config.H2", CPL_TYPE_DOUBLE,
                                "upper height limit [km]", "lblrtm_config",
                                H2_default);
    cpl_parameter_set_double(p, H2);
    cpl_parameterlist_append(plist, p);
    /*
     *  for ITYPE = 2, H2 is the end point altitude [km]
     *      ITYPE = 3, H2 is the tangent height [km] for H2 .GT. 0.
     *                 if H2 = 0. ANGLE determines tangent height
     */

    /* zenith angle at H1 [degrees] */
    p = cpl_parameter_new_range("lblrtm_config.ANGLE", CPL_TYPE_DOUBLE,
                                "zenith angle at H1 [degrees] [0.-90.]",
                                "lblrtm_config", ANGLE_default, 0., 90.);
    cpl_parameter_set_double(p, ANGLE);
    cpl_parameterlist_append(plist, p);

    /* length of a straight path from H1 to H2 [km] */
    p = cpl_parameter_new_value("lblrtm_config.RANGE", CPL_TYPE_DOUBLE,
                                "length of a straight path from H1 to H2 "
                                "[km]", "lblrtm_config", RANGE_default);
    cpl_parameter_set_double(p, RANGE);
    cpl_parameterlist_append(plist, p);

    /* earth centered angle from H1 to H2 [degrees] */
    p = cpl_parameter_new_value("lblrtm_config.BETA", CPL_TYPE_DOUBLE,
                                "earth centered angle from H1 to H2 "
                                "[degrees]", "lblrtm_config", BETA_default);
    cpl_parameter_set_double(p, BETA);
    cpl_parameterlist_append(plist, p);

    /* path length */
    p = cpl_parameter_new_range("lblrtm_config.LEN", CPL_TYPE_INT,
                                "path length [0,1]",
                                "lblrtm_config", LEN_default, 0, 1);
    cpl_parameter_set_int(p, LEN);
    cpl_parameterlist_append(plist, p);
    /*
     *  = 0  short path (default)
     *  = 1  long path through a tangent height
     *
     *  LEN is only used for H1 > H2 (ANGLE > 90`)
     *
     *  for ITYPE = 2, only 3 of the first 5 parameters are required to
     *                 specify the path, e.g., H1, H2, ANGLE or H1, H2 and
     *                 RANGE
     *
     *  for ITYPE = 3, H1 = observer altitude must be specified. Either
     *                 H2 = tangent height or ANGLE must be specified.
     *                 Other parameters are ignored.
     */

    /* Height of observer */
    p = cpl_parameter_new_value("lblrtm_config.HOBS", CPL_TYPE_DOUBLE,
                                "Height of observer", "lblrtm_config",
                                HOBS_default);
    cpl_parameter_set_double(p, HOBS);
    cpl_parameterlist_append(plist, p);
    /*
     *  Height of observer, used only for informational purposes in
     *  satellite-type simulations when computing output geometry
     *  above 120 km.
     */

    /* maximum Voigt width ratio across a layer (if zero, default=1.5) */
    p = cpl_parameter_new_value("lblrtm_config.AVTRAT", CPL_TYPE_DOUBLE,
                                "maximum Voigt width ratio across a layer",
                                "lblrtm_config", AVTRAT_default);
    cpl_parameter_set_double(p, AVTRAT);
    cpl_parameterlist_append(plist, p);

    /* maximum layer temperature difference at ALTD1 (if zero, default=5 K) */
    p = cpl_parameter_new_value("lblrtm_config.TDIFF1", CPL_TYPE_DOUBLE,
                                "maximum layer temperature difference at "
                                "ALTD1 [K]", "lblrtm_config", TDIFF1_default);
    cpl_parameter_set_double(p, TDIFF1);
    cpl_parameterlist_append(plist, p);

    /* maximum layer temperature difference at ALTD2 (if zero, default=8 K) */
    p = cpl_parameter_new_value("lblrtm_config.TDIFF2", CPL_TYPE_DOUBLE,
                                "maximum layer temperature difference at "
                                "ALTD1 [K]", "lblrtm_config", TDIFF2_default);
    cpl_parameter_set_double(p, TDIFF2);
    cpl_parameterlist_append(plist, p);

    /* altitude of TDIFF1 (if zero, default = 0 Km) */
    p = cpl_parameter_new_value("lblrtm_config.ALTD1", CPL_TYPE_DOUBLE,
                                "altitude of TDIFF [km]", "lblrtm_config",
                                ALTD1_default);
    cpl_parameter_set_double(p, ALTD1);
    cpl_parameterlist_append(plist, p);

    /* altitude of TDIFF2 (if zero, default = 100 Km) */
    p = cpl_parameter_new_value("lblrtm_config.ALTD2", CPL_TYPE_DOUBLE,
                                "altitude of TDIFF2 [km]", "lblrtm_config",
                                ALTD2_default);
    cpl_parameter_set_double(p, ALTD2);
    cpl_parameterlist_append(plist, p);

    /* number of wavenumbers [cm-1] per major division */
    p = cpl_parameter_new_value("lblrtm_config.DELV", CPL_TYPE_DOUBLE,
                                "number of wavenumbers [cm-1] per major "
                                "division", "lblrtm_config", DELV_default);
    cpl_parameter_set_double(p, DELV);
    cpl_parameterlist_append(plist, p);

//    if (VBAR_success == 0) {
//        cpl_msg_info(cpl_func, "Default value taken for VBAR");
//    }
    if (ICNTNM_success == 0 ||
        IAERSL_success == 0 || MPTS_success == 0 || NPTS_success == 0 ||
        V1_success == 0 || V2_success == 0 || SAMPLE_success == 0 ||
        ALFAL0_success == 0 || AVMASS_success == 0 || DPTMIN_success == 0 ||
        DPTFAC_success == 0 || TBOUND_success == 0 || SREMIS1_success == 0 ||
        SREMIS2_success == 0 || SREMIS3_success == 0 ||
        SRREFL1_success == 0 || SRREFL2_success == 0 ||
        SRREFL3_success == 0 || MODEL_success == 0 || ITYPE_success == 0 ||
        NOZERO_success == 0 || NOPRNT_success == 0 || IPUNCH_success == 0 ||
        RE_success == 0 || HSPACE_success == 0 ||
        REF_LAT_success == 0 || H1_success == 0 || H2_success == 0 ||
        ANGLE_success == 0 || RANGE_success == 0 || BETA_success == 0 ||
        LEN_success == 0 || HOBS_success == 0 || AVTRAT_success == 0 ||
        TDIFF1_success == 0 || TDIFF2_success == 0 || ALTD1_success == 0 ||
        ALTD2_success == 0 || DELV_success == 0) {

        /* check whether all parameters were successfully read */
        if (ICNTNM_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: ICNTNM");
        }
        if (IAERSL_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: IAERSL");
        }
        if (MPTS_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: MPTS");
        }
        if (NPTS_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: NPTS");
        }
        if (V1_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: V1");
        }
        if (V2_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: V2");
        }
        if (SAMPLE_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: SAMPLE");
        }
        if (ALFAL0_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: ALFAL0");
        }
        if (AVMASS_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: AVMASS");
        }
        if (DPTMIN_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: DPTMIN");
        }
        if (DPTFAC_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: DPTFAC");
        }
        if (TBOUND_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: TBOUND");
        }
        if (SREMIS1_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: SREMIS1");
        }
        if (SREMIS2_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: SREMIS2");
        }
        if (SREMIS3_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: SREMIS3");
        }
        if (SRREFL1_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: SRREFL1");
        }
        if (SRREFL2_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: SRREFL2");
        }
        if (SRREFL3_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: SRREFL3");
        }
        if (MODEL_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: MODEL");
        }
        if (ITYPE_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: ITYPE");
        }
        if (NOZERO_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: NOZERO");
        }
        if (NOPRNT_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: NOPRNT");
        }
        if (IPUNCH_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: IPUNCH");
        }
        if (RE_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: RE");
        }
        if (HSPACE_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: HSPACE");
        }
        if (REF_LAT_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: REF_LAT");
        }
        if (H1_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: H1");
        }
        if (H2_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: H2");
        }
        if (ANGLE_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: ANGLE");
        }
        if (RANGE_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: RANGE");
        }
        if (BETA_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: BETA");
        }
        if (LEN_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: LEN");
        }
        if (HOBS_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: HOBS");
        }
        if (AVTRAT_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: AVTRAT");
        }
        if (TDIFF1_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: TDIFF1");
        }
        if (TDIFF2_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: TDIFF2");
        }
        if (ALTD1_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: ALTD1");
        }
        if (ALTD2_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: ALTD2");
        }
        if (DELV_success == 0) {
            sprintf(errtxt, "syntax error in LBLRTM setup file: DELV");
        }

        return cpl_error_set_message(cpl_func, MF_ERROR_BADUSERINPUT,
                                     "%s", errtxt);
    } else {
        return CPL_ERROR_NONE;
    }
}


cpl_error_code mf_lblrtm_start1(const mfdrv *drvpar,
                                cpl_parameterlist *lblrtm_setup,
                                const cpl_table *prof, double V1, double V2,
                                const int range)
{
    /*!
     * \brief
     *   Create TAPE5 file and start LBLRTM for a single wavelength interval.
     *
     * This function parses the \c CPL_PARAMETERLIST \em lblrtm_setup and
     * creates a TAPE5 file in the working directory (specified in
     * \em lblrtm_setup) using the atmospheric profile \em prof. Then it
     * starts LBLRTM for the wavenumber interval \em V1 - \em V2. Additional
     * parameters are taken from the driver file \em drvpar.
     *
     * \b INPUT:
     * \param drvpar        ::mfdrv parameter structure
     * \param lblrtm_setup  LBLRTM setup \c CPL_PARAMETERLIST
     * \param prof          atmospheric profile
     * \param V1            starting wavenumber
     * \param V2            ending wavenumber
     * \param range         range number
     *
     * \b OUTPUT:
     * \return              CPL_ERROR_NONE on success,
     *                      MF_ERROR_BADUSERINPUT,
     *                      MF_ERROR_CD,
     *                      CPL_ERROR_FILE_NOT_FOUND,
     *                      CPL_ERROR_FILE_NOT_CREATED,
     *                      MF_ERROR_FOPEN,
     *                      CPL_ERROR_FILE_IO
     *                      MF_ERROR_SUBROUTINE else
     */

    cpl_parameter *p;

    cpl_errorstate err_state;

    char errtxt[MF_MAXLEN];

    int int1, int2, int3, int4, d;

    double dbl1, dbl2, dbl3, dbl4, dbl5, dbl6, dbl7;

    /* atmospheric model */
    int MODEL = 0;

    /* number of atmospheric profile boundaries */
    int immax = 0; /* number taken from length of merged profile */

    cpl_array *atm_molecs = NULL,     /* array for all LBLRTM molecules */
              *allmolecs = NULL;
    FILE *fp;                         /* file pointer for output TAPE5 */

    char sys[MF_MAXLEN],           /* string for various system commands */
         basedir[MF_MAXLEN],
         bindir[MF_MAXLEN],
         confdir[MF_MAXLEN],
         outdir[MF_MAXLEN],
         wdir[MF_MAXLEN],
         tape3[MF_MAXLEN];

    const char *molec;                /* string for a single molecule */
    const char *molec_string;         /* molecule list */

    int level,                        /* counter for current height level */
        allmol,                       /* counter for LBLRTM molecules */
        last_mol,
        mol,                          /* counter for molecules in list */
        success=1;                    /* flag recording successful writing of
                                         TAPE5 file */

    double tw;

//    p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.MODEL");
//    MODEL = cpl_parameter_get_int(p);

    /* get base directory */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));

    /* get bin directory */
    p = cpl_parameterlist_find(drvpar->parlist, "bindir");
    mf_basic_abspath(bindir, cpl_parameter_get_string(p), basedir);

    /* get config directory */
    p = cpl_parameterlist_find(drvpar->parlist, "configdir");

    char sharedir[MF_MAXLEN];
    sprintf(sharedir, "%s/..", mf_get_datadir());
    mf_basic_abspath(confdir, cpl_parameter_get_string(p), sharedir);

    /* get output directory */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);

    /* working directory */
    sprintf(wdir, "%swdir", outdir);

    /* TAPE3 file */
    sprintf(tape3, "%s", cpl_table_get_string(drvpar->rangetab, "lnfl",
                                              range-1));

    if (V1 > 0 && V2 - V1 > 2020) {
        return cpl_error_set_message(cpl_func, MF_ERROR_BADUSERINPUT,
                                     "|V2 - V1| > 2020 cm-1");;
    }

    /* change working directory */
    if (chdir(wdir)){
        return cpl_error_set_message(cpl_func, MF_ERROR_CD, MF_ERROR_CD_TXT);
    }

    /* remove existing TAPE3 in working directory
     * unless to be used TAPE3 is in working directory */
    sprintf(sys, "%s/TAPE3", wdir);
    err_state = cpl_errorstate_get();
    if (strcmp(sys, tape3) != 0 &&
            mf_basic_access(sys, F_OK) == CPL_ERROR_NONE) {
        if ((d = system("rm TAPE3"))){};
    }

    /* create symbolic link to TAPE3 in working directory */
    if (mf_basic_access(tape3, F_OK) != CPL_ERROR_NONE) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                                     "Could not find TAPE3");
    }
    if (mf_basic_access("TAPE3", F_OK) != CPL_ERROR_NONE) {
        sprintf(sys, "ln -s %s TAPE3", tape3);
        if (system(sys)) {
            return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                         "Could not create link to TAPE3");
        }
    }

    /* remove existing TAPE5, TAPE6, TAPE9, TAPE10, TAPE11, TAPE12, TAPE27 &
     * TAPE28 */
    sprintf(sys, "%s/TAPE5", wdir);
    if (mf_basic_access(sys, F_OK) == CPL_ERROR_NONE) {
        remove("TAPE5");
    }
    sprintf(sys, "%s/TAPE6", wdir);
    if (mf_basic_access(sys, F_OK) == CPL_ERROR_NONE) {
        remove("TAPE6");
    }
    sprintf(sys, "%s/TAPE9", wdir);
    if (mf_basic_access(sys, F_OK) == CPL_ERROR_NONE) {
        remove("TAPE9");
    }
    sprintf(sys, "%s/TAPE10", wdir);
    if (mf_basic_access(sys, F_OK) == CPL_ERROR_NONE) {
        remove("TAPE10");
    }
    sprintf(sys, "%s/TAPE11", wdir);
    if (mf_basic_access(sys, F_OK) == CPL_ERROR_NONE) {
        remove("TAPE11");
    }
    sprintf(sys, "%s/TAPE12", wdir);
    if (mf_basic_access(sys, F_OK) == CPL_ERROR_NONE) {
        remove("TAPE12");
    }
    sprintf(sys, "%s/TAPE27", wdir);
    if (mf_basic_access(sys, F_OK) == CPL_ERROR_NONE) {
        remove("TAPE27");
    }
    sprintf(sys, "%s/TAPE28", wdir);
    if (mf_basic_access(sys, F_OK) == CPL_ERROR_NONE) {
        remove("TAPE28");
    }
    cpl_errorstate_set(err_state);

    /* create TAPE5 file */
    fp = fopen("TAPE5", "w");
    if (fp == NULL) {
        return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                     "Could not open TAPE5 for writing");
    } else {
        /* get names of molecules in atmospheric profile*/
        atm_molecs = cpl_table_get_column_names(prof);

        /* store number of levels of merged profile and number of LBLRTM
         * molecules */
        immax = cpl_table_get_nrow(prof);

        /* get string with fit molecules */
        p = cpl_parameterlist_find(drvpar->parlist, "lbl_molecs");
        molec_string = cpl_parameter_get_string(p);

        /* get all molecule identifiers */
        allmolecs = mf_lblrtm_allmolecs();

        /* find last molecule in molec_string that is activated */
        last_mol = strrchr(molec_string, '1') - molec_string;

        /* LBLRTM record 1.1 */
        if (fputs("$ created by lblrtm_start1\n", fp) == EOF) {
            success = 0;
        }

        /* LBLRTM record 1.2 */
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.ICNTNM");
        int1 = cpl_parameter_get_int(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.IAERSL");
        int2 = cpl_parameter_get_int(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.MPTS");
        int3 = cpl_parameter_get_int(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.NPTS");
        int4 = cpl_parameter_get_int(p);
        // 4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1, 4X,I1, 4X,I1, 4X,I1, 3X,A2, 4X,I1, 4X,I1,  4X,I1, 1X,I4, 1X,I4,  4X,I1,4X,I1
        // 4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1, 4X,I1, 4X,I1, 4X,I1, 3X,A2, 4X,I1, 4X,I1,  4X,I1, 1X,I4, 1X,I4
        if (fprintf(fp, "    %1i    %1i    %1i    %1i    %1i    %1i    %1i"
                    "    %1i    %1i    %1i   %2s    %1i    %1i    %1i %4i %4i    %1i    %1i"
                    "\n", 1, 1, int1, int2, 1, 0, 0, 0, 0, 1, " 0", 0, 0,
                    0, int3, int4, 0, 0) < 0) {
            success = 0;
        }

        /* LBLRTM record 1.3 */
        // E10.3,  E10.3,    E10.3,   E10.3,   E10.3,    E10.3,    E10.3,    E10.3,    4X,I1,  5X,E10.3,       3x,I2
        // E10.3,  E10.3,    E10.3,   E10.3,   E10.3,    E10.3,    E10.3,    E10.3,    4X,I1,  5X,E10.3,       3x,I2
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.SAMPLE");
        int1 = cpl_parameter_get_int(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.ALFAL0");
        dbl4 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.AVMASS");
        dbl5 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.DPTMIN");
        dbl6 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.DPTFAC");
        dbl7 = cpl_parameter_get_double(p);
        if (fprintf(fp, "%10.3e%10.3e%10i%10.3e%10.3e%10.3e%10.3e%10.3e    "
                    "%1i     %10s   %2s\n", V1, V2, int1, 0., dbl4, dbl5,
                    dbl6, dbl7, 0, "", "") < 0) {
            success = 0;
        }

        /* LBLRTM record 1.4 */
        // E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3    4X,1A
        // E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3    4X,1A
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.TBOUND");
        dbl1 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.SREMIS1");
        dbl2 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.SREMIS2");
        dbl3 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.SREMIS3");
        dbl4 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.SRREFL1");
        dbl5 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.SRREFL2");
        dbl6 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.SRREFL3");
        dbl7 = cpl_parameter_get_double(p);
        if (fprintf(fp, "%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e    %c\n",
                    dbl1, dbl2, dbl3, dbl4, dbl5, dbl6, dbl7, 's') < 0) {
            success = 0;
        }

        /* LBLRTM record 3.1 */
        // I5,     I5,    I5,      I5,      I5,    I5,     I5,     I2,   1X, I2, F10.3,  F10.3, F10.3,   10x, F10.3
        // I5,     I5,    I5,      I5,      I5,    I5,     I5,     I2,   1X, I2, F10.3,  F10.3, F10.3,   10x, F10.3
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.ITYPE");
        int1 = cpl_parameter_get_int(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.NOZERO");
        int2 = cpl_parameter_get_int(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.NOPRNT");
        int3 = cpl_parameter_get_int(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.IPUNCH");
        int4 = cpl_parameter_get_int(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.RE");
        dbl1 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.HSPACE");
        dbl2 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.VBAR");
        dbl3 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.REF_LAT");
        dbl4 = cpl_parameter_get_double(p);
        if (fprintf(fp, "%5i%5i%5i%5i%5i%5i%5i%2i %2i%10.3e%10.3e%10.3e"
                    "          %10.3e\n", MODEL, int1, 0, int2, int3,
                    last_mol+1, int4, 0, 0, dbl1, dbl2, dbl3, dbl4) < 0) {
            success = 0;
        }

        /* LBLRTM record 3.2 */
        // F10.3, F10.3,   F10.3,   F10.3,  F10.3,    I5, 5X,F10.3
        // F10.3, F10.3,   F10.3,   F10.3,  F10.3,    I5, 5X,F10.3
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.H1");
        dbl1 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.H2");
        dbl2 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.ANGLE");
        dbl3 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.RANGE");
        dbl4 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.BETA");
        dbl5 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.LEN");
        int1 = cpl_parameter_get_int(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.HOBS");
        dbl6 = cpl_parameter_get_double(p);
        if (fprintf(fp, "%10.3e%10.3e%10.3e%10.3e%10.3e%5i     %10.3e\n",
                    dbl1, dbl2, dbl3, dbl4, dbl5, int1, dbl6) < 0) {
            success = 0;
        }

        /* LBLRTM record 3.3a */
        // F10.3,  F10.3,  F10.3, F10.3, F10.3
        // F10.3,  F10.3,  F10.3, F10.3, F10.3
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.AVTRAT");
        dbl1 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.TDIFF1");
        dbl2 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.TDIFF2");
        dbl3 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.ALTD1");
        dbl4 = cpl_parameter_get_double(p);
        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.ALTD2");
        dbl5 = cpl_parameter_get_double(p);
        if (fprintf(fp, "%10.3e%10.3e%10.3e%10.3e%10.3e\n",
                    dbl1, dbl2, dbl3, dbl4, dbl5)   < 0) {
            success = 0;
        }

        /* LBLRTM record 3.4 */
        // I5,    3A8
        // I5,    3A8
        if (fprintf(fp, "%5i%8s%8s%8s\n", immax, "", "", "") < 0) {
            success = 0;
        }

        /* loop over all levels */
        for (level = 0; level < immax; level++) {
            /* LBLRTM record 3.5 */
            // E10.3, E10.3, E10.3,   5x,  A1,     A1,  1x, A1,     1x,    39A1
            // E10.3, E10.3, E10.3,   5x,  A1,     A1,  1x, A1,     1x,    39A1
            if (fprintf(fp, "%10.3e%10.3e%10.3e     AA   ",
                        cpl_table_get_float(prof, "HGT", level, NULL),
                        cpl_table_get_float(prof, "PRE", level, NULL),
                        cpl_table_get_float(prof, "TEM", level, NULL)) < 0) {
                success = 0;
            }

            /* loop over required molecules in molec_string/allmolecs */
            for (allmol = 0, mol = 0; allmol <= last_mol; allmol++) {
                molec = cpl_array_get_string(allmolecs, allmol);
                fputc('A', fp);
            }

            if (fputc('\n', fp) == EOF) {
                success = 0;
            }

            /* LBLRTM record 3.6.1 ... 3.6.n */
            /* loop over all LBLRTM molecules */
            // 8E10.3
            // 8E10.3
            for (allmol = 0, mol = 0; allmol <= last_mol; allmol++) {
                molec = cpl_array_get_string(allmolecs, allmol);
                /* either write out value from table or 0. otherwise */
                if (cpl_table_has_column(prof, molec) &&
                         strncmp(molec_string+allmol, "1", 1) == 0) {
                    if (fprintf(fp, "%10.3e",
                                cpl_table_get_float(prof, molec, level,
                                                    NULL)) < 0) {
                        success = 0;
                    }
                    mol++;
                } else {
                    if (fprintf(fp, "%10s", "0.000e+00") < 0) {
                        success = 0;
                    }
                }
                if (mol > last_mol) {
                    break;
                }
                /* 8 molecules per line */
                if ((allmol + 1) % 8 == 0) {
                    if (fputc('\n', fp) == EOF) {
                        success = 0;
                    }
                }
            }
            if ((last_mol+1) % 8 != 0) {
                if (fputc('\n', fp) == EOF) {
                    success = 0;
                }
            }
        }

        if (fputs("-1\n", fp) == EOF) {
            success = 0;
        }
        if (fprintf(fp, "$ Transfer to ASCII plotting data\n") < 0 ||
            fprintf(fp, " HI=0 F4=0 CN=0 AE=0 EM=0 SC=0 FI=0 PL=1 TS=0 AM=0 "
                    "MG=0 LA=0 MS=0 XS=0    0    0\n") < 0 ||
            fprintf(fp, "# Plot title not used\n") < 0) {
            success = 0;
        }

        p = cpl_parameterlist_find(lblrtm_setup, "lblrtm_config.DELV");
        dbl1 = cpl_parameter_get_double(p);

        p = cpl_parameterlist_find(drvpar->parlist, "trans");
        if (cpl_parameter_get_int(p) == 0) {
            /* TAPE 27 corresponds to trans=0 */
            if (fprintf(fp, "%10.4e%10.4e%10.4e%10.4e%5i%5i%5i%5i%10.3e%2i"
                        "%3i%5i\n", V1, V2, 10.2, dbl1, 1, 0, 12, 0, 1., 0,
                        0, 0) < 0) {
                success = 0;
            }
            if (fprintf(fp, "%10.4g%10.4g%10.3e%10.3e%5i%5i%5i%5i%5i%5i%2i"
                        "   %2i%3i\n", 0., 1.2, 7.02, 0.2, 4, 0, 1, 1, 0, 0,
                        1, 3, 27) < 0) {
                success = 0;
            }
        } else {
            /* TAPE 28 corresponds to trans!=0 */
            if (fprintf(fp, "%10.4e%10.4e%10.4e%10.4e%5i%5i%5i%5i%10.3e%2i"
                        "%3i%5i\n", V1, V2, 10.2, dbl1, 1, 0, 12, 0, 1., 0,
                        0, 0) < 0) {
                success = 0;
            }
            if (fprintf(fp, "%10.4g%10.4g%10.3e%10.3e%5i%5i%5i%5i%5i%5i%2i"
                        "   %2i%3i\n", 0., 1.2, 7.02, 0.2, 4, 0, 1, 0, 0, 0,
                        1, 3, 28) < 0) {
                success = 0;
            }
        }
        if (fputs("-1\n", fp) == EOF) {
            success = 0;
        }
        if (fputs("% created by lblrtm_start1\n", fp) == EOF) {
            success = 0;
        }

        cpl_array_delete(atm_molecs);
    }
    fclose(fp);

    cpl_array_delete(allmolecs);

    if (success == 0) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                     "Problem writing TAPE5 file");;
    }

    /* run LBLRTM with parameters: */
//    cpl_msg_info(cpl_func, "SETUP VARIABLES LBLRTM:");
//    cpl_msg_info(cpl_func, "tape3: %s", tape3);
//    cpl_msg_info(cpl_func, "working directory: %s", wdir);
//    cpl_msg_info(cpl_func, "molecs: %s", molec_string);
//    cpl_msg_info(cpl_func, "V1: %f", V1);
//    cpl_msg_info(cpl_func, "V2: %f\n\n", V2);

    /* show info message */
    cpl_msg_info(cpl_func, "Run LBLRTM for %.4g - %.4g µm", 1e4 / V2,
                 1e4 / V1);

    sprintf(sys, "%slblrtm", bindir);
    if (mf_basic_access(sys, F_OK) != CPL_ERROR_NONE) {
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                                     "could not find LBLRTM executable");
    }

    // DEBUG
    strcat(sys, " > /dev/null 2> /dev/null");

    /* run lblrtm */
    tw = cpl_test_get_walltime();
    d = system(sys);
    t_code += cpl_test_get_walltime() - tw;

    /* Check for signal SIGINT and stop programme if present */
    signal(SIGINT, SIG_DFL);

    /* Check LBLRTM return value */
    if (d != 0) {
        sprintf(errtxt, "%s: LBLRTM failed", MF_ERROR_EIS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_EIS, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


cpl_array *mf_lblrtm_allmolecs(void)
{
    /*!
     * \brief
     *   Return an array with all molecules accepted by LBLRTM.
     *
     * This function returns a \c CPL_ARRAY with all molecules included in
     * LBLRTM. At the time of this writing these are:
     *
     * <table class="ec" align="center">
     *     <tr>
     *       <td class="ecl">0</td>
     *       <td class="ecr">H2O</td>
     *       <td class="ecl">1</td>
     *       <td class="ecr">CO2</td>
     *       <td class="ecl">2</td>
     *       <td class="ecr">O3</td>
     *       <td class="ecl">3</td>
     *       <td class="ecr">N2O</td>
     *       <td class="ecl">4</td>
     *       <td class="ecr">CO</td>
     *     </tr>
     *     <tr>
     *       <td class="ecl">5</td>
     *       <td class="ecr">CH4</td>
     *       <td class="ecl">6</td>
     *       <td class="ecr">O2</td>
     *       <td class="ecl">7</td>
     *       <td class="ecr">NO</td>
     *       <td class="ecl">8</td>
     *       <td class="ecr">SO2</td>
     *       <td class="ecl">9</td>
     *       <td class="ecr">NO2</td>
     *     </tr>
     *     <tr>
     *       <td class="ecl">10</td>
     *       <td class="ecr">NH3</td>
     *       <td class="ecl">11</td>
     *       <td class="ecr">HNO3</td>
     *       <td class="ecl">12</td>
     *       <td class="ecr">OH</td>
     *       <td class="ecl">13</td>
     *       <td class="ecr">HF</td>
     *       <td class="ecl">14</td>
     *       <td class="ecr">HCL</td>
     *     </tr>
     *     <tr>
     *       <td class="ecl">15</td>
     *       <td class="ecr">HBR</td>
     *       <td class="ecl">16</td>
     *       <td class="ecr">HI</td>
     *       <td class="ecl">17</td>
     *       <td class="ecr">CLO</td>
     *       <td class="ecl">18</td>
     *       <td class="ecr">OCS</td>
     *       <td class="ecl">19</td>
     *       <td class="ecr">H2CO</td>
     *     </tr>
     *     <tr>
     *       <td class="ecl">20</td>
     *       <td class="ecr">HOCL</td>
     *       <td class="ecl">21</td>
     *       <td class="ecr">N2</td>
     *       <td class="ecl">22</td>
     *       <td class="ecr">HCN</td>
     *       <td class="ecl">23</td>
     *       <td class="ecr">CH3CL</td>
     *       <td class="ecl">24</td>
     *       <td class="ecr">H2O2</td>
     *     </tr>
     *     <tr>
     *       <td class="ecl">25</td>
     *       <td class="ecr">C2H2</td>
     *       <td class="ecl">26</td>
     *       <td class="ecr">C2H6</td>
     *       <td class="ecl">27</td>
     *       <td class="ecr">PH3</td>
     *       <td class="ecl">28</td>
     *       <td class="ecr">COF2</td>
     *       <td class="ecl">29</td>
     *       <td class="ecr">SF6</td>
     *     </tr>
     *     <tr>
     *       <td class="ecl">30</td>
     *       <td class="ecr">H2S</td>
     *       <td class="ecl">31</td>
     *       <td class="ecr">HCOOH</td>
     *       <td class="ecl">32</td>
     *       <td class="ecr">HO2</td>
     *       <td class="ecl">33</td>
     *       <td class="ecr">O</td>
     *       <td class="ecl">34</td>
     *       <td class="ecr">CLONO2</td>
     *     </tr>
     *     <tr>
     *       <td class="ecl">35</td>
     *       <td class="ecr">NO+</td>
     *       <td class="ecl">36</td>
     *       <td class="ecr">HOBR</td>
     *       <td class="ecl">37</td>
     *       <td class="ecr">C2H4</td>
     *       <td class="ecl">38</td>
     *       <td class="ecr">CH3OH</td>
     *     </tr>
     *   </table>
     *
     * \b OUTPUT:
     * \return  string \c CPL_ARRAY \em molecs of type \c CPL_TYPE_STRING
     */

    /*
     *  List of valid molecules:
     *   1:H2O     2:CO2     3:O3      4:N2O     5:CO      6:CH4     7:O2
     *   8:NO      9:SO2    10:NO2    11:NH3    12:HNO3   13:OH     14:HF
     *  15:HCL    16:HBR    17:HI     18:CLO    19:OCS    20:H2CO   21:HOCL
     *  22:N2     23:HCN    24:CH3CL  25:H2O2   26:C2H2   27:C2H6   28:PH3
     *  29:COF2   30:SF6    31:H2S    32:HCOOH  33:HO2    34:O      35:CLONO2
     *  36:NO+    37:HOBR   38:C2H4   39:CH3OH  40:CH3Br  41:CH3CN  42:CF4
     *  43:C4H2   44:HC3N   45:H2     46:CS     47:SO3
     */

    static cpl_array *molecs;

    molecs = cpl_array_new(47, CPL_TYPE_STRING);
    cpl_array_set_string(molecs, 0, "H2O");
    cpl_array_set_string(molecs, 1, "CO2");
    cpl_array_set_string(molecs, 2, "O3");
    cpl_array_set_string(molecs, 3, "N2O");
    cpl_array_set_string(molecs, 4, "CO");
    cpl_array_set_string(molecs, 5, "CH4");
    cpl_array_set_string(molecs, 6, "O2");
    cpl_array_set_string(molecs, 7, "NO");
    cpl_array_set_string(molecs, 8, "SO2");
    cpl_array_set_string(molecs, 9, "NO2");
    cpl_array_set_string(molecs, 10, "NH3");
    cpl_array_set_string(molecs, 11, "HNO3");
    cpl_array_set_string(molecs, 12, "OH");
    cpl_array_set_string(molecs, 13, "HF");
    cpl_array_set_string(molecs, 14, "HCL");
    cpl_array_set_string(molecs, 15, "HBR");
    cpl_array_set_string(molecs, 16, "HI");
    cpl_array_set_string(molecs, 17, "CLO");
    cpl_array_set_string(molecs, 18, "OCS");
    cpl_array_set_string(molecs, 19, "H2CO");
    cpl_array_set_string(molecs, 20, "HOCL");
    cpl_array_set_string(molecs, 21, "N2");
    cpl_array_set_string(molecs, 22, "HCN");
    cpl_array_set_string(molecs, 23, "CH3CL");
    cpl_array_set_string(molecs, 24, "H2O2");
    cpl_array_set_string(molecs, 25, "C2H2");
    cpl_array_set_string(molecs, 26, "C2H6");
    cpl_array_set_string(molecs, 27, "PH3");
    cpl_array_set_string(molecs, 28, "COF2");
    cpl_array_set_string(molecs, 29, "SF6");
    cpl_array_set_string(molecs, 30, "H2S");
    cpl_array_set_string(molecs, 31, "HCOOH");
    cpl_array_set_string(molecs, 32, "HO2");
    cpl_array_set_string(molecs, 33, "O");
    cpl_array_set_string(molecs, 34, "CLONO2");
    cpl_array_set_string(molecs, 35, "NO+");
    cpl_array_set_string(molecs, 36, "HOBR");
    cpl_array_set_string(molecs, 37, "C2H4");
    cpl_array_set_string(molecs, 38, "CH3OH");
    cpl_array_set_string(molecs, 39, "CH3Br");
    cpl_array_set_string(molecs, 40, "CH3CN");
    cpl_array_set_string(molecs, 41, "CF4");
    cpl_array_set_string(molecs, 42, "C4H2");
    cpl_array_set_string(molecs, 43, "HC3N");
    cpl_array_set_string(molecs, 44, "H2");
    cpl_array_set_string(molecs, 45, "CS");
    cpl_array_set_string(molecs, 46, "SO3");

    return molecs;
}


cpl_error_code mf_lblrtm_renameoutput(const char *wdir, const int a,
                                      const char mode)
{
    /*!
     * \brief
     *   Rename LBLRTM output files: TAPE27/TAPE28.
     *
     * This private function reads the working directory from the LBLRTM setup
     * parameterlist \em lblrtm_setup and renames the TAPE27 / TAPE28 files to
     * TAPE27_XX / TAPE28_XX, respectively, where "XX" is the number specified
     * in the parameter \em a.
     *
     * \b INPUT:
     * \param wdir  working directory
     * \param a     number of the section
     *
     * \b OUTPUT:
     * \param mode  radiance ('r') or transmission ('t')
     *
     * \return      CPL_ERROR_NONE on success,
     *              MF_ERROR_CD or CPL_ERROR_FILE_IO else
     */

    char ptr[MF_MAXLEN],
         range[255],
         str[MF_MAXLEN*3];

    int i;

    char mode_str[3] = "27"
            ;
    cpl_errorstate err_state;

    sprintf(range, "%i", a);

    /* change working directory */
    if (chdir(wdir)){
        return cpl_error_set_message(cpl_func, MF_ERROR_CD, MF_ERROR_CD_TXT);
    }

    if (mode == 't') {
        sprintf(mode_str, "28");
    }

    /* rename existing TAPE27 & TAPE28 */
    sprintf(ptr, "%s/TAPE%s", wdir, mode_str);
    sprintf(str, "mv %s/TAPE%s %s/TAPE%s_%s",
            wdir, mode_str, wdir, mode_str, range);

    err_state = cpl_errorstate_get();
    if (mf_basic_access(ptr, F_OK) == CPL_ERROR_NONE) {
        if ((i = system(str))){};
    } else {
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                     "Could not rename file: TAPE27");
    }
    cpl_errorstate_set(err_state);

    return CPL_ERROR_NONE;
}

/**@}*/
