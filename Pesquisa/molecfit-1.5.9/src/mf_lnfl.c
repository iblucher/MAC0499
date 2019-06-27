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
 * \file mf_lnfl.c
 *
 * Routines for running LNFL (line file creation)
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 * \since  20 Jun 2012
 * \date   28 Aug 2013
 */


/*****************************************************************************
 *                                INCLUDES                                   *
 ****************************************************************************/

#include "mf_molecfit.h"
#include <mf_lnfl.h>


/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code mf_lnfl(const mfdrv *drvpar)
{
    /*!
     * Calls LNFL for the different wavelength ranges saved in the ::mfdrv
     * parameter structure.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * -
     *
     * \b ERRORS:
     * - Invalid object structure
     * - see ::mf_lnfl_call
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    char errtxt[MF_MAXLEN];
    int nrange = 0, singlespec = 0, j = 0;

    /* Get number of ranges */
    p = cpl_parameterlist_find(drvpar->parlist, "nrange");
    nrange = cpl_parameter_get_int(p);

    /* Calculate one molecular spectrum for full wavelength range? */
    p = cpl_parameterlist_find(drvpar->parlist, "singlespec");
    singlespec = cpl_parameter_get_int(p);
    if (singlespec == 1) {
        nrange = 1;
    }

    /* Check size of range table */
    if (singlespec == 0 && cpl_table_get_nrow(drvpar->rangetab) != nrange) {
        sprintf(errtxt, "%s: mfdrv *drvpar (unexpected size of rangetab)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Call LNFL for the different range numbers */
    for (j = 0; j < nrange; j++) {
        status = mf_lnfl_call(drvpar, j+1);
        if (status != CPL_ERROR_NONE) {
            return status;
        }
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_lnfl_call(const mfdrv *drvpar, const int range)
{
    /*!
     * \brief
     *   Run LNFL.
     *
     * This function checks for the existence of a TAPE3 file for the given
     * wavelength range (indicated by a range number). If the required file
     * does not exist, it runs LNFL using the following input data in the tree
     * under \em dir:
     *
     * <table class="ec" align="center">
     *     <tr>
     *       <td class="ecl">/data/hitran/aer_v_3.2</td>
     *       <td class="ecr">HITRAN data base</td>
     *     </tr>
     *     <tr>
     *       <td class="ecl">/output/lnfl_TAPE3</td>
     *       <td class="ecr">TAPE3 file</td>
     *     </tr>
     *     <tr>
     *       <td class="ecl">/output/wdir</td>
     *       <td class="ecr">working directory (will be deleted upon exit)</td>
     *     </tr>
     *     <tr>
     *       <td class="ecl">/bin/lnfl</td>
     *       <td class="ecr">LNFL executable</td>
     *     </tr>
     * </table>
     *
     * All available molecules are used.
     *
     * \note
     *   - If the file ./output/lnfl_TAPE3 exists under \em dir it gets
     *     overwritten.
     *   - Computation time of LBLRTM depends on wavelength range in LNFL:
     *     300-33333 takes 2x as long as 2981-3056
     *   - Computation time of LNFL depends on number of molecules:
     *     3 molecules take only 1/3 of the time all 39 molecules require
     *
     *
     * \b INPUT:
     * \param drvpar  driver parameter structure
     * \param range   range number
     *
     * \b OUTPUT:
     * \return        CPL_ERROR_NONE on success,
     *                CPL_ERROR_FILE_NOT_CREATED,
     *                MF_ERROR_CD,
     *                MF_ERROR_FOPEN,
     *                MF_ERROR_UFS on error.
     */

    /*
     * AVAILABLE MOLECULAR SPECIES:
     *  1) H2O   2) CO2   3)   O3   4)  N2O   5)  CO   6) CH4   7)    O2
     *  8)  NO   9) SO2  10)  NO2  11)  NH3  12)HNO3  13)  OH  14)    HF
     * 15) HCL  16) HBR  17)   HI  18)  CLO  19) OCS  20)H2CO  21)  HOCL
     * 22)  N2  23) HCN  24)CH3CL  25) H2O2  26)C2H2  27)C2H6  28)   PH3
     * 29)COF2  30) SF6  31)  H2S  32)HCOOH  33) HO2  34)   O  35)CLONO2
     * 36) NO+  37)HOBR  38) C2H4  39)CH3OH  40)CH3Br 41)CH3CN 42)   CF4
     * 43) C4H2 44) HC3N 45)    H2 46)    CS 47)   SO3
     */

    FILE *fp;                          /* file pointer */
    char line[MF_MAXLEN];              /* line in setup file */
    char confdir[MF_MAXLEN];
    char setup_file[MF_MAXLEN];
    char str[MF_MAXLEN],
         *cptr;
    char line_db[MF_MAXLEN];
    int  line_db_fmt = 0;

    float wn_start, wn_end;            /* start & end wavenumber in [cm-1] */

    char ptr[MF_MAXLEN],
         rmwdir[MF_MAXLEN];

    char tape1[MF_MAXLEN],
         tape3[MF_MAXLEN],
         tape5[MF_MAXLEN],
         basedir[MF_MAXLEN],
         outdir[MF_MAXLEN],
         datadir[MF_MAXLEN],
         lnfldir[MF_MAXLEN],
         bindir[MF_MAXLEN],
         wdir[MF_MAXLEN],
         lnfl_exe[MF_MAXLEN];
    const char *molecs = NULL;

    int d = 0;

    cpl_parameter *p;
    cpl_errorstate err_state;
    cpl_error_code err_code;

    struct stat buf;

    /* get molecules */
    p = cpl_parameterlist_find(drvpar->parlist, "lbl_molecs");
    molecs = cpl_parameter_get_string(p);

    /* get wavenumber range */
    wn_start = cpl_table_get(drvpar->rangetab, "wn_start", range-1, NULL);
    wn_end = cpl_table_get(drvpar->rangetab, "wn_end", range-1, NULL);

    /* exit if range is empty */
    if (wn_start == 0. && wn_end == 0.) {
        return CPL_ERROR_NONE;
    }

    /* get base directory */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));

    /* get data directory */
    p = cpl_parameterlist_find(drvpar->parlist, "datadir");
    // mf_basic_abspath(datadir, cpl_parameter_get_string(p), ".");
    strcpy(datadir, cpl_parameter_get_string(p));

    /* get lnfl directory */
    sprintf(lnfldir, "%slnfl/", datadir);

    /* create lnfl directory */
    err_state = cpl_errorstate_get();
    if (mf_basic_access(lnfldir, F_OK) != CPL_ERROR_NONE) {
        if (mkdir(lnfldir, 0755)) {
            return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                         "Could not create lnfl "
                                         "directory: %s", lnfldir);
        }
    }
    cpl_errorstate_set(err_state);

    /* get bin directory */
    p = cpl_parameterlist_find(drvpar->parlist, "bindir");
    mf_basic_abspath(bindir, cpl_parameter_get_string(p), basedir);

    /* get output directory */
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(outdir, cpl_parameter_get_string(p), curdir);

    /* get line list / line list format ------------------------------------*/
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

            if (strcmp(str, "_line_db") == 0){
                cptr = strtok(line, "=");
                cptr = strtok(NULL, "=");
                mf_basic_strtrim_inplace(cptr);    /* trim line */
                sprintf(line_db, "%s", cptr);

            }

            if (strcmp(str, "_line_db_fmt") == 0){
                cptr = strtok(line, "=");
                cptr = strtok(NULL, "=");
                mf_basic_strtrim_inplace(cptr);    /* trim line */
                if ( mf_basic_isnumber(cptr) != CPL_TRUE ){
                    return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                "Wrong line database format: '%s'\nMust be "
                                "a number. Please check parameter "
                                "_line_db_fmt in file %s\n", cptr,
                                setup_file);
                }
                if (strpbrk(cptr, ".,") != NULL){
                    return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                "Wrong line database format: '%s'\nMust be "
                                "a positive integer number. Please check "
                                "parameter _line_db_fmt in file %s\n", cptr,
                                setup_file);
                }
                if ( atoi(cptr) <= 0 ){
                    return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                "Wrong line database format: '%s'\nMust be "
                                "a positive integer number. Please check "
                                "parameter _line_db_fmt in file %s\n", cptr,
                                setup_file);
                }
                line_db_fmt = atoi(cptr);
            }

        }
    }
    fclose(fp);

    /* initialise setup structure */
    sprintf(tape1, "%shitran/%s", datadir, line_db);
    sprintf(wdir, "%swdir", outdir);
    sprintf(lnfl_exe, "%slnfl", bindir);

    if ((err_code = mf_lnfl_files(tape3, tape5, lnfldir,
                                  wn_start, wn_end, molecs)) !=
                                         CPL_ERROR_NONE) {
        if (err_code == CPL_ERROR_FILE_NOT_CREATED) {
            // no free slot found
            return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                         "Could not find slot for LNFL "
                                         "output files. Consider deleting "
                                         "files from %s", lnfldir);
        } else {
            return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                         MF_ERROR_SUBROUTINE_TXT);
        }
    }

    /* store LNFL TAPE3 filename for later use */
    cpl_table_set_string(drvpar->rangetab, "lnfl", range-1, tape3);

    /* run LNFL with parameters: */
    //cpl_msg_info(cpl_func, "SETUP VARIABLES LNFL:");
    //cpl_msg_info(cpl_func, "tape1: %s", tape1);
    //cpl_msg_info(cpl_func, "tape3: %s", tape3);
    //cpl_msg_info(cpl_func, "working directory: %s", wdir);
    //cpl_msg_info(cpl_func, "lnfl_exe: %s", lnfl_exe);
    //cpl_msg_info(cpl_func, "molecs: %s", molecs);
    //cpl_msg_info(cpl_func, "wn_start: %f", wn_start);
    //cpl_msg_info(cpl_func, "wn_end: %f\n\n", wn_end);

    /* do not overwrite existing lnfl_TAPE3 of non-zero length */
    err_state = cpl_errorstate_get();
    if (mf_basic_access(tape3, F_OK) == CPL_ERROR_NONE) {
        d = stat(tape3, &buf);
        if (d == 0 && buf.st_size != 0) {
           cpl_msg_info(cpl_func, "LNFL line file for %.4g - %.4g µm already "
                        "exists", 1e4 / wn_end, 1e4 / wn_start);
           return CPL_ERROR_NONE;
        }
    }
    cpl_errorstate_set(err_state);

    /* create working directory */
    err_state = cpl_errorstate_get();
    sprintf(wdir, "%swdir", outdir);
    if (mf_basic_access(wdir, F_OK) != CPL_ERROR_NONE) {
        if (mkdir(wdir, 0755)) {
            return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                         "Could not create working "
                                         "directory.");
        }
    }
    cpl_errorstate_set(err_state);

    /* command for removal of working directory */
    sprintf(rmwdir, "rm -rf %s", wdir);

    /* remove existing TAPE1 & TAPE3 */
    err_state = cpl_errorstate_get();
    sprintf(ptr, "%s/TAPE1", wdir);
    if (mf_basic_access(ptr, F_OK) == CPL_ERROR_NONE) {
        remove(ptr);
    }
    sprintf(ptr, "%s/TAPE3", wdir);
    if (mf_basic_access(ptr, F_OK) == CPL_ERROR_NONE) {
        remove(ptr);
    }
    cpl_errorstate_set(err_state);

    /* create symbolic link to HITRAN in working directory */
    sprintf(ptr, "%s/TAPE1", wdir);
    if (mf_basic_access(ptr, F_OK) != CPL_ERROR_NONE) {
        sprintf(ptr, "ln -s %s %s/TAPE1", tape1, wdir);
        if (system(ptr)) {
            d = system(rmwdir);
            return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                         "could not create symbolic link");
        }
        cpl_errorstate_set(err_state);
    }

    /* change to working directory */
    if (chdir(wdir)){
        d = system(rmwdir);
        return cpl_error_set_message(cpl_func, MF_ERROR_CD, MF_ERROR_CD_TXT);
    }

    /* create TAPE5 file */
    /*  1st line: 72 chars user info (e.g.: "$ f100 format")
     *  2nd line: lower & upper wavenumber (25cm-1 wider than required)
     *            format: F10.3,  F10.3 (e.g.: "   300.      3500.")
     *  3rd line: molecule indicator & hollerith indicator
     *            format: 39I1,3X,     A40
     *   (e.g.: "111111111111111111111111111111111111111   NBLK1 LNOUT")
     *  4th line: fixed: "%%%%%%%%%%%%%%%%%%"
     *  5th line: fixed: "1234567890123456789012345678901234567890"
     *             cont. "1234567890123456789012345678901234567890"
     */
    if ((fp = fopen("TAPE5", "w")) == NULL) {
        d = chdir(basedir);
        d = system(rmwdir);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                     MF_ERROR_FOPEN_TXT);
    } else {
        /* 1st line */
        if (fprintf(fp, "$ created by lnfl: %s\n", tape3) == EOF) {
            fclose(fp);
            d = chdir(basedir);
            d = system(rmwdir);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS,
                                         MF_ERROR_UFS_TXT);
        }
        /* 2nd line */
        if (fprintf(fp, "%10.3f%10.3f\n", wn_start, wn_end) == EOF) {
            fclose(fp);
            d = chdir(basedir);
            d = system(rmwdir);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS,
                                         MF_ERROR_UFS_TXT);
        }
        /* 3rd line */
        if (fprintf(fp, "%47s    LNOUT F%i\n", molecs, line_db_fmt) == EOF) {

//         if (fprintf(fp, "%39s   LNOUT F160\n", molecs) == EOF) {
            fclose(fp);
            d = chdir(basedir);
            d = system(rmwdir);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS,
                                         MF_ERROR_UFS_TXT);
        }
        /* 4th line */
        if (fputs("%%%%%%%%%%%%%%%%%%\n", fp) == EOF) {
            fclose(fp);
            d = chdir(basedir);
            d = system(rmwdir);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS,
                                         MF_ERROR_UFS_TXT);
        }
        /* 5th line */
        if (fputs("1234567890123456789012345678901234567890"
                  "1234567890123456789012345678901234567890\n", fp) == EOF) {
            fclose(fp);
            d = chdir(basedir);
            d = system(rmwdir);
            return cpl_error_set_message(cpl_func, MF_ERROR_UFS,
                                         MF_ERROR_UFS_TXT);
        }
    }

    fclose(fp);

    /* show info message */
    cpl_msg_info(cpl_func, "Run LNFL for %.4g - %.4g µm", 1e4 / wn_end,
                 1e4 / wn_start);

    /* run LNFL */
    if (mf_basic_access(lnfl_exe, F_OK) != CPL_ERROR_NONE) {
        d = chdir(basedir);
        d = system(rmwdir);
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                                     "File not found.");

    }

    // DEBUG
    strcat(lnfl_exe, " > /dev/null 2> /dev/null");
    
    if ((d = system(lnfl_exe))){};

    /* copy TAPE3 file */
    sprintf(ptr, "cp TAPE3 %s", tape3);
    d = system(ptr);

    /* copy TAPE5 file */
    sprintf(ptr, "cp TAPE5 %s", tape5);
    d = system(ptr);

    /* change to base directory */
    if (chdir(basedir)){
        d = system(rmwdir);
        return cpl_error_set_message(cpl_func, MF_ERROR_CD, MF_ERROR_CD_TXT);
    }

    /* remove working directory */
    d = system(rmwdir);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_lnfl_files(char *tape3_name, char *tape5_name,
                             const char *dir, const double wn_start,
                             const double wn_end, const char *molecs)
{
    /*!
     * \brief
     *   Find TAPE3/5 files for LNFL.
     *
     * This function reads all TAPE5 files in data/lnfl. All TAPE5 files are
     * checked for consistency with the setup of the current run. If a
     * matching TAPE5 file is found, the corresponding TAPE3/5 filenames are
     * returned. If no matching TAPE5 file was found, the first free "slot"
     * for a new TAPE3/5 file combination is being searched for. All TAPE3/5
     * files have the same naming scheme:
     *
     * \c lnfl_TAPE3_# and \c lnfl_TAPE5_#
     *
     * with "#" being an arbitrary \c long number.
     *
     * \note
     *   Currently, 10.000 slots are allowed. If all slots are filled, an
     *   error message is produced, with the suggestion, to clean-up the
     *   data/lnfl directory.
     *
     * \note
     *   The molecule parameter \em molecs is a string containing all
     *   molecules represented by "1" or "0", depending on whether it is to be
     *   included in the fit or not, respectively. Currently, LNFL allows for
     *   39 molecules. Thus, an example using all molecules would be:
     *   "111111111111111111111111111111111111111"
     *
     * \b INPUT:
     * \param dir        directory with LNFL files (default: data/lnfl)
     * \param wn_start   starting wavenumber [cm-1]
     * \param wn_end     ending wavenumber [cm-1]
     * \param molecs     molecule list ("binary" string)
     *
     * \b OUTPUT:
     * \param tape3_name TAPE3 file name
     * \param tape5_name TAPE5 file name
     * \return           CPL_ERROR_NONE on success,
     *                   MF_ERROR_FOPEN,
     *                   CPL_ERROR_FILE_NOT_FOUND,
     *                   CPL_ERROR_FILE_NOT_CREATED on error
     */

    cpl_errorstate err_state;

    char lnfldir[MF_MAXLEN] = "", file_name[MF_MAXLEN] = "";
    char molecs_file[MF_MAXLEN] = "";
    char wn_limits[MF_MAXLEN] = "0 0";
    char *wn_start_str = NULL, *wn_end_str = NULL;
    double wn_start_file = 0., wn_end_file = 0.;
    long i = 0, n_slots = 10000;

    FILE *fp;
    DIR *dp;
    struct dirent *dir_ent;

    /* Make sure that a slash is at the end of the directory string */
    sprintf(lnfldir, "%s", dir);
    mf_basic_dirslash(lnfldir);

    /* Find match in directory */
    dp = opendir(dir);

    if (dp != NULL) {
        /* Loop through files */
        while ((dir_ent = readdir (dp)) != NULL) {
            /* Select TAPE5 files */
            if (strncmp(dir_ent->d_name, "lnfl_TAPE5_", 11) == 0) {
                /* Set TAPE5 filename including path */
                sprintf(tape5_name, "%s%s", lnfldir, dir_ent->d_name);

                /* Open file & check content: wn_start, wn_end, molecs */
                if ((fp = fopen(tape5_name, "r")) == NULL) {
                    (void) closedir(dp);
                    return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                                 MF_ERROR_FOPEN_TXT);
                } else {
                    if (fgets(molecs_file, MF_MAXLEN-1, fp)){};
                    if (fgets(wn_limits, MF_MAXLEN-1, fp)){};
                    wn_start_str = strtok(wn_limits, " ");
                    wn_end_str = strtok(NULL, " ");
                    wn_start_file = strtod(wn_start_str, NULL);
                    wn_end_file = strtod(wn_end_str, NULL);
                    if (fscanf(fp, "%s", molecs_file)){};
                    fclose(fp);

                    if ((long) round(wn_start_file*1000.) ==
                            (long) round(wn_start*1000.) &&
                            (long) round(wn_end_file*1000.) ==
                                    (long) round(wn_end*1000.) &&
                            strcmp(molecs_file, molecs) == 0) {
                        /* Match found -> set TAPE3 filename including path */
                        sprintf(tape3_name, "%slnfl_TAPE3_%s", lnfldir,
                                (dir_ent->d_name+11));

                        (void) closedir(dp);

                        return CPL_ERROR_NONE;
                    }
                }
                /* Not compatible -> next */
            }
            /* Not a TAPE5 file -> next */
        }
        /* No matching TAPE5 files found */
        (void) closedir(dp);
    }
    else {
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_FOUND,
                                     "LNFL directory not found: %s", dir);
    }

    /* Find first free filename slot */
    for (i = 0; i < n_slots; i++) {
        sprintf(file_name, "%slnfl_TAPE3_%li", lnfldir, i);
        sprintf(tape5_name, "%slnfl_TAPE5_%li", lnfldir, i);

        err_state = cpl_errorstate_get();
        if (mf_basic_access(file_name, F_OK) == CPL_ERROR_NONE) {
            /* Slot used -> next */
            cpl_errorstate_set(err_state);
        } else {
            /* First free slot found -> set TAPE3/5 filename including path */
            sprintf(tape3_name, "%s", file_name);
            sprintf(tape5_name, "%slnfl_TAPE5_%li", lnfldir, i);

            cpl_errorstate_set(err_state);

            return CPL_ERROR_NONE;
        }
    }

    mf_basic_initstring(tape3_name, MF_MAXLEN);
    mf_basic_initstring(tape5_name, MF_MAXLEN);

    return CPL_ERROR_FILE_NOT_CREATED;
}

/**@}*/
