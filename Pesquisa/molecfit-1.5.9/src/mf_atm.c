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
 * \file mf_atm.c
 *
 * Routines for derivation and handling of atmospheric profiles
 *
 * \author Marco Barden, Stefan Noll, & ESO In-Kind Team Innsbruck
 * \since  20 Jun 2012
 * \date   27 Jul 2014
 */


/*****************************************************************************
 *                                 INCLUDES                                  *
 ****************************************************************************/

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif
#include "mf_molecfit.h"
#include <mf_atm.h>
#include <stdio.h>
static const char * base_gdas_url =
    "ftp://ftp.eso.org/pub/dfs/pipelines/skytools/molecfit/gdas/";

/*****************************************************************************
 *                                  CODE                                     *
 ****************************************************************************/

cpl_error_code mf_atm_createatm(cpl_table *out_profile, const mfdrv *drvpar)
{
    /*!
     * Creates a combined atmospheric profile based on GDAS data, which will
     * be downloaded from a webserver if necessary, and a standard atmospheric
     * profile (from MIPAS) based on the parameters specified in the ::mfdrv
     * parameter structure. If the parameter \e gdas_prof is set to the
     * default "auto", this profile building mode is performed. The selection
     * "none" will only use the standard atmospheric profile. In all other
     * cases, the file path and name of a specific GDAS-like profile is
     * expected. This will then be merged with the standard profile.
     *
     * \b INPUT:
     * \param drvpar       ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param out_profile  CPL table with GDAS & MIPAS combined profile
     *
     * \b ERRORS:
     * - Invalid object value(s)
     * - Error in subroutine
     * - see subroutines
     */

    cpl_error_code err_code = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_array *all_molecs, *atm_molecs, *timestamps, *molecs = NULL;
    cpl_table *atm_profile, *gdas_profile1, *gdas_profile2, *gdas_profile;
    char **mol_arr = NULL;
    char errtxt[MF_MAXLEN], basedir[MF_MAXLEN], filename[MF_MAXLEN] = "\0";
    char gdas_file[MF_MAXLEN] = "\0";
    cpl_boolean exmol = CPL_FALSE;
    int nmolec = 0, nrow_all = 0, j = 0, i = 0, nrow = 0;

    /* Get molecules from parameter list */
    p = cpl_parameterlist_find(drvpar->parlist, "nmolec");
    nmolec = cpl_parameter_get_int(p);
    mol_arr = cpl_table_get_data_string(drvpar->molectab, "list_molec");

    /* Get LBLRTM molecules */
    all_molecs = mf_lblrtm_allmolecs();
    nrow_all = cpl_array_get_size(all_molecs);

    /* Check whether LBLRTM can handle the selected molecules */
    for (j = 0; j < nmolec; j++) {
        for (exmol = CPL_FALSE, i = 0; i < nrow_all; i++) {
            if (strcmp(cpl_array_get_string(all_molecs, i), mol_arr[j])
                == 0) {
                exmol = CPL_TRUE;
                break;
            }
        }
        if (exmol == CPL_FALSE) {
            cpl_array_delete(all_molecs);
            sprintf(errtxt, "%s: cpl_table drvpar->molectab (molecule %s "
                    "cannot be handled)", MF_ERROR_IOV_TXT, mol_arr[j]);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s",
                                         errtxt);
        }
    }

    /* Delete temporary array */
    cpl_array_delete(all_molecs);

    /* Get base directory */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));

    /* Get file name of standard profile */
    p = cpl_parameterlist_find(drvpar->parlist, "ref_atm");
    sprintf(filename, "%s/profiles/mipas/", mf_get_datadir());
    strcat(filename, cpl_parameter_get_string(p));

    /* Write info message */
    cpl_msg_info(cpl_func, "Read standard profile %s", filename);

    /* Read standard profile */
    atm_profile = cpl_table_new(1);
    if (mf_atm_readatm(atm_profile, filename) != CPL_ERROR_NONE) {
        cpl_table_delete(atm_profile);
        return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                     "Could not read reference profile: %s",
                                     filename);
    }

    /* Get column labels from standard profile */
    atm_molecs = cpl_table_get_column_names(atm_profile);
    nrow = cpl_array_get_size(atm_molecs);

    /* Check existence of selected molecules in standard profile */
    for (j = 0; j < nmolec; j++) {
        for (exmol = CPL_FALSE, i = 0; i < nrow; i++) {
            if (strcmp(cpl_array_get_string(atm_molecs, i), mol_arr[j])
                == 0) {
                exmol = CPL_TRUE;
                break;
            }
        }
        if (exmol == CPL_FALSE) {
            cpl_table_delete(atm_profile);
            cpl_array_delete(atm_molecs);
            sprintf(errtxt, "%s: cpl_table drvpar->molectab (molecule %s not "
                    "found in %s)", MF_ERROR_IOV_TXT, mol_arr[j], filename);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s",
                                         errtxt);
        }
    }

    /* Delete temporary array */
    cpl_array_delete(atm_molecs);

    /* Copy selected molecules to CPL array */
    molecs = cpl_array_new(0, CPL_TYPE_STRING);
    mf_basic_col2arr(molecs, drvpar->molectab, "list_molec");

    /* Get name of specific GDAS profile if available */
    p = cpl_parameterlist_find(drvpar->parlist, "gdas_prof");
    strncpy(gdas_file, cpl_parameter_get_string(p), MF_MAXLEN);

    /* Get GDAS profile(s) depending on parameter gdas_prof and merge it with
       standard atmosphere */

    if (strncmp(gdas_file, "none", 4) == 0) {

        /* Take standard profile as output profile (no GDAS profile) */
        cpl_msg_info(cpl_func, "Do not consider GDAS profiles");
        cpl_table_delete(out_profile);
        out_profile = cpl_table_duplicate(atm_profile);

         /* Delete temporary structures */
        cpl_array_delete(molecs);
        cpl_table_delete(atm_profile);

   } else if (strncmp(gdas_file, "auto", 4) == 0) {

        /* Automatic retrieval of GDAS profiles */

        /* Initialise GDAS profile tables and tag array */
        gdas_profile1 = cpl_table_new(1);
        gdas_profile2 = cpl_table_new(1);
        timestamps = cpl_array_new(3, CPL_TYPE_FLOAT);

        /* Get GDAS profiles from local TAR archive or via grib */
        if ((err_code = mf_atm_getgdas_auto(&gdas_profile1, &gdas_profile2,
                                            timestamps, drvpar)) !=
            CPL_ERROR_NONE) {
            cpl_array_delete(molecs);
            cpl_table_delete(atm_profile);
            cpl_table_delete(gdas_profile1);
            cpl_table_delete(gdas_profile2);
            cpl_array_delete(timestamps);
            return err_code;
        }

        /* Initialise table for interpolated GDAS profiles */
        gdas_profile = cpl_table_new(1);

        /* Interpolate profiles linearly between two points in time */
        err_code = mf_atm_interpolprofile(gdas_profile, gdas_profile1,
                                          gdas_profile2, timestamps);

        /* Delete temporary structures */
        cpl_array_delete(timestamps);
        cpl_table_delete(gdas_profile1);
        cpl_table_delete(gdas_profile2);

        /* Check for errors */
        if (err_code != CPL_ERROR_NONE) {
            cpl_array_delete(molecs);
            cpl_table_delete(atm_profile);
            cpl_table_delete(gdas_profile);
            return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                         "Error in linear interpolation");
        }

        /* Merge GDAS and standard profiles */
        err_code = mf_atm_convertgdas(out_profile, atm_profile,
                                      gdas_profile, molecs, drvpar);

        /* Delete temporary structures */
        cpl_array_delete(molecs);
        cpl_table_delete(atm_profile);
        cpl_table_delete(gdas_profile);

        /* Check for errors */
        if (err_code != CPL_ERROR_NONE) {
            return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                         "Could not merge profiles");
        }

    } else {

        /* Take user-provided GDAS-like profile */

        /* Initialise GDAS profile table */
        gdas_profile = cpl_table_new(1);

        /* Get GDAS-like profile */
        if ((err_code = mf_atm_getgdas_user(gdas_profile, drvpar)) !=
            CPL_ERROR_NONE) {
            cpl_array_delete(molecs);
            cpl_table_delete(atm_profile);
            cpl_table_delete(gdas_profile);
            return err_code;
        }

        /* Merge GDAS and standard profile */
        err_code = mf_atm_convertgdas(out_profile, atm_profile, gdas_profile,
                                      molecs, drvpar);

        /* Delete temporary structures */
        cpl_array_delete(molecs);
        cpl_table_delete(atm_profile);
        cpl_table_delete(gdas_profile);

        /* Check for errors */
        if (err_code != CPL_ERROR_NONE) {
            return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                         "Could not merge profiles");
        }

    }

    return CPL_ERROR_NONE;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Hour closest in the file.
 *
 * @param filestr            File name including path.
 * @param found_hour         output: Hour found.
 * @param gdasdir            Path of the GDAS DB.
 * @param lon                Longitude.
 * @param lat                Latitude.
 * @param year               Year.
 * @param month              Month.
 * @param day                Day.
 * @param utch               Is Hour in UTC or not.
 * @param direction          Direction of the position.
 *
 * @return error code        1 on success, 0 on failure.
 *
 */
/* ---------------------------------------------------------------------------*/
static int mf_gdas_get_closest_file(
    char                     *filestr,
    double                   *found_hour,
    const char               *gdasdir,
    const double             lon,
    const double             lat,
    const int                year,
    const int                month,
    const int                day,
    int                      utch,
    int                      direction)
{
    cpl_errorstate err_state = cpl_errorstate_get();

    direction = direction > 0 ? 1 : -1;

    for (int i = 0; i < 6; i++) {

        char filebasestr[MF_MAXLEN];

        char date[100];
        char hstr[10];

        int lyear  = year;
        int lmonth = month;
        int lday   = day;

        if (utch < 0 || utch > 23) {

            long tmpday;
            mf_basic_greg2jd(&tmpday, year,    month, day   );

            tmpday += direction;
            mf_basic_jd2greg(&lyear,  &lmonth, &lday, tmpday);

            utch += -24 * direction;
        }

        sprintf(hstr, "%02i", utch);
        sprintf(date, "%4i%02i%02i", lyear, lmonth, lday);

        mf_atm_concatgdasname(filestr,     gdasdir, lon, lat, date, hstr);
        mf_atm_concatgdasname(filebasestr, "",      lon, lat, date, hstr);

        if (mf_basic_access(filestr, F_OK) != CPL_ERROR_NONE) {

            char gdas_filename[  MF_MAXLEN];
            char sys_tar_local[  MF_MAXLEN];
            char sys_tar_general[MF_MAXLEN];
            char sys_mkdir[      MF_MAXLEN];
            char sys_mv[         MF_MAXLEN];

            sprintf(gdas_filename,   "gdas_profiles_C%+.1f%+.1f.tar.gz",                   lon, lat);
            sprintf(sys_tar_local,   "cd %s 2>/dev/null && tar -zxf %s/%s %s 2>/dev/null", gdasdir, gdasdir, gdas_filename, filebasestr);
            sprintf(sys_tar_general, "tar -zxf %s/profiles/gdas/%s %s 2>/dev/null",        mf_get_datadir(), gdas_filename, filebasestr);
            sprintf(sys_mkdir,       "mkdir -p %s",                                        gdasdir);
            sprintf(sys_mv,          "mv %s %s",                                           filebasestr, filestr);

            /* Check temporary path */
            if (system(sys_tar_local) != 0) {
              /* Check general path */
              if (system(sys_tar_general) == 0) {
                  cpl_msg_info(cpl_func, "(mf_atm      ) %s", sys_tar_general);
                  cpl_msg_info(cpl_func, "(mf_atm      ) %s", sys_mv         );
                  system(sys_mkdir);
                  system(sys_mv);
              }
            }
        }

        if (mf_basic_access(filestr, F_OK) == CPL_ERROR_NONE) {
            *found_hour = utch;
            cpl_errorstate_set(err_state);
            return 1;
        }

        utch += direction;
    }

    cpl_errorstate_set(err_state);

    return 0;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Get number of entries in gdas tarball.
 *
 * @param path               File path of tarball.
 *
 * @return long              Number of files in tarball, 0 on error.
 *
 */
/* ---------------------------------------------------------------------------*/
static long mf_gdas_get_tarball_nentries(
    const char               *path)
{
    FILE * stream;
    char * endptr;
    char sys[MF_MAXLEN];
    long nentries = 0;
    sprintf(sys, "tar -tf \"%s\"  2>/dev/null| wc -l", path);
    stream = popen(sys, "r");
    if (stream == NULL) {
        return 0;
    }
    if (fread(sys, 1, MF_MAXLEN, stream) > 0) {
        sys[MF_MAXLEN - 1] = 0;
        nentries = strtol(sys, &endptr, 10);
        if (endptr == sys) {
            nentries = 0;
        }
    }
    pclose(stream);
    return nentries;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Attempt to update gdas tarball from eso ftp server.
 *
 * @param lon                Longitude.
 * @param lat                Latitude.
 *
 * @return error code        CPL_ERROR_NONE or a specific error code in other case.
 *
 */
/* ---------------------------------------------------------------------------*/
static cpl_error_code update_local_db(
    const char               *gdasdir,
    double                   lon,
    double                   lat)
{
  cpl_error_ensure(gdasdir, CPL_ERROR_NULL_INPUT, return CPL_ERROR_NULL_INPUT, "GDAS directory NULL");

  char gdas_filename[  MF_MAXLEN];
  char tarball_path[   MF_MAXLEN];
  char tarball_path_dw[MF_MAXLEN];
  char ftp_path[       MF_MAXLEN];
  char sys_backup[     MF_MAXLEN];
  char sys_restore[    MF_MAXLEN];
  char sys_mkdir[      MF_MAXLEN];
  char sys_download[   MF_MAXLEN];
  char sys_copy[       MF_MAXLEN];

  sprintf(gdas_filename,   "gdas_profiles_C%+.1f%+.1f.tar.gz",   lon,                 lat          );
  sprintf(tarball_path,    "%s/profiles/gdas/%s",                mf_get_datadir(),    gdas_filename);
  sprintf(tarball_path_dw, "%s/%s",                              gdasdir,             gdas_filename);
  sprintf(ftp_path,        "%s/%s",                              base_gdas_url,       gdas_filename);

  sprintf(sys_backup,      "mv %s     %s.old 2> /dev/null",      tarball_path,        tarball_path );
  sprintf(sys_restore,     "mv %s.old %s     2> /dev/null",      tarball_path,        tarball_path );
  sprintf(sys_mkdir,       "mkdir -p %s",                        gdasdir                           );
  sprintf(sys_download,    "cd %s && curl -m 600 -O %s",         gdasdir,             ftp_path     );
  sprintf(sys_copy,        "cp %s %s",                           tarball_path_dw,     tarball_path );

  long nentries_old = -1;
  long nentries_new;

  cpl_errorstate err_state = cpl_errorstate_get();

  /* Check if exist the *.tar.gz */
  if (mf_basic_access(tarball_path, F_OK) == CPL_ERROR_NONE) {

      /* Get nentries */
      nentries_old = mf_gdas_get_tarball_nentries(tarball_path);
      cpl_msg_info(cpl_func, "(mf_atm      ) Backup old tarball (contains %ld entries)", nentries_old);

      /* Backup existing tarball, overwrites old backups */
      cpl_msg_info(cpl_func, "(mf_atm      ) %s", sys_backup);
      if (system(sys_backup) != 0) {
          cpl_msg_warning(cpl_func, "(mf_atm      ) Failed GDAS backup.");
          return CPL_ERROR_FILE_NOT_FOUND;
      }

  } else {

      /* Clean the error */
      cpl_errorstate_set(err_state);
  }

  cpl_msg_info(cpl_func, "(mf_atm      ) Attempting to update local GDAS database");

  cpl_msg_info(cpl_func, "(mf_atm      ) %s", sys_mkdir);
  system(sys_mkdir);

  cpl_msg_info(cpl_func, "(mf_atm      ) %s", sys_download);
  if (system(sys_download) != 0) {

      cpl_msg_warning(cpl_func, "(mf_atm      ) Download and update the new GDAS data failed.");
      cpl_msg_warning(cpl_func, "(mf_atm      )   GDAS data could not be transferred in their entirety within 10 min.");
      cpl_msg_warning(cpl_func, "(mf_atm      )   Consider to restart the fit to complete the transfer.");

      /* Restore database */
      if (nentries_old >= 0) {
          cpl_msg_info(cpl_func, "(mf_atm      ) %s", sys_restore);
          if (system(sys_restore) != 0) {
              cpl_msg_error(cpl_func,
                            "(mf_atm      ) Restore of GDAS backup failed, please restore %s from %s.old manually",
                            tarball_path, tarball_path);
          }
      }

      return CPL_ERROR_FILE_NOT_CREATED;

  } else {

      cpl_msg_info(cpl_func, "(mf_atm      ) %s", sys_copy);
      if (system(sys_copy) != 0) {

          cpl_msg_warning(cpl_func, "(mf_atm      ) Copy newer GDAS database to the general data failed.");

          /* Restore database */
          if (nentries_old >= 0) {
              cpl_msg_info(cpl_func, "(mf_atm      ) %s", sys_restore);
              if (system(sys_restore) != 0) {
                  cpl_msg_error(cpl_func,
                                "(mf_atm      ) Restore of GDAS backup failed, please restore %s from %s.old manually",
                                tarball_path, tarball_path);
              }
          }
      }

      /* Not throw and error */
  }

  /* Report on changes on success */
  nentries_new = mf_gdas_get_tarball_nentries(tarball_path_dw);
  if (nentries_new > nentries_old) {
      cpl_msg_info(cpl_func, "(mf_atm      ) New gdas data contains %ld new entries", nentries_new - CPL_MAX(nentries_old, 0));
  } else {
      cpl_msg_info(cpl_func, "(mf_atm      ) Current database is up to date");
  }

  return CPL_ERROR_NONE;
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief Check local GDAS database for two closest in time data files.
 *
 * @param gdasdir            Path of the GDAS DB.
 * @param lon                Longitude.
 * @param lat                Latitude.
 * @param year               Year.
 * @param month              Month.
 * @param day                Day.
 * @param utch               Is Hour in UTC or not.
 * @param timestamps         output: Updated array of time stamps of found files.
 * @param gdas_files         output: Updated array of found filenames, empty on failure.
 *
 * @return error code        1 on success, 0 on failure.
 *
 */
/* ---------------------------------------------------------------------------*/
static int check_local_db(
    const char               *gdasdir,
    const double             lon,
    const double             lat,
    const int                year,
    const int                month,
    const int                day,
    int                      utch,
    cpl_array                *timestamps,
    cpl_array                *gdas_files)
{

    double t1;
    double t2;

    char   str1[MF_MAXLEN];
    char   str2[MF_MAXLEN];

    if (   mf_gdas_get_closest_file(str1, &t1, gdasdir, lon, lat, year, month, day, utch,     -1)
        && mf_gdas_get_closest_file(str2, &t2, gdasdir, lon, lat, year, month, day, utch + 1, +1) ){

        cpl_array_set_float(timestamps, 0, t1);
        cpl_array_set_float(timestamps, 1, t2);

        cpl_msg_info(cpl_func, "(mf_atm      ) GDAS files:");
        cpl_msg_info(cpl_func, "(mf_atm      ) 1. %s", str1);
        cpl_msg_info(cpl_func, "(mf_atm      ) 2. %s", str2);

        if (cpl_array_get_float(timestamps, 1, NULL) < cpl_array_get_float(timestamps, 0, NULL)) {
            cpl_array_set_float(timestamps, 1, cpl_array_get_float(timestamps, 1, NULL) + 24.);
        }

        cpl_array_set_string(gdas_files, 0, str1);
        cpl_array_set_string(gdas_files, 1, str2);

        return 1;

    } else {

        return 0;
    }
}


cpl_error_code mf_atm_getgdas_auto(cpl_table **gdas_profile1,
                                   cpl_table **gdas_profile2,
                                   cpl_array *timestamps, const mfdrv *drvpar)
{
    /*!
     * Gets the two GDAS profiles that are closest in requested time. The
     * profiles are either taken from a local TAR archive or downloaded from
     * a webserver if required. If even the latter fails, an average profile
     * is provided. The routines also provides the time stamps for the two
     * output GDAS profiles and the requested time.
     *
     * \b INPUT:
     * \param drvpar         ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param gdas_profile1  CPL table with first GDAS profile
     * \param gdas_profile2  CPL table with second GDAS profile
     * \param timestamps     time stamps for GDAS profiles and requested time
     *
     * \b ERRORS:
     * - Illegal input
     * - File I/O
     * - Unexpected file structure
     * - see subroutines
     */

    cpl_errorstate err_state = cpl_errorstate_get();
    cpl_error_code err_code1 = CPL_ERROR_NONE,
                   err_code2 = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_array *gdas_files;
    char xtract_res[MF_MAXLEN] = "\0";
    char hh[] = "00", mm[] = "00", ss[] = "00", /* hour, minute, second */
         date[] = "00000000..",                 /* YYYYMMDD */
         basedir[MF_MAXLEN],
         gdasdir[MF_MAXLEN] = "\0",
         filename[MF_MAXLEN] = "\0",
         errtxt[MF_MAXLEN];
    int year, month, day;
    long h, m, s;
    double utch,     /* time of observation */
           date_obs, /* observing date */
           lon, lat, /* geolocation of telescope (longitude,latitude) */
           geoelev,  /* geoelevation of telescope */
           dummy;

    /* Get base directory */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));

    /* Get GDAS directory */
    p = cpl_parameterlist_find(drvpar->parlist, "gdas_dir");
    mf_basic_abspath(gdasdir, cpl_parameter_get_string(p), mf_get_datadir());

    /* Set name of output file for extract_grib */
    strcpy(xtract_res, gdasdir);
    strcat(xtract_res, "extract_grib.tmp");

    /* Get UTC in h */
    p = cpl_parameterlist_find(drvpar->parlist, "utc");
    utch = cpl_parameter_get_double(p) / 3600.;

    /* Decompose UTC into h, m, and s */

    h = (long) utch;
    if (h < 10) {
        sprintf(hh, "0%1li", h);
    } else if (h < 25) {
        sprintf(hh, "%2li", h);
    } else {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "UTC out of range");
    }

    m = (long)((utch - (double)h) * 60.);
    if (m < 10) {
        sprintf(mm, "0%1li", m);
    } else if (m < 60) {
        sprintf(mm, "%2li", m);
    } else {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "UTC out of range");
    }

    s = (long)(((utch - (double)h) * 60. - (double)m) * 60.);
    if (s < 10) {
        sprintf(ss, "0%1li", s);
    } else if (s < 60) {
        sprintf(ss, "%2li", s);
    } else {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "UTC out of range");
    }

    /* Extract longitude / latitude from parlist */
    p = cpl_parameterlist_find(drvpar->parlist, "longitude");
    lon = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find(drvpar->parlist, "latitude");
    lat = cpl_parameter_get_double(p);

    /* Get date */
    p = cpl_parameterlist_find(drvpar->parlist, "obsdate");
    date_obs = cpl_parameter_get_double(p);
    if ((dummy = mf_basic_fracyear2date(&year, &month, &day, NULL, NULL, NULL,
                                        &date_obs))){};
    sprintf(date, "%4i%02i%02i", year, month, day);

    /* Build GDAS file names and time stamps */

    cpl_array_set_float(timestamps, 0, utch);
    cpl_array_set_float(timestamps, 1, utch);
    cpl_array_set_float(timestamps, 2, utch);

    gdas_files = cpl_array_new(2, CPL_TYPE_STRING);

    /* Check if expected GDAS profiles already exist at expected location */
    cpl_msg_info(cpl_func, "Checking GDAS files for: "
                 "C%+3.1f%+3.1f %4i%02i%02i %g",
                 lon, lat, year, month, day, utch);

    if (!check_local_db(gdasdir, lon, lat, year, month, day, utch,
                   timestamps, gdas_files)) {
        cpl_msg_warning(cpl_func, "Files not present locally, searching for "
                        "newer version of gdas database on %s", base_gdas_url);
        update_local_db(gdasdir, lon, lat);
        check_local_db(gdasdir, lon, lat, year, month, day, utch,
                        timestamps, gdas_files);
    }

    err_state = cpl_errorstate_get();

    /* Read GDAS profiles */

    if (cpl_array_get_string(gdas_files, 0) != NULL) {
        err_code1 = mf_atm_readgdas(*gdas_profile1,
                                    cpl_array_get_string(gdas_files, 0), "m");
    } else {
        err_code1 = MF_ERROR_IO;
    }
    if (cpl_array_get_string(gdas_files, 1) != NULL) {
        err_code2 = mf_atm_readgdas(*gdas_profile2,
                                    cpl_array_get_string(gdas_files, 1), "m");
    } else {
        err_code2 = MF_ERROR_IO;
    }
    cpl_errorstate_set(err_state);

    if (err_code1 != CPL_ERROR_NONE && err_code2 != CPL_ERROR_NONE) {
        cpl_msg_warning(cpl_func, "Could not find GDAS profile on server");
        cpl_msg_warning(cpl_func, "Using average profile");

        /* No GDAS profile was found by grib
           -> Retrieve best matching profile from library */
        int lmonth = month / 2;
        lmonth++;
        if (lmonth > 6) {
            lmonth = 1;
        }
        mf_basic_initstring(filename, MF_MAXLEN);
        sprintf(filename, "%s/profiles/lib/GDAS_t0_s%1i.atm",
                mf_get_datadir(), lmonth);

        if (mf_atm_readgdas(*gdas_profile1, filename, "km") !=
            CPL_ERROR_NONE) {
            cpl_array_delete(gdas_files);
            return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                         "Could not find profile: %s",
                                         filename);
        }

        cpl_table_delete(*gdas_profile2);
        *gdas_profile2 = cpl_table_duplicate(*gdas_profile1);
    } else if (err_code1 != CPL_ERROR_NONE && err_code2 == CPL_ERROR_NONE) {
        /* gdas_profile1 does not exist */
        cpl_table_delete(*gdas_profile1);
        *gdas_profile1 = cpl_table_duplicate(*gdas_profile2);
    } else if (err_code1 == CPL_ERROR_NONE && err_code2 != CPL_ERROR_NONE) {
        /* gdas_profile2 does not exist */
        cpl_table_delete(*gdas_profile2);
        *gdas_profile2 = cpl_table_duplicate(*gdas_profile1);
    }

    /* Get geoelevation in m */
    p = cpl_parameterlist_find(drvpar->parlist, "geoelev");
    geoelev = cpl_parameter_get_double(p);

    /* Compare lowest height of GDAS profiles and geoelevation */
    if (cpl_table_get(*gdas_profile1, "height", 0, NULL) > geoelev) {
        sprintf(errtxt, "%s: %s (lowest height %g > geoelev %g of mfdrv "
                "*drvpar)", MF_ERROR_UFS_TXT,
                cpl_array_get_string(gdas_files, 0),
                cpl_table_get(*gdas_profile1, "height", 0, NULL), geoelev);
        cpl_array_delete(gdas_files);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }
    if (cpl_table_get(*gdas_profile2, "height", 0, NULL) > geoelev) {
        sprintf(errtxt, "%s: %s (lowest height %g > geoelev %g of mfdrv "
                "*drvpar)", MF_ERROR_UFS_TXT,
                cpl_array_get_string(gdas_files, 1),
                cpl_table_get(*gdas_profile2, "height", 0, NULL), geoelev);
        cpl_array_delete(gdas_files);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    /* Delete temporary array */
    cpl_array_delete(gdas_files);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_atm_getgdas_user(cpl_table *gdas_profile,
                                   const mfdrv *drvpar)
{
    /*!
     * Gets a user-provided GDAS-like profile, which has to consist of four
     * columns (P[hPa] HGT[m] T[K] RELHUM[%]) and a possible header. The file
     * path and name has to be provided by the parameter \e gdas_prof.
     *
     * \b INPUT:
     * \param drvpar        ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param gdas_profile  CPL table with GDAS-like profile
     *
     * \b ERRORS:
     * - Unexpected file structure
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    char filename[MF_MAXLEN] = "", basedir[MF_MAXLEN] = "";
    char relfilename[MF_MAXLEN] = "", errtxt[MF_MAXLEN];
    double geoelev = 0.;

    /* Get name of GDAS-like profile */
    p = cpl_parameterlist_find(drvpar->parlist, "gdas_prof");
    strncpy(filename, cpl_parameter_get_string(p), MF_MAXLEN);
    if (filename[0] != '/') {
        p = cpl_parameterlist_find(drvpar->parlist, "basedir");
        strncpy(basedir, cpl_parameter_get_string(p), MF_MAXLEN);
        sprintf(relfilename, "%s/%s", basedir, filename);
        mf_basic_absfile(filename, relfilename);
    }

    /* Write info message */
    cpl_msg_info(cpl_func, "Read user-provided GDAS-like profile %s",
                 filename);

    /* Read profile */
    if ((status = mf_atm_readgdas(gdas_profile, filename, "m")) !=
        CPL_ERROR_NONE) {
        return status;
    }

    /* Get geoelevation in m and add tolerance */
    p = cpl_parameterlist_find(drvpar->parlist, "geoelev");
    geoelev = cpl_parameter_get_double(p);
    geoelev += 0.1 * (cpl_table_get(gdas_profile, "height", 1, NULL) -
                      cpl_table_get(gdas_profile, "height", 0, NULL));

    /* Compare lowest height of GDAS profiles and geoelevation */
    if (cpl_table_get(gdas_profile, "height", 0, NULL) > geoelev) {
        sprintf(errtxt, "%s: %s (lowest height > geoelev of mfdrv "
                "*drvpar)", MF_ERROR_UFS_TXT, filename);
        return cpl_error_set_message(cpl_func, MF_ERROR_UFS, "%s", errtxt);
    }

    return CPL_ERROR_NONE;
}


void mf_atm_concatgdasname(char *gdasname, const char *gdasdir,
                           const double lon, const double lat,
                           const char *date, const char *hour)
{
    /*!
     * Concatenates a GDAS file name depending on the input path, longitude,
     * latitude, date string (8 characters), and hour string (2 characters)
     *
     * \b INPUT:
     * \param gdasdir   path to GDAS file
     * \param lon       geographic longitude in deg
     * \param lat       geographic latitude in deg
     * \param date      date string consisting of year, month, and day
     * \param hour      hour string
     *
     * \b OUTPUT:
     * \param gdasname  GDAS file name including path
     */

    sprintf(gdasname, "%sC%+.1f%+.1fD", gdasdir, lon, lat);
    strncat(gdasname, date, 4);
    strcat(gdasname, "-");
    strncat(gdasname, date + 4, 2);
    strcat(gdasname, "-");
    strncat(gdasname, date + 6, 2);
    strcat(gdasname, "T");
    strncat(gdasname, hour, 2);
    strcat(gdasname, ".gdas");
}

cpl_error_code mf_atm_readatm_fromFits(cpl_table **atm_profile, const char *atm_file)
{

    char* fitsfilename = cpl_sprintf("%s", atm_file);

    if(!(*atm_profile = cpl_table_load(fitsfilename, 1, 0))) {
        cpl_error_code status = cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                     "Could not open fits file: %s",
                                     fitsfilename);;
        cpl_free(fitsfilename);
        return status;
    }

    cpl_free(fitsfilename);

    return CPL_ERROR_NONE;

}

cpl_error_code mf_atm_readatm(cpl_table *atm_profile, const char *atm_file)
{
    /*!
     * \brief
     *   Read atmospheric standard profile.
     *
     * This function reads an atmospheric standard profile into a
     * \c CPL_TABLE. The output \c CPL_TABLE must exist before calling this
     * routine. It will get resized and overwritten.
     *
     * \b INPUT:
     * \param atm_file     filename of an atmospheric standard profile.
     *
     * \b OUTPUT:
     * \param atm_profile  \c CPL_TABLE with standard profile.
     * \return             CPL_ERROR_NONE on success,
     *                     MF_ERROR_FOPEN on error.
     */

    FILE *fp;                        /* file pointer */

    char line[MF_MAXLEN];         /* line in setup file */
    char str[MF_MAXLEN];
    char *ptr, *d;

    long nrows = 0;                  /* number of rows in profile */

    long i = 036;

    mf_basic_initstring(line, MF_MAXLEN);
    mf_basic_initstring(str, MF_MAXLEN);

    /* open file for reading */
    fp = fopen(atm_file, "r");
    if (fp == NULL) {
        return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                     "Could not open atm_file: %s", atm_file);
    } else {
        while (fgets(line, MF_MAXLEN - 1, fp) != NULL) {
            mf_basic_strtrim_inplace(line);                    /* trim line */
            mf_basic_rmcntrl_inplace(line);    /* remove control characters */

            if (line[0] == '!') {
                continue;                             /* skip over comments */
            }

            /* read number of profile levels */
            /*ptr = strstr(line, "! Profile Levels");*/
            ptr = strstr(line, "Levels");
            if (ptr != NULL) {
                strncpy(str, line, ptr-&line[0]);
                mf_basic_strtrim_inplace(str);
                nrows = atol(str);
                break;                /* number of levels is known -> break */
            }
            if ((d = fgets(line, MF_MAXLEN - 1, fp) )){};
        }

        cpl_table_set_size(atm_profile, nrows);

        while (fgets(line, MF_MAXLEN - 1, fp) != NULL) {
            mf_basic_strtrim_inplace(line);                    /* trim line */
            mf_basic_rmcntrl_inplace(line);    /* remove control characters */

            if (line[0] == '*') {                  /* block with new column */
                ptr = strchr(line, ' ');        /* find first blank in line */
                if (ptr == NULL) {
                    break;                   /* no blank found -> new block */
                }

                i = ptr-&line[1];    /* extract column name (separate unit) */
                strncpy(str, &line[1], i);
                str[i] = '\0';
                /* ensure that all characters in string are upper case */
                for (i = 0; str[i]; i++) {
                    str[i] = toupper(str[i]);
                }
                cpl_table_new_column(atm_profile, str, CPL_TYPE_FLOAT);

                i = 0;                            /* read new block of data */
                while (fgets(line, MF_MAXLEN - 1, fp) != NULL) {
                    mf_basic_strtrim_inplace(line);            /* trim line */
                    mf_basic_rmcntrl_inplace(line); /* remove control char. */

                    /* values are separated by blank or tab */
                    ptr = strtok(line, " \t");
                    while (ptr != 0) {
                        /* extract float values & write to table */
                        cpl_table_set_float(atm_profile, str, i, atof(ptr));
                        ptr = strtok(NULL, " \t");
                        i++;
                        if (i >= nrows) {
                            break;          /* last value in block -> break */
                        }
                    }
                    if (i >= nrows) {
                        break;              /* last value in block -> break */
                    }
                }
            }
        }

        fclose(fp);

        return CPL_ERROR_NONE;
    }
}


cpl_error_code mf_atm_readgdas(cpl_table *gdas_profile, const char *gdas_file,
                               const char *hgt_units)
{
    /*!
     * \brief
     *   Read GDAS profile.
     *
     * This function reads an atmospheric GDAS profile into a \c CPL_TABLE.
     * The output \c CPL_TABLE must exist before calling this routine. It will
     * get resized and overwritten.
     *
     * \b INPUT:
     * \param gdas_file     filename of an atmospheric GDAS profile.
     * \param hgt_units     height unit ("m" or "km")
     *
     * \b OUTPUT:
     * \param gdas_profile  \c CPL_TABLE with GDAS profile.
     * \return              CPL_ERROR_NONE on success,
     *                      MF_ERROR_FOPEN on error.
     */

    FILE *fp;                  /* file pointer */

    char line[MF_MAXLEN];   /* line in setup file */
    char *ptr, *d;

    long nrows = 0, nhead = 0;

    int i = 0;

    float vals[4] = { 0., 0., 0., 0. };

    mf_basic_initstring(line, MF_MAXLEN);

    /* open file for reading */
    fp = fopen(gdas_file, "r");
    if (fp == NULL) {
        return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                     "Could not open gdas_file: %s",
                                     gdas_file);
    } else {
        /* Find header (non-data lines at beginning of file) */
        while (fgets(line, MF_MAXLEN - 1, fp) != NULL) {
            mf_basic_rmcntrl_inplace(line);
            mf_basic_strtrim_inplace(line);
            ptr = strtok(line, " \t");
            for (i = 0; i < 4; i++) {
                if (ptr == NULL) break;
                vals[i] = atof(ptr);
                ptr = strtok(NULL, " \t");
            }
            if (vals[0] + vals[1] + vals[2] + vals[3] == 0) {
                nhead++;
            } else {
                // first data line found
                break;
            }
        }
        rewind(fp);

        /* Skip header lines */
        for (i = 0; i < nhead; i++) {
            if ((d = fgets(line, MF_MAXLEN - 1, fp))){};
        }

        /* Count data lines (excluding comments) */
        while (fgets(line, MF_MAXLEN - 1, fp) != NULL) {
            mf_basic_rmcntrl_inplace(line);
            mf_basic_strtrim_inplace(line);
            if (line[0] != '#') {
                nrows++;
            }
        }

        cpl_table_set_size(gdas_profile, nrows);

        cpl_table_new_column(gdas_profile, "press", CPL_TYPE_FLOAT);
        cpl_table_new_column(gdas_profile, "height", CPL_TYPE_FLOAT);
        cpl_table_new_column(gdas_profile, "temp", CPL_TYPE_FLOAT);
        cpl_table_new_column(gdas_profile, "relhum", CPL_TYPE_FLOAT);

        rewind(fp);

        /* Skip header lines */
        for (i = 0; i < nhead; i++) {
            d = fgets(line, MF_MAXLEN - 1, fp);
        }

        i = 0;
        while (fgets(line, MF_MAXLEN - 1, fp) != NULL) {
            mf_basic_strtrim_inplace(line);

            /* skip comments */
            if (line[0] == '#') {
                continue;
            }

            ptr = strtok(line, " \t");
            cpl_table_set_float(gdas_profile, "press", i, atof(ptr));

            /* height in km */
            ptr = strtok(NULL, " \t");
            if (strcmp(hgt_units, "m") == 0) {
                cpl_table_set_float(gdas_profile, "height", i,
                                    atof(ptr)/1000.);
            } else {
                cpl_table_set_float(gdas_profile, "height", i, atof(ptr));
            }

            ptr = strtok(NULL, " \t");
            cpl_table_set_float(gdas_profile, "temp", i, atof(ptr));

            ptr = strtok(NULL, " \t");
            cpl_table_set_float(gdas_profile, "relhum", i, atof(ptr));
            i++;
        }

        fclose(fp);

        return CPL_ERROR_NONE;
    }
}


cpl_error_code mf_atm_interpolprofile(cpl_table *outprofile,
                                      const cpl_table *profile1,
                                      const cpl_table *profile2,
                                      const cpl_array *timestamps)
{
    /*!
     * \brief
     *   Interpolate profiles linearly between two points in time.
     *
     * This function calculates an interpolated profile from two input
     * profiles taken at different times. \em profile1 and \em profile2 are
     * expected to be sampled over the same geoelevation grid and to contain
     * the same data columns (molecules). The output profile \em outprofile is
     * linearly interpolated in time. If the number of layers in both profiles
     * differs, the profile closest to the requested time is taken as output
     * profile. In this case, a warning message is printed.
     *
     * \note
     *   No checks are performed whether the input profiles contain the same
     *   column information.
     *
     * \b INPUT:
     * \param profile1    first input atmospheric profile.
     * \param profile2    second input atmospheric profile.
     * \param timestamps  time stamps for input and output profiles.
     *
     * \b OUTPUT:
     * \param outprofile  interpolated output atmospheric profile.
     *
     * \b ERRORS:
     * - Invalid object value(s)
     * - Invalid object structure
     */

    cpl_array *cols, *col1, *col2;
    char errtxt[MF_MAXLEN] = "";
    const char *col_name;
    int i = 0;
    long nrow1 = 0, nrow2 = 0, ncol = 0;
    float t1 = 0., t2 = 0., tout = 0., dt = 0.;

    /* Get time values */
    t1 = cpl_array_get_float(timestamps, 0, NULL);
    t2 = cpl_array_get_float(timestamps, 1, NULL);
    tout = cpl_array_get_float(timestamps, 2, NULL);

    /* Check time values */
    if (t1 > t2 || tout < t1 || tout > t2) {
        sprintf(errtxt, "%s: cpl_array *timestamps (invalid time(s))",
                MF_ERROR_IOV_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOV, "%s", errtxt);
    }

    /* Ensure that input tables have the same structure */
    if (cpl_table_compare_structure(profile1, profile2) == 1) {
        sprintf(errtxt, "%s: cpl_table *profile1 != *profile2 (columns)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get relative position ([0,1]) of output profile time */
    if (t1 == t2) {
        dt = 0;
    } else {
        dt = tout - t1;
        dt /= t2 - t1;
    }

    /* Get number of rows in both profiles */
    nrow1 = cpl_table_get_nrow(profile1);
    nrow2 = cpl_table_get_nrow(profile2);

    /* Check number of rows -> Difference: Do not interpolate and take profile
       closest to requested time */
    if (nrow1 != nrow2) {
        cpl_msg_warning(cpl_func, "GDAS profiles differ in # of rows -> "
                        "no interpolation -> take profile closest in time");
        cpl_table_delete(outprofile);
        if (dt < 0.5) {
            outprofile = cpl_table_duplicate(profile1);
        } else {
            outprofile = cpl_table_duplicate(profile2);
        }
        return CPL_ERROR_NONE;
    }

    /* Get number and names of columns */
    ncol = cpl_table_get_ncol(profile1);
    cols = cpl_table_get_column_names(profile1);

    /* Prepare output table structure */
    cpl_table_copy_structure(outprofile, profile1);
    cpl_table_set_size(outprofile, nrow1);

    /* Create temporary arrays for column values */
    col1 = cpl_array_new(nrow1, CPL_TYPE_FLOAT);
    col2 = cpl_array_new(nrow1, CPL_TYPE_FLOAT);

    /* Interpolate data values of each column depending on time stamps */
    for (i = 0; i < ncol; i++) {
        col_name = cpl_array_get_string(cols, i);
        cpl_array_copy_data_float(col1,
                                  cpl_table_get_data_float_const(profile1,
                                                                 col_name));
        cpl_array_copy_data_float(col2,
                                  cpl_table_get_data_float_const(profile2,
                                                                 col_name));
        cpl_array_subtract(col2, col1);
        cpl_array_multiply_scalar(col2, dt);
        cpl_array_add(col2, col1);
        cpl_table_copy_data_float(outprofile, col_name,
                                  cpl_array_get_data_float(col2));
    }

    /* Delete temporary arrays */
    cpl_array_delete(cols);
    cpl_array_delete(col1);
    cpl_array_delete(col2);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_atm_convertgdas(cpl_table *merged_profile,
                                  const cpl_table *atm_profile,
                                  const cpl_table *gdas_profile,
                                  const cpl_array *molecs,
                                  const mfdrv *drvpar)
{
    /*!
     * Merges an atmospheric standard and a GDAS-like profile. Depending on
     * the parameter \e layers, either a fixed (= 1) or a natural (= 0) grid
     * of heights is used for the profile.
     *
     * \b INPUT:
     * \param atm_profile     atmospheric standard profile
     * \param gdas_profile    GDAS profile
     * \param molecs          list of molecules
     * \param drvpar          ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param merged_profile  merged profile
     *
     * \b ERRORS:
     * - see subroutines
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    int layers = 0;
    double geoelev = 0.;

    /* Get flag for grid of layers */
    p = cpl_parameterlist_find(drvpar->parlist, "layers");
    layers = cpl_parameter_get_int(p);

    /* Get geoelevation of observing site in km */
    p = cpl_parameterlist_find(drvpar->parlist, "geoelev");
    geoelev = cpl_parameter_get_double(p) / 1000.;

    /* Perform profile merging */
    if (layers == 1) {
        /* Merging for fixed grid of heights */
        cpl_msg_info(cpl_func, "Take fixed grid of layers for merged "
                     "profile");
        status = mf_atm_convertgdas_fixed(merged_profile, atm_profile,
                                          gdas_profile, molecs, geoelev);
    } else {
        /* Get natural grid of heights */
        mf_atm_mergelayers(merged_profile, atm_profile, gdas_profile,
                           drvpar);
        /* Merging for natural grid of heights */
        cpl_msg_info(cpl_func, "Take natural grid of layers for merged "
                     "profile");
        status = mf_atm_mergegdas(merged_profile, atm_profile, gdas_profile,
                                  molecs);
    }

    /* Treat errors */
    if (status != CPL_ERROR_NONE) {
        return status;
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_atm_convertgdas_fixed(cpl_table *merged_profile,
                                        const cpl_table *atm_profile,
                                        const cpl_table *gdas_profile,
                                        const cpl_array *molecs,
                                        const double geoelev)
{
    /*!
     * \brief
     *   Merge atmospheric standard and GDAS profiles.
     *
     * Information contained in GDAS profiles is incorporated into the
     * standard atmospheric profile. The resulting merged profile contains all
     * columns from the standard profile and has pressure, temperature and
     * H2O columns replaced with the GDAS values. The merged profile is
     * interpolated to a new irregular height grid with 50 levels, the first
     * 29 of which come from GDAS. The four height levels from 20 to 26 km are
     * a weighted mix of GDAS and standard profile. The influence of the GDAS
     * profile decreases with increasing height: 80%, 60%, 40%, 20% at 20km,
     * 22km, 24km, 26km. Beyond 26km no GDAS information is available.
     *
     * \note
     *   GDAS does not contain values for molecules other than H2O. At heights
     *   <26km only standard profile information is available for all other
     *   molecules.
     *
     * \b INPUT:
     * \param atm_profile     atmospheric standard profile
     * \param gdas_profile    GDAS profile
     * \param molecs          list of molecules
     * \param geoelev         geoelevation in km
     *
     * \b OUTPUT:
     * \param merged_profile  merged profile
     * \return                CPL_ERROR_NONE on success,
     *                        MF_ERROR_SUBROUTINE on error.
     */

    cpl_error_code err_code;

    int gdas_nrows, merged_nrows,            /* number of rows in profiles */
        i, j,                                /* loop variables */
        n_molecs,                            /* number of molecules */
        n_gdas_max;                          /* max number of GDAS levels */

    float step;                              /* spacing of upper levels */

    int n_low_hgt = 33,                      /* number of levels in low_hgt */
        i0,                                  /* number of lowest layer */
        exh2o = 0;                           /* existence of H2O column */

    double low_hgt[] = {  0, 0.5,   1, 1.5,
                          2, 2.5,   3, 3.5,   4, 4.5,   5, 5.5,   6, 6.5,
                          7, 7.5,   8, 8.5,   9, 9.5,  10,  11,  12,  13,
                         14,  15,  16,  17,  18,  20,  22,  24,  26},
           *hgt_levels,
           dh = 1.,                          /* height interval */
           t,                                /* temperature */
           p,                                /* pressure */
           h,                                /* height */
           rel_hum,                          /* relative humidity */
           ppmv = 0,                         /* volume mixing ratio */
           val_overlap = 0, val_merged = 0;

    char *mol;                               /* string for molecule */

    cpl_table *overlap_region,    /* to store GDAS values in overlap region */
              *tmp_gdas;                  /* temporarily store GDAS rel_hum */
    /*
     *  all information is in the ATM profile
     *  to introduce some variability at lower altitudes GDAS overwrites
     *  the following parameters:
     *  height, pressure, temperature, H2O via relative humidity
     */

    /*
     *  column names in GDAS:
     *  press height temp relhum
     *  column names in ATM:
     *  HGT PRE TEM N2 O2 CO2 O3 H2O CH4 N2O HNO3 OH CO NO2 N2O5 ClO HOCl
     *  ClONO2 NO HNO4 HCN NH3 F11 F12 F14 F22 CCl4 COF2 H2O2 C2H2 C2H6
     *  OCS SO2 SF6
     */

    /* find lowest hgt level */
    for (i0 = 0, i = 0; i < n_low_hgt; i++) {
         if (low_hgt[i] > geoelev - dh) {
            i0 = i;
            break;
        }
    }

    /* get the number of molecules */
    n_molecs = cpl_array_get_size(molecs);

    /* string for a single molecule */
    mol = (char *)malloc(MF_MAXLEN);

    /* get numbers of rows for GDAS, ATM and merged profile */
    gdas_nrows = cpl_table_get_nrow(gdas_profile);
    merged_nrows = cpl_table_get_nrow(merged_profile);

    /* ensure that output table has (at least) 50 rows for Cerro Paranal */
    if (merged_nrows < 54-i0) {
        cpl_table_set_size(merged_profile, 54-i0);
        merged_nrows = 54-i0;
    }

    /* define array for height levels */
    hgt_levels = (double *) malloc(merged_nrows * sizeof(double));

    /* fill hgt vector with values for GDAS range */
    for (i = 0; i < n_low_hgt-i0; i++) {
        hgt_levels[i] = low_hgt[i+i0];
    }

    /* calculate step size for higher hgt levels */
    step = (120 - low_hgt[n_low_hgt-1]) / (merged_nrows - (n_low_hgt-i0));

    /* fill hgt vector with higher level values */
    for (i = n_low_hgt-i0; i < merged_nrows; i++) {
        hgt_levels[i] = low_hgt[n_low_hgt-1] + (i-n_low_hgt+i0+1) * step;
    }

    /* get number of GDAS levels below height of 20 km in new grid */
    n_gdas_max = n_low_hgt-i0-4;

    /* insert new columns including unit and format */

    /* height in [km] */
    cpl_table_new_column(merged_profile, "HGT", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(merged_profile, "HGT", "km");
    cpl_table_set_column_format(merged_profile, "HGT", "%10.3e");

    /* pressure in [mb] */
    cpl_table_new_column(merged_profile, "PRE", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(merged_profile, "PRE", "mb");
    cpl_table_set_column_format(merged_profile, "PRE", "%10.3e");

    /* temperature in [k] */
    cpl_table_new_column(merged_profile, "TEM", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(merged_profile, "TEM", "K");
    cpl_table_set_column_format(merged_profile, "TEM", "%10.3e");

    /* loop over all molecules */
    for (i = 0; i < n_molecs; i++) {
        strcpy(mol, cpl_array_get_string(molecs, i));

        cpl_table_new_column(merged_profile, mol, CPL_TYPE_FLOAT);
        cpl_table_set_column_unit(merged_profile, mol, "ppmv");
        cpl_table_set_column_format(merged_profile, mol, "%10.3e");

        /* existence of H2O column */
        if (strcmp(mol, "H2O") == 0) {
            exh2o = 1;
        }
    }

    /* write height column */
    for (i = 0; i < merged_nrows; i++) {
        cpl_table_set_float(merged_profile, "HGT", i, *(hgt_levels + i));
    }

    free(hgt_levels);

    err_code = CPL_ERROR_NONE;

    /* interpolate pressure column (GDAS range) */
    err_code += mf_basic_interpolcolumn(gdas_profile, "height", "press",
                                        merged_profile, "HGT", "PRE", 0,
                                        n_gdas_max+4);

    /* interpolate temperature column (GDAS range) */
    err_code += mf_basic_interpolcolumn(gdas_profile, "height", "temp",
                                        merged_profile, "HGT", "TEM", 0,
                                        n_gdas_max+4);

    if (err_code != CPL_ERROR_NONE) {
        free(mol);
        return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                     "Problem in interpolation");
    }

    /* interpolate H2O column (GDAS range) */

    if (exh2o == 1) {

        /* new temporary table for calculating the overlap region */
        tmp_gdas = cpl_table_new(gdas_nrows);

        /* height in [km] */
        cpl_table_new_column(tmp_gdas, "height", CPL_TYPE_FLOAT);
        cpl_table_set_column_unit(tmp_gdas, "height", "km");
        cpl_table_set_column_format(tmp_gdas, "height", "%10.3e");

        /* H2O */
        cpl_table_new_column(tmp_gdas, "H2O", CPL_TYPE_FLOAT);
        cpl_table_set_column_unit(tmp_gdas, "H2O", "ppmv");
        cpl_table_set_column_format(tmp_gdas, "H2O", "%10.3e");

        /* convert relative humidity from GDAS to ppmv (stored in tmp_gdas) */
        for (i = 0; i < gdas_nrows; i++) {
            h = cpl_table_get_float(gdas_profile, "height", i, NULL);
            rel_hum = cpl_table_get_float(gdas_profile, "relhum", i, NULL);
            t = cpl_table_get_float(gdas_profile, "temp", i, NULL);
            p = cpl_table_get_float(gdas_profile, "press", i, NULL);
            mf_basic_rhum2ppmv(&ppmv, &t, &p, &rel_hum);
            cpl_table_set_float(tmp_gdas, "H2O", i, ppmv);
            cpl_table_set_float(tmp_gdas, "height", i, h);
        }

        /* interpolate H2O column (overlap region) */
        err_code = mf_basic_interpolcolumn(tmp_gdas, "height", "H2O",
                                           merged_profile, "HGT", "H2O", 0,
                                           n_gdas_max+4);

        cpl_table_delete(tmp_gdas);

        if (err_code != CPL_ERROR_NONE) {
            free(mol);
            return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                         "Problem in interpolation");
        }

    }

    /* new temporary table for calculating the overlap region */
    overlap_region = cpl_table_new(merged_nrows);

    /* pressure in [mb] */
    cpl_table_new_column(overlap_region, "PRE", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(overlap_region, "PRE", "mb");
    cpl_table_set_column_format(overlap_region, "PRE", "%10.3e");

    /* temperature in [K] */
    cpl_table_new_column(overlap_region, "TEM", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(overlap_region, "TEM", "K");
    cpl_table_set_column_format(overlap_region, "TEM", "%10.3e");

    /* H2O */
    if (exh2o == 1) {
        cpl_table_new_column(overlap_region, "H2O", CPL_TYPE_FLOAT);
        cpl_table_set_column_unit(overlap_region, "H2O", "ppmv");
        cpl_table_set_column_format(overlap_region, "H2O", "%10.3e");
    }

    /* set TEM, PRE, & H2O columns */
    for (i = 0; i < n_gdas_max+4; i++) {
        p = cpl_table_get_float(merged_profile, "PRE", i, NULL);
        cpl_table_set_float(overlap_region, "PRE", i, p);
        t = cpl_table_get_float(merged_profile, "TEM", i, NULL);
        cpl_table_set_float(overlap_region, "TEM", i, t);
        if (exh2o == 1) {
            ppmv = cpl_table_get_float(merged_profile, "H2O", i, NULL);
            cpl_table_set_float(overlap_region, "H2O", i, ppmv);
        }
    }

    /* interpolate pressure column (ATM range) */
    err_code += mf_basic_interpolcolumn(atm_profile, "HGT", "PRE",
                                        merged_profile, "HGT", "PRE",
                                        n_gdas_max, merged_nrows);

    /* interpolate temperature column (ATM range) */
    err_code += mf_basic_interpolcolumn(atm_profile, "HGT", "TEM",
                                        merged_profile, "HGT", "TEM",
                                        n_gdas_max, merged_nrows);

    /* insert molecules if requested */
    for (i = 0; i < n_molecs; i++) {
        strcpy(mol, cpl_array_get_string(molecs, i));

        /* insert H2O only above GDAS data, rest from 0km */
        if (strcmp(mol, "H2O") == 0) {
            err_code += mf_basic_interpolcolumn(atm_profile, "HGT", "H2O",
                                                merged_profile, "HGT", "H2O",
                                                n_gdas_max, merged_nrows);
        } else {
            err_code += mf_basic_interpolcolumn(atm_profile, "HGT", mol,
                                                merged_profile, "HGT", mol,
                                                0, merged_nrows);
        }
    }

    free(mol);

    /* insert overlap region */
    for (i = n_gdas_max, j = 1; i < n_gdas_max+4; i++, j++) {
        val_merged = cpl_table_get_float(merged_profile, "PRE", i, NULL);
        val_overlap = cpl_table_get_float(overlap_region, "PRE", i, NULL);
        cpl_table_set_float(merged_profile, "PRE", i,
                            val_merged*0.2*j+val_overlap*(1-0.2*j));
        val_merged = cpl_table_get_float(merged_profile, "TEM", i, NULL);
        val_overlap = cpl_table_get_float(overlap_region, "TEM", i, NULL);
        cpl_table_set_float(merged_profile, "TEM", i,
                            val_merged*0.2*j+val_overlap*(1-0.2*j));
        if (exh2o == 1) {
            val_merged = cpl_table_get_float(merged_profile, "H2O", i, NULL);
            val_overlap = cpl_table_get_float(overlap_region, "H2O", i, NULL);
            cpl_table_set_float(merged_profile, "H2O", i,
                                val_merged*0.2*j+val_overlap*(1-0.2*j));
        }
    }
    cpl_table_delete(overlap_region);

    if (err_code != CPL_ERROR_NONE) {
        return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                     "Problem in interpolation");
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_atm_mergelayers(cpl_table *merged_profile,
                                  const cpl_table *atm_profile,
                                  const cpl_table *gdas_profile,
                                  const mfdrv *drvpar)
{
    /*!
     * Merges the height layers of an atmospheric standard and a GDAS-like
     * profile. All layers of both profiles are combined and written into the
     * "HGT" column of the output profile. If local meteo data are considered,
     * the geoelevation will added as well.
     *
     * \b INPUT:
     * \param atm_profile     atmospheric standard profile
     * \param gdas_profile    GDAS profile
     * \param drvpar          ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param merged_profile  CPL table with grid of layers in km in column
     *                        "HGT"
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    int n_atm = 0, n_gdas = 0, n_max = 0, i_atm = 0, i_gdas = 0, i_obs = 0;
    int i = 0, j = 0;
    const float *hgt_atm, *hgt_gdas;
    float *hgt_merged;
    double geoelev = 0., emix = 0., dh = 1., h_min = 0., h0 = 0., h = 0.;

    /* Check existence of HGT column in output table */
    if (cpl_table_has_column(merged_profile, "HGT") != 1) {
        cpl_table_new_column(merged_profile, "HGT", CPL_TYPE_FLOAT);
        cpl_table_set_column_unit(merged_profile, "HGT", "km");
        cpl_table_set_column_format(merged_profile, "HGT", "%10.3e");
    }

    /* Get number of layers */
    n_atm = cpl_table_get_nrow(atm_profile);
    n_gdas = cpl_table_get_nrow(gdas_profile);

    /* Get maximum number of layers */
    n_max = n_atm + n_gdas;

    /* Get geoelevation of observing site in km */
    p = cpl_parameterlist_find(drvpar->parlist, "geoelev");
    geoelev = cpl_parameter_get_double(p) / 1000.;

    /* Get upper mixing height for consideration of meteo station
       parameters */
    p = cpl_parameterlist_find(drvpar->parlist, "emix");
    emix = cpl_parameter_get_double(p);

    /* If emix > gelev, add geoelevation for local meteo data */
    if (emix > geoelev) {
        n_max += 1;
    }

    /* Set size of output table */
    cpl_table_set_size(merged_profile, n_max);

    /* Get lowest layer height */
    h_min = geoelev - dh;
    h0 = h_min;

    /* Get pointers to height columns */
    hgt_atm = cpl_table_get_data_float_const(atm_profile, "HGT");
    hgt_gdas = cpl_table_get_data_float_const(gdas_profile, "height");

    /* Initialise output grid column and get pointer to it */
    cpl_table_fill_column_window_float(merged_profile, "HGT", 0, n_max, 0.);
    hgt_merged = cpl_table_get_data_float(merged_profile, "HGT");

    /* Build grid */
    for (j = 0, i_atm = 0, i_gdas = 0, i_obs = 0, i = 0; i < n_max; i++) {
        if (emix > geoelev && i_obs == 0 && geoelev < hgt_atm[i_atm] &&
            geoelev < hgt_gdas[i_gdas]) {
            h = geoelev;
            i_obs++;
        } else {
            if (i_atm < n_atm &&
                (i_gdas >= n_gdas || hgt_atm[i_atm] < hgt_gdas[i_gdas])) {
                h = hgt_atm[i_atm];
                i_atm++;
            } else {
                h = hgt_gdas[i_gdas];
                i_gdas++;
            }
        }
        if (h > h0) {
            hgt_merged[j] = h;
            h0 = h;
            j++;
        }
    }

    /* Resize output table */
    cpl_table_set_size(merged_profile, j);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_atm_mergegdas(cpl_table *merged_profile,
                                const cpl_table *atm_profile,
                                const cpl_table *gdas_profile,
                                const cpl_array *molecs)
{
    /*!
     * Merges an atmospheric standard and a GDAS-like profile for a given grid
     * of layers (column "HGT"). The resulting merged profile contains all
     * columns from the standard profile and has pressure, temperature, and
     * H2O columns replaced with the GDAS values. Above the uppermost GDAS
     * layer, the relative deviation between the standard and the GDAS profile
     * (as measured for the highest GDAS layer) is gradually decreased to be
     * zero at a 1 + ::MF_MERGEFRAC times higher altitude. Above that limiting
     * height the unmodified standard profile is written into the output
     * profile.
     *
     * \b INPUT:
     * \param merged_profile  CPL table with grid of layers in km in column
     *                        "HGT"
     * \param atm_profile     atmospheric standard profile
     * \param gdas_profile    GDAS profile
     * \param molecs          list of molecules
     *
     * \b OUTPUT:
     * \param merged_profile  merged profile
     *
     * \b ERRORS:
     * - Invalid object structure
     * - Error in subroutine
     */

    /*
     *  Column names in GDAS:
     *  press height temp relhum
     *  Column names in ATM:
     *  HGT PRE TEM N2 O2 CO2 O3 H2O CH4 N2O HNO3 OH CO NO2 N2O5 ClO HOCl
     *  ClONO2 NO HNO4 HCN NH3 F11 F12 F14 F22 CCl4 COF2 H2O2 C2H2 C2H6
     *  OCS SO2 SF6
     */

    cpl_error_code err_code = CPL_ERROR_NONE;
    cpl_table *tmp_gdas, *overlap_region;
    char *mol;
    char errtxt[MF_MAXLEN];
    int gdas_nrows = 0, merged_nrows = 0, n_gdas_max = 0, n_merge = 0, i = 0;
    int n_molecs = 0, exh2o = 0;
    float *hgt_levels;
    double rel_hum = 0., h_gdas_max = 0., h_merge = 0., dh = 0., h = 0.;
    double p = 0., t = 0., ppmv = 0., hfrac = 0., val_merged = 0;
    double val_overlap = 0, dp = 0., dt = 0, dh2o = 0.;

    /* Check existence of HGT column */
    if (cpl_table_has_column(merged_profile, "HGT") != 1) {
        sprintf(errtxt, "%s: cpl_table *merged_profile (no HGT column)",
                MF_ERROR_IOS_TXT);
        return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s", errtxt);
    }

    /* Get numbers of rows for GDAS and merged profile */
    gdas_nrows = cpl_table_get_nrow(gdas_profile);
    merged_nrows = cpl_table_get_nrow(merged_profile);

    /* Get pointer to HGT column */
    hgt_levels = cpl_table_get_data_float(merged_profile, "HGT");

    /* Get critical heights for profile merging */
    h_gdas_max = cpl_table_get_float(gdas_profile, "height", gdas_nrows-1,
                                     NULL);
    h_merge = h_gdas_max * (1. + MF_MERGEFRAC);
    dh = h_merge - h_gdas_max;

    /* Loop over new grid and find critical rows */
    for (n_gdas_max = -1, n_merge = -1, i = 0; i < merged_nrows; i++) {
        if (n_gdas_max < 0 && hgt_levels[i] > h_gdas_max) {
            n_gdas_max = i;
        }
        if (n_merge < 0 && hgt_levels[i] > h_merge) {
            n_merge = i;
        }
    }

    /* Insert new columns including unit and format */

    /* height in [km] (must exist) */
    cpl_table_set_column_unit(merged_profile, "HGT", "km");
    cpl_table_set_column_format(merged_profile, "HGT", "%10.3e");

    /* pressure in [mb] */
    cpl_table_new_column(merged_profile, "PRE", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(merged_profile, "PRE", "mb");
    cpl_table_set_column_format(merged_profile, "PRE", "%10.3e");

    /* temperature in [k] */
    cpl_table_new_column(merged_profile, "TEM", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(merged_profile, "TEM", "K");
    cpl_table_set_column_format(merged_profile, "TEM", "%10.3e");

    /* Get the number of molecules */
    n_molecs = cpl_array_get_size(molecs);

    /* String for a single molecule */
    mol = (char *) malloc(MF_MAXLEN);

    /* Loop over all molecules */
    for (i = 0; i < n_molecs; i++) {
        strcpy(mol, cpl_array_get_string(molecs, i));
        cpl_table_new_column(merged_profile, mol, CPL_TYPE_FLOAT);
        cpl_table_set_column_unit(merged_profile, mol, "ppmv");
        cpl_table_set_column_format(merged_profile, mol, "%10.3e");
        /* Existence of H2O column */
        if (strcmp(mol, "H2O") == 0) {
            exh2o = 1;
        }
    }

    /* Interpolate pressure column (GDAS range) */
    err_code += mf_basic_interpolcolumn(gdas_profile, "height", "press",
                                        merged_profile, "HGT", "PRE", 0,
                                        n_gdas_max);

    /* Interpolate temperature column (GDAS range) */
    err_code += mf_basic_interpolcolumn(gdas_profile, "height", "temp",
                                        merged_profile, "HGT", "TEM", 0,
                                        n_gdas_max);

    /* Check for errors */
    if (err_code != CPL_ERROR_NONE) {
        free(mol);
        return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                     "Problem in interpolation");
    }

    /* Interpolate H2O column (GDAS range) */

    if (exh2o == 1) {

        /* New temporary table for calculating H2O in ppmv */
        tmp_gdas = cpl_table_new(gdas_nrows);
        cpl_table_new_column(tmp_gdas, "height", CPL_TYPE_FLOAT); // km
        cpl_table_new_column(tmp_gdas, "H2O", CPL_TYPE_FLOAT); // ppmv

        /* Convert relative humidity from GDAS to ppmv (stored in tmp_gdas) */
        for (i = 0; i < gdas_nrows; i++) {
            h = cpl_table_get_float(gdas_profile, "height", i, NULL);
            rel_hum = cpl_table_get_float(gdas_profile, "relhum", i, NULL);
            t = cpl_table_get_float(gdas_profile, "temp", i, NULL);
            p = cpl_table_get_float(gdas_profile, "press", i, NULL);
            mf_basic_rhum2ppmv(&ppmv, &t, &p, &rel_hum);
            cpl_table_set_float(tmp_gdas, "H2O", i, ppmv);
            cpl_table_set_float(tmp_gdas, "height", i, h);
        }

        /* Interpolate H2O column (GDAS range) */
        err_code = mf_basic_interpolcolumn(tmp_gdas, "height", "H2O",
                                           merged_profile, "HGT", "H2O", 0,
                                           n_gdas_max);

        /* Delete temporary table */
        cpl_table_delete(tmp_gdas);

        /* Check for errors */
        if (err_code != CPL_ERROR_NONE) {
            free(mol);
            return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                         "Problem in interpolation");
        }

    }

    /* New temporary table for storing the GDAS values and calculating the
       overlap region */
    overlap_region = cpl_table_new(merged_nrows);
    cpl_table_new_column(overlap_region, "PRE", CPL_TYPE_FLOAT); // mb
    cpl_table_new_column(overlap_region, "TEM", CPL_TYPE_FLOAT); // K
    if (exh2o == 1) {
        cpl_table_new_column(overlap_region, "H2O", CPL_TYPE_FLOAT); // ppmv
    }

    /* Set TEM, PRE, & H2O columns */
    for (i = 0; i < n_gdas_max; i++) {
        p = cpl_table_get_float(merged_profile, "PRE", i, NULL);
        cpl_table_set_float(overlap_region, "PRE", i, p);
        t = cpl_table_get_float(merged_profile, "TEM", i, NULL);
        cpl_table_set_float(overlap_region, "TEM", i, t);
        if (exh2o == 1) {
            ppmv = cpl_table_get_float(merged_profile, "H2O", i, NULL);
            cpl_table_set_float(overlap_region, "H2O", i, ppmv);
        }
    }

    /* Interpolate pressure column (ATM range) */
    err_code += mf_basic_interpolcolumn(atm_profile, "HGT", "PRE",
                                        merged_profile, "HGT", "PRE",
                                        n_gdas_max-1, merged_nrows);

    /* Interpolate temperature column (ATM range) */
    err_code += mf_basic_interpolcolumn(atm_profile, "HGT", "TEM",
                                        merged_profile, "HGT", "TEM",
                                        n_gdas_max-1, merged_nrows);

    /* Insert molecules (Insert H2O only above GDAS data, rest from 0km) */
    for (i = 0; i < n_molecs; i++) {
        strcpy(mol, cpl_array_get_string(molecs, i));
        if (strcmp(mol, "H2O") == 0) {
            err_code += mf_basic_interpolcolumn(atm_profile, "HGT", "H2O",
                                                merged_profile, "HGT", "H2O",
                                                n_gdas_max-1, merged_nrows);
        } else {
            err_code += mf_basic_interpolcolumn(atm_profile, "HGT", mol,
                                                merged_profile, "HGT", mol,
                                                0, merged_nrows);
        }
    }

    /* Free memory */
    free(mol);

    /* Insert overlap region */
    for (i = n_gdas_max-1; i < n_merge; i++) {
        hfrac = (h_merge - hgt_levels[i]) / dh;
        val_merged = cpl_table_get_float(merged_profile, "PRE", i, NULL);
        if (i == n_gdas_max-1) {
            val_overlap = cpl_table_get_float(overlap_region, "PRE", i, NULL);
            dp = val_overlap / val_merged - 1;
        }
        cpl_table_set_float(merged_profile, "PRE", i,
                            val_merged * (1 + dp * hfrac));
        val_merged = cpl_table_get_float(merged_profile, "TEM", i, NULL);
        if (i == n_gdas_max-1) {
            val_overlap = cpl_table_get_float(overlap_region, "TEM", i, NULL);
            dt = val_overlap / val_merged - 1;
        }
        cpl_table_set_float(merged_profile, "TEM", i,
                            val_merged * (1 + dt * hfrac));
        if (exh2o == 1) {
            val_merged = cpl_table_get_float(merged_profile, "H2O", i, NULL);
            if (i == n_gdas_max-1) {
                val_overlap = cpl_table_get_float(overlap_region, "H2O", i,
                                                  NULL);
                dh2o = val_overlap / val_merged - 1;
            }
            cpl_table_set_float(merged_profile, "H2O", i,
                                val_merged * (1 + dh2o * hfrac));
        }
    }

    /* Delete temporary table */
    cpl_table_delete(overlap_region);

    /* Check for errors */
    if (err_code != CPL_ERROR_NONE) {
        return cpl_error_set_message(cpl_func, MF_ERROR_SUBROUTINE,
                                     "Problem in interpolation");
    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_atm_adaptenv(cpl_table *profile, const mfdrv *drvpar)
{
    /*!
     * \brief
     *   Routine for running ::mf_atm_adaptenv_basic with parameters from the
     *   molecfit driver file
     *
     * This function extracts the parameters gelev (geoelevation in [m]),
     * gpres (geopressure in [mbar]), gtemp (geotemperature in [deg C]), ghum
     * (geohumidity in [ppmv]) from the molecfit driver file. Then it starts
     * ::mf_atm_adaptenv_basic with these parameters and converts the input
     * profile. If the upper mixing height \e emix is lower than the
     * geoelevation, the input atmospheric profile remains unchanged.
     *
     * \note
     *   The input profile gets overwritten on output
     *
     * \b INPUT:
     * \param profile  input atmospheric profile
     * \param drvpar   ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param profile  atmospheric profile adapted to local conditions
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    double gelev = 0., emix = 0., gpres = 0., gtemp = 0., temp = 0., hum = 0.;
    double ghum = 0.;

    /* Get geoelevation of observing site in km */
    p = cpl_parameterlist_find(drvpar->parlist, "geoelev");
    gelev = cpl_parameter_get_double(p) / 1000.;

    /* Get upper mixing height for consideration of meteo station
       parameters */
    p = cpl_parameterlist_find(drvpar->parlist, "emix");
    emix = cpl_parameter_get_double(p);

    /* If emix <= gelev, do not consider data of local meteo station */
    if (emix <= gelev) {
        cpl_msg_info(cpl_func, "Do not consider local meteo data");
        return CPL_ERROR_NONE;
    }

    /* Get data of local meteo station */
    p = cpl_parameterlist_find(drvpar->parlist, "pres");
    gpres = cpl_parameter_get_double(p);
    p = cpl_parameterlist_find(drvpar->parlist, "temp");
    gtemp = cpl_parameter_get_double(p);
    temp = gtemp + 273.15;
    p = cpl_parameterlist_find(drvpar->parlist, "rhum");
    hum = cpl_parameter_get_double(p);

    /* Convert relative humdity into ppmv */
    mf_basic_rhum2ppmv(&ghum, &temp, &gpres, &hum);

    /* Write info message */
    cpl_msg_info(cpl_func, "Consider local meteo data");

    /* Modify atmospheric profile */
    mf_atm_adaptenv_basic(profile, gelev, gpres, gtemp, ghum, emix);

    return CPL_ERROR_NONE;
}


void mf_atm_adaptenv_basic(cpl_table *profile, const double gelev,
                           const double gpres, const double gtemp,
                           const double ghum, const double emix)
{
    /*!
     * \brief
     *   Adapt temperature/pressure/humidity values to local conditions.
     *
     * This routine adapts the pressure, temperature, and humidity to local
     * on-sight measurements:
     * Up to \e gelev, the ratio of the on-site value and the profile value
     * at \e gelev is used as scaling factor. Profile values above \e emix are
     * left untouched. In between, a smooth transition is created by linear
     * interpolation.
     *
     * \b INPUT:
     * \param gelev    geoelevation in [km]
     * \param gpres    geopressure in [mbar]
     * \param gtemp    geotemperature in [deg C]
     * \param ghum     geohumidity in [ppmv]
     * \param emix     upper mixing height in [km]
     *
     * \b OUTPUT:
     * \param profile  atmospheric profile
     *
     * \b ERRORS:
     * - none
     */

    cpl_boolean exh2o = CPL_TRUE;
    int i = 0, nint = 0;
    double pres = 0., temp = 0., hum = 0., dpres = 0., dtemp = 0., dhum = 0.;
    double frac = 0., hgt = 0., p = 0., t = 0., ppmv = 0., rhum = 0.;

    /* Find starting height */
    for (i = 1; gelev >= cpl_table_get(profile, "HGT", i, NULL); i++);

    /* Check existence of H2O column */
    if (cpl_table_has_column(profile, "H2O") != 1) {
        exh2o = CPL_FALSE;
    }

    /* Interpolate pressure / temperature / humidity value at gelev */
    pres = mf_atm_interpollog(&gelev, profile, i, "PRE");
    temp = mf_atm_interpollog(&gelev, profile, i, "TEM");
    if (exh2o == CPL_TRUE) {
        hum = mf_atm_interpollog(&gelev, profile, i, "H2O");
        if (hum < MF_TOL) {
            hum = ghum;
            cpl_msg_warning(cpl_func, "Do not consider local rel. humidity "
                          "due to zero value in model profile at altitude of "
                          "meteo station");
        }
    } else {
        hum = ghum;
    }

    /* Get the number nint of height levels up to emix */
    for (nint = i; emix >= cpl_table_get(profile, "HGT", nint, NULL); nint++);

    /* Relative pressure / temperature / humidity shifts at gelev */
    dpres = gpres / pres - 1;
    dtemp = (gtemp + 273.15) / temp - 1;
    dhum = ghum / hum - 1;
    //printf("d: %g %g %g | %g\n", dpres, dtemp, dhum, gelev);

    /* Adapt pressure / temperature / humidity profiles to local on-site
       measurements */

    for (i = 0; i < nint; i++) {

        /* Calculate change of scaling factor depending on height */
        hgt = cpl_table_get(profile, "HGT", i, NULL);
        frac = (hgt - gelev) / (emix - gelev);
        if (frac > 1) {
            frac = 1;
        } else if (frac < 0) {
            frac = 0;
        }

        /* Read pressure / temperature / humidity */
        p = cpl_table_get(profile, "PRE", i, NULL);
        t = cpl_table_get(profile, "TEM", i, NULL);
        if (exh2o == CPL_TRUE) {
            ppmv = cpl_table_get(profile, "H2O", i, NULL);
        }

        /* Check for reasonable relative humidity (no correction so far) */
        if (exh2o == CPL_TRUE) {
            mf_basic_ppmv2rhum(&rhum, &t, &p, &ppmv);
            if (rhum * (1 + dhum * (1 - frac)) > 100.) {
                cpl_msg_info(cpl_func, "rel. humidity > 100%% for layer %d",
                             i);
            }
        }

        /* Modify pressure / temperature / humidity */
        cpl_table_set(profile, "PRE", i, p * (1 + dpres * (1 - frac)));
        if (hgt < gelev) {
            /* Constant values below gelev improve LBLRTM stability */
            cpl_table_set(profile, "TEM", i, gtemp + 273.15);
            if (exh2o == CPL_TRUE) {
                cpl_table_set(profile, "H2O", i, ghum);
            }
        } else {
            cpl_table_set(profile, "TEM", i, t * (1 + dtemp * (1 - frac)));
            if (exh2o == CPL_TRUE) {
                cpl_table_set(profile, "H2O", i,
                              ppmv * (1 + dhum * (1 - frac)));
            }
        }

        /*
        p = cpl_table_get(profile, "PRE", i, NULL);
        t = cpl_table_get(profile, "TEM", i, NULL);
        ppmv = cpl_table_get(profile, "H2O", i, NULL);
        mf_basic_ppmv2rhum(&rhum, &t, &p, &ppmv);
        printf("%g %g %g %g %g\n", hgt, p, t, ppmv, rhum);
        */

    }
}


double mf_atm_interpollog(const double *x, cpl_table *profile,
                          const int start, const char *varstr)
{
    /*!
     * \brief
     *   Logarithmic interpolation.
     *
     * This routine reads the two values supposedly surrounding \em x, as
     * defined by the starting index \em start and logarithmically
     * interpolates the value at \em x. The values, which are to be
     * interpolated are taken from \em profile.
     *
     * \note
     *   The x-values reference column is hard-coded to "HGT"
     *
     * \b INPUT:
     * \param x        x-value for interpolation
     * \param profile  reference table
     * \param start    index of lower value in \em profile surrounding \em x
     * \param varstr   column-id of y-values for interpolation
     *
     * \b OUTPUT:
     * \return y       interpolated y-value
     *
     * \b ERRORS:
     * - none
     */

    int j;
    double y = 0;
    double xref[2] = {0, 0}, yref[2] = {0, 0};

    /* interpolate pressure value at gelev */
    for (j = 0; j < 2; j++) {
        yref[j] = (double)cpl_table_get_float(profile, varstr, start-1+j,
                                              NULL);
        /* ensure that yref >0 */
        if (yref[j] <= 0.) {
            yref[j] = DBL_MIN * 2;
            cpl_msg_warning(cpl_func, "%s must be larger than 0", varstr);
        }
        yref[j] = log(yref[j]);
        xref[j] = (double)cpl_table_get_float(profile, "HGT", start-1+j,
                                              NULL);
    }
    mf_basic_interpollin(x, &y, 1, xref, yref, 2);

    /* avoid overflow */
    return y >= DBL_MAX_EXP ? DBL_MAX : exp(y);
}


cpl_error_code mf_atm_writeatm(cpl_table *atm_profile,
    cpl_table **atm_profile_out,
    const char *atm_file)
{
    /*!
     * \brief
     *   Write atmospheric profile in FASCOD format.
     *
     * This function writes an atmospheric standard profile into a
     * \c CPL_TABLE.
     *
     * \b INPUT:
     * \param atm_file     filename of an atmospheric profile.
     *
     * \b OUTPUT:
     * \param atm_profile  \c CPL_TABLE with standard profile.
     * \return             CPL_ERROR_NONE on success,
     *                     CPL_ERROR_FILE_NOT_FOUND on error.
     */

    cpl_array *col_ids;

    FILE *fp;                        /* file pointer */

    long nrows, ncols;

    int i, j, cr;

    char molec[MF_MAXLEN];

    mf_basic_initstring(molec, MF_MAXLEN);

    /* number of rows & cols in atm_profile */
    nrows = cpl_table_get_nrow(atm_profile);
    ncols = cpl_table_get_ncol(atm_profile);

    /* get names of columns in atm_profile */
    col_ids = cpl_table_get_column_names(atm_profile);

    fp = fopen(atm_file, "w");
    if (fp == NULL) {
        cpl_array_delete(col_ids);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                     "Could not open atm_file: %s", atm_file);
    } else {
        /* fileheader */
        fprintf(fp, "! FASCOD\n");
        fprintf(fp, "          %li  ! No.Levels in profiles\n", nrows);

        /* HEIGHT */
        fprintf(fp, "*%s [km]\n", cpl_array_get_string(col_ids, 0));

        for (i = 0, cr = 0; i < nrows; i++, cr++) {
            if (cr == 5) {
                fprintf(fp, "\n");
                cr = 0;
            }
            fprintf(fp, " %9.3f",
                    cpl_table_get_float(atm_profile, "HGT", i, NULL));
        }

        /* PRESSURE */
        fprintf(fp, "\n*%s [mbar]\n", cpl_array_get_string(col_ids, 1));
        for (i = 0, cr = 0; i < nrows; i++, cr++) {
            if (cr == 5) {
                fprintf(fp, "\n");
                cr = 0;
            }
            if (i == 0) {
                /* Trick to avoid error by radiative transfer code */
                fprintf(fp, " %9.3E",
                        cpl_table_get_float(atm_profile, "PRE", i, NULL)+0.1);
            } else {
                fprintf(fp, " %9.3E",
                        cpl_table_get_float(atm_profile, "PRE", i, NULL));
            }
        }

        /* TEMPERATURE */
        fprintf(fp, "\n*%s [K]\n", cpl_array_get_string(col_ids, 2));
        for (i = 0, cr = 0; i < nrows; i++, cr++) {
            if (cr == 5) {
                fprintf(fp, "\n");
                cr = 0;
            }
            fprintf(fp, " %9.2f",
                    cpl_table_get_float(atm_profile, "TEM", i, NULL));
        }

        /* MOLECULES */
        for (i = 3; i < ncols; i++) {
            strcpy(molec, cpl_array_get_string(col_ids, i));
            fprintf(fp, "\n*%s [ppmv]\n", molec);
            for (j = 0, cr = 0; j < nrows; j++, cr++) {
                if (cr == 5) {
                    fprintf(fp, "\n");
                    cr = 0;
                }
                fprintf(fp, " %9.3E",
                        cpl_table_get_float(atm_profile, molec, j, NULL));
            }
        }

        fprintf(fp, "\n*END\n");
        fclose(fp);

        // Also save the table as a fits file
        char* fitsfilename = cpl_sprintf("%s.fits", atm_file);
        // Create output table
        // If not the first call, rewrite it.
        // Otherwise, we end up with a memory leak.
        if (*atm_profile_out) {
            cpl_table_delete(*atm_profile_out);
        }
        *atm_profile_out = cpl_table_duplicate(atm_profile);
        // Hack for the first value of the pressure
        cpl_table_set_float(*atm_profile_out, "PRE", 0, cpl_table_get_float(atm_profile, "PRE", 0, NULL) + 0.1);
        cpl_error_code status = cpl_table_save(*atm_profile_out, NULL, NULL, fitsfilename, CPL_IO_CREATE);

        cpl_free(fitsfilename);
        if (status) {
            cpl_array_delete(col_ids);
            return cpl_error_set_message(cpl_func, MF_ERROR_FOPEN,
                                     "Could not save atm_file: %s.fits", atm_file);
        }

        cpl_array_delete(col_ids);

        return CPL_ERROR_NONE;
    }
}


cpl_error_code mf_atm_scaletopwv(cpl_table *profile, const mfdrv *drvpar)
{
    /*!
     * Scales the water vapour profile of an atmospheric profile table to a
     * PWV value in mm given by the parameter \e pwv. If \e pwv is set to -1,
     * no scaling is performed.
     *
     * \b INPUT:
     * \param profile  input atmospheric profile
     * \param drvpar   ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param profile  atmospheric profile with scaled water vapour if
     *                 requested
     *
     * \b ERRORS:
     * - none
     */

    cpl_parameter *p;
    double outpwv = 0., gelev = 0., inpwv = 0.;

    /* Get output PWV */
    p = cpl_parameterlist_find(drvpar->parlist, "pwv");
    outpwv = cpl_parameter_get_double(p);

    /* Return if scaling is not requested (PWV <= 0) */
    if (outpwv <= 0) {
        return CPL_ERROR_NONE;
    }

    /* Get altitude of observer */
    p = cpl_parameterlist_find(drvpar->parlist, "geoelev");
    gelev = cpl_parameter_get_double(p) / 1000.;

    /* Calculate PWV for input profile */
    mf_atm_calpwv(&inpwv, profile, &gelev);

    /* Return if no H2O column exists (PWV = 0) */
    if (inpwv == 0.) {
        return CPL_ERROR_NONE;
    }

    /* Scale H2O profile */
    cpl_table_multiply_scalar(profile, "H2O", outpwv / inpwv);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_atm_calpwv(double *pwv, cpl_table *prof,
                             const double *geoelev)
{
    /*!
     * Calculates water vapour column pwv in [mm] from profile in [ppmv].
     * The input profiles are interpreted as functions.
     * The starting height is the altitude of the observing site.
     * The profile information is read from the output profiles file.
     *
     * \b INPUT:
     * \param prof     atmospheric profile
     * \param geoelev  geoelevation of observing site in [km]
     *
     * \b OUTPUT:
     * \param pwv      water vapour column in [mm]
     *
     * \b ERRORS:
     * - none
     */

    int nlayer = 0, i = 0;
    double uhgt = 0., lhgt = 0., dlayer = 0., pressure = 0., temp = 0.;
    double molcol = 0., ch = 0., nmol = 0.;
    const double km2m = 1e3, mbar2Pa = 1e2;
    const double molmass = 0.0182; /* mol mass of H2O in kg */

    /* Default water vapour column */
    *pwv = 0.;

    /* Check for H2O column in profile table and return if not present */
    if (cpl_table_has_column(prof, "H2O") != 1) {
        return CPL_ERROR_NONE;
    }

    /* Number of layers, initialisation of upper layer height, and total
       column height */
    nlayer = cpl_table_get_nrow(prof);
    uhgt = cpl_table_get(prof, "HGT", 0, NULL);

    /* Sum up water vapour columns of all layers */

    for (i = 0; i < nlayer-1; i++) {

        /* Lower and upper limit of layer */
        lhgt = uhgt;
        uhgt = cpl_table_get(prof, "HGT", i+1, NULL);

        /* Skip layers below height of observing site */
        if (uhgt <= *geoelev) continue;

        /* Thickness of layer in m */
        if (*pwv == 0. && uhgt > *geoelev) {
            dlayer = (uhgt - *geoelev) * km2m;
        } else {
            dlayer = (uhgt - lhgt) * km2m;
        }

        /* Average pressure and temperature for layer */
        pressure = (cpl_table_get(prof, "PRE", i, NULL)
                    + cpl_table_get(prof, "PRE", i+1, NULL)) / 2 ;
        temp = (cpl_table_get(prof, "TEM", i, NULL)
                + cpl_table_get(prof, "TEM", i+1, NULL)) / 2;
        if (temp <= 0) {
            temp = MF_TOL;
        }

        /* Average ppmv of H2O for layer */
        molcol = (cpl_table_get(prof, "H2O", i, NULL)
                  + cpl_table_get(prof, "H2O", i+1, NULL)) / 2;

        /* H2O column in mm */
        /* Column height [m] for layer */
        ch = 1e-6 * molcol * dlayer;
        /* Number of mols per unit area [mol m^-2] for layer */
        nmol = ch * pressure * mbar2Pa / (MF_R * temp);
        /* Mass per unit area [kg m^-2] */
        *pwv += nmol * molmass;
        /* Same value for column in mm since density of water is
           10^3 kg/m^3 and m to mm is 10^3 as well */

    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_atm_calpwv_histo(double *pwv, cpl_table *prof,
                                   const double *geoelev)
{
    /*!
     * Calculates water vapour column pwv in [mm] from profile in [ppmv].
     * The input profiles are interpreted as histograms.
     * The starting height is the altitude of the observing site.
     * The profile information is read from the output profiles file.
     *
     * \b INPUT:
     * \param prof     atmospheric profile
     * \param geoelev  geoelevation of observing site in [km]
     *
     * \b OUTPUT:
     * \param pwv      water vapour column in [mm]
     *
     * \b ERRORS:
     * - none
     */

    int nlayer = 0, i = 0;
    double uhgt = 0., hgtmax = 0., lhgt = 0., dlayer = 0., pressure = 0.;
    double temp = 0., molcol = 0., ch = 0., nmol = 0.;
    const double km2m = 1e3, mbar2Pa = 1e2;
    const double molmass = 0.0182; /* mol mass of H2O in kg */

    /* Default water vapour column */
    *pwv = 0.;

    /* Check for H2O column in profile table and return if not present */
    if (cpl_table_has_column(prof, "H2O") != 1) {
        return CPL_ERROR_NONE;
    }

    /* Number of layers, initialisation of upper layer height, and total
       column height */
    nlayer = cpl_table_get_nrow(prof);
    uhgt = 1.5 * cpl_table_get(prof, "HGT", 0, NULL) -
           0.5 * cpl_table_get(prof, "HGT", 1, NULL);
    hgtmax = 1.5 * cpl_table_get(prof, "HGT", nlayer-1, NULL) -
             0.5 * cpl_table_get(prof, "HGT", nlayer-2, NULL);

    /* Sum up water vapour columns of all layers */

    for (i = 0; i < nlayer-1; i++) {

        /* Lower and upper limit of layer */
        lhgt = uhgt;
        if (i == nlayer - 1) {
            uhgt = hgtmax;
        } else {
            uhgt = (cpl_table_get(prof, "HGT", i, NULL) +
                    cpl_table_get(prof, "HGT", i+1, NULL)) / 2;
        }

        /* Skip layers below height of observing site */
        if (uhgt <= *geoelev) continue;

        /* Thickness of layer in m */
        if (*pwv == 0. && uhgt > *geoelev) {
            dlayer = (uhgt - *geoelev) * km2m;
        } else {
            dlayer = (uhgt - lhgt) * km2m;
        }

        /* Average pressure and temperature for layer */
        pressure = cpl_table_get(prof, "PRE", i, NULL);
        temp = cpl_table_get(prof, "TEM", i, NULL);
        if (temp <= 0) {
            temp = MF_TOL;
        }

        /* ppmv of H2O for layer height */
        molcol = cpl_table_get(prof, "H2O", i, NULL);

        /* H2O column in mm */
        /* Column height [m] for layer */
        ch = 1e-6 * molcol * dlayer;
        /* Number of mols per unit area [mol m^-2] for layer */
        nmol = ch * pressure * mbar2Pa / (MF_R * temp);
        /* Mass per unit area [kg m^-2] */
        *pwv += nmol * molmass;
        /* Same value for column in mm since density of water is
           10^3 kg/m^3 and m to mm is 10^3 as well */

    }

    return CPL_ERROR_NONE;
}


cpl_error_code mf_atm_calcol(double *h2ocol, mfdrv *drvpar)
{
    /*!
     * Calculates water vapour column in mm from profile in ppmv.
     * The input profiles are interpreted as functions.
     * The starting height is the altitude of the observing site.
     * The profile information is read from the output profiles file.
     * Moreover, the ppmv of the entire atmospheric column is provided for
     * all molecules provided by the ::mfdrv parameter structure.
     * This information is written into a special ppmv column of the table of
     * molecules of ::mfdrv.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param h2ocol  water vapour column in mm
     * \param drvpar  ::mfdrv parameter structure supplemented by ppmv of
     *                provided molecules
     *
     * \b ERRORS:
     * - File opening failed
     * - Invalid object structure
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_table *prof;
    char basedir[MF_MAXLEN], output_dir[MF_MAXLEN];
    char output_name[MF_MAXLEN], atmfile[MF_MAXLEN];
    char errtxt[MF_MAXLEN];
    char **mol;
    int nmolec = 0, j = 0, nlayer = 0, i = 0;
    double *ppmv, geoelev = 0., uhgt = 0., nmol0 = 0., lhgt = 0., dlayer = 0.;
    double pressure = 0., temp = 0., molcol = 0., ch = 0., nmol = 0.;
    const double km2m = 1e3, mbar2Pa = 1e2;
    const double molmass = 0.0182; /* mol mass of H2O in kg */

    /* Default water vapour column */
    *h2ocol = 0.;

    /* Read profiles from file */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(output_dir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strcpy(output_name, cpl_parameter_get_string(p));
    // sprintf(atmfile, "%s%s_fit.atm", output_dir, output_name);
    // prof = cpl_table_new(0);
    // status = mf_atm_readatm(prof, atmfile);
    sprintf(atmfile, "%s%s_fit.atm.fits", output_dir, output_name);
    status = mf_atm_readatm_fromFits(&prof, atmfile);
    if (status != CPL_ERROR_NONE) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, atmfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Create column for ppmv of molecules in mfdrv structure
       if not present */
    if (cpl_table_has_column(drvpar->molectab, "ppmv") != 1) {
        cpl_table_new_column(drvpar->molectab, "ppmv", CPL_TYPE_DOUBLE);
        nmolec = cpl_table_get_nrow(drvpar->molectab);
    }

    /* Set initial ppmv values to 0 */
    nmolec = cpl_table_get_nrow(drvpar->molectab);
    cpl_table_fill_column_window(drvpar->molectab, "ppmv", 0, nmolec, 0.);

    /* Get pointer to names of molecules in driver parameter structure */
    mol = cpl_table_get_data_string(drvpar->molectab, "list_molec");

    /* Get pointer to ppmv of molecules in driver parameter structure */
    ppmv = cpl_table_get_data_double(drvpar->molectab, "ppmv");

    /* Check existence of selected molecular columns in profile table */
    for (j = 0; j < nmolec; j++) {
        if (cpl_table_has_column(prof, mol[j]) != 1) {
            sprintf(errtxt, "%s: cpl_table *prof (no %s column)",
                    MF_ERROR_IOS_TXT, mol[j]);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }
    }

    /* Altitude of observing site */
    p = cpl_parameterlist_find(drvpar->parlist, "geoelev");
    geoelev = cpl_parameter_get_double(p) / km2m;

    /* Number of layers, initialisation of upper layer height, and total
       column height */
    nlayer = cpl_table_get_nrow(prof);
    uhgt = cpl_table_get(prof, "HGT", 0, NULL);

    /* Sum up molecular columns of all layers */

    for (nmol0 = 0., i = 0; i < nlayer-1; i++) {

        /* Lower and upper limit of layer */
        lhgt = uhgt;
        uhgt = cpl_table_get(prof, "HGT", i+1, NULL);

        /* Skip layers below height of observing site */
        if (uhgt <= geoelev) continue;

        /* Thickness of layer in m */
        if (*h2ocol == 0. && uhgt > geoelev) {
            dlayer = (uhgt - geoelev) * km2m;
        } else {
            dlayer = (uhgt - lhgt) * km2m ;
        }

        /* Average pressure and temperature for layer */
        pressure = (cpl_table_get(prof, "PRE", i, NULL)
                    + cpl_table_get(prof, "PRE", i+1, NULL)) / 2 ;
        temp = (cpl_table_get(prof, "TEM", i, NULL)
                + cpl_table_get(prof, "TEM", i+1, NULL)) / 2;
        if (temp <= 0.) {
            temp = MF_TOL;
        }

        /* Number of mols per unit area for each molecule and PWV */
        for (j = 0; j < nmolec; j++) {
            /* Average ppmv of molecule for layer */
            molcol = (cpl_table_get(prof, mol[j], i, NULL)
                      + cpl_table_get(prof, mol[j], i+1, NULL)) / 2;
            /* Column height [m] of molecule for layer */
            ch = 1e-6 * molcol * dlayer;
            /* Number of mols per unit area [mol m^-2] for layer */
            nmol = ch * pressure * mbar2Pa / (MF_R * temp);
            /* Number of mols per unit area for atmosphere */
            ppmv[j] += nmol;
            /* Number of mols per unit area for air */
            if (j == 0) {
                nmol0 += dlayer * pressure * mbar2Pa / (MF_R * temp);
            }
            /* H2O column in mm (PWV) */
            if (strncmp(mol[j], "H2O", 3) == 0) {
                /* Mass per unit area [kg m^-2] */
                *h2ocol += nmol * molmass;
                /* Same value for column in mm since density of water is
                   10^3 kg/m^3 and m to mm is 10^3 as well */
            }
        }

    }

    /* Volume mixing ratio for each molecule */
    for (j = 0; j < nmolec; j++) {
        if (nmol0 <= 0.) {
            ppmv[j] = 0.;
        } else {
            ppmv[j] *= 1e6 / nmol0;
        }
    }

    /* Free allocated memory */
    cpl_table_delete(prof);

    return CPL_ERROR_NONE;
}


cpl_error_code mf_atm_calcol_histo(double *h2ocol, mfdrv *drvpar)
{
    /*!
     * Calculates water vapour column in mm from profile in ppmv.
     * The input profiles are interpreted as histograms.
     * The starting height is the altitude of the observing site.
     * The profile information is read from the output profiles file.
     * Moreover, the ppmv of the entire atmospheric column is provided for
     * all molecules provided by the ::mfdrv parameter structure.
     * This information is written into a special ppmv column of the table of
     * molecules of ::mfdrv.
     *
     * \b INPUT:
     * \param drvpar  ::mfdrv parameter structure
     *
     * \b OUTPUT:
     * \param h2ocol  water vapour column in mm
     * \param drvpar  ::mfdrv parameter structure supplemented by ppmv of
     *                provided molecules
     *
     * \b ERRORS:
     * - File opening failed
     * - Invalid object structure
     */

    cpl_error_code status = CPL_ERROR_NONE;
    cpl_parameter *p;
    cpl_table *prof;
    char basedir[MF_MAXLEN], output_dir[MF_MAXLEN];
    char output_name[MF_MAXLEN], atmfile[MF_MAXLEN];
    char errtxt[MF_MAXLEN];
    char **mol;
    int nmolec = 0, j = 0, nlayer = 0, i = 0;
    double *ppmv, geoelev = 0., uhgt = 0., hgtmax = 0., nmol0 = 0., lhgt = 0.;
    double dlayer = 0., pressure = 0., temp = 0., molcol = 0., ch = 0.;
    double nmol = 0.;
    const double km2m = 1e3, mbar2Pa = 1e2;
    const double molmass = 0.0182; /* mol mass of H2O in kg */

    /* Default water vapour column */
    *h2ocol = 0.;

    /* Read profiles from file */
    p = cpl_parameterlist_find(drvpar->parlist, "basedir");
    strcpy(basedir, cpl_parameter_get_string(p));
    char curdir[MF_MAXLEN];
    p = cpl_parameterlist_find(drvpar->parlist, "curdir");
    strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
    p = cpl_parameterlist_find(drvpar->parlist, "output_dir");
    mf_basic_abspath(output_dir, cpl_parameter_get_string(p), curdir);
    p = cpl_parameterlist_find(drvpar->parlist, "output_name");
    strcpy(output_name, cpl_parameter_get_string(p));
    // sprintf(atmfile, "%s%s_fit.atm", output_dir, output_name);
    // prof = cpl_table_new(0);
    // status = mf_atm_readatm(prof, atmfile);
    sprintf(atmfile, "%s%s_fit.atm.fits", output_dir, output_name);
    status = mf_atm_readatm_fromFits(&prof, atmfile);
    if (status != CPL_ERROR_NONE) {
        sprintf(errtxt, "%s: %s", MF_ERROR_FOF_TXT, atmfile);
        return cpl_error_set_message(cpl_func, MF_ERROR_FOF, "%s", errtxt);
    }

    /* Create column for ppmv of molecules in mfdrv structure
       if not present */
    if (cpl_table_has_column(drvpar->molectab, "ppmv") != 1) {
        cpl_table_new_column(drvpar->molectab, "ppmv", CPL_TYPE_DOUBLE);
        nmolec = cpl_table_get_nrow(drvpar->molectab);
    }

    /* Set initial ppmv values to 0 */
    nmolec = cpl_table_get_nrow(drvpar->molectab);
    cpl_table_fill_column_window(drvpar->molectab, "ppmv", 0, nmolec, 0.);

    /* Get pointer to names of molecules in driver parameter structure */
    mol = cpl_table_get_data_string(drvpar->molectab, "list_molec");

    /* Get pointer to ppmv of molecules in driver parameter structure */
    ppmv = cpl_table_get_data_double(drvpar->molectab, "ppmv");

    /* Check existence of selected molecular columns in profile table */
    for (j = 0; j < nmolec; j++) {
        if (cpl_table_has_column(prof, mol[j]) != 1) {
            sprintf(errtxt, "%s: cpl_table *prof (no %s column)",
                    MF_ERROR_IOS_TXT, mol[j]);
            return cpl_error_set_message(cpl_func, MF_ERROR_IOS, "%s",
                                         errtxt);
        }
    }

    /* Altitude of observing site */
    p = cpl_parameterlist_find(drvpar->parlist, "geoelev");
    geoelev = cpl_parameter_get_double(p) / km2m;

    /* Number of layers, initialisation of upper layer height, and total
       column height */
    nlayer = cpl_table_get_nrow(prof);
    uhgt = 1.5 * cpl_table_get(prof, "HGT", 0, NULL) -
           0.5 * cpl_table_get(prof, "HGT", 1, NULL);
    hgtmax = 1.5 * cpl_table_get(prof, "HGT", nlayer-1, NULL) -
             0.5 * cpl_table_get(prof, "HGT", nlayer-2, NULL);

    /* Sum up molecular columns of all layers */

    for (nmol0 = 0., i = 0; i < nlayer; i++) {

        /* Lower and upper limit of layer */
        lhgt = uhgt;
        if (i == nlayer - 1) {
            uhgt = hgtmax;
        } else {
            uhgt = (cpl_table_get(prof, "HGT", i, NULL) +
                    cpl_table_get(prof, "HGT", i+1, NULL)) / 2;
        }

        /* Skip layers below height of observing site */
        if (uhgt <= geoelev) continue;

        /* Thickness of layer in m */
        if (*h2ocol == 0. && uhgt > geoelev) {
            dlayer = (uhgt - geoelev) * km2m;
        } else {
            dlayer = (uhgt - lhgt) * km2m;
        }

        /* Average pressure and temperature for layer */
        pressure = cpl_table_get(prof, "PRE", i, NULL);
        temp = cpl_table_get(prof, "TEM", i, NULL);
        if (temp <= 0.) {
            temp = MF_TOL;
        }

        /* Number of mols per unit area for each molecule and PWV */
        for (j = 0; j < nmolec; j++) {
            /* ppmv of molecule for layer height */
            molcol = cpl_table_get(prof, mol[j], i, NULL);
            /* Column height [m] of molecule for layer */
            ch = 1e-6 * molcol * dlayer;
            /* Number of mols per unit area [mol m^-2] for layer */
            nmol = ch * pressure * mbar2Pa / (MF_R * temp);
            /* Number of mols per unit area for atmosphere */
            ppmv[j] += nmol;
            /* Number of mols per unit area for air */
            if (j == 0) {
                nmol0 += dlayer * pressure * mbar2Pa / (MF_R * temp);
            }
            /* H2O column in mm (PWV) */
            if (strncmp(mol[j], "H2O", 3) == 0) {
                /* Mass per unit area [kg m^-2] */
                *h2ocol += nmol * molmass;
                /* Same value for column in mm since density of water is
                   10^3 kg/m^3 and m to mm is 10^3 as well */
            }
        }

    }

    /* Volume mixing ratio for each molecule */
    for (j = 0; j < nmolec; j++) {
        if (nmol0 <= 0.) {
            ppmv[j] = 0.;
        } else {
            ppmv[j] *= 1e6 / nmol0;
        }
    }

    /* Free allocated memory */
    cpl_table_delete(prof);

    return CPL_ERROR_NONE;
}

/**@}*/
