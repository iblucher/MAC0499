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
 * \defgroup molecfit Atmospheric Molecular Fitting
 *
 * This module provides functions for MOLECFIT.
 */

/**@{*/

/*!
 * \file mf_basic.h
 *
 * Basic MOLECFIT header
 *
 * \author Stefan Noll & ESO In-Kind Team Innsbruck
 * \since  27 Apr 2010
 * \date   03 Nov 2014
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

/* C standard libraries */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <signal.h>
#include <unistd.h>
#include <sys/dir.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

/* CPL library */

#include <cpl.h>

/*****************************************************************************
 *                                DEFINES                                    *
 ****************************************************************************/

#ifndef MF_BASIC_H
#define MF_BASIC_H

/* Basic definitions */

/*! Enumeration structure for identifying integer and double type */
#ifndef MF_TYPE
#define MF_TYPE
typedef enum MF_TYPE {t_int, t_double} mftype;
#endif

/*! Macro for finding the minimum of two values */
#define MF_MIN(a, b) ((a) > (b) ? (b) : (a))
/*! Macro for finding the maximum of two values */
#define MF_MAX(a, b) ((a) < (b) ? (b) : (a))

/* Definition of constants */

/*! Maximum number of characters per line */
#define MF_LENLINE 510
/*! Maximum number of string characters */
#define MF_MAXLEN 1024
/*! Required relative accuracy for data comparisons */
#define MF_TOL 1e-7
/*! Maximum number of space- or tab-separated strings per line
    (see ::mf_basic_readline) */
#define MF_MAXPAR 100

/* Definition of error codes and corresponding standard error messages */

/*! Enumeration structure for sky model related error codes */

typedef enum _mf_error_code_ {
    MF_ERROR_FOPEN = CPL_ERROR_EOL + 11, // "File opening failed"
    MF_ERROR_UFS = CPL_ERROR_EOL + 12, // "Unexpected file structure"
    MF_ERROR_IFE = CPL_ERROR_EOL + 13, // "Invalid file name extension"
    MF_ERROR_NDA = CPL_ERROR_EOL + 20, // "No data"
    MF_ERROR_INSUFF_DATA = CPL_ERROR_EOL + 21, // "Insufficient data points"
    MF_ERROR_IDG = CPL_ERROR_EOL + 22, // "Inconsistent data grids"
    MF_ERROR_IDR = CPL_ERROR_EOL + 23, // "Invalid data range"
    MF_ERROR_IOD = CPL_ERROR_EOL + 24, // "Invalid order of data points"
    MF_ERROR_IIP = CPL_ERROR_EOL + 30, // "Invalid input parameter(s)"
    MF_ERROR_IOV = CPL_ERROR_EOL + 31, // "Invalid object value(s)"
    MF_ERROR_IOS = CPL_ERROR_EOL + 32, // "Invalid object structure"
    MF_ERROR_SUBROUTINE = CPL_ERROR_EOL + 40, // "Error in subroutine"
    // ERRORS RELATED TO access(), mkdir(), chdir()
    MF_ERROR_ACCES       = CPL_ERROR_EOL + 50, // access(), mkdir(), chdir()
    MF_ERROR_LOOP        = CPL_ERROR_EOL + 51, // access(), mkdir(), chdir()
    MF_ERROR_NAMETOOLONG = CPL_ERROR_EOL + 52, // access(), mkdir(), chdir()
    MF_ERROR_NOENT       = CPL_ERROR_EOL + 53, // access(), mkdir(), chdir()
    MF_ERROR_NOTDIR      = CPL_ERROR_EOL + 54, // access(), mkdir(), chdir()
    MF_ERROR_ROFS        = CPL_ERROR_EOL + 55, // access(), mkdir()
    MF_ERROR_FAULT       = CPL_ERROR_EOL + 56, // access(), mkdir(), chdir()
    MF_ERROR_INVAL       = CPL_ERROR_EOL + 57, // access()
    MF_ERROR_IO          = CPL_ERROR_EOL + 58, // access(),          chdir()
    MF_ERROR_NOMEM       = CPL_ERROR_EOL + 59, // access(), mkdir(), chdir()
    MF_ERROR_TXTBSY      = CPL_ERROR_EOL + 60, // access()
    MF_ERROR_EXIST       = CPL_ERROR_EOL + 61, //           mkdir()
    MF_ERROR_NOSPC       = CPL_ERROR_EOL + 62, //           mkdir()
    MF_ERROR_PERM        = CPL_ERROR_EOL + 63, //           mkdir()
    // General errors
    MF_ERROR_BADUSERINPUT = CPL_ERROR_EOL + 70,
    MF_ERROR_LINK        = CPL_ERROR_EOL + 71,
    MF_ERROR_CD          = CPL_ERROR_EOL + 72, // could not change directory
    MF_ERROR_GETCWD      = CPL_ERROR_EOL + 73, // could not get current
                                               // working directory
    MF_ERROR_UNDEF = CPL_ERROR_EOL + 80
} mf_error_code;

/* Aliases for sky model related error codes */

#define MF_ERROR_FOF MF_ERROR_FOPEN
#define MF_ERROR_ISM MF_ERROR_NOMEM
#define MF_ERROR_EIS MF_ERROR_SUBROUTINE
#define MF_ERROR_ISD MF_ERROR_INSUFF_DATA

/* Standard messages for sky model related errors */

#define MF_ERROR_FOPEN_TXT "File opening failed"
#define MF_ERROR_UFS_TXT "Unexpected file structure"
#define MF_ERROR_IFE_TXT "Invalid file name extension"
#define MF_ERROR_BDR_TXT "Bad directory"
#define MF_ERROR_NDA_TXT "No data"
#define MF_ERROR_ISD_TXT "Insufficient data points"
#define MF_ERROR_IDG_TXT "Inconsistent data grids"
#define MF_ERROR_IDR_TXT "Invalid data range"
#define MF_ERROR_IOD_TXT "Invalid order of data points"
#define MF_ERROR_IIP_TXT "Invalid input parameter(s)"
#define MF_ERROR_IOV_TXT "Invalid object value(s)"
#define MF_ERROR_IOS_TXT "Invalid object structure"
#define MF_ERROR_SUBROUTINE_TXT "Error in subroutine"
// ERRORS RELATED TO access(), mkdir(), chdir()
#define MF_ERROR_ACCES_TXT       "Permission denied"
#define MF_ERROR_LOOP_TXT        "Too many symbolic links"
#define MF_ERROR_NAMETOOLONG_TXT "Pathname too long"
#define MF_ERROR_NOENT_TXT       "File/dir does not exist"
#define MF_ERROR_NOTDIR_TXT      "Component used as directory in pathname " \
                                 "is not a directory"
#define MF_ERROR_ROFS_TXT        "Write permission requested for file/dir " \
                                 "on read-only file system"
#define MF_ERROR_FAULT_TXT       "Pathname points outside accessible " \
                                 "address space"
#define MF_ERROR_INVAL_TXT       "Mode was incorrectly specified"
#define MF_ERROR_IO_TXT          "I/O error occurred"
#define MF_ERROR_NOMEM_TXT       "Insufficient memory"
#define MF_ERROR_TXTBSY_TXT      "Write access requested to executable " \
                                 "which is being executed"
#define MF_ERROR_EXIST_TXT       "File/dir already exists"
#define MF_ERROR_NOSPC_TXT       "No space left on device"
#define MF_ERROR_PERM_TXT        "File system does not support creation of " \
                                 "directories"
#define MF_ERROR_LINK_TXT        "Could not create symbolic link"
#define MF_ERROR_UNDEF_TXT       "Undefined error"
#define MF_ERROR_CD_TXT          "Could not change directory"
#define MF_ERROR_GETCWD_TXT      "Could not get current working directory"

#define MF_ERROR_FOF_TXT MF_ERROR_FOPEN_TXT
#define MF_ERROR_ISM_TXT MF_ERROR_NOMEM_TXT
#define MF_ERROR_EIS_TXT MF_ERROR_SUBROUTINE_TXT

/*****************************************************************************
 *                                TYPEDEF                                    *
 ****************************************************************************/

/* Definition of structures */

/*!
 * Structure for a parameter value in different data types
 *
 * \param c  string (maximum length ::MF_LENLINE)
 * \param i  integer number
 * \param d  float number of double precision
 */

typedef struct _mfpar_ {
    char c[MF_LENLINE+2];
    int i;
    double d;
} mfpar;

/*****************************************************************************
 *                               PROTOTYPES                                  *
 ****************************************************************************/

/* Declaration of functions */

cpl_error_code mf_basic_readline(FILE *stream, mfpar par[], int *npar);
cpl_error_code mf_basic_joinpar(char *str, mfpar par[], int i1, int i2);
double mf_basic_mjd2fracyear(double mjd);
double mf_basic_fracyear2date(int *year, int *month, int *day,
                              int *hh, int *mm, double *ss,
                              const double *fracyear);
cpl_error_code mf_basic_rebin(cpl_table *outspec, const char *outlam,
                              const char *outflux, const cpl_table *inspec,
                              const char *inlam, const char *influx);
cpl_error_code mf_basic_convolve(cpl_table *spec, const char *colname,
                                 const cpl_array *kernel);
cpl_error_code mf_basic_convolvewindow(cpl_array *convflux,
                                       const cpl_array *flux,
                                       const int range[2],
                                       const cpl_array *kernel);
cpl_error_code mf_basic_convolvewindow_inv(cpl_array *convflux,
                                           const cpl_array *flux,
                                           const int range[2],
                                           const cpl_array *kernel,
                                           const int kerncen);
cpl_error_code mf_basic_convolvekernels(cpl_array *outkernel,
                                        const cpl_array *inkernel1,
                                        const cpl_array *inkernel2);
cpl_error_code mf_basic_linecount(long *n_lines, FILE *fp);
void mf_basic_ppmv2rhum(double *hum, const double *tem, const double *p,
                        const double *ppmv);
void mf_basic_rhum2ppmv(double *ppmv, const double *tem, const double *p,
                        const double *hum);
cpl_error_code mf_basic_rhum2ppmv_old(const cpl_array *temp,
                                      const cpl_array *pres,
                                      const cpl_array *rhum,
                                      cpl_array *ppmv);
cpl_error_code mf_basic_planck(cpl_array *bb, const cpl_array *wavelength,
                               const double temp);
void mf_basic_dirslash(char *dir);
void mf_basic_abspath(char *out, const char *dir, const char *cwd);
cpl_error_code mf_basic_access(const char *pathname, const int mode);
cpl_error_code mf_basic_greg2jd(long *jd, const int year, const int month,
                                const int day);
cpl_error_code mf_basic_jd2greg(int *year, int *month, int *day,
                                const long jd);
cpl_error_code mf_basic_absfile(char *absfilename, const char *filename);
cpl_boolean mf_get_tempdir(char * tmpdir_);
char * mf_get_cwd(void);
cpl_boolean mf_basic_parameterlists_same(cpl_parameterlist *list1,
                                         cpl_parameterlist *list2);
void mf_basic_initstring(char *str, const long n);
cpl_boolean mf_basic_isnumber(char *str);
void mf_basic_lowercase(char *inputstr);
char *mf_basic_replacestring(char *instring, char *oldsubstr,
                             char *newsubstr);
char *mf_basic_rmcntrl(char *str);
void mf_basic_rmcntrl_inplace(char *str);
char *mf_basic_strtrim(char *str);
void mf_basic_strtrim_inplace(char *str);
void mf_basic_terminatestring(char *str);
cpl_error_code mf_basic_interpollin(const double *x_out, double *y_out,
                                    const long n_out, const double *x_ref,
                                    const double *y_ref, const long n_ref);
cpl_error_code mf_basic_interpolcolumn(const cpl_table *in_tab,
                                       const char *in_str_x,
                                       const char *in_str_y,
                                       cpl_table *out_tab,
                                       const char *out_str_x,
                                       const char *out_str_y,
                                       const int start,
                                       const int end);
void mf_basic_fillvector(double *v, const int nrows, const cpl_table *tab,
                         const char *str);
cpl_boolean mf_basic_syntaxok(char *line);
void mf_basic_processvariable(const char *str, const char *var_str,
                              char *line, void *var, int *success,
                              const mftype v_type);
void mf_basic_processrange(const char *str, const char *var_str, char *line,
                           void *var, int *success, const mftype v_type,
                           const double lo, const double hi);
cpl_error_code mf_basic_getfilename(char *dir, char *filename, char *suffix,
                                    const char *path);
cpl_error_code mf_basic_copytable(cpl_table *outtab, const cpl_table *intab);
cpl_error_code mf_basic_copyvector(cpl_vector *outvec,
                                   const cpl_vector *invec);
cpl_error_code mf_basic_getmaskval_vector(double maskval[2],
                                          const cpl_vector *vector);
cpl_error_code mf_basic_getmaskval_image(double maskval[2],
                                         const cpl_image *image);
cpl_error_code mf_basic_col2arr(cpl_array *arr, const cpl_table *tab,
                                const char *colname);
cpl_error_code mf_basic_arr2col(cpl_table *tab, const char *colname,
                                const cpl_array *arr);

#endif /* MF_BASIC_H */

#ifdef __cplusplus
}
#endif

/**@}*/
