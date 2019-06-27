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
 * \callgraph
 *
 * \file test_mf_plot.c
 *
 * test routine for plotting library
 *
 * \author Wolfgang Kausch & ESO In-Kind Team Innsbruck
 *
 */

#include "mf_plot.h"

#include <cpl_test.h>

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    /* Initialize */
    cpl_errorstate errstate = cpl_errorstate_get();

    /* Get time */
    double cs = cpl_test_get_walltime();


    /* Reading parameter file */

    mfdrv drvpar;
    mf_par_initall(&drvpar);

    const char *parfile = MOLECFIT_SOURCEDIR"/config/molecfit_test_crires.par";
    mf_par_readfile(&drvpar, parfile, 0);

    char line[MF_MAXLEN];

    /*** Single + double plot */

    cpl_table *spectrum1 = cpl_table_new(0);
    cpl_table_new_column(spectrum1,"LAMBDA",CPL_TYPE_DOUBLE);
    cpl_table_new_column(spectrum1,"RADIANCE",CPL_TYPE_DOUBLE);

    cpl_table *spectrum2 = cpl_table_new(0);
    cpl_table_new_column(spectrum2,"LAMBDA",CPL_TYPE_DOUBLE);
    cpl_table_new_column(spectrum2,"RADIANCE",CPL_TYPE_DOUBLE);

    int loop      = 0;                   /* used for loops                  */
    int linecount = 0;                   /* for determining length of input */

    /* Determining cpl table sizes */

    const char *input_filename1 = MOLECFIT_SOURCEDIR"/input/testspec1.dat";
    cpl_msg_info(cpl_func, "%s", input_filename1);
    FILE *specfile1 = fopen(input_filename1, "r");
    while (fgets(line, MF_MAXLEN - 1, specfile1) != NULL) {
        linecount++;
    }
    fclose(specfile1);
    cpl_table_set_size(spectrum1,linecount);
    linecount=0;

    const char *input_filename2 = MOLECFIT_SOURCEDIR"/input/testspec2.dat";
    cpl_msg_info(cpl_func, "%s", input_filename1);
    FILE *specfile2 = fopen(input_filename2,"r");
    while (fgets(line, MF_MAXLEN - 1, specfile2) != NULL) {
        linecount++;
    }
    fclose(specfile2);
    cpl_table_set_size(spectrum2,linecount);

    char *entry;                         /* data point in spectrum file     */
    float lambda=0;                      /* wavelength                      */
    double data_point=0;                 /* data point                      */

    /* Reading spectrum 1 */
    specfile1 = fopen(input_filename1,"r");
    while (fgets(line, MF_MAXLEN - 1, specfile1) != NULL) {

        mf_basic_terminatestring(line);
        entry = strtok(line,"\t");
        lambda=1e4/atof(entry);

        entry = strtok(NULL, "\t" );
        data_point=atof(entry);
        cpl_table_set_double(spectrum1, "LAMBDA",   loop, lambda    );
        cpl_table_set_double(spectrum1, "RADIANCE", loop, data_point);
        loop++;
    }
    fclose(specfile1);
    loop=0;

    /* Reading spectrum 2 */
    specfile2 = fopen(input_filename2,"r");
    while (fgets(line, MF_MAXLEN - 1, specfile2) != NULL) {

        mf_basic_terminatestring(line);
        entry = strtok(line,"\t");
        lambda=1e4/atof(entry);

        entry = strtok(NULL, "\t" );
        data_point=atof(entry);
        cpl_table_set_double(spectrum2, "LAMBDA",   loop, lambda    );
        cpl_table_set_double(spectrum2, "RADIANCE", loop, data_point);
        loop++;
    }
    fclose(specfile2);


    /*** PLOTS ***/

    const char *title1 = "DEMO: mf plot single";
    mf_plot_single(spectrum1, NULL, NULL, title1, &drvpar);
    cpl_test_error(CPL_ERROR_NONE);

    mf_plot_double(spectrum1, spectrum2, &drvpar);
    cpl_test_error(CPL_ERROR_NONE);

    /* Cleanup */
    cpl_table_delete(spectrum1);
    cpl_table_delete(spectrum2);
    cpl_test_error(CPL_ERROR_NONE);


    /*** xyplot ***/
    const char *input_filename = MOLECFIT_SOURCEDIR"/input/testdat_plot_xy.dat";   /* Reading input data: testdata is f(x)=sin(x)/x */
    cpl_msg_info(cpl_func, "%s", input_filename1);

    cpl_table *testdat = cpl_table_new(0);
    cpl_table_new_column(testdat, "x", CPL_TYPE_DOUBLE);
    cpl_table_new_column(testdat, "y", CPL_TYPE_DOUBLE);

    linecount = 0;
    loop      = 0;

    /* Determining cpl table sizes */
    FILE *specfile = fopen(input_filename,"r");
    while (fgets(line, MF_MAXLEN - 1, specfile) != NULL) {
        linecount++;
    }
    fclose(specfile);
    cpl_table_set_size(testdat,linecount);

    specfile = fopen(input_filename, "r");
    while (fgets(line, MF_MAXLEN - 1, specfile) != NULL) {

        mf_basic_terminatestring(line);
        entry = strtok(line,"\t");
        lambda=atof(entry);

        entry = strtok(NULL, "\t" );
        data_point=atof(entry);
        cpl_table_set_double(testdat, "x", loop, lambda);
        cpl_table_set_double(testdat, "y", loop, data_point);
        loop++;
    }
    fclose(specfile);


    /*** PLOT ***/
    const char *title2   = "DEMO: mf plot single";
    const char *x_label2 = "x axis";
    const char *y_label2 = "y axis";
    mf_plot_xy(testdat, x_label2, y_label2, title2, &drvpar);
    cpl_test_error(CPL_ERROR_NONE);

    /* Cleanup */
    cpl_table_delete(testdat);
    cpl_test_error(CPL_ERROR_NONE);


    /*** Cleanup ***/
    mf_par_deleteall(&drvpar);
    cpl_test_error(CPL_ERROR_NONE);

    /* Show time */
    double ce = cpl_test_get_walltime();
    cpl_msg_info(cpl_func, "test_mf_model() -> Run time: %g min\n", (ce - cs) / 60.);

    /* Show errors and return */
    cpl_errorstate_dump(errstate, CPL_TRUE, cpl_errorstate_dump_one);
    return !cpl_errorstate_is_equal(errstate);
}
