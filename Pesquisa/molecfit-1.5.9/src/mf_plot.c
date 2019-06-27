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
 * \file mf_plot.c
 *
 * plotting library for molecfit
 *
 * \author Wolfgang Kausch & ESO In-Kind Team Innsbruck
 *
 */

/*****************************************************************************
 *                                  INCLUDES                                 *
 ****************************************************************************/

#include <mf_plot.h>

/*****************************************************************************
 *                                  DEFINES                                  *
 ****************************************************************************/

#define EXIST W_OK          /* check for file existence, unistd.h, access() */

/*****************************************************************************
 *                                 FUNCTIONS                                 *
 ****************************************************************************/

/*****************************************************************************
 *                                    CODE                                   *
 ****************************************************************************/

/*****************************************************************************
 *                                                                           *
 *  Routine to plot a single spectrum                                        *
 *                                                                           *
 ****************************************************************************/


cpl_error_code mf_plot_single(
    const cpl_table *spec,
    const char      *user_x_label,
    const char      *user_y_label,
    const char      *plot_title,
    const mfdrv     *drvpar)
{
/*!
 * \callgraph
 *
 * This program plots a single spectrum from a two-column CPL table. The first
 * column should contain wavelength information, the second one gives either
 * the flux or the transmission.
 * The second, third and fourth input parameter offer the possibility to give
 * particular x- and y-labels, and a plot title, respectively. In case of NULL
 * vectors the name of the cpl table columns is used for axes labelling.
 *
 * \b INPUT:
 * \param spec          2 column cpl table containing spectrum
 *                      (wavelength [micron], radiance flux/transmission)
 * \param user_x_label  string used for labelling of x-axis
 * \param user_y_label  string used for labelling of y-axis
 * \param plot_title    plot title
 * \param drvpar        structure containing information of the parameter file
 *
 * \b ERRORS:
 * - CPL_ERROR_NONE:           no error occurred
 * - CPL_ERROR_ILLEGAL_INPUT:  # of columns in table not equal two
 *
 */

    FILE *specfile;                 /* output ASCII spec file ptr           */
    FILE *gnufile;                  /* gnuplot driverfile                   */
    char filename[MF_MAXLEN]="sci_single_spec.dat";            /* filename  */
    char gnuname1[MF_MAXLEN]="single_gnufile_wxt.gnu";
    char gnuname2[MF_MAXLEN]="single_gnufile_ps.gnu";
    char gnuname3[MF_MAXLEN]="single_gnufile_x11.gnu";
    char path[MF_MAXLEN];           /* full path to output file             */
    char ps_filename[MF_MAXLEN];

    char tmpdir[MF_MAXLEN];
    char tmpfilename[MF_MAXLEN];
    char spectype[MF_MAXLEN];
    char plot_type[MF_MAXLEN];      /* plot type selection (plot_creation)  */

    char system_call[MF_MAXLEN];    /* system call                          */

    int len=0;                      /* table length                         */
    int run=0;                      /* runnning variable                    */
    int ncol=0;                     /* number of columns in inputspectable  */
    int dummy=0;                    /* dummy return value for system calls  */
    double lambda=0;                /* wavelength of spectrum               */
    double y_value=0;               /* y-value of plot (RAD/TRA)            */
    int dir_exist_flag=0;           /* Checking for existence of tmp dir    */
    /* Plot limits */
    double plot_xmin=0., plot_xmax=0.;

    /* plot labels */
    char x_label[MF_MAXLEN],y_label[MF_MAXLEN];                   /* labels */
    char title[MF_MAXLEN];                                        /* title  */

    cpl_array *column_names;        /* Column names of input table          */

    char err_msg[MF_MAXLEN];        /* error message to be returned         */

    /* directory + filename for output ps file */
    cpl_parameter *outdirpar, *filenamepar, *par;

/*--------------------------------------------------------------------------*/
/* ------------------------------ INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/

    /* Checking table (=input spectra) properties */
    len=cpl_table_get_nrow(spec);   /* length of spectrum                   */
    ncol=cpl_table_get_ncol(spec);  /* # of columns                         */
    column_names = cpl_table_get_column_names(spec);  /* Column names       */
    plot_xmin=cpl_table_get_double(spec,
                            cpl_array_get_string(column_names,0),0,NULL);
    plot_xmax=cpl_table_get_double(spec,
                            cpl_array_get_string(column_names,0),len-1,NULL);

    /* Checking plot title existence */
    if (plot_title == NULL)
    {
        sprintf(title," ");
    }
    else
    {
        sprintf(title,"%s",plot_title);
    }

    /* extract axes labels / title */
    if (ncol == 2)
    {
        if (user_x_label == NULL)
        {
            sprintf(x_label,"{/Symbol l} [micron]");
        }
        else
        {
            sprintf(x_label,"%s",user_x_label);
        }
        if ( !strcmp(cpl_array_get_string(column_names,1),"RADIANCE") )
        {
            if (user_y_label == NULL)
            {
                sprintf(y_label,"Radiance");
            }
            else
            {
                sprintf(y_label,"%s",user_y_label);
            }
            sprintf(spectype,"RADIANCE");
        }
        else
        {
            if (user_y_label == NULL)
            {
                sprintf(y_label,"Transmission");
            }
            else
            {
                sprintf(y_label,"%s",user_y_label);
            }
            sprintf(spectype,"TRANSMISSION");
        }
    }
    else
    {
        sprintf(err_msg,"Number of columns not equal 2 in input spectrum.");
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "%s", err_msg);
    }

    /* Creating tmp-directory */
    sprintf(tmpdir,"__tmpDIRtmp__");
    dir_exist_flag=access(tmpdir,EXIST);
    if (dir_exist_flag == 0){
        cpl_msg_warning(cpl_func,"Directory %s already exists!",tmpdir);
    }
    else
    {
        if ((dummy=mkdir(tmpdir,0777))){};
    }

    /* Checking plot options */
    par=cpl_parameterlist_find(drvpar->parlist, "plot_creation");
    sprintf(plot_type, "%s", cpl_parameter_get_string(par));

    /* writing .dat file containing spectrum information */
    sprintf(tmpfilename,"%s/%s",tmpdir,filename);
    specfile = fopen(tmpfilename,"w");
    for (run=0;run<len;run++)
    {
        lambda=cpl_table_get_double(spec,
                    cpl_array_get_string(column_names,0),run,NULL);
        y_value=cpl_table_get_double(spec,spectype,run,NULL);
        fprintf(specfile,"%5.6g\t%5.6g\n",lambda,y_value);
    }
    fclose(specfile);

    /* Creating wxt terminal gnuplot driver file */
    if (strchr(plot_type, 'W') != NULL) {
        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term wxt\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 13\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set title \"%s\"\n",title);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"set ylabel \"%s\"\n",y_label);
        fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",tmpdir,filename);
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname1);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        dummy=remove(tmpfilename);
    }

    /* Creating postscript terminal gnuplot driver file */
    if (strchr(plot_type, 'P') != NULL) {
        char curdir[MF_MAXLEN];
        cpl_parameter *p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        outdirpar = cpl_parameterlist_find(drvpar->parlist, "output_dir");
        mf_basic_abspath(path, cpl_parameter_get_string(outdirpar), curdir);
        filenamepar = cpl_parameterlist_find(drvpar->parlist, "output_name");

        sprintf(ps_filename,"%s/%s_singleplot.ps", path,
                cpl_parameter_get_string(filenamepar));

        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term postscript enhanced color\n");
        fprintf(gnufile,"set output \"%s\"\n",ps_filename);
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 13\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set title \"%s\"\n",title);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"set ylabel \"%s\"\n",y_label);
        fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",tmpdir,filename);
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname2);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
        dummy=remove(tmpfilename);
    }

    /* Creating x11 terminal gnuplot driver file */
    if (strchr(plot_type, 'X') != NULL) {
        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term x11\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 13\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set title \"%s\"\n",title);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"set ylabel \"%s\"\n",y_label);
        fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",tmpdir,
                                                              filename);
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname3);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        dummy=remove(tmpfilename);
    }

    /* Cleaning */
    mf_basic_initstring(tmpfilename,MF_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename);
    dummy=remove(tmpfilename);
    dummy=rmdir(tmpdir);

    cpl_array_delete(column_names);

    return CPL_ERROR_NONE;
}


/*****************************************************************************
 *                                                                           *
 *  Routine to plot two spectra                                              *
 *                                                                           *
 ****************************************************************************/

cpl_error_code mf_plot_double(const cpl_table *spec1, const cpl_table *spec2,
                              const mfdrv *drvpar)
{
/*!
 * \callgraph
 *
 * This program plots two spectra, e.g. radiance and transmission. Input
 * parameters must be two-column CPL tables, containing wavelength information
 * in col#1 and flux/transmission in col#2. Wavelength info might be in
 * micron or in wavenumber. This must be indicated by appropriate CPL column
 * names, either LAMBDA or WAVENUMBER. If column name #2 is not "RADIANCE",
 * a transmission spectrum is assumed. Flux must be in molecfit compliant
 * physical units [phot/s/m^2/micron/srqsec].
 *
 * \b INPUT:
 * \param spec1   spectrum 1
 * \param spec2   spectrum 2
 * \param drvpar  structure containing information of the parameter file
 *
 * \b ERRORS:
 * - CPL_ERROR_NONE:           no error occurred
 * - CPL_ERROR_ILLEGAL_INPUT:  error in input spectra
 */

    cpl_errorstate err_state;
    cpl_error_code err_code=CPL_ERROR_NONE;

    FILE *specfile;                 /* output ASCII spec file ptr           */
    FILE *gnufile;                  /* gnuplot driverfile                   */
    char filename1[MF_MAXLEN]="plot_spec1.dat";     /* input spec 1         */
    char filename2[MF_MAXLEN]="plot_spec2.dat";     /* input spec 2         */
    char gnuname1[MF_MAXLEN] ="double_gnufile_wxt.gnu";
    char gnuname2[MF_MAXLEN] ="double_gnufile_ps.gnu";
    char gnuname3[MF_MAXLEN] ="double_gnufile_x11.gnu";

    char tmpdir[MF_MAXLEN];
    char tmpfilename[MF_MAXLEN];
    char path[MF_MAXLEN];           /* full path to output file             */
    char ps_filename[MF_MAXLEN];    /* Name of postscript file              */
    char plot_type[MF_MAXLEN];      /* plot type selection (plot_creation)  */

    char spectype1[MF_MAXLEN], spectype2[MF_MAXLEN];       /* Rad or Tra    */
    char wavetype1[MF_MAXLEN], wavetype2[MF_MAXLEN];       /* lambda or wn  */
    char system_call[MF_MAXLEN];

    int len1=0, len2=0;             /* table lengtes                        */
    int ncol1=0, ncol2=0;           /* number of columns in inputspectables */
    int dummy=0;                    /* dummy return value for system calls  */
    int run=0;                      /* runnning variable                    */
    int gridcheckflag=0;            /* flag for checking lambda grid        */
    int dir_exist_flag=0;           /* Checking for existence of tmp dir    */

    /* Plot limits */
    double plot_xmin1=0., plot_xmax1=0., plot_xmin2=0., plot_xmax2=0.;

    double lambda=0;                /* wavelength of spectrum               */
    double y_value=0;               /* y-value of plot (RAD/TRA)            */

    char title1[MF_MAXLEN];                                  /* title       */
    char x_label1[MF_MAXLEN],y_label1[MF_MAXLEN];            /* labels      */
    char title2[MF_MAXLEN];                                  /* title       */
    char x_label2[MF_MAXLEN],y_label2[MF_MAXLEN];            /* labels      */

    cpl_array *column_names1, *column_names2; /* Col. names of input tables */

    char err_msg[MF_MAXLEN];        /* error message to be returned         */

    /* directory + filename for output ps file */
    cpl_parameter *outdirpar, *filenamepar, *par;

/*--------------------------------------------------------------------------*/
/* ------------------------------ INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/

    /* checking error state */
    if ((err_state = cpl_errorstate_get())){};

    /* CHECKING INPUT ------------------------------------------------------*/
    /* Input spectra will be checked against
         - same size (#cols, #rows)
         - same CPL table structure
         - same wavelength grid
         - type of input specta
    */
    /* Checking table (=input spectra) properties */
    len1=cpl_table_get_nrow(spec1);   /* length of spectrum 1               */
    ncol1=cpl_table_get_ncol(spec1);  /* # of columns spec 1                */
    column_names1 = cpl_table_get_column_names(spec1);  /* Col names spec 1 */

    len2=cpl_table_get_nrow(spec2);   /* length of spectrum 2               */
    ncol2=cpl_table_get_ncol(spec2);  /* # of columns spec 2                */
    column_names2 = cpl_table_get_column_names(spec2);  /* Col names spec 2 */

    plot_xmin1=cpl_table_get_double(spec1,"LAMBDA",0,NULL);
    plot_xmax1=cpl_table_get_double(spec1,"LAMBDA",len1-1,NULL);
    plot_xmin2=cpl_table_get_double(spec2,"LAMBDA",0,NULL);
    plot_xmax2=cpl_table_get_double(spec2,"LAMBDA",len2-1,NULL);

    /* Checking length of input CPL tables */
    if (len1 != len2)
    {
        cpl_msg_warning(cpl_func,"Input spectra do not have the same size.");
    }

    /* Checking structure of input CPL tables */
    if ( cpl_table_compare_structure(spec1,spec2) )
    {
        cpl_msg_warning(cpl_func,"Input spectra do not have the same "
                                 "structure.");
    }

    /* Checking type of input spectra and set appropriate labelling tags    */

    /* Checking for # columns in inputspec1 + determing type */
    if (ncol1 == 2)
    {
        sprintf(x_label1,"{/Symbol l} [micron]");

        if ( !strcmp(cpl_array_get_string(column_names1,1),"RADIANCE") )
        {
            sprintf(spectype1,"RADIANCE");
            sprintf(y_label1,"Radiance");
            sprintf(title1,"Radiance");
        }
        else
        {
            sprintf(spectype1,"TRANSMISSION");
            sprintf(y_label1,"Transmission");
            sprintf(title1,"Transmission");
        }
        if ( !strcmp(cpl_array_get_string(column_names1,0),"WAVENUMBER") )
        {
            sprintf(wavetype1,"WAVENUMBER");
        }
        else
        {
            sprintf(wavetype1,"LAMBDA");
        }
    }
    else
    {
        sprintf(err_msg,"Number of columns not equal 2 in input spectrum 1.");
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "%s", err_msg);
    }

    /* Checking for # columns in inputspec2 + determing type */
   if (ncol2 == 2)
    {
        sprintf(x_label2,"wavelength [micron]");
        if ( !strcmp(cpl_array_get_string(column_names2,1),"RADIANCE") )
        {
            sprintf(spectype2,"RADIANCE");
            sprintf(y_label2,"Radiance");
            sprintf(title2,"Radiance");
        }
        else
        {
            sprintf(spectype2,"TRANSMISSION");
            sprintf(y_label2,"Transmission");
            sprintf(title2,"Transmission");
        }
        if ( !strcmp(cpl_array_get_string(column_names2,0),"WAVENUMBER") )
        {
            sprintf(wavetype2,"WAVENUMBER");
        }
        else
        {
            sprintf(wavetype2,"LAMBDA");
        }

    }
    else
    {
        sprintf(err_msg,"Number of columns not equal 2 in input spectrum 2.");
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "%s", err_msg);
    }

    /* Checking if spectrum types are identical */
    if ( strcmp(spectype1,spectype2) )
    {
        cpl_msg_warning(cpl_func,"Input spectra are of different types "
                                 "(RAD/TRA).");
    }

    /* Checking if wavelength grid is identical */
    run=0;
    for (run=0;run<len1;run++)
    {
        if (cpl_table_get_double(spec1,
                cpl_array_get_string(column_names1,0),run,NULL) !=
            cpl_table_get_double(spec2,
                cpl_array_get_string(column_names2,0),run,NULL) )
        {
            gridcheckflag=1;
        }
    }
    if (gridcheckflag)
    {
        cpl_msg_warning(cpl_func,"Input spectra have differing wavelength "
                                 "grid.");
    }

    /* Checking plot options */
    par=cpl_parameterlist_find(drvpar->parlist, "plot_creation");
    sprintf(plot_type, "%s", cpl_parameter_get_string(par));

    /* Checking / Creating tmp-directory */
    sprintf(tmpdir,"__tmpDIRtmp__");
    dir_exist_flag=access(tmpdir,EXIST);
    if (dir_exist_flag == 0){
        cpl_msg_warning(cpl_func,"Directory %s already exists!",tmpdir);
    }
    else
    {
        if ((dummy=mkdir(tmpdir,0777))){};
    }

    /* writing .dat files containing spectrum information */
    mf_basic_initstring(tmpfilename,MF_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename1);
    specfile = fopen(tmpfilename,"w"); /* Science spectrum */
    for (run=0;run<len1;run++)
    {
        if ( !strcmp(wavetype1,"WAVENUMBER") )
        {
            lambda=1e4/cpl_table_get_double(spec1,"WAVENUMBER",run,NULL);
            lambda=1e4/cpl_table_get_double(spec1,
                   cpl_array_get_string(column_names1,0),run,NULL);
            y_value=cpl_table_get_double(spec1,spectype1,run,NULL);
            fprintf(specfile,"%5.6g\t%5.6g\n",lambda,y_value);
        }
        else
        {
            lambda=cpl_table_get_double(spec1,
                   cpl_array_get_string(column_names1,0),run,NULL);
            y_value=cpl_table_get_double(spec1,spectype1,run,NULL);
            fprintf(specfile,"%5.6g\t%5.6g\n",lambda,y_value);

        }
    }
    fclose(specfile);

    mf_basic_initstring(tmpfilename,MF_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename2);
    specfile = fopen(tmpfilename,"w"); /* Fit spectrum */
    for (run=0;run<len2;run++)
    {
        if ( !strcmp(wavetype1,"WAVENUMBER") )
        {
            lambda=1e4/cpl_table_get_double(spec2,
                   cpl_array_get_string(column_names2,0),run,NULL);
            y_value=cpl_table_get_double(spec2,spectype2,run,NULL);
            fprintf(specfile,"%5.6g\t%5.6g\n",lambda,y_value);
        }
        else
        {
            lambda=cpl_table_get_double(spec2,
                   cpl_array_get_string(column_names2,0),run,NULL);
            y_value=cpl_table_get_double(spec2,spectype2,run,NULL);
            fprintf(specfile,"%5.6g\t%5.6g\n",lambda,y_value);

        }
    }
    fclose(specfile);

/****************************************************************************/

    /* Creating wxt terminal gnuplot file */
    if (strchr(plot_type, 'W') != NULL) {
        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term wxt\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"# Plotting\n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 12\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin1,plot_xmax1);
        fprintf(gnufile,"set termoption font \"Times,7\"\n");
        fprintf(gnufile,"set multiplot layout 2,1 \n");
        fprintf(gnufile,"set title \"Spectrum 1\"\n");
        fprintf(gnufile,"set xlabel \"{/Symbol l} [micron]\"\n");
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label1);
        fprintf(gnufile,"set style data boxes\n");
        fprintf(gnufile,"plot '%s/%s' using 1:2 title "
                        " \"spectrum 1\" with lines lt 8\n",tmpdir,filename1);
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label2);
        fprintf(gnufile,"set xlabel \"{/Symbol l} [micron]\"\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin2,plot_xmax2);
        fprintf(gnufile,"set title \"Spectrum 2\"\n");
        fprintf(gnufile,"plot '%s/%s' "
                        "using 1:2 title \"spectrum 2\" with lines lt 8\n",
                        tmpdir,filename2);
        fprintf(gnufile,"unset multiplot    \n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname1);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        dummy=remove(tmpfilename);
    }

    /* Creating postscript terminal gnuplot driver file */
    if (strchr(plot_type, 'P') != NULL) {
        char curdir[MF_MAXLEN];
        cpl_parameter *p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        outdirpar = cpl_parameterlist_find(drvpar->parlist, "output_dir");
        mf_basic_abspath(path, cpl_parameter_get_string(outdirpar), curdir);
        filenamepar = cpl_parameterlist_find(drvpar->parlist, "output_name");

        sprintf(ps_filename,"%s/%s_doubleplot.ps", path,
                cpl_parameter_get_string(filenamepar));

        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
        gnufile = fopen(tmpfilename,"w");

        fprintf(gnufile,"set term postscript enhanced color\n");
        fprintf(gnufile,"set output \"%s\"\n", ps_filename);
        fprintf(gnufile,"# Plotting\n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 12\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin1,plot_xmax1);
        fprintf(gnufile,"set termoption font \"Times,12\"\n");
        fprintf(gnufile,"set multiplot layout 2,1 \n");
        fprintf(gnufile,"set title \"Spectrum 1\"\n");
        fprintf(gnufile,"set xlabel \"{/Symbol l} [micron]\"\n");
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label1);
        fprintf(gnufile,"set style data boxes\n");
        fprintf(gnufile,"plot '%s/%s' using 1:2 title "
                        " \"spectrum 1\" with lines ls 1\n",tmpdir,filename1);
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label2);
        fprintf(gnufile,"set xlabel \"{/Symbol l} [micron]\"\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin2,plot_xmax2);
        fprintf(gnufile,"set title \"Spectrum 2\"\n");
        fprintf(gnufile,"plot '%s/%s' "
                        "using 1:2 title \"spectrum 2\" with lines ls 1\n",
                        tmpdir,filename2);
        fprintf(gnufile,"unset multiplot    \n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname2);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
        dummy=remove(tmpfilename);
    }

    /* Creating x11 terminal gnuplot file */
    if (strchr(plot_type, 'X') != NULL) {
        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term x11\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"# Plotting\n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 12\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin1,plot_xmax1);

        fprintf(gnufile,"set termoption font \"Times,7\"\n");
        fprintf(gnufile,"set multiplot layout 2,1 \n");
        fprintf(gnufile,"set title \"Spectrum 1\"\n");
        fprintf(gnufile,"set xlabel \"{/Symbol l} [micron]\"\n");
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label1);
        fprintf(gnufile,"set style data boxes\n");
        fprintf(gnufile,"plot '%s/%s' using 1:2 title "
                        " \"spectrum 1\" with lines lt 8\n",tmpdir,filename1);
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set ylabel \"%s\" offset 1,0\n",y_label2);
        fprintf(gnufile,"set xlabel \"{/Symbol l} [micron]\"\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin2,plot_xmax2);

        fprintf(gnufile,"set title \"Spectrum 2\"\n");
        fprintf(gnufile,"plot '%s/%s' "
                        "using 1:2 title \"spectrum 2\" with lines lt 8\n",
                        tmpdir,filename2);
        fprintf(gnufile,"unset multiplot    \n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname3);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        dummy=remove(tmpfilename);
    }

    /* Cleaning */
    mf_basic_initstring(tmpfilename,MF_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename1);
    dummy=remove(tmpfilename);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename2);
    dummy=remove(tmpfilename);
    dummy=rmdir(tmpdir);

    cpl_array_delete(column_names1);
    cpl_array_delete(column_names2);

    return err_code;
}


/*****************************************************************************
 *                                                                           *
 *  Routine to plot a simple x vs. y plot                                    *
 *                                                                           *
 ****************************************************************************/

cpl_error_code mf_plot_xy(
    const cpl_table *spec,
    const char      *user_x_label,
    const char      *user_y_label,
    const char      *plot_title,
    const mfdrv     *drvpar)
{
/*!
 * \callgraph
 *
 * This program creates a simple x vs y plot from a two-column CPL table. The
 * first column corresponds to the x-axis, the second column to the y axis.
 * The second, third and fourth input parameter offer the possibility to give
 * particular x- and y-labels, and a title, respectively. In case of NULL
 * vectors the name of the cpl table columns is used for labelling or no
 * title is added.
 *
 * \b INPUT:
 * \param spec          2 column cpl table
 * \param user_x_label  string used for labelling of x-axis
 * \param user_y_label  string used for labelling of y-axis
 * \param plot_title    string used for plot title
 * \param drvpar        structure containing information of the parameter file
 *
 * \b ERRORS:
 * - CPL_ERROR_NONE: no error occurred
 * - CPL_ERROR_ILLEGAL_INPUT
 */

    FILE *specfile;                 /* output ASCII spec file ptr           */
    FILE *gnufile;                  /* gnuplot driverfile                   */
    char filename[MF_MAXLEN]="xy_plot.dat";     /* data filename            */
    char gnuname1[MF_MAXLEN]="xy_plot_gnufile_wxt.gnu";
    char gnuname2[MF_MAXLEN]="xy_plot_gnufile_ps.gnu";
    char gnuname3[MF_MAXLEN]="xy_plot_gnufile_x11.gnu";
    char path[MF_MAXLEN];           /* full path to output file             */
    char ps_filename[MF_MAXLEN];    /* Name of postscript file              */
    char tmpdir[MF_MAXLEN];
    char tmpfilename[MF_MAXLEN];
    char plot_type[MF_MAXLEN];      /* plot type selection (plot_creation)  */

    char system_call[MF_MAXLEN];    /* system call                          */
    int len=0;                      /* table length                         */
    int run=0;                      /* runnning variable                    */
    int ncol=0;                     /* number of columns in inputspectable  */
    int dummy=0;                    /* dummy return value for system calls  */
    double x_value=0;               /* wavelength of spectrum               */
    double y_value=0;               /* y-value of plot (RAD/TRA)            */
    int dir_exist_flag=0;           /* Checking for existence of tmp dir    */

    double plot_xmin=0., plot_xmax=0., plot_ymin=0., plot_ymax=0.;

    /* plot labels + title */
    char x_label[MF_MAXLEN],y_label[MF_MAXLEN],title[MF_MAXLEN];

    cpl_array *column_names;        /* Column names of input table          */
    char err_msg[MF_MAXLEN];     /* error message to be returned            */

    /* directory + filename for output ps file */
    cpl_parameter *outdirpar, *filenamepar, *par;

/*--------------------------------------------------------------------------*/
/* ------------------------------ INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/

    /* Checking table (=input spectra) properties */
    len=cpl_table_get_nrow(spec);   /* length of spectrum                   */
    ncol=cpl_table_get_ncol(spec);  /* # of columns                         */
    column_names = cpl_table_get_column_names(spec);  /* Column names       */

    plot_xmin=cpl_table_get_double(spec,
                            cpl_array_get_string(column_names,0),0,NULL);
    plot_xmax=cpl_table_get_double(spec,
                            cpl_array_get_string(column_names,0),len-1,NULL);
    plot_ymin=cpl_table_get_column_min(spec,
                            cpl_array_get_string(column_names,1));
    plot_ymax=cpl_table_get_column_max(spec,
                            cpl_array_get_string(column_names,1));

    mf_basic_initstring(x_label,MF_MAXLEN);
    mf_basic_initstring(y_label,MF_MAXLEN);
    /* extract axes labels / title */
    if (plot_title == NULL)
    {
        sprintf(title," ");
        cpl_msg_warning(cpl_func,"No title for plot given.");
    }
    else
    {
        sprintf(title,"%s",plot_title);
    }

    if (ncol == 2)
    {
        if (user_x_label == NULL)                      // checking title
        {
            sprintf(y_label,"%s",cpl_array_get_string(column_names,0));
            cpl_msg_warning(cpl_func,"No x-label for plot given, will "
                                     "use CPL columnname #1 instead.");
        }
        else
        {
            sprintf(x_label,"%s",user_x_label);
        }

        if (user_y_label == NULL)                      // checking y label
        {
            sprintf(y_label,"%s",cpl_array_get_string(column_names,1));
            cpl_msg_warning(cpl_func,"No y-label for plot given, will "
                                     "use CPL columnname #2 instead.");
        }
        else
        {
            sprintf(y_label,"%s",user_y_label);
        }

    }
    else
    {
        sprintf(err_msg,"Number of columns not equal 2 in input CPL "
                        "table...");
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "%s", err_msg);
    }

    /* Checking plot options */
    par=cpl_parameterlist_find(drvpar->parlist, "plot_creation");
    sprintf(plot_type, "%s", cpl_parameter_get_string(par));

    /* Checking / Creating tmp-directory */
    sprintf(tmpdir,"__tmpDIRtmp__");
    dir_exist_flag=access(tmpdir,EXIST);
    if (dir_exist_flag == 0){
        cpl_msg_warning(cpl_func,"Directory %s already exists!",tmpdir);
    }
    else
    {
        if ((dummy=mkdir(tmpdir,0777))){};
    }

    mf_basic_initstring(tmpfilename,MF_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename);
    specfile = fopen(tmpfilename,"w"); /* Science spectrum */
    for (run=0;run<len;run++)
    {
        x_value=cpl_table_get_double(spec,"x",run,NULL);
        y_value=cpl_table_get_double(spec,"y",run,NULL);
        fprintf(specfile,"%5.6g\t%5.6g\n",x_value,y_value);
    }
    fclose(specfile);

    /* Creating wxt terminal gnuplot driver file  */
    /* writing .dat file containing spectrum information */
    if (strchr(plot_type, 'W') != NULL) {
        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);

        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term wxt\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 11\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);

        fprintf(gnufile,"set title \"%s\"\n",title);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"set ylabel \"%s\"\n",y_label);
        fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",tmpdir,filename);
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname1);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        dummy=remove(tmpfilename);
    }

    /* Creating postscript terminal gnuplot driver file */
    if (strchr(plot_type, 'P') != NULL) {
        char curdir[MF_MAXLEN];
        cpl_parameter *p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        outdirpar = cpl_parameterlist_find(drvpar->parlist, "output_dir");
        mf_basic_abspath(path, cpl_parameter_get_string(outdirpar), curdir);
        filenamepar = cpl_parameterlist_find(drvpar->parlist, "output_name");

        sprintf(ps_filename,"%s/%s_plot_xy.ps", path,
                cpl_parameter_get_string(filenamepar));

        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term postscript enhanced color\n");
        fprintf(gnufile,"set output \"%s\"\n",ps_filename);
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 11\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);

        fprintf(gnufile,"set title \"%s\"\n",title);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"set ylabel \"%s\"\n",y_label);
        fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",tmpdir,filename);
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname2);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
        dummy=remove(tmpfilename);
    }

    /* Creating x11 terminal gnuplot driver file  */
    if (strchr(plot_type, 'X') != NULL) {
        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term x11\n");
        fprintf(gnufile,"set termoption enhanced  \n");
        fprintf(gnufile,"set nokey\n");
        fprintf(gnufile,"set tmargin 2\n");
        fprintf(gnufile,"set bmargin 5\n");
        fprintf(gnufile,"set lmargin 11\n");
        fprintf(gnufile,"set rmargin 3\n");
        fprintf(gnufile,"set title \"%s\"\n",title);
        fprintf(gnufile,"set xlabel \"%s\"\n",x_label);
        fprintf(gnufile,"set ylabel \"%s\"\n",y_label);
        fprintf(gnufile,"set xrange [ %g : %g ]\n",plot_xmin,plot_xmax);
        fprintf(gnufile,"set yrange [ %g : %g ]\n",plot_ymin,plot_ymax);

        fprintf(gnufile,"plot '%s/%s' using 1:2 with lines\n",tmpdir,filename);
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname3);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        dummy=remove(tmpfilename);
}

    /* Cleaning */
    mf_basic_initstring(tmpfilename,MF_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename);
    dummy=remove(tmpfilename);
    dummy=rmdir(tmpdir);

    cpl_array_delete(column_names);

    return CPL_ERROR_NONE;
}


/*****************************************************************************
 *                                                                           *
 *  Routine to plot a histogram                                              *
 *                                                                           *
 ****************************************************************************/

cpl_error_code mf_plot_hist(const cpl_table *histdat, char *x_col,
                            char *y_col, char *plot_title,
                            const mfdrv *drvpar)
{
/*!
 * \callgraph
 *
 * This program creates a histogram from a CPL table. This table must have
 * two columns labelled as "bins" and "counts".
 *
 * \b INPUT:
 * \param histdat     2 column cpl table
 * \param x_col       string used for labelling of x-axis
 * \param y_col       string used for labelling of y-axis
 * \param plot_title  string used for plot title
 * \param drvpar      structure containing information of the parameter file
 *
 * \b ERRORS:
 * - CPL_ERROR_NONE: no error occurred
 * - CPL_ERROR_ILLEGAL_INPUT
 */


    cpl_errorstate err_state;
    cpl_error_code err_code=CPL_ERROR_NONE;

    FILE *specfile;                 /* output ASCII spec file ptr           */
    FILE *gnufile;                  /* gnuplot driverfile                   */
    char filename1[MF_MAXLEN]="hist_data.dat";     /* data for histogram    */
    char gnuname1[MF_MAXLEN]="hist_plot_wxt.gnu";
    char gnuname2[MF_MAXLEN]="hist_plot_ps.gnu";
    char gnuname3[MF_MAXLEN]="hist_plot_x11.gnu";
    char path[MF_MAXLEN];           /* full path to output file             */
    char ps_filename[MF_MAXLEN];    /* Name of postscript file              */
    char tmpdir[MF_MAXLEN];
    char tmpfilename[MF_MAXLEN];
    char system_call[MF_MAXLEN];    /* string for system call               */
    char plot_type[MF_MAXLEN];      /* plot type selection (plot_creation)  */

    int  dir_exist_flag=0;          /* Checking for existence of tmp dir    */
    int  dummy=0;                   /* dummy return value for system calls  */
    int run=0;                      /* runnning variable                    */

//     char x_label[MF_MAXLEN],y_label[MF_MAXLEN];                /* labels */
    char title[MF_MAXLEN];                                        /* title  */
    char err_msg[MF_MAXLEN];        /* error message to be returned         */

    int  len1=0, ncol1=0;           /* Dimension of cpl-table 'histdat'     */

    /* directory + filename for output ps file */
    cpl_parameter *outdirpar, *filenamepar, *par;

    cpl_array *column_names1; /* Col. names of input tables */


/*--------------------------------------------------------------------------*/
/* ------------------------------ INITIALISING -----------------------------*/
/*--------------------------------------------------------------------------*/

    /* checking error state */
    if ((err_state = cpl_errorstate_get())){};

    /* CHECKING INPUT ------------------------------------------------------*/
    /* Checking table (=input data) properties                           */
    len1=cpl_table_get_nrow(histdat);   /* # of bins                        */
    ncol1=cpl_table_get_ncol(histdat);  /* # of columns spectrum            */
    column_names1 = cpl_table_get_column_names(histdat); /* Col names       */

    /* Checking whether 2 columns for input file exist */
    if (ncol1 != 2)
    {
        sprintf(err_msg,"Number of columns not equal 2 in data file...");
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                            "%s", err_msg);
    }

    /* Checking plot title / axes label existence */
    if (plot_title == NULL)
    {
        sprintf(title,"Histogram");
    }
    else
    {
        sprintf(title,"%s",plot_title);
    }

    if (x_col == NULL)
    {
        sprintf(title,"bins");
    }
    else
    {
        sprintf(title,"%s",x_col);
    }

    if (y_col == NULL)
    {
        sprintf(title,"counts per bin");
    }
    else
    {
        sprintf(title,"%s",y_col);
    }

    /* Checking / Creating tmp-directory */
    sprintf(tmpdir,"__tmpDIRtmp__");
    dir_exist_flag=access(tmpdir,EXIST);
    if (dir_exist_flag == 0){
        cpl_msg_warning(cpl_func,"Directory %s already exists!",tmpdir);
    }
    else
    {
        if ((dummy=mkdir(tmpdir,0777))){};
    }

    /* writing .dat files containing spectrum information                   */
    /* Science spectrum */
    mf_basic_initstring(tmpfilename,MF_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename1);
    specfile = fopen(tmpfilename,"w");
    for (run=0;run<len1;run++)
    {
        fprintf(specfile,"%5.3g\t%i\n",
                cpl_table_get_double(histdat,"bins",run,NULL),
                cpl_table_get_int(histdat,"counts",run,NULL));
    }
    fclose(specfile);

    /* Checking plot options */
    par=cpl_parameterlist_find(drvpar->parlist, "plot_creation");
    sprintf(plot_type, "%s", cpl_parameter_get_string(par));

    /* Creating wxt terminal gnuplot driver file */
    if (strchr(plot_type, 'W') != NULL) {
        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term wxt\n");
        fprintf(gnufile,"set termoption enhanced\n");
        fprintf(gnufile,"set title \"%s\"\n",plot_title);
        fprintf(gnufile,"set termoption font \"Times,9\"\n");
        fprintf(gnufile,"set xlabel '%s'\n",x_col);
        fprintf(gnufile,"set ylabel '%s'\n",y_col);
        fprintf(gnufile,"set boxwidth 0.75 absolute\n");
        fprintf(gnufile,"set style fill solid 1.00 border -1\n");
        fprintf(gnufile,"set style histogram rowstacked\n");
        fprintf(gnufile,"set style data histograms\n");
        fprintf(gnufile,"plot '%s/%s' u 1:2 smooth frequency with histeps"
                        " t \"%s\"\n",tmpdir,filename1,y_col);

        fprintf(gnufile,"unset xlabel\n");
        fprintf(gnufile,"unset ylabel\n");
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname1);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname1);
        dummy=remove(tmpfilename);
    }

    /* Creating postscript terminal gnuplot driver file */
    if (strchr(plot_type, 'P') != NULL) {
        char curdir[MF_MAXLEN];
        cpl_parameter *p = cpl_parameterlist_find(drvpar->parlist, "curdir");
        strncpy(curdir, cpl_parameter_get_string(p), MF_MAXLEN);
        outdirpar = cpl_parameterlist_find(drvpar->parlist, "output_dir");
        mf_basic_abspath(path, cpl_parameter_get_string(outdirpar), curdir);
        filenamepar = cpl_parameterlist_find(drvpar->parlist, "output_name");

        sprintf(ps_filename,"%s/%s_histogram.ps", path,
                cpl_parameter_get_string(filenamepar));

        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term postscript enhanced color\n");
        fprintf(gnufile,"set output \"%s\"\n",ps_filename);
        fprintf(gnufile,"set title \"%s\"\n",plot_title);
        fprintf(gnufile,"set termoption font \"Times,9\"\n");
        fprintf(gnufile,"set xlabel '%s'\n",x_col);
        fprintf(gnufile,"set ylabel '%s'\n",y_col);
        fprintf(gnufile,"set boxwidth 0.9 absolute\n");
        fprintf(gnufile,"set style fill solid 1.00 border -1\n");
        fprintf(gnufile,"set style data histograms\n");
        fprintf(gnufile,"plot '%s/%s' u 1:2 smooth frequency with histeps"
                        " t \"%s\"\n",tmpdir,filename1,y_col);
        fprintf(gnufile,"unset xlabel\n");
        fprintf(gnufile,"unset ylabel\n");
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname2);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname2);
        dummy=remove(tmpfilename);
    }

    /* Creating x11 terminal gnuplot driver file */
    if (strchr(plot_type, 'X') != NULL) {
        mf_basic_initstring(tmpfilename,MF_MAXLEN);
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        gnufile = fopen(tmpfilename,"w");
        fprintf(gnufile,"set term x11\n");
        fprintf(gnufile,"set termoption enhanced\n");

        fprintf(gnufile,"set title \"%s\"\n",plot_title);
        fprintf(gnufile,"set termoption font \"Times,9\"\n");
        fprintf(gnufile,"set xlabel '%s'\n",x_col);
        fprintf(gnufile,"set ylabel '%s'\n",y_col);
        fprintf(gnufile,"set boxwidth 0.75 absolute\n");
        fprintf(gnufile,"set style fill solid 1.00 border -1\n");
        fprintf(gnufile,"set style histogram rowstacked\n");
        fprintf(gnufile,"set style data histograms\n");
        fprintf(gnufile,"plot '%s/%s' u 1:2 smooth frequency with histeps"
                        " t \"%s\"\n",tmpdir,filename1,y_col);

        fprintf(gnufile,"unset xlabel\n");
        fprintf(gnufile,"unset ylabel\n");
        fprintf(gnufile,"\n");
        fclose(gnufile);

        /* Calling gnuplot */
        sprintf(system_call,"gnuplot -persist %s/%s",tmpdir,gnuname3);
        dummy=system(system_call);

        /* Cleaning */
        sprintf(tmpfilename,"%s/%s",tmpdir,gnuname3);
        dummy=remove(tmpfilename);
    }


    /* Cleaning */
    mf_basic_initstring(tmpfilename,MF_MAXLEN);
    sprintf(tmpfilename,"%s/%s",tmpdir,filename1);
    dummy=remove(tmpfilename);
    dummy=rmdir(tmpdir);

    cpl_array_delete(column_names1);

    return err_code;
}

/**@}*/
