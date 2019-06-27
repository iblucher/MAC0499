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

/*
 * test_mf_createatm.c
 *
 *  Created on: Jul 8, 2010
 *      Author: barden
 */

#include <cpl_test.h>
#include <mf_basic.h>

int main(void) {
    double date[] = {2000.,
                     2002.00273972602739726,
                     2002.0794520547945206,
                     2002.0821917808219178,
                     2002.0849315068493151,
                     2002.0876712328767123,
                     2002.158904109589041,
                     2002.161643835616438,
                     2002.243835616438356,
                     2002.246575342465753,
                     2002.328767123287671,
                     2002.413698630136986,
                     2002.495890410958904,
                     2002.580821917808219,
                     2002.997260273972603,
                     2002.000000000000001,
                     2002.00273,
                     2002.07945,
                     2002.082195,
                     2002.08493,
                     2002.08768,
                     2002.15890,
                     2002.16164,
                     2002.243835,
                     2002.24657,
                     2002.32877,
                     2002.41369,
                     2002.49589,
                     2002.58083,
                     2002.99726};

//    double rest;

    int day = 0, month = 0, year = 0;
    int hh = 0, mm = 0;
    double ss = 0;

    int i = 0;

    /************************************************************************
     * variables for mf_basic_rhum2ppmv
     ************************************************************************/
    cpl_table *profile;

    cpl_array *temp, *pres, *rhum, *ppmv;

    double d_temp, d_pres, d_rhum, d_ppmv, val;

    int n = 210, tlow = 123; //, col;

    /************************************************************************
     * variables for mf_basic_planck
     ************************************************************************/
//    double planck_true[] = {0.00000000000000000000,
//                            0.00560776004567742348,
//                            0.96991318464279174805};
    double planck_temp = 280;

    cpl_array *planck_wave, *planck_res;

    /************************************************************************
     * variables for mf_basic_greg2jd & mf_basic_jd2greg
     ************************************************************************/

    const long jd1 = 2455170, // 2009 12 4
               jd2 = 2454526, // 2008 2 29
               jd3 = 2415080, // 1900 3 1
               jd4 = 2454497; // 2008 1 31
    long jd = 0;

    /************************************************************************
     * variables for mf_basic_abspath
     ************************************************************************/

    char cwd[MF_MAXLEN],
         cwd_org[MF_MAXLEN],
         dir0[MF_MAXLEN] = "/dir1/./.././../dir2/../dir3/dir4/dir5/../..",
         dir1[MF_MAXLEN] = "",
         dir2[MF_MAXLEN] = ".",
         dir3[MF_MAXLEN] = "./",
         dir4[MF_MAXLEN] = "..",
         dir5[MF_MAXLEN] = "../",
         dir6[MF_MAXLEN] = "../..",
         dir7[MF_MAXLEN] = "../../",
         dir8[MF_MAXLEN] = "/",
         dir9[MF_MAXLEN] = "config",
         dir10[MF_MAXLEN] = "config/",
         dir11[MF_MAXLEN] = "../../config",
         dir12[MF_MAXLEN] = "../../config/",
         dir13[MF_MAXLEN] = "/config/",
         truedir[MF_MAXLEN],
         outdir[MF_MAXLEN];

    /************************************************************************
     * variables for mf_basic_dirslash
     ************************************************************************/

    char noslash[MF_MAXLEN] = "bla",
         slash[MF_MAXLEN] = "bla/";

    /************************************************************************
     * variables for mf_basic_absfile
     ************************************************************************/

    // already defined for mf_basic_abspath
    // char cwd[MF_MAXLEN],
    //      cwd_org[MF_MAXLEN]
    char file1[MF_MAXLEN] = "/test/filename",
         file2[MF_MAXLEN] = ".././filename",
         file3[MF_MAXLEN] = "/filename",
         file4[MF_MAXLEN] = "filename",
         truefile[MF_MAXLEN],
         outfile[MF_MAXLEN];

    /************************************************************************
     * other variables
     ************************************************************************/

    // already defined for mf_basic_fracyear2date
    //int year = 0, month = 0, day = 0;
    const int yr1 = 2009, mn1 = 12, dy1 =  4,
              yr2 = 2008, mn2 =  2, dy2 = 29,
              yr3 = 1900, mn3 =  3, dy3 = 1,
              yr4 = 2008, mn4 =  1, dy4 = 31;

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /************************************************************************
     * testing mf_basic_fracyear2date
     ************************************************************************/
    for (i = 0; i < 30; ++i) {
        if (mf_basic_fracyear2date(&year, &month, &day, &hh,
                                      &mm, &ss, &(date[i]))) {}
//        rest = mf_basic_fracyear2date(&year, &month, &day, &hh,
//                                      &mm, &ss, &(date[i]));
//        printf("%lf: %2i/%2i/%4i %2i:%2i:%5.2lf + %20.15lf\n",
//               date[i], day, month, year, hh, mm, ss, rest);
    }

    /************************************************************************
     * testing mf_basic_rhum2ppmv
     ************************************************************************/
    temp = cpl_array_new(3, CPL_TYPE_DOUBLE);
    pres = cpl_array_new(3, CPL_TYPE_DOUBLE);
    rhum = cpl_array_new(3, CPL_TYPE_DOUBLE);
    ppmv = cpl_array_new(3, CPL_TYPE_DOUBLE);

    cpl_array_set_double(temp, 0, 300);
    cpl_array_set_double(temp, 1, 290);
    cpl_array_set_double(temp, 2, 280);

    cpl_array_set_double(pres, 0, 900);
    cpl_array_set_double(pres, 1, 800);
    cpl_array_set_double(pres, 2, 700);

    cpl_array_set_double(rhum, 0, 5);
    cpl_array_set_double(rhum, 1, 40);
    cpl_array_set_double(rhum, 2, 20);

    if (mf_basic_rhum2ppmv_old(temp, pres, rhum, ppmv) != CPL_ERROR_NONE) {
        cpl_msg_warning(cpl_func, "mf_basic_rhum2ppmv() "
                                  "encountered a problem.");
    } else {
//        printf("\n\n");
//        printf("temp: %lf %lf %lf\n", cpl_array_get_double(temp, 0, NULL),
//                                      cpl_array_get_double(temp, 1, NULL),
//                                      cpl_array_get_double(temp, 2, NULL));
//        printf("pres: %lf %lf %lf\n", cpl_array_get_double(pres, 0, NULL),
//                                      cpl_array_get_double(pres, 1, NULL),
//                                      cpl_array_get_double(pres, 2, NULL));
//        printf("rhum: %lf %lf %lf\n", cpl_array_get_double(rhum, 0, NULL),
//                                      cpl_array_get_double(rhum, 1, NULL),
//                                      cpl_array_get_double(rhum, 2, NULL));
//        printf("ppmv: %lf %lf %lf\n", cpl_array_get_double(ppmv, 0, NULL),
//                                      cpl_array_get_double(ppmv, 1, NULL),
//                                      cpl_array_get_double(ppmv, 2, NULL));
        cpl_msg_info("mf_basic_rhum2ppmv()", "Conversion of relative humidity"
                     " to ppmv successful");
    }

    profile = cpl_table_new(n);
    cpl_table_new_column(profile, "temp", CPL_TYPE_DOUBLE);
    cpl_table_new_column(profile, "pres", CPL_TYPE_DOUBLE);
    cpl_table_new_column(profile, "rhum", CPL_TYPE_DOUBLE);
    cpl_table_new_column(profile, "ppmv_old", CPL_TYPE_DOUBLE);
    cpl_table_new_column(profile, "ppmv_new", CPL_TYPE_DOUBLE);
    cpl_table_new_column(profile, "ppmv_ratio", CPL_TYPE_DOUBLE);
    cpl_table_new_column(profile, "ppmv_diff", CPL_TYPE_DOUBLE);

    cpl_array_set_size(temp, n);
    cpl_array_set_size(pres, n);
    cpl_array_set_size(rhum, n);
    cpl_array_set_size(ppmv, n);

    for (i = 0; i < n; i++) {
        cpl_array_set_double(temp, i, tlow + i);
        cpl_array_set_double(pres, i, 800);
        cpl_array_set_double(rhum, i, 20);

        cpl_table_set_double(profile, "temp", i, tlow+i);
        cpl_table_set_double(profile, "pres", i, 800);
        cpl_table_set_double(profile, "rhum", i, 20);
    }

    mf_basic_rhum2ppmv_old(temp, pres, rhum, ppmv);

    for (i = 0; i < n; i++) {
        d_temp = cpl_array_get_double(temp, i, NULL);
        d_pres = cpl_array_get_double(pres, i, NULL);
        d_rhum =  cpl_array_get_double(rhum, i, NULL);
        mf_basic_rhum2ppmv(&d_ppmv, &d_temp, &d_pres, &d_rhum);

        val = cpl_array_get_double(ppmv, i, NULL);
//        printf("%5i %15.10lf %15.10lf %15.10lf\n",
//               tlow+i, d_ppmv, val, (d_ppmv / val - 1) * 100);

        cpl_table_set_double(profile, "ppmv_old", i, val);
        cpl_table_set_double(profile, "ppmv_new", i, d_ppmv);
        cpl_table_set_double(profile, "ppmv_ratio", i, val/d_ppmv-1);
        cpl_table_set_double(profile, "ppmv_diff", i, (val-d_ppmv));
    }

    /* for IDL plot */
//    col = 0;
//    printf("temp = [");
//    for (i = 0, col = 1; i < n; i++, col++) {
//        printf("%5i,", tlow+i);
//        if (col == 10) {
//            col = 0;
//            printf(" $\n");
//        }
//    }
//    printf("\b]\nppmv_old = [");
//    for (i = 0, col = 1; i < n; i++, col++) {
//        printf("%10.5e,", cpl_table_get_double(profile, "ppmv_old", i, NULL));
//        if (col == 5) {
//            col = 0;
//            printf(" $\n");
//        }
//    }
//    printf("\b]\nppmv_new = [");
//    for (i = 0, col = 1; i < n; i++, col++) {
//        printf("%10.5e,", cpl_table_get_double(profile, "ppmv_new", i, NULL));
//        if (col == 5) {
//            col = 0;
//            printf(" $\n");
//        }
//    }
//    printf("\b]\n");

//    cpl_plot_column(NULL, NULL, NULL, profile, "temp", "ppmv_diff");

    cpl_table_delete(profile);

    cpl_array_delete(temp);
    cpl_array_delete(pres);
    cpl_array_delete(rhum);
    cpl_array_delete(ppmv);


    /************************************************************************
     * testing mf_basic_planck
     ************************************************************************/
    planck_wave = cpl_array_new(3, CPL_TYPE_DOUBLE);
    planck_res = cpl_array_new(3, CPL_TYPE_DOUBLE);

    cpl_array_set_double(planck_wave, 0, 0.4);
    cpl_array_set_double(planck_wave, 1, 3);
    cpl_array_set_double(planck_wave, 2, 20);

    if (mf_basic_planck(planck_res, planck_wave,
                        planck_temp) != CPL_ERROR_NONE) {
        cpl_msg_warning(cpl_func, "mf_basic_planck() encountered a problem.");
    } else {
//        printf("%25.20lf %25.20lf %25.20lf\n",
//               cpl_array_get_double(planck_res, 0, NULL),
//               planck_true[0],
//               cpl_array_get_double(planck_res, 0, NULL) - planck_true[0]);
//        printf("%25.20lf %25.20lf %25.20lf\n",
//               cpl_array_get_double(planck_res, 1, NULL),
//               planck_true[1],
//               cpl_array_get_double(planck_res, 1, NULL) - planck_true[1]);
//        printf("%25.20lf %25.20lf %25.20lf\n",
//               cpl_array_get_double(planck_res, 2, NULL),
//               planck_true[2],
//               cpl_array_get_double(planck_res, 2, NULL) - planck_true[2]);
        cpl_msg_info("mf_basic_planck()", "mf_basic_planck() SUCCESSFUL");
    }

    cpl_array_delete(planck_wave);
    cpl_array_delete(planck_res);


    /************************************************************************
     * testing mf_basic_greg2jd & mf_basic_jd2greg
     ************************************************************************/
    mf_basic_greg2jd(&jd, yr1, mn1, dy1);
//    printf("%li %li %li\n", jd, jd1, jd-jd1);
    cpl_test_rel(jd1, jd, 1e-7);
    mf_basic_jd2greg(&year, &month, &day, jd1);
//    printf("%i %i %i --- %i %i %i\n", yr1, mn1, dy1, year, month, day);
    cpl_test_rel(yr1, year, 1e-7);
    cpl_test_rel(mn1, month, 1e-7);
    cpl_test_rel(dy1, day, 1e-7);

    mf_basic_greg2jd(&jd, yr2, mn2, dy2);
//    printf("%li %li %li\n", jd, jd2, jd-jd2);
    cpl_test_rel(jd2, jd, 1e-7);
    mf_basic_jd2greg(&year, &month, &day, jd2);
//    printf("%i %i %i --- %i %i %i\n", yr2, mn2, dy2, year, month, day);
    cpl_test_rel(yr2, year, 1e-7);
    cpl_test_rel(mn2, month, 1e-7);
    cpl_test_rel(dy2, day, 1e-7);

    mf_basic_greg2jd(&jd, yr3, mn3, dy3);
//    printf("%li %li %li\n", jd, jd3, jd-jd3);
    cpl_test_rel(jd3, jd, 1e-7);
    mf_basic_jd2greg(&year, &month, &day, jd3);
//    printf("%i %i %i --- %i %i %i\n", yr3, mn3, dy3, year, month, day);
    cpl_test_rel(yr3, year, 1e-7);
    cpl_test_rel(mn3, month, 1e-7);
    cpl_test_rel(dy3, day, 1e-7);

    mf_basic_greg2jd(&jd, yr4, mn4, dy4);
//    printf("%li %li %li\n", jd, jd4, jd-jd4);
    cpl_test_rel(jd4, jd, 1e-7);
    mf_basic_jd2greg(&year, &month, &day, jd4);
//    printf("%i %i %i --- %i %i %i\n", yr4, mn4, dy4, year, month, day);
    cpl_test_rel(yr4, year, 1e-7);
    cpl_test_rel(mn4, month, 1e-7);
    cpl_test_rel(dy4, day, 1e-7);

    /************************************************************************
     * testing mf_basic_abspath
     ************************************************************************/
    if (getcwd(cwd_org, sizeof(cwd_org))) {}
    mf_basic_dirslash(cwd_org);
    strcpy(cwd, cwd_org);

    // testing "/dir1/./.././../dir2/../dir3/dir4/dir5/../.."
    mf_basic_abspath(outdir, dir0, cwd);
    sprintf(truedir, "/dir3/");
    cpl_test_eq_string(outdir, truedir);

    // testing ""
    mf_basic_abspath(outdir, dir1, cwd);
    sprintf(truedir, "%s", cwd);
    cpl_test_eq_string(outdir, truedir);

    // testing  "."
    mf_basic_abspath(outdir, dir2, cwd);
    sprintf(truedir, "%s", cwd);
    cpl_test_eq_string(outdir, truedir);

    // testing  "./"
    mf_basic_abspath(outdir, dir3, cwd);
    sprintf(truedir, "%s", cwd);
    cpl_test_eq_string(outdir, truedir);

    // testing  ".."
    mf_basic_abspath(outdir, dir4, cwd);
    if (chdir("..")) {}
    if (getcwd(cwd, sizeof(cwd))) {}
    sprintf(truedir, "%s/", cwd);
    cpl_test_eq_string(outdir, truedir);
    mf_basic_initstring(cwd, strlen(cwd));
    strcpy(cwd, cwd_org);

    // testing  "../"
    mf_basic_abspath(outdir, dir5, cwd);
    sprintf(truedir, "%s/", cwd_org);
    if (chdir(truedir)) {}
    if (chdir("..")) {}
    if (getcwd(cwd, sizeof(cwd))) {}
    sprintf(truedir, "%s/", cwd);
    cpl_test_eq_string(outdir, truedir);
    mf_basic_initstring(cwd, strlen(cwd));
    strcpy(cwd, cwd_org);

    // testing  "../.."
    mf_basic_abspath(outdir, dir6, cwd);
    sprintf(truedir, "%s/", cwd_org);
    if (chdir(truedir)) {}
    if (chdir("../..")) {}
    if (getcwd(cwd, sizeof(cwd))) {}
    sprintf(truedir, "%s/", cwd);
    cpl_test_eq_string(outdir, truedir);
    mf_basic_initstring(cwd, strlen(cwd));
    strcpy(cwd, cwd_org);

    // testing  "../../"
    mf_basic_abspath(outdir, dir7, cwd);
    sprintf(truedir, "%s/", cwd_org);
    if (chdir(truedir)) {}
    if (chdir("../..")) {}
    if (getcwd(cwd, sizeof(cwd))) {}
    sprintf(truedir, "%s/", cwd);
    cpl_test_eq_string(outdir, truedir);
    mf_basic_initstring(cwd, strlen(cwd));
    strcpy(cwd, cwd_org);

    // testing  "/"
    mf_basic_abspath(outdir, dir8, cwd);
    sprintf(truedir, "/");
    cpl_test_eq_string(outdir, truedir);

    // testing  "config"
    mf_basic_abspath(outdir, dir9, cwd);
    sprintf(truedir, "%sconfig/", cwd_org);
    cpl_test_eq_string(outdir, truedir);

    // testing  "config/"
    mf_basic_abspath(outdir, dir10, cwd);
    sprintf(truedir, "%sconfig/", cwd_org);
    cpl_test_eq_string(outdir, truedir);

    // testing  "../../config"
    mf_basic_abspath(outdir, dir11, cwd);
    sprintf(truedir, "%s/", cwd_org);
    if (chdir(truedir)) {}
    if (chdir("../..")) {}
    if (getcwd(cwd, sizeof(cwd))) {}
    sprintf(truedir, "%s/config/", cwd);
    cpl_test_eq_string(outdir, truedir);
    mf_basic_initstring(cwd, strlen(cwd));
    strcpy(cwd, cwd_org);

    // testing  "../../config/"
    mf_basic_abspath(outdir, dir12, cwd);
    sprintf(truedir, "%s/", cwd_org);
    if (chdir(truedir)) {}
    if (chdir("../..")) {}
    if (getcwd(cwd, sizeof(cwd))) {}
    sprintf(truedir, "%s/config/", cwd);
    cpl_test_eq_string(outdir, truedir);
    mf_basic_initstring(cwd, strlen(cwd));
    strcpy(cwd, cwd_org);

    // testing  "/config/"
    mf_basic_abspath(outdir, dir13, cwd);
    sprintf(truedir, "/config/");
    cpl_test_eq_string(outdir, truedir);

   /************************************************************************
     * testing mf_basic_dirslash
     ************************************************************************/
    mf_basic_dirslash(noslash);
    cpl_test_eq_string(noslash, "bla/");

    mf_basic_dirslash(slash);
    cpl_test_eq_string(slash, "bla/");

    /************************************************************************
     * testing mf_basic_absfile
     ************************************************************************/

    mf_basic_initstring(cwd, strlen(cwd));
    strcpy(cwd, cwd_org);

    // testing "/test/filename"
    mf_basic_absfile(outfile, file1);
    sprintf(truefile, "/test/filename");
    cpl_test_eq_string(outfile, truefile);

    // testing ".././filename"
    sprintf(cwd, "%s/", cwd_org);
    if (chdir(cwd)) {}
    mf_basic_absfile(outfile, file2);
    if (chdir("..")) {}
    if (getcwd(cwd, sizeof(cwd))) {}
    sprintf(truefile, "%s/filename", cwd);
    cpl_test_eq_string(outfile, truefile);

    // testing "/filename"
    mf_basic_absfile(outfile, file3);
    sprintf(truefile, "/filename");
    cpl_test_eq_string(outfile, truefile);

    // testing "filename"
    sprintf(cwd, "%s/", cwd_org);
    if (chdir(cwd)) {}
    mf_basic_absfile(outfile, file4);
    if (getcwd(cwd, sizeof(cwd))) {}
    sprintf(truefile, "%s/filename", cwd);
    cpl_test_eq_string(outfile, truefile);

    return cpl_test_end(0);
}
