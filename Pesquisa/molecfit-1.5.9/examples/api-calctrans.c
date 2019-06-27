#include <molecfit.h>
#include <cpl.h>
#include <string.h>


int main(int argc, char * argv[])
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate cleanstate = cpl_errorstate_get();
    cpl_boolean cmp_mode = CPL_FALSE;
    /* use equal inputs to crires parameter file and store result
     * can be used to compare results of API usage vs parameter file usage */
    if (argc > 1 && strcmp(argv[1], "crires") == 0) {
        /* result should compare equal to molecfit_crires_tac.fits from
         * examples/config/molecfit_crires.par */
        cpl_msg_info(cpl_func, "CRIRES COMPARISON MODE");
        cmp_mode = CPL_TRUE;
    }

    cpl_parameterlist * parlist = cpl_parameterlist_new();


    /* load crires data, multiple chips spawning different wavelength in multi
     * extension fits tables */
    cpl_size nspec = 4;
    cpl_table * spec[nspec];
    const char * fn_fits = "examples/input/crires_spec_jitter_extracted_0000.fits";
    for (cpl_size i = 0; i < nspec; i++) {
        spec[i] = cpl_table_load(fn_fits, i + 1, 0);
    }

    /* load the primary header of crires to obtain all of the TEL/INS keywords
     * needed by molecfit (lat/lon, temparature etc */
    cpl_propertylist * plist = cpl_propertylist_load(fn_fits, 0);


    /* set column meanings in crires table data */
    cpl_parameter * p;
    p = cpl_parameter_new_value("col_lam", CPL_TYPE_STRING, "", "recipename", "Wavelength");
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("col_flux", CPL_TYPE_STRING, "", "recipename", "Extracted_OPT");
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("col_dflux", CPL_TYPE_STRING, "", "recipename", "Error_OPT");
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("col_mask", CPL_TYPE_STRING, "", "recipename", "NULL");
    cpl_parameterlist_append(parlist, p);

    /* update parameters (see examples/config/molecfit_crires.par) */
    p = cpl_parameter_new_value("ftol", CPL_TYPE_DOUBLE, "", "recipename", 0.01);
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("xtol", CPL_TYPE_DOUBLE, "", "recipename", 0.01);
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_range("fit_back", CPL_TYPE_INT, "", "", 0, 0, 1);
    cpl_parameterlist_append(parlist, p);

    p = cpl_parameter_new_value("cont_n", CPL_TYPE_INT, "", "", 3);
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("wlc_n", CPL_TYPE_INT, "", "", 3);
    cpl_parameterlist_append(parlist, p);

    p = cpl_parameter_new_range("fit_res_lorentz", CPL_TYPE_INT, "", "", 0, 0,
                                1);
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("res_lorentz", CPL_TYPE_DOUBLE, "", "", .5);
    cpl_parameterlist_append(parlist, p);

    p = cpl_parameter_new_range("fit_res_box", CPL_TYPE_INT, "", "", 0, 0,
                                1);
    cpl_parameterlist_append(parlist, p);
    p = cpl_parameter_new_value("relres_box", CPL_TYPE_DOUBLE, "", "", 0.);
    cpl_parameterlist_append(parlist, p);

    p = cpl_parameter_new_value("kernfac", CPL_TYPE_DOUBLE, "", "", 300.);
    cpl_parameterlist_append(parlist, p);

    p = cpl_parameter_new_value("wlgtomicron", CPL_TYPE_DOUBLE, "", "recipename", 0.001);
    cpl_parameterlist_append(parlist, p);

    /* define molecules to fit */
    cpl_table * molectab = cpl_table_new(3);
    cpl_table_new_column(molectab, "list_molec", CPL_TYPE_STRING);
    cpl_table_new_column(molectab, "fit_molec",  CPL_TYPE_INT);
    cpl_table_new_column(molectab, "relcol",     CPL_TYPE_DOUBLE);
    cpl_table_set_string(molectab, "list_molec", 0, "H2O");
    cpl_table_set_string(molectab, "list_molec", 1, "O3");
    cpl_table_set_string(molectab, "list_molec", 2, "CH4");
    cpl_table_set_int(   molectab, "fit_molec",  0, 1);
    cpl_table_set_int(   molectab, "fit_molec",  1, 1);
    cpl_table_set_int(   molectab, "fit_molec",  2, 1);
    cpl_table_set_double(molectab, "relcol",     0, 1.);
    cpl_table_set_double(molectab, "relcol",     1, 1.);
    cpl_table_set_double(molectab, "relcol",     2, 1.);

    /* Setup inclusion and exclusion values
     * TODO merge excludes into single table ?
     * pixeclude could be done via invalid in input spec (wlinclude could also
     * be via selection) */

    /* not used in crires */
    cpl_table * wlexclude = cpl_table_new(1);
    cpl_table_new_column(wlexclude, "llim", CPL_TYPE_DOUBLE);
    cpl_table_new_column(wlexclude, "ulim", CPL_TYPE_DOUBLE);

    /* not used in crires */
    cpl_table * wlinclude = cpl_table_new(1);
    cpl_table_new_column(wlinclude, "llim", CPL_TYPE_DOUBLE);
    cpl_table_new_column(wlinclude, "ulim", CPL_TYPE_DOUBLE);

    cpl_table * pixexclude = cpl_table_new(5);
    cpl_table_new_column(pixexclude, "llim", CPL_TYPE_DOUBLE);
    cpl_table_new_column(pixexclude, "ulim", CPL_TYPE_DOUBLE);
    cpl_table_set(pixexclude, "llim", 0, 1);
    cpl_table_set(pixexclude, "ulim", 0, 20);
    cpl_table_set(pixexclude, "llim", 1, 1005);
    cpl_table_set(pixexclude, "ulim", 1, 1044);
    cpl_table_set(pixexclude, "llim", 2, 2029);
    cpl_table_set(pixexclude, "ulim", 2, 2068);
    cpl_table_set(pixexclude, "llim", 3, 3053);
    cpl_table_set(pixexclude, "ulim", 3, 3112);
    cpl_table_set(pixexclude, "llim", 4, 4057);
    cpl_table_set(pixexclude, "ulim", 4, 4096);

    /* state variable to transfer internal data from molecfit -> calctrans ->
     * calctrans -> ...  e.g. the working directory that can contains molecfit
     * results for be used in calctrans or to avoid recomputing LNFL data */

    mf_calctrans_state * state = NULL;
    state = mf_init_calctrans_state();

    cpl_table* atmprof = NULL;
    cpl_table* res_table = NULL;
    cpl_table* spec_table = NULL;

    mf_run_molecfit(parlist, plist, spec, nspec, molectab,
                 /* wlinclude */ NULL,
                 /* wlexclude */ NULL,
                 pixexclude,
                 /*kernel */ NULL,
                 &atmprof,
                 &res_table,
                 &spec_table,
                 &state);

    const char *atmfitsname = "testATM.fits";
    cpl_error_code status = cpl_table_save(atmprof, NULL, NULL, atmfitsname, CPL_IO_CREATE);
    printf("Saved ATM table: %s, status = %d\n", atmfitsname, status);

    const char *resfitsname = "testRES.fits";
    status = cpl_table_save(res_table, NULL, NULL, resfitsname, CPL_IO_CREATE);
    printf("Saved RES table: %s, status = %d\n", resfitsname, status);

    const char *specfitsname = "testFIT.fits";
    status = cpl_table_save(spec_table, NULL, NULL, specfitsname, CPL_IO_CREATE);
    printf("Saved FIT table: %s, status = %d\n", specfitsname, status);

    cpl_table_delete(spec_table);

    if (!cpl_errorstate_is_equal(cleanstate))
        cpl_errorstate_dump(cleanstate, CPL_FALSE, NULL);


    /* dummy kernel, one kernel per input spectra pixel */
    cpl_matrix * kernel = cpl_matrix_new(4096, 3);
    for (cpl_size r = 0; r < cpl_matrix_get_nrow(kernel); r++) {
        for (cpl_size c = 0; c < cpl_matrix_get_ncol(kernel); c++) {
            cpl_matrix_set(kernel, r, c, c);
        }
    }

    /* extended range test */
    double wl_start = 3274.33 - 1e2;
    double wl_end = 3291.52 + 1e2;
    if (cmp_mode) {
        wl_start = -1;
        wl_end = -1;
        cpl_matrix_delete(kernel);
        kernel = NULL;
    }

    cpl_table *res = NULL;
    mf_run_calctrans_lblrtm(      parlist,
                                  plist,
                                  spec,
                                  nspec,
                                  molectab,
                                  wl_start,
                                  wl_end,
                                  atmprof,
                                  res_table,
                                  &state);

    mf_run_calctrans_convolution( parlist,
                                  plist,
                                  spec,
                                  nspec,
                                  kernel,
                                  wl_start,
                                  wl_end,
                                  atmprof,
                                  res_table,
                                  &state,
                                  &res);
    if (res) cpl_table_delete(res);


    cpl_msg_info(cpl_func, "second run");

    mf_run_calctrans_lblrtm(      parlist,
                                  plist,
                                  spec,
                                  nspec,
                                  molectab,
                                  wl_start,
                                  wl_end,
                                  atmprof,
                                  res_table,
                                  &state);

    mf_run_calctrans_convolution( parlist,
                                  plist,
                                  spec,
                                  nspec,
                                  kernel,
                                  wl_start,
                                  wl_end,
                                  atmprof,
                                  res_table,
                                  &state,
                                  &res);

    cpl_table_delete(atmprof);
    cpl_table_delete(res_table);

    if (cmp_mode) {
        cpl_table_save(res, NULL, NULL, "api-crires-output.fits", CPL_IO_CREATE);
    }

    if (res) {
        cpl_table_delete(res);
    }

    /* cleanup temporary data */
    mf_cleanup(state);


    if (molectab) {
        cpl_table_delete(molectab);
    }
    if (wlexclude) {
        cpl_table_delete(wlexclude);
    }
    if (pixexclude) {
        cpl_table_delete(pixexclude);
    }
    if (wlinclude) {
        cpl_table_delete(wlinclude);
    }

    if (parlist) {
        cpl_parameterlist_delete(parlist);
    }
    if (plist) {
        cpl_propertylist_delete(plist);
    }
    if (kernel) {
        cpl_matrix_delete(kernel);
    }
    for (cpl_size i = 0; i < nspec; i++) {
        cpl_table_delete(spec[i]);
    }

    if (!cpl_errorstate_is_equal(cleanstate))
        cpl_errorstate_dump(cleanstate, CPL_FALSE, NULL);
    cpl_end();

    return  0 ;
}
