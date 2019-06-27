#include "api-input.h"

cpl_parameterlist* fillParList(void) {

    cpl_parameterlist * parlist = cpl_parameterlist_new();

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

    return parlist;

}

void loadData(cpl_table** spec, cpl_size nspec, char* datapath) {

    for (cpl_size i = 0; i < nspec; i++) {
        spec[i] = cpl_table_load(datapath, i + 1, 0);
    }

}

cpl_propertylist* loadPropList(char* datapath) {

    return cpl_propertylist_load(datapath, 0);

}

cpl_table * fillMolecules(void) {

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

    return molectab;

}

cpl_table* fillExcludedPixels(void) {

    cpl_table* pixexclude = cpl_table_new(5);
    cpl_table_new_column(pixexclude, "llim", CPL_TYPE_DOUBLE);
    cpl_table_new_column(pixexclude, "ulim", CPL_TYPE_DOUBLE);
    cpl_table_set(       pixexclude, "llim", 0, 1);
    cpl_table_set(       pixexclude, "ulim", 0, 20);
    cpl_table_set(       pixexclude, "llim", 1, 1005);
    cpl_table_set(       pixexclude, "ulim", 1, 1044);
    cpl_table_set(       pixexclude, "llim", 2, 2029);
    cpl_table_set(       pixexclude, "ulim", 2, 2068);
    cpl_table_set(       pixexclude, "llim", 3, 3053);
    cpl_table_set(       pixexclude, "ulim", 3, 3112);
    cpl_table_set(       pixexclude, "llim", 4, 4057);
    cpl_table_set(       pixexclude, "ulim", 4, 4096);

    return pixexclude;

}
