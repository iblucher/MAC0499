#include <molecfit.h>
#include <cpl.h>
#include <string.h>
#include "api-input.h"

int main(int argc, char * argv[]) {

    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate cleanstate = cpl_errorstate_get();
    char datapath[1000];

    if (argc > 1) {
        strcpy(datapath, argv[1]);
    }
    else {
        strcpy(datapath, ".");
    }

    cpl_parameterlist *parlist = fillParList();

    cpl_size nspec = 4;
    cpl_table * spec[nspec];
    strcat(datapath, "/crires_spec_jitter_extracted_0000.fits");
    cpl_msg_info(cpl_func, "Reading input spectra: %s", datapath);
    loadData(spec, nspec, datapath);

    /* load the primary header of crires to obtain all of the TEL/INS keywords
     * needed by molecfit (lat/lon, temparature etc */
    cpl_propertylist * plist = loadPropList(datapath);

    cpl_table * molectab = fillMolecules();

    mf_calctrans_state * state = NULL;
    state = mf_init_calctrans_state();

    cpl_table* atmprof = NULL;
    cpl_table* res_table = NULL;
    cpl_table* spec_table = NULL;

    // Load the required results from disk.

    const char *atmfitsname = "testATM.fits";
    atmprof = cpl_table_load(atmfitsname, 1, 0);

    const char *resfitsname = "testRES.fits";
    res_table = cpl_table_load(resfitsname, 1, 0);

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

    // Can put -1 for extended range
    double wl_start = 3274.33 - 1e2;
    double wl_end = 3291.52 + 1e2;

    cpl_table * res = NULL;
    cpl_msg_info(cpl_func, "Running mf_run_calctrans_lblrtm");
    mf_run_calctrans_lblrtm(parlist, plist, spec, nspec, molectab,
                            wl_start, wl_end,
                            atmprof, res_table,
                            &state);
    cpl_msg_info(cpl_func, "Running mf_run_calctrans_convolution");
    mf_run_calctrans_convolution(parlist, plist, spec, nspec,
                                 kernel,
                                 wl_start, wl_end,
                                 atmprof, res_table,
                                 &state,
                                 &res);

    cpl_table_delete(atmprof);
    cpl_table_delete(res_table);

    const char *resultfitsname = "testRESULTS.fits";
    cpl_msg_info(cpl_func, "Saving RESULTS table: %s\n", resultfitsname);
    cpl_table_save(res, NULL, NULL, resultfitsname, CPL_IO_CREATE);

    if (res) {
        cpl_table_delete(res);
    }

    /* cleanup temporary data */
    mf_cleanup(state);

    if (molectab) {
        cpl_table_delete(molectab);
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
