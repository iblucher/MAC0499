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

    cpl_table * pixexclude = fillExcludedPixels();

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
    cpl_msg_info(cpl_func, "Saving ATM table: %s, status = %d\n", atmfitsname, status);

    const char *resfitsname = "testRES.fits";
    status = cpl_table_save(res_table, NULL, NULL, resfitsname, CPL_IO_CREATE);
    cpl_msg_info(cpl_func, "Saving RES table: %s, status = %d\n", resfitsname, status);

    const char *specfitsname = "testFIT.fits";
    status = cpl_table_save(spec_table, NULL, NULL, specfitsname, CPL_IO_CREATE);
    cpl_msg_info(cpl_func, "Saving FIT table: %s, status = %d\n", specfitsname, status);

    cpl_table_delete(spec_table);

    if (!cpl_errorstate_is_equal(cleanstate))
        cpl_errorstate_dump(cleanstate, CPL_FALSE, NULL);

    cpl_table_delete(atmprof);
    cpl_table_delete(res_table);

    /* cleanup temporary data */
    mf_cleanup(state);

    if (molectab) {
        cpl_table_delete(molectab);
    }
    if (pixexclude) {
        cpl_table_delete(pixexclude);
    }

    if (parlist) {
        cpl_parameterlist_delete(parlist);
    }
    if (plist) {
        cpl_propertylist_delete(plist);
    }

    for (cpl_size i = 0; i < nspec; i++) {
        cpl_table_delete(spec[i]);
    }

    if (!cpl_errorstate_is_equal(cleanstate))
        cpl_errorstate_dump(cleanstate, CPL_FALSE, NULL);
    cpl_end();

    return  0 ;
}
