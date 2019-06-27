#include <molecfit.h>
#include <cpl.h>
#include <string.h>
#include "mf_interface_fits.h"

int main (int argc, char * argv[]) {

    cpl_init(CPL_INIT_DEFAULT);
    cpl_boolean cmp_mode = CPL_FALSE;
    /* use equal inputs to crires parameter file and store result
     * can be used to compare results of API usage vs parameter file usage */
    if (argc > 1 && strcmp(argv[1], "crires") == 0) {
        /* result should compare equal to molecfit_crires_tac.fits from
         * examples/config/molecfit_crires.par */
        cpl_msg_info(cpl_func, "CRIRES COMPARISON MODE");
        cmp_mode = CPL_TRUE;
    }

    const char *conf_fits = "examples/config/molecfit_crires.fits";
    const char *fn_fits   = "examples/input/crires_spec_jitter_extracted_0000.fits";
    
    interface_fits(conf_fits, fn_fits, cmp_mode);

    cpl_end();
    return  0 ;

}

