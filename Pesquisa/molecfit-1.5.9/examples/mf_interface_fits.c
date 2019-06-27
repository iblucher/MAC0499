
#define _DEFAULT_SOURCE
#define _BSD_SOURCE

#include <molecfit.h>

#include "mf_interface_fits.h"

#include <string.h>

// gcc -o cplread_fits cplread_fits.c -I$HOME/workspace/INTROOT/include -L$HOME/workspace/INTROOT/lib -lcplcore -lcplui -O0 -g

void julians_interface(
    const char        *fn_fits,
    cpl_parameterlist *parlist,
    cpl_table         *molectab,
    cpl_table         *wlinclude,
    cpl_table         *wlexclude,
    cpl_table         *pixexclude,
    cpl_boolean       cmp_mode)
{

    cpl_errorstate cleanstate = cpl_errorstate_get();

    /* state variable to transfer internal data from molecfit -> calctrans ->
     * calctrans -> ...  e.g. the working directory that can contains molecfit
     * results for be used in calctrans or to avoid recomputing LNFL data */
    mf_calctrans_state *state;

    /* load crires data, multiple chips spawning different wavelength in multi extension fits tables */
    printf("Loading data...\n");
    cpl_size nspec = 4;
    cpl_table * inspec[nspec];
    for (cpl_size i = 0; i < nspec; i++) {
        inspec[i] = cpl_table_load(fn_fits, i + 1, 0);
    }

    /* load the primary header of crires to obtain all of the TEL/INS keywords
     * needed by molecfit (lat/lon, temparature etc */
    cpl_propertylist *plist = cpl_propertylist_load(fn_fits, 0);

    /* Initialize variables for molecfit execution */
    cpl_table *prof_out = NULL;
    cpl_table *res_out  = NULL;
    cpl_table *spec_out = NULL;

    printf("Running molecfit...\n");
    mf_run_molecfit( parlist,
                     plist,
                     inspec,
                     nspec,
                     molectab,
                     wlinclude,
                     wlexclude,
                     pixexclude,
                     NULL,          /* kernel */
                     &prof_out,
                     &res_out,
                     &spec_out,
                     &state);

    if (!cpl_errorstate_is_equal(cleanstate)) {
        cpl_errorstate_dump(cleanstate, CPL_FALSE, NULL);
    }

    /* dummy kernel, one kernel per input spectra pixel */
    cpl_matrix * kernel = cpl_matrix_new(4096, 3);
    for (cpl_size r = 0; r < cpl_matrix_get_nrow(kernel); r++) {
        for (cpl_size c = 0; c < cpl_matrix_get_ncol(kernel); c++) {
            cpl_matrix_set(kernel, r, c, c);
        }
    }

    /* extended range test */
    double wl_start = 3274.33 - 1e2;
    double wl_end   = 3291.52 + 1e2;
    if (cmp_mode) {
        wl_start = -1;
        wl_end   = -1;
        cpl_matrix_delete(kernel);
        kernel = NULL;
    }

    cpl_table *res = NULL;
    mf_run_calctrans( parlist,
                      plist,
                      inspec,
                      nspec,
                      molectab,
                      kernel,
                      wl_start,
                      wl_end,
                      prof_out,
                      res_out,
                      &state,
                      &res);
    cpl_table_delete(res);

    cpl_msg_info(cpl_func, "second run");
    mf_run_calctrans( parlist,
                      plist,
                      inspec,
                      nspec,
                      molectab,
                      kernel,
                      wl_start,
                      wl_end,
                      prof_out,
                      res_out,
                      &state,
                      &res);
    if (cmp_mode) {
        cpl_table_save(res, NULL, NULL, "api-crires-output.fits", CPL_IO_CREATE);
    }
    cpl_table_delete(res);

    if (prof_out) cpl_table_delete(prof_out);
    if (res_out ) cpl_table_delete(res_out );
    if (spec_out) cpl_table_delete(spec_out);

    /* cleanup temporary data */
    mf_cleanup(state);
    
    if (!cpl_errorstate_is_equal(cleanstate)) {
        cpl_errorstate_dump(cleanstate, CPL_FALSE, NULL);
    }
}

// This function is used to tokenize a string with a specified delimiter
// It makes use of strsep which is thread safe.
tokenizer_t tokenize_string(const char* input, const char* delimiter) {

    tokenizer_t result;
    result.n = 0;
    result.array = NULL;

    char* input_duplicate = cpl_strdup(input);
    char* input_memory    = input_duplicate;

    char *token = strsep(&input_duplicate, delimiter);
    while (token) {

        char **tmp = NULL;
        if (result.n) {
            // Duplicate what we already have
            tmp = cpl_malloc(result.n * sizeof(char*));
            for (unsigned index = 0; index < result.n; ++index) {
                tmp[index] = cpl_strdup(result.array[index]);
            }
            // Free the current array
            for (unsigned index = 0; index < result.n; ++index) {
                cpl_free(result.array[index]);
            }
            cpl_free(result.array);
        }
        // Allocate the new array
        result.array = cpl_malloc((result.n + 1) * sizeof(char*));
        // Copy back what we had before
        for (unsigned index = 0; index < result.n; ++index) {
            result.array[index] = cpl_strdup(tmp[index]);
        }
        // Free the backup
        for (unsigned index = 0; index < result.n; ++index) {
            cpl_free(tmp[index]);
        }
        if (tmp) {
            cpl_free(tmp);
        }
        // Add the new token
        result.array[result.n] = cpl_strdup(token);
        // Increment the counter
        ++(result.n);

        token = strsep(&input_duplicate, delimiter);
    }
    cpl_free(input_memory);

    return result;

}

// Find parameter name in the input array
int find_param(const cpl_array* colnames, const char* name) {

    for (unsigned index = 0; index < cpl_array_get_size(colnames); ++index) {
        if (!strcmp(cpl_array_get_string(colnames, index), name)) {
            return index;
        }
    }

    return -1;

}

// Process a parameter range
void process_range(const cpl_table* parameter_tab,
        const cpl_array* colnames,
        const char* name,
        cpl_parameterlist* plist,
        char* todo) {

    // Column type
    cpl_type coltype = cpl_table_get_column_type(parameter_tab, name);
    const char* strcoltype = cpl_type_get_name(coltype);
    // Duplicate the string
    unsigned idx;
    // Get the parameter name
    const char delim = '-';
    unsigned ids = 6;
    for (idx = ids; idx < strlen(name); ++idx) {
        if (name[idx] == delim) {
            break;
        }
    }
    unsigned nchar = idx - ids;
    char* parname = malloc((nchar + 1) * sizeof(char));
    strncpy(parname, &(name[ids]), nchar);
    parname[nchar] = '\0';
    printf("parname= %s\n", parname);
    // Find the other corresponding parameter
    char defname[100];
    char staname[100];
    char endname[100];
    sprintf(defname, "range-%s-default", parname);
    sprintf(staname, "range-%s-start", parname);
    sprintf(endname, "range-%s-end", parname);
    printf("names: %s\t%s\t%s\n", defname, staname, endname);
    int indexdef = find_param(colnames, defname);
    int indexsta = find_param(colnames, staname);
    int indexend = find_param(colnames, endname);
    printf("indexes: %d\t%d\t%d\n", indexdef, indexsta, indexend);
    todo[indexdef] = 0;
    todo[indexsta] = 0;
    todo[indexend] = 0;
    // Add parameter to the list
    cpl_parameter* p = NULL;
    if (coltype == CPL_TYPE_INT) {
        int null = 0;
        int def = cpl_table_get_int(parameter_tab, defname, 0, &null);
        int start = cpl_table_get_int(parameter_tab, staname, 0, &null);
        int end = cpl_table_get_int(parameter_tab, endname, 0, &null);
        p = cpl_parameter_new_range(parname, coltype, "", "", def, start, end);
    }
    else if (coltype == CPL_TYPE_DOUBLE) {
        int null = 0;
        double def = cpl_table_get_double(parameter_tab, defname, 0, &null);
        double start = cpl_table_get_double(parameter_tab, staname, 0, &null);
        double end = cpl_table_get_double(parameter_tab, endname, 0, &null);
        p = cpl_parameter_new_range(parname, coltype, "", "", def, start, end);
    }
    else {
        printf("ERROR: type %s not supported.", strcoltype);
    }
    if (p) {
        cpl_parameterlist_append(plist, p);
    }
    free(parname);

}

// Process a parameter which must be added to the molecule table
void process_molecule(const cpl_table* parameter_tab,
        const char* name,
        arraystr_t* molec_names,
        array_t* molec_fit,
        array_t* molec_relcol) {

    const char* value = cpl_table_get_string(parameter_tab, name, 0);
    tokenizer_t tokens = tokenize_string(value, " ");
    if (!strncmp(&(name[6]), "list_molec", 10)) {  // list_molec array
        molec_names->n = tokens.n;
        molec_names->array = malloc(tokens.n * sizeof(char*));
        for (unsigned idxtoken = 0; idxtoken < tokens.n; ++idxtoken) {
            printf("list_molec, appending %s\n", tokens.array[idxtoken]);
            molec_names->array[idxtoken] = cpl_strdup(tokens.array[idxtoken]);
        }
    }
    else if (!strncmp(&(name[6]), "fit_molec", 9)) {  // fit_molec array)
        molec_fit->n = tokens.n;
        molec_fit->array = malloc(tokens.n * sizeof(double));
        for (unsigned idxtoken = 0; idxtoken < tokens.n; ++idxtoken) {
            printf("fit_molec, appending %s\n", tokens.array[idxtoken]);
            molec_fit->array[idxtoken] = atof(tokens.array[idxtoken]);
        }
    }
    else {  // relcol array
        molec_relcol->n = tokens.n;
        molec_relcol->array = malloc(tokens.n * sizeof(double));
        for (unsigned idxtoken = 0; idxtoken < tokens.n; ++idxtoken) {
            printf("relcol, appending %s\n", tokens.array[idxtoken]);
            molec_relcol->array[idxtoken] = atof(tokens.array[idxtoken]);
        }
    }
    for (unsigned idxtoken = 0; idxtoken < tokens.n; ++idxtoken) {
        free(tokens.array[idxtoken]);
    }
    free(tokens.array);

}

void fill_molecules(arraystr_t molec_names,
        array_t molec_fit,
        array_t molec_relcol,
        cpl_table** molectab) {

    if (molec_names.n || molec_fit.n || molec_relcol.n) {
        if (molec_names.n == molec_fit.n && molec_names.n == molec_relcol.n) {
            *molectab = cpl_table_new(molec_names.n);
            cpl_table_new_column(*molectab, "list_molec", CPL_TYPE_STRING);
            cpl_table_new_column(*molectab, "fit_molec", CPL_TYPE_INT);
            cpl_table_new_column(*molectab, "relcol", CPL_TYPE_DOUBLE);
            for (unsigned index = 0; index < molec_names.n; ++index) {
                cpl_table_set_string(*molectab, "list_molec", index, molec_names.array[index]);
                cpl_table_set_int(*molectab, "fit_molec", index, molec_fit.array[index]);
                cpl_table_set_double(*molectab, "relcol", index, molec_relcol.array[index]);
            }    
        }
        else {
            printf("Error: arrays of molecules do not have the same size.");
        }
        if (molec_names.n) {
            for (unsigned index = 0; index < molec_names.n; ++index) {
                free(molec_names.array[index]);
            }
            free(molec_names.array);
        }
        if (molec_fit.n) {
            free(molec_fit.array);
        }
        if (molec_relcol.n) {
            free(molec_relcol.array);
        }
    }

}

// Process a parameter which must be added to an array
void process_array(const cpl_table* parameter_tab,
        const char* name,
        array_t* wlinclude_l,
        array_t* wlinclude_u,
        array_t* wlexclude_l,
        array_t* wlexclude_u,
        array_t* pixexclude_l,
        array_t* pixexclude_u) {

    const char* value = cpl_table_get_string(parameter_tab, name, 0);
    if (!strncmp(value, "none", 4)) {
        printf("contains nothing...\n");
        return;
    }
    tokenizer_t tokens = tokenize_string(value, " ");
    double* current_pointer = malloc(tokens.n * sizeof(double));
    printf("Values: ");
    for (unsigned idxtoken = 0; idxtoken < tokens.n; ++idxtoken) {
        double v = atof(tokens.array[idxtoken]);
        current_pointer[idxtoken] = v;
        printf("\t%lf", v);
    }
    for (unsigned idxtoken = 0; idxtoken < tokens.n; ++idxtoken) {
        free(tokens.array[idxtoken]);
    }
    free(tokens.array);
    if (!strncmp(&(name[7]), "include", 7)) {  // This is wrange_include
        if (!strncmp(&(name[15]), "llim", 4))  {  // This is wrange_include_llim
            printf("\tfor wrange_include_llim\n");
            wlinclude_l->n = tokens.n;
            wlinclude_l->array = current_pointer;
        }
        else {  // This is wrange_include_ulim
            printf("\tfor wrange_include_ulim\n");
            wlinclude_u->n = tokens.n;
            wlinclude_u->array = current_pointer;
        }
    }
    else if (!strncmp(name, "p", 1)) {  // This is prange_exclude
        if (!strncmp(&(name[15]), "llim", 4))  {  // This is prange_exclude_llim
            printf("\tfor prange_exclude_llim\n");
            pixexclude_l->n = tokens.n;
            pixexclude_l->array = current_pointer;
        }
        else {  // This is prange_exclude_ulim
            printf("\tfor prange_exclude_ulim\n");
            pixexclude_u->n = tokens.n;
            pixexclude_u->array = current_pointer;
        }
    }
    else {  // This is wrange_exclude
        if (!strncmp(&(name[15]), "llim", 4))  {  // This is wrange_exclude_llim
            printf("\tfor wrange_exclude_llim\n");
            wlexclude_l->n = tokens.n;
            wlexclude_l->array = current_pointer;
        }
        else {  // This is wrange_exclude_ulim
            printf("\tfor wrange_exclude_ulim\n");
            wlexclude_u->n = tokens.n;
            wlexclude_u->array = current_pointer;
        }
    }

}

// Process a scalar
void process_scalar(const cpl_table* parameter_tab,
        const char* name,
        cpl_parameterlist* plist) {

    // Column type
    cpl_type coltype = cpl_table_get_column_type(parameter_tab, name);
    const char* strcoltype = cpl_type_get_name(coltype);

    cpl_parameter* p = NULL;
    if (coltype == CPL_TYPE_INT) {
        int null = 0;
        int value = cpl_table_get_int(parameter_tab, name, 0, &null);
        p = cpl_parameter_new_value(name, coltype, "", "", value);
    }
    else if (coltype == CPL_TYPE_DOUBLE) {
        int null = 0;
        double value = cpl_table_get_double(parameter_tab, name, 0, &null);
        p = cpl_parameter_new_value(name, coltype, "", "", value);
    }
    else if (coltype == CPL_TYPE_STRING) {
        const char* value = cpl_table_get_string(parameter_tab, name, 0);
        p = cpl_parameter_new_value(name, coltype, "", "", value);
    }
    else {
        printf("ERROR: type %s not supported.", strcoltype);
    }
    if (p) {
        cpl_parameterlist_append(plist, p);
    }

}

void fill_array(array_t array_l,
        array_t array_u,
        cpl_table** table) {

    if (array_l.n || array_u.n) {
        if (array_l.n == array_u.n) {
            *table = cpl_table_new(array_u.n);
            cpl_table_new_column(*table, "llim", CPL_TYPE_DOUBLE);
            cpl_table_new_column(*table, "ulim", CPL_TYPE_DOUBLE);
            // Fill in the cpl_table
            for (unsigned index = 0; index < array_u.n; ++index) {
                cpl_table_set(*table, "llim", index, array_l.array[index]);
                cpl_table_set(*table, "ulim", index, array_u.array[index]);
            }
        }
        else {
            printf("Error: array of lower limits and array of upper limits do not have the same size.");
        }
        free(array_l.array);
        free(array_u.array);
    }

}

void interface_fits(
    const char* conf_filepath,
    const char* spec_filepath,
    cpl_boolean cmp_mode)
{

    cpl_table* parameter_tab = cpl_table_load(conf_filepath, 1, 0);

    unsigned ncol = cpl_table_get_ncol(parameter_tab);
    printf("Found %u columns\n", ncol);

    char* todo = malloc(ncol * sizeof(char));
    for (unsigned index = 0; index < ncol; ++index) {
        todo[index] = 1;
    }

    cpl_parameterlist* plist = cpl_parameterlist_new();

    cpl_table* molectab = NULL;

    cpl_array* colnames = cpl_table_get_column_names(parameter_tab);

    array_t wlinclude_l;
    array_t wlinclude_u;
    array_t wlexclude_l;
    array_t wlexclude_u;
    array_t pixexclude_l;
    array_t pixexclude_u;

    wlinclude_l.n = 0;
    wlinclude_l.array = NULL;
    wlinclude_u.n = 0;
    wlinclude_u.array = NULL;

    wlexclude_l.n = 0;
    wlexclude_l.array = NULL;
    wlexclude_u.n = 0;
    wlexclude_u.array = NULL;

    pixexclude_l.n = 0;
    pixexclude_l.array = NULL;
    pixexclude_u.n = 0;
    pixexclude_u.array = NULL;

    arraystr_t molec_names;
    array_t molec_fit;
    array_t molec_relcol;

    molec_names.n = 0;
    molec_names.array = NULL;

    molec_fit.n = 0;
    molec_fit.array = NULL;

    molec_relcol.n = 0;
    molec_relcol.array = NULL;

    for (unsigned index = 0; index < ncol; ++index) {

        const char* name = cpl_array_get_string(colnames, index);

        cpl_type coltype = cpl_table_get_column_type(parameter_tab, name);
        const char* strcoltype = cpl_type_get_name(coltype);

        printf("%u\t%s\t%s\n", index, name, strcoltype);

        int parindex = find_param(colnames, name);

        if (!todo[parindex]) {
            continue;
        }

        if (!strncmp(name, "range", 5)) {
            printf("Found a range\n");
            process_range(parameter_tab, colnames, name, plist, todo);
        }
        else if (!strncmp(name, "molec", 5)) {  // Molecules
            printf("Found a molecule element\n");
            process_molecule(parameter_tab, name,
                    &molec_names, &molec_fit, &molec_relcol);
            // process_molecule(parameter_tab, name, molectab);
        }
        else if (!strncmp(&(name[1]), "range_", 6) &&
                !(strncmp(name, "w", 1) && strncmp(name, "p", 1))) {  // wrange_include/exclude, prange_exclude
            printf("Found an include/exclude range: %s\n", name);
            process_array( parameter_tab,
                           name,
                           &wlinclude_l,
                           &wlinclude_u,
                           &wlexclude_l,
                           &wlexclude_u,
                           &pixexclude_l,
                           &pixexclude_u);
        }
        else{
            printf("Found a scalar\n");
            process_scalar(parameter_tab, name, plist);
        }

    }

    free(todo);

    cpl_table* wlinclude = NULL;
    cpl_table* wlexclude = NULL;
    cpl_table* pixexclude = NULL;

    printf("Filling molecules...\n");
    fill_molecules(molec_names, molec_fit, molec_relcol, &molectab);

    printf("Filling arrays...\n");
    fill_array(wlinclude_l,  wlinclude_u,  &wlinclude);
    fill_array(wlexclude_l,  wlexclude_u,  &wlexclude);
    fill_array(pixexclude_l, pixexclude_u, &pixexclude);

    printf("Julian's interface now...\n");
    julians_interface( spec_filepath,
                       plist,
                       molectab,
                       wlinclude,
                       wlexclude,
                       pixexclude,
                       cmp_mode);

    /* Cleanup */
    if (wlinclude ) cpl_table_delete(wlinclude );
    if (wlexclude ) cpl_table_delete(wlexclude );
    if (pixexclude) cpl_table_delete(pixexclude);
    if (molectab  ) cpl_table_delete(molectab  );

    cpl_array_delete(colnames);
    cpl_table_delete(parameter_tab);
    cpl_parameterlist_delete(plist);
}

