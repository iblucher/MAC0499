#ifndef MF_INTERFACE_FITS_H
#define MF_INTERFACE_FITS_H

#include <cpl.h>

typedef struct {
    unsigned n;
    double* array;
} array_t;

typedef struct {
    unsigned n;
    char** array;
} tokenizer_t;

typedef tokenizer_t arraystr_t;

void julians_interface(const char*,
        cpl_parameterlist*,
        cpl_table*,
        cpl_table*,
        cpl_table*,
        cpl_table*,
        cpl_boolean);

tokenizer_t tokenize_string(const char*, const char*);

int find_param(const cpl_array*, const char*);

void process_range(const cpl_table*,
        const cpl_array*,
        const char*,
        cpl_parameterlist*,
        char*);

void process_molecule(const cpl_table*,
        const char*,
        arraystr_t*,
        array_t*,
        array_t*);

void process_array(const cpl_table*,
        const char*,
        array_t*,
        array_t*,
        array_t*,
        array_t*,
        array_t*,
        array_t*);

void process_scalar(const cpl_table*,
        const char*,
        cpl_parameterlist* plist);

void fill_array(array_t,
        array_t,
        cpl_table**);

void interface_fits(const char*,
        const char*,
        cpl_boolean);

#endif

