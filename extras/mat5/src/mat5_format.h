/* Element-level layer of the MAT-file reader; internal to extras/mat5. */
#ifndef ARC_MAT5_FORMAT_H
#define ARC_MAT5_FORMAT_H

#include <arclp/mat5.h>

#include <stdint.h>

typedef struct
{
    const uint8_t *p;
    const uint8_t *end;
} mat_cursor;

typedef struct
{
    uint32_t type;
    uint32_t bytes;
    const uint8_t *data;
} mat_element;

enum
{
    MI_INT8 = 1, MI_UINT8, MI_INT16, MI_UINT16, MI_INT32, MI_UINT32,
    MI_SINGLE, MI_DOUBLE = 9, MI_INT64 = 12, MI_UINT64, MI_MATRIX,
    MI_COMPRESSED
};

enum { MX_STRUCT = 2, MX_SPARSE = 5 };

/* Reads the next tagged element; returns 0 on success. */
int mat_next(mat_cursor *cur, mat_element *el);

/* Converts `count` numbers of any numeric MAT type to double. */
int mat_to_doubles(const mat_element *el, double *out, size_t count);

/* Number of entries a numeric element holds. */
size_t mat_count(const mat_element *el);

/* Parsed header of a miMATRIX element. */
typedef struct
{
    int mx_class;
    int rows, cols;
    mat_cursor body;      /* remaining sub-elements after the name */
} mat_array;

int mat_open_array(const mat_element *el, mat_array *arr);

/* Nonzero when el is a miMATRIX element whose variable name is `name`. */
int mat_array_named(const mat_element *el, const char *name);

/* For a 1x1 struct: finds the field `name`; returns 0 when found. */
int mat_struct_field(const mat_array *arr, const char *name, mat_array *out);

#endif
