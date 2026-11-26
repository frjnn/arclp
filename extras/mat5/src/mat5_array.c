/* miMATRIX headers and struct field lookup. */
#include "mat5_format.h"

#include <string.h>

static int read_int(const mat_element *el, size_t k, int *out)
{
    double v[2];
    if (k >= 2 || mat_to_doubles(el, v, k + 1) != 0)
    {
        return -1;
    }
    *out = (int)v[k];
    return 0;
}

int mat_open_array(const mat_element *el, mat_array *arr)
{
    memset(arr, 0, sizeof(*arr));
    if (el->type != MI_MATRIX)
    {
        return -1;
    }
    mat_cursor cur = {el->data, el->data + el->bytes};
    arr->body = cur;
    if (el->bytes == 0) /* empty array */
    {
        return 0;
    }
    mat_element flags, dims, name;
    if (mat_next(&cur, &flags) != 0 || flags.bytes < 8
        || mat_next(&cur, &dims) != 0 || mat_next(&cur, &name) != 0
        || read_int(&dims, 0, &arr->rows) != 0
        || read_int(&dims, 1, &arr->cols) != 0)
        {
        return -1;
    }
    arr->mx_class = flags.data[0];
    arr->body = cur;
    return 0;
}

int mat_array_named(const mat_element *el, const char *name)
{
    mat_cursor cur = {el->data, el->data + el->bytes};
    mat_element flags, dims, label;
    if (el->type != MI_MATRIX || mat_next(&cur, &flags) != 0
        || mat_next(&cur, &dims) != 0 || mat_next(&cur, &label) != 0)
        {
        return 0;
    }
    return strlen(name) == label.bytes
        && memcmp(label.data, name, label.bytes) == 0;
}

int mat_struct_field(const mat_array *arr, const char *name, mat_array *out)
{
    if (arr->mx_class != MX_STRUCT)
    {
        return -1;
    }
    mat_cursor cur = arr->body;
    mat_element len_el, names, field;
    int len = 0;
    if (mat_next(&cur, &len_el) != 0 || read_int(&len_el, 0, &len) != 0
        || len <= 0 || mat_next(&cur, &names) != 0)
        {
        return -1;
    }
    const size_t fields = names.bytes / (size_t)len;
    for (size_t f = 0; f < fields; f++)
    {
        if (mat_next(&cur, &field) != 0)
        {
            return -1;
        }
        const char *field_name = (const char *)names.data + f * (size_t)len;
        if (strncmp(field_name, name, (size_t)len) == 0)
        {
            return mat_open_array(&field, out);
        }
    }
    return -1;
}
