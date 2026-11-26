/* Extraction of an LP (A, b, c) from a SuiteSparse "Problem" struct. */
#include "mat5_format.h"

#include "mat5_file.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fail(char *err, size_t len, const char *what)
{
    snprintf(err, len, "%s", what);
    return -1;
}

static int read_vector(const mat_array *problem, const char *path0,
                       const char *path1, double **out, int *length)
{
    mat_array arr, inner;
    if (mat_struct_field(problem, path0, &arr) != 0)
    {
        return -1;
    }
    if (path1 != NULL)
    {
        if (mat_struct_field(&arr, path1, &inner) != 0)
        {
            return -1;
        }
        arr = inner;
    }
    mat_element data;
    *length = arr.rows * arr.cols;
    *out = malloc(sizeof(double) * (size_t)(*length + 1));
    if (*out == NULL || mat_next(&arr.body, &data) != 0)
    {
        return -1;
    }
    return mat_to_doubles(&data, *out, (size_t)*length);
}

static int to_ints(const mat_element *el, int *out, size_t count)
{
    double *tmp = malloc(sizeof(double) * (count + 1));
    const int status = tmp == NULL ? -1 : mat_to_doubles(el, tmp, count);
    for (size_t k = 0; status == 0 && k < count; k++)
    {
        out[k] = (int)tmp[k];
    }
    free(tmp);
    return status;
}

/* Sparse arrays store row indices (ir), column pointers (jc), values (pr). */
static int read_sparse(const mat_array *problem, mat_lp *lp)
{
    mat_array a;
    mat_element ir, jc, pr;
    if (mat_struct_field(problem, "A", &a) != 0 || a.mx_class != MX_SPARSE
        || mat_next(&a.body, &ir) != 0 || mat_next(&a.body, &jc) != 0
        || mat_next(&a.body, &pr) != 0)
        {
        return -1;
    }
    lp->m = a.rows;
    lp->n = a.cols;
    lp->col_ptr = malloc(sizeof(int) * (size_t)(lp->n + 1));
    if (lp->col_ptr == NULL
        || to_ints(&jc, lp->col_ptr, (size_t)lp->n + 1) != 0)
        {
        return -1;
    }
    const size_t nnz = (size_t)lp->col_ptr[lp->n];
    lp->row_idx = malloc(sizeof(int) * (nnz + 1));
    lp->values = malloc(sizeof(double) * (nnz + 1));
    if (lp->row_idx == NULL || lp->values == NULL
        || to_ints(&ir, lp->row_idx, nnz) != 0)
        {
        return -1;
    }
    return mat_to_doubles(&pr, lp->values, nnz);
}

int mat_lp_load(const char *path, mat_lp *lp, char *err, size_t err_len)
{
    memset(lp, 0, sizeof(*lp));
    mat_file file;
    if (mat_file_find(path, "Problem", &file, err, err_len) != 0)
    {
        return -1;
    }
    mat_array problem;
    int m = 0, n = 0;
    int status = mat_open_array(&file.variable, &problem);
    if (status == 0 && read_sparse(&problem, lp) != 0)
    {
        status = fail(err, err_len, "cannot read Problem.A");
    }
    if (status == 0 && read_vector(&problem, "b", NULL, &lp->b, &m) != 0)
    {
        status = fail(err, err_len, "cannot read Problem.b");
    }
    if (status == 0
        && read_vector(&problem, "aux", "c", &lp->c, &n) != 0)
        {
        status = fail(err, err_len, "cannot read Problem.aux.c");
    }
    if (status == 0 && (m != lp->m || n != lp->n))
    {
        status = fail(err, err_len, "inconsistent dimensions");
    }
    mat_file_close(&file);
    if (status != 0)
    {
        mat_lp_free(lp);
    }
    return status;
}

void mat_lp_free(mat_lp *lp)
{
    free(lp->col_ptr);
    free(lp->row_idx);
    free(lp->values);
    free(lp->b);
    free(lp->c);
    memset(lp, 0, sizeof(*lp));
}

int mat_lp_load_named(const char *dir, const char *name, mat_lp *lp)
{
    char path[4096], err[256];
    snprintf(path, sizeof(path), "%s/lp_%s.mat", dir, name);
    const int status = mat_lp_load(path, lp, err, sizeof(err));
    if (status != 0)
    {
        fprintf(stderr, "skipping %s: %s\n", name, err);
    }
    return status;
}
