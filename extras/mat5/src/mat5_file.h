/* File layer of the MAT reader: header check, zlib, variable lookup. */
#ifndef ARC_MAT5_FILE_H
#define ARC_MAT5_FILE_H

#include "mat5_format.h"

typedef struct {
    uint8_t *bytes;        /* whole file              */
    uint8_t *inflated;     /* decompressed variable   */
    mat_element variable;  /* the miMATRIX element    */
} mat_file;

int mat_file_find(const char *path, const char *name, mat_file *file,
                  char *err, size_t err_len);

void mat_file_close(mat_file *file);

#endif
