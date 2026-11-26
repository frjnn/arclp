#include "mat5_file.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

enum { HEADER_BYTES = 128 };

static uint8_t *read_all(const char *path, size_t *size)
{
    FILE *f = fopen(path, "rb");
    if (f == NULL)
    {
        return NULL;
    }
    uint8_t *data = NULL;
    if (fseek(f, 0, SEEK_END) == 0)
    {
        const long length = ftell(f);
        data = length > 0 ? malloc((size_t)length) : NULL;
        rewind(f);
        if (data != NULL && fread(data, 1, (size_t)length, f) != (size_t)length)
        {
            free(data);
            data = NULL;
        }
        *size = length > 0 ? (size_t)length : 0;
    }
    fclose(f);
    return data;
}

/* Inflates a miCOMPRESSED payload into a freshly allocated buffer. */
static uint8_t *inflate_all(const mat_element *el, size_t *size)
{
    size_t capacity = 4 * (size_t)el->bytes + 1024;
    uint8_t *out = malloc(capacity);
    z_stream z;
    memset(&z, 0, sizeof(z));
    if (out == NULL || inflateInit(&z) != Z_OK)
    {
        free(out);
        return NULL;
    }
    z.next_in = (Bytef *)(uintptr_t)el->data;
    z.avail_in = el->bytes;
    int status = Z_OK;
    while (status == Z_OK)
    {
        if (z.total_out == capacity)
        {
            uint8_t *bigger = realloc(out, 2 * capacity);
            if (bigger == NULL)
            {
                break;
            }
            out = bigger;
            capacity *= 2;
        }
        z.next_out = out + z.total_out;
        z.avail_out = (uInt)(capacity - z.total_out);
        status = inflate(&z, Z_NO_FLUSH);
    }
    *size = z.total_out;
    inflateEnd(&z);
    if (status != Z_STREAM_END)
    {
        free(out);
        return NULL;
    }
    return out;
}

static int consider(mat_file *file, const mat_element *el, const char *name)
{
    if (el->type == MI_MATRIX)
    {
        file->variable = *el;
        return mat_array_named(el, name);
    }
    if (el->type != MI_COMPRESSED)
    {
        return 0;
    }
    size_t size = 0;
    uint8_t *raw = inflate_all(el, &size);
    mat_cursor inner = {raw, raw + size};
    if (raw != NULL && mat_next(&inner, &file->variable) == 0
        && mat_array_named(&file->variable, name))
        {
        free(file->inflated);
        file->inflated = raw;
        return 1;
    }
    free(raw);
    return 0;
}

int mat_file_find(const char *path, const char *name, mat_file *file,
                  char *err, size_t err_len)
{
    memset(file, 0, sizeof(*file));
    size_t size = 0;
    file->bytes = read_all(path, &size);
    if (file->bytes == NULL || size < HEADER_BYTES
        || memcmp(file->bytes + 126, "IM", 2) != 0)
        {
        snprintf(err, err_len, "%s: not a little-endian MAT-file", path);
        mat_file_close(file);
        return -1;
    }
    mat_cursor cur = {file->bytes + HEADER_BYTES, file->bytes + size};
    mat_element el;
    while (mat_next(&cur, &el) == 0)
    {
        if (consider(file, &el, name))
        {
            return 0;
        }
    }
    snprintf(err, err_len, "%s: no variable named %s", path, name);
    mat_file_close(file);
    return -1;
}

void mat_file_close(mat_file *file)
{
    free(file->bytes);
    free(file->inflated);
    memset(file, 0, sizeof(*file));
}
