/* Tagged-element layer of the level-5 MAT format (little-endian files). */
#include "mat5_format.h"

#include <string.h>

static uint32_t read_u32(const uint8_t *p)
{
    return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16
         | (uint32_t)p[3] << 24;
}

int mat_next(mat_cursor *cur, mat_element *el)
{
    if (cur->end - cur->p < 8)
    {
        return -1;
    }
    const uint32_t word = read_u32(cur->p);
    if (word >> 16 != 0) /* small element: data in the tag */
    {
        el->type = word & 0xffffu;
        el->bytes = word >> 16;
        el->data = cur->p + 4;
        cur->p += 8;
        return el->bytes <= 4 ? 0 : -1;
    }
    el->type = word;
    el->bytes = read_u32(cur->p + 4);
    el->data = cur->p + 8;
    if ((size_t)(cur->end - el->data) < el->bytes)
    {
        return -1;
    }
    const size_t padded = el->type == MI_COMPRESSED
                              ? el->bytes : (el->bytes + 7u) & ~7u;
    const size_t left = (size_t)(cur->end - el->data);
    cur->p = el->data + (padded < left ? padded : left);
    return 0;
}

static size_t width(uint32_t type)
{
    switch (type)
    {
    case MI_INT8: case MI_UINT8: return 1;
    case MI_INT16: case MI_UINT16: return 2;
    case MI_INT32: case MI_UINT32: case MI_SINGLE: return 4;
    case MI_DOUBLE: case MI_INT64: case MI_UINT64: return 8;
    default: return 0;
    }
}

size_t mat_count(const mat_element *el)
{
    const size_t w = width(el->type);
    return w == 0 ? 0 : el->bytes / w;
}

static double convert(uint32_t type, const uint8_t *p)
{
    int8_t i8; uint8_t u8; int16_t i16; uint16_t u16; int32_t i32;
    uint32_t u32; float f; double d; int64_t i64; uint64_t u64;
    switch (type)
    {
    case MI_INT8: memcpy(&i8, p, 1); return i8;
    case MI_UINT8: memcpy(&u8, p, 1); return u8;
    case MI_INT16: memcpy(&i16, p, 2); return i16;
    case MI_UINT16: memcpy(&u16, p, 2); return u16;
    case MI_INT32: memcpy(&i32, p, 4); return i32;
    case MI_UINT32: memcpy(&u32, p, 4); return u32;
    case MI_SINGLE: memcpy(&f, p, 4); return f;
    case MI_INT64: memcpy(&i64, p, 8); return (double)i64;
    case MI_UINT64: memcpy(&u64, p, 8); return (double)u64;
    default: memcpy(&d, p, 8); return d;
    }
}

int mat_to_doubles(const mat_element *el, double *out, size_t count)
{
    const size_t w = width(el->type);
    if (w == 0 || mat_count(el) < count)
    {
        return -1;
    }
    for (size_t k = 0; k < count; k++)
    {
        out[k] = convert(el->type, el->data + k * w);
    }
    return 0;
}
