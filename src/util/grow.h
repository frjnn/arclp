/* Growable arrays. */
#ifndef ARC_GROW_H
#define ARC_GROW_H

#include <stdlib.h>
#include <string.h>

/* Smallest doubling of capacity (from `start`) that holds `needed`. */
static inline int arc_grown_capacity(int capacity, int needed, int start)
{
    int grown = capacity > 0 ? 2 * capacity : start;
    while (grown < needed)
    {
        grown *= 2;
    }
    return grown;
}

/* Resizes the array whose pointer is at `address` (any object pointer type)
 * to count elements; on failure (-1) it is unchanged.  memcpy keeps this
 * valid for every pointer type. */
static inline int arc_resize(void *address, size_t count, size_t size)
{
    void *old;
    memcpy(&old, address, sizeof(old));
    void *moved = realloc(old, count * size);
    if (moved == NULL)
    {
        return -1;
    }
    memcpy(address, &moved, sizeof(moved));
    return 0;
}

#endif
