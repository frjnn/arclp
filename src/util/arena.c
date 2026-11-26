#include "arena.h"

#include <stdlib.h>
#include <string.h>

enum
{
    ARENA_ALIGNMENT = 64
};

static size_t round_up(size_t bytes)
{
    return (bytes + ARENA_ALIGNMENT - 1) & ~(size_t)(ARENA_ALIGNMENT - 1);
}

static void measure(arc_arena *arena)
{
    arena->base = NULL;
    arena->capacity = 0;
    arena->used = 0;
}

static int commit(arc_arena *arena)
{
    size_t bytes = arena->used > 0 ? arena->used : ARENA_ALIGNMENT;
    arena->base = aligned_alloc(ARENA_ALIGNMENT, round_up(bytes));
    if (arena->base == NULL)
    {
        return -1;
    }
    memset(arena->base, 0, round_up(bytes));
    arena->capacity = round_up(bytes);
    arena->used = 0;
    return 0;
}

int arc_arena_build(arc_arena *arena, arc_arena_layout layout, void *owner)
{
    measure(arena);
    layout(arena, owner);
    if (commit(arena) != 0)
    {
        return -1;
    }
    layout(arena, owner);
    return 0;
}

void *arc_arena_take(arc_arena *arena, size_t count, size_t size)
{
    size_t bytes = round_up(count * size);
    void *block = NULL;
    if (arena->base != NULL && arena->used + bytes <= arena->capacity)
    {
        block = arena->base + arena->used;
    }
    arena->used += bytes;
    return block;
}

void arc_arena_take_vectors(arc_arena *arena, double **const *vectors,
                            size_t count, size_t length)
{
    for (size_t k = 0; k < count; k++)
    {
        *vectors[k] = arc_arena_take(arena, length, sizeof(double));
    }
}

void arc_arena_release(arc_arena *arena)
{
    free(arena->base);
    measure(arena);
}
