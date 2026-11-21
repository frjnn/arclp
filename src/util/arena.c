#include "arena.h"

#include <stdlib.h>
#include <string.h>

enum { ARENA_ALIGNMENT = 64 };

static size_t round_up(size_t bytes)
{
    return (bytes + ARENA_ALIGNMENT - 1) & ~(size_t)(ARENA_ALIGNMENT - 1);
}

void arc_arena_measure(arc_arena *arena)
{
    arena->base = NULL;
    arena->capacity = 0;
    arena->used = 0;
}

int arc_arena_commit(arc_arena *arena)
{
    size_t bytes = arena->used > 0 ? arena->used : ARENA_ALIGNMENT;
    arena->base = aligned_alloc(ARENA_ALIGNMENT, round_up(bytes));
    if (arena->base == NULL) {
        return -1;
    }
    memset(arena->base, 0, round_up(bytes));
    arena->capacity = round_up(bytes);
    arena->used = 0;
    return 0;
}

void *arc_arena_take(arc_arena *arena, size_t count, size_t size)
{
    size_t bytes = round_up(count * size);
    void *block = NULL;
    if (arena->base != NULL && arena->used + bytes <= arena->capacity) {
        block = arena->base + arena->used;
    }
    arena->used += bytes;
    return block;
}

void arc_arena_release(arc_arena *arena)
{
    free(arena->base);
    arena->base = NULL;
    arena->capacity = 0;
    arena->used = 0;
}
