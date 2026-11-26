/*
 * Single-allocation arena.  A component describes its buffers once in a
 * layout function; that function runs first against a measuring arena
 * (which only adds up sizes) and then against a real one, so every
 * component performs exactly one heap allocation for all its workspace.
 */
#ifndef ARC_ARENA_H
#define ARC_ARENA_H

#include <stddef.h>

typedef struct
{
    unsigned char *base; /* NULL while measuring */
    size_t capacity;
    size_t used;
} arc_arena;

/* Takes every buffer of `owner` from the arena. */
typedef void (*arc_arena_layout)(arc_arena *arena, void *owner);

/* Runs layout to measure, allocates, and runs it again to hand out the
 * memory; returns 0 on success. */
int arc_arena_build(arc_arena *arena, arc_arena_layout layout, void *owner);

/* Called by layouts.  Returns `count` zero-initialized elements of `size` bytes, aligned to
 * a cache line.  Returns NULL while measuring. */
void *arc_arena_take(arc_arena *arena, size_t count, size_t size);

/* Hands out `count` vectors of `length` doubles. */
void arc_arena_take_vectors(arc_arena *arena, double **const *vectors,
                            size_t count, size_t length);

void arc_arena_release(arc_arena *arena);

#endif
