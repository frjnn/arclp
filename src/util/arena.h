/*
 * Single-allocation arena.  A component describes its buffers once in a
 * layout function; that function runs first against a measuring arena
 * (which only adds up sizes) and then against a real one, so every
 * component performs exactly one heap allocation for all its workspace.
 */
#ifndef ARC_ARENA_H
#define ARC_ARENA_H

#include <stddef.h>

typedef struct {
    unsigned char *base;   /* NULL while measuring */
    size_t capacity;
    size_t used;
} arc_arena;

/* Starts a measuring pass: arc_arena_take returns NULL and counts bytes. */
void arc_arena_measure(arc_arena *arena);

/* Allocates the measured number of bytes; returns 0 on success. */
int arc_arena_commit(arc_arena *arena);

/* Returns `count` zero-initialized elements of `size` bytes, aligned to
 * a cache line.  Returns NULL while measuring. */
void *arc_arena_take(arc_arena *arena, size_t count, size_t size);

void arc_arena_release(arc_arena *arena);

#endif
