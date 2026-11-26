#include "step_rule.h"

#include <stdlib.h>

arc_step_rule *arc_step_rule_create(const arc_step_rule_ops *ops, int m,
                                    int n)
{
    arc_step_rule *rule = calloc(1, ops->size);
    if (rule == NULL)
    {
        return NULL;
    }
    rule->ops = ops;
    rule->m = m;
    rule->n = n;
    if (arc_arena_build(&rule->arena, ops->layout, rule) != 0)
    {
        free(rule);
        return NULL;
    }
    if (ops->init != NULL)
    {
        ops->init(rule);
    }
    return rule;
}

void arc_step_rule_destroy(arc_step_rule *rule)
{
    if (rule != NULL)
    {
        arc_arena_release(&rule->arena);
        free(rule);
    }
}
