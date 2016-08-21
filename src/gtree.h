#ifndef GTREE_H
#define GTREE_H

/** gtree.h
 * encapsulate central methods for manipulating the core gtree data structure
 * used by our sequence alignment algorithm.
 */

#include "types.h"
#include "consts.h"


/**
 * malloc's a new gtree node and initializes its "next" array to all null
 * pointers. the "locs" array is left as garbage uninitialized values.
 *
 * @return:
 *      a pointer to a newly initialized gtree node
 */
gtree_t *init_gtree_node();

/**
 * destroy the gtree rooted at "node" by free'ing all nodes. Note this will
 * not free any memory allocated for strings in "locs". These strings must
 * be managed separately.
 *
 * @args:
 *      node - a pointer to the root node of the gtree to be free'd
 * @return:
 *      0        on success
 *      errcode  otherwise
 */
void destroy_gtree( gtree_t *node );

/**
 * prunes the gtree of extraneous nodes. That is:
 *
 * if a node has only one child, and that child has the same number of hits
 * as the parent, then no additional information is carried in that child
 * node and it can be pruned from the gtree.
 *
 * @args:
 *      node - a pointer to a node from which to prune all subtrees.
 */
void prune_gtree ( gtree_t *node );

/**
 * recursively count the numer of nodes in a gtree
 *
 * @args:
 *      node - a pointer to a node from which to count
 */
long count_gtree_nodes( gtree_t *node );

#endif
