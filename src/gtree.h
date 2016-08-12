#ifndef GTREE_H
#define GTREE_H

/** gtree.h
 * encapsulate central methods for manipulating the core gtree data structure
 * used by our sequence alignment algorithm.
 */

#define MAX_LOCS_PER_NODE 4
#define MAX_DESC_LEN 100
#define MAX_WINDOW_SIZE 10

typedef enum bp {
    A = 0,
    C = 1,
    T = 2,
    G = 3,
} bp_t;

typedef struct loc {
    char *desc;
    long pos;
} loc_t;

typedef struct gtree {
    unsigned int too_full : 1;        // bit flag to ignore intermediate matching
    unsigned int n_matches: 7;        // number of matches (is_match if >0)

    struct gtree *next[4];            // core gtree lookup array

    loc_t locs [MAX_LOCS_PER_NODE];   // get number of locs to match per contig
} gtree_t;

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

#endif
