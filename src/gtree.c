/**
 * basic C implementation of the gtree builder
 */

#include "gtree.h"
#include "build_gtree.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

gtree_t *init_gtree_node() { 
    gtree_t *node = malloc(sizeof(gtree_t));
    node->too_full  = 0;
    node->n_matches = 0;
    
    int i;
    for (i = 0; i < 4; i++) {
        node->next[i] = NULL;
    }

    return node;
}

void destroy_gtree( gtree_t *node ) { 
    gtree_t *prev[MAX_WINDOW_SIZE + 1];
    int depth = 0;

    prev[0] = node;

    while (depth >= 0) {
        // check all children of the node.
        int i;
        int dirty = 0;
        for (i = 0; i < 4; i++) {
            if (prev[depth]->next[i] != NULL) {
                prev[depth + 1] = prev[depth]->next[i];
                prev[depth]->next[i] = NULL;
                depth++;
                dirty = 1;
                break; 
            }
        }
        if (!dirty) {
            // no children were found, so we may safely free this node
            free(prev[depth]);
            depth--;
        }
    }
}

int _fewest_matches_in_subtree( gtree_t *node ) {

    int lowest = node->too_full ? INT_MAX : node->n_matches;

    int i, candidate;
    for (i = 0; i < 4; i++) {
        if (node->next[i] != NULL) {
            candidate = _fewest_matches_in_subtree(node->next[i]);
            if (candidate < lowest) {
                lowest = candidate;
            }
        }
    }

    return lowest;
}

void prune_gtree( gtree_t *node ) {

    int i, n_children = 0;
    gtree_t *children[4];
    int nextpos[4];

    for (i = 0; i < 4; i++) {
        if (node->next[i] != NULL) {
            children[n_children] = node->next[i];
            nextpos[n_children] = i;
            n_children++;
        }
    }

    if (n_children == 0) {
        // nothing to be done
        return;
    }
    else if (n_children > 1) {
        // recursively prune
        for (i = 0; i < n_children; i++) {
            prune_gtree(children[i]);
        }
        return;
    }
    // else 1 child
    if (node->too_full) {
        prune_gtree(children[0]);
        return;
    }

    // check if number of matches equal
    if (node->n_matches == children[0]->n_matches) {
        if (children[0]->n_matches == _fewest_matches_in_subtree(children[0])) {
            // then we may prune the child
            node->next[nextpos[0]] = NULL;
            destroy_gtree(children[0]);
        }
    }

    prune_gtree(children[0]);

    return;
}

long count_gtree_nodes( gtree_t *node ) {
    if (node == NULL) {
        return 0;
    }

    long sum = 1;
    int i;
    for (i = 0; i < 4; i++) {
        sum += count_gtree_nodes(node->next[i]);
    }
    return sum;
}
