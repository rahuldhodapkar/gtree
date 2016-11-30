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


void prune_gtree_helper( gtree_t *node, gtree_t *parent, int parent_ix) {

    // recursively prune all potential children
    int i;
    for (i = 0; i < 4; i++) {
        if (node->next[i] != NULL) {
            prune_gtree_helper(node->next[i], node, i);
        }
    }

    int n_children = 0;
    for (i = 0; i < 4; i++) {
        if (node->next[i] != NULL) {
            n_children++;
        }
    }

    if (n_children > 0) {
        return;
    }
    else if (node->too_full && parent != NULL) {
        parent->next[parent_ix] = NULL;
        destroy_gtree(node);
    }
}

void prune_gtree( gtree_t *node ) {
    prune_gtree_helper( node, NULL, -1 );
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
