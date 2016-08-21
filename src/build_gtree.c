/** build_gtree.c
 * encapsulate methods for building a gtree from a reference genome file.
 */
#include "build_gtree.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// global pushback buffer defs for getc wrappers
char PUSHBACK_BUFFER[MAX_WINDOW_SIZE];
int PUSHBACK_POS = -1;

char bufgetc(FILE *stream) {
    if (PUSHBACK_POS < 0) {
        return getc(stream);
    }
    else {
        return PUSHBACK_BUFFER[PUSHBACK_POS--];
    }
}

int bufungetc(char c) { 
    PUSHBACK_POS++;
    PUSHBACK_BUFFER[PUSHBACK_POS] = c;
    return 0;
}

int read_desc ( FILE *in, char *desc ) {
    char c;
    int head = 0;
    while ((c = fgetc(in)) != '\n') {
        desc[head++] = c;
    }
    desc[head] = '\0';
    return 0;
}

int process_base(bp_t base, gtree_t **cur_node_ref, long pos, char *desc) {
    gtree_t *cur_node = *cur_node_ref;

    if (cur_node->next[base] == NULL) { 
        cur_node->next[base] = init_gtree_node(); 
    }
    // move to the next node
    cur_node = cur_node->next[base];
    
    if (cur_node->n_matches < MAX_LOCS_PER_NODE) {
        cur_node->locs[cur_node->n_matches].desc = desc;
        cur_node->locs[cur_node->n_matches].pos = pos;
        cur_node->n_matches++;
    }
    else if (! cur_node->too_full && cur_node->n_matches == MAX_LOCS_PER_NODE) {
        cur_node->too_full = 1;
    } 

    *cur_node_ref = cur_node;
    return 0;
}

int build_gtree( char *ix_file,
                 gtree_t **gtree_root, 
                 char ***desc_strings,
                 unsigned int *n_descs ) { 
    printf("Building gtree on FASTA input %s.\n", ix_file);

    // declare local copies for readability
    char **descs = *desc_strings;
    gtree_t *root = *gtree_root;

    FILE *in = fopen(ix_file, "r");

    // define current position values
    char *cur_desc = malloc(MAX_DESC_LEN); 
    long cur_pos = 0; 
    *n_descs = 0;

    // define rolling window for build
    long cur_window_size = 0;
    char window_buffer[MAX_WINDOW_SIZE];
    gtree_t *cur_node = root;

    char force_window_rewind = 0;
    long iter = 0;

    char c;
    while ((c = bufgetc(in)) != EOF) {
        if (c == '>' && cur_window_size == 0) { 
            read_desc(in, cur_desc);
            descs = realloc(descs, sizeof(char *)*(*n_descs + 1));
            descs[*n_descs] = malloc(strlen(cur_desc) + 1);
            strcpy(descs[*n_descs], cur_desc);
            *n_descs = *n_descs + 1;
            cur_pos = 0;
            cur_node = root;
            continue;
        }
        else if (c == '>') {
            //rewind window and move forward.
            force_window_rewind = 1;
        }
        else if (c == '\n') {
            continue;
        }
        else if ( (c == 'N' || c == 'n') && cur_window_size == 0) {
            continue;
        }
        else if (c == 'N' || c == 'n') {
            //rewind window and move forward.
            force_window_rewind = 1;
        }
        
        // add current character to buffer
        window_buffer[cur_window_size] = c;
        cur_window_size++;
        cur_pos++;
        
        if (force_window_rewind || cur_window_size == MAX_WINDOW_SIZE) {
            // update window size
            int i;
            for (i = cur_window_size - 1; i > 0; i--) {
                // push back onto buffer
                bufungetc(window_buffer[i]);
            }
            cur_pos = cur_pos - cur_window_size + 1;
            cur_node = root;
            cur_window_size = 0; // delete the first char to move forward
            force_window_rewind = 0;

            continue;
        }

        switch (c) { 
            case 'a':
            case 'A':
                process_base(A, &cur_node, 
                             cur_pos - cur_window_size, descs[*n_descs - 1]);
                break;
                
            case 't':
            case 'T':
                process_base(T, &cur_node, 
                             cur_pos - cur_window_size, descs[*n_descs - 1]);
                break;
                
            case 'c':
            case 'C':
                process_base(C, &cur_node, 
                             cur_pos - cur_window_size, descs[*n_descs - 1]);
                break;

            case 'g':
            case 'G':
                process_base(G, &cur_node, 
                             cur_pos - cur_window_size, descs[*n_descs - 1]);
                break; 
            default:
                printf("ERROR - encountered illegal character [%c|%d] in %s:%ld",
                            c, c, descs[*n_descs - 1], cur_pos - cur_window_size); 
                break;
        }

        if (iter % 1000000 == 0) {
            iter = 0;
            printf("Working at desc:%s, pos:%ld, window:%ld\n",
                    descs[*n_descs - 1], cur_pos, cur_window_size);
        }
        iter++; 

    } 

    *desc_strings = descs;
    free(cur_desc);
    fclose(in);

    return 0;
}

gbuild_t *init_gbuild_node() {
    gbuild_t *node = malloc(sizeof(gbuild_t));
    node->too_full  = 0;
    node->n_matches = 0;
    
    int i;
    for (i = 0; i < 4; i++) {
        node->next[i] = NULL;
    }
    return node;
}

int process_gbase(bp_t base, gbuild_t **cur_node_ref) {
    gbuild_t *cur_node = *cur_node_ref;

    if (cur_node->next[base] == NULL) { 
        cur_node->next[base] = init_gbuild_node(); 
    }
    // move to the next node
    cur_node = cur_node->next[base];
    
    if (cur_node->n_matches < MAX_LOCS_PER_NODE) {
        cur_node->n_matches++;
    }
    else if (! cur_node->too_full && cur_node->n_matches == MAX_LOCS_PER_NODE) {
        cur_node->too_full = 1;
    } 

    *cur_node_ref = cur_node;
    return 0;
}

int build_gbuild_tree( char *ix_file, gbuild_t **gbuild_root ) {
    printf("sizeof(gbuild_t) = %lu\n", sizeof(gbuild_t));

    FILE *in = fopen(ix_file, "r");

    gbuild_t *root = init_gbuild_node();

    // define current position values
    char *cur_desc = malloc(MAX_DESC_LEN); 
    long cur_pos = 0; 

    // define rolling window for build
    long cur_window_size = 0;
    char window_buffer[MAX_WINDOW_SIZE];
    gbuild_t *cur_node = root;

    char force_window_rewind = 0;
    long iter = 0;

    char c;
    while ((c = bufgetc(in)) != EOF) {
        if (c == '>' && cur_window_size == 0) { 
            read_desc(in, cur_desc);
            cur_pos = 0;
            cur_node = root;
            continue;
        }
        else if (c == '>') {
            //rewind window and move forward.
            force_window_rewind = 1;
        }
        else if (c == '\n') {
            continue;
        }
        else if ( (c == 'N' || c == 'n') && cur_window_size == 0) {
            continue;
        }
        else if (c == 'N' || c == 'n') {
            //rewind window and move forward.
            force_window_rewind = 1;
        }
        
        // add current character to buffer
        window_buffer[cur_window_size] = c;
        cur_window_size++;
        cur_pos++;
        
        if (force_window_rewind || cur_window_size == MAX_WINDOW_SIZE) {
            // update window size
            int i;
            for (i = cur_window_size - 1; i > 0; i--) {
                // push back onto buffer
                bufungetc(window_buffer[i]);
            }
            cur_pos = cur_pos - cur_window_size + 1;
            cur_node = root;
            cur_window_size = 0; // delete the first char to move forward
            force_window_rewind = 0;

            continue;
        }

        switch (c) { 
            case 'a':
            case 'A':
                process_gbase(A, &cur_node);
                break;
                
            case 't':
            case 'T':
                process_gbase(T, &cur_node);
                break;
                
            case 'c':
            case 'C':
                process_gbase(C, &cur_node);
                break;

            case 'g':
            case 'G':
                process_gbase(G, &cur_node);
                break; 
            default:
                printf("ERROR - encountered illegal character [%c|%d] in %s:%ld",
                            c, c, cur_desc, cur_pos - cur_window_size); 
                break;
        }

        if (iter % 100000000 == 0) {
            iter = 0;
            printf("Pre-building at desc:%s, pos:%ld, window:%ld\n",
                    cur_desc, cur_pos, cur_window_size);
        }
        iter++; 

    } 

    *gbuild_root = root;
    free(cur_desc);
    fclose(in);

    return 0;
}

void destroy_gbuild( gbuild_t *node ) { 
    gbuild_t *prev[MAX_WINDOW_SIZE + 1];
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

void prune_gtree_builder( gbuild_t *node ) {

    int i, n_children = 0;
    gbuild_t *children[4];
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
            prune_gtree_builder(children[i]);
        }
        return;
    }
    // else 1 child
    if (node->too_full) {
        prune_gtree_builder(children[0]);
        return;
    }

    // check if number of matches equal
    if (node->n_matches == children[0]->n_matches) {
        // then we may prune the child
        node->next[nextpos[0]] = NULL;
        destroy_gbuild(children[0]);
    }
    return;
}

long count_builder_nodes( gbuild_t *node ) {
    if (node == NULL) {
        return 0;
    }

    long sum = 1;
    int i;
    for (i = 0; i < 4; i++) {
        sum += count_builder_nodes(node->next[i]);
    }
    return sum;
}

ix_t *build_ix_from_ref_seq( char *ref_filename ) {
    ix_t *ix = init_ix();
    //build_gtree(ref_filename, &(ix->root), &(ix->descs), &(ix->n_descs));
    gbuild_t *build_root;
    build_gbuild_tree(ref_filename, &build_root);
    printf("intermediary tree built with %ld nodes\n", 
                                    count_builder_nodes(build_root));
    prune_gtree_builder(build_root);
    printf("%ld nodes after pruning\n",
                                    count_builder_nodes(build_root));
    return ix;
}
