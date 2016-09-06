/** index.c
 * define and manipulate a gtree index
 */

#include "index.h"
#include "gtree.h"

#include "debug.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

ix_t *init_ix() {
    ix_t *ix = malloc(sizeof(ix_t));
    ix->root = init_gtree_node();
    ix->root->too_full = 1;
    ix->n_descs = 0;
    ix->descs = malloc( sizeof(char *) );
    return ix;
}

int destroy_ix( ix_t *ix ) {
    destroy_gtree(ix->root);
    
    int i;
    for (i = 0; i < ix->n_descs; i++) {
        free(*(ix->descs + i));
    }

    free(ix->descs);
    free(ix);
    
    return 0;
}

int _serialize_gtree( gtree_t *node, FILE *out, ix_t *ix ) {

    char has_data;

    if (node == NULL) {
        has_data = 0;
        fwrite(&has_data, sizeof(char), 1, out);
        return 0;
    }

    has_data = 1;
    fwrite(&has_data, sizeof(char), 1, out);

    char too_full = node->too_full ? 1 : 0;
    fwrite(&too_full, sizeof(char), 1, out);

    char n_matches = node->n_matches;
    fwrite(&n_matches, sizeof(char), 1, out);

    // write gtree nodes
    int i;
    for (i = 0; i < 4; i++) {
        _serialize_gtree(node->next[i], out, ix);
    }
    
    // write locs matches
    for (i = 0; i < node->n_matches; i++) {
        // write loc structure
        int j, matchpos = -1;
        for (j = 0; j < ix->n_descs; j++) {
            if (ix->descs[j] == node->locs[i].desc) {
                matchpos = j;
                break;
            }
        }

        if (matchpos == -1 && node->locs[i].desc != NULL) {
            //assert() somehow
            DIE("attempting to serialize corrupted gtree");
        }

        // write loc structure
        fwrite(&matchpos, sizeof(int), 1, out);
        fwrite(&(node->locs[i].pos), sizeof(long), 1, out);
    }

    return 0;
}

int serialize_ix( ix_t *ix, char *outfile ) {
    // needs to be done iteratively to avoid stackoverflow

    FILE *out = fopen(outfile, "w+");
    
    // write header
    // EMPTY right now

    // write n_desc_strings
    fwrite(&(ix->n_descs), sizeof(unsigned int), 1, out);

    // write strings
    int i;
    for (i = 0; i < ix->n_descs; i++) {
        // INT_N_LEN
        size_t desclen = strlen(ix->descs[i]);
        fwrite(&desclen, sizeof(size_t), 1, out);
        // DESC_STRING
        fwrite(ix->descs[i], sizeof(char), strlen(ix->descs[i]) + 1, out);
    }

    // write gtree
    _serialize_gtree(ix->root, out, ix);

    fclose(out);

    return 0;
}

gtree_t *_deserialize_gtree( FILE *in, ix_t *ix ) {

    char has_data;
    fread(&has_data, sizeof(char), 1, in);
    
    if ( !has_data ) {
        return NULL;
    }
    
    gtree_t *node = init_gtree_node();

    char too_full;
    fread(&too_full, sizeof(char), 1, in);
    if (too_full) {
        node->too_full = 1;
    }
    else {
        node->too_full = 0;
    }

    char n_matches;
    fread(&n_matches, sizeof(char), 1, in);
    node->n_matches = n_matches;

    int i;
    for (i = 0; i < 4; i++) {
        node->next[i] = _deserialize_gtree( in, ix );
    }

    for (i = 0; i < node->n_matches; i++) {
        // read loc structure
        int descpos;
        fread(&descpos, sizeof(int), 1, in);

        node->locs[i].desc = descpos < 0 ? NULL
                                         : ix->descs[descpos];

        fread(&(node->locs[i].pos), sizeof(long), 1, in);
    }

    return node;
}

ix_t *deserialize_ix( char *ixfile ) {

    FILE *in = fopen(ixfile, "r");

    ix_t *ix = init_ix();

    // read header
    // EMPTY right now

    // read n_desc_strings
    fread(&(ix->n_descs), sizeof(unsigned int), 1, in);
    ix->descs = realloc(ix->descs, sizeof(char *) * ix->n_descs);

    // read in desc strings
    int i;
    for (i = 0; i < ix->n_descs; i++) {
        // INT_N_LEN
        size_t desclen;
        fread(&desclen, sizeof(size_t), 1, in);

        // DESC_STRING
        ix->descs[i] = malloc(sizeof(char) * (desclen + 1));
        fread(ix->descs[i], sizeof(char), desclen + 1, in);
    }

    // read gtree
    free(ix->root);     // required since init_ix() alloc's a node
    ix->root = _deserialize_gtree(in, ix);

    fclose(in);

    return ix;
}

void print_ix_info( ix_t *ix ) {
    printf("printing index info:\n");
    printf("number of nodes: %ld\n", count_gtree_nodes(ix->root));
    printf("n_descs: %u\n", ix->n_descs);

    int i;
    for (i = 0; i < ix->n_descs; i++) {
        printf("    desc[%d]: %s\n", i, ix->descs[i]);
    }

    printf("done printing info\n");
}
