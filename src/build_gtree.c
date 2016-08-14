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
        iter++; 

        if (iter % 10 == 0) {
            iter = 0;
            printf("Working at desc:%s, pos:%ld, window:%ld\n",
                    descs[*n_descs - 1], cur_pos, cur_window_size);
        }

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
            cur_pos = cur_pos - cur_window_size;
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
    } 

    *desc_strings = descs;
    free(cur_desc);
    fclose(in);

    return 0;
}

ix_t *build_ix_from_ref_seq( char *ref_filename ) {
    ix_t *ix = init_ix();
    build_gtree(ref_filename, &(ix->root), &(ix->descs), &(ix->n_descs));
    return ix;
}

ix_t *init_ix() {
    ix_t *ix = malloc(sizeof(ix_t));
    ix->root = init_gtree_node();
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

    fclose(out);

    return 0;
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

    fclose(in);

    return ix;
}

void print_ix_info( ix_t *ix ) {
    printf("printing index info:\n");
    printf("n_descs: %u\n", ix->n_descs);

    int i;
    for (i = 0; i < ix->n_descs; i++) {
        printf("    desc[%d]: %s\n", i, ix->descs[i]);
    }

    printf("done printing info\n");
}

