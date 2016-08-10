/**
 * basic C implementation of the gtree builder
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LOCS_PER_NODE 4
#define MAX_DESC_LEN 100
#define MAX_WINDOW_SIZE 10

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
 * copy a description string from "in" to desc (terminated by a newline
 * character)
 *
 * ASSUME:
 *      - no description strings are greater than MAX_DESC_LEN chars long.
 *
 * @args:
 *      in - FILE to read from 
 *      desc - string to read data into
 *
 * @return:
 *      0        on success
 *      errcode  otherwise
 *
 */
int read_desc ( FILE *in, char *desc ) {
    char c;
    int head = 0;
    while ((c = fgetc(in)) != '\n') {
        desc[head++] = c;
    }
    desc[head] = '\0';
    return 0;
}

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

/**
 * process a new base and advance the current node
 *
 * @args:
 *      base - the base pair to process
 *      cur_node_ref - a reference to the current node, which will be updated
 *                     with "base".
 *      desc - a description of the sequence being indexed.
 *      pos - the current position in sequence being indexed.
 */
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

/**
 * build gtree index, assign to "gtree_root".
 *
 * ASSUME:
 *      - no description strings are greater than MAX_DESC_LEN chars long.
 *
 * @args:
 *      ix_file - FASTA file to 
 *      gtree_root - pointer to to assign gtree_root to
 *      desc_strings - pointer to assign description strings to for free later
 *      n_descs - number of description strings in desc_strings
 *
 * @return:
 *      0        on success
 *      errcode  otherwise
 *
 */
int build_gtree( char *ix_file,
                 gtree_t **gtree_root, 
                 char ***desc_strings,  // pointer to an array of strings
                 unsigned int *n_descs ) { 
    printf("Building gtree on FASTA input %s.\n", ix_file);

    *desc_strings = malloc(sizeof(char *)); 
    *gtree_root = init_gtree_node(); 

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

/**
 * destroy gtree index
 */
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

/**
 * destroy reserved description strings
 */ 
void destroy_strings( int n_strings, char **strings ) {
    int i;
    for (i = 0; i < n_strings; i++) {
        free(*(strings + i));
    }
    free(strings);
}

/**
 * main function to build gtree from fasta file
 */
int main(int argc, char *argv[]) {
    char *fgenome = NULL, *freads = NULL;

    if (argc <= 1) {
        printf("Usage: gtree [genome_to_index] [read_file]\n");
    }

    if (argc > 1) {
        printf("input file [%s] passed.\n", argv[1]);
        fgenome = argv[1];
    }

    if (argc > 2) {
        printf("reads file [%s] passed.\n", argv[2]);
        freads = argv[2];
    }

    gtree_t *ix_root = NULL;
    char **desc_strings = NULL;
    unsigned int n_descs = 0;
    
    printf("Building!\n");
    build_gtree(fgenome, &ix_root, &desc_strings, &n_descs);

    printf("Cleaning up!\n");
    destroy_gtree(ix_root);
    printf("destroy %d strings\n", n_descs);
    destroy_strings(n_descs, desc_strings);

    printf("finished running!\n");
    return 0;
}

