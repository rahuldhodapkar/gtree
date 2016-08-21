#ifndef TYPES_H
#define TYPES_H

/** types.h
 * define all types to be used in aligner
 */

 #include "consts.h"

typedef struct args {
    int exec_mode;
    int verbosity;
    char *ref_fasta_fn;
    char *ix_fn;
    char *out_fn;
    int out_format;
    char *in_fn;        // reads input file for first pair or unpaired
    char *in_fn2;       // reads input file for second pair
} args_t;

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

typedef struct gtreeix {
    gtree_t *root;           // root of gtree index
    unsigned int n_descs;    // number of description strings in gtree
    char **descs;            // access to all description strings in gtree
} ix_t;

#endif
