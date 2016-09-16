#ifndef TYPES_H
#define TYPES_H

/** types.h
 * define all types to be used in aligner
 */

#include "consts.h"

#include <stdio.h> 

typedef struct args {
    int exec_mode;
    int verbosity;
    char *ref_fasta_fn;
    char *ix_fn;
    char *out_fn;
    int out_format;
    char *in_fn;        // reads input file for first pair or unpaired
    char *in_fn2;       // reads input file for second pair
    char print_num;     // print num flag for `gtree ix stat`
} args_t;

typedef enum bp {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4,
    NOBP = 5
} bp_t;

typedef struct loc {
    char *desc;
    unsigned long pos;
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

typedef struct read {
    char template_id[MAX_DESC_LEN];       // unique id of read template
                                          // - empty string if anonymous
    bp_t *seq;                            // base pair sequence of read
    char *read_seq;                       // char sequence for read
    char *phred;                          // phred quality of read
    int len;                              // number of base pairs in read
    int malloc_len;                       // malloc'ed seq length
} read_t;

typedef struct gtreematch {
    unsigned int match_len;
    unsigned int n_matches;
    loc_t *locs;
} gtmatch_t;

/**
 * ***NOTE*** all alignments are reported on the "forward" strand
 */
typedef struct alignment {
    char *template_id;   // unique read template id
    char *desc;          // reference sequence name
    unsigned long pos;  // position of alignment start in reference
    bp_t *seq;           // sequence of the template
    int seq_len;         // length of the template sequence
    int align_len;       // length of the aligned sequence
    char cigar[MAX_CIGAR_STR_LEN];  // concise idiosyncratic gapped alignment report
} align_t;

typedef struct alignment_result {
    align_t alns[MAX_NUM_SEEDS + 1];
    int n_alns;
} alnres_t;

typedef struct desc_loc_map {
    FILE *ref_file;
    char ref_filename[FILENAME_MAX];
    char desc[MAX_DESC_LEN];
    unsigned long contig_start;
    unsigned long contig_len;
    struct desc_loc_map *next;
} desc_loc_map_t;

typedef struct reference {
    desc_loc_map_t *desc_map;
    char cur_desc[MAX_DESC_LEN];
    unsigned long cur_pos;
} ref_t;

#endif
