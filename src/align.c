/** align.c
 * perform genome alignment against gtree index.
 */

#include "align.h"
#include "ssw.h"

/* This table is used to transform nucleotide letters into numbers. */
static const int8_t _NT_TABLE[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int _get_longest_exact_match(bp_t *bp, int max_len, gtree_t *node,
                                                    gtmatch_t *res) {
    int pos = 0;

    res->match_len = 0;
    res->n_matches = node->too_full ? 0    : node->n_matches;
    res->locs =      node->too_full ? NULL : node->locs;

    while (pos < max_len && node->next[ bp[pos] ] != NULL) {
        node = node->next[ bp[pos] ];
        pos++;

        res->match_len++;
        res->n_matches = node->too_full ? 0    : node->n_matches;
        res->locs =      node->too_full ? NULL : node->locs;
    }

    return 0;
}

int seed_matches(read_t *read, ix_t *ix, alnres_t *res) {

    // define a heap in an array.
    align_t seeds[MAX_NUM_SEEDS + 1];
    align_t swp;        // for swap space;

    gtmatch_t cur_match;
    int n_seeds = 0;

    // want to use a heap here to maintain this.
    int pos, i;
    for (pos = 0; pos < read->len; pos++) {
        _get_longest_exact_match(read->seq + pos, read->len - pos,
                                                ix->root, &cur_match);
        if (cur_match.match_len < MIN_SEED_LEN) {
            continue;
        }
        for (i = 0; i < cur_match.n_matches; i++) {

            seeds[n_seeds].template_id = read->template_id;
            seeds[n_seeds].ref = cur_match.locs[i].desc;
            seeds[n_seeds].pos = cur_match.locs[i].pos;
            seeds[n_seeds].seq = read->seq + pos;
            seeds[n_seeds].seq_len = read->len;
            seeds[n_seeds].align_len = cur_match.match_len;
            sprintf(seeds[n_seeds].cigar, "%dM", cur_match.match_len);

            if (n_seeds < MAX_NUM_SEEDS) {
                n_seeds++;
            }
            else {
                // restrict heap
                int heap_pos = n_seeds;
                while (seeds[heap_pos].align_len 
                            > seeds[heap_pos - 1].align_len) {
                    swp = seeds[heap_pos - 1];
                    seeds[heap_pos - 1] = seeds[heap_pos];
                    seeds[heap_pos] = swp;
                    heap_pos--;
                }
            }
        }
    }

    return 0;
    //TODO: copy seeds into alignment result as seed alignments
}

int extend_matches(read_t *read, ix_t *ix, alnres_t *res) {

    // read in bp_t* form
    // read->seq;

    int8_t read_num[read->len];

    int i;
    for (i = 0; i < read->len; i++) {
        read_num[i] = read->seq[i];
    }

    // reference in parsable form
    int8_t ref_num[read->len + 2 * REF_PADDING_LEN];
    
    // reference as characters.
    char ref_seq[read->len + 2 * REF_PADDING_LEN];
    for (i = 0; i < read->len + 2 * REF_PADDING_LEN; i++)  ref_seq[i] = '\0';

    int32_t l, m, k, 
            match = 2, 
            mismatch = 2, 
            gap_open = 3, 
            gap_extension = 1;

    s_align* result;
    s_profile* profile;

    int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
    for (l = k = 0; l < 4; ++l) {
        for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : - mismatch;    
                                          /* weight_match : -weight_mismatch */
        mat[k++] = 0; // ambiguous base: no penalty
    }
    for (m = 0; m < 5; ++m) mat[k++] = 0;

    profile = ssw_init(read_num, 15, mat, 5, 2);
    result = ssw_align (profile, ref_num, 39, gap_open, gap_extension, 1, 0, 0, 15);

    ssw_write(result, ref_seq, read->read_seq, _NT_TABLE);

    align_destroy(result);
    init_destroy(profile);
    free(mat);
    return 0;
    // reference seq
}
