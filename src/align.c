/** align.c
 * perform genome alignment against gtree index.
 */

#include "align.h"

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

    return 0;
    // reference seq
}
