/** align.c
 * perform genome alignment against gtree index.
 */

#include "align.h"
#include "ssw.h"
#include "debug.h"
#include "ref.h"

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

bp_t char_to_bp(char c) {
    switch (c) {
        case 'a':
        case 'A':
            return A;
        case 'c':
        case 'C':
            return C;
        case 'g':
        case 'G':
            return G;
        case 't':
        case 'T':
            return T;
        default:
            return N;
    }
}

char bp_to_char (bp_t b) {
    switch (b) {
        case A:
            return 'A';
        case C:
            return 'C';
        case G:
            return 'G';
        case T:
            return 'T';
        case N:
            return 'N';
        default:
            DIE("Unable to marshall bp_t -> char");
            return '\0';
    }
}

int _get_longest_exact_match(bp_t *bp, int max_len, gtree_t *node,
                                                    gtmatch_t *res) {
    int pos = 0;

    res->match_len = 0;
    res->n_matches = node->too_full ? 0    : node->n_matches;
    res->locs =      node->too_full ? NULL : node->locs;

    while (pos < max_len 
            && bp[pos] != N 
            && (node->next[ bp[pos] ] != NULL) )  {
        node = node->next[ bp[pos] ];
        pos++;

        res->match_len++;
        res->n_matches = node->too_full ? 0    : node->n_matches;
        res->locs =      node->too_full ? NULL : node->locs;
    }

    return 0;
}

int get_next_read(FILE *read_file, read_t *read) {
    char c;
    char fastq_line = 1;
    int cur_line_pos = 0;

    while ( (c = fgetc(read_file)) != EOF) {
        if (c == '\n') {
            if (fastq_line == 1) {
                read->template_id[cur_line_pos] = '\0';
            }

            if (fastq_line == 2) {
                read->len = cur_line_pos;
            }

            fastq_line++;
            cur_line_pos = 0;
        }
        if (fastq_line > 4) {
            // read processing completed.
            return 1;
        }

        if (fastq_line == 1) {
            // header line
            if (cur_line_pos == 0) {
                // skip initial character
                if (c != '@') DIE("Malformed FASTQ file");
                cur_line_pos++;
                continue;
            }
            read->template_id[cur_line_pos - 1] = c;
        } else if ( fastq_line == 2 ) {
            if (cur_line_pos >= read->malloc_len) {
                read->malloc_len++;
                // realloc key pointers
                read->seq = realloc(read->seq, sizeof(bp_t) * read->malloc_len);
                read->read_seq = realloc(read->read_seq, sizeof(char) * read->malloc_len);
                read->phred = realloc(read->phred, sizeof(char) * read->malloc_len);
            }

            read->seq[cur_line_pos] = char_to_bp(c);
            read->read_seq[cur_line_pos] = toupper(c);

        } else if (fastq_line == 4) {
            read->phred[cur_line_pos] = c;
        }

        cur_line_pos++;
    }

    return 0;
}

int seed_matches(read_t *read, ix_t *ix, alnres_t *res) {

    // define a heap in an array.
    align_t swp;        // for swap space;

    gtmatch_t cur_match;
    res->n_alns = 0;

    // want to use a heap here to maintain this.
    int pos, i;
    for (pos = 0; pos < read->len; pos++) {
        _get_longest_exact_match(read->seq + pos, read->len - pos,
                                                ix->root, &cur_match);
        if (cur_match.match_len < MIN_SEED_LEN) {
            continue;
        }
        for (i = 0; i < cur_match.n_matches; i++) {

            res->alns[res->n_alns].template_id = read->template_id;
            res->alns[res->n_alns].desc = cur_match.locs[i].desc;
            res->alns[res->n_alns].pos = cur_match.locs[i].pos;
            res->alns[res->n_alns].seq = read->seq + pos;
            res->alns[res->n_alns].seq_len = read->len;
            res->alns[res->n_alns].align_len = cur_match.match_len;
            sprintf(res->alns[res->n_alns].cigar, "%dM", cur_match.match_len);

            if (res->n_alns < MAX_NUM_SEEDS) {
                res->n_alns++;
            }
            else {
                // restrict heap
                int heap_pos = res->n_alns;
                while (res->alns[heap_pos].align_len 
                            > res->alns[heap_pos - 1].align_len) {
                    swp = res->alns[heap_pos - 1];
                    res->alns[heap_pos - 1] = res->alns[heap_pos];
                    res->alns[heap_pos] = swp;
                    heap_pos--;
                }
            }
        }
    }

    return 0;
    //TODO: copy seeds into alignment result as seed alignments
}

int _extend_single_match(read_t *read, ix_t *ix, ref_t *ref, char *desc, long pos) {

    printf("starting to extend a single match\n");
    // read in bp_t* form
    // read->seq;
    int32_t l, m, k, i, 
        match = 2, 
        mismatch = 2, 
        gap_open = 3, 
        gap_extension = 1;

    // reference in parsable form
    unsigned long ref_bp_len = read->len + 2 * REF_PADDING_LEN;
    bp_t ref_bp_string[ref_bp_len];
    refcpy(ref, desc, pos, ref_bp_len, ref_bp_string, &ref_bp_len);

    int8_t ref_num[ref_bp_len];
    char ref_seq[ref_bp_len + 1]; // +1 for '\0'

    for (i = 0; i < ref_bp_len; i++) {
        ref_num[i] = ref_bp_string[i];
        ref_seq[i] = bp_to_char(ref_bp_string[i]);
    }
    ref_seq[ref_bp_len] = '\0';

    printf("DEBUG: read sequence [%s]\n", read->read_seq);
    printf("DEBUG: ref  sequence [%s]\n", ref_seq);

    int8_t read_num[read->len];
    for (i = 0; i < read->len; i++) read_num[i] = read->seq[i];

    // reference as characters.
    s_align *result;
    s_profile *profile;

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

int align_single_read(read_t *read, ix_t *ix, ref_t *ref, alnres_t *res) {
    seed_matches(read, ix, res);

    printf("DEBUG: align read [%s] : [%s]\n", read->read_seq, read->phred);
    printf("DEBUG: found %d seed matches from index\n", res->n_alns);

    int i;
    for (i = 0; i < res->n_alns; i++) {
        _extend_single_match(read, ix, ref, res->alns[i].desc, res->alns[i].pos);
    }

    return 0;
}
