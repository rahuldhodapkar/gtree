/** align.h
 * perform genome alignment against gtree index.
 */

/**
 * seed matches of single read against reference sequence with a smaller
 * exact match
 *
 * @args:
 *      read - a single quality-bearing read.
 *      ix - a gtree index built on the targeted region
 *      res - an alignment result structure into which to write results
 *
 * @return:
 *      0        on succcess
 *      errcode  otherwise
 *
 */
int seed_matches(read_t *read, ix_t *ix, alnres_t *res);

/**
 * extend seed matches against a reference genome, perhaps by using
 * Smith-Waterman or some other local alignment algorithm
 *
 * ***TODO*** not yet implemented
 *
 * @args:
 *      [alignments]
 *      ref - a set of file pointers to reference sequences.
 *
 * @return:
 *      0        on succcess
 *      errcode  otherwise
 */
int extend_matches(read_t *read, ix_t *ix, alnres_t *res);

/**
 * align single read against reference sequence
 * 
 * @args: 
 *      read - a single quality-bearing read
 *      ix - a gtree index built on the targeted region
 *      res - an alignment result structure into which to write results
 *
 * @return:
 *      0        on succcess
 *      errcode  otherwise
 */
int align_single_read(read_t *read, ix_t *ix, alnres_t *res);
