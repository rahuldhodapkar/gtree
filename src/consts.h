#ifndef CONSTS_H
#define CONSTS_H

/** consts.h
 * define all constants to be used in aligner
 */

// maximum length for a sequence description in a FASTA file
#define MAX_DESC_LEN 100

// maximum window size to search during gtree index construction.
#define MAX_WINDOW_SIZE 32

// maximum number of hits per node before declaring "too_full"
#define MAX_LOCS_PER_NODE 4

#endif
