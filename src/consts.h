#ifndef CONSTS_H
#define CONSTS_H

/** consts.h
 * define all constants to be used in aligner
 */

// verbosity of execution for debug logging / profiling
#define VERBOSITY_LEVEL_QUIET 0
#define VERBOSITY_LEVEL_DEBUG 1

// execution modes for aligner used in args_t
#define EXEC_MODE_IX_BUILD 0
#define EXEC_MODE_IX_MASK 1
#define EXEC_MODE_IX_PRUNE 2
#define EXEC_MODE_IX_STAT 3

#define EXEC_MODE_ALN 100

// output format modes
#define OUTPUT_FORMAT_SAM 0
#define OUTPUT_FORMAT_BAM 1

// maximum length for a sequence description in a FASTA file
#define MAX_DESC_LEN 100

// maximum window size to search during gtree index construction.
#define MAX_WINDOW_SIZE 32

// maximum number of hits per node before declaring "too_full"
#define MAX_LOCS_PER_NODE 4

#endif
