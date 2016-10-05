#ifndef CONSTS_H
#define CONSTS_H

/** consts.h
 * define all constants to be used in aligner
 *
 * also includes basic limits from <limits.h>
 */

#include <limits.h>


/////////////////////////////////////////////////////////////////////////
// ERROR CODES + DEBUG
/////////////////////////////////////////////////////////////////////////

#define ERRCODE_GENERAL_ERROR 1
#define ERRCODE_FILE_NOT_FOUND 2

#define DEBUG_BT_BUFFER_SIZE 100

#define VERBOSITY_LEVEL 5

// verbosity of execution for debug logging / profiling
#define VERBOSITY_LEVEL_QUIET 0
#define VERBOSITY_LEVEL_INFO 1
#define VERBOSITY_LEVEL_WARN 2
#define VERBOSITY_LEVEL_DEBUG 3

/////////////////////////////////////////////////////////////////////////
// ARUGMENT PARSING
/////////////////////////////////////////////////////////////////////////

// execution modes for aligner used in args_t
#define EXEC_MODE_IX_BUILD 0
#define EXEC_MODE_IX_MASK 1
#define EXEC_MODE_IX_PRUNE 2
#define EXEC_MODE_IX_STAT 3
#define EXEC_MODE_ALN 100

// output format modes
#define OUTPUT_FORMAT_SAM 0
#define OUTPUT_FORMAT_BAM 1

/////////////////////////////////////////////////////////////////////////
// INDEX BUILD
/////////////////////////////////////////////////////////////////////////

// maximum length for a sequence description in a FASTA file
#define MAX_DESC_LEN 100

// maximum window size to search during gtree index construction.
#define MAX_WINDOW_SIZE 32

// maximum number of hits per node before declaring "too_full"
#define MAX_LOCS_PER_NODE 4

/////////////////////////////////////////////////////////////////////////
// ALIGNMENT
/////////////////////////////////////////////////////////////////////////

// minimum length for a seed to be considered
#define MIN_SEED_LEN 20

#define MAX_NUM_SEEDS 10

#define MAX_CIGAR_STR_LEN 200

// define how much padding should be done when pulling reference
// sequence from a seeded alignment
#define REF_PADDING_LEN 50

#endif
