#ifndef INDEX_H
#define INDEX_H

#include "types.h"
#include "consts.h"

/** index.h
 * definition of gtree index structure and associated methods to manipulate 
 * them.
 */

/**
 * malloc all structures to be used by the gtree index. Specifically, any
 * additional descriptions will require a realloc() of the descs pointer
 *
 * @return:
 *      pointer to an allocated an initialized ix_t structure.
 */
ix_t *init_ix();

/**
 * free all structures used by "ix"
 *
 * @args:
 *      ix - the gtree index to be free'd
 * @return:
 *      0        on success
 *      errcode  otherwise
 */
int destroy_ix( ix_t *ix );

/**
 * serialize a gtree index into a standard on-disk representation format 
 * that can be ported to avoid the cost of re-computing an index from source.
 * 
 * on-disk format:
 * 
 * IX_SER := HEADER
 *           INT_N_DESC_STRINGS
 *           DESC_STRING (x INT_N_DESC_STRINGS)
 *           GTREE_NODE
 *
 * DESC_STRING := INT_N_LEN
 *                CHAR (x INT_N_LEN)
 *
 * GTREE_NODE := NULL             # no data
 *             | HAS_DATA         # non-zero flag to indicate data
 *               INT_TOO_FULL
 *               INT_N_MATCHES    # determines number of locs in serialization
 *               GTREE_NODE       # A
 *               GTREE_NODE       # C
 *               GTREE_NODE       # T
 *               GTREE_NODE       # G
 *               LOC_STRUCT (x INT_N_MATCHES)
 *
 * @args:
 *      ix - a pointer to the index to be serialized
 *      outfile - the name of the file to serialize the tree to
 * @return:
 *      0        on succcess
 *      errcode  otherwise
 */
int serialize_ix( ix_t *ix, char *outfile );

/**
 * deserialize a gtree index stored with "serialize_gtree" into an in-memory
 * representation.
 * 
 * @args:
 *      ixfile - the name of the file to read an index from.
 *
 * @return:
 *      a pointer to the deserialized index. When finished, this pointer
 *          should be freed with "destroy_ix"
 *      NULL if there is an error during index deserialization
 */
ix_t *deserialize_ix( char *ixfile );

/**
 * prints some information about the gtree index supplied to STDOUT
 *
 * @args:
 *      ix - an index pointer to print information about.
 */
void print_ix_info( ix_t *ix );

#endif
