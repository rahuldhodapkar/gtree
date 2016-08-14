#ifndef BUILD_GTREE_H
#define BUILD_GTREE_H

/** build_gtree.h
 * encapsulate methods for building a gtree from a reference genome file.
 * 
 * Currently, only FASTA-formatted reference genomes in the style of those
 * provided by the UCSC genome browser portal are supported.
 */

#include "gtree.h"

#include <stdlib.h>
#include <stdio.h>

// maximum length for a sequence description in a FASTA file
#define MAX_DESC_LEN 100

// maximum window size to search during gtree index construction.
#define MAX_WINDOW_SIZE 40

typedef struct gtreeix {
    gtree_t *root;      // root of gtree index
    unsigned int n_descs;        // number of description strings in gtree
    char **descs;       // access to all description strings in gtree
} ix_t;

/**
 * "getc()"-like interface to a file, with MAX_WINDOW_SIZE allowable chars of
 * pushback in a global buffer, allowing for efficient and easy-to-understand
 * index contruction logic.
 *
 * @args:
 *      stream - a FILE pointer from which to get the next char
 * @return
 *      the last character pushed back with "bufungetc" or the next character
 *      in "stream" if the buffer is empty.
 */
char bufgetc(FILE *stream);

/**
 * push back character to buffer, will return an error code if buffer is full.
 *
 * @args:
 *      c - the character to push back
 *
 * @return:
 *      0         on success
 *      errcode   otherwise
 */
int bufungetc(char c);

/**
 * copy a description string from "in" to desc (terminated by a newline
 * character)
 *
 * ASSUME:
 *      - no description strings are greater than MAX_DESC_LEN chars long.
 *
 * @args:
 *      in - FILE to read from 
 *      desc - string to read data into
 *
 * @return:
 *      0        on success
 *      errcode  otherwise
 *
 */
int read_desc ( FILE *in, char *desc );

/**
 * process a new base and advance the current node
 *
 * @args:
 *      base - the base pair to process
 *      cur_node_ref - a reference to the current node, which will be updated
 *                     with "base".
 *      desc - a description of the sequence being indexed.
 *      pos - the current position in sequence being indexed.
 */
int process_base(bp_t base, gtree_t **cur_node_ref, long pos, char *desc);

/**
 * build gtree index, assign to "gtree_root".
 *
 * ASSUME:
 *      - no description strings are greater than MAX_DESC_LEN chars long.
 *
 * @args:
 *      ix_file - FASTA file to 
 *      gtree_root - pointer to to assign gtree_root to
 *      desc_strings - pointer to an array of strings;
 *                     assign description strings to for free later
 *      n_descs - number of description strings in desc_strings
 *
 * @return:
 *      0        on success
 *      errcode  otherwise
 *
 */
int build_gtree( char *ix_file,
                 gtree_t **gtree_root, 
                 char ***desc_strings,
                 unsigned int *n_descs );

/**
 * Construct an index from a basic FASTA file
 *
 * @args:
 *      ref_filename - filename for the reference sequence
 *
 */
ix_t *build_ix_from_ref_seq( char *ref_filename );

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
