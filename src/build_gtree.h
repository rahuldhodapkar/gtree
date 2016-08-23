#ifndef BUILD_GTREE_H
#define BUILD_GTREE_H

/** build_gtree.h
 * encapsulate methods for building a gtree from a reference genome file.
 * 
 * Currently, only FASTA-formatted reference genomes in the style of those
 * provided by the UCSC genome browser portal are supported.
 */

#include "index.h"
#include "gtree.h"

#include "types.h"
#include "consts.h"

#include <stdlib.h>
#include <stdio.h>

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
 *      ix_file - FASTA file to build index from
 *      ix - gtree index to build
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
 * tests a gtree index for uniqueness against a reference FASTA file by
 * continuing a build, but not allocating new nodes for new sequences.
 * 
 * allows for building gene-specific index targets with small memory footprint
 * that are still selective against the genome as a whole!
 *
 * ***NOTE*** this function will modify the index "ix" passed in.
 *
 * @args:
 *      mask_file - FASTA file to run against index
 *      ix - gtree index to test
 *
 * @return:
 *      0        on success
 *      errcode  otherwise
 */
int mask_gtree( char *mask_file, ix_t *ix );

/**
 * Construct an index from a basic FASTA file
 *
 * @args:
 *      ref_filename - filename for the reference sequence
 *
 */
ix_t *build_ix_from_ref_seq( char *ref_filename );

#endif
