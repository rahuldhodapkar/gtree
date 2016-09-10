#ifndef REF_H
#define REF_H

/** ref.h
 * provide access interface to reference genome for local alignment.
 */

#include "types.h"
#include "consts.h"

/**
 * allocate a ref_t and initialize associated structures.
 *
 * @return:
 *      a pointer to a ref_t, must be free'd with destroy_ref() when done
 *      NULL if init failed
 */
ref_t *init_ref();

/**
 * destroy a ref_t created by init_ref() and all associated structures.
 *
 * @args:
 *      ref - a reference structure to destroy
 *
 * @return:
 *      0        on success
 *      errcode  otherwise
 */
int destroy_ref(ref_t *ref);

/**
 * Load a reference sequence for use by the aligner. Looks for:
 *
 *      <ref_file_base> + ".fa"     // reference FASTA file with no newlines
 *                                  // in sequence lines
 * 
 *      <ref_file_base> + ".refix"  // tab-delimited map of reference file and
 *                                  // fseek() positions by description
 *
 * and will automatically create <ref_file_base>.refix if it does not exist
 *
 * @args:
 *      ref_file_base - a path from which to load all associated reference files
 *
 * @return:
 *      a pointer to a loaded ref_t, must be free'd with destroy_ref() when done
 *      NULL if load failed
 */
ref_t *load_ref(char *ref_file_base);

/**
 * Print a description of a reference structure to STDOUT, including a list
 * of the description map and file seek positions
 *
 * @args:
 *      ref - a pointer to the ref structure to be printed
 */
 void print_ref_info(ref_t *ref);

/**
 * Copy bp string of "len" from "desc", "pos" to "refstr"
 *
 * if (desc, pos, len) results in a requested sequence longer than available,
 * the string written to "refstr" may be shorter than "len". In such a case,
 * the string will be padded with "N" wherever a reference character is not
 * available.
 *
 * @args:
 *      ref - reference from which to copy string
 *      desc - description of reference sequence to fetch
 *      pos - position in reference to fetch from
 *      len - length of string to copy
 *      refstr - bp string loc to copy to
 *      refstrlen - 
 *
 * @return:
 *      0        on success
 *      errcode  otherwise
 */
int refcpy( ref_t *ref,
            char *desc, unsigned long pos, unsigned int len, 
            bp_t *refstr, unsigned long *refstrlen);

#endif
