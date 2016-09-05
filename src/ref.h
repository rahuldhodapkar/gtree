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
 * Load a reference sequence for use by the aligner.
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
 * Copy bp string of "len" from "desc", "pos" to "refstr"
 *
 * @args:
 *      desc - description of reference sequence to fetch
 *      pos - position in reference to fetch from
 *      len - length of string to copy
 *      refstr - bp string loc to copy to
 *
 * @return:
 *      0        on success
 *      errcode  otherwise
 */
int refcpy(char *desc, unsigned long pos, unsigned int len, bp_t *refstr);

#endif
