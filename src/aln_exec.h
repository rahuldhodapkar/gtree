#ifndef ALN_EXEC_H
#define ALN_EXEC_H

/** aln_exec.h
 * methods for doing actual alignment
 */

/**
 * argument parsing function to call all other functions.
 *
 * @args:
 *      argc - the number of arguments passed to gtree
 *      argv - array of length "argc" containing command-line args
 *
 * @return:
 *      0        on success
 *      errcode  otherwise
 */
int gtree_aln(int argc, char *argv[]);

#endif
