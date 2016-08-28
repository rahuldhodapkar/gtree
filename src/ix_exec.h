#ifndef IX_EXEC_H
#define IX_EXEC_H

/** ix_exec.h
 * methods for doing building, masking, and otherwise manipulate gtree indexes.
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
int gtree_ix(int argc, char *argv[]);

#endif
