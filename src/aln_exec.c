/** aln_exec.c
 * core driver handling actual alignment
 */

#include "aln_exec.h"

#include "debug.h"

#include <stdio.h>
#include <stdlib.h>

#define GTREE_ALN_HELP_MESSAGE \
"Usage: gtree aln [options]\n"\
"    -v                         verbose mode\n"\
"    -r [path]                  reference sequence FASTA filename\n"\
"    -ix [path]                 prebuilt index path for alignment\n"\
"    -o [path]                  output file for alignment results\n"\
"    -of [format]               output file format - choose either SAM or BAM\n"\
"                               ** feature not yet supported **\n"\
"    -i                         input FASTQ file, implies single reads\n"\
"                               ** feature not yet supported **\n"\
"    -pe [path1] [path2]        input FASTQ files for paired-end alignment\n"\
"                               ** feature not yet supported **\n"\
"    -h                         print this message and quit\n"\
"\n"

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
int gtree_aln(int argc, char *argv[]) {
    if (argc <= 2) {
        printf(GTREE_ALN_HELP_MESSAGE);
        exit(EXIT_SUCCESS);
    } else {
        printf("ERROR: 'gtree aln' functions not yet implemented\n");
        printf(GTREE_ALN_HELP_MESSAGE);
        DIE("Invalid command line options");
    }
    return 0;
}
