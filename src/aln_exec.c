/** aln_exec.c
 * core driver handling actual alignment
 */

#include "aln_exec.h"

#include "debug.h"
#include "index.h"
#include "ref.h"

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define GTREE_ALN_HELP_MESSAGE \
"Usage: gtree aln [options]\n"\
"    -v                         verbose mode\n"\
"    -r [path]                  reference sequence FASTA filename\n"\
"    -ix [path]                 prebuilt index path for alignment\n"\
"    -o [path]                  output file for alignment results\n"\
"    -of [format]               output file format - choose either SAM or BAM\n"\
"                               ** feature not yet supported **\n"\
"    -i [path]                  input FASTQ file, implies single reads\n"\
"                               ** feature not yet supported **\n"\
"    -pe [path1] [path2]        input FASTQ files for paired-end alignment\n"\
"                               ** feature not yet supported **\n"\
"    -h                         print this message and quit\n"\
"\n"

int validate_aln_args(args_t *args) {
    return 0;
}

int aln_simple(args_t *args) {
    // use POSIX functions for timing harness
    struct timeval tval_before, tval_after, tval_result;
    ix_t *ix;
    ref_t *ref;

    /////////////////////////////////////////////////////////////////////////
    //  LOAD INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Loading index...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    ix = deserialize_ix(args->ix_fn);
    print_ix_info(ix);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    /////////////////////////////////////////////////////////////////////////
    //  LOAD REFERENCE
    /////////////////////////////////////////////////////////////////////////
    printf("Loading reference...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    ref = load_ref(args->ref_fasta_fn);
    print_ref_info(ref);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
    
    return 0;
}

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

    args_t args;
    // define default args
    args.exec_mode = -1;
    args.verbosity = VERBOSITY_LEVEL_QUIET;
    args.ref_fasta_fn = 
        args.ix_fn = 
        args.out_fn = 
        args.in_fn = 
        args.in_fn2 = NULL;
    args.out_format = OUTPUT_FORMAT_SAM;

    if (argc <= 2) {
        printf(GTREE_ALN_HELP_MESSAGE);
        exit(EXIT_SUCCESS);
    }

    int i = 2;
    while (i < argc) {
        if (strcmp("-h", argv[i]) == 0) {
            printf(GTREE_ALN_HELP_MESSAGE);
            exit(EXIT_SUCCESS);
        } else if (strcmp("-v", argv[i]) == 0) {
            args.verbosity = VERBOSITY_LEVEL_DEBUG;
        } else if (strcmp("-r", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no ref sequence passed with '-r'\n");
                DIE("Invalid command line options");
            }

            args.ref_fasta_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-ix", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no index filename passed with '-ix'\n");
                DIE("Invalid command line options");
            }

            args.ix_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-o", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no output file passed with '-o'\n");
                DIE("Invalid command line options");
            }

            args.out_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-of", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no output format passed with '-of'\n");
                DIE("Invalid command line options");
            }

            if (strcmp(argv[i+1], "SAM") == 0) {
                args.out_format = OUTPUT_FORMAT_SAM;
            } else if (strcmp(argv[i+1], "BAM") == 0) {
                args.out_format = OUTPUT_FORMAT_BAM;
            } else {
                printf("ERROR: invalid output format %s passed, " 
                       "choose 'SAM' or 'BAM'\n", argv[i+1]);
                DIE("Invalid command line options");
            }

            i++;
        } else if (strcmp("-i", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no input file passed with '-i'\n");
                DIE("Invalid command line options");
            }

            args.in_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-pe", argv[i]) == 0) {
            if ( i + 2 >= argc ) {
                printf("ERROR: not enough inputs with '-pe'; 2 required\n");
                DIE("Invalid command line options");
            }

            args.in_fn = argv[i+1]; 
            args.in_fn2 = argv[i+2];
            i += 2;
        }
        i++;
    }

    validate_aln_args(&args);

    aln_simple(&args);

    return 0;
}
