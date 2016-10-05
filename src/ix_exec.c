/** ix_exec.c
 * core driver for index build, stat, and manipulation
 */

#include "ix_exec.h"

#include "consts.h"
#include "types.h"
#include "debug.h"

// workhorse functions
#include "gtree.h"
#include "build_gtree.h"

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define GTREE_IX_HELP_MESSAGE \
"Usage: gtree ix [options]\n"\
"\n"\
"# INDEX BUILD \n"\
"    Usage: gtree ix build\n"\
"        -r [path]                 reference sequence FASTA filename\n"\
"        -o [path]                 prebuilt index path for alignment\n"\
"\n"\
"# INDEX MASK \n"\
"    Usage: gtree ix mask\n"\
"        -ix [path]                pre-built index to be masked for\n"\
"                                  selectivity\n"\
"        -r [path]                 reference sequence FASTA filename\n"\
"        -o [path]                 prebuilt index path for alignment\n"\
"\n"\
"# INDEX PRUNE \n"\
"    Usage: gtree ix prune\n"\
"        -ix [path]                pre-built index to be masked for\n"\
"                                  selectivity\n"\
"        -o [path]                 prebuilt index path for alignment\n"\
"\n"\
"# INDEX STATS \n"\
"    Usage: gtree ix stat\n"\
"        -ix [path]                pre-built index to be masked printed\n"\
"        -n                        print # of nodes in gtree to report on\n"\
"\n"\
"\n"

int validate_ix_args(args_t *args) {
    if (args->exec_mode < 0) {
        printf("ERROR: no execution mode chosen, use build or align\n");
        DIE("Invalid command line options, no exec mode chosen", 0);
    }
    return 0;
}

int ix_prune(args_t *args) {
    // use POSIX functions for timing harness
    struct timeval tval_before, tval_after, tval_result;
    ix_t *ix;

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
    //  PRUNE INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Pruning index...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    prune_gtree(ix->root);
    print_ix_info(ix);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
    
    /////////////////////////////////////////////////////////////////////////
    //  SERIALIZE INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Serializing...\n");
    gettimeofday(&tval_before, NULL);

    // call to time
    serialize_ix(ix, args->out_fn);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Serializing done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Destroying built index...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    destroy_ix(ix);
    // 
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Destroying done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    return 0;
}

int ix_build(args_t *args) {

    // use POSIX functions for timing harness
    struct timeval tval_before, tval_after, tval_result;
    ix_t *ix;

    /////////////////////////////////////////////////////////////////////////
    //  BUILD INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Building...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    ix = build_ix_from_ref_seq(args->ref_fasta_fn);
    print_ix_info(ix);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Building done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    /////////////////////////////////////////////////////////////////////////
    //  SERIALIZE INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Serializing...\n");
    gettimeofday(&tval_before, NULL);

    // call to time
    serialize_ix(ix, args->out_fn);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Serializing done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Destroying built index...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    destroy_ix(ix);
    // 
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Destroying done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    return 0;
}

int ix_mask(args_t *args) {

    // use POSIX functions for timing harness
    struct timeval tval_before, tval_after, tval_result;
    ix_t *ix;

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
    //  MASK INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Pruning index...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    mask_gtree(args->ref_fasta_fn, ix);
    print_ix_info(ix);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    /////////////////////////////////////////////////////////////////////////
    //  SERIALIZE INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Serializing...\n");
    gettimeofday(&tval_before, NULL);

    // call to time
    serialize_ix(ix, args->out_fn);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Serializing done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Destroying built index...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    destroy_ix(ix);
    // 
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Destroying done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    return 0;
}

int ix_stat(args_t *args) {

    // use POSIX functions for timing harness
    struct timeval tval_before, tval_after, tval_result;
    ix_t *ix;

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
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Destroying loaded index...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    destroy_ix(ix);
    // 
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Destroying done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    return 0;
}

int gtree_ix(int argc, char *argv[]) {
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
        printf(GTREE_IX_HELP_MESSAGE);
        DIE("Invalid command line options", 0);
    }

    if (strcmp(argv[2], "build") == 0) {
        args.exec_mode = EXEC_MODE_IX_BUILD;
    } else if (strcmp(argv[2], "mask") == 0) {
        args.exec_mode = EXEC_MODE_IX_MASK;
    } else if (strcmp(argv[2], "prune") == 0) {
        args.exec_mode = EXEC_MODE_IX_PRUNE;
    } else if (strcmp(argv[2], "stat") == 0) {
        args.exec_mode = EXEC_MODE_IX_STAT;
    }

    int i = 3;
    while (i < argc) {
        if (strcmp("-h", argv[i]) == 0) {
            printf(GTREE_IX_HELP_MESSAGE);
            exit(EXIT_SUCCESS);
        } else if (strcmp("-v", argv[i]) == 0) {
            args.verbosity = VERBOSITY_LEVEL_DEBUG;
        } else if (strcmp("-r", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no ref sequence passed with '-r'\n");
                DIE("Invalid command line option usage - %s", "'-r'");
            }

            args.ref_fasta_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-ix", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no index filename passed with '-ix'\n");
                DIE("Invalid command line option usage - %s", "'-ix'");
            }

            args.ix_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-o", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no output file passed with '-o'\n");
                DIE("Invalid command line option usage - %s", "'-o'");
            }

            args.out_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-of", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no output format passed with '-of'\n");
                DIE("Invalid command line option usage - %s", "'-of'");
            }

            if (strcmp(argv[i+1], "SAM") == 0) {
                args.out_format = OUTPUT_FORMAT_SAM;
            } else if (strcmp(argv[i+1], "BAM") == 0) {
                args.out_format = OUTPUT_FORMAT_BAM;
            } else {
                printf("ERROR: invalid output format %s passed, " 
                       "choose 'SAM' or 'BAM'\n", argv[i+1]);
                DIE("Invalid output format %s passed", argv[i+1]);
            }

            i++;
        }
        i++;
    }

    validate_ix_args(&args);

    if (args.exec_mode == EXEC_MODE_IX_BUILD) {
        ix_build(&args);
    } else if (args.exec_mode == EXEC_MODE_IX_MASK) {
        ix_mask(&args); 
    } else if (args.exec_mode == EXEC_MODE_IX_PRUNE) {
        ix_prune(&args); 
    } else if (args.exec_mode == EXEC_MODE_IX_STAT) {
        ix_stat(&args);
    } else {
        printf("ERROR: unknown exec_mode option '%d', passed\n", args.exec_mode);
        DIE("Invalid exec_mode option '%d' passed", args.exec_mode);
    }

    printf("finished running!\n");
    return 0;
}

