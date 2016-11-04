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
        DIE("Invalid command line options, no exec mode chosen;"
                                          " use build or align\n", 0);
    }
    return 0;
}

int ix_prune(args_t *args) {
    // use POSIX functions for timing harness

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    struct timeval tval_before, tval_after, tval_result;
#endif

    ix_t *ix;

    /////////////////////////////////////////////////////////////////////////
    //  LOAD INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Loading index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    ix = deserialize_ix(args->ix_fn);
    print_ix_info(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  PRUNE INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Pruning index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    prune_gtree(ix->root);
    print_ix_info(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif
    
    /////////////////////////////////////////////////////////////////////////
    //  SERIALIZE INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Serializing...\n");
    gettimeofday(&tval_before, NULL);
#endif

    serialize_ix(ix, args->out_fn);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Serializing done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Destroying built index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    destroy_ix(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Destroying done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    return 0;
}

int ix_build(args_t *args) {

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    struct timeval tval_before, tval_after, tval_result;
#endif
    ix_t *ix;

    /////////////////////////////////////////////////////////////////////////
    //  BUILD INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Building...\n");
    gettimeofday(&tval_before, NULL);
#endif

    ix = build_ix_from_ref_seq(args->ref_fasta_fn);
    print_ix_info(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Building done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  SERIALIZE INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Serializing...\n");
    gettimeofday(&tval_before, NULL);
#endif

    serialize_ix(ix, args->out_fn);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Serializing done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Destroying built index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    destroy_ix(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Destroying done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    return 0;
}

int ix_mask(args_t *args) {

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    struct timeval tval_before, tval_after, tval_result;
#endif
    ix_t *ix;

    /////////////////////////////////////////////////////////////////////////
    //  LOAD INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Loading index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    ix = deserialize_ix(args->ix_fn);
    print_ix_info(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  MASK INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Pruning index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    mask_gtree(args->ref_fasta_fn, ix);
    print_ix_info(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  SERIALIZE INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Serializing...\n");
    gettimeofday(&tval_before, NULL);
#endif

    serialize_ix(ix, args->out_fn);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Serializing done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Destroying built index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    destroy_ix(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Destroying done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    return 0;
}

int ix_stat(args_t *args) {

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    struct timeval tval_before, tval_after, tval_result;
#endif
    ix_t *ix;

    /////////////////////////////////////////////////////////////////////////
    //  LOAD INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Loading index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    ix = deserialize_ix(args->ix_fn);
    print_ix_info(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Destroying loaded index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    destroy_ix(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Destroying done in %ld.%06ld secs\n\n", 
                                        (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

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
        DIE("Invalid command line options\n", 0);
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
                DIE("Invalid command line option usage - %s\n", "'-r'");
            }

            args.ref_fasta_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-ix", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                DIE("Invalid command line option usage - %s\n", "'-ix'");
            }

            args.ix_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-o", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                DIE("Invalid command line option usage - %s\n", "'-o'");
            }

            args.out_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-of", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                DIE("Invalid command line option usage - %s\n", "'-of'");
            }

            if (strcmp(argv[i+1], "SAM") == 0) {
                args.out_format = OUTPUT_FORMAT_SAM;
            } else if (strcmp(argv[i+1], "BAM") == 0) {
                args.out_format = OUTPUT_FORMAT_BAM;
            } else {
                DIE("Invalid output format %s passed, "
                       "choose 'SAM' or 'BAM'\n", argv[i+1]);
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
        DIE("Invalid exec_mode option '%d' passed\n", args.exec_mode);
    }

    INFO("finished running!\n");
    return 0;
}

