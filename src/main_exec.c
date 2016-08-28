/** main.c
 * core driver handling actual alignment
 */

#include "ix_exec.h"
#include "aln_exec.h"

#include "consts.h"
#include "types.h"

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define GTREE_HELP_MESSAGE \
"Usage: gtree [func]\n"\
"    aln       core alignment algorithms\n"\
"    ix        index build, manipulation and info utilities\n"\
"\n"

int main(int argc, char *argv[]) {
    int rcode = 0;

    if (argc <= 1) {
        printf(GTREE_HELP_MESSAGE);
        exit(EXIT_SUCCESS);
    }


    if (strcmp(argv[1], "aln") == 0) {
        rcode = gtree_aln(argc, argv);
    } 
    else if (strcmp(argv[1], "ix") == 0) {
        rcode = gtree_ix(argc, argv);
    }
    else if (strcmp(argv[1], "-h") == 0) {
        printf(GTREE_HELP_MESSAGE);
        exit(EXIT_SUCCESS);
    }
    else {
        printf(GTREE_HELP_MESSAGE);
        exit(EXIT_FAILURE);
    }

    return rcode;
}

/*
int validate_args(args_t *args) {
    if (args->exec_mode < 0) {
        printf("ERROR: no execution mode chosen, use build or align\n");
        exit(EXIT_FAILURE);
    }
    return 0;
}

int prune_ix(args_t *args) {

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

    return 0;

}

int build_ix(args_t *args) {

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

int mask_ix(args_t *args) {

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

int align_seq(args_t *args) {

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

int main(int argc, char *argv[]) {
    args_t args;

    args.exec_mode = -1;        // invalid initial exec mode
    args.verbosity = VERBOSITY_LEVEL_QUIET;
    args.ref_fasta_fn = 
        args.ix_fn = 
        args.out_fn = 
        args.in_fn = 
        args.in_fn2 = NULL;
    args.out_format = OUTPUT_FORMAT_SAM;

    args.exec_mode = -1;
 
    if (argc <= 1) {
        printf(HELPMSG);
        exit(EXIT_SUCCESS);
    }

    int i = 1;
    while (i < argc) {
        if (strcmp("-h", argv[i]) == 0) {
            printf(HELPMSG);
            exit(EXIT_SUCCESS);
        } else if (strcmp("-b", argv[i]) == 0) {
            args.exec_mode = EXEC_MODE_BUILD_IX;
        } else if (strcmp("-a", argv[i]) == 0) {
            args.exec_mode = EXEC_MODE_ALIGN;
        } else if (strcmp("-maskix", argv[i]) == 0) {
            args.exec_mode = EXEC_MODE_MASK_IX;
        } else if (strcmp("-pruneix", argv[i]) == 0) {
            args.exec_mode = EXEC_MODE_PRUNE_IX;
        } else if (strcmp("-v", argv[i]) == 0) {
            args.verbosity = VERBOSITY_LEVEL_DEBUG;
        } else if (strcmp("-r", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no ref sequence passed with '-r'\n");
                exit(EXIT_FAILURE);
            }

            args.ref_fasta_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-ix", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no index filename passed with '-ix'\n");
                exit(EXIT_FAILURE);
            }

            args.ix_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-o", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no output file passed with '-o'\n");
                exit(EXIT_FAILURE);
            }

            args.out_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-of", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no output format passed with '-of'\n");
                exit(EXIT_FAILURE);
            }

            if (strcmp(argv[i+1], "SAM") == 0) {
                args.out_format = OUTPUT_FORMAT_SAM;
            } else if (strcmp(argv[i+1], "BAM") == 0) {
                args.out_format = OUTPUT_FORMAT_BAM;
            } else {
                printf("ERROR: invalid output format %s passed, " 
                       "choose 'SAM' or 'BAM'\n", argv[i+1]);
                exit(EXIT_FAILURE);
            }

            i++;
        } else if (strcmp("-i", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                printf("ERROR: no input sequence passed with '-r'\n");
                exit(EXIT_FAILURE);
            }

            args.in_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-pe", argv[i]) == 0) {
            if ( i + 2 >= argc ) {
                printf("ERROR: too few input seqs passed with '-pe'\n");
                exit(EXIT_FAILURE);
            }

            args.in_fn = argv[i+1]; 
            args.in_fn2 = argv[i+2]; 
            i += 2;
        }
        i++;
    }

    validate_args(&args);

    if (args.exec_mode == EXEC_MODE_BUILD_IX) {
        build_ix(&args);
    } else if (args.exec_mode == EXEC_MODE_ALIGN) {
        align_seq(&args);
    } else if (args.exec_mode == EXEC_MODE_MASK_IX) {
        mask_ix(&args); 
    } else if (args.exec_mode == EXEC_MODE_PRUNE_IX) {
        prune_ix(&args); 
    } else {
        printf("ERROR: unknown exec_mode option '%d', passed\n", args.exec_mode);
        exit(EXIT_FAILURE);
    }

    printf("finished running!\n");
    return 0;
}

*/
