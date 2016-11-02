/** aln_exec.c
 * core driver handling actual alignment
 */

#include "aln_exec.h"

#include "debug.h"
#include "index.h"
#include "ref.h"
#include "align.h"

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
"    -rl [read_fasta]           literal FASTA single read to align\n"\
"                               ** feature not yet supported **\n"\
"    -pe [path1] [path2]        input FASTQ files for paired-end alignment\n"\
"                               ** feature not yet supported **\n"\
"    -h                         print this message and quit\n"\
"\n"

read_t *_init_read() {
    read_t *ret = malloc(sizeof(read_t));

    ret->malloc_len = 0;
    ret->seq = NULL;
    ret->read_seq = NULL;
    ret->phred = NULL;
    ret->len = 0;

    return ret;
}

int _destroy_read(read_t *read) {

    if (read->seq != NULL)       free(read->seq);
    if (read->read_seq != NULL)  free(read->read_seq);
    if (read->phred != NULL)     free(read->phred);
    free(read);

    return 0;
}


int validate_aln_args(args_t *args) {
    return 0;
}

int aln_simple(args_t *args) {

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    struct timeval tval_before, tval_after, tval_result;
#endif

    ix_t *ix;
    ref_t *ref;

    /////////////////////////////////////////////////////////////////////////
    //  LOAD INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Loading index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    ix = deserialize_ix(args->ix_fn);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    print_ix_info(ix);
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  LOAD REFERENCE
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Loading reference...\n");
    gettimeofday(&tval_before, NULL);
#endif

    ref = load_ref(args->ref_fasta_fn);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    print_ref_info(ref);
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  ALIGN READS
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Aligning reads...\n");
    gettimeofday(&tval_before, NULL);
#endif

    FILE *reads_file = fopen(args->in_fn, "r");
    read_t *read = _init_read();
    alnres_t aln;

    print_sequence_headers(ref);
    while (get_next_read(reads_file, read)) {
        DEBUG("read length - [%d]\n", read->len);
        align_single_read(read, ix, ref, &aln);
    }

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    _destroy_read(read);

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Destroying index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    destroy_ix(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY REFERENCE
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Destroying ref...\n");
    gettimeofday(&tval_before, NULL);
#endif

    destroy_ref(ref);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    return 0;
}

/**
 * align only a single seed sequence against the reference.
 */
int aln_seed_seq(args_t *args) {
    INFO("Beginning single read literal alignment\n");

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    struct timeval tval_before, tval_after, tval_result;
#endif

    ix_t *ix;
    ref_t *ref;

    /////////////////////////////////////////////////////////////////////////
    //  LOAD INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Loading index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    ix = deserialize_ix(args->ix_fn);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    print_ix_info(ix);
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  LOAD REFERENCE
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Loading reference...\n");
    gettimeofday(&tval_before, NULL);
#endif

    ref = load_ref(args->ref_fasta_fn);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    print_ref_info(ref);
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  ALIGN READS
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Aligning reads...\n");
    gettimeofday(&tval_before, NULL);
#endif

    read_t *read = _init_read();
    alnres_t aln;

    size_t read_len = strlen(args->in_read_literal);

    read->template_id[0] = '\0'; // empty string for anon template
    read->read_seq = malloc(sizeof(char) * read_len + 1);
    read->seq = malloc(sizeof(bp_t) * read_len);
    read->phred = malloc(sizeof(char) * read_len + 1);
    read->len = read_len;
    read->malloc_len = read_len;

    int i;
    for (i = 0; i  < read_len; i++) {
        read->read_seq[i] = args->in_read_literal[i];
        read->seq[i] = char_to_bp(args->in_read_literal[i]);
        read->phred[i] = '~';   // max quality
    }
    read->phred[read_len]    = '\0';   // NULL terminate strings
    read->read_seq[read_len] = '\0';

    align_single_read(read, ix, ref, &aln);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    _destroy_read(read);

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Destroying index...\n");
    gettimeofday(&tval_before, NULL);
#endif

    destroy_ix(ix);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

    /////////////////////////////////////////////////////////////////////////
    //  DESTROY REFERENCE
    /////////////////////////////////////////////////////////////////////////

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    INFO("Destroying ref...\n");
    gettimeofday(&tval_before, NULL);
#endif

    destroy_ref(ref);

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    INFO("Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);
#endif

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
        args.in_read_literal =
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
        } else if (strcmp("-rl", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                DIE("Invalid command line option usage - %s\n", "'-rl'");
            }

            args.in_read_literal = argv[i+1]; 
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
        } else if (strcmp("-i", argv[i]) == 0) {
            if ( i + 1 >= argc ) {
                DIE("Invalid command line option usage - %s\n", "'-i'");
            }

            args.in_fn = argv[i+1]; 
            i++;
        } else if (strcmp("-pe", argv[i]) == 0) {
            if ( i + 2 >= argc ) {
                DIE("Invalid command line option usage - %s\n", "'-pe'");
            }

            args.in_fn = argv[i+1]; 
            args.in_fn2 = argv[i+2];
            i += 2;
        }
        i++;
    }

    validate_aln_args(&args);

    if (args.in_read_literal == NULL) {
        aln_simple(&args);
    }
    else {
        aln_seed_seq(&args);
    }

    return 0;
}
