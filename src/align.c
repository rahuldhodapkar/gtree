/** align.c
 * core driver handling actual alignment
 */
#include "gtree.h"
#include "build_gtree.h"

#include <time.h>
#include <sys/time.h>

int main(int argc, char *argv[]) {
    char *fgenome = NULL, *ix_fname = NULL;

    if (argc <= 1) {
        printf("Usage: gtree [genome_to_index] [ix_fname]\n");
    }

    if (argc > 1) {
        printf("input file [%s] passed.\n", argv[1]);
        fgenome = argv[1];
    }

    if (argc > 2) {
        printf("reads file [%s] passed.\n", argv[2]);
        ix_fname = argv[2];
    }

    // use POSIX functions for timing harness
    struct timeval tval_before, tval_after, tval_result;
    ix_t *ix;

    /////////////////////////////////////////////////////////////////////////
    //  BUILD INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Building...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    ix = build_ix_from_ref_seq(fgenome);
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
    serialize_ix(ix, ix_fname);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Serializing done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
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
    printf("INFO: Destroying done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    /////////////////////////////////////////////////////////////////////////
    //  LOAD INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Loading index...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    ix = deserialize_ix(ix_fname);
    print_ix_info(ix);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Loading done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    /////////////////////////////////////////////////////////////////////////
    //  PRUNE INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Loading index...\n");
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
    //  DESTROY INDEX
    /////////////////////////////////////////////////////////////////////////
    printf("Destroying loaded index...\n");
    gettimeofday(&tval_before, NULL);
    // call to time
    destroy_ix(ix);
    // 
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf("INFO: Destroying done in %ld.%06ld secs\n\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    printf("finished running!\n");
    return 0;
}
