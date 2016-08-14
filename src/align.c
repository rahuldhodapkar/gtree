/** align.c
 * core driver handling actual alignment
 */
#include "gtree.h"
#include "build_gtree.h"

#include <time.h>
#include <sys/time.h>

int main(int argc, char *argv[]) {
    char *fgenome = NULL, *freads = NULL;

    if (argc <= 1) {
        printf("Usage: gtree [genome_to_index] [read_file]\n");
    }

    if (argc > 1) {
        printf("input file [%s] passed.\n", argv[1]);
        fgenome = argv[1];
    }

    if (argc > 2) {
        printf("reads file [%s] passed.\n", argv[2]);
        freads = argv[2];
    }

    // use POSIX functions for timing harness
    struct timeval tval_before, tval_after, tval_result;
    
    printf("Building...");
    gettimeofday(&tval_before, NULL);
    // call to time
    ix_t *ix = build_ix_from_ref_seq(fgenome);
    //
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf(" done in %ld.%06ld secs\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    printf("Cleaning up...");
    gettimeofday(&tval_before, NULL);
    // call to time
    destroy_ix(ix);
    // 
    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    printf(" done in %ld.%06ld secs\n", (long int)tval_result.tv_sec, 
                                        (long int)tval_result.tv_usec);

    printf("finished running!\n");
    return 0;
}
