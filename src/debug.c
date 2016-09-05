
/** debug.c
 * define all constants to be used in aligner
 */

#include "debug.h"

#include <stdlib.h>
#include <execinfo.h>
#include <stdio.h>

void die(char *errmsg, int errcode) {
    void *bt_array[DEBUG_BT_BUFFER_SIZE];
    size_t n_frames;
    char **frames;

    n_frames = backtrace(bt_array, DEBUG_BT_BUFFER_SIZE);
    frames = backtrace_symbols(bt_array, n_frames);

    fprintf(stderr, "%s\n", errmsg);
    int i;
    for (i = 0; i < n_frames; i++) {
        fprintf(stderr, "%s\n", frames[i]);
    }
    free(frames);

    exit(errcode);
}
