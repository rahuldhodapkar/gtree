
/** debug.c
 * define all constants to be used in aligner
 */

#include "debug.h"

#include <stdlib.h>
#include <execinfo.h>
#include <stdio.h>
#include <stdarg.h>

void die(int errcode, char *msg, ...) {

    void *bt_array[DEBUG_BT_BUFFER_SIZE];
    size_t n_frames;
    char **frames;

    n_frames = backtrace(bt_array, DEBUG_BT_BUFFER_SIZE);
    frames = backtrace_symbols(bt_array, n_frames);

    va_list argptr;
    va_start(argptr, msg);
    vfprintf(stderr, msg, argptr);
    va_end(argptr);

    int i;
    for (i = 0; i < n_frames; i++) {
        fprintf(stderr, "%s\n", frames[i]);
    }
    free(frames);

    exit(errcode);
}
