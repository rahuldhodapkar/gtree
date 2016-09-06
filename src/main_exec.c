/** main.c
 * core driver handling actual alignment
 */

#include "ix_exec.h"
#include "aln_exec.h"

#include "consts.h"
#include "types.h"
#include "debug.h"

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
