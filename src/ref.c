/** ref.c
 * provide access interface to reference genome for local alignment.
 */

#include "ref.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

ref_t *init_ref() {
    ref_t *ret = malloc(sizeof(ref_t));

    int i;
    for (i = 0; i < MAX_DESC_LEN; i++) {
        ret->cur_desc[i] = '\0';
    }
    ret->cur_pos = 0;
    ret->desc_map = NULL;

    return ret;
}

desc_loc_map_t *_init_desc_loc_map() {
    desc_loc_map_t *ret = malloc(sizeof(desc_loc_map_t));

    ret->ref_file = NULL;
    int i;
    for (i = 0; i < MAX_DESC_LEN; i++) {
        ret->desc[i] = '\0';
    }
    ret->pos = 0;
    ret->next = NULL;

    return ret;
}

int _destroy_desc_loc_map(desc_loc_map_t *desc_map) {
    if (desc_map == NULL) {
        return 0;
    }

    _destroy_desc_loc_map(desc_map->next);
    free(desc_map);

    return 0;
}

int destroy_ref(ref_t *ref) {
    _destroy_desc_loc_map(ref->desc_map);
    free(ref);
    return 0;
}

ref_t *load_ref(char *ref_file_base) {

    ref_t *ret = init_ref();

    char refix_filename[strlen(ref_file_base) + 10];
    char ref_filename[strlen(ref_file_base) + 10];

    strcpy(refix_filename, ref_file_base);
    strcat(refix_filename, ".refix");

    strcpy(ref_filename, ref_file_base);
    strcat(ref_filename, ".fa");

    if( access(refix_filename, F_OK) == -1 ) {
        // create the reference index file
        // ERROR and exit with stacktrace.
    }

    // read the reference into memory, opening necessary file pointers
    // FILE *refix = fopen(ref_file_base, "r");

    return ret;
}
