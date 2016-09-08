/** ref.c
 * provide access interface to reference genome for local alignment.
 */

#include "ref.h"

#include "debug.h"
#include "build_gtree.h"

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

    if (desc_map->ref_file != NULL) {
        fclose(desc_map->ref_file);
    }
    free(desc_map);

    return 0;
}

/**
 * ***NOTE*** assumes no newlines in sequence strings and single byte
 *            representation of all characters (ASCII).
 */
desc_loc_map_t *_build_refix_from_ref(char *ref_filename) {

    char c;
    unsigned long pos = 0;
    desc_loc_map_t *ret = NULL,
                   *swp = NULL,
                   *cur = NULL;
    FILE *ref = fopen(ref_filename, "r");

    char desc[MAX_DESC_LEN];

    while ( (c = fgetc(ref)) != EOF ) {

        pos++;

        if (c == '>') {
            read_desc(ref, desc);
            pos = pos + strlen(desc) + 1; // ***NOTE*** + 1 for '\n' consumed
            swp = cur;

            cur = _init_desc_loc_map();
            strcpy(cur->ref_filename, ref_filename);
            strcpy(cur->desc, desc);
            cur->ref_file = fopen(cur->ref_filename, "r");
            cur->pos = pos;

            if (swp != NULL) {
                swp->next = cur;
            }

            if (ret == NULL) {
                ret = cur;
            }
        }
    }

    return ret;
}

int _serialize_refix(desc_loc_map_t *refix, char *refix_filename) {

    FILE *out = fopen(refix_filename, "w+");

    desc_loc_map_t *cur = refix;

    int n_nodes = 0;
    while (cur != NULL) {
        n_nodes++;
        cur = cur->next;
    }

    fwrite(&n_nodes, sizeof(int), 1, out);

    int desc_len, ref_filename_len;
    cur = refix;

    while (cur != NULL) {
        ref_filename_len = strlen(cur->ref_filename) + 1;   // +1 for '\0'
        desc_len = strlen(cur->desc) + 1;

        fwrite(&ref_filename_len, sizeof(int), 1, out);
        fwrite(&(cur->ref_filename), sizeof(char), ref_filename_len, out);
        fwrite(&desc_len, sizeof(int), 1, out);
        fwrite(&(cur->desc), sizeof(char), desc_len, out);   
        fwrite(&(cur->pos), sizeof(unsigned long), 1, out);

        cur = cur->next;
    }

    fclose(out);

    return 0;
}

desc_loc_map_t *_deserialize_refix(char *refix_filename) {

    FILE *in = fopen(refix_filename, "r");

    desc_loc_map_t *ret = NULL, *cur = NULL, *swp = NULL;

    int n_nodes;
    fread(&n_nodes, sizeof(int), 1, in);


    int desc_len, ref_filename_len;
    int i;
    for (i = 0; i < n_nodes; i++) {

        swp = cur;

        cur = _init_desc_loc_map();

        fread(&ref_filename_len, sizeof(int), 1, in);
        fread(&(cur->ref_filename), sizeof(char), ref_filename_len, in);
        fread(&desc_len, sizeof(int), 1, in);
        fread(&(cur->desc), sizeof(char), desc_len, in);
        fread(&(cur->pos), sizeof(unsigned long), 1, in);

        cur->ref_file = fopen(cur->ref_filename, "r");

        if (swp != NULL) {
            swp->next = cur;
        }

        if (ret == NULL) {
            ret = cur;
        }
    }

    return ret;
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

    if ( access(ref_filename, F_OK) == -1 ) {
        printf("ERROR: unable to access reference file %s\n", ref_filename);
        DIE("Unable to access file");
    }

    desc_loc_map_t *refix;

    if( access( refix_filename, F_OK) == -1 ) {
        refix = _build_refix_from_ref(ref_filename);
        _serialize_refix(refix, refix_filename);
        _destroy_desc_loc_map(refix);
    }

    refix = _deserialize_refix(refix_filename);

    ret->desc_map = refix;

    return ret;
}

void print_ref_info(ref_t *ref) {
    desc_loc_map_t *cur = ref->desc_map;
    int n_entries = 0;

    while (cur != NULL) {
        printf("%s\t%s\t%lu\n", cur->ref_filename, cur->desc, cur->pos);
        cur = cur->next;
        n_entries++;
    }

    printf("Reference contains [%d] seqs\n", n_entries);
}

