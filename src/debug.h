#ifndef DEBUG_H
#define DEBUG_H

#include "consts.h"

/**
 * convenience function for die(errmsg, ERRCODE_GENERAL_ERROR)
 *
 * @args:
 *      errmsg - message to print to STDERR
 *
 */
#define DIE(errmsg) die(errmsg, ERRCODE_GENERAL_ERROR)

/**
 * print error message to standard error and terminate process with errcode
 * 
 * @args:
 *      errmsg - message to print to STDERR
 *      errcode - code to exit with
 */
void die(char *errmsg, int errcode);

#endif
