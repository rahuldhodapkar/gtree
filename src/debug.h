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
#define DIE(msg, ...) die(ERRCODE_GENERAL_ERROR, msg, __VA_ARGS__)

/**
 * print error message to standard error and terminate process with errcode
 * 
 * @args:
 *      errcode - code to exit with
 *      msg - the message to exit with as a template string
 *      ... - must be at least of length 1, first argument passed as template 
 *            string to printf, remaining arguments passed as templating args
 */
void die(int errcode, char *msg, ...);

#endif
