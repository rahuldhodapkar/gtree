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

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_INFO
#define INFO(...) printf("INFO: "); printf(__VA_ARGS__)
#else
#define INFO(...) do {} while (0)
#endif

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_WARN
#define WARN(...) printf("WARN: "); printf(__VA_ARGS__)
#else
#define WARN(...) do {} while (0)
#endif

#if VERBOSITY_LEVEL >= VERBOSITY_LEVEL_DEBUG
#define DEBUG(...) printf("DEBUG: "); printf(__VA_ARGS__)
#else
#define DEBUG(...) do {} while (0)
#endif

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
