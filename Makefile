

CC=gcc
CFLAGS=-std=c99 -Wall -pedantic -O2


gtree: gtree.c
	$(CC) $(CFLAGS) $^ -o $@



