CC=gcc
CFLAGS=-Wall -pedantic -std=c99
DEBUG=-ggdb

gtree: src/gtree.c
	$(CC) $(CFLAGS) $^ -o $@

gtree-debug: src/gtree.c
	$(CC) $(CFLAGS) $(DEBU) $^ -o $@

.PHONY: clean debug

CLEAN_TARGETS=gtree gtree-debug *.dSYM
CLEAN_FLAGS=-rf

clean:
	rm $(CLEAN_FLAGS) $(CLEAN_TARGETS)

debug: gtree-debug
