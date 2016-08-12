CC=gcc
CFLAGS=-Wall -pedantic -std=c99
DEBUG=-ggdb

align: src/align.c gtree.o build_gtree.o
	$(CC) $(CFLAGS) $^ -o $@

build_gtree.o: src/build_gtree.c
	$(CC) $(CFLAGS) $^ -c -o $@

gtree.o: src/gtree.c
	$(CC) $(CFLAGS) $^ -c -o $@

.PHONY: clean debug

CLEAN_TARGETS=gtree gtree-debug *.dSYM *.o
CLEAN_FLAGS=-rf

clean:
	rm $(CLEAN_FLAGS) $(CLEAN_TARGETS)

debug: gtree-debug
