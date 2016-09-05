CC=gcc
CFLAGS=-Wall -pedantic -std=c99 -DTRACE -D_BSD_SOURCE \
		-fno-common
DEBUG=-ggdb

all: gtree

# Target-specific variable values for "debug"
debug-symbols: CFLAGS += -ggdb
debug-symbols: gtree

gtree: src/main_exec.c gtree.o build_gtree.o index.o \
                       ref.o align.o \
                       ix_exec.o aln_exec.o \
					   debug.o
	$(CC) $(CFLAGS) $^ -o $@

debug.o: src/debug.c
	$(CC) $(CFLAGS) $^ -c -o $@

ix_exec.o: src/ix_exec.c
	$(CC) $(CFLAGS) $^ -c -o $@

aln_exec.o: src/aln_exec.c
	$(CC) $(CFLAGS) $^ -c -o $@

ref.o: src/ref.c
	$(CC) $(CFLAGS) $^ -c -o $@

align.o: src/align.c
	$(CC) $(CFLAGS) $^ -c -o $@

build_gtree.o: src/build_gtree.c
	$(CC) $(CFLAGS) $^ -c -o $@

gtree.o: src/gtree.c
	$(CC) $(CFLAGS) $^ -c -o $@

index.o: src/index.c
	$(CC) $(CFLAGS) $^ -c -o $@

.PHONY: clean test test-all

CLEAN_TARGETS=gtree gtree-debug *.dSYM *.o
CLEAN_FLAGS=-rf

clean:
	rm $(CLEAN_FLAGS) $(CLEAN_TARGETS)

test:
	prove test

test-all:
	prove test test/heavy

