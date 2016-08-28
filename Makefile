CC=gcc
CFLAGS=-Wall -pedantic -std=c99 -DTRACE -D_BSD_SOURCE \
		-fno-common
DEBUG=-ggdb

all: gtree

# Target-specific variable values for "debug"
debug: CFLAGS += -ggdb
debug: gtree

gtree: src/main_exec.c gtree.o build_gtree.o index.o \
					   ix_exec.o aln_exec.o
	$(CC) $(CFLAGS) $^ -o $@

ix_exec.o: src/ix_exec.c
	$(CC) $(CFLAGS) $^ -c -o $@

aln_exec.o: src/aln_exec.c
	$(CC) $(CFLAGS) $^ -c -o $@

build_gtree.o: src/build_gtree.c
	$(CC) $(CFLAGS) $^ -c -o $@

gtree.o: src/gtree.c
	$(CC) $(CFLAGS) $^ -c -o $@

index.o: src/index.c
	$(CC) $(CFLAGS) $^ -c -o $@

.PHONY: clean test

CLEAN_TARGETS=gtree gtree-debug *.dSYM *.o
CLEAN_FLAGS=-rf

clean:
	rm $(CLEAN_FLAGS) $(CLEAN_TARGETS)

test:
	prove test

