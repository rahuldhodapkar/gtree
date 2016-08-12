CC=gcc
CFLAGS=-Wall -pedantic -std=c99
DEBUG=-ggdb

all: align

# Target-specific variable values for "debug"
debug: CFLAGS += -ggdb
debug: align

align: src/align.c gtree.o build_gtree.o
	$(CC) $(CFLAGS) $^ -o $@

build_gtree.o: src/build_gtree.c
	$(CC) $(CFLAGS) $^ -c -o $@

gtree.o: src/gtree.c
	$(CC) $(CFLAGS) $^ -c -o $@

.PHONY: clean

CLEAN_TARGETS=align align-debug *.dSYM *.o
CLEAN_FLAGS=-rf

clean:
	rm $(CLEAN_FLAGS) $(CLEAN_TARGETS)

