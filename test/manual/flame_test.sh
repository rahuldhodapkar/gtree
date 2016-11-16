#!/usr/bin/bash
## flame_test.sh
#
# run flame test analysis of `gtree` alignment software for profiling.
# uses `dtrace` to obtain stack traces and FlameGraph to consolidate
# and plot them.
#
# REQUIRES:
#   system which supports Linux perf
#   FlameGraph perl scripts in $PATH

set -x
set -e

# configure lambda macrophage test
bash test/manual/basic_align.sh

# run longreads alignment
./gtree aln -r lambda_virus_cleaned -ix lambda_virus_cleaned.gt -i data/longreads.fq > lambda_aligned.sam &

# 1. collect stack traces from longreads alignment
sudo perf record -F 99 -p $! -g -- sleep 10
perf script > out.perf
#dtrace -x ustackframes=100 -n "profile-97 /pid == $!/ { @[ustack()] = count(); } tick-60s { exit(0); }" -o gtree.user_stacks

# 2. fold stacks
stackcollapse-perf.pl out.perf > out.folded
#stackcollapse.pl gtree.user_stacks > gtree.user_folded

# 3. render flamegraph

flamegraph.pl out.folded > gtree.perf.svg
#flamegraph.pl gtree.user_folded > gtree.perf.svg

