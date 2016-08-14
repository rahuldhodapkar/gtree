#!/bin/sh

# ===== b2run.sh 
#
# simple test script to run for setting up a basic bowtie2 alignment with the
# material provided in the git repository.
#
# assumes the following relative directory structure:
#       b2run.sh
#       ix (initially empty)
#       out (initially empty)
#       reads
#       |  longreads.fq
#       |  reads_1.fq
#       |  reads_2.fq
#       |  simulate.pl
#       ref 
#       |  lambda_virus.fa

set -x
set -e

mkdir -p ix
mkdir -p out
mkdir -p stat

bowtie2-build ref/lambda_virus.fa ix/lambda_virus

for f in reads/*.fq
do
    b=$(basename "$f" .txt)
    bowtie2 -x ix/lambda_virus -U $f -S out/$b.sam 2>stat/$b.runstats
done


