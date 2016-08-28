#!/bin/sh

# ===== bwarun.sh 
#
# simple test script to run for setting up a basic bwa alignment with the
# material provided in the git repository.
#
# assumes the following relative directory structure:
#       bwarun.sh
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

rm -rf ix out stat
mkdir -p ix out stat

bwa index -p ix/lambda_virus ref/lambda_virus.fa

# single
bwa mem ix/lambda_virus reads/reads_1.fq \
		 >out/reads_1.sam 2>stat/reads_1.runstats

# longreads
bwa mem ix/lambda_virus reads/longreads.fq \
		 >out/longreads.sam 2>stat/longreads.runstats

# paired-end
bwa mem ix/lambda_virus \
			reads/reads_1.fq reads/reads_2.fq \
		 	>out/reads-pe.sam 2>stat/reads-pe.runstats

# do post-processing to get BAM files for samtools
for f in out/*.sam
do
    b=$(basename "$f" .sam)
    samtools view -bT ref/lambda_virus.fa out/$b.sam | samtools sort - out/$b.sorted
    samtools flagstat out/$b.sorted.bam >stat/$b.flagstat
done

