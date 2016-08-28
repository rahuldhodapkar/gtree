#!/bin/sh

# ===== bt2run.sh 
#
# simple test script to run for setting up a basic bowtie2 alignment with the
# material provided in the git repository.
#
# assumes the following relative directory structure:
#       bt2run.sh
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

bowtie2-build ref/lambda_virus.fa ix/lambda_virus

# single
bowtie2 -x ix/lambda_virus -U reads/reads_1.fq \
		-S out/reads_1.sam 2>stat/reads_1.runstats

# longreads
bowtie2 -x ix/lambda_virus -U reads/longreads.fq \
		-S out/longreads.sam 2>stat/longreads.runstats

# paired-end
bowtie2 -x ix/lambda_virus -1 reads/reads_1.fq -2 reads/reads_2.fq \
		-S out/reads-pe.sam 2>stat/reads-pe.runstats

# do post-processing to get BAM files for samtools
for f in out/*.sam
do
    b=$(basename "$f" .sam)
    samtools view -bT ref/lambda_virus.fa out/$b.sam | samtools sort - out/$b.sorted
    samtools flagstat out/$b.sorted.bam >stat/$b.flagstat
done

