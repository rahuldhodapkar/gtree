#!/usr/bin/bash
## basic_align.sh
#
# Create a basic alignment scenario from the lambda phage data.
# Performs all necessary cleaning / staging / etc of files and can be
# tuned as necessary.
# 
# NOTES:
#   should be run from the base directory of the git repository

# set -x
set -e

echo "===== BEGIN BASIC ALIGNMENT SETUP ====="

# clean FASTA references
utils/rm-fa-newlines -i data/lambda_virus.fa -o lambda_virus_cleaned.fa

# build index
./gtree ix build -r lambda_virus_cleaned.fa -o lambda_virus_cleaned.gt

# prune index
./gtree ix prune -ix lambda_virus_cleaned.gt -o lambda_virus_cleaned.msk.gt

echo "===== END BASIC ALIGNMENT SETUP ====="
cat <<HERE

Your alignment setup is now complete, and all related files can be found
under the "./tmp" directory

./gtree aln -r lambda_virus_cleaned -ix lambda_virus_cleaned.gt -i \\
        data/longreads.fq

or 

./gtree aln -r lambda_virus_cleaned -ix lambda_virus_cleaned.msk.gt \\
        -rl GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCATAACTTAATGTTTTTATTTAAAATACCCTCTGAAAAGAAAGGAAACGACAGGTGCTGAAAGCGAGGCTTTTTGGCCTCT

HERE

