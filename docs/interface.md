# Interface Specification for 'gtree'

## Abstract
This document defines an interface with which to interact with the 'gtree'
binary, developed for short read alignment.

The keywords "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD",
"SHOULD NOT", "RECOMMENDED", "MAY", and "OPTIONAL" in this document are to be
interpreted as in [rfc-2119](https://www.ietf.org/rfc/rfc2119.txt).

## Interface

Overview:

    # core alignment 
    gtree aln                               # core alignment algorithm, takes
                                            # a set of reads, either single or
                                            # paired-end, and aligns them
                                            # against an indexed reference. 

    # index stats
    gtree ix stat                           # print statistics about an index

    # index build
    gtree ix build                          # build an index that can match
                                            # against a reference sequence.

    # index manipulation + tuning
    gtree ix mask                           # mask an index for selectivity
                                            # against another sequence

    gtree ix prune                          # prune an index of nodes which do
                                            # not add additional information.

## User Stories
#### Unpaired read alignment against an entire reference genome "ref.fa"
1. Build a gtree index from the entire reference sequence

    gtree ix build -ref <ref.fa> -o <refix.gt>

2. Prune gtree index of unncessary nodes

    gtree ix prune -ix <refix.gt> -o <refix.pruned.gt>

3. Align reads against pruned index and reference sequence

    gtree aln -ix <refix.pruned.gt> -ref <ref.fa> -reads <reads.fq> \
                -o <aligned.sam>

#### Paired-end read alignment against an entire reference genome "ref.fa"
1. Build a gtree index from the entire reference sequence

    gtree ix build -ref <ref.fa> -o <refix.gt>

2. Prune gtree index of unncessary nodes

    gtree ix prune -ix <refix.gt> -o <refix.pruned.gt>

3. Align reads against pruned index and reference sequence

    gtree aln -ix <refix.pruned.gt> -ref <ref.fa> \
                -readsPE <reads1.fq> <reads2.fq> \
                -o <aligned.sam>

#### Get statistics on a gtree index

- Print number of nodes in a gtree index to STDOUT

    gtree ix stat -numNodes -ix <refix.gt>

- Print coverage of an index against a reference sequence to STDOUT
    
    gtree ix stat -cov -ref <ref.fa> -ix <refix.gt>

  ***NOTE*** sequence descriptions in <ref.fa> MUST match those used to build
             the index <refix.gt>.


