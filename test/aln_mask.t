#!/usr/bin/perl -w

# === aln1.t
#
# test core alignment functions
#
# @author rahuldhodapkar
# @version 2016-08-27
# @copyright Rahul Dhodapkar

use strict;
use warnings;

use Test::More tests => 11;
use POSIX qw(mkfifo);

my @test_files = qw/.ti0.fa .ti0.far \
                    .ti1.fa .ti1.far \
                    .ti0.refix \
                    .to0.ix \
                    .to0.msk.ix \
                    .reads.fq \
                    .reads_truly_selective.fq \
                    .aligned.sam /;

my $out;

my $utils_root_path = 'utils/';

####################################################
## GENERATE INPUT FILES
####################################################

open(FILE, '>', '.ti0.far') or die $!;
# 500 bp FASTA ref
# each line is 50 bp long
print FILE <<"HERE";
>chr1
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
TTATTATACCATAGGAGAGATACGTACGTGCTGACTGAGCATGCAGCTGC
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
HERE
close(FILE);

`cat .ti0.far | $utils_root_path/rm-fa-newlines > .ti0.fa`;

open(FILE, '>', '.ti1.far') or die $!;
# 500 bp FASTA ref
# each line is 50 bp long
print FILE <<"HERE";
>chr1
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
HERE
close(FILE);

`cat .ti1.far | $utils_root_path/rm-fa-newlines > .ti1.fa`;

open(FILE, '>', '.reads.fq') or die $!;
# 150 bp FASTA ref
print FILE <<"HERE";
\@r1
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
+
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HERE
close(FILE);

open(FILE, '>', '.reads_truly_selective.fq') or die $!;
# 150 bp FASTA ref
print FILE <<"HERE";
\@r2
TTATTATACCATAGGAGAGATACGTACGTGCTGACTGAGCATGCAGCTGC
+
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HERE
close(FILE);

####################################################
## BUILD INDEX AND MASKS
####################################################

$out = `./gtree ix build -r .ti0.fa -o .to0.ix`;
ok( $out !~ /ERROR/, 'error state detected' );

$out = `./gtree ix mask -ix .to0.ix -r .ti1.fa -o .to0.msk.ix`;
ok( $out !~ /ERROR/, 'error state detected' );

####################################################
## TEST ALIGNMENT
####################################################

#
# test unmasked index against locally and 
# globally selective sequences
#

$out = `./gtree aln -r .ti0 -ix .to0.ix -i .reads.fq > .aligned.sam`;
ok( $out !~ /ERROR/, 'error state detected' );

$out = `grep 'r1' .aligned.sam | wc -l`;
$out =~ s/^\s*(\d+)\s*$/$1/;
is($out, "1", "simple alignment gets expected number of matches");

$out = `./gtree aln -r .ti0 -ix .to0.ix -i .reads_truly_selective.fq > .aligned.sam`;
ok( $out !~ /ERROR/, 'error state detected' );

$out = `grep 'r2' .aligned.sam | wc -l`;
$out =~ s/^\s*(\d+)\s*$/$1/;
is($out, "1", "simple alignment gets expected number of matches");

#
# test masked index against locally selective sequence 
#

$out = `./gtree aln -r .ti0 -ix .to0.msk.ix -i .reads.fq > .aligned.sam`;
ok( $out !~ /ERROR/, 'error state detected' );

$out = `grep 'r1' .aligned.sam | wc -l`;
$out =~ s/^\s*(\d+)\s*$/$1/;
is($out, "0", "masked alignment gets expected number of matches");

#
# test masked index against globally selective sequence
#

$out = `./gtree aln -r .ti0 -ix .to0.msk.ix -i .reads_truly_selective.fq > .aligned.sam`;
ok( $out !~ /ERROR/, 'error state detected' );

$out = `grep 'r1' .aligned.sam | wc -l`;
$out =~ s/^\s*(\d+)\s*$/$1/;
is($out, "1", "masked alignment on selective sequence matches");

# clean up test files
unlink( @test_files );

ok( ! -e @test_files, 'fully cleaned up' );

