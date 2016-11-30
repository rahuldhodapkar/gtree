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

use Test::More tests => 8;
use POSIX qw(mkfifo);

my @test_files = qw/.ti0.fa .ti1.fa \
                    .ti0.far .ti1.far \
                    .ti0.refix \
                    .ti1.refix \
                    .ti1 .ti2 \
                    .to0.ix .to1.ix \
                    .reads.fq .same_aligned.sam \
                    .reads_same_seed.fq /;

my $out;

my $utils_root_path = 'utils/';

####################################################
## GENERATE INPUT FILES
####################################################

open(FILE, '>', '.ti0.far') or die $!;
# 150 bp FASTA ref
print FILE <<"HERE";
>chr1
GGGGGAAAA
GGGGGGGGG
GGGGGGGGG
>chr2
GGGGGGGGG
GGGGGGGGG
GGAAAAGGG
HERE
close(FILE);

`cat .ti0.far | $utils_root_path/rm-fa-newlines > .ti0.fa`;

open(FILE, '>', '.ti1.far') or die $!;
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
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
HERE
close(FILE);

`cat .ti1.far | $utils_root_path/rm-fa-newlines > .ti1.fa`;

open(FILE, '>', '.reads.fq') or die $!;
# 150 bp FASTA ref
print FILE <<"HERE";
\@r1
GGGGGAAAA
+
~~~~~~~~~
\@r2
GGGGGGGGG
+ 
~~~~~~~~~
HERE
close(FILE);

open(FILE, '>', '.reads_same_seed.fq') or die $!;
# 150 bp FASTA ref
print FILE <<"HERE";
\@r1
ATCACCCCACCCTATATAACCCCACACTTGGTTTAAAACACAAACCCATA
+
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HERE
close(FILE);

####################################################
## TEST REFERENCE SEEK
####################################################

$out = `./gtree ix build -r .ti0.fa -o .to0.ix`;
ok( $out !~ /ERROR/, 'error state detected' );

$out = `./gtree aln -r .ti0 -ix .to0.ix -rl GGGGGG`;
ok( $out !~ /ERROR/, 'error state detected' );

$out = `./gtree aln -r .ti0 -ix .to0.ix -i .reads.fq`;
ok( $out !~ /ERROR/, 'error state detected' );

####################################################
## TEST SAME_SEED_IGNORE_DIST
####################################################

$out = `./gtree ix build -r .ti1.fa -o .to1.ix`;
ok( $out !~ /ERROR/, 'error state detected' );

$out = `./gtree aln -r .ti1 -ix .to1.ix -i .reads_same_seed.fq > .same_aligned.sam`;
ok( $out !~ /ERROR/, 'error state detected' );

$out = `grep 'r1' .same_aligned.sam | wc -l`;
$out =~ s/^\s*(\d+)\s*$/$1/;
is($out, "1", "simple alignment gets expected number of matches");

$out = `cat .same_aligned.sam`;
ok( $out =~ /chr1\t251/, "simple alignment matches correct position");

# clean up test files
unlink( @test_files );

ok( ! -e @test_files, 'fully cleaned up' );

