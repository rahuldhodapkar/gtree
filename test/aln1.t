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

use Test::Simple tests => 5;
use POSIX qw(mkfifo);

my @test_files = qw/.ti0.fa .ti0.refix \
                    .ti1 .ti2 \
                    .to0.ix /;

my $out;

####################################################
## GENERATE INPUT FILES
####################################################

open(FILE, '>', '.ti0.fa') or die $!;
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

####################################################
## TEST REFERENCE SEEK
####################################################

$out = `./gtree ix build -r .ti0.fa -o .to0.ix`;
ok( $out !~ /ERROR/, 'error state detected' );

$out = `./gtree aln -r .ti0 -ix .to0.ix`;
ok( $out =~ /chr1\t6/, 'check "chr1" position' );
ok( $out =~ /chr2\t42/, 'check "chr2" position' );
ok( $out !~ /ERROR/, 'error state detected' );

# clean up test files
unlink( @test_files );

ok( ! -e @test_files, 'fully cleaned up' );

