#!/usr/bin/perl -w

# === core.t
#
# core test suite to ensure that things are working properly at a basic level.
#
# NOTE:
#   qw(.ti0 .ti1 .to0)
#
# @author rahuldhodapkar
# @version 2016-08-27
# @copyright Rahul Dhodapkar

use strict;
use warnings;

use Test::Simple tests => 10;
use POSIX qw(mkfifo);

my @test_files = qw/.ti0 .ti1 .ti2 \
                    .to0 .to1 .to2 \
                    .to0.prn .to1.prn .to2.prn \
                    .to0.msk .to1.prn .to2.prn \
                    .to0.msk.prn .to1.msk.prn .to2.msk.prn /;
my $out;

####################################################
## GENERATE INPUT FILES
####################################################

open(FILE, '>', '.ti0') or die $!;
# 150 bp FASTA ref
print FILE <<"HERE";
>chr1
gggggggggggggggggggggggggggggggg
HERE
close(FILE);

open(FILE, '>', '.ti1') or die $!;
print FILE <<"HERE";
>chr1
ggggggggggggggggggggggggggggggggg
HERE
close(FILE);

open(FILE, '>', '.ti2') or die $!;
print FILE <<"HERE";
>chr1
gggggggggggggggggggggggggggggggga
HERE
close(FILE);

####################################################
## TEST INDEX BUILD
####################################################

$out = `./gtree ix build -r .ti1 -o .to0`;
ok( $out =~ /nodes: 32/, 'build single index window' );

$out = `./gtree ix build -r .ti1 -o .to1`;
ok( $out =~ /nodes: 32/, 'build index window with duplicate nodes' );

$out = `./gtree ix build -r .ti2 -o .to2`;
ok( $out =~ /nodes: 33/, 'build index with branching' );

####################################################
## TEST INDEX LOAD
####################################################

$out = `./gtree ix stat -n -ix .to0`;
ok( $out =~ /nodes: 32/, 'stat single index window' );

$out = `./gtree ix stat -n -ix .to1`;
ok( $out =~ /nodes: 32/, 'stat index window with duplicate nodes' );

$out = `./gtree ix stat -n -ix .to2`;
ok( $out =~ /nodes: 33/, 'stat index with branching' );

####################################################
## TEST INDEX PRUNE
####################################################

$out = `./gtree ix prune -ix .to0 -o .to0.prn`;
ok( $out =~ /nodes: 2/, 'prune single index window' );

$out = `./gtree ix prune -ix .to1 -o .to1.prn`;
ok( $out =~ /nodes: 2/, 'prune index window with duplicate nodes' );

$out = `./gtree ix prune -ix .to2 -o .to2.prn`;
ok( $out =~ /nodes: 33/, 'prune index with branching' );

####################################################
## TEST INDEX LOAD
####################################################


# clean up test files
unlink( @test_files );

ok( ! -e @test_files, 'fully cleaned up' );

