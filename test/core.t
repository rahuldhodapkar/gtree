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

use Test::Simple tests => 4;
use POSIX qw(mkfifo);

my @test_files = qw/.ti0 .ti1 .ti2 .to0/;

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
my $out;

$out = `./gtree ix build -r .ti1 -o .to0`;
ok( $out =~ /nodes: 32/, 'single index window' );

$out = `./gtree ix build -r .ti1 -o .to0`;
ok( $out =~ /nodes: 32/, 'index window with duplicate nodes' );

$out = `./gtree ix build -r .ti2 -o .to0`;
ok( $out =~ /nodes: 33/, 'index with branching' );

# clean up test files
unlink( @test_files );

ok( ! -e @test_files, 'fully cleaned up' );

