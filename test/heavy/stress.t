#!/usr/bin/perl -w

# === stress.t
#
# "stress" tests to ensure that nothing disastrous happens when running
# against slightly more real data sets.
#
# @REQUIRES:
#       dm6.fa installed in test/ directory
#       dm6_min.fa installed in test/ directory
#
# @author rahuldhodapkar
# @version 2016-08-27
# @copyright Rahul Dhodapkar

use strict;
use warnings;

use Test::Simple tests => 6;
use POSIX qw(mkfifo);

my @test_files = qw/.minix     .minix.msk \
                    .minix.prn .minix.msk.prn/;
my $out;

####################################################
## ASSERT REQUIRED FILES
####################################################

if (! -e 'test/data/dm6.fa') {
    die "unable to find required 'dm6.fa' file\n";
}

if (! -e 'test/data/dm6_min.fa') {
    die "unable to find required 'dm6_min.fa' file\n";
}

####################################################
## TEST INDEX BUILD
####################################################

$out = `./gtree ix build -r test/data/dm6_min.fa -o .minix`;
ok( $? == 0, 'build index on dm6_min.fa' );

####################################################
## TEST INDEX LOAD
####################################################

$out = `./gtree ix stat -n -ix .minix`;
ok( $? == 0, 'stat index on dm6_min.fa' );

####################################################
## TEST INDEX PRUNE
####################################################

$out = `./gtree ix prune -ix .minix -o .minix.prn`;
ok( $? == 0, 'prune index on dm6_min.fa' );

####################################################
## TEST INDEX MASK
####################################################

$out = `./gtree ix mask -r test/data/dm6.fa \\
                            -ix .minix -o .minix.msk`;
ok( $? == 0, 'mask index on dm6_min.fa' );

$out = `./gtree ix prune -ix .minix.msk -o .minix.msk.prn`;
ok( $? == 0, 'prune masked index on dm6_min.fa' );

# clean up test files
unlink( @test_files );

ok( ! -e @test_files, 'fully cleaned up' );

