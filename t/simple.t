#!/usr/bin/perl

use strict;
use Test::More tests=>20;
use File::Basename;
use File::Spec;

chdir(dirname(__FILE__));

# set of very simple tests to check options and fastq i/o
my $dsrc=File::Spec->catfile("..", "bin", "dsrc");

for my $m (qw(m0 m1)) {
for my $f (qw(test1 test2)) {

    # remove prev out, if any
    unlink "$f-$m.fq.out.dz", "$f-$m.fq.out";

    # ****5 tests per loop*****


    # standard test
    sys("$dsrc c -m$m $f.fq $f-$m.fq.out.dz", "compress $m $f");
    sys("$dsrc d $f-$m.fq.out.dz $f-$m.fq.out", "decompress $m $f");
    same("$f.fq", "$f-$m.fq.out"); 

    #  stream decompress
    sys("$dsrc d -s $f-$m.fq.out.dz > $f-$m-s.fq.out", "decompress stream $m $f");
    same("$f.fq", "$f-$m-s.fq.out"); 
}
}

sub sys {
    if (!ok(!system($_[0]), $_[1])) {
        diag($_[0]);
    }
}

sub slurp {
    my ($f) = @_;
    local $/= undef;
    open S, $f;
    my $r=<S>;
    close S;
    return $r;
}

sub same {
    my ($f1, $f2) = @_;
    my $r1=slurp($f1);
    my $r2=slurp($f2);
    ok($r1 && $r2 && ($r1 eq $r2), "$f1 = $f2");
}

