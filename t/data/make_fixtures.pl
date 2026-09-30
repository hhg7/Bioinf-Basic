#!/usr/bin/env perl
# Regenerates query.fa, subjects.fa and hits.json in this directory; the tests
# only read what it wrote. Re-run from the distribution root, with Alien::Bioinf
# installed, as
#
#   perl t/data/make_fixtures.pl
#
# DEG20010421.fa is ~/identify.target/cryptococcus/fa/clustal.DEG20010421.input.fa
# (copied 2026-09-29): S.cerevisiae DEG20010421 and its three C.neoformans hits,
# with a few gaps left in from an earlier alignment, which are stripped here.
# hits.json is blastp 2.17.0+ (-outfmt 15) of S.cerevisiae against a database of
# the three C.neoformans sequences plus two made-up subjects that exercise
# get_best_alignment_hit's awkward cases:
#   - JEC21's sequence again under H99's defline: two hits with one title,
#     which the old hash-keyed version silently merged;
#   - S.cerevisiae residues 1-150 and 301-450 joined by 60 random residues
#     (srand 20260929), which blastp reports as one hit with two hsps.
use strict;
use warnings;
use FindBin;
use Alien::Bioinf;

chdir $FindBin::Bin or die "$FindBin::Bin: $!";
my (%seq, $name);
open my $in, '<', 'DEG20010421.fa' or die $!;
while (<$in>) {
	chomp;
	if (/^>(.+)/) { $name = $1; next }
	(my $s = $_) =~ tr/-//d;
	$seq{$name} .= $s;
}
close $in;
srand 20260929;
my @aa = split //, 'ACDEFGHIKLMNPQRSTVWY';
my $filler = join '', map { $aa[int rand @aa] } 1 .. 60;
my $q = $seq{'S.cerevisiae'};
open my $s, '>', 'subjects.fa' or die $!;
print {$s} ">$_\n$seq{$_}\n" foreach qw(C.neoformans.B.3501A C.neoformans.grubii.H99 C.neoformans.JEC21);
print {$s} ">C.neoformans.grubii.H99\n$seq{'C.neoformans.JEC21'}\n";
print {$s} '>S.cerevisiae.chimera' . "\n" . substr($q, 0, 150) . $filler . substr($q, 300, 150) . "\n";
close $s;
open my $qf, '>', 'query.fa' or die $!;
print {$qf} ">S.cerevisiae\n$q\n";
close $qf;
my $db = "$FindBin::Bin/subjects.db";
system(Alien::Bioinf->blast('makeblastdb'), '-in', 'subjects.fa', '-dbtype', 'prot', '-out', $db) == 0 or die 'makeblastdb failed';
system(Alien::Bioinf->blast('blastp'), '-query', 'query.fa', '-db', $db, '-outfmt', 15, '-out', 'hits.json') == 0 or die 'blastp failed';
unlink glob "$db.*";
