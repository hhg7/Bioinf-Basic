#!/usr/bin/env perl
# The simplest use of pgs_score: a polygenic score for every individual in a
# VCF, as a hash of individual => score. PLINK 2 does the scoring, from the
# files, so neither is read into perl.
#
#   perl pgs.score.KISS.pl [weights] [vcf]
#
# The score of a VCF that lists only the sites where an individual carries an
# ALT allele, as most single-genome VCFs do, leaves out every variant where the
# individual matches the reference, so it can't be compared with anyone else's;
# pgs_score warns of it, and pgs.score.pl fills those sites in and compares the
# score with 1000 Genomes.

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use FindBin;
use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch"; # this directory's build of Cn::Bioinf, once made
use Cn::Bioinf 'pgs_score';
use Cn::Util;

my $weights_file = shift // 'PGS000138.txt.gz';
my $vcf_file     = shift // '/home/con/Documents/Genomics/genome.vcf.gz';

my $s = pgs_score(vcf => $vcf_file, weights => $weights_file);
my %prs = map { $_ => $s->{$_}{score} } keys %{ $s };
p %prs;
