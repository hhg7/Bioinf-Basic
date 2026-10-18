#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
#use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
#use Devel::Confess 'color';
use Stats::LikeR;
use Bioinf::Basic 'pgs_score';
use Util::Cn;

my $s = pgs_score(
	vcf     => read_table('/home/con/Documents/Genomics/genome.vcf.gz'),
	weights => read_table('PGS000138.txt.gz', sep => "\t"),
);
p $s;
