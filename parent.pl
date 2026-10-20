#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Cn::Bioinf qw(:all);
use Cn::Util;

my $father   = "$ENV{HOME}/Documents/Genomics/father.vcf.gz";
my $daughter = "$ENV{HOME}/Documents/Genomics/daughter.vcf.gz";
my $mother   = "$ENV{HOME}/Documents/Genomics/mother.vcf.gz";
my ($is, $stats) = parent(parent => $father, child => $daughter);
p $is;
p $stats;
($is, $stats) = parent(parent => $mother, child => $father);
p $is; p $stats;
