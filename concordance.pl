#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
#use Stats::LikeR;
use Cn::Util;
use Cn::Bioinf;

# The test VCF and then the truth VCF may be given on the command line; the
# defaults are the links in ~/Documents/Genomics, as parent.pl's are, so that
# no sample's ID is written here.
my $test  = shift // "$ENV{HOME}/Documents/Genomics/genome.vcf.gz";
my $truth = shift // "$ENV{HOME}/Documents/Genomics/father.vcf.gz";
my $c = concordance(test => $test, truth => $truth);
p $c;
