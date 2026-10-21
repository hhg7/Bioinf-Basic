#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
#use Stats::LikeR;
use Cn::Util;
use Cn::Bioinf;

my $c = concordance(test => '/home/con/Documents/Genomics/56001801066929_WGZ.snp.vcf.gz', truth => '/home/con/Documents/Genomics/3025181.vcf.gz');
p $c;
