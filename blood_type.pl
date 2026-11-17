#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Cn::Util;
use Cn::Bioinf 'blood_type';

foreach my $vcf ('father', 'mother', 'daughter') {
	my $bt = blood_type("/home/con/Documents/Genomics/$vcf.vcf.gz");
	p $bt;
}
