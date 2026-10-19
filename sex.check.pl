#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Devel::Confess 'color';
use Cn::Bioinf qw(:all);
use Cn::Util;

foreach my $vcf (dir(dir => '/home/con/Documents/Genomics', regex => '^\d+\.vcf\.gz$')) {
	my $sc = sex_check($vcf);
	p $sc;
}
