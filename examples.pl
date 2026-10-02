#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Devel::Confess 'color';
use Bioinf::Basic qw(:all);

my $h = fasta2hash('t/data/DEG20010421.fa');

# clustal_view_residues aligns t/data/DEG20010421.fa itself, since its sequences differ in length
plot_msa(
	fasta      => 't/data/DEG20010421.fa',
	filename   => 'DEG20010421.msa.svg',
);
clustal_view_residues(
	'msa.file'        => 't/data/DEG20010421.fa',
	'output.tex.file' => 'view_residues.tab.tex'
);
clustal_view_residues(
	'msa.file'        => 't/data/DEG20010421.fa',
	'output.tex.file' => 'ex.tex',
	'color.residues'  => {
     # residue numbers are 1-based in S.cerevisiae's own sequence, not alignment columns
     'S.cerevisiae' => { 3 => 'red', 8 => 'blue', 9 => [0, 0.6, 0] },   # xcolor name, or [r,
	},
	track       => 'S.cerevisiae',   # adds a row under it with those residue numbers
	'row.width' => 7,
	caption     => 'Catalytic residues of \textit{S. cerevisiae}',
	label       => 'active',
);
plot_phylo(
	fasta => 't/data/DEG20010421.fa',
	filename => 'phylo.svg'
);   # phylo.svg

