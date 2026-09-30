require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use FindBin;
use File::Temp qw(tempdir);
use File::Spec::Functions qw(catfile);
use Bioinf::Basic qw(:all);
# These run the clustalo, blastp and Python that Alien::Bioinf installed, on
# t/data/DEG20010421.fa (see t/data/make_fixtures.pl). The clustal_view_residues
# expectations are worked out by hand from the 2-sequence alignment below.

my $dir = tempdir(CLEANUP => 1);
my $fa = "$FindBin::Bin/data/DEG20010421.fa";
my $seqs = fasta2hash($fa);
tr/-//d foreach values %{ $seqs };
sub png { open my $fh, '<:raw', $_[0] or return ''; read $fh, my $b, 8; $b }
sub slurp { open my $fh, '<:raw', $_[0] or return ''; local $/; <$fh> }
# the Creator each image's metadata carries: an SVG's <dc:title>, a PNG's tEXt chunk
sub creator { my ($sub, $by) = @_; qr/\Q$FindBin::RealScript\E called using "$sub" in \S+Basic\.pm version \Q$Bioinf::Basic::VERSION\E, drawn by $by/ }
my $by_py = qr/\S+msa_plot\.py with matplotlib [\d.]+/;

# ---- plot_msa and plot_phylo -----------------------------------------------

my %out = map { $_ => catfile($dir, $_) } qw(msa.svg tree.png aln.fa t.newick);
my $r = plot_msa(
	fasta => $fa, filename => $out{'msa.svg'}, 'msa.file' => $out{'aln.fa'}, 'tree.file' => $out{'t.newick'},
	title => 'DEG20010421', 'active.site.aa' => { Lys100 => 100 }, query => 'S.cerevisiae',
	labels => { 'S.cerevisiae' => '$\it{S. cerevisiae}$' },
);
is_deeply $r, { filename => $out{'msa.svg'}, 'msa.file' => $out{'aln.fa'}, 'tree.file' => $out{'t.newick'} },
	'plot_msa, a file in: every output named';
my $aln = fasta2hash($out{'aln.fa'});
is_deeply [sort keys %{ $aln }], [sort keys %{ $seqs }], 'the alignment has every sequence';
is scalar(keys %{ { map { length $_ => 1 } values %{ $aln } } }), 1, 'all of one length';
is_deeply { map { (my $s = $aln->{$_}) =~ tr/-//d; $_ => uc $s } keys %{ $aln } }, $seqs, 'and gapped copies of the input';
ok -s $out{'msa.svg'} && do { open my $fh, '<', $out{'msa.svg'}; local $/; <$fh> =~ /<svg/ }, 'an SVG alignment image';
like slurp($out{'msa.svg'}), creator('plot_msa', $by_py), 'whose Creator names this script, the sub, the versions and matplotlib';
open my $nw, '<', $out{'t.newick'} or die $!;
like do { local $/; <$nw> }, qr/S\.cerevisiae:[\d.]+/, 'the guide tree is newick with the sequence names';
close $nw;

$r = plot_phylo('tree.file' => $out{'t.newick'}, filename => $out{'tree.png'}, title => 'DEG20010421',
	labels => { 'S.cerevisiae' => '$\it{S. cerevisiae}$' });
is_deeply $r, { filename => $out{'tree.png'}, 'tree.file' => $out{'t.newick'} }, 'plot_phylo draws the tree plot_msa kept';
is png($out{'tree.png'}), "\x89PNG\r\n\x1a\n", 'a PNG tree image';
like slurp($out{'tree.png'}), creator('plot_phylo', $by_py), 'and its Creator';

$r = plot_msa(fasta => $seqs, filename => catfile($dir, 'h.png'), order => ['S.cerevisiae', 'C.neoformans.JEC21']);
is_deeply [sort keys %{ $r }], ['filename', 'msa.file'], 'plot_msa, a hash in: the image and a temporary alignment';
is png($r->{filename}), "\x89PNG\r\n\x1a\n", 'a PNG, drawing only the sequences in "order"';

$r = plot_phylo(fasta => $seqs, filename => catfile($dir, 'p.png'));
is_deeply [sort keys %{ $r }], ['filename', 'msa.file'], 'plot_phylo, a hash in: it aligns, and the tree is temporary';
is png($r->{filename}), "\x89PNG\r\n\x1a\n", 'a PNG tree image from the sequences';

foreach my $f (\&plot_msa, \&plot_phylo) {
	my @w;
	local $SIG{__WARN__} = sub { push @w, $_[0] };
	is_deeply $f->(fasta => { a => 'MKV' }, filename => catfile($dir, 'one.png')), {}, 'one sequence: nothing made';
	like $w[0], qr/clustalo needs at least 2 sequences, and there is only 1/, 'and a warning says why';
}
my %two = (fasta => $seqs, filename => catfile($dir, 'x.png'));
foreach my $bad (
	[{ %two, fasta => [] }, qr/"fasta" must be a FASTA file name or a hash ref/],
	[{ %two, colour => 1 }, qr/doesn't know "colour"/],
	[{ fasta => $seqs }, qr/needs "filename"/],
	[{ filename => 'x.png' }, qr/needs "fasta"/],
	[{ %two, 'active.site.aa' => { X => 1 } }, qr/needs "query" to place "active.site.aa"/],
	[{ %two, 'active.site.aa' => { X => 1 }, query => 'nope' }, qr/the query "nope" isn't in the alignment/],
	[{ %two, 'active.site.aa' => { X => 9999 }, query => 'S.cerevisiae' }, qr/active site X \(9999\) is past the end of "S.cerevisiae", which has 713 residues/],
	[{ %two, order => ['nope'] }, qr/"order" names sequences that aren't in the alignment: nope/],
) {
	eval { plot_msa(%{ $bad->[0] }) };
	like $@, $bad->[1], "plot_msa dies: $bad->[1]";
}
foreach my $bad (
	[{ filename => 'x.png' }, qr/needs "fasta" to align, or "tree.file" to draw/],
	[{ filename => 'x.png', 'tree.file' => catfile($dir, 'none.newick') }, qr/none\.newick" doesn't exist/],
	[{ filename => 'x.png', 'tree.file' => $out{'t.newick'}, threads => 2 }, qr/was given "threads" but no "fasta" to align/],
	[{ %two, order => [] }, qr/doesn't know "order"/],
	[{ fasta => $seqs }, qr/needs "filename"/],
) {
	eval { plot_phylo(%{ $bad->[0] }) };
	like $@, $bad->[1], "plot_phylo dies: $bad->[1]";
}
foreach my $f (\&plot_msa, \&plot_phylo, \&msa_quality_table, \&clustal_view_residues) {
	eval { $f->('x') };
	like $@, qr/takes name => value pairs, and was given an odd number of arguments/, 'dies on an odd number of arguments';
	eval { $f->({ %two }) };
	like $@, qr/takes name => value pairs, not a hash ref/, 'dies on a hash ref';
}

# ---- msa_quality_table ------------------------------------------------------

my $unaligned = hash2fasta_file($seqs, catfile($dir, 'u.fa'));
my $json = catfile($dir, 'all.json');
my $img = catfile($dir, 'q.png');
is msa_quality_table('alignment.json' => $json, 'unaligned.fa' => $unaligned, filename => $img, metric => 'bit_score'), $img,
	'blastp is run when "alignment.json" does not exist yet';
ok -s $json, 'and its report kept';
is png($img), "\x89PNG\r\n\x1a\n", 'a PNG table';
like slurp($img), creator('msa_quality_table', $by_py), 'and its Creator';
unlink $img;
msa_quality_table('alignment.json' => $json, filename => $img, normalize => 1, order => [sort keys %{ $seqs }]);
is png($img), "\x89PNG\r\n\x1a\n", 'the existing report is read';
unlink $img;
msa_quality_table('alignment.json' => $json, filename => $img, normalize => 1, cblogscale => 1);
is png($img), "\x89PNG\r\n\x1a\n", 'normalize and cblogscale together: a log scale from the smallest value above 0 to 1';
foreach my $bad (
	[{ filename => $img }, qr/needs "alignment.json"/],
	[{ filename => $img, 'alignment.json' => $json, metric => 'hseq' }, qr/"hseq" isn't one of the metrics/],
	[{ filename => $img, 'alignment.json' => catfile($dir, 'new.json') }, qr/"unaligned.fa" must be given/],
	[{ filename => $img, 'alignment.json' => $json, order => ['nodot'] }, qr/can't get a genus and species from "nodot"/],
) {
	eval { msa_quality_table(%{ $bad->[0] }) };
	like $@, $bad->[1], "dies: $bad->[1]";
}

# ---- clustal_view_residues --------------------------------------------------

# a: M K - L V   residues 1 2 . 3 4
# b: M - Q L V
# Residue 2 of a is column 1 and residue 3 is column 3; 3 columns per block.
my $tex = catfile($dir, 'r.tex');
my $ab = hash2fasta_file({ a => 'MK-LV', b => 'M-QLV' }, catfile($dir, 'ab.aln.fa'));
my %view = ('msa.file' => $ab, 'output.tex.file' => $tex, 'row.width' => 3, track => 'a',
	'color.residues' => { a => { 2 => 'red', 3 => [1, 0, 0] } });
sub body { open my $fh, '<', $_[0] or die $!; <$fh>; local $/; <$fh> } # all but the provenance line
my $rows = <<'EOT';
\textit{a} & \texttt{M\textcolor{red}{K}-} & 2
\\
\textit{a} track & \texttt{-2-} & 2
\\
\textit{b} & \texttt{M\textcolor{red}{-}Q} & 2
\\
\hline
\textit{a} & \texttt{{\color[rgb]{1,0,0}L}V} & 4
\\
\textit{a} track & \texttt{3-} & 4
\\
\textit{b} & \texttt{{\color[rgb]{1,0,0}L}V} & 4
\\
\hline
EOT
my $head = "\\begin{table}[htp]\\footnotesize\n\\begin{tabular}{|c|l|c|} \\hline\n\\textbf{Track} & \\textbf{Sequence} & \\textbf{Length}\\\\ \\hline\n";
is clustal_view_residues(%view), $tex, 'returns the file name';
is body($tex), "$head$rows\\end{tabular}\n\\caption{}\n\\end{table} \\FloatBarrier\n", 'one table, coloured and tracked';
my @rows = split /(?<=\\hline\n)/, $rows;
clustal_view_residues(%view, split => 1, caption => 'C', label => 'L');
is body($tex), "$head$rows[0]\\end{tabular}\n\\caption{C}\n\\label{tab:L0}\n\\end{table} \\FloatBarrier\n"
	. "$head$rows[1]\\end{tabular}\n\\caption{C (continued)}\n\\label{tab:L1}\n\\end{table} \\FloatBarrier\n", 'split into one block per table';
clustal_view_residues(%view, 'msa.file' => $out{'aln.fa'}, 'color.residues' => {}, track => undef, 'row.width' => 100);
like body($tex), qr/^\\textit\{S\.cerevisiae\} & \\texttt\{\Q${\ substr($aln->{'S.cerevisiae'}, 0, 100)}\E\} & \d+$/m,
	'plot_msa\'s alignment is shown as it is';
foreach my $bad (
	[{ %view, track => 'z' }, qr/tracker "z" isn't in the alignment/],
	[{ %view, 'color.residues' => { z => {} } }, qr/names proteins that aren't in the alignment: z/],
	[{ %view, 'color.residues' => { a => { 9 => 'red' } } }, qr/a has no residue 9/],
	[{ %view, 'color.residues' => { a => { 1 => [1, 0] } } }, qr/must be a name or 3 numbers/],
	[{ %view, order => ['z'] }, qr/"order" names proteins that aren't in the alignment: z/],
	[{ %view, 'msa.file' => $unaligned }, qr/u\.fa isn't aligned: its sequences are of different lengths/],
	[{ %view, 'msa.file' => catfile($dir, 'none.fa') }, qr/"msa.file" \S+none\.fa doesn't exist/],
	[{ %view, fasta => $ab }, qr/doesn't know "fasta"/],
	[{ %view, alignment => '\\centering' }, qr/doesn't know "alignment"/],
	[{ 'output.tex.file' => $tex }, qr/needs "msa.file"/],
) {
	eval { clustal_view_residues(%{ $bad->[0] }) };
	like $@, $bad->[1], "dies: $bad->[1]";
}
done_testing;
