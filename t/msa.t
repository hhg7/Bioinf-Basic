require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use FindBin;
use File::Temp qw(tempdir);
use File::Spec::Functions qw(catfile);
use Sys::Hostname ();
use Bioinf::Basic qw(:all);
# These run the clustalo, blastp and Python that Alien::Bioinf installed, on
# t/data/DEG20010421.fa (see t/data/make_fixtures.pl). The clustal_view_residues
# expectations are worked out by hand from the 2-sequence alignment below.
# Alien::Bioinf is only recommended, so without a working one this is skipped;
# the argument checks that need none of its programs are in t/checks.t.
unless (eval { require Alien::Bioinf; !grep { !-x $_ } Alien::Bioinf->clustalo, Alien::Bioinf->blast('blastp'), Alien::Bioinf->python }) {
	# the reason up to "(@INC contains: ...", which would fill the screen
	my ($why) = $@ ? $@ =~ /\A(.*?)(?: \(| at .+ line \d|\n|\z)/ : 'no clustalo, blastp or python';
	plan skip_all => "Alien::Bioinf can't run its tools here: $why";
}

my $dir = tempdir(CLEANUP => 1);
my $fa = "$FindBin::Bin/data/DEG20010421.fa";
my $seqs = fasta2hash($fa);
tr/-//d foreach values %{ $seqs };
sub png { open my $fh, '<:raw', $_[0] or return ''; read $fh, my $b, 8; $b }
sub slurp { open my $fh, '<:raw', $_[0] or return ''; local $/; <$fh> }
# the Creator each image's metadata carries, its whole provenance on one line:
# an SVG's <dc:title> in its <dc:creator>, a PNG's tEXt chunk
sub creator { my ($sub, $by) = @_; qr/\Q$FindBin::RealScript\E called using "$sub" in .+Basic\.pm version \Q$Bioinf::Basic::VERSION\E with Alien::Bioinf [\d.]+ by user \S.* on host \Q${\ Sys::Hostname::hostname()}\E \(\Q$^O\E\) with Perl [\d.]+ \(\Q$^X\E\), drawn by $by/ }
my $by_py = qr/.+msa_plot\.py with Python [\d.]+\S* \(.+?\), matplotlib [\d.]+/;
my $sha = qr/\(SHA-256 [0-9a-f]{64}\)/;
# what _align_provenance says of an alignment of $fa, by a run of clustalo
# that matched $cmd
sub aligned { my ($cmd) = @_; qr/; from the FASTA file \Q$fa\E $sha of 4 sequences, written with its gaps stripped to [^;]+ $sha, which Clustal Omega [\d.]+ from Alien::Bioinf [\d.]+ aligned into [^;]+ $sha, run as: [^;]*clustalo --in $cmd/ }
# an SVG's metadata, which should hold nothing but the date, format, type and Creator
sub svg_dc { my ($svg) = @_; my ($md) = $svg =~ m{<metadata>(.*?)</metadata>}s; [sort $md =~ /<dc:(\w+)/g] }

# ---- plot_msa and plot_phylo -----------------------------------------------

my %out = map { $_ => catfile($dir, $_) } qw(msa.svg tree.png aln.fa t.newick);
my $r = plot_msa(
	fasta => $fa, filename => $out{'msa.svg'}, 'msa_file' => $out{'aln.fa'}, 'tree_file' => $out{'t.newick'},
	title => 'DEG20010421', 'active_site_aa' => { Lys100 => 100 }, query => 'S.cerevisiae',
	labels => { 'S.cerevisiae' => '$\it{S. cerevisiae}$' },
);
is_deeply $r, { filename => $out{'msa.svg'}, 'msa_file' => $out{'aln.fa'}, 'tree_file' => $out{'t.newick'} },
	'plot_msa, a file in: every output named';
my $aln = fasta2hash($out{'aln.fa'});
is_deeply [sort keys %{ $aln }], [sort keys %{ $seqs }], 'the alignment has every sequence';
is scalar(keys %{ { map { length $_ => 1 } values %{ $aln } } }), 1, 'all of one length';
is_deeply { map { (my $s = $aln->{$_}) =~ tr/-//d; $_ => uc $s } keys %{ $aln } }, $seqs, 'and gapped copies of the input';
ok -s $out{'msa.svg'} && do { open my $fh, '<', $out{'msa.svg'}; local $/; <$fh> =~ /<svg/ }, 'an SVG alignment image';
my $msa_svg = slurp($out{'msa.svg'});
like $msa_svg, creator('plot_msa', $by_py . qr/; titled "DEG20010421"/ . aligned(qr/.*--guidetree-out=\Q$out{'t.newick'}\E/)
	. qr/; kept with it: msa_file \Q$out{'aln.fa'}\E, tree_file \Q$out{'t.newick'}\E<\/dc:title>/),
	'whose Creator names this script, the sub, the versions, the user, the host, Python, matplotlib, the title, the alignment and the files kept';
is_deeply svg_dc($msa_svg), [qw(creator date format title type)], 'and that is all its metadata, as from Matplotlib::Simple';
open my $nw, '<', $out{'t.newick'} or die $!;
like do { local $/; <$nw> }, qr/S\.cerevisiae:[\d.]+/, 'the guide tree is newick with the sequence names';
close $nw;

$r = plot_phylo('tree_file' => $out{'t.newick'}, 'output_file' => $out{'tree.png'}, title => 'DEG20010421',
	labels => { 'S.cerevisiae' => '$\it{S. cerevisiae}$' });
is_deeply $r, { 'output_file' => $out{'tree.png'}, 'tree_file' => $out{'t.newick'} }, 'plot_phylo draws the tree plot_msa kept';
is png($out{'tree.png'}), "\x89PNG\r\n\x1a\n", 'a PNG tree image';
my $tree_png = slurp($out{'tree.png'});
like $tree_png, creator('plot_phylo', $by_py . qr/, Biopython [\d.]+, NumPy [\d.]+; titled "DEG20010421"; from the newick file \Q$out{'t.newick'}\E $sha/
	. qr/; 4 tips: [^;]*S\.cerevisiae \(shown as \$\\it\{S\. cerevisiae\}\$\)[^;]*; the tree as drawn, in newick: \(\S*S\.cerevisiae:[\d.]+\S*\);/),
	'and its Creator, which names Biopython and NumPy, the title, the newick file with its digest, the tips with their labels and the newick';
unlike $tree_png, qr/(?:Source|Description|Title)\0/, 'and no other text chunks of provenance';

$r = plot_msa(fasta => $seqs, filename => catfile($dir, 'h.png'), order => ['S.cerevisiae', 'C.neoformans.JEC21']);
is_deeply [sort keys %{ $r }], ['filename', 'msa_file'], 'plot_msa, a hash in: the image and a temporary alignment';
is png($r->{filename}), "\x89PNG\r\n\x1a\n", 'a PNG, drawing only the sequences in "order"';

$r = plot_phylo(fasta => $seqs, 'output_file' => catfile($dir, 'p.png'));
is_deeply [sort keys %{ $r }], ['msa_file', 'output_file'], 'plot_phylo, a hash in: it aligns, and the tree is temporary';
is png($r->{'output_file'}), "\x89PNG\r\n\x1a\n", 'a PNG tree image from the sequences';

# "output_file" defaults to phylo.svg, in the working directory
{
	require Cwd;
	require Digest::SHA;
	my $was = Cwd::getcwd();
	chdir $dir or die "can't chdir to $dir: $!";
	my $keep = catfile($dir, 'kept.newick');
	$r = plot_phylo(fasta => $fa, 'tree_file' => $keep, title => 'DEG20010421');
	chdir $was or die "can't chdir back to $was: $!";
	is_deeply [sort keys %{ $r }], ['msa_file', 'output_file', 'tree_file'], 'plot_phylo with no "output_file"';
	is $r->{'output_file'}, 'phylo.svg', 'reports phylo.svg';
	my $svg = slurp(catfile($dir, 'phylo.svg'));
	like $svg, qr/<svg/, 'and writes it';
	my $nw_sha = Digest::SHA->new(256)->addfile($keep, 'b')->hexdigest;
	like $svg, creator('plot_phylo', $by_py . qr/, Biopython [\d.]+, NumPy [\d.]+; titled "DEG20010421"/ . aligned(qr/.*--guidetree-out=\Q$keep\E/)
		. qr/; drawing the guide tree it wrote, \Q$keep\E \(SHA-256 $nw_sha\); kept with it: tree_file \Q$keep\E; 4 tips: [^;]+; the tree as drawn, in newick: \(\S+\);<\/dc:title>/),
		'whose Creator names the FASTA, the ungapped copy, the alignment and the tree, each with its digest, the clustalo command, the tree kept, the tips and the newick';
	is_deeply svg_dc($svg), [qw(creator date format title type)], 'and that is all its metadata';
	unlike $svg, qr{<dc:title>DEG20010421</dc:title>}, 'with no Title beside the Creator';
}

foreach my $f ([\&plot_msa, 'filename'], [\&plot_phylo, 'output_file']) {
	my @w;
	local $SIG{__WARN__} = sub { push @w, $_[0] };
	is_deeply $f->[0]->(fasta => { a => 'MKV' }, $f->[1] => catfile($dir, 'one.png')), {}, 'one sequence: nothing made';
	like $w[0], qr/clustalo needs at least 2 sequences, and there is only 1/, 'and a warning says why';
}
my %two = (fasta => $seqs, filename => catfile($dir, 'x.png'));
foreach my $bad (
	[{ %two, fasta => [] }, qr/"fasta" must be a FASTA file name or a hash ref/],
	[{ %two, colour => 1 }, qr/doesn't know "colour"/],
	[{ fasta => $seqs }, qr/needs "filename"/],
	[{ filename => 'x.png' }, qr/needs "fasta"/],
	[{ %two, 'active_site_aa' => { X => 1 } }, qr/needs "query" to place "active_site_aa"/],
	[{ %two, 'active_site_aa' => { X => 1 }, query => 'nope' }, qr/the query "nope" isn't in the alignment/],
	[{ %two, 'active_site_aa' => { X => 9999 }, query => 'S.cerevisiae' }, qr/active site X \(9999\) is past the end of "S.cerevisiae", which has 713 residues/],
	[{ %two, order => ['nope'] }, qr/"order" names sequences that aren't in the alignment: nope/],
	[{ %two, labels => { 'C.neoformans.JEC21' => 'S.cerevisiae' } }, qr/"C.neoformans.JEC21", "S.cerevisiae" would all be drawn as "S.cerevisiae"/],
) {
	eval { plot_msa(%{ $bad->[0] }) };
	like $@, $bad->[1], "plot_msa dies: $bad->[1]";
}
foreach my $bad (
	[{ 'output_file' => 'x.png' }, qr/needs "fasta" to align, or "tree_file" to draw/],
	[{}, qr/needs "fasta" to align, or "tree_file" to draw/],
	[{ 'output_file' => 'x.png', 'tree_file' => catfile($dir, 'none.newick') }, qr/none\.newick" doesn't exist/],
	[{ 'output_file' => 'x.png', 'tree_file' => $out{'t.newick'}, threads => 2 }, qr/was given "threads" but no "fasta" to align/],
	[{ fasta => $seqs, 'output_file' => 'x.png', order => [] }, qr/doesn't know "order"/],
	[{ fasta => $seqs, filename => 'x.png' }, qr/doesn't know "filename"/],
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
is msa_quality_table(fasta => $fa, filename => $img), $img, 'a FASTA file in: blastp is run, into a temporary report';
is png($img), "\x89PNG\r\n\x1a\n", 'a PNG table';
unlink $img;
msa_quality_table(fasta => $aln, filename => $img);
is png($img), "\x89PNG\r\n\x1a\n", 'an aligned hash in: its gaps are stripped for blastp';
unlink $img;
is msa_quality_table('alignment_json' => $json, 'unaligned_fa' => $unaligned, filename => $img, metric => 'bit_score'), $img,
	'blastp is run when "alignment_json" does not exist yet, on "unaligned_fa", the old name for "fasta"';
ok -s $json, 'and its report kept';
is png($img), "\x89PNG\r\n\x1a\n", 'a PNG table';
like slurp($img), creator('msa_quality_table', $by_py . qr/; untitled; showing each pair's "bit_score"; from the FASTA file \Q$unaligned\E $sha of \d+ sequences, written with its gaps stripped to [^;]+ $sha, which blastp [\d.]+ from Alien::Bioinf [\d.]+ compared all against all into \Q${\ File::Spec->rel2abs($json)}\E $sha, run as: [^;]*blastp[^;]*; kept with it: alignment_json \Q${\ File::Spec->rel2abs($json)}\E/),
	'and its Creator, with the FASTA, the blastp command and the report kept';
unlink $img;
msa_quality_table('alignment_json' => $json, filename => $img, normalize => 1, order => [sort keys %{ $seqs }]);
is png($img), "\x89PNG\r\n\x1a\n", 'the existing report is read';
unlink $img;
msa_quality_table('alignment_json' => $json, filename => $img, normalize => 1, cblogscale => 1);
is png($img), "\x89PNG\r\n\x1a\n", 'normalize and cblogscale together: a log scale from the smallest value above 0 to 1';
foreach my $bad (
	[{ filename => $img }, qr/needs "fasta" to align, or an existing "alignment_json"/],
	[{ filename => $img, fasta => $fa, 'unaligned_fa' => $fa }, qr/both "fasta" and "unaligned_fa"/],
	[{ filename => $img, fasta => [] }, qr/"fasta" must be a FASTA file name or a hash ref/],
	[{ filename => $img, 'alignment_json' => $json, metric => 'hseq' }, qr/"hseq" isn't one of the metrics/],
	[{ filename => $img, 'alignment_json' => catfile($dir, 'new.json') }, qr/needs "fasta" to align, or an existing "alignment_json" \(.+new\.json doesn't exist yet\)/],
	[{ filename => $img, 'alignment_json' => $json, order => ['nodot'] }, qr/can't get a genus and species from "nodot"/],
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
my %view = ('msa_file' => $ab, 'output_tex_file' => $tex, 'row_width' => 3, track => 'a',
	'color_residues' => { a => { 2 => 'red', 3 => [1, 0, 0] } });
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
clustal_view_residues(%view, 'msa_file' => $out{'aln.fa'}, 'color_residues' => {}, track => undef, 'row_width' => 100);
like body($tex), qr/^\\textit\{S\.cerevisiae\} & \\texttt\{\Q${\ substr($aln->{'S.cerevisiae'}, 0, 100)}\E\} & \d+$/m,
	'plot_msa\'s alignment is shown as it is';
my $before = slurp($unaligned);
is clustal_view_residues(%view, 'msa_file' => $unaligned, 'color_residues' => {}, track => undef, 'row_width' => 1e5), $tex, # 1 block, so every row is the full width
	'sequences of different lengths are aligned rather than dying';
my @shown = body($tex) =~ /^\\textit\{[^}]+\} & \\texttt\{([^}]*)\}/mg;
ok @shown && (grep { /-/ } @shown) && !(grep { length $_ != length $shown[0] } @shown), 'into gapped rows of one width';
is slurp($unaligned), $before, 'and "msa_file" is left as it was';
foreach my $bad (
	[{ %view, track => 'z' }, qr/tracker "z" isn't in the alignment/],
	[{ %view, 'color_residues' => { z => {} } }, qr/names proteins that aren't in the alignment: z/],
	[{ %view, 'color_residues' => { a => { 9 => 'red' } } }, qr/a has no residue 9/],
	[{ %view, 'color_residues' => { a => { 1 => [1, 0] } } }, qr/must be a name or 3 numbers/],
	[{ %view, order => ['z'] }, qr/"order" names proteins that aren't in the alignment: z/],
	[{ %view, 'msa_file' => catfile($dir, 'none.fa') }, qr/"msa_file" .+none\.fa doesn't exist/],
	[{ %view, fasta => $ab }, qr/doesn't know "fasta"/],
	[{ %view, alignment => '\\centering' }, qr/doesn't know "alignment"/],
	[{ 'output_tex_file' => $tex }, qr/needs "msa_file"/],
) {
	eval { clustal_view_residues(%{ $bad->[0] }) };
	like $@, $bad->[1], "dies: $bad->[1]";
}
done_testing;
