require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use File::Temp qw(tempdir);
use File::Spec::Functions qw(catfile);
use Bioinf::Basic qw(:all);
# Checks that need none of Alien::Bioinf's programs, so that they run wherever
# the module builds: each one dies, or is decided, before clustalo, blastp or
# Python would be started. Every case is a bug that was found by review on
# 2026-10-01 and fixed then.

my $dir = tempdir(CLEANUP => 1);
sub spew {
	my ($name, $content) = @_;
	my $f = catfile($dir, $name);
	open my $fh, '>:raw', $f or die $!;
	print {$fh} $content;
	close $fh;
	$f;
}
sub dies_like {
	my ($code, $re, $name) = @_;
	eval { $code->(); 1 } and return fail($name);
	like $@, $re, $name;
}
sub slurp { open my $fh, '<', $_[0] or die "$_[0]: $!"; local $/; <$fh> }
my @inc = map { "-I$_" } grep { !ref } @INC;

# ---- loading ----------------------------------------------------------------

{
	package Plain;
	Bioinf::Basic->import;
	::ok !defined &Plain::fasta2hash && !defined &Plain::plot_msa, 'nothing is exported by default';
}
is system($^X, @inc, '-e', 'BEGIN { $0 = q{/no/such/script.pl} } require Bioinf::Basic'), 0,
	'the module loads when $0 is not a file, which FindBin would croak on';
# Alien::Bioinf hidden from require, as on a machine that could not install it
my $hide = 'BEGIN { unshift @INC, sub { die qq{hidden\n} if $_[1] eq q{Alien/Bioinf.pm}; return } }';
my $err = catfile($dir, 'err.txt');
# $err is passed in @ARGV, not pasted into the code, which a Windows path's
# backslashes or a '}' in it would break
system($^X, @inc, '-e', "$hide use Bioinf::Basic q{plot_msa}; open STDERR, q{>}, \$ARGV[0]; plot_msa(fasta => { a => q{MK}, b => q{MV} }, filename => q{x.png})", $err);
like slurp($err), qr/^_alien: plot_msa needs Alien::Bioinf, for Clustal Omega, BLAST\+ and the plotting Python, and it can't be loaded: hidden/,
	'without Alien::Bioinf a plotting function says what it needs';

# ---- _run -------------------------------------------------------------------

# Up to 2026-10-07 _run read only the exit code, the high byte of $?, so a
# program killed by a signal was reported as "exit 0".
dies_like sub { Bioinf::Basic::_run($^X, '-e', 'exit 3') }, qr/^_run: .+ failed: exit 3 at /, 'a non-zero exit is reported';
SKIP: {
	skip 'MSWin32 has no POSIX signals', 1 if $^O eq 'MSWin32';
	# set to DEFAULT in the child, which may have inherited TERM ignored
	dies_like sub { Bioinf::Basic::_run($^X, '-e', '$SIG{TERM} = q{DEFAULT}; kill q{TERM}, $$; sleep 5') },
		qr/^_run: .+ failed: killed by signal 15/, 'a death by signal is reported as one, not as "exit 0"';
}

# ---- undef file names -------------------------------------------------------

# Up to 2026-10-07 the message interpolated the undef, which warned
# "uninitialized value" -- fatally, under a caller's FATAL warnings, which then
# hid the message itself.
foreach my $case (
	[fasta2hash             => sub { fasta2hash(undef) },             qr/^_open_fasta: undef isn't a readable file/],
	[get_best_alignment_hit => sub { get_best_alignment_hit(undef) }, qr/^get_best_alignment_hit: undef isn't a readable file/],
) {
	my ($name, $code, $re) = @{ $case };
	my @w;
	local $SIG{__WARN__} = sub { push @w, $_[0] };
	dies_like $code, $re, "$name(undef) dies naming the undef";
	is "@w", '', 'and with no warning';
}

# ---- hash2fasta_file --------------------------------------------------------

# Up to 2026-10-07 the width went straight to the XS as a UV: -1 wrapped round
# to one line per sequence, and 'x' warned with no function name.
foreach my $w (-1, 1.5, 'x') {
	dies_like sub { hash2fasta_file({ a => 'MK' }, catfile($dir, 'w.fa'), undef, $w) },
		qr/^hash2fasta_file: 4th argument must be a whole number, 0 or more, or undef, not "\Q$w\E"/, "width $w dies";
}

# ---- get_best_alignment_hit -------------------------------------------------

my $report = '{"report":{"results":{"search":{"query_title":"q","hits":[]}}}}';
dies_like sub { get_best_alignment_hit(spew('twice.json', qq({"BlastOutput2":[$report,$report]}))) },
	qr/^get_best_alignment_hit: two queries in .+twice\.json are both titled "q"/, 'two queries of one title die rather than one being lost';

# ---- plot_msa ---------------------------------------------------------------

my %two = (fasta => { a => 'MK', b => 'MV' }, filename => catfile($dir, 'x.png'), query => 'a');
foreach my $n (0, -1, 1.5, 'two') {
	dies_like sub { plot_msa(%two, 'active_site_aa' => { X => $n }) },
		qr/^plot_msa: active site X must be a residue number, 1 or more, not "\Q$n\E"/, "active site $n dies, rather than counting from the end";
}
dies_like sub { plot_msa(%two, 'active_site_aa' => [1]) }, qr/^plot_msa: "active_site_aa" must be a hash ref/, 'active_site_aa must be a hash ref';

# ---- msa_quality_table ------------------------------------------------------

# Two sequences, every pair with an e-value of 0, as BLAST reports strong hits
my $hit = sub { { description => [{ title => $_[0] }], hsps => [{ evalue => 0 }] } };
my $blast = { BlastOutput2 => [map {
	my $q = $_;
	{ report => { results => { bl2seq => [{ query_title => $q, hits => [map { $hit->($_) } 'A.b', 'C.d'] }] } } }
} 'A.b', 'C.d'] };
my %q = ('alignment_json' => $blast, filename => catfile($dir, 'q.png'), metric => 'evalue');
dies_like sub { msa_quality_table(%q, normalize => 1) },
	qr/^msa_quality_table: can't normalize "evalue", since its largest value is 0, not above 0/, 'normalizing all-zero e-values dies, rather than dividing by 0';
dies_like sub { msa_quality_table(%q, order => ['E.f']) },
	qr/^msa_quality_table: no pair of the sequences in "order" has a "evalue" value/, 'an order with nothing to draw dies';

# ---- clustal_view_residues --------------------------------------------------

my $ab = hash2fasta_file({ 'a_b^c~d' => 'MK-LV', 'e{f}' => 'M-QLV' }, catfile($dir, 'ab.aln.fa'));
my $tex = catfile($dir, 'r.tex');
my %view = ('msa_file' => $ab, 'output_tex_file' => $tex);
foreach my $k ('row_width', 'split') {
	foreach my $v (0, -3, 'x') {
		dies_like sub { clustal_view_residues(%view, $k => $v) }, qr/^clustal_view_residues: "\Q$k\E" must be a whole number, 1 or more, not "\Q$v\E"/,
			"$k $v dies" . ($k eq 'row_width' && $v eq '0' ? ', rather than looping forever' : '');
	}
}
foreach my $n (0, -1) {
	dies_like sub { clustal_view_residues(%view, 'color_residues' => { 'e{f}' => { $n => 'red' } }) },
		qr/^clustal_view_residues: a residue of e\{f\} in "color_residues" must be a residue number, 1 or more, not "$n"/,
		"colouring residue $n dies, rather than colouring the last";
}
{
	open my $saved, '>&', \*STDOUT or die $!;
	open STDOUT, '>', catfile($dir, 'out.txt') or die $!;
	clustal_view_residues(%view);
	open STDOUT, '>&', $saved or die $!;
}
is slurp(catfile($dir, 'out.txt')), "wrote $tex\n", 'the "wrote" line has no colour when STDOUT is not a terminal';
my $body = slurp($tex);
like $body, qr/^\\textit\{a\\_b\\textasciicircum\{\}c\\textasciitilde\{\}d\} & /m, 'a name\'s _, ^ and ~ print as themselves';
like $body, qr/^\\textit\{e\\\{f\\\}\} & /m, 'and so do its braces';

done_testing;
