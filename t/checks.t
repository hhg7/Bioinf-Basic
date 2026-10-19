require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use File::Temp qw(tempdir);
use File::Spec::Functions qw(catdir catfile rel2abs);
use Cwd ();
use Cn::Bioinf qw(:all);
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
	Cn::Bioinf->import;
	::ok !defined &Plain::fasta2hash && !defined &Plain::plot_msa, 'nothing is exported by default';
}
is system($^X, @inc, '-e', 'BEGIN { $0 = q{/no/such/script.pl} } require Cn::Bioinf'), 0,
	'the module loads when $0 is not a file, which FindBin would croak on';
# Alien::Bioinf hidden from require, as on a machine that could not install it
my $hide = 'BEGIN { unshift @INC, sub { die qq{hidden\n} if $_[1] eq q{Alien/Bioinf.pm}; return } }';
my $err = catfile($dir, 'err.txt');
# $err is passed in @ARGV, not pasted into the code, which a Windows path's
# backslashes or a '}' in it would break
system($^X, @inc, '-e', "$hide use Cn::Bioinf q{plot_msa}; open STDERR, q{>}, \$ARGV[0]; plot_msa(fasta => { a => q{MK}, b => q{MV} }, filename => q{x.png})", $err);
like slurp($err), qr/^_alien: plot_msa needs Alien::Bioinf, for Clustal Omega, BLAST\+ and the plotting Python, and it can't be loaded: hidden/,
	'without Alien::Bioinf a plotting function says what it needs';

# Up to 2026-10-17 _creator joined the working directory to the script's file
# name, so a script run from elsewhere was recorded at a path that isn't it.
# Alien::Bioinf is stood in for by a stub, since _creator only reads its
# version; it is put in %INC rather than found through @INC, which an @INC hook
# in PERL5OPT could hide it from.
{
	my $scripts = catdir($dir, 'scripts');
	mkdir $scripts or die "$scripts: $!";
	my $script = catfile($scripts, 'creator.pl');
	open my $fh, '>', $script or die $!;
	print {$fh} 'BEGIN { $Alien::Bioinf::VERSION = q{9.99}; $INC{q{Alien/Bioinf.pm}} = __FILE__ } use Cn::Bioinf;'
		. ' open STDOUT, q{>}, $ARGV[0] or die $!; print Cn::Bioinf::_creator(q{x});';
	close $fh;
	my $out = catfile($dir, 'creator.txt');
	my $was = Cwd::getcwd();
	chdir $dir or die "$dir: $!";
	# absolute, since the relative entries prove -b gives would not survive the chdir
	system($^X, (map { '-I' . rel2abs(substr $_, 2) } @inc), $script, $out);
	chdir $was or die "$was: $!";
	my $real = catfile(Cwd::realpath($scripts), 'creator.pl');
	like slurp($out), qr/^\Q$real\E called using "x" in .+ with Alien::Bioinf 9\.99 by user /,
		"the Creator names the script's own path, not the working directory's";
}

# ---- _run -------------------------------------------------------------------

# Up to 2026-10-07 _run read only the exit code, the high byte of $?, so a
# program killed by a signal was reported as "exit 0".
dies_like sub { Cn::Bioinf::_run($^X, '-e', 'exit 3') }, qr/^_run: .+ failed: exit 3 at /, 'a non-zero exit is reported';
SKIP: {
	skip 'MSWin32 has no POSIX signals', 1 if $^O eq 'MSWin32';
	# set to DEFAULT in the child, which may have inherited TERM ignored
	dies_like sub { Cn::Bioinf::_run($^X, '-e', '$SIG{TERM} = q{DEFAULT}; kill q{TERM}, $$; sleep 5') },
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

# ---- pgs_score --------------------------------------------------------------

# Added with pgs_score on 2026-10-17; t/pgs.t has the scores themselves. Every
# case that dies does so before PLINK 2 starts; the ones that score need it.
{
	my $v = { 'chr1:5:A:G' => { CHROM => 'chr1', POS => 5, ID => 'rs1', REF => 'A', ALT => 'G', 'S.GT' => '0/1' } };
	my $wrow = { rsID => 'rs1', chr_name => '1', chr_position => 5, effect_allele => 'G', other_allele => 'A', effect_weight => '0.5' };
	my $w = [$wrow];
	my $plink2 = eval { require Alien::Bioinf; -x Alien::Bioinf->plink2 };
	SKIP: {
		skip 'PLINK 2, which Alien::Bioinf installs, is not installed', 1 unless $plink2;
		is pgs_score(vcf => $v, weights => $w)->{S}{score}, 0.5, 'pgs_score: the arguments below differ from these only where they say';
	}
	my $head = "rsID\tchr_name\tchr_position\teffect_allele\teffect_weight\n";
	my @cases = (
		[[vcf => $v], qr/^pgs_score: needs "weights"/, 'no weights'],
		[[{ vcf => $v, weights => $w }], qr/^pgs_score: takes name => value pairs, not a hash ref/, 'a hash ref'],
		[[vcf => $v, weights => $w, matches => 'id'], qr/^pgs_score: doesn't know "matches"; it takes build, match, memory, threads, vcf, weights/, 'a misspelled option'],
		[[vcf => $v, weights => $w, match => 'rsid'], qr/^pgs_score: "match" must be "position" or "id", not "rsid"/, 'match "rsid"'],
		[[vcf => $v, weights => $w, build => 'GRCh37'], qr/^pgs_score: "build" must be hg19, b37, hg38 or b38, not "GRCh37"/, 'build "GRCh37"'],
		[[vcf => $v, weights => $w, memory => 639], qr/^pgs_score: "memory" must be 640 or more \(MiB\), the least PLINK 2 takes, not 639/, 'memory 639'],
		[[vcf => $v, weights => $w, memory => '2G'], qr/^pgs_score: "memory" must be a whole number, 1 or more, not "2G"/, 'memory "2G"'],
		[[vcf => $v, weights => $w, threads => 0], qr/^pgs_score: "threads" must be a whole number, 1 or more, not "0"/, 'threads 0'],
		[[vcf => \'calls', weights => $w], qr/^pgs_score: "vcf" must be a VCF, BCF or PLINK 2 \.pgen file, or a data frame, not a SCALAR ref/, 'a scalar ref'],
		[[vcf => catfile($dir, 'none.vcf.gz'), weights => $w], qr/^pgs_score: "vcf" file ".+none\.vcf\.gz" doesn't exist or isn't readable/, 'a VCF that is not there'],
		[[vcf => spew('calls.txt', ''), weights => $w], qr/^pgs_score: "vcf" file .+calls\.txt must end in \.vcf, \.vcf\.gz, \.vcf\.bgz, \.bcf or \.pgen/, 'a VCF named .txt'],
		[[vcf => spew('lone.pgen', ''), weights => $w], qr/^pgs_score: "vcf" .+lone\.pgen has no .+lone\.pvar \(or \.pvar\.zst\) beside it/, 'a .pgen with no .pvar'],
		[[vcf => $v, weights => catfile($dir, 'none.txt.gz')], qr/^pgs_score: "weights" file ".+none\.txt\.gz" doesn't exist or isn't readable/, 'a scoring file that is not there'],
		[[vcf => $v, weights => spew('comments.txt', "#pgs_id=X\n#genome_build=NR\n")], qr/^pgs_score: "weights" file .+comments\.txt has no header line/, 'a scoring file of comments alone'],
		[[vcf => $v, weights => spew('header.txt', $head)], qr/^pgs_score: "weights" has no rows/, 'a scoring file of a header alone'],
		[[vcf => $v, weights => spew('bad.txt', $head . "rs1\t1\t5\tG\t0.5\nrs2\t1\t6\tG\tx\n")], qr/^pgs_score: row 2 of "weights" has effect_weight "x", which isn't a number/, 'a bad weight in a scoring file'],
		[[vcf => $v, weights => [{ %{ $wrow }, chr_name => 'X' }]], qr/^pgs_score: "weights" has variants on X, which PLINK 2 reads only when "build" says where X's pseudo-autosomal regions are: hg19 \(b37\) or hg38 \(b38\)/, 'a variant on X without build'],
		[[vcf => $v, weights => []], qr/^pgs_score: "weights" has no rows/, 'an empty array'],
		[[vcf => {}, weights => $w], qr/^pgs_score: "vcf" has no rows/, 'an empty hash'],
		[[vcf => $v, weights => [['effect_allele']]], qr/^pgs_score: "weights" has a header and no rows/, 'an array of arrays of only a header'],
		[[vcf => $v, weights => [$wrow, undef]], qr/^pgs_score: row 2 of "weights" must be a hash ref, not undef/, 'an undef row'],
		[[vcf => $v, weights => [[qw(a b)], 'x']], qr/^pgs_score: row 1 of "weights" must be an array ref, as its header is, not "x"/, 'a string row of an aoa'],
		# the shape is taken from whichever column perl's hash order gives first
		[[vcf => $v, weights => { a => [1], b => 2 }],
			qr/^pgs_score: (?:column "b" of "weights" must be an array ref, as "a" is, not "2"|"weights" must be .+, not a hash whose "b" is "2")/,
			'a string column of a hoa'],
		[[vcf => $v, weights => [1]], qr/^pgs_score: "weights" must be .+, not an array whose first row is "1"/, 'an array of strings'],
		[[vcf => { a => 1 }, weights => $w], qr/^pgs_score: "vcf" must be .+, not a hash whose "a" is "1"/, 'a hash of strings'],
		[[vcf => $v, weights => [{ map { $_ => $wrow->{$_} } qw(rsID effect_allele effect_weight) }]],
			qr/^pgs_score: "weights" has no "chr_name", "chr_position" columns, which a PGS Catalog scoring file has/, 'no positions'],
		[[vcf => $v, weights => [{ map { $_ => $wrow->{$_} } qw(chr_name chr_position effect_allele effect_weight) }], match => 'id'],
			qr/^pgs_score: "weights" has no "rsID" column/, 'no rsID to match by ID'],
		[[vcf => $v, weights => [{ %{ $wrow }, effect_weight => 'NaN' }]], qr/^pgs_score: row 1 of "weights" has effect_weight "NaN", which isn't a number/, 'a NaN weight'],
		[[vcf => $v, weights => [$wrow, { %{ $wrow }, effect_weight => undef }]], qr/^pgs_score: row 2 of "weights" has effect_weight undef, which isn't a number/, 'an undef weight'],
		[[vcf => $v, weights => { k => { %{ $wrow }, effect_weight => '0,5' } }], qr/^pgs_score: row "k" of "weights" has effect_weight "0,5", which isn't a number/, 'a decimal comma, in a hoh'],
		[[vcf => $v, weights => [{ %{ $wrow }, effect_allele => '' }]], qr/^pgs_score: row 1 of "weights" has no effect_allele/, 'an empty effect allele'],
		[[vcf => $v, weights => [{ %{ $wrow }, chr_position => '1e5' }]], qr/^pgs_score: row 1 of "weights" has chr_position "1e5", which isn't a position/, 'a position in exponent form'],
		[[vcf => $v, weights => [{ %{ $wrow }, hm_chr => '1', hm_pos => 0 }]], qr/^pgs_score: row 1 of "weights" has hm_pos "0", which isn't a position/, 'position 0, in hm_pos'],
		[[vcf => $v, weights => [{ %{ $wrow }, is_dominant => 'True' }]], qr/^pgs_score: row 1 of "weights" has is_dominant "True", and only additive weights can be scored/, 'a dominant weight'],
		[[vcf => $v, weights => [{ %{ $wrow }, is_haplotype => 'TRUE' }]], qr/^pgs_score: row 1 of "weights" has is_haplotype "TRUE", and only/, 'a haplotype'],
		[[vcf => $v, weights => [{ %{ $wrow }, dosage_1_weight => '0.2' }]], qr/^pgs_score: row 1 of "weights" has a dosage_1_weight, and only additive weights can be scored/, 'a dosage weight'],
		[[vcf => { k => { CHROM => 'chr1', POS => 5, REF => 'A', ALT => 'G', FORMAT => 'GT', S => '0/1' } }, weights => $w],
			qr/^pgs_score: "vcf" has no <sample>\.GT column; read the VCF with read_table's default of explode => 1/, 'a VCF read with explode => 0'],
		[[vcf => { k => { CHROM => 'chr1', POS => 5, ALT => 'G', 'S.GT' => '0/1' } }, weights => $w], qr/^pgs_score: "vcf" has no "REF" column/, 'no REF'],
		[[vcf => { k => { CHROM => 'chr1', POS => 5, REF => 'A', ALT => 'G', 'S.GT' => '0/1' } }, weights => $w, match => 'id'], qr/^pgs_score: "vcf" has no "ID" column/, 'no ID to match by ID'],
	);
	foreach my $c (@cases) {
		my ($args, $re, $name) = @{ $c };
		dies_like sub { pgs_score(@{ $args }) }, $re, "pgs_score: $name dies";
	}
	SKIP: {
		skip 'PLINK 2, which Alien::Bioinf installs, is not installed', 3 unless $plink2;
		is pgs_score(vcf => $v, weights => $w, match => 'id')->{S}{score}, 0.5, 'pgs_score: rs1 is found by ID';
		is pgs_score(vcf => $v, weights => [{ %{ $wrow }, chr_name => undef, chr_position => undef }], match => 'id')->{S}{score}, 0.5,
			'pgs_score: matching by ID needs no position';
		is pgs_score(vcf => $v, weights => [{ %{ $wrow }, rsID => undef }])->{S}{score}, 0.5, 'pgs_score: matching by position needs no rsID';
	}
}

# ---- sex_check --------------------------------------------------------------

# Added with sex_check on 2026-10-18; t/sex.t has the counts themselves, and the
# errors that come from reading the file.
{
	my $v = spew('sex.vcf', "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchrX\t3000000\t.\tA\tG\t.\t.\t.\tGT\t0/1\n");
	is sex_check(vcf => $v, report => 0)->{S}{sex}, 'female', 'sex_check: the arguments below differ from these only where they say';
	is_deeply sex_check($v, report => 0), sex_check(vcf => $v, report => 0), 'sex_check: the file may come first, unnamed, before the options';
	{
		open my $save, '>&', \*STDOUT or die $!;
		open STDOUT, '>', File::Spec->devnull or die $!;
		my $r = eval { sex_check($v) };
		open STDOUT, '>&', $save or die $!;
		is_deeply $r, { S => sex_check(vcf => $v, report => 0)->{S} }, 'sex_check: or alone, as sex_check($file), rather than dying on an odd number of arguments';
	}
	my @cases = (
		[[], qr/^sex_check: needs "vcf"/, 'no vcf'],
		[[{ vcf => $v }], qr/^sex_check: takes name => value pairs, not a hash ref/, 'a hash ref'],
		[[vcf => $v, depth => 7], qr/^sex_check: doesn't know "depth"; it takes cutoff, min_depth, report, vcf/, 'a misspelled option'],
		[[vcf => $v, cutoff => 0], qr/^sex_check: "cutoff" must be a number above 0, not "0"/, 'cutoff 0'],
		[[vcf => $v, cutoff => 'NaN'], qr/^sex_check: "cutoff" must be a number above 0, not "NaN"/, 'cutoff NaN'],
		[[vcf => $v, cutoff => [0.6]], qr/^sex_check: "cutoff" must be a number above 0, not an ARRAY ref/, 'cutoff an array ref'],
		[[vcf => $v, min_depth => -1], qr/^sex_check: "min_depth" must be a whole number, 0 or more, not "-1"/, 'min_depth -1'],
		[[vcf => $v, min_depth => 7.5], qr/^sex_check: "min_depth" must be a whole number, 0 or more, not "7\.5"/, 'min_depth 7.5'],
		[[vcf => $v, report => 2], qr/^sex_check: "report" must be 0 or 1, not "2"/, 'report 2'],
		[['cutoff', 0.5, 'min_depth'], qr/^sex_check: takes name => value pairs, and was given an odd number of arguments/, 'an odd list that starts with an option'],
		[[$v, 'report'], qr/^sex_check: ".+sex\.vcf" is the file, so the options after it must be name => value pairs, and are an odd number/, 'the file and an option with no value, which is read as a pair'],
		[[vcf => $v, report => 'yes'], qr/^sex_check: "report" must be 0 or 1, not "yes"/, 'report "yes"'],
		[[vcf => \'calls'], qr/^sex_check: "vcf" must be a file name, not a SCALAR ref/, 'a scalar ref'],
		[[vcf => catfile($dir, 'none.vcf.bz2')], qr/^sex_check: "vcf" file ".+none\.vcf\.bz2" doesn't exist or isn't readable/, 'a VCF that is not there'],
		[[vcf => $dir], qr/^sex_check: "vcf" file ".+" doesn't exist or isn't readable/, 'a directory'],
		[[vcf => spew('sex.bcf', '')], qr/^sex_check: "vcf" file .+sex\.bcf must end in \.vcf, \.vcf\.gz, \.vcf\.bgz or \.vcf\.bz2/, 'a BCF'],
		[[vcf => spew('sex.vcf.zst', '')], qr/^sex_check: "vcf" file .+sex\.vcf\.zst must end in/, 'a .vcf.zst'],
	);
	foreach my $c (@cases) {
		my ($args, $re, $name) = @{ $c };
		dies_like sub { sex_check(@{ $args }) }, $re, "sex_check: $name dies";
	}
}

# ---- plot_msa ---------------------------------------------------------------

my %two = (fasta => { a => 'MK', b => 'MV' }, filename => catfile($dir, 'x.png'), query => 'a');
foreach my $n (0, -1, 1.5, 'two') {
	dies_like sub { plot_msa(%two, 'active_site_aa' => { X => $n }) },
		qr/^plot_msa: active site X must be a residue number, 1 or more, not "\Q$n\E"/, "active site $n dies, rather than counting from the end";
}
dies_like sub { plot_msa(%two, 'active_site_aa' => [1]) }, qr/^plot_msa: "active_site_aa" must be a hash ref/, 'active_site_aa must be a hash ref';

# Up to 2026-10-17 a wrong kind of reference died with perl's own "Not an
# ARRAY reference" or "Can't use string", naming no function.
my %two_noq = (fasta => { a => 'MK', b => 'MV' }, filename => catfile($dir, 'x.png'));
foreach my $case (
	[sub { plot_msa(%two_noq, order => 'a') },        qr/^plot_msa: "order" must be an array ref, not "a"/],
	[sub { plot_msa(%two_noq, labels => ['a']) },     qr/^plot_msa: "labels" must be a hash ref, not an ARRAY ref/],
	[sub { plot_msa(%two_noq, 'clustal_args' => '--iter=2') }, qr/^plot_msa: "clustal_args" must be an array ref, not "--iter=2"/],
	[sub { plot_phylo(fasta => { a => 'MK', b => 'MV' }, labels => 'a') }, qr/^plot_phylo: "labels" must be a hash ref, not "a"/],
	[sub { plot_phylo(fasta => { a => 'MK', b => 'MV' }, 'clustal_args' => {}) }, qr/^plot_phylo: "clustal_args" must be an array ref, not a HASH ref/],
) {
	dies_like $case->[0], $case->[1], "dies: $case->[1]";
}

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
dies_like sub { msa_quality_table(%q, order => 'A.b') }, qr/^msa_quality_table: "order" must be an array ref, not "A\.b"/, 'an order that is not an array ref dies';
dies_like sub { msa_quality_table(%q, 'alignment_json' => []) },
	qr/^msa_quality_table: "alignment_json" must be a file name or a hash ref, not an ARRAY ref/, 'and so does an alignment_json that is neither';
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
dies_like sub { clustal_view_residues(%view, order => 'e{f}') }, qr/^clustal_view_residues: "order" must be an array ref, not "e\{f\}"/,
	'an order that is not an array ref dies';
dies_like sub { clustal_view_residues(%view, 'color_residues' => { 'e{f}' => [2] }) },
	qr/^clustal_view_residues: "color_residues" must map e\{f\} to a hash ref of residue number => colour, not an ARRAY ref/,
	'and so does a protein in color_residues that maps to anything but a hash ref';
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
