require 5.012;
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

#  loading

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

#  _run

# Up to 2026-10-07 _run read only the exit code, the high byte of $?, so a
# program killed by a signal was reported as "exit 0".
dies_like sub { Cn::Bioinf::_run($^X, '-e', 'exit 3') }, qr/^_run: .+ failed: exit 3 at /, 'a non-zero exit is reported';
SKIP: {
	skip 'MSWin32 has no POSIX signals', 1 if $^O eq 'MSWin32';
	# set to DEFAULT in the child, which may have inherited TERM ignored
	dies_like sub { Cn::Bioinf::_run($^X, '-e', '$SIG{TERM} = q{DEFAULT}; kill q{TERM}, $$; sleep 5') },
		qr/^_run: .+ failed: killed by signal 15/, 'a death by signal is reported as one, not as "exit 0"';
}

#  _plink2

# A plink2 that never starts is reported with the reason system gave. Up to
# 2026-11-02 $! was read after STDOUT and STDERR had been reopened, and the
# reason given was theirs ("Inappropriate ioctl for device"). The reason
# expected is the one perl gives for the same failed system, which is how this
# perl words it.
SKIP: {
	my $nowhere = catfile($dir, 'no such dir', 'plink2');
	my $want = do { no warnings 'exec'; system($nowhere) == -1 ? "$!" : undef };
	skip "this perl's system does not return -1 for a program that is not there", 1 unless defined $want;
	dies_like sub { Cn::Bioinf::_plink2('pgs_score', $nowhere, $dir, '--version') },
		qr/^pgs_score: .+ failed \(it didn't start: \Q$want\E\)/, '_plink2: a plink2 that does not start is reported with the reason it did not';
}

#  undef file names

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

#  hash2fasta_file 

# Up to 2026-10-07 the width went straight to the XS as a UV: -1 wrapped round
# to one line per sequence, and 'x' warned with no function name.
foreach my $w (-1, 1.5, 'x') {
	dies_like sub { hash2fasta_file({ a => 'MK' }, catfile($dir, 'w.fa'), undef, $w) },
		qr/^hash2fasta_file: 4th argument must be a whole number, 0 or more, or undef, not "\Q$w\E"/, "width $w dies";
}

#  get_best_alignment_hit 

my $report = '{"report":{"results":{"search":{"query_title":"q","hits":[]}}}}';
dies_like sub { get_best_alignment_hit(spew('twice.json', qq({"BlastOutput2":[$report,$report]}))) },
	qr/^get_best_alignment_hit: two queries in .+twice\.json are both titled "q"/, 'two queries of one title die rather than one being lost';
# A file that is not BLAST's JSON used to die with JSON::MaybeXS's message or
# perl's "Not a HASH reference", naming no function, or, for {}, return no
# queries as though BLAST had found none.
foreach my $case (
	['', qr/^get_best_alignment_hit: .+empty\.json isn't JSON: /, 'an empty file'],
	['{not json', qr/^get_best_alignment_hit: .+bad\.json isn't JSON: /, 'a file that is not JSON'],
	['[1]', qr/^get_best_alignment_hit: .+array\.json isn't a BLAST report in JSON \(-outfmt 15\): it has no "BlastOutput2" array/, 'a JSON array'],
	['{}', qr/^get_best_alignment_hit: .+object\.json isn't a BLAST report in JSON \(-outfmt 15\): it has no "BlastOutput2" array/, 'a JSON object of no report'],
	['{"BlastOutput2":[1]}', qr/^get_best_alignment_hit: report 1 of .+nosearch\.json has no report\.results\.search/, 'a report with no search'],
	['{"BlastOutput2":[{"report":{"results":{"search":{"hits":[]}}}}]}', qr/^get_best_alignment_hit: report 1 of .+noquery\.json names no query: it has no query_title or query_id/,
		'a report with no query'],
) {
	my ($json, $re, $name) = @{ $case };
	my ($file) = $re =~ /\.\+(\w+)\\\.json/;
	dies_like sub { get_best_alignment_hit(spew("$file.json", $json)) }, $re, "get_best_alignment_hit: $name dies";
}
# A query with no title is named by its query_id, where it used to be named
# "" with an "uninitialized" warning.
is_deeply [keys %{ get_best_alignment_hit(spew('untitled.json', '{"BlastOutput2":[{"report":{"results":{"search":{"query_id":"Query_1","hits":[]}}}}]}')) }],
	['Query_1'], 'get_best_alignment_hit: a query with no title is named by its query_id';

#  pgs_score 

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

#  sex_check 

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

#  relatedness 

# Added with relatedness on 2026-10-19; t/relatedness.t has the counts
# themselves, and the errors that come from reading the files. Each case here
# dies before any file is opened.
{
	my $hdr = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	my $two = spew('two.vcf', "$hdr\tA\tB\n1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\n");
	my $a = spew('a.vcf', "$hdr\tA\n1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1\n");
	my $b = spew('b.vcf', "$hdr\tB\n1\t6\t.\tA\tG\t.\t.\t.\tGT\t0/1\n");
	my @w;
	{
		local $SIG{__WARN__} = sub { push @w, $_[0] };
		is_deeply [sort keys %{ relatedness($two, report => 0) }], ['A', 'B'], 'relatedness: the arguments below differ from these only where they say';
		is relatedness($a, $b, report => 0)->{A}{B}{sites}, 2, 'relatedness: and two own VCFs';
	}
	is scalar(@w), 2, 'relatedness: which are undecided, with too few sites, and say so';
	my $none = catfile($dir, 'none.vcf');
	my @cases = (
		[[], qr/^relatedness: needs a VCF of several samples, or several VCFs of one/, 'nothing'],
		[[report => 0], qr/^relatedness: needs a VCF/, 'options and no VCF'],
		[[{ vcf => $two }], qr/^relatedness: takes name => value pairs, not a hash ref/, 'a hash ref'],
		[[$two, 'report'], qr/^relatedness: takes name => value pairs, and was given an odd number of arguments/, 'the file and an option with no value'],
		[[$two, depth => 7], qr/^relatedness: doesn't know "depth"; it takes max_ibs0, min_depth, report, vcf/, 'a misspelled option'],
		[[$a, vcf => $b], qr/^relatedness: the VCFs are given either first, unnamed, or as "vcf", not both/, 'files both first and as "vcf"'],
		[[vcf => []], qr/^relatedness: needs a VCF/, 'an empty array as "vcf"'],
		[[vcf => { a => $a }], qr/^relatedness: "vcf" must be a file name, not a HASH ref/, 'a hash ref as "vcf"'],
		[[vcf => [$a, undef]], qr/^relatedness: "vcf" must be a file name, not undef/, 'an undef among the files'],
		[[vcf => [$a, \'b']], qr/^relatedness: "vcf" must be a file name, not a SCALAR ref/, 'a scalar ref among the files'],
		[[$a, $none], qr/^relatedness: "vcf" file "\Q$none\E" doesn't exist or isn't readable/, 'a file that is not there'],
		[[$a, $dir], qr/^relatedness: "vcf" file ".+" doesn't exist or isn't readable/, 'a directory'],
		[[$a, spew('b.bcf', '')], qr/^relatedness: "vcf" file .+b\.bcf must end in \.vcf, \.vcf\.gz, \.vcf\.bgz or \.vcf\.bz2/, 'a BCF'],
		[[$a, $b, $a], qr/^relatedness: "\Q$a\E" is given twice/, 'one file twice'],
		[[$two, min_depth => -1], qr/^relatedness: "min_depth" must be a whole number, 0 or more, not "-1"/, 'min_depth -1'],
		[[$two, min_depth => 7.5], qr/^relatedness: "min_depth" must be a whole number, 0 or more, not "7\.5"/, 'min_depth 7.5'],
		[[$two, max_ibs0 => 0], qr/^relatedness: "max_ibs0" must be a number above 0, not "0"/, 'max_ibs0 0'],
		[[$two, max_ibs0 => 'NaN'], qr/^relatedness: "max_ibs0" must be a number above 0, not "NaN"/, 'max_ibs0 NaN'],
		[[$two, max_ibs0 => [0.01]], qr/^relatedness: "max_ibs0" must be a number above 0, not an ARRAY ref/, 'max_ibs0 an array ref'],
		[[$two, report => 'yes'], qr/^relatedness: "report" must be 0 or 1, not "yes"/, 'report "yes"'],
	);
	foreach my $c (@cases) {
		my ($args, $re, $name) = @{ $c };
		dies_like sub { relatedness(@{ $args }) }, $re, "relatedness: $name dies";
	}
}

# The XSUB's own checks of what perl hands it, which relatedness never gets
# wrong; they named XS_Cn__Bioinf__relatedness_own, __func__ in an XSUB, up to
# 2026-11-02.
{
	dies_like sub { Cn::Bioinf::_relatedness_own([sub { undef }], [], 0) },
		qr/^_relatedness_own: 1 readers for 0 files/, '_relatedness_own: names itself when its arrays differ in length';
	my @sparse;
	$sparse[1] = 'b.vcf'; # element 0 does not exist
	dies_like sub { Cn::Bioinf::_relatedness_own([sub { undef }, sub { undef }], \@sparse, 0) },
		qr/^_relatedness_own: file 0 has no reader/, '_relatedness_own: and when a file is missing';
}

#  concordance 

# Added with concordance on 2026-10-19; t/concordance.t has the counts
# themselves, and the errors that come from reading the files. Each case here
# dies before any file is opened.
{
	my $hdr = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	my $truth = spew('conc.truth.vcf', "$hdr\tS\n1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1\n");
	my $test = spew('conc.test.vcf', "$hdr\tS\n1\t5\t.\tA\tG\t.\t.\t.\tGT\t1/1\n");
	is concordance(truth => $truth, test => $test, report => 0)->{S}{counts}{discordant_het_ref_to_hom_alt}, 1,
		'concordance: the arguments below differ from these only where they say';
	my %ok = (truth => $truth, test => $test);
	my $none = catfile($dir, 'none.vcf');
	my @cases = (
		[[], qr/^concordance: needs "test", "truth"/, 'nothing'],
		[[truth => $truth], qr/^concordance: needs "test"/, 'no test'],
		[[test => $test], qr/^concordance: needs "truth"/, 'no truth'],
		[[{ %ok }], qr/^concordance: takes name => value pairs, not a hash ref/, 'a hash ref'],
		[[%ok, 'report'], qr/^concordance: takes name => value pairs, and was given an odd number of arguments/, 'an option with no value'],
		[[%ok, depth => 7], qr/^concordance: doesn't know "depth"; it takes discordances, min_depth, report, test, truth/, 'a misspelled option'],
		[[%ok, truth => [$truth]], qr/^concordance: "truth" must be a file name, not an ARRAY ref/, 'an array ref as "truth"'],
		[[%ok, test => \'x'], qr/^concordance: "test" must be a file name, not a SCALAR ref/, 'a scalar ref as "test"'],
		[[%ok, test => $none], qr/^concordance: "test" file "\Q$none\E" doesn't exist or isn't readable/, 'a test that is not there'],
		[[%ok, truth => $dir], qr/^concordance: "truth" file ".+" doesn't exist or isn't readable/, 'a directory as "truth"'],
		[[%ok, truth => spew('conc.bcf', '')], qr/^concordance: "truth" file .+conc\.bcf must end in \.vcf, \.vcf\.gz, \.vcf\.bgz or \.vcf\.bz2/, 'a BCF'],
		[[%ok, min_depth => -1], qr/^concordance: "min_depth" must be a whole number, 0 or more, not "-1"/, 'min_depth -1'],
		[[%ok, min_depth => 7.5], qr/^concordance: "min_depth" must be a whole number, 0 or more, not "7\.5"/, 'min_depth 7.5'],
		[[%ok, min_depth => [7]], qr/^concordance: "min_depth" must be a whole number, 0 or more, not an ARRAY ref/, 'min_depth an array ref'],
		[[%ok, discordances => 2], qr/^concordance: "discordances" must be 0 or 1, not "2"/, 'discordances 2'],
		[[%ok, discordances => 'yes'], qr/^concordance: "discordances" must be 0 or 1, not "yes"/, 'discordances "yes"'],
		[[%ok, report => 'yes'], qr/^concordance: "report" must be 0 or 1, not "yes"/, 'report "yes"'],
	);
	foreach my $c (@cases) {
		my ($args, $re, $name) = @{ $c };
		dies_like sub { concordance(@{ $args }) }, $re, "concordance: $name dies";
	}
}

#  plot_msa 

my %two = (fasta => { a => 'MK', b => 'MV' }, filename => catfile($dir, 'x.png'), query => 'a');
foreach my $n (0, -1, 1.5, 'two') {
	dies_like sub { plot_msa(%two, 'active_site_aa' => { X => $n }) },
		qr/^plot_msa: active site X must be a residue number, 1 or more, not "\Q$n\E"/, "active site $n dies, rather than counting from the end";
}
dies_like sub { plot_msa(%two, 'active_site_aa' => [1]) }, qr/^plot_msa: "active_site_aa" must be a hash ref/, 'active_site_aa must be a hash ref';
dies_like sub { plot_msa(%two, tooltips => 'yes') }, qr/^plot_msa: "tooltips" must be 0 or 1, not "yes"/, 'tooltips "yes" dies';
dies_like sub { plot_msa(%two, tooltips => 1) }, qr/^plot_msa: "tooltips" needs an SVG "filename", the only format that carries them, not "\Q$two{filename}\E"/,
	'tooltips in a PNG dies';

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
# A sequence that is undef, a reference or empty is refused before clustalo
# runs. Up to 2026-11-02 undef warned "uninitialized" and then died in
# write_fasta, a reference was written as "ARRAY(0x...)" for clustalo to
# align, and an empty one died in clustalo ("Sequence 's2' has 0 residues"),
# none of them naming plot_msa; and a fasta that was neither a file nor a hash
# ref named _ungapped_fasta.
foreach my $case (
	[{ a => 'MK', b => undef }, qr/^plot_msa: "fasta" sequence "b" must be a string of residues, not undef/, 'undef'],
	[{ a => 'MK', b => [1] }, qr/^plot_msa: "fasta" sequence "b" must be a string of residues, not an ARRAY ref/, 'an ARRAY ref'],
	[{ a => 'MK', b => '' }, qr/^plot_msa: "fasta" sequence "b" has no residues/, 'empty'],
	[{ a => 'MK', b => '--.-' }, qr/^plot_msa: "fasta" sequence "b" has no residues/, 'only gaps'],
	[[], qr/^plot_msa: "fasta" must be a FASTA file name or a hash ref of name => sequence, not an ARRAY ref/, 'an ARRAY ref for fasta'],
) {
	my ($fasta, $re, $name) = @{ $case };
	my $warned = '';
	local $SIG{__WARN__} = sub { $warned .= $_[0] };
	dies_like sub { plot_msa(%two_noq, fasta => $fasta) }, $re, "plot_msa: a sequence that is $name dies";
	is $warned, '', "plot_msa: with no warning first, for a sequence that is $name";
}
dies_like sub { msa_quality_table(fasta => { a => 'MK', b => '' }, filename => catfile($dir, 'q.png'), metric => 'evalue') },
	qr/^msa_quality_table: "fasta" sequence "b" has no residues/, 'msa_quality_table: so does an empty sequence';

# "fh" and "execute", as Matplotlib::Simple's plt takes them, are checked
# before clustalo would run.
open my $plain, '>', catfile($dir, 'plain.py') or die $!;
foreach my $f ([\&plot_msa, 'plot_msa', %two_noq], [\&plot_phylo, 'plot_phylo', fasta => $two_noq{fasta}]) {
	my ($code, $name, @args) = @{ $f };
	foreach my $case (
		[[execute => 'yes'], qr/^$name: "execute" must be 0 or 1, not "yes"/],
		[[execute => 2],     qr/^$name: "execute" must be 0 or 1, not "2"/],
		[[execute => []],    qr/^$name: "execute" must be 0 or 1, not an ARRAY ref/],
		[[fh => $plain],     qr/^$name: "fh" must be a File::Temp object, as for Matplotlib::Simple, not a GLOB ref/],
		[[fh => 'x.py'],     qr/^$name: "fh" must be a File::Temp object, as for Matplotlib::Simple, not "x\.py"/],
	) {
		dies_like sub { $code->(@args, @{ $case->[0] }) }, $case->[1], "dies: $case->[1]";
	}
}
close $plain;

#  msa_quality_table 

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
# a file of JSON that is not a BLAST report used to die with perl's "Not a HASH
# reference", and one that is not JSON at all with JSON::MaybeXS's message
dies_like sub { msa_quality_table(%q, 'alignment_json' => spew('q_array.json', '[1]')) },
	qr/^msa_quality_table: "alignment_json" isn't a BLAST report in JSON \(-outfmt 15\): it has no "BlastOutput2" array/, 'an alignment_json of a JSON array dies';
dies_like sub { msa_quality_table(%q, 'alignment_json' => spew('q_bad.json', '{not json')) },
	qr/^msa_quality_table: .+q_bad\.json isn't JSON: /, 'and one that is not JSON';

#  clustal_view_residues 

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
{
	# Up to 2026-11-02 |, < and > were left as they were, and print under
	# LaTeX's default OT1 fonts as an em dash and an inverted ! and ?; and the
	# residues were not escaped at all, so a % in one ended the line.
	local *STDOUT; # for the "wrote" lines, which are not TAP
	open STDOUT, '>', catfile($dir, 'quiet.txt') or die "can't send STDOUT to quiet.txt: $!";
	my $odd = hash2fasta_file({ 'sp|P1|<x>' => 'MK%L#' }, catfile($dir, 'odd.aln.fa'));
	clustal_view_residues('msa_file' => $odd, 'output_tex_file' => $tex);
	like slurp($tex), qr/^\\textit\{sp\\textbar\{\}P1\\textbar\{\}\\textless\{\}x\\textgreater\{\}\} & \\texttt\{MK\\%L\\#\} & 3$/m,
		'a name\'s |, < and > print as themselves, and so do a residue\'s % and #';
}
{
	# A number in a block's last columns used to be left off the track row,
	# and so did every number after it in that block; it now ends at the
	# block's end. 9 and 10 in the second case cannot both fit, so 10 is left
	# out and 9 kept.
	local *STDOUT; # for the "wrote" lines, which are not TAP
	open STDOUT, '>', catfile($dir, 'quiet.txt') or die "can't send STDOUT to quiet.txt: $!";
	my $one = hash2fasta_file({ a => 'ACDEFGHIKLMN', b => 'ACDEFGHIKLMN' }, catfile($dir, 'track.aln.fa'));
	my %t = ('msa_file' => $one, 'output_tex_file' => $tex, 'row_width' => 10, track => 'a');
	clustal_view_residues(%t, 'color_residues' => { a => { 3 => 'red', 5 => 'blue', 10 => 'green', 11 => 'red' } });
	my @track = slurp($tex) =~ /^\\textit\{a\} track & \\texttt\{([^}]*)\}/mg;
	is_deeply \@track, ['--3-5---10', '11'], 'clustal_view_residues: a residue number in the last column of a block ends at the block\'s end';
	clustal_view_residues(%t, 'color_residues' => { a => { 9 => 'red', 10 => 'green' } });
	@track = slurp($tex) =~ /^\\textit\{a\} track & \\texttt\{([^}]*)\}/mg;
	is_deeply \@track, ['--------9-', '--'], 'clustal_view_residues: and one that would overlap the number before it is left out';
	# A tracked protein that "order" leaves out has no track row, and used to be
	# counted as having one, so a table held more blocks than "split" said.
	clustal_view_residues(%t, track => 'b', order => ['a'], 'row_width' => 4, split => 2);
	my $tables = () = slurp($tex) =~ /\\begin\{table\}/g;
	is $tables, 2, 'clustal_view_residues: three blocks at 2 a table are two tables when the tracked protein is not shown';
	# a residue number that is not a number used to warn "isn't numeric in
	# sort" before it was refused
	my $warned = '';
	local $SIG{__WARN__} = sub { $warned .= $_[0] };
	dies_like sub { clustal_view_residues(%t, 'color_residues' => { a => { abc => 'red', 2 => 'blue' } }) },
		qr/^clustal_view_residues: a residue of a in "color_residues" must be a residue number, 1 or more, not "abc"/,
		'clustal_view_residues: a residue number "abc" dies';
	is $warned, '', 'clustal_view_residues: with no warning from sorting it first';
}

# ---- hirisplex_genotypes and trait_predict
#
# Neither needs any of Alien::Bioinf's programs, so every one of these runs
# wherever the module builds. The VCF below is a well-formed one that holds
# none of the panel's markers; what is being checked is the arguments, which
# are refused before it is opened.
{
	my $vcf = spew('checks.vcf', "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
		. join("\t", '#CHROM', qw(POS ID REF ALT QUAL FILTER INFO FORMAT S)) . "\n"
		. join("\t", 1, 100, 'rs1', 'A', 'G', '.', '.', '.', 'GT', '0/1') . "\n");
	my $bam = spew('checks.bam', "not a VCF\n");
	my $model = spew('checks_model.tsv', "trait\tcategory\tmarker\tbeta\neye\tblue\t(intercept)\t1\neye\tbrown\t(intercept)\t0\n");
	my @bad = (
		[['nowhere.vcf'], qr/^hirisplex_genotypes: "vcf" file "nowhere\.vcf" doesn't exist or isn't readable/, 'a missing file'],
		[[vcf => $vcf, match => 'rsid'], qr/^hirisplex_genotypes: "match" must be "position" or "id", not "rsid"/, 'an unknown match'],
		[[vcf => $vcf, build => 'hg18'], qr/^hirisplex_genotypes: "build" must be hg19, b37, hg38 or b38, not "hg18"/, 'an unknown build'],
		[[vcf => $vcf, min_depth => -1], qr/^hirisplex_genotypes: "min_depth" must be a whole number, 0 or more, not "-1"/, 'a negative min_depth'],
		[[vcf => $vcf, min_depth => '7.5'], qr/^hirisplex_genotypes: "min_depth" must be a whole number, 0 or more, not "7\.5"/, 'a fractional min_depth'],
		[[vcf => $vcf, report => 2], qr/^hirisplex_genotypes: "report" must be 0 or 1, not "2"/, 'a report of 2'],
		[[vcf => $vcf, file => []], qr/^hirisplex_genotypes: "file" must be a file name, not an ARRAY ref/, 'a file that is a ref'],
		[[vcf => $vcf, colour => 1], qr/^hirisplex_genotypes: doesn't know "colour"; it takes build, file, match, min_depth, report, vcf/, 'an unknown option'],
		[[vcf => $bam], qr/^hirisplex_genotypes: "vcf" file \Q$bam\E must end in \.vcf, \.vcf\.gz, \.vcf\.bgz or \.vcf\.bz2/, 'a file whose name says it is not a VCF'],
	);
	foreach my $case (@bad) {
		my ($args, $re, $name) = @{ $case };
		dies_like sub { hirisplex_genotypes(@{ $args }) }, $re, "hirisplex_genotypes: $name dies";
	}
	# A VCF whose header gives a length for none of the chromosomes the build
	# is told by, and one whose lengths are neither build's.
	my $nolen = spew('nolen.vcf', "##fileformat=VCFv4.2\n" . join("\t", '#CHROM', qw(POS ID REF ALT QUAL FILTER INFO FORMAT S)) . "\n");
	dies_like sub { hirisplex_genotypes(vcf => $nolen, report => 0) },
		qr/^hirisplex_genotypes: "vcf" .* has no ##contig line giving a length for chromosome 1, 15, 16, 20/,
		'hirisplex_genotypes: a header with no length to tell the build by dies';
	my $odd = spew('oddlen.vcf', "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n"
		. join("\t", '#CHROM', qw(POS ID REF ALT QUAL FILTER INFO FORMAT S)) . "\n");
	dies_like sub { hirisplex_genotypes(vcf => $odd, report => 0) },
		qr/^hirisplex_genotypes: "vcf" .* gives chromosome lengths that are neither GRCh37's nor GRCh38's \(1 is 1000\)/,
		'hirisplex_genotypes: a header whose lengths are neither build\'s dies';
	# Matching by rs number needs no build, so a header with no lengths is
	# enough for it; the file holds no marker, so every one is missing.
	{
		my $g;
		my $warned = '';
		local $SIG{__WARN__} = sub { $warned .= $_[0] };
		$g = hirisplex_genotypes(vcf => $vcf, match => 'id', report => 0);
		is $g->{S}{called}, 0, 'hirisplex_genotypes: a VCF with none of the markers calls none of them';
		is $g->{S}{markers}, 41, 'hirisplex_genotypes: and still counts the whole panel';
		is $g->{S}{eye}, 0, 'hirisplex_genotypes: with no HERC2 marker, no eye prediction is allowed';
		like $warned, qr/^hirisplex_genotypes: "vcf" .* holds none of the 41 HIrisPlex-S markers/,
			'hirisplex_genotypes: and it says so, rather than returning an empty answer quietly';
	}
	# The panel's insertion is checked to count one base, as every other marker
	# is; with none it used to be stored as base 0 and never counted.
	dies_like sub { Cn::Bioinf::_hirisplex_dosage(sub { undef }, 'x.vcf', 0, [[16, 89985751, '', 1, 'rs312262906', 0]], 1, 0) },
		qr/^hps_panel: entry 0 of the panel counts "", which isn't one base/, 'hirisplex_genotypes: a panel insertion of no base dies';
	dies_like sub { Cn::Bioinf::_hirisplex_dosage(sub { undef }, 'x.vcf', 0, [[16, 89985751, 'A', 1, 'rs312262906', 0]], 1, 0) },
		qr/^vcf_read: x\.vcf has no #CHROM line/, 'hirisplex_genotypes: and one of a base gets as far as reading the VCF';

	# a model that uses one marker, for checking the dosages given for it
	my $dose_model = [{ trait => 'eye', category => 'blue', marker => 'rs12913832_T', beta => 1 },
		{ trait => 'eye', category => 'brown', marker => '(intercept)', beta => 0 }];
	my @bad_model = (
		[[genotypes => 'x', model => $model], qr/^trait_predict: "genotypes" must be a hash ref of sample => genotypes, as hirisplex_genotypes returns, not "x"/, 'genotypes that are a string'],
		[[genotypes => {}, model => $model], qr/^trait_predict: "genotypes" has no samples/, 'no samples'],
		[[genotypes => { S => 'x' }, model => $model], qr/^trait_predict: "genotypes" sample "S" must be a hash ref, not "x"/, 'a sample that is a string'],
		[[genotypes => { S => {} }, model => 'nowhere.tsv'], qr/^trait_predict: "model" file "nowhere\.tsv" doesn't exist or isn't readable/, 'a missing model'],
		[[genotypes => { S => {} }, model => [{ trait => 'eye' }]], qr/^trait_predict: "model" has no "category", "marker", "beta" columns; a model table has trait, category, marker and beta/, 'a model missing columns'],
		[[genotypes => { S => {} }, model => [{ trait => 'eye', category => 'blue', marker => 'm', beta => 'x' }]],
			qr/^trait_predict: row 1 of "model" has beta "x", which isn't a number/, 'a beta that is not a number'],
		[[genotypes => { S => {} }, model => [{ trait => 'eye', category => 'blue', marker => '', beta => 1 }]],
			qr/^trait_predict: row 1 of "model" has no marker/, 'a row with no marker'],
		[[genotypes => { S => {} }, model => [{ trait => 'eye', category => 'blue', marker => 'm', beta => 1 },
			{ trait => 'eye', category => 'blue', marker => 'm', beta => 2 }]],
			qr/^trait_predict: row 2 of "model" gives eye category "blue" a second beta for m/, 'a marker given two betas'],
		[[genotypes => { S => {} }, model => $model, force => 2], qr/^trait_predict: "force" must be 0 or 1, not "2"/, 'a force of 2'],
		# a dosage used to be multiplied in whatever it held, so that "NA" was a
		# called genotype of 0 copies, with only perl's "isn't numeric" warning
		[[genotypes => { S => { 'rs12913832_T' => 'NA' } }, model => $dose_model],
			qr/^trait_predict: "genotypes" sample "S" has rs12913832_T "NA", which isn't a dosage from 0 to 2; a missing genotype is undef/,
			'a dosage of "NA", rather than counting it 0 copies'],
		[[genotypes => { S => { genotype => { 'rs12913832_T' => 3 } } }, model => $dose_model],
			qr/^trait_predict: "genotypes" sample "S" has rs12913832_T "3", which isn't a dosage from 0 to 2/, 'a dosage of 3'],
		[[genotypes => { S => { 'rs12913832_T' => -1 } }, model => $dose_model],
			qr/^trait_predict: "genotypes" sample "S" has rs12913832_T "-1", which isn't a dosage/, 'a dosage of -1'],
		[[genotypes => { S => { 'rs12913832_T' => [] } }, model => $dose_model],
			qr/^trait_predict: "genotypes" sample "S" has rs12913832_T an ARRAY ref, which isn't a dosage/, 'a dosage that is a ref'],
		[[genotypes => { S => {} }, model => $model, shade => 1], qr/^trait_predict: doesn't know "shade"; it takes force, genotypes, model, report/, 'an unknown option'],
	);
	foreach my $case (@bad_model) {
		my ($args, $re, $name) = @{ $case };
		dies_like sub { trait_predict(@{ $args }) }, $re, "trait_predict: $name dies";
	}
	# A model of one category is still a model: its probability is 1.
	my $one = trait_predict(genotypes => { S => { m => 2 } },
		model => [{ trait => 't', category => 'only', marker => '(intercept)', beta => '-3' }], report => 0);
	is $one->{S}{t}{p}{only}, 1, 'trait_predict: a model of one category gives it a probability of 1';
	is $one->{S}{t}{call}, 'only', 'trait_predict: and calls it';
	# A dosage between 0 and 2, as imputation gives, is a dosage; one the model
	# does not use is not looked at. The probability differs from the formula's
	# by 0 on perl 5.44.0; _softmax takes the larger predictor out first, which
	# can move the last bit, and 1e-12 is t/hirisplex.t's bound for the same.
	my $frac = trait_predict(genotypes => { S => { 'rs12913832_T' => '1.5', other => 'NA' } }, model => $dose_model, report => 0);
	cmp_ok abs($frac->{S}{eye}{p}{blue} - exp(1.5) / (exp(1.5) + 1)), '<', 1e-12, 'trait_predict: a dosage of 1.5 is multiplied in';
	is $frac->{S}{eye}{called}, 1, 'trait_predict: and is a called genotype';
}

# ---- pdb_search
#
# Each of these is refused before RCSB is asked, so _rcsb_post is replaced by
# one that dies, which would show in the message if it were reached.
{
	no warnings 'redefine';
	local *Cn::Bioinf::_rcsb_post = sub { die "the network was used\n" };
	my $seq = 'MVLSPADKTNVKAAWGKV';
	my @bad = (
		[[], qr/^pdb_search: needs "seq"/, 'no sequence'],
		[[seq => ['MK']], qr/^pdb_search: "seq" must be a sequence, as a string, not an ARRAY ref/, 'a sequence that is a ref'],
		[[seq => " \n*"], qr/^pdb_search: "seq" has no residues/, 'a sequence of only whitespace and a stop'],
		[[seq => 'MVLS-PADK'], qr/^pdb_search: "seq" has "-" at residue 5, which isn't a residue letter/, 'a gap in the sequence'],
		[[seq => 'MV*LSPADK'], qr/^pdb_search: "seq" has "\*" at residue 3, which isn't a residue letter/, 'a stop that is not at the end'],
		[[seq => $seq, evalue => 0], qr/^pdb_search: "evalue" must be a number above 0, not "0"/, 'an evalue of 0'],
		[[seq => $seq, evalue => 'x'], qr/^pdb_search: "evalue" must be a number above 0, not "x"/, 'an evalue that is not a number'],
		[[seq => $seq, identity => 101], qr/^pdb_search: "identity" must be a percentage, from 0 to 100, not "101"/, 'an identity above 100'],
		[[seq => $seq, identity => -1], qr/^pdb_search: "identity" must be a percentage, from 0 to 100, not "-1"/, 'a negative identity'],
		[[seq => $seq, max_hits => 0], qr/^pdb_search: "max_hits" must be a whole number from 1 to 10000, RCSB's limit, not "0"/, 'a max_hits of 0'],
		[[seq => $seq, max_hits => 10001], qr/^pdb_search: "max_hits" must be a whole number from 1 to 10000, RCSB's limit, not "10001"/, 'a max_hits past RCSB\'s limit'],
		[[seq => $seq, max_hits => 2.5], qr/^pdb_search: "max_hits" must be a whole number from 1 to 10000, RCSB's limit, not "2\.5"/, 'a fractional max_hits'],
		[[seq => $seq, computed => 'yes'], qr/^pdb_search: "computed" must be 0 or 1, not "yes"/, 'a computed of "yes"'],
		[[seq => $seq, type => 'DNA'], qr/^pdb_search: "type" must be protein, dna or rna, not "DNA"/, 'an unknown type'],
		[[seq => $seq, report => 2], qr/^pdb_search: "report" must be 0 or 1, not "2"/, 'a report of 2'],
		[[seq => $seq, timeout => 0], qr/^pdb_search: "timeout" must be a number of seconds above 0, not "0"/, 'a timeout of 0'],
		[[seq => $seq, fasta => 'x.fa'], qr/^pdb_search: doesn't know "fasta"; it takes computed, evalue, identity, max_hits, report, seq, timeout, type/, 'an unknown option'],
	);
	foreach my $case (@bad) {
		my ($args, $re, $name) = @{ $case };
		dies_like sub { pdb_search(@{ $args }) }, $re, "pdb_search: $name dies";
	}
	# the positive sentinel: a good call does get as far as asking RCSB
	dies_like sub { pdb_search(seq => $seq) }, qr/^the network was used/, 'pdb_search: a good call reaches RCSB';
}

#  the XSUBs that read a VCF, when the reader grows perl's stack
#
# Up to 2026-10-26 each XSUB that reads a VCF through a perl callback pushed
# its results with the stack pointer it took on entry. A callback that grew
# perl's stack moved it, and the results went to freed memory: _vcf_samples
# returned nothing instead of its samples, found by review. Each XSUB is run
# in a child perl, where the stack is still small, once with a reader that
# calls a sub with 2,000,000 arguments (16 MB of stack on a 64-bit perl, which
# glibc's malloc cannot grow in place from its heap, so the stack moves), and
# once with one that does not. The expected value is the second run's: what
# each XSUB counts is tested in t/sex.t, t/parent.t, t/relatedness.t,
# t/concordance.t and t/hirisplex.t, and this checks only that a moved stack
# does not change it.
{
	my $two = spew('stack.two.vcf', join '', map { join("\t", @{ $_ }) . "\n" }
		['##fileformat=VCFv4.2'],
		['#CHROM', qw(POS ID REF ALT QUAL FILTER INFO FORMAT A B)],
		[qw(1 100 rs1 A G . PASS . GT:DP 0/1:30 1/1:30)],
		[qw(1 200 rs2 C T . PASS . GT:DP 0/0:30 0/1:30)],
		[qw(X 5000000 rs3 G A . PASS . GT:DP 1/1:30 0/1:30)]);
	my $one_a = spew('stack.a.vcf', "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n"
		. "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n");
	my $one_b = spew('stack.b.vcf', "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\n"
		. "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t1/1\n");
	# The child's reader hands over 16 bytes at a time, so that lines span
	# blocks; $grown is 1 once its first call has grown the stack, for the
	# sentinel.
	my $child = spew('stack.pl', <<'CHILD');
use strict;
use warnings FATAL => 'all';
use Storable ();
use Cn::Bioinf ();
my ($case, $grow, $out, $two, $one_a, $one_b) = @ARGV;
my $grown = 0;
sub deep { scalar @_ }
sub reader {
	my ($file) = @_;
	open my $fh, '<:raw', $file or die "$file: $!";
	sub {
		if ($grow && !$grown) {
			deep((1) x 2_000_000);
			$grown = 1;
		}
		my $n = read $fh, my $block, 16;
		die "$file: $!" unless defined $n;
		$n ? $block : undef;
	};
}
my %call = (
	vcf_samples      => sub { Cn::Bioinf::_vcf_samples(reader($two), $two) },
	sex_tally        => sub { Cn::Bioinf::_sex_tally(reader($two), $two, 0) },
	parent_one       => sub { Cn::Bioinf::_parent_one(reader($two), $two, 0, 'A', 'B') },
	parent_own       => sub { Cn::Bioinf::_parent_own(reader($one_a), $one_a, reader($one_b), $one_b, 0) },
	relatedness_one  => sub { Cn::Bioinf::_relatedness_one(reader($two), $two, 0) },
	relatedness_own  => sub { Cn::Bioinf::_relatedness_own([reader($one_a), reader($one_b)], [$one_a, $one_b], 0) },
	concordance      => sub { Cn::Bioinf::_concordance(reader($two), $two, reader($two), $two, ['A', 'B'], 0, 1) },
	hirisplex_dosage => sub { Cn::Bioinf::_hirisplex_dosage(reader($two), $two, 0, [[1, 100, 'G', 0, 'rs1', 0], [1, 200, 'T', 0, 'rs2', 1]], 2, 0) },
);
my @result = $call{$case}->();
Storable::nstore({ grown => $grown, result => \@result }, $out);
CHILD
	require Storable;
	# each XSUB and how many values it returns
	my @cases = ([vcf_samples => 1], [sex_tally => 3], [parent_one => 12], [parent_own => 12],
		[relatedness_one => 3], [relatedness_own => 3], [concordance => 2], [hirisplex_dosage => 3]);
	foreach my $case (@cases) {
		my ($name, $values) = @{ $case };
		my %run; # 0 = the stack left alone, 1 = grown by the reader
		foreach my $grow (0, 1) {
			my $out = catfile($dir, "stack.$name.$grow");
			my $status = system($^X, @inc, $child, $name, $grow, $out, $two, $one_a, $one_b);
			is $status, 0, "_$name: the child perl exits 0 " . ($grow ? 'with' : 'without') . ' the stack grown';
			$run{$grow} = $status == 0 ? Storable::retrieve($out) : { grown => -1, result => [] };
		}
		is $run{1}{grown}, 1, "_$name: the reader grew the stack";
		is scalar @{ $run{0}{result} }, $values, "_$name: returns $values values";
		is_deeply $run{1}{result}, $run{0}{result}, "_$name: returns the same when the reader grows the stack, rather than losing its results";
	}
}

done_testing;
