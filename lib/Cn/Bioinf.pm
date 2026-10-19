package Cn::Bioinf;
# ABSTRACT: FASTA I/O in XS, BLAST hit ranking, and alignment plots and tables
#
# Taken from ~/Scripts/bioinf.pm (as it stood on 2026-08-19): fasta2hash,
# get_best_alignment_hit, hash2fasta_file, msa_quality_table, msa_phylo_plot and
# clustal_view_residues. fasta2hash and hash2fasta_file are now XS (Bioinf.xs).
# msa_phylo_plot is now plot_msa and plot_phylo, and no longer takes a BLAST
# report and a proteome database: each is given the sequences, aligns them with
# Clustal Omega and draws the result itself, with share/msa_plot.py (from
# ~/Scripts/HitList/scripts/x.my.align.py) in place of that script and of the
# R/ggtree tree. msa_quality_table's table is drawn by that script too, where
# bioinf.pm had Matplotlib::Simple's colored_table draw it. clustalo, blastp and
# that Python all come from Alien::Bioinf, never from PATH. Alien::Bioinf is only
# recommended, not required, since it cannot be installed everywhere (NCBI
# builds BLAST+ for few platforms) and fasta2hash, hash2fasta_file and
# get_best_alignment_hit do not need it; the functions that do load it through
# _alien, which says what is missing.
require 5.010;
use strict;
use warnings;
use Carp qw(croak carp);
use File::Spec;
use Sys::Hostname 'hostname';
use Exporter 'import';
use XSLoader;
our $VERSION = '0.01';
our @EXPORT_OK = qw(clustal_view_residues fasta2hash get_best_alignment_hit hash2fasta_file msa_quality_table parent pgs_score plot_msa plot_phylo sex_check);
our %EXPORT_TAGS = (all => \@EXPORT_OK);
XSLoader::load('Cn::Bioinf', $VERSION);
# The calling script's full path, for _creator. It is looked up at load time,
# as FindBin always does, but FindBin croaks when $0 is not a file on disk
# (after $0 is assigned to, or under an embedder), and that must not stop this
# module loading, so $0 made absolute stands in then. Up to 2026-10-17 this was
# the script's file name alone, which _creator joined to the working directory,
# so a script run from another directory was recorded at a path that wasn't it.
my $Script = eval { require FindBin; File::Spec->catfile($FindBin::RealBin, $FindBin::RealScript) } // File::Spec->rel2abs($0);

# ---- private helpers

# The name => value pairs in @$list as a hash ref. Dies unless they hold every
# @$required key and nothing outside @$required and @$optional. Its messages
# are prefixed with $sub, the public function whose arguments they are, not
# with "_args".
sub _args {
	my ($list, $sub, $required, $optional) = @_;
	croak "$sub: takes name => value pairs, not a hash ref" if @{ $list } == 1 && ref $list->[0] eq 'HASH';
	croak "$sub: takes name => value pairs, and was given an odd number of arguments" if @{ $list } % 2;
	my $args = { @{ $list } };
	my @missing = grep { !defined $args->{$_} } @{ $required };
	croak "$sub: needs " . join(', ', map { "\"$_\"" } @missing) if @missing;
	my %ok = map { $_ => 1 } @{ $required }, @{ $optional };
	my @bad = sort grep { !$ok{$_} } keys %{ $args };
	croak "$sub: doesn't know " . join(', ', map { "\"$_\"" } @bad) . '; it takes ' . join(', ', sort keys %ok) if @bad;
	$args;
}

# $v as a message names it: "a HASH ref", "an ARRAY ref", "undef", or the
# string quoted.
sub _what {
	my ($v) = @_;
	return 'undef' unless defined $v;
	ref $v ?(ref($v) =~ /\A[AEIOU]/ ? 'an ' : 'a ') . ref($v) . ' ref' : "\"$v\"";
}

# Dies unless each of @keys of $args that was given is a $type ('ARRAY' or
# 'HASH') ref, naming $sub, the public function whose arguments they are. Up to
# 2026-10-17 a wrong one died with perl's own "Not an ARRAY reference", naming
# no function, and often only once clustalo had run.
sub _refs {
	my ($sub, $args, $type, @keys) = @_;
	foreach my $k (@keys) {
		next if !defined $args->{$k} || ref $args->{$k} eq $type;
		croak "$sub: \"$k\" must be " . ($type eq 'ARRAY' ? 'an array' : 'a hash') . ' ref, not ' . _what($args->{$k});
	}
}

# $s as characters. A string perl holds as characters is taken as it is; one it
# holds as bytes is taken to be UTF-8 if it decodes as such, which is how the
# command line used to pass it to Python, and Latin-1 otherwise.
sub _chars {
	my $s = "$_[0]";
	utf8::decode($s) unless utf8::is_utf8($s);
	$s;
}

# $s, read as _chars reads it, as UTF-8 bytes, which is what msa_plot.py reads
# a FASTA or newick file as.
sub _utf8 {
	my $s = _chars($_[0]);
	utf8::encode($s);
	$s;
}

# The start of the "Creator" written into an image's metadata, after the one
# Matplotlib::Simple 0.319 writes: which script, run where, called which sub of
# which version of this module (and of Alien::Bioinf, whose tools and Python
# draw every image), by whom, on which computer (its hostname and OS), with
# which perl and where it is. Like Matplotlib::Simple's, the Creator is the image's
# whole provenance, on one line, and the only metadata written: msa_plot.py
# appends the Python, matplotlib and Biopython versions as it runs, since only
# it knows them, and then each clause it is given with "--p" (what the image
# was made from, with digests, and the commands that made it). Matplotlib::Simple has to py_str() its
# Creator, which it pastes into a Python literal; this one reaches Python as a
# JSON string through _python's --argfile, so backslashes and quotes in a path
# need no escaping here.
sub _creator {
	my ($sub) = @_;
	_alien($sub);
	"$Script called using \"$sub\" in " . __FILE__ . " version $VERSION with Alien::Bioinf $Alien::Bioinf::VERSION"
		. ' by user ' . _user() . ' on host ' . hostname() . " ($^O) with Perl " . sprintf('%vd', $^V) . " ($^X)";
}

# The name of the user running this. getpwuid is unimplemented on MSWin32,
# where it dies, and getlogin is undef without a controlling terminal (under
# cron, or a batch scheduler), so the environment is the last resort.
sub _user {
	my $user = eval { (getpwuid $<)[0] };
	$user = getlogin() unless defined $user && length $user;
	$user = $ENV{USERNAME} // $ENV{USER} unless defined $user && length $user;
	defined $user && length $user ? $user : 'unknown';
}

sub _json_file {
	my ($file) = @_;
	require JSON::MaybeXS;
	open my $fh, '<:raw', $file or croak "_json_file: can't read $file: $!";
	local $/;
	JSON::MaybeXS::decode_json(<$fh>);
}

# @cmd run with list-form system, which needs no shell. $? is -1 when the
# program never started; otherwise its low 7 bits are the signal that killed
# it, and must be read before the shift, which leaves only the exit code.
sub _run {
	my @cmd = @_;
	system(@cmd) == 0 and return;
	my $why = $? == -1 ? "it didn't start: $!" : $? & 127 ? 'killed by signal ' . ($? & 127) : 'exit ' . ($? >> 8);
	croak "_run: \"@cmd\" failed: $why";
}

# Alien::Bioinf, loaded for $sub, or a message saying that $sub needs it.
sub _alien {
	my ($sub) = @_;
	eval { require Alien::Bioinf; 1 }
		or croak "_alien: $sub needs Alien::Bioinf, for Clustal Omega, BLAST+ and the plotting Python, and it can't be loaded: $@";
	'Alien::Bioinf';
}

# share/msa_plot.py run with @args under Alien::Bioinf's Python, for $sub. The
# arguments go in a JSON file rather than on the command line: perl's
# system(LIST) on Windows does not escape a '"' inside an argument, and the
# Creator and the JSON options are full of them. Each argument is read as
# _chars reads it; encode_json alone would read UTF-8 bytes as Latin-1. -I
# keeps PYTHONPATH, PYTHONHOME and the user's own site-packages out, so that
# nothing but Alien::Bioinf's venv is imported.
sub _python {
	my ($sub, @args) = @_;
	require JSON::MaybeXS;
	my @chars = map { _chars($_) } @args;
	my $argfile = _tmp('.json');
	open my $fh, '>:raw', $argfile or croak "_python: can't write $argfile: $!";
	print {$fh} JSON::MaybeXS::encode_json(\@chars);
	close $fh or croak "_python: can't write $argfile: $!";
	_run(_alien($sub)->python, '-I', _script(), '--argfile', $argfile);
}

# $data, whose strings the caller has made characters with _chars, as JSON
# text: characters, not UTF-8 bytes, since it is one argument that _python
# encodes along with the rest. It is held as UTF-8, so that _python takes it as
# it is: decoding it whole, as up to 2026-10-17, turned a byte string inside it
# to mojibake whenever another string in it held characters.
sub _json_text {
	my ($data) = @_;
	require JSON::MaybeXS;
	my $text = JSON::MaybeXS->new(utf8 => 0, canonical => 1)->encode($data);
	utf8::upgrade($text);
	$text;
}

# "wrote $file" to STDOUT, as LikeR's write_table prints it: $file in black
# (30) on $bg (an ANSI background, 43 = yellow or 46 = cyan) and then reset
# (0), written out rather than loaded from Term::ANSIColor; the reset comes
# before the newline so the colour doesn't run on into the next line. The
# colour is left out when STDOUT is not a terminal, where the escapes would
# only be noise in a log or a test harness.
sub _wrote {
	my ($file, $bg) = @_;
	print STDOUT -t STDOUT ? "wrote \e[30;${bg}m$file\e[0m\n" : "wrote $file\n";
}

sub _tmp {
	my ($suffix) = @_;
	require File::Temp;
	my ($fh, $name) = File::Temp::tempfile(SUFFIX => $suffix, UNLINK => 1, TMPDIR => 1);
	close $fh;
	$name;
}

sub _script {
	require File::ShareDir;
	File::ShareDir::dist_file('Cn-Bioinf', 'msa_plot.py');
}

# The FASTA file's records and their deflines in file order.
sub _fasta_ordered {
	my ($file) = @_;
	my @order;
	my $fh = _open_fasta($file);
	my $h = _read_fasta($fh, undef, $file, \@order);
	close $fh;
	($h, \@order);
}

sub _open_fasta {
	my ($file) = @_;
	croak "_open_fasta: undef isn't a readable file" unless defined $file;
	croak "_open_fasta: \"$file\" doesn't exist or isn't a readable file" unless -f $file && -r _;
	my $fh;
	if ($file =~ /\.gz\z/) {
		# in an eval, since a perl without list-form piped open (MSWin32 before
		# 5.22) dies in open itself, with no function name in the message
		eval { open $fh, '-|', 'gzip', '-dc', $file } or croak "_open_fasta: can't run gzip on $file: " . ($@ || $!);
	} else {
		open $fh, '<:raw', $file or croak "_open_fasta: can't read $file: $!";
	}
	$fh;
}

sub _first_letter {
	my ($string) = @_;
	return $1 if $string =~ /^([A-Za-z])/;
	croak "_first_letter: \"$string\" doesn't start with a letter, so it has no first letter";
}

# Narrow but unambiguous labels, e.g. "$\it{C.al.}$", for [genus, species,
# strain] triples ('' for no strain). One letter of the species cannot tell
# C.albicans from C.auris, so this takes the shortest species prefix that
# separates every entry; the strain is always kept, since it is all that
# separates the three C.neoformans.
sub _abbreviated_labels {
	my ($names) = @_;
	my $max = 0;
	foreach my $n (@{ $names }) {
		$max = length $n->[1] if length $n->[1] > $max;
	}
	foreach my $n (1 .. $max) {
		my @labels = map {
			'$\it{' . ucfirst(_first_letter($_->[0])) . '.' . lc(substr($_->[1], 0, $n)) . '.}$' . ($_->[2] ne '' ? " $_->[2]" : '')
		} @{ $names };
		my %seen;
		return \@labels unless grep { $seen{$_}++ } @labels;
	}
	croak "_abbreviated_labels: can't tell " . join(', ', map { join '.', @{ $_ } } @{ $names }) . ' apart, even with the full species name and strain';
}

# ---- FASTA

sub fasta2hash {
	my ($file, $key) = @_;
	my $fh = _open_fasta($file);
	my $r = _read_fasta($fh, $key, $file);
	# A gzip that fails part-way, on a truncated or corrupt .gz, has merely
	# stopped writing, so its exit status is the only sign of it. That is
	# checked unless $key was found, when reading may have stopped early and
	# the closed pipe have killed gzip with a SIGPIPE.
	unless (close $fh) {
		my $why = $! ? "$!" : $? & 127 ? 'gzip was killed by signal ' . ($? & 127) : 'gzip exited ' . ($? >> 8);
		croak "fasta2hash: couldn't read all of $file: $why" unless defined $key && defined $r;
	}
	if (defined $key) {
		return $r if defined $r;
		croak "fasta2hash: couldn't find \"$key\" in $file";
	}
	croak "fasta2hash: couldn't find any sequences in $file" unless %{ $r };
	$r;
}

sub hash2fasta_file {
	my ($hash, $filename, $order, $width) = @_;
	croak "hash2fasta_file: 1st argument must be a hash ref" unless ref $hash eq 'HASH';
	croak "hash2fasta_file: 2nd argument must be a file name" unless defined $filename && ref $filename eq '';
	croak "hash2fasta_file: 3rd argument must be an array ref or undef" if defined $order && ref $order ne 'ARRAY';
	# checked here because the XS takes it as a UV, where -1 would wrap round
	croak "hash2fasta_file: 4th argument must be a whole number, 0 or more, or undef, not \"$width\"" if defined $width && $width !~ /\A[0-9]+\z/;
	open my $fh, '>:raw', $filename or croak "hash2fasta_file: can't write $filename: $!";
	_write_fasta($fh, $hash, $order // [sort keys %{ $hash }], $width // 80);
	close $fh or croak "hash2fasta_file: can't write $filename: $!";
	$filename;
}

# ---- BLAST

sub get_best_alignment_hit {
	my ($json_file, $sort_criterion) = @_;
	$sort_criterion //= 'evalue';
	# For each query in one BLAST JSON file (-outfmt 15), that query's hits as
	# an array sorted best-first on $sort_criterion, one element per hit: the
	# hit's single best hsp plus the fields that say which subject it is.
	#
	# The hits are an array, not a hash on the hit title, because two hits can
	# share a title (one protein deposited under two accessions): in
	# identify.target/potato/new/f.sambunicum.json 76 hits carry 44 distinct
	# titles, and a hash kept 44.
	croak "get_best_alignment_hit: undef isn't a readable file" unless defined $json_file;
	croak "get_best_alignment_hit: \"$json_file\" doesn't exist or isn't a readable file" unless -f $json_file && -r _;
	# Which way is better for each numeric hsp field: -1 = smaller, 1 = larger.
	# The other hsp keys (hseq, qseq, midline, num, *_from/*_to) are content
	# or position, not quality.
	my %better = (align_len => 1, bit_score => 1, evalue => -1, gaps => -1, identity => 1, positive => 1, score => 1);
	my $direction = $better{$sort_criterion}
		// croak "get_best_alignment_hit: \"$sort_criterion\" isn't a sortable hsp field; use one of " . join(', ', sort keys %better);
	my $blast = _json_file($json_file);
	my %alignments;
	foreach my $report (@{ $blast->{BlastOutput2} }) {
		my @hits;
		foreach my $hit (@{ $report->{report}{results}{search}{hits} }) {
			my $best; # undef until an hsp is seen; a hit with no hsps has nothing to record
			foreach my $hsp (@{ $hit->{hsps} }) {
				croak "get_best_alignment_hit: an hsp of \"$hit->{description}[0]{title}\" in $json_file has no \"$sort_criterion\" to sort on"
					unless defined $hsp->{$sort_criterion};
				unless (defined $best) {
					$best = $hsp;
					next;
				}
				my $cmp = $direction * ($hsp->{$sort_criterion} <=> $best->{$sort_criterion});
				# BLAST rounds any strong e-value to 0, so ties are the rule; the bit score breaks them
				$cmp = $hsp->{bit_score} <=> $best->{bit_score} if $cmp == 0 && $sort_criterion ne 'bit_score';
				$best = $hsp if $cmp > 0;
			}
			next unless defined $best;
			# A shallow copy (an hsp holds only numbers and strings), so these
			# fields are not written into $blast. The hit's own "num" is not
			# merged, as it would shadow the hsp's; its "len" is the subject's
			# full length, so it arrives as "hit_len" beside the hsp's "align_len".
			# A hit has one description for each identical sequence in a
			# non-redundant database such as nr; it is named by the first, as
			# BLAST lists it, and every one is in "descriptions". Up to
			# 2026-10-17 a hit with more than one died.
			push @hits, { %{ $best },
				accession => $hit->{description}[0]{accession},
				descriptions => [map { +{ %{ $_ } } } @{ $hit->{description} }],
				hit_len   => $hit->{len},
				id        => $hit->{description}[0]{id},
				title     => $hit->{description}[0]{title},
			};
		}
		# Indices are sorted so the last tiebreak can be BLAST's own order; the
		# array runs best-first, hence $b before $a.
		my @order = sort {
			   ($direction * ($hits[$b]{$sort_criterion} <=> $hits[$a]{$sort_criterion}))
			|| ($hits[$b]{bit_score} <=> $hits[$a]{bit_score})
			|| ($a <=> $b)
		} 0 .. $#hits;
		# Two queries of one title would leave only the second's hits, so that
		# is refused rather than lost silently.
		my $query = $report->{report}{results}{search}{query_title};
		croak "get_best_alignment_hit: two queries in $json_file are both titled \"$query\"" if exists $alignments{$query};
		$alignments{$query} = [ @hits[@order] ];
	}
	\%alignments;
}

# ---- polygenic scores

# The rows of $frame, given to $sub as $what, as an array ref of hash refs, in
# whichever of Stats::LikeR's four shapes it comes: an array of hashes
# (read_table's default), a hash of hashes (its default for a VCF), a hash of
# arrays, or an array of arrays whose first row is the header. A row of an AoH or
# HoH is the frame's own hash, and is only read; one of a HoA or AoA is made
# here. With $labelled, also returns how a message names each row: "row 3",
# counting data rows from 1, or 'row "key"' for a hash of hashes, whose rows are
# then in key order so that the first bad one is always the same one.
sub _frame_rows {
	my ($sub, $what, $frame, $labelled) = @_;
	my $shapes = 'an array of hashes or of arrays, or a hash of hashes or of arrays, as Stats::LikeR\'s read_table returns';
	my (@rows, @labels);
	if (ref $frame eq 'ARRAY') {
		croak "$sub: \"$what\" has no rows" unless @{ $frame };
		if (ref $frame->[0] eq 'HASH') {
			@rows = @{ $frame };
		} elsif (ref $frame->[0] eq 'ARRAY') {
			my ($head, @data) = @{ $frame };
			croak "$sub: \"$what\" has a header and no rows" unless @data;
			foreach my $i (0 .. $#data) {
				croak "$sub: row " . ($i + 1) . " of \"$what\" must be an array ref, as its header is, not " . _what($data[$i])
					unless ref $data[$i] eq 'ARRAY';
				my %row;
				@row{ @{ $head } } = @{ $data[$i] };
				push @rows, \%row;
			}
		} else {
			croak "$sub: \"$what\" must be $shapes, not an array whose first row is " . _what($frame->[0]);
		}
		@labels = map { "row $_" } 1 .. @rows if $labelled;
	} elsif (ref $frame eq 'HASH') {
		my ($key, $first) = each %{ $frame };
		keys %{ $frame }; # resets the iterator that each left part-way
		croak "$sub: \"$what\" has no rows" unless defined $key;
		if (ref $first eq 'HASH') {
			if ($labelled) {
				my @keys = sort keys %{ $frame };
				@rows = @{ $frame }{@keys};
				@labels = map { "row \"$_\"" } @keys;
			} else {
				@rows = values %{ $frame };
			}
		} elsif (ref $first eq 'ARRAY') {
			my @cols = keys %{ $frame };
			foreach my $c (@cols) {
				croak "$sub: column \"$c\" of \"$what\" must be an array ref, as \"$key\" is, not " . _what($frame->{$c})
					unless ref $frame->{$c} eq 'ARRAY';
			}
			foreach my $i (0 .. $#{ $first }) {
				push @rows, { map { $_ => $frame->{$_}[$i] } @cols };
			}
			croak "$sub: \"$what\" has no rows" unless @rows;
			@labels = map { "row $_" } 1 .. @rows if $labelled;
		} else {
			croak "$sub: \"$what\" must be $shapes, not a hash whose \"$key\" is " . _what($first);
		}
	} else {
		croak "$sub: \"$what\" must be $shapes, not " . _what($frame);
	}
	foreach my $i (0 .. $#rows) {
		croak "$sub: row " . ($i + 1) . " of \"$what\" must be a hash ref, not " . _what($rows[$i]) unless ref $rows[$i] eq 'HASH';
	}
	$labelled ? (\@rows, \@labels) : \@rows;
}

# A chromosome's name as a VCF, PLINK 2 and a scoring file are compared on: in
# upper case, without UCSC's "chr" prefix (GRCh37's own names, and the PGS
# Catalog's, have none), with PLINK's numeric codes for the sex and
# mitochondrial chromosomes and UCSC's "M" named as Ensembl names them, and with
# X's pseudo-autosomal regions, which PLINK 2 names PAR1 and PAR2 (or XY, its
# 25) once --split-par has split them off, taken as X.
my %Chrom_alias = (23 => 'X', 24 => 'Y', 25 => 'X', 26 => 'MT', M => 'MT', PAR1 => 'X', PAR2 => 'X', XY => 'X');
sub _chrom {
	my $c = uc $_[0];
	$c =~ s/\ACHR//;
	$Chrom_alias{$c} // $c;
}

# A decimal number, as an effect_weight must be. Infinity and NaN, which
# looks_like_number would take, are left out, as either would make every score
# the same non-number.
my $Number = qr/\A[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?\z/;

# The rows of "weights" one at a time, each a hash ref with how a message names
# it, and the column names: from a data frame (see _frame_rows), or read from a
# scoring file a line at a time, so that a score of millions of variants is
# never held whole. An empty field of a file is '', where read_table gives
# undef; every check below treats the two alike.
sub _weight_rows {
	my ($sub, $weights) = @_;
	if (ref $weights) {
		my ($rows, $labels) = _frame_rows($sub, 'weights', $weights, 1);
		my $i = 0;
		return (sub { $i < @{ $rows } ? ($rows->[$i], $labels->[$i++]) : () }, [keys %{ $rows->[0] }]);
	}
	croak "$sub: \"weights\" file \"$weights\" doesn't exist or isn't readable" unless -f $weights && -r _;
	my $fh;
	if ($weights =~ /\.gz\z/i) {
		require IO::Uncompress::Gunzip;
		$fh = IO::Uncompress::Gunzip->new($weights, MultiStream => 1)
			or croak "$sub: can't read $weights: $IO::Uncompress::Gunzip::GunzipError";
	} else {
		open $fh, '<:raw', $weights or croak "$sub: can't read $weights: $!";
	}
	my $head;
	while (defined(my $line = <$fh>)) {
		$line =~ s/\r?\n\z//;
		next if $line eq '' || $line =~ /\A#/;
		$head = $line;
		last;
	}
	croak "$sub: \"weights\" file $weights has no header line" unless defined $head;
	my @cols = split /\t/, $head, -1;
	my $n = 0;
	my $next = sub {
		while (defined(my $line = <$fh>)) {
			$line =~ s/\r?\n\z//;
			next if $line eq '';
			my %row;
			@row{@cols} = split /\t/, $line, -1;
			return (\%row, 'row ' . ++$n);
		}
		return;
	};
	($next, \@cols);
}

# The data frame $frame written to $file as a VCF of GT alone, for PLINK 2 to
# read. Its rows are sorted by chromosome and position, as PLINK 2 needs them
# to be, under the names _chrom gives the chromosomes, so that "chr1" and "1"
# are one chromosome. A row with no CHROM, POS, REF or ALT is left out, and so
# matches nothing, and a missing or empty GT is written "./.". Returns how many
# rows were written.
sub _frame_vcf {
	my ($sub, $frame, $file, $match) = @_;
	my $rows = _frame_rows($sub, 'vcf', $frame);
	# read_table names each sample's genotype column "<sample>.GT" when it splits
	# the sample columns by FORMAT, as it does by default for a VCF
	my @samples = sort map { /\A(.+)\.GT\z/s ? $1 : () } keys %{ $rows->[0] };
	croak "$sub: \"vcf\" has no <sample>.GT column; read the VCF with read_table's default of explode => 1" unless @samples;
	my @absent = grep { !exists $rows->[0]{$_} } qw(CHROM POS REF ALT), $match eq 'id' ? 'ID' : ();
	croak "$sub: \"vcf\" has no " . join(', ', map { "\"$_\"" } @absent) . ' column' . (@absent > 1 ? 's' : '') if @absent;
	my @keep = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] }
		map { [_chrom($_->{CHROM}), $_->{POS}, $_] }
		grep { defined $_->{CHROM} && length $_->{CHROM} && defined $_->{POS} && $_->{POS} =~ /\A[1-9][0-9]*\z/
			&& defined $_->{REF} && length $_->{REF} && defined $_->{ALT} && length $_->{ALT} } @{ $rows };
	open my $fh, '>:raw', $file or croak "_frame_vcf: can't write $file: $!";
	print {$fh} "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
		join("\t", '#CHROM', qw(POS ID REF ALT QUAL FILTER INFO FORMAT), @samples), "\n";
	foreach my $k (@keep) {
		my $r = $k->[2];
		my $id = defined $r->{ID} && length $r->{ID} ? $r->{ID} : '.';
		print {$fh} join("\t", $k->[0], $k->[1], $id, $r->{REF}, $r->{ALT}, '.', '.', '.', 'GT',
			map { my $gt = $r->{"$_.GT"}; defined $gt && length $gt ? $gt : './.' } @samples), "\n";
	}
	close $fh or croak "_frame_vcf: can't write $file: $!";
	scalar @keep;
}

# plink2 run with @args for $sub, its output in $dir/plink2.out rather than on
# the caller's terminal; a failure dies with the lines of it that say why.
sub _plink2 {
	my ($sub, $plink2, $dir, @args) = @_;
	my $log = File::Spec->catfile($dir, 'plink2.out');
	open my $out, '>&', \*STDOUT or croak "_plink2: can't save STDOUT: $!";
	open my $err, '>&', \*STDERR or croak "_plink2: can't save STDERR: $!";
	open STDOUT, '>', $log or croak "_plink2: can't write $log: $!";
	open STDERR, '>&', \*STDOUT or croak "_plink2: can't send STDERR to $log: $!";
	my $rc = system($plink2, @args);
	my $status = $?;
	open STDOUT, '>&', $out or croak "_plink2: can't restore STDOUT: $!";
	open STDERR, '>&', $err or croak "_plink2: can't restore STDERR: $!";
	return if $rc == 0;
	my $why = $status == -1 ? "it didn't start: $!" : $status & 127 ? 'killed by signal ' . ($status & 127) : 'exit ' . ($status >> 8);
	# what PLINK 2 says from its "Error:" on, or else its last line; it
	# rewrites its progress lines with a CR
	my @said;
	if (open my $fh, '<', $log) {
		@said = grep { /\S/ } map { (my $l = $_) =~ s/\s+\z//; $l } map { split /\r/ } <$fh>;
	}
	my ($first) = grep { $said[$_] =~ /\AError/ } 0 .. $#said;
	my $said = join ' ', defined $first ? @said[$first .. $#said] : @said ? $said[-1] : 'nothing';
	croak "$sub: \"$plink2 @args\" failed ($why): $said";
}

sub pgs_score {
	my $sub = 'pgs_score';
	my $args = _args(\@_, $sub, [qw(vcf weights)], [qw(build match memory threads)]);
	my %r;
	$r{match} = $args->{match} // 'position'; # 'position' = CHROM and POS against chr_name and chr_position; 'id' = ID against rsID
	croak "$sub: \"match\" must be \"position\" or \"id\", not " . _what($r{match})
		unless !ref $r{match} && ($r{match} eq 'position' || $r{match} eq 'id');
	$r{build} = $args->{build}; # undef = the sex chromosomes are left out; hg19, b37, hg38 or b38 = PLINK 2 reads them, with X's pseudo-autosomal regions placed for that build
	croak "$sub: \"build\" must be hg19, b37, hg38 or b38, not " . _what($r{build})
		if defined $r{build} && (ref $r{build} || $r{build} !~ /\A(?:hg19|b37|hg38|b38)\z/);
	foreach my $k (qw(memory threads)) {
		# undef = PLINK 2's own default: half the machine's memory, every CPU
		croak "$sub: \"$k\" must be a whole number, 1 or more, not " . _what($args->{$k})
			if defined $args->{$k} && (ref $args->{$k} || $args->{$k} !~ /\A[1-9][0-9]*\z/);
		$r{$k} = $args->{$k};
	}
	# PLINK 2 v2.0.0-a.7.13LM refuses less: "Invalid --memory argument '256' (minimum 640)"
	croak "$sub: \"memory\" must be 640 or more (MiB), the least PLINK 2 takes, not $r{memory}" if defined $r{memory} && $r{memory} < 640;
	my $vcf = $args->{vcf};
	my @input; # plink2's arguments that read the genotypes; filled below
	if (!ref $vcf) {
		croak "$sub: \"vcf\" file \"$vcf\" doesn't exist or isn't readable" unless -f $vcf && -r _;
		if ($vcf =~ /\A(.*)\.pgen\z/is) {
			my $prefix = $1;
			croak "$sub: \"vcf\" $vcf has no $prefix.pvar (or .pvar.zst) beside it" unless grep { -f "$prefix.$_" } qw(pvar pvar.zst);
			croak "$sub: \"vcf\" $vcf has no $prefix.psam beside it" unless -f "$prefix.psam";
			@input = ('--pfile', $prefix);
		} elsif ($vcf =~ /\.bcf\z/i) {
			@input = ('--bcf', $vcf);
		} elsif ($vcf =~ /\.vcf(?:\.b?gz)?\z/i) {
			@input = ('--vcf', $vcf);
		} else {
			croak "$sub: \"vcf\" file $vcf must end in .vcf, .vcf.gz, .vcf.bgz, .bcf or .pgen, as PLINK 2 reads it by its name";
		}
	} elsif (ref $vcf ne 'ARRAY' && ref $vcf ne 'HASH') {
		croak "$sub: \"vcf\" must be a VCF, BCF or PLINK 2 .pgen file, or a data frame, not " . _what($vcf);
	}

	my ($next, $cols) = _weight_rows($sub, $args->{weights});
	my %have = map { $_ => 1 } @{ $cols };
	# A harmonized (hmPOS) scoring file has its variants on the build its name
	# gives in the hm_ columns, and on the author's build in chr_name,
	# chr_position and rsID, so the hm_ columns are the ones read when present.
	my @key_cols = $r{match} eq 'id' ? ($have{hm_rsID} ? 'hm_rsID' : 'rsID')
		: $have{hm_pos} ? qw(hm_chr hm_pos) : qw(chr_name chr_position);
	my $chr_col = $have{hm_chr} ? 'hm_chr' : $have{chr_name} ? 'chr_name' : undef;
	my @absent = grep { !$have{$_} } @key_cols, 'effect_allele', 'effect_weight';
	croak "$sub: \"weights\" has no " . join(', ', map { "\"$_\"" } @absent) . ' column' . (@absent > 1 ? 's' : '')
		. ', which a PGS Catalog scoring file has' if @absent;
	# Only an additive weight is a weight times the dosage. These columns of the
	# PGS Catalog's format 2.0 mark a variant that is scored some other way, which
	# this function would get wrong, so a variant marked so is refused.
	my @flags = grep { $have{$_} } qw(is_haplotype is_diplotype is_interaction is_dominant is_recessive);
	my @dosage_weights = grep { $have{$_} } qw(dosage_0_weight dosage_1_weight dosage_2_weight);
	# removed as pgs_score returns, or, after an error, as the program ends;
	# File::Temp's newdir, which would do both, is newer than 5.10.0's File::Temp
	require File::Temp;
	my $dir = File::Temp::tempdir(CLEANUP => 1);
	my $regions = File::Spec->catfile($dir, 'regions.txt');
	open my $rfh, '>', $regions or croak "$sub: can't write $regions: $!";
	# Each row of weights is kept as "effect\tother\tweight", the alleles in upper
	# case and the weight as written, which a few million rows can afford.
	my @variant;
	my %at;     # key (CHROM:POS or rsID) => the rows of weights at it, joined by commas
	my %chroms; # chromosome, as _chrom names it => 1, for each of the rows' chromosomes
	while (my ($w, $label) = $next->()) {
		foreach my $c (@flags) {
			croak "$sub: $label of \"weights\" has $c \"$w->{$c}\", and only additive weights can be scored"
				if defined $w->{$c} && $w->{$c} =~ /\A(?:true|1)\z/i;
		}
		foreach my $c (@dosage_weights) {
			croak "$sub: $label of \"weights\" has a $c, and only additive weights can be scored"
				if defined $w->{$c} && length $w->{$c};
		}
		my ($e, $o, $weight) = @{ $w }{qw(effect_allele other_allele effect_weight)};
		croak "$sub: $label of \"weights\" has no effect_allele" unless defined $e && length $e;
		croak "$sub: $label of \"weights\" has effect_weight " . _what($weight) . ', which isn\'t a number'
			unless defined $weight && !ref $weight && $weight =~ $Number;
		my $i = @variant;
		push @variant, join "\t", uc $e, defined $o && length $o ? uc $o : '', $weight;
		my $chr = defined $chr_col && defined $w->{$chr_col} && length $w->{$chr_col} ? _chrom($w->{$chr_col}) : undef;
		my $key;
		if ($r{match} eq 'id') {
			$key = $w->{ $key_cols[0] };
		} else {
			my $pos = $w->{ $key_cols[1] };
			# A harmonized file leaves hm_chr and hm_pos empty for a variant it
			# could not place, and that variant is then not found, rather than an error.
			if (defined $chr && defined $pos && length $pos) {
				croak "$sub: $label of \"weights\" has $key_cols[1] \"$pos\", which isn't a position"
					unless $pos =~ /\A[1-9][0-9]*\z/;
				$key = "$chr:$pos";
				# X's positions are given under PAR1 and PAR2 too, which --split-par
				# may have moved them to
				print {$rfh} map { "$_\t$pos\t$pos\n" } $chr, $chr eq 'X' ? qw(PAR1 PAR2) : ();
			}
		}
		$chroms{$chr} = 1 if defined $chr;
		$at{$key} = defined $at{$key} ? "$at{$key},$i" : $i if defined $key && length $key;
	}
	close $rfh or croak "$sub: can't write $regions: $!";
	croak "$sub: \"weights\" has no rows" unless @variant;
	my @sex = grep { $chroms{$_} } qw(X Y);
	croak "$sub: \"weights\" has variants on " . join(' and ', @sex) . ', which PLINK 2 reads only when "build" says where X\'s'
		. ' pseudo-autosomal regions are: hg19 (b37) or hg38 (b38)' if @sex && !defined $r{build};
	_frame_vcf($sub, $vcf, File::Spec->catfile($dir, 'frame.vcf'), $r{match}) if ref $vcf;
	@input = ('--vcf', File::Spec->catfile($dir, 'frame.vcf')) if ref $vcf;

	my $plink2 = eval { _alien($sub)->plink2 };
	unless (defined $plink2) {
		(my $why = $@) =~ s/ at \S+ line \d+\.?\n?\z//;
		croak "$sub: needs PLINK 2, which Alien::Bioinf installs, and hasn't here: $why";
	}
	my @common = ('--allow-extra-chr', (defined $r{memory} ? ('--memory', $r{memory}) : ()), (defined $r{threads} ? ('--threads', $r{threads}) : ()));
	push @input, '--vcf-half-call', 'missing' unless $input[0] eq '--pfile';
	# The samples, each of unknown sex, written by PLINK 2 without reading the
	# sex chromosomes, which it would refuse to without each one's sex.
	my $samples = File::Spec->catfile($dir, 'samples');
	my @psam = ('--not-chr', 'X,Y,XY', '--make-just-psam', '--out', $samples);
	unless (eval { _plink2($sub, $plink2, $dir, @input, @common, @psam); 1 }) {
		# PLINK 2 refuses a single-sample VCF with no 0/0 among its first 1,000
		# called genotypes, as one that lists only the sites where the sample
		# carries an ALT allele; scored, every other variant is missing.
		croak $@ unless $@ =~ /contain at least one ALT allele/;
		carp "$sub: \"vcf\" lists only the sites where its sample carries an ALT allele, so every variant where it matches"
			. ' the reference is missing, and its score is not comparable with one from a VCF of every site';
		push @input, '--vcf-allow-no-nonvar';
		_plink2($sub, $plink2, $dir, @input, @common, @psam);
	}
	my @import = (@input, @common);
	# Only the chromosomes the score uses are read, or, where it gives none
	# (matching by ID alone), all but the sex chromosomes, unless "build" says
	# where X's pseudo-autosomal regions are. PLINK 2 needs that, and each
	# sample's sex, before it reads X at all; the sex is given as unknown, which
	# counts a haploid call on X as homozygous.
	if (%chroms) {
		push @import, '--chr', join ',', sort(keys %chroms), $chroms{X} ? qw(PAR1 PAR2) : ();
	} elsif (!defined $r{build}) {
		push @import, '--not-chr', 'X,Y,XY';
	}
	push @import, '--psam', "$samples.psam", '--split-par', $r{build} if defined $r{build} && $input[0] ne '--pfile';
	push @import, '--extract', 'bed1', $regions if $r{match} eq 'position';
	my $g = File::Spec->catfile($dir, 'g');
	# PLINK 2 exits 13 when none of the score's positions is in the VCF, which
	# is a result here: every variant is missing.
	my $none = eval { _plink2($sub, $plink2, $dir, @import, '--make-pgen', '--out', $g); 0 }
		// ($@ =~ /\(exit 13\): Error: No variants remaining/ ? 1 : croak $@);

	my @samples; # in PLINK 2's order
	open my $psam, '<', "$samples.psam" or croak "$sub: can't read $samples.psam: $!";
	my $iid; # the column of the IID
	while (my $line = <$psam>) {
		chomp $line;
		my @f = split /\t/, $line;
		if ($line =~ /\A#/) {
			($iid) = grep { $f[$_] =~ /\A#?IID\z/ } 0 .. $#f;
			next;
		}
		push @samples, $f[$iid // 0];
	}
	close $psam;

	# The .pvar a record at a time, matched against the rows of weights, and
	# written again with each record's ID made v1, v2, ..., which --score finds
	# it by, since a VCF's own IDs are often "." and need not be unique.
	my $found = ''; # bit $i set: row $i of weights is at a position (or ID) of the VCF
	my %use;        # row $i => "record number\tA1 as the record writes it\tits CHROM:POS:REF:ALT\t1, or 0 when the record lacks the effect allele", for the record whose alleles agree with it
	my ($in, $out);
	unless ($none) {
		open $in, '<', "$g.pvar" or croak "$sub: can't read $g.pvar: $!";
		open $out, '>', "$g.ids.pvar" or croak "$sub: can't write $g.ids.pvar: $!";
	}
	my ($n, %col) = (0);
	while (!$none && defined(my $line = <$in>)) {
		if ($line =~ /\A##/) {
			print {$out} $line;
			next;
		}
		chomp $line;
		my @f = split /\t/, $line, -1;
		if ($line =~ /\A#CHROM/) {
			%col = map { $f[$_] => $_ } 0 .. $#f;
			print {$out} "$line\n";
			next;
		}
		$n++;
		my ($chrom, $pos, $id, $ref, $alt) = @f[ @col{'#CHROM', 'POS', 'ID', 'REF', 'ALT'} ];
		$f[ $col{ID} ] = "v$n";
		print {$out} join("\t", @f), "\n";
		my @keys = $r{match} eq 'id' ? grep { length && $_ ne '.' } split(/;/, $id) : (_chrom($chrom) . ":$pos");
		my @here = map { split /,/, $at{$_} // '' } @keys;
		next unless @here;
		# An ALT of "." is a site with no alternate allele, which a VCF of every
		# site writes where the sample matches the reference.
		my @alleles = ($ref, $alt eq '.' ? () : split /,/, $alt);
		my %named; # allele in upper case => as the record writes it
		$named{ uc $_ } //= $_ foreach @alleles;
		my $sort = join ':', $chrom, $pos, $ref, $alt;
		foreach my $i (@here) {
			vec($found, $i, 1) = 1;
			my ($e, $o) = split /\t/, $variant[$i];
			# The alleles agree when both of the variant's are the record's, or,
			# for a record with no ALT, when either is its REF.
			next unless @alleles == 1
				? $e eq uc $ref || (length $o && $o eq uc $ref)
				: defined $named{$e} && (!length $o || defined $named{$o});
			# Of two records that agree, as a multi-allelic site split in two
			# can give, the first in CHROM:POS:REF:ALT order is used, so that the
			# order of the file does not decide it.
			next if defined $use{$i} && (split /\t/, $use{$i})[2] le $sort;
			# An effect allele the record lacks (its REF being the other allele)
			# has no copies; it is given as the REF with a weight of 0, which
			# counts the sample's alleles there and adds nothing.
			$use{$i} = join "\t", $n, $named{$e} // $ref, $sort, defined $named{$e} ? 1 : 0;
		}
	}
	unless ($none) {
		close $in;
		close $out or croak "$sub: can't write $g.ids.pvar: $!";
	}

	my %score = map { $_ => { alleles => 0, mismatched => 0, missing => 0, score => 0, variants => 0 } } @samples;
	my $agree = keys %use;
	my $disagree = unpack('%32b*', $found) - $agree;
	# --score takes one line for each record, so a record that two rows of
	# weights fall on is scored in a second pass, and so on.
	# Pass $p's score file and list of records are written as the rows are
	# read, rather than gathered first, so that a score of millions of variants
	# is not held twice.
	my %layer; # record => how many rows have taken it so far
	my @pass;  # pass => [its file's prefix, score file handle, ID list handle]
	foreach my $i (sort { $a <=> $b } keys %use) {
		my ($rec, $a1, undef, $has) = split /\t/, $use{$i};
		my $p = $layer{$rec}++;
		unless ($pass[$p]) {
			my $file = File::Spec->catfile($dir, "score$p");
			open my $fh, '>', "$file.txt" or croak "$sub: can't write $file.txt: $!";
			open my $ids, '>', "$file.ids" or croak "$sub: can't write $file.ids: $!";
			print {$fh} "ID\tA1\tWEIGHT\n";
			$pass[$p] = [$file, $fh, $ids];
		}
		print { $pass[$p][1] } "v$rec\t$a1\t", ($has ? (split /\t/, $variant[$i])[2] : 0), "\n";
		print { $pass[$p][2] } "v$rec\n";
	}
	undef %layer;
	foreach my $pass (@pass) {
		my ($file, $fh, $ids) = @{ $pass };
		close $fh or croak "$sub: can't write $file.txt: $!";
		close $ids or croak "$sub: can't write $file.ids: $!";
		_plink2($sub, $plink2, $dir, '--pgen', "$g.pgen", '--pvar', "$g.ids.pvar", '--psam', "$g.psam", @common,
			'--extract', "$file.ids", '--score', "$file.txt", 1, 2, 3, 'header', 'no-mean-imputation', 'cols=+scoresums',
			'--missing', 'sample-only', '--out', $file);
		# .smiss's OBS_CT is how many genotypes there were to see, and
		# MISSING_CT how many of those were missing, so a sample's variants are
		# the one less the other.
		foreach my $part (['sscore', ALLELE_CT => ['alleles', 1], SCORE1_SUM => ['score', 1]],
			['smiss', OBS_CT => ['variants', 1], MISSING_CT => ['variants', -1]]) {
			my ($ext, %field) = @{ $part };
			open my $res, '<', "$file.$ext" or croak "$sub: can't read $file.$ext: $!";
			chomp(my $head = <$res>);
			my @h = split /\t/, $head;
			my %at_col = map { $h[$_] => $_ } 0 .. $#h;
			my $id_col = $at_col{'#IID'} // $at_col{IID};
			while (my $line = <$res>) {
				chomp $line;
				my @f = split /\t/, $line;
				$score{ $f[$id_col] }{ $field{$_}[0] } += $field{$_}[1] * $f[ $at_col{$_} ] foreach keys %field;
			}
			close $res;
		}
	}
	foreach my $s (values %score) {
		$s->{missing} = @variant - $s->{variants};
		$s->{mismatched} = $disagree;
		$s->{average} = $s->{alleles} ? $s->{score} / $s->{alleles} : undef;
	}
	# On 2026-10-17, PGS000138 (whose header says genome_build=NR) against a
	# GRCh37 genome had 12,015 of its 22,274 variants at their positions, and
	# the alleles agreed at all 12,015; with every position moved on one base,
	# as a stand-in for the wrong build, they agreed at 9 of the 58 found.
	# Disagreement at more than half the variants found is the second case.
	if (!$agree) {
		carp "$sub: none of the " . scalar(@variant) . ' variants of "weights" was found in "vcf" with its alleles; check that the two'
			. ($r{match} eq 'id' ? ' use the same rsIDs' : ' are on the same genome build and name chromosomes alike');
	} elsif ($disagree > $agree) {
		carp "$sub: the alleles disagree at $disagree of the " . ($agree + $disagree) . ' variants of "weights" found in "vcf"; check that'
			. ' the two are on the same ' . ($r{match} eq 'id' ? 'strand' : 'genome build');
	}
	require File::Path;
	File::Path::rmtree($dir);
	\%score;
}

# ---- VCF genotypes: sex and parentage

# The bytes sex_check and parent read at a time, from a file or through
# IO::Uncompress, for Bioinf.xs to count. On 2026-10-18 sex_check read
# genome.vcf.gz (155 MB, 3.5 million records) in 4.84 s with
# blocks of 16 KiB, 64 KiB or 256 KiB, and 4.90 s and 4.86 s with 1 MiB and 4
# MiB (best of three, peak RSS 47 MB each), so 64 KiB, as for FASTA. The time
# is gzip's: the same file uncompressed (928 MB) took 0.26 s and 11 MB, and
# zcat alone 3.7 s; as .vcf.bz2 it took 19.4 s, and bzip2 -dc alone 18.5 s.
# IO::Uncompress::Gunzip's readline, which is perl, took 14.3 s just to read
# the lines, and the read loop alone in 1 MiB blocks 5.7 s.
my $Vcf_block = 1 << 16;

# The VCF $file, given to $sub as its option $key, opened: returns a sub that
# gives the next block of it and then undef, and a sub that closes it. Its
# messages are prefixed with $sub, whose argument the file is.
sub _vcf_open {
	my ($sub, $key, $file) = @_;
	croak "$sub: \"$key\" must be a file name, not " . _what($file) if ref $file;
	croak "$sub: \"$key\" file \"$file\" doesn't exist or isn't readable" unless -f $file && -r _;
	my ($fh, $error); # $error: what says why a read failed
	# Transparent => 0, or IO::Uncompress would pass a file that isn't compressed
	# through as it is, and the error would be that it isn't a VCF. For a file too
	# short to hold a header, new() gives undef and leaves the error empty in
	# IO::Uncompress::Gunzip 2.223 (perl 5.44.0), though 2.020 (perl 5.10.1) says
	# "Header Error: Minimum header size is 10 bytes".
	if ($file =~ /\.vcf\.b?gz\z/i) {
		require IO::Uncompress::Gunzip;
		$fh = IO::Uncompress::Gunzip->new($file, MultiStream => 1, Transparent => 0)
			or croak "$sub: can't read $file: " . ($IO::Uncompress::Gunzip::GunzipError || 'it isn\'t gzip-compressed');
		$error = sub { $IO::Uncompress::Gunzip::GunzipError };
	} elsif ($file =~ /\.vcf\.bz2\z/i) {
		require IO::Uncompress::Bunzip2;
		$fh = IO::Uncompress::Bunzip2->new($file, MultiStream => 1, Transparent => 0)
			or croak "$sub: can't read $file: " . ($IO::Uncompress::Bunzip2::Bunzip2Error || 'it isn\'t bzip2-compressed');
		$error = sub { $IO::Uncompress::Bunzip2::Bunzip2Error };
	} elsif ($file =~ /\.vcf\z/i) {
		open $fh, '<:raw', $file or croak "$sub: can't read $file: $!";
		$error = sub { "$!" };
	} else {
		croak "$sub: \"$key\" file $file must end in .vcf, .vcf.gz, .vcf.bgz or .vcf.bz2";
	}
	my $next = sub {
		my $block;
		my $n = read $fh, $block, $Vcf_block;
		# IO::Uncompress's read gives a negative number on an error, perl's undef
		croak "$sub: can't read $file: " . $error->() if !defined $n || $n < 0;
		$n ? $block : undef;
	};
	my $close = sub { close $fh or croak "$sub: can't read $file: " . $error->() };
	($next, $close);
}

# sex_check's report on $vcf: what was counted and the rule, then a table of
# the samples in $order, the VCF's, each with its sex and the counts it was told
# from. Columns are as wide as their widest entry, numbers right-aligned.
sub _sex_report {
	my ($vcf, $r, $tally, $order) = @_;
	my @head = ('sample', 'sex', 'het/hom-ALT', 'X het', 'X hom-ALT', 'X haploid', 'X missing', 'Y called');
	my @rows = map {
		my $s = $tally->{$_};
		[$_, $s->{sex} // 'unknown', defined $s->{het_ratio} ? sprintf('%.3f', $s->{het_ratio}) : '-',
			@{ $s }{qw(x_het x_hom_alt x_haploid x_missing y_called)}]
	} @{ $order };
	my @w = map { length $head[$_] } 0 .. $#head;
	foreach my $row (@rows) {
		foreach my $i (0 .. $#head) {
			$w[$i] = length $row->[$i] if length $row->[$i] > $w[$i];
		}
	}
	my $line = sub { join('  ', map { $_ < 2 ? sprintf('%-*s', $w[$_], $_[$_]) : sprintf('%*s', $w[$_], $_[$_]) } 0 .. $#head) . "\n" };
	my @unknown = grep { !defined $tally->{$_}{sex} } @{ $order };
	"sex_check: $vcf\n"
		. "male when het/hom-ALT, over SNVs on X outside its pseudo-autosomal regions, is below $r->{cutoff}; genotypes of DP below $r->{min_depth} not counted\n"
		. $line->(@head) . join('', map { $line->(@{ $_ }) } @rows)
		. (@unknown ? scalar(@unknown) . ' of ' . scalar(@rows) . " unknown: no het or hom-ALT call on X\n" : '');
}

sub sex_check {
	my $sub = 'sex_check';
	# sex_check($file) and sex_check($file, %options) take the file first; an
	# odd count whose first is an option's name is left for _args to refuse
	unshift @_, 'vcf' if @_ % 2 && defined $_[0] && !ref $_[0] && $_[0] !~ /\A(?:vcf|cutoff|min_depth|report)\z/;
	croak "$sub: \"$_[0]\" is the file, so the options after it must be name => value pairs, and are an odd number"
		if @_ && $_[0] ne 'vcf' && defined $_[0] && !ref $_[0] && -f $_[0];
	my $args = _args(\@_, $sub, ['vcf'], [qw(cutoff min_depth report)]);
	my %r;
	# peddy 0.4.3's defaults (sex_check() in peddy/peddy.py): male below a het
	# to hom-ALT ratio of 0.6 on non-PAR X, and genotypes of depth 7 or more
	$r{cutoff} = $args->{cutoff} // 0.6;
	croak "$sub: \"cutoff\" must be a number above 0, not " . _what($r{cutoff})
		unless !ref $r{cutoff} && $r{cutoff} =~ $Number && $r{cutoff} > 0;
	$r{min_depth} = $args->{min_depth} // 7; # 0 = count every genotype, whatever its DP
	croak "$sub: \"min_depth\" must be a whole number, 0 or more, not " . _what($r{min_depth})
		unless !ref $r{min_depth} && $r{min_depth} =~ /\A[0-9]+\z/;
	$r{report} = $args->{report} // 1; # 1 = print a table of the samples to STDOUT; 0 = only return them
	croak "$sub: \"report\" must be 0 or 1, not " . _what($r{report}) unless !ref $r{report} && ($r{report} eq '0' || $r{report} eq '1');
	my $vcf = $args->{vcf};
	my ($next, $close) = _vcf_open($sub, 'vcf', $vcf);
	my ($tally, $x_sites, $order) = _sex_tally($next, $vcf, $r{min_depth});
	$close->();
	foreach my $s (values %{ $tally }) {
		# peddy divides by a hom-ALT count of 0 as if it were 1, and so calls a
		# sample with no calls on X male; here it is undef, as is its sex
		my $calls = $s->{x_het} + $s->{x_hom_alt};
		$s->{het_ratio} = $calls ? $s->{x_het} / ($s->{x_hom_alt} || 1) : undef;
		$s->{sex} = !$calls ? undef : $s->{het_ratio} < $r{cutoff} ? 'male' : 'female';
	}
	if ($r{report}) {
		print STDOUT _sex_report($vcf, \%r, $tally, $order) or croak "$sub: can't print the report: $!";
	}
	carp "$sub: \"vcf\" $vcf has no SNV on X outside its pseudo-autosomal regions, so no sex can be told" unless $x_sites;
	$tally;
}

# KING's bounds on the kinship coefficient, from its manual
# (www.kingrelatedness.com/manual.shtml, read 2026-10-18): above 0.354 the same
# person or identical twins, and then first-, second- and third-degree
# relatives down to 0.177, 0.0884 and 0.0442, below which a pair is unrelated.
# They are 2^-1.5, 2^-2.5, 2^-3.5 and 2^-4.5, the geometric means of the
# kinships of neighbouring degrees (1/2, 1/4, 1/8, 1/16, 1/32), and the
# rounded values are used, as KING prints them.
my @King = ([0.354, 'the same person or identical twins'], [0.177, 'first-degree relatives'],
	[0.0884, 'second-degree relatives'], [0.0442, 'third-degree relatives']);

# The fewest sites parent will decide from. At the IBS0 rate of the unrelated
# pair measured (0.081, in the README), 1,000 sites hold 81 IBS0 sites on
# average, with a binomial standard deviation of 8.6, so that falling below
# the 12 of the default max_ibs0 is 8 standard deviations away. An exome
# gives about 150,000 and a genome 3 million, so only a small panel, or
# files that hardly overlap, comes near it.
my $Parent_sites = 1000;

# parent's report, given its sample or file names, the resolved options, and
# its result.
sub _parent_report {
	my ($child, $parent, $vcf, $r, $s) = @_;
	my $kin = sprintf '%.4f', $s->{kinship};
	my $rate = sprintf '%.4f', $s->{ibs0_rate};
	my $band = (grep { $s->{kinship} > $King[$_][0] } 0 .. $#King)[0];
	my $kin_says = !defined $band ? "kinship $kin is below 0.0442, as for unrelated people"
		: $band == 0 ? "kinship $kin is above 0.354, as for $King[0][1]"
		: "kinship $kin is from $King[$band][0] to $King[$band - 1][0], as for $King[$band][1]";
	my $ibs0_says = "the IBS0 rate, $rate, is " . ($s->{ibs0_rate} < $r->{max_ibs0} ? 'below' : 'not below') . " $r->{max_ibs0}";
	my $verdict = $s->{parent} ? "parent and child (1): $ibs0_says, and $kin_says"
		: "not parent and child (0): $kin_says, and $ibs0_says";
	$verdict .= ', so full siblings' if $s->{relationship} eq 'full siblings';
	my @rows = (["sites compared", $s->{sites}, "autosomal biallelic SNVs where either carries the ALT, DP $r->{min_depth} or more"]);
	push @rows, (['  in both VCFs', $s->{in_both}, ''], ["  in the child's alone", $s->{child_only}, 'the parent taken as hom-REF'],
		["  in the parent's alone", $s->{parent_only}, 'the child taken as hom-REF']) unless defined $vcf;
	push @rows, (['IBS0', $s->{ibs0}, 'opposite homozygotes'],
		['  child hom-ALT', $s->{child_alt_parent_ref}, 'and the parent hom-REF'], ['  child hom-REF', $s->{child_ref_parent_alt}, 'and the parent hom-ALT'],
		['IBS1', $s->{ibs1}, 'one allele shared'], ['IBS2', $s->{ibs2}, 'both alleles shared'],
		['het in both', $s->{het_het}, ''], ['het', $s->{child_het}, "the child's, and $s->{parent_het} the parent's"],
		['hom-ALT', $s->{child_hom_alt}, "the child's, and $s->{parent_hom_alt} the parent's"],
		['kinship', $kin, 'KING-robust, between-family']);
	my ($w0, $w1) = (0, 0);
	foreach my $row (@rows) {
		$w0 = length $row->[0] if length $row->[0] > $w0;
		$w1 = length $row->[1] if length $row->[1] > $w1;
	}
	"parent: child " . (defined $vcf ? "$child, parent $parent, in $vcf" : "$child, parent $parent") . "\n"
		. "$verdict\n"
		. join('', map { my $l = sprintf '%-*s  %*s', $w0, $_->[0], $w1, $_->[1]; $l .= "  $_->[2]" if length $_->[2]; "$l\n" } @rows);
}

sub parent {
	my $sub = 'parent';
	my $args = _args(\@_, $sub, [qw(child parent)], [qw(vcf min_depth max_ibs0 report)]);
	my %r;
	# Measured on three exomes of one family on 2026-10-18 (README.md): with
	# every genotype counted the two parent-child pairs had IBS0 rates of 0.254
	# and 0.215, from off-target calls of DP 1 to 3, and with DP 10 or more
	# 0.0030 and 0.0020, against 0.081 for the unrelated pair; DP 7 gave 0.0053
	# and 0.0037, half the margin below max_ibs0, for 8% more sites.
	$r{min_depth} = $args->{min_depth} // 10; # 0 = count every genotype, whatever its DP
	croak "$sub: \"min_depth\" must be a whole number, 0 or more, not " . _what($r{min_depth})
		unless !ref $r{min_depth} && $r{min_depth} =~ /\A[0-9]+\z/;
	# peddy 0.4.3's (ped_check() in peddy/peddy.py: predicted_parents is
	# ibs0 / n < 0.012), where it sits between the measured parent-child rates
	# above and the 0.02 or so expected of full siblings, who share no allele
	# by descent at a quarter of sites and so have a quarter of the unrelated
	# pair's 0.081
	$r{max_ibs0} = $args->{max_ibs0} // 0.012;
	croak "$sub: \"max_ibs0\" must be a number above 0, not " . _what($r{max_ibs0})
		unless !ref $r{max_ibs0} && $r{max_ibs0} =~ $Number && $r{max_ibs0} > 0;
	$r{report} = $args->{report} // 1; # 1 = print the report to STDOUT; 0 = only return
	croak "$sub: \"report\" must be 0 or 1, not " . _what($r{report}) unless !ref $r{report} && ($r{report} eq '0' || $r{report} eq '1');
	my ($child, $parent, $vcf) = @{ $args }{qw(child parent vcf)}; # $vcf: undef = child and parent are each a VCF of its own
	my @t; # sites by the child's ALT count and then the parent's, 0 to 2, flattened; then the informative sites in both files, the child's alone, the parent's alone
	if (defined $vcf) {
		foreach my $key (qw(child parent)) {
			croak "$sub: \"$key\" must be a sample name of \"vcf\", not " . _what($args->{$key}) if ref $args->{$key} || !length $args->{$key};
		}
		croak "$sub: \"child\" and \"parent\" are both \"$child\"" if $child eq $parent;
		my ($next, $close) = _vcf_open($sub, 'vcf', $vcf);
		@t = _parent_one($next, $vcf, $r{min_depth}, $child, $parent);
		$close->();
	} else {
		croak "$sub: \"child\" and \"parent\" are both \"$child\"" if !ref $child && !ref $parent && $child eq $parent;
		my ($child_next, $child_close) = _vcf_open($sub, 'child', $child);
		my ($parent_next, $parent_close) = _vcf_open($sub, 'parent', $parent);
		@t = _parent_own($child_next, $child, $parent_next, $parent, $r{min_depth});
		$child_close->();
		$parent_close->();
	}
	my %s;
	@s{qw(in_both child_only parent_only)} = @t[9 .. 11];
	$s{sites} = $s{in_both} + $s{child_only} + $s{parent_only};
	@s{qw(child_ref_parent_alt child_alt_parent_ref)} = @t[2, 6];
	$s{ibs0} = $t[2] + $t[6];
	$s{ibs1} = $t[1] + $t[3] + $t[5] + $t[7];
	$s{ibs2} = $t[4] + $t[8];
	$s{het_het} = $t[4];
	$s{child_het} = $t[3] + $t[4] + $t[5];
	$s{parent_het} = $t[1] + $t[4] + $t[7];
	$s{child_hom_alt} = $t[6] + $t[7] + $t[8];
	$s{parent_hom_alt} = $t[2] + $t[5] + $t[8];
	my $what = defined $vcf ? "\"$child\" and \"$parent\" of $vcf" : "$child and $parent";
	# Without a hom-REF call, IBS0 cannot be seen, and an unrelated pair would
	# look like parent and child. A VCF merged from VCFs that list only their own
	# sample's variant sites is like that, as bcftools merge writes it without
	# -0: each sample is "./." wherever it had no variant. It is checked before
	# the number of sites, which such a VCF often has too few of as well, since
	# this says what to do about it.
	croak "$sub: neither of $what has a hom-REF call where the other carries the ALT, so no opposite homozygote can be seen;"
		. ' a VCF merged from VCFs that list only their own variants is like that, so give each sample\'s own VCF as "child" and "parent"'
		if $s{sites} && !($t[1] + $t[2] + $t[3] + $t[6]);
	croak "$sub: $what have only $s{sites} autosomal biallelic SNVs of DP $r{min_depth} or more where either carries the ALT,"
		. " and $Parent_sites are needed" . (defined $vcf ? '' : '; if the files share few sites, check that they are of the same genome build')
		if $s{sites} < $Parent_sites;
	my $min_het = $s{child_het} < $s{parent_het} ? $s{child_het} : $s{parent_het};
	croak "$sub: " . ($s{child_het} ? 'the parent' : 'the child') . " of $what has no het call at the sites compared, so no kinship can be taken" unless $min_het;
	# KING-robust between-family kinship (Manichaikul et al. 2010,
	# Bioinformatics 26:2867), in the form PLINK 2's --make-king-table gives
	# it: t/parent.t's header says how the two were found to agree
	$s{kinship} = 0.5 - ($s{child_het} + $s{parent_het} - 2 * $s{het_het} + 4 * $s{ibs0}) / (4 * $min_het);
	$s{ibs0_rate} = $s{ibs0} / $s{sites};
	my $band = (grep { $s{kinship} > $King[$_][0] } 0 .. $#King)[0];
	$s{relationship} = !defined $band ? 'unrelated' : $band == 0 ? 'duplicate'
		: $band == 1 ? ($s{ibs0_rate} < $r{max_ibs0} ? 'parent-child' : 'full siblings')
		: $band == 2 ? 'second-degree' : 'third-degree';
	$s{parent} = $s{relationship} eq 'parent-child' ? 1 : 0;
	if ($r{report}) {
		print STDOUT _parent_report($child, $parent, $vcf, \%r, \%s) or croak "$sub: can't print the report: $!";
	}
	# Two VCFs of one pipeline have about as many sites each alone: on
	# 2026-10-18 the three pairs of one family's exomes had ratios of 1.13 at
	# most, where an exome of another capture kit had 3.0 to 13 against them,
	# and a genome against an exome 58, which is why 2.
	my ($more, $fewer) = $s{child_only} > $s{parent_only} ? ('child', 'parent') : ('parent', 'child');
	my ($n_more, $n_fewer) = @s{"${more}_only", "${fewer}_only"};
	carp "$sub: the ${more}'s VCF has $n_more sites that the ${fewer}'s lacks, and the ${fewer}'s only $n_fewer that the ${more}'s lacks,"
		. ' so the two were probably not sequenced over the same regions (an exome and a genome, or two capture kits);'
		. ' a site that one lacks is taken as hom-REF in it, so IBS0 and the kinship are not to be trusted'
		if !defined $vcf && $n_more > 2 * $n_fewer;
	wantarray ? ($s{parent}, \%s) : $s{parent};
}

# ---- alignment plots and tables

# Dies unless $n, given to $sub as $what, is a 1-based residue number. Without
# this, 0 and the negative numbers index @col from its end in plot_msa and
# clustal_view_residues, and silently mean the last residues.
sub _residue_number {
	my ($sub, $what, $n) = @_;
	croak "$sub: $what must be a residue number, 1 or more, not " . (defined $n ? "\"$n\"" : 'undef')
		unless defined $n && $n =~ /\A[1-9][0-9]*\z/;
}

# "fasta" -- a FASTA file name or a hash ref of name => sequence -- written to a
# temporary FASTA file with its gaps stripped, so that an input that is already
# aligned is aligned afresh. clustalo's --dealign does the same, but prints a
# "FORCED DEBUG" warning about them first; blastp would take a gap for a residue.
# Returns that file and the sequence names in input order (sorted for a hash).
# With $numbered, each sequence is written under its number in that order, "s1",
# "s2" and so on, rather than its name, for clustalo: it cuts a name at its
# first space in the guide tree and writes "(),:" in one unquoted, which drew 4
# sequences as 5 tips up to 2026-10-17. _unnumber_fasta and _unnumber_newick put
# the names back.
sub _ungapped_fasta {
	my ($fasta, $numbered) = @_;
	my ($seqs, @names);
	if (ref $fasta eq 'HASH') {
		$seqs = $fasta;
		@names = sort keys %{ $seqs };
	} elsif (ref $fasta eq '') {
		my $order;
		($seqs, $order) = _fasta_ordered($fasta);
		@names = @{ $order };
	} else {
		croak "_ungapped_fasta: \"fasta\" must be a FASTA file name or a hash ref of name => sequence, not " . _what($fasta);
	}
	my @written = $numbered ? map { "s$_" } 1 .. @names : @names;
	my %ungapped;
	@ungapped{@written} = map { (my $s = $seqs->{$_}) =~ tr/-.//d; $s } @names;
	(hash2fasta_file(\%ungapped, _tmp('.fa'), \@written), \@names);
}

# The name that clustalo's "s<n>" in $file stands for, from @$names, the list
# _ungapped_fasta numbered.
sub _unnumber {
	my ($id, $names, $file) = @_;
	croak "_unnumber: clustalo wrote \"$id\" in $file, which isn't one of the s1 to s" . scalar(@{ $names }) . ' it was given'
		unless $id =~ /\As([1-9][0-9]*)\z/ && $1 <= @{ $names };
	$names->[$1 - 1];
}

# The alignment clustalo wrote to $file under _ungapped_fasta's numbers,
# rewritten in place under the names in @$names, as UTF-8, in clustalo's order
# and 60-column lines. Returns the aligned sequences as a hash ref keyed by
# those names as perl was given them, which a name's UTF-8 may not equal.
sub _unnumber_fasta {
	my ($file, $names) = @_;
	my ($numbered, $ids) = _fasta_ordered($file);
	my (%aln, %utf8);
	foreach my $id (@{ $ids }) {
		my $name = _unnumber($id, $names, $file);
		$aln{$name} = $utf8{ _utf8($name) } = $numbered->{$id};
	}
	hash2fasta_file(\%utf8, $file, [map { _utf8(_unnumber($_, $names, $file)) } @{ $ids }], 60);
	\%aln;
}

# The guide tree clustalo wrote to $file under _ungapped_fasta's numbers,
# rewritten in place under the names in @$names, as UTF-8. A name holding
# whitespace or any of "()[]':;," is quoted as newick quotes one: in single
# quotes, each of its own doubled.
sub _unnumber_newick {
	my ($file, $names) = @_;
	open my $fh, '<:raw', $file or croak "_unnumber_newick: can't read $file: $!";
	my $tree = do { local $/; <$fh> };
	close $fh or croak "_unnumber_newick: can't read $file: $!";
	$tree =~ s{(?<=[(,\s])(s[0-9]+)(?=\s*[:,);])}{
		my $name = _utf8(_unnumber($1, $names, $file));
		$name =~ /[\s()\[\]':;,]/ ? "'" . join("''", split /'/, $name, -1) . "'" : $name;
	}ge;
	open $fh, '>:raw', $file or croak "_unnumber_newick: can't write $file: $!";
	print {$fh} $tree;
	close $fh or croak "_unnumber_newick: can't write $file: $!";
}

# Aligns "fasta" with clustalo for $sub. The alignment goes to "msa_file" and
# the guide tree to "tree_file" when those are given, and to temporary files
# otherwise; the tree is not asked for at all unless it is kept or $want_tree.
# Returns the hash ref of files to report (the alignment always, the tree only
# when it was named), the tree's file, the sequence names in input order
# (sorted for a hash), an array ref of the clustalo command as it was run, with
# the ungapped FASTA it read as its "--in", and the aligned sequences as a hash
# ref keyed by those names; or nothing, after a warning, for fewer than 2
# sequences. clustalo is given the sequences numbered, and the alignment and
# tree it writes get their names back.
sub _align {
	my ($args, $sub, $want_tree) = @_;
	my ($in, $names) = _ungapped_fasta($args->{fasta}, 1);
	if (@{ $names } < 2) {
		carp "$sub: clustalo needs at least 2 sequences, and there " . (@{ $names } ? 'is only 1' : 'are none');
		return;
	}
	my %r = ('msa_file' => $args->{'msa_file'} // _tmp('.aln.fa'));
	my $tree = $args->{'tree_file'} // ($want_tree ? _tmp('.newick') : undef);
	$r{'tree_file'} = $tree if defined $args->{'tree_file'};
	my @cmd = (_alien($sub)->clustalo, '--in', $in, '--out', $r{'msa_file'}, '--outfmt=fa', '--force',
		'--threads', $args->{threads} // 1, (defined $tree ? "--guidetree-out=$tree" : ()), @{ $args->{'clustal_args'} // [] });
	_run(@cmd);
	my $aln = _unnumber_fasta($r{'msa_file'}, $names);
	_unnumber_newick($tree, $names) if defined $tree;
	(\%r, $tree, $names, \@cmd, $aln);
}

sub plot_msa {
	my $sub = 'plot_msa';
	my $args = _args(\@_, $sub, [qw(fasta filename)], [qw(active_site_aa clustal_args labels msa_file order query threads title
		tree_file xlabel ylabel)]);
	_refs($sub, $args, 'ARRAY', qw(clustal_args order));
	_refs($sub, $args, 'HASH', 'labels');
	my $sites = $args->{'active_site_aa'};
	if (defined $sites) {
		croak "$sub: \"active_site_aa\" must be a hash ref of label => residue number" unless ref $sites eq 'HASH';
		croak "$sub: needs \"query\" to place \"active_site_aa\"" unless defined $args->{query};
		_residue_number($sub, "active site $_", $sites->{$_}) foreach sort keys %{ $sites };
	}
	my ($r, undef, $names, $cmd, $aln) = _align($args, $sub, 0) or return {};
	my @order = @{ $args->{order} // $names };
	my @unknown = grep { !defined $aln->{$_} } @order;
	croak "$sub: \"order\" names sequences that aren't in the alignment: @unknown" if @unknown;
	my $labels = $args->{labels} // {};
	# Each name as msa_plot.py is to show it, in UTF-8, which is what it reads.
	# Up to 2026-10-17 the name was written as perl held it, Latin-1 or UTF-8,
	# and Python died on the Latin-1.
	my %shown_as = map { $_ => _utf8($labels->{$_} // $_) } @order;
	# Two sequences shown under one label would be one row of the image, fed
	# two deflines.
	my %by_label;
	push @{ $by_label{ $shown_as{$_} } }, $_ foreach @order;
	foreach my $label (sort keys %by_label) {
		croak "$sub: " . join(', ', map { "\"$_\"" } @{ $by_label{$label} }) . " would all be drawn as \"$label\""
			if @{ $by_label{$label} } > 1;
	}
	my @segments;
	if (defined $sites) {
		my $q = $aln->{ $args->{query} } // croak "$sub: the query \"$args->{query}\" isn't in the alignment";
		# 1-based residue number => 0-based alignment column
		my (@col, %sites);
		while ($q =~ /[^-]/g) {
			push @col, pos($q) - 1;
		}
		foreach my $label (sort keys %{ $sites }) {
			my $n = $sites->{$label};
			$sites{$label} = $col[$n - 1] // croak "$sub: active site $label ($n) is past the end of \"$args->{query}\", which has " . scalar(@col) . ' residues';
		}
		@segments = ('--s', _json_text({ map { _chars($_) => 0 + $sites{$_} } keys %sites }));
	}
	my %shown = map { $shown_as{$_} => $aln->{$_} } @order;
	my $fa = hash2fasta_file(\%shown, _tmp('.fa'), [@shown_as{@order}], 0);
	_python($sub, '--f', $fa, '--o', $args->{filename}, '--c', _creator($sub), '--quiet',
		'--p', _json_text([map { _chars($_) } _align_provenance($args, $sub, $cmd, $names, $r->{'msa_file'}), _kept($args, qw(msa_file tree_file))]),
		(defined $args->{title} ? ('--t', $args->{title}) : ()),
		(defined $args->{xlabel} ? ('--x', $args->{xlabel}) : ()),
		(defined $args->{ylabel} ? ('--y', $args->{ylabel}) : ()), @segments);
	_wrote($args->{filename}, 43); # yellow
	$r->{filename} = $args->{filename};
	$r;
}

# "path (SHA-256 hex)" of an existing file, for the provenance of an image.
sub _file_id {
	my ($file) = @_;
	require Digest::SHA;
	require File::Spec;
	File::Spec->rel2abs($file) . ' (SHA-256 ' . Digest::SHA->new(256)->addfile($file, 'b')->hexdigest . ')';
}

# "fasta" as an image's Creator names it: a file by its path and digest, a hash
# ref by its size, since its contents were never on disk until
# _ungapped_fasta wrote them.
sub _fasta_id {
	my ($fasta, $n) = @_;
	ref $fasta ? "a hash ref of $n sequences" : 'the FASTA file ' . _file_id($fasta) . " of $n sequences";
}

# The files named in @keys of $args that were written and kept beside an
# image, as a clause of its Creator, or nothing when none was.
sub _kept {
	my ($args, @keys) = @_;
	require File::Spec;
	my @kept = map { defined $args->{$_} && !ref $args->{$_} ? "$_ " . File::Spec->rel2abs($args->{$_}) : () } @keys;
	@kept ? 'kept with it: ' . join(', ', @kept) : ();
}

# A command as it could be pasted into a shell, for an image's Creator.
sub _cmd_text { join ' ', map { /[^\w\/.,:=+-]/ ? "'$_'" : $_ } @_ }

# How _align made an alignment, as one clause of an image's Creator: what
# Clustal Omega read, which version ran, what it wrote and the command. $cmd,
# $names and $msa are what _align returned.
sub _align_provenance {
	my ($args, $sub, $cmd, $names, $msa) = @_;
	my $alien = _alien($sub);
	# $cmd->[2] is clustalo's "--in", the ungapped copy _ungapped_fasta wrote
	'from ' . _fasta_id($args->{fasta}, scalar @{ $names }) . ', written with its gaps stripped and the sequences numbered s1 to s'
		. scalar(@{ $names }) . ' to ' . _file_id($cmd->[2]) . ', which Clustal Omega ' . $alien->versions->{clustalo}
		. " from Alien::Bioinf $Alien::Bioinf::VERSION aligned into " . _file_id($msa) . ' (with the names put back), run as: '
		. _cmd_text(@{ $cmd });
}

sub plot_phylo {
	my $sub = 'plot_phylo';
	my $args = _args(\@_, $sub, [], [qw(clustal_args fasta labels msa_file output_file threads title tree_file)]);
	_refs($sub, $args, 'ARRAY', 'clustal_args');
	_refs($sub, $args, 'HASH', 'labels');
	$args->{output_file} //= 'phylo.svg';
	my ($r, $tree, $names, $cmd);
	if (defined $args->{fasta}) {
		($r, $tree, $names, $cmd) = _align($args, $sub, 1) or return {};
	} else {
		# no alignment to make, so "tree_file" is the tree to draw, not where to keep one
		$tree = $args->{'tree_file'} // croak "$sub: needs \"fasta\" to align, or \"tree_file\" to draw";
		my @moot = grep { defined $args->{$_} } qw(clustal_args msa_file threads);
		croak "$sub: was given " . join(', ', map { "\"$_\"" } @moot) . ' but no "fasta" to align' if @moot;
		croak "$sub: \"$tree\" doesn't exist or isn't a readable file" unless -f $tree && -r _;
		$r = { 'tree_file' => $tree };
	}
	my $labels = $args->{labels} // {};
	my @prov = defined $cmd
		? (_align_provenance($args, $sub, $cmd, $names, $r->{'msa_file'}), 'drawing the guide tree it wrote, ' . _file_id($tree),
			_kept($args, qw(msa_file tree_file)))
		: ('from the newick file ' . _file_id($tree));
	# the tips are named in UTF-8 by _unnumber_newick, which msa_plot.py decodes
	_python($sub, '--tree', $tree, '--o', $args->{'output_file'}, '--c', _creator($sub), '--quiet', '--p', _json_text([map { _chars($_) } @prov]),
		(defined $args->{title} ? ('--t', $args->{title}) : ()),
		(%{ $labels } ? ('--l', _json_text({ map { _chars($_) => _chars($labels->{$_}) } keys %{ $labels } })) : ()));
	$r->{'output_file'} = $args->{'output_file'};
	_wrote($args->{output_file}, 43); # yellow
	$r;
}

sub msa_quality_table {
	my $sub = 'msa_quality_table';
	my $args = _args(\@_, $sub, ['filename'], [qw(alignment_json cb_label cb_max cb_min cblogscale default_undefined fasta
		logscale_add metric msa_file normalize order show_numbers title unaligned_fa)]);
	_refs($sub, $args, 'ARRAY', 'order');
	my %metric = map { $_ => 1 } qw(num bit_score score evalue identity positive align_len);
	my $metric = $args->{metric} // 'score';
	croak "$sub: \"$metric\" isn't one of the metrics: " . join(', ', sort keys %metric) unless $metric{$metric};
	# "unaligned_fa" is the old name for "fasta"
	croak "$sub: was given both \"fasta\" and \"unaligned_fa\", which are the same thing" if defined $args->{fasta} && defined $args->{'unaligned_fa'};
	my $fasta = $args->{fasta} // $args->{'unaligned_fa'};
	# All-against-all blastp of the sequences: given as the parsed report, read
	# from "alignment_json", or -- when that file does not exist yet, or none is
	# named -- made by running blastp on "fasta", and kept in "alignment_json"
	# for next time if that is named.
	my $aj = $args->{'alignment_json'};
	croak "$sub: \"alignment_json\" must be a file name or a hash ref, not " . _what($aj) if ref $aj && ref $aj ne 'HASH';
	my ($blast, $prov);
	if (ref $aj eq 'HASH') {
		$blast = $aj;
		$prov = 'from a BLAST report given as a hash ref';
	} elsif (defined $aj && -f $aj) {
		$blast = _json_file($aj);
		$prov = 'from the BLAST report ' . _file_id($aj);
	} else {
		croak "$sub: needs \"fasta\" to align, or an existing \"alignment_json\"" . (defined $aj ? " ($aj doesn't exist yet)" : '')
			unless defined $fasta;
		my ($in, $names) = _ungapped_fasta($fasta);
		$aj //= _tmp('.json');
		my $alien = _alien($sub);
		my @cmd = ($alien->blast('blastp'), '-query', $in, '-subject', $in, '-out', $aj, '-outfmt', 15);
		_run(@cmd);
		$blast = _json_file($aj);
		$prov = 'from ' . _fasta_id($fasta, scalar @{ $names }) . ', written with its gaps stripped to ' . _file_id($in)
			. ', which blastp ' . $alien->versions->{blast} . " from Alien::Bioinf $Alien::Bioinf::VERSION compared all against all into "
			. _file_id($aj) . ', run as: ' . _cmd_text(@cmd);
		$prov = join '; ', $prov, _kept($args, 'alignment_json');
	}
	my (%data, $max);
	foreach my $query (@{ $blast->{BlastOutput2} }) {
		foreach my $hit_list (@{ $query->{report}{results}{bl2seq} }) {
			foreach my $hit (@{ $hit_list->{hits} }) {
				my $value = $hit->{hsps}[0]{$metric} // next;
				$data{ $hit_list->{query_title} }{ $hit->{description}[0]{title} } = $value;
				$max = $value if !defined $max || $value > $max;
			}
		}
	}
	croak "$sub: no \"$metric\" values in the BLAST report" unless defined $max;
	my @order = @{ $args->{order} // [sort { lc $a cmp lc $b } keys %data] };
	my (@row_labels, @names);
	foreach my $key (@order) {
		my ($genus, $species) = $key =~ /^([^.]+)\.(.+)/ or croak "$sub: can't get a genus and species from \"$key\"";
		my ($sp, @strain) = split /\./, $species; # C.neoformans.B.3501A
		push @names, [$genus, $sp, "@strain"];
		push @row_labels, '$\it{' . ucfirst(_first_letter($genus)) . ". $sp}\$" . (@strain ? " @strain" : '');
	}
	my @col_labels = @{ _abbreviated_labels(\@names) };
	my $add = $args->{'logscale_add'} // 0;
	my $norm = ($args->{normalize} // 0) > 0;
	# Normalized, every value is divided by the largest, with "logscale_add"
	# added to both so that the largest is still 1. That needs a largest above
	# 0, which "evalue" does not have when BLAST has rounded every e-value to 0.
	croak "$sub: can't normalize \"$metric\", since its largest value" . ($add ? ' plus "logscale_add"' : '') . ' is ' . ($max + $add) . ', not above 0'
		if $norm && $max + $add <= 0;
	my $log = $args->{cblogscale};
	my (@cells, $lo, $hi, $lo_positive);
	foreach my $i (0 .. $#order) {
		foreach my $j (0 .. $#order) {
			# a pair with no hit stays out of the table, rather than claiming a
			# score of 0, and is drawn grey
			my $value = $data{ $order[$i] }{ $order[$j] } // $args->{default_undefined};
			if (defined $value) {
				$value += $add;
				$value /= $max + $add if $norm;
				$lo = $value if !defined $lo || $value < $lo;
				$hi = $value if !defined $hi || $value > $hi;
				$lo_positive = $value if $value > 0 && (!defined $lo_positive || $value < $lo_positive);
			}
			$cells[$i][$j] = $value;
		}
	}
	croak "$sub: no pair of the sequences in \"order\" has a \"$metric\" value, or a \"default_undefined\", to draw" unless defined $hi;
	# A normalized scale runs from 0 to 1, but a log scale can't start at 0 or
	# below, so it starts at the smallest value above 0 instead, as it does
	# when not normalized.
	@{ $args }{qw(cb_min cb_max)} = ($log ? undef : 0, 1) if $norm;
	$lo = $args->{cb_min} // ($log ? $lo_positive : $lo);
	$hi = $args->{cb_max} // $hi;
	croak "$sub: \"cblogscale\" needs a value above 0 to start the scale at, or a \"cb_min\" above 0" if $log && (!defined $lo || $lo <= 0);
	require JSON::MaybeXS;
	my $table = _tmp('.json');
	open my $fh, '>:raw', $table or croak "$sub: can't write $table: $!";
	# The strings go through _chars, since encode_json reads a byte string as
	# Latin-1: up to 2026-10-17 a UTF-8 title "caf\xc3\xa9" was drawn as
	# "caf\xc3\x83\xc2\xa9".
	print {$fh} JSON::MaybeXS::encode_json({
		cells => \@cells, cols => [map { _chars($_) } @col_labels], rows => [map { _chars($_) } @row_labels], vmin => $lo + 0, vmax => $hi + 0,
		log => $log ? JSON::MaybeXS::true() : JSON::MaybeXS::false(),
		numbers => $args->{'show_numbers'} ? JSON::MaybeXS::true() : JSON::MaybeXS::false(),
		title => _chars($args->{title} // ''), cblabel => (defined $args->{cb_label} ? _chars($args->{cb_label}) : undef),
	});
	close $fh or croak "$sub: can't write $table: $!";
	_python($sub, '--table', $table, '--o', $args->{filename}, '--c', _creator($sub), '--p', _json_text([map { _chars($_) } "showing each pair's \"$metric\"", $prov]));
	$args->{filename};
}

# $text with every character that LaTeX treats specially written so that it
# prints as itself. "\^" would be a circumflex accent, not a caret, and "\~" a
# tilde accent, so those two (and the backslash) are spelled as text commands.
sub _latex_text {
	my ($text) = @_;
	my %special = ('\\' => '\\textbackslash{}', '^' => '\\textasciicircum{}', '~' => '\\textasciitilde{}');
	$text =~ s/([\\^~#\$%&_{}])/$special{$1} \/\/ "\\$1"/ge;
	$text;
}

sub clustal_view_residues {
	my $sub = 'clustal_view_residues';
	# msa_file: a FASTA file, aligned (such as plot_msa's) or not; color_residues:
	# {protein => {1-based residue number => colour}}, the colour an xcolor
	# name or an [r, g, b] ref; order: proteins top to bottom; row_width:
	# columns per block (100); split: blocks per LaTeX table (4); track: a
	# protein whose coloured residue numbers get a row of their own; threads and
	# clustal_args: for clustalo, if msa_file has to be aligned.
	my $args = _args(\@_, $sub, ['msa_file', 'output_tex_file'], [qw(caption clustal_args color_residues label
		order row_width split table_text_size threads track)]);
	_refs($sub, $args, 'ARRAY', qw(clustal_args order));
	my $color = $args->{'color_residues'} // {};
	croak "$sub: \"color_residues\" must be a hash ref" unless ref $color eq 'HASH';
	foreach my $protein (sort keys %{ $color }) {
		croak "$sub: \"color_residues\" must map $protein to a hash ref of residue number => colour, not " . _what($color->{$protein})
			unless ref $color->{$protein} eq 'HASH';
	}
	# 0 would never end the loop over blocks, and divide by 0 in the tables
	foreach my $k ('row_width', 'split') {
		croak "$sub: \"$k\" must be a whole number, 1 or more, not \"$args->{$k}\"" if defined $args->{$k} && $args->{$k} !~ /\A[1-9][0-9]*\z/;
	}
	croak "$sub: \"msa_file\" $args->{'msa_file'} doesn't exist" unless -e $args->{'msa_file'};
	my $data = fasta2hash($args->{'msa_file'});
	croak "$sub: has no sequences to show in $args->{'msa_file'}" unless %{ $data };
	my %len = map { length $_ => 1 } values %{ $data };
	if (keys %len > 1) {
		# Sequences of different lengths can't be an alignment, so align them,
		# into a temporary file: "msa_file" is the input, and is never written.
		my %clustal = map { $_ => $args->{$_} } grep { defined $args->{$_} } qw(clustal_args threads);
		(undef, undef, undef, undef, $data) = _align({ %clustal, fasta => $args->{'msa_file'} }, $sub, 0);
		%len = map { length $_ => 1 } values %{ $data };
	}
	my $track = $args->{track};
	croak "$sub: tracker \"$track\" isn't in the alignment" if defined $track && !defined $data->{$track};
	my ($aln_len) = keys %len;
	my @undef = grep { !defined $data->{$_} } sort keys %{ $color };
	croak "$sub: \"color_residues\" names proteins that aren't in the alignment: @undef" if @undef;
	my @proteins = @{ $args->{order} // [sort { lc $a cmp lc $b } keys %{ $data }] };
	my @bad = grep { !defined $data->{$_} } @proteins;
	croak "$sub: \"order\" names proteins that aren't in the alignment: @bad" if @bad;
	my $width = $args->{'row_width'} // 100;
	# Alignment column => colour, and for the tracked protein column => its
	# residue number. A coloured column is coloured in every protein.
	my (%col_color, %col_number);
	foreach my $protein (sort keys %{ $color }) {
		my @col;
		while ($data->{$protein} =~ /[A-Za-z]/g) {
			push @col, pos($data->{$protein}) - 1;
		}
		foreach my $n (sort { $a <=> $b } keys %{ $color->{$protein} }) {
			my $c = $color->{$protein}{$n};
			_residue_number($sub, "a residue of $protein in \"color_residues\"", $n);
			my $col = $col[$n - 1] // croak "$sub: $protein has no residue $n";
			croak "$sub: the colour of $protein residue $n must be a name or 3 numbers" if ref $c && (ref $c ne 'ARRAY' || @{ $c } != 3);
			$col_color{$col} = $c;
			$col_number{$col} = $n if defined $track && $protein eq $track;
		}
	}
	my (@table, %count);
	for (my $start = 0; $start < $aln_len; $start += $width) {
		my $w = $aln_len - $start < $width ? $aln_len - $start : $width;
		foreach my $protein (@proteins) {
			my @seq = split //, substr($data->{$protein}, $start, $w);
			$count{$protein} += grep { /[A-Za-z]/ } @seq;
			foreach my $col (grep { $_ >= $start && $_ < $start + $w } keys %col_color) {
				my ($c, $i) = ($col_color{$col}, $col - $start);
				$seq[$i] = ref $c ? '{\color[rgb]{' . join(',', @{ $c }) . "}$seq[$i]}" : "\\textcolor{$c}{$seq[$i]}";
			}
			my $name = _latex_text($protein);
			push @table, ["\\textit{$name}", '\texttt{' . join('', @seq) . '}', $count{$protein}];
			next unless defined $track && $protein eq $track;
			# Each residue number written from its own column rightwards, one
			# digit per column so the row stays aligned, and moved right past
			# the end of the previous number rather than overwriting it.
			my @row = ('-') x $w;
			my $free = 0; # first column not yet taken by a number
			foreach my $col (sort { $a <=> $b } grep { $_ >= $start && $_ < $start + $w } keys %col_number) {
				my $at = $col - $start < $free ? $free : $col - $start;
				my @digits = split //, $col_number{$col};
				last if $at + @digits > $w;
				@row[$at .. $at + $#digits] = @digits;
				$free = $at + @digits + 1;
			}
			push @table, ["\\textit{$name} track", '\texttt{' . join('', @row) . '}', $count{$protein}];
		}
		push @table, ['\hline'];
	}
	my $per_table = ($args->{split} // 4) * (1 + @proteins + (defined $track ? 1 : 0));
	my $size = $args->{'table_text_size'} // '\footnotesize';
	my $caption = $args->{caption} // '';
	open my $tex, '>', $args->{'output_tex_file'} or croak "$sub: can't write $args->{'output_tex_file'}: $!";
	print {$tex} "%written by $0, calling $sub in " . __FILE__ . "\n";
	my $n_tables = int((@table + $per_table - 1) / $per_table);
	foreach my $t (0 .. $n_tables - 1) {
		my $end = ($t + 1) * $per_table - 1;
		$end = $#table if $end > $#table;
		print {$tex} "\\begin{table}[htp]$size\n\\begin{tabular}{|c|l|c|} \\hline\n",
			'\textbf{Track} & \textbf{Sequence} & \textbf{Length}\\\\ \hline', "\n";
		print {$tex} join(' & ', @{ $_ }), ($_->[-1] eq '\hline' ? "\n" : "\n\\\\\n") foreach @table[$t * $per_table .. $end];
		print {$tex} "\\end{tabular}\n\\caption{$caption", ($t > 0 && $caption ne '' ? ' (continued)' : ''), "}\n";
		print {$tex} "\\label{tab:$args->{label}", ($n_tables > 1 ? $t : ''), "}\n" if defined $args->{label};
		print {$tex} "\\end{table} \\FloatBarrier\n";
	}
	close $tex or croak "$sub: can't write $args->{'output_tex_file'}: $!";
	_wrote($args->{'output_tex_file'}, 46); # cyan
	$args->{'output_tex_file'};
}

1;
__END__

=encoding utf8

=head1 SYNOPSIS

 use Cn::Bioinf ':all';

 my $seqs = fasta2hash('proteome.fa.gz');        # { defline => sequence }
 my $one  = fasta2hash('proteome.fa', 'P12345'); # just that sequence
 hash2fasta_file($seqs, 'copy.fa');

 my $hits = get_best_alignment_hit('blast.json', 'bit_score');

 # one clustalo run: plot_msa keeps the guide tree, and plot_phylo draws it
 plot_msa(
     fasta       => 'orthologs.fa',     # or { name => sequence }
     filename    => 'msa.svg',          # .png, .pdf, ... too
     'tree_file' => 'orthologs.newick',
     title       => 'EF-3',
 );
 plot_phylo('tree_file' => 'orthologs.newick', 'output_file' => 'tree.svg', title => 'EF-3');

 # a polygenic score for each sample of a VCF, from a PGS Catalog scoring file
 my $scores = pgs_score(vcf => 'genome.vcf.gz', weights => 'PGS000138.txt.gz');

 # each sample's sex, from its heterozygosity on X
 my $sex = sex_check(vcf => 'cohort.vcf.bz2'); # { sample => { sex => 'female', ... } }

 # whether two people are parent and child: 1 or 0
 my $is = parent(child => 'child.vcf.gz', parent => 'father.vcf.gz');

=head1 DESCRIPTION

FASTA I/O in XS, BLAST hit ranking, and multiple-sequence-alignment plots and
LaTeX tables, taken from the maintainer's C<bioinf.pm>; polygenic scores
from a VCF and a PGS Catalog scoring file; and, from VCFs, each sample's sex
and whether two people are parent and child.

Nothing is exported by default; ask for functions by name or with C<:all>.
Every function dies (via C<Carp::croak>) on bad arguments, with a message that
starts with the name of the function that raised it, as in
C<fasta2hash: couldn't find "four" in x.fa>. C<pgs_score>, C<sex_check>,
C<parent>, C<plot_msa>, C<plot_phylo>, C<msa_quality_table> and C<clustal_view_residues>
take C<< name =E<gt> value >> pairs, not a hash ref.

=head2 Alien::Bioinf

Clustal Omega, BLAST+, PLINK 2 and the Python that draws the plots come from
C<Alien::Bioinf>, which installs them module-locally and never uses C<PATH>. It
is recommended rather than required, because it can only be installed where
NCBI builds BLAST+ (Linux and macOS on x86_64 and aarch64, and Windows on
x86_64) and it needs the network; PLINK 2 has no build for Linux on aarch64,
so there it installs the rest. C<fasta2hash>, C<hash2fasta_file> and
C<get_best_alignment_hit> work without it; C<pgs_score>, C<plot_msa>,
C<plot_phylo>, C<msa_quality_table>, and C<clustal_view_residues> given an
unaligned file, die saying that they need it. A CPAN client does not install a recommended module
unless asked, so install it first:

 cpanm Alien::Bioinf
 cpanm Cn::Bioinf

or, from a clone of L<the repository|https://github.com/hhg7/Bioinf-Basic>,
where it is C<alien/Alien-Bioinf>:

 cpanm ./alien/Alien-Bioinf
 cpanm .

That downloads nothing it can reuse: archives are cached (C<ALIEN_BIOINF_CACHE>,
default C<~/.cache/alien-bioinf>), a clustalo or BLAST+ already on the machine at
the newest version is copied instead (C<ALIEN_BIOINF_CLUSTALO>,
C<ALIEN_BIOINF_BLAST>, or C<PATH>), and pip reuses the wheels in its own cache.
The Python packages go into a venv of their own, which nothing installed into
the base Python can break. To check for and apply updates:

 perl -MAlien::Bioinf -MData::Dumper -e 'print Dumper(Alien::Bioinf->check_updates)'
 perl -MAlien::Bioinf -e 'print "$_\n" for Alien::Bioinf->update'

=head2 Provenance in the images

Every PNG, SVG, PDF, PS or EPS image these functions draw carries its whole
provenance on one line, as its C<Creator> metadata and nothing else, just as
Matplotlib::Simple writes it. The line names the calling script (its full
path), the function, this file and its
version, the Alien::Bioinf version, the user who ran it, the computer it ran
on (hostname and operating system), the perl version and path, and what drew
it: C<msa_plot.py>, the Python version and path, the matplotlib version, and
for a tree the Biopython and NumPy versions too. The line has no date; an SVG
has the C<< E<lt>dc:dateE<gt> >> matplotlib writes beside it, taken from
C<SOURCE_DATE_EPOCH> where that is set. After a C<;> come the image's title
(or "untitled") and what it was made from and how, with the full path and
SHA-256 of each file:

=over

=item * C<plot_msa>: the FASTA file and how many sequences it holds (or the number
of sequences in a hash ref), the ungapped copy of it that Clustal Omega
read, with the sequences named C<s1>, C<s2> and so on, the alignment Clustal
Omega wrote, with the names put back, the Clustal Omega version, the exact
command, and the C<msa_file> and C<tree_file> kept, if any.

=item * C<plot_phylo>: the same when it aligns, with the guide tree it drew;
otherwise the newick file it drew. Then how many negative branch lengths
were drawn as 0, if any; the tips, with the label each was shown as; and the
whole newick tree that was drawn.

=item * C<msa_quality_table>: the metric shown, and either the BLAST report it read
or the FASTA file, its ungapped copy, the report blastp wrote, the blastp
version, the exact command, and the C<alignment_json> kept, if any.

=back

For example:

 /home/me/work/run.pl called using "plot_phylo" in /.../Cn/Bioinf.pm
 version 0.01 with Alien::Bioinf 0.01 by user me on host myhost (linux)
 with Perl 5.44.0 (/usr/bin/perl), drawn by /.../msa_plot.py with
 Python 3.14.2 (/.../venv/bin/python),
 matplotlib 3.11.2, Biopython 1.87, NumPy 2.4.6; titled "DEG20010421";
 from the newick file /home/me/work/t.newick (SHA-256 f0e8...); 4 tips:
 S.cerevisiae, ...; the tree as drawn, in newick: (S.cerevisiae:0.389085,...);

An SVG holds it in C<< E<lt>dc:creatorE<gt> >>; C<exiftool> or C<identify -verbose> shows it
in a PNG or PDF.

=head1 FUNCTIONS

=head2 fasta2hash($file, $key)

Reads a FASTA file (gzip-compressed if its name ends in C<.gz>, which needs a
C<gzip> program on C<PATH>). Returns a hash ref of defline (without the C<< E<gt> >>) =>
sequence, or, with C<$key>, just the sequence of that defline, reading no
further than the record after it. A defline that appears twice is warned about
and its sequences concatenated; with C<$key>, only a repeat of C<$key> is looked
for. Line endings may be C<\n> or C<\r\n>. Whitespace (space, tab, C<\r>, C<\f>
and C<\x0b>) is dropped from sequences and from the end of deflines, as
Biopython drops it, so a defline of C<x> and a space is found as C<x>. Deflines and sequences are bytes,
just as they are in the file. A C<.gz> that gzip cannot read to its end, such as
a truncated one, dies rather than returning the part that was read.

 my $h = fasta2hash('DEG20010421.fa');

=head2 hash2fasta_file($hash, $filename, $order, $width)

Writes C<$hash> as FASTA: the keys in C<@$order> (default: sorted), sequences
wrapped at C<$width> columns (default 80; 0 for one line each). Returns
C<$filename>. Keys and sequences are written as bytes, as C<print> to a C<:raw>
handle writes them, so a string that perl holds as UTF-8 is written as the
characters in it. A character above 255 has no single byte and dies; encode
such a string first, as with C<Encode::encode('UTF-8', $name)>.

=head2 get_best_alignment_hit($json_file, $sort_criterion)

For a BLAST C<-outfmt 15> JSON report, a hash ref of query title => array ref
of that query's hits, best first. Each hit is its best hsp plus C<accession>,
C<hit_len>, C<id> and C<title>, and C<descriptions>, an array ref of every
description BLAST gave the hit. A hit against a non-redundant database such as
nr has one for each accession of the same sequence; C<accession>, C<id> and
C<title> are those of the first, which is the one BLAST lists first. C<$sort_criterion> is an hsp field (default
C<evalue>): C<align_len>, C<bit_score>, C<evalue>, C<gaps>, C<identity>,
C<positive> or C<score>. Ties are broken on the bit score, then on BLAST's
order. A report with two queries of one title dies, since only one of them
could be returned.

=head2 pgs_score(%args)

Computes a polygenic score for every sample of a VCF from a scoring file of
the L<PGS Catalog|https://www.pgscatalog.org/downloads/>: for each variant,
the number of copies of the effect allele in the sample's genotype times
C<effect_weight>, summed over the variants. L<PLINK 2|https://www.cog-genomics.org/plink/2.0/>,
which C<Alien::Bioinf> installs, reads the genotypes and does the arithmetic
(C<--score> with C<no-mean-imputation>), so the genotypes are never held in perl:

 use Cn::Bioinf 'pgs_score';

 my $scores = pgs_score(vcf => 'genome.vcf.gz', weights => 'PGS000138.txt.gz');
 # { sample => { score => 1.6742, average => 6.97e-05, variants => 12015,
 #               alleles => 24030, missing => 10259, mismatched => 0 } }

Both arguments may also be data frames, as Stats::LikeR's C<read_table>
returns them, which this module takes without loading Stats::LikeR; a VCF
given as a frame is written to a temporary VCF for PLINK 2. Files cost less
memory, since a frame is all in perl: on 2026-10-17, PGS000907's 1,773,528
variants were scored from files against a 3.5-million-record genome in 25 s,
with perl peaking at 0.8 GB.

=over

=item * C<vcf> (required): a C<.vcf>, C<.vcf.gz>, C<.vcf.bgz> or C<.bcf> file, a PLINK 2
fileset given by its C<.pgen> (with its C<.pvar> and C<.psam> beside it), or a
data frame in any of C<read_table>'s four shapes whose genotypes are columns
named C<< E<lt>sampleE<gt>.GT >>, as C<read_table> gives a VCF by default; one read with
C<< explode =E<gt> 0 >> has none, and dies. A frame needs C<CHROM>, C<POS>, C<REF> and
C<ALT>, and C<ID> too with C<< match =E<gt> 'id' >>.

=item * C<weights> (required): a scoring file (C<.txt> or C<.txt.gz>), read a line at a
time, or a data frame of one. It needs C<effect_allele> and C<effect_weight>,
and either C<chr_name> and C<chr_position> or, with C<< match =E<gt> 'id' >>, C<rsID>.
C<other_allele> is used when present. A harmonized file of the PGS Catalog
(C<..._hmPOS_GRCh38.txt.gz>) is read by C<hm_chr> and C<hm_pos> (or C<hm_rsID>)
instead, which are on the build its name gives. A variant marked TRUE in
C<is_haplotype>, C<is_diplotype>, C<is_interaction>, C<is_dominant> or
C<is_recessive>, or given a C<dosage_0_weight>, C<dosage_1_weight> or
C<dosage_2_weight>, dies, since only additive weights are scored.

=item * C<match>: C<position> (the default) finds each variant by chromosome and
position; C<id> finds it by rsID, in a VCF whose C<ID> column holds them. An
C<ID> of C<rs1;rs2> is found by either.

=item * C<build>: C<hg19> (or C<b37>) or C<hg38> (or C<b38>). PLINK 2 reads X only once
told where its pseudo-autosomal regions are, so a score with variants on X
or Y dies without it. Every sample is taken to be of unknown sex, which PLINK
2 counts a haploid call on X of as homozygous: 2 copies, 2 alleles.

=item * C<memory>, C<threads>: the MiB of memory (640 or more) and the CPUs PLINK 2
may use. By default it reserves half the machine's memory and uses every
CPU.

=back

Returns a hash ref keyed by sample, a hash of hashes that C<write_table> can
write, with these fields for each sample:

=over

=item * C<score>: the sum of the weights times the dosages, PLINK 2's C<SCORE1_SUM>.
PLINK 2 writes it to six significant figures, so it is good to about one
part in a million.

=item * C<average>: C<score> divided by C<alleles>, PLINK 2's C<SCORE1_AVG>, or undef
when C<alleles> is 0.

=item * C<variants>: the variants of C<weights> that were scored for this sample.

=item * C<alleles>: the alleles counted in their genotypes, PLINK 2's C<ALLELE_CT>:
2 for each diploid call.

=item * C<missing>: the variants that were not scored, because the VCF has no
record of them whose alleles agree, or because the sample's genotype is
missing (C<./.>, or partly missing, such as C<./1>). C<variants> and C<missing>
add up to the number of rows of C<weights>.

=item * C<mismatched>: of C<missing>, those at a position (or ID) the VCF has, whose
alleles disagree with the record's there.

=back

PLINK 2 by default fills a missing genotype in with the effect allele's
frequency in the samples; C<pgs_score> leaves it out.

How a variant is matched, which perl decides, a C<.pvar> line at a time,
before PLINK 2 scores:

=over

=item * Chromosomes are compared without UCSC's C<chr> prefix and in any case, and
PLINK's C<23>, C<24> and C<26> and UCSC's C<M> are taken as C<X>, C<Y> and C<MT>,
so C<chr1> in a VCF finds C<1> in a scoring file.

=item * A variant agrees with a record when its effect allele and its other allele
(if it has one) are both among the record's C<REF> and C<ALT> alleles, in any
case. A multi-allelic record counts the copies of whichever C<ALT> is the
effect allele, so a C<1/2> genotype of C<G> to C<A,T> holds one C<T>. Where two
records agree, as a multi-allelic site split in two can, the first in
C<CHROM:POS:REF:ALT> order is used.

=item * A record whose C<ALT> is C<.>, which a VCF of every site writes where the
sample matches the reference, agrees with a variant either of whose alleles
is its C<REF>; a C<0/0> there is two copies of the effect allele if it is the
C<REF>, and none otherwise.

=item * Alleles are not complemented, since a VCF and the PGS Catalog's files are
both on the forward strand. A variant written on the other strand
disagrees, and A/T and C/G variants, which read the same on either strand,
are matched as they are written.

=item * Only C<GT> is read: a dosage (C<DS>) or a genotype likelihood is not used.

=back

The positions in the scoring file have to be on the same genome build as the
VCF's, since a variant's position in GRCh37 (hg19) and in GRCh38 (hg38) differ;
rsIDs do not, so C<< match =E<gt> 'id' >> needs no build. Many scoring files say
C<genome_build=NR> (not reported), so C<pgs_score> cannot compare the builds.
The alleles say what it cannot: it warns when none of the variants was found
with alleles that agree, and when the alleles disagree at more than half of
the variants found. On 2026-10-17, PGS000138, whose header says
C<genome_build=NR>, was scored against a GRCh37 genome: 12,015 of its 22,274
variants were found, and the alleles agreed at every one, so its positions
are on GRCh37. With every position moved on one base, the alleles agreed at 9
of the 58 found.

A VCF that lists only the variants a sample carries (as most single-sample
VCFs from a sequencing service do: every genotype holds an C<ALT>) has no
record where the sample matches the reference, so those variants are
C<missing>, not scored as zero copies. PLINK 2 refuses such a file; C<pgs_score>
warns that its score is not comparable and scores it anyway. That genome above
was scored at 12,015 of PGS000138's 22,274 variants, and an exome at 701. A
score is comparable only with another scored on the same variants, so check
C<variants> before comparing samples; a VCF called jointly over many samples,
or one of every site, has far fewer missing. The reference blocks of a gVCF,
which cover a range of positions with one record, are not expanded.

PLINK 2 converts the whole of a VCF or BCF before it keeps the score's
positions, which for the 12 GB, 81-million-record BCF of 1000 Genomes took 294
s and 6.2 GB on 2026-10-17. Cut the positions out first with C<bcftools view
-R> (16 s), or convert such a reference once with C<plink2 --make-pgen> and
give C<pgs_score> the C<.pgen>: from a C<.pgen> of those positions, all 2,504
people were scored in under a second.

The PGS Catalog's search pages export a table of scores with a link to each
one's scoring file (C<depression.pgs_scores_data.csv>, for one). That table
has no weights, so C<pgs_score> cannot use it, but it can say which files to
download and score:

 use Stats::LikeR 'read_table';
 my $catalog = read_table('depression.pgs_scores_data.csv');
 foreach my $pgs (@{ $catalog }) {
     my ($file) = $pgs->{'Scoring File (FTP Link)'} =~ m{([^/]+)\z};
     next unless -f $file;
     my $s = pgs_score(vcf => 'genome.vcf.gz', weights => $file);
     printf "%s\t%s\t%g\t%d variants\n", $pgs->{'Polygenic Score ID & Name'}, $_, @{ $s->{$_} }{qw(score variants)}
         foreach sort keys %{ $s };
 }

=head2 sex_check(%args)

Tells each sample of a VCF as male or female from its genotypes on X, by the
method of L<peddy|https://github.com/brentp/peddy>'s C<sex_check> (0.4.3): of
the SNVs on X outside its pseudo-autosomal regions, the ratio of heterozygous
calls to homozygous-ALT ones, which is near 0 for a male, whose one X can be
heterozygous only by error, and above 1 for a female. The file is
decompressed by perl and its lines counted in XS, which holds nothing but
seven counters for each sample, so that a whole genome costs about what
decompressing it does. It prints a report to STDOUT and returns the same
results as a hash:

 use Cn::Bioinf 'sex_check';

 my $sex = sex_check(vcf => 'genome.vcf.gz');
 my $same = sex_check('genome.vcf.gz');    # the file may come first, unnamed

prints

 sex_check: genome.vcf.gz
 male when het/hom-ALT, over SNVs on X outside its pseudo-autosomal regions, is below 0.6; genotypes of DP below 7 not counted
 sample              sex   het/hom-ALT  X het  X hom-ALT  X haploid  X missing  Y called
 genome  male        0.029   1990      68035          0          0      2746

with one row for each sample, in the VCF's order, and a last line counting
the samples whose sex is unknown, if any. It returns

 # { sample => { sex => 'male', het_ratio => 0.0292, x_het => 1990,
 #               x_hom_alt => 68035, x_hom_ref => 0, x_haploid => 0,
 #               x_missing => 0, y_called => 2746, y_alt => 2746 } }

On 2026-10-18 a 155 MB C<.vcf.gz> of a whole genome (3.5 million records)
took 4.8 s and 47 MB, where C<zcat> alone takes 3.7 s; the same file as a
C<.vcf.bz2> took 19.4 s, nearly all of it bzip2's.

=over

=item * C<vcf> (required): a C<.vcf>, C<.vcf.gz>, C<.vcf.bgz> or C<.vcf.bz2> file, the
suffix in any case. It may also be given first, without its name, as
C<sex_check($file)> or C<< sex_check($file, report =E<gt> 0) >>. A VCF may name its chromosomes C<chrX>, C<X> or C<23>
(and C<chrY>, C<Y> or C<24>); records on any other are skipped unread.

=item * C<cutoff>: a sample whose C<het_ratio> is below this is male, and otherwise
female. The default, 0.6, is peddy's. On five genomes the three males' ratios
were 0.029 (a whole genome), 0.103 and 0.201 (exomes), and the two
females' 1.660 and 1.712 (exomes); an exome's off-target calls make more
false heterozygotes than a genome's.

=item * C<min_depth>: a genotype whose C<DP> is below this is not counted. The
default, 7, is peddy's; 0 counts every genotype. A genotype with no C<DP>,
or a C<DP> of C<.>, is counted.

=item * C<report>: 1 (the default) prints the report; 0 prints nothing, and only
returns the results.

=back

What is counted, and where C<sex_check> differs from peddy:

=over

=item * X's pseudo-autosomal regions are taken as X:10000-2781479 and
X:154931044-156030895, peddy's bounds, which cover GRCh37's and GRCh38's,
so no build need be given. Y's are both builds' too.

=item * Only biallelic SNVs are counted: a record whose C<REF> or C<ALT> is not a
single base (an indel, a multi-allelic site, or an C<ALT> of C<.>, C<*> or
C<< E<lt>NON_REFE<gt> >>) is skipped, where peddy keeps an insertion. A record with no
C<GT> in its C<FORMAT> is skipped too.

=item * Every site is counted, where peddy stops at 10,000, and peddy's filters on
a site's call rate and ALT frequency across the cohort, which mean nothing
for one genome, are left out.

=item * A haploid call, as a caller that knows a sample is male writes it on X,
counts as homozygous: C<1> as hom-ALT, C<0> as hom-REF.

=item * A sample with no het or hom-ALT call on X has C<sex> and C<het_ratio> undef.
peddy divides by a hom-ALT count of 0 as if it were 1, and so would call
such a sample male. When the whole VCF has no SNV on X outside the
pseudo-autosomal regions, C<sex_check> also warns.

=back

Returns a hash ref keyed by sample, a hash of hashes that C<write_table> can
write, with these fields for each sample:

=over

=item * C<sex>: C<male>, C<female>, or undef.

=item * C<het_ratio>: C<x_het> divided by C<x_hom_alt> (by 1 when C<x_hom_alt> is 0),
or undef when both are 0.

=item * C<x_het>, C<x_hom_alt>, C<x_hom_ref>: the sample's heterozygous,
homozygous-ALT and homozygous-REF calls at the SNVs counted on X. A VCF
that lists only the sites where its sample carries an ALT allele has no
hom-REF calls.

=item * C<x_haploid>: of those, the haploid calls.

=item * C<x_missing>: the sample's calls there that were missing (C<./.>, C<.>, or
partly missing, such as C<./1>), or not a call of one or two alleles of the
record's two. A genotype left out by C<min_depth> is in no count.

=item * C<y_called>, C<y_alt>: the sample's calls at SNVs on Y outside its
pseudo-autosomal regions, and those of them that carry the C<ALT>. They
decide nothing: an exome of a woman had 530 called SNVs on Y, and one of a
man 659, since a short read from X often maps to Y. They show at most
whether Y was called at all.

=back

A BCF or a C<.pgen> is not read; convert it with C<bcftools view -Oz>.

=head2 parent(%args)

Tells whether two people are parent and child from their genotypes, and
returns 1 if they are and 0 if not. A child has one allele from each parent at
every site, so wherever a parent is homozygous the child carries that allele,
and the two are never opposite homozygotes (one hom-REF and the other
hom-ALT, which is called IBS0) except by a genotyping error. Genotypes cannot
tell which of the two is the parent, so C<child> and C<parent> may be given
either way round, with the same result. Paternity and maternity are this
together with the parent's sex, which C<sex_check> tells.

Over the autosomal biallelic SNVs where either of the two carries the ALT, a
pair is parent and child when the IBS0 rate (the share of those sites that are
IBS0) is below C<max_ibs0>, and the
L<KING|https://www.kingrelatedness.com/manual.shtml>-robust kinship
coefficient is from 0.177 to 0.354, KING's range for first-degree relatives.
Full siblings have that kinship too, but they share no allele by descent at a
quarter of sites, so they have opposite homozygotes there, and the IBS0 rate
tells them apart. One person sequenced twice, or identical twins, have a
kinship above 0.354, and are not parent and child.

The file is decompressed by perl and counted in XS, as for C<sex_check>. The
two people are either two samples of one VCF, or each the one sample of a VCF
of their own, as a sequencing company delivers it:

 use Cn::Bioinf 'parent';

 my $is = parent(child => 'child.vcf.gz', parent => 'father.vcf.gz');
 my ($same, $stats) = parent(vcf => 'family.vcf.gz', child => 'C', parent => 'F');

The first prints

 parent: child child.vcf.gz, parent father.vcf.gz
 parent and child (1): the IBS0 rate, 0.0030, is below 0.012, and kinship 0.2587 is from 0.177 to 0.354, as for first-degree relatives
 sites compared           150970  autosomal biallelic SNVs where either carries the ALT, DP 10 or more
   in both VCFs            98828
   in the child's alone    27645  the parent taken as hom-REF
   in the parent's alone   24497  the child taken as hom-REF
 IBS0                        458  opposite homozygotes
   child hom-ALT             305  and the parent hom-REF
   child hom-REF             153  and the parent hom-ALT
 IBS1                      76262  one allele shared
 IBS2                      74250  both alleles shared
 het in both               43727
 het                       82821  the child's, and 80895 the parent's
 hom-ALT                   43652  the child's, and 42430 the parent's
 kinship                  0.2587  KING-robust, between-family

and returns 1. On 2026-10-18 the exomes of one family, all three sequenced and
called together, gave:

 child, parent    sites compared  IBS0 rate  kinship  returned
 child, father           150,970     0.0030   0.2587         1
 child, mother           158,616     0.0020   0.2604         1
 father, mother          171,033     0.0810   0.0174         0

Each pair took 1.6 s. A genome against an exome took 7.3 s and 89 MB, and is
not to be trusted, for the reason under "Each sample's own VCF" below.

=over

=item * C<child>, C<parent> (required): without C<vcf>, the two VCF files, each of one
sample (C<.vcf>, C<.vcf.gz>, C<.vcf.bgz> or C<.vcf.bz2>, as for C<sex_check>);
with C<vcf>, the names of two of its samples.

=item * C<vcf>: a VCF that holds both samples.

=item * C<min_depth>: a genotype whose C<DP> is below this is not counted. The
default is 10; 0 counts every genotype. A genotype with no C<DP>, or a C<DP>
of C<.>, is counted. On the family above, counting every genotype gave the
parent and child pairs IBS0 rates of 0.254 and 0.215, as high as an
unrelated pair's, from calls of C<DP> 1 to 3 made off the exome's targets. A
C<min_depth> of 7 gave 0.0053 and 0.0037, and 10 gave 0.0030 and 0.0020,
keeping 92% of the sites of 7.

=item * C<max_ibs0>: the IBS0 rate a parent and child must be below. The default,
0.012, is L<peddy|https://github.com/brentp/peddy> 0.4.3's (C<ped_check>).
It lies between the parent and child rates above (0.002 to 0.005) and the
0.02 or so expected of full siblings, a quarter of the unrelated pair's 0.081.

=item * C<report>: 1 (the default) prints the report to STDOUT; 0 prints nothing.

=back

What is counted:

=over

=item * Only the autosomes, named C<1> to C<22> or C<chr1> to C<chr22>; X is left out,
since a father gives his X to a daughter and not to a son.

=item * Only biallelic SNVs: a record whose C<REF> or C<ALT> is not a single base, or
that has no C<GT>, is skipped.

=item * A call counts when it is diploid and of the record's C<REF> and C<ALT>, and
is not below C<min_depth>. A missing or partly missing call (C<./.>, C<./1>), a
haploid one, or one with an allele the site lacks leaves that site out.
So do sites where both are hom-REF, which say nothing about relatedness; a
cohort's VCF has a great many of them.

=back

Each sample's own VCF lists only the sites where its sample carries an ALT
allele, so a site that one file lacks is taken as hom-REF in that person, and
the two files are compared site by site, holding both in memory (8 bytes a
site, so about 30 MB for a genome). A position is left out when either file has more
than one SNV there (as after C<bcftools norm -m->), a call that does not count,
or a record of more than one C<ALT>, or when the two files give it different
C<REF> or C<ALT> alleles. The files may name and order their chromosomes
differently.

Taking a missing site as hom-REF is right only when the other file was
sequenced and called there too, so two VCFs from one pipeline are what this
compares well. A genome's VCF has millions of sites that an exome's lacks, and
two exome capture kits target different regions; either makes IBS0 out of
sites that one of the two never had sequenced. The report gives the sites in
each file alone for that reason, and C<parent> warns when one file has more
than twice as many as the other: the three pairs of the family above had
ratios of 1.13 at most, an exome of another kit had 3.0 to 13 against them,
and a genome against an exome 58 (with a kinship of -19.8). A VCF
that leaves out filtered calls has holes that look hom-REF too.

C<parent> dies, rather than return a 0 or 1 that it has no grounds for, when

=over

=item * fewer than 1,000 sites are compared, which at the unrelated pair's IBS0 rate
would hold 81 opposite homozygotes on average, 8 standard deviations above
the 12 that C<max_ibs0> allows (an exome gives about 150,000 and a genome 3
million);

=item * neither has a hom-REF call where the other carries the ALT, so that IBS0
cannot be seen. A VCF merged from VCFs that list only their own sample's
variants, as C<bcftools merge> writes it without C<-0>, is like that, with
C<./.> wherever a sample had no variant; give each sample's own VCF instead;

=item * either has no het call at the sites compared, so that no kinship can be
taken.

=back

In scalar context C<parent> returns 1 or 0. In list context it returns that and
a hash ref of what it was decided from, with the same fields on every return:

=over

=item * C<parent>: 1 or 0, as returned.

=item * C<relationship>: C<parent-child>, C<full siblings>, C<duplicate> (one person or
identical twins), C<second-degree>, C<third-degree> or C<unrelated>, by
KING's kinship ranges, with the IBS0 rate telling parent and child from
full siblings.

=item * C<kinship>: the KING-robust between-family estimate (Manichaikul et al. 2010,
I<Bioinformatics> 26:2867), which is what PLINK 2's C<--make-king-table>
gives, and which agreed with it exactly on the VCF of two of the exomes
above merged by C<bcftools merge -0>.

=item * C<ibs0_rate>: C<ibs0> divided by C<sites>.

=item * C<sites>: the sites compared, and of them C<in_both>, C<child_only> and
C<parent_only>: those in both own VCFs, in the child's alone and in the
parent's alone. With C<vcf>, every site is C<in_both>.

=item * C<ibs0>, C<ibs1>, C<ibs2>: the sites where the two share no allele, one, or
both; and C<ibs0> split into C<child_alt_parent_ref> and
C<child_ref_parent_alt>.

=item * C<het_het>: the sites where both are het; C<child_het>, C<parent_het>,
C<child_hom_alt> and C<parent_hom_alt>: each one's het and hom-ALT calls at
the sites compared.

=back

A BCF is not read; convert it with C<bcftools view -Oz>.

=head2 plot_msa(%args)

Aligns sequences with Clustal Omega and draws the alignment. Returns a hash ref
of the files made: C<filename>, C<msa_file>, and C<tree_file> if it was given.
Once the image is written it prints C<wrote> and its file name to STDOUT, the
name in black on yellow when STDOUT is a terminal.

 plot_msa(
     fasta      => 't/data/DEG20010421.fa',
     filename   => 'DEG20010421.msa.svg',
 );

=over

=item * C<fasta>, C<filename> (required): the sequences, as a FASTA file name or a
hash ref of name => sequence (the function tells the two apart by whether it
is a reference); and the image to draw, whose extension picks the format.
A name may hold spaces and punctuation: Clustal Omega is given the sequences
as C<s1>, C<s2> and so on, and the alignment and tree it writes get their
names back, in UTF-8. In the tree, a name with whitespace or any of
C<()[]':;,> in it is quoted as newick quotes one, as C<'it''s (x:1)'>.

=item * C<msa_file>, C<tree_file>: where to keep clustalo's alignment (FASTA; a
temporary file otherwise) and its guide tree (newick; not made otherwise).
Keep the tree to draw it with C<plot_phylo> without aligning a second time.

=item * C<order>: names, first to last (default: the input order; for a hash,
sorted). Only these are drawn, and the first is drawn at the bottom.

=item * C<labels>: a hash ref of name => label to show instead. matplotlib mathtext
works: C<< 'C.albicans' =E<gt> '$\it{C. albicans}$' >>. Two sequences drawn under one
label die, since they would be one row of the image.

=item * C<active_site_aa>, C<query>: C<< { His395 =E<gt> 395, ... } >>, a dashed vertical line
at each of these 1-based residue numbers of the sequence named by C<query>.
A number below 1 dies.

=item * C<title>, C<xlabel>, C<ylabel>, C<threads>, C<clustal_args>: the plot title; the
axis labels (default "Amino Acid Residue" and "Protein & Species");
clustalo threads (default 1); and an array ref of extra clustalo arguments.

=back

With fewer than two sequences it warns and returns an empty hash ref.

=head2 plot_phylo(%args)

Draws a guide tree, with Biopython's C<Bio.Phylo>, from a FASTA (aligned with
Clustal Omega first) or from a newick file that C<plot_msa> kept. Returns a
hash ref of the files made or used: C<output_file>, and C<tree_file> and
C<msa_file> as below.

 plot_phylo(fasta => 't/data/DEG20010421.fa');   # writes phylo.svg
 plot_phylo(
     'tree_file'   => 'DEG20010421.newick',
     'output_file' => 'DEG20010421.tree.png',
 );

=over

=item * C<output_file>: the image to draw (default C<phylo.svg>, in the working
directory); the extension picks the format.

=item * C<fasta>, C<tree_file>: with C<fasta> (as for C<plot_msa>), the sequences are
aligned with Clustal Omega and the guide tree drawn; C<tree_file> and
C<msa_file> then say where to keep the tree and alignment, and C<threads> and
C<clustal_args> are as for C<plot_msa>. Without C<fasta>, C<tree_file> is an
existing newick file to draw, such as one C<plot_msa> kept, and no alignment
is made.

=item * C<labels>, C<title>: a hash ref of name => label for the tips, as for
C<plot_msa>; and the plot title.

=back

With C<fasta> of fewer than two sequences it warns and returns an empty hash
ref.

=head2 msa_quality_table(%args)

Draws an all-against-all BLAST score table with matplotlib, and returns
C<filename>. The simplest call is

 msa_quality_table(
     fasta      => 't/data/DEG20010421.fa',  # or { name => sequence }, aligned or not
     filename   => 'DEG20010421.scores.png',
 );

=over

=item * C<fasta>, C<filename> (required): the sequences, as a FASTA file name or a
hash ref of name => sequence, as for C<plot_msa>; they are aligned all
against all with C<blastp>. Gaps are stripped first, so an aligned FASTA,
such as the one C<plot_msa> keeps, will do. Names must look like
C<Genus.species[.strain]>. C<fasta> is not needed when an existing
C<alignment_json> is given. C<unaligned_fa> is its old name.

=item * C<alignment_json>: the C<blastp -outfmt 15> report, as a parsed hash ref or a
file name. An existing file is read, and nothing is aligned; otherwise
C<blastp> is run on C<fasta> and its report kept there for next time. Without
it, the report is a temporary file.

=item * C<metric>: the hsp field to show (default C<score>).

=item * C<normalize>: divide every value by the largest, after adding
C<logscale_add> to both, so the scale runs to 1. A largest value of 0 or
less, as C<evalue> has when BLAST has rounded every e-value to 0, dies.

=item * C<order>, C<logscale_add>, C<default_undefined>, C<title>, C<cb_label>, C<cb_min>,
C<cb_max>, C<cblogscale>, C<show_numbers>: the sequences to show, in order; a
number added to every value; the value of a pair with no hit (otherwise
drawn grey); the title; the colour bar's label, lower and upper ends, and
whether it is logarithmic; and whether each cell shows its number.

=back

C<msa_file> is accepted and ignored, for old callers.

=head2 clustal_view_residues(%args)

Writes an alignment as LaTeX tables with chosen residues coloured, and returns
C<output_tex_file>, printing C<wrote> and that file name to STDOUT, the name in
black on cyan when STDOUT is a terminal. Protein names are written so that
LaTeX prints them as they are, C<_>, C<^>, C<{> and the like included.

=over

=item * C<msa_file> (required): a FASTA file. If its sequences are all one length
(such as the alignment C<plot_msa> keeps) it is shown as it is; if not, it is
first aligned with Clustal Omega into a temporary file, and C<msa_file>
itself is never written.

=item * C<output_tex_file> (required): the LaTeX file to write, meant to be
C<\input> into a document whose preamble loads C<xcolor> (for the colours)
and C<placeins> (for the C<\FloatBarrier> after each table).

=item * C<color_residues>: C<< { protein =E<gt> { residue number =E<gt> colour } } >>, where
residue numbers are 1-based and a colour is an xcolor name or
C<[r, g, b]>; a coloured column is coloured in every protein.

=item * C<track>: a protein that gets a row under it showing its coloured residue
numbers.

=item * C<order>: an array ref of the proteins to show, top to bottom (default:
sorted, ignoring case).

=item * C<row_width>: alignment columns per block (default 100).

=item * C<split>: blocks per LaTeX table (default 4); further tables are captioned
"(continued)".

=item * C<caption>: the table caption (default empty).

=item * C<label>: written as C<\label{tab:label}>, or as C<tab:label0>, C<tab:label1>,
... when there is more than one table.

=item * C<table_text_size>: the LaTeX size command put at the start of each table
(default C<\footnotesize>).

=item * C<threads>, C<clustal_args>: clustalo threads (default 1) and an array ref of
further clustalo arguments, when C<msa_file> has to be aligned.

=back

=head1 Thanks

A lot of this work (not all!) used Claude AI, which was paid for by the University of Idaho's IMCI
