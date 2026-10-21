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
#
# concordance is ~/Scripts/old/x.low_ram_local_concordance.pl (as it stood on
# 2026-10-19) made a function of two VCFs: its 45 kinds of comparison and its
# metrics are counted in Bioinf.xs, whose comments above CONC_NAME and CONC_TP
# say where it departs from that script; its XLSX, JSON and ACMG59 output are
# not made.
require 5.012;
use strict;
use warnings;
use Carp qw(croak carp);
use File::Spec;
use Sys::Hostname 'hostname';
use Exporter 'import';
use XSLoader;
our $VERSION = '0.01';
our @EXPORT_OK = qw(clustal_view_residues concordance fasta2hash get_best_alignment_hit hash2fasta_file hirisplex_genotypes msa_quality_table parent pdb_search pgs_score plot_msa plot_phylo relatedness sex_check trait_predict);
our @EXPORT = @EXPORT_OK;
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
	my ($file, $bg, $verb) = @_;
	$verb //= 'wrote'; # or 'will write', for an image whose script has not run yet
	print STDOUT -t STDOUT ? "$verb \e[30;${bg}m$file\e[0m\n" : "$verb $file\n";
}

# A new temporary file's name, with $suffix. It is removed when perl exits
# unless $keep is true.
sub _tmp {
	my ($suffix, $keep) = @_;
	require File::Temp;
	my ($fh, $name) = File::Temp::tempfile(SUFFIX => $suffix, UNLINK => $keep ? 0 : 1, TMPDIR => 1);
	close $fh;
	$name;
}

# "fh" and "execute" of $sub's $args, resolved and checked before anything
# runs, as Matplotlib::Simple's plt takes them: "fh" is a File::Temp that the
# Python drawing the image is appended to, so that one script, run once, can
# draw several images (and those of Matplotlib::Simple, which takes the same
# File::Temp); "execute" says whether to run that script now.
sub _drawing {
	my ($sub, $args) = @_;
	my %r = (fh => $args->{fh});
	$r{execute} = $args->{execute} // 1; # 1 = run the script now; 0 = only write it
	croak "$sub: \"execute\" must be 0 or 1, not " . _what($r{execute})
		unless !ref $r{execute} && ($r{execute} eq '0' || $r{execute} eq '1');
	if (defined $r{fh}) {
		require Scalar::Util;
		croak "$sub: \"fh\" must be a File::Temp object, as for Matplotlib::Simple, not " . _what($r{fh})
			unless Scalar::Util::blessed($r{fh}) && $r{fh}->isa('File::Temp');
	}
	\%r;
}

# share/msa_plot.py run with @args for $sub, as Python appended to the script
# that $drawing (from _drawing) names, or to a new temporary one; the script is
# run under Alien::Bioinf's Python if $drawing->{execute}. Returns the script's
# file name. A new script is kept after perl exits only when it is not run
# here, since only then is it of use afterwards.
#
# The script is the whole of what it needs, so that it can be run later, after
# perl has removed its temporary files: each file passed with --f, --tree or
# --table is carried in it, and written out again, to a temporary directory,
# as the script runs. Nothing perl was given is pasted into the Python source:
# the arguments, the files and msa_plot.py's path are one JSON object, base64
# encoded, as Matplotlib::Simple's write_data passes data, since a Windows
# path's "\b" is a backspace in a Python literal and an apostrophe would close
# it. msa_plot.py is run by runpy in a namespace of its own, so that its
# globals meet neither those of another image in the script nor those of
# Matplotlib::Simple's code. -I keeps PYTHONPATH, PYTHONHOME and the user's own
# site-packages out, as for _python.
sub _draw {
	my ($sub, $drawing, @args) = @_;
	require JSON::MaybeXS;
	require MIME::Base64;
	my %file = map { $_ => 1 } qw(--f --tree --table); # the options whose value is a file msa_plot.py reads
	my @json;
	for (my $i = 0; $i < @args; $i++) {
		push @json, _chars($args[$i]);
		next unless $file{ $args[$i] } && $i < $#args;
		my $in = $args[++$i];
		open my $fh, '<:raw', $in or croak "_draw: can't read $in: $!";
		my $bytes = do { local $/; <$fh> };
		close $fh or croak "_draw: can't read $in: $!";
		# named by its place, so that no two collide, and with its own suffix
		push @json, { name => 'input' . scalar(@json) . ($in =~ /(\.[A-Za-z]+)\z/ ? $1 : ''), b64 => MIME::Base64::encode_base64($bytes // '', '') };
	}
	my $b64 = MIME::Base64::encode_base64(JSON::MaybeXS::encode_json({ script => _chars(_script()), args => \@json }), '');
	my $py = <<"PYTHON";
# $sub, from Cn::Bioinf $VERSION: msa_plot.py run on the arguments below
import base64, json, os, runpy, sys, tempfile
cn_bioinf = json.loads(base64.b64decode('$b64').decode('utf-8'))
with tempfile.TemporaryDirectory() as cn_bioinf_dir:
	cn_bioinf_argv = []
	for cn_bioinf_arg in cn_bioinf['args']:
		if isinstance(cn_bioinf_arg, dict):
			cn_bioinf_path = os.path.join(cn_bioinf_dir, cn_bioinf_arg['name'])
			with open(cn_bioinf_path, 'wb') as cn_bioinf_fh:
				cn_bioinf_fh.write(base64.b64decode(cn_bioinf_arg['b64']))
			cn_bioinf_arg = cn_bioinf_path
		cn_bioinf_argv.append(cn_bioinf_arg)
	sys.argv = [cn_bioinf['script']] + cn_bioinf_argv
	runpy.run_path(cn_bioinf['script'], run_name='__main__')
PYTHON
	my $script;
	if (defined $drawing->{fh}) {
		my $fh = $drawing->{fh};
		$script = $fh->filename;
		print {$fh} $py or croak "$sub: can't write to \"fh\", $script: $!";
		$fh->flush or croak "$sub: can't write to \"fh\", $script: $!";
	} else {
		$script = _tmp('.py', !$drawing->{execute});
		open my $fh, '>:raw', $script or croak "$sub: can't write $script: $!";
		print {$fh} $py or croak "$sub: can't write $script: $!";
		close $fh or croak "$sub: can't write $script: $!";
	}
	_run(_alien($sub)->python, '-I', $script) if $drawing->{execute};
	$script;
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

# The rows of the option $what ("weights" for pgs_score, "model" for
# trait_predict) one at a time, each a hash ref with how a message names it,
# and the column names: from a data frame (see _frame_rows), or read from a
# tab-separated file a line at a time, so that a score of millions of variants
# is never held whole. An empty field of a file is '', where read_table gives
# undef; every check below treats the two alike.
sub _weight_rows {
	my ($sub, $weights, $what) = @_;
	$what //= 'weights';
	if (ref $weights) {
		my ($rows, $labels) = _frame_rows($sub, $what, $weights, 1);
		my $i = 0;
		return (sub { $i < @{ $rows } ? ($rows->[$i], $labels->[$i++]) : () }, [keys %{ $rows->[0] }]);
	}
	croak "$sub: \"$what\" file \"$weights\" doesn't exist or isn't readable" unless -f $weights && -r _;
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
	croak "$sub: \"$what\" file $weights has no header line" unless defined $head;
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
# a whole genome, genome.vcf.gz (155 MB, 3.5 million records), in 4.84 s with
# blocks of 16 KiB, 64 KiB or 256 KiB, and 4.90 s and 4.86 s with 1 MiB and 4
# MiB (best of three, peak RSS 47 MB each), so 64 KiB, as for FASTA. The time
# is gzip's: the same file uncompressed (928 MB) took 0.26 s and 11 MB, and
# zcat alone 3.7 s; as .vcf.bz2 it took 19.4 s, and bzip2 -dc alone 18.5 s.
# IO::Uncompress::Gunzip's readline, which is perl, took 14.3 s just to read
# the lines, and the read loop alone in 1 MiB blocks 5.7 s.
my $Vcf_block = 1 << 16;

# Dies unless the VCF $file, given to $sub as its option $key, is a readable
# file with a VCF's name. Returns how it is compressed: 'gz', 'bz2', or '' for
# not at all. Its messages are prefixed with $sub, whose argument the file is.
sub _vcf_check {
	my ($sub, $key, $file) = @_;
	croak "$sub: \"$key\" must be a file name, not " . _what($file) if ref $file || !defined $file;
	croak "$sub: \"$key\" file \"$file\" doesn't exist or isn't readable" unless -f $file && -r _;
	return 'gz' if $file =~ /\.vcf\.b?gz\z/i;
	return 'bz2' if $file =~ /\.vcf\.bz2\z/i;
	return '' if $file =~ /\.vcf\z/i;
	croak "$sub: \"$key\" file $file must end in .vcf, .vcf.gz, .vcf.bgz or .vcf.bz2";
}

# The VCF $file, given to $sub as its option $key, opened: returns a sub that
# gives the next block of it and then undef, and a sub that closes it. Its
# messages are prefixed with $sub, whose argument the file is.
sub _vcf_open {
	my ($sub, $key, $file) = @_;
	my $packed = _vcf_check($sub, $key, $file);
	my ($fh, $error); # $error: what says why a read failed
	# Transparent => 0, or IO::Uncompress would pass a file that isn't compressed
	# through as it is, and the error would be that it isn't a VCF. For a file too
	# short to hold a header, new() gives undef and leaves the error empty in
	# IO::Uncompress::Gunzip 2.223 (perl 5.44.0), though 2.020 (perl 5.10.1) says
	# "Header Error: Minimum header size is 10 bytes".
	if ($packed eq 'gz') {
		require IO::Uncompress::Gunzip;
		$fh = IO::Uncompress::Gunzip->new($file, MultiStream => 1, Transparent => 0)
			or croak "$sub: can't read $file: " . ($IO::Uncompress::Gunzip::GunzipError || 'it isn\'t gzip-compressed');
		$error = sub { $IO::Uncompress::Gunzip::GunzipError };
	} elsif ($packed eq 'bz2') {
		require IO::Uncompress::Bunzip2;
		$fh = IO::Uncompress::Bunzip2->new($file, MultiStream => 1, Transparent => 0)
			or croak "$sub: can't read $file: " . ($IO::Uncompress::Bunzip2::Bunzip2Error || 'it isn\'t bzip2-compressed');
		$error = sub { $IO::Uncompress::Bunzip2::Bunzip2Error };
	} else {
		open $fh, '<:raw', $file or croak "$sub: can't read $file: $!";
		$error = sub { "$!" };
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

# The header of the VCF $file, given to $sub as its option $key, read as far as
# its first line that does not start with "##": a hash ref whose "contig" is an
# array ref of the ##contig lines' [ID, length, assembly], in the header's
# order, and whose "reference" is the first ##reference line's value, or
# undef. An ID loses its "chr" prefix, in any case, as Bioinf.xs's
# conc_chrom() drops it from CHROM; a length or an assembly that the line
# lacks is undef. Whatever is wrong with the file beyond its ## lines is left
# for Bioinf.xs to find, as it reads it.
sub _vcf_meta {
	my ($sub, $key, $file) = @_;
	my ($next, $close) = _vcf_open($sub, $key, $file);
	my %meta = (contig => [], reference => undef);
	my $carry = ''; # the start of a line that the last block ended inside
	my $done = 0;   # 1 = a line not starting with "##" has been read
	while (!$done && defined(my $block = $next->())) {
		$carry .= $block;
		my $end = rindex $carry, "\n";
		next if $end < 0;
		foreach my $line (split /\n/, substr($carry, 0, $end)) {
			$line =~ s/\r\z//;
			if (substr($line, 0, 2) ne '##') {
				$done = 1;
				last;
			}
			if ($line =~ /\A##contig=<(.*)>\z/) {
				my $attr = $1;
				my ($id) = $attr =~ /(?:\A|,)ID=([^,]+)/;
				next unless defined $id;
				$id =~ s/\Achr(?=.)//i;
				my ($length) = $attr =~ /(?:\A|,)length=([0-9]+)/;
				my ($assembly) = $attr =~ /(?:\A|,)assembly="?([^,"]+)/;
				push @{ $meta{contig} }, [$id, $length, $assembly];
			} elsif ($line =~ /\A##reference=(.+)\z/) {
				$meta{reference} //= $1;
			}
		}
		$carry = substr($carry, $end + 1);
	}
	$close->();
	\%meta;
}

# Dies unless the VCFs $test and $truth, given to $sub, can be on one genome
# build: every contig that both headers give a length for must have the same
# length in each. A chromosome's length differs from one build to the next, and
# so does a site's POS, so that two builds share few sites, and only by chance.
# A header with no lengths cannot be checked, and is not. Contigs are matched by
# name less "chr", as sites are, so hg19 and b37 (chr1 and 1, both 249,250,621
# bases by the headers of ~/Documents/Genomics/3025181.vcf.gz and of 1000
# Genomes' hs37d5 BCF, read 2026-10-20) pass.
sub _same_build {
	my ($sub, $test, $truth) = @_;
	my %meta = (test => _vcf_meta($sub, 'test', $test), truth => _vcf_meta($sub, 'truth', $truth));
	my %truth_length; # contig => its length in the truth, for those that have one
	foreach my $c (@{ $meta{truth}{contig} }) {
		$truth_length{ $c->[0] } //= $c->[1] if defined $c->[1];
	}
	my ($shared, @differ) = (0); # @differ: [contig, test length, truth length], in the test's order
	my %seen;
	foreach my $c (@{ $meta{test}{contig} }) {
		my ($id, $length) = @{ $c };
		next if !defined $length || !defined $truth_length{$id} || $seen{$id}++;
		$shared++;
		push @differ, [$id, $length, $truth_length{$id}] if $length != $truth_length{$id};
	}
	return unless @differ;
	# what each header says its build is: the assembly of a ##contig line that
	# differs, or else its ##reference
	my %says;
	foreach my $key (qw(test truth)) {
		my ($asm) = map { $_->[2] } grep { defined $_->[2] && defined $_->[1] && $_->[0] eq $differ[0][0] } @{ $meta{$key}{contig} };
		$says{$key} = $asm // $meta{$key}{reference};
	}
	my $listed = 3; # contigs named in the message; the rest are counted
	my @named = map { "$_->[0] as $_->[1] bases long in the test and $_->[2] in the truth" } @differ[0 .. ($#differ < $listed - 1 ? $#differ : $listed - 1)];
	croak "$sub: \"test\" $test and \"truth\" $truth are on different genome builds: their ##contig lines give " . join(', ', @named)
		. (@differ > $listed ? ', and ' . (@differ - $listed) . ' more' : '') . ' (' . scalar(@differ) . " of the $shared contigs that both give a length for differ)"
		. (defined $says{test} && defined $says{truth} ? "; the test's header says \"$says{test}\" and the truth's \"$says{truth}\"" : '')
		. ', so a POS in one is not the same place in the other. Lift one over to the other\'s build first';
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

# peddy 0.4.3's (sex_check() in peddy/peddy.py): male below a het to hom-ALT
# ratio of 0.6 on non-PAR X. sex_check's default "cutoff", and what relatedness
# tells sex by.
my $Sex_cutoff = 0.6;

# The het_ratio and sex of one sample's counts, $s, from Bioinf.xs's
# sex_count(), set in $s: male when the ratio is below $cutoff. peddy divides
# by a hom-ALT count of 0 as if it were 1, and so calls a sample with no calls
# on X male; here it is undef, as is its sex.
sub _sex_of {
	my ($s, $cutoff) = @_;
	my $calls = $s->{x_het} + $s->{x_hom_alt};
	$s->{het_ratio} = $calls ? $s->{x_het} / ($s->{x_hom_alt} || 1) : undef;
	$s->{sex} = !$calls ? undef : $s->{het_ratio} < $cutoff ? 'male' : 'female';
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
	$r{cutoff} = $args->{cutoff} // $Sex_cutoff;
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
	_sex_of($_, $r{cutoff}) foreach values %{ $tally };
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

# The KING-robust between-family kinship (Manichaikul et al. 2010,
# Bioinformatics 26:2867) of two samples with $het_a and $het_b het calls at the
# sites compared, $het_het sites where both are het and $ibs0 opposite
# homozygotes, in the form PLINK 2's --make-king-table gives it (t/parent.t's
# header says how the two were found to agree); undef when either has no het
# call.
sub _king {
	my ($het_a, $het_b, $het_het, $ibs0) = @_;
	my $min_het = $het_a < $het_b ? $het_a : $het_b;
	return undef unless $min_het;
	0.5 - ($het_a + $het_b - 2 * $het_het + 4 * $ibs0) / (4 * $min_het);
}

# The relationship that a $kinship and an IBS0 rate say, by KING's ranges, with
# an IBS0 rate below $max_ibs0 telling parent and child from full siblings.
sub _relationship {
	my ($kinship, $ibs0_rate, $max_ibs0) = @_;
	my $band = (grep { $kinship > $King[$_][0] } 0 .. $#King)[0];
	!defined $band ? 'unrelated' : $band == 0 ? 'duplicate'
		: $band == 1 ? ($ibs0_rate < $max_ibs0 ? 'parent-child' : 'full siblings')
		: $band == 2 ? 'second-degree' : 'third-degree';
}

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
	$s{kinship} = _king(@s{qw(child_het parent_het het_het ibs0)});
	croak "$sub: " . ($s{child_het} ? 'the parent' : 'the child') . " of $what has no het call at the sites compared, so no kinship can be taken" unless defined $s{kinship};
	$s{ibs0_rate} = $s{ibs0} / $s{sites};
	$s{relationship} = _relationship($s{kinship}, $s{ibs0_rate}, $r{max_ibs0});
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

# The words for a male and a female relative, in relatedness's categories:
# the parent, the child, the sibling.
my %Kin_word = (male => [qw(father son brother)], female => [qw(mother daughter sister)]);

# relatedness's report, given its files, the resolved options, the samples in
# order, each one's sex counts, and $each_pair, which gives each pair in order
# to a sub. Columns are as wide as their widest entry, numbers right-aligned.
# The pairs' rows are made twice, once for the widths and once for the lines,
# rather than held: 1,000 samples have 499,500.
sub _relatedness_report {
	my ($files, $r, $names, $tally, $each_pair) = @_;
	# the lines of the rows that $rows gives to a sub, under @head, with the
	# columns in %$right right-aligned
	my $table = sub {
		my ($right, $rows, @head) = @_;
		my @w = map { length } @head;
		$rows->(sub {
			foreach my $i (0 .. $#_) {
				$w[$i] = length $_[$i] if length $_[$i] > $w[$i];
			}
		});
		my $text = '';
		my $line = sub {
			(my $l = join '  ', map { $right->{$_} ? sprintf('%*s', $w[$_], $_[$_]) : sprintf('%-*s', $w[$_], $_[$_]) } 0 .. $#_) =~ s/ +\z//;
			$text .= "$l\n";
		};
		$line->(@head);
		$rows->($line);
		$text;
	};
	my $several = @{ $files } > 1;
	my $samples = sub {
		my ($row) = @_;
		foreach my $x (@{ $names }) {
			my $s = $tally->{$x};
			$row->($x, $s->{sex} // 'unknown', defined $s->{het_ratio} ? sprintf('%.3f', $s->{het_ratio}) : '-', $several ? $s->{file} : ());
		}
	};
	my $pairs = sub {
		my ($row) = @_;
		$each_pair->(sub {
			my ($x, $y, $s) = @_;
			$row->($x, $y, $s->{sites}, defined $s->{ibs0_rate} ? sprintf('%.4f', $s->{ibs0_rate}) : '-', defined $s->{kinship} ? sprintf('%.4f', $s->{kinship}) : '-',
				defined $s->{undecided} ? "undecided: $s->{undecided}" : $s->{category} . (defined $s->{parent} ? " ($s->{parent} the parent)" : ''));
		});
	};
	'relatedness: ' . join(', ', @{ $files }) . "\n"
		. "KING-robust kinship over autosomal biallelic SNVs where either carries the ALT, genotypes of DP below $r->{min_depth} not counted;"
		. " first-degree from 0.177 to 0.354, and of those, parent and child when the IBS0 rate is below $r->{max_ibs0};"
		. " male when het/hom-ALT on X outside its pseudo-autosomal regions is below $Sex_cutoff\n"
		. $table->({ 2 => 1 }, $samples, 'sample', 'sex', 'X het/hom-ALT', $several ? 'VCF' : ())
		. $table->({ 2 => 1, 3 => 1, 4 => 1 }, $pairs, 'sample', 'sample', 'sites', 'IBS0 rate', 'kinship', 'relationship');
}

sub relatedness {
	my $sub = 'relatedness';
	# relatedness($file, ...) takes the VCFs first, unnamed, up to the first word:
	# a VCF's name has a dot, as its suffix must, and an option's name, or a
	# misspelling of one, has none
	my @vcf;
	push @vcf, shift @_ while @_ && defined $_[0] && !ref $_[0] && $_[0] !~ /\A\w+\z/;
	my $args = _args(\@_, $sub, [], [qw(vcf min_depth max_ibs0 report)]);
	if (defined $args->{vcf}) {
		croak "$sub: the VCFs are given either first, unnamed, or as \"vcf\", not both" if @vcf;
		@vcf = ref $args->{vcf} eq 'ARRAY' ? @{ $args->{vcf} } : $args->{vcf};
	}
	croak "$sub: needs a VCF of several samples, or several VCFs of one, as relatedness(\$file) or relatedness(\$file, \$file2)" unless @vcf;
	my %r;
	$r{min_depth} = $args->{min_depth} // 10; # 0 = count every genotype, whatever its DP; 10 is parent's, for parent's reasons
	croak "$sub: \"min_depth\" must be a whole number, 0 or more, not " . _what($r{min_depth})
		unless !ref $r{min_depth} && $r{min_depth} =~ /\A[0-9]+\z/;
	$r{max_ibs0} = $args->{max_ibs0} // 0.012; # parent's, and peddy 0.4.3's
	croak "$sub: \"max_ibs0\" must be a number above 0, not " . _what($r{max_ibs0})
		unless !ref $r{max_ibs0} && $r{max_ibs0} =~ $Number && $r{max_ibs0} > 0;
	$r{report} = $args->{report} // 1; # 1 = print the report to STDOUT; 0 = only return
	croak "$sub: \"report\" must be 0 or 1, not " . _what($r{report}) unless !ref $r{report} && ($r{report} eq '0' || $r{report} eq '1');
	# every file checked before any is read, which can take minutes
	my %given;
	foreach my $file (@vcf) {
		_vcf_check($sub, 'vcf', $file);
		croak "$sub: \"$file\" is given twice" if $given{$file}++;
	}
	my ($names, $counts, $tally); # the samples in order; REL_T counts a pair, as Bioinf.xs describes them, packed as UVs; each sample's sex counts
	if (@vcf == 1) {
		my ($next, $close) = _vcf_open($sub, 'vcf', $vcf[0]);
		($names, $counts, $tally) = _relatedness_one($next, $vcf[0], $r{min_depth});
		$close->();
	} else {
		# each opened as Bioinf.xs comes to it, and closed at its end, so that only
		# one is open at a time
		my @next = map {
			my $file = $_;
			my ($next, $close);
			sub {
				($next, $close) = _vcf_open($sub, 'vcf', $file) unless $next;
				my $block = $next->();
				$close->() unless defined $block;
				$block;
			};
		} @vcf;
		($names, $counts, $tally) = _relatedness_own(\@next, \@vcf, $r{min_depth});
	}
	foreach my $i (0 .. $#{ $names }) {
		my $s = $tally->{ $names->[$i] };
		_sex_of($s, $Sex_cutoff);
		$s->{file} = @vcf == 1 ? $vcf[0] : $vcf[$i];
	}
	my %out = map { $_ => {} } @{ $names };
	# each pair, in the samples' order, given to $code as (first, second, its
	# hash); walked when needed rather than held as a list, which for 1,000
	# samples would be 499,500 arrays, each with a copy of both names
	my $each_pair = sub {
		my ($code) = @_;
		foreach my $i (0 .. $#{ $names } - 1) {
			my $x = $names->[$i];
			$code->($x, $_, $out{$x}{$_}) foreach @{ $names }[$i + 1 .. $#{ $names }];
		}
	};
	my $pair = 11 * length pack 'J', 0; # bytes of a pair's counts in $counts
	my $k = 0; # where the pair's counts start in $counts
	foreach my $i (0 .. $#{ $names } - 1) {
		foreach my $j ($i + 1 .. $#{ $names }) {
			my ($x, $y) = @{ $names }[$i, $j];
			# t[3 * x's ALT count + y's], as parent's, and then the sites in both VCFs,
			# in x's alone and in y's alone
			my @t = (0, unpack "x$k J11", $counts);
			my @in = splice @t, 9;
			$k += $pair;
			# Flat, but for het and only: on 1,000 samples (2026-10-19) the pairs took
			# 560 MB with no field a hash, and 1,110 MB with four, which were these two
			# and each sample's hom-ALT count and share of ibs0.
			my %s;
			$s{sites} = 0;
			$s{sites} += $_ foreach @t;
			@in = ($s{sites}, 0, 0) if @vcf == 1;
			$s{in_both} = $in[0];
			$s{only} = { $x => $in[1], $y => $in[2] };
			$s{ibs0} = $t[2] + $t[6];
			$s{ibs1} = $t[1] + $t[3] + $t[5] + $t[7];
			$s{ibs2} = $t[4] + $t[8];
			$s{het_het} = $t[4];
			$s{het} = { $x => $t[3] + $t[4] + $t[5], $y => $t[1] + $t[4] + $t[7] };
			$s{kinship} = _king($s{het}{$x}, $s{het}{$y}, $s{het_het}, $s{ibs0});
			$s{ibs0_rate} = $s{sites} ? $s{ibs0} / $s{sites} : undef;
			# parent's reasons for dying, in parent's order
			$s{undecided} = $s{sites} && !($t[1] + $t[2] + $t[3] + $t[6])
				? 'neither has a hom-REF call where the other carries the ALT, so no opposite homozygote can be seen, as in a VCF merged from VCFs that list only their own variants; give each sample\'s own VCF instead'
				: $s{sites} < $Parent_sites ? "$s{sites} sites compared, and $Parent_sites are needed"
				: !defined $s{kinship} ? ($s{het}{$x} ? $y : $x) . ' has no het call at the sites compared'
				: undef;
			$s{relationship} = defined $s{undecided} ? undef : _relationship($s{kinship}, $s{ibs0_rate}, $r{max_ibs0});
			$s{parent} = $s{category} = undef; # below, once every pair is known
			$out{$x}{$y} = $out{$y}{$x} = \%s;
		}
	}
	undef $counts;
	# Which of a parent and child is the parent, where a third sample says, since
	# two people's genotypes cannot: x is y's parent when x has another child who
	# is y's full sibling, or when y has another parent who is unrelated to x, the
	# two being y's parents. Two of y's parent-child partners who are
	# second-degree relatives are y's parent and y's child, and say nothing about
	# which. "Third-degree" is not taken as unrelated: in the CEPH 1463 pedigree
	# (peddy's ceph1463.peddy.vcf.gz), the grandfather NA12889 and grandson
	# NA12893 had a kinship of 0.0839, below KING's 0.0884, and would have made
	# NA12893's father NA12877 his son too.
	my %pc; # sample => its parent-child partners
	$each_pair->(sub {
		my ($x, $y, $s) = @_;
		return unless ($s->{relationship} // '') eq 'parent-child';
		push @{ $pc{$x} }, $y;
		push @{ $pc{$y} }, $x;
	});
	my $is_parent = sub { # 1 or 0, and never a list, which would spill into the caller's next variable
		my ($x, $y) = @_;
		(grep { $_ ne $y && ($out{$y}{$_}{relationship} // '') eq 'full siblings' } @{ $pc{$x} })
			|| (grep { $_ ne $x && ($out{$x}{$_}{relationship} // '') eq 'unrelated' } @{ $pc{$y} }) ? 1 : 0;
	};
	$each_pair->(sub {
		my ($x, $y, $s) = @_;
		my ($rel, $sx, $sy) = ($s->{relationship}, $tally->{$x}{sex}, $tally->{$y}{sex});
		return unless defined $rel;
		if ($rel eq 'parent-child') {
			my ($x_is, $y_is) = ($is_parent->($x, $y), $is_parent->($y, $x));
			$s->{parent} = $x_is && !$y_is ? $x : $y_is && !$x_is ? $y : undef; # undef when neither is said, or both
			if (defined $s->{parent}) {
				my ($ps, $cs) = $s->{parent} eq $x ? ($sx, $sy) : ($sy, $sx);
				$s->{category} = ($ps ? $Kin_word{$ps}[0] : 'parent') . '-' . ($cs ? $Kin_word{$cs}[1] : 'child');
			} else { # two men are father and son, whichever is which
				$s->{category} = $sx && $sy && $sx eq $sy ? "$Kin_word{$sx}[0]-$Kin_word{$sx}[1]" : 'parent-child';
			}
		} elsif ($rel eq 'full siblings') {
			$s->{category} = !$sx || !$sy ? 'full siblings' : $sx ne $sy ? 'brother-sister' : $sx eq 'male' ? 'brothers' : 'sisters';
		} else {
			$s->{category} = $rel;
		}
	});
	if ($r{report}) {
		print STDOUT _relatedness_report(\@vcf, \%r, $names, $tally, $each_pair) or croak "$sub: can't print the report: $!";
	}
	my ($undecided, $first) = (0); # the pairs undecided, and the first of them as [first, second, why]
	my @regions; # as parent warns, and for its reasons: the pairs whose VCFs differ in their regions
	$each_pair->(sub {
		my ($x, $y, $s) = @_;
		if (defined $s->{undecided}) {
			$first //= [$x, $y, $s->{undecided}];
			$undecided++;
		}
		my ($more, $fewer) = $s->{only}{$x} > $s->{only}{$y} ? ($x, $y) : ($y, $x);
		push @regions, "$tally->{$more}{file} against $tally->{$fewer}{file} ($s->{only}{$more} to $s->{only}{$fewer})"
			if @vcf > 1 && $s->{only}{$more} > 2 * $s->{only}{$fewer};
	});
	my $n_pairs = @{ $names } * (@{ $names } - 1) / 2;
	carp "$sub: $undecided of $n_pairs pairs are undecided, and their relationship undef, such as $first->[0] and $first->[1]: $first->[2]" if $undecided;
	carp "$sub: the VCFs of " . scalar(@regions) . ' pair' . (@regions == 1 ? ' were' : 's were') . ' probably not sequenced over the same regions (an exome and a genome, or two capture kits),'
		. ' since one has more than twice as many sites that the other lacks as the other has that it lacks: '
		. join(', ', @regions[0 .. ($#regions < 2 ? $#regions : 2)]) . (@regions > 3 ? ', and ' . (@regions - 3) . ' more' : '')
		. '; a site that one lacks is taken as hom-REF in it, so their IBS0 and kinship are not to be trusted'
		if @regions;
	wantarray ? (\%out, $tally) : \%out;
}

# ---- genotype concordance

# concordance's report, given its files, the resolved options, its result,
# and the test's samples in its file's order. Columns are
# as wide as their widest entry, numbers right-aligned.
sub _concordance_report {
	my ($test, $truth, $r, $out, $test_names) = @_;
	my @pairs = map { $out->{$_} } grep { exists $out->{$_} } @{ $test_names }; # not a slice, which grep would autovivify into $out
	my $renamed = grep { $_->{test_sample} ne $_->{truth_sample} } @pairs; # the truth's sample is shown only when one is named otherwise
	my @head = ('sample', $renamed ? 'truth sample' : (), 'sites', 'in truth', 'concordant', 'discordant', 'gained', 'lost',
		'concordance', 'precision', 'recall', 'specificity');
	my $ratio = sub { defined $_[0] ? sprintf('%.4f', $_[0]) : '-' };
	my @rows = map {
		my ($c, $m) = @{ $_ }{qw(counts metrics)};
		my ($agree, $differ) = (0, 0);
		foreach my $kind (keys %{ $c }) {
			$agree += $c->{$kind} if $kind =~ /\Aconcordant_/;
			$differ += $c->{$kind} if $kind =~ /\Adiscordant_/;
		}
		[$_->{test_sample}, $renamed ? $_->{truth_sample} : (), $_->{sites}, $_->{in_truth}, $agree, $differ, @{ $c }{qw(gained_site lost_site)},
			map { $ratio->($m->{$_}) } qw(genotype_concordance precision_genotype recall_genotype specificity_genotype)]
	} @pairs;
	my @w = map { length } @head;
	foreach my $row (@rows) {
		foreach my $i (0 .. $#head) {
			$w[$i] = length $row->[$i] if length $row->[$i] > $w[$i];
		}
	}
	my $text_cols = $renamed ? 2 : 1;
	my $line = sub { join('  ', map { $_ < $text_cols ? sprintf('%-*s', $w[$_], $_[$_]) : sprintf('%*s', $w[$_], $_[$_]) } 0 .. $#head) . "\n" };
	my @alone = grep { !exists $out->{$_} } @{ $test_names };
	"concordance: test $test against truth $truth\n"
		. 'each test site compared with the truth\'s of the same CHROM, POS, REF and ALT, a site the truth lacks being a no-call there'
		. ($r->{min_depth} ? ", and a genotype of DP below $r->{min_depth} a no-call" : '')
		. "; precision, recall and specificity are of genotypes\n"
		. $line->(@head) . join('', map { $line->(@{ $_ }) } @rows)
		. (@alone ? scalar(@alone) . ' of the test\'s ' . scalar(@{ $test_names }) . ' samples are not in the truth: '
			. join(', ', @alone[0 .. ($#alone < 4 ? $#alone : 4)]) . (@alone > 5 ? ', and ' . (@alone - 5) . ' more' : '') . "\n" : '');
}

sub concordance {
	my $sub = 'concordance';
	my $args = _args(\@_, $sub, [qw(test truth)], [qw(discordances min_depth report)]);
	my %r;
	# x.low_ram_local_concordance.pl, which this is taken from, had no depth filter
	$r{min_depth} = $args->{min_depth} // 0; # 0 = every genotype compared, whatever its DP; otherwise a genotype of DP below it is a no-call
	croak "$sub: \"min_depth\" must be a whole number, 0 or more, not " . _what($r{min_depth})
		unless !ref $r{min_depth} && $r{min_depth} =~ /\A[0-9]+\z/;
	$r{discordances} = $args->{discordances} // 1; # 1 = list each discordant, gained and lost site; 0 = only count them, and "discordances" is undef
	croak "$sub: \"discordances\" must be 0 or 1, not " . _what($r{discordances})
		unless !ref $r{discordances} && ($r{discordances} eq '0' || $r{discordances} eq '1');
	$r{report} = $args->{report} // 1; # 1 = print a table of the samples to STDOUT; 0 = only return them
	croak "$sub: \"report\" must be 0 or 1, not " . _what($r{report}) unless !ref $r{report} && ($r{report} eq '0' || $r{report} eq '1');
	my ($test, $truth) = @{ $args }{qw(test truth)};
	# both checked before either is read, which can take minutes
	_vcf_check($sub, 'test', $test);
	_vcf_check($sub, 'truth', $truth);
	_same_build($sub, $test, $truth); # reads only the headers
	# the truth's samples first, so that the test's reading keeps the calls of
	# only those samples it will compare
	my ($next, $close) = _vcf_open($sub, 'truth', $truth);
	my $truth_names = _vcf_samples($next, $truth);
	$close->();
	my ($test_next, $test_close) = _vcf_open($sub, 'test', $test);
	my ($truth_next, $truth_close) = _vcf_open($sub, 'truth', $truth);
	my ($out, $test_names) = _concordance($test_next, $test, $truth_next, $truth, $truth_names, $r{min_depth}, $r{discordances});
	$test_close->();
	$truth_close->();
	if ($r{report}) {
		print STDOUT _concordance_report($test, $truth, \%r, $out, $test_names) or croak "$sub: can't print the report: $!";
	}
	my ($any) = values %{ $out }; # every pair has the same sites
	carp "$sub: \"test\" $test has no records, so nothing was compared" unless $any->{sites};
	carp "$sub: none of the $any->{sites} sites of \"test\" $test is in \"truth\" $truth, so every call of the test is gained or a no-call;"
		. ' check that the two are of the same genome build'
		if $any->{sites} && !$any->{in_truth};
	$out;
}

# ---- trait prediction

# GRCh37's and GRCh38's lengths for four of the chromosomes the HIrisPlex-S
# panel sits on, which _hirisplex_build tells one build from the other by.
# From Ensembl's /info/assembly/homo_sapiens, rest.ensembl.org for GRCh38 and
# grch37.rest.ensembl.org for GRCh37, read 2026-10-20. hg19 and b37 differ
# only in what they call their chromosomes, and hg38 and b38 likewise, so a
# length separates the two pairs but not the two within a pair, which the
# panel's positions do not depend on.
my %Build_length = (
	b37 => { 1 => 249250621, 15 => 102531392, 16 => 90354753, 20 => 63025520 },
	b38 => { 1 => 248956422, 15 => 101991189, 16 => 90338345, 20 => 64444167 },
);

# The HIrisPlex-S marker panel, share/hirisplex_s_panel.tsv, read once into an
# array ref of hash refs in the panel's own order. share/make_hirisplex_panel.pl
# writes that file and its header comment says where every column comes from.
my $Panel; # undef until the file has been read

sub _hirisplex_panel {
	my ($sub) = @_;
	return $Panel if $Panel;
	my $name = 'hirisplex_s_panel.tsv';
	# Installed, the panel is in the distribution's share directory, which is
	# what File::ShareDir finds. Two fallbacks follow it: the same layout looked
	# for in @INC, which is where blib puts the file under "make test" and keeps
	# this working on a perl that has not got File::ShareDir, and the share
	# directory beside lib/, which is what lets a working copy be used before
	# the module is installed.
	my $file = eval { require File::ShareDir; File::ShareDir::dist_file('Cn-Bioinf', $name) };
	unless (defined $file && -f $file) {
		my @try = map { File::Spec->catfile($_, 'auto', 'share', 'dist', 'Cn-Bioinf', $name) } grep { !ref } @INC;
		my (undef, $dir) = File::Spec->splitpath($INC{'Cn/Bioinf.pm'} // '');
		push @try, File::Spec->catfile($dir, File::Spec->updir, File::Spec->updir, 'share', $name) if length $dir;
		($file) = grep { -f $_ && -r _ } @try;
	}
	croak "$sub: can't find the marker panel $name, which this distribution installs in its share directory"
		unless defined $file && -f $file && -r _;
	open my $fh, '<', $file or croak "$sub: can't read $file: $!";
	my $head = <$fh>;
	croak "$sub: marker panel $file is empty" unless defined $head;
	$head =~ s/\r?\n\z//;
	$head =~ s/\A#//;
	my @cols = split /\t/, $head, -1;
	my @panel;
	while (defined(my $line = <$fh>)) {
		$line =~ s/\r?\n\z//;
		next if $line eq '';
		my %m;
		@m{@cols} = split /\t/, $line, -1;
		push @panel, \%m;
	}
	close $fh or croak "$sub: can't read $file: $!";
	croak "$sub: marker panel $file has no markers" unless @panel;
	foreach my $m (@panel) {
		foreach my $c (qw(rsid gene input_allele count_allele kind)) {
			croak "$sub: marker " . ($m->{rsid} // '?') . " of $file has no $c" unless defined $m->{$c} && length $m->{$c};
		}
		# Every marker of this panel is autosomal, and the reader counts a
		# diploid call only, so a sex chromosome here would be counted wrongly
		# rather than refused.
		foreach my $b (qw(b37 b38)) {
			my ($c, $p) = @{ $m }{"${b}_chrom", "${b}_pos"};
			next unless defined $c && length $c;
			croak "$sub: marker $m->{rsid} of $file is on $b chromosome \"$c\", and the panel is autosomal"
				unless $c =~ /\A(?:[1-9]|1[0-9]|2[0-2])\z/;
			croak "$sub: marker $m->{rsid} of $file has $b position " . _what($p) . ", which isn't a position"
				unless defined $p && $p =~ /\A[1-9][0-9]*\z/;
		}
	}
	$Panel = \@panel;
}

# Which build's coordinates the panel is matched on, read from the ##contig
# lines of the VCF $file, given to $sub as "vcf": "b37" or "b38". A header
# that gives no length for any of the chromosomes above, or lengths that are
# neither build's, cannot be told apart and dies asking for "build".
sub _hirisplex_build {
	my ($sub, $file) = @_;
	my $meta = _vcf_meta($sub, 'vcf', $file);
	my %length; # contig => the first length its header gives
	foreach my $c (@{ $meta->{contig} }) {
		$length{ $c->[0] } //= $c->[1] if defined $c->[1];
	}
	my @named = grep { defined $length{$_} } sort { $a <=> $b } keys %{ $Build_length{b37} };
	croak "$sub: \"vcf\" $file has no ##contig line giving a length for chromosome "
		. join(', ', sort { $a <=> $b } keys %{ $Build_length{b37} }) . ', so its genome build can\'t be told from its header;'
		. ' give "build" as hg19, b37, hg38 or b38' unless @named;
	my @fits = grep { my $b = $_; !grep { $length{$_} != $Build_length{$b}{$_} } @named } sort keys %Build_length;
	croak "$sub: \"vcf\" $file gives chromosome lengths that are neither GRCh37's nor GRCh38's ("
		. join(', ', map { "$_ is $length{$_}" } @named) . '), so its genome build can\'t be told from its header;'
		. ' give "build" as hg19, b37, hg38 or b38' unless @fits == 1;
	$fits[0];
}

# Whether each prediction can be made at all from $d, one genotype for each
# marker of $panel and undef where the sample had no call. The webtool refuses
# some outright, whatever the genotypes say, because the loss in AUC is too
# great; from the HIrisPlex-S user manual (version 2.0, 2018,
# https://hirisplex.erasmusmc.nl/pdf/hirisplex.erasmusmc.nl.pdf, read
# 2026-10-20), page 2: "Inputs missing the HERC2 SNP (rs12913832 = NA) will
# not provide an eye colour prediction result, as the loss in AUC is too
# great. Inputs missing all 11 MC1R SNPs will not provide a hair colour
# prediction result. Inputs missing a combination of HERC2-SLC45A2-IRF4 will
# not produce and eye or hair prediction result." The manual gives no such
# rule for skin, and says only that missing any marker "is not advised", so
# skin is refused here only when nothing at all was called.
sub _hirisplex_possible {
	my ($panel, $d) = @_;
	my %called; # rs number => 1 for each marker the sample has a genotype at
	my ($mc1r, $any) = (0, 0);
	foreach my $k (0 .. $#{ $panel }) {
		next unless defined $d->[$k];
		$called{ $panel->[$k]{rsid} } = 1;
		$mc1r++ if $panel->[$k]{gene} eq 'MC1R';
		$any++;
	}
	# HERC2 rs12913832, SLC45A2 rs16891982 and IRF4 rs12203592: eye needs the
	# first, so it needs one of the three already, and hair needs one of them
	my $three = grep { $called{$_} } qw(rs12913832 rs16891982 rs12203592);
	{
		eye  => $called{rs12913832} ? 1 : 0,
		hair => $mc1r && $three ? 1 : 0,
		skin => $any ? 1 : 0,
	};
}

# The genotypes of $out written to $file as the CSV the HIrisPlex-S webtool
# takes: a sampleid column and one for each marker, headed
# "<rs number>_<input allele>" as the webtool's own template heads it, holding
# the count of that allele or NA. The samples are written in $order, the VCF's.
sub _hirisplex_csv {
	my ($sub, $file, $panel, $out, $order) = @_;
	my @head = map { "$_->{rsid}_$_->{input_allele}" } @{ $panel };
	# The webtool reads the file as plain CSV, so a sample whose name holds a
	# comma, a quote or a newline would move or split its row rather than
	# being quoted into it.
	my @bad = grep { /[",\r\n]/ } @{ $order };
	croak "$sub: can't write $file: sample" . (@bad > 1 ? 's ' : ' ') . join(', ', map { "\"$_\"" } @bad)
		. (@bad > 1 ? ' have' : ' has') . ' a comma, a quote or a newline in the name, which a CSV column can\'t hold' if @bad;
	open my $fh, '>', $file or croak "$sub: can't write $file: $!";
	print {$fh} join(',', 'sampleid', @head), "\n" or croak "$sub: can't write $file: $!";
	foreach my $s (@{ $order }) {
		my $g = $out->{$s}{genotype};
		print {$fh} join(',', $s, map { defined $g->{$_} ? $g->{$_} : 'NA' } @head), "\n" or croak "$sub: can't write $file: $!";
	}
	close $fh or croak "$sub: can't write $file: $!";
}

# hirisplex_genotypes's report on $vcf: what was read and how, then a table of
# the samples in $order, the VCF's, each with how many markers were called and
# which predictions its genotypes allow.
sub _hirisplex_report {
	my ($vcf, $r, $panel, $out, $order, $absent) = @_;
	my @head = ('sample', 'called', 'missing', 'eye', 'hair', 'skin');
	my @rows = map {
		my $s = $out->{$_};
		[$_, $s->{called}, scalar @{ $s->{missing} }, map { $s->{$_} ? 'yes' : 'no' } qw(eye hair skin)]
	} @{ $order };
	my @w = map { length $head[$_] } 0 .. $#head;
	foreach my $row (@rows) {
		foreach my $i (0 .. $#head) {
			$w[$i] = length $row->[$i] if length $row->[$i] > $w[$i];
		}
	}
	my $line = sub { join('  ', map { $_ == 0 ? sprintf('%-*s', $w[$_], $_[$_]) : sprintf('%*s', $w[$_], $_[$_]) } 0 .. $#head) . "\n" };
	"hirisplex_genotypes: $vcf\n"
		. scalar(@{ $panel }) . ' HIrisPlex-S markers, matched by '
		. ($r->{match} eq 'id' ? 'rs number' : "position on $r->{build}")
		. "; genotypes of DP below $r->{min_depth} not counted\n"
		. (@{ $absent } ? scalar(@{ $absent }) . ' marker' . (@{ $absent } > 1 ? 's' : '') . ' not in the file at all: '
			. join(', ', @{ $absent }) . "\n" : '')
		. $line->(@head) . join('', map { $line->(@{ $_ }) } @rows)
		. "\"eye\", \"hair\" and \"skin\" are whether the HIrisPlex-S manual's rules allow that prediction from these markers\n";
}

sub hirisplex_genotypes {
	my $sub = 'hirisplex_genotypes';
	# hirisplex_genotypes($file) and hirisplex_genotypes($file, %options) take
	# the file first, as sex_check does; an odd count whose first is an
	# option's name is left for _args to refuse
	unshift @_, 'vcf' if @_ % 2 && defined $_[0] && !ref $_[0] && $_[0] !~ /\A(?:vcf|build|file|match|min_depth|report)\z/;
	croak "$sub: \"$_[0]\" is the file, so the options after it must be name => value pairs, and are an odd number"
		if @_ && $_[0] ne 'vcf' && defined $_[0] && !ref $_[0] && -f $_[0];
	my $args = _args(\@_, $sub, ['vcf'], [qw(build file match min_depth report)]);
	my %r;
	$r{match} = $args->{match} // 'position'; # 'position' = CHROM and POS against the panel's coordinates; 'id' = the ID column against its rs numbers
	croak "$sub: \"match\" must be \"position\" or \"id\", not " . _what($r{match})
		unless !ref $r{match} && ($r{match} eq 'position' || $r{match} eq 'id');
	$r{min_depth} = $args->{min_depth} // 7; # 0 = count every genotype, whatever its DP; 7 is sex_check's, for sex_check's reasons
	croak "$sub: \"min_depth\" must be a whole number, 0 or more, not " . _what($r{min_depth})
		unless !ref $r{min_depth} && $r{min_depth} =~ /\A[0-9]+\z/;
	$r{report} = $args->{report} // 1; # 1 = print a table of the samples to STDOUT; 0 = only return them
	croak "$sub: \"report\" must be 0 or 1, not " . _what($r{report}) unless !ref $r{report} && ($r{report} eq '0' || $r{report} eq '1');
	$r{file} = $args->{file}; # undef = no CSV is written
	croak "$sub: \"file\" must be a file name, not " . _what($r{file}) if defined $r{file} && ref $r{file};
	croak "$sub: \"build\" must be hg19, b37, hg38 or b38, not " . _what($args->{build})
		if defined $args->{build} && (ref $args->{build} || $args->{build} !~ /\A(?:hg19|b37|hg38|b38)\z/);
	my $vcf = $args->{vcf};
	_vcf_check($sub, 'vcf', $vcf);
	my $panel = _hirisplex_panel($sub);
	# which of the panel's two sets of coordinates the markers are looked for
	# at; undef when matching by rs number, which needs neither
	$r{build} = $r{match} eq 'id' ? undef
		: defined $args->{build} ? ($args->{build} =~ /\A(?:hg19|b37)\z/ ? 'b37' : 'b38')
		: _hirisplex_build($sub, $vcf);

	# What the reader searches for: one entry for each marker, and for matching
	# by rs number a second under the number that marker's has been merged
	# into, since a VCF annotated from either dbSNP release should match.
	my @entries;
	foreach my $k (0 .. $#{ $panel }) {
		my $m = $panel->[$k];
		my $insertion = $m->{kind} eq 'insertion' ? 1 : 0; # 1 = count the ALT longer than REF, 0 = count count_allele
		if ($r{match} eq 'id') {
			my %seen;
			foreach my $id (grep { defined $_ && length $_ && !$seen{$_}++ } $m->{rsid}, $m->{current_rsid}) {
				push @entries, [0, 0, $m->{count_allele}, $insertion, $id, $k];
			}
		} else {
			my ($c, $p) = @{ $m }{"$r{build}_chrom", "$r{build}_pos"};
			# a marker the panel could not place on this build is searched for
			# nowhere, and comes back as one the file did not hold
			next unless defined $c && length $c;
			push @entries, [$c, $p, $m->{count_allele}, $insertion, $m->{rsid}, $k];
		}
	}
	croak "$sub: no marker of the panel has a " . ($r{build} // 'known') . ' position' unless @entries;

	my ($next, $close) = _vcf_open($sub, 'vcf', $vcf);
	my ($dose, $names, $seen) = _hirisplex_dosage($next, $vcf, $r{min_depth}, \@entries, scalar @{ $panel }, $r{match} eq 'id' ? 1 : 0);
	$close->();

	my @column = map { "$_->{rsid}_$_->{input_allele}" } @{ $panel };
	my @absent = map { $panel->[$_]{rsid} } grep { !$seen->[$_] } 0 .. $#{ $panel };
	my %out;
	foreach my $s (keys %{ $dose }) {
		my $d = $dose->{$s};
		my (%genotype, @missing);
		foreach my $k (0 .. $#{ $panel }) {
			$genotype{ $column[$k] } = $d->[$k];
			push @missing, $panel->[$k]{rsid} unless defined $d->[$k];
		}
		my $can = _hirisplex_possible($panel, $d);
		$out{$s} = {
			genotype => \%genotype,
			missing  => \@missing,
			markers  => scalar @{ $panel },
			called   => scalar(@{ $panel }) - scalar(@missing),
			%{ $can },
		};
	}
	_hirisplex_csv($sub, $r{file}, $panel, \%out, $names) if defined $r{file};
	if ($r{report}) {
		print STDOUT _hirisplex_report($vcf, \%r, $panel, \%out, $names, \@absent) or croak "$sub: can't print the report: $!";
	}
	carp "$sub: \"vcf\" $vcf holds none of the " . scalar(@{ $panel }) . ' HIrisPlex-S markers'
		. ($r{match} eq 'id' ? '; it was matched by rs number, and a VCF whose ID column is "." matches nothing that way'
			: ", matched by position on $r{build}; check that it is the build the file is on")
		if @absent == @{ $panel };
	\%out;
}

# The marker that carries a category's intercept, in place of an rs number, in
# a "model" table. A marker can't be called this, so nothing is ambiguous.
my $Intercept = '(intercept)';

# The rows of "model" read into trait => category => { marker => beta }, with
# the intercept under $Intercept, and the traits and their categories in the
# order the table first names them, which is the order the report prints.
sub _trait_model {
	my ($sub, $model) = @_;
	my ($next, $cols) = _weight_rows($sub, $model, 'model');
	my %have = map { $_ => 1 } @{ $cols };
	my @absent = grep { !$have{$_} } qw(trait category marker beta);
	croak "$sub: \"model\" has no " . join(', ', map { "\"$_\"" } @absent) . ' column' . (@absent > 1 ? 's' : '')
		. '; a model table has trait, category, marker and beta' if @absent;
	my (%beta, @traits, %category);
	while (my ($row, $label) = $next->()) {
		my ($t, $c, $m, $b) = @{ $row }{qw(trait category marker beta)};
		foreach my $f (['trait', $t], ['category', $c], ['marker', $m]) {
			croak "$sub: $label of \"model\" has no $f->[0]" unless defined $f->[1] && length $f->[1];
		}
		croak "$sub: $label of \"model\" has beta " . _what($b) . ", which isn't a number"
			unless defined $b && !ref $b && $b =~ $Number;
		push @traits, $t unless $beta{$t};
		push @{ $category{$t} }, $c unless $beta{$t} && $beta{$t}{$c};
		croak "$sub: $label of \"model\" gives $t category \"$c\" a second beta for $m"
			if exists $beta{$t}{$c}{$m};
		$beta{$t}{$c}{$m} = $b;
	}
	croak "$sub: \"model\" has no rows" unless @traits;
	(\%beta, \@traits, \%category);
}

# exp() of each of @$eta, divided by their sum: the probabilities a multinomial
# logistic regression gives its categories. The largest is taken out first, so
# that a linear predictor big enough to overflow exp() (about 710 on an IEEE
# double, and larger on a long-double or quadmath perl) gives 1 and the rest 0
# rather than inf/inf, which is NaN.
sub _softmax {
	my ($eta) = @_;
	my $max = $eta->[0];
	foreach my $e (@{ $eta }) { $max = $e if $e > $max }
	my @exp = map { exp($_ - $max) } @{ $eta };
	my $sum = 0;
	$sum += $_ foreach @exp;
	[map { $_ / $sum } @exp];
}

# trait_predict's report on $out: a table of each sample's probabilities, one
# block for each trait, the categories in the model's order.
sub _trait_report {
	my ($r, $out, $order, $traits, $category) = @_;
	my $text = "trait_predict: " . scalar(@{ $order }) . ' sample' . (@{ $order } > 1 ? 's' : '')
		. ', ' . scalar(@{ $traits }) . ' trait' . (@{ $traits } > 1 ? 's' : '')
		. "; a marker the sample has no genotype for counts as 0 copies"
		. ($r->{force} ? ', and a prediction the HIrisPlex-S rules refuse is made anyway ("force")' : '') . "\n";
	foreach my $t (@{ $traits }) {
		my @head = ('sample', 'call', @{ $category->{$t} }, 'called');
		my @rows = map {
			my $p = $out->{$_}{$t};
			[$_, $p->{call} // '-', (map { defined $p->{p}{$_} ? sprintf('%.3f', $p->{p}{$_}) : '-' } @{ $category->{$t} }),
				"$p->{called}/$p->{markers}"]
		} @{ $order };
		my @w = map { length $head[$_] } 0 .. $#head;
		foreach my $row (@rows) {
			foreach my $i (0 .. $#head) {
				$w[$i] = length $row->[$i] if length $row->[$i] > $w[$i];
			}
		}
		my $line = sub { join('  ', map { $_ == 0 ? sprintf('%-*s', $w[$_], $_[$_]) : sprintf('%*s', $w[$_], $_[$_]) } 0 .. $#head) . "\n" };
		$text .= "\n$t\n" . $line->(@head) . join('', map { $line->(@{ $_ }) } @rows);
		my @refused = grep { defined $out->{$_}{$t}{refused} } @{ $order };
		$text .= join '', map { "  $_: no $t prediction: $out->{$_}{$t}{refused}\n" } @refused;
	}
	$text;
}

sub trait_predict {
	my $sub = 'trait_predict';
	my $args = _args(\@_, $sub, [qw(genotypes model)], [qw(force report)]);
	my %r;
	$r{force} = $args->{force} // 0; # 0 = a trait the genotypes say can't be predicted isn't; 1 = predict it anyway
	croak "$sub: \"force\" must be 0 or 1, not " . _what($r{force}) unless !ref $r{force} && ($r{force} eq '0' || $r{force} eq '1');
	$r{report} = $args->{report} // 1; # 1 = print a table of the samples to STDOUT; 0 = only return them
	croak "$sub: \"report\" must be 0 or 1, not " . _what($r{report}) unless !ref $r{report} && ($r{report} eq '0' || $r{report} eq '1');
	my $g = $args->{genotypes};
	croak "$sub: \"genotypes\" must be a hash ref of sample => genotypes, as hirisplex_genotypes returns, not " . _what($g)
		unless ref $g eq 'HASH';
	croak "$sub: \"genotypes\" has no samples" unless keys %{ $g };
	my ($beta, $traits, $category) = _trait_model($sub, $args->{model});

	# Which markers each trait's categories use between them: the same for every
	# sample, so it is worked out once rather than once per sample.
	my %uses;
	foreach my $t (@{ $traits }) {
		foreach my $c (@{ $category->{$t} }) {
			$uses{$t}{$_} = 1 foreach grep { $_ ne $Intercept } keys %{ $beta->{$t}{$c} };
		}
	}

	my %out;
	foreach my $s (sort keys %{ $g }) {
		my $rec = $g->{$s};
		croak "$sub: \"genotypes\" sample \"$s\" must be a hash ref, not " . _what($rec) unless ref $rec eq 'HASH';
		# Either what hirisplex_genotypes returns, whose genotypes are under
		# "genotype" and which also says which predictions its markers allow,
		# or a plain hash of marker => dosage.
		my $dose = ref $rec->{genotype} eq 'HASH' ? $rec->{genotype} : $rec;
		foreach my $t (@{ $traits }) {
			my (@eta, @missing);
			my $called = 0;
			my $marker = $uses{$t};
			foreach my $m (sort keys %{ $marker }) {
				if (defined $dose->{$m}) { $called++ } else { push @missing, $m }
			}
			foreach my $c (@{ $category->{$t} }) {
				my $eta = $beta->{$t}{$c}{$Intercept} // 0;
				foreach my $m (keys %{ $beta->{$t}{$c} }) {
					next if $m eq $Intercept;
					# A marker with no genotype drops out of the sum, which is
					# the same as counting it 0 copies of the counted allele.
					$eta += $beta->{$t}{$c}{$m} * $dose->{$m} if defined $dose->{$m};
				}
				push @eta, $eta;
			}
			# hirisplex_genotypes sets eye, hair and skin to 0 where the
			# HIrisPlex-S manual refuses the prediction outright; a model whose
			# trait it says nothing about is predicted whatever it holds.
			my $refused = !$r{force} && defined $rec->{$t} && !ref $rec->{$t} && $rec->{$t} eq '0'
				? 'the markers it needs are missing, and the HIrisPlex-S rules refuse it' : undef;
			my $p = $refused ? {} : do {
				my $soft = _softmax(\@eta);
				my %p;
				@p{ @{ $category->{$t} } } = @{ $soft };
				\%p;
			};
			my ($call) = $refused ? () : sort { $p->{$b} <=> $p->{$a} || $a cmp $b } keys %{ $p };
			$out{$s}{$t} = {
				p       => $p,
				call    => $call,
				markers => scalar keys %{ $marker },
				called  => $called,
				missing => \@missing,
				refused => $refused,
			};
		}
	}
	if ($r{report}) {
		print STDOUT _trait_report(\%r, \%out, [sort keys %out], $traits, $category) or croak "$sub: can't print the report: $!";
	}
	\%out;
}

# ---- PDB structures
#
# pdb_search asks two RCSB APIs, as their documentation at
# https://search.rcsb.org/ and https://data.rcsb.org/ described them on
# 2026-10-21: the search API finds the polymer entities whose sequence
# matches, by MMseqs2, and the data API's GraphQL says what each one's
# structure is.

my $Rcsb_search = 'https://search.rcsb.org/rcsbsearch/v2/query';
my $Rcsb_data   = 'https://data.rcsb.org/graphql';
# The search API's own limit on "rows", which it refuses past with "numeric
# instance is greater than the required maximum (maximum: 10000" (2026-10-21).
my $Rcsb_rows = 10000;
# How many entities each GraphQL request asks about. 500 came back in one
# request in 2.45 s, 1.6 MB of JSON, on 2026-10-21, so 200 leaves headroom
# for entities with more annotations than those had.
my $Rcsb_batch = 200;
my $Rcsb_entities = <<'GRAPHQL';
query ($ids: [String!]!) { polymer_entities(entity_ids: $ids) {
	rcsb_id
	entity_poly { rcsb_entity_polymer_type pdbx_strand_id }
	rcsb_polymer_entity { pdbx_description }
	rcsb_entity_source_organism { scientific_name ncbi_taxonomy_id }
	rcsb_polymer_entity_container_identifiers { reference_sequence_identifiers { database_name database_accession } }
	rcsb_polymer_entity_annotation { type annotation_id name }
	entry {
		struct { title }
		rcsb_entry_info { experimental_method resolution_combined polymer_composition structure_determination_methodology }
		rcsb_accession_info { deposit_date initial_release_date }
		rcsb_ma_qa_metric_global { ma_qa_metric_global { type value } }
		assemblies { rcsb_struct_symmetry { kind oligomeric_state stoichiometry symbol } }
		nonpolymer_entities { nonpolymer_comp { chem_comp { id name } } }
	}
} }
GRAPHQL

# $body POSTed to $url as JSON, for $sub, and the JSON it answers with,
# decoded; undef when it answers 204, which is how the search API says that
# nothing matched. A failure dies with RCSB's own "message" where it sent one.
# The tests replace this sub, so that they run with no network.
sub _rcsb_post {
	my ($sub, $url, $body, $timeout) = @_;
	require HTTP::Tiny;
	require JSON::MaybeXS;
	# verify_SSL is given, since it is the default only from HTTP::Tiny 0.083
	my $http = HTTP::Tiny->new(agent => "Cn-Bioinf/$VERSION ", timeout => $timeout, verify_SSL => 1);
	my ($ssl, $why) = $http->can_ssl;
	unless ($ssl) {
		$why =~ s/\s+\z//;
		croak "$sub: can't reach RCSB, which needs HTTPS: $why";
	}
	my $res = $http->post($url, {
		headers => { 'Content-Type' => 'application/json', Accept => 'application/json' },
		content => JSON::MaybeXS::encode_json($body),
	});
	return undef if $res->{status} == 204;
	my $json = eval { JSON::MaybeXS::decode_json($res->{content}) };
	unless ($res->{success}) {
		# 599 is HTTP::Tiny's own status, for a request that got no answer at
		# all, and its content says why
		my $msg = $res->{status} == 599 ? $res->{content}
			: ref $json eq 'HASH' && defined $json->{message} ? $json->{message} : $res->{reason};
		$msg =~ s/\s+\z//;
		croak "$sub: $url answered $res->{status}: $msg";
	}
	croak "$sub: $url answered with something that isn't JSON: " . substr($res->{content}, 0, 200) unless ref $json;
	$json;
}

# The alignment of one search hit, from the "match_context" RCSB gave it,
# under the names get_best_alignment_hit gives a BLAST hsp. RCSB's
# sequence_identity is truncated to 3 places, so the identity is counted from
# the aligned sequences instead, which is how BLAST counts it.
sub _pdb_alignment {
	my ($m) = @_;
	my ($q, $h) = @{ $m }{qw(query_aligned_seq subject_aligned_seq)};
	my ($identity, $gaps) = (0, 0);
	foreach my $i (0 .. length($q) - 1) {
		my ($x, $y) = (substr($q, $i, 1), substr($h, $i, 1));
		if ($x eq '-' || $y eq '-') { $gaps++ } elsif ($x eq $y) { $identity++ }
	}
	return {
		evalue      => $m->{evalue},
		bit_score   => $m->{bitscore},
		align_len   => $m->{alignment_length},
		identity    => $identity,
		pident      => 100 * $identity / $m->{alignment_length},
		mismatches  => $m->{mismatches},
		gaps        => $gaps,
		gap_opens   => $m->{gaps_opened},
		query_from  => $m->{query_beg},
		query_to    => $m->{query_end},
		query_len   => $m->{query_length},
		query_cover => 100 * ($m->{query_end} - $m->{query_beg} + 1) / $m->{query_length},
		hit_from    => $m->{subject_beg},
		hit_to      => $m->{subject_end},
		hit_len     => $m->{subject_length},
		qseq        => $q,
		hseq        => $h,
	};
}

# What the data API said of one polymer entity, $e, which is undef when it
# said nothing, as the fields pdb_search returns: every key is there either
# way, undef or an empty array when RCSB has no value.
sub _pdb_structure {
	my ($e) = @_;
	$e //= {};
	my $entry = $e->{entry} // {};
	my $info = $entry->{rcsb_entry_info} // {};
	my $dates = $entry->{rcsb_accession_info} // {};
	# The global symmetry of the first assembly, which is the one the
	# depositors gave first, and the one RCSB's own page shows.
	my ($sym) = grep { ($_->{kind} // '') eq 'Global Symmetry' }
		map { @{ $_->{rcsb_struct_symmetry} // [] } } grep { ref $_ eq 'HASH' } @{ $entry->{assemblies} // [] }[0];
	$sym //= {};
	my ($plddt) = map { $_->{value} } grep { ($_->{type} // '') eq 'pLDDT' }
		map { @{ $_->{ma_qa_metric_global} // [] } } @{ $entry->{rcsb_ma_qa_metric_global} // [] };
	my @annotation = @{ $e->{rcsb_polymer_entity_annotation} // [] };
	my (%seen, @ligands);
	foreach my $c (map { $_->{nonpolymer_comp}{chem_comp} // () } @{ $entry->{nonpolymer_entities} // [] }) {
		push @ligands, { id => $c->{id}, name => $c->{name} } unless $seen{ $c->{id} }++;
	}
	my $chains = ($e->{entity_poly} // {})->{pdbx_strand_id};
	return {
		title               => ($entry->{struct} // {})->{title},
		description         => ($e->{rcsb_polymer_entity} // {})->{pdbx_description},
		method              => $info->{experimental_method},
		methodology         => $info->{structure_determination_methodology},
		resolution          => ($info->{resolution_combined} // [])->[0],
		plddt               => $plddt,
		polymer_type        => ($e->{entity_poly} // {})->{rcsb_entity_polymer_type},
		polymer_composition => $info->{polymer_composition},
		oligomeric_state    => $sym->{oligomeric_state},
		stoichiometry       => $sym->{stoichiometry} ? join('', @{ $sym->{stoichiometry} }) : undef,
		symmetry            => $sym->{symbol},
		chains              => [defined $chains ? split(/,/, $chains) : ()],
		organisms           => [map { +{ name => $_->{scientific_name}, taxonomy_id => $_->{ncbi_taxonomy_id} } } @{ $e->{rcsb_entity_source_organism} // [] }],
		uniprot             => [map { $_->{database_accession} } grep { ($_->{database_name} // '') eq 'UniProt' }
			@{ ($e->{rcsb_polymer_entity_container_identifiers} // {})->{reference_sequence_identifiers} // [] }],
		pfam                => [map { +{ id => $_->{annotation_id}, name => $_->{name} } } grep { ($_->{type} // '') eq 'Pfam' } @annotation],
		interpro            => [map { +{ id => $_->{annotation_id}, name => $_->{name} } } grep { ($_->{type} // '') eq 'InterPro' } @annotation],
		ligands             => \@ligands,
		deposited           => defined $dates->{deposit_date} ? substr($dates->{deposit_date}, 0, 10) : undef,
		released            => defined $dates->{initial_release_date} ? substr($dates->{initial_release_date}, 0, 10) : undef,
	};
}

# How a hit's structure was determined, as the report shows it: RCSB's
# experimental method, which a computed model has none of. Its values, and
# those of the methodology, are the enums of rcsb_entry_info in
# https://data.rcsb.org/rest/v1/schema/entry (checked 2026-10-21).
sub _pdb_method {
	my ($h) = @_;
	return $h->{method} if defined $h->{method};
	($h->{methodology} // '') eq 'computational' ? 'model' : '-';
}

# pdb_search's report on @$hits, $total of which RCSB found: one line a hit,
# best first.
sub _pdb_report {
	my ($r, $hits, $total) = @_;
	my $text = "pdb_search: " . ($total == 0 ? 'no hits' : scalar(@{ $hits }) . " of $total hit" . ($total > 1 ? 's' : ''))
		. " in RCSB's " . ($r->{computed} ? 'PDB and computed models' : 'PDB')
		. " for a $r->{len}-residue $r->{type} sequence, at an E-value below $r->{evalue}"
		. ($r->{identity} > 0 ? " and $r->{identity}% identity or more" : '') . "\n";
	return $text unless @{ $hits };
	my @head = ('id', 'method', 'resolution', 'identity', 'cover', 'evalue', 'bits', 'oligomer', 'description');
	my @rows = map {
		[$_->{id}, _pdb_method($_),
			defined $_->{resolution} ? sprintf('%.2f A', $_->{resolution}) : defined $_->{plddt} ? "pLDDT $_->{plddt}" : '-',
			sprintf('%.1f%%', $_->{pident}), sprintf('%.0f%%', $_->{query_cover}), sprintf('%.2g', $_->{evalue}), $_->{bit_score},
			$_->{oligomeric_state} // '-', $_->{description} // '-']
	} @{ $hits };
	my @w = map { length $head[$_] } 0 .. $#head;
	foreach my $row (@rows) {
		foreach my $i (0 .. $#head) {
			$w[$i] = length $row->[$i] if length $row->[$i] > $w[$i];
		}
	}
	# left-aligned text, right-aligned numbers; the last column is not padded
	my %left = map { $_ => 1 } 0, 1, 7, 8;
	my $line = sub { join('  ', map { $_ == $#head ? $_[$_] : $left{$_} ? sprintf('%-*s', $w[$_], $_[$_]) : sprintf('%*s', $w[$_], $_[$_]) } 0 .. $#head) . "\n" };
	$text . $line->(@head) . join('', map { $line->(@{ $_ }) } @rows);
}

sub pdb_search {
	my $sub = 'pdb_search';
	my $args = _args(\@_, $sub, [qw(seq)], [qw(computed evalue identity max_hits report timeout type)]);
	my %r;
	# The cutoffs are made numbers with 0 +, since JSON::XS writes a string as
	# a string, and RCSB refuses a string where its schema wants a number.
	$r{evalue} = $args->{evalue} // 0.1; # hits at this E-value or above are left out; BLAST's default of 10 would let in hits expected 10 times by chance
	croak "$sub: \"evalue\" must be a number above 0, not " . _what($r{evalue})
		unless !ref $r{evalue} && $r{evalue} =~ $Number && $r{evalue} > 0;
	$r{evalue} = 0 + $r{evalue};
	$r{identity} = $args->{identity} // 0; # percent; hits with less identity than this are left out
	croak "$sub: \"identity\" must be a percentage, from 0 to 100, not " . _what($r{identity})
		unless !ref $r{identity} && $r{identity} =~ $Number && $r{identity} >= 0 && $r{identity} <= 100;
	$r{identity} = 0 + $r{identity};
	$r{max_hits} = $args->{max_hits} // 100; # the most hits returned, best first
	croak "$sub: \"max_hits\" must be a whole number from 1 to $Rcsb_rows, RCSB's limit, not " . _what($r{max_hits})
		unless !ref $r{max_hits} && $r{max_hits} =~ /\A[0-9]+\z/ && $r{max_hits} >= 1 && $r{max_hits} <= $Rcsb_rows;
	$r{max_hits} = 0 + $r{max_hits};
	$r{computed} = $args->{computed} // 0; # 0 = experimental structures only; 1 = also computed models (AlphaFold DB, ModelArchive)
	croak "$sub: \"computed\" must be 0 or 1, not " . _what($r{computed}) unless !ref $r{computed} && ($r{computed} eq '0' || $r{computed} eq '1');
	$r{type} = $args->{type} // 'protein'; # protein, dna or rna: which of RCSB's sequence databases to search
	croak "$sub: \"type\" must be protein, dna or rna, not " . _what($r{type}) unless !ref $r{type} && $r{type} =~ /\A(?:protein|dna|rna)\z/;
	$r{report} = $args->{report} // 1; # 1 = print a table of the hits to STDOUT; 0 = only return them
	croak "$sub: \"report\" must be 0 or 1, not " . _what($r{report}) unless !ref $r{report} && ($r{report} eq '0' || $r{report} eq '1');
	$r{timeout} = $args->{timeout} // 60; # seconds for each request, HTTP::Tiny's own default
	croak "$sub: \"timeout\" must be a number of seconds above 0, not " . _what($r{timeout})
		unless !ref $r{timeout} && $r{timeout} =~ $Number && $r{timeout} > 0;

	my $seq = $args->{seq};
	croak "$sub: \"seq\" must be a sequence, as a string, not " . _what($seq) if ref $seq;
	# A sequence copied out of a FASTA file keeps its line breaks, and one
	# translated from DNA can end in a stop; RCSB refuses both, and lower case.
	$seq =~ s/\s+//g;
	$seq = uc $seq;
	$seq =~ s/\*\z//;
	croak "$sub: \"seq\" has no residues" unless length $seq;
	croak "$sub: \"seq\" has \"$1\" at residue " . ($-[0] + 1) . ', which isn\'t a residue letter' if $seq =~ /([^A-Z])/;
	$r{len} = length $seq;

	my $found = _rcsb_post($sub, $Rcsb_search, {
		query => { type => 'terminal', service => 'sequence', parameters => {
			evalue_cutoff => $r{evalue}, identity_cutoff => $r{identity} / 100, sequence_type => $r{type}, value => $seq,
		} },
		request_options => {
			results_verbosity    => 'verbose', # which gives each hit's alignment
			results_content_type => [$r{computed} ? ('experimental', 'computational') : 'experimental'],
			scoring_strategy     => 'sequence',
			paginate             => { start => 0, rows => $r{max_hits} },
		},
		return_type => 'polymer_entity',
	}, $r{timeout});
	my $total = $found ? $found->{total_count} : 0;
	my @hits;
	foreach my $res (@{ $found ? $found->{result_set} : [] }) {
		my $best; # undef until a match is seen
		foreach my $m (map { @{ $_->{match_context} // [] } } map { @{ $_->{nodes} // [] } } @{ $res->{services} // [] }) {
			# none of the 1382 hits for human haemoglobin alpha had more than
			# one match on 2026-10-21, but the format allows it
			$best = $m if !defined $best || $m->{bitscore} > $best->{bitscore};
		}
		croak "$sub: RCSB gave no alignment for $res->{identifier}" unless defined $best;
		# An identifier is entry_entity: 1BZ1_1, or AF_AFP69905F1_1 for a computed model
		my ($entry, $entity) = $res->{identifier} =~ /\A(.+)_([^_]+)\z/
			or croak "$sub: RCSB returned \"$res->{identifier}\", which isn't an entry_entity identifier";
		push @hits, { id => $res->{identifier}, entry => $entry, entity => $entity,
			url => "https://www.rcsb.org/structure/$entry", %{ _pdb_alignment($best) } };
	}

	my %entity;
	for (my $i = 0; $i < @hits; $i += $Rcsb_batch) {
		my @ids = map { $_->{id} } @hits[$i .. ($i + $Rcsb_batch > @hits ? $#hits : $i + $Rcsb_batch - 1)];
		my $data = _rcsb_post($sub, $Rcsb_data, { query => $Rcsb_entities, variables => { ids => \@ids } }, $r{timeout});
		my $list = ref $data eq 'HASH' && ref $data->{data} eq 'HASH' ? $data->{data}{polymer_entities} : undef;
		unless (ref $list eq 'ARRAY') {
			my @why = map { $_->{message} // () } @{ ref $data eq 'HASH' && ref $data->{errors} eq 'ARRAY' ? $data->{errors} : [] };
			croak "$sub: $Rcsb_data said nothing of $ids[0]" . (@ids > 1 ? ' and the rest' : '') . (@why ? ': ' . join('; ', @why) : '');
		}
		$entity{ $_->{rcsb_id} } = $_ foreach grep { ref $_ eq 'HASH' } @{ $list };
	}
	foreach my $h (@hits) {
		# the search index can be ahead of the data API for a new release
		carp "$sub: $Rcsb_data had nothing on $h->{id}, so its structure fields are empty" unless $entity{ $h->{id} };
		%{ $h } = (%{ $h }, %{ _pdb_structure($entity{ $h->{id} }) });
	}
	# Best first, as get_best_alignment_hit sorts on "evalue": RCSB rounds a
	# strong E-value to 0, so the bit score breaks ties, and then RCSB's order.
	my @order = sort { $hits[$a]{evalue} <=> $hits[$b]{evalue} || $hits[$b]{bit_score} <=> $hits[$a]{bit_score} || $a <=> $b } 0 .. $#hits;
	@hits = @hits[@order];
	if ($r{report}) {
		my $text = _pdb_report(\%r, \@hits, $total);
		# RCSB's text is characters, and can hold some above 255 (a Greek
		# letter in a title), which print would warn of as a "Wide character",
		# fatally in a caller under FATAL warnings, unless STDOUT encodes.
		utf8::encode($text) unless grep { /\A(?:utf8|encoding)/ } PerlIO::get_layers(*STDOUT);
		print STDOUT $text or croak "$sub: can't print the report: $!";
	}
	wantarray ? (\@hits, $total) : \@hits;
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
	my $args = _args(\@_, $sub, [qw(fasta filename)], [qw(active_site_aa clustal_args execute fh labels msa_file order query threads title
		tree_file xlabel ylabel)]);
	_refs($sub, $args, 'ARRAY', qw(clustal_args order));
	_refs($sub, $args, 'HASH', 'labels');
	my $drawing = _drawing($sub, $args);
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
	$r->{'python_file'} = _draw($sub, $drawing, '--f', $fa, '--o', $args->{filename}, '--c', _creator($sub), '--quiet',
		'--p', _json_text([map { _chars($_) } _align_provenance($args, $sub, $cmd, $names, $r->{'msa_file'}), _kept($args, qw(msa_file tree_file))]),
		(defined $args->{title} ? ('--t', $args->{title}) : ()),
		(defined $args->{xlabel} ? ('--x', $args->{xlabel}) : ()),
		(defined $args->{ylabel} ? ('--y', $args->{ylabel}) : ()), @segments);
	_wrote($args->{filename}, 43, $drawing->{execute} ? 'wrote' : 'will write'); # yellow
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
	my $args = _args(\@_, $sub, [], [qw(clustal_args execute fasta fh labels msa_file output_file threads title tree_file)]);
	_refs($sub, $args, 'ARRAY', 'clustal_args');
	_refs($sub, $args, 'HASH', 'labels');
	my $drawing = _drawing($sub, $args);
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
	$r->{'python_file'} = _draw($sub, $drawing, '--tree', $tree, '--o', $args->{'output_file'}, '--c', _creator($sub), '--quiet',
		'--p', _json_text([map { _chars($_) } @prov]),
		(defined $args->{title} ? ('--t', $args->{title}) : ()),
		(%{ $labels } ? ('--l', _json_text({ map { _chars($_) => _chars($labels->{$_}) } keys %{ $labels } })) : ()));
	$r->{'output_file'} = $args->{'output_file'};
	_wrote($args->{output_file}, 43, $drawing->{execute} ? 'wrote' : 'will write'); # yellow
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

 # the PDB structures that match a sequence, how well, and what kind they are
 my $pdb = pdb_search(seq => $one);
 print "$pdb->[0]{id} $pdb->[0]{pident}% $pdb->[0]{method} $pdb->[0]{resolution}\n";

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

 # how every two people are related: father-son, brother-sister, ...
 my $kin = relatedness('mom.vcf.gz', 'dad.vcf.gz', 'kid.vcf.gz');
 print $kin->{'kid'}{'dad'}{category}; # by the sample names in the VCFs

 # how well a test VCF's calls agree with a truth VCF's, sample by sample
 my $conc = concordance(truth => 'NA12878.truth.vcf.gz', test => 'NA12878.run42.vcf.gz');
 print $conc->{'NA12878'}{metrics}{genotype_concordance};

 # the 41 HIrisPlex-S markers, as the webtool's CSV wants them
 my $g = hirisplex_genotypes(vcf => 'genome.vcf.gz', file => 'hirisplexs.csv');

 # ... and a prediction from a multinomial logistic model you supply
 my $p = trait_predict(genotypes => $g, model => 'my_model.tsv');
 print $p->{'NA12878'}{eye}{call};

=head1 DESCRIPTION

FASTA I/O in XS, BLAST hit ranking, a search of the PDB's structures for a
sequence, and multiple-sequence-alignment plots and LaTeX tables; polygenic scores from a VCF and a PGS Catalog scoring file; and, from VCFs, each sample's sex, whether two people are parent and child, how every two people are related, how well one VCF's calls agree with another's, and the genotypes of the HIrisPlex-S eye, hair and skin colour panel, with a trait predictor for any multinomial logistic model you supply.

Nothing is exported by default; ask for functions by name or with C<:all>.
Every function dies (via C<Carp::croak>) on bad arguments, with a message that
starts with the name of the function that raised it, as in
C<fasta2hash: couldn't find "four" in x.fa>. C<pdb_search>, C<pgs_score>, C<sex_check>,
C<parent>, C<relatedness>, C<concordance>, C<hirisplex_genotypes>, C<trait_predict>, C<plot_msa>, C<plot_phylo>, C<msa_quality_table> and C<clustal_view_residues>
take C<< name =E<gt> value >> pairs, not a hash ref.

=head2 Alien::Bioinf

Clustal Omega, BLAST+, PLINK 2 and the Python that draws the plots come from
C<Alien::Bioinf>, which installs them module-locally and never uses C<PATH>. It
is recommended rather than required, because it can only be installed where
NCBI builds BLAST+ (Linux and macOS on x86_64 and aarch64, and Windows on
x86_64) and it needs the network; PLINK 2 has no build for Linux on aarch64,
so there it installs the rest. C<fasta2hash>, C<hash2fasta_file>,
C<get_best_alignment_hit>, C<sex_check>, C<parent>, C<relatedness>, C<concordance>,
C<hirisplex_genotypes> and C<trait_predict> work without it; C<pgs_score>, C<plot_msa>,
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

=head2 pdb_search(%args)

Searches the structures of the L<RCSB PDB|https://www.rcsb.org/> for a
sequence, and returns the matching ones, best first, each with how well it
aligns and what kind of structure it is:

 use Cn::Bioinf qw(pdb_search);

 my $hits = pdb_search(seq => 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF...');
 my $best = $hits->[0];
 print "$best->{id}: $best->{pident}% identical over $best->{query_cover}% of the query, ",
     "$best->{method} at $best->{resolution} A, $best->{oligomeric_state}\n";
 # 1BZ1_1: 100% identical over 100% of the query, X-ray at 1.59 A, Hetero 4-mer

 my ($hits, $total) = pdb_search(seq => $seq, max_hits => 10); # and how many RCSB found

It needs the network, and HTTPS, which HTTP::Tiny gets from C<IO::Socket::SSL>
and C<Net::SSLeay>; without them it dies saying which is missing. It asks two
of RCSB's APIs: the L<search API|https://search.rcsb.org/> finds the polymer
entities (one distinct chain of a structure) whose sequence matches, and the
L<data API|https://data.rcsb.org/> says what each one's structure is.

RCSB aligns with MMseqs2, not BLAST, so the numbers are MMseqs2's: its E-value
and bit score are on BLAST's scale but are not what C<blastp> would give for the
same pair. They are returned under the names C<get_best_alignment_hit> gives a
BLAST hsp.

=head3 Arguments

=over

=item * C<seq> (required): the sequence, as a string. Whitespace is removed, lower
case is made upper case, and a final C<*> (a stop) is dropped, so a sequence
pasted from a FASTA file can be given as it is; any other character that is
not a letter dies. RCSB itself refuses a sequence shorter than 10 residues.

=item * C<type>: C<protein> (the default), C<dna> or C<rna>.

=item * C<evalue>: hits with this E-value or above are left out. The default is 0.1,
not BLAST's 10, which in a search of the whole PDB lets in hits expected 10
times over by chance.

=item * C<identity>: hits with less identity than this percentage are left out
(default 0).

=item * C<max_hits>: the most hits returned (default 100). RCSB allows up to 10000.

=item * C<computed>: 0 (the default) searches experimental structures only; 1 also
searches the computed structure models RCSB holds, from AlphaFold DB and
ModelArchive.

=item * C<report>: 1 (the default) prints a table of the hits to STDOUT; 0 returns
them only.

=item * C<timeout>: seconds to wait for each request to RCSB (default 60).

=back

=head3 What it returns

An array ref of hits, best first: by E-value, then bit score, then RCSB's own
order, as C<get_best_alignment_hit> sorts by C<evalue>. Many hits tie: every
structure of the query's own protein aligns the same way. In list context, the
number of hits RCSB found is returned as well, which can be more than
C<max_hits>. Nothing matching is an empty array ref, not an error.

Each hit is a hash ref with every one of these keys, C<undef> or an empty array
ref where RCSB has no value:

What it is:

=over

=item * C<id>: the polymer entity, as C<1BZ1_1>: entry C<1BZ1>, entity C<1>.

=item * C<entry>, C<entity>: the two halves of C<id>. A computed model's entry is like
C<AF_AFP69905F1>.

=item * C<url>: the entry's page at RCSB.

=back

How well it aligns:

=over

=item * C<evalue>, C<bit_score>: MMseqs2's.

=item * C<align_len>: the length of the alignment, gaps included.

=item * C<identity>: the number of aligned columns where the two agree, and C<pident>,
that as a percentage of C<align_len>. These are counted from the aligned
sequences, since RCSB gives the identity only truncated to 3 places.

=item * C<mismatches>, C<gaps> (gap columns), C<gap_opens>.

=item * C<query_from>, C<query_to>, C<query_len> and C<query_cover>, the percentage of
the query in the alignment.

=item * C<hit_from>, C<hit_to>, C<hit_len>: the same for the hit's sequence.

=item * C<qseq>, C<hseq>: the aligned query and hit, with C<-> for gaps.

=back

What kind of structure it is:

=over

=item * C<title>: the entry's title, and C<description>, the entity's name.

=item * C<method>: how the structure was determined: C<X-ray>, C<EM>, C<NMR>,
C<Neutron>, C<Integrative>, C<Multiple methods> or C<Other>; C<undef> for a
computed model.

=item * C<methodology>: C<experimental>, C<computational> or C<integrative>.

=item * C<resolution>: in Å, for X-ray and EM structures.

=item * C<plddt>: a computed model's mean pLDDT, 0 to 100.

=item * C<polymer_type>: C<Protein>, C<DNA>, C<RNA>, C<NA-hybrid> or C<Other>.

=item * C<polymer_composition>: the entry as a whole, such as C<homomeric protein>,
C<heteromeric protein> or C<protein/NA>.

=item * C<oligomeric_state>, C<stoichiometry> and C<symmetry>: the global symmetry of
the first biological assembly, such as C<Hetero 4-mer>, C<A2B2> and C<C2>.

=item * C<chains>: an array ref of the chains that are this entity.

=item * C<organisms>: an array ref of C<{ name, taxonomy_id }>.

=item * C<uniprot>: an array ref of UniProt accessions.

=item * C<pfam>, C<interpro>: array refs of C<{ id, name }>.

=item * C<ligands>: an array ref of C<{ id, name }>, one for each kind of small
molecule in the entry, such as C<HEM>.

=item * C<deposited>, C<released>: dates, as C<YYYY-MM-DD>.

=back

A hit that the data API has no record of (one so new that RCSB's search index
has it first) is warned of, and is returned with its alignment and no structure.

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
 sample  sex   het/hom-ALT  X het  X hom-ALT  X haploid  X missing  Y called
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

=head2 relatedness(@vcf, %args)

Tells how every two people are related, from their genotypes: each pair's
KING-robust kinship and its IBS0, IBS1 and IBS2 counts, as C<parent> takes them
for one pair, and from those a relationship (parent and child, full siblings,
one person twice, second- or third-degree relatives, or unrelated). Each
person's sex is told from X in the same reading, as C<sex_check> tells it, so
that the relationship can be named: C<father-son>, C<mother-daughter>,
C<brother-sister>, C<brothers>, C<sisters>, and so on. Genotypes alone cannot
tell which of a parent and child is the parent, but a third person often can,
and C<relatedness> says which when one does.

The people are either every sample of one VCF, or the one sample of each of
several VCFs, as a sequencing company delivers them:

 use Cn::Bioinf 'relatedness';

 my $family = relatedness('family.vcf.gz');              # every pair of its samples
 my $same   = relatedness('mom.vcf.gz', 'dad.vcf.gz', 'kid.vcf.gz', report => 0);
 my ($pairs, $samples) = relatedness(vcf => ['mom.vcf.gz', 'dad.vcf.gz']);

 print $family->{js}{dc}{category}; # father-son
 print $family->{dc}{js}{category}; # the same hash, so the same

The second prints

 relatedness: mom.vcf.gz, dad.vcf.gz, kid.vcf.gz
 KING-robust kinship over autosomal biallelic SNVs where either carries the ALT, genotypes of DP below 10 not counted; first-degree from 0.177 to 0.354, and of those, parent and child when the IBS0 rate is below 0.012; male when het/hom-ALT on X outside its pseudo-autosomal regions is below 0.6
 sample  sex     X het/hom-ALT  VCF
 mom     female          1.690  mom.vcf.gz
 dad     male            0.206  dad.vcf.gz
 kid     female          1.698  kid.vcf.gz
 sample  sample   sites  IBS0 rate  kinship  relationship
 mom     dad     171033     0.0810   0.0174  unrelated
 mom     kid     158616     0.0020   0.2604  mother-daughter (mom the parent)
 dad     kid     150970     0.0030   0.2587  father-daughter (dad the parent)

for the exomes of the family under C<parent>, and returns the same. Every
pair's counts are those C<parent> gives for it.

The VCFs are decompressed by perl and counted in XS. One VCF of several
samples is read a record at a time, keeping nothing but each pair's counters,
and each pair's counts are taken 64 sites at a time from bit masks of the
samples' genotypes; several VCFs of one sample each are held in memory, 8 bytes
a site, as C<parent> holds two, and opened one at a time. On 2026-10-19 the
three exomes above took 2.5 s and 30 MB, where C<parent> takes 1.6 s for one
pair of them. The 17 samples of the CEPH 1463 pedigree (peddy's
C<ceph1463.peddy.vcf.gz>, 20,000 sites) took 0.27 s. A made-up VCF of 1,000
samples and 20,000 sites took 5.3 s and 890 MB: the counting took 1.7 s and 96
MB, and nearly all the rest is the 499,500 hashes returned.

=over

=item * The VCFs are given first, as C<relatedness($file)> or
C<relatedness($file, $file2, ...)>, before any option, or as C<vcf>. Each is a
C<.vcf>, C<.vcf.gz>, C<.vcf.bgz> or C<.vcf.bz2> file, as for C<sex_check>. One VCF
must have two samples or more; each of several must have one, and no two the
same sample name, which is what each sample is known by.

=item * C<vcf>: a VCF, or an array ref of them, in place of giving them first.

=item * C<min_depth>: a genotype whose C<DP> is below this is not counted, on X as on
the autosomes. The default, 10, is C<parent>'s, for C<parent>'s reasons; 0
counts every genotype.

=item * C<max_ibs0>: the IBS0 rate that a parent and child must be below, which tells
them from full siblings. The default, 0.012, is C<parent>'s.

=item * C<report>: 1 (the default) prints the report to STDOUT; 0 prints nothing.

=back

What is counted is what C<parent> counts, for every pair: the autosomal
biallelic SNVs where either of the two carries the ALT, a site that one own VCF
lacks taken as hom-REF in that person, and a site left out of a pair for the
reasons given under C<parent>. Each person's sex is from SNVs on X outside its
pseudo-autosomal regions, as C<sex_check> tells it with its default C<cutoff>
(0.6), but with C<min_depth> as given here.

A pair's relationship is decided as C<parent> decides it: by KING's kinship
ranges, with the IBS0 rate telling parent and child from full siblings. Its
category names it from the two sexes:

=over

=item * A parent and child is C<father-son>, C<father-daughter>, C<mother-son> or
C<mother-daughter> when it is known which is the parent and both sexes are
known, with C<parent> or C<child> in place of a word whose sex is unknown
(C<parent-son>). When it is not known which is the parent, two men are
C<father-son> and two women C<mother-daughter> all the same, and a man and a
woman are C<parent-child>.

=item * Full siblings are C<brothers>, C<sisters> or C<brother-sister>, or
C<full siblings> when either's sex is unknown.

=item * Any other relationship is its own category: C<duplicate>, C<second-degree>,
C<third-degree> or C<unrelated>.

=back

Which of a parent and child is the parent is told by a third person: x is y's
parent when x has another child who is y's full sibling, or when y has another
parent who is unrelated to x (the two being y's parents). When neither is
said, or both, as in a pedigree that this reading cannot untangle, C<parent> is
undef. On CEPH 1463, three generations of 17 people, C<relatedness> agreed with
the pedigree on all 17 sexes and on 134 of the 136 pairs, including every
father, mother, son, daughter, brother and sister, and which of each parent
and child was the parent. The other two pairs were the grandson NA12893 with
his grandfather and with his grandmother, whose kinships of 0.0839 and 0.0799
are just under KING's 0.0884 for second-degree relatives, and so are
C<third-degree>.

Where C<parent> would die, C<relatedness> decides nothing about that pair and
goes on to the next: its C<relationship> and C<category> are undef and
C<undecided> says why, and C<relatedness> warns once, counting such pairs. That
is when fewer than 1,000 sites are compared, when neither has a hom-REF call
where the other carries the ALT (as in a VCF merged without C<-0>), or when one
has no het call, which leaves C<kinship> undef too. When several VCFs are
given, C<relatedness> also warns once about the pairs whose VCFs look
sequenced over different regions, as C<parent> warns about one.

In scalar context C<relatedness> returns a hash of hashes keyed by each sample
and then by each other sample, so that C<< $r-E<gt>{js}{dc} >> and C<< $r-E<gt>{dc}{js} >> are
the same pair, and the same hash: changing one changes the other. In list
context it returns that and a hash ref of the samples. A pair's hash has these
fields on every return:

=over

=item * C<relationship>: as C<parent> gives it (C<parent-child>, C<full siblings>,
C<duplicate>, C<second-degree>, C<third-degree> or C<unrelated>), or undef when
the pair is undecided.

=item * C<category>: the relationship named by sex, as above, or undef.

=item * C<parent>: the name of the one that is the parent, for a parent and child
where a third person says which; otherwise undef.

=item * C<undecided>: why the pair is undecided, or undef.

=item * C<kinship>: the KING-robust between-family estimate, or undef when either has
no het call. On C<t/data/parent.vcf> it agreed with PLINK 2's
C<--make-king-table> for all 21 pairs.

=item * C<ibs0_rate>: C<ibs0> divided by C<sites>, or undef when there are no sites.

=item * C<sites>: the sites compared, of which C<in_both> were in both VCFs, and
C<only>, a hash of each of the two samples, those in that sample's VCF alone.
With one VCF, every site is C<in_both>.

=item * C<ibs0>, C<ibs1>, C<ibs2>: the sites where the two share no allele, one, or
both; C<het_het>: those where both are het; C<het>, a hash of each of the two
samples: its het calls at the sites compared.

=back

C<parent>'s split of C<ibs0>, and each one's hom-ALT count, are left out, since
each hash of two numbers adds a third to a pair's memory; C<parent> gives them
for any one pair.

The samples' hash has, for each sample, what C<sex_check> returns for it
(C<sex>, C<het_ratio>, C<x_het>, C<x_hom_alt> and the rest), and C<file>, its VCF.

A BCF is not read; convert it with C<bcftools view -Oz>.

=head2 concordance(%args)

Compares the genotype calls of a VCF under test with those of a truth VCF,
sample by sample, and counts each kind of agreement and disagreement, as
C<x.low_ram_local_concordance.pl> did:

 use Cn::Bioinf 'concordance';

 my $r = concordance(truth => 'truth.vcf.gz', test => 'test.vcf.gz');
 print $r->{NA12878}{metrics}{genotype_concordance};
 foreach my $site (@{ $r->{NA12878}{discordances} }) {
     print "$site->{chrom}:$site->{pos} $site->{type}: test $site->{test}, truth $site->{truth}\n";
 }

This is based off of https://link.springer.com/article/10.1186/s12920-022-01199-8

Each of the test's samples is compared with the truth's sample of the same
name, and a sample in only one of the files is left out. When each file has
exactly one sample, the two are compared whatever they are called, so a run
named C<run42> can be compared with its truth C<NA12878>. When no sample can be
paired, C<concordance> dies and names some of each file's samples.

Every record of the test is compared with the truth's record of the same
CHROM, POS, REF and ALT. A C<chr> prefix on CHROM, in any case, is ignored, and
so is ID. When the truth has no such record, its call there is a no-call. A
site that the truth has and the test lacks is not compared at all, so the
test's sites are the ones measured, as for a genotyping panel against a
genome. A site whose ALTs are listed in another order, or split over two
records, is a different site.

The two files have to be on the same genome build, since a site's POS in
GRCh37 (hg19) is not the same place in GRCh38 (hg38). Before either file is
read, C<concordance> compares the lengths that their headers' C<##contig> lines
give, matching contigs by name without C<chr>, as sites are matched. It dies if
any contig that both give a length for has different lengths, and names the
contigs and the build each header gives (its C<assembly>, or else its
C<##reference>). hg19 and b37, whose chromosomes are the same lengths, pass. A
header with no C<##contig> lengths cannot be checked this way. When none of the
test's sites is in the truth, C<concordance> warns that the builds may differ.

On C<t/data/concordance.test.vcf> against C<t/data/concordance.truth.vcf> it
prints

 concordance: test concordance.test.vcf against truth concordance.truth.vcf
 each test site compared with the truth's of the same CHROM, POS, REF and ALT, a site the truth lacks being a no-call there; precision, recall and specificity are of genotypes
 sample  sites  in truth  concordant  discordant  gained  lost  concordance  precision  recall  specificity
 B          46        45          34          10       2     0       0.7674     0.9412  0.8889       0.9375
 A          46        45          11          32       2     1       0.2381     0.4286  0.4286       0.1111
 1 of the test's 3 samples are not in the truth: D

When the paired samples' names differ, a C<truth sample> column follows
C<sample>.

The VCFs are decompressed by perl and compared in XS. The test is read first,
and its sites are kept with their calls for the samples compared, about 100
bytes a site for one sample. Then the truth is read a record at a time and
looked up among them, so memory grows with the test, which is usually the
smaller file. On 2026-10-19, one exome compared against another (274,879 test
sites) took 1.9 s and 211 MB, or 1.7 s and 47 MB with C<< discordances =E<gt> 0 >>,
because most of that memory is the 163,376 discordant and gained sites
returned. A genome compared with itself (3.5 million sites) took 14 s and 362
MB, most of the time being gzip.

=over

=item * C<truth>, C<test>: the two VCFs, each a C<.vcf>, C<.vcf.gz>, C<.vcf.bgz> or
C<.vcf.bz2> file, as for C<sex_check>.

=item * C<min_depth>: a genotype whose C<DP> is below this is a no-call, in either
file. The default, 0, compares every genotype, as the script did. A C<DP>
that is absent or C<.> leaves the genotype compared.

=item * C<discordances>: 1 (the default) lists each discordant, gained and lost site
in C<discordances>; 0 only counts them, and leaves C<discordances> undef.

=item * C<report>: 1 (the default) prints the report to STDOUT; 0 prints nothing.

=back

A call is compared by its alleles, in either order and phased or not, so
C<1/0>, C<0|1> and C<0/1> are the same call. Each call is one of hom-REF
(C<0/0>), het-REF (C<0/1>), hom-ALT (C<1/1>, C<2/2>), het-ALT (C<1/2>),
haploid-REF (C<0>), haploid-ALT (C<1>), or a no-call. A missing call (C<.> or
C<./.>), a partly missing one (C<./1>), one of more than two alleles, and one
where FORMAT has no C<GT> are all no-calls. Each comparison, truth first, is one
of 45 kinds, the script's:

=over

=item * C<concordant_> and the class, when the two calls are the same:
C<concordant_hom_ref>, C<concordant_het_ref>, C<concordant_hom_alt>,
C<concordant_het_alt>, C<concordant_haploid_ref>, C<concordant_haploid_alt> and
C<concordant_no_call>.

=item * C<concordant_diploid_ref_to_haploid_ref> and
C<concordant_haploid_ref_to_diploid_ref>, for C<0/0> against C<0> and C<0>
against C<0/0>; C<concordant_hom_alt_to_diploid_alt> and
C<concordant_haploid_alt_to_hom_alt>, for C<1/1> against C<1> and C<1> against
C<1/1>.

=item * C<lost_site>, when the truth has a call and the test a no-call, and
C<gained_site>, when the test has a call and the truth a no-call.

=item * C<discordant_> and the two classes otherwise, as in
C<discordant_het_ref_to_hom_alt> (C<0/1> in the truth, C<1/1> in the test).
Two different calls of the same class are discordant too, as in
C<discordant_het_ref_to_het_ref> (C<0/1> against C<0/2>).

=back

C<concordance> returns a hash keyed by the test's name for each sample
compared. Each sample's hash has these fields on every return:

=over

=item * C<test_sample>, C<truth_sample>: the sample's name in each file.

=item * C<sites>: the test's records, every one compared; C<in_truth>: how many of them
the truth has.

=item * C<counts>: a hash of the 45 kinds, each with its count, 0 included.

=item * C<metrics>: a hash of the script's metrics, taken from C<counts>:

=over

=item * C<tp_site>: the matched calls that carry an ALT: C<concordant_het_ref>,
C<concordant_het_alt>, C<concordant_hom_alt>, C<concordant_haploid_alt>,
C<concordant_hom_alt_to_diploid_alt> and
C<concordant_haploid_alt_to_hom_alt>. C<tp_genotype> is the same number.

=item * C<fp_site>: C<gained_site>; C<fn_site>: C<lost_site>.

=item * C<tn_site>: the matched reference calls: C<concordant_hom_ref>,
C<concordant_haploid_ref>, C<concordant_diploid_ref_to_haploid_ref> and
C<concordant_haploid_ref_to_diploid_ref>. C<tn_genotype> is the same
number.

=item * C<fp_genotype>: the discordances where the test's call carries an ALT,
whether the truth's is a reference call or another ALT-carrying call;
C<fn_genotype>: those where the truth's call carries an ALT, whether the
test's is a reference call or another ALT-carrying call. A discordance of
two ALT-carrying calls, such as C<0/1> against C<1/1>, is thus counted in
both, since the test's call is wrong and the truth's is missed.

=item * C<genotype_concordance>: the concordant comparisons, less
C<concordant_no_call>, over every comparison but C<concordant_no_call>,
C<gained_site> and C<lost_site>.

=item * C<precision_site>: C<tp_site / (tp_site + fp_site)>; C<recall_site>:
C<tp_site / (tp_site + fn_site)>; C<precision_genotype>:
C<tp_genotype / (tp_genotype + fp_genotype)>; C<recall_genotype>:
C<tp_genotype / (tp_genotype + fn_genotype)>; C<specificity_genotype>:
C<tn_genotype> over C<tn_genotype> and those of C<fp_genotype> where the
truth is a reference call, since only those truths are negatives. Each
ratio is undef when its
denominator is 0, where the script gave the string C<undefined>.

=back

=item * C<discordances>: an array ref of the discordant, gained and lost sites, in the
test's order, each a hash of C<chrom>, C<pos>, C<id>, C<ref> and C<alt> as the
test has them, C<type> (the kind), and C<test> and C<truth>, the two calls as
GT writes them, or C<.> for a no-call; or undef with C<< discordances =E<gt> 0 >>.

=back

Four of the metrics are corrected from the script's. Its C<tp_site> pattern,
C</^conc\w+(?:alt|hom)/>, matched the C<hom> in C<concordant_hom_ref>, so a
hom-REF site was a true positive as well as a true negative, and it missed
C<concordant_het_ref>, the commonest true positive. Its C<tn_site> missed
C<concordant_diploid_ref_to_haploid_ref>, though it took the reverse. And its
C<tn_genotype> was C<concordant_hom_ref> alone, though its false positives
counted a haploid-REF truth called wrong, so C<specificity_genotype> gave no
credit for one called right, as on a male's X. And it counted a discordance
of two ALT-carrying calls in neither C<fp_genotype> nor C<fn_genotype>, so that
a wrong genotype at a variant site changed neither C<precision_genotype> nor
C<recall_genotype>.

C<concordance> departs from the script where the script failed. The script
took ID as part of a site, and dropped C<chr> from the truth's CHROM alone, so
two files that both had it shared no site. It compared GT strings, and so died
on C<1/0>, C<0|1>, C<./1> and calls of three alleles. And it made a no-call
against C<.> C<discordant_no_call_to_no_call>, outside its own 45 kinds. Its XLSX workbook,
JSON files and ACMG59 column are left to the caller, who has everything they
were made from.

C<concordance> warns when the test has no records, and when none of the test's
sites is in the truth, which usually means two genome builds, or two naming
schemes for the chromosomes other than C<chr>.

A BCF is not read; convert it with C<bcftools view -Oz>.

=head2 hirisplex_genotypes(%args)

Reads the 41 markers of the HIrisPlex-S eye, hair and skin colour panel out of
a VCF, counting for each sample how many copies of each marker's model allele
it carries, and can write them as the CSV the HIrisPlex-S webtool uploads:

 use Cn::Bioinf 'hirisplex_genotypes';

 my $g = hirisplex_genotypes(vcf => 'genome.vcf.gz', file => 'hirisplexs.csv');
 print $g->{NA12878}{genotype}{'rs12913832_T'};  # 0, 1, 2, or undef
 print scalar @{ $g->{NA12878}{missing} }, " markers missing\n";

B<< This distribution ships no prediction model, and C<hirisplex_genotypes> makes
no prediction. >> The HIrisPlex-S regression coefficients are not published:
Chaitanya et al. (2018) cite them by reference rather than tabulating them,
the webtool offers no download and says that no part of it may be reproduced
without consent, and part of the knowledge behind it is covered by European
patent EP2195448. What is public, and what this function uses, is the panel
itself -- which markers, in which order, and which allele each model counts --
taken from the webtool's own input form. So C<hirisplex_genotypes> prepares the
input, and the prediction is made either at
L<the webtool|https://hirisplex.erasmusmc.nl/> or by C<trait_predict> from a
model table you hold yourself.

=head3 The panel

C<share/hirisplex_s_panel.tsv> holds the panel: each marker's rs number, gene,
the allele the model counts, its positions on GRCh37 and GRCh38, and whether
it belongs to the 6-marker IrisPlex and 24-marker HIrisPlex subsets.
C<share/make_hirisplex_panel.pl> writes that file from the webtool's form and
Ensembl, and its header comment says how to re-run it; the committed file is
what is used, and the generator is never run by the module.

The webtool asks for a count of the allele on the assay's strand, which is not
always the strand a VCF is written on, and the panel records the forward-strand
allele to count for each marker so that nothing has to be guessed at run time.
Fourteen of the 41 are counted as the complement of the allele the webtool
names. One, ASIP C<rs6119471>, carries the named allele on the forward strand
I<and> still has to be flipped -- the HIrisPlex-S manual singles it out as "the
only SNP that may cause confusion and therefore must be converted (from NCBI's
forward orientation)" -- so a tool that simply matched the named allele would
silently count the wrong one. MC1R C<rs312262906> is not a SNP but the 86insA
frameshift, merged in dbSNP into C<rs796296176>; it is counted as the allele
longer than REF, and matching by rs number looks for both numbers.

=head3 Arguments

=over

=item * C<vcf> (required): a C<.vcf>, C<.vcf.gz>, C<.vcf.bgz> or C<.vcf.bz2> file. As
C<sex_check> does, it may be given first and on its own:
C<hirisplex_genotypes('genome.vcf.gz')>.

=item * C<match>: C<position> (the default) matches CHROM and POS against the panel's
coordinates; C<id> matches the ID column against the markers' rs numbers,
which finds nothing in a VCF whose ID column is C<.>.

=item * C<build>: C<hg19>, C<b37>, C<hg38> or C<b38>, saying which of the panel's two sets
of coordinates to match on. By default it is read from the VCF's C<##contig>
lengths, and C<hirisplex_genotypes> dies asking for it when the header gives
no length it knows, or lengths that are neither build's.

=item * C<min_depth>: a genotype whose FORMAT C<DP> is below this is not counted.
Defaults to 7, C<sex_check>'s default; 0 counts every call whatever its depth.
A genotype with no C<DP> is always counted.

=item * C<file>: write the webtool's CSV here. Its header is that of the webtool's own
template, column for column, and a marker with no genotype is written C<NA>.

=item * C<report>: 1 (the default) prints a table of the samples to STDOUT; 0 returns
them only.

=back

=head3 What it returns

A hash ref keyed by sample name, each holding:

=over

=item * C<genotype>: a hash ref keyed as the webtool's CSV heads its columns,
C<< E<lt>rs numberE<gt>_E<lt>input alleleE<gt> >>, holding 0, 1 or 2 copies of that allele, or
C<undef> where the sample has no usable call.

=item * C<missing>: an array ref of the rs numbers with no genotype, in the panel's
order.

=item * C<markers> and C<called>: how many markers the panel holds, and how many this
sample has a genotype for.

=item * C<eye>, C<hair>, C<skin>: 1 or 0, whether the HIrisPlex-S rules allow that
prediction from these markers at all. The manual refuses an eye prediction
without HERC2 C<rs12913832>, a hair prediction when all 11 MC1R markers are
missing, and both when HERC2, SLC45A2 and IRF4 are all missing, because the
loss in AUC is too great. It sets no such rule for skin, so C<skin> is 0 only
when nothing at all was called. C<trait_predict> honours these.

=back

A genotype is counted only when the call is diploid and complete: a C<./.>, a
C<./1>, a haploid call (the panel is autosomal) or one below C<min_depth> leaves
the marker without a genotype. A site split into one record per ALT is counted
from whichever record carries the counted allele, in either order, and a record
at a marker's position whose REF is longer than one base is some other variant
and is passed over.

C<hirisplex_genotypes> warns when the file holds none of the 41 markers, which
usually means the wrong build or an ID column of C<.>.

Reading the markers is done in XS, in one pass, keeping one byte per sample per
marker, so a whole-genome VCF costs no more memory than a small one. A
275,879-record exome of 22 MB gzipped took 0.86 s (perl 5.44.0, 2026-10-20).

Most of the panel is intronic or intergenic, so an exome carries little of it:
that same exome held 9 of the 41 markers. A whole-genome VCF, a genotyping
array that covers the panel, or the SNaPshot assay the panel was designed for
is what the prediction expects, and C<missing> and C<called> are returned so that
a prediction made from a handful of markers can be seen for what it is.

=head2 trait_predict(%args)

Applies a multinomial logistic regression -- which is what the HIrisPlex-S,
HIrisPlex and IrisPlex models are -- to genotypes, and returns a probability
for each category of each trait:

 use Cn::Bioinf qw(hirisplex_genotypes trait_predict);

 my $g = hirisplex_genotypes(vcf => 'genome.vcf.gz', report => 0);
 my $p = trait_predict(genotypes => $g, model => 'my_model.tsv');
 print $p->{NA12878}{eye}{call};          # the most probable category
 print $p->{NA12878}{eye}{p}{blue};       # its probability

The model is yours to supply: see C<hirisplex_genotypes> above for why this
distribution ships none. Any model written in the form below works, so this is
not limited to pigmentation.

=head3 The model table

A tab-separated file (optionally gzipped), or a data frame as
C<Stats::LikeR>'s C<read_table> returns, with four columns:

 trait    category    marker    beta
 eye    blue    (intercept)    3.94
 eye    blue    rs12913832_T    -4.81
 eye    brown    (intercept)    0

C<marker> is a key of the C<genotype> hash C<hirisplex_genotypes> returns, which
is C<< E<lt>rs numberE<gt>_E<lt>input alleleE<gt> >>, or the literal C<(intercept)> for the
category's intercept. A category whose only row is an intercept of 0 is the
reference category, which is how a model of I<k> categories is written with
I<k> - 1 sets of betas; every category that should appear in the result needs at
least one row. Traits and their categories come out in the order the table
first names them.

For each category the linear predictor is its intercept plus the sum of its
betas times the genotypes, and the probabilities are C<exp> of those over their
sum. The largest is taken out before C<exp>, so a linear predictor large enough
to overflow gives 1 rather than a NaN.

=head3 Arguments

=over

=item * C<genotypes> (required): what C<hirisplex_genotypes> returned, or a plain hash
ref of sample => { marker => 0, 1 or 2 }.

=item * C<model> (required): the table above, as a file name or a data frame.

=item * C<force>: 0 (the default) makes no prediction for a trait that the genotypes
say cannot be predicted -- that is, where C<hirisplex_genotypes> set C<eye>,
C<hair> or C<skin> to 0; 1 predicts it anyway. A trait the genotypes say
nothing about is always predicted.

=item * C<report>: 1 (the default) prints a table of the samples to STDOUT; 0 returns
them only.

=back

=head3 What it returns

A hash ref keyed by sample and then by trait, each holding:

=over

=item * C<p>: a hash ref of category => probability, which sum to 1; empty when the
prediction was refused.

=item * C<call>: the most probable category, or C<undef> when it was refused. A tie is
broken by the category's name, so the answer does not depend on hash order.

=item * C<markers> and C<called>: how many markers the trait's categories use between
them, and how many this sample has a genotype for.

=item * C<missing>: an array ref of the markers the model names and the sample has no
genotype for.

=item * C<refused>: C<undef>, or why no prediction was made.

=back

A marker the sample has no genotype for drops out of the sum, which is the same
as counting it 0 copies of the counted allele. That is a choice, not something
the HIrisPlex-S authors specify, and it biases the result towards whatever 0
copies means for that marker; C<missing> and C<called> are returned so that a
prediction made from few markers can be reported as such. The webtool reports a
loss in AUC for exactly this reason.

=head2 plot_msa(%args)

Aligns sequences with Clustal Omega and draws the alignment. Returns a hash ref
of the files made: C<filename>, C<msa_file>, C<python_file> (see C<fh> and
C<execute> below), and C<tree_file> if it was given. Once the image is written it
prints C<wrote> and its file name to STDOUT, the name in black on yellow when
STDOUT is a terminal; with C<< execute =E<gt> 0 >> it prints C<will write> instead.

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

=item * C<fh>, C<execute>: as for Matplotlib::Simple's C<plt>. The image is drawn by a
Python script, and C<fh> is a C<File::Temp> object to append that script to;
without it, the script goes in a new temporary file. C<execute> is 1 (the
default) to run the script now, under Alien::Bioinf's Python, or 0 only to
write it. So several calls with one C<fh> and C<< execute =E<gt> 0 >>, then a last one
with C<< execute =E<gt> 1 >>, draw every image in one run of Python. The C<fh> may be
one that Matplotlib::Simple's C<plt> also writes to, since the Python each
function appends stands alone. Clustal Omega runs either way; only the
drawing waits.

C<python_file> is the script's file name: C<fh>'s, or the temporary file's,
which is kept after perl exits only when C<execute> is 0. The script carries
the alignment and everything else it reads, so it can be run later, with
C<< perl -MAlien::Bioinf -e 'print Alien::Bioinf-E<gt>python' >> giving the Python to
run it with, after perl's own temporary files have gone.

=back

With fewer than two sequences it warns and returns an empty hash ref.

=head2 plot_phylo(%args)

Draws a guide tree, with Biopython's C<Bio.Phylo>, from a FASTA (aligned with
Clustal Omega first) or from a newick file that C<plot_msa> kept. Returns a
hash ref of the files made or used: C<output_file>, C<python_file>, and
C<tree_file> and C<msa_file> as below.

 plot_phylo(fasta => 't/data/DEG20010421.fa');   # writes phylo.svg
 plot_phylo(
     tree_file   => 'DEG20010421.newick',
     output_file => 'DEG20010421.tree.png',
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

=item * C<fh>, C<execute>: as for C<plot_msa>, which can share an C<fh> with this; the
script carries the newick tree it draws. C<python_file> is the script's file
name.

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
