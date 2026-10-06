package Bioinf::Basic;
# ABSTRACT: FASTA I/O in XS, BLAST hit ranking, and alignment plots and tables
#
# Taken from ~/Scripts/bioinf.pm (as it stood on 2026-08-19): fasta2hash,
# get_best_alignment_hit, hash2fasta_file, msa_quality_table, msa_phylo_plot and
# clustal_view_residues. fasta2hash and hash2fasta_file are now XS (Basic.xs).
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
use Cwd 'getcwd';
use Sys::Hostname 'hostname';
use Exporter 'import';
use XSLoader;
our $VERSION = '0.01';
our @EXPORT_OK = qw(clustal_view_residues fasta2hash get_best_alignment_hit hash2fasta_file msa_quality_table plot_msa plot_phylo);
our %EXPORT_TAGS = (all => \@EXPORT_OK);
XSLoader::load('Bioinf::Basic', $VERSION);
# The calling script's name, for _creator. It is looked up at load time, as
# FindBin always does, but FindBin croaks when $0 is not a file on disk (after
# $0 is assigned to, or under an embedder), and that must not stop this module
# loading, so $0 itself stands in then.
my $RealScript = eval { require FindBin; $FindBin::RealScript } // $0;

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
	getcwd() . "/$RealScript called using \"$sub\" in " . __FILE__ . " version $VERSION with Alien::Bioinf $Alien::Bioinf::VERSION"
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

sub _run {
	my @cmd = @_;
	system(@cmd) == 0 or croak "_run: \"@cmd\" failed: " . ($? == -1 ? $! : 'exit ' . ($? >> 8));
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
# Creator and the JSON options are full of them. A string that perl holds as
# bytes is taken to be UTF-8 if it decodes as such, which is how the command
# line used to pass it on; encode_json alone would read such bytes as Latin-1.
sub _python {
	my ($sub, @args) = @_;
	require JSON::MaybeXS;
	my @chars = map { my $s = "$_"; utf8::decode($s) unless utf8::is_utf8($s); $s } @args;
	my $argfile = _tmp('.json');
	open my $fh, '>:raw', $argfile or croak "_python: can't write $argfile: $!";
	print {$fh} JSON::MaybeXS::encode_json(\@chars);
	close $fh or croak "_python: can't write $argfile: $!";
	_run(_alien($sub)->python, _script(), '--argfile', $argfile);
}

# $data as JSON text: characters, not UTF-8 bytes, since it is one argument
# that _python encodes along with the rest. encode_json here would read a byte
# string's UTF-8 as Latin-1, so that _python's decoding came too late.
sub _json_text {
	my ($data) = @_;
	require JSON::MaybeXS;
	JSON::MaybeXS->new(utf8 => 0, canonical => 1)->encode($data);
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
	File::ShareDir::dist_file('Bioinf-Basic', 'msa_plot.py');
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
	croak "_open_fasta: \"$file\" doesn't exist or isn't a readable file" unless defined $file && -f $file && -r _;
	my $fh;
	if ($file =~ /\.gz\z/) {
		open $fh, '-|', 'gzip', '-dc', $file or croak "_open_fasta: can't run gzip on $file: $!";
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
	croak "get_best_alignment_hit: \"$json_file\" doesn't exist or isn't a readable file" unless defined $json_file && -f $json_file && -r _;
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
			croak "get_best_alignment_hit: a hit has more than 1 description" if @{ $hit->{description} } > 1;
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
			push @hits, { %{ $best },
				accession => $hit->{description}[0]{accession},
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
sub _ungapped_fasta {
	my ($fasta) = @_;
	my ($seqs, @names);
	if (ref $fasta eq 'HASH') {
		$seqs = $fasta;
		@names = sort keys %{ $seqs };
	} elsif (ref $fasta eq '') {
		my $order;
		($seqs, $order) = _fasta_ordered($fasta);
		@names = @{ $order };
	} else {
		croak "_ungapped_fasta: \"fasta\" must be a FASTA file name or a hash ref of name => sequence, not a " . ref($fasta) . ' ref';
	}
	my %ungapped = map { (my $s = $seqs->{$_}) =~ tr/-.//d; $_ => $s } @names;
	(hash2fasta_file(\%ungapped, _tmp('.fa'), \@names), \@names);
}

# Aligns "fasta" with clustalo for $sub. The alignment goes to "msa_file" and
# the guide tree to "tree_file" when those are given, and to temporary files
# otherwise; the tree is not asked for at all unless it is kept or $want_tree.
# Returns the hash ref of files to report (the alignment always, the tree only
# when it was named), the tree's file, the sequence names in input order
# (sorted for a hash), and an array ref of the clustalo command as it was run,
# with the ungapped FASTA it read as its "--in"; or nothing, after a warning,
# for fewer than 2 sequences.
sub _align {
	my ($args, $sub, $want_tree) = @_;
	my ($in, $names) = _ungapped_fasta($args->{fasta});
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
	(\%r, $tree, $names, \@cmd);
}

sub plot_msa {
	my $sub = 'plot_msa';
	my $args = _args(\@_, $sub, [qw(fasta filename)], [qw(active_site_aa clustal_args labels msa_file order query threads title
		tree_file xlabel ylabel)]);
	my $sites = $args->{'active_site_aa'};
	if (defined $sites) {
		croak "$sub: \"active_site_aa\" must be a hash ref of label => residue number" unless ref $sites eq 'HASH';
		croak "$sub: needs \"query\" to place \"active_site_aa\"" unless defined $args->{query};
		_residue_number($sub, "active site $_", $sites->{$_}) foreach sort keys %{ $sites };
	}
	my ($r, undef, $names, $cmd) = _align($args, $sub, 0) or return {};
	my $aln = fasta2hash($r->{'msa_file'});
	my @order = @{ $args->{order} // $names };
	my @unknown = grep { !defined $aln->{$_} } @order;
	croak "$sub: \"order\" names sequences that aren't in the alignment: @unknown" if @unknown;
	my $labels = $args->{labels} // {};
	# Two sequences shown under one label would be one row of the image, fed
	# two deflines.
	my %by_label;
	push @{ $by_label{ $labels->{$_} // $_ } }, $_ foreach @order;
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
		@segments = ('--s', _json_text(\%sites));
	}
	my %shown = map { ($labels->{$_} // $_) => $aln->{$_} } @order;
	my $fa = hash2fasta_file(\%shown, _tmp('.fa'), [map { $labels->{$_} // $_ } @order], 0);
	_python($sub, '--f', $fa, '--o', $args->{filename}, '--c', _creator($sub), '--quiet',
		'--p', _json_text([_align_provenance($args, $sub, $cmd, $names, $r->{'msa_file'}), _kept($args, qw(msa_file tree_file))]),
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
	'from ' . _fasta_id($args->{fasta}, scalar @{ $names }) . ', written with its gaps stripped to ' . _file_id($cmd->[2])
		. ', which Clustal Omega ' . $alien->versions->{clustalo} . " from Alien::Bioinf $Alien::Bioinf::VERSION aligned into "
		. _file_id($msa) . ', run as: ' . _cmd_text(@{ $cmd });
}

sub plot_phylo {
	my $sub = 'plot_phylo';
	my $args = _args(\@_, $sub, [], [qw(clustal_args fasta labels msa_file output_file threads title tree_file)]);
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
	_python($sub, '--tree', $tree, '--o', $args->{'output_file'}, '--c', _creator($sub), '--quiet', '--p', _json_text(\@prov),
		(defined $args->{title} ? ('--t', $args->{title}) : ()),
		(%{ $labels } ? ('--l', _json_text($labels)) : ()));
	$r->{'output_file'} = $args->{'output_file'};
	_wrote($args->{output_file}, 43); # yellow
	$r;
}

sub msa_quality_table {
	my $sub = 'msa_quality_table';
	my $args = _args(\@_, $sub, ['filename'], [qw(alignment_json cb_label cb_max cb_min cblogscale default_undefined fasta
		logscale_add metric msa_file normalize order show_numbers title unaligned_fa)]);
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
	print {$fh} JSON::MaybeXS::encode_json({
		cells => \@cells, cols => \@col_labels, rows => \@row_labels, vmin => $lo + 0, vmax => $hi + 0,
		log => $log ? JSON::MaybeXS::true() : JSON::MaybeXS::false(),
		numbers => $args->{'show_numbers'} ? JSON::MaybeXS::true() : JSON::MaybeXS::false(),
		title => $args->{title} // '', cblabel => $args->{cb_label},
	});
	close $fh or croak "$sub: can't write $table: $!";
	_python($sub, '--table', $table, '--o', $args->{filename}, '--c', _creator($sub), '--p', _json_text(["showing each pair's \"$metric\"", $prov]));
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
	my $color = $args->{'color_residues'} // {};
	croak "$sub: \"color_residues\" must be a hash ref" unless ref $color eq 'HASH';
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
		my ($r) = _align({ %clustal, fasta => $args->{'msa_file'} }, $sub, 0);
		$data = fasta2hash($r->{'msa_file'});
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

 use Bioinf::Basic ':all';

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

=head1 DESCRIPTION

FASTA I/O in XS, BLAST hit ranking, and multiple-sequence-alignment plots and
LaTeX tables, taken from the maintainer's C<bioinf.pm>.

Nothing is exported by default; ask for functions by name or with C<:all>.
Every function dies (via C<Carp::croak>) on bad arguments, with a message that
starts with the name of the function that raised it, as in
C<fasta2hash: couldn't find "four" in x.fa>. C<plot_msa>, C<plot_phylo>,
C<msa_quality_table> and C<clustal_view_residues> take C<< name =E<gt> value >> pairs,
not a hash ref.

=head2 Alien::Bioinf

Clustal Omega, BLAST+ and the Python that draws the plots come from
C<Alien::Bioinf>, which installs them module-locally and never uses C<PATH>. It
is recommended rather than required, because it can only be installed where
NCBI builds BLAST+ (Linux and macOS on x86_64 and aarch64, and Windows on
x86_64) and it needs the network. C<fasta2hash>, C<hash2fasta_file> and
C<get_best_alignment_hit> work without it; C<plot_msa>, C<plot_phylo>,
C<msa_quality_table>, and C<clustal_view_residues> given an unaligned file, die
saying that they need it. From this repository, install it first:

 cpanm ./alien/Alien-Bioinf
 cpanm .

That downloads nothing it can reuse: archives are cached (C<ALIEN_BIOINF_CACHE>,
default C<~/.cache/alien-bioinf>), a clustalo or BLAST+ already on the machine at
the newest version is copied instead (C<ALIEN_BIOINF_CLUSTALO>,
C<ALIEN_BIOINF_BLAST>, or C<PATH>), and pip installs only what the base Python
lacks. To check for and apply updates:

 perl -MAlien::Bioinf -MData::Dumper -e 'print Dumper(Alien::Bioinf->check_updates)'
 perl -MAlien::Bioinf -e 'print "$_\n" for Alien::Bioinf->update'

=head2 Provenance in the images

Every PNG, SVG, PDF, PS or EPS image these functions draw carries its whole
provenance on one line, as its C<Creator> metadata and nothing else, just as
Matplotlib::Simple writes it. The line names the calling script (as the
working directory plus the script's name), the function, this file and its
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
read, the alignment Clustal Omega wrote, the Clustal Omega version, the
exact command, and the C<msa_file> and C<tree_file> kept, if any.

=item * C<plot_phylo>: the same when it aligns, with the guide tree it drew;
otherwise the newick file it drew. Then how many negative branch lengths
were drawn as 0, if any; the tips, with the label each was shown as; and the
whole newick tree that was drawn.

=item * C<msa_quality_table>: the metric shown, and either the BLAST report it read
or the FASTA file, its ungapped copy, the report blastp wrote, the blastp
version, the exact command, and the C<alignment_json> kept, if any.

=back

For example:

 /home/me/work/run.pl called using "plot_phylo" in /.../Bioinf/Basic.pm
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

Reads a FASTA file (gzip-compressed if its name ends in C<.gz>). Returns a hash
ref of defline (without the C<< E<gt> >>) => sequence, or, with C<$key>, just the
sequence of that defline, reading no further than the record after it. A
defline that appears twice is warned about and its sequences concatenated;
with C<$key>, only a repeat of C<$key> is looked for. Line endings may be C<\n>
or C<\r\n>. A C<.gz> that gzip cannot read to its end, such as a truncated one,
dies rather than returning the part that was read.

 my $h = fasta2hash('DEG20010421.fa');

=head2 hash2fasta_file($hash, $filename, $order, $width)

Writes C<$hash> as FASTA: the keys in C<@$order> (default: sorted), sequences
wrapped at C<$width> columns (default 80; 0 for one line each). Returns
C<$filename>.

=head2 get_best_alignment_hit($json_file, $sort_criterion)

For a BLAST C<-outfmt 15> JSON report, a hash ref of query title => array ref
of that query's hits, best first. Each hit is its best hsp plus C<accession>,
C<hit_len>, C<id> and C<title>. C<$sort_criterion> is an hsp field (default
C<evalue>): C<align_len>, C<bit_score>, C<evalue>, C<gaps>, C<identity>,
C<positive> or C<score>. Ties are broken on the bit score, then on BLAST's
order. A report with two queries of one title dies, since only one of them
could be returned.

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
C<\input> into a document.

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
