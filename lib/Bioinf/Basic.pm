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
# that Python all come from Alien::Bioinf, never from PATH.
require 5.010;
use strict;
use warnings;
use Carp qw(croak carp);
use Cwd 'getcwd';
use FindBin '$RealScript';
use Sys::Hostname 'hostname';
use Exporter 'import';
use XSLoader;
our $VERSION = '0.01';
our @EXPORT_OK = qw(clustal_view_residues fasta2hash get_best_alignment_hit hash2fasta_file msa_quality_table plot_msa plot_phylo);
our @EXPORT = @EXPORT_OK;
our %EXPORT_TAGS = (all => \@EXPORT_OK);
XSLoader::load('Bioinf::Basic', $VERSION);

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

# The "Creator" written into an image's metadata, after the one
# Matplotlib::Simple writes:
# which script, run where, called which sub of which version of this module,
# on which computer (its hostname, OS and perl).
sub _creator {
	my ($sub) = @_;
	getcwd() . "/$RealScript called using \"$sub\" in " . __FILE__ . " version $VERSION on host " . hostname()
		. " ($^O, perl " . sprintf('%vd', $^V) . ')';
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
	# Not checked: with $key given, reading stops early, and a gzip pipe
	# closed then reports the SIGPIPE it got as a failure.
	close $fh;
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
		$alignments{ $report->{report}{results}{search}{query_title} } = [ @hits[@order] ];
	}
	\%alignments;
}

# ---- alignment plots and tabless

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

# Aligns "fasta" with clustalo for $sub. The alignment goes to "msa.file" and
# the guide tree to "tree.file" when those are given, and to temporary files
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
	my %r = ('msa.file' => $args->{'msa.file'} // _tmp('.aln.fa'));
	my $tree = $args->{'tree.file'} // ($want_tree ? _tmp('.newick') : undef);
	$r{'tree.file'} = $tree if defined $args->{'tree.file'};
	require Alien::Bioinf;
	my @cmd = (Alien::Bioinf->clustalo, '--in', $in, '--out', $r{'msa.file'}, '--outfmt=fa', '--force',
		'--threads', $args->{threads} // 1, (defined $tree ? "--guidetree-out=$tree" : ()), @{ $args->{'clustal.args'} // [] });
	_run(@cmd);
	(\%r, $tree, $names, \@cmd);
}

sub plot_msa {
	my $sub = 'plot_msa';
	my $args = _args(\@_, $sub, [qw(fasta filename)], [qw(active.site.aa clustal.args labels msa.file order query threads title
		tree.file xlabel ylabel)]);
	croak "$sub: needs \"query\" to place \"active.site.aa\"" if defined $args->{'active.site.aa'} && !defined $args->{query};
	my ($r, undef, $names) = _align($args, $sub, 0) or return {};
	my $aln = fasta2hash($r->{'msa.file'});
	my @order = @{ $args->{order} // $names };
	my @unknown = grep { !defined $aln->{$_} } @order;
	croak "$sub: \"order\" names sequences that aren't in the alignment: @unknown" if @unknown;
	my @segments;
	if (defined $args->{'active.site.aa'}) {
		my $q = $aln->{ $args->{query} } // croak "$sub: the query \"$args->{query}\" isn't in the alignment";
		# 1-based residue number => 0-based alignment column
		my (@col, %sites);
		while ($q =~ /[^-]/g) {
			push @col, pos($q) - 1;
		}
		foreach my $label (sort keys %{ $args->{'active.site.aa'} }) {
			my $n = $args->{'active.site.aa'}{$label};
			$sites{$label} = $col[$n - 1] // croak "$sub: active site $label ($n) is past the end of \"$args->{query}\", which has " . scalar(@col) . ' residues';
		}
		require JSON::MaybeXS;
		@segments = ('--s', JSON::MaybeXS::encode_json(\%sites));
	}
	my $labels = $args->{labels} // {};
	my %shown = map { ($labels->{$_} // $_) => $aln->{$_} } @order;
	my $fa = hash2fasta_file(\%shown, _tmp('.fa'), [map { $labels->{$_} // $_ } @order], 0);
	require Alien::Bioinf;
	_run(Alien::Bioinf->python, _script(), '--f', $fa, '--o', $args->{filename}, '--c', _creator($sub), '--quiet',
		(defined $args->{title} ? ('--t', $args->{title}) : ()),
		(defined $args->{xlabel} ? ('--x', $args->{xlabel}) : ()),
		(defined $args->{ylabel} ? ('--y', $args->{ylabel}) : ()), @segments);
	# Black (30) on yellow (43), then reset (0), written out rather than loaded
	# from Term::ANSIColor; the reset comes before the newline so the colour
	# doesn't run on into the next line.
	print STDOUT "\e[30;43mwrote $args->{filename}\e[0m\n";
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

# Everything plot_phylo knows about where its tree came from, as the Dublin
# Core fields msa_plot.py writes into the image beside the Creator: what was
# read (Source), how the tree was made (Description), the files kept with it
# (Relation), and the Perl-side software that took part (Contributor). $cmd,
# $names and $msa are what _align returned, and undef when nothing was aligned.
# msa_plot.py adds the date, the newick it drew and its digest, the tip names,
# and its own Python, matplotlib, Biopython and NumPy.
sub _phylo_provenance {
	my ($args, $tree, $cmd, $names, $msa) = @_;
	require Alien::Bioinf;
	require File::Spec;
	my $v = Alien::Bioinf->versions;
	my %p = (Title => $args->{title} // 'Phylogenetic tree');
	my @contrib = ("Bioinf::Basic $VERSION (" . __FILE__ . ')', "Alien::Bioinf $Alien::Bioinf::VERSION",
		"perl " . sprintf('%vd', $^V) . " ($^X)");
	if (defined $cmd) {
		my $fasta = $args->{fasta};
		# $cmd->[2] is clustalo's "--in", the ungapped copy _ungapped_fasta wrote
		$p{Source} = (ref $fasta ? 'a hash ref of ' . scalar(@{ $names }) . ' sequences' : 'the FASTA file ' . _file_id($fasta))
			. ', written with its gaps stripped to ' . _file_id($cmd->[2]) . ', which Clustal Omega aligned into '
			. _file_id($msa);
		$p{Description} = "The guide tree Clustal Omega $v->{clustalo} wrote while aligning "
			. scalar(@{ $names }) . ' sequences, run as: ' . join(' ', map { /[^\w\/.,:=+-]/ ? "'$_'" : $_ } @{ $cmd }) . '.';
		unshift @contrib, "Clustal Omega $v->{clustalo} ($cmd->[0])";
		my @kept = map { defined $args->{$_} ? "$_ " . File::Spec->rel2abs($args->{$_}) : () } qw(msa.file tree.file);
		$p{Relation} = 'kept with it: ' . join('; ', @kept) if @kept;
	} else {
		$p{Source} = 'the newick file ' . _file_id($tree);
		$p{Description} = 'Drawn from an existing newick file; no alignment was made.';
	}
	$p{Contributor} = \@contrib;
	\%p;
}

sub plot_phylo {
	my $sub = 'plot_phylo';
	my $args = _args(\@_, $sub, [], [qw(clustal.args fasta labels msa.file output.file threads title tree.file)]);
	$args->{'output.file'} //= 'phylo.svg';
	my ($r, $tree, $names, $cmd);
	if (defined $args->{fasta}) {
		($r, $tree, $names, $cmd) = _align($args, $sub, 1) or return {};
	} else {
		# no alignment to make, so "tree.file" is the tree to draw, not where to keep one
		$tree = $args->{'tree.file'} // croak "$sub: needs \"fasta\" to align, or \"tree.file\" to draw";
		my @moot = grep { defined $args->{$_} } qw(clustal.args msa.file threads);
		croak "$sub: was given " . join(', ', map { "\"$_\"" } @moot) . ' but no "fasta" to align' if @moot;
		croak "$sub: \"$tree\" doesn't exist or isn't a readable file" unless -f $tree && -r _;
		$r = { 'tree.file' => $tree };
	}
	my $labels = $args->{labels} // {};
	require Alien::Bioinf;
	require JSON::MaybeXS;
	_run(Alien::Bioinf->python, _script(), '--tree', $tree, '--o', $args->{'output.file'}, '--c', _creator($sub),
		'--meta', JSON::MaybeXS::encode_json(_phylo_provenance($args, $tree, $cmd, $names, $r->{'msa.file'})),
		(defined $args->{title} ? ('--t', $args->{title}) : ()),
		(%{ $labels } ? ('--l', JSON::MaybeXS::encode_json($labels)) : ()));
	$r->{'output.file'} = $args->{'output.file'};
	$r;
}

sub msa_quality_table {
	my $sub = 'msa_quality_table';
	my $args = _args(\@_, $sub, ['filename'], [qw(alignment.json cb_label cb_max cb_min cblogscale default_undefined fasta
		logscale.add metric msa.file normalize order show.numbers title unaligned.fa)]);
	my %metric = map { $_ => 1 } qw(num bit_score score evalue identity positive align_len);
	my $metric = $args->{metric} // 'score';
	croak "$sub: \"$metric\" isn't one of the metrics: " . join(', ', sort keys %metric) unless $metric{$metric};
	# "unaligned.fa" is the old name for "fasta"
	croak "$sub: was given both \"fasta\" and \"unaligned.fa\", which are the same thing" if defined $args->{fasta} && defined $args->{'unaligned.fa'};
	my $fasta = $args->{fasta} // $args->{'unaligned.fa'};
	# All-against-all blastp of the sequences: given as the parsed report, read
	# from "alignment.json", or -- when that file does not exist yet, or none is
	# named -- made by running blastp on "fasta", and kept in "alignment.json"
	# for next time if that is named.
	my $aj = $args->{'alignment.json'};
	my $blast;
	if (ref $aj eq 'HASH') {
		$blast = $aj;
	} elsif (defined $aj && -f $aj) {
		$blast = _json_file($aj);
	} else {
		croak "$sub: needs \"fasta\" to align, or an existing \"alignment.json\"" . (defined $aj ? " ($aj doesn't exist yet)" : '')
			unless defined $fasta;
		my ($in) = _ungapped_fasta($fasta);
		$aj //= _tmp('.json');
		require Alien::Bioinf;
		_run(Alien::Bioinf->blast('blastp'), '-query', $in, '-subject', $in, '-out', $aj, '-outfmt', 15);
		$blast = _json_file($aj);
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
	my $add = $args->{'logscale.add'} // 0;
	my $norm = ($args->{normalize} // 0) > 0;
	my $log = $args->{cblogscale};
	my (@cells, $lo, $hi, $lo_positive);
	foreach my $i (0 .. $#order) {
		foreach my $j (0 .. $#order) {
			# a pair with no hit stays out of the table, rather than claiming a
			# score of 0, and is drawn grey
			my $value = $data{ $order[$i] }{ $order[$j] } // $args->{default_undefined};
			if (defined $value) {
				$value += $add;
				$value /= $max if $norm;
				$lo = $value if !defined $lo || $value < $lo;
				$hi = $value if !defined $hi || $value > $hi;
				$lo_positive = $value if $value > 0 && (!defined $lo_positive || $value < $lo_positive);
			}
			$cells[$i][$j] = $value;
		}
	}
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
		numbers => $args->{'show.numbers'} ? JSON::MaybeXS::true() : JSON::MaybeXS::false(),
		title => $args->{title} // '', cblabel => $args->{cb_label},
	});
	close $fh or croak "$sub: can't write $table: $!";
	require Alien::Bioinf;
	_run(Alien::Bioinf->python, _script(), '--table', $table, '--o', $args->{filename}, '--c', _creator($sub));
	$args->{filename};
}

sub clustal_view_residues {
	my $sub = 'clustal_view_residues';
	# msa.file: a FASTA file, aligned (such as plot_msa's) or not; color.residues:
	# {protein => {1-based residue number => colour}}, the colour an xcolor
	# name or an [r, g, b] ref; order: proteins top to bottom; row.width:
	# columns per block (100); split: blocks per LaTeX table (4); track: a
	# protein whose coloured residue numbers get a row of their own; threads and
	# clustal.args: for clustalo, if msa.file has to be aligned.
	my $args = _args(\@_, $sub, ['msa.file', 'output.tex.file'], [qw(caption clustal.args color.residues label
		order row.width split table.text.size threads track)]);
	my $color = $args->{'color.residues'} // {};
	croak "$sub: \"color.residues\" must be a hash ref" unless ref $color eq 'HASH';
	croak "$sub: \"msa.file\" $args->{'msa.file'} doesn't exist" unless -e $args->{'msa.file'};
	my $data = fasta2hash($args->{'msa.file'});
	croak "$sub: has no sequences to show in $args->{'msa.file'}" unless %{ $data };
	my %len = map { length $_ => 1 } values %{ $data };
	if (keys %len > 1) {
		# Sequences of different lengths can't be an alignment, so align them,
		# into a temporary file: "msa.file" is the input, and is never written.
		my %clustal = map { $_ => $args->{$_} } grep { defined $args->{$_} } qw(clustal.args threads);
		my ($r) = _align({ %clustal, fasta => $args->{'msa.file'} }, $sub, 0);
		$data = fasta2hash($r->{'msa.file'});
		%len = map { length $_ => 1 } values %{ $data };
	}
	my $track = $args->{track};
	croak "$sub: tracker \"$track\" isn't in the alignment" if defined $track && !defined $data->{$track};
	my ($aln_len) = keys %len;
	my @undef = grep { !defined $data->{$_} } sort keys %{ $color };
	croak "$sub: \"color.residues\" names proteins that aren't in the alignment: @undef" if @undef;
	my @proteins = @{ $args->{order} // [sort { lc $a cmp lc $b } keys %{ $data }] };
	my @bad = grep { !defined $data->{$_} } @proteins;
	croak "$sub: \"order\" names proteins that aren't in the alignment: @bad" if @bad;
	my $width = $args->{'row.width'} // 100;
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
			(my $name = $protein) =~ s/([#&^_%])/\\$1/g;
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
	my $size = $args->{'table.text.size'} // '\footnotesize';
	my $caption = $args->{caption} // '';
	open my $tex, '>', $args->{'output.tex.file'} or croak "$sub: can't write $args->{'output.tex.file'}: $!";
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
	close $tex or croak "$sub: can't write $args->{'output.tex.file'}: $!";
	# LikeR's write_table confirmation line: black (30) on cyan (46), then
	# reset (0), written out rather than loaded from Term::ANSIColor.
	print STDOUT "wrote \e[30;46m$args->{'output.tex.file'}\e[0m\n";
	$args->{'output.tex.file'};
}

1;
__END__

=head1 SYNOPSIS

 use Bioinf::Basic ':all';

 my $seqs = fasta2hash('proteome.fa.gz');      # { defline => sequence }
 my $one  = fasta2hash('proteome.fa', 'P12345'); # just that sequence
 hash2fasta_file($seqs, 'copy.fa');

 my $hits = get_best_alignment_hit('blast.json', 'bit_score');

 # one clustalo run: plot_msa keeps the guide tree, and plot_phylo draws it
 plot_msa(
 	fasta       => 'orthologs.fa',     # or { name => sequence }
 	filename    => 'msa.svg',          # .png, .pdf, ... too
 	'tree.file' => 'orthologs.newick',
 	title       => 'EF-3',
 );
 plot_phylo('tree.file' => 'orthologs.newick', filename => 'tree.svg', title => 'EF-3');

=head1 DESCRIPTION

Nothing is exported by default; ask for functions by name or with C<:all>.
Every function dies (via L<Carp/croak>) on bad arguments, with a message that
starts with the name of the function that raised it, as in
C<fasta2hash: couldn't find "four" in x.fa>. C<plot_msa>,
C<plot_phylo>, C<msa_quality_table> and C<clustal_view_residues> take
C<< name => value >> pairs, not a hash ref.

Every PNG, SVG, PDF, PS or EPS image these functions draw carries its
provenance as C<Creator> metadata: the calling script (as the working
directory plus the script's name, like Matplotlib::Simple), the function,
this file and its version, the computer it ran on (hostname, operating system
and perl version), and what drew it: C<msa_plot.py> and the matplotlib
version, and for a tree the Biopython version too. For example:

  /home/me/work/run.pl called using "plot_msa" in /.../Bioinf/Basic.pm
  version 0.01 on host myhost (linux, perl 5.44.0), drawn by /.../msa_plot.py
  with matplotlib 3.10.0

An SVG holds it in C<< <dc:creator> >>; C<exiftool> or C<identify -verbose>
shows it in a PNG or PDF.

Clustal Omega, BLAST+ and the Python that draws the plots are the ones
L<Alien::Bioinf> installed alongside this module. See there to check them for
updates or to update them.

=head1 FUNCTIONS

=head2 fasta2hash($file, $key)

Reads a FASTA file (gzip-compressed if its name ends in C<.gz>). Returns a hash
ref of defline (without the C<< > >>) => sequence, or, with C<$key>, just the
sequence of that defline, reading no further than the record after it. A
defline that appears twice is warned about and its sequences concatenated.
Line endings may be C<\n> or C<\r\n>.

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
order.

=head2 plot_msa(%args)

Aligns sequences with Clustal Omega and draws the alignment. Returns a hash ref
of the files made: C<filename>, C<msa.file>, and C<tree.file> if it was given.
Once the image is written it prints C<wrote> and its file name to STDOUT, in
black on yellow.

=over

=item fasta, filename (required)

The sequences, as a FASTA file name or a hash ref of name => sequence (the
function tells the two apart by whether it is a reference); and the image to
draw, whose extension picks the format.

=item msa.file, tree.file

Where to keep clustalo's alignment (FASTA; a temporary file otherwise) and its
guide tree (newick; not made otherwise). Keep the tree to draw it with
C<plot_phylo> without aligning a second time.

=item order

Names, first to last (default: the input order; for a hash, sorted). Only
these are drawn, and the first is drawn at the bottom.

=item labels

A hash ref of name => label to show instead. matplotlib mathtext works:
C<< 'C.albicans' => '$\it{C. albicans}$' >>.

=item active.site.aa, query

C<< { His395 => 395, ... } >>: a dashed vertical line at each of these
1-based residue numbers of the sequence named by C<query>.

=item title, xlabel, ylabel, threads, clustal.args

Plot title; axis labels (default "Amino Acid Residue" and "Protein &
Species"); clustalo threads (default 1); an array ref of extra clustalo
arguments.

=back

With fewer than two sequences it warns and returns an empty hash ref.

=head2 plot_phylo(%args)

Draws a guide tree, with Biopython's C<Bio.Phylo>. Returns a hash ref of the
files made or used: C<output.file>, and C<tree.file> and C<msa.file> as below.

 plot_phylo(fasta => 'orthologs.fa');   # writes phylo.svg

=over

=item output.file

The image to draw (default F<phylo.svg>, in the working directory); the
extension picks the format.

=item fasta, tree.file

With C<fasta> (as for C<plot_msa>), the sequences are aligned with Clustal
Omega and the guide tree drawn; C<tree.file> and C<msa.file> then say where to
keep the tree and alignment, and C<threads> and C<clustal.args> are as for
C<plot_msa>. Without C<fasta>, C<tree.file> is an existing newick file to
draw, such as one C<plot_msa> kept, and no alignment is made.

=item labels, title

A hash ref of name => label for the tips, as for C<plot_msa>; plot title.

=back

With C<fasta> of fewer than two sequences it warns and returns an empty hash
ref.

Besides the C<Creator> every image carries, a tree records where it came from
in the image's own metadata, so the file can be traced and checked without
the script that made it. An SVG holds all of this as Dublin Core in its
C<< <metadata> >>, and a PNG as text chunks:

=over

=item Title, Date

The title (default "Phylogenetic tree"), and when it was drawn, to the second
and with the UTC offset. Where C<SOURCE_DATE_EPOCH> is set, matplotlib's date
from it is kept instead, so a reproducible build stays reproducible.

=item Source

The input, with the full path and SHA-256 of each file: the FASTA file (or the
number of sequences in a hash ref), the ungapped copy of it that Clustal Omega
read, and the alignment Clustal Omega wrote; or the newick file drawn.

=item Description

How the tree was made: the Clustal Omega version and the exact command line;
how many negative branch lengths were drawn as 0; the tips, with the label each
was shown as; and the whole newick tree that was drawn.

=item Identifier, Relation

The SHA-256 of the newick file drawn, and the C<msa.file> and C<tree.file>
kept with the image, if any.

=item Contributor

Each piece of software that took part, with its version and, where it runs as
a program, its path: Clustal Omega, Bioinf::Basic, Alien::Bioinf, perl,
Python, matplotlib, Biopython and NumPy.

=item Keywords

"phylogenetic tree", "newick" and the tip names.

=back

A PDF keeps the Title, the Description (as its Subject) and the Keywords; PS
and EPS keep only the C<Creator>, and other formats nothing.

=head2 msa_quality_table(%args)

Draws an all-against-all BLAST score table with matplotlib, and returns
C<filename>. The simplest call is

 msa_quality_table(fasta => 'orthologs.fa', filename => 'scores.png');

=over

=item fasta, filename (required)

The sequences, as a FASTA file name or a hash ref of name => sequence, as for
C<plot_msa>; they are aligned all against all with C<blastp>. Gaps are stripped
first, so an aligned FASTA, such as the one C<plot_msa> keeps, will do. Names
must look like C<Genus.species[.strain]>. C<fasta> is not needed when an
existing C<alignment.json> is given. C<unaligned.fa> is its old name.

=item alignment.json

The C<blastp -outfmt 15> report, as a parsed hash ref or a file name. An
existing file is read, and nothing is aligned; otherwise C<blastp> is run on
C<fasta> and its report kept there for next time. Without it, the report is a
temporary file.

=item metric, order, normalize, logscale.add, default_undefined, title, cb_label, cb_min, cb_max, cblogscale, show.numbers

C<metric> is an hsp field (default C<score>). C<msa.file> is accepted and
ignored, for old callers.

=back

=head2 clustal_view_residues(%args)

Writes an alignment as LaTeX tables with chosen residues coloured, and returns
C<output.tex.file>, printing C<wrote> and that file name (on cyan) to STDOUT
once it is written. C<msa.file> (required) is a FASTA file. One whose
sequences are all one length, such as the one C<plot_msa> keeps, is shown as it
is; one whose sequences differ in length is first aligned with Clustal Omega,
into a temporary file (C<msa.file> is never written), using C<threads> and
C<clustal.args> as C<plot_msa> does.
C<color.residues> is C<< { protein => { residue number => colour } } >>, where
residue numbers are 1-based and a colour is an xcolor name or C<[r, g, b]>; a
coloured column is coloured in every protein. C<track> adds a row under that
protein with its coloured residue numbers. Also C<order>, C<row.width>
(default 100 columns), C<split> (blocks per table, default 4), C<caption>,
C<label> and C<table.text.size> (default C<\footnotesize>).

=cut
