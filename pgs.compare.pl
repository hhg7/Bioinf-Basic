#!/usr/bin/env perl
# Polygenic scores for the samples of one or more VCFs, each placed in a
# reference panel scored the same way: its z-score and percentile within each
# population of the panel. A raw score means nothing on its own -- PGS000138's
# weights allow anything from -776 to +774 -- only its place among people
# scored on the same variants does.
#
#   perl pgs.compare.pl --weights PGS000138.txt.gz \
#   	--vcf /home/con/Documents/Genomics/genome.vcf.gz --fill
#   perl pgs.compare.pl --weights PGS000138.txt.gz --weights PGS000907.txt.gz \
#   	--vcf a.vcf.gz --vcf cohort.bcf --group pop --out depression.tsv
#   perl pgs.compare.pl --weights PGS000138.txt.gz --vcf cohort.vcf.gz \
#   	--reference cohort.bcf --panel cohort.groups.tsv   # any other panel
#   perl pgs.compare.pl --weights PGS000138.txt.gz --vcf a.vcf.gz --no-reference
#
# The default reference is the 2,504 people of 1000 Genomes phase 3 (GRCh37),
# and the default panel is that release's sample table,
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel,
# which is downloaded to the default path if it is not there. A panel is any
# tab-separated table with a header naming a "sample" column and the --group
# column; with --panel '' every reference sample is in one group, ALL.
#
# A single-sample VCF from a sequencing service usually lists only the sites
# where the sample carries an ALT allele, so a variant it doesn't list is
# missing rather than counted as zero copies, and the score is not comparable
# with the panel's. --fill writes those sites in as homozygous for the
# reference's REF allele, which is right only where the sample was sequenced.
# On 2026-10-17, genome.vcf.gz scored 1.674 on 12,015 of
# PGS000138's 22,274 variants unfilled, and 1.037 on 22,270 filled, the
# same 22,270 that each 1000 Genomes sample was scored on.
#
# The panel's scores are kept in --cache, one file per scoring file, reference
# and match, and reused; --refresh recomputes them. Scoring 1000 Genomes on
# PGS000138 took 95 s and 1.7 GB in 500-sample groups on 2026-10-17; the
# reference is read one chromosome and one --batch of samples at a time, so
# that a scoring file of a million variants still fits in memory.
#
# The Cn::Bioinf in this directory's blib/ is used if it has been built
# (perl Makefile.PL && make), and the installed one otherwise.
use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use FindBin;
use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use File::Basename 'basename';
use File::Path 'make_path';
use File::Spec;
use File::Temp 'tempdir';
use Getopt::Long;
use IO::Uncompress::Gunzip;
use Stats::LikeR qw(mean median quantile read_table sd);
use Cn::Bioinf 'pgs_score';

my $panel_url = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel';
my %r = (
	batch     => 250,       # reference samples read at a time
	bcftools  => 'bcftools',
	cache     => 'pgs_cache',
	fill      => 0,         # 1 = sites a target VCF doesn't list are 0/0, with the reference's REF and ALT
	group     => 'super_pop', # the panel column whose groups the scores are compared within
	match     => 'position',  # passed to pgs_score: 'position' or 'id'
	out       => 'pgs_scores.tsv',
	panel     => '/home/con/Documents/Genomics/integrated_call_samples_v3.20130502.ALL.panel',
	reference => '/home/con/Documents/Genomics/1000genomes.bcf',
	refresh   => 0,         # 1 = recompute the reference's scores rather than read --cache
	threads   => 4,         # bcftools' compression threads
);
my $default_panel = $r{panel};
GetOptions(\%r, 'weights=s@', 'vcf=s@', 'reference=s', 'panel=s', 'group=s', 'match=s', 'fill!', 'out=s', 'cache=s',
	'refresh!', 'batch=i', 'threads=i', 'bcftools=s', 'no-reference', 'help')
	or die "pgs.compare.pl: bad options; see the comment at the top of $0\n";
if ($r{help} || !$r{weights}) {
	open my $self, '<', $0;
	while (<$self>) { last if /^use /; print s/^# ?//r if /^#[^!]/ }
	exit !$r{help};
}
die "pgs.compare.pl: --batch must be 1 or more, not $r{batch}\n" if $r{batch} < 1;
die "pgs.compare.pl: --fill needs a --reference to take REF and ALT from\n" if $r{fill} && $r{'no-reference'};
$r{vcf} //= [];
my $tmp = tempdir(CLEANUP => 1);

# ---- running bcftools

sub run (@cmd) {
	system(@cmd) == 0 and return;
	my $why = $? == -1 ? "it didn't start: $!" : $? & 127 ? 'killed by signal ' . ($? & 127) : 'exit ' . ($? >> 8);
	die "run: \"@cmd\" failed: $why\n";
}

sub lines (@cmd) {
	open my $fh, '-|', @cmd;
	chomp(my @lines = <$fh>);
	close $fh;
	@lines;
}

# 'chr' if $file names its chromosomes as UCSC does (chr1), '' if as Ensembl
# and the PGS Catalog do (1), judged by its first ##contig line, or else by its
# first record.
sub chr_prefix ($file) {
	foreach (lines($r{bcftools}, 'view', '-h', $file)) {
		return /\A##contig=<ID=chr/i ? 'chr' : '' if /\A##contig=<ID=/;
	}
	open my $fh, '-|', $r{bcftools}, 'view', '-H', $file;
	my $first = <$fh> // '';
	# closed early on purpose, so bcftools dies of SIGPIPE, and close says so
	eval { close $fh };
	$first =~ /\Achr/i ? 'chr' : '';
}

# ---- the scoring file

# The '#key=value' lines at the top of a PGS Catalog scoring file, which
# read_table skips as comments: pgs_id, trait_reported, genome_build and so on.
sub pgs_header ($file) {
	my $fh = IO::Uncompress::Gunzip->new($file, Transparent => 1) or die "pgs_header: can't read $file: $IO::Uncompress::Gunzip::GunzipError\n";
	my %h;
	while (my $line = <$fh>) {
		last unless $line =~ /\A#/;
		$h{$1} = $2 if $line =~ /\A#([^#=]+)=(.*?)\s*\z/;
	}
	close $fh;
	\%h;
}

# The columns pgs_score finds a variant by, which it prefers in the same order.
sub position_cols ($weights) {
	exists $weights->[0]{hm_pos} ? ('hm_chr', 'hm_pos') : ('chr_name', 'chr_position');
}

# A chromosome as pgs_score compares them, so the scoring file's chromosomes
# and a VCF's can be grouped together.
sub norm_chrom ($c) {
	$c = uc $c;
	$c =~ s/\ACHR//;
	{ 23 => 'X', 24 => 'Y', 26 => 'MT', M => 'MT' }->{$c} // $c;
}

# A regions file of every position of $weights, named with $prefix as the VCF
# it is for names them. A variant with no position, such as one a harmonized
# file could not place, is left out, and pgs_score counts it as missing.
sub write_regions ($weights, $prefix, $name) {
	my ($cc, $pc) = position_cols($weights);
	die "write_regions: the scoring file has no $cc and $pc columns; use the PGS Catalog's harmonized file (..._hmPOS_GRCh37.txt.gz)\n"
		unless exists $weights->[0]{$cc} && exists $weights->[0]{$pc};
	my %seen;
	my @at;
	foreach my $w (@{ $weights }) {
		my ($c, $p) = @{ $w }{$cc, $pc};
		next unless defined $c && length $c && defined $p && $p =~ /\A[1-9][0-9]*\z/;
		$c = norm_chrom($c);
		push @at, [$c, $p] unless $seen{"$c:$p"}++;
	}
	my $file = File::Spec->catfile($tmp, $name);
	open my $fh, '>', $file;
	foreach my $at (sort { ($a->[0] =~ /\A[0-9]+\z/ ? $a->[0] : 1e9) <=> ($b->[0] =~ /\A[0-9]+\z/ ? $b->[0] : 1e9) || $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @at) {
		my $c = $prefix && $at->[0] eq 'MT' ? 'M' : $at->[0];
		say {$fh} "$prefix$c\t$at->[1]";
	}
	close $fh;
	$file;
}

# ---- the reference panel

# The reference's records at the positions of $weights, all samples, as an
# indexed BCF made once per scoring file.
my %sites;
sub reference_sites ($id, $weights) {
	return $sites{$id} if $sites{$id};
	my $regions = write_regions($weights, chr_prefix($r{reference}), "$id.ref.regions");
	my $bcf = File::Spec->catfile($tmp, "$id.ref.bcf");
	say "extracting $id's positions from $r{reference}";
	# --regions-overlap 0 keeps the records that start at a position; the
	# default, 1, would also keep a deletion that starts before it
	run($r{bcftools}, 'view', '-R', $regions, '--regions-overlap', 0, '--threads', $r{threads}, '-Ob', '-o', $bcf, $r{reference});
	run($r{bcftools}, 'index', $bcf);
	$sites{$id} = $bcf;
}

# sample => its group in the panel's --group column
sub read_panel () {
	return {} if $r{panel} eq '';
	if (!-e $r{panel} && $r{panel} eq $default_panel) {
		require HTTP::Tiny;
		say "downloading $panel_url to $r{panel}";
		my $res = HTTP::Tiny->new->mirror($panel_url, $r{panel});
		die "read_panel: couldn't download $panel_url: $res->{status} $res->{reason}\n" unless $res->{success};
	}
	open my $fh, '<', $r{panel};
	chomp(my $head = <$fh>);
	my @cols = split /\t/, $head;
	my ($s) = grep { $cols[$_] eq 'sample' } 0 .. $#cols;
	my ($g) = grep { $cols[$_] eq $r{group} } 0 .. $#cols;
	die "read_panel: $r{panel} has no \"sample\" column\n" unless defined $s;
	die "read_panel: $r{panel} has no \"$r{group}\" column; it has " . join(', ', grep { length } @cols) . "\n" unless defined $g;
	my %group;
	while (<$fh>) {
		chomp;
		my @f = split /\t/;
		$group{ $f[$s] } = $f[$g] if defined $f[$g] && length $f[$g];
	}
	close $fh;
	\%group;
}

my @fields = qw(score variants alleles missing mismatched);

# sample => { group, and pgs_score's fields } for every sample of the
# reference, from --cache or else scored a chromosome and a batch at a time and
# summed, which a score is, since it is a sum over variants.
sub reference_scores ($id, $weights, $panel) {
	(my $ref = basename($r{reference})) =~ s/\.(?:bcf|vcf\.gz|vcf\.bgz|vcf)\z//;
	my $cache = File::Spec->catfile($r{cache}, "$id.$ref.$r{match}.tsv");
	my %score;
	if (-e $cache && !$r{refresh}) {
		my $t = read_table($cache, sep => "\t", 'output_type' => 'hoh', 'row_names' => 'sample');
		$score{$_} = { %{ $t->{$_} }, group => $panel->{$_} // 'ALL' } foreach keys %{ $t };
		say "read the reference's scores from $cache";
		return \%score;
	}
	my $bcf = reference_sites($id, $weights);
	my @samples = lines($r{bcftools}, 'query', '-l', $bcf);
	my @chroms = map { (split /\t/)[0] } lines($r{bcftools}, 'index', '-s', $bcf);
	my ($cc) = position_cols($weights);
	my %by_chrom; # normalised chromosome => the rows of $weights on it
	foreach my $w (@{ $weights }) {
		push @{ $by_chrom{ norm_chrom($w->{$cc}) } }, $w if defined $w->{$cc} && length $w->{$cc};
	}
	$score{$_} = { map { $_ => 0 } @fields } foreach @samples;
	my $batches = int((@samples + $r{batch} - 1) / $r{batch});
	foreach my $b (0 .. $batches - 1) {
		my @batch = @samples[$b * $r{batch} .. ($b == $batches - 1 ? $#samples : ($b + 1) * $r{batch} - 1)];
		my $list = File::Spec->catfile($tmp, 'batch.txt');
		open my $fh, '>', $list;
		say {$fh} $_ foreach @batch;
		close $fh;
		printf "scoring %s in the reference: samples %d to %d of %d\n", $id, $b * $r{batch} + 1, $b * $r{batch} + @batch, scalar @samples;
		foreach my $chrom (@chroms) {
			my $w = $by_chrom{ norm_chrom($chrom) } or next;
			my $chunk = File::Spec->catfile($tmp, 'chunk.vcf.gz');
			run($r{bcftools}, 'view', '-r', $chrom, '-S', $list, '--threads', $r{threads}, '-Oz', '-o', $chunk, $bcf);
			my $frame = read_table($chunk);
			next unless %{ $frame };
			my $s = do {
				local $SIG{__WARN__} = sub { warn "[$id, reference chromosome $chrom] $_[0]" };
				pgs_score(vcf => $frame, weights => $w, match => $r{match});
			};
			foreach my $sample (keys %{ $s }) {
				$score{$sample}{$_} += $s->{$sample}{$_} foreach qw(score variants alleles mismatched);
			}
		}
	}
	# Variants on a chromosome the reference doesn't have were never given to
	# pgs_score, so missing is counted from the total, as pgs_score counts it.
	$_->{missing} = @{ $weights } - $_->{variants} foreach values %score;
	make_path($r{cache});
	open my $out, '>', $cache;
	say {$out} join "\t", 'sample', @fields;
	say {$out} join "\t", $_, @{ $score{$_} }{@fields} foreach sort keys %score;
	close $out;
	say "wrote $cache";
	$score{$_}{group} = $panel->{$_} // 'ALL' foreach keys %score;
	\%score;
}

# ---- the VCFs to score

# pgs_score's result for every sample of $vcf, read at the positions of
# $weights only, and, with --fill, given a 0/0 record at every position of the
# reference that $vcf has no record at.
sub target_scores ($vcf, $id, $weights) {
	my $regions = write_regions($weights, chr_prefix($vcf), "$id.target.regions");
	my $part = File::Spec->catfile($tmp, 'target.vcf.gz');
	# -T streams, and so needs no index for $vcf
	run($r{bcftools}, 'view', '-T', $regions, '--threads', $r{threads}, '-Oz', '-o', $part, $vcf);
	my @samples = lines($r{bcftools}, 'query', '-l', $vcf);
	my $frame = read_table($part);
	my $homref = grep { my $rec = $_; grep { ($rec->{"$_.GT"} // '') =~ m{\A0[/|]0\z} } @samples } values %{ $frame };
	if ($r{fill}) {
		my $prefix = chr_prefix($vcf);
		my %have = map { (norm_chrom($_->{CHROM}) . ":$_->{POS}" => 1) } values %{ $frame };
		my $added = 0;
		foreach (lines($r{bcftools}, 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\n', reference_sites($id, $weights))) {
			my ($c, $p, $ref, $alt) = split /\t/;
			my $n = norm_chrom($c);
			next if $have{"$n:$p"};
			$frame->{"fill:$n:$p:$ref:$alt"} = { CHROM => $prefix . ($prefix && $n eq 'MT' ? 'M' : $n), POS => $p, ID => '.',
				REF => $ref, ALT => $alt, map { ("$_.GT" => '0/0') } @samples };
			$added++;
		}
		say "filled in $added positions of $id that " . basename($vcf) . ' does not list, as 0/0';
	} elsif (%{ $frame } && !$homref) {
		say basename($vcf) . " has no 0/0 genotype at $id's positions, so it probably lists only the sites where a sample"
			. ' carries an ALT, and its scores are not comparable with the reference\'s; --fill writes the rest in as 0/0';
	}
	unless (%{ $frame }) {
		say basename($vcf) . " has no record at any of $id's positions; check that it is on the same genome build";
		return { map { $_ => { score => 0, variants => 0, alleles => 0, missing => scalar @{ $weights }, mismatched => 0 } } @samples };
	}
	local $SIG{__WARN__} = sub { warn '[' . basename($vcf) . "] $_[0]" };
	pgs_score(vcf => $frame, weights => $weights, match => $r{match});
}

# ---- comparison

# n, mean, sd, 5% and 95% quantiles (R's type 7) and sorted scores of each
# group of the reference, and of ALL. sd is undef for a group of one, where
# Stats::LikeR's sd dies, as R's gives NA.
sub group_stats ($ref) {
	my %by;
	foreach my $s (values %{ $ref }) {
		push @{ $by{ $s->{group} } }, $s->{score};
		push @{ $by{ALL} }, $s->{score} unless $s->{group} eq 'ALL';
	}
	my %stats;
	foreach my $g (keys %by) {
		my @v = sort { $a <=> $b } @{ $by{$g} };
		my $q = quantile(x => \@v, probs => [0.05, 0.95]);
		$stats{$g} = { n => scalar @v, mean => mean(\@v), sd => @v > 1 ? sd(\@v) : undef, q05 => $q->{'5%'}, q95 => $q->{'95%'}, sorted => \@v };
	}
	\%stats;
}

# The percentage of @$sorted below $x, counting a tie as half below.
# Stats::LikeR has no ECDF to evaluate at a value outside the sample.
sub percentile ($x, $sorted) {
	my ($below, $tied) = (0, 0);
	foreach (@{ $sorted }) {
		last if $_ > $x;
		$_ < $x ? $below++ : $tied++;
	}
	100 * ($below + $tied / 2) / @{ $sorted };
}

# ---- main

my $panel = $r{'no-reference'} ? {} : read_panel();
open my $out, '>', $r{out};
say {$out} join "\t", qw(pgs_id vcf sample score variants missing mismatched comparable group ref_n ref_mean ref_sd z percentile);
foreach my $file (@{ $r{weights} }) {
	my $h = pgs_header($file);
	my $id = $h->{pgs_id} // basename($file) =~ s/\..*//r;
	my $weights = read_table($file, sep => "\t");
	say "\n$id", (defined $h->{pgs_name} ? " ($h->{pgs_name})" : ''), (defined $h->{trait_reported} ? ": $h->{trait_reported}" : ''),
		'; genome build ', $h->{genome_build} // 'not given', (defined $h->{HmPOS_build} ? ", harmonized to $h->{HmPOS_build}" : ''),
		'; ', scalar @{ $weights }, ' variants';
	my ($ref, $stats, $ref_median) = ({}, {}, 0);
	unless ($r{'no-reference'}) {
		$ref = reference_scores($id, $weights, $panel);
		$stats = group_stats($ref);
		$ref_median = median([map { $_->{variants} } values %{ $ref }]);
		say "reference: ", basename($r{reference}), ", $ref_median of the variants scored in its median sample";
		printf "  %-6s %5s %8s %6s %8s %8s\n", qw(group n mean sd 5% 95%);
		foreach my $g (sort keys %{ $stats }) {
			my $s = $stats->{$g};
			printf "  %-6s %5d %8.3f %6s %8.3f %8.3f\n", $g, $s->{n}, $s->{mean}, defined $s->{sd} ? sprintf('%.3f', $s->{sd}) : 'NA',
				@{ $s }{qw(q05 q95)};
		}
	}
	foreach my $vcf (@{ $r{vcf} }) {
		my $t = target_scores($vcf, $id, $weights);
		foreach my $sample (sort keys %{ $t }) {
			my $s = $t->{$sample};
			# On 2026-10-17 the filled genome was scored on all 22,270 of the
			# variants a 1000 Genomes sample was, the unfilled one on 54% and an
			# exome on 3%, so anything short of 99% is a different set of variants.
			my $comparable = $r{'no-reference'} ? '' : $s->{variants} >= 0.99 * $ref_median ? 1 : 0; # '' = no reference to compare with
			printf "%s %s: score %.4f on %d variants (%d missing, %d with alleles that disagree)%s\n", basename($vcf), $sample,
				@{ $s }{qw(score variants missing mismatched)},
				$r{'no-reference'} ? '' : $comparable ? '' : " -- NOT comparable: the reference's median sample was scored on $ref_median";
			my @row = ($id, $vcf, $sample, @{ $s }{qw(score variants missing mismatched)}, $comparable);
			unless (%{ $stats }) {
				say {$out} join "\t", @row, ('') x 6;
				next;
			}
			foreach my $g (sort keys %{ $stats }) {
				my $st = $stats->{$g};
				# '' = no z, for a group whose scores have no spread
				my $z = $st->{sd} ? ($s->{score} - $st->{mean}) / $st->{sd} : '';
				my $pct = percentile($s->{score}, $st->{sorted});
				printf "  %-6s z %6s, above %5.1f%% of %d\n", $g, $z eq '' ? 'NA' : sprintf('%.2f', $z), $pct, $st->{n};
				say {$out} join "\t", @row, $g, $st->{n}, $st->{mean}, $st->{sd} // '', $z, $pct;
			}
		}
	}
}
close $out;
say "\nwrote $r{out}";
