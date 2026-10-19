#!/usr/bin/env perl
# A polygenic score for each sample of a VCF, compared with each population of
# 1000 Genomes phase 3 scored on the same weights:
#
#   perl pgs.score.pl [weights] [vcf]
#
# weights is a PGS Catalog scoring file (default PGS000138.txt.gz), which has
# to be downloaded first, from the link in the catalog's table of scores (such
# as depression.pgs_scores_data.csv), and vcf a .vcf or .vcf.gz on GRCh37
# (default the WGZ genome). The table, one row for
# each of the 26 populations and then each of the 5 super-populations ("all"),
# is printed and written to pgs.<score>.<vcf>.tsv: n, mean, sd, 5% and 95% of
# the population's scores, and each sample's z and percentile among them.
#
# 1000genomes.ancestry.tsv (sample, pop, super_pop, population) was made on
# 2026-10-17 from 1000 Genomes' phase 3 panel and population list,
#   https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
#   https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv
# because read_table refuses the panel itself, whose header has two more tabs
# than its rows.
#
# A VCF with no 0/0 genotype at the score's positions lists only the sites
# where its samples carry an ALT allele. The other positions are then written
# in as 0/0, with 1000 Genomes' REF and ALT, so that the samples are scored on
# the same variants as 1000 Genomes; that is right only where the sample was
# sequenced, as in a whole genome, not an exome. On 2026-10-17 the WGZ genome
# scored 1.037 on 22,270 of PGS000138's variants this way, as each 1000
# Genomes sample was, and 1.674 on 12,015 without.
#
# PLINK 2 does the scoring, through pgs_score, and bcftools cuts the score's
# positions out of the files, so that neither the genotypes nor the weights
# are ever held in perl. bcftools is the one on PATH.

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
#use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
#use Devel::Confess 'color';
use FindBin;
use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch"; # this directory's build of Cn::Bioinf, once made
use File::Basename 'basename';
use File::Temp 'tempdir';
use IO::Uncompress::Gunzip;
use Stats::LikeR;
use Cn::Bioinf 'pgs_score';
use Cn::Util;

my $weights_file  = shift // 'PGS000138.txt.gz';
my $vcf_file      = shift // '/home/con/Documents/Genomics/genome.vcf.gz';
my $kg_file       = '/home/con/Documents/Genomics/1000genomes.bcf';
my $ancestry_file = '/home/con/Documents/Genomics/1000genomes.ancestry.tsv';
my $tmp = tempdir(CLEANUP => 1);

sub run (@cmd) {
	system(@cmd) == 0 and return;
	die "\"@cmd\" failed: " . ($? == -1 ? "it didn't start: $!" : $? & 127 ? 'killed by signal ' . ($? & 127) : 'exit ' . ($? >> 8)) . "\n";
}

# the percentage of @v below $x, a tie counting as half below
sub percentile ($x, @v) {
	100 * (grep({ $_ < $x } @v) + grep({ $_ == $x } @v) / 2) / @v;
}

(my $id = basename($weights_file)) =~ s/\..*//;
die "pgs.score.pl: there is no $weights_file; the PGS Catalog has it at\n"
	. "  https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/$id/ScoringFiles/$id.txt.gz\n" unless -f $weights_file;

# The score's positions, read a line at a time from the scoring file, by the
# columns pgs_score reads: hm_chr and hm_pos in a harmonized file, chr_name and
# chr_position otherwise. 1000 Genomes names chromosomes 1, 2, ..., as the PGS
# Catalog does; the VCF's names (chr1 or 1) are written into a file of their own.
my $vcf_chr = (grep { /\A##contig=<ID=chr/ } map { chomp; $_ } do { open my $h, '-|', 'bcftools', 'view', '-h', $vcf_file; <$h> }) ? 'chr' : '';
my $n_weights = 0;
{
	my $in = IO::Uncompress::Gunzip->new($weights_file, Transparent => 1) or die "$weights_file: $IO::Uncompress::Gunzip::GunzipError\n";
	my ($head, @cols);
	open my $kg_regions, '>', "$tmp/kg.regions";
	open my $vcf_regions, '>', "$tmp/vcf.regions";
	while (my $line = <$in>) {
		next if $line =~ /\A#/;
		chomp $line;
		my @f = split /\t/, $line, -1;
		unless (@cols) {
			my %col = map { $f[$_] => $_ } 0 .. $#f;
			@cols = exists $col{hm_pos} ? @col{qw(hm_chr hm_pos)} : @col{qw(chr_name chr_position)};
			next;
		}
		$n_weights++;
		my ($chr, $pos) = @f[@cols];
		next unless defined $pos && length $pos;
		$chr =~ s/\Achr//i;
		say {$kg_regions} "$chr\t$pos";
		say {$vcf_regions} "$vcf_chr$chr\t$pos";
	}
	close $kg_regions;
	close $vcf_regions;
}

# 1000 Genomes' scores are kept in pgs.<score>.1000genomes.tsv and read back on
# the next run of the same score; delete that file to compute them again. The
# score's positions are cut out of 1000 Genomes by bcftools, which reads its
# index (--regions-overlap 0 keeps the records that start at a position, and
# not a deletion that spans one), and PLINK 2 scores all 2,504 people at once.
my %about; # pop => a row of the ancestry table with it, for its names
my %pop_of = map { $about{ $_->{pop} } = $_; ($_->{sample} => $_->{pop}) } @{ read_table($ancestry_file, sep => "\t") };
run('bcftools', 'view', '-R', "$tmp/kg.regions", '--regions-overlap', 0, '-Ob', '-o', "$tmp/kg.bcf", $kg_file);
my $kg_file_scores = "pgs.$id.1000genomes.tsv";
my $kg; # sample => { pop, super_pop, score, variants }
if (-f $kg_file_scores) {
	$kg = read_table($kg_file_scores, sep => "\t", 'output_type' => 'hoh', 'row_names' => 'sample');
	say "read 1000 Genomes' scores from $kg_file_scores";
} else {
	my $s = pgs_score(vcf => "$tmp/kg.bcf", weights => $weights_file);
	$kg = { map { $_ => { pop => $pop_of{$_}, super_pop => $about{ $pop_of{$_} }{super_pop}, %{ $s->{$_} }{qw(score variants)} } } keys %{ $s } };
	write_table($kg, $kg_file_scores, sep => "\t", 'row_names' => 'sample', 'col_names' => [qw(pop super_pop score variants)]);
}
my %scores; # pop or super_pop => the scores of its people
foreach my $person (values %{ $kg }) {
	push @{ $scores{ $person->{pop} } }, $person->{score};
	push @{ $scores{ $person->{super_pop} } }, $person->{score};
}
my $kg_variants = median([map { $_->{variants} } values %{ $kg }]); # how many variants a 1000 Genomes sample was scored on

# The VCF's records at the score's positions, streamed by bcftools (-T reads
# the file through, and so needs no index), and whether any genotype among
# them is 0/0. A VCF with none lists only the sites with an ALT allele, and is
# scored with 1000 Genomes' records at every other position written in as 0/0;
# any other VCF is scored as it is.
open my $part, '>', "$tmp/vcf.vcf";
{
	open my $h, '-|', 'bcftools', 'view', '-h', $vcf_file;
	print {$part} <$h>;
	close $h;
}
my (%have, $homref, @samples);
{
	open my $body, '-|', 'bcftools', 'view', '-H', '-T', "$tmp/vcf.regions", $vcf_file;
	while (my $line = <$body>) {
		print {$part} $line;
		my @f = split /\t/, $line;
		$have{"$f[0]:$f[1]"} = 1;
		$homref ||= grep { m{\A0[/|]0(?::|\s*\z)} } @f[9 .. $#f];
	}
	close $body;
}
open my $names, '-|', 'bcftools', 'query', '-l', $vcf_file;
chomp(@samples = <$names>);
close $names;
my $scored = $vcf_file;
unless ($homref) {
	my $added = 0;
	open my $sites, '-|', 'bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\n', "$tmp/kg.bcf";
	while (my $line = <$sites>) {
		chomp $line;
		my ($c, $p, $ref, $alt) = split /\t/, $line;
		next if $have{"$vcf_chr$c:$p"};
		say {$part} join "\t", "$vcf_chr$c", $p, '.', $ref, $alt, '.', '.', '.', 'GT', ('0/0') x @samples;
		$added++;
	}
	close $sites;
	close $part;
	run('bcftools', 'sort', '-T', $tmp, '-Oz', '-o', "$tmp/filled.vcf.gz", "$tmp/vcf.vcf");
	$scored = "$tmp/filled.vcf.gz";
	say basename($vcf_file) . " lists only sites with an ALT allele, so $added records of 1000 Genomes were written in as 0/0";
}
my $target = pgs_score(vcf => $scored, weights => $weights_file);
say "$_: score $target->{$_}{score} on $target->{$_}{variants} of $id\'s $n_weights variants (1000 Genomes: $kg_variants)"
	foreach @samples;

my @rows;
my %super = map { $about{$_}{super_pop} => 1 } keys %about;
foreach my $sp (sort keys %super) {
	foreach my $group ((sort grep { $about{$_}{super_pop} eq $sp } keys %about), $sp) {
		my @v = @{ $scores{$group} };
		my ($mean, $sd) = (mean(\@v), sd(\@v));
		my $q = quantile(x => \@v, probs => [0.05, 0.95]);
		my %row = (
			super_pop  => $sp,
			pop        => $group eq $sp ? 'all' : $group,
			population => $group eq $sp ? "all $sp" : $about{$group}{population},
			n          => scalar @v,
			mean       => sprintf('%.3f', $mean),
			sd         => sprintf('%.3f', $sd),
			'5%'       => sprintf('%.3f', $q->{'5%'}),
			'95%'      => sprintf('%.3f', $q->{'95%'}),
		);
		foreach my $s (@samples) {
			$row{"$s z"} = sprintf '%.2f', ($target->{$s}{score} - $mean) / $sd;
			$row{"$s percentile"} = sprintf '%.1f', percentile($target->{$s}{score}, @v);
		}
		push @rows, \%row;
	}
}
my @cols = (qw(super_pop pop population n mean sd 5% 95%), map { ("$_ z", "$_ percentile") } @samples);
(my $vcf_name = basename($vcf_file)) =~ s/\.vcf(?:\.gz)?\z//;
write_table(\@rows, "pgs.$id.$vcf_name.tsv", sep => "\t", col_names => \@cols);
view(\@rows, n => scalar @rows, cols => \@cols, max_width => 40);
