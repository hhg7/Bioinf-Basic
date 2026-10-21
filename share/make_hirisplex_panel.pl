#!/usr/bin/env perl
# Writes share/hirisplex_s_panel.tsv, the HIrisPlex-S marker panel that
# hirisplex_genotypes() reads. Run it from the distribution's root:
#
#	perl share/make_hirisplex_panel.pl > share/hirisplex_s_panel.tsv
#
# Cn::Bioinf never runs this, and [PruneFiles] keeps it out of the tarball:
# the committed .tsv is what ships and what is
# used. Re-run it only to refresh the panel against its sources, and read the
# diff before committing, since a changed coordinate or allele silently changes
# every genotype the module extracts.
#
# Two sources, both fetched over the network:
#
# - https://hirisplex.erasmusmc.nl/ -- the webtool's own input form, which is
#   the authority for which markers are in the panel, in which order, under
#   which gene, and which allele each model counts. The page lays out three
#   tables, HIrisPlex-S (41 markers), then HIrisPlex (24), then IrisPlex (6),
#   so the membership flags come from the page too rather than from a list
#   typed out here.
# - Ensembl's REST API, grch37.rest.ensembl.org and rest.ensembl.org, for each
#   marker's position and its alleles on the forward strand of GRCh37 and
#   GRCh38.
#
# The webtool asks for the count of an allele on the assay's strand, which is
# not always the forward strand a VCF is written on. The HIrisPlex-S user
# manual (https://hirisplex.erasmusmc.nl/pdf/hirisplex.erasmusmc.nl.pdf,
# version 2.0, 2018, read 2026-10-20) puts it this way: "The only SNP that may
# cause confusion and therefore must be converted (from NCBI's forward
# orientation) is rs6119471. All other G/C, A/T SNPs are in the correct
# orientation for input or are opposite alleles i.e. SNP G/A, input requests
# C/T." So count_allele below is worked out as:
#
# - rs6119471, the documented exception, is complemented. Its forward alleles
#   are C/A/G, so its input allele C is also a forward allele, and matching the
#   input allele against the forward alleles would silently count the wrong one.
#   This is the single marker where that happens, and the reason this script
#   does not simply match.
# - otherwise the input allele when the forward alleles include it, and its
#   complement when they do not.
#
# Every result is then checked against the forward alleles, so that a refreshed
# Ensembl record that no longer carries the counted allele stops the run rather
# than writing a panel that counts nothing.
use 5.044;
use warnings FATAL => 'all';
use JSON::MaybeXS qw(decode_json encode_json);
use HTTP::Tiny;

my $Form   = 'https://hirisplex.erasmusmc.nl/';
my %Rest   = (b38 => 'https://rest.ensembl.org', b37 => 'https://grch37.rest.ensembl.org');
my $Expect = 71; # the page's three tables: 41 + 24 + 6

# rs6119471 is complemented whatever its forward alleles say; see above.
my %Flip = (rs6119471 => 1);

my %Complement = (A => 'T', C => 'G', G => 'C', T => 'A');

my $http = HTTP::Tiny->new(agent => 'Cn-Bioinf make_hirisplex_panel.pl ', timeout => 120);

sub get {
	my ($url) = @_;
	my $r = $http->get($url);
	die "GET $url failed: $r->{status} $r->{reason}\n" unless $r->{success};
	$r->{content};
}

# Each table row of the form is "<n> <gene> <rsID> <allele>".
my $html = get($Form);
my @row;
while ($html =~ m{<td[^>]*>\s*([0-9]+)\s*</td>\s*<td>\s*<b><em>([^<]+)</em></b>\s*</td>\s*<td>\s*(rs[0-9]+)\s*</td>\s*<td>\s*([ACGT])\s*</td>}gs) {
	push @row, { n => $1, gene => $2, rsid => $3, input_allele => $4 };
}
die "$Form gave " . scalar(@row) . " marker rows, not $Expect; its layout has changed\n" unless @row == $Expect;
my @panel = @row[0 .. 40];
my %hirisplex = map { $_->{rsid} => 1 } @row[41 .. 64];
my %irisplex  = map { $_->{rsid} => 1 } @row[65 .. 70];
die "the 24-marker table is not the first 24 of the 41\n"
	if grep { !$hirisplex{ $panel[$_]{rsid} } } 0 .. 23;
my %in_panel = map { $_->{rsid} => 1 } @panel;
die "the 6-marker table is not a subset of the 41\n" if grep { !$in_panel{$_} } keys %irisplex;

# Ensembl, both builds at once. A marker whose rs number has been merged comes
# back under its current name with the asked-for one among its synonyms, so the
# reply is keyed back to the panel by synonym as well as by name.
my %var;
foreach my $build (sort keys %Rest) {
	my $body = encode_json({ ids => [map { $_->{rsid} } @panel] });
	my $r = $http->request('POST', "$Rest{$build}/variation/homo_sapiens?content-type=application/json",
		{ content => $body, headers => { 'Content-Type' => 'application/json', Accept => 'application/json' } });
	die "POST $Rest{$build} failed: $r->{status} $r->{reason}\n" unless $r->{success};
	my $got = decode_json($r->{content});
	foreach my $name (keys %{ $got }) {
		my $v = $got->{$name};
		my ($m) = grep { $_->{coord_system} eq 'chromosome' } @{ $v->{mappings} || [] };
		next unless $m;
		foreach my $key ($name, @{ $v->{synonyms} || [] }) {
			$var{$key}{$build} = { %{ $m }, name => $name, var_class => $v->{var_class} };
		}
	}
}

print "#n\tgene\trsid\tinput_allele\tcount_allele\tkind\tb37_chrom\tb37_pos\tb38_chrom\tb38_pos\tb38_alleles\tcurrent_rsid\tirisplex\thirisplex\n";
foreach my $p (@panel) {
	my ($rsid, $in) = @{ $p }{qw(rsid input_allele)};
	my ($m37, $m38) = ($var{$rsid}{b37}, $var{$rsid}{b38});
	die "$rsid is not in Ensembl on either build\n" unless $m37 || $m38;
	my $alleles = $m38 ? $m38->{allele_string} : $m37->{allele_string};
	my $class   = $m38 ? $m38->{var_class}     : $m37->{var_class};
	my @forward = split m{/}, $alleles;

	# "allele" = count count_allele among the called alleles; "insertion" =
	# count the longer of REF and ALT, for the one marker that is an indel.
	my ($kind, $count) = ("allele", $in);
	if ($class ne 'SNP') {
		$kind = "insertion";
		$count = $in;
	} elsif ($Flip{$rsid} || !grep { $_ eq $in } @forward) {
		$count = $Complement{$in};
		die "$rsid: the counted allele $count is not among its forward alleles $alleles\n"
			unless grep { $_ eq $count } @forward;
	}
	printf "%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n",
		$p->{n}, $p->{gene}, $rsid, $in, $count, $kind,
		($m37 ? ($m37->{seq_region_name}, $m37->{start}) : ('', '')),
		($m38 ? ($m38->{seq_region_name}, $m38->{start}) : ('', '')),
		$alleles, ($m38 ? $m38->{name} : $m37->{name}),
		($irisplex{$rsid} ? 1 : 0), ($hirisplex{$rsid} ? 1 : 0);
}
