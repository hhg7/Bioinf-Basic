#!/usr/bin/env perl
# Writes share/blood_type_panel.tsv, the markers that blood_type() reads. Run
# it from the distribution's root:
#
#	perl share/make_blood_type_panel.pl > share/blood_type_panel.tsv
#
# Cn::Bioinf never runs this, and [PruneFiles] keeps it out of the tarball:
# the committed .tsv is what ships and what is used. It needs bcftools on PATH
# (1.21 was used) and the network, and took 10 minutes on 2026-11-10. Read the
# diff before committing, since a changed coordinate or allele silently
# changes every blood type the module gives.
#
# The panel has three parts.
#
# ABO. Three markers, each named by its change in the coding sequence, which
# is on the reverse strand, so the allele counted here is its complement:
#
#   rs8176719   c.261delG, the frameshift of the common O alleles (Yamamoto
#               et al., Nature 345:229-233, 1990). GRCh37 and GRCh38 both
#               carry the deletion, so a VCF writes the working gene as an
#               insertion of C, and that insertion is what is counted: a
#               chromosome with it is A or B, and one without it is O.
#   rs8176746   c.796C>A (p.Leu266Met), one of the four amino-acid changes
#               that make the B transferase (the same paper); T is counted.
#   rs41302905  c.802G>A (p.Gly268Arg), the O2 (O03) allele: an O with no
#               deletion, which would otherwise be read as A; T is counted.
#
# Their positions and forward-strand alleles come from Ensembl's REST API,
# grch37.rest.ensembl.org and rest.ensembl.org, and the run dies when Ensembl
# no longer lists the allele counted. Ensembl's VEP (rest.ensembl.org
# /vep/human/region, GRCh38, 2026-11-10) gave a frameshift, Leu to Met and Gly
# to Arg for these three alleles, at codons 86, 265 and 267 of its transcript
# ENST00000611156, which numbers them one lower than the names above do.
#
# RHD, the gene. Its span on each build, from Ensembl's lookup of the symbol,
# which is what a record of a deletion has to cover to be read as the RHD
# deletion.
#
# RHD, the sites. Most people who are RhD negative have no RHD gene on either
# chromosome, and a VCF of SNVs does not record a missing gene. What it does
# record is that somebody who has the gene differs from the reference at
# certain sites almost without fail, because the reference carries a rare
# allele there. Those sites are found here rather than taken from anywhere:
#
#   truth      1000 Genomes phase 3 genotyped the deletion as a copy-number
#              variant, 1:25592642 <CN0> (END=25661222, GRCh37). A sample's
#              count of <CN0> is how many chromosomes lack the gene; the 11
#              samples carrying <CN2> at that record are left out.
#   calls      The New York Genome Center's 30x resequencing of the same
#              people, jointly genotyped (GRCh38), in which a sample with no
#              gene is written "./." or with a depth of 0 inside it.
#   a site     is a biallelic SNV inside the gene at which at least
#              $Min_carriers of the samples with a gene have a heterozygous
#              or homozygous-ALT call of depth $Min_depth or more, and no
#              sample without one has.
#
# On 2026-11-10 that gave 15 sites from 2,332 SNVs and 2,493 samples (150 with
# no gene, 730 with one, 1,613 with two). The 15th had a call in 94.3% of the
# samples with a gene and the next best in 82.4%, so the cutoff of 0.90 falls
# in a gap. No sample without the gene had a call at any of the 15; of the
# 2,343 with it, 5 had a call at none of them and 2 at only one. Choosing the
# sites on a random half of the samples and counting on the other half, each
# way round, chose the same 15 and gave the same counts, so the result is not
# an artefact of choosing and counting on the same people.
#
# Each site's GRCh37 position is Ensembl's mapping of its GRCh38 one, and its
# rs number is the variant Ensembl has there with both alleles.
#
# blood_type() reads a marker that a VCF has no record of as no copies of the
# counted allele, which is only right when that allele is not the reference's
# own. So the run dies unless every counted allele differs from the reference
# base on both builds: ref_allele is the reference's on GRCh38 ("-" for the
# insertion), and each site's GRCh37 base is fetched from Ensembl and must be
# the same.
use 5.044;
use warnings FATAL => 'all';
use JSON::MaybeXS qw(decode_json);
use HTTP::Tiny;
use File::Spec;
use File::Temp qw(tempdir);

# htslib saves the index of a remote VCF in the working directory, so bcftools
# is run from a temporary one, which is left before File::Temp removes it.
my $Tmp = tempdir(CLEANUP => 1);
chdir $Tmp or die "can't chdir to $Tmp: $!\n";
END { chdir File::Spec->rootdir if defined $Tmp }

my %Rest = (b37 => 'https://grch37.rest.ensembl.org', b38 => 'https://rest.ensembl.org');
my $Phase3 = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz';
my $Deep = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.recalibrated_variants.vcf.gz';
my $Deletion = 25592642;  # POS of phase 3's <CN0> record over RHD, on GRCh37
my $Min_carriers = 0.90;  # the share of samples with the gene that must have a call at a site; see the header for the gap it falls in
my $Min_depth = 7;        # blood_type's default min_depth, which is sex_check's and peddy's

# [rs number, role, allele counted on the forward strand]
my @Abo = (['rs8176719', 'functional', 'C'], ['rs8176746', 'B', 'T'], ['rs41302905', 'O2', 'T']);

my $http = HTTP::Tiny->new(agent => 'Cn-Bioinf make_blood_type_panel.pl ', timeout => 120);

# Ensembl answers with 429 or a timeout when it is busy, so a GET is tried up
# to 5 times, a second longer apart each time.
sub get {
	my ($url) = @_;
	my $r;
	foreach my $try (1 .. 5) {
		$r = $http->get($url, { headers => { 'Content-Type' => 'application/json' } });
		return decode_json($r->{content}) if $r->{success};
		sleep $try;
	}
	die "GET $url failed: $r->{status} $r->{reason}\n";
}

# The records of $region of the VCF at $url, a line at a time, through
# bcftools; list-form, so that nothing goes through a shell.
sub records {
	my ($url, $region, $each) = @_;
	open my $fh, '-|', 'bcftools', 'view', '-r', $region, $url or die "can't run bcftools: $!\n";
	while (defined(my $line = <$fh>)) {
		next if $line =~ /\A##/;
		chomp $line;
		$each->([split /\t/, $line]);
	}
	close $fh or die "bcftools view -r $region $url failed: exit " . ($? >> 8) . "\n";
}

# ---- ABO

my @rows;
foreach my $marker (@Abo) {
	my ($rs, $role, $allele) = @{ $marker };
	my %at;
	foreach my $build (sort keys %Rest) {
		my $v = get("$Rest{$build}/variation/human/$rs");
		my ($m) = grep { $_->{seq_region_name} eq '9' } @{ $v->{mappings} || [] };
		die "$rs has no place on chromosome 9 of $build\n" unless $m;
		my ($ref, @other) = split m{/}, $m->{allele_string}; # Ensembl names the reference's allele first
		die "$rs: the counted allele $allele is not among its forward alleles $m->{allele_string} on $build\n"
			unless grep { $_ eq $allele } @other;
		# an insertion lies between two bases, which Ensembl gives as start =
		# end + 1; a VCF anchors it on the base before, which is "end"
		$at{$build} = { pos => $m->{start} > $m->{end} ? $m->{end} : $m->{start}, insertion => $m->{start} > $m->{end} ? 1 : 0, ref => $ref };
	}
	die "$rs is an insertion on one build and not on the other\n" if $at{b37}{insertion} != $at{b38}{insertion};
	die "$rs has reference allele $at{b37}{ref} on GRCh37 and $at{b38}{ref} on GRCh38\n" if $at{b37}{ref} ne $at{b38}{ref};
	push @rows, [$rs, 'ABO', $role, $at{b38}{ref}, $allele, $at{b37}{insertion} ? 'insertion' : 'allele', 9, $at{b37}{pos}, '', 9, $at{b38}{pos}, '', '', ''];
}

# ---- RHD, the gene

my %gene = map { $_ => get("$Rest{$_}/lookup/symbol/homo_sapiens/RHD") } sort keys %Rest;
foreach my $build (sort keys %gene) {
	die "RHD is not on chromosome 1 of $build\n" unless $gene{$build}{seq_region_name} eq '1';
}
push @rows, ['RHD', 'RHD', 'gene', '', '', 'region', 1, @{ $gene{b37} }{qw(start end)}, 1, @{ $gene{b38} }{qw(start end)}, '', ''];

# ---- RHD, the sites

# sample => how many of its chromosomes lack the gene; absent for a sample
# with a duplication
my %lacks;
{
	my @sample;
	records($Phase3, "1:$Deletion-$Deletion", sub {
		my ($f) = @_;
		if ($f->[0] eq '#CHROM') {
			@sample = @{ $f }[9 .. $#{ $f }];
			return;
		}
		my @alt = split /,/, $f->[4];
		my ($del) = grep { $alt[$_] eq '<CN0>' } 0 .. $#alt;
		return unless $f->[1] == $Deletion && defined $del;
		foreach my $i (0 .. $#sample) {
			my ($gt) = split /:/, $f->[9 + $i];
			my @allele = split m{[/|]}, $gt;
			next if grep { $_ ne '0' && $_ != $del + 1 } @allele;
			$lacks{ $sample[$i] } = grep { $_ ne '0' } @allele;
		}
	});
	die "phase 3 has no <CN0> record at 1:$Deletion\n" unless %lacks;
}

my @site;   # [GRCh38 position, REF, ALT, carriers with a call, samples with no gene with a call]
my ($carriers, $negatives, $snvs) = (0, 0, 0);
{
	my @class; # one per sample of the 30x file: 1 = has the gene, 0 = has none, undef = not in phase 3's truth
	records($Deep, "chr1:$gene{b38}{start}-$gene{b38}{end}", sub {
		my ($f) = @_;
		if ($f->[0] eq '#CHROM') {
			@class = map { defined $lacks{$_} ? ($lacks{$_} == 2 ? 0 : 1) : undef } @{ $f }[9 .. $#{ $f }];
			$carriers = grep { defined && $_ } @class;
			$negatives = grep { defined && !$_ } @class;
			return;
		}
		return unless $f->[3] =~ /\A[ACGT]\z/ && $f->[4] =~ /\A[ACGT]\z/;
		$snvs++;
		my @key = split /:/, $f->[8];
		my ($gt) = grep { $key[$_] eq 'GT' } 0 .. $#key;
		my ($dp) = grep { $key[$_] eq 'DP' } 0 .. $#key;
		return unless defined $gt;
		my @with = (0, 0); # samples with a call: [of those with no gene, of those with one]
		foreach my $i (0 .. $#class) {
			next unless defined $class[$i];
			my @v = split /:/, $f->[9 + $i];
			next unless defined $v[$gt] && $v[$gt] =~ m{\A(?:0[/|]1|1[/|]0|1[/|]1)\z};
			next if defined $dp && defined $v[$dp] && $v[$dp] =~ /\A[0-9]+\z/ && $v[$dp] < $Min_depth;
			$with[ $class[$i] ]++;
		}
		push @site, [$f->[1], $f->[3], $f->[4], $with[1], $with[0]] if $with[1] >= $Min_carriers * $carriers && !$with[0];
	});
}
die "no site of RHD passed; the 30x file's layout may have changed\n" unless @site;
warn "$snvs SNVs in RHD, $carriers samples with the gene and $negatives without: " . scalar(@site) . " sites chosen\n";

foreach my $s (@site) {
	my ($pos38, $ref, $alt, $with, $without) = @{ $s };
	my $map = get("$Rest{b38}/map/human/GRCh38/1:$pos38..$pos38:1/GRCh37");
	my ($m) = @{ $map->{mappings} || [] };
	die "1:$pos38 does not map to GRCh37\n" unless $m && $m->{mapped}{seq_region_name} eq '1' && $m->{mapped}{start} == $m->{mapped}{end};
	my $known = get("$Rest{b38}/overlap/region/human/1:$pos38-$pos38?feature=variation");
	my ($rs) = map { $_->{id} } grep {
		my %has = map { $_ => 1 } @{ $_->{alleles} || [] };
		$_->{id} =~ /\Ars[0-9]+\z/ && $_->{start} == $pos38 && $_->{end} == $pos38 && $has{$ref} && $has{$alt}
	} @{ $known };
	die "Ensembl has no rs number for 1:$pos38 $ref>$alt\n" unless defined $rs;
	my $pos37 = $m->{mapped}{start};
	my $base37 = get("$Rest{b37}/sequence/region/human/1:$pos37..$pos37:1")->{seq};
	die "$rs: GRCh37 has $base37 at 1:$pos37, where GRCh38 has $ref at 1:$pos38\n" unless uc $base37 eq $ref;
	push @rows, [$rs, 'RHD', 'present', $ref, $alt, 'allele', 1, $pos37, '', 1, $pos38, '', "$with/$carriers", "$without/$negatives"];
}

print "#name\tsystem\trole\tref_allele\tcount_allele\tkind\tb37_chrom\tb37_pos\tb37_end\tb38_chrom\tb38_pos\tb38_end\tcarriers\tnegatives\n";
print join("\t", @{ $_ }), "\n" foreach @rows;
