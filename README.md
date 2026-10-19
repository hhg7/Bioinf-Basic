# SYNOPSIS

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

# DESCRIPTION

FASTA I/O in XS, BLAST hit ranking, and multiple-sequence-alignment plots and
LaTeX tables, taken from the maintainer's `bioinf.pm`; polygenic scores
from a VCF and a PGS Catalog scoring file; and, from VCFs, each sample's sex
and whether two people are parent and child.

Nothing is exported by default; ask for functions by name or with `:all`.
Every function dies (via `Carp::croak`) on bad arguments, with a message that
starts with the name of the function that raised it, as in
`fasta2hash: couldn't find "four" in x.fa`. `pgs_score`, `sex_check`,
`parent`, `plot_msa`, `plot_phylo`, `msa_quality_table` and `clustal_view_residues`
take `name => value` pairs, not a hash ref.

## Alien::Bioinf

Clustal Omega, BLAST+, PLINK 2 and the Python that draws the plots come from
`Alien::Bioinf`, which installs them module-locally and never uses `PATH`. It
is recommended rather than required, because it can only be installed where
NCBI builds BLAST+ (Linux and macOS on x86\_64 and aarch64, and Windows on
x86\_64) and it needs the network; PLINK 2 has no build for Linux on aarch64,
so there it installs the rest. `fasta2hash`, `hash2fasta_file` and
`get_best_alignment_hit` work without it; `pgs_score`, `plot_msa`,
`plot_phylo`, `msa_quality_table`, and `clustal_view_residues` given an
unaligned file, die saying that they need it. A CPAN client does not install a recommended module
unless asked, so install it first:

    cpanm Alien::Bioinf
    cpanm Cn::Bioinf

or, from a clone of [the repository](https://github.com/hhg7/Bioinf-Basic),
where it is `alien/Alien-Bioinf`:

    cpanm ./alien/Alien-Bioinf
    cpanm .

That downloads nothing it can reuse: archives are cached (`ALIEN_BIOINF_CACHE`,
default `~/.cache/alien-bioinf`), a clustalo or BLAST+ already on the machine at
the newest version is copied instead (`ALIEN_BIOINF_CLUSTALO`,
`ALIEN_BIOINF_BLAST`, or `PATH`), and pip reuses the wheels in its own cache.
The Python packages go into a venv of their own, which nothing installed into
the base Python can break. To check for and apply updates:

    perl -MAlien::Bioinf -MData::Dumper -e 'print Dumper(Alien::Bioinf->check_updates)'
    perl -MAlien::Bioinf -e 'print "$_\n" for Alien::Bioinf->update'

## Provenance in the images

Every PNG, SVG, PDF, PS or EPS image these functions draw carries its whole
provenance on one line, as its `Creator` metadata and nothing else, just as
Matplotlib::Simple writes it. The line names the calling script (its full
path), the function, this file and its
version, the Alien::Bioinf version, the user who ran it, the computer it ran
on (hostname and operating system), the perl version and path, and what drew
it: `msa_plot.py`, the Python version and path, the matplotlib version, and
for a tree the Biopython and NumPy versions too. The line has no date; an SVG
has the `<dc:date>` matplotlib writes beside it, taken from
`SOURCE_DATE_EPOCH` where that is set. After a `;` come the image's title
(or "untitled") and what it was made from and how, with the full path and
SHA-256 of each file:

- `plot_msa`: the FASTA file and how many sequences it holds (or the number
  of sequences in a hash ref), the ungapped copy of it that Clustal Omega
  read, with the sequences named `s1`, `s2` and so on, the alignment Clustal
  Omega wrote, with the names put back, the Clustal Omega version, the exact
  command, and the `msa_file` and `tree_file` kept, if any.
- `plot_phylo`: the same when it aligns, with the guide tree it drew;
  otherwise the newick file it drew. Then how many negative branch lengths
  were drawn as 0, if any; the tips, with the label each was shown as; and the
  whole newick tree that was drawn.
- `msa_quality_table`: the metric shown, and either the BLAST report it read
  or the FASTA file, its ungapped copy, the report blastp wrote, the blastp
  version, the exact command, and the `alignment_json` kept, if any.

For example:

    /home/me/work/run.pl called using "plot_phylo" in /.../Cn/Bioinf.pm
    version 0.01 with Alien::Bioinf 0.01 by user me on host myhost (linux)
    with Perl 5.44.0 (/usr/bin/perl), drawn by /.../msa_plot.py with
    Python 3.14.2 (/.../venv/bin/python),
    matplotlib 3.11.2, Biopython 1.87, NumPy 2.4.6; titled "DEG20010421";
    from the newick file /home/me/work/t.newick (SHA-256 f0e8...); 4 tips:
    S.cerevisiae, ...; the tree as drawn, in newick: (S.cerevisiae:0.389085,...);

An SVG holds it in `<dc:creator>`; `exiftool` or `identify -verbose` shows it
in a PNG or PDF.

# FUNCTIONS

## fasta2hash($file, $key)

Reads a FASTA file (gzip-compressed if its name ends in `.gz`, which needs a
`gzip` program on `PATH`). Returns a hash ref of defline (without the `>`) =>
sequence, or, with `$key`, just the sequence of that defline, reading no
further than the record after it. A defline that appears twice is warned about
and its sequences concatenated; with `$key`, only a repeat of `$key` is looked
for. Line endings may be `\n` or `\r\n`. Whitespace (space, tab, `\r`, `\f`
and `\x0b`) is dropped from sequences and from the end of deflines, as
Biopython drops it, so a defline of `x` and a space is found as `x`. Deflines and sequences are bytes,
just as they are in the file. A `.gz` that gzip cannot read to its end, such as
a truncated one, dies rather than returning the part that was read.

    my $h = fasta2hash('DEG20010421.fa');

## hash2fasta_file($hash, $filename, $order, $width)

Writes `$hash` as FASTA: the keys in `@$order` (default: sorted), sequences
wrapped at `$width` columns (default 80; 0 for one line each). Returns
`$filename`. Keys and sequences are written as bytes, as `print` to a `:raw`
handle writes them, so a string that perl holds as UTF-8 is written as the
characters in it. A character above 255 has no single byte and dies; encode
such a string first, as with `Encode::encode('UTF-8', $name)`.

## get_best_alignment_hit($json_file, $sort_criterion)

For a BLAST `-outfmt 15` JSON report, a hash ref of query title => array ref
of that query's hits, best first. Each hit is its best hsp plus `accession`,
`hit_len`, `id` and `title`, and `descriptions`, an array ref of every
description BLAST gave the hit. A hit against a non-redundant database such as
nr has one for each accession of the same sequence; `accession`, `id` and
`title` are those of the first, which is the one BLAST lists first. `$sort_criterion` is an hsp field (default
`evalue`): `align_len`, `bit_score`, `evalue`, `gaps`, `identity`,
`positive` or `score`. Ties are broken on the bit score, then on BLAST's
order. A report with two queries of one title dies, since only one of them
could be returned.

## pgs_score(%args)

Computes a polygenic score for every sample of a VCF from a scoring file of
the [PGS Catalog](https://www.pgscatalog.org/downloads/): for each variant,
the number of copies of the effect allele in the sample's genotype times
`effect_weight`, summed over the variants. [PLINK 2](https://www.cog-genomics.org/plink/2.0/),
which `Alien::Bioinf` installs, reads the genotypes and does the arithmetic
(`--score` with `no-mean-imputation`), so the genotypes are never held in perl:

    use Cn::Bioinf 'pgs_score';

    my $scores = pgs_score(vcf => 'genome.vcf.gz', weights => 'PGS000138.txt.gz');
    # { sample => { score => 1.6742, average => 6.97e-05, variants => 12015,
    #               alleles => 24030, missing => 10259, mismatched => 0 } }

Both arguments may also be data frames, as Stats::LikeR's `read_table`
returns them, which this module takes without loading Stats::LikeR; a VCF
given as a frame is written to a temporary VCF for PLINK 2. Files cost less
memory, since a frame is all in perl: on 2026-10-17, PGS000907's 1,773,528
variants were scored from files against a 3.5-million-record genome in 25 s,
with perl peaking at 0.8 GB.

- `vcf` (required): a `.vcf`, `.vcf.gz`, `.vcf.bgz` or `.bcf` file, a PLINK 2
  fileset given by its `.pgen` (with its `.pvar` and `.psam` beside it), or a
  data frame in any of `read_table`'s four shapes whose genotypes are columns
  named `<sample>.GT`, as `read_table` gives a VCF by default; one read with
  `explode => 0` has none, and dies. A frame needs `CHROM`, `POS`, `REF` and
  `ALT`, and `ID` too with `match => 'id'`.
- `weights` (required): a scoring file (`.txt` or `.txt.gz`), read a line at a
  time, or a data frame of one. It needs `effect_allele` and `effect_weight`,
  and either `chr_name` and `chr_position` or, with `match => 'id'`, `rsID`.
  `other_allele` is used when present. A harmonized file of the PGS Catalog
  (`..._hmPOS_GRCh38.txt.gz`) is read by `hm_chr` and `hm_pos` (or `hm_rsID`)
  instead, which are on the build its name gives. A variant marked TRUE in
  `is_haplotype`, `is_diplotype`, `is_interaction`, `is_dominant` or
  `is_recessive`, or given a `dosage_0_weight`, `dosage_1_weight` or
  `dosage_2_weight`, dies, since only additive weights are scored.
- `match`: `position` (the default) finds each variant by chromosome and
  position; `id` finds it by rsID, in a VCF whose `ID` column holds them. An
  `ID` of `rs1;rs2` is found by either.
- `build`: `hg19` (or `b37`) or `hg38` (or `b38`). PLINK 2 reads X only once
  told where its pseudo-autosomal regions are, so a score with variants on X
  or Y dies without it. Every sample is taken to be of unknown sex, which PLINK
  2 counts a haploid call on X of as homozygous: 2 copies, 2 alleles.
- `memory`, `threads`: the MiB of memory (640 or more) and the CPUs PLINK 2
  may use. By default it reserves half the machine's memory and uses every
  CPU.

Returns a hash ref keyed by sample, a hash of hashes that `write_table` can
write, with these fields for each sample:

- `score`: the sum of the weights times the dosages, PLINK 2's `SCORE1_SUM`.
  PLINK 2 writes it to six significant figures, so it is good to about one
  part in a million.
- `average`: `score` divided by `alleles`, PLINK 2's `SCORE1_AVG`, or undef
  when `alleles` is 0.
- `variants`: the variants of `weights` that were scored for this sample.
- `alleles`: the alleles counted in their genotypes, PLINK 2's `ALLELE_CT`:
  2 for each diploid call.
- `missing`: the variants that were not scored, because the VCF has no
  record of them whose alleles agree, or because the sample's genotype is
  missing (`./.`, or partly missing, such as `./1`). `variants` and `missing`
  add up to the number of rows of `weights`.
- `mismatched`: of `missing`, those at a position (or ID) the VCF has, whose
  alleles disagree with the record's there.

PLINK 2 by default fills a missing genotype in with the effect allele's
frequency in the samples; `pgs_score` leaves it out.

How a variant is matched, which perl decides, a `.pvar` line at a time,
before PLINK 2 scores:

- Chromosomes are compared without UCSC's `chr` prefix and in any case, and
  PLINK's `23`, `24` and `26` and UCSC's `M` are taken as `X`, `Y` and `MT`,
  so `chr1` in a VCF finds `1` in a scoring file.
- A variant agrees with a record when its effect allele and its other allele
  (if it has one) are both among the record's `REF` and `ALT` alleles, in any
  case. A multi-allelic record counts the copies of whichever `ALT` is the
  effect allele, so a `1/2` genotype of `G` to `A,T` holds one `T`. Where two
  records agree, as a multi-allelic site split in two can, the first in
  `CHROM:POS:REF:ALT` order is used.
- A record whose `ALT` is `.`, which a VCF of every site writes where the
  sample matches the reference, agrees with a variant either of whose alleles
  is its `REF`; a `0/0` there is two copies of the effect allele if it is the
  `REF`, and none otherwise.
- Alleles are not complemented, since a VCF and the PGS Catalog's files are
  both on the forward strand. A variant written on the other strand
  disagrees, and A/T and C/G variants, which read the same on either strand,
  are matched as they are written.
- Only `GT` is read: a dosage (`DS`) or a genotype likelihood is not used.

The positions in the scoring file have to be on the same genome build as the
VCF's, since a variant's position in GRCh37 (hg19) and in GRCh38 (hg38) differ;
rsIDs do not, so `match => 'id'` needs no build. Many scoring files say
`genome_build=NR` (not reported), so `pgs_score` cannot compare the builds.
The alleles say what it cannot: it warns when none of the variants was found
with alleles that agree, and when the alleles disagree at more than half of
the variants found. On 2026-10-17, PGS000138, whose header says
`genome_build=NR`, was scored against a GRCh37 genome: 12,015 of its 22,274
variants were found, and the alleles agreed at every one, so its positions
are on GRCh37. With every position moved on one base, the alleles agreed at 9
of the 58 found.

A VCF that lists only the variants a sample carries (as most single-sample
VCFs from a sequencing service do: every genotype holds an `ALT`) has no
record where the sample matches the reference, so those variants are
`missing`, not scored as zero copies. PLINK 2 refuses such a file; `pgs_score`
warns that its score is not comparable and scores it anyway. That genome above
was scored at 12,015 of PGS000138's 22,274 variants, and an exome at 701. A
score is comparable only with another scored on the same variants, so check
`variants` before comparing samples; a VCF called jointly over many samples,
or one of every site, has far fewer missing. The reference blocks of a gVCF,
which cover a range of positions with one record, are not expanded.

PLINK 2 converts the whole of a VCF or BCF before it keeps the score's
positions, which for the 12 GB, 81-million-record BCF of 1000 Genomes took 294
s and 6.2 GB on 2026-10-17. Cut the positions out first with `bcftools view
-R` (16 s), or convert such a reference once with `plink2 --make-pgen` and
give `pgs_score` the `.pgen`: from a `.pgen` of those positions, all 2,504
people were scored in under a second.

The PGS Catalog's search pages export a table of scores with a link to each
one's scoring file (`depression.pgs_scores_data.csv`, for one). That table
has no weights, so `pgs_score` cannot use it, but it can say which files to
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

## sex_check(%args)

Tells each sample of a VCF as male or female from its genotypes on X, by the
method of [peddy](https://github.com/brentp/peddy)'s `sex_check` (0.4.3): of
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

On 2026-10-18 a 155 MB `.vcf.gz` of a whole genome (3.5 million records)
took 4.8 s and 47 MB, where `zcat` alone takes 3.7 s; the same file as a
`.vcf.bz2` took 19.4 s, nearly all of it bzip2's.

- `vcf` (required): a `.vcf`, `.vcf.gz`, `.vcf.bgz` or `.vcf.bz2` file, the
  suffix in any case. It may also be given first, without its name, as
  `sex_check($file)` or `sex_check($file, report => 0)`. A VCF may name its chromosomes `chrX`, `X` or `23`
  (and `chrY`, `Y` or `24`); records on any other are skipped unread.
- `cutoff`: a sample whose `het_ratio` is below this is male, and otherwise
  female. The default, 0.6, is peddy's. On five genomes the three males' ratios
  were 0.029 (a whole genome), 0.103 and 0.201 (exomes), and the two
  females' 1.660 and 1.712 (exomes); an exome's off-target calls make more
  false heterozygotes than a genome's.
- `min_depth`: a genotype whose `DP` is below this is not counted. The
  default, 7, is peddy's; 0 counts every genotype. A genotype with no `DP`,
  or a `DP` of `.`, is counted.
- `report`: 1 (the default) prints the report; 0 prints nothing, and only
  returns the results.

What is counted, and where `sex_check` differs from peddy:

- X's pseudo-autosomal regions are taken as X:10000-2781479 and
  X:154931044-156030895, peddy's bounds, which cover GRCh37's and GRCh38's,
  so no build need be given. Y's are both builds' too.
- Only biallelic SNVs are counted: a record whose `REF` or `ALT` is not a
  single base (an indel, a multi-allelic site, or an `ALT` of `.`, `*` or
  `<NON_REF>`) is skipped, where peddy keeps an insertion. A record with no
  `GT` in its `FORMAT` is skipped too.
- Every site is counted, where peddy stops at 10,000, and peddy's filters on
  a site's call rate and ALT frequency across the cohort, which mean nothing
  for one genome, are left out.
- A haploid call, as a caller that knows a sample is male writes it on X,
  counts as homozygous: `1` as hom-ALT, `0` as hom-REF.
- A sample with no het or hom-ALT call on X has `sex` and `het_ratio` undef.
  peddy divides by a hom-ALT count of 0 as if it were 1, and so would call
  such a sample male. When the whole VCF has no SNV on X outside the
  pseudo-autosomal regions, `sex_check` also warns.

Returns a hash ref keyed by sample, a hash of hashes that `write_table` can
write, with these fields for each sample:

- `sex`: `male`, `female`, or undef.
- `het_ratio`: `x_het` divided by `x_hom_alt` (by 1 when `x_hom_alt` is 0),
  or undef when both are 0.
- `x_het`, `x_hom_alt`, `x_hom_ref`: the sample's heterozygous,
  homozygous-ALT and homozygous-REF calls at the SNVs counted on X. A VCF
  that lists only the sites where its sample carries an ALT allele has no
  hom-REF calls.
- `x_haploid`: of those, the haploid calls.
- `x_missing`: the sample's calls there that were missing (`./.`, `.`, or
  partly missing, such as `./1`), or not a call of one or two alleles of the
  record's two. A genotype left out by `min_depth` is in no count.
- `y_called`, `y_alt`: the sample's calls at SNVs on Y outside its
  pseudo-autosomal regions, and those of them that carry the `ALT`. They
  decide nothing: an exome of a woman had 530 called SNVs on Y, and one of a
  man 659, since a short read from X often maps to Y. They show at most
  whether Y was called at all.

A BCF or a `.pgen` is not read; convert it with `bcftools view -Oz`.

## parent(%args)

Tells whether two people are parent and child from their genotypes, and
returns 1 if they are and 0 if not. A child has one allele from each parent at
every site, so wherever a parent is homozygous the child carries that allele,
and the two are never opposite homozygotes (one hom-REF and the other
hom-ALT, which is called IBS0) except by a genotyping error. Genotypes cannot
tell which of the two is the parent, so `child` and `parent` may be given
either way round, with the same result. Paternity and maternity are this
together with the parent's sex, which `sex_check` tells.

Over the autosomal biallelic SNVs where either of the two carries the ALT, a
pair is parent and child when the IBS0 rate (the share of those sites that are
IBS0) is below `max_ibs0`, and the
[KING](https://www.kingrelatedness.com/manual.shtml)-robust kinship
coefficient is from 0.177 to 0.354, KING's range for first-degree relatives.
Full siblings have that kinship too, but they share no allele by descent at a
quarter of sites, so they have opposite homozygotes there, and the IBS0 rate
tells them apart. One person sequenced twice, or identical twins, have a
kinship above 0.354, and are not parent and child.

The file is decompressed by perl and counted in XS, as for `sex_check`. The
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

- `child`, `parent` (required): without `vcf`, the two VCF files, each of one
  sample (`.vcf`, `.vcf.gz`, `.vcf.bgz` or `.vcf.bz2`, as for `sex_check`);
  with `vcf`, the names of two of its samples.
- `vcf`: a VCF that holds both samples.
- `min_depth`: a genotype whose `DP` is below this is not counted. The
  default is 10; 0 counts every genotype. A genotype with no `DP`, or a `DP`
  of `.`, is counted. On the family above, counting every genotype gave the
  parent and child pairs IBS0 rates of 0.254 and 0.215, as high as an
  unrelated pair's, from calls of `DP` 1 to 3 made off the exome's targets. A
  `min_depth` of 7 gave 0.0053 and 0.0037, and 10 gave 0.0030 and 0.0020,
  keeping 92% of the sites of 7.
- `max_ibs0`: the IBS0 rate a parent and child must be below. The default,
  0.012, is [peddy](https://github.com/brentp/peddy) 0.4.3's (`ped_check`).
  It lies between the parent and child rates above (0.002 to 0.005) and the
  0.02 or so expected of full siblings, a quarter of the unrelated pair's 0.081.
- `report`: 1 (the default) prints the report to STDOUT; 0 prints nothing.

What is counted:

- Only the autosomes, named `1` to `22` or `chr1` to `chr22`; X is left out,
  since a father gives his X to a daughter and not to a son.
- Only biallelic SNVs: a record whose `REF` or `ALT` is not a single base, or
  that has no `GT`, is skipped.
- A call counts when it is diploid and of the record's `REF` and `ALT`, and
  is not below `min_depth`. A missing or partly missing call (`./.`, `./1`), a
  haploid one, or one with an allele the site lacks leaves that site out.
  So do sites where both are hom-REF, which say nothing about relatedness; a
  cohort's VCF has a great many of them.

Each sample's own VCF lists only the sites where its sample carries an ALT
allele, so a site that one file lacks is taken as hom-REF in that person, and
the two files are compared site by site, holding both in memory (8 bytes a
site, so about 30 MB for a genome). A position is left out when either file has more
than one SNV there (as after `bcftools norm -m-`), a call that does not count,
or a record of more than one `ALT`, or when the two files give it different
`REF` or `ALT` alleles. The files may name and order their chromosomes
differently.

Taking a missing site as hom-REF is right only when the other file was
sequenced and called there too, so two VCFs from one pipeline are what this
compares well. A genome's VCF has millions of sites that an exome's lacks, and
two exome capture kits target different regions; either makes IBS0 out of
sites that one of the two never had sequenced. The report gives the sites in
each file alone for that reason, and `parent` warns when one file has more
than twice as many as the other: the three pairs of the family above had
ratios of 1.13 at most, an exome of another kit had 3.0 to 13 against them,
and a genome against an exome 58 (with a kinship of -19.8). A VCF
that leaves out filtered calls has holes that look hom-REF too.

`parent` dies, rather than return a 0 or 1 that it has no grounds for, when

- fewer than 1,000 sites are compared, which at the unrelated pair's IBS0 rate
  would hold 81 opposite homozygotes on average, 8 standard deviations above
  the 12 that `max_ibs0` allows (an exome gives about 150,000 and a genome 3
  million);
- neither has a hom-REF call where the other carries the ALT, so that IBS0
  cannot be seen. A VCF merged from VCFs that list only their own sample's
  variants, as `bcftools merge` writes it without `-0`, is like that, with
  `./.` wherever a sample had no variant; give each sample's own VCF instead;
- either has no het call at the sites compared, so that no kinship can be
  taken.

In scalar context `parent` returns 1 or 0. In list context it returns that and
a hash ref of what it was decided from, with the same fields on every return:

- `parent`: 1 or 0, as returned.
- `relationship`: `parent-child`, `full siblings`, `duplicate` (one person or
  identical twins), `second-degree`, `third-degree` or `unrelated`, by
  KING's kinship ranges, with the IBS0 rate telling parent and child from
  full siblings.
- `kinship`: the KING-robust between-family estimate (Manichaikul et al. 2010,
  *Bioinformatics* 26:2867), which is what PLINK 2's `--make-king-table`
  gives, and which agreed with it exactly on the VCF of two of the exomes
  above merged by `bcftools merge -0`.
- `ibs0_rate`: `ibs0` divided by `sites`.
- `sites`: the sites compared, and of them `in_both`, `child_only` and
  `parent_only`: those in both own VCFs, in the child's alone and in the
  parent's alone. With `vcf`, every site is `in_both`.
- `ibs0`, `ibs1`, `ibs2`: the sites where the two share no allele, one, or
  both; and `ibs0` split into `child_alt_parent_ref` and
  `child_ref_parent_alt`.
- `het_het`: the sites where both are het; `child_het`, `parent_het`,
  `child_hom_alt` and `parent_hom_alt`: each one's het and hom-ALT calls at
  the sites compared.

A BCF is not read; convert it with `bcftools view -Oz`.

## plot_msa(%args)

Aligns sequences with Clustal Omega and draws the alignment. Returns a hash ref
of the files made: `filename`, `msa_file`, and `tree_file` if it was given.
Once the image is written it prints `wrote` and its file name to STDOUT, the
name in black on yellow when STDOUT is a terminal.

    plot_msa(
    	fasta      => 't/data/DEG20010421.fa',
    	filename   => 'DEG20010421.msa.svg',
    );

- `fasta`, `filename` (required): the sequences, as a FASTA file name or a
  hash ref of name => sequence (the function tells the two apart by whether it
  is a reference); and the image to draw, whose extension picks the format.
  A name may hold spaces and punctuation: Clustal Omega is given the sequences
  as `s1`, `s2` and so on, and the alignment and tree it writes get their
  names back, in UTF-8. In the tree, a name with whitespace or any of
  `()[]':;,` in it is quoted as newick quotes one, as `'it''s (x:1)'`.
- `msa_file`, `tree_file`: where to keep clustalo's alignment (FASTA; a
  temporary file otherwise) and its guide tree (newick; not made otherwise).
  Keep the tree to draw it with `plot_phylo` without aligning a second time.
- `order`: names, first to last (default: the input order; for a hash,
  sorted). Only these are drawn, and the first is drawn at the bottom.
- `labels`: a hash ref of name => label to show instead. matplotlib mathtext
  works: `'C.albicans' => '$\it{C. albicans}$'`. Two sequences drawn under one
  label die, since they would be one row of the image.
- `active_site_aa`, `query`: `{ His395 => 395, ... }`, a dashed vertical line
  at each of these 1-based residue numbers of the sequence named by `query`.
  A number below 1 dies.
- `title`, `xlabel`, `ylabel`, `threads`, `clustal_args`: the plot title; the
  axis labels (default "Amino Acid Residue" and "Protein & Species");
  clustalo threads (default 1); and an array ref of extra clustalo arguments.

With fewer than two sequences it warns and returns an empty hash ref.

## plot_phylo(%args)

Draws a guide tree, with Biopython's `Bio.Phylo`, from a FASTA (aligned with
Clustal Omega first) or from a newick file that `plot_msa` kept. Returns a
hash ref of the files made or used: `output_file`, and `tree_file` and
`msa_file` as below.

    plot_phylo(fasta => 't/data/DEG20010421.fa');   # writes phylo.svg
    plot_phylo(
    	'tree_file'   => 'DEG20010421.newick',
    	'output_file' => 'DEG20010421.tree.png',
    );

- `output_file`: the image to draw (default `phylo.svg`, in the working
  directory); the extension picks the format.
- `fasta`, `tree_file`: with `fasta` (as for `plot_msa`), the sequences are
  aligned with Clustal Omega and the guide tree drawn; `tree_file` and
  `msa_file` then say where to keep the tree and alignment, and `threads` and
  `clustal_args` are as for `plot_msa`. Without `fasta`, `tree_file` is an
  existing newick file to draw, such as one `plot_msa` kept, and no alignment
  is made.
- `labels`, `title`: a hash ref of name => label for the tips, as for
  `plot_msa`; and the plot title.

With `fasta` of fewer than two sequences it warns and returns an empty hash
ref.

## msa_quality_table(%args)

Draws an all-against-all BLAST score table with matplotlib, and returns
`filename`. The simplest call is

    msa_quality_table(
    	fasta      => 't/data/DEG20010421.fa',  # or { name => sequence }, aligned or not
    	filename   => 'DEG20010421.scores.png',
    );

- `fasta`, `filename` (required): the sequences, as a FASTA file name or a
  hash ref of name => sequence, as for `plot_msa`; they are aligned all
  against all with `blastp`. Gaps are stripped first, so an aligned FASTA,
  such as the one `plot_msa` keeps, will do. Names must look like
  `Genus.species[.strain]`. `fasta` is not needed when an existing
  `alignment_json` is given. `unaligned_fa` is its old name.
- `alignment_json`: the `blastp -outfmt 15` report, as a parsed hash ref or a
  file name. An existing file is read, and nothing is aligned; otherwise
  `blastp` is run on `fasta` and its report kept there for next time. Without
  it, the report is a temporary file.
- `metric`: the hsp field to show (default `score`).
- `normalize`: divide every value by the largest, after adding
  `logscale_add` to both, so the scale runs to 1. A largest value of 0 or
  less, as `evalue` has when BLAST has rounded every e-value to 0, dies.
- `order`, `logscale_add`, `default_undefined`, `title`, `cb_label`, `cb_min`,
  `cb_max`, `cblogscale`, `show_numbers`: the sequences to show, in order; a
  number added to every value; the value of a pair with no hit (otherwise
  drawn grey); the title; the colour bar's label, lower and upper ends, and
  whether it is logarithmic; and whether each cell shows its number.

`msa_file` is accepted and ignored, for old callers.

## clustal_view_residues(%args)

Writes an alignment as LaTeX tables with chosen residues coloured, and returns
`output_tex_file`, printing `wrote` and that file name to STDOUT, the name in
black on cyan when STDOUT is a terminal. Protein names are written so that
LaTeX prints them as they are, `_`, `^`, `{` and the like included.

- `msa_file` (required): a FASTA file. If its sequences are all one length
  (such as the alignment `plot_msa` keeps) it is shown as it is; if not, it is
  first aligned with Clustal Omega into a temporary file, and `msa_file`
  itself is never written.
- `output_tex_file` (required): the LaTeX file to write, meant to be
  `\input` into a document whose preamble loads `xcolor` (for the colours)
  and `placeins` (for the `\FloatBarrier` after each table).
- `color_residues`: `{ protein => { residue number => colour } }`, where
  residue numbers are 1-based and a colour is an xcolor name or
  `[r, g, b]`; a coloured column is coloured in every protein.
- `track`: a protein that gets a row under it showing its coloured residue
  numbers.
- `order`: an array ref of the proteins to show, top to bottom (default:
  sorted, ignoring case).
- `row_width`: alignment columns per block (default 100).
- `split`: blocks per LaTeX table (default 4); further tables are captioned
  "(continued)".
- `caption`: the table caption (default empty).
- `label`: written as `\label{tab:label}`, or as `tab:label0`, `tab:label1`,
  ... when there is more than one table.
- `table_text_size`: the LaTeX size command put at the start of each table
  (default `\footnotesize`).
- `threads`, `clustal_args`: clustalo threads (default 1) and an array ref of
  further clustalo arguments, when `msa_file` has to be aligned.

# Thanks

A lot of this work (not all!) used Claude AI, which was paid for by the University of Idaho's IMCI
