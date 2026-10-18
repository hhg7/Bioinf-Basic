# SYNOPSIS

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

    # a polygenic score for each sample of a VCF, from a PGS Catalog scoring file
    use Stats::LikeR 'read_table';
    my $scores = pgs_score(
    	vcf     => read_table('genome.vcf.gz'),
    	weights => read_table('PGS000138.txt.gz', sep => "\t"),
    );

# DESCRIPTION

FASTA I/O in XS, BLAST hit ranking, and multiple-sequence-alignment plots and
LaTeX tables, taken from the maintainer's `bioinf.pm`, and polygenic scores
from a VCF and a PGS Catalog scoring file.

Nothing is exported by default; ask for functions by name or with `:all`.
Every function dies (via `Carp::croak`) on bad arguments, with a message that
starts with the name of the function that raised it, as in
`fasta2hash: couldn't find "four" in x.fa`. `pgs_score`, `plot_msa`,
`plot_phylo`, `msa_quality_table` and `clustal_view_residues` take
`name => value` pairs, not a hash ref.

## Alien::Bioinf

Clustal Omega, BLAST+ and the Python that draws the plots come from
`Alien::Bioinf`, which installs them module-locally and never uses `PATH`. It
is recommended rather than required, because it can only be installed where
NCBI builds BLAST+ (Linux and macOS on x86\_64 and aarch64, and Windows on
x86\_64) and it needs the network. `fasta2hash`, `hash2fasta_file`,
`get_best_alignment_hit` and `pgs_score` work without it; `plot_msa`,
`plot_phylo`, `msa_quality_table`, and `clustal_view_residues` given an
unaligned file, die saying that they need it. A CPAN client does not install a recommended module
unless asked, so install it first:

    cpanm Alien::Bioinf
    cpanm Bioinf::Basic

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

    /home/me/work/run.pl called using "plot_phylo" in /.../Bioinf/Basic.pm
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
`effect_weight`, summed over the variants. Both arguments are data frames, as
Stats::LikeR's `read_table` returns them, so this module never reads either
file itself and does not need Stats::LikeR:

    use Stats::LikeR 'read_table';
    use Bioinf::Basic 'pgs_score';

    my $weights = read_table('PGS000138.txt.gz', sep => "\t");
    my $vcf = read_table('genome.vcf.gz');
    my $scores = pgs_score(vcf => $vcf, weights => $weights);
    # { sample => { score => 1.674, average => 6.97e-05, variants => 12015,
    #               alleles => 24030, missing => 10259, mismatched => 0 } }

A scoring file needs `sep => "\t"`: `read_table` takes a tab from the name of
a `.tsv` or `.vcf`, but a `.txt.gz` would be read as commas. Its `#` header
lines are skipped as comments.

- `vcf` (required): the VCF, in any of `read_table`'s four shapes, as long as
  each sample's genotype is a column of its own named `<sample>.GT`.
  `read_table` gives a `.vcf`, `.vcf.gz` or `.vcf.bgz` those columns by
  default; a VCF read with `explode => 0` has none, and dies. It needs
  `CHROM`, `POS`, `REF` and `ALT` (or `ID`, `REF` and `ALT` with
  `match => 'id'`).
- `weights` (required): the scoring file, in any of the four shapes. It needs
  `effect_allele` and `effect_weight`, and either `chr_name` and
  `chr_position` or, with `match => 'id'`, `rsID`. `other_allele` is used
  when present. A harmonized file of the PGS Catalog (`..._hmPOS_GRCh38.txt.gz`)
  is read by `hm_chr` and `hm_pos` (or `hm_rsID`) instead, which are on the
  build its name gives. A variant marked TRUE in `is_haplotype`,
  `is_diplotype`, `is_interaction`, `is_dominant` or `is_recessive`, or given
  a `dosage_0_weight`, `dosage_1_weight` or `dosage_2_weight`, dies, since
  only additive weights are scored.
- `match`: `position` (the default) finds each variant by chromosome and
  position; `id` finds it by rsID, in a VCF whose `ID` column holds them. An
  `ID` of `rs1;rs2` is found by either.

Returns a hash ref keyed by sample, a hash of hashes that `write_table` can
write, with these fields for each sample:

- `score`: the sum of the weights times the dosages.
- `average`: `score` divided by `alleles`, or undef when `alleles` is 0.
- `variants`: the variants of `weights` that were scored for this sample.
- `alleles`: the alleles counted in their genotypes, 2 for each diploid call.
- `missing`: the variants that were not scored, because the VCF has no
  record of them whose alleles agree, or because the sample's genotype is
  missing (`./.`, or partly missing, such as `./1`). `variants` and `missing`
  add up to the number of rows of `weights`.
- `mismatched`: of `missing`, those at a position (or ID) the VCF has, whose
  alleles disagree with the record's there.

These are the `SCORE1_SUM`, `SCORE1_AVG` and `ALLELE_CT` that PLINK 2 writes
with `--score ... no-mean-imputation cols=+scoresums`. PLINK 2 by default
fills a missing genotype in with the effect allele's frequency in the samples
instead; `pgs_score` leaves it out.

How a variant is matched:

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
- Only `GT` is read: a dosage (`DS`) or a genotype likelihood is not used, and
  a haploid call such as `1` counts its one allele.

The positions in the scoring file have to be on the same genome build as the
VCF's, since a variant's position in GRCh37 (hg19) and in GRCh38 (hg38) differ;
rsIDs do not, so `match => 'id'` needs no build. `read_table` skips the
`#genome_build=` header line, and many files say `NR` (not reported) there
anyway, so `pgs_score` cannot compare the builds. The alleles say what it
cannot: it warns when none of the variants was found with alleles that agree,
and when the alleles disagree at more than half of the variants found. On
2026-10-17, PGS000138, whose header says `genome_build=NR`, was scored against
a GRCh37 genome: 12,015 of its 22,274 variants were found, and the alleles
agreed at every one, so its positions are on GRCh37. With every position
moved on one base, the alleles agreed at 9 of the 58 found.

A VCF that lists only the variants a sample carries (as most single-sample
VCFs from a sequencing service do: every genotype holds an `ALT`) has no
record where the sample matches the reference, so those variants are
`missing`, not scored as zero copies. That genome above was scored at 12,015
of PGS000138's 22,274 variants, and an exome at 701. A score is comparable
only with another scored on the same variants, so check `variants` before
comparing samples; a VCF called jointly over many samples, or one of every
site, has far fewer missing. The reference blocks of a gVCF, which cover a
range of positions with one record, are not expanded.

A whole-genome VCF read into memory takes several GB. `read_table`'s `filter`
can keep only the positions the score needs while the file is read; the genome
above then took 9 s and under 1 GB:

    my %pos = map { $_->{chr_position} => 1 } @{ $weights };
    my $vcf = read_table('genome.vcf.gz', filter => { POS => sub { $pos{$_} } });

The PGS Catalog's search pages export a table of scores with a link to each
one's scoring file (`depression.pgs_scores_data.csv`, for one). That table
has no weights, so `pgs_score` cannot use it, but it can say which files to
download and score:

    # $vcf read without a filter, since each score has variants of its own
    my $catalog = read_table('depression.pgs_scores_data.csv');
    foreach my $pgs (@{ $catalog }) {
    	my ($file) = $pgs->{'Scoring File (FTP Link)'} =~ m{([^/]+)\z};
    	next unless -f $file;
    	my $s = pgs_score(vcf => $vcf, weights => read_table($file, sep => "\t"));
    	printf "%s\t%s\t%g\t%d variants\n", $pgs->{'Polygenic Score ID & Name'}, $_, @{ $s->{$_} }{qw(score variants)}
    		foreach sort keys %{ $s };
    }

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
