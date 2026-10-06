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

# DESCRIPTION

FASTA I/O in XS, BLAST hit ranking, and multiple-sequence-alignment plots and
LaTeX tables, taken from the maintainer's `bioinf.pm`.

Nothing is exported by default; ask for functions by name or with `:all`.
Every function dies (via `Carp::croak`) on bad arguments, with a message that
starts with the name of the function that raised it, as in
`fasta2hash: couldn't find "four" in x.fa`. `plot_msa`, `plot_phylo`,
`msa_quality_table` and `clustal_view_residues` take `name => value` pairs,
not a hash ref.

## Alien::Bioinf

Clustal Omega, BLAST+ and the Python that draws the plots come from
`Alien::Bioinf`, which installs them module-locally and never uses `PATH`. It
is recommended rather than required, because it can only be installed where
NCBI builds BLAST+ (Linux and macOS on x86\_64 and aarch64, and Windows on
x86\_64) and it needs the network. `fasta2hash`, `hash2fasta_file` and
`get_best_alignment_hit` work without it; `plot_msa`, `plot_phylo`,
`msa_quality_table`, and `clustal_view_residues` given an unaligned file, die
saying that they need it. From this repository, install it first:

    cpanm ./alien/Alien-Bioinf
    cpanm .

That downloads nothing it can reuse: archives are cached (`ALIEN_BIOINF_CACHE`,
default `~/.cache/alien-bioinf`), a clustalo or BLAST+ already on the machine at
the newest version is copied instead (`ALIEN_BIOINF_CLUSTALO`,
`ALIEN_BIOINF_BLAST`, or `PATH`), and pip installs only what the base Python
lacks. To check for and apply updates:

    perl -MAlien::Bioinf -MData::Dumper -e 'print Dumper(Alien::Bioinf->check_updates)'
    perl -MAlien::Bioinf -e 'print "$_\n" for Alien::Bioinf->update'

## Provenance in the images

Every PNG, SVG, PDF, PS or EPS image these functions draw carries its whole
provenance on one line, as its `Creator` metadata and nothing else, just as
Matplotlib::Simple writes it. The line names the calling script (as the
working directory plus the script's name), the function, this file and its
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
  read, the alignment Clustal Omega wrote, the Clustal Omega version, the
  exact command, and the `msa_file` and `tree_file` kept, if any.
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

Reads a FASTA file (gzip-compressed if its name ends in `.gz`). Returns a hash
ref of defline (without the `>`) => sequence, or, with `$key`, just the
sequence of that defline, reading no further than the record after it. A
defline that appears twice is warned about and its sequences concatenated;
with `$key`, only a repeat of `$key` is looked for. Line endings may be `\n`
or `\r\n`. A `.gz` that gzip cannot read to its end, such as a truncated one,
dies rather than returning the part that was read.

    my $h = fasta2hash('DEG20010421.fa');

## hash2fasta_file($hash, $filename, $order, $width)

Writes `$hash` as FASTA: the keys in `@$order` (default: sorted), sequences
wrapped at `$width` columns (default 80; 0 for one line each). Returns
`$filename`.

## get_best_alignment_hit($json_file, $sort_criterion)

For a BLAST `-outfmt 15` JSON report, a hash ref of query title => array ref
of that query's hits, best first. Each hit is its best hsp plus `accession`,
`hit_len`, `id` and `title`. `$sort_criterion` is an hsp field (default
`evalue`): `align_len`, `bit_score`, `evalue`, `gaps`, `identity`,
`positive` or `score`. Ties are broken on the bit score, then on BLAST's
order. A report with two queries of one title dies, since only one of them
could be returned.

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
  `\input` into a document.
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
