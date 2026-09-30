# Synopsis

FASTA I/O in XS, BLAST hit ranking, and multiple-sequence-alignment plots and
LaTeX tables, taken from `~/Scripts/bioinf.pm`. See the POD in
`lib/Bioinf/Basic.pm`.

Clustal Omega, BLAST+ and the Python that draws the plots are installed with
it, module-locally, by `Alien::Bioinf` in `alien/Alien-Bioinf`, which has to be
installed first:

    cpanm ./alien/Alien-Bioinf
    cpanm .

That downloads nothing it can reuse: archives are cached (`ALIEN_BIOINF_CACHE`,
default `~/.cache/alien-bioinf`), a clustalo or BLAST+ already on the machine at
the newest version is copied instead (`ALIEN_BIOINF_CLUSTALO`,
`ALIEN_BIOINF_BLAST`, or `PATH`), and pip installs only what the base Python
lacks. To check for and apply updates:

    perl -MAlien::Bioinf -MData::Dumper -e 'print Dumper(Alien::Bioinf->check_updates)'
    perl -MAlien::Bioinf -e 'print "$_\n" for Alien::Bioinf->update'

# Functions/Subroutines

## clustal_view_residues

View sequences aligned in a latex file

### Arguments

| Argument        | Default | Meaning |
|-----------------|---------|-------------------------------------------------------------------------|
| `msa.file`        | (required) | A FASTA file. If its sequences are all one length (such as the alignment `plot_msa` keeps) it is shown as it is; if not, it is first aligned with Clustal Omega into a temporary file, and `msa.file` itself is never written |
| `output.tex.file` | (required) | The LaTeX file to write, meant to be `\input` into a document; it is also the return value |
| `color.residues`  | none       | `{ protein => { residue number => colour } }`: residue numbers are 1-based, and a colour is an xcolor name or `[r, g, b]`; a coloured column is coloured in every protein |
| `track`           | none       | A protein that gets a row under it showing its coloured residue numbers |
| `order`           | sorted, ignoring case | An array ref of the proteins to show, top to bottom |
| `row.width`       | 100        | Alignment columns per block |
| `split`           | 4          | Blocks per LaTeX table; further tables are captioned "(continued)" |
| `caption`         | empty      | The table caption |
| `label`           | none       | Written as `\label{tab:label}`, or as `tab:label0`, `tab:label1`, ... when there is more than one table |
| `table.text.size` | `\footnotesize` | The LaTeX size command put at the start of each table |
| `threads`         | 1          | Threads for clustalo, when `msa.file` has to be aligned |
| `clustal.args`    | none       | An array ref of further clustalo arguments, when `msa.file` has to be aligned |


## fasta2hash

Read a fasta file to a hash

    my $h = fasta2hash('DEG20010421.fa');

## get_best_alignment_hit

## hash2fasta_file

Write a FASTA to a file.

## msa_quality_table

    msa_quality_table(
    	fasta      => 't/data/DEG20010421.fa',  # or { name => sequence }, aligned or not
    	filename   => 'DEG20010421.scores.png',
    );

## plot_msa

    plot_msa(
    	fasta      => 't/data/DEG20010421.fa',
    	filename   => 'DEG20010421.msa.svg',
    );

## plot_phylo

Draw a guide tree with Biopython's `Bio.Phylo`, from a FASTA (aligned with
Clustal Omega first) or from a newick file that `plot_msa` kept. The image is
`output.file`, `phylo.svg` by default:

    plot_phylo(fasta => 't/data/DEG20010421.fa');   # writes phylo.svg
    plot_phylo(
    	'tree.file'   => 'DEG20010421.newick',
    	'output.file' => 'DEG20010421.tree.png',
    );

The image records its own provenance: in an SVG's Dublin Core `<metadata>` (or
a PNG's text chunks) are the date to the second, the input files and the
alignment with each one's SHA-256, the exact clustalo command, the newick tree
drawn and its SHA-256, and every program and library that took part, with its
version.

# COPYRIGHT AND LICENSE

This software is free.  It is licensed under the same terms as Perl itself

# Thanks

A lot of this work (not all!) used Claude AI, which was paid for by the University of Idaho's IMCI
