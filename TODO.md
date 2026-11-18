# TODO

## Housekeeping

These two items were here before 2026-11-11.

- Make `dist.ini` like `~/Scripts/stats/dist.ini`, with all dependencies.
- Alter `Bioinf.xs` to use the `restrict` keyword where it will help.

## Basic techniques that Cn::Bioinf does not have

This is the list of 2026-11-11, less its first group (sequences held as
strings), which was written that day: `complement` to `motif_search`,
`fasta_index` and `fasta_region`. Each item names a package whose own tests
could supply the expected values, as `CLAUDE.md` asks under "Expected values
come from a reference and are frozen".

### Left over from the sequence functions

- Canonical k-mers for `kmer_count`: a k-mer and its reverse complement
  counted as one, as `jellyfish count -C` does.
- A minus strand for `subseq` and `fasta_region`, as `samtools faidx -i` gives.
- `fasta_region` on a bgzip-compressed FASTA, with its `.gzi` index (htslib
  `test/faidx` has the case).
- An index and regions for FASTQ, which `samtools fqidx` writes with a sixth
  column.
- Codon usage and the codon adaptation index (Biopython `test_SeqUtils.py`).
- The other ProtParam numbers: aromaticity, instability index, GRAVY and the
  molar extinction coefficient (Biopython `test_ProtParam.py`).
- The melting temperature of an oligonucleotide (Biopython
  `test_SeqUtils.py`, `MeltingTemp`).

### Other file formats

- FASTQ reading and writing, quality trimming and adapter removal (Biopython
  `SeqIO` tests; seqkit and cutadapt for trimming).
- GenBank and EMBL records (Biopython `test_GenBank.py`, BioPerl `t/SeqIO`).
- GFF3, GTF and BED parsing, and interval arithmetic: overlap, merge and
  intersect (bedtools' test suite).
- SAM, BAM and CRAM reading, and pileup or depth of coverage (samtools and
  htslib tests).
- General VCF work: filtering, merging, normalization and annotation. The
  shared `vcf_*` reader serves only the module's own analyses (bcftools'
  test suite).

### Alignment and search

- Pairwise alignment, global (Needleman-Wunsch) and local (Smith-Waterman),
  with substitution matrices (Biopython `test_pairwise_aligner.py`, EMBOSS
  `needle` and `water`).
- Running BLAST as a public function. `get_best_alignment_hit` ranks the hits
  of a report that already exists, and Alien::Bioinf already installs BLAST+.
- Read mapping, by running bwa or minimap2.
- Profile HMM searches against Pfam, by running HMMER.
- Primer design, by running Primer3.

### Phylogenetics and evolution

- Tree inference: neighbor joining, maximum likelihood and Bayesian, with
  bootstrap support. `plot_phylo` draws only Clustal Omega's guide tree.
- Distance matrices from an alignment (Biopython `test_TreeConstruction.py`).
- Substitution models, and dN/dS (Biopython `test_codonalign.py`, PAML).

### Genomics pipelines

- Variant calling from aligned reads.
- Genome assembly.
- Gene prediction and annotation.
- RNA-seq quantification and differential expression.
- Enrichment analysis against GO and KEGG.
- Population genetics was done on 2026-11-17: `pca`, `fst`, `hardy_weinberg`,
  `ld`, `ld_prune` and `gwas` by PLINK 2, and `admixture` in Bioinf.xs. Left:
  X in all of them, which needs each sample's sex; `ld_prune`'s list applied
  to the next function without a filtered VCF; ADMIXTURE's cross-validation
  error for choosing K (`--cv`), and its standard errors by bootstrap.

### Structure and databases

- Reading PDB and mmCIF files, RMSD, superposition and secondary structure
  assignment. `pdb_search` finds the structures for a sequence and does not
  read them (Biopython `test_PDB_*.py`).
- Fetching records from NCBI Entrez, UniProt and Ensembl, tested offline by
  replaying saved responses, as `t/pdb.t` replays `t/data/pdb_search.json`.
