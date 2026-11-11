# Instructions for Claude in ~/Scripts/bioinf

These add to `~/.claude/CLAUDE.md`, and take precedence over it where they
differ.

Cn::Bioinf is FASTA reading and writing in XS (`Bioinf.xs`), BLAST hit
ranking, alignment plots and tables drawn by running Clustal Omega, BLAST+
and `share/msa_plot.py` under the Python that Alien::Bioinf (the sibling
distribution in `alien/Alien-Bioinf`) installs, and polygenic scores
(`pgs_score`) computed by the PLINK 2 that Alien::Bioinf also installs (since
2026-10-17), and each sample's sex (`sex_check`), whether two people are
parent and child (`parent`), how every two people are related
(`relatedness`), and how well a test VCF's calls agree with a truth VCF's
(`concordance`) from VCFs, counted in XS by one shared VCF reader, the `vcf_*`
functions (`sex_check` and `parent` since 2026-10-18, `relatedness` and
`concordance` since 2026-10-19), and the genotypes of the HIrisPlex-S eye,
hair and skin colour panel (`hirisplex_genotypes`, read by the same VCF
reader) with a predictor for any multinomial logistic model the caller
supplies (`trait_predict`), both since 2026-10-20, and each sample's ABO and
RhD blood type (`blood_type`, since 2026-11-10, read by the same VCF reader
from the markers in `share/blood_type_panel.tsv`, which
`share/make_blood_type_panel.pl` chose from 1000 Genomes), and a search of
the RCSB PDB's structures for a sequence (`pdb_search`, since 2026-10-21),
which asks RCSB's search and data APIs over HTTPS with HTTP::Tiny and is
tested offline by replaying `t/data/pdb_search.json`. Since 2026-11-11 it
also has the everyday operations on a sequence held as a string (`complement`,
`reverse_complement`, `transcribe`, `back_transcribe`, `translate`,
`translate_six_frames`, `find_orfs`, `gc_content`, `kmer_count`,
`assembly_stats`, `molecular_weight`, `isoelectric_point`, `protein_charge`,
`subseq` and `motif_search`) and samtools's FASTA index (`fasta_index`,
`fasta_region`). Of those only the codon loop and the k-mer count are in XS;
the genetic codes are `share/genetic_codes.tsv`, which
`share/make_genetic_codes.pl` writes from NCBI's gc.prt. It is meant for CPAN,
so it has to work on perls and platforms that are not this one.

## No prediction model ships

`hirisplex_genotypes` prepares the HIrisPlex-S webtool's input and makes no
prediction, and `trait_predict` applies only a model the caller supplies.
That is deliberate and must stay so: the HIrisPlex-S coefficients are not
published (Chaitanya et al. 2018 cites them by reference; the webtool offers
no download and says no part of it may be reproduced without consent), and
part of the knowledge behind them is covered by European patent EP2195448. Do
not add a coefficient table to this distribution, and never write one from
memory -- an unsourced number here would be a defect twice over. What is
public, and what `share/hirisplex_s_panel.tsv` holds, is the panel itself,
taken from the webtool's own input form by
`share/make_hirisplex_panel.pl`.

On 2026-10-07 the rules below were merged in from the `CLAUDE.md` of three
sibling distributions: `~/Scripts/stats` (Stats::LikeR, XS),
`~/Scripts/SimpleFlow` (pure perl, runs commands), and
`~/Scripts/python/matplotlib/MatPlotLib-Simple` (perl that runs generated
Python). A section taken from one of them says which. Where a sibling's rule
was changed to fit this repo, the section gives the reason. When a sibling
learns something new, check whether it belongs here too.

## Error messages name the function that raised them

Every error or warning raised in `Bioinf.xs` or `lib/Cn/Bioinf.pm` starts
with the name of the function that raised it, a colon and a space:
`fasta2hash: couldn't find "four" in x.fa`. This covers `croak`, `carp`, `die`
and `warn` in perl and `croak()` and `warn()` in C, including new ones.

- In perl, a public function that keeps its own name in `my $sub = '...'`
  writes `"$sub: ..."`; any other function writes its name literally.
- A private helper names itself (`_open_fasta: ...`), except one that is handed
  the caller's `$sub` to report on that caller's arguments (`_args`, `_align`),
  which prefixes `$sub` instead, since the fault is in the caller's arguments.
- In C, the prefix is `__func__` (`croak("%s: ...", __func__, ...)`), which
  names the C function (`read_fasta`, `write_fasta`) and cannot drift out of
  step with it.
- Name the offending value in the message: the key, the index, the file, the
  exit code. This rule comes from SimpleFlow. `plot_msa: active site X must be a
  residue number, 1 or more, not "0"` is the model.
- Tests match messages with `like`, so a new prefix should not break them; if a
  test pins the start of a message, update it to include the prefix.
- Call every XSUB through `_located(sub { ... })` (since 2026-11-10). C's
  `croak()` gives the line of `lib/Cn/Bioinf.pm` that called the XSUB, and
  `_located` raises the error again at the caller's line, as `croak` in perl
  does. A test that pins the location allows for no period after the line
  number (`line \d+\.?\n`), since Carp before 1.25, perl 5.12.5's, writes none.

## Generated files: edit the source, never the copy

### README.md is the source of the POD

`md2pod.pl` replaces everything after the `1;` line of `lib/Cn/Bioinf.pm`
with `__END__` and POD converted from README.md, then checks the module with
`pod_file_ok()` and `Changes` with `changes_file_ok()`. It is a copy of
`~/Scripts/stats/md2pod.pl`, made on 2026-09-30, so a fix to the conversion
belongs in both copies.

- Document a function in README.md and run `perl md2pod.pl`; do not edit the
  POD in `lib/Cn/Bioinf.pm` directly (no `Edit`, `Write`, `sed -i` or patch),
  since the next run discards it. Fix a POD problem in the Markdown it came
  from.
- On 2026-09-30 the POD in `lib/Cn/Bioinf.pm` was still written by hand and
  said far more than README.md, which has empty sections for
  `get_best_alignment_hit` and `plot_phylo`. Running `md2pod.pl` then would
  have thrown that POD away. Move it into README.md first, and check with
  `git diff lib/Cn/Bioinf.pm` after a run that nothing was lost.
- `md2pod.pl` is an author-only helper that `[PruneFiles]` keeps out of the
  tarball. Like the helper scripts in `~/Scripts/stats`, it may use modern perl
  (`use 5.044`), whereas the module and `t/` are held to 5.12.

### Other build output

- `Cn-Bioinf-*/` and `Cn-Bioinf-*.tar.gz` are `dzil build` output, and
  `read.me.pod` is written by `md2pod.pl`. All are in `.gitignore`. Do not edit
  them, and leave them out of any repo-wide search-and-replace.
- `cover_db/` and the `*.gcov`, `*.gcda` and `*.gcno` files in the root are
  coverage output and are tracked in git, but are stale (see "Coverage" below).
  Regenerate them; never edit them.
- The `*.pl` and `*.tex` files in the root, other than the helpers named in
  this file, are the maintainer's scratch and example scripts. `[PruneFiles]`
  keeps them out of the tarball.

## Editing `Changes`

`Changes` is written by hand and is the only copy of the release notes; README.md
does not carry them, and nothing regenerates it, so a bad edit is a lost
release note. These rules are taken from `~/Scripts/stats/CLAUDE.md` and
`~/Scripts/SimpleFlow/CLAUDE.md`:

- Add a release by prepending a section, and never re-word a release that has
  already shipped.
- A version line is `<version> <date> <tz>`, as in `0.01 2026-09-30 CDT`.
  `changes_file_ok()` in `md2pod.pl` fails on a release with no date or a
  version that does not parse, so run it after an edit, or load the file
  through `CPAN::Changes`.
- Entries go under ` [Bracketed headings]`, with one leading space before the
  heading. Bullets are written ` - text`, continuation lines are indented three
  spaces, and lines wrap at the width the 0.01 section uses. The 0.01 section
  is the model to follow.
- An entry says what was wrong and what the user will now observe, not "fixed a
  bug". When it is about a platform, name the platform and the report that
  found it. This rule comes from Matplotlib::Simple.
- `$VERSION` in `lib/Cn/Bioinf.pm` is the version in progress. Work done
  while it reads a given version belongs under that version's section, so do
  not open a section for the next one. Check what has actually shipped (on
  PAUSE, not in the local tarballs) before assuming a version is released.

## Perl in `lib/` and `t/`

### House style

- Indent with tabs. Perl comments are `#`. A block comment is a run of `#`
  lines directly above the code, and a single fact about one line trails it.
  The `/* … */` layout in the global `CLAUDE.md` is for C only; the rest of
  that comment doctrine applies to perl too.
- Comment every enumerated value where it is declared, as the global doctrine
  says. For an option key, that means its type and meaning.
- Resolve each option exactly once, then read only the resolved copy, e.g.
  `$r{'threads'} = $args->{'threads'} // 1;`. Never test the raw
  `$args->{...}` again for an option that has a default. SimpleFlow 0.14
  silently ignored a documented default that way.
- Booleans returned to a caller are 0 or 1. They are never undef and never the
  caller's own string.
- Prefer a core module to a new dependency. `_wrote` writes its ANSI escapes
  itself rather than loading `Term::ANSIColor`. A new prereq must be justified,
  added to `dist.ini`, and installable on `perl-5.12.5`.

### Callers run under `use warnings FATAL => 'all'`

Every test does, and so do many callers, which turns each "uninitialized
value" warning from a nuisance into a crash in the caller's program. SimpleFlow
shipped three fixes of exactly this kind (0.13, 0.14 and 0.16).

- Never let a filetest, `length` or a comparison see an undef. Validate a name
  before it is filetested. `-s $file` on a missing file is undef, so default it
  or test for existence first.
- A returned record has the same keys on every return path. A field that
  exists on only one path is a bug.

### No `autodie` in the module; check every builtin yourself

The module does not use `autodie`, and that rule comes from Matplotlib::Simple,
where removing it uncovered a lost Python traceback. Every `open`, `close`,
`binmode`, `mkdir`, `unlink`, `rename` and `system` checks its own result and
`croak`s with the function name and `$!`, as `_json_file`, `_python` and `_run`
do. Never use `autodie ':all'` in a test either: it loads `IPC::System::Simple`
at compile time, which is not core and not a prereq.

### Running external programs

clustalo, blastp, plink2 and the plotting Python come from Alien::Bioinf, never
from PATH, and are run with list-form `system`: by `_run`, or for plink2 by
`_plink2`, which keeps PLINK 2's output off the caller's terminal and puts its
`Error:` lines in the message. `pgs_score` hands PLINK 2 files and reads its
results a line at a time, so that genotypes are never held in perl. Some of
the rules below come from Matplotlib::Simple and SimpleFlow, and some describe
what `lib/Cn/Bioinf.pm` already does:

- Never use the one-argument `system`, backticks or `qx//` to run anything. Each
  of those goes through a shell (`cmd.exe` on Windows), so a path containing a
  space or a backslash breaks the command.
- Keep double quotes out of every argument. Perl on Win32 wraps an argument
  that contains a space in double quotes, but does not escape a double quote
  already inside it. That is why `_python` passes its arguments in a JSON
  `--argfile`; anything new that goes to Python goes the same way.
- Never paste a filesystem path or caller-supplied text into Python source: a
  Windows path's `\b` is a backspace in a Python literal, and an apostrophe
  closes it. Pass such values as data through the argfile.
- Decode `$?` fully. `-1` means the program never started (`$!` says why), the
  low 7 bits give the signal, and the high byte gives the exit code. Read the
  signal before shifting.
- A feature that needs `fork`, process groups or signals is POSIX-only. It must
  refuse on `MSWin32` with a message that says why, and must not quietly do
  nothing. Never wrap `system` in `$SIG{ALRM}`, because that leaves the child
  running as an orphan.

### JSON::MaybeXS picks whichever backend is installed

This section comes from Matplotlib::Simple, whose 0.318 broke on it. Cpanel::JSON::XS
and JSON::PP write a dualvar (a scalar that is both a number and a string, such as
`scalar @empty`) as a number, but JSON::XS 4.04 writes it as a string. A number
bound for `encode_json` or `_json_text` must therefore be built as a number
(`0 + $x`). Run any change to what reaches the encoder under JSON::XS too,
with Cpanel::JSON::XS hidden:

    echo 'unshift @INC, sub { die "hidden\n" if $_[1] eq "Cpanel/JSON/XS.pm" }; 1;' > /path/HideXS.pm
    PERL5OPT="-I/path -MHideXS" prove -Ilib t/

### Compatibility of the interface

- The functions' argument keys and the fields they return are the public API.
  New behaviour that would change what an existing caller gets defaults to off.
  Removing, renaming, or changing the type of a key or returned field is an
  incompatible change, and the reply must say so, so that it reaches `Changes`.
- `$VERSION` stays a quoted string. A bare `0.20` stringifies as `0.2`, which
  CPAN would treat as older than `0.15`. `dist.ini` takes the version from the
  module (`[VersionFromModule]`), so the module is the one place to set it.

### Prereqs

`[Prereqs]` in `dist.ini` lists what `lib/Cn/Bioinf.pm` loads, core modules
included. Anything only `t/` needs goes under `[Prereqs / TestRequires]`, and
anything only the `Makefile.PL` probe needs goes under
`[Prereqs / ConfigureRequires]`. Alien::Bioinf stays under `RuntimeRecommends`,
for the reason the comment there gives. The author-only scripts declare
nothing.

## C in `Bioinf.xs`

This section is taken from `~/Scripts/stats/CLAUDE.md`, whose `LikeR.xs` is
the house style. The only floating-point arithmetic in `Bioinf.xs` is
`concordance`'s ratios (`conc_ratio()`, since 2026-10-19), one `NV` division
each, which `t/concordance.t` compares within a relative 1e-15 so that
long-double and quadmath perls pass; follow the `NV` rules for anything
more. The sequence functions' floating point (`gc_content`,
`molecular_weight`, `isoelectric_point`, `protein_charge`) is in perl on
purpose (2026-11-11): perl's own `**` and `/` are already the width of an
`NV`, so none of the `nv_*` libm block was needed.

### Types match the value's real domain

Plain `int` is not a default.

- Two states only → `bool`, assigned `TRUE`/`FALSE`, which are perl's own (do
  not add `<stdbool.h>`).
- Cannot be negative → unsigned and wide: `size_t` for sizes, lengths, counts
  and indices; `STRLEN` for a length crossing the perl string API; `UV` (or
  `IV` when signed) for an integer to or from perl; `uint64_t`/`uint32_t` when
  the algorithm depends on the width; and `unsigned char` for a byte.
- A small enumerated set → `short int` (signed when `-1` means "not yet
  decided"), with each value commented at the declaration.
- A loop counter bounded by a small literal → `unsigned short int`; bounded
  by a runtime count → `size_t`; bounded by `av_len()`/`AvFILLp()` →
  `SSize_t`.
- Floating point → `NV`, never `double`. Use `NV_INF`, `NV_NAN` and
  `NV_EPSILON`, never call libm bare (copy the `nv_*` macro block from
  `LikeR.xs` together with its link probe), and format with
  `my_snprintf` and `NVgf`.

Changing a width must not silently change a perl-visible conversion: fix the
`Sv*` calls and the format strings in the same edit, and keep the build free of
warnings, including mixed-sign comparisons.

### `restrict`, and when not to use it

Pointer parameters and locals get `restrict`, with `const` when the pointee is
not written, unless they may alias. The header comment of `Bioinf.xs` explains
why only `read_fasta()`'s block buffer carries it: every other pointer there is
perl's, and a defline and its sequence can be COW copies of one buffer. Leave
it off, with a short comment saying why, for overlapping buffers, re-pointed
pointers, anything reachable by another route, perl-managed memory, and
anything that escapes. Keep the `#if !defined(__cplusplus) &&
!defined(restrict)` block below the includes. Do not churn existing
declarations just to add or remove `restrict`.

### Loop variables live inside their loop

Declare a counter in the `for` and a per-iteration temporary in the innermost
block that uses it. Each loop gets its own counter. The exception is a value
the code reads after the loop, which is declared before it with a comment
saying so. In perl the same rule is `for my $i` and a `my` inside the loop
body. Apply this to new code and to functions already being edited; do not
churn the rest.

### Portable C

There is no local Windows, Solaris or BSD perl, so this is discipline applied
while writing; no test run here will catch it.

- Allocate with `Newx`/`Newxz`/`Renew`/`Safefree`, never `malloc`/`free`. Reach
  for the perl API (`PerlIO`, `my_snprintf`, `sortsv`, `strEQ`) before libc.
- No GNU extensions: no VLAs, nested functions, statement expressions,
  `typeof`, unguarded `__builtin_*`, zero-length arrays, case ranges, or `void
  *` arithmetic. Check with `-std=c99`, not `gnu99`, when in doubt.
- No POSIX-only headers or functions (`<unistd.h>`, `strcasecmp`, `fork`), and
  no glibc-only or BSD-only ones (`qsort_r`, `memmem`, `strndup`, `getline`,
  `asprintf`, `strchrnul`, `reallocarray`, `arc4random`).
- Printf lengths are perl's: `%" UVuf "`, `%" IVdf "`, `%" NVgf "`, never
  `%zu`, `%lld` or `%llu`.
- Never cast a `char *` to a wider pointer and dereference it, because SPARC
  and some ARM systems fault on unaligned access.
- Before calling any libc function listed in `reentr.h`, check what perl's core
  does with it. On a threaded perl, the core gets the `_r` form while an XS
  module gets the plain one, which put Stats::LikeR 0.316's `srand` out of step.
- If a Solaris or POSIX declaration ever goes missing under strict C99, add
  the `_GNU_SOURCE`/`__EXTENSIONS__` block from the top of `LikeR.xs`. If a
  perl API newer than 5.12 is used, add `ppport.h`; `Bioinf.xs` has needed
  neither so far.
- The C99 flag is probed by the `[MakeMaker::Awesome]` header in `dist.ini`.
  Change the probe there, not only in the hand-kept `Makefile.PL`.

`perl xs.check.pl` runs XS::Check and the house rules over `Bioinf.xs`; run it
after any edit there.

## Tests

### Expected values come from a reference and are frozen

This section adapts the rule in `~/Scripts/stats/CLAUDE.md` that tests come
from R's and SciPy's own suites. Here the references are
bioinf.pm's original pure-perl functions (the source of `t/fasta.t` and
`t/blast.t`), and the programs and libraries the function imitates:
Biopython's `SeqIO` for FASTA (installed under
`/home/con/.pyenv/versions/3.14.2/lib/python3.14/site-packages/Bio`, although
its `Tests/` directory is not installed), and BLAST+ and Clustal Omega
themselves.

For the sequence functions (2026-11-11) the references are other packages'
own test suites, fetched from their repositories at a tag, and the programs
themselves: Biopython 1.87 (`Tests/test_seq.py`, `test_translate.py`,
`test_SeqUtils.py`, `test_ProtParam.py`, from the `biopython-187` tag, since
`Tests/` is not installed), scikit-bio 0.7.4, BioPerl (bioperl-live at commit
4a43b46), htslib 1.25 (`test/faidx`), assembly-stats 1.0.1, EMBOSS 6.6.0 and
samtools 1.19.2. `t/seq.t`, `t/orf.t` and `t/faidx.t` say which value is
whose. EMBOSS and samtools are not installed here: for
`t/data/emboss_reference.pl` and `t/data/faidx_reference.pl` they were unpacked
from Ubuntu's packages into a scratch directory, with no root needed
(`apt-get download emboss emboss-lib emboss-data libhpdf-2.3.0`, then `dpkg -x`
each, and `LD_LIBRARY_PATH`, `EMBOSS_ACDROOT` and `EMBOSS_DATA` pointed into
the unpacked tree; `samtools libhts3t64 libhtscodecs2 libdeflate0` likewise).
A test of a new sequence function takes its cases from one of those suites
before anything is made up for it.

- Expected values are literals in the `.t` file, or are fixtures in `t/data/`.
  A generator that is worth keeping sits next to its output, with a comment
  saying how to re-run it, as `t/data/make_fixtures.pl` does. The tests never
  call it.
- Each file's header comment records where every expected value came from:
  which reference and version, which function, and the date it was checked.
  "what bioinf.pm's pure-perl fasta2hash returned for the same input (checked
  2026-09-29)" is the standard.
- Where the module deliberately differs from its reference, the test says so
  at that case.
- Choose a tolerance or limit from a measurement and say what was measured.
  Never widen one to make a failure go away.

### What a smoker has, and what it may skip

A CPAN smoker has the prereqs in `dist.ini` and nothing else: no network, no
Alien::Bioinf (which is not on CPAN), and often no gzip and no Python.

- Anything that needs Alien::Bioinf's tools is skipped when they are missing,
  as `t/msa.t` and `t/pgs.t` do with `plan skip_all`. A missing dependency must
  produce a SKIP and never a FAIL.
- Everything that does not need those tools must still run there. Argument
  checks, and anything else that is decided before clustalo, blastp, plink2 or
  Python starts, go in `t/checks.t`, because that is the part a Windows smoker
  actually exercises. `t/seq.t`, `t/orf.t` and `t/faidx.t` need none of the
  tools, so they run in full on every smoker.
- `t/data/faidx.fa` is htslib's file byte for byte, CRLF lines included, and
  `fasta_index`'s offsets are bytes: `t/faidx.t` checks its SHA-256 first, and
  `.gitattributes` marks it `-text` so that a checkout cannot rewrite it.
- Never skip a test because a reference implementation (R, Python, Biopython)
  is missing, since expected values are frozen and a skipped cross-check never
  runs on a smoker.
- A test must never build a shell command line. Use list-form `system($^X,
  ...)`, put the program in a file or pass it with `-e`, and pass paths as
  arguments. `t/fasta.t` makes its `.gz` fixtures with `IO::Compress::Gzip`
  for this reason, and skips its `.gz` tests only when `fasta2hash` reports
  that it can't run the gzip program.
- No hardcoded `/tmp`, no `/` pasted into a path, and no `which`, `ls` or `cp`.
  Use `File::Temp` and `File::Spec`, as the existing tests do. Matplotlib::Simple
  0.312 failed 131 of 166 Windows subtests because of one `/tmp`.
- A child perl that loads the module gets this working copy's `@INC`, as
  `t/checks.t` passes `@inc`, and never a bare `-MCn::Bioinf`.
- Never depend on mtimes or `time` advancing during a `sleep`: a VM's clock can
  lag, and mtimes can be whole seconds. Check contents, or set times apart with
  `utime`.

### Regression tests

The rules in this section come from SimpleFlow.

- Write the test first, and confirm that it fails against the code that had the
  bug (`git stash` the fix, or unpack the previous tarball) before the fix goes
  in. Use only arguments that the old code accepted.
- Name the old behaviour in the test's description, as `t/checks.t` does with
  "dies, rather than counting from the end", so that a re-break can be
  recognised from the failure message alone.
- Never assert only that something is empty or that no warning was emitted.
  First assert a positive sentinel showing that the probe ran.

### Leak tests and old perls

`t/fasta.t` uses `Test::LeakTrace` where it is installed. A perl's own bug
broke Stats::LikeR's leak tests (another, on 5.10.0, is below this module's
minimum):

- On 5.32, regcomp leaks on every run-time-compiled pattern. Never interpolate
  a variable into a match that runs more than once per value in a leak-tested
  function; cache the compiled `qr//` instead.

### Coverage

Aim to exercise every path: each call form, each option and its default, every
`croak`/`carp` message, every validation branch, and each returned field. The
committed `cover_db/` and the `.gcov` files are stale, so regenerate them with
`Devel::Cover` (and `gcov` for `Bioinf.xs`) before quoting a number. Gate
genuinely slow cases behind `EXTENDED_TESTING` rather than dropping them.

## Every change must hold across the support matrix

`dist.ini` declares `perl = 5.012`, the Makefile.PL's `MIN_PERL_VERSION` is
`5.012`, and `lib/Cn/Bioinf.pm` says `require 5.012`. The minimum was 5.10
until 2026-10-21, when the maintainer raised it to match Alien::Bioinf, which
cannot be installed on 5.10 (see below). `test.all.perls.pl` reads it from
`dist.ini` and leaves older perls out.

### Back to perl 5.12

- 5.12 syntax only: no signatures, no postfix dereference, no `s///r`, no
  `package NAME BLOCK`, no `keys`/`values`/`each` on a reference, no `__SUB__`,
  `fc`, `isa` or lexical subs. `//` and `state` are 5.10, and `package NAME
  VERSION`, `...` and `keys`/`values`/`each` on an array are 5.12, but `say`
  needs `use feature 'say'`.
- Tests keep the existing header: `require 5.012; use strict; use warnings
  FATAL => 'all';`. `md2pod.pl`, `test.all.perls.pl`, `xs.check.pl` and the
  other helpers do not ship and may use modern perl.
- Raising the minimum is the maintainer's decision, never a way to make an
  error go away. If a change truly needs a newer perl, say so in the reply and
  stop.

### The perls installed here

On 2026-10-07, `/home/con/perl5/perlbrew/perls/` held `perl-5.44.0` (the
default), `perl-5.42.3` (threaded), `perl-5.32.1`, `perl-5.12.5` (long double),
`perl-5.10.1`, `5.44.0-quadmath`, `5.44.0-i686` (32-bit, so `IV_MAX` is
2147483647), and `5.16.3-thr-ld` (threaded, long double, older than 5.20). By
2026-11-02 `perl-5.18.4` had been added; it has no Alien::Bioinf, so
`test.all.perls.pl` will try to install it there. On 2026-11-10
`5.44.0-i686-64int` was there too, and neither it nor `perl-5.18.4` had
File::ShareDir::Install, so both failed at `perl Makefile.PL` under
`./test.all.perls.pl --no-install --no-alien`; `--deps` would install it.
On 2026-11-11 the same command gave the same result: seven perls passed
(2,342 tests where Alien::Bioinf is installed, 1,760 on `5.44.0-i686`) and
those two failed at `perl Makefile.PL` for that reason.

`./test.all.perls.pl` builds and tests against all of them, installing
Alien::Bioinf into any perl that lacks it, and with `-p` against chosen ones.
Run it for anything that touches `Bioinf.xs` or `lib/Cn/Bioinf.pm`.

On 2026-10-21 every perl had Alien::Build, File::ShareDir::Install and the
rest of the CPAN prereqs (installed by `--deps`), and six had Alien::Bioinf:
`perl-5.44.0`, `perl-5.42.3`, `perl-5.32.1`, `5.44.0-quadmath`, `perl-5.12.5`
and `5.16.3-thr-ld`. The last two needed Alien::Bioinf's Makefile.PL to ask for
HTTP::Tiny 0.018 (5.16.3's core 0.017 got "599 SSL connection failed" from
GitHub) and IO::Uncompress::Unzip 2.030 (5.12.5's core 2.024 never ended the
loop over the plink2 zip, and the install was killed after 30 minutes); its
comment there gives the details. The other two cannot install it here:

- `perl-5.10.1`: Alien::CPython3 needs File::Copy::Recursive, whose tests need
  Test::Deep, and Test::Deep 1.205 requires perl 5.12. So since 2026-10-21
  Alien::Bioinf requires 5.012 (`use 5.012` heads its Makefile.PL, so a 5.10
  smoker reports NA), and so does Cn::Bioinf. `test.all.perls.pl` skips this
  perl, which stays installed only as a perl below the minimum.
- `5.44.0-i686`: Net::SSLeay does not compile, since this machine has no
  32-bit OpenSSL headers (`openssl/opensslconf.h`).

Test `5.44.0-i686` with `./test.all.perls.pl --no-alien`; `t/msa.t` and
`t/pgs.t` then skip, as on a CPAN smoker.

A failed Alien::Bioinf install keeps its build tree in `.build/multiperl/work/`,
and one that got as far as its venv holds symlinks to the same `python3`.
MakeMaker's kit check walks the whole tree and dies on them, silently, with
exit 255, so every later `perl Makefile.PL` in the root fails, on every perl,
until that tree is removed. The parallel mode builds in private copies that
leave `.build/` out, and is not affected.

Because the IV width and the NV width are independent, any cast from a size or
an `NV` to `IV` must be gated on `IV_MAX`, and `5.44.0-i686` is the perl that
shows when it is not.

When a smoker report fails, find out what was different about that smoker
before reaching for another perl. SimpleFlow's 0.193 and 0.194 failures came
from signals inherited ignored and from a lagging clock, not from the perl
version.

### CI

`.github/workflows/test.yml` was untracked on 2026-10-07. It builds the
tarball with dzil and tests it on Linux, macOS and Strawberry Perl on Windows,
without Alien::Bioinf and with a space in `TMPDIR`. Until it is committed and
running, the only Windows signal is a CPAN Testers report after the release.

## Releasing

`sh dzil.sh` runs `md2pod.pl`, **commits every changed tracked file with
`git commit -am`**, then runs `dzil clean` and `dzil build`, and fails if a
build product reached the tarball. It stops before `dzil release`. So:

- Do not run `sh dzil.sh` without being asked, since it commits.
- Never run `dzil release`. Do not commit or push unless asked.
- Anything added to the root that should not ship must be kept out by the
  `[PruneFiles]` pattern in `dist.ini`; `[Git::GatherDir]` takes every tracked
  file otherwise.

## Keep the facts in this file true

The installed perls, which prereqs each perl has, and the dates and counts
quoted here change without anyone editing this file. Check a fact on the
machine before relying on it in a recommendation. If it is wrong, correct it
in the same session and give the date it was checked.
