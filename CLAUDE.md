# Instructions for Claude in ~/Scripts/bioinf

These add to `~/.claude/CLAUDE.md`, and take precedence over it where they
differ.

## Error messages name the function that raised them

Every error or warning raised in `Basic.xs` or `lib/Bioinf/Basic.pm` starts
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
- Tests match messages with `like`, so a new prefix should not break them; if a
  test pins the start of a message, update it to include the prefix.

## README.md is the source of the POD

`md2pod.pl` replaces everything after the `1;` line of `lib/Bioinf/Basic.pm`
with `__END__` and POD converted from README.md, then checks the module with
`pod_file_ok()` and `Changes` with `changes_file_ok()`. It is a copy of
`~/Scripts/stats/md2pod.pl`, made on 2026-09-30, so a fix to the conversion
belongs in both copies.

- Document a function in README.md and run `perl md2pod.pl`; do not edit the
  POD in `lib/Bioinf/Basic.pm` directly, since the next run discards it.
- On 2026-09-30 the POD in `lib/Bioinf/Basic.pm` was still written by hand and
  said far more than README.md, which has empty sections for
  `get_best_alignment_hit` and `plot_phylo`. Running `md2pod.pl` then would
  have thrown that POD away. Move it into README.md first, and check with
  `git diff lib/Bioinf/Basic.pm` after a run that nothing was lost.
- `md2pod.pl` is an author-only helper that `[PruneFiles]` keeps out of the
  tarball. Like the helper scripts in `~/Scripts/stats`, it may use modern perl
  (`use 5.044`), whereas the module and `t/` are held to 5.10.

## Editing `Changes`

`Changes` is written by hand and is the only copy of the release notes; README.md
does not carry them. These rules are taken from `~/Scripts/stats/CLAUDE.md`:

- Add a release by prepending a section, and never re-word a release that has
  already shipped.
- A version line is `<version> <date> <tz>`, as in `0.01 2026-09-30 CDT`.
  `changes_file_ok()` in `md2pod.pl` fails on a release with no date or a
  version that does not parse, so run it after an edit.
- `$VERSION` in `lib/Bioinf/Basic.pm` is the version in progress. Work done
  while it reads a given version belongs under that version's section, so do
  not open a section for the next one.
