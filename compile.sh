#!/bin/sh
# Builds, tests and installs Bioinf::Basic, installing Alien::Bioinf (clustalo,
# BLAST+ and the plotting Python) from alien/Alien-Bioinf first if this perl does
# not have it yet. Once it is installed it is left alone, since reinstalling it
# rebuilds clustalo; update its tools with Alien::Bioinf->update instead.
#
# OPTIMIZE *replaces* perl's own $Config{optimize} (-O2 here) rather than adding
# to it, so bare -Wall builds at -O0. For Basic.xs that is the difference
# between 0.85 s and 0.81 s to read a 928 MB, 400,000-record FASTA (best of
# three, perl-5.44.0, 2026-09-29) -- small, since the reading is mostly
# memchr() and memcpy() inside libc, but not nothing. Keep -O2 ahead of the
# warning flags.
set -e
cd "$(dirname "$0")"
perl -MAlien::Bioinf -e1 2>/dev/null || cpanm ./alien/Alien-Bioinf
[ -f Makefile ] && make clean
perl Makefile.PL OPTIMIZE='-O2 -Wall' && make && make test && make install
