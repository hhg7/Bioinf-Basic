#!/bin/sh
# Commit, build and push the version in progress, then stop short of uploading:
# `dzil release` is left to be run by hand once the tarball has been checked.
#
# Adapted from ~/Scripts/stats/dzil.sh (commit b463c94 there) on 2026-11-18;
# until then this copy committed tracked files alone, under "Update generated
# docs", and so left new tests out of the tarball. Like md2pod.pl it is an
# author-only helper, and dist.ini's [PruneFiles] keeps it out of the tarball.
#
#   sh dzil.sh           # for a stable release, then: dzil release
#   sh dzil.sh --trial   # for a trial, then: dzil release --trial
set -e
# Every path below is relative to the repo root, so run from there whatever the
# caller's directory.
cd "$(dirname "$0")"

TRIAL=''  # '' = a stable release; '--trial' = a trial, its tarball named -TRIAL
case "${1-}" in
	'') ;;
	--trial) TRIAL='--trial' ;;
	*)
		echo "dzil.sh: the only option is --trial, not '$1'" >&2
		exit 1
		;;
esac

perl md2pod.pl

# The version in progress, used as the commit message, as ~/Scripts/stats does.
VERSION=$(perl -ne 'if (/\$VERSION\s*=\s*[\x27"]([^\x27"]+)/) { print $1; exit }' lib/Cn/Bioinf.pm)
if [ -z "$VERSION" ]; then
	echo "dzil.sh: could not read \$VERSION from lib/Cn/Bioinf.pm" >&2
	exit 1
fi

# md2pod.pl only checks that Changes parses. A newest section that names some
# other version still parses, and would ship as this release's notes.
CHANGES_VERSION=$(perl -ne 'if (/^(\d\S*)\s/) { print $1; exit }' Changes)
if [ "$CHANGES_VERSION" != "$VERSION" ]; then
	echo "dzil.sh: \$VERSION is $VERSION but the newest section of Changes is '$CHANGES_VERSION'" >&2
	exit 1
fi

# [Git::GatherDir] builds the dist only from files git tracks, so a new test left
# untracked is silently missing from the tarball, as t/popgen.t and
# t/admixture.t were from the first build of 0.03. Add the shipped paths by name
# rather than `git add *`, and exclude patch(1)'s *.orig and *.rej: PruneFiles
# keeps anything under lib/, t/ and share/, so a committed backup or reject
# there would ship, as lib/Stats/LikeR.pm.rej did in Stats::LikeR 0.3422.
# .gitattributes ships nothing, but it is what keeps a checkout from rewriting
# the CRLF lines of t/data/faidx.fa, whose bytes t/faidx.t checks.
git add -- lib t share Bioinf.xs Changes README.md CONTRIBUTING.md SECURITY.md \
	Makefile.PL dist.ini MANIFEST .gitattributes ':!*.orig' ':!*.rej'
# -a also picks up modified and deleted tracked files outside those paths.
# Commit only when something changed, rather than `git commit || true`, which
# also swallowed a commit that failed for a real reason.
if ! git diff-index --quiet HEAD --; then
	git commit -am "$VERSION"
fi

# `dzil clean` removes all of .build/, including the logs and work trees that
# test.all.perls.pl keeps in .build/multiperl, and every Cn-Bioinf-*/ and
# tarball of an earlier build. That is intentional: a release starts with them
# gone.
dzil clean
dzil build $TRIAL

TARBALL="Cn-Bioinf-$VERSION${TRIAL:+-TRIAL}.tar.gz"
# the directory every file of the tarball is under, which has no -TRIAL even in
# a trial's (Dist::Zilla 6, checked on 2026-11-18)
TOP="Cn-Bioinf-$VERSION"
if [ ! -f "$TARBALL" ]; then
	echo "dzil.sh: dzil build did not write $TARBALL" >&2
	exit 1
fi
# Build products, coverage output and editor/patch leftovers must never ship.
# The top-level directory (listed with no trailing slash) is dropped and its
# prefix stripped first, so a Cn-Bioinf-* left is a nested one: a dzil build
# that gathered an earlier build into this one.
BAD=$(tar tzf "$TARBALL" | sed -e "/^$TOP\/*\$/d" -e "s|^$TOP/||" | grep -E \
	'\.(c|o|obj|so|dll|dylib|bs|def|exp|gcda|gcno|gcov|orig|rej|swp)$|(^|/)(blib|_alien|cover_db|\.build|Cn-Bioinf-[^/]*)(/|$)|^(Makefile|Makefile\.old|MYMETA\.(json|yml)|pm_to_blib)$' \
	|| true)	# grep exits 1 when nothing matches, which is the good case
if [ -n "$BAD" ]; then
	echo "dzil.sh: $TARBALL contains files that must not ship:" >&2
	echo "$BAD" >&2
	exit 1
fi

# Pushed only once the build and the tarball check have both passed.
git push
echo "==== $TARBALL contents ===="
tar tzf "$TARBALL"
echo "If that looks right, run: dzil release${TRIAL:+ --trial}"
