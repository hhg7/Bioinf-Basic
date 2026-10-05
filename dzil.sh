#!/bin/sh
# Build the Bioinf-Basic release tarball, for inspection before `dzil release`.
#
# Adapted on 2026-10-05 from ~/Scripts/stats/dzil.sh (commit a7869a7 there),
# which does the same for Stats-LikeR. This copy also fails, rather than only
# asking to be looked at, when a build product is in the tarball.
#
# Run it from this directory as `sh dzil.sh`. Like md2pod.pl it is an
# author-only helper, and dist.ini's [PruneFiles] keeps it out of the tarball.
set -e
cd "$(dirname "$0")"
perl md2pod.pl
# [Git::GatherDir] takes only what is committed, so the regenerated POD has to
# be committed before the build sees it. -a commits every tracked file that has
# changed, not only lib/Bioinf/Basic.pm, and never adds untracked ones (-A
# would); `|| true` because there is nothing to commit when the POD is current.
git commit -am "Update generated docs" || true
dzil clean
dzil build
tarball=$(ls Bioinf-Basic-*.tar.gz)
echo "==== $tarball ===="
tar tzf "$tarball"
# Products of compile.sh and of make (see .gitignore), coverage data, and a
# nested Bioinf-Basic-*/ from an earlier build gathered into this one. Each
# entry is <dist>-<version>/<path>, so the first component is stripped first.
if tar tzf "$tarball" | sed 's|^[^/]*/||' |
	grep -E '\.(c|o|obj|bs|dll|so|gcda|gcno)$|(^|/)(blib|_alien|\.build|cover_db|Bioinf-Basic-[^/]*)/|(^|/)(pm_to_blib|MYMETA\.[a-z]+)$'
then
	echo "dzil.sh: the build products above are in $tarball; it is not fit to release" >&2
	exit 1
fi
echo "If the rest of that list looks right, run: dzil release"
