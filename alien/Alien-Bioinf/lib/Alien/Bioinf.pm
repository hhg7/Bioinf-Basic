package Alien::Bioinf;
# Module-local Clustal Omega, NCBI BLAST+, PLINK 2, bcftools, Prodigal,
# kallisto, SPAdes and a Python venv for Cn::Bioinf.
#
# Where each tool comes from, as checked on 2026-09-29:
#
#   BLAST+     https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/,
#              ncbi-blast-<v>+-<platform>.tar.gz beside a <file>.md5 that is
#              checked before anything is unpacked. 2.17.0 was current.
#   clustalo   https://www.clustal.org/omega/ binaries: clustalo-<v>-Ubuntu-x86_64
#              and -Ubuntu-32-bit, clustal-omega-<v>-macosx, and
#              clustal-omega-<v>-win64.zip -- 1.2.4, 1.2.4, 1.2.3 and 1.2.2.
#              The site sits behind Cloudflare, which answered 403 to every
#              request from the machine this was written on, so wherever that
#              page cannot be read, or has no binary for the platform, the
#              version comes from the tags of github.com/GSLBiotech/clustal-omega
#              (the upstream repository; 1.2.4 was the newest) and it is built
#              from that tag's tarball against argtable2-13, the release
#              clustalo's configure asks for, fetched from SourceForge.
#   PLINK 2    https://www.cog-genomics.org/plink/2.0/ (added 2026-10-17), the
#              zips in the page's stable column, which link to
#              s3.amazonaws.com/plink2-assets/alpha<n>/plink2_<platform>_<date>.zip:
#              alpha 7.13 of 15 Oct 2026 was current. Each zip holds plink2,
#              vcf_subset and Intel's licence; only plink2 is kept. Nothing on
#              the page gives a checksum, so the zip is checked only by
#              unpacking it and by the build date plink2 --version reports.
#              The builds are linux_x86_64 and linux_i686, mac and mac_arm64,
#              and win64 and win32; there is none for Linux on aarch64. The
#              AVX2 builds are faster but stop on a CPU without AVX2, so the
#              plain ones are used. A version such as 2.0.0-a.7.13 doesn't
#              compare as a number, so the build date stands for it.
#   Python     Alien::CPython3 supplies the interpreter; a venv on top of it
#              gets numpy, matplotlib, adjustText and biopython -- what
#              share/msa_plot.py and Matplotlib::Simple import -- and their
#              dependencies (scipy, for adjustText) from PyPI, and, since
#              2026-11-13, pydeseq2 (0.5.4 then, which brings pandas, anndata
#              and scikit-learn), for Cn::Bioinf's differential_expression.
#              The newest release for check_updates() is the highest
#              "Python 3.x.y" in
#              https://www.python.org/api/v2/downloads/release/.
#
# The four below were added on 2026-11-13, for Cn::Bioinf's call_variants,
# predict_genes, rnaseq_quant and assemble. Each is released on GitHub, and
# the newest release is the tag that https://github.com/<repository>/releases/latest
# redirects to, which is read from the redirect rather than from GitHub's API,
# since the API answers only 60 requests an hour from one address. None of the
# four publishes a checksum beside its files, so each is checked by the version
# the installed program reports.
#
#   bcftools   github.com/samtools/bcftools: bcftools-<v>.tar.bz2 (1.25 was
#              current), which holds HTSlib too, and is built here, as there
#              are no binaries: sh ./configure --disable-libcurl, then make,
#              which needs a C compiler and zlib's headers, and bzip2's and
#              liblzma's for the CRAM files that use them. Where those two are
#              missing it is built without them. There is no Windows build.
#   Prodigal   github.com/hyattpd/Prodigal: prodigal.linux (x86_64),
#              prodigal.osx.10.9.5 (x86_64) and prodigal.windows.exe, of
#              v2.6.3, from 2016, and the newest. Elsewhere, and should a newer
#              release name its files otherwise, it is built from the tag's
#              source with make, which needs only a C compiler.
#   kallisto   github.com/pachterlab/kallisto: kallisto_linux-v<v>.tar.gz
#              (x86_64) and kallisto_mac-v<v>.tar.gz, which in v0.52.0 is an
#              arm64 binary (file(1), 2026-11-13). Releases up to v0.50.1 had
#              a Windows zip; v0.52.0 has none, and there is no source build
#              here, since that needs CMake and autoconf.
#   SPAdes     github.com/ablab/spades: SPAdes-<v>-Linux.tar.gz (x86_64),
#              -Darwin-x86_64 and -Darwin-arm64. 4.3.0's Linux tarball is 196
#              MB and unpacks to 684 MB, of which 499 MB is programs that
#              spades.py does not run to assemble reads (pathracer, spaligner,
#              binspreader, splitter, and the spades-g*, -kmer*, -read-filter
#              and -convert-bin-to-fasta tools), the SARS-CoV-2 tables and
#              gtest. What is kept is 185 MB: bin/spades.py and
#              spades_init.py, spades-core, spades-hammer, spades-ionhammer,
#              spades-corrector-core and spades-bwa, and share/spades without
#              sewage/ and coronaspades_hmms/. spades.py is run under this
#              distribution's own Python, so it is installed after that.
#
# Nothing is downloaded twice. Every archive is kept in cache_dir() and reused
# from there; a clustalo or BLAST already on this machine at the wanted version
# is copied instead of fetched (see _local_copy); and pip reuses the wheels in
# its own cache.
#
# The venv does not see the base interpreter's packages. Up to 0.01 it was made
# with --system-site-packages, so pip skipped whatever the base Python already
# had, and the venv then broke whenever that changed: on 2026-10-13 a numpy
# 1.26.4 installed into this machine's base Python left its scipy 1.18.1, which
# needs numpy 2, unable to import, and adjustText with it.
#
# ALIEN_INSTALL_NETWORK=0, the Alien::Build setting, stops every download, pip's
# included (see _http and _pip_index).
use strict;
use warnings;
use parent 'Alien::Base';
use Carp qw(croak);
use Config ();
use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path remove_tree);
use File::Copy qw(copy);
use File::Basename qw(basename dirname);
our $VERSION = '0.02';

my @PY_PACKAGES = qw(numpy matplotlib adjustText biopython pydeseq2);
# the name each of @PY_PACKAGES is imported by, where that is not its PyPI name
my %PY_IMPORT = (biopython => 'Bio');
my @TOOLS = qw(clustalo blast plink2 bcftools prodigal kallisto spades python python-packages);
# The tools that are built for fewer platforms than BLAST+ is. One that cannot
# be had is left out, and the function of Cn::Bioinf that needs it says so.
my %OPTIONAL = map { $_ => 1 } qw(plink2 bcftools prodigal kallisto spades);
# For each tool released on GitHub: its repository; the arguments that make it
# print its version; and the pattern whose first group is that version, in
# bcftools's "bcftools 1.25", Prodigal's "Prodigal V2.6.3: February, 2016"
# (which it writes to STDERR), kallisto's "kallisto, version 0.52.0" and
# spades.py's "SPAdes genome assembler v4.3.0".
my %GITHUB = (
	bcftools => { repo => 'samtools/bcftools', args => ['--version'], says => qr/^bcftools (\d[\w.]*)/m },
	prodigal => { repo => 'hyattpd/Prodigal', args => ['-v'], says => qr/^Prodigal V(\d[\w.]*)/m },
	kallisto => { repo => 'pachterlab/kallisto', args => ['version'], says => qr/^kallisto, version (\d[\w.]*)/m },
	spades => { repo => 'ablab/spades', args => ['--version'], says => qr/^SPAdes genome assembler v(\d[\w.]*)/m },
);
my $BLAST_URL = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/';
my $CLUSTAL_URL = 'https://www.clustal.org/omega/';
my $CLUSTAL_GIT = 'https://api.github.com/repos/GSLBiotech/clustal-omega/tags';
my $ARGTABLE_URL = 'https://downloads.sourceforge.net/project/argtable/argtable/argtable-2.13/argtable2-13.tar.gz';
my $PY_URL = 'https://www.python.org/api/v2/downloads/release/?is_published=true&pre_release=false';
my $PLINK2_URL = 'https://www.cog-genomics.org/plink/2.0/';
my $WIN = $^O eq 'MSWin32';
my $EXE = $WIN ? '.exe' : '';

# ---- paths to the installed tools -------------------------------------------

# absolute, since a relative @INC entry gives a relative dist_dir, and a venv
# python run by a relative path warns that its sys.prefix is unexpected
sub _root { File::Spec->rel2abs($_[1] // $_[0]->dist_dir) }

sub clustalo { my ($c, $root) = @_; catfile($c->_root($root), 'clustalo', $c->_ver('clustalo', $root), 'bin', "clustalo$EXE") }

sub blast {
	my ($c, $prog, $root) = @_;
	catfile($c->_root($root), 'blast', $c->_ver('blast', $root), 'bin', ($prog // 'blastp') . $EXE);
}

sub plink2 { my ($c, $root) = @_; catfile($c->_root($root), 'plink2', $c->_ver('plink2', $root), 'bin', "plink2$EXE") }

sub bcftools { my ($c, $root) = @_; catfile($c->_root($root), 'bcftools', $c->_ver('bcftools', $root), 'bin', "bcftools$EXE") }

sub prodigal { my ($c, $root) = @_; catfile($c->_root($root), 'prodigal', $c->_ver('prodigal', $root), 'bin', "prodigal$EXE") }

sub kallisto { my ($c, $root) = @_; catfile($c->_root($root), 'kallisto', $c->_ver('kallisto', $root), 'bin', "kallisto$EXE") }

# spades.py, which is run under python(): it is a script, and its "#!" line
# would find whatever python3 is first in PATH
sub spades { my ($c, $root) = @_; catfile($c->_root($root), 'spades', $c->_ver('spades', $root), 'bin', 'spades.py') }

sub python { my ($c, $root) = @_; catfile($c->_root($root), 'venv', $WIN ? ('Scripts', 'python.exe') : ('bin', 'python')) }

# an optional tool only where it was installed, which plink2, say, is not where
# PLINK 2 has no build
sub bin_dir {
	my $c = shift;
	my $v = $c->versions;
	map { dirname($_) } $c->clustalo, $c->blast, (map { $c->$_ } grep { $v->{$_} } sort keys %OPTIONAL), $c->python;
}

sub versions { my ($c, $root) = @_; my $s = $c->_state($root); +{ map { $_ => $s->{$_}{version} } keys %{ $s } } }

sub _ver {
	my ($c, $tool, $root) = @_;
	$c->_state($root)->{$tool}{version} // croak "$tool is not installed under " . $c->_root($root);
}

# tools.json in the share dir records which version of each tool is current. It
# is the only record: update() rewrites it, so the version Alien::Build wrote
# into its own runtime properties at install time goes stale after an update.
sub _state {
	my ($c, $root) = @_;
	my $f = catfile($c->_root($root), 'tools.json');
	return {} unless -f $f;
	require JSON::PP;
	open my $fh, '<', $f or croak "$f: $!";
	local $/;
	JSON::PP::decode_json(<$fh>);
}

sub _save {
	my ($c, $root, $tool, $version) = @_;
	my $s = $c->_state($root);
	$s->{$tool} = { version => $version };
	require JSON::PP;
	my $f = catfile($c->_root($root), 'tools.json');
	unlink $f; # ExtUtils::Install leaves it read-only; the directory is writable
	open my $fh, '>', $f or croak "$f: $!";
	print {$fh} JSON::PP->new->canonical->pretty->encode($s);
	close $fh or croak "$f: $!";
}

# ---- update checking --------------------------------------------------------

sub check_updates {
	my ($c, @want) = @_;
	@want = @TOOLS unless @want;
	my ($have, %r) = ($c->versions);
	foreach my $tool (@want) {
		if ($tool eq 'python-packages') {
			my $out = _run_out($c->python, qw(-m pip list --outdated --format=json --disable-pip-version-check));
			require JSON::PP;
			my %old = map { lc $_->{name} => $_ } @{ JSON::PP::decode_json($out) };
			# A package that was added to @PY_PACKAGES after this venv was made
			# is not outdated, but absent, and update() installs it too; pip
			# names what it lists as the package was installed, in any case.
			my %here = map { lc $_->{name} => 1 }
				@{ JSON::PP::decode_json(_run_out($c->python, qw(-m pip list --format=json --disable-pip-version-check))) };
			foreach my $p (@PY_PACKAGES) {
				$r{$tool}{$p} = { installed => undef, latest => undef, update => 1 } unless $here{lc $p};
				my $o = $old{lc $p} or next;
				$r{$tool}{$p} = { installed => $o->{version}, latest => $o->{latest_version}, update => 1 };
			}
			next;
		}
		my $l = eval { $c->latest($tool) };
		unless ($l) {
			# an optional tool with no build for this platform is not behind
			die $@ unless $OPTIONAL{$tool};
			(my $why = $@) =~ s/ at \S+ line \d+\.?\n?\z//;
			$r{$tool} = { installed => $have->{$tool}, latest => undef, update => 0, unavailable => $why };
			next;
		}
		my $latest = $l->{version};
		$r{$tool} = { installed => $have->{$tool}, latest => $latest,
			update => (!defined $have->{$tool} || _vcmp($latest, $have->{$tool}) > 0) ? 1 : 0 };
	}
	\%r;
}

sub latest {
	my ($c, $tool) = @_;
	return _blast_latest() if $tool eq 'blast';
	return _clustalo_latest() if $tool eq 'clustalo';
	return _plink2_latest() if $tool eq 'plink2';
	return _github_latest($tool) if $GITHUB{$tool};
	if ($tool eq 'python') {
		require JSON::PP;
		my ($v) = sort { _vcmp($b, $a) } map { $_->{name} =~ /^Python (3\.\d+\.\d+)$/ ? $1 : () } @{ JSON::PP::decode_json(_get($PY_URL)) };
		return { version => $v };
	}
	croak "latest() knows " . join(', ', grep { $_ ne 'python-packages' } @TOOLS) . ", not \"$tool\"";
}

# Install whatever check_updates() says is behind, then drop the version it
# replaced. Returns the names of the tools that changed.
sub update {
	my ($c, @want) = @_;
	my $root = $c->dist_dir;
	my $todo = $c->check_updates(@want);
	my @done;
	foreach my $tool (sort keys %{ $todo }) {
		if ($tool eq 'python-packages') {
			next unless %{ $todo->{$tool} };
			_run($c->python, qw(-m pip install --upgrade --disable-pip-version-check), _pip_index(), sort keys %{ $todo->{$tool} });
			_drop_scripts(catdir($root, 'venv'));
			push @done, $tool;
			next;
		}
		next unless $todo->{$tool}{update};
		if ($tool eq 'python') {
			# A venv cannot change interpreter, so a newer Python means a new
			# venv -- from whatever Alien::CPython3 now provides, which may be
			# no newer than before.
			my $base = _base_python_version();
			if (_vcmp($base, $todo->{$tool}{installed} // 0) <= 0) {
				warn "Python $todo->{$tool}{latest} is out, but Alien::CPython3 still provides $base; install a newer Python there first\n";
				next;
			}
			remove_tree(catdir($root, 'venv'));
		}
		my $old = $todo->{$tool}{installed};
		$c->install($tool, $root);
		remove_tree(catdir($root, $tool, $old)) if $tool ne 'python' && defined $old && -d catdir($root, $tool, $old);
		push @done, $tool;
	}
	@done;
}

# ---- installation -----------------------------------------------------------

sub install {
	my ($c, $tool, $root) = @_;
	make_path($root);
	if ($tool eq 'blast') {
		my $l = _blast_latest();
		my $bin = catdir($root, 'blast', $l->{version}, 'bin');
		remove_tree(catdir($root, 'blast', $l->{version})); # installed files are read-only, so copy() cannot overwrite them
		_local_copy('blastp', $l->{version}, $bin, 1) or _blast_unpack($l, $bin);
		$c->_save($root, blast => $l->{version});
	} elsif ($tool eq 'clustalo') {
		my $l = _clustalo_latest();
		my $bin = catdir($root, 'clustalo', $l->{version}, 'bin');
		remove_tree(catdir($root, 'clustalo', $l->{version}));
		_local_copy('clustalo', $l->{version}, $bin, 0) or ($l->{url} ? _clustalo_binary($l, $bin) : _clustalo_source($l, $bin));
		$c->_save($root, clustalo => $l->{version});
	} elsif ($tool eq 'plink2') {
		my $l = _plink2_latest();
		my $bin = catdir($root, 'plink2', $l->{version}, 'bin');
		remove_tree(catdir($root, 'plink2', $l->{version}));
		_local_copy('plink2', $l->{version}, $bin, 0) or _plink2_unpack($l, $bin);
		$c->_save($root, plink2 => $l->{version});
	} elsif ($tool eq 'bcftools' || $tool eq 'prodigal' || $tool eq 'kallisto') {
		my $l = _github_latest($tool);
		my $bin = catdir($root, $tool, $l->{version}, 'bin');
		remove_tree(catdir($root, $tool, $l->{version}));
		unless (_local_tool($tool, $l->{version}, $bin)) {
			$tool eq 'bcftools' ? _bcftools_source($l, $bin) : $tool eq 'prodigal' ? _prodigal_get($l, $bin) : _kallisto_unpack($l, $bin);
		}
		my $says = _says($tool, catfile($bin, "$tool$EXE")) // 'no version';
		$says eq $l->{version} or croak "the $tool installed in $bin reports $says, not $l->{version}";
		$c->_save($root, $tool => $l->{version});
	} elsif ($tool eq 'spades') {
		my $l = _github_latest($tool);
		my $to = catdir($root, 'spades', $l->{version});
		remove_tree($to);
		_spades_unpack($l, $to);
		# under its own Python, with -E and -s as Cn::Bioinf runs it, which keep
		# PYTHONPATH and the user's site-packages out. -I would do that too, but
		# also keeps the script's directory out of sys.path, and spades.py
		# imports spades_init.py from there.
		my $says = _says($tool, $c->python($root), '-E', '-s', catfile($to, 'bin', 'spades.py')) // 'no version';
		$says eq $l->{version} or croak "the spades.py installed in $to reports $says, not $l->{version}";
		$c->_save($root, spades => $l->{version});
	} elsif ($tool eq 'python') {
		my $venv = catdir($root, 'venv');
		_run(_base_python(), '-m', 'venv', $venv);
		_tidy_venv($venv) unless $WIN;
		my $py = $c->python($root);
		_run($py, qw(-m pip install --disable-pip-version-check), _pip_index(), @PY_PACKAGES);
		_drop_scripts($venv);
		_run($py, '-I', '-c', _py_imports());
		$c->_save($root, python => _line($py, '-c', 'import platform; print(platform.python_version())'));
	} else {
		croak 'install() knows ' . join(', ', sort grep { $_ ne 'python-packages' } @TOOLS) . ", not \"$tool\"";
	}
}

# A venv's bin/ is symlinks back to the base interpreter, and ExtUtils::Install
# copies what a link points at, so each of python, python3, python3.x (and on
# 3.14 a "U+1D70B thon" joke link) would become a copy of the interpreter, and
# lib64 -> lib a second copy of site-packages. Keep one link, straight to the
# base binary, and nothing else that is a link.
sub _tidy_venv {
	my ($venv) = @_;
	my $bin = catdir($venv, 'bin');
	opendir my $dh, $bin or croak "$bin: $!";
	foreach my $f (readdir $dh) {
		my $p = catfile($bin, $f);
		unlink $p if -l $p;
	}
	unlink catfile($venv, 'lib64') if -l catfile($venv, 'lib64');
	my $base = _line(_base_python(), '-c', 'import sys; print(sys.executable)');
	symlink($base, catfile($bin, 'python')) or croak "symlink $base: $!";
}

# Every file in the venv's bin/ (Scripts\ on Windows) but the interpreter: pip,
# activate and the console scripts of the packages. The venv is built in the
# staging directory and copied into place, so each of them names the staging
# copy, which is gone once the build is cleaned; python -m pip is what this
# module runs.
sub _drop_scripts {
	my ($venv) = @_;
	my $bin = catdir($venv, $WIN ? 'Scripts' : 'bin');
	opendir my $dh, $bin or croak "$bin: $!";
	foreach my $f (readdir $dh) {
		my $p = catfile($bin, $f);
		next if !-f $p || ($WIN ? $f =~ /\Apythonw?\.exe\z/i : $f eq 'python');
		unlink $p or croak "can't delete $p: $!";
	}
	closedir $dh;
}

# Python code that imports every one of @PY_PACKAGES, which install() runs to
# check the venv. Up to 0.01 it imported each by its PyPI name, and
# "import biopython" fails: Biopython is imported as Bio.
sub _py_imports { 'import ' . join(',', map { $PY_IMPORT{$_} // $_ } @PY_PACKAGES) }

# --no-index when ALIEN_INSTALL_NETWORK forbids downloads, so that pip installs
# only from what is already here, or fails saying what it can't find.
sub _pip_index { _network_ok() ? () : '--no-index' }

# The interpreter itself, not a pyenv shim or the like: a venv made through a
# shim records the shim, which picks a Python by the working directory.
sub _base_python {
	require Alien::CPython3;
	my $exe = Alien::CPython3->exe;
	my ($dir) = Alien::CPython3->bin_dir;
	$exe = catfile($dir, $exe) if defined $dir;
	_line($exe, '-c', 'import sys; print(sys.executable)');
}

sub _base_python_version { _line(_base_python(), '-c', 'import platform; print(platform.python_version())') }

# A copy of $prog already on this machine at exactly $version, copied into
# $bin instead of downloading. ALIEN_BIOINF_CLUSTALO / ALIEN_BIOINF_BLAST /
# ALIEN_BIOINF_PLINK2 name one explicitly (the clustalo or plink2 binary itself,
# or BLAST's bin directory); failing that, PATH is searched. BLAST is a
# directory of ~30 programs, so a PATH hit is only copied from a directory that
# is plainly an unpacked NCBI tarball (.../ncbi-blast-<v>+/bin) -- never, say,
# all of /usr/bin. A plink2's version is its build date, as _plink2_date reads it.
my %LOCAL_ENV = (blastp => 'ALIEN_BIOINF_BLAST', clustalo => 'ALIEN_BIOINF_CLUSTALO', plink2 => 'ALIEN_BIOINF_PLINK2');
sub _local_copy {
	my ($prog, $version, $bin, $whole_dir) = @_;
	my $env = $ENV{ $LOCAL_ENV{$prog} };
	my @cand = defined $env ? ($whole_dir ? catfile($env, "$prog$EXE") : $env) : ();
	push @cand, map { catfile($_, "$prog$EXE") } File::Spec->path unless @cand;
	foreach my $exe (grep { -f $_ && -x _ } @cand) {
		next if $whole_dir && !defined $env && dirname($exe) !~ /ncbi-blast-\Q$version\E\+[\/\\]bin\z/;
		my $out = eval { _run_out($exe, $prog eq 'blastp' ? '-version' : '--version') } // next;
		next unless $prog eq 'plink2' ? (_plink2_date($out) // '') eq $version : $out =~ /^\s*(?:blastp:\s*)?\Q$version\E\+?\s*$/m;
		make_path($bin);
		my @files = $whole_dir ? _dir_files(dirname($exe)) : ($exe);
		foreach my $f (@files) {
			my $to = catfile($bin, $whole_dir ? basename($f) : "$prog$EXE");
			copy($f, $to) or croak "copy $f: $!";
			chmod 0755, $to;
		}
		print "copied $prog $version from ", dirname($exe), " instead of downloading it\n";
		return 1;
	}
	0;
}

sub _dir_files { my ($d) = @_; opendir my $dh, $d or croak "$d: $!"; grep { -f $_ } map { catfile($d, $_) } readdir $dh }

# NCBI's name for this machine's build, e.g. "x64-linux".
sub _blast_platform {
	my $cpu = _cpu();
	my %p = (linux => { x64 => 'x64-linux', aarch64 => 'aarch64-linux' },
		darwin => { x64 => 'x64-macosx', aarch64 => 'aarch64-macosx' },
		MSWin32 => { x64 => 'x64-win64' });
	$p{$^O}{$cpu} // croak "NCBI publishes no BLAST+ build for $^O/$cpu; set ALIEN_BIOINF_BLAST to the bin directory of one built here";
}

sub _blast_latest {
	my $plat = _blast_platform();
	my $page = _get($BLAST_URL);
	$page =~ /href="(ncbi-blast-(\d[\d.]*)\+-\Q$plat\E\.tar\.gz)"/ or croak "no $plat tarball listed at $BLAST_URL";
	{ version => $2, url => "$BLAST_URL$1", file => $1 };
}

sub _blast_unpack {
	my ($l, $bin) = @_;
	my $tgz = _fetch($l->{url}, $l->{file});
	my ($md5) = _get("$l->{url}.md5") =~ /^([0-9a-f]{32})/ or croak "no md5 at $l->{url}.md5";
	require Digest::MD5;
	open my $fh, '<:raw', $tgz or croak "$tgz: $!";
	my $got = Digest::MD5->new->addfile($fh)->hexdigest;
	if ($got ne $md5) {
		unlink $tgz;
		croak "$l->{file}: md5 $got, but NCBI says $md5; the cached copy has been deleted, so a re-run downloads it again";
	}
	make_path($bin);
	# iter() streams the archive; read() would hold all of its ~250 MB of
	# programs in memory at once. Only bin/ is kept: doc/ and the licence
	# text are not needed to run anything.
	require Archive::Tar;
	my $next = Archive::Tar->iter($tgz, 1) or croak "can't read $tgz";
	while (my $f = $next->()) {
		next unless $f->is_file && $f->full_path =~ m{/bin/([^/]+)\z};
		$f->extract(catfile($bin, $1)) or croak "can't extract $1 from $tgz";
		chmod 0755, catfile($bin, $1);
	}
	-f catfile($bin, "blastp$EXE") or croak "$tgz had no bin/blastp";
}

# clustal.org's binary for this platform: [regex, file extension or '' for a
# bare executable]. Anything missing here is built from source.
sub _clustalo_pattern {
	my $cpu = _cpu();
	return ['clustalo-(\d[\d.]*)-Ubuntu-x86_64', ''] if $^O eq 'linux' && $cpu eq 'x64';
	return ['clustalo-(\d[\d.]*)-Ubuntu-32-bit', ''] if $^O eq 'linux' && $cpu eq 'x86';
	return ['clustal-omega-(\d[\d.]*)-macosx', ''] if $^O eq 'darwin';
	return ['clustal-omega-(\d[\d.]*)-win64\.zip', 'zip'] if $WIN;
	undef;
}

sub _clustalo_latest {
	my $pat = _clustalo_pattern();
	if ($pat) {
		my $page = eval { _get($CLUSTAL_URL) } // '';
		my ($best, $file);
		while ($page =~ /href="(?:[^"]*\/)?($pat->[0])"/g) {
			($best, $file) = ($2, $1) if !defined $best || _vcmp($2, $best) > 0;
		}
		return { version => $best, url => "$CLUSTAL_URL$file", file => $file, zip => $pat->[1] } if defined $best;
		croak "$CLUSTAL_URL could not be read and there is no Windows source build; set ALIEN_BIOINF_CLUSTALO to a clustalo.exe" if $WIN;
	}
	require JSON::PP;
	my ($v) = sort { _vcmp($b, $a) } grep { /^\d+(?:\.\d+)*\z/ } map { $_->{name} } @{ JSON::PP::decode_json(_get($CLUSTAL_GIT)) };
	defined $v or croak "no release tags at $CLUSTAL_GIT";
	{ version => $v, src => "https://github.com/GSLBiotech/clustal-omega/archive/refs/tags/$v.tar.gz" };
}

sub _clustalo_binary {
	my ($l, $bin) = @_;
	my $got = _fetch($l->{url}, $l->{file});
	make_path($bin);
	my $to = catfile($bin, "clustalo$EXE");
	if ($l->{zip}) {
		# tar is bsdtar on Windows 10 and later, which reads zip; the zip's
		# clustalo.exe needs the DLLs beside it, so everything is kept.
		my $tmp = _tmpdir();
		_run('tar', '-xf', $got, '-C', $tmp);
		require File::Find;
		File::Find::find(sub { copy($_, catfile($bin, $_)) if -f $_ && ($_ =~ /\.dll\z/i || lc $_ eq 'clustalo.exe') }, $tmp);
	} else {
		copy($got, $to) or croak "copy $got: $!";
	}
	chmod 0755, $to;
	_run_out($to, '--version') =~ /\Q$l->{version}\E/ or croak "$to does not report version $l->{version}";
}

sub _clustalo_source {
	my ($l, $bin) = @_;
	print "no clustal.org binary for this platform could be had; building clustalo $l->{version} from source\n";
	my $tmp = _tmpdir();
	my $at = catdir($tmp, 'argtable');
	require Archive::Tar;
	require Cwd;
	my $here = Cwd::getcwd();
	my $make = $Config::Config{gmake} || $Config::Config{make} || 'make';
	# argtable2-13's config.guess dates from 2010 and cannot name an aarch64
	# machine -- the very one with no clustal.org binary -- so tell configure
	# what the compiler targets instead of letting it guess.
	my $triplet = eval { _line($Config::Config{cc}, '-dumpmachine') };
	my @build = $triplet ? ("--build=$triplet") : ();
	eval {
		chdir $tmp or croak "$tmp: $!";
		Archive::Tar->extract_archive(_fetch($ARGTABLE_URL, 'argtable2-13.tar.gz'), 1) or croak Archive::Tar->error;
		Archive::Tar->extract_archive(_fetch($l->{src}, "clustal-omega-$l->{version}.tar.gz"), 1) or croak Archive::Tar->error;
		# Static argtable2, so the clustalo that comes out needs nothing that
		# is not already on any machine with a C++ compiler.
		chdir catdir($tmp, 'argtable2-13') or croak "argtable2-13: $!";
		_run('sh', './configure', @build, "--prefix=$at", '--disable-shared');
		_run($make);
		_run($make, 'install');
		chdir catdir($tmp, "clustal-omega-$l->{version}") or croak "clustal-omega-$l->{version}: $!";
		_run('sh', './configure', @build, "--prefix=" . catdir($tmp, 'co'), '--disable-shared',
			"CPPFLAGS=-I" . catdir($at, 'include'), "LDFLAGS=-L" . catdir($at, 'lib'));
		_run($make);
		_run($make, 'install');
		1;
	} or do { my $e = $@; chdir $here; die $e };
	chdir $here;
	make_path($bin);
	copy(catfile($tmp, 'co', 'bin', 'clustalo'), catfile($bin, 'clustalo')) or croak "copy clustalo: $!";
	chmod 0755, catfile($bin, 'clustalo');
}

# The platform in the name of PLINK 2's zip for this machine, or undef where it
# publishes none.
sub _plink2_platform {
	my $cpu = _cpu();
	my %p = (linux => { x64 => 'linux_x86_64', x86 => 'linux_i686' },
		darwin => { x64 => 'mac', aarch64 => 'mac_arm64' },
		MSWin32 => { x64 => 'win64', x86 => 'win32' });
	$p{$^O}{$cpu};
}

my %MONTH = (Jan => 1, Feb => 2, Mar => 3, Apr => 4, May => 5, Jun => 6, Jul => 7, Aug => 8, Sep => 9, Oct => 10, Nov => 11, Dec => 12);

# The build date in plink2 --version's "PLINK v2.0.0-a.7.13LM 64-bit Intel (15
# Oct 2026)", as 20261015, which is how its zips are named; undef if there is none.
sub _plink2_date {
	my ($out) = @_;
	$out =~ /^PLINK v2\S* .*\((\d{1,2}) ([A-Z][a-z]{2}) (\d{4})\)/m && $MONTH{$2} or return undef;
	sprintf '%04d%02d%02d', $3, $MONTH{$2}, $1;
}

# The newest zip in the stable column of $PLINK2_URL for this platform. Where
# there is none, a plink2 named by ALIEN_BIOINF_PLINK2 is used at whatever
# build it is.
sub _plink2_latest {
	my $plat = _plink2_platform();
	unless (defined $plat) {
		my $exe = $ENV{ALIEN_BIOINF_PLINK2}
			// croak "PLINK 2 publishes no build for $^O/" . _cpu() . "; set ALIEN_BIOINF_PLINK2 to a plink2 built here";
		my $date = _plink2_date(_run_out($exe, '--version')) // croak "$exe --version gives no build date";
		return { version => $date };
	}
	my $page = _get($PLINK2_URL);
	my ($best, $url);
	while ($page =~ m{href="(https://s3\.amazonaws\.com/plink2-assets/alpha\d+/plink2_\Q$plat\E_(\d{8})\.zip)"}g) {
		($best, $url) = ($2, $1) if !defined $best || $2 > $best;
	}
	defined $best or croak "no plink2_${plat}_<date>.zip in the stable column of $PLINK2_URL";
	{ version => $best, url => $url, file => basename($url) };
}

sub _plink2_unpack {
	my ($l, $bin) = @_;
	my $zip = _fetch($l->{url}, $l->{file});
	require IO::Uncompress::Unzip;
	my $z = IO::Uncompress::Unzip->new($zip) or croak "can't read $zip: $IO::Uncompress::Unzip::UnzipError";
	my $to = catfile($bin, "plink2$EXE");
	make_path($bin);
	my $status = 1;
	while ($status > 0) {
		if ($z->getHeaderInfo->{Name} eq "plink2$EXE") {
			open my $fh, '>:raw', $to or croak "$to: $!";
			my $buf;
			while (($status = $z->read($buf)) > 0) {
				print {$fh} $buf;
			}
			croak "can't unzip plink2$EXE from $zip: $IO::Uncompress::Unzip::UnzipError" if $status < 0;
			close $fh or croak "$to: $!";
		}
		$status = $z->nextStream;
	}
	croak "can't read $zip: $IO::Uncompress::Unzip::UnzipError" if $status < 0;
	-s $to or croak "$zip had no plink2$EXE";
	chmod 0755, $to;
	my $date = _plink2_date(_run_out($to, '--version')) // '';
	$date eq $l->{version} or croak "$to reports build date \"$date\", not $l->{version}, as $l->{file} is named";
}

# ---- the tools released on GitHub ---------------------------------------------

# The version that the program @cmd, $tool or the interpreter and script that
# run it, reports, or undef if it does not run or says something else. STDERR
# is read too, which is where Prodigal writes its version.
sub _says {
	my ($tool, @cmd) = @_;
	require Capture::Tiny;
	my $out = eval { Capture::Tiny::capture_merged(sub { system @cmd, @{ $GITHUB{$tool}{args} } }) } // return undef;
	$out =~ $GITHUB{$tool}{says} ? $1 : undef;
}

# The newest release of $repo on GitHub, as its version and its tag, which is
# the version with a "v" before it for all but bcftools: the tag that
# /releases/latest redirects to. HTTP::Tiny follows the redirect and gives the
# last URL it asked for.
sub _github_tag {
	my ($repo) = @_;
	my $url = "https://github.com/$repo/releases/latest";
	my $r = _http()->request('HEAD', $url);
	croak "HEAD $url: $r->{status} $r->{reason}" unless $r->{success};
	my @release = _github_tag_of($r->{url} // '');
	@release or croak "$url leads to \"" . ($r->{url} // '') . '", which names no release';
	@release;
}

# The version and the tag in the URL of a release's page on GitHub, or an empty
# list for any other URL.
sub _github_tag_of { $_[0] =~ m{/releases/tag/(v?(\d[\w.]*))\z} ? ($2, $1) : () }

# The file of each tool's release for this machine, as the part of its name
# that says so; a machine left out has none. Prodigal's are whole file names.
my %ASSET = (
	prodigal => { linux => { x64 => 'prodigal.linux' }, darwin => { x64 => 'prodigal.osx.10.9.5' }, MSWin32 => { x64 => 'prodigal.windows.exe' } },
	kallisto => { linux => { x64 => 'linux' }, darwin => { aarch64 => 'mac' } },
	spades => { linux => { x64 => 'Linux' }, darwin => { x64 => 'Darwin-x86_64', aarch64 => 'Darwin-arm64' } },
);
my %ENV_OF = (bcftools => 'ALIEN_BIOINF_BCFTOOLS', prodigal => 'ALIEN_BIOINF_PRODIGAL', kallisto => 'ALIEN_BIOINF_KALLISTO');

# The newest release of $tool and where its file for this machine is: {version,
# url, file}, "url" being left out for one that is only built from "src".
# Where there is neither a file nor a way to build it, the program that
# ALIEN_BIOINF_<TOOL> names is used at whatever version it reports.
sub _github_latest {
	my ($tool) = @_;
	my $repo = $GITHUB{$tool}{repo};
	my $asset = $ASSET{$tool}{$^O}{ _cpu() };
	my $buildable = $tool eq 'bcftools' ? !$WIN : $tool eq 'prodigal' ? !$WIN : 0; # 1 = it can be built from source here
	unless (defined $asset || $buildable) {
		my $exe = defined $ENV_OF{$tool} ? $ENV{ $ENV_OF{$tool} } : undef;
		croak "$repo publishes no build of $tool for $^O/" . _cpu()
			. (defined $ENV_OF{$tool} ? "; set $ENV_OF{$tool} to a $tool built here" : '') unless defined $exe;
		my $v = _says($tool, $exe) // croak "$exe does not report a $tool version";
		return { version => $v };
	}
	my ($v, $tag) = _github_tag($repo);
	my $base = "https://github.com/$repo/releases/download/$tag";
	return { version => $v, url => "$base/bcftools-$v.tar.bz2", file => "bcftools-$v.tar.bz2" } if $tool eq 'bcftools';
	if ($tool eq 'prodigal') {
		# the cached file carries the version, which the release's own name lacks
		return { version => $v, src => "https://github.com/$repo/archive/refs/tags/$tag.tar.gz", srcfile => "Prodigal-$v.tar.gz",
			(defined $asset ? (url => "$base/$asset", file => "prodigal-$v." . ($asset =~ /\Aprodigal\.(.+)\z/)[0]) : ()) };
	}
	my $file = $tool eq 'kallisto' ? "kallisto_$asset-$tag.tar.gz" : "SPAdes-$v-$asset.tar.gz";
	{ version => $v, url => "$base/$file", file => $file };
}

# A $tool already on this machine at exactly $version, copied into $bin instead
# of being downloaded or built: the one ALIEN_BIOINF_<TOOL> names, or else the
# first in PATH. As _local_copy, for the tools that are one program.
sub _local_tool {
	my ($tool, $version, $bin) = @_;
	my $env = $ENV{ $ENV_OF{$tool} };
	my @cand = defined $env ? ($env) : map { catfile($_, "$tool$EXE") } File::Spec->path;
	foreach my $exe (grep { -f $_ && -x _ } @cand) {
		next unless (_says($tool, $exe) // '') eq $version;
		make_path($bin);
		my $to = catfile($bin, "$tool$EXE");
		copy($exe, $to) or croak "copy $exe: $!";
		chmod 0755, $to;
		print "copied $tool $version from ", dirname($exe), " instead of ", ($tool eq 'bcftools' ? 'building' : 'downloading'), " it\n";
		return 1;
	}
	0;
}

# The archive $tgz unpacked into $dir, which Archive::Tar does only from the
# working directory.
sub _untar {
	my ($tgz, $dir) = @_;
	require Archive::Tar;
	require Cwd;
	croak "$tgz is bzip2-compressed, which Archive::Tar reads only with IO::Uncompress::Bunzip2"
		if $tgz =~ /\.bz2\z/ && !Archive::Tar->has_bzip2_support;
	my $here = Cwd::getcwd();
	chdir $dir or croak "$dir: $!";
	my $ok = Archive::Tar->extract_archive($tgz);
	my $why = Archive::Tar->error;
	chdir $here or croak "$here: $!";
	$ok or croak "can't unpack $tgz: $why";
}

# @cmd run in $dir.
sub _run_in {
	my ($dir, @cmd) = @_;
	require Cwd;
	my $here = Cwd::getcwd();
	chdir $dir or croak "$dir: $!";
	my $ok = eval { _run(@cmd); 1 };
	my $why = $@;
	chdir $here or croak "$here: $!";
	$ok or die $why;
}

# bcftools built from its release tarball, which carries the HTSlib it links
# statically, so the one program is all there is to install. libcurl is left
# out: it is for reading files over the network, which Cn::Bioinf does not ask
# of it, and it would tie the program to the libcurl of the day it was built.
# HTSlib's configure stops where bzip2's or liblzma's headers are missing, and
# is then run again without the two, which costs only the CRAM files whose
# blocks are compressed with them.
sub _bcftools_source {
	my ($l, $bin) = @_;
	croak 'bcftools has no Windows build, and is not built here; set ALIEN_BIOINF_BCFTOOLS to a bcftools.exe' if $WIN;
	my $tmp = _tmpdir();
	_untar(_fetch($l->{url}, $l->{file}), $tmp);
	my $src = catdir($tmp, "bcftools-$l->{version}");
	my $make = $Config::Config{gmake} || $Config::Config{make} || 'make';
	print "building bcftools $l->{version} from source\n";
	unless (eval { _run_in($src, 'sh', './configure', '--disable-libcurl'); 1 }) {
		print "configuring bcftools again without bzip2 and liblzma, which CRAM files may need: $@";
		_run_in($src, 'sh', './configure', '--disable-libcurl', '--disable-bz2', '--disable-lzma');
	}
	# the "bcftools" target alone: the default one builds its 40 plugins and its
	# test programs as well, none of which is installed
	_run_in($src, $make, 'bcftools');
	make_path($bin);
	copy(catfile($src, 'bcftools'), catfile($bin, 'bcftools')) or croak "copy bcftools: $!";
	chmod 0755, catfile($bin, 'bcftools');
}

# Prodigal's program for this machine, or, where it has none or the release
# does not hold the file by the name v2.6.3 gave it, one built from source.
sub _prodigal_get {
	my ($l, $bin) = @_;
	my $to = catfile($bin, "prodigal$EXE");
	make_path($bin);
	if ($l->{url}) {
		my $got = eval { _fetch($l->{url}, $l->{file}) };
		if (defined $got) {
			copy($got, $to) or croak "copy $got: $!";
			chmod 0755, $to;
			return;
		}
		(my $why = $@) =~ s/\s+\z//;
		croak "$why; Prodigal is not built from source on Windows, so set ALIEN_BIOINF_PRODIGAL to a prodigal.exe" if $WIN;
		print "no Prodigal binary for this platform could be had ($why); building it from source\n";
	}
	my $tmp = _tmpdir();
	_untar(_fetch($l->{src}, $l->{srcfile}), $tmp);
	# its Makefile names gcc, and uses GNU make's $(shell)
	_run_in(catdir($tmp, "Prodigal-$l->{version}"), $Config::Config{gmake} || $Config::Config{make} || 'make', "CC=$Config::Config{cc}");
	copy(catfile($tmp, "Prodigal-$l->{version}", 'prodigal'), $to) or croak "copy prodigal: $!";
	chmod 0755, $to;
}

# kallisto's program, from the release's tarball; its README, licence and test
# reads are left there.
sub _kallisto_unpack {
	my ($l, $bin) = @_;
	my $tgz = _fetch($l->{url}, $l->{file});
	require Archive::Tar;
	my $next = Archive::Tar->iter($tgz, 1) or croak "can't read $tgz";
	make_path($bin);
	my $to = catfile($bin, "kallisto$EXE");
	while (my $f = $next->()) {
		next unless $f->is_file && $f->full_path =~ m{\Akallisto/kallisto\z};
		$f->extract($to) or croak "can't extract kallisto from $tgz";
	}
	-s $to or croak "$tgz had no kallisto/kallisto";
	chmod 0755, $to;
}

# The part of SPAdes's tarball that spades.py runs to assemble reads, unpacked
# into $to; the comment at the top of this file says what is left out and how
# much room that saves. iter() streams the archive, which read() would hold
# whole in memory. The links from metaspades.py and the like to spades.py are
# left out too: spades.py takes --meta for the same thing, and
# ExtUtils::Install would install each as another copy.
sub _spades_unpack {
	my ($l, $to) = @_;
	my $tgz = _fetch($l->{url}, $l->{file});
	require Archive::Tar;
	my $next = Archive::Tar->iter($tgz, 1) or croak "can't read $tgz";
	my $programs = qr{\Abin/(?:spades\.py|spades_init\.py|spades-core|spades-hammer|spades-ionhammer|spades-corrector-core|spades-bwa)\z};
	while (my $f = $next->()) {
		next unless $f->is_file;
		(my $path = $f->full_path) =~ s{\A[^/]+/}{}; # less "SPAdes-4.3.0-Linux/"
		next unless $path =~ $programs || ($path =~ m{\Ashare/spades/} && $path !~ m{\Ashare/spades/(?:sewage|coronaspades_hmms)/});
		my $file = catfile($to, split m{/}, $path);
		make_path(dirname($file));
		$f->extract($file) or croak "can't extract $path from $tgz";
		chmod 0755, $file if $path =~ m{\Abin/};
	}
	-s catfile($to, 'bin', 'spades.py') && -s catfile($to, 'bin', 'spades-core') or croak "$tgz had no bin/spades.py and bin/spades-core";
}

# ---- helpers ----------------------------------------------------------------

# x64, aarch64, x86 or other. $Config{archname} says nothing about the CPU on
# darwin ("darwin-2level"), so ask uname, and on Windows the environment.
sub _cpu {
	my $m = $WIN ? ($ENV{PROCESSOR_ARCHITEW6432} // $ENV{PROCESSOR_ARCHITECTURE} // '') : do { require POSIX; (POSIX::uname())[4] };
	return 'x64' if $m =~ /^(?:x86_64|amd64|x64)\z/i;
	return 'aarch64' if $m =~ /^(?:aarch64|arm64)\z/i;
	return 'x86' if $m =~ /^(?:i[3-6]86|x86)\z/i;
	'other';
}

sub cache_dir {
	my $base = $ENV{ALIEN_BIOINF_CACHE} // catdir($ENV{XDG_CACHE_HOME} // ($WIN ? $ENV{LOCALAPPDATA} : undef) // catdir($ENV{HOME} // $ENV{USERPROFILE}, '.cache'), 'alien-bioinf');
	make_path($base);
	$base;
}

sub _tmpdir { require File::Temp; File::Temp::tempdir(CLEANUP => 1) }

# False when ALIEN_INSTALL_NETWORK is set to a false value, which is how
# Alien::Build is told not to download anything.
sub _network_ok { !defined $ENV{ALIEN_INSTALL_NETWORK} || $ENV{ALIEN_INSTALL_NETWORK} ? 1 : 0 }

sub _http {
	croak 'ALIEN_INSTALL_NETWORK is ' . ($ENV{ALIEN_INSTALL_NETWORK} eq '' ? 'empty' : $ENV{ALIEN_INSTALL_NETWORK}) . ', so nothing may be downloaded'
		unless _network_ok();
	require HTTP::Tiny;
	HTTP::Tiny->new(agent => "Alien-Bioinf/$VERSION ", timeout => 120, verify_SSL => 1);
}

sub _get {
	my ($url) = @_;
	my $r = _http()->get($url);
	croak "GET $url: $r->{status} $r->{reason}" unless $r->{success};
	$r->{content};
}

# $url saved as $name in cache_dir(), unless it is there already. The download
# goes to a .part file that is renamed only once complete, so an interrupted
# one is never mistaken for a cached archive.
sub _fetch {
	my ($url, $name) = @_;
	my $f = catfile(cache_dir(), $name);
	if (-s $f) {
		print "using cached $f\n";
		return $f;
	}
	print "downloading $url\n";
	open my $fh, '>:raw', "$f.part" or croak "$f.part: $!";
	my $r = _http()->get($url, { data_callback => sub { print {$fh} $_[0] } });
	close $fh or croak "$f.part: $!";
	unless ($r->{success}) {
		unlink "$f.part";
		croak "GET $url: $r->{status} $r->{reason}";
	}
	rename "$f.part", $f or croak "rename $f.part: $!";
	$f;
}

sub _run { system(@_) == 0 or croak "\"@_\" failed: " . ($? == -1 ? $! : 'exit ' . ($? >> 8)) }

sub _run_out {
	my @cmd = @_;
	require Capture::Tiny;
	my ($out, $err, $rc) = Capture::Tiny::capture(sub { system @cmd });
	croak "\"@cmd\" failed: $err" if $rc != 0;
	$out;
}

sub _line { my $s = _run_out(@_); $s =~ s/\s+\z//; $s }

# Dotted-number comparison: 2.17.0 > 2.9.1, which "gt" gets wrong.
sub _vcmp {
	my @a = split /\./, $_[0];
	my @b = split /\./, $_[1];
	while (@a || @b) {
		my $c = (shift(@a) // 0) <=> (shift(@b) // 0);
		return $c if $c;
	}
	0;
}

1;
__END__

=head1 NAME

Alien::Bioinf - module-local Clustal Omega, BLAST+, PLINK 2, bcftools, Prodigal, kallisto, SPAdes and Python for Cn::Bioinf

=head1 SYNOPSIS

 use Alien::Bioinf;
 system Alien::Bioinf->clustalo, '--in', 'seqs.fa', '--out', 'aln.fa';
 system Alien::Bioinf->blast('blastp'), '-query', 'q.fa', '-subject', 's.fa';
 system Alien::Bioinf->plink2, '--vcf', 'calls.vcf.gz', '--make-pgen', '--out', 'calls';
 system Alien::Bioinf->bcftools, 'view', '-h', 'calls.vcf.gz';
 system Alien::Bioinf->prodigal, '-i', 'contigs.fa', '-f', 'gff', '-o', 'genes.gff';
 system Alien::Bioinf->kallisto, 'index', '-i', 'tx.idx', 'transcripts.fa';
 system Alien::Bioinf->python, '-E', '-s', Alien::Bioinf->spades, '-1', 'r1.fq.gz', '-2', 'r2.fq.gz', '-o', 'assembly';
 system Alien::Bioinf->python, 'script.py';

 my $behind = Alien::Bioinf->check_updates;          # every tool
 Alien::Bioinf->update('blast', 'python-packages');  # just these

From the shell:

 perl -MAlien::Bioinf -MData::Dumper -e 'print Dumper(Alien::Bioinf->check_updates)'
 perl -MAlien::Bioinf -e 'print "$_\n" for Alien::Bioinf->update'

=head1 DESCRIPTION

Installing this distribution puts the newest Clustal Omega, BLAST+, PLINK 2,
bcftools, Prodigal, kallisto and SPAdes for the running operating system, and
a Python venv holding numpy, matplotlib, adjustText, Biopython and PyDESeq2,
into its own share directory. None of it touches PATH.

Clustal Omega, BLAST+ and the Python are always installed. The other five are
installed where they can be, and left out, with a warning, where they cannot:

=over

=item PLINK 2

has no build for Linux on aarch64, unless C<ALIEN_BIOINF_PLINK2> names a plink2
built there.

=item bcftools

is built from source, with the HTSlib its tarball carries, which needs a C
compiler, make and zlib's headers; with bzip2's and liblzma's headers too it
reads every CRAM file, and without them all but those compressed with either.
It is not built on Windows. C<ALIEN_BIOINF_BCFTOOLS> may name a bcftools.

=item Prodigal

has binaries for Linux and macOS on x86_64 and for Windows, and is built from
source elsewhere but on Windows. C<ALIEN_BIOINF_PRODIGAL> may name a prodigal.

=item kallisto

has binaries for Linux on x86_64 and macOS on arm64, and is not built from
source. C<ALIEN_BIOINF_KALLISTO> may name a kallisto.

=item SPAdes

has binaries for Linux on x86_64 and for macOS. Only what C<spades.py> runs to
assemble reads is kept (185 MB of the 684 MB that 4.3.0 unpacks to), and
C<spades.py> is run under this distribution's Python.

=back

A tool added by a newer version of this module is installed by C<update>, as
one that is behind is.

=head2 Nothing is downloaded twice

Downloads are kept in L</cache_dir> and reused from there by later installs
and updates. Before downloading clustalo, BLAST, PLINK 2, bcftools, Prodigal or
kallisto, the installer looks for a copy already on the machine at the version
it wants and copies that instead: C<ALIEN_BIOINF_CLUSTALO>,
C<ALIEN_BIOINF_PLINK2>, C<ALIEN_BIOINF_BCFTOOLS>, C<ALIEN_BIOINF_PRODIGAL> and
C<ALIEN_BIOINF_KALLISTO> may each name the program, and C<ALIEN_BIOINF_BLAST>
the C<bin> directory of an unpacked BLAST+, and otherwise C<PATH> is searched.
PLINK 2's version is its build date, such as C<20261015>.
pip reuses the wheels in its own cache. The venv does not see the base
interpreter's packages, so a change to those cannot break it.

Setting C<ALIEN_INSTALL_NETWORK> to 0 stops every download, pip's included.

=head1 METHODS

=over

=item clustalo, blast($program), plink2, bcftools, prodigal, kallisto, python

Full paths to the installed programs. C<blast> defaults to C<blastp>.
C<plink2>, C<bcftools>, C<prodigal> and C<kallisto> die where the program was
not installed.

=item spades

The full path to C<spades.py>, to be run under C<python>, with C<-E -s> to
keep C<PYTHONPATH> and the user's own packages out (C<-I> would also hide
C<spades_init.py>, which sits beside it). It dies where SPAdes was not
installed.

=item versions

A hash reference of tool => installed version.

=item check_updates(@tools)

For each of C<clustalo>, C<blast>, C<plink2>, C<bcftools>, C<prodigal>,
C<kallisto>, C<spades>, C<python> and C<python-packages> (all nine by
default), what is installed and what the upstream site offers. Each tool maps
to C<< { installed, latest, update } >>; one that has no build for this
platform has C<latest> undef, C<update> 0 and the reason in C<unavailable>.
C<python-packages> maps each outdated package to the same, and each package
that is missing to C<< { installed => undef, latest => undef, update => 1 } >>.

=item update(@tools)

Installs whatever C<check_updates> reports as behind or missing and returns
the names of the tools that changed. A newer Python needs a newer interpreter
from L<Alien::CPython3> first, since a venv cannot change interpreter.

=item latest($tool)

The newest upstream version of any tool but C<python-packages>.

=item cache_dir

C<ALIEN_BIOINF_CACHE>, else C<alien-bioinf> under C<XDG_CACHE_HOME>, else
under C<~/.cache> (C<%LOCALAPPDATA%> on Windows).

=back

=cut
