package Alien::Bioinf;
# Module-local Clustal Omega, NCBI BLAST+ and a Python venv for Bioinf::Basic.
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
#   Python     Alien::CPython3 supplies the interpreter; a venv on top of it
#              gets numpy, matplotlib and adjustText -- what share/msa_plot.py
#              and Matplotlib::Simple import -- from PyPI. The newest release
#              for check_updates() is the highest "Python 3.x.y" in
#              https://www.python.org/api/v2/downloads/release/.
#
# Nothing is downloaded twice. Every archive is kept in cache_dir() and reused
# from there; a clustalo or BLAST already on this machine at the wanted version
# is copied instead of fetched (see _local_copy); and the venv is made with
# --system-site-packages, so pip installs only what the base Python lacks.
use strict;
use warnings;
use parent 'Alien::Base';
use Carp qw(croak);
use Config ();
use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path remove_tree);
use File::Copy qw(copy);
use File::Basename qw(basename dirname);
our $VERSION = '0.01';

my @PY_PACKAGES = qw(numpy matplotlib adjustText);
my @TOOLS = qw(clustalo blast python python-packages);
my $BLAST_URL = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/';
my $CLUSTAL_URL = 'https://www.clustal.org/omega/';
my $CLUSTAL_GIT = 'https://api.github.com/repos/GSLBiotech/clustal-omega/tags';
my $ARGTABLE_URL = 'https://downloads.sourceforge.net/project/argtable/argtable/argtable-2.13/argtable2-13.tar.gz';
my $PY_URL = 'https://www.python.org/api/v2/downloads/release/?is_published=true&pre_release=false';
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

sub python { my ($c, $root) = @_; catfile($c->_root($root), 'venv', $WIN ? ('Scripts', 'python.exe') : ('bin', 'python')) }

sub bin_dir { my $c = shift; map { dirname($_) } $c->clustalo, $c->blast, $c->python }

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
			foreach my $p (@PY_PACKAGES) {
				my $o = $old{lc $p} or next;
				$r{$tool}{$p} = { installed => $o->{version}, latest => $o->{latest_version}, update => 1 };
			}
			next;
		}
		my $latest = $c->latest($tool)->{version};
		$r{$tool} = { installed => $have->{$tool}, latest => $latest,
			update => (!defined $have->{$tool} || _vcmp($latest, $have->{$tool}) > 0) ? 1 : 0 };
	}
	\%r;
}

sub latest {
	my ($c, $tool) = @_;
	return _blast_latest() if $tool eq 'blast';
	return _clustalo_latest() if $tool eq 'clustalo';
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
			_run($c->python, qw(-m pip install --upgrade --disable-pip-version-check), sort keys %{ $todo->{$tool} });
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
	} elsif ($tool eq 'python') {
		my $venv = catdir($root, 'venv');
		_run(_base_python(), '-m', 'venv', '--system-site-packages', $venv);
		_tidy_venv($venv) unless $WIN;
		my $py = $c->python($root);
		_run($py, qw(-m pip install --disable-pip-version-check), @PY_PACKAGES);
		_run($py, '-c', 'import ' . join(',', @PY_PACKAGES));
		$c->_save($root, python => _line($py, '-c', 'import platform; print(platform.python_version())'));
	} else {
		croak "install() knows blast, clustalo and python, not \"$tool\"";
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
# $bin instead of downloading. ALIEN_BIOINF_CLUSTALO / ALIEN_BIOINF_BLAST name
# one explicitly (the clustalo binary itself, or BLAST's bin directory); failing
# that, PATH is searched. BLAST is a directory of ~30 programs, so a PATH hit is
# only copied from a directory that is plainly an unpacked NCBI tarball
# (.../ncbi-blast-<v>+/bin) -- never, say, all of /usr/bin.
sub _local_copy {
	my ($prog, $version, $bin, $whole_dir) = @_;
	my $env = $ENV{ $prog eq 'blastp' ? 'ALIEN_BIOINF_BLAST' : 'ALIEN_BIOINF_CLUSTALO' };
	my @cand = defined $env ? ($whole_dir ? catfile($env, "$prog$EXE") : $env) : ();
	push @cand, map { catfile($_, "$prog$EXE") } File::Spec->path unless @cand;
	foreach my $exe (grep { -f $_ && -x _ } @cand) {
		next if $whole_dir && !defined $env && dirname($exe) !~ /ncbi-blast-\Q$version\E\+[\/\\]bin\z/;
		my $out = eval { _run_out($exe, $prog eq 'blastp' ? '-version' : '--version') } // next;
		next unless $out =~ /^\s*(?:blastp:\s*)?\Q$version\E\+?\s*$/m;
		make_path($bin);
		my @files = $whole_dir ? _dir_files(dirname($exe)) : ($exe);
		foreach my $f (@files) {
			my $to = catfile($bin, $whole_dir ? basename($f) : "clustalo$EXE");
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

sub _http {
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

Alien::Bioinf - module-local Clustal Omega, NCBI BLAST+ and Python for Bioinf::Basic

=head1 SYNOPSIS

 use Alien::Bioinf;
 system Alien::Bioinf->clustalo, '--in', 'seqs.fa', '--out', 'aln.fa';
 system Alien::Bioinf->blast('blastp'), '-query', 'q.fa', '-subject', 's.fa';
 system Alien::Bioinf->python, 'script.py';

 my $behind = Alien::Bioinf->check_updates;          # every tool
 Alien::Bioinf->update('blast', 'python-packages');  # just these

From the shell:

 perl -MAlien::Bioinf -MData::Dumper -e 'print Dumper(Alien::Bioinf->check_updates)'
 perl -MAlien::Bioinf -e 'print "$_\n" for Alien::Bioinf->update'

=head1 DESCRIPTION

Installing this distribution puts the newest Clustal Omega and BLAST+ for the
running operating system, and a Python venv holding numpy, matplotlib and
adjustText, into its own share directory. None of it touches PATH.

=head2 Nothing is downloaded twice

Downloads are kept in L</cache_dir> and reused from there by later installs
and updates. Before downloading clustalo or BLAST, the installer looks for a
copy already on the machine at the version it wants and copies that instead:
C<ALIEN_BIOINF_CLUSTALO> may name a clustalo binary and C<ALIEN_BIOINF_BLAST>
the C<bin> directory of an unpacked BLAST+, and otherwise C<PATH> is searched.
The venv sees the base interpreter's packages, so pip fetches only what that
lacks.

=head1 METHODS

=over

=item clustalo, blast($program), python

Full paths to the installed programs. C<blast> defaults to C<blastp>.

=item versions

A hash reference of tool => installed version.

=item check_updates(@tools)

For each of C<clustalo>, C<blast>, C<python> and C<python-packages> (all four
by default), what is installed and what the upstream site offers. Each tool
maps to C<< { installed, latest, update } >>; C<python-packages> maps each
outdated package to the same.

=item update(@tools)

Installs whatever C<check_updates> reports as behind and returns the names of
the tools that changed. A newer Python needs a newer interpreter from
L<Alien::CPython3> first, since a venv cannot change interpreter.

=item latest($tool)

The newest upstream version of C<clustalo>, C<blast> or C<python>.

=item cache_dir

C<ALIEN_BIOINF_CACHE>, else C<alien-bioinf> under C<XDG_CACHE_HOME>, else
under C<~/.cache> (C<%LOCALAPPDATA%> on Windows).

=back

=cut
