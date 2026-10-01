#!/usr/bin/env perl
# test.all.perls.pl - build & test Bioinf::Basic against every perlbrew perl.
#
# Adapted on 2026-09-30 from ~/Scripts/stats/test.all.perls.pl, which does the
# same for Stats::LikeR; a fix to the shared machinery (discovery, the private
# trees, the fork/reap loop, the summary) belongs in both copies.  What differs:
#
#   * Bioinf::Basic needs Alien::Bioinf, the sibling distribution in
#     alien/Alien-Bioinf, and every perl has its own site_perl to install it
#     into.  So, as compile.sh does, a perl that cannot load Alien::Bioinf gets
#     it installed with cpanm first (see the 'alien' step in build_one).
#   * The x87 row, the quadmath-first run order and the 32-bit-IV advice are
#     gone.  They exist in Stats::LikeR because its XS does floating-point
#     arithmetic whose result depends on how wide an NV is and how it is
#     computed; Basic.xs does none, so those axes have nothing to catch here.
#
# Automates the manual loop of:
#     perlbrew use perl-5.10.1 && ./compile.sh
#     perlbrew use perl-5.12.5 && ./compile.sh
#     perlbrew use perl-5.44.0 && ./compile.sh
#
# Each version is built in its own environment (no `perlbrew use` needed: the
# perl binary is invoked directly and PATH is rewritten for the child), the
# tree is cleaned between versions, and a pass/fail summary is printed at the
# end.  Exit status is non-zero if any version failed.
#
# By default four perls are built at once (-P 4).  They cannot share the
# distribution root: there is one Makefile, one Basic.c, one Basic.o and one
# blib, and a Basic.o built by another perl links without complaint and then
# dies at load.  So each parallel child builds in a private copy of the tree
# under the log directory and reports its result back to the parent, which
# leaves the root's Makefile and blib exactly as they were.  Nothing forces
# Parallel::ForkManager here: the fork/throttle/reap loop is a few dozen lines
# and this script already forks in run_cmd, so the helper keeps needing nothing
# but core modules.
#
# -P 1 is the old behaviour: one perl at a time, in the distribution root,
# which is left built and installed against the perl that runs last -- the
# newest plain one, see the run order below.

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use Getopt::Long 'GetOptions';
# -P (how many perls at once) and -p (which perl) are different options, so the
# default case-folding of single-letter aliases has to go.
Getopt::Long::Configure('no_ignore_case');
use File::Spec;
use File::Copy 'copy';
use File::Path 'remove_tree';
use Fcntl ':flock';
use IO::Handle;
use Cwd 'getcwd';
use Data::Dumper;
use POSIX 'strftime';
use Time::HiRes 'time';

my $PERLBREW_ROOT = $ENV{PERLBREW_ROOT} || File::Spec->catdir($ENV{HOME}, 'perl5', 'perlbrew');

# How many cores this machine has, or 0 when that cannot be determined.  Only
# an author-side helper runs this, but the probes stay portable so a run on a
# BSD or Solaris box degrades to a sensible default instead of dying.
sub cpu_count {
	if (open my $fh, '<', '/proc/cpuinfo') {           # linux
		my $n = grep { /^processor\s*:/ } <$fh>;
		close $fh;
		return $n if $n;
	}
	for my $probe (['getconf', '_NPROCESSORS_ONLN'],   # glibc, solaris, aix
		['sysctl', '-n', 'hw.ncpu']) {             # *bsd, darwin
		my $out = `@$probe 2>/dev/null`;
		return $1 if defined $out && $out =~ /(\d+)/ && $1 > 0;
	}
	return 0;
}

my ($help, $list, $install, $deps, $clean, $stop, $quiet, $jobs, $log_dir, $optimize, @only);
my ($par, $keep_work);
# On by default, as in compile.sh: without Alien::Bioinf, t/blast.t and
# t/msa.t fail, so a perl that lacks it would only report a missing prerequisite.
my $alien = 1;
my $ncpu  = cpu_count();
$install  = 1; # make install, as compile.sh does
$clean    = 1;
# OPTIMIZE= replaces perl's own $Config{optimize} rather than adding to it, so
# a bare '-Wall' here would build every perl in the matrix at -O0.  The same
# string as compile.sh.
$optimize = '-O2 -Wall';
$log_dir  = File::Spec->catdir('.build', 'multiperl');
# Parallel by default, leaving 4 cores for the rest of the machine: a full
# matrix run is minutes of `make` and a test harness, not something that should
# need remembering to ask for.
my $budget = $ncpu > 5 ? $ncpu - 4 : 1;
# That budget is split two ways: several perls at once, each with its own
# parallel make and harness.  -P x -j multiply, so -j follows from -P rather
# than claiming the budget twice.
$par = 4;
my $par_default  = $par;
my $jobs_default = jobs_for($par);

sub jobs_for {
	my $p = shift;
	my $j = int($budget / ($p || 1));
	return $j > 1 ? $j : 1;
}

GetOptions(
	'perl|p=s@'     => \@only,
	'install!'      => \$install,
	'deps!'         => \$deps,
	'alien!'        => \$alien,
	'clean!'        => \$clean,
	'stop-on-fail!' => \$stop,
	'jobs|j=i'      => \$jobs,
	'parallel|P:i'  => \$par,
	'keep-work!'    => \$keep_work,
	'optimize=s'    => \$optimize,
	'log-dir=s'     => \$log_dir,
	'quiet|q'       => \$quiet,
	'list|l'        => \$list,
	'help|h'        => \$help,
) or usage(1);
usage(0) if $help;

sub usage {
	my $rc = shift;
	my $cpus = $ncpu ? "$ncpu CPUs - 4" : 'CPU count unknown';
	print STDERR <<"END";
usage: $0 [options]

Builds and tests the distribution in the current directory against each
perl installed under $PERLBREW_ROOT.

Run order: oldest first, then the newest plain perl (double NV, no ithreads)
last, so a serial run leaves the tree built and installed against the reference
build rather than an outlier.

A perl that cannot load Alien::Bioinf has it installed first, with cpanm, from
a fresh copy of alien/Alien-Bioinf.  That downloads clustalo, BLAST+ and a
Python into ~/.cache/alien-bioinf the first time; installs are run one at a
time, even under -P, because the perls share that cache.

options:
  -p, --perl VERSION   only this perl (repeatable); accepts "5.10.1",
                       "perl-5.10.1" or an exact directory name such as
                       "5.44.0-quadmath".  default: every installed perl
  -l, --list           list the perls that would be tested, in run order, with
                       each one's NV width and threading, then exit
      --no-install     skip "make install" (build + test only)
      --no-clean       skip "make clean" before each version
      --deps           cpanm any missing CONFIGURE_REQUIRES, PREREQ_PM or
                       TEST_REQUIRES for that perl first
      --no-alien       do not install Alien::Bioinf where it is missing
                       (t/blast.t and t/msa.t will then fail on that perl)
      --stop-on-fail   abort at the first version that fails (with -P: launch
                       no more perls; the running ones finish)
  -P, --parallel [N]   build and test N perls at once, each in a private copy of
                       the tree under --log-dir/work, leaving the distribution
                       root untouched.  bare -P, or -P 0, uses one child per
                       perl, capped at the CPU count; -P 1 builds serially in
                       the distribution root and leaves it built against the
                       newest plain perl.  default: $par_default
      --keep-work      keep the private trees of a perl that passed (those of
                       one that failed are always kept, for inspection)
  -j, --jobs N         parallel make, and HARNESS_OPTIONS=j<N> for the tests;
                       0 builds and tests serially.  -P and -j multiply, so the
                       default is the CPU budget ($cpus) divided by -P
                       (default: $jobs_default at -P $par_default, $budget at -P 1)
      --optimize STR   OPTIMIZE= passed to Makefile.PL (default: $optimize)
      --log-dir DIR    where per-version logs go (default: $log_dir)
  -q, --quiet          only write logs; do not echo build output.  implied by
                       -P > 1, where interleaved output would be unreadable
  -h, --help           this message

exit status: 0 if every perl built, tested and installed cleanly, else 1.
END
	exit $rc;
}

# ---------------------------------------------------------------- discovery --

my $perls_dir = File::Spec->catdir($PERLBREW_ROOT, 'perls');
die "$0: no perlbrew perls directory at $perls_dir\n" unless -d $perls_dir;

opendir my $dh, $perls_dir or die "$0: cannot read $perls_dir: $!\n";
my @installed = grep { -x File::Spec->catfile($perls_dir, $_, 'bin', 'perl') }
	grep { !/^\.\.?$/ } readdir $dh;
closedir $dh;

# numeric sort, oldest first, which is where the run order below starts from.
sub vkey {
	(my $v = shift) =~ s/^perl-//;
	my @p = ($v =~ /(\d+)/g);
	push @p, 0 while @p < 3;
	# the whole name breaks a tie, so 5.44.0, 5.44.0-i686 and 5.44.0-quadmath
	# keep one order from run to run
	return sprintf('%05d%05d%05d', @p[0 .. 2]) . $v;
}
@installed = sort { vkey($a) cmp vkey($b) } @installed;

my @targets = @installed;
if (@only) {
	my %have = map { $_ => 1 } @installed;
	my (@want, @missing);
	for my $arg (map { split /,/ } @only) {
		# an exact directory name wins, so a build installed with no perl-
		# prefix, such as "5.44.0-quadmath", is selectable.
		my ($match) = grep { $have{$_} } $arg, "perl-$arg";
		if (defined $match) { push @want, $match }
		else                { push @missing, $arg }
	}
	die "$0: not installed under perlbrew: @missing\n(installed: @installed)\n" if @missing;
	my %seen;
	@targets = sort { vkey($a) cmp vkey($b) } grep { !$seen{$_}++ } @want;
}
die "$0: no perls found in $perls_dir\n" unless @targets;

# ------------------------------------------------------------------ run order --

# What each perl actually is.  The answer has to come from the interpreter:
# "5.44.0-quadmath" is a local naming habit rather than a promise, and the
# long-double build is only identifiable as "perl-5.12.5" plus an archname
# nobody parses.  A perl that will not answer is described as unknown, never as
# plain, so a broken interpreter cannot become the reference build the tree is
# left standing on.
#
# ivsize decides nothing about the run order.  It is reported because
# 5.44.0-i686 is a 32-bit perl, the only one here where STRLEN and SSize_t --
# what Basic.xs's FASTA reader and writer count lengths in -- are 32 bits wide,
# and the label makes that row look different.  It is also not plain, so it is
# never the build a serial run leaves the tree standing on.
my %facts;
sub facts {
	my $version = shift;
	return $facts{$version} if $facts{$version};
	my $perl = File::Spec->catfile($perls_dir, $version, 'bin', 'perl');
	my %f = (nv => '?', iv => 0, threads => 0, known => 0);
	if (open my $fh, '-|', $perl, '-MConfig', '-e',
			'print "$Config{nvtype}\t", ($Config{useithreads} ? 1 : 0),'
			. '"\t$Config{ivsize}"') {
		my $line = <$fh>;
		close $fh;
		if (defined $line && $line =~ /^(\S[^\t]*)\t([01])\t(\d+)/) {
			%f = (nv => $1, threads => $2, iv => $3, known => 1);
		}
	}
	# plain: the double-NV, unthreaded build that the long-double, quadmath and
	# threaded configurations are each a variation on.  Unknown does not count.
	$f{plain} = $f{known} && $f{nv} eq 'double' && !$f{threads} && $f{iv} == 8;
	return $facts{$version} = \%f;
}

sub nv_label {
	my $f = shift;
	my %short = ('double' => 'double', 'long double' => 'long-double',
		'__float128' => 'quadmath');
	return ($short{ $f->{nv} } || $f->{nv}) . ($f->{threads} ? '-thr' : '')
		. ($f->{known} && $f->{iv} != 8 ? "/iv$f->{iv}" : '');
}

# Run order, given @targets oldest first: everything else, still oldest first,
# then the newest plain perl last, because a serial run leaves the distribution
# root built and installed against whichever perl went last, and that should be
# the ordinary double-NV reference build rather than a long-double, threaded,
# quadmath or 32-bit outlier.
sub order_targets {
	my @rest = @_;
	my ($plain) = grep { facts($_)->{plain} } reverse @rest;
	@rest = grep { $_ ne $plain } @rest if defined $plain;
	return (@rest, defined $plain ? $plain : ());
}
@targets = order_targets(@targets);

if ($list) {
	printf "%-20s %-16s %s\n", $_, nv_label(facts($_)),
		File::Spec->catfile($perls_dir, $_, 'bin', 'perl') for @targets;
	exit 0;
}

die "$0: no Makefile.PL in " . getcwd() . " - run this from the distribution root\n"
	unless -f 'Makefile.PL';
my $alien_src = File::Spec->catdir('alien', 'Alien-Bioinf');
die "$0: no $alien_src in " . getcwd() . " (or pass --no-alien)\n"
	if $alien && !-f File::Spec->catfile($alien_src, 'Makefile.PL');

# ------------------------------------------------------------------ logging --

mkdirp($log_dir);
my $stamp = strftime '%Y%m%d-%H%M%S', localtime;

# ------------------------------------------------------------- parallelism --

# bare -P (or -P 0) means "all of them", within reason: more compilers than
# cores only makes every perl slower.
$par = 0 if $par < 0;
$par ||= do {
	my $cap = $ncpu || 4;   # an unknown CPU count stays conservative rather
	                        # than fanning out over however many perls exist
	@targets < $cap ? scalar @targets : $cap;
};
$par = @targets if $par > @targets;
$jobs = jobs_for($par) unless defined $jobs;   # -j follows -P unless asked for

# Children chdir into their own tree, and the alien step chdirs nowhere but
# hands cpanm a path, so the log paths are made absolute in every mode.
$log_dir = File::Spec->rel2abs($log_dir);
my $work_root   = File::Spec->catdir($log_dir, 'work');
my $alien_lock  = File::Spec->catfile($log_dir, 'alien.lock');
my $in_parallel = 0;   # set by the dispatch below, once it knows there is
                       # more than one perl to run

sub mkdirp {
	my @parts = File::Spec->splitdir(shift);
	my $path;   # undef, not '': splitdir gives an absolute path a leading '',
	            # and catdir('', 'home') is '/home' where 'home' would be a
	            # directory of that name in the cwd.
	for my $p (@parts) {
		$path = defined $path ? File::Spec->catdir($path, $p) : $p;
		next if !length($path) || -d $path;
		mkdir $path or die "$0: mkdir $path: $!\n";
	}
}

# ------------------------------------------------------------- child runner --

# Run @cmd with STDERR folded into STDOUT, echoing to the terminal and to
# $logfh.  Returns ($exit_code, \@lines).
sub run_cmd {
	my ($cmd, $logfh) = @_;
	print $logfh "\n\$ @$cmd\n";
	print "\$ @$cmd\n" unless $quiet;

	my $pid = open my $fh, '-|';
	die "$0: fork failed: $!\n" unless defined $pid;
	if (!$pid) {                              # child
		open STDERR, '>&', \*STDOUT or die "$0: dup STDERR: $!\n";
		$| = 1;
		{ exec { $cmd->[0] } @$cmd; }
		print "exec @$cmd failed: $!\n";
		exit 127;
	}

	my @lines;
	while (defined(my $line = <$fh>)) {
		push @lines, $line;
		print $logfh $line;
		print $line unless $quiet;
	}
	close $fh;
	my $status = $?;
	my $code = $status == -1  ? -1
		: ($status & 127) ? 128 + ($status & 127)
		: ($status >> 8);
	return ($code, \@lines);
}

# Module names from the CONFIGURE_REQUIRES, PREREQ_PM and TEST_REQUIRES blocks
# of Makefile.PL, so this stays in sync with it.  CONFIGURE_REQUIRES is included
# because Makefile.PL itself loads File::ShareDir::Install, which a fresh perl
# lacks.  Alien::Bioinf is left out: it is not on CPAN, and the 'alien' step
# installs it from alien/Alien-Bioinf instead.
sub prereqs {
	open my $fh, '<', 'Makefile.PL' or return ();
	local $/;
	my $src = <$fh>;
	close $fh;
	my %seen = ('Alien::Bioinf' => 1);
	my @mods;
	while ($src =~ /\b(?:CONFIGURE_REQUIRES|PREREQ_PM|TEST_REQUIRES)\s*=>\s*\{(.*?)\}/sg) {
		push @mods, grep { !$seen{$_}++ } ($1 =~ /['"]([\w:]+)['"]\s*=>/g);
	}
	return @mods;
}
my @prereqs = prereqs();

sub cpanm_for {
	my $bin = shift;
	my $own = File::Spec->catfile($bin, 'cpanm');
	return -x $own ? $own : File::Spec->catfile($PERLBREW_ROOT, 'bin', 'cpanm');
}

# ------------------------------------------------------------ private trees --

# What must not be copied into a private build tree: the build products of
# whichever perl last used the source tree (copying them would recreate exactly
# the stale-Basic.o hazard the private trees exist to avoid), the repository, a
# dzil build directory, and the log/work tree itself.  The names are
# .gitignore's, which also covers alien/Alien-Bioinf: its _alien/ holds
# Alien::Build's state for the perl that last built it, and reusing that state
# under another perl is the same hazard one level down.
my %SKIP_DIR  = map { $_ => 1 } qw(.git .build blib _alien);
my @SKIP_FILE = (qr/\.(?:o|a|so|bs|dylib|dll)$/, qr/^Basic\.c$/,
	qr/^Makefile(?:\.old)?$/, qr/^pm_to_blib$/, qr/^blibdirs$/,
	qr/^MYMETA\./, qr/^Bioinf-Basic-.*\.tar\.gz$/);

sub copy_tree {
	my ($src, $dst) = @_;
	mkdirp($dst);
	opendir my $dh, $src or die "$0: cannot read $src: $!\n";
	my @entries = grep { !/^\.\.?$/ } readdir $dh;
	closedir $dh;
	for my $e (@entries) {
		my $s = File::Spec->catfile($src, $e);
		my $d = File::Spec->catfile($dst, $e);
		if (-d $s) {
			next if $SKIP_DIR{$e} || $e =~ /^Bioinf-Basic-[\d.]+$/;
			# a --log-dir inside the tree would copy itself forever
			my $abs = File::Spec->rel2abs($s);
			next if $abs eq $work_root || $abs eq $log_dir;
			copy_tree($s, $d);
			next;
		}
		next if grep { $e =~ $_ } @SKIP_FILE;
		copy($s, $d) or die "$0: copy $s -> $d: $!\n";
		chmod((stat $s)[2] & 07777, $d);
	}
}

# ------------------------------------------------------------ build one perl --

# Install Alien::Bioinf into $perl if it cannot load it, and return the step
# record, or nothing when it was already there.  cpanm builds a local directory
# in place, so it is handed a fresh copy of alien/Alien-Bioinf rather than the
# one in the source tree, which may hold another perl's Makefile and _alien/.
#
# The lock serialises the installs, in every mode: Alien::Bioinf downloads
# into one cache directory shared by every perl, and two children fetching the
# same archive would write the same .part file at once.  Only the install is
# serialised; the perls build and test in parallel around it.
sub install_alien {
	my ($version, $perl, $bin, $logfh) = @_;
	my ($have) = run_cmd([$perl, '-MAlien::Bioinf', '-e', '1'], $logfh);
	return if $have == 0;

	my $t = time;
	my $work = File::Spec->catdir($work_root, "alien.$version.$stamp");
	remove_tree($work) if -d $work;
	copy_tree($alien_src, $work);

	open my $lock, '>', $alien_lock or die "$0: cannot write $alien_lock: $!\n";
	unless (flock $lock, LOCK_EX | LOCK_NB) {
		print $logfh "-- waiting for another perl's Alien::Bioinf install\n";
		flock $lock, LOCK_EX or die "$0: flock $alien_lock: $!\n";
	}
	my ($code) = run_cmd([$perl, cpanm_for($bin), $work], $logfh);
	close $lock;   # releases the lock

	if ($code == 0 && !$keep_work) { remove_tree($work) }
	else { print $logfh "-- Alien::Bioinf build tree kept at $work\n" }
	return { label => 'alien', code => $code, seconds => time - $t };
}

# Build, test and install the distribution in the current directory with one
# perl.  Returns the result hashref, and prints its own progress unless
# $silent: a parallel child is silent and the parent reports for it.
sub build_one {
	my ($version, $silent) = @_;
	my $root = File::Spec->catdir($perls_dir, $version);
	my $bin  = File::Spec->catdir($root, 'bin');
	my $perl = File::Spec->catfile($bin, 'perl');

	my $log = File::Spec->catfile($log_dir, "$version.$stamp.log");
	# Re-create the log directory rather than trusting that it survived.  A
	# matrix run is minutes long, and anything that removes .build meanwhile
	# turns a perl that was building fine into a failure with no log to read.
	mkdirp($log_dir);
	open my $logfh, '>', $log or die "$0: cannot write $log: $!\n";
	$logfh->autoflush(1);
	STDOUT->autoflush(1);

	unless ($silent) {
		print "\n", '=' x 72, "\n";
		printf "== %s   (log: %s)\n", $version, $log;
		print '=' x 72, "\n";
	}
	print $logfh "== $version at " . strftime('%F %T', localtime)
		. ' in ' . getcwd() . "\n";

	# Emulate `perlbrew use $version` for the children: this perl's bin first,
	# every other perlbrew perl stripped out, and local::lib / PERL5LIB
	# leftovers from the calling shell removed so nothing bleeds across
	# versions.
	local %ENV = %ENV;
	my @path = grep { index($_, File::Spec->catdir($perls_dir, '')) != 0 }
		split /:/, ($ENV{PATH} || '/usr/bin:/bin');
	$ENV{PATH}          = join ':', $bin, @path;
	$ENV{PERLBREW_ROOT} = $PERLBREW_ROOT;
	$ENV{PERLBREW_PERL} = $version;
	$ENV{PERLBREW_PATH} = $bin;
	delete @ENV{qw(PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MM_OPT PERL_MB_OPT
		PERLBREW_LIB PERL_MM_USE_DEFAULT)};
	$ENV{HARNESS_OPTIONS} = "j$jobs" if $jobs;

	my $t0 = time;
	my %r = (version => $version, log => $log, steps => [], warnings => 0);

	# sanity: the interpreter really is the version we think it is, and say so
	# when it is a threaded build -- an XS bug can be invisible on an
	# unthreaded perl and a hard compile error under MULTIPLICITY, so the
	# summary has to make that coverage visible rather than implied by a name.
	my ($vc, $vout) = run_cmd([$perl, '-MConfig', '-e',
		'printf "%vd%s\n", $^V, $Config{useithreads} ? "-thr" : ""'], $logfh);
	$r{reported} = $vc == 0 && @$vout ? do { my $s = $vout->[0]; chomp $s; $s } : '?';

	my $failed;
	my @steps;
	if ($deps && @prereqs) {
		my @missing;
		for my $mod (@prereqs) {
			my ($c) = run_cmd([$perl, "-M$mod", '-e', '1'], $logfh);
			push @missing, $mod if $c != 0;
		}
		push @steps, ['deps', [$perl, cpanm_for($bin), '--notest', @missing]] if @missing;
	}

	# The deps step has to run before the Alien::Bioinf install, which is not a
	# plain command and so is not in @steps: run it here, ahead of the rest.
	for my $step (splice @steps) {
		my ($label, $cmd) = @$step;
		my $t_step = time;
		my ($code) = run_cmd($cmd, $logfh);
		push @{ $r{steps} }, { label => $label, code => $code, seconds => time - $t_step };
		if ($code != 0) { $failed = $label; last }
	}
	if ($alien && !$failed) {
		if (my $s = install_alien($version, $perl, $bin, $logfh)) {
			push @{ $r{steps} }, $s;
			$failed = 'alien' if $s->{code} != 0;
		}
	}

	if ($clean && -f 'Makefile') {
		push @steps, ['clean', ['make', 'clean'], 1];   # 1 = failure tolerated
	}

	push @steps, ['Makefile.PL', [$perl, 'Makefile.PL', "OPTIMIZE=$optimize"]];
	push @steps, ['make',        ['make', $jobs ? ("-j$jobs") : ()]];
	push @steps, ['make test',   ['make', 'test']];
	push @steps, ['make install',['make', 'install']] if $install;

	for my $step ($failed ? () : @steps) {
		my ($label, $cmd, $soft) = @$step;
		my $t_step = time;
		my ($code, $lines) = run_cmd($cmd, $logfh);
		my $secs = time - $t_step;
		printf $logfh "-- step '%s' exited %d after %.1fs\n", $label, $code, $secs;

		if ($label eq 'make test') {
			for my $l (@$lines) {
				$r{files}  = $1 if $l =~ /^(Files=\d+.*)/;
				$r{result} = $1 if $l =~ /^Result:\s*(\S+)/;
				$r{passed} = 1  if $l =~ /^All tests successful/;
			}
		}
		# residual build warnings are worth surfacing even on success
		$r{warnings} += grep { /: warning:/ } @$lines if $label eq 'make';

		push @{ $r{steps} }, { label => $label, code => $code, seconds => $secs };
		next if $code == 0 || $soft;
		$failed = $label;
		last;
	}

	if ($clean && !$failed && !$in_parallel) {
		# leave a clean tree behind before the next perl takes over.  in
		# parallel mode the whole private tree goes instead, so this would
		# only be a `make clean` nobody reads.
		run_cmd(['make', 'clean'], $logfh) unless $version eq $targets[-1];
	}

	$r{seconds} = time - $t0;
	$r{failed}  = $failed;
	close $logfh;

	report_one(\%r) unless $silent;
	return \%r;
}

# the two progress lines a finished perl prints, from the parent in either mode
sub report_one {
	my $r = shift;
	printf "-- %s: %s in %.1fs%s\n", $r->{version},
		($r->{failed} ? "FAILED at '$r->{failed}'" : 'ok'),
		$r->{seconds},
		(defined $r->{files} ? " ($r->{files})" : '');
	print '-- ', join('  ', map { sprintf '%s %.1fs', $_->{label}, $_->{seconds} }
		@{ $r->{steps} }), "\n";
}

# The child's result has to cross a fork, and %r is plain data, so Dumper out /
# eval in beats any IPC here: it also leaves the numbers next to the log when
# something needs explaining afterwards.
sub write_result {
	my ($file, $r) = @_;
	# Every caller writes into $log_dir, and this one runs *after* the whole
	# build: a log directory that went away mid-run must cost the report, not
	# the minutes of work the report is about.
	mkdirp($log_dir);
	open my $fh, '>', $file or die "$0: cannot write $file: $!\n";
	local $Data::Dumper::Indent   = 0;
	local $Data::Dumper::Sortkeys = 1;
	print $fh Data::Dumper->Dump([$r], ['R']);
	close $fh or die "$0: close $file: $!\n";
}

sub read_result {
	my $file = shift;
	open my $fh, '<', $file or return undef;
	local $/;
	my $src = <$fh>;
	close $fh;
	my $R;
	eval $src;                      # our own Dumper output, nobody else's
	return ref $R eq 'HASH' ? $R : undef;
}

# ------------------------------------------------------------------ drivers --

sub run_serial {
	my @queue = @_;
	my @out;
	for my $version (@queue) {
		my $r = build_one($version);
		push @out, $r;
		next unless $r->{failed} && $stop;
		my %done = map { $_->{version} => 1 } @out;
		my @rest = grep { !$done{$_} } @queue;
		print "-- --stop-on-fail: skipping @rest\n" if @rest;
		last;
	}
	return @out;
}

# Fork up to $par children, each in its own copy of the tree, reaping them as
# they finish and starting the next perl in the freed slot.
sub run_parallel {
	my @order = @_;
	my @queue = @order;
	my $src   = getcwd();
	my %kid;                        # pid => { version, work, result, log }
	my (@out, $halt);

	my $reaper = sub {
		my ($pid, $status) = @_;
		my $kid = delete $kid{$pid} or return;
		my $r   = read_result($kid->{result});
		if (!$r) {
			# The child died without reporting: exec failure, signal, OOM, or
			# a log directory that disappeared under it.  Naming a log file
			# that was never created sends the reader looking for evidence
			# that is not there, so distinguish the two.
			my $log = -e $kid->{log} ? $kid->{log}
				: "$kid->{log} (never written)";
			$r = { version => $kid->{version}, log => $log,
				reported => '?', steps => [], warnings => 0, seconds => 0,
				failed => sprintf('child exited %d%s', $status >> 8,
					($status & 127) ? ' on signal ' . ($status & 127) : '') };
		}
		push @out, $r;
		report_one($r);
		$halt = 1 if $r->{failed} && $stop;
		if ($r->{failed} || $keep_work) {
			print "-- $kid->{version}: build tree kept at $kid->{work}\n";
		}
		else {
			remove_tree($kid->{work});
		}
	};

	local $SIG{INT} = local $SIG{TERM} = sub {
		my $sig = shift;
		print "\n-- $sig: stopping " . keys(%kid) . " running build(s)\n";
		# each child leads its own process group (see the fork below), so one
		# signal per group takes its make, its compilers and its harness with
		# it; signalling the child perl alone would orphan those.
		kill 'TERM', map { -$_ } keys %kid;
		sleep 1;
		kill 'KILL', map { -$_ } keys %kid;
		exit 130;
	};

	while (@queue || %kid) {
		while (@queue && keys(%kid) < $par && !$halt) {
			my $version = shift @queue;
			my $work    = File::Spec->catdir($work_root, "$version.$stamp");
			my $result  = File::Spec->catfile($log_dir, "$version.$stamp.result");
			my $log     = File::Spec->catfile($log_dir, "$version.$stamp.log");

			remove_tree($work) if -d $work;
			printf "-- %-16s starting (tree: %s)\n", $version, $work;
			copy_tree($src, $work);

			my $pid = fork;
			die "$0: fork failed: $!\n" unless defined $pid;
			if (!$pid) {                              # child
				$SIG{$_} = 'DEFAULT' for qw(INT TERM);
				setpgrp 0, 0;   # so ^C reaches this whole build, once, via
						# the parent's handler and not the terminal
				# A die in here would reach the parent as nothing but an exit
				# status, reported as a bare 'child exited N' with no steps
				# and no elapsed time -- indistinguishable from a build that
				# never started.  Catch it and put the reason where the parent
				# actually looks, which is the result file.
				my $r = eval {
					chdir $work or die "chdir $work: $!\n";
					build_one($version, 1);
				};
				unless ($r) {
					my $why = $@ || 'build_one returned nothing';
					chomp $why;   # $r->{failed} is printed inside quotes
					$r = { version => $version, log => $log,
						reported => '?', steps => [], warnings => 0,
						seconds => 0, failed => $why };
				}
				# last resort: the result file is unwritable, so the only
				# place left to say why is the terminal.
				eval { write_result($result, $r); 1 }
					or print STDERR "$0: $version: $@";
				exit 0;
			}
			$kid{$pid} = { version => $version, work => $work,
				result => $result, log => $log };
		}
		last unless %kid;
		my $pid = waitpid -1, 0;
		last if $pid <= 0;
		$reaper->($pid, $?);
	}

	print "-- --stop-on-fail: skipping @queue\n" if $halt && @queue;
	# report in version order, not in the order they happened to finish
	my %by = map { $_->{version} => $_ } @out;
	return map { $by{$_} } grep { $by{$_} } @order;
}

# --------------------------------------------------------------- main loop --

my @results;
my $t_all = time;

if ($par > 1 && @targets > 1) {
	$in_parallel = 1;
	$quiet       = 1;   # N interleaved build logs on one terminal is noise
	printf "-- %d perl(s), %d at a time%s; per-version output goes to the logs\n",
		scalar @targets, $par, ($jobs ? " with make -j$jobs each" : '');
	print "-- note: --deps has the children sharing one ~/.cpanm; if a "
		. "prerequisite install misbehaves, run once with -P 1 --deps first\n"
		if $deps;
	@results = run_parallel(@targets);
}
else {
	@results = run_serial(@targets);
}

# ----------------------------------------------------------------- summary --

my $bad = grep { $_->{failed} } @results;
print "\n", '=' x 72, "\n";
printf "%-16s %-11s %-8s %-7s %-6s %s\n",
	qw(PERL REPORTED STATUS TIME WARN TESTS);
print '-' x 72, "\n";
for my $r (@results) {
	printf "%-16s %-11s %-8s %6.1fs %-6s %s\n",
		$r->{version},
		$r->{reported},
		($r->{failed} ? 'FAIL' : 'PASS'),
		$r->{seconds},
		($r->{warnings} || 0),
		# Send the reader to the log only when there is one to read: if
		# whatever removed the log directory mid-run took the log with it,
		# "see <path>" is an invitation to hunt for a file that is not there.
		($r->{failed}
			? "failed at '$r->{failed}' - "
				. (-e $r->{log} ? "see $r->{log}"
					: "log gone: $r->{log}")
			: ($r->{result} ? "Result: $r->{result}" : 'no test summary')),
		;
}
print '-' x 72, "\n";
my $skipped = @targets - @results;
printf "%d/%d perl(s) passed%s in %.1fs.  Logs in %s\n",
	scalar(@results) - $bad, scalar(@targets),
	($skipped ? " ($skipped not run)" : ''),
	time - $t_all, $log_dir;
# every perl installed into its own site_perl, but nothing was built here: say
# so, because a serial run leaves the root built against the newest plain perl.
print "-- built in private trees; this directory is untouched (-P 1 to build here)\n"
	if $in_parallel;

exit($bad || $skipped ? 1 : 0);
