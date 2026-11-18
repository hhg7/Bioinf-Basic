use strict;
use warnings;
use ExtUtils::MakeMaker;
use File::ShareDir::Install;
use Config;
use File::Temp ();
use File::Spec ();

install_share dist => 'share';

# Bioinf.xs declares variables in for() and after statements, which is C99. The
# flag that asks for it is per-vendor and cannot be told from $Config{cc}
# (that is plain "cc" for Solaris Studio and HP-UX aCC), so each candidate is
# trial-compiled -- the same probe as ~/Scripts/stats. MSVC needs none and
# only warns on flags it does not know, which would fool the probe.
my $std = '';
unless ($Config{cc} =~ /(?:^|[\/\\])cl(?:\.exe)?$/i) {
	my $dir = File::Temp::tempdir(CLEANUP => 1);
	my $src = File::Spec->catfile($dir, 'c99.c');
	open my $fh, '>', $src or die "$src: $!";
	print {$fh} "int main(void) { int s = 0; for (int i = 0; i < 3; i++) s += i; int t = s; return t - 3; }\n";
	close $fh;
	my $null = File::Spec->devnull;
	foreach my $flag ('-std=gnu99', '-std=c99', '-xc99=all', '-qlanglvl=extc99', '-AC99', '') {
		# quoted, since a temporary directory can have a space in it (a Windows user name)
		next if system("$Config{cc} $flag -c \"$src\" -o \"" . File::Spec->catfile($dir, 'c99.o') . "\" > $null 2> $null") != 0;
		$std = " $flag" if $flag ne '';
		last;
	}
}

# On a long-double perl, admixture's fit in Bioinf.xs wants logl() and sqrtl(),
# so that an NV is not computed at a double's precision. They are C99, but some
# libms (the thinner BSD ones especially) lack them, so they are link-tested,
# and without them Bioinf.xs uses log() and sqrt(). __float128 needs no probe:
# quadmath.h and -lquadmath come with the quadmath perl itself. This is the
# probe of ~/Scripts/stats, cut down to the two functions used.
my $define = '';
if ($Config{nvtype} eq 'long double') {
	my $dir = File::Temp::tempdir(CLEANUP => 1);
	my $src = File::Spec->catfile($dir, 'ldmath.c');
	my $bin = File::Spec->catfile($dir, 'ldmath');
	if (open my $ld_fh, '>', $src) {
		print {$ld_fh} "#include <math.h>\nint main(void) { volatile long double x = 2.0L; volatile long double r = logl(x) + sqrtl(x); return (int)r * 0; }\n";
		close $ld_fh;
		my $quiet = File::Spec->devnull;
		$define = '-DBIOINF_HAVE_LONG_DOUBLE_MATH'
			if system("$Config{cc}$std -o \"$bin\" \"$src\" -lm > $quiet 2> $quiet") == 0;
	}
}

WriteMakefile(
	NAME             => 'Cn::Bioinf',
	DISTNAME         => 'Cn-Bioinf',
	VERSION_FROM     => 'lib/Cn/Bioinf.pm',
	# Only the module. Without PM, MakeMaker also installs every *.pm and *.pl in
	# the root, which holds the maintainer's scratch scripts: on 2026-11-02 a
	# build here put md2pod.pl and nine others in blib/lib/Cn/. The tarball has
	# none, since [PruneFiles] leaves them out.
	PM               => { 'lib/Cn/Bioinf.pm' => '$(INST_LIB)/Cn/Bioinf.pm' },
	ABSTRACT         => 'FASTA I/O in XS, BLAST hit ranking, and alignment plots and tables', # the POD has no NAME section; PodWeaver makes one from '# ABSTRACT:'
	AUTHOR           => 'David E. Condon <dcondon@uidaho.edu>',
	LICENSE          => 'perl_5',
	MIN_PERL_VERSION => '5.012',
	CCFLAGS          => $Config{ccflags} . $std,
	DEFINE           => $define,
	LIBS             => ['-lm'],
	CONFIGURE_REQUIRES => { 'ExtUtils::MakeMaker' => '6.64', 'File::ShareDir::Install' => '0.06' },
	PREREQ_PM        => { 'File::ShareDir' => '1.00', 'HTTP::Tiny' => '0.056', 'JSON::MaybeXS' => '0' },
	TEST_REQUIRES    => { 'Test::More' => '0.88' },
	# recommended, not required: see the [Prereqs / RuntimeRecommends] comments in dist.ini
	META_MERGE       => { 'meta-spec' => { version => 2 }, prereqs => { runtime => { recommends => { 'Alien::Bioinf' => '0', 'IO::Socket::SSL' => '1.968', 'Net::SSLeay' => '1.49' } } } },
);

package MY;
use File::ShareDir::Install 'postamble';
