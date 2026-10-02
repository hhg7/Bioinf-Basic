use strict;
use warnings;
use ExtUtils::MakeMaker;
use File::ShareDir::Install;
use Config;
use File::Temp ();
use File::Spec ();

install_share dist => 'share';

# Basic.xs declares variables in for() and after statements, which is C99. The
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
		next if system("$Config{cc} $flag -c $src -o " . File::Spec->catfile($dir, 'c99.o') . " > $null 2> $null") != 0;
		$std = " $flag" if $flag ne '';
		last;
	}
}

WriteMakefile(
	NAME             => 'Bioinf::Basic',
	DISTNAME         => 'Bioinf-Basic',
	VERSION_FROM     => 'lib/Bioinf/Basic.pm',
	ABSTRACT         => 'FASTA I/O in XS, BLAST hit ranking, and alignment plots and tables', # the POD has no NAME section; PodWeaver makes one from '# ABSTRACT:'
	AUTHOR           => 'David E. Condon <dcondon@uidaho.edu>',
	LICENSE          => 'perl_5',
	MIN_PERL_VERSION => '5.010',
	CCFLAGS          => $Config{ccflags} . $std,
	CONFIGURE_REQUIRES => { 'ExtUtils::MakeMaker' => '6.64', 'File::ShareDir::Install' => '0.06' },
	PREREQ_PM        => { 'File::ShareDir' => '1.00', 'JSON::MaybeXS' => '0' },
	TEST_REQUIRES    => { 'Test::More' => '0.88' },
	# recommended, not required: see the [Prereqs / RuntimeRecommends] comment in dist.ini
	META_MERGE       => { 'meta-spec' => { version => 2 }, prereqs => { runtime => { recommends => { 'Alien::Bioinf' => '0' } } } },
);

package MY;
use File::ShareDir::Install 'postamble';
