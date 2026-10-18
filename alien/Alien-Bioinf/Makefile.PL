use strict;
use warnings;
use ExtUtils::MakeMaker;
use Alien::Build::MM;

# Where NCBI builds no BLAST+ this cannot install, so it says so and exits 0
# without writing a Makefile, which a CPAN smoker reports as NA rather than
# FAIL. ALIEN_BIOINF_BLAST, a BLAST+ built locally, gets past that.
use lib 'lib';
require Alien::Bioinf;
unless (defined $ENV{ALIEN_BIOINF_BLAST} || eval { Alien::Bioinf::_blast_platform(); 1 }) {
	print "OS unsupported: $@";
	exit 0;
}

my $abmm = Alien::Build::MM->new;
WriteMakefile($abmm->mm_args(
	NAME             => 'Alien::Bioinf',
	DISTNAME         => 'Alien-Bioinf',
	VERSION_FROM     => 'lib/Alien/Bioinf.pm',
	ABSTRACT_FROM    => 'lib/Alien/Bioinf.pm',
	AUTHOR           => 'David E. Condon <dcondon@uidaho.edu>',
	LICENSE          => 'perl_5',
	MIN_PERL_VERSION => '5.010',
	# parent: lib/Alien/Bioinf.pm is loaded above, and perl 5.10.0 lacks it
	CONFIGURE_REQUIRES => { 'Alien::Build::MM' => '0', 'Alien::Build' => '2.00', 'ExtUtils::MakeMaker' => '6.52', 'parent' => '0' },
	BUILD_REQUIRES   => { 'Alien::CPython3' => '0' },
	PREREQ_PM        => {
		'Alien::Base' => '0', 'Alien::CPython3' => '0', 'Archive::Tar' => '0', 'Capture::Tiny' => '0',
		'Digest::MD5' => '0', 'HTTP::Tiny' => '0.014', 'IO::Socket::SSL' => '1.56', 'JSON::PP' => '0',
		'Net::SSLeay' => '1.49', 'parent' => '0',
	},
	TEST_REQUIRES    => { 'Test::More' => '0.88', 'Test::Alien' => '0' },
));

sub MY::postamble { $abmm->mm_postamble(@_) }
sub MY::install { $abmm->mm_install(@_) }
