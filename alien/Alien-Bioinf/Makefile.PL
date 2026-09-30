use strict;
use warnings;
use ExtUtils::MakeMaker;
use Alien::Build::MM;

my $abmm = Alien::Build::MM->new;
WriteMakefile($abmm->mm_args(
	NAME             => 'Alien::Bioinf',
	DISTNAME         => 'Alien-Bioinf',
	VERSION_FROM     => 'lib/Alien/Bioinf.pm',
	ABSTRACT_FROM    => 'lib/Alien/Bioinf.pm',
	AUTHOR           => 'David E. Condon <dcondon@uidaho.edu>',
	LICENSE          => 'perl_5',
	MIN_PERL_VERSION => '5.010',
	CONFIGURE_REQUIRES => { 'Alien::Build::MM' => '0', 'Alien::Build' => '2.00', 'ExtUtils::MakeMaker' => '6.52' },
	BUILD_REQUIRES   => { 'Alien::CPython3' => '0' },
	PREREQ_PM        => {
		'Alien::Base' => '0', 'Alien::CPython3' => '0', 'Archive::Tar' => '0', 'Capture::Tiny' => '0',
		'Digest::MD5' => '0', 'HTTP::Tiny' => '0.014', 'IO::Socket::SSL' => '1.56', 'JSON::PP' => '0',
		'Net::SSLeay' => '1.49',
	},
	TEST_REQUIRES    => { 'Test::More' => '0.88', 'Test::Alien' => '0' },
));

sub MY::postamble { $abmm->mm_postamble(@_) }
sub MY::install { $abmm->mm_install(@_) }
