# 5.012, not the 5.010 Cn::Bioinf holds to: Alien::CPython3 needs
# File::Copy::Recursive, whose tests need Test::Deep, and Test::Deep 1.205
# requires 5.12 (found on perl-5.10.1, 2026-10-21). This comes first, ahead of
# "use Alien::Build::MM", so that a 5.10 perl stops on its version, which a
# CPAN smoker reports as NA, rather than on a missing module.
use 5.012;
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
	MIN_PERL_VERSION => '5.012',
	# parent: lib/Alien/Bioinf.pm is loaded above; it is core in every perl
	# from 5.10.1, and listed so that the list is the whole truth
	CONFIGURE_REQUIRES => { 'Alien::Build::MM' => '0', 'Alien::Build' => '2.00', 'ExtUtils::MakeMaker' => '6.52', 'parent' => '0' },
	BUILD_REQUIRES   => { 'Alien::CPython3' => '0' },
	# HTTP::Tiny 0.018 added the verify_SSL option that _http passes (its
	# Changes, 2012-04-18), and the 0.017 core in perl 5.16.3 ignores it and got
	# "599 SSL connection failed" from api.github.com, where 0.018 and 0.019 got
	# 200 (checked on 5.16.3-thr-ld, 2026-10-21). IO::Uncompress::Unzip 2.030
	# "Fixed problem with nextStream not returning 0 when there is no next
	# stream" (IO-Compress Changes, 2010-07-22): with 2.020 and 2.024, the core
	# copies in perl 5.10.1 and 5.12.5, _plink2_unpack's nextStream loop never
	# ended on the 3-file plink2 zip, and the install ran 30 minutes until it was
	# killed. 2.030 stopped after the 3 files on 5.12.5 (2026-10-21).
	PREREQ_PM        => {
		'Alien::Base' => '0', 'Alien::CPython3' => '0', 'Archive::Tar' => '0', 'Capture::Tiny' => '0',
		'Digest::MD5' => '0', 'HTTP::Tiny' => '0.018', 'IO::Uncompress::Unzip' => '2.030', 'IO::Socket::SSL' => '1.56', 'JSON::PP' => '0',
		'Net::SSLeay' => '1.49', 'parent' => '0',
	},
	TEST_REQUIRES    => { 'Test::More' => '0.88', 'Test::Alien' => '0' },
));

sub MY::postamble { $abmm->mm_postamble(@_) }
sub MY::install { $abmm->mm_install(@_) }
