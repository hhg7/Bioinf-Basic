require 5.012;
use strict;
use warnings FATAL => 'all';
use Test::More;
use Test::Alien;
use File::Spec;
use File::Basename qw(dirname);
use Alien::Bioinf;

alien_ok 'Alien::Bioinf';
my $v = Alien::Bioinf->versions;
ok defined $v->{$_}, "$_ has a recorded version" foreach qw(clustalo blast python);
run_ok([Alien::Bioinf->clustalo, '--version'])->success->out_like(qr/^\Q$v->{clustalo}\E\s*$/);
run_ok([Alien::Bioinf->blast('blastp'), '-version'])->success->out_like(qr/blastp: \Q$v->{blast}\E\+/);
SKIP: {
	skip "PLINK 2 has no build for $^O on this CPU, and none was given in ALIEN_BIOINF_PLINK2", 2
		unless defined $v->{plink2} || defined Alien::Bioinf::_plink2_platform() || defined $ENV{ALIEN_BIOINF_PLINK2};
	like $v->{plink2} // '', qr/\A20[0-9]{6}\z/, 'plink2 has a recorded build date';
	run_ok([Alien::Bioinf->plink2, '--version'])->success->out_like(qr/^PLINK v2\S* .*\(\d{1,2} [A-Z][a-z]{2} \d{4}\)/);
}
# The four tools added on 2026-11-13, each of which is installed only where it
# has a build or can be built: one with no recorded version was left out, with
# a warning, and is skipped here. What each prints for its version was read
# from bcftools 1.25, Prodigal 2.6.3 (which writes it to STDERR), kallisto
# 0.52.0 and SPAdes 4.3.0 that day.
my %ver = map { $_ => quotemeta($v->{$_} // '') } qw(bcftools prodigal kallisto);
foreach my $t (['bcftools', ['--version'], qr/^bcftools $ver{bcftools}$/m, 'out'],
	['prodigal', ['-v'], qr/^Prodigal V$ver{prodigal}:/m, 'err'],
	['kallisto', ['version'], qr/^kallisto, version $ver{kallisto}$/m, 'out']) {
	my ($tool, $args, $says, $stream) = @{ $t };
	SKIP: {
		skip "$tool was not installed on $^O", 2 unless defined $v->{$tool};
		my $run = run_ok([Alien::Bioinf->$tool, @{ $args }])->success;
		$stream eq 'err' ? $run->err_like($says) : $run->out_like($says);
	}
}
SKIP: {
	skip "SPAdes was not installed on $^O", 4 unless defined $v->{spades};
	# -E and -s, not -I, which would keep spades_init.py out of sys.path
	run_ok([Alien::Bioinf->python, '-E', '-s', Alien::Bioinf->spades, '--version'])->success
		->out_like(qr/^SPAdes genome assembler v\Q$v->{spades}\E$/m);
	my $bin = dirname(Alien::Bioinf->spades);
	ok -x File::Spec->catfile($bin, 'spades-core'), 'spades-core is beside spades.py';
	ok !-e File::Spec->catfile($bin, 'pathracer'), 'and pathracer, which spades.py does not run to assemble reads, is not';
}
is_deeply [Alien::Bioinf::_github_tag_of('https://github.com/hyattpd/Prodigal/releases/tag/v2.6.3')], ['2.6.3', 'v2.6.3'],
	'a release\'s version is its tag less the "v"';
is_deeply [Alien::Bioinf::_github_tag_of('https://github.com/samtools/bcftools/releases/tag/1.25')], ['1.25', '1.25'], 'and a tag with no "v" is the version';
is_deeply [Alien::Bioinf::_github_tag_of('https://github.com/samtools/bcftools/releases')], [], 'and a URL that names no release gives nothing';
# the build date of PLINK v2.0.0-a.7.13LM, downloaded on 2026-10-17
is Alien::Bioinf::_plink2_date("PLINK v2.0.0-a.7.13LM 64-bit Intel (15 Oct 2026)\n"), '20261015', '_plink2_date reads plink2 --version';
is Alien::Bioinf::_plink2_date("PLINK v2.0.0-a.6.9LM 64-bit Intel (5 Jan 2025)\n"), '20250105', 'and pads a one-digit day';
is Alien::Bioinf::_plink2_date('vcf_subset v1.1.0 64-bit (16 May 2026)'), undef, 'and gives undef for anything else';
run_ok([Alien::Bioinf->python, '-c', 'import numpy, matplotlib, adjustText, Bio, pydeseq2; import platform; print(platform.python_version())'])
	->success->out_like(qr/^\Q$v->{python}\E\s*$/);
# Up to 0.01 the venv was checked with "import biopython", which always fails
like Alien::Bioinf::_py_imports(), qr/\Aimport numpy,matplotlib,adjustText,Bio,pydeseq2\z/, 'the venv check imports Biopython as Bio, not as biopython';
# Up to 0.01 the venv saw the base Python's packages, and broke when they changed
my $cfg = File::Spec->catfile(dirname(dirname(Alien::Bioinf->python)), 'pyvenv.cfg');
open my $fh, '<', $cfg or die "$cfg: $!";
my $text = do { local $/; <$fh> };
like $text, qr/^home\s*=/m, 'pyvenv.cfg is read';
like $text, qr/^include-system-site-packages\s*=\s*false\s*$/m, "and the venv doesn't see the base Python's packages";
opendir my $dh, dirname(Alien::Bioinf->python) or die $!;
my @bin = sort grep { -f File::Spec->catfile(dirname(Alien::Bioinf->python), $_) } readdir $dh;
ok scalar(grep { /\Apython(?:w?\.exe)?\z/i } @bin), "the venv's interpreter is there";
is_deeply [grep { !/\Apython(?:w?\.exe)?\z/i } @bin], [], 'and nothing else, such as a pip that names the staging directory';
{
	local $ENV{ALIEN_INSTALL_NETWORK} = 0;
	ok !eval { Alien::Bioinf::_get('https://www.python.org/'); 1 }, 'ALIEN_INSTALL_NETWORK=0 stops a download';
	like $@, qr/ALIEN_INSTALL_NETWORK is 0, so nothing may be downloaded/, 'and says why';
	is_deeply [Alien::Bioinf::_pip_index()], ['--no-index'], 'and keeps pip off the network';
}
{
	local $ENV{ALIEN_INSTALL_NETWORK};
	is_deeply [Alien::Bioinf::_pip_index()], [], 'pip may use the network by default';
}
is Alien::Bioinf::_vcmp('2.17.0', '2.9.1'), 1, '_vcmp compares numerically, not as strings';
is Alien::Bioinf::_vcmp('1.2', '1.2.0'), 0, '_vcmp pads the shorter version with zeros';
done_testing;
