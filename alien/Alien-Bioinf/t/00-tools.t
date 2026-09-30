require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use Test::Alien;
use Alien::Bioinf;

alien_ok 'Alien::Bioinf';
my $v = Alien::Bioinf->versions;
ok defined $v->{$_}, "$_ has a recorded version" foreach qw(clustalo blast python);
run_ok([Alien::Bioinf->clustalo, '--version'])->success->out_like(qr/^\Q$v->{clustalo}\E\s*$/);
run_ok([Alien::Bioinf->blast('blastp'), '-version'])->success->out_like(qr/blastp: \Q$v->{blast}\E\+/);
run_ok([Alien::Bioinf->python, '-c', 'import numpy, matplotlib, adjustText; import platform; print(platform.python_version())'])
	->success->out_like(qr/^\Q$v->{python}\E\s*$/);
is Alien::Bioinf::_vcmp('2.17.0', '2.9.1'), 1, '_vcmp compares numerically, not as strings';
is Alien::Bioinf::_vcmp('1.2', '1.2.0'), 0, '_vcmp pads the shorter version with zeros';
done_testing;
