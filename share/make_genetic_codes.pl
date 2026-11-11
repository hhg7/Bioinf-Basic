#!/usr/bin/env perl
# Writes share/genetic_codes.tsv, the NCBI genetic codes that translate(),
# translate_six_frames() and find_orfs() read. Run it from the distribution's
# root:
#
#	perl share/make_genetic_codes.pl > share/genetic_codes.tsv
#
# or, to read a copy of gc.prt that is already on disk rather than fetch it:
#
#	perl share/make_genetic_codes.pl gc.prt > share/genetic_codes.tsv
#
# Cn::Bioinf never runs this, and [PruneFiles] keeps it out of the tarball:
# the committed .tsv is what ships and what is used. Re-run it only to refresh
# the codes against NCBI, and read the diff before committing, since a changed
# letter silently changes every translation made with that code.
#
# One source: https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt, NCBI's own
# table of genetic codes, which is what https://www.ncbi.nlm.nih.gov/Taxonomy/
# Utils/wprintgc.cgi displays. The committed file was written from version 4.6
# of it, fetched 2026-11-11. Each code there is a block of
#
#	name "Standard" ,
#	name "SGC0" ,
#	id 1 ,
#	ncbieaa  "FFLLSSSSYY**CC*WLLLL...",
#	sncbieaa "---M------**--*----M...",
#
# where ncbieaa gives the amino acid of each of the 64 codons, with * for a
# stop, and sncbieaa has M at each codon that can start a protein and * at each
# that can end one. The codons run TTT, TTC, TTA, TTG, TCT .. GGG: the first
# base changes slowest and each base goes through T, C, A, G, which the Base1,
# Base2 and Base3 comment lines of every block spell out and this script checks.
# A codon with an amino acid in ncbieaa and * in sncbieaa is one that the
# organism reads as either (TGA in code 27, for one).
#
# The output has one line for each code: its id, its names joined with "; ",
# ncbieaa and sncbieaa.
use 5.044;
use warnings FATAL => 'all';
use autodie ':default';

my $url = 'https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt';
my $text;
if (@ARGV) {
	open my $fh, '<', $ARGV[0];
	local $/;
	$text = <$fh>;
	close $fh;
} else {
	require HTTP::Tiny;
	my $r = HTTP::Tiny->new(timeout => 60)->get($url);
	die "couldn't fetch $url: $r->{status} $r->{reason}\n" unless $r->{success};
	$text = $r->{content};
}
my ($version) = $text =~ /^--\s+Version\s+(\S+)/m or die "gc.prt gives no version\n";

# the three comment lines that give each codon's bases, which must be these for
# the 64 letters to mean what the module takes them to mean
my @bases = (
	'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG',
	'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG',
	'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG',
);
my ($body) = $text =~ /Genetic-code-table ::= \{(.*)\}\s*\z/s or die "gc.prt has no Genetic-code-table\n";
my @rows;
foreach my $block ($body =~ /\{(.*?)\n \}/sg) {
	my @names = $block =~ /name "([^"]*)"/g;
	s/\s+/ /g for @names; # a long name is wrapped over two lines
	my ($id) = $block =~ /^\s*id (\d+)\s*,/m or die "a code has no id:\n$block\n";
	my ($aas) = $block =~ /ncbieaa\s+"([^"]{64})"/ or die "code $id has no ncbieaa of 64 letters\n";
	my ($starts) = $block =~ /sncbieaa\s+"([^"]{64})"/ or die "code $id has no sncbieaa of 64 letters\n";
	foreach my $i (1 .. 3) {
		my ($b) = $block =~ /-- Base$i\s+(\S+)/ or die "code $id has no Base$i line\n";
		die "code $id has its codons in another order: Base$i is $b\n" unless $b eq $bases[$i - 1];
	}
	die "code $id has no name\n" unless @names;
	die "code $id has a letter in ncbieaa that is not an amino acid or *: $aas\n" unless $aas =~ /\A[A-Z*]+\z/;
	die "code $id has a letter in sncbieaa that is not M, * or -: $starts\n" unless $starts =~ /\A[M*-]+\z/;
	# every stop of ncbieaa is a stop in sncbieaa too; the reverse does not hold
	foreach my $c (0 .. 63) {
		die "code $id: codon $c is * in ncbieaa and not in sncbieaa\n"
			if substr($aas, $c, 1) eq '*' && substr($starts, $c, 1) ne '*';
	}
	push @rows, [$id, join('; ', @names), $aas, $starts];
}
die "found only " . @rows . " codes\n" if @rows < 27; # version 4.6 has 27
say "## NCBI genetic codes: gc.prt version $version, written by share/make_genetic_codes.pl";
say join "\t", '#id', 'names', 'aas', 'starts';
say join "\t", @{ $_ } for sort { $a->[0] <=> $b->[0] } @rows;
