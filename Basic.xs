/*FASTA reading and writing for Bioinf::Basic.

fasta2hash() and hash2fasta_file() used to be pure perl in ~/Scripts/bioinf.pm:
one readline, chomp and string append per line in, and one substr-sized SV per
80-column line out. Here a file is read in large blocks and every line of a
record is appended straight onto that record's SV, so no per-line SV is ever
made, and each finished sequence is shrunk to its length, which is the RAM a
proteome or genome actually needs rather than whatever the last buffer growth
left. On the file described at FASTA_BLOCK below, bioinf.pm read it in 3.83 s
with a peak RSS of 1,098 MB, and this in 0.80 s and 967 MB, 900 MB of which is
the sequences themselves; writing it back took 44.44 s and takes 0.97 s, and
the two files are byte-identical. Dropping the whitespace in sequences, which
began on 2026-10-17, costs about 15%: a made-up file of the same shape (938 MB,
400,000 records of 2,280 residues in 60-column lines, in the page cache) took
0.85 s to read before and 0.98 s after, best of five timed inside perl. The
perl side (lib/Bioinf/Basic.pm) opens and closes the files.

Every croak() and warn() here starts with C99's __func__, so a message names
the C function that raised it (read_fasta or write_fasta) and stays right if
one is renamed.*/
#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include <string.h>
#include <stdint.h>
/*`restrict` is spelled bare, which is C99 syntax.

It is carried only by read_fasta()'s block buffer, which this file allocates
itself and reaches through nothing but that pointer and the ones derived from
it. Every other pointer here is perl's: an SV, AV or HV, or the PV inside one,
which the perl API may alias behind the code's back -- in write_fasta() a
defline and its sequence can be COW copies of the one buffer -- so the
qualifier would be a false promise. perl is built with -fno-strict-aliasing,
where `restrict` is the only disambiguation there is, but it earns nothing even
on the buffer: the loop over it stores only through memchr() and sv_catpvn(),
which are opaque calls the compiler must assume write anywhere. The generated
code for this file was byte-identical with and without it (gcc -O2 -std=gnu99,
2026-09-30), so the one annotation records the contract rather than buying
speed.

Makefile.PL probes for a C99 flag but cannot always find one. MSVC accepts the
keyword only under /std:c11 or later and otherwise spells it __restrict, older
gcc spells it __restrict__, and a strict C89 compiler has no equivalent. Map it
where a spelling exists and define it away where none does -- the same block as
~/Scripts/stats/LikeR.xs. It sits below every #include so that it cannot
rewrite a parameter named `restrict` inside a system header.*/
#if !defined(__cplusplus) && !defined(restrict)
#  if defined(_MSC_VER)
#    define restrict __restrict
#  elif !defined(__STDC_VERSION__) || __STDC_VERSION__ < 199901L
#    if defined(__GNUC__)
#      define restrict __restrict__
#    else
#      define restrict /*no spelling available; annotation only*/
#    endif
#  endif
#endif
/*PerlIO_read() request size. Reading a 928 MB FASTA of 400,000 records in
60-column lines, on this machine on 2026-09-29, took 0.84 s with 16 KiB blocks,
0.80 s with 64 KiB, 0.81 s with 256 KiB and 0.82 s with 1 MiB (best of three
each, file in the page cache), so nothing is gained past 64 KiB. Build with
DEFINE=-DFASTA_BLOCK=<bytes> to try another.*/
#ifndef FASTA_BLOCK
#define FASTA_BLOCK 65536
#endif

/*Whether c is whitespace as Python has it, which is what Biopython strips: ' '
or one of '\t' (9), '\v' (11), '\f' (12) and '\r' (13), which is the range
9..13 less the '\n' (10) that ends every line, and so is never tested. perl's
isSPACE() is not used, since that took in '\v' only from 5.18.*/
#define FASTA_SPACE(c) ((c) == ' ' || (unsigned char)((c) - '\t') < 5)

/*Whether p..e holds any FASTA_SPACE() byte, 8 bytes at a time.

With a byte-at-a-time loop, the 938 MB file described at the top of this file
took 1.46 s to read rather than 0.85 s (best of five each, 2026-10-17): gcc
-O2 does not vectorise a loop of unknown length. This takes 8 bytes at once: a space is
a zero byte of w ^ 0x2020...20, found by Alan Mycroft's haszero(), and a byte
in 9..13 is found by hasbetween(w, 8, 14), both from Sean Eron Anderson's "Bit
Twiddling Hacks" (graphics.stanford.edu/~seander/bithacks.html, sections
ZeroInWord and HasBetweenInWord). Both are exact, and work byte by byte with no
carry between bytes, so the byte order does not matter.*/
static bool any_space(const char *p, const char *e) {
	const uint64_t ones = UINT64_C(0x0101010101010101), low7 = ones * 0x7f, high = ones * 0x80;
	for (; e - p >= 8; p += 8) {
		uint64_t w;
		memcpy(&w, p, 8); // not a cast to uint64_t *, since SPARC and some ARM fault on an unaligned load
		const uint64_t sp = w ^ (ones * ' ');
		if ((sp - ones) & ~sp & high) return TRUE;
		if ((ones * (127 + 14) - (w & low7)) & ~w & ((w & low7) + ones * (127 - 8)) & high) return TRUE;
	}
	for (; p < e; p++) {
		if (FASTA_SPACE(*p)) return TRUE;
	}
	return FALSE;
}

/*Drop the whitespace at the end of a defline, the '\r' of a CRLF file
included, as Biopython does. It is chopped when the line ends rather than as it
is read, because a block can end inside it.*/
static void chop_space(SV *sv) {
	STRLEN n = SvCUR(sv);
	while (n && FASTA_SPACE(SvPVX(sv)[n - 1])) n--;
	SvCUR_set(sv, n);
}

/*Append the residues in p..e to seq, leaving out the whitespace. Biopython
drops the spaces in a sequence line and the whitespace at its end; a tab inside
a line is dropped here too. A line with none, which is nearly every line, is
appended whole; otherwise each run between whitespace is one sv_catpvn().*/
static void cat_residues(pTHX_ SV *seq, const char *p, const char *e) {
	if (!any_space(p, e)) {
		sv_catpvn(seq, p, (STRLEN)(e - p));
		return;
	}
	while (p < e) {
		const char *run = p;
		while (p < e && !FASTA_SPACE(*p)) p++;
		if (p > run) sv_catpvn(seq, run, (STRLEN)(p - run));
		while (p < e && FASTA_SPACE(*p)) p++;
	}
}

/*Whether p..e holds nothing but whitespace.*/
static bool all_space(const char *p, const char *e) {
	for (; p < e; p++) {
		if (!FASTA_SPACE(*p)) return FALSE;
	}
	return TRUE;
}

/*Every record of fh as a hash of defline => sequence, or with key given only
that record's sequence (undef if it is absent), reading no further than the
defline after it. name is the file's name, for messages. order, when not NULL,
gets each defline in file order.

A repeated defline is warned about and its sequences concatenated, which is
what bioinf.pm did: it is the file that is wrong, and nothing here can say
which of the two records was meant. With key given, only a repeat of key is
looked for. Finding the others would mean keeping every defline read on the way
to key: a hash of 400,000 keys, for one record of a proteome.*/
static SV *read_fasta(pTHX_ PerlIO *fh, SV *key, const char *name, AV *order) {
	HV *out = (HV*)sv_2mortal((SV*)newHV());
	SV *hdr = sv_2mortal(newSVpvs(""));
	SV *seq = NULL; // record being filled; NULL before the first defline and while skipping
	bool have_def = FALSE, in_hdr = FALSE, bol = TRUE, done = FALSE, fed = FALSE; // fed: the '\n' that ends an unterminated last line has been supplied
	UV line = 0; // 1-based, as $. would count it
	char *restrict buf;
	Newx(buf, FASTA_BLOCK, char);
	SAVEFREEPV(buf); // freed at the caller's LEAVE, croak or not
	while (!done) {
		SSize_t n = PerlIO_read(fh, buf, FASTA_BLOCK);
		if (n < 0 || PerlIO_error(fh)) croak("%s: error reading %s: %s", __func__, name, Strerror(errno));
		if (n == 0) {
			if (bol || fed) break;
			buf[0] = '\n';
			n = 1;
			fed = TRUE;
		}
		const char *p = buf, *end = buf + n;
		while (p < end) {
			if (bol) {
				line++;
				bol = FALSE;
				if (*p == '>') {
					in_hdr = TRUE;
					SvCUR_set(hdr, 0);
					p++;
					continue;
				}
			}
			const char *nl = (const char*)memchr(p, '\n', (size_t)(end - p));
			const char *e = nl ? nl : end;
			if (in_hdr) sv_catpvn(hdr, p, (STRLEN)(e - p));
			else if (seq) cat_residues(aTHX_ seq, p, e);
			else if (!have_def && !all_space(p, e))
				croak("%s: %s line %" UVuf " is sequence, but no defline (\">...\") has come before it", __func__, name, line);
			if (!nl) break;
			p = nl + 1;
			bol = TRUE;
			if (!in_hdr) continue;
			in_hdr = FALSE;
			chop_space(hdr);
			if (!SvCUR(hdr)) croak("%s: %s line %" UVuf " is a defline with no name", __func__, name, line);
			if (seq) SvPV_shrink_to_cur(seq);
			bool is_key = key && sv_eq(hdr, key);
			if (key && !is_key && hv_exists_ent(out, key, 0)) {
				done = TRUE; // the wanted record is complete
				break;
			}
			have_def = TRUE;
			if (key && !is_key) {
				seq = NULL;
				continue;
			}
			if (hv_exists_ent(out, hdr, 0))
				warn("%s: \"%s\" appears more than once in %s (line %" UVuf "); its sequences will be concatenated", __func__, SvPV_nolen(hdr), name, line);
			else if (order) av_push(order, newSVsv(hdr));
			seq = HeVAL(hv_fetch_ent(out, hdr, 1, 0));
			if (!SvOK(seq)) sv_setpvs(seq, "");
		}
	}
	if (seq) SvPV_shrink_to_cur(seq);
	if (!key) return newRV_inc((SV*)out);
	HE *he = hv_fetch_ent(out, key, 0, 0);
	return he ? SvREFCNT_inc_simple_NN(HeVAL(he)) : newSV(0); // not a copy: out is mortal, so this becomes the only reference
}

/*The bytes of sv for a FASTA file, as print to a :raw handle writes them: a
string perl holds as UTF-8 is written as the characters it holds, so that how
perl happened to store it cannot change the file. Up to 2026-10-17 it was
written as its UTF-8, and "caf\xe9" came out as one byte or two. scratch holds
that downgraded copy. A character above 255 has no byte, and gets NULL, for
the caller to die naming the string; print would write UTF-8 and warn.*/
static const char *fasta_bytes(pTHX_ SV *sv, SV *scratch, STRLEN *len) {
	if (SvUTF8(sv)) {
		sv_setsv(scratch, sv);
		if (!sv_utf8_downgrade(scratch, TRUE)) return NULL;
		sv = scratch;
	}
	return SvPV_const(sv, *len);
}

/*Each key of h named in order, as ">key" and then its sequence broken into
lines of width characters; width 0 puts each sequence on one line.*/
static void write_fasta(pTHX_ PerlIO *fh, HV *h, AV *order, STRLEN width) {
	SV *kbuf = sv_newmortal(), *sbuf = sv_newmortal(); // fasta_bytes()'s scratch copies, of the key and of its sequence
	for (SSize_t i = 0; i <= av_len(order); i++) {
		SV **k = av_fetch(order, i, 0);
		if (!k || !SvOK(*k)) croak("%s: element %" IVdf " of the order is undefined", __func__, (IV)i);
		HE *he = hv_fetch_ent(h, *k, 0, 0);
		STRLEN kl, sl;
		// kp and sp not restrict: a defline and its sequence may share one COW buffer
		const char *kp = fasta_bytes(aTHX_ *k, kbuf, &kl);
		if (!kp) croak("%s: \"%s\" has a character above 255, so it has to be encoded (with Encode::encode('UTF-8', ...), say) before it is written", __func__, SvPV_nolen(*k));
		if (!he || !SvOK(HeVAL(he))) croak("%s: \"%s\" has no sequence in the hash", __func__, kp);
		const char *sp = fasta_bytes(aTHX_ HeVAL(he), sbuf, &sl);
		if (!sp) croak("%s: the sequence of \"%s\" has a character above 255, so it has to be encoded before it is written", __func__, kp);
		PerlIO_putc(fh, '>');
		PerlIO_write(fh, kp, kl);
		PerlIO_putc(fh, '\n');
		STRLEN w = width ? width : sl;
		for (STRLEN o = 0; o < sl; o += w) {
			PerlIO_write(fh, sp + o, sl - o < w ? sl - o : w);
			PerlIO_putc(fh, '\n');
		}
	}
	if (PerlIO_error(fh)) croak("%s: error writing FASTA: %s", __func__, Strerror(errno));
}

MODULE = Bioinf::Basic	PACKAGE = Bioinf::Basic
PROTOTYPES: DISABLE

SV *
_read_fasta(fh, key, name, order = NULL)
	PerlIO *fh
	SV *key
	const char *name
	AV *order
CODE:
	RETVAL = read_fasta(aTHX_ fh, SvOK(key) ? key : NULL, name, order);
OUTPUT:
	RETVAL

void
_write_fasta(fh, h, order, width)
	PerlIO *fh
	HV *h
	AV *order
	UV width
CODE:
	write_fasta(aTHX_ fh, h, order, (STRLEN)width);
