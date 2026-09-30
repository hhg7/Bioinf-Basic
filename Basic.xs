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
the two files are byte-identical. The perl side (lib/Bioinf/Basic.pm) opens and
closes the files.*/
#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include <string.h>
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

/*Drop a trailing '\r' left by a CRLF file. It is chopped when the line ends
rather than when it is read, because a block can end between the '\r' and
its '\n'.*/
static void chop_cr(SV *sv) {
	if (SvCUR(sv) && SvPVX(sv)[SvCUR(sv) - 1] == '\r') SvCUR_set(sv, SvCUR(sv) - 1);
}

/*Every record of fh as a hash of defline => sequence, or with key given only
that record's sequence (undef if it is absent), reading no further than the
defline after it. name is the file's name, for messages. order, when not NULL,
gets each defline in file order.

A repeated defline is warned about and its sequences concatenated, which is
what bioinf.pm did: it is the file that is wrong, and nothing here can say
which of the two records was meant.*/
static SV *read_fasta(pTHX_ PerlIO *fh, SV *key, const char *name, AV *order) {
	HV *out = (HV*)sv_2mortal((SV*)newHV());
	HV *seen = key ? (HV*)sv_2mortal((SV*)newHV()) : out; // other deflines met so far; with every record kept, out already is that
	SV *hdr = sv_2mortal(newSVpvs(""));
	SV *seq = NULL; // record being filled; NULL before the first defline and while skipping
	bool have_def = FALSE, in_hdr = FALSE, bol = TRUE, done = FALSE, fed = FALSE; // fed: the '\n' that ends an unterminated last line has been supplied
	UV line = 0; // 1-based, as $. would count it
	char *restrict buf;
	Newx(buf, FASTA_BLOCK, char);
	SAVEFREEPV(buf); // freed at the caller's LEAVE, croak or not
	while (!done) {
		SSize_t n = PerlIO_read(fh, buf, FASTA_BLOCK);
		if (n < 0 || PerlIO_error(fh)) croak("error reading %s: %s", name, Strerror(errno));
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
			else if (seq) sv_catpvn(seq, p, (STRLEN)(e - p));
			else if (!have_def && e > p && !(e - p == 1 && *p == '\r'))
				croak("%s line %" UVuf " is sequence, but no defline (\">...\") has come before it", name, line);
			if (!nl) break;
			p = nl + 1;
			bol = TRUE;
			if (!in_hdr) {
				if (seq) chop_cr(seq);
				continue;
			}
			in_hdr = FALSE;
			chop_cr(hdr);
			if (!SvCUR(hdr)) croak("%s line %" UVuf " is a defline with no name", name, line);
			if (seq) SvPV_shrink_to_cur(seq);
			bool is_key = key && sv_eq(hdr, key);
			if (key && !is_key && hv_exists_ent(out, key, 0)) {
				done = TRUE; // the wanted record is complete
				break;
			}
			have_def = TRUE;
			if (hv_exists_ent(is_key ? out : seen, hdr, 0))
				warn("\"%s\" appears more than once in %s (line %" UVuf "); its sequences will be concatenated", SvPV_nolen(hdr), name, line);
			else if (order) av_push(order, newSVsv(hdr));
			if (key && !is_key) {
				(void)hv_store_ent(seen, hdr, newSV(0), 0);
				seq = NULL;
				continue;
			}
			seq = HeVAL(hv_fetch_ent(out, hdr, 1, 0));
			if (!SvOK(seq)) sv_setpvs(seq, "");
		}
	}
	if (seq) SvPV_shrink_to_cur(seq);
	if (!key) return newRV_inc((SV*)out);
	HE *he = hv_fetch_ent(out, key, 0, 0);
	return he ? SvREFCNT_inc_simple_NN(HeVAL(he)) : newSV(0); // not a copy: out is mortal, so this becomes the only reference
}

/*Each key of h named in order, as ">key" and then its sequence broken into
lines of width characters; width 0 puts each sequence on one line.*/
static void write_fasta(pTHX_ PerlIO *fh, HV *h, AV *order, STRLEN width) {
	for (SSize_t i = 0; i <= av_len(order); i++) {
		SV **k = av_fetch(order, i, 0);
		if (!k || !SvOK(*k)) croak("element %" IVdf " of the order is undefined", (IV)i);
		HE *he = hv_fetch_ent(h, *k, 0, 0);
		STRLEN kl, sl;
		// kp and sp not restrict: a defline and its sequence may share one COW buffer
		const char *kp = SvPV_const(*k, kl);
		if (!he || !SvOK(HeVAL(he))) croak("\"%s\" has no sequence in the hash", kp);
		const char *sp = SvPV_const(HeVAL(he), sl);
		PerlIO_putc(fh, '>');
		PerlIO_write(fh, kp, kl);
		PerlIO_putc(fh, '\n');
		STRLEN w = width ? width : sl;
		for (STRLEN o = 0; o < sl; o += w) {
			PerlIO_write(fh, sp + o, sl - o < w ? sl - o : w);
			PerlIO_putc(fh, '\n');
		}
	}
	if (PerlIO_error(fh)) croak("error writing FASTA: %s", Strerror(errno));
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
