/*FASTA reading and writing for Cn::Bioinf, and the genotype counts behind its
sex_check and parent.

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
perl side (lib/Cn/Bioinf.pm) opens and closes the files.

Every croak() and warn() here starts with C99's __func__, so a message names
the C function that raised it (read_fasta, write_fasta, or one of the vcf_,
sex_ and parent_ ones) and stays right if one is renamed.*/
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

/*Reading a VCF, for sex_check() and parent().

The file is read by perl, which knows gzip and bzip2 through IO::Uncompress
(linking zlib and libbz2 here would be a build dependency on every platform),
and handed over a block at a time. vcf_read() cuts the blocks into lines and
gives each record to the reader's on_record, keeping nothing itself but the
one line a block can end inside; what each caller keeps is its own.*/

typedef struct vcf_reader vcf_reader;
struct vcf_reader {
	const char *name; // the file's, for messages
	UV line;          // 1-based, as $. would count it
	UV min_depth;     // a genotype whose DP is below this is not counted
	size_t n;         // samples; 0 until the #CHROM line
	bool header;      // the #CHROM line has been read
	HV *samples;      // sample => its 0-based column among the samples
	AV *names;        // the samples in file order
	void (*on_header)(pTHX_ vcf_reader *r);                               // once the #CHROM line has been read
	void (*on_record)(pTHX_ vcf_reader *r, const char *p, const char *e); // each record, its '\n' gone
};

/*The end of the tab-separated field that starts at p, before e.*/
static const char *vcf_field_end(const char *p, const char *e) {
	const char *t = (const char*)memchr(p, '\t', (size_t)(e - p));
	return t ? t : e;
}

/*Sub-field k (0-based) of the colon-separated p..e, in *s..*f, or FALSE when
there are fewer.*/
static bool vcf_subfield(const char *p, const char *e, size_t k, const char **s, const char **f) {
	for (; k; k--) {
		p = (const char*)memchr(p, ':', (size_t)(e - p));
		if (!p) return FALSE;
		p++;
	}
	*s = p;
	const char *c = (const char*)memchr(p, ':', (size_t)(e - p));
	*f = c ? c : e;
	return TRUE;
}

/*The chromosome that CHROM p..f names, numbered as PLINK numbers them: 1 to
22, 23 for X and 24 for Y, and 0 for any other (MT, an unplaced contig). A
"chr" prefix, in any case, is dropped.*/
static unsigned short int vcf_chrom(const char *p, const char *f) {
	if (f - p > 3 && (p[0] == 'c' || p[0] == 'C') && (p[1] == 'h' || p[1] == 'H') && (p[2] == 'r' || p[2] == 'R')) p += 3;
	if (f - p == 1 && (*p == 'X' || *p == 'x')) return 23;
	if (f - p == 1 && (*p == 'Y' || *p == 'y')) return 24;
	if (f - p < 1 || f - p > 2 || *p < '1' || *p > '9') return 0; // "01" is not chromosome 1
	unsigned short int c = (unsigned short int)(*p - '0');
	if (f - p == 2) {
		if (p[1] < '0' || p[1] > '9') return 0;
		c = (unsigned short int)(c * 10 + (p[1] - '0'));
	}
	return c <= 24 ? c : 0;
}

/*The index of base c, 0 to 3 for A, C, G and T in either case, or 4 for
anything else.*/
static unsigned short int vcf_base_index(char c) {
	switch (c) {
		case 'A': case 'a': return 0;
		case 'C': case 'c': return 1;
		case 'G': case 'g': return 2;
		case 'T': case 't': return 3;
		default: return 4;
	}
}

/*Whether p..e is one base, as a SNV's REF and ALT are. An ALT of ".", "*" or
"<NON_REF>", or with a comma, is not.*/
static bool vcf_base(const char *p, const char *e) {
	return e - p == 1 && vcf_base_index(*p) < 4;
}

/*Where each of the 9 fixed columns of the record p..e starts, in col, given
that CHROM ends at f; returns where FORMAT ends, so that col[i + 1] - 1 ends
column i.*/
static const char *vcf_columns(pTHX_ const vcf_reader *r, const char *p, const char *f, const char *e, const char **col) {
	col[0] = p;
	for (unsigned short int i = 1; i < 9; i++) {
		if (f >= e) croak("%s: %s line %" UVuf " has %u columns, and a record needs 9 and a sample", __func__, r->name, r->line, (unsigned)i);
		col[i] = f + 1;
		f = vcf_field_end(col[i], e);
	}
	return f;
}

/*POS, the field p..f.*/
static UV vcf_pos(pTHX_ const vcf_reader *r, const char *p, const char *f) {
	if (p == f) croak("%s: %s line %" UVuf " has an empty POS", __func__, r->name, r->line);
	UV pos = 0;
	for (const char *q = p; q < f; q++) {
		if (*q < '0' || *q > '9' || pos > (UV_MAX - 9) / 10) croak("%s: %s line %" UVuf " has POS \"%.*s\", which isn't a position", __func__, r->name, r->line, (int)(f - p), p);
		pos = pos * 10 + (UV)(*q - '0');
	}
	return pos;
}

#define VCF_ABSENT ((size_t)-1) // a key that FORMAT doesn't name

/*GT's and DP's places in FORMAT, the field p..e, in *gt and *dp, or
VCF_ABSENT.*/
static void vcf_format(const char *p, const char *e, size_t *gt, size_t *dp) {
	*gt = *dp = VCF_ABSENT;
	size_t k = 0;
	for (const char *s = p; s <= e; k++) {
		const char *c = (const char*)memchr(s, ':', (size_t)(e - s));
		const char *ce = c ? c : e;
		if (ce - s == 2 && s[0] == 'G' && s[1] == 'T') *gt = k;
		else if (ce - s == 2 && s[0] == 'D' && s[1] == 'P') *dp = k;
		s = ce + 1;
	}
}

/*What vcf_call() found.*/
enum {
	VCF_SHALLOW, // a DP below min_depth: the genotype is not counted at all
	VCF_NO_CALL, // missing, partly missing ("./1"), more than two alleles, or an allele the site lacks
	VCF_CALL     // one or two alleles, each REF or the one ALT
};

/*The call in the sample field p..e of a biallelic site, whose FORMAT has GT
at gt and DP at dp: one of the values above, and for VCF_CALL the number of
alleles in *alleles and how many of them are the ALT in *alt. A DP that is
absent, "." or not a whole number leaves the genotype counted.*/
static short int vcf_call(const vcf_reader *r, const char *p, const char *e, size_t gt, size_t dp, unsigned short int *alleles, unsigned short int *alt) {
	const char *s, *sf;
	if (dp != VCF_ABSENT && vcf_subfield(p, e, dp, &s, &sf) && s < sf && *s != '.') {
		UV depth = 0;
		bool number = TRUE;
		for (; s < sf; s++) {
			if (*s < '0' || *s > '9' || depth > (UV_MAX - 9) / 10) {
				number = FALSE;
				break;
			}
			depth = depth * 10 + (UV)(*s - '0');
		}
		if (number && depth < r->min_depth) return VCF_SHALLOW;
	}
	*alleles = *alt = 0;
	if (!vcf_subfield(p, e, gt, &s, &sf) || s == sf) return VCF_NO_CALL;
	while (s < sf) {
		if (++*alleles > 2 || *s == '.') return VCF_NO_CALL;
		if (*s == '1') ++*alt;
		else if (*s != '0') return VCF_NO_CALL; // no other allele exists at a biallelic site
		s++;
		if (s < sf && *s != '/' && *s != '|') return VCF_NO_CALL;
		if (s < sf) s++;
	}
	return VCF_CALL;
}

/*The end of sample i's field, which starts at p; dies when the record ends
before it.*/
static const char *vcf_sample(pTHX_ const vcf_reader *r, const char *p, const char *e, size_t i) {
	if (p > e) croak("%s: %s line %" UVuf " has %" UVuf " samples, and the #CHROM line names %" UVuf, __func__, r->name, r->line, (UV)i, (UV)r->n);
	return vcf_field_end(p, e);
}

/*Dies when anything is left of the record after its last sample, which would
start at p.*/
static void vcf_samples_end(pTHX_ const vcf_reader *r, const char *p, const char *e) {
	if (p <= e) croak("%s: %s line %" UVuf " has more samples than the %" UVuf " its #CHROM line names", __func__, r->name, r->line, (UV)r->n);
}

/*The #CHROM line p..e: each sample after FORMAT gets its column in samples
and its place in names, and then the reader's on_header is called.*/
static void vcf_header(pTHX_ vcf_reader *r, const char *p, const char *e) {
	if (r->header) croak("%s: %s line %" UVuf " is a second #CHROM line", __func__, r->name, r->line);
	r->header = TRUE;
	unsigned short int col = 0; // columns passed; the samples start after the 9th, FORMAT
	while (p <= e && col < 9) {
		p = vcf_field_end(p, e) + 1;
		col++;
	}
	if (p > e) croak("%s: %s has no samples: its #CHROM line (line %" UVuf ") ends at %s", __func__, r->name, r->line, col < 9 ? "INFO or before" : "FORMAT");
	while (p <= e) {
		const char *f = vcf_field_end(p, e);
		SV *name = newSVpvn(p, (STRLEN)(f - p));
		if (hv_exists_ent(r->samples, name, 0)) {
			// mortal on a line of its own: SvPV_nolen() is a macro that evaluates its
			// argument twice on perl 5.10.1 to 5.32.1, which then freed the name twice
			// ("Attempt to free unreferenced scalar"); 5.42.3 and 5.44.0 were unaffected
			sv_2mortal(name);
			const char *dup = SvPV_nolen(name);
			croak("%s: %s names sample \"%s\" twice in its #CHROM line", __func__, r->name, dup);
		}
		(void)hv_store_ent(r->samples, name, newSVuv((UV)r->n), 0);
		av_push(r->names, name);
		r->n++;
		p = f + 1;
	}
	if (r->on_header) r->on_header(aTHX_ r);
}

/*The line p..e, its '\n' gone.*/
static void vcf_line(pTHX_ vcf_reader *r, const char *p, const char *e) {
	r->line++;
	if (e > p && e[-1] == '\r') e--;
	if (p == e) return;
	if (*p != '#') {
		if (!r->header) croak("%s: %s line %" UVuf " is a record, but no #CHROM line has come before it", __func__, r->name, r->line);
		r->on_record(aTHX_ r, p, e);
	} else if (e - p >= 6 && memcmp(p, "#CHROM", 6) == 0) {
		vcf_header(aTHX_ r, p, e);
	}
}

/*The VCF r->name, which next, a perl sub, returns a block of at a time and
then undef. r's name, min_depth, on_header and on_record are the caller's to
set; the rest is set here.*/
static void vcf_read(pTHX_ SV *next, vcf_reader *r) {
	r->line = 0;
	r->n = 0;
	r->header = FALSE;
	r->samples = (HV*)sv_2mortal((SV*)newHV());
	r->names = (AV*)sv_2mortal((SV*)newAV());
	SV *carry = sv_2mortal(newSVpvs("")); // the start of a line that the last block ended inside
	/*The block, copied out of the scope that next's own temporaries are freed
	with, so that whatever on_header saves for freeing (with SAVEFREEPV, say)
	lands in the caller's scope rather than that one; from perl 5.20 the copy
	is copy-on-write, and moves no bytes.*/
	SV *block = sv_newmortal();
	for (;;) {
		dSP;
		ENTER;
		SAVETMPS;
		PUSHMARK(SP);
		PUTBACK;
		I32 got = call_sv(next, G_SCALAR);
		SPAGAIN;
		sv_setsv(block, got ? POPs : &PL_sv_undef);
		PUTBACK;
		FREETMPS;
		LEAVE;
		if (!SvOK(block)) break;
		STRLEN n;
		const char *b = SvPV_const(block, n), *end = b + n; // not restrict: perl's buffer
		while (b < end) {
			const char *nl = (const char*)memchr(b, '\n', (size_t)(end - b));
			if (!nl) {
				sv_catpvn(carry, b, (STRLEN)(end - b));
				break;
			}
			if (SvCUR(carry)) {
				sv_catpvn(carry, b, (STRLEN)(nl - b));
				vcf_line(aTHX_ r, SvPVX_const(carry), SvPVX_const(carry) + SvCUR(carry));
				SvCUR_set(carry, 0);
			} else {
				vcf_line(aTHX_ r, b, nl);
			}
			b = nl + 1;
		}
	}
	if (SvCUR(carry)) vcf_line(aTHX_ r, SvPVX_const(carry), SvPVX_const(carry) + SvCUR(carry));
	if (!r->header) croak("%s: %s has no #CHROM line, so it isn't a VCF", __func__, r->name);
}

/*The genotype counts behind sex_check().

The method is peddy 0.4.3's sex_check() (peddy/peddy.py): on X outside its
pseudo-autosomal regions, count each sample's heterozygous and homozygous-ALT
calls at biallelic sites, leaving out a genotype whose depth is below a
minimum. Its PAR bounds, X:10000-2781479 and X:154931044-156030895, cover both
GRCh37's and GRCh38's, so no build need be named. It departs from peddy in
three ways. Every site is read, where peddy stops at 10,000. Only SNVs are
counted, where peddy keeps an insertion (its test is len(REF) > 1 alone). And
peddy's cohort filters (a call rate of at least 0.5, an ALT frequency of at
least 0.01) are left out, since they mean nothing for a VCF of one genome,
which is the usual input here.

Non-PAR Y is counted too, for the caller to see, but decides nothing: on
2026-10-18 the exome mother, whose X is plainly female (2,592 het to 2,673
hom-ALT), had 530 called SNVs on non-PAR Y, and the male exome father had
659. Y's PARs are GRCh37's (Y:10001-2649520, Y:59034050-59363566) and
GRCh38's (Y:10001-2781479, Y:56887903-57217415), from the GRC's assembly
reports.

Nothing is kept but a fixed set of counters for each sample.*/

/*The counters kept for each sample, in this order in sex_tally's count; their
names in the hash returned to perl are SEX_NAME's.*/
enum { SEX_X_HET, SEX_X_HOM_ALT, SEX_X_HOM_REF, SEX_X_HAPLOID, SEX_X_MISSING, SEX_Y_CALLED, SEX_Y_ALT, SEX_FIELDS };
static const char *const SEX_NAME[SEX_FIELDS] = { "x_het", "x_hom_alt", "x_hom_ref", "x_haploid", "x_missing", "y_called", "y_alt" };

typedef struct {
	vcf_reader r;  // first, so that the vcf_reader * that sex_header() and sex_record() are given is this
	UV x_sites;    // non-PAR X SNVs read, for perl to say when there were none
	UV *count;     // r.n * SEX_FIELDS counters, sample by sample
} sex_tally;

/*Whether pos on chrom (23 for X, 24 for Y) lies in a pseudo-autosomal
region; see the comment above SEX_NAME for the bounds.*/
static bool sex_par(unsigned short int chrom, UV pos) {
	if (chrom == 23) return (pos >= 10000 && pos <= 2781479) || (pos >= 154931044 && pos <= 156030895);
	return (pos >= 10001 && pos <= 2781479) || (pos >= 56887903 && pos <= 57217415) || (pos >= 59034050 && pos <= 59363566);
}

/*Each sample's counters, once the #CHROM line has said how many samples
there are.*/
static void sex_header(pTHX_ vcf_reader *r) {
	sex_tally *t = (sex_tally*)r;
	Newxz(t->count, r->n * SEX_FIELDS, UV);
	SAVEFREEPV(t->count); // freed at the caller's LEAVE, croak or not; see vcf_read()'s block
}

/*The record p..e.*/
static void sex_record(pTHX_ vcf_reader *r, const char *p, const char *e) {
	sex_tally *t = (sex_tally*)r;
	const char *f = vcf_field_end(p, e);
	unsigned short int chrom = vcf_chrom(p, f); // only 23 (X) and 24 (Y) are read
	if (chrom != 23 && chrom != 24) return;
	const char *col[9];
	const char *format_end = vcf_columns(aTHX_ r, p, f, e, col);
	UV pos = vcf_pos(aTHX_ r, col[1], col[2] - 1);
	if (sex_par(chrom, pos) || !vcf_base(col[3], col[4] - 1) || !vcf_base(col[4], col[5] - 1)) return;
	size_t gt, dp;
	vcf_format(col[8], format_end, &gt, &dp);
	if (gt == VCF_ABSENT) return;
	if (chrom == 23) t->x_sites++;
	p = format_end + 1;
	for (size_t i = 0; i < r->n; i++) {
		const char *se = vcf_sample(aTHX_ r, p, e, i), *sample = p;
		p = se + 1;
		unsigned short int alleles, alt;
		short int call = vcf_call(r, sample, se, gt, dp, &alleles, &alt);
		if (call == VCF_SHALLOW) continue;
		UV *c = t->count + i * SEX_FIELDS;
		if (chrom == 24) {
			if (call == VCF_NO_CALL) continue;
			c[SEX_Y_CALLED]++;
			if (alt) c[SEX_Y_ALT]++;
		} else if (call == VCF_NO_CALL) {
			c[SEX_X_MISSING]++;
		} else {
			if (alleles == 1) c[SEX_X_HAPLOID]++;
			// a haploid call is hom-ALT or hom-REF, as peddy's cyvcf2 counts it
			c[alt == 0 ? SEX_X_HOM_REF : alt == alleles ? SEX_X_HOM_ALT : SEX_X_HET]++;
		}
	}
	vcf_samples_end(aTHX_ r, p, e);
}

/*The counts of every sample of the VCF name, which next returns a block of
at a time. Returns sample => hash of counts, and sets *x_sites and *names,
the samples in file order, for the report.*/
static HV *sex_tally_vcf(pTHX_ SV *next, const char *name, UV min_depth, UV *x_sites, AV **names) {
	sex_tally t;
	Zero(&t, 1, sex_tally);
	t.r.name = name;
	t.r.min_depth = min_depth;
	t.r.on_header = sex_header;
	t.r.on_record = sex_record;
	vcf_read(aTHX_ next, &t.r);
	HV *out = t.r.samples;
	for (size_t i = 0; i < t.r.n; i++) {
		SV **sample = av_fetch(t.r.names, (SSize_t)i, 0);
		if (!sample) croak("%s: sample %" UVuf " of %s has lost its name", __func__, (UV)i, name); // cannot happen: vcf_header() pushed every one
		HV *row = newHV();
		for (unsigned short int k = 0; k < SEX_FIELDS; k++) (void)hv_store(row, SEX_NAME[k], (I32)strlen(SEX_NAME[k]), newSVuv(t.count[i * SEX_FIELDS + k]), 0);
		(void)hv_store_ent(out, *sample, newRV_noinc((SV*)row), 0); // in place of the column vcf_header() stored
	}
	*x_sites = t.x_sites;
	*names = t.r.names;
	return out;
}

/*The genotype counts behind parent(): how many autosomal biallelic SNVs the
child and the parent have with each pair of ALT counts, 0 to 2. perl takes
IBS0, IBS1, IBS2 and the KING-robust kinship from those nine numbers.

The two samples are either columns of one VCF, which is read a record at a
time and nothing kept but the nine counts, or each the one sample of its own
VCF. A VCF of one genome usually lists only the sites where its sample
carries an ALT allele, so that a site absent from it stands for hom-REF; the
two are then compared site by site, which takes each file's sites in memory,
8 bytes a site, sorted (by C89's qsort(), since the two files need not list
their chromosomes in the same order). A position is left out when either file
has more than one SNV there, a call that is not a diploid call of REF and
ALT, a DP below min_depth, or more than one ALT, or when the files give it
different REF or ALT alleles. A file that leaves out filtered sites, or
off-target ones, makes the comparison see hom-REF where the sample was never
called; perl's report gives the sites each file had alone, for that reason.*/

#define PARENT_NONE 3 // the code of a sample with no diploid call of REF and ALT, or a DP below min_depth

/*A site of an own VCF is one uint64_t, packed so that sorting the numbers
sorts the sites by chromosome and then position, and the sites at one
position share every bit from bit 7 up:

chromosome (1-22) << 39 | POS (32 bits) << 7 | REF << 5 | ALT << 2 | code

where REF and ALT are vcf_base_index()'s 0 to 3, ALT is 4 for a record of
more than one ALT, and code is the ALT count or PARENT_NONE.*/
#define PARENT_POS_SHIFT 7
#define PARENT_ALLELES(k) (((k) >> 2) & 0x1f) // REF and ALT together
#define PARENT_ALT(k) (((k) >> 2) & 7)
#define PARENT_CODE(k) ((unsigned short int)((k) & 3))

/*An own VCF's sites, in a holder that the scope's destructor frees, since
Renew() moves the array.*/
typedef struct {
	uint64_t *site;
	size_t n, room; // sites used and allocated
} parent_sites;

static void parent_sites_free(pTHX_ void *p) {
	parent_sites *s = (parent_sites*)p;
	Safefree(s->site);
	Safefree(s);
}

typedef struct {
	vcf_reader r;        // first, so that the vcf_reader * that the callbacks are given is this
	SV *want[2];         // one VCF: the child's and the parent's sample names; unused for own VCFs
	size_t col[2];       // one VCF: their columns among the samples
	UV t[9];             // one VCF: t[3 * child + parent], sites by the two codes
	parent_sites *sites; // own VCF: the file's sites
} parent_reader;

/*The code of the call in the sample field p..e: its ALT count, or
PARENT_NONE.*/
static unsigned short int parent_code(const vcf_reader *r, const char *p, const char *e, size_t gt, size_t dp) {
	unsigned short int alleles, alt;
	return vcf_call(r, p, e, gt, dp, &alleles, &alt) == VCF_CALL && alleles == 2 ? alt : PARENT_NONE;
}

/*One VCF: where the child's and the parent's columns are.*/
static void parent_header_one(pTHX_ vcf_reader *r) {
	parent_reader *t = (parent_reader*)r;
	for (unsigned short int k = 0; k < 2; k++) {
		HE *he = hv_fetch_ent(r->samples, t->want[k], 0, 0);
		if (he) {
			t->col[k] = (size_t)SvUV(HeVAL(he));
			continue;
		}
		SV *some = sv_2mortal(newSVpvs("")); // the first 5 samples, quoted, which is enough to show how they are named
		for (SSize_t i = 0; i <= av_len(r->names) && i < 5; i++) {
			SV **s = av_fetch(r->names, i, 0);
			sv_catpvf(some, "%s\"%s\"", i ? ", " : "", s ? SvPV_nolen(*s) : "");
		}
		if (r->n > 5) sv_catpvf(some, " and %" UVuf " more", (UV)(r->n - 5));
		const char *want = SvPV_nolen(t->want[k]), *list = SvPV_nolen(some);
		croak("%s: %s has no sample \"%s\"; its samples are %s", __func__, r->name, want, list);
	}
}

/*One VCF: the record p..e.*/
static void parent_record_one(pTHX_ vcf_reader *r, const char *p, const char *e) {
	parent_reader *t = (parent_reader*)r;
	const char *f = vcf_field_end(p, e);
	unsigned short int chrom = vcf_chrom(p, f);
	if (chrom < 1 || chrom > 22) return;
	const char *col[9];
	const char *format_end = vcf_columns(aTHX_ r, p, f, e, col);
	if (!vcf_base(col[3], col[4] - 1) || !vcf_base(col[4], col[5] - 1)) return;
	size_t gt, dp;
	vcf_format(col[8], format_end, &gt, &dp);
	if (gt == VCF_ABSENT) return;
	unsigned short int code[2] = { PARENT_NONE, PARENT_NONE }; // the child's and the parent's
	p = format_end + 1;
	for (size_t i = 0; i < r->n; i++) {
		const char *se = vcf_sample(aTHX_ r, p, e, i);
		for (unsigned short int k = 0; k < 2; k++) {
			if (i == t->col[k]) code[k] = parent_code(r, p, se, gt, dp);
		}
		p = se + 1;
	}
	vcf_samples_end(aTHX_ r, p, e);
	if (code[0] != PARENT_NONE && code[1] != PARENT_NONE) t->t[3 * code[0] + code[1]]++;
}

/*Own VCF: one sample, or it isn't.*/
static void parent_header_own(pTHX_ vcf_reader *r) {
	if (r->n != 1) croak("%s: %s has %" UVuf " samples, and a VCF given as \"child\" or \"parent\" must have one; give a VCF of several as \"vcf\", and the samples by name", __func__, r->name, (UV)r->n);
}

/*Own VCF: the record p..e, kept as a site.*/
static void parent_record_own(pTHX_ vcf_reader *r, const char *p, const char *e) {
	parent_reader *t = (parent_reader*)r;
	const char *f = vcf_field_end(p, e);
	unsigned short int chrom = vcf_chrom(p, f);
	if (chrom < 1 || chrom > 22) return;
	const char *col[9];
	const char *format_end = vcf_columns(aTHX_ r, p, f, e, col);
	if (!vcf_base(col[3], col[4] - 1)) return;
	unsigned short int alt; // vcf_base_index()'s, or 4 for more than one ALT
	if (vcf_base(col[4], col[5] - 1)) alt = vcf_base_index(*col[4]);
	else if (memchr(col[4], ',', (size_t)(col[5] - 1 - col[4]))) alt = 4;
	else return; // an indel, or an ALT of "." or "<NON_REF>"
	UV pos = vcf_pos(aTHX_ r, col[1], col[2] - 1);
	if (pos > UINT32_MAX) croak("%s: %s line %" UVuf " has POS %" UVuf ", beyond the 4294967295 that parent can compare", __func__, r->name, r->line, pos);
	size_t gt, dp;
	vcf_format(col[8], format_end, &gt, &dp);
	p = format_end + 1;
	const char *se = vcf_sample(aTHX_ r, p, e, 0);
	vcf_samples_end(aTHX_ r, se + 1, e);
	unsigned short int code = alt == 4 || gt == VCF_ABSENT ? PARENT_NONE : parent_code(r, p, se, gt, dp);
	parent_sites *s = t->sites;
	if (s->n == s->room) {
		s->room = s->room ? 2 * s->room : 65536; // 512 KiB, a sixth of an exome's 280,000 records
		Renew(s->site, s->room, uint64_t);
	}
	s->site[s->n++] = (uint64_t)chrom << 39 | (uint64_t)pos << PARENT_POS_SHIFT | (uint64_t)vcf_base_index(*col[3]) << 5 | (uint64_t)alt << 2 | code;
}

static int parent_cmp(const void *a, const void *b) {
	const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
	return (x > y) - (x < y);
}

/*Own VCF: the sites of the VCF name, which next returns a block of at a
time, sorted.*/
static parent_sites *parent_sites_vcf(pTHX_ SV *next, const char *name, UV min_depth) {
	parent_reader t;
	Zero(&t, 1, parent_reader);
	Newxz(t.sites, 1, parent_sites);
	SAVEDESTRUCTOR_X(parent_sites_free, t.sites); // at the caller's LEAVE, croak or not
	t.r.name = name;
	t.r.min_depth = min_depth;
	t.r.on_header = parent_header_own;
	t.r.on_record = parent_record_own;
	vcf_read(aTHX_ next, &t.r);
	if (t.sites->n) qsort(t.sites->site, t.sites->n, sizeof *t.sites->site, parent_cmp);
	return t.sites;
}

/*The sites of the child's VCF, a, and of the parent's, b, both sorted, as t
counts them; in[0], in[1] and in[2] get the sites where either carries the
ALT that were in both files, in the child's alone, and in the parent's alone.
The rules are in the comment above PARENT_NONE.*/
static void parent_merge(const parent_sites *a, const parent_sites *b, UV *t, UV *in) {
	size_t i = 0, j = 0;
	while (i < a->n || j < b->n) {
		const uint64_t pa = i < a->n ? a->site[i] >> PARENT_POS_SHIFT : UINT64_MAX, pb = j < b->n ? b->site[j] >> PARENT_POS_SHIFT : UINT64_MAX;
		const uint64_t at = pa < pb ? pa : pb; // the position to compare
		size_t na = 0, nb = 0;                 // its sites in each file
		while (i + na < a->n && a->site[i + na] >> PARENT_POS_SHIFT == at) na++;
		while (j + nb < b->n && b->site[j + nb] >> PARENT_POS_SHIFT == at) nb++;
		const uint64_t ka = na ? a->site[i] : 0, kb = nb ? b->site[j] : 0;
		i += na;
		j += nb;
		if (na > 1 || nb > 1) continue;
		if (na && (PARENT_CODE(ka) == PARENT_NONE || PARENT_ALT(ka) == 4)) continue;
		if (nb && (PARENT_CODE(kb) == PARENT_NONE || PARENT_ALT(kb) == 4)) continue;
		if (na && nb && PARENT_ALLELES(ka) != PARENT_ALLELES(kb)) continue;
		const unsigned short int c = na ? PARENT_CODE(ka) : 0, p = nb ? PARENT_CODE(kb) : 0; // a site absent from a file is hom-REF there
		t[3 * c + p]++;
		if (c || p) in[na && nb ? 0 : na ? 1 : 2]++;
	}
}

MODULE = Cn::Bioinf	PACKAGE = Cn::Bioinf
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

void
_sex_tally(next, name, min_depth)
	SV *next
	const char *name
	UV min_depth
PPCODE:
	UV x_sites;
	AV *names;
	HV *out = sex_tally_vcf(aTHX_ next, name, min_depth, &x_sites, &names);
	EXTEND(SP, 3);
	PUSHs(sv_2mortal(newRV_inc((SV*)out)));
	PUSHs(sv_2mortal(newSVuv(x_sites)));
	PUSHs(sv_2mortal(newRV_inc((SV*)names)));

void
_parent_one(next, name, min_depth, child, parent)
	SV *next
	const char *name
	UV min_depth
	SV *child
	SV *parent
PPCODE:
	parent_reader t;
	Zero(&t, 1, parent_reader);
	t.r.name = name;
	t.r.min_depth = min_depth;
	t.r.on_header = parent_header_one;
	t.r.on_record = parent_record_one;
	t.want[0] = child;
	t.want[1] = parent;
	vcf_read(aTHX_ next, &t.r);
	UV in = 0; // the sites where either carries the ALT, every one of them in both, since there is one file
	for (unsigned short int k = 1; k < 9; k++) in += t.t[k];
	EXTEND(SP, 12);
	for (unsigned short int k = 0; k < 9; k++) PUSHs(sv_2mortal(newSVuv(t.t[k])));
	PUSHs(sv_2mortal(newSVuv(in)));
	PUSHs(sv_2mortal(newSVuv(0)));
	PUSHs(sv_2mortal(newSVuv(0)));

void
_parent_own(child_next, child_name, parent_next, parent_name, min_depth)
	SV *child_next
	const char *child_name
	SV *parent_next
	const char *parent_name
	UV min_depth
PPCODE:
	const parent_sites *a = parent_sites_vcf(aTHX_ child_next, child_name, min_depth);
	const parent_sites *b = parent_sites_vcf(aTHX_ parent_next, parent_name, min_depth);
	UV t[9] = { 0 }, in[3] = { 0 }; // as parent_merge() fills them
	parent_merge(a, b, t, in);
	EXTEND(SP, 12);
	for (unsigned short int k = 0; k < 9; k++) PUSHs(sv_2mortal(newSVuv(t[k])));
	for (unsigned short int k = 0; k < 3; k++) PUSHs(sv_2mortal(newSVuv(in[k])));
