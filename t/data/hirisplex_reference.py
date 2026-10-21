# hirisplex_genotypes's dosages, worked out independently with cyvcf2, for
# t/hirisplex.t's expected values. It reads share/hirisplex_s_panel.tsv and a
# VCF and prints, for each sample, the count of each marker's counted allele
# or "NA", in the panel's order, as the webtool's CSV holds them. The tests
# never run it. Run it from the distribution root as
#
#   python -m venv venv && venv/bin/pip install cyvcf2
#   venv/bin/python -I t/data/hirisplex_reference.py t/data/hirisplex.vcf 7
#
# (the second argument is min_depth). On 2026-10-20, with cyvcf2 0.34.0 under
# Python 3.14.2, it printed what t/hirisplex.t expects.
#
# It is a second implementation of the same rules, not a copy of the first:
# cyvcf2 parses the VCF, and the rules it applies are those Bioinf.xs's
# comment above hps_site gives, namely
#
# - the counted allele is looked for among REF and the ALTs, so it may be REF;
# - a record whose REF is longer than one base is passed over, except at the
#   one marker that is an indel, where what is counted is the ALT longer than
#   REF;
# - where a site is split into one record for each ALT, the record carrying
#   the counted allele is the one that counts, whichever order they come in;
# - a genotype that is missing, partly missing, not diploid, or of a DP below
#   min_depth is no genotype at all.
import sys
from cyvcf2 import VCF

vcf_path, min_depth = sys.argv[1], int(sys.argv[2])

panel = []
with open('share/hirisplex_s_panel.tsv') as fh:
    cols = fh.readline().lstrip('#').rstrip('\n').split('\t')
    for line in fh:
        line = line.rstrip('\n')
        if line:
            panel.append(dict(zip(cols, line.split('\t'))))

at = {}  # (chrom, pos) -> list of marker indexes
for k, m in enumerate(panel):
    if m['b37_chrom']:
        at.setdefault((m['b37_chrom'], int(m['b37_pos'])), []).append(k)

vcf = VCF(vcf_path)
samples = list(vcf.samples)
dose = [["NA"] * len(panel) for _ in samples]
carried = [False] * len(panel)
seen = [False] * len(panel)

for v in vcf:
    chrom = v.CHROM[3:] if v.CHROM.lower().startswith('chr') else v.CHROM
    for k in at.get((chrom, v.POS), []):
        m = panel[k]
        ref, alts = v.REF, list(v.ALT)
        insertion = m['kind'] == 'insertion'
        if not insertion and len(ref) != 1:
            continue  # some other variant sharing the position
        idx = -1
        if insertion:
            for i, a in enumerate(alts):
                if len(a) > len(ref):
                    idx = i + 1
                    break
        elif ref == m['count_allele']:
            idx = 0
        else:
            for i, a in enumerate(alts):
                if a == m['count_allele']:
                    idx = i + 1
                    break
        if seen[k] and (carried[k] or idx < 0):
            continue
        seen[k] = True
        if idx >= 0:
            carried[k] = True
        depths = v.format('DP') if 'DP' in (v.FORMAT or []) else None
        for s, gt in enumerate(v.genotypes):
            alleles = gt[:-1]  # the last element says whether the call is phased
            if depths is not None and depths[s][0] >= 0 and depths[s][0] < min_depth:
                dose[s][k] = "NA"
            elif len(alleles) != 2 or any(a < 0 for a in alleles):
                dose[s][k] = "NA"
            else:
                dose[s][k] = str(sum(1 for a in alleles if a == idx)) if idx >= 0 else "0"

print(','.join(['sampleid'] + [m['rsid'] + '_' + m['input_allele'] for m in panel]))
for s, name in enumerate(samples):
    print(','.join([name] + dose[s]))
print('# markers with no record in the file: '
      + (', '.join(panel[k]['rsid'] for k in range(len(panel)) if not seen[k]) or 'none'))
