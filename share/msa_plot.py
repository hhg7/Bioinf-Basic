"""Draw an aligned FASTA, a newick tree, or a coloured table of scores, as an
image in any format matplotlib writes (the output file's extension picks it:
.svg, .png, .pdf, ...).

Run by Bioinf::Basic::plot_msa, plot_phylo and msa_quality_table with the Python
from Alien::Bioinf. The
alignment drawing is ~/Scripts/HitList/scripts/x.my.align.py, cut down to what
it uses; that script was itself derived from CIAlign (Tumescheit, Firth & Brown,
PeerJ 2022, MIT licence), whose colour-blind-safe "CBS" palette and RasMol-based
residue colouring are copied below as flat tables. The tree drawing is new: it
replaces the R/ggtree script bioinf.pm used to write. The table follows the
script Matplotlib::Simple 0.317's colored_table wrote for msa_quality_table
before this took its place: a gist_rainbow-coloured matplotlib table over a
hidden imshow that carries the colorbar, with empty cells grey.
"""
import argparse, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, LogNorm, Normalize

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--f', help='aligned FASTA to draw')
ap.add_argument('--tree', help='newick tree to draw instead')
ap.add_argument('--o', required=True, help='output image')
ap.add_argument('--t', default='', help='title')
ap.add_argument('--x', default='Amino Acid Residue', help='x-axis label')
ap.add_argument('--y', default='Protein & Species', help='y-axis label')
ap.add_argument('--s', help='JSON {label: 0-based alignment column} of vertical lines to draw')
ap.add_argument('--l', help='JSON {tip name: label} for --tree')
ap.add_argument('--table', help='JSON file of the table to draw instead: rows, cols, cells (null for none), vmin, vmax, log, numbers, title, cblabel')
ap.add_argument('--c', default='Bioinf::Basic', help='Creator metadata: the script and sub that called this')
a = ap.parse_args()

# CIAlign getAAColours('CBS') and getNtColours('CBS'), in the order CIAlign
# assigns them to the colour map.
AA = dict(D='#a22c49', E='#a22c49', C='#c9c433', M='#c9c433', K='#0038a2', R='#0038a2', S='#e57700',
	T='#e57700', F='#589aab', Y='#589aab', N='#50d3cb', Q='#50d3cb', G='#eae2ea', L='#56ae6c', V='#56ae6c',
	I='#56ae6c', A='#888988', W='#89236a', H='#e669ca', P='#ffc4a9', X='#000000', B='#936e23', Z='#936e23',
	J='#936e23', U='#936e23', O='#936e23', **{'-': '#FFFFFF', '*': '#FFFFFF'})
NT = dict(A='#56ae6c', G='#c9c433', T='#a22c49', C='#0038a2', N='#6979d3', U='#a22c49',
	**{c: '#6979d3' for c in 'RYSWKMBDHVX'}, **{'-': '#FFFFFF'})

def save(fig):
	# 'Creator' is metadata only these backends take; jpg and the rest refuse it
	creator = a.c + ', drawn by ' + __file__ + ' with matplotlib ' + matplotlib.__version__
	md = {'Creator': creator} if a.o.rsplit('.', 1)[-1].lower() in ('png', 'svg', 'pdf', 'ps', 'eps') else None
	fig.savefig(a.o, bbox_inches='tight', pad_inches=0.1, metadata=md)
	print('wrote ' + a.o)

def read_fasta(path):
	names, seqs = [], []
	with open(path) as fh:
		for line in fh:
			line = line.strip()
			if line.startswith('>'):
				names.append(line[1:])
				seqs.append([])
			elif line:
				if not names:
					raise SystemExit(path + ' is not FASTA: sequence before the first ">" line')
				seqs[-1].append(line.upper())
	seqs = [''.join(s) for s in seqs]
	if len({len(s) for s in seqs}) > 1:
		raise SystemExit(path + ': the sequences are not all the same length, so they are not aligned')
	return names, seqs

def draw_alignment():
	names, seqs = read_fasta(a.f)
	# CIAlign's seqType(): nucleotide when every character is a nucleotide code
	chars = set(''.join(seqs))
	pal = NT if chars <= set(NT) else AA
	if not chars <= set(pal):
		raise SystemExit('not IUPAC amino acid or nucleotide codes: ' + ' '.join(sorted(chars - set(pal))))
	# rows flipped, as CIAlign draws them: the first sequence is the bottom row
	arr = np.flip(np.array([list(s) for s in seqs]), axis=0)
	present = [c for c in pal if c in chars]
	code = {c: i for i, c in enumerate(present)}
	num = np.vectorize(code.get)(arr)
	height, width = arr.shape
	f = plt.figure(figsize=(12, 3))
	ax = f.add_subplot(1, 1, 1)
	ax.set_xlim(-0.5, width)
	for side in ('right', 'top', 'left'):
		ax.spines[side].set_visible(False)
	f.suptitle(a.t, fontsize=11)
	ax.imshow(num, cmap=ListedColormap([pal[c] for c in present]), aspect='auto', interpolation='nearest')
	ax.set_yticks(list(range(height))[::-1])
	ax.set_yticklabels(names, fontsize=10)
	ax.set_xlabel(a.x)
	ax.set_ylabel(a.y)
	if a.s:
		try:
			from adjustText import adjust_text
		except ImportError:
			raise SystemExit('--s needs the adjustText module, which Alien::Bioinf installs')
		sites = json.loads(a.s)
		texts = [ax.text(col, -1, label) for label, col in sites.items()]
		ax.add_collection(LineCollection([[(col, height), (col, 0)] for col in sites.values()],
			colors='black', linewidth=3, alpha=0.5, linestyle='dashed'))
		adjust_text(texts, only_move={'text': 'x', 'static': 'x', 'explode': 'x', 'pull': 'x'})
	save(f)

def parse_newick(s):
	"""[name, branch length, children]; clustalo writes no quoted names."""
	s = ''.join(s.split()).rstrip(';')
	i = 0
	def node():
		nonlocal i
		kids = []
		if s[i] == '(':
			i += 1
			while True:
				kids.append(node())
				i += 1
				if s[i - 1] == ')':
					break
		j = i
		while i < len(s) and s[i] not in ',():':
			i += 1
		name, length = s[j:i], 0.0
		if i < len(s) and s[i] == ':':
			i += 1
			j = i
			while i < len(s) and s[i] not in ',()':
				i += 1
			length = float(s[j:i])
		return [name, length, kids]
	return node()

def draw_tree():
	labels = json.loads(a.l) if a.l else {}
	with open(a.tree) as fh:
		root = parse_newick(fh.read())
	tips, segs = [], []
	def place(n, x):
		# clustalo can write a slightly negative length; drawn as 0, as the
		# ggtree script's ignore.negative.edge=TRUE did
		x += max(n[1], 0.0)
		if not n[2]:
			tips.append((n[0], x))
			return x, len(tips) - 1
		ys = [place(k, x) for k in n[2]]
		for kx, ky in ys:
			segs.append([(x, ky), (kx, ky)])
		segs.append([(x, ys[0][1]), (x, ys[-1][1])])
		return x, (ys[0][1] + ys[-1][1]) / 2
	place(root, 0.0)
	right = max(x for _, x in tips)
	f, ax = plt.subplots(figsize=(10, 0.3 * len(tips) + 1))
	ax.add_collection(LineCollection(segs, colors='black', linewidth=1))
	# tip labels aligned at the right, joined to their tips by dotted lines
	ax.add_collection(LineCollection([[(x, y), (right, y)] for y, (_, x) in enumerate(tips)],
		colors='grey', linewidth=0.5, linestyle='dotted'))
	for y, (name, _) in enumerate(tips):
		ax.text(right * 1.02, y, labels.get(name, name), va='center')
	ax.set_xlim(0, right * 1.02 if right else 1)
	ax.set_ylim(len(tips) - 0.5, -0.5)
	ax.set_yticks([])
	for side in ('right', 'top', 'left'):
		ax.spines[side].set_visible(False)
	ax.set_xlabel('substitutions per site')
	ax.set_title(a.t)
	save(f)

def draw_table():
	with open(a.table) as fh:
		t = json.load(fh)
	d = np.array([[np.nan if v is None else v for v in row] for row in t['cells']], dtype=float)
	cmap = matplotlib.colormaps['gist_rainbow'].copy()
	cmap.set_bad('gray')
	norm = (LogNorm if t['log'] else Normalize)(vmin=t['vmin'], vmax=t['vmax'])
	f, ax = plt.subplots(1, 1, layout='constrained')
	# the table can't carry a colorbar, so an invisible image of the same
	# cells does
	img = ax.imshow(d, cmap=cmap, norm=norm)
	f.colorbar(img, label=t['cblabel'] or '')
	img.set_visible(False)
	text = [['' if v is None else str(v) for v in row] for row in t['cells']] if t['numbers'] else None
	ax.table(cellText=text, rowLabels=t['rows'], colLabels=t['cols'], cellColours=cmap(norm(d)), loc='center', bbox=[0, 0, 1, 1])
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_title(t['title'])
	save(f)

if a.table:
	draw_table()
elif a.tree:
	draw_tree()
elif a.f:
	draw_alignment()
else:
	raise SystemExit('give --f (an alignment), --tree (a newick file) or --table (a JSON table)')
