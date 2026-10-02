"""Draw an aligned FASTA, a newick tree, or a coloured table of scores, as an
image in any format matplotlib writes (the output file's extension picks it:
.svg, .png, .pdf, ...).

Run by Bioinf::Basic::plot_msa, plot_phylo and msa_quality_table with the Python
from Alien::Bioinf. The
alignment drawing is ~/Scripts/HitList/scripts/x.my.align.py, cut down to what
it uses; that script was itself derived from CIAlign (Tumescheit, Firth & Brown,
PeerJ 2022, MIT licence), whose colour-blind-safe "CBS" palette and RasMol-based
residue colouring are copied below as flat tables. The tree is read and drawn by
Biopython's Bio.Phylo (Cock et al., Bioinformatics 2009), onto matplotlib axes;
it replaces the R/ggtree script bioinf.pm used to write. The table follows the
script Matplotlib::Simple 0.317's colored_table wrote for msa_quality_table
before this took its place: a gist_rainbow-coloured matplotlib table over a
hidden imshow that carries the colorbar, with empty cells grey.
"""
import argparse, datetime, hashlib, json, os, platform, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, LogNorm, Normalize

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--f', help='aligned FASTA to draw')
ap.add_argument('--tree', help='newick tree to draw instead')
ap.add_argument('--o', help='output image (required)')
ap.add_argument('--t', default='', help='title')
ap.add_argument('--x', default='Amino Acid Residue', help='x-axis label')
ap.add_argument('--y', default='Protein & Species', help='y-axis label')
ap.add_argument('--s', help='JSON {label: 0-based alignment column} of vertical lines to draw')
ap.add_argument('--l', help='JSON {tip name: label} for --tree')
ap.add_argument('--table', help='JSON file of the table to draw instead: rows, cols, cells (null for none), vmin, vmax, log, numbers, title, cblabel')
ap.add_argument('--meta', help='JSON of further Dublin Core fields for --tree: Title, Source, Description, Relation, Contributor')
ap.add_argument('--c', default='Bioinf::Basic', help='Creator metadata: the script and sub that called this')
ap.add_argument('--quiet', action='store_true', help="don't print 'wrote' and the output file; the caller prints its own")
ap.add_argument('--argfile', help='a JSON array of every other argument, in place of them on the command line')
a = ap.parse_args()
# Bioinf::Basic passes everything through --argfile, because perl's system(LIST)
# on Windows does not escape a '"' inside an argument, and the JSON options and
# the Creator ("... called using \"plot_msa\" ...") are full of them.
if a.argfile:
	with open(a.argfile, encoding='utf-8') as fh:
		a = ap.parse_args(json.load(fh))
# checked here rather than with required=True, which would refuse --argfile alone
if not a.o:
	ap.error('--o is required')

# CIAlign getAAColours('CBS') and getNtColours('CBS'), in the order CIAlign
# assigns them to the colour map.
AA = dict(D='#a22c49', E='#a22c49', C='#c9c433', M='#c9c433', K='#0038a2', R='#0038a2', S='#e57700',
	T='#e57700', F='#589aab', Y='#589aab', N='#50d3cb', Q='#50d3cb', G='#eae2ea', L='#56ae6c', V='#56ae6c',
	I='#56ae6c', A='#888988', W='#89236a', H='#e669ca', P='#ffc4a9', X='#000000', B='#936e23', Z='#936e23',
	J='#936e23', U='#936e23', O='#936e23', **{'-': '#FFFFFF', '*': '#FFFFFF'})
NT = dict(A='#56ae6c', G='#c9c433', T='#a22c49', C='#0038a2', N='#6979d3', U='#a22c49',
	**{c: '#6979d3' for c in 'RYSWKMBDHVX'}, **{'-': '#FFFFFF'})

# The Dublin Core fields each backend takes, after matplotlib's savefig(): SVG
# writes them all into <metadata>, PNG writes any key as a tEXt chunk, PDF takes
# only its own Info keys (Description going in as Subject), and PS and EPS
# only Creator. jpg and the rest take none, and refuse any.
DC = ('Title', 'Date', 'Source', 'Description', 'Identifier', 'Relation', 'Contributor', 'Keywords')

def save(fig, also=None, dc=None):
	"""also: a further library that drew the image, as "name version"; dc:
	Dublin Core fields to write beside the Creator, where the format allows."""
	creator = a.c + ', drawn by ' + __file__ + ' with matplotlib ' + matplotlib.__version__
	if also:
		creator += ' and ' + also
	fmt = a.o.rsplit('.', 1)[-1].lower()
	dc = dc or {}
	if fmt == 'svg':
		md = {'Creator': creator, **dc}
	elif fmt == 'png':
		# tEXt values are plain strings
		md = {'Creator': creator, **{k: '; '.join(v) if isinstance(v, list) else str(v) for k, v in dc.items()}}
	elif fmt == 'pdf':
		md = {'Creator': creator, **{k: v for k, v in (('Title', dc.get('Title')), ('Subject', dc.get('Description')),
			('Keywords', '; '.join(dc.get('Keywords', [])) or None)) if v}}
	elif fmt in ('ps', 'eps'):
		md = {'Creator': creator}
	else:
		md = None
	fig.savefig(a.o, bbox_inches='tight', pad_inches=0.1, metadata=md)
	if not a.quiet:
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

def draw_tree():
	try:
		import Bio
		from Bio import Phylo
	except ImportError:
		raise SystemExit('--tree needs Biopython, which Alien::Bioinf installs')
	labels = json.loads(a.l) if a.l else {}
	with open(a.tree, 'rb') as fh:
		raw = fh.read()
	newick = raw.decode()
	tree = Phylo.read(a.tree, 'newick')
	# clustalo can write a slightly negative length; drawn as 0, as the
	# ggtree script's ignore.negative.edge=TRUE did
	clamped = 0
	for clade in tree.find_clades():
		if clade.branch_length is not None and clade.branch_length < 0:
			clade.branch_length = 0.0
			clamped += 1
	tips = [c.name for c in tree.get_terminals()]
	f, ax = plt.subplots(figsize=(10, 0.3 * len(tips) + 1))
	# Phylo.draw labels every clade for which label_func is not None, and only
	# the tips have names
	Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
		label_func=lambda c: labels.get(c.name, c.name) if c.is_terminal() else None)
	ax.set_yticks([])
	ax.set_ylabel('')
	for side in ('right', 'top', 'left'):
		ax.spines[side].set_visible(False)
	ax.set_xlabel('substitutions per site')
	ax.set_title(a.t)
	dc = json.loads(a.meta) if a.meta else {}
	# to the second, with the UTC offset: matplotlib's own SVG Date is the day
	# alone. SOURCE_DATE_EPOCH, where set, is left to matplotlib instead, so
	# that a reproducible build stays reproducible.
	if 'SOURCE_DATE_EPOCH' not in os.environ:
		dc['Date'] = datetime.datetime.now().astimezone().isoformat(timespec='seconds')
	dc['Identifier'] = 'sha256:' + hashlib.sha256(raw).hexdigest() + ' of the newick file drawn'
	desc = dc.get('Description', '')
	if clamped:
		desc += ' ' + str(clamped) + ' negative branch length' + ('s were' if clamped > 1 else ' was') + ' drawn as 0.'
	shown = [t + (' (shown as ' + labels[t] + ')' if t in labels else '') for t in tips]
	dc['Description'] = (desc + ' ' + str(len(tips)) + ' tips: ' + ', '.join(shown)
		+ '. The tree as drawn, in newick: ' + ''.join(newick.split())).strip()
	dc['Keywords'] = ['phylogenetic tree', 'newick'] + tips
	dc['Contributor'] = dc.get('Contributor', []) + ['Python ' + platform.python_version() + ' (' + sys.executable + ')',
		'matplotlib ' + matplotlib.__version__, 'Biopython ' + Bio.__version__, 'NumPy ' + np.__version__]
	save(f, 'Biopython ' + Bio.__version__, {k: dc[k] for k in DC if k in dc})

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
