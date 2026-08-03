#!/usr/bin/env python3
"""Draw an SVG chart of the Adams-Novikov E2 page from mr_BP's output tables.

Run this in the directory holding a finished mr_BP run:

    ./anss_chart.py 25                  # -> 25_anss_E2.svg

The only inputs are the text tables mr_BP already writes; nothing is
recomputed and nothing needs to be rebuilt.

    <halfT>_BPAANSS_table.txt    classes, and the algebraic Novikov differentials
    <halfT>_BPAANSS_a0.txt       multiplication by 3
    <halfT>_BPAANSS_h0.txt       multiplication by h0 = alpha_1

Grading
-------
Default (--grading anss) plots the Adams-Novikov bidegree (t-s, s): the stem
on the x-axis and the homological degree s on the y-axis, one dot per
generator, no 3-towers.  This is the convention used in the published ANSS
charts, and it is NOT the pair mr_BP prints.  mr_BP's "|deg=(a,b)" has
b = s + i, where i is the algebraic Novikov filtration (algNov.cpp:118), so
alpha_2 = v1^1[1-0] is printed at height 2 but belongs at height 1.  The
s used here is read off the "[s-n]" bracket in the class name instead; the
printed stem a is used as-is.

--grading algnov plots mr_BP's own (a, b) pair instead, keeps the 3-multiples,
and draws them as vertical 3-towers plus the algebraic Novikov differentials.
That is the right picture for debugging a run; it is not an ANSS chart.

Which classes are on the E2 page
--------------------------------
The algebraic Novikov SS converges to Ext_{BP*BP} = the ANSS E2 page, so a
class belongs on the chart exactly when it survives that spectral sequence:

  * Lines of the form "cycle <- tag |dr" record a differential.  Both the
    cycle and the tag die, so both are dropped.  (The tag never gets a line
    of its own -- SS_table::output, SS.h:139-146, only emits untagged
    entries -- so dropping the "<-" lines removes both.)
  * A class that appears as a *target* in the multiplication-by-3 table is
    3 times another class, i.e. a step up a 3-tower rather than a generator,
    so it is dropped too.  Note this is a data question, not a syntactic one:
    v0^1[1-1] carries a v0 in its name but is not a 3-multiple of anything
    surviving (its would-be predecessor [1-1] supports a differential), and
    it is exactly alpha_3 in stem 11.  Filtering on the name would delete it.
"""

import argparse
import os
import re
import sys
from collections import defaultdict

# ---------------------------------------------------------------- parsing

NAME = r'(?:v\d+\^\d+)*\[\d+-\d+\]'
CLASS_RE = re.compile(
    rf'^(?P<cyc>{NAME})'
    rf'(?:\t<-\t(?P<tag>{NAME})\t\|d(?P<dr>\d+))?'
    r'\t\|deg=\((?P<stem>-?\d+),(?P<filt>-?\d+)\)$')
MULT_RE = re.compile(rf'^(?P<src>{NAME})\t->\t(?P<targets>.*?)\+?o$')
BRACKET = re.compile(r'\[(\d+)-(\d+)\]$')


def parse_table(path):
    """Return (classes, differentials).

    classes maps name -> {name, stem, algnov, s, gen, killed}
    differentials is a list of (source_name, target_name, r).
    """
    classes, diffs, bad = {}, [], 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            m = CLASS_RE.match(line)
            if not m:
                bad += 1
                sys.stderr.write(f'warning: unparsed line: {line!r}\n')
                continue
            cyc, stem, filt = m['cyc'], int(m['stem']), int(m['filt'])
            s, gen = map(int, BRACKET.search(cyc).groups())
            classes[cyc] = dict(name=cyc, stem=stem, algnov=filt, s=s, gen=gen,
                                killed=bool(m['tag']))
            if m['tag']:
                # The source of a differential is a tagged entry and so never
                # gets a line of its own; reconstruct it.  A d_r on this chart
                # runs (stem, f) -> (stem - 1, f + r).
                tag, dr = m['tag'], int(m['dr'])
                ts, tgen = map(int, BRACKET.search(tag).groups())
                classes.setdefault(tag, dict(name=tag, stem=stem + 1,
                                             algnov=filt - dr, s=ts, gen=tgen,
                                             killed=True))
                diffs.append((tag, cyc, dr))
    if bad:
        sys.stderr.write(f'warning: {bad} unparsed line(s) in {path}\n')
    return classes, diffs


def parse_mult(path):
    """Return the list of (source, target) pairs in a multiplication table.

    "-> o" means the product is zero: the trailing 'o' is a terminator
    (multiplication.cpp:159-172), not a class.
    """
    edges = []
    if not os.path.exists(path):
        sys.stderr.write(f'note: {path} not found, skipping those lines\n')
        return edges
    with open(path) as fh:
        for line in fh:
            m = MULT_RE.match(line.rstrip('\n'))
            if not m:
                continue
            for target in m['targets'].split('+'):
                if target:
                    edges.append((m['src'], target))
    return edges


# ------------------------------------------------------------- selection

def e2_page(classes, a0_edges):
    """Restrict to what survives to the ANSS E2 page (see module docstring)."""
    three_multiples = {t for _, t in a0_edges}
    return {n: c for n, c in classes.items()
            if not c['killed'] and n not in three_multiples}


# --------------------------------------------------------------- drawing

# Colours sampled from the reference chart.
C_GRID = '#e6e6e6'
C_AXIS = '#9a9a9a'
C_DOT = '#7dded9'
C_DOT_TOWER = '#979797'
C_STRUCT = '#e8d9c5'
C_DIFF = '#c294f0'
C_TEXT = '#666666'


def esc(s):
    return (s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))


def render(classes, diffs, struct, towers, opts):
    U = opts.unit
    pad_l, pad_b, pad_t, pad_r = 52, 42, 34, 24
    ykey = opts.ykey

    xs = [c['stem'] for c in classes.values()] or [0]
    ys = [c[ykey] for c in classes.values()] or [0]
    x0, x1 = min(min(xs), 0), max(xs)
    y0, y1 = min(min(ys), 0), max(ys)
    W = (x1 - x0 + 1) * U + pad_l + pad_r
    H = (y1 - y0 + 1) * U + pad_t + pad_b

    def X(stem):
        return pad_l + (stem - x0 + 0.5) * U

    def Y(f):
        return H - pad_b - (f - y0 + 0.5) * U

    # Several classes can share a bidegree: spread them within the cell.
    cells = defaultdict(list)
    for c in classes.values():
        cells[(c['stem'], c[ykey])].append(c)
    for group in cells.values():
        group.sort(key=lambda c: (c['s'], c['gen'], c['name']))
        n = len(group)
        step = min(0.30 * U, (0.62 * U) / n)
        for i, c in enumerate(group):
            c['px'] = X(c['stem']) + (i - (n - 1) / 2) * step
            c['py'] = Y(c[ykey])

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" '
         f'height="{H:.0f}" viewBox="0 0 {W:.0f} {H:.0f}" '
         f'font-family="Georgia,serif">',
         f'<rect width="100%" height="100%" fill="#ffffff"/>']

    tick = opts.tick or (2 if x1 - x0 <= 40 else 10)
    o.append(f'<g stroke="{C_GRID}" stroke-width="0.8">')
    for stem in range(x0, x1 + 1):
        if stem % 2 == 0:
            o.append(f'<line x1="{X(stem):.1f}" y1="{pad_t:.1f}" '
                     f'x2="{X(stem):.1f}" y2="{H - pad_b:.1f}"/>')
    for f in range(y0, y1 + 1):
        if f % 2 == 0:
            o.append(f'<line x1="{pad_l:.1f}" y1="{Y(f):.1f}" '
                     f'x2="{W - pad_r:.1f}" y2="{Y(f):.1f}"/>')
    o.append('</g>')

    o.append(f'<g font-size="10" fill="{C_TEXT}">')
    for stem in range(x0, x1 + 1):
        if stem % tick == 0:
            o.append(f'<text x="{X(stem):.1f}" y="{H - pad_b + 15:.1f}" '
                     f'text-anchor="middle">{stem}</text>')
    for f in range(y0, y1 + 1):
        if f % tick == 0 or tick > 2:
            if f % (2 if tick <= 2 else tick) == 0:
                o.append(f'<text x="{pad_l - 9:.1f}" y="{Y(f) + 3.5:.1f}" '
                         f'text-anchor="end">{f}</text>')
    o.append('</g>')

    o.append(f'<line x1="{pad_l:.1f}" y1="{H - pad_b:.1f}" x2="{W - pad_r:.1f}" '
             f'y2="{H - pad_b:.1f}" stroke="{C_AXIS}" stroke-width="1"/>')
    o.append(f'<line x1="{pad_l:.1f}" y1="{pad_t:.1f}" x2="{pad_l:.1f}" '
             f'y2="{H - pad_b:.1f}" stroke="{C_AXIS}" stroke-width="1"/>')

    def line(a, b, **kw):
        A, B = classes.get(a), classes.get(b)
        if not A or not B or 'px' not in A or 'px' not in B:
            return None
        attrs = ' '.join(f'{k.replace("_", "-")}="{v}"' for k, v in kw.items())
        return (f'<line x1="{A["px"]:.1f}" y1="{A["py"]:.1f}" '
                f'x2="{B["px"]:.1f}" y2="{B["py"]:.1f}" {attrs}/>')

    for a, b in towers:
        ln = line(a, b, stroke=C_DOT_TOWER, stroke_width=1.0)
        if ln:
            o.append(ln)
    for a, b in struct:
        ln = line(a, b, stroke=C_STRUCT, stroke_width=1.5)
        if ln:
            o.append(ln)
    for a, b, r in diffs:
        ln = line(a, b, stroke=C_DIFF, stroke_width=1.3, stroke_dasharray='5,3')
        if ln:
            o.append(ln)

    r = opts.dot_radius or max(1.8, 0.075 * U)
    for c in sorted(classes.values(), key=lambda c: (c['stem'], c[ykey])):
        fill = C_DOT_TOWER if c.get('is_tower') else C_DOT
        o.append(f'<circle cx="{c["px"]:.1f}" cy="{c["py"]:.1f}" r="{r:.2f}" '
                 f'fill="{fill}" stroke="#3f8f8b" stroke-width="0.5">'
                 f'<title>{esc(c["name"])}  (stem {c["stem"]}, s={c["s"]}, '
                 f'algNov {c["algnov"]})</title></circle>')

    if opts.title:
        o.append(f'<text x="{pad_l:.1f}" y="{pad_t - 12:.1f}" font-size="14" '
                 f'fill="#222">{esc(opts.title)}</text>')
    o.append(f'<text x="{W - pad_r:.1f}" y="{H - 8:.1f}" font-size="10" '
             f'text-anchor="end" fill="{C_TEXT}">stem  t-s</text>')
    o.append('</svg>')
    return '\n'.join(o)


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description='Chart the Adams-Novikov E2 page from mr_BP output.',
        epilog='Run in the directory containing a finished mr_BP run.')
    ap.add_argument('halfT', help="mr_BP's first argument, e.g. 25")
    ap.add_argument('-d', '--dir', default='.', help='where the tables live')
    ap.add_argument('-o', '--output', help='default <halfT>_anss_E2.svg')
    ap.add_argument('--grading', choices=['anss', 'algnov'], default='anss',
                    help='anss (default): (t-s, s), one dot per generator. '
                         'algnov: mr_BP\'s printed (t-s, s+i) with 3-towers.')
    ap.add_argument('--omit-stem0', action='store_true',
                    help='drop stem 0 (Ext^0 = Z_(3)), as the published '
                         'charts do')
    ap.add_argument('--max-stem', type=int, help='crop the chart')
    ap.add_argument('--max-filt', type=int, help='crop the chart')
    ap.add_argument('--unit', type=float, default=31.2,
                    help='px per lattice step (default 31.2, as in the '
                         'published charts)')
    ap.add_argument('--dot-radius', type=float)
    ap.add_argument('--tick', type=int, help='label every N stems')
    ap.add_argument('--title', default=None)
    ap.add_argument('--no-title', action='store_true')
    a = ap.parse_args()

    base = os.path.join(a.dir, f'{a.halfT}_BP')
    table = base + 'AANSS_table.txt'
    if not os.path.exists(table):
        sys.exit(f'error: {table} not found.\n'
                 f'Run ./mr_st {a.halfT} <s+1> && ./BPtab {a.halfT} && '
                 f'./mr_BP {a.halfT} <s> first, or pass --dir.')

    classes, diffs = parse_table(table)
    a0 = parse_mult(base + 'AANSS_a0.txt')
    h0 = parse_mult(base + 'AANSS_h0.txt')

    if a.grading == 'anss':
        a.ykey = 's'
        classes = e2_page(classes, a0)
        diffs, towers = [], []
        struct = [(x, y) for x, y in h0 if x in classes and y in classes]
    else:
        a.ykey = 'algnov'
        towers = [(x, y) for x, y in a0 if x in classes and y in classes]
        struct = [(x, y) for x, y in h0 if x in classes and y in classes]
        for _, t in towers:
            if t in classes:
                classes[t]['is_tower'] = True

    if a.omit_stem0:
        classes = {n: c for n, c in classes.items() if c['stem'] != 0}
    if a.max_stem is not None:
        classes = {n: c for n, c in classes.items() if c['stem'] <= a.max_stem}
    if a.max_filt is not None:
        classes = {n: c for n, c in classes.items() if c[a.ykey] <= a.max_filt}
    keep = set(classes)
    struct = [e for e in struct if e[0] in keep and e[1] in keep]
    towers = [e for e in towers if e[0] in keep and e[1] in keep]
    diffs = [d for d in diffs if d[0] in keep and d[1] in keep]

    if not classes:
        sys.exit('error: no classes to draw (is the run empty, or the crop too tight?)')

    if a.no_title:
        a.title = None
    elif a.title is None:
        grading = ('Adams-Novikov E2, p=3' if a.grading == 'anss'
                   else 'algebraic Novikov table, p=3')
        a.title = f'{grading}   (halfT={a.halfT}, {len(classes)} classes)'

    svg = render(classes, diffs, struct, towers, a)
    out = a.output or os.path.join(a.dir, f'{a.halfT}_anss_E2.svg')
    with open(out, 'w') as fh:
        fh.write(svg)

    stems = [c['stem'] for c in classes.values()]
    print(f'{len(classes)} classes (stems {min(stems)}-{max(stems)}), '
          f'{len(struct)} alpha_1 lines, {len(towers)} 3-tower lines, '
          f'{len(diffs)} differentials -> {out}')
    print('note: classes near the top stems are cut off by the degree bound, '
          'not absent.')


if __name__ == '__main__':
    main()
