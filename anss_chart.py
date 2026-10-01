#!/usr/bin/env python3
"""Draw an SVG chart of the Adams-Novikov E2 page from mr_BP's output tables.

Run this in the directory holding a finished mr_BP run:

    ./anss_chart.py 25                  # -> 25_anss_E2.svg

The only inputs are the text tables mr_BP already writes; nothing is
recomputed and nothing needs to be rebuilt.

    <halfT>_BPAANSS_table.txt    classes, and the algebraic Novikov differentials
    <halfT>_BPAANSS_a0.txt       multiplication by 3  (characteristic 0 only;
                                 a run over BP_*/I_n writes AANSS_a<n>.txt,
                                 multiplication by v_n, instead)
    <halfT>_BPAANSS_h0.txt       multiplication by h0 = alpha_1
    <halfT>_BPrun_info.txt       which ring the run was done over

Any comodule resolved by mr_BP_comod works the same way via -c/--comodule;
nothing below is specific to the sphere.

Grading
-------
Default (--grading anss) plots the Adams-Novikov bidegree (t-s, s): the stem
on the x-axis and the homological degree s on the y-axis.  This is the
convention used in the published ANSS charts, and it is NOT the pair mr_BP
prints.  mr_BP's "|deg=(a,b)" has b = s + i, where i is the algebraic Novikov
filtration (algNov.cpp:118), so alpha_2 = v1^1[1-0] is printed at height 2
but belongs at height 1.  The s used here is read off the "[s-n]" bracket in
the class name instead; the printed stem a is used as-is.

--grading algnov plots mr_BP's own (a, b) pair instead, one dot per class,
draws the algebraic Novikov differentials, and joins each class to its
multiples by the bottom generator of the base ring's maximal invariant ideal
-- p = v_0 over BP_*, which makes vertical towers, or v_n over BP_*/I_n,
which does not (v_n raises the stem by |v_n| as well as the filtration).
That is the right picture for debugging a run; it is not an ANSS chart.

What one mark means (--grading anss)
------------------------------------
One mark per cyclic summand, not per F_3 basis element: the classes in the
table are a basis of gr Ext, so a Z/9 shows up as two of them and would be
two dots on a naive chart.  Multiplication by 3 is exactly what the a0 table
records, so the summands are recovered by chaining a0 edges.  The glyphs
follow Belmont's published 3-primary ANSS chart:

    filled square, marked oo    Z_(3)
    filled dot                  Z/3
    dot inside n-1 rings        Z/3^n
    dashed outermost ring       the tower runs off the end of the computed
                                range, so the order shown is a lower bound

Several summands in one bidegree are drawn side by side, each labelled with
the name of its generator (the class the rest of the tower is 3 times).
Labels are the table's own names, typeset as v_0v_1^3[1-1] rather than
v0^1v1^3[1-1]; --raw-labels keeps the raw text and --no-labels drops them.

Filtration 0 is decided by mathematics rather than by the table, but WHICH
mathematics depends on the ring the run was done over, which it reads from
<prefix>run_info.txt (written by BPInit):

    characteristic 0   the base ring is BP_*, M is free over it, so
                       Ext^0 = Prim(M) is a submodule of a free BP_*-module,
                       hence torsion-free: every summand there is Z_(3), a box.
    characteristic p   the base ring is BP_*/I_n (n >= 1), which is an
                       F_p-algebra, so p kills M and every Ext group -- Ext^0
                       included -- is an F_p-vector space: every summand
                       anywhere is Z/p, a dot.

Getting that wrong is silent, which is why it is read and not assumed; pass
--height for a run whose run_info.txt is missing.  Above filtration 0 in
characteristic 0 the order comes from the a0 chain, and where that chain hits
the truncation (b = s + i is cut off at the resolution length, see
docs/CHARTS.md) the mark gets the dashed ring.

Which classes are on the E2 page
--------------------------------
The algebraic Novikov SS converges to Ext_{BP*BP} = the ANSS E2 page, so a
class belongs on the chart exactly when it survives that spectral sequence:
lines of the form "cycle <- tag |dr" record a differential, and both the
cycle and the tag die, so both are dropped.  (The tag never gets a line of
its own -- SS_table::output, SS.h:139-146, only emits untagged entries -- so
dropping the "<-" lines removes both.)

Unlike the previous version of this script, 3-multiples are NOT dropped:
they are what determines each summand's isomorphism type, and they are now
drawn as part of one mark.
"""

import argparse
import os
import re
import sys
from collections import defaultdict

# ---------------------------------------------------------------- parsing

# this is the p=3 fork; the prime shows up in glyph names ("Z/3") and in the
# characteristic of the base ring of a height > 0 run
PRIME = 3

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


def parse_run_info(path):
    """Read the run-info file BPInit writes beside its tables.

    One "key value" per line, '#' starts a comment.  Returns None if the file
    is absent, which is the case for output produced before runs started
    recording this.
    """
    if not os.path.exists(path):
        return None
    info = {}
    with open(path) as fh:
        for line in fh:
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            key, _, value = line.partition(' ')
            info[key.strip()] = value.strip()
    return info


def looks_characteristic_p(classes):
    """True if the table itself is evidence that p is 0 in the base ring.

    `v0` in a class name is the p-adic valuation of that class's coefficient
    (algNov.cpp:28), not a polynomial generator. So a run over BP_* is full of
    them -- Ext^{0,0} = Z_(3) alone contributes a whole v_0-tower, and in
    practice half the classes carry one -- while a run over BP_*/I_n, where p
    is 0, cannot produce a single one.

    This is evidence, not proof: a height-0 run truncated before any v_0
    appeared would look the same. So the caller uses it to refuse to guess,
    not to decide.
    """
    return bool(classes) and not any('v0^' in name for name in classes)


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


def parse_mult_detail(path):
    """Return source -> [targets] for every source the table states.

    A source missing from the result had no line at all, which is *unknown*
    rather than zero: the tables are pruned by homological degree
    (multiplication.cpp:169).  A source mapping to [] is a stated zero.
    """
    if not os.path.exists(path):
        return None
    products = {}
    with open(path) as fh:
        for line in fh:
            m = MULT_RE.match(line.rstrip('\n'))
            if not m:
                continue
            products[m['src']] = [t for t in m['targets'].split('+') if t]
    return products


# ------------------------------------------------------- cyclic summands

def build_marks(classes, a0_products, ext0_torsion_free=True):
    """Group the surviving classes into cyclic summands, one mark each.

    Follows multiplication-by-3 chains: a class that is not 3 times any
    surviving class generates a summand, and the chain of its 3-multiples
    gives the order.  Returns a list of marks and a name -> mark index map.

    ext0_torsion_free says whether the base ring has characteristic 0; see
    make_mark, which is where it matters.
    """
    surv = {n: c for n, c in classes.items() if not c['killed']}
    maxb = max((c['algnov'] for c in surv.values()), default=0)
    if a0_products is None:
        # No multiplication-by-3 table: the isomorphism types are simply not
        # knowable, so draw every class as a plain dot rather than inventing
        # towers (the caller has already warned about this).
        marks = [make_mark(c, [n], True, ext0_torsion_free)
                 for n, c in surv.items()]
        return marks, {m['name']: i for i, m in enumerate(marks)}
    prod = a0_products

    # only edges out of a surviving class count; a 3-multiple of something
    # that died is a generator in its own right (CHARTS.md section 2)
    is_multiple = {t for src, ts in prod.items() if src in surv
                   for t in ts if t in surv}

    marks, owner, seen = [], {}, set()
    order = sorted(surv, key=lambda n: (surv[n]['stem'], surv[n]['s'],
                                        surv[n]['algnov'], n))
    for name in order:
        if name in seen or name in is_multiple:
            continue
        chain, cur, bounded = [name], name, None
        seen.add(name)
        while True:
            ts = prod.get(cur)
            if ts is None:                       # no line: top not visible
                bounded = False
                break
            if not ts:                           # stated "-> o": tower ends
                bounded = True
                break
            nxt = ts[0]
            if len(ts) != 1 or nxt not in surv or nxt in seen:
                bounded = False                  # not a simple chain
                break
            chain.append(nxt)
            seen.add(nxt)
            cur = nxt
        if bounded and surv[chain[-1]]['algnov'] >= maxb:
            bounded = False                      # ends at the truncation edge
        marks.append(make_mark(surv[name], chain, bounded, ext0_torsion_free))
        for member in chain:
            owner[member] = len(marks) - 1

    for name in order:                           # defensive: nothing orphaned
        if name not in owner:
            marks.append(make_mark(surv[name], [name], True, ext0_torsion_free))
            owner[name] = len(marks) - 1
    return marks, owner


def make_mark(gen, chain, bounded, ext0_torsion_free=True):
    """Classify one summand.

    s = 0 is decided by mathematics, not by the tables -- but only when that
    mathematics applies.  Over BP_* the comodule is free, so Ext^0 = Prim(M)
    is a submodule of a free BP_*-module, hence torsion-free, and every
    summand there is Z_(3) -- a box, whatever the a0 table can or cannot see.
    That is what ext0_torsion_free records, and it is FALSE for a run over
    BP_*/I_n: there the base ring is an F_p-algebra, p kills M, and Ext^0 is
    an F_p-vector space like every other Ext group, so its summands are Z/p
    dots.  Asserting the box unconditionally is how this used to draw the
    Ext^0 = F_3[v_1] of S/p as a row of Z_(3)'s.

    Above filtration 0 the order is read off the a0 chain.  If the chain runs
    into the truncation (its top has no a0 line, or sits in the last
    filtration the run prints) the visible length is only a lower bound, so
    the mark is drawn with a dashed outer ring and reported as ">=".
    """
    n = len(chain)
    if gen['s'] == 0 and ext0_torsion_free:
        kind, group, trunc = 'box', 'Z_(3)', False
    elif bounded:
        kind, group, trunc = ('dot' if n == 1 else 'circles'), f'Z/3^{n}', False
        if n == 1:
            group = 'Z/3'
    else:
        kind, group, trunc = 'circles', f'>= Z/3^{n} (truncated)', True
    return dict(name=gen['name'], stem=gen['stem'], s=gen['s'],
                algnov=gen['algnov'], kind=kind, n=n, group=group,
                truncated=trunc, members=chain)


# --------------------------------------------------------------- drawing

# Colours sampled from the reference chart (Belmont's 3-primary ANSS chart).
C_GRID = '#e6e6e6'
C_AXIS = '#9a9a9a'
C_DOT = '#808080'
C_EDGE = '#6e6e6e'
C_DOT_TOWER = '#979797'
C_STRUCT = '#e8d9c5'
C_DIFF = '#c294f0'
C_TEXT = '#666666'
C_LABEL = '#444444'

SUB = str.maketrans('0123456789', '₀₁₂₃₄'
                                  '₅₆₇₈₉')
SUP = str.maketrans('0123456789', '⁰¹²³⁴'
                                  '⁵⁶⁷⁸⁹')


def esc(s):
    return (s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))


def pretty(name):
    """v0^1v1^3[1-1] -> v(0)v(1)^3[1-1] with real sub/superscripts.

    Same names the tables use, typeset the way the published charts do: the
    generator index drops to a subscript and an exponent of 1 is left off.
    """
    def repl(m):
        idx, exp = m.group(1), m.group(2)
        out = 'v' + idx.translate(SUB)
        return out if exp == '1' else out + exp.translate(SUP)
    return re.sub(r'v(\d+)\^(\d+)', repl, name)


def mark_glyph(m, r, opts):
    """SVG for one summand.

    Filled square (with an infinity sign) = Z_(3); filled dot = Z/3; a dot
    inside n-1 concentric rings = Z/3^n.  A dashed outer ring means the tower
    runs past the end of the computed range, so the order shown is a lower
    bound.
    """
    x, y = m['px'], m['py']
    note = ''
    if m.get('cohabitants', 1) > 1:
        note = (f' -- splitting NOT determined: {m["cohabitants"]} summands '
                f'drawn in this bidegree, and a hidden extension could merge '
                f'them into a bigger cyclic group')
    title = (f'<title>{esc(pretty(m["name"]))}  {esc(m["group"])}  '
             f'(stem {m["stem"]}, s={m["s"]}){esc(note)}</title>')
    if m['kind'] == 'box':
        side = 2 * r
        body = (f'<rect x="{x - r:.1f}" y="{y - r:.1f}" width="{side:.1f}" '
                f'height="{side:.1f}" fill="{C_DOT}" stroke="{C_EDGE}" '
                f'stroke-width="0.6"/>'
                f'<text x="{x:.1f}" y="{y - r - 1.0:.1f}" font-size="{2.2 * r:.1f}" '
                f'text-anchor="middle" fill="{C_EDGE}">∞</text>')
    elif m['kind'] == 'dot':
        body = (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.2f}" fill="{C_DOT}" '
                f'stroke="none"/>')
    else:
        solid, dashed = rings_of(m, opts)
        body = (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.2f}" '
                f'fill="{C_DOT}" stroke="none"/>')
        for k in range(1, solid + dashed + 1):
            dash = (' stroke-dasharray="1.6,1.3"'
                    if dashed and k == solid + dashed else '')
            body += (f'<circle cx="{x:.1f}" cy="{y:.1f}" '
                     f'r="{r + opts.ring_gap * k:.2f}" fill="none" '
                     f'stroke="{C_EDGE}" stroke-width="0.9"{dash}/>')
        if m['n'] - 1 > solid:           # too many rings to draw legibly
            R = m['radius']
            body += (f'<text x="{x + R + 1:.1f}" y="{y - R:.1f}" font-size="5" '
                     f'text-anchor="start" fill="{C_EDGE}">{m["n"]}</text>')
    return f'<g>{title}{body}</g>'


def rings_of(m, opts):
    """(solid rings, dashed rings) around the central dot of a Z/3^n mark."""
    solid = min(m['n'] - 1, opts.max_rings - 1)
    return solid, (1 if m['truncated'] else 0)


def mark_radius(m, r, opts):
    """Outer radius of a mark, for label placement and line endpoints."""
    if m['kind'] != 'circles':
        return r
    solid, dashed = rings_of(m, opts)
    return r + opts.ring_gap * (solid + dashed)


def render(marks, diffs, struct, towers, opts):
    U = opts.unit
    pad_l, pad_b, pad_t, pad_r = 52, 46, 34, 30
    ykey = opts.ykey

    xs = [m['stem'] for m in marks] or [0]
    ys = [m[ykey] for m in marks] or [0]
    x0, x1 = min(min(xs), 0), max(xs)
    y0, y1 = min(min(ys), 0), max(ys)
    W = (x1 - x0 + 1) * U + pad_l + pad_r
    H = (y1 - y0 + 1) * U + pad_t + pad_b

    def X(stem):
        return pad_l + (stem - x0 + 0.5) * U

    def Y(f):
        return H - pad_b - (f - y0 + 0.5) * U

    base_r = opts.dot_radius or max(2.2, 0.075 * U)
    for m in marks:
        m['radius'] = mark_radius(m, base_r, opts)

    # Several summands can share a bidegree: spread them within the cell.
    cells = defaultdict(list)
    for m in marks:
        cells[(m['stem'], m[ykey])].append(m)
    for group in cells.values():
        group.sort(key=lambda m: (m['s'], m['algnov'], m['name']))
        n = len(group)
        step = min(0.34 * U, (0.72 * U) / n)
        for i, m in enumerate(group):
            m['px'] = X(m['stem']) + (i - (n - 1) / 2) * step
            m['py'] = Y(m[ykey])
            m['slot'] = i

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

    by_name = {m['name']: m for m in marks}

    def line(a, b, **kw):
        A, B = by_name.get(a), by_name.get(b)
        if not A or not B or 'px' not in A or 'px' not in B or A is B:
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

    if opts.bound and ykey == 's':
        # The run computes internal degree t = stem + s up to the degree
        # bound, so what it cuts off is a DIAGONAL, not the top of the
        # y-axis: (31,1) is in range at t=32 while (30,6) is not, at t=36.
        T = opts.bound
        xa, xb = max(x0, T - y1), min(x1, T - y0)
        if xa <= xb:
            o.append(f'<line x1="{X(xa):.1f}" y1="{Y(T - xa):.1f}" '
                     f'x2="{X(xb):.1f}" y2="{Y(T - xb):.1f}" '
                     f'stroke="{C_AXIS}" stroke-width="0.8" '
                     f'stroke-dasharray="4,4" opacity="0.55"/>')
            o.append(f'<text x="{X(xb) - 4:.1f}" y="{Y(T - xb) - 5:.1f}" '
                     f'font-size="8" text-anchor="end" fill="{C_TEXT}" '
                     f'opacity="0.8">t = {T} (degree bound)</text>')

    for m in sorted(marks, key=lambda m: (m['stem'], m[ykey])):
        o.append(mark_glyph(m, base_r, opts))

    if opts.labels:
        o.append(f'<g font-size="{opts.label_size}" fill="{C_LABEL}" '
                 f'font-family="monospace">')
        for m in marks:
            dy = m['radius'] + opts.label_size + 1 + (m['slot'] % 2) * \
                (opts.label_size + 0.5)
            text = m['name'] if opts.raw_labels else pretty(m['name'])
            o.append(f'<text x="{m["px"]:.1f}" y="{m["py"] + dy:.1f}" '
                     f'text-anchor="middle">{esc(text)}</text>')
        o.append('</g>')

    if opts.title:
        o.append(f'<text x="{pad_l:.1f}" y="{pad_t - 12:.1f}" font-size="14" '
                 f'fill="#222">{esc(opts.title)}</text>')
    if opts.legend and ykey == 's':
        o.append(legend_svg(pad_l, H - 9, char_p=getattr(opts, 'char_p', 0)))
    o.append(f'<text x="{W - pad_r:.1f}" y="{H - 8:.1f}" font-size="10" '
             f'text-anchor="end" fill="{C_TEXT}">stem  t-s</text>')
    o.append('</svg>')
    return '\n'.join(o)


def legend_svg(x, y, r=2.6, gap=2.2, char_p=0):
    """One line explaining the glyphs, drawn under the x-axis.

    Over a base ring of characteristic p every summand is Z/p, so the box and
    the rings cannot occur and advertising them would be misleading.
    """
    parts, cx = [], x + r
    def label(text, cx, pad=0.0):
        parts.append(f'<text x="{cx + r + 3 + pad:.1f}" y="{y + 3:.1f}" '
                     f'font-size="8" fill="{C_TEXT}">{text}</text>')
        return cx + r + 3 + pad + 5.2 * len(text) + 16
    if char_p:
        parts.append(f'<circle cx="{cx:.1f}" cy="{y:.1f}" r="{r:.1f}" '
                     f'fill="{C_DOT}" stroke="none"/>')
        label(f'ℤ/{char_p}  (the base ring is an F_{char_p}-algebra, so every '
              f'summand is ℤ/{char_p})', cx)
        return '<g>' + ''.join(parts) + '</g>'
    parts.append(f'<rect x="{cx - r:.1f}" y="{y - r:.1f}" width="{2*r:.1f}" '
                 f'height="{2*r:.1f}" fill="{C_DOT}" stroke="{C_EDGE}" '
                 f'stroke-width="0.6"/>')
    cx = label('ℤ_(3)', cx) + r
    parts.append(f'<circle cx="{cx:.1f}" cy="{y:.1f}" r="{r:.1f}" '
                 f'fill="{C_DOT}" stroke="none"/>')
    cx = label('ℤ/3', cx) + r + gap
    parts.append(f'<circle cx="{cx:.1f}" cy="{y:.1f}" r="{r:.1f}" '
                 f'fill="{C_DOT}" stroke="none"/>')
    parts.append(f'<circle cx="{cx:.1f}" cy="{y:.1f}" r="{r + gap:.1f}" '
                 f'fill="none" stroke="{C_EDGE}" stroke-width="0.9"/>')
    cx = label('ℤ/3ⁿ  (n-1 rings)', cx, pad=gap) + r + gap
    parts.append(f'<circle cx="{cx:.1f}" cy="{y:.1f}" r="{r:.1f}" '
                 f'fill="{C_DOT}" stroke="none"/>')
    parts.append(f'<circle cx="{cx:.1f}" cy="{y:.1f}" r="{r + gap:.1f}" '
                 f'fill="none" stroke="{C_EDGE}" stroke-width="0.9" '
                 f'stroke-dasharray="1.6,1.3"/>')
    label('order cut off by the range (lower bound)', cx, pad=gap)
    return '<g>' + ''.join(parts) + '</g>'


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description='Chart the Adams-Novikov E2 page from mr_BP output.',
        epilog='Run in the directory containing a finished mr_BP run.')
    ap.add_argument('halfT', help="mr_BP's first argument, e.g. 25")
    ap.add_argument('-d', '--dir', default='.', help='where the tables live')
    ap.add_argument('-c', '--comodule', default='',
                    help="chart an ./mr_BP_comod run for this comodule (e.g. "
                         "alpha_1), which writes <halfT>_<comodule>BP... "
                         "instead of <halfT>_BP.... Default: plain mr_BP output.")
    ap.add_argument('-o', '--output', help='default <halfT>_anss_E2.svg')
    ap.add_argument('--grading', choices=['anss', 'algnov'], default='anss',
                    help='anss (default): (t-s, s), one mark per cyclic '
                         'summand. algnov: mr_BP\'s printed (t-s, s+i), one '
                         'dot per class, with differentials and the '
                         'multiplication-by-p (or by v_n) lines.')
    ap.add_argument('--no-labels', dest='labels', action='store_false',
                    help='omit the class names next to the marks')
    ap.add_argument('--no-legend', dest='legend', action='store_false',
                    help='omit the box/dot/circles key under the axis')
    ap.add_argument('--no-bound', dest='bound', action='store_false',
                    help='omit the dashed t = <halfT> line marking where the '
                         'degree bound cuts the page off')
    ap.add_argument('--label-size', type=float, default=5.5,
                    help='font size for the class names (default 5.5)')
    ap.add_argument('--raw-labels', action='store_true',
                    help='label with the table text verbatim (v1^3[1-0]) '
                         'instead of typesetting it as v_1^3[1-0]')
    ap.add_argument('--max-rings', type=int, default=4,
                    help='draw at most this many rings; beyond that the '
                         'order is written beside the mark (default 4)')
    ap.add_argument('--ring-gap', type=float, default=2.0,
                    help='px between concentric rings (default 2.0)')
    ap.add_argument('--omit-stem0', action='store_true',
                    help='drop stem 0, as the published charts do')
    ap.add_argument('--height', type=int, default=None,
                    help='the height n of the run, i.e. the base ring is '
                         'BP_*/I_n (0 = BP_* itself). Normally read from '
                         '<prefix>run_info.txt; pass it only for output '
                         'produced before runs recorded it.')
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
    if a.bound:                                  # the degree bound to draw
        try:
            a.bound = int(a.halfT)
        except ValueError:
            a.bound = None

    base = os.path.join(a.dir, f'{a.halfT}_{a.comodule}BP')
    table = base + 'AANSS_table.txt'
    if not os.path.exists(table) and not a.comodule:
        # With no -c the prefix is <halfT>_BP, which is what a plain mr_BP run
        # writes. But `mr_BP_comod <halfT> <s>` with its comodule argument
        # omitted also computes the sphere, and writes <halfT>_sphereBP
        # instead -- so fall back to that rather than making the two tools
        # disagree about what "no comodule named" means.
        sphere = os.path.join(a.dir, f'{a.halfT}_sphereBP')
        if os.path.exists(sphere + 'AANSS_table.txt'):
            a.comodule, base = 'sphere', sphere
            table = base + 'AANSS_table.txt'
            sys.stderr.write(f'note: no {a.halfT}_BP... run here; charting '
                             f'the mr_BP_comod sphere run instead\n')
    if not os.path.exists(table):
        if a.comodule:
            sys.exit(f'error: {table} not found.\n'
                     f'Run ./BPtab {a.halfT} && '
                     f'./mr_BP_comod {a.halfT} <s> {a.comodule} first, '
                     f'or pass --dir.')
        sys.exit(f'error: {table} not found, and no '
                 f'{a.halfT}_sphereBPAANSS_table.txt to fall back on.\n'
                 f'Run ./BPtab {a.halfT} && ./mr_BP_comod {a.halfT} <s>, or '
                 f'./mr_st {a.halfT} <s+1> && ./BPtab {a.halfT} && '
                 f'./mr_BP {a.halfT} <s> first, or pass --dir.')

    classes, diffs = parse_table(table)
    h0 = parse_mult(base + 'AANSS_h0.txt')

    # Which ring was this run done over?  It decides what the tables MEAN --
    # over BP_* a filtration-0 class generates a Z_(3), over BP_*/I_n a Z/p --
    # so it is read from the run rather than assumed.
    info = parse_run_info(base + 'run_info.txt')
    if a.height is not None:                     # explicit override wins
        height = a.height
        char_p = PRIME if height > 0 else 0
    elif info is not None and 'height' in info:
        height = int(info['height'])
        # the run publishes the characteristic as well, and that is the fact
        # that actually matters here, so prefer it over re-deriving it
        char_p = int(info.get('characteristic', PRIME if height > 0 else 0))
    else:
        # No statement of the ring, so the height has to be assumed. That is
        # safe for output produced before runs recorded it -- the height
        # feature did not exist, so every such run really was over BP_* --
        # but NOT for a height-n table that has been separated from its
        # run_info.txt, where assuming BP_* draws filtration 0 as Z_(3) boxes
        # and hedges on bidegrees holding several summands, both wrong. The
        # table says which case this is, so check before assuming.
        if looks_characteristic_p(classes):
            sys.exit(
                f'error: cannot tell which ring this run was done over.\n'
                f'  {base}run_info.txt is missing, and no class in\n'
                f'  {table}\n'
                f'  carries a v0 -- which is what a run over BP_*/I_n looks '
                f'like, since p is 0 there.\n'
                f'Assuming BP_* would draw filtration 0 as Z_(3) boxes and '
                f'hedge on the bidegrees\nholding several summands, and both '
                f'would be wrong for such a run.\n'
                f'Pass --height n (0 for BP_*, 1 for BP_*/p, 2 for '
                f'BP_*/(p,v_1), ...), or put the\nrun_info.txt back beside '
                f'the tables.')
        height, char_p = 0, 0
        sys.stderr.write(
            f'note: {base}run_info.txt not found; assuming a run over BP_* '
            f'itself (height 0). Pass --height if it is not.\n')

    # multiplication by the bottom generator of the base ring's maximal
    # invariant ideal: p = v_0 over BP_*, v_n over BP_*/I_n
    a0_name = f'AANSS_a{height}.txt'

    if a.grading == 'anss':
        a.ykey = 's'
        if char_p:
            # p is 0 in the base ring, so every Ext group is an F_p-vector
            # space and every class is its own Z/p summand.  That is a fact
            # about the ring, not a shrug at a missing table: the a<n> table
            # records multiplication by v_n, which does not decompose Ext into
            # cyclic summands the way multiplication by p does.
            a0 = None
        else:
            a0 = parse_mult_detail(base + a0_name)
            if a0 is None:
                sys.stderr.write(
                    f'warning: {base}{a0_name} not found, so the 3-tower '
                    f'structure is unknown; every class is drawn as its own '
                    f'Z/3, which UNDERSTATES any Z_(3) or Z/9 above '
                    f'filtration 0\n')
        a.char_p = char_p                        # render() reads this too
        marks, owner = build_marks(classes, a0, ext0_torsion_free=not char_p)
        # A bidegree holding more than one chain is the one case where the
        # data does not settle the group: multiplication by 3 is recorded
        # only up to algebraic Novikov filtration, so a product that jumps
        # filtration could join two of the chains drawn here.
        #
        # That argument needs p to be a non-zero element to multiply BY. In
        # characteristic p it is zero, the whole cobar complex -- hence every
        # Ext group -- is an F_p-vector space, and several summands in a
        # bidegree just means dimension > 1: no hidden extension can merge
        # them into a Z/p^2. Flagging them would hedge on the one case that
        # IS proved, so above height 0 the count stays at 1.
        per_cell = defaultdict(int)
        for m in marks:
            per_cell[(m['stem'], m['s'])] += 1
        for m in marks:
            m['cohabitants'] = 1 if char_p else per_cell[(m['stem'], m['s'])]
        diffs, towers = [], []
        # h0-multiplication between summands: redirect each class to its mark
        struct, seen = [], set()
        for x, y in h0:
            if x in owner and y in owner and owner[x] != owner[y]:
                e = (marks[owner[x]]['name'], marks[owner[y]]['name'])
                if e not in seen:
                    seen.add(e)
                    struct.append(e)
    else:
        a.ykey = 'algnov'
        # same file, per height: multiplication by p, or by v_n
        a0 = parse_mult(base + a0_name)
        marks = [dict(m, kind='dot', n=1, group='', members=[m['name']])
                 for m in classes.values()]
        names = {m['name'] for m in marks}
        towers = [(x, y) for x, y in a0 if x in names and y in names]
        struct = [(x, y) for x, y in h0 if x in names and y in names]

    if a.omit_stem0:
        marks = [m for m in marks if m['stem'] != 0]
    if a.max_stem is not None:
        marks = [m for m in marks if m['stem'] <= a.max_stem]
    if a.max_filt is not None:
        marks = [m for m in marks if m[a.ykey] <= a.max_filt]
    keep = {m['name'] for m in marks}
    struct = [e for e in struct if e[0] in keep and e[1] in keep]
    towers = [e for e in towers if e[0] in keep and e[1] in keep]
    diffs = [d for d in diffs if d[0] in keep and d[1] in keep]

    if not marks:
        sys.exit('error: nothing to draw (is the run empty, or the crop too tight?)')

    if a.no_title:
        a.title = None
    elif a.title is None:
        grading = ('Adams-Novikov E2, p=3' if a.grading == 'anss'
                   else 'algebraic Novikov table, p=3')
        subject = f' of {a.comodule}' if a.comodule else ''
        a.title = (f'{grading}{subject}   '
                   f'(halfT={a.halfT}, {len(marks)} '
                   f'{"summands" if a.grading == "anss" else "classes"})')

    svg = render(marks, diffs, struct, towers, a)
    out = a.output or os.path.join(a.dir, f'{a.halfT}_{a.comodule}anss_E2.svg')
    with open(out, 'w') as fh:
        fh.write(svg)

    stems = [m['stem'] for m in marks]
    if a.grading == 'anss':
        kinds = defaultdict(int)
        for m in marks:
            kinds[m['kind']] += 1
        trunc = sum(1 for m in marks if m['truncated'])
        ring = f'BP_*/I_{height}' if height else 'BP_*'
        if char_p:
            breakdown = (f'{kinds["dot"]} Z/{char_p} (dot) -- every summand, '
                         f'since the base ring is an F_{char_p}-algebra')
        else:
            breakdown = (f'{kinds["box"]} Z_(3) (box), {kinds["dot"]} Z/3 '
                         f'(dot), {kinds["circles"]} Z/3^n (rings), of which '
                         f'{trunc} have the top of the tower cut off by the '
                         f'range (dashed)')
        print(f'{len(marks)} summands over {ring} '
              f'(stems {min(stems)}-{max(stems)}): {breakdown}; '
              f'{len(struct)} alpha_1 lines -> {out}')
        cells = defaultdict(list)
        for m in marks:
            cells[(m['stem'], m['s'])].append(m)
        if char_p:
            # Both of the caveats below are about reading orders off the a0
            # chain. In characteristic p there is no chain to read: every Ext
            # group is an F_p-vector space, so the class count settles it and
            # every bidegree is exact.
            print(f'groups: all {len(cells)} bidegree(s) pinned down exactly '
                  f'-- over an F_{char_p}-algebra the class count IS the '
                  f'group, see CHARTS.md section 2')
        else:
            shared = sum(1 for ms in cells.values() if len(ms) > 1)
            cut = sum(1 for ms in cells.values()
                      if len(ms) == 1 and ms[0]['truncated'])
            print(f'groups: {len(cells) - shared - cut} bidegree(s) pinned down '
                  f'exactly (a single 3-tower, terminating inside the range), '
                  f'{cut} known only as a lower bound (tower truncated), '
                  f'{shared} holding more than one summand, where a hidden '
                  f'extension could merge them -- see CHARTS.md section 2')
    else:
        print(f'{len(marks)} classes (stems {min(stems)}-{max(stems)}), '
              f'{len(struct)} alpha_1 lines, {len(towers)} '
              f'{"v_" + str(height) if height else "3"}-multiplication lines, '
              f'{len(diffs)} differentials -> {out}')
    print(f'note: this run computes internal degree t = stem + s up to '
          f'{a.halfT}, so the page is cut off along a diagonal (drawn '
          f'dashed), not at the top of either axis. Filtration is capped '
          f'separately by the resolution length: s + i < that argument.')


if __name__ == '__main__':
    main()
