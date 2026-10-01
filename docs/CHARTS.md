# Charting the E2 page

`anss_chart.py` turns a finished `mr_BP` run into an SVG chart of the
Adams–Novikov E2 page. It reads only the text tables `mr_BP` already writes,
so it never recomputes anything and needs no rebuild — you can point it at
old output, including output from other machines.

```sh
./mr_st 35 31 && ./BPtab 35 && ./mr_BP 35 30    # the usual pipeline
./anss_chart.py 35                              # -> 35_anss_E2.svg
```

Requires Python 3 and nothing else. The argument is the same `<halfT>` you
gave the pipeline; that's how it finds the files. Run it in the directory
holding the run, or pass `--dir`.

### Charting a comodule other than the sphere

`mr_BP_comod` (see [`GENERAL_COMODULES.md`](GENERAL_COMODULES.md)) writes its
tables with a `<halfT>_<comodule>BP...` prefix rather than `<halfT>_BP...`.
Pass `-c/--comodule` to chart one:

```sh
./BPtab 20 && ./mr_BP_comod 20 4 alpha_1   # resolve S/alpha_1
./anss_chart.py 20 -c alpha_1              # -> 20_alpha_1anss_E2.svg
```

The comodule name lands in the output filename and the chart title, so
charts for different complexes don't overwrite each other. Nothing in the
script is specific to the sphere: the glyphs, groups and structure lines
below are all derived from whichever run's tables it is pointed at.

With no `-c` the prefix is `<halfT>_BP...`, what a plain `mr_BP` run writes.
Since `mr_BP_comod` also defaults to the sphere but writes
`<halfT>_sphereBP...`, the script falls back to that prefix when there is no
`<halfT>_BP...` run in the directory, and says so on stderr. The fallback is
exactly equivalent to passing `-c sphere`, output filename and title
included; an explicit `-c` never falls back.

One caveat when reading such a chart: `mr_BP_comod` runs `mult_table()` (so
the `α₁` structure lines are drawn) but deliberately **not** `mult_theta()`,
which is specific to the Moore spectrum and carries a hardcoded table of
theta degrees. The `theta_i` tables are Bockstein-side and aren't drawn on an
ANSS chart anyway (see §4), so this costs nothing here.

## 1. Grading convention

The default (`--grading anss`) plots the Adams–Novikov bidegree:

- **x** = stem = `t - s`
- **y** = `s`, the homological degree

**This is not the pair `mr_BP` prints.** `algNov.cpp:118` emits
`|deg=(t-s, s+i)`, where `i` is the algebraic Novikov filtration, so
`v1^1[1-0]` (which is α₂, in stem 7 of Ext¹) is printed at height **2** and
belongs at height **1**. The chart takes the stem from the printed pair but
reads `s` off the `[s-n]` bracket in the class name instead. The internal
degree is recoverable as `t = x + y`.

`--grading algnov` plots `mr_BP`'s printed pair as-is, one dot per class,
draws the algebraic Novikov differentials, and joins each class to its
multiples by the bottom generator of the base ring's maximal invariant ideal
— `p = v₀` over `BP_*`, which makes vertical towers, or `v_n` over
`BP_*/I_n`, which does not (`v_n` raises the stem as well as the
filtration). That view is useful for checking a run against the raw tables;
it is not an ANSS chart.

### What a class name means

A name like `v0^1v1^3[1-1]` has two parts:

```
<monomial in the base ring>[<homological degree s>-<index of the generator of V_s>]
```

For a minimal resolution `0 → M → F_0 → F_1 → …` with `F_s = Γ ⊗_A V_s`, the
complex of primitives is `V_0 → V_1 → …`, with each `V_s` free over the base
ring `A` — `BP_*`, or `BP_*/I_n` for a height `n` run. So `[s-g]` names
generator `g` of `V_s`, and the monomial is that generator's coefficient in
`A`. `algNov.cpp`'s `naming`/`output` build it: slot 1 of a name is the
`v₀`-exponent and slots 2… are the `v₁…v₅` exponents, printed as `v<i-1>^e`.
Note `v₀` is `p` — it is the `p`-adic valuation of the coefficient
(`algNov.cpp:28`), not a polynomial generator — so no name carries a `v0` at
height ≥ 1, where `p` is 0.

For `mod_p` the zero line therefore reads `[0-0]`, `v1^1[0-0]`, `v1^2[0-0]`,
…: `BP_*/p` has rank 1 in degree 0, so `V_0` has a single generator,
corresponding to `1 ∈ BP_*/p`, and `v1^n[0-0]` is `v₁ⁿ · 1` in stem `4n`
(`|v₁| = 4`). That is right on the nose, since `Ext⁰ = Prim(M)` and the
invariants of `BP_*/p` are exactly `F₃[v₁]`. The sphere's zero line is the
same `[0-0]` with prefix `v0^j`, spelling out `Ext⁰(BP_*) = Z₍₃₎`.

**The monomial is a leading term, not the class.** The tables list a basis of
`gr Ext` for the algebraic Novikov filtration, so in general `v1^1[2-0]`
means "leading term `v₁` times generator 0 of `V_2`" — a statement about the
associated graded, which is why §2 warns against reading `v0` in a name as
"divisible by 3". In filtration 0 the distinction collapses whenever `Ext⁰`
consists of honest invariants of `A`, as it does for `BP_*/p`. The `g` in
`[s-g]` is just a position in the generator list this particular resolution
chose; it carries no further meaning.

## 2. One mark per cyclic summand

The classes in the table are an F₃-basis of `gr Ext`, so a `Z/9` appears as
*two* of them. Drawing one dot each would misrepresent the group, so the
ANSS view groups them: multiplication by 3 is exactly what
`<prefix>AANSS_a0.txt` records, and chaining those edges recovers each cyclic
summand. A class that is not 3 times any *surviving* class generates a
summand; the chain of its 3-multiples gives the order.

The glyphs follow Belmont's published 3-primary ANSS chart:

| glyph | meaning |
|---|---|
| filled square, marked ∞ | `Z_(3)` |
| filled dot | `Z/3` |
| dot inside `n-1` concentric rings | `Z/3^n` |
| dashed outermost ring | the tower runs off the end of the computed range, so the order drawn is a **lower bound** |

Several summands in one bidegree are drawn side by side. Each is labelled
with the name of its generator, typeset the way the published charts do —
`v₀v₁³[1-1]` for the table's `v0^1v1^3[1-1]`. `--raw-labels` keeps the
table's exact text, `--no-labels` drops labels entirely, and `--no-legend`
drops the key under the axis.

Two rules decide the glyph, and they come from different places:

- **Filtration 0 is `Z_(3)` when the base ring has characteristic 0**, by
  mathematics rather than by the tables: over `BP_*` the comodule is free, so
  `Ext^0 = Prim(M)` is a submodule of a free `BP_*`-module, hence
  torsion-free. (This is why `S/α₁` has a box in stem 4: `Ext^{0,4}` there is
  `Z_(3){v_1x_0 + 3x_4}` — see [`SES_CHECK.md`](SES_CHECK.md) §4.)
- **Above filtration 0** the order comes from the `a0` chain. Where that
  chain hits the truncation the mark gets the dashed ring: the printed
  second coordinate `s + i` is cut off at the resolution length (§6), and the
  `a0` table is pruned one step earlier still (`multiplication.cpp:169`), so
  near the top of the filtration range the top of a tower is simply not
  visible. At `halfT=185, L=14`, 60 of 344 summands are in that situation.

### Over `BP_*/I_n` every summand is `Z/3`

Both rules above are about `BP_*`. A run at height `n ≥ 1` (see
[`GENERAL_COMODULES.md`](GENERAL_COMODULES.md)) is over `BP_*/I_n`, which is
an `F_3`-algebra: `p` kills the comodule, so **every** `Ext` group —
filtration 0 included — is an `F_3`-vector space and every summand is a
`Z/3` dot. There are no boxes and no rings, and the legend says so.

This also settles the group in *every* bidegree, not just filtration 0, so
the uncertainty discussed under
["What the glyph proves, and what it assumes"](#what-the-glyph-proves-and-what-it-assumes)
below does not arise: several dots in one bidegree mean dimension `> 1`, and
no hidden extension can merge them into a `Z/9`. The run's summary line says
so, and no tooltip hedges.

The script does not guess this. Each run writes `<prefix>run_info.txt`:

```
# what this run computed. Written by BPInit; read by anss_chart.py.
height 1
characteristic 3
max_degree 30
resolution_length 8
```

and the chart reads the characteristic from it. Output produced before runs
recorded this has no such file; the script then assumes height 0, says so,
and `--height n` overrides it. Getting this wrong is silent rather than an
error — asserting `Z_(3)` unconditionally is what once drew the
`Ext^0 = F_3[v_1]` of `S/p` as a row of boxes — which is why it is published
by the run rather than inferred from the tables.

That assumption is sound for output predating `run_info.txt`, since the
height did not exist then, but not for a height-`n` table that has been
*separated* from its `run_info.txt` — copying one into an examples directory
is enough. So before assuming, the script cross-checks the table against
itself: `v0` in a class name is the `p`-adic valuation of its coefficient, so
a run over `BP_*` is full of them (`Ext^{0,0} = Z_(3)` alone contributes a
whole `v_0`-tower; in practice about half the classes carry one) and a run
over `BP_*/I_n` cannot produce a single one. If the metadata is missing *and*
no class carries a `v0`, the script refuses to draw rather than guess:

```
error: cannot tell which ring this run was done over.
  …run_info.txt is missing, and no class in
  …AANSS_table.txt
  carries a v0 -- which is what a run over BP_*/I_n looks like, since p is 0 there.
```

It is evidence, not proof — a height-0 run truncated before any `v_0`
appeared would look the same — which is why it stops and asks for `--height`
instead of flipping the assumption on its own.

One consequence for the `a0` chain: multiplication by `p` is the zero map at
height `n ≥ 1`, so there is nothing to chain, and the summands are singletons
by mathematics rather than for want of a table. What a height-`n` run writes
as `<prefix>AANSS_a<n>.txt` is multiplication by `v_n` — the bottom generator
of the base ring's maximal invariant ideal, the true analogue of `a0` — and
that is structure data, not summand data: the `--grading algnov` view draws
those lines, and the ANSS view does not use them.

A class killed by an algebraic Novikov differential is dropped, along with
the tag that killed it: lines containing `<-` record a differential and both
ends die. (The tag never gets a line of its own — `SS_table::output`,
`SS.h:139-146`, only emits untagged entries — so dropping the `<-` lines
removes both.)

Note that 3-multiples must be found from the `a0` table, **not** by looking
for `v0` in the name. `v0^1[1-1]` carries a `v0` but is not 3 times anything
that survives — its predecessor `[1-1]` supports a d₂ — and it is the
generator of the `Z/9` in stem 11. A syntactic filter would hide it.

If the `a0` file is missing in a characteristic-0 run the script warns and
draws every class as a plain dot, since without it no isomorphism type is
knowable — which understates any `Z_(3)` or `Z/9` above filtration 0.

### What the glyph proves, and what it assumes

All of this concerns a run over `BP_*`. Over `BP_*/I_n` there is no `a0`
chain to read — `p` is `0` there — and the class count settles every group on
its own, so none of the caveats below apply; see
[Over `BP_*/I_n` every summand is `Z/3`](#over-bp_i_n-every-summand-is-z3)
above.

The two halves of the evidence behave differently, so it is worth being
precise about which conclusions are theorems:

- **A non-zero `a0` entry is sound.** It says the leading term of `3x` is
  non-zero, hence `3x ≠ 0` in `Ext`. So a chain of length `n` exhibits an
  element of order `3^n`: a **lower** bound.
- **A `-> o` entry is not.** It says the leading term of `3x` vanishes, which
  means either `3x = 0` or `3x` has jumped further up the algebraic Novikov
  filtration than the table can see — a hidden extension.
- **The class count is exact.** The surviving classes are an F₃-basis of
  `gr Ext`, so a bidegree with `n` of them has length exactly `n`.

Put together: **when one chain accounts for every class in a bidegree, the
group is pinned down**, because a group of length `n` containing an element
of order `3^n` is cyclic. That is how `S/α₁` gets `Z/9` in (7,1) and `Z/27`
in (11,1) rather than `Z/3 ⊕ Z/3` and `Z/9 ⊕ Z/3` — those have no element of
the required order, and the `a0` chain exhibits one.

When a bidegree holds **more than one** chain, the splitting drawn is only
the finest possibility: a hidden extension could merge two of the marks into
one bigger cyclic group. The script says how many bidegrees are in each
state after every run, and the tooltip of any mark sharing its bidegree says
so explicitly. At `halfT=35` all 15 of the sphere's bidegrees and all 14 of
`S/α₁`'s are pinned down exactly; at `halfT=185, L=14` the sphere has 168
pinned down, 37 known only as lower bounds (truncated towers) and 64 holding
several summands.

The repo does contain an independent record of this: the Bockstein tables
`<prefix>BocSS_table.txt` are the `v_0`-Bockstein spectral sequence, whose
`d_r` differentials are exactly the `p`-divisibility data, with
`<prefix>B2A_table.txt` translating Bockstein names to algebraic Novikov
ones. For `S/α₁` those tables confirm both orders above — `v0^2[1-0]` dies by
a `d2` and `v0^3[1-1]` by a `d3` — and they are the natural place to look for
the ambiguous bidegrees. The chart does not read them, for the grading
reason in §4.

## 3. Structure lines

Solid tan lines are multiplication by `h0 = (η_R(v1) - η_L(v1))/p` = α₁,
read from `<prefix>AANSS_h0.txt`, so they run `(+3, +1)` — slope 1/3, as in
the published charts. Products between classes inside the same summand are
suppressed; each line joins two marks.

Multiplication by α₂ (`(+7, +1)` lines) and β₁-divisibility (which Belmont's
chart indicates by colouring the classes) are not currently available:
`BPInit::mult_table` (`BP_init.cpp:180-187`) only computes the `h0` table for
the algebraic Novikov side, and the `theta_i` tables `mult_theta` writes are
Bockstein-side. Adding either means calling `mult_table(<class>, <degree>,
"<name>.txt")` with the appropriate `BPBP` element (`BP_oper.thetas()` has
β₁) and teaching this script the extra file — the parsing side already
handles any multiplication table. Massey products, which that chart also
draws, are out of reach here entirely.

## 4. Bockstein tables are not supported

`Boc.cpp:37-38` computes `(degree*2 - k, k + filtration)` but prints only the
first coordinate, and that coordinate is `2t - s`, not the stem — for α₁ it
prints `7` where the algebraic Novikov table prints stem `3`. So the
Bockstein `.txt` files are missing the second coordinate entirely and use a
different first one. Charting them needs either a fix at that line or a
reconstruction from the class names (`Boc_table::v_valuation`, `Boc.cpp:4-6`,
counts only powers of `v0`, so `y = s + v0-exponent`).

## 5. Verification

At `halfT=35, s=30` the chart's marks are

```
(0,0) box        (3,1) (7,1) (10,2) (11,1)°° (13,3) (15,1) (19,1) (20,4)
(23,1)°° (23,5) (26,2) (27,1) (29,3) (31,1)
```

i.e. 1, α₁, α₂, β₁, α_{3/2}, α₁β₁, α₄, α₅, β₁², α_{6/2}, α₁β₁², β₂, α₇,
α₁β₂, α₈ — agreeing with the published p=3 ANSS charts on every stem in
range. The two double-ringed marks (°°) are stems 11 and 23, which are
exactly the two α-family classes of order 9 in that range: `α_{i/j}` has
order `3^j`, and `j = 2` precisely when `3 | i`.

For `S/α₁` the same chart has boxes in stems 0 **and 4** and nothing in
stems 3, 13, 23, 29 — see [`SES_CHECK.md`](SES_CHECK.md), which checks the
whole page against the sphere's automatically.

### Worked examples you can re-run

`example_data/` holds three finished runs at `t = 30, length = 8` — the
sphere, `S/α₁` and `S/p` — each as the tables the chart reads plus the SVG it
produced. They are complete, so each redraws itself byte for byte:

```sh
./anss_chart.py 30 -c mod_p -d example_data -o /tmp/check.svg
cmp /tmp/check.svg example_data/30_mod_panss_E2.svg    # identical
```

That makes them a regression test on the chart as well as an illustration: a
change that alters any of the three shows up as a `cmp` failure. Each run
contributes `AANSS_table.txt` (the classes), `AANSS_h0.txt` (α₁ lines),
`AANSS_a0.txt` or `AANSS_a1.txt` (the summand structure — `a<n>` per the
height) and `run_info.txt` (the ring). Drop any one of them and the redraw
degrades: without `h0` the α₁ lines vanish, without `a0` the towers are
understated, and without `run_info.txt` the `S/p` example refuses to draw at
all rather than guess (§2).

## 6. Truncation: the page is cut off along a diagonal

Neither axis is what limits the picture, and the y-axis is not capped at all
— it grows to fit whatever the run produced. **What cuts the page off is the
degree bound, and since `t = stem + s` that boundary is a diagonal.** The
chart draws it as a dashed line labelled `t = <halfT>` (`--no-bound` to hide
it); everything beyond it is missing because it was never computed.

This is worth internalising, because it is not where the eye expects the
edge. At `halfT=35` the sphere has a class in stem **31** (`v1^7[1-0]`, `s=1`,
so `t = 32 ≤ 35`) but *not* in stem **30** — β₁³ = `[6-0]` lives at `s = 6`,
so `t = 36 > 35`. Re-running with `./mr_BP_comod 40 12 sphere` produces it,
along with α₁β₁³ = `[7-0]` at `(33,7)`, and the y-axis extends to 7 by
itself. So:

- a missing class high up the chart usually means **raise the first
  argument** (the internal-degree bound), not the second;
- the second argument still matters independently: filtration is capped by
  `s + i <` resolution length, which is what the dashed rings of §2 are
  about, and a class at filtration `s` needs a resolution at least that long.

Near the boundary `mr_BP` also prints things like `out of range for beta1`
when a needed class falls outside the budget, so treat the last stems before
the line as unreliable too, and crop with `--max-stem` when showing the chart
to anyone.

## 7. Options

| Flag | Effect |
|---|---|
| `-d, --dir` | directory holding the tables (default `.`) |
| `-c, --comodule` | chart an `mr_BP_comod` run for this comodule (e.g. `alpha_1`); default is a plain `mr_BP` run |
| `-o, --output` | output file (default `<halfT>_<comodule>anss_E2.svg`) |
| `--grading anss\|algnov` | see §1 |
| `--no-labels` | drop the class names |
| `--raw-labels` | label with the table text verbatim (`v1^3[1-0]`) |
| `--label-size` | font size for the names (default 5.5) |
| `--no-legend` | drop the glyph key under the axis |
| `--no-bound` | drop the dashed `t = <halfT>` degree-bound diagonal (§6) |
| `--max-rings` | most rings drawn before the order is written beside the mark (default 4) |
| `--ring-gap` | px between rings (default 2.0) |
| `--omit-stem0` | drop `Ext^0` in stem 0, as the published charts do |
| `--max-stem`, `--max-filt` | crop |
| `--unit` | px per lattice step (default 31.2, matching published charts) |
| `--dot-radius`, `--tick`, `--title`, `--no-title` | cosmetics |

Every mark carries an SVG `<title>`, so hovering it in a browser shows the
class name, the group, and both gradings.
