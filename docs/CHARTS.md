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
charts for different complexes don't overwrite each other. Without the flag
the script behaves exactly as before, reading a plain `mr_BP` run.

One caveat when reading such a chart: `mr_BP_comod` runs `mult_table()` (so
the `α₁` structure lines are drawn) but deliberately **not** `mult_theta()`,
which is specific to the Moore spectrum and carries a hardcoded table of
theta degrees. The `theta_i` tables are Bockstein-side and aren't drawn on an
ANSS chart anyway (see §3), so this costs nothing here.

## 1. Grading convention

The default (`--grading anss`) plots the Adams–Novikov bidegree:

- **x** = stem = `t - s`
- **y** = `s`, the homological degree
- one dot per generator; no 3-towers

**This is not the pair `mr_BP` prints.** `algNov.cpp:118` emits
`|deg=(t-s, s+i)`, where `i` is the algebraic Novikov filtration, so
`v1^1[1-0]` (which is α₂, in stem 7 of Ext¹) is printed at height **2** and
belongs at height **1**. The chart takes the stem from the printed pair but
reads `s` off the `[s-n]` bracket in the class name instead.

`--grading algnov` plots `mr_BP`'s printed pair as-is, keeps the 3-multiples,
draws them as vertical 3-towers, and draws the algebraic Novikov
differentials. That view is useful for checking a run against the raw tables;
it is not an ANSS chart.

## 2. Which classes are on the chart

The algebraic Novikov SS converges to `Ext_{BP_*BP} =` the ANSS E2 page, so a
class is drawn exactly when it survives that spectral sequence. Two filters
implement this:

1. **Lines containing `<-` record a differential**, and both ends die. The
   source never appears as a line of its own — `SS_table::output`
   (`SS.h:139-146`) only emits untagged entries — so dropping the `<-` lines
   removes both ends.
2. **Classes appearing as a target in `<halfT>_BPAANSS_a0.txt`** are 3 times
   another class, i.e. a rung on a 3-tower rather than a generator.

Filter 2 has to be driven by the multiplication-by-3 table, **not** by
looking for `v0` in the name. `v0^1[1-1]` carries a `v0` but is not 3 times
anything that survives — its predecessor `[1-1]` supports a d₂ — and it is
exactly α₃ in stem 11. A syntactic filter deletes it and puts a hole in the
chart.

## 3. Structure lines

Solid tan lines are multiplication by `h0 = (η_R(v1) - η_L(v1))/p` = α₁,
read from `<halfT>_BPAANSS_h0.txt`, so they run `(+3, +1)`.

Multiplication by α₂ (`(+7, +1)` lines, which published charts also draw) is
not currently available: `BPInit::mult_table` (`BP_init.cpp:180-187`) only
computes the `h0` table for the algebraic Novikov side. Adding it means
calling `mult_table(<class>, <degree>, "alpha2.txt")` with the appropriate
`BPBP` element, then teaching this script the extra file — the parsing side
already handles any multiplication table.

The `theta_i` tables that `mr_BP` writes are Bockstein-side
(`<halfT>_BPBocSS_theta*.txt`), so they are not drawn on an ANSS chart.

## 4. Bockstein tables are not supported

`Boc.cpp:37-38` computes `(degree*2 - k, k + filtration)` but prints only the
first coordinate, and that coordinate is `2t - s`, not the stem — for α₁ it
prints `7` where the algebraic Novikov table prints stem `3`. So the
Bockstein `.txt` files are missing the second coordinate entirely and use a
different first one. Charting them needs either a fix at that line or a
reconstruction from the class names (`Boc_table::v_valuation`, `Boc.cpp:4-6`,
counts only powers of `v0`, so `y = s + v0-exponent`).

## 5. Verification

At `halfT=35, s=30` the chart's dots are

```
(3,1) (7,1) (10,2) (11,1) (13,3) (15,1) (19,1) (20,4) (23,1) (23,5)
(26,2) (27,1) (29,3) (31,1)
```

i.e. α₁, α₂, β₁, α₃, α₁β₁, α₄, α₅, β₁², α₆, α₁β₁², … — agreeing with the
published p=3 ANSS charts on every stem in range.

## 6. Truncation

Classes near the top stems are missing because of the degree bound, not
because they are absent — `mr_BP` prints things like `out of range for beta1`
when a needed class falls outside the budget. The chart cannot tell the
difference, so treat the right-hand edge as unreliable and crop it with
`--max-stem` when showing the chart to anyone.

## 7. Options

| Flag | Effect |
|---|---|
| `-d, --dir` | directory holding the tables (default `.`) |
| `-c, --comodule` | chart an `mr_BP_comod` run for this comodule (e.g. `alpha_1`); default is a plain `mr_BP` run |
| `-o, --output` | output file (default `<halfT>_<comodule>anss_E2.svg`) |
| `--grading anss\|algnov` | see §1 |
| `--omit-stem0` | drop `Ext^0 = Z_(3)`, as the published charts do |
| `--max-stem`, `--max-filt` | crop |
| `--unit` | px per lattice step (default 31.2, matching published charts) |
| `--dot-radius`, `--tick`, `--title`, `--no-title` | cosmetics |

Every dot carries an SVG `<title>`, so hovering it in a browser shows the
class name and both gradings.
