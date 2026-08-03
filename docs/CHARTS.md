# Charting the E2 page

`anss_chart.py` turns a finished `mr_BP` run into an SVG chart of the
Adams–Novikov E2 page. It reads only the text tables `mr_BP` already writes,
so it recomputes nothing, needs no rebuild, and can be pointed at old output —
including output produced on another machine.

```sh
./anss_chart.py 35            # -> 35_anss_E2.svg
```

The argument is the same `<halfT>` you gave the pipeline; that is how it finds
the files. Run it in the directory holding the run, or pass `--dir`. Python 3
is the only requirement — no packages to install.

Open the result in a browser: every dot carries an SVG `<title>`, so hovering
it shows the class name and both gradings.

## 1. Running the whole thing from scratch

### 1.1 Prerequisites

**Linux:**

```sh
sudo apt-get install g++ libgmp-dev python3
```

**macOS:** Apple's `clang` rejects `-fopenmp`, which every build script here
uses, so the stock compiler will not work:

```sh
brew install gcc gmp
```

Then edit `st_compiling`, `BPtable_compile` and `BP_compile` to call `g++-14`
(or whichever version Homebrew installed) instead of `g++`, and add
`-I/opt/homebrew/include -L/opt/homebrew/lib` to `BPtable_compile` so it finds
GMP. Nothing else needs changing. (Alternatively `brew install libomp` and
compile with `-Xpreprocessor -fopenmp -lomp`, but switching compilers is less
fiddly.)

### 1.2 Build

```sh
sh st_compiling        # -> mr_st
sh BPtable_compile     # -> BPtab
sh BP_compile          # -> mr_BP
```

The top-level README misspells the second one as `BPtable_complile`; the file
is `BPtable_compile`.

`st_compiling` builds with `-O0`, which makes `mr_st` far slower than it needs
to be. Changing it to `-O2` is safe.

### 1.3 Run

All four programs read and write files in the **current** directory, and the
repo has no `.gitignore`, so work in a scratch subdirectory rather than
dumping dozens of output files into your checkout:

```sh
mkdir -p run35 && cd run35
../mr_st 35 31
../BPtab 35
../mr_BP 35 30
../anss_chart.py 35 --omit-stem0      # -> 35_anss_E2.svg
```

Do a small run first — `20 15 / 20 / 20 14` finishes in seconds and confirms
the toolchain before you commit to a long one.

Three things that bite:

- `argv[1]` must be **identical** across all four commands, and `mr_st`'s
  second argument must be at least one more than `mr_BP`'s. Nothing validates
  this; a mismatch segfaults or allocates wildly rather than erroring cleanly
  (see [`BUILD_AND_RUN.md`](BUILD_AND_RUN.md) for the mechanism).
- Don't mix runs in one directory. The `<halfT>_` filename prefix is the only
  linkage between the programs, so stale files from a different `halfT` are
  picked up silently.
- Classes near the top stems are missing because of the degree bound, not
  because they are absent (§7).

## 2. Grading convention

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

## 3. Which classes are on the chart

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

## 4. Structure lines

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

## 5. Bockstein tables are not supported

`Boc.cpp:37-38` computes `(degree*2 - k, k + filtration)` but prints only the
first coordinate, and that coordinate is `2t - s`, not the stem — for α₁ it
prints `7` where the algebraic Novikov table prints stem `3`. So the
Bockstein `.txt` files are missing the second coordinate entirely and use a
different first one. Charting them needs either a fix at that line or a
reconstruction from the class names (`Boc_table::v_valuation`, `Boc.cpp:4-6`,
counts only powers of `v0`, so `y = s + v0-exponent`).

## 6. Verification

At `halfT=35, s=30` the chart's dots are

```
(3,1) (7,1) (10,2) (11,1) (13,3) (15,1) (19,1) (20,4) (23,1) (23,5)
(26,2) (27,1) (29,3) (31,1)
```

i.e. α₁, α₂, β₁, α₃, α₁β₁, α₄, α₅, β₁², α₆, α₁β₁², … — agreeing with the
published p=3 ANSS charts on every stem in range. That check exercises the
parsing, both filters in §3, and the regrading in §2 together, so it is the
one to re-run after changing any of them.

## 7. Truncation

Classes near the top stems are missing because of the degree bound, not
because they are absent — `mr_BP` prints things like `out of range for beta1`
when a needed class falls outside the budget. The chart cannot tell the
difference, so treat the right-hand edge as unreliable and crop it with
`--max-stem` when showing the chart to anyone.

## 8. Options

| Flag | Effect |
|---|---|
| `-d, --dir` | directory holding the tables (default `.`) |
| `-o, --output` | output file (default `<halfT>_anss_E2.svg`) |
| `--grading anss\|algnov` | see §2 |
| `--omit-stem0` | drop `Ext^0 = Z_(3)`, as the published charts do |
| `--max-stem`, `--max-filt` | crop |
| `--unit` | px per lattice step (default 31.2, matching published charts) |
| `--dot-radius`, `--tick`, `--title`, `--no-title` | cosmetics |

## 9. How the conventions were pinned down

The grading in §2 and the class-selection rules in §3 were not guessed. They
were derived by measuring a published p=3 ANSS chart — extracting its dot
positions and structure-line slopes — and then requiring that the tool
reproduce that dot set from `mr_BP`'s tables. Both §2 and §3 are corrections
that this comparison forced:

- plotting `mr_BP`'s printed second coordinate puts α₂ at height 2 instead of
  1, and every other `v`-multiple with it;
- filtering 3-multiples by name instead of by the `a0` table deletes α₃.

Either mistake produces a chart that looks plausible and is wrong, which is
why §6 exists.
