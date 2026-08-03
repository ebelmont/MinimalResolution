# Build and Run

Concrete instructions for building and running every executable in the repo,
consolidated from the individual pipeline docs. See
[`ARCHITECTURE.md`](ARCHITECTURE.md) for how these fit together and the p=2
vs p=3 caveat (§5 there) before you pick which pipeline to run.

Requires GCC with C++11/C++14 support, OpenMP, and (for `BPtab` only) the GNU
Multiple Precision library (`libgmp`, `libgmpxx`). Charting the result with
`anss_chart.py` additionally needs Python 3 — see [`CHARTS.md`](CHARTS.md),
including the macOS notes, since Apple's `clang` rejects `-fopenmp`.

All executables take degree parameters as `argv[1]`/`argv[2]` and read/write
files in the **current working directory**, prefixed with the first
argument (e.g. `25_etaL`, `25_maps`). Run everything for a given computation
from the same directory, and don't mix output files from runs with different
degree arguments.

## Quick reference: build scripts → executables

| Script | Executable | Pipeline |
|---|---|---|
| `sh st_compiling` | `mr_st` | Steenrod |
| `sh kos_compile` | `kos` | Steenrod (Koszul variant) |
| `sh BPtable_compile` | `BPtab` | BP |
| `sh BP_compile` | `mr_BP` | BP |
| `sh ex_compile` | `mr_ex` | Ex/Ctau |
| `sh e2p_compile` | `e2p` | Ex/Ctau (helper for motivic) |
| `sh mottable_compile` | `motTab` | Motivic (helper) |
| `sh mot_compile` | `mr_mot` | Motivic |
| `sh mot_comb_compile` | `mot_comb` | Motivic |
| `sh tau_mult_compile` | `mot_mult` | Motivic |
| `sh tauboc_compile` | `tauBoc` | Motivic |

Build everything:
```sh
for s in st_compiling kos_compile BPtable_compile BP_compile ex_compile \
         e2p_compile mottable_compile mot_compile mot_comb_compile \
         tau_mult_compile tauboc_compile; do
  sh "$s"
done
```

## Pipeline 1: BP (the repo's main deliverable, p=3)

```sh
sh BPtable_compile   # -> BPtab
sh BP_compile        # -> mr_BP
sh st_compiling      # -> mr_st (produces the seed generators mr_BP needs)

./mr_st 25 21        # writes 25_gens_data, among other 25_-prefixed files
./BPtab 25           # writes 25_etaL, 25_R2L, 25_delta (structure maps)
./mr_BP 25 20        # loads all of the above; writes 25_BP*-prefixed output
./anss_chart.py 25   # draws 25_anss_E2.svg from the tables mr_BP just wrote
```

- `argv[1]` (here `25`) **must be the same value** across all three
  invocations — it's used directly as the max internal degree cutoff and
  determines the filenames each program looks for.
- The README calls `argv[1]` "half of `t`" — the max degree by which
  everything is truncated. Per `docs/pipelines/BP.md` §2.2, the *code itself*
  doesn't visibly halve anything again; treat this as the project's grading
  convention for choosing what value to pass, not a further transformation
  the programs perform.
- `mr_BP`'s second argument (`20` above) is the resolution length `s`; the
  README requires it to be **at least one less than** the `s` used for
  `mr_st`'s second argument (`mr_st`'s `21` above), matching the BP/I ↔ BP
  relationship described in `ARCHITECTURE.md` §4.
- **Getting the order or the degree argument wrong doesn't fail cleanly** —
  `docs/pipelines/BP.md` §2.2 traces this to an unchecked bounds read in the
  binary matrix loader (`matrices/4.h:30-37`↔`modules/5.h:40-49`): a
  mismatched `<halfT>` desyncs the file's read cursor from its record
  boundaries, and the resulting garbage record-length can trigger a huge
  allocation or an out-of-bounds read — this is the mechanism behind the
  README's "usually a break-down of the program such as a segmentation
  error" warning. There is no argument validation to catch this earlier.
- `mr_BP`'s output includes the actual algebraic-Novikov E2-page tables
  (`<halfT>_BPAANSS_table_binary`/`.txt`), the Bockstein-SS tables
  (`<halfT>_BPBocSS_table_binary`/`.txt`), and multiplicative structure
  (`h0`, `theta_i` tables) — see `docs/pipelines/BP.md` §4 for the full
  output-file inventory and what each one is.
- `anss_chart.py` post-processes those tables into an SVG chart and needs no
  rebuild, so it can be re-run over an existing run's output as often as you
  like. It plots `(t-s, s)`, which is **not** the pair the tables print —
  [`CHARTS.md`](CHARTS.md) §2 explains the fix-up and §3 which classes are
  drawn.

## Pipeline 2: Classical Steenrod Ext (p=3, standalone)

```sh
sh st_compiling      # -> mr_st
./mr_st 25 21        # <halfT> <resolution_length>
```

Self-contained; computes `Ext_{A_*}(F_3,F_3)` (a truncation: only the
polynomial part of the dual Steenrod algebra, `xi_1..xi_5`, no exterior
generators — see `docs/pipelines/STEENROD.md` §1). Also the generator seed
for the BP pipeline above (`ARCHITECTURE.md` §4).

## Pipeline 3: Koszul cross-check (p=2 hardcoded — likely legacy)

```sh
sh kos_compile       # -> kos
./kos 25 21          # same argument shape as mr_st
```

Reuses the classical Steenrod machinery but resolves a different ("doubled")
comodule and self-checks with `check_splitting`. **Hardcoded to p=2**
(`kosul.cpp:47`) regardless of arguments — see `docs/pipelines/STEENROD.md`
§5/§6 for why this looks like unported p=2 legacy code rather than a
p=3-valid cross-check.

## Pipeline 4: Ex/Ctau (p=2 hardcoded)

```sh
sh ex_compile        # -> mr_ex
./mr_ex 25 21        # <maxdeg> <resolution_length>
```

Self-generates its coproduct table on first run. Outputs (`_gens_data_ctau`,
`_extables`) are consumed by the motivic pipeline below. **Hardcoded to
p=2** (`ex_main.cpp:11`) — see `docs/pipelines/EX_CTAU.md` §1/§6.

## Pipeline 5: Motivic (p=2, tau-graded)

Run order matters and spans three separate compiled tools before `mr_mot`
can run:

```sh
sh ex_compile            # -> mr_ex   (if not already built/run above)
sh e2p_compile            # -> e2p
sh mottable_compile        # -> motTab
sh mot_compile              # -> mr_mot
sh mot_comb_compile          # -> mot_comb   (optional, see below)
sh tau_mult_compile           # -> mot_mult
sh tauboc_compile              # -> tauBoc

./mr_ex 25 21          # writes 25_gens_data_ctau, 25_extables
./e2p 25               # writes 25_ex2poly_index
./motTab 25            # writes 25_mot_deltas, 25_poly_exponents
./mr_mot 25 21         # loads all of the above; writes 25_mot_gens, 25_mot_res, ...
./mot_mult 25 21       # loads 25_mot_gens/25_mot_res; writes h0.txt..h3.txt
./tauBoc 25 21         # loads 25_mot_gens/25_mot_res; writes tau_bockstein.txt etc.
```

- `mot_comb` is **not** a normal pipeline step — it re-runs only `mr_mot`'s
  final combine/compose stage from an interrupted or partial run, without
  re-resolving. Only needed if a `mr_mot` run was cut short after producing
  per-step files but before finishing. See `docs/pipelines/MOTIVIC.md` §2.
- Every tool here is p=2 in the actual arithmetic implemented (see
  `ARCHITECTURE.md` §5), whatever the degree arguments — there is currently
  no p=3 motivic path.

## Degree-parameter conventions, summarized

- All primary executables: `argv[1]` = max internal degree cutoff (the
  README calls it "half of `t`" as a project convention — see caveats
  above), `argv[2]` = resolution length `s` (number of resolution steps),
  except `BPtab` and `motTab`/`e2p`, which take only the degree argument
  (they don't run a resolution themselves).
- Always use the **same degree argument** across every tool in a pipeline
  run — this is how they find each other's output files (no other linkage
  mechanism exists).
- There is no runtime validation of matching degree arguments anywhere in
  the codebase (see the segfault mechanism above) — get this right by
  convention, not because the program will catch you if you don't.
