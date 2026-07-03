# The Motivic Pipeline (`mr_mot` / `mot_comb` / `mot_mult` / `tauBoc`)

Motivic layer of the codebase, built on the generic template framework
described in `docs/FRAMEWORK.md` (`matrix<ring>`, `ModuleOp`,
`comodule_generic`/`cofree_comodule`, `Hopf_Algebroid<ring,algebroid>`,
`curtis_table`, `polynomial<ring>`). Covers: `mot_steenrod.h/.cpp`,
`mot_main.cpp`, `mot_combine.cpp`, `mot_coactions.cpp`, `tao_bockstein.h/.cpp`,
`taubocmain.cpp`, `tao_boc.cpp`, `mot_mult.cpp`.

## 1. Purpose

- The ring being resolved is the **motivic dual Steenrod algebra**
  `A_*^mot = F_p[tau][xi_1, xi_2, ...]` (`motSteenrod = polynomial<tauPoly>`,
  `mot_steenrod.h:43`), a polynomial ring on `xi_i` with coefficients not in
  `F_p` but in `F_p[tau]` — `tau` a formal variable of its own bidegree
  `(0, weight)` (topological degree 0, motivic weight tracked separately,
  `MotDegree` class, `mot_steenrod.h:58-73`). Every class is bigraded by
  `(deg, weight) = (t, w)`.
- `tauOper` (`mot_steenrod.h:12-38`) implements arithmetic on (a sparse
  encoding of) monomials in `tau`: a `tauPoly` is just the exponent of `tau`
  packed into an `int16_t`, with a sentinel `internal_zero = 1<<15`
  (`mot_steenrod.h:14`) standing for the zero polynomial — i.e. this models
  `F_p[tau]` restricted to *monomials* `c·tau^k` (`c` implicitly absorbed:
  `add` of two nonzero same-valued terms returns `internal_zero`,
  `mot_steenrod.cpp:7-12`), not general polynomials. **TODO(math):** the
  comment "if x and y are both nontrivial, ... then x=y in F3[tau]"
  (`mot_steenrod.cpp:10`) suggests coefficients are only ever tracked
  mod 2 (so `x+x=0`), consistent with `unit(n)`'s parity check
  (`mot_steenrod.cpp:31-34`) and `lift()`'s `x%2!=0` test
  (`mot_steenrod.cpp:270-273`) — i.e. despite the `F3`/`F_p` naming
  throughout, the actual coefficient arithmetic implemented is mod **2**,
  not mod 3. See §6 and the parallel finding in `EX_CTAU.md` §1 — nothing in
  this file states why a repo billed as "p=3 fork" computes this pipeline's
  coefficients mod 2.
- `MotSteenrodOp` (`mot_steenrod.h:83-140`) is the `Hopf_Algebroid<tauPoly,
  motSteenrod>` structure on this ring: a cofree-comodule coaction table
  (`cofree_coaction`, loaded from disk, `generate_cofree_coaction`,
  `mot_steenrod.cpp:191-213`) built from a **separately precomputed**
  comultiplication table (see `mot_coactions.cpp`, §2). `mr_mot` then
  resolves the trivial comodule (§4) to compute (bigraded) `Ext` groups
  `Ext_{A_*^mot}^{s,(t,w)}(F_p, F_p)` — the intended **E2 page of the
  motivic/algebraic-Novikov-style Adams spectral sequence** for the motivic
  sphere. **TODO(math):** no file states the precise SS this feeds (motivic
  Adams vs. algebraic-Novikov-for-BP-motivic); it is inferred purely from
  the bigrading and the `tau`-Bockstein tooling built on top (see below).
- `tao_bockstein.{h,cpp}` / `taubocmain.cpp` (`tauBoc`) then treat the
  resulting minimal resolution as a **complex over `F_p[tau]`**
  (`motComplex`, `tao_bockstein.h:19-26`) and run a **tau-Bockstein spectral
  sequence**: `tau_table`/`tau_table_entry` (`tao_bockstein.h:28-74`)
  implement an algorithm that, for each resolution differential (viewed as
  a matrix over `F_p[tau]`), finds a "leading term" by `tau`-adic valuation
  (`tau_table::leading_term`, `tao_bockstein.cpp:72-78`) and iteratively
  cancels cycles against boundaries divisible by increasing powers of
  `tau` (`make_pretable`, `tao_bockstein.cpp:81-145`) — the standard
  Bockstein-SS bookkeeping (tag/cycle pairs with a "differential length" =
  power of `tau`) applied to the `tau`-adic filtration relating the motivic
  and classical (tau-inverted / tau=0) Adams `E2` pages.
  **TODO(math):** the code never states which classical target this
  Bockstein SS converges to (presumably the classical p-primary Adams `E2`
  page via `tau=1`, dual to setting weight aside) — it's inferred from the
  name and structure only.
- `mot_mult.cpp` (`mot_mult`) additionally computes **Massey-product-free
  multiplications** by fixed classes `h_i = MOP.hi(i)` (`mot_steenrod.h:124`,
  constructed as the monomial `x_{2^i}` with coefficient `tau^{-2^{i-1}}`,
  `mot_steenrod.cpp:230-231`) on the `E2`/associated-graded page, using the
  tau-Bockstein cycle data to express products in terms of chosen cycle
  representatives (`multiplication_table`, `tao_bockstein.cpp:304-336`).

## 2. Executables

All four share the `<maxdeg>_`-prefixed file-naming convention
(`directory = argv[1] + "_"`, e.g. `mot_main.cpp:10-12`).

### `mr_mot` (from `mot_main.cpp`)

- Build (`mot_compile:1`): `g++ -O2 -omr_mot Fp.cpp streams.cpp mot_main.cpp
  mot_steenrod.cpp exponents.cpp -fopenmp -Wall -Wfatal-errors -I./`
- Usage: `./mr_mot <maxdeg> <resolution_length>` (`mot_main.cpp:8-9`).
- **Inputs** (must already exist, all `<maxdeg>_`-prefixed):
  - `ex2poly_index` — read by `MOP.init_mon_array` (`mot_main.cpp:22`); this
    is the monomial list, produced by **`e2p`/`ex2poly.cpp`** (Part B, see
    `EX_CTAU.md`) — i.e. the motivic monomial ordering is inherited from the
    ex-pipeline's monomial index, not generated independently
    (`ex2poly.cpp:39-44` writes exactly `<maxdeg>_ex2poly_index`).
  - `mot_deltas`, `poly_exponents` — read by
    `MOP.generate_cofree_coaction` (`mot_main.cpp:25`); produced by
    **`motTab`/`mot_coactions.cpp`** (see §2 below; not one of the 4
    assigned executables but a required upstream step run first, via
    `mottable_compile`).
  - `gens_data_ctau` — read by `MOP.load_gens` (`mot_main.cpp:54`);
    produced by **`mr_ex`/`ex_main.cpp`** (`ex_main.cpp:17`,
    `SteenrodInit::save_gens`) — the motivic resolution's choice of
    generators at each stage is **modeled on** the classical ex/ctau
    resolution's generators (see `EX_CTAU.md` for that pipeline). This is
    the key cross-pipeline dependency: **`mr_ex` must run before `mr_mot`**
    for the same `<maxdeg>`.
  - `extables` — read inside `pre_resolution_modeled` (`mot_main.cpp:57`,
    consumed via `hopf_algebroid/6.h:41,48,78`); produced by `mr_ex`'s
    `SteenrodInit::resolve` → `pre_resolution_tab(..., director +
    "extables")` (`ex_steenrod_init.cpp:57`). This is a `curtis_table<Fp>`
    (over the *ex/ctau* coefficient ring `Fp`, not `tauPoly`) — the
    `tfm`/`transformer` lambda (`mot_main.cpp:44-45`, `MOP.lift`) converts
    `Fp` vectors from this table into `tauPoly` vectors on the fly
    (`mot_steenrod.cpp:276-281`).
- **Outputs** (all `<maxdeg>_`-prefixed): `mot_maps` (injection/quotient
  matrices per step), `mot_gens` (generator data per step, combined by
  `gens_file_combiner`, `mot_main.cpp:60`), `mot_res` (composed resolution,
  `mot_main.cpp:64`), intermediate `motinj`/`motqut`/`motindj`/`motnewm`,
  and `back<i>` backup files (`pre_resolution_modeled`'s
  `back_up_file_name`, `mot_main.cpp:57`). Also writes `maps.txt` (fixed
  name, cwd-relative, `mot_main.cpp:63`) and prints the monomial list to
  stdout (`MOP.output_monomials()`, `mot_main.cpp:23`).
- Control flow: construct `MotSteenrodOp` → load monomials → build cofree
  coaction from the `motTab` output → set comodule to trivial (`F_p[tau]`
  in degree `(0,0)`, `mot_main.cpp:32-33`) → load `ex`-modeled generators →
  `pre_resolution_modeled` (the actual resolution loop, modeled on the
  `mr_ex` curtis tables) → `gens_file_combiner` → `resolution()` (composes
  consecutive injection/quotient maps into the final differentials).

### `mot_comb` (from `mot_combine.cpp`)

- Build (`mot_comb_compile:1`): `g++ -O2 -omot_comb Fp.cpp streams.cpp
  mot_combine.cpp mot_steenrod.cpp exponents.cpp -fopenmp -Wall
  -Wfatal-errors -I./`
- Usage: `./mot_comb <maxdeg> <resolution_length>`
  (`mot_combine.cpp:8-9`) — same argument meaning as `mr_mot`.
- Purpose: **re-run only the combining/output stage** of `mr_mot` without
  re-resolving — it constructs `MotSteenrodOp` with `cofa=NULL`
  (`mot_combine.cpp:19`, no coaction matrix, so it cannot resolve, only
  reload), reloads the monomial list and the (already-computed) trivial
  comodule, then calls the same `gens_file_combiner` +
  `resolution()` pair as the tail of `mr_mot` (`mot_combine.cpp:35,39`
  vs. `mot_main.cpp:60,64`).
- Inputs: `ex2poly_index`, `mot_coac_d` (a `matrix_file<motSteenrod>` for the
  comodule, `mot_combine.cpp:25`), and the **per-step** `mot_gens<i>` /
  `mot_maps<i>` files left behind by a `mr_mot` run — i.e. this is meant to
  be run **after a `mr_mot` run that produced per-step but not yet combined
  output** (e.g. a run interrupted before its final combine/resolve step,
  or one where those two steps are deliberately split out).
- Outputs: `mot_gens` (combined), `mot_res`, `maps.txt` — identical output
  set/names to `mr_mot`'s tail.
- **TODO(math)/TODO(code):** no comment states why this split executable
  exists (e.g. to avoid recomputing the expensive resolution loop when only
  the combine/compose step needs re-running, or to combine partial output
  from a crashed/resumed `mr_mot` run) — inferred purely from the code
  overlap with `mot_main.cpp`'s tail.

### `motTab` (from `mot_coactions.cpp`) — upstream helper, not one of the 4 named executables but load-bearing

- Build (`mottable_compile:1`): `g++ -std=c++11 -O2 exponents.cpp
  mot_steenrod.cpp mot_coactions.cpp mon_index.cpp -omotTab -I./ -fopenmp
  -Wall -Wfatal-errors`. Note it links `mon_index.cpp`'s **classical**
  `monomial_index` (built on `exponents.h`, no exterior part) — not
  `ex_index.cpp`'s variant — since `mot_coactions.cpp:3` includes
  `mon_index.h` directly.
- Usage: `./motTab <maxdeg>`; writes `<maxdeg>_mot_deltas` (comultiplication
  table) and `<maxdeg>_poly_exponents` (monomial list in `motTab`'s own
  ordering) (`mot_coactions.cpp:56-71`).
- Computes the **motivic Milnor coproduct**: for generator `x_n`,
  `Δ(x_n) = x_n⊗1 + 1⊗x_n + Σ_{i=1}^{n-1} τ^{-2^{i-1}}[x_i|x_{n-i}^{2^i}]`
  (`mot_coactions.cpp:20-27`), with the comment "`xi_i = τ^{-1}τ_i^2`" at
  line 23 — i.e. the motivic `xi_i` here is being expressed via a "square
  root" generator `τ_i` exactly as in the ex/ctau pipeline (see
  `EX_CTAU.md` §1), extended to a genuine coproduct formula over
  `F_p[tau]` using negative powers of `tau`. **TODO(math):** the formula's
  relationship to the classical/Milnor coproduct at `p=2` vs. the intended
  motivic coproduct of Voevodsky/Hoyois-Kelly-Østvær is not stated in the
  comments — it's presented as a direct, unexplained generalization.
- This produces the `mot_deltas`/`poly_exponents` files consumed by
  `mr_mot`'s `generate_cofree_coaction` (`mot_main.cpp:25`) — **must run
  before `mr_mot`**, alongside `mr_ex` (for `gens_data_ctau`/`extables`) and
  `e2p` (for `ex2poly_index`). So the full pipeline order is:
  `mr_ex` (and/or `e2p`) → `motTab` → `mr_mot` → (`mot_comb` if needed) →
  `mot_mult` / `tauBoc`.

### `mot_mult` (from `mot_mult.cpp`)

- Build (`tau_mult_compile:1`): `g++ -std=c++11 -O2 mot_mult.cpp
  tao_bockstein.cpp mot_steenrod.cpp exponents.cpp -fopenmp -I./ -omot_mult`
- Usage: `./mot_mult <maxdeg> <resolution_length>` (`mot_mult.cpp:9-10`).
- Inputs: `ex2poly_index`, `mot_gens`, `mot_res` — i.e. **runs after
  `mr_mot`/`mot_comb`** have produced the combined resolution.
- Control flow (`mot_mult.cpp:24-47`): load the resolution as a
  `motComplex` (`Complex.load`) → `make_table` (build the tau-Bockstein
  tables, same as `tauBoc`, see below) → `make_cycle_tables` (extract
  cycle representatives per resolution degree) → for `i=0..3`, call
  `multiplication_table(MOP.hi(i), 2^i, ...)` (`mot_mult.cpp:44-47`) to
  compute multiplication-by-`h_i` maps on cycles, writing
  `h0.txt`..`h3.txt`.
- Outputs: `tau_bockstein.txt` (same as `tauBoc`, recomputed redundantly —
  `mot_mult` does not reuse a `tauBoc` run's output, it recomputes the
  tables itself), `new_generators` (cycle data dump,
  `mot_mult.cpp:36-41`), `h0.txt`..`h3.txt` (multiplication-by-`h_i`
  tables).

### `tauBoc` (from `taubocmain.cpp`)

- Build (`tauboc_compile:1`): `g++ -std=c++11 -O2 taubocmain.cpp
  tao_bockstein.cpp mot_steenrod.cpp exponents.cpp -fopenmp -I./ -otauBoc`
- Usage: `./tauBoc <maxdeg> <resolution_length>`
  (`taubocmain.cpp:7-8`) — same inputs as `mot_mult` (`ex2poly_index`,
  `mot_gens`, `mot_res`), also produced by `mr_mot`/`mot_comb`.
- Outputs: `mot_gen.txt` (human-readable generator degrees,
  `output_generators`, `taubocmain.cpp:23`), `mot_maps.txt` (human-readable
  differentials, `Complex.output`, `taubocmain.cpp:32`),
  `tau_bockstein.txt` (the Bockstein table, `taubocmain.cpp:38-39`).
- This is the **read-only / reporting** counterpart of `mot_mult`: it
  builds and dumps the tau-Bockstein tables and human-readable resolution
  data but does not compute any `h_i`-multiplications.

## 3. Class catalog

| Class | File | Template Params | Base Classes | Purpose |
|---|---|---|---|---|
| `tauPoly` | `mot_steenrod.h:10` | `typedef int16_t` | — | Encodes a single monomial `c·tau^k` in `F_p[tau]` as its `tau`-exponent; sentinel `internal_zero=1<<15` for `0` (`mot_steenrod.h:14`) |
| `tauOper` | `mot_steenrod.h:12` | — | `virtual RingOp<tauPoly>` | Arithmetic on (monomials of) `F_p[tau]`; `add` of distinct-valuation or equal nonzero terms collapses per mod-2 rule (`mot_steenrod.cpp:7-12`, see §1 TODO); `tau_valuation`, `shift`, `power_tau` for Bockstein bookkeeping |
| `motSteenrod` | `mot_steenrod.h:43` | `typedef polynomial<tauPoly>` | — | The motivic dual Steenrod algebra `F_p[tau][xi_1,xi_2,...]` |
| `mStmSt` | `mot_steenrod.h:45` | `typedef polynomial<motSteenrod>` | — | `A_*^mot ⊗ A_*^mot`, used for the coproduct |
| `MotSteenrodRingOp` | `mot_steenrod.h:48` | — | `virtual PolynomialOp<tauPoly>` | Ring/module ops on `motSteenrod`; wires `output_term` (`mot_steenrod.cpp:77-80`) |
| `MotDegree` | `mot_steenrod.h:58` | — | — | Bigrading `(deg, weight)`; `varWeight(n) = xnDeg(n)/2` (`mot_steenrod.cpp:113`) — motivic weight of `xi_n` is half its topological degree |
| `FreeMotCoMod` | `mot_steenrod.h:76` | `typedef cofree_comodule<motSteenrod,MotDegree>` | — | Cofree `A_*^mot`-comodule bigraded by `MotDegree` |
| `MotSteenrodOp` | `mot_steenrod.h:83` | — | `virtual Hopf_Algebroid<tauPoly,motSteenrod>` | The Hopf algebroid structure on `A_*^mot`; owns `cofree_coaction`, `mon_array`/`mon_index`, `etaL`/`etaR`, `delta`, `hi(i)`, `lift()` (F_p→F_p[tau] and vector lift), `algebroid2vector[St]`/`vector2algebroid` |
| `motComplex` | `tao_bockstein.h:19` | — | `complex<MotDegree,tauPoly>` | Loads a resolution's differentials as matrices over `tauPoly` from `mot_gens`/`mot_res` (`tao_bockstein.cpp:24-61`) |
| `cycle_name` | `tao_bockstein.h:28` | — | — | Single private `int i`; unused elsewhere in the read files — likely a stub/leftover |
| `tau_table_entry` | `tao_bockstein.h:33` | — | — | One row of the tau-Bockstein table: `tag`/`cycle` indices, full vectors, `diff_length` (power of `tau` of the differential) |
| `cycle_data` | `tao_bockstein.h:49` | `typedef std::map<int,vectors<matrix_index,tauPoly>>` | — | Map from generator index to its chosen cycle representative |
| `tau_table` | `tao_bockstein.h:52` | — | — | The Bockstein table at one resolution degree: `make_pretable` (core recursive cancellation algorithm, `tao_bockstein.cpp:81-145`), `get_cycles`/`get_tags`, `pot_maker` |

## 4. Data/control flow per executable

- **`mr_mot`**: `MotSteenrodOp(coa, maxdeg)` → `init_mon_array` (from
  `ex2poly_index`) → `generate_cofree_coaction` (from `motTab`'s
  `mot_deltas`/`poly_exponents`) → `set_to_trivial` (comodule = `F_p[tau]`
  at `(0,0)`) → `load_gens` (from `mr_ex`'s `gens_data_ctau`) →
  `pre_resolution_modeled` (loop over resolution steps, each guided by the
  `mr_ex` curtis table `extables` via the `Fp→tauPoly` lift transformer) →
  `gens_file_combiner` → `resolution()` (compose maps, write `mot_res`).
- **`mot_comb`**: identical tail to `mr_mot` (`gens_file_combiner` +
  `resolution()`) without the resolving step — reloads existing per-step
  `mot_gens<i>`/`mot_maps<i>` files.
- **`motTab`**: `monomial_index(maxdeg)` (classical, non-exterior) →
  `init_mon_array` → `make_delta_table` (recursive motivic Milnor coproduct
  via `mon_index.substitution_table`) → write `mot_deltas` +
  `poly_exponents`.
- **`mot_mult`**: load `motComplex` from `mot_gens`/`mot_res` →
  `make_table` (tau-Bockstein) → `make_cycle_tables` → for each `h_i`
  (`i=0..3`), `multiplication_table` (builds a multiplication-by-`h_i`
  matrix via `make_multiplication_table`, pushes each generator's cycle
  through it, re-expresses the result in the next degree's chosen cycles
  via `find_cycle`) → write `h<i>.txt`.
- **`tauBoc`**: load `motComplex` → `output_generators` (`mot_gen.txt`) →
  `Complex.output` (`mot_maps.txt`) → `make_table` → `output_tables`
  (`tau_bockstein.txt`).

## 5. Open math questions

- **TODO(math):** Why do `tauOper`'s arithmetic (`mot_steenrod.cpp:7-12`),
  `unit(n)`'s parity test, and `MotSteenrodOp::lift`'s `x%2!=0` test all
  operate mod 2, when the surrounding type is named `F3` (`mot_steenrod.h:7,
  typedef Fp F3`) and the repo is framed as a "p=3 fork"? Is the entire
  motivic pipeline still computing at `p=2` (unported), or is there a
  mod-3-compatible reading of "add of equal terms → 0" that isn't literally
  mod-2 arithmetic?
- **TODO(math):** What is the precise motivic coproduct formula
  `mot_coactions.cpp:20-27` computing, relative to a standard reference for
  the motivic dual Steenrod algebra (e.g. Voevodsky) — in particular the
  role of the negative powers of `tau` (`τ^{-2^{i-1}}`) and whether the
  resulting `Ext` groups are literally `Ext_{A_*^mot}(F_p,F_p)` (motivic
  Adams `E2`) or something else bigraded compatibly.
- **TODO(math):** What SS the tau-Bockstein tooling (`tao_bockstein.*`,
  `tauBoc`) is computing the differentials of — is it the motivic
  analogue of the classical mod-`p`-Bockstein SS (converging the mod-`p^k`
  Adams `E2` pages), or the "tau-Bockstein SS" relating motivic and
  classical Adams `E2` pages via inverting `tau`? No file states the target.
- **TODO(math)/TODO(code):** Why does `mot_comb` exist as a separate
  executable from `mr_mot`'s tail rather than a flag/mode of `mr_mot` —
  purely a resumability/recomputation-avoidance tool, or something else?
- **TODO(math):** The exact justification for "modeling" the motivic
  resolution's generator choices and cell structure on the ex/ctau
  pipeline's `F_2` resolution (`pre_resolution_modeled`, `extables`,
  `gens_data_ctau`) — presumably because `A_*^mot` reduces to (a
  square-root cover of) the classical mod-2 dual Steenrod algebra when
  `tau` is inverted/forgotten, so the same minimal generating set works,
  but no comment states this identification explicitly.
