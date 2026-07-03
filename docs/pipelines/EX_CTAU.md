# The Ex/Ctau Pipeline (`mr_ex`)

An "exterior ⊗ polynomial" reformulation of a mod-2 dual-Steenrod-algebra
computation, built on the same generic template framework as the classical
and motivic pipelines (`docs/FRAMEWORK.md`, `docs/pipelines/STEENROD.md`,
`docs/pipelines/MOTIVIC.md`). Covers: `ex_exponents.h/.cpp`, `expoly.h`,
`ex_index.h/.cpp`, `ex_steenrod_init.h/.cpp`, `ctau_steenrod.h/.cpp`,
`ex_main.cpp`, `ex_test.cpp`, `ex2poly.cpp`. This document's most actionable
content is §3 (duplication) and §5 (dead code) — read those first.

## 1. Purpose (best-effort inference — speculative)

- `ctau_steenrod.h` defines the dual object `P = poly<ex_poly,Fp>`
  (`ctau_steenrod.h:10`): a ring with **two kinds of generators**,
  polynomial `x_i` (`xnDegs = {0,2,6,14,30,62,126,254}`, `ex_exponents.cpp:46`)
  and **exterior** `z_i` (`znDegs = {0,1,3,7,15,31,63,127}`,
  `ex_exponents.cpp:51`, with `z_i^2=0` enforced by `common_ex`,
  `ex_exponents.cpp:15-19`, and `ExPolyOp_Para::mon_multiply`,
  `ex_exponents.h:158-172`). Note **`znDegs[i] == xnDegs[i]/2` for every
  `i`** — i.e. `x_i` and `z_i` are set up with exactly the degree relation
  `|x_i| = 2|z_i|`, the same relation `mot_coactions.cpp:23`'s comment
  states for the motivic algebra (`"ξ_i = τ^{-1}τ_i^2"`, see
  `MOTIVIC.md` §2). **Inference:** `P` here is very likely the classical,
  non-motivic, **mod-2** dual Steenrod algebra `A_* = F_2[xi_1,xi_2,...]`
  represented via an auxiliary "square-root" presentation `x_i = z_i^2`
  (with `z_i` an odd/exterior formal square root, playing the role `tau_i`
  plays in the motivic algebra with `tau` set to a constant/forgotten) —
  hence the name "**c-tau**" ≈ "constant/classical tau": the same
  square-root trick used motivically, but with the `tau`-grading collapsed
  (no separate weight, no formal `tau` variable) so ordinary `Hopf_Algebroid<
  Fp,P>` machinery applies directly. This would also explain why this
  variant needs its own `monomial_index` with an exterior half (`ex_index.h`,
  §3) that the plain classical pipeline's `mon_index.h` doesn't have.
  **TODO(math):** this reading is not stated anywhere in the code/comments;
  confirm against `ctau_steenrod.cpp:52-77`'s coproduct formulas (`deltaXn`,
  `deltaZn`) matching the classical odd-primary Milnor coproduct on
  `ξ_j`/`τ_j` specialized/doubled to `p=2`.
  - Corroborating: `ctau_steenrod.cpp:54-59`'s `deltaXn(j)` = `x_j⊗1 + 1⊗x_j
    + Σ_{k=1}^{j-1} x_k⊗x_{j-k}^{power_p(k)}` and `deltaZn(j)`
    (`ctau_steenrod.cpp:62-74`) = `z_j⊗1 + 1⊗z_j + Σ_{k=1}^{j-1}
    z_k⊗x_{j-k}^{power_p(k-1)}` are exactly the shape of the classical
    Milnor coproduct on `ξ_j` and `τ_j` (odd-primary dual Steenrod
    algebra), with `power_p(i)=prime^i` — and `ex_main.cpp:11` hardcodes
    `prime=2` regardless of `argv`, so `power_p(i)=2^i`.
  - This directly answers (with a plausible mechanism, not proof) the open
    question left in `docs/pipelines/STEENROD.md` §6 about why the
    classical pipeline's `P = polynomial<Fp>` (`steenrod.h:10`) has *no*
    exterior `tau_i` part: **the exterior part is handled here, in this
    separate ex/ctau pipeline**, not in `steenrod.h`.
- Relation to the motivic pipeline (`MOTIVIC.md`): `mr_mot` explicitly
  consumes this pipeline's output as a *model* for its own resolution
  (`gens_data_ctau`, `extables`, both written by `ex_main.cpp`/
  `ex_steenrod_init.cpp` and read by `mot_main.cpp:54,57` — see
  `MOTIVIC.md` §2). So concretely, independent of the speculative math
  above: **`mr_ex`'s resolution structure (its curtis tables and chosen
  generators) is reused verbatim to drive the motivic resolution.** This is
  the strongest, code-verified relationship between the two pipelines.
- Relation to the classical pipeline (`steenrod.h`/`mr_st`, `p=3`): same
  class names (`P_Op`, `PP_Op`, `Steenrod_Op`, `ComodInit`, `SteenrodInit`)
  and same overall driver shape, but a structurally different ring (`P =
  poly<ex_poly,Fp>` here vs. `polynomial<Fp>` there) and — per `ex_main.cpp:
  11` — a hardcoded **`p=2`**, whereas `stmain.cpp:11` hardcodes `p=3`
  (see `STEENROD.md` §2). **TODO(math):** is `mr_ex` intentionally a
  `p=2`-only auxiliary tool (independent of the `p=3` fork's main target),
  or a leftover from the pre-fork p=2 codebase not (yet) ported to p=3?

## 2. Executable

### `mr_ex` (from `ex_main.cpp`, header comment still says `stmain.cpp`)

- Build (`ex_compile:1`): `g++ -omr_ex ex_steenrod_init.cpp
  ex_exponents.cpp ex_index.cpp ex_main.cpp ctau_steenrod.cpp Fp.cpp
  streams.cpp -I./ -std=c++14 -fopenmp -O2 -Wfatal-errors`
- Usage: `./mr_ex <maxdeg> <resolution_length>` (`ex_main.cpp:8-9`).
  `ex_main.cpp:1` still carries the comment `//stmain.cpp` — this file is a
  copy of the classical pipeline's `stmain.cpp` (diffed: only the include,
  the hardcoded prime `3→2`, the constructor's extra `director` argument,
  and output filenames differ; see §3).
- Prime is **hardcoded to 2**: `SteenrodInit st(2, maxdeg, length, ...)`
  (`ex_main.cpp:11`), regardless of any "p=3 fork" framing elsewhere in the
  repo.
- Inputs: none required upfront — `<maxdeg>_ctau_steenrod_coaction.data` is
  self-generated on first run (constructor default `IorO=false`,
  `ex_steenrod_init.cpp:6-20`), same pattern as the classical `SteenrodInit`.
  Note: this pipeline does **not** consume `ex2poly_index`/`e2p`'s output —
  `ex2poly.cpp`/`e2p` instead **depends on this pipeline's monomial
  ordering** (`ex2poly.cpp:41-44` builds its own `monomial_index` from
  `ex_index.h` and just re-derives the doubled exponent list; it does not
  read any `mr_ex` output file). The dependency direction for the *motivic*
  chain is: `mr_ex` (or a bare `monomial_index`) and `e2p` both derive
  independently from the same `ex_index.h` monomial enumeration, and both
  feed `mr_mot` (see `MOTIVIC.md` §2).
- Outputs (all `<maxdeg>_`-prefixed):
  - `ctau_steenrod_coaction.data` — the coproduct table (self-generated).
  - `maps`, `gens` — per `SteenrodInit::resolve` (`ex_steenrod_init.cpp:56-58`,
    calls `Steenrod_Op::pre_resolution_tab`).
  - `extables` — the `tablename` argument to `pre_resolution_tab`
    (`ex_steenrod_init.cpp:57`) — **consumed by `mr_mot`**
    (`mot_main.cpp:57`, `MOTIVIC.md` §2).
  - `ResTables_ctau` — `SteenrodInit::saveResolutionTables`
    (`ex_main.cpp:16`).
  - `gens_data_ctau` — `SteenrodInit::save_gens` (`ex_main.cpp:17`) —
    **consumed by `mr_mot`** (`mot_main.cpp:54`).
  - stdout: `Steenrod_Op::output_resolution` rank table
    (`ex_main.cpp:19`, shared implementation with the classical pipeline's
    `steenrod.cpp`, but here in `ctau_steenrod.cpp:156-191`).
  - Intermediate matrix files: `inex`, `indjex`, `qutex`, `mewexp`
    (injection/quotient/new-map scratch, `ex_steenrod_init.cpp:5`) and
    `excac`/`excoamat` (the comodule's backing coaction matrix,
    `ex_steenrod_init.cpp:5,52`) — note these are **`matrix_file<...>`**
    (disk-backed), not `matrix_mem<...>` as in the classical pipeline (§3).

## 3. Duplication analysis

This pipeline is built almost entirely from **near-duplicates of files used
elsewhere** in the repo. Every pair below was diffed directly; findings are
concrete, not inferred.

| Pair | Locations | Verdict |
|---|---|---|
| `ex_poly`/`ExPolyOp_Para` in **`ex_exponents.h`** vs. **`expoly.h`** | `ex_exponents.h:1-177` vs. `expoly.h:1-177` | **Byte-for-byte identical** (`diff` returns no output). `expoly.h:1` even still carries the header comment `//ex_exponents.h`. Nothing in the repo `#include`s `expoly.h` (checked all `.h`/`.cpp` includes) — it is a **fully orphaned, unreferenced exact copy** of `ex_exponents.h`, not compiled or used anywhere. Whoever/whatever created it appears to have duplicated the file (perhaps as a WIP fork point) and never diverged or wired it in. |
| `monomial_index` in **`ex_index.h`** vs. **`mon_index.h`** | `ex_index.h:10-99` vs. `mon_index.h:10-95` | **Genuinely divergent, but clearly copy-paste-derived.** `ex_index.h:1` still carries the header comment `//mon_index.h` (a leftover from copying `mon_index.h` as the starting point). Real differences: `ex_index.h` includes `ex_exponents.h` (maxVar=7, has exterior generators) vs. `mon_index.h` includes `exponents.h` (maxVar=5, no exterior generators, `exponents.h:3,7`); `ex_index.h`'s `substitution_table` takes an extra `exvalues` callback and has a `max_ex` member (`ex_index.h:42-43,80-89,36`) to handle the exterior half of `ex_poly` that `mon_index.h`'s version has no concept of. So the divergence is functionally necessary (mon_index.h's monomial enumeration has nothing to enumerate for exterior generators), but the file was clearly forked from `mon_index.h` wholesale rather than written fresh or factored via inheritance/templates — the non-exterior 80% of the two files is line-for-line parallel. |
| `ComodInit`/`SteenrodInit` in **`ex_steenrod_init.h`** vs. **`steenrod_init.h`** | `ex_steenrod_init.h:1-56` vs. `steenrod_init.h:1-50` | **Near-duplicate class definitions, confirmed.** Same two class names, same member list shape. Differences: `ex_steenrod_init.h` includes `ctau_steenrod.h` (not `steenrod.h`) and `matrices_stream.h` (not present in `steenrod_init.h`); its `ComodInit` backs onto `matrix_file<P>` (disk-backed, `ex_steenrod_init.h:11`) vs. `steenrod_init.h`'s `matrix_mem<P>` (in-memory, `steenrod_init.h:8`); likewise `inj`/`indj`/`qut`/`new_map` are `matrix_file<Fp>` here (`ex_steenrod_init.h:30-31`) vs. `matrix_mem<Fp>` there (`steenrod_init.h:27`); the `SteenrodInit` constructor takes an extra `string director` parameter here (`ex_steenrod_init.h:46` vs. `steenrod_init.h:40`). This looks like a deliberate storage-strategy variant (disk-backed matrices, presumably for a larger `p=2` computation) layered onto an otherwise copy-paste-identical class structure — **genuinely divergent in intent (memory vs. disk), but structurally copy-drift** (the two files could plausibly have been unified via a template parameter on storage backend). |
| `ex_main.cpp` vs. classical `stmain.cpp` | `ex_main.cpp:1-22` vs. `stmain.cpp:1-22` | **Direct copy**, diff is 3 hunks: include (`ex_steenrod_init.h` vs `steenrod_init.h`), hardcoded prime (`2` vs `3`) and an extra `filename` ctor arg, and output filenames (`..._ctau` suffix vs. none). `ex_main.cpp:1` still says `//stmain.cpp`. |

**Overall assessment:** this pipeline was very likely bootstrapped by
copying four files from the classical pipeline (`mon_index.h`→`ex_index.h`,
`steenrod_init.h`→`ex_steenrod_init.h`, `steenrod.h`→`ctau_steenrod.h`,
`stmain.cpp`→`ex_main.cpp`) and adapting each for the exterior-generator /
disk-backed-matrix / `p=2` variant, plus one stray exact duplicate
(`expoly.h`) that appears to be an abandoned copy of `ex_exponents.h` with
no callers at all. None of this is inherently wrong, but it means any bugfix
or improvement to the classical `steenrod_init.h`/`mon_index.h`/`steenrod.h`
machinery will silently **not** propagate to this pipeline's copies, and
vice versa.

## 4. Class catalog

| Class | File | Template Params | Base Classes | Purpose |
|---|---|---|---|---|
| `ex_poly` | `ex_exponents.h:16` | — | — | Packed exponent for `P = F[x_i]⊗Λ(z_i)`: `exponent e` packs polynomial powers in low bits + exterior bitmask in high bits (`ex_exponents.cpp:5-19`) |
| `ExPolyOp_Para<ring>` | `ex_exponents.h:155` | `ring` | `PolyOp_Para<ex_poly,ring>` | Multiplication respecting `z_i^2=0` (`mon_multiply`, `ex_exponents.h:157-172`, skips via `common_ex`) |
| `monomial_index` (ex variant) | `ex_index.h:10` | — | — | Enumerates monomials of `P` (poly ⊗ exterior) up to `max_degree`; `substitution_table` takes both `values` (poly gens) and `exvalues` (exterior gens) callbacks (`ex_index.h:42-43`) |
| `P` | `ctau_steenrod.h:10` | `typedef poly<ex_poly,Fp>` | — | The "ctau" dual object: `Fp[x_1..x_7] ⊗ Λ(z_1..z_7)` |
| `PP` | `ctau_steenrod.h:12` | `typedef poly<ex_poly,P>` | — | `P ⊗ P`, for the coproduct |
| `P_Op` | `ctau_steenrod.h:18` | — | `ExPolyOp_Para<Fp>` | Ring ops on `P`; `construct(k,i)`=`x_k^i`, `construct_ex(k)`=`z_k` (`ctau_steenrod.cpp:126-133`) |
| `PP_Op` | `ctau_steenrod.h:29` | — | `virtual ExPolyOp_Para<P>` | Ops on `P⊗P`; `construct(k1,i1,k2,i2)`, `construct_expoly(k1,k2,i2)`=`z_{k1}⊗x_{k2}^{i2}` (`ctau_steenrod.cpp:136-143`) |
| `Steenrod_Op` (ctau) | `ctau_steenrod.h:42` | — | `virtual Hopf_Algebroid<Fp,P>`, `virtual Fp_Op` | Hopf algebroid on `P`; `make_delta` builds `deltaXn`/`deltaZn` coproduct formulas (`ctau_steenrod.cpp:52-77`, see §1) |
| `ComodInit` (ex) | `ex_steenrod_init.h:10` | — | `SteenrodCoMod_generic` | Backs onto **`matrix_file<P>`** (disk), ctor takes a directory string (`ex_steenrod_init.h:11-15`) |
| `SteenrodInit` (ex) | `ex_steenrod_init.h:18` | — | (none) | Driver; `inj`/`indj`/`qut`/`new_map` are `matrix_file<Fp>` (disk); commented-out `curtisTable_stream`/`unique_ptr` fields (`ex_steenrod_init.h:36-37`) suggest an abandoned streaming-curtis-table variant |
| `FreeSteenrodCoMod`/`SteenrodCoMod_generic` (ctau) | `ctau_steenrod.h:103,106` | `typedef cofree_comodule<P,int>` / `comodule_generic<P,int>` | — | Same shape as classical pipeline's identically-named typedefs (`steenrod.h:99,102`) |

## 5. Dead code candidates

Cross-checked against **every** compile script in the repo (`*_compile`,
10 scripts total), not just the 5 named in this task:

| File | Referenced by any compile script? | Verdict |
|---|---|---|
| `expoly.h` | No — and not even `#include`d by any `.cpp`/`.h` | **Dead**: an exact orphaned duplicate of `ex_exponents.h` (§3). Safe to delete once confirmed unused by any out-of-tree/WIP work. |
| `ex_test.cpp` | No | **Dead/standalone dev utility**: `#include`s only `exponents.h` (the *classical* one, maxVar=5, no exterior part — not `ex_exponents.h`!) and has a bare `main(){ ... }` (no return type, pre-C++ style, `ex_test.cpp:3`) that builds a `monomial_index` over degree 142 and dumps degrees/monomials to stdout. Reads like a scratch test of the classical `monomial_index`, not actually exercising anything `ex`-specific despite its name and location among the `ex_*` files. |
| `tao_boc.cpp` | No (only `tao_bockstein.cpp` is referenced, by `tau_mult_compile`/`tauboc_compile`) | **Dead — orphaned older revision.** Confirmed by direct inspection: `tao_boc.cpp` uses namespaces/types that do not match the current headers at all — `TauBockstein::`, `Matrices::matrix_array`, `Modules::vectors`, `HopfAlgebroid::cofree_comodule`, `MotSteentrod::tauPoly`, `motSt`, `tauP`, `new_vector`, and a `cycle_data` used as a *class* with `.insert()`/`.search()` methods (`tao_boc.cpp:240-260`) — whereas the current `tao_bockstein.h:49` defines `cycle_data` as a plain `std::map` typedef with no methods. This is a **pre-refactor snapshot** of what is now `tao_bockstein.cpp`, kept in the tree but never updated or wired into any build. |
| `ex2poly.cpp` | **Yes** — built into `e2p` via `e2p_compile:1` (`g++ ... ex_exponents.cpp ex_index.cpp ex2poly.cpp ... -oe2p`) | **Not dead** — it's a legitimate small standalone tool, just not one of the 5 executables named in this task's scope. Included here only because the task asked to check; it *is* wired into the motivic chain as the producer of `<maxdeg>_ex2poly_index` (`MOTIVIC.md` §2). |
| `mot_coactions.cpp` (cross-reference) | **Yes** — built into `motTab` via `mottable_compile` (out of this doc's Part-B scope, documented in `MOTIVIC.md` §2) | Not dead; flagged here only to close the loop since the task listed it as "check if it has its own main." |

**Recommendation for the repo owner:** `expoly.h`, `ex_test.cpp`, and
`tao_boc.cpp` can very likely be deleted outright — none are referenced by
any of the 10 compile scripts in the repo, and `tao_boc.cpp` in particular
is clearly superseded by `tao_bockstein.cpp` (same filename family, totally
different/older API).

## 6. Open math questions

- **TODO(math):** Confirm or refute the central speculative claim of §1 —
  that `ctau_steenrod.h`'s `P = poly<ex_poly,Fp>` is the ordinary mod-2 dual
  Steenrod algebra presented via square-root generators `x_i=z_i^2`,
  by checking `ctau_steenrod.cpp:52-77`'s coproduct formulas against a
  known reference for the Milnor coproduct at `p=2`.
- **TODO(math):** Why is the prime hardcoded to `2` in `ex_main.cpp:11`
  in a repository whose stated purpose is a `p=3` fork? Is `mr_ex` an
  intentionally `p=2`-specific auxiliary computation feeding the motivic
  pipeline (which itself appears to also be stuck at `p=2`, see
  `MOTIVIC.md` §1/§6), or unported legacy code?
- **TODO(math):** What, precisely, does "modeling" the motivic resolution
  on this pipeline's curtis tables/generator choices (`gens_data_ctau`,
  `extables`, consumed by `mr_mot`) mathematically justify — i.e. why
  should the minimal generating set and cell structure of `Ext_P(F_2,F_2)`
  for this `ex_poly` presentation coincide with (or lift to) the minimal
  generating set for `Ext_{A_*^mot}(F_p,F_p)`? No comment states this.
- **TODO(math):** The commented-out `curtisTable_stream`/streaming variant
  in `ex_steenrod_init.h:36-37`/`ex_steenrod_init.cpp:40-41` suggests an
  abandoned attempt at a disk-streamed (rather than fully in-memory per
  table) curtis-table implementation for this pipeline specifically —
  unclear whether this was abandoned for correctness or performance
  reasons, or simply unfinished.
