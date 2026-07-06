# A Second Hopf Algebroid: Resolving F_p over Gamma = F_p[x]/(x^n - 1)

This page documents a second, independent test of the `Hopf_Algebroid<ring,
algebroid>` framework — separate from the `BP_*BP`/dual-Steenrod-algebra
pipelines the rest of `docs/` covers. It's the smallest nontrivial instance
of the general "construct a minimal cofree resolution of a comodule"
algorithm described in `MinimalResolution.pdf` §5 (the same algorithm
`docs/CODE_WALKTHROUGH.md` §4.0 traces through the code for the sphere), and
it exists mainly to exercise `Hopf_Algebroid`'s resolution machinery
end-to-end on a Hopf algebroid that has nothing to do with `BP`.

## The Hopf algebroid

- **Base ring**: `F_p` (a field), for a chosen prime `p`.
- **Algebroid**: `Gamma = F_p[x]/(x^GammaDim - 1)` — the truncated polynomial
  ring on one generator `x` of degree 1, truncated above degree
  `GammaDim-1` (`trunc_hopf.h` fixes `GammaDim=5`, i.e. `Gamma` has basis
  `1, x, x^2, x^3, x^4`; edit that constant to try other heights).
- **Left/right units**: both equal to the standard inclusion `F_p ->
  Gamma`, `c |-> c*1` (`TruncHopf_Op::etaL`/`etaR`, `trunc_hopf.cpp`).
- **Comultiplication**: determined by `Delta(x) = 1 (x) x + x (x) 1` (`x` is
  primitive) extended multiplicatively, i.e.
  `Delta(x^n) = sum_{i=0}^n binom(n,i) x^i (x) x^{n-i}`
  (`TruncHopf_Op::delta`, `trunc_hopf.cpp` — computed directly from the
  binomial-coefficient formula, not from a coproduct table).

Because the base ring `F_p` is already a field, `pre_resolution_tab` (the
from-scratch resolver — see `docs/GENERAL_COMODULES.md`'s explanation of why
that function needs a field) applies **directly**, exactly as it does for
the dual Steenrod algebra in `mr_st`. There is no `BP_*`-style non-field
base ring here, so none of the mod-`I`-reduce-then-lift machinery in
`BP_generic_init.*`/`Steenrod_generic_init.*`/`BP_mod_I.*` is needed — this
example is architecturally simpler than the general-`BP_*BP`-comodule case.

## Files

- **`trunc_hopf.h`/`.cpp`** — `TruncHopf_Op : public Hopf_Algebroid<Fp,Gamma>,
  public Fp_Op`, mirroring `Steenrod_Op`'s exact wiring
  (`Gamma_opers`/`FpMod_opers`/`GammaMod_opers` composed the same way
  `Steenrod_Op` composes `P_opers`/`FpMod_opers`/`PMod_opers`). Since
  `Gamma`'s basis element `x^n` is already indexed by its own exponent `n`
  (a single variable, no `exponents.h`-style packing needed the way the
  multi-variable dual Steenrod algebra requires), `algebroid2vector`/
  `vector2algebroid` are just the identity map on indices, shifted.
- **`trunc_init.h`/`.cpp`** — `TruncInit`, mirroring `SteenrodInit` exactly:
  builds the trivial comodule `F_p` via the inherited `set_to_trivial`, then
  calls the inherited `pre_resolution_tab`.
- **`mr_trunc.cpp`** — the executable, mirroring `stmain.cpp`.
- **`trunc_compile`** — build script, mirroring `st_compiling`.

## How to use it

```
sh trunc_compile
./mr_trunc <max_degree> <resolution_length> <p>
```

`p` defaults to 7 if omitted. `max_degree` plays the same role as `mr_st`'s
first argument: it must be set comfortably larger than the highest internal
degree the requested `resolution_length` will reach (see "the `maxDeg`
subtlety" below).

The program prints, for each resolution step (most recent first), one row
of generator counts indexed by `t-s` (internal degree minus resolution
step) — the same `dims[degree][step]` layout `Steenrod_Op::output_resolution`
uses for the sphere.

## Verifying this against the paper's own worked example

`MinimalResolution.pdf`'s example (the one this file is "inspired from")
works out, by hand, the *first* step of exactly this resolution for
`Gamma = F_p[x_1]/x_1^5` (`p>=5`): the cokernel of `F_p -> Gamma` is
`F_p{x_1,x_1^2,x_1^3,x_1^4}`, and solving `Delta-bar(y) = ...` degree by
degree finds a **single** new cogenerator, in degree 1 (`Y={x_1}`), with
`f(x_1^2)=2x_1[x_1]`, `f(x_1^3)=3x_1^2[x_1]`, `f(x_1^4)=4x_1^3[x_1]` — i.e.
**no further generators** in degrees 2, 3, or 4.

Running `./mr_trunc 20 4 7` reproduces exactly this: the row for step 1 is
`1  0  0  0  0  ...` (one generator at `t-s=0`, i.e. absolute degree 1;
nothing else through degree 4). This is the same "check against a
known-correct case, not just a compile" discipline `docs/GENERAL_COMODULES.md`
describes for the `BP_*BP` case.

## A pre-existing bug this example exposed: `Fp_Op::inverse`

`Fp_Op::inverse` (`Fp.cpp`) previously special-cased only `p=2` and `p=3` —
the only primes the rest of this repo ever runs at — and printed
`"not implemented!"` and silently returned `0` for any other prime.
Since `p=7` is exactly the point of this example, that bug was hit
immediately (as a `std::out_of_range` crash a few steps later, once the
wrong "inverse" corrupted a Gaussian-elimination pivot). It's now
implemented generally via Fermat's little theorem, `x^(p-2) mod p`, which
is algebraically identical to the old code's special cases for `p=2` (`x^0
=1`) and `p=3` (`x^1=x`), so this is not a behavior change for any
existing (`p=3`-only) pipeline — confirmed by re-running the sphere's own
`mr_st` after the fix.

## Caveats

- Same "no comodule-axiom verification" caveat as `GENERAL_COMODULES.md`:
  nothing checks that the `delta` formula above is actually coassociative
  for a given `GammaDim`/`p` combination (it is, by construction, for any
  `p` and `GammaDim`, since it's the standard truncated-polynomial-Hopf-
  algebra comultiplication — but if you edit `delta` to try a different
  algebroid, you're back to that caveat).
- **The `maxDeg` field is a global degree truncation, not `Gamma`'s own top
  degree.** `Hopf_Algebroid::maxDeg` is used only as `ranksBelowDeg(maxDeg -
  deg_x)` (`hopf_algebroid/9.h`, `/11.h`) to bound how many algebroid basis
  elements to consider — for `Gamma`, `ranksBelowDeg` already saturates at
  `GammaDim` regardless of its argument, so the only requirement on
  `maxDeg` is that it stay at least as large as the highest degree the
  resolution reaches (exactly the same convention `mr_st`'s `max_degree`
  CLI argument follows for the infinite-dimensional dual Steenrod algebra).
  Since `maxDeg` is `unsigned`, setting it too small relative to
  `resolution_length` risks the subtraction underflowing rather than
  cleanly erroring — pick it generously.
