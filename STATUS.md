# Current Status

_Last updated: 2026-07-12_

The target-cofree rewrite of `lift.h` (per the plan below) is **implemented, fully
verified, and closed out**. `TRACE_SQUARE_CHECK` — the direct chain-map correctness test
— now runs clean across **every genuine class at every homological-degree level of the
entire resolution** (`bs=0..29`, 161/161 genuine `(bs,bb)` pairs, zero failures), not
just the one `{6-9}` chain originally tested. See `debug_notes.md` for the full
blow-by-blow; this file gives the coarse summary.

## Exhaustive verification sweep (2026-07-12)

Ran `TRACE_SQUARE_CHECK` (full level coverage) for every **genuine** class (excluding
tau-Bockstein boundary/killed elements — see "A separate, unrelated finding" below) at
every filtration `bs=0` through `bs=29` — **161 `(bs,bb)` pairs total, 161 clean, 0
genuine failures**. This is a far more conclusive confirmation than the single beta
chain checked earlier: the `φ` chain-map construction (all 3 bugs fixed this session)
holds exactly across the entire genuine-class space of this resolution.

## Active plan

`~/.claude/plans/look-at-the-top-level-smooth-llama.md` — **implemented and verified**.
Three additional real bugs were found and fixed along the way (beyond what the plan
text anticipated — see "Bugs found and fixed" below). `PHI_EXTENSION_ISSUE.md` is
**RESOLVED**.

## The resolution (math) — unchanged from before, now implemented and verified

The bug was using the **wrong universal property**. The old code stored `φ_k` by its
values on the **source's** cogenerators and reconstructed via `φ(a[y]) = a·lg(y)` — the
property of maps *out of* the cofree source, which does not determine the map. The fix:
use the **target's** cofree property. Represent `φ_k` by `φ̄_k := (ε⊗1)∘φ_k : G_k →
Z_{k+bs}` and reconstruct `φ_k = (1⊗φ̄_k)∘ρ_{G_k}`. `φ̄_k` is forced on `image(inj_k)`
(`φ̄_k(inj_k(x)) = (ε⊗1)(im_gk) =: prescribed(x)`) and free (=0) elsewhere, read off at a
τ⁰-pivot echelon reduction of `{inj_k(x)}`. Full derivation in `CLAUDE.md`.

## Bugs found and fixed this session

1. **`gens_k` (the module-generating subset) is not all of `X_{k+1}`.** `inj_{k+1}` is
   defined, and the chain condition must be imposed, on **all** of `X_{k+1}.rank()`
   (`R[k+1].Xrank`), not just the small subset used to build `G_{k+1}`'s cofree
   structure. Fixed: `lift_one_step_sum` now takes `X_rank` and loops over all of it
   (`yoneda2.cpp:596-598`).
2. **The echelon reduction's right-hand side (`prescribed`) must be transformed by the
   same row operations as the matrix, not left untouched.** Standard Gauss-Jordan
   requires reducing the augmented system `[M | b]` together; only `M` was being
   reduced. `echelon_pivots_tau0` now takes an optional `rhs` parameter and updates it
   in lockstep (`lift.h`). Also recorded as **CLAUDE.md pitfall #13**.
3. **A single forward pass over known pivots isn't enough — switched to true
   Gauss-Jordan.** Eliminating one already-known pivot column from a row can reintroduce
   nonzero content at a *different* already-known column the same pass already visited.
   Fixed by (a) looping each row's own reduction to a fixed point, and (b) retroactively
   eliminating each brand-new pivot from every earlier-finalized row too. This also let
   the earlier "backward pass" workaround be deleted — after true RREF, `prescribed[x]`
   directly *is* `φ̄(pivot[x])`.

(Note: while chasing bug 3, a `TRACE_LIFT_SUM` "duplicate pivot" observation turned out
to be a trace-scoping artifact — `lift_one_step_sum`'s `trace_budget` is a
function-static counter shared across *all* calls for the whole program run, so naive
`x=`-based `grep`s can silently mix unrelated calls. `X_rank=` was added to the trace
line to disambiguate. Bug 3 itself is a real, legitimate fix, just not for the reason
first suspected — see `debug_notes.md`.)

## The `p=112` investigation — resolved as a diagnostic-tool limitation, not a phi bug

After the three fixes above, `TRACE_SQUARE_CHECK` still reported a failure at level
k=0, p=112. Deep tracing (coaction dumps, primitivity checks, exact degree
verification) found the actual cause: **the resolution is deliberately degree-truncated
at `maxDeg=40`** (`adjoint()`, `hopf_algebroid/9.h:10`:
`total_rank = ranksBelowDeg(maxDeg - deg(generator))`), and position 112 in G₀
corresponds to the monomial `x_1^19` (degree 19), which — when it needs to act on
cogenerator 9 of G₆ (degree 22) — requires total degree 41, **one past `maxDeg`**.
Confirmed exactly: `deg(x_1^18)+deg(cogen9)=40` (boundary, valid), `deg(x_1^19)+deg(cogen9)=41`
(one past it, and cogenerator 9's block is exactly 112 positions long — the precise
boundary). This is the resolution's own data genuinely running out, not a construction
bug — the truncated model is self-consistent, and the exhaustive position sweep is
bound to eventually probe one degree past where any data exists.

**`TRACE_SQUARE_CHECK` was hardened accordingly** (`yoneda2.cpp:631+`): it now computes
`deg_beta` once and, for every failure candidate, checks whether
`deg(p) + deg_beta > maxDeg` — if so, it's a truncation-boundary artifact, logged and
skipped (scan continues) rather than treated as a stopping failure. Re-running with
this hardening, `TRACE_SQUARE_CHECK=24` (covering **every** level of the `bs=6` chain)
completes with **zero genuine failures** (20,802 boundary artifacts correctly
identified and skipped, spot-checked for sanity).

## A separate, unrelated finding: `cyc[bs]` includes tau-Bockstein boundary elements

While sweeping many `(bs,bb)` pairs, some genuinely failed at trivially low degree
(e.g. `bs=3,bb=2`, `deg(p)+deg(beta)=8`, nowhere near `maxDeg`). Traced via the
existing `TRACE_TAG_SEM` hook: these are classes where the tau-Bockstein table
(`tao_bockstein.cpp`) registers the index as a **boundary** (killed/paired element,
`tag != Invalid` in the OTHER level's table), not a genuine survivor — user confirmed
independently that e.g. `{3-2}` is mathematically not a cycle. `get_cycles()`
(`tao_bockstein.cpp:231-245`) includes both genuine cycles AND corrected boundaries in
its result map with no distinguishing marker; this is very likely intentional for
`find_cycle`/`find_cycle_sum`'s reduction algorithm (which needs entries for boundary
indices too, to cancel them via substitution), but means `yoneda2.cpp:546`'s
`cyc[bs].count(bb)` check can't tell a real class from a boundary. **Not a bug in the
resolution or in anything fixed this session** — just a narrow gap where `yoneda2.cpp`
doesn't validate that a user-supplied `bb` is a genuine survivor before using it. Safe
to ignore as long as `(bs,bb)` pairs are chosen from an independently-known-correct
list of classes (as the exhaustive sweep above did, via `TRACE_TAG_SEM`'s
`tag=INVALID/genuine-cycle` marker).

## Verification status (all PASS)

1. `φ_1(τ_1[1-0]) = {7-10}` — **PASS**.
2. `{4-6}*{6-9}` vs `{6-9}*{4-6}` — **PASS**, both cleanly `t^0{10-25}`.
3. `./yoneda2 40 30 1 1` regression — **PASS**, unchanged/clean.
4. `TRACE_SQUARE_CHECK=24 ./yoneda2 40 30 4 6 6 9` — **PASS**, clean across all 24
   levels of this beta chain (previously stopped at k=0 due to bugs 1-3, then at the
   `p=112` boundary artifact until the hardening above).

## Currently open

Nothing outstanding for this plan. Possible follow-ups, not blocking:
- The `TRACE_SQUARE_CHECK` hardening only reasons about *this* beta chain's own
  `deg_beta`; if reused for a different beta/product, the same logic applies
  automatically (recomputes `deg_beta` from that chain's own `beta_rep_pos`), so no
  further generalization should be needed.
- Temporary debug instrumentation added this session (`TRACE_ECHELON_ROW`/
  `TRACE_ECHELON_N` in `lift.h`, `TRACE_COACTION_POS`/`TRACE_COACTION_SRC_RANK` in
  `lift.h`) is harmless (env-gated) and left in place — useful enough to keep, but could
  be removed for cleanliness if desired.

## Dead ends / ruled out

- Comparing the two product directions structurally (they build different chain maps).
- Solving `φ_k`'s extension by local/degree-by-degree reasoning on `G_k`'s own coaction
  (genuinely underdetermined — the whole reason the target-cofree reframing was needed).
- Back-substitution being unnecessary ("forward elimination alone gives a valid basis,
  so it should be enough") — true for basis validity, **false** for value assignment;
  this was an actual bug, now fixed (superseded by full Gauss-Jordan, see above).
- "phi_0 preimage-independence gap" (a suspected math-level issue) — retracted, the
  actual cause of that failure (p=27) was the missing rhs-propagation bug (#2 above).
- "Duplicate pivot 1058" (x=116 and x=578 both claiming it) — retracted, a
  trace-scoping artifact, not a real bug.
- "Missing degree-truncation check causes a real phi bug at p=112" — retracted; it's a
  genuine boundary of the truncated resolution's own data, not a bug, confirmed via
  exact degree arithmetic.
