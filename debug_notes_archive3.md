## Current status (new session, 2026-07-10) — fresh investigation per handoff doc
## `~/.claude/plans/compiled-hopping-truffle.md`; three parallel read-only passes done,
## now implementing the confirmed tag_index fix and re-testing.

**Repro (from handoff doc, to reconfirm today)**: `./yoneda2 40 30 6 9 | grep '{4-6}'` and
`./yoneda2 40 30 4 6 | grep '{6-9}'` should give the same nonzero value; the former is
reportedly silently wrong, the latter hard-aborts on a `tauOper::add` "t^i + t^j"
INVARIANT BROKEN.

**Ruled in/confirmed this session:**
- `recover_inv_ind`'s prior tau^0-preference fix (lift.h:56-93) IS correctly in place —
  re-read in full, logic checks out (freezes on first clean/tau^0 singleton found per
  target column, via `inv_is_clean[k]` guard at line 68/70).
- The "no clean preimage" diagnostic in `recover_inv_ind` (lift.h:79-90) is STILL a soft
  `std::cerr` warning, not a hard abort, exactly as the handoff doc suspected. Not yet
  changed this session (user said to do plan items 2-3 first, hardening this is item 1 —
  holding off).
- **NEW confirmed bug** (separate from the tau-preimage bug fixed last session):
  `tables[s].tag_index`'s keys are level-`(s-1)` indices, not level-`s` (traced through
  `tao_bockstein.cpp`'s `make_table`/`make_pretable`/`insert`, ~lines 189-201, 85-153).
  But `yoneda2.cpp:567`'s output filter does
  `if (tables[s].tag_index.count(a)) continue;` where `a` is a level-`s` index — comparing
  mismatched levels. Correct fix: check `tables[s+1].tag_index` instead (keys = level-`s`
  "tags" that die feeding into level s+1), guarding `s+1 >= tables.size()` (top level: no
  data means nothing can kill it, so don't filter). There's already a `TRACE_TAG_SEM`
  diagnostic block (yoneda2.cpp:342-363) independently confirming this exact mismatch and
  flagging line 567 as "(possibly buggy)".
- Located the mechanical crash site for `4 6`: `lift_one_step` (lift.h:190-244) builds
  `main_term` divided by `tau_inv = inverse(inv_tau_prev[x])` (line 198-199) and a
  `correction` summed from other cogenerators in the same `inj_src` row (lines 202-216),
  then `add`s them (line 218) — this is where `tauOper::add`'s "t^i+t^j" abort fires if two
  contributions land at the same index with different tau exponents.

**Hypothesis to check today**: the `4 6`/`6 9` failures are downstream fallout from the
tag_index bug (a spurious "surviving" class like `{4-5}` participating in the resolution
with a tau-exponent that doesn't belong, corrupting `inj_src`/`phi_prev` entries that later
collide in `lift_one_step`'s correction-term sum) — NOT a separate remaining bug in
`lift_one_step` itself. Not yet verified; plan is to fix tag_index, rebuild, clear caches,
rerun the two repros, and see if they resolve. If not, trace `lift_one_step` directly with
`TRACE_LIFT` to find the exact colliding tau-exponents/positions.

**Next concrete action**: implement the `tables[s+1].tag_index` fix in yoneda2.cpp around
line 565-572, rebuild, clear `<N>_yoneda2_*` caches, rerun repros.

## Update: tag_index fix implemented, does NOT fix the `4 6` abort by itself

Implemented the fix at yoneda2.cpp:565-573 (`tables[s+1].tag_index.count(a)` with an
`s+1 < tables.size()` guard). Rebuilt cleanly. Cleared `40_yoneda2_4_6_phi*` and
`40_yoneda2_6_9_phi*` caches, reran `./build/yoneda2 40 30 4 6`: **still hard-aborts** with
the identical `INVARIANT BROKEN: add(t^2, t^1)` message. So the tag_index bug is real and
now fixed, but it is NOT (solely) the cause of the `4 6` abort — the earlier hypothesis
("downstream fallout from tag_index") is RULED OUT as the sole cause; there's a genuinely
separate issue in `lift_one_step`'s correction-term math.

**Localized the crash precisely**: `phi_beta[0..6]` load from stale-but-valid cache
(`phi_exists` true), so the only `lift_one_step` call that actually executes fresh is
`k=6` (building `phi_beta[7]` = the level-7 lift). Ran with `TRACE_LIFT=500`
(lift.h:187/224-242): the crash happens partway through cogenerator `j=12` of level 7's
`ncog=26` cogens (trace budget still had ~488 left, so it's not a budget cutoff) — the
last successfully-traced cogenerator is `j=11` (`x=1116, cog_pos=2089`, `inj_row.size=1`,
i.e. inj_row was PURE for j=11, no correction needed there). So the abort is inside j=12's
correction-term loop (lift.h:206-216), before that iteration's trace print at line 224-237
fires — need per-term tracing inside the loop itself (not just the after-the-fact summary)
to see which two tau-exponents collide. **In progress**: adding targeted std::cerr prints
inside the correction loop (lift.h:206-216), gated on the existing `trace_budget>0`
condition, showing each `tm.ind`/`tm.coeficient`/`phi_pos_cogens` term and its running
`correction` value before/after each `add()` call.

## Update: exact collision fully pinpointed with concrete trace data

Added `TRACE_LIFT_PRE`/`TRACE_LIFT_CORR`/`TRACE_LIFT_FINALADD` prints (lift.h, gated on
existing `trace_budget>0`) to `lift_one_step`, right before/inside/after the
correction-term loop (lift.h:202-231ish). Rebuilt, reran with `TRACE_LIFT=500`. Full
concrete trace of the crashing cogenerator, level 7 (k=6 -> building phi_beta[7]),
`j=12`, `x=1146` (X_7 index), `cog_pos=2141` (its position in G_7):

```
full_inj_row: 2091^t1 2141^t0        (inj_7(x=1146): self term at 2141^t0, correction term at pos 2091^t1)
main_term(im_gk, pre-tau-div): 965^t2   tau_inv=t^0   (inv_tau_prev[1146]=0, i.e. a CLEAN
                                          preimage was chosen -- no twisted division happened here)
  -> main_term_cogens (reindexed): 22^t2
TRACE_LIFT_CORR: term.ind=2091 term.coef=t^1  scaled: 21^t1 22^t1
  (phi_pos = cofree_adjoint_row(lg,...,2091,...) reindexed to cogens gives {21^t0, 22^t0}
   BEFORE scaling; inj_7's t^1 correction coefficient scales it up to {21^t1, 22^t1})
FINALADD: main_term_cogens: 22^t2  |  correction: 21^t1 22^t1
-> add(main_term_cogens, minus(correction)) at index 22: t^2 + t^1  -->  ABORT
```

So concretely: two independent derivations of the coefficient at target cogenerator 22
disagree — the direct chain-condition path (`im_gk` via position 965, reindexed) says
`t^2`, while the correction-term path (`phi_beta[7]` at position 2091, scaled by inj_7's
own `t^1` correction coefficient) says `t^1`. Since `tau_inv=t^0` here, this is NOT an
artifact of the tau-division-on-torsion issue fixed last session (no twisted preimage was
used) — it's a genuine disagreement between the "main" and "correction" halves of the
`lift_one_step` formula itself, i.e. a distinct bug from either previously-identified one.

**Ruled out**: this is not the tag_index bug's fallout (confirmed above: fix applied,
abort persists identically). Not a tau-division-on-torsion artifact (tau_inv=t^0 here,
no division occurred).

**Not yet determined**: which side (main-term chain condition, or the correction-term's
`phi_beta[7]`-at-position-2091 computation) is actually wrong, or whether both are
"correct" under their own local logic but the formula combining them is structurally
incomplete (e.g. missing a further reindexing/twist step when `inj_src`'s correction
coefficient itself carries a tau twist that should interact multiplicatively with
`phi_beta`'s own internal tau-bookkeeping in a way not currently modeled). This needs
deeper case-by-case math verification (e.g. hand-computing what cogen 22's value SHOULD
be here) rather than more mechanical tracing — stopping here to report to the user and
get direction before proposing/attempting any fix, per standing instructions not to
change strategy without asking.

**Instrumentation added this session, all still in tree, env-gated on `TRACE_LIFT`**:
`TRACE_LIFT_PRE`, `TRACE_LIFT_CORR`, `TRACE_LIFT_FINALADD` (lift.h, inside
`lift_one_step`, all print only when the existing `trace_budget>0`).

## New direction (user-proposed, plan approved): tauPoly is the wrong data model —
## switch phi_beta's coefficients to a genuine F2[tau] polynomial type

User's insight: Ext over the C-motivic Steenrod algebra is genuinely a module over
`F2[tau]`, so a basis coefficient legitimately CAN be a real sum like `tau+tau^2` — the
`t^2 vs t^1` clash found above may not be a bug at all, just `tauPoly`'s "single monomial
only" simplification breaking down. Checked `tao_bockstein.cpp`/`.h` per the user's
suggestion: **no dedicated multi-term polynomial type exists there** — every coefficient
is still plain `tauPoly`. BUT the exact reusable template already exists elsewhere in this
codebase: `polynomial<base_ring>`/`PolynomialOp<base_ring>` (`polynomial/2.h:60-65`,
generic `poly<exponent,base_ring>=vectors<exponent,base_ring>` + `PolyOp` with real
`add`/`multiply`/`invertible`/`inverse`) is ALREADY instantiated with `base_ring=tauPoly`
to build `motSteenrod` (`mot_steenrod.h:43`) — instantiating the same template with
`base_ring=F2` gives genuine `F2[tau]` for free.

**Key de-risking finding**: `PolyOp::invertible`/`inverse` (polynomial/6.h:26-33) only
support the degree-0/constant case — which matches the ALREADY-established invariant that
`tau_inv` in `lift_one_step` is always `t^0` in a correct computation (clean preimage
always chosen/exists, per last session's `recover_inv_ind` fix). So no Laurent/negative-
exponent generalization is needed for the new type to work as a drop-in for the "divide by
tau^n" step.

**Approved plan** (full plan at
`~/.claude/plans/look-at-claude-plans-compiled-hopping-tr-synthetic-puffin.md`, scope
confirmed with user as "phi_beta only", NOT touching resolution's own inj/qut/differential
matrices):
1. Harden `recover_inv_ind`'s diagnostic to a hard abort (lift.h:79-90) — still not done,
   now doubly important as a safety net for the tau_inv==t^0 invariant.
2. Add `typedef polynomial<F2> tauPolySum;` + `PolynomialOp<F2> tauPolySum_oper` (mirroring
   `Z2.h:22`'s `Fp_Op F2_opers` pattern) in mot_steenrod.h/.cpp, plus `tauPoly->tauPolySum`
   lift helpers mirroring the existing `lift(F2 x)->tauPoly` pattern (mot_steenrod.h:135-139).
3. Change `lg`/`phi`/`main_term`/`correction`/`combined` in lift.h's `cofree_adjoint_row`/
   `lift_first_step`/`lift_one_step` to `tauPolySum`, with small bridging helpers wherever a
   `tauPoly`-valued resolution matrix (`qut`,`inj`) is applied to a now-`tauPolySum`-valued
   vector (apply termwise per input monomial, re-sum with `tauPolySum_oper.add`).
4. Change `phi_beta`/`M_beta`'s type in yoneda2.cpp, update `fmt()` to print genuine sums,
   lift `beta_rep`/`beta_rep_pos` to `tauPolySum` before use.
5. Clear ALL `<N>_yoneda2_*_phi*` caches (binary format changes), rebuild, rerun `1 1`
   (should be unchanged) and the `4 6`/`6 9` pair (should now agree, possibly as a genuine
   multi-term value — that would be an expected correct outcome, not evidence of a new bug).

**Step 1 done**: hardened `recover_inv_ind`'s diagnostic to a hard `abort()` (lift.h:79-93,
matching `tauOper::add`'s `fprintf`+`abort` style). Added `<cstdio>`/`<cstdlib>` includes.
Rebuilt cleanly (`cmake --build build --target yoneda2`).

**Step 2 done**: added `typedef polynomial<F2> tauPolySum;` + `extern PolynomialOp<F2>
tauPolySum_oper;` + `liftToPolySum(tauPoly)`/`liftToPolySum(vectors<...,tauPoly>)` helpers
to mot_steenrod.h (declarations) and mot_steenrod.cpp (definitions, backed by a new
`Fp_Op tauPolySum_F2_ops(2)` instance, mirroring `Z2.h`'s `Fp_Op F2_opers` pattern).
Needed to add `Fp.cpp` to the `yoneda2` CMake target's sources (CMakeLists.txt:240) since
it previously didn't need `Fp_Op`'s definitions (link failure caught this immediately:
`undefined reference to Fp_Op::unit(int)`/`Fp_Op::Fp_Op(int)`). Rebuilt clean.

**Step 3 in progress**: added new `_sum` variants in lift.h (`apply_ring_matrix_to_sum`,
`cofree_adjoint_row_sum`, `lift_first_step_sum`, `lift_one_step_sum`) that use `tauPolySum`
for the accumulator while resolution data (`qut`/`inj`) stays `tauPoly`. Also had to add
`find_cycle_sum`/`leading_term_sum`/`tau_valuation_sum` to tao_bockstein.h/.cpp (not
originally called out in the plan by name, but necessary: `M_beta`'s cogen-projected
output and the final `img`/`find_cycle` step are downstream of phi_beta and also need to
become tauPolySum-valued for the pipeline to type-check end to end).

**Important side effect found**: `mot_steenrod.cpp` is shared by several OTHER CMake
targets (`motTab`, `mot_mult`, `dump_gens`, `test_lift`, `tauBoc`, and `yoneda` -- the
out-of-scope one) that didn't previously link `Fp.cpp`. Adding `Fp_Op` usage there broke
their link step. Fixed by adding `Fp.cpp` to all of those targets' CMakeLists.txt sources
EXCEPT `yoneda` (left alone since editing anything related to yoneda.cpp, even its build
config, felt too close to the "don't touch yoneda.cpp" hard constraint -- flagging this
to the user rather than deciding unilaterally; `yoneda`'s build is presumably now broken
until either Fp.cpp is added to its CMake entry or the user says otherwise).

**Step 3+4 complete, builds cleanly end to end.** Finished updating all yoneda2.cpp call
sites (`lift_one_step_sum`, `M_beta`/`phi_beta` retyped to `tauPolySum`, `fwd_product`
uses `find_cycle_sum`, all TRACE_* print sites fixed to use `fmtCoef()`). Rebuilt clean
(only pre-existing/benign warnings). Cleared all 134 stale `<N>_yoneda2_*_phi*` cache
files (binary format changed) and reran fresh.

**Results**:
- `./yoneda2 40 30 1 1` (previously-passing regression check): UNCHANGED, still all clean
  single-term `t^0{...}` results, e.g. `{2-1} -> t^0{3-1}`. No regression.
- `./yoneda2 40 30 4 6 | grep '{6-9}'`: **no longer aborts** -- now gives
  `{6-9} -> t^0{10-25}` (a clean single-term result, not even a genuine multi-term sum in
  this case).
- `./yoneda2 40 30 6 9 | grep '{4-6}'`: still gives `{4-6} -> 0`.
- **Confirmed via direct single-product mode too** (not a table-dump filtering artifact):
  `./yoneda2 40 30 4 6 6 9` => `{4-6}*{6-9} = 0`; `./yoneda2 40 30 6 9 4 6` =>
  `{6-9}*{4-6} = t^0{10-25}`. These should be equal (graded-commutativity, no signs in
  char 2) but still disagree -- **the original bug (asymmetric zero vs nonzero) is NOT
  yet fully fixed**, though the hard-abort crash that was blocking `4 6` entirely IS fixed,
  and the `tag_index` level-mismatch bug is a real, confirmed, independent fix.

**Status**: reporting progress to user now rather than continuing to guess. Three
confirmed/fixed bugs this session (tag_index mismatch, recover_inv_ind hard-abort
hardening, the tau-sum representation gap that caused the `INVARIANT BROKEN` crash), but
the core "{s-a}*{t-b} should equal {t-b}*{s-a}" asymmetry from the original bug report
persists in a new, crash-free form. Next step (not yet started): investigate why
`{4-6}*{6-9}` computes as exactly 0 via `TRACE_LIFT`/`TRACE_PROD_S`/`TRACE_PROD_A`
tracing on the `bs=6,bb=9` computation, analogous to the investigation done earlier this
session for the `4 6` abort.

## Result of raw-phi-with-offsets dump: no offsets actually appear in this range

Ran the raw-phi dump (both programs, `TRACE_PHI_TABLE=12`): in this range (hom degree 0-4,
internal degree<=12, beta={1-1}), phi's raw output is ALWAYS either exactly 0 or a single
PURE cogenerator term in both implementations — never an offset term. So the earlier
`M_beta`-projected tables weren't hiding anything; raw==projected everywhere here. The
`{2-1}` divergence (yoneda2=0 vs yoneda-old=t^0{3-1}) stands, still pure on both sides.

## ROOT CAUSE CONFIRMED (2026-07-06): the real bug is a divide-by-tau on torsion, caused by
## recover_inv_ind picking a tau^n-twisted preimage when a guaranteed tau^0 one exists

User's diagnosis (confirmed correct): the chosen preimage `{1-0}+4` gives
`qut_1 = t^1 * x467`, so `lift_one_step`'s formula computes `phi({2-1}) = tau^{-1} *
phi_prev(applied through the chain)`, i.e. `tau^{-1} * (t^1 * {3-1}) = {3-1}` IF phi were
linear over the torsion-free part — but `{3-1}` is itself `tau`-torsion (killed by
multiplication by tau, since it's a genuine Ext class with a nontrivial tau-Bockstein
differential nearby, per the `tau_table` structure investigated earlier this session), so
the correct value is silently divided away to 0. This is a real, previously-unidentified
bug in the MAIN strategy of `lift_one_step`/`recover_inv_ind` (`lift.h`) — not the
correction-terms issue this session started with, and not something the earlier "ruled
out" appendix note (#2, checking for `inv_tau!=0 AND im_prev nonempty` cases) would have
caught, since that check was scoped to a different product's dependency chain only.

**Checked `matrices/6.h`'s `make_quotient`/`quot_index` (the code that actually builds
`qut` inside `resolvor_modeled`) to answer the user's question: does the construction
ALWAYS guarantee a tau^0 (untwisted) preimage exists for every X_{i+1} target?** Answer:
**YES, unconditionally, by construction** — `quot_index` (lines 3-23) partitions G_i's
positions into "pivot" positions (`inj_index`, the leading terms of the echelon-reduced
`indj`) and "complement" positions (`inverse_ind`, everything else, numbered `0,1,2,...`
in scan order to become X_{i+1}'s basis order). `make_quotient`'s SECOND loop (lines 40-41):
`result->insert(inverse_ind[i], moduleOper->singleton(i))` — this sets, for EVERY
complement position, `qut(complement_position) = singleton(i) = EXACTLY tau^0 * e_i`, no
twisting, no exceptions. Only the PIVOT positions (first loop, lines 31-38) can get
arbitrary tau-twisted rows (from the echelon reduction of `indj`). So **every X_{i+1}
basis element has a guaranteed clean tau^0 preimage: its own designated complement
position** — the twisted ones (like position 4, a pivot-adjacent artifact) are just
ADDITIONAL (non-canonical) singleton rows that happen to exist alongside the guaranteed
clean one, from unrelated echelon-reduction coincidences.

**The bug is entirely in `recover_inv_ind`'s selection heuristic** (`lift.h`): it scans
positions 0..F_rank-1 in raw order and keeps the FIRST singleton row found per target,
with no preference for `inv_tau==0`. Since twisted rows can appear at lower raw positions
than the guaranteed-clean complement position (as happened here: position 4 < position
976), `recover_inv_ind` can and does pick a torsion-twisted preimage even though a clean
one always exists elsewhere. **Fix direction (not yet implemented, awaiting user go-ahead
per the "ask before changing lift.h" standing instruction)**: change `recover_inv_ind` to
prefer/require `inv_tau==0` rows — e.g. two-pass scan (first pass only accepts tau^0
singletons; only fall back to a twisted one if literally none exists, which per the above
should never actually happen) — or, more robustly, thread the `inverse_ind`/complement-
position data out of `resolvor_modeled` directly instead of rediscovering it by scanning
`qut`. Reported this full finding to the user; have not modified lift.h or recover_inv_ind
yet.

## MAJOR FINDING: recover_inv_ind's "pick the first/lowest singleton position" heuristic
## looks like it picked the WRONG preimage for {2-1}, contradicting the earlier "proven
## safe by exactness" conclusion (2026-07-06)

User's hypothesis: {2-1} "comes from {1-0}" (per the inv_inj table: preimage =
`t^1 * {1-0}+4`) but should intuitively come from {1-1} instead. Added `TRACE_QUT_SWEEP_*`
(env-gated: `_LEVEL`, `_COG`, `_TARGET_X`) sweeping every position in a given cogenerator's
own cofree summand, checking `qut_h`'s row for a target X index. Swept `{1-1}`'s summand
(positions 974..1870) for x=467 (the X_2 index behind `{2-1}`): **FOUND A HIT** —
`qut_1({1-1}+2) = t^0 * x467` — a PURE tau^0 singleton, directly in `{1-1}`'s own summand
(at monomial offset 2 = x_1^2). This is a perfectly valid alternative preimage of x=467,
and arguably the "expected" one (cleaner: tau^0, and living in {1-1}'s summand as the user
intuited) — but `recover_inv_ind` (`lift.h`) iterates positions 0..F_rank-1 in raw order
and keeps only the FIRST singleton row found for each target; since position 4 (in
`{1-0}`'s summand, giving `t^1 * x467`) comes before position 976 (`{1-1}+2`), position 4
is the one actually used by `lift_one_step` for the `{2-1}` computation — NOT position 976.

**This directly contradicts finding #6 from the original plan's appendix** ("recover_inv_
ind's 'any singleton row' preimage choice — proven false [as a concern], both
mathematically and empirically... any two valid singleton preimages differ by an element
of image(inj), and composing with the next qut in the chain condition kills that
difference exactly"). That earlier "proof" was checked only for the `{3-1}*{1-1}`
dependency chain, not for this specific case, and now we have direct empirical reason to
suspect it does NOT hold here — using position 976 as the preimage instead of position 4
could plausibly give a nonzero downstream `lg[{2-1}]` instead of 0. **Not yet verified
computationally** — next concrete step: directly compute and compare
`phi_beta[0].find(4)` vs `phi_beta[0].find(976)`, and what `lift_one_step`'s formula would
give for `lg[{2-1}]` under EACH choice, to see whether they actually agree (as claimed) or
differ (which would locate a real bug in `recover_inv_ind`'s tie-breaking, not just in the
correction-terms code touched so far). Reported the sweep finding to the user; this
computation is in progress.

## offset -> monomial name table done, results reported to user (2026-07-06)

Ran `TRACE_MON_NAMES=40`. Uses `::output(e)`'s hardcoded default name "x" (not "xi" as I'd
planned — `MOP.output_monomials()` calls the free function without a name override, and I
didn't modify shared code to change that), so names are `x_1^k`, `x_2^k`, `x_3^k` etc.
(presumably Milnor basis generators xi_1, xi_2, xi_3, ...). Applying this to the earlier
inv_inj/inj-correction table: e.g. `{3-2}`'s correction `t^1 * {3-1}+2` means "t^1 times
x_1^2 acting on cogenerator {3-1}"; `{4-1}`'s preimage `t^2 * {3-0}+8` means "t^2 times
x_1^2x_2^1 acting on cogenerator {3-0}". Reported to user; standing by for next direction.

## New ask: offset -> human-readable monomial name table (2026-07-06)

"Offset" = local index within a cofree comodule summand = index into `MotSteenrodOp`'s
private `mon_array` (confirmed via `mot_steenrod.cpp:166,176`: `algebroid2vector` places a
monomial `e` at absolute position `mon_index[e]+shift`, i.e. `mon_index[e]` IS the offset,
and `mon_array[offset]` is the corresponding `exponent`). `mon_array`/`mon_index` are
PRIVATE members of `MotSteenrodOp` — no direct external access — but `MOP.output_monomials()`
is already public and prints one line per `mon_array` entry in index order (via `::output(e)`
from exponents.cpp, format like `x_1^2x_3^1`), so offset->name can be recovered by parsing
that string's lines by position, with NO source changes to mot_steenrod.h needed. This
table is universal (one shared dual-Steenrod-algebra basis for every cofree summand, not
per-cogenerator), so only needs to be printed once, covering the offsets actually seen so
far (max offset seen = 38, from `{4-5}`'s preimage). Implementing as `TRACE_MON_NAMES` in
yoneda2.cpp (added `#include <sstream>` already). Using name "xi" instead of default "x"
for interpretability as dual-Steenrod-algebra Milnor generators.

## inv_inj + inj-corrections dump complete, results reported to user (2026-07-06)

Implemented and ran `TRACE_INV_INJ=12` in yoneda2.cpp. Full output captured (s=1..4, deg<=12,
raw resolution data only — not program-specific). Key facts for the record:
- s=1,2: ALL preimages are pure singletons of the form `{s-1 - 0}+<offset>` (i.e. every
  X_s basis element's canonical qut_{s-1}-preimage lives in cogenerator 0's own summand at
  levels 1,2) and every inj_s row is pure (no corrections) in this range.
- s=3: `{3-2}` (x=776) is the first impure inj row: correction = `t^1 * {3-1}+2`. Its own
  preimage is pure: `t^0 * {2-0}+9`.
- s=4: `{4-2}` (x=685) impure, correction `t^1 * {4-1}+2`; `{4-5}` (x=1195) impure,
  correction `t^1 * {4-4}+2`. Preimages for both are pure singletons.
- Notably `{4-1}` (x=395) has preimage `t^2 * {3-0}+8` — a tau^2 twist, the highest tau
  exponent seen in this range (matches the appendix note about recurring inv_tau!=0
  "curiosity" from earlier sessions).
Reported the full table to user; no further action taken yet, standing by for direction.

## New ask: dump inv_inj (canonical qut-preimage) + inj corrections, not phi (2026-07-06)

User wants, for each cogenerator {s-i} of level s: (1) `recover_inv_ind`'s canonical
preimage in G_{s-1} (the position p with `qut_{s-1}(p) = tau^n * e_x`, x=gens[s][i]),
described in `{cog}+offset` notation; and (2) `inj_s(x)`'s correction terms (excluding the
self term), same notation — i.e. exactly the two raw ingredients `lift_one_step` combines,
laid out explicitly per cogenerator. Implementing as `TRACE_INV_INJ` (env-gated) in
yoneda2.cpp only (raw resolution data + `recover_inv_ind`, both from lift.h/R[], nothing
program-specific — shared ground truth for both implementations, no need to duplicate in
yoneda.cpp). Reuses the `describe_pos`-style helper already used for inj/qut and phi dumps.

## New ask: dump RAW phi (with offset terms, not just cogen-projected M_beta) (2026-07-06)

User wants the raw phi values (like the inj/qut dump's `{h-cog}+offset` notation) rather
than the cogen-only-projected `M_beta`/`build_M` table from before. Key realization:
`build_M`(yoneda2)/`build_M_at`(yoneda old) both DROP any term of raw phi landing on a
non-cogenerator position (`find_index`-filtered) when producing the M_beta table — so the
earlier phi-table dump could have silently hidden real "offset" terms. For yoneda2,
`phi_beta[h].find(pos)` IS already the raw position-space row (no change needed, just print
with the same `describe_pos` naming used for inj/qut). For yoneda(old), `psi_beta` is
ALREADY collapsed to cogen-index space at storage time (via `proj_cogens` inside
`lift_step`), so to see offsets there need the PRE-`proj_cogens` value, i.e. call
`apply_chainmap(singleton(pos), R[h].F, psi_beta[h], R[h+bs].F, MOP)` directly (this is
exactly `fa` inside `build_M_at`, before its `proj_cogens` call). Implementing: reuse/
duplicate the `describe_pos` helper (currently scoped inside `TRACE_INJ_QUT`'s block in
yoneda2.cpp) into the `TRACE_PHI_TABLE` blocks of both files, printing the raw row
alongside the existing cogen-projected row for each {h-a} in range.

## inj/qut dump implemented and run (2026-07-06)

Added `TRACE_INJ_QUT=<maxdeg>` (env-gated, `yoneda2.cpp`, right after "resolution loaded")
dumping `inj_h`/`qut_h` for h=0..4, restricted to cogenerators/X_h basis elements with
internal degree<=maxdeg, human-readable (`describe_pos` names a G_h position as `{h-cog}`
if it IS a cogenerator, else `{h-owner}+<local offset>` for a non-cogen element inside
owner's cofree summand; qut's X_{h+1} targets are labeled with their cogen name via a
`gens[h+1]` reverse-lookup when known). This is raw resolution data, identical in both
programs (same on-disk files) — dumped once from yoneda2.cpp only, not duplicated.

Confirmed inj_3/inj_4's known impure rows read exactly as expected, e.g.
`inj_4({4-2}) = t^1 * {4-1}+2  +  t^0 * {4-2}` — a genuine A_C (dual-Steenrod) monomial
(local offset 2 within cogenerator {4-1}'s own cofree summand) contributing at tau^1,
plus the tau^0 self term — i.e. exactly the "A_C elements multiplied by lower-degree
generators" pattern the user expected to see. Reported full inj_0..inj_4/qut_0..qut_4
dump to user. Not yet asked to act on anything further; standing by.

## Progress: corrected phi-table comparison found a REAL divergence at {2-1}; now adding
## inj/qut dumps (2026-07-06)

Redid the phi comparison properly (both programs projected to target COGEN-INDEX space via
`fmt()`, not mixing position-space vs cogen-index-space as before): added `TRACE_PHI_TABLE`
(env-gated, in both `yoneda2.cpp` after `M_beta` is built, and `yoneda.cpp` right after its
`M_beta` build using a RAW `apply_chainmap`+`proj_cogens` — no cycle correction — to match
yoneda2's uncorrected `M_beta`/`build_M`). Ran both at bs=1,bb=1 (beta={1-1}),
`TRACE_PHI_TABLE=12` (internal degree cutoff), homological degrees 0-4. **Both agree
through `{2-0}`, then diverge at `{2-1}`: yoneda2 gives `0`, yoneda(old) gives
`t^0{3-1}`.** Everything below that in yoneda2's table is also wrongly zero (consistent
with a single upstream failure cascading forward, matching the earlier `qut_1(pos 975)==
empty` trace from before the position/cogen-index mixup was caught). This is confirmed as
a REAL divergence (not an artifact) since both dumps now use the identical representation.
User confirmed this is testing `Hom_A(F2[tau],-)` applied to beta={1-1}'s lift, and asked
to also dump `inj_h`/`qut_h` (raw resolution matrices) in the same range (hom degree 0-4,
internal degree <=12) in human-readable form, expecting to see some A_C elements at lower
degrees. **In progress**: adding an env-gated dump of `R[h].inj`/`R[h].qut` restricted to
cogenerator positions/X_h basis elements with degree<=12, to `yoneda2.cpp` (this is raw
resolution data — loaded identically in both programs, so a single dump suffices; not
duplicating into yoneda.cpp). Not implemented yet as of this note.

## CORRECTION (2026-07-06): the "first divergence" finding below was a false alarm —
## comparing apples to oranges, not a real bug in yoneda.cpp. User flagged (correctly)
## that I might have misunderstood yoneda.cpp's conventions before I reported this as a
## bug; re-reading confirmed the user was right to be suspicious.

**yoneda.cpp's `psi_beta` and yoneda2.cpp's `phi_beta` store values in DIFFERENT spaces**,
and I compared them as if they were the same:
- `yoneda2.cpp`'s `phi_beta[lev]` stores values as **absolute positions** in the target
  cofree comodule (built via `cofree_adjoint_row` -> `HA.algebroid2vector(..., F_tgt.
  position_of_gens[ind])`, i.e. real G-positions).
- `yoneda.cpp`'s `psi_beta[lev]` stores values as **cogenerator INDICES** of the target
  (see `apply_chainmap`'s own comment: "psik = psi.find(tm.ind); // psi([k]) : tauPoly vec
  over cogens(Ftgt)", and `lift_step`'s `psi[i+1].insert(x, proj_cogens(dfv, ...))` —
  `proj_cogens` explicitly converts to cogen-index space before storing).

So `psi_beta[0].find(0) = {cogen-index 1: tau^0}` is actually EXACTLY CORRECT under
yoneda.cpp's own convention — it's directly `beta_rep` (already in cogen-index space),
inserted with no conversion needed, because psi never stores raw positions. My earlier
claim ("yoneda.cpp forgot to convert cog-index to position, that's the bug") compared
yoneda.cpp's cogen-index value against yoneda2's position value and called the mismatch a
bug — but the mismatch is expected/correct given the two programs use different internal
representations for the same mathematical quantity. **Retracted.** No conclusion should be
drawn from that comparison as previously stated.

## New task in progress: print phi's action on each cogenerator (hom degree 0-4, internal
## degree <=12) for BOTH programs in a common, human-readable format, to properly compare
## apples to apples this time (target values expressed as target COGENERATOR names, not
## raw internal representation) — reusing the existing per-program cogen-index projection
## each program already has (`M_beta`/`build_M` in yoneda2.cpp uses position->cogen-index
## reindexing already; yoneda.cpp's `build_M_at` uses `apply_chainmap`+`proj_cogens` but
## ALSO applies an extra "cycle correction" step beyond raw phi — for a fair comparison of
## raw phi action, need the UNCORRECTED apply_chainmap+proj_cogens result from yoneda.cpp,
## not the corrected M_beta). Internal degree per cogenerator = `generators.degree[a].deg
## - h` (h = homological/filtration level), per the convention already used in
## `output_generators` (tao_bockstein.cpp). Using the same bs=1,bb=1 test case as before
## for continuity. Implementation in progress: adding env-gated dump blocks to both files.

## FOUND: first phi/psi divergence, and it's a clear bug in yoneda.cpp (2026-07-06)

Dumped `phi_beta`/`psi_beta` (via `TRACE_PHI_DUMP`, see below) for both programs at
bs=1,bb=1, joined by (level,position) key, sorted ascending. **First divergence is at
level=0, position=0** (the very first value either program ever inserts):
- `yoneda2` (new): `phi_beta[0].find(0) = {974: tau^0}` — CORRECT. `beta_rep = cyc[1].at(1)
  = {cogen-index 1: tau^0}` and yoneda2 explicitly converts cog-index -> position via
  `beta_rep_pos.push({R[bs].F.position_of_gens[tm.ind], tm.coeficient})` before use
  (`yoneda2.cpp` ~line 218-220), giving position 974 = `R[1].F.position_of_gens[1]`
  (matches independently-confirmed cog1-of-G_1 position from earlier TRACE_LIFT output).
- `yoneda` (old): `psi_beta[0].find(0) = {1: tau^0}` — **WRONG**. `yoneda.cpp`'s
  initialization (~line 435-438) does:
  `psi_beta[0].insert(p, (p==cog0) ? beta_rep : ...)` — it inserts `beta_rep` **directly**,
  never converting it from cogenerator-INDEX space (`{cogen 1: tau^0}`) to POSITION space.
  So it stores raw cogen-index `1` as if it were position `1` in G_1 — but position 1 is a
  1coaction-generated (non-cogenerator) element in the SAME summand as G_1's cogenerator 0,
  totally unrelated to cogenerator 1 (real position 974). This is a plain reindexing bug:
  yoneda.cpp is missing the exact conversion step yoneda2.cpp does correctly.

This single wrong seed value (`beta` planted at the wrong position before the very first
lift step) propagates through yoneda.cpp's entire computation from level 0 onward, and is
the most likely explanation for both its wrong answers on this case AND its later
`INVARIANT BROKEN: add(t^2, t^1)` crash (garbage position data eventually colliding with
mismatched tau exponents in its own tau-echelon solver). **This is a bug in yoneda.cpp,
not in yoneda2.cpp/lift.h** — yoneda2's level-0 value is exactly the expected/correct
`beta_rep_pos` conversion. Reported to user; per explicit instruction NOT debugging
yoneda.cpp further or changing the plan — this was purely the requested diff task.

**Debug/dump methodology** (all instrumentation still in tree, env-gated, harmless):
added `TRACE_PHI_DUMP=<maxlev>` to both `yoneda2.cpp` (dumps `phi_beta[lev]`, called after
each level is built) and `yoneda.cpp` (dumps `psi_beta[lev]`, same pattern), identical
output format `PHIDUMP lev=<L> pos=<P> :: <ind>^t<coef> ...` via unbuffered `std::cerr` (so
yoneda.cpp's later crash doesn't lose earlier dumps). Ran both at bs=1,bb=1 with
`TRACE_PHI_DUMP=3`, cleared caches (`40_yoneda2_1_1_phi*`, `40_yoneda_te_1_1_phi*`); joined
the two dumps by (level,position) key in Python and diffed in ascending key order to find
the first mismatch (see above). yoneda2 has vastly more nonempty positions (3396 dumped
rows vs. yoneda's 96 by level 3) — expected, since `cofree_adjoint_row`-based extension in
lift.h reaches many more positions; not itself evidence of a bug.

## Progress update: TRACE_PHI_DUMP instrumentation added to both files, about to build+run

Added `TRACE_PHI_DUMP=<maxlev>` (env-gated) to both `yoneda2.cpp` (dumps `phi_beta[lev]`,
called `dump_phi(0)` right after level-0 build and `dump_phi(k+1)` at the end of each
k-loop iteration) and `yoneda.cpp` (dumps `psi_beta[lev]`, same call pattern around its
`lift_step` loop). Both print identical format: `PHIDUMP lev=<L> pos=<P> :: <ind>^t<coef>
...` for every nonempty row, using `std::cerr` (unbuffered, survives yoneda.cpp's later
crash). Since both programs load the same on-disk `<pre>mot_gens{i}`/`<pre>mot_maps{i}`
files, level and position numbers are directly comparable between the two dumps. Next:
build both, clear phi/psi caches, run both with `bs=1 bb=1` (the crash case) and
`TRACE_PHI_DUMP` set to a small level cutoff, then diff the two dumps line-by-line to find
the first (lev,pos) where their stored values disagree.

## Current status (session 6 continued #2, 2026-07-06) — user has decided yoneda.cpp is
## simply wrong (its own tau^i+tau^j crash is taken as evidence of that, not something to
## chase), and does NOT want yoneda.cpp debugged. New ask: dump phi (yoneda2's phi_beta)
## and psi (yoneda.cpp's psi_beta) side by side in a low level range and find the FIRST
## position where their values actually differ -- a direct empirical diff, nothing more.
## Explicitly told not to change the plan or go out of scope beyond this dump+diff task.

Plan: add a small env-gated dump (e.g. `TRACE_PHI_DUMP`) to BOTH yoneda2.cpp (dumping
`phi_beta[k]` rows) and yoneda.cpp (dumping `psi_beta[i+1]` rows) for levels 0..3ish,
using IDENTICAL position/level numbering (both load the same on-disk `<pre>mot_gens{i}`/
`<pre>mot_maps{i}` files, so positions are directly comparable). Since yoneda.cpp crashes
partway through building psi_beta for bs=1,bb=1 (session-6-continued finding above), use
std::cerr (unbuffered) so early-level dumps survive even if it later aborts. Run both with
bs=1,bb=1, cleared caches, diff the dumps level-by-level/position-by-position, report the
first mismatch (or report "no mismatch found in the range checked" if none). This is a
pure read/diagnostic task -- not modifying lift.h or the plan itself.

## Current status (session 6 continued, 2026-07-06) — user reports {3-1}*{1-1} STILL zero
## after fresh rebuild, PLUS new more-basic failures ({2-1} as alpha, both bs=1 and bs=2).
## Traced root cause in detail; found a SEPARATE, more serious issue in the process.

**Re-checked lift.h's tau handling carefully — found no logic bug in the new correction-term
code.** Added `TRACE_LIFT` (env-gated, `lift.h` inside `lift_one_step`) printing c_self,
inj_row contents, im_prev/im_ck/im_gk/correction/combined for every cogenerator processed.
Ran with `TRACE_LIFT=10..60` at bs=1: **c_self was `t^0` in EVERY case observed** (levels 1
through 6, ncog up to 21) — including the one previously-known impure row (level 3, j=2,
x=776, `inj_row.size=2`) — so the "general c_self^{-1}" code path is never actually exercised
on a non-unit value in this data; matches the plan's expectation ("every case checked so far
has c_self = tau^0"). No evidence of a c_self-related bug.

**Traced the {2-1} (s=2,a=1, bs=1) zero to its exact root cause, step by step** (all via
TRACE_LIFT/TRACE_PROD_S+A/TRACE_QUT_POS, all still in the tree, env-gated):
- `cyc[2].at(1)` = pure `{cog1^t0}` (single term, not a sum) → M_beta[2].maps_to(rep) is
  zero because `phi_beta[2].find(897)` (897 = position of G_2's cogen 1) is empty.
- That's empty because in `lift_one_step` (k=1, cogen j=1, x=467): `im_prev =
  phi_beta[1].find(g_prev_pos=4)` is EMPTY.
- That's empty because `phi_beta[1].find(4)` depends (via `cofree_adjoint_row`, since
  position 4 is a non-cogenerator position inside G_1's cogen-0 summand) entirely on
  `lg[position 0]` (G_1's cogen 0), which was computed to be EMPTY at k=0.
- That's empty because at k=0, cogen j=0 (x=0) of G_1: `im_prev = phi_beta[0].find(1) =
  {p975^t0}` (NONEMPTY — position 975 = 974+1, i.e. beta's coaction-shifted copy one step
  into G_1 cogen-1's own summand; this itself is a legitimate, expected value, NOT a bug —
  974 is cog1 of G_1's position, so 975 is the first non-cogen offset within that summand).
  But `im_ck = R[1].qut.maps_to(im_prev)` is EMPTY.
- **Directly verified against the raw resolution data** (added `TRACE_QUT_POS`/`TRACE_QUT_LVL`
  env vars in `yoneda2.cpp`, right after "bockstein tables loaded"): `R[1].qut.find(975)` is
  genuinely, truly EMPTY in the on-disk resolution matrices — this is a raw fact about the
  data, not something computed/miscomputed by `lift_one_step`. So `lift_one_step`'s math is
  behaving exactly as its (correct, per the plan) formula dictates given real inputs — the
  zero is not fabricated by a coding error in my rewrite.

**NEW, more serious/orthogonal finding — cross-checked against the OLD `yoneda.cpp`
(previously reported to get `{3-1}*{1-1}`-type products right) using the SAME freshly
rebuilt resolution** (after this session's to_del-permanent fix + user's full mr_mot
rebuild): **`./yoneda 40 30 1 1` now CRASHES**: `INVARIANT BROKEN: add(t^2, t^1) --
aborting` (a `tauOper::add` call being asked to add two different-exponent tau monomials
at the same vector index — structurally should never happen for a well-formed sparse
vector; see `mot_steenrod.cpp` and the plan's Appendix note #1 on `tauOper::inverse`).
Meanwhile **`./yoneda 40 30` (the solver self-check, verifying `qut_i`/`inj_i` invert
correctly at every level) still reports "SOLVER CHECK: success"** on this same resolution.
So the resolution's basic per-level invariants still hold, but yoneda.cpp's own (more
elaborate, Gaussian-elimination-based) product algorithm now hits a code path it apparently
never hit before, on a product that used to work, on the SAME freshly-rebuilt resolution
that yoneda2.cpp is also using. This raises real doubt about whether the "to_del fix has no
downside" conclusion from earlier this session (verified only via exactness-violation counts
and RANK comparisons, never by exercising downstream product code) was actually complete —
or alternatively this could be a preexisting yoneda.cpp bug that just never got exercised
before. **Have not yet determined which.** This is orthogonal to my lift.h changes (I did not
touch yoneda.cpp, and yoneda.cpp's crash is in its own from-scratch cycle-correction code,
not anything shared with lift.h), but it means the ground truth I'd normally cross-check
{2-1}/{3-1} against is currently unusable/crashing on this resolution.

**Instrumentation added this pass (all env-gated, still in the tree, not yet cleaned up)**:
`TRACE_LIFT` (lift.h, prints per-cogenerator c_self/im_prev/im_ck/im_gk/inj_row/
correction/combined, budget via `TRACE_LIFT=<N>`), `TRACE_PROD_S`/`TRACE_PROD_A` (yoneda2.cpp,
dumps a cycle representative and M_beta application for a given (s,a)), `TRACE_QUT_POS`/
`TRACE_QUT_LVL` (yoneda2.cpp, dumps a raw qut row at a given level/position).

**Next step**: report both findings to the user (lift.h math checks out on the data it's
given; but the {2-1} zero traces back to a raw `qut_1.find(975)==empty` fact about the
resolution itself, AND yoneda.cpp itself now crashes on the same fresh resolution) and ask
how to proceed — e.g. whether to investigate the yoneda.cpp crash / to_del side effect next,
or independently verify whether `qut_1(position 975)` SHOULD be nonzero mathematically.

## Current status (session 6, 2026-07-06) — implementation done, compiles, NOT YET VERIFIED

Implemented the plan at `~/.claude/plans/zany-baking-steele.md` in full:

- `lift.h`'s `lift_one_step` rewritten: takes a new `matrix<ring> &inj_src` parameter
  (`INJ_{k+1}: X_{k+1} -> G_{k+1}`), builds `lg` via an explicit sequential loop over
  cogenerators `j=0..rank-1` (ascending position order, per `position_of_gens` already
  being ascending) instead of the old single non-sequential `matrix_mem::construct` call.
  For each cogenerator: reads `inj_src.find(x)`, separates the self term (`c_self` at
  `cog_pos`) from correction terms at other (earlier) positions, computes
  `phi_{k+1}(pos_i)` via `cofree_adjoint_row(lg, ...)` on the partially-built `lg`,
  reindexes both the main term and each correction term to `G_tgt` cogen-index space via
  the same `cr = G_tgt.find_index` rule, subtracts (`add` + `minus`, correction UNSCALED
  relative to the already-tau^{-n}-divided main term, per the plan's order-of-operations
  warning), then multiplies by `c_self^{-1}` if `c_self != tau^0`, and `lg.insert(cog_pos,
  combined)`. Final `phi.construct(...)` step unchanged (`lift.h:200-204`).
- `yoneda2.cpp`'s `lift_one_step(...)` call site (~line 245-255): added `R[k+1].inj` as
  the new `inj_src` argument, correctly distinct from the two other `inj`/`qut` args
  already passed (`R[k+bs].qut`, `R[k+bs+1].inj`).
- Stripped all now-superseded debug instrumentation exactly as the plan listed: in
  `yoneda2.cpp` removed TMPCHECKNEG/TMPCHECK0/TMPCHECK-exactness (right after "resolution
  loaded"), TRACE_CYC block, and TRACE_LIVE_DIV/TRACE_IMPURE_HIT/TRACEB inside the k-loop.
  (Left TMPCENSUS alone — not in the plan's explicit removal list.) In
  `hopf_algebroid/6.h` removed TRACE_RM/TRACE_FRANK/TRACE_P1/TRACE_P2 and
  TRACE_PIVOT_SCAN/TRACE_FRANK_SCAN, and made the `to_del`-skip (the `to_del` bugfix)
  unconditional/default by deleting the `EXPERIMENT_NO_DELETE` env-var gate entirely —
  `to_del` is now simply always empty, matching what the plan asked for ("make it the
  unconditional default ... rather than leaving it opt-in").
- **Build verified**: `cmake --build build --target yoneda2` and `--target mr_mot` both
  compile and link cleanly (only pre-existing unrelated sign-compare warnings elsewhere).
  Editor/clangd diagnostics seen during editing (missing `matrices.h`, undeclared
  `matrix_index`, etc.) are standalone-file include-path noise from clangd, not real
  build errors — confirmed by the actual cmake build succeeding.

**NOT done this pass** (deliberately, per user's token-budget instruction to just get it
compiling and stop): rebuilding the motivic resolution fresh, clearing `40_yoneda2_*_phi*`
caches, and re-running the verification steps from the plan (exactness re-check,
`{3-1}*{1-1}`, `{6-9}*{4-6}`/`{4-6}*{6-9}`, regression checks on `{1-1}*{3-1}` and the
bs=6 case). **Next step next session: run the plan's "Verification" section end-to-end**
— rebuild `mr_mot`/`yoneda2`, regenerate resolution data fresh, clear stale phi caches
before testing (stale caches have caused false readings repeatedly this debugging effort),
then check both the originally-broken and previously-correct products.

# Debugging Notes: Yoneda Product Zero Bug (session 5 start, 2026-07-06)

Full history from earlier sessions is in `debug_notes_archive1.md` and
`debug_notes_archive2.md`. This file is the current, active summary — keep it short and
correct; append detail only as it becomes load-bearing for the next step.

## What we're chasing

Yoneda products in the C-motivic Adams E2 page that should be nonzero are computing as 0
via `yoneda2.cpp` (the "clean" chain-map implementation we're supposed to get working,
per the plan at `~/.claude/plans/tidy-yawning-dongarra.md`). Smallest known repro:
`./build/yoneda2 40 30 3 1 1 1` gives `{3-1}*{1-1} = 0` (should be nonzero; the mirror
order `{1-1}*{3-1}` correctly gives `t^0{4-1}`).

## Two real, confirmed bugs found and fixed so far

1. **`to_del` in `resolvor_modeled` (hopf_algebroid/6.h)** — deletes columns from `indj`
   (scratch matrix used to build `qut` via Gaussian elimination) for "complement
   positions not needed as generators at the next level," but does this BEFORE using
   `indj` to compute the CURRENT level's `qut` — corrupting `qut` even for OTHER rows
   that still needed the deleted column as a cross-reference. Fully traced with a
   concrete minimal example (md=10 build): `qut(pos 10)` was wrongly empty because its
   pivot row's needed cross-reference to position 48 got deleted by `to_del` before
   Gaussian elimination ran.
   - **Fix**: env-var-gated experiment (`EXPERIMENT_NO_DELETE=1`) that skips building
     `to_del` entirely (leaves it empty). Verified at BOTH md=10 and md=40: (a) zero
     exactness violations (`qut_i(inj_i(y))==0` for every y, every level checked) with
     the fix, vs. widespread violations (up to 96% of rows) without it; (b) ranks
     IDENTICAL between fixed and broken versions at every level — no structural/Ext
     change at all, purely a correctness fix with no apparent downside found so far.
   - **This fix alone did NOT fix the actual Yoneda products** — `{3-1}*{1-1}` etc. are
     still 0 even with `to_del` disabled. So there's at least one more issue.

2. Two theories investigated and RULED OUT with exhaustive (not just anecdotal) checks:
   - Gaussian-elimination `unify()`/`tauOper::inverse()` corruption (matrices_mem/4.h):
     `inverse()` doesn't check `invertible()`, structurally looked risky, but scanned
     ALL ~9051 pivot self-coefficients across 7 levels — zero non-invertible cases found.
   - `lift_one_step`'s own `tau^{-n}` division (lift.h): scanned every cogen at every
     level of the `{3-1}*{1-1}` computation for a case where `inv_tau!=0` AND `im_prev`
     is simultaneously nonzero (the only situation where this division could matter) —
     zero matches found anywhere. The division is never numerically exercised on a
     nonzero value in this computation.
   - `tauVal`/`algebroid2vector` (mot_steenrod.cpp) — the per-monomial tau-bookkeeping
     system (`tauVal(e) = sum floor(exponent_i/2)`, encoding the standard "doubling"
     trick for A_C/tau's associated-graded structure) — read and understood, looks like
     sound, well-established, pervasively-used infrastructure; no bug found, not pursued
     further without a more specific lead.

## Most promising NEW lead — impure inj rows starting at level 3

The plan's `lift_one_step` formula assumes `inj_k(e_j)` (the injection map row at a
cogenerator) is always a PURE singleton (`tau^0` at the cogenerator's own position, no
other terms) — session 3 verified this held for the `{6-9}`/`{4-6}` beta case (bs=6).
**Just re-checked this assumption for the CURRENT beta (`{1-1}`, bs=1) and it FAILS**
starting at level 3:
```
lvl=1: impure_inj_rows=0/6      (pure, as assumed)
lvl=2: impure_inj_rows=0/15     (pure, as assumed)
lvl=3: impure_inj_rows=1/18     -- j=2, x=776, cog_pos=1461, row={p827^t1, p1461^t0}
lvl=4: impure_inj_rows=5/21     -- 5 rows have an extra tau^1 term beyond the self singleton
```
So for levels 3+, `inj_k(e_j)` is genuinely NOT a pure singleton — it has an EXTRA term
(always at `tau^1` in the samples seen) beyond the cogenerator's own position. This
directly means the plan's simplifying assumption (`lift.h`'s comment: "No subtraction of
other inj-row terms is needed") is WRONG for this case, and `lift_one_step`'s formula is
silently dropping the contribution from these extra terms. This looks like the most
likely remaining root cause — NOT yet confirmed as the actual fix, and NOT yet acted on.

**Important**: the session-2 attempt to add "correction terms" for exactly this situation
was buggy and was reverted earlier this session — but the REASON for reverting was that
it was unjustified/undocumented scope creep AND didn't fix anything at the time (before
the to_del fix). Now that to_del is fixed AND we have concrete evidence that impure rows
really exist, correction-term handling may need to be reintroduced — but carefully, and
NOT without the user's explicit sign-off (see constraints below).

## Constraints (from the user, updated 2026-07-06)

- User has SLIGHTLY relaxed the "don't touch yoneda.cpp" rule: may now read yoneda.cpp
  (both committed and uncommitted-history versions) SPECIFICALLY to compare its
  RESOLUTION DATA against yoneda2's resolution — both versions of yoneda.cpp reportedly
  get `{3-1}*{1-1}`-type products RIGHT. Explicit instructions:
  - Compare the RESOLUTION (data/structure), not the code/algorithm — "don't take
    inspiration from the code."
  - Apply the to_del fix when generating/comparing yoneda.cpp's resolution if relevant.
  - **DO NOT adopt any strategy changes without asking the user first** — this is a
    read-only investigation task right now, not a license to start modifying lift.h/
    yoneda2.cpp based on what's found.
- Keep using the debug-notes-update hook (fires every 5 tool calls) — confirmed working
  well, keep going without being reminded.
- Still must not read git history of yoneda.cpp without separately asking (the relaxation
  was specifically about reading yoneda.cpp's current file and comparing resolution data;
  re-confirm before diffing against old commits if that becomes relevant).

## Instrumentation currently in the tree (all still present, all env-var-gated/harmless)

- `yoneda2.cpp`: TMPCHECK0, TMPCHECK, TMPCHECKNEG (exactness/dimension checks, right
  after "resolution loaded"), TRACEB (per-cogen trace in the k-loop, bs==1||bs==3),
  TRACE_CYC (env `TRACE_CYC`, cycle/bidegree/inj-purity inspection), TRACE_LIVE_DIV (env
  `TRACE_LIVE_DIV`, scans for numerically-live tau divisions).
- `hopf_algebroid/6.h` (`resolvor_modeled`): TRACE_RM/TRACE_FRANK/TRACE_P1/TRACE_P2 (env-
  gated position tracing), TRACE_PIVOT_SCAN/TRACE_FRANK_SCAN (pivot self-coefficient
  scan), `EXPERIMENT_NO_DELETE` (the to_del fix toggle — currently the way to get the
  "fixed" resolution; not yet made permanent/default).
- Scratchpad tools (session scratchpad dir): `check_ex_exactness.cpp` (classical F2
  resolution loader/checker, supports text dump mode), `dump_mot_mod_tau.cpp` (motivic
  mod-tau dumper), `check_mot_exactness_small.cpp` (small motivic exactness checker with
  full per-row detail).
- `ex_backup_fresh/`, `mot_backup_baseline/`, `exp_nodel/` — backup/comparison
  directories from earlier in the session; safe to keep around for now.

## Resolution-comparison request — resolved quickly, redirected

Checked `yoneda.cpp`'s file loading (`load_F`/`load_maps`, lines 65-75): it reads the
EXACT SAME `<pre>mot_gens{i}`/`<pre>mot_maps{i}` files as yoneda2.cpp — same filenames,
same format, no separate resolution-building step inside yoneda.cpp. So there is no
"yoneda.cpp resolution" distinct from "yoneda2.cpp resolution" to diff — both consume
identical on-disk data. This means the difference must be in the CHAIN-MAP ALGORITHM,
not the resolution. Reported this to user; user asked me to instead keep chasing the
`inj_k(e_j)` impurity finding directly (without touching yoneda.cpp), so redirected.

## Chasing the impure inj_k(e_j) finding

Added `TRACE_IMPURE_HIT` (yoneda2.cpp, in the k-loop): scans every cogen at every level
for a case where `im_ck` (computed from a genuinely nonempty `im_prev`) has a term
landing on an X-index whose `inj_{k+bs+1}` row is IMPURE (extra terms beyond the self
singleton) — i.e. checking whether the plan's "no correction terms needed" assumption is
actually exercised on a nonzero value anywhere in the `{3-1}*{1-1}` computation. This is
the direct test of whether the impure-inj-row finding (levels 3-4 have impure rows, see
above) is actually ON THE CRITICAL PATH for this product, or just structurally present
but never touched (like the earlier to_del-adjacent findings that turned out inert).

**Ran TRACE_IMPURE_HIT — but had mis-targeted the check at first, then corrected it.**
Initial run found one hit (`k=1 j=14 x=2355`, `im_ck` landing on X-index 1776 whose
`R[3].inj` row is impure) — but on reflection this was checking the WRONG inj: I'd
checked `R[k+bs+1].inj` (the TARGET-level injection used via `im_gk =
inj_tgt.maps_to(im_ck)`), which handles impure rows CORRECTLY regardless of purity
(`maps_to` sums the full row, no special casing needed there).

**Re-derived where impurity actually matters** (worth recording carefully, since I got
this wrong once already): the plan's claim is about `inj_{k+1}`'s OWN row (`inj_k(e_j)`
in the original notation) — the SOURCE-level injection, used only IMPLICITLY via
`cofree_adjoint_row`, never passed explicitly to `lift_one_step`. The actual chain
condition is `phi_k(inj_k(e_j)) = im_gk/tau^n`. Since `phi_k` is linear and
`inj_k(e_j) = cog_j + (correction terms at OTHER, earlier cogens)`, this expands to
`lg[cog_j] + phi_k(correction terms) = im_gk/tau^n`, i.e.
`lg[cog_j] = im_gk/tau^n - phi_k(correction terms)`. The CURRENT simplified formula
just sets `lg[cog_j] = im_gk/tau^n`, silently DROPPING the correction whenever
`inj_{k+1}`'s row for that cogen is impure. This is exactly the check I'd already done
via `TRACE_CYC`'s `impure_inj_rows` scan (R[lvl].inj, SAME level as the cogens) — found
impure rows starting at level 3.

**Confirmed on our exact example**: level 3, cogen j=2 (x=776) has impure
`inj_3(776) = {p827^t1, p1461^t0}` (1461=self, 827=correction term, presumably belonging
to an earlier cogen of G_3 since 827<1461). Re-checked `TRACEB` for k=2, j=2: `im_prev`
is EMPTY, so the CURRENT formula gives `lg[cog_2 of G_3] = im_gk/tau^n = 0`. But per the
correct formula, `lg[cog_2] = 0 - phi_3(position 827) = -phi_3(position 827)` — which is
NONZERO if `phi_beta[3]` at position 827 (an earlier cogen) is itself nonzero. **The
current code cannot even compute this correctly as structured** — `lift_one_step`
builds ALL of `lg` via a single non-sequential `construct()` call using only `phi_prev`
(level k, already complete); it has no mechanism to reference `phi_beta[k+1]`'s OWN
(earlier-cogen) values while still building `phi_beta[k+1]` itself. Doing so requires
SEQUENTIAL processing (cogen by cogen, using already-computed earlier lg entries for
corrections) — which is structurally what the reverted session-2 "correction terms"
code attempted, though it was reverted for good reasons at the time (undocumented scope
creep, didn't fix anything before the to_del fix existed).

**Conclusion so far**: this looks like a real, well-justified, mathematically necessary
gap in `lift_one_step`'s simplified formula — not a hunch, a derived consequence of
`inj_{k+1}` being genuinely impure at levels 3+. Re-introducing SOME form of sequential
correction-term handling (done carefully, unlike the buggy session-2 attempt, and now
justified by concrete evidence rather than a vague worry) looks like the most promising
concrete next action. **This would be a real strategy change to lift.h/yoneda2.cpp — have
NOT implemented anything, per the user's explicit instruction to ask first.**

**Next step**: report this finding clearly to the user and ask permission before
touching lift.h. If given the go-ahead, the fix would need: (a) process cogens of
G_{k+1} in position order within `lift_one_step`, (b) for each cogen j, after computing
the main `im_gk/tau^n` term, subtract `phi_{k+1}(correction term positions)` using
ALREADY-SET `lg`/`phi` values for earlier (lower-position) cogens via
`cofree_adjoint_row`, (c) verify carefully against BOTH this failing case AND the
previously-passing cases (e.g. `{1-1}*{3-1}`, and the bs=6 case from session 3 where
`inj` rows were pure and corrections should be a no-op) to avoid repeating session-2's
mistake of breaking things silently.
