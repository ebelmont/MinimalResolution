# Debugging Notes: Yoneda Product Zero Bug (session 4 start, 2026-07-05)

Full history of dead ends and detailed traces from earlier sessions is in
`debug_notes_archive1.md`. This file is the current, active summary — keep it short and
correct rather than exhaustive; append detail only as it becomes load-bearing for the
next step.

## What we're chasing

Yoneda products in the C-motivic Adams E2 page that should be nonzero are computing as 0.
Concrete failing examples (md=40 len=30, via `./build/yoneda2 40 30 <a_s> <a_a> <b_s>
<b_b>`):
- `{6-9} * {4-6}` and `{4-6} * {6-9}` (both orders) — both give 0, should be nonzero.
- Smaller repro: `{3-1} * {1-1}` gives 0 (should be nonzero); the mirror order
  `{1-1} * {3-1}` gives the correct nonzero answer `t^0{4-1}`.

## What's ruled out (high confidence, with evidence)

- **yoneda2.cpp / lift.h (the chain-map construction we're supposed to fix per the plan
  at `~/.claude/plans/tidy-yawning-dongarra.md`) is NOT the bug.** Reverted an earlier
  session's unjustified "correction terms" deviation from the plan back to the clean
  singleton-based formula. Verified directly: (a) any valid qut-preimage gives the same
  phi value (exactness-based invariant, proven by hand and confirmed empirically); (b)
  every `inj_row` checked is a pure `tau^0` singleton, matching the plan's simplifying
  assumption exactly, no hidden correction terms needed.
- **The classical (F2, no tau) resolution built by `mr_ex` is exact.** Directly checked
  `qut_i(inj_i(y)) == 0` for every y at every level 0-7 on a **freshly rebuilt**
  `mr_ex 40 32` resolution (the old data had a corrupted 0-byte `40_gens` file — always
  regenerate fresh before trusting these files). Zero violations at every level. This
  resolution-construction machinery (Gaussian elimination, curtis tables) works
  correctly, at least at F2.
- **"tau=1" collapse doesn't explain it.** Collapsing tau-exponent distinctions (XOR by
  position only, ignoring exponent) on the motivic exactness violations gives the SAME
  nonzero support as the tau-respecting computation. Not a tau-exponent-mismatch-prevents-
  cancellation issue.
- **"mod tau" (tau->0) does NOT make the motivic inj/qut match the classical ones either**
  — see below, this was the most recent hypothesis and it's now falsified too.

## What's confirmed broken

**The motivic (tauPoly) resolution violates qut_i(inj_i(y))=0 pervasively** — 78-96% of
rows at levels 1+ (level 0 is trivially clean, rank 1). This was checked directly on the
loaded `mot_maps<i>` data used by yoneda2.cpp.

**This points at the "modeled" F2->F2[tau] lifting stage specifically**, not the
resolution-construction machinery in general:
- `mot_main.cpp` calls `pre_resolution_modeled` / `resolvor_modeled`
  (hopf_algebroid/6.h), NOT the plain `resolvor`/`pre_resolution_tab`
  (hopf_algebroid/5.h) that mr_ex uses and that's proven exact.
- `resolvor_modeled` reuses a classical (F2) "model" resolution (built by `mr_ex`,
  loaded via `gens[]`/curtis tables) and reconstructs the tauPoly inj/qut from it via
  `embed2cofree_modeled` (hopf_algebroid/7.h:176-198) and a special `cycle_matrix`
  overload (others.h:266-282 / matrices/9.h:186-200) that takes a `transformer` +
  `Map` argument.
- In `mot_main.cpp`, `transformer = [&MOP](vectors<matrix_index,F2> const&v){ return
  MOP.lift(v); }`, and `MOP.lift()` (mot_steenrod.cpp:293-304) maps every nonzero F2
  coefficient to exactly `tau^0`, always.
- User's insight (still unconfirmed in detail): the lift is NOT simply "apply the trivial
  F2->F2[tau] inclusion" — `resolvor_modeled` genuinely recomputes inj/qut using data from
  the F2 model, and this recomputation is the suspected site of a **misunderstanding**
  (self-consistent-but-undocumented internal representation), not necessarily a "bug" in
  the naive sense. User remains confident the resolution is NOT simply wrong 80% of the
  time — something about how we're reading/interpreting it is off.

**Latest test (mod-tau comparison, just completed, and it's informative but not
conclusive)**: hypothesis was that reducing the motivic inj/qut "mod tau" (keep only
tau^0 entries, drop tau^n n!=0) should reproduce the classical F2 inj/qut exactly. Built
this comparison (scratchpad/check_ex_exactness.cpp dump mode + scratchpad/
dump_mot_mod_tau.cpp, diffed as plain text). Result: they do NOT match, in two different
ways:
1. `qut`: 30-45% of rows have entries in the classical version "missing" from the
   motivic mod-tau reduction (i.e. the motivic entry exists but has nonzero tau power).
   Scale is large enough that this might just reflect genuine/expected tau-grading,
   inconclusive on its own.
2. `inj`: a handful of rows per level (4-19) have EXTRA columns in the motivic mod-tau
   version that don't exist in the classical version at all — e.g. level 2 row 1009:
   classical `{2365,3370,3371}`, motivic-mod-tau `{2364,2365,3370,3371}`. This is
   structurally different from a re-grading — it's new support. Not yet investigated
   further (haven't looked at the raw, non-reduced tauPoly value at these positions).

## Data/tooling notes

- **Always rebuild `mr_ex`/`mr_mot` fresh before trusting `40_*` derived files** — found
  and fixed a corrupted 0-byte `40_gens` mid-session. Backup of a known-good fresh
  classical resolution is in `ex_backup_fresh/` (40_maps, 40_gens, 40_inex, 40_qutex,
  40_indjex, 40_extables, 40_gens_data_ctau, 40_ResTables_ctau,
  40_ctau_steenrod_coaction.data, 40_excacexcoamat).
- Scratchpad tools (in the session's scratchpad dir, not in the repo):
  `check_ex_exactness.cpp` (loads classical `<pre>maps`/`<pre>gens` combined-file format
  via `SteenrodInit(...,IorO=true)`; supports a dump-to-text mode), `dump_mot_mod_tau.cpp`
  (loads motivic per-level `mot_maps<i>` files, dumps mod-tau reduction to matching text
  format). Note: classical (`ex_exponents.h`) and motivic (`exponents.h`) headers define
  `maxVar`/`unpack` incompatibly and CANNOT be included in the same translation unit —
  keep these as separate small programs, compare via text dump + diff.
- Current yoneda2.cpp has temporary instrumentation still in place (TMPCHECK0, TMPCHECK,
  TMPCHECKNEG right after "resolution loaded"; TRACEB in the k-loop) — harmless (gated by
  specific `bs`/`k` values or always-run read-only checks) but should be cleaned up before
  considering this done.

## Constraints (from the user, still in force)

- Do not read `yoneda.cpp` or its git history without asking first — it's a previous
  incorrect attempt, deliberately out of scope.
- May read any file linked/included by yoneda2.cpp.
- Keep using the debug-notes-update hook (fires every 5 tool calls) — user confirmed
  this workflow is working well, keep it up without being reminded.

## FOUND: smallest concrete example, pinned down to one missing matrix entry

Built a small max_deg=10 resolution from scratch (`e2p 10; motTab 10; mr_ex 10 12;
mr_mot 10 10` — clean slate, no stale 10_* files existed). Classical resolution at md=10
is confirmed exact (0 violations, levels 0-2, via scratchpad/check_ex_exactness.cpp).
Motivic resolution at md=10 has violations starting at level 1 (tiny — inj.rank=30,
qut.rank=60 — small enough to read by hand). Wrote scratchpad/check_mot_exactness_small.cpp
to print full per-row detail for every violation.

**The minimal example — level 1, y=7** (X_1 index 7's image under inj_1):
- Motivic `inj(7) = {p28^t0, p29^t0, p46^t0}` — MATCHES classical `inj(7) = {28,29,46}`
  exactly (confirmed via scratchpad/cl10_inj1.txt row 8). inj is NOT the problem here.
- Classical `qut(28) = {12, 20}`, `qut(29) = {12}`, `qut(46) = {20}` (scratchpad/
  cl10_qut1.txt rows 29,30,47). Classical `qut(inj(7))` = {12,20}+{12}+{20} = 0 exactly
  (12 cancels between rows 28&29, 20 cancels between rows 28&46) — CORRECT, exact.
- Motivic (RAW, not mod-tau-reduced) `qut.find(28) = {x20^t0}` ONLY — **completely
  missing the entry to x12** that the classical structure has! `qut.find(29) = {x12^t0}`
  (matches classical), `qut.find(46) = {x20^t0}` (matches classical). So motivic
  `qut(inj(7))` = {x20^t0} + {x12^t0} + {x20^t0} = {x12^t0} (x20 cancels, but x12 has no
  partner to cancel it since row 28 is missing its x12 component) — NONZERO, the
  violation.

**Root cause pinned to one specific defect**: motivic `qut` row 28 is missing a
component (to X-index 12) that both the classical model AND exactness require. This is
not a tau-grading/exponent issue at all — the entry is COMPLETELY ABSENT, not present
with a wrong exponent. This is the single smallest reproducible instance of the bug.

**User has zero confidence in the classical-vs-motivic comparison approach** — dropping
that angle entirely. Refocused purely on finding the smallest `qut(inj(y)) != 0` failure
and tracing the ACTUAL motivic construction code from the beginning for that one example,
without reference to the classical model.

**Even simpler example found, y=12 (level 1)**: `inj(12) = {p10^t0, p48^t0}` (just 2
terms). `qut.find(10) = {}` (completely empty), `qut.find(48) = {x21^t0}`. Sum = `{x21^t0}`
— nonzero, no cancellation partner at all (not even a "missing one of two" case like
y=7 — position 10 contributes NOTHING). This is the new primary trace target: simplest
possible violation, 2 terms, one of which is entirely zero.

(y=18 is equally simple: `inj(18)={p15,p51}`, `qut(15)={}`, `qut(51)={x23}` — same
shape, kept as a backup/cross-check example.)

**Instrumented resolvor_modeled (hopf_algebroid/6.h)**: added TRACE_RM/TRACE_FRANK/
TRACE_P1/TRACE_P2 env-var-gated debug prints from the top of the function, tracing inj,
inj_ind, indj (before/after gaussian), quot_inds, to_del, and final qut — all at G-1
positions 10 and 48. Gated by `F.rank() == TRACE_FRANK` so it only fires for the ONE
relevant call (level 1, F.rank()=60), not every level.

**First attempt crashed + had a real bug in the trace itself, now fixed**: initially
looked up `indj->find(10)`/`indj->find(48)` directly — WRONG. `indj`'s rows are indexed
by X-TAG (0..X.rank()-1), not by G-position, unlike `inj`/`qut` which ARE G-position-
indexed. Also the trace wasn't gated by level, so it fired (and crashed via out-of-bounds
access) at level 0 where these G-1-specific positions don't even exist. Fixed both: gate
by `F.rank()==TRACE_FRANK`, and look up indj via the TAG i such that `inj_ind[i]==P1/P2`
(if such a tag exists at all — a position might be a "complement" position, with no
corresponding indj row).

## ROOT CAUSE FOUND — fully traced, not just inferred

Ran `TRACE_RM=1 TRACE_FRANK=60 TRACE_P1=10 TRACE_P2=48 ./build/mr_mot 10 10`. Full trace
for level 1 (F.rank()=60, the G_1 that contains our target positions):

```
inj.find(11) contains pos 48
inj.find(12) contains pos 10          <- confirms inj(12) = {10, 48}, our target row
TAG for P1(10) = 12   TAG for P2(48) = -1
indj.find(tag=12 ->pos10) [BEFORE gaussian]: 10^t0 48^t0     <- has BOTH 10 and 48!
P1=10 is PIVOT
P2=48 is COMPLEMENT (quot idx 21)
to_del contains P1? no   to_del contains P2? YES   (to_del.size()=24)
indj.find(tag=12 ->pos10) [AFTER del_and_gaussian]: 10^t0     <- the "48" entry is GONE
qut.find(10): (empty)
qut.find(48): x21^t0
```

**The causal chain, fully explained**:
1. Position 10 in G_1 is a PIVOT (tag 12 — i.e. it's X_1's own basis element 12's
   "leading term" position). Position 48 is a COMPLEMENT position (destined to become
   X_2's generator 21).
2. `indj`'s raw row for tag 12 (BEFORE any deletion/elimination) correctly contains BOTH
   `10^t0` (its own pivot term) AND `48^t0` — this second term is exactly what SHOULD let
   `qut(10)` correctly end up as `{x21}` after `make_quotient`'s pivot formula
   (`-row(i) + singleton(pivot)`, which needs the `48` term to survive as the nonzero
   "correction" after canceling the pivot term against itself).
3. **But position 48 is in `to_del`** — `to_del` starts as ALL complement positions
   (`quot_inds.first`) and then REMOVES (`erase`s) only those needed as `nex_gens`
   (gens[2], the classically-known list of X_2's "real" generators). Position 48
   (-> X_2 generator 21) apparently isn't in `nex_gens`, so it stays marked for deletion.
4. `del_and_gaussian(gs, to_del)` calls `del_cols(to_del)` FIRST — this deletes column 48
   from EVERY row of `indj`, INCLUDING tag 12's pivot row, wiping out the `48^t0` term
   that was needed.
5. `make_quotient` then builds `qut(10)` from this now-incomplete row — with the `48`
   term gone, the pivot formula's `-row(i) + singleton(pivot)` cancels to exactly zero.
   `qut(48)` itself is unaffected (complement positions get their identity assignment
   directly, not from `indj`'s rows) — hence `qut(48) = {x21}` survives untouched, with
   no partner left to cancel it in `qut(inj(12)) = qut(10) + qut(48)`.

**Root cause, one sentence**: `resolvor_modeled` (hopf_algebroid/6.h) deletes columns
for "complement positions not needed as generators at the next level" (`to_del`) from
`indj` BEFORE using `indj` to compute the CURRENT level's `qut` — but those columns can
still be load-bearing for correctly computing OTHER rows' `qut` values at the CURRENT
level (like pivot row 10's dependence on column 48), even if position 48 itself isn't
needed as a generator going forward. The deletion is happening at the wrong point in the
pipeline relative to what it's allowed to assume is safe to discard.

**Not yet done**: haven't identified what `to_del`/this whole pruning step is FOR (i.e.
why `resolvor_modeled` wants to delete "unneeded future generators" from indj at all —
some kind of intentional model-driven pruning, possibly to keep resolution size closer to
the classical model's, but as constructed it silently breaks exactness). Haven't proposed
or attempted a fix. This seems like exactly the kind of architectural finding that
warrants stopping to report back before touching resolution-construction code further.

## User's follow-up questions (2026-07-05, answered/in progress)

- **inj vs indj**: `inj` is the real, saved `INJ_i: X_i->G_i` (used everywhere
  downstream). `indj` is throwaway scratch, rebuilt from curtis-table/model data purely
  to drive the Gaussian elimination that produces `qut` — never saved, not equal to
  `inj`'s values.
- **Bidegree of the class in question**: added `F.degree(P1)`/`F.degree(P2)` to the
  trace. Both G_1 positions 10 and 48 are at `(deg=7, weight=3)` (internal degree 7,
  motivic weight 3). Position 48 becomes X_2's generator 21 — internal grading is
  preserved by the differential, so this is an Ext class at filtration s=2, internal
  degree t=7, weight 3 (stem t-s=5).
- **Experiment DONE — very clean, encouraging result**: added an env-var-gated override
  (`EXPERIMENT_NO_DELETE=1`) in hopf_algebroid/6.h that skips building `to_del` entirely
  (leaves it empty, so `del_and_gaussian` doesn't delete any columns). Reran the whole
  md=10 resolution from scratch in a scratch subdirectory (`exp_nodel/`, copying the
  shared inputs — ex2poly_index/mot_deltas/poly_exponents/gens_data_ctau/ResTables_ctau/
  extables — then running `EXPERIMENT_NO_DELETE=1 ../build/mr_mot 10 10` from inside it
  so outputs don't clobber the known-bad baseline).
  - **Exactness check: ZERO violations across levels 0-2** (previously the baseline had
    violations starting at level 1). Specifically `qut(inj(12))` is now empty (correctly
    cancels), exactly as predicted.
  - **Ranks are IDENTICAL to the baseline at every level 0-4** (inj.rank/qut.rank exactly
    match between the deletion and no-deletion runs). This means `to_del` has ZERO effect
    on which positions become generators / the resolution's size — that's governed
    entirely by `nex_gens`/`gens[i+1]` at the NEXT level's `embed2cofree_modeled` call,
    independent of `to_del`. `to_del`'s ONLY observed effect is corrupting `qut`'s VALUES
    by deleting columns from `indj` before Gaussian elimination, for no apparent benefit.
  - **Answer to "what happens to Ext if un-deleted"**: nothing structural — no spurious
    new generators appear, no existing generators disappear, same ranks throughout. The
    only change is that `qut` becomes exact (correct) instead of silently broken. This
    looks like a very low-risk fix: the `to_del` mechanism appears to be unnecessary
    entirely, at least for this md=10 test range.

## Scaling up the experiment to md=40 (per user request)

Testing at the ORIGINAL failing scale (md=40) to check (1) does removing `to_del` ever
change Ext (ranks/generators), (2) does it fix the qut*inj issue broadly, (3) do the
originally-broken products ({3-1}*{1-1}, {6-9}*{4-6}, {4-6}*{6-9}) become correctly
nonzero.

- Backed up baseline (with-deletion) 40_mot_gens*/40_mot_maps*/40_mot_res to
  `mot_backup_baseline/` before touching anything.
- Deleted the old 40_mot_gens*/40_mot_maps*/40_mot_res*, reran
  `EXPERIMENT_NO_DELETE=1 ./build/mr_mot 40 30` fresh — completed successfully.
- Classical (40_maps/40_gens/etc, from mr_ex) is untouched by any of this — `to_del`
  only exists in the motivic `resolvor_modeled` path.

**Results at md=40, levels 0-19 — very strong confirmation**:
- **Ranks IDENTICAL between baseline (with `to_del`) and no-deletion, at every one of the
  20 levels checked** (e.g. level 7: inj=1275/qut=2360 both ways; level 15: inj=197/
  qut=355 both ways). No structural change to the resolution/Ext generators at all,
  across the full range that matters for our failing examples.
- **ZERO exactness violations across all 20 levels** with `to_del` disabled (checked via
  scratchpad/check_mot_exactness_small, full qut_i(inj_i(y))==0 check for every row at
  every level 0-19). Previously (baseline) this had violations starting at level 1,
  78-96% of rows at several levels.

This is a very clean result: at both the small (md=10) and full (md=40) scale, removing
`to_del` entirely (a) never changes Ext (same ranks/generators throughout), (b) fully
fixes the qut∘inj exactness property everywhere checked.

**SURPRISING NEGATIVE RESULT — products are STILL zero despite exactness now being
perfect**: rebuilt yoneda2, cleared stale `40_yoneda2_*_phi*` caches (confirmed fresh
regeneration, not stale-cache artifacts), reran against the no-deletion resolution:
```
{3-1} * {1-1}  =  0        (still wrong, should be nonzero)
{1-1} * {3-1}  =  t^0{4-1} (control case, still correctly nonzero, unchanged)
{6-9} * {4-6}  =  0        (still wrong)
{4-6} * {6-9}  =  0        (still wrong)
```
So: `to_del` removal fixed qut∘inj=0 exactness completely (verified: 0 violations, all
20 levels, md=40) AND changed nothing structurally (identical ranks) — but the actual
broken Yoneda products are COMPLETELY UNCHANGED, still exactly 0. This means either (a)
the exactness bug we found and fixed, while real, is NOT the (or not the only) cause of
the wrong-zero products — there's a separate issue — or (b) fixing exactness at the
qut/inj level doesn't automatically propagate to fixing phi_beta's construction for some
other reason not yet understood.

Recall from earlier in this session (see debug_notes_archive1.md): phi_beta (the chain
map built by lift.h/yoneda2.cpp) was found to collapse entirely to zero by level 2 for
the {6-9} case, and this collapse was independently judged "mathematically legitimate
given the qut/inj data available at the time" — but that judgment was made using the
OLD (exactness-violating) data. Worth re ekhecking: does the SAME collapse-to-zero
pattern still happen with the FIXED (exact) qut/inj data? If yes, that confirms the
exactness bug wasn't the (sole) cause and there's a distinct issue in phi_beta's
construction or elsewhere. If the collapse pattern changed (e.g. phi_beta is less sparse
now) but the FINAL answer is still 0, that would point to a bug further downstream
(build_M, find_cycle, or the cycle representative itself).

## Continuing to debug {3-1}*{1-1} against the FIXED resolution (per user: "maybe there
were multiple issues")

**The to_del fix DID have a real, substantial effect on phi_beta — good sign, real
progress**: re-ran the level-by-level nonzero census (TMPCENSUS, re-added to yoneda2.cpp
before build_M) for bs=1 against the fixed (no-to_del) resolution:
```
                    OLD (broken qut/inj)     NEW (fixed qut/inj)
level 0: cogens 1/1     ALLPOS 585/1056   -> cogens 1/1   ALLPOS 1056/1056 (ALL nonzero!)
level 1: cogens ?/6     ALLPOS 416/3414   -> cogens 4/6   ALLPOS 1683/3414 (4x denser)
level 2: cogens 0/15    ALLPOS 0/4633     -> cogens 4/15  ALLPOS 657/4633  (was TOTAL zero, now has content!)
level 3: cogens 0/18    ALLPOS 0/4105     -> cogens 0/18  ALLPOS 0/4105    (still total collapse)
level 4-6: still 0                       -> still 0
```
So the fix pushed the "total collapse to zero" point from level 2 to level 3 — genuine,
substantial improvement, but NOT a full fix. Something else causes a complete collapse
between level 2 and level 3. Since alpha={3-1} needs phi_beta[3], this remaining bug is
exactly what's still causing the product to be wrong.

**Traced k=2->3 transition**: of the 18 cogens of G_3, only ONE (j=13, x=2243,
g_prev_pos=3575) has nonzero im_prev = `{p3040^t0}` (phi_beta[2] at position 3575). But
`im_ck = R[3].qut.maps_to({p3040^t0})` is EMPTY — i.e. `qut_3(3040) = 0`. All other 17
cogens have im_prev genuinely empty (unremarkable, given phi_beta[2] is only ~14% dense:
657/4633).

**Important: this can't be another to_del-type bug** — `EXPERIMENT_NO_DELETE=1` makes
`to_del` unconditionally EMPTY at every level, including level 3's own resolvor_modeled
call. We already independently verified qut_k(inj_k(y))=0 holds for ALL y at ALL 20
levels checked (0 violations) — so `qut_3(3040)=0` is very likely a genuinely CORRECT
value now (position 3040 legitimately in ker(qut_3)=image(inj_3)), not a data bug of the
same kind we just fixed. Getting only 1-of-18 canonical preimages to land in phi_beta[2]'s
~14%-dense nonzero support isn't statistically shocking on its own (expected ~2.5 hits;
got 1, and that 1 happened to land in qut_3's kernel).

**Reassessment**: the remaining collapse at level 3 may not indicate a SECOND distinct
"data corruption" bug like `to_del` — it may just mean phi_beta[1]/phi_beta[2] still
aren't as rich as they mathematically should be, for a DIFFERENT reason (not an inj/qut
data bug, since that's now verified exact). Candidates: (a) `lift_first_step`/beta_rep's
seeding at level 0 (even though level 0 is now 100% dense, maybe the SPECIFIC values are
still not quite right), (b) something in how `recover_inv_ind`/`lift_one_step` choose
which preimage to use — even though ANY valid preimage should give the same answer
in principle (proven earlier, and empirically confirmed once already this session), worth
re-verifying that invariant still holds now that the underlying data has changed, (c) a
genuinely separate bug elsewhere not yet identified.

**User's response**: surprised the F2-linear (classical) qut∘inj was actually exact
(expected it might be broken too). No concrete lead other than a strong hunch this is
tau-related, specifically thinks the remaining failures correspond to tau-TORSION
elements in Ext, but expects this to be hard to see at the low (raw matrix) level —
asked me to keep tracing carefully.

**Traced position 3040 in G_3 directly (TRACE_RM, TRACE_FRANK=4105, TRACE_P1=P2=3040)**:
- 3040 is a PIVOT, tag=1653 (X_3 element 1653's own leading term). Bidegree (deg=36,
  weight=18) — note deg=2*weight exactly, i.e. "chow degree" (deg-2*weight) is 0, which
  is the generic/expected relation for non-torsion motivic classes (a weak heuristic, not
  proof either way).
- `indj` row for tag=1653, BEFORE gaussian: `{3040^t0, 4070^t0}` — two terms: its own
  pivot (3040) plus another position, 4070.
- **`to_del` is completely EMPTY this run** (`to_del.size()=0`, confirmed —
  EXPERIMENT_NO_DELETE=1 is fully in effect). Yet AFTER `del_and_gaussian`, the row is
  JUST `{3040^t0}` — the `4070` term is GONE, even though nothing was explicitly
  deleted!
- This means the `4070` term was eliminated by ORDINARY Gaussian elimination itself
  (`gaussian(gs)`, canceling row 1653 against position 4070's OWN pivot row, since 4070
  must also be a pivot for some other tag) — NOT by any deletion mechanism. This is
  normal/expected elimination behavior IF the ring is well-behaved (a field) — but
  F2[tau] is NOT a field, and elimination (subtracting a multiple of one pivot row from
  another) is only guaranteed to preserve the correct quotient structure when the ring
  behaves like a field throughout. If pivot 4070 is itself connected to a tau-torsion
  relation, using it to eliminate another row's entry could be exactly the kind of
  operation that's invalid in this non-field setting — matching the user's tau-torsion
  hunch closely.
- Also notable: `sort_deg(gs, X.base_module.degree)` (hopf_algebroid/6.h:25, right
  before `del_and_gaussian`) explicitly SORTS the elimination order by degree — meaning
  the pivot elimination sequence is NOT arbitrary, it's degree-ordered. If this ordering
  interacts badly with torsion relations (e.g. eliminating using a "later"/torsion-linked
  pivot before a class's own relation is fully resolved), that's a very plausible
  mechanism for a subtle, hard-to-spot bug of exactly the kind being hunted.

**Traced position 4070 too — reconsidering, this may NOT be a bug after all**: 4070 is
ALSO a pivot (tag=2182), and its OWN raw indj row (before AND after gaussian) is just
`{4070^t0}` — a pure self-term, nothing else. `qut.find(4070)` is empty too. 26 different
X_3 rows (649, 652, 655, ... 2246) reference position 4070 in their inj expansion —
it's a heavily-shared, ordinary-looking pivot.

Given `make_quotient`'s pivot formula is `-row(i) + singleton(pivot)`, a pivot whose row
reduces to PURELY its own self-term legitimately gives qut=0 — this is completely
expected/normal for "dependent" pivot positions in general (most pivots SHOULD map to
zero under qut; only combinations reaching the true complement/generator space survive).
Canceling row 3040's "4070" term using 4070's own (valid, self-contained) pivot row looks
like textbook-correct Gaussian elimination, not obviously wrong.

**Reassessment**: qut(3040)=0 may simply be the mathematically CORRECT value — position
3040 (a pivot, not a generator) isn't guaranteed to be nonzero under qut at all. This
specific row may not be a bug; the real deficiency may be upstream, in phi_beta[1]/
phi_beta[2] not being AS RICH as they should be for some other reason (still only 14%
dense at level 2) — i.e. this specific (im_prev) draw for cogen j=13 of G_3 legitimately
lands on a pivot that's legitimately zero, and the actual missing richness is elsewhere
in phi_beta[2], not in this qut computation.

**tau_table_entry has a `diff_length` field (tao_bockstein.h:40)** — exactly the
Bockstein-differential-length concept for identifying torsion classes (a class supporting
or receiving a nontrivial tau-Bockstein differential is torsion; presumably some sentinel
value means "permanent"/non-torsion). This is the right tool to directly check whether
the ACTUAL Ext classes we care about ({3-1}, {1-1}, or whatever cyc[2] index corresponds
to G_2's cogen j=13) are torsion, rather than inferring from bidegree heuristics.

## MAJOR FINDING: direct link to a real tau-Bockstein differential

Built `tauBoc` and ran it fresh (`./build/tauBoc 40 30`, consistent with the to_del-fixed
resolution) to get `40_tau_bockstein.txt`, the actual tau-Bockstein differential table.
Checked our exact classes of interest:
```
{1-1}	:cycle       <- permanent (non-torsion)
{3-1}	:cycle       <- permanent (non-torsion)
t^1{4-1}	<-	{3-2}   <- {3-2} is TORSION, supports a Bockstein differential hitting {4-1}!
```
Neither `{1-1}` nor `{3-1}` themselves are torsion. BUT the target class of our WORKING
example — `{1-1}*{3-1} = t^0{4-1}` — is `{4-1}`, and `{4-1}` is EXACTLY the target of a
real tau-Bockstein differential from the torsion class `{3-2}` (`t^1{4-1} <- {3-2}`).
This is a genuine, concrete torsion phenomenon sitting exactly at the bidegree connecting
our working and broken examples (both `{1-1}*{3-1}` and `{3-1}*{1-1}` should land on
`{4-1}` at s=4). This strongly supports the user's hypothesis: reaching `{4-1}` via the
OPPOSITE product order (`{3-1}*{1-1}`) apparently requires correctly navigating this
torsion relationship (i.e. correctly constructing the chain map/resolution structure
around the {3-2}-torsion class), which the current construction fails to do — while the
WORKING order manages to avoid needing that navigation.

**Next step**: connect this directly to the resolution data — find which G_3/G_4
position(s) correspond to the torsion class `{3-2}` (cyc[3] index 2) and its target
`{4-1}` (cyc[4] index 1), and check whether OUR traced positions (3040, 4070, or the
G_2->G_3 chain we were following for cogen 13) relate to this torsion class's own
generator/pivot structure in the resolution. Also worth reading tao_boc.cpp around
line ~130-150 to fully understand what `diff_length`/the Bockstein table construction
implies about how such a class's resolution data should look, to know what "correct"
should look like for comparison.

**User's hypothesis to keep in mind while tracing (their words, speculative, "don't take
too seriously" but worth watching for)**: could the repeated `tau^{-k}` divisions in
lift_one_step's chain condition (multiplying by inv_tau's inverse at each step) be
implicitly introducing an EXTRA factor of tau into the computed product — i.e. instead of
computing the true class X, the code computes `tau^m * X` for some accumulated m from the
chain of divisions — and if X is itself tau-torsion (killed by tau^k for k<=m), then
`tau^m * X = 0` even though the true product X is nonzero?

**DIRECT CONFIRMATION that the resolution data itself correctly encodes this torsion
relationship** — added TRACE_CYC env-var-gated debug prints (yoneda2.cpp) to inspect
`cyc[3].at(2)`, `cyc[4].at(1)`, and the actual qut_3/inj_4 matrix entries at those exact
positions:
```
cyc[3].at(2) = cog2 (of G_3), cyc[4].at(1) = cog1 (of G_4)   -- both simple singletons
G_3 cogen2 position = 1461, G_4 cogen1 position = 757
qut_3.find(1461) = {x395^t1}      <- TAU EXPONENT 1, exactly matching the Bockstein
                                      table's "t^1{4-1} <- {3-2}"!
inj_4.find(395)  = {p757^t0}      <- pure self-term (G_4's own cogen1), tau^0
```
This is a genuinely exciting, concrete confirmation: the qut_3 matrix ITSELF directly
encodes the real tau^1 Bockstein relationship at exactly the position connecting these
two classes — the resolution's raw data is NOT wrong here, it correctly captures the
torsion structure. This is precisely the kind of place where `recover_inv_ind` (which
scans qut for singleton rows to find canonical preimages) would find `inv_tau=1` for
X_4 index 395 if it picks position 1461 as the preimage — triggering lift_one_step's
`tau^{-1}` division exactly at this torsion-connected step. This matches the user's
hypothesis structurally: dividing by tau^1 here is only valid/meaningful if what's being
divided is genuinely a multiple of tau^1 in a torsion-FREE sense; if the true chain-map
value at this point is itself connected to the {3-2} torsion class in a way that doesn't
admit clean division, this specific step is a strong candidate for silently producing 0
(or a wrong value) instead of the correct nonzero answer.

**Correction — traced the wrong level relationship at first**: realized `qut_3.find(1461)
={x395^t1}` (G_3 cogen2 -> X_4 index 395) is actually used when building phi_beta[4] from
phi_beta[3] (recover_inv_ind on R[3].qut, at loop iteration k=3) — but {3-1}*{1-1} only
needs phi_beta[3] (alpha's own filtration s=3 = maxlev needed, per fwd_product using
M_beta[s]). So that specific tau^1 relationship, while a real and exciting torsion
signature, is NOT on the critical path for THIS product. Real target: G_3's cogen INDEX 1
directly (= alpha's own class {3-1} = cyc[3].at(1) = cog1) — in the lift_one_step loop
building phi_beta[3] (k=2), the loop variable `j` IS the cogen index, so **j=1 IS
alpha's own class**.

**Checked j=1 directly (forced print, bypassing the "interesting rows" filter for k=2,
after clearing stale phi cache which was silently skipping the whole TRACEB block on
cache-hit runs — worth remembering: TRACEB is inside the `else` branch of `phi_exists`,
so it's silently skipped whenever cached phi files already exist)**:
```
TRACEB bs=1 k=2 j=1 x=430 g_prev_pos=5 inv_tau=1 im_prev: im_ck: im_gk:
```
**j=1 (alpha's own cogen) has inv_tau=1 — a NONZERO tau exponent right on alpha's own
chain condition!** This is the SECOND instance this session of a nonzero tau exponent
sitting exactly on/near the classes of interest (first was the {3-2}->{4-1} Bockstein
connection). However `im_prev` (phi_beta[2] at G_2 position 5) is already EMPTY here, so
the `tau^{-1}` division is numerically inert this time (0 divided by anything is still
0) — meaning the immediate zero at j=1 is not YET demonstrating the "torsion collapse"
mechanism directly; the real deficiency is upstream: WHY is phi_beta[2] empty at G_2
position 5?

**Traced backward through all levels — full chain now visible, but no numerically-visible
"tau collapse" yet**. G_2 position 5 belongs to cogen 0 of G_2. Forced full (unfiltered)
printing for k=1 and k=2 (note: TRACEB is inside the `phi_exists` cache-miss branch, so
it's silently SKIPPED on any run where `40_yoneda2_1_1_phi*` is already cached — must
clear that cache every time before re-running a trace, easy to forget). Full picture for
{3-1}*{1-1}, cogen-index-1 chain across all levels:
```
k=0 (building phi_beta[1]): j=1 x=1  g_prev_pos=2 inv_tau=0  im_prev=p976  im_ck=x467  im_gk=p897   -> NONZERO, survives
k=1 (building phi_beta[2]): j=0 x=0  g_prev_pos=1 inv_tau=0  im_prev=EMPTY                          -> zero (trivial)
k=1 (building phi_beta[2]): j=1 x=467 g_prev_pos=4 inv_tau=1 im_prev=EMPTY                          -> zero, but inv_tau=1!
k=2 (building phi_beta[3]): j=1 x=430 g_prev_pos=5 inv_tau=1 im_prev=EMPTY (pos 5 = G_2 cogen 0)    -> zero, but inv_tau=1!
```
**Notable recurring pattern**: `inv_tau=1` shows up at cogen-index-1's own chain
condition at BOTH k=1 and k=2 (two consecutive levels) — i.e. cogen 1 specifically (not
cogen 0, 2, etc.) seems to systematically involve a tau^1 relationship at multiple
levels. This is suggestive of a real structural/torsion pattern (not coincidence) but in
EVERY case so far `im_prev` is already empty at the point of division, so the division
itself is numerically inert (0 * tau^{-1} = 0) — we have NOT yet caught the "tau^m * X=0
for nonzero X" mechanism red-handed. The empty `im_prev` values themselves trace back to
`im_ck` being empty at k=0's j=0 and j=2 (i.e. qut_1 annihilating specific positions),
which is the SAME "qut(pivot)=0" pattern already established as looking legitimate.

**Status / open question**: we have a fully traced dependency chain from alpha's class
all the way back to level 0, and a recurring tau^1 signature at cogen-index-1 across
multiple levels, but have not yet found a point where this actually multiplies a
genuinely nonzero value by a torsion-killing tau power (the smoking-gun moment the user's
hypothesis predicts). Candidates for why: (a) the real collapse happens even earlier /
elsewhere not yet checked (e.g. beta_rep's own seeding, or the {1-1} side of the product
rather than the {3-1} side), (b) the recurring inv_tau=1 pattern is a red herring / just
ordinary bigrading and not actually torsion-related, (c) need to check a DIFFERENT
cogen's chain (not cogen 1) since alpha's OWN j=1 at k=2 might not be the only relevant
path — build_M/M_beta[3] might combine MULTIPLE cogens' phi_beta[3] values, not just j=1
in isolation.

## Reopening the unify()/inverse() Gaussian-elimination theory — found a real gap in my
## earlier reasoning

User asked whether the tau^{-1} investigation connects back to the EARLIER `unify()`/
`tauOper::inverse()` theory (matrices_mem/4.h, used by `gaussian`/`row_reduction` during
`del_and_gaussian`'s Gaussian elimination) — the one I'd previously (session-early)
argued was a no-op in practice because `embed2cofree`'s `symplify_to_led` only ever
accepts pivots with an ALREADY-invertible (tau^0) coefficient before inserting them into
the table.

**That reasoning has a real gap for the MODELED path specifically**: it's true for plain
`resolvor`/`embed2cofree` (F2 or tauPoly built directly with its own `symplify_to_led`
check). But `resolvor_modeled` is different: `embed2cofree_modeled` does NOT use
`symplify_to_led` at all (builds inj directly from a fixed `gens` list via `adjoint()`).
The pivot COLUMN (`inj_ind[i]`) comes from the CLASSICAL (F2) model table's recorded
`itm.cycle`, guaranteed invertible only in the F2, classical sense, at the time `mr_ex`
built it. But the actual MOTIVIC VALUE at that column, in the freshly-reconstructed
`indj` row (`new_cyc = inj->maps_to(transformer(itm.full_tag))`), is a completely
independent computation — nothing guarantees `new_cyc`'s coefficient AT column
`inj_ind[i]` is tau^0! It's whatever `inj`'s (genuinely tau-graded) rows produce. If it's
NOT tau^0, then `unify(inj_ind[i], indj row i)` calls `tauOper::inverse()` on a truly
non-tau^0 value — and since `inverse()` doesn't check `invertible()` (returns `-x`
unconditionally), this WOULD produce a bogus scaling factor and corrupt the whole row by
scaling it with a meaningless "tau^{-n}" (n!=0) multiplier. This is a live, previously
un-verified gap — not yet proven to occur, but structurally very plausible and untested.

**RULED OUT, definitively, with a complete scan (not just one example)**: added a
systematic `[PIVOTSCAN]` check in `resolvor_modeled` (hopf_algebroid/6.h, env-gated
`TRACE_PIVOT_SCAN`/`TRACE_FRANK_SCAN`) that checks, for EVERY tag `i` at a given level,
whether `indj`'s row `i` has an invertible (tau^0) coefficient at its own pivot column
`inj_ind[i]` — exactly the precondition needed for `unify()`'s buggy `inverse()` call to
be safe. Ran across ALL 7 levels relevant to {3-1}*{1-1} (F.rank = 1056, 3414, 4633,
4105, 3669, 3531, 2968, covering resolution levels 0-6):
```
every level: non-invertible-self-coefficient count = 0    self-column-missing count = 0
```
**Zero exceptions across 2274+2359+1055+1831+1838+1693+1 = ~9051 total tags checked.**
This DEFINITIVELY rules out the `unify()`/Gaussian-elimination tau-corruption theory —
not just for the one example traced by hand, but exhaustively for every pivot at every
level that matters for this product. The "master plan" tau-division-during-elimination
mechanism is confirmed NOT the source of the bug.

**So the remaining live leads are**: (a) the cogen-index-1 `inv_tau=1` pattern recurring
across levels (still not caught in the act multiplying a nonzero value), or (b) some
other, not-yet-identified "different kind of tau-trickery" per the user's suggestion —
possibly in `lift_one_step`'s own `tau^{-n}` division (lift.h, NOT the Gaussian
elimination), or in `cofree_adjoint_row`'s handling of the coaction, or somewhere in
`embed2cofree_modeled`/`adjoint()` itself (not yet examined in this much detail).

**User chose to explore the LESS-explored angle first**: not `lift_one_step`'s own
tau^{-n} division (already partially traced via TRACEB), but the OTHER, unexamined
candidate: `adjoint()` — the function used both by `embed2cofree_modeled` (to build inj
rows directly from the fixed `gens` list, WITHOUT any `symplify_to_led` check) and by
`cofree_adjoint_row` (to expand the coaction when constructing phi_beta). This is the
core function translating the classical model's combinatorial structure into actual
tauPoly-graded values — the natural place for a genuine tau-bookkeeping bug distinct
from both the (ruled out) Gaussian-elimination theory and the (partially explored)
lift_one_step division.

Found two overloads (hopf_algebroid/9.h + declarations in hopf_algebroid/4.h:43,47):
- `cofree_comodule<...> adjoint(const CoModule *X, int n, matrix<ring> *adjoint_map, int shift)` (9.h:4)
- `vectors<matrix_index,ring> adjoint(const CoModule *X, std::vector<int> const& gens, std::vector<uint32_t> const& pos, int i, int shift=0)` (9.h:25)

**Read adjoint() (hopf_algebroid/9.h)**: both overloads just EXTRACT existing algebroid
coefficients from `X->coaction(i)` via `component()` and convert via `algebroid2vector()`
— they don't themselves compute/assign tau exponents from scratch, they faithfully
propagate whatever X's own coaction structure already has. Not obviously buggy by
inspection; the real tau-assignment logic must be further down, in `algebroid2vector`
itself.

## FOUND THE CORE TAU-BOOKKEEPING MECHANISM — `tauVal`

Read `MotSteenrodOp::algebroid2vector` (mot_steenrod.cpp:173-179) and
`MotSteenrodOp::tauVal` (mot_steenrod.cpp:157-162) — this is the actual, pervasive
"different kind of tau trickery" the user hinted at, distinct from both the (ruled out)
Gaussian-elimination theory and the (partially explored) lift_one_step division:

```cpp
// The tau valuation of the monomial e. x_i^2 can be divided by tau, so has valuation 1
int MotSteenrodOp::tauVal(exponent e) {
    int res=0;
    for(int i=1;i<=maxVar;i++)
        res += xnVal(e,i)/2;     // integer division — only EVEN part of each exponent counts
    return res;
}

vectors<matrix_index,tauPoly> MotSteenrodOp::algebroid2vector(motSteenrod const &x,int shift){
    ... return std::pair<matrix_index,tauPoly>(mon_index[e]+shift,
             this->ringOper->multiply(tauPoly(tauVal(e)),r)); ...
}
```

**What this reveals**: the classical monomial basis (built by `e2p`'s `doubling_ex`,
which converts $z_i x_j$-style classical exponents into "doubled" `x_i` exponents) is
NOT in 1-1 correspondence with tauPoly-graded motivic basis elements directly — each
classical monomial `e` carries an IMPLICIT tau power, `tauVal(e)`, computed by counting
`floor(exponent_i / 2)` per variable and summing. Converting a classical algebra element
into its motivic vector representation (`algebroid2vector`) MULTIPLIES the given
coefficient by `tau^{tauVal(e)}` — i.e. this is exactly where "extra" tau powers get
introduced into the computation, driven by a per-monomial COMBINATORIAL formula, not by
any genuine Adams/motivic-weight bookkeeping tied to torsion. The inverse conversion
(`vector2algebroid`, mot_steenrod.cpp:189-195) undoes this via `tau^{-tauVal(e)}`, and
`generate_cofree_coaction` (mot_steenrod.cpp:221) applies YET ANOTHER
`tau^{-tauVal(mon_array[n])}` scaling when building the base coaction table.

This is a genuinely dense, easy-to-misread piece of bookkeeping — exactly matching the
user's description of the "opaque" tau tricks. Have NOT yet found a concrete bug in it,
but this is now the strongest remaining candidate for where a torsion-specific mistake
could hide, since it's a per-monomial combinatorial rule (based purely on exponent
parity) rather than anything computed FROM the actual resolution structure — if this
rule ever disagrees with the TRUE tau valuation for a class that's genuinely torsion
(where the "well-defined tau power" question is subtler than for a free/generic class),
that mismatch would silently propagate through every downstream computation.

**Found and read `xnVal`/`unpack`/`doubling_ex`**: `tauVal(e) = sum_i floor(xnVal(e,i)/2)`.
Cross-referencing with `ex2poly.cpp`'s `doubling_ex` (`res[i] = eps.first[i]*2 +
eps.second[i]`, converting classical exponents into motivic ones): a monomial built this
way has `xnVal(e,i) = 2*eps.first[i] + eps.second[i]` (eps.second[i] in {0,1}), so
`floor(xnVal(e,i)/2) = eps.first[i]` exactly — i.e. `tauVal(e)` is precisely the sum of
the classical "doubled" (z_i-style) exponents. This is a STANDARD, well-established
technique for encoding A_C/tau's associated-graded structure (each doubled generator
genuinely corresponds to one hidden factor of tau in the real motivic dual Steenrod
algebra) — this is foundational infrastructure used throughout ALL motivic computations
in this codebase, not something specific to resolvor_modeled. Given it's used pervasively
and other computations work, and given we ALREADY independently confirmed `qut_3`
correctly encodes the REAL tau^1 Bockstein relationship (`{3-2}->{4-1}`) using this exact
system, this mechanism looks sound — not proven bug-free, but not an obvious culprit
either. Not pursuing further without a more specific lead into this system.

**DEFINITIVE NEGATIVE RESULT on the "tau^m * X = 0" hypothesis for lift_one_step's own
division**: added a full scan (`TRACE_LIVE_DIV` env var) across ALL cogens at ALL levels
(0 through maxlev) for the {3-1}*{1-1} computation, looking for ANY case where
`inv_tau != 0` AND `im_prev` is simultaneously non-empty (the only situation where the
`tau^{-n}` division in lift_one_step could numerically matter). Result: **zero matches,
anywhere in the entire computation**. The division is NEVER exercised on a nonzero value
for this product — every time `inv_tau` is nonzero, `im_prev` happens to already be
empty (inert 0 * tau^{-n}), and every time `im_prev` is nonzero, `inv_tau` happens to be
0 (no division applied). This DEFINITIVELY rules out `lift_one_step`'s own tau^{-n}
division as the cause, for this example — not just the one path traced by hand, but
exhaustively across the whole computation.

**Status check — both concrete "tau-division" theories now ruled out**: (1) the
Gaussian-elimination `unify()`/`inverse()` theory (ruled out: every pivot self-coefficient
is invertible, zero exceptions across ~9051 tags), (2) `lift_one_step`'s own tau^{-n}
division (ruled out: never applied to a nonzero value anywhere in this computation). The
recurring `inv_tau=1` pattern at cogen-index-1 (levels 1,2) is real but NUMERICALLY INERT
in both cases — interesting structurally (probably reflects genuine torsion nearby) but
not itself the mechanism causing the wrong zero.

**Reassessment**: the actual deficiency must be elsewhere — most likely simply that
`phi_beta` isn't propagating enough nonzero content from level 0 forward (a sparsity/
richness problem, not a sign-or-scaling corruption). Given `im_ck` (qut applied to a
nonzero `im_prev`) was EMPTY at k=0 for cogens 0 and 2 of G_1 (positions 975, 979) — this
was earlier judged "looks like an ordinary pivot, probably legitimate" but was NOT
checked with the same rigor as the `to_del` and pivot-scan investigations. Given how many
"looks legitimate" assumptions turned out to be exactly where bugs were hiding this
session (to_del, and almost the Gaussian elimination pivot check), this deserves the same
level of scrutiny: is qut_1(975)=0 and qut_1(979)=0 actually forced by real exactness, or
could these ALSO be symptomatic of some other still-undiscovered issue specific to how
level 1's resolvor_modeled handles a nearby torsion class?

**Traced positions 975 and 979 directly (TRACE_RM, TRACE_FRANK=3414)**:
- 975 (deg=3,weight=1): pure self-only pivot (tag=3), `indj` row is just `{975^t0}`
  before AND after gaussian — completely ordinary, qut(975)=0 looks legitimate, same
  pattern as the earlier (verified-OK) 4070 case.
- 979 (deg=6,weight=3): pivot (tag=8), raw row BEFORE gaussian = `{979^t0, 1873^t0}` (two
  terms), AFTER gaussian = `{979^t0}` only (1873 eliminated via some other pivot, to_del
  confirmed empty/uninvolved) — same shape as the earlier 3040/4070 case, which we
  concluded was ordinary valid elimination, not a bug.
Both look unremarkable — no fresh evidence of a to_del-style or elimination-style defect
here specifically.

**Next step**: re-verify a foundational assumption for THIS beta (`{1-1}`) specifically:
session 3 confirmed `inj_row`s at G_1's cogens are pure `tau^0` singletons for the
`{6-9}`/`{4-6}` case, and the plan's `lift_one_step` formula relies on this (no
correction-term handling). Given beta is now `{1-1}` (bs=1, different resolution
sub-structure engaged via `R[k+bs]`/`R[k+bs+1]` indices), re-check this assumption holds
here too — if `inj_1` (or `inj` at other relevant levels for bs=1) ever has a NON-pure
cogen row (extra terms beyond the self singleton), that would directly justify bringing
back SOME form of correction-term handling (though carefully, not the buggy session-2
version) — this was previously verified only for bs=6, never for bs=1.

## Instrumentation currently in the tree

`hopf_algebroid/6.h`'s `resolvor_modeled` has TRACE_RM/TRACE_FRANK/TRACE_P1/TRACE_P2
env-var-gated debug prints (harmless when the env vars are unset — default off). Should
be removed once the investigation concludes. `yoneda2.cpp` also still has TMPCHECK0/
TMPCHECK/TMPCHECKNEG/TRACEB from earlier in the session (see debug_notes_archive1.md for
what these do) — also still present, also currently harmless/gated.

## Next step (per user, 2026-07-05)

Find the SMALLEST concrete example of `qut_i(inj_i(y)) != 0` in the motivic resolution
and investigate that one instance in full detail (rather than aggregate statistics across
thousands of rows). Likely candidate: level 1 (the first level where mismatches appear —
level 0 is rank 1 and trivially clean) has 1010 qut violations out of 3414; the SMALLEST
row index among those violations is a good starting candidate. May also be worth
rebuilding at a much smaller max_deg (e.g. 10 or 12) to get a small enough resolution to
trace completely by hand, rather than picking one row out of a large md=40 resolution.
