# Debugging Notes: Yoneda Product Bugs in yoneda2.cpp (fresh cycle)

Full history through this point is archived in `debug_notes_archive1.md` through
`debug_notes_archive4.md` (chronological). This file is the current, active summary — keep
it short; append detail only as it becomes load-bearing. For a higher-level (less
granular) summary of where things stand, see `STATUS.md`.

**Hard constraint, unchanged across all sessions (still the standing default)**: do not
read, edit, or run `yoneda.cpp` (no "2") for any reason. Only `yoneda2.cpp` and its
supporting headers are in scope. **One-time exception granted this session**: the user
explicitly asked to check whether `yoneda.cpp` (both working-tree and git `HEAD`) gets the
`φ_1(τ_1[1-0])` case right — see "Checked yoneda.cpp" below. This was a one-off,
explicitly requested check, not a standing change to the rule above; treat the
prohibition as back in force unless the user says otherwise again.

## Current status (2026-07-11 — IMPLEMENTING the approved plan)

**The root cause is understood, fix designed and approved.** Now implementing
`~/.claude/plans/look-at-the-top-level-smooth-llama.md` (rewrite of `lift_one_step_sum`/
`lift_first_step_sum` in `lift.h` to use the target's cofree universal property instead of
the source's). See `STATUS.md` for the math summary.

**This session so far:** re-read the plan and `lift.h` in full. Confirmed
`lift_first_step_sum`/`lift_one_step_sum` (the `tauPolySum` variants, `lift.h:120-236`) are
the ones actually invoked by `yoneda2.cpp` (`yoneda2.cpp:579,596`) — the plain
`lift_first_step`/`lift_one_step` (`lift.h:308-476`, `tauPoly`-only) appear unused by
`yoneda2.cpp` and are lower priority (only rewrite if still referenced elsewhere; check
before touching). Confirmed `gaussian` (`matrices_mem/4.h:33`) and `symplify_to_led`
(`matrices/9.h:34,68`, dense/sparse variants) exist and are the intended reuse targets for
the τ⁰-pivot echelon reduction of `inj_src` rows (Step 2 of the plan).

**Major design clarification worked out this session (resolves earlier confusion, no
code yet — see below for the concrete next action):**

- Re-derived from scratch what `cofree_comodule::coaction(i)` (`hopf_algebroid/13.h:63`)
  and `cofree_adjoint_row`/`cofree_adjoint_row_sum` (`lift.h:15,79`) actually compute.
  `coaction(i)` returns RAW F_src positions (same-owner-shifted local offsets from the
  algebroid's own `cofree_coaction`), NOT cogenerator indices. `cofree_adjoint_row`'s
  original doc comment ("lg values are cogenerator indices of F_tgt") confirms `lg`'s
  DOMAIN is all raw F_src positions and its RANGE is target-cogenerator-index space —
  i.e. `cofree_adjoint_row` **already implements the target-cofree cofree-up formula**
  `phi = (1⊗phi_bar)∘ρ_src` from the plan, exactly. So `cofree_adjoint_row_sum` does NOT
  need restructuring at the formula level — only its VALUE representation needs fixing:
  the `_sum` variant was generalized to let `lg` hold raw F_tgt positions (via
  `vector2algebroid` decomposition, see its own comment at `lift.h:69-77`), which is
  exactly the piece the plan says to delete (plan "Files to change" item 3).
- **Real bug isolated to just two things**: (1) `lift_one_step_sum`/`lift_first_step_sum`
  populate `lg` (=phi_bar) only at SOURCE COGENERATOR positions (`cog_pos`,
  `lift.h:174,229`), never at the true pivot positions of `image(inj_k)` (which the plan
  notes are "generally non-cogenerator offsets") — so `phi_bar` is pinned at the wrong
  domain points entirely; (2) the correction-term/`c_self` machinery (`lift.h:184-219`)
  is a leftover attempt to patch this via impure-row subtraction, order-dependent and
  the actual source of the underdetermination bug.
- **Checked whether level-0 (`lift_first_step_sum`) needs the same pivot-finding fix**:
  traced `beta_rep_pos`'s construction (`yoneda2.cpp:547-549`) — it's built as
  `position_of_gens[tm.ind]` for cogenerator indices `tm.ind`, i.e. it is ALREADY
  pure-cogenerator-supported (only ever touches cogenerator positions of `G_bs`). Since
  the comodule-map condition forces `phi_0(generator)` to be primitive
  (`ρ_tgt(phi(gen))=1⊗phi(gen)`), which for this algebroid's cofree comodules means
  living purely in cogenerator space (per `PHI_EXTENSION_ISSUE.md`'s resolved
  discussion), `(ε⊗1)v = v` here is a no-op — **no pivot/echelon-finding is needed at
  level 0**; only `lg`'s value representation needs to change (raw position `v` →
  cogenerator-index space via `filtered_reindex`-style reindexing, matching what the
  plain non-sum `lift_first_step` (`lift.h:308-335`) already does, just lifted to
  `tauPolySum`). Confirmed NOT a case requiring user clarification.

**Next concrete action:** implement, in order:
1. `cofree_adjoint_row_sum` (`lift.h:78-116`): change `lg`'s value convention from raw
   F_tgt positions to cogenerator-index space (drop the `vector2algebroid` decomposition
   path entirely), matching original `cofree_adjoint_row`'s convention.
2. `lift_first_step_sum` (`lift.h:120-140`): reindex `v` to cogenerator-index space
   before storing into `lg` (safe/lossless per the above), keep "every cogenerator of
   G_src gets this value" logic.
3. `lift_one_step_sum` (`lift.h:145-236`): add a local tau^0-pivot echelon reduction of
   `{inj_src.find(gens_k[j])}` (need to check `ModuleOp`/`vectors` API for a clean way to
   do this — next tool calls), assign `phi_bar[pivot_j] = prescribed(x_j)` (counit-project
   `im_gk` via cogenerator reindex, no `tau_inv` scaling needed beyond what pivot-unify
   requires), delete the `c_self`/correction-term loop entirely, add a hard abort if no
   tau^0 pivot is found for some row (mirroring `recover_inv_ind`'s abort style).
No code modified yet this session; this is a corrected/refined understanding of the
plan's math superseding my earlier (wrong) worry that `cofree_adjoint_row_sum`'s formula
itself needed to change or that level-0 needed a real pivot search.

## Implementation done this session (2026-07-11, later)

Implemented all three code changes above in `lift.h`:
- `cofree_adjoint_row_sum` (`lift.h:80`) simplified to the pure-cogenerator-value path
  (dropped `vector2algebroid`/local-offset decomposition).
- `lift_first_step_sum` (`lift.h:117`) now reindexes `v` into cogenerator-index space
  (hard-aborts if `v` ever has non-cogenerator support, which shouldn't happen given how
  `beta_rep_pos` is built in `yoneda2.cpp:547-549`).
- New `echelon_pivots_tau0` helper (`lift.h:165`) does forward-only Gauss-Jordan
  elimination (pivot identity only, no back-substitution needed) on
  `{inj_src.find(gens_k[j])}` to find tau^0 pivot positions.
- `lift_one_step_sum` (`lift.h:212`) rewritten: pins `lg[pivot[j]] = prescribed(x)` (no
  more `cog_pos`, `c_self`, or correction-term loop). Added a hard-abort check that
  `inv_tau_prev[x]` is always tau^0 (replacing the old blind `tau_inv` division).

`cmake --build build --target yoneda2` succeeded clean (only pre-existing, unrelated
warnings: undefined inline `tauOper::multiply`/`isZero` -- those must be defined
elsewhere and linked, not a new issue; unused `fmt` function in yoneda2.cpp).

**Next concrete action**: per CLAUDE.md pitfall #11, must do all three before testing:
(1) already rebuilt; (2) `\cp -f build/yoneda2 ./yoneda2` (note: plain `\cp`/`\rm` still
hit the interactive alias/function in this shell despite the leading backslash -- had to
use `/bin/rm -f` / literal path to actually bypass it, contradicting CLAUDE.md pitfall
#12's claim that `\cmd` suffices; worth revisiting that note); (3) cleared phi caches.

## MAJOR FINDING: verification surfaced a real gap in the just-implemented rewrite

`TRACE_SQUARE_CHECK=30 ./yoneda2 40 30 4 6 6 9` still fails at the SAME spot as before
the rewrite: level k=0, position p=6, DIFF={7-10} missing from LHS. Traced it down:

- `TRACE_PHI_AT_LVL=1 TRACE_PHI_AT_POS=4` (=τ_1[1-0]) still gives **0**, not the expected
  `{7-10}` from `PHI_EXTENSION_ISSUE.md`.
- Root cause found by dumping `inj_1`'s rows for the 6 entries of `gens[1]`
  (`TRACE_INJ_AT_LVL=1 TRACE_INJ_AT_X={0,1,4,15,71,482}`): **every one of them is a pure
  singleton** `t^0*{1-j}` -- no impurity at all among `gens[1]`.
- But `PHI_EXTENSION_ISSUE.md:288` traces a DIFFERENT, impure row:
  `TRACE_INJ_AT_LVL=1 TRACE_INJ_AT_X=5` (raw X_1 index **5**, literal) gives
  `t^0*{1-0}+4 + t^1*{1-1}+2` -- reconfirmed exactly with today's build/args
  (`40 30 4 6 6 9`). And `TRACE_QUT_LVL=0 TRACE_QUT_POS=6` gives `x5^t0`, matching
  `PHI_EXTENSION_ISSUE.md`'s `qut_0(6)=x5` and explaining `TRACE_SQUARE_CHECK`'s failure
  at exactly position p=6.
- **The X_1 index 5 is NOT one of `gens[1]`'s six values (0,1,4,15,71,482) at all.**
  `gens[k+1]` (passed as `gens_k` to `lift_one_step_sum`) is the MODULE-GENERATING
  subset used by `embed2cofree_modeled` to build `G_{k+1}` (one cofree summand per
  generator, `hopf_algebroid/7.h:160-169`) -- it is NOT the full basis of `X_{k+1}`.
  `inj_{k+1}` (`R[k+1].inj`), however, is defined and meaningful on **all** of
  `X_{k+1}.rank()` (`inj->set2zero(X->rank())`, `hopf_algebroid/7.h:154`), including
  non-generator indices like `x=5`.
- **This means my rewritten `lift_one_step_sum` (and the old one before it) only ever
  processes the `gens_k` subset (6 rows here) when it needs to process ALL of
  `X_{k+1}`'s rank** to correctly pin `phi_bar` on the whole of `image(inj_{k+1})` --
  CLAUDE.md's own worked example ("tau_1[1-0] is the pivot of inj_1(x5)") literally names
  a non-generator index, which I'd misread earlier as "gens_k[5]" -- confirmed wrong by
  directly checking both queries give different answers.

**Ruled OUT**: this is not a bug in the echelon-pivot logic itself (it correctly found
pivot=cog_pos for each of the 6 pure rows it was given -- correct given pure input). Also
not a bug in `phi_0`/level-0 handling (`phi_beta[0].find(1) = 2384^t^0`, nonzero and
presumably correct) or in `qut_6`/`inj_7` (confirmed `qut_6.find(2384)=0` directly,
consistent, not obviously wrong on its own -- though this chain needs re-checking once
X_1 index 5 is actually processed, since prescribed(x=5) is likely what should carry the
`{7-10}` answer, not prescribed(x=0)).

**Next concrete action**: extend `lift_one_step_sum` (and its pivot-finding) to run over
**all** of `X_{k+1}`'s rank, not just `gens_k`. Need to find how to obtain
`X_{k+1}.rank()` inside `lift.h`/at the call site in `yoneda2.cpp` (candidates: a `Xrank`
field already used elsewhere in `yoneda2.cpp`, e.g. `R[k+1].Xrank` per
`recover_inv_ind(R[k].qut, R[k].F.total_rank, (unsigned)R[k+1].Xrank)` at
`yoneda2.cpp:592-593` -- this is exactly the quantity needed, already available at the
call site). Plan: add an `unsigned X_rank` parameter to `lift_one_step_sum`, loop
`x=0..X_rank-1` for both the echelon-pivot pass and the prescribed/lg-insert pass
(instead of looping over `gens_k`), and pass `R[k+1].Xrank` from `yoneda2.cpp`. Was
mid-investigation when checking `matrix<ring>`'s own rank field
(`matrices.h`/`matrices/*.h`) as an alternative source of truth for `X_rank` -- not yet
resolved which is more reliable; check next.

## X_rank fix implemented and real progress made (2026-07-11, still later)

Implemented: `lift_one_step_sum` now takes `unsigned X_rank` instead of `gens_k`, loops
`x=0..X_rank-1` for both the echelon-pivot pass and the prescribed/lg pass
(`lift.h:211-284`). Call site updated (`yoneda2.cpp:596-598`) to pass
`(unsigned)R[k+1].Xrank` instead of `gens[k+1]`. Rebuilt, redeployed
(`/bin/cp -f build/yoneda2 ./yoneda2` -- confirmed `\cp -f` still doesn't bypass the
alias in this shell, must use `/bin/cp`/`/bin/rm` with the literal path), caches cleared
(`/bin/rm -f *_yoneda2_*_phi*`).

**Result: real, substantive change.** `TRACE_SQUARE_CHECK=30 ./yoneda2 40 30 4 6 6 9`:
the previous failure at level k=0, p=6 is GONE (no longer reported). The product
`{4-6}*{6-9}` is no longer trivially `0` -- it now gives `t^0{10-25}+t^0{10-26}`
(nonzero, and `{10-25}` matches the long-expected correct answer from earlier sessions'
notes, e.g. `STATUS.md`'s "both {10-25}"). **But** a NEW `TRACE_SQUARE_CHECK` failure
now appears further along, at level k=0, p=14:
```
qut_0(14): 13^t0
LHS = phi_1(inj_1(qut_0(14))): 1951^t^0({7-9}+7) 2141^t^0({7-12}) 2178^t^0({7-13})
RHS = inj_7(qut_6(phi_0(2397))):                 1951^t^0({7-9}+7)
DIFF: 2141^t^0({7-12}) 2178^t^0({7-13})   (LHS has 2 EXTRA terms RHS lacks)
```
Different shape of failure than before (LHS is a superset of RHS, not missing a term) --
suggests either phi_1 is now producing spurious extra terms at some position (possibly
an echelon/pivot-order artifact, or a genuine remaining gap in the construction), or
there's a second, distinct bug. Also gives an extra `{10-26}` term in the product beyond
the expected `{10-25}` -- consistent with "some phi_1 value has spurious extra content."

**Reverse-order product checked**: `./yoneda2 40 30 6 9 4 6` gives `t^0{10-25}` cleanly
(matches the historical "known correct" value, no spurious term) -- so the bug is
specific to this direction's `phi_beta` build (or at least not (yet) manifesting in the
reverse one).

**Also reconfirmed the ORIGINAL known-bug case is now fixed**:
`TRACE_PHI_AT_LVL=1 TRACE_PHI_AT_POS=4` (τ_1[1-0]) now gives `2028^t^0` = `{7-10}`,
exactly matching `PHI_EXTENSION_ISSUE.md`'s hand-derived expected answer. Real progress.

## ROOT CAUSE FOUND for the p=14 failure: back-substitution is NOT optional after all

Traced fully: `qut_0(14)=13`ᵗ⁰(sing.), `inj_1.find(13) = t^0*{1-1}+7 (pos 981) + t^0*{1-2}+4 (pos 1875)`.
Separately, `qut_0(15)=14`ᵗ⁰(sing.), and **`inj_1.find(14) = t^0*{1-2}+4 (pos 1875) ALONE`**
-- i.e. row 14 is a PURE singleton at the SAME position (1875) that also appears in row
13. `echelon_pivots_tau0` processes rows in x-order (13 before 14), so when row 13 is
processed, position 1875 hasn't been claimed yet -- row 13 picks pivot=981 (the smaller
index) without ever being reduced against row 14's (not-yet-known) pivot at 1875. Then
row 14 -- forced, no choice -- claims pivot=1875.

**This is the bug**: I earlier argued (see the `echelon_pivots_tau0` comment,
`lift.h:150-163`) that skipping back-substitution was fine because forward-only
elimination already gives an invertible (triangular) change-of-basis matrix -- true for
BASIS VALIDITY, but WRONG for correctly assigning **values**. The actual constraint is
`phi_bar(inj_1(x)) = prescribed(x)` for the ORIGINAL (unreduced) row. Row 13's original
content is `e_981 + e_1875`, so the constraint is
`phi_bar(981) + 1*phi_bar(1875) = prescribed(13)`, i.e.
**`phi_bar(981) = prescribed(13) - phi_bar(1875) = prescribed(13) - prescribed(14)`**,
NOT `prescribed(13)` alone (what the code currently stores). Confirmed numerically:
`phi_1(981)=1951` (only, matches RHS exactly) but `phi_1(1875)=2141+2178` (the exact
spurious DIFF terms) -- i.e. `lg[981]` is currently `prescribed(13)` unadjusted, and the
missing subtraction of `phi_bar(1875)=prescribed(14)` is exactly what's leaking through
as the extra terms in `phi_1(inj_1(13)) = phi_1(981)+phi_1(1875)`.

**Fix (not yet implemented)**: after the forward pass in `echelon_pivots_tau0` (which
already leaves `rows[x]` reduced -- zero at every EARLIER-claimed pivot, by
construction, so the only "unknown" leftover terms in a forward-reduced row are LATER
rows' pivots or genuinely-free non-pivot positions), do a SECOND, BACKWARD pass in
`lift_one_step_sum` (highest x down to 0): `phi_bar(pivot[x]) = prescribed(x) - sum over
(q,coeff) in the FORWARD-REDUCED row (q != pivot[x], q happens to be pivot[x'] for some
x') of coeff * phi_bar(pivot[x'])` -- valid because x' is guaranteed > x (earlier pivots
were already eliminated from the forward-reduced row), so `phi_bar(pivot[x'])` is already
computed by the time the backward pass reaches x. Requires keeping the FORWARD-REDUCED
rows (currently discarded after `echelon_pivots_tau0` returns just the pivot list) and a
`pivot_of_col`-style reverse map from column to owning row index.

**Next concrete action**: implement the backward substitution pass, rebuild, redeploy,
reclear caches, re-run the full verification suite (`TRACE_SQUARE_CHECK=30`, the τ_1[1-0]
check, both product orders, and the `40 30 1 1` regression).

## Back-substitution implemented -- more real progress, one more failure surfaced

Implemented per the plan above (`lift.h:245-303`): compute `prescribed(x)` for all
`x` independently first, then a BACKWARD pass (`xi` from `X_rank-1` down to `0`) using
the forward-reduced `inj_rows` (kept, not discarded) and an `owner_of_pivot` reverse map,
subtracting already-resolved LATER pivots' `phi_bar` contributions. Rebuilt, redeployed,
cleared caches.

**Result: the p=14 failure is GONE, and `{4-6}*{6-9}` now gives the clean, correct
`t^0{10-25}`** (previously `t^0{10-25}+t^0{10-26}`) -- matching the reverse-order product
`{6-9}*{4-6}` exactly (both now `t^0{10-25}`). **This is likely the actual fix the user
cares about for this specific product pair.**

`TRACE_SQUARE_CHECK=30` still finds a failure, now further along at level k=0, p=27:
```
qut_0(27)=26(t^0);  inj_1(26) = t^0{1-0}+22  t^1{1-1}+16  t^0{1-1}+18
                                t^1{1-2}+9   t^0{1-2}+11  t^1{1-3}+2
LHS: 1960^t1({7-9}+16) 1962^t0({7-9}+18) 2037^t0({7-10}+9) 2271^t1({7-16})
RHS: 1960^t1({7-9}+16) 1962^t0({7-9}+18) 2037^t0({7-10}+9)
DIFF: 2271^t1({7-16})   -- from phi_1(pos 1880 = {1-2}+9), which equals prescribed(x=25)
```
**Traced fully**: pos 1880 is pivot of X_1 index **25** (row `inj_1(25) = t^0{1-2}+9
t^0{1-2}+10 t^0{1-3}+2`, pivot=1880=smallest). The row's OTHER two positions (1881,
2630) were checked against all ~1055 pivots (`TRACE_LIFT=2000`, covers full X_rank) and
are NOT anyone's pivot -- genuinely free, so `phi_bar[25] = prescribed(25)` exactly, NO
back-substitution correction applies here. So `phi_bar[25]={7-16}` is exactly what
`prescribed(25)` says, and prescribed/phi_bar bookkeeping is internally consistent --
**this isn't a bug in the echelon/back-substitution machinery itself.**

**Open question (ruled out simple explanations, not yet resolved)**: is `prescribed(25)`
*itself* wrong, i.e. is `im_gk` (built from `phi_0` at X_1 index 25's recovered
preimage in G_0) not actually the value the chain condition requires? Live hypothesis:
`lift_first_step_sum`'s level-0 construction ("every cogenerator of G_0 maps to
`beta_rep_pos`, cofree-up the rest") has NO chain condition to satisfy (there is no
level `-1`), so nothing currently forces `phi_0(p1) - phi_0(p2)` to land in
`image(inj_bs)` whenever `qut_0(p1) = qut_0(p2)` for two DIFFERENT G_0 positions p1,p2 --
but exactly this property is implicitly REQUIRED for `prescribed(x)` to be
preimage-independent (the `d∘d=0` argument in CLAUDE.md implicitly assumes phi_{k-1} is
"already a genuine chain map," which for k=0 needs its own justification, not just
inherited for free). Haven't yet confirmed whether this is the actual gap, or whether
`recover_inv_ind`'s "first singleton wins" preimage choice for X_1 index 25 is itself
inconsistent with what TRACE_SQUARE_CHECK's direct `qut_0(27)=26` route implies -- these
involve DIFFERENT X_1 indices (25 vs 26) so may not even be directly comparable; need to
re ORIENT/recheck this reasoning fresh next session rather than trust this half-formed
hypothesis.

**This is a math-level open question, not a mechanical implementation bug** -- per
CLAUDE.md's "stop and explain math-level plan changes" directive, this should be
reported to the user rather than silently patched further. **Status to report**: real,
verified progress (two genuine bugs found and fixed: missing non-generator X indices,
and missing back-substitution), the originally-tested product pair `{4-6}*{6-9}` /
`{6-9}*{4-6}` now agrees and matches the historical expected value, but
`TRACE_SQUARE_CHECK` still finds a residual failure further into the resolution
(level 0, p=27) whose root cause is now narrowed down to either (a) `phi_0`'s
preimage-independence not actually holding, or (b) something in `recover_inv_ind`'s
preimage choice -- not yet distinguished. `./yoneda2 40 30 1 1` regression and the
`phi_1(τ_1[1-0])={7-10}` check have NOT been re-verified against this latest build yet
either -- do that first thing next session, before diving back into the p=27 question.

**Update**: ran both. `./yoneda2 40 30 1 1` still looks clean (all self-products `{n-0}
-> 0` at high filtration, matching the historical "all clean" regression pattern).
`phi_1(τ_1[1-0])` still correctly gives `{7-10}`. So of the plan's 4 verification items,
3/4 pass cleanly (τ_1[1-0] case, both product orders agreeing, the 1-1 regression);
only the full `TRACE_SQUARE_CHECK=30` sweep still finds the level-0/p=27 residual issue
described above.

## p=27 failure: ROOT CAUSE FOUND -- retracts the "math-level phi_0 gap" hypothesis above

The "phi_0 preimage-independence" hypothesis above was WRONG (or at least not needed to
explain this failure) -- **this is a plain third implementation bug**, same species as
the back-substitution bug, not a math-level issue. Full derivation:

**The missing piece: `prescribed`/rhs must be forward-eliminated TOGETHER with the
matrix rows, not just at the end.** `echelon_pivots_tau0`'s forward pass correctly
updates `rows[j]` when eliminating an earlier pivot column (`rows[j] -= coeff *
rows[owner]`), but `lift_one_step_sum` never does the SAME operation on `prescribed[j]`
during that forward pass -- `prescribed` is computed completely independently (from
`phi_prev`/`im_gk`) and consulted only in the backward pass as if untouched by row
26's own the forward reduction of `rows[26]` and NOT rhs). Standard Gauss-Jordan
requires transforming the augmented matrix `[M | b]` together; I only transformed `M`.

**Concrete confirmation, worked by hand and cross-checked against actual traced
values**: row 26 (`inj_1.find(26)` = `t^0{1-0}+22 t^1{1-1}+16(=990) t^0{1-1}+18(=992)
t^1{1-2}+9(=1880) t^0{1-2}+11(=1882) t^1{1-3}+2(=2630)`) gets column 1880 eliminated
using row 25 (`inj_1.find(25) = t^0*1880 t^0*1881 t^0*2630`, pivot=1880, unified coeff
1) during the forward pass: `rows[26] -= t^1 * rows[25]`, giving `rows[26]_reduced =
{22:t0, 990:t1, 992:t0, 1882:t0, 1881:t1}` (1880 and 2630 both cancel exactly, 1881
newly appears as a side effect of the subtraction). **The forward elimination step
implicitly used `t^1 * rows[25]`'s contribution -- but the CORRESPONDING `t^1 *
prescribed(25)` was never subtracted from `prescribed(26)`.** The correct equation
(solving `Σ_q M_orig[26][q]·φ̄(q) = prescribed(26)` with the now-known
`φ̄(1880)=prescribed(25)={16:t^0}` and `φ̄(2630)=0`, `φ̄(1881)`/`φ̄(992)`/`φ̄(1882)`
confirmed free/unclaimed, `φ̄(990)=prescribed(24)=0`) is:
```
φ̄(22) = prescribed(26) - t^1·prescribed(25) - t^1·prescribed(24)
       = 0 - t^1·{16:t^0} - 0  =  {16: t^1}
```
but the CURRENT code gives `phi_bar[26] = prescribed(26) = 0` (confirmed via
`TRACE_LIFT`: `x=26 pivot=22, prescribed: (empty), phi_bar: (empty)`) -- **missing
exactly the `{16:t^1}` term**.

**Verified this exactly explains the observed DIFF**: if `phi_bar(22)` correctly
included `{16:t^1}`, then `phi_1(22) = cofree_adjoint_row_sum(...)` would gain a NEW
direct term `{2271:t^1}` (2271 = `G_7.position_of_gens[16]` = `{7-16}`) on top of its
current confirmed value `{2037,2144,2145,2181,2182}` (all t^0, reconfirmed via
`TRACE_PHI_AT_LVL=1 TRACE_PHI_AT_POS=22`). Since row 26's own coefficient at position 22
is `t^0`, this new term propagates into LHS as `2271:t^1` -- landing on TOP of the
*already-present* `2271:t^1` contributed by position 1880's term in the same LHS sum
(`t^1 * phi_1(1880) = t^1*{2271:t^0} = {2271:t^1}`, the term that IS the current DIFF).
**Two `t^1` contributions at the same position cancel mod 2** -- so fixing `phi_bar(22)`
would make the spurious `2271^t^1({7-16})` term in LHS disappear entirely, resolving the
p=27 `TRACE_SQUARE_CHECK` failure. This is a clean, mechanistic, fully-traced
explanation, not a guess.

**Fix needed (not yet implemented, pending user go-ahead)**: track a parallel `rhs[x]`
(seeded as `prescribed[x]`, `tauPolySum`-valued) alongside `rows[x]` inside
`echelon_pivots_tau0`'s forward loop, updating `rhs[j] -= coeff * rhs[owner]` every time
`rows[j] -= coeff * rows[owner]` happens (owner's `rhs` is always already-finalized at
that point, since forward processing order = pivot-claim order). Then the BACKWARD pass
should start from this forward-adjusted `rhs[x]` (not the raw `prescribed[x]`) before
subtracting later-pivot corrections. Requires either extending `echelon_pivots_tau0`'s
signature to also carry/return the adjusted rhs (parametrized over the value type, since
rows are `ring`-valued but rhs is `tauPolySum`-valued -- can't literally reuse the same
loop body without genericizing over both types), or duplicating the forward-elimination
loop inline in `lift_one_step_sum` with both `rows` and `rhs` tracked together.

## rhs-propagation fix implemented -- moved the failure much further (p=27 -> p=112), found the REAL remaining bug

Implemented per the plan: `echelon_pivots_tau0` now takes an optional
`rhs` (`prescribed`) parameter and applies the identical row operation to it whenever a
row is reduced (`lift.h:167-236`); `lift_one_step_sum` computes `prescribed` BEFORE the
echelon call so it gets forward-adjusted in lockstep (`lift.h:255-303`). Rebuilt,
redeployed, cleared caches.

**Big improvement**: `TRACE_SQUARE_CHECK` failure moved from p=27 all the way to
**p=112** (product `{4-6}*{6-9}` still correctly `t^0{10-25}`). But NOT fully clean yet.

**New failure traced**: `qut_0(112)=111`, `inj_1(111) = t^0{1-0}+97(pos97)
t^0{1-1}+84(pos1058) t^0{1-4}+3(pos3161)`. `phi_1(97)={2089:t1}`,
`phi_1(1058)={2028:t0, 2089:t1}`, `phi_1(3161)=0`. Sum cancels the `2089:t1` parts,
leaving `{2028:t0}` = the DIFF. Traced `pos1058`'s owner via a full-budget
`TRACE_LIFT` dump and found **TWO different X_1 indices both reporting pivot=1058**:
`x=116` AND `x=578`(!). Pivots are supposed to be globally distinct (each raw position
claimed by exactly one row) -- this is a genuine invariant violation in
`echelon_pivots_tau0` itself, not a math-level issue.

**Root cause of THIS bug**: the forward-elimination loop does only a SINGLE PASS over
`pivot_of_col` per row. When eliminating pivot column `col` (owned by row `owner`) from
row `j`, the subtraction `rows[j] -= coeff*rows[owner]` can REINTRODUCE nonzero content
at a DIFFERENT already-claimed column `col2` if `rows[owner]` itself still has leftover
content there (`owner` was finalized before `col2`'s owning row existed, so `rows[owner]`
was never reduced against `col2`). If `col2 < col` in the map's iteration order, the
single forward sweep has already passed `col2` and won't revisit it -- so row `j` can
end up choosing an ALREADY-CLAIMED column as if it were free (exactly what happened:
row 578 apparently picked column 1058, already claimed by row 116, because some
intermediate elimination reintroduced content there after the sweep had passed it).
This is the exact same "byproduct of subtraction reintroduces a term" mechanism already
seen with position 1881 in the row 25/26 example -- but there it only affected a FREE
(non-pivot) position, harmlessly; here it corrupts pivot uniqueness itself.

**Fix (not yet implemented)**: switch to TRUE Gauss-Jordan (full RREF), not
forward-only elimination. Whenever row `j`'s own new pivot `chosen` is established,
immediately eliminate that column from EVERY other already-finalized row `j2 < j` too
(propagating the same operation to `rhs[j2]`), not just eliminate earlier pivots FROM
row `j`. This maintains the invariant that every finalized row has zero at every OTHER
row's pivot (both earlier- and later-claimed) at all times, which prevents the
reintroduction issue by construction. As a bonus, **this should make the separate
backward-substitution pass unnecessary**: after full RREF (with rhs correctly
propagated throughout), each row's equation reduces to exactly
`phi_bar(pivot_x)*1 = rhs[x]` (every other term is either an eliminated pivot column or
a genuinely free column contributing 0), so `phi_bar[x] = prescribed[x]` directly,
post-reduction, no backward loop needed. Complexity becomes O(n^2) in the number of rows
per level (~1000 rows at level 1 -> ~10^6 small-vector operations, should be fine).

## Full Gauss-Jordan implemented -- moved failure p=27->p=112, but SAME p=112 persists after the fix (unexpected)

Implemented full Gauss-Jordan (`lift.h`: retroactively eliminates each brand-new pivot
column from every earlier-finalized row too, not just forward elimination), with `rhs`
propagated through both directions. Simplified `lift_one_step_sum` to drop the
now-theoretically-unnecessary backward pass (`prescribed[x]` should directly BE
`phi_bar(pivot[x])` after true RREF). Rebuilt, redeployed, cleared caches.

**Result: unexpectedly, EXACTLY the same p=112 failure persists, unchanged.** Re-checked
pivot ownership: **both x=116 and x=578 still report `pivot=1058`**, even after switching
to full Gauss-Jordan (which should make this structurally impossible if implemented
correctly). This means either (a) the Gauss-Jordan fix has a bug of its own, or (b) the
p=112 failure was never actually caused by the pivot-duplication bug in the first place
(coincidentally unaffected by the fix) and I'm chasing two separate things.

**New data point**: dumped RAW (unreduced) `inj_1.find(578)` -- it does **NOT** contain
position 1058 at all (its positions are all offsets 479/481/390/395/396/248/259/261/262
/87/93/95 within owners 1-4, nowhere near offset 84 within owner 1 = 1058). So if row578
ends up with pivot=1058, that value can ONLY have arrived via elimination
"reintroduction" from some intermediate row's stale content -- confirming the mechanism
is real, but should have been prevented by the full-Gauss-Jordan retroactive step (per
my own re-derivation: every finalized row should stay continuously clean against every
OTHER established pivot, forward or backward, so no row should ever be able to "pick up"
an already-claimed column from another row used in elimination). This re-derivation
hasn't found the actual flaw yet -- added TEMPORARY env-gated debug output
(`TRACE_ECHELON_ROW=<row index>`, prints that row's full content after every processing
step from when it's first touched, `lift.h`, right after the retroactive elimination
loop) to directly watch row 578 across the whole sweep and see exactly which step (and
via which other row) introduces column 1058. **Not yet run** -- next action: run with
`TRACE_ECHELON_ROW=578`, grep for `1058` in the output to find the exact step it first
appears, then cross-reference which row was being finalized at that step.

**This instrumentation is temporary and should be removed once the bug is found** (it's
inside `echelon_pivots_tau0`, gated on `TRACE_ECHELON_ROW` env var).

## CORRECTION: the "duplicate pivot 1058" finding was a trace-scoping artifact, not a real bug

**Important gotcha discovered**: `lift_one_step_sum`'s `trace_budget` (`lift.h`, `static
int trace_budget = getenv("TRACE_LIFT")...`) is a **function-static variable, shared
across ALL calls to `lift_one_step_sum` for the entire program run** -- it does NOT
reset per level or per beta chain. `yoneda2.cpp`'s invocation `40 30 4 6 6 9` builds
**two separate `phi_beta` chains** (one for beta={4,6}, one for beta={6,9}), each
calling `lift_one_step_sum` once per resolution level (k=0,1,2,...) -- and **both
chains' level-1 step has the SAME `X_rank=1055`** (X_1's rank doesn't depend on which
beta is being lifted, only `inj_1` does, which is shared). So a `TRACE_LIFT=N` dump
mixes rows from **both beta chains' level-1 calls** indistinguishably when just grepping
by `x=`.

This is exactly what happened: `TRACE_LIFT_SUM x=578 pivot=1058` was NOT the same call
as `TRACE_LIFT_SUM x=116 pivot=1058` -- one belonged to the beta={4,6} chain's phi_1
build, the other to a *different* call (turned out to be a coincidence of two unrelated
dumps both landing on row-index 578 at various points, not actually the relevant
computation). Added `X_rank=` to the `TRACE_LIFT_SUM` print line (`lift.h`) and an
`TRACE_ECHELON_N` env-var (must match `X_rank` exactly) alongside `TRACE_ECHELON_ROW` to
properly scope the `echelon_pivots_tau0` debug dumps to a single call. **Re-verified,
properly scoped this time**: `x=116 pivot=1058` and `x=578 pivot=1453` -- **no
duplicate**. The pivot-uniqueness "bug" was never real; the fixed-point-elimination fix
and the row-578 investigation that prompted it were chasing a phantom caused by this
trace-scoping mistake.

**The full-Gauss-Jordan and fixed-point-elimination changes made while chasing this are
NOT reverted** -- they're legitimate defensive strengthening (the original single-pass
forward elimination genuinely COULD have this bug in principle, even though it turned
out not to be triggered here) and did not break anything (`TRACE_SQUARE_CHECK`, the
τ_1[1-0] check, and both product orders all still pass/match). Keeping them.

**Lesson for future sessions**: when using `TRACE_LIFT` (or any trace gated by a
function-static counter) to inspect a SPECIFIC resolution level/beta-chain's
`lift_one_step_sum` call, always cross-check against the `X_rank=` value now printed
alongside `x=`/`pivot=` -- do not assume two lines with the same `x=` value belong to
the same call just because a single `grep` run shows them near each other.

## Back to the ACTUAL p=112 failure (still open, correctly scoped this time)

Re-ran `phi_1` at positions 97/1058/3161 (the raw content of `inj_1(111)`, the actual
X_1 index under test at `TRACE_SQUARE_CHECK`'s p=112 failure) -- **unchanged from
before**: `phi_1(97)={2089:t1}`, `phi_1(1058)={2028:t0, 2089:t1}`, `phi_1(3161)=0`, sum
`={2028:t0}` = the DIFF. Properly-scoped check confirms **pivot 1058 belongs solely to
x=116** (no duplicate), with `phi_bar(1058) = prescribed(116) = {11:t^1}` (cogenerator
11 of `G_7`, coefficient t^1) -- i.e. `phi_bar(1058)` does NOT itself contain
cogenerator 10 (`{7-10}`=position 2028) at all. So the `2028:t0` term in `phi_1(1058)`
must arise from `cofree_adjoint_row_sum`'s expansion: position 1058 has a nontrivial
local offset (84) within owner-1's block, so `G_src.coaction(1058)` decomposes into
MULTIPLE `(algebroid, position)` pairs within that same block (not just the identity
term reaching `lg[1058]` itself) -- one of which apparently also reaches position 974
(owner-1's own base/cogenerator, where `lg[974]=prescribed(x=1)={9:t0}`, established
much earlier in the session) with some algebroid coefficient that, when expanded via
`algebroid2vector`, happens to land exactly on cogenerator 10's own base (position
2028) rather than merely re-scaling cogenerator 9. **Not yet fully traced through the
actual algebroid arithmetic** -- this requires either dumping `G_src.coaction(1058)`'s
raw terms directly (via a targeted trace hook, not yet added) or inspecting
`cofree_coaction`'s behavior at local offset 84, which goes beyond generic
linear-algebra bug-hunting into the specific Steenrod-algebroid comultiplication
structure. This is the next concrete thing to trace, but is a deeper, more
algebroid-specific investigation than the previous three bugs (which were all pure
linear-algebra/bookkeeping issues in `echelon_pivots_tau0`/`lift_one_step_sum`).

**Final re-verification after all three fixes (rhs-propagation, full Gauss-Jordan,
X_rank scoping)**: regression `40 30 1 1` still clean, `phi_1(tau_1[1-0])` still `{7-10}`
(correct), reverse-order product `{6-9}*{4-6}` still `t^0{10-25}` (matches forward
order). Only the full `TRACE_SQUARE_CHECK` sweep still finds the p=112 issue described
above. Temporary debug instrumentation (`TRACE_ECHELON_ROW`/`TRACE_ECHELON_N`,
env-gated, harmless when unset) left in `echelon_pivots_tau0` for next session's use —
should be removed once the p=112 root cause is found and fixed.

## p=112: full coaction trace + several hypotheses ruled out — genuinely deep now

Added a new temporary trace (`TRACE_COACTION_POS=<pos>`, env-gated, in
`cofree_adjoint_row_sum`, `lift.h`) dumping `F_src.coaction(i)`'s raw
`(algebroid, position)` terms with human-readable algebroid printing (via the global
`motSteenrod_oper.output(...)`, `mot_steenrod.h:69`), plus what each term contributes.
**Same trace-scoping gotcha as before applies** — `TRACE_COACTION_POS=1058` fires for
every call across the whole program where position 1058 occurs in ANY comodule, so the
first (small, 4-term) block in the output is the relevant one; a later 84-term block
belongs to an unrelated call and should be ignored.

**Full trace of the relevant call** (`coaction(1058)`, level 1, building phi_1):
```
COACTION_DUMP i=1058: 4 term(s)
  coaction term: pos=974  algebroid=o+t^-8x_1^17
  coaction term: pos=975  algebroid=o+t^-8x_1^16
  coaction term: pos=1046 algebroid=o+1x_1^1
  coaction term: pos=1058 algebroid=o+1
  lg.find(pos=974)  [phi_bar] = cogen9^[t^0]      -> contributes 2028^t0 ({7-10})
  lg.find(pos=1058) [phi_bar] = cogen11^[t^1]     -> contributes 2089^t1 ({7-11})
```
(positions 975, 1046 have no `lg` value -- genuinely free, correctly contribute
nothing). This exactly explains `phi_1(1058) = {2028:t0, 2089:t1}` mechanically: the
`{7-10}` term comes from expanding `phi_bar(974) = prescribed(x=1) = {9:t0}` through the
algebroid element `t^-8 * x_1^17` via `algebroid2vector` (`mot_steenrod.cpp:195`).

**Checked `algebroid2vector`'s arithmetic directly** (`mot_steenrod.cpp:195-201`): it
computes `position = mon_index[e] + shift` and `coefficient = tauVal(e) * r` (`r` being
the polynomial's own stored coefficient, here `t^-8`). The `t^-8` is NOT a red flag by
itself -- it's designed to cancel against the monomial's own intrinsic `tauVal(x_1^17)`
weight (evidently `+8`) to produce the correctly-normalized `t^0` result. This looks like
deliberate, correct bookkeeping, not an overflow/sign bug -- **could not find a code bug
in this arithmetic**.

**Ruled out three alternative hypotheses for where the "real" bug might be**, all
concretely checked rather than assumed:
1. **Preimage ambiguity for X_1 index 1** (would only matter if `recover_inv_ind`'s
   "first found" choice differs from some other equally-valid preimage): checked via
   `TRACE_QUT_ALL_LVL=0 TRACE_QUT_ALL_TARGET=1` -- **exactly ONE** preimage exists
   (position 2). No ambiguity possible here.
2. **G_0 having multiple cogenerators** (would undermine the "no chain condition needed
   at level 0" argument from earlier in this file): checked via `TRACE_INJ_QUT=3` --
   **G_0 has exactly one cogenerator** (`{0-0}`), confirming the earlier reasoning.
3. **`im_gk` for x=1 silently dropping content**: traced the FULL chain by hand --
   `phi_0(2)=2385` (single term) -> `qut_6(2385)=x1038` (clean singleton) ->
   `inj_7(1038)={7-9}` (PURE, single term, nothing else) -- so
   `prescribed(1)={9:t^0}` is not just correct but has literally nothing else that could
   have been dropped. Airtight.

**Where this leaves things**: every individual computation checked out consistent and
internally correct (the row/pivot logic, the prescribed values feeding in, the algebroid
arithmetic's own tau-bookkeeping). The failure is real (`RHS` is definitively `0` since
`qut_6(2495)=0` is directly confirmed, while `LHS` is definitively `{7-10}` via a fully
traced, self-consistent mechanical computation) but its root cause has NOT been
isolated to a specific bug in this codebase -- it may require either (a) independently
verifying the actual dual-Steenrod-algebra comultiplication value at this specific
bidegree (whether `x_1^17` really should act on cogenerator 9 in a way that lands
exactly on cogenerator 10's own base), which is a genuine math check beyond generic
code-reading, or (b) reconsidering whether pinning `phi_bar` independently "one x at a
time" (the core simplifying claim of the whole target-cofree construction) can actually
break down when two DIFFERENT rows' pivots (like x=1's 974 and x=116's 1058) interact
through a SHARED algebroid pathway in a way the per-x `prescribed` derivation doesn't
account for -- this second possibility, if true, would be a genuine **math-level**
finding requiring discussion before any further code change, not something to guess at
silently.

**Status**: reporting this to the user rather than continuing to guess further, per
CLAUDE.md's directive to stop and discuss math-level uncertainties rather than push
forward silently.

## NEW LEAD (2026-07-12): per-generator degree truncation, found via user's questions

User asked to check primitivity of offset 112 (the local monomial index within G_0's
one generator that starts the p=112 failure). **Confirmed NOT primitive**:
`coaction(112)` [scoped to `F_src.total_rank=1056`=G_0] has **8** terms (positions
0,1,2,3,72,84,97,112, all algebroid coefficients pure powers of `x_1`), not the 2 a
primitive element would have. This matches `Δ(x_1^19) = Σ_{i submask of 19} x_1^i ⊗
x_1^(19-i)` exactly (19 = `10011`₂ has 3 set bits → 8 submasks) -- strong confirmation
`mon_array[112] = x_1^19`.

Only ONE of these 8 terms has a populated `φ̄` (position 0, since G_0's only cogenerator
is at 0): the "`x_1^19 ⊗ 1`" component. This computes
`algebroid2vector(x_1^19 (twisted), shift=base_of_cogen9_in_G_6=2383)` = position
`mon_index[x_1^19] + 2383 = 112 + 2383 = 2495`, which is **exactly cogenerator 10's own
base** in G_6 (confirmed: `2495` is itself a cogenerator, and `2383+112=2495` matches
cogenerator 9's block being exactly 112 positions wide).

**Traced this to `adjoint()` (`hopf_algebroid/9.h:10`)**:
`result.total_rank = ranksBelowDeg(maxDeg - underlyingDeg(generator's own degree))` --
each generator's cofree block is explicitly degree-truncated relative to the GLOBAL
`maxDeg` (=40, from argv[1] `"40 30 4 6 6 9"`). This is a strong candidate explanation:
`x_1^19` may be exactly the first monomial whose degree pushes `deg(x_1^19) +
deg(cogen9)` past `maxDeg=40`, meaning "x_1^19 · cogen9" was never allocated a
legitimate position in G_6 at all -- and `algebroid2vector`'s blind `mon_index[e]+shift`
arithmetic, with no per-generator bounds check, spills into cogenerator 10's
numerically-adjacent (but semantically unrelated) position, because generators are
packed contiguously (`position_of_gens`).

**Important correction found while re-reading `mot_steenrod.cpp:260-301`**: `ranksbelow[d]`
(what `ranksBelowDeg` returns) is a **count** of monomials with `expoDeg(e) <= d`, built
by scanning `mon_array` in file order and incrementing `ranksbelow[d]` for every
monomial of degree `<= d`, for every `d` from 0 to `maxDeg`. This is a per-degree
histogram/count, **not** inherently tied to `mon_array`'s ORDER unless `mon_array` is
itself sorted by ascending degree. The whole "generator's block = positions
`0..total_rank-1` are exactly its in-budget monomials" scheme silently **assumes**
`mon_array` (loaded from an external file, built by the separate `e2p`/`motTab` tools)
is degree-sorted ascending -- not yet independently verified, though the "8 submask
terms, matching Δ(x_1^19)'s exact structure" finding is strong indirect evidence that
`mon_index`/`mon_array` are behaving sensibly at least combinatorially.

**Not yet confirmed numerically**: the actual degree values (`deg(x_1^19)`,
`deg(cogenerator 9 in G_6)`) to directly verify `deg(x_1^19) + deg(cogen9) > maxDeg=40`.
Was about to add a temporary trace hook (extending the existing `TRACE_DESCRIBE_POS`
block in `yoneda2.cpp:513-527`) to print `MotDegree` (via `.degree(pos)`, which has an
`.output()` giving `"(deg,weight)"`) alongside the existing cogenerator/offset
description, to get concrete numbers. **Next action**: add this print, rebuild, and
directly check the degree arithmetic.

**Open framing question from the user, not yet resolved**: even once the degree
inequality is confirmed, there are two very different explanations consistent with it:
(1) the truncated model is self-consistent and this term is *correctly* absent (in which
case `TRACE_SQUARE_CHECK` at `p=112` may be testing outside where the theory promises
exactness), or (2) `maxDeg=40` is simply too small for what `φ`'s chain-map construction
needs to express (a real, if different, bug/limitation). Need more investigation to
distinguish these -- e.g. checking whether the RHS side of the failing check is
ALSO subject to the same kind of degree constraint in a way that's consistent with (1),
or whether there's reason to think φ inherently needs more room than the resolution's
own generators do.

## Degree numbers confirmed the boundary hypothesis exactly

Added a temporary degree print to the `TRACE_DESCRIBE_POS` hook (`yoneda2.cpp`) and got:
`deg(x_1^18)=18`, `deg(x_1^19)=19`, `deg(cogenerator 9 in G_6)=22`, `deg(cogenerator 10
in G_6)=23`. `maxDeg=40` (argv[1]). `18+22=40` (exactly at the truncation boundary,
allowed); `19+22=41` (one past it). This is an exact, not approximate, match to the
"generator 9's block is exactly 112 positions wide" finding -- offset 111 (`x_1^18`) is
the last monomial that fits within cogenerator 9's degree budget, offset 112 (`x_1^19`)
is the first that doesn't.

**Conclusion**: the truncated model is self-consistent; `TRACE_SQUARE_CHECK`'s
exhaustive per-position sweep is bound to eventually test degree-41 data against a
degree-40-truncated resolution, and failing there reflects the sweep asking a question
beyond where the resolution has any data, not a bug in `φ`'s construction. Supporting
evidence: the actual product `{4-6}*{6-9}` still comes out as the clean, correct
`t^0{10-25}` despite this artifact being present.

## TRACE_SQUARE_CHECK hardened to recognize this and skip past it (2026-07-12)

Per user request ("avoid false positives... check higher homological degrees"),
modified `TRACE_SQUARE_CHECK` (`yoneda2.cpp:631+`) to compute `deg_beta` (the beta
class's own degree, via `R[bs].F.degree(beta_rep_pos.dataArray[0].ind).deg`) once, and
for every failure candidate, compute `d_p = R[k].F.degree(p).deg` and skip (continue
scanning, incrementing a `n_boundary_skipped` counter) whenever `d_p + deg_beta >
MOP.maxDeg` -- exactly mirroring the `adjoint()` truncation formula
(`hopf_algebroid/9.h:10`) and the SAME check pattern already used elsewhere in the
codebase (`tao_bockstein.cpp:382`, `F0.generators.degree[i].deg + x_deg <=
ms_oper.maxDeg`). Genuine (non-boundary) failures still stop the scan and print full
detail as before; the summary line reports how many boundary artifacts were skipped.
`TRACE_SQUARE_CHECK_VERBOSE_SKIPS` env var lists each skipped one if set.

**Next action**: rebuild, redeploy, clear caches, and run the improved check with a
large `checkmax` to scan as many homological-degree levels (k=1,2,3,...) as possible,
to check for genuine (non-boundary) bugs at higher levels that were previously
unreachable because the scan got stuck at the k=0 boundary artifact.

## RESULT: clean scan across ALL 24 levels of this beta chain -- no genuine failures

`TRACE_SQUARE_CHECK=24 ./yoneda2 40 30 4 6 6 9` (24 = `maxlev` for `bs=6,
resolution_length=30`, i.e. every level this beta chain has) now runs to completion with
**zero genuine failures** -- `20802` boundary artifacts correctly identified and
skipped, `no failure found up to level 24`. Spot-checked a sample of the skips with
`TRACE_SQUARE_CHECK_VERBOSE_SKIPS=1`: every one genuinely has `d_p + deg_beta > maxDeg`,
and the excess grows sensibly at higher levels (e.g. k=3 shows excess up to 22, i.e.
`deg(p)=40, deg(beta)=22, sum=62` vs `maxDeg=40`) -- consistent with deeper resolution
levels sitting at higher base degrees, not a classifier that's spuriously permissive.

Re-ran the full regression suite after this change (diagnostic-only, shouldn't affect
`phi` itself, but verified anyway): `40 30 1 1` regression clean, `phi_1(tau_1[1-0])`
still `{7-10}`, both product orders still `t^0{10-25}`.

**Conclusion for this session**: the `φ` chain-map construction (the three real bugs
fixed: X_rank scoping, rhs-propagation through Gauss-Jordan, full Gauss-Jordan instead
of forward-only elimination) is now verified correct across the ENTIRE homological
range of this beta chain, not just levels 0-1. `TRACE_SQUARE_CHECK` itself is hardened
to distinguish genuine chain-map bugs from expected degree-truncation-boundary
artifacts going forward, so it can be reused directly (with a large `checkmax`) as a
real correctness check on other beta/product computations without manual boundary
triage each time.

**Cleanup still to do** (deferred, low priority, purely cosmetic): the temporary debug
instrumentation added while chasing bugs 2-4 this session --
`TRACE_ECHELON_ROW`/`TRACE_ECHELON_N` (`lift.h`, `echelon_pivots_tau0`) and
`TRACE_COACTION_POS`/`TRACE_COACTION_SRC_RANK` (`lift.h`, `cofree_adjoint_row_sum`) --
are harmless (env-gated, no-op when unset) and still present. They were useful enough
this session that they may be worth keeping rather than removing; up to the user.

## NEW: exhaustive multi-(bs,bb) sweep per user request finds a DIFFERENT, non-degree-related failure

User asked to run `TRACE_SQUARE_CHECK` exhaustively across many `(bs,bb)` pairs (not
just `6,9` and `1,1`), warning there might be OTHER kinds of bounds-related false
positives near degree 40. Swept `bs=1..8`, several `bb` per `bs`.

**Result so far**: `bs=1` (all 6 valid classes) and `bs=2` (all 15 valid classes) are
ALL clean. Starting at `bs=3`, exactly ONE class out of 18 valid ones fails:
**`bb=2`** -- and `bb=2` ALSO fails at `bs=4,5,6,7,8` (every one tested), while `bb=0,1,
3,4,...,17` at `bs=3` (tested exhaustively) are ALL clean. So this is NOT the
degree-truncation artifact (confirmed: `deg(p)=0, deg(beta)=8` for `bs=3,bb=2` -- tiny,
nowhere near `maxDeg=40`) -- it's something else, and suspiciously specific to index
"2" recurring across several different filtrations.

**The failure itself** (`bs=3,bb=2`, smallest case): fails at the very first test,
`k=0,p=0` (`{0-0}`, G_0's one cogenerator). `qut_0(0)=0` always (independent of beta --
G_0 is shared/beta-independent), so `LHS=phi_1(inj_1(0))=0` trivially. The check
reduces to: **is `v` (phi_0's value at the cogenerator, i.e. the class's OWN chosen
representative) itself a cycle under `qut_bs`** (`d_bs(v)=0`)? Found `RHS =
inj_4(qut_3(v)) = {4-1}` (NONZERO) for `bs=3,bb=2` -- i.e. `v` does NOT die under
`qut_3`, meaning `cyc[3][2]` (as chosen by `tao_bockstein.cpp`'s Bockstein/cycle
tracking, NOT anything in `lift.h`) may not actually be a genuine cycle in the
RESOLUTION's own `ker(qut_bs)` sense. This is a **necessary condition for phi to be a
valid chain map at all** (forced since `d_0(cogenerator_0)=0` unconditionally) -- so if
it fails, `v` itself is a bad choice of representative, not a `lift.h` bug.

**Not yet resolved**: whether this reflects (a) a genuine bug in `tao_bockstein.cpp`'s
cycle selection for this specific class, (b) a misunderstanding on my part of what
`cyc[bs][bb]` is guaranteed to satisfy (maybe it's a cycle w.r.t a DIFFERENT
differential, e.g. the tau-Bockstein `d_tau`, not literally `qut_bs`, and my assumption
that it must also die under `qut_bs` is wrong), or (c) something else. Added
`TRACE_CYC_TERMS`/`TRACE_CYC_BB` (env-gated, `yoneda2.cpp`) to dump `cyc[s][bb]`'s raw
cogenerator-index terms directly, to check whether `bb=2` is structurally different
(e.g. a genuine multi-term Bockstein correction) from the passing classes. **Next
action**: rebuild, compare `cyc[3][2]`'s terms against a passing class like `cyc[3][0]`
or `cyc[3][1]`, and investigate further before concluding anything.

**Result**: `cyc[3][0]={0:t^0}`, `cyc[3][1]={1:t^0}`, `cyc[3][2]={2:t^0}`,
`cyc[3][3]={3:t^0}` -- ALL FOUR are bare, uncorrected singletons (`cyc[3][bb] =
{bb:t^0}` exactly), including the failing `bb=2`! This rules out "multi-term Bockstein
correction vs simple singleton" as the distinguishing factor entirely -- `bb=2` looks
STRUCTURALLY IDENTICAL (a trivial, uncorrected cogenerator) to the passing `bb=0,1,3`.
So whatever makes cogenerator 2 of `X_3` fail `d_bs(v)=0` isn't visible in `cyc`'s own
data at all -- it must be a genuine fact about `qut_3`/`inj_4`'s structure specific to
that one cogenerator, unrelated to how `cyc` built its representative (since `cyc` did
nothing special here -- it just passed the raw cogenerator through unchanged).

**This now looks like either**: (a) cogenerator 2 of `X_3` is genuinely NOT a permanent
cycle in the resolution's own `ker(qut_3)` sense, and `cyc[3]`/`tao_bockstein.cpp`'s
validity criterion (whatever it actually checks -- likely Bockstein/`d_tau` survival)
is NOT the same thing as "is a genuine resolution-level cycle," making this class an
invalid input to `lift_first_step_sum`/`phi_beta`'s WHOLE construction in the first
place (a scope/precondition question about which `(bs,bb)` pairs are even legitimate
to test this way, not a `lift.h` bug), or (b) my own understanding of what property `v`
needs (`d_bs(v)=0`, forced by `d_0(cogenerator_0)=0`) is subtly wrong. Not yet resolved
-- reporting this to the user now rather than digging further into `tao_bockstein.cpp`
(a different subsystem) without checking in first, especially since it may reflect a
real issue there rather than in anything touched this session.

## ROOT CAUSE FOUND for the {3-2}-is-not-a-cycle anomaly -- in tao_bockstein.cpp, confirmed by direct trace

User confirmed independently that `{3-2}` is mathematically NOT a cycle, and asked why
`cyc[3]` still lists it. Traced via the EXISTING `TRACE_TAG_SEM` hook
(`yoneda2.cpp:370`, no new code needed):

```
TRACE_TAG_SEM level s=3 rank=18
  tables[3].tag_index keys (level 2 indices!): -1
  tables[3].cycle_index keys (level 3 indices) with tag value: 0(INVALID) 1(INVALID)
    3(INVALID) 4(INVALID) ... 17(INVALID)          <-- note: "2" is ABSENT here
  tables[4].tag_index keys (level 3 indices -- i.e. level-3 SOURCES feeding level 4): -1 2
```

**Conclusive**: `tables[3].cycle_index` (level 3's own table of entries confirmed as
genuine, uncorrected cycles via `tag==Invalid`) has NO entry for cycle=2 -- level-3
index 2 was never independently confirmed as a cycle there. Meanwhile
`tables[4].tag_index` DOES contain "2" as a key -- level-3 index 2 is registered as the
SOURCE of a real (non-Invalid) tau-Bockstein differential feeding level 4, i.e. it is a
**boundary** (killed/paired element), not a survivor.

**The actual bug is in `tau_table::get_cycles()` (`tao_bockstein.cpp:231-245`)**:
```cpp
cycle_data tau_table::get_cycles() const{
    cycle_data res;
    for(auto &itm : table){
        if(itm.tag == tau_table_entry::Invalid)
            res.emplace(itm.cycle,itm.full_cycle);
        else{
            //modify the boundary by dividing powers of tau
            tauPoly fc = tau_oper.power_tau(-itm.diff_length);
            auto nt = tau_module_oper.scalor_mult(fc, itm.full_cycle);
            res.emplace(itm.cycle, nt);      // <-- BOUNDARIES included here too!
        }
    }
    return res;
}
```
The `else` branch's own comment says "modify **the boundary**" -- it knows
`itm.tag != Invalid` means this entry is a boundary/killed element, yet it STILL
inserts it into the returned `cycle_data` map (with a tau-power correction),
indistinguishable from a genuine `tag==Invalid` survivor. `make_cycle_tables`
(`tao_bockstein.cpp:267-276`) then merges `tab[i].get_tags()` (the level-(i-1)-indexed
SOURCES of these same pairings) into `result[i-1]` via `std::map::insert(range)`, which
**does not overwrite existing keys** -- so even where "this index is a differential
source, hence non-surviving" information exists separately, it silently loses to
whatever `get_cycles()` already inserted first, with no reconciliation.

**Net effect**: `cyc[bs]`'s keys conflate two structurally different things --
genuinely-surviving, uncorrected cycles, and boundary/killed elements that happen to
get a "corrected" representative computed (presumably useful elsewhere for REDUCTION
purposes, e.g. `find_cycle`/`find_cycle_sum`, to rewrite arbitrary vectors into normal
form) -- with no marker in `cyc[bs]` itself to distinguish which is which.
`yoneda2.cpp:546`'s `cyc[bs].count(bb)` check (deciding whether `bb` is a legitimate
class to test) can't tell the difference, so it happily picks up boundaries like
`{3-2}` as if they were real classes.

**This is a genuine bug, in `tao_bockstein.cpp`, not in anything fixed this session
(`lift.h`/`yoneda2.cpp`'s phi-construction code)** -- reported to the user rather than
fixed unilaterally, since it's a different subsystem and the right fix (exclude
boundaries from `get_cycles()`'s result entirely? mark them so callers can
distinguish? something else?) depends on how `cyc[]`'s other consumers
(`find_cycle_sum` etc.) rely on the current mixed behavior.

## Refined verdict: NOT a bug in get_cycles() itself; a narrower gap in yoneda2.cpp's usage -- SAFE to ignore if you already know which classes are real

User asked "is this a real bug or can I safely ignore it." On reflection: `get_cycles()`
including boundary entries (with a tau-power correction) is very likely INTENTIONAL and
CORRECT for its OTHER consumer, `find_cycle`/`find_cycle_sum` (`tao_bockstein.cpp:279-`,
used by `yoneda2.cpp`'s product-reduction code at lines ~842,866-868): that's a
Gaussian-elimination-style rewrite algorithm (subtract `table.at(leading_term)` 
repeatedly until the remainder vanishes) which NEEDS an entry for every basis index
that could appear as a leading term during reduction -- including boundaries -- so it
can correctly cancel them via substitution and end up with only genuine survivors in
the final result. Removing the `else` branch would likely break THAT algorithm.

So the real, narrow issue is only in `yoneda2.cpp:546`'s specific usage:
`cyc[bs].count(bb)` is used to decide "is `bb` a legitimate class the user can query
directly," but it can't distinguish a genuine survivor from a boundary-with-correction
-- both are present as keys. This only matters if someone queries `./yoneda2 <bs> <bb>`
for a `bb` that's ACTUALLY a boundary index without realizing it; the tool won't warn
you, it'll just silently compute a non-representative answer for a non-class.

**Verified all three classes already tested this session are genuine, not boundaries**
(via `TRACE_TAG_SEM`, checking `tag=INVALID/genuine-cycle` in each's own
`cycle_index`): `{4-6}` genuine, `{6-9}` genuine, `{1-1}` genuine. So none of this
session's verification work is affected.

**Verdict given to the user**: safe to ignore in practice, AS LONG AS you (or whoever
picks `(bs,bb)`) already know from an independent, correct source (e.g. an Ext chart)
which classes are genuine survivors, and don't query boundary indices. It is a real,
if narrow, usability gap (`yoneda2.cpp` doesn't proactively validate/warn), not a
correctness bug in the resolution or in `get_cycles()`'s primary purpose. Not fixed
(no action requested beyond the diagnosis).

## Exhaustive sweep across all genuine (bs,bb) pairs launched (2026-07-12)

Per user request: run `TRACE_SQUARE_CHECK` exhaustively across many `(bs,bb)`, skipping
the boundary/non-genuine classes found above (via `TRACE_TAG_SEM`'s
`tag=INVALID/genuine-cycle` marker in each level's own `cycle_index`).

Sized the sweep first (cheap `TRACE_TAG_SEM` calls, no `TRACE_SQUARE_CHECK`): genuine
class counts per `bs=0..29` range from 1 (high `bs`, e.g. 21-29) up to 17 (bs=3,5),
totaling **161 genuine `(bs,bb)` pairs**. At ~3-5s/run (each reloads the full 30-level
resolution), this is ~10-15 min total -- decided to run it myself in the background
rather than handing it to the user, per the "if long-running, tell the user how to run
it" framing (161 pairs is long but still very much automatable).

Script: `/tmp/.../scratchpad/sqcheck/full_sweep.sh` -- for each `bs=0..29`: get genuine
`bb` list via `TRACE_TAG_SEM=$bs`, grep for `\d+(?=\(tag=INVALID)`; for each genuine
`bb`, run `TRACE_SQUARE_CHECK=$((30-bs))` (full level coverage for that beta chain) and
classify as OK / GENUINE FAILURE / UNEXPECTED. Launched via Bash `run_in_background`,
job id `bactdb4y7`, log at
`/tmp/.../scratchpad/sqcheck/full_sweep_results.log`.

Hit and fixed one shell portability snag while writing the script: `readarray`/
`mapfile` (bash 4+ builtins) aren't available in this environment's shell (zsh, or a
restricted bash) -- switched to a POSIX `while IFS= read -r bb; do ... done <
tempfile` pattern instead, which works reliably. Also added `</dev/null` to every
`./yoneda2` invocation in the loop as a defensive measure against accidental stdin
consumption inside loops.

**RESULT: sweep complete, 161/161 genuine `(bs,bb)` pairs clean, ZERO failures.**
`grep -c "^OK " full_sweep_results.log` = 161 (matches the sized total exactly, so no
pairs were skipped/dropped by the script); `grep "GENUINE FAILURE\|UNEXPECTED"` = empty.
Covered `bs=0..29`, every genuine class at each, each checked across its FULL level
range (`TRACE_SQUARE_CHECK=30-bs`). This is a much broader, more conclusive
confirmation than the single `{6-9}` chain checked earlier in the session -- the `φ`
chain-map construction (the 3 bugs fixed this session) holds exactly across the entire
genuine-class space of this resolution, not just one hand-picked example.

## Session closed out (2026-07-12)

`STATUS.md` rewritten to reflect the final, conclusive state: plan implemented, three
real bugs fixed, `TRACE_SQUARE_CHECK` hardened against degree-truncation false
positives and passing clean across all 24 levels of the `{6-9}` chain, full regression
suite green. Nothing outstanding for this plan. Task list was empty at session end (no
open TaskCreate items to reconcile).

**This retracts the earlier "phi_0 preimage-independence" math-level hypothesis** --
that was speculation before finding this much more concrete, mechanistically-verified
explanation; no evidence currently supports the phi_0-gap theory, and this bug alone
appears sufficient to explain the observed failure. Should still double check after the
fix that no OTHER failures crop up further along (same pattern as p=6 -> p=14 -> p=27:
each fix has so far revealed the SAME kind of bug one level further in, so it's worth
bracing for a possible p=NN four should this fix only address one instance of the
missing-rhs-propagation class of bug -- though this fix, unlike the previous two, is a
FULL general fix (applies to every row uniformly during the forward pass), so it's
plausible this is the last one).

**Hypothesis (now upgraded to "understood"):** the lift used the wrong universal property
(maps *out of* the cofree source, `φ(a[y])=a·lg(y)`, which is underdetermined). Fix: use
the **target's** cofree property — represent `φ_k` by `φ̄_k=(ε⊗1)φ_k: G_k→Z_{k+bs}`,
reconstruct `φ_k=(1⊗φ̄_k)ρ`. Comodule-map condition becomes free; only the counit-projected
chain condition needs imposing (full square then follows via `E=φ_k∘d−d∘φ_{k-1}` being a
comodule map into a cofree target with `(ε⊗1)E=0 ⇒ E=0`). `φ̄_k` forced on `image(inj_k)`
(=`(ε⊗1)(im_gk)`), free (=0) elsewhere; read off at τ⁰-pivots of an echelon-reduced `inj_k`.

**Ruled IN this session (evidence):**
- Leading terms / pivots are gated to **τ⁰ only** — `symplify_to_led` `matrices/9.h:62`
  (dense) & `:95` (sparse) accept a leading term iff `invertible`, and
  `tauOper::invertible` is τ⁰-only (`mot_steenrod.cpp:44`). So `unify` (`matrices_mem/4.h:3`)
  divides only by τ⁰ ⇒ **no negative-τ powers enter `inj`/`qut`** ⇒ the new F2[τ]-free
  basis change is sound. (Answers the user's negative-τ worry.)
- `im_gk` (`lift.h:179-182`) is exactly the prescribed chain value; only its counit
  projection (offset-0 terms) is needed for `φ̄_k`.
- `cofree_adjoint_row_sum` (`lift.h:78`) already IS the cofree-up formula; with pure-
  cogenerator `φ̄_k` values the `vector2algebroid` raw-position path is unnecessary.

**Next concrete action:** implement the plan — rewrite `lift_one_step_sum` /
`lift_first_step_sum`, simplify `cofree_adjoint_row_sum`, add τ⁰-pivot assert, delete the
correction machinery; then verify `TRACE_SQUARE_CHECK=30 ./yoneda2 40 30 4 6 6 9` is clean
from level 0, `φ_1(τ_1[1-0])={7-10}`, product directions agree, `40 30 1 1` regression
unchanged (rebuild + `\cp -f` binary + `\rm -f *_yoneda2_*_phi*` first).

---

### (Superseded) prior "current status" — the `vector2algebroid` fix cycle

The prior active plan's corrected design (using
`vector2algebroid`/`algebroid2vector` in `cofree_adjoint_row_sum`, `lift.h`) was
implemented and verified, but did NOT resolve the core disagreement — it was necessary but
not sufficient, and is now subsumed by the target-cofree construction above:
- Build is clean (`cmake --build build --target yoneda2`).
- **Gotcha hit during verification**: `./yoneda2` in the repo root was a *stale* binary
  from an earlier build; `cp`/`rm` are aliased/wrapped (`cp='/bin/cp -ia'`, `rm` is a
  shell-function, both print interactive-looking prompts) in this shell — use `\cp -f` /
  `\rm -f` to bypass and actually overwrite/delete. First round of "verification" runs
  silently used the old binary and old phi caches and looked like nothing had changed;
  re-ran after `\cp -f build/yoneda2 ./yoneda2` and `\rm -f *_yoneda2_*_phi*` and got a
  real result.
- Regression check `./yoneda2 40 30 1 1` unchanged (matches prior clean single-term
  output).
- `TRACE_LIFT=100000 ./yoneda2 40 30 4 6 6 9` (fresh phi caches) confirms `phi_beta[3]`'s
  cogenerator `j=8` (`{3-8}`, `x=2025`, `cog_pos=3583`) now gets the correct genuine
  5-term nonzero `combined` value (`1252^[t^1] 1345^[t^0] 1386^[t^1] 1387^[t^0]
  1476^[t^1]`) instead of the previously-dropped `0`. Confirms the plan's fix is
  correctly implemented.

**However, the core disagreement is NOT resolved**: `./yoneda2 40 30 4 6 6 9` still gives
`{4-6} * {6-9} = 0`, while `./yoneda2 40 30 6 9 4 6` still gives `{6-9} * {4-6} =
t^0{10-25}` (fresh caches both times, confirmed via `\rm -f *_yoneda2_*_phi*` before each
run). So the plan's fix was real and necessary but not sufficient — there is at least one
more bug.

## Current thread: re-confirmed and sharpened the "different preimages give different
## downstream phi" lead from `debug_notes_archive4.md` (still open after today's fix)

Re-ran the exact `2442` vs `2444` check from archive4 (X-index 1416 behind `{4-6}`, level
`G_3`, beta=`{6-9}`) with **today's `vector2algebroid` fix in place** (fresh phi caches,
via `TRACE_PHI_AT_LVL=3 TRACE_PHI_AT_POS=<pos>`):
- `phi_beta[3].find(2442) = {1552:t^0, 1568:t^0}` — **byte-identical** to before today's
  fix (expected: this value was already "pure"/base-position-only, so the
  `vector2algebroid` fix never touched it).
- `phi_beta[3].find(2444) = {1554:t^0, 1570:t^0}` — also unchanged.

**Directly verified the exactness-based invariant is violated** (per user's framing:
`2442 - 2444 ∈ ker(qut_3) = image(inj_3)` since both are confirmed clean `tau^0`
singleton preimages of the same X-index 1416, so `qut_9(phi_beta[3](2442) -
phi_beta[3](2444))` MUST be `0` for `phi_beta[3]` to be a genuine chain map — checked via
`TRACE_QUT_LVL=9 TRACE_QUT_POS=<p>` on all 4 relevant positions):
- `qut_9(1552) = t^1·x663`, `qut_9(1568) = t^1·x663` — equal, cancel.
- `qut_9(1554) = t^1·x664 + t^0·x724`, `qut_9(1570) = t^1·x664 + t^0·x725` — the `x664`
  terms cancel between these two, but `x724`/`x725` do NOT cancel against anything.
- Net: `qut_9(phi_beta[3](2442) - phi_beta[3](2444)) = t^0·x724 + t^0·x725 ≠ 0`.
  **Confirmed nonzero — the chain-map identity is genuinely violated**, not a
  representation artifact.

**User's exact diagnostic request (not yet done)**: trace precisely WHERE this breaks,
using the 3-part chain: (1) is `2442 - 2444` really in `image(inj_3)`? (2) does
`phi_beta[3]` actually send `image(inj_3)` elements correctly (i.e. does `phi_3 ∘ inj_3`
match what the chain condition requires)? (3) where does `qut_9(...)` end up nonzero
instead of 0 — i.e. which of (1)/(2) fails. Next concrete step: find the explicit
`z ∈ X_3` with `inj_3(z) = e_2442 + e_2444` (mod 2) by brute-force scanning `inj_3`'s rows
(no existing trace tool does this reverse lookup — need to add one, e.g. a temporary
`TRACE_FIND_INJ_TARGETS` env var near the existing `TRACE_INJ_AT_X` block,
yoneda2.cpp:429-446), then check whether `phi_beta[3](inj_3(z))` matches what
`lift_one_step_sum`'s formula would require it to be, to localize the bug to either (2)
`cofree_adjoint_row_sum`'s Part-2 extension formula at these specific offsets, or
somewhere in how `inj_3`/`z` interact with the level-2→3 construction.

**Ground rule for this thread (from user)**: investigate freely and add tracing without
asking permission; only check in once something interesting is found or a new fix idea
emerges — do NOT implement a fix without checking in first (per the CLAUDE.md math-level-
change guardrail, now explicitly confirmed to still apply here even though tracing itself
is pre-approved).

## Found z with inj_3(z) touching {2442,2444}; confirmed Part-2 extension is mechanically
## exact; found the SAME preimage-ambiguity failure recurs ONE LEVEL UP too (level m=2)

Added `TRACE_FIND_INJ_LVL`/`TRACE_FIND_INJ_TARGETS` (yoneda2.cpp, brute-force reverse
lookup: which `z` have an `inj_lvl(z)` row touching given target positions, flagging
whether the target is the row's OWN leading/smallest-index term). Used it to hunt for
`z` with `inj_3(z) = e_2442+e_2444`; got partial matches (`z=1017`, `z=1018` etc.) but a
full manual Gaussian-elimination reconstruction by hand turned out to be a rabbit hole
(residual terms like `{3-6}+1`/pos3187 that are genuinely in `ker(qut_3)` per direct
check but never anyone's own leading term in a first-pass sweep) — **abandoned this
specific manual-reconstruction approach as too unreliable/slow; did not find explicit
z**. This code (`TRACE_FIND_INJ_TARGETS`) is still in the tree (yoneda2.cpp ~line
448-480) and works, just wasn't the most direct path to an answer.

**Confirmed Part 2 (`cofree_adjoint_row_sum`'s coaction-based extension) is mechanically
exact, not buggy**: `lg[{3-4}]` (cog_pos=2429) `= {1539:t^0, 1555:t^0}` (call it `v`).
Position 2442=`{3-4}+13` and 2444=`{3-4}+15` are exactly `x₁₃·v` and `x₁₅·v` (local
monomial at that offset acting on `v`): `1539+13=1552`, `1555+13=1568` ✓ matches observed
`φ(2442)`; `1539+15=1554`, `1555+15=1570` ✓ matches observed `φ(2444)`. So both φ values
are exactly what the forced cofree-adjoint formula requires given `v` — no bug in this
step.

**Found the SAME kind of preimage-ambiguity failure ONE LEVEL UP, at level m=2** (i.e.
building `phi_3` from `G_2`, for `x=1006`, `{3-4}`'s own `X_3` index): `qut_2` has FOUR
clean singleton preimages of `x=1006`: `912, 1663, 1664, 2351`. `recover_inv_ind` picks
`912`. `phi_beta[2].find(912) = {1619:t^0}` (nonzero) but `phi_beta[2].find(1663) =
find(1664) = find(2351) = 0` (all empty) — i.e. picking any of the OTHER three instead
would make `lg[{3-4}]` come out `0` instead of `{1539,1555}`, an even starker
disagreement (zero vs. nonzero) than the original 2442-vs-2444 case. Directly verified
the exactness invariant is violated here too: `qut_8(phi_beta[2](912) -
phi_beta[2](1663)) = qut_8(1619) - 0 = {834:t^0,844:t^0} ≠ 0`.

**Conclusion so far**: this is not a one-off; the same failure (multiple valid clean
qut-preimages giving genuinely different downstream phi, violating exactness) recurs at
multiple levels (confirmed at m=2, building phi_3; original finding at the m=3 level
building phi_4). Reported this to the user with full evidence; user's response: this is
about whether "phi ∘ (inj∘qut)" commutes with "(inj∘qut) ∘ phi" — i.e. the FULL chain-map
square, not just the ker(qut)-restricted weaker form — and asked to **start checking
systematically from level (resolution degree) 0 upward and find the FIRST level where
this fails**, rather than continuing to spot-check higher levels first.

**Just implemented (not yet run/verified)**: `TRACE_CHAIN_CHECK=<maxlev>` (yoneda2.cpp,
right after the phi_beta build loop, ~line 611-663). For each level `m=0..min(maxlev-1,
checkmax)` in order, for each cogenerator `x` of `X_{m+1}` (via `gens[m+1]`), sweeps
`qut_m` for ALL clean (`tau^0`) singleton preimages of `x`, and if there are 2+
candidates, computes `phi_beta[m]` at each, pushes through `qut_{m+bs}`
(`apply_ring_matrix_to_sum`), and checks whether all candidates agree (pairwise diff via
`tauPolySum_module_oper.add`, checking for exactly zero). Stops and reports full detail
at the FIRST (smallest `m`) failure found. Added `#include <map>` (was missing,
needed for `candidates_by_x`).

**Next step**: build, copy binary, clear phi_beta caches for the `6_9` beta (needed for
`4 6 6 9`), run `TRACE_CHAIN_CHECK=<maxlev> ./yoneda2 40 30 4 6 6 9` and report the FIRST
failing level found (expect it to be at or before `m=2`, confirming or sharpening the
level-2 finding above — but let the tool find the true first failure rather than
assuming).

## Ran TRACE_CHAIN_CHECK; found+diagnosed a clean minimal counterexample at level m=2

`TRACE_CHAIN_CHECK`'s systematic sweep (level 0 upward) confirmed level 0 and 1 are
clean; FIRST failure at level m=2, cogenerator `{3-2}` (X_3 index x=776, candidate
preimages in G_2: 9,10,903,1657; `recover_inv_ind` picks 9). Hand-traced the FULL
arithmetic for both preimage 9 (currently used) and preimage 903 (alternate): preimage 9
gives `main_term=0` (since `lg[{2-0}]=0`) so `combined = -correction = {1295:t^1}` (a
junk, non-cogenerator position). Preimage 903 gives `main_term={1295:t^1,1418:t^0}`
(since `lg[{2-1}]={1604:t^0}`, `qut_8(1610)={776:t^0}` as an X_9 index — coincidental
numeral reuse, unrelated to the X_3 index of the same value —, `inj_9(776)=t^1*pos1295 +
t^0*pos1418`), so `combined = {1418:t^0}` = a CLEAN cogenerator hit (`{9-12}`). Reported
this to the user as the first failure found via preimage-ambiguity sweep.

**User clarified this was NOT what they asked for**: they want a literal, direct test of
the commutative square `phi_{k+1} . (inj_{k+1} . qut_k) == (inj_{k+bs+1} . qut_{k+bs}) .
phi_k` as maps `G_k -> G_{k+bs+1}`, evaluated at EVERY position of `G_k` (not just
cogenerators, not just via alternate-preimage comparison), scanning levels `k=0,1,2,...`
in order and stopping at the first failing `(k,p)`. This is a materially different,
more direct/exhaustive check than `TRACE_CHAIN_CHECK`.

**Just implemented (not yet run)**: `TRACE_SQUARE_CHECK=<maxlev>` (yoneda2.cpp, inserted
right before the `TRACE_CHAIN_CHECK` block). For each level `k=0..min(maxlev-1,
checkmax)` and position `p=0..R[k].F.total_rank-1` (in order), computes:
- `LHS = phi_beta[k+1]` applied to `inj_{k+1}(qut_k(p))` (via a new local `applyPhi`
  lambda that composes a `tauPolySum`-valued matrix with a `tauPolySum`-valued vector,
  since `apply_ring_matrix_to_sum` in lift.h only handles `tauPoly`-valued matrices).
- `RHS = inj_{k+bs+1}(qut_{k+bs}(phi_beta[k].find(p)))` (reusing `apply_ring_matrix_to_sum`
  twice, same matrices `lift_one_step_sum` itself uses for `im_ck`/`im_gk`).
Stops and dumps full detail (`qut_k(p)`, LHS, `phi_k(p)`, RHS, DIFF, all with
`{lvl-cogen}`/`{lvl-owner}+offset` naming) at the first nonzero `LHS-RHS`. Build is clean
(verified). **Not yet run** — next action: copy binary, clear `6_9` phi caches, run
`TRACE_SQUARE_CHECK=30 ./yoneda2 40 30 4 6 6 9` and report the first failure (or
"no failure found up to level N" if none, which would itself be an important, surprising
data point contradicting the current active-bug hypothesis).

## Instrumentation currently in the tree (all env-gated, harmless when unset)

- `yoneda2.cpp`: `TRACE_QUT_SWEEP_LEVEL/_COG/_TARGET_X`, `TRACE_MON_NAMES`,
  `TRACE_INJ_QUT`, `TRACE_INV_INJ`, `TRACE_TAG_SEM`, `TRACE_QUT_POS/_LVL`,
  `TRACE_QUT_ALL_LVL/_TARGET` (full sweep of every valid preimage of a target, not just
  the first), `TRACE_PHI_DUMP`, `TRACE_PHI_TABLE`, `TRACE_PROD_S`/`_A`,
  `TRACE_DESCRIBE_LVL`/`_POS`, `TRACE_PHI_AT_LVL`/`_POS`, `TRACE_INJ_AT_LVL`/`_X`,
  `TRACE_EXACTNESS_LVL` (checks `qut(inj(x))==0` for every basis element at a level, plus
  a dimension-count sanity print).
- `lift.h`: `TRACE_LIFT` (per-cogenerator trace budget in `lift_one_step_sum`), plus
  `TRACE_LIFT_SUM_PRE`/`_CORR` sub-prints and a per-cogenerator summary line, all gated
  on the same budget.
- `tao_bockstein.cpp`: `TRACE_FIND_CYCLE_SUM` (traces `find_cycle_sum`'s reduction steps).

## Key facts established

- `tauPoly` = single tau-monomial (`int16_t` exponent); `tauPolySum` = genuine `F2[tau]`
  polynomial (`polynomial<F2>`), used only for `phi_beta`/`M_beta`'s accumulator — the
  resolution's own `inj`/`qut`/differential data stays `tauPoly`-valued throughout.
- Exactness of the resolution (`ker(qut)=image(inj)`) has been independently re-verified
  (sweep + dimension counting) at every level checked — treat it as solid ground.
- `lg`/`phi`'s values must be raw target-comodule positions, not projected down to
  "combinations of target cogenerators only" — see `CLAUDE.md` for the general lesson.
- `.coaction()` computes the full comultiplication (can be a large sum) — not a simple
  "which monomial, which generator" lookup. Use `vector2algebroid` for that instead.
- Always verify empirically with concrete printed values before concluding a root cause.

## Full math write-up: PHI_EXTENSION_ISSUE.md (this session's main output)

The root cause understanding from this session (why `φ_k`'s extension from cogenerators
to the rest of `G_k` is underdetermined, with a fully-worked counterexample
`φ_1(τ_1[1-0])` needing to be `{7-10}` not `0`) is written up in full at
`PHI_EXTENSION_ISSUE.md` (repo root) — read that file first, it's self-contained.
**Important correction already folded into that file**: an earlier draft wrongly claimed
the comodule-map condition alone uniquely determines `φ_k` from cogenerator values (cited
injectivity of `Δ_target` via the counit). User caught this: injectivity doesn't imply the
equation `Δ_target(w)=1⊗w` has a unique solution — *every* cogenerator combination `w`
satisfies it trivially (`Δ_A(1)=1⊗1`). So the comodule condition alone is genuinely
underdetermined; only the chain condition (via a diagram chase using a *different*,
higher-degree source element) picks out the right answer. Don't re-derive this by hand
again without reading the corrected writeup — the derivation is subtle and easy to
re-get-wrong the same way.

## Checked yoneda.cpp (one-time exception, see top of file)

User asked to check whether `yoneda.cpp` (working tree AND git `HEAD`, which use an
identical `TauEchelon`/`invInj`/`secQut` algorithm — HEAD's version is just inlined where
the working tree extracted it into a `lift_step` lambda) gets `φ_1(τ_1[1-0])` right.
**It does not** — same wrong answer (`0`) as `yoneda2.cpp`. Mechanism (traced via new
`TRACE_LIFT_STEP_LVL`/`_X` and `TRACE_PSI_AT_LVL`/`_POS` hooks added to `yoneda.cpp`):
position `x=4` (`{1-0}+4`) is neither in `image(inj_1)` (per the `TauEchelon` solver) nor
itself a cogenerator, so `yoneda.cpp` falls into its final fallback branch and **explicitly
forces it to `0`** with no attempt at any comodule/coaction-based extension at all —
strictly less sophisticated than `yoneda2.cpp`'s (still-insufficient) attempt. Confirms
this is a genuinely open algorithmic problem, not something already solved in the older
code.

**Build note**: `yoneda` target's CMakeLists.txt was missing `Fp.cpp` (link failure,
`Fp_Op::unit`/`Fp_Op::Fp_Op` undefined) — this is the exact pre-existing gap flagged in
`CLAUDE.md` pitfall #10 (tauPolySum's `Fp_Op` dependency was added to `mot_steenrod.cpp`
without updating the `yoneda` target). Fixed in `CMakeLists.txt` (added `Fp.cpp` to the
`yoneda` target's sources) since it was needed just to compile and check this. HEAD's
version was checked via a separate manual build (`g++ -std=c++11 -DEXPONENT_WIDTH=32`,
not through CMake) with the same trace hooks patched in, in `/tmp/yoneda_head_check/` —
not part of the repo, just a scratch comparison copy.
