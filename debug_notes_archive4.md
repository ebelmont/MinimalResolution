# Debugging Notes: Yoneda Product Bugs in yoneda2.cpp (fresh cycle, 2026-07-10)

Full history through this point is archived in `debug_notes_archive1.md`,
`debug_notes_archive2.md`, `debug_notes_archive3.md` (chronological). This file is the
current, active summary — keep it short; append detail only as it becomes load-bearing.

**Hard constraint, unchanged across all sessions**: do not read, edit, or run `yoneda.cpp`
(no "2") for any reason. Only `yoneda2.cpp` and its supporting headers are in scope.

## Bugs fixed so far this cycle

1. **`tag_index` level mismatch** (yoneda2.cpp:567, confirmed real): `tables[s].tag_index`
   holds level-`(s-1)` indices, not level-`s`. Fixed to check `tables[s+1].tag_index`
   instead (with an `s+1 < tables.size()` bounds guard).
2. **`recover_inv_ind`'s diagnostic hardened to a hard abort** (lift.h:79-93), per
   standing user request from a prior session.
3. **`tauPolySum` genuine multi-term F2[tau] type** introduced for `phi_beta`'s
   accumulator (lift.h's `_sum` function variants, `mot_steenrod.h/.cpp`,
   `tao_bockstein.h/.cpp`'s `find_cycle_sum`, `yoneda2.cpp`'s `phi_beta`/`M_beta`/`fmt`).
   This fixed the `INVARIANT BROKEN: add(t^2, t^1)` hard abort on
   `./yoneda2 40 30 4 6` — the crash is gone, `1 1` regression-checks clean.
   Side effect: had to add `Fp.cpp` to several CMake targets that share `mot_steenrod.cpp`
   (`motTab`, `mot_mult`, `dump_gens`, `test_lift`, `tauBoc`) — NOT done for the `yoneda`
   target (out of scope); that target's build is likely broken until fixed or the user
   says otherwise.

## Bug still open (current focus)

`{4-6}*{6-9}` computes as `0`; `{6-9}*{4-6}` computes as `t^0{10-25}`. User has confirmed
the nonzero value is correct and the zero is the bug. These must be equal (graded
commutativity, no signs mod 2).

Repro:
```
./build/yoneda2 40 30 4 6 6 9   # => {4-6} * {6-9}  =  0        (BUG)
./build/yoneda2 40 30 6 9 4 6   # => {6-9} * {4-6}  =  t^0{10-25}   (correct)
```

User feedback: don't expect a "diff the two directions" strategy to be fruitful, since
they're structurally different computations (different beta, different chain map
levels/target). Just trace the failing direction thoroughly and follow whatever's found.

## Progress tracing `{4-6}*{6-9}=0` (this cycle)

**Found and fixed a real gap**: `lift_one_step_sum` (lift.h) was written without the
`TRACE_LIFT` prints the original `lift_one_step` had — since `yoneda2.cpp` now calls the
`_sum` variant exclusively, `TRACE_LIFT` was silently producing zero output. Added
equivalent tracing as `TRACE_LIFT_SUM_PRE`/`_CORR` and a per-cogenerator summary line
(still gated on `TRACE_LIFT`, lift.h). Also added a small new `TRACE_DESCRIBE_LVL`/
`TRACE_DESCRIBE_POS` env var (yoneda2.cpp) to name an arbitrary raw position as
`{lvl-cogen}` or `{lvl-owner}+offset`, reusing the same describe-pos logic already used
elsewhere.

**IMPORTANT gotcha hit and corrected**: `cog_pos=N` values in the trace are per-LEVEL
position numbers (each `G_k` renumbers positions from 0), so grepping for a raw number
like `cog_pos=912` across the WHOLE trace log matches unrelated levels too — always
narrow to the specific level's line range first (bounded by consecutive `j=0` trace
lines) before trusting a match.

**Trace chain for the failing computation** (`bs=6,bb=9`, i.e. beta=`{6-9}`; alpha=
`{4-6}`, i.e. `s=4,a=6`; single-product command:
`./build/yoneda2 40 30 4 6 6 9`, env `TRACE_PROD_S=4 TRACE_PROD_A=6 TRACE_LIFT=100000`):

1. `cyc[4].at(6)` (alpha's cycle rep) = pure `cog6^t0` (single term).
2. `phi_beta[4].find(pos=2711)` (cog6's own position) = **empty** — this is the direct
   cause of the final `0`.
3. Traced `lift_one_step_sum`'s build of this exact cogenerator (level 4, `j=6, x=1416,
   cog_pos=2711`): `im_prev = phi_beta[3].find(g_prev_pos=2442) = {1552:t^0, 1568:t^0}`
   (nonzero, 2 terms) — but `qut_9.find(1552) = qut_9.find(1568) = {663: t^1}`
   (**identical** target+coefficient for both), so `apply_ring_matrix_to_sum` correctly
   sums two copies of the same tau^1 contribution to x663 → **cancels to exactly 0 mod
   2** (`im_ck` empty). `inj_row` for this cogenerator is pure (self-term only, no
   correction), so `combined` is also empty. This is a genuine mod-2 cancellation, not a
   representation/truncation bug — the mechanism itself is mathematically sound *given*
   `im_prev`'s two terms.
4. Traced where `im_prev`'s two terms come from: position 1552 = `{9-17}+13`, position
   1568 = `{9-18}+13` (SAME local offset 13 in two different cogenerators' summands).
   These come from `phi_beta[3].find(cog4's own position 2429) = main_term_cogens =
   {17:t^0, 18:t^0}` — i.e. `phi_3({3-4}) = cog17 + cog18` (a genuine 2-term value AT
   THE COGENERATOR LEVEL, level 3, under beta=`{6-9}`), extended to offset 13 via
   `cofree_adjoint_row_sum`'s coaction-based formula (same coaction term applied to both
   summands, naturally landing at the same local offset in each).
5. Traced level 3's build of `{3-4}`: `im_prev = phi_beta[2].find(g_prev_pos=912) =
   {1619:t^0}` (single, clean term) → `qut_8.find(1619) = {834:t^0, 844:t^0}` (a
   genuinely 2-term row — position 1619 must be a PIVOT position in `qut_8`'s
   Gaussian-elimination-based construction, which CAN legitimately have multiple terms,
   per `make_quotient`'s documented invariant that only non-pivot positions are
   guaranteed pure singletons). `inj_9` then sends `{834,844}` to `{1539:t^0,1555:t^0}` =
   cog17/cog18's own positions. `inj_row` here is also pure, so no correction;
   `main_term_cogens = combined = {17:t^0,18:t^0}` directly.
6. Traced position 912 (level 2, owner `{2-1}+15`): `phi_beta[2]`'s cogenerator value
   `lg[{2-1}]` (cog_pos=897) = pure `cog9^t0` (single, clean). Position 912's value is
   then a deterministic coaction-based extension of this single clean value — nothing
   suspicious found here; this looks like solid ground.
7. **Hypothesis considered and RETRACTED**: initially thought cogenerators 17 and 18 (of
   `G_9`) might be a direct tau-Bockstein pair (one tagging the other), which would make
   their sum-then-cancel pattern "expected" via the Bockstein reduction rather than a raw
   qut collision. Checked `tables[9].cycle_index` via `TRACE_TAG_SEM=9`: entry
   `17(tag=18)` exists, BUT per this session's own earlier-confirmed finding,
   `tables[9].tag_index`'s keys (and hence the `.tag` field read off table rows) are
   **level-8** indices, not level-9 — so "tag=18" here almost certainly refers to an
   unrelated level-8 object, not level-9's cogenerator 18. This looked like a tempting
   lead but is very likely a numeric coincidence, not a real relationship. **Do not reuse
   this reasoning without re-deriving it properly** (e.g. by checking what level-8 object
   tag 18 actually names, and whether IT has any bearing on cogen 18 of level 9 — probably
   not, but not independently re-verified since abandoning this thread).

## Current state / open questions

- The raw resolution data (`qut_9`, `qut_8`, `inj_9`) is being treated as ground truth
  (out of scope per the earlier "phi_beta only" scope decision) — everything traced so
  far is consistent with `lift_one_step_sum`/`cofree_adjoint_row_sum` correctly
  implementing the intended formula on that data. The cancellation at step 3 looks like a
  mathematically real consequence of the resolution's actual `qut_9` values, not an
  artifact of the tauPolySum machinery.
- **Not yet resolved**: if the mechanism is sound and the resolution data is ground truth,
  why does the final answer come out wrong? Candidate remaining suspects, none yet
  checked: (a) `cyc[4].at(6)`'s cycle representative (`pure cog6^t0`) might itself be an
  incomplete/wrong representative for class `{4-6}` — i.e. the classical (`tauPoly`-only)
  Bockstein cycle-selection machinery in `tao_bockstein.cpp` might have its own
  representational gap (analogous to the one just fixed for `phi_beta`, but for CYCLE
  REPRESENTATIVES needing to be a genuine sum of multiple cogenerators, not explored yet);
  (b) something about how `g_prev_pos`/`inv_ind_prev` (canonical preimages, via
  `recover_inv_ind`) is chosen could differ subtly between the two product directions in
  a way not yet identified; (c) the resolution data itself, while presumed ground truth,
  hasn't been independently cross-checked for this specific corner (e.g. verifying
  `qut_9`'s pivot-row corrections at positions 1552/1568 by hand against the raw
  Gaussian-elimination process in `hopf_algebroid/6.h`).
- **Checked and ruled out (partially)**: `tau_table::get_cycles()` (tao_bockstein.cpp:231-
  245) stores `full_cycle` directly for genuine surviving classes (`itm.tag ==
  tau_table_entry::Invalid`, which `{4-6}` is, confirmed by it appearing un-filtered in
  table dumps) with NO reduction/simplification applied — so `cyc[4].at(6) = cog6^t0`
  being a pure singleton is presumably an accurate reflection of whatever `full_cycle`
  was recorded for this entry, not something `get_cycles()` itself truncates. Have NOT
  yet checked further upstream (`make_pretable`, tao_bockstein.cpp:81-140ish) to verify
  `full_cycle` was itself constructed completely/correctly for this entry — that's a
  deeper, not-yet-explored possibility.
- Have NOT yet independently hand-verified `qut_9`'s pivot-row values at positions
  1552/1568 against the Gaussian-elimination construction in `hopf_algebroid/6.h`
  (`make_quotient`/echelon reduction) — this is ground-truth resolution data that's been
  assumed correct (per the earlier scope decision) but not actually re-derived by hand for
  this specific case.
- Reported full findings to the user; awaiting direction on which of these remaining
  threads to pull next.

## User pushback: "it can't be an accident that the product lands in tau-torsion" —
## re-examining every tau-related step at high effort

User strongly suspects a tau-representation/choice bug given the target class is
tau-torsion, and asked to re-examine every step for wrong tau additions/cancellations.

**Re-analyzed `apply_ring_matrix_to_sum`/`cofree_adjoint_row_sum` line by line**: the
"cog17+cog18" value is just an ordinary MULTI-ENTRY sparse vector (two different target
indices, each with a plain single-tau-monomial coefficient) — this is NOT a case that
needed the tauPolySum fix at all (no same-index collision), and would have worked
identically in the old tauPoly-only code. So the level-3→4 cancellation is not a
tauPolySum-introduced artifact — it's a property of raw `qut_9` data, same in old or new
code. Re-confirmed `cofree_adjoint_row_sum` calls `algebroid2vector` with the exact same
(unchanged) function for each vterm, so no tau cross-contamination between the two loop
iterations (each stays within its own target cogenerator's summand).

**New, promising lead**: `recover_inv_ind` picks the FIRST tau^0 singleton preimage found
in `qut_3` for X-index `1416` (position `2442` = `{3-4}+13`) — but a **full sweep**
(`TRACE_QUT_ALL_LVL=3 TRACE_QUT_ALL_TARGET=1416`) shows THREE tau^0 singleton candidates:
`pos=2442 ({3-4}+13)`, `pos=2444 ({3-4}+15)`, and **`pos=2837 ({3-5}+11)`** — the third one
in a totally different cogenerator, `{3-5}`. `recover_inv_ind`'s docstring already flags
this exact ambiguity ("multiple positions can have singleton qut rows for the same
target... only the FIRST is kept"), historically justified by an old (already once-
disproven-elsewhere-this-session) assumption that any two valid singleton preimages
differ by something that cancels out downstream. Given the current bug's shape (matches
exactly this same failure pattern from earlier in the project), this looks like the most
promising lead yet: **check whether using preimage `2837` (inside `{3-5}`'s summand)
instead of `2442` avoids the collision** — i.e. is `phi_beta[3].find(2837)` (under
beta=`{6-9}`) something that does NOT collide under `qut_9`, unlike the `2442` choice?

**Next concrete action**: query `phi_beta[3].find(2837)` (adding a small env-gated
query hook to yoneda2.cpp, `TRACE_PHI_AT_LVL`/`TRACE_PHI_AT_POS`, since no existing tool
dumps an arbitrary phi_beta position on demand) and see whether it leads to a nonzero,
non-cancelling result at level 4.

## ROOT CAUSE CONFIRMED: recover_inv_ind's arbitrary tie-break among MULTIPLE clean
## (tau^0) singleton preimages picks the wrong one

Added `TRACE_QUT_ALL_LVL`/`TRACE_QUT_ALL_TARGET` (full sweep of every position whose
`qut` row is a singleton hitting a given X-index, not just the first found like
`recover_inv_ind`), `TRACE_PHI_AT_LVL`/`TRACE_PHI_AT_POS` (dump `phi_beta[lvl].find(pos)`
on demand), and `TRACE_INJ_AT_LVL`/`TRACE_INJ_AT_X` (dump `inj_lvl.find(x)` on demand) —
all in yoneda2.cpp, small additions to the existing describe/trace pattern.

**Swept `qut_3` for X-index 1416** (the X-index behind cogenerator `{4-6}`, whose
canonical preimage `recover_inv_ind` needs when building level 4 under beta=`{6-9}`):
found **THREE** clean (`tau^0`) singleton preimages, not one:
- `pos=2442` (`{3-4}+13`) — the one `recover_inv_ind` actually picks (lowest position).
- `pos=2444` (`{3-4}+15`) — same owner cogenerator, different offset.
- `pos=2837` (`{3-5}+11`) — a different cogenerator entirely.

**Traced what each choice leads to**:
- `2442` (current code's choice): `im_prev={1552:t^0,1568:t^0}`; `qut_9.find(1552) =
  qut_9.find(1568) = {663:t^1}` — IDENTICAL, so both cancel completely mod 2 →
  `im_ck` empty → final answer **0** (the bug).
- `2837`: owner cogenerator `{3-5}` itself has `phi_beta[3]` completely empty here (lg
  itself is 0 under this beta) → also gives 0, but for an unrelated, uninteresting
  reason (no rescue).
- **`2444` (the fix)**: `im_prev = phi_beta[3].find(2444) = {1554:t^0, 1570:t^0}`.
  `qut_9.find(1554) = {664:t^1, 724:t^0}`, `qut_9.find(1570) = {664:t^1, 725:t^0}` —
  these are NOT identical (unlike the 2442 case): the shared `664:t^1` term cancels, but
  `724`/`725` are each unique and SURVIVE → `im_ck = {724:t^0, 725:t^0}`. Pushing this
  through `inj_10`: `inj_10.find(724) = t^1*{10-22}+2 + t^0*{10-25}`,
  `inj_10.find(725) = t^1*{10-22}+2 + t^0*{10-26}` — the shared `{10-22}+2 @ t^1` term
  cancels between the two, leaving **`main_term_cogens = {25:t^0, 26:t^0}`** — this is
  EXACTLY the same `img` that the independently-verified-correct `{6-9}*{4-6}`
  computation produces, which `find_cycle_sum` correctly reduces to `t^0{10-25}` (via
  `table.at(25) = {25:t^0, 26:t^0}`, a genuine, table-driven cancellation of the `26`
  term). **Manually confirmed: using preimage 2444 instead of 2442 gives exactly the
  correct, expected answer.**

**Conclusion**: `recover_inv_ind`'s documented assumption — "any two singleton qut
preimages of the same target differ by an element of ker(qut) that cancels out
automatically downstream" — is FALSE here, exactly reproducing (in a new form) the same
category of bug already flagged (but only partially fixed, for the clean-vs-twisted case)
earlier this project. Picking the WRONG one of several valid-looking clean preimages
causes a spurious full cancellation instead of the correct partial cancellation that
reveals the right answer. This is very likely THE bug the user suspected — it is entirely
about which qut-preimage is used, i.e. exactly a "wrong choice interacting with a
tau-torsion-adjacent structure" issue, matching the user's instinct.

**Not yet resolved (superseded by a deeper finding below)**: a general, principled fix for
WHY `2444` works and `2442` doesn't. Superseded once the actual structural cause was found.

## ROOT CAUSE, take 2: user's diagram-chase argument nailed it — filtered_reindex silently
## discards genuine chain-map information

User gave a precise math argument: if x,y ∈ G_{s-1} both map to the same target under
qut_{s-1}, then x-y ∈ ker(qut_{s-1}) = image(inj_{s-1}) (exactness), so x-y = inj_{s-1}(z)
for some z; then phi(x-y) = phi(inj_{s-1}(z)), which BY THE CHAIN CONDITION equals a value
built via inj_{s'-1} at the shifted level, hence lands in image(inj_{s'-1}), hence dies
under qut_{s'} (exactness one level up). This REQUIRES two things: (1) exactness
(ker(qut)=image(inj)) at both levels, and (2) that phi actually satisfies the chain
condition EXACTLY on inj's image, not just approximately.

**Re-verified exactness directly** (new `TRACE_EXACTNESS_LVL` sweep, yoneda2.cpp: checks
`qut_lvl(inj_lvl.find(x)) == 0` for every X_lvl basis element x, i.e. image(inj)⊆ker(qut)):
**0 violations at levels 2, 3, 4, 8, 9, 10** (thousands of basis elements each). Also
checked the reverse inclusion via dimension counting (new print in the same block):
`total_rank(G_3) = 4105 = |X_3| (2274) + |X_4| (1831)` EXACTLY — combined with 0
violations, this forces `ker(qut_3) = image(inj_3)` exactly (not just ⊆). **Exactness of
the resolution itself is solid** — this rules out the resolution DATA being at fault
(confirmed freshly, not just assumed, learning from this project's past
misunderstandings here).

**So the bug must be in (2): phi doesn't actually satisfy the chain condition exactly.**
Directly computed `qut_9(phi_beta[3].find(2442) - phi_beta[3].find(2444))` by hand from
already-traced values: `qut_9(1552)=t^1x663`, `qut_9(1568)=t^1x663` (cancel),
`qut_9(1554)=t^1x664+t^0x724`, `qut_9(1570)=t^1x664+t^0x725` (the `x664` terms cancel
between these two, but `x724`/`x725` are unique to each and DON'T cancel) → net result
`{724:t^0, 725:t^0}`, confirmed NONZERO. So `phi_beta[3]` does NOT send this specific
`ker(qut_3)` element to `ker(qut_9)` — a direct, confirmed violation of the chain
condition that SHOULD hold by the user's argument.

**Found the actual mechanism**: scanned every level-3 cogenerator's `TRACE_LIFT_SUM_PRE`
data for cases where `main_term_cogens` (the FILTERED, cogenerator-only reindexing of
`im_gk`) has FEWER terms than raw `im_gk` — i.e., cases where `filtered_reindex(cr, ...)`
silently drops nonzero terms of `im_gk` that don't land exactly on a target cogenerator's
own position. **Found one immediately**: cogenerator `j=8` (X-index 2025, some class
`{3-8}`, cog_pos=3583): `im_gk` = 5 genuinely nonzero terms (`1252^t1 1345^t0 1386^t1
1387^t0 1476^t1`, confirmed via `TRACE_DESCRIBE_POS` that NONE of these 5 positions are
cogenerators — e.g. `1252` = `{9-8}+11`, `1345` = `{9-10}+8`, all "offset"/coaction-
descendant positions) — but `main_term_cogens` (after filtering) has **ZERO** terms,
because none of the 5 land on an actual cogenerator position! Since `inj_row` here is
pure (self-term only, no correction), `combined = main_term_cogens = 0` — i.e.
`phi_beta[3]({3-8}) = 0` even though the chain condition genuinely requires it to be a
5-term nonzero value. **This is a hard, silent information loss**, not a subtle
tau-representation issue: `lift.h`'s whole design (stated explicitly in its own header
comment: "matrix lg ... values are cogenerator indices of F_tgt") assumes `phi` restricted
to a source generator can ALWAYS be expressed purely as a combination of TARGET
generators — i.e. that `im_gk` will always land purely on cogenerator positions. That
assumption is FALSE in this data (confirmed directly), and whenever it's false, the
current code just throws away the non-cogenerator part instead of representing it.

**Mathematical framing**: for a cofree comodule `cof(M)`, a comodule map `cof(M) -> Y` (Y
= `cof(N)`) is freely determined by an arbitrary LINEAR map `M -> (underlying vector space
of Y)` — i.e. generators of the SOURCE can map to ANY element of the TARGET comodule,
including non-generator ("coaction descendant") positions, not just elements of `N` itself.
`lg`'s restriction to "values are cogenerator indices of F_tgt" is consequently NOT
general enough — it's an assumption that happens to hold often (many cogenerators checked
this session had `im_gk` landing purely on cogenerators) but is false in general, and
silently corrupts the result when it fails, exactly the kind of bug that would show up as
"answer depends on an arbitrary choice that shouldn't matter" and "answer incorrectly
zero for a torsion-adjacent class" — matching the user's instinct precisely.

**User correction (important, simplifies the fix)**: my framing of this as a deep
"category-theory gap" was overcomplicated. The actual point: a free/cofree comodule's
elements are all linear combinations of `a[x]` (`a` a dual-Steenrod-algebra element, `x` a
cogenerator) — and a RAW POSITION in the comodule already IS exactly one such `a[x]` pair.
So `phi`'s value at a generator legitimately can be `a[x]` for non-unit `a` (a
non-generator position) — nothing exotic. `filtered_reindex`'s restriction to "value must
be `c[x]` for `c` a bare F2[tau] scalar and `x` itself a generator" is simply too narrow
and wrong, not a sign of a deeper design problem needing a rethink.

**Approved plan** (full text at
`~/.claude/plans/look-at-claude-plans-compiled-hopping-tr-synthetic-puffin.md`): in
lift.h's `_sum` functions only (dead-code plain versions left alone):
1. `lift_one_step_sum`: remove `cr`/both `filtered_reindex` calls; use raw
   `main_term`/`phi_pos` directly; `lg.insert(cog_pos, combined)` stores a raw
   G_tgt-position-space vector.
2. `cofree_adjoint_row_sum`: change `algebroid2vector(contribution,
   F_tgt.position_of_gens[vterm.ind])` to `algebroid2vector(contribution, vterm.ind)`
   since `vterm.ind` is now already the raw position.
3. `lift_first_step_sum`: remove `cogen_rule`/`filtered_reindex(cogen_rule, v)`; use `v`
   (already raw) directly.
4. Update `TRACE_LIFT_SUM_*` prints for the renamed/unfiltered variables.
5. Leave `M_beta`/`build_M` (yoneda2.cpp) unchanged — its own final, output-facing
   cogenerator projection is legitimate and different from the bug.
6. Clear all phi caches again, rebuild, verify `{4-6}*{6-9}` now equals `{6-9}*{4-6}`
   (both `t^0{10-25}`), plus `1 1` regression check.

**Status**: implemented (lift.h: removed `filtered_reindex`/`cr` from `lift_one_step_sum`
and `lift_first_step_sum`; changed `cofree_adjoint_row_sum` to use raw positions). Builds
clean, `1 1` regression check unchanged AND now shows MORE nonzero cogenerators
(TMPCENSUS counts went up, e.g. level 3: 4/18 -> 5/18) confirming previously-lost
information (like `{3-8}`) is now preserved.

**BUT the target bug persists unchanged**: `{4-6}*{6-9}` still `0`; `phi_beta[3].find(2442)`
and `find(2444)` are BYTE-IDENTICAL to before this fix (`{3-4}`'s own lg value was already
"pure"/base-position-only, so this fix never touched this specific chain).

**Follow-on bug found in my own fix, then fixed**: `cofree_adjoint_row_sum`'s
`algebroid2vector(contribution, vterm.ind)` treated `vterm.ind` (now an arbitrary raw
F_tgt position, since lg is no longer cogen-restricted) as if it were a generator's own
base position — but `algebroid2vector`'s second arg MUST be a genuine base (positions are
computed as `mon_index[e]+shift`, only meaningful if `shift` anchors a real generator's
summand). Fixed by decomposing `vterm.ind` via `F_tgt.coaction(vterm.ind)` first (same
coaction machinery already used on the source side) into (algebra element, generator)
pairs, THEN anchoring at each generator's base — this is necessary for correctness
whenever `lg` holds a non-base-position value (like the now-fixed `{3-8}`) that later gets
read back in as input to another extension. User confirmed this should be fixed before
continuing.

**Next**: rebuild, clear caches, retest `{4-6}*{6-9}` vs `{6-9}*{4-6}` and the `1 1`
regression check with BOTH fixes in place.

## Instrumentation currently in the tree (all env-gated, harmless when unset)

- `yoneda2.cpp`: `TRACE_QUT_SWEEP_LEVEL/_COG/_TARGET_X`, `TRACE_MON_NAMES`,
  `TRACE_INJ_QUT`, `TRACE_INV_INJ`, `TRACE_TAG_SEM`, `TRACE_QUT_POS/_LVL`,
  `TRACE_PHI_DUMP`, `TRACE_PHI_TABLE`, `TRACE_PROD_S`/`_A`.
- `lift.h`: `TRACE_LIFT` (per-cogenerator trace budget in `lift_one_step`/
  `lift_one_step_sum`), plus `TRACE_LIFT_PRE`/`_CORR`/`_FINALADD` sub-prints gated on the
  same budget.

## Key facts established (still true)

- `tauPoly` = single tau-monomial (`int16_t` exponent); `tauPolySum` = genuine `F2[tau]`
  polynomial (`polynomial<F2>`), used only for `phi_beta`/`M_beta`'s accumulator — the
  resolution's own `inj`/`qut`/differential data stays `tauPoly`-valued throughout.
- A clean (`tau^0`) preimage is always structurally guaranteed by `make_quotient`
  (`hopf_algebroid/6.h`) and `recover_inv_ind` now hard-aborts if that invariant is ever
  violated — this is load-bearing for `tauPolySum`'s inverse (constant-only) being valid.
- Two representations of chain-map values coexist: "cogenerator-projected" (`M_beta`,
  `fmt()`) vs. "raw position" (`phi_beta` directly) — don't compare across them naively.
- Always verify empirically with concrete printed values before concluding a root cause.
