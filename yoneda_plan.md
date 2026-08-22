## Session handoff files

This file (`CLAUDE.md`) is for durable knowledge that doesn't go stale as the specific bug
being chased changes — math, strategy, and pitfalls every session should know. It should
NOT describe the current specific bug or the exact next debugging step; that's what
`STATUS.md` and any other linked plans are for.

There are two persistent tracking documents, at different levels of granularity, kept
alongside this file:

- **`STATUS.md`** — coarse-grained, updated only when something material changes (a bug
  gets fixed, a new one is found, the active plan changes). Read this FIRST in any new
  session on this problem — it links to the currently-active plan file (the concrete next
  implementation step) and gives a short summary of what's fixed vs. still open. Keep it
  updated whenever status changes; don't let it go stale.
- **`debug_notes.md`** (plus `debug_notes_archive*.md` for history) — fine-grained,
  blow-by-blow trace log. Cycle it (archive current content to the next
  `debug_notes_archiveN.md`, start a short fresh file) when it gets long, the same way
  past cycles are already archived.

## ⚠️ Stop and explain math-level plan changes before proceeding ⚠️

If, while implementing or debugging, you conclude that the **active plan itself** needs a
math-level change (not just an implementation detail) — e.g. a different formula, a
different notion of what needs to be true, a new understanding of what's actually
required mathematically — **stop and explain the issue to the user in terms of the math**
before changing the plan or continuing to implement. Let them approve or edit the revised
plan rather than silently proceeding on your own revised understanding, even if you're
confident in it. You can investigate, add tracing, and write tests without
asking for explicit permission.

## ⚠️ DO NOT TOUCH `yoneda.cpp` ⚠️

**`yoneda.cpp` is explicitly and permanently out of scope. Do not read it, edit it, run it,
or compare against it, even implicitly ("let me just check how yoneda.cpp handles this").**
It is an older, independently-buggy implementation that the user does not want debugged or
used as a reference. `yoneda2.cpp` is the actively maintained rewrite and the only program
that matters for this work. If you find yourself about to `cat`, `grep`, or `Read` a path
containing `yoneda.cpp`, stop — that is almost certainly a mistake. If genuinely unsure
whether something is in scope, ask the user first rather than opening the file "just to
check."

## The math

This is a computation of Ext groups over the C-motivic Steenrod algebra (Adams spectral
sequence machinery), specifically **Yoneda products** of Ext classes, computed via a
**minimal free resolution** of `F2` (or a comodule) over the motivic Steenrod algebroid.

- The resolution is a sequence of "cofree comodules" `G_0 -> G_1 -> G_2 -> ...` with
  differentials `d: G_k -> G_{k+1}`, each `G_k` free on a set of "cogenerators" in various
  internal degrees. Ext classes at resolution-level `s` correspond to cogenerators of `G_s`
  (named things like `{s-i}`, the i-th cogenerator at level s), possibly further reduced by
  a "tau-Bockstein" homology computation (`tables`, `cyc`, `tag_index`/`cycle_index`) since
  not every cogenerator survives as an actual Ext class — some die in a `d_tau`-type
  differential.
- Coefficients live in `F2[tau]`: a polynomial ring over F2 in one variable tau (motivic
  weight). **Ext over `F2[tau]` is a module over the whole polynomial ring**, so a
  coefficient can genuinely be a nontrivial sum like `tau + tau^2`, not just a single
  monomial `tau^n` — keep this in mind before assuming any "single monomial only"
  restriction is a law of the math rather than a simplification of a particular data
  structure (see the `tauPoly`/`tauPolySum` split below, which exists precisely because an
  earlier simplification turned out to be too narrow).
- A Yoneda product `{s2-a} * {s1-b}` (product of two Ext classes) is computed by:
  1. Building a **chain map** `phi_beta: G_k -> G_{k+s1}` lifting multiplication by the
     "beta" class `{s1-b}` through the resolution, degree level by level (`lift_first_step`,
     then `lift_one_step` iteratively for k=1,2,...). See "The lifting strategy" below for
     the construction (represent by the counit projection `phibar_k`, force it on
     `image(inj_k)`, cofree-up to the rest).
  2. Applying `phi_beta` (as a table `M_beta`, cogenerator-index -> position in target
     comodule) to a cycle representative of the other class `{s2-a}`, and reading off which
     target cogenerators appear.

## The lifting strategy, in detail

_(Rewritten 2026-07-11. This is the intended construction; the code in `lift.h` does NOT
yet implement it — see "How this differs from the current implementation" at the end. Plan
file: `~/.claude/plans/look-at-the-top-level-smooth-llama.md`.)_

`lift_one_step` builds `phi_k: G_k -> G_{k+s1}` given `phi_{k-1}`. `phi_k` must satisfy
**two** conditions:

1. **Comodule map**: `phi_k` intertwines the coactions.
2. **Chain square**: `phi_k ∘ d_{k-1} = d_{k-1+s1} ∘ phi_{k-1}` (where `d = inj∘qut`),
   required on **all** of `G_{k-1}`.

The whole trick is to use the universal property of the **TARGET** being cofree (not the
source). Everything below follows from that one choice.

### Step 0 — represent `phi_k` by its counit projection `phibar_k`

Write the target `G_{k+s1} = A ⊗ Z_{k+s1}`, with `Z_{k+s1}` the free `F2[tau]`-module on
the target cogenerators. Set

  `phibar_k := (eps⊗1) ∘ phi_k : G_k -> Z_{k+s1}`

— an **arbitrary** `F2[tau]`-linear map to target cogenerators. Any such `phibar_k`
reconstructs a comodule map by the cofree-up formula

  `phi_k(e_i) = (1 ⊗ phibar_k)(rho_{G_k}(e_i)) = Σ_j a_j · phibar_k(p_j)`

where `rho_{G_k}(e_i) = coaction(i) = Σ_j a_j ⊗ p_j`. **Condition 1 is automatic** for any
`phibar_k`, and `phi_k` is the unique comodule map with counit projection `phibar_k`. So
there is no "extension ambiguity" — `phibar_k` is genuine free data and Condition 1 is
solved by construction.

### Step 1 — only the counit-projected chain condition needs solving (`E = 0` argument)

Consider `E := phi_k ∘ d_{k-1} − d_{k-1+s1} ∘ phi_{k-1} : G_{k-1} -> G_{k+s1}`. It is a
composite of comodule maps, hence a comodule map into the cofree target, hence
`E = (1⊗((eps⊗1)E))∘rho`. Therefore **imposing `(eps⊗1)E = 0` forces `E = 0`** — the full
chain square commutes automatically once its counit projection does. No separate
commutativity check is ever needed.

`(eps⊗1)E = 0` pins `phibar_k` exactly on `image(inj_k)`:

  `phibar_k(inj_k(x)) = prescribed(x) := (eps⊗1)(im_gk)`,

where `im_gk = inj_{k+s1}(qut_{k-1+s1}(phi_{k-1}(preimage of x under qut_{k-1})))`. This is
the quantity the lift already computes. `prescribed` is independent of the preimage choice
because `d∘d = 0` and `phi_{k-1}` is (inductively) a genuine chain map.

- **Clean `tau^0` preimage under `qut_{k-1}`.** `recover_inv_ind` (`lift.h`) finds, for each
  target cogenerator, a position whose `qut` row is a `tau^0` singleton. Such a clean
  preimage is guaranteed to exist by construction of `make_quotient` (`matrices/6.h`); using
  a `tau^i`-twisted one would force a divide-by-`tau^i` that silently zeroes tau-torsion
  values. `recover_inv_ind` hard-aborts if no clean preimage exists.
- **`(eps⊗1)` = keep the offset-0 terms.** `(eps⊗1)` of a target element keeps exactly the
  cogenerator (offset-0) positions and drops everything at positive offset. So
  `prescribed(x)` is just `im_gk` restricted to its offset-0 terms, reindexed to cogenerator
  indices.

### Step 2 — read `phibar_k` off the `inj_k` pivots (a free `F2[tau]`-basis change)

`phibar_k` is forced on `image(inj_k)` and **free everywhere else** (chain-homotopy
freedom). To turn "known on the submodule `image(inj_k)`, free on a complement" into a
per-position table, use the echelon structure of `inj_k`:

- The rows `inj_k(x)` are impure (`inj_k(x) = c_x·e_{pivot_x} + Σ_{q∈N} b_{x,q}·e_q`, where
  `N` = non-pivot positions), so reduce `{inj_k(x)}` to **reduced echelon form** — reusing
  the resolution's own `gaussian`/`symplify_to_led` machinery (`matrices_mem/4.h`,
  `matrices/9.h`). Leading terms / pivots are **gated to `tau^0`** (`symplify_to_led`
  `matrices/9.h:62`,`:95`, via `tauOper::invertible` = "`tau^0` only",
  `mot_steenrod.cpp:44`); a tau-torsion leading term is skipped, never a pivot. So `c_x` is
  a `tau^0` unit and `unify` never divides by a positive tau-power — **no negative-tau
  powers ever enter `inj`/`qut`**, and `{inj_k(x)} ∪ {e_q : q∈N}` is a genuine
  `F2[tau]`-FREE basis of `G_k`.
- Then set, per position:

    `phibar_k(pivot_x) = c_x^{-1} · prescribed(x)`   (`c_x = tau^0`, so no real division)
    `phibar_k(p)       = 0`   for every non-pivot position `p`.

- **Add a hard assert** (à la `recover_inv_ind`) that each `inj_k(x)` pivot entry is a
  single `tau^0` monomial, so any violation of the guarantee aborts instead of silently
  corrupting the basis change.

Note the pivot position is generally a **non-cogenerator offset** (e.g. `tau_1[1-0]` is the
pivot of `inj_1(x5)`), which is exactly why `phibar_k` must be stored at pivots, not at
cogenerator positions.

### Step 3 — cofree-up to all of `G_k`

Apply the cofree-up formula (`cofree_adjoint_row_sum`, `lift.h`) with the full-support,
cogenerator-valued `phibar_k` as `lg`. Because `phibar_k`'s values are **pure cogenerators**
(offset 0), each `a_j · phibar_k(p_j)` is just an algebra element acting on a generator, i.e.
`algebroid2vector(a_j · lambda, position_of_gens[c])` — no `vector2algebroid` decomposition
needed. (`coaction()` on an arbitrary position still computes the full comultiplication and
can blow up combinatorially — keep using `vector2algebroid`/`algebroid2vector` for any
"position ↔ (owner generator, local monomial)" needs, never `coaction()`.)

### Why any valid choice is fine

Different free choices (the `phibar_k = 0` on non-pivots, and the preimage) yield different
but chain-homotopic `phi_k`, hence the same Yoneda product. Exactness (`ker(qut_k) =
image(inj_k)`, re-verified repeatedly — see next section) is what makes `prescribed`
well-defined and the whole scheme consistent.

### How this differs from the current implementation (what to change in `lift.h`)

The current `lift_one_step_sum`/`cofree_adjoint_row_sum` implement the **wrong** universal
property (maps *out of* the cofree **source**: it stores `lg` only on **source
cogenerators** with **raw-position** values and reconstructs `phi(a[y]) = a·lg(y)`, which
does not determine `phi_k`). Concretely, to move to the construction above:

1. **Counit-project `im_gk`** to get `prescribed(x)` (keep offset-0 terms, reindex to
   cogenerators) — instead of storing the full raw `im_gk`.
2. **Store `phibar_k` at `inj_k` pivot positions, 0 elsewhere** — instead of storing only at
   cogenerator positions (`cog_pos`). This requires an echelon reduction of `inj_k`
   (reuse `gaussian`/`symplify_to_led`), or persisting `inj_ind` from the resolution build.
3. **`lg`/`phibar_k` values are pure cogenerators**, so drop the `vector2algebroid`
   raw-position-decomposition path in `cofree_adjoint_row_sum`.
4. **Delete the correction-term machinery** (`subtract Σ c_i·phi(pos_i)`, the
   process-cogenerators-in-order logic) — it was an artifact of solving `phi` at
   cogenerators via impure `inj` rows in the source-out construction, and has no place here.
5. Add the `tau^0`-pivot assert (Step 2).

## Verifying exactness (a good sanity check, not usually where bugs are)

Exactness of the resolution — `ker(qut_k) = image(inj_k)` at every level — is the fact that
makes the "any valid preimage choice gives the same answer" guarantee above hold; it also
gives a genuinely useful debugging tool independent of any specific product computation: if
`x, y` in `G_k` map to the same target under `qut_k`, then `x - y ∈ image(inj_k)`, so for
any genuine chain map `phi`, `phi(x - y)` must land in `image(inj)` one level up, and hence
die under the *next* `qut`. This gives a way to directly test whether a constructed `phi` is
behaving as a true chain map at an arbitrary point, not just by cross-checking two different
orders of a product against each other.

Exactness itself has been checked two independent ways and holds cleanly everywhere tested
so far: (a) a direct sweep confirming `qut_k(inj_k(x)) == 0` for every basis element `x` at
several levels (zero violations), and (b) dimension counting, confirming
`total_rank(G_k) == |X_k| + |X_{k+1}|` exactly. Treat the resolution's own data (`inj`/`qut`/
differentials) as solid ground unless you have comparably strong evidence otherwise — the
bugs found in past sessions have consistently been in how `phi` is built on top of that data,
not in the data itself.

## The algorithm-level machinery (`lift.h`, `matrices/6.h`, `hopf_algebroid/{5,6}.h`)

- Each resolution level `k` stores, alongside `G_k` (`F` in code), two auxiliary matrices
  `inj_k` and `qut_k` built by `quot_index`/`make_quotient` in `matrices/6.h`. These
  represent an injective map (echelon-reduced pivot rows) and its "quotient" complement.
  Crucially — and this was hand-verified by reading `make_quotient` — **every position that
  is NOT a pivot gets `qut(pos) = singleton(i)`, i.e. an EXACT `tau^0` entry, unconditionally,
  by construction.** This is used identically by both `resolvor` (`hopf_algebroid/5.h`) and
  `resolvor_modeled` (`hopf_algebroid/6.h`), so it holds for every `qut` matrix in the whole
  codebase, not just one code path. **This guarantee — a clean tau^0 preimage always exists
  for every quotient target — is the load-bearing fact behind the clean-preimage rule in
  "The lifting strategy" (Step 1) above.**
- `recover_inv_ind(qut, F_rank, C_rank)` (in `lift.h`): for each target index `k` in the
  quotient space, find a position `j` in the previous level whose `qut` row is a `tau^0`
  singleton hitting `k` — i.e. a clean "canonical preimage" of cogenerator/target `k`. This
  preimage is used to compute `im_gk = inj∘qut∘phi_{k-1}(preimage)`, whose counit projection
  is `prescribed(x)` (Step 1). Choosing the clean `tau^0` preimage means **no division by a
  positive tau-power is ever needed** — sidestepping the tau-torsion hazard entirely.
- `tauPoly` (`int16_t`, a bare tau-exponent) represents a single monomial and is used
  throughout the resolution's own data (`inj`/`qut`/differentials). `tauPolySum`
  (`polynomial<F2>`, a genuine multi-term `F2[tau]` polynomial) is used specifically for
  `phi_beta`'s accumulated chain-map values, because — per "The math" above — a correct
  `phi` value sometimes really is a nontrivial sum of tau-powers, and forcing it into a
  single monomial silently corrupts the answer. These two representations legitimately
  coexist; conversions between them (`liftToPolySum` and friends, `mot_steenrod.h`/`.cpp`)
  are a normal, expected part of the code, not a sign that something is half-migrated.

## Pitfalls and lessons learned (read before touching this code again)

1. **Don't assume `tauPoly`'s single-monomial restriction is a mathematical law** — it's a
   simplification that's valid for the resolution's own data but not for `phi_beta`'s
   accumulated values (see above). If you're working with `phi_beta`/`lift_one_step`, you
   should be looking at `tauPolySum`, not `tauPoly`.
2. **Dividing by `tau^n` is dangerous when the dividend might be tau-torsion.** Many Ext
   classes at low degree are killed by tau, so `x / tau^n` is not a safe general operation —
   it can silently produce 0 for a class that should survive. Always prefer a preimage
   choice that avoids the division entirely (see "The lifting strategy" Step 1, and the
   `tau^0`-pivot gate in Step 2) rather than dividing and hoping the tau-torsion issue
   doesn't bite.
3. **A clean (`tau^0`) preimage is always structurally available** for every quotient
   target, *by construction* of `make_quotient` — verified by reading the actual
   matrix-construction code, not assumed. Any place that ends up using a twisted preimage
   when an untwisted one exists is picking the wrong one, and any place where NO clean
   preimage is found for a valid target indicates something is actually broken upstream (a
   genuine invariant violation) — `recover_inv_ind` hard-aborts on this rather than warning.
4. **Keep the two objects straight: `phibar_k` (free data, cogenerator-valued) vs `phi_k`
   (the reconstruction, raw-position-valued).** In the correct construction (see "The
   lifting strategy"), `phibar_k = (eps⊗1)phi_k` is `F2[tau]`-linear into the *target
   cogenerators* — so `lg`/`phibar_k`'s stored values ARE pure cogenerators (offset 0), and
   they are stored on `inj_k` **pivot positions** (generally non-cogenerator offsets), not
   only on cogenerator positions. The **cofree-up output `phi_k`, by contrast, is a general
   raw-position value** (`a[x]` pairings at arbitrary offsets), and must NOT be projected
   down to "cogenerators only" — doing that to `phi_k`'s final value is lossy and silently
   drops real terms. (This is the corrected version of an earlier note that wrongly said
   `lg`'s *own* values must be raw positions — that was the old source-out representation.)
5. **`coaction()` computes the full comultiplication, not a simple "which monomial, which
   generator" lookup** — using it to decompose an arbitrary raw target position causes a
   combinatorial blowup. Use `vector2algebroid` (the actual inverse of `algebroid2vector`)
   for that instead.
6. **Two different value spaces coexist and are NOT directly comparable**: the
   "cogenerator" form (target cogenerator indices — the natural home of `phibar_k` and of
   `prescribed(x) = (eps⊗1)(im_gk)`) versus the "raw position" form (positions in the
   target cofree comodule, including non-cogenerator offsets — the natural home of the
   cofree-up output `phi_k` and of `im_gk` itself). Converting between them is exactly
   `(eps⊗1)` one way (keep offset-0 terms) and the cofree-up the other way; they need the
   right projection/formatting helper before comparison (`fmt()` for cogen-space, a
   `describe_pos`-style helper for raw position + offset). **When debugging phi tables,
   always double check which space you're looking at before drawing conclusions from a
   diff** — e.g. `prescribed` lives in cogen-space, a `phi_k` row lives in raw-position
   space.
7. **Debug output should be human-readable and directly tied to named Ext classes/monomials**,
   not raw indices — dumping `inj`/`qut` matrices restricted to relevant degree ranges with
   a `describe_pos`-style naming scheme (`{level-cogen}` or `{level-owner}+offset`), plus an
   offset→monomial-name table (parsed from `MOP.output_monomials()`'s string output, since
   `mon_array` itself is private), has repeatedly been the most effective way to make sense
   of what's actually happening. Prior instrumentation (grep `yoneda2.cpp` for `TRACE_`) is
   env-var-gated and harmless when unset — check what's already there before adding new
   hooks; there is a lot to reuse (`TRACE_LIFT`, `TRACE_INV_INJ`, `TRACE_INJ_QUT`,
   `TRACE_PHI_TABLE`, `TRACE_MON_NAMES`, `TRACE_QUT_POS`/`_LVL`, `TRACE_QUT_SWEEP_*`,
   `TRACE_QUT_ALL_LVL`/`_TARGET` — a full sweep for every valid preimage of a target, not
   just the first, `TRACE_TAG_SEM`, `TRACE_PROD_S`/`_A`, `TRACE_DESCRIBE_LVL`/`_POS`,
   `TRACE_PHI_AT_LVL`/`_POS`, `TRACE_INJ_AT_LVL`/`_X`, `TRACE_EXACTNESS_LVL` — the
   `qut(inj(x))==0` sweep described above).
8. **Always verify empirically, not by inference/speculation**, especially when reasoning
   about which preimage/representation is "the" canonical one. Multiple bugs in this project
   have only been found because someone insisted on directly checking `qut`/`inj` values and
   tracing through by hand rather than trusting a plausible-sounding story. Back claims with
   printed, hand-checkable values, not just "this should work."
9. **`phi_beta`'s on-disk caches (`<N>_yoneda2_*_phi*`) are binary and tied to its current
   value representation.** Whenever that representation changes, *all* cached phi files
   become stale and must be cleared — not just the specific repro cases you're testing —
   since loading an old-format cache silently gives wrong (not erroring) results.
10. **`mot_steenrod.cpp` is compiled into many CMake targets** besides `yoneda2`
    (`motTab`, `mot_mult`, `dump_gens`, `test_lift`, `tauBoc`, and `yoneda`). Adding a new
    dependency to it (e.g. needing `Fp_Op`) can silently break those other targets' links
    unless their CMakeLists sources are updated too — check for this whenever touching a
    widely-shared file.
11. **Before testing ANY code change, do all three of the following, in order — skipping
    any one of them silently reproduces the OLD behavior and looks exactly like "the fix
    didn't work":**
    1. Rebuild (`cmake --build build --target yoneda2`, or whichever target changed).
    2. Copy the freshly-built binary from `build/` over the repo-root copy actually being
       invoked (e.g. `\cp -f build/yoneda2 ./yoneda2`) — `./yoneda2` at the repo root is a
       separate, stale copy that does NOT get updated by the build step, and running it
       silently tests old code with no error or warning.
    3. Delete the on-disk phi caches (`\rm -f *_yoneda2_*_phi*`, matching `<N>_yoneda2_*_phi*`
       e.g. `40_...`, `20_...`) — see point 9 above; a stale cache silently loads old-format
       or old-logic results instead of recomputing.

    A session already burned real time on this: after implementing a real, correct fix, the
    first several verification runs silently used the stale root-level `./yoneda2` binary
    (built before the fix) and produced output that looked unchanged, appearing to show the
    fix had no effect. Always re-run all three steps above, in order, before drawing any
    conclusion from a test run — and if a change "had no effect," suspect this before
    suspecting the fix itself.
12. **This shell's `rm` and `cp`/`mv` are not the plain coreutils.** `rm` is a shell
    function (likely a trash/safe-delete wrapper) and `cp`/`mv` are aliased to `-i`
    (interactive-confirm) variants. Both print what looks like an interactive
    "Actually remove? (Y/n)" or "overwrite?" prompt when invoked non-interactively (e.g.
    from a tool call with no attached TTY) — in that context the prompt is typically
    harmless/bypassed rather than actually blocking, but it makes the true outcome
    (removed? overwritten? left alone?) ambiguous from the visible output alone.
    **`\rm -f ...` / `\cp -f ...` does NOT reliably bypass this** (re-confirmed in a
    2026-07-11 session: the interactive prompt still appeared despite the leading
    backslash) — use the **literal absolute path** instead, e.g. `/bin/rm -f ...` /
    `/bin/cp -f ...`, which is what actually bypasses the shell function/alias. This
    matters especially for clearing phi caches and copying build output (see point 11).
    After any such command, it's worth double-checking the result (e.g.
    `ls | grep -c ...`) rather than trusting the "Going to remove the following..."
    banner alone.
13. **When solving a linear system via Gaussian elimination for a right-hand-side value
    (not just to find a valid pivot basis), the right-hand side must be transformed by
    the SAME row operations as the matrix, during the SAME pass — not consulted
    afterward as if untouched.** Found in `lift_one_step_sum`'s echelon/pivot
    construction (2026-07-11): the forward pass correctly eliminated earlier pivot
    columns from each row (`rows[j] -= coeff*rows[owner]`) to find a valid pivot
    assignment, but the corresponding `prescribed(j) -= coeff*prescribed(owner)` was
    never applied — silently leaving `prescribed` inconsistent with the reduced rows.
    This is the textbook augmented-matrix `[M | b]` requirement: reducing `M` alone
    finds a valid *basis* (any invertible change of basis is fine for that), but
    reducing `M` while leaving `b` alone breaks the *system* `Mx=b` you're actually
    trying to solve. Symptom was subtle and easy to misattribute: it looked exactly
    like a math-level gap (see `debug_notes.md`'s retracted "phi_0 preimage-independence"
    hypothesis) rather than an arithmetic bookkeeping bug, until traced through by hand
    with concrete positions and coefficients. Whenever a pivot/echelon reduction exists
    only to determine a basis (no rhs to solve for), forward-only elimination is fine;
    the moment there's an rhs vector whose VALUE matters, track and reduce it in lockstep
    with the matrix, every step, not as an afterthought.
