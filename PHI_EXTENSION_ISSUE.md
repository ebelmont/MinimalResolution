# ✅ RESOLVED — kept for the record; see resolution below

**Status (2026-07-11): SOLVED.** The obstacle described in this file is real, but it was
being attacked from the wrong side. The resolution, worked out and stress-tested with the
user this session, is in `~/.claude/plans/look-at-the-top-level-smooth-llama.md` (and
summarized in `STATUS.md`). Short version:

- The extension is underdetermined **because the code used the wrong universal property** —
  maps *out of* the cofree **source** (`φ(a[y])=a·lg(y)`), which genuinely doesn't
  determine `φ_k`.
- Use the **target's** cofree property instead: represent `φ_k` by
  `φ̄_k := (ε⊗1)∘φ_k : G_k → Z_{k+bs}` (arbitrary F2[τ]-linear), reconstruct
  `φ_k=(1⊗φ̄_k)∘ρ`. Then the comodule-map condition is automatic, and the FULL chain
  square follows from just its counit projection (`E=φ_k∘d−d∘φ_{k-1}` is a comodule map
  into a cofree target with `(ε⊗1)E=0`, so `E=0`).
- `φ̄_k` is forced on `image(inj_k)` (= `(ε⊗1)(im_gk)`) and free elsewhere: reduce
  `{inj_k(x)}` to τ⁰-pivot echelon form, set `φ̄_k(pivot_x)=prescribed(x)`, 0 on
  non-pivots. This is a direct construction (no linear solve, no τ-division), and it
  reproduces `φ_1(τ_1[1-0])={7-10}` — the exact case this document said was dropped.

Everything below is the ORIGINAL open-problem write-up. It correctly identifies the
underdetermination and (in the "actual obstacle" section) why *local reasoning on `G_k`'s
own coaction* can't fix it — which is precisely why the target-cofree reframing was
needed. It is retained as context; the "Open question / Suggested next step" sections at
the bottom are now answered by the plan above.

---

# (ORIGINAL) Open issue: `φ`'s extension from cogenerators to non-cogenerator positions is
# underdetermined by local (degree-by-degree) reasoning alone

## Context / where this fits

This is a write-up of a specific mathematical obstacle found while debugging why
`{4-6}*{6-9}` and `{6-9}*{4-6}` (Yoneda products in `yoneda2.cpp`) disagree. See
`STATUS.md` / `debug_notes.md` / `CLAUDE.md` in this repo for the surrounding project
context (C-motivic Adams spectral sequence, minimal free resolution over the motivic
Steenrod algebroid, Yoneda products via chain maps `φ_beta`). This file is meant to be
readable on its own — it re-derives everything needed from scratch with concrete,
independently-verified numbers.

**This is a math-level finding, not yet a fix.** Per this project's standing rule (see
`CLAUDE.md`), any actual code change implementing a new algorithm should be proposed and
discussed before being implemented — this write-up is the "propose/discuss" artifact.

## Background: the two conditions φ must satisfy

`φ_k: G_k → G_{k+bs}` (built for a fixed "beta" class `{bs-bb}`) must satisfy **two**
independent conditions:

1. **Comodule map**: `φ_k` intertwines the coactions, i.e. `Δ_target ∘ φ_k = (1⊗φ_k) ∘
   Δ_{G_k}`.
2. **Chain condition** (the actual thing we care about): `φ` commutes with the
   resolution's differential, factored via `qut`/`inj`. Concretely, for consecutive
   levels, this is the square
   ```
   G_{k-1} --qut_{k-1}--> X_k --inj_k--> G_k
      |                                   |
    φ_{k-1}                              φ_k
      |                                   |
      v                                   v
   G_{(k-1)+bs} --qut_{(k-1)+bs}--> X_{k+bs} --inj_{k+bs}--> G_{k+bs}
   ```
   i.e. `φ_k ∘ (inj_k ∘ qut_{k-1}) = (inj_{k+bs} ∘ qut_{(k-1)+bs}) ∘ φ_{k-1}`, required to
   hold on **every element of `G_{k-1}`**, not just cogenerators.

**⚠️ Corrected claim (an earlier version of this document got this wrong — flagged and
fixed by the user)**: condition (1) alone does **NOT** uniquely determine `φ_k` from its
values on cogenerators, even though `G_k` is cofree. It's true that an arbitrary choice
of generator values freely extends to *at least one* valid comodule map (via the standard
`φ(a[x]) = (1⊗φ)(Δ_{G_k}(a[x]))` formula, evaluated using only the generator values). But
that extension is **not unique** — there exist other, genuinely different comodule maps
agreeing on every cogenerator. Concretely (worked out below): whenever the induction
solving for `φ_k(a[x])` from lower-degree data reduces to an equation of the form
`Δ_target(w) = 1⊗w` (no other terms), **every** `F2[τ]`-combination of target cogenerators
satisfies this equation, not just `w=0` — because `Δ_A(1)=1⊗1` trivially, so
`Δ_target(cogenerator) = 1⊗cogenerator` exactly for *any* cogenerator. So the comodule-map
condition, applied locally, leaves this choice completely open; it takes the *chain*
condition (2) — via information from a different, unrelated part of the diagram — to pin
it down. (The original mistake: treating injectivity of `Δ_target`, via the counit, as
implying the equation `Δ_target(w) = 1⊗w` has a unique solution. Injectivity only rules
out two *different* values of `Δ_target` colliding; it says nothing about how many `w`
satisfy one *specific* target value, and here many do.)

**The current code (`lift.h`, `lift_one_step_sum`/`cofree_adjoint_row_sum`) only ever
solves condition (2) at the cogenerators of `G_k`** (via `recover_inv_ind`'s single
canonical preimage per generator), then extends to the rest of `G_k` using ONLY condition
(1), via the naive formula `a[x] ↦ a·φ_k(x)` (valid when `a` is primitive, see below). Per
the correction above, this produces *a* valid comodule map, but — since condition (1)
doesn't uniquely determine the extension — there's no reason to expect it's *the* one that
also satisfies condition (2) everywhere, and indeed it isn't. The rest of this document is
a fully-worked, independently verified counterexample.

## Notation

- `{k-j}` = the `j`-th cogenerator of `G_k`.
- `{k-j}+n` = the position at local offset `n` within cogenerator `j`'s summand of `G_k`
  (offset `n` indexes a basis monomial of the dual motivic Steenrod algebra, via
  `mon_array`/`mon_index` in `mot_steenrod.cpp`).
- **Offset-to-generator dictionary** (motivic dual Steenrod algebra basis, low degree):
  `offset 0 = 1` (unit), `offset 1 = τ_0`, `offset 2 = ξ_1`, `offset 3 = τ_0ξ_1`,
  `offset 4 = τ_1`, `offset 6 = τ_0τ_1`. **Do not trust the generic `x_1^k`/`x_2^k`
  labels printed by the existing `TRACE_MON_NAMES` instrumentation in `yoneda2.cpp`** —
  those are placeholder names assigned in whatever order `mon_array` happens to list
  algebra elements, and do **not** correspond to the actual `τ_i`/`ξ_i` generator
  identities. (E.g. `TRACE_MON_NAMES` calls offset 2 `x_1^2` and offset 4 `x_2^1` — this
  session initially took that literally and got the coproduct structure wrong as a
  result. The correct names come from motivic Milnor basis theory, not from that tool.)
- `τ` (no subscript) is the *coefficient ring* variable (`F2[τ]`, the polynomial ring
  Ext is a module over) — unrelated to the algebra generators `τ_0, τ_1` above, which are
  elements of the dual Steenrod algebra itself. This is a genuine, confusing notational
  collision baked into the field; keep the two `τ`s mentally separate.

## The concrete example

Working beta = `{6-9}` (`bs=6, bb=9`), i.e. `./yoneda2 40 30 4 6 6 9`'s beta side. All
values below are computed by (already-existing or session-added) `TRACE_*` env-var
instrumentation in `yoneda2.cpp`/`lift.h` and cross-checked by hand; see "How to
reproduce" below.

### What's already correctly solved (cogenerators of `G_1`)

Via the existing `lift_one_step_sum`/`recover_inv_ind` machinery (condition (2), solved
correctly for cogenerators):
```
φ_1({1-0}) = 0
φ_1({1-1}) = t^0·{7-9}          (cog_pos of {7-9} = position 1944 in G_7)
```

### The problem: extending to non-cogenerator positions of `{1-0}`'s summand

We need `φ_1(ξ_1[1-0])` (offset 2) and `φ_1(τ_1[1-0])` (offset 4). **The correct
answers are**:
```
φ_1(ξ_1[1-0]) = 0
φ_1(τ_1[1-0]) = t^0·{7-10}
```
The current code gets the first one right (by accident, via the naive rule below) but
gets the second one **wrong** (naively computes `0`, when it must be `{7-10}`).

### Why the naive "comodule map only" formula fails

`cofree_adjoint_row_sum`'s current formula for `φ(a[x])`, when `a` is a *primitive*
element (`Δ(a) = a⊗1 + 1⊗a`, no middle terms), reduces to the naive rule `φ(a[x]) =
a·φ(x)`. `τ_0` and `ξ_1` are primitive, so this rule is exactly right for offsets 1 and 2
in isolation. **But `τ_1` is NOT primitive**: the motivic Milnor coproduct gives
```
Δ(τ_1) = τ_1⊗1  +  ξ_1⊗τ_0  +  1⊗τ_1
```
— a genuine middle term `ξ_1⊗τ_0`. The current code's `cofree_adjoint_row_sum` only ever
populates its lookup table (`lg`) at cogenerator (offset-0) positions, so when it
processes this middle term, it looks up `lg` at a *non-cogenerator* position (offset 1,
i.e. `{1-0}+1`, matching the `τ_0` piece) and finds nothing — silently dropping real
information. That is the immediate mechanism of the bug, but see the next section for why
naively "fixing" this by populating `lg` at offset 1 too **isn't enough** either.

### The diagram chase that gives the actual right answer

Independently of any coaction reasoning on `G_1` itself, look at position `{0-0}+6`
(`= τ_0τ_1`, in `G_0`) and push it through the actual chain-condition square:
```
qut_0({0-0}+6) = t^0 · x5                     (X_1 basis index 5)
inj_1(x5)      = t^0 · {1-0}+4  +  t^1 · {1-1}+2      (= τ_1[1-0] + t^1·ξ_1[1-1])
```
So the LHS of the square, `φ_1(inj_1(qut_0({0-0}+6)))`, expands to
```
φ_1(τ_1[1-0])  +  t^1 · φ_1(ξ_1[1-1])
```
The second summand is already fully known (via the naive/correct rule, since `ξ_1` is
primitive and `φ_1({1-1}) = {7-9}`): `φ_1(ξ_1[1-1]) = {7-9}+2 = 1946`. Meanwhile the RHS
of the square (computed straight from `φ_0`, independently, via `inj_7(qut_6(φ_0({0-0}+6)))`)
comes out to `t^1·1946 + t^0·2028` where `2028 = {7-10}`. Requiring LHS = RHS:
```
φ_1(τ_1[1-0]) + t^1·1946  =  t^1·1946 + t^0·2028
⟹ φ_1(τ_1[1-0]) = t^0·2028 = t^0·{7-10}
```
This is a **fully independent, externally-forced** equation for `φ_1(τ_1[1-0])`, coming
from `G_0`'s structure — not from anything internal to `G_1`'s own coaction — and it
gives the nonzero answer directly. (`φ_1(ξ_1[1-0]) = 0` also needs justifying, but that
one *is* consistent with the naive rule; the point of concern is that you can't tell
*from local reasoning on `{1-0}`'s summand alone* that offset 2 is fine and offset 4 is
not — see next section.)

## The actual obstacle: local reasoning is underdetermined

The natural next idea is: don't just apply the naive rule blindly — instead, for each
non-cogenerator position `a[x]`, write out the *full* coaction `Δ_{G_1}(a[x])`, subtract
off the known lower-degree pieces (already solved, by induction on degree of `a`), and
solve the resulting equation `Δ_target(φ(a[x])) = (\text{known RHS})` for the unknown
`φ(a[x])`.

**This does not work as a self-contained induction, and here is exactly why.** Take
`w = φ_1(ξ_1[1-0])` (offset 2; `ξ_1` is primitive, `Δ(ξ_1)=ξ_1⊗1+1⊗ξ_1`, no middle term).
Writing out `Δ_{G_1}(ξ_1[1-0]) = ξ_1⊗(1[1-0]) + 1⊗(ξ_1[1-0])` and requiring
`Δ_target(φ_1(ξ_1[1-0])) = (1⊗φ_1)(\text{that})`, the first term contributes
`ξ_1⊗φ_1({1-0}) = ξ_1⊗0 = 0` (using the already-known `φ_1({1-0})=0`), leaving:
```
Δ_target(w) = 1⊗w        (exactly — no other terms)
```
This equation is satisfied by `w=0` — but **also** by `w` equal to *any* `F2[τ]`-linear
combination of `G_7`'s own cogenerators: for a pure cogenerator `c` (offset 0 in its own
summand), `Δ_target(c) = Δ_A(1)⊗c = (1⊗1)⊗c = 1⊗c` exactly, trivially, with zero cross
terms — so `Δ_target(w)=1⊗w` holds for `w=c`, and (by linearity) for any combination of
cogenerators, not just `w=0`. So the equation does **not** pin down `w`; it only tells you
`w` must be *some* cogenerator combination (possibly the zero one). Nothing internal to
this equation, or to processing offsets in increasing degree order within `{1-0}`'s own
summand, resolves which combination is correct.

The only way found so far to actually resolve it is the diagram chase above — pushing a
**higher-degree** source position (`{0-0}+6`, degree of `τ_0τ_1`) through the chain
condition, which ties together `φ_1(τ_1[1-0])` and `φ_1(ξ_1[1-1])` in one equation with a
fully-known RHS. Note the source position used (`{0-0}+6`) has *higher* degree than the
unknown being solved for (`{1-0}+4`) — so this is not a simple "process by increasing
degree" fix either; the informative equation comes from a different, higher-degree
element of the *source* level entirely.

## Practical scoping observation (not a theorem — a simplification worth exploiting)

In this example, the trouble is entirely confined to `{1-0}`'s summand (`φ_1({1-0})=0`):
`φ_1(a[1-0])` needs real solving for at least one `a` (`a=τ_1`). `{1-1}`'s summand
(`φ_1({1-1})={7-9}≠0`) is *not* a problem here — `φ_1(ξ_1[1-1])` computed via the naive
rule (`ξ_1·φ_1({1-1})`) already matches the value required by the diagram chase, with no
correction needed. **This is not a proven general fact** — the same "homogeneous
ambiguity" derived above (any target-cogenerator combination solves the local equation)
applies just as much when `φ_k(cogen)≠0` as when it's `0`; nothing in the coaction
equation itself rules out a nonzero correction in the `{1-1}` case either, it just so
happens empirically that none was needed there. Nonetheless, **if a solution ends up using
something like the naive rule as a fallback/starting point, it should specifically target
its solving effort at cogenerators with `φ_k(cogen)=0`** (the naive rule gives the
uninteresting, easily-*wrong* answer `φ_k(a[cogen])=0` for every `a` there, whereas for
`φ_k(cogen)≠0` the naive rule has, at least in every case checked so far, matched the
required answer) — narrowing the "must actually solve" set to zero-valued cogenerators
would be a substantial practical simplification, worth building into any implementation
even without a proof that it's always sufficient. Worth checking empirically (e.g. via
`TRACE_SQUARE_CHECK`) whether a nonzero-cogenerator counterexample ever actually occurs
before relying on this too heavily.

**Open question**: is there a systematic (not ad hoc / case-by-case diagram-chase)
algorithm that correctly determines `φ_k` on all of `G_k`, given `φ_{k-1}` and the
resolution's own `qut`/`inj` data? Ideas considered and their status:

- *Solve `φ_k` on `image(inj_k)` first (as a whole subspace, not just cogenerators), by
  sweeping every position `p` of `G_{k-1}` and treating `φ_k(inj_k(qut_{k-1}(p))) =
  inj_{k+bs}(qut_{k-1+bs}(φ_{k-1}(p)))` as one linear equation per `p`, then solve the
  resulting system simultaneously.* — Not yet tried; this is the current leading idea
  (raised in this session, un-implemented). It sidesteps the "which local candidate is
  right" ambiguity by using the *actual* forcing equations (from `G_{k-1}`) instead of
  guessing from `G_k`'s internal coalgebra structure. Open question: is this system
  always exactly-determined (enough independent equations, no leftover freedom), or can
  it *also* leave some positions of `G_k` undetermined (in `ker(qut_k)`'s complement or
  similar), requiring a further step?
- *Reuse `curtis_table`'s row-reduction (`symplify_to_led`) machinery, as used to build
  `inj`/`qut` in the first place.* — Considered and set aside: `curtis_table` (via
  `resolvor_modeled`/`embed2cofree_modeled`) is specifically used to find **new
  primitives** during incremental resolution construction, working in an intermediate
  "F2 first, lift `τ` back in afterward" representation. That's a different task
  (detecting new generators) from ours (extracting a `τ`-correct value for a *known*
  target subspace), and reusing its F2-first approach risks silently mishandling `τ`
  exactly the way the next bullet warns about.
- **`τ`-divide danger, explicitly flagged by the user**: any step that concludes
  "`φ_k(τ·y)` is known to be `X`, therefore `φ_k(y) = X/τ`" is unsound in general —
  `τ`-multiplication is not injective (there is real `τ`-torsion), so knowing the image
  of `τ·y` under `φ_k` does not determine the image of `y`. This exact mistake caused a
  previous bug in this project; it was only safe there because the object being divided
  was *provably* a cogenerator with a guaranteed clean (`τ^0`) preimage (see
  `CLAUDE.md`, `recover_inv_ind`). No equivalent safety guarantee is currently known for
  the general non-cogenerator-position case. Any proposed algorithm needs an explicit
  answer to: when (if ever) does solving require dividing by `τ`, and how do we know it's
  legitimate?

## How to reproduce every number above

All via existing (or this-session-added) env-gated trace hooks in `yoneda2.cpp`/`lift.h`.
**Important**: `./yoneda2` at repo root must be freshly rebuilt (`cmake --build build
--target yoneda2`, then `\cp -f build/yoneda2 ./yoneda2` — the aliased `cp`/`rm` in this
shell print interactive-looking prompts even non-interactively, use `\cp -f`/`\rm -f`) and
all `*_yoneda2_*_phi*` cache files cleared (`\rm -f *_yoneda2_*_phi*`) before trusting any
run — see `CLAUDE.md` pitfall #11/#12.

```
# cogenerator values (via existing TRACE_PHI_AT_LVL/_POS, after building phi_beta for 6_9):
./yoneda2 40 30 4 6 6 9   # builds and caches phi_beta for beta={6-9}
TRACE_PHI_AT_LVL=1 TRACE_PHI_AT_POS=0   ./yoneda2 40 30 4 6 6 9   # -> phi_1({1-0}) = 0
TRACE_PHI_AT_LVL=1 TRACE_PHI_AT_POS=974 ./yoneda2 40 30 4 6 6 9   # -> phi_1({1-1}) = {1944:t^0} = {7-9}

# qut_0 / inj_1 raw data:
TRACE_QUT_LVL=0 TRACE_QUT_POS=6      ./yoneda2 40 30 1 1   # -> qut_0(6) = x5
TRACE_INJ_AT_LVL=1 TRACE_INJ_AT_X=5  ./yoneda2 40 30 1 1   # -> inj_1(5) = t^0*{1-0}+4 + t^1*{1-1}+2

# offset naming (CAUTION: labels are placeholder names, not real generator names -- see above):
TRACE_MON_NAMES=10 ./yoneda2 40 30 1 1

# the full square check that originally surfaced this (added this session, yoneda2.cpp):
TRACE_SQUARE_CHECK=30 ./yoneda2 40 30 4 6 6 9
  # reports FIRST failure at level k=0, position p=6 ({0-0}+6), with full LHS/RHS/DIFF dump
```

`TRACE_SQUARE_CHECK` (new this session) directly implements the square-commutativity test
described above (`φ_{k+1}(inj_{k+1}(qut_k(p)))` vs `inj_{k+bs+1}(qut_{k+bs}(φ_k(p)))`),
scanning `k=0,1,2,...` and `p=0..G_k.total_rank-1` in order, and is the tool that found
this specific counterexample as the very first failure in the whole computation (level 0,
before any of the higher-level preimage-ambiguity symptoms found earlier this session,
which are likely downstream consequences of this same root issue rather than independent
bugs — not yet confirmed).

## Suggested next step for a new session

Start from the "solve `φ_k` on `image(inj_k)` as one linear system swept over all of
`G_{k-1}`" idea above; work out on paper (or by extending this same `{0-0}` example a
little further, e.g. also including offset 3 = `τ_0ξ_1` and offset 5) whether that system
is always exactly-determined, or whether it can also leave residual freedom — and if so,
what *additional* principle (beyond "solve the linear system") pins down the remainder.
