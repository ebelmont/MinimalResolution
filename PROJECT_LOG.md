# Project log — general comodules, Sept 2026

A record of what was added to this fork during the September 2026 sessions, why,
and the things that were painful to find out. Written both as an overview and as
working context for a later session (including a later Claude session) that has
no access to the original chat transcripts.

It points at the real documentation rather than repeating it:
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) is the map,
[`docs/GENERAL_COMODULES.md`](docs/GENERAL_COMODULES.md) is the reference for
everything below.

## Orientation

The repo computes the **Adams–Novikov E2 page** (`Ext_{BP_*BP}(BP_*, M)`) via the
**algebraic Novikov spectral sequence**, at `p = 3`. Originally it did this for the
**sphere only** — `mr_BP` resolves the trivial comodule `BP_*` itself. The work
below generalised it to an arbitrary comodule `M`, given as a rank, per-generator
degrees, and a coaction matrix.

Two companion references now live in the repo: `MinimalResolution.pdf` (Guozhen
Wang's write-up, which the `docs/` tree is grounded in) and `ravenel3rd.pdf`
(Ravenel's green book, digital 3rd ed.).

## What was added, in order

1. **General-comodule machinery.** `BP_generic_init.*`, `Steenrod_generic_init.*`,
   `BP_mod_I.*`. Lets you resolve a user-supplied comodule instead of the sphere.
2. **Documentation tree.** `docs/` — `ARCHITECTURE`, `FRAMEWORK`, `BUILD_AND_RUN`,
   `CLASSES`, `GLOSSARY`, `CODE_WALKTHROUGH`, `pipelines/*`, plus styled HTML twins
   and `docs/index.html`. Grounded in `MinimalResolution.pdf` with file:line
   citations.
3. **A second Hopf algebroid**, as an independent test of the framework:
   `trunc_hopf.*`, `trunc_init.*`, `mr_trunc.cpp`, `trunc_compile` — the minimal
   cofree resolution of `F_7` over `Γ = F_7[x]/(x^5)`. See
   [`docs/TRUNCATED_HOPF_EXAMPLE.md`](docs/TRUNCATED_HOPF_EXAMPLE.md).
4. **A comodule registry with CLI selection.** `comodules.h`/`.cpp`,
   `mr_BP_comod.cpp`, `BP_comod_compile`. Replaced two single-purpose drivers
   (`mr_BP_generic_example.cpp`, `mr_BP_S_alpha1.cpp`) that each needed their own
   compile script and a source edit to change comodule.
5. **Comodules**: `sphere`, `alpha_1` (`S/α₁`), `triv_01` (`S ∨ S¹`), and later
   `mod_p`, `mod_p_v1`, `alpha_1_mod_p` with a `height` field for modules free only
   over `BP_*/I_n`.
6. **Charting and checking** (`anss_chart.py`, `ses_check.py`, `example_data/`) —
   largely later work; see `docs/CHARTS.md` and `docs/SES_CHECK.md`.

## The traps

These cost the most time. Each is documented where it bites, but collected here
because they are not guessable.

### 1. `BPBP`'s exponent slots are reversed

`BPBP = polynomial<BP>` is **not** "polynomial in the `t_i` with `BP_*`
coefficients". Per `BP.cpp:69` — *"the right unit, vn is in the outer"*:

- the **outer** exponent indexes the `v_i`, via the **right** unit `η_R`
- the **inner** (coefficient) `BP`'s exponent indexes the `t_i`

So `BPBP_opers.monomial(singleVar(1,1), unit(1))` is **`η_R(v_1)`, not `t_1`**.
For `t_1`, use `BP_Op::h0()`, which computes `(η_R(v_1) − η_L(v_1))/p` from the
loaded structure tables and cannot be inverted.

Getting it backwards fails **silently**: `η_R(v_1)` is not in the augmentation
ideal, so the coaction violates counitality, and nothing checks the comodule
axioms. Full warning at the top of `comodules.cpp`.

**This trap produced a real bug that shipped.** `reduce_BPBP_mod_I` was written
for the suggested-but-wrong layout and read both slots backwards, reducing every
off-diagonal coaction entry to **zero** — so phase 2 was handed a *split* comodule
(the mod-`I` model of `S⁰ ∨ S⁴` instead of `S/α₁`). Fixed in commit `35c6370`.
Two lessons:

- The documented `sphere`-vs-`mr_BP` regression **could never have caught it**,
  because the sphere's coaction is `1`, which reduces correctly either way. A
  byte-identical regression on a degenerate case proves less than it appears to.
- It was found by reading the warning above, which had been written for a
  different reason. Writing the trap down was what surfaced the bug.

### 2. `pre_resolution_tab` needs the base ring to be a *field*

`PolyOp::invertible` (`polynomial/6.h:27-28`) only recognises bare constants with
invertible coefficients. Over `BP_* = Z_(3)[v_1,…]` almost nothing a real coaction
produces qualifies, so `curtis_table`'s pivot search cannot work there. This is
why resolving directly over `BP_*` compiles cleanly and silently produces an
**empty** page.

Hence the **two-phase architecture**: reduce mod `I = (p,v_1,…)` to get a comodule
over the field-based `P = BP_*BP/I`, resolve *that* with `pre_resolution_tab`,
then lift via `pre_resolution_modeled`. This is not new for general comodules —
`mr_st`+`mr_BP` already did exactly this for the sphere, with `set_to_trivial`
(`hopf_algebroid/12.h:1-12`) as the degenerate one-line mod-`I` reduction. See
`docs/CODE_WALKTHROUGH.md` §4.0.

### 3. `Fp_Op::inverse` only worked for p = 2 and 3

It special-cased those two primes and silently returned `0` otherwise (printing
`"not implemented!"`). Only ever exercised at `p=3`, so it went unnoticed until
the `Γ = F_7[x]/(x^5)` example. Now general, via Fermat's little theorem
(`x^(p-2) mod p`), which is algebraically identical to the old code at `p=2,3` —
so no behaviour change for existing pipelines.

### 4. Degrees are full topological degrees

`exponents.cpp:31`: `xnDegs = {0,4,16,52,160,484}`, i.e. `|v_n| = |t_n| = 2(3ⁿ−1)`,
so `|v_1| = |t_1| = 4`. The first CLI argument is **`t`, not half of `t`** — the
"halfT"/"half of t" wording was inherited from Guozhen Wang's `p=2` original and
does not hold in this fork. Corrected in `89348f1`.

Also: only generators with `|v_n| ≤ argv[1]` are built at all, so `v_2`/`t_2`
first appear at 16 and `v_3`/`t_3` at 52 — very small runs look empty rather than
wrong.

### 5. `Hopf_Algebroid::maxDeg` is a global truncation, not the algebroid's top degree

Used only as `ranksBelowDeg(maxDeg − deg_x)` (`hopf_algebroid/9.h`, `/11.h`). It
must stay at least as large as the highest degree the resolution reaches; it is
`unsigned`, so too small a value **underflows** rather than erroring cleanly.

### 6. Private-member shadowing in the `*Init` classes

`BPComodInit`/`ComodInit` each declare their own private `coaction_matrix` that
shadows the inherited public `comodule_generic::coaction_matrix` pointer. Reaching
the base member needs explicit qualification:
`comod.BPCoMod_generic::coaction_matrix->construct(...)`.

### 7. Operational gotchas

- **`BPtab <t>` must be run first** for anything on the BP side. Forgetting it
  gives `std::bad_alloc` while "loading delta table", not a clear error.
- `mr_BP` also needs `mr_st` first; **`mr_BP_comod` does not** — it builds its own
  model resolution internally.
- `BPtab` needs **GMP** (`libgmp-dev`); without it you get
  `gmpxx.h: No such file or directory`.
- Compile scripts name the executable with `-o`, which differs from the source
  file name (`mr_BP_comod.cpp` → `./mr_BP_comod`).
- `mult_table()` (multiplication by `h0 = t_1`) is valid for **any** comodule,
  since `Ext(BP_*,M)` is a module over `Ext(BP_*,BP_*)`. `mult_theta()` is **not** —
  it is Moore-spectrum-specific with hardcoded θ degrees (`BP_init.cpp:152`), and
  `mr_BP_comod` deliberately does not call it.

## Verification: what is actually established

The discipline that emerged, after the `pre_resolution_tab`-over-`BP` attempt
compiled cleanly and produced a silently wrong answer: **never trust a change here
because it builds — check it against a case whose answer is already known.**

Established by measurement:

- **`sphere` through the generic path is byte-identical to a real
  `mr_st`+`BPtab`+`mr_BP` run**, across all eight output tables. This is the
  plumbing regression. Note trap 1: it is structurally blind to coaction bugs.
- **`alpha_1` kills α₁**: the sphere has exactly one class at `deg=(3,1)`, `S/α₁`
  has none. Coning off α₁ must kill it.
- **`triv_01` splits exactly**: `Ext(BP_* ⊕ ΣBP_*)` must be two copies of the
  sphere's page, one shifted a step along the stem. Verified as a degree multiset —
  11 classes → 22, matching term for term, single `d2` doubling to two.
- **`trunc_hopf` reproduces `MinimalResolution.pdf`'s hand-worked example**: one new
  generator in degree 1, none in degrees 2–4, for the first resolution step.
- **The S/α₁ mathematics was checked against Ravenel** (see below).
- Later work added `ses_check.py` and `docs/SES_CHECK.md`, which checks `S/α₁`'s E2
  page against the sphere's through the cofiber long exact sequence — the
  independent check that was flagged as missing at the time.

Ravenel citations, all confirmed against `ravenel3rd.pdf`:

| claim | where |
|---|---|
| `(A,Γ) = (BP_*, BP_*BP)`, `dim v_i = dim t_i = 2(pⁱ−1)` | Thm 4.1.19 |
| `α_t = δ_0(v_1^t) ∈ E_2^{1,qt}`, `q = 2p−2` | Def 1.3.10 |
| `Ext^{1,t}(BP_*) = Ext^{0,t}(M¹)` for `t>0`, `p` odd; `= Z/p^{i+1}` at `t = sp^i q` | Thm 5.2.6(i),(ii) |
| `ᾱ_1 = h_{10}`, and `h_{1,0}` is represented by `[ξ_1]` (BP analogue `[t_1]`) | lines 9781, 4647, 4979 |
| geometric boundary theorem, hypothesis `E_*(h) = 0` | Thm 2.3.4 |
| Hazewinkel `η_R(v_1) = v_1 + p·t_1` | A2.2.1 (explicit as `v_1 + 2t_1` at `p=2`) |

At `p=3`: `t = q = 4` gives `s=1, i=0`, so `Ext^{1,4}(BP_*) = Z/3`. Note that
`BP_Op::h0()` computes exactly Ravenel's `δ_0(v_1)` — the code's naming matches
the book.

**Still not independently verified**: nothing checks the comodule axioms for a
hand-supplied coaction. Freeness over `BP_*` (or `BP_*/I_n`) is assumed, not
checked. The `height` field is likewise unchecked and silently wrong in both
directions — see the comment on `ComoduleSpec::height` in `comodules.h`.

## How the S/α₁ comodule was derived

The method generalises, so it is worth stating rather than just the answer.

1. **Underlying module from the cofiber LES.** `BP_*` is evenly graded, so
   `BP_3 = 0`, so `(α₁)_* = 0` on BP-homology and the LES collapses to
   `0 → BP_* → BP_*(S/α₁) → Σ⁴BP_* → 0`. Rank 2, degrees 0 and 4.
2. **Extension class = the attaching map's E2 class** (geometric boundary theorem,
   Ravenel Thm 2.3.4, whose hypothesis is exactly step 1's conclusion). Here that
   is `[t_1] ∈ Ext^{1,4} = Z/3`.
3. **Write the matrix**: diagonal `1`, off-diagonal the cobar cocycles.

$$\text{sphere: } \begin{pmatrix} 1 \end{pmatrix} \qquad S/\alpha_1: \begin{pmatrix} 1 & 0 \\ t_1 & 1\end{pmatrix}$$

The answer is forced up to a unit: in degree 4 the reduced part of `BP_*BP` over
`BP_*` is spanned by `t_1` alone, so the only content is that the extension is
**nonsplit** — i.e. `α₁ ≠ 0`.

Do not confuse `S/α₁` (cofiber of a map of **spheres**) with the cofiber of the
Adams self-map `Σ⁴V(0) → V(0)`, which is `V(1) = S/(3,v_1)`.

## Quick reference

```sh
sh BPtable_compile && sh BP_comod_compile    # needs libgmp-dev
./BPtab 30                                   # once per degree bound
./mr_BP_comod --list
./mr_BP_comod 30 10 alpha_1                  # <t> <resolution_length> [comodule]
python3 anss_chart.py 30 -c alpha_1          # -> 30_alpha_1anss_E2.svg
```

Adding a comodule is two steps in `comodules.cpp`: write a builder, add a row to
`comodule_table` (with its `height`). Nothing else changes.

The original sphere output, before any of this, is `<t>_BPAANSS_table.txt` from a
plain `mr_st`+`BPtab`+`mr_BP` run. Output files are not committed; `sr` holds the
author's own reference invocation.
