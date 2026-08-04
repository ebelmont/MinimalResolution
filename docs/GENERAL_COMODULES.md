# Computing the E2 Page for a General Comodule

By default, this repo computes the algebraic Novikov E2 page only for the
**sphere**: `mr_BP` resolves the trivial `BP_*BP`-comodule `BP_*` itself (see
`docs/CODE_WALKTHROUGH.md` §3 for exactly where `set_to_trivial` builds
that). This page documents a new, additive capability — `BP_generic_init.h`/
`.cpp`, `Steenrod_generic_init.h`/`.cpp`, `BP_mod_I.h`/`.cpp`,
`comodules.h`/`.cpp`, and `mr_BP_comod.cpp` — for computing the same kind of
E2 page for **any** finitely generated `BP_*BP`-comodule `M` that is free
over `BP_*` (e.g. the BP-homology of a finite complex), given as a rank, a
per-generator degree, and a coaction matrix.

Comodules live in a **registry** (`comodules.cpp`) and are selected by name
on the command line:

```
./mr_BP_comod <halfT> <resolution_length> [comodule]   # default: sphere
./mr_BP_comod --list                                    # show what's available
```

Nothing in the existing pipeline is modified; this is purely additive.

## Why this isn't a one-line change

The obvious approach — build a `comodule_generic<BPBP,int>` from your own
data instead of calling `set_to_trivial`, then call
`Hopf_Algebroid::pre_resolution_tab` directly — **does not work**, and the
reason is worth understanding before using this:

`PolyOp::invertible(x)` (`polynomial/6.h:27-28`), which `curtis_table`'s
pivot-finding (`symplify_to_led`) relies on to decide whether a leading term
is "already accounted for," only recognizes an element as invertible if it's
a bare constant *and* that constant is itself invertible in the base ring.
For `BP = BP_* = Z_(3)[v_1,v_2,...]`, that means an element is invertible
only if it has **no `v_i`'s in it at all**. Almost anything a real coaction
produces (`η_R(v_1)`, anything genuinely involving the `v_i`'s) is never
invertible by this test — so the pivot search that `pre_resolution_tab`
relies on can't reliably find pivots when the ring is `BP_*` itself. This
isn't a corner case: it's why the *existing*, shipped `mr_BP` never calls
`pre_resolution_tab` with `ring=BP` at all — it exclusively uses the
"modeled" path (`pre_resolution_modeled`), which never searches for pivots
over `BP_*`; it reuses pivots already found on the `F_3` side (a genuine
field, where every nonzero element is invertible) and only *lifts* the
resulting values into `BP_*BP` via `BP_Op::lift`.

(This was found by testing, not just reasoning: an earlier version of this
feature that called `pre_resolution_tab` directly over `BP` compiled cleanly
but produced an empty algebraic Novikov page for the sphere, where the known
answer has several classes — a concrete illustration of why every nontrivial
change here should be checked against a case with an already-known answer,
not just compiled.)

## Is "reduce mod `I`, then resolve, then lift" new? No

It's tempting to read the two-phase approach below as a new architecture
introduced for the general case. It isn't — **the sphere's own E2-page
computation already works exactly this way**; see
`docs/CODE_WALKTHROUGH.md` §4.0 for the precise, function-by-function
trace. In brief: `BP_*` is the simplest possible `BP_*BP`-comodule (one
generator, coaction "1 times itself"), so its reduction mod `I` has a
one-line answer that doesn't need computing — that one-liner *is*
`set_to_trivial` (`hopf_algebroid/12.h:1-12`), called on the `P` side by
`mr_st` (`steenrod_init.cpp:32`). `mr_st` resolving that reduction and
`mr_BP` lifting it (`BP_init.cpp:43-54`) are phases 2 and 3 of the exact
same pattern used below. The only genuinely new piece for a general `M` is
`BP_mod_I.cpp`'s reduction functions — a real computation, because a general
`M`'s mod-`I` reduction isn't simple enough to hardcode by hand the way
`set_to_trivial` does for `BP_*` itself.

## The two-phase approach

1. **Reduce `M` mod `I = (p, v_1, v_2, ...)`** to get a comodule over
   `P = BP_*BP/I` — a genuine field-based Hopf algebra (`F_3[t_1,t_2,...]`).
   Since `M` is free over `BP_*`, this reduction is always well-defined, and
   it's exactly the general `M`/`M/I` language `MinimalResolution.pdf` §5
   already uses (not specific to `M = BP_*`). `BP_mod_I.h`'s
   `reduce_coaction_rows_mod_I` derives this automatically from your
   `BP_*BP`-valued coaction — you only enter your comodule's data once.
2. **Resolve that reduction with the classical (Steenrod-style) machinery**
   — `SteenrodGenericInit::set_comodule` + the inherited
   `SteenrodInit::resolve()`, i.e. `pre_resolution_tab` over `ring=Fp`. This
   *does* work correctly, since `F_3` is a field.
3. **Lift that resolution to `BP_*BP` via `pre_resolution_modeled`** —
   `BPGenericInit::set_comodule` + `BPGenericInit::resolve(model_gens_data,
   model_table_file)`, pointed at the files phase 2 wrote, exactly mirroring
   how `mr_BP` already lifts `mr_st`'s resolution for the sphere.

```mermaid
flowchart LR
    A["Your data:<br/>rank, degree,<br/>BP_*BP-valued coaction"] --> B["BPGenericInit::set_comodule"]
    A --> C["reduce_coaction_rows_mod_I<br/>(BP_mod_I.h)"]
    C --> D["SteenrodGenericInit::set_comodule<br/>(comodule over P = BP_*BP/I)"]
    D --> E["SteenrodInit::resolve()<br/>(pre_resolution_tab, correct: F_3 is a field)"]
    E --> F["model gens/table files"]
    B --> G["BPGenericInit::resolve(...)<br/>(pre_resolution_modeled, lifts the model)"]
    F --> G
    G --> H["BPInit::resolution() / make_algNov() / make_Boc()<br/>(unchanged -- only ever read the resolution's own output)"]
```

## How to use it

Build once, then pick a comodule by name:

```
sh BP_comod_compile
./BPtab 20                      # same prerequisite as mr_BP
./mr_BP_comod 20 4              # the sphere (default)
./mr_BP_comod 20 4 alpha_1      # S/alpha_1
./mr_BP_comod --list            # what's available
```

Output is prefixed `<halfT>_<comodule>BP...` for the final resolution and its
`algNov`/`Boc` tables, and `<halfT>_<comodule>P...` for the mod-`I` model
resolution — so different comodules, and a plain `mr_st`/`mr_BP` run, all
coexist in one directory without clobbering each other.

### Shipped comodules

| name | complex | rank | degrees | coaction |
|---|---|---|---|---|
| `sphere` (default) | `S` | 1 | `0` | `ψ(x_0) = 1 ⊗ x_0` |
| `alpha_1` | `S/α₁ = cofib(S³ → S⁰)` | 2 | `0, 4` | `ψ(x_0) = 1 ⊗ x_0`, `ψ(x_4) = 1 ⊗ x_4 + t_1 ⊗ x_0` |
| `triv_01` | `S ∨ S¹` | 2 | `0, 1` | identity: `ψ(x_i) = 1 ⊗ x_i` |

`triv_01` is split (`BP_*(S ∨ S¹) = BP_* ⊕ ΣBP_*`), which makes it a useful
self-test of the rank > 1 machinery: since `Ext` takes direct sums to direct
sums, its E2 page must be *exactly* two copies of the sphere's, one shifted
one step along the stem — see "A sanity check with mathematical content"
below. It also exercises odd-degree bookkeeping, which neither the sphere nor
`alpha_1` touches (everything in `BP_*` sits in even degrees).

### Adding your own

Two steps, both in `comodules.cpp`, and nothing else in the program changes:

1. Write a builder function setting `rank`, `degree`, and `coaction_rows`.
   `coaction_rows(i)` returns generator `i`'s coaction as a sparse
   `vectors<matrix_index,BPBP>` — a list of `(j, c)` pairs, `j` the index of
   another generator (`0..rank-1`), `c` a `BP_*BP` element built with
   `BP_oper.BPBP_opers`'s ring operations (`.unit(1)`, `.add(x,y)`,
   `.multiply(x,y)`, …). Read the warning directly below before writing any
   coefficient by hand.
2. Add one row to `comodule_table` at the bottom of the file.

Degrees are **full topological degrees**, the same units `exponents.cpp`'s
`xnDegs` uses: `|v_n| = |t_n| = 2(3ⁿ−1)`, so `|v_1| = |t_1| = 4` at `p=3`.

`BP_oper`'s structure tables are already loaded when your builder runs, so
`BP_oper.h0()` (`= t_1`) and `BP_oper.thetas()` (`= β₁`, …) are available —
these are the standard connecting-homomorphism constructions
(`α₁ = δ(v_1/p)`, `β₁ = δ(v_2/v_1)`) and are the safest way to get hold of
the cobar cocycles you'll need as off-diagonal entries.

### Warning: `BPBP`'s two exponent slots are not what you'd guess

`BPBP = polynomial<BP>` is **not** "polynomial in the `t_i` with `BP_*`
coefficients" in the layout the type name suggests. Per `BP.cpp:69` — *"the
right unit, `vn` is in the outer"* — it is the other way around:

- the **outer** exponent indexes the `v_i` (included via the **right** unit
  `η_R`), and
- the **inner** (coefficient) `BP`'s exponent indexes the `t_i`.

So `BPBP_opers.monomial(singleVar(1,1), unit(1))` — the obvious-looking way
to write `t_1` — is actually **`η_R(v_1)`**. The correct `t_1` is

```cpp
// t_1: outer exponent 0 (no v's), inner coefficient carrying t-exponent 1
BP   inner = BP_oper.monomial(singleVar(1,1), BP_oper.Z3_oper->unit(1));
BPBP t1    = BP_oper.BPBP_opers.monomial(0, inner);
```

Better still, for `t_1` specifically, just call **`BP_oper.h0()`**, which the
repo already provides: it computes `(η_R(v_1) - η_L(v_1))/p` from the loaded
structure tables, which is exactly `t_1`, and it can't be gotten backwards.
(Verified: `h0()` and the two-line construction above are bitwise equal, and
both differ from `BPBP_opers.monomial(singleVar(1,1),unit(1))`.)

This matters because getting it backwards fails *silently* in a particularly
nasty way: `η_R(v_1)` is not in the augmentation ideal, so a coaction using it
violates counitality — and, per the caveats below, nothing checks that.

## Verifying a change here

Because resolving the trivial comodule through this path should produce
*exactly* the same answer as the existing `mr_BP` (they're computing the
same thing, `Ext_{BP_*BP}(BP_*,BP_*)`, just via different code paths), the
`sphere` comodule — the default — is a real regression test, not a smoke
test. Run both and diff:

```
./mr_st 20 5 && ./BPtab 20 && ./mr_BP 20 4     # reference
./mr_BP_comod 20 4                              # same thing, generic path
diff 20_BPAANSS_table.txt 20_sphereBPAANSS_table.txt
diff 20_BPBocSS_table.txt 20_sphereBPBocSS_table.txt
```

Both are currently byte-identical (as are the `_binary` and `_a0.txt`
variants), at `halfT=20`/`length=4` and at two smaller ranges including one
with a nontrivial differential. This comparison is how the original
from-scratch-over-`BP` approach was caught as incorrect. If you change
anything in `BP_generic_init.*`, `Steenrod_generic_init.*`, or `BP_mod_I.*`,
re-run it before trusting the result on a real complex.

### A sanity check with mathematical content

Diffing against the sphere only tests the plumbing. For a check that the
*comodule data itself* is right, look for something the topology predicts.
For `alpha_1`: coning off `α₁` should kill it, so the class the sphere has
in stem 3, filtration 1 must be absent from `S/α₁`:

```
grep -c "deg=(3,1)" 20_sphereBPAANSS_table.txt    # 1  -- alpha_1 itself
grep -c "deg=(3,1)" 20_alpha_1BPAANSS_table.txt   # 0  -- killed, as it must be
```

For `triv_01` the prediction is sharper still, and needs no topology input:
`Ext` takes direct sums to direct sums, so its page must be exactly the
sphere's page plus a copy of it shifted one step along the stem — for every
sphere class at `deg=(a,b)`, classes at `(a,b)` and `(a+1,b)`, and nothing
else. At `halfT=20, length=4` that holds exactly: 11 sphere classes → 22,
with the degree multiset matching term for term, and the single `d2`
doubling to two.

Find the analogous prediction for your own complex before trusting its
output; it is the only check that can catch a wrong coaction matrix, since
nothing verifies the comodule axioms.

## Caveats

- **No comodule-axiom verification.** Nothing here checks that the coaction
  you supply actually satisfies coassociativity or counitality. A mistake
  won't crash — it will silently produce a wrong `Ext` computation.
  `set_comodule` only checks that your data is structurally well-formed
  (indices in range, rows sorted).
- **Same degree-range and generator-count limits as the rest of the
  pipeline** — bounded by `max_degree`/`monomial_index`, and by the shared
  `exponents.h` table's generator cap (`v_1..v_5`/`t_1..t_5`), exactly as for
  the sphere computation.
- **"Free over `BP_*`" is assumed, not checked.** This matches the
  framework's blanket assumption throughout (`comodule_generic` has no
  notion of torsion), per `MinimalResolution.pdf` Definition 1.
