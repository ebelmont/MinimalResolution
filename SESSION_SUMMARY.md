# Session summary — September 2026

A narrative record of the September 2026 working sessions: what was asked for, what
got built, and the two things most worth knowing afterwards. Companion to
[`PROJECT_LOG.md`](PROJECT_LOG.md), which is the reference version — the full trap
list with citations, the verification status, and the Ravenel confirmations.

## Overview of what happened

**The goal:** the code computed the Adams–Novikov E2 page for the **sphere only**.
It now does it for an arbitrary `BP_*BP`-comodule.

**Built, in order:**

1. General-comodule machinery (`BP_generic_init.*`, `Steenrod_generic_init.*`,
   `BP_mod_I.*`)
2. The `docs/` tree, grounded in `MinimalResolution.pdf` with file:line citations
3. A second Hopf algebroid as an independent test — `F_7` over `Γ = F_7[x]/(x^5)`
4. A comodule **registry with CLI selection** — `./mr_BP_comod <t> <len> [comodule]`,
   replacing two single-purpose drivers
5. The comodules `sphere`, `alpha_1`, `triv_01`

## Two things worth flagging

### The repo has moved well past where these sessions left it

`prime3` is roughly twenty commits ahead with work from later sessions: three more
comodules (`mod_p`, `mod_p_v1`, `alpha_1_mod_p`), a `height` mechanism for modules
free only over `BP_*/I_n`, `ses_check.py` + `docs/SES_CHECK.md`, `example_data/`, a
`.gitignore`, and a much-expanded `anss_chart.py`. `PROJECT_LOG.md` reflects the
current state, not only the part described here.

### A bug introduced here shipped, and was fixed later

Commit `35c6370`. `reduce_BPBP_mod_I` read `BPBP`'s exponent slots backwards,
reducing every off-diagonal coaction entry to **zero** — so phase 2 received a
*split* comodule (the mod-`I` model of `S⁰ ∨ S⁴` instead of `S/α₁`). Two things
worth carrying forward:

- The `sphere`-vs-`mr_BP` regression that this work leaned on **could not possibly
  have caught it**. The sphere's coaction is `1`, which reduces correctly either
  way. A byte-identical regression on a degenerate case proves less than it looks
  like it does.
- It was found by reading the `BPBP`-slot warning written here for an unrelated
  reason. Writing the trap down is what surfaced it.

The damage turned out to be indexing only — generator labels shifted, while
degrees, summand structure and α₁-multiplication edges were unchanged — because the
BP side always used the true coaction. That is luck, not design.

## The traps, in brief

Full versions, with file:line citations, are in `PROJECT_LOG.md`. The ones that cost
the most time:

1. **`BPBP`'s exponent slots are reversed** — outer indexes `v_i` via `η_R`, inner
   indexes `t_i`. So `monomial(singleVar(1,1), unit(1))` is `η_R(v_1)`, *not* `t_1`.
   Use `BP_Op::h0()`. Fails silently: it breaks counitality, and nothing checks the
   comodule axioms.
2. **`pre_resolution_tab` needs a field base ring** — over `BP_*` it compiles and
   returns an *empty* page. Hence the reduce-mod-`I` → resolve → lift architecture,
   which the sphere computation already used.
3. **Degrees are full topological degrees**, and `argv[1]` is `t`, not half of `t`
   (the "halfT" wording was inherited from the p=2 original and is wrong here).
4. **`Fp_Op::inverse` was p=2,3 only** — silently returned `0` otherwise. Now
   general.
5. **`BPtab` must run first**, or you get `std::bad_alloc`. `mr_BP_comod` does not
   need `mr_st`; `mr_BP` does.

## Confirmed against Ravenel

Once `ravenel3rd.pdf` was added to the repo, the S/α₁ mathematics was checked
against it — all six claims hold: Def 1.3.10 for `α_t = δ_0(v_1^t)`, Thm 5.2.6 for
`Ext^{1,4} = Z/3`, Thm 2.3.4 for the geometric boundary theorem, Thm 4.1.19 for the
degrees. The citation table is in `PROJECT_LOG.md`. One nice detail: `BP_Op::h0()`
computes exactly Ravenel's `δ_0(v_1)`, so the code's naming matches the book.

**Still unverified:** nothing checks the comodule axioms, freeness, or the `height`
field — which is silently wrong in *both* directions.
