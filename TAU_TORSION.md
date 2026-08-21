# Tau-torsion: the kernel of multiplication by tau on Ext

This note explains the math behind `tautorsion.cpp`: what it computes, why the
tau-Bockstein machinery it leans on works the way it does, and what's genuinely new versus
already available.

## Usage

`tautorsion` finds tau-torsion in the motivic Ext groups this codebase resolves: it looks
at every filtration/stem/weight, and reports every linearly independent combination of Ext
classes that becomes zero after one multiplication by `tau` — including combinations of
classes that aren't individually torsion on their own. Build and run it against an
already-computed resolution:

```
cmake --build build --target tautorsion
./tautorsion <md> <len>
```

`<md>` and `<len>` are the same resolution-identifying arguments used by `yoneda2` and
`tauBoc` (half the max topological degree, and the resolution length), and the resolution
must already have been built for that `<md>`/`<len>` (via `e2p`, `motTab`, `mr_ex`,
`mr_mot`). Output looks like:

```
tau-torsion  s=8 t=8 w=7  (single-generator):  tau^0{8-1}
tau-torsion  s=3 t=31 w=11  (combination):      tau^4{3-13}+tau^3{3-14}
```

Each line is one basis vector of a kernel at a given `(s,t,w)`: `single-generator` means
one already-tau-torsion class; `combination` means a genuinely new relation among classes
that were not individually known to be torsion. Setting `TRACE_GROUP_SIZES=1` also prints
every bidegree where more than one class is active (useful for seeing where a nontrivial
kernel is even possible), whether or not a kernel vector was found there.

## 0. Setup, briefly

We're computing `Ext_A^{s,t,w}(F2,F2)` for the C-motivic Steenrod algebroid `A`, via a
minimal free resolution `G_0, G_1, G_2, ...` of `F2` (cofree comodules, one per filtration
`s`), with coefficients in `F2[tau]` (`tau` a formal variable of motivic weight). Each
`G_s` is free on a finite set of **cogenerators** in various bidegrees `(t,w)`; how the
resolution itself is built doesn't matter for what follows — just that it exists, and that
each `G_s` comes with a differential `d: G_s -> G_{s+1}` and a known list of cogenerators
with recorded bidegrees.

Fixing `s` and `t`, `Ext_s^{t,*}` — everything in that filtration and stem, as the
motivic weight `w` varies — is a finite-dimensional graded `F2`-vector space, one
dimension-slot per weight. Multiplication by `tau` is a genuine **linear map**

```
    mu_w : Ext_s^{t,w} -> Ext_s^{t,w-1}
```

(`tau` lowers weight by 1, never touches `s` or `t`). "Tau-torsion" means a nonzero class
`x` with `mu_w(x) = 0` — an element of `ker(mu_w)`.

## 1. Why cogenerators alone aren't Ext, and what the tau-Bockstein computation is for

Over a *field*, a minimal free resolution is special precisely because its cogenerators
already *are* a basis of Ext — "minimal" means the differential, written as a matrix, has
no invertible (unit) entries, and over a field every nonzero scalar is a unit, so
"minimal" forces the cogenerator-to-cogenerator part of the differential to vanish
identically. No further homology computation is needed; you just read off the cogenerators.

Over `F2[tau]`, this shortcut breaks. `F2[tau]` has plenty of nonzero non-units — `tau,
tau^2, tau^3, ...` — so "minimal" only rules out *unit* (`tau^0`) entries; a cogenerator's
differential is allowed to be a genuine nonzero multiple of `tau`. So the cogenerators of
`G_s` are only a *candidate generating set* for `Ext_s`, not the answer itself: you still
have to compute the homology of the "cogenerator-restricted differential"

```
    d_cog : cogens(s) -> cogens(s+1)
```

(the actual resolution differential `d: G_s -> G_{s+1}`, restricted to land only on
cogenerator positions of `G_{s+1}` — this restriction is a genuine chain complex in its
own right, since `d^2=0` for the resolution implies `d_cog^2 = 0` too). `Ext_s = ker(d_cog
\text{ at level } s) / im(d_cog \text{ from level } s-1)`. This is the actual, complete
description of what needs computing; the rest of this note is about how.

## 2. The tau-Bockstein idea: Gaussian elimination that tracks tau-valuation

Ordinary Gaussian elimination over a field, applied to `d_cog`, would say: "this
cogenerator's image has a nonzero (hence invertible) leading entry on some target
cogenerator — cancel it, that target is a boundary, done." Over `F2[tau]`, the leading
entry of `d_cog(cogenerator)` isn't automatically invertible; it's `tau^v` times a unit,
for some **valuation** `v >= 0` (`v=0` would mean invertible, ruled out by minimality — so
in fact `v >= 1` always, for any nonzero entry). The elimination has to track this
valuation explicitly instead of just "zero vs. nonzero":

- Process cogenerators one at a time (roughly: level by level, lowest filtration first).
  For a cogenerator `cur`, compute `d_cog(cur)` and look at its current leading term.
- If that leading term's target cogenerator hasn't been claimed by anything yet, register
  the pair: **target = boundary of `cur`, killed at exactly valuation `v`** — record `v`
  as this pairing's `diff_length`.
- If the target *has* already been claimed by an earlier, lower-valuation pairing, subtract
  off the appropriate `tau`-power multiple of that known relation and continue reducing
  (exactly the "eliminate using an already-known pivot" step of ordinary Gaussian
  elimination, just with an extra `tau`-power bookkeeping step).
- If `d_cog(cur)` fully reduces to zero this way, `cur` is a genuine cocycle — nothing to
  pair it with, it's simply a surviving class.

This is "tau-Bockstein" because it's exactly the associated-graded/Bockstein philosophy
applied to the `tau`-adic filtration: instead of asking "is this exactly zero," you ask
"how many powers of `tau` divide it," and read off relations order by order in that
valuation. (The classical motivic tau-Bockstein spectral sequence — relating motivic Ext
to classical Ext via the short exact sequence `0 -> F2[tau] --tau--> F2[tau] -> F2 -> 0` —
is the same idea one level up, applied to Ext groups that are already known; here it's run
at the chain level, directly on the cogenerator complex, as part of *computing* those Ext
groups in the first place.)

Running this elimination on every level classifies each cogenerator `a` into exactly one
of three outcomes:

- **GENUINE**: `d_cog(a)` reduces to zero — a genuine cocycle, never paired with anything.
  Its tower (see below) is infinite: `tau^k*{s-a} != 0` for every `k >= 0`.
- **BOUNDARY-WITH-`diff_length D`**: `a` is the target that got paired with some `cur` one
  level down, at valuation `D`. This means `tau^D * {s-a} = d_cog(cur)/tau^D \cdot tau^D`
  is exactly a boundary — i.e. `tau^D * {s-a} = 0` in `Ext_s`. `{s-a}` is still a genuine,
  *nonzero* Ext class (it's not itself a boundary — it's the class a boundary happens to
  hit), just one with known finite tau-order `D`.
- **NONCYCLE**: `a` itself is a `cur` whose own `d_cog(a)` never reduces to zero (it stays
  nonzero, itself becoming the source of some *other* cogenerator's pairing one level up).
  Such `a` is not in `ker(d_cog)` at all — not a valid Ext class at any weight, full stop.

`yoneda2.cpp` (the Yoneda-product code) already relies on exactly this classification to
decide which cogenerator indices name real Ext classes worth reporting; `tautorsion.cpp`
reuses the identical classification for a different purpose.

## 3. From cogenerators to a basis of `Ext_s^{t,w}`: towers

A cogenerator is not "a basis vector at one weight" — it's the *top* of a tower. A GENUINE
or BOUNDARY-WITH-`diff_length` cogenerator `a` sits at one specific weight `w0(a)` (its
own recorded bidegree), and `tau * {s-a}, tau^2*{s-a}, ...` are further classes at weights
`w0(a)-1, w0(a)-2, ...` — the *same* generator, `tau`-shifted, not new cogenerators. A
GENUINE tower never dies (infinite descending sequence of nonzero classes); a
BOUNDARY-WITH-`diff_length D` tower is exactly `D` classes long (`k=0,...,D-1`), then dies
(`tau^D*{s-a}=0`). So:

```
    basis of Ext_s^{t,w}  =  { tau^k * {s-a} : a a GENUINE-or-BOUNDARY cogenerator with
                                stem t, own weight w0(a), k = w0(a) - w >= 0,
                                and (a GENUINE, or k < diff_length(a)) }
```

## 4. What's actually new: from "trivially multiply by tau" to a real matrix

Multiplying a class by `tau` is definitionally trivial: `tau^k{s-a}` becomes
`tau^{k+1}{s-a}`. The real content is *re-expressing* that shifted class as a combination
of named towers, because the raw shifted vector — viewed as a combination of raw
cogenerator positions, before you know which towers it touches — can genuinely involve
*several different* towers once you allow for the corrections `d_cog`'s own construction
forces (the same kind of correction that produces a `diff_length` pairing in the first
place). Reducing an arbitrary cogenerator-space vector against the known tower
representatives (already-existing code, used elsewhere for a different purpose — decomposing
Yoneda products) can legitimately return something like `tau^3{s-b} + tau^1{s-c}` for two
different generators `b,c`, not just the original tower reasserting itself. (It's equally
normal for it to come back as just the same tower shifted by one — e.g. the unit class
generates a free, infinite tower with nothing else to mix with, so `tau*{unit}` really is
just `tau*{unit}`, not an error.)

One correction has to be applied before trusting this reduced expression: it doesn't know
about `diff_length` cutoffs on its own. A term `tau^e{s-b}` in the raw reduction is
*already zero* if `b` is BOUNDARY-WITH-`diff_length D` and `e >= D` — that's exactly
section 2's pairing fact, applied after the fact to clean up an expression instead of
during the original elimination. Dropping these already-dead terms is the "similar filter"
to the NONCYCLE exclusion above: that one throws away terms that were never valid classes;
this one throws away terms that are valid classes but have already tau-died by the
exponent in question. Only after this drop does a surviving term's presence become a plain
0/1 fact — the exponent on any surviving `b` is forced by bidegree homogeneity (it must
land at the one weight being asked about), so there's no remaining freedom to track.

## 5. The linear algebra

Fix `(s,t,w)`. Collect the domain classes active at weight `w+1` — every `(a,k)` pair with
`w0(a)-k = w+1` — and, for each, its filtered, reduced image at weight `w` (a 0/1 row over
the basis of classes active at `w`, from section 3). This is exactly the matrix of
`mu_{w+1}` restricted to this bidegree, and `ker(mu_{w+1})` is its ordinary `F2` null
space — computed by a small augmented-matrix Gaussian elimination (`[A | I]`, reduce `A`,
track the same row operations on `I`; any row of `A` that goes fully to zero has its
`I`-row as a certified kernel vector). If the null space has dimension `d`, the tool
reports `d` independent kernel vectors forming a basis for it — not a single "torsion
exists here" flag. A kernel vector touching one `(a,k)` reproduces already-known
single-generator torsion (section 2's `diff_length` pairing, rediscovered); one touching
several is the genuinely new "individually-non-obviously-torsion classes summing to a
torsion class" phenomenon this tool exists to find.

## 6. Two scope caveats worth stating plainly

- **Tower depth is capped, not unbounded.** GENUINE towers are infinite in principle; the
  tool only walks each one `D_max` weights past its own top (`D_max` = the largest known
  `diff_length` at that level, or 4 if there is none). Generous and cheap, but not a
  completeness proof — a combination relation requiring a deeper reach would be missed.
- **The last loaded resolution level is skipped.** Classifying a cogenerator as GENUINE
  vs. NONCYCLE needs the *next* level's elimination data, which doesn't exist for the last
  level of a truncated resolution — an expected truncation artifact, not a bug.

## 7. Empirical status so far

Run against `md=40, len=30`, multi-dimensional bidegrees do occur (e.g. `s=3, t=31, w=12`
has two active classes, `tau^4{3-13}` and `tau^3{3-14}`), and in every such case checked so
far, `mu_w` turned out to be full rank — no new combination torsion, only the already-known
single-generator kind. That's a real (if so far negative) answer, not a limitation of the
method: the matrix/kernel computation is genuinely exercised on rank>1 data, it simply
hasn't found a nontrivial combination kernel vector yet in the range checked.
