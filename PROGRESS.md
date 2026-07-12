# Yoneda Product Implementation Progress

## Goal

Compute the Yoneda product `{s-i} * {s'-i'}` in the motivic Adams spectral sequence,
expressed as a linear combination of generators `tau^k * {a-b}` in Ext^{s+s'}.

### Mathematical description

- `{s'-i'}` = i'-th generator in filtration s' (cogenerator i' of G_{s'})
- `{s-i}` = i-th generator in filtration s (cogenerator i of G_s)
- Product lives in Ext^{s+s'}
- Computation: pick out generator i of G_s via the map F_2[tau] → G_s (unit → position_of_gens[i]),
  compose with the lifted chain map phi_s: G_s → G_{s+s'}, project result to cogenerators of G_{s+s'}.

### Lift algorithm (existing, in lift.h)

For a fixed element {s'-i'} (with v = singleton(G_{s'}.position_of_gens[i'])):
- phi_0: G_0 → G_{s'} via lift_first_step(G_0, v, G_{s'})
- phi_k: G_k → G_{s'+k} via lift_one_step(G_k, gens[k], inv_ind_{k-1}, phi_{k-1}, QUT_{s'+k-1}, INJ_{s'+k}, G_{s'+k})

To compute {s-i} * {s'-i'}, we only need phi_s (not all phi_0..phi_{res_len-s'}).

## Program Interface (CURRENT)

```
./yoneda max_deg res_len s' i'
```

- `max_deg`  : maximum degree (same as mr_mot); `dir` is inferred as `"<max_deg>_"`
- `res_len`  : resolution length
- `s'`       : filtration of the Ext element being lifted
- `i'`       : generator index within filtration s' (i'-th cogenerator of G_{s'})

Output (stdout): one tab-separated line per generator {s-i}, for all s from 0 to res_len-s'
and all i, in the form:
  `{s-i}\t->\tt^a{s+s'-j}+t^b{s+s'-k}+...`  or  `{s-i}\t->\t0`

Reads: `{dir}mot_gens{k}`, `{dir}mot_maps{k}`, `{dir}ex2poly_index`,
       `{dir}mot_deltas`, `{dir}poly_exponents`, `{dir}gens_data_ctau`

Cache files: `{dir}yoneda_{s'}_{i'}_phi{k}` for k = 0..s
  (phi_k: G_k → G_{s'+k}, stored as matrix_mem<tauPoly>)

Output (stdout): one line of the form
  `{s-i} * {s'-i'} = tau^a * {s+s'-j1} + tau^b * {s+s'-j2} + ...`
  or `{s-i} * {s'-i'} = 0` if the product is zero.

## Cache Strategy

Before computing phi_k, check if `{dir}yoneda_{s'}_{i'}_phi{k}` exists.
- If yes: load from file (skip computation).
- If no:  compute from phi_{k-1} (or from v if k=0), then save.

Invariant: if phi_k file exists, then phi_0..phi_{k-1} files also exist
(ensured by always saving phi_k before advancing to k+1).

## Files

- `lift.h`      : lift_first_step, lift_one_step, recover_inv_ind, cofree_adjoint_row (unchanged)
- `yoneda.cpp`  : main program (to be rewritten with new interface)
- `mot_steenrod.h/.cpp`: provides tauPoly, motSteenrod, MotSteenrodOp, FreeMotCoMod
- `PROGRESS.md` : this file

## Implementation Steps

- [x] Write PROGRESS.md (this file)
- [x] Rewrite yoneda.cpp with new 7-argument interface
- [x] Test: compute {1-0} * {1-0} = {2-0} (h₀² at stem 0) ✓
- [x] Test: verify cache reuse (second call loads from cache, same result) ✓
- [x] Test: {1-1}*{1-1} = {2-1} (h₁² at stem 2) ✓
- [x] Test: {1-2}*{1-2} = {2-3} (h₂² at stem 6) ✓
- [x] Test: {1-0}*{1-1} = 0, {1-1}*{1-0} = 0 (commutativity, both orders agree) ✓

## Current State (COMPLETE)

- `lift.h` is complete and tested (no changes needed).
- `yoneda.cpp` has been rewritten with the new interface.
  - Old interface: 8 args (max_deg res_len s' stem weight gen_idx dir), computed all phi_0..phi_{res_len-s'}, saved as `{dir}yoneda_maps{k}`, element identified by (s', t', w', sub-index).
  - New interface: 8 args (max_deg res_len s' i' s i dir), computes only phi_0..phi_s (as needed), saves as `{dir}yoneda_{s'}_{i'}_phi{k}`, element identified by (s', i') directly.

## Key Code Notes

- `tauPoly` = int16_t; value x means tau^x; value (1<<15) = internal_zero means 0.
  Display: tau_ring_oper.output(r) → "1" / "t^k" / "0".
- `G.position_of_gens[j]` = position in G of j-th cogenerator.
- `G.find_index(pos)` = cogenerator index if pos is a cogenerator position, else -1.
- `G.generators.rank` = number of cogenerators of G.
- For a cogenerator input pos = G_s.position_of_gens[i]:
  phi_s.find(pos) has entries only at cogenerator positions of G_{s+s'}
  (because the cofree coaction of the cogenerator is trivial: 1 ⊗ e_{pos}).

## Constraint

Requires s + s' <= res_len (so that G_{s+s'} = mot_gens{s+s'} exists).
