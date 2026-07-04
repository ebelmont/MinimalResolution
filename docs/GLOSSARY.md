# Glossary: Math ↔ Code

A term index mapping mathematical concepts to the identifiers that implement
them, and vice versa. See the linked doc section for full context on each
entry — this page is a lookup table, not an explanation. Definitions/
Propositions cited by number refer to `MinimalResolution.pdf` (Guozhen Wang,
*Computations of the Adams-Novikov E2-term*, April 2020), added to the repo
after this glossary was first written — see `docs/ARCHITECTURE.md` §2 for a
guided summary and §9 for the full reference list.

## Math concept → code identifier

| Math concept | Code identifier(s) | Where |
|---|---|---|
| Abelian group (abstract) | `AbGroupOp<A>` | `FRAMEWORK.md` §3 |
| Ring (abstract) | `RingOp<R>` | `FRAMEWORK.md` §3 |
| Module over a ring | `ModuleOp<index,R>`, elements are `vectors<index,R>` | `FRAMEWORK.md` §3 |
| Polynomial ring | `polynomial<R>` (= `vectors<exponent,R>`), `PolyOp`/`PolyOp_Para` | `FRAMEWORK.md` §3 |
| Comodule over a Hopf algebroid | `CoModule<algebroid,degree_type>`, concretely `comodule_generic` | `FRAMEWORK.md` §3 |
| Cofree comodule = direct sum of copies of Γ (and shifts) | `cofree_comodule<algebroid,degree_type>` (Definition, paper §2) | `FRAMEWORK.md` §3/§5 |
| Hopf algebroid `(A, Γ)` | `Hopf_Algebroid<ring, algebroid>` (`ring`=A, `algebroid`=Γ) | `FRAMEWORK.md` §3/§5 |
| Left unit `η_L : A → Γ` | `Hopf_Algebroid::etaL` | `FRAMEWORK.md` §3, `BP.md` §3 |
| Right unit `η_R : A → Γ` | `Hopf_Algebroid::etaR` | `FRAMEWORK.md` §3, `BP.md` §3 |
| Comultiplication `Δ : Γ → Γ ⊗_A Γ` | `Hopf_Algebroid::delta` | `FRAMEWORK.md` §3 |
| `Prim(M)`, the primitive elements of a comodule `M` | `BPComplex`/`primitive_data` computes the complex of primitives of a resolution (paper §2, §7) | `BP.md` §1/§3, `FRAMEWORK.md` §5 |
| Strong injection/surjection (split injection/surjection of underlying `BP_*`-modules) | criterion `matrix<R>::gaussian`/pivot structure implements; checkable mod `I` (Remark 3.1, Proposition 3) | `FRAMEWORK.md` §5 |
| Cogenerators of a comodule `M` (an `X` with `M → BP_*BP ⊗ X` strongly injective) | `Hopf_Algebroid::adjoint`'s target (Definitions/Prop. 8) | `FRAMEWORK.md` §5 |
| Minimal (cofree) resolution — reduction mod `I` is a minimal `P`-comodule resolution | `Hopf_Algebroid::embed2cofree`'s degree-by-degree bijectivity check (Definitions 6–7) | `FRAMEWORK.md` §5 |
| Minimal relative-injective resolution, one step | `Hopf_Algebroid::resolvor` (embed + Gaussian-eliminate + quotient) | `FRAMEWORK.md` §5 |
| Embedding a comodule into a cofree comodule | `Hopf_Algebroid::embed2cofree` | `FRAMEWORK.md` §5 |
| Curtis table (entries `a` = surviving cycle, `a → b` = differential; Curtis–Goerss–Mahowald–Milgram [1]) | `curtis_table<ring>` (Proposition 9–10 give the cycle/differential correspondence) | `FRAMEWORK.md` §3/§5 |
| "Compute the minimal resolution of `M/I` first, use it as a model to lift a resolution of `M`" (paper §5, "Optimization of the process") | `resolvor_modeled` / `pre_resolution_modeled` (a `transformer` converts the model's ring elements) | `FRAMEWORK.md` §5, `ARCHITECTURE.md` §5, `MOTIVIC.md` §2, `EX_CTAU.md` §1 |
| `BP_* = Z_(3)[v_1,v_2,...]` | `typedef polynomial<Z3> BP` | `BP.md` §1/§3 |
| `BP_*BP = BP_*[t_1,t_2,...]` | `typedef polynomial<BP> BPBP` | `BP.md` §1/§3 |
| `BP_*⊗Q` (rational BP, used to derive structure maps) | `BPQ = polynomial<Qp>`, `BPQ_Op` | `BP.md` §2/§3 |
| The invariant ideal `I = (p, v1, v2, ...)` | tracked via `algNov_table`'s cycle name `(filtration, v0-valuation, v1-exp,...,v5-exp)` | `BP.md` §1/§3, `ARCHITECTURE.md` §1 |
| `P := BP_*BP/I = F_p[t1,t2,...]`, a sub-Hopf-algebra of the dual Steenrod algebra (= the whole thing, doubled degrees, at `p=2`) | `P = polynomial<Fp>`, `Steenrod_Op` (`steenrod.h`) | `STEENROD.md` §1, `ARCHITECTURE.md` §1 (paper §2) |
| Algebraic Adams-Novikov filtration: `v0^{i0}v1^{i1}···vk^{ik}·a` has filtration `i0+i1+...+ik` (`v0 = p`) | `algNov_table` / `algNov_tables` (paper §7) | `BP.md` §1/§3, `ARCHITECTURE.md` §2 |
| Bockstein spectral sequence = same construction, ordered lexicographically instead (paper §7, Remark 7.1) | `Boc_table` / `Boc_tables` (= `algNov_table` with `v_valuation` counting only powers of `p`) | `BP.md` §3/§6, `ARCHITECTURE.md` §2 |
| Complex of primitives of the BP resolution | `BPComplex` / `primitive_data` / `prim_entry` (**not** "BP/I" — see `BP.md` §1) | `BP.md` §1/§3 |
| Atiyah-Hirzebruch differential = multiplication by the extension class `h` of a two-step filtration (Proposition 12) | `BP_Op::h0()`, `BP_Op::thetas()`, `multiplication.cpp` | `BP.md` §4/§6, `ARCHITECTURE.md` §2 |
| 3-adic integers (truncated; `p=2` original documented as `Z2.h` in paper §9.14) | `Z3` (alias for `uint64_t`), `Z3_Op` | `BP.md` §3/§6 |
| `p`-local rationals (exact, GMP-backed; original `Qp`/Hazewinkel-generator construction, paper §9.12-13) | `Qp` struct, `Qp_Op`/`Q3_Op`/`Qp_int`/`Q3_int` | `BP.md` §3 |
| Milnor coproduct `Δ(ξ_j) = Σ ξ_{j-i}^{p^i} ⊗ ξ_i` | `Steenrod_Op::make_delta` | `STEENROD.md` §4 |
| `Ext_P(F_3, F_3)` (the model resolution `mr_BP` lifts, paper §5) | output of `mr_st` / `SteenrodInit::resolve` | `STEENROD.md` §1/§4, `ARCHITECTURE.md` §5 |
| Mod-2 dual Steenrod algebra via square-root presentation `x_i = z_i²` (best-effort inference, not paper-confirmed; a *different* ring from `steenrod.h`'s `P`, used by the independent p=2 cross-check) | `P = poly<ex_poly,Fp>`, `Steenrod_Op` (in `ctau_steenrod.h`) | `EX_CTAU.md` §1/§4 |
| Exterior generator `z_i` with `z_i² = 0` | `ex_poly`'s packed exterior bitmask, enforced in `ExPolyOp_Para::mon_multiply` | `EX_CTAU.md` §4 |
| Motivic dual Steenrod algebra `F_p[τ][ξ_1,ξ_2,...]` | `motSteenrod = polynomial<tauPoly>`, `MotSteenrodOp` | `MOTIVIC.md` §1/§3 |
| Motivic bidegree `(t, weight)` | `MotDegree` | `MOTIVIC.md` §3 |
| The τ variable, encoded as (a monomial in) its own tiny "ring" | `tauPoly` (packed exponent, `int16_t`), `tauOper` | `MOTIVIC.md` §1/§3 |
| τ-Bockstein spectral sequence (τ-adic filtration on the motivic resolution; converges to the algebraic Novikov SS via Gheorghe–Wang–Xu [2]'s cofiber-of-τ isomorphism) | `tau_table` / `tau_table_entry`, `tao_bockstein.*` | `MOTIVIC.md` §1/§3, `ARCHITECTURE.md` §6 |
| Multiplication by `h_i` on the E2/associated-graded page | `MotSteenrodOp::hi(i)`, `multiplication_table` | `MOTIVIC.md` §1/§4 |
| Monomial enumeration/indexing scheme | `monomial_index` (two independent implementations: `mon_index.h` for non-exterior rings, `ex_index.h` for exterior-augmented rings) | `FRAMEWORK.md` §3, `EX_CTAU.md` §3 |

## Code identifier → math concept (reverse index)

For quickly answering "what is this class actually *for*" when reading code
top-down.

| Identifier | Meaning |
|---|---|
| `AbGroupOp<A>` | Abstract abelian group interface |
| `RingOp<R>` | Abstract ring interface (extends `AbGroupOp`) |
| `ModuleOp<index,R>` | Operations on free-module sparse vectors |
| `vectors<index,R>` | A sparse vector (the universal "row/column" datatype) |
| `matrix<R>` / `matrix_mem` / `matrix_stream` / `matrix_file` | Matrix over a ring; three storage backends |
| `polynomial<R>` | A polynomial ring, represented as `vectors<exponent,R>` |
| `PolyOp` / `PolyOp_Para` | Ring structure on `polynomial<R>` (serial / OpenMP-parallel) |
| `CoModule` / `comodule_generic` / `cofree_comodule` | Comodule interface / explicit-coaction impl / cofree-summand impl |
| `Hopf_Algebroid<ring,algebroid>` | **The resolution engine** — structure maps + generic resolution algorithm |
| `curtis_table` (+ `_mem`, `curtisTable_stream`) | Leading-term reduction table for one resolution step |
| `SS_table` / `SS_entry` | Generic spectral-sequence-page storage, subclassed by `algNov_table` |
| `monomial_index` | Enumerates/indexes monomials of a graded polynomial ring |
| `con_streams` / `con_fstreams` | Thread-safe stream/file wrapper |
| `Fp_Op` | Arithmetic on `F_p` |
| `Z3_Op` | Arithmetic on (truncated) 3-adic integers |
| `Qp_Op` / `Q3_Op` / `Qp_int` / `Q3_int` | Arithmetic on p-local rationals (GMP-backed) |
| `BP_Op` / `BPBP_Op` / `BPBPBP_Op` | The `BP_*`/`BP_*BP` Hopf algebroid structure and ring ops |
| `BPQ_Op` / `BPBPQ_Op` / `BPBPBPQ_Op` | The rational (`⊗Q`) version, used to derive `BP_Op`'s structure maps |
| `BPComplex` / `primitive_data` / `prim_entry` | Complex of primitives of the BP resolution |
| `algNov_table` / `algNov_tables` | Algebraic-Novikov E1/E2-page table(s) |
| `Boc_table` / `Boc_tables` | Bockstein-SS variant of the above |
| `multiplication` / `multiplication_table_entry` | Multiplicative structure on algNov/Boc tables |
| `Steenrod_Op` (in `steenrod.h`) | Structure on `P = BP_*BP/I` (p=3) — the model resolution for `mr_BP` |
| `Steenrod_Op` (in `ctau_steenrod.h`) | Structure on a *different* ring, `poly<ex_poly,Fp>` (p=2, poly+exterior) — the independent motivic-cross-check pipeline, unrelated to `steenrod.h`'s `P` despite the shared class name |
| `SteenrodInit` / `ComodInit` (either variant) | Driver: builds the coproduct table, resolves the trivial comodule |
| `MotSteenrodOp` | Motivic dual Steenrod algebra structure |
| `motComplex` | A resolution loaded as a complex over `F_p[τ]` for Bockstein analysis |

See [`CLASSES.md`](CLASSES.md) for the full class-by-class catalog with file
locations, template parameters, and base classes.
