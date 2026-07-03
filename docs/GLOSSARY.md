# Glossary: Math ↔ Code

A term index mapping mathematical concepts to the identifiers that implement
them, and vice versa. See the linked doc section for full context on each
entry — this page is a lookup table, not an explanation.

## Math concept → code identifier

| Math concept | Code identifier(s) | Where |
|---|---|---|
| Abelian group (abstract) | `AbGroupOp<A>` | `FRAMEWORK.md` §3 |
| Ring (abstract) | `RingOp<R>` | `FRAMEWORK.md` §3 |
| Module over a ring | `ModuleOp<index,R>`, elements are `vectors<index,R>` | `FRAMEWORK.md` §3 |
| Polynomial ring | `polynomial<R>` (= `vectors<exponent,R>`), `PolyOp`/`PolyOp_Para` | `FRAMEWORK.md` §3 |
| Comodule over a Hopf algebroid | `CoModule<algebroid,degree_type>`, concretely `comodule_generic` | `FRAMEWORK.md` §3 |
| Cofree (extended/induced) comodule | `cofree_comodule<algebroid,degree_type>` | `FRAMEWORK.md` §3/§5 |
| Hopf algebroid `(A, Γ)` | `Hopf_Algebroid<ring, algebroid>` (`ring`=A, `algebroid`=Γ) | `FRAMEWORK.md` §3/§5 |
| Left unit `η_L : A → Γ` | `Hopf_Algebroid::etaL` | `FRAMEWORK.md` §3, `BP.md` §3 |
| Right unit `η_R : A → Γ` | `Hopf_Algebroid::etaR` | `FRAMEWORK.md` §3, `BP.md` §3 |
| Comultiplication `Δ : Γ → Γ ⊗_A Γ` | `Hopf_Algebroid::delta` | `FRAMEWORK.md` §3 |
| Minimal relative-injective resolution, one step | `Hopf_Algebroid::resolvor` (embed + Gaussian-eliminate + quotient) | `FRAMEWORK.md` §5 |
| Embedding a comodule into a cofree comodule | `Hopf_Algebroid::embed2cofree` | `FRAMEWORK.md` §5 |
| Gröbner/leading-term reduction table for one resolution step | `curtis_table<ring>` | `FRAMEWORK.md` §3/§5 |
| "Resolution transported from a related Hopf algebroid" | `resolvor_modeled` / `pre_resolution_modeled` (a `transformer` converts the model's ring elements) | `FRAMEWORK.md` §5, `MOTIVIC.md` §2, `EX_CTAU.md` §1 |
| `BP_* = Z_(3)[v_1,v_2,...]` | `typedef polynomial<Z3> BP` | `BP.md` §1/§3 |
| `BP_*BP = BP_*[t_1,t_2,...]` | `typedef polynomial<BP> BPBP` | `BP.md` §1/§3 |
| `BP_*⊗Q` (rational BP, used to derive structure maps) | `BPQ = polynomial<Qp>`, `BPQ_Op` | `BP.md` §2/§3 |
| The invariant ideal `I = (p, v1, v2, ...)` | tracked via `algNov_table`'s cycle name `(filtration, v0-valuation, v1-exp,...,v5-exp)` | `BP.md` §1/§3 |
| Algebraic Novikov filtration / associated-graded Ext | `algNov_table` / `algNov_tables` | `BP.md` §1/§3 |
| Complex of primitives of the BP resolution | `BPComplex` / `primitive_data` / `prim_entry` (**not** "BP/I" — see `BP.md` §1) | `BP.md` §1/§3 |
| Bockstein spectral sequence (mod-3 reduction) | `Boc_table` / `Boc_tables` (= `algNov_table` with `v_valuation` counting only powers of `p`) | `BP.md` §3/§6 |
| `h0` (the class `α_1`-ish element from `η_R(v1)-η_L(v1)`) | `BP_Op::h0()` | `BP.md` §4 |
| "Theta" elements (Greek-letter/Massey-product-like classes) | `BP_Op::thetas()` | `BP.md` §4/§6 |
| 3-adic integers (truncated) | `Z3` (alias for `uint64_t`), `Z3_Op` | `BP.md` §3 |
| `p`-local rationals (exact, GMP-backed) | `Qp` struct, `Qp_Op`/`Q3_Op`/`Qp_int`/`Q3_int` | `BP.md` §3 |
| Dual Steenrod algebra `A_* = F_p[ξ_1,ξ_2,...]` (polynomial part only, `p=3`, ≤5 generators) | `P = polynomial<Fp>`, `Steenrod_Op` | `STEENROD.md` §1/§3 |
| Milnor coproduct `Δ(ξ_j) = Σ ξ_{j-i}^{p^i} ⊗ ξ_i` | `Steenrod_Op::make_delta` | `STEENROD.md` §4 |
| `Ext_{A_*}(F_p, F_p)` (input to algebraic Novikov E1) | output of `mr_st` / `SteenrodInit::resolve` | `STEENROD.md` §1/§4 |
| Dual Steenrod algebra with exterior part, `F_2[x_i] ⊗ Λ(z_i)` (speculative: `x_i = z_i²` square-root presentation) | `P = poly<ex_poly,Fp>`, `Steenrod_Op` (in `ctau_steenrod.h`) | `EX_CTAU.md` §1/§4 |
| Exterior generator `z_i` with `z_i² = 0` | `ex_poly`'s packed exterior bitmask, enforced in `ExPolyOp_Para::mon_multiply` | `EX_CTAU.md` §4 |
| Motivic dual Steenrod algebra `F_p[τ][ξ_1,ξ_2,...]` | `motSteenrod = polynomial<tauPoly>`, `MotSteenrodOp` | `MOTIVIC.md` §1/§3 |
| Motivic bidegree `(t, weight)` | `MotDegree` | `MOTIVIC.md` §3 |
| The τ variable, encoded as (a monomial in) its own tiny "ring" | `tauPoly` (packed exponent, `int16_t`), `tauOper` | `MOTIVIC.md` §1/§3 |
| τ-Bockstein spectral sequence (τ-adic filtration on the motivic resolution) | `tau_table` / `tau_table_entry`, `tao_bockstein.*` | `MOTIVIC.md` §1/§3 |
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
| `Steenrod_Op` (in `steenrod.h`) | Classical (p=3, polynomial-part-only) dual Steenrod algebra structure |
| `Steenrod_Op` (in `ctau_steenrod.h`) | Ex/ctau (p=2, poly+exterior) dual Steenrod algebra structure |
| `SteenrodInit` / `ComodInit` (either variant) | Driver: builds the coproduct table, resolves the trivial comodule |
| `MotSteenrodOp` | Motivic dual Steenrod algebra structure |
| `motComplex` | A resolution loaded as a complex over `F_p[τ]` for Bockstein analysis |

See [`CLASSES.md`](CLASSES.md) for the full class-by-class catalog with file
locations, template parameters, and base classes.
