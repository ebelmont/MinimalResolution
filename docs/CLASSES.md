# Class Catalog

Every class/struct in the repository, one row each, grouped by layer/pipeline
matching `ARCHITECTURE.md` §2. This is the flat index for "where is X
defined and what does it do" — for template parameters, base classes, and
full method references, follow the "Details" link into the relevant doc.

## Layer 2: Generic framework (`docs/FRAMEWORK.md`)

| Class | File | Purpose |
|---|---|---|
| `AbGroupOp<A>` | `algebra/2.h` | Abstract abelian-group interface |
| `RingOp<R>` | `algebra/3.h` | Abstract ring interface |
| `vectors<index,R>` | `modules/2.h` | Sparse vector — the universal row/column datatype |
| `ModuleOp<index,R>` | `modules/3.h`–`9.h` | Operations on the free module of sparse vectors |
| `matrix<R>` | `matrices/2.h`–`8.h` | Abstract matrix over a ring |
| `matrix_mem<R>` | `matrices_mem/2.h`, `4.h` | In-memory matrix backend |
| `matrix_stream<R>` | `matrices_stream.h` | Stream-backed matrix backend |
| `matrix_file<R>` | `matrices_stream.h` | File-backed matrix backend (owns a `con_fstreams`) |
| `curtis_table<ring>` | `matrices/9.h` | **Live** leading-term reduction table for one resolution step |
| `curtis_table_mem<ring>` | `matrices_mem/3.h` | In-memory `curtis_table` backend |
| `curtisTable_stream<ring>` | `curtis.h` | Disk-streamed `curtis_table` backend (with write cache) |
| `curtis_table<R>` (others.h) | `others.h` | **Dead** — unrelated second definition, not `#include`d anywhere |
| `matrix_array<R>` | `others.h` | Dead code (see above) |
| `quasi_table<R>` | `others.h` | Dead code (see above) |
| `poly<exponent_type,base_ring>` | `polynomial/1.h` | Type alias: polynomial = sparse vector keyed by exponent |
| `PolyOp<exponent_type,base_ring>` | `polynomial/1.h`, `2.h`, `6.h` | Ring structure on `poly<...>` |
| `PolyOp_Para<exponent_type,base_ring>` | `polynomial/3.h`, `4.h` | OpenMP-parallel `PolyOp` variant |
| `modules<degree_type>` | `modules/10.h` | A graded free module's shape (rank + per-generator degrees) |
| `monomial_index` | `mon_index.h`/`.cpp` | Enumerates/indexes monomials up to a max degree (non-exterior rings) |
| `con_streams` | `streams.h`/`.cpp` | Thread-safe stream wrapper |
| `con_fstreams` | `streams.h`/`.cpp` | `con_streams` owning a real file |
| `CoModule<algebroid,degree_type>` | `hopf_algebroid/2.h` | Abstract comodule interface |
| `comodule_generic<algebroid,degree_type>` | `hopf_algebroid/2.h` | Comodule with explicit coaction matrix |
| `cofree_comodule<algebroid,degree_type>` | `hopf_algebroid/3.h`, `9.h`, `11.h` | Direct sum of cofree comodule summands |
| `Hopf_Algebroid<ring,algebroid>` | `hopf_algebroid/4.h`–`13.h` | **The resolution engine** — see `FRAMEWORK.md` §5 |
| `complex<degree_type,ring>` | `hopf_algebroid/6.h` | Minimal chain-complex container |
| `SS_entry<cycle_name,ring>` | `SS.h` | One spectral-sequence-page entry |
| `SS_table<cycle_name,ring>` | `SS.h` | Generic spectral-sequence-page storage/algorithm |

Not classes but worth knowing: `inverse.h` (`inverse()`, `check_splitting()`
— partial-inverse/splitting computation from a `curtis_table`) and `lift.h`
(`lift_resolvor`, `resolution_lift` — **looks unfinished/non-compiling**, see
`FRAMEWORK.md` §3).

## Layer 3: base rings

| Class | File | Purpose |
|---|---|---|
| `Fp_Op` | `Fp.h` | Arithmetic on `F_p` (`inverse()` only implemented for p=2,3) |
| `Z3_Op` | `Z3.h` | Arithmetic on truncated 3-adic integers (`Z3` = `uint64_t` alias) |
| `Qp_Op` | `Qp.h` | Arbitrary-precision p-local rationals (GMP), `prime()` pure virtual |
| `Q3_Op` | `Qp.h` | `Qp_Op` fixed to `prime()=3` |
| `Qp_int` | `Qp.h` | `Qp_Op` variant that serializes only the truncated integer part |
| `Q3_int` | `Qp.h` | `Qp_int` + `Q3_Op` combined |

## Layer 4a: Classical Steenrod pipeline (`docs/pipelines/STEENROD.md`)

| Class | File | Purpose |
|---|---|---|
| `exponent` | `polynomial/2.h` | Packed monomial exponent vector (`typedef uint32_t`) |
| `exponentArry` | `exponents.h` | Unpacked per-variable form of `exponent` |
| `P_Op` | `steenrod.h` | Ring/module ops on `P = polynomial<Fp>` |
| `PP_Op` | `steenrod.h` | Ring/module ops on `P⊗P` |
| `Steenrod_Op` | `steenrod.h` | The Hopf algebra structure on `P` (dual Steenrod algebra, polynomial part, p=3) |
| `ComodInit` | `steenrod_init.h` | Backing comodule (in-memory `matrix_mem<P>`) |
| `SteenrodInit` | `steenrod_init.h` | Top-level driver: builds coproduct table, resolves trivial comodule |
| `FreeSteenrodCoMod` | `steenrod.h` (typedef) | `cofree_comodule<P,int>` |
| `SteenrodCoMod_generic` | `steenrod.h` (typedef) | `comodule_generic<P,int>` |

## Layer 4b: BP pipeline (`docs/pipelines/BP.md`) — the repo's main deliverable

| Class | File | Purpose |
|---|---|---|
| `BPQ_Op` | `BPQ.h` | The **rational** Hopf algebroid `BP_*⊗Q`, used to derive structure maps |
| `BPBPQ_Op` | `BPQ.h` | Ring ops on `BP_*BP⊗Q` |
| `BPBPBPQ_Op` | `BPQ.h` | Ring ops on `BP_*BP⊗BP_*BP⊗Q` |
| `BP_Op` | `BP.h` | The **integral** Hopf algebroid `(BP_*, BP_*BP)` — central class of the repo |
| `BPBP_Op` | `BP.h` | Ring ops on `BP_*BP` |
| `BPBPBP_Op` | `BP.h` | Ring ops on `BP_*BP⊗_{BP_*}BP_*BP` |
| `BPComodInit` | `BP_init.h` | Loads a generic `BP_*BP`-comodule from file |
| `BPInit` | `BP_init.h` | Driver owning every operator/table for one `mr_BP` run |
| `prim_entry` | `BPcomplex.h` | One primitive basis element |
| `primitive_data` | `BPcomplex.h` | Full primitive basis of one resolution term |
| `BPComplex` | `BPcomplex.h` | Chain complex of primitives (**not** "BP/I" — see `BP.md` §1) |
| `algNov_table` | `algNov.h` | One homological degree's algebraic-Novikov E1/E2 table |
| `algNov_tables` | `algNov.h` | Owns one `algNov_table` per homological degree |
| `Boc_table` | `Boc.h` | Bockstein-SS variant of `algNov_table` |
| `Boc_tables` | `Boc.h` | Owns one `Boc_table` per homological degree |
| `multiplication_table_entry<cycle_name>` | `multiplication.h` | One row of a multiplication table |
| `multiplication` | `multiplication.h` | Computes multiplicative structure (`h0`, `theta_i`, mult-by-3) |

## Layer 4c: Ex/Ctau pipeline (`docs/pipelines/EX_CTAU.md`) — p=2, largely copy-derived from Steenrod

| Class | File | Purpose |
|---|---|---|
| `ex_poly` | `ex_exponents.h` | Packed exponent for `F[x_i] ⊗ Λ(z_i)` (poly + exterior) |
| `ExPolyOp_Para<ring>` | `ex_exponents.h` | Multiplication respecting `z_i²=0` |
| `monomial_index` (ex variant) | `ex_index.h` | Monomial enumeration including the exterior half — **separate type from the framework's `mon_index.h` version** |
| `P` | `ctau_steenrod.h` (typedef) | `poly<ex_poly,Fp>` — the "ctau" dual object |
| `PP` | `ctau_steenrod.h` (typedef) | `P⊗P` |
| `P_Op` (ctau) | `ctau_steenrod.h` | Ring ops on `P` |
| `PP_Op` (ctau) | `ctau_steenrod.h` | Ring ops on `P⊗P` |
| `Steenrod_Op` (ctau) | `ctau_steenrod.h` | Hopf algebroid structure on `P` |
| `ComodInit` (ex) | `ex_steenrod_init.h` | Backing comodule (**disk-backed** `matrix_file<P>`, unlike classical pipeline's in-memory version) |
| `SteenrodInit` (ex) | `ex_steenrod_init.h` | Driver (disk-backed matrices throughout) |

Also: `ex_poly`/`ExPolyOp_Para` are **duplicated byte-for-byte** in
`expoly.h` (dead, unreferenced file) — see `EX_CTAU.md` §3.

## Layer 4d: Motivic pipeline (`docs/pipelines/MOTIVIC.md`) — p=2 (despite `F3` naming)

| Class | File | Purpose |
|---|---|---|
| `tauPoly` | `mot_steenrod.h` (typedef `int16_t`) | Encodes a monomial `c·τ^k` as its τ-exponent |
| `tauOper` | `mot_steenrod.h` | Arithmetic on (monomials of) `F_p[τ]` |
| `motSteenrod` | `mot_steenrod.h` (typedef) | `polynomial<tauPoly>` — the motivic dual Steenrod algebra |
| `mStmSt` | `mot_steenrod.h` (typedef) | `polynomial<motSteenrod>`, used for the coproduct |
| `MotSteenrodRingOp` | `mot_steenrod.h` | Ring/module ops on `motSteenrod` |
| `MotDegree` | `mot_steenrod.h` | Bigrading `(degree, weight)` |
| `FreeMotCoMod` | `mot_steenrod.h` (typedef) | `cofree_comodule<motSteenrod,MotDegree>` |
| `MotSteenrodOp` | `mot_steenrod.h` | The Hopf algebroid structure on `motSteenrod` |
| `motComplex` | `tao_bockstein.h` | Loads a resolution's differentials as matrices over `tauPoly` |
| `cycle_name` | `tao_bockstein.h` | Single `int` — appears to be an unused stub |
| `tau_table_entry` | `tao_bockstein.h` | One row of the τ-Bockstein table |
| `cycle_data` | `tao_bockstein.h` (typedef) | `std::map<int, vectors<matrix_index,tauPoly>>` |
| `tau_table` | `tao_bockstein.h` | The τ-Bockstein table at one resolution degree |

## Inheritance overview (selected, load-bearing relationships only)

```mermaid
classDiagram
    AbGroupOp <|-- RingOp
    RingOp <|-- ModuleOp : (A = vectors)
    ModuleOp <|-- PolyOp
    RingOp <|-- PolyOp
    PolyOp <|-- PolyOp_Para
    matrix <|-- matrix_mem
    matrix <|-- matrix_stream
    matrix_stream <|-- matrix_file
    curtis_table <|-- curtis_table_mem
    curtis_table <|-- curtisTable_stream
    CoModule <|-- comodule_generic
    CoModule <|-- cofree_comodule
    "std::vector~SS_entry~" <|-- SS_table
    SS_table <|-- algNov_table
    algNov_table <|-- Boc_table

    RingOp <|-- Fp_Op
    RingOp <|-- Z3_Op
    RingOp <|-- Qp_Op
    Qp_Op <|-- Q3_Op
    Qp_Op <|-- Qp_int
    Qp_int <|-- Q3_int

    PolyOp_Para <|-- BPBP_Op
    PolyOp_Para <|-- BPBPBP_Op
    Hopf_Algebroid <|-- BP_Op
    PolyOp_Para <|-- BP_Op

    Hopf_Algebroid <|-- Steenrod_Op_classical["Steenrod_Op (steenrod.h)"]
    Fp_Op <|-- Steenrod_Op_classical

    Hopf_Algebroid <|-- Steenrod_Op_ctau["Steenrod_Op (ctau_steenrod.h)"]
    Fp_Op <|-- Steenrod_Op_ctau

    Hopf_Algebroid <|-- MotSteenrodOp

    comodule_generic <|-- BPComodInit
    comodule_generic <|-- ComodInit_classical["ComodInit (steenrod_init.h)"]
    comodule_generic <|-- ComodInit_ex["ComodInit (ex_steenrod_init.h)"]
```

Full per-subsystem diagrams (with methods/fields shown) are in each detail
doc: `FRAMEWORK.md` §4, `BP.md` §5.
