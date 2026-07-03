# Generic Framework Layer

This document covers the math-agnostic template infrastructure in
`MinimalResolution` — the layer that the Steenrod, BP, and motivic pipelines
are all built on top of. It does **not** cover any concrete Hopf algebroid,
comodule, or spectral-sequence pipeline; those are documented separately.

Source files covered: `algebra.h` (+`algebra/1-3.h`), `matrices.h`
(+`matrices/1-11.h`), `matrices_mem.h` (+`matrices_mem/1-4.h`),
`matrices_stream.h`, `polynomial.h` (+`polynomial/1-6.h`), `modules.h`
(+`modules/1-10.h`), `hopf_algebroid.h` (+`hopf_algebroid/1-13.h`),
`mon_index.h/.cpp`, `streams.h`, `others.h`, `inverse.h`, `lift.h`, `curtis.h`,
`SS.h`.

## 1. Overview

This layer supplies the reusable machinery that every concrete pipeline
(Steenrod algebra at p=3, BP, motivic) instantiates with its own ring and
comodule structure: generic abelian-group/ring interfaces (`algebra.h`),
generic sparse-vector and matrix types over an arbitrary ring
(`modules.h`, `matrices.h`, `matrices_mem.h`, `matrices_stream.h`),
a generic polynomial ring builder (`polynomial.h`), and — the centerpiece —
`Hopf_Algebroid<ring,algebroid>` (`hopf_algebroid.h`), which implements the
comodule-resolution algorithm (embedding into a cofree comodule, Gaussian
elimination, quotienting) purely in terms of abstract ring/module operations
and a caller-supplied Hopf-algebroid structure (`etaL`, `etaR`, `delta`, …).
Separating this from the concrete math means the same resolution engine,
matrix backend (in-memory or disk-streamed, for handling resolutions too
large to fit in RAM), and polynomial arithmetic are shared verbatim across
completely different Hopf algebroids (dual Steenrod algebra, BP*BP, motivic
variants) — only the small algebroid-specific classes (structure maps,
concrete comodules) differ per pipeline.

## 2. The numbered-fragment pattern

Several "big" classes/headers are not written as a single file. Instead the
top-level header (e.g. `matrices.h`) is a pure `#include` manifest:

```cpp
// matrices.h
#pragma once
#include"matrices/1.h"
#include"matrices/2.h"
#include"matrices/3.h"
...
#include"matrices/9.h"
```

The **first** numbered fragment (`matrices/1.h`, `modules/1.h`,
`hopf_algebroid/1.h`, `algebra/1.h`, `matrices_mem/1.h`) contains only the
`#pragma once` guard and the `#include` of dependency headers (e.g.
`matrices/1.h` pulls in `modules.h`). The **second** fragment (`matrices/2.h`,
`modules/2.h`, `hopf_algebroid/2.h`...) contains the actual class *definition*
with its full body (data members + declarations), often declaring some methods
`virtual`/pure and leaving others merely declared (not defined). Every
**subsequent** fragment (`matrices/3.h` … `matrices/9.h`, etc.) is **not** a
class body at all — it is a sequence of out-of-line `template<typename R> ...
matrix<R>::method(...) { ... }` definitions for the methods declared in
fragment 2, grouped by theme (IO in one file, composition operators in
another, quotient-construction in another, etc.). Because C++ template
member functions can be defined anywhere after the class declaration is
visible, this works as long as fragment 2 is included before the later
fragments — which the manifest header guarantees.

Concretely, for `matrices.h`:
- `matrices/1.h` — dependencies only (`#include "modules.h"`, `typedef ... matrix_index`).
- `matrices/2.h` — the `matrix<R>` class body (pure-virtual interface + a few inline helper methods like `filter`, `del_cols`, `construct`, `merge`, `direct_sum`), plus the static member definition `ModuleOp<matrix_index,R> *matrix<R>::moduleOper;` at the bottom.
- `matrices/3.h` — `set2unit`, `set2zero`, `equal` definitions.
- `matrices/4.h` — `output`, `save`, `load`, `load_modify` definitions.
- `matrices/5.h` — `maps_to`, `compose`, and free-function `maps_to`/`maps_to_p` template overloads (parallel map-to-generators).
- `matrices/6.h` — `quot_index`, `make_quotient` definitions.
- `matrices/7.h` — `row_operation`, `filter` (two overloads) definitions.
- `matrices/8.h` — the free function `sort_deg` (degree-based stable sort of `(row,col)` pairs).
- `matrices/9.h` — **a second, independent class**, `curtis_table<ring>` (full definition, not a fragment of `matrix`!). See §3 and the note on redundancy.
- `matrices/11.h` — empty (0 bytes of content, just a blank line). `matrices/10.h` does not exist — the numbering has a gap.

The same three-part shape (deps-only fragment → class-body fragment → N
method-body fragments) recurs in `modules/*.h`, `hopf_algebroid/*.h`, and
`matrices_mem/*.h`, though `hopf_algebroid/*.h` interleaves multiple classes
(`CoModule`, `comodule_generic`, `cofree_comodule`, `Hopf_Algebroid`, and even
a `complex` struct tacked onto the end of `hopf_algebroid/6.h`) rather than
splitting one class across all fragments — see §3 for exactly which fragment
holds which class.

**Reader takeaway:** if you open a fragment file like `matrices/5.h` or
`hopf_algebroid/10.h` and see a bare `template<typename R> return_type
matrix<R>::method(...) { ... }` with no enclosing `class { ... }`, this is
normal — it is an out-of-line definition for a class declared in an earlier,
lower-numbered fragment of the same directory. Grep the whole numbered
directory (not just the file you're looking at) to find the class body and
all its method definitions.

## 3. Class catalog

| Class/Struct | File(s) | Template Parameters | Base Classes | Purpose |
|---|---|---|---|---|
| `AbGroupOp` | `algebra/2.h` | `<A>` | (none, pure interface) | Abstract interface for abelian-group operations on an opaque element type `A`: `add`, `zero`, `isZero`, `minus`, a destructive/move `add`, plus `output`/`save`/`load` for serialization. Every ring, module, and vector-of-vectors type in the codebase is driven through some concrete subclass of this. |
| `RingOp` | `algebra/3.h` | `<R>` | `virtual public AbGroupOp<R>` | Adds `multiply`, `unit(int)` (ring hom from ℤ), `invertible`, `inverse`, and a default `power` (binary exponentiation) on top of the abelian-group interface. This is the ring-level contract that `ModuleOp`, `PolyOp`, and `Hopf_Algebroid` all consume via a `RingOp<R>*`. |
| `vectors<index,R>` | `modules/2.h` | `<index, R>` | (none — plain data class) | A sparse vector: `dataArray` is a sorted `std::vector<term>` where `term{ind, coeficient}`. Provides `direct_sum`/`un_direct_sum` (block sum / split, used to build direct sums of comodules), dense/sparse conversion (`toDense`/`deDense`/`add2Dense`). This is the fundamental "column"/"row" datatype used everywhere: rows of matrices, elements of modules, coefficients of a comodule coaction. |
| `ModuleOp<index,R>` | `modules/3.h`–`modules/9.h` | `<index, R>` | `virtual public AbGroupOp<vectors<index,R>>` | Operations on the free `R`-module of sparse vectors `vectors<index,R>`: `add`, `scalor_mult`, `singleton`, `re_index`/`filtered_reindex`, `filter`, `termwise_operation`, `component` (binary-search lookup by index, since `dataArray` is sorted), `next_non_trvial_entry`, `first_invertible[_index]`, `last_entry`, and a parallel `sum` over an array of vectors. Holds a `RingOp<R> *ringOper`. Nearly every generic algorithm (Gaussian elimination, resolution) is expressed in terms of a `ModuleOp` instance rather than raw arithmetic. |
| `matrix<R>` | `matrices/2.h`–`matrices/8.h` | `<R>` | (abstract; concrete backends below implement it) | A matrix over ring `R`, represented as an array of row vectors (`vectors<matrix_index,R>`), addressed by `matrix_index` (`uint32_t`). Pure-virtual storage primitives (`clear`, `find`, `insert`, `set_rank`, `update_all`, `gaussian`, `del_and_gaussian`) are implemented by backends; everything else (`compose`, `maps_to`, `quot_index`, `make_quotient`, `direct_sum`, IO) is written generically against those primitives. Holds a `static ModuleOp<matrix_index,R> *moduleOper` shared across all instances of `matrix<R>`. |
| `matrix_mem<R>` | `matrices_mem/2.h`, `matrices_mem/4.h` | `<R>` | `public matrix<R>` | In-memory backend: rows stored in a `std::vector<vectors<matrix_index,R>>`. Implements `gaussian`/`row_reduction` (parallel row reduction via OpenMP) and helper `unify`/`take_away` used by Gaussian elimination. |
| `matrix_stream<R>` | `matrices_stream.h` | `<R>` | `public matrix<R>` | Backend that stores each row at a byte offset in a `con_streams*` (concurrent stream wrapper), with an in-memory index (`datapos`) of stream positions. `gaussian`/`del_and_gaussian` are implemented by round-tripping through a temporary `matrix_mem<R>` (load, reduce in memory, write back) rather than operating on the stream directly. |
| `matrix_file<R>` | `matrices_stream.h` | `<R>` | `public matrix_stream<R>` | `matrix_stream<R>` bound to an owned `con_fstreams` (a real file), for resolutions too large to keep in RAM. |
| `curtis_table<ring>` | `matrices/9.h` | `<ring>` | (abstract) | **The** live/used curtis-table interface — see §5 for what it stores. Abstract storage primitives (`clear`, `is_member`, `search`/`search_ref`, `num_entries`, `insert`, `run_through`) plus generic algorithms built on them: `symplify_to_led`/`simplify_to_led` (reduce a vector against the table, tracking a homotopy), `cycle_matrix` (extract the matrix of stored cycles/tags), `save`/`load`. Holds `static ModuleOp<matrix_index,ring> *ModOper`. |
| `curtis_table_mem<ring>` | `matrices_mem/3.h` | `<ring>` | `public curtis_table<ring>` | In-memory backend for `curtis_table`: entries kept in a `std::map<matrix_index, entry>`. |
| `curtisTable_stream<ring>` | `curtis.h` | `<ring>` | `public curtis_table<ring>` | Disk/stream-backed `curtis_table`: keeps a small in-memory write cache (`curtis_table_mem<ring> cached`) plus a `std::map<matrix_index,std::streampos>` index into a `con_streams`; `flush()` drains the cache to the stream. Used for resolutions whose Curtis table is too large for memory. |
| `curtis_table<R>` **(others.h)** | `others.h` (lines ~159–283) | `<R>` | none listed (independent hierarchy) | **A second, unrelated, unused class of the same name.** See "Redundancy" note below. |
| `matrix_array<R>` | `others.h` | `<R>` | `virtual public matrix<R>` (of the `others.h` `matrix<R>`, not the real one) | Dead code — see redundancy note. |
| `quasi_table<R>` | `others.h` | `<R>` | `public curtis_table<R>` (the `others.h` one) | Dead code — see redundancy note. |
| `poly<exponent_type,base_ring>` | `polynomial/1.h` | alias | `= vectors<exponent_type,base_ring>` | Type alias: a polynomial is just a sparse vector keyed by exponent. |
| `PolyOp<exponent_type,base_ring>` | `polynomial/1.h`, `.../2.h`, `.../6.h` | `<exponent_type, base_ring>` | `virtual public ModuleOp<exponent_type,base_ring>`, `virtual public RingOp<poly<...>>` | Turns the additive `ModuleOp` on sparse vectors into a full ring by adding `multiply` (divide-and-conquer polynomial multiplication via repeated `mon_multiply` + `add`), `unit`, `invertible`/`inverse` (only degree-0 invertible constants), `constant`, `monomial`. |
| `PolyOp_Para<exponent_type,base_ring>` | `polynomial/3.h`, `.../4.h` | `<exponent_type, base_ring>` | `public PolyOp<exponent_type,base_ring>` | Parallel-algorithm override of `mon_multiply`/`multiply` (OpenMP `parallel for` / recursive `sum`), otherwise identical semantics to `PolyOp`. |
| `modules<degree_type>` | `modules/10.h` | `<degree_type>` | (plain data class) | A graded free module's shape only: `rank` plus a `degree` array (one degree per generator). Used as the "generators" bookkeeping structure inside `cofree_comodule`. Has `direct_sum`, `save`/`load`, `output`. |
| `monomial_index` | `mon_index.h`, `mon_index.cpp` | (non-template) | none | Enumerates/indexes monomials (type `exponent`, from `exponents.h`) up to a max degree: builds `mon_array` (all monomials sorted by degree), `mon_index` (monomial → `matrix_index`), `ranksBelow` (generator count below each degree). Provides `substitution_table` (evaluate all monomials at given values of generators, streaming the results out) and `poly2vec` (convert a `polynomial<base_ring>` into a `vectors<matrix_index,base_ring>` using the index). Math-agnostic in the sense that it works for any graded polynomial ring on generators with `exponent`-encoded exponents (used concretely for the dual Steenrod algebra generators, but the machinery itself doesn't know that). |
| `con_streams` | `streams.h`, `streams.cpp` | (non-template) | none | Thread-safe wrapper around a `std::iostream*`: `access`/`read`/`write` all take a `std::lock_guard` on an internal `std::mutex` before touching the stream, so multiple threads can safely read/write to disjoint regions of one stream (positions returned by `write` are handed back to the caller for later `read`). Used by `matrix_stream`/`matrix_file` and `curtisTable_stream`. |
| `con_fstreams` | `streams.h`, `streams.cpp` | (non-template) | `public con_streams` | `con_streams` that owns a real `std::fstream` opened on a given filename; adds `fclear()` to truncate/reopen. |
| `CoModule<algebroid,degree_type>` | `hopf_algebroid/2.h` | `<algebroid, degree_type>` | (abstract) | The comodule interface consumed by `Hopf_Algebroid`: `rank()`, `coaction(int i)` (returns the coaction of the i-th base generator, as a `vectors<matrix_index,algebroid>`), `degree(int i)`. Assumes the comodule is free as a module over the base ring with a specified basis. |
| `comodule_generic<algebroid,degree_type>` | `hopf_algebroid/2.h` | `<algebroid, degree_type>` | `virtual public CoModule<algebroid,degree_type>` | The generic (non-cofree) implementation: coaction stored explicitly as a `matrix<algebroid> *coaction_matrix`, generator data in a `modules<degree_type> base_module`. This is the type used to hold "the comodule currently being resolved" at each stage of the resolution (it gets overwritten in place, stage by stage, to become successive connected covers/quotients — see §5). |
| `cofree_comodule<algebroid,degree_type>` | `hopf_algebroid/3.h`, `hopf_algebroid/9.h`, `hopf_algebroid/11.h` | `<algebroid, degree_type>` | `virtual public CoModule<algebroid,degree_type>` | Represents a direct sum of cofree (extended/induced) comodules, one per co-generator, without materializing the coaction matrix explicitly. Static class-wide data (`cofree_coaction`, `modOpers`, `cofree_degree` — set once via `Hopf_Algebroid::init_cofree_data`) describes the coaction/degree of a single cofree summand co-generated in degree 0; instance data (`generators`, `position_of_gens`, `total_rank`) describes how many summands and where each starts. `coaction(i)`/`degree(i)` look up which summand `i` falls in (`findPos`, binary search over `position_of_gens`) and re-index/shift accordingly. Also provides `multiply_using_table` (apply a precomputed multiplication-table matrix to compute the algebroid action on generator `i`) and `direct_sum`. |
| `Hopf_Algebroid<ring,algebroid>` | `hopf_algebroid/4.h`, `.../5.h`, `.../6.h`, `.../7.h`, `.../8.h`, `.../9.h`, `.../10.h`, `.../12.h`, `.../13.h` | `<ring, algebroid>` | (abstract; holds `RingOp`/`ModuleOp` pointers, no base class) | **The central class.** Bundles ring/module operations for the base ring and for the algebroid, and declares the Hopf-algebroid structure maps as pure virtuals (`etaL`, `etaR`, `algebroid2vector`, `vector2algebroid`, `delta`, `ranksBelowDeg`) that a concrete pipeline must supply. On top of these it implements, generically, the entire resolution algorithm: `adjoint` (built from projection + structure maps), `embed2cofree`/`embed2cofree_modeled`, `resolvor`/`resolvor_modeled` (one resolution step), `pre_resolution_tab`/`pre_resolution_modeled` (iterate resolution steps to a fixed length, with file-backed persistence), `quotient`/`quotient_p`, `resolution` (splice per-step short exact sequences into one chain complex), `make_multiplication_table`. See §5 for the resolution semantics. |
| `complex<degree_type,ring>` | `hopf_algebroid/6.h` (tail of file, after `resolution`) | `<degree_type, ring>` | (plain data class) | A chain complex: parallel arrays `terms` (one `modules<degree_type>` per stage) and `maps` (one `matrix<ring>*` per differential). Minimal container, not used elsewhere in the files read. |
| `SS_entry<cycle_name,ring>` | `SS.h` | `<cycle_name, ring>` | (plain data class) | One entry of a spectral-sequence page: leading-term names `tag`/`cycle` (of type `cycle_name`, chosen by the concrete pipeline) plus the corresponding full vectors `full_tag`/`full_cycle` (type `vectors<matrix_index,ring>`). |
| `SS_table<cycle_name,ring>` | `SS.h` | `<cycle_name, ring>` | `public std::vector<SS_entry<cycle_name,ring>>` (abstract — several pure virtuals) | Generic storage/algorithm for one page of a spectral sequence built from a Curtis-table-like reduction process: indexes entries by `tag`/`cycle` name, and provides `simplify` (reduce a cycle modulo known boundaries, recording a homotopy — directly analogous to `curtis_table::symplify_to_led`), `make_table` (build a table from a candidate pool of tag names, ordered by filtration, calling a caller-supplied differential matrix `M`), and `name_of_cycle` (classify a vector's leading terms as tags/boundaries/cycles, optionally consulting a `next_table` for the following page). Concrete pipelines (BP's `algNov_table` — see the BP-pipeline doc) subclass this and supply `filtration`, `naming`, `leading_term`, `tagged`, `invalid`, `get_tag`, `cycle_pot`, and IO. Generic in the sense that it doesn't know what a spectral sequence *of* — it operates purely on `vectors<matrix_index,ring>` and caller-provided naming/filtration functions. |

### Note: the `curtis_table` / `matrix_array` / `quasi_table` redundancy

`others.h` (397 lines) contains **a second, independent, and incompatible**
definition of `curtis_table<R>`, `matrix_array<R>`, `matrix_array_data<R>`,
`quasi_table<R>`, plus free functions (`BP_Op::Q2Z` overloads,
`make_tables_gens`, a stray `matrix<R>::direct_sum`/`del_cols`). This code:

- Is **not `#include`d anywhere** in the repository (`grep` for `others.h`
  across all `.h`/`.cpp` files returns no hits).
- Uses namespaces (`Modules::`, `Algebra::`, `MapTable::`, `Paralell::`) that
  do not exist anywhere else in the codebase — the rest of the code is
  entirely un-namespaced.
- References a header `filedmatrices.h` (line 379) that **does not exist**
  in the repository.
- Has a `curtis_table` constructor signature
  (`curtis_table(ModuleOp*, MapTable::data_base*)`, `others.h:170`) completely
  incompatible with the real `curtis_table<ring>` in `matrices/9.h` (which has
  no constructor and a `static ModuleOp* ModOper`).

This is dead legacy code from an earlier refactor (probably the point where
the codebase moved from a namespaced `Modules::`/`Algebra::` design to the
flat, unnamespaced design seen in `algebra.h`/`modules.h`/`matrices.h`), left
in the tree but never wired into any build. **`others.h` should not be treated
as part of the live framework** — the framework's actual, used
`curtis_table<ring>` is the one in `matrices/9.h`, with backends
`curtis_table_mem<ring>` (`matrices_mem/3.h`, in-memory) and
`curtisTable_stream<ring>` (`curtis.h`, disk-streamed with an in-memory write
cache that gets `flush()`ed). Concrete pipelines pick whichever backend fits
their memory budget (e.g. `BP_init.h:58-59` uses `curtis_table_mem<F3>`;
`ex_steenrod_init.h:34-35` keeps both a `curtis_table<Fp>*` pointer and a
`curtis_table_mem<Fp>` vector side by side).

`inverse.h` and `lift.h` are small, genuinely-used utility headers, not
fragments of a bigger class:
- `inverse.h` computes a partial inverse of an injective matrix from its
  Curtis table (`inverse()`, `inverse.h:6`), plus a `resolution()` overload
  and `check_splitting()` sanity-checker. `resolution_length`-scale batch
  version at `inverse.h:36`.
- `lift.h` implements lifting a map between comodules to a lift of one step of
  their resolutions (`lift_resolvor`, `resolution_lift`) — this is dead-ish
  code too in parts (e.g. `adjoint` at `lift.h:19` references
  `F.pos`/`HA_oper->` with `->` on a reference, and `HA_oper->algebroidModuleOper->component(i)`
  passes only one argument to a two-argument `component` — this file does not
  look like it compiles as written; treat it as unfinished/unverified).

## 4. Key relationships

```mermaid
classDiagram
    class AbGroupOp~A~{
        <<interface>>
        +add(A,A) A
        +zero() A
        +isZero(A) bool
        +minus(A) A
        +save(A, stream)
        +load(stream) A
    }
    class RingOp~R~{
        <<interface>>
        +multiply(R,R) R
        +unit(int) R
        +invertible(R) bool
        +inverse(R) R
        +power(R,n) R
    }
    AbGroupOp <|-- RingOp : R = R

    class ModuleOp~index,R~{
        +ringOper RingOp~R~*
        +add(vec,vec) vec
        +scalor_mult(R,vec) vec
        +singleton(index,R) vec
        +component(index,vec) R
        +re_index(rule,vec) vec
    }
    AbGroupOp <|-- ModuleOp : A = vectors~index,R~

    class PolyOp~exp,R~{
        +multiply(poly,poly) poly
        +unit(int) poly
        +monomial(exp,R) poly
    }
    ModuleOp <|-- PolyOp : index=exp
    RingOp <|-- PolyOp : R = poly~exp,R~
    class PolyOp_Para~exp,R~{
        +mon_multiply(...) poly
        +multiply(poly,poly) poly
    }
    PolyOp <|-- PolyOp_Para

    class vectors~index,R~{
        +dataArray vector~term~
        +direct_sum(vec, rank)
        +toDense(rank) vector~R~
    }
    ModuleOp ..> vectors : operates on

    class matrix~R~{
        <<abstract>>
        +rank unsigned
        +moduleOper ModuleOp~matrix_index,R~$
        +find(i) vec
        +insert(i,vec)
        +gaussian(rowcols)*
        +compose(...)
        +quot_index(...)$
        +make_quotient(...)
    }
    matrix ..> ModuleOp : static moduleOper
    class matrix_mem~R~{
        +gaussian(rowcols)
        +row_reduction(row,col)
    }
    matrix <|-- matrix_mem
    class matrix_stream~R~{
        +datas con_streams*
        +find(i) vec
        +gaussian(rowcols)
    }
    matrix <|-- matrix_stream
    class matrix_file~R~{
        +files con_fstreams
    }
    matrix_stream <|-- matrix_file

    class curtis_table~ring~{
        <<abstract>>
        +ModOper ModuleOp~matrix_index,ring~$
        +entry cycle,tag,full_cycle,full_tag
        +symplify_to_led(vec,...) matrix_index
        +cycle_matrix(...) vector
        +insert(...)*
        +search(i)*
    }
    curtis_table ..> ModuleOp : static ModOper
    curtis_table ..> matrix : cycle_matrix builds matrix
    class curtis_table_mem~ring~{
        +data map~matrix_index,entry~
    }
    curtis_table <|-- curtis_table_mem
    class curtisTable_stream~ring~{
        +cached curtis_table_mem~ring~
        +stream_data con_streams
    }
    curtis_table <|-- curtisTable_stream

    class CoModule~algebroid,degree_type~{
        <<abstract>>
        +rank() int*
        +coaction(i) vec*
        +degree(i) degree_type*
    }
    class comodule_generic~algebroid,degree_type~{
        +coaction_matrix matrix~algebroid~*
        +base_module modules~degree_type~
    }
    CoModule <|-- comodule_generic
    class cofree_comodule~algebroid,degree_type~{
        +cofree_coaction fn$
        +modOpers ModuleOp~matrix_index,algebroid~*$
        +generators modules~degree_type~
        +position_of_gens vector~uint32~
        +findPos(n) unsigned
    }
    CoModule <|-- cofree_comodule

    class Hopf_Algebroid~ring,algebroid~{
        <<abstract>>
        +ringOper RingOp~ring~*
        +algebroidRingOper RingOp~algebroid~*
        +moduleOper ModuleOp~matrix_index,ring~*
        +algebroidModuleOper ModuleOp~matrix_index,algebroid~*
        +etaL(ring) algebroid*
        +etaR(ring) algebroid*
        +delta(i) vec*
        +embed2cofree(...) cofree_comodule
        +resolvor(...) cofree_comodule
        +pre_resolution_tab(...)
        +quotient(...)
    }
    Hopf_Algebroid --> RingOp : uses (base + algebroid)
    Hopf_Algebroid --> ModuleOp : uses (base + algebroid)
    Hopf_Algebroid ..> CoModule : resolves
    Hopf_Algebroid ..> cofree_comodule : embed2cofree produces
    Hopf_Algebroid ..> curtis_table : resolvor/embed2cofree use as scratch structure
    Hopf_Algebroid ..> matrix : inj/quot/indj matrices

    class SS_entry~cycle_name,ring~{
        +tag cycle_name
        +cycle cycle_name
        +full_tag vec
        +full_cycle vec
    }
    class SS_table~cycle_name,ring~{
        <<abstract>>
        +Modop ModuleOp~matrix_index,ring~*
        +simplify(...) cycle_type
        +make_table(pric, M, T)
    }
    SS_table --|> "std::vector~SS_entry~" : extends
    SS_table ..> ModuleOp : uses
    SS_table ..> curtis_table : structurally analogous (simplify ~ symplify_to_led)
```

## 5. The resolution algorithm in this framework's terms

The framework computes a **minimal relative-injective (cofree) resolution**
of a comodule `X` over a Hopf algebroid, one step at a time. Everything below
is read directly off `hopf_algebroid/*.h`; the *mathematical* justification
(why this specific construction gives a minimal resolution, why the "cofree"
comodules used are the relative injectives for this category, etc.) is
external (`MinimalResolution.pdf`, not present in this repo) — flagged with
`TODO(math)` below wherever the code doesn't itself explain the "why".

### `curtis_table<ring>`: what it stores

A `curtis_table` (`matrices/9.h:1-263`) stores **one Gröbner-basis-like
reduction table for one resolution step**: a set of `entry { cycle, tag,
full_cycle, full_tag }` records where

- `cycle` is the index of a basis element of the *source* comodule `X` (or,
  during `embed2cofree`, of the cofree comodule being built) that is a
  "leading term",
- `tag` is the index of the newly introduced generator (in the cofree
  comodule / the next term of the resolution) that this cycle corresponds to,
- `full_cycle` is the entire vector (over the base ring) expressing how that
  generator's image cancels down to leading term `cycle`,
- `full_tag` is the corresponding "homotopy" vector recording the coefficient
  data needed to invert/split the map later (consumed by `inverse.h`'s
  `inverse()` to build a partial inverse/splitting map).

`symplify_to_led`/`simplify_to_led` (`matrices/9.h:34-140`) is the core
reduction routine: given a vector `x` (a would-be cycle) and a `homotopy`
accumulator, it walks term-by-term; if the current leading term's index is
already a `cycle` in the table, it subtracts off `coefficient * full_cycle`
(accumulating the matching multiple of `full_tag` into `homotopy`) and
recurses on the remainder; if the leading term is **invertible** in the ring
and not in the table, that index becomes the new leading term returned
(`Boundary` — actually `matrix_index Boundary = -1`, `matrices/9.h:31` — is
returned only once the vector reduces to zero). This is literally
row-reduction against a partial basis, done term-by-term rather than as a
single matrix operation, and is what makes the resolution "minimal": a
generator is only introduced (in `embed2cofree`, see below) when reduction
against everything already known fails to reach the zero vector.

`cycle_matrix` (`matrices/9.h:159-204`) turns the table into two matrices: one
whose row `tag` is `full_cycle` (used as an injective/embedding matrix into
the cofree comodule) and one whose row `tag` is `full_tag`.

### `cofree_comodule<algebroid,degree_type>`: what it represents

A cofree comodule here is a direct sum of **cofree (extended/induced)
comodules**, one per co-generator — i.e. comodules of the shape
"algebroid co-generated freely in one degree", assembled via the static class
data `cofree_coaction`/`cofree_degree`/`modOpers` set once via
`Hopf_Algebroid::init_cofree_data` (`hopf_algebroid/10.h:14-25`), which wires
`cofree_coaction` to the Hopf algebroid's own comultiplication `delta`
(`hopf_algebroid/4.h:24`). Each `cofree_comodule` instance just records how
many summands (`generators`, a `modules<degree_type>`), where each summand's
block starts (`position_of_gens`), and the total underlying rank
(`total_rank`). `coaction(i)`/`degree(i)` locate which summand index `i` falls
in (`findPos`, binary search, `hopf_algebroid/11.h:40-54`) and shift/reindex
the static cofree data accordingly (`hopf_algebroid/11.h:62-75`).

**TODO(math):** the code assumes (via `adjoint`, see below) that
`Hom_{comod}(X, cofree-on-degree-d-generator) ≅ (X in degree d, as a
base-ring vector)`, i.e. an adjunction between "cofree on one generator" and
"evaluate the coaction at one coordinate". This adjunction (and why cofree
comodules are the relative injectives needed for a resolution in this
category) is not derived in the code — it is simply implemented as the
`adjoint` method. The precise categorical statement should come from
`MinimalResolution.pdf`.

### `adjoint`: the embedding-detection map

`Hopf_Algebroid::adjoint(X, n, adjoint_map, shift)` (`hopf_algebroid/8.h:1-20`)
builds, for a single base generator `n` of `X`, the cofree comodule
co-generated by that one generator, and (if `adjoint_map` is non-null) the
matrix of the adjoint map `X → cofree(n)`, whose `i`-th row is
`algebroid2vector(component(n, X.coaction(i)), shift)` — i.e. project the
coaction of `X`'s `i`-th generator onto the `n`-th algebroid factor, then
re-express that algebroid element as a base-ring vector. The overload at
`hopf_algebroid/8.h:23-35` does the same for several generators at once
(`gens`), packing results side-by-side via `direct_sum` — this is the map
`irow` used inside `embed2cofree`.

### `embed2cofree`: one resolution half-step (the embedding)

`Hopf_Algebroid::embed2cofree` (`hopf_algebroid/9.h:1-83`) builds a cofree
comodule `F` and an injective map `inj : X → F`, one basis element of `X` at
a time, processed **in increasing underlying degree**
(`std::stable_sort` by `underlyingDeg`, line 23) — a comment-free but clear
sign that the algorithm is degree-by-degree, consistent with computing a
resolution "up to `maxDeg`" one internal degree at a time so that lower
degrees are fully resolved (and hence usable in reduction) before higher ones
are attempted:

1. Compute `irow = adjoint(X, gens, pos_of_gens, i, 0)` — the image of basis
   element `i` in the cofree comodule spanned by generators chosen *so far*.
2. `table->simplify_to_led(...)` reduces `irow` against the current
   `curtis_table`, tracking the homotopy `sc` (initialized to
   `singleton(i)`).
3. If reduction hits a non-boundary leading term (`pos != Boundary`) — i.e.
   `irow`, after subtracting known relations, still has an invertible leading
   coefficient — then element `i` is *already* hit by the existing cofree
   summands: record `table->insert(pos, i, irow, sc)` and `inj->insert(i,
   irow)`. No new generator is introduced.
4. If reduction instead collapses `irow` to the zero vector (`Boundary`),
   element `i` is **not** in the image of the current partial map — a new
   cofree summand (`adjoint(X, i, NULL, old_rank)`) must be added, `gens`
   grows by one, and the table is updated with this new generator as its own
   cycle/tag.

This is precisely "add generators only when the current map fails to be
surjective onto the needed target in this degree" — the minimality condition.

**TODO(math):** why "reduces to a non-boundary" (step 3) is the correct
criterion for "already covered", as opposed to some other test, and why
processing by increasing degree suffices for correctness/termination (rather
than needing degree-by-degree completion certificates), is not spelled out in
comments — likely covered by the minimal-resolution construction in
`MinimalResolution.pdf`.

### `resolvor`: one full resolution step (embed + quotient)

`Hopf_Algebroid::resolvor` (`hopf_algebroid/5.h:1-41`) is one full step:

1. `F = embed2cofree(X, inj, table, gens, ...)` — build the cofree comodule
   `F` and injection `inj : X ↪ F` as above.
2. `table->cycle_matrix(*indj, X.rank())` — extract from the table the
   injective matrix `indj` (rows = table cycles, i.e. one row per element of
   `X` mapped through the embedding) together with `inj_ind`, the leading-term
   (pivot) column index for each row.
3. Gaussian-eliminate `indj` using those pivot positions
   (`indj->gaussian(gs)` where `gs` pairs row `i` with pivot column
   `inj_ind[i]`) — this puts the embedding matrix into reduced echelon form.
4. `matrix<ring>::quot_index(inj_ind, F.rank())` computes which columns of
   `F` are *not* hit by any pivot — these indices span the **quotient**
   comodule `F / X`.
5. `indj->make_quotient(...)` builds the explicit quotient map
   `F → F/X` (transposing `(I,A)` into `(-A,I)`-style, per the comment at
   `matrices/6.h:25`).
6. `quotient_p(&F, quot, quot_inds.first, X)` recomputes `X`'s coaction *as
   the quotient comodule* `F/X` (using the quotient map applied on the right,
   with `etaR` twisted in via `right_scalor_mult`) — this **overwrites** the
   input `X` in place with the next comodule in the resolution.

So one `resolvor` call realizes the short exact sequence
`0 → X → F → X' → 0` (with `X'` the new value written into the `X` argument),
which is exactly one step `Ω⁻¹` of a cobar/relative-injective resolution: `F`
is the `s`-th term of the resolution, `inj` the `s`-th differential (into the
cofree term), and `X'` becomes the input to the next call.

**TODO(math):** the precise sense in which `F` is "the injective hull" (i.e.
why this particular cofree comodule, built exactly as in `embed2cofree`, is
guaranteed to be relative-injective and why the resulting resolution is
*minimal* in the technical sense used for the Adams-Novikov / algebraic
Novikov spectral sequence) is not present in code comments and should be
looked up in `MinimalResolution.pdf`.

### `pre_resolution_tab`: iterating steps into a resolution

`Hopf_Algebroid::pre_resolution_tab` (`hopf_algebroid/4.h:49-96`) just calls
`resolvor` `resolution_length+1` times in a loop, each time re-using (and
overwriting) `resolved` as the running comodule, saving `inj`/`qut` matrices
to `filename_maps` and the cofree generator data (`F.save`) to
`filename_generators`, and consuming one `curtis_table<ring>*` per step from
a caller-supplied vector `result` (so callers can retain/reuse the tables,
e.g. to later invert maps via `inverse.h`, or to reuse a table as a "model"
for `resolvor_modeled` on a related Hopf algebroid). `resolution()`
(`hopf_algebroid/6.h:1-69`) is the separate post-processing pass that splices
consecutive steps' `inj`/`qut` matrices into a single composed differential
per stage of the resulting chain complex.

`resolvor_modeled`/`pre_resolution_modeled`/`embed2cofree_modeled`
(`hopf_algebroid/3.h`, `.../9.h:85-135`) are variants that skip the
degree-by-degree cycle search and instead reuse a previously-computed
`curtis_table<table_type>` (from a *different*, "model" Hopf algebroid, via a
`transformer` function converting the model's ring elements into the current
ring) — i.e. transporting a known resolution shape from one setting to
another rather than recomputing it from scratch. **TODO(math):** the
conditions under which a resolution computed for one Hopf algebroid can be
validly transported ("modeled") onto another are not stated in the code.

## 6. Key method reference

### `matrix<R>` (`matrices/2.h` unless noted)

| Method | Description | Location |
|---|---|---|
| `find(matrix_index n) const` | Fetch row `n` as a sparse vector | `matrices/2.h:15` |
| `insert(i, x)` | Store row `i` | `matrices/2.h:18` |
| `gaussian(row_cols)` | Gaussian elimination using given (row, pivot-col) pairs | `matrices/2.h:124` (pure virtual, impl in `matrices_mem/4.h:33`) |
| `del_and_gaussian(row_cols, to_del)` | Delete columns then Gaussian-eliminate | `matrices/2.h:127` |
| `compose(rows, res_rank, result)` | Compose this matrix with a source given as a row-generating function | `matrices/5.h:12` |
| `maps_to(v)` | Right-multiply this matrix by row vector `v` (apply the map) | `matrices/5.h:3` |
| `quot_index(inj_index, rank)` (static) | Compute quotient-space indices/leading-term map from an injection's pivots | `matrices/6.h:3` |
| `make_quotient(inverse_ind, quot_ind, inj_index, result)` | Build the explicit quotient map from an echelon-form embedding | `matrices/6.h:27` |
| `direct_sum(Y, tar_rank)` | Block-sum with another matrix | `matrices/2.h:141` |
| `save`/`load` | Binary serialize/deserialize | `matrices/4.h:22`, `:41` |

### `ModuleOp<index,R>` (`modules/3.h` unless noted)

| Method | Description | Location |
|---|---|---|
| `add(x,y)` | Merge-add two sorted sparse vectors | `modules/4.h:3` |
| `scalor_mult(r,x)` | Scale a vector by a ring element | `modules/3.h:45` |
| `singleton(e,r)` | Build a one-term vector | `modules/3.h:73` |
| `component(n,x)` | Binary-search the coefficient at index `n` | `modules/6.h:29` |
| `re_index(rule,x)` / `filtered_reindex(rule,x,invalid)` | Reindex (optionally dropping invalid terms) | `modules/5.h:4`, `:33` |
| `filter(rule,x)` | Keep only terms matching a predicate | `modules/5.h:48` |
| `termwise_operation(rule,x)` | Map a function over every (index,coefficient) or coefficient | `modules/5.h:73`, `:88` |
| `next_non_trvial_entry(n,x)` / `first_invertible(x)` / `last_entry(x)` | Locate structurally-significant terms (used by `curtis_table` reduction) | `modules/6.h:34`,`:61`,`:81` |
| `sum(summands,start,end)` | Parallel (OpenMP) reduction-tree sum of a vector array | `modules/9.h:28` |

### `Hopf_Algebroid<ring,algebroid>` (`hopf_algebroid/4.h` unless noted)

| Method | Description | Location |
|---|---|---|
| `etaL`/`etaR` | Left/right unit maps `ring → algebroid` (pure virtual) | `hopf_algebroid/4.h:15-16` |
| `delta(matrix_index)` | Comultiplication on `i`-th base element (pure virtual) | `hopf_algebroid/4.h:24` |
| `adjoint(X, n, adjoint_map, shift)` | Build cofree-comodule-on-one-generator + adjoint map from `X` | `hopf_algebroid/8.h:4` |
| `adjoint(X, gens, pos, i, shift)` | Compute row `i` of the adjoint map for several generators at once | `hopf_algebroid/8.h:25` |
| `embed2cofree(X, inj, table, gens, ...)` | Embed `X` into a cofree comodule, minimally (degree order) | `hopf_algebroid/9.h:4` |
| `resolvor(X, inj, indj, quot, table, ...)` | One full resolution step: embed, Gaussian-eliminate, quotient | `hopf_algebroid/5.h:4` |
| `pre_resolution_tab(resolved, ..., resolution_length, result, gens, ...)` | Iterate `resolvor` for `resolution_length+1` steps, persisting to disk | `hopf_algebroid/4.h:51` |
| `quotient`/`quotient_p` | Compute the coaction of the quotient comodule (serial / OpenMP-parallel) | `hopf_algebroid/12.h:4`, `:30` |
| `resolution(...)` | Splice per-step short exact sequences into a chain complex | `hopf_algebroid/6.h:4` |
| `make_multiplication_table(x, deg_x, result)` | Tabulate right-multiplication-by-`x` on all basis elements below a degree | `hopf_algebroid/13.h:3` |
| `set_to_trivial(X, deg)` | Initialize `X` as the rank-1 trivial comodule (resolution starting point) | `hopf_algebroid/13.h:4` |
| `init_cofree_data(cofree_degree)` | Wire up `cofree_comodule`'s static coaction/degree functions to this algebroid's `delta` | `hopf_algebroid/13.h:17` |

### `curtis_table<ring>` (`matrices/9.h`)

| Method | Description | Location |
|---|---|---|
| `is_member(matrix_index)` / `search(matrix_index)` / `search_ref` | Table lookup by cycle index (pure virtual, backend-specific) | `matrices/9.h:19-23` |
| `symplify_to_led(x, homotopy)` | Reduce `x` against the table, accumulating a homotopy | `matrices/9.h:108` |
| `simplify_to_led(target_rank, source_rank, x, homotopy)` | Dense-vector variant of the above (converts to `std::vector<ring>` for speed) | `matrices/9.h:124` |
| `insert(cycle, tag, full_cycle, full_tag)` | Add a new table entry (pure virtual) | `matrices/9.h:156` |
| `cycle_matrix(cyc_mat, tag_mat, source_rank)` | Extract cycle/tag matrices + leading-term array from the table | `matrices/9.h:159` |
| `cycle_matrix(result, source_rank, transformer, Map)` | Transform-and-reinsert variant, used by the "modeled" resolution path | `matrices/9.h:191` |
| `save`/`load` | Binary serialize/deserialize the whole table | `matrices/9.h:215`, `:236` |

### `comodule_generic<algebroid,degree_type>` (`hopf_algebroid/2.h`)

| Method | Description | Location |
|---|---|---|
| `rank()` | Number of base generators | `hopf_algebroid/2.h:27` |
| `coaction(i)` | Coaction on generator `i`, read from `coaction_matrix` | `hopf_algebroid/2.h:31` |
| `degree(i)` | Degree of generator `i`, read from `base_module` | `hopf_algebroid/2.h:35` |

## Notably confusing / redundant (flag for top-level architecture doc)

- **`others.h` (397 lines) is dead code**: it defines a second, incompatible
  `curtis_table<R>`/`matrix_array<R>`/`quasi_table<R>` hierarchy under
  namespaces (`Modules::`, `Algebra::`, `MapTable::`, `Paralell::`) that
  appear nowhere else in the repo, and it `#include`s a nonexistent
  `filedmatrices.h`. It is not included from any other file. The real,
  live `curtis_table<ring>` lives in `matrices/9.h`, with backends
  `curtis_table_mem` (`matrices_mem/3.h`) and `curtisTable_stream`
  (`curtis.h`).
- **The numbered-fragment convention is easy to misread**: opening a
  mid-range fragment (e.g. `matrices/6.h`, `hopf_algebroid/10.h`) shows bare
  out-of-line template method bodies with no visible class declaration —
  the declaration lives in a lower-numbered fragment in the same directory,
  stitched together only via the parent manifest header's `#include` order.
  `matrices/9.h` is a trap in this scheme: unlike its siblings, it is not a
  fragment of `matrix<R>` at all but an entirely separate class
  (`curtis_table<ring>`) that only happens to live in the same directory.
- **`lift.h` looks unfinished/possibly non-compiling**: e.g. `adjoint()` at
  `lift.h:19-34` calls `HA_oper->algebroidModuleOper->component(i)` with a
  single argument on a `component(index, vec)` that needs two, and accesses
  `F.pos` where `cofree_comodule` has no member `pos` (only
  `position_of_gens`). Treat as an unverified/legacy utility, not a reference
  implementation.
