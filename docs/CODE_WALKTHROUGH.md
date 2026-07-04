# Code Walkthrough: From Math Objects to C++ Classes

This document is for readers who know *what* a comodule, a resolution, or a
Hopf algebroid is mathematically (see [`ARCHITECTURE.md`](ARCHITECTURE.md)
§1–2), but aren't fluent in reading C++, and want to know concretely: **which
class is the comodule, which member variable holds the coaction, and which
line of code actually runs the resolution?**

Everything here is a slower, more explicit version of material that's also
in [`FRAMEWORK.md`](FRAMEWORK.md), [`pipelines/BP.md`](pipelines/BP.md), and
[`pipelines/STEENROD.md`](pipelines/STEENROD.md) — those documents assume
you can read C++ fluently; this one doesn't. Every code snippet below is
quoted directly from the repository (with a file:line citation) and then
explained piece by piece.

## 0. A 60-second primer on the C++ used here

You'll see five constructs repeatedly. Here's what each one means, so the
rest of this document doesn't have to stop and explain them every time.

- **Templates — "fill-in-the-blank" classes.** A line like
  `template<typename algebroid, typename degree_type> class comodule_generic
  { ... }` declares a class that isn't finished yet — it has two blanks,
  named `algebroid` and `degree_type`. Nothing gets built until some other
  code "fills in the blanks" by writing e.g. `comodule_generic<BPBP, int>`.
  That fills `algebroid = BPBP` and `degree_type = int`, and the result is a
  genuine, concrete, usable class — as concrete as if you'd written a
  separate class by hand for exactly that combination. The whole point of
  this codebase's design is that the *math* (a comodule, a resolution
  algorithm) is written once as a fill-in-the-blank template, and then
  reused by filling in different rings/algebroids for the Steenrod, BP, and
  motivic pipelines.
- **Pointers (`T*`) — "the address of a `T`", not a `T` itself.** A member
  like `matrix<algebroid> *coaction_matrix` doesn't directly contain a
  matrix; it holds a reference to one that lives somewhere else. This is
  what lets the same class work whether that matrix is stored in memory or
  streamed from a file on disk — the pointer doesn't care, only the code on
  the other end of it does.
- **Inheritance and `virtual` — "is a kind of," with some blanks the child
  must fill in.** `class comodule_generic : virtual public CoModule<...>`
  means "a `comodule_generic` is a kind of `CoModule`." When a base class
  writes `virtual int rank() const = 0;` (note the `= 0`), it's declaring a
  *contract*: "every concrete comodule type must provide its own `rank()`
  function," without saying yet how any particular one computes it.
- **`typedef` — a nickname.** `typedef comodule_generic<P,int>
  SteenrodCoMod_generic;` just gives a short name to a long, filled-in
  template type, purely for readability. `SteenrodCoMod_generic` and
  `comodule_generic<P,int>` are exactly the same type.
- **Lambdas (`[this](int i){ return ...; }`) — a small, throwaway function
  defined inline**, right where it's used, instead of being given its own
  name elsewhere. `[this]` means "this little function is allowed to use
  the surrounding object's own data."

## 1. How a comodule is represented

Three classes matter here, all declared in `hopf_algebroid/2.h`. They form a
small hierarchy: one is a *contract* (interface), and two are different,
interchangeable ways of *fulfilling* that contract.

### 1.1 The contract: `CoModule<algebroid, degree_type>`

```cpp
template<typename algebroid, typename degree_type>
class CoModule{
public:
	virtual int rank() const=0;
	virtual vectors<matrix_index, algebroid> coaction(int i) const=0;
	virtual degree_type degree(int i) const=0;
};
```
(`hopf_algebroid/2.h:1-13`)

Read this as: "a comodule (over some algebroid, graded by some
`degree_type`) is anything that can answer three questions:
`rank()` — how many generators does it have (as a module over the base
ring)?
`coaction(i)` — what is the coaction on the `i`-th generator, expressed as a
sparse vector of algebroid elements?
`degree(i)` — what internal degree is the `i`-th generator in?
The `=0` on each means `CoModule` itself doesn't know how to compute any of
these — it's a promise that some other class will."

### 1.2 The concrete, matrix-backed comodule: `comodule_generic<algebroid, degree_type>`

This is the class you guessed at: "some class with a coaction matrix." Here
it is in full:

```cpp
template<typename algebroid, typename degree_type>
class comodule_generic : virtual public CoModule<algebroid, degree_type>{
public:
	matrix<algebroid> *coaction_matrix;
	modules<degree_type> base_module;

	comodule_generic(matrix<algebroid> *coactor) { coaction_matrix = coactor; }

	int rank() const{
		return base_module.rank; }

	vectors<matrix_index, algebroid> coaction(int i) const{
		return coaction_matrix->find(i); }

	degree_type degree(int i) const{
		return base_module.degree[i]; }
};
```
(`hopf_algebroid/2.h:15-37`)

Reading the template line first: `comodule_generic<algebroid, degree_type>`
has two blanks. `algebroid` is filled in with whatever ring the coaction
lands in — for the BP pipeline that's `BPBP` (i.e. `BP_*BP`), for the
classical Steenrod pipeline it's `P` (i.e. `BP_*BP/I`). `degree_type` is
filled in with whatever type is used to record a generator's grading — in
every pipeline documented here that's just a plain `int` (a single integer
degree), except the motivic pipeline, which uses a small `MotDegree` class
(bidegree `(t, weight)`) instead.

Now the members, concretely:

- **`base_module.rank`** — a plain integer: how many generators this
  comodule has. `rank()` just returns it.
- **`base_module.degree`** — an array with one entry per generator, giving
  each generator's degree. `degree(i)` just looks up entry `i`.
- **`coaction_matrix`** — a *pointer* to a `matrix<algebroid>` (a matrix
  whose entries are elements of the algebroid ring). Row `i` of this matrix
  *is* the coaction on generator `i`, written out as a sparse linear
  combination — `coaction(i)` is literally `coaction_matrix->find(i)`, "go
  fetch row `i`." So: **the coaction map, as data, is nothing more than a
  matrix, one row per generator.** Because it's a pointer rather than an
  owned object, the matrix backing a `comodule_generic` can be an in-memory
  matrix (fast, but limited by RAM) or a disk-file-backed matrix (slower,
  but scales to computations too large to fit in memory) — the
  `comodule_generic` code above doesn't change either way; only which
  concrete matrix class gets pointed to changes. This is exactly the
  difference between the classical Steenrod pipeline (in-memory,
  `matrix_mem<P>`) and the BP pipeline (disk-backed, `matrix_file<BPBP>`) —
  see the table in §2 below.

This is the class that represents "the comodule currently being resolved,"
at every stage of every pipeline's resolution loop.

### 1.3 The concrete, formula-backed comodule: `cofree_comodule<algebroid, degree_type>`

The *other* place comodules show up is as the resolution's own terms — the
`F_0, F_1, F_2, ...` in a resolution `M → F_0 → F_1 → ...`. These are always
**cofree** (a direct sum of copies of the algebroid itself, with degree
shifts — see `ARCHITECTURE.md` §2 point 1 for the mathematical definition).
Storing an explicit coaction matrix for one of these would be wasteful,
since a cofree comodule's coaction is *always* just the algebroid's own
comultiplication, shifted around — so `cofree_comodule` computes it on
demand instead of storing it:

```cpp
cofree_comodule<algebroid,degree_type>::cofree_coaction;      // a shared function
cofree_comodule<algebroid,degree_type>::modOpers;             // shared module operations
cofree_comodule<algebroid,degree_type>::cofree_degree;        // a shared degree function
```
(static members, `hopf_algebroid/12.h:29-33`)

The word **`static`** here means: these three aren't three different values
per `cofree_comodule` object — they're **one shared value for every
`cofree_comodule` in the whole program** (of a given `<algebroid,
degree_type>` combination), like a single shared blackboard rather than each
object having its own notebook. They get written onto that blackboard once,
via:

```cpp
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::init_cofree_data(std::function<degree_type(matrix_index)> cofree_degree){
	cofree_comodule<algebroid,degree_type>::cofree_coaction = [this] (matrix_index n){
		return delta(n); };
	cofree_comodule<algebroid,degree_type>::modOpers = algebroidModuleOper;
	cofree_comodule<algebroid,degree_type>::cofree_degree = cofree_degree;
}
```
(`hopf_algebroid/12.h:14-25`)

In plain English: "the shared coaction formula for any cofree comodule is
*literally the algebroid's own comultiplication*, `delta`." So a
`cofree_comodule`'s `coaction(i)` (`hopf_algebroid/11.h:62-68`) doesn't read
a matrix at all — it works out which cofree summand index `i` falls into,
then calls this shared `delta`-based formula and shifts the result into the
right place. This is *why* the framework can build resolution terms that
are, in principle, enormous (a big direct sum of copies of the whole
algebroid) without ever materializing that whole thing as a matrix.

## 2. The Hopf algebroid itself: `Hopf_Algebroid<ring, algebroid>`

Same fill-in-the-blank idea, one level up. `Hopf_Algebroid<ring, algebroid>`
(`hopf_algebroid/4.h`) has two blanks: `ring` (the base ring `A`) and
`algebroid` (the bigger ring `Γ`), and it's what actually declares
`etaL`, `etaR`, `delta`, and the whole resolution algorithm (`embed2cofree`,
`resolvor`, `pre_resolution_tab`, ... — see `FRAMEWORK.md` §5). Each concrete
pipeline fills in the blanks differently:

| Pipeline | Concrete class | `ring` (=`A`) | `algebroid` (=`Γ`) | Comodule type used (§1.2) | Coaction matrix backend |
|---|---|---|---|---|---|
| BP (`mr_BP`) | `class BP_Op : Hopf_Algebroid<BP, BPBP>` (`BP.h:33`) | `BP` = `polynomial<Z3>` (i.e. `BP_*`) | `BPBP` = `polynomial<BP>` (i.e. `BP_*BP`) | `BPCoMod_generic` = `comodule_generic<BPBP,int>` (`BP.h:122`) | `matrix_file<BPBP>` — disk-backed (`BP_init.h:11-12`) |
| Classical Steenrod (`mr_st`) | `class Steenrod_Op : Hopf_Algebroid<Fp, P>` (`steenrod.h:38`) | `Fp` (i.e. `F_3`) | `P` = `polynomial<Fp>` (i.e. `BP_*BP/I`) | `SteenrodCoMod_generic` = `comodule_generic<P,int>` (`steenrod.h:102`) | `matrix_mem<P>` — in memory (`steenrod_init.h:7-8`) |

Concretely, this is why `docs/pipelines/BP.md` and `docs/pipelines/
STEENROD.md` describe two structurally near-identical driver classes
(`BPInit`/`SteenrodInit`, `BPComodInit`/`ComodInit`) — they're the *same*
generic pattern (one `Hopf_Algebroid` instantiation + one matrix-backed
comodule + a resolution loop), filled in with different rings and a
different storage backend for the comodule's matrix.

## 3. The trivial comodule: `set_to_trivial`

Every pipeline needs a starting point: resolve the *trivial* comodule (rank
1, one generator, everything else derived from that). This is one shared
generic function, not something re-implemented per pipeline:

```cpp
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::set_to_trivial(comodule_generic<algebroid,degree_type> &X, degree_type deg) {
	X.base_module.rank = 1;
	X.base_module.degree.resize(1);
	X.base_module.degree[0] = deg;

	std::function<vectors<matrix_index,algebroid>(int)> rw = [this] (int i){
		return algebroidModuleOper->singleton(0,algebroidRingOper->unit(1)); };
	X.coaction_matrix->construct(1,rw);
}
```
(`hopf_algebroid/12.h:1-12`)

Line by line, in plain English:

1. `X.base_module.rank = 1;` — this comodule has exactly one generator.
2. `X.base_module.degree[0] = deg;` — that one generator sits in whatever
   degree the caller asked for (every pipeline documented here passes
   `deg = 0`).
3. The `rw` lambda is the "recipe for computing row `i` of the coaction
   matrix" (only `i=0` will ever be asked for, since there's one
   generator). Its recipe: `algebroidRingOper->unit(1)` is the algebroid's
   own multiplicative identity element (its "1"); `singleton(0, ...)` builds
   a one-term sparse vector saying "this coefficient, at position 0." So
   the coaction on the single generator is "1 times itself" — the trivial
   comodule's coaction is as simple as it can possibly be, and doesn't
   depend on which concrete algebroid `Γ` you've plugged in.
4. `X.coaction_matrix->construct(1, rw)` — build the (one-row) coaction
   matrix using that recipe, and store it via the pointer described in §1.2.

Every pipeline's initialization code calls this *exact* function, right
before starting its resolution loop:

- `steenrod_init.cpp:32` — `steenrod_oper.set_to_trivial(comod, 0);`
- `BP_init.cpp:23` — `BP_oper.set_to_trivial(comod, 0);`
- `mot_main.cpp:33` and `mot_combine.cpp:29` — `MOP.set_to_trivial(comod, zo);`
- `ex_steenrod_init.cpp:34` — `steenrod_oper.set_to_trivial(comod, 0);`

There is no other "trivial comodule" logic anywhere in the codebase — it's
this one routine, reused everywhere.

## 4. Where the BP/I and BP resolutions actually run

This traces the exact call path from each executable's `main()` down to the
line where the resolution loop executes.

### 4.1 The BP/I resolution (`mr_st`)

```cpp
// stmain.cpp:11-15
SteenrodInit st(3, maxdeg, length, filename + "steenrod_coaction.data");
st.resolve(filename);
```

`SteenrodInit`'s constructor (`steenrod_init.cpp:5-39`) does the setup:
builds the coproduct table, then at line 32 calls `set_to_trivial` (§3
above) to initialize `comod` — a `ComodInit` (`steenrod_init.h:7-12`), which
is a `SteenrodCoMod_generic` (§1.2/§2) backed by its own in-memory
`matrix_mem<P> coaction_matrix`.

`st.resolve(filename)` (`steenrod_init.cpp:46-48`) is one line:

```cpp
steenrod_oper.pre_resolution_tab(comod, director + "maps", director + "gens",
                                  resolution_length, resolutionTables, gens,
                                  &inj, &indj, &qut, &new_map, director + "BPtables");
```

`pre_resolution_tab` (declared `hopf_algebroid/4.h:51`, defined
`hopf_algebroid/5.h:51-96`) is the actual resolution loop: it calls
`resolvor` (one embed-and-quotient step, `FRAMEWORK.md` §5) repeatedly,
overwriting `comod` each time with the next comodule in the resolution, and
writing the maps/generators/curtis-tables to disk as it goes. **This loop is
where `Ext_P(F_3, F_3)` actually gets computed**, one resolution step at a
time.

### 4.2 The BP resolution (`mr_BP`)

```cpp
// BPmain.cpp:18-27
BPInit BPoper(max_degree, resolution_length, filename0+"etaL", filename0+"R2L", filename0+"delta", filename);
BPoper.load_gens(filename0 + "gens_data");   // the mr_st output — see ARCHITECTURE.md §5
BPoper.resolve();
```

`BPInit`'s constructor (`BP_init.cpp:9-40`) mirrors `SteenrodInit`'s: it
calls `BP_oper.set_to_trivial(comod, 0)` at line 23 to initialize `comod` —
a `BPComodInit` (`BP_init.h:11-16`), a `BPCoMod_generic` (§1.2/§2) backed by
a *disk-file* `matrix_file<BPBP> coaction_matrix` this time, because a full
BP resolution can be far too large to keep in memory.

`BPoper.resolve()` (`BP_init.cpp:43-54`) is:

```cpp
BP_oper.pre_resolution_modeled(comod, director+"maps", director+"gens",
                                resolution_length, director+"tables", &ctable, gens,
                                tfm, &inj, &qut, &indj, &new_map, director+"back");
```

This is **the one line where the actual BP resolution is computed** —
`pre_resolution_modeled` (declared `hopf_algebroid/4.h:83`, defined
`hopf_algebroid/6.h:41-...`) is the "lift a known model resolution" variant
of `pre_resolution_tab` (`ARCHITECTURE.md` §5). Two of its arguments are
worth spelling out:

- **`&ctable`** — a pointer to the *already-computed* curtis table from
  `mr_st`'s run (loaded from the `mr_st`-produced files); this is the
  "model" being lifted, playing the role `resolutionTables` played in
  §4.1's `pre_resolution_tab` call, except reused rather than recomputed.
- **`tfm`** — a small lambda (`BP_init.cpp:47-48`) defined as
  `[this](const vectors<matrix_index,Fp>& v){ return BP_oper.lift(v); }`:
  "given a vector over `F_p` (the model's ring), lift it to a vector over
  `BP` (the ring being actually resolved)." This is the concrete "transport
  a resolution from one Hopf algebroid to another" step —
  `pre_resolution_modeled`'s generic `transformer` parameter
  (`hopf_algebroid/4.h:83`), filled in here with `BP_Op::lift`.

## 5. Cheat sheet

| Math concept | C++ syntax | Plain-English meaning |
|---|---|---|
| A comodule (abstractly) | `CoModule<algebroid, degree_type>` | "Anything with a `rank()`, `coaction(i)`, `degree(i)`" — a contract, no data |
| A comodule with an explicit coaction | `comodule_generic<algebroid, degree_type>` | `rank`/`degree` are a plain int + array; `coaction(i)` is "row `i` of a matrix" |
| A cofree (relative-injective) comodule | `cofree_comodule<algebroid, degree_type>` | Coaction computed on demand from the algebroid's own comultiplication, never stored as a matrix |
| The Hopf algebroid `(A, Γ)` itself, plus the resolution algorithm | `Hopf_Algebroid<ring, algebroid>` | `ring`=`A`, `algebroid`=`Γ`; declares `etaL`/`etaR`/`delta` and implements `embed2cofree`/`resolvor`/`pre_resolution_tab` generically |
| "Fill in the blanks" for BP | `BP_Op : Hopf_Algebroid<BP, BPBP>` | `A = BP_* `, `Γ = BP_*BP` |
| "Fill in the blanks" for BP/I | `Steenrod_Op : Hopf_Algebroid<Fp, P>` | `A = F_3`, `Γ = BP_*BP/I` |
| The trivial (rank-1) comodule | `Hopf_Algebroid::set_to_trivial(X, deg)` | One generator, coaction = "1 times itself," reused by every pipeline |
| Where a resolution loop actually runs | `pre_resolution_tab` (from scratch) / `pre_resolution_modeled` (lifted from a model) | One function call in each pipeline's `Init::resolve()` |
