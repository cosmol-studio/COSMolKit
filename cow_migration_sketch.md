# COSMolKit Copy-on-Write Migration Sketch

## Goal

Introduce a clear pandas 3.0 / Polars-style public molecule model without paying full deep-copy cost on every value-returning API call.

This is primarily meant to solve the API ambiguity that exists in RDKit's mixed inplace/non-inplace style, while preserving RDKit-compatible behavior in the internal algorithm paths where compatibility is the goal.

COSMolKit's public model should be deterministic and predictable: normal molecule transformations return new molecule objects and never modify the input. This mirrors pandas 3.0 Copy-on-Write and Polars-style explicit dataflow: ordinary transformations produce new values, while mutation belongs behind clearly named editing APIs.

The target model is:

- public `Molecule` is cheap to clone
- user-facing transforms return a new molecule by default
- input molecules are unchanged after normal transformations
- only mutated storage is copied
- complex RDKit-aligned algorithms still run on explicit mutable working state when needed

This is an API/data-layout change, not a chemistry-semantics change.

## Non-Goal

Do not force every internal algorithm to operate directly on copy-on-write storage.

For parity-critical code paths such as:

- kekulization
- SMILES writer canonicalization
- stereochemistry assignment
- distance-geometry bounds
- molblock/SDF write pipelines

it is still acceptable, and often preferable, to materialize an explicit mutable working copy and translate RDKit source logic onto that state directly.

## Why COW Fits COSMolKit

The intended COSMolKit Python/Rust API is deliberately aligned with modern dataframe API design:

- pandas 3.0-style Copy-on-Write behavior at the public API boundary
- Polars-style explicit transformation chains that return new values
- explicit editing contexts
- no `inplace=True`-style API ambiguity

That API style maps naturally to copy-on-write:

- `mol2 = mol.with_hydrogens()` is logically "return a new molecule"
- internally we can share data first and only detach the modified storage on write

This avoids the two bad extremes:

1. RDKit-style implicit mutation ambiguity
2. naive full deep-copy on every returned molecule

## Recommended Storage Split

Do not start with one giant `Arc<MolData>`.

That is simple, but too coarse for chemistry workloads where topology, coordinates, and properties often change independently.

Recommended split:

```rust
pub struct Molecule {
    topology: Arc<TopologyData>,
    conformers: Arc<ConformerStore>,
    props: Arc<PropertyStore>,
}
```

### `TopologyData`

Contains graph and graph-derived persistent state:

```rust
pub struct TopologyData {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    adjacency: AdjacencyList,
}
```

Future candidates:

- cached ring info
- fragment membership
- stable stereo-supporting cached topology data

Topology should own only data that is conceptually tied to graph identity.

### `ConformerStore`

Contains coordinates only:

```rust
pub struct ConformerStore {
    conformers: Vec<Conformer>,
}
```

This lets operations such as:

- 2D coordinate generation
- translation/rotation
- conformer replacement

detach only coordinates, without cloning atoms/bonds.

### `PropertyStore`

Contains user or IO metadata:

```rust
pub struct PropertyStore {
    sdf_data_fields: Vec<(String, String)>,
    name: Option<String>,
    // later: typed property map if needed
}
```

This isolates cheap metadata updates from graph/coordinate duplication.

## Public API Shape

The public `Molecule` API should remain value-like:

```rust
let mol_h = mol.with_hydrogens()?;
let mol_2d = mol_h.with_2d_coords()?;
let mol_named = mol_2d.with_name("ligand");
```

Implementation idea:

```rust
impl Molecule {
    pub fn topology_mut(&mut self) -> &mut TopologyData {
        Arc::make_mut(&mut self.topology)
    }

    pub fn conformers_mut(&mut self) -> &mut ConformerStore {
        Arc::make_mut(&mut self.conformers)
    }

    pub fn props_mut(&mut self) -> &mut PropertyStore {
        Arc::make_mut(&mut self.props)
    }
}
```

This gives the desired semantics:

- `clone()` is cheap
- first mutation detaches the touched block only
- untouched blocks remain shared

## Internal Working-State Layer

Do not use copy-on-write as the only mutable model.

For RDKit-aligned algorithm ports, add explicit working-state types.

Examples:

```rust
pub struct EditableMolecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    adjacency: AdjacencyList,
    conformers: Vec<Conformer>,
    props: PropertyStore,
}
```

or narrower task-specific forms:

```rust
pub struct SmilesWriteState<'a> {
    base: &'a Molecule,
    atom_chiral_tags: Vec<ChiralTag>,
    bond_directions: Vec<BondDirection>,
    bond_stereos: Vec<BondStereo>,
    bond_orders: Vec<BondOrder>,
    atom_aromatic_flags: Vec<bool>,
    bond_aromatic_flags: Vec<bool>,
    mol_stack: Vec<MolStackElem>,
    atom_visit_orders: Vec<usize>,
    bond_visit_orders: Vec<usize>,
}
```

This separation is important:

- `Molecule` is the persistent shared value object
- working state is explicit, mutable, and algorithm-local

That matches RDKit source translation much better than trying to express every intermediate writer/sanitizer mutation through `Arc::make_mut`.

## When Each Block Should Detach

### Topology detaches

Operations that change graph identity or graph-derived state:

- add/remove hydrogens
- parse from SMILES/SDF
- sanitize / aromaticity / valence updates if persisted on atoms/bonds
- kekulize if persisted on atoms/bonds
- bond order changes
- stereo flags stored on atoms/bonds

### Conformers detach

Operations that change only coordinates:

- compute 2D coordinates
- future 3D embedding
- transform/translate/rotate conformers
- replace conformer positions

### Properties detach

Operations that change only metadata:

- set molecule name
- set SDF data fields
- later arbitrary user properties

## Immediate Design Constraint from Current Codebase

Current `cosmolkit-core` code frequently assumes direct ownership of:

- `mol.atoms`
- `mol.bonds`
- `mol.adjacency`

So a full immediate conversion of `Molecule` internals to `Arc<...>` would ripple through:

- `smiles.rs`
- `smiles_write.rs`
- `canon_smiles.rs`
- `kekulize.rs`
- `stereo.rs`
- `distgeom/*`
- `io/molblock.rs`
- `io/sdf.rs`
- Python bindings

This is too much churn to mix with active parity debugging.

Therefore migration should be staged.

## Recommended Migration Stages

### Stage 0: Keep parity work stable

Before touching storage, continue landing parity fixes against the current plain-owned `Molecule`.

Reason:

- parity failures should reflect chemistry/logic differences
- not data-layout churn

### Stage 1: Introduce internal storage types without changing public behavior

Add:

- `TopologyData`
- `ConformerStore`
- `PropertyStore`

But keep `Molecule` layout unchanged for one short transition period if needed.

Purpose:

- isolate responsibilities
- reduce future diff size

Possible transitional form:

```rust
pub struct Molecule {
    topology: TopologyData,
    conformers: ConformerStore,
    props: PropertyStore,
}
```

This already improves structure before adding `Arc`.

### Stage 2: Switch persistent blocks to `Arc`

After Stage 1 is stable:

```rust
pub struct Molecule {
    topology: Arc<TopologyData>,
    conformers: Arc<ConformerStore>,
    props: Arc<PropertyStore>,
}
```

Then add:

- `topology_mut()`
- `conformers_mut()`
- `props_mut()`

and keep algorithm-facing helper methods so existing code does not immediately touch `Arc` internals everywhere.

### Stage 3: Add explicit editable/working types

Add:

- `EditableMolecule`
- `SmilesWriteState`
- later possibly `SanitizeWorkingState`

Then move high-mutation internal code away from mutating `Molecule` directly when that improves parity readability.

This should start with the most mutation-heavy or RDKit-structured code:

1. SMILES writer
2. kekulization
3. stereo assignment

### Stage 4: Align Python API with immutable-first semantics

Once COW is in place, Python bindings can expose clearer semantics without paying full clone cost:

- `mol.with_hydrogens()`
- `mol.without_hydrogens()`
- `mol.with_2d_coords()`
- `mol.with_name(...)`

Optional explicit mutation surface:

- `mol.edit()`
- `editor.commit()`

This is the right time to expose more "modern API" behavior, because storage cost will already be controlled.

## Suggested First Refactor Targets

These are the safest first candidates after Stage 0 parity stabilization.

### 1. Conformer split

This is the most mechanically isolated split.

Reason:

- coordinates are conceptually separate already
- upcoming 2D/3D APIs benefit immediately
- fewer parity-critical graph algorithms depend on conformer storage layout

### 2. Property split

Also low risk.

Reason:

- minimal chemistry impact
- helps SDF/Python metadata API

### 3. Topology split

Do this after the previous two so the graph-bearing piece is the only hard part left.

## Interaction with Current SMILES Writer Design

The current writer path already exposed the need for an explicit mutable writer state.

Recent fixes translated RDKit `Canon.cpp` double-bond direction canonicalization and applied it on a temporary molecule, then copied final results into traversal overrides.

That is acceptable as an intermediate parity-preserving step.

Longer-term, it should evolve into a proper writer-local mutable state object instead of accumulating more `*_overrides` fields on `FragmentTraversal`.

Recommended future direction:

- keep `FragmentTraversal` for traversal structure only
- move temporary mutable atom/bond output state into `SmilesWriteState`

## Risks

### Risk 1: Mixing parity logic changes with storage changes

Avoid this.

If a parity regression appears during migration, it must be obvious whether the cause is:

- chemistry logic
- data-layout / ownership refactor

### Risk 2: Overusing COW in internal hot paths

Do not force every internal algorithm to bounce through `Arc::make_mut`.

For long mutation sequences, convert to a dedicated mutable working object first.

### Risk 3: Letting shared storage leak mutability semantics into API

The user should still see a crisp semantic distinction:

- default methods return a new molecule
- explicit editing APIs are clearly named

Copy-on-write is an implementation detail, not a semantic excuse for hidden mutation.

## Minimal Example

Desired behavior:

```rust
let mol = Molecule::from_smiles("CCO")?;
let mol2 = mol.with_2d_coords()?;
let mol3 = mol2.with_name("ethanol");
```

Possible detach pattern:

- `mol2` shares topology + props with `mol`, detaches conformers only
- `mol3` shares topology + conformers with `mol2`, detaches props only

This is exactly the kind of persistent value behavior we want.

## Recommendation

Adopt copy-on-write, but not as an immediate repo-wide rewrite.

The recommended end state is:

- persistent public `Molecule` with split COW blocks
- explicit mutable working states for parity-critical algorithms
- immutable-first Python/Rust API on top

That gives COSMolKit a cleaner semantic model than RDKit while keeping selected behavior compatible with RDKit where compatibility is required.

## Concrete Execution Plan

Use the staged design above, but split implementation into seven reviewable
steps so parity and data-layout regressions are easy to separate.

### Step 1: Freeze Current Parity Baseline

Status: effectively complete before storage migration.

Required checks:

- `cargo test --workspace`
- `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- `.venv/bin/pytest`
- targeted RDKit parity tests when a touched module has its own golden suite

No COW storage changes should land in this step.

### Step 2: Introduce Owned Storage Blocks

Goal: add the data blocks without `Arc` yet.

Implementation target:

```rust
pub struct Molecule {
    topology: TopologyData,
    conformers: ConformerStore,
    props: PropertyStore,
}
```

For the current codebase, this step should be transition-friendly:

- `TopologyData` owns `atoms`, `bonds`, and `adjacency`.
- `ConformerStore` owns `coords_2d`, `conformers_3d`, and
  `source_coordinate_dim`.
- `PropertyStore` starts small and can remain empty until molecule-level
  metadata APIs are added.
- Keep compatibility helpers so parity-heavy code does not need a full rewrite
  in the same commit.

Validation:

- No public behavior change.
- No RDKit parity golden updates.
- Benchmark before and after the step using the debug Python extension.

### Step 3: Replace Direct Storage Access With Helper Methods

Goal: prepare for real `Arc::make_mut` without touching chemistry behavior.

Add and use helper methods such as:

- `atoms()`
- `atoms_mut()`
- `bonds()`
- `bonds_mut()`
- `topology_mut()`
- `conformers_mut()`
- `props_mut()`

This step is mostly mechanical and should avoid algorithmic edits.

### Step 4: Switch Storage Blocks To `Arc`

Goal: make public `Molecule::clone()` cheap.

Implementation target:

```rust
pub struct Molecule {
    topology: Arc<TopologyData>,
    conformers: Arc<ConformerStore>,
    props: Arc<PropertyStore>,
}
```

Mutation helpers should detach only the block being modified:

- topology-changing operations detach topology
- coordinate-only operations detach conformers
- metadata-only operations detach props

### Step 5: Move Public Transform APIs Onto COW

Goal: preserve immutable-first semantics while avoiding full deep copies.

Primary APIs:

- `with_hydrogens()`
- `without_hydrogens()`
- `with_2d_coords()`
- `with_kekulized_bonds()`
- future `with_name(...)`

Each method should clone the public value cheaply and detach only the needed
storage block.

### Step 6: Add Explicit Working-State Types

Goal: keep RDKit source-level algorithm ports readable.

Do not force high-mutation parity code through persistent COW storage.

Priority:

1. SMILES writer state
2. kekulization state
3. stereo assignment state

This is where temporary `*_overrides` fields should be retired when a clearer
working-state object exists.

### Step 7: Align Python Docs, Examples, And Tests

Goal: make COW semantics visible and tested at the API boundary.

Add tests for:

- value transforms do not mutate input molecules
- coordinate-only operations preserve topology sharing after Step 4
- metadata-only operations preserve topology and conformer sharing after Step 4
- editing APIs remain explicit

Update Python docs and examples after behavior is implemented, not before.

## Benchmark Protocol

Before each storage-migration step, run:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Use the debug Python extension for now. Do not switch to release numbers until
the COW design is settled.

The benchmark currently tracks:

- SDF read + atom/bond count
- SDF read + canonical isomeric SMILES
- cached molecule to canonical isomeric SMILES
- cached add hydrogens only
- cached explicit-H molecule to canonical isomeric SMILES
- cached remove hydrogens only
- cached add/remove hydrogens only
- cached add hydrogens + canonical isomeric SMILES
- cached remove hydrogens + canonical isomeric SMILES
- cached add/remove hydrogens + canonical isomeric SMILES
- COW/value-semantics add hydrogens preserving the input molecule
- COW/value-semantics 2D coordinate generation preserving the input molecule

Record the debug relative speeds in this document after each step.

### Benchmark Log

#### Before Step 2

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.535 ms | 20.151 ms | 37.685x |
| `read_sdf_to_smiles` | yes | 0.788 ms | 78.817 ms | 100.080x |
| `cached_mol_to_smiles` | yes | 0.223 ms | 58.548 ms | 262.488x |
| `cached_add_remove_hs_to_smiles` | yes | 1.312 ms | 127.852 ms | 97.437x |

#### After Step 2

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.505 ms | 21.114 ms | 41.838x |
| `read_sdf_to_smiles` | yes | 0.794 ms | 81.295 ms | 102.406x |
| `cached_mol_to_smiles` | yes | 0.234 ms | 60.470 ms | 257.993x |
| `cached_add_remove_hs_to_smiles` | yes | 1.324 ms | 133.467 ms | 100.820x |

#### Before Step 3

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension from Step 2
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.514 ms | 20.507 ms | 39.863x |
| `read_sdf_to_smiles` | yes | 0.794 ms | 79.895 ms | 100.658x |
| `cached_mol_to_smiles` | yes | 0.287 ms | 59.655 ms | 207.681x |
| `cached_add_remove_hs_to_smiles` | yes | 1.126 ms | 132.308 ms | 117.487x |

#### After Step 3

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.699 ms | 20.619 ms | 29.504x |
| `read_sdf_to_smiles` | yes | 0.849 ms | 77.656 ms | 91.445x |
| `cached_mol_to_smiles` | yes | 0.226 ms | 58.089 ms | 256.720x |
| `cached_add_remove_hs_to_smiles` | yes | 1.377 ms | 127.107 ms | 92.322x |

#### Before Step 4

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension from Step 3
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.537 ms | 20.662 ms | 38.469x |
| `read_sdf_to_smiles` | yes | 0.814 ms | 78.354 ms | 96.272x |
| `cached_mol_to_smiles` | yes | 0.221 ms | 58.346 ms | 264.602x |
| `cached_add_remove_hs_to_smiles` | yes | 1.124 ms | 129.398 ms | 115.112x |

#### After Step 4

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.741 ms | 21.646 ms | 29.203x |
| `read_sdf_to_smiles` | yes | 0.885 ms | 78.788 ms | 88.980x |
| `cached_mol_to_smiles` | yes | 0.227 ms | 59.099 ms | 260.599x |
| `cached_add_remove_hs_to_smiles` | yes | 1.144 ms | 128.083 ms | 111.935x |

#### Before Step 5

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension from Step 4
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.639 ms | 20.164 ms | 31.573x |
| `read_sdf_to_smiles` | yes | 1.006 ms | 78.278 ms | 77.831x |
| `cached_mol_to_smiles` | yes | 0.249 ms | 58.219 ms | 233.495x |
| `cached_add_remove_hs_to_smiles` | yes | 1.121 ms | 127.138 ms | 113.428x |

#### After Step 5

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.513 ms | 20.245 ms | 39.459x |
| `read_sdf_to_smiles` | yes | 0.813 ms | 77.440 ms | 95.301x |
| `cached_mol_to_smiles` | yes | 0.224 ms | 57.254 ms | 255.349x |
| `cached_add_remove_hs_to_smiles` | yes | 1.138 ms | 127.580 ms | 112.079x |

#### Before Step 6

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension from Step 5
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.519 ms | 20.606 ms | 39.674x |
| `read_sdf_to_smiles` | yes | 0.834 ms | 76.907 ms | 92.175x |
| `cached_mol_to_smiles` | yes | 0.276 ms | 56.827 ms | 205.757x |
| `cached_add_remove_hs_to_smiles` | yes | 1.342 ms | 128.090 ms | 95.443x |

#### After Step 6A

Step 6A introduced a SMILES writer-local `SmilesWriteState` so atom SMILES
emission reuses one RDKit-like valence assignment per `mol_to_smiles()` call
instead of recomputing valence for each emitted atom field.

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.519 ms | 20.116 ms | 38.757x |
| `read_sdf_to_smiles` | yes | 0.786 ms | 38.129 ms | 48.520x |
| `cached_mol_to_smiles` | yes | 0.228 ms | 18.014 ms | 78.995x |
| `cached_add_remove_hs_to_smiles` | yes | 1.105 ms | 90.133 ms | 81.532x |

#### Before Step 6B

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension from Step 6A
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.508 ms | 20.318 ms | 40.001x |
| `read_sdf_to_smiles` | yes | 0.786 ms | 37.994 ms | 48.343x |
| `cached_mol_to_smiles` | yes | 0.236 ms | 18.110 ms | 76.876x |
| `cached_add_remove_hs_to_smiles` | yes | 1.163 ms | 87.859 ms | 75.529x |

#### After Step 6B

Step 6B moved canonical double-bond direction/stereo mutation off a cloned
`Molecule` and into a writer-local `CanonBondState`. This is primarily a
storage-boundary cleanup: the benchmark is expected to be near-neutral because
Step 6A already removed the dominant repeated valence assignment cost.

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.509 ms | 20.441 ms | 40.126x |
| `read_sdf_to_smiles` | yes | 0.905 ms | 38.775 ms | 42.826x |
| `cached_mol_to_smiles` | yes | 0.224 ms | 18.103 ms | 80.842x |
| `cached_add_remove_hs_to_smiles` | yes | 1.398 ms | 91.582 ms | 65.519x |

#### Before Step 6C

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension from Step 6B
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.554 ms | 20.047 ms | 36.207x |
| `read_sdf_to_smiles` | yes | 0.889 ms | 38.226 ms | 42.977x |
| `cached_mol_to_smiles` | yes | 0.226 ms | 17.977 ms | 79.523x |
| `cached_add_remove_hs_to_smiles` | yes | 1.248 ms | 91.668 ms | 73.450x |

#### After Step 6C

Step 6C made canonical ranking, fragment canonicalization, and atom emission
share the same `SmilesWriteState` valence assignment on the ordinary
non-kekule `mol_to_smiles()` path. The benchmark remained near-neutral, which
suggests the current writer hot path is dominated by later stereo/ranking work
rather than the remaining top-level valence assignment calls.

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.511 ms | 20.507 ms | 40.097x |
| `read_sdf_to_smiles` | yes | 0.793 ms | 38.034 ms | 47.962x |
| `cached_mol_to_smiles` | yes | 0.230 ms | 18.486 ms | 80.473x |
| `cached_add_remove_hs_to_smiles` | yes | 1.497 ms | 89.173 ms | 59.575x |

#### Before Step 6D

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension from Step 6C
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.535 ms | 20.229 ms | 37.826x |
| `read_sdf_to_smiles` | yes | 0.893 ms | 38.083 ms | 42.641x |
| `cached_mol_to_smiles` | yes | 0.220 ms | 17.862 ms | 81.030x |
| `cached_add_remove_hs_to_smiles` | yes | 1.252 ms | 86.815 ms | 69.362x |

#### After Step 6D

Step 6D cached `find_chiral_atom_special_cases()` once per molecule state
inside the SMILES writer path. Canonical ranking and fragment canonicalization
now share that derived ring-stereo data when they operate on the same molecule;
the `do_kekule` path still computes separate data for the kekulized temporary
molecule.

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.614 ms | 20.557 ms | 33.461x |
| `read_sdf_to_smiles` | yes | 0.804 ms | 33.015 ms | 41.073x |
| `cached_mol_to_smiles` | yes | 0.222 ms | 13.045 ms | 58.822x |
| `cached_add_remove_hs_to_smiles` | yes | 1.174 ms | 84.144 ms | 71.683x |

#### Before Step 6E

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension from Step 6D
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.553 ms | 20.454 ms | 36.969x |
| `read_sdf_to_smiles` | yes | 1.037 ms | 33.116 ms | 31.922x |
| `cached_mol_to_smiles` | yes | 0.277 ms | 12.927 ms | 46.668x |
| `cached_add_remove_hs_to_smiles` | yes | 1.178 ms | 83.112 ms | 70.549x |

#### After Step 6E

Step 6E added a lazy `RingQueryCache` inside
`find_chiral_atom_special_cases()`. The cache replaces repeated
`all_cycle_candidates()` calls only after the original cheap gates
(`chiral_tag` present and no CIP code) indicate that ring-stereo special-case
logic is actually needed. This avoids precomputing ring state for ordinary
non-ring stereocenters.

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.530 ms | 20.461 ms | 38.581x |
| `read_sdf_to_smiles` | yes | 0.786 ms | 33.970 ms | 43.232x |
| `cached_mol_to_smiles` | yes | 0.277 ms | 13.002 ms | 46.955x |
| `cached_add_remove_hs_to_smiles` | yes | 1.158 ms | 86.783 ms | 74.929x |

#### Expanded Baseline Before Step 6F

The benchmark was expanded before touching hydrogen handling so the
add/remove-H cost is visible separately from SMILES writing. Two small-molecule
COW/value-semantics cases were also added: one for topology-changing
`with_hydrogens()` and one for coordinate-only `with_2d_coords()`. Both compare
against RDKit's explicit copied-molecule workflow and assert that the original
input remains unchanged.

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py --hide-output
```

Environment:

- debug Python extension from Step 6E
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.524 ms | 20.400 ms | 38.967x |
| `read_sdf_to_smiles` | yes | 0.895 ms | 33.713 ms | 37.683x |
| `cached_mol_to_smiles` | yes | 0.233 ms | 13.692 ms | 58.840x |
| `cached_add_hs_only` | yes | 0.190 ms | 0.974 ms | 5.111x |
| `cached_remove_hs_only` | yes | 0.625 ms | 72.644 ms | 116.269x |
| `cached_add_remove_hs_only` | yes | 0.788 ms | 72.562 ms | 92.118x |
| `cached_add_hs_to_smiles` | yes | 0.894 ms | 44.613 ms | 49.913x |
| `cached_remove_hs_to_smiles` | yes | 1.076 ms | 85.926 ms | 79.856x |
| `cached_add_remove_hs_to_smiles` | yes | 1.227 ms | 85.938 ms | 70.045x |
| `cow_add_hs_preserves_input` | yes | 0.006 ms | 0.022 ms | 3.373x |
| `cow_2d_coords_preserves_input` | yes | 0.007 ms | 0.048 ms | 7.287x |

Conclusion: the next useful optimization target is `without_hydrogens()`, not
`with_hydrogens()`. On this peptide input, adding hydrogens costs about
`0.974 ms`, while removing hydrogens from the explicit-H molecule costs about
`72.644 ms` before SMILES writing.

#### After Step 6F

Step 6F optimized `without_hydrogens()` while keeping the RDKit source execution
model: hydrogens are still selected first and removed in descending atom-index
order. The hot path was repeated valence assignment during each individual H
removal. RDKit's `removeHs()` calls `updatePropertyCache(false)` once before
the removal loop and then `molRemoveH()` reads the cached total valence; COSMolKit
now mirrors that by computing total valence once and passing it through the
per-H removal path.

The physical deletion path also got a narrow trailing atom/bond fast path for
the common `AddHs()` then `RemoveHs()` layout, and the precomputed H
neighbor/bond table is reused as a validated hint before falling back to the
old bond scan.

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py --hide-output
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.541 ms | 20.388 ms | 37.699x |
| `read_sdf_to_smiles` | yes | 0.827 ms | 34.110 ms | 41.248x |
| `cached_mol_to_smiles` | yes | 0.248 ms | 13.777 ms | 55.643x |
| `cached_add_hs_only` | yes | 0.203 ms | 1.048 ms | 5.158x |
| `cached_remove_hs_only` | yes | 0.634 ms | 1.447 ms | 2.284x |
| `cached_add_remove_hs_only` | yes | 0.849 ms | 2.002 ms | 2.357x |
| `cached_add_hs_to_smiles` | yes | 0.814 ms | 44.908 ms | 55.175x |
| `cached_remove_hs_to_smiles` | yes | 1.245 ms | 14.857 ms | 11.937x |
| `cached_add_remove_hs_to_smiles` | yes | 1.291 ms | 15.456 ms | 11.973x |
| `cow_add_hs_preserves_input` | yes | 0.006 ms | 0.022 ms | 3.428x |
| `cow_2d_coords_preserves_input` | yes | 0.006 ms | 0.048 ms | 7.718x |

Key deltas from the expanded baseline:

- `cached_remove_hs_only`: `72.644 ms -> 1.447 ms`
- `cached_add_remove_hs_only`: `72.562 ms -> 2.002 ms`
- `cached_remove_hs_to_smiles`: `85.926 ms -> 14.857 ms`
- `cached_add_remove_hs_to_smiles`: `85.938 ms -> 15.456 ms`

#### Expanded Baseline Before Step 6G

The benchmark was expanded again before optimizing explicit-H SMILES writing.
The new `cached_explicit_h_mol_to_smiles` case precomputes an explicit-H
molecule and times only canonical isomeric SMILES generation. This separates
writer cost from `with_hydrogens()` cost.

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py --hide-output
```

Environment:

- debug Python extension from Step 6F
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.540 ms | 20.338 ms | 37.686x |
| `read_sdf_to_smiles` | yes | 0.890 ms | 33.842 ms | 38.043x |
| `cached_mol_to_smiles` | yes | 0.238 ms | 13.786 ms | 57.807x |
| `cached_add_hs_only` | yes | 0.199 ms | 0.987 ms | 4.955x |
| `cached_explicit_h_mol_to_smiles` | yes | 0.648 ms | 44.619 ms | 68.847x |
| `cached_remove_hs_only` | yes | 0.774 ms | 1.448 ms | 1.870x |
| `cached_add_remove_hs_only` | yes | 0.804 ms | 1.650 ms | 2.053x |
| `cached_add_hs_to_smiles` | yes | 1.050 ms | 44.961 ms | 42.819x |
| `cached_remove_hs_to_smiles` | yes | 1.034 ms | 14.630 ms | 14.147x |
| `cached_add_remove_hs_to_smiles` | yes | 1.285 ms | 15.062 ms | 11.717x |
| `cow_add_hs_preserves_input` | yes | 0.007 ms | 0.026 ms | 3.542x |
| `cow_2d_coords_preserves_input` | yes | 0.008 ms | 0.057 ms | 7.612x |

#### After Step 6G

Step 6G optimized canonical SMILES graph traversal for explicit-H molecules
without changing traversal semantics:

- `adjacency(mol)` now borrows cached adjacency instead of cloning it on every
  call, and still materializes an owned adjacency only when the molecule lacks
  a cache.
- Removed an unused canonical-ranking precheck that computed whether any bond
  was in a ring but did not feed the result into later logic.
- Added a traversal-local `BondCycleCache` so `bond_in_any_cycle()` results are
  shared between cycle discovery and stack construction for the same fragment.
- Reused the already-borrowed adjacency in chiral traversal post-processing and
  double-bond controlling-atom setup.

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py --hide-output
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.522 ms | 20.633 ms | 39.546x |
| `read_sdf_to_smiles` | yes | 0.805 ms | 27.706 ms | 34.423x |
| `cached_mol_to_smiles` | yes | 0.301 ms | 6.971 ms | 23.174x |
| `cached_add_hs_only` | yes | 0.183 ms | 1.054 ms | 5.760x |
| `cached_explicit_h_mol_to_smiles` | yes | 0.683 ms | 15.321 ms | 22.421x |
| `cached_remove_hs_only` | yes | 0.642 ms | 1.767 ms | 2.754x |
| `cached_add_remove_hs_only` | yes | 0.838 ms | 1.789 ms | 2.135x |
| `cached_add_hs_to_smiles` | yes | 0.828 ms | 15.523 ms | 18.742x |
| `cached_remove_hs_to_smiles` | yes | 1.044 ms | 8.285 ms | 7.935x |
| `cached_add_remove_hs_to_smiles` | yes | 1.172 ms | 8.656 ms | 7.385x |
| `cow_add_hs_preserves_input` | yes | 0.006 ms | 0.021 ms | 3.315x |
| `cow_2d_coords_preserves_input` | yes | 0.007 ms | 0.048 ms | 7.142x |

Key deltas from the Step 6G baseline:

- `cached_mol_to_smiles`: `13.786 ms -> 6.971 ms`
- `cached_explicit_h_mol_to_smiles`: `44.619 ms -> 15.321 ms`
- `cached_add_hs_to_smiles`: `44.961 ms -> 15.523 ms`
- `cached_remove_hs_to_smiles`: `14.630 ms -> 8.285 ms`
- `cached_add_remove_hs_to_smiles`: `15.062 ms -> 8.656 ms`

#### After Step 6H

Step 6H split the benchmark surface and optimized the SDF/3D stereo path:

- Added Python `Molecule.num_atoms()` / `Molecule.num_bonds()` so signature
  benchmarks no longer pay for Python `atoms()` feature construction.
- Added record-text, non-isomeric SMILES, explicit-H non-isomeric SMILES, and
  COW cases to separate parser, writer, stereo, and COW costs.
- SDF V2000/V3000 parsing now fills a new molecule's atom/bond tables directly
  instead of calling `add_atom()` / `add_bond()` for every parsed line.
- Molfile explicit valence is computed once and reused by total-valence and
  aromaticity setup.
- 3D MolBlock stereo finalization and legacy stereo property calculation now
  use local adjacency and bond-cycle caches instead of repeated full-bond scans
  and repeated ring BFS. The RDKit-aligned decision branches are unchanged.
- `rdkit_legacy_stereo_atom_props(false)` no longer computes potential stereo
  centers/bonds that only affect the `flag_possible_stereo_centers=true` output.

Command:

```bash
.venv/bin/python benchmark_sdf_to_smiles.py --hide-output
```

Environment:

- debug Python extension rebuilt with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
- input: `环肽_9hyt.sdf`
- iterations: 120
- warmup: 10

Results:

| Case | Same Output | RDKit Mean | COSMolKit Mean | COSMolKit / RDKit Mean |
| --- | --- | ---: | ---: | ---: |
| `read_sdf_signature` | yes | 0.541 ms | 10.938 ms | 20.202x |
| `read_sdf_record_from_str_signature` | yes | 0.646 ms | 10.689 ms | 16.538x |
| `read_sdf_to_smiles` | yes | 0.816 ms | 15.045 ms | 18.440x |
| `read_sdf_to_nonisomeric_smiles` | yes | 0.868 ms | 12.467 ms | 14.366x |
| `cached_mol_to_smiles` | yes | 0.240 ms | 3.867 ms | 16.097x |
| `cached_mol_to_nonisomeric_smiles` | yes | 0.212 ms | 1.578 ms | 7.444x |
| `cached_add_hs_only` | yes | 0.184 ms | 0.258 ms | 1.405x |
| `cached_explicit_h_mol_to_smiles` | yes | 0.838 ms | 11.959 ms | 14.272x |
| `cached_explicit_h_mol_to_nonisomeric_smiles` | yes | 0.555 ms | 5.018 ms | 9.044x |
| `cached_remove_hs_only` | yes | 0.641 ms | 1.144 ms | 1.786x |
| `cached_add_remove_hs_only` | yes | 0.786 ms | 1.632 ms | 2.075x |
| `cached_add_hs_to_smiles` | yes | 0.843 ms | 12.332 ms | 14.628x |
| `cached_remove_hs_to_smiles` | yes | 1.097 ms | 4.936 ms | 4.499x |
| `cached_add_remove_hs_to_smiles` | yes | 1.453 ms | 5.323 ms | 3.665x |
| `cow_add_hs_preserves_input` | yes | 0.007 ms | 0.006 ms | 0.966x |
| `cow_2d_coords_preserves_input` | yes | 0.006 ms | 0.049 ms | 7.626x |

Key deltas from the expanded Step 6H pre-optimization baseline:

- `read_sdf_signature`: `20.500 ms -> 10.938 ms`
- `read_sdf_to_smiles`: `27.159 ms -> 15.045 ms`
- `read_sdf_to_nonisomeric_smiles`: `21.922 ms -> 12.467 ms`
- `cached_mol_to_smiles`: `6.744 ms -> 3.867 ms`
- `cached_remove_hs_to_smiles`: `8.083 ms -> 4.936 ms`
- `cached_add_remove_hs_to_smiles`: `8.280 ms -> 5.323 ms`

Current interpretation:

- `AddHs()` / `RemoveHs()` count-only paths are close enough in debug mode that
  further work should be justified by a concrete profile.
- The remaining SDF parser gap is mostly baseline parsing/aromaticity/valence
  plus a smaller 3D stereo cost: forcing `coordinate_dim="2d"` reads this input
  at roughly `8.2-9.1 ms` median while `auto/3d` is roughly `9.7-10.3 ms`.
- Isomeric writer cost is still higher than non-isomeric writer cost, but the
  legacy stereo graph-query portion is no longer the dominant avoidable cost.
