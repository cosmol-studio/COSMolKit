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
