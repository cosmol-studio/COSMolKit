# Topology Operation Abstraction Sketch

This document describes a lightweight implementation strategy for handling topology-related molecule operations in COSMolKit.

The goal is not to create a large framework.  
The goal is to make topology-changing operations safer, more consistent, and easier to test, while keeping runtime overhead low.

---

## 1. Core Problem

Some molecule operations change atom/bond topology or chemical graph state.

Examples:

```text
without_hydrogens
with_hydrogens
remove_atom
renumber_atoms
kekulize
sanitize
set_bond_type
```

These operations may affect dependent data:

```text
coordinates
stereochemistry
atom-level properties
bond-level properties
molecule-level properties
adjacency cache
ring cache
valence cache
drawing cache
fingerprint cache
```

The main risk is not that the atom list itself is wrong.

The main risk is that dependent data becomes stale or misaligned.

Example:

```text
Atom was removed,
but coordinate row was not removed.

Atom index changed,
but stereo still references old index.

Bond state changed,
but cache still stores old graph state.
```

So we want a lightweight internal protocol:

```text
Strong operation:
  index system may change -> use mapping

Weak operation:
  index system stays stable -> invalidate/recompute derived state
```

---

## 2. Strong vs Weak Topology Operations

### Strong topology-changing operation

A strong operation changes one or more of:

```text
atom count
bond count
atom table ordering
bond table ordering
atom index mapping
bond index mapping
```

Examples:

```text
without_hydrogens
with_hydrogens
remove_atom
add_atom
remove_bond
add_bond
renumber_atoms
fragment
combine
```

Strong operations need explicit index mapping.

---

### Weak topology-state operation

A weak operation keeps atom and bond indices stable, but changes chemical graph state.

Examples:

```text
kekulize
sanitize
set_bond_type
set_formal_charge
set_aromaticity
update_implicit_hydrogens
```

Weak operations usually do not need atom mapping, but they must invalidate or recompute affected derived state.

---

## 3. IndexMapping

For strong operations, define a lightweight internal mapping type.

This mapping records how old indices correspond to new indices.

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct IndexMapping {
    old_to_new: Vec<Option<usize>>,
    new_to_old: Vec<Option<usize>>,
}
```

Use `Option<usize>` because strong operations may delete or create atoms.

For deletion:

```text
old atom exists, but no new atom
old_to_new[old_idx] = None
```

For creation:

```text
new atom exists, but no old atom
new_to_old[new_idx] = None
```

Example: `without_hydrogens`

```text
old atoms:
0 C
1 C
2 O
3 H
4 H
5 H

new atoms:
0 C <- old 0
1 C <- old 1
2 O <- old 2
```

Mapping:

```rust
old_to_new = vec![
    Some(0),
    Some(1),
    Some(2),
    None,
    None,
    None,
];

new_to_old = vec![
    Some(0),
    Some(1),
    Some(2),
];
```

Example: `with_hydrogens`

```text
old atoms:
0 C
1 O

new atoms:
0 C <- old 0
1 O <- old 1
2 H <- newly created
3 H <- newly created
4 H <- newly created
```

Mapping:

```rust
old_to_new = vec![
    Some(0),
    Some(1),
];

new_to_old = vec![
    Some(0),
    Some(1),
    None,
    None,
    None,
];
```

---

## 4. Minimal IndexMapping API

Keep this type small.

Do not over-abstract it into a large trait framework.

```rust
impl IndexMapping {
    pub(crate) fn new(
        old_to_new: Vec<Option<usize>>,
        new_to_old: Vec<Option<usize>>,
    ) -> Self {
        Self {
            old_to_new,
            new_to_old,
        }
    }

    pub(crate) fn old_to_new(&self, old_idx: usize) -> Option<usize> {
        self.old_to_new.get(old_idx).copied().flatten()
    }

    pub(crate) fn new_to_old(&self, new_idx: usize) -> Option<usize> {
        self.new_to_old.get(new_idx).copied().flatten()
    }

    pub(crate) fn old_len(&self) -> usize {
        self.old_to_new.len()
    }

    pub(crate) fn new_len(&self) -> usize {
        self.new_to_old.len()
    }

    pub(crate) fn is_deleted_old(&self, old_idx: usize) -> bool {
        self.old_to_new(old_idx).is_none()
    }

    pub(crate) fn is_newly_created(&self, new_idx: usize) -> bool {
        self.new_to_old(new_idx).is_none()
    }
}
```

A helper for filter-style operations:

```rust
impl IndexMapping {
    pub(crate) fn from_kept_indices(old_len: usize, kept_old_indices: &[usize]) -> Self {
        let mut old_to_new = vec![None; old_len];
        let mut new_to_old = Vec::with_capacity(kept_old_indices.len());

        for (new_idx, &old_idx) in kept_old_indices.iter().enumerate() {
            old_to_new[old_idx] = Some(new_idx);
            new_to_old.push(Some(old_idx));
        }

        Self {
            old_to_new,
            new_to_old,
        }
    }
}
```

This is useful for operations like:

```text
without_hydrogens
remove_atom
filter_atoms
fragment extraction
```

---

## 5. AtomMapping and BondMapping

For readability, we can use type aliases:

```rust
pub(crate) type AtomMapping = IndexMapping;
pub(crate) type BondMapping = IndexMapping;
```

If atom mapping and bond mapping later need different behavior, they can become separate structs.

For now, aliases are enough.

---

## 6. TopologyOperationResult

Strong operations may internally return a result that includes the molecule and mapping.

```rust
#[derive(Debug)]
pub(crate) struct TopologyOperationResult {
    pub molecule: Molecule,
    pub atom_mapping: Option<AtomMapping>,
    pub bond_mapping: Option<BondMapping>,
}
```

This does not need to be public Python API.

It is mainly useful internally and for tests.

---

## 7. Strong Operation Implementation Pipeline

Strong operations must not hand-code the state-migration pipeline as free
steps inside each operation body. The operation body decides chemistry
selection, but compacting topology edits must go through a dedicated `OpParts`
API such as `remove_atoms(...)`.

The dedicated API owns this pipeline:

```text
1. Decide which atoms/bonds are kept, removed, created, or reordered.
2. Build AtomMapping and BondMapping.
3. Rebuild topology using the mapping.
4. Remap conformers using AtomMapping.
5. Remap stereochemistry using AtomMapping/BondMapping.
6. Remap, drop, or reject atom-level and bond-level properties.
7. Preserve molecule-level properties unless operation documents otherwise.
8. Invalidate all affected caches.
9. Return new Molecule.
```

Operation bodies should look like:

```rust
#[mol_op_body(without_hydrogens, parts)]
fn without_hydrogens_impl(sanitize: bool) -> Result<OpOutcome, OperationError> {
    let removable = find_removable_hydrogens(parts.topology(), sanitize)?;
    if removable.is_empty() {
        return Ok(OpOutcome::NoOp { reason: "no removable hydrogens" });
    }

    parts.remove_atoms(&removable)?;
    parts.recompute_valence()?;
    parts.recompute_aromaticity()?;
    parts.recompute_stereo()?;
    Ok(OpOutcome::Changed)
}
```

The operation body must not manually interleave:

```text
topology_mut().remove...
coordinates_mut().remap...
properties_mut().remap...
set mapping report...
invalidate...
```

Those steps are registry-controlled state migration and belong in `OpParts`.

---

## 8. Remapping Conformers

Conformer remapping should use `new_to_old`.

That makes the logic direct:

```rust
pub(crate) fn remap_conformer_rows<T: Clone>(
    old_rows: &[T],
    atom_mapping: &AtomMapping,
) -> Vec<T> {
    let mut new_rows = Vec::with_capacity(atom_mapping.new_len());

    for new_idx in 0..atom_mapping.new_len() {
        match atom_mapping.new_to_old(new_idx) {
            Some(old_idx) => {
                new_rows.push(old_rows[old_idx].clone());
            }
            None => {
                // Newly created atom.
                // The caller must define how to initialize its coordinate.
                panic!("newly created atom requires coordinate initialization");
            }
        }
    }

    new_rows
}
```

For `without_hydrogens`, every new atom has an old source.

For `with_hydrogens`, newly created H atoms need explicit coordinate initialization.

So add-H should not blindly call this helper unless it also provides coordinates for new atoms.

A more general API:

```rust
pub(crate) fn remap_conformer_rows_with_new<T: Clone>(
    old_rows: &[T],
    atom_mapping: &AtomMapping,
    mut init_new: impl FnMut(usize) -> T,
) -> Vec<T> {
    let mut new_rows = Vec::with_capacity(atom_mapping.new_len());

    for new_idx in 0..atom_mapping.new_len() {
        match atom_mapping.new_to_old(new_idx) {
            Some(old_idx) => new_rows.push(old_rows[old_idx].clone()),
            None => new_rows.push(init_new(new_idx)),
        }
    }

    new_rows
}
```

---

## 9. Remapping Stereo

Stereo is high-risk.

A stereo object that references atom indices must either be:

```text
remapped
recomputed
dropped explicitly
or rejected with an error
```

Never leave stereo pointing to deleted atoms.

Example skeleton:

```rust
pub(crate) fn remap_tetrahedral_stereo(
    stereo: &TetrahedralStereo,
    atom_mapping: &AtomMapping,
) -> Option<TetrahedralStereo> {
    let new_center = atom_mapping.old_to_new(stereo.center)?;

    let mut new_ligands = Vec::with_capacity(stereo.ligands.len());

    for old_ligand in &stereo.ligands {
        match old_ligand {
            Some(old_idx) => {
                let new_idx = atom_mapping.old_to_new(*old_idx)?;
                new_ligands.push(Some(new_idx));
            }
            None => {
                new_ligands.push(None);
            }
        }
    }

    Some(TetrahedralStereo {
        center: new_center,
        ligands: new_ligands,
        // Other fields must be preserved or recomputed carefully.
        ..stereo.clone()
    })
}
```

If remapping would make stereo invalid, return `None` or return an explicit `StereoError`, depending on operation policy.

For high-risk cases, prefer explicit error until behavior is verified against RDKit.

---

## 10. Weak Operation Implementation Pipeline

Weak operations do not need index mapping.

They should follow this pipeline:

```text
1. Clone molecule using value-style behavior.
2. Mutate topology state in the new molecule.
3. Keep atom/bond indices stable.
4. Preserve coordinates unless documented otherwise.
5. Preserve molecule-level properties unless documented otherwise.
6. Invalidate or recompute all affected derived state.
7. Return new Molecule.
```

Example: `kekulize`

```rust
pub(crate) fn kekulize_impl(mol: &Molecule) -> Result<Molecule> {
    let mut out = mol.clone();

    {
        let topology = out.make_topology_mut();
        kekulize_topology_in_place(topology)?;

        topology.invalidate_derived_state(
            DerivedState::AROMATICITY
                | DerivedState::VALENCE
                | DerivedState::DRAWING
                | DerivedState::FINGERPRINT,
        );
    }

    Ok(out)
}
```

The exact invalidated states should match the actual cache system.

---

## 11. DerivedState

Use a lightweight bitflag to describe derived state affected by operations.

```rust
bitflags::bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub(crate) struct DerivedState: u32 {
        const ADJACENCY   = 1 << 0;
        const RINGS       = 1 << 1;
        const VALENCE     = 1 << 2;
        const AROMATICITY = 1 << 3;
        const STEREO      = 1 << 4;
        const DRAWING     = 1 << 5;
        const FINGERPRINT = 1 << 6;
    }
}
```

Each operation declares what it invalidates.

Example:

```rust
const WITHOUT_HYDROGENS_INVALIDATES: DerivedState =
    DerivedState::ADJACENCY
        .union(DerivedState::RINGS)
        .union(DerivedState::VALENCE)
        .union(DerivedState::STEREO)
        .union(DerivedState::DRAWING)
        .union(DerivedState::FINGERPRINT);
```

Or if using bitflags operators:

```rust
const WITHOUT_HYDROGENS_INVALIDATES: DerivedState =
    DerivedState::ADJACENCY
        | DerivedState::RINGS
        | DerivedState::VALENCE
        | DerivedState::STEREO
        | DerivedState::DRAWING
        | DerivedState::FINGERPRINT;
```

---

## 12. TopologyOpSpec

Define operation metadata.

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TopologyOpKind {
    Strong,
    Weak,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DependentData {
    Coordinates,
    Stereo,
    MoleculeProps,
    AtomProps,
    BondProps,
    AdjacencyCache,
    RingCache,
    ValenceCache,
    AromaticityState,
    DrawingState,
    FingerprintState,
}

#[derive(Debug)]
pub(crate) struct TopologyOpSpec {
    pub name: &'static str,
    pub kind: TopologyOpKind,
    pub affects: &'static [DependentData],
    pub invalidates: DerivedState,
    pub support: SupportStatus,
    pub parity: ParityPolicy,
    pub io_roundtrip: bool,
}
```

This document predates the macro-controlled operation registry. In the current
architecture, topology operation specs should be generated from
`molecule_ops!`, not hand-maintained as scattered static values. The fields
above map to the generated `MoleculeOpSpec` shape:

```text
kind -> MoleculeOpKind
affects -> may_mutate / must_handle
invalidates -> invalidates
topology edit style -> topology_edit
required remaps -> auto_remap
mapping/report contract -> requires_mapping / report
support -> SupportStatus from FeatureSpec
parity -> ParityPolicy
```

`rdkit_parity: bool` is no longer allowed because it conflates "should match
RDKit once supported" with "currently claims RDKit parity". Use:

```text
ParityPolicy::NotApplicable
ParityPolicy::RequiredWhenSupported
ParityPolicy::RequiredNow
```

Example:

```rust
pub(crate) static WITHOUT_HYDROGENS_SPEC: TopologyOpSpec = TopologyOpSpec {
    name: "without_hydrogens",
    kind: TopologyOpKind::Strong,
    affects: &[
        DependentData::Coordinates,
        DependentData::Stereo,
        DependentData::MoleculeProps,
        DependentData::AdjacencyCache,
        DependentData::RingCache,
        DependentData::ValenceCache,
    ],
    invalidates: DerivedState::ADJACENCY
        | DerivedState::RINGS
        | DerivedState::VALENCE
        | DerivedState::STEREO
        | DerivedState::DRAWING
        | DerivedState::FINGERPRINT,
    support: SupportStatus::Unsupported {
        reason: "hydrogen removal has not been ported",
    },
    parity: ParityPolicy::RequiredWhenSupported,
    io_roundtrip: true,
};
```

Example weak operation:

```rust
pub(crate) static KEKULIZE_SPEC: TopologyOpSpec = TopologyOpSpec {
    name: "kekulize",
    kind: TopologyOpKind::Weak,
    affects: &[
        DependentData::AromaticityState,
        DependentData::ValenceCache,
        DependentData::DrawingState,
        DependentData::FingerprintState,
    ],
    invalidates: DerivedState::AROMATICITY
        | DerivedState::VALENCE
        | DerivedState::DRAWING
        | DerivedState::FINGERPRINT,
    support: SupportStatus::Unsupported {
        reason: "kekulization has not been ported",
    },
    parity: ParityPolicy::RequiredWhenSupported,
    io_roundtrip: true,
};
```

---

## 13. Macro for Operation Spec

A macro can reduce boilerplate.

This macro should only generate metadata.

It should not hide core molecule transformation logic.

```rust
macro_rules! topology_op_spec {
    (
        $vis:vis static $ident:ident = {
            name: $name:literal,
            kind: $kind:ident,
            affects: [$($affects:ident),* $(,)?],
            invalidates: $invalidates:expr,
            support: $support:expr,
            parity: $parity:expr,
            io_roundtrip: $io:expr $(,)?
        };
    ) => {
        $vis static $ident: TopologyOpSpec = TopologyOpSpec {
            name: $name,
            kind: TopologyOpKind::$kind,
            affects: &[$(DependentData::$affects),*],
            invalidates: $invalidates,
            support: $support,
            parity: $parity,
            io_roundtrip: $io,
        };
    };
}
```

Usage:

```rust
topology_op_spec! {
    pub(crate) static WITHOUT_HYDROGENS_SPEC = {
        name: "without_hydrogens",
        kind: Strong,
        affects: [
            Coordinates,
            Stereo,
            MoleculeProps,
            AdjacencyCache,
            RingCache,
            ValenceCache,
        ],
        invalidates: DerivedState::ADJACENCY
            | DerivedState::RINGS
            | DerivedState::VALENCE
            | DerivedState::STEREO
            | DerivedState::DRAWING
            | DerivedState::FINGERPRINT,
        support: SupportStatus::Unsupported {
            reason: "hydrogen removal has not been ported",
        },
        parity: ParityPolicy::RequiredWhenSupported,
        io_roundtrip: true,
    };
}
```

---

## 14. Operation Wrapper

Use a generic wrapper for debug-only invariant checks.

This avoids dynamic dispatch.

```rust
#[inline]
pub(crate) fn apply_topology_op<F>(
    mol: &Molecule,
    spec: &'static TopologyOpSpec,
    f: F,
) -> Result<Molecule>
where
    F: FnOnce(&Molecule) -> Result<Molecule>,
{
    let out = f(mol)?;

    #[cfg(debug_assertions)]
    {
        crate::testing::check_topology_op_contract(mol, &out, spec)?;
    }

    Ok(out)
}
```

Example:

```rust
pub fn without_hydrogens(&self) -> Result<Molecule> {
    apply_topology_op(self, &WITHOUT_HYDROGENS_SPEC, without_hydrogens_public_impl)
}
```

For strong operations that internally return mappings:

```rust
#[inline]
pub(crate) fn apply_strong_topology_op<F>(
    mol: &Molecule,
    spec: &'static TopologyOpSpec,
    f: F,
) -> Result<Molecule>
where
    F: FnOnce(&Molecule) -> Result<TopologyOperationResult>,
{
    debug_assert_eq!(spec.kind, TopologyOpKind::Strong);

    let result = f(mol)?;
    let out = result.molecule;

    #[cfg(debug_assertions)]
    {
        crate::testing::check_strong_topology_op_contract(
            mol,
            &out,
            spec,
            result.atom_mapping.as_ref(),
            result.bond_mapping.as_ref(),
        )?;
    }

    Ok(out)
}
```

Usage:

```rust
pub fn without_hydrogens(&self) -> Result<Molecule> {
    apply_strong_topology_op(self, &WITHOUT_HYDROGENS_SPEC, without_hydrogens_impl)
}
```

---

## 15. Optional Macro for Operation Wrapper

If the wrapper boilerplate becomes repetitive, use a small macro.

```rust
macro_rules! define_strong_topology_method {
    (
        $vis:vis fn $name:ident(&self) -> Result<Molecule>
        spec = $spec:ident;
        body = $body:ident;
    ) => {
        #[inline]
        $vis fn $name(&self) -> Result<Molecule> {
            apply_strong_topology_op(self, &$spec, $body)
        }
    };
}
```

Usage:

```rust
impl Molecule {
    define_strong_topology_method! {
        pub fn without_hydrogens(&self) -> Result<Molecule>
        spec = WITHOUT_HYDROGENS_SPEC;
        body = without_hydrogens_impl;
    }
}
```

This macro only removes repetitive wrapper code.

It should not contain chemistry logic.

---

## 16. Why Not Pure Macro for Everything?

Do not put the actual molecule transformation logic inside macros.

Avoid this:

```rust
strong_topology_op! {
    without_hydrogens {
        // 200 lines of atom filtering, bond rebuilding, stereo handling
    }
}
```

Reasons:

```text
harder to debug
worse compiler errors
harder to unit test pieces
harder to profile
harder for future contributors to read
```

Preferred split:

```text
normal Rust functions:
  build mapping
  rebuild topology
  remap coordinates
  remap stereo
  invalidate caches

macros:
  declare operation metadata
  generate repetitive wrappers
  generate test matrix boilerplate
```

---

## 17. Performance Notes

The recommended design avoids expensive runtime abstraction.

Avoid:

```rust
Vec<Box<dyn TopologyOperation>>
```

Avoid unnecessary dynamic dispatch.

Prefer:

```text
static specs
generic wrappers
#[inline]
debug-only invariant checks
macro-generated boilerplate
normal optimized Rust loops
```

`IndexMapping` has runtime cost, but it represents real runtime information.

This cost cannot be removed by macros because mapping depends on the input molecule.

Mapping should be temporary and internal:

```rust
let mapping = build_mapping(...);
let out = rebuild_with_mapping(...);
// mapping is dropped here
```

The expected cost is O(n_atoms) or O(n_bonds), which is appropriate for topology-changing operations.

---

## 18. Possible Optimization Later

If `Vec<Option<usize>>` becomes a benchmarked bottleneck, use a sentinel-based representation.

```rust
const ABSENT: usize = usize::MAX;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CompactIndexMapping {
    old_to_new: Vec<usize>,
    new_to_old: Vec<usize>,
}
```

Then:

```rust
impl CompactIndexMapping {
    pub(crate) fn old_to_new(&self, old_idx: usize) -> Option<usize> {
        let value = self.old_to_new[old_idx];
        if value == ABSENT {
            None
        } else {
            Some(value)
        }
    }
}
```

Do not start with this unless profiling shows mapping overhead matters.

The simple `Option<usize>` version is easier to read and safer.

---

## 19. Test Matrix Macro

Test boilerplate can be generated with macros.

Example:

```rust
macro_rules! topology_op_invariant_tests {
    (
        $(
            $test_name:ident:
            spec = $spec:ident,
            corpus = $corpus:expr,
            apply = $apply:expr;
        )*
    ) => {
        $(
            #[test]
            fn $test_name() {
                for case in $corpus {
                    let mol = Molecule::from_smiles(case.smiles).unwrap();

                    let snapshot = MoleculeSnapshot::capture(&mol);

                    force_all_relevant_caches(&mol);

                    let out = $apply(&mol).unwrap();

                    assert_source_unchanged(&snapshot, &mol);
                    assert_topology_op_contract(&mol, &out, &$spec).unwrap();
                }
            }
        )*
    };
}
```

Usage:

```rust
topology_op_invariant_tests! {
    invariant_without_hydrogens:
        spec = WITHOUT_HYDROGENS_SPEC,
        corpus = topology_core_cases(),
        apply = |mol: &Molecule| mol.with_hydrogens()?.without_hydrogens();

    invariant_kekulize:
        spec = KEKULIZE_SPEC,
        corpus = aromatic_cases(),
        apply = |mol: &Molecule| mol.with_kekulized_bonds();
}
```

The macro should not know chemistry.

It should only run the operation and invoke invariant checks.

---

## 20. Recommended File Layout

Suggested internal structure:

```text
crates/cosmolkit-core/src/
  topology/
    mod.rs
    mapping.rs
    ops.rs
    remap.rs
    derived_state.rs

  testing/
    topology_contract.rs
    invariants.rs
```

Or, if you want to keep it smaller at first:

```text
crates/cosmolkit-core/src/
  topology_mapping.rs
  topology_ops.rs
```

Test-side support:

```text
crates/cosmolkit-core/tests/
  support/
    topology_cases.rs
    topology_invariants.rs
    topology_ops.rs

  topology_operation_invariants.rs
```

---

## 21. Minimal First Implementation

Do not implement everything at once.

Start with this minimum:

```text
1. IndexMapping
2. TopologyOpKind
3. TopologyOpSpec
4. DerivedState
5. WITHOUT_HYDROGENS_SPEC
6. KEKULIZE_SPEC
7. apply_strong_topology_op()
8. apply_topology_op()
9. basic invariant tests for without_hydrogens and kekulize
```

Initial checks:

For strong operation:

```text
source molecule unchanged
graph valid
coordinate rows match atom count
stereo indices valid
```

For weak operation:

```text
source molecule unchanged
graph valid
atom count unchanged
bond count unchanged
coordinates unchanged unless documented otherwise
stereo indices valid
```

Add cache and property checks after the basic runner is stable.

---

## 22. Key Design Rule

The mapping is not the abstraction goal.

The mapping is the shared factual record of index migration.

The real goal is:

```text
All dependent data must update from the same index migration record.
```

This prevents:

```text
topology using one index mapping
coordinates using another
stereo forgetting the mapping
props silently misaligning
cache using old graph state
```

---

## 23. Summary

Recommended design:

```text
IndexMapping:
  lightweight internal runtime data structure

Strong operation:
  build mapping -> rebuild topology -> remap dependent data -> invalidate caches

Weak operation:
  mutate chemical graph state -> keep indices stable -> invalidate/recompute derived state

Macros:
  declare operation specs
  generate wrapper boilerplate
  generate invariant test boilerplate

Do not:
  use dyn operation objects
  expose mapping to Python
  hide chemistry logic inside macros
  over-design trait frameworks early
```

Short version:

```text
strong = index mapping
weak = cache invalidation / derived-state refresh

mapping is runtime fact
macros are compile-time boilerplate reduction
```
