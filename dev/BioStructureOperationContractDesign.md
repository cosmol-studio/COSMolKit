# COSMolKit BioStructure Operation Contract Design

## Draft Status

This document is a design contract draft. Sections marked **Provisional**
(暂定) are expected to be revisited when the first BioStructure parser and
operation implementation lands. Provisional sections are still binding for
agents during draft-driven work: do not silently implement a conflicting
design. If an implementation needs to diverge, update this document or ask the
human author to confirm the design exception first.

## Design Goal

`BioStructure` is the unified structural representation for:

- proteins
- DNA
- RNA
- protein–ligand complexes
- protein–DNA/RNA complexes
- biological assemblies
- multi-model structures

The operation system must guarantee:

- hierarchy consistency
- coordinate alignment
- assembly correctness
- altloc/disorder correctness
- stable row-id semantics
- future extensibility toward DNA/RNA-specific operations

while remaining:

- cache-friendly
- ndarray-friendly
- GPU-friendly
- graph-friendly
- compatible with molecule extraction workflows

---

# Core Philosophy

```text
MoleculeOps:
    chemical graph correctness

BioStructureOps:
    structural hierarchy correctness
```

`BioStructure` is NOT a giant `Molecule`.

Protein/DNA/RNA structures are primarily:

- hierarchy objects
- coordinate objects
- assembly objects

rather than fully sanitized chemistry graphs.

---

# Top-Level Object

```rust
pub struct BioStructure {
    atoms: AtomBlock,
    residues: ResidueBlock,
    chains: ChainBlock,
    entities: EntityBlock,
    models: ModelBlock,

    coordinates: CoordinateBlock,
    bonds: BondBlock,

    assemblies: AssemblyBlock,
    annotations: AnnotationBlock,

    derived: BioDerivedCacheBlock,
    props: PropertyBlock,
}
```

---

# Structural Hierarchy

```text
Model
  └── Chain
        └── Residue
              └── Atom
```

All hierarchy relationships must use stable row-id references.

---

# Flat Hierarchy Storage Contract (Provisional)

COSMolKit should store BioStructure hierarchy as flat row blocks, not nested
vectors. This is intentionally different from Gemmi's source model
(`Structure -> Model -> Chain -> Residue -> Atom`) and is chosen for stable ids,
ndarray export, GPU-friendly access, and registry-controlled row migration.

Provisional internal shape:

```rust
pub struct AtomRow {
    pub residue_id: ResidueId,
    pub name: AtomName,
    pub element: Element,
    pub altloc: Option<AltLocLabel>,
    pub source: AtomSourceIds,
}

pub struct ResidueRow {
    pub chain_id: ChainId,
    pub atom_span: RowSpan<AtomId>,
    pub name: ResidueName,
    pub kind: ResidueKind,
    pub source: ResidueSourceIds,
}

pub struct ChainRow {
    pub model_id: ModelId,
    pub entity_id: Option<EntityId>,
    pub residue_span: RowSpan<ResidueId>,
    pub kind: ChainKind,
    pub source: ChainSourceIds,
}

pub struct ModelRow {
    pub chain_span: RowSpan<ChainId>,
    pub source_model_number: Option<i32>,
}
```

Provisional hierarchy invariant:

```text
atoms are grouped by residue
residues are grouped by chain
chains are grouped by model
row spans are contiguous and valid
child.parent_id agrees with the parent span that contains it
```

If a future implementation chooses non-contiguous parent membership for
streaming or incremental editing, it must update this section and add an index
cache contract before implementation. Do not mix contiguous spans and
non-contiguous membership implicitly.

---

# Stable Row IDs

```rust
AtomId
ResidueId
ChainId
EntityId
ModelId
AssemblyId
BondId
AltLocGroupId
```

Core invariant:

```text
within one immutable BioStructure snapshot:
    row index == stable id
```

Stable ids are snapshot-local identities. A strong operation creates a new
row-id space when it removes, inserts, reorders, expands, splits, or merges
rows. The relationship between the old and new id spaces must be represented by
the operation mapping.

Persistent PDB/mmCIF identifiers such as author chain id, label asym id,
residue sequence number, insertion code, segment id, atom serial, entity id, and
assembly id are source identifiers or annotations. They must not be used as
COSMolKit row ids.

---

# Polymer Types

This taxonomy is provisional. It is intended to be good enough for the first
BioStructure parser and selection APIs, but residue/polymer/chain classification
will likely need refinement when CCD, modified residues, branched glycans, and
hybrid complexes are implemented.

## ResidueKind

```rust
pub enum ResidueKind {
    AminoAcid,
    DNA,
    RNA,
    Saccharide,
    Water,
    Ligand,
    Ion,
    Cofactor,
    Unknown,
}
```

## PolymerKind

```rust
pub enum PolymerKind {
    Peptide,
    DNA,
    RNA,
    PeptideLike,
    NucleicAcidHybrid,
    Saccharide,
    NonPolymer,
    Water,
    Unknown,
}
```

## ChainKind

```rust
pub enum ChainKind {
    Protein,
    DNA,
    RNA,
    ProteinDNAComplex,
    ProteinRNAComplex,
    LigandOnly,
    WaterOnly,
    Mixed,
    Unknown,
}
```

Protein/DNA/RNA are annotations, NOT separate structs.

---

# Source Identifier Semantics

Gemmi is the planned source-level reproduction target for PDB/mmCIF
macromolecular parsing. COSMolKit's internal representation should be flat,
typed, and row-id based, but it must preserve the source identifiers that make
PDB/mmCIF structures addressable.

The following are not COSMolKit row ids:

```text
PDB atom serial
PDB chain id
PDB segment id
auth_asym_id
label_asym_id / subchain id
auth_seq_id
label_seq_id
insertion code
label_entity_id / entity id
assembly id
operator id
connection id
```

They should be represented as typed fields, typed annotations, or explicit
source-provenance metadata. They should not be flattened into generic string
properties if later operations need to validate, remap, or compare them.

The parser contract must distinguish:

```text
COSMolKit row id       stable id inside the current snapshot
source identifier      author/spec/reference identifier from PDB/mmCIF
selection address      user-facing address assembled from source identifiers
```

Gemmi-compatible parsing tests should compare a documented structured dump of
these fields rather than relying on incidental internal row ordering alone.

---

# Core Blocks

```rust
pub struct BioBlockSet(u32);

impl BioBlockSet {
    pub const ATOMS: Self;
    pub const RESIDUES: Self;
    pub const CHAINS: Self;
    pub const ENTITIES: Self;
    pub const MODELS: Self;

    pub const COORDINATES: Self;
    pub const BONDS: Self;

    pub const ASSEMBLIES: Self;
    pub const ANNOTATIONS: Self;

    pub const DERIVED_CACHE: Self;
    pub const PROPERTIES: Self;
}
```

---

# Coordinate and Multi-Model Contract (Provisional)

The initial BioStructure design should treat atom rows as model-scoped rows:
each atom belongs to exactly one residue, chain, and model, and each atom has at
most one primary coordinate row in the active coordinate block.

Provisional coordinate invariant:

```text
coordinates.len() == atoms.len()
coordinates[atom_id] belongs to atom_id
atom.model_id is derived through residue -> chain -> model
```

This matches Gemmi's practical model where each `Model` owns its own chains,
residues, and atoms. It avoids forcing different models to share one hierarchy
when PDB/mmCIF inputs may have model-specific atom presence, disorder, waters,
or ligands.

Operations that compare or align models should use explicit model-aware
selection APIs rather than assuming that atom row `N` in model 1 corresponds to
atom row `N` in model 2.

Alternative designs, such as shared hierarchy plus coordinate tensor
`[model][atom]`, remain possible but are not the current default. Switching to
that model would require revising:

```text
keep_model
split_by_model
coordinate row alignment
assembly expansion
Gemmi structured dumps
BioStructureMapping
```

---

# BioStructureOpSpec

```rust
pub struct BioStructureOpSpec {
    pub method: &'static str,
    pub impl_fn: &'static str,

    pub domain: BioOpDomain,
    pub kind: BioOpKind,
    pub edit_kind: BioEditKind,

    pub may_mutate: BioBlockSet,
    pub auto_remap: BioBlockSet,

    pub must_handle: BioStateSet,
    pub needs_update: BioDerivedState,

    pub requires_mapping: MappingRequirement,
    pub report: BioOperationReportSet,

    pub allows_noop: bool,
    pub support: SupportStatus,

    pub parity: BioParityPolicy,
    pub io_roundtrip: bool,
}
```

`BioStructureOpSpec` should mirror the current `MoleculeOpSpec` runtime shape:
it describes the operation contract, permissions, support state, and parity
policy. Invariant and parity profiles are declared in `bio_structure_ops!`
registry entries and emitted into generated matrices rather than stored directly
as runtime spec fields.

Generated source-of-truth tables should follow the molecule operation pattern:

```text
BIO_STRUCTURE_OPS
BIO_SUPPORT_MATRIX
BIO_OPERATION_INVARIANT_MATRIX
BIO_PARITY_MATRIX
```

Do not maintain separate hand-written BioStructure operation lists unless they
are generated from the registry.

---

# BioOpDomain

The initial domain list is provisional. New domains should be added only when
they change operation contract behavior, testing profiles, or capability
boundaries; do not add a domain only to mirror a user-facing module name.

```rust
pub enum BioOpDomain {
    Selection,
    Hierarchy,
    Coordinate,
    Assembly,
    Annotation,
    Bonding,
    Polymer,
    ChemistryBridge,
}
```

---

# BioOpKind

```rust
pub enum BioOpKind {
    Weak,
    Strong,
}
```

## Weak

Does not fundamentally change row identity.

Examples:

- coordinate transforms
- cache recomputation
- annotations
- sequence assignment

## Strong

Changes row topology or hierarchy identity.

Examples:

- remove residues
- remove chains
- choose altloc
- assembly expansion
- merge structures

Strong ops require mappings.

---

# BioEditKind

```rust
pub enum BioEditKind {
    None,
    Local,
    Compacting,
    Expanding,
    Renumbering,
    Splitting,
    Merging,
    Transforming,
}
```

---

# BioStateSet

```rust
pub struct BioStateSet(u32);

impl BioStateSet {
    pub const HIERARCHY: Self;

    pub const RESIDUE_SPANS: Self;
    pub const CHAIN_SPANS: Self;
    pub const MODEL_SPANS: Self;

    pub const COORDINATE_ALIGNMENT: Self;

    pub const ENTITY_MAPPING: Self;

    pub const ALTLOC_GROUPS: Self;
    pub const ASSEMBLY_REFERENCES: Self;

    pub const BOND_REFERENCES: Self;

    pub const SELECTION_PROVENANCE: Self;

    pub const POLYMER_ANNOTATION: Self;
    pub const SECONDARY_STRUCTURE: Self;
}
```

---

# BioDerivedState

```rust
pub struct BioDerivedState(u64);

impl BioDerivedState {
    pub const ATOM_INDEX: Self;
    pub const RESIDUE_INDEX: Self;
    pub const CHAIN_INDEX: Self;
    pub const ENTITY_INDEX: Self;

    pub const SEQUENCE_CACHE: Self;
    pub const POLYMER_CACHE: Self;

    pub const ALTLOC_CACHE: Self;
    pub const ASSEMBLY_CACHE: Self;

    pub const BOND_CACHE: Self;

    pub const BACKBONE_GEOMETRY: Self;
    pub const SIDECHAIN_GEOMETRY: Self;

    pub const NUCLEIC_GEOMETRY: Self;

    pub const SECONDARY_STRUCTURE: Self;

    pub const CONTACT_MAP: Self;
    pub const GRAPH_CACHE: Self;
}
```

---

# Mapping System

```rust
pub struct BioStructureMapping {
    pub atoms: BioRowMapping<AtomId>,
    pub residues: BioRowMapping<ResidueId>,
    pub chains: BioRowMapping<ChainId>,
    pub entities: BioRowMapping<EntityId>,
    pub models: BioRowMapping<ModelId>,
    pub bonds: BioRowMapping<BondId>,
    pub assemblies: BioRowMapping<AssemblyId>,
}
```

```rust
pub struct BioSourceId(u32);

pub struct BioRowRef<Id> {
    pub source: BioSourceId,
    pub row: Id,
}

pub struct BioRowMapping<Id> {
    pub old_to_new: Vec<RowTargets<Id>>,
    pub new_to_old: Vec<RowSource<Id>>,
}

pub enum RowTargets<Id> {
    Deleted,
    One(Id),
    Many(Vec<Id>),
}

pub enum RowSource<Id> {
    Preserved(BioRowRef<Id>),
    Copied {
        source: BioRowRef<Id>,
        copy_index: u32,
    },
    Merged(Vec<BioRowRef<Id>>),
    Created,
}
```

Strong hierarchy-changing operations MUST produce mappings.

The simpler molecule-style shape:

```rust
old_to_new: Vec<Option<Id>>
new_to_old: Vec<Id>
```

is sufficient for deletion and simple renumbering, but it is not sufficient for
BioStructure operations. BioStructure mappings must be able to express:

- compacting deletion, such as `without_waters`
- renumbering and row reordering
- assembly expansion, where one source atom/residue/chain may produce many rows
- structure merge, where a new row may originate from a different source
  structure
- structure split, where one source structure produces multiple output
  structures
- altloc/disorder collapse, where multiple candidate rows may be merged into or
  represented by one selected row

For single-source operations, `BioSourceId(0)` represents the source structure.
For multi-source operations, mappings must record source structure identity in
addition to row identity.

This mapping shape is provisional, but source-indexed row references are not:
BioStructure mappings must be able to represent multi-source merge and
multi-output split provenance without changing public report semantics later.

Operations that cannot preserve a dependent block through the mapping must
explicitly drop it, recompute it, or return a structured unsupported/invalid
operation error.

## Mapping Requirement Profiles (Provisional)

The molecule-level `MappingRequirement::{None, Identity, Required}` shape is
the starting point, but BioStructure operations need profile-specific
obligations:

```text
none                  no mapping expected
identity              mapping may be identity and should be cheap to verify
required_single       changed single-output operation must report mappings
required_per_output   split operation must report one mapping per output
required_multi_source merge operation must report source-indexed mappings
required_expansion    expansion must report one-to-many source mappings
```

This may become a Bio-specific enum, or it may remain a profile-level rule
attached to `MappingRequirement::Required`. The implementation should not expose
an operation that relies on split, merge, or expansion mappings until this is
settled in code and tests.

---

# BioOpParts

```rust
pub struct BioOpParts<'a> {
    spec: &'static BioStructureOpSpec,

    working: BioStructure,

    mapping: Option<BioStructureMapping>,

    _source: PhantomData<&'a BioStructure>,

    #[cfg(feature = "bio-op-contracts")]
    trace: BioOperationTrace,
}
```

---

# Access API

```rust
impl BioOpParts<'_> {
    pub fn atoms(&self) -> &AtomBlock;
    pub fn atoms_mut(&mut self) -> &mut AtomBlock;

    pub fn residues_mut(&mut self) -> &mut ResidueBlock;
    pub fn chains_mut(&mut self) -> &mut ChainBlock;
    pub fn models_mut(&mut self) -> &mut ModelBlock;

    pub fn coordinates_mut(&mut self) -> &mut CoordinateBlock;
    pub fn bonds_mut(&mut self) -> &mut BondBlock;

    pub fn assemblies_mut(&mut self) -> &mut AssemblyBlock;
    pub fn annotations_mut(&mut self) -> &mut AnnotationBlock;

    pub fn derived_cache_mut(&mut self) -> &mut BioDerivedCacheBlock;

    pub fn remove_atoms(
        &mut self,
        atoms: &[AtomId],
    ) -> Result<&BioStructureMapping>;

    pub fn remove_residues(
        &mut self,
        residues: &[ResidueId],
    ) -> Result<&BioStructureMapping>;

    pub fn remove_chains(
        &mut self,
        chains: &[ChainId],
    ) -> Result<&BioStructureMapping>;

    pub fn remap_required_blocks(
        &mut self,
        mapping: &BioStructureMapping,
    ) -> Result<()>;

    pub fn clear_cache(
        &mut self,
        states: BioDerivedState,
    );

    pub fn recompute_hierarchy(&mut self) -> Result<()>;

    pub fn recompute_polymer_annotation(
        &mut self,
    ) -> Result<()>;

    pub fn recompute_altloc_groups(
        &mut self,
    ) -> Result<()>;

    pub fn finish(
        self,
        outcome: BioOpOutcome,
    ) -> Result<BioStructure>;

    pub fn finish_many(
        self,
        outcome: BioOpOutcome,
    ) -> Result<Vec<BioStructure>>;
}
```

Most BioStructure operations are single-output value-style transforms and use
`finish()`. Splitting operations such as `split_by_chain()` and
`split_by_model()` are multi-output transforms and must use an explicit
multi-output finish path. Multi-output operations must define source-structure
provenance, per-output mappings, and property/annotation propagation policy.

---

# Operation Body Style

```rust
#[bio_op_body(remove_waters, parts)]
fn remove_waters_impl() -> Result<BioOpOutcome, BioOperationError> {
    let waters = parts.find_residues_by_kind(
        ResidueKind::Water
    );

    parts.remove_residues(&waters)?;

    parts.recompute_hierarchy()?;

    parts.clear_cache(
        BioDerivedState::GRAPH_CACHE |
        BioDerivedState::CONTACT_MAP
    );

    Ok(BioOpOutcome::Changed)
}
```

---

# Registry Example

```rust
bio_structure_ops! {
    op remove_waters() {
        method: without_waters,
        impl_fn: remove_waters_impl,

        domain: selection,
        kind: strong,
        edit_kind: compacting,

        may_mutate: [
            atoms,
            residues,
            chains,
            models,
            coordinates,
            bonds,
            derived_cache,
        ],

        auto_remap: [
            coordinates,
            bonds,
            assemblies,
            annotations,
        ],

        must_handle: [
            hierarchy,
            residue_spans,
            chain_spans,
            model_spans,
            bond_references,
        ],

        needs_update: [
            atom_index,
            residue_index,
            chain_index,
            sequence_cache,
            graph_cache,
            contact_map,
        ],

        requires_mapping: required,

        report: [
            atom_mapping,
            residue_mapping,
            chain_mapping,
        ],

        allows_noop: true,

        feature: BIO_SELECTION_FEATURE,

        parity: gemmi_when_applicable,

        io_roundtrip: true,

        invariant_profile:
            "compacting_hierarchy_with_coordinates",

        parity_profile:
            "remove_waters_gemmi_like",
    }
}
```

The registry macro must generate all public operation metadata and test-driver
matrices. Tests should discover BioStructure operation obligations from:

```text
BIO_STRUCTURE_OPS
BIO_SUPPORT_MATRIX
BIO_OPERATION_INVARIANT_MATRIX
BIO_PARITY_MATRIX
```

`BIO_OPERATION_INVARIANT_MATRIX` includes unsupported and experimental
operations. Unsupported operations should still be tested for structured
unsupported errors, source immutability, and absence of partially mutated
placeholder output.

`BIO_PARITY_MATRIX` includes operations whose `BioParityPolicy` requires a
source-level comparison profile. A `GemmiWhenApplicable` operation may appear
before it is fully supported; until support is claimed, tests should assert
explicit unsupported behavior rather than silently skipping the operation.

---

# Operation Categories

This operation list is provisional. Public operations must still be registered
before implementation, and the registry entry is the source of truth if this
sketch and generated metadata diverge.

## Selection / Filtering

```rust
without_waters()
without_hydrogens()

keep_model(model)
keep_chain(chain)

keep_entity(entity)

keep_polymer()
keep_protein()
keep_nucleic_acid()

keep_ligands()

extract_ligands()

extract_chain(chain)
extract_model(model)
```

---

# AltLoc / Disorder

```rust
choose_altloc(policy)

keep_altloc(label)

collapse_disorder(policy)
```

## AltLocPolicy

```rust
pub enum AltLocPolicy {
    First,
    HighestOccupancy,
    Prefer(char),
    KeepAll,
    ErrorIfDisordered,
}
```

## AltLoc / Disorder Contract (Provisional)

Altloc/disorder handling must be typed state, not a loose atom property.

Provisional grouping rule:

```text
AltLocGroup = model + chain + residue identity + atom name group
blank altloc participates in every compatible conformer
named altlocs represent alternative conformers
occupancy is part of selection policy input
```

`choose_altloc(policy)` and `collapse_disorder(policy)` are strong operations
when they remove atom rows or residue rows. They must report mappings for atom,
residue, bond, coordinate, annotation, and assembly-dependent state as required
by their registry profile.

Tie-breaking is provisional and must be documented per policy:

```text
First              source order after parsing
HighestOccupancy   highest occupancy, then source order
Prefer(label)      requested label if present, otherwise explicit error unless
                   the policy says fallback is allowed
KeepAll            no compaction; may only normalize/cache altloc groups
ErrorIfDisordered  returns structured invalid-input error if disorder exists
```

Gemmi-compatible parsing tests should preserve raw altloc labels, blank-altloc
behavior, occupancy, and source order so later policy decisions can reproduce
or intentionally diverge from Gemmi.

---

# Hierarchy Editing

```rust
rename_chain(old, new)

renumber_residues(policy)

normalize_residue_ids()

split_by_chain()
split_by_model()

merge_structures(other)
```

---

# Coordinate Operations

```rust
center()

translate(vec3)

rotate(rotation)

transform(matrix4)

superpose_onto(
    reference,
    atom_selection,
)
```

Most coordinate ops are weak operations.

---

# Assembly Operations

```rust
list_assemblies()

apply_assembly(id)

expand_symmetry()
```

Assembly expansion is a strong expanding op.

## Assembly Expansion Contract (Provisional)

Assembly expansion is not ordinary compacting topology migration. It is a
strong expanding operation that copies hierarchy rows, transforms coordinates,
and records assembly provenance.

`apply_assembly(id)` / `expand_symmetry()` must define:

```text
source assembly id or pseudo-assembly id
operator ids applied to each copied row
chain/subchain selection used by each generator
chain naming policy for copied chains
whether original unexpanded rows are preserved or replaced
whether assembly metadata is consumed, preserved, or rewritten
how bonds/connections across copied rows are created or preserved
```

Required mapping behavior:

```text
old atom/residue/chain rows may map to many new rows
new rows record source row and applied operator
coordinate rows are transformed by the same operator provenance
annotations are preserved, remapped, dropped, or explicitly marked stale
```

Gemmi's assembly utilities include chain naming and copy/merge behavior. COSMolKit
may reproduce a pinned Gemmi profile, but the chosen chain naming policy must be
explicit in the operation profile and test output schema.

---

# Annotation Operations

```rust
assign_residue_kinds()

assign_polymer_kinds()

assign_sequences()

assign_secondary_structure()

assign_backbone_atoms()

assign_nucleic_acid_atoms()
```

---

# Bond Operations

```rust
infer_standard_bonds()

read_conect_bonds()

infer_disulfide_bonds()

infer_covalent_ligand_bonds()

clear_bonds()
```

## BondSource

```rust
pub enum BondSource {
    ExplicitConect,
    StandardResidueTemplate,
    DistanceInferred,
    UserProvided,
    Unknown,
}
```

BioStructure bonds should remain provenance-aware.

BioStructure bonds are structural connectivity records. They are not a promise
that full small-molecule chemistry perception has been performed.

Provisional boundary:

```text
BioStructure bond:
    atom endpoints
    source/provenance
    optional structural connection type
    optional order when explicitly provided or template-derived

Molecule bond:
    chemistry graph edge
    valence/aromaticity/stereo/sanitize semantics
```

BioStructure operations must not infer valence, aromaticity, or sanitized
small-molecule chemistry unless they explicitly enter the chemistry bridge and
return a `Molecule`. Protein CONECT records, mmCIF `struct_conn`, residue
template bonds, inferred disulfides, and ligand covalent links should remain
distinguishable by `BondSource`.

---

# Chemistry Bridge Operations

```rust
extract_ligand_as_molecule()

extract_residue_as_molecule()

extract_chain_as_polymer_molecule()

attach_ligand(molecule, pose)
```

Once converted into `Molecule`,
the chemistry contract system becomes active.

---

# Protein-Specific Extensions

```rust
compute_phi_psi()

compute_chi_angles()

select_backbone()

select_sidechains()

mutate_residue()

standardize_residue()
```

These should live in:

```rust
ProteinView
ProteinOps
```

rather than polluting BioStructure core.

---

# DNA/RNA Extension Hooks

```rust
compute_backbone_torsions()

compute_sugar_pucker()

detect_base_pairs()

detect_base_stacking()

select_bases()

select_sugar_phosphate_backbone()
```

These are future extension points.

---

# Parser Import Contract (Provisional)

PDB/mmCIF parsing constructs a new `BioStructure`; it should not be modeled as
a normal value-style mutation operation on an existing structure.

Parser entry points should instead have an explicit import contract:

```text
input format
parser options
source reference profile
structured unsupported behavior
post-parse invariant profile
Gemmi comparison profile when applicable
raw metadata preservation policy
```

The parser must return structured parse/import errors. It must not emit a
plausible-looking partial structure for unsupported required fields unless the
selected parser mode explicitly allows preserved-only behavior and records it.

Parser import may populate typed raw/provenance fields that are not yet used by
operations. Such fields should be marked preserved-only, not interpreted.

---

# Annotation and Property Remap Policy (Provisional)

BioStructure annotations are not generic molecule properties. Row-attached
annotations must declare remap behavior just like coordinates and bonds.

Examples:

```text
secondary structure ranges
SIFTS / UniProt cross references
TLS groups
entity_poly_seq records
source organism metadata
crystal cell and spacegroup
raw remarks / audit metadata
assembly annotations
validation annotations
```

Single-source filtering operations should preserve structure-level metadata by
default, but row-attached annotations must be remapped, recomputed, dropped, or
reported unsupported. Split and merge operations must define how metadata is
copied, namespaced, merged, or rejected on conflict.

Provisional rule:

```text
annotations tied to AtomId/ResidueId/ChainId/ModelId must not remain attached
to stale ids after any strong operation
```

Do not use `annotations` or `props` as an escape hatch for typed PDB/mmCIF state
that affects selection, validation, writeback, or Gemmi comparison.

---

# Invariant Profiles

Profile names are provisional until `bio_structure_ops!` exists. The important
contract is that profiles are generated from registry metadata and remain the
source of truth for tests.

```text
readonly

weak_annotation_update

coordinate_transform

compacting_hierarchy_with_coordinates

altloc_compaction

assembly_expansion

structure_merge

structure_split

bond_inference

chemistry_bridge

parser_import

source_identifier_preservation

gemmi_parse_reproduction
```

---

# Core Invariant Checks

## Hierarchy

```text
atom ids continuous
residue ids continuous
chain ids continuous
model ids continuous

atom.residue_id valid
residue.chain_id valid
chain.entity_id valid
```

## Spans

```text
residue atom spans valid
chain residue spans valid
model spans valid
```

## Coordinates

```text
coordinate rows aligned with atoms
```

## AltLoc

```text
altloc groups valid
```

## Assembly

```text
assembly references valid
assembly generator chains/subchains valid
assembly operators valid
```

## Bonds

```text
bond atom ids valid
bond provenance valid
```

## Source Identifiers

```text
source identifiers preserved or intentionally normalized
selection addresses resolve to valid rows when applicable
auth/label/subchain/entity fields remain internally consistent
```

## Mapping

```text
required mappings are present
mapping source and target ids are valid
one-to-many expansion mappings are valid
multi-source merge mappings identify their source structure
split outputs carry per-output mappings and provenance
```

## Source

```text
source object unchanged after operation
```

---

# Parity Policy

```rust
pub enum BioParityPolicy {
    NotApplicable,

    GemmiWhenApplicable,

    BiopythonWhenApplicable,

    PdbSpecRequired,

    RequiredNow,
}
```

BioStructure parity compares externally observable structural behavior, not
internal implementation shape. It must remain separate from invariant tests.

Gemmi is the pinned source-level reproduction target for future PDB/mmCIF
parsing. Gemmi comparison profiles should record:

```text
pinned Gemmi commit
input format and parser options
case id and input reference
structured output schema
normalization rules, if any
known-failure expectation, if any
```

Recommended Gemmi structured dump fields:

```text
models and model numbers
chains and chain/source ids
subchains / label_asym_id
entities and entity types
residues with auth id, label id, insertion code, name, segment
atoms with name, element, altloc, serial, occupancy, b_iso, charge
coordinates aligned to atom rows
connections / CONECT / struct_conn with provenance
assemblies, generators, operators, selected chains/subchains
secondary-structure annotations when parsed
raw source metadata that affects writeback or later parsing behavior
```

Do not compare vague dumps or incidental debug formatting. Every Gemmi parity
profile must define exactly which fields are compared and which source-level
differences are intentionally normalized.

`PdbSpecRequired` is for behavior where the wwPDB/PDBx specification is the
oracle even if Gemmi and Biopython differ. Such cases must cite the relevant
format rule in the test or profile documentation.

## Gemmi Reproduction Contract (Provisional)

`GemmiWhenApplicable` means Gemmi is the pinned reproduction target for the
declared parser or operation profile. It does not mean Gemmi is the permanent
semantic truth for every BioStructure feature.

Gemmi may use heuristics for:

```text
PDB element inference
subchain assignment
entity/type inference
chain naming during assembly expansion
connection completion
altloc/conformer iteration helpers
format-specific tolerated malformed input
```

When COSMolKit follows Gemmi, the test profile must identify the compared
behavior and the pinned Gemmi commit. When COSMolKit intentionally follows the
PDB/mmCIF specification or a COSMolKit-native policy instead, the operation or
parser profile must use `PdbSpecRequired` or `NotApplicable` and document the
difference.

Allowed responses to Gemmi mismatch:

```text
fix COSMolKit reproduction
update the profile with source evidence
mark behavior unsupported
add executable known-failure
document intentional divergence
```

Disallowed responses:

```text
skip mismatching fields without a profile rule
compare only debug text
guess a heuristic and call it Gemmi-compatible
silently normalize source identifiers
```

---

# Final Principle

```text
BioStructure prioritizes:

- hierarchy correctness
- coordinate correctness
- assembly semantics
- structural provenance
- polymer annotations

rather than full chemistry correctness.
```

```text
MoleculeOps:
    atom/bond/valence/aromaticity/stereo

BioStructureOps:
    atom/residue/chain/model/entity/
    assembly/altloc/coordinates/polymer
```
