# RDKit-Compatible PDB Molecule Conversion Plan

This plan defines how COSMolKit should implement a `Molecule` result that is
compatible with RDKit `Chem.MolFromPDBBlock()` while preserving the current
BioStructure-first architecture.

## Architecture

PDB/mmCIF structural parsing remains owned by the Gemmi-aligned `BioStructure`
path. Do not add a second public PDB parser.

The target flow is:

```text
PDB text
-> BioStructure
-> RDKit-compatible BioStructure-to-Molecule conversion profile
-> Molecule
```

The user-facing convenience API may be:

```python
mol = Molecule.from_pdb_block(
    pdb,
    sanitize=True,
    remove_hs=True,
    proximity_bonding=True,
)
```

Internally this must reuse the BioStructure parse, then apply RDKit-compatible
molecule conversion and post-processing.

`Protein` is not the correct parity boundary for `Chem.MolFromPDBBlock()`.
`Protein` is a protein-only projection and may drop ligands, waters, ions,
nucleic acids, HETATM records, and connection context. A future
`Protein.to_molecule()` may exist, but it must be documented as a protein-subset
projection and must not claim RDKit PDB parser equivalence.

## Public Rust API Shape

Add an explicit conversion profile:

```rust
pub struct RdkitPdbMolProfile {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub flavor: u32,
    pub proximity_bonding: bool,
}
```

Add BioStructure conversion:

```rust
impl BioStructure {
    pub fn to_rdkit_pdb_molecule(
        &self,
        profile: RdkitPdbMolProfile,
    ) -> Result<Molecule, PdbMoleculeConversionError>;
}
```

Add Molecule convenience construction:

```rust
impl Molecule {
    pub fn from_pdb_block_with_options(
        text: &str,
        profile: RdkitPdbMolProfile,
    ) -> Result<Self, PdbMoleculeConversionError>;
}
```

The implementation should live in a narrow module such as:

```text
crates/cosmolkit-core/src/io/pdb_molecule.rs
```

## RDKit Source Anchors

Follow `dev/source_reproduction_protocol.md`. Relevant RDKit source blocks must
be copied verbatim into corresponding Rust functions with two-axis markers.

Primary source anchors:

```text
RDKit Code/GraphMol/FileParsers/PDBParser.cpp::PDBAtomLine
RDKit Code/GraphMol/FileParsers/PDBParser.cpp::PDBBondLine
RDKit Code/GraphMol/FileParsers/PDBParser.cpp::parsePdbBlock
RDKit Code/GraphMol/FileParsers/PDBParser.cpp::BasicPDBCleanup
RDKit Code/GraphMol/FileParsers/PDBParser.cpp::StandardPDBResidueChirality
RDKit Code/GraphMol/FileParsers/ProximityBonds.cpp::ConnectTheDots
RDKit Code/GraphMol/FileParsers/ProximityBonds.cpp::IsBonded
RDKit Code/GraphMol/FileParsers/ProximityBonds.cpp::IsBlacklistedPair
RDKit Code/GraphMol/FileParsers/ProximityBonds.cpp::StandardPDBResidueBondOrders
RDKit Code/GraphMol/MonomerInfo.h::AtomPDBResidueInfo
```

## RDKit-Compatible Retention And Filtering

The conversion profile should match `Chem.MolFromPDBBlock()` for modeled
Molecule state.

Retain:

- `ATOM  ` records as atoms.
- `HETATM` records as atoms.
- Protein, nucleic acid, ligand, water, ion, and other PDB atoms unless RDKit
  would filter them.

Default filtering when `(flavor & 1) == 0`:

- Ignore alternate locations where altLoc is not space, `A`, or `1`.
- Ignore XPLOR pseudo atoms with coordinate fields
  `9999.0009999.0009999.000`.
- Ignore NMR pseudo atoms where PDB atom-name columns match RDKit's
  `ptr[12] == ' '` and `ptr[13] == 'Q'` condition.
- Ignore PDB dummy residues with residue name `DUM`.

When `(flavor & 1) != 0`, do not apply the above filtering, matching RDKit.

## Atom State

Element inference must follow RDKit order:

1. Atomic symbol columns 76-77.
2. PDB atom name fallback.

The conversion must preserve:

- Atomic number, including `D` and `T` as hydrogen isotopes 2 and 3.
- Formal charge from PDB charge columns 78-79 using RDKit rules.
- One 3D conformer from PDB coordinates.
- Conformer `is_3d` set if any z coordinate is non-zero.
- Atom-level PDB residue info.

`AtomPdbResidueInfo` fields should map RDKit `AtomPDBResidueInfo`:

```text
atom name
serial number
altLoc
residue name
residue number
chain id
insertion code
occupancy
temperature factor
isHeteroAtom
```

`HETATM` must set `isHeteroAtom=true`; `ATOM  ` must set it to false.

## Bond State

Explicit `CONECT` handling must follow RDKit `PDBBondLine()`:

- Ignore source or destination serials that do not map to retained atoms.
- Skip self bonds.
- Add explicit single bonds by default.
- Use zero bonds for RDKit blacklisted explicit pairs.
- Use repeated CONECT entries to upgrade single bonds to double, triple, or
  quadruple according to RDKit's byte bitmap logic.

`proximity_bonding=True` must run RDKit-compatible proximity bonding:

- Use RDKit distance thresholds: `EXTDIST`, `MAXDIST`, `MINDIST2`,
  `MAXDIST2`.
- Use covalent radii and the same `IsBonded` rule.
- Use RDKit's blacklist behavior for cross-residue metals/noble gases/halogens
  and water.
- Implement hydrogen cleanup for multivalent hydrogens according to
  `ConnectTheDots`.

After proximity bonding and explicit CONECT processing, apply:

- `StandardPDBResidueBondOrders()`.
- `BasicPDBCleanup()`.
- `StandardPDBResidueChirality()` when the corresponding stereochemistry state
  is modeled.

## Suggested Function-Level Breakdown

```text
bio_structure_to_rdkit_pdb_molecule(...)
include_atom_like_rdkit(...)
atom_spec_from_bio_atom_like_rdkit(...)
pdb_residue_info_from_bio_rows_like_rdkit(...)
add_conformer_from_bio_coordinates(...)
apply_conect_records_like_rdkit(...)
pdb_bond_line_from_bio_conect_like_rdkit(...)
connect_the_dots_like_rdkit(...)
is_bonded_like_rdkit(...)
same_pdb_residue_like_rdkit(...)
is_blacklisted_pair_like_rdkit(...)
standard_pdb_residue_bond_orders_like_rdkit(...)
standard_pdb_double_bond_like_rdkit(...)
basic_pdb_cleanup_like_rdkit(...)
standard_pdb_residue_chirality_like_rdkit(...)
finish_rdkit_pdb_molecule_processing(...)
```

## Python API Shape

Expose Python APIs matching the Rust behavior:

```python
Molecule.from_pdb_block(
    text: str,
    *,
    sanitize: bool = True,
    remove_hs: bool = True,
    flavor: int = 0,
    proximity_bonding: bool = True,
) -> Molecule
```

Optional BioStructure conversion:

```python
BioStructure.to_molecule(
    profile: str = "rdkit_pdb",
    *,
    sanitize: bool = True,
    remove_hs: bool = True,
    flavor: int = 0,
    proximity_bonding: bool = True,
) -> Molecule
```

Do not expose legacy compatibility aliases.

## Test Plan

Use RDKit-generated fixtures where RDKit is available. Compare modeled state:

```text
atom count
bond count
atomic numbers
isotopes
formal charges
bond begin/end/order
3D coordinates
conformer is_3d
PDB residue info fields
SMILES when sanitization succeeds
```

Required focused tests:

- `ATOM  ` and `HETATM` are both retained.
- Protein, ligand, water, ion, and nucleic-acid atoms are retained unless RDKit
  filters them.
- Default altLoc filtering matches RDKit.
- `flavor & 1` disables alternate/dummy/pseudo filtering.
- `D` and `T` become hydrogen isotopes.
- PDB charge columns map to formal charges.
- `HETATM` sets `is_hetero_atom`.
- Coordinates create a conformer and set `is_3d` correctly.
- `CONECT` creates explicit bonds.
- Repeated `CONECT` records upgrade bond order.
- Blacklisted explicit pairs produce RDKit-compatible zero bonds.
- `proximity_bonding=True` and `False` differ as RDKit does.
- Water cross-residue proximity blacklisting matches RDKit.
- Multivalent hydrogen cleanup matches RDKit.
- Standard residue double-bond assignment matches RDKit.
- Four-valent neutral nitrogen cleanup matches RDKit.
- `sanitize` and `remove_hs` parameter order matches RDKit.

## Completion Criteria

`Molecule.from_pdb_block(pdb)` may be documented as RDKit-compatible only for
the modeled state once the tests above pass against RDKit fixtures and all
source-backed functions carry source-reproduction markers.
