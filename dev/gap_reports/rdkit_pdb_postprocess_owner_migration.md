# RDKit PDB Postprocess Owner Migration

## Outcome

`cosmolkit-io` is now the single production behavior owner for the
chemistry-owned tail of RDKit PDB molecule construction:

- `ConnectTheDots()` / `ConnectTheDots_Large()` hash-grid proximity bonding;
- `IsBonded()` distance and covalent-radius rules;
- cross-residue metal, halogen, noble-gas, and water blacklisting;
- multivalent-hydrogen cleanup;
- `StandardPDBResidueBondOrders()` with the complete source table;
- `BasicPDBCleanup()` four-coordinate neutral-nitrogen correction; and
- `StandardPDBResidueChirality()` filtering after 3D stereo assignment.

The owner API is implemented in
`crates/cosmolkit-io/src/pdb_chemistry.rs` as
`postprocess_pdb_detached()` and
`apply_standard_pdb_residue_chirality_detached()`.

## Runtime Boundary

The former `crates/cosmolkit/src/runtime/io/pdb_molecule.rs` copy has been
deleted with the old mixed runtime. The canonical `cosmolkit` public adapter
still needs to perform the BioStructure-to-detached-row projection, call the
detached owner, apply operation policy, validate the result, and construct the
live `Molecule`. No PDB chemistry implementation may be added to `cosmolkit`.

The required public flow is:

```text
BioStructure rows
-> detached atom/bond/conformer assembly
-> cosmolkit_io::postprocess_pdb_detached
-> runtime sanitize/remove-H policy
-> 3D stereo assignment
-> cosmolkit_io::apply_standard_pdb_residue_chirality_detached
```

This ordering keeps live `Molecule`, cache, and operation policy outside
`cosmolkit-io` while leaving all PDB chemistry in the detached owner.

## Source-Reproduction Status

The moved functions retain the verbatim RDKit source anchors and two-axis
markers from:

- `Code/GraphMol/FileParsers/ProximityBonds.cpp`;
- `Code/GraphMol/FileParsers/PDBParser.cpp`; and
- the pinned periodic-table covalent-radius data used by RDKit.

`BasicPDBCleanup()` now delegates explicit-valence evaluation to
`cosmolkit_core::calculate_explicit_valence_for_topology()` instead of using
the former local integer bond-order approximation. This closes the source call
edge without creating a second valence implementation.

Three conservative `RDKit❗✔️` markers remain:

1. the inline source excerpt shows the first periodic-table rows while the
   complete pinned `rCov` column is stored in the adjacent Rust table;
2. the `ConnectTheDots_Large()` coordinator points to the separately anchored
   cleanup helper body; and
3. the residue-table source excerpt points to the complete adjacent Rust tuple
   table.

These are source-framing disclosures, not silent behavior fallbacks. No
first-axis unsupported marker remains in the detached PDB chemistry owner.

## Parity Evidence

Focused expected values were checked against the repository RDKit environment
with `sanitize=false` and `removeHs=false`:

| Case | RDKit result | Detached gate |
|---|---|---|
| C/O at 1.2 Å | one single bond, stored O-to-C | exact endpoints/order |
| H/H at 0.7 Å | no bond | exact count |
| ALA `C`/`O`, flavor 8, explicit CONECT | one double bond | exact order |
| four-coordinate neutral N | formal charge `+1` | exact charge |
| H between C at 1.0 Å and O at 1.1 Å | retains only H-C | exact retained endpoints |

The standard-residue chirality unit gate additionally checks that an ALA alpha
carbon retains its assigned tetrahedral/CIP state, a non-allowlisted ALA atom
is cleared, and a hetero-residue atom is not filtered.

Historical validation recorded before the clean-break removal:

```text
cargo test -p cosmolkit-io --release
44 passed, 0 failed

cargo test -p cosmolkit --features op-contracts-strict pdb_molecule
12 passed, 0 failed

cargo test -p cosmolkit --release --features op-contracts-strict pdb_molecule
12 passed, 0 failed
```

The `cosmolkit` rows above are no longer current gates because the copied
runtime adapter was deleted. New validation must exercise the canonical public
adapter after it is implemented.

## Remaining PDB / IO Boundaries

- The detached PDB writer still returns structured unsupported for aromatic
  input requiring RDKit Kekulize preprocessing; Kekulize is a topology
  operation and must not be approximated inside the writer.
- The BioStructure projection and direct fixed-column detached reader are
  still separate input adapters. Their chemistry tail is shared, but parser
  ownership cleanup remains part of the later runtime/IO consolidation.
- Concrete `TopologyBlock` cannot represent all V2000/V3000 query, collection,
  OBJ3D, and SGroup source states; those branches continue to fail closed until
  the corresponding model/lowering contracts exist.
