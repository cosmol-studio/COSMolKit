# RDKit CIPLabeler Gap Report For Morgan Chirality

## Trigger

Command:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_morgan_fingerprint_parity -- --nocapture
```

Observed failing test:

```text
morgan_fingerprint_matches_rdkit_golden_across_param_branches
```

Observed failing branch:

```text
row 83
SMILES: O=C(NC[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2)[C@@H]1C[C@H]2c3ccccc3[C@@H]1c1ccccc12
branch: r2_n2048_chiral_true_bonds_true
```

COSMolKit currently fails closed:

```text
unsupported Morgan fingerprint option includeChirality:
MorganEnvGenerator requires RDKit CIPLabeler parity before recomputing missing _CIPCode labels
```

This is the correct temporary behavior. Returning Morgan bits from a partial or
legacy CIP implementation would violate the source-reproduction protocol and
would produce chemically meaningful-looking fingerprint output from unsupported
state.

## RDKit Source Dependency

The Morgan chirality path depends on RDKit's new CIPLabeler path, not only on
legacy `Chirality.cpp` CIP rank/code helpers.

Relevant Morgan call sites:

- `third_party/rdkit/Code/GraphMol/Fingerprints/FingerprintGenerator.cpp`
  - `FingerprintGenerator<OutputType>::getFingerprintHelper`
  - When `df_includeChirality` is true and `_StereochemDone` is absent, RDKit
    clones the molecule and calls `MolOps::assignStereochemistry`.
- `third_party/rdkit/Code/GraphMol/Fingerprints/MorganGenerator.cpp`
  - `MorganEnvGenerator<OutputType>::getEnvironments`
  - When `df_includeChirality` is true, legacy stereo perception is disabled,
    and `_CIPComputed` is absent, RDKit calls
    `CIPLabeler::assignCIPLabels`.
  - Chiral atom hashing then consumes `_CIPCode` values `R` and `S`.
- `third_party/rdkit/Code/GraphMol/Fingerprints/MorganGenerator.cpp`
  - `MorganBondInvGenerator::getBondInvariants`
  - For chiral double bonds, RDKit also calls `CIPLabeler::assignCIPLabels`
    when `_CIPComputed` is absent, then consumes bond `_CIPCode` values `E`
    and `Z`.

Relevant CIPLabeler source surface:

- `third_party/rdkit/Code/GraphMol/CIPLabeler/CIPLabeler.cpp`
  - `findConfigs`
  - `labelAux`
  - `label`
  - `assignCIPLabels(ROMol &, const boost::dynamic_bitset<> &, const boost::dynamic_bitset<> &, unsigned int)`
  - `assignCIPLabels(ROMol &, unsigned int)`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/CIPMol.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/Digraph.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/Node.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/Edge.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/Sort.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/Priority.h`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/Descriptor.h`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/configs/Configuration.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/configs/Tetrahedral.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/configs/Sp2Bond.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/configs/AtropisomerBond.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/SequenceRule.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Rules.h`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Rule1a.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Rule1b.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Rule2.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Rule3.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Rule4a.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Rule4b.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Rule4c.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Rule5New.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Rule6.{h,cpp}`
- `third_party/rdkit/Code/GraphMol/CIPLabeler/rules/Pairlist.h`

The source surface is approximately 6k lines across graph construction,
sequence-rule comparison, auxiliary descriptor propagation, and configuration
labeling. It is not a narrow Morgan-local helper.

## Why Existing COSMolKit CIP Is Not Enough

COSMolKit already has legacy CIP-related code:

- `_CIPRank` assignment in `crates/cosmolkit-core/src/chemistry/stereo.rs`
- legacy-style `_CIPCode` assignment for some atom stereochemistry paths
- SMILES writer stereo cleanup and assignment paths
- drawing and graph-feature consumers of existing `_CIPCode`

Those paths are not equivalent to RDKit `CIPLabeler::assignCIPLabels`.
The existing code is marked partial in multiple source anchors and does not
model the CIPLabeler digraph/rule system, auxiliary descriptors, Rule4/Rule5/Rule6
behavior, `Sp2Bond` E/Z labeling, or `AtropisomerBond` labeling.

Using the existing legacy path to unblock Morgan would be a heuristic shortcut:

- it would satisfy only some tetrahedral cases;
- it would not reproduce RDKit's `_CIPComputed` semantics;
- it would not reproduce bond `_CIPCode` for E/Z double bonds;
- it would not reproduce pseudoasymmetric lower-case descriptors;
- it would hide unsupported behavior behind exact-looking fingerprint bits.

## Required Closure

Morgan include-chirality parity must remain blocked until COSMolKit has a
source-backed private RDKit CIPLabeler port and Morgan calls it at the same
source call sites where RDKit does.

The follow-on plan steps split that work by RDKit function group. A step is not
complete unless it includes copied source anchors, Rust implementation, tests
at the relevant boundary, and no placeholder or compatibility-only branch.
