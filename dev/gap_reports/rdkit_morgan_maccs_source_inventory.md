# RDKit Morgan/MACCS Source Inventory

## Scope

This report records the exact boundary of the Morgan and MACCS implementation
against pinned RDKit `2026.03.1`. It is not a completion claim for every class,
wrapper, overload, input-state preparation branch, or fingerprint family under
`GraphMol/Fingerprints`.

The source inventory is anchored to:

- `FingerprintGenerator.h` / `FingerprintGenerator.cpp`
- `FingerprintUtil.h` / `FingerprintUtil.cpp`
- `MorganGenerator.h` / `MorganGenerator.cpp`
- `MorganFingerprints.h` / `MorganFingerprints.cpp`
- `MACCS.h` / `MACCS.cpp`
- the dependent CIPLabeler, SMARTS/substructure, ring, valence, and periodic
  table paths reached by the exposed branch matrix

## Morgan: Implemented And Validated Boundary

The following exposed outputs have source-backed implementations and strict
RDKit golden comparisons:

- sparse count fingerprint
- sparse bit fingerprint
- hashed count fingerprint
- explicit bit fingerprint
- count simulation and count bounds
- connectivity and default feature atom invariants
- explicit Morgan bond invariants
- `fromAtoms`, RDKit's currently unused `ignoreAtoms` behavior, custom atom
  invariants, and custom bond invariants
- redundant-environment and `numBitsPerFeature` branches
- `AdditionalOutput` atom counts, atom-to-bits, bit-info map, and atoms-per-bit
- covered chirality/CIPLabeler branches in the committed Morgan matrix

This boundary is validated by
`crates/cosmolkit-core/tests/rdkit_morgan_fingerprint_parity.rs`. It does not
mean every RDKit generator object API or every molecule preparation state is
available. In particular:

- fingerprint generation with `includeChirality=true` and missing
  `_StereochemDone` fails closed until the exact RDKit
  `MolOps::assignStereochemistry` preparation branch is reproduced there
- supplied feature SMARTS generator patterns remain unsupported
- any branch marked open in the source file remains unsupported even when an
  adjacent branch shares the same output type

Those branches return `FingerprintError`; they do not use legacy CIP,
alternative SMARTS classification, or another local fingerprint as fallback.

## MACCS: Implemented And Validated Boundary

The exposed MACCS path follows `MACCS.cpp::GenerateFP`:

- source SMARTS patterns and raw RDKit bit numbering
- direct element/count/ring/fragment key logic
- raw 167-bit vector with bit 0 unused
- COSMolKit public 166-bit projection from raw bits 1 through 166

Strict tests compare both vectors on targeted per-key fixtures and the selected
SMILES profile. Non-166 public output sizes return `FingerprintError`. SMARTS
matcher construction is fallible and its error is propagated; it cannot panic
or silently omit a key.

## Other Fingerprint Families

`RDKFingerprintMol` and Avalon are not part of the completed Morgan/MACCS
boundary.

- `RDKFingerprintMol` depends on the RDKitFP generator path, branched subgraph
  enumeration, random bit generation, density folding, atom invariants,
  atom-bit output, and bit-path output.
- Avalon depends on the external Avalon/reaccs implementation, including
  `bitFlags`, `isQuery`, `resetVect`, and byte-rounded vector semantics.

The previous local DFS/path-hash implementations were not source-equivalent and
have been removed. The retained public entry points return structured
unsupported errors until those complete source paths are ported and exact-bit
tests pass.

## Acceptance Rule

Only exact equality on the enumerated source-backed boundary is accepted.
Similarity correlation, structurally similar hashes, partial key sets, and
99.9% bit agreement are failures, not compatibility states.
