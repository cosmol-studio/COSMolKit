# RDKit Morgan/MACCS Source Inventory Gap Report

## Scope

This report inventories the pinned RDKit source under
`third_party/rdkit/Code/GraphMol/Fingerprints/` for the Morgan and MACCS
surfaces targeted by `dev/rdkit_morgan_maccs_full_port_plan.md`, and compares
that surface to the currently exposed COSMolKit fingerprint API.

## RDKit Source Inventory

### `FingerprintGenerator.h` / `FingerprintGenerator.cpp`

Current RDKit responsibilities:

- `AdditionalOutput` allocation holders and reinitialization semantics
- `FingerprintArguments` construction, `commonArgumentsString`, `toJSON`,
  `fromJSON`
- `FingerprintGenerator` construction, ownership, `infoString`, `toJSON`,
  `fromJSON`
- `getFingerprintHelper` orchestration
- count-simulation helper bits and temporary additional-output duplication
- `getSparseCountFingerprint`, `getSparseFingerprint`, `getCountFingerprint`,
  `getFingerprint`

Gap status in COSMolKit:

- Only a partial Morgan-specific inline path exists.
- AdditionalOutput is represented, but RDKit allocation/reinitialization
  semantics are not yet source-ported as a generator-owned object model.
- The shared generator/argument pipeline is not yet fully ported.

### `FingerprintUtil.h` / `FingerprintUtil.cpp`

Current RDKit responsibilities:

- `MorganFingerprints::ss_matcher`
- `MorganFingerprints::defaultFeatureSmarts`
- `MorganFingerprints::getConnectivityInvariants`
- `MorganFingerprints::getFeatureInvariants`
- `RDKitFPUtils::buildDefaultRDKitFingerprintAtomInvariants`
- path enumeration and bond-hash helpers used by other fingerprint families

Gap status in COSMolKit:

- Connectivity invariants exist in partial form.
- Feature invariants are currently approximated locally instead of being
  source-backed through SMARTS/substructure matching.
- `ss_matcher` and `defaultFeatureSmarts` are not source-ported.

### `MorganGenerator.h` / `MorganGenerator.cpp`

Current RDKit responsibilities:

- `MorganAtomInvGenerator`
- `MorganFeatureAtomInvGenerator`
- `MorganBondInvGenerator`
- `MorganArguments`
- `MorganAtomEnv`
- `MorganEnvGenerator`
- `getMorganGenerator`
- `getAtomInvariants`, `getBondInvariants`, `getEnvironments`
- info-string, JSON, clone, result-size, and additional-output semantics

Gap status in COSMolKit:

- Morgan fingerprint behavior exists only as a partial inline implementation in
  `crates/cosmolkit-core/src/properties/fingerprint.rs`.
- Generator ownership, JSON wiring, and exact environment generation are not
  source-complete.
- Feature and bond invariant generators are not yet ported as RDKit-shaped
  generator objects.

### `MorganFingerprints.h` / `MorganFingerprints.cpp`

Current RDKit responsibilities:

- `getFingerprint`
- `getHashedFingerprint`
- `getFingerprintAsBitVect`
- sparse/count/bit-vector wrapper behavior over `FingerprintGenerator`

Gap status in COSMolKit:

- Public Morgan wrappers exist, but they still call a simplified inline path.
- Wrapper behavior around exact generator plumbing and additional-output
  propagation is not yet source-complete.

### `MACCS.h` / `MACCS.cpp`

Current RDKit responsibilities:

- `Patterns` table initialization
- `GenerateFP`
- `MACCSFingerprints::getFingerprintAsBitVect`
- 167-bit internal vector behavior with bit 0 unused and public 166-bit
  projection

Gap status in COSMolKit:

- MACCS is present only as a heuristic/local rule set.
- RDKit SMARTS patterns, bit numbering, and exact key logic are not
  source-ported.
- The current implementation is not acceptable as parity closure for the
  pinned RDKit oracle.

## Current COSMolKit Public Fingerprint Surface

Exposed from `crates/cosmolkit-core/src/lib.rs`:

- `Fingerprint`
- `FingerprintError`
- `MorganAdditionalOutput`
- `MorganAtomInvariantsGenerator`
- `MorganBondInvariantsGenerator`
- `MorganFingerprintOutput`
- `MorganFingerprintParams`

Exposed from `crates/cosmolkit-core/src/properties/fingerprint.rs`:

- `morgan_fingerprint`
- `morgan_fingerprint_with_output`
- `topological_fingerprint`
- `maccs_fingerprint`
- `TopologicalFingerprintParams`
- `MaccsFingerprintParams`

Policy surface:

- `crates/cosmolkit-core/src/support.rs` currently marks the Morgan fingerprint
  feature as experimental and explicitly unfinished.
- The public docs and tests still describe Morgan parity as pending strict
  bit-identical coverage.

## Gap Summary

1. The RDKit generator pipeline is not yet source-complete in COSMolKit.
2. Morgan feature invariants still rely on local approximation instead of
   source-backed SMARTS/substructure machinery.
3. MACCS is still heuristic and not source-ported.
4. The public surface currently exposes unfinished fingerprint behavior that
   must be replaced or narrowed as the RDKit port progresses.

## Immediate Follow-Ups Implied by the Plan

- Port `AdditionalOutput` allocation and storage semantics.
- Port `FingerprintArguments` and Morgan argument/generator objects.
- Replace feature-invariant heuristics with source-backed SMARTS/substructure
  behavior.
- Port MACCS patterns and exact key logic.
