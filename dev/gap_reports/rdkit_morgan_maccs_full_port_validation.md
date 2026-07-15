# RDKit Morgan/MACCS Exposed-Branch Validation

## Scope

This validation closes `dev/rdkit_morgan_maccs_full_port_plan.md` only for the
enumerated COSMolKit Morgan and MACCS branches below. It is not a blanket claim
for every RDKit fingerprint class, overload, molecule-preparation state, or
other fingerprint family.

The source oracle is the pinned RDKit tree under:

- `third_party/rdkit/Code/GraphMol/Fingerprints/FingerprintGenerator.*`
- `third_party/rdkit/Code/GraphMol/Fingerprints/FingerprintUtil.*`
- `third_party/rdkit/Code/GraphMol/Fingerprints/MorganGenerator.*`
- `third_party/rdkit/Code/GraphMol/Fingerprints/MorganFingerprints.*`
- `third_party/rdkit/Code/GraphMol/Fingerprints/MACCS.*`
- dependent RDKit CIPLabeler and PeriodicTable source used by Morgan chirality
  and isotope branches.

No compatibility threshold is accepted for this feature. A corpus result with
one differing bit is a failure. This report therefore records exact command
success and covered conditions, not pass percentages.

## Source Marker Closure

Morgan fingerprint behavior is source-backed for the exposed COSMolKit paths:

- RDKit fingerprint argument/default/JSON shape.
- Morgan argument wiring and generator ownership shape.
- connectivity atom invariants, feature atom invariants, and bond invariants.
- Morgan environment generation, duplicate suppression, from-atoms handling,
  ignore-atoms behavior matching the RDKit source, redundant-environment mode,
  count simulation, sparse count, sparse bit, hashed count, explicit bit vector,
  and additional output.
- include-chirality paths that require missing `_CIPCode` or `_CIPComputed`
  now route through the source-backed RDKit CIPLabeler port instead of a legacy
  shortcut.
- isotope mass and most-common-isotope data used by invariant and CIP paths are
  generated from RDKit `atomic_data.cpp`, not maintained as local fallback
  tables.

MACCS behavior is source-backed for the exposed COSMolKit path:

- RDKit `MACCS.cpp::Patterns` table bit numbering and SMARTS strings.
- substructure semantics required by the MACCS pattern table.
- RDKit direct-key logic for keys 1 through 166.
- RDKit 167-bit internal vector semantics with bit 0 unused.
- COSMolKit public 166-bit projection from RDKit bits 1 through 166.

Unsupported behavior policy remains unchanged: any option outside the modeled
source-backed surface returns a structured error instead of emitting
chemically meaningful-looking fingerprint bits. Morgan chirality requiring a
missing `_StereochemDone` preparation branch and supplied feature SMARTS are
examples of branches outside this validated boundary.

## Golden Conditions

Morgan golden tests are generated from pinned RDKit APIs and assert exact field
equality for:

- sparse-count output
- sparse-bit output
- hashed-count output
- explicit-bit output
- count simulation
- additional-output fields
- custom atom invariants
- custom bond invariants
- from-atoms behavior
- chirality behavior, including the CIPLabeler recomputation branch

MACCS golden tests are generated from `MACCSkeys.GenMACCSKeys` and assert exact
equality for:

- RDKit raw 167-bit vectors
- COSMolKit public 166-bit projected vectors
- targeted per-key fixtures
- daily small-corpus coverage
- strict 5000-row corpus coverage when `COSMOLKIT_PARITY_PROFILE=smiles_5000`
  is selected

## Validation Commands

Focused fingerprint parity:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_graph_feature_parity --test rdkit_morgan_fingerprint_parity --test rdkit_maccs_fingerprint_parity -- --nocapture
```

Result:

```text
rdkit_graph_feature_parity: 11 passed
rdkit_maccs_fingerprint_parity: 2 passed
rdkit_morgan_fingerprint_parity: 9 passed
```

Strict 5000-row Morgan/MACCS profile:

```bash
COSMOLKIT_PARITY_PROFILE=smiles_5000 cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_morgan_fingerprint_parity --test rdkit_maccs_fingerprint_parity -- --nocapture
```

Result:

```text
rdkit_maccs_fingerprint_parity: 2 passed
rdkit_morgan_fingerprint_parity: 9 passed
```

Full strict core validation:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict
```

Result:

```text
library tests: 2510 passed; 45 ignored
integration tests: passed
doc tests: 6 passed
```

Formatting:

```bash
cargo fmt --all
```

Result: completed successfully.

## No Overclaim

The Morgan and MACCS exposed paths must be described as exact RDKit parity for
the covered source-backed behavior or as unsupported. They must not be
described as approximate, correlation-compatible, mostly-compatible, or
percentage-compatible.

The accepted evidence is exact equality against pinned RDKit golden outputs and
source-backed implementation markers. `99.9%` agreement is not a passing state
for these fingerprint features.

This report does not cover `RDKFingerprintMol`/topological or Avalon. Their
former local path-hash approximations have been removed, and the public calls
return structured unsupported errors until the complete RDKit/Avalon source
paths are ported.
