# RDKit Fingerprint SMARTS/Substructure Boundary Audit

## Purpose

This audit records the SMARTS and substructure behavior reached by the
validated Morgan feature-invariant and MACCS branches. It must not be read as
a blanket closure of the public substructure API.

The source targets are:

- `third_party/rdkit/Code/GraphMol/Fingerprints/FingerprintUtil.cpp`
- `third_party/rdkit/Code/GraphMol/Fingerprints/MACCS.cpp`
- `third_party/rdkit/Code/GraphMol/Substruct/SubstructMatch.cpp`
- `third_party/rdkit/Code/GraphMol/Substruct/vf2.hpp`

## Validated Fingerprint Boundary

Morgan uses the six RDKit default feature SMARTS patterns for donor, acceptor,
aromatic, halogen, basic, and acidic invariants. MACCS uses the pinned RDKit
pattern table and direct-key logic. The committed Morgan and MACCS parity tests
compare exact output bits, raw/public projections, and the covered provenance
fields on targeted fixtures and the selected SMILES profiles.

The Rust implementation constructs these matchers through the source-backed
SMARTS/query path. Matcher construction is fallible. A parse failure is
returned as `FingerprintError::InvalidSmartsPattern`; an empty matcher is not
silently skipped; no local element/degree/ring heuristic is used as fallback.

## Deliberately Unfinished Boundaries

The following are separate from the validated fingerprint branch and remain
unfinished:

- the full public SMARTS-query API
- every `SubstructMatchParams` branch, including marker-open stereo/query
  behavior
- recursive SMARTS/query primitives not covered by the fingerprint pattern
  matrix
- user-supplied Morgan feature SMARTS generator patterns

Those branches must remain explicitly unsupported or marker-open. Passing the
Morgan/MACCS fingerprint matrix does not upgrade them to general substructure
parity.

## Policy

The previous report text describing Morgan feature invariants or MACCS as local
heuristics is obsolete. The current implementation either uses the validated
source-backed pattern path or returns a structured error. No similarity
correlation, partial key set, or 99.9% bit agreement is accepted as parity.
