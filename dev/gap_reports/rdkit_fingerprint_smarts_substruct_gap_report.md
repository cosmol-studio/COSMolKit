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

## Canonical SMARTS Boundary

All ordinary-molecule consumers now compile through `mol_from_smarts` into the
same query-bearing `Molecule` model. Fingerprints, descriptors, force fields,
coordinate templates, Marvin SMARTS lines, and SMARTSQ SGroups do not retain a
consumer-local parser, fallback decoder, or compatibility query graph.

The following remain deliberately outside this port:

- reaction SMARTS parsing and reaction-container matching
- cartridge/database query containers and serialization
- any `SubstructMatchParams` branch not proven by the pinned SMARTS corpus

Those excluded surfaces require separate source inventories and parity claims.

## Policy

The previous report text describing Morgan feature invariants or MACCS as local
heuristics is obsolete. The current implementation either uses the validated
source-backed pattern path or returns a structured error. No similarity
correlation, partial key set, or 99.9% bit agreement is accepted as parity.
