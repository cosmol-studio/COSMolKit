# RDKit Fingerprint SMARTS/Substructure Gap Report

## Purpose

This report covers the SMARTS parsing and substructure-matching capability
required by:

- `MorganFeatureAtomInvGenerator::getAtomInvariants`
- `FingerprintUtil.cpp::defaultFeatureSmarts`
- `FingerprintUtil.cpp::getFeatureInvariants`
- `MACCS.cpp::Patterns`
- `MACCS.cpp::GenerateFP`

The target is exact RDKit source-backed behavior, not heuristic feature
classification.

## RDKit Behaviors That Must Be Covered

### Morgan feature invariants

RDKit uses the six default feature SMARTS patterns from
`FingerprintUtil.cpp::defaultFeatureSmarts`:

- Donor
- Acceptor
- Aromatic
- Halogen
- Basic
- Acidic

These are not local element buckets. They depend on SMARTS parsing and
substructure matching over query trees.

### MACCS key patterns

RDKit MACCS keys require:

- the full `Patterns` table from `MACCS.cpp`
- recursive SMARTS support through `$()`
- ring closures, branch syntax, and ring-bond queries
- atom and bond negation
- OR alternatives
- aromatic and aliphatic atom tokens
- hydrogen-count and degree constraints
- bond direction / aromatic / ring semantics where present in the pattern

## Current COSMolKit Support

### SMARTS parser

`crates/cosmolkit-core/src/search/smarts_parse.rs` already supports a broad
subset of SMARTS syntax:

- bracket atoms
- organic atoms and aromatic atoms
- AND / OR / NOT query trees
- ring closures
- recursive SMARTS labeling
- bond tokens such as `-`, `=`, `#`, `:`, `~`, `@`, `/`, `\\`

This is sufficient to parse many RDKit pattern strings, including the Morgan
default patterns and much of MACCS.

### Query predicates

`crates/cosmolkit-core/src/search/query.rs` already evaluates many RDKit-style
atom and bond predicates:

- atomic number and atom type
- aromaticity
- formal charge
- isotopes
- hydrogen counts
- explicit degree / total degree / substitution count
- ring membership / ring bond counts
- hybridization
- basic bond order and bond ring predicates

### Substructure matching

`crates/cosmolkit-core/src/search/substruct.rs` already provides VF2-based
matching with recursive SMARTS caching and query-aware predicate evaluation.

## Gaps Remaining

### 1. Morgan feature invariants are still heuristic in the fingerprint module

`crates/cosmolkit-core/src/properties/fingerprint.rs` still computes Morgan
feature invariants with local donor/acceptor/aromatic/halogen/basic/acidic
rules. That is not a source-backed replacement for RDKit
`defaultFeatureSmarts` + `getFeatureInvariants`.

### 2. RDKit recursive SMARTS and query semantics are not yet closed as a
parity claim

The parser and query engine can represent and evaluate many SMARTS forms, but
the current codebase still marks recursive SMARTS, selected molfile/query-code
branches, and parts of query evaluation as unfinished or explicitly unsupported.
That is acceptable for fail-closed behavior, but it is not yet a closed parity
surface for all Morgan/MACCS pattern strings.

### 3. MACCS requires exact pattern-table parity and exact match semantics

The current MACCS implementation in `fingerprint.rs` is a local heuristic rule
set. It does not provide:

- exact RDKit `Patterns` initialization
- exact 167-bit internal semantics
- exact key-by-key `GenerateFP` behavior
- exact first-match / pattern-match behavior for the source pattern table

### 4. Unsupported query branches must stay explicit

Any query primitive that is not fully reproduced must continue to fail closed
instead of producing a chemically plausible fallback. That includes recursive
SMARTS and any query code not yet mapped to exact RDKit semantics.

## Conclusion

COSMolKit has enough parser and query infrastructure to continue the port, but
the Morgan feature path and MACCS path are not source-complete yet. The next
required work is to replace heuristic feature classification with exact
RDKit-pattern matching and to close any remaining unsupported query branches
explicitly.
