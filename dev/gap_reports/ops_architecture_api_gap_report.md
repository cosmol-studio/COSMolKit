# Ops Architecture API Gap Report

## Confirmed Policy Gaps

- The redesigned core still exposes a public mutable bypass at
  `chemistry::stereo::assign_atom_cip_ranks_in_place(&mut Molecule)`.
- That bypass writes `_CIPRank` directly into atom properties through
  `topology_block_mut()` instead of going through a registered operation.
- This violates the redesigned-core policy that public mutation of existing
  molecules must go through registered operations.

## Architecture Decisions For This Pass

- Public mutable CIP-rank writeback is removed from the public API surface in
  this pass by demoting `assign_atom_cip_ranks_in_place()` to crate visibility.
- No public operation-backed replacement for CIP-rank property writeback is
  added in this pass.
- `_CIPRank` property writeback remains an internal helper path for depiction
  and SMILES-writing code until the project explicitly models persistent
  stereochemistry metadata as approved public molecule state.

## Convenience API Decisions

- Already-supported fragment behaviors should be reachable from `Molecule`
  through `fragments()` and `largest_fragment()`.
- Already-supported scaffold behaviors should be reachable from `Molecule`
  through `murcko_scaffold()` and `net_scaffold()`.
- Already-supported fingerprint behaviors should be reachable from `Molecule`
  through `avalon_fingerprint()`, `topological_fingerprint()`, and
  `maccs_fingerprint()`.
- Existing hash helpers should be reachable from `Molecule` through `hash()`
  and `hash_with_ranks()`.
- Existing PDB export should be reachable from `Molecule` through
  `to_pdb_block(conf_id, flavor)`.

## Stereochemistry Naming Decisions

- The current read-only stereochemistry entry point should be exposed under the
  clearer public name `perceive_stereochemistry()`.
- The existing `assign_stereochemistry()` symbol should remain temporarily as a
  deprecated compatibility wrapper with unchanged read-only behavior.
