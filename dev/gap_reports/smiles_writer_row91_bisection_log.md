# SMILES Writer Row 91 Bisection Log

## Probe 1

- Checked: `raw_parse`, `post_update_property_cache`, `post_symmetrize_sssr`, `post_assign_stereochemistry`, final writer output
- Last matching boundary: `post_assign_stereochemistry`
- First differing boundary: final plain non-isomeric output
- Reduced range: writer-side bond-direction cleanup / emission path
- Next probe: inspect traversal, then helper cleanup after `canonicalize_double_bond_directions_for_writer()`

## Probe 2

- Checked: traversal stack, post-traversal bond directions, post-helper bond directions, post-helper SMILES
- Last matching boundary: post-traversal stack capture
- First differing boundary: `canonicalize_double_bond_directions_for_writer()`
- Reduced range: helper cleanup of bond-dir state for row 91
- Next probe: compare helper cleanup against RDKit `clearBondDirs` / `removeRedundantBondDirSpecs`

## Probe 3

- Checked: RDKit source `clearBondDirs`, `removeUnwantedBondDirSpecs`, `removeRedundantBondDirSpecs`
- Last matching boundary: `clearDirection` counter updates
- First differing boundary: Rust `clear_bond_dirs_from_atom_for_writer()` decision logic
- Reduced range: `removeRedundantBondDirSpecs` atom-side clearing branch
- Next probe: align Rust helper with RDKit `clearBondDirs` exactly and lock row 91 with a dedicated regression test
# Probe 4
- Checked: Rust `after_double_bonds`, `after_remove_unwanted`, `after_remove_redundant`, `post_write`
- Last matching boundary: `after_remove_unwanted`
- First differing boundary: `after_remove_redundant`
- Reduced range: RDKit/Rust `removeRedundantBondDirSpecs` cleanup branch
- Next probe: upstream `post_canonicalize_fragment` on the same row 91 input

## Probe 5

- Checked: RDKit `removeRedundantBondDirSpecs` / `clearBondDirs` source vs Rust `remove_redundant_bond_dir_specs_for_writer()`
- Last matching boundary: `canHaveDirection(*tBond) && bondDirCounts[tBond->getIdx()]`
- First differing boundary: Rust adds an extra `any(other stereo double bond)` gate before calling `clear_bond_dirs_from_atom_for_writer()`
- Reduced range: `removeRedundantBondDirSpecs` outer trigger condition, not `clearBondDirs` itself
- Next probe: remove the extra Rust gate, rerun row 91 regression, and add a focused unit test that proves `clearBondDirsFromAtom` is entered without a neighboring stereo double bond

## Probe 6

- Checked: row 91 regression after removing the extra Rust gate
- Last matching boundary: `remove_redundant_bond_dir_specs_for_writer()` entry
- First differing boundary: `clear_bond_dirs_from_atom_for_writer()` still does not produce the RDKit cleanup effect
- Reduced range: `clearBondDirs` helper internals, not the outer `removeRedundantBondDirSpecs` gate
- Next probe: instrument `clear_bond_dirs_from_atom_for_writer()` with the row 91 working copy and inspect which bond/atom pair is skipped or cleared

## Probe 7

- Checked: manual `clear_bond_dirs_from_atom_for_writer()` calls on the row 91 working copy
- Last matching boundary: helper entry and the `atomDirCounts[fromAtom] >= 2` guard
- First differing boundary: the first eligible incident bond picked by the helper on the relevant atom
- Reduced range: `clearBondDirs` candidate selection / adjacency iteration order on the row 91 working copy
- Next probe: compare the row 91 atom-bond iteration order against RDKit working-copy state, then patch only if the order or the candidate-selection rule is source-different

## Probe 8

- Checked: row 91 working-copy bond 7 stereo state and CIP ranks before writer double-bond direction canonicalization
- Last matching boundary: `clearBondDirs` and `removeRedundantBondDirSpecs` source logic
- First differing boundary: Rust entered writer double-bond canonicalization with stale `Trans` stereo on bond 7 even though one double-bond end has two equal-rank methyl substituents
- Source basis: RDKit `Chirality::findAtomNeighborDirHelper` clears the neighbor list when two substituents on one double-bond end have identical CIP ranks, so no double-bond stereochemistry is assigned for that center
- Fix: reject such indistinguishable double-bond stereo state before writer direction canonicalization and clear it as non-stereo; row 91 now emits the RDKit plain non-isomeric string without slash/backslash bond directions
