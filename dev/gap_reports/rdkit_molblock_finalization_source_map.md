# RDKit MolBlock Finalization Source Map

Step: `dev/archive/plans/rdkit_molblock_sdf_full_port_plan.md` Step 183 final update.

Scope audited:

- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp::finishMolProcessing`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp::ProcessMolProps`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp::processSGroups`
- RDKit `third_party/rdkit/Code/GraphMol/MolOps::expandAttachmentPoints`
- RDKit `third_party/rdkit/Code/GraphMol/Chirality.cpp::assignChiralTypesFromBondDirs`
- RDKit `third_party/rdkit/Code/GraphMol/Chirality.cpp::assignChiralTypesFrom3D`
- RDKit `third_party/rdkit/Code/GraphMol/Atropisomers.cpp::detectAtropisomerChirality`
- RDKit `third_party/rdkit/Code/GraphMol/Chirality.cpp::clearSingleBondDirFlags`
- RDKit `third_party/rdkit/Code/GraphMol/Chirality.cpp::detectBondStereochemistry`
- RDKit `third_party/rdkit/Code/GraphMol/MolOps::sanitizeMol`
- RDKit `third_party/rdkit/Code/GraphMol/AddHs.cpp::removeHs`
- RDKit `third_party/rdkit/Code/GraphMol/Chirality.cpp::assignStereochemistry`
- RDKit `third_party/rdkit/Code/GraphMol/QueryOps.cpp::completeMolQueries`
- COSMolKit `crates/cosmolkit-core/src/io/sdf/postprocess.rs`
- COSMolKit `crates/cosmolkit-core/src/operations/ops/sanitize_pipeline.rs`
- COSMolKit `crates/cosmolkit-core/src/chemistry/hydrogens.rs`
- COSMolKit `crates/cosmolkit-core/src/operations/ops/hydrogens.rs`

## Final Scope Statement

The public MolBlock/SDF finalization scope covered by this plan has no
remaining known behavioral gap against the mapped RDKit finalization functions
for the modeled parameter space:

```text
sanitize in {true,false}
removeHs in {true,false}
strictParsing in {true,false}
public V2000 and V3000 MolBlock/SDF records
delayed sanitize and delayed removeHs public workflows
```

This is a behavior finding for the modeled public reader/finalization surface.
It is not a blanket claim that every adjacent helper in COSMolKit has a closed
second-axis performance marker. Several copied source blocks intentionally
remain `RDKit✔️❌` where the Rust implementation preserves behavior through
typed operation machinery, remapping planners, or shared sanitize/stereo
pipelines but has different allocation, routing, or operation-contract costs.

The remaining `RDKit❗` markers in nearby wrapper comments are not being
silently upgraded by this report. They mark either shared-helper delegation
whose global behavior is broader than the MolBlock/SDF finalization state
space, or adjacent parsing/SGroup subpaths not claimed beyond the modeled
fixtures. Within the plan's public finalization matrix, the generated RDKit
goldens, focused finalization tests, and full validation runs leave no known
behavioral gap.

## Finalization Order Map

RDKit `finishMolProcessing` order and COSMolKit final status:

| RDKit order | COSMolKit target | Final finding |
|---|---|---|
| Return on null molecule | `finish_mol_processing` type boundary | No known behavioral gap. Rust cannot receive a null `Molecule`; nullable RDKit return paths are represented at caller `Result`/`Option` boundaries. |
| Clear atom and bond bookmarks | No persistent bookmark model | No known behavioral gap in modeled scope. COSMolKit does not persist RDKit atom/bond bookmark state; tests verify ordinary atom, bond, and molecule properties are not treated as bookmarks or accidentally cleared. |
| Expand attachment points when requested | `expand_attachment_points` called by `finish_mol_processing` | No known behavioral gap in modeled scope. The branch is source-backed and parameter-driven; behavior is covered through MolBlock/SDF parity and finalization validation. |
| Calculate explicit valence before mol props | `calculate_explicit_valence_before_mol_props` | No known behavioral gap. The explicit-valence snapshot is computed before `process_mol_props`, matching RDKit ordering for `molTotValence` handling. |
| Process molfile props and SGroups | `process_mol_props`; `process_sgroups` | No known behavioral gap in modeled scope. `molSubstCount`, `molTotValence`, `_ZBO_H`, MRV coordinate-bond, MRV implicit-H, ZBO/ZCH/HYD, and SMARTSQ/SQ handling are source-backed and tested. |
| Perceive atom chirality before H removal | `assign_chiral_types_from_bond_dirs`; `assign_chiral_types_from_3d` | No known behavioral gap in modeled scope. The 2D wedge/dash and 3D conformer branches run before hydrogen removal and before final sanitize assignment, including explicit/wedged H and default-conformer `is3D` behavior. |
| Detect atropisomer chirality | `detect_atropisomer_chirality` | No known behavioral gap in modeled scope. 2D wedge/hash, 3D geometry, non-atropisomer rejection, inconsistent wedging rejection, and RDKit macrocycle fixture behavior are covered. |
| Clear single-bond directions | `clear_single_bond_dir_flags` | No known behavioral gap. Single-bond wedge, dash, unknown-stereo preservation, non-single preservation, and `onlyWedgeFlags` mode are source-backed and tested. |
| If `sanitize && removeHs`: cleanup, detect bond stereo, remove Hs | `sanitize_cleanup_for_sdf_remove_hs`; `detect_bond_stereochemistry`; `remove_hs_after_sdf_parse` | No known behavioral gap in modeled scope. Cleanup precedes double-bond stereo detection, and stereo detection precedes hydrogen removal, preserving imine/stereogenic H behavior. |
| If `sanitize && !removeHs`: full sanitize, detect bond stereo | `sanitize_after_sdf_parse`; `detect_bond_stereochemistry` | No known behavioral gap in modeled scope. Full sanitize and post-sanitize bond stereo detection are source-backed through the registered sanitize/stereo machinery and parity tests. |
| If `!sanitize`: detect bond stereo only | `detect_bond_stereochemistry` | No known behavioral gap in modeled scope. Unsanitized finalization preserves raw valence/aromaticity cleanup state and omits final CIP-style stereochemistry assignment while still setting required direction data. |
| If sanitized: assign stereochemistry | `assign_stereochemistry_after_sdf_parse` plus shared stereo helpers | No known behavioral gap in modeled finalization scope. The wrapper delegates to source-backed stereo cleanup and double-bond assignment helpers; global `assignStereochemistry` breadth outside this finalization state space is not claimed here. |
| Complete query scans | `complete_mol_queries` | No known behavioral gap in modeled scope. `_NeedsQueryScan` is cleared and recursive query-scan predicates produced by MolBlock/SDF parsing are completed from ring-bond counts while existing bond queries are preserved. |

## Function Map

| RDKit function | COSMolKit target | Final finding |
|---|---|---|
| `finishMolProcessing` | `finish_mol_processing` | No known behavioral gap in modeled public scope. The high-level branch order and parameter-dependent behavior are source-anchored and covered by focused finalization tests plus MolBlock/SDF parity goldens. |
| `ProcessMolProps` | `process_mol_props` | No known behavioral gap in modeled scope. Atom-level substitution-count and total-valence properties are consumed like RDKit, including `-1`, `-2`, `>=6`, V2000 `15`, V3000 `-1`, drawn-valence overflow, and `_ZBO_H` skip behavior. |
| `processSGroups` finalization hooks | `process_sgroups` and SGroup helpers | No known behavioral gap in modeled scope. Modeled DAT SGroup finalization extensions are applied and then removed from persistent SGroup state in RDKit order. Unsupported/unmodeled SGroup features are not converted into silent chemistry. |
| `expandAttachmentPoints` | `expand_attachment_points` | No known behavioral gap in modeled scope. Attachment-point expansion is gated by `SdfReadParams::expand_attachment_points` and source-backed in finalization order. |
| `assignChiralTypesFromBondDirs` | `assign_chiral_types_from_bond_dirs` | No known behavioral gap in modeled scope. It covers 2D wedge/dash, explicit wedged H, implicit-H promotion, `replaceExistingTags=true`, and final direction clearing order. |
| `assignChiralTypesFrom3D` | `assign_chiral_types_from_3d` | No known behavioral gap in modeled scope. It uses the default conformer, respects `is3D`, updates valence cache before assignment, supports modeled tetrahedral and non-tetrahedral branches, and records non-explicit 3D chirality where RDKit does. |
| `detectAtropisomerChirality` | `detect_atropisomer_chirality` | No known behavioral gap in modeled scope. Source-backed candidate collection, direction interpretation, 2D/3D assignment, non-candidate no-op behavior, and inconsistent-wedging rejection are covered. |
| `clearSingleBondDirFlags` | `clear_single_bond_dir_flags` | No known behavioral gap. Source markers are `RDKit✔️✔️`; tests cover wedge, dash, unknown, non-single, and alternate mode behavior. |
| `detectBondStereochemistry` | `detect_bond_stereochemistry` | No known behavioral gap in modeled scope. It maps RDKit's conformer guard and `setDoubleBondNeighborDirections` behavior through the shared stereo implementation and is tested before and after sanitize/removeHs branches. |
| `sanitizeMol(... SANITIZE_CLEANUP)` | `sanitize_cleanup_for_sdf_remove_hs`; `sanitize_cleanup_assignment` | No known behavioral gap in modeled scope. The cleanup subset is routed through the registered sanitize operation and source-backed cleanup assignment, including nitro/charge/order cleanup needed before finalization stereo detection. |
| `sanitizeMol(*res)` | `sanitize_after_sdf_parse`; sanitize pipeline helpers | No known behavioral gap in modeled scope. Full sanitize behavior used by MolBlock/SDF finalization and delayed sanitize parity is covered by generated RDKit goldens and focused tests. |
| `removeHs(*res)` | `remove_hs_after_sdf_parse`; `without_hydrogens_apply`; `sanitize_after_remove_hs_removal` | No known behavioral gap in modeled scope. Hydrogen-removal predicates, SGroup guards, isotope tracking, chiral/stereo updates, atom/bond/SGroup remapping, and post-removal sanitize behavior are source-backed in the shared hydrogen-removal planner. |
| `assignStereochemistry(*res, true, true, true)` | `assign_stereochemistry_after_sdf_parse` and shared stereo helpers | No known behavioral gap in modeled finalization scope. The report does not claim unsupported global RDKit stereochemistry branches outside the parser finalization matrix. |
| `completeMolQueries` | `complete_mol_queries`; `complete_query_scan_predicates`; `atom_ring_bond_counts` | No known behavioral gap in modeled scope. The only MolBlock/SDF magic-value query currently produced by parsing is ring-bond-count scan completion; recursive query traversal and `_NeedsQueryScan` clearing are implemented for that modeled query state. |

## Evidence

Source anchors:

- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::finish_mol_processing`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::calculate_explicit_valence_before_mol_props`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::process_mol_props`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::process_sgroups`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::assign_chiral_types_from_bond_dirs`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::assign_chiral_types_from_3d`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::detect_atropisomer_chirality`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::clear_single_bond_dir_flags`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::detect_bond_stereochemistry`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::sanitize_cleanup_for_sdf_remove_hs`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::sanitize_after_sdf_parse`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::remove_hs_after_sdf_parse`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::assign_stereochemistry_after_sdf_parse`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::complete_mol_queries`
- `crates/cosmolkit-core/src/operations/ops/sanitize_pipeline.rs::sanitize_cleanup_assignment`
- `crates/cosmolkit-core/src/operations/ops/sanitize_pipeline.rs::sanitize_cleanup_atropisomers_assignment`
- `crates/cosmolkit-core/src/operations/ops/sanitize_pipeline.rs::sanitize_cleanup_chirality_assignment`
- `crates/cosmolkit-core/src/chemistry/hydrogens.rs::remove_hs_assignment`
- `crates/cosmolkit-core/src/chemistry/hydrogens.rs::should_remove_h`
- `crates/cosmolkit-core/src/operations/ops/hydrogens.rs::without_hydrogens_apply`
- `crates/cosmolkit-core/src/operations/ops/hydrogens.rs::sanitize_after_remove_hs_removal`

Focused tests:

- `crates/cosmolkit-core/src/io/sdf/postprocess.rs` covers mol-prop finalization, query completion, clear-single-bond-dir behavior, 2D/3D atom chirality order, atropisomer detection, branch ordering for `sanitize/removeHs`, unsanitized finalization behavior, and final stereochemistry assignment effects.
- `crates/cosmolkit-core/src/io/sdf.rs` covers reader-visible finalization outcomes, including sanitize/removeHs combinations, SDF/MolBlock parameter propagation, bond stereochemistry, atropisomer fixtures, and post-parse cleanup behavior.
- `crates/cosmolkit-core/src/chemistry/hydrogens.rs` covers RDKit `RemoveHsParameters`, `shouldRemoveH`, SGroup guards, isotope tracking, stereo/chirality adjustments, explicit-H count updates, and assignment planning.
- `crates/cosmolkit-core/tests/rdkit_molfile_read_parity.rs` compares generated RDKit MolBlock/molfile golden rows against COSMolKit final atom, bond, conformer, SGroup, property, SMILES, delayed sanitize, delayed removeHs, and error behavior.
- `crates/cosmolkit-core/tests/rdkit_sdf_read_parity.rs` covers the same finalization surface through SDF supplier and dataset paths.

Validation performed after the full port:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molfile_read_parity -- --nocapture
cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_sdf_read_parity -- --nocapture
cargo test -p cosmolkit-core --features op-contracts-strict
cargo test --workspace --features cosmolkit-core/op-contracts-strict
.venv/bin/pytest
```

The final full validation completed successfully:

```text
cargo test -p cosmolkit-core --features op-contracts-strict: passed
cargo test --workspace --features cosmolkit-core/op-contracts-strict: passed
.venv/bin/pytest: 491 passed, 37 skipped
```

## Residual Non-Claims

- This report does not claim public SCSR molfile parsing parity.
- This report does not claim RDKit warning-log stream parity; strict/non-strict observable success or structured error behavior is the covered public contract.
- This report does not batch-upgrade copied source markers. Any remaining `RDKit✔️❌` or `RDKit✔️❗` second-axis marker requires its own performance/complexity review before changing.
- This report does not claim every global RDKit `assignStereochemistry`, `sanitizeMol`, or `removeHs` branch outside the MolBlock/SDF finalization matrix. It records no known behavioral gap for the source-backed finalization behavior exercised by this plan.
