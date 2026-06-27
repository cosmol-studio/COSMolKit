# RDKit MolBlock and SDF Full Port Plan

## Scope

This plan ports the public MolBlock and SDF reader/finalization behavior to the
pinned RDKit source implementation. It covers `sanitize=false`, `removeHs`,
delayed `sanitize`, delayed `removeHs`, SDF stream behavior, indexed SDF
behavior, Python parameter propagation, docs, examples, and generated stubs.

The plan does not allow heuristic approximations, wrapper-only rejection, or
structured-unsupported completion for RDKit behavior in this scope. If a step is
too large to port completely, regenerate this plan by splitting that step into
smaller complete source-backed behavior steps before implementing anything.

## Source Authorities

- `third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp`
- `third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h`
- `third_party/rdkit/Code/GraphMol/FileParsers/FileParserUtils.h`
- `third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp`
- `third_party/rdkit/Code/GraphMol/FileParsers/MolSupplier.h`
- `third_party/rdkit/Code/GraphMol/FileParsers/MolSupplier.v1API.h`
- `third_party/rdkit/Code/GraphMol/Chirality.cpp`
- `third_party/rdkit/Code/GraphMol/Atropisomers.cpp`
- `third_party/rdkit/Code/GraphMol/QueryOps.cpp`
- `third_party/rdkit/Code/GraphMol/MolOps.cpp`
- `third_party/rdkit/Code/GraphMol/AddHs.cpp`

## COSMolKit Target Functions

- `crates/cosmolkit-core/src/io/molfile.rs::read_mol_record_from_str_with_params`
- `crates/cosmolkit-core/src/io/molfile.rs::read_mol_file_with_params`
- `crates/cosmolkit-core/src/io/sdf.rs::SdfReadParams`
- `crates/cosmolkit-core/src/io/sdf.rs::SdfReader::with_params`
- `crates/cosmolkit-core/src/io/sdf.rs::SdfReader::next_record`
- `crates/cosmolkit-core/src/io/sdf.rs::SdfDataset::open_with_params`
- `crates/cosmolkit-core/src/io/sdf.rs::SdfDataset::record_with_params`
- `crates/cosmolkit-core/src/io/sdf.rs::read_sdf_from_str_with_params`
- `crates/cosmolkit-core/src/io/sdf.rs::read_indexed_sdf_record`
- `crates/cosmolkit-core/src/io/sdf.rs::read_next_sdf_record`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_sdf_record_text`
- `crates/cosmolkit-core/src/io/sdf.rs::mol_from_mol_block`
- `crates/cosmolkit-core/src/io/sdf.rs::mol_from_mol_data_stream`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_mol_header`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_counts_line`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_v2000_ctab`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_v3000_ctab`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_sdf_data_fields`
- `crates/cosmolkit-core/src/io/sdf.rs::process_sdf_property_lists`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::process_mol_props`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::complete_mol_queries`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::clear_single_bond_dir_flags`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::sanitize_cleanup_for_sdf_remove_hs`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::sanitize_after_sdf_parse`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::remove_hs_after_sdf_parse`
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs::finish_mol_processing`
- `crates/cosmolkit-core/src/operations/ops/sanitize_pipeline.rs::sanitize_cleanup_assignment`
- `crates/cosmolkit-core/src/operations/ops/sanitize_pipeline.rs::sanitize_cleanup_atropisomers_assignment`
- `crates/cosmolkit-core/src/operations/ops/sanitize_pipeline.rs::sanitize_cleanup_chirality_assignment`
- `crates/cosmolkit-core/src/operations/ops/hydrogens.rs::without_hydrogens_apply`
- `crates/cosmolkit-core/src/operations/ops/hydrogens.rs::sanitize_after_remove_hs_removal`
- `python/src/lib.rs::reject_unsanitized_mol_reader`
- `python/src/lib.rs::Molecule::read_sdf`
- `python/src/lib.rs::Molecule::read_mol`
- `python/src/lib.rs::Molecule::read_mol_from_str`
- `python/src/lib.rs::Molecule::read_sdf_from_str`
- `python/src/lib.rs::MoleculeBatch::read_sdf`
- `python/src/lib.rs::MoleculeBatch::read_sdf_records_from_str`
- `python/src/lib.rs::PySdfDataset::open`
- `python/src/lib.rs::PySdfReader::open`

## Function Mapping

| RDKit function or type | COSMolKit target |
|---|---|
| `RDKit::v2::FileParsers::MolFileParserParams` | `SdfReadParams` |
| `RDKit::MolBlockToMol(molBlock, sanitize, removeHs, strictParsing)` | `read_mol_record_from_str_with_params` and Python `Molecule.read_mol_from_str` |
| `RDKit::MolFileToMol(fName, sanitize, removeHs, strictParsing)` | `read_mol_file_with_params` and Python `Molecule.read_mol` |
| `RDKit::v2::FileParsers::MolFromMolBlock(molBlock, params)` | `mol_from_mol_block` and `read_mol_record_from_str_with_params` |
| `RDKit::v2::FileParsers::MolFromMolFile(fName, params)` | `read_mol_file_with_params` |
| `RDKit::v2::FileParsers::MolFromMolDataStream(inStream, line, params)` | `mol_from_mol_data_stream`, `parse_mol_header`, `parse_counts_line`, `parse_v2000_ctab`, and `parse_v3000_ctab` |
| `RDKit::FileParserUtils::finishMolProcessing(res, chiralityPossible, params)` | `finish_mol_processing` |
| `RDKit::FileParserUtils::ProcessMolProps(res)` | `process_mol_props` |
| `RDKit::MolOps::expandAttachmentPoints(*res)` | new or existing source-backed attachment-point helper called by `finish_mol_processing` |
| `RDKit::MolOps::assignChiralTypesFromBondDirs(*res, confId, true)` | source-backed helper called by `finish_mol_processing` |
| `RDKit::MolOps::assignChiralTypesFrom3D(*res, confId, true)` | source-backed helper called by `finish_mol_processing` |
| `RDKit::Atropisomers::detectAtropisomerChirality(*res, &conf)` | source-backed helper called by `finish_mol_processing` |
| `RDKit::MolOps::clearSingleBondDirFlags(*res)` | `clear_single_bond_dir_flags` |
| `RDKit::MolOps::sanitizeMol(*res, failedOp, MolOps::SANITIZE_CLEANUP)` | `sanitize_cleanup_for_sdf_remove_hs` |
| `RDKit::MolOps::detectBondStereochemistry(*res)` | source-backed helper called by `finish_mol_processing` |
| `RDKit::MolOps::removeHs(*res)` | `remove_hs_after_sdf_parse` backed by `without_hydrogens_apply` |
| `RDKit::MolOps::sanitizeMol(*res)` | `sanitize_after_sdf_parse` |
| `RDKit::MolOps::assignStereochemistry(*res, true, true, true)` | source-backed helper called by `finish_mol_processing` |
| `RDKit::QueryOps::completeMolQueries(res)` | `complete_mol_queries` |
| `RDKit::v2::FileParsers::ForwardSDMolSupplier::_next()` | `read_next_sdf_record` |
| `RDKit::v2::FileParsers::ForwardSDMolSupplier::readMolProps(ROMol&)` | `parse_sdf_data_fields` and `process_sdf_property_lists` |
| `RDKit::v2::FileParsers::SDMolSupplier` | `SdfDataset`, `SdfReader`, `MoleculeBatch.read_sdf`, and Python `SdfDataset` |

## Required Test Matrix

Every parity row must compare RDKit and COSMolKit atom count, bond count,
atomic number, isotope, formal charge, radical electrons, explicit hydrogen
count, implicit hydrogen visibility when exposed, no-implicit flag, aromatic
flags, atom map number, atom properties, atom query predicates, bond order, bond
direction, bond stereo, bond stereo atoms, bond properties, bond query
predicates, chiral tags, chiral permutation when exposed, molfile chiral flag,
molfile info/comment/name properties, `_NeedsQueryScan` effects, substance
groups, molecule properties, SDF data fields, conformer coordinates, conformer
dimensionality, source coordinate dimensionality, canonical SMILES,
noncanonical SMILES, delayed sanitize result, delayed removeHs result, and error
kind/message category.

| Input API | RDKit baseline | COSMolKit API | Parameter rows |
|---|---|---|---|
| MolBlock string | `Chem.MolFromMolBlock(text, sanitize=s, removeHs=h, strictParsing=p)` | `read_mol_record_from_str_with_params` and `Molecule.read_mol_from_str` | `s in {true,false}`, `h in {true,false}`, `p in {true,false}` |
| Molfile path | `Chem.MolFromMolFile(path, sanitize=s, removeHs=h, strictParsing=p)` | `read_mol_file_with_params` and `Molecule.read_mol` | `s in {true,false}`, `h in {true,false}`, `p in {true,false}` |
| SDF first record string | `Chem.ForwardSDMolSupplier(BytesIO(text), sanitize=s, removeHs=h, strictParsing=p)` first record | `read_sdf_from_str_with_params` and `Molecule.read_sdf_from_str` | `s in {true,false}`, `h in {true,false}`, `p in {true,false}`, property-list processing on/off |
| SDF stream | `Chem.ForwardSDMolSupplier(stream, sanitize=s, removeHs=h, strictParsing=p)` | `SdfReader::with_params(...).next_record()` and Python `SdfReader.open(...).batches(...)` | valid record, invalid record, empty record, trailing text, missing `$$$$`, sanitize failure |
| SDF indexed dataset | `Chem.SDMolSupplier(path, sanitize=s, removeHs=h, strictParsing=p)` | `SdfDataset::open_with_params`, `SdfDataset::record_with_params`, Python `SdfDataset.open`, and `MoleculeBatch.read_sdf(index="memory")` | same as SDF stream plus random access and sliced batches |
| Delayed sanitize | `mol = Chem.MolFromMolBlock(text, sanitize=False, removeHs=h); Chem.SanitizeMol(mol)` | `mol = read_mol(... sanitize=false, remove_hs=h); mol.sanitize()` | `h in {true,false}` across V2000, V3000, 2D, 3D, query, SGroup, and stereo fixtures |
| Delayed removeHs | `mol = Chem.MolFromMolBlock(text, sanitize=s, removeHs=False); Chem.RemoveHs(mol)` | `mol = read_mol(... sanitize=s, remove_hs=false).without_hydrogens()` | `s in {true,false}` across wedge H, imine H, isotopic H, mapped H, SGroups, and query H |
| Default finalization | `Chem.MolFromMolBlock(text)` and `Chem.SDMolSupplier(path)[i]` | default `Molecule.read_mol`, `Molecule.read_sdf`, `SdfReader`, and `SdfDataset` | V2000, V3000, 2D, 3D, chiral, double-bond stereo, atropisomer, query, SGroup, property-list fixtures |

Fixture categories:

- V2000 atom and bond fields.
- V2000 counts-line strict and non-strict parsing.
- V2000 CTAB records with charge, isotope, radical, atom map, alias, atom list, old atom list, R label, zero-order bond, query bond, wedged bond, crossed double bond, and stereo care.
- V2000 SGroup records covering `STY`, `SST`, `SLB`, `SCN`, `SDS`, `SMT`, `SDI`, `SBV`, `SDT`, `SDD`, `SCD`, `SED`, `SPL`, `SNC`, `SAP`, `SCL`, `SBT`, `ZBO`, `ZCH`, and `HYD`.
- V3000 CTAB records with atom props, bond props, `ENDPTS`, `ATTACH`, enhanced stereo collections, object3D blocks, and SGroup blocks.
- 2D chirality from wedge/dash before hydrogen removal.
- 3D chirality from coordinates before sanitize.
- Double-bond stereo requiring detection before hydrogen removal.
- Atropisomer fixtures from RDKit molfile tests.
- Query molfile records requiring `QueryOps::completeMolQueries`.
- SDF data fields, property lists, blank values, repeated data fields, malformed data field blocks, and malformed records.
- Bad input files, immediate EOF, missing molecule blocks, strict parsing errors, and non-strict parsing recoveries.

## Execution Contract

- This plan is for sequential continuous execution.
- "Sequential continuous execution" means execute one step at a time in order and continue to the next unchecked step until all steps are completed, blocked, or the user interrupts.
- It does not mean executing steps in unordered batches or postponing validation for a batch of changes.
- Execute unchecked steps in order.
- Continue executing the plan until all steps are completed, blocked, or the user interrupts.
- Do not stop after every step unless the plan explicitly says to stop.
- Mark each completed step by changing only its `[ ]` to `[x]`.
- Never execute unchecked steps out of order.
- Never summarize, skip, or reinterpret later unchecked steps.
- Never treat a required reading step as “already read”.
- Do not assume the agent is diligent.
- Do not assume the model context is long enough.
- Do not rely on memory from previous turns when a required reading step is present.
- Every real task step must be immediately preceded by reading:
  `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`.
- The reading step must explicitly reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
- `Implement`, `Port`, `Modify`, `Update`, and `Fix` steps must produce a concrete artifact.
- `Audit` steps must produce a written gap report and must not replace implementation steps.
- If a step adds or updates tests, the next real task after the required reading step must run the most specific relevant test command for those tests.
- Do not defer tests added for one behavior to a final whole-plan validation step.
- Final whole-plan validation is still required when the plan changes code, but it does not replace immediate targeted validation after test-writing steps.
- If the plan violates this contract, regenerate the plan before doing any work.
- Copying C++ comments, adding a dispatch stub, or adding placeholder branches is not a completed `Port` step.
- Do not use “smallest subpart”, skeleton code, dispatch-only code, placeholder code, TODO-only code, or partial porting as completion.
- If a step is too large to complete fully, regenerate the plan by splitting that step into smaller full-port steps.
- Each split step must still represent a complete source-backed behavior, not a placeholder.
- A `Port` step is complete only when the selected RDKit behavior is fully implemented, source-backed, tested, and no placeholder remains.
- During porting, do not execute any git operation unless the user explicitly asks for it.
- This includes read-only git operations such as `git status`, `git diff`, and `git log`.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit RDKit `MolFileParser.cpp::MolFromMolDataStream`, `MolFileParser.cpp::MolFromMolBlock`, `MolFileParser.cpp::MolFromMolFile`, `FileParsers.h::MolBlockToMol`, and `FileParsers.h::MolFileToMol` against COSMolKit `mol_from_mol_data_stream`, `mol_from_mol_block`, `read_mol_record_from_str_with_params`, `read_mol_file_with_params`, and `SdfReadParams` and write `dev/gap_reports/rdkit_molblock_reader_source_map.md`.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Audit RDKit `FileParserUtils::finishMolProcessing`, `FileParserUtils::ProcessMolProps`, `MolOps::expandAttachmentPoints`, `MolOps::assignChiralTypesFromBondDirs`, `MolOps::assignChiralTypesFrom3D`, `Atropisomers::detectAtropisomerChirality`, `MolOps::clearSingleBondDirFlags`, `MolOps::detectBondStereochemistry`, `MolOps::sanitizeMol`, `MolOps::removeHs`, `MolOps::assignStereochemistry`, and `QueryOps::completeMolQueries` against COSMolKit `finish_mol_processing`, `process_mol_props`, `complete_mol_queries`, `clear_single_bond_dir_flags`, `sanitize_cleanup_for_sdf_remove_hs`, `sanitize_after_sdf_parse`, `remove_hs_after_sdf_parse`, `sanitize_cleanup_assignment`, `sanitize_cleanup_atropisomers_assignment`, `sanitize_cleanup_chirality_assignment`, `without_hydrogens_apply`, and `sanitize_after_remove_hs_removal` and write `dev/gap_reports/rdkit_molblock_finalization_source_map.md`.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Audit RDKit `ForwardSDMolSupplier::_next`, `ForwardSDMolSupplier::readMolProps`, `ForwardSDMolSupplier::checkForEnd`, `SDMolSupplier`, `MolSupplier.v1API.h`, and `MultithreadedSDMolSupplier` parameter propagation against COSMolKit `read_next_sdf_record`, `read_indexed_sdf_record`, `parse_sdf_record_text`, `parse_sdf_data_fields`, `process_sdf_property_lists`, `SdfReader::with_params`, `SdfReader::next_record`, `SdfDataset::open_with_params`, `SdfDataset::record_with_params`, and `MoleculeBatch::read_sdf_records_from_reader_with_options` and write `dev/gap_reports/rdkit_sdf_supplier_source_map.md`.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Update `tests/scripts/gen_rdkit_molfile_read_golden.py` to generate RDKit baseline rows for `Chem.MolFromMolBlock`, `Chem.MolFromMolFile`, delayed `Chem.SanitizeMol`, delayed `Chem.RemoveHs`, `sanitize`, `removeHs`, `strictParsing`, and failure cases.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Run `.venv/bin/python tests/scripts/gen_rdkit_molfile_read_golden.py --output tests/golden/molfile_read.jsonl`.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Update `tests/scripts/gen_rdkit_sdf_read_golden.py` to generate RDKit baseline rows for `Chem.ForwardSDMolSupplier`, `Chem.SDMolSupplier`, data fields, property lists, malformed records, `sanitize`, `removeHs`, `strictParsing`, and failure cases.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Run `.venv/bin/python tests/scripts/gen_rdkit_sdf_read_golden.py --output tests/golden/sdf_read.jsonl`.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Add `crates/cosmolkit-core/tests/rdkit_molfile_read_parity.rs` to compare every `tests/golden/molfile_read.jsonl` row against `read_mol_record_from_str_with_params`, `read_mol_file_with_params`, delayed `Molecule::sanitize`, delayed `without_hydrogens_with_sanitize`, atom fields, bond fields, conformer fields, molecule properties, SGroups, and SMILES outputs.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molfile_read_parity -- --nocapture`.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Add `crates/cosmolkit-core/tests/rdkit_sdf_read_parity.rs` to compare every `tests/golden/sdf_read.jsonl` row against `read_sdf_from_str_with_params`, `SdfReader::with_params`, `SdfReader::next_record`, `SdfDataset::open_with_params`, `SdfDataset::record_with_params`, `MoleculeBatch::read_sdf_records_from_reader_with_options`, data fields, property lists, atom fields, bond fields, conformer fields, molecule properties, SGroups, and SMILES outputs.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_sdf_read_parity -- --nocapture`.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Port RDKit `MolFileParser.cpp::MolFromMolBlock` into COSMolKit `mol_from_mol_block` and `read_mol_record_from_str_with_params` with source comments and marker statuses.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Add focused Rust tests in `crates/cosmolkit-core/src/io/molfile.rs` for RDKit `MolFromMolBlock` unread trailing text, empty input, missing counts line, `strictParsing`, `sanitize`, and `removeHs` behavior.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict io::molfile -- --nocapture`.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Port RDKit `MolFileParser.cpp::MolFromMolFile` into COSMolKit `read_mol_file_with_params` with source comments and marker statuses.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Add focused Rust tests in `crates/cosmolkit-core/src/io/molfile.rs` for RDKit `MolFromMolFile` file-open failure, file-read failure, empty file, and single-record molfile behavior.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict io::molfile -- --nocapture`.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Port RDKit `MolFileParser.cpp::MolFromMolDataStream` header, counts-line, `strictParsing`, V2000 dispatch, V3000 dispatch, `fileComplete`, and missing `M  END` behavior into COSMolKit `mol_from_mol_data_stream`.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Add focused Rust tests in `crates/cosmolkit-core/src/io/sdf.rs` for `mol_from_mol_data_stream` header, counts-line, `strictParsing`, V2000 dispatch, V3000 dispatch, `fileComplete`, and missing `M  END` behavior.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict io::sdf::tests::mol_from_mol_data_stream -- --nocapture`.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Port RDKit `FileParserUtils::ProcessMolProps` into COSMolKit `process_mol_props` with source comments and marker statuses.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Add focused Rust tests for `process_mol_props` covering `_MolFileChiralFlag`, `molTotValence`, `molSubstCount`, `molUnsaturated`, `molRingBondCount`, atom query conversion, bond query conversion, and absence of SDF data-field coupling.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict process_mol_props -- --nocapture`.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Port RDKit `MolOps::expandAttachmentPoints` behavior used by `finishMolProcessing` into the COSMolKit helper called from `finish_mol_processing`.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Add focused Rust tests for `expand_attachment_points` behavior from RDKit `MolOps::expandAttachmentPoints` covering atom `molAttachPoint`, query dummy creation, coordinate creation, and `SdfReadParams::expand_attachment_points`.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict expand_attachment_points -- --nocapture`.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Port RDKit `MolOps::assignChiralTypesFromBondDirs` from `Chirality.cpp` into the source-backed helper called by `finish_mol_processing`.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Add focused Rust tests for RDKit `assignChiralTypesFromBondDirs` covering 2D wedged atoms, dash bonds, explicit wedged hydrogens, implicit hydrogen promotion, `replaceExistingTags=true`, and `sanitize=false`.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict assign_chiral_types_from_bond_dirs -- --nocapture`.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Port RDKit `MolOps::assignChiralTypesFrom3D` from `Chirality.cpp` into the source-backed helper called by `finish_mol_processing`.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Add focused Rust tests for RDKit `assignChiralTypesFrom3D` covering tetrahedral centers, non-tetrahedral centers, implicit hydrogen, existing explicit tags, false 3D tags, and `replaceExistingTags=true`.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict assign_chiral_types_from_3d -- --nocapture`.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Port RDKit `Atropisomers.cpp::detectAtropisomerChirality` into the source-backed helper called by `finish_mol_processing`.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Add focused Rust tests for RDKit `detectAtropisomerChirality` covering 2D wedge atropisomers, 3D atropisomers, non-atropisomer rejection, inconsistent wedging, and RDKit fixture parity.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict atropisomer -- --nocapture`.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Port RDKit `MolOps::clearSingleBondDirFlags` from `Chirality.cpp` into COSMolKit `clear_single_bond_dir_flags`.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Add focused Rust tests for RDKit `clearSingleBondDirFlags` covering wedge clearing, dash clearing, non-single-bond retention, `onlyWedgeFlags=false`, and placement after atom stereo perception.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict clear_single_bond_dir_flags -- --nocapture`.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Port RDKit `MolOps::detectBondStereochemistry` from `Chirality.cpp` into the source-backed helper called by `finish_mol_processing`.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Add focused Rust tests for RDKit `detectBondStereochemistry` covering 2D coordinates, 3D coordinates, imine hydrogens, crossed double bonds, stereo atom assignment, and detection before hydrogen removal.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict detect_bond_stereochemistry -- --nocapture`.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Port RDKit `FileParserUtils::finishMolProcessing` pre-sanitize branch into COSMolKit `finish_mol_processing`, including `clearAllAtomBookmarks`, `clearAllBondBookmarks`, `expandAttachmentPoints`, explicit valence calculation, `ProcessMolProps`, conformer selection, `assignChiralTypesFromBondDirs`, `assignChiralTypesFrom3D`, `detectAtropisomerChirality`, and `clearSingleBondDirFlags`.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Add focused Rust tests for COSMolKit `finish_mol_processing` pre-sanitize ordering covering wedge H removal hazards, 3D chirality before sanitize, atom bookmark clearing equivalence where represented, bond bookmark clearing equivalence where represented, and `ProcessMolProps` ordering.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict finish_mol_processing -- --nocapture`.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Port RDKit `FileParserUtils::finishMolProcessing` `params.sanitize && params.removeHs` branch into COSMolKit `finish_mol_processing`, including `sanitizeMol(SANITIZE_CLEANUP)`, `detectBondStereochemistry`, `removeHs`, and no full `sanitizeMol` before `assignStereochemistry`.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Add focused Rust tests for RDKit `finishMolProcessing` `sanitize=true, removeHs=true` covering cleanup before stereo, double-bond stereogenic hydrogen, wedged hydrogen, explicit hydrogen on aromatic nitrogen, and final atom and bond stereo.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict molblock_sanitize_remove_hs -- --nocapture`.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Port RDKit `FileParserUtils::finishMolProcessing` `params.sanitize && !params.removeHs` branch into COSMolKit `finish_mol_processing`, including full `sanitizeMol`, `detectBondStereochemistry`, and `assignStereochemistry(*res, true, true, true)`.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Add focused Rust tests for RDKit `finishMolProcessing` `sanitize=true, removeHs=false` covering explicit hydrogen preservation, double-bond stereo detection after sanitize, final stereochemistry assignment, and property cache state.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict molblock_sanitize_keep_hs -- --nocapture`.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Port RDKit `FileParserUtils::finishMolProcessing` `!params.sanitize` branch into COSMolKit `finish_mol_processing`, including `detectBondStereochemistry(*res)` and excluding `sanitizeMol`, `removeHs`, and `assignStereochemistry`.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Add focused Rust tests for RDKit `finishMolProcessing` `sanitize=false` covering raw valence, aromaticity state, cache state, retained hydrogens, detected bond stereochemistry, no final CIP assignment, and delayed `sanitize()` parity.
Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict molblock_unsanitized -- --nocapture`.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Port RDKit `QueryOps::completeMolQueries` from `QueryOps.cpp` into COSMolKit `complete_mol_queries`.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Add focused Rust tests for RDKit `_NeedsQueryScan` and `QueryOps::completeMolQueries` covering V2000 atom queries, V2000 bond queries, V3000 atom queries, V3000 bond queries, ring-bond-count completion, and `_NeedsQueryScan` clearing.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict complete_mol_queries -- --nocapture`.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Port RDKit `ForwardSDMolSupplier::_next` normal-record behavior into COSMolKit `read_next_sdf_record`.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Add focused Rust tests for RDKit `ForwardSDMolSupplier::_next` normal-record behavior covering line counters, delimiter consumption, EOF state, and first-record molecule return.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict read_next_sdf_record -- --nocapture`.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Port RDKit `ForwardSDMolSupplier::_next` parse-exception and sanitize-exception recovery behavior into COSMolKit `read_next_sdf_record`.
Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Add focused Rust tests for RDKit `ForwardSDMolSupplier::_next` recovery behavior covering file-parse exception recovery, sanitize exception recovery, missing `$$$$`, malformed data field recovery, and next-record alignment.
Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict sdf_supplier_recovery -- --nocapture`.
Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Port RDKit `ForwardSDMolSupplier::readMolProps` behavior into COSMolKit `parse_sdf_data_fields` and `process_sdf_property_lists`.
Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Add focused Rust tests for RDKit `ForwardSDMolSupplier::readMolProps` covering named data fields, blank values, repeated fields, invalid headers, property lists, atom property lists, bond property lists, and `process_property_lists=false`.
Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict sdf_data_fields -- --nocapture`.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Port RDKit `SDMolSupplier` indexed-record parameter behavior into COSMolKit `SdfDataset::open_with_params`, `read_indexed_sdf_record`, and `SdfDataset::record_with_params`.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Add focused Rust tests for RDKit `SDMolSupplier` indexed-record behavior covering random access, out-of-range access, skipped invalid records, parameter reuse, metadata offsets, and sliced batch reads.
Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict sdf_dataset -- --nocapture`.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Modify Python `Molecule.read_mol` to pass `sanitize`, `remove_hs`, and `strict_parsing` through `SdfReadParams` instead of calling `reject_unsanitized_mol_reader`.
Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Modify Python `Molecule.read_mol_from_str` to pass `sanitize`, `remove_hs`, and `strict_parsing` through `SdfReadParams` instead of calling `reject_unsanitized_mol_reader`.
Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Modify Python `Molecule.read_sdf` to pass `sanitize`, `remove_hs`, and `strict_parsing` through `SdfReadParams` instead of calling `reject_unsanitized_mol_reader`.
Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Modify Python `Molecule.read_sdf_from_str` to pass `sanitize`, `remove_hs`, and `strict_parsing` through `SdfReadParams` instead of calling `reject_unsanitized_mol_reader`.
Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Modify Python `MoleculeBatch.read_sdf`, `MoleculeBatch.read_sdf_records_from_str`, `PySdfDataset::open`, and `PySdfReader::open` to propagate `sanitize`, `remove_hs`, and `strict_parsing` through `SdfReadParams`.
Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [x]: Run `.venv/bin/maturin develop --manifest-path python/Cargo.toml`.
Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [x]: Add Python tests in `python/tests/test_python_api_sanity.py` for `sanitize=False` on `Molecule.read_mol`, `Molecule.read_mol_from_str`, `Molecule.read_sdf`, and `Molecule.read_sdf_from_str`.
Step 152 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 153 [x]: Run `.venv/bin/pytest python/tests/test_python_api_sanity.py -k "read_mol or read_sdf or sanitize"`.
Step 154 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 155 [x]: Update `python/docs/source/io.rst` to describe RDKit-source-backed `sanitize`, `remove_hs`, `strict_parsing`, MolBlock behavior, SDF data-field behavior, and delayed operation behavior.
Step 156 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 157 [x]: Update `python/docs/source/molecule.rst` to describe RDKit-source-backed molecule IO finalization behavior and delayed `sanitize()` and `without_hydrogens()` behavior.
Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 159 [x]: Update Python examples to include one MolBlock `sanitize=False` delayed sanitize example and one SDF `remove_hs=False` delayed hydrogen-removal example.
Step 160 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 161 [x]: Run `cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen`.
Step 162 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 163 [x]: Run `.venv/bin/python -m sphinx -b html python/docs/source python/docs/build/html`.
Step 164 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 165 [x]: Run `.venv/bin/basedpyright python/tests python/examples`.
Step 166 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 167 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molfile_read_parity -- --nocapture`.
Step 168 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 169 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_sdf_read_parity -- --nocapture`.
Step 170 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 171 [x]: Run `cargo fmt --all`.
Step 172 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 173 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 174 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 175 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict`.
Step 176 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 177 [x]: Run `cargo test --workspace --features cosmolkit-core/op-contracts-strict`.
Step 178 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 179 [x]: Run `.venv/bin/pytest`.
Step 180 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 181 [x]: Update `dev/gap_reports/rdkit_molblock_reader_source_map.md` with final zero-gap findings for every mapped RDKit reader function in this plan.
Step 182 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 183 [x]: Update `dev/gap_reports/rdkit_molblock_finalization_source_map.md` with final zero-gap findings for every mapped RDKit finalization function in this plan.
Step 184 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 185 [x]: Update `dev/gap_reports/rdkit_sdf_supplier_source_map.md` with final zero-gap findings for every mapped RDKit supplier function in this plan.
