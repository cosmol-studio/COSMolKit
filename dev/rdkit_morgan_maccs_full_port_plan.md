# RDKit Morgan And MACCS Full Source Port Plan

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
- Never treat a required reading step as "already read".
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
- Do not use "smallest subpart", skeleton code, dispatch-only code, placeholder code, TODO-only code, or partial porting as completion.
- If a step is too large to complete fully, regenerate the plan by splitting that step into smaller full-port steps.
- Each split step must still represent a complete source-backed behavior, not a placeholder.
- A `Port` step is complete only when the selected RDKit behavior is fully implemented, source-backed, tested, and no placeholder remains.
- During porting, do not execute any git operation unless the user explicitly asks for it.
- This includes read-only git operations such as `git status`, `git diff`, and `git log`.

## Scope And Non-Negotiable Standard

- Scope is only Morgan fingerprints and MACCS keys exposed by COSMolKit.
- The target oracle is pinned RDKit source under `third_party/rdkit/Code/GraphMol/Fingerprints/`.
- Morgan must reproduce the RDKit `FingerprintGenerator` pipeline, `MorganGenerator`, `MorganFingerprints` wrappers, invariants, environment generation, count simulation, sparse/count/bit-vector behavior, and additional outputs.
- MACCS must reproduce RDKit `MACCS.cpp` `Patterns`, `GenerateFP`, and `MACCSFingerprints::getFingerprintAsBitVect`.
- SMARTS and substructure behavior required by Morgan feature invariants or MACCS keys must be source-backed instead of replaced with local element/degree/ring heuristics.
- Exact bit-vector parity is required on every covered corpus row and branch.
- No pass percentage, similarity correlation, compatible hashing, structural approximation, or partial branch support may be described as parity.
- Unsupported behavior must fail explicitly until ported; it must not emit chemically meaningful-looking fingerprint bits.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit RDKit fingerprint source functions and write `dev/gap_reports/rdkit_morgan_maccs_source_inventory.md` covering `FingerprintGenerator.{h,cpp}`, `FingerprintUtil.{h,cpp}`, `MorganGenerator.{h,cpp}`, `MorganFingerprints.{h,cpp}`, `MACCS.{h,cpp}`, and every currently exposed COSMolKit fingerprint API.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Audit COSMolKit SMARTS and substructure support needed by RDKit `MorganFeatureAtomInvGenerator::getAtomInvariants` and `MACCS.cpp::Patterns` and write `dev/gap_reports/rdkit_fingerprint_smarts_substruct_gap_report.md`.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Port RDKit `AdditionalOutput` allocation and storage semantics from `FingerprintGenerator.h::AdditionalOutput` into COSMolKit without placeholder fields.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Add unit tests for COSMolKit `AdditionalOutput` allocation, empty-state, atom-to-bits, bit-info-map, atom-counts, atoms-per-bit, and reinitialization behavior.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict fingerprint additional_output -- --nocapture`.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Port RDKit `FingerprintArguments::FingerprintArguments`, `FingerprintArguments::commonArgumentsString`, `FingerprintArguments::toJSON`, and `FingerprintArguments::fromJSON` into COSMolKit.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Add unit tests for `FingerprintArguments` defaults, count bounds, fingerprint size, number of bits per feature, chirality flag, and JSON round-trip behavior.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict fingerprint_arguments -- --nocapture`.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Port RDKit `MorganArguments::MorganArguments`, `MorganArguments::infoString`, `MorganArguments::toJSON`, and `MorganArguments::fromJSON` into COSMolKit.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Add unit tests for `MorganArguments` radius, count simulation, chirality, only-nonzero invariants, count bounds, fp size, redundant-environment flag, bond-type flag, and JSON round-trip behavior.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict morgan_arguments -- --nocapture`.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Port RDKit `MorganAtomInvGenerator::MorganAtomInvGenerator`, `MorganAtomInvGenerator::getAtomInvariants`, `MorganAtomInvGenerator::infoString`, `MorganAtomInvGenerator::toJSON`, `MorganAtomInvGenerator::fromJSON`, and `MorganAtomInvGenerator::clone`.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Add RDKit golden tests for `MorganAtomInvGenerator::getAtomInvariants` covering atomic number, total degree, total hydrogens, formal charge, isotope delta mass, and ring-membership branches.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict morgan_atom_invariants -- --nocapture`.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Port RDKit `FingerprintUtil.cpp::ss_matcher::ss_matcher`, `FingerprintUtil.cpp::ss_matcher::ss_matcher(const std::string &)`, and `FingerprintUtil.cpp::ss_matcher::getMatcher` for Morgan feature invariant pattern handling.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Port RDKit `FingerprintUtil.cpp::defaultFeatureSmarts` into a source-backed COSMolKit feature-pattern table without manual donor/acceptor classification.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Port RDKit SMARTS parser behavior required by `defaultFeatureSmarts` into COSMolKit `parse_smarts` and query-node construction.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Add targeted SMARTS parser tests for every RDKit `defaultFeatureSmarts` pattern string.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict smarts default_feature_smarts -- --nocapture`.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Port RDKit `SubstructMatch` behavior required by `FingerprintUtil.cpp::getFeatureInvariants`, including recursive SMARTS, atom OR, atom negation, atom hydrogen-count, aromaticity, ring-membership, valence, and bond query semantics used by `defaultFeatureSmarts`.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Add RDKit parity tests for feature SMARTS substructure matches and matched atom indices for every `defaultFeatureSmarts` pattern.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict substruct feature_smarts -- --nocapture`.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Port RDKit `FingerprintUtil.cpp::getFeatureInvariants` into COSMolKit using parsed SMARTS matchers and source-backed substructure matching only.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Port RDKit `MorganFeatureAtomInvGenerator::MorganFeatureAtomInvGenerator`, `MorganFeatureAtomInvGenerator::~MorganFeatureAtomInvGenerator`, `MorganFeatureAtomInvGenerator::cleanUpPatterns`, `MorganFeatureAtomInvGenerator::getAtomInvariants`, `MorganFeatureAtomInvGenerator::infoString`, `MorganFeatureAtomInvGenerator::toJSON`, `MorganFeatureAtomInvGenerator::fromJSON`, and `MorganFeatureAtomInvGenerator::clone`.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Add RDKit golden tests for `MorganFeatureAtomInvGenerator::getAtomInvariants` covering default patterns, supplied patterns, overlapping matches, no-match molecules, and invalid pattern failure behavior.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict morgan_feature_invariants -- --nocapture`.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Port RDKit `MorganBondInvGenerator::MorganBondInvGenerator`, `MorganBondInvGenerator::getBondInvariants`, `MorganBondInvGenerator::infoString`, `MorganBondInvGenerator::toJSON`, `MorganBondInvGenerator::fromJSON`, and `MorganBondInvGenerator::clone`.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Add RDKit golden tests for `MorganBondInvGenerator::getBondInvariants` covering all RDKit bond types, aromatic bonds, dative bonds, no-bond-type mode, and bond chirality mode.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict morgan_bond_invariants -- --nocapture`.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Port RDKit `MorganAtomEnv<OutputType>::getBitId` and `MorganAtomEnv<OutputType>::updateAdditionalOutput`.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Add unit tests for `MorganAtomEnv` bit-id identity behavior and additional-output updates for atom-to-bits, bit-info-map, atom-counts, and atoms-per-bit.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict morgan_atom_env -- --nocapture`.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Port RDKit `MorganEnvGenerator<OutputType>::getEnvironments` exactly, including atom order, dead atom handling, ignore atoms, from atoms, neighborhood bond sets, redundant environment handling, chirality hashing from CIP properties, duplicate environment suppression, and layer semantics.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Add RDKit golden tests for `MorganEnvGenerator<OutputType>::getEnvironments` covering radius 0, radius 1, radius 2, isolated atoms, duplicate neighborhoods, redundant-environment true and false, fromAtoms, ignoreAtoms, chirality, and bond-type branches.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict morgan_env_generator -- --nocapture`.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Port RDKit `MorganEnvGenerator<OutputType>::infoString`, `MorganEnvGenerator<OutputType>::toJSON`, `MorganEnvGenerator<OutputType>::fromJSON`, and `MorganEnvGenerator<OutputType>::getResultSize`.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Add unit tests for `MorganEnvGenerator` info string, JSON behavior, and result-size behavior.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict morgan_env_generator_metadata -- --nocapture`.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Port RDKit `MorganFingerprint::getMorganGenerator<OutputType>` ownership, default generator, custom atom-invariant generator, custom bond-invariant generator, count simulation, chirality, only-nonzero, redundant-environment, and bond-type wiring.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Add unit tests for `getMorganGenerator` option wiring and generator ownership behavior.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict get_morgan_generator -- --nocapture`.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Port RDKit `FingerprintGenerator<OutputType>::getFingerprintHelper`, including input invariant validation, generator invariant selection, environment generation, bit-id extraction, hash-results handling, count accumulation, and additional-output update order.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Add RDKit golden tests for `FingerprintGenerator<OutputType>::getFingerprintHelper` on Morgan environments covering custom atom invariants, custom bond invariants, fromAtoms, ignoreAtoms, and additional output.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict fingerprint_generator_helper morgan -- --nocapture`. Note: Cargo rejected the two-filter form as invalid syntax; the equivalent single-filter command `cargo test -p cosmolkit-core --features op-contracts-strict fingerprint_generator_helper_morgan -- --nocapture` was run and passed 5 tests.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Port RDKit `FingerprintGenerator.cpp::duplicateAdditionalOutputBit` and `FingerprintGenerator.cpp::setupTempAdditionalOutput`.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Add unit tests for additional-output bit duplication and temporary additional-output setup during count simulation.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict fingerprint_count_simulation_additional_output -- --nocapture`.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Port RDKit `FingerprintGenerator<OutputType>::getSparseCountFingerprint`, `FingerprintGenerator<OutputType>::getSparseFingerprint`, `FingerprintGenerator<OutputType>::getCountFingerprint`, and `FingerprintGenerator<OutputType>::getFingerprint`.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Add RDKit golden tests for Morgan sparse-count, sparse-bit, hashed-count, and explicit-bit outputs across radius, fp size, count simulation, and additional output.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict morgan_fingerprint_generator_outputs -- --nocapture`.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Port RDKit `MorganFingerprints.cpp::getFingerprint`.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Port RDKit `MorganFingerprints.cpp::getHashedFingerprint`.
Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Port RDKit `MorganFingerprints.cpp::getFingerprintAsBitVect`.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Add RDKit parity tests for public Morgan API wrappers covering sparse count, hashed count, explicit bit vector, atomsSettingBits, custom invariants, fromAtoms, chirality, and count simulation.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_morgan_fingerprint_parity -- --nocapture`.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Port RDKit SMARTS parser behavior required by every `MACCS.cpp::Patterns` SMARTS string into COSMolKit `parse_smarts` and query-node construction.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Add targeted parser tests for every `MACCS.cpp::Patterns` SMARTS string, including recursive SMARTS, ring query atoms, ring bonds, non-ring bonds, wildcard atoms, negation, OR alternatives, hydrogen-count queries, and branch/ring closure syntax.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict smarts maccs_patterns -- --nocapture`. Note: Cargo rejected the two-filter form as invalid syntax; the equivalent filter commands `cargo test -p cosmolkit-core --features op-contracts-strict maccs_patterns -- --nocapture` and `cargo test -p cosmolkit-core --features op-contracts-strict smarts -- --nocapture` were run and passed.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Port RDKit `SubstructMatch` behavior required by every `MACCS.cpp::Patterns` SMARTS string.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Add RDKit parity tests for MACCS pattern substructure match truth values and first-match atom maps for every `MACCS.cpp::Patterns` SMARTS string.
Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict substruct maccs_patterns -- --nocapture`. Note: Cargo rejected the two-filter form as invalid syntax; the equivalent single-filter commands `cargo test -p cosmolkit-core --features op-contracts-strict maccs_patterns -- --nocapture` and `cargo test -p cosmolkit-core --features op-contracts-strict substruct -- --nocapture` were run and passed 4 and 20 tests respectively.
Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Port RDKit `MACCS.cpp::Patterns` initialization into a source-backed COSMolKit pattern table preserving RDKit bit numbering and bit-zero-unused semantics.
Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Add tests proving the COSMolKit MACCS pattern table contains every RDKit `Patterns` field with the exact SMARTS string and bit number.
Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict maccs pattern_table -- --nocapture`. Note: Cargo rejected the two-filter form as invalid syntax; the equivalent single-filter commands `cargo test -p cosmolkit-core --features op-contracts-strict pattern_table -- --nocapture` and `cargo test -p cosmolkit-core --features op-contracts-strict maccs -- --nocapture` were run and passed 1 and 5 tests respectively.
Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Port the element and atom-property direct-key blocks from RDKit `MACCS.cpp::GenerateFP` for keys 1 through 40.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Add RDKit parity tests for MACCS keys 1 through 40 using targeted molecules that independently trigger and independently avoid each key.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict maccs_keys_001_040 -- --nocapture`.
Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Port RDKit `MACCS.cpp::GenerateFP` keys 41 through 80.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Add RDKit parity tests for MACCS keys 41 through 80 using targeted molecules that independently trigger and independently avoid each key.
Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict maccs_keys_041_080 -- --nocapture`.
Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Port RDKit `MACCS.cpp::GenerateFP` keys 81 through 120.
Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Add RDKit parity tests for MACCS keys 81 through 120 using targeted molecules that independently trigger and independently avoid each key.
Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict maccs_keys_081_120 -- --nocapture`.
Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [x]: Port RDKit `MACCS.cpp::GenerateFP` keys 121 through 166.
Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [x]: Add RDKit parity tests for MACCS keys 121 through 166 using targeted molecules that independently trigger and independently avoid each key.
Step 152 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 153 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict maccs_keys_121_166 -- --nocapture`.
Step 154 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 155 [x]: Port RDKit `MACCSFingerprints::getFingerprintAsBitVect` including 167-bit internal vector semantics, bit 0 unused behavior, and COSMolKit public 166-bit projection.
Step 156 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 157 [x]: Add RDKit parity tests for full MACCS bit vectors covering internal RDKit bit numbering, public COSMolKit bit projection, empty molecules, single atoms, salts, aromatic molecules, heterocycles, charged molecules, isotopes, and all targeted key fixtures.
Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 159 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict maccs_fingerprint -- --nocapture`.
Step 160 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 161 [x]: Update `tests/scripts/gen_rdkit_morgan_fingerprint_golden.py` to generate RDKit sparse-count, sparse-bit, hashed-count, explicit-bit, additional-output, atom-invariant, bond-invariant, fromAtoms, ignoreAtoms, count-simulation, and chirality baselines from the pinned RDKit API.
Step 162 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 163 [x]: Add or update the Rust Morgan golden parity test to assert exact RDKit equality for every generated Morgan output field without percentages or compatibility thresholds.
Step 164 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 165 [x]: Write `dev/gap_reports/rdkit_ciplabeler_for_morgan_gap_report.md` documenting why the Step 165 Morgan parity failure requires RDKit `CIPLabeler::assignCIPLabels` rather than legacy COSMolKit CIP shortcuts, with exact source files and failing golden row.
Step 166 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 167 [x]: Port RDKit `CIPLabeler/Descriptor.h::Descriptor`, `CIPLabeler/Descriptor.h::to_string`, and descriptor ordering constants into a private COSMolKit CIPLabeler module.
Step 168 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 169 [x]: Add unit tests for RDKit `Descriptor` string mapping and descriptor ordering used by CIPLabeler sequence rules.
Step 170 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 171 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_descriptor -- --nocapture`.
Step 172 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 173 [x]: Port RDKit `CIPMol.cpp::CIPMol`, `CIPMol.cpp::getFractionalAtomicNum`, `CIPMol.cpp::getNumAtoms`, `CIPMol.cpp::getNumBonds`, `CIPMol.cpp::getAtom`, `CIPMol.cpp::getBond`, `CIPMol.cpp::getBonds`, `CIPMol.cpp::getNeighbors`, `CIPMol.cpp::isInRing`, and `CIPMol.cpp::getBondOrder` into the COSMolKit CIPLabeler module.
Step 174 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 175 [x]: Add RDKit golden tests for `CIPMol` fractional atomic numbers, kekulized bond orders, ring-bond checks, atom access, bond access, bond iteration, and neighbor iteration.
Step 176 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 177 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_cipmol -- --nocapture`.
Step 178 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 179 [x]: Port RDKit `Node.{h,cpp}` functions `Node::Node`, `Node::duplicate`, `Node::newChild`, `Node::getEdges`, `Node::getEdges(Atom *)`, `Node::getDistance`, `Node::getAtom`, `Node::getDigraph`, `Node::isDuplicate`, `Node::setAux`, `Node::getAux`, `Node::set`, `Node::isSet`, and `Node::getFlags` into the COSMolKit CIPLabeler module.
Step 180 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 181 [x]: Port RDKit `Edge.{h,cpp}` functions `Edge::Edge`, `Edge::getBeg`, `Edge::getEnd`, `Edge::getOther`, `Edge::isBeg`, and `Edge::isEnd` into the COSMolKit CIPLabeler module.
Step 182 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 183 [x]: Add unit tests for RDKit `Node` and `Edge` duplicate-node flags, aux descriptor storage, child creation, edge access, and opposite-end lookup.
Step 184 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 185 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_node_edge -- --nocapture`. Note: the concrete tests are named with separate `ciplabeler_node` and `ciplabeler_edge` filters; both equivalent focused commands were run and passed.
Step 186 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 187 [x]: Port RDKit `Digraph.{h,cpp}` functions `Digraph::Digraph`, `Digraph::init`, `Digraph::expand`, `Digraph::getOriginalRoot`, `Digraph::getCurrentRoot`, `Digraph::changeRoot`, `Digraph::getNodes`, `Digraph::setRule6Ref`, `Digraph::getRule6Ref`, and `Digraph::getMol` into the COSMolKit CIPLabeler module.
Step 188 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 189 [x]: Add RDKit golden tests for `Digraph` root construction, recursive expansion, duplicate-node creation, root changes, focus-node lookup, and Rule6 reference storage.
Step 190 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 191 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_digraph -- --nocapture`.
Step 192 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 193 [x]: Port RDKit `Priority.h::Priority`, `Sort.cpp::Sort::prioritize`, `Sort.cpp::Sort::getGroups`, and `Sort.cpp::Sort::getSorter` into the COSMolKit CIPLabeler module.
Step 194 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 195 [x]: Add unit tests for CIPLabeler priority uniqueness, pseudoasymmetric priority, edge sorting order, and partition grouping.
Step 196 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 197 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_sort_priority -- --nocapture`.
Step 198 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 199 [x]: Port RDKit `rules/SequenceRule.{h,cpp}` functions `SequenceRule::compare`, `SequenceRule::getSorter`, `SequenceRule::recursiveCompare`, and `SequenceRule::sort` into the COSMolKit CIPLabeler module.
Step 200 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 201 [x]: Port RDKit `rules/Rules.h::Rules` constructor, `Rules::compare`, `Rules::getNumSubRules`, `Rules::getSorter`, and `Rules::sort` into the COSMolKit CIPLabeler module.
Step 202 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 203 [x]: Port RDKit `rules/Rule1a.cpp::compare`, `Rule1b.cpp::compare`, `Rule2.cpp::compare`, and `Rule3.cpp::compare` into the COSMolKit CIPLabeler module.
Step 204 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 205 [x]: Port RDKit `rules/Pairlist.h` pair-list construction, `PairList::ref`, `PairList::getRefDescriptor`, `PairList::add`, `PairList::addAll`, `PairList::getPairing`, `PairList::compareTo`, `PairList::operator<`, `PairList::toString`, and `PairList::addAndPair` into the COSMolKit CIPLabeler module.
Step 206 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 207 [x]: Port RDKit `rules/Rule4a.cpp::compare`, `Rule4b.cpp` functions `Rule4b::compare`, `Rule4b::getReferenceDescriptors`, `Rule4b::hasDescriptors`, `Rule4b::getReference`, `Rule4b::initialLevel`, `Rule4b::getNextLevel`, `Rule4b::toNodeList`, `Rule4b::newPairLists`, `Rule4b::fillPairs`, `Rule4b::comparePairs`, `Rule4b::getRefSorter`, and `Rule4c.cpp::compare` into the COSMolKit CIPLabeler module with focused Rule4 tests.
Step 208 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 209 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_rules_rule4 -- --nocapture`.
Step 210 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 211 [x]: Port RDKit `rules/Rule5New.cpp::compare` and `Rule6.cpp::compare` into the COSMolKit CIPLabeler module with focused Rule5New and Rule6 tests.
Step 212 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 213 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_rules_rule5new_rule6 -- --nocapture`.
Step 214 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 215 [x]: Port RDKit `configs/Configuration.{h,cpp}` functions `Configuration::Configuration`, `Configuration::getFocus`, `Configuration::getFoci`, `Configuration::setCarriers`, `Configuration::getCarriers`, `Configuration::getDigraph`, `Configuration::parity4`, and `Configuration::label(Node *, Digraph &, const Rules &)` into the COSMolKit CIPLabeler module with focused Configuration tests.
Step 216 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 217 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_configuration -- --nocapture`.
Step 218 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 219 [x]: Port RDKit `configs/Tetrahedral.cpp` functions `Tetrahedral::Tetrahedral`, `Tetrahedral::setPrimaryLabel`, `Tetrahedral::hasPrimaryLabel`, `Tetrahedral::resetPrimaryLabel`, `Tetrahedral::label(const Rules &)`, `Tetrahedral::label(Node *, Digraph &, const Rules &)`, and `Tetrahedral::label(Node *, const Rules &)` into the COSMolKit CIPLabeler module with focused Tetrahedral tests.
Step 220 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 221 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_tetrahedral -- --nocapture`.
Step 222 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 223 [x]: Port RDKit `configs/Sp2Bond.cpp` functions `Sp2Bond::Sp2Bond`, `Sp2Bond::setPrimaryLabel`, `Sp2Bond::hasPrimaryLabel`, `Sp2Bond::resetPrimaryLabel`, and all `Sp2Bond::label` overloads into the COSMolKit CIPLabeler module with focused Sp2Bond tests.
Step 224 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 225 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_sp2bond -- --nocapture`.
Step 226 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 227 [x]: Port RDKit `configs/AtropisomerBond.cpp` functions `AtropisomerBond::AtropisomerBond`, `AtropisomerBond::setPrimaryLabel`, `AtropisomerBond::hasPrimaryLabel`, `AtropisomerBond::resetPrimaryLabel`, and all `AtropisomerBond::label` overloads into the COSMolKit CIPLabeler module with focused AtropisomerBond tests.
Step 228 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 229 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_atropisomerbond -- --nocapture`.
Step 230 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 231 [x]: Add RDKit golden tests for Configuration parity, Tetrahedral R/S/r/s labeling, Sp2Bond E/Z labeling, AtropisomerBond M/P labeling, primary-label reset, and CIP neighbor-order writeback.
Step 232 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 233 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_configs -- --nocapture`.
Step 234 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 235 [x]: Port RDKit `Chirality.cpp::findHighestCIPNeighbor` into `CipSp2Bond::new` so RDKit E/Z bonds without explicit stereo atoms construct source-backed carriers instead of returning unsupported.
Step 236 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 237 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_sp2bond -- --nocapture`.
Step 238 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 239 [x]: Port RDKit `CIPLabeler.cpp::findConfigs`, `CIPLabeler.cpp::labelAux`, `CIPLabeler.cpp::label`, `CIPLabeler.cpp::assignCIPLabels(ROMol &, const boost::dynamic_bitset<> &, const boost::dynamic_bitset<> &, unsigned int)`, and `CIPLabeler.cpp::assignCIPLabels(ROMol &, unsigned int)` into the COSMolKit CIPLabeler module.
Step 240 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 241 [x]: Add RDKit golden tests for `assignCIPLabels` on explicit tetrahedral centers, implicit-hydrogen tetrahedral centers, pseudoasymmetric centers, E/Z double bonds, atropisomer bonds, molecule-level `_CIPComputed`, subset atom/bond bitsets, and recursive-iteration limit behavior.
Step 242 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 243 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ciplabeler_assign -- --nocapture`.
Step 244 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 245 [x]: Modify Morgan `FingerprintGenerator<OutputType>::getFingerprintHelper`, `MorganEnvGenerator<OutputType>::getEnvironments`, and `MorganBondInvGenerator::getBondInvariants` to call the ported RDKit CIPLabeler branch instead of returning structured unsupported for missing `_CIPCode` or `_CIPComputed`.
Step 246 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 247 [x]: Add Morgan regression tests proving includeChirality recomputes missing atom and bond `_CIPCode` through the ported CIPLabeler path and preserves exact sparse-count, sparse-bit, hashed-count, explicit-bit, and additional-output golden fields for the former failing row.
Step 248 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 249 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_morgan_fingerprint_parity -- --nocapture`.
Step 250 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 251 [x]: Add `tests/scripts/gen_rdkit_maccs_fingerprint_golden.py` to generate RDKit MACCS 167-bit raw vectors, public 166-bit projected vectors, and per-key targeted fixture expectations from `MACCSkeys.GenMACCSKeys`.
Step 252 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 253 [x]: Add `crates/cosmolkit-core/tests/rdkit_maccs_fingerprint_parity.rs` to assert exact RDKit MACCS bit-vector parity on targeted fixtures, the small daily corpus, and the strict 5000-row corpus when explicitly selected.
Step 254 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 255 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_maccs_fingerprint_parity -- --nocapture`.
Step 256 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 257 [x]: Modify public Rust and Python fingerprint APIs so Morgan and MACCS either return exact source-ported results or return structured unsupported errors for any option not reproduced from RDKit.
Step 258 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 259 [x]: Add Python parity tests for exposed Morgan and MACCS APIs that compare exact RDKit bit vectors and exact additional-output fields where exposed.
Step 260 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 261 [x]: Run `.venv/bin/pytest python/tests/test_fingerprint_rdkit_bit_parity.py -q`.
Step 262 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 263 [x]: Remove or rewrite every Morgan and MACCS heuristic helper, approximate comment, unfinished marker, and compatibility claim in `crates/cosmolkit-core/src/properties/fingerprint.rs`, `crates/cosmolkit-core/src/support.rs`, `tests/README.md`, and `tests/parity_scope.md` according to the actual source-backed status.
Step 264 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 265 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict fingerprint -- --nocapture`.
Step 266 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 267 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_morgan_fingerprint_parity --test rdkit_maccs_fingerprint_parity -- --nocapture`.
Step 268 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 269 [x]: Run `COSMOLKIT_PARITY_PROFILE=smiles_5000 cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_morgan_fingerprint_parity --test rdkit_maccs_fingerprint_parity -- --nocapture`.
Step 270 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 271 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict`.
Step 272 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 273 [x]: Run `cargo fmt --all`.
Step 274 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 275 [x]: Write `dev/gap_reports/rdkit_morgan_maccs_full_port_validation.md` recording source marker closure, exact test commands, exact RDKit corpus/golden conditions, remaining unsupported branches if any, and why no compatibility threshold or subset claim remains.
