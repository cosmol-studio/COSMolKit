# RDKit MMFF Parity Completion Plan

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
- Every real task step must be immediately preceded by reading `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`.
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

## Parity Gates

- Every RDKit-successful, deterministically embedded force-field row in the active profile must execute initial-energy, initial-gradient, single-conformer, and multi-conformer comparisons.
- The normal `smiles_small` profile must compare every eligible row at the declared `1e-6` energy, gradient, and coordinate tolerances.
- The strict force-field audit must compare parameter availability, explicit-H behavior, MMFF atom types, and MMFF charges on every `smiles_5000` row. Initial energy, full gradient, and single-/multi-conformer optimizer parity remain exhaustive on every eligible `smiles_small` row and on every curated regression. This is a declared two-tier boundary, not a hidden skip or reduced optimizer subset presented as 5,000-row optimizer coverage.
- ChEMBL audit failures must become committed curated corpus rows or smallest-boundary regression fixtures before their production fix is marked complete.
- RDKit-success/COSMolKit-unsupported is a parity failure, not an accepted optimizer result.
- Embedding failures and RDKit parameter unavailability must be recorded and asserted explicitly rather than silently skipped.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Modify the force-field parity integration tests so every eligible active-profile row executes all declared comparison surfaces and explicit-hydrogen parameter availability is compared.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Run the focused force-field parity integration tests.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 52 [x]: Regenerate the `smiles_small` force-field expected-data family through `tools/testdata/rdkit/generate_all.py`.
Step 53 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 54 [x]: Run the focused force-field parity integration tests against the regenerated expected data.
Step 55 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 56 [x]: Audit the complete focused-test failures and write `dev/gap_reports/rdkit_mmff_parity_failures.md`.
Step 57 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 58 [x]: Fix force-field parity assertions to model RDKit no-force-field sentinels and unchanged coordinates without skipping any corpus row.
Step 59 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 60 [x]: Run the focused force-field parity integration tests after the sentinel assertion fix.
Step 61 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Instrument pinned RDKit and COSMolKit bond-stretch parameter lookup boundaries and write the first-divergence probe report to `dev/gap_reports/rdkit_mmff_bond_empirical_probe.md`.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Port the complete RDKit MMFF bond-stretch empirical-rule behavior into the corresponding Rust source block.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Add smallest-boundary bond empirical-rule tests and public optimizer regressions for every audited bond-triggering case.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Run the focused bond empirical-rule and force-field parity tests.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Instrument pinned RDKit and COSMolKit angle-bend parameter lookup boundaries and write the first-divergence probe report to `dev/gap_reports/rdkit_mmff_angle_empirical_probe.md`.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Port the complete RDKit MMFF angle-bend empirical-rule behavior into the corresponding Rust source block.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Add smallest-boundary angle empirical-rule tests and public optimizer regressions for every audited angle-triggering case.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run the focused angle empirical-rule and force-field parity tests.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Instrument pinned RDKit and COSMolKit torsion parameter lookup boundaries and write the first-divergence probe report to `dev/gap_reports/rdkit_mmff_torsion_empirical_probe.md`.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Port the complete RDKit MMFF torsion empirical-rule behavior into the corresponding Rust source block.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Add smallest-boundary torsion empirical-rule tests and public optimizer regressions for every audited torsion-triggering case.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Run the focused torsion empirical-rule and force-field parity tests.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Instrument pinned RDKit and COSMolKit explicit-hydrogen MMFF typing boundaries and write the first-divergence probe report to `dev/gap_reports/rdkit_mmff_hydrogen_typing_probe.md`.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Port the complete RDKit MMFF hydrogen atom-typing behavior for every modeled neighbor element and neighbor atom type.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Add smallest-boundary hydrogen-typing tests and explicit-hydrogen parameter-availability regressions for every audited failure class.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Run the focused hydrogen-typing and force-field parity tests.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Instrument pinned RDKit and COSMolKit force-field term construction, gradient, minimizer iteration, and termination boundaries and write the first-divergence probe report to `dev/gap_reports/rdkit_mmff_optimizer_probe.md`.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Fix the complete source-backed MMFF optimizer behavior identified by the first-divergence report.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Add smallest-boundary optimizer-state tests and public coordinate regressions for every audited optimizer failure class.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Instrument the remaining MMFF minimizer divergence after exact initial energy and gradient parity and write `dev/gap_reports/rdkit_mmff_minimizer_second_divergence.md`.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Fix the complete source-backed minimizer or contribution-gradient behavior identified by the second-divergence report.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Run the focused second-divergence regressions and record the remaining row-103 public optimizer mismatch without weakening its assertion.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Instrument pinned RDKit and COSMolKit after the contribution-gradient evaluation-order fix, locate the next bit-level optimizer divergence, and write `dev/gap_reports/rdkit_mmff_minimizer_third_divergence.md`.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Fix the complete source-backed MMFF behavior identified by the third-divergence report.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Add smallest-boundary regressions for every behavior fixed from the third-divergence report.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Run the focused third-divergence smallest-boundary and row-103 regressions, then record the remaining row-123 public optimizer mismatch without weakening its assertion.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Instrument pinned RDKit and COSMolKit after the torsion match-order fix, locate the next bit-level optimizer divergence, and write `dev/gap_reports/rdkit_mmff_minimizer_fourth_divergence.md`.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Fix the complete source-backed MMFF behavior identified by the fourth-divergence report.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Add smallest-boundary regressions for every behavior fixed from the fourth-divergence report.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Run the focused fourth-divergence regressions, charge and 1-4 regressions, angle regressions, and complete MMFF optimizer tests.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Run the focused MMFF optimizer and complete `smiles_small` force-field parity suites, then record the remaining UFF row-34 multi-conformer and row-47 single-conformer mismatches without weakening their assertions.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Instrument pinned RDKit and COSMolKit at the UFF construction, contribution, gradient, and minimizer boundaries for rows 34 and 47, locate each first divergence, and write `dev/gap_reports/rdkit_uff_smiles_small_divergences.md`.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Fix the complete source-backed UFF behaviors identified by the row-34 and row-47 divergence report.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Add smallest-boundary UFF regressions for every behavior fixed from the row-34 and row-47 divergence report.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Run the focused UFF regressions and complete `smiles_small` force-field parity suite.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Instrument pinned RDKit and COSMolKit after the UFF central-bond match-order fix, compare complete ordered torsion contribution tuples for rows 34 and 47, and extend `dev/gap_reports/rdkit_uff_smiles_small_divergences.md` with the next source-level divergence.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Fix the complete source-backed UFF behavior identified by the full torsion-contribution order probe.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Extend the smallest-boundary and public optimizer regressions for every additional UFF behavior fixed by the full torsion-contribution order probe.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Run the focused UFF regressions after the full torsion-contribution order fix and record the remaining row-47 optimizer mismatch without weakening its assertion.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Instrument pinned RDKit and COSMolKit across every row-47 UFF contribution and minimizer state, locate the next bit-level source divergence, and extend `dev/gap_reports/rdkit_uff_smiles_small_divergences.md`.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Isolate every force-field golden surface from RDKit molecule mutation and regenerate the `smiles_small` force-field expected data.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Add a smallest-boundary golden-generator state-isolation regression while retaining the complete row-47 public UFF regression.
Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Run the focused row-47 UFF regression and complete `smiles_small` force-field parity suite, recording the remaining row-81 multi-conformer and row-113 single-conformer coordinate mismatches without weakening their assertions.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Extend the complete public UFF optimizer regression to rows 81 and 113 so both remaining failures are independently executable.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Run the focused rows-81-and-113 UFF regression and record both unchanged `1e-6` failures.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Instrument pinned RDKit and COSMolKit across complete UFF contribution and minimizer-state boundaries for rows 81 and 113, locate the first source-level divergence, and extend `dev/gap_reports/rdkit_uff_smiles_small_divergences.md`.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Fix the complete source-backed UFF behavior identified by the rows-81-and-113 probe.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Add smallest-boundary regressions for every UFF behavior fixed by the rows-81-and-113 probe.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Run the focused rows-81-and-113 regressions and complete `smiles_small` force-field parity suite.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Run the lightweight strict force-field audit: all 5,000 rows for parameter availability, explicit-H behavior, MMFF atom types, and charges, plus the complete `smiles_small` energy, gradient, and single-/multi-conformer optimizer suite; record any mismatch as an executable regression before continuing.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Run the deterministic ChEMBL MMFF audit and add every new failure class to the committed curated corpus.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Run formatting, strict core checks, strict release core tests, and strict release workspace tests.
