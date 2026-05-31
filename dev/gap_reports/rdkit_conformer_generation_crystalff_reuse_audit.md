# RDKit Conformer Generation CrystalFF Reuse Audit

## Scope

- Plan step: `dev/rdkit_conformer_generation_full_port_plan.md` Step 75.
- RDKit sources audited:
  - `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.cpp`
  - `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleContribs.cpp`
  - `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleM6.cpp`
- Rust modules audited:
  - `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/torsion_preferences.rs`
  - `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/torsion_angle_contribs.rs`
  - `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/torsion_angle_m6.rs`

## Reusable Implemented Surface

- `TorsionAngleContribs` already has source-backed construction, `add_contrib`,
  energy, and gradient paths corresponding to
  `CrystalFF::TorsionAngleContribs`.
- `TorsionAngleContribM6` already has source-backed construction, energy, and
  gradient paths corresponding to `CrystalFF::TorsionAngleContribM6`.
- `calc_torsion_energy` already implements the six-term M6 energy polynomial
  used by both RDKit `calcTorsionEnergyM6` implementations.
- `get_experimental_torsions` and `get_experimental_torsions_without_bonds`
  already model the RDKit overload shape and fill `CrystalFFDetails` fields
  used by the DistGeom force-field constructors.
- CrystalFF parameter data is already vendored through generated constants and
  parsed into `ExpTorsionAngleCollection`.

## Required Follow-Up Gaps

- `calc_torsion_energy` is named generically and carries `RDKit✔️❌` markers
  while the RDKit source symbol is `calcTorsionEnergyM6`; Step 77 should either
  add a source-named `calc_torsion_energy_m6` wrapper or rename/re-export the
  helper, then upgrade markers only after local complexity review.
- `TorsionAngleContribs::{new, add_contrib, get_energy, get_grad}` still carry
  `RDKit✔️❌` markers even though the behavior appears source-shaped. Step 81
  must inspect allocation, loop structure, indexing, owner-pointer refresh, and
  gradient math before upgrading markers or fixing real gaps.
- `TorsionAngleContribM6::{new, get_energy, get_grad}` still carry
  `RDKit✔️❌` markers. Step 85 must perform the same local review and targeted
  tests for the single-contribution path.
- `get_experimental_torsions` main overload is marked `RDKit✔️❗`; Step 89 must
  close the remaining proof gap against bridged-ring exclusion, constrained
  atoms, SMARTS matching, ETversion selection, small-ring torsions, macrocycle
  torsions, improper atom generation, flat-ring torsions, and basic-knowledge
  behavior.
- Verbose logging in `getExperimentalTorsions` is intentionally not observable
  in the Rust public path today; Step 89 must either prove it is not part of the
  modeled state space or document it as unsupported without changing chemistry.

## Reuse Decision

- Reuse the existing CrystalFF modules as the implementation base.
- Do not introduce a parallel CrystalFF implementation.
- Do not add RDKit interop.
- Do not replace the current code with heuristic torsion handling.
- The next steps should be marker upgrades plus targeted source-backed tests
  where parity is already implemented, and source-backed fixes only where
  tests or local review reveal an actual mismatch.
