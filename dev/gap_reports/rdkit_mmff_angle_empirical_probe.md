# RDKit MMFF Angle Empirical-Rule Probe

## Scope

- RDKit source commit: `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- RDKit runtime oracle: `2026.03.1` from the project `.venv`.
- Upstream functions: `MMFFMolProperties::getMMFFAngleBendParams` and
  `getMMFFAngleBendEmpiricalRuleParams` in
  `Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp`.
- Rust function: `MmffMolProperties::get_mmff_angle_bend_params`.
- First active-profile failure: row 17 of `smiles_small`,
  `[13CH3:7][C@H](F)Cl`.

## Probe 1

- Enumerated every RDKit angle in the embedded row and recorded atom indices,
  MMFF atom types, both bond-stretch parameters, and the public angle result.
- Added a temporary environment-gated probe at the corresponding Rust angle
  table boundary and removed it after collecting the result.
- Last matching boundary: topology, angle indices, atom types, angle type,
  tabulated angle row, central-atom properties, bond parameters, and ring state.
- First differing boundary: RDKit invokes the angle empirical helper when the
  tabulated row is absent or has a zero force constant; Rust returned
  `UnsupportedFeature` at that branch.
- Reduced source range: `AtomTyper.cpp:2731-2872` plus the dispatch in
  `AtomTyper.cpp:3520-3562`.

Exact first boundary:

| Field | Value |
|---|---|
| Angle indices | `2-1-3` |
| Atom types | `11/1/12` |
| Atomic numbers | `9/6/17` |
| Angle type | `0` |
| Tabulated angle | `ka=0.0`, `theta0=108.9` |
| Central property | `crd=4`, `val=4`, `mltb=0`, `linh=0` |
| First bond | `kb=6.011`, `r0=1.36` |
| Second bond | `kb=2.974`, `r0=1.773` |
| Ring size 3 or 4 | none |
| RDKit final angle | `ka=1.2566039721725888`, `theta0=108.9` |
| Rust result before port | `UnsupportedFeature` |

The retained `theta0=108.9` proves this case takes the upstream
`oldMMFFAngleParams` branch. Equation (20) supplies only the force constant.
The complete port must also cover absent old rows, central coordination and
element branches, and three- and four-membered ring scaling.

## Conclusion

The angle failures are incomplete-port failures, not unsupported RDKit input.
The source-backed repair is the complete empirical helper and its existing
tabulated-or-zero dispatch. No molecule, row, atom-type, or element-specific
shortcut is justified.
