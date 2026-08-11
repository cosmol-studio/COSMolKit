# RDKit MMFF Explicit-Hydrogen Typing Probe

## Scope

- RDKit source commit: `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- RDKit runtime oracle: `2026.03.1` from the project `.venv`.
- Upstream function: `MMFFMolProperties::setMMFFHydrogenType` in
  `Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:2083-2373`.
- Rust function: `set_mmff_hydrogen_type` in
  `crates/cosmolkit-core/src/chemistry/forcefield/mmff/mol_properties.rs`.
- First active-profile failure: explicit-hydrogen row 1 of `smiles_small`.

## Probe Method

- Read row 1 from the pinned `forcefield_params.jsonl` golden.
- Ran pinned RDKit on `Chem.AddHs(Chem.MolFromSmiles("C=C"))` and recorded
  atom indices, neighbor indices/elements, MMFF atom types, and
  `MMFFHasAllMoleculeParams`.
- Added a temporary environment-gated probe at the corresponding Rust
  hydrogen neighbor/type boundary. It recorded molecule size, hydrogen and
  neighbor indices, neighbor atomic number, and the already-assigned heavy
  atom type. The probe was removed after collecting the result.
- Traced the result through the pinned C++ source and the Rust source block.

## First Boundary: Row 1

Input SMILES:

```text
C=C
```

| Field | RDKit | COSMolKit |
|---|---:|---:|
| AddHs atom count | `6` | `6` |
| Heavy-atom types | `[2, 2]` | `[2, 2]` |
| First hydrogen index | `2` | `2` |
| First hydrogen neighbor | atom `0`, carbon, type `2` | atom `0`, carbon, type `2` |
| Explicit-H atom types | `[5, 5, 5, 5]` | no result |
| Parameter coverage | `true` | `false` |

The topology and heavy-atom typing agree before the first hydrogen is typed.
RDKit enters the carbon/silicon branch and assigns atom type 5 unconditionally:

```cpp
case 6:
case 14:
  atomType = 5;
  break;
```

The corresponding Rust lines remain marked `RDKit❌❌` and the current
catch-all branch returns `UnsupportedFeature` for carbon. The public
`mmff_has_all_molecule_params` wrapper intentionally converts that structured
unsupported result to `false`, which is the observed coverage mismatch.

## Remaining Source Surface

The same Rust function is only partially ported. A complete repair must port
the whole source function rather than only the row-1 carbon case:

- carbon and silicon neighbors -> type 5;
- nitrogen neighbor-type groups -> types 23, 36, 27, or 28;
- oxygen types 49, 51, and 70 -> types 50, 52, and 31;
- oxygen type 6 environment classification -> types 24, 29, 33, or 21;
- other oxygen types -> type 21;
- phosphorus and sulfur neighbors -> type 71;
- source-compatible zero/untyped validity behavior.

The oxygen type-6 branch depends on two nested neighbor traversals and bond
orders. It cannot be replaced by a neighbor-element shortcut without changing
RDKit behavior.

## Conclusion

The first explicit-H failure is an incomplete-port failure, not an AddHs
topology mismatch or an unsupported RDKit molecule. The next implementation
step should reproduce all branches of `setMMFFHydrogenType`, preserving the
existing two-phase heavy-then-hydrogen typing order and adding branch-level
tests plus full-corpus explicit-H coverage/type comparisons.
