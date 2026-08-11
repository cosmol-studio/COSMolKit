# RDKit MMFF Bond Empirical-Rule Probe

## Scope

- RDKit source commit: `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- RDKit runtime oracle: `2026.03.1` from the project `.venv`.
- Upstream functions: `MMFFMolProperties::getMMFFBondStretchParams` and
  `MMFFMolProperties::getMMFFBondStretchEmpiricalRuleParams` in
  `Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp`.
- Rust function: `MmffMolProperties::get_mmff_bond_stretch_params`.
- Corpus failures: rows 103 and 123 of `smiles_small`.

## Probe 1

- Checked: parsed atom elements, MMFF atom types, MMFF bond type, public
  `GetMMFFBondStretchParams` result, and Rust force-field construction result.
- Last matching boundary: parsed topology, atom indices, element numbers,
  atom types, and bond type.
- First differing boundary: RDKit returns a bond parameter while Rust returns
  `UnsupportedFeature` during force-field construction.
- Reduced range: tabulated bond lookup through the empirical fallback in
  `getMMFFBondStretchParams`.
- Next probe: record whether the same bond key hits the tabulated MMFF bond
  collection on both sides.

## Probe 2

- Checked: every bond lookup preceding the failure, including atom indices,
  atomic numbers, atom types, bond type, and tabulated hit/miss status.
- Last matching boundary: the tabulated `MMFFBondCollection` lookup misses on
  the same triggering bond.
- First differing boundary: RDKit calls
  `getMMFFBondStretchEmpiricalRuleParams`; Rust returns the explicit unported
  error at that branch.
- Reduced range: `AtomTyper.cpp:2573-2729`, the complete bond-stretch
  empirical-rule helper.
- Next probe: compare its covalent-radius/electronegativity inputs, Bndk
  reference row, equation (18) `r0`, and equation (19) `kb`.

Triggering boundaries:

| Corpus row | Bond | Atomic numbers | Atom types | Bond type | Tabulated |
|---:|---|---|---|---:|---|
| 103 | 21-22 | 6-6 | 63-22 | 0 | miss |
| 123 | 0-1 | 6-8 | 60-6 | 0 | miss |

## Probe 3

- Checked: `MMFFCovRadPauEle`, `MMFFBndk`, equation (18), equation (19), and
  the final upstream public parameter tuple.
- Last matching boundary: the default parameter rows and all formula inputs.
- First differing boundary: Rust has no structures, default collections, or
  execution of the empirical helper after the tabulated miss.
- Reduced range: the missing Rust reproduction of the complete helper and its
  three supporting default parameter collections; no downstream builder or
  optimizer behavior is involved in this failure.
- Next probe: none; the first divergent source block is identified.

Exact source-backed calculations:

| Row | `r0_i` | `chi_i` | `c` | Bndk `(r0, kb)` | Empirical `r0` | Empirical `kb` |
|---:|---|---|---:|---|---:|---:|
| 103 | 0.77, 0.77 | 2.50, 2.50 | 0.085 | 1.512, 3.80 | 1.54 | 3.403846905179782 |
| 123 | 0.77, 0.72 | 2.50, 3.50 | 0.085 | 1.393, 5.40 | 1.405 | 5.129115902527102 |

The values follow the upstream source without special cases:

```text
r0 = r0_i + r0_j - c * abs(chi_i - chi_j)^1.4
kb = bndk.kb * (bndk.r0 / r0)^6
```

The Herschbach-Laurie branch remains part of the required complete port even
though these two corpus triggers use a Bndk reference row.

## Conclusion

The unsupported results are incomplete-port failures, not expected RDKit
behavior. The source-backed repair is to reproduce the full empirical helper,
including Bndk scaling and Herschbach-Laurie fallback, and to preserve the
existing tabulated-first dispatch. No molecule-specific or atom-type-specific
patch is justified.
