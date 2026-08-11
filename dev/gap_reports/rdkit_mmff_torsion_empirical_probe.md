# RDKit MMFF Torsion Empirical-Rule Probe

## Scope

- RDKit source commit: `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- RDKit runtime oracle: `2026.3.1` from the project `.venv`.
- Upstream functions: `MMFFMolProperties::getMMFFTorsionParams` and
  `MMFFMolProperties::getMMFFTorsionEmpiricalRuleParams` in
  `Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp`.
- Rust function: `MmffMolProperties::get_mmff_torsion_params`.
- First active-profile table-miss failure after the bond and angle ports: row
  103 of `smiles_small`.

## Probe Method

- Queried pinned RDKit's public torsion parameter boundary for the exact
  four-atom paths and both endpoint orders.
- Added a temporary environment-gated probe at the corresponding Rust torsion
  table-miss boundary. It recorded atom indices and types, atomic numbers,
  principal/secondary torsion types, central atom properties, and central-bond
  state. The probe was removed after collecting the results.
- Traced the observed results through the pinned source helper at
  `AtomTyper.cpp:2874-3067` and its dispatch at `AtomTyper.cpp:3629-3670`.

## First Boundary: Row 103

SMILES:

```text
O=C(C1=C(C2CC2)N(C3=C4C=CC=NC4=CC=C3)N=C1)NC(N)=N.[H]Cl.[H]O[H]
```

| Field | Value |
|---|---|
| Torsion indices | `10-21-22-23` |
| Atom types | `39/63/22/22` |
| Atomic numbers | `7/6/6/6` |
| Torsion type pair | `0/0` |
| Tabulated lookup | absent |
| J property | `crd=3`, `val=4`, `pilp=0`, `mltb=2`, `linh=0` |
| K property | `crd=4`, `val=4`, `pilp=0`, `mltb=0`, `linh=0` |
| Central bond | single, non-aromatic |
| RDKit forward/reverse result | `None` / `None` |
| Rust result before port | `UnsupportedFeature` |

The helper takes rule (f). J has `crd=3` and `val=4`, so RDKit sets all three
torsion coefficients to zero. `getMMFFTorsionParams` then returns false and the
builder excludes this torsion contribution. A table miss is therefore not an
operation-wide parameter-availability failure.

## Independent Boundary: Row 123

| Field | Value |
|---|---|
| Torsion indices | `3-5-6-7` |
| Atom types | `2/4/4/22` |
| Atomic numbers | `6/6/6/6` |
| Torsion type pair | `1/0` |
| Tabulated lookup | absent |
| J/K properties | both `crd=2`, `val=4`, `mltb=3`, `linh=1` |
| Central bond | single, non-aromatic |
| RDKit forward/reverse result | `None` / `None` |
| Rust result before port | `UnsupportedFeature` |

This case takes rule (a), which also sets all coefficients to zero because both
central atoms are linear. It independently locks the linear path rather than
only the row-103 coordination path.

## Conclusion

The failures are incomplete-port failures, not unsupported RDKit inputs. The
source-backed repair is the complete helper, including element constants and
rules (a) through (h), followed by the existing all-zero exclusion check. It
must not special-case the two observed zero-result molecules: other table
misses can produce nonzero `V1`, `V2`, or `V3` through the same helper.

The row-19 and row-21 numerical energy mismatches occur after successful
torsion lookup and are not attributed by this probe. They remain executable
parity failures for the later complete term/gradient/optimizer investigation.
