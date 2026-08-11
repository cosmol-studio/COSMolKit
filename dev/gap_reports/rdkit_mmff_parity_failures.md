# RDKit MMFF Parity Failure Inventory

## Scope

- Reference: RDKit `2026.03.1`.
- COSMolKit: current workspace source and Python extension.
- Corpus: all 150 rows in `testdata/smiles/corpus/smiles_small.smi`.
- Coordinates: the deterministic RDKit CXSMILES and multi-conformer coordinates in `forcefield_params.jsonl`.
- Thresholds: `1e-6` absolute energy and coordinate error, matching the integration-test contract.
- MMFF profile: `MMFF94`, 200 optimizer iterations, non-bonded threshold 100, and ignored inter-fragment interactions.

## Test-Coverage Finding

The former initial-energy, initial-gradient, single-conformer, and
multi-conformer tests selected only the first successful embedded row. That row
was `C=C`. The tests now execute every RDKit-successful row for the applicable
surface. The golden schema now also records explicit-hydrogen UFF/MMFF
parameter availability and MMFF atom types.

The first all-row test run passed only the three implicit-topology checks and
the golden row-count check. All optimizer/gradient surfaces and the new
explicit-hydrogen surface exposed failures.

## MMFF Inventory

Of 150 corpus rows, 138 parse successfully in RDKit and 135 also have a
successful deterministic embedded MMFF record.

| Surface | Match | Numeric/status mismatch | COSMolKit error | RDKit embedding failure |
|---|---:|---:|---:|---:|
| Explicit-H `has_all` | 24 | 114 | 0 | n/a |
| Initial energy | 49 | 67 | 19 | 3 |
| Single-conformer optimization | 47 | 69 | 19 | 3 |
| Two-conformer optimization | 44 | 72 | 19 | 3 |

The 19 COSMolKit errors are deterministic source-port gaps:

- 17 rows require the unported RDKit MMFF angle-bend empirical rule.
- 2 rows require the unported RDKit MMFF bond-stretch empirical rule.
- No row in this 150-row profile reaches the unported torsion empirical rule
  before another comparison outcome is recorded.

Representative failures:

| Row | SMILES | First observed failure |
|---:|---|---|
| 1 | `C=C` | RDKit explicit-H parameters true, COSMolKit false |
| 17 | `[13CH3:7][C@H](F)Cl` | angle-bend empirical rule unsupported |
| 19 | `C[C@H](O)[C@@H](C)O` | initial energy error `-5.373040665`; optimized coordinate error `0.003461206` |
| 20 | `N[C@@H](C)C(=O)O` | optimized energy error `-14.713178081`; coordinate error `0.139677563` |
| 103 | `O=C(C1=C(C2CC2)N(C3=C4C=CC=NC4=CC=C3)N=C1)NC(N)=N.[H]Cl.[H]O[H]` | bond-stretch empirical rule unsupported |

## Unsupported Rows

Angle-bend empirical-rule rows:

```text
17, 18, 21, 22, 48, 49, 50, 51, 72, 73, 76, 78, 79, 90, 99, 104, 129
```

Bond-stretch empirical-rule rows:

```text
103, 123
```

These are incomplete COSMolKit ports. RDKit succeeds on every one of these
rows, so `UnsupportedFeature` is not a parity result.

## Oracle Sentinel Case

`O=O` is not an optimizer-success row. RDKit reports no complete MMFF
parameters and returns its documented force-field sentinel:

```text
needs_more = -1
energy = -1.0
optimized coordinates = null for the single-force-field golden
```

The all-row test initially diagnosed the null coordinates as a malformed
golden. That is a test bug. The parity test must instead assert the `-1`
sentinel and unchanged coordinates. The initial-gradient test must likewise
distinguish an absent RDKit force field from a successful force field with a
gradient.

## UFF Findings Exposed By The Same Test Repair

The stronger shared force-field tests also exposed UFF mismatches that the
first-row selector hid:

- Single-conformer row 47, `c1cc[se]c1`, differs by approximately
  `0.023599 A` at the first reported coordinate.
- Multi-conformer row 29, benzene, differs by approximately
  `1.055e-6 A`, just outside the declared tolerance.

These remain separate from the MMFF source-port sequence but must not be
removed or tolerance-relaxed.

## Required Repair Order

1. Correct the RDKit no-force-field sentinel assertions without skipping the row.
2. Instrument and port the complete bond-stretch empirical rule.
3. Instrument and port the complete angle-bend empirical rule.
4. Instrument and port the complete torsion empirical rule even though this corpus does not reach it first.
5. Instrument and port complete explicit-hydrogen MMFF typing.
6. Bisect term construction, gradients, minimizer iterations, and termination for the remaining numeric mismatches.
7. Run the same complete schema on `smiles_5000` and the deterministic ChEMBL audit corpus.
