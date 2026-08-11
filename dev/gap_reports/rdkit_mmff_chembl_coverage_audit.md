# RDKit MMFF ChEMBL Coverage Audit

## Scope

- Reference: RDKit `2026.03.1`.
- Source corpus: ChEMBL 37 chemical representations.
- Deterministic corpus selection: the existing 1-of-16 sample in
  `tmp/parity-audit/chembl_sample_1of16.jsonl`.
- Source sample SHA-256:
  `fe2fac7fcc0bf4aa913bac0a0fb2da0076d7092c79d6529f4573d850e092d86a`.
- Resource boundary: RDKit-parseable molecules with at most 80 atoms.
- Scanned records: 181,306.
- Selected records: 178,531.
- RDKit parse failures: 0.
- Records skipped above the atom boundary: 2,775.

The audit uses the declared lightweight large-corpus force-field boundary. It
compares implicit- and explicit-hydrogen UFF/MMFF parameter availability, MMFF
atom types, and MMFF formal and partial charges. It does not perform conformer
embedding or optimization and does not claim large-corpus optimizer parity.

The selected SMILES file and RDKit oracle remain generated audit artifacts:

```text
tmp/parity-audit/chembl_mmff_coverage_le80.smi
SHA-256 7aecaa0da97456858b821fb3336dc710a75f5d0415697836cd9f0f11c33ebabf

tmp/parity-audit/chembl_mmff_coverage_le80.rdkit.jsonl
SHA-256 387c1535a4f037cdf5991bf6d2a08597156cfd1457c8daf164604ac525655acf
```

## Failure Classes

### Aromatic five-ring unlisted-element fallthrough

- Filtered audit row: 22,671.
- ChEMBL source row: 367,491.
- ChEMBL id: `CHEMBL362680`.
- Trigger: a phosphorus atom in a five-member MMFF-aromatic ring.
- RDKit result: complete MMFF parameters; the ring phosphorus is atom type 75.
- COSMolKit result before the fix: `UnsupportedFeature` from the aromatic
  five-ring switch.

RDKit `MMFFMolProperties::setMMFFHeavyAtomType()` has no `default` case in the
five-member aromatic-ring element switch. Elements other than C, N, O, and S
retain `atomType == 0` and continue to the element-specific aliphatic switch.
The phosphorus then reaches the source `getTotalDegree() == 2` branch and is
assigned type 75. COSMolKit incorrectly returned from the first switch.

The source-level fix restores the fallthrough for every unlisted element. The
full molecule is retained in `smiles_small.smi`, and the focused regression
checks the complete atom-type vector plus phosphorus partial charges.

### Degree-zero hydrogen component in AddHs

- Filtered audit row: 148,579.
- ChEMBL source row: 2,412,464.
- ChEMBL id: `CHEMBL3209673`.
- Trigger: a disconnected `[H+]` component.
- RDKit result: `AddHs` preserves the proton and adds hydrogens to eligible
  heavy atoms.
- COSMolKit result before the fix: the operation-level
  `hydrogen_ownership_represented` precondition rejected the molecule.

RDKit `MolOps::addHs()` has no ownership precondition. A degree-zero hydrogen
has zero explicit and implicit hydrogens to materialize, so it remains
unchanged while the normal heavy-atom loop continues. COSMolKit had both an
operation-registry precondition and an assignment-planner rejection that were
absent from the source.

Both rejection paths were removed. The smallest-boundary `C.[H+]` regression
checks that four carbon hydrogens are added while the charged degree-zero
hydrogen remains disconnected. The original ChEMBL row is also retained in
`smiles_small.smi`.

## Commands

Prepare the bounded audit corpus and oracle:

```bash
.venv/bin/python tools/testdata/rdkit/audit_forcefield_coverage.py \
  --input-jsonl tmp/parity-audit/chembl_sample_1of16.jsonl \
  --smiles-output tmp/parity-audit/chembl_mmff_coverage_le80.smi \
  --golden-output tmp/parity-audit/chembl_mmff_coverage_le80.rdkit.jsonl \
  --max-atoms 80
```

Run the Rust comparison:

```bash
COSMOLKIT_PARITY_SMILES="$PWD/tmp/parity-audit/chembl_mmff_coverage_le80.smi" \
COSMOLKIT_FORCEFIELD_COVERAGE_GOLDEN="$PWD/tmp/parity-audit/chembl_mmff_coverage_le80.rdkit.jsonl" \
cargo test -p cosmolkit-core --test rdkit_forcefield_coverage_parity
```

The post-fix continuation from row 148,579 through row 178,531 passed in
232.69 seconds. The final row-1 audit then passed all 178,531 selected records
in four deterministic shards with zero mismatches. The shard runtimes were
329.71, 331.31, 331.65, and 332.27 seconds.
