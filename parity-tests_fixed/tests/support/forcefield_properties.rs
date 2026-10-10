use cosmolkit::Molecule;
use serde::Deserialize;

const CHARGE_TOLERANCE: f64 = 1.0e-12;

#[derive(Debug, Deserialize)]
struct ForcefieldResult {
    ok: bool,
    has_all: Option<bool>,
    #[serde(default)]
    atom_types: Option<Vec<u8>>,
    #[serde(default)]
    formal_charges: Option<Vec<f64>>,
    #[serde(default)]
    partial_charges: Option<Vec<f64>>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct ForcefieldCoverageRecord {
    smiles: String,
    rdkit_ok: bool,
    uff: ForcefieldResult,
    mmff: ForcefieldResult,
    uff_explicit_h: ForcefieldResult,
    mmff_explicit_h: ForcefieldResult,
    error: Option<String>,
}

fn load_golden() -> Vec<ForcefieldCoverageRecord> {
    let snapshot = cosmolkit_parity_tests_fixed::special_regression::preflight(
        "forcefield_properties",
        &cosmolkit_parity_tests_fixed::expected(),
    )
    .unwrap();
    snapshot
        .rows
        .into_iter()
        .map(|row| serde_json::from_value(row).unwrap())
        .collect()
}

fn assert_uff_coverage(
    row: usize,
    smiles: &str,
    mol: &Molecule,
    expected: &ForcefieldResult,
    surface: &str,
) {
    assert!(
        expected.ok,
        "row {row} ({smiles}) has RDKit {surface} UFF error: {:?}",
        expected.error
    );
    let expected_has_all = expected.has_all.unwrap_or_else(|| {
        panic!("row {row} ({smiles}) has RDKit {surface} UFF result without has_all")
    });
    let actual_has_all = mol.uff_has_all_molecule_params().unwrap_or_else(|err| {
        panic!("row {row} ({smiles}) COSMolKit {surface} UFF coverage errored: {err}")
    });
    assert_eq!(
        actual_has_all, expected_has_all,
        "row {row} ({smiles}) {surface} UFF parameter coverage mismatch"
    );
}

fn assert_mmff_coverage(
    row: usize,
    smiles: &str,
    mol: &Molecule,
    expected: &ForcefieldResult,
    surface: &str,
) {
    assert!(
        expected.ok,
        "row {row} ({smiles}) has RDKit {surface} MMFF error: {:?}",
        expected.error
    );
    let expected_has_all = expected.has_all.unwrap_or_else(|| {
        panic!("row {row} ({smiles}) has RDKit {surface} MMFF result without has_all")
    });
    let actual_has_all = mol.mmff_has_all_molecule_params().unwrap_or_else(|err| {
        panic!("row {row} ({smiles}) COSMolKit {surface} MMFF coverage errored: {err}")
    });
    assert_eq!(
        actual_has_all, expected_has_all,
        "row {row} ({smiles}) {surface} MMFF parameter coverage mismatch"
    );

    let Some(expected_atom_types) = expected.atom_types.as_ref() else {
        assert!(
            !expected_has_all,
            "row {row} ({smiles}) has RDKit {surface} MMFF parameters without atom types"
        );
        assert!(expected.formal_charges.is_none());
        assert!(expected.partial_charges.is_none());
        return;
    };

    let props = mol.mmff_properties().unwrap_or_else(|err| {
        panic!("row {row} ({smiles}) COSMolKit {surface} MMFF properties errored: {err}")
    });
    let actual_atom_types = (0..mol.num_atoms())
        .map(|idx| {
            props.atom_type(idx).unwrap_or_else(|err| {
                panic!(
                    "row {row} ({smiles}) COSMolKit {surface} MMFF atom {idx} type errored: {err}"
                )
            })
        })
        .collect::<Vec<_>>();
    assert_eq!(
        actual_atom_types, *expected_atom_types,
        "row {row} ({smiles}) {surface} MMFF atom type mismatch"
    );

    let expected_formal = expected.formal_charges.as_ref().unwrap_or_else(|| {
        panic!("row {row} ({smiles}) RDKit {surface} MMFF result has no formal charges")
    });
    let expected_partial = expected.partial_charges.as_ref().unwrap_or_else(|| {
        panic!("row {row} ({smiles}) RDKit {surface} MMFF result has no partial charges")
    });
    assert_eq!(expected_formal.len(), mol.num_atoms());
    assert_eq!(expected_partial.len(), mol.num_atoms());
    for idx in 0..mol.num_atoms() {
        let actual_formal = props.formal_charge(idx).unwrap();
        let actual_partial = props.partial_charge(idx).unwrap();
        assert!(
            (actual_formal - expected_formal[idx]).abs() <= CHARGE_TOLERANCE,
            "row {row} ({smiles}) {surface} MMFF formal charge mismatch at atom {idx}: actual={actual_formal} expected={}",
            expected_formal[idx]
        );
        assert!(
            (actual_partial - expected_partial[idx]).abs() <= CHARGE_TOLERANCE,
            "row {row} ({smiles}) {surface} MMFF partial charge mismatch at atom {idx}: actual={actual_partial} expected={}",
            expected_partial[idx]
        );
    }
}

#[test]
fn forcefield_coverage_matches_rdkit_for_every_active_profile_row() {
    let records = load_golden();
    assert_eq!(
        records.len(),
        152,
        "force-field coverage golden row count must match the active corpus"
    );

    for (row_idx, record) in records.iter().enumerate() {
        let row = row_idx + 1;
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {row} ({}) is RDKit-not-ok without an error",
                record.smiles
            );
            continue;
        }

        let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
            panic!(
                "row {row} ({}) COSMolKit parse errored: {err}",
                record.smiles
            )
        });
        assert_uff_coverage(row, &record.smiles, &mol, &record.uff, "implicit-H");
        assert_mmff_coverage(row, &record.smiles, &mol, &record.mmff, "implicit-H");

        let explicit_h_mol = mol.with_hydrogens().unwrap_or_else(|err| {
            panic!(
                "row {row} ({}) COSMolKit explicit-H construction errored: {err}",
                record.smiles
            )
        });
        // The detached runtime deliberately invalidates valence after AddHs.
        // The source UFF typer reads a prepared property cache. Reproduce that
        // input boundary through the canonical registered valence operation;
        // no topology, expected vector, option, or tolerance is changed.
        let explicit_h_mol = explicit_h_mol
            .with_assigned_valence()
            .unwrap_or_else(|err| {
                panic!(
                    "row {row} ({}) explicit-H property cache preparation errored: {err}",
                    record.smiles
                )
            });
        assert_uff_coverage(
            row,
            &record.smiles,
            &explicit_h_mol,
            &record.uff_explicit_h,
            "explicit-H",
        );
        assert_mmff_coverage(
            row,
            &record.smiles,
            &explicit_h_mol,
            &record.mmff_explicit_h,
            "explicit-H",
        );
    }
}
