use cosmolkit::{
    MmffConformerOptimizationParams, MmffEvaluationParams, MmffOptimizationParams,
    MmffOptimizeMoleculeConfsResult, MmffOptimizeMoleculeResult, Molecule, OperationError,
    UffConformerOptimizationParams, UffEvaluationParams, UffOptimizationParams,
};
use serde::Deserialize;
use std::fs::File;
use std::io::{BufRead, BufReader};
#[path = "mmff_test_support.rs"]
mod parity_data;
const ENERGY_TOLERANCE: f64 = 1.0e-6;
const GRADIENT_TOLERANCE: f64 = 1.0e-6;
const COORDINATE_TOLERANCE: f64 = 1.0e-6;
const FORCEFIELD_PARITY_NONBONDED_THRESH: f64 = 100.0;
const FORCEFIELD_OPT_MAX_ITERS: usize = 200;
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
struct ForcefieldInitialEnergyResult {
    ok: bool,
    needs_more: Option<i32>,
    energy: Option<f64>,
    gradient: Option<Vec<f64>>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct ForcefieldOptimizedResult {
    ok: bool,
    needs_more: Option<i32>,
    energy: Option<f64>,
    coords: Option<Vec<[f64; 3]>>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct ForcefieldMultiOptimizedResult {
    ok: bool,
    conformer_results: Option<Vec<ForcefieldOptimizedResult>>,
    initial_coords: Option<Vec<Vec<[f64; 3]>>>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct EmbeddedForcefieldRecord {
    ok: bool,
    cxsmiles: Option<String>,
    coords: Option<Vec<[f64; 3]>>,
    uff: ForcefieldInitialEnergyResult,
    mmff: ForcefieldInitialEnergyResult,
    uff_optimized: Option<ForcefieldOptimizedResult>,
    mmff_optimized: Option<ForcefieldOptimizedResult>,
    uff_multi_optimized: Option<ForcefieldMultiOptimizedResult>,
    mmff_multi_optimized: Option<ForcefieldMultiOptimizedResult>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct ForcefieldParamsRecord {
    smiles: String,
    rdkit_ok: bool,
    uff: ForcefieldResult,
    mmff: ForcefieldResult,
    uff_explicit_h: ForcefieldResult,
    mmff_explicit_h: ForcefieldResult,
    embedded: Option<EmbeddedForcefieldRecord>,
    error: Option<String>,
}

fn load_golden() -> Vec<ForcefieldParamsRecord> {
    let path = parity_data::golden_path("forcefield_params.jsonl");
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate RDKit goldens with `{}`: {err}",
            path.display(),
            parity_data::regenerate_command()
        )
    });
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(idx, line)| {
            let line = line.unwrap_or_else(|err| {
                panic!("failed to read {} line {}: {err}", path.display(), idx + 1)
            });
            serde_json::from_str(&line).unwrap_or_else(|err| {
                panic!("failed to parse {} line {}: {err}", path.display(), idx + 1)
            })
        })
        .collect()
}

#[test]
fn forcefield_params_golden_has_one_record_per_smiles_library_entry() {
    let expected = parity_data::count_smiles_rows();
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        if let Some(embedded) = &record.embedded {
            assert!(
                embedded.ok || embedded.error.is_some(),
                "row {} ({}) has failed embedded forcefield golden without error detail",
                row_idx + 1,
                record.smiles
            );
        }
    }
    assert_eq!(
        records.len(),
        expected,
        "forcefield params golden row count must match the active parity corpus"
    );
}

#[test]
fn uff_has_all_molecule_params_matches_rdkit_golden() {
    let records = load_golden();
    let row_filter = std::env::var("COSMOLKIT_FORCEFIELD_ROW_FILTER")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());

    for (row_idx, record) in records.iter().enumerate() {
        if let Some(filter) = row_filter {
            if row_idx + 1 != filter {
                continue;
            }
        }

        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}) is rdkit not ok but has no error",
                row_idx + 1,
                record.smiles
            );
            continue;
        }

        let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
            panic!(
                "cosmolkit failed to parse row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });

        assert!(
            record.uff.ok,
            "row {} ({}) has RDKit UFF error: {:?}",
            row_idx + 1,
            record.smiles,
            record.uff.error
        );
        let expected_uff = record.uff.has_all.unwrap_or_else(|| {
            panic!(
                "row {} ({}) has RDKit UFF ok result without has_all",
                row_idx + 1,
                record.smiles
            )
        });
        let actual_uff = uff_has_all_molecule_params(&mol).unwrap_or_else(|err| {
            panic!(
                "COSMolKit UFF parameter coverage errored at row {} ({}), expected RDKit comparison: {err}",
                row_idx + 1,
                record.smiles
            )
        });
        assert_eq!(
            actual_uff,
            expected_uff,
            "UFF parameter coverage mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
    }
}

#[test]
fn forcefield_initial_energy_matches_rdkit_golden_for_all_embedded_rows() {
    let records = load_golden();
    let row_filter = std::env::var("COSMOLKIT_FORCEFIELD_ROW_FILTER")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());
    let mut uff_comparisons = 0;
    let mut mmff_comparisons = 0;

    for (row_idx, record) in records.iter().enumerate() {
        if let Some(filter) = row_filter
            && row_idx + 1 != filter
        {
            continue;
        }
        let Some(embedded) = &record.embedded else {
            continue;
        };
        if !embedded.ok {
            assert!(
                embedded.error.is_some(),
                "row {} ({}) has failed RDKit embedding without error detail",
                row_idx + 1,
                record.smiles
            );
            continue;
        }
        let cxsmiles = embedded.cxsmiles.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) has embedded forcefield golden without CXSMILES",
                row_idx + 1,
                record.smiles
            )
        });
        let expected_coords = embedded.coords.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) has embedded forcefield golden without coordinates",
                row_idx + 1,
                record.smiles
            )
        });
        let mol = Molecule::from_smiles(cxsmiles).unwrap_or_else(|err| {
            panic!(
                "COSMolKit failed to parse embedded CXSMILES at row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });
        assert_eq!(
            mol.conformers_3d().len(),
            1,
            "row {} ({}) must parse exactly one 3D conformer from RDKit CXSMILES",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            mol.conformers_3d()[0].coordinates(),
            expected_coords.as_slice(),
            "row {} ({}) initial coordinates must be the RDKit golden coordinates",
            row_idx + 1,
            record.smiles
        );

        if embedded.uff.ok {
            uff_comparisons += 1;
            assert_initial_energy_matches(
                row_idx + 1,
                &record.smiles,
                "UFF",
                embedded.uff.needs_more,
                embedded.uff.energy,
                embedded.uff.error.as_deref(),
                uff_optimize_molecule_confs(&mol, 1, 0, FORCEFIELD_PARITY_NONBONDED_THRESH, true)
                    .map(|result| result.conformer_results),
            );
        } else {
            assert!(embedded.uff.error.is_some());
        }
        if embedded.mmff.ok {
            mmff_comparisons += 1;
            assert_initial_energy_matches(
                row_idx + 1,
                &record.smiles,
                "MMFF",
                embedded.mmff.needs_more,
                embedded.mmff.energy,
                embedded.mmff.error.as_deref(),
                mmff_optimize_molecule_confs(
                    &mol,
                    1,
                    0,
                    "MMFF94",
                    FORCEFIELD_PARITY_NONBONDED_THRESH,
                    true,
                )
                .map(|result| result.conformer_results),
            );
        } else {
            assert!(embedded.mmff.error.is_some());
        }
    }

    assert!(uff_comparisons > 0, "golden must contain UFF energy rows");
    assert!(mmff_comparisons > 0, "golden must contain MMFF energy rows");
}

#[test]
fn forcefield_initial_gradient_matches_rdkit_golden_for_all_embedded_rows() {
    let records = load_golden();
    let mut uff_comparisons = 0;
    let mut mmff_comparisons = 0;

    for (row_idx, record) in records.iter().enumerate() {
        let Some(embedded) = &record.embedded else {
            continue;
        };
        if !embedded.ok {
            assert!(embedded.error.is_some());
            continue;
        }
        let cxsmiles = embedded.cxsmiles.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) has embedded forcefield golden without CXSMILES",
                row_idx + 1,
                record.smiles
            )
        });
        let mol = Molecule::from_smiles(cxsmiles).unwrap_or_else(|err| {
            panic!(
                "COSMolKit failed to parse embedded CXSMILES at row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });

        if embedded.uff.ok {
            uff_comparisons += 1;
            if let Some(expected_gradient) = embedded.uff.gradient.as_deref() {
                assert_initial_gradient_matches(
                    row_idx + 1,
                    &record.smiles,
                    "UFF",
                    expected_gradient,
                    embedded.uff.error.as_deref(),
                    uff_initial_gradient_for_parity(
                        &mol,
                        FORCEFIELD_PARITY_NONBONDED_THRESH,
                        -1,
                        true,
                    ),
                );
            } else {
                assert_unavailable_forcefield_result(
                    row_idx + 1,
                    &record.smiles,
                    "UFF",
                    embedded.uff.needs_more,
                    embedded.uff.energy,
                    embedded.uff.error.as_deref(),
                    mol.conformers_3d()[0].coordinates(),
                    uff_optimize_molecule_confs(
                        &mol,
                        1,
                        0,
                        FORCEFIELD_PARITY_NONBONDED_THRESH,
                        true,
                    )
                    .map(|result| {
                        (
                            result.conformer_results[0].needs_more,
                            result.conformer_results[0].energy,
                            result.molecule.conformers_3d()[0].coordinates().to_vec(),
                        )
                    }),
                );
            }
        } else {
            assert!(embedded.uff.error.is_some());
        }
        if embedded.mmff.ok {
            mmff_comparisons += 1;
            if let Some(expected_gradient) = embedded.mmff.gradient.as_deref() {
                assert_initial_gradient_matches(
                    row_idx + 1,
                    &record.smiles,
                    "MMFF",
                    expected_gradient,
                    embedded.mmff.error.as_deref(),
                    mmff_initial_gradient_for_parity(
                        &mol,
                        "MMFF94",
                        FORCEFIELD_PARITY_NONBONDED_THRESH,
                        -1,
                        true,
                    ),
                );
            } else {
                assert_unavailable_forcefield_result(
                    row_idx + 1,
                    &record.smiles,
                    "MMFF",
                    embedded.mmff.needs_more,
                    embedded.mmff.energy,
                    embedded.mmff.error.as_deref(),
                    mol.conformers_3d()[0].coordinates(),
                    mmff_optimize_molecule_confs(
                        &mol,
                        1,
                        0,
                        "MMFF94",
                        FORCEFIELD_PARITY_NONBONDED_THRESH,
                        true,
                    )
                    .map(|result| {
                        (
                            result.conformer_results[0].needs_more,
                            result.conformer_results[0].energy,
                            result.molecule.conformers_3d()[0].coordinates().to_vec(),
                        )
                    }),
                );
            }
        } else {
            assert!(embedded.mmff.error.is_some());
        }
    }

    assert!(uff_comparisons > 0, "golden must contain UFF gradient rows");
    assert!(
        mmff_comparisons > 0,
        "golden must contain MMFF gradient rows"
    );
}

#[test]
fn uff_single_conformer_final_coordinates_match_rdkit_golden_for_all_embedded_rows() {
    let records = load_golden();
    let eligible = records.iter().enumerate().filter(|(_, record)| {
        record
            .embedded
            .as_ref()
            .and_then(|embedded| embedded.uff_optimized.as_ref())
            .is_some_and(|uff_optimized| uff_optimized.ok)
    });
    let mut comparisons = 0;
    for (row_idx, record) in eligible {
        comparisons += 1;
        let embedded = record
            .embedded
            .as_ref()
            .expect("selected row must have embedded forcefield golden");
        let cxsmiles = embedded.cxsmiles.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) has embedded forcefield golden without CXSMILES",
                row_idx + 1,
                record.smiles
            )
        });
        let expected = embedded.uff_optimized.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) has no UFF optimized golden",
                row_idx + 1,
                record.smiles
            )
        });
        assert!(
            expected.error.is_none(),
            "row {} ({}) has RDKit UFF optimized-coordinate error: {:?}",
            row_idx + 1,
            record.smiles,
            expected.error
        );

        let mol = Molecule::from_smiles(cxsmiles).unwrap_or_else(|err| {
            panic!(
                "COSMolKit failed to parse embedded CXSMILES at row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });
        let actual = uff_optimize_molecule(
            &mol,
            FORCEFIELD_OPT_MAX_ITERS,
            FORCEFIELD_PARITY_NONBONDED_THRESH,
            -1,
            true,
        )
        .unwrap_or_else(|err| {
            panic!(
                "COSMolKit UFF optimized-coordinate parity errored at row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });

        let expected_needs_more = expected.needs_more.unwrap_or_else(|| {
            panic!(
                "row {} ({}) RDKit UFF optimized result has no needs_more",
                row_idx + 1,
                record.smiles
            )
        });
        assert_eq!(
            actual.needs_more,
            expected_needs_more,
            "row {} ({}) UFF optimized result-code mismatch",
            row_idx + 1,
            record.smiles
        );

        let expected_energy = expected.energy.unwrap_or_else(|| {
            panic!(
                "row {} ({}) RDKit UFF optimized result has no energy",
                row_idx + 1,
                record.smiles
            )
        });
        assert!(
            (actual.energy - expected_energy).abs() <= ENERGY_TOLERANCE,
            "row {} ({}) UFF optimized final-energy mismatch: actual={} expected={}",
            row_idx + 1,
            record.smiles,
            actual.energy,
            expected_energy
        );

        let actual_coords = actual.molecule.conformers_3d()[0].coordinates();
        if let Some(expected_coords) = expected.coords.as_ref() {
            assert_coordinate_matrix_close(
                row_idx + 1,
                &record.smiles,
                "UFF optimized final",
                actual_coords,
                expected_coords,
            );
        } else {
            assert_eq!(expected_needs_more, -1);
            assert_eq!(expected_energy, -1.0);
            assert_coordinate_matrix_close(
                row_idx + 1,
                &record.smiles,
                "UFF unavailable force field preserves initial coordinates",
                actual_coords,
                embedded
                    .coords
                    .as_deref()
                    .expect("embedded row has initial coordinates"),
            );
        }
    }
    assert!(comparisons > 0, "golden must contain UFF optimized rows");
}

#[test]
fn uff_multi_conformer_final_coordinates_match_rdkit_golden_for_all_embedded_rows() {
    let records = load_golden();
    let eligible = records.iter().enumerate().filter(|(_, record)| {
        record
            .embedded
            .as_ref()
            .and_then(|embedded| embedded.uff_multi_optimized.as_ref())
            .is_some_and(|uff_multi_optimized| uff_multi_optimized.ok)
    });
    let mut comparisons = 0;
    for (row_idx, record) in eligible {
        comparisons += 1;
        let embedded = record
            .embedded
            .as_ref()
            .expect("selected row must have embedded forcefield golden");
        let cxsmiles = embedded.cxsmiles.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) has embedded forcefield golden without CXSMILES",
                row_idx + 1,
                record.smiles
            )
        });
        let expected = embedded.uff_multi_optimized.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) has no UFF multi-conformer optimized golden",
                row_idx + 1,
                record.smiles
            )
        });
        assert!(
            expected.error.is_none(),
            "row {} ({}) has RDKit UFF multi-conformer optimized error: {:?}",
            row_idx + 1,
            record.smiles,
            expected.error
        );
        let expected_initial_coords = expected.initial_coords.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) RDKit UFF multi-conformer result has no initial coordinates",
                row_idx + 1,
                record.smiles
            )
        });
        let expected_results = expected.conformer_results.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) RDKit UFF multi-conformer result has no conformer results",
                row_idx + 1,
                record.smiles
            )
        });

        let topology_mol = Molecule::from_smiles(cxsmiles).unwrap_or_else(|err| {
            panic!(
                "COSMolKit failed to parse embedded CXSMILES at row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });
        let mol = molecule_with_3d_conformers(&topology_mol, expected_initial_coords)
        .unwrap_or_else(|err| {
            panic!(
                "COSMolKit failed to build UFF multi-conformer parity molecule at row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });
        assert_eq!(
            mol.conformers_3d().len(),
            expected_initial_coords.len(),
            "row {} ({}) UFF multi-conformer initial coordinate count mismatch",
            row_idx + 1,
            record.smiles
        );
        for (conf_idx, expected_coords) in expected_initial_coords.iter().enumerate() {
            assert_coordinate_matrix_close(
                row_idx + 1,
                &record.smiles,
                &format!("UFF multi-conformer initial conformer {conf_idx}"),
                mol.conformers_3d()[conf_idx].coordinates(),
                expected_coords,
            );
        }

        let actual = uff_optimize_molecule_confs(
            &mol,
            1,
            FORCEFIELD_OPT_MAX_ITERS,
            FORCEFIELD_PARITY_NONBONDED_THRESH,
            true,
        )
        .unwrap_or_else(|err| {
            panic!(
                "COSMolKit UFF multi-conformer optimized parity errored at row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });
        assert_eq!(
            actual.conformer_results.len(),
            expected_results.len(),
            "row {} ({}) UFF multi-conformer optimized result count mismatch",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            actual.molecule.conformers_3d().len(),
            expected_results.len(),
            "row {} ({}) UFF multi-conformer optimized coordinate count mismatch",
            row_idx + 1,
            record.smiles
        );

        for (conf_idx, (actual_result, expected_result)) in actual
            .conformer_results
            .iter()
            .zip(expected_results)
            .enumerate()
        {
            assert_optimized_result_matches(
                row_idx + 1,
                &record.smiles,
                &format!("UFF multi-conformer final conformer {conf_idx}"),
                actual_result.needs_more,
                actual_result.energy,
                actual.molecule.conformers_3d()[conf_idx].coordinates(),
                expected_result,
            );
        }
    }
    assert!(
        comparisons > 0,
        "golden must contain UFF multi-conformer optimized rows"
    );
}

#[test]
fn uff_torsion_match_order_optimizer_regressions_rows_34_and_47_match_rdkit() {
    for row in [34, 47] {
        assert_uff_optimizer_regression(row);
    }
}

#[test]
fn uff_optimizer_regression_row_81_matches_rdkit() {
    assert_uff_optimizer_regression(81);
}

#[test]
fn uff_optimizer_regression_row_113_matches_rdkit() {
    assert_uff_optimizer_regression(113);
}

fn assert_uff_optimizer_regression(row: usize) {
    let records = load_golden();
    let record = records
        .get(row - 1)
        .unwrap_or_else(|| panic!("force-field golden has no row {row}"));
    let embedded = record
        .embedded
        .as_ref()
        .unwrap_or_else(|| panic!("row {row} ({}) has no embedded golden", record.smiles));
    assert!(
        embedded.ok,
        "row {row} ({}) failed RDKit embedding",
        record.smiles
    );
    assert!(
        embedded.uff.ok,
        "row {row} ({}) failed RDKit UFF",
        record.smiles
    );

    let cxsmiles = embedded
        .cxsmiles
        .as_ref()
        .unwrap_or_else(|| panic!("row {row} ({}) has no embedded CXSMILES", record.smiles));
    let mol = Molecule::from_smiles(cxsmiles).unwrap_or_else(|err| {
        panic!(
            "COSMolKit failed to parse embedded CXSMILES at row {row} ({}): {err}",
            record.smiles
        )
    });

    let actual_has_all = uff_has_all_molecule_params(&mol).unwrap_or_else(|err| {
        panic!(
            "COSMolKit UFF parameter coverage errored at row {row} ({}): {err}",
            record.smiles
        )
    });
    assert_eq!(
        Some(actual_has_all),
        record.uff.has_all,
        "UFF parameter coverage mismatch at row {row} ({})",
        record.smiles
    );

    assert_initial_energy_matches(
        row,
        &record.smiles,
        "UFF torsion-match-order regression",
        embedded.uff.needs_more,
        embedded.uff.energy,
        embedded.uff.error.as_deref(),
        uff_optimize_molecule_confs(&mol, 1, 0, FORCEFIELD_PARITY_NONBONDED_THRESH, true)
            .map(|result| result.conformer_results),
    );

    let expected_gradient = embedded.uff.gradient.as_deref().unwrap_or_else(|| {
        panic!(
            "row {row} ({}) RDKit UFF result has no initial gradient",
            record.smiles
        )
    });
    assert_initial_gradient_matches(
        row,
        &record.smiles,
        "UFF torsion-match-order regression",
        expected_gradient,
        embedded.uff.error.as_deref(),
        uff_initial_gradient_for_parity(&mol, FORCEFIELD_PARITY_NONBONDED_THRESH, -1, true),
    );

    let expected_single = embedded.uff_optimized.as_ref().unwrap_or_else(|| {
        panic!(
            "row {row} ({}) has no RDKit single-conformer UFF result",
            record.smiles
        )
    });
    assert!(expected_single.ok);
    let actual_single = uff_optimize_molecule(
        &mol,
        FORCEFIELD_OPT_MAX_ITERS,
        FORCEFIELD_PARITY_NONBONDED_THRESH,
        -1,
        true,
    )
    .unwrap_or_else(|err| {
        panic!(
            "COSMolKit UFF single-conformer regression errored at row {row} ({}): {err}",
            record.smiles
        )
    });
    assert_optimized_result_matches(
        row,
        &record.smiles,
        "UFF torsion-match-order single-conformer final",
        actual_single.needs_more,
        actual_single.energy,
        actual_single.molecule.conformers_3d()[0].coordinates(),
        expected_single,
    );

    let expected_multi = embedded.uff_multi_optimized.as_ref().unwrap_or_else(|| {
        panic!(
            "row {row} ({}) has no RDKit multi-conformer UFF result",
            record.smiles
        )
    });
    assert!(expected_multi.ok);
    let expected_initial_coords = expected_multi.initial_coords.as_ref().unwrap_or_else(|| {
        panic!(
            "row {row} ({}) RDKit multi-conformer result has no initial coordinates",
            record.smiles
        )
    });
    let expected_results = expected_multi
        .conformer_results
        .as_ref()
        .unwrap_or_else(|| {
            panic!(
                "row {row} ({}) RDKit multi-conformer result has no conformer results",
                record.smiles
            )
        });
    let multi_mol = molecule_with_3d_conformers(&mol, expected_initial_coords).unwrap_or_else(|err| {
        panic!(
            "COSMolKit failed to build UFF multi-conformer regression molecule at row {row} ({}): {err}",
            record.smiles
        )
    });
    let actual_multi = uff_optimize_molecule_confs(
        &multi_mol,
        1,
        FORCEFIELD_OPT_MAX_ITERS,
        FORCEFIELD_PARITY_NONBONDED_THRESH,
        true,
    )
    .unwrap_or_else(|err| {
        panic!(
            "COSMolKit UFF multi-conformer regression errored at row {row} ({}): {err}",
            record.smiles
        )
    });
    assert_eq!(actual_multi.conformer_results.len(), expected_results.len());
    for (conf_idx, (actual_result, expected_result)) in actual_multi
        .conformer_results
        .iter()
        .zip(expected_results)
        .enumerate()
    {
        assert_optimized_result_matches(
            row,
            &record.smiles,
            &format!("UFF torsion-match-order multi-conformer final {conf_idx}"),
            actual_result.needs_more,
            actual_result.energy,
            actual_multi.molecule.conformers_3d()[conf_idx].coordinates(),
            expected_result,
        );
    }
}

fn assert_initial_energy_matches<T>(
    row: usize,
    smiles: &str,
    forcefield: &str,
    expected_needs_more: Option<i32>,
    expected_energy: Option<f64>,
    expected_error: Option<&str>,
    actual: Result<Vec<T>, impl std::fmt::Display>,
) where
    T: InitialEnergyResult,
{
    assert!(
        expected_error.is_none(),
        "row {row} ({smiles}) has RDKit {forcefield} initial-energy error: {expected_error:?}"
    );
    let actual = actual.unwrap_or_else(|err| {
        panic!(
            "COSMolKit {forcefield} initial-energy parity errored at row {row} ({smiles}): {err}"
        )
    });
    assert_eq!(
        actual.len(),
        1,
        "row {row} ({smiles}) {forcefield} must return one conformer result"
    );
    let expected_needs_more = expected_needs_more.unwrap_or_else(|| {
        panic!("row {row} ({smiles}) RDKit {forcefield} result has no needs_more")
    });
    let expected_energy = expected_energy
        .unwrap_or_else(|| panic!("row {row} ({smiles}) RDKit {forcefield} result has no energy"));
    assert_eq!(
        actual[0].needs_more(),
        expected_needs_more,
        "row {row} ({smiles}) {forcefield} max_iters=0 result-code mismatch"
    );
    let actual_energy = actual[0].energy();
    assert!(
        (actual_energy - expected_energy).abs() <= ENERGY_TOLERANCE,
        "row {row} ({smiles}) {forcefield} initial-energy mismatch: actual={actual_energy} expected={expected_energy}"
    );
}

fn assert_initial_gradient_matches(
    row: usize,
    smiles: &str,
    forcefield: &str,
    expected_gradient: &[f64],
    expected_error: Option<&str>,
    actual: Result<Vec<f64>, impl std::fmt::Display>,
) {
    assert!(
        expected_error.is_none(),
        "row {row} ({smiles}) has RDKit {forcefield} gradient error: {expected_error:?}"
    );
    let actual = actual.unwrap_or_else(|err| {
        panic!("COSMolKit {forcefield} gradient parity errored at row {row} ({smiles}): {err}")
    });
    assert_eq!(
        actual.len(),
        expected_gradient.len(),
        "row {row} ({smiles}) {forcefield} gradient length mismatch"
    );
    for (axis_idx, (actual_value, expected_value)) in
        actual.iter().zip(expected_gradient.iter()).enumerate()
    {
        assert!(
            (actual_value - expected_value).abs() <= GRADIENT_TOLERANCE,
            "row {row} ({smiles}) {forcefield} gradient mismatch at flat index {axis_idx}: actual={actual_value} expected={expected_value}"
        );
    }
}

fn assert_unavailable_forcefield_result(
    row: usize,
    smiles: &str,
    forcefield: &str,
    expected_needs_more: Option<i32>,
    expected_energy: Option<f64>,
    expected_error: Option<&str>,
    initial_coords: &[[f64; 3]],
    actual: Result<(i32, f64, Vec<[f64; 3]>), impl std::fmt::Display>,
) {
    assert!(
        expected_error.is_none(),
        "row {row} ({smiles}) has RDKit {forcefield} unavailable-force-field error: {expected_error:?}"
    );
    assert_eq!(
        expected_needs_more,
        Some(-1),
        "row {row} ({smiles}) RDKit {forcefield} missing gradient must be explained by the no-force-field sentinel"
    );
    assert_eq!(
        expected_energy,
        Some(-1.0),
        "row {row} ({smiles}) RDKit {forcefield} no-force-field energy sentinel mismatch"
    );
    let (actual_needs_more, actual_energy, actual_coords) = actual.unwrap_or_else(|err| {
        panic!(
            "COSMolKit {forcefield} unavailable-force-field parity errored at row {row} ({smiles}): {err}"
        )
    });
    assert_eq!(actual_needs_more, -1);
    assert_eq!(actual_energy, -1.0);
    assert_coordinate_matrix_close(
        row,
        smiles,
        &format!("{forcefield} unavailable force field preserves initial coordinates"),
        &actual_coords,
        initial_coords,
    );
}

fn assert_optimized_result_matches(
    row: usize,
    smiles: &str,
    label: &str,
    actual_needs_more: i32,
    actual_energy: f64,
    actual_coords: &[[f64; 3]],
    expected: &ForcefieldOptimizedResult,
) {
    assert!(
        expected.error.is_none(),
        "row {row} ({smiles}) RDKit {label} error: {:?}",
        expected.error
    );
    let expected_needs_more = expected
        .needs_more
        .unwrap_or_else(|| panic!("row {row} ({smiles}) RDKit {label} result has no needs_more"));
    assert_eq!(
        actual_needs_more, expected_needs_more,
        "row {row} ({smiles}) {label} result-code mismatch"
    );
    let expected_energy = expected
        .energy
        .unwrap_or_else(|| panic!("row {row} ({smiles}) RDKit {label} result has no energy"));
    assert!(
        (actual_energy - expected_energy).abs() <= ENERGY_TOLERANCE,
        "row {row} ({smiles}) {label} energy mismatch: actual={actual_energy} expected={expected_energy}"
    );
    let expected_coords = expected
        .coords
        .as_ref()
        .unwrap_or_else(|| panic!("row {row} ({smiles}) RDKit {label} result has no coordinates"));
    assert_coordinate_matrix_close(row, smiles, label, actual_coords, expected_coords);
}

fn assert_coordinate_matrix_close(
    row: usize,
    smiles: &str,
    label: &str,
    actual: &[[f64; 3]],
    expected: &[[f64; 3]],
) {
    assert_eq!(
        actual.len(),
        expected.len(),
        "row {row} ({smiles}) {label} coordinate row count mismatch"
    );
    for (atom_idx, (actual_coord, expected_coord)) in actual.iter().zip(expected).enumerate() {
        for axis in 0..3 {
            assert!(
                (actual_coord[axis] - expected_coord[axis]).abs() <= COORDINATE_TOLERANCE,
                "row {row} ({smiles}) {label} coordinate mismatch at atom {atom_idx} axis {axis}: actual={} expected={}",
                actual_coord[axis],
                expected_coord[axis]
            );
        }
    }
}

fn molecule_with_3d_conformers(
    template: &Molecule,
    conformer_coords: &[Vec<[f64; 3]>],
) -> Result<Molecule, OperationError> {
    Molecule::from_parts(
        template.topology().clone(),
        cosmolkit::CoordinateBlock {
            conformers_3d: conformer_coords
                .iter()
                .enumerate()
                .map(|(id, coords)| cosmolkit::Conformer3D::new(id, coords.clone(), true))
                .collect(),
            ..Default::default()
        },
        template.properties().clone(),
    )
}

fn mmff_has_all_molecule_params(mol: &Molecule) -> Result<bool, cosmolkit::MmffMolPropertiesError> {
    mol.mmff_has_all_molecule_params()
}
fn mmff_optimize_molecule(
    mol: &Molecule,
    variant: &str,
    max: usize,
    threshold: f64,
    id: i32,
    ignore: bool,
) -> Result<MmffOptimizeMoleculeResult, OperationError> {
    assert!(
        id >= -1,
        "this original Rust oracle selects -1 or nonnegative conformer IDs"
    );
    mol.with_mmff_optimized_with_params(&MmffOptimizationParams {
        mmff_variant: variant.into(),
        max_iterations: i32::try_from(max).unwrap(),
        non_bonded_threshold: threshold,
        conformer_id: (id >= 0).then_some(id as usize),
        ignore_interfragment_interactions: ignore,
    })
}
fn mmff_optimize_molecule_confs(
    mol: &Molecule,
    threads: i32,
    max: usize,
    variant: &str,
    threshold: f64,
    ignore: bool,
) -> Result<MmffOptimizeMoleculeConfsResult, OperationError> {
    mol.with_mmff_optimized_conformers_with_params(&MmffConformerOptimizationParams {
        num_threads: threads,
        max_iterations: i32::try_from(max).unwrap(),
        mmff_variant: variant.into(),
        non_bonded_threshold: threshold,
        ignore_interfragment_interactions: ignore,
    })
}
fn mmff_initial_gradient_for_parity(
    mol: &Molecule,
    variant: &str,
    threshold: f64,
    id: i32,
    ignore: bool,
) -> Result<Vec<f64>, OperationError> {
    assert!(id >= -1);
    let result = mol.mmff_energy_gradient_with_params(&MmffEvaluationParams {
        mmff_variant: variant.into(),
        non_bonded_threshold: threshold,
        conformer_id: (id >= 0).then_some(id as usize),
        ignore_interfragment_interactions: ignore,
    })?;
    Ok(result
        .expect("original gradient-positive oracle row must have an available force field")
        .gradient()
        .to_vec())
}
trait InitialEnergyResult {
    fn needs_more(&self) -> i32;
    fn energy(&self) -> f64;
}
impl InitialEnergyResult for cosmolkit::MmffOptimizeMoleculeConfResult {
    fn needs_more(&self) -> i32 {
        self.needs_more
    }
    fn energy(&self) -> f64 {
        self.energy
    }
}

// Test-only transport of original raw result records from the canonical APIs;
// every expected value/assertion above remains the original source condition.
struct OriginalUffSingleRecord {
    molecule: Molecule,
    needs_more: i32,
    energy: f64,
}
struct OriginalUffConformerRecord {
    needs_more: i32,
    energy: f64,
}
struct OriginalUffMultiRecord {
    molecule: Molecule,
    conformer_results: Vec<OriginalUffConformerRecord>,
}
impl InitialEnergyResult for OriginalUffConformerRecord {
    fn needs_more(&self) -> i32 {
        self.needs_more
    }
    fn energy(&self) -> f64 {
        self.energy
    }
}
fn uff_has_all_molecule_params(mol: &Molecule) -> Result<bool, cosmolkit::UffParameterQueryError> {
    mol.with_assigned_valence()
        .map_err(cosmolkit::UffParameterQueryError::Cache)?
        .uff_has_all_molecule_params()
}
fn uff_optimize_molecule(
    mol: &Molecule,
    max: usize,
    threshold: f64,
    id: i32,
    ignore: bool,
) -> Result<OriginalUffSingleRecord, OperationError> {
    assert!(id >= -1);
    let value = mol
        .with_assigned_valence()?
        .with_uff_optimized_with_params(&UffOptimizationParams {
            max_iterations: i32::try_from(max).unwrap(),
            vdw_threshold: threshold,
            conformer_id: (id >= 0).then_some(id as usize),
            ignore_interfragment_interactions: ignore,
        })?;
    Ok(OriginalUffSingleRecord {
        molecule: value.molecule,
        needs_more: value.status,
        energy: value.energy,
    })
}
fn uff_optimize_molecule_confs(
    mol: &Molecule,
    threads: i32,
    max: usize,
    threshold: f64,
    ignore: bool,
) -> Result<OriginalUffMultiRecord, OperationError> {
    let value = mol
        .with_assigned_valence()?
        .with_uff_optimized_conformers_with_params(&UffConformerOptimizationParams {
            num_threads: threads,
            max_iterations: i32::try_from(max).unwrap(),
            vdw_threshold: threshold,
            ignore_interfragment_interactions: ignore,
        })?;
    Ok(OriginalUffMultiRecord {
        molecule: value.molecule,
        conformer_results: value
            .conformers
            .into_iter()
            .map(|row| OriginalUffConformerRecord {
                needs_more: row.status,
                energy: row.energy,
            })
            .collect(),
    })
}
fn uff_initial_gradient_for_parity(
    mol: &Molecule,
    threshold: f64,
    id: i32,
    ignore: bool,
) -> Result<Vec<f64>, OperationError> {
    assert!(id >= -1);
    Ok(mol
        .with_assigned_valence()?
        .uff_energy_gradient_with_params(&UffEvaluationParams {
            vdw_threshold: threshold,
            conformer_id: (id >= 0).then_some(id as usize),
            ignore_interfragment_interactions: ignore,
        })?
        .gradient()
        .to_vec())
}
