use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit_core::Molecule;
use cosmolkit_core::properties::descriptors::{
    CrippenDescriptorValues, DescriptorError, DescriptorResult, NumRotatableBondsOptions,
    calc_crippen_descriptors, calc_exact_mol_wt, calc_fraction_csp3, calc_mol_formula, calc_mol_wt,
    calc_num_aromatic_rings, calc_num_hba, calc_num_hbd, calc_num_rotatable_bonds, calc_qed,
    calc_tpsa,
};
use serde::Deserialize;

mod common;
use common::parity_data;

#[derive(Debug, Deserialize)]
struct DescriptorRecord {
    smiles: String,
    rdkit_ok: bool,
    descriptors: Option<DescriptorSet>,
    descriptor_bits: Option<DescriptorBits>,
    descriptor_options: Option<DescriptorOptionMatrix>,
    descriptor_option_bits: Option<DescriptorOptionBits>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct DescriptorOptionMatrix {
    mol_wt_only_heavy: f64,
    exact_mol_wt_only_heavy: f64,
    crippen: BTreeMap<String, CrippenOptionValues>,
    tpsa: BTreeMap<String, f64>,
}

#[derive(Debug, Deserialize)]
struct DescriptorOptionBits {
    mol_wt_only_heavy: String,
    exact_mol_wt_only_heavy: String,
    crippen: BTreeMap<String, CrippenOptionBits>,
    tpsa: BTreeMap<String, String>,
}

#[derive(Debug, Deserialize)]
struct CrippenOptionValues {
    logp: f64,
    molar_refractivity: f64,
}

#[derive(Debug, Deserialize)]
struct CrippenOptionBits {
    logp: String,
    molar_refractivity: String,
}

#[derive(Debug, Deserialize)]
struct DescriptorSet {
    mol_wt: f64,
    exact_mol_wt: f64,
    formula: String,
    formula_separate_isotopes: String,
    formula_separate_isotopes_no_h_abbrev: String,
    num_hbd: u32,
    num_hba: u32,
    fraction_csp3: f64,
    crippen_logp: f64,
    crippen_mr: f64,
    tpsa: f64,
    tpsa_include_sandp: f64,
    num_aromatic_rings: u32,
    num_rotatable_bonds_default: u32,
    num_rotatable_bonds_non_strict: u32,
    num_rotatable_bonds_strict: u32,
    num_rotatable_bonds_strict_linkages: u32,
    qed: f64,
}

#[derive(Debug, Deserialize)]
struct DescriptorBits {
    mol_wt: String,
    exact_mol_wt: String,
    fraction_csp3: String,
    crippen_logp: String,
    crippen_mr: String,
    tpsa: String,
    tpsa_include_sandp: String,
    qed: String,
}

fn load_golden() -> Vec<DescriptorRecord> {
    let path = parity_data::golden_path("molecular_descriptors.jsonl");
    assert!(
        path.exists(),
        "missing RDKit molecular descriptor golden: {}. Regenerate goldens with {}",
        path.display(),
        parity_data::regenerate_command()
    );
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!("failed to open {}: {err}", path.display());
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

fn record_context(record: &DescriptorRecord, row_idx: usize) -> String {
    format!("row {} ({})", row_idx + 1, record.smiles)
}

#[test]
fn molecular_descriptor_golden_has_one_record_per_smiles() {
    let records = load_golden();
    assert_eq!(
        records.len(),
        parity_data::count_smiles_rows(),
        "molecular descriptor golden row count must match the selected SMILES profile"
    );
}

#[test]
fn split_source_isotope_mass_row_matches_rdkit_descriptors() {
    let molecule = Molecule::from_smiles("[47Ca+2].[Cl-].[Cl-]").unwrap();

    assert_eq!(
        calc_mol_wt(&molecule, false).unwrap().to_bits(),
        0x405d_7713_2f87_ad08
    );
    assert_eq!(
        calc_mol_wt(&molecule, true).unwrap().to_bits(),
        0x405d_7713_2f87_ad08
    );
    assert_eq!(
        calc_exact_mol_wt(&molecule, false).unwrap().to_bits(),
        0x405d_391a_a572_c0bd
    );
    assert_eq!(
        calc_exact_mol_wt(&molecule, true).unwrap().to_bits(),
        0x405d_391a_a572_c0bd
    );
    assert_eq!(
        calc_qed(&molecule).unwrap().to_bits(),
        0x3fd1_bfb4_805a_d56c
    );
}

#[test]
fn crippen_force_and_include_hs_follow_rdkit_computed_property_cache() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let structural_snapshot = molecule.clone();

    let without_hs = calc_crippen_descriptors(&molecule, false, true).unwrap();
    assert_eq!(without_hs.logp.to_bits(), 0xbfd6_5119_ce07_5f70);
    assert_eq!(
        without_hs.molar_refractivity.to_bits(),
        0x4018_51b7_1758_e21a
    );

    let cached = calc_crippen_descriptors(&molecule, true, false).unwrap();
    assert_eq!(cached, without_hs);

    let forced_with_hs = calc_crippen_descriptors(&molecule, true, true).unwrap();
    assert_eq!(forced_with_hs.logp.to_bits(), 0xbf56_f006_8db8_bb00);
    assert_eq!(
        forced_with_hs.molar_refractivity.to_bits(),
        0x4029_8504_816f_006a
    );
    assert_eq!(molecule, structural_snapshot);
}

#[test]
fn crippen_computed_property_cache_is_copied_but_not_shared() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let without_hs = calc_crippen_descriptors(&molecule, false, true).unwrap();
    let clone = molecule.clone();

    let clone_with_hs = calc_crippen_descriptors(&clone, true, true).unwrap();
    assert_ne!(clone_with_hs, without_hs);
    assert_eq!(
        calc_crippen_descriptors(&molecule, true, false).unwrap(),
        without_hs
    );
    assert_eq!(
        calc_crippen_descriptors(&clone, false, false).unwrap(),
        clone_with_hs
    );
}

#[test]
fn rdkit_computed_property_clearing_operations_invalidate_the_crippen_cache() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let without_hs = calc_crippen_descriptors(&molecule, false, true).unwrap();

    let with_hydrogens = molecule.with_hydrogens().unwrap();
    let explicit_h_result = calc_crippen_descriptors(&with_hydrogens, false, false).unwrap();

    assert_ne!(explicit_h_result, without_hs);
    assert_eq!(explicit_h_result.logp.to_bits(), 0xbf56_f006_8db8_bb00);
    assert_eq!(
        explicit_h_result.molar_refractivity.to_bits(),
        0x4029_8504_816f_006a
    );

    let sanitized = molecule
        .sanitize_with_ops(cosmolkit_core::SanitizeOps::NONE)
        .unwrap();
    let sanitized_result = calc_crippen_descriptors(&sanitized, true, false).unwrap();
    assert_ne!(sanitized_result, without_hs);
    assert_eq!(sanitized_result.logp.to_bits(), 0xbf56_f006_8db8_bb00);

    let without_hydrogens = with_hydrogens
        .without_hydrogens_with_sanitize(false)
        .unwrap();
    let removed_h_result = calc_crippen_descriptors(&without_hydrogens, false, false).unwrap();
    assert_ne!(removed_h_result, explicit_h_result);
    assert_eq!(removed_h_result.logp.to_bits(), 0xbfd6_5119_ce07_5f70);
}

#[test]
fn rdkit_topology_state_operations_preserve_the_crippen_cache() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let without_hs = calc_crippen_descriptors(&molecule, false, true).unwrap();

    let kekulized = molecule.with_kekulized_bonds(false).unwrap();
    assert_eq!(
        calc_crippen_descriptors(&kekulized, true, false).unwrap(),
        without_hs
    );

    let with_radicals = molecule.with_assigned_radicals().unwrap();
    assert_eq!(
        calc_crippen_descriptors(&with_radicals, true, false).unwrap(),
        without_hs
    );

    let with_aromaticity = molecule.with_assigned_aromaticity().unwrap();
    assert_eq!(
        calc_crippen_descriptors(&with_aromaticity, true, false).unwrap(),
        without_hs
    );
}

#[test]
fn molecular_descriptors_match_rdkit_golden_for_supported_properties() {
    let records = load_golden();
    let mut failures = Vec::<String>::new();
    let mut by_field = BTreeMap::<&'static str, usize>::new();

    for (row_idx, record) in records.iter().enumerate() {
        let context = record_context(record, row_idx);
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "RDKit-not-ok molecular descriptor record has no error in {context}"
            );
            assert!(
                record.descriptors.is_none()
                    && record.descriptor_bits.is_none()
                    && record.descriptor_options.is_none()
                    && record.descriptor_option_bits.is_none(),
                "RDKit-not-ok row unexpectedly carries descriptor output in {context}"
            );
            assert!(
                Molecule::from_smiles(&record.smiles).is_err(),
                "RDKit rejected {context}, but COSMolKit accepted the same descriptor input"
            );
            continue;
        }
        let expected = record.descriptors.as_ref().unwrap_or_else(|| {
            panic!("RDKit-ok molecular descriptor record missing descriptor set in {context}")
        });
        let expected_bits = record.descriptor_bits.as_ref().unwrap_or_else(|| {
            panic!(
                "RDKit-ok molecular descriptor record missing descriptor_bits in {context}; regenerate goldens with {}",
                parity_data::regenerate_command()
            )
        });
        let expected_options = record.descriptor_options.as_ref().unwrap_or_else(|| {
            panic!("RDKit-ok molecular descriptor record missing option matrix in {context}")
        });
        let expected_option_bits = record.descriptor_option_bits.as_ref().unwrap_or_else(|| {
            panic!("RDKit-ok molecular descriptor record missing option bits in {context}")
        });
        let mol = Molecule::from_smiles(&record.smiles)
            .unwrap_or_else(|err| panic!("COSMolKit failed to parse {context}: {err}"));

        check_f64(
            "mol_wt_only_heavy",
            &context,
            calc_mol_wt(&mol, true),
            expected_options.mol_wt_only_heavy,
            expected_option_bits.mol_wt_only_heavy.as_str(),
            &mut failures,
            &mut by_field,
        );
        check_f64(
            "exact_mol_wt_only_heavy",
            &context,
            calc_exact_mol_wt(&mol, true),
            expected_options.exact_mol_wt_only_heavy,
            expected_option_bits.exact_mol_wt_only_heavy.as_str(),
            &mut failures,
            &mut by_field,
        );

        for include_hs in [false, true] {
            for force in [false, true] {
                let key = format!("include_hs_{include_hs}_force_{force}");
                let expected_branch = expected_options
                    .crippen
                    .get(&key)
                    .unwrap_or_else(|| panic!("missing Crippen branch {key} in {context}"));
                let expected_branch_bits = expected_option_bits
                    .crippen
                    .get(&key)
                    .unwrap_or_else(|| panic!("missing Crippen bit branch {key} in {context}"));
                // The golden generator starts each option branch from a fresh
                // RDKit molecule so one branch's computed-property cache cannot
                // mask the next branch's includeHs argument. Keep that branch
                // isolation here; cache-order parity has separate focused tests.
                let branch_mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
                    panic!("COSMolKit failed to parse Crippen branch {context}: {err}")
                });
                match calc_crippen_descriptors(&branch_mol, include_hs, force) {
                    Ok(CrippenDescriptorValues {
                        logp,
                        molar_refractivity,
                    }) => {
                        compare_f64(
                            "crippen_option_logp",
                            &format!("{context} {key}"),
                            logp,
                            expected_branch.logp,
                            expected_branch_bits.logp.as_str(),
                            &mut failures,
                            &mut by_field,
                        );
                        compare_f64(
                            "crippen_option_mr",
                            &format!("{context} {key}"),
                            molar_refractivity,
                            expected_branch.molar_refractivity,
                            expected_branch_bits.molar_refractivity.as_str(),
                            &mut failures,
                            &mut by_field,
                        );
                    }
                    Err(err) => {
                        record_error(
                            "crippen_option_logp",
                            &format!("{context} {key}"),
                            err.clone(),
                            &mut failures,
                            &mut by_field,
                        );
                        record_error(
                            "crippen_option_mr",
                            &format!("{context} {key}"),
                            err,
                            &mut failures,
                            &mut by_field,
                        );
                    }
                }
            }
        }

        for force in [false, true] {
            for include_sandp in [false, true] {
                let key = format!("force_{force}_include_sandp_{include_sandp}");
                let expected_value = *expected_options
                    .tpsa
                    .get(&key)
                    .unwrap_or_else(|| panic!("missing TPSA branch {key} in {context}"));
                let expected_bits = expected_option_bits
                    .tpsa
                    .get(&key)
                    .unwrap_or_else(|| panic!("missing TPSA bit branch {key} in {context}"));
                check_f64(
                    "tpsa_option",
                    &format!("{context} {key}"),
                    calc_tpsa(&mol, force, include_sandp),
                    expected_value,
                    expected_bits,
                    &mut failures,
                    &mut by_field,
                );
            }
        }

        check_f64(
            "mol_wt",
            &context,
            calc_mol_wt(&mol, false),
            expected.mol_wt,
            expected_bits.mol_wt.as_str(),
            &mut failures,
            &mut by_field,
        );
        check_f64(
            "exact_mol_wt",
            &context,
            calc_exact_mol_wt(&mol, false),
            expected.exact_mol_wt,
            expected_bits.exact_mol_wt.as_str(),
            &mut failures,
            &mut by_field,
        );
        check_string(
            "formula",
            &context,
            calc_mol_formula(&mol, false, true),
            &expected.formula,
            &mut failures,
            &mut by_field,
        );
        check_string(
            "formula_separate_isotopes",
            &context,
            calc_mol_formula(&mol, true, true),
            &expected.formula_separate_isotopes,
            &mut failures,
            &mut by_field,
        );
        check_string(
            "formula_separate_isotopes_no_h_abbrev",
            &context,
            calc_mol_formula(&mol, true, false),
            &expected.formula_separate_isotopes_no_h_abbrev,
            &mut failures,
            &mut by_field,
        );
        check_u32(
            "num_hbd",
            &context,
            calc_num_hbd(&mol),
            expected.num_hbd,
            &mut failures,
            &mut by_field,
        );
        check_u32(
            "num_hba",
            &context,
            calc_num_hba(&mol),
            expected.num_hba,
            &mut failures,
            &mut by_field,
        );
        check_f64(
            "fraction_csp3",
            &context,
            calc_fraction_csp3(&mol),
            expected.fraction_csp3,
            expected_bits.fraction_csp3.as_str(),
            &mut failures,
            &mut by_field,
        );
        match calc_crippen_descriptors(&mol, true, false) {
            Ok(CrippenDescriptorValues {
                logp,
                molar_refractivity,
            }) => {
                compare_f64(
                    "crippen_logp",
                    &context,
                    logp,
                    expected.crippen_logp,
                    expected_bits.crippen_logp.as_str(),
                    &mut failures,
                    &mut by_field,
                );
                compare_f64(
                    "crippen_mr",
                    &context,
                    molar_refractivity,
                    expected.crippen_mr,
                    expected_bits.crippen_mr.as_str(),
                    &mut failures,
                    &mut by_field,
                );
            }
            Err(err) => {
                record_error(
                    "crippen_logp",
                    &context,
                    err.clone(),
                    &mut failures,
                    &mut by_field,
                );
                record_error("crippen_mr", &context, err, &mut failures, &mut by_field);
            }
        }
        check_f64(
            "tpsa",
            &context,
            calc_tpsa(&mol, false, false),
            expected.tpsa,
            expected_bits.tpsa.as_str(),
            &mut failures,
            &mut by_field,
        );
        check_f64(
            "tpsa_include_sandp",
            &context,
            calc_tpsa(&mol, true, true),
            expected.tpsa_include_sandp,
            expected_bits.tpsa_include_sandp.as_str(),
            &mut failures,
            &mut by_field,
        );
        check_u32(
            "num_aromatic_rings",
            &context,
            calc_num_aromatic_rings(&mol),
            expected.num_aromatic_rings,
            &mut failures,
            &mut by_field,
        );
        check_u32(
            "num_rotatable_bonds_default",
            &context,
            calc_num_rotatable_bonds(&mol, NumRotatableBondsOptions::Default),
            expected.num_rotatable_bonds_default,
            &mut failures,
            &mut by_field,
        );
        check_u32(
            "num_rotatable_bonds_non_strict",
            &context,
            calc_num_rotatable_bonds(&mol, NumRotatableBondsOptions::NonStrict),
            expected.num_rotatable_bonds_non_strict,
            &mut failures,
            &mut by_field,
        );
        check_u32(
            "num_rotatable_bonds_strict",
            &context,
            calc_num_rotatable_bonds(&mol, NumRotatableBondsOptions::Strict),
            expected.num_rotatable_bonds_strict,
            &mut failures,
            &mut by_field,
        );
        check_u32(
            "num_rotatable_bonds_strict_linkages",
            &context,
            calc_num_rotatable_bonds(&mol, NumRotatableBondsOptions::StrictLinkages),
            expected.num_rotatable_bonds_strict_linkages,
            &mut failures,
            &mut by_field,
        );
        check_f64(
            "qed",
            &context,
            calc_qed(&mol),
            expected.qed,
            expected_bits.qed.as_str(),
            &mut failures,
            &mut by_field,
        );
    }

    if !failures.is_empty() {
        let field_summary = by_field
            .iter()
            .map(|(field, count)| format!("{field}: {count}"))
            .collect::<Vec<_>>()
            .join(", ");
        let sample = failures
            .iter()
            .take(24)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n");
        panic!(
            "molecular descriptor strict parity failed with {} field failures; by field: {field_summary}\nfirst failures:\n{sample}",
            failures.len()
        );
    }
}

#[test]
fn rotatable_descriptor_modes_match_rdkit_golden() {
    let records = load_golden();
    let mut shown = 0usize;
    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let expected = record.descriptors.as_ref().unwrap();
        let mol = Molecule::from_smiles(&record.smiles)
            .unwrap_or_else(|err| panic!("COSMolKit failed to parse row {}: {err}", row_idx + 1));
        let actual_default = calc_num_rotatable_bonds(&mol, NumRotatableBondsOptions::Default);
        let actual_non_strict = calc_num_rotatable_bonds(&mol, NumRotatableBondsOptions::NonStrict);
        let actual_strict = calc_num_rotatable_bonds(&mol, NumRotatableBondsOptions::Strict);
        if actual_default != Ok(expected.num_rotatable_bonds_default)
            || actual_non_strict != Ok(expected.num_rotatable_bonds_non_strict)
            || actual_strict != Ok(expected.num_rotatable_bonds_strict)
        {
            eprintln!(
                "row {} {} default {:?}/{} non {:?}/{} strict {:?}/{}",
                row_idx + 1,
                record.smiles,
                actual_default,
                expected.num_rotatable_bonds_default,
                actual_non_strict,
                expected.num_rotatable_bonds_non_strict,
                actual_strict,
                expected.num_rotatable_bonds_strict
            );
            shown += 1;
            if shown >= 24 {
                break;
            }
        }
    }
    assert_eq!(shown, 0, "rotatable descriptor differences were printed");
}

fn check_f64(
    field: &'static str,
    context: &str,
    actual: DescriptorResult<f64>,
    expected: f64,
    expected_bits: &str,
    failures: &mut Vec<String>,
    by_field: &mut BTreeMap<&'static str, usize>,
) {
    match actual {
        Ok(actual) => compare_f64(
            field,
            context,
            actual,
            expected,
            expected_bits,
            failures,
            by_field,
        ),
        Err(err) => record_error(field, context, err, failures, by_field),
    }
}

fn check_u32(
    field: &'static str,
    context: &str,
    actual: DescriptorResult<u32>,
    expected: u32,
    failures: &mut Vec<String>,
    by_field: &mut BTreeMap<&'static str, usize>,
) {
    match actual {
        Ok(actual) if actual == expected => {}
        Ok(actual) => record_failure(
            field,
            format!("{context}: {field} mismatch actual={actual} expected={expected}"),
            failures,
            by_field,
        ),
        Err(err) => record_error(field, context, err, failures, by_field),
    }
}

fn check_string(
    field: &'static str,
    context: &str,
    actual: DescriptorResult<String>,
    expected: &str,
    failures: &mut Vec<String>,
    by_field: &mut BTreeMap<&'static str, usize>,
) {
    match actual {
        Ok(actual) if actual == expected => {}
        Ok(actual) => record_failure(
            field,
            format!("{context}: {field} mismatch actual={actual:?} expected={expected:?}"),
            failures,
            by_field,
        ),
        Err(err) => record_error(field, context, err, failures, by_field),
    }
}

fn compare_f64(
    field: &'static str,
    context: &str,
    actual: f64,
    expected: f64,
    expected_bits: &str,
    failures: &mut Vec<String>,
    by_field: &mut BTreeMap<&'static str, usize>,
) {
    let expected_bits = parse_expected_f64_bits(field, context, expected_bits);
    if actual.to_bits() == expected_bits {
        return;
    }
    record_failure(
        field,
        format!(
            "{context}: {field} mismatch actual={actual:?} bits={:#018x} expected_json={expected:?} expected_bits={:#018x}",
            actual.to_bits(),
            expected_bits
        ),
        failures,
        by_field,
    );
}

fn parse_expected_f64_bits(field: &'static str, context: &str, bits: &str) -> u64 {
    u64::from_str_radix(bits, 16).unwrap_or_else(|err| {
        panic!("{context}: invalid {field} descriptor_bits value {bits:?}: {err}")
    })
}

fn record_error(
    field: &'static str,
    context: &str,
    err: DescriptorError,
    failures: &mut Vec<String>,
    by_field: &mut BTreeMap<&'static str, usize>,
) {
    record_failure(
        field,
        format!("{context}: {field} unsupported or failed closed: {err}"),
        failures,
        by_field,
    );
}

fn record_failure(
    field: &'static str,
    message: String,
    failures: &mut Vec<String>,
    by_field: &mut BTreeMap<&'static str, usize>,
) {
    *by_field.entry(field).or_default() += 1;
    failures.push(message);
}
