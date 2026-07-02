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
    error: Option<String>,
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
        let mol = Molecule::from_smiles(&record.smiles)
            .unwrap_or_else(|err| panic!("COSMolKit failed to parse {context}: {err}"));

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
