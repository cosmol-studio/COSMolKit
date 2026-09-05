#![recursion_limit = "256"]

use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use cosmolkit_core::{
    Molecule, calc_chi_0, calc_chi_0n, calc_chi_0v, calc_chi_1, calc_chi_1n, calc_chi_1v, calc_chi_2n, calc_chi_2v,
    calc_chi_3n, calc_chi_3v, calc_chi_4n, calc_chi_4v, calc_chi_nn, calc_chi_nv, calc_fraction_csp3,
    calc_hall_kier_alpha, calc_hall_kier_alpha_with_contributions, calc_kappa_1, calc_kappa_2, calc_kappa_3,
    calc_labute_asa, calc_labute_asa_contributions, calc_lipinski_hba, calc_lipinski_hbd, calc_mqns,
    calc_num_aliphatic_carbocycles, calc_num_aliphatic_heterocycles, calc_num_aliphatic_rings, calc_num_amide_bonds,
    calc_num_aromatic_carbocycles, calc_num_aromatic_heterocycles, calc_num_aromatic_rings,
    calc_num_atom_stereo_centers, calc_num_atoms, calc_num_bridgehead_atoms, calc_num_hba, calc_num_hbd,
    calc_num_heavy_atoms, calc_num_heteroatoms, calc_num_heterocycles, calc_num_rings, calc_num_saturated_carbocycles,
    calc_num_saturated_heterocycles, calc_num_saturated_rings, calc_num_spiro_atoms,
    calc_num_unspecified_atom_stereo_centers, calc_phi, calc_slogp_vsa, calc_slogp_vsa_with_bins, calc_smr_vsa,
    calc_smr_vsa_with_bins,
};
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Debug, Deserialize)]
struct Record {
    smiles: String,
    rdkit_ok: bool,
    error: Option<String>,
    high_feasibility_descriptor_bits: Option<Value>,
    high_feasibility_contribution_bits: Option<Value>,
    high_feasibility_cache_profile_bits: Option<Value>,
}

fn bit(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn bits(values: impl IntoIterator<Item = f64>) -> Vec<String> {
    values.into_iter().map(bit).collect()
}

fn molecule(smiles: &str) -> Molecule {
    Molecule::from_smiles(smiles)
        .unwrap_or_else(|error| panic!("failed to parse focused descriptor row {smiles:?}: {error}"))
}

fn descriptor_bits(smiles: &str) -> Value {
    let mol = molecule(smiles);
    let chi_nv = (0..=6)
        .map(|order| bit(calc_chi_nv(&mol, order, false).unwrap()))
        .collect::<Vec<_>>();
    let chi_nn = (0..=6)
        .map(|order| bit(calc_chi_nn(&mol, order, false).unwrap()))
        .collect::<Vec<_>>();

    json!({
        "chi_0": bit(calc_chi_0(&mol)),
        "chi_1": bit(calc_chi_1(&mol)),
        "chi_nv_orders_0_6": chi_nv,
        "chi_nn_orders_0_6": chi_nn,
        "chi_0v": bit(calc_chi_0v(&mol, false).unwrap()),
        "chi_1v": bit(calc_chi_1v(&mol, false).unwrap()),
        "chi_2v": bit(calc_chi_2v(&mol, false).unwrap()),
        "chi_3v": bit(calc_chi_3v(&mol, false).unwrap()),
        "chi_4v": bit(calc_chi_4v(&mol, false).unwrap()),
        "chi_0n": bit(calc_chi_0n(&mol, false).unwrap()),
        "chi_1n": bit(calc_chi_1n(&mol, false).unwrap()),
        "chi_2n": bit(calc_chi_2n(&mol, false).unwrap()),
        "chi_3n": bit(calc_chi_3n(&mol, false).unwrap()),
        "chi_4n": bit(calc_chi_4n(&mol, false).unwrap()),
        "hall_kier_alpha": bit(calc_hall_kier_alpha(&mol)),
        "kappa_1": bit(calc_kappa_1(&mol)),
        "kappa_2": bit(calc_kappa_2(&mol).unwrap()),
        "kappa_3": bit(calc_kappa_3(&mol).unwrap()),
        "phi": bit(calc_phi(&mol).unwrap()),
        "lipinski_hba": calc_lipinski_hba(&mol).unwrap(),
        "lipinski_hbd": calc_lipinski_hbd(&mol).unwrap(),
        "num_hba": calc_num_hba(&mol).unwrap(),
        "num_hbd": calc_num_hbd(&mol).unwrap(),
        "num_heteroatoms": calc_num_heteroatoms(&mol).unwrap(),
        "num_amide_bonds": calc_num_amide_bonds(&mol).unwrap(),
        "num_heavy_atoms": calc_num_heavy_atoms(&mol).unwrap(),
        "num_atoms": calc_num_atoms(&mol).unwrap(),
        "num_rings": calc_num_rings(&mol).unwrap(),
        "num_heterocycles": calc_num_heterocycles(&mol).unwrap(),
        "num_aromatic_rings": calc_num_aromatic_rings(&mol).unwrap(),
        "num_saturated_rings": calc_num_saturated_rings(&mol).unwrap(),
        "num_aliphatic_rings": calc_num_aliphatic_rings(&mol).unwrap(),
        "num_aromatic_heterocycles": calc_num_aromatic_heterocycles(&mol).unwrap(),
        "num_aromatic_carbocycles": calc_num_aromatic_carbocycles(&mol).unwrap(),
        "num_aliphatic_heterocycles": calc_num_aliphatic_heterocycles(&mol).unwrap(),
        "num_aliphatic_carbocycles": calc_num_aliphatic_carbocycles(&mol).unwrap(),
        "num_saturated_heterocycles": calc_num_saturated_heterocycles(&mol).unwrap(),
        "num_saturated_carbocycles": calc_num_saturated_carbocycles(&mol).unwrap(),
        "num_spiro_atoms": calc_num_spiro_atoms(&mol).unwrap(),
        "num_bridgehead_atoms": calc_num_bridgehead_atoms(&mol).unwrap(),
        "num_atom_stereo_centers": calc_num_atom_stereo_centers(&mol).unwrap(),
        "num_unspecified_atom_stereo_centers": calc_num_unspecified_atom_stereo_centers(&mol).unwrap(),
        "fraction_csp3": bit(calc_fraction_csp3(&mol).unwrap()),
        "mqns": calc_mqns(&mol).unwrap().into_iter().collect::<Vec<_>>(),
        "labute_asa_include_hs_false": bit(calc_labute_asa(&mol, false, false)),
        "labute_asa_include_hs_true": bit(calc_labute_asa(&mol, true, true)),
        "slogp_vsa": bits(calc_slogp_vsa(&mol, false).unwrap()),
        "smr_vsa": bits(calc_smr_vsa(&mol, false).unwrap()),
    })
}

fn contribution_bits(smiles: &str) -> Value {
    let alpha = calc_hall_kier_alpha_with_contributions(&molecule(smiles));
    let labute_false = calc_labute_asa_contributions(&molecule(smiles), false, false);
    let labute_true = calc_labute_asa_contributions(&molecule(smiles), true, false);
    json!({
        "hall_kier_alpha": {
            "value": bit(alpha.alpha),
            "atom_contributions": bits(alpha.atom_contributions),
        },
        "labute_asa": {
            "include_hs_false": {
                "asa": bit(labute_false.asa),
                "atom_contributions": bits(labute_false.atom_contributions),
                "hydrogen_contribution": bit(labute_false.hydrogen_contribution),
            },
            "include_hs_true": {
                "asa": bit(labute_true.asa),
                "atom_contributions": bits(labute_true.atom_contributions),
                "hydrogen_contribution": bit(labute_true.hydrogen_contribution),
            },
        },
    })
}

fn chi_cache_profile(smiles: &str, valence: bool) -> Value {
    let mol = molecule(smiles);
    let calculate = |order, force| {
        if valence {
            calc_chi_nv(&mol, order, force)
        } else {
            calc_chi_nn(&mol, order, force)
        }
        .unwrap()
    };
    json!({
        "cold": bits((0..=6).map(|order| calculate(order, false))),
        "warm": bits((0..=6).map(|order| calculate(order, false))),
        "forced": bits((0..=6).map(|order| calculate(order, true))),
    })
}

fn cache_profile_bits(smiles: &str) -> Value {
    let labute_mol = molecule(smiles);
    let labute_calls = [(false, false), (true, false), (true, true), (false, false)]
        .into_iter()
        .map(|(include_hs, force)| {
            json!({
                "include_hs": include_hs,
                "force": force,
                "value": bit(calc_labute_asa(&labute_mol, include_hs, force)),
            })
        })
        .collect::<Vec<_>>();
    let custom_bins = [-0.2, 0.0, 0.25, 0.25, 0.8];
    let vsa_mol = molecule(smiles);
    let slogp_cold = calc_slogp_vsa(&vsa_mol, false).unwrap();
    let slogp_warm = calc_slogp_vsa(&vsa_mol, false).unwrap();
    let slogp_custom = calc_slogp_vsa_with_bins(&vsa_mol, &custom_bins, true).unwrap();
    let smr_warm = calc_smr_vsa(&vsa_mol, false).unwrap();
    let smr_custom = calc_smr_vsa_with_bins(&vsa_mol, &custom_bins, true).unwrap();

    json!({
        "chi_nv": chi_cache_profile(smiles, true),
        "chi_nn": chi_cache_profile(smiles, false),
        "labute_asa_sequence": labute_calls,
        "slogp_vsa_default_cold": bits(slogp_cold),
        "slogp_vsa_default_warm": bits(slogp_warm),
        "slogp_vsa_custom_forced": bits(slogp_custom),
        "smr_vsa_default_warm": bits(smr_warm),
        "smr_vsa_custom_forced": bits(smr_custom),
        "custom_bins": bits(custom_bins),
    })
}

fn load_profile(profile: &str) -> Vec<Record> {
    let path = cosmolkit_test_support::expected_path_for_profile(
        "descriptors",
        "rdkit",
        profile,
        "molecular_descriptors.jsonl",
    );
    BufReader::new(File::open(&path).unwrap())
        .lines()
        .enumerate()
        .map(|(index, line)| {
            serde_json::from_str(&line.unwrap())
                .unwrap_or_else(|error| panic!("failed to parse {} line {}: {error}", path.display(), index + 1))
        })
        .collect()
}

#[test]
fn focused_high_feasibility_descriptors_match_pinned_rdkit_exactly() {
    let records = load_profile("descriptors_focused");
    assert_eq!(records.len(), 69);
    compare_records(&records);
}

#[test]
fn active_corpus_high_feasibility_descriptors_match_pinned_rdkit_exactly() {
    let profile = cosmolkit_test_support::profile_name();
    let records = load_profile(&profile);
    assert_eq!(records.len(), cosmolkit_test_support::count_smiles_rows());
    compare_records(&records);
}

fn compare_records(records: &[Record]) {
    for (index, record) in records.iter().enumerate() {
        let context = format!("row {} ({:?})", index + 1, record.smiles);
        if !record.rdkit_ok {
            assert!(record.error.is_some(), "RDKit rejection lacks an error in {context}");
            assert!(
                Molecule::from_smiles(&record.smiles).is_err(),
                "RDKit rejected {context}, but COSMolKit accepted it"
            );
            continue;
        }
        assert_eq!(
            descriptor_bits(&record.smiles),
            *record.high_feasibility_descriptor_bits.as_ref().unwrap(),
            "descriptor mismatch in {context}"
        );
        assert_eq!(
            contribution_bits(&record.smiles),
            *record.high_feasibility_contribution_bits.as_ref().unwrap(),
            "contribution mismatch in {context}"
        );
        assert_eq!(
            cache_profile_bits(&record.smiles),
            *record.high_feasibility_cache_profile_bits.as_ref().unwrap(),
            "cache-profile mismatch in {context}"
        );
    }
}
