//! Independent libtest functions; Cargo parallelizes them by default.
use cosmolkit_parity_tests::testing::run_registered;

macro_rules! tests {
    ($($key:ident),+ $(,)?) => {
        const KEYS: &[&str] = &[$(stringify!($key)),+];
        $(#[test] fn $key() { run_registered(stringify!($key), KEYS).unwrap(); })+
    };
}

tests!(
    bio_pdb_output_pdb,
    bio_pdb_output_cif,
    fuzzy_and_fingerprint_pairs,
    fuzzy_or_fingerprint_pairs,
    smiles_read_smiles,
    sanitize_smiles,
    kekulize_smiles,
    molecular_weight_smiles,
    exact_molecular_weight_smiles,
    molecular_formula_smiles,
    num_heavy_atoms_smiles,
    total_atom_count_smiles,
    lipinski_hba_smiles,
    lipinski_hbd_smiles,
    fraction_csp3_smiles,
    num_heteroatoms_smiles,
    num_hba_smiles,
    num_hbd_smiles,
    add_hydrogens_smiles,
    remove_hydrogens_smiles,
    coordinates_2d_smiles,
    svg_smiles,
    distance_matrix_smiles,
    num_rings_smiles,
    num_heterocycles_smiles,
    num_aromatic_rings_smiles,
    num_saturated_rings_smiles,
    num_aliphatic_rings_smiles,
    num_aromatic_heterocycles_smiles,
    num_aromatic_carbocycles_smiles,
    num_aliphatic_heterocycles_smiles,
    num_aliphatic_carbocycles_smiles,
    num_saturated_heterocycles_smiles,
    num_saturated_carbocycles_smiles,
    uff_has_all_molecule_params_smiles,
    uff_optimize_smiles,
    uff_optimize_conformers_smiles,
    morgan_fingerprint_smiles,
    morgan_sparse_fingerprint_smiles,
    morgan_count_fingerprint_smiles,
    morgan_sparse_count_fingerprint_smiles,
    chi_0_smiles,
    chi_1_smiles,
    hall_kier_alpha_smiles,
    hall_kier_alpha_with_contributions_smiles,
    kappa_1_smiles,
    kappa_2_smiles,
    kappa_3_smiles,
    phi_smiles,
    mqns_smiles,
    chi_0_v_smiles,
    chi_1_v_smiles,
    chi_2_v_smiles,
    chi_3_v_smiles,
    chi_4_v_smiles,
    chi_0_n_smiles,
    chi_1_n_smiles,
    chi_2_n_smiles,
    chi_3_n_smiles,
    chi_4_n_smiles,
    chi_n_v_smiles,
    chi_n_n_smiles,
);

// Preserve the existing BIO smoke checks alongside the corpus reference entrypoints.
#[test]
fn bio_pdb_output_pdb_registered_reference() {
    let tasks = cosmolkit_parity_tests::registry::select(Some("bio_pdb_output_pdb"))
        .expect("task is registered");
    assert_eq!(tasks.len(), 1);
    let task = tasks[0];
    assert_eq!(
        task.corpus_type,
        cosmolkit_parity_tests::registry::CorpusType::Pdb
    );
    assert_eq!(
        task.operation,
        cosmolkit_parity_tests::registry::Operation::BioPdbOutput
    );
    // Native Gemmi reference identity is bound to the task's generator.
    assert!(task.generator.starts_with("generate_bio_pdb_output"));

    let input = cosmolkit_parity_tests::registry::Input::BioPdbOutput {
        case: cosmolkit_parity_tests::registry::BioPdbCase {
            id: "ref_pdb_01".to_string(),
            text: "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nTER       2      ALA A   1                                                      \nEND".to_string(),
            format: cosmolkit_parity_tests::registry::BioPdbCorpusFormat::Pdb,
        },
        profile: cosmolkit_parity_tests::registry::BioPdbOutputProfile::ALL[31],
    };

    let record = cosmolkit_parity_tests::execute::run(&input)
        .expect("framework executes registered PDB reference");
    match &record.output {
        cosmolkit_parity_tests::registry::Value::BioPdbOutput(v) => {
            assert!(v.error.is_none(), "no error for valid PDB: {:?}", v.error);
            assert!(!v.text.is_empty(), "nonempty output");
            assert!(v.text.contains("ATOM"));
        }
        other => panic!("expected BioPdbOutput, got {other:?}"),
    }
    assert!(
        task.validate_reference(&input, &input, &record.output)
            .is_ok()
    );
}

#[test]
fn bio_pdb_output_cif_registered_reference() {
    let tasks = cosmolkit_parity_tests::registry::select(Some("bio_pdb_output_cif"))
        .expect("task is registered");
    assert_eq!(tasks.len(), 1);
    let task = tasks[0];
    assert_eq!(
        task.corpus_type,
        cosmolkit_parity_tests::registry::CorpusType::Cif
    );

    // Use actual mmCIF text for the CIF corpus (not PDB substitution).
    // A minimal mmCIF with one ATOM_SITE record.
    let mmcif_text = r#"data_test
_entry.id test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA A 1 1.000 2.000 3.000 1.00 20.00
"#;

    let input = cosmolkit_parity_tests::registry::Input::BioPdbOutput {
        case: cosmolkit_parity_tests::registry::BioPdbCase {
            id: "ref_cif_01".to_string(),
            text: mmcif_text.to_string(),
            format: cosmolkit_parity_tests::registry::BioPdbCorpusFormat::Cif,
        },
        profile: cosmolkit_parity_tests::registry::BioPdbOutputProfile::ALL[31],
    };

    let record = cosmolkit_parity_tests::execute::run(&input)
        .expect("framework executes registered CIF reference");
    match &record.output {
        cosmolkit_parity_tests::registry::Value::BioPdbOutput(v) => {
            // mmCIF parse may succeed (producing text) or produce a typed
            // Parse error — both are valid reference outcomes.
            match &v.error {
                Some(e) => {
                    assert!(
                        matches!(
                            e,
                            cosmolkit_parity_tests::registry::BioPdbOutputError::Parse { .. }
                        ),
                        "CIF parse error is stage-typed: {e:?}"
                    );
                }
                None => {
                    assert!(!v.text.is_empty(), "nonempty output when successful");
                }
            }
        }
        other => panic!("expected BioPdbOutput, got {other:?}"),
    }
    assert!(
        task.validate_reference(&input, &input, &record.output)
            .is_ok()
    );
}
