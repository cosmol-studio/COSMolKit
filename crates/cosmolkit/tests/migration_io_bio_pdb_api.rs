#![cfg(feature = "cap-bio")]

use cosmolkit::{
    BINDING_CONTRACT, BioCoordinateFormat, BioPdbReadParams, BioStructure, FunctionStatus, Protein,
    ProteinReadError,
};
use std::error::Error;

const ATOM: &str =
    "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n";
const WATER: &str =
    "HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 20.00           O  \n";

#[test]
fn public_pdb_preserves_structure_and_explicitly_projects_protein() {
    let text = format!("{ATOM}{WATER}END\n");
    let structure = BioStructure::from_pdb(&text).unwrap();
    structure.validate().unwrap();
    assert_eq!(structure.input_format(), BioCoordinateFormat::Pdb);
    assert_eq!(structure.atoms().len(), 2);
    assert_eq!(structure.atoms()[0].name().as_bytes(), b" CA ");
    assert_eq!(
        structure.coordinates().positions(),
        &[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    );
    let protein = Protein::from_pdb(&text).unwrap();
    assert_eq!(protein, structure.protein().unwrap());
    assert_eq!(protein.num_atoms(), 1);
    assert_eq!(protein.atoms().next().unwrap().name().as_str(), " CA ");
    assert_eq!(protein.atoms().next().unwrap().position(), [1.0, 2.0, 3.0]);
}

#[test]
fn public_pdb_defaults_and_explicit_options_are_forwarded() {
    let text = format!("{ATOM}TER\n{WATER}END\n");
    let defaults = BioPdbReadParams::default();
    assert_eq!(defaults.max_line_length, 0);
    assert!(
        !defaults.check_non_ascii
            && !defaults.ignore_ter
            && !defaults.split_chain_on_ter
            && !defaults.skip_remarks
    );
    assert_eq!(
        BioStructure::from_pdb(&text).unwrap(),
        BioStructure::from_pdb_with_params(&text, &defaults).unwrap()
    );
    assert_eq!(
        Protein::from_pdb(&text).unwrap(),
        Protein::from_pdb_with_params(&text, &defaults).unwrap()
    );
    for params in [
        BioPdbReadParams {
            max_line_length: 80,
            ..defaults
        },
        BioPdbReadParams {
            check_non_ascii: true,
            ..defaults
        },
        BioPdbReadParams {
            ignore_ter: true,
            ..defaults
        },
        BioPdbReadParams {
            split_chain_on_ter: true,
            ..defaults
        },
        BioPdbReadParams {
            skip_remarks: true,
            ..defaults
        },
    ] {
        // A local forwarding regression, not a second parser or corpus oracle.
        let expected = BioStructure::from_parts(
            cosmolkit_io::read_pdb_bio_structure(&text, "<string>", &params)
                .unwrap()
                .into_parts(),
        )
        .unwrap();
        assert_eq!(
            BioStructure::from_pdb_with_params(&text, &params).unwrap(),
            expected
        );
        assert_eq!(
            Protein::from_pdb_with_params(&text, &params).unwrap(),
            expected.protein().unwrap()
        );
    }
}

#[test]
fn public_pdb_errors_preserve_stage_line_record_and_source() {
    // Unknown records are skipped by the pinned dispatcher; retain this
    // first attempted error fixture as an explicit successful-empty regression.
    assert!(
        BioStructure::from_pdb("{\"not\":\"pdb\"}")
            .unwrap()
            .atoms()
            .is_empty()
    );
    // Gemmi pdb.cpp rejects an ATOM record shorter than the coordinate fields.
    let bad = "ATOM  \n";
    let direct =
        cosmolkit_io::read_pdb_bio_structure(bad, "<string>", &BioPdbReadParams::default())
            .unwrap_err();
    let error = BioStructure::from_pdb(bad).unwrap_err();
    assert_eq!(error.stage(), direct.stage());
    assert_eq!(error.line_number(), direct.line_number());
    assert_eq!(error.record_tag(), direct.record_tag());
    assert_eq!(error.to_string(), direct.to_string());
    assert!(error.source().is_some());
    let wrapped = Protein::from_pdb(bad).unwrap_err();
    assert!(matches!(&wrapped, ProteinReadError::Pdb(_)));
    assert!(wrapped.source().is_some());
}

#[test]
fn public_pdb_registry_has_real_experimental_signatures() {
    for target in ["BioStructure", "Protein"] {
        for name in ["from_pdb", "from_pdb_with_params"] {
            let id = format!("{target}.{name}");
            let entries: Vec<_> = BINDING_CONTRACT
                .iter()
                .filter(|e| e.semantic_id == id)
                .collect();
            assert_eq!(entries.len(), 1);
            let entry = entries[0];
            assert_eq!(entry.status, FunctionStatus::Experimental);
            assert_eq!(entry.feature, "cap-bio");
            assert_eq!(entry.python_name, name);
            let callable = entry.callable.unwrap();
            assert_eq!(callable.operation_semantic_id, None);
            assert_eq!(
                callable.parameters.len(),
                if name.ends_with("with_params") { 2 } else { 1 }
            );
            assert!(callable.error_type.is_some());
        }
    }
}
