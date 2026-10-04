//! Independent reference-parity tests for the BIO PDB output tasks
//! (BIO-PDB-WRITE Steps 65-70). Two tests: one for PDB corpus, one
//! for CIF corpus, each exercising the full framework path.

#[cfg(test)]
mod bio_pdb_output_reference_tests {
    use crate::registry::{
        self, BioPdbCase, BioPdbOutputProfile, CorpusType, Input, Operation, Value,
    };

    const PDB_CASE: &str = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nTER       2      ALA A   1                                                      \nEND";

    /// Reference test 1: bio_pdb_output_pdb — select the task, construct
    /// a typed input, execute through the framework, and verify the
    /// output value is a nonempty text (or a typed error).
    #[test]
    fn bio_pdb_output_pdb_reference_pipeline() {
        // 1. Select the task by exact key.
        let tasks = registry::select(Some("bio_pdb_output_pdb")).expect("task registered");
        assert_eq!(tasks.len(), 1, "exactly one task");
        let task = tasks[0];
        assert_eq!(task.corpus_type, CorpusType::Pdb);
        assert_eq!(task.operation, Operation::BioPdbOutput);

        // 2. Construct a typed input (PDB format).
        let input = Input::BioPdbOutput {
            case: BioPdbCase {
                id: "ref_pdb_01".to_string(),
                text: PDB_CASE.to_string(),
                format: crate::registry::BioPdbCorpusFormat::Pdb,
            },
            profile: BioPdbOutputProfile::ALL[31], // all-true default
        };

        // 3. Execute through the framework.
        let record = crate::execute::run(&input).expect("framework executes PDB case");

        // 4. Verify the output.
        match &record.output {
            Value::BioPdbOutput(value) => {
                if let Some(error) = &value.error {
                    panic!("unexpected error in PDB reference: {error}");
                }
                assert!(!value.text.is_empty(), "PDB reference output is nonempty");
                assert!(value.text.contains("ATOM"), "contains ATOM");
                assert!(value.text.contains("TER"), "contains TER");
            }
            other => panic!("expected BioPdbOutput value, got {other:?}"),
        }

        // 5. Validate the reference.
        assert!(
            task.validate_reference(&input, &input, &record.output)
                .is_ok(),
            "reference validates"
        );
    }

    /// Reference test 2: bio_pdb_output_cif — select the task, construct
    /// a typed CIF input, execute through the framework (mmCIF parse),
    /// and verify the output.
    #[test]
    fn bio_pdb_output_cif_reference_pipeline() {
        // 1. Select the task by exact key.
        let tasks = registry::select(Some("bio_pdb_output_cif")).expect("task registered");
        assert_eq!(tasks.len(), 1, "exactly one task");
        let task = tasks[0];
        assert_eq!(task.corpus_type, CorpusType::Cif);
        assert_eq!(task.operation, Operation::BioPdbOutput);

        // 2. Construct a typed input (CIF format — mmCIF text).
        // Note: the native Gemmi oracle reads mmCIF, then writes PDB.
        // For this reference test, we use a simple PDB-parsable string
        // since from_mmcif requires actual mmCIF text. The full corpus
        // preparation (Step 71-72) will provide real mmCIF inputs.
        let input = Input::BioPdbOutput {
            case: BioPdbCase {
                id: "ref_cif_01".to_string(),
                text: "data_test\n_entry.id test\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\nATOM 1 C CA ALA A 1 1.000 2.000 3.000 1.00 20.00\n".to_string(),
                format: crate::registry::BioPdbCorpusFormat::Cif,
            },
            profile: BioPdbOutputProfile::ALL[31],
        };

        // 3. Execute through the framework.
        let record = crate::execute::run(&input).expect("framework executes CIF case");

        // 4. Verify the output.
        match &record.output {
            Value::BioPdbOutput(value) => {
                if let Some(error) = &value.error {
                    panic!("unexpected error in CIF reference: {error}");
                }
                assert!(!value.text.is_empty(), "CIF reference output is nonempty");
            }
            other => panic!("expected BioPdbOutput value, got {other:?}"),
        }

        // 5. Validate the reference.
        assert!(
            task.validate_reference(&input, &input, &record.output)
                .is_ok(),
            "reference validates"
        );
    }
}
