//! Rust-native cheminformatics and structural biology toolkit.
//!
//! COSMolKit provides molecular graph operations, SMILES/SMARTS processing,
//! molecular file IO, fingerprints and descriptors, 2D depiction, native 3D
//! conformer generation, UFF/MMFF optimization, InChI, substructure search,
//! batch workflows, and protein structure APIs.
//!
//! For supported cheminformatics operations, RDKit-compatible behavior is
//! treated as the correctness floor and validated through explicit parity
//! tests.
//!
//! This crate is the public Rust facade that re-exports COSMolKit core modules.
//!
//! # Examples
//!
//! ```
//! use cosmolkit::Molecule;
//!
//! let mol = Molecule::from_smiles("CCO").unwrap();
//! let with_hydrogens = mol.with_hydrogens().unwrap();
//!
//! assert_eq!(mol.num_atoms(), 3);
//! assert!(with_hydrogens.num_atoms() > mol.num_atoms());
//! ```
//!
//! ```
//! use cosmolkit::Molecule;
//!
//! let mut mol = Molecule::from_smiles("CCO").unwrap();
//! mol.add_hydrogens_().unwrap();
//! mol.compute_2d_coordinates_().unwrap();
//!
//! assert!(mol.num_atoms() > 3);
//! assert_eq!(mol.coordinates_2d().unwrap().len(), mol.num_atoms());
//! ```

pub use cosmolkit_core as core;
pub use cosmolkit_core::bio;
pub use cosmolkit_core::io;
pub use cosmolkit_core::*;

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn methane() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        builder.build().expect("trusted methane")
    }

    fn assert_generation_type(_: &MolToInchiOutput) {}
    fn assert_molecule_key_type(_: &MolToInchiKeyOutput) {}
    fn assert_direct_key_type(_: &InchiToInchiKeyOutput) {}
    fn assert_parse_type(_: &MolFromInchiOutput) {}
    fn assert_error_type(_: &InchiError) {}

    #[test]
    fn inchi_rust_facade_four_scalar_apis_match_exact_methane_results_and_preserve_source() {
        let source = methane()
            .with_name("source methane")
            .with_prop("record", "kept");
        let before = source.clone();

        let generated = mol_to_inchi(&source, None).expect("methane InChI");
        assert_generation_type(&generated);
        assert_eq!(generated.inchi, b"InChI=1S/CH4/h1H4");
        assert_eq!(generated.return_values.return_code, 0);
        assert!(generated.return_values.message.is_empty());
        assert_eq!(
            generated.return_values.log,
            concat!(
                "Generating standard InChI\n",
                "Input format: MOLfile\n",
                "Output format: Plain text\n",
                "Full Aux. info\n",
                "No timeout\n",
                "Up to 1024 atoms per structure"
            )
            .as_bytes()
        );
        assert_eq!(
            generated.return_values.aux_info,
            b"AuxInfo=1/0/N:1/rA:1C/rB:/rC:;"
        );
        assert!(generated.diagnostics.is_empty());

        let molecule_key = mol_to_inchi_key(&source, None).expect("methane InChIKey");
        assert_molecule_key_type(&molecule_key);
        assert_eq!(molecule_key.key, b"VNWKTOKETHGBQD-UHFFFAOYSA-N");
        assert!(molecule_key.diagnostics.is_empty());

        let direct_key = inchi_to_inchi_key(&generated.inchi).expect("direct methane InChIKey");
        assert_direct_key_type(&direct_key);
        assert_eq!(direct_key.key, molecule_key.key);
        assert!(direct_key.diagnostics.is_empty());

        let parsed = mol_from_inchi(&generated.inchi, false, false).expect("parsed methane");
        assert_parse_type(&parsed);
        assert_eq!(parsed.return_values.return_code, 0);
        assert!(parsed.diagnostics.is_empty());
        let parsed_molecule = parsed.molecule.expect("source returned methane graph");
        assert_eq!(parsed_molecule.num_atoms(), 1);
        assert_eq!(parsed_molecule.num_bonds(), 0);
        assert_eq!(parsed_molecule.atoms()[0].atomic_number(), 6);
        assert_eq!(parsed_molecule.atoms()[0].explicit_hydrogens(), 4);

        assert_eq!(source, before);
        assert_eq!(source.properties().name(), Some("source methane"));
        assert_eq!(source.prop("record"), Some("kept"));
    }

    #[test]
    fn inchi_rust_facade_preserves_source_defined_stereo_isotope_and_diagnostics() {
        let parsed = mol_from_inchi(b"InChI=1S/CHBrClF/c2-1(3)4/t1-/m0/s1/i1+1", false, false)
            .expect("parsed isotope and stereo fixture");
        assert_eq!(parsed.return_values.return_code, 1);
        let molecule = parsed.molecule.expect("source returned a graph");
        assert_eq!(molecule.num_atoms(), 4);
        assert_eq!(molecule.num_bonds(), 3);
        assert_eq!(molecule.atoms()[0].atomic_number(), 6);
        assert_eq!(molecule.atoms()[0].isotope(), Some(13));
        assert!(matches!(
            molecule.atoms()[0].chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ));

        let invalid_key = inchi_to_inchi_key(b"").expect("source-defined key diagnostic");
        assert!(invalid_key.key.is_empty());
        assert_eq!(invalid_key.diagnostics.len(), 1);
        assert_eq!(
            invalid_key.diagnostics[0].level,
            InchiDiagnosticLevel::Error
        );
        assert_eq!(
            invalid_key.diagnostics[0].message,
            "Invalid InChI prefix in generating InChI Key\n"
        );
    }

    #[test]
    fn inchi_rust_facade_preserves_structured_unsupported_and_allocation_error_categories() {
        let unsupported = mol_to_inchi(&Molecule::new(), None).expect_err("untrusted topology");
        assert_error_type(&unsupported);
        assert_eq!(unsupported.operation, "mol_to_inchi");
        assert_eq!(unsupported.kind, InchiErrorKind::UnsupportedState);
        assert!(unsupported.detail.contains("trusted"));

        let mut query_builder = MoleculeBuilder::new();
        query_builder.add_atom(
            AtomSpec::new(Element::C).with_query(QueryNode::predicate(AtomQueryPredicate::Any)),
        );
        let query = query_builder.build().expect("query graph");
        let unsupported = mol_to_inchi_key(&query, None).expect_err("query atom");
        assert_eq!(unsupported.operation, "mol_to_inchi_key");
        assert_eq!(unsupported.kind, InchiErrorKind::UnsupportedState);
        assert!(unsupported.detail.contains("query atoms"));

        // The source/core focused tests force the official-undefined allocation
        // path. This facade assertion ensures its public category is not erased.
        let allocation = InchiError {
            operation: "mol_to_inchi",
            kind: InchiErrorKind::AllocationFailed,
            detail: "AllocationFailed".to_owned(),
        };
        assert_error_type(&allocation);
        assert_eq!(allocation.kind, InchiErrorKind::AllocationFailed);
        assert_eq!(allocation.detail, "AllocationFailed");
    }

    #[test]
    fn conformer_generation_public_api_is_reexported_from_facade() {
        let molecule = Molecule::from_smiles("CC").expect("ethane");

        let mut params = EmbedParameters::etkdg_v3();
        params.random_seed = 0xf00d;
        params.num_threads = 1;
        let generated = molecule
            .with_3d_conformer_with_params(params.clone())
            .expect("value-style conformer generation");
        assert!(molecule.conformers_3d().is_empty());
        assert_eq!(generated.conformers_3d().len(), 1);

        let direct_result =
            embed_molecule_result(&molecule, &mut params).expect("embed_molecule_result re-export");
        assert!(direct_result.ok());
        assert_eq!(direct_result.conf_id, 0);
        assert_eq!(direct_result.molecule.conformers_3d().len(), 1);

        let (embedded, conf_id) =
            embed_molecule(&molecule, &mut params).expect("direct embed_molecule re-export");
        assert_eq!(conf_id, 0);
        assert_eq!(embedded.conformers_3d().len(), 1);

        let mut multi_params = EmbedParameters::etkdg();
        multi_params.random_seed = 0xf00d;
        multi_params.num_threads = 1;
        let generated_multi = molecule
            .with_3d_conformers_with_params(2, multi_params.clone())
            .expect("value-style multi-conformer generation");
        assert_eq!(generated_multi.conformers_3d().len(), 2);

        let multi_result = embed_multiple_confs_result(&molecule, 2, &mut multi_params)
            .expect("embed_multiple_confs_result re-export");
        assert_eq!(multi_result.conf_ids, vec![0, 1]);
        assert_eq!(multi_result.generated_count(), 2);
        assert_eq!(multi_result.molecule.conformers_3d().len(), 2);

        let (embedded_multi, ids) = embed_multiple_confs(&molecule, 2, &mut multi_params)
            .expect("direct embed_multiple_confs re-export");
        assert_eq!(ids, vec![0, 1]);
        assert_eq!(embedded_multi.conformers_3d().len(), 2);
    }
}
