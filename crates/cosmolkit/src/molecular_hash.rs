//! Thin read-only projection of the sole native hashing owner.
#[cfg(all(test, feature = "cap-smiles", feature = "cap-rings"))]
mod scaffold_tests {
    use crate::*;
    type Transform = fn(&Molecule) -> Result<Molecule, OperationError>;
    const TRANSFORMS: [Transform; 3] = [
        Molecule::murcko_scaffold,
        Molecule::net_scaffold,
        Molecule::murcko_decompose,
    ];
    #[test]
    fn scaffold_rdkit_2026_03_6_regressions_and_source_cow() {
        // Fixed literal outputs verified with pinned rdMolHash.MolHash and
        // Chem.MurckoDecompose, not expectations inferred from CK.
        for (input, expected) in [
            ("", ["", "", ""]),
            ("CCO", ["", "", ""]),
            ("Cc1ccccc1", ["c1ccccc1", "*c1ccccc1", "c1ccccc1"]),
            ("O=C1CCCCC1", ["C1CCCCC1", "*=C1CCCCC1", "O=C1CCCCC1"]),
            ("Cn1cccc1", ["c1cc[nH]c1", "*n1cccc1", "c1cc[nH]c1"]),
            (
                "C1CCCCC1CCC2CCCCC2",
                [
                    "C1CCC(CCC2CCCCC2)CC1",
                    "C1CCC(CCC2CCCCC2)CC1",
                    "C1CCC(CCC2CCCCC2)CC1",
                ],
            ),
            ("C[C@H]1CCCCO1", ["C1CCOCC1", "*[C@H]1CCCCO1", "C1CCOCC1"]),
            ("CC1=CC=CC=C1.Cl", ["c1ccccc1", "*c1ccccc1", "c1ccccc1"]),
            (
                "[13CH3:7]C1CCCCC1",
                ["C1CCCCC1", "C1CCC([*:7])CC1", "C1CCCCC1"],
            ),
        ] {
            let source = Molecule::from_smiles(input).unwrap();
            let peer = source.clone();
            let original_cache = source.derived_cache_arc_runtime();
            for (transform, want) in TRANSFORMS.into_iter().zip(expected) {
                let output = transform(&source).unwrap();
                assert_eq!(
                    output.to_smiles().unwrap().as_bytes(),
                    want.as_bytes(),
                    "{input}"
                );
                assert_eq!(source, peer, "non-mutating scaffold changed source {input}");
                assert!(std::ptr::eq(source.topology(), peer.topology()));
                assert!(std::sync::Arc::ptr_eq(
                    &original_cache,
                    &source.derived_cache_arc_runtime()
                ));
            }
        }
        // MolHash's force=true stereo assignment only computes Fast rings
        // when the source ring state is absent/insufficient. No deletion must
        // retain the existing SymmSSSR type, not downgrade it to Fast rings.
        let source = Molecule::from_smiles("c1ccccc1").unwrap();
        let ring_type = source
            .derived_cache_runtime()
            .valid_ring_info()
            .unwrap()
            .find_type();
        for transform in [
            Molecule::murcko_scaffold as Transform,
            Molecule::net_scaffold,
        ] {
            let output = transform(&source).unwrap();
            assert_eq!(
                output
                    .derived_cache_runtime()
                    .valid_ring_info()
                    .unwrap()
                    .find_type(),
                ring_type
            );
        }
    }
    #[test]
    fn scaffold_mapping_preserves_all_conformers_and_ordinary_metadata() {
        let parsed = Molecule::from_smiles("CCOc1ccccc1").unwrap();
        let mut graph = parsed.topology().clone();
        for (i, atom) in graph.atoms.iter_mut().enumerate() {
            atom.set_atom_map(Some(i as u32 + 1));
            atom.set_prop("label", format!("original-{i}")).unwrap();
        }
        let n = graph.atoms.len();
        let coordinates = CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(5, (0..n).map(|i| [i as f64, -0.0]).collect())
                    .with_prop("kind", "2d"),
            ],
            conformers_3d: vec![
                Conformer3D::new(7, (0..n).map(|i| [i as f64, -0.0, 3.0]).collect(), true)
                    .with_prop("kind", "first"),
                Conformer3D::new(11, (0..n).map(|i| [i as f64, 5.0, -0.0]).collect(), false)
                    .with_prop("kind", "second"),
            ],
            source_coordinate_dim: None,
            source_conformer_order: Some(vec![
                CoordinateDimension::TwoD,
                CoordinateDimension::ThreeD,
                CoordinateDimension::ThreeD,
            ]),
        };
        let properties = MoleculeProperties::default()
            .with_prop("kept", "ordinary")
            .unwrap()
            .with_computed_prop("discarded", 9_i32)
            .unwrap();
        // from_parts intentionally has no prepared runtime valence cache;
        // prepare it explicitly, just as native MolHash requires cached Hs.
        let source = Molecule::from_parts(graph, coordinates, properties)
            .unwrap()
            .with_assigned_valence()
            .unwrap()
            .with_assigned_symm_sssr()
            .unwrap();
        let before = source.clone();
        for transform in TRANSFORMS {
            let output = transform(&source).unwrap();
            assert_eq!(output.property("kept"), source.property("kept"));
            assert!(output.property("discarded").is_none());
            assert_eq!(output.conformers_3d().len(), 2);
            let out_coords = output.coordinate_block_runtime();
            assert_eq!(
                out_coords.source_conformer_order,
                source.coordinate_block_runtime().source_conformer_order
            );
            for (new, atom) in output.atoms().iter().enumerate() {
                let old = atom.atom_map().unwrap() as usize - 1;
                assert_eq!(atom.prop("label"), source.atoms()[old].prop("label"));
                for (actual, original) in out_coords
                    .conformers_2d
                    .iter()
                    .zip(&source.coordinate_block_runtime().conformers_2d)
                {
                    assert_eq!(actual.id(), original.id());
                    assert_eq!(actual.props(), original.props());
                    assert_eq!(
                        actual.coordinates()[new].map(f64::to_bits),
                        original.coordinates()[old].map(f64::to_bits)
                    );
                }
                for (actual, original) in
                    out_coords.conformers_3d.iter().zip(source.conformers_3d())
                {
                    assert_eq!(actual.id(), original.id());
                    assert_eq!(actual.is_3d(), original.is_3d());
                    assert_eq!(actual.props(), original.props());
                    assert_eq!(
                        actual.coordinates()[new].map(f64::to_bits),
                        original.coordinates()[old].map(f64::to_bits)
                    );
                }
            }
            assert_eq!(source, before);
        }
    }
    #[test]
    fn scaffold_errors_preserve_source_and_typed_cause() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        let source = builder.build().unwrap();
        let before = source.clone();
        for transform in [
            Molecule::murcko_scaffold as Transform,
            Molecule::net_scaffold,
        ] {
            assert!(matches!(transform(&source), Err(OperationError::Scaffold(
                cosmolkit_core::ScaffoldError::MissingImplicitHydrogens { atom }
            )) if atom == AtomId::new(0)));
            assert_eq!(source, before);
        }
        assert!(matches!(
            source.murcko_decompose(),
            Err(OperationError::Scaffold(
                cosmolkit_core::ScaffoldError::MissingRingInfo
            ))
        ));
        assert_eq!(source, before);
    }
}
use crate::{Molecule, MoleculeHashError};
impl Molecule {
    /// Compute the original native 64-bit hash in legacy CIP rank order.
    /// Requires already prepared valence; absent ordinary rings remain absent.
    pub fn molecular_hash(&self) -> Result<u64, MoleculeHashError> {
        let cache = self.derived_cache_runtime();
        cosmolkit_fingerprints::molecule_hash(
            self.topology(),
            cache.valence_assignment(),
            cache.valid_ring_info(),
        )
    }
    /// Compute the same native hash with exactly one supplied rank per atom.
    pub fn molecular_hash_with_ranks(&self, ranks: &[u32]) -> Result<u64, MoleculeHashError> {
        cosmolkit_fingerprints::molecule_hash_with_ranks(
            self.topology(),
            self.derived_cache_runtime().valid_ring_info(),
            ranks,
        )
    }
}

#[cfg(all(test, feature = "cap-smiles"))]
mod tests {
    use super::*;
    fn sanitized_mol(smiles: &str) -> Molecule {
        // Canonical construction already runs sanitization and retains rings.
        Molecule::from_smiles(smiles).expect("valid SMILES")
    }
    #[test]
    fn test_mol_hash_benzene() {
        let mol = sanitized_mol("c1ccccc1");
        let hash = mol.molecular_hash().expect("hash");
        assert!(hash != 0);
    }

    #[test]
    fn test_mol_hash_deterministic() {
        let mol = sanitized_mol("c1ccccc1");
        let hash1 = mol.molecular_hash().expect("hash1");
        let hash2 = mol.molecular_hash().expect("hash2");
        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_mol_hash_different_molecules() {
        let mol1 = sanitized_mol("c1ccccc1");
        let mol2 = sanitized_mol("CCO");
        let hash1 = mol1.molecular_hash().expect("hash1");
        let hash2 = mol2.molecular_hash().expect("hash2");
        assert_ne!(hash1, hash2);
    }

    #[test]
    fn test_mol_hash_empty_error() {
        let mol = Molecule::new();
        assert!(mol.molecular_hash().is_err());
    }

    #[test]
    fn test_mol_hash_with_ranks() {
        let mol = sanitized_mol("c1ccccc1");
        let ranks = cosmolkit_core::assign_atom_cip_ranks(
            mol.topology(),
            mol.derived_cache_runtime().valence_assignment().unwrap(),
        )
        .expect("ranks");
        let hash = mol.molecular_hash_with_ranks(&ranks).expect("hash");
        assert!(hash != 0);
    }

    #[test]
    fn rank_length_and_missing_prepared_valence_are_structured() {
        let mol = sanitized_mol("CCO");
        for actual in [0, 2, 4] {
            assert!(
                matches!(mol.molecular_hash_with_ranks(&vec![0; actual]), Err(MoleculeHashError::RankCount {actual: n, atom_count: 3}) if n == actual)
            );
        }
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let raw = builder.build().unwrap();
        assert!(matches!(
            raw.molecular_hash(),
            Err(MoleculeHashError::MissingPreparedValence)
        ));
        assert!(raw.molecular_hash_with_ranks(&[0]).is_ok());
    }
}
