#![cfg(all(
    feature = "cap-smiles",
    feature = "cap-io",
    feature = "cap-descriptors"
))]

use cosmolkit::{MolBlockWriteParams, Molecule, PropertyText, SdfFormat};

#[test]
fn default_mol_and_sdf_writers_preserve_stereo_without_input_coordinates() {
    // Pinned MolFileWriter.cpp prepareMol generates a temporary 2D conformer
    // when includeStereo is true and the source has no conformers.
    for (smiles, expected) in [
        ("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"),
        ("N[C@H](C)C(=O)O", "C[C@@H](N)C(=O)O"),
        ("F/C=C/F", "F/C=C/F"),
        ("F/C=C\\F", "F/C=C\\F"),
        ("C[C@H](O)[C@@H](O)C", "C[C@H](O)[C@H](C)O"),
    ] {
        let expected = PropertyText::from(expected);
        let molecule = Molecule::from_smiles(smiles).unwrap();
        let before = molecule.clone();
        assert!(molecule.coordinates_2d().is_none());
        assert_eq!(molecule.num_3d_conformers(), 0);
        for format in [SdfFormat::V2000, SdfFormat::V3000] {
            let params = MolBlockWriteParams {
                format,
                ..Default::default()
            };
            let mol = molecule.to_mol_with_params(&params).unwrap();
            let sdf = molecule.to_sdf_with_params(&params).unwrap();
            assert_eq!(
                Molecule::from_mol(&mol).unwrap().to_smiles().unwrap(),
                expected,
                "{smiles}: MOL {format:?}"
            );
            assert_eq!(
                Molecule::from_sdf(&sdf).unwrap().to_smiles().unwrap(),
                expected,
                "{smiles}: SDF {format:?}"
            );
            assert_eq!(molecule, before, "writer must not change {smiles}");
        }
        assert_eq!(
            Molecule::from_mol(&molecule.to_mol().unwrap())
                .unwrap()
                .to_smiles()
                .unwrap(),
            expected
        );
        assert_eq!(
            Molecule::from_sdf(&molecule.to_sdf().unwrap())
                .unwrap()
                .to_smiles()
                .unwrap(),
            expected
        );
    }
}

#[test]
fn disabling_stereo_does_not_generate_writer_coordinates() {
    let molecule = Molecule::from_smiles("N[C@@H](C)C(=O)O").unwrap();
    let before = molecule.clone();
    let mol = molecule
        .to_mol_with_params(&MolBlockWriteParams {
            include_stereo: false,
            ..Default::default()
        })
        .unwrap();
    for row in mol.lines().skip(4).take(molecule.num_atoms()) {
        assert_eq!(&row[..30], "    0.0000    0.0000    0.0000");
    }
    assert_eq!(molecule, before);
    assert!(molecule.coordinates_2d().is_none());
}

#[test]
fn public_parser_rejects_invalid_bond_and_separator_sequences() {
    for smiles in ["C==C", "C#", "C..C"] {
        assert!(Molecule::from_smiles(smiles).is_err(), "{smiles}");
    }
}

#[test]
fn public_empty_molecule_qed_matches_pinned_rdkit() {
    let molecule = Molecule::from_smiles("").unwrap();
    assert_eq!(molecule.qed().unwrap().to_bits(), 0x3fd5b91db87826fc);
    assert_eq!(molecule.num_atoms(), 0);
    assert_eq!(molecule.to_smiles().unwrap(), PropertyText::from(""));
}
