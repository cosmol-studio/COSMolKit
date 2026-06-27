use cosmolkit_core::{
    ChiralTag, LigandRef, SdfReadParams,
    io::{
        mol2::{Mol2ReadParams, read_mol2_from_str_with_params},
        sdf::read_sdf_from_str_with_params,
    },
};

const ONE_AID_SDF: &str = include_str!("../../../tests/fixtures/regression/1aid_ligand.sdf");
const ONE_AID_MOL2: &str = include_str!("../../../tests/fixtures/regression/1aid_ligand.mol2");

fn assert_1aid_tetrahedral_stereo_is_preserved(molecule: &cosmolkit_core::Molecule, source: &str) {
    assert_eq!(
        molecule.atoms()[0].chiral_tag(),
        ChiralTag::TetrahedralCw,
        "{source} should preserve RDKit 3D chirality at atom 0 after reader removeHs"
    );
    assert_eq!(
        molecule.atoms()[3].chiral_tag(),
        ChiralTag::TetrahedralCw,
        "{source} should preserve RDKit 3D chirality at atom 3 after reader removeHs"
    );

    let stereo = molecule.tetrahedral_stereo().expect("tetrahedral stereo");
    let center0 = stereo
        .iter()
        .find(|entry| entry.center.index() == 0)
        .unwrap_or_else(|| panic!("{source} should report tetrahedral stereo at atom 0"));
    let center3 = stereo
        .iter()
        .find(|entry| entry.center.index() == 3)
        .unwrap_or_else(|| panic!("{source} should report tetrahedral stereo at atom 3"));
    assert!(
        center0.ligands.contains(&LigandRef::ImplicitHydrogen),
        "{source} atom 0 should retain the removed hydrogen as an implicit ligand"
    );
    assert!(
        !center3.ligands.contains(&LigandRef::ImplicitHydrogen),
        "{source} atom 3 should remain four-explicit-ligand tetrahedral stereo"
    );
}

#[test]
fn sdf_reader_preserves_1aid_3d_chirality_through_sanitize_remove_hs_like_rdkit() {
    let molecule = read_sdf_from_str_with_params(
        ONE_AID_SDF,
        SdfReadParams {
            sanitize: true,
            remove_hs: true,
            ..Default::default()
        },
    )
    .expect("1aid SDF should parse")
    .molecule;

    assert_1aid_tetrahedral_stereo_is_preserved(&molecule, "SDF");
}

#[test]
fn mol2_reader_preserves_1aid_3d_chirality_through_sanitize_remove_hs_like_rdkit() {
    let molecule = read_mol2_from_str_with_params(
        ONE_AID_MOL2,
        Mol2ReadParams {
            sanitize: true,
            remove_hs: true,
            ..Default::default()
        },
    )
    .expect("1aid MOL2 should parse")
    .expect("1aid MOL2 should produce a molecule")
    .molecule;

    assert_1aid_tetrahedral_stereo_is_preserved(&molecule, "MOL2");
}
