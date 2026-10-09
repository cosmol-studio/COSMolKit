#![cfg(feature = "cap-io")]
//! Source-extracted writer regressions; proposals pending independent review.
use cosmolkit::{MolBlockWriteParams, Molecule, SdfReadParams};
fn bond_rows(mol: &Molecule, text: &str) -> Vec<(usize, usize, usize)> {
    text.lines()
        .skip(4 + mol.num_atoms())
        .take(mol.num_bonds())
        .filter_map(|line| {
            let stereo = line[9..12].trim().parse::<usize>().unwrap();
            (stereo != 0).then(|| {
                (
                    line[..3].trim().parse::<usize>().unwrap() - 1,
                    line[3..6].trim().parse::<usize>().unwrap() - 1,
                    stereo,
                )
            })
        })
        .collect()
}
#[test]
fn io44_sdf_row196_wedges_use_final_symmetrized_ring_carrier() {
    let mol = Molecule::from_mol(include_str!(
        "../../../testdata/molblock/fixtures/io44-sdf-row196.mol"
    ))
    .unwrap();
    let peer = mol.clone();
    let before = mol
        .to_mol_with_params(&MolBlockWriteParams {
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        })
        .unwrap();
    let rows: Vec<_> = before
        .lines()
        .skip(4 + mol.num_atoms())
        .take(mol.num_bonds())
        .enumerate()
        .filter_map(|(i, l)| (l[9..12].trim() != "0").then_some(i))
        .collect();
    assert_eq!(rows, vec![19, 30, 34, 55]);
    assert_eq!(mol.atoms(), peer.atoms());
    assert_eq!(mol.bonds(), peer.bonds());
    assert_eq!(
        mol.to_mol_with_params(&MolBlockWriteParams {
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        })
        .unwrap(),
        before
    );
}
#[test]
fn io44_sdf_row3500_equal_score_source_sort_controls_wedge_centers() {
    let mol = Molecule::from_mol_with_params(
        include_str!("../../../testdata/molblock/fixtures/io44-sdf-row3500.mol"),
        &SdfReadParams {
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    let peer = mol.clone();
    let text = mol
        .to_mol_with_params(&MolBlockWriteParams {
            include_stereo: false,
            kekulize: false,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        bond_rows(&mol, &text),
        vec![
            (24, 1, 6),
            (26, 25, 6),
            (27, 26, 1),
            (28, 27, 1),
            (29, 28, 1),
            (30, 29, 6),
            (31, 32, 1),
            (25, 30, 6)
        ]
    );
    assert_eq!(mol.atoms(), peer.atoms());
    assert_eq!(mol.bonds(), peer.bonds());
}
