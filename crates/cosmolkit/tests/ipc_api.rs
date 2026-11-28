#![cfg(feature = "cap-descriptors")]
#[cfg(feature = "descriptors")]
use cosmolkit::{AtomId, IpcDistanceMatrix, IpcError, Molecule, SmilesParseParams};

#[test]
fn checked_rows_work_without_smiles_or_dense_matrix_support() {
    use cosmolkit::{AtomSpec, Element, MoleculeBuilder};

    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    let mol = builder.build().unwrap();
    let loop_row = [vec![1.0]];

    assert_eq!(mol.ipc().unwrap(), 0.0);
    assert_eq!(
        mol.ipc_with_distance_matrix(loop_row.as_slice().into())
            .unwrap(),
        2.0
    );
    assert_eq!(
        mol.avg_ipc_with_distance_matrix(loop_row.as_slice().into())
            .unwrap(),
        1.0
    );
}

#[cfg(feature = "descriptors")]
#[test]
fn default_graph_uses_all_vertices_and_needs_no_preparation() {
    for smiles in ["CC", "C=C", "[H]O"] {
        let mol = Molecule::from_smiles_with_params(
            smiles,
            &SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        let edited = mol
            .with_atoms_renumbered(&[AtomId::new(1), AtomId::new(0)])
            .unwrap();
        assert_eq!(
            (edited.ipc().unwrap(), edited.avg_ipc().unwrap()),
            (2.0, 1.0)
        );
    }
    let oxygen = Molecule::from_smiles("O").unwrap();
    assert_eq!(oxygen.ipc().unwrap(), 0.0);
    let two_edges = Molecule::from_smiles("CC.CC").unwrap();
    assert_eq!(
        (two_edges.ipc().unwrap(), two_edges.avg_ipc().unwrap()),
        (6.0, 1.5)
    );
}

#[cfg(feature = "descriptors")]
#[test]
fn path_information_uses_the_released_float_polynomial() {
    let mol = Molecule::from_smiles("CCCC").unwrap();
    assert_eq!(
        mol.ipc().unwrap().to_bits(),
        6.854752972273345_f64.to_bits()
    );
    assert_eq!(
        mol.avg_ipc().unwrap().to_bits(),
        1.370950594454669_f64.to_bits()
    );
}

#[cfg(feature = "descriptors")]
#[test]
fn borrowed_matrices_replace_the_current_graph_without_changing_it() {
    let mol = Molecule::from_smiles("CC").unwrap();
    let other = Molecule::from_smiles("C.C")
        .unwrap()
        .distance_matrix()
        .unwrap();
    let score = |matrix: IpcDistanceMatrix<'_>| {
        (
            mol.ipc_with_distance_matrix(matrix).unwrap(),
            mol.avg_ipc_with_distance_matrix(matrix).unwrap(),
        )
    };

    assert_eq!(mol.ipc().unwrap(), 2.0);
    assert_eq!(score((&other).into()), (0.0, 0.0));
    assert_eq!(
        score([vec![0.0, 1.0], vec![0.0, 0.0]].as_slice().into()),
        (0.0, 0.0)
    );
    assert_eq!(
        score([vec![0.0, f64::NAN], vec![1.0, 0.0]].as_slice().into()),
        (0.0, 0.0)
    );
    assert_eq!(score((&other).into()), (0.0, 0.0));
    assert_eq!(mol.ipc().unwrap(), 2.0);

    let carbon = Molecule::from_smiles("C").unwrap();
    assert_eq!(
        carbon
            .avg_ipc_with_distance_matrix([vec![1.0]].as_slice().into())
            .unwrap(),
        1.0
    );
    assert_eq!(
        carbon
            .avg_ipc_with_distance_matrix([vec![0.999]].as_slice().into())
            .unwrap(),
        0.0
    );
    assert_eq!(
        carbon
            .ipc_with_distance_matrix((&other).into())
            .unwrap_err(),
        IpcError::InvalidDistanceMatrixShape {
            row: None,
            actual: 2,
            expected: 1
        }
    );
    assert_eq!(
        mol.avg_ipc_with_distance_matrix([vec![0.0, 0.0], vec![1.0]].as_slice().into())
            .unwrap_err(),
        IpcError::InvalidDistanceMatrixShape {
            row: Some(1),
            actual: 1,
            expected: 2
        }
    );
}
