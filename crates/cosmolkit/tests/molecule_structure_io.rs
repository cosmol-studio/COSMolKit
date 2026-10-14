#![cfg(feature = "cap-io")]
//! Ordinary behavior regressions for molecular structural text IO.
use cosmolkit::*;
use std::error::Error;

#[cfg(feature = "cap-bio")]
const PDB: &str = "HETATM    1  C1  LIG A   7       1.000   2.000   3.000  0.50 12.00           C  \nHETATM    2  H1  LIG A   7       1.800   2.000   3.000  1.00  4.00           H  \nCONECT    1    2\nEND\n";
#[cfg(feature = "cap-bio")]
const MMCIF: &str = "data_lig\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\nHETATM 1 C C1 . LIG A A 7 1.000 2.000 3.000 0.50 12.00\n";

fn carbon(coordinates: CoordinateBlock) -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
        vec![],
        vec![],
        vec![],
    )
    .unwrap();
    Molecule::from_parts(topology, coordinates, MoleculeProperties::default()).unwrap()
}
fn xyz(text: &str) -> Vec<[f64; 3]> {
    text.lines()
        .filter(|l| l.starts_with("ATOM  ") || l.starts_with("HETATM"))
        .map(|l| {
            [
                l[30..38].trim().parse().unwrap(),
                l[38..46].trim().parse().unwrap(),
                l[46..54].trim().parse().unwrap(),
            ]
        })
        .collect()
}

#[test]
fn pdb_auto_zero_or_unique_and_explicit_dimension_ids() {
    assert_eq!(
        xyz(&carbon(CoordinateBlock::default()).to_pdb().unwrap()),
        vec![[0.0; 3]]
    );
    let two = Conformer2D::new(7, vec![[1.25, -2.5]]);
    let three = Conformer3D::new(7, vec![[4.0, 5.0, 6.0]], true);
    for coordinates in [
        CoordinateBlock {
            conformers_2d: vec![two.clone()],
            ..Default::default()
        },
        CoordinateBlock {
            conformers_3d: vec![three.clone()],
            ..Default::default()
        },
    ] {
        let molecule = carbon(coordinates);
        assert_eq!(xyz(&molecule.to_pdb().unwrap()).len(), 1);
        let all = molecule
            .to_pdb_with_params(&PdbWriteParams {
                all_conformers: true,
                ..Default::default()
            })
            .unwrap();
        assert!(all.starts_with("MODEL        1\n"));
        assert!(all.ends_with("ENDMDL\nEND\n"));
    }
    let mol = carbon(CoordinateBlock {
        conformers_2d: vec![two],
        conformers_3d: vec![three],
        ..Default::default()
    });
    let before = mol.clone();
    assert!(matches!(
        mol.to_pdb(),
        Err(MolecularIoError::PdbWrite(
            PdbWriteError::AmbiguousCoordinates {
                two_d: 1,
                three_d: 1
            }
        ))
    ));
    for (selection, expected) in [
        (MolCoordinateSelection::TwoD { id: 7 }, [1.25, -2.5, 0.0]),
        (MolCoordinateSelection::ThreeD { id: 7 }, [4.0, 5.0, 6.0]),
    ] {
        assert_eq!(
            xyz(&mol
                .to_pdb_with_params(&PdbWriteParams {
                    coordinate_selection: selection,
                    ..Default::default()
                })
                .unwrap()),
            vec![expected]
        );
    }
    let missing = mol
        .to_pdb_with_params(&PdbWriteParams {
            coordinate_selection: MolCoordinateSelection::ThreeD { id: 0 },
            ..Default::default()
        })
        .unwrap_err();
    assert!(matches!(
        missing,
        MolecularIoError::PdbWrite(PdbWriteError::MissingCoordinate {
            dimension: CoordinateDimension::ThreeD,
            id: 0
        })
    ));
    assert!(missing.source().unwrap().is::<PdbWriteError>());
    assert_eq!(mol, before);
}

#[test]
fn pdb_all_enumerates_noncontiguous_ids_without_provenance_or_mutation() {
    let mol = carbon(CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(90, vec![[1.0, 2.0]])],
        conformers_3d: vec![
            Conformer3D::new(7, vec![[3.0, 4.0, 5.0]], true),
            Conformer3D::new(42, vec![[6.0, 7.0, 8.0]], true),
        ],
        ..Default::default()
    });
    let before = mol.clone();
    let text = mol
        .to_pdb_with_params(&PdbWriteParams {
            all_conformers: true,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        xyz(&text),
        vec![[1.0, 2.0, 0.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]
    );
    assert_eq!(
        text.lines()
            .filter(|l| l.starts_with("MODEL"))
            .collect::<Vec<_>>(),
        ["MODEL        1", "MODEL        2", "MODEL        3"]
    );
    assert_eq!(text.matches("ENDMDL\n").count(), 3);
    assert!(text.ends_with("END\n"));
    assert!(matches!(
        mol.to_pdb_with_params(&PdbWriteParams {
            all_conformers: true,
            coordinate_selection: MolCoordinateSelection::TwoD { id: 90 },
            ..Default::default()
        }),
        Err(MolecularIoError::PdbWrite(
            PdbWriteError::ConflictingCoordinateSelection
        ))
    ));
    assert_eq!(mol, before);
}

#[cfg(feature = "cap-smiles")]
#[test]
fn pdb_aromatic_preprocessing_is_local_and_writes_kekule_conect() {
    let mol = Molecule::from_smiles("c1ccccc1").unwrap();
    let before = mol.clone();
    let text = mol.to_pdb().unwrap();
    // PDBWriter.cpp Kekulizes its copy before writing bond multiplicity.
    // Fixed expected CONECT rows checked against pinned RDKit 2026.09.1;
    // canonical Kekulization places the first double bond on closure 1-6.
    assert_eq!(
        text.lines()
            .filter(|l| l.starts_with("CONECT"))
            .collect::<Vec<_>>(),
        [
            "CONECT    1    2    6    6",
            "CONECT    2    3    3",
            "CONECT    3    4",
            "CONECT    4    5    5",
            "CONECT    5    6",
        ]
    );
    assert_eq!(mol, before);
    assert!(mol.topology().bonds.iter().all(|b| b.is_aromatic()));
}

#[cfg(feature = "cap-bio")]
#[test]
fn molecular_readers_reuse_bio_options_and_keep_residue_metadata() {
    let molecule = Molecule::from_pdb(PDB).unwrap();
    assert_eq!(molecule.num_atoms(), 1);
    assert_eq!(molecule.num_bonds(), 0);
    assert_eq!(
        molecule.conformers_3d()[0].coordinates(),
        &[[1.0, 2.0, 3.0]]
    );
    let info = molecule.topology().atoms[0].pdb_residue_info().unwrap();
    assert_eq!(info.residue_number(), 7);
    assert_eq!(info.chain_id().as_bytes(), b"A");
    assert_eq!(info.occupancy(), 0.5);
    assert_eq!(mol_from_pdb(PDB).unwrap(), molecule);
    let params = BioMoleculeParams {
        remove_hs: false,
        proximity_bonding: false,
        ..Default::default()
    };
    let retained = Molecule::from_pdb_with_params(PDB, &params).unwrap();
    assert_eq!(retained.num_atoms(), 2);
    assert_eq!(retained.num_bonds(), 1);
    let proximity_only = PDB.replace("CONECT    1    2\n", "");
    assert_eq!(
        Molecule::from_pdb_with_params(&proximity_only, &params)
            .unwrap()
            .num_bonds(),
        0
    );
    assert_eq!(
        Molecule::from_pdb_with_params(
            &proximity_only,
            &BioMoleculeParams {
                proximity_bonding: true,
                ..params
            }
        )
        .unwrap()
        .num_bonds(),
        1
    );
    let cif = Molecule::from_mmcif(MMCIF).unwrap();
    assert_eq!(cif.num_atoms(), 1);
    assert_eq!(cif.conformers_3d()[0].coordinates(), &[[1.0, 2.0, 3.0]]);
    assert_eq!(
        cif.topology().atoms[0]
            .pdb_residue_info()
            .unwrap()
            .residue_number(),
        7
    );
    assert_eq!(mol_from_mmcif(MMCIF).unwrap(), cif);
    assert_eq!(
        Molecule::from_mmcif_with_params(MMCIF, &params)
            .unwrap()
            .num_atoms(),
        1
    );
}

#[cfg(feature = "cap-bio")]
#[test]
fn structural_reader_failure_retains_typed_cause() {
    let error = Molecule::from_mmcif("not a CIF block").unwrap_err();
    assert!(matches!(error, MolecularIoError::MmcifRead(_)));
    assert!(error.source().unwrap().is::<BioMmcifReadError>());
}
