#![cfg(all(feature = "cap-io", feature = "cap-batch"))]
use cosmolkit::*;
use std::io::Cursor;

#[test]
fn all_public_sdf_paths_receive_one_finalization_before_lists_and_metadata() {
    let text = concat!(
        "source-order\n  COSMolKit         2D\n\n",
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n",
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "  1  2  1  0  0  0  0\nM  END\n",
        ">  <atom.iprop.rank>\n7\n\n",
        ">  <_StereochemDone>\nafter-finalization\n\n$$$$\n"
    );
    let assert_molecule = |molecule: &Molecule| {
        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(
            molecule.atom(AtomId::new(0)).unwrap().prop("rank"),
            Some(&PropertyValue::Int(7))
        );
        assert_eq!(
            molecule.property("_StereochemDone"),
            Some(&PropertyValue::from("after-finalization"))
        );
        assert_eq!(
            molecule.properties().sdf_data_fields(),
            [
                (
                    PropertyText::from("atom.iprop.rank"),
                    PropertyText::from("7")
                ),
                (
                    PropertyText::from("_StereochemDone"),
                    PropertyText::from("after-finalization")
                )
            ]
        );
    };
    let direct = SdfRecord::from_sdf(text).unwrap();
    assert_molecule(direct.molecule().unwrap());
    assert_molecule(&Molecule::from_sdf(text).unwrap());
    // MolFromMolBlock finalizes the same CTAB but leaves trailing SD fields unread.
    let mol = Molecule::from_mol(text).unwrap();
    assert_eq!(mol.num_atoms(), 1);
    assert_eq!(mol.atom(AtomId::new(0)).unwrap().prop("rank"), None);
    assert!(mol.properties().sdf_data_fields().is_empty());
    let unremoved = Molecule::from_mol_with_params(
        text,
        &SdfReadParams {
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(unremoved.num_atoms(), 2);
    let mut stream = SdfRecordStream::new(Cursor::new(text.as_bytes()));
    assert_molecule(stream.next_record().unwrap().unwrap().molecule().unwrap());
    assert!(stream.next_record().unwrap().is_none());
    let path = std::env::temp_dir().join(format!(
        "cosmolkit-sdf-source-order-{}.sdf",
        std::process::id()
    ));
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)
        .unwrap();
    std::io::Write::write_all(&mut file, text.as_bytes()).unwrap();
    drop(file);
    let dataset = SdfDataset::open(path.to_str().unwrap()).unwrap();
    assert_eq!(dataset.len(), 1);
    assert_molecule(dataset.record(0).unwrap().molecule().unwrap());
    let selected =
        MoleculeBatch::from_dataset_indices(&dataset, &[0, 0], BatchErrorMode::Strict).unwrap();
    assert_eq!(selected.len(), 2);
    for record in selected.records() {
        let BatchRecord::Molecule(molecule) = record else {
            panic!("retained supported row");
        };
        assert_molecule(molecule);
    }
    let batch = SdfRecordStream::new(Cursor::new(text.as_bytes()))
        .batches(1, BatchErrorMode::Strict, Some(1))
        .unwrap()
        .next_batch()
        .unwrap()
        .unwrap();
    let BatchRecord::Molecule(molecule) = batch.get(0).unwrap() else {
        panic!("retained supported row");
    };
    assert_molecule(molecule);
    std::fs::remove_file(path).unwrap();
}
