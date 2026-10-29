//! Native canonical BIO writer re-exports and destination preservation.
use crate::{BioMmcifWriteParams, BioPdbWriteError, BioPdbWriteParams, BioStructure};
#[test]
fn bio_writers_reexports_text_file_and_prewrite_failure_atomicity() {
    let pdb =
        "ATOM      7  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nEND\n";
    let cif = "data_hierarchy\n_entity.id e1\n_entity.type polymer\n_entity_poly.entity_id e1\n_entity_poly.type \"polypeptide(L)\"\n_struct_asym.id X\n_struct_asym.entity_id e1\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\ne1 1 ALA\ne1 2 GLY\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_alt_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\nATOM 19 C CA ALA . X e1 1 A -4 1 2 3 1 20\n";

    let s = BioStructure::from_pdb(pdb).unwrap();
    let p = BioPdbWriteParams::default();
    let c = BioMmcifWriteParams::default();
    assert_eq!(s.to_pdb().unwrap(), s.to_pdb_with_params(&p).unwrap());
    assert_eq!(s.to_mmcif().unwrap(), s.to_mmcif_with_params(&c).unwrap());
    let dir = std::env::temp_dir().join(format!(
        "cosmolkit-bio-writer-binding-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("coordinates.pdb");
    s.write_pdb(&path).unwrap();
    assert_eq!(std::fs::read_to_string(&path).unwrap(), s.to_pdb().unwrap());
    s.write_pdb_with_params(&path, &p).unwrap();
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        s.to_pdb_with_params(&p).unwrap()
    );
    s.write_mmcif(&path).unwrap();
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        s.to_mmcif().unwrap()
    );
    s.write_mmcif_with_params(&path, &c).unwrap();
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        s.to_mmcif_with_params(&c).unwrap()
    );
    let bad = BioStructure::from_mmcif(&cif.replace("ATOM 19", "ATOM -1")).unwrap();
    let preserve = BioPdbWriteParams {
        preserve_serial: true,
        ..p
    };
    let before = std::fs::read(&path).unwrap();
    assert!(matches!(
        bad.write_pdb_with_params(&path, &preserve),
        Err(BioPdbWriteError::NegativeSerial { serial: -1 })
    ));
    assert_eq!(std::fs::read(&path).unwrap(), before);
    assert_eq!(bad.num_atoms(), 1);
}
