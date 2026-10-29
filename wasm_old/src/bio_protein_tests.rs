//! Native Protein re-export behavior and public iterator contracts.
use crate::{
    BioCoordinateFormat, BioPdbReadParams, BioReadParams, BioSelection, BioStructure, Protein,
};
#[test]
fn bio_protein_reexports_projection_iterators_cow_and_file_read() {
    let pdb = "ATOM      1  CA AALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nHETATM    2  SE  MSE B   2       4.000   5.000   6.000  1.00 20.00          SE  \nHETATM    3  O   HOH C   3       7.000   8.000   9.000  1.00 20.00           O  \nEND\n";
    let cif = "data_hierarchy\n_entity.id e1\n_entity.type polymer\n_entity_poly.entity_id e1\n_entity_poly.type \"polypeptide(L)\"\n_struct_asym.id X\n_struct_asym.entity_id e1\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\ne1 1 ALA\ne1 2 GLY\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_alt_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\nATOM 19 C CA ALA . X e1 1 A -4 1 2 3 1 20\n";

    let s = BioStructure::from_pdb(pdb).unwrap();
    let p = s.protein().unwrap();
    assert_eq!(s.num_atoms(), 3);
    assert_eq!(p.num_atoms(), 2);
    assert_eq!(p.selection_summary().atoms, 2);
    assert_eq!(p.chains().len(), 2);
    assert_eq!(p.residues().len(), 2);
    assert_eq!(p.atoms().len(), 2);
    let mut it = p.atoms();
    assert_eq!(it.next().unwrap().id().value(), 0);
    assert_eq!(it.len(), 1);
    assert_eq!(it.next().unwrap().id().value(), 1);
    assert!(it.next().is_none());
    assert!(it.next().is_none());
    let snapshot = p.clone();
    let mut changed = p.with_translated_coordinates([1., 0., 0.]).unwrap();
    assert_eq!(changed.atoms().next().unwrap().position(), [2., 2., 3.]);
    assert_eq!(snapshot.atoms().next().unwrap().position(), [1., 2., 3.]);
    let select = BioSelection::from_cid("/1/A/1").unwrap();
    changed.retain_selection_(&select).unwrap();
    assert_eq!(changed.num_atoms(), 1);
    assert_eq!(snapshot.num_atoms(), 2);
    assert_eq!(p.clone().into_bio_structure().num_atoms(), 2);
    assert_eq!(Protein::from_pdb(pdb).unwrap(), p);
    assert_eq!(
        Protein::from_pdb_with_params(pdb, &BioPdbReadParams::default()).unwrap(),
        p
    );
    assert_eq!(Protein::from_text(pdb).unwrap().num_atoms(), 2);
    assert_eq!(
        Protein::from_text_with_params(
            pdb,
            &BioReadParams {
                format: BioCoordinateFormat::Pdb,
                source_name: "fixture".into()
            }
        )
        .unwrap()
        .num_atoms(),
        2
    );
    assert_eq!(Protein::from_mmcif(cif).unwrap().num_atoms(), 1);
    let dir = std::env::temp_dir().join(format!(
        "cosmolkit-protein-binding-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("protein.pdb");
    std::fs::write(&path, pdb).unwrap();
    assert_eq!(Protein::read(&path).unwrap().num_atoms(), 2);
    assert_eq!(
        Protein::read_with_format(&path, BioCoordinateFormat::Pdb)
            .unwrap()
            .num_atoms(),
        2
    );
}
