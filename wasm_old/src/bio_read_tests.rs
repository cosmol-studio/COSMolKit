//! Native checks of the binding crate's canonical BIO re-export.
use crate::{BioCoordinateFormat, BioPdbReadParams, BioReadParams, BioStructure};
const PDB: &str =
    "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nEND\n";
const INCOMPLETE_CIF: &str = "data_test\n_entry.id test\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\nATOM 1 C CA ALA A 1 1.000 2.000 3.000 1.00 20.00\n";
const CIF: &str = "data_test\n_entry.id test\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_alt_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\nATOM 1 C CA ALA . A A 1 1.000 2.000 3.000 1.00 20.00\n";
#[test]
fn all_seven_structural_reader_routes_preserve_counts_formats_and_file_values() {
    assert_eq!(
        BioStructure::from_mmcif(INCOMPLETE_CIF)
            .unwrap()
            .num_atoms(),
        0
    );
    let values = [
        BioStructure::from_pdb(PDB).unwrap(),
        BioStructure::from_pdb_with_params(PDB, &BioPdbReadParams::default()).unwrap(),
        BioStructure::from_mmcif(CIF).unwrap(),
        BioStructure::from_text(PDB).unwrap(),
        BioStructure::from_text_with_params(
            CIF,
            &BioReadParams {
                format: BioCoordinateFormat::Mmcif,
                source_name: "named.cif".into(),
            },
        )
        .unwrap(),
    ];
    for v in &values {
        assert_eq!(
            (
                v.num_models(),
                v.num_chains(),
                v.num_residues(),
                v.num_atoms()
            ),
            (1, 1, 1, 1)
        );
    }
    assert_eq!(values[0], values[1]);
    assert_eq!(values[0].input_format(), BioCoordinateFormat::Pdb);
    assert_eq!(values[2].input_format(), BioCoordinateFormat::Mmcif);
    let executable = std::env::current_exe().unwrap();
    let directory = tempfile::Builder::new()
        .prefix("bio-reader-fixtures-")
        .tempdir_in(executable.parent().unwrap())
        .unwrap();
    let path = directory.path().join("fixture.pdb");
    std::fs::write(&path, PDB).unwrap();
    for v in [
        BioStructure::read(&path).unwrap(),
        BioStructure::read_with_format(&path, BioCoordinateFormat::Pdb).unwrap(),
    ] {
        assert_eq!(v.num_atoms(), 1);
        assert_eq!(
            v.atom_position(crate::BioAtomId::new(0)),
            Some([1., 2., 3.])
        );
    }
    let error = BioStructure::from_pdb("ATOM  \n").unwrap_err();
    assert_eq!(error.stage(), crate::BioPdbReadStage::Record);
    assert_eq!(error.line_number(), Some(1));
    assert_eq!(error.record_tag(), Some(*b"ATOM"));
    assert!(
        BioStructure::from_pdb("{\"not\":\"pdb\"}")
            .unwrap()
            .num_atoms()
            == 0
    );
}
