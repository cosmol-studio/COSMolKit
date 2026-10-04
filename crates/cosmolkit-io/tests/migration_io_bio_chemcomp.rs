use cosmolkit_bio::{BioCoordinateFormat, BioStructureData};
use cosmolkit_io::read_bio_structure_file;

const COORDS: [(&str, &str, &str); 3] = [
    ("x", "y", "z"),
    ("model_Cartn_x", "model_Cartn_y", "model_Cartn_z"),
    (
        "pdbx_model_Cartn_x_ideal",
        "pdbx_model_Cartn_y_ideal",
        "pdbx_model_Cartn_z_ideal",
    ),
];

fn document(shape: usize, mask: u8) -> String {
    // mmcif.hpp selects exactly the CCD, two-block or global/comp_list shape.
    let mut text = match shape {
        0 => "data_comp_TST\n".to_owned(),
        1 => "data_comp_list\n_chem_comp.id LIST\ndata_comp_TST\n".to_owned(),
        2 => "global_\ndata_comp_list\n_chem_comp.id LIST\ndata_comp_TST\n".to_owned(),
        _ => unreachable!(),
    };
    text.push_str("_chem_comp.id 'TST'\n");
    text.push_str("loop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.comp_id\n_chem_comp_atom.charge\n");
    for (bit, tags) in COORDS.iter().enumerate() {
        if mask & (1 << bit) != 0 {
            for tag in [tags.0, tags.1, tags.2] {
                text.push_str(&format!("_chem_comp_atom.{tag}\n"));
            }
        }
    }
    for (index, atom) in ["C2", "N1"].iter().enumerate() {
        text.push_str(&format!(
            "{atom} {} TST {}",
            ["C", "N"][index],
            ["1.6", "."][index]
        ));
        for bit in 0..3 {
            if mask & (1 << bit) != 0 {
                let base = bit * 10 + index * 3;
                text.push_str(&format!(" {} {} {}", base + 1, base + 2, base + 3));
            }
        }
        text.push('\n');
    }
    text
}

fn check(data: &BioStructureData, mask: u8) {
    data.validate().unwrap();
    assert_eq!(data.input_format(), BioCoordinateFormat::ChemComp);
    assert_eq!(data.source_state().name, "'TST'");
    let count = mask.count_ones() as usize;
    assert_eq!(
        (
            data.models().len(),
            data.chains().len(),
            data.residues().len(),
            data.atoms().len()
        ),
        (count, count, count, count * 2)
    );
    for (index, bit) in (0..3).filter(|bit| mask & (1 << bit) != 0).enumerate() {
        let model = data.models()[index];
        assert_eq!(model.source_model_number(), Some(index as i32 + 1));
        assert_eq!(
            (model.chain_span().start(), model.chain_span().len()),
            (index as u32, 1)
        );
        let chain = &data.chains()[index];
        assert_eq!(
            (chain.residue_span().start(), chain.residue_span().len()),
            (index as u32, 1)
        );
        let residue = &data.residues()[index];
        assert_eq!(
            (residue.atom_span().start(), residue.atom_span().len()),
            ((index * 2) as u32, 2)
        );
        assert_eq!(residue.name().as_str(), "TST");
        assert_eq!(residue.source().seq_id().unwrap().seq_num(), 1);
        for row in 0..2 {
            let atom = &data.atoms()[index * 2 + row];
            assert_eq!(atom.name().as_str(), ["C2", "N1"][row]);
            assert_eq!(atom.element().symbol(), ["C", "N"][row]);
            assert_eq!(atom.formal_charge(), [2, 0][row]);
            assert_eq!(atom.source().serial().unwrap().value(), 0);
            assert_eq!(atom.occupancy(), 1.0);
            assert_eq!(atom.b_iso(), 20.0);
            let base = bit * 10 + row * 3;
            assert_eq!(
                data.coordinates().positions()[index * 2 + row],
                [(base + 1) as f64, (base + 2) as f64, (base + 3) as f64]
            );
        }
    }
}

#[test]
fn chemcomp_detached_three_recognizers_and_all_coordinate_presences() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("component.cif");
    for shape in 0..3 {
        for mask in 0..8 {
            std::fs::write(&path, document(shape, mask)).unwrap();
            let data = read_bio_structure_file(&path, BioCoordinateFormat::ChemComp).unwrap();
            check(&data, mask);
        }
    }
}

#[test]
fn chemcomp_detached_pair_columns_and_null_only_coordinate_are_not_present() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("pair.cif");
    std::fs::write(&path, "data_comp_TST\n_chem_comp.id 'TST'\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.comp_id TST\n_chem_comp_atom.x .\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n_chem_comp_atom.model_Cartn_x 4.5\n_chem_comp_atom.model_Cartn_y 6\n_chem_comp_atom.model_Cartn_z 7\n").unwrap();
    let data = read_bio_structure_file(&path, BioCoordinateFormat::ChemComp).unwrap();
    assert_eq!(data.source_state().name, "'TST'");
    assert_eq!(data.models().len(), 1);
    assert_eq!(
        (
            data.residues()[0].atom_span().start(),
            data.residues()[0].atom_span().len()
        ),
        (0, 1)
    );
    assert_eq!(data.atoms()[0].name().as_str(), "C1");
    assert_eq!(data.coordinates().positions(), &[[4.5, 6.0, 7.0]]);
}
