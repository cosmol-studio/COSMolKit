use cosmolkit_bio::{BioAsu, BioConnectionKind, BioCoordinateFormat, EntityKind};
use cosmolkit_io::{BioPdbReadParams, read_pdb_bio_structure};

fn fixed_record(tag: &[u8; 6], payload: &[u8]) -> String {
    let mut line = [b' '; 80];
    line[..6].copy_from_slice(tag);
    let end = 10 + payload.len();
    assert!(end <= line.len());
    line[10..end].copy_from_slice(payload);
    let mut text = String::from_utf8(line.to_vec()).unwrap();
    text.push('\n');
    text
}

fn pdb_atom(
    tag: &[u8; 6],
    serial: i32,
    atom_name: &[u8; 4],
    altloc: u8,
    residue: &[u8; 3],
    sequence: i32,
    position: [f64; 3],
    element: &[u8; 2],
) -> String {
    let mut line = [b' '; 80];
    line[..6].copy_from_slice(tag);
    line[6..11].copy_from_slice(format!("{serial:>5}").as_bytes());
    line[12..16].copy_from_slice(atom_name);
    line[16] = altloc;
    line[17..20].copy_from_slice(residue);
    line[20..22].copy_from_slice(b"A ");
    line[22..26].copy_from_slice(format!("{sequence:>4}").as_bytes());
    line[30..38].copy_from_slice(format!("{:>8.3}", position[0]).as_bytes());
    line[38..46].copy_from_slice(format!("{:>8.3}", position[1]).as_bytes());
    line[46..54].copy_from_slice(format!("{:>8.3}", position[2]).as_bytes());
    line[54..60].copy_from_slice(b"  1.00");
    line[60..66].copy_from_slice(b" 20.00");
    line[76..78].copy_from_slice(element);
    let mut text = String::from_utf8(line.to_vec()).unwrap();
    text.push('\n');
    text
}

fn pdb_cryst1_p21() -> String {
    let mut line = [b' '; 80];
    line[..6].copy_from_slice(b"CRYST1");
    line[6..15].copy_from_slice(b"   10.000");
    line[15..24].copy_from_slice(b"   10.000");
    line[24..33].copy_from_slice(b"   10.000");
    line[33..40].copy_from_slice(b"  90.00");
    line[40..47].copy_from_slice(b"  90.00");
    line[47..54].copy_from_slice(b"  90.00");
    line[55..66].copy_from_slice(b"P 1 21 1   ");
    let mut text = String::from_utf8(line.to_vec()).unwrap();
    text.push('\n');
    text
}

fn pdb_ssbond_without_symop_columns() -> String {
    let mut line = [b' '; 60];
    line[..6].copy_from_slice(b"SSBOND");
    line[7..10].copy_from_slice(b"  1");
    line[11..14].copy_from_slice(b"CYS");
    line[14..16].copy_from_slice(b"A ");
    line[17..22].copy_from_slice(b"   1 ");
    line[25..28].copy_from_slice(b"CYS");
    line[28..30].copy_from_slice(b"A ");
    line[31..36].copy_from_slice(b"   2 ");
    let mut text = String::from_utf8(line.to_vec()).unwrap();
    text.push('\n');
    text
}

fn parse(text: &str) -> cosmolkit_bio::BioStructureData {
    read_pdb_bio_structure(text, "fixture.pdb", &BioPdbReadParams::default()).unwrap()
}

#[test]
fn p42_finalizes_default_model_and_ter_entity_rollback() {
    let empty = parse("");
    assert_eq!(empty.input_format(), BioCoordinateFormat::Pdb);
    assert_eq!(empty.models().len(), 1);
    assert!(empty.atoms().is_empty());
    empty.validate().unwrap();

    let mut water = pdb_atom(
        b"HETATM",
        1,
        b" O  ",
        b' ',
        b"HOH",
        1,
        [0.0, 0.0, 0.0],
        b" O",
    );
    water.push_str("TER\n");
    let rolled_back = parse(&water);
    assert_eq!(rolled_back.source_state().ter_status, b'e');
    assert_eq!(rolled_back.residues().len(), 1);
    assert_eq!(rolled_back.residues()[0].entity_kind(), EntityKind::Unknown);
    rolled_back.validate().unwrap();
}

#[test]
fn p42_prepares_cell_images_before_ssbond_altloc_completion() {
    let mut input = pdb_cryst1_p21();
    input.push_str(&pdb_ssbond_without_symop_columns());
    // For the pinned P 1 21 1 screw operation, pair A/A is zero-distance
    // through a crystallographic image; without that image B/B is closer.
    // Gemmi calls setup_cell_images before process_conn/complete_ssbond.
    input.push_str(&pdb_atom(
        b"ATOM  ",
        1,
        b" SG ",
        b'A',
        b"CYS",
        1,
        [1.0, 0.0, 1.0],
        b" S",
    ));
    input.push_str(&pdb_atom(
        b"ATOM  ",
        2,
        b" SG ",
        b'B',
        b"CYS",
        1,
        [0.0, 0.0, 0.0],
        b" S",
    ));
    input.push_str(&pdb_atom(
        b"ATOM  ",
        3,
        b" SG ",
        b'A',
        b"CYS",
        2,
        [-1.0, -5.0, -1.0],
        b" S",
    ));
    input.push_str(&pdb_atom(
        b"ATOM  ",
        4,
        b" SG ",
        b'B',
        b"CYS",
        2,
        [0.1, 0.0, 0.0],
        b" S",
    ));

    let structure = parse(&input);
    let crystal = structure.crystal().unwrap();
    // Gemmi's P 1 21 1 group has two operations including identity, while
    // set_cell_images_from_groupops stores only the non-identity operation.
    assert_eq!(crystal.cs_count(), 1);
    assert_eq!(crystal.symmetry_images().len(), 1);
    assert_eq!(structure.connections().len(), 1);
    let connection = &structure.connections()[0];
    assert_eq!(connection.kind, BioConnectionKind::Disulf);
    assert_eq!(connection.asu, BioAsu::Any);
    assert_eq!(connection.partner1.altloc(), b'A');
    assert_eq!(connection.partner2.altloc(), b'A');
    structure.validate().unwrap();
}

#[test]
fn p42_finalizes_author_continuations_and_optional_remarks() {
    let mut input = String::new();
    input.push_str(&fixed_record(b"AUTHOR", b"Doe, A.-"));
    input.push_str(&fixed_record(b"AUTHOR", b"B.Smith, C.Jones"));
    input.push_str("REMARK   2 RESOLUTION.    2.50 ANGSTROMS.\n");

    let interpreted = parse(&input);
    assert_eq!(
        interpreted.metadata().authors,
        ["Doe", "Smith, A.-B.", "Jones, C."]
    );
    assert_eq!(interpreted.source_state().raw_remarks.len(), 1);
    assert_eq!(interpreted.source_state().resolution, 2.5);
    interpreted.validate().unwrap();

    let skipped = read_pdb_bio_structure(
        &input,
        "fixture.pdb",
        &BioPdbReadParams {
            skip_remarks: true,
            ..BioPdbReadParams::default()
        },
    )
    .unwrap();
    assert_eq!(skipped.metadata().authors, interpreted.metadata().authors);
    assert_eq!(skipped.source_state().raw_remarks.len(), 1);
    assert_eq!(skipped.source_state().resolution, 0.0);
    skipped.validate().unwrap();
}

#[test]
fn p42_restores_short_hetname_codes_after_complete_dispatch() {
    let mut alias = [b' '; 80];
    alias[..6].copy_from_slice(b"HETNAM");
    alias[11..14].copy_from_slice(b"AAA");
    alias[71..74].copy_from_slice(b"BBB");
    let mut input = String::from_utf8(alias.to_vec()).unwrap();
    input.push('\n');
    input.push_str(&pdb_atom(
        b"HETATM",
        1,
        b" CA ",
        b' ',
        b"AAA",
        1,
        [0.0, 0.0, 0.0],
        b" C",
    ));

    let structure = parse(&input);
    assert_eq!(structure.residues().len(), 1);
    assert_eq!(structure.residues()[0].name().as_str(), "BBB");
    structure.validate().unwrap();
}
