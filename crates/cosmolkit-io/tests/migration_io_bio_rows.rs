//! BIO-ROWS R12: end-to-end privately constructed selector cases against
//! detached PDB/mmCIF reader results, through the owning IO crate.
//!
//! Expected row IDs and logical names are derived independently from the
//! pinned reader semantics (file order defines the row tables; PDB names
//! trim through read_string; mmCIF names stay verbatim; auth chain
//! identity is preferred over label; D projects to H+Some(2)).

use cosmolkit_bio::__bio_rows_probe::{ProbeHit, ProbeSelection};
use cosmolkit_bio::BioCoordinateFormat;
use cosmolkit_io::{BioReadParams, read_bio_structure};

/// Column-exact PDB ATOM line builder (1-based columns per the format).
fn pdb_atom(
    serial: u32,
    name: &str,
    altloc: char,
    resname: &str,
    chain: char,
    resseq: i32,
    occ: f64,
    b: f64,
    element: &str,
) -> String {
    format!(
        "ATOM  {serial:>5} {name:<4}{altloc}{resname:>3} {chain}{resseq:>4}    {x:>8.3}{y:>8.3}{z:>8.3}{occ:>6.2}{b:>6.2}          {element:>2}",
        serial = serial,
        name = name,
        altloc = altloc,
        resname = resname,
        chain = chain,
        resseq = resseq,
        x = 0.0,
        y = 0.0,
        z = 0.0,
        occ = occ,
        b = b,
        element = element,
    ) + "\n"
}

fn pdb_fixture() -> String {
    // Single implicit model (number 1 via size+1), chain A with ALA 1 and
    // GLY 2, chain B with HOH 1; the deuterium atom carries element D.
    [
        pdb_atom(1, "N", ' ', "ALA", 'A', 1, 1.00, 20.0, "N"),
        pdb_atom(2, "CA", ' ', "ALA", 'A', 1, 1.00, 20.0, "C"),
        pdb_atom(3, "D", ' ', "ALA", 'A', 1, 1.00, 20.0, "D"),
        pdb_atom(4, "N", 'A', "GLY", 'A', 2, 0.50, 20.0, "N"),
        pdb_atom(5, "O", ' ', "HOH", 'B', 1, 0.50, 20.0, "O"),
        "TER\nEND\n".to_string(),
    ]
    .concat()
}

fn mmcif_fixture() -> String {
    // Same logical structure: label asym ids L1/L2 with AUTH asym ids A/B
    // (differing on purpose), one deuterium (type_symbol D), altLoc A,
    // occupancy 0.5 rows, single implicit model.
    [
        "data_r12\n".to_string(),
        "#\nloop_\n".to_string(),
        "_atom_site.group_PDB\n".to_string(),
        "_atom_site.id\n".to_string(),
        "_atom_site.type_symbol\n".to_string(),
        "_atom_site.label_atom_id\n".to_string(),
        "_atom_site.label_alt_id\n".to_string(),
        "_atom_site.label_comp_id\n".to_string(),
        "_atom_site.label_asym_id\n".to_string(),
        "_atom_site.label_entity_id\n".to_string(),
        "_atom_site.label_seq_id\n".to_string(),
        "_atom_site.Cartn_x\n".to_string(),
        "_atom_site.Cartn_y\n".to_string(),
        "_atom_site.Cartn_z\n".to_string(),
        "_atom_site.occupancy\n".to_string(),
        "_atom_site.B_iso_or_equiv\n".to_string(),
        "_atom_site.auth_seq_id\n".to_string(),
        "_atom_site.auth_comp_id\n".to_string(),
        "_atom_site.auth_asym_id\n".to_string(),
        "ATOM 1 N N . ALA L1 1 1 0.0 0.0 0.0 1.00 20.0 1 ALA A\n".to_string(),
        "ATOM 2 C CA . ALA L1 1 1 0.0 0.0 0.0 1.00 20.0 1 ALA A\n".to_string(),
        "ATOM 3 D D . ALA L1 1 1 0.0 0.0 0.0 1.00 20.0 1 ALA A\n".to_string(),
        "ATOM 4 N N A GLY L1 1 2 0.0 0.0 0.0 0.50 20.0 2 GLY A\n".to_string(),
        "ATOM 5 O O . HOH L2 2 1 0.0 0.0 0.0 0.50 20.0 1 HOH B\n".to_string(),
        "#\n".to_string(),
    ]
    .concat()
}

fn read(text: &str, format: BioCoordinateFormat) -> cosmolkit_bio::BioStructureData {
    read_bio_structure(
        text,
        &BioReadParams {
            format,
            source_name: "r12".to_string(),
        },
    )
    .expect("reader fixture must parse")
}

fn wildcard() -> ProbeSelection {
    ProbeSelection::default()
}

#[test]
fn bio_rows_r12_same_logical_structure_across_formats() {
    // Independently derived table layout (file order): one model, chain A
    // (residues ALA, GLY), chain B (HOH); atoms N, CA, D, N(A), O.
    let expected_first = ProbeHit {
        model: 0,
        chain: 0,
        residue: 0,
        atom: 0,
    };
    for (data, format_name) in [
        (read(&pdb_fixture(), BioCoordinateFormat::Pdb), "pdb"),
        (read(&mmcif_fixture(), BioCoordinateFormat::Mmcif), "mmcif"),
    ] {
        assert_eq!(
            cosmolkit_bio::__bio_rows_probe::probe_first_rows(&wildcard(), &data, b' ', b' '),
            Ok(Some(expected_first)),
            "{format_name}"
        );
        assert_eq!(
            cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(
                &wildcard(),
                &data,
                b' ',
                b' '
            ),
            Ok(vec![0u32, 1, 2, 3, 4]),
            "{format_name}"
        );
    }
}

#[test]
fn bio_rows_r12_auth_label_chain_identity() {
    let pdb = read(&pdb_fixture(), BioCoordinateFormat::Pdb);
    let cif = read(&mmcif_fixture(), BioCoordinateFormat::Mmcif);

    // Auth identity matches in both formats.
    let mut auth_a = wildcard();
    auth_a.chain_list = Some("A".to_string());
    assert_eq!(
        cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&auth_a, &pdb, b' ', b' '),
        Ok(vec![0u32, 1, 2, 3])
    );
    assert_eq!(
        cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&auth_a, &cif, b' ', b' '),
        Ok(vec![0u32, 1, 2, 3])
    );

    // Label asym id differs from auth in the mmCIF fixture: the reader
    // stores both, and the canonical identity is AUTH — "L1"/"L2" never
    // match; the PDB fixture has no label identity at all.
    let mut label = wildcard();
    label.chain_list = Some("L1,L2".to_string());
    assert_eq!(
        cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&label, &pdb, b' ', b' '),
        Ok(vec![])
    );
    assert_eq!(
        cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&label, &cif, b' ', b' '),
        Ok(vec![])
    );

    // The label identity is still recorded: chain B by auth hits water.
    let mut auth_b = wildcard();
    auth_b.chain_list = Some("B".to_string());
    for data in [&pdb, &cif] {
        assert_eq!(
            cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&auth_b, data, b' ', b' '),
            Ok(vec![4u32])
        );
    }
}

#[test]
fn bio_rows_r12_deuterium_altloc_occupancy() {
    let pdb = read(&pdb_fixture(), BioCoordinateFormat::Pdb);
    let cif = read(&mmcif_fixture(), BioCoordinateFormat::Mmcif);

    // Deuterium: the reader projects source D to H + Some(2); the
    // selection keeps D's distinct ordinal, so D (El ordinal 119) selects
    // exactly the deuterium atom and H (ordinal 1) selects nothing.
    // BIO-CID C02: the probe takes explicit ordinals (no CID parsing in
    // BIO; the lexical parser moved to the IO owner).
    let mut deu = wildcard();
    deu.element_ordinals = Some(vec![119]);
    for data in [&pdb, &cif] {
        assert_eq!(
            cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&deu, data, b' ', b' '),
            Ok(vec![2u32])
        );
    }
    let mut hyd = wildcard();
    hyd.element_ordinals = Some(vec![1]);
    for data in [&pdb, &cif] {
        assert_eq!(
            cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&hyd, data, b' ', b' '),
            Ok(vec![])
        );
    }

    // Altloc: only the GLY N carries label A in both fixtures.
    let mut alt_a = wildcard();
    alt_a.altloc_list = Some("A".to_string());
    for data in [&pdb, &cif] {
        assert_eq!(
            cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&alt_a, data, b' ', b' '),
            Ok(vec![3u32])
        );
    }

    // Occupancy: strict q > 0.5 keeps the full-occupancy rows only.
    let mut occ = wildcard();
    occ.occ_gt = Some(0.5);
    for data in [&pdb, &cif] {
        assert_eq!(
            cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&occ, data, b' ', b' '),
            Ok(vec![0u32, 1, 2])
        );
    }

    // Names: PDB stored columns trim (CA member matches stored " CA ");
    // mmCIF label_atom_id "CA" stays verbatim. Both yield the CA atom.
    let mut ca = wildcard();
    ca.atom_names = Some("CA".to_string());
    assert_eq!(
        cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&ca, &pdb, b' ', b' '),
        Ok(vec![1u32])
    );
    assert_eq!(
        cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&ca, &cif, b' ', b' '),
        Ok(vec![1u32])
    );
}

#[test]
fn bio_rows_r12_no_hit_and_seqid_bounds() {
    let pdb = read(&pdb_fixture(), BioCoordinateFormat::Pdb);
    let cif = read(&mmcif_fixture(), BioCoordinateFormat::Mmcif);

    // No hit.
    let mut none = wildcard();
    none.atom_names = Some("ZZ".to_string());
    for data in [&pdb, &cif] {
        assert_eq!(
            cosmolkit_bio::__bio_rows_probe::probe_first_rows(&none, data, b' ', b' '),
            Ok(None)
        );
        assert_eq!(
            cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&none, data, b' ', b' '),
            Ok(vec![])
        );
    }

    // Sequence bounds: auth seq 2..2 selects only GLY's atom.
    let mut band = wildcard();
    band.from_seqid = Some((2, b'*'));
    band.to_seqid = Some((2, b'*'));
    for data in [&pdb, &cif] {
        assert_eq!(
            cosmolkit_bio::__bio_rows_probe::probe_selected_atom_ids(&band, data, b' ', b' '),
            Ok(vec![3u32])
        );
    }
}
