use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Write};
use std::path::PathBuf;

use cosmolkit_core::{
    BatchErrorMode, BatchRecord, BondDirection, BondOrder, BondStereo, ChiralTag, Molecule,
    MoleculeBatch, SdfPropertyListTarget, SmilesWriteParams,
    io::sdf::{
        SdfCoordinateMode, SdfDataset, SdfReadParams, SdfReader, SdfRecord,
        read_sdf_from_str_with_params,
    },
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct SdfReadRecord {
    smiles: Option<String>,
    case_id: String,
    dimension: String,
    format: String,
    stereo_markers: String,
    #[serde(default)]
    api: Option<String>,
    #[serde(default)]
    operation: Option<String>,
    #[serde(default)]
    sanitize: Option<bool>,
    #[serde(default)]
    remove_hs: Option<bool>,
    #[serde(default)]
    strict_parsing: Option<bool>,
    #[serde(default)]
    process_property_lists: Option<bool>,
    rdkit_ok: bool,
    sdf: Option<String>,
    atoms: Option<Vec<AtomRecord>>,
    bonds: Option<Vec<BondRecord>>,
    positions: Option<Vec<[f64; 3]>>,
    chiral_tags: Option<Vec<String>>,
    smiles_out: Option<SmilesOut>,
    properties: Option<BTreeMap<String, String>>,
    atom_properties: Option<Vec<BTreeMap<String, String>>>,
    bond_properties: Option<Vec<BTreeMap<String, String>>>,
    error: Option<String>,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
struct AtomRecord {
    atomic_num: u8,
    isotope: Option<u16>,
    formal_charge: i8,
    is_aromatic: bool,
    atom_map_num: Option<u32>,
    #[serde(default)]
    chiral_tag: Option<String>,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
struct BondRecord {
    begin: usize,
    end: usize,
    bond_type: String,
    is_aromatic: bool,
    direction: String,
    stereo: String,
    stereo_atoms: Vec<usize>,
}

#[derive(Debug, Deserialize)]
struct SmilesOut {
    canonical: String,
    noncanonical: String,
}

fn golden_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/sdf_read.jsonl")
}

fn ensure_golden_exists() {
    let path = golden_path();
    assert!(
        path.exists(),
        "missing RDKit SDF read golden: {}. Generate it before running tests:\n\
         uv sync --group dev && .venv/bin/python tests/scripts/gen_all_rdkit_goldens.py --python .venv/bin/python --clean --jobs 4",
        path.display()
    );
}

fn load_golden() -> Vec<SdfReadRecord> {
    ensure_golden_exists();
    let path = golden_path();
    let file = File::open(&path).expect("should read SDF read golden");
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(idx, line)| {
            let line = line.unwrap_or_else(|error| {
                panic!(
                    "failed to read {} line {}: {error}",
                    path.display(),
                    idx + 1
                )
            });
            serde_json::from_str(&line).unwrap_or_else(|error| {
                panic!(
                    "failed to parse {} line {}: {error}",
                    path.display(),
                    idx + 1
                )
            })
        })
        .collect()
}

fn record_label(record: &SdfReadRecord) -> &str {
    record.smiles.as_deref().unwrap_or("<no-smiles>")
}

fn read_params(record: &SdfReadRecord) -> SdfReadParams {
    SdfReadParams {
        sanitize: record.sanitize.unwrap_or(true),
        remove_hs: record.remove_hs.unwrap_or(true),
        strict_parsing: record.strict_parsing.unwrap_or(true),
        process_property_lists: record.process_property_lists.unwrap_or(true),
        ..Default::default()
    }
}

fn bond_order_name(order: BondOrder) -> &'static str {
    order.rdkit_name()
}

fn bond_direction_name(direction: BondDirection) -> &'static str {
    match direction {
        BondDirection::None => "NONE",
        BondDirection::BeginWedge => "BEGINWEDGE",
        BondDirection::BeginDash => "BEGINDASH",
        BondDirection::EndUpRight => "ENDUPRIGHT",
        BondDirection::EndDownRight => "ENDDOWNRIGHT",
        BondDirection::EitherDouble => "EITHERDOUBLE",
        BondDirection::Unknown => "UNKNOWN",
    }
}

fn bond_stereo_name(stereo: BondStereo) -> &'static str {
    match stereo {
        BondStereo::None => "STEREONONE",
        BondStereo::Any => "STEREOANY",
        BondStereo::Z => "STEREOZ",
        BondStereo::E => "STEREOE",
        BondStereo::Cis => "STEREOZ",
        BondStereo::Trans => "STEREOE",
        BondStereo::AtropCw => "STEREOATROPCW",
        BondStereo::AtropCcw => "STEREOATROPCCW",
    }
}

fn chiral_tag_name(tag: ChiralTag) -> &'static str {
    tag.rdkit_name()
}

fn direct_sdf_result(record: &SdfReadRecord) -> Result<SdfRecord, String> {
    let sdf = record
        .sdf
        .as_ref()
        .ok_or_else(|| "golden row has no SDF text".to_string())?;
    read_sdf_from_str_with_params(sdf, read_params(record)).map_err(|error| error.to_string())
}

fn forward_reader_result(record: &SdfReadRecord) -> Result<SdfRecord, String> {
    let sdf = record
        .sdf
        .as_ref()
        .ok_or_else(|| "golden row has no SDF text".to_string())?;
    let mut reader = SdfReader::with_params(Cursor::new(sdf.as_bytes()), read_params(record));
    reader
        .next_record()
        .map_err(|error| error.to_string())?
        .ok_or_else(|| "SdfReader returned no records".to_string())
}

fn indexed_dataset_result(record: &SdfReadRecord) -> Result<SdfRecord, String> {
    let sdf = record
        .sdf
        .as_ref()
        .ok_or_else(|| "golden row has no SDF text".to_string())?;
    let params = read_params(record);
    let mut temp = tempfile::NamedTempFile::new().expect("should create temp SDF file");
    temp.write_all(sdf.as_bytes())
        .expect("should write temp SDF file");
    let dataset = SdfDataset::open_with_params(temp.path(), params).map_err(|e| e.to_string())?;
    dataset
        .record_with_params(0, params)
        .map_err(|error| error.to_string())
}

fn assert_result_matches_rdkit(
    actual: Result<SdfRecord, String>,
    record: &SdfReadRecord,
    row_idx: usize,
    api_under_test: &str,
) {
    if !record.rdkit_ok {
        assert!(
            record.error.is_some(),
            "row {} {} ({}) is rdkit not ok but has no error",
            row_idx + 1,
            record.case_id,
            record_label(record)
        );
        if record.sdf.is_some() {
            assert!(
                actual.is_err(),
                "{api_under_test} unexpectedly parsed RDKit failure row {} case {} ({})",
                row_idx + 1,
                record.case_id,
                record_label(record)
            );
        }
        return;
    }

    let parsed = actual.unwrap_or_else(|error| {
        panic!(
            "{api_under_test} failed at row {} case {} ({}): {error}",
            row_idx + 1,
            record.case_id,
            record_label(record)
        )
    });
    assert_record_matches_rdkit(&parsed, record, row_idx, api_under_test);
}

fn assert_record_matches_rdkit(
    actual: &SdfRecord,
    record: &SdfReadRecord,
    row_idx: usize,
    api_under_test: &str,
) {
    assert_molecule_matches_rdkit(&actual.molecule, record, row_idx, api_under_test);
    assert_sdf_data_fields_match(actual, record, row_idx, api_under_test);
    assert_sdf_property_lists_match(&actual.molecule, record, row_idx, api_under_test);
}

fn assert_molecule_matches_rdkit(
    molecule: &Molecule,
    record: &SdfReadRecord,
    row_idx: usize,
    api_under_test: &str,
) {
    let expected_atoms = record.atoms.as_ref().expect("rdkit_ok row has atoms");
    let expected_bonds = record.bonds.as_ref().expect("rdkit_ok row has bonds");

    let actual_atoms = molecule
        .atoms()
        .iter()
        .map(|atom| AtomRecord {
            atomic_num: atom.atomic_number(),
            isotope: atom.isotope(),
            formal_charge: atom.formal_charge(),
            is_aromatic: atom.is_aromatic(),
            atom_map_num: atom.atom_map(),
            chiral_tag: Some(chiral_tag_name(atom.chiral_tag()).to_owned()),
        })
        .collect::<Vec<_>>();
    assert_eq!(
        actual_atoms,
        *expected_atoms,
        "{api_under_test} atom field mismatch at row {} case {} ({})",
        row_idx + 1,
        record.case_id,
        record_label(record)
    );

    let actual_bonds = molecule
        .bonds()
        .iter()
        .map(|bond| BondRecord {
            begin: bond.begin().index(),
            end: bond.end().index(),
            bond_type: bond_order_name(bond.order()).to_owned(),
            is_aromatic: bond.is_aromatic(),
            direction: bond_direction_name(bond.direction()).to_owned(),
            stereo: bond_stereo_name(bond.stereo()).to_owned(),
            stereo_atoms: bond
                .stereo_atoms()
                .map(|atoms| atoms.into_iter().map(|atom| atom.index()).collect())
                .unwrap_or_default(),
        })
        .collect::<Vec<_>>();
    assert_eq!(
        actual_bonds,
        *expected_bonds,
        "{api_under_test} bond field mismatch at row {} case {} ({})",
        row_idx + 1,
        record.case_id,
        record_label(record)
    );

    assert_coordinates_match_rdkit(molecule, record, row_idx, api_under_test);
    assert_chirality_matches_rdkit(molecule, record, row_idx, api_under_test);
    assert_smiles_matches_rdkit(molecule, record, row_idx, api_under_test);
    assert_molfile_properties_match_rdkit(molecule, record, row_idx, api_under_test);
    assert_public_atom_bond_properties_match_rdkit(molecule, record, row_idx, api_under_test);
}

fn assert_coordinates_match_rdkit(
    molecule: &Molecule,
    record: &SdfReadRecord,
    row_idx: usize,
    api_under_test: &str,
) {
    let expected_positions = record
        .positions
        .as_ref()
        .expect("rdkit_ok row should have positions");
    if record.dimension == "2D" {
        let coords = molecule.coordinates_2d().unwrap_or_else(|| {
            panic!(
                "{api_under_test} row {} {} should preserve 2D coords",
                row_idx + 1,
                record.case_id
            )
        });
        assert_eq!(coords.len(), expected_positions.len());
        for (atom_idx, (actual, expected)) in coords.iter().zip(expected_positions).enumerate() {
            assert!(
                (actual[0] - expected[0]).abs() <= 1e-12
                    && (actual[1] - expected[1]).abs() <= 1e-12
                    && expected[2].abs() <= 1e-12,
                "{api_under_test} 2D coordinate mismatch at row {} case {} atom {} ({})",
                row_idx + 1,
                record.case_id,
                atom_idx,
                record_label(record)
            );
        }
    } else if record.dimension == "3D" {
        let coords = molecule
            .conformers_3d()
            .first()
            .map(|c| c.coordinates())
            .unwrap_or_else(|| {
                panic!(
                    "{api_under_test} row {} {} should preserve 3D coords",
                    row_idx + 1,
                    record.case_id
                )
            });
        assert_eq!(coords.len(), expected_positions.len());
        for (atom_idx, (actual, expected)) in coords.iter().zip(expected_positions).enumerate() {
            assert!(
                (actual[0] - expected[0]).abs() <= 1e-12
                    && (actual[1] - expected[1]).abs() <= 1e-12
                    && (actual[2] - expected[2]).abs() <= 1e-12,
                "{api_under_test} 3D coordinate mismatch at row {} case {} atom {} ({})",
                row_idx + 1,
                record.case_id,
                atom_idx,
                record_label(record)
            );
        }
    }
}

fn assert_chirality_matches_rdkit(
    molecule: &Molecule,
    record: &SdfReadRecord,
    row_idx: usize,
    api_under_test: &str,
) {
    let expected_tags = record
        .chiral_tags
        .as_ref()
        .expect("rdkit_ok row should have chiral tags");
    let actual_tags = molecule
        .atoms()
        .iter()
        .map(|atom| chiral_tag_name(atom.chiral_tag()).to_owned())
        .collect::<Vec<_>>();
    assert_eq!(
        actual_tags,
        *expected_tags,
        "{api_under_test} chirality mismatch at row {} case {} ({})",
        row_idx + 1,
        record.case_id,
        record_label(record)
    );
}

fn assert_smiles_matches_rdkit(
    molecule: &Molecule,
    record: &SdfReadRecord,
    row_idx: usize,
    api_under_test: &str,
) {
    let expected = record
        .smiles_out
        .as_ref()
        .expect("rdkit_ok row should have SMILES output");

    match molecule.to_smiles_with_params(&SmilesWriteParams {
        canonical: true,
        do_isomeric_smiles: true,
        ..Default::default()
    }) {
        Ok(canonical) => assert_eq!(
            canonical,
            expected.canonical,
            "{api_under_test} canonical SMILES mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record_label(record)
        ),
        Err(error) => assert!(
            error.to_string().contains("unsupported path"),
            "{api_under_test} canonical SMILES write failed without explicit unsupported error at row {} case {} ({}): {error}",
            row_idx + 1,
            record.case_id,
            record_label(record)
        ),
    }

    let noncanonical = molecule
        .to_smiles_with_params(&SmilesWriteParams {
            canonical: false,
            do_isomeric_smiles: true,
            ..Default::default()
        })
        .unwrap_or_else(|error| {
            panic!(
                "{api_under_test} noncanonical SMILES write failed at row {} case {} ({}): {error}",
                row_idx + 1,
                record.case_id,
                record_label(record)
            )
        });
    assert_eq!(
        noncanonical,
        expected.noncanonical,
        "{api_under_test} noncanonical SMILES mismatch at row {} case {} ({})",
        row_idx + 1,
        record.case_id,
        record_label(record)
    );
}

fn assert_molfile_properties_match_rdkit(
    molecule: &Molecule,
    record: &SdfReadRecord,
    row_idx: usize,
    api_under_test: &str,
) {
    let Some(expected_props) = record.properties.as_ref() else {
        return;
    };
    let expected_name = expected_props.get("_Name").map(String::as_str);
    assert_eq!(
        molecule.properties().name(),
        expected_name,
        "{api_under_test} _Name mismatch at row {} case {} ({})",
        row_idx + 1,
        record.case_id,
        record_label(record)
    );
    for (rdkit_key, cosmolkit_key) in [
        ("_MolFileComments", "_MolFileComments"),
        ("_MolFileChiralFlag", "_MolFileChiralFlag"),
        ("_MolFileInfo", "_MolFileInfoLine"),
    ] {
        if let Some(expected) = expected_props.get(rdkit_key) {
            assert_eq!(
                molecule.prop(cosmolkit_key),
                Some(expected.as_str()),
                "{api_under_test} property {rdkit_key}/{cosmolkit_key} mismatch at row {} case {} ({})",
                row_idx + 1,
                record.case_id,
                record_label(record)
            );
        }
    }
}

fn assert_sdf_data_fields_match(
    actual: &SdfRecord,
    record: &SdfReadRecord,
    row_idx: usize,
    api_under_test: &str,
) {
    let Some(expected_props) = record.properties.as_ref() else {
        return;
    };
    let expected_fields = expected_sdf_fields(record);
    assert_eq!(
        actual.data_fields,
        expected_fields,
        "{api_under_test} SDF data fields mismatch at row {} case {} ({})",
        row_idx + 1,
        record.case_id,
        record_label(record)
    );
    assert_eq!(
        actual.molecule.properties().sdf_data_fields(),
        expected_fields.as_slice(),
        "{api_under_test} molecule SDF data field storage mismatch at row {} case {} ({})",
        row_idx + 1,
        record.case_id,
        record_label(record)
    );
    let mut last_field_value_by_key = BTreeMap::<String, String>::new();
    for (key, expected_value) in expected_fields {
        last_field_value_by_key.insert(key, expected_value);
    }
    for (key, expected_value) in last_field_value_by_key {
        assert_eq!(
            expected_props.get(&key).map(String::as_str),
            Some(expected_value.as_str()),
            "RDKit golden properties do not mirror final data field value {key} at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record_label(record)
        );
        assert_eq!(
            actual.molecule.prop(&key),
            Some(expected_value.as_str()),
            "{api_under_test} molecule property for final SDF data field value {key} mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record_label(record)
        );
    }
}

fn expected_sdf_fields(record: &SdfReadRecord) -> Vec<(String, String)> {
    let Some(sdf) = record.sdf.as_deref() else {
        return Vec::new();
    };
    let Some(data_start) = sdf.find("\n>") else {
        return Vec::new();
    };
    let mut fields = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_lines: Vec<String> = Vec::new();
    for raw_line in sdf[data_start + 1..].lines() {
        if raw_line == "$$$$" {
            break;
        }
        if raw_line.trim_start().starts_with('>') {
            if let Some(name) = current_name.take() {
                fields.push((name, current_lines.join("\n")));
                current_lines.clear();
            }
            current_name = extract_sdf_data_label(raw_line);
        } else if raw_line.is_empty() {
            if let Some(name) = current_name.take() {
                fields.push((name, current_lines.join("\n")));
                current_lines.clear();
            }
        } else if current_name.is_some() {
            current_lines.push(raw_line.strip_suffix('\r').unwrap_or(raw_line).to_string());
        }
    }
    if let Some(name) = current_name {
        fields.push((name, current_lines.join("\n")));
    }
    fields
}

fn extract_sdf_data_label(line: &str) -> Option<String> {
    let after_arrow = line.trim_start().strip_prefix('>')?;
    let begin = after_arrow.find('<')?;
    let end = after_arrow.rfind('>')?;
    (end > begin + 1).then(|| after_arrow[begin + 1..end].to_string())
}

fn assert_public_atom_bond_properties_match_rdkit(
    molecule: &Molecule,
    record: &SdfReadRecord,
    row_idx: usize,
    api_under_test: &str,
) {
    if let Some(expected_atom_props) = record.atom_properties.as_ref() {
        assert_eq!(expected_atom_props.len(), molecule.num_atoms());
        for (atom_idx, expected_props) in expected_atom_props.iter().enumerate() {
            for (key, expected_value) in public_modeled_props(expected_props) {
                assert_eq!(
                    molecule.atoms()[atom_idx].prop(key),
                    Some(expected_value),
                    "{api_under_test} atom property {key} mismatch at row {} case {} atom {} ({})",
                    row_idx + 1,
                    record.case_id,
                    atom_idx,
                    record_label(record)
                );
            }
        }
    }
    if let Some(expected_bond_props) = record.bond_properties.as_ref() {
        assert_eq!(expected_bond_props.len(), molecule.num_bonds());
        for (bond_idx, expected_props) in expected_bond_props.iter().enumerate() {
            for (key, expected_value) in public_modeled_props(expected_props) {
                assert_eq!(
                    molecule.bonds()[bond_idx].prop(key),
                    Some(expected_value),
                    "{api_under_test} bond property {key} mismatch at row {} case {} bond {} ({})",
                    row_idx + 1,
                    record.case_id,
                    bond_idx,
                    record_label(record)
                );
            }
        }
    }
}

fn public_modeled_props(props: &BTreeMap<String, String>) -> impl Iterator<Item = (&str, &str)> {
    props.iter().filter_map(|(key, value)| {
        (!key.starts_with('_') && key != "__computedProps")
            .then_some((key.as_str(), value.as_str()))
    })
}

fn assert_sdf_property_lists_match(
    molecule: &Molecule,
    record: &SdfReadRecord,
    row_idx: usize,
    api_under_test: &str,
) {
    let expected = expected_property_lists_from_rdkit(record);
    let actual = molecule
        .properties()
        .sdf_property_lists()
        .iter()
        .map(|list| {
            let target = match list.target() {
                SdfPropertyListTarget::Atom => "atom",
                SdfPropertyListTarget::Bond => "bond",
            };
            (
                target.to_string(),
                list.name().to_string(),
                list.values().to_vec(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        actual,
        expected,
        "{api_under_test} SDF property-list mismatch at row {} case {} ({})",
        row_idx + 1,
        record.case_id,
        record_label(record)
    );
}

fn expected_property_lists_from_rdkit(
    record: &SdfReadRecord,
) -> Vec<(String, String, Vec<Option<String>>)> {
    if !record.process_property_lists.unwrap_or(true) {
        return Vec::new();
    }
    let mut out = Vec::new();
    for (name, _) in expected_sdf_fields(record) {
        if let Some(prop_name) = name.strip_prefix("atom.prop.") {
            push_expected_list(
                &mut out,
                "atom",
                prop_name,
                record.atom_properties.as_deref(),
            );
        } else if let Some(prop_name) = name.strip_prefix("atom.iprop.") {
            push_expected_list(
                &mut out,
                "atom",
                prop_name,
                record.atom_properties.as_deref(),
            );
        } else if let Some(prop_name) = name.strip_prefix("atom.dprop.") {
            push_expected_list(
                &mut out,
                "atom",
                prop_name,
                record.atom_properties.as_deref(),
            );
        } else if let Some(prop_name) = name.strip_prefix("atom.bprop.") {
            push_expected_list(
                &mut out,
                "atom",
                prop_name,
                record.atom_properties.as_deref(),
            );
        } else if let Some(prop_name) = name.strip_prefix("bond.prop.") {
            push_expected_list(
                &mut out,
                "bond",
                prop_name,
                record.bond_properties.as_deref(),
            );
        } else if let Some(prop_name) = name.strip_prefix("bond.iprop.") {
            push_expected_list(
                &mut out,
                "bond",
                prop_name,
                record.bond_properties.as_deref(),
            );
        } else if let Some(prop_name) = name.strip_prefix("bond.dprop.") {
            push_expected_list(
                &mut out,
                "bond",
                prop_name,
                record.bond_properties.as_deref(),
            );
        } else if let Some(prop_name) = name.strip_prefix("bond.bprop.") {
            push_expected_list(
                &mut out,
                "bond",
                prop_name,
                record.bond_properties.as_deref(),
            );
        }
    }
    out
}

fn push_expected_list(
    out: &mut Vec<(String, String, Vec<Option<String>>)>,
    target: &str,
    prop_name: &str,
    props: Option<&[BTreeMap<String, String>]>,
) {
    let Some(props) = props else {
        return;
    };
    if !props.iter().any(|row| row.contains_key(prop_name)) {
        return;
    }
    out.push((
        target.to_string(),
        prop_name.to_string(),
        props
            .iter()
            .map(|row| row.get(prop_name).cloned())
            .collect(),
    ));
}

fn assert_case_matrix(records: &[SdfReadRecord]) {
    let expected_cases = [
        "2d_v2000_with_markers",
        "2d_v2000_coords_only",
        "2d_v3000_with_markers",
        "2d_v3000_coords_only",
        "3d_v2000_with_markers",
        "3d_v2000_coords_only",
        "3d_v3000_with_markers",
        "3d_v3000_coords_only",
    ];
    assert!(
        records.iter().any(|record| record.rdkit_ok),
        "SDF read golden should include at least one successful RDKit case"
    );
    for case_id in expected_cases {
        assert!(
            records.iter().any(|record| record.case_id == case_id),
            "SDF read golden is missing case {case_id}"
        );
    }
    for api in ["ForwardSDMolSupplier", "SDMolSupplier"] {
        assert!(
            records
                .iter()
                .any(|record| record.api.as_deref() == Some(api)),
            "SDF read golden is missing API {api}"
        );
    }
    for format in ["V2000", "V3000"] {
        assert!(
            records.iter().any(|record| record.format == format),
            "SDF read golden is missing format {format}"
        );
    }
    for sanitize in [true, false] {
        for remove_hs in [true, false] {
            for strict_parsing in [true, false] {
                for process_property_lists in [true, false] {
                    assert!(
                        records.iter().any(|record| {
                            record.sanitize == Some(sanitize)
                                && record.remove_hs == Some(remove_hs)
                                && record.strict_parsing == Some(strict_parsing)
                                && record.process_property_lists == Some(process_property_lists)
                        }),
                        "SDF read golden is missing parameter row sanitize={sanitize} remove_hs={remove_hs} strict_parsing={strict_parsing} process_property_lists={process_property_lists}"
                    );
                }
            }
        }
    }
    assert!(
        records.iter().any(|record| record
            .sdf
            .as_deref()
            .is_some_and(|sdf| sdf.contains(">  <ID>"))),
        "SDF read golden must cover data fields"
    );
    assert!(
        records.iter().any(|record| record
            .sdf
            .as_deref()
            .is_some_and(|sdf| sdf.contains(">  <atom.prop.note>"))),
        "SDF read golden must cover atom property lists"
    );
    assert!(
        records
            .iter()
            .any(|record| record.operation.as_deref() == Some("malformed")),
        "SDF read golden must cover malformed records"
    );
    assert!(
        records.iter().any(|record| record.dimension == "3D"
            && record.stereo_markers == "coords_only"
            && record
                .chiral_tags
                .as_ref()
                .is_some_and(|tags| tags.iter().any(|tag| tag != "CHI_UNSPECIFIED"))),
        "SDF read golden must cover coordinate-inferred 3D chirality"
    );
}

#[test]
fn sdf_read_golden_covers_expected_case_matrix() {
    let records = load_golden();
    assert_case_matrix(&records);
}

#[test]
fn sdf_read_representative_apis_and_parameters_match_rdkit() {
    let records = load_golden();
    for (row_idx, record) in quick_sdf_parity_rows(&records) {
        assert_result_matches_rdkit(
            direct_sdf_result(record),
            record,
            row_idx,
            "read_sdf_from_str_with_params",
        );
        match record.api.as_deref() {
            Some("ForwardSDMolSupplier") => assert_result_matches_rdkit(
                forward_reader_result(record),
                record,
                row_idx,
                "SdfReader::with_params::next_record",
            ),
            Some("SDMolSupplier") => assert_result_matches_rdkit(
                indexed_dataset_result(record),
                record,
                row_idx,
                "SdfDataset::record_with_params",
            ),
            None => {}
            Some(other) => panic!(
                "unsupported SDF golden API {other} at row {} case {}",
                row_idx + 1,
                record.case_id
            ),
        }
    }
}

#[test]
#[ignore = "full RDKit SDF reader API/parameter matrix is expensive; run explicitly for exhaustive parity"]
fn sdf_read_all_apis_and_parameters_match_rdkit() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        assert_result_matches_rdkit(
            direct_sdf_result(record),
            record,
            row_idx,
            "read_sdf_from_str_with_params",
        );
        match record.api.as_deref() {
            Some("ForwardSDMolSupplier") => assert_result_matches_rdkit(
                forward_reader_result(record),
                record,
                row_idx,
                "SdfReader::with_params::next_record",
            ),
            Some("SDMolSupplier") => assert_result_matches_rdkit(
                indexed_dataset_result(record),
                record,
                row_idx,
                "SdfDataset::record_with_params",
            ),
            None => {}
            Some(other) => panic!(
                "unsupported SDF golden API {other} at row {} case {}",
                row_idx + 1,
                record.case_id
            ),
        }
    }
}

#[test]
fn molecule_batch_sdf_reader_and_dataset_paths_match_rdkit_defaults() {
    let records = load_golden();
    let forward_rows = default_success_rows(&records, "ForwardSDMolSupplier", 32);
    let forward_sdf = forward_rows
        .iter()
        .map(|(_, record)| record.sdf.as_deref().unwrap())
        .collect::<String>();
    let batch = MoleculeBatch::read_sdf_records_from_reader_with_options(
        Cursor::new(forward_sdf.as_bytes()),
        SdfCoordinateMode::Preserve,
        BatchErrorMode::KeepErrors,
        Some(1),
        Some(false),
    )
    .expect("batch reader should parse default ForwardSDMolSupplier rows");
    assert_batch_matches_rows(
        &batch,
        &forward_rows,
        "MoleculeBatch::read_sdf_records_from_reader_with_options",
    );

    let dataset_rows = default_success_rows(&records, "SDMolSupplier", 32);
    let dataset_sdf = dataset_rows
        .iter()
        .map(|(_, record)| record.sdf.as_deref().unwrap())
        .collect::<String>();
    let mut temp = tempfile::NamedTempFile::new().expect("should create temp SDF dataset");
    temp.write_all(dataset_sdf.as_bytes())
        .expect("should write temp SDF dataset");
    let dataset = SdfDataset::open_with_params(temp.path(), SdfReadParams::default())
        .expect("should open temp SDF dataset");
    let batch = MoleculeBatch::read_sdf_dataset_with_params_and_options(
        &dataset,
        SdfReadParams::default(),
        BatchErrorMode::KeepErrors,
        Some(1),
        Some(false),
    )
    .expect("batch dataset reader should parse default SDMolSupplier rows");
    assert_batch_matches_rows(
        &batch,
        &dataset_rows,
        "MoleculeBatch::read_sdf_dataset_with_params_and_options",
    );
}

fn quick_sdf_parity_rows(records: &[SdfReadRecord]) -> Vec<(usize, &SdfReadRecord)> {
    let mut rows = Vec::<(usize, &SdfReadRecord)>::new();
    let mut seen = std::collections::BTreeSet::<usize>::new();
    let mut push_first = |predicate: &dyn Fn(&SdfReadRecord) -> bool| {
        if let Some((idx, record)) = records
            .iter()
            .enumerate()
            .find(|(idx, record)| !seen.contains(idx) && predicate(record))
        {
            seen.insert(idx);
            rows.push((idx, record));
        }
    };

    for api in [None, Some("ForwardSDMolSupplier"), Some("SDMolSupplier")] {
        push_first(&|record| {
            record.rdkit_ok
                && record.api.as_deref() == api
                && read_params(record) == SdfReadParams::default()
        });
    }
    for sanitize in [false, true] {
        for remove_hs in [false, true] {
            push_first(&|record| {
                record.rdkit_ok
                    && record.sanitize == Some(sanitize)
                    && record.remove_hs == Some(remove_hs)
                    && record.strict_parsing == Some(true)
                    && record.process_property_lists == Some(true)
            });
        }
    }
    push_first(&|record| record.rdkit_ok && record.process_property_lists == Some(false));
    push_first(&|record| record.rdkit_ok && record.strict_parsing == Some(false));
    push_first(&|record| {
        record
            .sdf
            .as_deref()
            .is_some_and(|sdf| sdf.contains(">  <atom.prop.note>"))
    });
    push_first(&|record| !record.rdkit_ok);

    assert!(
        rows.len() >= 8,
        "quick SDF parity should cover default APIs, parameter variants, property lists, and failures"
    );
    rows
}

fn default_success_rows<'a>(
    records: &'a [SdfReadRecord],
    api: &str,
    limit: usize,
) -> Vec<(usize, &'a SdfReadRecord)> {
    records
        .iter()
        .enumerate()
        .filter(|(_, record)| {
            record.rdkit_ok
                && record.api.as_deref() == Some(api)
                && read_params(record) == SdfReadParams::default()
                && record.sdf.is_some()
        })
        .take(limit)
        .collect()
}

fn assert_batch_matches_rows(
    batch: &MoleculeBatch,
    rows: &[(usize, &SdfReadRecord)],
    api_under_test: &str,
) {
    assert_eq!(
        batch.len(),
        rows.len(),
        "{api_under_test} batch length mismatch"
    );
    for ((row_idx, record), batch_record) in rows.iter().zip(batch.iter()) {
        let BatchRecord::Molecule(molecule) = batch_record else {
            panic!(
                "{api_under_test} missing valid molecule at row {} case {} ({})",
                row_idx + 1,
                record.case_id,
                record_label(record)
            );
        };
        assert_molecule_matches_rdkit(molecule, record, *row_idx, api_under_test);
    }
}
