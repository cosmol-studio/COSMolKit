//! Special source regression: all 77 fixed structure-tag cases, not a SMILES corpus.
//! Preparation is explicit; comparison never invokes RDKit or writes expectations.
use std::{
    collections::{BTreeMap, BTreeSet},
    io::{Read, Write},
    path::PathBuf,
    process::{Command, Stdio},
};

use cosmolkit_core::{
    StereoError, StructureTagParams, ValenceAssignment, assign_chiral_tags_from_structure,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer3D, CoordinateBlock, PropertyValue,
    TopologyBlock, ordered_atom_properties,
};
use cosmolkit_types::{BondDirection, BondOrder, ChiralTag, Element};
use serde::{Deserialize, Serialize};
use serde_json::Value;

const TASK: &str = "structure_tags";
const CHILD_CASE_ENV: &str = "COSMOLKIT_STRUCTURE_TAG_ORACLE_CHILD_CASE";
const NONTETRAHEDRAL_ENV: &str = "RDK_ENABLE_NONTETRAHEDRAL_STEREO";

#[derive(Debug, Deserialize, Serialize)]
struct OracleEnvironment {
    mode: String,
    #[serde(default)]
    value: Option<String>,
}

#[derive(Debug, Deserialize, Serialize)]
struct OracleRow {
    case_id: String,
    environment: OracleEnvironment,
    conf_id: i32,
    replace_existing_tags: bool,
    status: String,
    #[serde(default)]
    error_type: Option<String>,
    #[serde(default)]
    error_text: Option<String>,
    #[serde(default)]
    selected_conformer_id: Option<usize>,
    before: OracleSnapshot,
    after: OracleSnapshot,
}

#[derive(Debug, Deserialize, Serialize)]
struct OracleSnapshot {
    atom_count: usize,
    bond_count: usize,
    atoms: Vec<OracleAtom>,
    bonds: Vec<OracleBond>,
    conformers: Vec<OracleConformer>,
    #[serde(default)]
    stereochem_done: Option<i64>,
}

#[derive(Debug, Deserialize, Serialize)]
struct OracleAtom {
    index: usize,
    atomic_number: u8,
    formal_charge: i8,
    explicit_hydrogens: u8,
    implicit_hydrogens: i32,
    chiral_tag: String,
    #[serde(default)]
    chiral_permutation: Option<u32>,
    #[serde(default)]
    non_explicit_3d_chirality: Option<i64>,
    #[serde(default)]
    props: BTreeMap<String, Value>,
}

#[derive(Debug, Deserialize, Serialize)]
struct OracleBond {
    index: usize,
    begin: usize,
    end: usize,
    #[serde(rename = "type")]
    order: String,
    direction: String,
    #[serde(default)]
    unknown_stereo: Option<i64>,
    #[serde(default)]
    props: BTreeMap<String, Value>,
}

#[derive(Debug, Deserialize, Serialize)]
struct OracleConformer {
    id: usize,
    is_3d: bool,
    coordinates: Vec<[String; 3]>,
    #[serde(default)]
    props: BTreeMap<String, Value>,
}

fn data() -> PathBuf {
    cosmolkit_parity_tests_fixed::expected()
}

fn checked_snapshot() -> cosmolkit_parity_tests_fixed::special_regression::Snapshot {
    cosmolkit_parity_tests_fixed::special_regression::preflight(TASK, &data())
        .expect("complete special regression preflight before CK calls")
}

fn oracle_rows() -> Vec<OracleRow> {
    checked_snapshot()
        .rows
        .into_iter()
        .map(|row| serde_json::from_value(row).expect("complete typed structure-tag row"))
        .collect()
}

fn property_string(value: &Value) -> String {
    match value {
        Value::String(value) => value.clone(),
        Value::Number(value) => value.to_string(),
        Value::Bool(value) => value.to_string(),
        other => panic!("unsupported fixture property value {other:?}"),
    }
}

fn atom_or_bond_property(key: &str, value: &Value) -> PropertyValue {
    // Fixture generator: SetUnsignedProp for _chiralPermutation, SetIntProp
    // for other numeric fields. Chirality.cpp:3242 declares unsigned int perm;
    // :3430/:3497 writes _NonExplicit3DChirality as int. Preserve those types.
    match value {
        Value::String(value) => PropertyValue::String(value.as_str().into()),
        Value::Number(value) if key == "_chiralPermutation" => PropertyValue::UInt(
            u32::try_from(value.as_u64().expect("reference unsigned integer"))
                .expect("reference source unsigned int range"),
        ),
        Value::Number(value) => PropertyValue::Int(
            i32::try_from(value.as_i64().expect("reference signed integer"))
                .expect("reference source int range"),
        ),
        Value::Bool(value) => PropertyValue::Bool(*value),
        other => panic!("unsupported fixture property value {other:?}"),
    }
}

fn chiral_tag(name: &str) -> ChiralTag {
    ChiralTag::from_rdkit_name(name).unwrap_or_else(|| panic!("unknown chiral tag {name}"))
}

fn bond_order(name: &str) -> BondOrder {
    BondOrder::from_rdkit_name(name).unwrap_or_else(|| panic!("unknown bond order {name}"))
}

fn bond_direction(name: &str) -> BondDirection {
    BondDirection::from_rdkit_name(name).unwrap_or_else(|| panic!("unknown bond direction {name}"))
}

fn topology_from_snapshot(
    snapshot: &OracleSnapshot,
    property_predecessor: Option<&TopologyBlock>,
) -> TopologyBlock {
    assert_eq!(snapshot.atom_count, snapshot.atoms.len());
    assert_eq!(snapshot.bond_count, snapshot.bonds.len());
    let atoms = snapshot
        .atoms
        .iter()
        .map(|row| {
            let element = Element::from_atomic_number(row.atomic_number)
                .unwrap_or_else(|| panic!("unmodeled atomic number {}", row.atomic_number));
            let mut spec = AtomSpec::new(element)
                .with_formal_charge(row.formal_charge)
                .with_explicit_hydrogens(row.explicit_hydrogens)
                .with_chiral_tag(chiral_tag(&row.chiral_tag));
            if let Some(permutation) = row.chiral_permutation {
                spec = spec.with_chiral_permutation(permutation);
            }
            if let Some(predecessor) = property_predecessor {
                let predecessor = &predecessor.atoms[row.index];
                let mut inserted = BTreeSet::new();
                for (key, _) in ordered_atom_properties(predecessor) {
                    let key = std::str::from_utf8(key.as_bytes())
                        .expect("structure-tag reference property key must be UTF-8");
                    let value = row.props.get(key).unwrap_or_else(|| {
                        panic!("atom {} lost predecessor property {key}", row.index)
                    });
                    spec = spec
                        .with_prop(key, atom_or_bond_property(key, value))
                        .expect("non-empty oracle atom property key");
                    inserted.insert(key.to_owned());
                }
                // The oracle JSON carrier is a BTreeMap and cannot represent
                // Dict insertion order. Reapply the exact source write order
                // after all predecessor properties: the non-tetrahedral branch
                // writes `_chiralPermutation` first, and both branches then
                // write `_NonExplicit3DChirality` when required.
                for key in ["_chiralPermutation", "_NonExplicit3DChirality"] {
                    if inserted.contains(key) {
                        continue;
                    }
                    if let Some(value) = row.props.get(key) {
                        spec = spec
                            .with_prop(key, atom_or_bond_property(key, value))
                            .expect("non-empty generated atom property key");
                        inserted.insert(key.to_owned());
                    }
                }
                assert_eq!(
                    inserted.len(),
                    row.props.len(),
                    "atom {} has an unaccounted structure-tag property transition",
                    row.index
                );
            } else {
                for (key, value) in &row.props {
                    spec = spec
                        .with_prop(key, atom_or_bond_property(key, value))
                        .expect("non-empty fixture atom property key");
                }
            }
            let prop_value = row
                .props
                .get("_NonExplicit3DChirality")
                .and_then(Value::as_i64);
            assert_eq!(
                row.non_explicit_3d_chirality, prop_value,
                "atom {}",
                row.index
            );
            Atom::from_spec(AtomId::new(row.index), spec)
        })
        .collect();
    let bonds = snapshot
        .bonds
        .iter()
        .map(|row| {
            let mut spec = BondSpec::new(
                AtomId::new(row.begin),
                AtomId::new(row.end),
                bond_order(&row.order),
            )
            .with_direction(bond_direction(&row.direction))
            .with_unknown_stereo(row.unknown_stereo.is_some_and(|value| value != 0));
            for (key, value) in &row.props {
                spec = spec
                    .with_prop(key, atom_or_bond_property(key, value))
                    .expect("non-empty fixture bond property key");
            }
            let prop_value = row.props.get("_UnknownStereo").and_then(Value::as_i64);
            assert_eq!(row.unknown_stereo, prop_value, "bond {}", row.index);
            Bond::from_spec(BondId::new(row.index), spec)
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])
        .expect("oracle topology is structurally valid")
}

fn parse_hex_float(token: &str) -> f64 {
    match token {
        "nan" => return f64::NAN,
        "inf" => return f64::INFINITY,
        "-inf" => return f64::NEG_INFINITY,
        _ => {}
    }
    let (negative, token) = token
        .strip_prefix('-')
        .map_or((false, token), |rest| (true, rest));
    let token = token
        .strip_prefix("0x")
        .unwrap_or_else(|| panic!("coordinate is not hexadecimal: {token}"));
    let (mantissa, exponent) = token
        .split_once('p')
        .unwrap_or_else(|| panic!("coordinate has no binary exponent: {token}"));
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits = format!("{whole}{fraction}");
    let significand = u64::from_str_radix(&digits, 16).expect("valid hexadecimal significand");
    let exponent: i32 = exponent.parse().expect("valid hexadecimal exponent");
    let value = (significand as f64) * 2.0_f64.powi(exponent - 4 * fraction.len() as i32);
    if negative { -value } else { value }
}

fn coordinates_from_snapshot(snapshot: &OracleSnapshot) -> CoordinateBlock {
    CoordinateBlock {
        conformers_3d: snapshot
            .conformers
            .iter()
            .map(|row| {
                let coordinates = row
                    .coordinates
                    .iter()
                    .map(|coordinate| {
                        [
                            parse_hex_float(&coordinate[0]),
                            parse_hex_float(&coordinate[1]),
                            parse_hex_float(&coordinate[2]),
                        ]
                    })
                    .collect();
                let mut conformer = Conformer3D::new(row.id, coordinates, row.is_3d);
                for (key, value) in &row.props {
                    conformer = conformer.with_prop(key, property_string(value));
                }
                conformer
            })
            .collect(),
        ..CoordinateBlock::default()
    }
}

fn valence_from_snapshot(snapshot: &OracleSnapshot) -> ValenceAssignment {
    ValenceAssignment {
        explicit_valence: vec![0; snapshot.atoms.len()],
        implicit_hydrogens: snapshot
            .atoms
            .iter()
            .map(|atom| atom.implicit_hydrogens)
            .collect(),
    }
}

fn assert_coordinates_bitwise_eq(left: &CoordinateBlock, right: &CoordinateBlock) {
    assert_eq!(left.conformers_2d, right.conformers_2d);
    assert_eq!(left.source_coordinate_dim, right.source_coordinate_dim);
    assert_eq!(left.conformers_3d.len(), right.conformers_3d.len());
    for (left, right) in left.conformers_3d.iter().zip(&right.conformers_3d) {
        assert_eq!(left.id(), right.id());
        assert_eq!(left.is_3d(), right.is_3d());
        assert_eq!(left.props(), right.props());
        assert_eq!(left.coordinates().len(), right.coordinates().len());
        for (left, right) in left.coordinates().iter().zip(right.coordinates()) {
            assert_eq!(left.map(f64::to_bits), right.map(f64::to_bits));
        }
    }
}

fn run_oracle_row(row: &OracleRow) {
    let topology = topology_from_snapshot(&row.before, None);
    let original_topology = topology.clone();
    let coordinates = coordinates_from_snapshot(&row.before);
    let original_coordinates = coordinates.clone();
    let valence = valence_from_snapshot(&row.before);
    let result = assign_chiral_tags_from_structure(
        &topology,
        &coordinates,
        &valence,
        &StructureTagParams {
            conformer_id: row.conf_id,
            replace_existing_tags: row.replace_existing_tags,
        },
    );

    assert_eq!(
        topology, original_topology,
        "{} mutated topology input",
        row.case_id
    );
    assert_coordinates_bitwise_eq(&coordinates, &original_coordinates);

    match row.status.as_str() {
        "ok" => {
            assert!(row.error_type.is_none(), "{}", row.case_id);
            assert!(row.error_text.is_none(), "{}", row.case_id);
            let assignment = result.unwrap_or_else(|error| panic!("{}: {error}", row.case_id));
            assert_eq!(
                assignment.selected_conformer_id, row.selected_conformer_id,
                "{} selected conformer",
                row.case_id
            );
            assert_eq!(
                assignment.clear_stereochem_done,
                row.before.stereochem_done != row.after.stereochem_done,
                "{} stereochemistry clear signal",
                row.case_id
            );
            assert_eq!(
                assignment.topology,
                topology_from_snapshot(&row.after, Some(&topology)),
                "{} topology result",
                row.case_id
            );
            assert_eq!(
                row.before.conformers.len(),
                row.after.conformers.len(),
                "{} conformer count changed in oracle",
                row.case_id
            );
            for (before, after) in row.before.conformers.iter().zip(&row.after.conformers) {
                assert_eq!(before.id, after.id, "{} conformer id", row.case_id);
                assert_eq!(before.is_3d, after.is_3d, "{} dimensionality", row.case_id);
                assert_eq!(
                    before.coordinates, after.coordinates,
                    "{} coordinate rows",
                    row.case_id
                );
                assert_eq!(before.props, after.props, "{} conformer props", row.case_id);
            }
        }
        "error" => {
            let expected = match row.case_id.as_str() {
                "missing_specific_conformer" => StereoError::ConformerNotFound { requested: 8 },
                "direction_length_below_zero_tolerance" | "direction_duplicate_coordinate" => {
                    StereoError::ZeroLengthVector {
                        center: AtomId::new(0),
                        neighbor: AtomId::new(1),
                    }
                }
                "partial_mutation_before_later_exception" => StereoError::ZeroLengthVector {
                    center: AtomId::new(4),
                    neighbor: AtomId::new(5),
                },
                other => panic!("unregistered oracle error row {other}"),
            };
            assert_eq!(result, Err(expected), "{} structured error", row.case_id);
            assert!(
                row.error_type.is_some(),
                "{} missing error type",
                row.case_id
            );
            assert!(
                row.error_text.is_some(),
                "{} missing error text",
                row.case_id
            );
        }
        status => panic!("{} has unknown oracle status {status}", row.case_id),
    }
}

fn spawn_oracle_case(row: &OracleRow) -> std::process::Output {
    let mut command = Command::new(std::env::current_exe().expect("current test executable"));
    command
        .arg("--exact")
        .arg("structure_tags::structure_tags_oracle_child")
        .arg("--ignored")
        .arg("--nocapture")
        .env(CHILD_CASE_ENV, &row.case_id)
        .env_remove(NONTETRAHEDRAL_ENV)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    match row.environment.mode.as_str() {
        "unset" => assert!(row.environment.value.is_none(), "{}", row.case_id),
        "set" => {
            command.env(
                NONTETRAHEDRAL_ENV,
                row.environment
                    .value
                    .as_deref()
                    .expect("frozen environment value"),
            );
        }
        mode => panic!("{} has unknown environment mode {mode}", row.case_id),
    }
    // Pass the already validated, owned snapshot to the child, not a mutable
    // filesystem path. Concurrent pipe writing avoids output/input deadlock.
    let payload = serde_json::to_vec(row).expect("encode validated oracle row");
    let mut child = command.spawn().expect("start isolated special regression");
    let mut stdin = child.stdin.take().expect("child input");
    let writer = std::thread::spawn(move || stdin.write_all(&payload));
    let output = child
        .wait_with_output()
        .expect("collect special regression child");
    writer
        .join()
        .expect("input writer")
        .expect("write frozen oracle row");
    output
}

fn compare_case(id: &str) {
    let rows = oracle_rows();
    let row = rows
        .iter()
        .find(|row| row.case_id == id)
        .expect("registered case present");
    let output = spawn_oracle_case(row);
    assert!(
        output.status.success(),
        "{id}: {}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}
include!(concat!(env!("OUT_DIR"), "/structure_cases.rs"));
#[test]
#[ignore = "isolated case worker; parent case tests invoke it explicitly"]
fn structure_tags_oracle_child() {
    let case_id =
        std::env::var(CHILD_CASE_ENV).expect("only the parent case test invokes this worker");
    let mut payload = String::new();
    std::io::stdin()
        .read_to_string(&mut payload)
        .expect("child snapshot");
    let row: OracleRow = serde_json::from_str(&payload).expect("typed child snapshot");
    assert_eq!(row.case_id, case_id);
    run_oracle_row(&row);
}
