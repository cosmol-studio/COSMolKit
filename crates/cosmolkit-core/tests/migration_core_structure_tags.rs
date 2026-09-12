use std::{
    collections::{BTreeMap, BTreeSet},
    process::Command,
};

use cosmolkit_core::{
    StereoError, StructureTagParams, ValenceAssignment, assign_chiral_tags_from_structure,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer3D, CoordinateBlock,
    TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, ChiralTag, Element};
use serde::Deserialize;
use serde_json::Value;

const FIXTURE: &str = include_str!(
    "../../../testdata/stereo/fixtures/assign_atom_chiral_tags_from_structure_cases.json"
);
const ORACLE: &str = include_str!(
    "../../../testdata/stereo/expected/rdkit/smiles_small/assign_atom_chiral_tags_from_structure.jsonl"
);
const CHILD_CASE_ENV: &str = "COSMOLKIT_STRUCTURE_TAG_ORACLE_CHILD_CASE";
const NUMERIC_CHILD_ENV: &str = "COSMOLKIT_STRUCTURE_TAG_NUMERIC_CHILD";
const NONTETRAHEDRAL_ENV: &str = "RDK_ENABLE_NONTETRAHEDRAL_STEREO";

#[derive(Debug, Deserialize)]
struct OracleEnvironment {
    mode: String,
    #[serde(default)]
    value: Option<String>,
}

#[derive(Debug, Deserialize)]
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

#[derive(Debug, Deserialize)]
struct OracleSnapshot {
    atom_count: usize,
    bond_count: usize,
    atoms: Vec<OracleAtom>,
    bonds: Vec<OracleBond>,
    conformers: Vec<OracleConformer>,
    #[serde(default)]
    stereochem_done: Option<i64>,
}

#[derive(Debug, Deserialize)]
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

#[derive(Debug, Deserialize)]
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

#[derive(Debug, Deserialize)]
struct OracleConformer {
    id: usize,
    is_3d: bool,
    coordinates: Vec<[String; 3]>,
    #[serde(default)]
    props: BTreeMap<String, Value>,
}

fn oracle_rows() -> Vec<OracleRow> {
    ORACLE
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| serde_json::from_str(line).expect("valid structure-tag oracle row"))
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

fn chiral_tag(name: &str) -> ChiralTag {
    ChiralTag::from_rdkit_name(name).unwrap_or_else(|| panic!("unknown chiral tag {name}"))
}

fn bond_order(name: &str) -> BondOrder {
    BondOrder::from_rdkit_name(name).unwrap_or_else(|| panic!("unknown bond order {name}"))
}

fn bond_direction(name: &str) -> BondDirection {
    BondDirection::from_rdkit_name(name).unwrap_or_else(|| panic!("unknown bond direction {name}"))
}

fn topology_from_snapshot(snapshot: &OracleSnapshot) -> TopologyBlock {
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
            for (key, value) in &row.props {
                spec = spec
                    .with_prop(key, property_string(value))
                    .expect("non-empty fixture atom property key");
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
                    .with_prop(key, property_string(value))
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
    let topology = topology_from_snapshot(&row.before);
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
                topology_from_snapshot(&row.after),
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

fn spawn_oracle_case(row: &OracleRow) {
    let mut command = Command::new(std::env::current_exe().expect("current test executable"));
    command
        .arg("--exact")
        .arg("structure_tags_oracle_child")
        .arg("--nocapture")
        .env(CHILD_CASE_ENV, &row.case_id)
        .env_remove(NONTETRAHEDRAL_ENV);
    match row.environment.mode.as_str() {
        "unset" => assert!(row.environment.value.is_none(), "{}", row.case_id),
        "set" => {
            command.env(
                NONTETRAHEDRAL_ENV,
                row.environment
                    .value
                    .as_deref()
                    .unwrap_or_else(|| panic!("{} missing environment value", row.case_id)),
            );
        }
        mode => panic!("{} has unknown environment mode {mode}", row.case_id),
    }
    let output = command.output().expect("run isolated structure-tag case");
    assert!(
        output.status.success(),
        "isolated case {} failed\nstdout:\n{}\nstderr:\n{}",
        row.case_id,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn simple_star(center_atomic_number: u8, explicit_hydrogens: u8) -> TopologyBlock {
    let center = AtomSpec::new(
        Element::from_atomic_number(center_atomic_number).expect("modeled center element"),
    )
    .with_explicit_hydrogens(explicit_hydrogens);
    let atoms = [
        center,
        AtomSpec::new(Element::F),
        AtomSpec::new(Element::CL),
        AtomSpec::new(Element::BR),
    ]
    .into_iter()
    .enumerate()
    .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
    .collect();
    let bonds = (1..4)
        .map(|neighbor| {
            Bond::from_spec(
                BondId::new(neighbor - 1),
                BondSpec::new(AtomId::new(0), AtomId::new(neighbor), BondOrder::Single),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn simple_coordinates(rows: Vec<[f64; 3]>) -> CoordinateBlock {
    CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(0, rows, true)],
        ..CoordinateBlock::default()
    }
}

#[test]
fn defaults_and_committed_fixture_inventory_are_exact() {
    assert_eq!(
        StructureTagParams::default(),
        StructureTagParams {
            conformer_id: -1,
            replace_existing_tags: true,
        }
    );

    let fixture: Value = serde_json::from_str(FIXTURE).expect("valid committed fixture");
    assert_eq!(fixture["schema_version"], 1);
    assert_eq!(fixture["reference"]["implementation"], "RDKit");
    assert_eq!(fixture["reference"]["version"], "2026.3.1");
    assert_eq!(fixture["defaults"]["conf_id"], -1);
    assert_eq!(fixture["defaults"]["replace_existing_tags"], true);
    let primary = fixture["cases"].as_array().expect("fixture cases");
    let switch = fixture["octahedral_switch_cases"]
        .as_array()
        .expect("octahedral cases");
    assert_eq!(primary.len(), 47);
    assert_eq!(switch.len(), 30);

    let fixture_ids = primary
        .iter()
        .chain(switch)
        .map(|case| {
            case["case_id"]
                .as_str()
                .expect("fixture case id")
                .to_owned()
        })
        .collect::<BTreeSet<_>>();
    let rows = oracle_rows();
    let oracle_ids = rows
        .iter()
        .map(|row| row.case_id.clone())
        .collect::<BTreeSet<_>>();
    assert_eq!(rows.len(), 77);
    assert_eq!(oracle_ids.len(), 77);
    assert_eq!(fixture_ids, oracle_ids);
}

#[test]
fn rdkit_oracle_all_rows_run_in_controlled_processes_without_filtering() {
    let rows = oracle_rows();
    assert_eq!(rows.len(), 77);
    for row in &rows {
        spawn_oracle_case(row);
    }
}

#[test]
fn structure_tags_oracle_child() {
    let Ok(case_id) = std::env::var(CHILD_CASE_ENV) else {
        return;
    };
    let rows = oracle_rows();
    let matching = rows
        .iter()
        .filter(|row| row.case_id == case_id)
        .collect::<Vec<_>>();
    assert_eq!(matching.len(), 1, "child case id must be unique");
    run_oracle_row(matching[0]);
}

#[test]
fn model_boundaries_return_precise_topology_valence_and_coordinate_errors() {
    let invalid_topology = TopologyBlock {
        atoms: vec![Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C))],
        adjacency: AdjacencyList::from_topology(1, &[]),
        ..TopologyBlock::default()
    };
    assert_eq!(
        assign_chiral_tags_from_structure(
            &invalid_topology,
            &CoordinateBlock::default(),
            &ValenceAssignment {
                explicit_valence: vec![0],
                implicit_hydrogens: vec![0],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            }
        ))
    );

    let topology = simple_star(6, 1);
    let coordinates = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]);
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 3],
                implicit_hydrogens: vec![0; 4],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::InvalidValence {
            field: "explicit_valence",
            actual: 3,
            atom_count: 4,
        })
    );
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0; 2],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::InvalidValence {
            field: "implicit_hydrogens",
            actual: 2,
            atom_count: 4,
        })
    );
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0, 0, -1, 0],
                implicit_hydrogens: vec![0; 4],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::InvalidValenceValue {
            field: "explicit_valence",
            atom: AtomId::new(2),
            value: -1,
        })
    );
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0, -2, 0, 0],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::InvalidValenceValue {
            field: "implicit_hydrogens",
            atom: AtomId::new(1),
            value: -2,
        })
    );

    let short_coordinates = simple_coordinates(vec![[0.0, 0.0, 0.0]; 3]);
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &short_coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0; 4],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::ConformerAtomCountMismatch {
            conformer: 0,
            rows: 3,
            atom_count: 4,
        })
    );

    let duplicate_coordinates = CoordinateBlock {
        conformers_3d: vec![
            Conformer3D::new(9, vec![[0.0, 0.0, 0.0]; 4], true),
            Conformer3D::new(9, vec![[1.0, 0.0, 0.0]; 4], true),
        ],
        ..CoordinateBlock::default()
    };
    assert_eq!(
        assign_chiral_tags_from_structure(
            &topology,
            &duplicate_coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0; 4],
            },
            &StructureTagParams::default(),
        ),
        Err(StereoError::DuplicateConformerId { id: 9 })
    );
}

#[test]
fn wedge_and_dash_explicitness_prevent_nonexplicit_property_insertion() {
    for direction in [BondDirection::BeginWedge, BondDirection::BeginDash] {
        let mut topology = simple_star(6, 1);
        topology.bonds[0].set_direction(direction);
        let before = topology.clone();
        let coordinates = simple_coordinates(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        let assignment = assign_chiral_tags_from_structure(
            &topology,
            &coordinates,
            &ValenceAssignment {
                explicit_valence: vec![0; 4],
                implicit_hydrogens: vec![0; 4],
            },
            &StructureTagParams::default(),
        )
        .unwrap();
        assert_eq!(
            assignment.topology.atoms[0].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(
            assignment.topology.atoms[0].prop("_NonExplicit3DChirality"),
            None
        );
        assert_eq!(topology, before);
    }
}

#[test]
fn numeric_extremes_run_in_a_process_with_source_default_environment() {
    let output = Command::new(std::env::current_exe().expect("current test executable"))
        .arg("--exact")
        .arg("structure_tags_numeric_extremes_child")
        .arg("--nocapture")
        .env(NUMERIC_CHILD_ENV, "1")
        .env_remove(NONTETRAHEDRAL_ENV)
        .output()
        .expect("run isolated numeric boundary test");
    assert!(
        output.status.success(),
        "numeric boundary child failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn structure_tags_numeric_extremes_child() {
    if std::env::var_os(NUMERIC_CHILD_ENV).is_none() {
        return;
    }
    let tetrahedral = simple_star(6, 1);
    let signed_zero = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -0.0],
    ]);
    let valence = ValenceAssignment {
        explicit_valence: vec![0; 4],
        implicit_hydrogens: vec![0; 4],
    };
    let assignment = assign_chiral_tags_from_structure(
        &tetrahedral,
        &signed_zero,
        &valence,
        &StructureTagParams::default(),
    )
    .unwrap();
    assert_eq!(
        assignment.topology.atoms[0].chiral_tag(),
        ChiralTag::Unspecified
    );

    let phosphorus = simple_star(15, 0);
    let minimum_subnormal = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [f64::from_bits(1), 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ]);
    assert_eq!(
        assign_chiral_tags_from_structure(
            &phosphorus,
            &minimum_subnormal,
            &valence,
            &StructureTagParams::default(),
        ),
        Err(StereoError::ZeroLengthVector {
            center: AtomId::new(0),
            neighbor: AtomId::new(1),
        })
    );

    let maximum_finite = simple_coordinates(vec![
        [0.0, 0.0, 0.0],
        [f64::MAX, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ]);
    let assignment = assign_chiral_tags_from_structure(
        &phosphorus,
        &maximum_finite,
        &valence,
        &StructureTagParams::default(),
    )
    .unwrap();
    assert_eq!(
        assignment.topology.atoms[0].chiral_tag(),
        ChiralTag::Unspecified
    );
}
