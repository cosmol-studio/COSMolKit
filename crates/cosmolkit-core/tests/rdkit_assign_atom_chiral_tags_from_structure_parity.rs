mod common;

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::{Mutex, OnceLock};

use common::parity_data;
use cosmolkit_core::{
    AtomId, AtomSpec, BondDirection, BondOrder, BondSpec, ChiralTag, Conformer3D, Element, Molecule, MoleculeBuilder,
    OperationError, SanitizeOps, StereoError, ValenceModel, assign_valence_with_options,
};
use serde::Deserialize;
use serde_json::{Map, Value, json};

const EXPECTED_RECORD_COUNT: usize = 77;
const NON_TETRAHEDRAL_ENV: &str = "RDK_ENABLE_NONTETRAHEDRAL_STEREO";

static ENVIRONMENT_LOCK: OnceLock<Mutex<()>> = OnceLock::new();

#[derive(Debug, Deserialize)]
struct Fixture {
    schema_version: u32,
    defaults: Value,
    cases: Vec<Value>,
    octahedral_switch_cases: Vec<Value>,
}

#[derive(Debug, Deserialize)]
struct OracleRecord {
    case_id: String,
    selection_reason: String,
    environment: Value,
    conf_id: i32,
    replace_existing_tags: bool,
    status: String,
    error_type: Option<String>,
    error_text: Option<String>,
    selected_conformer_id: Option<usize>,
    before: Value,
    after: Value,
}

struct EnvironmentGuard {
    original: Option<String>,
}

impl EnvironmentGuard {
    fn apply(config: &Value) -> Self {
        let original = std::env::var(NON_TETRAHEDRAL_ENV).ok();
        match string_field(config, "mode") {
            "unset" => {
                // SAFETY: this integration target has one test and serializes all environment
                // access through ENVIRONMENT_LOCK.
                unsafe { std::env::remove_var(NON_TETRAHEDRAL_ENV) };
            }
            "set" => {
                let value = string_field(config, "value");
                // SAFETY: this integration target has one test and serializes all environment
                // access through ENVIRONMENT_LOCK.
                unsafe { std::env::set_var(NON_TETRAHEDRAL_ENV, value) };
            }
            mode => panic!("unsupported environment mode {mode:?}"),
        }
        Self { original }
    }
}

impl Drop for EnvironmentGuard {
    fn drop(&mut self) {
        match self.original.as_deref() {
            Some(value) => {
                // SAFETY: EnvironmentGuard remains under ENVIRONMENT_LOCK until it is dropped.
                unsafe { std::env::set_var(NON_TETRAHEDRAL_ENV, value) };
            }
            None => {
                // SAFETY: EnvironmentGuard remains under ENVIRONMENT_LOCK until it is dropped.
                unsafe { std::env::remove_var(NON_TETRAHEDRAL_ENV) };
            }
        }
    }
}

fn fixture_path() -> std::path::PathBuf {
    parity_data::repo_root().join("testdata/stereo/fixtures/assign_atom_chiral_tags_from_structure_cases.json")
}

fn load_fixture() -> Fixture {
    let path = fixture_path();
    let bytes = std::fs::read(&path).unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()));
    serde_json::from_slice(&bytes).unwrap_or_else(|error| panic!("failed to parse {}: {error}", path.display()))
}

fn load_oracle() -> Vec<OracleRecord> {
    let path = parity_data::golden_path("assign_atom_chiral_tags_from_structure.jsonl");
    let file = File::open(&path).unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()));
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(line_index, line)| {
            let line = line
                .unwrap_or_else(|error| panic!("failed to read {} line {}: {error}", path.display(), line_index + 1));
            serde_json::from_str(&line)
                .unwrap_or_else(|error| panic!("failed to parse {} line {}: {error}", path.display(), line_index + 1))
        })
        .collect()
}

fn deep_merge(base: &Value, update: &Value) -> Value {
    match (base, update) {
        (Value::Object(base), Value::Object(update)) => {
            let mut result = base.clone();
            for (key, value) in update {
                let merged = result
                    .get(key)
                    .map_or_else(|| value.clone(), |base_value| deep_merge(base_value, value));
                result.insert(key.clone(), merged);
            }
            Value::Object(result)
        }
        (_, update) => update.clone(),
    }
}

fn octahedral_case(raw: &Value) -> Value {
    let axes = HashMap::from([
        ("px", json!([1.0, 0.0, 0.0])),
        ("nx", json!([-1.0, 0.0, 0.0])),
        ("py", json!([0.0, 1.0, 0.0])),
        ("ny", json!([0.0, -1.0, 0.0])),
        ("pz", json!([0.0, 0.0, 1.0])),
        ("nz", json!([0.0, 0.0, -1.0])),
    ]);
    let ligands = array_field(raw, "neighbor_order")
        .iter()
        .map(|name| axes[string_value(name)].clone())
        .collect::<Vec<_>>();
    json!({
        "case_id": string_field(raw, "case_id"),
        "selection_reason": format!(
            "OctahedralPermFrom3D {} with {} VOLTEST",
            string_field(raw, "branch"),
            string_field(raw, "volume_sign")
        ),
        "center": {"atomic_number": 15},
        "ligands": ligands,
        "conformers": [{"id": 0, "is_3d": true}],
    })
}

fn object(value: &Value) -> &Map<String, Value> {
    value.as_object().expect("fixture value must be an object")
}

fn array_field<'a>(value: &'a Value, field: &str) -> &'a [Value] {
    object(value)[field]
        .as_array()
        .unwrap_or_else(|| panic!("fixture field {field:?} must be an array"))
}

fn string_field<'a>(value: &'a Value, field: &str) -> &'a str {
    string_value(&object(value)[field])
}

fn string_value(value: &Value) -> &str {
    value.as_str().expect("fixture value must be a string")
}

fn i64_field(value: &Value, field: &str) -> i64 {
    object(value)[field]
        .as_i64()
        .unwrap_or_else(|| panic!("fixture field {field:?} must be an integer"))
}

fn bool_field(value: &Value, field: &str) -> bool {
    object(value)[field]
        .as_bool()
        .unwrap_or_else(|| panic!("fixture field {field:?} must be a boolean"))
}

fn property_string(value: &Value) -> String {
    match value {
        Value::String(value) => value.clone(),
        Value::Bool(value) => value.to_string(),
        Value::Number(value) => value.to_string(),
        other => panic!("unsupported fixture property value {other}"),
    }
}

fn decode_number(value: &Value) -> f64 {
    if let Some(value) = value.as_f64() {
        return value;
    }
    match string_value(value) {
        "negative_zero" => -0.0,
        "min_subnormal" => f64::from_bits(1),
        "max_finite" => f64::MAX,
        "nextafter_0.1_down" => f64::from_bits(0x3fb9_9999_9999_9999),
        "nextafter_0.1_up" => f64::from_bits(0x3fb9_9999_9999_999b),
        "nextafter_negative_0.1_down" => f64::from_bits(0xbfb9_9999_9999_999b),
        "nextafter_negative_0.1_up" => f64::from_bits(0xbfb9_9999_9999_9999),
        "nextafter_zero_tolerance_down" => f64::from_bits(0x3c9c_d2b2_97d8_89bb),
        "zero_tolerance" => f64::from_bits(0x3c9c_d2b2_97d8_89bc),
        "nextafter_zero_tolerance_up" => f64::from_bits(0x3c9c_d2b2_97d8_89bd),
        "cos_100_degrees" => f64::from_bits(0xbfc6_3a1a_7e0b_7388),
        "sin_100_degrees" => f64::from_bits(0x3fef_838b_8c81_1c17),
        "nan" => f64::NAN,
        "positive_infinity" => f64::INFINITY,
        "negative_infinity" => f64::NEG_INFINITY,
        token => panic!("unknown numeric token {token:?}"),
    }
}

fn decode_point(value: &Value) -> [f64; 3] {
    let values = value.as_array().expect("coordinate must be an array");
    assert_eq!(values.len(), 3, "coordinate must have exactly three values");
    [
        decode_number(&values[0]),
        decode_number(&values[1]),
        decode_number(&values[2]),
    ]
}

fn configure_atom(config: &Value) -> AtomSpec {
    let atomic_number = object(config).get("atomic_number").and_then(Value::as_u64).unwrap_or(0) as u8;
    let mut atom = AtomSpec::new(
        Element::from_atomic_number(atomic_number)
            .unwrap_or_else(|| panic!("unsupported atomic number {atomic_number}")),
    )
    .with_formal_charge(object(config).get("formal_charge").and_then(Value::as_i64).unwrap_or(0) as i8)
    .with_explicit_hydrogens(
        object(config)
            .get("explicit_hydrogens")
            .and_then(Value::as_u64)
            .unwrap_or(0) as u8,
    )
    .with_chiral_tag(
        ChiralTag::from_rdkit_name(
            object(config)
                .get("chiral_tag")
                .and_then(Value::as_str)
                .unwrap_or("CHI_UNSPECIFIED"),
        )
        .expect("fixture chiral tag must be supported"),
    );
    if let Some(permutation) = object(config).get("chiral_permutation").and_then(Value::as_u64) {
        atom = atom.with_chiral_permutation(permutation as u32);
    }
    if let Some(props) = object(config).get("props").and_then(Value::as_object) {
        for (name, value) in props {
            atom = atom.with_prop(name, property_string(value));
        }
    }
    atom
}

fn bond_order(name: &str) -> BondOrder {
    match name {
        "ZERO" => BondOrder::Zero,
        name => BondOrder::from_rdkit_name(name).unwrap_or_else(|| panic!("unsupported fixture bond type {name:?}")),
    }
}

fn bond_direction(name: &str) -> BondDirection {
    match name {
        "NONE" => BondDirection::None,
        "BEGINWEDGE" => BondDirection::BeginWedge,
        "BEGINDASH" => BondDirection::BeginDash,
        "ENDUPRIGHT" => BondDirection::EndUpRight,
        "ENDDOWNRIGHT" => BondDirection::EndDownRight,
        "EITHERDOUBLE" => BondDirection::EitherDouble,
        "UNKNOWN" => BondDirection::Unknown,
        other => panic!("unsupported fixture bond direction {other:?}"),
    }
}

fn add_bond(builder: &mut MoleculeBuilder, config: &Value) {
    let begin = AtomId::new(i64_field(config, "begin") as usize);
    let end = AtomId::new(i64_field(config, "end") as usize);
    let mut bond = BondSpec::new(
        begin,
        end,
        bond_order(object(config).get("type").and_then(Value::as_str).unwrap_or("SINGLE")),
    )
    .with_direction(bond_direction(
        object(config)
            .get("direction")
            .and_then(Value::as_str)
            .unwrap_or("NONE"),
    ));
    if let Some(unknown) = object(config).get("unknown_stereo").and_then(Value::as_i64) {
        bond = bond.with_prop("_UnknownStereo", unknown.to_string());
    }
    if let Some(props) = object(config).get("props").and_then(Value::as_object) {
        for (name, value) in props {
            bond = bond.with_prop(name, property_string(value));
        }
    }
    builder.add_bond(bond).expect("fixture bond must reference valid atoms");
}

fn build_molecule(case: &Value) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    match string_field(case, "builder") {
        "explicit" => {
            for atom in array_field(case, "atoms") {
                builder.add_atom(configure_atom(atom));
            }
            for bond in array_field(case, "bonds") {
                add_bond(&mut builder, bond);
            }
        }
        "single_center" => {
            builder.add_atom(configure_atom(&object(case)["center"]));
            let ligand_elements = array_field(case, "ligand_atomic_numbers");
            for (ligand_index, _) in array_field(case, "ligands").iter().enumerate() {
                let atomic_number = ligand_elements[ligand_index]
                    .as_u64()
                    .expect("ligand atomic number must be unsigned") as u8;
                builder.add_atom(AtomSpec::new(
                    Element::from_atomic_number(atomic_number).expect("ligand atomic number must be supported"),
                ));
            }

            let overrides = object(case)
                .get("bond_overrides")
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
                .map(|item| (i64_field(item, "ligand") as usize, item))
                .collect::<HashMap<_, _>>();
            for ligand_index in 0..array_field(case, "ligands").len() {
                let mut config = json!({
                    "begin": 0,
                    "end": ligand_index + 1,
                    "type": string_field(case, "bond_type"),
                    "direction": string_field(case, "bond_direction"),
                    "unknown_stereo": object(case)["unknown_stereo"].clone(),
                });
                if let Some(override_config) = overrides.get(&ligand_index) {
                    config = deep_merge(&config, override_config);
                }
                let config_object = config.as_object_mut().expect("bond config must be an object");
                config_object.remove("ligand");
                if config_object
                    .remove("reverse")
                    .and_then(|value| value.as_bool())
                    .unwrap_or(false)
                {
                    let begin = config_object["begin"].clone();
                    let end = config_object["end"].clone();
                    config_object.insert("begin".to_string(), end);
                    config_object.insert("end".to_string(), begin);
                }
                add_bond(&mut builder, &config);
            }
        }
        kind => panic!("unsupported fixture builder {kind:?}"),
    }

    if let Some(props) = object(case).get("molecule_props").and_then(Value::as_object) {
        for (name, value) in props {
            builder = builder.with_property(name, property_string(value));
        }
    }
    for conformer_config in array_field(case, "conformers") {
        let coordinates =
            if let Some(coordinates) = object(conformer_config).get("coordinates").and_then(Value::as_array) {
                coordinates.iter().map(decode_point).collect()
            } else {
                let center = object(conformer_config)
                    .get("center_position")
                    .unwrap_or(&object(case)["center_position"]);
                let ligands = object(conformer_config)
                    .get("ligands")
                    .and_then(Value::as_array)
                    .map(Vec::as_slice)
                    .unwrap_or_else(|| array_field(case, "ligands"));
                std::iter::once(decode_point(center))
                    .chain(ligands.iter().map(decode_point))
                    .collect()
            };
        let mut conformer = Conformer3D::new(
            i64_field(conformer_config, "id") as usize,
            coordinates,
            bool_field(conformer_config, "is_3d"),
        );
        if let Some(props) = object(conformer_config).get("props").and_then(Value::as_object) {
            for (name, value) in props {
                conformer = conformer.with_prop(name, property_string(value));
            }
        }
        builder
            .add_conformer(conformer)
            .expect("fixture conformer rows must match atom count");
    }

    builder
        .build()
        .expect("fixture molecule must satisfy structural invariants")
        .sanitize_with_ops(SanitizeOps::NONE)
        .expect("RDKit UpdatePropertyCache(strict=False) fixture setup must succeed")
}

fn typed_props(props: &BTreeMap<String, String>, typed_integer_keys: &[&str]) -> Value {
    let mut result = Map::new();
    for (name, value) in props {
        if typed_integer_keys.contains(&name.as_str()) {
            result.insert(
                name.clone(),
                json!(
                    value
                        .parse::<i64>()
                        .unwrap_or_else(|error| { panic!("typed property {name:?} must contain an integer: {error}") })
                ),
            );
        } else {
            result.insert(name.clone(), Value::String(value.clone()));
        }
    }
    Value::Object(result)
}

fn python_float_hex(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value == f64::INFINITY {
        return "inf".to_string();
    }
    if value == f64::NEG_INFINITY {
        return "-inf".to_string();
    }
    let bits = value.to_bits();
    let sign = if bits >> 63 == 0 { "" } else { "-" };
    let exponent_bits = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & 0x000f_ffff_ffff_ffff;
    if exponent_bits == 0 {
        if fraction == 0 {
            format!("{sign}0x0.0p+0")
        } else {
            format!("{sign}0x0.{fraction:013x}p-1022")
        }
    } else {
        let exponent = exponent_bits - 1023;
        format!("{sign}0x1.{fraction:013x}p{exponent:+}")
    }
}

fn snapshot_molecule(molecule: &Molecule) -> Value {
    let valence = assign_valence_with_options(molecule, ValenceModel::RdkitLike, false)
        .expect("snapshot valence must match fixture preparation");
    let atoms = molecule
        .atoms()
        .iter()
        .enumerate()
        .map(|(index, atom)| {
            let mut props = match typed_props(atom.props(), &["_NonExplicit3DChirality"]) {
                Value::Object(props) => props,
                _ => unreachable!(),
            };
            if let Some(permutation) = atom.chiral_permutation() {
                props.insert("_chiralPermutation".to_string(), json!(permutation));
            }
            json!({
                "index": index,
                "atomic_number": atom.atomic_number(),
                "formal_charge": atom.formal_charge(),
                "explicit_hydrogens": atom.explicit_hydrogens(),
                "implicit_hydrogens": valence.implicit_hydrogens[index],
                "chiral_tag": atom.chiral_tag().rdkit_name(),
                "chiral_permutation": atom.chiral_permutation(),
                "non_explicit_3d_chirality": atom.prop("_NonExplicit3DChirality")
                    .map(|value| value.parse::<i64>().expect("3D chirality property must be integer")),
                "props": props,
            })
        })
        .collect::<Vec<_>>();
    let bonds = molecule
        .bonds()
        .iter()
        .enumerate()
        .map(|(index, bond)| {
            json!({
                "index": index,
                "begin": bond.begin().index(),
                "end": bond.end().index(),
                "type": bond.order().rdkit_name(),
                "direction": bond.direction().rdkit_name(),
                "unknown_stereo": bond.prop("_UnknownStereo")
                    .map(|value| value.parse::<i64>().expect("unknown stereo property must be integer")),
                "props": typed_props(bond.props(), &["_UnknownStereo"]),
            })
        })
        .collect::<Vec<_>>();
    let conformers = molecule
        .conformers_3d()
        .iter()
        .map(|conformer| {
            let coordinates = conformer
                .coordinates()
                .iter()
                .map(|point| point.map(python_float_hex))
                .collect::<Vec<_>>();
            json!({
                "id": conformer.id(),
                "is_3d": conformer.is_3d(),
                "coordinates": coordinates,
                "props": typed_props(conformer.props(), &[]),
            })
        })
        .collect::<Vec<_>>();
    json!({
        "atom_count": molecule.num_atoms(),
        "bond_count": molecule.num_bonds(),
        "atoms": atoms,
        "bonds": bonds,
        "conformers": conformers,
        "molecule_props": typed_props(molecule.properties().props(), &["_StereochemDone"]),
        "stereochem_done": molecule.prop("_StereochemDone")
            .map(|value| value.parse::<i64>().expect("stereochem-done property must be integer")),
    })
}

fn selected_conformer_id(molecule: &Molecule, conformer_id: i32) -> Option<usize> {
    if conformer_id < 0 {
        molecule.conformers_3d().first().map(Conformer3D::id)
    } else {
        molecule
            .conformers_3d()
            .iter()
            .find(|conformer| conformer.id() == conformer_id as usize)
            .map(Conformer3D::id)
    }
}

fn rdkit_error(error: &OperationError) -> (&'static str, &'static str) {
    match error {
        OperationError::Stereo {
            source: StereoError::ConformerNotFound { .. },
            ..
        } => ("ValueError", "Bad Conformer Id"),
        OperationError::Stereo {
            source: StereoError::ZeroLengthVector,
            ..
        } => ("RuntimeError", "Cannot normalize a zero length vector"),
        other => panic!("unexpected structured Rust error: {other:?}"),
    }
}

#[test]
fn assign_atom_chiral_tags_from_structure_matches_pinned_rdkit_for_every_fixture_row() {
    let _environment_lock = ENVIRONMENT_LOCK
        .get_or_init(|| Mutex::new(()))
        .lock()
        .expect("environment lock must not be poisoned");
    let fixture = load_fixture();
    assert_eq!(fixture.schema_version, 1);
    let mut cases = fixture.cases;
    cases.extend(fixture.octahedral_switch_cases.iter().map(octahedral_case));
    let oracle = load_oracle();
    assert_eq!(cases.len(), EXPECTED_RECORD_COUNT);
    assert_eq!(oracle.len(), EXPECTED_RECORD_COUNT);

    for (row_index, (raw_case, expected)) in cases.iter().zip(&oracle).enumerate() {
        let case = deep_merge(&fixture.defaults, raw_case);
        assert_eq!(string_field(&case, "case_id"), expected.case_id);
        assert_eq!(string_field(&case, "selection_reason"), expected.selection_reason);
        assert_eq!(object(&case)["environment"], expected.environment);
        assert_eq!(i64_field(&case, "conf_id") as i32, expected.conf_id);
        assert_eq!(
            bool_field(&case, "replace_existing_tags"),
            expected.replace_existing_tags
        );

        let _environment = EnvironmentGuard::apply(&expected.environment);
        let molecule = build_molecule(&case);
        assert_eq!(
            snapshot_molecule(&molecule),
            expected.before,
            "{} row {} fixture translation differs from pinned RDKit before state",
            expected.case_id,
            row_index + 1
        );
        assert_eq!(
            selected_conformer_id(&molecule, expected.conf_id),
            expected.selected_conformer_id,
            "{} selected conformer differs",
            expected.case_id
        );

        match molecule.with_chiral_tags_from_structure(expected.conf_id, expected.replace_existing_tags) {
            Ok(actual) => {
                assert_eq!(expected.status, "ok", "{} status differs", expected.case_id);
                assert_eq!(expected.error_type, None, "{} error type differs", expected.case_id);
                assert_eq!(expected.error_text, None, "{} error text differs", expected.case_id);
                assert_eq!(
                    snapshot_molecule(&actual),
                    expected.after,
                    "{} row {} after state differs from pinned RDKit",
                    expected.case_id,
                    row_index + 1
                );
            }
            Err(error) => {
                assert_eq!(expected.status, "error", "{} status differs", expected.case_id);
                let (error_type, error_text) = rdkit_error(&error);
                assert_eq!(expected.error_type.as_deref(), Some(error_type));
                assert_eq!(expected.error_text.as_deref(), Some(error_text));
                assert_eq!(
                    snapshot_molecule(&molecule),
                    expected.before,
                    "{} failed value operation mutated its source",
                    expected.case_id
                );
            }
        }
    }
}
