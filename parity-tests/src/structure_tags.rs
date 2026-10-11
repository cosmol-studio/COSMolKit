//! Fixed stereo construction recipes. Process-local environment is isolated.
use crate::{Result, tasks, value::Value};
use cosmolkit::*;
use serde_json::Value as Json;
use std::{
    collections::BTreeMap,
    io::{Read, Write},
    process::{Command, Stdio},
};

fn property(key: &str, v: &Json) -> Result<PropertyValue> {
    Ok(match v {
        Json::Bool(v) => PropertyValue::Bool(*v),
        Json::String(v) => PropertyValue::String(v.as_str().into()),
        Json::Number(v) if key == "_chiralPermutation" => {
            PropertyValue::UInt(u32::try_from(v.as_u64().ok_or("invalid permutation")?)?)
        }
        Json::Number(v) if v.is_i64() => PropertyValue::Int(i32::try_from(v.as_i64().unwrap())?),
        Json::Number(v) => PropertyValue::Double(v.as_f64().ok_or("invalid numeric property")?),
        _ => return Err("unsupported fixed property input".into()),
    })
}
fn properties(p: &BTreeMap<PropertyText, PropertyValue>) -> Value {
    Value::List(
        p.iter()
            .map(|(key, v)| {
                Value::List(vec![
                    Value::Bytes(key.as_bytes().to_vec()),
                    match v {
                        PropertyValue::String(v) => Value::Text(
                            String::from_utf8(v.as_bytes().to_vec()).expect("fixture UTF-8"),
                        ),
                        PropertyValue::Int(v) => Value::Int(*v as i64),
                        PropertyValue::UInt(v) => Value::Int(*v as i64),
                        PropertyValue::Bool(v) => Value::Bool(*v),
                        PropertyValue::Double(v) => Value::float(*v),
                        PropertyValue::IntVector(v) => {
                            Value::List(v.iter().map(|v| Value::Int(*v as i64)).collect())
                        }
                        PropertyValue::StringVector(v) => Value::List(
                            v.iter()
                                .map(|v| Value::Bytes(v.as_bytes().to_vec()))
                                .collect(),
                        ),
                    },
                ])
            })
            .collect(),
    )
}
fn molecule(row: &Json) -> Result<Molecule> {
    let atoms = row["atoms"]
        .as_array()
        .ok_or("missing stereo atoms")?
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let mut spec = AtomSpec::new(
                Element::from_atomic_number(u8::try_from(
                    a["atomic_number"].as_u64().ok_or("missing atomic number")?,
                )?)
                .ok_or("unmodeled element")?,
            )
            .with_formal_charge(i8::try_from(a["formal_charge"].as_i64().unwrap_or(0))?)
            .with_explicit_hydrogens(u8::try_from(a["explicit_hydrogens"].as_u64().unwrap_or(0))?)
            .with_chiral_tag(
                ChiralTag::from_rdkit_name(a["chiral_tag"].as_str().unwrap_or("CHI_UNSPECIFIED"))
                    .ok_or("invalid chiral tag")?,
            );
            if let Some(v) = a["chiral_permutation"].as_u64() {
                spec = spec.with_chiral_permutation(u32::try_from(v)?);
            }
            if let Some(props) = a["props"].as_object() {
                for (k, v) in props {
                    spec = spec.with_prop(k, property(k, v)?)?;
                }
            }
            Ok(Atom::from_spec(AtomId::new(i), spec))
        })
        .collect::<Result<Vec<_>>>()?;
    let bonds = row["bonds"]
        .as_array()
        .ok_or("missing stereo bonds")?
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let mut spec = BondSpec::new(
                AtomId::new(b["begin"].as_u64().ok_or("missing bond begin")? as usize),
                AtomId::new(b["end"].as_u64().ok_or("missing bond end")? as usize),
                BondOrder::from_rdkit_name(b["type"].as_str().unwrap_or("SINGLE"))
                    .ok_or("invalid bond type")?,
            )
            .with_direction(
                BondDirection::from_rdkit_name(b["direction"].as_str().unwrap_or("NONE"))
                    .ok_or("invalid bond direction")?,
            );
            if let Some(v) = b["unknown_stereo"].as_i64() {
                spec = spec.with_prop("_UnknownStereo", PropertyValue::Int(i32::try_from(v)?))?;
            }
            if let Some(props) = b["props"].as_object() {
                for (k, v) in props {
                    spec = spec.with_prop(k, property(k, v)?)?;
                }
            }
            Ok(Bond::from_spec(BondId::new(i), spec))
        })
        .collect::<Result<Vec<_>>>()?;
    let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])?;
    let conformers = row["conformers"]
        .as_array()
        .ok_or("missing stereo conformers")?
        .iter()
        .map(|c| {
            let bits: Vec<[u64; 3]> = serde_json::from_value(c["coordinate_bits"].clone())?;
            Ok(Conformer3D::new(
                c["id"].as_u64().ok_or("missing conformer id")? as usize,
                bits.into_iter().map(|p| p.map(f64::from_bits)).collect(),
                c["is_3d"].as_bool().ok_or("missing dimension")?,
            ))
        })
        .collect::<Result<_>>()?;
    let mut props = MoleculeProperties::default();
    for (k, v) in row["molecule_props"]
        .as_object()
        .ok_or("missing stereo molecule properties")?
    {
        props = props.with_prop(k, property(k, v)?)?;
    }
    Ok(Molecule::from_parts(
        topology,
        CoordinateBlock {
            conformers_3d: conformers,
            ..Default::default()
        },
        props,
    )?
    .with_assigned_valence_with_params(&ValenceParams {
        strict: false,
        ..Default::default()
    })?)
}
fn snapshot(m: &Molecule) -> Result<Value> {
    let metadata = m.atom_metadata(false)?;
    let atoms = m
        .topology()
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            // The public model stores this RDKit property in its typed slot.
            let mut atom_props = a.props().clone();
            if let Some(permutation) = a.chiral_permutation() {
                atom_props.insert(
                    "_chiralPermutation".into(),
                    PropertyValue::UInt(permutation),
                );
            }
            Value::List(vec![
                Value::Int(a.atomic_number() as i64),
                Value::Int(a.formal_charge() as i64),
                Value::Int(a.explicit_hydrogens() as i64),
                Value::Int(metadata[i].implicit_hydrogens as i64),
                Value::Int(a.chiral_tag().rdkit_code()),
                properties(&atom_props),
            ])
        })
        .collect();
    let bonds = m
        .topology()
        .bonds
        .iter()
        .map(|b| {
            Value::List(vec![
                Value::Int(b.begin().index() as i64),
                Value::Int(b.end().index() as i64),
                Value::Int(b.order().rdkit_code()),
                Value::Int(b.direction().rdkit_code()),
                properties(b.props()),
            ])
        })
        .collect();
    Ok(Value::List(vec![
        Value::List(atoms),
        Value::List(bonds),
        crate::geometry::coordinates(m),
        properties(m.properties().props()),
    ]))
}
pub fn child() -> Result<()> {
    let mut input = String::new();
    std::io::stdin()
        .take(1024 * 1024)
        .read_to_string(&mut input)?;
    let row: Json = serde_json::from_str(&input)?;
    let m = molecule(&row)?;
    let before = snapshot(&m)?;
    let params = StructureTagParams {
        conformer_id: i32::try_from(
            row["conf_id"]
                .as_i64()
                .ok_or("missing conformer selection")?,
        )?,
        replace_existing_tags: row["replace_existing_tags"]
            .as_bool()
            .ok_or("missing replacement policy")?,
    };
    let (status, after) = match m.with_chiral_tags_from_structure_with_params(&params) {
        Ok(after) => (Value::Bool(true), snapshot(&after)?),
        Err(OperationError::Stereo(StereoError::ConformerNotFound { requested })) => (
            Value::List(vec![
                Value::Text("conformer_not_found".into()),
                Value::Int(requested as i64),
            ]),
            Value::Null,
        ),
        Err(OperationError::Stereo(StereoError::ZeroLengthVector { .. })) => (
            Value::List(vec![Value::Text("zero_length_vector".into())]),
            Value::Null,
        ),
        Err(e) => return Err(e.into()),
    };
    if before != snapshot(&m)? {
        return Err("stereo assignment modified source state".into());
    }
    std::io::stdout().write_all(&tasks::success(Value::List(vec![status, before, after])))?;
    Ok(())
}
pub fn observe(row: &Json) -> Result<Vec<u8>> {
    let mut command = Command::new(std::env::current_exe()?);
    command
        .arg("structure-case")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    match row["environment"]["mode"].as_str() {
        Some("unset") => {
            command.env_remove("RDK_ENABLE_NONTETRAHEDRAL_STEREO");
        }
        Some("set") => {
            command.env(
                "RDK_ENABLE_NONTETRAHEDRAL_STEREO",
                row["environment"]["value"]
                    .as_str()
                    .ok_or("missing environment value")?,
            );
        }
        _ => return Err("invalid stereo environment mode".into()),
    };
    let mut child = command.spawn()?;
    serde_json::to_writer(child.stdin.take().ok_or("missing child stdin")?, row)?;
    let output = child.wait_with_output()?;
    if !output.status.success() {
        return Err(format!(
            "stereo case process {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }
    if output.stdout.is_empty() {
        return Err("stereo case produced no outcome".into());
    }
    Ok(output.stdout)
}
