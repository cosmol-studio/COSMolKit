//! Import the existing, checksummed descriptor reference through the parity
//! framework's preparation stage. This does not run CK or regenerate RDKit data.
//! The original manifest (including its malformed revision spelling) remains
//! explicit provenance; it is not promoted to a freshly authenticated oracle.
use crate::molecular::Outcome;
use crate::registry::{self, Corpus, Input, Record, Task, Value, molecule_plan::Profile};
use crate::{Result, digest, read, root};
use serde_json::{Value as Json, json};

const INPUT: &str = "testdata/smiles/corpus/smiles_5000.smi";
const DIRECTORY: &str = "testdata/descriptors/expected/rdkit/smiles_5000";
const INPUT_SHA: &str = "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849";
const GOLDEN_SHA: &str = "9b3aa0b3c5ca537388503dc4d3d495a2d04bebe0d64bf09fda0225ce452274c8";
const MANIFEST_SHA: &str = "7da152d6f7a454d2dddecb7356e847fd61367690cdab6c91cae69896212452aa";

pub fn handles(task: &Task) -> bool {
    use crate::registry::{Operation, molecule_plan::TaskId::*};
    matches!(
        task.operation,
        Operation::Molecular(
            Chi0 | Chi1
                | HallKierAlpha
                | HallKierAlphaWithContributions
                | Kappa1
                | Kappa2
                | Kappa3
                | Phi
                | Mqns
                | Chi0V
                | Chi1V
                | Chi2V
                | Chi3V
                | Chi4V
                | ChiNV
                | Chi0N
                | Chi1N
                | Chi2N
                | Chi3N
                | Chi4N
                | ChiNN
        )
    )
}

pub fn provenance(task: &Task) -> Option<Json> {
    handles(task).then(|| json!({
        "kind": "existing_reference_import",
        "path": format!("{DIRECTORY}/molecular_descriptors.jsonl"),
        "sha256": GOLDEN_SHA,
        "manifest_sha256": MANIFEST_SHA,
        "input_sha256": INPUT_SHA,
        "recorded_source_revision": "351f8f378f8ad6bbd517980c38896e66bf907af8c",
        "revision_status": "original 41-character spelling retained; not a fresh source-identity claim",
        "reference_version": "2026.03.1",
        "records": 5000
    }))
}

fn verified(path: &std::path::Path, expected: &str) -> Result<Vec<u8>> {
    let bytes = read(path)?;
    if digest(&bytes) != expected {
        return Err(format!(
            "{}: original reference checksum mismatch",
            path.display()
        ));
    }
    Ok(bytes)
}

fn bits(value: &Json) -> Result<u64> {
    let text = value.as_str().ok_or("missing reference f64 bits")?;
    if text.len() != 16 {
        return Err("invalid reference f64 bit width".into());
    }
    let value = u64::from_str_radix(text, 16).map_err(|e| e.to_string())?;
    if !f64::from_bits(value).is_finite() {
        return Err("nonfinite descriptor reference".into());
    }
    Ok(value)
}

fn observation(profile: Profile, row: &Json) -> Result<Outcome> {
    use Profile::*;
    let scalar = &row["high_feasibility_descriptor_bits"];
    let key = match profile {
        Chi0 => "chi_0",
        Chi1 => "chi_1",
        HallKierAlpha => "hall_kier_alpha",
        Kappa1 => "kappa_1",
        Kappa2 => "kappa_2",
        Kappa3 => "kappa_3",
        Phi => "phi",
        Chi0V => "chi_0v",
        Chi1V => "chi_1v",
        Chi2V => "chi_2v",
        Chi3V => "chi_3v",
        Chi4V => "chi_4v",
        Chi0N => "chi_0n",
        Chi1N => "chi_1n",
        Chi2N => "chi_2n",
        Chi3N => "chi_3n",
        Chi4N => "chi_4n",
        ChiNV { order } | ChiNN { order } => {
            let field = if matches!(profile, ChiNV { .. }) {
                "chi_nv_orders_0_6"
            } else {
                "chi_nn_orders_0_6"
            };
            let values = scalar[field].as_array().ok_or("missing Chi order vector")?;
            if values.len() != 7 {
                return Err("Chi order vector must contain orders 0 through 6".into());
            }
            return Ok(Outcome::Float64Bits(bits(
                values.get(order as usize).ok_or("unregistered Chi order")?,
            )?));
        }
        HallKierAlphaWithContributions => {
            let contribution = &row["high_feasibility_contribution_bits"]["hall_kier_alpha"];
            let atoms = contribution["atom_contributions"]
                .as_array()
                .ok_or("missing Hall atom contributions")?;
            // Original generator: alpha_contributions = [0.0] * alpha_mol.GetNumAtoms()
            // Its separate CalcNumAtoms descriptor includes implicit H and
            // cannot validate this vector length. Exact ordered comparison
            // below retains every contribution, including explicit H rows.
            return Ok(Outcome::Float64ContributionsBits {
                value: bits(&contribution["value"])?,
                atom_contributions: atoms.iter().map(bits).collect::<Result<_>>()?,
            });
        }
        Mqns { .. } => {
            // The pinned source calcMQNs explicitly ignores force. Exercise
            // both public arguments against the same 42 reference counts.
            let values = row["high_feasibility_descriptors"]["mqns"]
                .as_array()
                .ok_or("missing MQNs")?;
            if values.len() != 42 {
                return Err("MQNs must contain 42 entries".into());
            }
            return Ok(Outcome::UnsignedVector(
                values
                    .iter()
                    .map(|value| {
                        value
                            .as_u64()
                            .and_then(|value| value.try_into().ok())
                            .ok_or("MQN entry is not u32".into())
                    })
                    .collect::<Result<_>>()?,
            ));
        }
        _ => return Err("profile is outside the descriptor reference import".into()),
    };
    Ok(Outcome::Float64Bits(bits(&scalar[key])?))
}

pub fn generate(task: &Task, cases: &Corpus) -> Result<Vec<Record>> {
    if !handles(task) {
        return Err("task is outside descriptor import scope".into());
    }
    let directory = root().join(DIRECTORY);
    let input = verified(&root().join(INPUT), INPUT_SHA)?;
    let original_cases = crate::molecular::read_corpus(&root().join(INPUT))?;
    if cases.molecules != original_cases || original_cases.len() != 5000 {
        return Err("existing descriptor reference requires its exact ordered 5000 input".into());
    }
    let manifest_bytes = verified(&directory.join("manifest.json"), MANIFEST_SHA)?;
    let manifest: Json = serde_json::from_slice(&manifest_bytes).map_err(|e| e.to_string())?;
    let golden = verified(&directory.join("molecular_descriptors.jsonl"), GOLDEN_SHA)?;
    if manifest["input"]["sha256"] != digest(&input)
        || manifest["outputs"][0]["sha256"] != digest(&golden)
        || manifest["outputs"][0]["records"] != 5000
        || manifest["reference_implementation"]["version"] != registry::RDKIT_VERSION
    {
        return Err("existing descriptor manifest identity mismatch".into());
    }
    let rows: Vec<Json> = golden
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| serde_json::from_slice(line).map_err(|e| e.to_string()))
        .collect::<Result<_>>()?;
    if rows.len() != original_cases.len() {
        return Err("existing descriptor reference row count mismatch".into());
    }
    for (case, row) in original_cases.iter().zip(&rows) {
        if row["smiles"] != case.smiles || row["rdkit_ok"] != true || !row["error"].is_null() {
            return Err(format!(
                "{}: invalid or reordered original descriptor reference",
                case.id
            ));
        }
    }
    // No chemistry execution occurs here; preparation validates every row.
    registry::expand(cases, task)
        .into_iter()
        .map(|input| {
            let Input::Molecular { case, profile } = &input else {
                return Err("expected molecular descriptor input".into());
            };
            let index = case
                .id
                .strip_prefix("line:")
                .and_then(|n| n.parse::<usize>().ok())
                .and_then(|n| n.checked_sub(1))
                .ok_or("invalid case identity")?;
            let row = rows.get(index).ok_or("case identity outside reference")?;
            let output = Value::Molecular(observation(*profile, row)?);
            task.validate_reference(&input, &input, &output)?;
            Ok(Record { input, output })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hall_reference_rows_do_not_use_total_atom_count_with_implicit_h() {
        // Ethanol has three stored atoms but CalcNumAtoms reports nine.
        // Reproduce the original generator's three-row contribution output.
        let row = json!({
            "high_feasibility_descriptors": {"num_atoms": 9},
            "high_feasibility_contribution_bits": {"hall_kier_alpha": {
                "value": "0000000000000000",
                "atom_contributions": ["0000000000000000", "0000000000000000", "0000000000000000"]
            }}
        });
        assert_eq!(
            observation(Profile::HallKierAlphaWithContributions, &row).unwrap(),
            Outcome::Float64ContributionsBits {
                value: 0,
                atom_contributions: vec![0, 0, 0]
            }
        );
    }
}
