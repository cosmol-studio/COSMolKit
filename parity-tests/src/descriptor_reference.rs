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
            NumAmideBonds
                | NumSpiroAtoms
                | NumBridgeheadAtoms
                | NumAtomStereoCenters
                | NumUnspecifiedAtomStereoCenters
                | NumRotatableBonds
                | CrippenDescriptors
                | LabuteAsa
                | LabuteAsaContributions
                | Tpsa
                | SlogpVsa
                | SmrVsa
                | SlogpVsa1
                | SlogpVsa2
                | SlogpVsa3
                | SlogpVsa4
                | SlogpVsa5
                | SlogpVsa6
                | SlogpVsa7
                | SlogpVsa8
                | SlogpVsa9
                | SlogpVsa10
                | SlogpVsa11
                | SlogpVsa12
                | SmrVsa1
                | SmrVsa2
                | SmrVsa3
                | SmrVsa4
                | SmrVsa5
                | SmrVsa6
                | SmrVsa7
                | SmrVsa8
                | SmrVsa9
                | SmrVsa10
                | Qed
                | Chi0
                | Chi1
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

fn bit_vector(value: &Json) -> Result<Vec<u64>> {
    value
        .as_array()
        .ok_or("missing reference f64 vector")?
        .iter()
        .map(bits)
        .collect()
}

fn observation(profile: Profile, row: &Json) -> Result<Outcome> {
    use Profile::*;
    let scalar = &row["high_feasibility_descriptor_bits"];
    let key = match profile {
        Chi0VWithParams { .. } => "chi_0v",
        Chi1VWithParams { .. } => "chi_1v",
        Chi2VWithParams { .. } => "chi_2v",
        Chi3VWithParams { .. } => "chi_3v",
        Chi4VWithParams { .. } => "chi_4v",
        Chi0NWithParams { .. } => "chi_0n",
        Chi1NWithParams { .. } => "chi_1n",
        Chi2NWithParams { .. } => "chi_2n",
        Chi3NWithParams { .. } => "chi_3n",
        Chi4NWithParams { .. } => "chi_4n",
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
        ChiNV { order }
        | ChiNN { order }
        | ChiNVWithParams { order, .. }
        | ChiNNWithParams { order, .. } => {
            let field = if matches!(profile, ChiNV { .. } | ChiNVWithParams { .. }) {
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
        NumAmideBonds => {
            return Ok(Outcome::Unsigned(
                row["high_feasibility_descriptors"]["num_amide_bonds"]
                    .as_u64()
                    .and_then(|v| v.try_into().ok())
                    .ok_or("missing u32 num_amide_bonds")?,
            ));
        }
        NumSpiroAtoms => {
            return Ok(Outcome::Unsigned(
                row["high_feasibility_descriptors"]["num_spiro_atoms"]
                    .as_u64()
                    .and_then(|v| v.try_into().ok())
                    .ok_or("missing u32 num_spiro_atoms")?,
            ));
        }
        NumBridgeheadAtoms => {
            return Ok(Outcome::Unsigned(
                row["high_feasibility_descriptors"]["num_bridgehead_atoms"]
                    .as_u64()
                    .and_then(|v| v.try_into().ok())
                    .ok_or("missing u32 num_bridgehead_atoms")?,
            ));
        }
        NumAtomStereoCenters => {
            return Ok(Outcome::Unsigned(
                row["high_feasibility_descriptors"]["num_atom_stereo_centers"]
                    .as_u64()
                    .and_then(|v| v.try_into().ok())
                    .ok_or("missing u32 num_atom_stereo_centers")?,
            ));
        }
        NumUnspecifiedAtomStereoCenters => {
            return Ok(Outcome::Unsigned(
                row["high_feasibility_descriptors"]["num_unspecified_atom_stereo_centers"]
                    .as_u64()
                    .and_then(|v| v.try_into().ok())
                    .ok_or("missing u32 num_unspecified_atom_stereo_centers")?,
            ));
        }
        SlogpVsa1 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][0])?)),
        SlogpVsa2 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][1])?)),
        SlogpVsa3 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][2])?)),
        SlogpVsa4 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][3])?)),
        SlogpVsa5 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][4])?)),
        SlogpVsa6 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][5])?)),
        SlogpVsa7 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][6])?)),
        SlogpVsa8 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][7])?)),
        SlogpVsa9 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][8])?)),
        SlogpVsa10 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][9])?)),
        SlogpVsa11 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][10])?)),
        SlogpVsa12 => return Ok(Outcome::Float64Bits(bits(&scalar["slogp_vsa"][11])?)),
        SmrVsa1 => return Ok(Outcome::Float64Bits(bits(&scalar["smr_vsa"][0])?)),
        SmrVsa2 => return Ok(Outcome::Float64Bits(bits(&scalar["smr_vsa"][1])?)),
        SmrVsa3 => return Ok(Outcome::Float64Bits(bits(&scalar["smr_vsa"][2])?)),
        SmrVsa4 => return Ok(Outcome::Float64Bits(bits(&scalar["smr_vsa"][3])?)),
        SmrVsa5 => return Ok(Outcome::Float64Bits(bits(&scalar["smr_vsa"][4])?)),
        SmrVsa6 => return Ok(Outcome::Float64Bits(bits(&scalar["smr_vsa"][5])?)),
        SmrVsa7 => return Ok(Outcome::Float64Bits(bits(&scalar["smr_vsa"][6])?)),
        SmrVsa8 => return Ok(Outcome::Float64Bits(bits(&scalar["smr_vsa"][7])?)),
        SmrVsa9 => return Ok(Outcome::Float64Bits(bits(&scalar["smr_vsa"][8])?)),
        SmrVsa10 => return Ok(Outcome::Float64Bits(bits(&scalar["smr_vsa"][9])?)),
        Qed => return Ok(Outcome::Float64Bits(bits(&row["descriptor_bits"]["qed"])?)),
        NumRotatableBonds { mode } => {
            use crate::registry::molecule_plan::RotatableBondMode as M;
            let mode = match mode {
                None | Some(M::Default) => "default",
                Some(M::NonStrict) => "non_strict",
                Some(M::Strict) => "strict",
                Some(M::StrictLinkages) => "strict_linkages",
            };
            return Ok(Outcome::Unsigned(
                row["descriptors"][format!("num_rotatable_bonds_{mode}")]
                    .as_u64()
                    .and_then(|v| v.try_into().ok())
                    .ok_or("missing rotatable count")?,
            ));
        }
        CrippenDescriptors {
            include_hydrogens,
            force,
        } => {
            let value = &row["descriptor_option_bits"]["crippen"][format!(
                "include_hs_{}_force_{}",
                include_hydrogens.unwrap_or(true),
                force
            )];
            return Ok(Outcome::Float64PairBits {
                first: bits(&value["logp"])?,
                second: bits(&value["molar_refractivity"])?,
            });
        }
        Tpsa {
            include_sulfur_phosphorus,
            force,
        } => {
            return Ok(Outcome::Float64Bits(bits(
                &row["descriptor_option_bits"]["tpsa"][format!(
                    "force_{force}_include_sandp_{}",
                    include_sulfur_phosphorus.unwrap_or(false)
                )],
            )?));
        }
        LabuteAsa {
            include_hydrogens, ..
        } => {
            return Ok(Outcome::Float64Bits(bits(
                &scalar[format!(
                    "labute_asa_include_hs_{}",
                    include_hydrogens.unwrap_or(true)
                )],
            )?));
        }
        LabuteAsaContributions {
            include_hydrogens, ..
        } => {
            let value = &row["high_feasibility_contribution_bits"]["labute_asa"]
                [format!("include_hs_{}", include_hydrogens.unwrap_or(true))];
            let atoms = bit_vector(&value["atom_contributions"])?;
            let hall =
                row["high_feasibility_contribution_bits"]["hall_kier_alpha"]["atom_contributions"]
                    .as_array()
                    .ok_or("missing stored-atom row identity")?;
            if atoms.len() != hall.len() {
                return Err("Labute contribution count differs from stored atoms".into());
            }
            return Ok(Outcome::LabuteAsaContributionsBits {
                asa: bits(&value["asa"])?,
                atom_contributions: atoms,
                hydrogen_contribution: bits(&value["hydrogen_contribution"])?,
            });
        }
        SlogpVsa { bins, .. } | SmrVsa { bins, .. } => {
            let family = if matches!(profile, SlogpVsa { .. }) {
                "slogp_vsa"
            } else {
                "smr_vsa"
            };
            let value = match bins {
                crate::registry::molecule_plan::VsaBins::Default => &scalar[family],
                crate::registry::molecule_plan::VsaBins::CustomDuplicates => {
                    &row["high_feasibility_cache_profile_bits"][format!("{family}_custom_forced")]
                }
            };
            return Ok(Outcome::Float64VectorBits(bit_vector(value)?));
        }
        LabuteAsaCacheSequence => {
            let rows = row["high_feasibility_cache_profile_bits"]["labute_asa_sequence"]
                .as_array()
                .ok_or("missing Labute sequence")?;
            let flags = [(false, false), (true, false), (true, true), (false, false)];
            if rows.len() != flags.len() {
                return Err("wrong Labute sequence length".into());
            }
            for (r, (include, force)) in rows.iter().zip(flags) {
                if r["include_hs"] != include || r["force"] != force {
                    return Err("wrong Labute sequence options".into());
                }
            }
            return Ok(Outcome::Float64VectorBits(
                rows.iter()
                    .map(|r| bits(&r["value"]))
                    .collect::<Result<_>>()?,
            ));
        }
        ChiNVCacheSequence | ChiNNCacheSequence => {
            let family = if matches!(profile, ChiNVCacheSequence) {
                "chi_nv"
            } else {
                "chi_nn"
            };
            let value = &row["high_feasibility_cache_profile_bits"][family];
            return Ok(Outcome::Float64VectorsBits(
                ["cold", "warm", "forced"]
                    .into_iter()
                    .map(|key| bit_vector(&value[key]))
                    .collect::<Result<_>>()?,
            ));
        }
        SlogpVsaCacheSequence | SmrVsaCacheSequence => {
            let keys: &[&str] = if matches!(profile, SlogpVsaCacheSequence) {
                &[
                    "slogp_vsa_default_cold",
                    "slogp_vsa_default_warm",
                    "slogp_vsa_custom_forced",
                ]
            } else {
                &["smr_vsa_default_warm", "smr_vsa_custom_forced"]
            };
            return Ok(Outcome::Float64VectorsBits(
                keys.iter()
                    .map(|key| bit_vector(&row["high_feasibility_cache_profile_bits"][key]))
                    .collect::<Result<_>>()?,
            ));
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

#[cfg(test)]
mod complete_descriptor_tests {
    use super::*;
    use crate::registry::molecule_plan::{Category, TaskId, VsaBins};

    #[test]
    fn complete_descriptor_registration_has_all78_and_preserves_individual_tasks() {
        let descriptors = registry::TASKS
            .iter()
            .filter(|task| {
                matches!(task.operation,
            registry::Operation::Molecular(id) if id.category()==Category::Descriptors)
            })
            .collect::<Vec<_>>();
        assert_eq!(descriptors.len(), 78);
        let selected = registry::select(Some("descriptors")).unwrap();
        assert_eq!(selected.len(), 78);
        assert!(
            selected
                .iter()
                .zip(&descriptors)
                .all(|(a, b)| std::ptr::eq(*a, *b))
        );
        assert_eq!(descriptors.iter().filter(|task| handles(task)).count(), 56);
        assert_eq!(TaskId::ChiNV.profiles().len(), 22);
        assert_eq!(TaskId::SlogpVsa.profiles().len(), 6);
        assert_eq!(TaskId::LabuteAsa.profiles().len(), 6);
    }

    #[test]
    fn full_labute_reference_requires_every_atom_and_hydrogen_field() {
        let mut row = json!({"high_feasibility_contribution_bits": {
            "hall_kier_alpha": {"atom_contributions":["0000000000000000","0000000000000000"]},
            "labute_asa": {"include_hs_false": {"asa":"4000000000000000", "atom_contributions":["3ff0000000000000","3ff0000000000000"], "hydrogen_contribution":"0000000000000000"}}}});
        let profile = Profile::LabuteAsaContributions {
            include_hydrogens: Some(false),
            force: true,
        };
        let expected = Outcome::LabuteAsaContributionsBits {
            asa: 2.0_f64.to_bits(),
            atom_contributions: vec![1.0_f64.to_bits(); 2],
            hydrogen_contribution: 0.0_f64.to_bits(),
        };
        assert_eq!(observation(profile, &row).unwrap(), expected);
        row["high_feasibility_contribution_bits"]["labute_asa"]["include_hs_false"]["atom_contributions"] =
            json!(["4000000000000000"]);
        assert!(observation(profile, &row).is_err());
        row["high_feasibility_contribution_bits"]["labute_asa"]["include_hs_false"]["hydrogen_contribution"] =
            Json::Null;
        assert!(observation(profile, &row).is_err());
    }

    #[test]
    fn vsa_reference_keeps_duplicate_bin_vectors_and_rejects_wrong_shape() {
        let profile = Profile::SlogpVsa {
            bins: VsaBins::CustomDuplicates,
            force: Some(true),
        };
        let row = json!({"high_feasibility_cache_profile_bits":{"slogp_vsa_custom_forced":["0000000000000000","0000000000000000","3ff0000000000000","0000000000000000","0000000000000000","0000000000000000"]}});
        let outcome = observation(profile, &row).unwrap();
        crate::molecular::validate_output(&profile, &outcome).unwrap();
        assert!(
            crate::molecular::validate_output(&profile, &Outcome::Float64VectorBits(vec![0; 5]))
                .is_err()
        );
        assert!(
            crate::molecular::validate_output(
                &profile,
                &Outcome::Float64VectorBits(vec![f64::NAN.to_bits(); 6])
            )
            .is_err()
        );
    }

    #[test]
    fn labute_sequence_preserves_option_order_before_any_execution() {
        let mut row = json!({"high_feasibility_cache_profile_bits":{"labute_asa_sequence":[
            {"include_hs":false,"force":false,"value":"0000000000000000"},
            {"include_hs":true,"force":false,"value":"0000000000000000"},
            {"include_hs":true,"force":true,"value":"0000000000000000"},
            {"include_hs":false,"force":false,"value":"0000000000000000"}]}});
        assert_eq!(
            observation(Profile::LabuteAsaCacheSequence, &row).unwrap(),
            Outcome::Float64VectorBits(vec![0; 4])
        );
        row["high_feasibility_cache_profile_bits"]["labute_asa_sequence"][2]["force"] =
            json!(false);
        assert!(observation(Profile::LabuteAsaCacheSequence, &row).is_err());
    }
}
