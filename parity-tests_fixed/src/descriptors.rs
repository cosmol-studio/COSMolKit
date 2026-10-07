//! Typed projection of freshly generated RDKit descriptor rows.
use crate::{
    Result,
    molecular::Outcome,
    registry::{Task, molecule_plan::Profile},
};
use serde_json::Value as Json;
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

pub(crate) fn observation(profile: Profile, row: &Json) -> Result<Outcome> {
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
