//! Deterministic optional parameter sampling, never chemistry randomness.
use crate::{Result, tasks};
use serde_json::{Value, json};

pub fn sample(task: &str, count: usize, seed: u64) -> Result<Vec<Value>> {
    if count == 0 || crate::special::input(task).is_some() {
        return Ok(Vec::new());
    }
    let mut state = seed;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state >> 32
    };
    let base = tasks::profiles(task)?;
    let mut sampled = Vec::new();
    for _ in 0..count {
        let mut p = base[(next() as usize) % base.len()].clone();
        match task {
            task if crate::special::input(task).is_some() => {}
            "bio_pdb_output" => {}
            "fingerprint_morgan" => {
                p["radius"] = json!(next() % 5);
                p["fp_size"] = json!([257, 1024, 2048, 4096][next() as usize % 4]);
                for flag in [
                    "include_chirality",
                    "use_bond_types",
                    "include_ring_membership",
                    "only_nonzero_invariants",
                    "include_redundant_environments",
                    "count_simulation",
                    "feature_invariants",
                    "from_first",
                    "ignore_first",
                    "custom_atoms",
                    "custom_bonds",
                ] {
                    p[flag] = json!(next() % 2 != 0);
                }
                p["bits_per_feature"] = json!(1 + next() % 3);
            }
            "fingerprint_atom_pair" => {
                p["min_distance"] = json!(1 + next() % 3);
                p["max_distance"] = json!(4 + next() % 27);
                p["fp_size"] = json!([257, 1024, 2048][next() as usize % 3]);
                for flag in [
                    "include_chirality",
                    "count_simulation",
                    "from_first",
                    "ignore_first",
                    "custom_atoms",
                    "additional_output",
                ] {
                    p[flag] = json!(next() % 2 != 0);
                }
            }
            "fingerprint_torsion" => {
                p["torsion_atom_count"] = json!(3 + next() % 3);
                p["fp_size"] = json!([257, 1024, 2048][next() as usize % 3]);
                for flag in [
                    "include_chirality",
                    "count_simulation",
                    "only_shortest_paths",
                    "custom_atoms",
                    "additional_output",
                ] {
                    p[flag] = json!(next() % 2 != 0);
                }
            }
            "fingerprint_topological" => {
                p["min_path"] = json!(1 + next() % 2);
                p["max_path"] = json!(3 + next() % 5);
                p["fp_size"] = json!([257, 1024, 2048][next() as usize % 3]);
                for flag in ["use_hs", "branched_paths", "use_bond_order", "custom_atoms"] {
                    p[flag] = json!(next() % 2 != 0);
                }
            }
            "fingerprint_layered" => {
                p["layers"] = json!(next() as u32);
                p["min_path"] = json!(1 + next() % 2);
                p["max_path"] = json!(3 + next() % 5);
                p["fp_size"] = json!([64, 257, 2048][next() as usize % 3]);
                p["branched_paths"] = json!(next() % 2 != 0);
            }
            "fingerprint_pattern" => {
                p["fp_size"] = json!([7, 127, 257, 2048][next() as usize % 4]);
                p["tautomeric"] = json!(next() % 2 != 0);
            }
            "fingerprint_avalon" => {
                p["fp_size"] = json!([256, 512, 1024][next() as usize % 3]);
                p["is_query"] = json!(next() % 2 != 0);
            }
            "fingerprint_fuzzy" => {
                p["radius"] = json!(next() % 5);
                p["fp_size"] = json!([257, 1024, 2048][next() as usize % 3]);
                p["wide"] = json!(next() % 2 != 0);
            }
            "smiles_write" => {
                for flag in [
                    "isomeric_smiles",
                    "kekule",
                    "canonical",
                    "clean_stereo",
                    "all_bonds_explicit",
                    "all_hydrogens_explicit",
                    "include_dative_bonds",
                    "ignore_atom_map_numbers",
                ] {
                    p[flag] = json!(next() % 2 != 0);
                }
            }
            // These have a finite enumerated parameter space; selecting a
            // historical profile is the complete available sampling operation.
            "cip" | "valence" | "fingerprint_maccs" => (),
            task if crate::descriptors::TASKS
                .iter()
                .any(|(name, _)| *name == task) =>
            {
                ()
            }
            task if crate::chemistry::TASKS
                .iter()
                .any(|(name, _)| *name == task) =>
            {
                ()
            }
            "substructure" => {
                p["compiled_default"] = json!(false);
                for flag in [
                    "use_chirality",
                    "uniquify",
                    "recursion_possible",
                    "use_enhanced_stereo",
                    "specified_stereo_query_matches_unspecified",
                    "aromatic_matches_conjugated",
                    "aromatic_matches_single_or_double",
                ] {
                    p[flag] = json!(next() % 2 != 0);
                }
                p["max_matches"] = json!([1, 3, 100, 1000][next() as usize % 4]);
            }
            "tautomer" => {
                for flag in [
                    "remove_sp3_stereo",
                    "remove_bond_stereo",
                    "remove_isotopic_hydrogens",
                    "reassign_stereo",
                ] {
                    p[flag] = json!(next() % 2 != 0);
                }
                p["max_tautomers"] = json!([1, 10, 1000][next() as usize % 3]);
                p["max_transforms"] = json!([1, 10, 1000][next() as usize % 3]);
            }
            "coordinates_2d" => {
                p["canonical"] = json!(next() % 2 != 0);
                p["samples"] = json!(4);
                p["flips"] = json!(2);
                p["seed"] = json!(next() % 10000);
                p["permute"] = json!(next() % 2 != 0);
                p["ring_templates"] = json!(next() % 2 != 0);
            }
            "svg" => {
                p["width"] = json!([200, 300, 400][next() as usize % 3]);
                p["height"] = json!([200, 300, 400][next() as usize % 3]);
            }
            "bounds_matrix" | "uff_available" | "mmff_available" => {
                p["add_hs"] = json!(next() % 2 != 0);
            }
            "embedding" => {
                p["seed"] = json!(next() % 1000000);
                p["random_coords"] = json!(next() % 2 != 0);
                p["count"] = json!(1 + next() % 3);
            }
            "alignment" => {
                p["reflect"] = json!(next() % 2 != 0);
                p["explicit_map"] = json!(next() % 2 != 0);
                p["symmetrize"] = json!(next() % 2 != 0);
                p["seed"] = json!(next() % 10000);
            }
            task if crate::geometry::TASKS.iter().any(|(name, _)| *name == task) => {
                p["ignore_interfragment"] = json!(next() % 2 != 0);
                p["seed"] = json!(next() % 10000);
                p["fixed_first"] = json!(task.ends_with("force_field") && next() % 2 != 0);
                p["threshold"] = json!([10.0, 100.0][next() as usize % 2]);
            }
            _ => return Err(format!("no parameter sampler for {task}").into()),
        }
        sampled.push(tasks::normalize_profile(task, p)?);
    }
    Ok(sampled)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn samples_reproduce_exactly_without_hidden_oracle_parameters() {
        for (task, _) in tasks::available() {
            let a = sample(task, 4, 20260306).unwrap();
            assert_eq!(a, sample(task, 4, 20260306).unwrap());
            assert_eq!(
                a.len(),
                if crate::special::input(task).is_some() {
                    0
                } else {
                    4
                }
            );
        }
        assert_ne!(
            sample("fingerprint_morgan", 4, 1).unwrap(),
            sample("fingerprint_morgan", 4, 2).unwrap()
        );
    }
}
