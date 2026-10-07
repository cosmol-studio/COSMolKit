//! Ordinary SMILES writer corpus parity, independent of batch tests.
use crate::registry::{Input, Record, SmilesCase, Value};
use cosmolkit::{AtomId, Molecule, SmilesWriteParams};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Root {
    None,
    First,
    Last,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Profile {
    pub do_isomeric_smiles: bool,
    pub do_kekule: bool,
    pub canonical: bool,
    pub clean_stereo: bool,
    pub all_bonds_explicit: bool,
    pub all_hs_explicit: bool,
    pub include_dative_bonds: bool,
    pub ignore_atom_map_numbers: bool,
    pub rooted_at_atom: Root,
}

/// Same ordered product as `_generate_smiles_writer_golden.py::iter_branches`:
/// eight booleans, then none/first/last. Random traversal is fixed off.
pub fn profiles() -> Vec<Profile> {
    (0..256)
        .flat_map(|bits| {
            [Root::None, Root::First, Root::Last].map(move |rooted_at_atom| Profile {
                do_isomeric_smiles: bits & 128 != 0,
                do_kekule: bits & 64 != 0,
                canonical: bits & 32 != 0,
                clean_stereo: bits & 16 != 0,
                all_bonds_explicit: bits & 8 != 0,
                all_hs_explicit: bits & 4 != 0,
                include_dative_bonds: bits & 2 != 0,
                ignore_atom_map_numbers: bits & 1 != 0,
                rooted_at_atom,
            })
        })
        .collect()
}

impl Profile {
    fn params(self, atoms: usize) -> SmilesWriteParams {
        let rooted_at_atom = match self.rooted_at_atom {
            Root::None => None,
            Root::First => (atoms != 0).then_some(AtomId::new(0)),
            Root::Last => atoms.checked_sub(1).map(AtomId::new),
        };
        SmilesWriteParams {
            do_isomeric_smiles: self.do_isomeric_smiles,
            do_kekule: self.do_kekule,
            canonical: self.canonical,
            clean_stereo: self.clean_stereo,
            all_bonds_explicit: self.all_bonds_explicit,
            all_hydrogens_explicit: self.all_hs_explicit,
            include_dative_bonds: self.include_dative_bonds,
            ignore_atom_map_numbers: self.ignore_atom_map_numbers,
            rooted_at_atom,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WriteInput {
    pub case: SmilesCase,
    pub profile: Profile,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Outcome {
    Smiles(String),
    Error(Stage),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Stage {
    Parse,
    Write,
}

pub fn run(input: &WriteInput) -> Result<Record, String> {
    let output = match Molecule::from_smiles(&input.case.smiles) {
        Err(_) => Outcome::Error(Stage::Parse),
        Ok(mol) => match mol.to_smiles_with_params(&input.profile.params(mol.num_atoms())) {
            Ok(text) => Outcome::Smiles(text),
            Err(_) => Outcome::Error(Stage::Write),
        },
    };
    Ok(Record {
        input: Input::SmilesWrite(input.clone()),
        output: Value::SmilesWrite(output),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_has_all_768_unique_ordered_profiles_and_exact_parameter_mapping() {
        let profiles = profiles();
        assert_eq!(profiles.len(), 768);
        for (i, p) in profiles.iter().enumerate() {
            assert!(!profiles[..i].contains(p));
            let values = [
                p.do_isomeric_smiles,
                p.do_kekule,
                p.canonical,
                p.clean_stereo,
                p.all_bonds_explicit,
                p.all_hs_explicit,
                p.include_dative_bonds,
                p.ignore_atom_map_numbers,
            ];
            for (axis, value) in values.into_iter().enumerate() {
                assert_eq!(value, (i / 3) & (1 << (7 - axis)) != 0);
            }
            assert_eq!(
                p.rooted_at_atom,
                [Root::None, Root::First, Root::Last][i % 3]
            );
            let params = p.params(5);
            assert_eq!(params.do_isomeric_smiles, p.do_isomeric_smiles);
            assert_eq!(params.do_kekule, p.do_kekule);
            assert_eq!(params.canonical, p.canonical);
            assert_eq!(params.clean_stereo, p.clean_stereo);
            assert_eq!(params.all_bonds_explicit, p.all_bonds_explicit);
            assert_eq!(params.all_hydrogens_explicit, p.all_hs_explicit);
            assert_eq!(params.include_dative_bonds, p.include_dative_bonds);
            assert_eq!(params.ignore_atom_map_numbers, p.ignore_atom_map_numbers);
            assert_eq!(
                params.rooted_at_atom.map(AtomId::index),
                [None, Some(0), Some(4)][i % 3]
            );
            assert_eq!(p.params(0).rooted_at_atom, None);
            assert_eq!(
                p.params(1).rooted_at_atom.map(AtomId::index),
                [None, Some(0), Some(0)][i % 3]
            );
        }
    }
}
