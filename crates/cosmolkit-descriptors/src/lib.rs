//! Detached molecular mass and formula descriptors.
//!
//! The descriptor boundary accepts only a validated [`TopologyBlock`]. It
//! cannot observe or mutate a live runtime molecule or derived-cache state.

use std::{cmp::Ordering, collections::BTreeMap};

use cosmolkit_core::{
    ValenceModel, assign_valence_with_options_for_topology, atomic_mass, most_common_isotope_mass,
    rdkit_element_symbol,
};
use cosmolkit_model::TopologyBlock;

const RDKIT_ELECTRON_MASS: f64 = 0.00054857991;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum DescriptorError {
    #[error("descriptor `{function}` cannot evaluate this topology: {detail}")]
    Unsupported {
        function: &'static str,
        detail: String,
    },
}

pub type DescriptorResult<T> = Result<T, DescriptorError>;

fn validate_topology(topology: &TopologyBlock, function: &'static str) -> DescriptorResult<()> {
    topology
        .validate()
        .map_err(|error| DescriptorError::Unsupported {
            function,
            detail: error.to_string(),
        })
}

fn valence(
    topology: &TopologyBlock,
    function: &'static str,
) -> DescriptorResult<cosmolkit_core::ValenceAssignment> {
    assign_valence_with_options_for_topology(topology, ValenceModel::RdkitLike, false).map_err(
        |error| DescriptorError::Unsupported {
            function,
            detail: error.to_string(),
        },
    )
}

/// Calculates average molecular weight using all atoms and implicit Hs.
#[must_use]
pub fn molecular_weight(topology: &TopologyBlock) -> DescriptorResult<f64> {
    molecular_weight_with_options(topology, false)
}

/// RDKit-compatible average molecular weight with an explicit heavy-atom mode.
#[must_use]
pub fn molecular_weight_with_options(
    topology: &TopologyBlock,
    only_heavy: bool,
) -> DescriptorResult<f64> {
    // RDKit✔️✔️: double calcAMW(const ROMol &mol, bool onlyHeavy) {
    // RDKit✔️✔️:   return MolOps::getAvgMolWt(mol, onlyHeavy);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double getAvgMolWt(const ROMol &mol, bool onlyHeavy) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   for (const auto &atom : mol.atoms()) {
    // RDKit✔️✔️:     if (!onlyHeavy || atom->getAtomicNum() != 1) res += atom->getMass();
    // RDKit✔️✔️:     if (!onlyHeavy) res += atom->getTotalNumHs() * table->getAtomicWeight(1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    validate_topology(topology, "molecular_weight")?;
    let assignment = valence(topology, "molecular_weight")?;
    let hydrogen_mass = atomic_mass(1, None);
    let mut result = 0.0;
    for (index, atom) in topology.atoms.iter().enumerate() {
        if !only_heavy || atom.atomic_number() != 1 {
            result += atomic_mass(atom.atomic_number(), atom.isotope());
        }
        if !only_heavy {
            let implicit = assignment
                .implicit_hydrogens
                .get(index)
                .copied()
                .unwrap_or(0)
                .max(0);
            result +=
                f64::from(u32::from(atom.explicit_hydrogens()) + implicit as u32) * hydrogen_mass;
        }
    }
    Ok(result)
}

/// Calculates exact molecular weight using the most common isotope mass.
#[must_use]
pub fn exact_molecular_weight(topology: &TopologyBlock) -> DescriptorResult<f64> {
    exact_molecular_weight_with_options(topology, false)
}

/// RDKit-compatible exact molecular weight with an explicit heavy-atom mode.
#[must_use]
pub fn exact_molecular_weight_with_options(
    topology: &TopologyBlock,
    only_heavy: bool,
) -> DescriptorResult<f64> {
    // RDKit✔️✔️: double calcExactMW(const ROMol &mol, bool onlyHeavy) {
    // RDKit✔️✔️:   return MolOps::getExactMolWt(mol, onlyHeavy);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: const double electronMass = 0.00054857991;
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (atNum != 1 || !onlyHeavy) res += table->getMostCommonIsotopeMass(atNum);
    // RDKit✔️✔️:     res -= constants::electronMass * atom->getFormalCharge();
    // RDKit✔️✔️:     if (!onlyHeavy) nHsToCount += atom->getTotalNumHs(false);
    // RDKit✔️✔️:   }
    validate_topology(topology, "exact_molecular_weight")?;
    let assignment = valence(topology, "exact_molecular_weight")?;
    let mut result = 0.0;
    let mut hydrogens_to_count = 0_i32;
    for (index, atom) in topology.atoms.iter().enumerate() {
        if atom.atomic_number() != 1 || !only_heavy {
            result += if atom.isotope().is_none() {
                most_common_isotope_mass(atom.atomic_number()).map_err(|error| {
                    DescriptorError::Unsupported {
                        function: "exact_molecular_weight",
                        detail: error.to_string(),
                    }
                })?
            } else {
                atomic_mass(atom.atomic_number(), atom.isotope())
            };
            result -= RDKIT_ELECTRON_MASS * f64::from(atom.formal_charge());
        }
        if !only_heavy {
            hydrogens_to_count += i32::from(atom.explicit_hydrogens())
                + assignment
                    .implicit_hydrogens
                    .get(index)
                    .copied()
                    .unwrap_or(0)
                    .max(0);
        }
    }
    if !only_heavy {
        result += f64::from(hydrogens_to_count)
            * most_common_isotope_mass(1).map_err(|error| DescriptorError::Unsupported {
                function: "exact_molecular_weight",
                detail: error.to_string(),
            })?;
    }
    Ok(result)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct FormulaKey {
    isotope: u32,
    symbol: &'static str,
}

fn hill_compare(left: &FormulaKey, right: &FormulaKey) -> Ordering {
    if left.symbol == "C" {
        return if right.symbol == "C" {
            left.isotope.cmp(&right.isotope)
        } else {
            Ordering::Less
        };
    }
    if right.symbol == "C" {
        return Ordering::Greater;
    }
    if left.symbol == "H" {
        return if right.symbol == "H" {
            left.isotope.cmp(&right.isotope)
        } else {
            Ordering::Less
        };
    }
    if right.symbol == "H" {
        return Ordering::Greater;
    }
    if left.symbol == "D" {
        return Ordering::Less;
    }
    if right.symbol == "D" {
        return Ordering::Greater;
    }
    if left.symbol == "T" {
        return Ordering::Less;
    }
    if right.symbol == "T" {
        return Ordering::Greater;
    }
    left.cmp(right)
}

/// Calculates the Hill-ordered molecular formula.
#[must_use]
pub fn molecular_formula(topology: &TopologyBlock) -> DescriptorResult<String> {
    molecular_formula_with_options(topology, false, false)
}

/// RDKit-compatible formula generation with isotope formatting controls.
#[must_use]
pub fn molecular_formula_with_options(
    topology: &TopologyBlock,
    separate_isotopes: bool,
    abbreviate_h_isotopes: bool,
) -> DescriptorResult<String> {
    // RDKit✔️✔️: std::string getMolFormula(const ROMol &mol, bool separateIsotopes,
    // RDKit✔️✔️:                           bool abbreviateHIsotopes) {
    // RDKit✔️✔️:   std::map<std::pair<unsigned int, std::string>, unsigned int> counts;
    // RDKit✔️✔️:   unsigned int nHs = 0;
    validate_topology(topology, "molecular_formula")?;
    let assignment = valence(topology, "molecular_formula")?;
    let mut counts = BTreeMap::<FormulaKey, u32>::new();
    let mut charge = 0_i32;
    let mut hydrogens = 0_u32;
    for (index, atom) in topology.atoms.iter().enumerate() {
        let atomic_number = atom.atomic_number();
        let mut key = FormulaKey {
            isotope: 0,
            symbol: rdkit_element_symbol(atomic_number).map_err(|error| {
                DescriptorError::Unsupported {
                    function: "molecular_formula",
                    detail: error.to_string(),
                }
            })?,
        };
        if separate_isotopes {
            let isotope = atom.isotope().map(u32::from).unwrap_or(0);
            if abbreviate_h_isotopes && atomic_number == 1 && (isotope == 2 || isotope == 3) {
                key.symbol = if isotope == 2 { "D" } else { "T" };
            } else {
                key.isotope = isotope;
            }
        }
        *counts.entry(key).or_insert(0) += 1;
        let total = i32::from(atom.explicit_hydrogens())
            + assignment
                .implicit_hydrogens
                .get(index)
                .copied()
                .unwrap_or(0);
        hydrogens += u32::try_from(total).map_err(|_| DescriptorError::Unsupported {
            function: "molecular_formula",
            detail: "negative total hydrogen count".to_owned(),
        })?;
        charge += i32::from(atom.formal_charge());
    }
    if hydrogens != 0 {
        *counts
            .entry(FormulaKey {
                isotope: 0,
                symbol: "H",
            })
            .or_insert(0) += hydrogens;
    }
    let mut keys = counts.keys().copied().collect::<Vec<_>>();
    keys.sort_by(hill_compare);
    let mut result = String::new();
    for key in keys {
        if key.isotope > 0 {
            result.push('[');
            result.push_str(&key.isotope.to_string());
            result.push_str(key.symbol);
            result.push(']');
        } else {
            result.push_str(key.symbol);
        }
        if let Some(count) = counts.get(&key).copied()
            && count > 1
        {
            result.push_str(&count.to_string());
        }
    }
    if charge > 0 {
        result.push('+');
        if charge > 1 {
            result.push_str(&charge.to_string());
        }
    } else if charge < 0 {
        result.push('-');
        if charge < -1 {
            result.push_str(&(-charge).to_string());
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_core::{AdjacencyList, Molecule};

    fn detached(smiles: &str) -> TopologyBlock {
        let molecule = Molecule::from_smiles(smiles).expect("reference SMILES must parse");
        TopologyBlock {
            atoms: molecule.atoms().to_vec(),
            bonds: molecule.bonds().to_vec(),
            adjacency: AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds()),
            ..TopologyBlock::default()
        }
    }

    #[test]
    fn detached_descriptors_match_legacy_core_for_ethanol() {
        let molecule = Molecule::from_smiles("CCO").expect("reference SMILES must parse");
        let topology = detached("CCO");
        let expected_average = cosmolkit_core::calc_mol_wt(&molecule, false).unwrap();
        let expected_exact = cosmolkit_core::calc_exact_mol_wt(&molecule, false).unwrap();
        let expected_formula = cosmolkit_core::calc_mol_formula(&molecule, false, false).unwrap();

        assert_eq!(molecular_formula(&topology).unwrap(), expected_formula);
        assert!((molecular_weight(&topology).unwrap() - expected_average).abs() < 1e-12);
        assert!((exact_molecular_weight(&topology).unwrap() - expected_exact).abs() < 1e-12);
    }
}
