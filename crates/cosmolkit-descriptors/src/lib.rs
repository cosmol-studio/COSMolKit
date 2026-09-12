//! Detached molecular mass, formula, and topology descriptor primitives.
//!
//! The descriptor boundary accepts only a validated [`TopologyBlock`]. It
//! cannot observe or mutate a live runtime molecule or derived-cache state.

use std::{cmp::Ordering, collections::BTreeMap};

use cosmolkit_core::{
    ValenceModel, assign_valence_with_options_for_topology, atomic_mass, most_common_isotope_mass,
    rdkit_element_symbol, total_hydrogen_count,
};
use cosmolkit_model::{Element, TopologyBlock};

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

fn descriptor_atomic_mass(
    element: Element,
    isotope: Option<u16>,
    function: &'static str,
) -> DescriptorResult<f64> {
    atomic_mass(element, isotope).map_err(|error| DescriptorError::Unsupported {
        function,
        detail: error.to_string(),
    })
}

/// Counts nitrogen and oxygen atoms using RDKit's direct Lipinski definition.
pub fn lipinski_hba(topology: &TopologyBlock) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBA
    // RDKit✔️✔️: unsigned int calcLipinskiHBA(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (ROMol::ConstAtomIterator iter = mol.beginAtoms(); iter != mol.endAtoms();
    // RDKit✔️✔️:        ++iter) {
    // RDKit✔️✔️:     if ((*iter)->getAtomicNum() == 7 || (*iter)->getAtomicNum() == 8) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBA
    validate_topology(topology, "lipinski_hba")?;
    u32::try_from(
        topology
            .atoms
            .iter()
            .filter(|atom| matches!(atom.atomic_number(), 7 | 8))
            .count(),
    )
    .map_err(|_| DescriptorError::Unsupported {
        function: "lipinski_hba",
        detail: "acceptor count exceeds the RDKit unsigned result model".to_owned(),
    })
}

/// Sums hydrogens on nitrogen and oxygen using RDKit's direct Lipinski definition.
pub fn lipinski_hbd(topology: &TopologyBlock) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBD
    // RDKit✔️✔️: unsigned int calcLipinskiHBD(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (ROMol::ConstAtomIterator iter = mol.beginAtoms(); iter != mol.endAtoms();
    // RDKit✔️✔️:        ++iter) {
    // RDKit✔️✔️:     if (((*iter)->getAtomicNum() == 7 || (*iter)->getAtomicNum() == 8)) {
    // RDKit✔️✔️:       res += (*iter)->getTotalNumHs(true);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBD
    validate_topology(topology, "lipinski_hbd")?;
    let assignment = valence(topology, "lipinski_hbd")?;
    let mut result = 0_u32;
    for atom in &topology.atoms {
        if matches!(atom.atomic_number(), 7 | 8) {
            let hydrogens =
                total_hydrogen_count(topology, &assignment, atom.id(), true).map_err(|error| {
                    DescriptorError::Unsupported {
                        function: "lipinski_hbd",
                        detail: error.to_string(),
                    }
                })?;
            result = result
                .checked_add(hydrogens)
                .ok_or_else(|| DescriptorError::Unsupported {
                    function: "lipinski_hbd",
                    detail: "donor hydrogen count exceeds the RDKit unsigned result model"
                        .to_owned(),
                })?;
        }
    }
    Ok(result)
}

/// Counts explicit atoms whose atomic number is greater than one.
pub fn num_heavy_atoms(topology: &TopologyBlock) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumHeavyAtoms
    // RDKit✔️✔️: unsigned int calcNumHeavyAtoms(const ROMol &mol) {
    // RDKit✔️✔️:   return mol.getNumHeavyAtoms();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumHeavyAtoms
    // BEGIN RDKIT CPP FUNCTION: RDKit::ROMol::getNumHeavyAtoms
    // RDKit✔️✔️: unsigned int ROMol::getNumHeavyAtoms() const {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto atom : atoms()) {
    // RDKit✔️✔️:     if (atom->getAtomicNum() > 1) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::ROMol::getNumHeavyAtoms
    validate_topology(topology, "num_heavy_atoms")?;
    u32::try_from(
        topology
            .atoms
            .iter()
            .filter(|atom| atom.atomic_number() > 1)
            .count(),
    )
    .map_err(|_| DescriptorError::Unsupported {
        function: "num_heavy_atoms",
        detail: "heavy-atom count exceeds the RDKit unsigned result model".to_owned(),
    })
}

/// Counts explicit atoms plus attached implicit/explicit atom-state hydrogens.
pub fn num_atoms(topology: &TopologyBlock) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAtoms
    // RDKit✔️✔️: unsigned int calcNumAtoms(const ROMol &mol) {
    // RDKit✔️✔️:   bool onlyExplicit = false;
    // RDKit✔️✔️:   return mol.getNumAtoms(onlyExplicit);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAtoms
    // BEGIN RDKIT CPP FUNCTION: RDKit::ROMol::getNumAtoms
    // RDKit✔️✔️: unsigned int ROMol::getNumAtoms(bool onlyExplicit) const {
    // RDKit✔️✔️:   int res = rdcast<int>(boost::num_vertices(d_graph));
    // RDKit✔️✔️:   if (!onlyExplicit) {
    // RDKit✔️✔️:     // if we are interested in hydrogens as well add them up from
    // RDKit✔️✔️:     // each
    // RDKit✔️✔️:     for (const auto atom : atoms()) {
    // RDKit✔️✔️:       res += atom->getTotalNumHs();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::ROMol::getNumAtoms
    validate_topology(topology, "num_atoms")?;
    let assignment = valence(topology, "num_atoms")?;
    let mut result =
        u32::try_from(topology.atoms.len()).map_err(|_| DescriptorError::Unsupported {
            function: "num_atoms",
            detail: "explicit atom count exceeds the RDKit unsigned result model".to_owned(),
        })?;
    for atom in &topology.atoms {
        let hydrogens =
            total_hydrogen_count(topology, &assignment, atom.id(), false).map_err(|error| {
                DescriptorError::Unsupported {
                    function: "num_atoms",
                    detail: error.to_string(),
                }
            })?;
        result = result
            .checked_add(hydrogens)
            .ok_or_else(|| DescriptorError::Unsupported {
                function: "num_atoms",
                detail: "total atom count exceeds the RDKit unsigned result model".to_owned(),
            })?;
    }
    Ok(result)
}

/// Fraction of carbon atoms whose total degree is four.
pub fn fraction_csp3(topology: &TopologyBlock) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcFractionCSP3
    // RDKit✔️✔️: double calcFractionCSP3(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int nCSP3 = 0;
    // RDKit✔️✔️:   unsigned int nC = 0;
    // RDKit✔️✔️:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️✔️:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️✔️:   while (atBegin != atEnd) {
    // RDKit✔️✔️:     const Atom *at = mol[*atBegin];
    // RDKit✔️✔️:     if (at->getAtomicNum() == 6) {
    // RDKit✔️✔️:       ++nC;
    // RDKit✔️✔️:       if (at->getTotalDegree() == 4) {
    // RDKit✔️✔️:         ++nCSP3;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atBegin;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!nC) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return static_cast<double>(nCSP3) / nC;
    // RDKit✔️✔️: }
    // BEGIN RDKIT CPP FUNCTION: RDKit::Atom::getTotalDegree
    // RDKit✔️✔️: unsigned int Atom::getTotalDegree() const {
    // RDKit✔️✔️:   unsigned int res = this->getTotalNumHs(false) + this->getDegree();
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Atom::getTotalDegree
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcFractionCSP3
    validate_topology(topology, "fraction_csp3")?;
    let assignment = valence(topology, "fraction_csp3")?;
    let mut carbon_count = 0_u32;
    let mut sp3_count = 0_u32;
    for atom in &topology.atoms {
        if atom.atomic_number() != 6 {
            continue;
        }
        carbon_count += 1;
        let hydrogens =
            total_hydrogen_count(topology, &assignment, atom.id(), false).map_err(|error| {
                DescriptorError::Unsupported {
                    function: "fraction_csp3",
                    detail: error.to_string(),
                }
            })?;
        let degree = u32::try_from(topology.adjacency.neighbors_of(atom.id().index()).len())
            .map_err(|_| DescriptorError::Unsupported {
                function: "fraction_csp3",
                detail: "atom degree exceeds the RDKit unsigned result model".to_owned(),
            })?;
        if degree.saturating_add(hydrogens) == 4 {
            sp3_count += 1;
        }
    }
    if carbon_count == 0 {
        return Ok(0.0);
    }
    Ok(f64::from(sp3_count) / f64::from(carbon_count))
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
    let hydrogen_mass = descriptor_atomic_mass(Element::H, None, "molecular_weight")?;
    let mut result = 0.0;
    for (index, atom) in topology.atoms.iter().enumerate() {
        if !only_heavy || atom.atomic_number() != 1 {
            result += descriptor_atomic_mass(atom.element(), atom.isotope(), "molecular_weight")?;
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
                most_common_isotope_mass(atom.element())
            } else {
                descriptor_atomic_mass(atom.element(), atom.isotope(), "exact_molecular_weight")?
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
        result += f64::from(hydrogens_to_count) * most_common_isotope_mass(Element::H);
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
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element,
    };

    fn detached_ethanol() -> TopologyBlock {
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            ),
        ];
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    fn smiles_topology(smiles: &str) -> TopologyBlock {
        cosmolkit_smiles::parse_smiles(smiles, &Default::default())
            .expect("parse detached descriptor fixture")
            .topology
    }

    #[test]
    fn detached_descriptors_match_rdkit_ethanol_fixture() {
        let topology = detached_ethanol();
        assert_eq!(molecular_formula(&topology).unwrap(), "C2H6O");
        assert!((molecular_weight(&topology).unwrap() - 46.069).abs() < 1e-9);
        assert!((exact_molecular_weight(&topology).unwrap() - 46.041864812).abs() < 1e-9);
    }

    #[test]
    fn detached_lipinski_and_atom_counts_match_pinned_rdkit_cases() {
        const CASES: [(&str, [u32; 4], f64); 6] = [
            ("CCO", [1, 1, 3, 9], 1.0),
            ("NC(=O)C", [2, 2, 4, 9], 0.5),
            ("NC(=O)N", [3, 4, 4, 8], 0.0),
            ("NCC(=O)O", [3, 3, 5, 10], 0.5),
            ("c1ncc[nH]1", [2, 1, 5, 9], 0.0),
            ("[Na+].[Cl-]", [0, 0, 2, 2], 0.0),
        ];

        for (smiles, expected, expected_fraction) in CASES {
            let topology = smiles_topology(smiles);
            assert_eq!(
                [
                    lipinski_hba(&topology).unwrap(),
                    lipinski_hbd(&topology).unwrap(),
                    num_heavy_atoms(&topology).unwrap(),
                    num_atoms(&topology).unwrap(),
                ],
                expected,
                "{smiles}"
            );
            assert_eq!(
                fraction_csp3(&topology).unwrap(),
                expected_fraction,
                "{smiles}"
            );
        }
    }

    #[test]
    fn detached_lipinski_donor_count_includes_explicit_hydrogen_neighbors() {
        let mut atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::N)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        atoms.extend(
            (2..4).map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::H))),
        );
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(2),
                BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
            ),
        ];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };

        assert_eq!(lipinski_hba(&topology), Ok(1));
        assert_eq!(lipinski_hbd(&topology), Ok(2));
        assert_eq!(num_heavy_atoms(&topology), Ok(2));
        assert_eq!(num_atoms(&topology), Ok(7));
        assert_eq!(fraction_csp3(&topology), Ok(1.0));
    }
}
