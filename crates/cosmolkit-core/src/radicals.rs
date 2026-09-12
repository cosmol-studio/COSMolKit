//! RDKit-aligned radical-electron assignment over detached topology values.

use cosmolkit_model::{Atom, AtomId, TopologyBlock, TopologyValidationError};

use crate::{periodic_table, valence};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RadicalAssignment {
    pub radical_electrons: Vec<u8>,
    pub diagnostics: Vec<RadicalDiagnostic>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RadicalDiagnostic {
    UnusualChargeClamped {
        atom: AtomId,
        atomic_number: u8,
        formal_charge: i8,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RadicalError {
    #[error("invalid topology: {source}")]
    InvalidTopology { source: TopologyValidationError },
    #[error("periodic-table field {field} is unavailable for atomic number {atomic_number}")]
    PeriodicTableLookup {
        atomic_number: u8,
        field: &'static str,
    },
    #[error(transparent)]
    Valence(#[from] valence::ValenceError),
    #[error("radical electron count out of range at atom {atom}: {count}")]
    RadicalCountOutOfRange { atom: AtomId, count: i32 },
}

/// Assign radical-electron counts in canonical atom-row order.
pub fn assign_radicals(topology: &TopologyBlock) -> Result<RadicalAssignment, RadicalError> {
    topology
        .validate()
        .map_err(|source| RadicalError::InvalidTopology { source })?;

    // BEGIN RDKIT CPP FUNCTION MolOps::assignRadicals
    // RDKit✔️❌: void assignRadicals(RWMol &mol) {
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     // we only automatically assign radicals to atoms that
    // RDKit✔️❌:     // don't have implicit Hs:
    // RDKit✔️❌:     if (!atom->getNoImplicit() || !atom->getAtomicNum()) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     const auto &valens =
    // RDKit✔️❌:         PeriodicTable::getTable()->getValenceList(atom->getAtomicNum());
    // RDKit✔️❌:     int chg = atom->getFormalCharge();
    // RDKit✔️❌:     int nOuter =
    // RDKit✔️❌:         PeriodicTable::getTable()->getNouterElecs(atom->getAtomicNum());
    // RDKit✔️❌:     if (valens.size() != 1 || valens[0] != -1) {
    // RDKit✔️❌:       double accum = 0.0;
    // RDKit✔️❌:       RWMol::OEDGE_ITER beg, end;
    // RDKit✔️❌:       boost::tie(beg, end) = mol.getAtomBonds(atom);
    // RDKit✔️❌:       while (beg != end) {
    // RDKit✔️❌:         accum += mol[*beg]->getValenceContrib(atom);
    // RDKit✔️❌:         ++beg;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       accum += atom->getNumExplicitHs();
    // RDKit✔️❌:       int totalValence = static_cast<int>(accum + 0.1);
    // RDKit✔️❌:       int baseCount = 8;
    // RDKit✔️❌:       if (atom->getAtomicNum() == 1 || atom->getAtomicNum() == 2) {
    // RDKit✔️❌:         baseCount = 2;
    // RDKit✔️❌:       }
    // RDKit✔️❌:
    // RDKit✔️❌:       // applies to later (more electronegative) elements:
    // RDKit✔️❌:       int numRadicals = baseCount - nOuter - totalValence + chg;
    // RDKit✔️❌:       if (numRadicals < 0) {
    // RDKit✔️❌:         numRadicals = 0;
    // RDKit✔️❌:         // can the atom be "hypervalent"?  (was github #447)
    // RDKit✔️❌:         const INT_VECT &valens =
    // RDKit✔️❌:             PeriodicTable::getTable()->getValenceList(atom->getAtomicNum());
    // RDKit✔️❌:         if (valens.size() > 1) {
    // RDKit✔️❌:           for (auto val : valens) {
    // RDKit✔️❌:             if (val - totalValence + chg >= 0) {
    // RDKit✔️❌:               numRadicals = val - totalValence + chg;
    // RDKit✔️❌:               break;
    // RDKit✔️❌:             }
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       // applies to earlier elements:
    // RDKit✔️❌:       int numRadicals2 = nOuter - totalValence - chg;
    // RDKit✔️❌:       if (numRadicals2 >= 0) {
    // RDKit✔️❌:         numRadicals = std::min(numRadicals, numRadicals2);
    // RDKit✔️❌:       }
    // RDKit✔️❌:       atom->setNumRadicalElectrons(numRadicals);
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       // #7122: if there's a bond to the metal center, then don't assign
    // RDKit✔️❌:       // radicals:
    // RDKit✔️❌:       if (atom->getDegree() > 0) {
    // RDKit✔️❌:         atom->setNumRadicalElectrons(0);
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         auto nValence = nOuter - chg;
    // RDKit✔️❌:         //  if this is an atom where we have no preferred valence info at all,
    // RDKit✔️❌:         //  e.g. for transition metals, then we shouldn't be guessing. This was
    // RDKit✔️❌:         //  #3330
    // RDKit✔️❌:         if (nValence < 0) {
    // RDKit✔️❌:           // this was github #5462
    // RDKit✔️❌:           nValence = 0;
    // RDKit✔️❌:           BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:               << "Unusual charge on atom " << atom->getIdx()
    // RDKit✔️❌:               << " number of radical electrons set to zero" << std::endl;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         atom->setNumRadicalElectrons(nValence % 2);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::assignRadicals
    // The detached API must allocate one ordered result row per atom, while
    // RDKit mutates existing atom fields in place. Both traverse atoms and
    // adjacency in O(V + E), but the required O(V) result allocation is a
    // material allocation difference, so the second marker remains `❌`.

    let mut radical_electrons = topology
        .atoms
        .iter()
        .map(Atom::radical_electrons)
        .collect::<Vec<_>>();
    let mut diagnostics = Vec::new();

    for atom in &topology.atoms {
        if !atom.no_implicit() || atom.atomic_number() == 0 {
            continue;
        }

        let valences = periodic_table::valences(atom.atomic_number()).ok_or(
            RadicalError::PeriodicTableLookup {
                atomic_number: atom.atomic_number(),
                field: "valences",
            },
        )?;
        let outer = periodic_table::outer_electrons(atom.atomic_number()).ok_or(
            RadicalError::PeriodicTableLookup {
                atomic_number: atom.atomic_number(),
                field: "outer_electrons",
            },
        )?;
        let charge = i32::from(atom.formal_charge());

        let radical_count = if valences.len() != 1 || valences[0] != -1 {
            let mut accumulation = 0.0;
            for bond in valence::incident(topology, atom.id())? {
                accumulation += valence::bond_valence_contrib(bond, atom.id())?;
            }
            accumulation += f64::from(atom.explicit_hydrogens());
            let total_valence = (accumulation + 0.1) as i32;
            let base_count = if matches!(atom.atomic_number(), 1 | 2) {
                2
            } else {
                8
            };

            let mut count = base_count - outer - total_valence + charge;
            if count < 0 {
                count = 0;
                if valences.len() > 1 {
                    for &allowed_valence in valences {
                        if allowed_valence - total_valence + charge >= 0 {
                            count = allowed_valence - total_valence + charge;
                            break;
                        }
                    }
                }
            }
            let alternate = outer - total_valence - charge;
            if alternate >= 0 {
                count = count.min(alternate);
            }
            count
        } else if !topology
            .adjacency
            .neighbors_of(atom.id().index())
            .is_empty()
        {
            0
        } else {
            let mut valence = outer - charge;
            if valence < 0 {
                diagnostics.push(RadicalDiagnostic::UnusualChargeClamped {
                    atom: atom.id(),
                    atomic_number: atom.atomic_number(),
                    formal_charge: atom.formal_charge(),
                });
                valence = 0;
            }
            valence % 2
        };

        radical_electrons[atom.id().index()] =
            u8::try_from(radical_count).map_err(|_| RadicalError::RadicalCountOutOfRange {
                atom: atom.id(),
                count: radical_count,
            })?;
    }

    Ok(RadicalAssignment {
        radical_electrons,
        diagnostics,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{AdjacencyList, AtomSpec};
    use cosmolkit_types::Element;

    #[test]
    fn detached_radical_assignment_handles_no_implicit_carbon() {
        let atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_radical_electrons(2),
        );
        let topology = TopologyBlock {
            atoms: vec![atom],
            adjacency: AdjacencyList::from_topology(1, &[]),
            ..TopologyBlock::default()
        };
        assert_eq!(
            assign_radicals(&topology).unwrap().radical_electrons,
            vec![4]
        );
    }
}
