// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use crate::{
    AdjacencyList, Atom, AtomId, Bond, BondOrder, BondQueryPredicate, Molecule, QueryNode,
    read_parts::MoleculeReadParts,
};

include!(concat!(
    env!("OUT_DIR"),
    "/rdkit_isotope_masses_generated.rs"
));

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValenceModel {
    RdkitLike,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceAssignment {
    pub explicit_valence: Vec<i32>,
    pub implicit_hydrogens: Vec<i32>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ValenceError {
    #[error("{message}")]
    InvalidValence {
        atom: crate::AtomId,
        atomic_number: u8,
        formal_charge: i8,
        message: String,
    },
    #[error("unsupported valence branch: {reason}")]
    UnsupportedBranch { reason: &'static str },
    #[error("radical electron count out of range at atom {atom}: {count}")]
    RadicalCountOutOfRange { atom: crate::AtomId, count: i32 },
    #[error("{message}")]
    BadBondType { message: String },
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AtomValenceState {
    explicit_valence: i32,
    implicit_valence: i32,
}

pub(crate) fn assign_valence_state_for_atom_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
) -> Result<(i32, i32), ValenceError> {
    let state = atom_update_property_cache(atoms, bonds, adjacency, atom_id, strict)?;
    Ok((state.explicit_valence, state.implicit_valence))
}

pub(crate) fn assign_implicit_valence_for_atom_from_parts_with_explicit_valence(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    explicit_valence: i32,
    strict: bool,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Atom::calcImplicitValence
    // RDKit✔️✔️: int Atom::calcImplicitValence(bool strict) {
    // RDKit✔️✔️:   if (d_explicitValence == -1) {
    // RDKit✔️✔️:     calcExplicitValence(strict);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool checkIt = false;
    // RDKit✔️✔️:   d_implicitValence = calculateImplicitValence(*this, strict, checkIt);
    // RDKit✔️✔️:   return d_implicitValence;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::calcImplicitValence
    atom_calc_implicit_valence(atoms, bonds, adjacency, atom_id, explicit_valence, strict)
}

fn atom_from_parts(atoms: &[Atom], atom_id: AtomId) -> Result<&Atom, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::getAtomWithIdx
    // RDKit✔️✔️: //! returns a pointer to a particular Atom
    // RDKit✔️✔️: Atom *getAtomWithIdx(unsigned int idx);
    // RDKit✔️✔️: //! \overload
    // RDKit✔️✔️: const Atom *getAtomWithIdx(unsigned int idx) const;
    // END RDKIT CPP FUNCTION ROMol::getAtomWithIdx
    atoms
        .get(atom_id.index())
        .ok_or(ValenceError::UnsupportedBranch {
            reason: "atom index out of range",
        })
}

fn incident_bonds_from_parts<'a>(
    atom_count: usize,
    bonds: &'a [Bond],
    adjacency: &'a AdjacencyList,
    atom_id: AtomId,
) -> Result<impl Iterator<Item = &'a Bond> + 'a, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::atomBonds
    // RDKit✔️✔️: CXXBondIterator<const MolGraph, Bond *const, MolGraph::out_edge_iterator>
    // RDKit✔️✔️: atomBonds(Atom const *at) const {
    // RDKit✔️✔️:   auto pr = getAtomBonds(at);
    // RDKit✔️✔️:   return {&d_graph, pr.first, pr.second};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ROMol::atomBonds
    if atom_id.index() >= atom_count {
        return Err(ValenceError::UnsupportedBranch {
            reason: "atom index out of range",
        });
    }
    let neighbors = adjacency.neighbors_of(atom_id.index());
    if neighbors
        .iter()
        .any(|neighbor| neighbor.bond.index() >= bonds.len())
    {
        return Err(ValenceError::UnsupportedBranch {
            reason: "topology adjacency bond index out of range",
        });
    }
    Ok(neighbors
        .iter()
        .map(move |neighbor| &bonds[neighbor.bond.index()]))
}

pub fn rdkit_valence_list(atomic_number: u8) -> Result<Option<&'static [i32]>, ValenceError> {
    rdkit_valence_list_for_atomic_number(atomic_number)
}

pub fn assign_radicals(molecule: &Molecule) -> Result<Vec<u8>, ValenceError> {
    assign_radicals_from_parts(
        molecule.atoms(),
        molecule.bonds(),
        &molecule.topology_block().adjacency,
    )
}

pub(crate) fn assign_radicals_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
) -> Result<Vec<u8>, ValenceError> {
    // COSMolKit uses value-style molecule state. This helper computes the
    // radical-electron values that RDKit writes into RWMol atoms; the registered
    // operation applies them atomically through OpParts.
    // BEGIN RDKIT CPP FUNCTION MolOps::assignRadicals
    // RDKit✔️✔️: void assignRadicals(RWMol &mol) {
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     // we only automatically assign radicals to atoms that
    // RDKit✔️✔️:     // don't have implicit Hs:
    // RDKit✔️✔️:     if (!atom->getNoImplicit() || !atom->getAtomicNum()) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    let mut radicals = atoms
        .iter()
        .map(Atom::radical_electrons)
        .collect::<Vec<_>>();

    for atom in atoms {
        if !atom.no_implicit() || atom.atomic_number() == 0 {
            continue;
        }

        // RDKit✔️✔️:     const auto &valens =
        // RDKit✔️✔️:         PeriodicTable::getTable()->getValenceList(atom->getAtomicNum());
        // RDKit✔️✔️:     int chg = atom->getFormalCharge();
        // RDKit✔️✔️:     int nOuter =
        // RDKit✔️✔️:         PeriodicTable::getTable()->getNouterElecs(atom->getAtomicNum());
        let valens = required_valence_list(atom.atomic_number())?;
        let chg = i32::from(atom.formal_charge());
        let n_outer = periodic_table_outer_electrons(atom.atomic_number())?;

        // RDKit✔️✔️:     if (valens.size() != 1 || valens[0] != -1) {
        let num_radicals = if valens.len() != 1 || valens[0] != -1 {
            // RDKit✔️✔️:       double accum = 0.0;
            // RDKit✔️✔️:       RWMol::OEDGE_ITER beg, end;
            // RDKit✔️✔️:       boost::tie(beg, end) = mol.getAtomBonds(atom);
            // RDKit✔️✔️:       while (beg != end) {
            // RDKit✔️✔️:         accum += mol[*beg]->getValenceContrib(atom);
            // RDKit✔️✔️:         ++beg;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       accum += atom->getNumExplicitHs();
            // RDKit✔️✔️:       int totalValence = static_cast<int>(accum + 0.1);
            let mut accum = 0.0;
            for bond in incident_bonds_from_parts(atoms.len(), bonds, adjacency, atom.id())? {
                accum += bond_valence_contrib(bond, atom.id())?;
            }
            accum += f64::from(atom.explicit_hydrogens());
            let total_valence = (accum + 0.1) as i32;

            // RDKit✔️✔️:       int baseCount = 8;
            // RDKit✔️✔️:       if (atom->getAtomicNum() == 1 || atom->getAtomicNum() == 2) {
            // RDKit✔️✔️:         baseCount = 2;
            // RDKit✔️✔️:       }
            // RDKit uses the two-electron shell for H/He and the octet count
            // for the later main-group branch.
            let base_count = if matches!(atom.atomic_number(), 1 | 2) {
                2
            } else {
                8
            };

            // RDKit✔️✔️:       // applies to later (more electronegative) elements:
            // RDKit✔️✔️:       int numRadicals = baseCount - nOuter - totalValence + chg;
            let mut num_radicals = base_count - n_outer - total_valence + chg;
            // RDKit✔️✔️:       if (numRadicals < 0) {
            // RDKit✔️✔️:         numRadicals = 0;
            // RDKit✔️✔️:         // can the atom be "hypervalent"?  (was github #447)
            // RDKit✔️✔️:         const INT_VECT &valens =
            // RDKit✔️✔️:             PeriodicTable::getTable()->getValenceList(atom->getAtomicNum());
            // RDKit✔️✔️:         if (valens.size() > 1) {
            // RDKit✔️✔️:           for (auto val : valens) {
            // RDKit✔️✔️:             if (val - totalValence + chg >= 0) {
            // RDKit✔️✔️:               numRadicals = val - totalValence + chg;
            // RDKit✔️✔️:               break;
            // RDKit✔️✔️:             }
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            if num_radicals < 0 {
                num_radicals = 0;
                if valens.len() > 1 {
                    for &valence in valens {
                        if valence - total_valence + chg >= 0 {
                            num_radicals = valence - total_valence + chg;
                            break;
                        }
                    }
                }
            }
            // RDKit✔️✔️:       // applies to earlier elements:
            // RDKit✔️✔️:       int numRadicals2 = nOuter - totalValence - chg;
            // RDKit✔️✔️:       if (numRadicals2 >= 0) {
            // RDKit✔️✔️:         numRadicals = std::min(numRadicals, numRadicals2);
            // RDKit✔️✔️:       }
            let num_radicals2 = n_outer - total_valence - chg;
            if num_radicals2 >= 0 {
                num_radicals = num_radicals.min(num_radicals2);
            }
            // RDKit✔️✔️:       atom->setNumRadicalElectrons(numRadicals);
            num_radicals
        } else {
            // RDKit✔️✔️:     } else {
            // RDKit✔️✔️:       // #7122: if there's a bond to the metal center, then don't assign
            // RDKit✔️✔️:       // radicals:
            // RDKit✔️✔️:       if (atom->getDegree() > 0) {
            // RDKit✔️✔️:         atom->setNumRadicalElectrons(0);
            if incident_bonds_from_parts(atoms.len(), bonds, adjacency, atom.id())?
                .next()
                .is_some()
            {
                0
            } else {
                // RDKit✔️✔️:       } else {
                // RDKit✔️✔️:         auto nValence = nOuter - chg;
                // RDKit✔️✔️:         //  if this is an atom where we have no preferred valence info at all,
                // RDKit✔️✔️:         //  e.g. for transition metals, then we shouldn't be guessing. This was
                // RDKit✔️✔️:         //  #3330
                let mut n_valence = n_outer - chg;
                // RDKit✔️✔️:         if (nValence < 0) {
                // RDKit✔️✔️:           // this was github #5462
                // RDKit✔️✔️:           nValence = 0;
                // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
                // RDKit✔️✔️:               << "Unusual charge on atom " << atom->getIdx()
                // RDKit✔️✔️:               << " number of radical electrons set to zero" << std::endl;
                // RDKit✔️✔️:         }
                if n_valence < 0 {
                    n_valence = 0;
                    log::warn!(
                        "Unusual charge on atom {} number of radical electrons set to zero",
                        atom.id()
                    );
                }
                // RDKit✔️✔️:         atom->setNumRadicalElectrons(nValence % 2);
                // RDKit✔️✔️:       }
                // RDKit✔️✔️:     }
                n_valence % 2
            }
        };

        radicals[atom.id().index()] =
            u8::try_from(num_radicals).map_err(|_| ValenceError::RadicalCountOutOfRange {
                atom: atom.id(),
                count: num_radicals,
            })?;
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::assignRadicals
    Ok(radicals)
}

pub fn assign_valence(
    molecule: &Molecule,
    model: ValenceModel,
) -> Result<ValenceAssignment, ValenceError> {
    assign_valence_with_options(molecule, model, true)
}

pub fn assign_valence_with_options(
    molecule: &Molecule,
    model: ValenceModel,
    strict: bool,
) -> Result<ValenceAssignment, ValenceError> {
    assign_valence_with_options_from_parts(
        molecule.atoms(),
        molecule.bonds(),
        &molecule.topology_block().adjacency,
        model,
        strict,
    )
}

pub(crate) fn assign_valence_with_options_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    model: ValenceModel,
    strict: bool,
) -> Result<ValenceAssignment, ValenceError> {
    match model {
        ValenceModel::RdkitLike => romol_update_property_cache(atoms, bonds, adjacency, strict),
    }
}

pub fn atom_has_valence_violation(
    molecule: &Molecule,
    atom_id: AtomId,
) -> Result<bool, ValenceError> {
    atom_has_valence_violation_from_parts(
        molecule.atoms(),
        molecule.bonds(),
        &molecule.topology_block().adjacency,
        atom_id,
    )
}

fn romol_update_property_cache(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    strict: bool,
) -> Result<ValenceAssignment, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::updatePropertyCache
    // RDKit✔️✔️: void ROMol::updatePropertyCache(bool strict) {
    // RDKit✔️✔️:   for (auto atom : atoms()) {
    // RDKit✔️✔️:     atom->updatePropertyCache(strict);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto bond : bonds()) {
    // RDKit✔️✔️:     bond->updatePropertyCache(strict);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ROMol::updatePropertyCache
    let mut explicit_valence = Vec::with_capacity(atoms.len());
    let mut implicit_hydrogens = Vec::with_capacity(atoms.len());
    for atom in atoms {
        let state = atom_update_property_cache(atoms, bonds, adjacency, atom.id(), strict)?;
        explicit_valence.push(state.explicit_valence);
        implicit_hydrogens.push(state.implicit_valence);
    }
    for bond in bonds {
        bond_update_property_cache(bond, strict);
    }
    Ok(ValenceAssignment {
        explicit_valence,
        implicit_hydrogens,
    })
}

fn atom_update_property_cache(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
) -> Result<AtomValenceState, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Atom::updatePropertyCache
    // RDKit✔️✔️: void Atom::updatePropertyCache(bool strict) {
    // RDKit✔️✔️:   calcExplicitValence(strict);
    // RDKit✔️✔️:   calcImplicitValence(strict);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::updatePropertyCache
    let explicit_valence = atom_calc_explicit_valence(atoms, bonds, adjacency, atom_id, strict)?;
    let implicit_valence =
        atom_calc_implicit_valence(atoms, bonds, adjacency, atom_id, explicit_valence, strict)?;
    Ok(AtomValenceState {
        explicit_valence,
        implicit_valence,
    })
}

fn bond_update_property_cache(_bond: &Bond, strict: bool) {
    // BEGIN RDKIT CPP FUNCTION Bond::updatePropertyCache
    // RDKit✔️✔️: void updatePropertyCache(bool strict = true) { (void)strict; }
    // END RDKIT CPP FUNCTION Bond::updatePropertyCache
    let _ = strict;
}

fn atom_calc_explicit_valence(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Atom::calcExplicitValence
    // RDKit✔️✔️: int Atom::calcExplicitValence(bool strict) {
    // RDKit✔️✔️:   bool checkIt = false;
    // RDKit✔️✔️:   d_explicitValence = calculateExplicitValence(*this, strict, checkIt);
    // RDKit✔️✔️:   return d_explicitValence;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::calcExplicitValence
    calculate_explicit_valence(atoms, bonds, adjacency, atom_id, strict, false)
}

pub(crate) fn assign_explicit_valence_for_atom_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
) -> Result<i32, ValenceError> {
    atom_calc_explicit_valence(atoms, bonds, adjacency, atom_id, strict)
}

fn atom_calc_implicit_valence(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    explicit_valence: i32,
    strict: bool,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Atom::calcImplicitValence
    // RDKit✔️✔️: int Atom::calcImplicitValence(bool strict) {
    // RDKit✔️✔️:   if (d_explicitValence == -1) {
    // RDKit✔️✔️:     calcExplicitValence(strict);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool checkIt = false;
    // RDKit✔️✔️:   d_implicitValence = calculateImplicitValence(*this, strict, checkIt);
    // RDKit✔️✔️:   return d_implicitValence;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::calcImplicitValence
    calculate_implicit_valence(
        atoms,
        bonds,
        adjacency,
        atom_id,
        explicit_valence,
        strict,
        false,
    )
}

fn calculate_explicit_valence(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    strict: bool,
    check_it: bool,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION calculateExplicitValence
    // RDKit✔️✔️: int calculateExplicitValence(const Atom &atom, bool strict, bool checkIt) {
    // RDKit✔️✔️:   // FIX: contributions of bonds to valence are being done at best
    // RDKit✔️✔️:   // approximately
    // RDKit✔️✔️:   double accum = 0;
    // RDKit✔️✔️:   for (const auto bnd : atom.getOwningMol().atomBonds(&atom)) {
    // RDKit✔️✔️:     accum += bnd->getValenceContrib(&atom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   accum += atom.getNumExplicitHs();
    let atom = atom_from_parts(atoms, atom_id)?;
    let mut accum = 0.0;
    for bond in incident_bonds_from_parts(atoms.len(), bonds, adjacency, atom_id)? {
        accum += bond_valence_contrib(bond, atom_id)?;
    }
    accum += f64::from(atom.explicit_hydrogens());

    // RDKit✔️✔️:   const auto &ovalens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(atom.getAtomicNum());
    // RDKit✔️✔️:   // if we start with an atom that doesn't have specified valences, we stick
    // RDKit✔️✔️:   // with that. otherwise we will use the effective valence
    // RDKit✔️✔️:   unsigned int effectiveAtomicNum = atom.getAtomicNum();
    // RDKit✔️✔️:   if (ovalens.size() > 1 || ovalens[0] != -1) {
    // RDKit✔️✔️:     effectiveAtomicNum = getEffectiveAtomicNum(atom, checkIt);
    // RDKit✔️✔️:   }
    let ovalens = required_valence_list(atom.atomic_number())?;
    let mut effective_atomic_num = atom.atomic_number();
    if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_atomic_num = get_effective_atomic_num(atom, check_it)?;
    }

    // RDKit✔️✔️:   unsigned int dv =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getDefaultValence(effectiveAtomicNum);
    // RDKit✔️✔️:   const auto &valens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(effectiveAtomicNum);
    let default_valence = rdkit_default_valence(effective_atomic_num)?;
    let valens = required_valence_list(effective_atomic_num)?;

    // RDKit✔️✔️:   if (accum > dv && isAromaticAtom(atom)) {
    // RDKit✔️✔️:     // this needs some explanation : if the atom is aromatic and
    // RDKit✔️✔️:     // accum > dv we assume that no hydrogen can be added
    // RDKit✔️✔️:     // to this atom.  We set x = (v + chr) such that x is the
    // RDKit✔️✔️:     // closest possible integer to "accum" but less than
    // RDKit✔️✔️:     // "accum".
    // RDKit✔️✔️:     //
    // RDKit✔️✔️:     // "v" here is one of the allowed valences. For example:
    // RDKit✔️✔️:     //    sulfur here : O=c1ccs(=O)cc1
    // RDKit✔️✔️:     //    nitrogen here : c1cccn1C
    if accum > f64::from(default_valence)
        && is_aromatic_atom_from_parts(atoms, bonds, adjacency, atom_id)?
    {
        // RDKit✔️✔️:     int pval = dv;
        // RDKit✔️✔️:     for (auto val : valens) {
        // RDKit✔️✔️:       if (val == -1) {
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (val > accum) {
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         pval = val;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mut pval = default_valence;
        for &val in valens {
            if val == -1 || f64::from(val) > accum {
                break;
            }
            pval = val;
        }
        // RDKit✔️✔️:     // if we're within 1.5 of the allowed valence, go ahead and take it.
        // RDKit✔️✔️:     // this reflects things like the N in c1cccn1C, which starts with
        // RDKit✔️✔️:     // accum of 4, but which can be kekulized to C1=CC=CN1C, where
        // RDKit✔️✔️:     // the valence is 3 or the bridging N in c1ccn2cncc2c1, which starts
        // RDKit✔️✔️:     // with a valence of 4.5, but can be happily kekulized down to a valence
        // RDKit✔️✔️:     // of 3
        // RDKit✔️✔️:     if (accum - pval <= 1.5) {
        // RDKit✔️✔️:       accum = pval;
        // RDKit✔️✔️:     }
        if accum - f64::from(pval) <= 1.5 {
            accum = f64::from(pval);
        }
    }
    // RDKit✔️✔️:   // despite promising to not to blame it on him - this a trick Greg
    // RDKit✔️✔️:   // came up with: if we have a bond order sum of x.5 (i.e. 1.5, 2.5
    // RDKit✔️✔️:   // etc) we would like it to round to the higher integer value --
    // RDKit✔️✔️:   // 2.5 to 3 instead of 2 -- so we will add 0.1 to accum.
    // RDKit✔️✔️:   // this plays a role in the number of hydrogen that are implicitly
    // RDKit✔️✔️:   // added. This will only happen when the accum is a non-integer
    // RDKit✔️✔️:   // value and less than the default valence (otherwise the above if
    // RDKit✔️✔️:   // statement should have caught it). An example of where this can
    // RDKit✔️✔️:   // happen is the following smiles:
    // RDKit✔️✔️:   //    C1ccccC1
    // RDKit✔️✔️:   // Daylight accepts this smiles and we should be able to Kekulize
    // RDKit✔️✔️:   // correctly.
    // RDKit✔️✔️:   accum += 0.1;
    // RDKit✔️✔️:   auto res = static_cast<int>(std::round(accum));
    accum += 0.1;
    let res = accum.round() as i32;

    // RDKit✔️✔️:   if (strict || checkIt) {
    // RDKit✔️✔️:     int maxValence = valens.back();
    // RDKit✔️✔️:     int offset = 0;
    if strict || check_it {
        let mut max_valence = *valens.last().expect("valence list is nonempty");
        let mut offset = 0;
        // RDKit✔️✔️:     // we have to include a special case here for negatively charged P, S, As,
        // RDKit✔️✔️:     // and Se, which all support "hypervalent" forms, but which can be
        // RDKit✔️✔️:     // isoelectronic to Cl/Ar or Br/Kr, which do not support hypervalent forms.
        // RDKit✔️✔️:     if (canBeHypervalent(atom, effectiveAtomicNum)) {
        // RDKit✔️✔️:       maxValence = ovalens.back();
        // RDKit✔️✔️:       offset -= atom.getFormalCharge();
        // RDKit✔️✔️:     }
        if can_be_hypervalent(atom, effective_atomic_num) {
            max_valence = *ovalens.last().expect("valence list is nonempty");
            offset -= i32::from(atom.formal_charge());
        }
        // RDKit✔️✔️:     // we have historically accepted two-coordinate [H-] as a valid atom. This
        // RDKit✔️✔️:     // is highly questionable, but changing it requires some thought. For now we
        // RDKit✔️✔️:     // will just keep accepting it
        // RDKit✔️✔️:     if (atom.getAtomicNum() == 1 && atom.getFormalCharge() == -1) {
        // RDKit✔️✔️:       maxValence = 2;
        // RDKit✔️✔️:     }
        if atom.atomic_number() == 1 && atom.formal_charge() == -1 {
            max_valence = 2;
        }
        // RDKit✔️✔️:     // maxValence == -1 signifies that we'll take anything at the high end
        // RDKit✔️✔️:     if (maxValence >= 0 && ovalens.back() >= 0 && (res + offset) > maxValence) {
        // RDKit✔️✔️:       // the explicit valence is greater than any
        // RDKit✔️✔️:       // allowed valence for the atoms
        // RDKit✔️✔️:       if (strict) {
        // RDKit✔️✔️:         // raise an error
        // RDKit✔️✔️:         std::ostringstream errout;
        // RDKit✔️✔️:         errout << "Explicit valence for atom # " << atom.getIdx() << " "
        // RDKit✔️✔️:                << PeriodicTable::getTable()->getElementSymbol(
        // RDKit✔️✔️:                       atom.getAtomicNum())
        // RDKit✔️✔️:                << ", " << res << ", is greater than permitted";
        // RDKit✔️✔️:         std::string msg = errout.str();
        // RDKit✔️✔️:         BOOST_LOG(rdErrorLog) << msg << std::endl;
        // RDKit✔️✔️:         throw AtomValenceException(msg, atom.getIdx());
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         return -1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        if max_valence >= 0
            && *ovalens.last().expect("valence list is nonempty") >= 0
            && (res + offset) > max_valence
        {
            if strict {
                return Err(invalid_valence_with_message(
                    atom,
                    format!(
                        "Explicit valence for atom # {} {}, {}, is greater than permitted",
                        atom.id(),
                        rdkit_element_symbol(atom.atomic_number())?,
                        res
                    ),
                ));
            }
            // RDKit keeps the non-strict/check-it path on the `-1` sentinel so
            // later implicit-valence logic can distinguish "too high" without
            // throwing.
            return Ok(-1);
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION calculateExplicitValence
    Ok(res)
}

fn calculate_implicit_valence(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
    explicit_valence: i32,
    strict: bool,
    check_it: bool,
) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION calculateImplicitValence
    // RDKit✔️✔️: int calculateImplicitValence(const Atom &atom, bool strict, bool checkIt) {
    // RDKit✔️✔️:   if (atom.df_noImplicit) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    let atom = atom_from_parts(atoms, atom_id)?;
    if atom.no_implicit() {
        return Ok(0);
    }
    // RDKit✔️✔️:   auto explicitValence = atom.d_explicitValence;
    // RDKit✔️✔️:   if (explicitValence == -1) {
    // RDKit✔️✔️:     explicitValence = calculateExplicitValence(atom, strict, checkIt);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // special cases
    // RDKit✔️✔️:   auto atomicNum = atom.d_atomicNum;
    // RDKit✔️✔️:   if (atomicNum == 0) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    let atomic_num = atom.atomic_number();
    if atomic_num == 0 {
        return Ok(0);
    }
    // RDKit✔️✔️:   for (const auto bnd : atom.getOwningMol().atomBonds(&atom)) {
    // RDKit✔️✔️:     if (QueryOps::hasComplexBondTypeQuery(*bnd)) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for bond in incident_bonds_from_parts(atoms.len(), bonds, adjacency, atom_id)? {
        if bond.query().is_some_and(has_complex_bond_type_query) {
            return Ok(0);
        }
    }
    // RDKit✔️✔️:   auto formalCharge = atom.d_formalCharge;
    // RDKit✔️✔️:   auto numRadicalElectrons = atom.d_numRadicalElectrons;
    // RDKit✔️✔️:   if (explicitValence == 0 && numRadicalElectrons == 0 && atomicNum == 1) {
    // RDKit✔️✔️:     if (formalCharge == 1 || formalCharge == -1) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     } else if (formalCharge == 0) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if (strict) {
    // RDKit✔️✔️:         std::ostringstream errout;
    // RDKit✔️✔️:         errout << "Unreasonable formal charge on atom # " << atom.getIdx()
    // RDKit✔️✔️:                << ".";
    // RDKit✔️✔️:         std::string msg = errout.str();
    // RDKit✔️✔️:         BOOST_LOG(rdErrorLog) << msg << std::endl;
    // RDKit✔️✔️:         throw AtomValenceException(msg, atom.getIdx());
    // RDKit✔️✔️:       } else if (checkIt) {
    // RDKit✔️✔️:         return -1;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         return 0;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if explicit_valence == 0 && atom.radical_electrons() == 0 && atomic_num == 1 {
        return match atom.formal_charge() {
            1 | -1 => Ok(0),
            0 => Ok(1),
            _ if strict => Err(invalid_valence_with_message(
                atom,
                format!("Unreasonable formal charge on atom # {}.", atom.id()),
            )),
            _ if check_it => Ok(-1),
            _ => Ok(0),
        };
    }
    // RDKit✔️✔️:   int explicitPlusRadV = atom.d_explicitValence + atom.d_numRadicalElectrons;
    let mut explicit_plus_rad_v = explicit_valence + i32::from(atom.radical_electrons());

    // RDKit✔️✔️:   const auto &ovalens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(atom.d_atomicNum);
    // RDKit✔️✔️:   // if we start with an atom that doesn't have specified valences, we stick
    // RDKit✔️✔️:   // with that. otherwise we will use the effective valence for the rest of
    // RDKit✔️✔️:   // this.
    // RDKit✔️✔️:   unsigned int effectiveAtomicNum = atom.d_atomicNum;
    // RDKit✔️✔️:   if (ovalens.size() > 1 || ovalens[0] != -1) {
    // RDKit✔️✔️:     effectiveAtomicNum = getEffectiveAtomicNum(atom, checkIt);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (effectiveAtomicNum == 0) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    let ovalens = required_valence_list(atomic_num)?;
    let mut effective_atomic_num = atomic_num;
    if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_atomic_num = get_effective_atomic_num(atom, check_it)?;
    }
    if effective_atomic_num == 0 {
        return Ok(0);
    }

    // RDKit✔️✔️:   // The d-block and f-block of the periodic table (i.e. transition metals,
    // RDKit✔️✔️:   // lanthanoids and actinoids) have no default valence.
    // RDKit✔️✔️:   int dv = PeriodicTable::getTable()->getDefaultValence(effectiveAtomicNum);
    // RDKit✔️✔️:   if (dv == -1) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    let default_valence = rdkit_default_valence(effective_atomic_num)?;
    if default_valence == -1 {
        return Ok(0);
    }

    // RDKit✔️✔️:   // here is how we are going to deal with the possibility of
    // RDKit✔️✔️:   // multiple valences
    // RDKit✔️✔️:   // - check the explicit valence "ev"
    // RDKit✔️✔️:   // - if it is already equal to one of the allowed valences for the
    // RDKit✔️✔️:   //    atom return 0
    // RDKit✔️✔️:   // - otherwise take return difference between next larger allowed
    // RDKit✔️✔️:   //   valence and "ev"
    // RDKit✔️✔️:   // if "ev" is greater than all allowed valences for the atom raise an
    // RDKit✔️✔️:   // exception
    // RDKit✔️✔️:   // finally aromatic cases are dealt with differently - these atoms are allowed
    // RDKit✔️✔️:   // only default valences
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // we have to include a special case here for negatively charged P, S, As,
    // RDKit✔️✔️:   // and Se, which all support "hypervalent" forms, but which can be
    // RDKit✔️✔️:   // isoelectronic to Cl/Ar or Br/Kr, which do not support hypervalent forms.
    // RDKit✔️✔️:   if (canBeHypervalent(atom, effectiveAtomicNum)) {
    // RDKit✔️✔️:     effectiveAtomicNum = atomicNum;
    // RDKit✔️✔️:     explicitPlusRadV -= atom.d_formalCharge;
    // RDKit✔️✔️:   }
    if can_be_hypervalent(atom, effective_atomic_num) {
        // RDKit switches back to the original atomic number here so the
        // allowed-valence scan uses the hypervalent atom's native table row.
        effective_atomic_num = atomic_num;
        explicit_plus_rad_v -= i32::from(atom.formal_charge());
    }
    // RDKit✔️✔️:   const auto &valens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(effectiveAtomicNum);
    let valens = required_valence_list(effective_atomic_num)?;

    // RDKit✔️✔️:   int res = 0;
    let res;
    // RDKit✔️✔️:   // if we have an aromatic case treat it differently
    // RDKit✔️✔️:   if (isAromaticAtom(atom)) {
    if is_aromatic_atom_from_parts(atoms, bonds, adjacency, atom_id)? {
        // RDKit✔️✔️:     if (explicitPlusRadV <= dv) {
        // RDKit✔️✔️:       res = dv - explicitPlusRadV;
        // RDKit✔️✔️:     } else {
        if explicit_plus_rad_v <= default_valence {
            res = default_valence - explicit_plus_rad_v;
        } else {
            // RDKit✔️✔️:       // As we assume when finding the explicitPlusRadValence if we are
            // RDKit✔️✔️:       // aromatic we should not be adding any hydrogen and already
            // RDKit✔️✔️:       // be at an accepted valence state,
            // RDKit✔️✔️:
            // RDKit✔️✔️:       // FIX: this is just ERROR checking and probably moot - the
            // RDKit✔️✔️:       // explicitPlusRadValence function called above should assure us that
            // RDKit✔️✔️:       // we satisfy one of the accepted valence states for the
            // RDKit✔️✔️:       // atom. The only diff I can think of is in the way we handle
            // RDKit✔️✔️:       // formal charge here vs the explicit valence function.
            // RDKit✔️✔️:       bool satis = false;
            // RDKit✔️✔️:       for (auto vi = valens.begin(); vi != valens.end() && *vi > 0; ++vi) {
            // RDKit✔️✔️:         if (explicitPlusRadV == *vi) {
            // RDKit✔️✔️:           satis = true;
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            let satisfied = valens
                .iter()
                .take_while(|&&val| val > 0)
                .any(|&val| explicit_plus_rad_v == val);
            // RDKit✔️✔️:       if (!satis && (strict || checkIt)) {
            // RDKit✔️✔️:         if (strict) {
            // RDKit✔️✔️:           std::ostringstream errout;
            // RDKit✔️✔️:           errout << "Explicit valence for aromatic atom # " << atom.getIdx()
            // RDKit✔️✔️:                  << " not equal to any accepted valence\n";
            // RDKit✔️✔️:           std::string msg = errout.str();
            // RDKit✔️✔️:           BOOST_LOG(rdErrorLog) << msg << std::endl;
            // RDKit✔️✔️:           throw AtomValenceException(msg, atom.getIdx());
            // RDKit✔️✔️:         } else {
            // RDKit✔️✔️:           return -1;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            if !satisfied && (strict || check_it) {
                if strict {
                    return Err(invalid_valence_with_message(
                        atom,
                        format!(
                            "Explicit valence for aromatic atom # {} not equal to any accepted valence\n",
                            atom.id()
                        ),
                    ));
                }
                return Ok(-1);
            }
            // RDKit✔️✔️:       res = 0;
            res = 0;
        }
    } else {
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     // non-aromatic case we are allowed to have non default valences
        // RDKit✔️✔️:     // and be able to add hydrogens
        // RDKit✔️✔️:     res = -1;
        // RDKit✔️✔️:     for (auto vi = valens.begin(); vi != valens.end() && *vi >= 0; ++vi) {
        // RDKit✔️✔️:       int tot = *vi;
        // RDKit✔️✔️:       if (explicitPlusRadV <= tot) {
        // RDKit✔️✔️:         res = tot - explicitPlusRadV;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mut candidate = -1;
        for &valence in valens.iter().take_while(|&&val| val >= 0) {
            if explicit_plus_rad_v <= valence {
                candidate = valence - explicit_plus_rad_v;
                break;
            }
        }
        // RDKit✔️✔️:     if (res < 0) {
        // RDKit✔️✔️:       if ((strict || checkIt) && valens.back() != -1 && ovalens.back() > 0) {
        // RDKit✔️✔️:         // this means that the explicit valence is greater than any
        // RDKit✔️✔️:         // allowed valence for the atoms
        // RDKit✔️✔️:         if (strict) {
        // RDKit✔️✔️:           // raise an error
        // RDKit✔️✔️:           std::ostringstream errout;
        // RDKit✔️✔️:           errout << "Explicit valence for atom # " << atom.getIdx() << " "
        // RDKit✔️✔️:                  << PeriodicTable::getTable()->getElementSymbol(atomicNum)
        // RDKit✔️✔️:                  << " greater than permitted";
        // RDKit✔️✔️:           std::string msg = errout.str();
        // RDKit✔️✔️:           BOOST_LOG(rdErrorLog) << msg << std::endl;
        // RDKit✔️✔️:           throw AtomValenceException(msg, atom.getIdx());
        // RDKit✔️✔️:         } else {
        // RDKit✔️✔️:           return -1;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         res = 0;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        if candidate < 0 {
            if (strict || check_it)
                && *valens.last().expect("valence list is nonempty") != -1
                && *ovalens.last().expect("valence list is nonempty") > 0
            {
                if strict {
                    return Err(invalid_valence_with_message(
                        atom,
                        format!(
                            "Explicit valence for atom # {} {} greater than permitted",
                            atom.id(),
                            rdkit_element_symbol(atomic_num)?
                        ),
                    ));
                }
                return Ok(-1);
            }
            candidate = 0;
        }
        res = candidate;
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION calculateImplicitValence
    Ok(res)
}

fn atom_has_valence_violation_impl(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
) -> Result<bool, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Atom::hasValenceViolation
    // RDKit✔️✔️: bool Atom::hasValenceViolation() const {
    // RDKit✔️✔️:   // Ignore dummy atoms, query atoms, or atoms attached to query bonds
    // RDKit✔️✔️:   auto bonds = getOwningMol().atomBonds(this);
    // RDKit✔️✔️:   auto is_query = [](auto b) { return b->hasQuery(); };
    // RDKit✔️✔️:   if (getAtomicNum() == 0 || hasQuery() ||
    // RDKit✔️✔️:       std::any_of(bonds.begin(), bonds.end(), is_query)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    let atom = atom_from_parts(atoms, atom_id)?;
    if atom.atomic_number() == 0
        || atom.query().is_some()
        || incident_bonds_from_parts(atoms.len(), bonds, adjacency, atom_id)?
            .any(|bond| bond.query().is_some())
    {
        return Ok(false);
    }

    // RDKit✔️✔️:   unsigned int effectiveAtomicNum;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     bool checkIt = true;
    // RDKit✔️✔️:     effectiveAtomicNum = getEffectiveAtomicNum(*this, checkIt);
    // RDKit✔️✔️:   } catch (const AtomValenceException &) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    let effective_atomic_num = match get_effective_atomic_num(atom, true) {
        // RDKit treats an out-of-range effective atomic number as an immediate
        // valence violation instead of propagating the exception shape.
        Ok(value) => value,
        Err(ValenceError::InvalidValence { .. }) => return Ok(true),
        Err(error) => return Err(error),
    };

    // RDKit✔️✔️:   // special case for H:
    // RDKit✔️✔️:   if (getAtomicNum() == 1) {
    // RDKit✔️✔️:     if (getFormalCharge() > 1 || getFormalCharge() < -1) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // Non-H checks for absurd charge values:
    // RDKit✔️✔️:     //   1. the formal charge is larger than the atomic number
    // RDKit✔️✔️:     //   2. the formal charge moves us to a different row of the periodic table
    // RDKit✔️✔️:     if (getFormalCharge() > getAtomicNum() ||
    // RDKit✔️✔️:         PeriodicTable::getTable()->getRow(d_atomicNum) !=
    // RDKit✔️✔️:             PeriodicTable::getTable()->getRow(effectiveAtomicNum)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if atom.atomic_number() == 1 {
        if atom.formal_charge() > 1 || atom.formal_charge() < -1 {
            return Ok(true);
        }
    } else if atom.formal_charge() > atom.atomic_number() as i8
        || periodic_table_row(atom.atomic_number()) != periodic_table_row(effective_atomic_num)
    {
        return Ok(true);
    }

    // RDKit✔️✔️:   bool strict = false;
    // RDKit✔️✔️:   bool checkIt = true;
    // RDKit✔️✔️:   if (calculateExplicitValence(*this, strict, checkIt) == -1 ||
    // RDKit✔️✔️:       calculateImplicitValence(*this, strict, checkIt) == -1) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::hasValenceViolation
    let strict = false;
    let check_it = true;
    let explicit = calculate_explicit_valence(atoms, bonds, adjacency, atom_id, strict, check_it)?;
    if explicit == -1 {
        return Ok(true);
    }
    let implicit =
        calculate_implicit_valence(atoms, bonds, adjacency, atom_id, explicit, strict, check_it)?;
    Ok(implicit == -1)
}

pub(crate) fn atom_has_valence_violation_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
) -> Result<bool, ValenceError> {
    atom_has_valence_violation_impl(atoms, bonds, adjacency, atom_id)
}

pub(crate) fn assign_valence_with_options_from_read_parts(
    read: MoleculeReadParts<'_>,
    model: ValenceModel,
    strict: bool,
) -> Result<ValenceAssignment, ValenceError> {
    assign_valence_with_options_from_parts(
        read.atoms(),
        read.bonds(),
        &read.topology().adjacency,
        model,
        strict,
    )
}

pub(crate) fn assign_radicals_from_read_parts(
    read: MoleculeReadParts<'_>,
) -> Result<Vec<u8>, ValenceError> {
    assign_radicals_from_parts(read.atoms(), read.bonds(), &read.topology().adjacency)
}

fn is_aromatic_atom_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_id: AtomId,
) -> Result<bool, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION isAromaticAtom
    // RDKit✔️✔️: bool isAromaticAtom(const Atom &atom) {
    // RDKit✔️✔️:   if (atom.getIsAromatic()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    let atom = atom_from_parts(atoms, atom_id)?;
    if atom.is_aromatic() {
        return Ok(true);
    }
    // RDKit✔️✔️:   if (atom.hasOwningMol()) {
    // RDKit✔️✔️:     for (const auto &bond : atom.getOwningMol().atomBonds(&atom)) {
    // RDKit✔️✔️:       if (bond->getIsAromatic() ||
    // RDKit✔️✔️:           bond->getBondType() == Bond::BondType::AROMATIC) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for bond in incident_bonds_from_parts(atoms.len(), bonds, adjacency, atom_id)? {
        if bond.is_aromatic() || bond.order() == BondOrder::Aromatic {
            return Ok(true);
        }
    }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isAromaticAtom
    Ok(false)
}

fn get_effective_atomic_num(atom: &Atom, check_value: bool) -> Result<u8, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION getEffectiveAtomicNum
    // RDKit✔️✔️: unsigned int getEffectiveAtomicNum(const Atom &atom, bool checkValue) {
    // RDKit✔️✔️:   auto effectiveAtomicNum = atom.getAtomicNum() - atom.getFormalCharge();
    let effective_atomic_num = i32::from(atom.atomic_number()) - i32::from(atom.formal_charge());
    // RDKit✔️✔️:   if (checkValue &&
    // RDKit✔️✔️:       (effectiveAtomicNum < 0 ||
    // RDKit✔️✔️:        effectiveAtomicNum >
    // RDKit✔️✔️:            static_cast<int>(PeriodicTable::getTable()->getMaxAtomicNumber()))) {
    // RDKit✔️✔️:     throw AtomValenceException("Effective atomic number out of range",
    // RDKit✔️✔️:                                atom.getIdx());
    // RDKit✔️✔️:   }
    if check_value && !(0..=118).contains(&effective_atomic_num) {
        return Err(invalid_valence_with_message(
            atom,
            "Effective atomic number out of range".to_string(),
        ));
    }
    // RDKit✔️✔️:   effectiveAtomicNum = std::clamp(
    // RDKit✔️✔️:       effectiveAtomicNum, 0,
    // RDKit✔️✔️:       static_cast<int>(PeriodicTable::getTable()->getMaxAtomicNumber()));
    // RDKit✔️✔️:   return static_cast<unsigned int>(effectiveAtomicNum);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getEffectiveAtomicNum
    // RDKit leaves unchecked callers on the clamped periodic-table range even
    // when the checked branch would have thrown.
    Ok(effective_atomic_num.clamp(0, 118) as u8)
}

fn can_be_hypervalent(atom: &Atom, effective_atomic_num: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION canBeHypervalent
    // RDKit✔️✔️: bool canBeHypervalent(const Atom &atom, unsigned int effectiveAtomicNum) {
    // RDKit✔️✔️:   return (effectiveAtomicNum > 16 &&
    // RDKit✔️✔️:           (atom.getAtomicNum() == 15 || atom.getAtomicNum() == 16)) ||
    // RDKit✔️✔️:          (effectiveAtomicNum > 34 &&
    // RDKit✔️✔️:           (atom.getAtomicNum() == 33 || atom.getAtomicNum() == 34));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION canBeHypervalent
    (effective_atomic_num > 16 && matches!(atom.atomic_number(), 15 | 16))
        || (effective_atomic_num > 34 && matches!(atom.atomic_number(), 33 | 34))
}

pub(crate) fn bond_valence_contrib(bond: &Bond, atom_id: AtomId) -> Result<f64, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION QueryBond::getValenceContrib
    // RDKit✔️✔️: double QueryBond::getValenceContrib(const Atom *atom) const {
    // RDKit✔️✔️:   if (!hasQuery() || !QueryOps::hasComplexBondTypeQuery(*getQuery())) {
    // RDKit✔️✔️:     return Bond::getValenceContrib(atom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION QueryBond::getValenceContrib
    if bond.query().is_some_and(has_complex_bond_type_query) {
        return Ok(0.0);
    }

    // BEGIN RDKIT CPP FUNCTION Bond::getValenceContrib
    // RDKit✔️✔️: double Bond::getValenceContrib(const Atom *atom) const {
    // RDKit✔️✔️:   if (atom != getBeginAtom() && atom != getEndAtom()) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    if bond.begin() != atom_id && bond.end() != atom_id {
        return Ok(0.0);
    }
    // RDKit✔️✔️:   double res;
    // RDKit✔️✔️:   if ((getBondType() == DATIVE || getBondType() == DATIVEONE) &&
    // RDKit✔️✔️:       atom->getIdx() != getEndAtomIdx()) {
    // RDKit✔️✔️:     res = 0.0;
    // RDKit✔️✔️:   } else if (getBondType() == HYDROGEN) {
    // RDKit✔️✔️:     res = 0.0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = getBondTypeAsDouble();
    // RDKit✔️✔️:   }
    let res = if matches!(bond.order(), BondOrder::Dative | BondOrder::DativeOne)
        && atom_id != bond.end()
    {
        0.0
    } else if bond.order() == BondOrder::Hydrogen {
        0.0
    } else {
        bond_type_as_double(bond.order())?
    };
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Bond::getValenceContrib
    Ok(res)
}

pub(crate) fn has_complex_bond_type_query(query: &QueryNode<BondQueryPredicate>) -> bool {
    // BEGIN RDKIT CPP FUNCTION QueryOps::hasComplexBondTypeQuery
    // RDKit✔️✔️: bool hasComplexBondTypeQueryHelper(
    // RDKit✔️✔️:     const Queries::Query<int, Bond const *, true> &qry, bool seenBondOrder) {
    // RDKit✔️✔️:   const auto df = qry.getDescription();
    // RDKit✔️✔️:   bool isBondOrder = (df == "BondOrder");
    // RDKit✔️✔️:   // is this a bond order query?
    // RDKit✔️✔️:   if (std::find(bondOrderQueryFunctions.begin(), bondOrderQueryFunctions.end(),
    // RDKit✔️✔️:                 df) != bondOrderQueryFunctions.end()) {
    // RDKit✔️✔️:     if (seenBondOrder || !isBondOrder || qry.getNegation()) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto &child :
    // RDKit✔️✔️:        boost::make_iterator_range(qry.beginChildren(), qry.endChildren())) {
    // RDKit✔️✔️:     if (hasComplexBondTypeQueryHelper(*child, seenBondOrder | isBondOrder)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (child->getDescription() == "BondOrder") {
    // RDKit✔️✔️:       seenBondOrder = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT bool hasComplexBondTypeQuery(
    // RDKit✔️✔️:     const Queries::Query<int, Bond const *, true> &qry) {
    // RDKit✔️✔️:   return hasComplexBondTypeQueryHelper(qry, false);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION QueryOps::hasComplexBondTypeQuery
    has_complex_bond_type_query_helper(query, false).0
}

pub(crate) fn has_bond_type_query(query: &QueryNode<BondQueryPredicate>) -> bool {
    // BEGIN RDKIT CPP FUNCTION QueryOps::hasBondTypeQuery
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT bool hasBondTypeQuery(
    // RDKit✔️✔️:     const Queries::Query<int, Bond const *, true> &qry) {
    // RDKit✔️✔️:   const auto df = qry.getDescription();
    // RDKit✔️✔️:   const auto dt = qry.getTypeLabel();
    // RDKit✔️✔️:   // is this a bond order query?
    // RDKit✔️✔️:   if (dt == "BondOrder" ||
    // RDKit✔️✔️:       (dt.empty() &&
    // RDKit✔️✔️:        std::find(bondOrderQueryFunctions.begin(), bondOrderQueryFunctions.end(),
    // RDKit✔️✔️:                  df) != bondOrderQueryFunctions.end())) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto &child :
    // RDKit✔️✔️:        boost::make_iterator_range(qry.beginChildren(), qry.endChildren())) {
    // RDKit✔️✔️:     if (hasBondTypeQuery(*child)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION QueryOps::hasBondTypeQuery
    if query_node_is_bond_order_function(query) {
        return true;
    }
    match query {
        QueryNode::And(children) | QueryNode::Or(children) => {
            children.iter().any(has_bond_type_query)
        }
        QueryNode::Not(child) => has_bond_type_query(child),
        QueryNode::Predicate(_) => false,
    }
}

fn has_complex_bond_type_query_helper(
    query: &QueryNode<BondQueryPredicate>,
    seen_bond_order: bool,
) -> (bool, bool) {
    if let QueryNode::Not(child) = query
        && query_node_is_bond_order_function(child)
    {
        return (true, false);
    }

    let is_bond_order = query_node_is_bond_order(query);
    if query_node_is_bond_order_function(query)
        && (seen_bond_order || !is_bond_order || matches!(query, QueryNode::Not(_)))
    {
        return (true, is_bond_order);
    }

    let mut seen_bond_order = seen_bond_order || is_bond_order;
    match query {
        QueryNode::And(children) | QueryNode::Or(children) => {
            for child in children {
                // RDKit threads the "already saw a bond-order term" state
                // across siblings so later recursive children can detect a
                // compound bond-order tree.
                let (complex, child_is_bond_order) =
                    has_complex_bond_type_query_helper(child, seen_bond_order);
                if complex {
                    return (true, is_bond_order);
                }
                if child_is_bond_order {
                    seen_bond_order = true;
                }
            }
        }
        QueryNode::Not(child) => {
            let (complex, _) = has_complex_bond_type_query_helper(child, seen_bond_order);
            if complex {
                return (true, is_bond_order);
            }
        }
        QueryNode::Predicate(_) => {}
    }

    (false, is_bond_order)
}

fn query_node_is_bond_order(query: &QueryNode<BondQueryPredicate>) -> bool {
    match query {
        QueryNode::Predicate(BondQueryPredicate::Order(_)) => true,
        QueryNode::Predicate(BondQueryPredicate::OrderIn(orders)) => orders.len() == 1,
        QueryNode::Predicate(_) | QueryNode::And(_) | QueryNode::Or(_) | QueryNode::Not(_) => false,
    }
}

fn query_node_is_bond_order_function(query: &QueryNode<BondQueryPredicate>) -> bool {
    matches!(
        query,
        QueryNode::Predicate(BondQueryPredicate::Order(_))
            | QueryNode::Predicate(BondQueryPredicate::OrderIn(_))
    )
}

pub(crate) fn bond_type_as_double(order: BondOrder) -> Result<f64, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Bond::getBondTypeAsDouble
    // RDKit✔️✔️: double Bond::getBondTypeAsDouble() const {
    // RDKit✔️✔️:   double res;
    // RDKit✔️✔️:   switch (getBondType()) {
    let value = match order {
        // RDKit✔️✔️:     case UNSPECIFIED:
        // RDKit✔️✔️:     case IONIC:
        // RDKit✔️✔️:     case ZERO:
        // RDKit✔️✔️:       res = 0;
        BondOrder::Unspecified | BondOrder::Null | BondOrder::Ionic | BondOrder::Zero => 0.0,
        // RDKit✔️✔️:     case SINGLE:
        // RDKit✔️✔️:       res = 1;
        BondOrder::Single => 1.0,
        // RDKit✔️✔️:     case DOUBLE:
        // RDKit✔️✔️:       res = 2;
        BondOrder::Double => 2.0,
        // RDKit✔️✔️:     case TRIPLE:
        // RDKit✔️✔️:       res = 3;
        BondOrder::Triple => 3.0,
        // RDKit✔️✔️:     case QUADRUPLE:
        // RDKit✔️✔️:       res = 4;
        BondOrder::Quadruple => 4.0,
        // RDKit✔️✔️:     case QUINTUPLE:
        // RDKit✔️✔️:       res = 5;
        BondOrder::Quintuple => 5.0,
        // RDKit✔️✔️:     case HEXTUPLE:
        // RDKit✔️✔️:       res = 6;
        BondOrder::Hextuple => 6.0,
        // RDKit✔️✔️:     case ONEANDAHALF:
        // RDKit✔️✔️:       res = 1.5;
        BondOrder::OneAndHalf => 1.5,
        // RDKit✔️✔️:     case TWOANDAHALF:
        // RDKit✔️✔️:       res = 2.5;
        BondOrder::TwoAndHalf => 2.5,
        // RDKit✔️✔️:     case THREEANDAHALF:
        // RDKit✔️✔️:       res = 3.5;
        BondOrder::ThreeAndHalf => 3.5,
        // RDKit✔️✔️:     case FOURANDAHALF:
        // RDKit✔️✔️:       res = 4.5;
        BondOrder::FourAndHalf => 4.5,
        // RDKit✔️✔️:     case FIVEANDAHALF:
        // RDKit✔️✔️:       res = 5.5;
        BondOrder::FiveAndHalf => 5.5,
        // RDKit✔️✔️:     case AROMATIC:
        // RDKit✔️✔️:       res = 1.5;
        BondOrder::Aromatic => 1.5,
        // RDKit✔️✔️:     case DATIVEONE:
        // RDKit✔️✔️:       res = 1.0;
        // RDKit✔️✔️:       break;  // FIX: this should probably be different
        // RDKit✔️✔️:     case DATIVE:
        // RDKit✔️✔️:       res = 1.0;
        // RDKit✔️✔️:       break;  // FIX: again probably wrong
        BondOrder::Dative | BondOrder::DativeOne => 1.0,
        // RDKit✔️✔️:     case HYDROGEN:
        // RDKit✔️✔️:       res = 0.0;
        BondOrder::Hydrogen => 0.0,
        // RDKit✔️✔️:     default:
        // RDKit✔️✔️:       UNDER_CONSTRUCTION("Bad bond type");
        BondOrder::DativeLeft
        | BondOrder::DativeRight
        | BondOrder::ThreeCenter
        | BondOrder::Other => {
            return Err(ValenceError::BadBondType {
                message: "Bad bond type".to_string(),
            });
        }
    };
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Bond::getBondTypeAsDouble
    Ok(value)
}

#[allow(dead_code)]
fn rdkit_periodic_table_atom_data_source_frame() {
    // BEGIN RDKIT CPP FUNCTION periodicTableAtomData / atomicData::atomicData
    // RDKit✔️✔️: const std::string periodicTableAtomData =
    // RDKit✔️✔️:     R"DAT(0	*	0	0	0	0	0	0	0	0	-1
    // RDKit✔️✔️: 1	H		1	0.31	0.33	1.2	1.008	1	1	1.007825032	1
    // RDKit✔️✔️: 2	He	1	0.28	0.7	1.4	4.003	2	4	4.002603254	0
    // RDKit✔️✔️: 3	Li	2	1.28	1.23	2.2	6.941	1	7	7.01600455	1	-1
    // RDKit✔️✔️: 4	Be	2	0.96	0.9	1.9	9.012	2	9	9.0121822	2
    // RDKit✔️✔️: 5	B	2	0.84	0.82	1.8	10.812	3	11	11.0093054	3
    // RDKit✔️✔️: 6	C	2	0.76	0.77	1.7	12.011	4	12	12	4
    // RDKit✔️✔️: 7	N	2	0.71	0.7	1.6	14.007	5	14	14.003074	3
    // RDKit✔️✔️: 8	O	2	0.66	0.66	1.55	15.999	6	16	15.99491462	2
    // RDKit✔️✔️: 9	F	2	0.57	0.611	1.5	18.998	7	19	18.99840322	1
    // RDKit✔️✔️: 10	Ne	2	0.58	0.7	1.54	20.18	8	20	19.99244018	0
    // RDKit✔️✔️: 11	Na	3	1.66	1.54	2.4	22.99	1	23	22.98976928	1	-1
    // RDKit✔️✔️: 12	Mg	3	1.41	1.36	2.2	24.305	2	24	23.9850417	2	-1
    // RDKit✔️✔️: 13	Al	3	1.21	1.18	2.1	26.982	3	27	26.98153863	3
    // RDKit✔️✔️: 14	Si	3	1.11	0.937	2.1	28.086	4	28	27.97692653	4
    // RDKit✔️✔️: 15	P	3	1.07	0.89	1.95	30.974	5	31	30.97376163	3	5
    // RDKit✔️✔️: 16	S	3	1.05	1.04	1.8	32.067	6	32	31.972071	2	4	6
    // RDKit✔️✔️: 17	Cl	3	1.02	0.997	1.8	35.453	7	35	34.96885268	1
    // RDKit✔️✔️: 18	Ar	3	1.06	1.74	1.88	39.948	8	40	39.96238312	0
    // RDKit✔️✔️: 19	K	4	2.03	2.03	2.8	39.098	1	39	38.96370668	1	-1
    // RDKit✔️✔️: 20	Ca	4	1.76	1.74	2.4	40.078	2	40	39.96259098	2	-1
    // RDKit✔️✔️: 21	Sc	4	1.70	1.44	2.3	44.956	3	45	44.9559119	-1
    // RDKit✔️✔️: 22	Ti	4	1.60	1.32	2.15	47.867	4	48	47.9479463	-1
    // RDKit✔️✔️: 23	V	4	1.52	1.22	2.05	50.944	5	51	50.9439595	-1
    // RDKit✔️✔️: 24	Cr	4	1.39	1.18	2.05	51.996	6	52	51.9405075	-1
    // RDKit✔️✔️: 25	Mn	4	1.39	1.17	2.05	54.938	7	55	54.9380451	-1)DAT"
    // RDKit✔️✔️:     R"DAT(
    // RDKit✔️✔️: 26	Fe	4	1.32	1.17	2.05	55.845	8	56	55.9349375	-1
    // RDKit✔️✔️: 27	Co	4	1.26	1.16	2.0	58.933	9	59	58.933195	-1
    // RDKit✔️✔️: 28	Ni	4	1.24	1.15	2.0	58.693	10	58	57.9353429	-1
    // RDKit✔️✔️: 29	Cu	4	1.32	1.17	2.0	63.546	11	63	62.9295975	-1
    // RDKit✔️✔️: 30	Zn	4	1.22	1.25	2.1	65.39	2	64	63.9291422	-1
    // RDKit✔️✔️: 31	Ga	4	1.22	1.26	2.1	69.723	3	69	68.9255736	3
    // RDKit✔️✔️: 32	Ge	4	1.20	1.188	2.1	72.61	4	74	73.9211778	4
    // RDKit✔️✔️: 33	As	4	1.19	1.2	2.05	74.922	5	75	74.9215965	3	5
    // RDKit✔️✔️: 34	Se	4	1.20	1.17	1.9	78.96	6	80	79.9165213	2	4	6
    // RDKit✔️✔️: 35	Br	4	1.20	1.167	1.9	79.904	7	79	78.9183371	1
    // RDKit✔️✔️: 36	Kr	4	1.16	1.91	2.02	83.8	8	84	83.911507	0
    // RDKit✔️✔️: 37	Rb	5	2.20	2.16	2.9	85.468	1	85	84.91178974	1	-1
    // RDKit✔️✔️: 38	Sr	5	1.95	1.91	2.55	87.62	2	88	87.9056121	2	-1
    // RDKit✔️✔️: 39	Y	5	1.90	1.62	2.4	88.906	3	89	88.9058483	-1
    // RDKit✔️✔️: 40	Zr	5	1.75	1.45	2.3	91.224	4	90	89.9047044	-1
    // RDKit✔️✔️: 41	Nb	5	1.64	1.34	2.15	92.906	5	93	92.9063781	-1
    // RDKit✔️✔️: 42	Mo	5	1.54	1.3	2.1	95.94	6	98	97.9054082	-1
    // RDKit✔️✔️: 43	Tc	5	1.47	1.27	2.05	98	7	97	96.906365	-1
    // RDKit✔️✔️: 44	Ru	5	1.46	1.25	2.05	101.07	8	102	101.9043493	-1
    // RDKit✔️✔️: 45	Rh	5	1.42	1.25	2.0	102.906	9	103	102.905504	-1
    // RDKit✔️✔️: 46	Pd	5	1.39	1.28	2.05	106.42	10	106	105.903486	-1
    // RDKit✔️✔️: 47	Ag	5	1.45	1.34	2.1	107.868	11	107	106.905097	-1
    // RDKit✔️✔️: 48	Cd	5	1.44	1.48	2.2	112.412	2	114	113.9033585	-1
    // RDKit✔️✔️: 49	In	5	1.42	1.44	2.2	114.818	3	115	114.903878	3
    // RDKit✔️✔️: 50	Sn	5	1.39	1.385	2.25	118.711	4	120	119.9021947	2	4)DAT"
    // RDKit✔️✔️:     R"DAT(
    // RDKit✔️✔️: 51	Sb	5	1.39	1.4	2.2	121.76	5	121	120.9038157	3	5
    // RDKit✔️✔️: 52	Te	5	1.38	1.378	2.1	127.6	6	130	129.9062244	2	4	6
    // RDKit✔️✔️: 53	I	5	1.39	1.387	2.1	126.904	7	127	126.904473	1	3	5
    // RDKit✔️✔️: 54	Xe	5	1.40	1.98	2.16	131.29	8	132	131.9041535	0	2	4	6
    // RDKit✔️✔️: 55	Cs	6	2.44	2.35	3.0	132.905	1	133	132.9054519	1
    // RDKit✔️✔️: 56	Ba	6	2.15	1.98	2.7	137.328	2	138	137.9052472	2	-1
    // RDKit✔️✔️: 57	La	6	2.07	1.69	2.5	138.906	3	139	138.9063533	-1
    // RDKit✔️✔️: 58	Ce	6	2.04	1.83	2.48	140.116	4	140	139.9054387	-1
    // RDKit✔️✔️: 59	Pr	6	2.03	1.82	2.47	140.908	3	141	140.9076528	-1
    // RDKit✔️✔️: 60	Nd	6	2.01	1.81	2.45	144.24	4	142	141.9077233	-1
    // RDKit✔️✔️: 61	Pm	6	1.99	1.8	2.43	145	5	145	144.912749	-1
    // RDKit✔️✔️: 62	Sm	6	1.98	1.8	2.42	150.36	6	152	151.9197324	-1
    // RDKit✔️✔️: 63	Eu	6	1.98	1.99	2.4	151.964	7	153	152.9212303	-1
    // RDKit✔️✔️: 64	Gd	6	1.96	1.79	2.38	157.25	8	158	157.9241039	-1
    // RDKit✔️✔️: 65	Tb	6	1.94	1.76	2.37	158.925	9	159	158.9253468	-1
    // RDKit✔️✔️: 66	Dy	6	1.92	1.75	2.35	162.5	10	164	163.9291748	-1
    // RDKit✔️✔️: 67	Ho	6	1.92	1.74	2.33	164.93	11	165	164.9303221	-1
    // RDKit✔️✔️: 68	Er	6	1.89	1.73	2.32	167.26	12	166	165.9302931	-1
    // RDKit✔️✔️: 69	Tm	6	1.90	1.72	2.3	168.934	13	169	168.9342133	-1
    // RDKit✔️✔️: 70	Yb	6	1.87	1.94	2.28	173.04	14	174	173.9388621	-1
    // RDKit✔️✔️: 71	Lu	6	1.87	1.72	2.27	174.967	15	175	174.9407718	-1
    // RDKit✔️✔️: 72	Hf	6	1.75	1.44	2.25	178.49	4	180	179.94655	-1
    // RDKit✔️✔️: 73	Ta	6	1.70	1.34	2.2	180.948	5	181	180.9479958	-1
    // RDKit✔️✔️: 74	W	6	1.62	1.3	2.1	183.84	6	184	183.9509312	-1
    // RDKit✔️✔️: 75	Re	6	1.51	1.28	2.05	186.207	7	187	186.9557531	-1)DAT"
    // RDKit✔️✔️:     R"DAT(
    // RDKit✔️✔️: 76	Os	6	1.44	1.26	2.0	190.23	8	192	191.9614807	-1
    // RDKit✔️✔️: 77	Ir	6	1.41	1.27	2.0	192.217	9	193	192.9629264	-1
    // RDKit✔️✔️: 78	Pt	6	1.36	1.3	2.05	195.078	10	195	194.9647911	-1
    // RDKit✔️✔️: 79	Au	6	1.36	1.34	2.1	196.967	11	197	196.9665687	-1
    // RDKit✔️✔️: 80	Hg	6	1.32	1.49	2.05	200.59	2	202	201.970643	-1
    // RDKit✔️✔️: 81	Tl	6	1.45	1.48	2.2	204.383	3	205	204.9744275	-1
    // RDKit✔️✔️: 82	Pb	6	1.46	1.48	2.3	207.2	4	208	207.9766521	2	4
    // RDKit✔️✔️: 83	Bi	6	1.48	1.45	2.3	208.98	5	209	208.9803987	3	5
    // RDKit✔️✔️: 84	Po	6	1.40	1.46	2.0	209	6	209	208.9824304	2	4	6
    // RDKit✔️✔️: 85	At	6	1.50	1.45	2.0	210	7	210	209.987148	1	3	5)DAT"
    // RDKit✔️✔️:     // the values for Ra and Rn are from:
    // RDKit✔️✔️:     // https://www.ciaaw.org/atomic-masses.htm
    // RDKit✔️✔️:     // their most common isotopes are from
    // RDKit✔️✔️:     // https://www.ciaaw.org/radioactive-elements.htm
    // RDKit✔️✔️:     R"DAT(
    // RDKit✔️✔️: 86	Rn	6	1.50	2.4	2.0	222	8	222	222.0175706	0
    // RDKit✔️✔️: 87	Fr	7	2.6	2	2.0	223	1	223	223.0197359	1
    // RDKit✔️✔️: 88	Ra	7	2.2	1.9	2.0	226	2	226	226.0254026	2	-1
    // RDKit✔️✔️: 89	Ac	7	2.15	1.88	2.0	227	3	227	227.0277521	-1
    // RDKit✔️✔️: 90	Th	7	2.06	1.79	2.4	232.038	4	232	232.0380553	-1
    // RDKit✔️✔️: 91	Pa	7	2.00	1.61	2.0	231.036	3	231	231.035884	-1
    // RDKit✔️✔️: 92	U	7	1.96	1.58	2.3	238.029	4	238	238.0507882	-1
    // RDKit✔️✔️: 93	Np	7	1.90	1.55	2.0	237	5	236	236.04657	-1
    // RDKit✔️✔️: 94	Pu	7	1.87	1.53	2.0	244	6	238	238.0495599	-1
    // RDKit✔️✔️: 95	Am	7	1.80	1.07	2.0	243	7	241	241.0568291	-1
    // RDKit✔️✔️: 96	Cm	7	1.69	0	2.0	247	8	243	243.0613891	-1
    // RDKit✔️✔️: 97	Bk	7	1.9	0	2.0	247	9	247	247.070307	-1
    // RDKit✔️✔️: 98	Cf	7	1.9	0	2.0	251	10	249	249.0748535	-1
    // RDKit✔️✔️: 99	Es	7	1.9	0	2.0	252	11	252	252.08298	-1
    // RDKit✔️✔️: 100	Fm	7	1.9	0	2.0	257	12	257	257.095105	-1)DAT"
    // RDKit✔️✔️:     R"DAT(
    // RDKit✔️✔️: 101	Md	7	1.9	0	2.0	258	13	258	258.098431	-1
    // RDKit✔️✔️: 102	No	7	1.9	0	2.0	259	14	259	259.10103	-1
    // RDKit✔️✔️: 103	Lr	7	1.9	0	2.0	262	15	262	262.10963	-1
    // RDKit✔️✔️: 104	Rf	7	1.9	0	2.0	267	2	267	267.12153	-1
    // RDKit✔️✔️: 105	Db	7	1.9	0	2.0	268	2		268	268.12545	-1
    // RDKit✔️✔️: 106	Sg	7	1.9	0	2.0	269	2	271	271.13347	-1
    // RDKit✔️✔️: 107	Bh	7	1.9	0	2.0	270	2	270	270.13362	-1
    // RDKit✔️✔️: 108	Hs	7	1.9	0	2.0	269	2	269	269.13406	-1
    // RDKit✔️✔️: 109	Mt	7	1.9	0	2.0	278	2	278	278.15481	-1
    // RDKit✔️✔️: 110	Ds	7	1.9	0	2.0	281	2	281	281.16206	-1
    // RDKit✔️✔️: 111	Rg	7	1.9	0	2.0	281	2	281	281.16537	-1
    // RDKit✔️✔️: 112	Cn	7	1.9	0	2.0	285	2	285	285.17411	-1)DAT"
    // RDKit✔️✔️:     // Values for the below elements are from the BODR.
    // RDKit✔️✔️:     // The Blue Obelisk Data Repository has no VdW radii for those:
    // RDKit✔️✔️:     // Ds Rg Cn Uut Fl Mc Uup Lv Ts Og; we use 2.0 for all of them
    // RDKit✔️✔️:     // ---
    // RDKit✔️✔️:     // added from BODR 30.10.2016
    // RDKit✔️✔️:     // atomic mass data from NIST
    // RDKit✔️✔️:     // we leave Uut and Uup in here for backwards
    // RDKit✔️✔️:     // compatibility. Nh and Mc (the entries appearing first
    // RDKit✔️✔️:     // for a particular atomic number) will be the values returned
    // RDKit✔️✔️:     // when looking an atomic symbol up using atomic number.
    // RDKit✔️✔️:     // In the event that symbols need to be corrected in the future,
    // RDKit✔️✔️:     // follow this same pattern to ensure backwards compatibility.
    // RDKit✔️✔️:     R"DAT(
    // RDKit✔️✔️: 113	Nh	7	1.36	0	2.0	284	2	284	284.17873	-1
    // RDKit✔️✔️: 113	Uut	7	1.36	0	2.0	284	2	284	284.17873	-1
    // RDKit✔️✔️: 114	Fl	7	1.43	0	2.0	289	2	289	289.19042	-1
    // RDKit✔️✔️: 115	Mc	7	1.62	0	2.0	288	2	288	288.19274	-1
    // RDKit✔️✔️: 115	Uup	7	1.62	0	2.0	288	2	288	288.19274	-1
    // RDKit✔️✔️: 116	Lv	7	1.75	0	2.0	293	2	293	293.20449	-1
    // RDKit✔️✔️: 117	Ts	7	1.65	0	2.0	292	2	292	292.20746	-1
    // RDKit✔️✔️: 118	Og	7	1.57	0	2.0	294	2	294	294.21392	-1
    // RDKit✔️✔️: )DAT";
    // RDKit✔️✔️: atomicData::atomicData(const std::string &dataLine) {
    // RDKit✔️✔️:   boost::char_separator<char> spaceSep(" \t");
    // RDKit✔️✔️:   tokenizer tokens(dataLine, spaceSep);
    // RDKit✔️✔️:   tokenizer::iterator token = tokens.begin();
    // RDKit✔️✔️:   std::istringstream istr;
    // RDKit✔️✔️:   istr.imbue(std::locale("C"));
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> anum;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️:   name = elementNames[anum];
    // RDKit✔️✔️:   symb = *token;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> row;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> rCov;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> rB0;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> rVdw;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> mass;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> nVal;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> commonIsotope;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> commonIsotopeMass;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️:   valence.clear();
    // RDKit✔️✔️:   while (token != tokens.end()) {
    // RDKit✔️✔️:     istr.clear();
    // RDKit✔️✔️:     istr.str(*token);
    // RDKit✔️✔️:     int tval;
    // RDKit✔️✔️:     istr >> tval;
    // RDKit✔️✔️:     valence.push_back(tval);
    // RDKit✔️✔️:     ++token;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION periodicTableAtomData / atomicData::atomicData
}

#[derive(Debug, Clone, Copy)]
struct RdkitPeriodicTableEntry {
    symbol: &'static str,
    row: u8,
    n_outer: i32,
    valence_list: &'static [i32],
}

const RDKIT_PERIODIC_TABLE: [RdkitPeriodicTableEntry; 119] = [
    RdkitPeriodicTableEntry {
        symbol: "*",
        row: 0,
        n_outer: 0,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "H",
        row: 1,
        n_outer: 1,
        valence_list: &[1],
    },
    RdkitPeriodicTableEntry {
        symbol: "He",
        row: 1,
        n_outer: 2,
        valence_list: &[0],
    },
    RdkitPeriodicTableEntry {
        symbol: "Li",
        row: 2,
        n_outer: 1,
        valence_list: &[1, -1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Be",
        row: 2,
        n_outer: 2,
        valence_list: &[2],
    },
    RdkitPeriodicTableEntry {
        symbol: "B",
        row: 2,
        n_outer: 3,
        valence_list: &[3],
    },
    RdkitPeriodicTableEntry {
        symbol: "C",
        row: 2,
        n_outer: 4,
        valence_list: &[4],
    },
    RdkitPeriodicTableEntry {
        symbol: "N",
        row: 2,
        n_outer: 5,
        valence_list: &[3],
    },
    RdkitPeriodicTableEntry {
        symbol: "O",
        row: 2,
        n_outer: 6,
        valence_list: &[2],
    },
    RdkitPeriodicTableEntry {
        symbol: "F",
        row: 2,
        n_outer: 7,
        valence_list: &[1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ne",
        row: 2,
        n_outer: 8,
        valence_list: &[0],
    },
    RdkitPeriodicTableEntry {
        symbol: "Na",
        row: 3,
        n_outer: 1,
        valence_list: &[1, -1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Mg",
        row: 3,
        n_outer: 2,
        valence_list: &[2, -1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Al",
        row: 3,
        n_outer: 3,
        valence_list: &[3],
    },
    RdkitPeriodicTableEntry {
        symbol: "Si",
        row: 3,
        n_outer: 4,
        valence_list: &[4],
    },
    RdkitPeriodicTableEntry {
        symbol: "P",
        row: 3,
        n_outer: 5,
        valence_list: &[3, 5],
    },
    RdkitPeriodicTableEntry {
        symbol: "S",
        row: 3,
        n_outer: 6,
        valence_list: &[2, 4, 6],
    },
    RdkitPeriodicTableEntry {
        symbol: "Cl",
        row: 3,
        n_outer: 7,
        valence_list: &[1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ar",
        row: 3,
        n_outer: 8,
        valence_list: &[0],
    },
    RdkitPeriodicTableEntry {
        symbol: "K",
        row: 4,
        n_outer: 1,
        valence_list: &[1, -1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ca",
        row: 4,
        n_outer: 2,
        valence_list: &[2, -1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Sc",
        row: 4,
        n_outer: 3,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ti",
        row: 4,
        n_outer: 4,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "V",
        row: 4,
        n_outer: 5,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Cr",
        row: 4,
        n_outer: 6,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Mn",
        row: 4,
        n_outer: 7,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Fe",
        row: 4,
        n_outer: 8,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Co",
        row: 4,
        n_outer: 9,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ni",
        row: 4,
        n_outer: 10,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Cu",
        row: 4,
        n_outer: 11,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Zn",
        row: 4,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ga",
        row: 4,
        n_outer: 3,
        valence_list: &[3],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ge",
        row: 4,
        n_outer: 4,
        valence_list: &[4],
    },
    RdkitPeriodicTableEntry {
        symbol: "As",
        row: 4,
        n_outer: 5,
        valence_list: &[3, 5],
    },
    RdkitPeriodicTableEntry {
        symbol: "Se",
        row: 4,
        n_outer: 6,
        valence_list: &[2, 4, 6],
    },
    RdkitPeriodicTableEntry {
        symbol: "Br",
        row: 4,
        n_outer: 7,
        valence_list: &[1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Kr",
        row: 4,
        n_outer: 8,
        valence_list: &[0],
    },
    RdkitPeriodicTableEntry {
        symbol: "Rb",
        row: 5,
        n_outer: 1,
        valence_list: &[1, -1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Sr",
        row: 5,
        n_outer: 2,
        valence_list: &[2, -1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Y",
        row: 5,
        n_outer: 3,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Zr",
        row: 5,
        n_outer: 4,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Nb",
        row: 5,
        n_outer: 5,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Mo",
        row: 5,
        n_outer: 6,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Tc",
        row: 5,
        n_outer: 7,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ru",
        row: 5,
        n_outer: 8,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Rh",
        row: 5,
        n_outer: 9,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Pd",
        row: 5,
        n_outer: 10,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ag",
        row: 5,
        n_outer: 11,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Cd",
        row: 5,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "In",
        row: 5,
        n_outer: 3,
        valence_list: &[3],
    },
    RdkitPeriodicTableEntry {
        symbol: "Sn",
        row: 5,
        n_outer: 4,
        valence_list: &[2, 4],
    },
    RdkitPeriodicTableEntry {
        symbol: "Sb",
        row: 5,
        n_outer: 5,
        valence_list: &[3, 5],
    },
    RdkitPeriodicTableEntry {
        symbol: "Te",
        row: 5,
        n_outer: 6,
        valence_list: &[2, 4, 6],
    },
    RdkitPeriodicTableEntry {
        symbol: "I",
        row: 5,
        n_outer: 7,
        valence_list: &[1, 3, 5],
    },
    RdkitPeriodicTableEntry {
        symbol: "Xe",
        row: 5,
        n_outer: 8,
        valence_list: &[0, 2, 4, 6],
    },
    RdkitPeriodicTableEntry {
        symbol: "Cs",
        row: 6,
        n_outer: 1,
        valence_list: &[1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ba",
        row: 6,
        n_outer: 2,
        valence_list: &[2, -1],
    },
    RdkitPeriodicTableEntry {
        symbol: "La",
        row: 6,
        n_outer: 3,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ce",
        row: 6,
        n_outer: 4,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Pr",
        row: 6,
        n_outer: 3,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Nd",
        row: 6,
        n_outer: 4,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Pm",
        row: 6,
        n_outer: 5,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Sm",
        row: 6,
        n_outer: 6,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Eu",
        row: 6,
        n_outer: 7,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Gd",
        row: 6,
        n_outer: 8,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Tb",
        row: 6,
        n_outer: 9,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Dy",
        row: 6,
        n_outer: 10,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ho",
        row: 6,
        n_outer: 11,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Er",
        row: 6,
        n_outer: 12,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Tm",
        row: 6,
        n_outer: 13,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Yb",
        row: 6,
        n_outer: 14,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Lu",
        row: 6,
        n_outer: 15,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Hf",
        row: 6,
        n_outer: 4,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ta",
        row: 6,
        n_outer: 5,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "W",
        row: 6,
        n_outer: 6,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Re",
        row: 6,
        n_outer: 7,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Os",
        row: 6,
        n_outer: 8,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ir",
        row: 6,
        n_outer: 9,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Pt",
        row: 6,
        n_outer: 10,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Au",
        row: 6,
        n_outer: 11,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Hg",
        row: 6,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Tl",
        row: 6,
        n_outer: 3,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Pb",
        row: 6,
        n_outer: 4,
        valence_list: &[2, 4],
    },
    RdkitPeriodicTableEntry {
        symbol: "Bi",
        row: 6,
        n_outer: 5,
        valence_list: &[3, 5],
    },
    RdkitPeriodicTableEntry {
        symbol: "Po",
        row: 6,
        n_outer: 6,
        valence_list: &[2, 4, 6],
    },
    RdkitPeriodicTableEntry {
        symbol: "At",
        row: 6,
        n_outer: 7,
        valence_list: &[1, 3, 5],
    },
    RdkitPeriodicTableEntry {
        symbol: "Rn",
        row: 6,
        n_outer: 8,
        valence_list: &[0],
    },
    RdkitPeriodicTableEntry {
        symbol: "Fr",
        row: 7,
        n_outer: 1,
        valence_list: &[1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ra",
        row: 7,
        n_outer: 2,
        valence_list: &[2, -1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ac",
        row: 7,
        n_outer: 3,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Th",
        row: 7,
        n_outer: 4,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Pa",
        row: 7,
        n_outer: 3,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "U",
        row: 7,
        n_outer: 4,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Np",
        row: 7,
        n_outer: 5,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Pu",
        row: 7,
        n_outer: 6,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Am",
        row: 7,
        n_outer: 7,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Cm",
        row: 7,
        n_outer: 8,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Bk",
        row: 7,
        n_outer: 9,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Cf",
        row: 7,
        n_outer: 10,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Es",
        row: 7,
        n_outer: 11,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Fm",
        row: 7,
        n_outer: 12,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Md",
        row: 7,
        n_outer: 13,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "No",
        row: 7,
        n_outer: 14,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Lr",
        row: 7,
        n_outer: 15,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Rf",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Db",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Sg",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Bh",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Hs",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Mt",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ds",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Rg",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Cn",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Nh",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Fl",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Mc",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Lv",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Ts",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
    RdkitPeriodicTableEntry {
        symbol: "Og",
        row: 7,
        n_outer: 2,
        valence_list: &[-1],
    },
];

fn rdkit_periodic_table_entry(atomic_number: u8) -> Option<&'static RdkitPeriodicTableEntry> {
    RDKIT_PERIODIC_TABLE.get(usize::from(atomic_number))
}

const RDKit_ATOMIC_WEIGHTS: [f64; 119] = [
    0.0, 1.008, 4.003, 6.941, 9.012, 10.812, 12.011, 14.007, 15.999, 18.998, 20.18, 22.99, 24.305,
    26.982, 28.086, 30.974, 32.067, 35.453, 39.948, 39.098, 40.078, 44.956, 47.867, 50.944, 51.996,
    54.938, 55.845, 58.933, 58.693, 63.546, 65.39, 69.723, 72.61, 74.922, 78.96, 79.904, 83.8,
    85.468, 87.62, 88.906, 91.224, 92.906, 95.94, 98.0, 101.07, 102.906, 106.42, 107.868, 112.412,
    114.818, 118.711, 121.76, 127.6, 126.904, 131.29, 132.905, 137.328, 138.906, 140.116, 140.908,
    144.24, 145.0, 150.36, 151.964, 157.25, 158.925, 162.5, 164.93, 167.26, 168.934, 173.04,
    174.967, 178.49, 180.948, 183.84, 186.207, 190.23, 192.217, 195.078, 196.967, 200.59, 204.383,
    207.2, 208.98, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.038, 231.036, 238.029, 237.0,
    244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 262.0, 267.0, 268.0, 269.0,
    270.0, 269.0, 278.0, 281.0, 281.0, 285.0, 284.179, 289.190, 288.193, 293.204, 292.207, 294.214,
];

#[must_use]
pub(crate) fn rdkit_atomic_mass(atomic_number: u8, isotope: Option<u16>) -> f64 {
    // BEGIN RDKIT CPP FUNCTION Atom::getMass / PeriodicTable::getAtomicWeight / PeriodicTable::getMassForIsotope
    // RDKit✔️✔️: double Atom::getMass() const {
    // RDKit✔️✔️:   if (d_isotope) {
    // RDKit✔️✔️:     return PeriodicTable::getTable()->getMassForIsotope(d_atomicNum, d_isotope);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return PeriodicTable::getTable()->getAtomicWeight(d_atomicNum);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double getAtomicWeight(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   double mass = byanum[atomicNumber].Mass();
    // RDKit✔️✔️:   return mass;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double getMassForIsotope(UINT atomicNumber, UINT isotope) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   const std::map<unsigned int, std::pair<double, double>> &m =
    // RDKit✔️✔️:       byanum[atomicNumber].d_isotopeInfoMap;
    // RDKit✔️✔️:   std::map<unsigned int, std::pair<double, double>>::const_iterator item =
    // RDKit✔️✔️:       m.find(isotope);
    // RDKit✔️✔️:   if (item == m.end()) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return item->second.first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getMass / PeriodicTable::getAtomicWeight / PeriodicTable::getMassForIsotope
    let atomic_weight = RDKit_ATOMIC_WEIGHTS
        .get(usize::from(atomic_number))
        .copied()
        .unwrap_or(0.0);
    if let Some(isotope) = isotope {
        RDKIT_ISOTOPE_MASSES
            .binary_search_by_key(&(atomic_number, isotope), |&(key, _)| key)
            .map(|index| RDKIT_ISOTOPE_MASSES[index].1)
            .unwrap_or(0.0)
    } else {
        atomic_weight
    }
}

#[must_use]
pub fn rdkit_atomic_number_from_symbol(symbol: &str) -> Option<u8> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/PeriodicTable.h :: getAtomicNumber
    // RDKit✔️🔝: int getAtomicNumber(const std::string &elementSymbol) const {
    // RDKit✔️🔝:   // this little optimization actually makes a measurable difference
    // RDKit✔️🔝:   // in molecule-construction time
    // RDKit✔️🔝:   int anum = -1;
    // RDKit✔️🔝:   if (elementSymbol == "C") {
    // RDKit✔️🔝:     anum = 6;
    // RDKit✔️🔝:   } else if (elementSymbol == "N") {
    // RDKit✔️🔝:     anum = 7;
    // RDKit✔️🔝:   } else if (elementSymbol == "O") {
    // RDKit✔️🔝:     anum = 8;
    // RDKit✔️🔝:   } else {
    // RDKit✔️🔝:     STR_UINT_MAP::const_iterator iter = byname.find(elementSymbol);
    // RDKit✔️🔝:     if (iter != byname.end()) {
    // RDKit✔️🔝:       anum = iter->second;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   POSTCONDITION(anum > -1, "Element '" + elementSymbol + "' not found");
    // RDKit✔️🔝:   return anum;
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/PeriodicTable.h :: getAtomicNumber
    //
    // BEGIN RDKIT CPP TYPEDEF third_party/rdkit/Code/RDGeneral/types.h :: STR_UINT_MAP
    // RDKit✔️🔝: typedef std::map<std::string, UINT> STR_UINT_MAP;
    // END RDKIT CPP TYPEDEF third_party/rdkit/Code/RDGeneral/types.h :: STR_UINT_MAP
    //
    // The symbol set is static and tiny. A direct string match preserves
    // RDKit's byname coverage, including legacy aliases retained in
    // atomic_data, while avoiding RDKit's O(log n) std::map traversal and
    // runtime table storage.
    match symbol {
        "*" => Some(0),
        "H" => Some(1),
        "He" => Some(2),
        "Li" => Some(3),
        "Be" => Some(4),
        "B" => Some(5),
        "C" => Some(6),
        "N" => Some(7),
        "O" => Some(8),
        "F" => Some(9),
        "Ne" => Some(10),
        "Na" => Some(11),
        "Mg" => Some(12),
        "Al" => Some(13),
        "Si" => Some(14),
        "P" => Some(15),
        "S" => Some(16),
        "Cl" => Some(17),
        "Ar" => Some(18),
        "K" => Some(19),
        "Ca" => Some(20),
        "Sc" => Some(21),
        "Ti" => Some(22),
        "V" => Some(23),
        "Cr" => Some(24),
        "Mn" => Some(25),
        "Fe" => Some(26),
        "Co" => Some(27),
        "Ni" => Some(28),
        "Cu" => Some(29),
        "Zn" => Some(30),
        "Ga" => Some(31),
        "Ge" => Some(32),
        "As" => Some(33),
        "Se" => Some(34),
        "Br" => Some(35),
        "Kr" => Some(36),
        "Rb" => Some(37),
        "Sr" => Some(38),
        "Y" => Some(39),
        "Zr" => Some(40),
        "Nb" => Some(41),
        "Mo" => Some(42),
        "Tc" => Some(43),
        "Ru" => Some(44),
        "Rh" => Some(45),
        "Pd" => Some(46),
        "Ag" => Some(47),
        "Cd" => Some(48),
        "In" => Some(49),
        "Sn" => Some(50),
        "Sb" => Some(51),
        "Te" => Some(52),
        "I" => Some(53),
        "Xe" => Some(54),
        "Cs" => Some(55),
        "Ba" => Some(56),
        "La" => Some(57),
        "Ce" => Some(58),
        "Pr" => Some(59),
        "Nd" => Some(60),
        "Pm" => Some(61),
        "Sm" => Some(62),
        "Eu" => Some(63),
        "Gd" => Some(64),
        "Tb" => Some(65),
        "Dy" => Some(66),
        "Ho" => Some(67),
        "Er" => Some(68),
        "Tm" => Some(69),
        "Yb" => Some(70),
        "Lu" => Some(71),
        "Hf" => Some(72),
        "Ta" => Some(73),
        "W" => Some(74),
        "Re" => Some(75),
        "Os" => Some(76),
        "Ir" => Some(77),
        "Pt" => Some(78),
        "Au" => Some(79),
        "Hg" => Some(80),
        "Tl" => Some(81),
        "Pb" => Some(82),
        "Bi" => Some(83),
        "Po" => Some(84),
        "At" => Some(85),
        "Rn" => Some(86),
        "Fr" => Some(87),
        "Ra" => Some(88),
        "Ac" => Some(89),
        "Th" => Some(90),
        "Pa" => Some(91),
        "U" => Some(92),
        "Np" => Some(93),
        "Pu" => Some(94),
        "Am" => Some(95),
        "Cm" => Some(96),
        "Bk" => Some(97),
        "Cf" => Some(98),
        "Es" => Some(99),
        "Fm" => Some(100),
        "Md" => Some(101),
        "No" => Some(102),
        "Lr" => Some(103),
        "Rf" => Some(104),
        "Db" => Some(105),
        "Sg" => Some(106),
        "Bh" => Some(107),
        "Hs" => Some(108),
        "Mt" => Some(109),
        "Ds" => Some(110),
        "Rg" => Some(111),
        "Cn" => Some(112),
        "Nh" | "Uut" => Some(113),
        "Fl" => Some(114),
        "Mc" | "Uup" => Some(115),
        "Lv" => Some(116),
        "Ts" => Some(117),
        "Og" => Some(118),
        _ => None,
    }
}

pub(crate) fn rdkit_default_valence(atomic_number: u8) -> Result<i32, ValenceError> {
    Ok(required_valence_list(atomic_number)?[0])
}

fn periodic_table_row(atomic_number: u8) -> Option<u8> {
    // BEGIN RDKIT CPP FUNCTION atomicData::Row / PeriodicTable::getRow
    // RDKit✔️✔️: unsigned int Row() const { return row; }
    // END RDKIT CPP FUNCTION atomicData::Row / PeriodicTable::getRow
    rdkit_periodic_table_entry(atomic_number).map(|entry| entry.row)
}

pub(crate) fn periodic_table_outer_electrons(atomic_number: u8) -> Result<i32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION atomicData::NumOuterShellElec / PeriodicTable::getNouterElecs
    // RDKit✔️✔️: int NumOuterShellElec() const { return nVal; }
    // RDKit✔️✔️: int getNouterElecs(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].NumOuterShellElec();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT extern const std::string periodicTableAtomData;
    // END RDKIT CPP FUNCTION atomicData::NumOuterShellElec / PeriodicTable::getNouterElecs
    rdkit_periodic_table_entry(atomic_number)
        .map(|entry| entry.n_outer)
        .ok_or(ValenceError::UnsupportedBranch {
            reason: "PeriodicTable outer electrons atomic number out of range",
        })
}

pub(crate) fn rdkit_most_common_isotope(atomic_number: u8) -> Result<i64, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getMostCommonIsotope / atomicData::MostCommonIsotope
    // RDKit✔️✔️: int MostCommonIsotope() const { return commonIsotope; }
    // RDKit✔️✔️: int getMostCommonIsotope(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].MostCommonIsotope();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: atomicData::atomicData(const std::string &dataLine) {
    // RDKit✔️✔️:   // most common isotope
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> commonIsotope;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT extern const std::string periodicTableAtomData;
    // END RDKIT CPP FUNCTION PeriodicTable::getMostCommonIsotope / atomicData::MostCommonIsotope
    RDKIT_MOST_COMMON_ISOTOPES
        .get(usize::from(atomic_number))
        .map(|isotope| i64::from(*isotope))
        .ok_or(ValenceError::UnsupportedBranch {
            reason: "PeriodicTable most-common isotope atomic number out of range",
        })
}

pub(crate) fn rdkit_most_common_isotope_mass(atomic_number: u8) -> Result<f64, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getMostCommonIsotopeMass / atomicData::MostCommonIsotopeMass
    // RDKit✔️✔️: double MostCommonIsotopeMass() const { return commonIsotopeMass; }
    // RDKit✔️✔️: double getMostCommonIsotopeMass(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].MostCommonIsotopeMass();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: atomicData::atomicData(const std::string &dataLine) {
    // RDKit✔️✔️:   // most common isotopic mass
    // RDKit✔️✔️:   istr.clear();
    // RDKit✔️✔️:   istr.str(*token);
    // RDKit✔️✔️:   istr >> commonIsotopeMass;
    // RDKit✔️✔️:   ++token;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT extern const std::string periodicTableAtomData;
    // END RDKIT CPP FUNCTION PeriodicTable::getMostCommonIsotopeMass / atomicData::MostCommonIsotopeMass
    RDKIT_MOST_COMMON_ISOTOPE_MASSES
        .get(usize::from(atomic_number))
        .copied()
        .ok_or(ValenceError::UnsupportedBranch {
            reason: "PeriodicTable most-common isotope mass atomic number out of range",
        })
}

pub(crate) fn periodic_table_more_electronegative(
    atomic_number_1: u8,
    atomic_number_2: u8,
) -> Result<bool, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::moreElectroNegative
    // RDKit✔️✔️: bool moreElectroNegative(UINT anum1, UINT anum2) const {
    // RDKit✔️✔️:   PRECONDITION(anum1 < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   PRECONDITION(anum2 < byanum.size(), "Atomic number not found");
    let outer_1 = periodic_table_outer_electrons(atomic_number_1)?;
    let outer_2 = periodic_table_outer_electrons(atomic_number_2)?;
    // RDKit✔️✔️:   // FIX: the atomic_data needs to have real electronegativity values
    // RDKit✔️✔️:   UINT ne1 = getNouterElecs(anum1);
    // RDKit✔️✔️:   UINT ne2 = getNouterElecs(anum2);
    // RDKit✔️✔️:   if (ne1 > ne2) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    if outer_1 > outer_2 {
        return Ok(true);
    }
    // RDKit✔️✔️:   if (ne1 == ne2) {
    // RDKit✔️✔️:     if (anum1 < anum2) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if outer_1 == outer_2 && atomic_number_1 < atomic_number_2 {
        return Ok(true);
    }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::moreElectroNegative
    Ok(false)
}

fn required_valence_list(atomic_number: u8) -> Result<&'static [i32], ValenceError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getValenceList
    // RDKit✔️✔️: //! returns a vector of all stable valences. For atoms where
    // RDKit✔️✔️: //! we really don't have any idea what a reasonable maximum
    // RDKit✔️✔️: //! valence is (like transition metals), the vector ends with -1
    // RDKit✔️✔️: const INT_VECT &getValenceList(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].ValenceList();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::getValenceList
    rdkit_valence_list_for_atomic_number(atomic_number)?.ok_or(ValenceError::UnsupportedBranch {
        reason: "PeriodicTable atomic number out of range",
    })
}

fn rdkit_valence_list_for_atomic_number(
    atomic_number: u8,
) -> Result<Option<&'static [i32]>, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION atomicData::ValenceList / PeriodicTable::getValenceList
    // RDKit✔️✔️: int DefaultValence() const { return valence.front(); }
    // RDKit✔️✔️: const INT_VECT &ValenceList() const { return valence; }
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT extern const std::string periodicTableAtomData;
    // END RDKIT CPP FUNCTION atomicData::ValenceList / PeriodicTable::getValenceList
    Ok(rdkit_periodic_table_entry(atomic_number).map(|entry| entry.valence_list))
}

fn invalid_valence_with_message(atom: &Atom, message: String) -> ValenceError {
    ValenceError::InvalidValence {
        atom: atom.id(),
        atomic_number: atom.atomic_number(),
        formal_charge: atom.formal_charge(),
        message,
    }
}

pub fn rdkit_element_symbol(atomic_number: u8) -> Result<&'static str, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getElementSymbol / atomicData::Symbol
    // RDKit✔️✔️: std::string getElementSymbol(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].Symbol();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: std::string Symbol() const { return symb; }
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT extern const std::string periodicTableAtomData;
    // END RDKIT CPP FUNCTION PeriodicTable::getElementSymbol / atomicData::Symbol
    rdkit_periodic_table_entry(atomic_number)
        .map(|entry| entry.symbol)
        .ok_or(ValenceError::UnsupportedBranch {
            reason: "PeriodicTable element symbol atomic number out of range",
        })
}

#[cfg(test)]
mod tests {
    use super::rdkit_element_symbol;
    use crate::{
        AtomId, AtomQueryPredicate, AtomSpec, BondOrder, BondQueryPredicate, BondSpec, Element,
        Molecule, MoleculeBuilder, QueryNode, ValenceAssignment, ValenceModel, assign_radicals,
        assign_valence, assign_valence_with_options, atom_has_valence_violation,
        rdkit_valence_list,
    };

    #[test]
    fn rdkit_valence_list_covers_periodic_table_entries() {
        assert_eq!(rdkit_valence_list(0).unwrap(), Some(&[-1][..]));
        assert_eq!(rdkit_valence_list(6).unwrap(), Some(&[4][..]));
        assert_eq!(rdkit_valence_list(16).unwrap(), Some(&[2, 4, 6][..]));
        assert_eq!(rdkit_valence_list(53).unwrap(), Some(&[1, 3, 5][..]));
        assert_eq!(rdkit_valence_list(118).unwrap(), Some(&[-1][..]));
        assert_eq!(rdkit_valence_list(119).unwrap(), None);
    }

    #[test]
    fn rdkit_element_symbol_matches_periodic_table_entries_used_in_valence_errors() {
        assert_eq!(rdkit_element_symbol(0).unwrap(), "*");
        assert_eq!(rdkit_element_symbol(6).unwrap(), "C");
        assert_eq!(rdkit_element_symbol(118).unwrap(), "Og");
        assert!(rdkit_element_symbol(119).is_err());
    }

    #[test]
    fn rdkit_atomic_number_from_symbol_covers_canonical_and_legacy_entries() {
        for atomic_number in 0..=118 {
            let symbol = rdkit_element_symbol(atomic_number).unwrap();
            assert_eq!(
                super::rdkit_atomic_number_from_symbol(symbol),
                Some(atomic_number)
            );
        }
        assert_eq!(super::rdkit_atomic_number_from_symbol("Uut"), Some(113));
        assert_eq!(super::rdkit_atomic_number_from_symbol("Uup"), Some(115));
    }

    #[test]
    fn valence_parts_atom_and_incident_bonds_follow_topology_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        assert_eq!(
            super::atom_from_parts(molecule.atoms(), AtomId::new(1))
                .unwrap()
                .atomic_number(),
            6
        );
        assert_eq!(
            super::incident_bonds_from_parts(
                molecule.num_atoms(),
                molecule.bonds(),
                &molecule.topology_block().adjacency,
                AtomId::new(1),
            )
            .unwrap()
            .map(|bond| bond.id().index())
            .collect::<Vec<_>>(),
            vec![0, 1]
        );
    }

    #[test]
    fn valence_parts_use_topology_adjacency_without_cached_derived_cache() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        assert_eq!(
            super::incident_bonds_from_parts(
                molecule.num_atoms(),
                molecule.bonds(),
                &molecule.topology_block().adjacency,
                AtomId::new(0),
            )
            .unwrap()
            .map(|bond| bond.id().index())
            .collect::<Vec<_>>(),
            vec![0]
        );
    }

    #[test]
    fn valence_parts_report_out_of_range_atom_access() {
        let molecule = Molecule::from_smiles_with_sanitize("C", false).unwrap();
        assert_eq!(
            super::atom_from_parts(molecule.atoms(), AtomId::new(1))
                .unwrap_err()
                .to_string(),
            "unsupported valence branch: atom index out of range"
        );
        match super::incident_bonds_from_parts(
            molecule.num_atoms(),
            molecule.bonds(),
            &molecule.topology_block().adjacency,
            AtomId::new(1),
        ) {
            Err(error) => assert_eq!(
                error.to_string(),
                "unsupported valence branch: atom index out of range"
            ),
            Ok(_) => panic!("expected out-of-range incident_bonds() to fail"),
        }
    }

    #[test]
    fn valence_parts_do_not_depend_on_cached_adjacency_from_registered_operation() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false)
            .unwrap()
            .with_assigned_valence()
            .unwrap();
        assert_eq!(
            super::incident_bonds_from_parts(
                molecule.num_atoms(),
                molecule.bonds(),
                &molecule.topology_block().adjacency,
                AtomId::new(1),
            )
            .unwrap()
            .map(|bond| bond.id().index())
            .collect::<Vec<_>>(),
            vec![0, 1]
        );
    }

    #[test]
    fn periodic_table_row_and_required_valence_list_match_rdkit_entries() {
        assert_eq!(super::periodic_table_row(0), Some(0));
        assert_eq!(super::periodic_table_row(1), Some(1));
        assert_eq!(super::periodic_table_row(92), Some(7));
        assert_eq!(super::periodic_table_row(118), Some(7));
        assert_eq!(super::periodic_table_row(119), None);

        assert_eq!(super::required_valence_list(0).unwrap(), &[-1]);
        assert_eq!(super::required_valence_list(16).unwrap(), &[2, 4, 6]);
        assert_eq!(super::required_valence_list(118).unwrap(), &[-1]);
        assert_eq!(
            super::required_valence_list(119).unwrap_err().to_string(),
            "unsupported valence branch: PeriodicTable atomic number out of range"
        );
    }

    #[test]
    fn get_effective_atomic_num_checks_then_clamps_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::O).with_formal_charge(9));
        builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(118).unwrap()).with_formal_charge(-10),
        );
        let molecule = builder.build().unwrap();
        let low_atom = &molecule.atoms()[0];
        let high_atom = &molecule.atoms()[1];

        let checked = super::get_effective_atomic_num(low_atom, true).unwrap_err();
        assert_eq!(checked.to_string(), "Effective atomic number out of range");
        assert_eq!(super::get_effective_atomic_num(low_atom, false).unwrap(), 0);
        assert!(super::get_effective_atomic_num(high_atom, true).is_err());
        assert_eq!(
            super::get_effective_atomic_num(high_atom, false).unwrap(),
            118
        );
        assert_eq!(
            super::rdkit_valence_list_for_atomic_number(119).unwrap(),
            None
        );
    }

    #[test]
    fn can_be_hypervalent_matches_rdkit_supported_element_pairs() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::from_atomic_number(15).unwrap()));
        builder.add_atom(AtomSpec::new(Element::from_atomic_number(16).unwrap()));
        builder.add_atom(AtomSpec::new(Element::from_atomic_number(33).unwrap()));
        builder.add_atom(AtomSpec::new(Element::from_atomic_number(34).unwrap()));
        builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let molecule = builder.build().unwrap();

        assert!(super::can_be_hypervalent(&molecule.atoms()[0], 17));
        assert!(super::can_be_hypervalent(&molecule.atoms()[1], 17));
        assert!(super::can_be_hypervalent(&molecule.atoms()[2], 35));
        assert!(super::can_be_hypervalent(&molecule.atoms()[3], 35));
        assert!(!super::can_be_hypervalent(&molecule.atoms()[0], 16));
        assert!(!super::can_be_hypervalent(&molecule.atoms()[4], 18));
    }

    #[test]
    fn has_bond_type_query_matches_rdkit_description_and_child_walk() {
        assert!(super::has_bond_type_query(&QueryNode::predicate(
            BondQueryPredicate::Order(BondOrder::Single)
        )));
        assert!(super::has_bond_type_query(&QueryNode::predicate(
            BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Aromatic])
        )));
        assert!(super::has_bond_type_query(&QueryNode::and(vec![
            QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        ])));
        assert!(super::has_bond_type_query(&QueryNode::not(
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single))
        )));
        assert!(!super::has_bond_type_query(&QueryNode::predicate(
            BondQueryPredicate::Any
        )));
        assert!(!super::has_bond_type_query(&QueryNode::predicate(
            BondQueryPredicate::IsInRing(true)
        )));
        assert!(!super::has_bond_type_query(&QueryNode::predicate(
            BondQueryPredicate::IsAromatic(true)
        )));
    }

    #[test]
    fn has_complex_bond_type_query_distinguishes_simple_and_complex_queries_like_rdkit() {
        assert!(!super::has_complex_bond_type_query(&QueryNode::predicate(
            BondQueryPredicate::Order(BondOrder::Single)
        )));
        assert!(!super::has_complex_bond_type_query(&QueryNode::predicate(
            BondQueryPredicate::OrderIn(vec![BondOrder::Single])
        )));
        assert!(super::has_complex_bond_type_query(&QueryNode::predicate(
            BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Double])
        )));
        assert!(super::has_complex_bond_type_query(&QueryNode::not(
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single))
        )));
        assert!(!super::has_complex_bond_type_query(&QueryNode::predicate(
            BondQueryPredicate::IsInRing(true)
        )));
    }

    #[test]
    fn assign_valence_updates_simple_aliphatic_atoms_like_rdkit_property_cache() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let assignment = assign_valence(&molecule, ValenceModel::RdkitLike).unwrap();
        assert_eq!(assignment.explicit_valence, vec![1, 2, 1]);
        assert_eq!(assignment.implicit_hydrogens, vec![3, 2, 1]);
    }

    #[test]
    fn assign_valence_counts_dative_valence_only_on_acceptor_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("N->O", false).unwrap();
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Dative);
        let assignment = assign_valence(&molecule, ValenceModel::RdkitLike).unwrap();
        assert_eq!(assignment.explicit_valence, vec![0, 1]);
        assert_eq!(assignment.implicit_hydrogens, vec![3, 1]);
    }

    #[test]
    fn assign_valence_rejects_dative_left_and_right_like_rdkit_default_bond_type_branch() {
        for order in [BondOrder::DativeLeft, BondOrder::DativeRight] {
            let mut builder = MoleculeBuilder::new();
            let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
            let oxygen = builder.add_atom(AtomSpec::new(Element::O));
            builder
                .add_bond(BondSpec::new(nitrogen, oxygen, order))
                .unwrap();
            let molecule = builder.build().unwrap();

            let error = assign_valence(&molecule, ValenceModel::RdkitLike).unwrap_err();
            assert_eq!(error.to_string(), "Bad bond type");
        }
    }

    #[test]
    fn assign_valence_handles_complex_query_bonds_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Single).with_query(
                QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                    BondOrder::Single,
                    BondOrder::Double,
                ])),
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let assignment = assign_valence(&molecule, ValenceModel::RdkitLike).unwrap();
        assert_eq!(assignment.explicit_valence, vec![0, 0]);
        assert_eq!(assignment.implicit_hydrogens, vec![0, 0]);
    }

    #[test]
    fn assign_valence_counts_simple_bond_order_queries_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Single).with_query(
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let assignment = assign_valence(&molecule, ValenceModel::RdkitLike).unwrap();
        assert_eq!(assignment.explicit_valence, vec![1, 1]);
        assert_eq!(assignment.implicit_hydrogens, vec![3, 1]);
    }

    #[test]
    fn assign_radicals_matches_rdkit_for_no_implicit_main_group_atoms() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        builder.add_atom(
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(3),
        );
        builder.add_atom(AtomSpec::new(Element::N).with_no_implicit(true));
        builder.add_atom(
            AtomSpec::new(Element::N)
                .with_no_implicit(true)
                .with_explicit_hydrogens(4)
                .with_formal_charge(1),
        );
        let molecule = builder.build().unwrap();
        let radicals = assign_radicals(&molecule).unwrap();
        assert_eq!(radicals, vec![4, 1, 3, 0]);
    }

    #[test]
    fn assign_radicals_matches_rdkit_for_no_preferred_valence_atoms() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(27).unwrap()).with_no_implicit(true),
        );
        builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(26).unwrap())
                .with_no_implicit(true)
                .with_formal_charge(9),
        );
        let bonded_fe = builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(26).unwrap()).with_no_implicit(true),
        );
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H).with_no_implicit(true));
        builder
            .add_bond(BondSpec::new(bonded_fe, hydrogen, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let radicals = assign_radicals(&molecule).unwrap();
        assert_eq!(radicals, vec![1, 0, 0, 0]);
    }

    #[test]
    fn assign_radicals_preserves_existing_radicals_for_skipped_atoms_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C).with_radical_electrons(2));
        builder.add_atom(
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_radical_electrons(2),
        );
        builder.add_atom(
            AtomSpec::new(Element::DUMMY)
                .with_no_implicit(true)
                .with_radical_electrons(3),
        );
        let molecule = builder.build().unwrap();

        let radicals = assign_radicals(&molecule).unwrap();
        assert_eq!(radicals, vec![2, 4, 3]);
    }

    #[test]
    fn assign_radicals_zeroes_bonded_no_preferred_valence_atoms_like_rdkit() {
        let iron = Element::from_atomic_number(26).unwrap();
        let mut builder = MoleculeBuilder::new();
        let fe = builder.add_atom(AtomSpec::new(iron).with_no_implicit(true));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H).with_no_implicit(true));
        builder
            .add_bond(BondSpec::new(fe, hydrogen, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let radicals = assign_radicals(&molecule).unwrap();
        assert_eq!(radicals, vec![0, 0]);
    }

    #[test]
    fn assign_radicals_uses_hypervalent_fallback_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let sulfur_element = Element::from_atomic_number(16).unwrap();
        let fluorine_element = Element::from_atomic_number(9).unwrap();
        let sulfur = builder.add_atom(
            AtomSpec::new(sulfur_element)
                .with_no_implicit(true)
                .with_formal_charge(-1),
        );
        for _ in 0..4 {
            let fluorine = builder.add_atom(AtomSpec::new(fluorine_element));
            builder
                .add_bond(BondSpec::new(sulfur, fluorine, BondOrder::Single))
                .unwrap();
        }
        let molecule = builder.build().unwrap();

        let radicals = assign_radicals(&molecule).unwrap();
        assert_eq!(radicals[0], 1);
    }

    #[test]
    fn assign_radicals_zeroes_unusual_charge_without_preferred_valence_like_rdkit() {
        let iron = Element::from_atomic_number(26).unwrap();
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(
            AtomSpec::new(iron)
                .with_no_implicit(true)
                .with_formal_charge(9),
        );
        let molecule = builder.build().unwrap();

        let radicals = assign_radicals(&molecule).unwrap();
        assert_eq!(radicals, vec![0]);
    }

    #[test]
    fn assign_valence_reports_invalid_strict_valence() {
        let molecule = Molecule::from_smiles_with_sanitize("C(=O)(=O)(=O)", false).unwrap();
        let error = assign_valence(&molecule, ValenceModel::RdkitLike).unwrap_err();
        assert_eq!(
            error.to_string(),
            "Explicit valence for atom # 0 C, 6, is greater than permitted"
        );
    }

    #[test]
    fn assign_valence_with_options_non_strict_keeps_property_cache_path_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("C(=O)(=O)(=O)", false).unwrap();

        let strict_error =
            assign_valence_with_options(&molecule, ValenceModel::RdkitLike, true).unwrap_err();
        assert_eq!(
            strict_error.to_string(),
            "Explicit valence for atom # 0 C, 6, is greater than permitted"
        );

        let non_strict =
            assign_valence_with_options(&molecule, ValenceModel::RdkitLike, false).unwrap();
        assert_eq!(non_strict.explicit_valence[0], 6);
    }

    #[test]
    fn calculate_explicit_valence_handles_hypervalent_anion_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("[S-](F)(F)(F)(F)F", false).unwrap();
        let explicit = super::calculate_explicit_valence(
            molecule.atoms(),
            molecule.bonds(),
            &molecule.topology_block().adjacency,
            AtomId::new(0),
            false,
            false,
        )
        .unwrap();
        assert_eq!(explicit, 5);
        assert_eq!(
            super::atom_calc_explicit_valence(
                molecule.atoms(),
                molecule.bonds(),
                &molecule.topology_block().adjacency,
                AtomId::new(0),
                true,
            )
            .unwrap(),
            5
        );
    }

    #[test]
    fn calculate_explicit_valence_returns_check_it_sentinel_for_overvalent_atom() {
        let molecule = Molecule::from_smiles_with_sanitize("C(=O)(=O)(=O)", false).unwrap();
        let explicit = super::calculate_explicit_valence(
            molecule.atoms(),
            molecule.bonds(),
            &molecule.topology_block().adjacency,
            AtomId::new(0),
            false,
            true,
        )
        .unwrap();
        assert_eq!(explicit, -1);
    }

    #[test]
    fn calculate_implicit_valence_handles_unreasonable_hydrogen_charge_by_strictness() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::H).with_formal_charge(2));
        let molecule = builder.build().unwrap();
        let strict_error = super::calculate_implicit_valence(
            molecule.atoms(),
            molecule.bonds(),
            &molecule.topology_block().adjacency,
            AtomId::new(0),
            0,
            true,
            false,
        )
        .unwrap_err();
        assert_eq!(
            strict_error.to_string(),
            "Unreasonable formal charge on atom # 0."
        );

        let non_strict = super::calculate_implicit_valence(
            molecule.atoms(),
            molecule.bonds(),
            &molecule.topology_block().adjacency,
            AtomId::new(0),
            0,
            false,
            false,
        )
        .unwrap();
        assert_eq!(non_strict, 0);

        let check_it = super::calculate_implicit_valence(
            molecule.atoms(),
            molecule.bonds(),
            &molecule.topology_block().adjacency,
            AtomId::new(0),
            0,
            false,
            true,
        )
        .unwrap();
        assert_eq!(check_it, -1);
    }

    #[test]
    fn calculate_implicit_valence_handles_hypervalent_anion_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("[S-](F)(F)(F)(F)F", false).unwrap();
        let explicit = super::calculate_explicit_valence(
            molecule.atoms(),
            molecule.bonds(),
            &molecule.topology_block().adjacency,
            AtomId::new(0),
            false,
            false,
        )
        .unwrap();

        let implicit = super::calculate_implicit_valence(
            molecule.atoms(),
            molecule.bonds(),
            &molecule.topology_block().adjacency,
            AtomId::new(0),
            explicit,
            false,
            false,
        )
        .unwrap();
        assert_eq!(implicit, 0);
        assert_eq!(
            super::atom_calc_implicit_valence(
                molecule.atoms(),
                molecule.bonds(),
                &molecule.topology_block().adjacency,
                AtomId::new(0),
                explicit,
                true,
            )
            .unwrap(),
            0
        );
    }

    #[test]
    fn assign_valence_rounds_half_aromatic_accumulation_up_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("C1ccccC1", false).unwrap();

        let assignment = assign_valence(&molecule, ValenceModel::RdkitLike).unwrap();
        assert_eq!(assignment.explicit_valence, vec![2, 3, 3, 3, 3, 2]);
        assert_eq!(assignment.implicit_hydrogens, vec![2, 1, 1, 1, 1, 2]);
    }

    #[test]
    fn assign_valence_accepts_hypervalent_negatively_charged_sulfur_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("[S-](F)(F)(F)(F)F", false).unwrap();

        let assignment = assign_valence(&molecule, ValenceModel::RdkitLike).unwrap();
        assert_eq!(assignment.explicit_valence, vec![5, 1, 1, 1, 1, 1]);
        assert_eq!(assignment.implicit_hydrogens, vec![0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn assign_valence_handles_hydrogen_charge_special_cases_like_rdkit() {
        let cases = [
            ("[H]", vec![0], vec![0]),
            ("[H+]", vec![0], vec![0]),
            ("[H-]", vec![0], vec![0]),
            ("[HH]", vec![1], vec![0]),
        ];

        for (smiles, explicit, implicit) in cases {
            let molecule = Molecule::from_smiles_with_sanitize(smiles, false).unwrap();
            let assignment = assign_valence(&molecule, ValenceModel::RdkitLike).unwrap();
            assert_eq!(assignment.explicit_valence, explicit, "{smiles}");
            assert_eq!(assignment.implicit_hydrogens, implicit, "{smiles}");
        }
    }

    #[test]
    fn atom_has_valence_violation_matches_rdkit_valid_and_invalid_examples() {
        for smiles in [
            "C",
            "C(C)(C)(C)C",
            "S(C)(C)(C)(C)(C)C",
            "O(C)C",
            "[H]",
            "[H+]",
            "[H-]",
            "[HH]",
            "[He]",
            "[CH3+]",
            "[CH3-]",
            "[NH4+]",
            "[Na]",
            "[Na][H]",
            "[Na]([H])[H]",
            "[Og][Og]([Og])([Og])([Og])([Og])([Og])[Og]",
            "[Lv-2]",
            "[Lv+4]",
            "[Lv+8]",
            "*",
        ] {
            let molecule = Molecule::from_smiles_with_sanitize(smiles, false).unwrap();
            for atom in molecule.atoms() {
                assert!(
                    !atom_has_valence_violation(&molecule, atom.id()).unwrap(),
                    "{smiles} atom {} should not have a valence violation",
                    atom.id()
                );
            }
        }

        for smiles in [
            "C(C)(C)(C)(C)C",
            "S(C)(C)(C)(C)(C)(C)C",
            "[C-](C)(C)(C)C",
            "O(C)=C",
            "[H+2]",
            "[O-3]",
            "[F-2]",
            "[Lv-4]",
        ] {
            let molecule = Molecule::from_smiles_with_sanitize(smiles, false).unwrap();
            assert!(
                atom_has_valence_violation(&molecule, AtomId::new(0)).unwrap(),
                "{smiles} atom 0 should have a valence violation"
            );
        }
    }

    #[test]
    fn atom_has_valence_violation_ignores_dummy_atoms_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::DUMMY));
        let molecule = builder.build().unwrap();

        assert!(!atom_has_valence_violation(&molecule, AtomId::new(0)).unwrap());
    }

    #[test]
    fn atom_has_valence_violation_flags_absurd_non_hydrogen_charge_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::O).with_formal_charge(9));
        let molecule = builder.build().unwrap();

        assert!(atom_has_valence_violation(&molecule, AtomId::new(0)).unwrap());
    }

    #[test]
    fn atom_has_valence_violation_flags_periodic_table_row_change_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let fluorine_element = Element::from_atomic_number(9).unwrap();
        builder.add_atom(AtomSpec::new(fluorine_element).with_formal_charge(-2));
        let molecule = builder.build().unwrap();

        assert!(atom_has_valence_violation(&molecule, AtomId::new(0)).unwrap());
    }

    #[test]
    fn atom_has_valence_violation_reports_non_strict_valence_sentinel_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("C(C)(C)(C)(C)C", false).unwrap();

        assert!(atom_has_valence_violation(&molecule, AtomId::new(0)).unwrap());
    }

    #[test]
    fn atom_has_valence_violation_flags_effective_atomic_number_overflow_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::H).with_formal_charge(-120));
        let molecule = builder.build().unwrap();
        assert!(atom_has_valence_violation(&molecule, AtomId::new(0)).unwrap());
        assert!(
            super::atom_has_valence_violation_impl(
                molecule.atoms(),
                molecule.bonds(),
                &molecule.topology_block().adjacency,
                AtomId::new(0),
            )
            .unwrap()
        );
    }

    #[test]
    fn atom_has_valence_violation_ignores_query_atoms_and_query_bonds_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let query = builder.add_atom(AtomSpec::new(Element::C).with_query(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumberIn(vec![6]),
        )));
        for _ in 0..5 {
            let carbon = builder.add_atom(AtomSpec::new(Element::C));
            builder
                .add_bond(BondSpec::new(query, carbon, BondOrder::Single))
                .unwrap();
        }
        let query_atom = builder.build().unwrap();
        assert!(!atom_has_valence_violation(&query_atom, AtomId::new(0)).unwrap());

        let mut builder = MoleculeBuilder::new();
        let left = builder.add_atom(AtomSpec::new(Element::C));
        let middle_left = builder.add_atom(AtomSpec::new(Element::C));
        let middle_right = builder.add_atom(AtomSpec::new(Element::C));
        let right = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(left, middle_left, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(
                BondSpec::new(middle_left, middle_right, BondOrder::Single)
                    .with_query(QueryNode::predicate(BondQueryPredicate::Any)),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(middle_right, right, BondOrder::Single))
            .unwrap();
        let query_bond = builder.build().unwrap();
        assert!(!atom_has_valence_violation(&query_bond, AtomId::new(1)).unwrap());
        assert!(!atom_has_valence_violation(&query_bond, AtomId::new(2)).unwrap());
    }

    #[test]
    fn assign_radicals_uses_light_element_base_count_two_branch_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(2).unwrap()).with_no_implicit(true),
        );
        builder.add_atom(AtomSpec::new(Element::H).with_no_implicit(true));
        let molecule = builder.build().unwrap();

        let radicals = assign_radicals(&molecule).unwrap();
        assert_eq!(radicals, vec![0, 1]);
    }

    #[test]
    fn has_complex_bond_type_query_marks_compound_bond_order_trees_like_rdkit() {
        assert!(super::has_complex_bond_type_query(&QueryNode::and(vec![
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        ])));
        assert!(super::has_complex_bond_type_query(&QueryNode::or(vec![
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
            ]),
        ])));
        assert!(!super::has_complex_bond_type_query(&QueryNode::and(vec![
            QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        ])));
    }

    #[test]
    fn with_assigned_valence_updates_derived_cache_through_registered_operation() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let original = molecule.clone();
        let result = molecule.with_assigned_valence().unwrap();

        assert_eq!(molecule, original);
        let adjacency = &result.topology_block().adjacency;
        assert_eq!(adjacency.neighbors_of(0).len(), 1);
        assert_eq!(adjacency.neighbors_of(1).len(), 2);
        assert_eq!(adjacency.neighbors_of(2).len(), 1);
        assert_eq!(
            result.derived_cache().valence,
            Some(ValenceAssignment {
                explicit_valence: vec![1, 2, 1],
                implicit_hydrogens: vec![3, 2, 1],
            })
        );
    }

    #[test]
    fn with_assigned_valence_fails_without_committing_invalid_valence_cache() {
        let molecule = Molecule::from_smiles_with_sanitize("C(=O)(=O)(=O)", false).unwrap();
        let original = molecule.clone();
        let error = molecule.with_assigned_valence().unwrap_err();

        assert_eq!(molecule, original);
        assert!(
            error
                .to_string()
                .contains("Explicit valence for atom # 0 C, 6, is greater than permitted")
        );
    }

    #[test]
    fn with_assigned_radicals_updates_topology_state_and_valence_cache_atomically() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(3),
        );
        let molecule = builder.build().unwrap();
        let original = molecule.clone();

        let result = molecule.with_assigned_radicals().unwrap();

        assert_eq!(molecule, original);
        assert_eq!(result.atoms()[0].radical_electrons(), 1);
        let adjacency = &result.topology_block().adjacency;
        assert_eq!(adjacency.neighbors_of(0).len(), 0);
        assert_eq!(
            result.derived_cache().valence,
            Some(ValenceAssignment {
                explicit_valence: vec![3],
                implicit_hydrogens: vec![0],
            })
        );
    }
}
