//! RDKit 2026.03.6 scaffold graph transformations over detached values.
//! MolHash and ChemTransforms have deliberately distinct scaffold semantics.
use crate::{RingFindType, RingInfo, ValenceAssignment};
use cosmolkit_model::{
    AtomId, BondOrder, ChiralTag, CoordinateBlock, Element, MoleculeProperties, TopologyBlock,
    TopologyEditError, TopologyMapping,
};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ScaffoldError {
    #[error(transparent)]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error(transparent)]
    Edit(#[from] TopologyEditError),
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
    #[error(transparent)]
    Rings(#[from] crate::RingFindingError),
    #[error(transparent)]
    Radicals(#[from] crate::RadicalError),
    #[error(transparent)]
    Stereo(#[from] crate::LegacyStereoError),
    #[error(transparent)]
    Property(#[from] cosmolkit_model::MoleculePropertyError),
    #[error(transparent)]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error(transparent)]
    BondProperty(#[from] cosmolkit_model::BondValueError),
    #[error(transparent)]
    UIntProperty(#[from] crate::PropertyUIntReadError),
    #[error(transparent)]
    TextProperty(#[from] crate::PropertyStringError),
    #[error(transparent)]
    Matrix(#[from] crate::MatrixError),
    #[error("atom {atom} has no calculated implicit hydrogen count")]
    MissingImplicitHydrogens { atom: AtomId },
    #[error("source valence {field} has {actual} rows, expected {expected}")]
    ValenceRowCount {
        field: &'static str,
        actual: usize,
        expected: usize,
    },
    #[error("MurckoDecompose requires initialized ring information")]
    MissingRingInfo,
    #[error("distance matrix dimension {atoms} overflows storage")]
    MatrixSize { atoms: usize },
}

/// Detached candidate; the facade alone validates and publishes its mapping.
pub struct ScaffoldResult {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
    pub mapping: TopologyMapping,
    pub valence: Option<ValenceAssignment>,
    pub rings: Option<RingInfo>,
}

fn result(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
) -> ScaffoldResult {
    let mapping = TopologyMapping::identity(topology.atoms.len(), topology.bonds.len());
    ScaffoldResult {
        topology,
        coordinates,
        properties,
        mapping,
        valence: None,
        rings: None,
    }
}

fn total_hydrogens(topology: &TopologyBlock, id: AtomId) -> Result<u32, ScaffoldError> {
    // Delegates the pinned Atom::getNumImplicitHs precondition and noImplicit
    // shortcut to its existing sole owner; explicit neighboring Hs are excluded.
    let atom = &topology.atoms[id.index()];
    let implicit = crate::valence::source_atom_implicit_hydrogens(atom)
        .map_err(|_| ScaffoldError::MissingImplicitHydrogens { atom: id })?;
    Ok(u32::from(atom.explicit_hydrogens()) + implicit)
}

fn sanitize_hydrogens(
    topology: &mut TopologyBlock,
    source_valence: Option<&ValenceAssignment>,
) -> Result<(), ScaffoldError> {
    // RDKit✔️✔️: void NMRDKitSanitizeHydrogens(RDKit::RWMol *mol) {
    // RDKit✔️✔️:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️✔️:   // Move all of the implicit Hs into one box
    // RDKit✔️✔️:   for (auto aptr : mol->atoms()) {
    // RDKit✔️✔️:     unsigned int hcount = aptr->getTotalNumHs();
    // RDKit✔️✔️:     aptr->setNoImplicit(true);
    // RDKit✔️✔️:     aptr->setNumExplicitHs(hcount);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     bool strict = false;
    // RDKit✔️✔️:     aptr->updatePropertyCache(
    // RDKit✔️✔️:         strict);  // or else the valence is reported incorrectly
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local cost review: one atom pass, source cache owner scans each adjacency;
    // no graph rebuild or second whole-topology clone.
    if let Some(valence) = source_valence {
        for (field, actual) in [
            ("explicit_valence", valence.explicit_valence.len()),
            ("implicit_hydrogens", valence.implicit_hydrogens.len()),
        ] {
            if actual != topology.atoms.len() {
                return Err(ScaffoldError::ValenceRowCount {
                    field,
                    actual,
                    expected: topology.atoms.len(),
                });
            }
        }
    }
    for i in 0..topology.atoms.len() {
        let id = AtomId::new(i);
        // CK stores prepared source cache scalars in the detached assignment;
        // this is the same getter value, not a new valence calculation.
        let count = if let Some(valence) = source_valence {
            let atom = &topology.atoms[i];
            let implicit = if atom.no_implicit() {
                0
            } else {
                valence.implicit_hydrogens[i]
            };
            let implicit = u32::try_from(implicit)
                .map_err(|_| ScaffoldError::MissingImplicitHydrogens { atom: id })?;
            u32::from(atom.explicit_hydrogens()) + implicit
        } else {
            total_hydrogens(topology, id)?
        };
        topology.atoms[i].set_no_implicit(true);
        // RDKit stores this unsigned count in its uint8_t atom field.
        topology.atoms[i].set_explicit_hydrogens(count as u8);
        crate::valence::update_source_atom_cache(topology, id, false)?;
    }
    Ok(())
}

fn bond_order(order: BondOrder) -> u32 {
    // RDKit✔️✔️: unsigned int NMRDKitBondGetOrder(const RDKit::Bond *bnd) {
    // RDKit✔️✔️:   PRECONDITION(bnd, "bad bond");
    // RDKit✔️✔️:   switch (bnd->getBondType()) {
    // RDKit✔️✔️:     case RDKit::Bond::AROMATIC:
    // RDKit✔️✔️:     case RDKit::Bond::SINGLE:
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     case RDKit::Bond::DOUBLE:
    // RDKit✔️✔️:       return 2;
    // RDKit✔️✔️:     case RDKit::Bond::TRIPLE:
    // RDKit✔️✔️:       return 3;
    // RDKit✔️✔️:     case RDKit::Bond::QUADRUPLE:
    // RDKit✔️✔️:       return 4;
    // RDKit✔️✔️:     case RDKit::Bond::QUINTUPLE:
    // RDKit✔️✔️:       return 5;
    // RDKit✔️✔️:     case RDKit::Bond::HEXTUPLE:
    // RDKit✔️✔️:       return 6;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    match order {
        BondOrder::Aromatic | BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Quintuple => 5,
        BondOrder::Hextuple => 6,
        _ => 0,
    }
}

fn remove_rows(value: &mut ScaffoldResult, removed: &[AtomId]) -> Result<(), ScaffoldError> {
    // The sole source batch kernel owns neighbor/stereo/SGroups/property/coordinate
    // transport. Moving the topology into it avoids the borrowed edit's clone.
    let topology = std::mem::take(&mut value.topology);
    let mut edit = topology.into_batch_edit()?;
    for id in removed {
        edit.remove_atom(*id)?;
    }
    let (topology, mapping) = edit.finish_source(
        &mut value.coordinates,
        &mut value.properties,
        &mut |v| crate::property_value_to_uint(v).map_err(ScaffoldError::from),
        &mut |v| crate::property_value_to_string(v).map_err(ScaffoldError::from),
    )?;
    // Compose original -> current -> next without rebuilding retained atom/bond
    // values. Both directions remain reciprocal and are validated at commit.
    for row in &mut value.mapping.atoms.old_to_new {
        *row = row.and_then(|id| mapping.atoms.old_to_new[id.index()]);
    }
    value.mapping.atoms.new_to_old = mapping
        .atoms
        .new_to_old
        .iter()
        .map(|id| id.and_then(|id| value.mapping.atoms.new_to_old[id.index()]))
        .collect();
    for row in &mut value.mapping.bonds.old_to_new {
        *row = row.and_then(|id| mapping.bonds.old_to_new[id.index()]);
    }
    value.mapping.bonds.new_to_old = mapping
        .bonds
        .new_to_old
        .iter()
        .map(|id| id.and_then(|id| value.mapping.bonds.new_to_old[id.index()]))
        .collect();
    value.topology = topology;
    Ok(())
}

fn finish_hash(
    value: &mut ScaffoldResult,
    source_rings: Option<&RingInfo>,
    prepared_rings: Option<RingInfo>,
) -> Result<(), ScaffoldError> {
    // RDKit✔️❌: MolOps::assignRadicals(*mol);
    // RDKit✔️❌: bool cleanIt = true;
    // RDKit✔️❌: bool force = true;
    // RDKit✔️❌: MolOps::assignStereochemistry(*mol, cleanIt, force);
    // This is the graph-state boundary preceding MolHash's string rendering:
    // CK returns the transformed graph, rather than parsing its SMILES back.
    // The existing detached radical owner allocates a result row per atom;
    // that known allocation difference is retained in the cost marker.
    let assignment = crate::assign_radicals(&value.topology)?;
    for (atom, count) in value
        .topology
        .atoms
        .iter_mut()
        .zip(assignment.radical_electrons)
    {
        atom.set_radical_electrons(count);
    }
    // assignStereochemistry updates ALL atom caches only if the source
    // molecule-level predicate finds at least one unprepared atom.
    if value
        .topology
        .atoms
        .iter()
        .any(crate::valence::source_atom_needs_cache_update)
    {
        for i in 0..value.topology.atoms.len() {
            crate::valence::update_source_atom_cache(&mut value.topology, AtomId::new(i), false)?;
        }
    }
    let valence = ValenceAssignment {
        explicit_valence: value
            .topology
            .atoms
            .iter()
            .map(|a| i32::from(a.source_valence_facts().explicit_valence))
            .collect(),
        implicit_hydrogens: value
            .topology
            .atoms
            .iter()
            .map(|a| {
                if a.no_implicit() {
                    0
                } else {
                    i32::from(a.source_valence_facts().implicit_valence)
                }
            })
            .collect(),
    };
    let found = if prepared_rings.is_some() {
        prepared_rings
    } else if source_rings
        .is_some_and(|r| r.is_initialized() && r.find_type() != RingFindType::OtherOrUnknown)
    {
        None
    } else {
        Some(crate::fast_find_rings(&value.topology)?)
    };
    let rings = found
        .as_ref()
        .or(source_rings)
        .expect("initialized source or newly found rings");
    let mut ring_update = None;
    crate::assign_legacy_stereochemistry_source(
        &mut value.topology,
        &valence,
        &rings,
        None,
        true,
        false,
        &mut ring_update,
    )?;
    value
        .properties
        .set_computed_prop("_StereochemDone", 1_i32)?;
    value.valence = Some(valence);
    // None means the caller's existing cache is retained; do not clone a
    // borrowed SymmSSSR just to reinstall it or replace it with Fast rings.
    value.rings = ring_update.or(found);
    Ok(())
}

/// MolHash::MurckoScaffold graph state (not ChemTransforms::MurckoDecompose).
pub fn murcko_scaffold(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    source_rings: Option<&RingInfo>,
    source_valence: Option<&ValenceAssignment>,
) -> Result<ScaffoldResult, ScaffoldError> {
    // Source: GraphMol/MolHash/hashfunctions.cpp, MurckoScaffoldHash.
    // RDKit✔️❌: std::string MurckoScaffoldHash(RWMol *mol, bool useCXSmiles,
    // RDKit✔️❌:                                unsigned cxFlagsToSkip = 0) {
    // RDKit✔️❌:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❌:   std::vector<Atom *> for_deletion;
    // RDKit✔️❌:   do {
    // RDKit✔️❌:     for_deletion.clear();
    // RDKit✔️❌:     for (auto aptr : mol->atoms()) {
    // RDKit✔️❌:       unsigned int deg = aptr->getDegree();
    // RDKit✔️❌:       if (deg < 2) {
    // RDKit✔️❌:         if (deg == 1) {  // i.e. not 0 and the last atom in the molecule
    // RDKit✔️❌:           for (const auto &nbri : boost::make_iterator_range(
    // RDKit✔️❌:                    aptr->getOwningMol().getAtomBonds(aptr))) {
    // RDKit✔️❌:             auto bptr = (aptr->getOwningMol())[nbri];
    // RDKit✔️❌:             Atom *nbr = bptr->getOtherAtom(aptr);
    // RDKit✔️❌:             unsigned int hcount = nbr->getTotalNumHs(false);
    // RDKit✔️❌:             nbr->setNumExplicitHs(hcount + NMRDKitBondGetOrder(bptr));
    // RDKit✔️❌:             nbr->setNoImplicit(true);
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:         for_deletion.push_back(aptr);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     mol->beginBatchEdit();
    // RDKit✔️❌:     for (auto &i : for_deletion) {
    // RDKit✔️❌:       mol->removeAtom(i);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     mol->commitBatchEdit();
    // RDKit✔️❌:   } while (!for_deletion.empty());
    // RDKit✔️❌:   MolOps::assignRadicals(*mol);
    // RDKit✔️❌:
    // RDKit✔️❌:   // we may have just destroyed some stereocenters/bonds
    // RDKit✔️❌:   // clean that up:
    // RDKit✔️❌:   bool cleanIt = true;
    // RDKit✔️❌:   bool force = true;
    // RDKit✔️❌:   MolOps::assignStereochemistry(*mol, cleanIt, force);
    // RDKit✔️❌:
    // RDKit❌❌:   std::string result;
    // RDKit❌❌:   result = convertToSmilesWithCXFlags(*mol, useCXSmiles, cxFlagsToSkip);
    // RDKit❌❌:   if (useCXSmiles) {
    // RDKit❌❌:     addCXExtensions(mol, result, cxFlagsToSkip | SmilesWrite::CX_RADICALS);
    // RDKit❌❌:   }
    // RDKit❌❌:   return result;
    // RDKit✔️❌: }
    // Source MolHash first normalizes implicit Hs; deletion uses the shared
    // batch-edit kernel. Mapping composition adds O(V+E) rows per round but
    // never clones the surviving graph. Radical-owner allocation is qualified.
    topology.validate()?;
    let mut value = result(topology, coordinates, properties);
    sanitize_hydrogens(&mut value.topology, source_valence)?;
    let mut deleted = false;
    loop {
        let mut removed = Vec::new();
        for i in 0..value.topology.atoms.len() {
            let neighbors = value.topology.adjacency.neighbors_of(i);
            if neighbors.len() < 2 {
                if let Some(neighbor) = neighbors.first() {
                    let id = AtomId::new(neighbor.atom_index);
                    let count = total_hydrogens(&value.topology, id)?
                        + bond_order(value.topology.bonds[neighbor.bond.index()].order());
                    value.topology.atoms[id.index()].set_explicit_hydrogens(count as u8);
                    value.topology.atoms[id.index()].set_no_implicit(true);
                }
                removed.push(AtomId::new(i));
            }
        }
        remove_rows(&mut value, &removed)?;
        if removed.is_empty() {
            break;
        }
        deleted = true;
    }
    finish_hash(&mut value, if deleted { None } else { source_rings }, None)?;
    Ok(value)
}

fn traverse_for_ring(
    topology: &TopologyBlock,
    atom: usize,
    rings: &RingInfo,
    visited: &mut [bool],
) -> bool {
    // RDKit✔️✔️: bool TraverseForRing(Atom *atom, unsigned char *visit) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom pointer");
    // RDKit✔️✔️:   PRECONDITION(visit, "bad pointer");
    // RDKit✔️✔️:   visit[atom->getIdx()] = 1;
    // RDKit✔️✔️:   for (auto nbri : boost::make_iterator_range(
    // RDKit✔️✔️:            atom->getOwningMol().getAtomNeighbors(atom))) {
    // RDKit✔️✔️:     auto nptr = atom->getOwningMol()[nbri];
    // RDKit✔️✔️:     if (visit[nptr->getIdx()] == 0) {
    // RDKit✔️✔️:       if (RDKit::queryIsAtomInRing(nptr)) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (TraverseForRing(nptr, visit)) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    visited[atom] = true;
    for neighbor in topology.adjacency.neighbors_of(atom) {
        let next = neighbor.atom_index;
        if !visited[next]
            && (rings.num_atom_rings(AtomId::new(next)) != 0
                || traverse_for_ring(topology, next, rings, visited))
        {
            return true;
        }
    }
    false
}

fn reaches_ring(topology: &TopologyBlock, root: usize, neighbor: usize, rings: &RingInfo) -> bool {
    // RDKit✔️✔️: bool DepthFirstSearchForRing(Atom *root, Atom *nbor, unsigned int maxatomidx) {
    // RDKit✔️✔️:   PRECONDITION(root, "bad atom pointer");
    // RDKit✔️✔️:   PRECONDITION(nbor, "bad atom pointer");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int natoms = maxatomidx;
    // RDKit✔️✔️:   std::vector<unsigned char> visit(natoms, 0);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   visit[root->getIdx()] = true;
    // RDKit✔️✔️:   return TraverseForRing(nbor, visit.data());
    // RDKit✔️✔️: }
    let mut visited = vec![false; topology.atoms.len()];
    visited[root] = true;
    traverse_for_ring(topology, neighbor, rings, &mut visited)
}

fn in_scaffold(topology: &TopologyBlock, atom: usize, rings: &RingInfo) -> bool {
    // RDKit✔️✔️: bool IsInScaffold(Atom *atom, unsigned int maxatomidx) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom pointer");
    // RDKit✔️✔️:   if (RDKit::queryIsAtomInRing(atom)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int count = 0;
    // RDKit✔️✔️:   for (auto nbri : boost::make_iterator_range(
    // RDKit✔️✔️:            atom->getOwningMol().getAtomNeighbors(atom))) {
    // RDKit✔️✔️:     auto nptr = atom->getOwningMol()[nbri];
    // RDKit✔️✔️:     if (DepthFirstSearchForRing(atom, nptr, maxatomidx)) {
    // RDKit✔️✔️:       ++count;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return count > 1;
    // RDKit✔️✔️: }
    if rings.num_atom_rings(AtomId::new(atom)) != 0 {
        return true;
    }
    topology
        .adjacency
        .neighbors_of(atom)
        .iter()
        .filter(|n| reaches_ring(topology, atom, n.atom_index, rings))
        .count()
        > 1
}

fn has_scaffold_neighbor(topology: &TopologyBlock, atom: usize, keep: &[bool]) -> bool {
    // RDKit✔️✔️: bool HasNbrInScaffold(Atom *aptr, unsigned char *is_in_scaffold) {
    // RDKit✔️✔️:   PRECONDITION(aptr, "bad atom pointer");
    // RDKit✔️✔️:   PRECONDITION(is_in_scaffold, "bad pointer");
    // RDKit✔️✔️:   for (auto nbri : boost::make_iterator_range(
    // RDKit✔️✔️:            aptr->getOwningMol().getAtomNeighbors(aptr))) {
    // RDKit✔️✔️:     auto nptr = aptr->getOwningMol()[nbri];
    // RDKit✔️✔️:     if (is_in_scaffold[nptr->getIdx()]) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    topology
        .adjacency
        .neighbors_of(atom)
        .iter()
        .any(|n| keep[n.atom_index])
}

/// MolHash::ExtendedMurcko graph state: retain the first substituent as dummy.
pub fn net_scaffold(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    source_rings: Option<&RingInfo>,
    source_valence: Option<&ValenceAssignment>,
) -> Result<ScaffoldResult, ScaffoldError> {
    // Source: GraphMol/MolHash/hashfunctions.cpp, ExtendedMurckoScaffold.
    // RDKit✔️❌: std::string ExtendedMurckoScaffold(RWMol *mol, bool useCXSmiles,
    // RDKit✔️❌:                                    unsigned cxFlagsToSkip = 0) {
    // RDKit✔️❌:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❌:   if (!mol->getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️❌:     MolOps::fastFindRings(*mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   unsigned int maxatomidx = mol->getNumAtoms();
    // RDKit✔️❌:   std::vector<unsigned char> is_in_scaffold(maxatomidx);
    // RDKit✔️❌:   for (auto aptr : mol->atoms()) {
    // RDKit✔️❌:     is_in_scaffold[aptr->getIdx()] = IsInScaffold(aptr, maxatomidx);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   std::vector<Atom *> for_deletion;
    // RDKit✔️❌:   for (auto aptr : mol->atoms()) {
    // RDKit✔️❌:     unsigned int aidx = aptr->getIdx();
    // RDKit✔️❌:     if (is_in_scaffold[aidx]) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (HasNbrInScaffold(aptr, is_in_scaffold.data())) {
    // RDKit✔️❌:       aptr->setAtomicNum(0);
    // RDKit✔️❌:       aptr->setFormalCharge(0);
    // RDKit✔️❌:       aptr->setNoImplicit(true);
    // RDKit✔️❌:       aptr->setNumExplicitHs(0);
    // RDKit✔️❌:       aptr->setIsotope(0);
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       for_deletion.push_back(aptr);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   mol->beginBatchEdit();
    // RDKit✔️❌:   for (auto &i : for_deletion) {
    // RDKit✔️❌:     mol->removeAtom(i);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   mol->commitBatchEdit();
    // RDKit✔️❌:   MolOps::assignRadicals(*mol);
    // RDKit✔️❌:
    // RDKit✔️❌:   // we may have just destroyed some stereocenters/bonds
    // RDKit✔️❌:   // clean that up:
    // RDKit✔️❌:   bool cleanIt = true;
    // RDKit✔️❌:   bool force = true;
    // RDKit✔️❌:   MolOps::assignStereochemistry(*mol, cleanIt, force);
    // RDKit✔️❌:
    // RDKit❌❌:   std::string result;
    // RDKit❌❌:   result = convertToSmilesWithCXFlags(*mol, useCXSmiles, cxFlagsToSkip);
    // RDKit❌❌:   if (useCXSmiles) {
    // RDKit❌❌:     addCXExtensions(mol, result, cxFlagsToSkip | SmilesWrite::CX_RADICALS);
    // RDKit❌❌:   }
    // RDKit❌❌:   return result;
    // RDKit✔️❌: }
    topology.validate()?;
    let mut value = result(topology, coordinates, properties);
    sanitize_hydrogens(&mut value.topology, source_valence)?;
    let found = if source_rings
        .is_some_and(|r| r.is_initialized() && r.find_type() != RingFindType::OtherOrUnknown)
    {
        None
    } else {
        Some(crate::fast_find_rings(&value.topology)?)
    };
    let rings = found
        .as_ref()
        .or(source_rings)
        .expect("initialized source or newly found rings");
    let keep: Vec<_> = (0..value.topology.atoms.len())
        .map(|i| in_scaffold(&value.topology, i, rings))
        .collect();
    let mut removed = Vec::new();
    for i in 0..keep.len() {
        if keep[i] {
            continue;
        }
        if has_scaffold_neighbor(&value.topology, i, &keep) {
            let atom = &mut value.topology.atoms[i];
            atom.set_element(Element::DUMMY);
            atom.set_formal_charge(0);
            atom.set_no_implicit(true);
            atom.set_explicit_hydrogens(0);
            atom.set_isotope(None);
        } else {
            removed.push(AtomId::new(i));
        }
    }
    remove_rows(&mut value, &removed)?;
    finish_hash(
        &mut value,
        if removed.is_empty() {
            source_rings
        } else {
            None
        },
        if removed.is_empty() { found } else { None },
    )?;
    Ok(value)
}

fn predecessor_matrix(topology: &TopologyBlock) -> Result<Vec<i32>, ScaffoldError> {
    // RDKit✔️❌: template <class T>
    // RDKit✔️❌: void FloydWarshall(int dim, T *adjMat, int *pathMat) {
    // RDKit✔️❌:   int k, i, j;
    // RDKit✔️❌:   T *currD, *lastD, *tTemp;
    // RDKit✔️❌:   int *currP, *lastP, *iTemp;
    // RDKit✔️❌:
    // RDKit✔️❌:   currD = new T[dim * dim];
    // RDKit✔️❌:   currP = new int[dim * dim];
    // RDKit✔️❌:   lastD = new T[dim * dim];
    // RDKit✔️❌:   lastP = new int[dim * dim];
    // RDKit✔️❌:
    // RDKit✔️❌:   memcpy(static_cast<void *>(lastD), static_cast<void *>(adjMat),
    // RDKit✔️❌:          dim * dim * sizeof(T));
    // RDKit✔️❌:
    // RDKit✔️❌:   // initialize the paths
    // RDKit✔️❌:   for (i = 0; i < dim; i++) {
    // RDKit✔️❌:     int itab = i * dim;
    // RDKit✔️❌:     for (j = 0; j < dim; j++) {
    // RDKit✔️❌:       if (i == j || adjMat[itab + j] == LOCAL_INF) {
    // RDKit✔️❌:         pathMat[itab + j] = -1;
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         pathMat[itab + j] = i;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   memcpy(static_cast<void *>(lastP), static_cast<void *>(pathMat),
    // RDKit✔️❌:          dim * dim * sizeof(int));
    // RDKit✔️❌:
    // RDKit✔️❌:   for (k = 0; k < dim; k++) {
    // RDKit✔️❌:     int ktab = k * dim;
    // RDKit✔️❌:     for (i = 0; i < dim; i++) {
    // RDKit✔️❌:       int itab = i * dim;
    // RDKit✔️❌:       for (j = 0; j < dim; j++) {
    // RDKit✔️❌:         T v1 = lastD[itab + j];
    // RDKit✔️❌:         T v2 = lastD[itab + k] + lastD[ktab + j];
    // RDKit✔️❌:         if (v1 <= v2) {
    // RDKit✔️❌:           currD[itab + j] = v1;
    // RDKit✔️❌:           currP[itab + j] = lastP[itab + j];
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           currD[itab + j] = v2;
    // RDKit✔️❌:           currP[itab + j] = lastP[ktab + j];
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     tTemp = currD;
    // RDKit✔️❌:     currD = lastD;
    // RDKit✔️❌:     lastD = tTemp;
    // RDKit✔️❌:
    // RDKit✔️❌:     iTemp = currP;
    // RDKit✔️❌:     currP = lastP;
    // RDKit✔️❌:     lastP = iTemp;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   memcpy(static_cast<void *>(adjMat), static_cast<void *>(lastD),
    // RDKit✔️❌:          dim * dim * sizeof(T));
    // RDKit✔️❌:   memcpy(static_cast<void *>(pathMat), static_cast<void *>(lastP),
    // RDKit✔️❌:          dim * dim * sizeof(int));
    // RDKit✔️❌:
    // RDKit✔️❌:   delete[] currD;
    // RDKit✔️❌:   delete[] currP;
    // RDKit✔️❌:   delete[] lastD;
    // RDKit✔️❌:   delete[] lastP;
    // RDKit✔️❌: }
    // getDistanceMat(mol, false, false, true): all bonds contribute exactly 1.
    // Preserve the source's <= tie decision and previous-k buffers. O(V^3)
    // work/O(V^2) storage; no BFS substitution that changes equal-length paths.
    // The existing matrix owner implements this verbatim source algorithm;
    // reuse it rather than maintaining another Floyd-Warshall implementation.
    // RDKit✔️❌: this shared owner has extra quadratic scratch copies and an
    // active-row vector compared with the all-atoms native overload above.
    let n = topology.atoms.len();
    let size = n
        .checked_mul(n)
        .ok_or(ScaffoldError::MatrixSize { atoms: n })?;
    let mut previous_d = vec![100_000_000.0; size];
    for i in 0..n {
        previous_d[i * n + i] = 0.0;
    }
    for bond in &topology.bonds {
        let i = bond.begin().index();
        let j = bond.end().index();
        previous_d[i * n + j] = 1.0;
        previous_d[j * n + i] = 1.0;
    }
    let active = (0..n).collect::<Vec<_>>();
    Ok(crate::matrices::floyd_warshall(
        n,
        &mut previous_d,
        &active,
    )?)
}

/// ChemTransforms::MurckoDecompose; preserve ring-exocyclic double bonds.
pub fn murcko_decompose(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    rings: Option<&RingInfo>,
) -> Result<ScaffoldResult, ScaffoldError> {
    // Source: GraphMol/ChemTransforms/ChemTransforms.cpp, MurckoDecompose.
    // RDKit✔️❌: ROMol *MurckoDecompose(const ROMol &mol) {
    // RDKit✔️❌:   auto *res = new RWMol(mol);
    // RDKit✔️❌:   unsigned int nAtoms = res->getNumAtoms();
    // RDKit✔️❌:   if (!nAtoms) {
    // RDKit✔️❌:     return res;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // start by getting the shortest paths matrix:
    // RDKit✔️❌:   MolOps::getDistanceMat(mol, false, false, true);
    // RDKit✔️❌:   boost::shared_array<int> pathMat;
    // RDKit✔️❌:   mol.getProp(common_properties::DistanceMatrix_Paths, pathMat);
    // RDKit✔️❌:
    // RDKit✔️❌:   boost::dynamic_bitset<> keepAtoms(nAtoms);
    // RDKit✔️❌:   const RingInfo *ringInfo = res->getRingInfo();
    // RDKit✔️❌:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️❌:     if (ringInfo->numAtomRings(i)) {
    // RDKit✔️❌:       keepAtoms[i] = 1;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   const VECT_INT_VECT &rings = ringInfo->atomRings();
    // RDKit✔️❌:   // std::cerr<<"  rings: "<<rings.size()<<std::endl;
    // RDKit✔️❌:   // now find the shortest paths between each ring system and mark the atoms
    // RDKit✔️❌:   // along each as being keepers:
    // RDKit✔️❌:   for (auto ringsItI = rings.begin(); ringsItI != rings.end(); ++ringsItI) {
    // RDKit✔️❌:     for (auto ringsItJ = ringsItI + 1; ringsItJ != rings.end(); ++ringsItJ) {
    // RDKit✔️❌:       int atomI = (*ringsItI)[0];
    // RDKit✔️❌:       int atomJ = (*ringsItJ)[0];
    // RDKit✔️❌:       // std::cerr<<atomI<<" -> "<<atomJ<<": ";
    // RDKit✔️❌:       while (atomI != atomJ) {
    // RDKit✔️❌:         keepAtoms[atomI] = 1;
    // RDKit✔️❌:         atomI = pathMat[atomJ * nAtoms + atomI];
    // RDKit✔️❌:         // test for the disconnected case:
    // RDKit✔️❌:         if (atomI < 0) {
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         // std::cerr<<atomI<<" ";
    // RDKit✔️❌:       }
    // RDKit✔️❌:       // std::cerr<<std::endl;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   boost::dynamic_bitset<> removedAtoms(nAtoms);
    // RDKit✔️❌:   res->beginBatchEdit();
    // RDKit✔️❌:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️❌:     if (!keepAtoms[i]) {
    // RDKit✔️❌:       Atom *atom = res->getAtomWithIdx(i);
    // RDKit✔️❌:       bool removeIt = true;
    // RDKit✔️❌:
    // RDKit✔️❌:       // check if the atom has a neighboring keeper:
    // RDKit✔️❌:       for (auto nbr : res->atomNeighbors(atom)) {
    // RDKit✔️❌:         if (keepAtoms[nbr->getIdx()]) {
    // RDKit✔️❌:           if (res->getBondBetweenAtoms(atom->getIdx(), nbr->getIdx())
    // RDKit✔️❌:                   ->getBondType() == Bond::DOUBLE) {
    // RDKit✔️❌:             removeIt = false;
    // RDKit✔️❌:             break;
    // RDKit✔️❌:           } else if (nbr->getIsAromatic() && nbr->getAtomicNum() != 6) {
    // RDKit✔️❌:             // fix aromatic heteroatoms:
    // RDKit✔️❌:             nbr->setNumExplicitHs(1);
    // RDKit✔️❌:           } else if (nbr->getIsAromatic() && nbr->getAtomicNum() == 6 &&
    // RDKit✔️❌:                      nbr->getFormalCharge() == 1) {
    // RDKit✔️❌:             // fix aromatic carbocations
    // RDKit✔️❌:             nbr->setNumExplicitHs(1);
    // RDKit✔️❌:           } else if (nbr->getNoImplicit() ||
    // RDKit✔️❌:                      nbr->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️❌:             nbr->setNoImplicit(false);
    // RDKit✔️❌:             nbr->setNumExplicitHs(0);
    // RDKit✔️❌:             nbr->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:
    // RDKit✔️❌:       if (removeIt) {
    // RDKit✔️❌:         res->removeAtom(atom);
    // RDKit✔️❌:         removedAtoms.set(atom->getIdx());
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   res->commitBatchEdit();
    // RDKit✔️❌:
    // RDKit✔️❌:   details::updateSubMolConfs(mol, *res, removedAtoms);
    // RDKit✔️❌:   res->clearComputedProps();
    // RDKit✔️❌:
    // RDKit✔️❌:   return (ROMol *)res;
    // RDKit✔️❌: }
    // Unlike the MolHash transforms this function does not normalize Hs,
    // assign radicals, sanitize or assign stereo after deletion.
    topology.validate()?;
    let mut value = result(topology, coordinates, properties);
    let n = value.topology.atoms.len();
    if n == 0 {
        return Ok(value);
    }
    let paths = predecessor_matrix(&value.topology)?;
    let rings = rings
        .filter(|r| r.is_initialized())
        .ok_or(ScaffoldError::MissingRingInfo)?;
    let mut keep: Vec<_> = (0..n)
        .map(|i| rings.num_atom_rings(AtomId::new(i)) != 0)
        .collect();
    for (i, first) in rings.atom_rings().iter().enumerate() {
        for second in &rings.atom_rings()[i + 1..] {
            let mut atom = first[0].index();
            let target = second[0].index();
            while atom != target {
                keep[atom] = true;
                let next = paths[target * n + atom];
                if next < 0 {
                    break;
                }
                atom = next as usize;
            }
        }
    }
    let mut removed = Vec::new();
    for i in 0..n {
        if keep[i] {
            continue;
        }
        let mut remove = true;
        for neighbor in value.topology.adjacency.neighbors_of(i) {
            let next = neighbor.atom_index;
            if keep[next] {
                if value.topology.bonds[neighbor.bond.index()].order() == BondOrder::Double {
                    remove = false;
                    break;
                }
                let atom = &mut value.topology.atoms[next];
                if atom.is_aromatic() && (atom.atomic_number() != 6 || atom.formal_charge() == 1) {
                    atom.set_explicit_hydrogens(1);
                } else if atom.no_implicit() || atom.chiral_tag() != ChiralTag::Unspecified {
                    atom.set_no_implicit(false);
                    atom.set_explicit_hydrogens(0);
                    atom.set_chiral_tag(ChiralTag::Unspecified);
                }
            }
        }
        if remove {
            removed.push(AtomId::new(i));
        }
    }
    remove_rows(&mut value, &removed)?;
    for atom in &mut value.topology.atoms {
        atom.clear_computed_props()?;
    }
    for bond in &mut value.topology.bonds {
        bond.clear_computed_props()?;
    }
    value.properties.clear_computed_props()?;
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, Bond, BondId, BondSpec};
    fn carbonyl_ring() -> TopologyBlock {
        let atoms = (0..7)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(if i == 6 { Element::O } else { Element::C }),
                )
            })
            .collect();
        let edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6)];
        let bonds = edges
            .into_iter()
            .enumerate()
            .map(|(i, (a, b))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(
                        AtomId::new(a),
                        AtomId::new(b),
                        if i == 6 {
                            BondOrder::Double
                        } else {
                            BondOrder::Single
                        },
                    ),
                )
            })
            .collect();
        let mut graph = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        for i in 0..7 {
            crate::valence::update_source_atom_cache(&mut graph, AtomId::new(i), false).unwrap();
        }
        graph
    }
    #[test]
    fn source_distinguishes_molhash_pruning_dummy_and_exocyclic_double_bond() {
        let graph = carbonyl_ring();
        let rings = crate::fast_find_rings(&graph).unwrap();
        let make = || {
            (
                graph.clone(),
                CoordinateBlock::default(),
                MoleculeProperties::default(),
            )
        };
        let (g, c, p) = make();
        let bare = murcko_scaffold(g, c, p, Some(&rings), None).unwrap();
        assert_eq!(bare.topology.atoms.len(), 6);
        assert_eq!(bare.topology.atoms[0].explicit_hydrogens(), 2);
        assert_eq!(bare.mapping.atoms.old_to_new[6], None);
        let (g, c, p) = make();
        let net = net_scaffold(g, c, p, Some(&rings), None).unwrap();
        assert_eq!(net.topology.atoms.len(), 7);
        assert_eq!(net.topology.atoms[6].element(), Element::DUMMY);
        assert_eq!(net.topology.bonds[6].order(), BondOrder::Double);
        let (g, c, p) = make();
        let decomposed = murcko_decompose(g, c, p, Some(&rings)).unwrap();
        assert_eq!(decomposed.topology, graph);
        assert!(decomposed.valence.is_none());
        assert!(decomposed.rings.is_none());
    }
    #[test]
    fn empty_source_behavior_and_unprepared_hydrogen_error() {
        for kind in 0..3 {
            let g = TopologyBlock::default();
            let c = CoordinateBlock::default();
            let p = MoleculeProperties::default();
            let out = match kind {
                0 => murcko_scaffold(g, c, p, None, None),
                1 => net_scaffold(g, c, p, None, None),
                _ => murcko_decompose(g, c, p, None),
            }
            .unwrap();
            assert!(out.topology.atoms.is_empty());
            assert_eq!(out.mapping, TopologyMapping::identity(0, 0));
        }
        let raw = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        assert!(
            matches!(murcko_scaffold(raw, CoordinateBlock::default(), MoleculeProperties::default(), None, None),
            Err(ScaffoldError::MissingImplicitHydrogens { atom }) if atom == AtomId::new(0))
        );
    }
    #[test]
    fn floyd_warshall_retains_first_equal_distance_predecessor() {
        let atoms = (0..4)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let bonds = [(0, 1), (0, 2), (1, 3), (2, 3)]
            .into_iter()
            .enumerate()
            .map(|(i, (a, b))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect();
        let graph = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        let paths = predecessor_matrix(&graph).unwrap();
        assert_eq!(paths[3], 1);
        assert_eq!(paths[12], 1);
        assert_eq!(paths[0], -1);
    }
}
