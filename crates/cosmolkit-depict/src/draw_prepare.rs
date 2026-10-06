//! Detached legacy drawing preparation. Inherited source markers are not
//! newly certified; the dimension-specific 2D selection is packet-authorized.
#[cfg(test)]
#[path = "draw_prepare_stage_probe.rs"]
mod draw_prepare_stage_probe;
#[cfg(test)]
#[path = "drawing_state_probe.rs"]
mod drawing_state_probe;
#[cfg(test)]
#[path = "drawing_svg_boundary_probe.rs"]
mod drawing_svg_boundary_probe;
#[cfg(test)]
#[path = "draw_prepare_tests.rs"]
mod tests;
use crate::draw::PreparedDrawingInput;
use crate::{Compute2DCoordinatesParams, DrawingError, DrawingInput, compute_2d_coordinates};
use cosmolkit_core::{
    AddHsParams, AtropisomerConformer, KekulizeAttempt, KekulizeParams, RingInfo,
    ValenceAssignment, ValenceModel, WedgeInfo, add_hydrogens_with_params,
    assign_valence_with_options_for_topology, determine_bond_wedge_state,
    kekulize_if_possible_with_query_state_and_ring_info,
    pick_bonds_to_wedge_with_existing_ring_info,
};
use cosmolkit_model::{
    AtomId, BondDirection, BondId, BondOrder, ChiralTag, CoordinateBlock, MoleculeProperties,
    SdfPropertyListTarget, TopologyBlock, TopologyMapping,
};

pub(crate) struct PreparedDrawing {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
    pub valence: ValenceAssignment,
    pub rings: RingInfo,
}

impl PreparedDrawing {
    pub fn borrow(&self) -> PreparedDrawingInput<'_> {
        PreparedDrawingInput {
            topology: &self.topology,
            layout: &self.coordinates.conformers_2d[0],
            properties: &self.properties,
            valence: &self.valence,
            rings: &self.rings,
        }
    }
}

fn check_ring_dimensions(rings: &RingInfo, atoms: usize, bonds: usize) -> Result<(), DrawingError> {
    for (field, actual, expected) in [
        ("ring atom memberships", rings.atom_row_count(), atoms),
        ("ring bond memberships", rings.bond_row_count(), bonds),
    ] {
        if actual != expected {
            return Err(DrawingError::StateRows {
                field,
                actual,
                expected,
            });
        }
    }
    Ok(())
}

fn check_properties(
    properties: &MoleculeProperties,
    atoms: usize,
    bonds: usize,
) -> Result<(), DrawingError> {
    for list in properties.sdf_property_lists() {
        let expected_rows = match list.target() {
            SdfPropertyListTarget::Atom => atoms,
            SdfPropertyListTarget::Bond => bonds,
        };
        if list.values().len() != expected_rows {
            return Err(cosmolkit_core::HydrogenError::InvalidPropertyList {
                target: list.target(),
                name: list.name().to_owned(),
                expected_rows,
                actual_rows: list.values().len(),
            }
            .into());
        }
    }
    Ok(())
}

pub(crate) fn prepare(input: DrawingInput<'_>) -> Result<PreparedDrawing, DrawingError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/MolDraw2D/MolDraw2DUtils.cpp :: prepareMolForDrawing
    // RDKit✔️✔️: void prepareMolForDrawing(RWMol &mol, bool kekulize, bool addChiralHs,
    // RDKit✔️✔️:                           bool wedgeBonds, bool forceCoords, bool wavyBonds) {
    // RDKit✔️✔️:   if (kekulize) {
    // RDKit✔️✔️:     RDLog::LogStateSetter blocker;
    // RDKit✔️✔️:     MolOps::KekulizeIfPossible(
    // RDKit✔️✔️:         mol, false);  // kekulize, but keep the aromatic flags!
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (addChiralHs) {
    // RDKit✔️✔️:     std::vector<unsigned int> chiralAts;
    // RDKit✔️✔️:     for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:       if (isAtomCandForChiralH(mol, atom)) {
    // RDKit✔️✔️:         chiralAts.push_back(atom->getIdx());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (chiralAts.size()) {
    // RDKit✔️✔️:       bool addCoords = false;
    // RDKit✔️✔️:       if (!forceCoords && mol.getNumConformers()) {
    // RDKit✔️✔️:         addCoords = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       MolOps::addHs(mol, false, addCoords, &chiralAts);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (forceCoords || !mol.getNumConformers()) {
    // RDKit✔️✔️:     const bool canonOrient = true;
    // RDKit✔️✔️:     RDDepict::compute2DCoords(mol, nullptr, canonOrient);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (wedgeBonds) {
    // RDKit✔️✔️:     Chirality::wedgeMolBonds(mol, &mol.getConformer());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (wavyBonds) {
    // RDKit✔️✔️:     addWavyBondsForStereoAny(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/MolDraw2D/MolDraw2DUtils.cpp :: prepareMolForDrawing
    input.topology.validate()?;
    input
        .coordinates
        .validate_for_atom_count(input.topology.atoms.len())?;
    check_properties(
        input.properties,
        input.topology.atoms.len(),
        input.topology.bonds.len(),
    )?;
    if let Some(valence) = input.valence {
        for (field, actual) in [
            ("explicit valence", valence.explicit_valence.len()),
            ("implicit hydrogens", valence.implicit_hydrogens.len()),
        ] {
            if actual != input.topology.atoms.len() {
                return Err(DrawingError::StateRows {
                    field,
                    actual,
                    expected: input.topology.atoms.len(),
                });
            }
        }
    }
    if let Some(rings) = input.rings {
        check_ring_dimensions(
            rings,
            input.topology.atoms.len(),
            input.topology.bonds.len(),
        )?;
    }
    // One detached source copy; never install working state in a live owner.
    let mut topology = input.topology.clone();
    let mut coordinates = input.coordinates.clone();
    let mut properties = input.properties.clone();
    let mut rings = input.rings.cloned();
    #[cfg(test)]
    draw_prepare_stage_probe::observe(
        "entry",
        &topology,
        &coordinates,
        &properties,
        rings.as_ref(),
    );
    if !topology.bonds.is_empty() {
        match kekulize_if_possible_with_query_state_and_ring_info(
            &topology,
            &KekulizeParams {
                mark_atoms_bonds: false,
                ..Default::default()
            },
            None,
            rings.as_ref(),
        )? {
            KekulizeAttempt::Applied(assignment) => {
                topology = assignment.topology;
                if let Some(update) = assignment.ring_update {
                    rings = Some(update);
                }
            }
            // Existing source-defined sanitize fallback; failure-side ring
            // transport remains qualified by the accepted core prerequisite.
            KekulizeAttempt::NotKekulizable {
                topology: fallback, ..
            } => topology = fallback,
        }
    }
    #[cfg(test)]
    draw_prepare_stage_probe::observe(
        "post_kekulize",
        &topology,
        &coordinates,
        &properties,
        rings.as_ref(),
    );
    if let Some(carrier) = rings.as_ref() {
        check_ring_dimensions(carrier, topology.atoms.len(), topology.bonds.len())?;
    }
    let chiral_atoms = rings
        .as_ref()
        .filter(|r| r.is_initialized())
        .map(|r| {
            topology
                .atoms
                .iter()
                .filter(|a| {
                    r.num_atom_rings(a.id()) > 1
                        && matches!(
                            a.chiral_tag(),
                            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
                        )
                })
                .map(|a| a.id())
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    if !chiral_atoms.is_empty() {
        let old_atoms = topology.atoms.iter().map(|a| a.id()).collect::<Vec<_>>();
        let old_bonds = topology
            .bonds
            .iter()
            .map(|b| (b.id(), b.begin(), b.end()))
            .collect::<Vec<_>>();
        let add_coords =
            !coordinates.conformers_2d.is_empty() || !coordinates.conformers_3d.is_empty();
        let result = add_hydrogens_with_params(
            topology,
            coordinates,
            properties,
            &AddHsParams {
                add_coords,
                only_on_atoms: Some(chiral_atoms),
                ..Default::default()
            },
        )?;
        // Validate the ACTUAL result before extending the retained memberships:
        // mappings alone cannot establish that new rows are terminal hydrogens.
        check_hydrogen_append(&old_atoms, &old_bonds, &result.topology, &result.mapping)?;
        if let Some(carrier) = rings.as_mut() {
            check_ring_dimensions(carrier, old_atoms.len(), old_bonds.len())?;
            carrier.preallocate(result.topology.atoms.len(), result.topology.bonds.len());
        }
        topology = result.topology;
        coordinates = result.coordinates;
        properties = result.properties;
    }
    #[cfg(test)]
    draw_prepare_stage_probe::observe(
        "post_chiral_h",
        &topology,
        &coordinates,
        &properties,
        rings.as_ref(),
    );
    // Approved dimension-specific default, regardless of stored 3D. ID zero
    // belongs only to the newly generated 2D layout; retained IDs stay intact.
    if coordinates.conformers_2d.is_empty() {
        let conformer = compute_2d_coordinates(
            &topology,
            &properties,
            &Compute2DCoordinatesParams {
                canonical_orientation: true,
                ..Default::default()
            },
        )?;
        coordinates
            .record_source_conformer_append(cosmolkit_model::CoordinateDimension::TwoD)
            .map_err(DrawingError::Coordinates)?;
        coordinates.conformers_2d.push(conformer);
    }
    #[cfg(test)]
    draw_prepare_stage_probe::observe(
        "post_compute_2d",
        &topology,
        &coordinates,
        &properties,
        rings.as_ref(),
    );
    let rings = wedge(&mut topology, &coordinates.conformers_2d[0], rings)?;
    #[cfg(test)]
    draw_prepare_stage_probe::observe(
        "post_wedge",
        &topology,
        &coordinates,
        &properties,
        Some(&rings),
    );
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    check_properties(&properties, topology.atoms.len(), topology.bonds.len())?;
    check_ring_dimensions(&rings, topology.atoms.len(), topology.bonds.len())?;
    let valence =
        assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, false)?;
    Ok(PreparedDrawing {
        topology,
        coordinates,
        properties,
        valence,
        rings,
    })
}

fn check_hydrogen_append(
    old_atoms: &[AtomId],
    old_bonds: &[(BondId, AtomId, AtomId)],
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
) -> Result<(), DrawingError> {
    let fail = |row, reason| DrawingError::HydrogenAppend { row, reason };
    let (a, b) = (old_atoms.len(), old_bonds.len());
    let (new_a, new_b) = (topology.atoms.len(), topology.bonds.len());
    topology.validate()?;
    mapping.validate_for_counts(a, new_a, b, new_b)?;
    if new_a < a || new_b < b || new_a - a != new_b - b {
        return Err(fail(None, "atom/bond counts are not a paired append"));
    }
    for (i, &id) in old_atoms.iter().enumerate() {
        if topology.atoms[i].id() != id
            || mapping.atoms.old_to_new[i] != Some(id)
            || mapping.atoms.new_to_old[i] != Some(id)
        {
            return Err(fail(Some(i), "retained atom IDs/mappings changed"));
        }
    }
    for (i, &(id, begin, end)) in old_bonds.iter().enumerate() {
        let bond = &topology.bonds[i];
        if (bond.id(), bond.begin(), bond.end()) != (id, begin, end)
            || mapping.bonds.old_to_new[i] != Some(id)
            || mapping.bonds.new_to_old[i] != Some(id)
        {
            return Err(fail(
                Some(i),
                "retained bond IDs/endpoints/mappings changed",
            ));
        }
    }
    for i in a..new_a {
        if mapping.atoms.new_to_old[i].is_some()
            || topology.atoms[i].atomic_number() != 1
            || topology.adjacency.neighbors_of(i).len() != 1
        {
            return Err(fail(
                Some(i),
                "appended atom is not a new terminal hydrogen",
            ));
        }
        let neighbor = topology.adjacency.neighbors_of(i)[0];
        if neighbor.atom_index >= a || neighbor.bond.index() < b {
            return Err(fail(
                Some(i),
                "appended hydrogen is not bonded to one retained atom by a new bond",
            ));
        }
    }
    for i in b..new_b {
        let bond = &topology.bonds[i];
        if mapping.bonds.new_to_old[i].is_some()
            || bond.order() != BondOrder::Single
            || !((bond.begin().index() < a && bond.end().index() >= a)
                || (bond.end().index() < a && bond.begin().index() >= a))
        {
            return Err(fail(
                Some(i),
                "appended bond is not a new single bond to a terminal hydrogen",
            ));
        }
    }
    Ok(())
}

fn wedge(
    topology: &mut TopologyBlock,
    layout: &cosmolkit_model::Conformer2D,
    rings: Option<RingInfo>,
) -> Result<RingInfo, DrawingError> {
    // BEGIN RDKIT CPP CALL third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: outputMolToMolBlock
    // RDKit❗✔️: auto wedgeBonds = Chirality::pickBondsToWedge(tmol, nullptr, conf);
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: wedgeMolBonds
    // RDKit❗✔️: auto wedgeBonds = Chirality::pickBondsToWedge(mol, params, conf);
    // RDKit❗✔️: for (const auto &[wbi, wedgeInfo] : wedgeBonds) {
    // RDKit❗✔️:   auto bond = mol.getBondWithIdx(wbi);
    // RDKit❗✔️:   auto dir =
    // RDKit❗✔️:       detail::determineBondWedgeState(bond, wedgeInfo->getIdx(), conf);
    // RDKit❗✔️:   if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
    // RDKit❗✔️:     bond->setBondDir(dir);
    // RDKit❗✔️:     if (static_cast<unsigned int>(wedgeInfo->getIdx()) !=
    // RDKit❗✔️:         bond->getBeginAtomIdx()) {
    // RDKit❗✔️:       auto tmp = bond->getBeginAtomIdx();
    // RDKit❗✔️:       bond->setBeginAtomIdx(bond->getEndAtomIdx());
    // RDKit❗✔️:       bond->setEndAtomIdx(tmp);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: wedgeMolBonds
    let conformer = Some(AtropisomerConformer::TwoD(layout));
    let (assignments, rings) =
        pick_bonds_to_wedge_with_existing_ring_info(topology, conformer, rings)?;
    let mut updates = Vec::new();
    for (id, info) in assignments.iter() {
        match *info {
            WedgeInfo::Chiral { center } => {
                let direction = determine_bond_wedge_state(topology, id, center, conformer)?;
                if matches!(
                    direction,
                    BondDirection::BeginWedge | BondDirection::BeginDash
                ) {
                    let bond = &topology.bonds[id.index()];
                    let end = if bond.begin() == center {
                        bond.end()
                    } else {
                        bond.begin()
                    };
                    updates.push((id, direction, center, end));
                }
            }
            WedgeInfo::Atropisomer { update } => {
                updates.push((id, update.direction, update.begin, update.end))
            }
        }
    }
    for (id, direction, begin, end) in updates {
        let bond = &mut topology.bonds[id.index()];
        bond.set_direction(direction);
        bond.set_endpoints(begin, end);
    }
    // Reversing endpoint orientation leaves the undirected adjacency and its
    // incident-bond order intact. Validation checks that exact model invariant.
    topology.validate()?;
    Ok(rings)
}
