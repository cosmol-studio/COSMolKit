use super::*;

#[mol_op_body(with_hydrogens, parts)]
pub(super) fn with_hydrogens_impl(
    params: crate::hydrogens::AddHsParams,
) -> Result<(), OperationError> {
    // RDKit✔️✔️:   mol.clearComputedProps(false);
    parts.clear_computed_properties();
    parts.with_topology_coordinates_properties_mut(
        |parts, topology, coordinates, properties| {
            let assignment = parts.with_block_read_parts(
                topology.clone(),
                coordinates.clone(),
                properties.clone(),
                |read| {
                    read.add_hs_assignment(&params)
                        .map_err(|source| OperationError::AddHydrogens {
                            operation: &WITH_HYDROGENS_SPEC,
                            source,
                        })
                },
            )?;

            apply_add_hs_assignment(parts, topology, coordinates, properties, &assignment)
        },
    )?;
    parts.prove_preserved(
        DerivedState::RINGS | DerivedState::RING_FAMILIES,
        PreservationProof::LeafAtomAppend,
    )?;

    Ok(())
}

pub(super) fn apply_add_hs_assignment(
    parts: &mut OpParts<'_>,
    topology: &mut TopologyBlock,
    coordinates: &mut CoordinateBlock,
    properties: &mut MoleculeProperties,
    assignment: &crate::hydrogens::AddHsAssignment,
) -> Result<bool, OperationError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::addHs(RWMol&, const AddHsParameters&, const UINT_VECT*)
    if assignment.clear_computed_properties {
        // RDKit✔️✔️:   RDProps::clearComputedProps();
        crate::chemistry::cip::clear_computed_source_properties(topology, properties);
    }
    // RDKit❗✔️:   ...
    // RDKit❗✔️:     newAt->clearComputedProps();
    // END RDKIT CPP FUNCTION MolOps::addHs(RWMol&, const AddHsParameters&, const UINT_VECT*)

    let old_atom_count = topology.atoms.len();
    let old_bond_count = topology.bonds.len();
    let added_atom_count = assignment.hydrogens_to_add.len();
    let added_bond_count = assignment.hydrogens_to_add.len();
    let changed = added_atom_count != 0
        || !assignment.atom_explicit_hydrogen_updates.is_empty()
        || !assignment.atom_pdb_residue_info_updates.is_empty()
        || !assignment.clear_isotopic_hydrogen_properties.is_empty();

    let (coords_2d_to_append, conformer_coords_to_append) = parts.with_block_read_parts(
        topology.clone(),
        coordinates.clone(),
        properties.clone(),
        |read| {
            let coords_2d_to_append = read
                .coordinates_2d()
                .map(|coords| add_hs_terminal_coords_2d(read, assignment, coords));
            let coords_2d_to_append = match coords_2d_to_append {
                Some(coords) => Some(coords?),
                None => None,
            };
            let conformer_coords_to_append = read
                .conformers_3d()
                .iter()
                .map(|conformer| {
                    add_hs_terminal_coords_3d(read, assignment, conformer.coordinates())
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok((coords_2d_to_append, conformer_coords_to_append))
        },
    )?;

    for update in &assignment.atom_explicit_hydrogen_updates {
        let Some(atom) = topology.atoms.get_mut(update.atom.index()) else {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs explicit hydrogen update references an out-of-range atom",
            });
        };
        atom.set_explicit_hydrogens(update.explicit_hydrogens);
    }
    for atom_id in &assignment.clear_isotopic_hydrogen_properties {
        let Some(atom) = topology.atoms.get_mut(atom_id.index()) else {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs isotope property cleanup references an out-of-range atom",
            });
        };
        atom.set_tracked_isotopic_hydrogens(Vec::new());
    }
    for update in &assignment.atom_pdb_residue_info_updates {
        let Some(atom) = topology.atoms.get_mut(update.atom.index()) else {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs PDB residue info update references an out-of-range atom",
            });
        };
        atom.set_pdb_residue_info(Some(update.pdb_residue_info.clone()));
    }

    for hydrogen in &assignment.hydrogens_to_add {
        if topology.atoms.get(hydrogen.heavy_atom.index()).is_none() {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs hydrogen references an out-of-range heavy atom",
            });
        }
        let atom_id = crate::AtomId::new(topology.atoms.len());
        let mut spec = crate::AtomSpec::new(crate::Element::H);
        if let Some(isotope) = hydrogen.isotope {
            spec = spec.with_isotope(isotope);
        }
        if hydrogen.is_implicit {
            spec = spec.with_implicit_hydrogen(true);
        }
        for (key, value) in &hydrogen.props {
            spec = spec.with_prop(key, value);
        }
        if let Some(info) = hydrogen.pdb_residue_info.clone() {
            spec = spec.with_pdb_residue_info(info);
        }
        topology.atoms.push(crate::Atom::from_spec(atom_id, spec));

        let bond_id = crate::BondId::new(topology.bonds.len());
        topology.bonds.push(crate::Bond::from_spec(
            bond_id,
            crate::BondSpec::new(hydrogen.heavy_atom, atom_id, crate::BondOrder::Single),
        ));
    }

    if let Some(coords) = coords_2d_to_append {
        let Some(existing) = coordinates.conformers_2d.first_mut() else {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs lost 2D coordinate block while appending hydrogen coordinates",
            });
        };
        for coord in coords {
            existing.push_coord(coord);
        }
    }
    if !conformer_coords_to_append.is_empty() {
        for (conformer, coords) in coordinates
            .conformers_3d
            .iter_mut()
            .zip(conformer_coords_to_append)
        {
            for coord in coords {
                conformer.push_coord(coord);
            }
        }
    }

    let mapping = TopologyMapping::with_appended(
        old_atom_count,
        old_bond_count,
        added_atom_count,
        added_bond_count,
    );
    properties.remap_topology(&mapping);

    let adjacency = crate::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
    let valence = crate::hydrogens::add_hs_valence_assignment_from_parts(
        &topology.atoms,
        &topology.bonds,
        &adjacency,
        old_atom_count,
        assignment,
    )
    .map_err(|source| OperationError::AddHydrogens {
        operation: &WITH_HYDROGENS_SPEC,
        source,
    })?;

    parts.record_topology_edit(TopologyEditKind::Appending)?;
    parts.record_topology_mapping(mapping);
    parts.clear_cache(WITH_HYDROGENS_SPEC.needs_update());
    parts.set_valence_cache(valence);
    Ok(changed)
}

pub(super) fn add_hs_terminal_coords_2d(
    read_parts: MoleculeReadParts<'_>,
    assignment: &crate::hydrogens::AddHsAssignment,
    coords: &[[f64; 2]],
) -> Result<Vec<[f64; 2]>, OperationError> {
    let molecule = read_parts;
    let coords_3d = coords
        .iter()
        .map(|coord| [coord[0], coord[1], 0.0])
        .collect::<Vec<_>>();
    Ok(
        add_hs_terminal_coords(molecule, assignment, &coords_3d, false)?
            .into_iter()
            .map(|coord| [coord[0], coord[1]])
            .collect(),
    )
}

pub(super) fn add_hs_terminal_coords_3d(
    read_parts: MoleculeReadParts<'_>,
    assignment: &crate::hydrogens::AddHsAssignment,
    coords: &[[f64; 3]],
) -> Result<Vec<[f64; 3]>, OperationError> {
    let molecule = read_parts;
    add_hs_terminal_coords(molecule, assignment, coords, true)
}

fn add_hs_terminal_coords(
    molecule: MoleculeReadParts<'_>,
    assignment: &crate::hydrogens::AddHsAssignment,
    coords: &[[f64; 3]],
    is_3d: bool,
) -> Result<Vec<[f64; 3]>, OperationError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: void setTerminalAtomCoords(ROMol &mol, unsigned int idx, unsigned int otherIdx)
    // RDKit✔️✔️: void setTerminalAtomCoords(ROMol &mol, unsigned int idx,
    // RDKit✔️✔️:                            unsigned int otherIdx) {
    // RDKit✔️✔️:   PRECONDITION(otherIdx != idx, "degenerate atoms");
    // RDKit✔️✔️:   Atom *atom = mol.getAtomWithIdx(idx);
    // RDKit✔️✔️:   PRECONDITION(mol.getAtomDegree(atom) == 1, "bad atom degree");
    // RDKit✔️✔️:   const Bond *bond = mol.getBondBetweenAtoms(otherIdx, idx);
    // RDKit✔️✔️:   PRECONDITION(bond, "no bond between atoms");
    // RDKit✔️✔️:   const Atom *otherAtom = mol.getAtomWithIdx(otherIdx);
    // RDKit✔️✔️:   double bondLength =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getRb0(1) +
    // RDKit✔️✔️:       PeriodicTable::getTable()->getRb0(otherAtom->getAtomicNum());
    // RDKit✔️✔️:   RDGeom::Point3D dirVect(0, 0, 0);
    // RDKit✔️✔️:   switch (otherAtom->getDegree()) {
    // COSMolKit computes the same sequential coordinate assignment here, then
    // writes the finished coordinates through OpParts after topology mutation.
    let old_atom_count = molecule.num_atoms();
    let mut working_coords = coords.to_vec();
    let mut adjacency = add_hs_virtual_adjacency(molecule, assignment.hydrogens_to_add.len());
    let mut appended = Vec::with_capacity(assignment.hydrogens_to_add.len());

    for (offset, hydrogen) in assignment.hydrogens_to_add.iter().enumerate() {
        let heavy_index = hydrogen.heavy_atom.index();
        let hydrogen_index = old_atom_count + offset;
        if coords.get(heavy_index).is_none() {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs coordinate append references an out-of-range heavy atom",
            });
        }
        adjacency[heavy_index].push((hydrogen_index, None));
        adjacency[hydrogen_index].push((heavy_index, None));

        let coord = if assignment.add_terminal_coordinates {
            add_hs_set_terminal_atom_coord(
                molecule,
                &adjacency,
                &working_coords,
                hydrogen_index,
                heavy_index,
                is_3d,
            )?
        } else {
            working_coords[heavy_index]
        };
        working_coords.push(coord);
        appended.push(coord);
    }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: void setTerminalAtomCoords(ROMol &mol, unsigned int idx, unsigned int otherIdx)
    Ok(appended)
}

fn add_hs_virtual_adjacency(
    molecule: MoleculeReadParts<'_>,
    hydrogens_to_add: usize,
) -> Vec<Vec<(usize, Option<crate::BondId>)>> {
    let mut adjacency = vec![Vec::new(); molecule.num_atoms() + hydrogens_to_add];
    for bond in molecule.bonds() {
        adjacency[bond.begin().index()].push((bond.end().index(), Some(bond.id())));
        adjacency[bond.end().index()].push((bond.begin().index(), Some(bond.id())));
    }
    adjacency
}

pub(crate) fn add_hs_set_terminal_atom_coord(
    molecule: MoleculeReadParts<'_>,
    adjacency: &[Vec<(usize, Option<crate::BondId>)>],
    coords: &[[f64; 3]],
    atom_index: usize,
    other_index: usize,
    is_3d: bool,
) -> Result<[f64; 3], OperationError> {
    if atom_index == other_index {
        return Err(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs terminal coordinate assignment received degenerate atoms",
        });
    }
    let atom_neighbors = adjacency
        .get(atom_index)
        .ok_or(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs terminal coordinate assignment references an out-of-range atom",
        })?;
    if atom_neighbors.len() != 1 {
        return Err(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs terminal coordinate assignment requires terminal degree one",
        });
    }
    if atom_neighbors[0].0 != other_index {
        return Err(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs terminal coordinate assignment requires a bond to the heavy atom",
        });
    }
    let other_atom = molecule
        .atoms()
        .get(other_index)
        .ok_or(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs coordinate append references an out-of-range heavy atom",
        })?;
    let other_pos = *coords
        .get(other_index)
        .ok_or(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs coordinate append references an out-of-range heavy atom",
        })?;
    let bond_length = if is_3d {
        add_hs_rdkit_rb0(1) + add_hs_rdkit_rb0(other_atom.atomic_number())
    } else {
        1.0
    };
    let neighbors = adjacency
        .get(other_index)
        .ok_or(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs coordinate append references an out-of-range heavy atom",
        })?
        .iter()
        .filter_map(|(neighbor, _)| (*neighbor != atom_index).then_some(*neighbor))
        .collect::<Vec<_>>();

    let mut dir: [f64; 3];
    match neighbors.len() + 1 {
        1 => {
            // RDKit✔️✔️:     case 1:
            // RDKit✔️✔️:       // No other atoms present.
            // RDKit✔️✔️: if ((*cfi)->is3D()) { dirVect.z = 1; } else { dirVect.x = 1; }
            dir = if is_3d {
                [0.0, 0.0, 1.0]
            } else {
                [1.0, 0.0, 0.0]
            };
        }
        2 => {
            // RDKit✔️✔️:     case 2:
            // RDKit✔️✔️:       // One other neighbor.
            // RDKit✔️✔️: nbr1Vect = nbr1Pos - otherPos; ... nbr1Vect *= -1;
            let nbr1 = neighbors[0];
            let Some(nbr1_vec) = add_hs_direction_away(coords, other_pos, nbr1) else {
                return Ok(other_pos);
            };
            match other_atom.hybridization() {
                crate::Hybridization::Sp3 => {
                    // RDKit✔️✔️: tform.SetRotation((180 - 109.471) * M_PI / 180., perpVect);
                    let perp = if is_3d {
                        add_hs_perpendicular(nbr1_vec)
                    } else {
                        [0.0, 0.0, 1.0]
                    };
                    dir = add_hs_rotate(nbr1_vec, perp, (180.0 - 109.471_f64).to_radians());
                }
                crate::Hybridization::Sp2 => {
                    // RDKit✔️✔️: default 3D position is an arbitrary perpendicular; for 2D z-normal.
                    let mut perp = if is_3d {
                        add_hs_perpendicular(nbr1_vec)
                    } else {
                        [0.0, 0.0, 1.0]
                    };
                    if adjacency[nbr1].len() > 1
                        && add_hs_bond_is_pi_like(molecule, adjacency, other_index, nbr1)
                        && let Some(nbr2) = adjacency[nbr1].iter().find_map(|(neighbor, _)| {
                            (*neighbor != other_index).then_some(*neighbor)
                        })
                        && let (Some(nbr1_pos), Some(nbr2_pos)) =
                            (coords.get(nbr1).copied(), coords.get(nbr2).copied())
                        && let Some(nbr2_vec) =
                            add_hs_normalized(add_hs_vec_sub(nbr2_pos, nbr1_pos))
                    {
                        let cross = add_hs_cross(nbr2_vec, nbr1_vec);
                        if add_hs_len_sq(cross) >= ADD_HS_SQ_DIST_ZERO_TOL {
                            perp = cross;
                        }
                    }
                    perp = add_hs_normalized(perp).unwrap_or(perp);
                    dir = add_hs_rotate(nbr1_vec, perp, 60.0_f64.to_radians());
                }
                crate::Hybridization::Sp => {
                    // RDKit✔️✔️: just lay the H along the vector.
                    dir = nbr1_vec;
                }
                _ => {
                    // RDKit✔️✔️: default branch lays the H along the vector.
                    dir = nbr1_vec;
                }
            }
        }
        3 => {
            // RDKit✔️✔️:     case 3:
            // RDKit✔️✔️:       // Two other neighbors.
            // RDKit✔️✔️: start along the average of the two vectors.
            let Some(nbr1_vec) = add_hs_direction_away(coords, other_pos, neighbors[0]) else {
                return Ok(other_pos);
            };
            let Some(nbr2_vec) = add_hs_direction_away(coords, other_pos, neighbors[1]) else {
                return Ok(other_pos);
            };
            let Some(mut avg_dir) = add_hs_normalized(add_hs_vec_add(nbr1_vec, nbr2_vec)) else {
                return Ok(other_pos);
            };
            if is_3d {
                match other_atom.hybridization() {
                    crate::Hybridization::Sp3 => {
                        // RDKit✔️✔️: rotate about the perpendicular-to-neighbors axis.
                        let nbr_perp = add_hs_cross(nbr1_vec, nbr2_vec);
                        let rotn_axis = add_hs_cross(nbr_perp, avg_dir);
                        if let Some(axis) = add_hs_normalized(rotn_axis) {
                            avg_dir =
                                add_hs_rotate(avg_dir, axis, (109.471_f64 / 2.0).to_radians());
                        }
                    }
                    crate::Hybridization::Sp2 => {}
                    _ => {}
                }
            }
            dir = avg_dir;
        }
        4 => {
            // RDKit✔️✔️:     case 4:
            // RDKit✔️✔️:       // Three other neighbors.
            // RDKit✔️✔️: three other neighbors.
            let Some(nbr1_vec) = add_hs_direction_away(coords, other_pos, neighbors[0]) else {
                return Ok(other_pos);
            };
            let Some(nbr2_vec) = add_hs_direction_away(coords, other_pos, neighbors[1]) else {
                return Ok(other_pos);
            };
            let Some(nbr3_vec) = add_hs_direction_away(coords, other_pos, neighbors[2]) else {
                return Ok(other_pos);
            };
            if is_3d {
                if add_hs_dot(nbr3_vec, add_hs_cross(nbr1_vec, nbr2_vec)).abs() < 0.1 {
                    dir = add_hs_cross(nbr1_vec, nbr2_vec);
                    if add_hs_len_sq(dir) < ADD_HS_SQ_DIST_ZERO_TOL {
                        dir = add_hs_vec_scale(add_hs_cross(nbr1_vec, nbr3_vec), -1.0);
                    }
                    if add_hs_len_sq(dir) < ADD_HS_SQ_DIST_ZERO_TOL {
                        dir = add_hs_cross(nbr2_vec, nbr3_vec);
                    }
                    if add_hs_len_sq(dir) < ADD_HS_SQ_DIST_ZERO_TOL {
                        return Ok(other_pos);
                    }
                    if other_atom.prop("_CIPCode").is_some() {
                        let v1 = add_hs_vec_sub(dir, nbr3_vec);
                        let v2 = add_hs_vec_sub(nbr1_vec, nbr3_vec);
                        let v3 = add_hs_vec_sub(nbr2_vec, nbr3_vec);
                        let vol = add_hs_dot(v1, add_hs_cross(v2, v3));
                        if (other_atom.chiral_tag() == crate::ChiralTag::TetrahedralCcw
                            && vol < 0.0)
                            || (other_atom.chiral_tag() == crate::ChiralTag::TetrahedralCw
                                && vol > 0.0)
                        {
                            dir = add_hs_vec_scale(dir, -1.0);
                        }
                    }
                } else {
                    dir = add_hs_vec_add(add_hs_vec_add(nbr1_vec, nbr2_vec), nbr3_vec);
                }
            } else {
                // RDKit✔️✔️: pick the bisector opposite the smallest accumulated outer angle.
                let angle12 = add_hs_angle(nbr1_vec, nbr2_vec);
                let angle13 = add_hs_angle(nbr1_vec, nbr3_vec);
                let angle23 = add_hs_angle(nbr2_vec, nbr3_vec);
                let accum1 = angle12 + angle13;
                let accum2 = angle12 + angle23;
                let accum3 = angle13 + angle23;
                dir = if accum1 <= accum2 && accum1 <= accum3 {
                    add_hs_pick_bisector(nbr1_vec, nbr2_vec, nbr3_vec)
                } else if accum2 <= accum1 && accum2 <= accum3 {
                    add_hs_pick_bisector(nbr2_vec, nbr1_vec, nbr3_vec)
                } else {
                    add_hs_pick_bisector(nbr3_vec, nbr1_vec, nbr2_vec)
                };
            }
        }
        _ => {
            // RDKit✔️✔️:     default:
            // RDKit✔️✔️:       atomPos = otherPos + dirVect * bondLength;
            // RDKit✔️✔️:       for (...) { (*cfi)->setAtomPos(idx, atomPos); }
            return Ok([0.0, 0.0, 0.0]);
        }
    }

    let dir = add_hs_normalized(dir).unwrap_or(dir);
    Ok(add_hs_vec_add(
        other_pos,
        add_hs_vec_scale(dir, bond_length),
    ))
}

const ADD_HS_SQ_DIST_ZERO_TOL: f64 = 1.0e-8;

fn add_hs_rdkit_rb0(atomic_number: u8) -> f64 {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getRb0 / atomicData::Rb0
    // RDKit✔️✔️: double Rb0() const { return rB0; }
    // RDKit✔️✔️: double getRb0(UINT atomicNumber) const { ... return byanum[atomicNumber].Rb0(); }
    const RB0: [f64; 113] = [
        0.0, 0.33, 0.7, 1.23, 0.9, 0.82, 0.77, 0.7, 0.66, 0.611, 0.7, 1.54, 1.36, 1.18, 0.937,
        0.89, 1.04, 0.997, 1.74, 2.03, 1.74, 1.44, 1.32, 1.22, 1.18, 1.17, 1.17, 1.16, 1.15, 1.17,
        1.25, 1.26, 1.188, 1.2, 1.17, 1.167, 1.91, 2.16, 1.91, 1.62, 1.45, 1.34, 1.3, 1.27, 1.25,
        1.25, 1.28, 1.34, 1.48, 1.44, 1.385, 1.4, 1.378, 1.387, 1.98, 2.35, 1.98, 1.69, 1.83, 1.82,
        1.81, 1.8, 1.8, 1.99, 1.79, 1.76, 1.75, 1.74, 1.73, 1.72, 1.94, 1.72, 1.44, 1.34, 1.3,
        1.28, 1.26, 1.27, 1.3, 1.34, 1.49, 1.48, 1.48, 1.45, 1.46, 1.45, 2.4, 2.0, 1.9, 1.88, 1.79,
        1.61, 1.58, 1.55, 1.53, 1.07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    RB0.get(usize::from(atomic_number)).copied().unwrap_or(0.0)
    // END RDKIT CPP FUNCTION PeriodicTable::getRb0 / atomicData::Rb0
}

fn add_hs_bond_is_pi_like(
    molecule: MoleculeReadParts<'_>,
    adjacency: &[Vec<(usize, Option<crate::BondId>)>],
    begin: usize,
    end: usize,
) -> bool {
    adjacency[begin]
        .iter()
        .find(|(neighbor, _)| *neighbor == end)
        .and_then(|(_, bond)| *bond)
        .and_then(|bond| molecule.bonds().get(bond.index()))
        .is_some_and(|bond| {
            bond.is_aromatic() || bond.order() == crate::BondOrder::Double || bond.is_conjugated()
        })
}

fn add_hs_direction_away(
    coords: &[[f64; 3]],
    other_pos: [f64; 3],
    neighbor: usize,
) -> Option<[f64; 3]> {
    let neighbor_pos = *coords.get(neighbor)?;
    add_hs_normalized(add_hs_vec_scale(
        add_hs_vec_sub(neighbor_pos, other_pos),
        -1.0,
    ))
}

fn add_hs_pick_bisector(nbr1_vec: [f64; 3], nbr2_vec: [f64; 3], nbr3_vec: [f64; 3]) -> [f64; 3] {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: RDGeom::Point3D pickBisector(const RDGeom::Point3D&, const RDGeom::Point3D&, const RDGeom::Point3D&)
    // RDKit✔️✔️: auto dirVect = nbr2Vect + nbr3Vect;
    let mut dir = add_hs_vec_add(nbr2_vec, nbr3_vec);
    if add_hs_len_sq(dir) < ADD_HS_SQ_DIST_ZERO_TOL {
        // RDKit✔️✔️: std::swap(dirVect.x, dirVect.y); dirVect.x *= -1;
        dir = [-nbr2_vec[1], nbr2_vec[0], nbr2_vec[2]];
    }
    if add_hs_dot(dir, nbr1_vec) < 0.0 {
        dir = add_hs_vec_scale(dir, -1.0);
    }
    dir
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: RDGeom::Point3D pickBisector(const RDGeom::Point3D&, const RDGeom::Point3D&, const RDGeom::Point3D&)
}

fn add_hs_vec_add(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
}

fn add_hs_vec_sub(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn add_hs_vec_scale(vec: [f64; 3], scale: f64) -> [f64; 3] {
    [vec[0] * scale, vec[1] * scale, vec[2] * scale]
}

fn add_hs_dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn add_hs_cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn add_hs_len_sq(vec: [f64; 3]) -> f64 {
    add_hs_dot(vec, vec)
}

fn add_hs_normalized(vec: [f64; 3]) -> Option<[f64; 3]> {
    let len_sq = add_hs_len_sq(vec);
    if len_sq < ADD_HS_SQ_DIST_ZERO_TOL {
        return None;
    }
    Some(add_hs_vec_scale(vec, 1.0 / len_sq.sqrt()))
}

fn add_hs_perpendicular(vec: [f64; 3]) -> [f64; 3] {
    let axis = if vec[0].abs() < vec[1].abs() {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    add_hs_normalized(add_hs_cross(vec, axis)).unwrap_or([0.0, 0.0, 1.0])
}

fn add_hs_rotate(vec: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let axis = add_hs_normalized(axis).unwrap_or([0.0, 0.0, 1.0]);
    let cos = angle.cos();
    let sin = angle.sin();
    let term1 = add_hs_vec_scale(vec, cos);
    let term2 = add_hs_vec_scale(add_hs_cross(axis, vec), sin);
    let term3 = add_hs_vec_scale(axis, add_hs_dot(axis, vec) * (1.0 - cos));
    add_hs_vec_add(add_hs_vec_add(term1, term2), term3)
}

fn add_hs_angle(left: [f64; 3], right: [f64; 3]) -> f64 {
    let denom = (add_hs_len_sq(left) * add_hs_len_sq(right)).sqrt();
    if denom < ADD_HS_SQ_DIST_ZERO_TOL {
        return 0.0;
    }
    (add_hs_dot(left, right) / denom).clamp(-1.0, 1.0).acos()
}

#[mol_op_body(without_hydrogens, parts)]
pub(super) fn without_hydrogens_impl(sanitize: bool) -> Result<(), OperationError> {
    let params = crate::hydrogens::RemoveHsParams::default();
    without_hydrogens_apply(parts, &params, sanitize, &WITHOUT_HYDROGENS_SPEC)
}

#[mol_op_body(without_hydrogens_with_params, parts)]
pub(super) fn without_hydrogens_with_params_impl(
    params: crate::hydrogens::RemoveHsParams,
    sanitize: bool,
) -> Result<(), OperationError> {
    without_hydrogens_apply(
        parts,
        &params,
        sanitize,
        &WITHOUT_HYDROGENS_WITH_PARAMS_SPEC,
    )
}

fn without_hydrogens_apply(
    parts: &mut OpParts<'_>,
    params: &crate::hydrogens::RemoveHsParams,
    sanitize: bool,
    operation: &'static MoleculeOpSpec,
) -> Result<(), OperationError> {
    parts.with_topology_coordinates_properties_mut(
        |parts, topology, coordinates, properties| {
            let assignment = parts.with_block_read_parts(
                topology.clone(),
                coordinates.clone(),
                properties.clone(),
                |read| {
                    read.remove_hs_assignment(params, sanitize)
                        .map_err(|source| OperationError::RemoveHydrogens { operation, source })
                },
            )?;

            let mut changed = apply_remove_hs_assignment(topology, &assignment, operation)?;
            let atoms_to_remove = assignment.atoms_to_remove.clone();
            let mapping = if atoms_to_remove.is_empty() {
                let atom_count = topology.atoms.len();
                let bond_count = topology.bonds.len();
                parts.record_topology_edit(TopologyEditKind::Compacting)?;
                TopologyMapping::identity(atom_count, bond_count)
            } else {
                let mapping = topology.remove_atoms_with_mapping(&atoms_to_remove);
                coordinates.remap_topology(&mapping);
                properties.remap_topology(&mapping);
                parts.record_topology_edit(TopologyEditKind::Compacting)?;
                changed = true;
                mapping
            };
            let valence = remap_remove_hs_valence_cache(&assignment, &mapping, operation)?;
            parts.record_topology_mapping(mapping);
            parts.clear_cache(
                WITHOUT_HYDROGENS_SPEC
                    .derived_effects
                    .recompute()
                    .union(WITHOUT_HYDROGENS_SPEC.derived_effects.invalidate()),
            );
            parts.set_valence_cache(valence);
            // RDKit✔️✔️:   mol.clearComputedProps(true);
            crate::chemistry::cip::clear_computed_source_properties(topology, properties);
            parts.clear_computed_properties();
            if assignment.sanitize_after_removal {
                changed |= sanitize_after_remove_hs_removal(
                    parts,
                    topology,
                    coordinates,
                    properties,
                    operation,
                )?;
            }
            Ok(changed)
        },
    )?;

    Ok(())
}

fn remap_remove_hs_valence_cache(
    assignment: &crate::hydrogens::RemoveHsAssignment,
    mapping: &TopologyMapping,
    operation: &'static MoleculeOpSpec,
) -> Result<crate::ValenceAssignment, OperationError> {
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     atom->updatePropertyCache(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ... mol.removeAtom(atom, false); ...
    // RDKit✔️❌:   mol.clearComputedProps(true);
    // Atom valence is stored outside RDProps in RDKit. Compacting removal
    // therefore drops removed rows and retains every surviving row unchanged.
    let new_atom_count = mapping.atoms().new_to_old().len();
    let mut explicit_valence = vec![None; new_atom_count];
    let mut implicit_hydrogens = vec![None; new_atom_count];
    for update in &assignment.atom_valence_cache_updates {
        let Some(new_atom) = mapping
            .atoms()
            .old_to_new()
            .get(update.atom.index())
            .and_then(|mapped| *mapped)
        else {
            continue;
        };
        explicit_valence[new_atom.index()] = Some(update.explicit_valence);
        implicit_hydrogens[new_atom.index()] = Some(update.implicit_hydrogens);
    }
    let explicit_valence = explicit_valence
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or(OperationError::InvalidInput {
            operation,
            message: "remove-H valence migration did not cover every surviving atom",
        })?;
    let implicit_hydrogens = implicit_hydrogens
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or(OperationError::InvalidInput {
            operation,
            message: "remove-H implicit-H migration did not cover every surviving atom",
        })?;
    Ok(crate::ValenceAssignment {
        explicit_valence,
        implicit_hydrogens,
    })
}

fn apply_remove_hs_assignment(
    topology: &mut TopologyBlock,
    assignment: &crate::hydrogens::RemoveHsAssignment,
    operation: &'static MoleculeOpSpec,
) -> Result<bool, OperationError> {
    let has_local_updates = !assignment.atom_explicit_hydrogen_updates.is_empty()
        || !assignment.atom_chirality_inversions.is_empty()
        || !assignment.atom_property_updates.is_empty()
        || !assignment.bond_direction_updates.is_empty()
        || !assignment.bond_stereo_updates.is_empty()
        || !assignment.bond_stereo_atom_replacements.is_empty()
        || !assignment.sgroup_updates.is_empty();

    if has_local_updates {
        for update in &assignment.atom_explicit_hydrogen_updates {
            let Some(atom) = topology.atoms.get_mut(update.atom.index()) else {
                return Err(OperationError::InvalidInput {
                    operation,
                    message: "remove-H explicit hydrogen update references an out-of-range atom",
                });
            };
            atom.set_explicit_hydrogens(update.explicit_hydrogens);
        }
        for atom_id in &assignment.atom_chirality_inversions {
            let Some(atom) = topology.atoms.get_mut(atom_id.index()) else {
                return Err(OperationError::InvalidInput {
                    operation,
                    message: "remove-H chirality update references an out-of-range atom",
                });
            };
            atom.set_chiral_tag(match atom.chiral_tag() {
                crate::ChiralTag::TetrahedralCw => crate::ChiralTag::TetrahedralCcw,
                crate::ChiralTag::TetrahedralCcw => crate::ChiralTag::TetrahedralCw,
                other => other,
            });
        }
        for update in &assignment.atom_property_updates {
            match update {
                crate::hydrogens::AtomPropertyUpdate::SetUnknownStereo { atom } => {
                    let Some(atom) = topology.atoms.get_mut(atom.index()) else {
                        return Err(OperationError::InvalidInput {
                            operation,
                            message: "remove-H atom property update references an out-of-range atom",
                        });
                    };
                    atom.set_unknown_stereo(true);
                }
                crate::hydrogens::AtomPropertyUpdate::SetIsotopicHydrogens { atom, isotopes } => {
                    let Some(atom) = topology.atoms.get_mut(atom.index()) else {
                        return Err(OperationError::InvalidInput {
                            operation,
                            message: "remove-H isotope property update references an out-of-range atom",
                        });
                    };
                    atom.set_tracked_isotopic_hydrogens(isotopes.clone());
                }
                crate::hydrogens::AtomPropertyUpdate::ClearIsotopicHydrogens { atom } => {
                    let Some(atom) = topology.atoms.get_mut(atom.index()) else {
                        return Err(OperationError::InvalidInput {
                            operation,
                            message: "remove-H isotope property cleanup references an out-of-range atom",
                        });
                    };
                    atom.set_tracked_isotopic_hydrogens(Vec::new());
                }
                crate::hydrogens::AtomPropertyUpdate::ClearExcessChiralExplicitHydrogens {
                    atom,
                } => {
                    let Some(atom) = topology.atoms.get_mut(atom.index()) else {
                        return Err(OperationError::InvalidInput {
                            operation,
                            message: "remove-H chiral explicit-H cleanup references an out-of-range atom",
                        });
                    };
                    atom.set_explicit_hydrogens(0);
                }
            }
        }
        for update in &assignment.bond_direction_updates {
            let Some(bond) = topology.bonds.get_mut(update.bond.index()) else {
                return Err(OperationError::InvalidInput {
                    operation,
                    message: "remove-H bond direction update references an out-of-range bond",
                });
            };
            bond.set_direction(update.direction);
        }
        for update in &assignment.bond_stereo_updates {
            let Some(bond) = topology.bonds.get_mut(update.bond.index()) else {
                return Err(OperationError::InvalidInput {
                    operation,
                    message: "remove-H bond stereo update references an out-of-range bond",
                });
            };
            bond.set_stereo_atoms(update.stereo_atoms);
            bond.set_stereo(update.stereo);
        }
        for update in &assignment.bond_stereo_atom_replacements {
            let Some(bond) = topology.bonds.get_mut(update.bond.index()) else {
                return Err(OperationError::InvalidInput {
                    operation,
                    message: "remove-H stereo atom replacement references an out-of-range bond",
                });
            };
            let Some([left, right]) = bond.stereo_atoms() else {
                continue;
            };
            let stereo_atoms = [
                if left == update.old_atom {
                    update.new_atom
                } else {
                    left
                },
                if right == update.old_atom {
                    update.new_atom
                } else {
                    right
                },
            ];
            bond.set_stereo_atoms(Some(stereo_atoms));
        }
        for update in &assignment.sgroup_updates {
            for substance_group in &mut topology.substance_groups {
                match update {
                    crate::hydrogens::SGroupRemoveHsUpdate::RemoveBond { bond } => {
                        substance_group.remove_bond(*bond);
                    }
                    crate::hydrogens::SGroupRemoveHsUpdate::RemoveAtom { atom } => {
                        substance_group.remove_atom(*atom);
                    }
                    crate::hydrogens::SGroupRemoveHsUpdate::RemoveParentAtom { atom } => {
                        substance_group.remove_parent_atom(*atom);
                    }
                    crate::hydrogens::SGroupRemoveHsUpdate::ClearAttachPointLeavingAtom {
                        atom,
                    } => {
                        substance_group.clear_attach_point_leaving_atom(*atom);
                    }
                }
            }
        }
    }
    Ok(has_local_updates)
}

fn sanitize_after_remove_hs_removal(
    parts: &mut OpParts<'_>,
    topology: &mut TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    operation: &'static MoleculeOpSpec,
) -> Result<bool, OperationError> {
    // BEGIN RDKIT CPP BODY MolOps::removeHs sanitize-after-removal
    // RDKit✔️✔️:   if (!atomsToRemove.empty() && ps.removeNonimplicit && sanitize) {
    // RDKit✔️✔️:     sanitizeMol(mol);
    // RDKit✔️✔️:   }
    let changed = run_sanitize_pipeline_on_topology(
        parts,
        topology,
        Some(coordinates),
        Some(properties),
        crate::SanitizeOps::ALL,
        operation,
    )?;
    // END RDKIT CPP BODY MolOps::removeHs sanitize-after-removal
    Ok(changed)
}
