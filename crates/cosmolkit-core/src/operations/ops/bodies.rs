//! Operation bodies for registry entries declared in `ops.rs`.
//!
//! These functions intentionally live outside the runtime module that owns
//! `OpParts::working`, so they can only use the public capability methods.

use super::*;

#[mol_op_body(with_kekulized_bonds, parts)]
pub(super) fn with_kekulized_bonds_impl(clear_aromatic_flags: bool) -> Result<(), OperationError> {
    // RDKit✔️✔️: void kekulizeMol(ROMol &mol, bool clearAromaticFlags = false,
    // RDKit✔️✔️:                  bool canonical = true) {
    // RDKit✔️✔️:   auto &wmol = static_cast<RWMol &>(mol);
    // RDKit✔️✔️:   MolOps::Kekulize(wmol, clearAromaticFlags, canonical);
    // RDKit✔️✔️: }
    let valence = parts.with_topology_mut(|parts, topology| {
        let rings = parts.with_topology_read_parts(topology.clone(), |read| {
            read.symmetrize_sssr().map_err(|source| OperationError::RingFinding {
                operation: &WITH_KEKULIZED_BONDS_SPEC,
                source,
            })
        })?;
        parts.set_rings_cache(rings);
        let assignment = parts.with_topology_read_parts(topology.clone(), |read| {
            let ring_info = read
                .derived_cache()
                .rings
                .as_ref()
                .expect("rings were recomputed immediately above")
                .clone();
            read.kekulize_assignment(Some(&ring_info), clear_aromatic_flags, true, 100)
                .map_err(|source| OperationError::Kekulize {
                    operation: &WITH_KEKULIZED_BONDS_SPEC,
                    source,
                })
        })?;

        crate::kekulize::apply_kekulize_assignment(topology, &assignment);
        let valence = parts.with_topology_read_parts(topology.clone(), |read| {
            read.assign_valence_with_options(crate::ValenceModel::RdkitLike, true)
                .map_err(|source| OperationError::Valence {
                    operation: &WITH_KEKULIZED_BONDS_SPEC,
                    source,
                })
        })?;
        Ok(valence)
    })?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.clear_cache(DerivedState::AROMATICITY);
    parts.set_valence_cache(valence);
    parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
    Ok(())
}

#[mol_op_body(assigned_valence, parts)]
pub(super) fn assigned_valence_impl(strict: bool) -> Result<(), OperationError> {
    let read = parts.begin_topology_read()?;
    let valence = read
        .assign_valence_with_options(crate::ValenceModel::RdkitLike, strict)
        .map_err(|source| OperationError::Valence {
            operation: &ASSIGNED_VALENCE_SPEC,
            source,
        })?;
    parts.set_valence_cache(valence);
    Ok(())
}

#[mol_op_body(assigned_rings, parts)]
pub(super) fn assigned_rings_impl() -> Result<(), OperationError> {
    let read = parts.begin_topology_read()?;
    let rings = read.symmetrize_sssr().map_err(|source| OperationError::RingFinding {
        operation: &ASSIGNED_RINGS_SPEC,
        source,
    })?;
    parts.set_rings_cache(rings);
    Ok(())
}

#[mol_op_body(assigned_ring_families, parts)]
pub(super) fn assigned_ring_families_impl() -> Result<(), OperationError> {
    let read = parts.begin_topology_read()?;
    let ring_families = read
        .find_ring_families(false, false)
        .map_err(|source| OperationError::RingFinding {
            operation: &ASSIGNED_RING_FAMILIES_SPEC,
            source,
        })?;
    parts.set_ring_families_cache(ring_families);
    Ok(())
}

#[mol_op_body(assigned_aromaticity, parts)]
pub(super) fn assigned_aromaticity_impl() -> Result<(), OperationError> {
    let valence = parts.with_topology_mut(|parts, topology| {
        let rings = parts.with_topology_read_parts(topology.clone(), |read| {
            read.symmetrize_sssr().map_err(|source| OperationError::RingFinding {
                operation: &ASSIGNED_AROMATICITY_SPEC,
                source,
            })
        })?;
        parts.set_rings_cache(rings);
        let assignment = parts.with_topology_read_parts(topology.clone(), |read| {
            read.set_aromaticity(crate::AromaticityModel::Default)
                .map_err(|source| OperationError::Aromaticity {
                    operation: &ASSIGNED_AROMATICITY_SPEC,
                    source,
                })
        })?;
        for (atom, is_aromatic) in topology.atoms.iter_mut().zip(assignment.atom_aromatic.iter().copied()) {
            atom.set_aromatic(is_aromatic);
        }
        for (bond, is_aromatic) in topology.bonds.iter_mut().zip(assignment.bond_aromatic.iter().copied()) {
            bond.set_aromatic(is_aromatic);
            if is_aromatic && matches!(bond.order(), crate::BondOrder::Single | crate::BondOrder::Double) {
                bond.set_order(crate::BondOrder::Aromatic);
            }
        }
        parts.with_topology_read_parts(topology.clone(), |read| {
            read.assign_valence_with_options(crate::ValenceModel::RdkitLike, true)
                .map_err(|source| OperationError::Valence {
                    operation: &ASSIGNED_AROMATICITY_SPEC,
                    source,
                })
        })
    })?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.set_valence_cache(valence);
    parts.mark_aromaticity_valid();
    parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
    Ok(())
}

#[mol_op_body(assigned_radicals, parts)]
pub(super) fn assigned_radicals_impl() -> Result<(), OperationError> {
    let valence = parts.with_topology_mut(|parts, topology| {
        let (radicals, changed) = parts.with_topology_read_parts(topology.clone(), |read| {
            let radicals = read.assign_radicals().map_err(|source| OperationError::Valence {
                operation: &ASSIGNED_RADICALS_SPEC,
                source,
            })?;
            let changed = read
                .atoms()
                .iter()
                .zip(radicals.iter().copied())
                .any(|(atom, radical)| atom.radical_electrons() != radical);
            Ok((radicals, changed))
        })?;

        if changed {
            for (atom, radical) in topology.atoms.iter_mut().zip(radicals) {
                atom.set_radical_electrons(radical);
            }
        }

        parts.with_topology_read_parts(topology.clone(), |read| {
            read.assign_valence_with_options(crate::ValenceModel::RdkitLike, true)
                .map_err(|source| OperationError::Valence {
                    operation: &ASSIGNED_RADICALS_SPEC,
                    source,
                })
        })
    })?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.set_valence_cache(valence);
    Ok(())
}

#[mol_op_body(with_chiral_tags_from_structure, parts)]
pub(super) fn with_chiral_tags_from_structure_impl(
    conformer_id: i32,
    replace_existing_tags: bool,
) -> Result<(), OperationError> {
    let (original_topology, original_properties, selected_conformer, implicit_hydrogens) = {
        let read = parts.begin_operation_read()?;
        let Some(selected_conformer) = (if read.conformers_3d().is_empty() {
            None
        } else if conformer_id < 0 {
            read.conformers_3d().first()
        } else {
            read.conformers_3d()
                .iter()
                .find(|conformer| conformer.id() == conformer_id as usize)
        }) else {
            if read.conformers_3d().is_empty() {
                return Ok(());
            }
            return Err(OperationError::Stereo {
                operation: &WITH_CHIRAL_TAGS_FROM_STRUCTURE_SPEC,
                source: crate::StereoError::ConformerNotFound { conformer_id },
            });
        };
        if !selected_conformer.is_3d() {
            return Ok(());
        }
        (
            read.topology().clone(),
            read.properties().clone(),
            selected_conformer.clone(),
            read.derived_cache()
                .valence
                .as_ref()
                .map(|valence| valence.implicit_hydrogens.clone()),
        )
    };

    let mut topology = original_topology.clone();
    let mut properties = original_properties.clone();
    crate::chemistry::stereo::assign_chiral_types_from_3d_kernel(
        &mut topology.atoms,
        &topology.bonds,
        &topology.adjacency,
        std::slice::from_ref(&selected_conformer),
        &mut properties,
        implicit_hydrogens.as_deref(),
        conformer_id,
        replace_existing_tags,
    )
    .map_err(|source| OperationError::Stereo {
        operation: &WITH_CHIRAL_TAGS_FROM_STRUCTURE_SPEC,
        source,
    })?;

    if topology == original_topology && properties == original_properties {
        return Ok(());
    }

    parts.with_topology_and_properties_mut(|_parts, working_topology, working_properties| {
        *working_topology = topology;
        *working_properties = properties;
        Ok(())
    })?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.clear_cache(DerivedState::STEREO | DerivedState::DRAWING | DerivedState::FINGERPRINT);
    Ok(())
}

#[mol_op_body(assigned_cip_labels, parts)]
pub(super) fn assigned_cip_labels_impl(options: crate::CipLabelOptions) -> Result<(), OperationError> {
    let assignment_error = parts.with_topology_and_properties_mut(|parts, topology, properties| {
        let transition = parts.with_borrowed_optional_block_read_parts(topology, None, Some(properties), |read| {
            crate::chemistry::ciplabeler::assign_cip_labels_from_read_parts(read, &options).map_err(|source| {
                OperationError::CipLabeler {
                    operation: &ASSIGNED_CIP_LABELS_SPEC,
                    source,
                }
            })
        })?;
        let (assigned_topology, assigned_properties, assignment_error) = transition.into_parts();
        *topology = assigned_topology;
        *properties = assigned_properties;
        Ok(assignment_error)
    })?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.clear_cache(DerivedState::STEREO | DerivedState::DRAWING | DerivedState::FINGERPRINT);
    match assignment_error {
        Some(source) => Err(OperationError::CipLabeler {
            operation: &ASSIGNED_CIP_LABELS_SPEC,
            source,
        }),
        None => Ok(()),
    }
}

#[mol_op_body(with_2d_coordinates, parts)]
pub(super) fn with_2d_coordinates_impl(params: crate::With2DCoordinatesParams) -> Result<(), OperationError> {
    let coords = parts.with_coordinate_update_read_parts(|read| {
        crate::coordinates::compute_2d_coords_with_properties_and_params(
            read.atoms(),
            read.bonds(),
            read.properties(),
            &params.as_compute_params(),
        )
        .map_err(|source| match source {
            crate::coordinates::Coordinate2DError::InvalidInput(message) => OperationError::InvalidInput {
                operation: &WITH_2D_COORDINATES_SPEC,
                message,
            },
            crate::coordinates::Coordinate2DError::UnsupportedFeature(_) => OperationError::UnsupportedFeature {
                operation: &WITH_2D_COORDINATES_SPEC,
                source: crate::UnsupportedFeatureError::from_spec(&crate::COORDINATE_2D_FEATURE),
            },
        })
    })?;
    parts.with_coordinates_mut(|_parts, coord_block| {
        if params.clear_confs {
            coord_block.conformers_2d.clear();
        }
        coord_block
            .conformers_2d
            .push(crate::Conformer2D::new(coord_block.conformers_2d.len(), coords));
        coord_block.source_coordinate_dim = Some(crate::CoordinateDimension::TwoD);
        Ok(())
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(())
}

#[mol_op_body(with_2d_coordinate_block, parts)]
pub(super) fn with_2d_coordinate_block_impl(coords: Vec<[f64; 2]>) -> Result<(), OperationError> {
    let atom_count = {
        let read = parts.begin_topology_read()?;
        read.num_atoms()
    };
    if coords.len() != atom_count {
        return Err(OperationError::InvalidInput {
            operation: &WITH_2D_COORDINATE_BLOCK_SPEC,
            message: "2D coordinate row count mismatch",
        });
    }

    parts.with_coordinates_mut(|_parts, coord_block| {
        coord_block.conformers_2d.clear();
        coord_block.conformers_2d.push(crate::Conformer2D::new(0, coords));
        coord_block.source_coordinate_dim = Some(crate::CoordinateDimension::TwoD);
        Ok(())
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(())
}

fn source_coordinate_dim_for_block(coord_block: &CoordinateBlock) -> Option<crate::CoordinateDimension> {
    if coord_block.conformers_3d.iter().any(crate::Conformer3D::is_3d) {
        Some(crate::CoordinateDimension::ThreeD)
    } else if !coord_block.conformers_2d.is_empty() || !coord_block.conformers_3d.is_empty() {
        Some(crate::CoordinateDimension::TwoD)
    } else {
        None
    }
}

#[mol_op_body(with_3d_conformer, parts)]
pub(super) fn with_3d_conformer_impl(mut params: crate::EmbedParameters) -> Result<(), OperationError> {
    parts.with_coordinates_mut(|parts, coord_block| {
        parts.with_coordinate_update_read_parts(|read| {
            crate::distgeom::embed_molecule_coordinate_update(read, coord_block, &mut params).map_err(|source| {
                OperationError::DistanceGeometry {
                    operation: &WITH_3D_CONFORMER_SPEC,
                    source,
                }
            })
        })
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(())
}

#[mol_op_body(with_3d_coordinates, parts)]
pub(super) fn with_3d_coordinates_impl(coords: Vec<[f64; 3]>, conformer_index: usize) -> Result<(), OperationError> {
    let atom_count = {
        let read = parts.begin_topology_read()?;
        read.num_atoms()
    };
    if coords.len() != atom_count {
        return Err(OperationError::InvalidInput {
            operation: &WITH_3D_COORDINATES_SPEC,
            message: "3D conformer row count mismatch",
        });
    }

    parts.with_coordinates_mut(|_parts, coord_block| {
        if conformer_index >= coord_block.conformers_3d.len() {
            return Err(OperationError::InvalidInput {
                operation: &WITH_3D_COORDINATES_SPEC,
                message: "3D conformer index out of range",
            });
        }
        let existing = &coord_block.conformers_3d[conformer_index];
        coord_block.conformers_3d[conformer_index] = crate::Conformer3D::new(existing.id(), coords, existing.is_3d());
        coord_block.source_coordinate_dim = source_coordinate_dim_for_block(coord_block);
        Ok(())
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(())
}

#[mol_op_body(with_aligned_conformers, parts)]
pub(super) fn with_aligned_conformers_impl(
    params: crate::ConformerAlignmentParameters,
) -> Result<crate::ConformerAlignmentReport, OperationError> {
    let atom_count = parts.begin_topology_read()?.num_atoms();
    let rmsds = parts.with_coordinates_mut(|_parts, coordinates| {
        crate::mol_align::align_conformers_in_coordinate_block(coordinates, atom_count, &params).map_err(|source| {
            OperationError::Alignment {
                operation: &WITH_ALIGNED_CONFORMERS_SPEC,
                source,
            }
        })
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(crate::ConformerAlignmentReport { rmsds })
}

#[mol_op_body(with_alignment_to, parts)]
pub(super) fn with_alignment_to_impl(
    reference: &crate::Molecule,
    params: &crate::AlignmentParameters,
) -> Result<crate::AlignmentResult, OperationError> {
    let result = parts.with_coordinate_update_read_parts(|read| {
        crate::mol_align::alignment_result_from_read_parts(read, reference, params).map_err(|source| {
            OperationError::Alignment {
                operation: &WITH_ALIGNMENT_TO_SPEC,
                source,
            }
        })
    })?;
    parts.with_coordinates_mut(|_parts, coordinates| {
        crate::mol_align::apply_alignment_result_to_coordinate_block(coordinates, params.probe_conformer_id, &result)
            .map_err(|source| OperationError::Alignment {
                operation: &WITH_ALIGNMENT_TO_SPEC,
                source,
            })
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(result)
}

#[mol_op_body(with_atom_position, parts)]
pub(super) fn with_atom_position_impl(
    atom: usize,
    position: [f64; 3],
    conformer_index: usize,
) -> Result<(), OperationError> {
    let atom_count = parts.begin_topology_read()?.num_atoms();
    if atom >= atom_count {
        return Err(OperationError::MolTransform {
            operation: &WITH_ATOM_POSITION_SPEC,
            source: crate::mol_transforms::MolTransformError::AtomIndexOutOfBounds {
                atom,
                num_atoms: atom_count,
            },
        });
    }
    parts.with_coordinates_mut(|_parts, coordinates| {
        crate::mol_transforms::set_atom_position_in_coordinate_block(coordinates, atom, position, conformer_index)
            .map_err(|source| OperationError::MolTransform {
                operation: &WITH_ATOM_POSITION_SPEC,
                source,
            })
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(())
}

#[mol_op_body(with_cleared_3d_conformers, parts)]
pub(super) fn with_cleared_3d_conformers_impl() -> Result<(), OperationError> {
    parts.with_coordinates_mut(|_parts, coord_block| {
        coord_block.conformers_3d.clear();
        coord_block.source_coordinate_dim = source_coordinate_dim_for_block(coord_block);
        Ok(())
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(())
}

#[mol_op_body(with_3d_conformers, parts)]
pub(super) fn with_3d_conformers_impl(
    num_confs: usize,
    mut params: crate::EmbedParameters,
) -> Result<(), OperationError> {
    let num_confs = u32::try_from(num_confs).map_err(|_| OperationError::InvalidInput {
        operation: &WITH_3D_CONFORMERS_SPEC,
        message: "num_confs does not fit in RDKit unsigned int parameter",
    })?;
    parts.with_coordinates_mut(|parts, coord_block| {
        parts.with_coordinate_update_read_parts(|read| {
            crate::distgeom::embed_multiple_confs_coordinate_update(read, coord_block, num_confs, &mut params).map_err(
                |source| OperationError::DistanceGeometry {
                    operation: &WITH_3D_CONFORMERS_SPEC,
                    source,
                },
            )
        })
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(())
}

#[mol_op_body(with_added_3d_conformer, parts)]
pub(super) fn with_added_3d_conformer_impl(coords: Vec<[f64; 3]>, is_3d: bool) -> Result<(), OperationError> {
    let atom_count = {
        let read = parts.begin_topology_read()?;
        read.num_atoms()
    };
    if coords.len() != atom_count {
        return Err(OperationError::InvalidInput {
            operation: &WITH_ADDED_3D_CONFORMER_SPEC,
            message: "3D conformer row count mismatch",
        });
    }

    parts.with_coordinates_mut(|_parts, coord_block| {
        let next_id = coord_block
            .conformers_3d
            .iter()
            .map(crate::Conformer3D::id)
            .max()
            .map_or(0, |max_id| max_id + 1);
        coord_block
            .conformers_3d
            .push(crate::Conformer3D::new(next_id, coords, is_3d));
        coord_block.source_coordinate_dim = source_coordinate_dim_for_block(coord_block);
        Ok(())
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(())
}

#[mol_op_body(with_only_3d_conformer, parts)]
pub(super) fn with_only_3d_conformer_impl(coords: Vec<[f64; 3]>, is_3d: bool) -> Result<(), OperationError> {
    let atom_count = {
        let read = parts.begin_topology_read()?;
        read.num_atoms()
    };
    if coords.len() != atom_count {
        return Err(OperationError::InvalidInput {
            operation: &WITH_ONLY_3D_CONFORMER_SPEC,
            message: "3D conformer row count mismatch",
        });
    }

    parts.with_coordinates_mut(|_parts, coord_block| {
        coord_block.conformers_3d.clear();
        coord_block
            .conformers_3d
            .push(crate::Conformer3D::new(0, coords, is_3d));
        coord_block.source_coordinate_dim = source_coordinate_dim_for_block(coord_block);
        Ok(())
    })?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(())
}
