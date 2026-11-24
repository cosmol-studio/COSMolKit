//! Component chemistry belongs to core; runtime validates every mapped output.
use crate::{AtomId, DerivedState, Molecule, OperationError};
use cosmolkit_macros::{mol_multi_op_body, mol_op_body};

/// Options for extracting connected components with their source atom indices.
#[derive(Clone, Debug, PartialEq, Eq, cosmolkit_macros::Configuration)]
#[configuration(scalar)]
pub struct FragmentExtractionParams {
    /// Sanitize each fragment before returning it (default: true).
    pub sanitize: bool,
}

impl Default for FragmentExtractionParams {
    fn default() -> Self {
        Self { sanitize: true }
    }
}

/// A connected component and its atom correspondence at extraction time.
#[derive(Clone, Debug, PartialEq)]
pub struct MappedFragment {
    molecule: Molecule,
    source_atom_indices: Vec<AtomId>,
}

impl MappedFragment {
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }

    /// Entry `i` is the index of local atom `i` in the input to extraction.
    pub fn source_atom_indices(&self) -> &[AtomId] {
        &self.source_atom_indices
    }
}

#[mol_op_body(with_bonds_cut, parts)]
pub(crate) fn with_bonds_cut_impl(cuts: &[crate::BondCut]) -> Result<(), OperationError> {
    parts.transform_mapped_blocks(|topology, coordinates, properties| {
        let result = cosmolkit_core::cut_selected_bonds(topology, coordinates, properties, cuts)
            .map_err(OperationError::BondCut)?;
        Ok((
            result.topology,
            result.coordinates,
            result.properties,
            result.mapping,
        ))
    })?;
    parts.clear_cache(
        DerivedState::VALENCE
            .union(DerivedState::RINGS)
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.apply_cip_policy()
}

#[mol_multi_op_body(fragments, parts)]
pub(crate) fn fragments_impl() -> Result<(), OperationError> {
    let fragments = cosmolkit_core::get_molecule_fragments(
        parts.topology()?,
        parts.coordinates()?,
        parts.properties()?,
        true,
        true,
    )
    .map_err(OperationError::Fragments)?;
    parts.emit_mapped(
        fragments
            .into_iter()
            .map(|fragment| fragment.into_mapped_parts())
            .collect(),
    )
}

#[mol_multi_op_body(fragments_with_atom_mapping_with_params, parts)]
pub(crate) fn fragments_with_atom_mapping_impl(
    params: &FragmentExtractionParams,
) -> Result<Vec<Vec<AtomId>>, OperationError> {
    let fragments = cosmolkit_core::get_molecule_fragments(
        parts.topology()?,
        parts.coordinates()?,
        parts.properties()?,
        params.sanitize,
        true,
    )
    .map_err(OperationError::Fragments)?;

    let (rows, candidates) = fragments
        .into_iter()
        .map(|f| f.into_indexed_parts())
        .unzip();
    parts.emit_mapped(candidates)?;
    Ok(rows)
}

pub(crate) fn assemble_mapped_fragments(
    molecules: Vec<Molecule>,
    rows: Vec<Vec<AtomId>>,
) -> Result<Vec<MappedFragment>, OperationError> {
    let mismatch = |field, actual, expected| OperationError::InvalidAlgorithmResult {
        operation: "fragments_with_atom_mapping_with_params",
        field,
        actual,
        expected,
    };
    if rows.len() != molecules.len() {
        return Err(mismatch("component count", rows.len(), molecules.len()));
    }

    molecules
        .into_iter()
        .zip(rows)
        .map(|(molecule, source_atom_indices)| {
            if source_atom_indices.len() != molecule.num_atoms() {
                return Err(mismatch(
                    "source atom indices",
                    source_atom_indices.len(),
                    molecule.num_atoms(),
                ));
            }
            Ok(MappedFragment {
                molecule,
                source_atom_indices,
            })
        })
        .collect()
}

#[mol_op_body(largest_fragment, parts)]
pub(crate) fn largest_fragment_impl() -> Result<(), OperationError> {
    let (mapping, prepared) =
        parts.with_mutable_candidate_blocks(|topology, coordinates, properties, _cache| {
            let fragment =
                cosmolkit_core::get_largest_molecule_fragment(topology, coordinates, properties)
                    .map_err(OperationError::Fragments)?
                    .ok_or(OperationError::EmptyFragments)?;
            let (t, c, p, mapping, prepared) = fragment.into_mapped_parts();
            *topology = t;
            *coordinates = c;
            *properties = p;
            let has_prepared = prepared.is_some();
            if let Some((valence, rings)) = prepared {
                _cache.install_valence_assignment(valence);
                _cache.install_ring_info(rings);
            }
            Ok((mapping, has_prepared))
        })?;
    parts.record_topology_mapping(mapping)?;
    parts.clear_cache(
        DerivedState::RING_FAMILIES
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    if prepared {
        parts.mark_cache_updated(DerivedState::VALENCE.union(DerivedState::RINGS))?;
    } else {
        parts.clear_cache(DerivedState::VALENCE.union(DerivedState::RINGS))?;
    }
    parts.apply_cip_policy()
}

#[cfg(test)]
mod tests {
    use super::{AtomId, OperationError, assemble_mapped_fragments};
    use crate::{AtomSpec, Element, MoleculeBuilder};

    #[test]
    fn inconsistent_source_atom_rows_reject_the_entire_fragment_result() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        let carbon = builder.build().unwrap();
        let [a, b, c] = [0, 1, 2].map(AtomId::new);

        for (rows, field, actual, expected) in [
            (vec![vec![a]], "component count", 1, 2),
            (vec![vec![a], vec![b], vec![c]], "component count", 3, 2),
            (vec![vec![a], vec![]], "source atom indices", 0, 1),
            (vec![vec![a], vec![b, c]], "source atom indices", 2, 1),
        ] {
            assert_eq!(
                assemble_mapped_fragments(vec![carbon.clone(); 2], rows),
                Err(OperationError::InvalidAlgorithmResult {
                    operation: "fragments_with_atom_mapping_with_params",
                    field,
                    actual,
                    expected,
                })
            );
        }
    }
}
