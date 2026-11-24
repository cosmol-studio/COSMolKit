//! Registered single-molecule edits on source-indexed graph rows.
use crate::{AtomId, BondId, DerivedState, Molecule, OperationError};
use cosmolkit_core::indexed_edits as core;
use cosmolkit_macros::mol_op_body;

/// Local-to-source atom and bond indices at the time of subset extraction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubgraphMapping {
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
}
impl SubgraphMapping {
    fn from_subset(mapping: &crate::TopologyMapping) -> Result<Self, OperationError> {
        fn source_rows<I: Copy>(
            rows: &[Option<I>],
            field: &'static str,
        ) -> Result<Vec<I>, OperationError> {
            rows.iter()
                .copied()
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| OperationError::InvalidAlgorithmResult {
                    operation: "extract_subgraph",
                    field,
                    actual: rows.iter().flatten().count(),
                    expected: rows.len(),
                })
        }
        Ok(Self {
            atoms: source_rows(&mapping.atoms.new_to_old, "source atoms")?,
            bonds: source_rows(&mapping.bonds.new_to_old, "source bonds")?,
        })
    }

    pub fn source_atom_indices(&self) -> &[AtomId] {
        &self.atoms
    }
    pub fn source_bond_indices(&self) -> &[BondId] {
        &self.bonds
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MappedSubgraph {
    molecule: Molecule,
    mapping: SubgraphMapping,
}
impl MappedSubgraph {
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
    pub fn mapping(&self) -> &SubgraphMapping {
        &self.mapping
    }
}
impl From<(Molecule, SubgraphMapping)> for MappedSubgraph {
    fn from((molecule, mapping): (Molecule, SubgraphMapping)) -> Self {
        Self { molecule, mapping }
    }
}

const GRAPH_DERIVED: DerivedState = DerivedState::VALENCE
    .union(DerivedState::RINGS)
    .union(DerivedState::RING_FAMILIES)
    .union(DerivedState::AROMATICITY)
    .union(DerivedState::STEREO)
    .union(DerivedState::COORDINATES)
    .union(DerivedState::DRAWING)
    .union(DerivedState::FINGERPRINT);

#[mol_op_body(with_atoms_renumbered, parts)]
pub(crate) fn with_atoms_renumbered_impl(order: &[crate::AtomId]) -> Result<(), OperationError> {
    parts.transform_mapped_blocks(|t, c, _p| {
        core::renumber_atoms(t, c, order).map_err(OperationError::IndexedGraph)
    })?;
    parts.clear_cache(GRAPH_DERIVED)?;
    parts.apply_cip_policy()
}

#[mol_op_body(without_atom, parts)]
pub(crate) fn without_atom_impl(atom: crate::AtomId) -> Result<(), OperationError> {
    parts.transform_mapped_blocks(|t, c, p| {
        core::remove_atom(t, c, p, atom).map_err(OperationError::IndexedGraph)
    })?;
    parts.clear_cache(GRAPH_DERIVED)?;
    parts.apply_cip_policy()
}

#[mol_op_body(with_atom_replaced, parts)]
pub(crate) fn with_atom_replaced_impl(
    atom: crate::AtomId,
    replacement: &crate::AtomSpec,
    preserve_props: bool,
) -> Result<(), OperationError> {
    parts.transform_mapped_blocks(|t, c, p| {
        core::replace_atom(t, c, p, atom, replacement, preserve_props)
            .map_err(OperationError::IndexedGraph)
    })?;
    parts.clear_cache(GRAPH_DERIVED)?;
    parts.apply_cip_policy()
}

#[mol_op_body(with_bond_replaced, parts)]
pub(crate) fn with_bond_replaced_impl(
    bond: crate::BondId,
    replacement: &crate::BondSpec,
    preserve_props: bool,
    keep_sgroups: bool,
) -> Result<(), OperationError> {
    parts.transform_mapped_blocks(|t, c, p| {
        core::replace_bond(t, c, p, bond, replacement, preserve_props, keep_sgroups)
            .map_err(OperationError::IndexedGraph)
    })?;
    parts.clear_cache(GRAPH_DERIVED)?;
    parts.apply_cip_policy()
}

#[mol_op_body(with_molecule_appended, parts)]
pub(crate) fn with_molecule_appended_impl(right: &crate::Molecule) -> Result<(), OperationError> {
    parts.transform_mapped_blocks(|t, c, p| {
        core::append_molecule(t, c, p, right.topology(), right.coordinate_block_runtime())
            .map_err(OperationError::IndexedGraph)
    })?;
    parts.clear_cache(GRAPH_DERIVED)?;
    parts.apply_cip_policy()
}

#[mol_op_body(extract_subgraph_with_params, parts)]
pub(crate) fn extract_subgraph_with_params_impl(
    selection: &crate::SubgraphSelection,
    params: &crate::SubgraphExtractionParams,
) -> Result<SubgraphMapping, OperationError> {
    let mapping = parts.transform_mapped_blocks(|t, c, p| {
        cosmolkit_core::extract_indexed_subgraph(t, c, p, selection, params)
            .map_err(|e| OperationError::IndexedGraph(e.into()))
    })?;
    let report = SubgraphMapping::from_subset(mapping)?;
    parts.clear_cache(GRAPH_DERIVED)?;
    parts.apply_cip_policy()?;
    Ok(report)
}
