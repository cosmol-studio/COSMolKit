//! Private generated reaction transport; the detached owner owns chemistry.

use std::borrow::Cow;

use crate::{
    DerivedState, Molecule, OperationError, Reaction, ReactionApplyParams, ReactionRunParams,
    ReactionSingleRunParams, TopologyEditKind,
};
use cosmolkit_macros::{mol_multi_op_body, mol_op_body};

#[mol_multi_op_body(reaction_products, parts)]
pub(crate) fn reaction_products_impl(
    reaction: &mut Reaction,
    reactant_template: usize,
    params: &ReactionSingleRunParams,
) -> Result<Vec<usize>, OperationError> {
    let input = parts.reconstruction_source()?;
    let sets = cosmolkit_reaction::run_reactant(reaction, input, reactant_template, params)
        .map_err(OperationError::ReactionRun)?;
    let lengths = sets.iter().map(Vec::len).collect();
    parts.emit_reconstructed(sets.into_iter().flatten().collect())?;
    Ok(lengths)
}

#[mol_multi_op_body(reaction_products_from_inputs, parts)]
pub(crate) fn reaction_products_from_inputs_impl(
    reaction: &mut Reaction,
    reactants: &[&Molecule],
    params: &ReactionRunParams,
) -> Result<Vec<usize>, OperationError> {
    let inputs = parts.reconstruction_inputs(reactants)?;
    let sets = cosmolkit_reaction::run_reactants(reaction, &inputs, params)
        .map_err(OperationError::ReactionRun)?;
    let lengths = sets.iter().map(Vec::len).collect();
    parts.emit_reconstructed(sets.into_iter().flatten().collect())?;
    Ok(lengths)
}

#[mol_op_body(apply_reaction, parts)]
pub(crate) fn apply_reaction_impl(
    reaction: &mut Reaction,
    params: &ReactionApplyParams,
) -> Result<bool, OperationError> {
    let (changes, rows) = {
        let input = parts.reaction_input()?;
        let rows = (input.topology.atoms.len(), input.topology.bonds.len());
        (
            cosmolkit_reaction::apply_reaction(reaction, input, params)
                .map_err(OperationError::ReactionApply)?,
            rows,
        )
    };
    let has_candidate = changes.change.is_some();
    let mapping = if let Some((topology, mapping)) = changes.change {
        // Stage the already owned detached graph; retain the original borrowed
        // properties without an extra graph/property clone or commit authority.
        parts.stage_topology_properties_cow(|_, properties, _| {
            Ok(((), Some((Cow::Owned(topology), properties))))
        })?;
        mapping
    } else {
        cosmolkit_model::TopologyMapping::identity(rows.0, rows.1)
    };
    parts.record_topology_edit(TopologyEditKind::Compacting)?;
    parts.record_topology_mapping(mapping)?;
    parts.apply_runtime_remap()?;
    if changes.clears_computed_properties {
        let mut properties = parts.checkout_properties()?;
        properties
            .clear_computed_props()
            .map_err(OperationError::InvalidProperty)?;
        parts.install_properties(properties)?;
    }
    if has_candidate {
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
    }
    parts.apply_cip_policy()?;
    Ok(changes.changed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, Conformer3D,
        CoordinateBlock, CoordinateDimension, Element, MoleculeProperties, TopologyBlock,
    };
    use std::sync::Arc;

    fn molecule(elements: &[Element], edges: &[(usize, usize)]) -> Molecule {
        let atoms = elements
            .iter()
            .enumerate()
            .map(|(row, &element)| Atom::from_spec(AtomId::new(row), AtomSpec::new(element)))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(row, &(a, b))| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect();
        Molecule::from_parts(
            TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .unwrap()
    }

    fn assert_shared(left: &Molecule, right: &Molecule) {
        assert_eq!(left, right);
        assert!(std::ptr::eq(left.topology(), right.topology()));
        assert!(std::ptr::eq(
            left.coordinate_block_runtime(),
            right.coordinate_block_runtime()
        ));
        assert!(std::ptr::eq(left.properties(), right.properties()));
        assert!(Arc::ptr_eq(
            &left.derived_cache_arc_runtime(),
            &right.derived_cache_arc_runtime()
        ));
    }

    #[test]
    fn reaction_private_pipeline_single_run_initializes_actual_reaction_and_keeps_source() {
        let source = molecule(&[Element::C], &[]);
        let observer = source.clone();
        let mut reaction = cosmolkit_reaction::parse_smirks("[C:1]>>[N:1]").unwrap();
        assert!(!reaction.is_initialized());
        let sets = source.reaction_products(&mut reaction, 0).unwrap();
        assert!(reaction.is_initialized());
        assert_eq!(sets.len(), 1);
        assert_eq!(sets[0].len(), 1);
        assert_eq!(sets[0][0].topology().atoms[0].element(), Element::N);
        assert_shared(&source, &observer);
    }

    #[test]
    fn reaction_private_pipeline_source_product_growth_accepts_new_rows() {
        let source = molecule(&[Element::C], &[]);
        let observer = source.clone();
        let mut reaction = cosmolkit_reaction::parse_smirks("[C:1]>>[C:1]C").unwrap();
        let sets = source.reaction_products(&mut reaction, 0).unwrap();
        assert_eq!(sets.len(), 1);
        assert_eq!(sets[0].len(), 1);
        assert_eq!(sets[0][0].topology().atoms.len(), 2);
        assert_eq!(sets[0][0].topology().bonds.len(), 1);
        assert_shared(&source, &observer);
    }

    #[test]
    fn reaction_private_pipeline_actual_duplicate_inputs_keep_both_product_slots() {
        let source = molecule(&[Element::C], &[]);
        let observer = source.clone();
        let mut reaction = cosmolkit_reaction::parse_smirks("[C:1].[C:2]>>[C:1].[C:2]").unwrap();
        let sets = source
            .reaction_products_from_inputs(
                &mut reaction,
                &[&source, &source],
                &ReactionRunParams::default(),
            )
            .unwrap();
        assert_eq!(sets.len(), 1);
        assert_eq!(sets[0].len(), 2);
        assert_eq!(sets[0][0].topology().atoms.len(), 1);
        assert_eq!(sets[0][1].topology().atoms.len(), 1);
        assert_shared(&source, &observer);
    }

    #[test]
    fn reaction_private_pipeline_no_match_keeps_live_cache_and_every_shared_block() {
        let source = molecule(&[Element::O], &[]);
        let mut cache = crate::molecule::DerivedCacheBlock::default();
        cache.install_valence_assignment(cosmolkit_core::ValenceAssignment {
            explicit_valence: vec![0],
            implicit_hydrogens: vec![2],
        });
        cache.mark_valid(DerivedState::VALENCE);
        let mut source = Molecule::from_runtime_parts(
            source.topology_arc_runtime(),
            source.coordinates_arc_runtime(),
            source.properties_arc_runtime(),
            Arc::new(cache),
        )
        .unwrap();
        let observer = source.clone();
        let mut reaction = cosmolkit_reaction::parse_smirks("[C:1]>>[N:1]").unwrap();
        let result = source.apply_reaction(&mut reaction).unwrap();
        assert!(!result.changed);
        assert_shared(&result.molecule, &observer);
        let changed: bool = source.apply_reaction_(&mut reaction).unwrap();
        assert!(!changed);
        assert_shared(&source, &observer);
    }

    #[test]
    fn reaction_private_pipeline_apply_compaction_remaps_all_actual_coordinate_frames() {
        let source = molecule(&[Element::C, Element::N], &[(0, 1)]);
        let coordinates = CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(9, vec![[1.0, 2.0], [3.0, 4.0]]),
                Conformer2D::new(2, vec![[5.0, 6.0], [7.0, 8.0]]),
            ],
            conformers_3d: vec![Conformer3D::new(
                12,
                vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                true,
            )],
            source_conformer_order: Some(vec![
                CoordinateDimension::TwoD,
                CoordinateDimension::ThreeD,
                CoordinateDimension::TwoD,
            ]),
            ..Default::default()
        };
        let mut source = Molecule::from_parts(
            source.topology().clone(),
            coordinates,
            MoleculeProperties::default(),
        )
        .unwrap();
        let observer = source.clone();
        let mut reaction = cosmolkit_reaction::parse_smirks("[C:1][N:2]>>[C:1]").unwrap();
        let result = source.apply_reaction(&mut reaction).unwrap();
        assert!(result.changed);
        assert_eq!(result.molecule.topology().atoms.len(), 1);
        assert_eq!(
            result.molecule.coordinate_block_runtime().conformers_2d[0].coordinates(),
            &[[1.0, 2.0]]
        );
        assert_eq!(
            result.molecule.coordinate_block_runtime().conformers_2d[1].coordinates(),
            &[[5.0, 6.0]]
        );
        assert_eq!(
            result.molecule.coordinate_block_runtime().conformers_3d[0].coordinates(),
            &[[1.0, 2.0, 3.0]]
        );
        assert_shared(&source, &observer);
        assert!(source.apply_reaction_(&mut reaction).unwrap());
        assert_eq!(source, result.molecule);
    }

    #[test]
    fn reaction_private_pipeline_failure_preserves_molecule_and_source_initialization_prefix() {
        let source = molecule(&[Element::C], &[]);
        let observer = source.clone();
        let mut reaction = cosmolkit_reaction::parse_smirks("[C:1]>>[N:1]").unwrap();
        assert!(matches!(
            source.reaction_products(&mut reaction, 99),
            Err(OperationError::ReactionRun(
                cosmolkit_reaction::ReactionRunError::ReactantTemplateIndex {
                    index: 99,
                    count: 1
                }
            ))
        ));
        assert!(reaction.is_initialized());
        assert_shared(&source, &observer);
    }

    #[test]
    fn reaction_private_pipeline_computed_properties_follow_actual_source_removals() {
        // ReactionRunner.cpp::run_Reactant(RWMol&) uses commitBatchEdit().
        // RWMol.cpp::commitBatchEdit returns without clearing when neither
        // deletion mask contains a removal; otherwise clearComputedProps(true)
        // clears computed molecule/atom/bond properties, retaining ordinary ones.
        for (smirks, clears, atom_count, bond_count) in [
            ("[C:1][C:2]>>[C:1]", true, 1, 0),
            ("[C:1][C:2]>>([C:1].[C:2])", true, 2, 0),
            ("[C:1][C:2]>>[C:1]=[C:2]", false, 2, 1),
            ("[C:1][C:2]>>[N:1][C:2]", false, 2, 1),
            ("[C:1][C:2]>>[C:1][C:2]", false, 2, 1),
            ("[O:1]>>[N:1]", false, 2, 1),
        ] {
            let source = molecule(&[Element::C, Element::C], &[(0, 1)]);
            let mut topology = source.topology().clone();
            for atom in &mut topology.atoms {
                atom.set_computed_prop("_CIPCode", "R").unwrap();
                atom.set_computed_prop("_CIPRank", 7_u32).unwrap();
                atom.set_prop("ordinary", "atom").unwrap();
                atom.set_prop("_CIPNeighborOrder", "retained").unwrap();
            }
            topology.bonds[0]
                .set_computed_prop("_CIPCode", "E")
                .unwrap();
            topology.bonds[0].set_prop("ordinary", "bond").unwrap();
            let mut properties = MoleculeProperties::default();
            properties.set_computed_prop("_CIPComputed", true).unwrap();
            properties
                .set_computed_prop("other_computed", 9_u32)
                .unwrap();
            properties.set_prop("ordinary", "molecule").unwrap();
            let mut source =
                Molecule::from_parts(topology, CoordinateBlock::default(), properties).unwrap();
            let observer = source.clone();
            let mut reaction = cosmolkit_reaction::parse_smirks(smirks).unwrap();
            let result = source.apply_reaction(&mut reaction).unwrap();
            let product = &result.molecule;
            assert_eq!(product.topology().atoms.len(), atom_count, "{smirks}");
            assert_eq!(product.topology().bonds.len(), bond_count, "{smirks}");
            for atom in &product.topology().atoms {
                assert_eq!(atom.prop("_CIPCode").is_none(), clears, "{smirks}");
                assert_eq!(atom.prop("_CIPRank").is_none(), clears, "{smirks}");
                assert_eq!(
                    atom.prop("ordinary"),
                    observer.topology().atoms[0].prop("ordinary")
                );
                assert_eq!(
                    atom.prop("_CIPNeighborOrder"),
                    observer.topology().atoms[0].prop("_CIPNeighborOrder")
                );
            }
            for bond in &product.topology().bonds {
                assert_eq!(bond.prop("_CIPCode").is_none(), clears, "{smirks}");
                assert_eq!(
                    bond.prop("ordinary"),
                    observer.topology().bonds[0].prop("ordinary")
                );
            }
            assert_eq!(
                product.properties().prop("_CIPComputed").is_none(),
                clears,
                "{smirks}"
            );
            assert_eq!(
                product.properties().prop("other_computed").is_none(),
                clears,
                "{smirks}"
            );
            assert_eq!(
                product.properties().prop("ordinary"),
                observer.properties().prop("ordinary")
            );
            assert_shared(&source, &observer);
            assert_eq!(
                source.apply_reaction_(&mut reaction).unwrap(),
                result.changed
            );
            assert_eq!(source, result.molecule, "{smirks}");
        }
    }

    #[test]
    fn reaction_private_pipeline_product_growth_retains_source_copied_atom_properties() {
        let source = molecule(&[Element::C, Element::O], &[(0, 1)]);
        let mut topology = source.topology().clone();
        topology.atoms[1]
            .set_computed_prop("_CIPCode", "R")
            .unwrap();
        topology.atoms[1].set_prop("ordinary", "spectator").unwrap();
        let source = Molecule::from_parts(
            topology,
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .unwrap();
        let observer = source.clone();
        let mut reaction = cosmolkit_reaction::parse_smirks("[C:1]>>[C:1]C").unwrap();
        let sets = source.reaction_products(&mut reaction, 0).unwrap();
        let product = &sets[0][0];
        assert_eq!(product.topology().atoms.len(), 3);
        // ReactionRunner.cpp::addMissingProductAtom copies the original Atom,
        // including its computed property membership, into the growing product.
        let oxygen = product
            .topology()
            .atoms
            .iter()
            .find(|atom| atom.element() == Element::O)
            .unwrap();
        assert_eq!(
            oxygen.prop("_CIPCode"),
            source.topology().atoms[1].prop("_CIPCode")
        );
        assert_eq!(oxygen.is_prop_computed("_CIPCode").unwrap(), true);
        assert_eq!(
            oxygen.prop("ordinary"),
            source.topology().atoms[1].prop("ordinary")
        );
        assert_shared(&source, &observer);
    }
}
