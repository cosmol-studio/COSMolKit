//! Query-preserving projection of source-ordered connected components.

use crate::SmartsParseError;
use cosmolkit_model::{
    AtomId, AtomMapping, BondId, BondMapping, QueryGraph, StereoGroup, SubstanceGroup,
    SubstanceGroupId, TopologyMapping, query_substance_groups, replace_query_substance_groups,
};

/// Fragment detached reaction-agent templates without passing through a live
/// molecule or converting arbitrary query identities into concrete elements.
#[doc(hidden)]
pub fn query_graph_fragments(graph: &QueryGraph) -> Result<Vec<QueryGraph>, SmartsParseError> {
    // RDKit❗❌:   int nFrags = getMolFrags(mol, *frags);
    // RDKit❗❌:   if (nFrags == 1) {
    // RDKit❗❌:     res.emplace_back(new RWMol(mol));
    // RDKit❗❌:   } else {
    // RDKit❗❌:     res.reserve(nFrags);
    // RDKit❗❌:     for (int i = 0; i < nFrags; ++i) {
    // The exact CORE connected-components implementation and existing SEARCH
    // removal-reference helpers are reused. Every query carrier and predicate
    // origin survives; the source's one-component full copy is a distinct path.
    // Complexity review: per-component source-row projection is O(K(V+E));
    // coordinate_block clones the source coordinates once, an extra O(CV)
    // input copy over source's borrowed positions, hence the cost marker.
    let components = cosmolkit_core::query_connected_components(graph)
        .map_err(|error| SmartsParseError::Parse(error.to_string()))?;
    if components.components.len() == 1 {
        return Ok(vec![graph.clone()]);
    }
    let coordinates = graph.coordinate_block(None);
    let mut result = Vec::with_capacity(components.components.len());
    for (component, source_atoms) in components.components.iter().enumerate() {
        let mut atom_map = vec![None; graph.num_atoms()];
        for (row, old) in source_atoms.iter().enumerate() {
            atom_map[old.index()] = Some(AtomId::new(row));
        }
        let mut bond_map = vec![None; graph.num_bonds()];
        let mut source_bonds = Vec::new();
        for bond in graph.bonds() {
            if atom_map[bond.begin().index()].is_some() && atom_map[bond.end().index()].is_some() {
                bond_map[bond.id().index()] = Some(BondId::new(source_bonds.len()));
                source_bonds.push(bond.id());
            }
        }
        let mapping = TopologyMapping {
            atoms: AtomMapping {
                old_to_new: atom_map,
                new_to_old: source_atoms.iter().copied().map(Some).collect(),
            },
            bonds: BondMapping {
                old_to_new: bond_map,
                new_to_old: source_bonds.iter().copied().map(Some).collect(),
            },
        };
        mapping
            .validate_for_counts(
                graph.num_atoms(),
                source_atoms.len(),
                graph.num_bonds(),
                source_bonds.len(),
            )
            .map_err(|e| SmartsParseError::Parse(e.to_string()))?;
        // RDKit❗✔️: if (comp.size() == 1 ||
        // RDKit❗✔️:     (nFrags > 3 && !fragmentHasChallengingFeatures(comp, atomsInFrag))) {
        // RDKit❗✔️:   SubsetOptions opts{.sanitize = sanitizeFrags,
        // RDKit❗✔️:                      .clearComputedProps = true,
        // RDKit❗✔️:                      .copyCoordinates = copyConformers,
        // RDKit❗✔️:                      .method = SubsetMethod::BONDS_BETWEEN_ATOMS};
        let challenging_atom = source_atoms.iter().any(|id| {
            !matches!(
                graph.atoms()[id.index()].chiral_tag(),
                cosmolkit_types::ChiralTag::Unspecified | cosmolkit_types::ChiralTag::Other
            )
        });
        let challenging_bond = source_bonds.iter().any(|id| {
            !matches!(
                graph.bonds()[id.index()].bond().stereo(),
                cosmolkit_types::BondStereo::None | cosmolkit_types::BondStereo::Any
            )
        });
        let challenging_sgroup = query_substance_groups(graph).iter().any(|group| {
            group
                .atoms()
                .iter()
                .chain(group.parent_atoms())
                .any(|id| mapping.atoms.old_to_new[id.index()].is_some())
        });
        let challenging_stereo = graph.stereo_groups().iter().any(|group| {
            group
                .atoms()
                .iter()
                .any(|id| mapping.atoms.old_to_new[id.index()].is_some())
                || group
                    .bonds()
                    .iter()
                    .any(|id| mapping.bonds.old_to_new[id.index()].is_some())
        });
        let fast_subset = source_atoms.len() == 1
            || (components.components.len() > 3
                && !(challenging_atom
                    || challenging_bond
                    || challenging_sgroup
                    || challenging_stereo));
        let mut atoms = Vec::with_capacity(source_atoms.len());
        for (row, old) in source_atoms.iter().enumerate() {
            let mut atom = graph.atoms()[old.index()].clone().with_id(AtomId::new(row));
            atom.remap_template_attachment_order(mapping.atoms.old_to_new())
                .map_err(|source| SmartsParseError::TemplateAttachmentRemap {
                    carrier: old.index(),
                    source,
                })?;
            // RDKit❗✔️: clearComputedProps(true);
            // Both a nonempty batch deletion and clearComputedProps in the
            // subset path clear every retained atom's computed properties.
            atom.clear_computed_props()?;
            atoms.push(atom);
        }
        let mut bonds = Vec::with_capacity(source_bonds.len());
        for (row, old) in source_bonds.iter().enumerate() {
            let mut bond = graph.bonds()[old.index()].clone();
            let carrier = bond.bond();
            let begin = mapping.atoms.old_to_new[carrier.begin().index()].ok_or_else(|| {
                SmartsParseError::Parse("retained fragment bond has absent begin".into())
            })?;
            let end = mapping.atoms.old_to_new[carrier.end().index()].ok_or_else(|| {
                SmartsParseError::Parse("retained fragment bond has absent end".into())
            })?;
            let mut stereo = carrier.stereo_atoms().and_then(|[a, b]| {
                Some([
                    mapping.atoms.old_to_new[a.index()]?,
                    mapping.atoms.old_to_new[b.index()]?,
                ])
            });
            if fast_subset && let Some(ref mut references) = stereo {
                // RDKit❗✔️: for (auto &atomidx : atoms) {
                // RDKit❗✔️:   auto map = atomMapping.find(atomidx);
                // RDKit❗✔️:   if (map != atomMapping.end()) {
                // RDKit❗✔️:     atomidx = map->second;
                // RDKit❗✔️:   }
                // RDKit❗✔️: }
                // Subset.cpp executes this after already mapping the pair.
                // Keep the second lookup against ORIGINAL IDs literally.
                for reference in references {
                    if let Some(mapped) = mapping
                        .atoms
                        .old_to_new
                        .get(reference.index())
                        .copied()
                        .flatten()
                    {
                        *reference = mapped;
                    }
                }
            }
            *bond.bond_mut() = carrier
                .clone()
                .remapped(BondId::new(row), begin, end, stereo);
            bond.bond_mut().clear_computed_props()?;
            bonds.push(bond);
        }
        let mut fragment_coordinates = coordinates.clone();
        fragment_coordinates
            .remap_topology(&source_atoms.iter().map(|id| id.index()).collect::<Vec<_>>());
        let groups = if fast_subset {
            subset_stereo_groups(graph.stereo_groups(), &mapping)?
        } else {
            crate::smarts_parse::remap_query_stereo_groups_after_removal(
                graph.stereo_groups(),
                &mapping,
            )?
        };
        let substance_groups = if fast_subset {
            subset_substance_groups(query_substance_groups(graph), &mapping)
        } else {
            crate::smarts_parse::remap_query_substance_groups_after_removal(
                query_substance_groups(graph),
                &mapping,
            )?
        };
        // RDKit❗✔️: auto extracted_mol = std::make_unique<RWMol>();
        // RDKit❗✔️: if (options.clearComputedProps) {
        // RDKit❗✔️:   extracted_mol->clearComputedProps();
        // The subset branch never updateProps(mol); the batch-removal branch
        // started from a complete RWMol copy.
        let properties: Vec<_> = if fast_subset {
            Vec::new()
        } else {
            graph
                .ordered_props()
                .map(|(key, value)| (key.clone(), value.clone()))
                .collect()
        };
        let mut fragment = QueryGraph::from_parts(
            atoms,
            bonds,
            properties,
            fragment_coordinates.conformers_2d,
            fragment_coordinates.conformers_3d,
            groups,
        )
        .map_err(|error| {
            SmartsParseError::Parse(format!("agent component {component}: {error}"))
        })?;
        replace_query_substance_groups(&mut fragment, substance_groups).map_err(|error| {
            SmartsParseError::Parse(format!("agent component {component}: {error}"))
        })?;
        result.push(fragment);
    }
    Ok(result)
}

fn subset_stereo_groups(
    groups: &[StereoGroup],
    mapping: &TopologyMapping,
) -> Result<Vec<StereoGroup>, SmartsParseError> {
    // RDKit❗✔️: return objects.empty() ||
    // RDKit❗✔️:        std::any_of(objects.begin(), objects.end(), [&](auto &object) {
    // RDKit❗✔️:          return selected_indices[object->getIdx()];
    // RDKit❗✔️:        });
    // RDKit❗✔️: return is_selected_component(stereo_group.getAtoms(),
    // RDKit❗✔️:                              selection_info.selectedAtoms) &&
    // RDKit❗✔️:        is_selected_component(stereo_group.getBonds(),
    // RDKit❗✔️:                              selection_info.selectedBonds);
    // Source requires selection in each NONEMPTY domain. A deletion instead
    // retains a group with either atom or bond membership remaining.
    let selected: Vec<_> = groups
        .iter()
        .filter(|group| {
            (group.atoms().is_empty()
                || group
                    .atoms()
                    .iter()
                    .any(|id| mapping.atoms.old_to_new[id.index()].is_some()))
                && (group.bonds().is_empty()
                    || group
                        .bonds()
                        .iter()
                        .any(|id| mapping.bonds.old_to_new[id.index()].is_some()))
        })
        .cloned()
        .collect();
    Ok(cosmolkit_model::merge_absolute_stereo_groups(
        crate::smarts_parse::remap_query_stereo_groups_after_removal(&selected, mapping)?,
    ))
}

fn subset_substance_groups(
    groups: &[SubstanceGroup],
    mapping: &TopologyMapping,
) -> Vec<SubstanceGroup> {
    // RDKit❗✔️: return indices.empty() ||
    // RDKit❗✔️:        std::all_of(indices.begin(), indices.end(), selection_test);
    // RDKit❗✔️: return is_selected_component(sgroup.getAtoms(), atom_test) &&
    // RDKit❗✔️:        is_selected_component(sgroup.getBonds(), bond_test) &&
    // RDKit❗✔️:        is_selected_component(sgroup.getParentAtoms(), atom_test);
    // RDKit❗✔️: SubstanceGroup extracted_sgroup(sgroup);
    // RDKit❗✔️: update_indices(extracted_sgroup, std::mem_fn(&SubstanceGroup::getAtoms),
    // RDKit❗✔️:                std::mem_fn(&SubstanceGroup::setAtoms), atomMapping);
    // RDKit❗✔️: update_indices(extracted_sgroup,
    // RDKit❗✔️:                std::mem_fn(&SubstanceGroup::getParentAtoms),
    // RDKit❗✔️:                std::mem_fn(&SubstanceGroup::setParentAtoms), atomMapping);
    // RDKit❗✔️: update_indices(extracted_sgroup, std::mem_fn(&SubstanceGroup::getBonds),
    // RDKit❗✔️:                std::mem_fn(&SubstanceGroup::setBonds), bondMapping);
    // RDKit❗✔️: addSubstanceGroup(extracted_mol, std::move(extracted_sgroup));
    // Source leaves parent-group, cstate and attachment references as copied.
    // Preserve them literally; model validation reports an invalid reference,
    // instead of silently repairing or dropping it.
    let mut result = Vec::new();
    for group in groups {
        let atoms: Option<Vec<_>> = group
            .atoms()
            .iter()
            .map(|id| mapping.atoms.old_to_new[id.index()])
            .collect();
        let bonds: Option<Vec<_>> = group
            .bonds()
            .iter()
            .map(|id| mapping.bonds.old_to_new[id.index()])
            .collect();
        let parents: Option<Vec<_>> = group
            .parent_atoms()
            .iter()
            .map(|id| mapping.atoms.old_to_new[id.index()])
            .collect();
        let (Some(atoms), Some(bonds), Some(parents)) = (atoms, bonds, parents) else {
            continue;
        };
        let mut copied = group
            .clone()
            .with_atoms(atoms)
            .with_bonds(bonds)
            .with_parent_atoms(parents);
        copied.set_id(SubstanceGroupId::new(result.len()));
        result.push(copied);
    }
    result
}
