// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use cosmolkit_core::{CanonicalRankError, CanonicalRankParams, rank_mol_atoms_with_params};
use cosmolkit_model::{AdjacencyList, AtomId, BondId, StereoGroup, TopologyBlock};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CanonicalRankPolicy {
    pub break_ties: bool,
    pub include_chirality: bool,
    pub include_isotopes: bool,
    pub include_atom_maps: bool,
    pub include_chiral_presence: bool,
    pub include_stereo_groups: bool,
    pub use_non_stereo_ranks: bool,
    pub include_ring_stereo: bool,
}

impl Default for CanonicalRankPolicy {
    fn default() -> Self {
        Self {
            break_ties: true,
            include_chirality: true,
            include_isotopes: true,
            include_atom_maps: true,
            include_chiral_presence: false,
            include_stereo_groups: true,
            use_non_stereo_ranks: false,
            include_ring_stereo: true,
        }
    }
}

impl CanonicalRankPolicy {
    fn core_params(self) -> CanonicalRankParams {
        let mut params = CanonicalRankParams::default();
        params.break_ties = self.break_ties;
        params.include_chirality = self.include_chirality;
        params.include_isotopes = self.include_isotopes;
        params.include_atom_maps = self.include_atom_maps;
        params.include_chiral_presence = self.include_chiral_presence;
        params.include_stereo_groups = self.include_stereo_groups;
        params.use_non_stereo_ranks = self.use_non_stereo_ranks;
        params.include_ring_stereo = self.include_ring_stereo;
        params
    }
}

pub(crate) fn rank_component_atoms(
    topology: &TopologyBlock,
    component_atoms: &[usize],
) -> Result<Vec<usize>, CanonicalRankError> {
    rank_component_atoms_with_policy(topology, component_atoms, CanonicalRankPolicy::default())
}

pub(crate) fn rank_component_atoms_with_policy(
    topology: &TopologyBlock,
    component_atoms: &[usize],
    policy: CanonicalRankPolicy,
) -> Result<Vec<usize>, CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment ranking boundary
    // RDKit✔️✔️:   auto mols =
    // RDKit✔️✔️:       MolOps::getMolFrags(mol, false, nullptr, &fragsMolAtomMapping, false);
    // RDKit✔️✔️:   for (unsigned fragIdx = 0; fragIdx < mols.size(); fragIdx++) {
    // RDKit✔️✔️:     ROMol *tmol = mols[fragIdx].get();
    // RDKit✔️✔️:     if (params.canonical) {
    // RDKit✔️✔️:       Canon::rankMolAtoms(*tmol, ranks, breakTies, includeChirality,
    // RDKit✔️✔️:                           includeIsotopes, includeAtomMaps,
    // RDKit✔️✔️:                           includeChiralPresence, includeStereoGroups,
    // RDKit✔️✔️:                           useNonStereoRanks);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment ranking boundary
    let fragment = build_ranking_fragment(topology, component_atoms)?;
    let local_ranks = rank_mol_atoms_with_params(&fragment, &policy.core_params())?;
    let mut ranks = vec![usize::MAX; topology.atoms.len()];
    for (new_index, &old_index) in component_atoms.iter().enumerate() {
        ranks[old_index] = local_ranks[new_index];
    }
    Ok(ranks)
}

/// Builds the detached renumbered fragment for one connected component.
/// Stays private: the fragment is a ranking intermediate, and the shared
/// attachment remapper is the only template-attachment maintenance path.
fn build_ranking_fragment(
    topology: &TopologyBlock,
    component_atoms: &[usize],
) -> Result<TopologyBlock, CanonicalRankError> {
    // Pass 1: validate every selected index and build the complete
    // old-to-new mapping before any atom is remapped. Forward references
    // (for example old atom 1 targeting the later retained old atom 2) must
    // resolve against the finished map, not a partially-populated one.
    let mut old_to_new = vec![None; topology.atoms.len()];
    for (new_index, &old_index) in component_atoms.iter().enumerate() {
        if old_index >= topology.atoms.len() {
            return Err(CanonicalRankError::ComponentAtomOutOfRange {
                atom_index: old_index,
                atom_count: topology.atoms.len(),
            });
        }
        if old_to_new[old_index].is_some() {
            return Err(CanonicalRankError::DuplicateComponentAtom {
                atom_index: old_index,
            });
        }
        old_to_new[old_index] = Some(AtomId::new(new_index));
    }
    // Pass 2: clone, renumber, and remap each selected atom through the one
    // shared primitive. A carrier whose target is genuinely excluded from
    // the component still fails under the established lost-target policy
    // instead of keeping a stale row.
    let mut atoms = Vec::with_capacity(component_atoms.len());
    for &old_index in component_atoms {
        let mut atom = topology.atoms[old_index]
            .clone()
            .with_id(old_to_new[old_index].expect("selected index mapped in the first pass"));
        atom.remap_template_attachment_order(&old_to_new)
            .map_err(|source| CanonicalRankError::TemplateAttachmentRemap {
                carrier: topology.atoms[old_index].id(),
                source,
            })?;
        atoms.push(atom);
    }
    let mut bonds = Vec::new();
    let mut bond_old_to_new = vec![None; topology.bonds.len()];
    for (old_bond_index, bond) in topology.bonds.iter().enumerate() {
        let (Some(begin), Some(end)) = (
            old_to_new[bond.begin().index()],
            old_to_new[bond.end().index()],
        ) else {
            continue;
        };
        let stereo_atoms = bond.stereo_atoms().and_then(|stereo_atoms| {
            Some([
                old_to_new[stereo_atoms[0].index()]?,
                old_to_new[stereo_atoms[1].index()]?,
            ])
        });
        let new_bond = BondId::new(bonds.len());
        bond_old_to_new[old_bond_index] = Some(new_bond);
        bonds.push(bond.clone().remapped(new_bond, begin, end, stereo_atoms));
    }
    let stereo_groups =
        remap_component_stereo_groups(&topology.stereo_groups, &old_to_new, &bond_old_to_new);
    Ok(TopologyBlock {
        adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
        atoms,
        bonds,
        substance_groups: Vec::new(),
        stereo_groups,
    })
}

fn remap_component_stereo_groups(
    groups: &[StereoGroup],
    atom_old_to_new: &[Option<AtomId>],
    bond_old_to_new: &[Option<BondId>],
) -> Vec<StereoGroup> {
    // BEGIN RDKIT CPP FUNCTION Subset.cpp::copySelectedStereoGroups
    // RDKit✔️✔️:   std::vector<StereoGroup> extracted_stereo_groups;
    // RDKit✔️✔️:   for (const auto &stereo_group : reference_mol.getStereoGroups()) {
    // RDKit✔️✔️:     if (!is_selected_stereo_group(stereo_group)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     std::vector<Atom *> atoms;
    // RDKit✔️✔️:     for (const auto &atom : stereo_group.getAtoms()) {
    // RDKit✔️✔️:       auto mapping = atomMapping.find(atom->getIdx());
    // RDKit✔️✔️:       if (mapping != atomMapping.end()) {
    // RDKit✔️✔️:         atoms.push_back(extracted_atoms[mapping->second]);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     std::vector<Bond *> bonds;
    // RDKit✔️✔️:     for (const auto &bond : stereo_group.getBonds()) {
    // RDKit✔️✔️:       auto mapping = bondMapping.find(bond->getIdx());
    // RDKit✔️✔️:       if (mapping != bondMapping.end()) {
    // RDKit✔️✔️:         bonds.push_back(extracted_bonds[mapping->second]);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     extracted_stereo_groups.push_back({stereo_group.getGroupType(),
    // RDKit✔️✔️:                                        std::move(atoms), std::move(bonds),
    // RDKit✔️✔️:                                        stereo_group.getReadId()});
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION Subset.cpp::copySelectedStereoGroups
    groups
        .iter()
        .filter_map(|group| {
            let atoms = group
                .atoms()
                .iter()
                .filter_map(|atom| atom_old_to_new.get(atom.index()).copied().flatten())
                .collect::<Vec<_>>();
            let bonds = group
                .bonds()
                .iter()
                .filter_map(|bond| bond_old_to_new.get(bond.index()).copied().flatten())
                .collect::<Vec<_>>();
            if (!group.atoms().is_empty() && atoms.is_empty())
                || (!group.bonds().is_empty() && bonds.is_empty())
            {
                return None;
            }
            let remapped = StereoGroup::new(group.kind(), atoms, bonds);
            Some(match group.id() {
                Some(id) => remapped.with_id(id),
                None => remapped,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_smiles;
    use cosmolkit_model::{
        Atom, AtomSpec, Bond, BondSpec, StereoGroupKind, TemplateAttachment,
        TemplateAttachmentOrder, TemplateAttachmentOrderError,
    };
    use cosmolkit_types::{BondOrder, Element};

    #[test]
    fn detached_rank_mol_port_matches_rdkit_canonical_rank_atoms() {
        for (smiles, expected) in [
            ("OCC", vec![1, 2, 0]),
            ("OC(C)C", vec![2, 3, 0, 1]),
            ("c1ccncc1", vec![0, 1, 3, 5, 4, 2]),
            ("C12(CCCCC1)CCCCC2", vec![10, 6, 2, 0, 3, 7, 8, 4, 1, 5, 9]),
            ("C1C2C3C1C2C3", vec![0, 2, 4, 3, 5, 1]),
            ("[13CH3]CO", vec![0, 2, 1]),
            ("[CH3:7]CO", vec![2, 1, 0]),
        ] {
            let record = parse_smiles(smiles, &Default::default()).unwrap();
            let atoms = (0..record.topology.atoms.len()).collect::<Vec<_>>();
            assert_eq!(
                rank_component_atoms(&record.topology, &atoms).unwrap(),
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn component_selection_rejects_invalid_rows_and_marks_nonmembers() {
        let record = parse_smiles("CC.O", &Default::default()).unwrap();
        assert_eq!(
            rank_component_atoms(&record.topology, &[2]).unwrap(),
            vec![usize::MAX, usize::MAX, 0]
        );
        assert_eq!(
            rank_component_atoms(&record.topology, &[]).unwrap(),
            vec![usize::MAX; 3]
        );
        assert_eq!(
            rank_component_atoms(&record.topology, &[0, 0]),
            Err(CanonicalRankError::DuplicateComponentAtom { atom_index: 0 })
        );
        assert_eq!(
            rank_component_atoms(&record.topology, &[3]),
            Err(CanonicalRankError::ComponentAtomOutOfRange {
                atom_index: 3,
                atom_count: 3,
            })
        );
    }

    #[test]
    fn component_stereo_groups_keep_selected_members_and_source_order() {
        let groups = vec![
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(0), AtomId::new(2)],
                vec![BondId::new(0), BondId::new(1)],
            )
            .with_id(7),
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(1)],
                vec![BondId::new(0)],
            )
            .with_id(9),
        ];
        let atom_map = vec![Some(AtomId::new(1)), None, Some(AtomId::new(0))];
        let bond_map = vec![None, Some(BondId::new(0))];

        let remapped = remap_component_stereo_groups(&groups, &atom_map, &bond_map);
        assert_eq!(remapped.len(), 1);
        assert_eq!(remapped[0].kind(), StereoGroupKind::Or);
        assert_eq!(remapped[0].id(), Some(7));
        assert_eq!(remapped[0].atoms(), &[AtomId::new(1), AtomId::new(0)]);
        assert_eq!(remapped[0].bonds(), &[BondId::new(0)]);
    }

    fn two_component_topology() -> TopologyBlock {
        // Component A: old atom 0. Component B: old atoms 1, 2. The carrier at
        // old atom 1 targets old atom 2 (inside its own component).
        let order =
            TemplateAttachmentOrder::new(vec![TemplateAttachment::new(AtomId::new(2), "port")])
                .unwrap();
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(
                AtomId::new(1),
                AtomSpec::new(Element::N).with_template_attachment_order(order),
            ),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
        )];
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    #[test]
    fn ranking_fragment_remaps_surviving_template_attachment_targets() {
        let topology = two_component_topology();
        let fragment =
            build_ranking_fragment(&topology, &[1, 2]).expect("component fragment builds");
        // The carrier at old atom 1 becomes new atom 0, and its target old
        // atom 2 becomes new atom 1: the fragment state is internally
        // consistent, not a stale copy of the source row.
        assert_eq!(fragment.atoms.len(), 2);
        let order = fragment.atoms[0]
            .template_attachment_order()
            .expect("carrier keeps its typed state in the fragment");
        assert_eq!(order.entries().len(), 1);
        assert_eq!(order.entries()[0].target(), AtomId::new(1));
        assert_eq!(order.entries()[0].label(), "port");
        fragment
            .validate()
            .expect("fragment stays a valid topology");
    }

    #[test]
    fn ranking_fragment_fails_when_a_carrier_loses_its_target() {
        // The same carrier in a fragment that excludes its target (old atom 2
        // lives in the other component) fails under the established
        // lost-target policy instead of surviving with a stale row.
        let topology = two_component_topology();
        assert!(matches!(
            build_ranking_fragment(&topology, &[0, 1]),
            Err(CanonicalRankError::TemplateAttachmentRemap {
                carrier,
                source: TemplateAttachmentOrderError::TargetRemoved { target, .. }
            }) if carrier == AtomId::new(1) && target == AtomId::new(2)
        ));
    }
}
