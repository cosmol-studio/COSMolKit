//! Detached source-backed modern RDKit CIP label assignment.

use cosmolkit_core::{atropisomer_carriers, find_double_bond_stereo_atoms};
use cosmolkit_model::{AtomId, BondId, MoleculeProperties, TopologyBlock};
use cosmolkit_types::{BondOrder, BondStereo, ChiralTag};

use crate::cip_graph::{
    CipDigraph, CipEdgeId, CipLabelerContext, CipLabelerError, CipNode, CipNodeId, CipRules,
    CipSequenceRule, Descriptor, cip_all_rules, cip_constitutional_rules, descriptor_to_string,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CipLabelOptions {
    atoms: Option<Vec<AtomId>>,
    bonds: Option<Vec<BondId>>,
    max_recursive_iterations: u32,
}

impl Default for CipLabelOptions {
    fn default() -> Self {
        Self {
            atoms: None,
            bonds: None,
            max_recursive_iterations: 0,
        }
    }
}

impl CipLabelOptions {
    #[must_use]
    pub fn with_atoms(mut self, atoms: impl IntoIterator<Item = AtomId>) -> Self {
        self.atoms = Some(atoms.into_iter().collect());
        self
    }

    #[must_use]
    pub fn with_bonds(mut self, bonds: impl IntoIterator<Item = BondId>) -> Self {
        self.bonds = Some(bonds.into_iter().collect());
        self
    }

    #[must_use]
    pub const fn with_max_recursive_iterations(mut self, limit: u32) -> Self {
        self.max_recursive_iterations = limit;
        self
    }

    #[must_use]
    pub fn atoms(&self) -> Option<&[AtomId]> {
        self.atoms.as_deref()
    }

    #[must_use]
    pub fn bonds(&self) -> Option<&[BondId]> {
        self.bonds.as_deref()
    }

    #[must_use]
    pub const fn max_recursive_iterations(&self) -> u32 {
        self.max_recursive_iterations
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CipLabelAssignment {
    topology: TopologyBlock,
    properties: MoleculeProperties,
}

impl CipLabelAssignment {
    #[must_use]
    pub const fn topology(&self) -> &TopologyBlock {
        &self.topology
    }

    #[must_use]
    pub const fn properties(&self) -> &MoleculeProperties {
        &self.properties
    }

    #[must_use]
    pub fn into_parts(self) -> (TopologyBlock, MoleculeProperties) {
        (self.topology, self.properties)
    }
}

// BEGIN RDKIT CPP FUNCTION findConfigs (CIPLabeler.cpp)
// RDKit✔️✔️: std::vector<std::unique_ptr<Configuration>> findConfigs(
// RDKit✔️✔️:     CIPMol &mol, const boost::dynamic_bitset<> &atoms,
// RDKit✔️✔️:     const boost::dynamic_bitset<> &bonds) {
// RDKit✔️✔️:   std::vector<std::unique_ptr<Configuration>> configs;
// RDKit✔️✔️:
// RDKit✔️✔️:   for (auto index = atoms.find_first(); index != boost::dynamic_bitset<>::npos;
// RDKit✔️✔️:        index = atoms.find_next(index)) {
// RDKit✔️✔️:     auto atom = mol.getAtom(index);
// RDKit✔️✔️:     auto chiraltag = atom->getChiralTag();
// RDKit✔️✔️:     if (chiraltag == Atom::CHI_TETRAHEDRAL_CW ||
// RDKit✔️✔️:         chiraltag == Atom::CHI_TETRAHEDRAL_CCW) {
// RDKit✔️✔️:       std::unique_ptr<Tetrahedral> cfg{new Tetrahedral(mol, atom)};
// RDKit✔️✔️:       configs.push_back(std::move(cfg));
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   for (auto index = bonds.find_first(); index != boost::dynamic_bitset<>::npos;
// RDKit✔️✔️:        index = bonds.find_next(index)) {
// RDKit✔️✔️:     auto bond = mol.getBond(index);
// RDKit✔️✔️:
// RDKit✔️✔️:     auto bond_cfg = bond->getStereo();
// RDKit✔️✔️:     switch (bond_cfg) {
// RDKit✔️✔️:       case Bond::STEREOE:
// RDKit✔️✔️:         bond_cfg = Bond::STEREOTRANS;
// RDKit✔️✔️:         break;
// RDKit✔️✔️:       case Bond::STEREOZ:
// RDKit✔️✔️:         bond_cfg = Bond::STEREOCIS;
// RDKit✔️✔️:         break;
// RDKit✔️✔️:       default:
// RDKit✔️✔️:         break;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     switch (bond_cfg) {
// RDKit✔️✔️:       case Bond::STEREOTRANS:
// RDKit✔️✔️:       case Bond::STEREOCIS: {
// RDKit✔️✔️:         std::unique_ptr<Sp2Bond> cfg(new Sp2Bond(
// RDKit✔️✔️:             mol, bond, bond->getBeginAtom(), bond->getEndAtom(), bond_cfg));
// RDKit✔️✔️:         configs.push_back(std::move(cfg));
// RDKit✔️✔️:       } break;
// RDKit✔️✔️:
// RDKit✔️✔️:       case Bond::STEREOATROPCCW:
// RDKit✔️✔️:       case Bond::STEREOATROPCW: {
// RDKit✔️✔️:         std::unique_ptr<AtropisomerBond> cfgAtrop(new AtropisomerBond(
// RDKit✔️✔️:             mol, bond, bond->getBeginAtom(), bond->getEndAtom(), bond_cfg));
// RDKit✔️✔️:         configs.push_back(std::move(cfgAtrop));
// RDKit✔️✔️:       } break;
// RDKit✔️✔️:
// RDKit✔️✔️:       default:
// RDKit✔️✔️:         break;
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   return configs;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION findConfigs
fn cip_find_configs<'a>(
    molecule: &'a TopologyBlock,
    atom_mask: &[bool],
    bond_mask: &[bool],
) -> Result<Vec<CipConfig<'a>>, CipLabelerError> {
    let mut configs = Vec::new();

    for (index, selected) in atom_mask.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        let atom = molecule
            .atoms
            .get(index)
            .ok_or(CipLabelerError::AtomIndexOutOfRange {
                index,
                atom_count: molecule.atoms.len(),
            })?;
        match atom.chiral_tag() {
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                configs.push(CipConfig::Tetrahedral(CipTetrahedral::new(
                    molecule, index,
                )?));
            }
            ChiralTag::Tetrahedral
            | ChiralTag::Allene
            | ChiralTag::SquarePlanar
            | ChiralTag::TrigonalBipyramidal
            | ChiralTag::Octahedral => {
                return Err(CipLabelerError::UnsupportedConfiguration {
                    atom: index,
                    tag: atom.chiral_tag(),
                });
            }
            ChiralTag::Unspecified | ChiralTag::Other => {}
        }
    }

    for (index, selected) in bond_mask.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        let bond = molecule
            .bonds
            .get(index)
            .ok_or(CipLabelerError::BondIndexOutOfRange {
                index,
                bond_count: molecule.bonds.len(),
            })?;
        let bond_cfg = match bond.stereo() {
            BondStereo::E => BondStereo::Trans,
            BondStereo::Z => BondStereo::Cis,
            other => other,
        };
        match bond_cfg {
            BondStereo::Trans | BondStereo::Cis => {
                configs.push(CipConfig::Sp2Bond(CipSp2Bond::new(
                    molecule,
                    index,
                    bond.begin().index(),
                    bond.end().index(),
                    bond_cfg,
                )?));
            }
            BondStereo::AtropCcw | BondStereo::AtropCw => {
                configs.push(CipConfig::AtropisomerBond(CipAtropisomerBond::new(
                    molecule,
                    index,
                    bond.begin().index(),
                    bond.end().index(),
                    bond_cfg,
                )?));
            }
            _ => {}
        }
    }

    Ok(configs)
}

fn cip_label_with_center_digraph(
    configs: &mut [CipConfig<'_>],
    center_idx: usize,
    target_idx: usize,
    node: CipNodeId,
    rules: &CipRules,
    context: &mut CipLabelerContext,
) -> Result<Descriptor, CipLabelerError> {
    debug_assert_ne!(center_idx, target_idx);
    if center_idx < target_idx {
        let (left, right) = configs.split_at_mut(target_idx);
        let center = &mut left[center_idx];
        let target = &mut right[0];
        let digraph = center.get_digraph_mut();
        target.label_with_external_digraph(node, digraph, rules, context)
    } else {
        let (left, right) = configs.split_at_mut(center_idx);
        let target = &mut left[target_idx];
        let center = &mut right[0];
        let digraph = center.get_digraph_mut();
        target.label_with_external_digraph(node, digraph, rules, context)
    }
}

fn cip_set_center_node_aux(
    configs: &mut [CipConfig<'_>],
    center_idx: usize,
    node: CipNodeId,
    desc: Descriptor,
) {
    let center = &mut configs[center_idx];
    center.get_digraph_mut().set_node_aux(node, desc);
}

// BEGIN RDKIT CPP FUNCTION labelAux (CIPLabeler.cpp)
// RDKit✔️✔️: bool labelAux(std::vector<std::unique_ptr<Configuration>> &configs,
// RDKit✔️✔️:               const Rules &rules,
// RDKit✔️✔️:               const std::unique_ptr<Configuration> &center) {
// RDKit✔️✔️:   using Node_Cfg_Pair = std::pair<Node *, Configuration *>;
// RDKit✔️✔️:   std::vector<Node_Cfg_Pair> aux;
// RDKit✔️✔️:
// RDKit✔️✔️:   auto &digraph = center->getDigraph();
// RDKit✔️✔️:   for (const auto &config : configs) {
// RDKit✔️✔️:     if (config == center) {
// RDKit✔️✔️:       continue;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     // FIXME: specific to each descriptor
// RDKit✔️✔️:     const auto &foci = config->getFoci();
// RDKit✔️✔️:     for (const auto &node : digraph.getNodes(foci[0])) {
// RDKit✔️✔️:       if (node->isDuplicate()) {
// RDKit✔️✔️:         continue;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       auto low = node;
// RDKit✔️✔️:       if (foci.size() == 2) {
// RDKit✔️✔️:         for (const auto &edge : node->getEdges(foci[1])) {
// RDKit✔️✔️:           const auto &other_node = edge->getOther(node);
// RDKit✔️✔️:           if (other_node->getDistance() < node->getDistance()) {
// RDKit✔️✔️:             low = other_node;
// RDKit✔️✔️:           }
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (!low->isDuplicate()) {
// RDKit✔️✔️:         aux.emplace_back(low, config.get());
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   auto farthest = [](const Node_Cfg_Pair &a, const Node_Cfg_Pair &b) {
// RDKit✔️✔️:     return a.first->getDistance() > b.first->getDistance();
// RDKit✔️✔️:   };
// RDKit✔️✔️:   std::sort(aux.begin(), aux.end(), farthest);
// RDKit✔️✔️:
// RDKit✔️✔️:   // Using a boost::unordered_map because it is more performant
// RDKit✔️✔️:   // than the STL version.
// RDKit✔️✔️:   boost::unordered_map<Node *, Descriptor> queue;
// RDKit✔️✔️:   int prev = std::numeric_limits<int>::max();
// RDKit✔️✔️:   for (const auto &e : aux) {
// RDKit✔️✔️:     const auto &node = e.first;
// RDKit✔️✔️:
// RDKit✔️✔️:     if (node->getDistance() < prev) {
// RDKit✔️✔️:       for (const auto &e2 : queue) {
// RDKit✔️✔️:         e2.first->setAux(e2.second);
// RDKit✔️✔️:       }
// RDKit✔️✔️:       queue.clear();
// RDKit✔️✔️:       prev = node->getDistance();
// RDKit✔️✔️:     }
// RDKit✔️✔️:     const auto &config = e.second;
// RDKit✔️✔️:     auto label = config->label(node, digraph, rules);
// RDKit✔️✔️:     queue.emplace(node, label);
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   for (const auto &e : queue) {
// RDKit✔️✔️:     e.first->setAux(e.second);
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   return true;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION labelAux
fn cip_label_aux(
    configs: &mut [CipConfig<'_>],
    rules: &CipRules,
    center_idx: usize,
    context: &mut CipLabelerContext,
) -> Result<bool, CipLabelerError> {
    let config_foci = configs
        .iter()
        .enumerate()
        .filter(|(idx, _)| *idx != center_idx)
        .map(|(idx, config)| (idx, config.get_foci().to_vec()))
        .collect::<Vec<_>>();
    let mut aux = Vec::<(CipNodeId, usize, i32)>::new();

    {
        let digraph = configs[center_idx].get_digraph_mut();
        for (config_idx, foci) in config_foci {
            for node in digraph.get_nodes(foci[0])? {
                if digraph.node(node).is_duplicate() {
                    continue;
                }
                let mut low = node;
                if foci.len() == 2 {
                    for edge in digraph.node_edges_for_atom(node, Some(foci[1]))? {
                        let other_node = digraph.edge(edge).get_other(edge, node)?;
                        if digraph.node(other_node).get_distance()
                            < digraph.node(node).get_distance()
                        {
                            low = other_node;
                        }
                    }
                }
                if !digraph.node(low).is_duplicate() {
                    aux.push((low, config_idx, digraph.node(low).get_distance()));
                }
            }
        }
    }

    aux.sort_by(|left, right| right.2.cmp(&left.2));

    let mut queue = Vec::<(CipNodeId, Descriptor)>::new();
    let mut prev = i32::MAX;
    for (node, config_idx, distance) in aux {
        if distance < prev {
            for (queued_node, desc) in queue.drain(..) {
                cip_set_center_node_aux(configs, center_idx, queued_node, desc);
            }
            prev = distance;
        }
        let label =
            cip_label_with_center_digraph(configs, center_idx, config_idx, node, rules, context)?;
        if !queue.iter().any(|(queued_node, _)| *queued_node == node) {
            queue.push((node, label));
        }
    }

    for (queued_node, desc) in queue {
        cip_set_center_node_aux(configs, center_idx, queued_node, desc);
    }

    Ok(true)
}

// BEGIN RDKIT CPP FUNCTION label (CIPLabeler.cpp)
// RDKit✔️✔️: void label(std::vector<std::unique_ptr<Configuration>> &configs,
// RDKit✔️✔️:            unsigned int maxRecursiveIterations) {
// RDKit✔️✔️:   // First, if the specified number of iterations allows it, run all centers
// RDKit✔️✔️:   // through a fast pass with the constitutional rules allow easy stuff to be
// RDKit✔️✔️:   // resolved.
// RDKit✔️✔️:   for (auto &conf : configs) {
// RDKit✔️✔️:     // Make sure this stereo center has no label
// RDKit✔️✔️:     conf->resetPrimaryLabel();
// RDKit✔️✔️:
// RDKit✔️✔️:     remainingCallCount = constitutionalRuleTimeout;
// RDKit✔️✔️:     try {
// RDKit✔️✔️:       auto desc = conf->label(constitutional_rules);
// RDKit✔️✔️:       if (desc != Descriptor::UNKNOWN) {
// RDKit✔️✔️:         conf->setPrimaryLabel(desc);
// RDKit✔️✔️:       }
// RDKit✔️✔️:     } catch (const MaxIterationsExceeded &) {
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   // Now, retry everything that hasn't been solved with a more generous
// RDKit✔️✔️:   // threshold
// RDKit✔️✔️:   if (maxRecursiveIterations != 0) {
// RDKit✔️✔️:     remainingCallCount = maxRecursiveIterations;
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     remainingCallCount = UINT_MAX;  // really big - will never be hit
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   // try again on everything that hasn't been resolved yet
// RDKit✔️✔️:   for (const auto &conf : configs) {
// RDKit✔️✔️:     if (conf->hasPrimaryLabel()) {
// RDKit✔️✔️:       // already resolved!
// RDKit✔️✔️:       continue;
// RDKit✔️✔️:     }
// RDKit✔️✔️:
// RDKit✔️✔️:     auto desc = conf->label(constitutional_rules);
// RDKit✔️✔️:     if (desc != Descriptor::UNKNOWN) {
// RDKit✔️✔️:       conf->setPrimaryLabel(desc);
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       if (labelAux(configs, all_rules, conf)) {
// RDKit✔️✔️:         desc = conf->label(all_rules);
// RDKit✔️✔️:
// RDKit✔️✔️:         if (desc != Descriptor::UNKNOWN) {
// RDKit✔️✔️:           conf->setPrimaryLabel(desc);
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION label
fn cip_label(
    configs: &mut [CipConfig<'_>],
    max_recursive_iterations: u32,
) -> Result<(), CipLabelerError> {
    let constitutional_rules = cip_constitutional_rules()?;
    for conf in configs.iter_mut() {
        conf.reset_primary_label();
        let mut context = CipLabelerContext::with_remaining_call_count(
            CipLabelerContext::CONSTITUTIONAL_RULE_TIMEOUT,
        );
        match conf.label(&constitutional_rules, &mut context) {
            Ok(desc) if desc != Descriptor::Unknown => conf.set_primary_label(desc)?,
            Ok(_) | Err(CipLabelerError::MaxIterationsExceeded) => {}
            Err(err) => return Err(err),
        }
    }

    let constitutional_rules = cip_constitutional_rules()?;
    let all_rules = cip_all_rules()?;
    let mut context = CipLabelerContext::new(max_recursive_iterations);
    for idx in 0..configs.len() {
        if configs[idx].has_primary_label() {
            continue;
        }

        let desc = configs[idx].label(&constitutional_rules, &mut context)?;
        if desc != Descriptor::Unknown {
            configs[idx].set_primary_label(desc)?;
        } else if cip_label_aux(configs, &all_rules, idx, &mut context)? {
            let desc = configs[idx].label(&all_rules, &mut context)?;
            if desc != Descriptor::Unknown {
                configs[idx].set_primary_label(desc)?;
            }
        }
    }
    Ok(())
}

fn cip_neighbor_order_value(indices: &[u32]) -> String {
    let body = indices
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join(",");
    format!("[{body}]")
}

fn cip_clear_selected_labels(molecule: &mut TopologyBlock, atom_mask: &[bool], bond_mask: &[bool]) {
    for (idx, selected) in atom_mask.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        if let Some(atom) = molecule.atoms.get_mut(idx) {
            if matches!(
                atom.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            ) {
                atom.clear_prop("_CIPCode");
            }
        }
    }
    for (idx, selected) in bond_mask.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        if let Some(bond) = molecule.bonds.get_mut(idx) {
            if matches!(
                bond.stereo(),
                BondStereo::E
                    | BondStereo::Z
                    | BondStereo::Cis
                    | BondStereo::Trans
                    | BondStereo::AtropCcw
                    | BondStereo::AtropCw
            ) {
                bond.clear_prop("_CIPCode");
            }
        }
    }
}

fn cip_apply_primary_labels(
    molecule: &mut TopologyBlock,
    labels: Vec<CipPrimaryLabel>,
) -> Result<(), CipLabelerError> {
    let topology = molecule;
    for label in labels {
        match label {
            CipPrimaryLabel::Atom(label) => {
                let atom_count = topology.atoms.len();
                let atom = topology.atoms.get_mut(label.atom_idx).ok_or(
                    CipLabelerError::AtomIndexOutOfRange {
                        index: label.atom_idx,
                        atom_count,
                    },
                )?;
                atom.set_prop("_CIPCode", label.cip_code)
                    .map_err(|source| CipLabelerError::AtomProperty {
                        atom: label.atom_idx,
                        source,
                    })?;
                atom.set_computed_prop(
                    "_CIPNeighborOrder",
                    cip_neighbor_order_value(&label.cip_neighbor_order),
                )
                .map_err(|source| CipLabelerError::AtomProperty {
                    atom: label.atom_idx,
                    source,
                })?;
            }
            CipPrimaryLabel::Bond(label) => {
                let bond_count = topology.bonds.len();
                let bond = topology.bonds.get_mut(label.bond_idx).ok_or(
                    CipLabelerError::BondIndexOutOfRange {
                        index: label.bond_idx,
                        bond_count,
                    },
                )?;
                bond.set_stereo_atoms(Some([
                    AtomId::new(label.stereo_atoms[0]),
                    AtomId::new(label.stereo_atoms[1]),
                ]));
                bond.set_stereo(label.stereo)
                    .map_err(|source| CipLabelerError::BondValue {
                        bond: label.bond_idx,
                        source,
                    })?;
                bond.set_prop("_CIPCode", label.cip_code)
                    .map_err(|source| CipLabelerError::BondValue {
                        bond: label.bond_idx,
                        source,
                    })?;
                bond.set_computed_prop(
                    "_CIPNeighborOrder",
                    cip_neighbor_order_value(&label.cip_neighbor_order),
                )
                .map_err(|source| CipLabelerError::BondValue {
                    bond: label.bond_idx,
                    source,
                })?;
            }
            CipPrimaryLabel::AtropisomerBond(label) => {
                let bond_count = topology.bonds.len();
                let bond = topology.bonds.get_mut(label.bond_idx).ok_or(
                    CipLabelerError::BondIndexOutOfRange {
                        index: label.bond_idx,
                        bond_count,
                    },
                )?;
                bond.set_prop("_CIPCode", label.cip_code)
                    .map_err(|source| CipLabelerError::BondValue {
                        bond: label.bond_idx,
                        source,
                    })?;
                bond.set_computed_prop(
                    "_CIPNeighborOrder",
                    cip_neighbor_order_value(&label.cip_neighbor_order),
                )
                .map_err(|source| CipLabelerError::BondValue {
                    bond: label.bond_idx,
                    source,
                })?;
            }
        }
    }
    Ok(())
}

// BEGIN RDKIT CPP FUNCTION assignCIPLabels selected overload (CIPLabeler.cpp)
// RDKit✔️✔️: void assignCIPLabels(ROMol &mol, const boost::dynamic_bitset<> &atoms,
// RDKit✔️✔️:                      const boost::dynamic_bitset<> &bonds,
// RDKit✔️✔️:                      unsigned int maxRecursiveIterations) {
// RDKit✔️✔️:   ControlCHandler::reset();
// RDKit✔️✔️:
// RDKit✔️✔️:   // reset the mark, for the case that this fails
// RDKit✔️✔️:   mol.clearProp(common_properties::_CIPComputed);
// RDKit✔️✔️:   CIPMol cipmol{mol};
// RDKit✔️✔️:   auto configs = findConfigs(cipmol, atoms, bonds);
// RDKit✔️✔️:
// RDKit✔️✔️:   try {
// RDKit✔️✔️:     label(configs, maxRecursiveIterations);
// RDKit❌❌:   } catch (const ControlCCaught &) {
// RDKit❌❌:   }
// RDKit❌❌:   if (ControlCHandler::getGotSignal()) {
// RDKit❌❌:     BOOST_LOG(rdWarningLog)
// RDKit❌❌:         << "Interrupted, cancelling CIP label calculation" << std::endl;
// RDKit❌❌:     return;
// RDKit❌❌:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   const bool computed = true;
// RDKit✔️✔️:   mol.setProp(common_properties::_CIPComputed, true, computed);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION assignCIPLabels selected overload
fn assign_cip_labels_for_masks(
    mut topology: TopologyBlock,
    mut properties: MoleculeProperties,
    atom_mask: &[bool],
    bond_mask: &[bool],
    max_recursive_iterations: u32,
) -> Result<CipLabelAssignment, CipLabelerError> {
    properties.clear_prop("_CIPComputed");
    cip_clear_selected_labels(&mut topology, atom_mask, bond_mask);

    let mut configs = cip_find_configs(&topology, atom_mask, bond_mask)?;
    cip_label(&mut configs, max_recursive_iterations)?;
    let labels = configs
        .iter()
        .filter_map(CipConfig::primary_label)
        .collect::<Vec<_>>();
    drop(configs);

    cip_apply_primary_labels(&mut topology, labels)?;
    topology.validate()?;
    properties.set_computed_prop("_CIPComputed", "true")?;
    Ok(CipLabelAssignment {
        topology,
        properties,
    })
}

// BEGIN RDKIT CPP FUNCTION assignCIPLabels all-molecule overload (CIPLabeler.cpp)
// RDKit✔️✔️: void assignCIPLabels(ROMol &mol, unsigned int maxRecursiveIterations) {
// RDKit✔️✔️:   boost::dynamic_bitset<> atoms(mol.getNumAtoms());
// RDKit✔️✔️:   boost::dynamic_bitset<> bonds(mol.getNumBonds());
// RDKit✔️✔️:   atoms.set();
// RDKit✔️✔️:   bonds.set();
// RDKit✔️✔️:   assignCIPLabels(mol, atoms, bonds, maxRecursiveIterations);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION assignCIPLabels all-molecule overload
pub fn assign_cip_labels(
    topology: TopologyBlock,
    properties: MoleculeProperties,
    options: &CipLabelOptions,
) -> Result<CipLabelAssignment, CipLabelerError> {
    topology.validate()?;
    if topology.atoms.len() > u32::MAX as usize {
        return Err(CipLabelerError::SourceIndexWidthExceeded {
            kind: "atom count",
            index: topology.atoms.len(),
        });
    }
    if topology.bonds.len() > u32::MAX as usize {
        return Err(CipLabelerError::SourceIndexWidthExceeded {
            kind: "bond count",
            index: topology.bonds.len(),
        });
    }

    let atom_selection = options.atoms.as_deref().unwrap_or_default();
    let bond_selection = options.bonds.as_deref().unwrap_or_default();
    for atom in atom_selection {
        if atom.index() >= topology.atoms.len() {
            return Err(CipLabelerError::AtomIndexOutOfRange {
                index: atom.index(),
                atom_count: topology.atoms.len(),
            });
        }
    }
    for bond in bond_selection {
        if bond.index() >= topology.bonds.len() {
            return Err(CipLabelerError::BondIndexOutOfRange {
                index: bond.index(),
                bond_count: topology.bonds.len(),
            });
        }
    }

    let select_all = atom_selection.is_empty() && bond_selection.is_empty();
    let mut atom_mask = vec![select_all; topology.atoms.len()];
    let mut bond_mask = vec![select_all; topology.bonds.len()];
    if !select_all {
        for atom in atom_selection {
            atom_mask[atom.index()] = true;
        }
        for bond in bond_selection {
            bond_mask[bond.index()] = true;
        }
    }
    assign_cip_labels_for_masks(
        topology,
        properties,
        &atom_mask,
        &bond_mask,
        options.max_recursive_iterations,
    )
}

pub(crate) struct CipConfiguration<'a> {
    foci: Vec<usize>,
    carriers: Vec<Option<usize>>,
    digraph: CipDigraph<'a>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CipAtomPrimaryLabel {
    atom_idx: usize,
    cip_code: &'static str,
    cip_neighbor_order: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CipBondPrimaryLabel {
    bond_idx: usize,
    stereo_atoms: [usize; 2],
    stereo: BondStereo,
    cip_code: &'static str,
    cip_neighbor_order: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CipAtropisomerBondPrimaryLabel {
    bond_idx: usize,
    cip_code: &'static str,
    cip_neighbor_order: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum CipPrimaryLabel {
    Atom(CipAtomPrimaryLabel),
    Bond(CipBondPrimaryLabel),
    AtropisomerBond(CipAtropisomerBondPrimaryLabel),
}

pub(crate) enum CipConfig<'a> {
    Tetrahedral(CipTetrahedral<'a>),
    Sp2Bond(CipSp2Bond<'a>),
    AtropisomerBond(CipAtropisomerBond<'a>),
}

impl<'a> CipConfig<'a> {
    fn get_foci(&self) -> &[usize] {
        match self {
            Self::Tetrahedral(config) => config.get_foci(),
            Self::Sp2Bond(config) => config.get_foci(),
            Self::AtropisomerBond(config) => config.get_foci(),
        }
    }

    fn get_digraph_mut(&mut self) -> &mut CipDigraph<'a> {
        match self {
            Self::Tetrahedral(config) => config.configuration.get_digraph(),
            Self::Sp2Bond(config) => config.configuration.get_digraph(),
            Self::AtropisomerBond(config) => config.configuration.get_digraph(),
        }
    }

    fn reset_primary_label(&mut self) {
        match self {
            Self::Tetrahedral(config) => config.reset_primary_label(),
            Self::Sp2Bond(config) => config.reset_primary_label(),
            Self::AtropisomerBond(config) => config.reset_primary_label(),
        }
    }

    fn has_primary_label(&self) -> bool {
        match self {
            Self::Tetrahedral(config) => config.has_primary_label(),
            Self::Sp2Bond(config) => config.has_primary_label(),
            Self::AtropisomerBond(config) => config.has_primary_label(),
        }
    }

    fn label(
        &mut self,
        rules: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        match self {
            Self::Tetrahedral(config) => config.label(rules, context),
            Self::Sp2Bond(config) => config.label(rules, context),
            Self::AtropisomerBond(config) => config.label(rules, context),
        }
    }

    fn label_with_external_digraph(
        &mut self,
        node: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        rules: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        match self {
            Self::Tetrahedral(config) => {
                config.label_with_external_digraph(node, digraph, rules, context)
            }
            Self::Sp2Bond(config) => {
                config.label_with_external_digraph(node, digraph, rules, context)
            }
            Self::AtropisomerBond(config) => {
                config.label_with_external_digraph(node, digraph, rules, context)
            }
        }
    }

    fn set_primary_label(&mut self, desc: Descriptor) -> Result<(), CipLabelerError> {
        match self {
            Self::Tetrahedral(config) => config.set_primary_label(desc),
            Self::Sp2Bond(config) => config.set_primary_label(desc),
            Self::AtropisomerBond(config) => config.set_primary_label(desc),
        }
    }

    fn primary_label(&self) -> Option<CipPrimaryLabel> {
        match self {
            Self::Tetrahedral(config) => config.primary_label().cloned().map(CipPrimaryLabel::Atom),
            Self::Sp2Bond(config) => config.primary_label().cloned().map(CipPrimaryLabel::Bond),
            Self::AtropisomerBond(config) => config
                .primary_label()
                .cloned()
                .map(CipPrimaryLabel::AtropisomerBond),
        }
    }
}

impl<'a> CipConfiguration<'a> {
    // BEGIN RDKIT CPP FUNCTION Configuration::parity4 (configs/Configuration.h)
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: static int parity4(const std::vector<T> &trg, const std::vector<T> &ref) {
    // RDKit✔️✔️:   if (ref.size() != 4 || trg.size() != ref.size()) {
    // RDKit✔️✔️:     throw std::runtime_error("Parity vectors must have size 4.");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (ref[0] == trg[0]) {
    // RDKit✔️✔️:     if (ref[1] == trg[1]) {
    // RDKit✔️✔️:       // a,b,c,d -> a,b,c,d
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> a,b,d,c
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[2]) {
    // RDKit✔️✔️:       // a,b,c,d -> a,c,b,d
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> a,c,d,b
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[3]) {
    // RDKit✔️✔️:       // a,b,c,d -> a,d,c,b
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> a,d,b,c
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (ref[0] == trg[1]) {
    // RDKit✔️✔️:     if (ref[1] == trg[0]) {
    // RDKit✔️✔️:       // a,b,c,d -> b,a,c,d
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> b,a,d,c
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[2]) {
    // RDKit✔️✔️:       // a,b,c,d -> b,c,a,d
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> b,c,d,a
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[3]) {
    // RDKit✔️✔️:       // a,b,c,d -> b,d,c,a
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> b,d,a,c
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (ref[0] == trg[2]) {
    // RDKit✔️✔️:     if (ref[1] == trg[1]) {
    // RDKit✔️✔️:       // a,b,c,d -> c,b,a,d
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> c,b,d,a
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[0]) {
    // RDKit✔️✔️:       // a,b,c,d -> c,a,b,d
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> c,a,d,b
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[3]) {
    // RDKit✔️✔️:       // a,b,c,d -> c,d,a,b
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> c,d,b,a
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (ref[0] == trg[3]) {
    // RDKit✔️✔️:     if (ref[1] == trg[1]) {
    // RDKit✔️✔️:       // a,b,c,d -> d,b,c,a
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> d,b,a,c
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[2]) {
    // RDKit✔️✔️:       // a,b,c,d -> d,c,b,a
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> d,c,a,b
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[0]) {
    // RDKit✔️✔️:       // a,b,c,d -> d,a,c,b
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> d,a,b,c
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // We should never hit this, but the compiler still complains
    // RDKit✔️✔️:   // about a missing return statement.
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Configuration::parity4
    pub(crate) fn parity4<T: PartialEq>(
        trg: &[T],
        reference: &[T],
    ) -> Result<i32, CipLabelerError> {
        if reference.len() != 4 || trg.len() != reference.len() {
            return Err(CipLabelerError::ParityVectorsMustHaveSize4);
        }

        let r = reference;
        let t = trg;
        if r[0] == t[0] {
            if r[1] == t[1] {
                if r[2] == t[2] && r[3] == t[3] {
                    return Ok(2);
                }
                if r[2] == t[3] && r[3] == t[2] {
                    return Ok(1);
                }
            } else if r[1] == t[2] {
                if r[2] == t[1] && r[3] == t[3] {
                    return Ok(1);
                }
                if r[2] == t[3] && r[3] == t[1] {
                    return Ok(2);
                }
            } else if r[1] == t[3] {
                if r[2] == t[2] && r[3] == t[1] {
                    return Ok(1);
                }
                if r[2] == t[1] && r[3] == t[2] {
                    return Ok(2);
                }
            }
        } else if r[0] == t[1] {
            if r[1] == t[0] {
                if r[2] == t[2] && r[3] == t[3] {
                    return Ok(1);
                }
                if r[2] == t[3] && r[3] == t[2] {
                    return Ok(2);
                }
            } else if r[1] == t[2] {
                if r[2] == t[0] && r[3] == t[3] {
                    return Ok(2);
                }
                if r[2] == t[3] && r[3] == t[0] {
                    return Ok(1);
                }
            } else if r[1] == t[3] {
                if r[2] == t[2] && r[3] == t[0] {
                    return Ok(2);
                }
                if r[2] == t[0] && r[3] == t[2] {
                    return Ok(1);
                }
            }
        } else if r[0] == t[2] {
            if r[1] == t[1] {
                if r[2] == t[0] && r[3] == t[3] {
                    return Ok(1);
                }
                if r[2] == t[3] && r[3] == t[0] {
                    return Ok(2);
                }
            } else if r[1] == t[0] {
                if r[2] == t[1] && r[3] == t[3] {
                    return Ok(2);
                }
                if r[2] == t[3] && r[3] == t[1] {
                    return Ok(1);
                }
            } else if r[1] == t[3] {
                if r[2] == t[0] && r[3] == t[1] {
                    return Ok(2);
                }
                if r[2] == t[1] && r[3] == t[0] {
                    return Ok(1);
                }
            }
        } else if r[0] == t[3] {
            if r[1] == t[1] {
                if r[2] == t[2] && r[3] == t[0] {
                    return Ok(1);
                }
                if r[2] == t[0] && r[3] == t[2] {
                    return Ok(2);
                }
            } else if r[1] == t[2] {
                if r[2] == t[1] && r[3] == t[0] {
                    return Ok(2);
                }
                if r[2] == t[0] && r[3] == t[1] {
                    return Ok(1);
                }
            } else if r[1] == t[0] {
                if r[2] == t[2] && r[3] == t[1] {
                    return Ok(2);
                }
                if r[2] == t[1] && r[3] == t[2] {
                    return Ok(1);
                }
            }
        }
        Ok(0)
    }

    // BEGIN RDKIT CPP FUNCTION Configuration constructors/getters/setCarriers/label (configs/Configuration.cpp)
    // RDKit✔️✔️: Configuration::Configuration(const CIPMol &mol, Atom *focus)
    // RDKit✔️✔️:     : d_foci{focus}, d_digraph{mol, focus} {};
    // RDKit✔️✔️:
    // RDKit✔️✔️: Configuration::Configuration(const CIPMol &mol, std::vector<Atom *> &&foci,
    // RDKit✔️✔️:                              bool atropisomerMode)
    // RDKit✔️✔️:     : d_foci{std::move(foci)}, d_digraph{mol, d_foci[0], atropisomerMode} {}
    // RDKit✔️✔️:
    // RDKit✔️✔️: Configuration::~Configuration() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: Atom *Configuration::getFocus() const { return d_foci[0]; }
    // RDKit✔️✔️:
    // RDKit✔️✔️: const std::vector<Atom *> &Configuration::getFoci() const { return d_foci; }
    // RDKit✔️✔️:
    // RDKit✔️✔️: const std::vector<Atom *> &Configuration::getCarriers() const {
    // RDKit✔️✔️:   return d_carriers;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: Digraph &Configuration::getDigraph() { return d_digraph; }
    // RDKit✔️✔️:
    // RDKit✔️✔️: Descriptor Configuration::label(Node *node, Digraph &digraph,
    // RDKit✔️✔️:                                 const Rules &comp) {
    // RDKit✔️✔️:   (void)node;
    // RDKit✔️✔️:   (void)digraph;
    // RDKit✔️✔️:   (void)comp;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return Descriptor::UNKNOWN;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Configuration constructors/getters/setCarriers/label
    pub(crate) fn new(molecule: &'a TopologyBlock, focus: usize) -> Result<Self, CipLabelerError> {
        Ok(Self {
            foci: vec![focus],
            carriers: Vec::new(),
            digraph: CipDigraph::new(molecule, focus, false)?,
        })
    }

    pub(crate) fn with_foci(
        molecule: &'a TopologyBlock,
        foci: Vec<usize>,
        atropisomer_mode: bool,
    ) -> Result<Self, CipLabelerError> {
        let focus = *foci
            .first()
            .ok_or(CipLabelerError::EmptyConfigurationFoci)?;
        Ok(Self {
            foci,
            carriers: Vec::new(),
            digraph: CipDigraph::new(molecule, focus, atropisomer_mode)?,
        })
    }

    pub(crate) fn get_focus(&self) -> usize {
        self.foci[0]
    }

    pub(crate) fn get_foci(&self) -> &[usize] {
        &self.foci
    }

    pub(crate) fn set_carriers(&mut self, carriers: Vec<Option<usize>>) {
        self.carriers = carriers;
    }

    pub(crate) fn get_carriers(&self) -> &[Option<usize>] {
        &self.carriers
    }

    pub(crate) fn get_digraph(&mut self) -> &mut CipDigraph<'a> {
        &mut self.digraph
    }

    pub(crate) fn label(
        &self,
        _node: CipNodeId,
        _digraph: &mut CipDigraph<'_>,
        _comp: &CipRules,
    ) -> Descriptor {
        Descriptor::Unknown
    }

    // BEGIN RDKIT CPP FUNCTION Configuration::findInternalEdge/isInternalEdge/removeInternalEdges (configs/Configuration.cpp)
    // RDKit✔️✔️: Edge *Configuration::findInternalEdge(const std::vector<Edge *> &edges,
    // RDKit✔️✔️:                                       Atom *f1, Atom *f2) {
    // RDKit✔️✔️:   for (const auto &edge : edges) {
    // RDKit✔️✔️:     if (edge->getBeg()->isDuplicate() || edge->getEnd()->isDuplicate()) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (isInternalEdge(edge, f1, f2)) {
    // RDKit✔️✔️:       return edge;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: bool Configuration::isInternalEdge(const Edge *edge, Atom *f1, Atom *f2) {
    // RDKit✔️✔️:   const auto &beg = edge->getBeg();
    // RDKit✔️✔️:   const auto &end = edge->getEnd();
    // RDKit✔️✔️:   if (f1 == beg->getAtom() && f2 == end->getAtom()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else if (f1 == end->getAtom() && f2 == beg->getAtom()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: void Configuration::removeInternalEdges(std::vector<Edge *> &edges, Atom *f1,
    // RDKit✔️✔️:                                         Atom *f2) {
    // RDKit✔️✔️:   std::vector<Edge *> new_edges;
    // RDKit✔️✔️:   for (auto &&e : edges) {
    // RDKit✔️✔️:     if (!isInternalEdge(e, f1, f2)) {
    // RDKit✔️✔️:       new_edges.push_back(std::move(e));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::swap(edges, new_edges);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Configuration::findInternalEdge/isInternalEdge/removeInternalEdges
    fn find_internal_edge(
        digraph: &CipDigraph<'_>,
        edges: &[CipEdgeId],
        f1: usize,
        f2: usize,
    ) -> Option<CipEdgeId> {
        edges.iter().copied().find(|edge_id| {
            let edge = digraph.edge(*edge_id);
            if digraph.node(edge.get_beg()).is_duplicate()
                || digraph.node(edge.get_end()).is_duplicate()
            {
                return false;
            }
            Self::is_internal_edge(digraph, *edge_id, f1, f2)
        })
    }

    fn is_internal_edge(
        digraph: &CipDigraph<'_>,
        edge_id: CipEdgeId,
        f1: usize,
        f2: usize,
    ) -> bool {
        let edge = digraph.edge(edge_id);
        let beg = digraph.node(edge.get_beg()).atom_idx();
        let end = digraph.node(edge.get_end()).atom_idx();
        (beg == Some(f1) && end == Some(f2)) || (beg == Some(f2) && end == Some(f1))
    }

    fn remove_internal_edges(
        digraph: &CipDigraph<'_>,
        edges: &mut Vec<CipEdgeId>,
        f1: usize,
        f2: usize,
    ) {
        edges.retain(|edge_id| !Self::is_internal_edge(digraph, *edge_id, f1, f2));
    }

    // BEGIN RDKIT CPP FUNCTION Configuration::isDuplicateOrHydrogenEdge/removeDuplicatesAndHs (configs/Configuration.cpp)
    // RDKit✔️✔️: bool Configuration::isDuplicateOrHydrogenEdge(const Edge *edge) {
    // RDKit✔️✔️:   return edge->getBeg()->isDuplicateOrH() || edge->getEnd()->isDuplicateOrH();
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: void Configuration::removeDuplicatesAndHs(std::vector<Edge *> &edges) {
    // RDKit✔️✔️:   std::vector<Edge *> new_edges;
    // RDKit✔️✔️:   for (auto &&e : edges) {
    // RDKit✔️✔️:     if (!isDuplicateOrHydrogenEdge(e)) {
    // RDKit✔️✔️:       new_edges.push_back(std::move(e));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   edges = std::move(new_edges);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Configuration::isDuplicateOrHydrogenEdge/removeDuplicatesAndHs
    fn is_duplicate_or_hydrogen_edge(digraph: &CipDigraph<'_>, edge_id: CipEdgeId) -> bool {
        let edge = digraph.edge(edge_id);
        digraph.node(edge.get_beg()).is_duplicate_or_h()
            || digraph.node(edge.get_end()).is_duplicate_or_h()
    }

    fn remove_duplicates_and_hs(digraph: &CipDigraph<'_>, edges: &mut Vec<CipEdgeId>) {
        edges.retain(|edge_id| !Self::is_duplicate_or_hydrogen_edge(digraph, *edge_id));
    }
}

pub(crate) struct CipTetrahedral<'a> {
    configuration: CipConfiguration<'a>,
    ranked_anchors: Vec<u32>,
    primary_label: Option<CipAtomPrimaryLabel>,
}

impl<'a> CipTetrahedral<'a> {
    // BEGIN RDKIT CPP FUNCTION Tetrahedral::Tetrahedral (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: Tetrahedral::Tetrahedral(const CIPMol &mol, Atom *focus)
    // RDKit✔️✔️:     : Configuration(mol, focus) {
    // RDKit✔️✔️:   CHECK_INVARIANT(focus, "bad atom")
    // RDKit✔️✔️:   CHECK_INVARIANT(focus->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW ||
    // RDKit✔️✔️:                       focus->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW,
    // RDKit✔️✔️:                   "bad config")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Atom *> carriers;
    // RDKit✔️✔️:   carriers.reserve(4);
    // RDKit✔️✔️:   for (auto &nbr : mol.getNeighbors(focus)) {
    // RDKit✔️✔️:     carriers.push_back(nbr);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (carriers.size() < 4) {
    // RDKit✔️✔️:     // Implicit H -- use the central atom instead of a dummy H
    // RDKit✔️✔️:     carriers.push_back(focus);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (carriers.size() < 4) {
    // RDKit✔️✔️:     // Trigonal pyramid centers with an implicit H need a phantom
    // RDKit✔️✔️:     // atom as fourth carrier. This one must be represented differently
    // RDKit✔️✔️:     // than the implicit H.
    // RDKit✔️✔️:     carriers.push_back(nullptr);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   POSTCONDITION(carriers.size() == 4, "configuration must have 4 carriers");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   setCarriers(std::move(carriers));
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION Tetrahedral::Tetrahedral
    pub(crate) fn new(molecule: &'a TopologyBlock, focus: usize) -> Result<Self, CipLabelerError> {
        let mut configuration = CipConfiguration::new(molecule, focus)?;
        let atom = configuration.digraph.mol().atom(focus)?;
        match atom.chiral_tag() {
            ChiralTag::TetrahedralCcw | ChiralTag::TetrahedralCw => {}
            ChiralTag::Unspecified
            | ChiralTag::Other
            | ChiralTag::Tetrahedral
            | ChiralTag::Allene
            | ChiralTag::SquarePlanar
            | ChiralTag::TrigonalBipyramidal
            | ChiralTag::Octahedral => return Err(CipLabelerError::BadTetrahedralConfig),
        }

        let mut carriers = Vec::with_capacity(4);
        for nbr in configuration.digraph.mol().neighbor_indices(focus)? {
            carriers.push(Some(nbr));
        }
        if carriers.len() < 4 {
            carriers.push(Some(focus));
        }
        if carriers.len() < 4 {
            carriers.push(None);
        }
        if carriers.len() != 4 {
            return Err(CipLabelerError::TetrahedralConfigurationMustHave4Carriers);
        }

        configuration.set_carriers(carriers);
        Ok(Self {
            configuration,
            ranked_anchors: Vec::new(),
            primary_label: None,
        })
    }

    pub(crate) fn get_focus(&self) -> usize {
        self.configuration.get_focus()
    }

    pub(crate) fn get_foci(&self) -> &[usize] {
        self.configuration.get_foci()
    }

    pub(crate) fn get_carriers(&self) -> &[Option<usize>] {
        self.configuration.get_carriers()
    }

    pub(crate) fn ranked_anchors(&self) -> &[u32] {
        &self.ranked_anchors
    }

    pub(crate) fn primary_label(&self) -> Option<&CipAtomPrimaryLabel> {
        self.primary_label.as_ref()
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::setPrimaryLabel (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: void Tetrahedral::setPrimaryLabel(Descriptor desc) {
    // RDKit✔️✔️:   switch (desc) {
    // RDKit✔️✔️:     case Descriptor::R:
    // RDKit✔️✔️:     case Descriptor::S:
    // RDKit✔️✔️:     case Descriptor::r:
    // RDKit✔️✔️:     case Descriptor::s: {
    // RDKit✔️✔️:       auto chiralAtom = getFocus();
    // RDKit✔️✔️:       chiralAtom->setProp(common_properties::_CIPCode, to_string(desc));
    // RDKit✔️✔️:       chiralAtom->setProp(common_properties::_CIPNeighborOrder,
    // RDKit✔️✔️:                           d_ranked_anchors, true);
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     case Descriptor::seqTrans:
    // RDKit✔️✔️:     case Descriptor::seqCis:
    // RDKit✔️✔️:     case Descriptor::E:
    // RDKit✔️✔️:     case Descriptor::Z:
    // RDKit✔️✔️:     case Descriptor::M:
    // RDKit✔️✔️:     case Descriptor::P:
    // RDKit✔️✔️:     case Descriptor::m:
    // RDKit✔️✔️:     case Descriptor::p:
    // RDKit✔️✔️:     case Descriptor::SP_4:
    // RDKit✔️✔️:     case Descriptor::TBPY_5:
    // RDKit✔️✔️:     case Descriptor::OC_6:
    // RDKit✔️✔️:       throw std::runtime_error(
    // RDKit✔️✔️:           "Received a Descriptor that is not supported for atoms");
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       throw std::runtime_error("Received an invalid Atom Descriptor");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::setPrimaryLabel
    pub(crate) fn set_primary_label(&mut self, desc: Descriptor) -> Result<(), CipLabelerError> {
        match desc {
            Descriptor::R | Descriptor::S | Descriptor::r | Descriptor::s => {
                self.primary_label = Some(CipAtomPrimaryLabel {
                    atom_idx: self.configuration.get_focus(),
                    cip_code: descriptor_to_string(desc),
                    cip_neighbor_order: self.ranked_anchors.clone(),
                });
                Ok(())
            }
            Descriptor::seqTrans
            | Descriptor::seqCis
            | Descriptor::E
            | Descriptor::Z
            | Descriptor::M
            | Descriptor::P
            | Descriptor::m
            | Descriptor::p
            | Descriptor::SP_4
            | Descriptor::TBPY_5
            | Descriptor::OC_6 => Err(CipLabelerError::DescriptorNotSupportedForAtoms),
            Descriptor::None | Descriptor::Unknown | Descriptor::ns => {
                Err(CipLabelerError::InvalidAtomDescriptor)
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::hasPrimaryLabel (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: bool Tetrahedral::hasPrimaryLabel() const {
    // RDKit✔️✔️:   return getFocus()->hasProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::hasPrimaryLabel
    pub(crate) fn has_primary_label(&self) -> bool {
        self.primary_label.is_some()
            || self
                .configuration
                .digraph
                .mol()
                .atom(self.configuration.get_focus())
                .is_ok_and(|atom| atom.prop("_CIPCode").is_some())
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::resetPrimaryLabel (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: void Tetrahedral::resetPrimaryLabel() const {
    // RDKit✔️✔️:   getFocus()->clearProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::resetPrimaryLabel
    pub(crate) fn reset_primary_label(&mut self) {
        self.primary_label = None;
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::label(const Rules &) (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: Descriptor Tetrahedral::label(const Rules &comp) {
    // RDKit✔️✔️:   auto &digraph = getDigraph();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto root = digraph.getOriginalRoot();
    // RDKit✔️✔️:   if (digraph.getCurrentRoot() != root) {
    // RDKit✔️✔️:     digraph.changeRoot(root);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return label(root, comp);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::label(const Rules &)
    pub(crate) fn label(
        &mut self,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let root = self.configuration.digraph.get_original_root();
        if self.configuration.digraph.get_current_root() != root {
            self.configuration.digraph.change_root(root)?;
        }
        self.label_node(root, comp, context)
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::label(Node *, Digraph &, const Rules &) (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: Descriptor Tetrahedral::label(Node *node, Digraph &digraph, const Rules &comp) {
    // RDKit✔️✔️:   digraph.changeRoot(node);
    // RDKit✔️✔️:   return label(node, comp);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::label(Node *, Digraph &, const Rules &)
    pub(crate) fn label_with_external_digraph(
        &mut self,
        node: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        digraph.change_root(node)?;
        self.label_node_in_digraph(node, digraph, comp, context)
    }

    fn label_node(
        &mut self,
        node: CipNodeId,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let focus = self.configuration.get_focus();
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            focus,
            &carriers,
            &mut self.ranked_anchors,
            node,
            &mut self.configuration.digraph,
            comp,
            context,
        )
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::label(Node *, const Rules &) (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: Descriptor Tetrahedral::label(Node *node, const Rules &comp) {
    // RDKit✔️✔️:   auto focus = getFocus();
    // RDKit✔️✔️:   auto edges = node->getEdges();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   d_ranked_anchors.clear();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // something not right!?! bad creation
    // RDKit✔️✔️:   if (edges.size() < 3) {
    // RDKit✔️✔️:     return Descriptor::ns;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto priority = comp.sort(node, edges);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool isUnique = priority.isUnique();
    // RDKit✔️✔️:   if (!isUnique && edges.size() == 4) {
    // RDKit✔️✔️:     if (comp.getNumSubRules() == 3) {
    // RDKit✔️✔️:       return Descriptor::UNKNOWN;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto partition = comp.getSorter()->getGroups(edges);
    // RDKit✔️✔️:     if (partition.size() == 2) {
    // RDKit✔️✔️:       node->getDigraph()->setRule6Ref(edges[1]->getEnd()->getAtom());
    // RDKit✔️✔️:       priority = comp.sort(node, edges);
    // RDKit✔️✔️:       node->getDigraph()->setRule6Ref(nullptr);
    // RDKit✔️✔️:     } else if (partition.size() == 1) {
    // RDKit✔️✔️:       // S4 symmetric case
    // RDKit✔️✔️:       node->getDigraph()->setRule6Ref(edges[0]->getEnd()->getAtom());
    // RDKit✔️✔️:       comp.sort(node, edges);
    // RDKit✔️✔️:       auto nbrs1 = std::vector<Edge *>(edges.begin(), edges.end());
    // RDKit✔️✔️:
    // RDKit✔️✔️:       node->getDigraph()->setRule6Ref(edges[1]->getEnd()->getAtom());
    // RDKit✔️✔️:       priority = comp.sort(node, edges);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       node->getDigraph()->setRule6Ref(nullptr);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (parity4(nbrs1, edges) == 1) {
    // RDKit✔️✔️:         return Descriptor::UNKNOWN;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!priority.isUnique()) {
    // RDKit✔️✔️:       return Descriptor::UNKNOWN;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (!isUnique) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto ordered = std::vector<Atom *>(4, nullptr);
    // RDKit✔️✔️:   int idx = 0;
    // RDKit✔️✔️:   d_ranked_anchors.reserve(4);
    // RDKit✔️✔️:   for (const auto &edge : edges) {
    // RDKit✔️✔️:     if (edge->getEnd()->isSet(Node::BOND_DUPLICATE) ||
    // RDKit✔️✔️:         edge->getEnd()->isSet(Node::IMPL_HYDROGEN)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto atom = edge->getEnd()->getAtom();
    // RDKit✔️✔️:     ordered[idx] = atom;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // In this case we don't worry about implicit H (see Sp2Bond
    // RDKit✔️✔️:     // and Atropisomer): chirality is positional, and we don't
    // RDKit✔️✔️:     // know where the implicit H may be ("before" or "after" a
    // RDKit✔️✔️:     // potential 1H with lower priority?), so we just ignore it
    // RDKit✔️✔️:     // in the ranked neighbors list
    // RDKit✔️✔️:     d_ranked_anchors.push_back(atom->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     ++idx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // if we are resolving a trigonal pyramid with an implicit H,
    // RDKit✔️✔️:   // the 4th carrier will be a nullptr: we need to add a phantom
    // RDKit✔️✔️:   // atom, which will always have the lowest priority, so that
    // RDKit✔️✔️:   // it must be different than the representation of the implicit H.
    // RDKit✔️✔️:   if (idx < 4) {
    // RDKit✔️✔️:     ordered[idx] = focus;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int parity = parity4(ordered, getCarriers());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (parity == 0) {
    // RDKit✔️✔️:     throw std::runtime_error("Could not calculate parity! Carrier mismatch");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto config = focus->getChiralTag();
    // RDKit✔️✔️:   if (parity == 1) {
    // RDKit✔️✔️:     if (config == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:       config = Atom::CHI_TETRAHEDRAL_CW;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       config = Atom::CHI_TETRAHEDRAL_CCW;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (config == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:     if (priority.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::s;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::S;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (config == Atom::CHI_TETRAHEDRAL_CW) {
    // RDKit✔️✔️:     if (priority.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::r;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::R;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return Descriptor::UNKNOWN;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::label(Node *, const Rules &)
    fn label_node_in_digraph(
        &mut self,
        node: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let focus = self.configuration.get_focus();
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            focus,
            &carriers,
            &mut self.ranked_anchors,
            node,
            digraph,
            comp,
            context,
        )
    }

    fn label_node_impl(
        focus: usize,
        carriers: &[Option<usize>],
        ranked_anchors: &mut Vec<u32>,
        node: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let mut edges = digraph.node_edges(node)?;
        ranked_anchors.clear();

        if edges.len() < 3 {
            return Ok(Descriptor::ns);
        }

        let mut priority = comp.sort(digraph, context, node, &mut edges, true)?;
        let is_unique = priority.is_unique();
        if !is_unique && edges.len() == 4 {
            if comp.get_num_sub_rules() == 3 {
                return Ok(Descriptor::Unknown);
            }
            let partition = comp.get_sorter().get_groups(digraph, context, &edges)?;
            if partition.len() == 2 {
                let ref_atom = digraph.node(digraph.edge(edges[1]).get_end()).atom_idx();
                digraph.set_rule6_ref(ref_atom)?;
                priority = comp.sort(digraph, context, node, &mut edges, true)?;
                digraph.set_rule6_ref(None)?;
            } else if partition.len() == 1 {
                let ref_atom = digraph.node(digraph.edge(edges[0]).get_end()).atom_idx();
                digraph.set_rule6_ref(ref_atom)?;
                comp.sort(digraph, context, node, &mut edges, true)?;
                let nbrs1 = edges.clone();

                let ref_atom = digraph.node(digraph.edge(edges[1]).get_end()).atom_idx();
                digraph.set_rule6_ref(ref_atom)?;
                priority = comp.sort(digraph, context, node, &mut edges, true)?;

                digraph.set_rule6_ref(None)?;

                if CipConfiguration::parity4(&nbrs1, &edges)? == 1 {
                    return Ok(Descriptor::Unknown);
                }
            }
            if !priority.is_unique() {
                return Ok(Descriptor::Unknown);
            }
        } else if !is_unique {
            return Ok(Descriptor::Unknown);
        }

        let mut ordered = vec![None; 4];
        let mut idx = 0_usize;
        ranked_anchors.reserve(4);
        for edge in &edges {
            let end = digraph.edge(*edge).get_end();
            if digraph.node(end).is_set(CipNode::BOND_DUPLICATE)
                || digraph.node(end).is_set(CipNode::IMPL_HYDROGEN)
            {
                continue;
            }

            let atom = digraph.node(end).atom_idx();
            if idx < 4 {
                ordered[idx] = atom;
            }
            if let Some(atom_idx) = atom {
                ranked_anchors.push(u32::try_from(atom_idx).map_err(|_| {
                    CipLabelerError::SourceIndexWidthExceeded {
                        kind: "atom",
                        index: atom_idx,
                    }
                })?);
            }
            idx += 1;
        }

        if idx < 4 {
            ordered[idx] = Some(focus);
        }

        let parity = CipConfiguration::parity4(&ordered, carriers)?;
        if parity == 0 {
            return Err(CipLabelerError::CarrierMismatch);
        }

        let mut config = digraph.mol().atom(focus)?.chiral_tag();
        if parity == 1 {
            config = match config {
                ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                _ => config,
            };
        }

        if config == ChiralTag::TetrahedralCcw {
            if priority.is_pseudo_asymetric() {
                Ok(Descriptor::s)
            } else {
                Ok(Descriptor::S)
            }
        } else if config == ChiralTag::TetrahedralCw {
            if priority.is_pseudo_asymetric() {
                Ok(Descriptor::r)
            } else {
                Ok(Descriptor::R)
            }
        } else {
            Ok(Descriptor::Unknown)
        }
    }
}

pub(crate) struct CipSp2Bond<'a> {
    configuration: CipConfiguration<'a>,
    bond_idx: usize,
    cfg: BondStereo,
    ranked_anchors: Vec<u32>,
    primary_label: Option<CipBondPrimaryLabel>,
}

impl<'a> CipSp2Bond<'a> {
    // BEGIN RDKIT CPP FUNCTION Sp2Bond::Sp2Bond (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: Sp2Bond::Sp2Bond(const CIPMol &mol, Bond *bond, Atom *startAtom, Atom *endAtom,
    // RDKit✔️✔️:                  Bond::BondStereo cfg)
    // RDKit✔️✔️:     : Configuration(mol, {startAtom, endAtom}), dp_bond{bond}, d_cfg{cfg} {
    // RDKit✔️✔️:   CHECK_INVARIANT(startAtom && endAtom, "bad foci")
    // RDKit✔️✔️:   CHECK_INVARIANT(d_cfg == Bond::STEREOTRANS || d_cfg == Bond::STEREOCIS,
    // RDKit✔️✔️:                   "bad config")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto stereo_atoms = Chirality::findStereoAtoms(bond);
    // RDKit✔️✔️:   CHECK_INVARIANT(stereo_atoms.size() == 2, "incorrect number of stereo atoms")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Atom *> anchors{
    // RDKit✔️✔️:       {mol.getAtom(stereo_atoms[0]), mol.getAtom(stereo_atoms[1])}};
    // RDKit✔️✔️:
    // RDKit✔️✔️:   setCarriers(std::move(anchors));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::Sp2Bond
    // BEGIN RDKIT CPP FUNCTION Chirality::findStereoAtoms (GraphMol/Chirality.cpp)
    // RDKit✔️✔️: INT_VECT findStereoAtoms(const Bond *bond) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:   PRECONDITION(bond->hasOwningMol(), "no mol");
    // RDKit✔️✔️:   PRECONDITION(bond->getBondType() == Bond::DOUBLE, "not double bond");
    // RDKit✔️✔️:   PRECONDITION(bond->getStereo() > Bond::BondStereo::STEREOANY,
    // RDKit✔️✔️:                "no defined stereo");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!bond->getStereoAtoms().empty()) {
    // RDKit✔️✔️:     return bond->getStereoAtoms();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bond->getStereo() == Bond::BondStereo::STEREOE ||
    // RDKit✔️✔️:       bond->getStereo() == Bond::BondStereo::STEREOZ) {
    // RDKit✔️✔️:     const Atom *startStereoAtom =
    // RDKit✔️✔️:         findHighestCIPNeighbor(bond->getBeginAtom(), bond->getEndAtom());
    // RDKit✔️✔️:     const Atom *endStereoAtom =
    // RDKit✔️✔️:         findHighestCIPNeighbor(bond->getEndAtom(), bond->getBeginAtom());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (startStereoAtom == nullptr || endStereoAtom == nullptr) {
    // RDKit✔️✔️:       return {};
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     int startStereoAtomIdx = static_cast<int>(startStereoAtom->getIdx());
    // RDKit✔️✔️:     int endStereoAtomIdx = static_cast<int>(endStereoAtom->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     return {startStereoAtomIdx, endStereoAtomIdx};
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog) << "Unable to assign stereo atoms for bond "
    // RDKit✔️✔️:                             << bond->getIdx() << std::endl;
    // RDKit✔️✔️:     return {};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::findStereoAtoms
    pub(crate) fn new(
        molecule: &'a TopologyBlock,
        bond_idx: usize,
        start_atom: usize,
        end_atom: usize,
        cfg: BondStereo,
    ) -> Result<Self, CipLabelerError> {
        let mut configuration =
            CipConfiguration::with_foci(molecule, vec![start_atom, end_atom], false)?;
        configuration.digraph.mol().atom(start_atom)?;
        configuration.digraph.mol().atom(end_atom)?;
        let bond = configuration.digraph.mol().bond(bond_idx)?;
        if bond.order() != BondOrder::Double {
            return Err(CipLabelerError::BadSp2BondFoci);
        }
        if !((bond.begin().index() == start_atom && bond.end().index() == end_atom)
            || (bond.begin().index() == end_atom && bond.end().index() == start_atom))
        {
            return Err(CipLabelerError::BadSp2BondFoci);
        }
        if !matches!(cfg, BondStereo::Trans | BondStereo::Cis) {
            return Err(CipLabelerError::BadSp2BondConfig);
        }

        let mut ranks = vec![0_u32; molecule.atoms.len()];
        if bond.stereo_atoms().is_none() {
            for focus in [start_atom, end_atom] {
                let skip = if focus == start_atom {
                    end_atom
                } else {
                    start_atom
                };
                for neighbor in molecule.adjacency.neighbors_of(focus) {
                    if neighbor.atom_index == skip {
                        continue;
                    }
                    let rank = molecule.atoms[neighbor.atom_index]
                        .prop("_CIPRank")
                        .and_then(|value| value.parse::<u32>().ok())
                        .ok_or(CipLabelerError::IncorrectNumberOfStereoAtoms)?;
                    ranks[neighbor.atom_index] = rank;
                }
            }
        }
        let stereo_atoms = find_double_bond_stereo_atoms(molecule, BondId::new(bond_idx), &ranks)?
            .ok_or(CipLabelerError::IncorrectNumberOfStereoAtoms)?
            .map(AtomId::index);
        configuration.digraph.mol().atom(stereo_atoms[0])?;
        configuration.digraph.mol().atom(stereo_atoms[1])?;
        configuration.set_carriers(vec![Some(stereo_atoms[0]), Some(stereo_atoms[1])]);

        Ok(Self {
            configuration,
            bond_idx,
            cfg,
            ranked_anchors: Vec::new(),
            primary_label: None,
        })
    }

    pub(crate) fn get_foci(&self) -> &[usize] {
        self.configuration.get_foci()
    }

    pub(crate) fn get_carriers(&self) -> &[Option<usize>] {
        self.configuration.get_carriers()
    }

    pub(crate) fn ranked_anchors(&self) -> &[u32] {
        &self.ranked_anchors
    }

    pub(crate) fn primary_label(&self) -> Option<&CipBondPrimaryLabel> {
        self.primary_label.as_ref()
    }

    // BEGIN RDKIT CPP FUNCTION Sp2Bond::setPrimaryLabel (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: void Sp2Bond::setPrimaryLabel(Descriptor desc) {
    // RDKit✔️✔️:   switch (desc) {
    // RDKit✔️✔️:     case Descriptor::seqTrans:
    // RDKit✔️✔️:     case Descriptor::E:
    // RDKit✔️✔️:     case Descriptor::seqCis:
    // RDKit✔️✔️:     case Descriptor::Z: {
    // RDKit✔️✔️:       auto carriers = getCarriers();
    // RDKit✔️✔️:       dp_bond->setStereoAtoms(carriers[0]->getIdx(), carriers[1]->getIdx());
    // RDKit✔️✔️:       dp_bond->setStereo(d_cfg);
    // RDKit✔️✔️:       dp_bond->setProp(common_properties::_CIPCode, to_string(desc));
    // RDKit✔️✔️:       dp_bond->setProp(common_properties::_CIPNeighborOrder, d_ranked_anchors,
    // RDKit✔️✔️:                        true);
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     case Descriptor::R:
    // RDKit✔️✔️:     case Descriptor::S:
    // RDKit✔️✔️:     case Descriptor::r:
    // RDKit✔️✔️:     case Descriptor::s:
    // RDKit✔️✔️:     case Descriptor::M:
    // RDKit✔️✔️:     case Descriptor::P:
    // RDKit✔️✔️:     case Descriptor::m:
    // RDKit✔️✔️:     case Descriptor::p:
    // RDKit✔️✔️:     case Descriptor::SP_4:
    // RDKit✔️✔️:     case Descriptor::TBPY_5:
    // RDKit✔️✔️:     case Descriptor::OC_6:
    // RDKit✔️✔️:       throw std::runtime_error(
    // RDKit✔️✔️:           "Received a Descriptor that is not supported for double bonds");
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       throw std::runtime_error("Received an invalid Bond Descriptor");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::setPrimaryLabel
    pub(crate) fn set_primary_label(&mut self, desc: Descriptor) -> Result<(), CipLabelerError> {
        match desc {
            Descriptor::seqTrans | Descriptor::E | Descriptor::seqCis | Descriptor::Z => {
                let carriers = self.configuration.get_carriers();
                let stereo_atoms = [
                    carriers
                        .first()
                        .and_then(|carrier| *carrier)
                        .ok_or(CipLabelerError::IncorrectNumberOfStereoAtoms)?,
                    carriers
                        .get(1)
                        .and_then(|carrier| *carrier)
                        .ok_or(CipLabelerError::IncorrectNumberOfStereoAtoms)?,
                ];
                self.primary_label = Some(CipBondPrimaryLabel {
                    bond_idx: self.bond_idx,
                    stereo_atoms,
                    stereo: self.cfg,
                    cip_code: descriptor_to_string(desc),
                    cip_neighbor_order: self.ranked_anchors.clone(),
                });
                Ok(())
            }
            Descriptor::R
            | Descriptor::S
            | Descriptor::r
            | Descriptor::s
            | Descriptor::M
            | Descriptor::P
            | Descriptor::m
            | Descriptor::p
            | Descriptor::SP_4
            | Descriptor::TBPY_5
            | Descriptor::OC_6 => Err(CipLabelerError::DescriptorNotSupportedForDoubleBonds),
            Descriptor::None | Descriptor::Unknown | Descriptor::ns => {
                Err(CipLabelerError::InvalidBondDescriptor)
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION Sp2Bond::hasPrimaryLabel (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: bool Sp2Bond::hasPrimaryLabel() const {
    // RDKit✔️✔️:   return dp_bond->hasProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::hasPrimaryLabel
    pub(crate) fn has_primary_label(&self) -> bool {
        self.primary_label.is_some()
            || self
                .configuration
                .digraph
                .mol()
                .bond(self.bond_idx)
                .is_ok_and(|bond| bond.prop("_CIPCode").is_some())
    }

    // BEGIN RDKIT CPP FUNCTION Sp2Bond::resetPrimaryLabel (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: void Sp2Bond::resetPrimaryLabel() const {
    // RDKit✔️✔️:   dp_bond->clearProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::resetPrimaryLabel
    pub(crate) fn reset_primary_label(&mut self) {
        self.primary_label = None;
    }

    // BEGIN RDKIT CPP FUNCTION Sp2Bond::label(const Rules &) (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: Descriptor Sp2Bond::label(const Rules &comp) {
    // RDKit✔️✔️:   auto &digraph = getDigraph();
    // RDKit✔️✔️:   auto root1 = digraph.getOriginalRoot();
    // RDKit✔️✔️:   if (digraph.getCurrentRoot() != root1) {
    // RDKit✔️✔️:     digraph.changeRoot(root1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return label(root1, digraph, comp);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::label(const Rules &)
    pub(crate) fn label(
        &mut self,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let root1 = self.configuration.digraph.get_original_root();
        if self.configuration.digraph.get_current_root() != root1 {
            self.configuration.digraph.change_root(root1)?;
        }
        self.label_node(root1, comp, context)
    }

    pub(crate) fn label_with_external_digraph(
        &mut self,
        root1: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let foci = self.configuration.get_foci().to_vec();
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            &foci,
            &carriers,
            self.cfg,
            &mut self.ranked_anchors,
            root1,
            digraph,
            comp,
            context,
        )
    }

    // BEGIN RDKIT CPP FUNCTION Sp2Bond::label(Node *, Digraph &, const Rules &) (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: Descriptor Sp2Bond::label(Node *root1, Digraph &digraph, const Rules &comp) {
    // RDKit✔️✔️:   const auto &focus1 = getFoci()[0];
    // RDKit✔️✔️:   const auto &focus2 = getFoci()[1];
    // RDKit✔️✔️:
    // RDKit✔️✔️:   d_ranked_anchors.clear();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto &internal = findInternalEdge(root1->getEdges(), focus1, focus2);
    // RDKit✔️✔️:   if (internal == nullptr) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto &root2 = internal->getOther(root1);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto edges1 = root1->getEdges();
    // RDKit✔️✔️:   auto edges2 = root2->getEdges();
    // RDKit✔️✔️:   removeInternalEdges(edges1, focus1, focus2);
    // RDKit✔️✔️:   removeInternalEdges(edges2, focus1, focus2);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto carriers = std::vector<Atom *>(getCarriers());
    // RDKit✔️✔️:   auto config = d_cfg;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (root1->getAtom() == focus2) {
    // RDKit✔️✔️:     std::swap(carriers[0], carriers[1]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   digraph.changeRoot(root1);
    // RDKit✔️✔️:   const auto &priority1 = comp.sort(root1, edges1);
    // RDKit✔️✔️:   if (!priority1.isUnique()) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // swap
    // RDKit✔️✔️:   if (edges1.size() > 1 && carriers[0] != edges1[0]->getEnd()->getAtom()) {
    // RDKit✔️✔️:     if (config == Bond::STEREOCIS) {
    // RDKit✔️✔️:       config = Bond::STEREOTRANS;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       config = Bond::STEREOCIS;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   digraph.changeRoot(root2);
    // RDKit✔️✔️:   const auto &priority2 = comp.sort(root2, edges2);
    // RDKit✔️✔️:   if (!priority2.isUnique()) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // swap
    // RDKit✔️✔️:   if (edges2.size() > 1 && carriers[1] != edges2[0]->getEnd()->getAtom()) {
    // RDKit✔️✔️:     if (config == Bond::STEREOCIS) {
    // RDKit✔️✔️:       config = Bond::STEREOTRANS;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       config = Bond::STEREOCIS;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   {
    // RDKit✔️✔️:     auto carrier1_idx = edges1[0]->getEnd()->getAtomIdx();
    // RDKit✔️✔️:     auto carrier2_idx = edges2[0]->getEnd()->getAtomIdx();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (edges1[0]->getBeg()->getAtom() == focus1) {
    // RDKit✔️✔️:       d_ranked_anchors.assign({carrier1_idx, carrier2_idx});
    // RDKit✔️✔️:     } else if (edges2[0]->getBeg()->getAtom() == focus1) {
    // RDKit✔️✔️:       d_ranked_anchors.assign({carrier2_idx, carrier1_idx});
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (config == Bond::STEREOCIS) {
    // RDKit✔️✔️:     if (priority1.isPseudoAsymetric() != priority2.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::seqCis;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::Z;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (config == Bond::STEREOTRANS) {
    // RDKit✔️✔️:     if (priority1.isPseudoAsymetric() != priority2.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::seqTrans;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::E;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return Descriptor::UNKNOWN;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::label(Node *, Digraph &, const Rules &)
    fn label_node(
        &mut self,
        root1: CipNodeId,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let foci = self.configuration.get_foci().to_vec();
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            &foci,
            &carriers,
            self.cfg,
            &mut self.ranked_anchors,
            root1,
            &mut self.configuration.digraph,
            comp,
            context,
        )
    }

    fn label_node_impl(
        foci: &[usize],
        carriers: &[Option<usize>],
        cfg: BondStereo,
        ranked_anchors: &mut Vec<u32>,
        root1: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let focus1 = foci[0];
        let focus2 = foci[1];
        ranked_anchors.clear();

        let root1_edges = digraph.node_edges(root1)?;
        let Some(internal) =
            CipConfiguration::find_internal_edge(digraph, &root1_edges, focus1, focus2)
        else {
            return Ok(Descriptor::Unknown);
        };
        let root2 = digraph.edge(internal).get_other(internal, root1)?;

        let mut edges1 = digraph.node_edges(root1)?;
        let mut edges2 = digraph.node_edges(root2)?;
        CipConfiguration::remove_internal_edges(digraph, &mut edges1, focus1, focus2);
        CipConfiguration::remove_internal_edges(digraph, &mut edges2, focus1, focus2);

        if edges1.is_empty() || edges2.is_empty() {
            return Ok(Descriptor::Unknown);
        }

        let mut carriers = carriers.to_vec();
        let mut config = cfg;
        if digraph.node(root1).atom_idx() == Some(focus2) {
            carriers.swap(0, 1);
        }

        digraph.change_root(root1)?;
        let priority1 = comp.sort(digraph, context, root1, &mut edges1, true)?;
        if !priority1.is_unique() {
            return Ok(Descriptor::Unknown);
        }
        if edges1.len() > 1
            && carriers[0] != digraph.node(digraph.edge(edges1[0]).get_end()).atom_idx()
        {
            config = match config {
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                _ => config,
            };
        }

        digraph.change_root(root2)?;
        let priority2 = comp.sort(digraph, context, root2, &mut edges2, true)?;
        if !priority2.is_unique() {
            return Ok(Descriptor::Unknown);
        }
        if edges2.len() > 1
            && carriers[1] != digraph.node(digraph.edge(edges2[0]).get_end()).atom_idx()
        {
            config = match config {
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                _ => config,
            };
        }

        let carrier1_idx = digraph
            .node(digraph.edge(edges1[0]).get_end())
            .get_atom_idx()?;
        let carrier2_idx = digraph
            .node(digraph.edge(edges2[0]).get_end())
            .get_atom_idx()?;
        if digraph.node(digraph.edge(edges1[0]).get_beg()).atom_idx() == Some(focus1) {
            ranked_anchors.extend([carrier1_idx, carrier2_idx]);
        } else if digraph.node(digraph.edge(edges2[0]).get_beg()).atom_idx() == Some(focus1) {
            ranked_anchors.extend([carrier2_idx, carrier1_idx]);
        }

        if config == BondStereo::Cis {
            if priority1.is_pseudo_asymetric() != priority2.is_pseudo_asymetric() {
                Ok(Descriptor::seqCis)
            } else {
                Ok(Descriptor::Z)
            }
        } else if config == BondStereo::Trans {
            if priority1.is_pseudo_asymetric() != priority2.is_pseudo_asymetric() {
                Ok(Descriptor::seqTrans)
            } else {
                Ok(Descriptor::E)
            }
        } else {
            Ok(Descriptor::Unknown)
        }
    }
}

pub(crate) struct CipAtropisomerBond<'a> {
    configuration: CipConfiguration<'a>,
    bond_idx: usize,
    cfg: BondStereo,
    ranked_anchors: Vec<u32>,
    primary_label: Option<CipAtropisomerBondPrimaryLabel>,
}

impl<'a> CipAtropisomerBond<'a> {
    // BEGIN RDKIT CPP FUNCTION AtropisomerBond::AtropisomerBond (configs/AtropisomerBond.cpp)
    // RDKit✔️✔️: AtropisomerBond::AtropisomerBond(const CIPMol &mol, Bond *bond, Atom *startAtom,
    // RDKit✔️✔️:                                  Atom *endAtom, Bond::BondStereo cfg)
    // RDKit✔️✔️:     : Configuration(mol, {startAtom, endAtom}, true),
    // RDKit✔️✔️:       dp_bond{bond},
    // RDKit✔️✔️:       d_cfg{cfg} {
    // RDKit✔️✔️:   CHECK_INVARIANT(startAtom && endAtom, "bad foci")
    // RDKit✔️✔️:   CHECK_INVARIANT(d_cfg == Bond::STEREOATROPCW || d_cfg == Bond::STEREOATROPCCW,
    // RDKit✔️✔️:                   "bad config")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Atropisomers::AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit✔️✔️:   if (!Atropisomers::getAtropisomerAtomsAndBonds(bond, atomAndBondVecs,
    // RDKit✔️✔️:                                                  bond->getOwningMol())) {
    // RDKit✔️✔️:     return;  // not an atropisomer
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atom1 = mol.getAtom(atomAndBondVecs[0].second[0]->getOtherAtomIdx(
    // RDKit✔️✔️:       atomAndBondVecs[0].first->getIdx()));
    // RDKit✔️✔️:   auto atom2 = mol.getAtom(atomAndBondVecs[1].second[0]->getOtherAtomIdx(
    // RDKit✔️✔️:       atomAndBondVecs[1].first->getIdx()));
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Atom *> anchors{atom1, atom2};
    // RDKit✔️✔️:
    // RDKit✔️✔️:   setCarriers(std::move(anchors));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtropisomerBond::AtropisomerBond
    pub(crate) fn new(
        molecule: &'a TopologyBlock,
        bond_idx: usize,
        start_atom: usize,
        end_atom: usize,
        cfg: BondStereo,
    ) -> Result<Self, CipLabelerError> {
        let mut configuration =
            CipConfiguration::with_foci(molecule, vec![start_atom, end_atom], true)?;
        configuration.digraph.mol().atom(start_atom)?;
        configuration.digraph.mol().atom(end_atom)?;
        let bond = configuration.digraph.mol().bond(bond_idx)?;
        if !((bond.begin().index() == start_atom && bond.end().index() == end_atom)
            || (bond.begin().index() == end_atom && bond.end().index() == start_atom))
        {
            return Err(CipLabelerError::BadAtropisomerBondFoci);
        }
        if !matches!(cfg, BondStereo::AtropCw | BondStereo::AtropCcw) {
            return Err(CipLabelerError::BadAtropisomerBondConfig);
        }

        if let Some(ends) = atropisomer_carriers(molecule, BondId::new(bond_idx))? {
            let mut carriers = Vec::with_capacity(2);
            for end in ends {
                let carrier_bond = *end
                    .carrier_bonds()
                    .first()
                    .ok_or(CipLabelerError::BadAtropisomerBondFoci)?;
                let bond = &molecule.bonds[carrier_bond.index()];
                let carrier = if bond.begin() == end.focus() {
                    bond.end()
                } else if bond.end() == end.focus() {
                    bond.begin()
                } else {
                    return Err(CipLabelerError::BondNotIncident {
                        bond: carrier_bond.index(),
                        atom: end.focus().index(),
                    });
                };
                carriers.push(Some(carrier.index()));
            }
            configuration.set_carriers(carriers);
        }

        Ok(Self {
            configuration,
            bond_idx,
            cfg,
            ranked_anchors: Vec::new(),
            primary_label: None,
        })
    }

    pub(crate) fn get_foci(&self) -> &[usize] {
        self.configuration.get_foci()
    }

    pub(crate) fn get_carriers(&self) -> &[Option<usize>] {
        self.configuration.get_carriers()
    }

    pub(crate) fn ranked_anchors(&self) -> &[u32] {
        &self.ranked_anchors
    }

    pub(crate) fn primary_label(&self) -> Option<&CipAtropisomerBondPrimaryLabel> {
        self.primary_label.as_ref()
    }

    // BEGIN RDKIT CPP FUNCTION AtropisomerBond::setPrimaryLabel (configs/AtropisomerBond.cpp)
    // RDKit✔️✔️: void AtropisomerBond::setPrimaryLabel(Descriptor desc) {
    // RDKit✔️✔️:   switch (desc) {
    // RDKit✔️✔️:     case Descriptor::M:
    // RDKit✔️✔️:     case Descriptor::P:
    // RDKit✔️✔️:     case Descriptor::m:
    // RDKit✔️✔️:     case Descriptor::p: {
    // RDKit✔️✔️:       dp_bond->setProp(common_properties::_CIPCode, to_string(desc));
    // RDKit✔️✔️:       dp_bond->setProp(common_properties::_CIPNeighborOrder, d_ranked_anchors,
    // RDKit✔️✔️:                        true);
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     case Descriptor::R:
    // RDKit✔️✔️:     case Descriptor::S:
    // RDKit✔️✔️:     case Descriptor::r:
    // RDKit✔️✔️:     case Descriptor::s:
    // RDKit✔️✔️:     case Descriptor::SP_4:
    // RDKit✔️✔️:     case Descriptor::TBPY_5:
    // RDKit✔️✔️:     case Descriptor::OC_6:
    // RDKit✔️✔️:     case Descriptor::seqTrans:
    // RDKit✔️✔️:     case Descriptor::E:
    // RDKit✔️✔️:     case Descriptor::seqCis:
    // RDKit✔️✔️:     case Descriptor::Z:
    // RDKit✔️✔️:       throw std::runtime_error(
    // RDKit✔️✔️:           "Received a Descriptor that is not supported for atropisomer bonds");
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       throw std::runtime_error("Received an invalid Bond Descriptor");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtropisomerBond::setPrimaryLabel
    pub(crate) fn set_primary_label(&mut self, desc: Descriptor) -> Result<(), CipLabelerError> {
        match desc {
            Descriptor::M | Descriptor::P | Descriptor::m | Descriptor::p => {
                self.primary_label = Some(CipAtropisomerBondPrimaryLabel {
                    bond_idx: self.bond_idx,
                    cip_code: descriptor_to_string(desc),
                    cip_neighbor_order: self.ranked_anchors.clone(),
                });
                Ok(())
            }
            Descriptor::R
            | Descriptor::S
            | Descriptor::r
            | Descriptor::s
            | Descriptor::SP_4
            | Descriptor::TBPY_5
            | Descriptor::OC_6
            | Descriptor::seqTrans
            | Descriptor::E
            | Descriptor::seqCis
            | Descriptor::Z => Err(CipLabelerError::DescriptorNotSupportedForAtropisomerBonds),
            Descriptor::None | Descriptor::Unknown | Descriptor::ns => {
                Err(CipLabelerError::InvalidBondDescriptor)
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION AtropisomerBond::hasPrimaryLabel/resetPrimaryLabel (configs/AtropisomerBond.cpp)
    // RDKit✔️✔️: bool AtropisomerBond::hasPrimaryLabel() const {
    // RDKit✔️✔️:   return dp_bond->hasProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: void AtropisomerBond::resetPrimaryLabel() const {
    // RDKit✔️✔️:   dp_bond->clearProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtropisomerBond::hasPrimaryLabel/resetPrimaryLabel
    pub(crate) fn has_primary_label(&self) -> bool {
        self.primary_label.is_some()
            || self
                .configuration
                .digraph
                .mol()
                .bond(self.bond_idx)
                .is_ok_and(|bond| bond.prop("_CIPCode").is_some())
    }

    pub(crate) fn reset_primary_label(&mut self) {
        self.primary_label = None;
    }

    // BEGIN RDKIT CPP FUNCTION AtropisomerBond::label(const Rules &) (configs/AtropisomerBond.cpp)
    // RDKit✔️✔️: Descriptor AtropisomerBond::label(const Rules &comp) {
    // RDKit✔️✔️:   auto &digraph = getDigraph();
    // RDKit✔️✔️:   auto root1 = digraph.getOriginalRoot();
    // RDKit✔️✔️:   if (digraph.getCurrentRoot() != root1) {
    // RDKit✔️✔️:     digraph.changeRoot(root1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return label(root1, digraph, comp);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtropisomerBond::label(const Rules &)
    pub(crate) fn label(
        &mut self,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let root1 = self.configuration.digraph.get_original_root();
        if self.configuration.digraph.get_current_root() != root1 {
            self.configuration.digraph.change_root(root1)?;
        }
        self.label_with_root(root1, comp, context)
    }

    pub(crate) fn label_with_external_digraph(
        &mut self,
        root1: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let foci = self.configuration.get_foci().to_vec();
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            &foci,
            &carriers,
            self.cfg,
            &mut self.ranked_anchors,
            root1,
            digraph,
            comp,
            context,
        )
    }

    // BEGIN RDKIT CPP FUNCTION AtropisomerBond::label(Node *, Digraph &, const Rules &) (configs/AtropisomerBond.cpp)
    // RDKit✔️✔️: Descriptor AtropisomerBond::label(Node *root1, Digraph &digraph,
    // RDKit✔️✔️:                                   const Rules &comp) {
    // RDKit✔️✔️:   const auto &focus1 = getFoci()[0];
    // RDKit✔️✔️:   const auto &focus2 = getFoci()[1];
    // RDKit✔️✔️:
    // RDKit✔️✔️:   d_ranked_anchors.clear();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto &internal = findInternalEdge(root1->getEdges(), focus1, focus2);
    // RDKit✔️✔️:   if (internal == nullptr) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto &root2 = internal->getOther(root1);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto edges1 = root1->getEdges();
    // RDKit✔️✔️:   auto edges2 = root2->getEdges();
    // RDKit✔️✔️:   removeInternalEdges(edges1, focus1, focus2);
    // RDKit✔️✔️:   removeInternalEdges(edges2, focus1, focus2);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   removeDuplicatesAndHs(edges1);
    // RDKit✔️✔️:   removeDuplicatesAndHs(edges2);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto carriers = std::vector<Atom *>(getCarriers());
    // RDKit✔️✔️:   auto config = d_cfg;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (root1->getAtom() == focus2) {
    // RDKit✔️✔️:     std::swap(carriers[0], carriers[1]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   digraph.changeRoot(root1);
    // RDKit✔️✔️:   const auto &priority1 = comp.sort(root1, edges1);
    // RDKit✔️✔️:   if (!priority1.isUnique()) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // swap
    // RDKit✔️✔️:   if (edges1.size() > 1 && carriers[0] == edges1[1]->getEnd()->getAtom()) {
    // RDKit✔️✔️:     if (config == Bond::STEREOATROPCCW) {
    // RDKit✔️✔️:       config = Bond::STEREOATROPCW;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       config = Bond::STEREOATROPCCW;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   digraph.changeRoot(root2);
    // RDKit✔️✔️:   const auto &priority2 = comp.sort(root2, edges2);
    // RDKit✔️✔️:   if (!priority2.isUnique()) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // swap
    // RDKit✔️✔️:   if (edges2.size() > 1 && carriers[1] == edges2[1]->getEnd()->getAtom()) {
    // RDKit✔️✔️:     if (config == Bond::STEREOATROPCCW) {
    // RDKit✔️✔️:       config = Bond::STEREOATROPCW;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       config = Bond::STEREOATROPCCW;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (config == Bond::STEREOATROPCCW) {
    // RDKit✔️✔️:     if (priority1.isPseudoAsymetric() || priority2.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::m;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::M;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (config == Bond::STEREOATROPCW) {
    // RDKit✔️✔️:     if (priority1.isPseudoAsymetric() || priority2.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::p;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::P;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return Descriptor::UNKNOWN;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtropisomerBond::label(Node *, Digraph &, const Rules &)
    fn label_with_root(
        &mut self,
        root1: CipNodeId,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let foci = self.configuration.get_foci().to_vec();
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            &foci,
            &carriers,
            self.cfg,
            &mut self.ranked_anchors,
            root1,
            &mut self.configuration.digraph,
            comp,
            context,
        )
    }

    fn label_node_impl(
        foci: &[usize],
        carriers: &[Option<usize>],
        cfg: BondStereo,
        ranked_anchors: &mut Vec<u32>,
        root1: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let focus1 = foci[0];
        let focus2 = foci[1];
        ranked_anchors.clear();
        if carriers.len() < 2 {
            return Ok(Descriptor::Unknown);
        }

        let root1_edges = digraph.node_edges(root1)?;
        let Some(internal) =
            CipConfiguration::find_internal_edge(digraph, &root1_edges, focus1, focus2)
        else {
            return Ok(Descriptor::Unknown);
        };
        let root2 = digraph.edge(internal).get_other(internal, root1)?;

        let mut edges1 = digraph.node_edges(root1)?;
        let mut edges2 = digraph.node_edges(root2)?;
        CipConfiguration::remove_internal_edges(digraph, &mut edges1, focus1, focus2);
        CipConfiguration::remove_internal_edges(digraph, &mut edges2, focus1, focus2);
        CipConfiguration::remove_duplicates_and_hs(digraph, &mut edges1);
        CipConfiguration::remove_duplicates_and_hs(digraph, &mut edges2);
        if edges1.is_empty() || edges2.is_empty() {
            return Ok(Descriptor::Unknown);
        }

        let mut carriers = carriers.to_vec();
        let mut config = cfg;
        if digraph.node(root1).atom_idx() == Some(focus2) {
            carriers.swap(0, 1);
        }

        digraph.change_root(root1)?;
        let priority1 = comp.sort(digraph, context, root1, &mut edges1, true)?;
        if !priority1.is_unique() {
            return Ok(Descriptor::Unknown);
        }
        if edges1.len() > 1
            && carriers[0] == digraph.node(digraph.edge(edges1[1]).get_end()).atom_idx()
        {
            config = match config {
                BondStereo::AtropCcw => BondStereo::AtropCw,
                BondStereo::AtropCw => BondStereo::AtropCcw,
                _ => config,
            };
        }

        digraph.change_root(root2)?;
        let priority2 = comp.sort(digraph, context, root2, &mut edges2, true)?;
        if !priority2.is_unique() {
            return Ok(Descriptor::Unknown);
        }
        if edges2.len() > 1
            && carriers[1] == digraph.node(digraph.edge(edges2[1]).get_end()).atom_idx()
        {
            config = match config {
                BondStereo::AtropCcw => BondStereo::AtropCw,
                BondStereo::AtropCw => BondStereo::AtropCcw,
                _ => config,
            };
        }

        let carrier1_idx = digraph
            .node(digraph.edge(edges1[0]).get_end())
            .get_atom_idx()?;
        let carrier2_idx = digraph
            .node(digraph.edge(edges2[0]).get_end())
            .get_atom_idx()?;
        if digraph.node(digraph.edge(edges1[0]).get_beg()).atom_idx() == Some(focus1) {
            ranked_anchors.extend([carrier1_idx, carrier2_idx]);
        } else if digraph.node(digraph.edge(edges2[0]).get_beg()).atom_idx() == Some(focus1) {
            ranked_anchors.extend([carrier2_idx, carrier1_idx]);
        }

        if config == BondStereo::AtropCcw {
            if priority1.is_pseudo_asymetric() || priority2.is_pseudo_asymetric() {
                Ok(Descriptor::m)
            } else {
                Ok(Descriptor::M)
            }
        } else if config == BondStereo::AtropCw {
            if priority1.is_pseudo_asymetric() || priority2.is_pseudo_asymetric() {
                Ok(Descriptor::p)
            } else {
                Ok(Descriptor::P)
            }
        } else {
            Ok(Descriptor::Unknown)
        }
    }
}
