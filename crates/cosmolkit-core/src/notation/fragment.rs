// RDKit marker convention defined in dev/source_reproduction_protocol.md.
//
// Port of RDKit MolOps fragmentation functions from:
//   third_party/rdkit/Code/GraphMol/MolOps.cpp
//
// The C++ source uses boost::connected_components for graph connectivity.
// Here we implement an equivalent DFS-based connected-components algorithm.

use crate::{
    AdjacencyList, Atom, AtomSpec, BondSpec, Molecule, MoleculeBuilder, error::MoleculeBuildError,
};

// ---------------------------------------------------------------------------
// FragmentError
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum FragmentError {
    #[error("molecule has no atoms")]
    EmptyMolecule,
    #[error("failed to build fragment molecule: {0}")]
    BuildError(#[from] MoleculeBuildError),
}

// ---------------------------------------------------------------------------
// find_connected_components — DFS-based graph traversal
//
// RDKit equivalent: boost::connected_components (used in getMolFrags)
// RDKit source: MolOps.cpp lines 852-857
// RDKit✔️✔️: unsigned int getMolFrags(const ROMol &mol, INT_VECT &mapping) {
// RDKit✔️✔️:   unsigned int natms = mol.getNumAtoms();
// RDKit✔️✔️:   mapping.resize(natms);
// RDKit✔️✔️:   return natms ? boost::connected_components(mol.getTopology(), &mapping[0])
// RDKit✔️✔️:                : 0;
// RDKit✔️✔️: };
//
// Boost connected_components assigns a component ID to each vertex.
// We re-implement with iterative DFS over CSR adjacency.
// ---------------------------------------------------------------------------
fn find_connected_components(mol: &Molecule) -> Vec<usize> {
    let num_atoms = mol.num_atoms();
    if num_atoms == 0 {
        return Vec::new();
    }

    // Build adjacency from bonds (like boost does from the graph topology)
    let adjacency = AdjacencyList::from_topology(num_atoms, mol.bonds());

    // RDKit (line 854): mapping.resize(natms);
    let mut component_of = vec![usize::MAX; num_atoms];
    let mut next_component: usize = 0;

    for start in 0..num_atoms {
        if component_of[start] != usize::MAX {
            continue;
        }

        // Iterative DFS (stack-based) — equivalent to boost's BFS traversal
        let mut stack: Vec<usize> = Vec::new();
        stack.push(start);
        component_of[start] = next_component;

        while let Some(atom_idx) = stack.pop() {
            for neighbor_ref in adjacency.neighbors_of(atom_idx) {
                let nbr = neighbor_ref.atom_index;
                if component_of[nbr] == usize::MAX {
                    component_of[nbr] = next_component;
                    stack.push(nbr);
                }
            }
        }

        next_component += 1;
    }

    component_of
}

// ---------------------------------------------------------------------------
// get_fragment_atom_mapping
//
// RDKit equivalent: getMolFrags(mol, mapping) returning frag assignment per atom
// RDKit source: MolOps.cpp lines 852-857
//
// Returns a Vec where index i contains the fragment index for atom i.
// ---------------------------------------------------------------------------
pub fn get_fragment_atom_mapping(mol: &Molecule) -> Vec<usize> {
    find_connected_components(mol)
}

// ---------------------------------------------------------------------------
// get_num_fragments
//
// Counts the number of connected components.
// RDKit equivalent: getMolFrags returning unsigned int count
// RDKit source: MolOps.cpp lines 852-857
// ---------------------------------------------------------------------------
pub fn get_num_fragments(mol: &Molecule) -> usize {
    let component_of = find_connected_components(mol);
    if component_of.is_empty() {
        return 0;
    }
    // The last component ID is one less than the count of unique IDs,
    // because IDs are assigned consecutively starting from 0.
    *component_of.iter().max().unwrap_or(&0) + 1
}

// ---------------------------------------------------------------------------
// group_atoms_by_fragment
//
// RDKit equivalent: getMolFrags(mol, VECT_INT_VECT &frags) — lines 859-878
// RDKit✔️✔️: unsigned int getMolFrags(const ROMol &mol, VECT_INT_VECT &frags) {
// RDKit✔️✔️:   frags.clear();
// RDKit✔️✔️:   INT_VECT mapping;
// RDKit✔️✔️:   getMolFrags(mol, mapping);
// RDKit✔️✔️:   INT_INT_VECT_MAP comMap;
// RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
// RDKit✔️✔️:     int mi = mapping[i];
// RDKit✔️✔️:     if (comMap.find(mi) == comMap.end()) {
// RDKit✔️✔️:       INT_VECT comp;
// RDKit✔️✔️:       comMap[mi] = comp;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     comMap[mi].push_back(i);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   for (INT_INT_VECT_MAP_CI mci = comMap.begin(); mci != comMap.end(); mci++) {
// RDKit✔️✔️:     frags.push_back((*mci).second);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return rdcast<unsigned int>(frags.size());
// RDKit✔️✔️: }
// ---------------------------------------------------------------------------
fn group_atoms_by_fragment(mol: &Molecule) -> Vec<Vec<usize>> {
    let component_of = find_connected_components(mol);
    if component_of.is_empty() {
        return Vec::new();
    }

    let num_fragments = get_num_fragments(mol);
    let mut frags: Vec<Vec<usize>> = vec![Vec::new(); num_fragments];
    for (atom_idx, frag_idx) in component_of.into_iter().enumerate() {
        frags[frag_idx].push(atom_idx);
    }
    frags
}

// ---------------------------------------------------------------------------
// Helper: build a Molecule for one fragment from atom indices
//
// RDKit equivalent: the per-fragment molecule construction inside getTheFrags()
// RDKit source: MolOps.cpp lines 725-814
//
// We use MoleculeBuilder to construct a new molecule from the atoms and bonds
// that belong to this fragment.
// ---------------------------------------------------------------------------
pub(crate) fn build_fragment_molecule(
    mol: &Molecule,
    fragment_atoms: &[usize],
    copy_conformers: bool,
) -> Result<Molecule, FragmentError> {
    let mut builder = MoleculeBuilder::new();

    let num_frag_atoms = fragment_atoms.len();

    // Build old-index -> new-index mapping for atoms in this fragment
    // RDKit equivalent: the atomsInFrag bitset + new atom index tracking
    // RDKit source: MolOps.cpp lines 727-733
    let mut old_to_new = vec![usize::MAX; mol.num_atoms()];
    let mut atom_mapping = vec![None; mol.num_atoms()];
    for (new_idx, old_idx) in fragment_atoms.iter().enumerate() {
        old_to_new[*old_idx] = new_idx;
        atom_mapping[*old_idx] = Some(AtomId_usize(new_idx));
    }
    let mut bond_mapping = vec![None; mol.num_bonds()];

    // Add each atom to the builder
    // RDKit equivalent: getTheFrags — atoms are either copied via copyMolSubset
    // or implicitly kept when batch-removing non-fragment atoms (lines 801-809).
    // Here we add atoms explicitly.
    for old_idx in fragment_atoms {
        let atom = mol
            .atoms()
            .get(*old_idx)
            .expect("fragment atom index out of range");
        let spec = atom_to_spec(atom);
        builder.add_atom(spec);
    }

    // Add bonds whose both endpoints are in this fragment
    // RDKit equivalent: getTheFrags — after removing atoms, bonds between
    // surviving atoms are kept automatically (lines 804-809).
    for bond in mol.bonds() {
        let begin_old = bond.begin().index();
        let end_old = bond.end().index();
        if old_to_new[begin_old] != usize::MAX && old_to_new[end_old] != usize::MAX {
            let begin_new = AtomId_usize(old_to_new[begin_old]);
            let end_new = AtomId_usize(old_to_new[end_old]);

            // Remap stereo atoms if present
            let stereo_atoms = bond.stereo_atoms().map(|[sa, sb]| {
                (
                    AtomId_usize(old_to_new[sa.index()]),
                    AtomId_usize(old_to_new[sb.index()]),
                )
            });

            let mut spec = BondSpec::new(begin_new, end_new, bond.order())
                .with_aromatic(bond.is_aromatic())
                .with_conjugated(bond.is_conjugated())
                .with_direction(bond.direction())
                .with_stereo(bond.stereo());
            if let Some((sa, sb)) = stereo_atoms {
                spec = spec.with_stereo_atoms(sa, sb);
            }
            // Copy bond properties
            for (key, value) in bond.props() {
                spec = spec.with_prop(key.clone(), value.clone());
            }
            let new_bond = builder.add_bond(spec)?;
            bond_mapping[bond.id().index()] = Some(new_bond);
        }
    }

    for stereo_group in mol.stereo_groups() {
        if let Some(fragment_group) = stereo_group.remapped(&atom_mapping, &bond_mapping) {
            builder.add_stereo_group(fragment_group)?;
        }
    }

    // Copy conformers if requested
    // RDKit equivalent: getTheFrags lines 817-821
    if copy_conformers {
        if let Some(coords_2d) = mol.coordinates_2d() {
            let frag_coords: Vec<[f64; 2]> = fragment_atoms
                .iter()
                .filter_map(|old_idx| coords_2d.get(*old_idx).copied())
                .collect();
            if frag_coords.len() == num_frag_atoms {
                builder.set_2d_coordinates(frag_coords)?;
            }
        }
        for conformer in mol.conformers_3d() {
            let frag_coords: Vec<[f64; 3]> = fragment_atoms
                .iter()
                .filter_map(|old_idx| conformer.coordinates().get(*old_idx).copied())
                .collect();
            if frag_coords.len() == num_frag_atoms {
                builder.add_3d_conformer(frag_coords)?;
            }
        }
    }

    Ok(builder.build()?)
}

// ---------------------------------------------------------------------------
// Helper: convert Atom to AtomSpec for building a new molecule
// ---------------------------------------------------------------------------
fn atom_to_spec(atom: &Atom) -> AtomSpec {
    let mut spec = AtomSpec::new(atom.element())
        .with_formal_charge(atom.formal_charge())
        .with_explicit_hydrogens(atom.explicit_hydrogens())
        .with_chiral_tag(atom.chiral_tag())
        .with_aromatic(atom.is_aromatic())
        .with_no_implicit(atom.no_implicit())
        .with_radical_electrons(atom.radical_electrons())
        .with_hybridization(atom.hybridization());
    if let Some(isotope) = atom.isotope() {
        spec = spec.with_isotope(isotope);
    }
    if let Some(atom_map) = atom.atom_map() {
        spec = spec.with_atom_map(atom_map);
    }
    if let Some(chiral_perm) = atom.chiral_permutation() {
        spec = spec.with_chiral_permutation(chiral_perm);
    }
    if atom.unknown_stereo() {
        spec = spec.with_unknown_stereo(true);
    }
    if let Some(mol_parity) = atom.mol_parity() {
        spec = spec.with_mol_parity(mol_parity);
    }
    if let Some(mol_inv) = atom.mol_inversion_flag() {
        spec = spec.with_mol_inversion_flag(mol_inv);
    }
    if atom.implicit_hydrogen() {
        spec = spec.with_implicit_hydrogen(true);
    }
    let tracked: Vec<u16> = atom.tracked_isotopic_hydrogens().to_vec();
    if !tracked.is_empty() {
        spec = spec.with_tracked_isotopic_hydrogens(tracked);
    }
    if let Some(query) = atom.query() {
        spec = spec.with_query(query.clone());
    }
    for (key, value) in atom.props() {
        spec = spec.with_prop(key.clone(), value.clone());
    }
    if let Some(pdb_info) = atom.pdb_residue_info() {
        spec = spec.with_pdb_residue_info(pdb_info.clone());
    }
    spec
}

// ---------------------------------------------------------------------------
// Helper: wrap a usize into AtomId for builder calls
// ---------------------------------------------------------------------------
fn AtomId_usize(index: usize) -> crate::AtomId {
    crate::AtomId::new(index)
}

// ---------------------------------------------------------------------------
// get_mol_frags
//
// Decompose a molecule into its connected components (fragments).
// Returns a list of Molecule, one per fragment.
//
// RDKit equivalent: the molecule-producing overload of getMolFrags
// RDKit source: MolOps.cpp lines 704-835 (getTheFrags) + lines 838-888 (public overloads)
//
// RDKit✔️✔️: std::vector<std::unique_ptr<ROMol>> getTheFrags(
// RDKit✔️✔️:     const ROMol &mol, bool sanitizeFrags, INT_VECT *frags,
// RDKit✔️✔️:     VECT_INT_VECT *fragsMolAtomMapping, bool copyConformers) {
// RDKit✔️✔️:   std::unique_ptr<INT_VECT> mappingStorage;
// RDKit✔️✔️:   if (!frags) {
// RDKit✔️✔️:     mappingStorage.reset(new INT_VECT);
// RDKit✔️✔️:     frags = mappingStorage.get();
// RDKit✔️✔️:   }
// RDKit✔️✔️:   int nFrags = getMolFrags(mol, *frags);
// RDKit✔️✔️:   std::vector<std::unique_ptr<RWMol>> res;
// RDKit✔️✔️:   if (nFrags == 1) {
// RDKit✔️✔️:     res.emplace_back(new RWMol(mol));
// RDKit✔️✔️:     if (fragsMolAtomMapping) {
// RDKit✔️✔️:       INT_VECT comp;
// RDKit✔️✔️:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
// RDKit✔️✔️:         comp.push_back(idx);
// RDKit✔️✔️:       }
// RDKit✔️✔️:       (*fragsMolAtomMapping).push_back(comp);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     res.reserve(nFrags);
// RDKit✔️✔️:     for (int i = 0; i < nFrags; ++i) {
// RDKit✔️✔️:       boost::dynamic_bitset<> atomsInFrag(mol.getNumAtoms());
// RDKit✔️✔️:       INT_VECT comp;
// RDKit✔️✔️:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
// RDKit✔️✔️:         if ((*frags)[idx] == i) {
// RDKit✔️✔️:           comp.push_back(idx);
// RDKit✔️✔️:           atomsInFrag.set(idx);
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:       auto fragmentHasChallengingFeatures =
// RDKit✔️✔️:           [&](const INT_VECT &comp,
// RDKit✔️✔️:               const boost::dynamic_bitset<> &atomsInFrag) -> bool {
// RDKit✔️✔️:         for (auto idx : comp) {
// RDKit✔️✔️:           const auto atom = mol.getAtomWithIdx(idx);
// RDKit✔️✔️:           if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED &&
// RDKit✔️✔️:               atom->getChiralTag() != Atom::ChiralType::CHI_OTHER) {
// RDKit✔️✔️:             return true;
// RDKit✔️✔️:           }
// RDKit✔️✔️:           for (auto bnd : mol.atomBonds(atom)) {
// RDKit✔️✔️:             if (atomsInFrag[bnd->getOtherAtomIdx(idx)]) {
// RDKit✔️✔️:               if (bnd->getStereo() != Bond::BondStereo::STEREONONE &&
// RDKit✔️✔️:                   bnd->getStereo() != Bond::BondStereo::STEREOANY) {
// RDKit✔️✔️:                 return true;
// RDKit✔️✔️:               }
// RDKit✔️✔️:             }
// RDKit✔️✔️:           }
// RDKit✔️✔️:         }
// RDKit✔️✔️:         for (auto sgroup : getSubstanceGroups(mol)) {
// RDKit✔️✔️:           for (auto aid : sgroup.getAtoms()) {
// RDKit✔️✔️:             if (atomsInFrag[aid]) { return true; }
// RDKit✔️✔️:           }
// RDKit✔️✔️:           for (auto aid : sgroup.getParentAtoms()) {
// RDKit✔️✔️:             if (atomsInFrag[aid]) { return true; }
// RDKit✔️✔️:           }
// RDKit✔️✔️:         }
// RDKit✔️✔️:         for (auto stereoGroup : mol.getStereoGroups()) {
// RDKit✔️✔️:           for (auto atom : stereoGroup.getAtoms()) {
// RDKit✔️✔️:             if (atomsInFrag[atom->getIdx()]) { return true; }
// RDKit✔️✔️:           }
// RDKit✔️✔️:           for (auto bond : stereoGroup.getBonds()) {
// RDKit✔️✔️:             if (atomsInFrag[bond->getBeginAtomIdx()] &&
// RDKit✔️✔️:                 atomsInFrag[bond->getEndAtomIdx()]) { return true; }
// RDKit✔️✔️:           }
// RDKit✔️✔️:         }
// RDKit✔️✔️:         return false;
// RDKit✔️✔️:       };
// RDKit✔️✔️:       if (comp.size() == 1 || (nFrags > 3 && !fragmentHasChallengingFeatures(comp, atomsInFrag))) {
// RDKit✔️✔️:         SubsetOptions opts{.sanitize = sanitizeFrags,
// RDKit✔️✔️:                            .clearComputedProps = true,
// RDKit✔️✔️:                            .copyCoordinates = copyConformers,
// RDKit✔️✔️:                            .method = SubsetMethod::BONDS_BETWEEN_ATOMS};
// RDKit✔️✔️:         std::vector<unsigned int> atoms{comp.begin(), comp.end()};
// RDKit✔️✔️:         auto submol = copyMolSubset(mol, atoms, info, opts);
// RDKit✔️✔️:         res.push_back(std::move(submol));
// RDKit✔️✔️:       } else {
// RDKit✔️✔️:         res.emplace_back(new RWMol(mol));
// RDKit✔️✔️:         frag->beginBatchEdit();
// RDKit✔️✔️:         for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
// RDKit✔️✔️:           if (!atomsInFrag[idx]) { frag->removeAtom(idx); }
// RDKit✔️✔️:         }
// RDKit✔️✔️:         frag->commitBatchEdit();
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (fragsMolAtomMapping) { (*fragsMolAtomMapping).push_back(comp); }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (!copyConformers) {
// RDKit✔️✔️:     for (auto &frag : res) { frag->clearConformers(); }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (sanitizeFrags) {
// RDKit✔️✔️:     for (auto &frag : res) { sanitizeMol(*frag); }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   std::vector<std::unique_ptr<ROMol>> finalRes;
// RDKit✔️✔️:   for (auto &r : res) {
// RDKit✔️✔️:     finalRes.emplace_back(r.get());
// RDKit✔️✔️:     r.release();
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return finalRes;
// RDKit✔️✔️: }
// ---------------------------------------------------------------------------
pub fn get_mol_frags(mol: &Molecule) -> Result<Vec<Molecule>, FragmentError> {
    if mol.num_atoms() == 0 {
        return Err(FragmentError::EmptyMolecule);
    }

    let frag_groups = group_atoms_by_fragment(mol);

    if frag_groups.len() == 1 {
        // RDKit (line 716): res.emplace_back(new RWMol(mol));
        // Single fragment — return a copy of the input molecule
        return Ok(vec![mol.clone()]);
    }

    let mut result = Vec::with_capacity(frag_groups.len());
    for group in &frag_groups {
        let fragment = build_fragment_molecule(mol, group, true)?;
        result.push(fragment);
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// get_largest_fragment
//
// Returns the fragment with the most atoms.
// RDKit does not have a direct equivalent; this is a convenience wrapper
// over get_mol_frags.
// ---------------------------------------------------------------------------
pub fn get_largest_fragment(mol: &Molecule) -> Result<Molecule, FragmentError> {
    let fragments = get_mol_frags(mol)?;
    fragments
        .into_iter()
        .max_by_key(|frag| frag.num_atoms())
        .ok_or(FragmentError::EmptyMolecule)
}
