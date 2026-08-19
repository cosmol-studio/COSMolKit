// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::{
    cmp::Ordering,
    collections::{BTreeSet, VecDeque},
};

use crate::{
    AdjacencyList, Atom, AtomId, Bond, BondId, BondOrder, BondStereo, ChiralTag, KekulizeError,
    Molecule, RingInfo, StereoGroupKind, ValenceAssignment, ValenceModel, molecule::TopologyBlock,
    read_parts::MoleculeReadParts, stereo::StereoGroup,
};

struct CanonRankReadView<'a> {
    atoms: &'a [Atom],
    bonds: &'a [Bond],
    adjacency: &'a AdjacencyList,
    stereo_groups: &'a [StereoGroup],
    rings: RingInfo,
    valence: ValenceAssignment,
}

impl<'a> CanonRankReadView<'a> {
    fn from_molecule(molecule: &'a Molecule) -> Result<Self, KekulizeError> {
        let rings = crate::fast_find_rings(molecule)?;
        let valence =
            crate::valence::assign_valence_with_options(molecule, ValenceModel::RdkitLike, false)?;
        Ok(Self {
            atoms: molecule.atoms(),
            bonds: molecule.bonds(),
            adjacency: &molecule.topology_block().adjacency,
            stereo_groups: molecule.stereo_groups(),
            rings,
            valence,
        })
    }

    fn from_read_parts(read: MoleculeReadParts<'a>) -> Result<Self, KekulizeError> {
        Self::from_parts(
            read.atoms(),
            read.bonds(),
            read.adjacency(),
            read.stereo_groups(),
        )
    }

    fn from_topology(
        topology: &'a TopologyBlock,
        stereo_groups: &'a [StereoGroup],
    ) -> Result<Self, KekulizeError> {
        Self::from_parts(
            &topology.atoms,
            &topology.bonds,
            &topology.adjacency,
            stereo_groups,
        )
    }

    fn from_parts(
        atoms: &'a [Atom],
        bonds: &'a [Bond],
        adjacency: &'a AdjacencyList,
        stereo_groups: &'a [StereoGroup],
    ) -> Result<Self, KekulizeError> {
        // RDKit✔️✔️:   bool clearRings = false;
        // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
        // RDKit✔️✔️:     MolOps::fastFindRings(mol);
        // RDKit✔️✔️:     clearRings = true;
        // RDKit✔️✔️:   }
        let rings = crate::rings::fast_find_rings_from_parts(atoms.len(), bonds, adjacency)?;
        let valence = crate::valence::assign_valence_with_options_from_parts(
            atoms,
            bonds,
            adjacency,
            ValenceModel::RdkitLike,
            false,
        )?;
        Ok(Self {
            atoms,
            bonds,
            adjacency,
            stereo_groups,
            rings,
            valence,
        })
    }

    fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    fn atom_degree(&self, atom: AtomId) -> usize {
        self.adjacency.neighbors_of(atom.index()).len()
    }

    fn atom_neighbors(&self, atom: AtomId) -> Vec<usize> {
        self.adjacency
            .neighbors_of(atom.index())
            .iter()
            .map(|neighbor| neighbor.atom_index)
            .collect()
    }

    fn bond_other_atom_index(&self, bond_id: BondId, atom_id: AtomId) -> Option<usize> {
        let bond = self.bonds.get(bond_id.index())?;
        if bond.begin() == atom_id {
            Some(bond.end().index())
        } else if bond.end() == atom_id {
            Some(bond.begin().index())
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct FragmentRankScope<'a> {
    pub(crate) atoms_to_use: &'a [bool],
    pub(crate) bonds_to_use: &'a [bool],
    pub(crate) atom_symbols: Option<&'a [String]>,
    pub(crate) bond_symbols: Option<&'a [String]>,
}

impl<'a> FragmentRankScope<'a> {
    pub(crate) const fn new(atoms_to_use: &'a [bool], bonds_to_use: &'a [bool]) -> Self {
        Self {
            atoms_to_use,
            bonds_to_use,
            atom_symbols: None,
            bond_symbols: None,
        }
    }

    #[allow(dead_code)]
    pub(crate) const fn with_atom_symbols(mut self, atom_symbols: &'a [String]) -> Self {
        self.atom_symbols = Some(atom_symbols);
        self
    }

    #[allow(dead_code)]
    pub(crate) const fn with_bond_symbols(mut self, bond_symbols: &'a [String]) -> Self {
        self.bond_symbols = Some(bond_symbols);
        self
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CanonicalRankOptions {
    pub(crate) break_ties: bool,
    pub(crate) include_chirality: bool,
    pub(crate) include_isotopes: bool,
    pub(crate) include_atom_maps: bool,
    pub(crate) include_chiral_presence: bool,
    pub(crate) include_stereo_groups: bool,
    pub(crate) use_non_stereo_ranks: bool,
    pub(crate) include_ring_stereo: bool,
    pub(crate) chirality_rings_use_ring_stereo: bool,
}

impl CanonicalRankOptions {
    pub(crate) const fn rank_mol_default() -> Self {
        Self {
            break_ties: true,
            include_chirality: true,
            include_isotopes: true,
            include_atom_maps: true,
            include_chiral_presence: false,
            include_stereo_groups: true,
            use_non_stereo_ranks: false,
            include_ring_stereo: true,
            chirality_rings_use_ring_stereo: true,
        }
    }

    pub(crate) const fn kekulize_default() -> Self {
        Self {
            break_ties: true,
            include_chirality: true,
            include_isotopes: true,
            include_atom_maps: true,
            include_chiral_presence: false,
            include_stereo_groups: false,
            use_non_stereo_ranks: false,
            include_ring_stereo: true,
            chirality_rings_use_ring_stereo: false,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct CanonRankFlags {
    use_isotopes: bool,
    use_chirality: bool,
    use_chirality_rings: bool,
    use_atom_maps: bool,
    use_non_stereo_ranks: bool,
    use_chiral_presence: bool,
    use_atom_maps_on_dummies: bool,
}

impl CanonRankFlags {
    const fn from_fragment_options(options: CanonicalRankOptions) -> Self {
        Self {
            use_isotopes: options.include_isotopes,
            use_chirality: options.include_chirality,
            use_chirality_rings: if options.chirality_rings_use_ring_stereo {
                options.include_chirality && options.include_ring_stereo
            } else {
                options.include_chirality
            },
            use_atom_maps: options.include_atom_maps,
            use_non_stereo_ranks: options.use_non_stereo_ranks,
            use_chiral_presence: options.include_chiral_presence,
            use_atom_maps_on_dummies: true,
        }
    }
}

pub(crate) fn rank_mol_atoms(molecule: &Molecule) -> Result<Vec<usize>, KekulizeError> {
    let view = CanonRankReadView::from_molecule(molecule)?;
    rank_mol_atoms_with_options_from_view(&view, CanonicalRankOptions::rank_mol_default())
}

pub(crate) fn rank_mol_atoms_from_read_parts(
    read: MoleculeReadParts<'_>,
) -> Result<Vec<usize>, KekulizeError> {
    let view = CanonRankReadView::from_read_parts(read)?;
    rank_mol_atoms_with_options_from_view(&view, CanonicalRankOptions::rank_mol_default())
}

pub(crate) fn rank_mol_atoms_with_options(
    molecule: &Molecule,
    options: CanonicalRankOptions,
) -> Result<Vec<usize>, KekulizeError> {
    let view = CanonRankReadView::from_molecule(molecule)?;
    rank_mol_atoms_with_options_from_view(&view, options)
}

fn rank_mol_atoms_with_options_from_view(
    view: &CanonRankReadView<'_>,
    options: CanonicalRankOptions,
) -> Result<Vec<usize>, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: void rankMolAtoms(const ROMol&, std::vector<unsigned int>&, bool, bool, bool, bool, bool, bool, bool, bool)
    // RDKit✔️✔️: void rankMolAtoms(const ROMol &mol, std::vector<unsigned int> &res,
    // RDKit✔️✔️:                   bool breakTies, bool includeChirality, bool includeIsotopes,
    // RDKit✔️✔️:                   bool includeAtomMaps, bool includeChiralPresence,
    // RDKit✔️✔️:                   bool includeStereoGroups, bool useNonStereoRanks,
    // RDKit✔️✔️:                   bool includeRingStereo) {
    // RDKit✔️✔️:   if (!mol.getNumAtoms()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if view.num_atoms() == 0 {
        return Ok(Vec::new());
    }

    // RDKit✔️✔️:   bool clearRings = false;
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:     clearRings = true;
    // RDKit✔️✔️:   }
    // COSMolKit rank computation is read-only, so ring info is computed ephemerally.
    // RDKit✔️✔️:   res.resize(mol.getNumAtoms());
    // RDKit✔️✔️:   std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
    // RDKit✔️✔️:   initCanonAtoms(mol, atoms, includeChirality, includeStereoGroups);
    let mut atoms = init_canon_atoms_for_rank_mol(
        view,
        options.include_chirality,
        options.include_stereo_groups,
    )?;

    // RDKit✔️✔️:   AtomCompareFunctor ftor(&atoms.front(), mol);
    // RDKit✔️✔️:   ftor.df_useIsotopes = includeIsotopes;
    // RDKit✔️✔️:   ftor.df_useChirality = includeChirality;
    // RDKit✔️✔️:   ftor.df_useChiralityRings = includeChirality && includeRingStereo;
    // RDKit✔️✔️:   ftor.df_useAtomMaps = includeAtomMaps;
    // RDKit✔️✔️:   ftor.df_useNonStereoRanks = useNonStereoRanks;
    // RDKit✔️✔️:   ftor.df_useChiralPresence = includeChiralPresence;
    let mut order = vec![0usize; view.num_atoms()];
    let flags = CanonRankFlags::from_fragment_options(options);
    // RDKit✔️✔️:   detail::rankWithFunctor(ftor, breakTies, order, true, includeChirality,
    // RDKit✔️✔️:                           includeRingStereo);
    rank_with_atom_compare_functor_for_kekulize(
        view,
        &mut atoms,
        options.break_ties,
        options.include_ring_stereo,
        flags,
        &mut order,
    )?;

    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:     res[order[i]] = atoms[order[i]].index;
    // RDKit✔️✔️:   }
    let mut res = vec![0usize; view.num_atoms()];
    for idx in 0..view.num_atoms() {
        res[order[idx]] = usize::try_from(atoms[order[idx]].index).unwrap_or(usize::MAX);
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: void rankMolAtoms(const ROMol&, std::vector<unsigned int>&, bool, bool, bool, bool, bool, bool, bool, bool)
    Ok(res)
}

pub(crate) fn rank_fragment_atoms_for_kekulize(
    topology: &TopologyBlock,
    stereo_groups: &[StereoGroup],
    scope: FragmentRankScope<'_>,
) -> Result<Vec<usize>, KekulizeError> {
    let view = CanonRankReadView::from_topology(topology, stereo_groups)?;
    rank_fragment_atoms_from_view(&view, scope, CanonicalRankOptions::kekulize_default())
}

pub(crate) fn rank_fragment_atoms(
    molecule: &Molecule,
    scope: FragmentRankScope<'_>,
    options: CanonicalRankOptions,
) -> Result<Vec<usize>, KekulizeError> {
    let view = CanonRankReadView::from_molecule(molecule)?;
    rank_fragment_atoms_from_view(&view, scope, options)
}

fn rank_fragment_atoms_from_view(
    view: &CanonRankReadView<'_>,
    scope: FragmentRankScope<'_>,
    options: CanonicalRankOptions,
) -> Result<Vec<usize>, KekulizeError> {
    let atoms_to_use = scope.atoms_to_use;
    let bonds_to_use = scope.bonds_to_use;
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: void rankFragmentAtoms(const ROMol&, std::vector<unsigned int>&, const boost::dynamic_bitset<>&, const boost::dynamic_bitset<>&, const std::vector<std::string>*, const std::vector<std::string>*, bool, bool, bool, bool, bool, bool)
    // RDKit✔️✔️: void rankFragmentAtoms(const ROMol &mol, std::vector<unsigned int> &res,
    // RDKit✔️✔️:                        const boost::dynamic_bitset<> &atomsInPlay,
    // RDKit✔️✔️:                        const boost::dynamic_bitset<> &bondsInPlay,
    // RDKit✔️✔️:                        const std::vector<std::string> *atomSymbols,
    // RDKit✔️✔️:                        const std::vector<std::string> *bondSymbols,
    // RDKit✔️✔️:                        bool breakTies, bool includeChirality,
    // RDKit✔️✔️:                        bool includeIsotopes, bool includeAtomMaps,
    // RDKit✔️✔️:                        bool includeChiralPresence, bool includeRingStereo) {
    // RDKit✔️✔️:   PRECONDITION(atomsInPlay.size() == mol.getNumAtoms(), "bad atomsInPlay size");
    // RDKit✔️✔️:   PRECONDITION(bondsInPlay.size() == mol.getNumBonds(), "bad bondsInPlay size");
    // RDKit✔️✔️:   PRECONDITION(!atomSymbols || atomSymbols->size() == mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atomSymbols size");
    // RDKit✔️✔️:   PRECONDITION(!bondSymbols || bondSymbols->size() == mol.getNumBonds(),
    // RDKit✔️✔️:                "bad bondSymbols size");
    // RDKit✔️✔️:   if (!mol.getNumAtoms()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool clearRings = false;
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:     clearRings = true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res.resize(mol.getNumAtoms());
    // RDKit✔️✔️:   std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
    // RDKit✔️✔️:   detail::initFragmentCanonAtoms(mol, atoms, includeChirality, atomSymbols,
    // RDKit✔️✔️:                                  bondSymbols, atomsInPlay, bondsInPlay, true);
    // RDKit✔️✔️:   AtomCompareFunctor ftor(&atoms.front(), mol, &atomsInPlay, &bondsInPlay);
    // RDKit✔️✔️:   ftor.df_useIsotopes = includeIsotopes;
    // RDKit✔️✔️:   ftor.df_useChirality = includeChirality;
    // RDKit✔️✔️:   ftor.df_useAtomMaps = includeAtomMaps;
    // RDKit✔️✔️:   ftor.df_useChiralityRings = includeChirality;
    // RDKit✔️✔️:   ftor.df_useChiralPresence = includeChiralPresence;
    // RDKit✔️✔️:   std::vector<int> order(mol.getNumAtoms());
    // RDKit✔️✔️:   detail::rankWithFunctor(ftor, breakTies, order, true, includeChirality,
    // RDKit✔️✔️:                           includeRingStereo, &atomsInPlay, &bondsInPlay);
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:     res[order[i]] = atoms[order[i]].index;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (clearRings) {
    // RDKit✔️✔️:     mol.getRingInfo()->reset();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: void rankFragmentAtoms(const ROMol&, std::vector<unsigned int>&, const boost::dynamic_bitset<>&, const boost::dynamic_bitset<>&, const std::vector<std::string>*, const std::vector<std::string>*, bool, bool, bool, bool, bool, bool)
    if atoms_to_use.len() != view.num_atoms() || bonds_to_use.len() != view.num_bonds() {
        return Err(KekulizeError::FragmentBitsetSizeMismatch {
            atoms: atoms_to_use.len(),
            bonds: bonds_to_use.len(),
        });
    }
    if let Some(atom_symbols) = scope.atom_symbols
        && atom_symbols.len() != view.num_atoms()
    {
        return Err(KekulizeError::CanonicalRankSymbolSizeMismatch {
            kind: "atomSymbols",
            expected: view.num_atoms(),
            actual: atom_symbols.len(),
        });
    }
    if let Some(bond_symbols) = scope.bond_symbols
        && bond_symbols.len() != view.num_bonds()
    {
        return Err(KekulizeError::CanonicalRankSymbolSizeMismatch {
            kind: "bondSymbols",
            expected: view.num_bonds(),
            actual: bond_symbols.len(),
        });
    }
    if view.num_atoms() == 0 {
        return Ok(Vec::new());
    }
    // RDKit✔️✔️: rankFragmentAtoms() has no includeStereoGroups parameter and
    // RDKit✔️✔️: detail::initFragmentCanonAtoms() does not assign stereo groups.
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        view,
        atoms_to_use,
        bonds_to_use,
        options.include_chirality,
        false,
        scope.atom_symbols,
        scope.bond_symbols,
    )?;
    let mut order = vec![0usize; view.num_atoms()];
    let flags = CanonRankFlags::from_fragment_options(options);
    rank_with_atom_compare_functor_for_kekulize(
        view,
        &mut atoms,
        options.break_ties,
        options.include_ring_stereo,
        flags,
        &mut order,
    )?;
    let mut res = vec![0usize; view.num_atoms()];
    for idx in 0..view.num_atoms() {
        res[order[idx]] = usize::try_from(atoms[order[idx]].index).unwrap_or(usize::MAX);
    }
    Ok(res)
}

#[derive(Debug, Clone)]
struct CanonBondHolder<'a> {
    bond_type: BondOrder,
    bond_stereo: u8,
    stype: BondStereo,
    controlling_atoms: [Option<usize>; 4],
    nbr_sym_class: usize,
    nbr_idx: usize,
    p_symbol: Option<&'a str>,
    // Source-aligned storage for RDKit's needsInit=false symbol update branch.
    // rankFragmentAtoms forces needsInit=true today, so this remains dormant.
    #[allow(dead_code)]
    bond_idx: usize,
}

#[derive(Debug, Clone)]
struct CanonAtom<'a> {
    index: i32,
    is_in_play: bool,
    degree: usize,
    atomic_number: u8,
    isotope: u16,
    atom_map: u32,
    canonical_ranking_number: i32,
    formal_charge: i8,
    chiral_tag: ChiralTag,
    total_num_hs: usize,
    is_ring_atom: bool,
    has_ring_nbr: bool,
    is_ring_stereo_atom: bool,
    which_stereo_group: usize,
    type_of_stereo_group: CanonStereoGroupType,
    neighbor_num: Vec<i32>,
    revisted_neighbors: Vec<i32>,
    all_nbr_ids: Vec<usize>,
    nbr_ids: Vec<usize>,
    p_symbol: Option<&'a str>,
    bonds: Vec<CanonBondHolder<'a>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum CanonStereoGroupType {
    Absolute = 0,
    Or = 1,
    And = 2,
}

#[derive(Debug, Clone, Copy)]
enum CanonCompareMode {
    Atom,
    SpecialChirality,
    SpecialSymmetry,
}

fn atropisomer_atoms_and_bonds_for_canonical_rank(
    view: &CanonRankReadView<'_>,
    bond_id: BondId,
) -> Option<[(AtomId, Vec<BondId>); 2]> {
    // BEGIN RDKIT CPP FUNCTION Atropisomers::getAtropisomerAtomsAndBonds
    // RDKit✔️✔️: bool getAtropisomerAtomsAndBonds(const Bond *bond,
    // RDKit✔️✔️:                                  AtropAtomAndBondVec atomsAndBondVects[2],
    // RDKit✔️✔️:                                  const ROMol &mol) {
    // RDKit✔️✔️:   PRECONDITION(bond, "no bond");
    let bond = view.bonds.get(bond_id.index())?;
    // RDKit✔️✔️:   atomsAndBondVects[0].first = bond->getBeginAtom();
    // RDKit✔️✔️:   atomsAndBondVects[1].first = bond->getEndAtom();
    let atoms = [bond.begin(), bond.end()];
    let mut result = [(atoms[0], Vec::new()), (atoms[1], Vec::new())];
    // RDKit✔️✔️:   for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    for bond_atom_index in 0..2 {
        // RDKit✔️✔️:     for (const auto nbrBond :
        // RDKit✔️✔️:          mol.atomBonds(atomsAndBondVects[bondAtomIndex].first)) {
        for neighbor_bond in view.bonds {
            if neighbor_bond.begin() != atoms[bond_atom_index]
                && neighbor_bond.end() != atoms[bond_atom_index]
            {
                continue;
            }
            // RDKit✔️✔️:       if (nbrBond == bond) {
            // RDKit✔️✔️:         continue;
            // RDKit✔️✔️:       }
            if neighbor_bond.id() == bond_id {
                continue;
            }
            // RDKit✔️✔️:       atomsAndBondVects[bondAtomIndex].second.push_back(nbrBond);
            result[bond_atom_index].1.push(neighbor_bond.id());
        }
        // RDKit✔️✔️:     if (atomsAndBondVects[bondAtomIndex].second.size() == 0) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        if result[bond_atom_index].1.is_empty() {
            return None;
        }
        // RDKit✔️✔️:     if (atomsAndBondVects[bondAtomIndex].second.size() == 2 &&
        // RDKit✔️✔️:         atomsAndBondVects[bondAtomIndex]
        // RDKit✔️✔️:                 .second[1]
        // RDKit✔️✔️:                 ->getOtherAtom(atomsAndBondVects[bondAtomIndex].first)
        // RDKit✔️✔️:                 ->getIdx() <
        // RDKit✔️✔️:             atomsAndBondVects[bondAtomIndex]
        // RDKit✔️✔️:                 .second[0]
        // RDKit✔️✔️:                 ->getOtherAtom(atomsAndBondVects[bondAtomIndex].first)
        // RDKit✔️✔️:                 ->getIdx()) {
        // RDKit✔️✔️:       std::swap(atomsAndBondVects[bondAtomIndex].second[0],
        // RDKit✔️✔️:                 atomsAndBondVects[bondAtomIndex].second[1]);
        // RDKit✔️✔️:     }
        if result[bond_atom_index].1.len() == 2 {
            let first_other =
                bond_other_atom_index(view, result[bond_atom_index].1[0], atoms[bond_atom_index])?;
            let second_other =
                bond_other_atom_index(view, result[bond_atom_index].1[1], atoms[bond_atom_index])?;
            if second_other < first_other {
                result[bond_atom_index].1.swap(0, 1);
            }
        }
    }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atropisomers::getAtropisomerAtomsAndBonds
    Some(result)
}

fn bond_other_atom_index(
    view: &CanonRankReadView<'_>,
    bond_id: BondId,
    atom_id: AtomId,
) -> Option<usize> {
    view.bond_other_atom_index(bond_id, atom_id)
}

fn count_swaps_to_interconvert<T: Copy + Eq>(reference: &[T], probe: Vec<T>) -> usize {
    crate::source_port_helpers::count_swaps_to_interconvert(reference, &probe).unwrap_or_else(
        |error| match error {
            crate::source_port_helpers::CountSwapsError::SizeMismatch => panic!("size mismatch"),
            crate::source_port_helpers::CountSwapsError::MissingProbeElement => {
                panic!("could not find probe element")
            }
        },
    )
}

fn empty_canon_atom_from_source_atom(atom: &Atom) -> CanonAtom<'static> {
    CanonAtom {
        index: i32::try_from(atom.id().index()).unwrap_or(i32::MAX),
        is_in_play: true,
        degree: 0,
        atomic_number: atom.atomic_number(),
        isotope: atom.isotope().unwrap_or(0),
        atom_map: atom.atom_map().unwrap_or(0),
        canonical_ranking_number: atom
            .prop("_CanonicalRankingNumber")
            .and_then(|value| value.parse::<i32>().ok())
            .unwrap_or(0),
        formal_charge: atom.formal_charge(),
        chiral_tag: atom.chiral_tag(),
        total_num_hs: 0,
        is_ring_atom: false,
        has_ring_nbr: false,
        is_ring_stereo_atom: false,
        which_stereo_group: 0,
        type_of_stereo_group: CanonStereoGroupType::Absolute,
        neighbor_num: Vec::new(),
        revisted_neighbors: Vec::new(),
        all_nbr_ids: Vec::new(),
        nbr_ids: Vec::new(),
        p_symbol: None,
        bonds: Vec::new(),
    }
}

fn init_canon_atoms_for_rank_mol(
    view: &CanonRankReadView<'_>,
    include_chirality: bool,
    include_stereo_groups: bool,
) -> Result<Vec<CanonAtom<'static>>, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: void initCanonAtoms(const ROMol&, std::vector<Canon::canon_atom>&, bool, bool)
    // RDKit✔️✔️: void initCanonAtoms(const ROMol &mol, std::vector<Canon::canon_atom> &atoms,
    // RDKit✔️✔️:                     bool includeChirality, bool includeStereoGroups) {
    let mut atoms = view
        .atoms
        .iter()
        .map(empty_canon_atom_from_source_atom)
        .collect::<Vec<_>>();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    for atom_idx in 0..view.num_atoms() {
        // RDKit✔️✔️:     basicInitCanonAtom(mol, atoms[i], i);
        // RDKit✔️✔️:   atom.atom = mol.getAtomWithIdx(idx);
        // RDKit✔️✔️:   atom.index = idx;
        // RDKit✔️✔️:   atom.p_symbol = nullptr;
        // RDKit✔️✔️:   atom.degree = atom.atom->getDegree();
        // RDKit✔️✔️:   atom.nbrIds = std::make_unique<int[]>(atom.degree);
        // RDKit✔️✔️:   getNbrs(mol, atom.atom, atom.nbrIds.get());
        let neighbors = view.adjacency.neighbors_of(atom_idx);
        let degree = neighbors.len();
        let nbr_ids = neighbors
            .iter()
            .map(|neighbor| neighbor.atom_index)
            .collect::<Vec<_>>();
        atoms[atom_idx].degree = degree;
        atoms[atom_idx].all_nbr_ids.clone_from(&nbr_ids);
        atoms[atom_idx].nbr_ids = nbr_ids;

        // RDKit✔️✔️:     advancedInitCanonAtom(mol, atoms[i], i);
        // RDKit✔️✔️:   atom.totalNumHs = atom.atom->getTotalNumHs();
        let atom = &view.atoms[atom_idx];
        atoms[atom_idx].total_num_hs = usize::from(atom.explicit_hydrogens())
            + usize::try_from(view.valence.implicit_hydrogens[atom_idx].max(0))
                .unwrap_or(usize::MAX);
        // RDKit✔️✔️:   atom.isRingStereoAtom =
        // RDKit✔️✔️:       (atom.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
        // RDKit✔️✔️:        atom.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) &&
        // RDKit✔️✔️:       atom.atom->hasProp(common_properties::_ringStereoAtoms);
        atoms[atom_idx].is_ring_atom = view.rings.num_atom_rings(atom.id()) > 0;
        atoms[atom_idx].is_ring_stereo_atom = matches!(
            atom.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) && atom.prop("_ringStereoAtoms").is_some();
        // RDKit✔️✔️:   atom.hasRingNbr = hasRingNbr(mol, atom.atom);
        atoms[atom_idx].has_ring_nbr = has_ring_nbr_for_kekulize(view, atom_idx);

        // RDKit✔️✔️:     atoms[i].bonds.reserve(atoms[i].degree);
        atoms[atom_idx].bonds.reserve(degree);
        // RDKit✔️✔️:     getBonds(mol, atoms[i].atom, atoms[i].bonds, includeChirality, atoms);
        for neighbor in neighbors {
            atoms[atom_idx].bonds.push(make_canon_bond_holder(
                view,
                neighbor.bond.index(),
                neighbor.atom_index,
                include_chirality,
            )?);
        }
        // RDKit✔️✔️:   std::sort(nbrs.begin(), nbrs.end(), bondholder::greater);
        let initial_ranks = canon_atom_rank_snapshot(&atoms);
        sort_canon_bonds_descending(&mut atoms[atom_idx].bonds, &initial_ranks);
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   if (includeChirality && includeStereoGroups) {
    // RDKit✔️✔️:     unsigned int sgidx = 1;
    // RDKit✔️✔️:     for (auto &sg : mol.getStereoGroups()) {
    // RDKit✔️✔️:       for (auto atom : sg.getAtoms()) {
    // RDKit✔️✔️:         atoms[atom->getIdx()].whichStereoGroup = sgidx;
    // RDKit✔️✔️:         atoms[atom->getIdx()].typeOfStereoGroup = sg.getGroupType();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ++sgidx;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if include_chirality && include_stereo_groups {
        for (group_idx, group) in view.stereo_groups.iter().enumerate() {
            for atom in group.atoms() {
                atoms[atom.index()].which_stereo_group = group_idx + 1;
                atoms[atom.index()].type_of_stereo_group =
                    canon_stereo_group_type_from_kind(group.kind());
            }
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: void initCanonAtoms(const ROMol&, std::vector<Canon::canon_atom>&, bool, bool)
    Ok(atoms)
}

fn init_fragment_canon_atoms_for_kekulize<'a>(
    view: &CanonRankReadView<'_>,
    atoms_to_use: &[bool],
    bonds_to_use: &[bool],
    include_chirality: bool,
    include_stereo_groups: bool,
    atom_symbols: Option<&'a [String]>,
    bond_symbols: Option<&'a [String]>,
) -> Result<Vec<CanonAtom<'a>>, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION detail::initFragmentCanonAtoms
    // RDKit✔️✔️: void initFragmentCanonAtoms(const ROMol &mol,
    // RDKit✔️✔️:                             std::vector<Canon::canon_atom> &atoms,
    // RDKit✔️✔️:                             bool includeChirality,
    // RDKit✔️✔️:                             const std::vector<std::string> *atomSymbols,
    // RDKit✔️✔️:                             const std::vector<std::string> *bondSymbols,
    // RDKit✔️✔️:                             const boost::dynamic_bitset<> &atomsInPlay,
    // RDKit✔️✔️:                             const boost::dynamic_bitset<> &bondsInPlay,
    // RDKit✔️✔️:                             bool needsInit) {
    // RDKit✔️✔️:   needsInit = true;
    // RDKit✔️✔️:   PRECONDITION(!atomSymbols || atomSymbols->size() == mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atom symbols");
    // RDKit✔️✔️:   PRECONDITION(!bondSymbols || bondSymbols->size() == mol.getNumBonds(),
    // RDKit✔️✔️:                "bad bond symbols");
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     auto i = atom->getIdx();
    // RDKit✔️✔️:     auto &atomsi = atoms[i];
    // RDKit✔️✔️:     atomsi.atom = atom;
    // RDKit✔️✔️:     atomsi.index = i;
    // RDKit✔️✔️:     atomsi.degree = 0;
    let mut atoms = view
        .atoms
        .iter()
        .map(|atom| CanonAtom {
            index: i32::try_from(atom.id().index()).unwrap_or(i32::MAX),
            is_in_play: false,
            degree: 0,
            atomic_number: atom.atomic_number(),
            isotope: atom.isotope().unwrap_or(0),
            atom_map: atom.atom_map().unwrap_or(0),
            canonical_ranking_number: atom
                .prop("_CanonicalRankingNumber")
                .and_then(|value| value.parse::<i32>().ok())
                .unwrap_or(0),
            formal_charge: atom.formal_charge(),
            chiral_tag: atom.chiral_tag(),
            total_num_hs: 0,
            is_ring_atom: false,
            has_ring_nbr: false,
            is_ring_stereo_atom: false,
            which_stereo_group: 0,
            type_of_stereo_group: CanonStereoGroupType::Absolute,
            neighbor_num: Vec::new(),
            revisted_neighbors: Vec::new(),
            all_nbr_ids: Vec::new(),
            nbr_ids: Vec::new(),
            p_symbol: None,
            bonds: Vec::new(),
        })
        .collect::<Vec<_>>();
    for bond in view.bonds {
        let begin = bond.begin().index();
        let end = bond.end().index();
        atoms[begin].all_nbr_ids.push(end);
        atoms[end].all_nbr_ids.push(begin);
    }
    // RDKit✔️✔️:     if (atomsInPlay[i]) {
    // RDKit✔️✔️:       if (atomSymbols) {
    // RDKit✔️✔️:         atomsi.p_symbol = &(*atomSymbols)[i];
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         atomsi.p_symbol = nullptr;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (needsInit) {
    // RDKit✔️✔️:         atomsi.nbrIds = std::make_unique<int[]>(atom->getDegree());
    // RDKit✔️✔️:         advancedInitCanonAtom(mol, atomsi, i);
    // RDKit✔️✔️:         atomsi.bonds.reserve(4);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for idx in 0..atoms.len() {
        atoms[idx].is_in_play = atoms_to_use[idx];
        if atoms_to_use[idx] {
            atoms[idx].p_symbol = atom_symbols.map(|symbols| symbols[idx].as_str());
        }
    }

    // RDKit✔️✔️:   if (needsInit) {
    // RDKit✔️✔️:     for (const auto bond : mol.bonds()) {
    // RDKit✔️✔️:       if (!bondsInPlay[bond->getIdx()] ||
    // RDKit✔️✔️:           !atomsInPlay[bond->getBeginAtomIdx()] ||
    // RDKit✔️✔️:           !atomsInPlay[bond->getEndAtomIdx()]) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       Canon::canon_atom &begAt = atoms[bond->getBeginAtomIdx()];
    // RDKit✔️✔️:       Canon::canon_atom &endAt = atoms[bond->getEndAtomIdx()];
    // RDKit✔️✔️:       begAt.nbrIds[begAt.degree++] = bond->getEndAtomIdx();
    // RDKit✔️✔️:       endAt.nbrIds[endAt.degree++] = bond->getBeginAtomIdx();
    // RDKit✔️✔️:       begAt.bonds.push_back(
    // RDKit✔️✔️:           makeBondHolder(bond, bond->getEndAtomIdx(), includeChirality, atoms));
    // RDKit✔️✔️:       endAt.bonds.push_back(makeBondHolder(bond, bond->getBeginAtomIdx(),
    // RDKit✔️✔️:                                            includeChirality, atoms));
    for bond in view.bonds {
        let bond_idx = bond.id().index();
        let begin = bond.begin().index();
        let end = bond.end().index();
        if !bonds_to_use[bond_idx] || !atoms_to_use[begin] || !atoms_to_use[end] {
            continue;
        }
        atoms[begin].nbr_ids.push(end);
        atoms[begin].degree += 1;
        atoms[end].nbr_ids.push(begin);
        atoms[end].degree += 1;
        let mut begin_bond = make_canon_bond_holder(view, bond_idx, end, include_chirality)?;
        let mut end_bond = make_canon_bond_holder(view, bond_idx, begin, include_chirality)?;
        if let Some(symbols) = bond_symbols {
            begin_bond.p_symbol = Some(symbols[bond_idx].as_str());
            end_bond.p_symbol = Some(symbols[bond_idx].as_str());
        }
        atoms[begin].bonds.push(begin_bond);
        atoms[end].bonds.push(end_bond);
    }
    // RDKit✔️✔️:       if (bondSymbols) {
    // RDKit✔️✔️:         begAt.bonds.back().p_symbol = &(*bondSymbols)[bond->getIdx()];
    // RDKit✔️✔️:         endAt.bonds.back().p_symbol = &(*bondSymbols)[bond->getIdx()];
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if (bondSymbols) {
    // RDKit✔️✔️:       for (auto &atom : atoms) {
    // RDKit✔️✔️:         for (auto &bond : atom.bonds) {
    // RDKit✔️✔️:           bond.p_symbol = &(*bondSymbols)[bond.bondIdx];
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }

    // RDKit✔️✔️:   for (size_t i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:     if (!atomsInPlay[i]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto &atomsi = atoms[i];
    // RDKit✔️✔️:     if (needsInit) {
    // RDKit✔️✔️:       atomsi.totalNumHs += (mol.getAtomWithIdx(i)->getDegree() - atomsi.degree);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::sort(atomsi.bonds.begin(), atomsi.bonds.end(), bondholder::greater);
    // RDKit✔️✔️:   }
    let initial_ranks = canon_atom_rank_snapshot(&atoms);
    for atom in view.atoms {
        let idx = atom.id().index();
        if !atoms_to_use[idx] {
            continue;
        }
        let original_degree = view.atom_degree(atom.id());
        atoms[idx].total_num_hs = usize::from(atom.explicit_hydrogens())
            + usize::try_from(view.valence.implicit_hydrogens[idx].max(0)).unwrap_or(usize::MAX)
            + original_degree.saturating_sub(atoms[idx].degree);
        sort_canon_bonds_descending(&mut atoms[idx].bonds, &initial_ranks);
    }
    advanced_init_canon_atoms_stereo_state_for_kekulize(
        view,
        &mut atoms,
        atoms_to_use,
        include_chirality,
        include_stereo_groups,
    )?;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION detail::initFragmentCanonAtoms
    Ok(atoms)
}

fn advanced_init_canon_atoms_stereo_state_for_kekulize(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
    atoms_to_use: &[bool],
    include_chirality: bool,
    include_stereo_groups: bool,
) -> Result<(), KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION advancedInitCanonAtom
    // RDKit✔️✔️: void advancedInitCanonAtom(const ROMol &mol, Canon::canon_atom &atom,
    // RDKit✔️✔️:                            const int &) {
    // RDKit✔️✔️:   atom.totalNumHs = atom.atom->getTotalNumHs();
    // RDKit✔️✔️:   atom.isRingStereoAtom =
    // RDKit✔️✔️:       (atom.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:        atom.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) &&
    // RDKit✔️✔️:       atom.atom->hasProp(common_properties::_ringStereoAtoms);
    // RDKit✔️✔️:   atom.hasRingNbr = hasRingNbr(mol, atom.atom);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION advancedInitCanonAtom
    for atom in view.atoms {
        let idx = atom.id().index();
        if !atoms_to_use[idx] {
            continue;
        }
        atoms[idx].is_ring_atom = view.rings.num_atom_rings(atom.id()) > 0;
        atoms[idx].is_ring_stereo_atom = matches!(
            atom.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) && atom.prop("_ringStereoAtoms").is_some();
        atoms[idx].has_ring_nbr = has_ring_nbr_for_kekulize(view, idx);
    }

    // BEGIN RDKIT CPP FUNCTION initCanonAtoms stereo group assignment
    // RDKit✔️✔️:   if (includeChirality && includeStereoGroups) {
    // RDKit✔️✔️:     unsigned int sgidx = 1;
    // RDKit✔️✔️:     for (auto &sg : mol.getStereoGroups()) {
    // RDKit✔️✔️:       for (auto atom : sg.getAtoms()) {
    // RDKit✔️✔️:         atoms[atom->getIdx()].whichStereoGroup = sgidx;
    // RDKit✔️✔️:         atoms[atom->getIdx()].typeOfStereoGroup = sg.getGroupType();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ++sgidx;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION initCanonAtoms stereo group assignment
    if include_chirality && include_stereo_groups {
        for (group_idx, group) in view.stereo_groups.iter().enumerate() {
            let rdkit_group_idx = group_idx + 1;
            for atom in group.atoms() {
                let atom_idx = atom.index();
                if atom_idx < atoms.len() {
                    atoms[atom_idx].which_stereo_group = rdkit_group_idx;
                    atoms[atom_idx].type_of_stereo_group =
                        canon_stereo_group_type_from_kind(group.kind());
                }
            }
        }
    }
    Ok(())
}

fn has_ring_nbr_for_kekulize(view: &CanonRankReadView<'_>, atom_idx: usize) -> bool {
    // BEGIN RDKIT CPP FUNCTION hasRingNbr
    // RDKit✔️✔️: bool hasRingNbr(const ROMol &mol, const Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad pointer");
    // RDKit✔️✔️:   for (const auto nbr : mol.atomNeighbors(at)) {
    // RDKit✔️✔️:     if ((nbr->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:          nbr->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) &&
    // RDKit✔️✔️:         nbr->hasProp(common_properties::_ringStereoAtoms)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION hasRingNbr
    view.adjacency
        .neighbors_of(atom_idx)
        .iter()
        .any(|neighbor_ref| {
            let neighbor = &view.atoms[neighbor_ref.atom_index];
            matches!(
                neighbor.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            ) && neighbor.prop("_ringStereoAtoms").is_some()
        })
}

const fn canon_stereo_group_type_from_kind(kind: StereoGroupKind) -> CanonStereoGroupType {
    match kind {
        StereoGroupKind::Absolute => CanonStereoGroupType::Absolute,
        StereoGroupKind::Or => CanonStereoGroupType::Or,
        StereoGroupKind::And => CanonStereoGroupType::And,
    }
}

fn make_canon_bond_holder<'a>(
    view: &CanonRankReadView<'_>,
    bond_idx: usize,
    other_idx: usize,
    include_chirality: bool,
) -> Result<CanonBondHolder<'a>, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION makeBondHolder
    // RDKit✔️✔️: bondholder makeBondHolder(const Bond *bond, unsigned int otherIdx,
    // RDKit✔️✔️:                           bool includeChirality,
    // RDKit✔️✔️:                           const std::vector<Canon::canon_atom> &atoms) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad pointer");
    // RDKit✔️✔️:   Bond::BondStereo stereo = Bond::STEREONONE;
    // RDKit✔️✔️:   if (includeChirality) {
    // RDKit✔️✔️:     stereo = bond->getStereo();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Bond::BondType bt =
    // RDKit✔️✔️:       bond->getIsAromatic() ? Bond::AROMATIC : bond->getBondType();
    // RDKit✔️✔️:   bondholder res(bt, stereo, otherIdx, 0, bond->getIdx());
    let bond = &view.bonds[bond_idx];
    let bond_type = if bond.is_aromatic() {
        BondOrder::Aromatic
    } else {
        bond.order()
    };
    let mut controlling_atoms = [None; 4];
    // RDKit✔️✔️:   if (includeChirality) {
    // RDKit✔️✔️:     res.stype = bond->getStereo();
    // RDKit✔️✔️:     if (res.stype == Bond::BondStereo::STEREOCIS ||
    // RDKit✔️✔️:         res.stype == Bond::BondStereo::STEREOTRANS) {
    // RDKit✔️✔️:       res.controllingAtoms[0] = &atoms[bond->getStereoAtoms()[0]];
    // RDKit✔️✔️:       res.controllingAtoms[2] = &atoms[bond->getStereoAtoms()[1]];
    if include_chirality && matches!(bond.stereo(), BondStereo::Cis | BondStereo::Trans) {
        let Some(stereo_atoms) = bond.stereo_atoms() else {
            return Err(KekulizeError::ProtocolDebt {
                branch: "makeBondHolder cis/trans stereo atoms",
                reason: "cis/trans bond stereo requires the two RDKit stereo atom references",
            });
        };
        controlling_atoms[0] = Some(stereo_atoms[0].index());
        controlling_atoms[2] = Some(stereo_atoms[1].index());
        // RDKit✔️✔️:       if (bond->getBeginAtom()->getDegree() > 2) {
        // RDKit✔️✔️:         for (const auto nbr :
        // RDKit✔️✔️:              bond->getOwningMol().atomNeighbors(bond->getBeginAtom())) {
        // RDKit✔️✔️:           if (nbr->getIdx() != bond->getEndAtomIdx() &&
        // RDKit✔️✔️:               nbr->getIdx() !=
        // RDKit✔️✔️:                   static_cast<unsigned int>(bond->getStereoAtoms()[0])) {
        // RDKit✔️✔️:             res.controllingAtoms[1] = &atoms[nbr->getIdx()];
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        if view.atom_degree(bond.begin()) > 2 {
            for neighbor in view.atom_neighbors(bond.begin()) {
                if neighbor != bond.end().index() && neighbor != stereo_atoms[0].index() {
                    controlling_atoms[1] = Some(neighbor);
                }
            }
        }
        // RDKit✔️✔️:       if (bond->getEndAtom()->getDegree() > 2) {
        // RDKit✔️✔️:         for (const auto nbr :
        // RDKit✔️✔️:              bond->getOwningMol().atomNeighbors(bond->getEndAtom())) {
        // RDKit✔️✔️:           if (nbr->getIdx() != bond->getBeginAtomIdx() &&
        // RDKit✔️✔️:               nbr->getIdx() !=
        // RDKit✔️✔️:                   static_cast<unsigned int>(bond->getStereoAtoms()[1])) {
        // RDKit✔️✔️:             res.controllingAtoms[3] = &atoms[nbr->getIdx()];
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        if view.atom_degree(bond.end()) > 2 {
            for neighbor in view.atom_neighbors(bond.end()) {
                if neighbor != bond.begin().index() && neighbor != stereo_atoms[1].index() {
                    controlling_atoms[3] = Some(neighbor);
                }
            }
        }
    }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (res.stype == Bond::BondStereo::STEREOATROPCCW ||
    // RDKit✔️✔️:         res.stype == Bond::BondStereo::STEREOATROPCW) {
    if include_chirality && matches!(bond.stereo(), BondStereo::AtropCcw | BondStereo::AtropCw) {
        // RDKit✔️✔️:       Atropisomers::AtropAtomAndBondVec atropAtomAndBondVecs[2];
        // RDKit✔️✔️:       CHECK_INVARIANT(Atropisomers::getAtropisomerAtomsAndBonds(
        // RDKit✔️✔️:                           bond, atropAtomAndBondVecs, bond->getOwningMol()),
        // RDKit✔️✔️:                       "Could not find atropisomer controlling atoms")
        let Some(atrop) = atropisomer_atoms_and_bonds_for_canonical_rank(view, bond.id()) else {
            return Err(KekulizeError::ProtocolDebt {
                branch: "Atropisomers::getAtropisomerAtomsAndBonds invariant",
                reason: "could not find atropisomer controlling atoms",
            });
        };
        // RDKit✔️✔️:       res.controllingAtoms[0] =
        // RDKit✔️✔️:           &atoms[atropAtomAndBondVecs[0]
        // RDKit✔️✔️:                      .second[0]
        // RDKit✔️✔️:                      ->getOtherAtom(atropAtomAndBondVecs[0].first)
        // RDKit✔️✔️:                      ->getIdx()];
        controlling_atoms[0] = Some(
            bond_other_atom_index(view, atrop[0].1[0], atrop[0].0)
                .expect("atropisomer neighbor bond must include focus atom"),
        );
        // RDKit✔️✔️:       res.controllingAtoms[2] =
        // RDKit✔️✔️:           &atoms[atropAtomAndBondVecs[1]
        // RDKit✔️✔️:                      .second[0]
        // RDKit✔️✔️:                      ->getOtherAtom(atropAtomAndBondVecs[1].first)
        // RDKit✔️✔️:                      ->getIdx()];
        controlling_atoms[2] = Some(
            bond_other_atom_index(view, atrop[1].1[0], atrop[1].0)
                .expect("atropisomer neighbor bond must include focus atom"),
        );
        // RDKit✔️✔️:       if (atropAtomAndBondVecs[0].second.size() > 1) {
        // RDKit✔️✔️:         res.controllingAtoms[1] =
        // RDKit✔️✔️:             &atoms[atropAtomAndBondVecs[0]
        // RDKit✔️✔️:                        .second[1]
        // RDKit✔️✔️:                        ->getOtherAtom(atropAtomAndBondVecs[0].first)
        // RDKit✔️✔️:                        ->getIdx()];
        // RDKit✔️✔️:       }
        if atrop[0].1.len() > 1 {
            controlling_atoms[1] = Some(
                bond_other_atom_index(view, atrop[0].1[1], atrop[0].0)
                    .expect("atropisomer neighbor bond must include focus atom"),
            );
        }
        // RDKit✔️✔️:       if (atropAtomAndBondVecs[1].second.size() > 1) {
        // RDKit✔️✔️:         res.controllingAtoms[3] =
        // RDKit✔️✔️:             &atoms[atropAtomAndBondVecs[1]
        // RDKit✔️✔️:                        .second[1]
        // RDKit✔️✔️:                        ->getOtherAtom(atropAtomAndBondVecs[1].first)
        // RDKit✔️✔️:                        ->getIdx()];
        // RDKit✔️✔️:     }
        if atrop[1].1.len() > 1 {
            controlling_atoms[3] = Some(
                bond_other_atom_index(view, atrop[1].1[1], atrop[1].0)
                    .expect("atropisomer neighbor bond must include focus atom"),
            );
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION makeBondHolder
    let stype = if include_chirality {
        bond.stereo()
    } else {
        BondStereo::None
    };
    Ok(CanonBondHolder {
        bond_type,
        bond_stereo: rdkit_bond_stereo_rank(stype),
        stype,
        controlling_atoms,
        nbr_sym_class: 0,
        nbr_idx: other_idx,
        p_symbol: None,
        bond_idx,
    })
}

fn rank_with_atom_compare_functor_for_kekulize(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
    break_ties: bool,
    include_ring_stereo: bool,
    flags: CanonRankFlags,
    order: &mut [usize],
) -> Result<(), KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION detail::rankWithFunctor
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: void rankWithFunctor(T &ftor, bool breakTies, std::vector<int> &order,
    // RDKit✔️✔️:                      bool useSpecial, bool useChirality, bool includeRingStereo,
    // RDKit✔️✔️:                      const boost::dynamic_bitset<> *atomsInPlay,
    // RDKit✔️✔️:                      const boost::dynamic_bitset<> *bondsInPlay) {
    // RDKit✔️✔️:   PRECONDITION(!order.empty(), "order should not be empty");
    // RDKit✔️✔️:   const ROMol &mol = *ftor.dp_mol;
    // RDKit✔️✔️:   canon_atom *atoms = ftor.dp_atoms;
    // RDKit✔️✔️:   const unsigned int nAts = mol.getNumAtoms();
    let n_atoms = view.num_atoms();
    // RDKit✔️✔️:   std::vector<int> count(nAts);
    // RDKit✔️✔️:   std::vector<int> next(nAts);
    // RDKit✔️✔️:   std::vector<int> changed(nAts, 1);
    // RDKit✔️✔️:   std::vector<char> touched(nAts, 0);
    // RDKit✔️✔️:   int activeset;
    let mut count = vec![0usize; n_atoms];
    let mut next = vec![-2isize; n_atoms];
    let mut changed = vec![true; n_atoms];
    let mut touched = vec![false; n_atoms];
    let mut active_set = -1isize;
    // RDKit✔️✔️:   CreateSinglePartition(nAts, order, count, atoms);
    create_single_partition_for_kekulize(n_atoms, order, &mut count, atoms);
    // RDKit✔️✔️:   ftor.df_useNbrs = true;
    // RDKit✔️✔️:   ActivatePartitions(nAts, order, count, activeset, next, changed);
    activate_partitions_for_kekulize(
        n_atoms,
        order,
        &count,
        &mut active_set,
        &mut next,
        &mut changed,
    );
    // RDKit✔️✔️:   RefinePartitions(mol, atoms, ftor, true, order, count, activeset, next,
    // RDKit✔️✔️:                    changed, touched);
    refine_partitions_for_kekulize(
        view,
        atoms,
        true,
        CanonCompareMode::Atom,
        flags,
        order,
        &mut count,
        &mut active_set,
        &mut next,
        &mut changed,
        &mut touched,
    );
    // RDKit✔️✔️:   bool ties = false;
    // RDKit✔️✔️:   for (unsigned i = 0; i < nAts; ++i) {
    // RDKit✔️✔️:     if (!count[i]) {
    // RDKit✔️✔️:       ties = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let ties = count.iter().any(|&value| value == 0);
    // RDKit✔️✔️:   if (useChirality && ties && includeRingStereo) {
    // RDKit✔️✔️:     SpecialChiralityAtomCompareFunctor scftor(atoms, mol, atomsInPlay,
    // RDKit✔️✔️:                                               bondsInPlay);
    // RDKit✔️✔️:     ActivatePartitions(nAts, order, count, activeset, next, changed);
    // RDKit✔️✔️:     RefinePartitions(mol, atoms, scftor, true, order, count, activeset, next,
    // RDKit✔️✔️:                      changed, touched);
    // RDKit✔️✔️:   }
    if flags.use_chirality && ties && include_ring_stereo {
        activate_partitions_for_kekulize(
            n_atoms,
            order,
            &count,
            &mut active_set,
            &mut next,
            &mut changed,
        );
        refine_partitions_for_kekulize(
            view,
            atoms,
            true,
            CanonCompareMode::SpecialChirality,
            flags,
            order,
            &mut count,
            &mut active_set,
            &mut next,
            &mut changed,
            &mut touched,
        );
    }
    // RDKit✔️✔️:   ties = false;
    // RDKit✔️✔️:   unsigned symRingAtoms = 0;
    // RDKit✔️✔️:   unsigned ringAtoms = 0;
    // RDKit✔️✔️:   bool branchingRingAtom = false;
    // RDKit✔️✔️:   RingInfo *ringInfo = mol.getRingInfo();
    let use_special_symmetry = special_symmetry_rank_refinement_required(view, order, &count);
    // RDKit✔️✔️:   if (useSpecial && ties && ringAtoms > 0 &&
    // RDKit✔️✔️:       static_cast<float>(symRingAtoms) / ringAtoms > 0.5 && branchingRingAtom) {
    // RDKit✔️✔️:     SpecialSymmetryAtomCompareFunctor sftor(atoms, mol, atomsInPlay,
    // RDKit✔️✔️:                                             bondsInPlay);
    // RDKit✔️✔️:     compareRingAtomsConcerningNumNeighbors(atoms, nAts, mol);
    // RDKit✔️✔️:     ActivatePartitions(nAts, order, count, activeset, next, changed);
    // RDKit✔️✔️:     RefinePartitions(mol, atoms, sftor, true, order, count, activeset, next,
    // RDKit✔️✔️:                      changed, touched);
    // RDKit✔️✔️:   }
    if use_special_symmetry {
        compare_ring_atoms_concerning_num_neighbors_for_kekulize(view, atoms);
        activate_partitions_for_kekulize(
            n_atoms,
            order,
            &count,
            &mut active_set,
            &mut next,
            &mut changed,
        );
        refine_partitions_for_kekulize(
            view,
            atoms,
            true,
            CanonCompareMode::SpecialSymmetry,
            flags,
            order,
            &mut count,
            &mut active_set,
            &mut next,
            &mut changed,
            &mut touched,
        );
    }
    // RDKit✔️✔️:   if (breakTies) {
    // RDKit✔️✔️:     BreakTies(mol, atoms, ftor, true, order, count, activeset, next, changed,
    // RDKit✔️✔️:               touched);
    // RDKit✔️✔️:   }
    if break_ties {
        break_ties_for_kekulize(
            view,
            atoms,
            true,
            CanonCompareMode::Atom,
            flags,
            order,
            &mut count,
            &mut active_set,
            &mut next,
            &mut changed,
            &mut touched,
        );
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION detail::rankWithFunctor
    Ok(())
}

fn create_single_partition_for_kekulize(
    n_atoms: usize,
    order: &mut [usize],
    count: &mut [usize],
    atoms: &mut [CanonAtom<'_>],
) {
    // BEGIN RDKIT CPP FUNCTION CreateSinglePartition
    // RDKit✔️✔️: void CreateSinglePartition(unsigned int nAtoms, std::vector<int> &order,
    // RDKit✔️✔️:                            std::vector<int> &count, canon_atom *atoms) {
    // RDKit✔️✔️:   PRECONDITION(!order.empty(), "order should not be empty");
    // RDKit✔️✔️:   PRECONDITION(!count.empty(), "count should not be empty");
    // RDKit✔️✔️:   PRECONDITION(atoms, "bad pointer");
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; i++) {
    // RDKit✔️✔️:     atoms[i].index = 0;
    // RDKit✔️✔️:     order[i] = i;
    // RDKit✔️✔️:     count[i] = 0;
    // RDKit✔️✔️:   }
    for i in 0..n_atoms {
        atoms[i].index = 0;
        order[i] = i;
        count[i] = 0;
    }
    // RDKit✔️✔️:   count[0] = nAtoms;
    // RDKit✔️✔️: }
    count[0] = n_atoms;
    // END RDKIT CPP FUNCTION CreateSinglePartition
}

fn activate_partitions_for_kekulize(
    n_atoms: usize,
    order: &[usize],
    count: &[usize],
    active_set: &mut isize,
    next: &mut [isize],
    changed: &mut [bool],
) {
    // BEGIN RDKIT CPP FUNCTION ActivatePartitions
    // RDKit✔️✔️: void ActivatePartitions(unsigned int nAtoms, std::vector<int> &order,
    // RDKit✔️✔️:                         std::vector<int> &count, int &activeset,
    // RDKit✔️✔️:                         std::vector<int> &next, std::vector<int> &changed) {
    // RDKit✔️✔️:   unsigned int i, j;
    // RDKit✔️✔️:   activeset = -1;
    *active_set = -1;
    // RDKit✔️✔️:   for (i = 0; i < nAtoms; i++) {
    // RDKit✔️✔️:     next[i] = -2;
    // RDKit✔️✔️:   }
    next.fill(-2);
    // RDKit✔️✔️:   i = 0;
    // RDKit✔️✔️:   do {
    // RDKit✔️✔️:     j = order[i];
    // RDKit✔️✔️:     if (count[j] > 1) {
    // RDKit✔️✔️:       next[j] = activeset;
    // RDKit✔️✔️:       activeset = j;
    // RDKit✔️✔️:       i += count[j];
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       i++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } while (i < nAtoms);
    let mut i = 0usize;
    while i < n_atoms {
        let j = order[i];
        if count[j] > 1 {
            next[j] = *active_set;
            *active_set = isize::try_from(j).unwrap_or(isize::MAX);
            i += count[j];
        } else {
            i += 1;
        }
    }
    // RDKit✔️✔️:   for (i = 0; i < nAtoms; i++) {
    // RDKit✔️✔️:     j = order[i];
    // RDKit✔️✔️:     int flag = 1;
    // RDKit✔️✔️:     changed[j] = flag;
    // RDKit✔️✔️:   }
    for &j in order.iter().take(n_atoms) {
        changed[j] = true;
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ActivatePartitions
}

#[allow(clippy::too_many_arguments)]
fn refine_partitions_for_kekulize(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
    mode: bool,
    compare_mode: CanonCompareMode,
    flags: CanonRankFlags,
    order: &mut [usize],
    count: &mut [usize],
    active_set: &mut isize,
    next: &mut [isize],
    changed: &mut [bool],
    touched_partitions: &mut [bool],
) {
    // BEGIN RDKIT CPP FUNCTION RefinePartitions
    // RDKit✔️✔️: template <typename CompareFunc>
    // RDKit✔️✔️: void RefinePartitions(const ROMol &mol, canon_atom *atoms, CompareFunc compar,
    // RDKit✔️✔️:                       int mode, std::vector<int> &order,
    // RDKit✔️✔️:                       std::vector<int> &count, int &activeset,
    // RDKit✔️✔️:                       std::vector<int> &next, std::vector<int> &changed,
    // RDKit✔️✔️:                       std::vector<char> &touchedPartitions) {
    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = view.num_atoms();
    // RDKit✔️✔️:   while (activeset != -1) {
    while *active_set != -1 {
        // RDKit✔️✔️:     partition = activeset;
        // RDKit✔️✔️:     activeset = next[partition];
        // RDKit✔️✔️:     next[partition] = -2;
        let partition = usize::try_from(*active_set).expect("active partition is non-negative");
        *active_set = next[partition];
        next[partition] = -2;
        // RDKit✔️✔️:     len = count[partition];
        // RDKit✔️✔️:     offset = atoms[partition].index;
        let len = count[partition];
        let offset = usize::try_from(atoms[partition].index).unwrap_or(usize::MAX);
        // RDKit✔️✔️:     auto start = std::span<int>(&order[offset], len);
        hanoi_sort_order_for_kekulize(
            order,
            offset,
            len,
            count,
            changed,
            atoms,
            compare_mode,
            flags,
        );
        // RDKit✔️✔️:     for (int k = 0; k < len; ++k) {
        // RDKit✔️✔️:       changed[start[k]] = 0;
        // RDKit✔️✔️:     }
        for k in 0..len {
            changed[order[offset + k]] = false;
        }
        // RDKit✔️✔️:     index = start[0];
        // RDKit✔️✔️:     for (i = count[index]; i < len; i++) {
        // RDKit✔️✔️:       index = start[i];
        // RDKit✔️✔️:       if (count[index]) {
        // RDKit✔️✔️:         symclass = offset + i;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       atoms[index].index = symclass;
        // RDKit✔️✔️:       for (unsigned j = 0; j < atoms[index].degree; ++j) {
        // RDKit✔️✔️:         changed[atoms[index].nbrIds[j]] = 1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mut index = order[offset];
        let mut sym_class = 0usize;
        let mut i = count[index];
        while i < len {
            index = order[offset + i];
            if count[index] > 0 {
                sym_class = offset + i;
            }
            atoms[index].index = i32::try_from(sym_class).unwrap_or(i32::MAX);
            for nbr in atoms[index].nbr_ids.iter().copied() {
                changed[nbr] = true;
            }
            i += 1;
        }
        // RDKit✔️✔️:     if (mode) {
        // RDKit✔️✔️:       index = start[0];
        // RDKit✔️✔️:       for (i = count[index]; i < len; i++) {
        // RDKit✔️✔️:         index = start[i];
        // RDKit✔️✔️:         for (unsigned j = 0; j < atoms[index].degree; ++j) {
        // RDKit✔️✔️:           unsigned int nbor = atoms[index].nbrIds[j];
        // RDKit✔️✔️:           touchedPartitions[atoms[nbor].index] = 1;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       for (unsigned int ii = 0; ii < nAtoms; ++ii) {
        // RDKit✔️✔️:         if (touchedPartitions[ii]) {
        // RDKit✔️✔️:           partition = order[ii];
        // RDKit✔️✔️:           if ((count[partition] > 1) && (next[partition] == -2)) {
        // RDKit✔️✔️:             next[partition] = activeset;
        // RDKit✔️✔️:             activeset = partition;
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:           touchedPartitions[ii] = 0;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        if mode {
            index = order[offset];
            let mut i = count[index];
            while i < len {
                index = order[offset + i];
                for nbr in atoms[index].nbr_ids.iter().copied() {
                    let partition_idx = usize::try_from(atoms[nbr].index).unwrap_or(usize::MAX);
                    if partition_idx < touched_partitions.len() {
                        touched_partitions[partition_idx] = true;
                    }
                }
                i += 1;
            }
            for ii in 0..n_atoms {
                if touched_partitions[ii] {
                    let partition = order[ii];
                    if count[partition] > 1 && next[partition] == -2 {
                        next[partition] = *active_set;
                        *active_set = isize::try_from(partition).unwrap_or(isize::MAX);
                    }
                    touched_partitions[ii] = false;
                }
            }
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RefinePartitions
}

fn break_ties_for_kekulize(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
    mode: bool,
    compare_mode: CanonCompareMode,
    flags: CanonRankFlags,
    order: &mut [usize],
    count: &mut [usize],
    active_set: &mut isize,
    next: &mut [isize],
    changed: &mut [bool],
    touched_partitions: &mut [bool],
) {
    // BEGIN RDKIT CPP FUNCTION BreakTies
    // RDKit✔️✔️: template <typename CompareFunc>
    // RDKit✔️✔️: void BreakTies(const ROMol &mol, canon_atom *atoms, CompareFunc compar,
    // RDKit✔️✔️:                int mode, std::vector<int> &order, std::vector<int> &count,
    // RDKit✔️✔️:                int &activeset, std::vector<int> &next,
    // RDKit✔️✔️:                std::vector<int> &changed,
    // RDKit✔️✔️:                std::vector<char> &touchedPartitions) {
    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = view.num_atoms();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; i++) {
    let mut i = 0usize;
    while i < n_atoms {
        // RDKit✔️✔️:     partition = order[i];
        // RDKit✔️✔️:     oldPart = atoms[partition].index;
        let partition = order[i];
        let old_part = atoms[partition].index;
        // RDKit✔️✔️:     while (count[partition] > 1) {
        while count[partition] > 1 {
            // RDKit✔️✔️:       len = count[partition];
            // RDKit✔️✔️:       offset = atoms[partition].index + len - 1;
            // RDKit✔️✔️:       index = order[offset];
            // RDKit✔️✔️:       atoms[index].index = offset;
            // RDKit✔️✔️:       count[partition] = len - 1;
            // RDKit✔️✔️:       count[index] = 1;
            let len = count[partition];
            let offset = usize::try_from(atoms[partition].index).unwrap_or(usize::MAX) + len - 1;
            let index = order[offset];
            atoms[index].index = i32::try_from(offset).unwrap_or(i32::MAX);
            count[partition] = len - 1;
            count[index] = 1;
            // RDKit✔️✔️:       if (atoms[index].degree < 1) {
            // RDKit✔️✔️:         continue;
            // RDKit✔️✔️:       }
            if atoms[index].degree < 1 {
                continue;
            }
            // RDKit✔️✔️:       for (unsigned j = 0; j < atoms[index].degree; ++j) {
            // RDKit✔️✔️:         unsigned int nbor = atoms[index].nbrIds[j];
            // RDKit✔️✔️:         touchedPartitions[atoms[nbor].index] = 1;
            // RDKit✔️✔️:         changed[nbor] = 1;
            // RDKit✔️✔️:       }
            for nbr in atoms[index].nbr_ids.iter().copied() {
                let partition_idx = usize::try_from(atoms[nbr].index).unwrap_or(usize::MAX);
                if partition_idx < touched_partitions.len() {
                    touched_partitions[partition_idx] = true;
                }
                changed[nbr] = true;
            }
            // RDKit✔️✔️:       for (unsigned int ii = 0; ii < nAtoms; ++ii) {
            // RDKit✔️✔️:         if (touchedPartitions[ii]) {
            // RDKit✔️✔️:           int npart = order[ii];
            // RDKit✔️✔️:           if ((count[npart] > 1) && (next[npart] == -2)) {
            // RDKit✔️✔️:             next[npart] = activeset;
            // RDKit✔️✔️:             activeset = npart;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           touchedPartitions[ii] = 0;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            for ii in 0..n_atoms {
                if touched_partitions[ii] {
                    let npart = order[ii];
                    if count[npart] > 1 && next[npart] == -2 {
                        next[npart] = *active_set;
                        *active_set = isize::try_from(npart).unwrap_or(isize::MAX);
                    }
                    touched_partitions[ii] = false;
                }
            }
            // RDKit✔️✔️:       RefinePartitions(mol, atoms, compar, mode, order, count, activeset, next,
            // RDKit✔️✔️:                        changed, touchedPartitions);
            refine_partitions_for_kekulize(
                view,
                atoms,
                mode,
                compare_mode,
                flags,
                order,
                count,
                active_set,
                next,
                changed,
                touched_partitions,
            );
        }
        // RDKit✔️✔️:     if (atoms[partition].index != oldPart) {
        // RDKit✔️✔️:       i -= 1;
        // RDKit✔️✔️:     }
        if atoms[partition].index != old_part {
            if i > 0 {
                i -= 1;
            }
        } else {
            i += 1;
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION BreakTies
}

fn hanoi_sort_order_for_kekulize(
    order: &mut [usize],
    offset: usize,
    len: usize,
    count: &mut [usize],
    changed: &[bool],
    atoms: &mut [CanonAtom<'_>],
    compare_mode: CanonCompareMode,
    flags: CanonRankFlags,
) {
    // BEGIN RDKIT CPP FUNCTION hanoisort
    // RDKit✔️✔️: template <typename CompareFunc>
    // RDKit✔️✔️: void hanoisort(std::span<int> &base, std::vector<int> &count,
    // RDKit✔️✔️:                std::vector<int> &changed, CompareFunc compar) {
    // RDKit✔️✔️:   std::vector<int> tempVec(base.size());
    let mut temp = vec![0usize; len];
    // RDKit✔️✔️:   if (detail::hanoi(base.data(), base.size(), tempVec.data(), count.data(),
    // RDKit✔️✔️:                     changed.data(), compar)) {
    if hanoi_order_for_kekulize(
        &mut order[offset..offset + len],
        &mut temp,
        count,
        changed,
        atoms,
        compare_mode,
        flags,
    ) {
        // RDKit✔️✔️:     std::copy(tempVec.begin(), tempVec.end(), base.begin());
        order[offset..offset + len].copy_from_slice(&temp);
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION hanoisort
}

fn hanoi_order_for_kekulize(
    base: &mut [usize],
    temp: &mut [usize],
    count: &mut [usize],
    changed: &[bool],
    atoms: &mut [CanonAtom<'_>],
    compare_mode: CanonCompareMode,
    flags: CanonRankFlags,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION detail::hanoi
    // RDKit✔️✔️: template <typename CompareFunc>
    // RDKit✔️✔️: bool hanoi(int *base, int nel, int *temp, int *count, int *changed,
    // RDKit✔️✔️:            CompareFunc compar) {
    // RDKit✔️✔️:   assert(base);
    // RDKit✔️✔️:   assert(temp);
    // RDKit✔️✔️:   assert(count);
    // RDKit✔️✔️:   assert(changed);
    debug_assert_eq!(base.len(), temp.len());
    // RDKit✔️✔️:   int *b1, *b2;
    // RDKit✔️✔️:   int *t1, *t2;
    // RDKit✔️✔️:   int *s1, *s2;
    // RDKit✔️✔️:   int n1, n2;
    // RDKit✔️✔️:   int result;
    // RDKit✔️✔️:   int *ptr;
    let nel = base.len();

    // RDKit✔️✔️:   if (nel == 1) {
    // RDKit✔️✔️:     count[base[0]] = 1;
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   } else if (nel == 2) {
    if nel == 1 {
        count[base[0]] = 1;
        return false;
    } else if nel == 2 {
        // RDKit✔️✔️:     n1 = base[0];
        // RDKit✔️✔️:     n2 = base[1];
        let n1 = base[0];
        let n2 = base[1];
        // RDKit✔️✔️:     int stat =
        // RDKit✔️✔️:         (/*!changed || */ changed[n1] || changed[n2]) ? compar(n1, n2) : 0;
        let stat = if changed[n1] || changed[n2] {
            compare_canon_atoms_for_kekulize(atoms, n1, n2, compare_mode, flags)
        } else {
            Ordering::Equal
        };
        // RDKit✔️✔️:     if (stat == 0) {
        // RDKit✔️✔️:       count[n1] = 2;
        // RDKit✔️✔️:       count[n2] = 0;
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     } else if (stat < 0) {
        // RDKit✔️✔️:       count[n1] = 1;
        // RDKit✔️✔️:       count[n2] = 1;
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     } else /* stat > 0 */ {
        // RDKit✔️✔️:       count[n1] = 1;
        // RDKit✔️✔️:       count[n2] = 1;
        // RDKit✔️✔️:       base[0] = n2; /* temp[0] = n2; */
        // RDKit✔️✔️:       base[1] = n1; /* temp[1] = n1; */
        // RDKit✔️✔️:       return false; /* return True;  */
        // RDKit✔️✔️:     }
        match stat {
            Ordering::Equal => {
                count[n1] = 2;
                count[n2] = 0;
            }
            Ordering::Less => {
                count[n1] = 1;
                count[n2] = 1;
            }
            Ordering::Greater => {
                count[n1] = 1;
                count[n2] = 1;
                base[0] = n2;
                base[1] = n1;
            }
        }
        return false;
    }

    // RDKit✔️✔️:   n1 = nel / 2;
    // RDKit✔️✔️:   n2 = nel - n1;
    let left_len = nel / 2;
    let right_len = nel - left_len;
    // RDKit✔️✔️:   b1 = base;
    // RDKit✔️✔️:   t1 = temp;
    // RDKit✔️✔️:   b2 = base + n1;
    // RDKit✔️✔️:   t2 = temp + n1;
    let (left_in_temp, right_in_temp) = {
        let (base_left, base_right) = base.split_at_mut(left_len);
        let (temp_left, temp_right) = temp.split_at_mut(left_len);

        // RDKit✔️✔️:   if (hanoi(b1, n1, t1, count, changed, compar)) {
        let left_in_temp = hanoi_order_for_kekulize(
            base_left,
            temp_left,
            count,
            changed,
            atoms,
            compare_mode,
            flags,
        );
        // RDKit✔️✔️:     if (hanoi(b2, n2, t2, count, changed, compar)) {
        let right_in_temp = hanoi_order_for_kekulize(
            base_right,
            temp_right,
            count,
            changed,
            atoms,
            compare_mode,
            flags,
        );
        (left_in_temp, right_in_temp)
    };
    // RDKit✔️✔️:     s1 = t1;
    // RDKit✔️✔️:     s1 = b1;
    // RDKit✔️✔️:       s2 = t2;
    // RDKit✔️✔️:       s2 = b2;
    // RDKit✔️✔️:     result = false;
    // RDKit✔️✔️:     ptr = base;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     result = true;
    // RDKit✔️✔️:     ptr = temp;
    // RDKit✔️✔️:   }
    let result = !left_in_temp;
    let base_ptr = base.as_ptr();
    let temp_ptr = temp.as_ptr();
    let left_source_ptr = if left_in_temp { temp_ptr } else { base_ptr };
    let right_source_ptr = if right_in_temp {
        // SAFETY: right half starts at left_len within the same backing buffer.
        unsafe { temp_ptr.add(left_len) }
    } else {
        // SAFETY: right half starts at left_len within the same backing buffer.
        unsafe { base_ptr.add(left_len) }
    };
    let output_ptr = if result {
        temp.as_mut_ptr()
    } else {
        base.as_mut_ptr()
    };
    let mut left_pos = 0usize;
    let mut right_pos = 0usize;
    let mut out_pos = 0usize;
    let mut remaining_left = left_len;
    let mut remaining_right = right_len;

    // RDKit✔️✔️:   while (true) {
    loop {
        let left_atom = unsafe { *left_source_ptr.add(left_pos) };
        let right_atom = unsafe { *right_source_ptr.add(right_pos) };
        // RDKit✔️✔️:     assert(*s1 != *s2);
        debug_assert_ne!(left_atom, right_atom);
        // RDKit✔️✔️:     int stat =
        // RDKit✔️✔️:         (/*!changed || */ changed[*s1] || changed[*s2]) ? compar(*s1, *s2) : 0;
        let stat = if changed[left_atom] || changed[right_atom] {
            compare_canon_atoms_for_kekulize(atoms, left_atom, right_atom, compare_mode, flags)
        } else {
            Ordering::Equal
        };
        // RDKit✔️✔️:     int len1 = count[*s1];
        // RDKit✔️✔️:     int len2 = count[*s2];
        // RDKit✔️✔️:     assert(len1 > 0);
        // RDKit✔️✔️:     assert(len2 > 0);
        let class_left = count[left_atom];
        let class_right = count[right_atom];
        debug_assert!(class_left > 0);
        debug_assert!(class_right > 0);
        // RDKit✔️✔️:     if (stat == 0) {
        if stat == Ordering::Equal {
            // RDKit✔️✔️:       count[*s1] = len1 + len2;
            // RDKit✔️✔️:       count[*s2] = 0;
            count[left_atom] = class_left + class_right;
            count[right_atom] = 0;
            // RDKit✔️✔️:       memmove(ptr, s1, len1 * sizeof(int));
            unsafe {
                std::ptr::copy(
                    left_source_ptr.add(left_pos),
                    output_ptr.add(out_pos),
                    class_left,
                );
            }
            // RDKit✔️✔️:       ptr += len1;
            // RDKit✔️✔️:       n1 -= len1;
            out_pos += class_left;
            remaining_left -= class_left;
            // RDKit✔️✔️:       if (n1 == 0) {
            // RDKit✔️✔️:         if (ptr != s2) {
            // RDKit✔️✔️:           memmove(ptr, s2, n2 * sizeof(int));
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         return result;
            // RDKit✔️✔️:       }
            if remaining_left == 0 {
                unsafe {
                    std::ptr::copy(
                        right_source_ptr.add(right_pos),
                        output_ptr.add(out_pos),
                        remaining_right,
                    );
                }
                return result;
            }
            // RDKit✔️✔️:       s1 += len1;
            left_pos += class_left;
            // RDKit✔️✔️:       memmove(ptr, s2, len2 * sizeof(int));
            unsafe {
                std::ptr::copy(
                    right_source_ptr.add(right_pos),
                    output_ptr.add(out_pos),
                    class_right,
                );
            }
            // RDKit✔️✔️:       ptr += len2;
            // RDKit✔️✔️:       n2 -= len2;
            out_pos += class_right;
            remaining_right -= class_right;
            // RDKit✔️✔️:       if (n2 == 0) {
            // RDKit✔️✔️:         memmove(ptr, s1, n1 * sizeof(int));
            // RDKit✔️✔️:         return result;
            // RDKit✔️✔️:       }
            if remaining_right == 0 {
                unsafe {
                    std::ptr::copy(
                        left_source_ptr.add(left_pos),
                        output_ptr.add(out_pos),
                        remaining_left,
                    );
                }
                return result;
            }
            // RDKit✔️✔️:       s2 += len2;
            right_pos += class_right;
            // RDKit✔️✔️:     } else if (stat < 0 && len1 > 0) {
        } else if stat == Ordering::Less {
            // RDKit✔️✔️:       memmove(ptr, s1, len1 * sizeof(int));
            unsafe {
                std::ptr::copy(
                    left_source_ptr.add(left_pos),
                    output_ptr.add(out_pos),
                    class_left,
                );
            }
            // RDKit✔️✔️:       ptr += len1;
            // RDKit✔️✔️:       n1 -= len1;
            out_pos += class_left;
            remaining_left -= class_left;
            // RDKit✔️✔️:       if (n1 == 0) {
            // RDKit✔️✔️:         if (ptr != s2) {
            // RDKit✔️✔️:           memmove(ptr, s2, n2 * sizeof(int));
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         return result;
            // RDKit✔️✔️:       }
            if remaining_left == 0 {
                unsafe {
                    std::ptr::copy(
                        right_source_ptr.add(right_pos),
                        output_ptr.add(out_pos),
                        remaining_right,
                    );
                }
                return result;
            }
            // RDKit✔️✔️:       s1 += len1;
            left_pos += class_left;
            // RDKit✔️✔️:     } else if (stat > 0 && len2 > 0) /* stat > 0 */ {
        } else {
            // RDKit✔️✔️:       memmove(ptr, s2, len2 * sizeof(int));
            unsafe {
                std::ptr::copy(
                    right_source_ptr.add(right_pos),
                    output_ptr.add(out_pos),
                    class_right,
                );
            }
            // RDKit✔️✔️:       ptr += len2;
            // RDKit✔️✔️:       n2 -= len2;
            out_pos += class_right;
            remaining_right -= class_right;
            // RDKit✔️✔️:       if (n2 == 0) {
            // RDKit✔️✔️:         memmove(ptr, s1, n1 * sizeof(int));
            // RDKit✔️✔️:         return result;
            // RDKit✔️✔️:       }
            if remaining_right == 0 {
                unsafe {
                    std::ptr::copy(
                        left_source_ptr.add(left_pos),
                        output_ptr.add(out_pos),
                        remaining_left,
                    );
                }
                return result;
            }
            // RDKit✔️✔️:       s2 += len2;
            right_pos += class_right;
            // RDKit✔️✔️:     } else {
            // RDKit✔️✔️:       assert(0);
            // RDKit✔️✔️:     }
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION detail::hanoi
}

fn compare_canon_atoms_for_kekulize(
    atoms: &mut [CanonAtom<'_>],
    left: usize,
    right: usize,
    mode: CanonCompareMode,
    flags: CanonRankFlags,
) -> Ordering {
    if matches!(mode, CanonCompareMode::SpecialSymmetry) {
        return compare_special_symmetry_atoms_for_kekulize(atoms, left, right);
    }
    if matches!(mode, CanonCompareMode::SpecialChirality) {
        return compare_special_chirality_atoms_for_kekulize(atoms, left, right);
    }
    if !atom_pair_has_any_in_play_for_kekulize(atoms, left, right) {
        return Ordering::Equal;
    }
    // RDKit✔️✔️:     int v = basecomp(i, j);
    // RDKit✔️✔️:     if (v) {
    // RDKit✔️✔️:       return v;
    // RDKit✔️✔️:     }
    let base_cmp = compare_canon_atom_base_for_kekulize(atoms, left, right, flags);
    if base_cmp != Ordering::Equal {
        return base_cmp;
    }
    // RDKit✔️✔️:     if (df_useNbrs) {
    // RDKit✔️✔️:       if (!dp_atomsInPlay || (*dp_atomsInPlay)[i]) {
    if atoms[left].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, left);
    }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!dp_atomsInPlay || (*dp_atomsInPlay)[j]) {
    if atoms[right].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, right);
    }
    // RDKit✔️✔️:       }
    let ranks = canon_atom_rank_snapshot(atoms);
    // RDKit✔️✔️:       for (unsigned int ii = 0;
    // RDKit✔️✔️:            ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size();
    // RDKit✔️✔️:            ++ii) {
    // RDKit✔️✔️:         int cmp =
    // RDKit✔️✔️:             bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
    // RDKit✔️✔️:         if (cmp) {
    // RDKit✔️✔️:           return cmp;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    for idx in 0..atoms[left].bonds.len().min(atoms[right].bonds.len()) {
        let cmp =
            compare_canon_bond_holder(&atoms[left].bonds[idx], &atoms[right].bonds[idx], &ranks);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    // RDKit✔️✔️:       if (dp_atoms[i].bonds.size() < dp_atoms[j].bonds.size()) {
    // RDKit✔️✔️:         return -1;
    // RDKit✔️✔️:       } else if (dp_atoms[i].bonds.size() > dp_atoms[j].bonds.size()) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return 0;
    atoms[left].bonds.len().cmp(&atoms[right].bonds.len())
}

fn compare_special_chirality_atoms_for_kekulize(
    atoms: &mut [CanonAtom<'_>],
    left: usize,
    right: usize,
) -> Ordering {
    // BEGIN RDKIT CPP FUNCTION SpecialChiralityAtomCompareFunctor::operator()
    // RDKit✔️✔️:   int operator()(int i, int j) const {
    // RDKit✔️✔️:     PRECONDITION(dp_atoms, "no atoms");
    // RDKit✔️✔️:     PRECONDITION(dp_mol, "no molecule");
    // RDKit✔️✔️:     PRECONDITION(i != j, "bad call");
    // RDKit✔️✔️:     if (dp_atomsInPlay && !((*dp_atomsInPlay)[i] || (*dp_atomsInPlay)[j])) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[i]) {
    // RDKit✔️✔️:       updateAtomNeighborIndex(dp_atoms, dp_atoms[i].bonds);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[j]) {
    // RDKit✔️✔️:       updateAtomNeighborIndex(dp_atoms, dp_atoms[j].bonds);
    // RDKit✔️✔️:     }
    if !atom_pair_has_any_in_play_for_kekulize(atoms, left, right) {
        return Ordering::Equal;
    }
    if atoms[left].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, left);
    }
    if atoms[right].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, right);
    }
    // RDKit✔️✔️:     for (unsigned int ii = 0;
    // RDKit✔️✔️:          ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size(); ++ii) {
    // RDKit✔️✔️:       int cmp =
    // RDKit✔️✔️:           bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
    // RDKit✔️✔️:       if (cmp) {
    // RDKit✔️✔️:         return cmp;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    let ranks = canon_atom_rank_snapshot(atoms);
    for idx in 0..atoms[left].bonds.len().min(atoms[right].bonds.len()) {
        let cmp =
            compare_canon_bond_holder(&atoms[left].bonds[idx], &atoms[right].bonds[idx], &ranks);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    // RDKit✔️✔️:     std::vector<std::pair<unsigned int, unsigned int>> swapsi;
    // RDKit✔️✔️:     std::vector<std::pair<unsigned int, unsigned int>> swapsj;
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[i]) {
    // RDKit✔️✔️:       updateAtomNeighborNumSwaps(dp_atoms, dp_atoms[i].bonds, i, swapsi);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[j]) {
    // RDKit✔️✔️:       updateAtomNeighborNumSwaps(dp_atoms, dp_atoms[j].bonds, j, swapsj);
    // RDKit✔️✔️:     }
    let swaps_left = if atoms[left].is_in_play {
        update_atom_neighbor_num_swaps_for_kekulize(atoms, left)
    } else {
        Vec::new()
    };
    let swaps_right = if atoms[right].is_in_play {
        update_atom_neighbor_num_swaps_for_kekulize(atoms, right)
    } else {
        Vec::new()
    };
    // RDKit✔️✔️:     for (unsigned int ii = 0; ii < swapsi.size() && ii < swapsj.size(); ++ii) {
    // RDKit✔️✔️:       int cmp = swapsi[ii].second - swapsj[ii].second;
    // RDKit✔️✔️:       if (cmp) {
    // RDKit✔️✔️:         return cmp;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    for idx in 0..swaps_left.len().min(swaps_right.len()) {
        let cmp = swaps_left[idx].1.cmp(&swaps_right[idx].1);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION SpecialChiralityAtomCompareFunctor::operator()
    Ordering::Equal
}

fn compare_special_symmetry_atoms_for_kekulize(
    atoms: &mut [CanonAtom<'_>],
    left: usize,
    right: usize,
) -> Ordering {
    // BEGIN RDKIT CPP FUNCTION SpecialSymmetryAtomCompareFunctor::operator()
    // RDKit✔️✔️:   int operator()(int i, int j) const {
    // RDKit✔️✔️:     PRECONDITION(dp_atoms, "no atoms");
    // RDKit✔️✔️:     PRECONDITION(dp_mol, "no molecule");
    // RDKit✔️✔️:     PRECONDITION(i != j, "bad call");
    // RDKit✔️✔️:     if (dp_atomsInPlay && !((*dp_atomsInPlay)[i] || (*dp_atomsInPlay)[j])) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (dp_atoms[i].neighborNum < dp_atoms[j].neighborNum) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if (dp_atoms[i].neighborNum > dp_atoms[j].neighborNum) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    if !atom_pair_has_any_in_play_for_kekulize(atoms, left, right) {
        return Ordering::Equal;
    }
    let neighbor_cmp = atoms[left].neighbor_num.cmp(&atoms[right].neighbor_num);
    if neighbor_cmp != Ordering::Equal {
        return neighbor_cmp;
    }
    // RDKit✔️✔️:     if (dp_atoms[i].revistedNeighbors < dp_atoms[j].revistedNeighbors) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if (dp_atoms[i].revistedNeighbors > dp_atoms[j].revistedNeighbors) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    let revisited_cmp = atoms[left]
        .revisted_neighbors
        .cmp(&atoms[right].revisted_neighbors);
    if revisited_cmp != Ordering::Equal {
        return revisited_cmp;
    }
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[i]) {
    // RDKit✔️✔️:       updateAtomNeighborIndex(dp_atoms, dp_atoms[i].bonds);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[j]) {
    // RDKit✔️✔️:       updateAtomNeighborIndex(dp_atoms, dp_atoms[j].bonds);
    // RDKit✔️✔️:     }
    if atoms[left].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, left);
    }
    if atoms[right].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, right);
    }
    let ranks = canon_atom_rank_snapshot(atoms);
    // RDKit✔️✔️:     for (unsigned int ii = 0;
    // RDKit✔️✔️:          ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size(); ++ii) {
    // RDKit✔️✔️:       int cmp =
    // RDKit✔️✔️:           bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
    // RDKit✔️✔️:       if (cmp) {
    // RDKit✔️✔️:         return cmp;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    for idx in 0..atoms[left].bonds.len().min(atoms[right].bonds.len()) {
        let cmp =
            compare_canon_bond_holder(&atoms[left].bonds[idx], &atoms[right].bonds[idx], &ranks);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    // RDKit✔️✔️:     if (dp_atoms[i].bonds.size() < dp_atoms[j].bonds.size()) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if (dp_atoms[i].bonds.size() > dp_atoms[j].bonds.size()) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION SpecialSymmetryAtomCompareFunctor::operator()
    atoms[left].bonds.len().cmp(&atoms[right].bonds.len())
}

fn compare_canon_atom_base_for_kekulize(
    atoms: &[CanonAtom<'_>],
    left: usize,
    right: usize,
    flags: CanonRankFlags,
) -> Ordering {
    // BEGIN RDKIT CPP FUNCTION AtomCompareFunctor::basecomp
    // RDKit✔️✔️:   ivi = dp_atoms[i].index;
    // RDKit✔️✔️:   ivj = dp_atoms[j].index;
    // RDKit✔️✔️:   if (df_useNonStereoRanks) {
    // RDKit✔️✔️:     int rankingNumber_i = 0;
    // RDKit✔️✔️:     int rankingNumber_j = 0;
    // RDKit✔️✔️:     dp_atoms[i].atom->getPropIfPresent(
    // RDKit✔️✔️:         common_properties::_CanonicalRankingNumber, rankingNumber_i);
    // RDKit✔️✔️:     dp_atoms[j].atom->getPropIfPresent(
    // RDKit✔️✔️:         common_properties::_CanonicalRankingNumber, rankingNumber_j);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (df_useAtomMaps || df_useAtomMapsOnDummies) {
    // RDKit✔️✔️:     int molAtomMapNumber_i = 0;
    // RDKit✔️✔️:     int molAtomMapNumber_j = 0;
    // RDKit✔️✔️:     if (df_useAtomMaps ||
    // RDKit✔️✔️:         (df_useAtomMapsOnDummies && dp_atoms[i].atom->getAtomicNum() == 0)) {
    // RDKit✔️✔️:       dp_atoms[i].atom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit✔️✔️:                                          molAtomMapNumber_i);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (df_useAtomMaps ||
    // RDKit✔️✔️:         (df_useAtomMapsOnDummies && dp_atoms[j].atom->getAtomicNum() == 0)) {
    // RDKit✔️✔️:       dp_atoms[j].atom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit✔️✔️:                                          molAtomMapNumber_j);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ivi = dp_atoms[i].degree;
    // RDKit✔️✔️:   ivj = dp_atoms[j].degree;
    // RDKit✔️✔️:   if (dp_atoms[i].p_symbol && dp_atoms[j].p_symbol) {
    // RDKit✔️✔️:     if (*(dp_atoms[i].p_symbol) < *(dp_atoms[j].p_symbol)) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if (*(dp_atoms[i].p_symbol) > *(dp_atoms[j].p_symbol)) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ivi = dp_atoms[i].atom->getAtomicNum();
    // RDKit✔️✔️:   ivj = dp_atoms[j].atom->getAtomicNum();
    // RDKit✔️✔️:   if (df_useIsotopes) {
    // RDKit✔️✔️:     ivi = dp_atoms[i].atom->getIsotope();
    // RDKit✔️✔️:     ivj = dp_atoms[j].atom->getIsotope();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ivi = dp_atoms[i].totalNumHs;
    // RDKit✔️✔️:   ivj = dp_atoms[j].totalNumHs;
    // RDKit✔️✔️:   ivi = dp_atoms[i].atom->getFormalCharge();
    // RDKit✔️✔️:   ivj = dp_atoms[j].atom->getFormalCharge();
    // RDKit✔️✔️:   if (df_useChiralPresence) {
    // RDKit✔️✔️:     ivi =
    // RDKit✔️✔️:         dp_atoms[i].atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED;
    // RDKit✔️✔️:     ivj =
    // RDKit✔️✔️:         dp_atoms[j].atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (df_useChirality) {
    // RDKit✔️✔️:     ivi = dp_atoms[i].whichStereoGroup;
    // RDKit✔️✔️:     ivj = dp_atoms[j].whichStereoGroup;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (df_useChiralityRings) {
    // RDKit✔️✔️:     ivi = getAtomRingNbrCode(i);
    // RDKit✔️✔️:     ivj = getAtomRingNbrCode(j);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION AtomCompareFunctor::basecomp
    let left_atom = &atoms[left];
    let right_atom = &atoms[right];
    let mut cmp = left_atom.index.cmp(&right_atom.index);
    if cmp != Ordering::Equal {
        return cmp;
    }

    if flags.use_non_stereo_ranks {
        cmp = left_atom
            .canonical_ranking_number
            .cmp(&right_atom.canonical_ranking_number);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }

    if flags.use_atom_maps || flags.use_atom_maps_on_dummies {
        let left_map = if flags.use_atom_maps
            || (flags.use_atom_maps_on_dummies && left_atom.atomic_number == 0)
        {
            left_atom.atom_map
        } else {
            0
        };
        let right_map = if flags.use_atom_maps
            || (flags.use_atom_maps_on_dummies && right_atom.atomic_number == 0)
        {
            right_atom.atom_map
        } else {
            0
        };
        cmp = left_map.cmp(&right_map);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }

    cmp = left_atom.degree.cmp(&right_atom.degree);
    if cmp != Ordering::Equal {
        return cmp;
    }

    if let (Some(left_symbol), Some(right_symbol)) = (left_atom.p_symbol, right_atom.p_symbol) {
        return left_symbol.cmp(right_symbol);
    }

    cmp = left_atom.atomic_number.cmp(&right_atom.atomic_number);
    if cmp != Ordering::Equal {
        return cmp;
    }

    if flags.use_isotopes {
        cmp = left_atom.isotope.cmp(&right_atom.isotope);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }

    cmp = left_atom.total_num_hs.cmp(&right_atom.total_num_hs);
    if cmp != Ordering::Equal {
        return cmp;
    }

    // RDKit basecomp stores comparison temporaries as unsigned int.
    // Preserve that C++ conversion behavior here so negative formal charges
    // wrap and order after non-negative charges (e.g. -1 > 0 in this step).
    let left_charge = left_atom.formal_charge as i32 as u32;
    let right_charge = right_atom.formal_charge as i32 as u32;
    cmp = left_charge.cmp(&right_charge);
    if cmp != Ordering::Equal {
        return cmp;
    }

    if flags.use_chiral_presence {
        cmp = chiral_presence_for_kekulize(left_atom.chiral_tag)
            .cmp(&chiral_presence_for_kekulize(right_atom.chiral_tag));
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    if flags.use_chirality {
        // RDKit✔️✔️:     ivi = dp_atoms[i].whichStereoGroup;
        // RDKit✔️✔️:     ivj = dp_atoms[j].whichStereoGroup;
        cmp = compare_stereo_group_state_for_kekulize(atoms, left, right);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    if flags.use_chirality_rings {
        // RDKit✔️✔️:     ivi = getAtomRingNbrCode(i);
        // RDKit✔️✔️:     ivj = getAtomRingNbrCode(j);
        cmp = get_atom_ring_nbr_code_for_kekulize(atoms, left)
            .cmp(&get_atom_ring_nbr_code_for_kekulize(atoms, right));
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    Ordering::Equal
}

fn atom_pair_has_any_in_play_for_kekulize(
    atoms: &[CanonAtom<'_>],
    left: usize,
    right: usize,
) -> bool {
    atoms[left].is_in_play || atoms[right].is_in_play
}

fn compare_stereo_group_state_for_kekulize(
    atoms: &[CanonAtom<'_>],
    left: usize,
    right: usize,
) -> Ordering {
    let left_group = atoms[left].which_stereo_group;
    let right_group = atoms[right].which_stereo_group;
    match (left_group, right_group) {
        (0, 0) => Ordering::Equal,
        (_, 0) => Ordering::Greater,
        (0, _) => Ordering::Less,
        _ => atoms[left]
            .type_of_stereo_group
            .cmp(&atoms[right].type_of_stereo_group)
            .then_with(|| {
                if left_group == right_group {
                    if atoms[left].type_of_stereo_group == CanonStereoGroupType::Absolute {
                        get_chiral_rank_for_kekulize(atoms, left)
                            .cmp(&get_chiral_rank_for_kekulize(atoms, right))
                    } else {
                        Ordering::Equal
                    }
                } else {
                    stereo_group_symmetry_set_for_kekulize(atoms, left_group)
                        .cmp(&stereo_group_symmetry_set_for_kekulize(atoms, right_group))
                }
            }),
    }
    .then_with(|| {
        if left_group == 0 && right_group == 0 {
            chiral_presence_for_kekulize(atoms[left].chiral_tag)
                .cmp(&chiral_presence_for_kekulize(atoms[right].chiral_tag))
                .then_with(|| {
                    if chiral_presence_for_kekulize(atoms[left].chiral_tag)
                        && chiral_presence_for_kekulize(atoms[right].chiral_tag)
                    {
                        get_chiral_rank_for_kekulize(atoms, left)
                            .cmp(&get_chiral_rank_for_kekulize(atoms, right))
                    } else {
                        Ordering::Equal
                    }
                })
        } else {
            Ordering::Equal
        }
    })
}

fn chiral_presence_for_kekulize(chiral_tag: ChiralTag) -> bool {
    chiral_tag != ChiralTag::Unspecified
}

fn get_chiral_rank_for_kekulize(atoms: &[CanonAtom<'_>], atom_idx: usize) -> u32 {
    // BEGIN RDKIT CPP FUNCTION getChiralRank
    // RDKit✔️✔️: unsigned int getChiralRank(const ROMol *dp_mol, canon_atom *dp_atoms,
    // RDKit✔️✔️:                            unsigned int i) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   std::vector<unsigned int> perm;
    // RDKit✔️✔️:   perm.reserve(dp_atoms[i].atom->getDegree());
    let mut res = 0u32;
    let mut perm = Vec::<i32>::with_capacity(atoms[atom_idx].all_nbr_ids.len());
    // RDKit✔️✔️:   for (const auto nbr : dp_mol->atomNeighbors(dp_atoms[i].atom)) {
    // RDKit✔️✔️:     auto rnk = dp_atoms[nbr->getIdx()].index;
    // RDKit✔️✔️:     if (std::find(perm.begin(), perm.end(), rnk) != perm.end()) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       perm.push_back(rnk);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for neighbor in atoms[atom_idx].all_nbr_ids.iter().copied() {
        let rank = atoms[neighbor].index;
        if perm.contains(&rank) {
            break;
        }
        perm.push(rank);
    }
    // RDKit✔️✔️:   if (perm.size() == dp_atoms[i].atom->getDegree()) {
    if perm.len() == atoms[atom_idx].all_nbr_ids.len() {
        // RDKit✔️✔️:     auto ctag = dp_atoms[i].atom->getChiralTag();
        // RDKit✔️✔️:     if (ctag == Atom::ChiralType::CHI_TETRAHEDRAL_CW ||
        // RDKit✔️✔️:         ctag == Atom::ChiralType::CHI_TETRAHEDRAL_CCW) {
        if matches!(
            atoms[atom_idx].chiral_tag,
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            // RDKit✔️✔️:       auto sortedPerm = perm;
            // RDKit✔️✔️:       std::sort(sortedPerm.begin(), sortedPerm.end());
            // RDKit✔️✔️:       auto nswaps = countSwapsToInterconvert(perm, sortedPerm);
            let mut sorted_perm = perm.clone();
            sorted_perm.sort_unstable();
            let swaps = count_swaps_to_interconvert(&perm, sorted_perm);
            // RDKit✔️✔️:       res = ctag == Atom::ChiralType::CHI_TETRAHEDRAL_CW ? 2 : 1;
            res = if atoms[atom_idx].chiral_tag == ChiralTag::TetrahedralCw {
                2
            } else {
                1
            };
            // RDKit✔️✔️:       if (nswaps % 2) {
            // RDKit✔️✔️:         res = res == 2 ? 1 : 2;
            // RDKit✔️✔️:       }
            if swaps % 2 == 1 {
                res = if res == 2 { 1 } else { 2 };
            }
        }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getChiralRank
    res
}

fn stereo_group_symmetry_set_for_kekulize(
    atoms: &[CanonAtom<'_>],
    group_idx: usize,
) -> BTreeSet<i32> {
    atoms
        .iter()
        .filter(|atom| atom.which_stereo_group == group_idx)
        .map(|atom| atom.index)
        .collect()
}

fn get_atom_ring_nbr_code_for_kekulize(atoms: &[CanonAtom<'_>], atom_idx: usize) -> i32 {
    // BEGIN RDKIT CPP FUNCTION AtomCompareFunctor::getAtomRingNbrCode
    // RDKit✔️✔️:   unsigned int getAtomRingNbrCode(unsigned int i) const {
    // RDKit✔️✔️:     if (!dp_atoms[i].hasRingNbr) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    if !atoms[atom_idx].has_ring_nbr {
        return 0;
    }
    // RDKit✔️✔️:     auto nbrs = dp_atoms[i].nbrIds.get();
    // RDKit✔️✔️:     unsigned int code = 0;
    // RDKit✔️✔️:     for (unsigned j = 0; j < dp_atoms[i].degree; ++j) {
    // RDKit✔️✔️:       if (dp_atoms[nbrs[j]].isRingStereoAtom) {
    // RDKit✔️✔️:         code += dp_atoms[nbrs[j]].index * 10000 + 1;  // j;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return code;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION AtomCompareFunctor::getAtomRingNbrCode
    atoms[atom_idx]
        .nbr_ids
        .iter()
        .filter(|&&neighbor| atoms[neighbor].is_ring_stereo_atom)
        .map(|&neighbor| {
            atoms[neighbor]
                .index
                .saturating_mul(10_000)
                .saturating_add(1)
        })
        .sum()
}

fn special_symmetry_rank_refinement_required(
    view: &CanonRankReadView<'_>,
    order: &[usize],
    count: &[usize],
) -> bool {
    let mut ties = false;
    let mut sym_ring_atoms = 0usize;
    let mut ring_atoms = 0usize;
    let mut branching_ring_atom = false;
    for &atom_idx in order.iter().take(view.num_atoms()) {
        let atom = AtomId::new(atom_idx);
        if view.rings.num_atom_rings(atom) > 0 {
            if count[atom_idx] > 2 {
                sym_ring_atoms += count[atom_idx];
            }
            ring_atoms += 1;
            if view.rings.num_atom_rings(atom) > 1 && count[atom_idx] > 1 {
                branching_ring_atom = true;
            }
        }
        if count[atom_idx] == 0 {
            ties = true;
        }
    }
    ties && ring_atoms > 0
        && (sym_ring_atoms as f32) / (ring_atoms as f32) > 0.5
        && branching_ring_atom
}

fn compare_ring_atoms_concerning_num_neighbors_for_kekulize(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
) {
    // BEGIN RDKIT CPP FUNCTION compareRingAtomsConcerningNumNeighbors
    // RDKit✔️✔️: void compareRingAtomsConcerningNumNeighbors(Canon::canon_atom *atoms,
    // RDKit✔️✔️:                                             unsigned int nAtoms,
    // RDKit✔️✔️:                                             const ROMol &mol) {
    // RDKit✔️✔️:   PRECONDITION(atoms, "bad pointer");
    // RDKit✔️✔️:   RingInfo *ringInfo = mol.getRingInfo();
    let n_atoms = view.num_atoms();
    // RDKit✔️✔️:   std::vector<char> visited(nAtoms);
    // RDKit✔️✔️:   std::vector<char> lastLevelNbrs(nAtoms);
    // RDKit✔️✔️:   std::vector<char> currentLevelNbrs(nAtoms);
    // RDKit✔️✔️:   std::vector<int> revisitedNeighbors(nAtoms);
    let mut visited = vec![false; n_atoms];
    let mut last_level_nbrs = vec![false; n_atoms];
    let mut current_level_nbrs = vec![false; n_atoms];
    let mut revisited_neighbors = vec![0i32; n_atoms];
    // RDKit✔️✔️:   std::vector<int> visitedIndices;
    // RDKit✔️✔️:   std::vector<int> lastLevelIndices;
    // RDKit✔️✔️:   std::vector<int> modifiedRevisited;
    let mut visited_indices = Vec::<usize>::new();
    let mut last_level_indices = Vec::<usize>::new();
    let mut modified_revisited = Vec::<usize>::new();
    // RDKit✔️✔️:   for (unsigned idx = 0; idx < nAtoms; ++idx) {
    // RDKit✔️✔️:     const Canon::canon_atom &a = atoms[idx];
    // RDKit✔️✔️:     if (!ringInfo->isInitialized() ||
    // RDKit✔️✔️:         ringInfo->numAtomRings(a.atom->getIdx()) < 1) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    for idx in 0..n_atoms {
        if view.rings.num_atom_rings(AtomId::new(idx)) < 1 {
            continue;
        }
        // RDKit✔️✔️:     std::deque<int> neighbors;
        // RDKit✔️✔️:     neighbors.push_back(idx);
        // RDKit✔️✔️:     unsigned currentRNIdx = 0;
        // RDKit✔️✔️:     atoms[idx].neighborNum.reserve(1000);
        // RDKit✔️✔️:     atoms[idx].revistedNeighbors.assign(1000, 0);
        let mut neighbors = VecDeque::new();
        neighbors.push_back(idx);
        let mut current_rn_idx = 0usize;
        atoms[idx].neighbor_num.clear();
        atoms[idx].neighbor_num.reserve(1000);
        atoms[idx].revisted_neighbors.clear();
        atoms[idx].revisted_neighbors.resize(1000, 0);
        // RDKit✔️✔️:     for (int i : visitedIndices) {
        // RDKit✔️✔️:       visited[i] = 0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     visitedIndices.clear();
        for i in visited_indices.drain(..) {
            visited[i] = false;
        }
        // RDKit✔️✔️:     for (int i : lastLevelIndices) {
        // RDKit✔️✔️:       lastLevelNbrs[i] = 0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     lastLevelIndices.clear();
        for i in last_level_indices.drain(..) {
            last_level_nbrs[i] = false;
        }
        // RDKit✔️✔️:     std::vector<int> nextLevelNbrs;
        let mut next_level_nbrs = Vec::<usize>::new();
        // RDKit✔️✔️:     while (!neighbors.empty()) {
        while !neighbors.is_empty() {
            // RDKit✔️✔️:       unsigned int numLevelNbrs = 0;
            // RDKit✔️✔️:       nextLevelNbrs.resize(0);
            let mut num_level_nbrs = 0i32;
            next_level_nbrs.clear();
            // RDKit✔️✔️:       while (!neighbors.empty()) {
            while let Some(nidx) = neighbors.pop_front() {
                // RDKit✔️✔️:         int nidx = neighbors.front();
                // RDKit✔️✔️:         neighbors.pop_front();
                // RDKit✔️✔️:         const Canon::canon_atom &atom = atoms[nidx];
                // RDKit✔️✔️:         if (!ringInfo->isInitialized() ||
                // RDKit✔️✔️:             ringInfo->numAtomRings(atom.atom->getIdx()) < 1) {
                // RDKit✔️✔️:           continue;
                // RDKit✔️✔️:         }
                // symmetrize_sssr() always returns initialized RingInfo here,
                // so only the ring-membership branch remains.
                if view.rings.num_atom_rings(AtomId::new(nidx)) < 1 {
                    continue;
                }
                // RDKit✔️✔️:         lastLevelNbrs[nidx] = 1;
                // RDKit✔️✔️:         lastLevelIndices.push_back(nidx);
                // RDKit✔️✔️:         visited[nidx] = 1;
                // RDKit✔️✔️:         visitedIndices.push_back(nidx);
                last_level_nbrs[nidx] = true;
                last_level_indices.push(nidx);
                visited[nidx] = true;
                visited_indices.push(nidx);
                // RDKit✔️✔️:         for (unsigned int j = 0; j < atom.degree; j++) {
                // RDKit✔️✔️:           int iidx = atom.nbrIds[j];
                // RDKit✔️✔️:           if (!visited[iidx]) {
                // RDKit✔️✔️:             currentLevelNbrs[iidx] = 1;
                // RDKit✔️✔️:             numLevelNbrs++;
                // RDKit✔️✔️:             visited[iidx] = 1;
                // RDKit✔️✔️:             visitedIndices.push_back(iidx);
                // RDKit✔️✔️:             nextLevelNbrs.push_back(iidx);
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                for iidx in atoms[nidx].nbr_ids.iter().copied() {
                    if !visited[iidx] {
                        current_level_nbrs[iidx] = true;
                        num_level_nbrs += 1;
                        visited[iidx] = true;
                        visited_indices.push(iidx);
                        next_level_nbrs.push(iidx);
                    }
                }
            }
            // RDKit✔️✔️:       for (int i : nextLevelNbrs) {
            // RDKit✔️✔️:         const Canon::canon_atom &natom = atoms[i];
            // RDKit✔️✔️:         for (unsigned int k = 0; k < natom.degree; k++) {
            // RDKit✔️✔️:           int jidx = natom.nbrIds[k];
            // RDKit✔️✔️:           if (currentLevelNbrs[jidx] || lastLevelNbrs[jidx]) {
            // RDKit✔️✔️:             if (revisitedNeighbors[jidx] == 0) {
            // RDKit✔️✔️:               modifiedRevisited.push_back(jidx);
            // RDKit✔️✔️:             }
            // RDKit✔️✔️:             revisitedNeighbors[jidx] += 1;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            for i in next_level_nbrs.iter().copied() {
                for jidx in atoms[i].nbr_ids.iter().copied() {
                    if current_level_nbrs[jidx] || last_level_nbrs[jidx] {
                        if revisited_neighbors[jidx] == 0 {
                            modified_revisited.push(jidx);
                        }
                        revisited_neighbors[jidx] += 1;
                    }
                }
            }
            // RDKit✔️✔️:       for (int i : lastLevelIndices) {
            // RDKit✔️✔️:         lastLevelNbrs[i] = 0;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       lastLevelIndices.clear();
            for i in last_level_indices.drain(..) {
                last_level_nbrs[i] = false;
            }
            // RDKit✔️✔️:       for (int i : nextLevelNbrs) {
            // RDKit✔️✔️:         lastLevelNbrs[i] = 1;
            // RDKit✔️✔️:         lastLevelIndices.push_back(i);
            // RDKit✔️✔️:       }
            for i in next_level_nbrs.iter().copied() {
                last_level_nbrs[i] = true;
                last_level_indices.push(i);
            }
            // RDKit✔️✔️:       for (int i : nextLevelNbrs) {
            // RDKit✔️✔️:         currentLevelNbrs[i] = 0;
            // RDKit✔️✔️:       }
            for i in next_level_nbrs.iter().copied() {
                current_level_nbrs[i] = false;
            }
            // RDKit✔️✔️:       std::vector<int> tmp;
            // RDKit✔️✔️:       tmp.reserve(30);
            // RDKit✔️✔️:       for (int i : modifiedRevisited) {
            // RDKit✔️✔️:         tmp.push_back(revisitedNeighbors[i]);
            // RDKit✔️✔️:       }
            let mut tmp = Vec::with_capacity(30);
            for i in modified_revisited.iter().copied() {
                tmp.push(revisited_neighbors[i]);
            }
            // RDKit✔️✔️:       std::sort(tmp.begin(), tmp.end());
            // RDKit✔️✔️:       tmp.push_back(-1);
            tmp.sort_unstable();
            tmp.push(-1);
            // RDKit✔️✔️:       for (int i : tmp) {
            // RDKit✔️✔️:         if (currentRNIdx >= atoms[idx].revistedNeighbors.size()) {
            // RDKit✔️✔️:           atoms[idx].revistedNeighbors.resize(
            // RDKit✔️✔️:               atoms[idx].revistedNeighbors.size() + 1000);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         atoms[idx].revistedNeighbors[currentRNIdx] = i;
            // RDKit✔️✔️:         currentRNIdx++;
            // RDKit✔️✔️:       }
            for value in tmp {
                if current_rn_idx >= atoms[idx].revisted_neighbors.len() {
                    atoms[idx]
                        .revisted_neighbors
                        .resize(atoms[idx].revisted_neighbors.len() + 1000, 0);
                }
                atoms[idx].revisted_neighbors[current_rn_idx] = value;
                current_rn_idx += 1;
            }
            // RDKit✔️✔️:       for (int i : modifiedRevisited) {
            // RDKit✔️✔️:         revisitedNeighbors[i] = 0;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       modifiedRevisited.clear();
            for i in modified_revisited.drain(..) {
                revisited_neighbors[i] = 0;
            }
            // RDKit✔️✔️:       atoms[idx].neighborNum.push_back(numLevelNbrs);
            // RDKit✔️✔️:       atoms[idx].neighborNum.push_back(-1);
            atoms[idx].neighbor_num.push(num_level_nbrs);
            atoms[idx].neighbor_num.push(-1);
            // RDKit✔️✔️:       neighbors.insert(neighbors.end(), nextLevelNbrs.begin(),
            // RDKit✔️✔️:                        nextLevelNbrs.end());
            // RDKit✔️✔️:     }
            neighbors.extend(next_level_nbrs.iter().copied());
        }
        // RDKit✔️✔️:     atoms[idx].revistedNeighbors.resize(currentRNIdx);
        atoms[idx].revisted_neighbors.truncate(current_rn_idx);
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION compareRingAtomsConcerningNumNeighbors
}

fn update_atom_neighbor_index_for_kekulize(atoms: &mut [CanonAtom<'_>], atom_idx: usize) {
    // BEGIN RDKIT CPP FUNCTION updateAtomNeighborIndex
    // RDKit✔️✔️: void updateAtomNeighborIndex(canon_atom *atoms, std::vector<bondholder> &nbrs) {
    // RDKit✔️✔️:   PRECONDITION(atoms, "bad pointer");
    // RDKit✔️✔️:   for (auto &nbr : nbrs) {
    // RDKit✔️✔️:     unsigned nbrIdx = nbr.nbrIdx;
    // RDKit✔️✔️:     unsigned newSymClass = atoms[nbrIdx].index;
    // RDKit✔️✔️:     nbr.nbrSymClass = newSymClass;
    // RDKit✔️✔️:   }
    let updates = atoms[atom_idx]
        .bonds
        .iter()
        .map(|bond| usize::try_from(atoms[bond.nbr_idx].index).unwrap_or(usize::MAX))
        .collect::<Vec<_>>();
    for (bond, nbr_sym_class) in atoms[atom_idx].bonds.iter_mut().zip(updates) {
        bond.nbr_sym_class = nbr_sym_class;
    }
    // RDKit✔️✔️:   std::sort(nbrs.begin(), nbrs.end(), bondholder::greater);
    // RDKit✔️✔️: }
    let ranks = canon_atom_rank_snapshot(atoms);
    sort_canon_bonds_descending(&mut atoms[atom_idx].bonds, &ranks);
    // END RDKIT CPP FUNCTION updateAtomNeighborIndex
}

fn update_atom_neighbor_num_swaps_for_kekulize(
    atoms: &[CanonAtom<'_>],
    atom_idx: usize,
) -> Vec<(usize, u32)> {
    // BEGIN RDKIT CPP FUNCTION updateAtomNeighborNumSwaps
    // RDKit✔️✔️: void updateAtomNeighborNumSwaps(
    // RDKit✔️✔️:     canon_atom *atoms, std::vector<bondholder> &nbrs, unsigned int atomIdx,
    // RDKit✔️✔️:     std::vector<std::pair<unsigned int, unsigned int>> &result) {
    // RDKit✔️✔️:   bool isRingAtom = queryIsAtomInRing(atoms[atomIdx].atom);
    let is_ring_atom = atoms[atom_idx].is_ring_atom;
    let mut result = Vec::<(usize, u32)>::new();
    // RDKit✔️✔️:   for (auto &nbr : nbrs) {
    for nbr in &atoms[atom_idx].bonds {
        // RDKit✔️✔️:     unsigned nbrIdx = nbr.nbrIdx;
        let nbr_idx = nbr.nbr_idx;
        // RDKit✔️✔️:     std::list<unsigned int> neighborsSeen;
        // RDKit✔️✔️:     bool tooManySimilarNbrs = false;
        let mut neighbors_seen = Vec::<i32>::new();
        let mut too_many_similar_neighbors = false;
        // RDKit✔️✔️:     if (isRingAtom && atoms[nbrIdx].atom->getChiralTag() != 0) {
        if is_ring_atom && atoms[nbr_idx].chiral_tag != ChiralTag::Unspecified {
            // RDKit✔️✔️:       std::vector<int> ref, probe;
            let mut reference = Vec::<usize>::new();
            // RDKit✔️✔️:       for (unsigned i = 0; i < atoms[nbrIdx].degree; ++i) {
            for nbr_nbr_id in atoms[nbr_idx].nbr_ids.iter().copied() {
                // RDKit✔️✔️:         auto nbrNbrId =
                // RDKit✔️✔️:             atoms[nbrIdx].nbrIds[i];
                // RDKit✔️✔️:         ref.push_back(nbrNbrId);
                reference.push(nbr_nbr_id);
                // RDKit✔️✔️:         if ((int)atomIdx != nbrNbrId) {
                if atom_idx != nbr_nbr_id {
                    // RDKit✔️✔️:           if ((std::find(neighborsSeen.begin(), neighborsSeen.end(),
                    // RDKit✔️✔️:                          atoms[nbrNbrId].index) != neighborsSeen.end())) {
                    // RDKit✔️✔️:             tooManySimilarNbrs = true;
                    // RDKit✔️✔️:           } else {
                    // RDKit✔️✔️:             neighborsSeen.push_back(atoms[nbrNbrId].index);
                    // RDKit✔️✔️:           }
                    let neighbor_rank = atoms[nbr_nbr_id].index;
                    if neighbors_seen.contains(&neighbor_rank) {
                        too_many_similar_neighbors = true;
                    } else {
                        neighbors_seen.push(neighbor_rank);
                    }
                }
            }
            // RDKit✔️✔️:       probe.push_back(atomIdx);
            let mut probe = vec![atom_idx];
            // RDKit✔️✔️:       for (auto &bond : atoms[nbrIdx].bonds) {
            // RDKit✔️✔️:         if (bond.nbrIdx != atomIdx) {
            // RDKit✔️✔️:           probe.push_back(bond.nbrIdx);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            for bond in &atoms[nbr_idx].bonds {
                if bond.nbr_idx != atom_idx {
                    probe.push(bond.nbr_idx);
                }
            }
            // RDKit✔️✔️:       if (tooManySimilarNbrs) {
            // RDKit✔️✔️:         result.emplace_back(nbr.nbrSymClass, 0);
            // RDKit✔️✔️:       } else {
            if too_many_similar_neighbors {
                result.push((nbr.nbr_sym_class, 0));
            } else {
                // RDKit✔️✔️:         int nSwaps = static_cast<int>(countSwapsToInterconvert(ref, probe));
                let swaps = count_swaps_to_interconvert(&reference, probe);
                // RDKit✔️✔️:         if (atoms[nbrIdx].atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW) {
                // RDKit✔️✔️:           if (nSwaps % 2) {
                // RDKit✔️✔️:             result.emplace_back(nbr.nbrSymClass, 2);
                // RDKit✔️✔️:           } else {
                // RDKit✔️✔️:             result.emplace_back(nbr.nbrSymClass, 1);
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         } else if (atoms[nbrIdx].atom->getChiralTag() ==
                // RDKit✔️✔️:                    Atom::CHI_TETRAHEDRAL_CCW) {
                // RDKit✔️✔️:           if (nSwaps % 2) {
                // RDKit✔️✔️:             result.emplace_back(nbr.nbrSymClass, 1);
                // RDKit✔️✔️:           } else {
                // RDKit✔️✔️:             result.emplace_back(nbr.nbrSymClass, 2);
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                match atoms[nbr_idx].chiral_tag {
                    ChiralTag::TetrahedralCw => {
                        result.push((nbr.nbr_sym_class, if swaps % 2 == 1 { 2 } else { 1 }));
                    }
                    ChiralTag::TetrahedralCcw => {
                        result.push((nbr.nbr_sym_class, if swaps % 2 == 1 { 1 } else { 2 }));
                    }
                    _ => {}
                }
            }
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     } else {
        } else {
            // RDKit✔️✔️:       result.emplace_back(nbr.nbrSymClass, 0);
            result.push((nbr.nbr_sym_class, 0));
        }
        // RDKit✔️✔️:     }
    }
    // RDKit✔️✔️:   sort(result.begin(), result.end());
    // RDKit✔️✔️: }
    result.sort_unstable();
    // END RDKIT CPP FUNCTION updateAtomNeighborNumSwaps
    result
}

fn canon_atom_rank_snapshot(atoms: &[CanonAtom<'_>]) -> Vec<i32> {
    atoms.iter().map(|atom| atom.index).collect()
}

fn sort_canon_bonds_descending(bonds: &mut [CanonBondHolder<'_>], atom_ranks: &[i32]) {
    bonds.sort_by(|left, right| compare_canon_bond_holder(right, left, atom_ranks));
}

fn compare_canon_bond_holder(
    left: &CanonBondHolder<'_>,
    right: &CanonBondHolder<'_>,
    atom_ranks: &[i32],
) -> Ordering {
    // BEGIN RDKIT CPP FUNCTION bondholder::compare
    // RDKit✔️✔️: static int compare(const bondholder &x, const bondholder &y,
    // RDKit✔️✔️:                    unsigned int div = 1) {
    // RDKit✔️✔️:   if (x.p_symbol && y.p_symbol) {
    // RDKit✔️✔️:     if ((*x.p_symbol) < (*y.p_symbol)) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if ((*x.p_symbol) > (*y.p_symbol)) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (x.bondType < y.bondType) {
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   } else if (x.bondType > y.bondType) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    if let (Some(left_symbol), Some(right_symbol)) = (left.p_symbol, right.p_symbol) {
        let symbol_cmp = left_symbol.cmp(right_symbol);
        if symbol_cmp != Ordering::Equal {
            return symbol_cmp;
        }
    }
    let base_cmp = rdkit_bond_order_rank(left.bond_type)
        .cmp(&rdkit_bond_order_rank(right.bond_type))
        // RDKit✔️✔️:   if (x.bondStereo < y.bondStereo) {
        // RDKit✔️✔️:     return -1;
        // RDKit✔️✔️:   } else if (x.bondStereo > y.bondStereo) {
        // RDKit✔️✔️:     return 1;
        // RDKit✔️✔️:   }
        .then_with(|| left.bond_stereo.cmp(&right.bond_stereo))
        // RDKit✔️✔️:   auto scdiv = x.nbrSymClass / div - y.nbrSymClass / div;
        // RDKit✔️✔️:   if (scdiv) {
        // RDKit✔️✔️:     return scdiv;
        // RDKit✔️✔️:   }
        .then_with(|| left.nbr_sym_class.cmp(&right.nbr_sym_class));
    if base_cmp != Ordering::Equal {
        return base_cmp;
    }
    // RDKit✔️✔️:   if (x.bondStereo && y.bondStereo) {
    // RDKit✔️✔️:     auto cs = x.compareStereo(y);
    // RDKit✔️✔️:     if (cs) {
    // RDKit✔️✔️:       return cs;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if left.bond_stereo != 0 && right.bond_stereo != 0 {
        let stereo_cmp = compare_canon_bond_stereo(left, right, atom_ranks);
        if stereo_cmp != Ordering::Equal {
            return stereo_cmp;
        }
    }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION bondholder::compare
    Ordering::Equal
}

fn compare_canon_bond_stereo(
    left: &CanonBondHolder<'_>,
    right: &CanonBondHolder<'_>,
    atom_ranks: &[i32],
) -> Ordering {
    // BEGIN RDKIT CPP FUNCTION bondholder::compareStereo
    // RDKit✔️✔️: int bondholder::compareStereo(const bondholder &o) const {
    // RDKit✔️✔️:   auto st1 = stype;
    // RDKit✔️✔️:   auto st2 = o.stype;
    let mut st1 = left.stype;
    let mut st2 = right.stype;
    // RDKit✔️✔️:   if (st1 == Bond::BondStereo::STEREONONE) {
    // RDKit✔️✔️:     if (st2 == Bond::BondStereo::STEREONONE) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if st1 == BondStereo::None {
        return if st2 == BondStereo::None {
            Ordering::Equal
        } else {
            Ordering::Less
        };
    }
    // RDKit✔️✔️:   if (st2 == Bond::BondStereo::STEREONONE) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    if st2 == BondStereo::None {
        return Ordering::Greater;
    }
    // RDKit✔️✔️:   if (st1 == Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:     if (st2 == Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if st1 == BondStereo::Any {
        return if st2 == BondStereo::Any {
            Ordering::Equal
        } else {
            Ordering::Less
        };
    }
    // RDKit✔️✔️:   if (st2 == Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    if st2 == BondStereo::Any {
        return Ordering::Greater;
    }
    // RDKit✔️✔️:   // we have some kind of specified stereo on both bonds, work is required
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // if both have absolute stereo labels we can compare them directly
    // RDKit✔️✔️:   if ((st1 == Bond::BondStereo::STEREOE || st1 == Bond::BondStereo::STEREOZ) &&
    // RDKit✔️✔️:       (st2 == Bond::BondStereo::STEREOE || st2 == Bond::BondStereo::STEREOZ)) {
    // RDKit✔️✔️:     if (st1 < st2) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if (st1 > st2) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    if matches!(st1, BondStereo::E | BondStereo::Z) && matches!(st2, BondStereo::E | BondStereo::Z)
    {
        return rdkit_bond_stereo_rank(st1).cmp(&rdkit_bond_stereo_rank(st2));
    }
    // RDKit✔️✔️:   // check to see if we need to flip the controlling atoms due to atom ranks
    // RDKit✔️✔️:   flipIfNeeded(st1, controllingAtoms);
    // RDKit✔️✔️:   flipIfNeeded(st2, o.controllingAtoms);
    st1 = flip_bond_stereo_if_needed(st1, left.controlling_atoms, atom_ranks);
    st2 = flip_bond_stereo_if_needed(st2, right.controlling_atoms, atom_ranks);
    // RDKit✔️✔️:   if (st1 < st2) {
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   } else if (st1 > st2) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION bondholder::compareStereo
    rdkit_bond_stereo_rank(st1).cmp(&rdkit_bond_stereo_rank(st2))
}

fn flip_bond_stereo_if_needed(
    mut stereo: BondStereo,
    controlling_atoms: [Option<usize>; 4],
    atom_ranks: &[i32],
) -> BondStereo {
    // BEGIN RDKIT CPP FUNCTION flipIfNeeded
    // RDKit✔️✔️: void flipIfNeeded(Bond::BondStereo &st1,
    // RDKit✔️✔️:                   const canon_atom *const *controllingAtoms) {
    // RDKit✔️✔️:   CHECK_INVARIANT(controllingAtoms[0], "missing controlling atom");
    // RDKit✔️✔️:   CHECK_INVARIANT(controllingAtoms[2], "missing controlling atom");
    let controlling_0 = controlling_atoms[0].expect("missing controlling atom");
    let controlling_2 = controlling_atoms[2].expect("missing controlling atom");
    // RDKit✔️✔️:   bool flip = false;
    let mut flip = false;
    // RDKit✔️✔️:   if (controllingAtoms[1] &&
    // RDKit✔️✔️:       controllingAtoms[1]->index > controllingAtoms[0]->index) {
    // RDKit✔️✔️:     flip = !flip;
    // RDKit✔️✔️:   }
    if controlling_atoms[1].is_some_and(|idx| atom_ranks[idx] > atom_ranks[controlling_0]) {
        flip = !flip;
    }
    // RDKit✔️✔️:   if (controllingAtoms[3] &&
    // RDKit✔️✔️:       controllingAtoms[3]->index > controllingAtoms[2]->index) {
    // RDKit✔️✔️:     flip = !flip;
    // RDKit✔️✔️:   }
    if controlling_atoms[3].is_some_and(|idx| atom_ranks[idx] > atom_ranks[controlling_2]) {
        flip = !flip;
    }
    // RDKit✔️✔️:   if (flip) {
    // RDKit✔️✔️:     if (st1 == Bond::BondStereo::STEREOCIS) {
    // RDKit✔️✔️:       st1 = Bond::BondStereo::STEREOTRANS;
    // RDKit✔️✔️:     } else if (st1 == Bond::BondStereo::STEREOTRANS) {
    // RDKit✔️✔️:       st1 = Bond::BondStereo::STEREOCIS;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if flip {
        stereo = match stereo {
            BondStereo::Cis => BondStereo::Trans,
            BondStereo::Trans => BondStereo::Cis,
            other => other,
        };
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION flipIfNeeded
    stereo
}

fn rdkit_bond_order_rank(order: BondOrder) -> u8 {
    match order {
        BondOrder::Null | BondOrder::Unspecified => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Quintuple => 5,
        BondOrder::Hextuple => 6,
        BondOrder::OneAndHalf => 7,
        BondOrder::TwoAndHalf => 8,
        BondOrder::ThreeAndHalf => 9,
        BondOrder::FourAndHalf => 10,
        BondOrder::FiveAndHalf => 11,
        BondOrder::Aromatic => 12,
        BondOrder::Ionic => 13,
        BondOrder::Hydrogen => 14,
        BondOrder::ThreeCenter => 15,
        BondOrder::DativeOne => 16,
        BondOrder::Dative => 17,
        BondOrder::DativeLeft => 18,
        BondOrder::DativeRight => 19,
        BondOrder::Other => 20,
        BondOrder::Zero => 21,
    }
}

fn rdkit_bond_stereo_rank(stereo: BondStereo) -> u8 {
    match stereo {
        BondStereo::None => 0,
        BondStereo::Any => 1,
        BondStereo::Z => 2,
        BondStereo::E => 3,
        BondStereo::Cis => 4,
        BondStereo::Trans => 5,
        BondStereo::AtropCw => 6,
        BondStereo::AtropCcw => 7,
    }
}

#[cfg(test)]
mod tests;
