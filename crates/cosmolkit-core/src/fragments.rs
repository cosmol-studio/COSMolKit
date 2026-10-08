//! Source-backed private fragment selection and construction.

use std::collections::BTreeMap;

use crate::{SanitizeError, SanitizeParams, sanitize_topology};
use cosmolkit_model::{
    Atom, AtomId, AtomMapping, Bond, BondId, BondMapping, BondStereo, ChiralTag, Conformer2D,
    Conformer3D, CoordinateBlock, CoordinateDimension, CoordinateValidationError,
    MappingValidationError, MoleculeProperties, PropertyText, StereoGroup, SubstanceGroup,
    SubstanceGroupId, TopologyBlock, TopologyEditError, TopologyMapping, TopologyValidationError,
};

/// Actual optional source annotations for detached fragment copying.
/// `None` means that independent source capability is not modeled by this input.
#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct FragmentSourceMetadataView<'a> {
    pub rings: Option<&'a crate::RingInfo>,
    pub atom_bookmarks: Option<&'a BTreeMap<i32, Vec<AtomId>>>,
    pub bond_bookmarks: Option<&'a BTreeMap<i32, Vec<BondId>>>,
}

impl FragmentSourceMetadataView<'_> {
    fn unmodeled() -> Self {
        Self {
            rings: None,
            atom_bookmarks: None,
            bond_bookmarks: None,
        }
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FragmentSourceMetadata {
    pub rings: Option<crate::RingInfo>,
    pub atom_bookmarks: Option<BTreeMap<i32, Vec<AtomId>>>,
    pub bond_bookmarks: Option<BTreeMap<i32, Vec<BondId>>>,
}

impl FragmentSourceMetadata {
    fn unmodeled() -> Self {
        Self {
            rings: None,
            atom_bookmarks: None,
            bond_bookmarks: None,
        }
    }
}

const MASK_WORD_BITS: usize = usize::BITS as usize;

#[derive(Debug, Clone, Copy)]
enum FragmentCoordinates3D<'a> {
    Packed(&'a [[f64; 3]]),
    KernelRows(&'a [&'a [f64]]),
}

impl FragmentCoordinates3D<'_> {
    fn len(self) -> usize {
        match self {
            Self::Packed(rows) => rows.len(),
            Self::KernelRows(rows) => rows.len(),
        }
    }

    fn get(self, row: usize) -> Option<[f64; 3]> {
        match self {
            Self::Packed(rows) => rows.get(row).copied(),
            Self::KernelRows(rows) => {
                let coordinates = rows.get(row)?;
                Some([
                    *coordinates.first()?,
                    *coordinates.get(1)?,
                    *coordinates.get(2)?,
                ])
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct FragmentConformer2D<'a> {
    id: usize,
    coordinates: &'a [[f64; 2]],
    props: &'a BTreeMap<PropertyText, PropertyText>,
}

#[derive(Debug, Clone, Copy)]
struct FragmentConformer3D<'a> {
    id: usize,
    coordinates: FragmentCoordinates3D<'a>,
    is_3d: bool,
    props: &'a BTreeMap<PropertyText, PropertyText>,
}

/// Borrowed conformer rows and metadata used by detached fragment copying.
///
/// The descriptor can substitute the selected force-field coordinate rows
/// while retaining that source conformer's ID, 3D flag, and properties.
#[derive(Debug, Clone)]
pub struct FragmentCoordinateView<'a> {
    conformers_2d: Vec<FragmentConformer2D<'a>>,
    conformers_3d: Vec<FragmentConformer3D<'a>>,
    source_coordinate_dim: Option<CoordinateDimension>,
    source_conformer_order: Option<Vec<CoordinateDimension>>,
}

/// Invalid selected-coordinate metadata for a borrowed fragment view.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum FragmentCoordinateViewError {
    #[error("selected 3D conformer {id} is absent")]
    SelectedConformerMissing { id: usize },
    #[error("selected 3D conformer {id} occurs more than once")]
    SelectedConformerDuplicated { id: usize },
    #[error("selected 3D conformer {id} has {actual} kernel coordinate rows, expected {expected}")]
    KernelCoordinateRowCount {
        id: usize,
        expected: usize,
        actual: usize,
    },
    #[error("selected 3D conformer {id} row {row} has {actual} values, expected 3")]
    KernelCoordinateRowWidth {
        id: usize,
        row: usize,
        actual: usize,
    },
}

impl<'a> FragmentCoordinateView<'a> {
    /// Drop the 3D descriptors from this detached borrowed copy view.
    /// Coordinate rows are never cloned or mutated.
    #[must_use]
    pub fn without_3d_conformers(mut self) -> Self {
        self.conformers_3d.clear();
        if let Some(order) = &mut self.source_conformer_order {
            order.retain(|d| *d != CoordinateDimension::ThreeD);
        }
        self
    }

    /// Borrow every conformer row, ID, flag, and property map in source order.
    ///
    /// Complexity: allocates two O(C) descriptor vectors for C conformers;
    /// coordinate rows and property maps remain borrowed without a data clone.
    #[must_use]
    pub fn from_coordinate_block(source: &'a CoordinateBlock) -> Self {
        Self {
            conformers_2d: source
                .conformers_2d
                .iter()
                .map(|conformer| FragmentConformer2D {
                    id: conformer.id(),
                    coordinates: conformer.coordinates(),
                    props: conformer.props(),
                })
                .collect(),
            conformers_3d: source
                .conformers_3d
                .iter()
                .map(|conformer| FragmentConformer3D {
                    id: conformer.id(),
                    coordinates: FragmentCoordinates3D::Packed(conformer.coordinates()),
                    is_3d: conformer.is_3d(),
                    props: conformer.props(),
                })
                .collect(),
            source_coordinate_dim: source.source_coordinate_dim,
            source_conformer_order: source.source_conformer_order.clone(),
        }
    }

    /// Borrow a selected kernel conformer between disjoint source slices.
    ///
    /// The selected conformer's property map is passed separately so callers
    /// can retain its source ID/properties while its position rows are borrowed
    /// mutably by the force-field kernel. The other source conformers and all
    /// 2D conformers remain borrowed in their original order.
    #[allow(clippy::too_many_arguments)]
    pub fn from_split_conformers(
        conformers_2d: &'a [Conformer2D],
        conformers_3d_before: &'a [Conformer3D],
        selected_id: usize,
        selected_is_3d: bool,
        selected_props: &'a BTreeMap<PropertyText, PropertyText>,
        selected_kernel_rows: &'a [&'a [f64]],
        conformers_3d_after: &'a [Conformer3D],
        source_coordinate_dim: Option<CoordinateDimension>,
        source_conformer_order: Option<&[CoordinateDimension]>,
    ) -> Result<Self, FragmentCoordinateViewError> {
        for (row, coordinates) in selected_kernel_rows.iter().enumerate() {
            if coordinates.len() != 3 {
                return Err(FragmentCoordinateViewError::KernelCoordinateRowWidth {
                    id: selected_id,
                    row,
                    actual: coordinates.len(),
                });
            }
        }

        let conformers_2d = conformers_2d
            .iter()
            .map(|conformer| FragmentConformer2D {
                id: conformer.id(),
                coordinates: conformer.coordinates(),
                props: conformer.props(),
            })
            .collect();
        let mut conformers_3d =
            Vec::with_capacity(conformers_3d_before.len() + conformers_3d_after.len() + 1);
        conformers_3d.extend(
            conformers_3d_before
                .iter()
                .map(|conformer| FragmentConformer3D {
                    id: conformer.id(),
                    coordinates: FragmentCoordinates3D::Packed(conformer.coordinates()),
                    is_3d: conformer.is_3d(),
                    props: conformer.props(),
                }),
        );
        conformers_3d.push(FragmentConformer3D {
            id: selected_id,
            coordinates: FragmentCoordinates3D::KernelRows(selected_kernel_rows),
            is_3d: selected_is_3d,
            props: selected_props,
        });
        conformers_3d.extend(
            conformers_3d_after
                .iter()
                .map(|conformer| FragmentConformer3D {
                    id: conformer.id(),
                    coordinates: FragmentCoordinates3D::Packed(conformer.coordinates()),
                    is_3d: conformer.is_3d(),
                    props: conformer.props(),
                }),
        );

        Ok(Self {
            conformers_2d,
            conformers_3d,
            source_coordinate_dim,
            source_conformer_order: source_conformer_order.map(<[CoordinateDimension]>::to_vec),
        })
    }

    /// Borrow the selected kernel's 3D rows in place of one canonical 3D row
    /// set while retaining all source conformer metadata and other dimensions.
    pub fn with_selected_3d_kernel_rows(
        source: &'a CoordinateBlock,
        selected_id: usize,
        kernel_rows: &'a [&'a [f64]],
    ) -> Result<Self, FragmentCoordinateViewError> {
        for (row, coordinates) in kernel_rows.iter().enumerate() {
            if coordinates.len() != 3 {
                return Err(FragmentCoordinateViewError::KernelCoordinateRowWidth {
                    id: selected_id,
                    row,
                    actual: coordinates.len(),
                });
            }
        }

        let mut view = Self::from_coordinate_block(source);
        let mut selected = None;
        for (index, conformer) in view.conformers_3d.iter_mut().enumerate() {
            if conformer.id == selected_id {
                if selected.replace(index).is_some() {
                    return Err(FragmentCoordinateViewError::SelectedConformerDuplicated {
                        id: selected_id,
                    });
                }
                if kernel_rows.len() != conformer.coordinates.len() {
                    return Err(FragmentCoordinateViewError::KernelCoordinateRowCount {
                        id: selected_id,
                        expected: conformer.coordinates.len(),
                        actual: kernel_rows.len(),
                    });
                }
                conformer.coordinates = FragmentCoordinates3D::KernelRows(kernel_rows);
            }
        }
        if selected.is_none() {
            return Err(FragmentCoordinateViewError::SelectedConformerMissing { id: selected_id });
        }
        Ok(view)
    }

    fn to_coordinate_block(&self) -> CoordinateBlock {
        let conformers_2d = self
            .conformers_2d
            .iter()
            .map(|conformer| {
                let mut copied = Conformer2D::new(conformer.id, conformer.coordinates.to_vec());
                for (key, value) in conformer.props {
                    copied = copied.with_prop(key.clone(), value.clone());
                }
                copied
            })
            .collect();
        let conformers_3d = self
            .conformers_3d
            .iter()
            .map(|conformer| {
                let mut copied = Conformer3D::new(
                    conformer.id,
                    (0..conformer.coordinates.len())
                        .map(|row| {
                            conformer
                                .coordinates
                                .get(row)
                                .expect("borrowed fragment coordinate row has width three")
                        })
                        .collect(),
                    conformer.is_3d,
                );
                for (key, value) in conformer.props {
                    copied = copied.with_prop(key.clone(), value.clone());
                }
                copied
            })
            .collect();
        CoordinateBlock {
            conformers_2d,
            conformers_3d,
            source_coordinate_dim: self.source_coordinate_dim,
            source_conformer_order: self.source_conformer_order.clone(),
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
struct SelectionMask {
    bit_count: usize,
    words: Vec<usize>,
}

impl SelectionMask {
    fn clear(&mut self) {
        self.bit_count = 0;
        self.words.clear();
    }

    fn resize(&mut self, bit_count: usize) {
        self.bit_count = bit_count;
        let word_count = bit_count / MASK_WORD_BITS + usize::from(bit_count % MASK_WORD_BITS != 0);
        self.words.resize(word_count, 0);
    }

    fn set(&mut self, index: usize) {
        debug_assert!(index < self.bit_count);
        self.words[index / MASK_WORD_BITS] |= 1 << (index % MASK_WORD_BITS);
    }

    fn contains(&self, index: usize) -> bool {
        index < self.bit_count
            && self.words[index / MASK_WORD_BITS] & (1 << (index % MASK_WORD_BITS)) != 0
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
struct FragmentSubsetInfo {
    selected_atoms: SelectionMask,
    selected_bonds: SelectionMask,
    atom_mapping: BTreeMap<AtomId, AtomId>,
    bond_mapping: BTreeMap<BondId, BondId>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
enum SelectedBondCopyError {
    #[error("copyMolSubset: subset bonds contain atoms not contained in subset atoms")]
    MissingEndpointMapping { bond: BondId, atom: AtomId },
}

fn get_subset_info_for_atom_path(
    topology: &TopologyBlock,
    path: &[AtomId],
    selection_info: &mut FragmentSubsetInfo,
) {
    // This is the caller's BONDS_BETWEEN_ATOMS specialization. The other
    // source branch remains separate from this atom-path helper.
    // BEGIN RDKIT CPP FUNCTION getSubsetInfo
    // RDKit✔️✔️: static void getSubsetInfo(SubsetInfo &selection_info, const RDKit::ROMol &mol,
    // RDKit✔️✔️:                           const std::vector<unsigned int> &path,
    // RDKit✔️✔️:                           const SubsetOptions &options) {
    // RDKit✔️✔️:   const auto num_atoms = mol.getNumAtoms();
    // RDKit✔️✔️:   const auto num_bonds = mol.getNumBonds();
    // RDKit✔️✔️:   selection_info.selectedAtoms.clear();
    // RDKit✔️✔️:   selection_info.selectedAtoms.resize(num_atoms);
    // RDKit✔️✔️:   selection_info.selectedBonds.clear();
    // RDKit✔️✔️:   selection_info.selectedBonds.resize(num_bonds);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   selection_info.atomMapping.clear();
    // RDKit✔️✔️:   selection_info.bondMapping.clear();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto &[selectedAtoms, selectedBonds, atomMapping, bondMapping] =
    // RDKit✔️✔️:       selection_info;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (options.method == SubsetMethod::BONDS_BETWEEN_ATOMS) {
    // RDKit✔️✔️:     for (const auto &atom_idx : path) {
    // RDKit✔️✔️:       if (atom_idx < num_atoms) {
    // RDKit✔️✔️:         selectedAtoms.set(atom_idx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto &bond : mol.bonds()) {
    // RDKit✔️✔️:       if (selectedAtoms[bond->getBeginAtomIdx()] &&
    // RDKit✔️✔️:           selectedAtoms[bond->getEndAtomIdx()]) {
    // RDKit✔️✔️:         selectedBonds.set(bond->getIdx());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit❌❌:   } else if (options.method == SubsetMethod::BONDS) {
    // RDKit❌❌:     for (const auto &bond_idx : path) {
    // RDKit❌❌:       if (bond_idx < num_bonds) {
    // RDKit❌❌:         selectedBonds.set(bond_idx);
    // RDKit❌❌:         const auto &bnd = mol.getBondWithIdx(bond_idx);
    // RDKit❌❌:         selectedAtoms.set(bnd->getBeginAtomIdx());
    // RDKit❌❌:         selectedAtoms.set(bnd->getEndAtomIdx());
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getSubsetInfo

    let num_atoms = topology.atoms.len();
    let num_bonds = topology.bonds.len();
    selection_info.selected_atoms.clear();
    selection_info.selected_atoms.resize(num_atoms);
    selection_info.selected_bonds.clear();
    selection_info.selected_bonds.resize(num_bonds);
    selection_info.atom_mapping.clear();
    selection_info.bond_mapping.clear();

    for atom in path {
        if atom.index() < num_atoms {
            selection_info.selected_atoms.set(atom.index());
        }
    }
    for (bond_index, bond) in topology.bonds.iter().enumerate() {
        if selection_info.selected_atoms.contains(bond.begin().index())
            && selection_info.selected_atoms.contains(bond.end().index())
        {
            selection_info.selected_bonds.set(bond_index);
        }
    }
}

fn copy_selected_atoms(
    reference: &TopologyBlock,
    selection_info: &mut FragmentSubsetInfo,
) -> Result<Vec<Atom>, cosmolkit_model::AtomPropertyError> {
    // The UFF-FRAG specialization fixes copyAsQuery=false and uses detached
    // concrete atoms; query lowering and live-molecule ownership are outside it.
    // BEGIN RDKIT CPP FUNCTION copySelectedAtomsAndBonds concrete atom loop
    // RDKit✔️✔️:   for (const auto &ref_atom : reference_mol.atoms()) {
    // RDKit✔️✔️:     if (!selectedAtoms[ref_atom->getIdx()]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     std::unique_ptr<Atom> extracted_atom{
    // RDKit✔️✔️:         options.copyAsQuery ? new QueryAtom(*ref_atom) : ref_atom->copy()};
    // RDKit✔️✔️:     extracted_atom->clearComputedProps();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     constexpr bool updateLabel = false;
    // RDKit✔️✔️:     constexpr bool takeOwnership = true;
    // RDKit✔️✔️:     atomMapping[ref_atom->getIdx()] = extracted_mol.addAtom(
    // RDKit✔️✔️:         extracted_atom.release(), updateLabel, takeOwnership);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION copySelectedAtomsAndBonds concrete atom loop
    // Behavior: selected source rows are scanned in order, copied as concrete
    // atoms, cleared unconditionally, appended, and then entered in the map.
    // Complexity: O(V + A log A + copied property bytes) time and O(A) output
    // storage; selection lookup is constant-time and the source map is ordered.
    let mut copied_atoms = Vec::new();
    for ref_atom in &reference.atoms {
        if !selection_info
            .selected_atoms
            .contains(ref_atom.id().index())
        {
            continue;
        }

        let new_id = AtomId::new(copied_atoms.len());
        let mut extracted_atom = ref_atom.clone().with_id(new_id);
        extracted_atom.clear_computed_props()?;

        copied_atoms.push(extracted_atom);
        selection_info.atom_mapping.insert(ref_atom.id(), new_id);
    }

    Ok(copied_atoms)
}

fn copy_selected_bonds(
    reference: &TopologyBlock,
    selection_info: &mut FragmentSubsetInfo,
) -> Result<Vec<Bond>, SelectedBondCopyError> {
    // UFF-FRAG uses copyAsQuery=false. Preserve this source bond loop as a
    // raw detached-row operation; validated topology construction happens at
    // a later owner boundary.
    // BEGIN RDKIT CPP FUNCTION copySelectedAtomsAndBonds selected bond loop
    // RDKit✔️✔️:   for (const auto &ref_bond : reference_mol.bonds()) {
    // RDKit✔️✔️:     if (!selectedBonds[ref_bond->getIdx()]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atomMapping.find(ref_bond->getBeginAtomIdx()) == atomMapping.end() ||
    // RDKit✔️✔️:         atomMapping.find(ref_bond->getEndAtomIdx()) == atomMapping.end()) {
    // RDKit✔️✔️:       throw ValueErrorException("copyMolSubset: subset bonds contain atoms not contained in subset atoms");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     std::unique_ptr<Bond> extracted_bond{
    // RDKit✔️✔️:         options.copyAsQuery ? new QueryBond(*ref_bond) : ref_bond->copy()};
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // Check the stereo atoms
    // RDKit✔️✔️:     auto &atoms = extracted_bond->getStereoAtoms();
    // RDKit✔️✔️:     if (atoms.size() == 2) {
    // RDKit✔️✔️:       auto map1 = atomMapping.find(atoms[0]);
    // RDKit✔️✔️:       auto map2 = atomMapping.find(atoms[1]);
    // RDKit✔️✔️:       if (map1 != atomMapping.end() && map2 != atomMapping.end()) {
    // RDKit✔️✔️:         atoms[0] = map1->second;
    // RDKit✔️✔️:         atoms[1] = map2->second;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         atoms.clear();  // We couldn't map the stereo atoms
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     for (auto &atomidx : atoms) {
    // RDKit✔️✔️:       auto map = atomMapping.find(atomidx);
    // RDKit✔️✔️:       if (map != atomMapping.end()) {
    // RDKit✔️✔️:         atomidx = map->second;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     extracted_bond->setBeginAtomIdx(atomMapping[ref_bond->getBeginAtomIdx()]);
    // RDKit✔️✔️:     extracted_bond->setEndAtomIdx(atomMapping[ref_bond->getEndAtomIdx()]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     constexpr bool takeOwnership = true;
    // RDKit✔️✔️:     auto num_bonds =
    // RDKit✔️✔️:         extracted_mol.addBond(extracted_bond.release(), takeOwnership);
    // RDKit✔️✔️:     bondMapping[ref_bond->getIdx()] = num_bonds - 1;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION copySelectedAtomsAndBonds selected bond loop
    // Behavior: selected source bond rows are copied in order. Endpoint map
    // checks fail before copying; stereo references get both source passes,
    // and clone/remap preserves all other typed and property state.
    // Complexity: O(B log A + C log B + copied property bytes) time and O(C)
    // output/map storage; BTreeMap matches source std::map lookup costs.
    let mut copied_bonds = Vec::new();
    for ref_bond in &reference.bonds {
        if !selection_info
            .selected_bonds
            .contains(ref_bond.id().index())
        {
            continue;
        }

        let begin = ref_bond.begin();
        let Some(&new_begin) = selection_info.atom_mapping.get(&begin) else {
            return Err(SelectedBondCopyError::MissingEndpointMapping {
                bond: ref_bond.id(),
                atom: begin,
            });
        };
        let end = ref_bond.end();
        let Some(&new_end) = selection_info.atom_mapping.get(&end) else {
            return Err(SelectedBondCopyError::MissingEndpointMapping {
                bond: ref_bond.id(),
                atom: end,
            });
        };

        let mut stereo_atoms = ref_bond.stereo_atoms();
        if let Some([first, second]) = stereo_atoms {
            stereo_atoms = match (
                selection_info.atom_mapping.get(&first),
                selection_info.atom_mapping.get(&second),
            ) {
                (Some(&new_first), Some(&new_second)) => Some([new_first, new_second]),
                _ => None,
            };
        }

        if let Some(atoms) = &mut stereo_atoms {
            for atom_index in atoms {
                if let Some(&mapped_index) = selection_info.atom_mapping.get(atom_index) {
                    *atom_index = mapped_index;
                }
            }
        }

        let new_id = BondId::new(copied_bonds.len());
        let copied_bond = ref_bond
            .clone()
            .remapped(new_id, new_begin, new_end, stereo_atoms);
        copied_bonds.push(copied_bond);
        selection_info.bond_mapping.insert(ref_bond.id(), new_id);
    }

    Ok(copied_bonds)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SubstanceGroupAtomList {
    Atoms,
    ParentAtoms,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
enum SelectedSubstanceGroupCopyError {
    #[error("subset SGroup {group:?} {list:?} atom {atom:?} has no atom mapping")]
    MissingAtomMapping {
        group: SubstanceGroupId,
        list: SubstanceGroupAtomList,
        atom: AtomId,
    },
    #[error(
        "subset SGroup {group:?} {list:?} atom {source_atom:?} maps to {mapped_atom:?}, outside atom count {atom_count}"
    )]
    AtomOutOfRange {
        group: SubstanceGroupId,
        list: SubstanceGroupAtomList,
        source_atom: AtomId,
        mapped_atom: AtomId,
        atom_count: usize,
    },
    #[error("subset SGroup {group:?} bonds bond {bond:?} has no bond mapping")]
    MissingBondMapping {
        group: SubstanceGroupId,
        bond: BondId,
    },
    #[error(
        "subset SGroup {group:?} bonds bond {source_bond:?} maps to {mapped_bond:?}, outside bond count {bond_count}"
    )]
    BondOutOfRange {
        group: SubstanceGroupId,
        source_bond: BondId,
        mapped_bond: BondId,
        bond_count: usize,
    },
}

fn is_selected_substance_group(
    group: &SubstanceGroup,
    selection_info: &FragmentSubsetInfo,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION isSelectedSGroup
    // RDKit✔️✔️: static bool isSelectedSGroup(const SubstanceGroup &sgroup,
    // RDKit✔️✔️:                              const SubsetInfo &selection_info) {
    // RDKit✔️✔️:   auto is_selected_component = [](auto &indices, auto &selection_test) {
    // RDKit✔️✔️:     return indices.empty() ||
    // RDKit✔️✔️:            std::all_of(indices.begin(), indices.end(), selection_test);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto atom_test = [&](int idx) { return selection_info.selectedAtoms[idx]; };
    // RDKit✔️✔️:   auto bond_test = [&](int idx) { return selection_info.selectedBonds[idx]; };
    // RDKit✔️✔️:   return is_selected_component(sgroup.getAtoms(), atom_test) &&
    // RDKit✔️✔️:          is_selected_component(sgroup.getBonds(), bond_test) &&
    // RDKit✔️✔️:          is_selected_component(sgroup.getParentAtoms(), atom_test);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isSelectedSGroup
    let atoms_selected = group.atoms().is_empty()
        || group
            .atoms()
            .iter()
            .all(|atom| selection_info.selected_atoms.contains(atom.index()));
    let bonds_selected = group.bonds().is_empty()
        || group
            .bonds()
            .iter()
            .all(|bond| selection_info.selected_bonds.contains(bond.index()));
    let parent_atoms_selected = group.parent_atoms().is_empty()
        || group
            .parent_atoms()
            .iter()
            .all(|atom| selection_info.selected_atoms.contains(atom.index()));
    atoms_selected && bonds_selected && parent_atoms_selected
}

fn map_substance_group_atoms(
    group: &SubstanceGroup,
    list: SubstanceGroupAtomList,
    atoms: &[AtomId],
    selection_info: &FragmentSubsetInfo,
) -> Result<Vec<AtomId>, SelectedSubstanceGroupCopyError> {
    atoms
        .iter()
        .map(|atom| {
            selection_info.atom_mapping.get(atom).copied().ok_or(
                SelectedSubstanceGroupCopyError::MissingAtomMapping {
                    group: group.id(),
                    list,
                    atom: *atom,
                },
            )
        })
        .collect()
}

fn validate_substance_group_atom_list(
    group: &SubstanceGroup,
    list: SubstanceGroupAtomList,
    source_atoms: &[AtomId],
    mapped_atoms: &[AtomId],
    atom_count: usize,
) -> Result<(), SelectedSubstanceGroupCopyError> {
    for (&source_atom, &mapped_atom) in source_atoms.iter().zip(mapped_atoms) {
        if mapped_atom.index() >= atom_count {
            return Err(SelectedSubstanceGroupCopyError::AtomOutOfRange {
                group: group.id(),
                list,
                source_atom,
                mapped_atom,
                atom_count,
            });
        }
    }
    Ok(())
}

fn map_substance_group_bonds(
    group: &SubstanceGroup,
    bonds: &[BondId],
    selection_info: &FragmentSubsetInfo,
) -> Result<Vec<BondId>, SelectedSubstanceGroupCopyError> {
    bonds
        .iter()
        .map(|bond| {
            selection_info.bond_mapping.get(bond).copied().ok_or(
                SelectedSubstanceGroupCopyError::MissingBondMapping {
                    group: group.id(),
                    bond: *bond,
                },
            )
        })
        .collect()
}

fn validate_substance_group_bonds(
    group: &SubstanceGroup,
    source_bonds: &[BondId],
    mapped_bonds: &[BondId],
    bond_count: usize,
) -> Result<(), SelectedSubstanceGroupCopyError> {
    for (&source_bond, &mapped_bond) in source_bonds.iter().zip(mapped_bonds) {
        if mapped_bond.index() >= bond_count {
            return Err(SelectedSubstanceGroupCopyError::BondOutOfRange {
                group: group.id(),
                source_bond,
                mapped_bond,
                bond_count,
            });
        }
    }
    Ok(())
}

fn copy_selected_substance_groups(
    reference: &[SubstanceGroup],
    selection_info: &FragmentSubsetInfo,
) -> Result<Vec<SubstanceGroup>, SelectedSubstanceGroupCopyError> {
    // BEGIN RDKIT CPP FUNCTION copySelectedSubstanceGroups
    // RDKit✔️❌: static void copySelectedSubstanceGroups(RWMol &extracted_mol,
    // RDKit✔️❌:                                         const RDKit::ROMol &reference_mol,
    // RDKit✔️❌:                                         const SubsetInfo &selection_info,
    // RDKit✔️❌:                                         const SubsetOptions &) {
    // RDKit✔️❌:   auto update_indices = [](auto &sgroup, auto getter, auto setter,
    // RDKit✔️❌:                            auto &mapping) {
    // RDKit✔️❌:     auto indices = getter(sgroup);
    // RDKit✔️❌:     std::for_each(indices.begin(), indices.end(),
    // RDKit✔️❌:                   [&](auto &idx) { idx = mapping.at(idx); });
    // RDKit✔️❌:     setter(sgroup, std::move(indices));
    // RDKit✔️❌:   };
    // RDKit✔️❌:
    // RDKit✔️❌:   const auto &[selectedAtoms, selectedBonds, atomMapping, bondMapping] =
    // RDKit✔️❌:       selection_info;
    // RDKit✔️❌:   for (const auto &sgroup : getSubstanceGroups(reference_mol)) {
    // RDKit✔️❌:     if (!isSelectedSGroup(sgroup, selection_info)) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     SubstanceGroup extracted_sgroup(sgroup);
    // RDKit✔️❌:     extracted_sgroup.setOwningMol(&extracted_mol);
    // RDKit✔️❌:
    // RDKit✔️❌:     update_indices(extracted_sgroup, std::mem_fn(&SubstanceGroup::getAtoms),
    // RDKit✔️❌:                    std::mem_fn(&SubstanceGroup::setAtoms), atomMapping);
    // RDKit✔️❌:     update_indices(extracted_sgroup,
    // RDKit✔️❌:                    std::mem_fn(&SubstanceGroup::getParentAtoms),
    // RDKit✔️❌:                    std::mem_fn(&SubstanceGroup::setParentAtoms), atomMapping);
    // RDKit✔️❌:     update_indices(extracted_sgroup, std::mem_fn(&SubstanceGroup::getBonds),
    // RDKit✔️❌:                    std::mem_fn(&SubstanceGroup::setBonds), bondMapping);
    // RDKit✔️❌:
    // RDKit✔️❌:     addSubstanceGroup(extracted_mol, std::move(extracted_sgroup));
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION copySelectedSubstanceGroups
    let mut copied_groups = Vec::new();
    let atom_count = selection_info.atom_mapping.len();
    let bond_count = selection_info.bond_mapping.len();

    for source_group in reference {
        if !is_selected_substance_group(source_group, selection_info) {
            continue;
        }

        let mapped_atoms = map_substance_group_atoms(
            source_group,
            SubstanceGroupAtomList::Atoms,
            source_group.atoms(),
            selection_info,
        )?;
        validate_substance_group_atom_list(
            source_group,
            SubstanceGroupAtomList::Atoms,
            source_group.atoms(),
            &mapped_atoms,
            atom_count,
        )?;
        let mapped_parent_atoms = map_substance_group_atoms(
            source_group,
            SubstanceGroupAtomList::ParentAtoms,
            source_group.parent_atoms(),
            selection_info,
        )?;
        validate_substance_group_atom_list(
            source_group,
            SubstanceGroupAtomList::ParentAtoms,
            source_group.parent_atoms(),
            &mapped_parent_atoms,
            atom_count,
        )?;
        let mapped_bonds =
            map_substance_group_bonds(source_group, source_group.bonds(), selection_info)?;
        validate_substance_group_bonds(
            source_group,
            source_group.bonds(),
            &mapped_bonds,
            bond_count,
        )?;

        // The target model stores CBONDS/XBONDS roles by BondId. Carry each
        // source member's typed role to its mapped member after cloning the
        // source's remaining fields unchanged. Rebuilding this sidecar map
        // adds an ordered-map pass/allocation relative to RDKit's topology-
        // derived bond type; retain this known target-representation cost.
        let mut extracted_group = source_group
            .clone()
            .with_atoms(mapped_atoms)
            .with_parent_atoms(mapped_parent_atoms)
            .with_bonds(Vec::new());
        for (source_bond, mapped_bond) in source_group.bonds().iter().copied().zip(mapped_bonds) {
            extracted_group.push_bond_with_role(mapped_bond, source_group.bond_role(source_bond));
        }
        extracted_group.set_id(SubstanceGroupId::new(copied_groups.len()));
        copied_groups.push(extracted_group);
    }

    Ok(copied_groups)
}

fn copy_selected_stereo_groups(
    reference: &[StereoGroup],
    selection_info: &FragmentSubsetInfo,
) -> Vec<StereoGroup> {
    // BEGIN RDKIT CPP FUNCTION copySelectedStereoGroups
    // RDKit✔️🔝: static void copySelectedStereoGroups(RWMol &extracted_mol,
    // RDKit✔️🔝:                                      const RDKit::ROMol &reference_mol,
    // RDKit✔️🔝:                                      const SubsetInfo &selection_info) {
    // RDKit✔️🔝:   auto is_selected_component = [](auto &objects, auto &selected_indices) {
    // RDKit✔️🔝:     return objects.empty() ||
    // RDKit✔️🔝:            std::any_of(objects.begin(), objects.end(), [&](auto &object) {
    // RDKit✔️🔝:              return selected_indices[object->getIdx()];
    // RDKit✔️🔝:            });
    // RDKit✔️🔝:   };
    // RDKit✔️🔝:
    // RDKit✔️🔝:   auto is_selected_stereo_group = [&](const auto &stereo_group) {
    // RDKit✔️🔝:     return is_selected_component(stereo_group.getAtoms(),
    // RDKit✔️🔝:                                  selection_info.selectedAtoms) &&
    // RDKit✔️🔝:            is_selected_component(stereo_group.getBonds(),
    // RDKit✔️🔝:                                  selection_info.selectedBonds);
    // RDKit✔️🔝:   };
    // RDKit✔️🔝:
    // RDKit✔️🔝:   std::vector<Atom *> extracted_atoms(extracted_mol.getNumAtoms());
    // RDKit✔️🔝:   for (const auto &atom : extracted_mol.atoms()) {
    // RDKit✔️🔝:     extracted_atoms[atom->getIdx()] = atom;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:
    // RDKit✔️🔝:   std::vector<Bond *> extracted_bonds(extracted_mol.getNumBonds());
    // RDKit✔️🔝:   for (const auto &bond : extracted_mol.bonds()) {
    // RDKit✔️🔝:     extracted_bonds[bond->getIdx()] = bond;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:
    // RDKit✔️🔝:   const auto &[selectedAtoms, selectedBonds, atomMapping, bondMapping] =
    // RDKit✔️🔝:       selection_info;
    // RDKit✔️🔝:   std::vector<StereoGroup> extracted_stereo_groups;
    // RDKit✔️🔝:   for (const auto &stereo_group : reference_mol.getStereoGroups()) {
    // RDKit✔️🔝:     if (!is_selected_stereo_group(stereo_group)) {
    // RDKit✔️🔝:       continue;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:
    // RDKit✔️🔝:     std::vector<Atom *> atoms;
    // RDKit✔️🔝:     for (const auto &atom : stereo_group.getAtoms()) {
    // RDKit✔️🔝:       auto mapping = atomMapping.find(atom->getIdx());
    // RDKit✔️🔝:       if (mapping != atomMapping.end()) {
    // RDKit✔️🔝:         atoms.push_back(extracted_atoms[mapping->second]);
    // RDKit✔️🔝:       }
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:
    // RDKit✔️🔝:     std::vector<Bond *> bonds;
    // RDKit✔️🔝:     for (const auto &bond : stereo_group.getBonds()) {
    // RDKit✔️🔝:       auto mapping = bondMapping.find(bond->getIdx());
    // RDKit✔️🔝:       if (mapping != bondMapping.end()) {
    // RDKit✔️🔝:         bonds.push_back(extracted_bonds[mapping->second]);
    // RDKit✔️🔝:       }
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:
    // RDKit✔️🔝:     extracted_stereo_groups.push_back({stereo_group.getGroupType(),
    // RDKit✔️🔝:                                        std::move(atoms),
    // RDKit✔️🔝:                                        std::move(bonds),
    // RDKit✔️🔝:                                        stereo_group.getReadId()});
    // RDKit✔️🔝:     extracted_stereo_groups.back().setWriteId(stereo_group.getWriteId());
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:
    // RDKit✔️🔝:   extracted_mol.setStereoGroups(std::move(extracted_stereo_groups));
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION copySelectedStereoGroups
    // Stable model row IDs are returned directly through the existing ordered
    // maps, avoiding RDKit's two temporary dense target-row pointer vectors.
    // The category/member scans and map lookup complexity remain source-shaped.
    let mut copied_groups = Vec::new();
    for source_group in reference {
        let atoms_selected = source_group.atoms().is_empty()
            || source_group
                .atoms()
                .iter()
                .any(|atom| selection_info.selected_atoms.contains(atom.index()));
        let bonds_selected = source_group.bonds().is_empty()
            || source_group
                .bonds()
                .iter()
                .any(|bond| selection_info.selected_bonds.contains(bond.index()));
        if !atoms_selected || !bonds_selected {
            continue;
        }

        let atoms = source_group
            .atoms()
            .iter()
            .filter_map(|atom| selection_info.atom_mapping.get(atom).copied())
            .collect();
        let bonds = source_group
            .bonds()
            .iter()
            .filter_map(|bond| selection_info.bond_mapping.get(bond).copied())
            .collect();
        let mut copied_group = StereoGroup::new(source_group.kind(), atoms, bonds);
        if let Some(read_id) = source_group.id() {
            copied_group = copied_group.with_id(read_id);
        }
        copied_groups.push(copied_group.with_write_id(source_group.write_id()));
    }
    copied_groups
}

fn copy_full_molecule_coordinates(source: &FragmentCoordinateView<'_>) -> CoordinateBlock {
    // BEGIN RDKIT CPP FUNCTION copyMolSubset whole-molecule copy and ROMol conformer copy
    // RDKit✔️❌:   if ((atoms.size() == natoms && bonds.size() == nbonds)) {
    // RDKit✔️❌:     // optimization to copy the entire thing
    // RDKit✔️❌:     return std::make_unique<RDKit::RWMol>(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌: RWMol(const RWMol &other) : ROMol(other) {}
    // RDKit✔️❌: ROMol(const ROMol &other, bool quickCopy = false, int confId = -1)
    // RDKit✔️❌:     : RDProps() {
    // RDKit✔️❌:   initFromOther(other, quickCopy, confId);
    // RDKit✔️❌: }
    // RDKit✔️❌: if (!quickCopy) {
    // RDKit✔️❌:   // copy conformations
    // RDKit✔️❌:   for (const auto &conf : other.d_confs) {
    // RDKit✔️❌:     if (confId < 0 || rdcast<int>(conf->getId()) == confId) {
    // RDKit✔️❌:       this->addConformer(new Conformer(*conf));
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // RDKit✔️❌: Conformer(const Conformer &other) = default;
    // RDKit✔️❌: unsigned int addConformer(Conformer *conf, bool assignId = false);
    // RDKit✔️❌: if (assignId) {
    // RDKit✔️❌:   int maxId = -1;
    // RDKit✔️❌:   for (auto cptr : d_confs) {
    // RDKit✔️❌:     maxId = std::max((int)(cptr->getId()), maxId);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   maxId++;
    // RDKit✔️❌:   conf->setId((unsigned int)maxId);
    // RDKit✔️❌: }
    // RDKit✔️❌: conf->setOwningMol(this);
    // RDKit✔️❌: CONFORMER_SPTR nConf(conf);
    // RDKit✔️❌: d_confs.push_back(nConf);
    // RDKit✔️❌: return conf->getId();
    // END RDKIT CPP FUNCTION copyMolSubset whole-molecule copy and ROMol conformer copy
    // Behavior: this branch copies every modeled dimension and conformer,
    // including IDs, dimensional flags, coordinate rows and conformer props.
    // Complexity: output rows and properties are cloned once as in the source;
    // the borrowed descriptor vectors were allocated upstream at O(C), but
    // this function creates no intermediate coordinate buffer.
    source.to_coordinate_block()
}

fn copy_single_full_molecule_component(
    source_topology: &TopologyBlock,
    source_coordinates: &CoordinateBlock,
    source_properties: &MoleculeProperties,
) -> FullCopyComponent {
    let coordinate_view = FragmentCoordinateView::from_coordinate_block(source_coordinates);
    copy_single_full_molecule_component_with_view(
        source_topology,
        &coordinate_view,
        source_properties,
    )
}

fn copy_single_full_molecule_component_with_view(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    source_properties: &MoleculeProperties,
) -> FullCopyComponent {
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags single-fragment full-copy branch
    // RDKit❗❌:   if (nFrags == 1) {
    // RDKit❗❌:     res.emplace_back(new RWMol(mol));
    // RDKit❗❌:     if (fragsMolAtomMapping) {
    // RDKit❗❌:       INT_VECT comp;
    // RDKit❗❌:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:         comp.push_back(idx);
    // RDKit❗❌:       }
    // RDKit❗❌:       (*fragsMolAtomMapping).push_back(comp);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION MolOps::getTheFrags single-fragment full-copy branch
    // Behavior: this private stage clones the complete represented topology,
    // coordinate block, and molecule properties and returns the ascending
    // identity row map. RingInfo and RDKit bookmark tables are not represented
    // by these detached model values. Conformer clearing and sanitation belong
    // to the later collection post-processing stage, not this copy branch.
    // Complexity: value cloning is linear in stored topology, coordinate and
    // property payload; identity mapping also allocates O(A+B) rows, and the
    // topology clone copies stored adjacency in addition to atom/bond rows.
    FullCopyComponent {
        topology: source_topology.clone(),
        coordinates: copy_full_molecule_coordinates(source_coordinates),
        molecule_properties: source_properties.clone(),
        mapping: TopologyMapping::identity(
            source_topology.atoms.len(),
            source_topology.bonds.len(),
        ),
        source_metadata: FragmentSourceMetadata::unmodeled(),
    }
}

#[derive(Debug, Clone, PartialEq)]
struct FullCopyComponent {
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    molecule_properties: MoleculeProperties,
    mapping: TopologyMapping,
    source_metadata: FragmentSourceMetadata,
}

#[derive(Debug, thiserror::Error)]
enum FullCopyComponentError {
    #[error(transparent)]
    SourceUInt(#[from] crate::PropertyUIntReadError),
    #[error(transparent)]
    SourceText(#[from] crate::PropertyStringError),

    #[error(transparent)]
    TopologyEdit(#[from] TopologyEditError),
    #[error(transparent)]
    MappingValidation(#[from] MappingValidationError),
    #[error(transparent)]
    CoordinateValidation(#[from] CoordinateValidationError),
    #[error("fragment atom mask has {actual} bits, expected {expected}")]
    AtomMaskSize { expected: usize, actual: usize },
    #[error("batch deletion produced atom row {row} without a source identity")]
    AddedAtomRow { row: usize },
    #[error("source atom row {source_row} is absent from conformer {conformer_id}")]
    MissingCoordinateRow {
        conformer_id: usize,
        source_row: usize,
    },
}

fn remap_full_copy_coordinate_rows<const DIMENSION: usize>(
    source_row: impl Fn(usize) -> Option<[f64; DIMENSION]>,
    atom_new_to_old: &[Option<AtomId>],
    conformer_id: usize,
) -> Result<Vec<[f64; DIMENSION]>, FullCopyComponentError> {
    let mut positions = Vec::with_capacity(atom_new_to_old.len());
    for (new_row, source_atom) in atom_new_to_old.iter().enumerate() {
        let Some(source_atom) = source_atom else {
            return Err(FullCopyComponentError::AddedAtomRow { row: new_row });
        };
        let Some(position) = source_row(source_atom.index()) else {
            return Err(FullCopyComponentError::MissingCoordinateRow {
                conformer_id,
                source_row: source_atom.index(),
            });
        };
        positions.push(position);
    }
    Ok(positions)
}

fn copy_full_copy_coordinates(
    source: &FragmentCoordinateView<'_>,
    atom_new_to_old: &[Option<AtomId>],
    source_atom_count: usize,
) -> Result<CoordinateBlock, FullCopyComponentError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::initFromOther conformer copy
    // RDKit❗❌:     for (const auto &conf : other.d_confs) {
    // RDKit❗❌:       if (confId < 0 || rdcast<int>(conf->getId()) == confId) {
    // RDKit❗❌:         this->addConformer(new Conformer(*conf));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌: Conformer(const Conformer &other) = default;
    // END RDKIT CPP FUNCTION ROMol::initFromOther conformer copy
    // BEGIN RDKIT CPP FUNCTION RWMol::batchRemoveAtoms conformer row deletion
    // RDKit❗❌:   for (auto conf : d_confs) {
    // RDKit❗❌:     RDGeom::POINT3D_VECT &positions = conf->getPositions();
    // RDKit❗❌:     RDGeom::POINT3D_VECT newPositions;
    // RDKit❗❌:     newPositions.reserve(getNumAtoms());
    // RDKit❗❌:     for (RDGeom::POINT3D_VECT::size_type i = 0; i < positions.size(); ++i) {
    // RDKit❗❌:       if (oldIndices[i] != nullptr) {
    // RDKit❗❌:         newPositions.push_back(positions[i]);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     CHECK_INVARIANT(newPositions.size() == getNumAtoms(), "Lost coordinates!");
    // RDKit❗❌:     positions.swap(newPositions);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION RWMol::batchRemoveAtoms conformer row deletion
    // Behavior: the full-copy route preserves each copied conformer's source
    // ID, dimensional flag, properties, and retained source rows in their
    // post-deletion order. `None` output rows are impossible because this
    // branch only removes source atoms; report them as a typed invariant error.
    // Source clone/delete transports IEEE values verbatim; it does not reject
    // non-finite rows. Preserve signed zero and NaN payload bits, while retaining
    // structural row-count and mapping checks before any indexed access.
    // Complexity: one output coordinate row pass per conformer and retained
    // atom, with one output vector per conformer; no discarded full-size
    // coordinate block is allocated before deletion.
    for conformer in &source.conformers_2d {
        if conformer.coordinates.len() != source_atom_count {
            return Err(FullCopyComponentError::CoordinateValidation(
                CoordinateValidationError::RowCount {
                    dimension: "2D",
                    conformer: conformer.id,
                    rows: conformer.coordinates.len(),
                    atom_count: source_atom_count,
                },
            ));
        }
    }
    for conformer in &source.conformers_3d {
        if conformer.coordinates.len() != source_atom_count {
            return Err(FullCopyComponentError::CoordinateValidation(
                CoordinateValidationError::RowCount {
                    dimension: "3D",
                    conformer: conformer.id,
                    rows: conformer.coordinates.len(),
                    atom_count: source_atom_count,
                },
            ));
        }
    }
    let conformers_2d = source
        .conformers_2d
        .iter()
        .map(|conformer| {
            let positions = remap_full_copy_coordinate_rows(
                |row| conformer.coordinates.get(row).copied(),
                atom_new_to_old,
                conformer.id,
            )?;
            let mut copied = Conformer2D::new(conformer.id, positions);
            for (key, value) in conformer.props {
                copied = copied.with_prop(key.clone(), value.clone());
            }
            Ok(copied)
        })
        .collect::<Result<Vec<_>, FullCopyComponentError>>()?;
    let conformers_3d = source
        .conformers_3d
        .iter()
        .map(|conformer| {
            let positions = remap_full_copy_coordinate_rows(
                |row| conformer.coordinates.get(row),
                atom_new_to_old,
                conformer.id,
            )?;
            let mut copied = Conformer3D::new(conformer.id, positions, conformer.is_3d);
            for (key, value) in conformer.props {
                copied = copied.with_prop(key.clone(), value.clone());
            }
            Ok(copied)
        })
        .collect::<Result<Vec<_>, FullCopyComponentError>>()?;

    Ok(CoordinateBlock {
        conformers_2d,
        conformers_3d,
        source_coordinate_dim: source.source_coordinate_dim,
        source_conformer_order: source.source_conformer_order.clone(),
    })
}

fn copy_full_molecule_remove_atoms_outside_component(
    source_topology: &TopologyBlock,
    source_coordinates: &CoordinateBlock,
    source_properties: &MoleculeProperties,
    atoms_in_fragment: &SelectionMask,
) -> Result<FullCopyComponent, FullCopyComponentError> {
    copy_full_molecule_remove_atoms_outside_component_with_view(
        source_topology,
        &FragmentCoordinateView::from_coordinate_block(source_coordinates),
        source_properties,
        atoms_in_fragment,
    )
}

fn copy_full_molecule_remove_atoms_outside_component_with_view(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    source_properties: &MoleculeProperties,
    atoms_in_fragment: &SelectionMask,
) -> Result<FullCopyComponent, FullCopyComponentError> {
    // BEGIN RDKIT CPP FUNCTION RWMol copy constructor and ROMol::initFromOther
    // RDKit❗❌: RWMol(const RWMol &other) : ROMol(other) {}
    // RDKit❗❌: ROMol(const ROMol &other, bool quickCopy = false, int confId = -1)
    // RDKit❗❌:     : RDProps() {
    // RDKit❗❌:   initFromOther(other, quickCopy, confId);
    // RDKit❗❌: }
    // RDKit❗❌: void ROMol::initFromOther(const ROMol &other, bool quickCopy, int confId) {
    // RDKit❗❌:   if (this == &other) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   numBonds = 0;
    // RDKit❗❌:   for (const auto oatom : other.atoms()) {
    // RDKit❗❌:     constexpr bool updateLabel = false;
    // RDKit❗❌:     constexpr bool takeOwnership = true;
    // RDKit❗❌:     addAtom(oatom->copy(), updateLabel, takeOwnership);
    // RDKit❗❌:   }
    // RDKit❗❌:   for (const auto obond : other.bonds()) {
    // RDKit❗❌:     addBond(obond->copy(), true);
    // RDKit❗❌:   }
    // RDKit❗❌:   d_stereo_groups.clear();
    // RDKit❗❌:   for (auto &otherGroup : other.d_stereo_groups) {
    // RDKit❗❌:     std::vector<Atom *> atoms;
    // RDKit❗❌:     for (auto &otherAtom : otherGroup.getAtoms()) {
    // RDKit❗❌:       atoms.push_back(getAtomWithIdx(otherAtom->getIdx()));
    // RDKit❗❌:     }
    // RDKit❗❌:     std::vector<Bond *> bonds;
    // RDKit❗❌:     for (auto &otherBond : otherGroup.getBonds()) {
    // RDKit❗❌:       bonds.push_back(getBondWithIdx(otherBond->getIdx()));
    // RDKit❗❌:     }
    // RDKit❗❌:     d_stereo_groups.emplace_back(otherGroup.getGroupType(), std::move(atoms),
    // RDKit❗❌:                                  std::move(bonds), otherGroup.getReadId());
    // RDKit❗❌:     d_stereo_groups.back().setWriteId(otherGroup.getWriteId());
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!quickCopy) {
    // RDKit❗❌:     for (const auto &conf : other.d_confs) {
    // RDKit❗❌:       if (confId < 0 || rdcast<int>(conf->getId()) == confId) {
    // RDKit❗❌:         this->addConformer(new Conformer(*conf));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     for (const auto &sg : getSubstanceGroups(other)) {
    // RDKit❗❌:       addSubstanceGroup(*this, sg);
    // RDKit❗❌:     }
    // RDKit❗❌:     d_props = other.d_props;
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RWMol copy constructor and ROMol::initFromOther
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags slow-copy branch
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res.emplace_back(new RWMol(mol));
    // RDKit❗❌:         auto &frag = res.back();
    // RDKit❗❌:         frag->beginBatchEdit();
    // RDKit❗❌:         for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:           if (!atomsInFrag[idx]) {
    // RDKit❗❌:             frag->removeAtom(idx);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         frag->commitBatchEdit();
    // RDKit❗❌:       }
    // END RDKIT CPP FUNCTION MolOps::getTheFrags slow-copy branch
    // BEGIN RDKIT CPP FUNCTION RWMol::commitBatchEdit changed/no-op order
    // RDKit❗❌: void RWMol::commitBatchEdit() {
    // RDKit❗❌:   if (!(dp_delBonds || dp_delAtoms)) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   } else if (dp_delBonds->none() && dp_delAtoms->none()) {
    // RDKit❗❌:     // no need to reset ring info & calculated properties,
    // RDKit❗❌:     // since nothing gets removed
    // RDKit❗❌:     dp_delBonds.reset();
    // RDKit❗❌:     dp_delAtoms.reset();
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   batchRemoveBonds();
    // RDKit❗❌:   batchRemoveAtoms();
    // RDKit❗❌:   dp_ringInfo->reset();
    // RDKit❗❌:   clearComputedProps(true);
    // RDKit❗❌:   dp_delBonds.reset();
    // RDKit❗❌:   dp_delAtoms.reset();
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RWMol::commitBatchEdit changed/no-op order
    // Behavior: this private helper corresponds only to the source-selected
    // slow route. It edits a detached, full topology copy through the existing
    // validated model batch editor; the original topology, coordinates, and
    // properties remain unchanged. Source-index order is used to schedule
    // atom removals; the editor performs the pinned descending commit work.
    // The unchanged batch early-return means computed values are cleared only
    // when at least one complementary atom is actually removed.
    // Complexity: the source copies once then deletes in place. The model
    // editor validates and clones both source and working topology, then
    // rebuilds remapped rows/adjacency, so this private bridge is materially
    // more allocation-heavy despite linear source-shaped row processing.
    if atoms_in_fragment.bit_count != source_topology.atoms.len() {
        return Err(FullCopyComponentError::AtomMaskSize {
            expected: source_topology.atoms.len(),
            actual: atoms_in_fragment.bit_count,
        });
    }

    let mut edit = source_topology.begin_batch_edit()?;
    let mut removed_atom = false;
    for atom_index in 0..source_topology.atoms.len() {
        if !atoms_in_fragment.contains(atom_index) {
            edit.remove_atom(AtomId::new(atom_index))?;
            removed_atom = true;
        }
    }
    let (mut topology, mapping) = edit.finish()?;
    mapping.validate_for_counts(
        source_topology.atoms.len(),
        topology.atoms.len(),
        source_topology.bonds.len(),
        topology.bonds.len(),
    )?;
    let coordinates = copy_full_copy_coordinates(
        source_coordinates,
        mapping.atoms().new_to_old(),
        source_topology.atoms.len(),
    )?;
    let mut molecule_properties = source_properties.clone();
    molecule_properties.remap_topology(mapping.atoms().new_to_old(), mapping.bonds().new_to_old());
    if removed_atom {
        clear_subset_computed_props(&mut topology, &mut molecule_properties)?;
    }

    Ok(FullCopyComponent {
        topology,
        coordinates,
        molecule_properties,
        mapping,
        source_metadata: FragmentSourceMetadata::unmodeled(),
    })
}

fn copy_full_molecule_remove_atoms_outside_component_with_source_metadata(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    source_properties: &MoleculeProperties,
    atoms_in_fragment: &SelectionMask,
    source_metadata: FragmentSourceMetadataView<'_>,
) -> Result<FullCopyComponent, FullCopyComponentError> {
    // BEGIN RDKIT CPP FUNCTION RWMol copy constructor and ROMol::initFromOther
    // RDKit❗❌: RWMol(const RWMol &other) : ROMol(other) {}
    // RDKit❗❌: ROMol(const ROMol &other, bool quickCopy = false, int confId = -1)
    // RDKit❗❌:     : RDProps() {
    // RDKit❗❌:   initFromOther(other, quickCopy, confId);
    // RDKit❗❌: }
    // RDKit❗❌: void ROMol::initFromOther(const ROMol &other, bool quickCopy, int confId) {
    // RDKit❗❌:   if (this == &other) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   numBonds = 0;
    // RDKit❗❌:   for (const auto oatom : other.atoms()) {
    // RDKit❗❌:     constexpr bool updateLabel = false;
    // RDKit❗❌:     constexpr bool takeOwnership = true;
    // RDKit❗❌:     addAtom(oatom->copy(), updateLabel, takeOwnership);
    // RDKit❗❌:   }
    // RDKit❗❌:   for (const auto obond : other.bonds()) {
    // RDKit❗❌:     addBond(obond->copy(), true);
    // RDKit❗❌:   }
    // RDKit❗❌:   d_stereo_groups.clear();
    // RDKit❗❌:   for (auto &otherGroup : other.d_stereo_groups) {
    // RDKit❗❌:     std::vector<Atom *> atoms;
    // RDKit❗❌:     for (auto &otherAtom : otherGroup.getAtoms()) {
    // RDKit❗❌:       atoms.push_back(getAtomWithIdx(otherAtom->getIdx()));
    // RDKit❗❌:     }
    // RDKit❗❌:     std::vector<Bond *> bonds;
    // RDKit❗❌:     for (auto &otherBond : otherGroup.getBonds()) {
    // RDKit❗❌:       bonds.push_back(getBondWithIdx(otherBond->getIdx()));
    // RDKit❗❌:     }
    // RDKit❗❌:     d_stereo_groups.emplace_back(otherGroup.getGroupType(), std::move(atoms),
    // RDKit❗❌:                                  std::move(bonds), otherGroup.getReadId());
    // RDKit❗❌:     d_stereo_groups.back().setWriteId(otherGroup.getWriteId());
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!quickCopy) {
    // RDKit❗❌:     for (const auto &conf : other.d_confs) {
    // RDKit❗❌:       if (confId < 0 || rdcast<int>(conf->getId()) == confId) {
    // RDKit❗❌:         this->addConformer(new Conformer(*conf));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     for (const auto &sg : getSubstanceGroups(other)) {
    // RDKit❗❌:       addSubstanceGroup(*this, sg);
    // RDKit❗❌:     }
    // RDKit❗❌:     d_props = other.d_props;
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RWMol copy constructor and ROMol::initFromOther
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags slow-copy branch
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res.emplace_back(new RWMol(mol));
    // RDKit❗❌:         auto &frag = res.back();
    // RDKit❗❌:         frag->beginBatchEdit();
    // RDKit❗❌:         for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:           if (!atomsInFrag[idx]) {
    // RDKit❗❌:             frag->removeAtom(idx);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         frag->commitBatchEdit();
    // RDKit❗❌:       }
    // END RDKIT CPP FUNCTION MolOps::getTheFrags slow-copy branch
    // BEGIN RDKIT CPP FUNCTION RWMol::commitBatchEdit changed/no-op order
    // RDKit❗❌: void RWMol::commitBatchEdit() {
    // RDKit❗❌:   if (!(dp_delBonds || dp_delAtoms)) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   } else if (dp_delBonds->none() && dp_delAtoms->none()) {
    // RDKit❗❌:     // no need to reset ring info & calculated properties,
    // RDKit❗❌:     // since nothing gets removed
    // RDKit❗❌:     dp_delBonds.reset();
    // RDKit❗❌:     dp_delAtoms.reset();
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   batchRemoveBonds();
    // RDKit❗❌:   batchRemoveAtoms();
    // RDKit❗❌:   dp_ringInfo->reset();
    // RDKit❗❌:   clearComputedProps(true);
    // RDKit❗❌:   dp_delBonds.reset();
    // RDKit❗❌:   dp_delAtoms.reset();
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RWMol::commitBatchEdit changed/no-op order
    // Behavior: the source slow route copies raw coordinates and real optional
    // source metadata, schedules complement atoms in ascending order, and
    // invokes the canonical MODEL source batch commit. No finite-coordinate
    // guard or final stereo-cardinality validation precedes native filtering.
    // Unmodeled metadata remains explicitly absent; represented metadata is
    // copied, deleted, remapped, and reset in source order. Typed molecule row
    // projections still add an explicitly recorded transport step below.
    // Complexity: one topology/coordinate/property copy followed by native
    // descending removal. MODEL adjacency rebuilding, per-row Vec erasure,
    // mapping creation, and typed property transport add material costs versus
    // native pointer graph deletion; no second detached topology clone occurs.
    if atoms_in_fragment.bit_count != source_topology.atoms.len() {
        return Err(FullCopyComponentError::AtomMaskSize {
            expected: source_topology.atoms.len(),
            actual: atoms_in_fragment.bit_count,
        });
    }

    let mut edit = source_topology.begin_batch_edit()?;
    for atom_index in 0..source_topology.atoms.len() {
        if !atoms_in_fragment.contains(atom_index) {
            edit.remove_atom(AtomId::new(atom_index))?;
        }
    }
    let (mut topology, mut atoms, mut bonds) = edit.into_source_batch_parts();
    let mut coordinates = copy_full_molecule_coordinates(source_coordinates);
    let mut molecule_properties = source_properties.clone();
    let mut metadata = clone_source_fragment_metadata(source_topology, source_metadata)?;
    let mut ring_reset = metadata.rings.as_mut().map(|ring| move || ring.reset());
    cosmolkit_model::commit_batch_edit_source(
        &mut topology,
        cosmolkit_model::SourceBatchCommitState {
            atoms: &mut atoms,
            bonds: &mut bonds,
            atom_bookmarks: metadata.atom_bookmarks.as_mut(),
            bond_bookmarks: metadata.bond_bookmarks.as_mut(),
            coordinates: &mut coordinates,
            properties: &mut molecule_properties,
            reset_ring: ring_reset.as_mut().map(|reset| reset as &mut dyn FnMut()),
            uint_reader: &mut |value| {
                crate::property_value_to_uint(value).map_err(FullCopyComponentError::from)
            },
            text_reader: &mut |value| {
                crate::property_value_to_string(value).map_err(FullCopyComponentError::from)
            },
        },
    )?;
    drop(ring_reset);
    let mut atom_old_to_new = vec![None; source_topology.atoms.len()];
    let mut atom_new_to_old = Vec::with_capacity(topology.atoms.len());
    for old in 0..source_topology.atoms.len() {
        if atoms_in_fragment.contains(old) {
            let new = AtomId::new(atom_new_to_old.len());
            atom_old_to_new[old] = Some(new);
            atom_new_to_old.push(Some(AtomId::new(old)));
        }
    }
    let mut bond_old_to_new = vec![None; source_topology.bonds.len()];
    let mut bond_new_to_old = Vec::with_capacity(topology.bonds.len());
    for bond in &source_topology.bonds {
        if atoms_in_fragment.contains(bond.begin().index())
            && atoms_in_fragment.contains(bond.end().index())
        {
            let new = BondId::new(bond_new_to_old.len());
            bond_old_to_new[bond.id().index()] = Some(new);
            bond_new_to_old.push(Some(bond.id()));
        }
    }
    let mapping = TopologyMapping {
        atoms: AtomMapping {
            old_to_new: atom_old_to_new,
            new_to_old: atom_new_to_old,
        },
        bonds: BondMapping {
            old_to_new: bond_old_to_new,
            new_to_old: bond_new_to_old,
        },
    };
    mapping.validate_for_counts(
        source_topology.atoms.len(),
        topology.atoms.len(),
        source_topology.bonds.len(),
        topology.bonds.len(),
    )?;
    // Existing typed row projections follow the returned structural mapping;
    // native generic property dictionaries remain the canonical property owner.
    // This extra modeled row transport remains an explicit source-boundary gap.
    molecule_properties.remap_topology(mapping.atoms().new_to_old(), mapping.bonds().new_to_old());

    Ok(FullCopyComponent {
        topology,
        coordinates,
        molecule_properties,
        mapping,
        source_metadata: metadata,
    })
}

fn copy_selected_coordinate_rows<const DIMENSION: usize>(
    source_positions: &[[f64; DIMENSION]],
    atom_mapping: &BTreeMap<AtomId, AtomId>,
) -> Vec<[f64; DIMENSION]> {
    let mut copied_positions = vec![[0.0; DIMENSION]; atom_mapping.len()];
    for (source_atom, target_atom) in atom_mapping {
        copied_positions[target_atom.index()] = source_positions[source_atom.index()];
    }
    copied_positions
}

fn copy_subset_coordinates(
    source: &CoordinateBlock,
    atom_mapping: &BTreeMap<AtomId, AtomId>,
) -> CoordinateBlock {
    let coordinate_view = FragmentCoordinateView::from_coordinate_block(source);
    copy_subset_coordinates_with_view(&coordinate_view, atom_mapping)
}

fn copy_subset_coordinates_with_view(
    source: &FragmentCoordinateView<'_>,
    atom_mapping: &BTreeMap<AtomId, AtomId>,
) -> CoordinateBlock {
    // BEGIN RDKIT CPP FUNCTION copyCoords
    // RDKit✔️✔️: void copyCoords(RDKit::RWMol &copy, const RDKit::RWMol &mol,
    // RDKit✔️✔️:                 SubsetInfo &subset_info) {
    // RDKit✔️✔️:   if (mol.getNumConformers()) {
    // RDKit✔️✔️:     // copy coordinates over:
    // RDKit✔️✔️:     for (auto confIt = mol.beginConformers(); confIt != mol.endConformers();
    // RDKit✔️✔️:          ++confIt) {
    // RDKit✔️✔️:       auto *conf = new Conformer(copy.getNumAtoms());
    // RDKit✔️✔️:       conf->set3D((*confIt)->is3D());
    // RDKit✔️✔️:       for (auto &mapping : subset_info.atomMapping) {
    // RDKit✔️✔️:         conf->setAtomPos(mapping.second, (*confIt)->getAtomPos(mapping.first));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       conf->setId((*confIt)->getId());
    // RDKit✔️✔️:       copy.addConformer(conf, false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION copyCoords
    // Behavior: build fresh conformers in dimension-vector order, retain each
    // source ID/3D flag, and write selected source rows into mapped target rows.
    // Fresh constructors intentionally leave subset conformer property maps
    // empty, matching the source's new Conformer(numAtoms) path.
    // Complexity: O(CK) coordinate work and O(CK) output storage for C
    // conformers and K selected atoms; the source and its conformer descriptors
    // stay borrowed, and no full-size coordinate vectors are cloned.
    let conformers_2d = source
        .conformers_2d
        .iter()
        .map(|conformer| {
            Conformer2D::new(
                conformer.id,
                copy_selected_coordinate_rows(conformer.coordinates, atom_mapping),
            )
        })
        .collect();
    let conformers_3d = source
        .conformers_3d
        .iter()
        .map(|conformer| {
            let mut coordinates = vec![[0.0; 3]; atom_mapping.len()];
            for (source_atom, target_atom) in atom_mapping {
                coordinates[target_atom.index()] = conformer
                    .coordinates
                    .get(source_atom.index())
                    .expect("source conformer contains each selected atom row");
            }
            Conformer3D::new(conformer.id, coordinates, conformer.is_3d)
        })
        .collect();

    CoordinateBlock {
        conformers_2d,
        conformers_3d,
        source_coordinate_dim: source.source_coordinate_dim,
        source_conformer_order: source.source_conformer_order.clone(),
    }
}

fn sanitize_subset_if_requested(
    topology: &mut TopologyBlock,
    sanitize: bool,
) -> Result<(), SanitizeError> {
    // BEGIN RDKIT CPP FUNCTION copyMolSubset sanitize option and MolOps::sanitizeMol default
    // RDKit✔️❌:   if (options.sanitize) {
    // RDKit✔️❌:     MolOps::sanitizeMol(*extracted_mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌: void sanitizeMol(RWMol &mol) {
    // RDKit✔️❌:   unsigned int failedOp = 0;
    // RDKit✔️❌:   sanitizeMol(mol, failedOp, SANITIZE_ALL);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION copyMolSubset sanitize option and MolOps::sanitizeMol default
    // Behavior: false skips sanitation and topology validation. True reuses
    // the existing default source-stage implementation; only a successful
    // returned topology is installed, while the original SanitizeError stage
    // and typed cause pass through unchanged.
    // Complexity: the disabled branch is O(1) and allocates nothing. Enabled
    // sanitation retains the detached pipeline's known extra topology/cloned
    // stage allocations versus RDKit's in-place mutation; this wrapper adds
    // no second topology clone.
    if sanitize {
        let sanitized = sanitize_topology(topology, &SanitizeParams::default())?;
        *topology = sanitized.topology;
    }
    Ok(())
}

fn clear_subset_computed_props(
    topology: &mut TopologyBlock,
    molecule_properties: &mut MoleculeProperties,
) -> Result<(), TopologyEditError> {
    // BEGIN RDKIT CPP FUNCTION copyMolSubset clearComputedProps=true and ROMol::clearComputedProps
    // RDKit✔️🔝: if (options.clearComputedProps) {
    // RDKit✔️🔝:   // this clears atom/bond and molecule computed props
    // RDKit✔️🔝:   extracted_mol->clearComputedProps();
    // RDKit✔️🔝: }
    // RDKit✔️🔝: void ROMol::clearComputedProps(bool includeRings) const {
    // RDKit✔️🔝:   // the SSSR information:
    // RDKit✔️🔝:   if (includeRings) {
    // RDKit✔️🔝:     this->dp_ringInfo->reset();
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:
    // RDKit✔️🔝:   RDProps::clearComputedProps();
    // RDKit✔️🔝:
    // RDKit✔️🔝:   for (auto atom : atoms()) {
    // RDKit✔️🔝:     atom->clearComputedProps();
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:
    // RDKit✔️🔝:   for (auto bond : bonds()) {
    // RDKit✔️🔝:     bond->clearComputedProps();
    // RDKit✔️🔝:   }
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION copyMolSubset clearComputedProps=true and ROMol::clearComputedProps
    // BEGIN RDKIT CPP FUNCTION RDProps::clearComputedProps
    // RDKit✔️🔝: void clearComputedProps() const {
    // RDKit✔️🔝:   STR_VECT compLst;
    // RDKit✔️🔝:   if (getPropIfPresent(RDKit::detail::computedPropName, compLst) &&
    // RDKit✔️🔝:       !compLst.empty()) {
    // RDKit✔️🔝:     for (const auto &sv : compLst) {
    // RDKit✔️🔝:       d_props.clearVal(sv);
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     compLst.clear();
    // RDKit✔️🔝:     d_props.setVal(RDKit::detail::computedPropName, compLst);
    // RDKit✔️🔝:   }
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION RDProps::clearComputedProps
    // BEGIN RDKIT CPP FUNCTION Dict::clearVal
    // RDKit✔️🔝: void clearVal(const std::string_view what) {
    // RDKit✔️🔝:   for (auto it = _data.begin(); it < _data.end(); ++it) {
    // RDKit✔️🔝:     if (it->key == what) {
    // RDKit✔️🔝:       if (_hasNonPodData) {
    // RDKit✔️🔝:         RDValue::cleanup_rdvalue(it->val);
    // RDKit✔️🔝:       }
    // RDKit✔️🔝:       _data.erase(it);
    // RDKit✔️🔝:       return;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   }
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION Dict::clearVal
    // Behavior: the detached representation has no persistent RingInfo cache
    // on TopologyBlock, so there is no attached ring assignment to reset.
    // Clear molecule, atom and bond computed values in source order while
    // preserving all unmarked properties and chemical fields. Existing
    // carrier helpers use computed-key membership, not property-name guesses.
    // Complexity: one atom and bond traversal plus each carrier's computed
    // keys. BTreeMap removal avoids RDKit Dict's per-key linear vector scan;
    // this preserves key membership semantics with logarithmic lookup.
    molecule_properties.clear_computed_props()?;
    for atom in &mut topology.atoms {
        atom.clear_computed_props()?;
    }
    for bond in &mut topology.bonds {
        bond.clear_computed_props()
            .map_err(TopologyEditError::InvalidBond)?;
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
struct AtomPathSubsetCopy {
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    molecule_properties: MoleculeProperties,
    mapping: TopologyMapping,
}

#[derive(Debug, thiserror::Error)]
enum AtomPathSubsetCopyError {
    #[error(transparent)]
    ComputedProperties(#[from] TopologyEditError),
    #[error("molecule property operation failed: {0}")]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
    #[error("bond property operation failed: {0}")]
    BondProperty(#[from] cosmolkit_model::BondValueError),
    #[error("atom property operation failed: {0}")]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error(transparent)]
    BondCopy(#[from] SelectedBondCopyError),
    #[error(transparent)]
    SubstanceGroupCopy(#[from] SelectedSubstanceGroupCopyError),
    #[error(transparent)]
    TopologyValidation(#[from] TopologyValidationError),
    #[error(transparent)]
    MappingValidation(#[from] MappingValidationError),
    #[error(transparent)]
    Sanitize(#[from] SanitizeError),
}

fn topology_mapping_from_subset_info(
    selection_info: &FragmentSubsetInfo,
    source_atom_count: usize,
    source_bond_count: usize,
) -> Result<TopologyMapping, MappingValidationError> {
    let mut atom_old_to_new = vec![None; source_atom_count];
    let mut atom_new_to_old = vec![None; selection_info.atom_mapping.len()];
    for (&old_atom, &new_atom) in &selection_info.atom_mapping {
        let Some(slot) = atom_old_to_new.get_mut(old_atom.index()) else {
            return Err(MappingValidationError::Length {
                entity: "atom",
                direction: "old_to_new",
                actual: old_atom.index() + 1,
                expected: source_atom_count,
            });
        };
        *slot = Some(new_atom);
        let Some(slot) = atom_new_to_old.get_mut(new_atom.index()) else {
            return Err(MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "old_to_new",
                row: old_atom.index(),
                mapped: new_atom.index(),
                target_count: atom_new_to_old.len(),
            });
        };
        *slot = Some(old_atom);
    }

    let mut bond_old_to_new = vec![None; source_bond_count];
    let mut bond_new_to_old = vec![None; selection_info.bond_mapping.len()];
    for (&old_bond, &new_bond) in &selection_info.bond_mapping {
        let Some(slot) = bond_old_to_new.get_mut(old_bond.index()) else {
            return Err(MappingValidationError::Length {
                entity: "bond",
                direction: "old_to_new",
                actual: old_bond.index() + 1,
                expected: source_bond_count,
            });
        };
        *slot = Some(new_bond);
        let Some(slot) = bond_new_to_old.get_mut(new_bond.index()) else {
            return Err(MappingValidationError::OutOfRange {
                entity: "bond",
                direction: "old_to_new",
                row: old_bond.index(),
                mapped: new_bond.index(),
                target_count: bond_new_to_old.len(),
            });
        };
        *slot = Some(old_bond);
    }

    let mapping = TopologyMapping {
        atoms: AtomMapping {
            old_to_new: atom_old_to_new,
            new_to_old: atom_new_to_old,
        },
        bonds: BondMapping {
            old_to_new: bond_old_to_new,
            new_to_old: bond_new_to_old,
        },
    };
    mapping.validate_for_counts(
        source_atom_count,
        selection_info.atom_mapping.len(),
        source_bond_count,
        selection_info.bond_mapping.len(),
    )?;
    Ok(mapping)
}

fn copy_mol_subset_atom_path(
    source_topology: &TopologyBlock,
    source_coordinates: &CoordinateBlock,
    path: &[AtomId],
    sanitize: bool,
    copy_coordinates: bool,
) -> Result<AtomPathSubsetCopy, AtomPathSubsetCopyError> {
    let coordinate_view = FragmentCoordinateView::from_coordinate_block(source_coordinates);
    copy_mol_subset_atom_path_with_view(
        source_topology,
        &coordinate_view,
        path,
        sanitize,
        copy_coordinates,
    )
}

fn copy_mol_subset_atom_path_with_view(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    path: &[AtomId],
    sanitize: bool,
    copy_coordinates: bool,
) -> Result<AtomPathSubsetCopy, AtomPathSubsetCopyError> {
    // This is the private getTheFrags BONDS_BETWEEN_ATOMS specialization:
    // copyAsQuery=false and clearComputedProps=true. It deliberately accepts
    // the atom path, so it cannot use the separate atom/bond-vector fast copy.
    // BEGIN RDKIT CPP FUNCTION copyMolSubset path overload and selection pipeline
    // RDKit✔️❌: std::unique_ptr<RDKit::RWMol> copyMolSubset(
    // RDKit✔️❌:     const RDKit::ROMol &mol, const std::vector<unsigned int> &path,
    // RDKit✔️❌:     SubsetInfo &selection_info, const SubsetOptions &options) {
    // RDKit✔️❌:   getSubsetInfo(selection_info, mol, path, options);
    // RDKit✔️❌:   auto res = copyMolSubset(mol, selection_info, options);
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // RDKit✔️❌: std::unique_ptr<RDKit::RWMol> copyMolSubset(const RDKit::ROMol &mol,
    // RDKit✔️❌:                                             SubsetInfo &selection_info,
    // RDKit✔️❌:                                             const SubsetOptions &options) {
    // RDKit✔️❌:   auto extracted_mol = std::make_unique<RWMol>();
    // RDKit✔️❌:   copySelectedAtomsAndBonds(*extracted_mol, mol, selection_info, options);
    // RDKit✔️❌:   copySelectedSubstanceGroups(*extracted_mol, mol, selection_info, options);
    // RDKit✔️❌:   copySelectedStereoGroups(*extracted_mol, mol, selection_info);
    // RDKit✔️❌:   if (options.copyCoordinates) {
    // RDKit✔️❌:     copyCoords(*extracted_mol, mol, selection_info);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (options.sanitize) {
    // RDKit✔️❌:     MolOps::sanitizeMol(*extracted_mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (options.clearComputedProps) {
    // RDKit✔️❌:     // this clears atom/bond and molecule computed props
    // RDKit✔️❌:     extracted_mol->clearComputedProps();
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   return extracted_mol;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION copyMolSubset path overload and selection pipeline
    // Behavior: each source stage is delegated once in pinned order. The
    // detached model validates the assembled topology and requires an owned
    // adjacency block; enabled sanitization also carries its known detached
    // allocation cost. These costs do not alter the source option branches.
    // Complexity: selection and row copying remain O(V+B) plus ordered-map
    // costs; copied coordinates are O(CK). Topology validation/adjacency and
    // the existing sanitized detached pipeline add work beyond RDKit's
    // in-place RWMol construction, hence the known performance loss marker.
    let mut selection_info = FragmentSubsetInfo::default();
    get_subset_info_for_atom_path(source_topology, path, &mut selection_info);

    let atoms = copy_selected_atoms(source_topology, &mut selection_info)?;
    let bonds = copy_selected_bonds(source_topology, &mut selection_info)?;
    let substance_groups =
        copy_selected_substance_groups(&source_topology.substance_groups, &selection_info)?;
    let stereo_groups =
        copy_selected_stereo_groups(&source_topology.stereo_groups, &selection_info);

    let mut topology =
        TopologyBlock::try_from_parts(atoms, bonds, substance_groups, stereo_groups)?;
    let mapping = topology_mapping_from_subset_info(
        &selection_info,
        source_topology.atoms.len(),
        source_topology.bonds.len(),
    )?;
    let coordinates = if copy_coordinates {
        copy_subset_coordinates_with_view(source_coordinates, &selection_info.atom_mapping)
    } else {
        CoordinateBlock::default()
    };
    sanitize_subset_if_requested(&mut topology, sanitize)?;

    // A new source RWMol starts with no copied molecule-level property map;
    // clearComputedProps=true then runs after optional sanitation.
    let mut molecule_properties = MoleculeProperties::default();
    clear_subset_computed_props(&mut topology, &mut molecule_properties)?;

    Ok(AtomPathSubsetCopy {
        topology,
        coordinates,
        molecule_properties,
        mapping,
    })
}

fn fragment_has_challenging_features(
    topology: &TopologyBlock,
    component: &[AtomId],
    atoms_in_fragment: &SelectionMask,
) -> bool {
    // The private caller supplies a validated topology: atom and bond IDs are
    // their row indices, and all group references are in range.
    // BEGIN RDKIT CPP FUNCTION getTheFrags::fragmentHasChallengingFeatures
    // RDKit✔️✔️: auto fragmentHasChallengingFeatures =
    // RDKit✔️✔️:     [&](const INT_VECT &comp,
    // RDKit✔️✔️:         const boost::dynamic_bitset<> &atomsInFrag) -> bool {
    // RDKit✔️✔️:   for (auto idx : comp) {
    // RDKit✔️✔️:     // check for atoms with stereochem:
    // RDKit✔️✔️:     const auto atom = mol.getAtomWithIdx(idx);
    // RDKit✔️✔️:     if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:         atom->getChiralTag() != Atom::ChiralType::CHI_OTHER) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (auto bnd : mol.atomBonds(atom)) {
    // RDKit✔️✔️:       if (atomsInFrag[bnd->getOtherAtomIdx(idx)]) {
    // RDKit✔️✔️:         if (bnd->getStereo() != Bond::BondStereo::STEREONONE &&
    // RDKit✔️✔️:             bnd->getStereo() != Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:           return true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto sgroup : getSubstanceGroups(mol)) {
    // RDKit✔️✔️:     for (auto aid : sgroup.getAtoms()) {
    // RDKit✔️✔️:       if (atomsInFrag[aid]) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (auto aid : sgroup.getParentAtoms()) {
    // RDKit✔️✔️:       if (atomsInFrag[aid]) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto stereoGroup : mol.getStereoGroups()) {
    // RDKit✔️✔️:     // doesn't seem like this should be necessary, but in case
    // RDKit✔️✔️:     // we ever need stereogroups where the atoms aren't marked
    // RDKit✔️✔️:     // with stereo...
    // RDKit✔️✔️:     for (auto atom : stereoGroup.getAtoms()) {
    // RDKit✔️✔️:       if (atomsInFrag[atom->getIdx()]) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // same check for stereo groups involving bonds:
    // RDKit✔️✔️:     for (auto bond : stereoGroup.getBonds()) {
    // RDKit✔️✔️:       if (atomsInFrag[bond->getBeginAtomIdx()] &&
    // RDKit✔️✔️:           atomsInFrag[bond->getEndAtomIdx()]) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION getTheFrags::fragmentHasChallengingFeatures
    // Behavior: preserve source atom, incident-bond, SGroup and stereo-group
    // order and short-circuit points; fixed F11 regressions cover every
    // represented chiral/stereo variant and each source group-member branch.
    // Complexity: component and incident-edge scans plus full group-member
    // scans match the source structure; mask tests and bond-row resolution are
    // constant-time, and this helper allocates no temporary collections.
    for atom_id in component {
        let atom_index = atom_id.index();
        let atom = &topology.atoms[atom_index];
        if !matches!(atom.chiral_tag(), ChiralTag::Unspecified | ChiralTag::Other) {
            return true;
        }

        for neighbor in topology.adjacency.neighbors_of(atom_index) {
            if atoms_in_fragment.contains(neighbor.atom_index) {
                let bond = &topology.bonds[neighbor.bond.index()];
                if !matches!(bond.stereo(), BondStereo::None | BondStereo::Any) {
                    return true;
                }
            }
        }
    }

    for substance_group in &topology.substance_groups {
        for atom in substance_group.atoms() {
            if atoms_in_fragment.contains(atom.index()) {
                return true;
            }
        }
        for atom in substance_group.parent_atoms() {
            if atoms_in_fragment.contains(atom.index()) {
                return true;
            }
        }
    }

    for stereo_group in &topology.stereo_groups {
        for atom in stereo_group.atoms() {
            if atoms_in_fragment.contains(atom.index()) {
                return true;
            }
        }
        for bond_id in stereo_group.bonds() {
            let bond = &topology.bonds[bond_id.index()];
            if atoms_in_fragment.contains(bond.begin().index())
                && atoms_in_fragment.contains(bond.end().index())
            {
                return true;
            }
        }
    }

    false
}

#[derive(Debug, Clone, PartialEq)]
struct OrderedFragmentCopy {
    component_atoms: Vec<AtomId>,
    copy: FullCopyComponent,
}

#[derive(Debug, thiserror::Error)]
enum OrderedFragmentBuildError {
    #[error(transparent)]
    ComponentLabels(#[from] crate::paths::PathError),
    #[error("fast subset construction failed for component {component_index}: {source}")]
    FastSubset {
        component_index: usize,
        #[source]
        source: AtomPathSubsetCopyError,
    },
    #[error("full-copy deletion failed for component {component_index}: {source}")]
    SlowFullCopy {
        component_index: usize,
        #[source]
        source: FullCopyComponentError,
    },
}

#[derive(Debug, thiserror::Error)]
enum MoleculeFragmentsFailure {
    #[error("computed property clearing failed for fragment {component_index}: {source}")]
    FinalComputedProperties {
        component_index: usize,
        #[source]
        source: cosmolkit_model::MoleculePropertyError,
    },
    #[error(transparent)]
    Build(#[from] OrderedFragmentBuildError),
    #[error("final sanitation failed for fragment {component_index}: {source}")]
    FinalSanitize {
        component_index: usize,
        #[source]
        source: SanitizeError,
    },
}

/// A source-ordered, detached connected component and its source-row mapping.
///
/// This core value is consumed by sibling algorithm crates. It does not wrap
/// or retain a live Molecule or runtime state.
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeFragment {
    component_atoms: Vec<AtomId>,
    copy: FullCopyComponent,
}

impl MoleculeFragment {
    /// Consume the detached fragment values, without cloning graph or coordinates.
    pub fn into_parts(self) -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
        (
            self.copy.topology,
            self.copy.coordinates,
            self.copy.molecule_properties,
        )
    }

    /// Source atom rows in the component, in ascending order.
    pub fn component_atoms(&self) -> &[AtomId] {
        &self.component_atoms
    }

    /// The detached component topology.
    pub fn topology(&self) -> &TopologyBlock {
        &self.copy.topology
    }

    /// The detached component coordinates and conformer metadata.
    pub fn coordinates(&self) -> &CoordinateBlock {
        &self.copy.coordinates
    }

    /// The detached molecule-level properties copied for this component.
    pub fn molecule_properties(&self) -> &MoleculeProperties {
        &self.copy.molecule_properties
    }

    /// Explicitly modeled copied source metadata; missing capabilities stay None.
    #[doc(hidden)]
    pub fn source_metadata(&self) -> &FragmentSourceMetadata {
        &self.copy.source_metadata
    }

    /// The validated source-to-component and component-to-source row maps.
    pub fn topology_mapping(&self) -> &TopologyMapping {
        &self.copy.mapping
    }
}

/// A typed component-copy or final-sanitization failure.
#[derive(Debug)]
pub struct MoleculeFragmentsError {
    failure: MoleculeFragmentsFailure,
}

impl MoleculeFragmentsError {
    /// Component index for component-local failures; component labeling errors
    /// occur before a component index exists.
    pub fn component_index(&self) -> Option<usize> {
        match &self.failure {
            MoleculeFragmentsFailure::Build(OrderedFragmentBuildError::ComponentLabels(_)) => None,
            MoleculeFragmentsFailure::Build(OrderedFragmentBuildError::FastSubset {
                component_index,
                ..
            })
            | MoleculeFragmentsFailure::Build(OrderedFragmentBuildError::SlowFullCopy {
                component_index,
                ..
            })
            | MoleculeFragmentsFailure::FinalComputedProperties {
                component_index, ..
            }
            | MoleculeFragmentsFailure::FinalSanitize {
                component_index, ..
            } => Some(*component_index),
        }
    }
}

impl std::fmt::Display for MoleculeFragmentsError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{}", self.failure)
    }
}

impl std::error::Error for MoleculeFragmentsError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.failure)
    }
}

impl From<OrderedFragmentBuildError> for MoleculeFragmentsError {
    fn from(failure: OrderedFragmentBuildError) -> Self {
        Self {
            failure: MoleculeFragmentsFailure::Build(failure),
        }
    }
}

fn build_ordered_fragment_copies(
    source_topology: &TopologyBlock,
    source_coordinates: &CoordinateBlock,
    source_properties: &MoleculeProperties,
    sanitize: bool,
    copy_conformers: bool,
) -> Result<Vec<OrderedFragmentCopy>, OrderedFragmentBuildError> {
    let coordinate_view = FragmentCoordinateView::from_coordinate_block(source_coordinates);
    build_ordered_fragment_copies_with_view(
        source_topology,
        &coordinate_view,
        source_properties,
        sanitize,
        copy_conformers,
    )
}

fn build_ordered_fragment_copies_with_view(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    source_properties: &MoleculeProperties,
    sanitize: bool,
    copy_conformers: bool,
) -> Result<Vec<OrderedFragmentCopy>, OrderedFragmentBuildError> {
    let connected = crate::paths::connected_components(source_topology)?;
    build_ordered_fragment_copies_from_components(
        source_topology,
        source_coordinates,
        source_properties,
        sanitize,
        copy_conformers,
        connected,
        None,
        FragmentSourceMetadataView::unmodeled(),
        false,
    )
}

fn build_ordered_fragment_copies_from_components(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    source_properties: &MoleculeProperties,
    sanitize: bool,
    copy_conformers: bool,
    connected: crate::paths::ConnectedComponents,
    mut component_output: Option<&mut Vec<Vec<i32>>>,
    source_metadata: FragmentSourceMetadataView<'_>,
    source_boundary: bool,
) -> Result<Vec<OrderedFragmentCopy>, OrderedFragmentBuildError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags ordered component construction
    // RDKit❗❌: int nFrags = getMolFrags(mol, *frags);
    // RDKit❗❌: std::vector<std::unique_ptr<RWMol>> res;
    // RDKit❗❌: if (nFrags == 1) {
    // RDKit❗❌:   res.emplace_back(new RWMol(mol));
    // RDKit❗❌:   if (fragsMolAtomMapping) {
    // RDKit❗❌:     INT_VECT comp;
    // RDKit❗❌:     for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:       comp.push_back(idx);
    // RDKit❗❌:     }
    // RDKit❗❌:     (*fragsMolAtomMapping).push_back(comp);
    // RDKit❗❌:   }
    // RDKit❗❌: } else {
    // RDKit❗❌:   res.reserve(nFrags);
    // RDKit❗❌:   for (int i = 0; i < nFrags; ++i) {
    // RDKit❗❌:     boost::dynamic_bitset<> atomsInFrag(mol.getNumAtoms());
    // RDKit❗❌:     INT_VECT comp;
    // RDKit❗❌:     for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:       if ((*frags)[idx] == i) {
    // RDKit❗❌:         comp.push_back(idx);
    // RDKit❗❌:         atomsInFrag.set(idx);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (comp.size() == 1 ||
    // RDKit❗❌:         (nFrags > 3 && !fragmentHasChallengingFeatures(comp, atomsInFrag))) {
    // RDKit❗❌:       SubsetOptions opts{.sanitize = sanitizeFrags,
    // RDKit❗❌:                          .clearComputedProps = true,
    // RDKit❗❌:                          .copyCoordinates = copyConformers,
    // RDKit❗❌:                          .method = SubsetMethod::BONDS_BETWEEN_ATOMS};
    // RDKit❗❌:       std::vector<unsigned int> atoms{comp.begin(), comp.end()};
    // RDKit❗❌:       SubsetInfo info;
    // RDKit❗❌:       auto submol = copyMolSubset(mol, atoms, info, opts);
    // RDKit❗❌:       res.push_back(std::move(submol));
    // RDKit❗❌:     } else {
    // RDKit❗❌:       res.emplace_back(new RWMol(mol));
    // RDKit❗❌:       auto &frag = res.back();
    // RDKit❗❌:       frag->beginBatchEdit();
    // RDKit❗❌:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:         if (!atomsInFrag[idx]) {
    // RDKit❗❌:           frag->removeAtom(idx);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       frag->commitBatchEdit();
    // RDKit❗❌:     }
    // RDKit❗❌:     if (fragsMolAtomMapping) {
    // RDKit❗❌:       (*fragsMolAtomMapping).push_back(comp);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION MolOps::getTheFrags ordered component construction
    // Behavior: source labels and component order are retained, and component
    // atoms are projected in ascending source-row order. Singleton components
    // short-circuit to subset copying; other fast fragments require more than
    // three components and no challenging feature. The path helper sanitizes
    // before the loop advances; F15 slow copies defer sanitation to collection
    // postprocessing. The source outer sanitation loop is not performed here.
    // Complexity: component labeling is O(V+E), then source-shaped per-fragment
    // atom scans and masks are O(FV). The validated model subset and batch paths
    // carry known additional allocations compared with RDKit's in-place rows.
    let fragment_count = connected.components.len();
    if fragment_count == 0 {
        return Ok(Vec::new());
    }
    if fragment_count == 1 {
        let component_atoms = connected.components.into_iter().next().unwrap();
        let mut copy = copy_single_full_molecule_component_with_view(
            source_topology,
            source_coordinates,
            source_properties,
        );
        copy.source_metadata = clone_source_fragment_metadata(source_topology, source_metadata)
            .map_err(|source| OrderedFragmentBuildError::SlowFullCopy {
                component_index: 0,
                source,
            })?;
        if let Some(output) = component_output.as_deref_mut() {
            output.push(
                component_atoms
                    .iter()
                    .map(|atom| atom.index() as i32)
                    .collect(),
            );
        }
        return Ok(vec![OrderedFragmentCopy {
            component_atoms,
            copy,
        }]);
    }

    let mut fragments = Vec::with_capacity(fragment_count);
    for component_index in 0..fragment_count {
        let mut atoms_in_fragment = SelectionMask::default();
        atoms_in_fragment.resize(source_topology.atoms.len());
        let mut component_atoms = Vec::new();
        for (atom_index, label) in connected.atom_to_component.iter().copied().enumerate() {
            if label == component_index {
                let atom = AtomId::new(atom_index);
                component_atoms.push(atom);
                atoms_in_fragment.set(atom_index);
            }
        }

        let copied = if component_atoms.len() == 1
            || (fragment_count > 3
                && !fragment_has_challenging_features(
                    source_topology,
                    &component_atoms,
                    &atoms_in_fragment,
                )) {
            let subset = copy_mol_subset_atom_path_with_view(
                source_topology,
                source_coordinates,
                &component_atoms,
                sanitize,
                copy_conformers,
            )
            .map_err(|source| OrderedFragmentBuildError::FastSubset {
                component_index,
                source,
            })?;
            FullCopyComponent {
                topology: subset.topology,
                coordinates: subset.coordinates,
                molecule_properties: subset.molecule_properties,
                mapping: subset.mapping,
                source_metadata: fresh_subset_source_metadata(source_metadata),
            }
        } else {
            let copied = if source_boundary {
                copy_full_molecule_remove_atoms_outside_component_with_source_metadata(
                    source_topology,
                    source_coordinates,
                    source_properties,
                    &atoms_in_fragment,
                    source_metadata,
                )
            } else {
                copy_full_molecule_remove_atoms_outside_component_with_view(
                    source_topology,
                    source_coordinates,
                    source_properties,
                    &atoms_in_fragment,
                )
            };
            copied.map_err(|source| OrderedFragmentBuildError::SlowFullCopy {
                component_index,
                source,
            })?
        };
        if let Some(output) = component_output.as_deref_mut() {
            output.push(
                component_atoms
                    .iter()
                    .map(|atom| atom.index() as i32)
                    .collect(),
            );
        }
        fragments.push(OrderedFragmentCopy {
            component_atoms,
            copy: copied,
        });
    }
    Ok(fragments)
}

/// Source shared-ownership fragment overload over detached domain values.
#[doc(hidden)]
pub fn get_shared_molecule_fragments_with_source_outputs(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    source_properties: &MoleculeProperties,
    sanitize_fragments: bool,
    copy_conformers: bool,
    label_output: Option<&mut Vec<i32>>,
    component_output: Option<&mut Vec<Vec<i32>>>,
    source_metadata: FragmentSourceMetadataView<'_>,
) -> Result<Vec<std::sync::Arc<MoleculeFragment>>, MoleculeFragmentsError> {
    // RDKit❗✔️: std::vector<ROMOL_SPTR> getMolFrags(const ROMol &mol, bool sanitizeFrags,
    // RDKit❗✔️:                                     INT_VECT *frags,
    // RDKit❗✔️:                                     VECT_INT_VECT *fragsMolAtomMapping,
    // RDKit❗✔️:                                     bool copyConformers) {
    // RDKit❗✔️:   auto upFrags = getTheFrags(mol, sanitizeFrags, frags, fragsMolAtomMapping,
    // RDKit❗✔️:                              copyConformers);
    // RDKit❗✔️:   std::vector<boost::shared_ptr<ROMol>> finalRes;
    // RDKit❗✔️:   for (auto &r : upFrags) {
    // RDKit❗✔️:     finalRes.emplace_back(r.get());
    // RDKit❗✔️:     r.release();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return finalRes;
    // RDKit❗✔️: }
    // Behavior: build every unique detached fragment before transferring each
    // into shared ownership. Graph, coordinates, properties, and mappings are
    // moved, not cloned; all kernel errors and partial output updates propagate.
    // Native allocation exceptions/pointer deletion and Rust allocation failure
    // remain distinct owner-state capabilities, hence the behavior marker.
    // Complexity: one result-vector pass and one shared allocation per fragment,
    // as in the native shared control-block construction; Arc moves the owned
    // detached value and does not clone its vectors or underlying blocks.
    let fragments = get_molecule_fragments_with_source_outputs(
        source_topology,
        source_coordinates,
        source_properties,
        sanitize_fragments,
        copy_conformers,
        label_output,
        component_output,
        source_metadata,
    )?;
    Ok(fragments.into_iter().map(std::sync::Arc::new).collect())
}

/// Source owned-output fragment overload over detached domain values.
#[doc(hidden)]
pub fn assign_molecule_fragments_with_source_outputs(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    source_properties: &MoleculeProperties,
    fragments_output: &mut Vec<MoleculeFragment>,
    sanitize_fragments: bool,
    copy_conformers: bool,
    label_output: Option<&mut Vec<i32>>,
    component_output: Option<&mut Vec<Vec<i32>>>,
    source_metadata: FragmentSourceMetadataView<'_>,
) -> Result<u32, MoleculeFragmentsError> {
    // RDKit❗✔️: unsigned int getMolFrags(const ROMol &mol,
    // RDKit❗✔️:                          std::vector<std::unique_ptr<ROMol>> &molFrags,
    // RDKit❗✔️:                          bool sanitizeFrags, std::vector<int> *frags,
    // RDKit❗✔️:                          std::vector<std::vector<int>> *fragsMolAtomMapping,
    // RDKit❗✔️:                          bool copyConformers) {
    // RDKit❗✔️:   molFrags = getTheFrags(mol, sanitizeFrags, frags, fragsMolAtomMapping,
    // RDKit❗✔️:                          copyConformers);
    // RDKit❗✔️:   return rdcast<unsigned int>(molFrags.size());
    // RDKit❗✔️: }
    // Behavior: evaluate the complete kernel before replacing the caller's
    // prior owned fragments. On error those fragments stay intact while actual
    // label/map updates performed by the kernel remain observable. Assignment
    // moves the complete vector before reading its source unsigned count.
    // RDKit❗✔️: #define rdcast static_cast
    // Ordinary pinned native builds use unsigned truncation for rdcast; native
    // RDDEBUG numeric_cast exceptions are a separate unmodeled build policy.
    // Complexity: no topology/coordinate/property clone; vector assignment
    // drops the old owned results and transfers the new vector in O(1), aside
    // from the same old-fragment destruction costs as native unique ownership.
    let fragments = get_molecule_fragments_with_source_outputs(
        source_topology,
        source_coordinates,
        source_properties,
        sanitize_fragments,
        copy_conformers,
        label_output,
        component_output,
        source_metadata,
    )?;
    *fragments_output = fragments;
    Ok(fragments_output.len() as u32)
}

/// Build connected components in source order from borrowed detached values.
///
/// The returned fragment values own their topology, coordinates, properties,
/// component atom rows, and validated topology mappings.
pub fn get_molecule_fragments(
    source_topology: &TopologyBlock,
    source_coordinates: &CoordinateBlock,
    source_properties: &MoleculeProperties,
    sanitize_fragments: bool,
    copy_conformers: bool,
) -> Result<Vec<MoleculeFragment>, MoleculeFragmentsError> {
    let coordinate_view = FragmentCoordinateView::from_coordinate_block(source_coordinates);
    get_molecule_fragments_with_coordinate_view(
        source_topology,
        &coordinate_view,
        source_properties,
        sanitize_fragments,
        copy_conformers,
    )
}

/// Build source-ordered fragment copies from a borrowed coordinate view.
///
/// This sibling algorithm entry accepts the selected kernel rows without
/// materializing another input coordinate block.
pub fn get_molecule_fragments_with_coordinate_view(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    source_properties: &MoleculeProperties,
    sanitize_fragments: bool,
    copy_conformers: bool,
) -> Result<Vec<MoleculeFragment>, MoleculeFragmentsError> {
    get_molecule_fragments_impl(
        source_topology,
        source_coordinates,
        source_properties,
        sanitize_fragments,
        copy_conformers,
        None,
        None,
        FragmentSourceMetadataView::unmodeled(),
        false,
    )
}

/// Complete source fragment dispatch with actual optional output buffers and metadata.
/// Independently unmodeled ring/bookmark fields remain `None`, never invented defaults.
#[doc(hidden)]
pub fn get_molecule_fragments_with_source_outputs(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    source_properties: &MoleculeProperties,
    sanitize_fragments: bool,
    copy_conformers: bool,
    label_output: Option<&mut Vec<i32>>,
    component_output: Option<&mut Vec<Vec<i32>>>,
    source_metadata: FragmentSourceMetadataView<'_>,
) -> Result<Vec<MoleculeFragment>, MoleculeFragmentsError> {
    get_molecule_fragments_impl(
        source_topology,
        source_coordinates,
        source_properties,
        sanitize_fragments,
        copy_conformers,
        label_output,
        component_output,
        source_metadata,
        true,
    )
}

fn get_molecule_fragments_impl(
    source_topology: &TopologyBlock,
    source_coordinates: &FragmentCoordinateView<'_>,
    source_properties: &MoleculeProperties,
    sanitize_fragments: bool,
    copy_conformers: bool,
    label_output: Option<&mut Vec<i32>>,
    component_output: Option<&mut Vec<Vec<i32>>>,
    source_metadata: FragmentSourceMetadataView<'_>,
    source_boundary: bool,
) -> Result<Vec<MoleculeFragment>, MoleculeFragmentsError> {
    // RDKit❗❌: std::vector<std::unique_ptr<ROMol>> getTheFrags(
    // RDKit❗❌:     const ROMol &mol, bool sanitizeFrags, INT_VECT *frags,
    // RDKit❗❌:     VECT_INT_VECT *fragsMolAtomMapping, bool copyConformers) {
    // RDKit❗❌:   std::unique_ptr<INT_VECT> mappingStorage;
    // RDKit❗❌:   if (!frags) {
    // RDKit❗❌:     mappingStorage.reset(new INT_VECT);
    // RDKit❗❌:     frags = mappingStorage.get();
    // RDKit❗❌:   }
    // RDKit❗❌:   int nFrags = getMolFrags(mol, *frags);
    // RDKit❗❌:   std::vector<std::unique_ptr<RWMol>> res;
    // RDKit❗❌:
    // RDKit❗❌:   if (nFrags == 1) {
    // RDKit❗❌:     res.emplace_back(new RWMol(mol));
    // RDKit❗❌:     if (fragsMolAtomMapping) {
    // RDKit❗❌:       INT_VECT comp;
    // RDKit❗❌:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:         comp.push_back(idx);
    // RDKit❗❌:       }
    // RDKit❗❌:       (*fragsMolAtomMapping).push_back(comp);
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     res.reserve(nFrags);
    // RDKit❗❌:     for (int i = 0; i < nFrags; ++i) {
    // RDKit❗❌:       boost::dynamic_bitset<> atomsInFrag(mol.getNumAtoms());
    // RDKit❗❌:       INT_VECT comp;
    // RDKit❗❌:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:         if ((*frags)[idx] == i) {
    // RDKit❗❌:           comp.push_back(idx);
    // RDKit❗❌:           atomsInFrag.set(idx);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       auto fragmentHasChallengingFeatures =
    // RDKit❗❌:           [&](const INT_VECT &comp,
    // RDKit❗❌:               const boost::dynamic_bitset<> &atomsInFrag) -> bool {
    // RDKit❗❌:         for (auto idx : comp) {
    // RDKit❗❌:           // check for atoms with stereochem:
    // RDKit❗❌:           const auto atom = mol.getAtomWithIdx(idx);
    // RDKit❗❌:           if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED &&
    // RDKit❗❌:               atom->getChiralTag() != Atom::ChiralType::CHI_OTHER) {
    // RDKit❗❌:             return true;
    // RDKit❗❌:           }
    // RDKit❗❌:           for (auto bnd : mol.atomBonds(atom)) {
    // RDKit❗❌:             if (atomsInFrag[bnd->getOtherAtomIdx(idx)]) {
    // RDKit❗❌:               if (bnd->getStereo() != Bond::BondStereo::STEREONONE &&
    // RDKit❗❌:                   bnd->getStereo() != Bond::BondStereo::STEREOANY) {
    // RDKit❗❌:                 return true;
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         for (auto sgroup : getSubstanceGroups(mol)) {
    // RDKit❗❌:           for (auto aid : sgroup.getAtoms()) {
    // RDKit❗❌:             if (atomsInFrag[aid]) {
    // RDKit❗❌:               return true;
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:           for (auto aid : sgroup.getParentAtoms()) {
    // RDKit❗❌:             if (atomsInFrag[aid]) {
    // RDKit❗❌:               return true;
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         for (auto stereoGroup : mol.getStereoGroups()) {
    // RDKit❗❌:           // doesn't seem like this should be necessary, but in case
    // RDKit❗❌:           // we ever need stereogroups where the atoms aren't marked
    // RDKit❗❌:           // with stereo...
    // RDKit❗❌:           for (auto atom : stereoGroup.getAtoms()) {
    // RDKit❗❌:             if (atomsInFrag[atom->getIdx()]) {
    // RDKit❗❌:               return true;
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:           // same check for stereo groups involving bonds:
    // RDKit❗❌:           for (auto bond : stereoGroup.getBonds()) {
    // RDKit❗❌:             if (atomsInFrag[bond->getBeginAtomIdx()] &&
    // RDKit❗❌:                 atomsInFrag[bond->getEndAtomIdx()]) {
    // RDKit❗❌:               return true;
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         return false;
    // RDKit❗❌:       };
    // RDKit❗❌:       if (comp.size() == 1 ||
    // RDKit❗❌:           (nFrags > 3 && !fragmentHasChallengingFeatures(comp, atomsInFrag))) {
    // RDKit❗❌:         // special case for a small, simple fragments when a bunch of fragments
    // RDKit❗❌:         // are present. The check on the number of fragments is purely
    // RDKit❗❌:         // empirical. This is mainly intended to catch situations like proteins
    // RDKit❗❌:         // where you have a bunch of single-atom fragments (waters); the
    // RDKit❗❌:         // standard approach below ends up being horribly inefficient there
    // RDKit❗❌:         SubsetOptions opts{.sanitize = sanitizeFrags,
    // RDKit❗❌:                            .clearComputedProps = true,
    // RDKit❗❌:                            .copyCoordinates = copyConformers,
    // RDKit❗❌:                            .method = SubsetMethod::BONDS_BETWEEN_ATOMS};
    // RDKit❗❌:         std::vector<unsigned int> atoms{comp.begin(), comp.end()};
    // RDKit❗❌:         SubsetInfo info;
    // RDKit❗❌:         auto submol = copyMolSubset(mol, atoms, info, opts);
    // RDKit❗❌:         res.push_back(std::move(submol));
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res.emplace_back(new RWMol(mol));
    // RDKit❗❌:         auto &frag = res.back();
    // RDKit❗❌:
    // RDKit❗❌:         frag->beginBatchEdit();
    // RDKit❗❌:         for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:           if (!atomsInFrag[idx]) {
    // RDKit❗❌:             frag->removeAtom(idx);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         frag->commitBatchEdit();
    // RDKit❗❌:       }
    // RDKit❗❌:       if (fragsMolAtomMapping) {
    // RDKit❗❌:         (*fragsMolAtomMapping).push_back(comp);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!copyConformers) {
    // RDKit❗❌:     for (auto &frag : res) {
    // RDKit❗❌:       frag->clearConformers();
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (sanitizeFrags) {
    // RDKit❗❌:     for (auto &frag : res) {
    // RDKit❗❌:       sanitizeMol(*frag);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<std::unique_ptr<ROMol>> finalRes;
    // RDKit❗❌:   for (auto &r : res) {
    // RDKit❗❌:     finalRes.emplace_back(r.get());
    // RDKit❗❌:     r.release();
    // RDKit❗❌:   }
    // RDKit❗❌:   return finalRes;
    // RDKit❗❌: }

    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags postprocessing
    // RDKit❗❌:   if (!copyConformers) {
    // RDKit❗❌:     for (auto &frag : res) {
    // RDKit❗❌:       frag->clearConformers();
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (sanitizeFrags) {
    // RDKit❗❌:     for (auto &frag : res) {
    // RDKit❗❌:       sanitizeMol(*frag);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: void sanitizeMol(RWMol &mol) {
    // RDKit❗❌:   unsigned int failedOp = 0;
    // RDKit❗❌:   sanitizeMol(mol, failedOp, SANITIZE_ALL);
    // RDKit❗❌: }
    // RDKit❗❌: void sanitizeMol(RWMol &mol, unsigned int &operationThatFailed,
    // RDKit❗❌:                  unsigned int sanitizeOps) {
    // RDKit❗❌:   mol.clearComputedProps();
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<std::unique_ptr<ROMol>> finalRes;
    // RDKit❗❌:   for (auto &r : res) {
    // RDKit❗❌:     finalRes.emplace_back(r.get());
    // RDKit❗❌:     r.release();
    // RDKit❗❌:   }
    // RDKit❗❌:   return finalRes;
    // RDKit❗❌: }
    // RDKit❗❌: std::vector<ROMOL_SPTR> getMolFrags(
    // RDKit❗❌:     const ROMol &mol, bool sanitizeFrags, INT_VECT *frags,
    // RDKit❗❌:     VECT_INT_VECT *fragsMolAtomMapping, bool copyConformers) {
    // RDKit❗❌:   auto upFrags = getTheFrags(mol, sanitizeFrags, frags,
    // RDKit❗❌:                              fragsMolAtomMapping, copyConformers);
    // RDKit❗❌:   std::vector<boost::shared_ptr<ROMol>> finalRes;
    // RDKit❗❌:   for (auto &r : upFrags) {
    // RDKit❗❌:     finalRes.emplace_back(r.get());
    // RDKit❗❌:     r.release();
    // RDKit❗❌:   }
    // RDKit❗❌:   return finalRes;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION MolOps::getTheFrags postprocessing
    // Behavior: all component copies finish before outer conformer clearing;
    // default sanitation then runs in fragment order, after which Rust moves
    // each completed detached value into the returned owner vector. Fast
    // subsets may already have had their source-required inner sanitize pass.
    // Clearing removes conformer rows while retaining the separate modeled
    // source-coordinate-dimension provenance field.
    // Complexity: clearing is linear in copied conformer records; optional
    // default sanitation visits each completed fragment once more. Moving the
    // already-built vector adds no second topology or coordinate clone.
    let connected = crate::paths::connected_components(source_topology)
        .map_err(OrderedFragmentBuildError::from)
        .map_err(MoleculeFragmentsError::from)?;
    if let Some(output) = label_output {
        output.clear();
        output.extend(
            connected
                .atom_to_component
                .iter()
                .map(|label| *label as i32),
        );
    }
    let mut fragments = build_ordered_fragment_copies_from_components(
        source_topology,
        source_coordinates,
        source_properties,
        sanitize_fragments,
        copy_conformers,
        connected,
        component_output,
        source_metadata,
        source_boundary,
    )
    .map_err(MoleculeFragmentsError::from)?;

    if !copy_conformers {
        for fragment in &mut fragments {
            fragment.copy.coordinates.clear_2d_conformers();
            fragment.copy.coordinates.clear_3d_conformers();
        }
    }

    if sanitize_fragments {
        for (component_index, fragment) in fragments.iter_mut().enumerate() {
            fragment
                .copy
                .molecule_properties
                .clear_computed_props()
                .map_err(|source| MoleculeFragmentsError {
                    failure: MoleculeFragmentsFailure::FinalComputedProperties {
                        component_index,
                        source,
                    },
                })?;
            let sanitized = sanitize_topology(&fragment.copy.topology, &SanitizeParams::default())
                .map_err(|source| MoleculeFragmentsError {
                    failure: MoleculeFragmentsFailure::FinalSanitize {
                        component_index,
                        source,
                    },
                })?;
            fragment.copy.topology = sanitized.topology;
            if fragment.copy.source_metadata.rings.is_some() {
                fragment.copy.source_metadata.rings = Some(match sanitized.final_rings {
                    Some(rings) => rings,
                    None => source_uninitialized_ring_info(),
                });
            }
        }
    }

    Ok(fragments
        .into_iter()
        .map(|fragment| MoleculeFragment {
            component_atoms: fragment.component_atoms,
            copy: fragment.copy,
        })
        .collect())
}

#[cfg(test)]
fn fixed_property_text(value: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(value.as_bytes()).expect("fixed fixture text is UTF8")
}

#[cfg(test)]
mod cf3d_frag_f06_tests {
    use super::{FragmentSubsetInfo, SelectionMask, copy_selected_stereo_groups};
    use cosmolkit_model::{AtomId, BondId, StereoGroup, StereoGroupKind};

    #[derive(Clone, Copy)]
    enum CategorySelection {
        Empty,
        NoneSelected,
        SomeSelected,
        AllSelected,
    }

    struct Case {
        atoms: CategorySelection,
        bonds: CategorySelection,
        retain_nonempty_groups: bool,
        expected_atoms: &'static [usize],
        expected_bonds: &'static [usize],
    }

    fn source_members(category: CategorySelection) -> Vec<usize> {
        match category {
            CategorySelection::Empty => vec![],
            CategorySelection::NoneSelected
            | CategorySelection::SomeSelected
            | CategorySelection::AllSelected => {
                vec![2, 0]
            }
        }
    }

    fn selected_members(category: CategorySelection) -> &'static [usize] {
        match category {
            CategorySelection::Empty | CategorySelection::NoneSelected => &[],
            CategorySelection::SomeSelected => &[2],
            CategorySelection::AllSelected => &[2, 0],
        }
    }

    fn group(
        kind: StereoGroupKind,
        atoms: Vec<AtomId>,
        bonds: Vec<BondId>,
        read_id: Option<u32>,
        write_id: u32,
    ) -> StereoGroup {
        let group = StereoGroup::new(kind, atoms, bonds);
        let group = match read_id {
            Some(id) => group.with_id(id),
            None => group,
        };
        group.with_write_id(write_id)
    }

    fn ids(indices: &[usize]) -> Vec<AtomId> {
        indices.iter().copied().map(AtomId::new).collect()
    }

    fn bond_ids(indices: &[usize]) -> Vec<BondId> {
        indices.iter().copied().map(BondId::new).collect()
    }

    #[test]
    fn cf3d_frag_f06_stereo_group_copy_covers_category_matrix_and_identity_order() {
        use CategorySelection::{AllSelected, Empty, NoneSelected, SomeSelected};

        let cases = [
            Case {
                atoms: Empty,
                bonds: Empty,
                retain_nonempty_groups: true,
                expected_atoms: &[],
                expected_bonds: &[],
            },
            Case {
                atoms: Empty,
                bonds: NoneSelected,
                retain_nonempty_groups: false,
                expected_atoms: &[],
                expected_bonds: &[],
            },
            Case {
                atoms: Empty,
                bonds: SomeSelected,
                retain_nonempty_groups: true,
                expected_atoms: &[],
                expected_bonds: &[1],
            },
            Case {
                atoms: Empty,
                bonds: AllSelected,
                retain_nonempty_groups: true,
                expected_atoms: &[],
                expected_bonds: &[1, 2],
            },
            Case {
                atoms: NoneSelected,
                bonds: Empty,
                retain_nonempty_groups: false,
                expected_atoms: &[],
                expected_bonds: &[],
            },
            Case {
                atoms: NoneSelected,
                bonds: NoneSelected,
                retain_nonempty_groups: false,
                expected_atoms: &[],
                expected_bonds: &[],
            },
            Case {
                atoms: NoneSelected,
                bonds: SomeSelected,
                retain_nonempty_groups: false,
                expected_atoms: &[],
                expected_bonds: &[],
            },
            Case {
                atoms: NoneSelected,
                bonds: AllSelected,
                retain_nonempty_groups: false,
                expected_atoms: &[],
                expected_bonds: &[],
            },
            Case {
                atoms: SomeSelected,
                bonds: Empty,
                retain_nonempty_groups: true,
                expected_atoms: &[3],
                expected_bonds: &[],
            },
            Case {
                atoms: SomeSelected,
                bonds: NoneSelected,
                retain_nonempty_groups: false,
                expected_atoms: &[],
                expected_bonds: &[],
            },
            Case {
                atoms: SomeSelected,
                bonds: SomeSelected,
                retain_nonempty_groups: true,
                expected_atoms: &[3],
                expected_bonds: &[1],
            },
            Case {
                atoms: SomeSelected,
                bonds: AllSelected,
                retain_nonempty_groups: true,
                expected_atoms: &[3],
                expected_bonds: &[1, 2],
            },
            Case {
                atoms: AllSelected,
                bonds: Empty,
                retain_nonempty_groups: true,
                expected_atoms: &[3, 1],
                expected_bonds: &[],
            },
            Case {
                atoms: AllSelected,
                bonds: NoneSelected,
                retain_nonempty_groups: false,
                expected_atoms: &[],
                expected_bonds: &[],
            },
            Case {
                atoms: AllSelected,
                bonds: SomeSelected,
                retain_nonempty_groups: true,
                expected_atoms: &[3, 1],
                expected_bonds: &[1],
            },
            Case {
                atoms: AllSelected,
                bonds: AllSelected,
                retain_nonempty_groups: true,
                expected_atoms: &[3, 1],
                expected_bonds: &[1, 2],
            },
        ];

        for (case_index, case) in cases.into_iter().enumerate() {
            let source_atoms = source_members(case.atoms)
                .into_iter()
                .map(AtomId::new)
                .collect::<Vec<_>>();
            let source_bonds = source_members(case.bonds)
                .into_iter()
                .map(BondId::new)
                .collect::<Vec<_>>();
            let reference = vec![
                group(
                    StereoGroupKind::Or,
                    source_atoms.clone(),
                    source_bonds.clone(),
                    Some(17),
                    9,
                ),
                group(StereoGroupKind::Absolute, vec![], vec![], None, 0),
                group(StereoGroupKind::And, source_atoms, source_bonds, Some(0), 4),
            ];
            let original_reference = reference.clone();

            let mut selection = FragmentSubsetInfo::default();
            selection.selected_atoms = SelectionMask::default();
            selection.selected_atoms.resize(4);
            selection.selected_bonds = SelectionMask::default();
            selection.selected_bonds.resize(3);
            for source_index in selected_members(case.atoms) {
                selection.selected_atoms.set(*source_index);
                let mapped = match source_index {
                    2 => 3,
                    0 => 1,
                    _ => unreachable!("the fixed F06 matrix uses only source rows 2 and 0"),
                };
                selection
                    .atom_mapping
                    .insert(AtomId::new(*source_index), AtomId::new(mapped));
            }
            for source_index in selected_members(case.bonds) {
                selection.selected_bonds.set(*source_index);
                let mapped = match source_index {
                    2 => 1,
                    0 => 2,
                    _ => unreachable!("the fixed F06 matrix uses only source rows 2 and 0"),
                };
                selection
                    .bond_mapping
                    .insert(BondId::new(*source_index), BondId::new(mapped));
            }

            let actual = copy_selected_stereo_groups(&reference, &selection);
            let mut expected = Vec::new();
            if case.retain_nonempty_groups {
                expected.push(group(
                    StereoGroupKind::Or,
                    ids(case.expected_atoms),
                    bond_ids(case.expected_bonds),
                    Some(17),
                    9,
                ));
            }
            expected.push(group(StereoGroupKind::Absolute, vec![], vec![], None, 0));
            if case.retain_nonempty_groups {
                expected.push(group(
                    StereoGroupKind::And,
                    ids(case.expected_atoms),
                    bond_ids(case.expected_bonds),
                    Some(0),
                    4,
                ));
            }

            assert_eq!(actual, expected, "fixed category case {case_index}");
            assert_eq!(
                reference, original_reference,
                "input case {case_index} mutated"
            );
        }
    }
}

#[cfg(test)]
mod cf3d_frag_f07_tests {
    use super::{FragmentCoordinateView, copy_full_molecule_coordinates, copy_subset_coordinates};
    use cosmolkit_model::{AtomId, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension};
    use std::collections::BTreeMap;

    #[test]
    fn cf3d_frag_f07_empty_and_single_conformer_cases_preserve_source_rules() {
        let empty = CoordinateBlock::default();
        let empty_view = FragmentCoordinateView::from_coordinate_block(&empty);
        assert_eq!(copy_subset_coordinates(&empty, &BTreeMap::new()), empty);
        assert_eq!(copy_full_molecule_coordinates(&empty_view), empty);

        let one = CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(41, vec![[0.0, 0.0], [1.0, 2.0], [3.0, 4.0]])
                    .with_prop("subset-only", "drop"),
            ],
            conformers_3d: vec![
                Conformer3D::new(
                    41,
                    vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                    false,
                )
                .with_prop("subset-only", "drop"),
            ],
            source_coordinate_dim: Some(CoordinateDimension::TwoD),
            source_conformer_order: None,
        };
        let original_one = one.clone();
        let one_view = FragmentCoordinateView::from_coordinate_block(&one);
        let atom_mapping = BTreeMap::from([(AtomId::new(2), AtomId::new(0))]);
        let expected_subset = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(41, vec![[3.0, 4.0]])],
            conformers_3d: vec![Conformer3D::new(41, vec![[4.0, 5.0, 6.0]], false)],
            source_coordinate_dim: Some(CoordinateDimension::TwoD),
            source_conformer_order: None,
        };

        assert_eq!(
            copy_subset_coordinates(&one, &atom_mapping),
            expected_subset
        );
        assert_eq!(copy_full_molecule_coordinates(&one_view), one);
        assert_eq!(one, original_one);
    }

    #[test]
    fn cf3d_frag_f07_subset_multiple_conformers_maps_gapped_rows_and_drops_props() {
        let source = CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(
                    19,
                    vec![
                        [0.0, 0.0],
                        [1.0, 10.0],
                        [2.0, 20.0],
                        [3.0, 30.0],
                        [4.0, 40.0],
                        [5.0, 50.0],
                    ],
                )
                .with_prop("kind", "2d-first"),
                Conformer2D::new(
                    103,
                    vec![
                        [10.0, 100.0],
                        [11.0, 110.0],
                        [12.0, 120.0],
                        [13.0, 130.0],
                        [14.0, 140.0],
                        [15.0, 150.0],
                    ],
                )
                .with_prop("kind", "2d-second"),
            ],
            conformers_3d: vec![
                Conformer3D::new(
                    19,
                    vec![
                        [0.0, 0.0, 0.0],
                        [1.0, 10.0, 100.0],
                        [2.0, 20.0, 200.0],
                        [3.0, 30.0, 300.0],
                        [4.0, 40.0, 400.0],
                        [5.0, 50.0, 500.0],
                    ],
                    true,
                )
                .with_prop("kind", "3d-first"),
                Conformer3D::new(
                    207,
                    vec![
                        [10.0, 100.0, 1000.0],
                        [11.0, 110.0, 1100.0],
                        [12.0, 120.0, 1200.0],
                        [13.0, 130.0, 1300.0],
                        [14.0, 140.0, 1400.0],
                        [15.0, 150.0, 1500.0],
                    ],
                    false,
                )
                .with_prop("kind", "3d-second"),
            ],
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            source_conformer_order: None,
        };
        let original_source = source.clone();
        let atom_mapping = BTreeMap::from([
            (AtomId::new(1), AtomId::new(0)),
            (AtomId::new(4), AtomId::new(1)),
            (AtomId::new(5), AtomId::new(2)),
        ]);
        let expected = CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(19, vec![[1.0, 10.0], [4.0, 40.0], [5.0, 50.0]]),
                Conformer2D::new(103, vec![[11.0, 110.0], [14.0, 140.0], [15.0, 150.0]]),
            ],
            conformers_3d: vec![
                Conformer3D::new(
                    19,
                    vec![[1.0, 10.0, 100.0], [4.0, 40.0, 400.0], [5.0, 50.0, 500.0]],
                    true,
                ),
                Conformer3D::new(
                    207,
                    vec![
                        [11.0, 110.0, 1100.0],
                        [14.0, 140.0, 1400.0],
                        [15.0, 150.0, 1500.0],
                    ],
                    false,
                ),
            ],
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            source_conformer_order: None,
        };

        assert_eq!(copy_subset_coordinates(&source, &atom_mapping), expected);
        assert_eq!(source, original_source);
    }

    #[test]
    fn cf3d_frag_f07_full_copy_keeps_order_ids_flags_props_and_owned_rows() {
        let source = CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(72, vec![[1.0, 2.0], [3.0, 4.0]])
                    .with_prop("copy", "preserved-2d"),
                Conformer2D::new(8, vec![[5.0, 6.0], [7.0, 8.0]])
                    .with_prop("copy", "preserved-second-2d"),
            ],
            conformers_3d: vec![
                Conformer3D::new(72, vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], true)
                    .with_prop("copy", "preserved-3d"),
                Conformer3D::new(9, vec![[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]], false)
                    .with_prop("copy", "preserved-second-3d"),
            ],
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            source_conformer_order: None,
        };
        let original_source = source.clone();
        let source_view = FragmentCoordinateView::from_coordinate_block(&source);

        let mut copied = copy_full_molecule_coordinates(&source_view);
        assert_eq!(copied, source);
        copied.conformers_2d[0].coordinates_mut()[0] = [-1.0, -2.0];
        assert_eq!(source, original_source);
        assert_eq!(source.conformers_2d[0].coordinates()[0], [1.0, 2.0]);
    }
}

#[cfg(test)]
mod cf3d_frag_f08_tests {
    use super::sanitize_subset_if_requested;
    use crate::{
        KekulizeError, PropertyCacheError, SanitizeError, SanitizeStage, ValenceError, ValencePhase,
    };
    use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
    use cosmolkit_types::{BondOrder, Element};

    fn topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
        let atoms = atom_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect();
        let bonds = bond_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    fn alternating_benzene() -> TopologyBlock {
        let atom_specs = vec![AtomSpec::new(Element::C); 6];
        let bond_specs = (0..6)
            .map(|index| {
                BondSpec::new(
                    AtomId::new(index),
                    AtomId::new((index + 1) % 6),
                    if index % 2 == 0 {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    },
                )
            })
            .collect();
        topology(atom_specs, bond_specs)
    }

    #[test]
    fn cf3d_frag_f08_disabled_skips_and_enabled_sanitizes_benzene() {
        let source = alternating_benzene();
        let original_source = source.clone();

        let mut disabled = source.clone();
        sanitize_subset_if_requested(&mut disabled, false).unwrap();
        assert_eq!(disabled, source);

        let mut enabled = source.clone();
        sanitize_subset_if_requested(&mut enabled, true).unwrap();
        assert_ne!(enabled, source);
        assert!(enabled.atoms.iter().all(Atom::is_aromatic));
        assert!(enabled.bonds.iter().all(Bond::is_aromatic));
        assert_eq!(source, original_source);
    }

    #[test]
    fn cf3d_frag_f08_preserves_valence_failure_and_input_for_both_options() {
        let mut atom_specs = vec![AtomSpec::new(Element::C)];
        atom_specs.extend((0..5).map(|_| AtomSpec::new(Element::H)));
        let bond_specs = (1..=5)
            .map(|hydrogen| BondSpec::new(AtomId::new(0), AtomId::new(hydrogen), BondOrder::Single))
            .collect();
        let source = topology(atom_specs, bond_specs);
        let original_source = source.clone();

        let mut disabled = source.clone();
        sanitize_subset_if_requested(&mut disabled, false).unwrap();
        assert_eq!(disabled, source);

        let mut enabled = source.clone();
        let error = sanitize_subset_if_requested(&mut enabled, true).unwrap_err();
        assert!(matches!(
            error,
            SanitizeError::Properties {
                stage: SanitizeStage::Properties,
                source: PropertyCacheError::Valence(ValenceError::InvalidValence {
                    atom,
                    atomic_number: 6,
                    phase: ValencePhase::Explicit,
                    ..
                }),
            } if atom == AtomId::new(0)
        ));
        assert_eq!(enabled, source);
        assert_eq!(source, original_source);
    }

    #[test]
    fn cf3d_frag_f08_preserves_kekulization_failure_and_input_for_both_options() {
        let source = topology(
            vec![
                AtomSpec::new(Element::C)
                    .with_aromatic(true)
                    .with_no_implicit(true),
            ],
            Vec::new(),
        );
        let original_source = source.clone();

        let mut disabled = source.clone();
        sanitize_subset_if_requested(&mut disabled, false).unwrap();
        assert_eq!(disabled, source);

        let mut enabled = source.clone();
        let error = sanitize_subset_if_requested(&mut enabled, true).unwrap_err();
        assert!(matches!(
            error,
            SanitizeError::Kekulize {
                stage: SanitizeStage::Kekulize,
                source: KekulizeError::AromaticAtomOutsideRing { atom },
            } if atom == AtomId::new(0)
        ));
        assert_eq!(enabled, source);
        assert_eq!(source, original_source);
    }
}

#[cfg(test)]
mod cf3d_frag_f09_tests {
    use super::{clear_subset_computed_props, sanitize_subset_if_requested};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, MoleculeProperties, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    fn ring_with_registered_properties() -> (TopologyBlock, MoleculeProperties) {
        let atom_specs = (0..6)
            .map(|index| {
                AtomSpec::new(Element::C)
                    .with_prop("ordinary-atom", format!("atom-{index}"))
                    .unwrap()
                    .with_prop("ring-cache-note", format!("atom-note-{index}"))
                    .unwrap()
                    .with_computed_prop("computed-ring-cache", format!("atom-cache-{index}"))
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let bond_specs = (0..6)
            .map(|index| {
                BondSpec::new(
                    AtomId::new(index),
                    AtomId::new((index + 1) % 6),
                    if index % 2 == 0 {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    },
                )
                .with_prop("ordinary-bond", format!("bond-{index}"))
                .unwrap()
                .with_prop("ring-cache-note", format!("bond-note-{index}"))
                .unwrap()
                .with_computed_prop("computed-ring-cache", format!("bond-cache-{index}"))
                .unwrap()
            })
            .collect::<Vec<_>>();
        let atoms = atom_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect();
        let bonds = bond_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
            .collect();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
        let molecule_properties = MoleculeProperties::default()
            .with_name("ring-cache-clear")
            .with_prop("ordinary-molecule", "retained")
            .unwrap()
            .with_prop("ring-cache-note", "ordinary value")
            .unwrap()
            .with_computed_prop("computed-ring-cache", "stale")
            .unwrap();

        (topology, molecule_properties)
    }

    #[test]
    fn cf3d_frag_f09_clears_registered_props_without_clearing_chemistry() {
        let (source_topology, source_properties) = ring_with_registered_properties();
        let original_topology = source_topology.clone();
        let original_properties = source_properties.clone();

        for sanitize in [false, true] {
            let mut topology = source_topology.clone();
            let mut molecule_properties = source_properties.clone();
            sanitize_subset_if_requested(&mut topology, sanitize).unwrap();

            let bond_orders_before_clear =
                topology.bonds.iter().map(Bond::order).collect::<Vec<_>>();
            if sanitize {
                assert!(topology.atoms.iter().all(Atom::is_aromatic));
                assert!(topology.bonds.iter().all(Bond::is_aromatic));
            } else {
                assert!(topology.atoms.iter().all(|atom| !atom.is_aromatic()));
                assert!(topology.bonds.iter().all(|bond| !bond.is_aromatic()));
                assert_eq!(
                    bond_orders_before_clear,
                    vec![
                        BondOrder::Double,
                        BondOrder::Single,
                        BondOrder::Double,
                        BondOrder::Single,
                        BondOrder::Double,
                        BondOrder::Single,
                    ]
                );
            }

            // The detached topology stores no persistent RingInfo cache.
            // Clearing computed carrier values must still preserve the
            // chemical ring flags and bond orders produced above.
            clear_subset_computed_props(&mut topology, &mut molecule_properties);

            assert_eq!(
                topology.bonds.iter().map(Bond::order).collect::<Vec<_>>(),
                bond_orders_before_clear
            );
            assert_eq!(
                topology
                    .atoms
                    .iter()
                    .map(Atom::is_aromatic)
                    .collect::<Vec<_>>(),
                vec![sanitize; 6]
            );
            assert_eq!(
                topology
                    .bonds
                    .iter()
                    .map(Bond::is_aromatic)
                    .collect::<Vec<_>>(),
                vec![sanitize; 6]
            );
            for (index, atom) in topology.atoms.iter().enumerate() {
                assert_eq!(
                    atom.prop("ordinary-atom"),
                    Some(&cosmolkit_model::PropertyValue::String(
                        (format!("atom-{index}")).into()
                    ))
                );
                assert_eq!(
                    atom.prop("ring-cache-note"),
                    Some(&cosmolkit_model::PropertyValue::String(
                        (format!("atom-note-{index}")).into()
                    ))
                );
                assert_eq!(atom.prop("computed-ring-cache"), None);
                assert!(!atom.is_prop_computed("computed-ring-cache").unwrap());
            }
            for (index, bond) in topology.bonds.iter().enumerate() {
                assert_eq!(
                    bond.prop("ordinary-bond"),
                    Some(&cosmolkit_model::PropertyValue::String(
                        (format!("bond-{index}")).into()
                    ))
                );
                assert_eq!(
                    bond.prop("ring-cache-note"),
                    Some(&cosmolkit_model::PropertyValue::String(
                        (format!("bond-note-{index}")).into()
                    ))
                );
                assert_eq!(bond.prop("computed-ring-cache"), None);
                assert!(!bond.is_prop_computed("computed-ring-cache").unwrap());
            }
            assert_eq!(
                molecule_properties.name().map(super::fixed_property_text),
                Some("ring-cache-clear")
            );
            assert_eq!(
                molecule_properties.prop("ordinary-molecule"),
                Some(&cosmolkit_model::PropertyValue::String("retained".into()))
            );
            assert_eq!(
                molecule_properties.prop("ring-cache-note"),
                Some(&cosmolkit_model::PropertyValue::String(
                    "ordinary value".into()
                ))
            );
            assert_eq!(molecule_properties.prop("computed-ring-cache"), None);
            assert!(
                molecule_properties
                    .computed_prop_names()
                    .unwrap()
                    .expect("cleared computed StringVector exists")
                    .is_empty()
            );
            assert_eq!(source_topology, original_topology);
            assert_eq!(source_properties, original_properties);
        }
    }
}

#[cfg(test)]
mod cf3d_frag_f02_tests {
    use super::{FragmentSubsetInfo, get_subset_info_for_atom_path};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, TopologyBlock,
    };
    use cosmolkit_types::Element;

    fn topology(atom_count: usize, edges: &[(usize, usize)]) -> TopologyBlock {
        let atoms = (0..atom_count)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed atom-path topology is valid")
    }

    fn assert_selection(
        selection: &FragmentSubsetInfo,
        atom_count: usize,
        bond_count: usize,
        expected_atoms: &[usize],
        expected_bonds: &[usize],
    ) {
        for atom in 0..atom_count {
            assert_eq!(
                selection.selected_atoms.contains(atom),
                expected_atoms.contains(&atom),
                "atom row {atom}"
            );
        }
        for bond in 0..bond_count {
            assert_eq!(
                selection.selected_bonds.contains(bond),
                expected_bonds.contains(&bond),
                "bond row {bond}"
            );
        }
        assert!(selection.atom_mapping.is_empty());
        assert!(selection.bond_mapping.is_empty());
    }

    #[test]
    fn cf3d_frag_f02_empty_path_clears_masks_and_maps() {
        let graph = topology(4, &[(0, 1), (1, 2)]);
        let mut selection = FragmentSubsetInfo::default();
        selection.selected_atoms.resize(4);
        selection.selected_atoms.set(3);
        selection.selected_bonds.resize(2);
        selection.selected_bonds.set(1);
        selection
            .atom_mapping
            .insert(AtomId::new(3), AtomId::new(0));
        selection
            .bond_mapping
            .insert(BondId::new(1), BondId::new(0));

        get_subset_info_for_atom_path(&graph, &[], &mut selection);

        assert_selection(&selection, 4, 2, &[], &[]);
    }

    #[test]
    fn cf3d_frag_f02_full_path_selects_every_valid_bond() {
        let graph = topology(5, &[(0, 1), (1, 2), (2, 3), (3, 0), (3, 4)]);
        let path = (0..5).map(AtomId::new).collect::<Vec<_>>();
        let mut selection = FragmentSubsetInfo::default();

        get_subset_info_for_atom_path(&graph, &path, &mut selection);

        assert_selection(&selection, 5, 5, &[0, 1, 2, 3, 4], &[0, 1, 2, 3, 4]);
    }

    #[test]
    fn cf3d_frag_f02_atom_path_permutations_have_same_selection() {
        let graph = topology(4, &[(0, 1), (1, 2), (2, 3)]);
        let mut first = FragmentSubsetInfo::default();
        let mut second = FragmentSubsetInfo::default();

        get_subset_info_for_atom_path(
            &graph,
            &[AtomId::new(2), AtomId::new(0), AtomId::new(1)],
            &mut first,
        );
        get_subset_info_for_atom_path(
            &graph,
            &[AtomId::new(1), AtomId::new(2), AtomId::new(0)],
            &mut second,
        );

        assert_selection(&first, 4, 3, &[0, 1, 2], &[0, 1]);
        assert_eq!(first, second);
    }

    #[test]
    fn cf3d_frag_f02_duplicate_and_out_of_range_rows_are_ignored_or_collapsed() {
        let graph = topology(3, &[(0, 1), (1, 2)]);
        let path = [
            AtomId::new(1),
            AtomId::new(1),
            AtomId::new(9),
            AtomId::new(0),
            AtomId::new(usize::MAX),
        ];
        let mut selection = FragmentSubsetInfo::default();

        get_subset_info_for_atom_path(&graph, &path, &mut selection);

        assert_selection(&selection, 3, 2, &[0, 1], &[0]);
    }

    #[test]
    fn cf3d_frag_f02_isolated_selected_atom_is_retained() {
        let graph = topology(4, &[(0, 1), (1, 2)]);
        let mut selection = FragmentSubsetInfo::default();

        get_subset_info_for_atom_path(&graph, &[AtomId::new(0), AtomId::new(3)], &mut selection);

        assert_selection(&selection, 4, 2, &[0, 3], &[]);
    }
}

#[cfg(test)]
mod cf3d_frag_f03_tests {
    use std::collections::BTreeMap;

    use super::{FragmentSubsetInfo, copy_selected_atoms};
    use cosmolkit_model::{
        Atom, AtomId, AtomPdbResidueInfo, AtomSpec, TemplateAttachment, TemplateAttachmentOrder,
        TopologyBlock,
    };
    use cosmolkit_types::{ChiralTag, Element, Hybridization};

    fn basic_atom(id: usize, element: Element) -> Atom {
        Atom::from_spec(AtomId::new(id), AtomSpec::new(element))
    }

    fn rich_atom(id: usize, attachment_target: usize) -> Atom {
        let residue = AtomPdbResidueInfo::new("N7", 73, "LIG", 12, "Q", true)
            .with_alt_loc("B")
            .with_insertion_code("C")
            .with_occupancy(0.75)
            .with_temp_factor(18.25)
            .with_secondary_structure(4)
            .with_segment_number(9)
            .with_monomer_class("fragment-fixture");
        let attachment_order = TemplateAttachmentOrder::new(vec![TemplateAttachment::new(
            AtomId::new(attachment_target),
            "ordered-port",
        )])
        .expect("one fixed attachment is valid");
        let spec = AtomSpec::new(Element::N)
            .with_formal_charge(-2)
            .with_explicit_hydrogens(3)
            .with_chiral_tag(ChiralTag::TetrahedralCw)
            .with_chiral_permutation(5)
            .with_unknown_stereo(true)
            .with_mol_parity(6)
            .with_mol_inversion_flag(7)
            .with_implicit_hydrogen(true)
            .with_tracked_isotopic_hydrogens(vec![2, 3])
            .with_aromatic(true)
            .with_isotope(15)
            .with_atom_map(42)
            .with_no_implicit(true)
            .with_radical_electrons(1)
            .with_hybridization(Hybridization::Sp2)
            .with_pdb_residue_info(residue)
            .with_template_attachment_order(attachment_order)
            .with_prop("ordinary", format!("source-{id}"))
            .expect("fixed ordinary property is valid")
            .with_computed_prop("_cache", format!("computed-{id}"))
            .expect("fixed computed property is valid");
        let mut atom = Atom::from_spec(AtomId::new(id), spec);
        atom.set_temporary_flags(if id == 1 {
            0x8123_4567_89AB_CDEF
        } else {
            u64::MAX
        });
        atom
    }

    fn topology(atoms: Vec<Atom>) -> TopologyBlock {
        TopologyBlock::try_from_parts(atoms, Vec::new(), Vec::new(), Vec::new())
            .expect("fixed atom-only topology is valid")
    }

    fn selection(atom_count: usize, rows: &[usize]) -> FragmentSubsetInfo {
        let mut selection = FragmentSubsetInfo::default();
        selection.selected_atoms.resize(atom_count);
        for row in rows {
            selection.selected_atoms.set(*row);
        }
        selection
    }

    #[test]
    fn cf3d_frag_f03_gapped_copy_preserves_all_atom_state_and_clears_only_computed_props() {
        let graph = topology(vec![
            basic_atom(0, Element::C),
            rich_atom(1, 3),
            basic_atom(2, Element::S),
            rich_atom(3, 1),
        ]);
        let original = graph.clone();
        let mut selection = selection(4, &[1, 3]);

        let copied = copy_selected_atoms(&graph, &mut selection).unwrap();

        assert_eq!(copied.len(), 2);
        assert_eq!(copied[0].id(), AtomId::new(0));
        assert_eq!(copied[1].id(), AtomId::new(1));
        assert_eq!(
            selection.atom_mapping,
            BTreeMap::from([
                (AtomId::new(1), AtomId::new(0)),
                (AtomId::new(3), AtomId::new(1)),
            ])
        );
        assert!(selection.bond_mapping.is_empty());

        let source = &graph.atoms[1];
        let copied_atom = &copied[0];
        assert_eq!(copied_atom.element(), Element::N);
        assert_eq!(copied_atom.formal_charge(), -2);
        assert_eq!(copied_atom.explicit_hydrogens(), 3);
        assert_eq!(copied_atom.chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(copied_atom.chiral_permutation(), Some(5));
        assert!(copied_atom.unknown_stereo());
        assert_eq!(copied_atom.mol_parity(), Some(6));
        assert_eq!(copied_atom.mol_inversion_flag(), Some(7));
        assert!(copied_atom.implicit_hydrogen());
        assert_eq!(copied_atom.tracked_isotopic_hydrogens(), &[2, 3]);
        assert!(copied_atom.is_aromatic());
        assert_eq!(copied_atom.isotope(), Some(15));
        assert_eq!(copied_atom.atom_map(), Some(42));
        assert!(copied_atom.no_implicit());
        assert_eq!(copied_atom.radical_electrons(), 1);
        assert_eq!(copied_atom.hybridization(), Hybridization::Sp2);
        assert_eq!(copied_atom.temporary_flags(), 0x8123_4567_89AB_CDEF);
        assert_eq!(copied_atom.pdb_residue_info(), source.pdb_residue_info());
        assert_eq!(
            copied_atom.template_attachment_order(),
            source.template_attachment_order()
        );
        assert_eq!(
            copied_atom.prop("ordinary"),
            Some(&cosmolkit_model::PropertyValue::String(
                ("source-1".to_owned()).into()
            ))
        );
        assert_eq!(copied_atom.prop("_cache"), None);
        assert!(!copied_atom.is_prop_computed("_cache").unwrap());

        assert_eq!(copied[1].temporary_flags(), u64::MAX);
        assert_eq!(
            copied[1].prop("ordinary"),
            Some(&cosmolkit_model::PropertyValue::String(
                ("source-3".to_owned()).into()
            ))
        );
        assert_eq!(copied[1].prop("_cache"), None);
        assert_eq!(graph, original, "copying must not mutate source atoms");
        assert_eq!(graph.atoms[1].temporary_flags(), 0x8123_4567_89AB_CDEF);
        assert_eq!(
            graph.atoms[1].prop("_cache"),
            Some(&cosmolkit_model::PropertyValue::String(
                ("computed-1".to_owned()).into()
            ))
        );
    }

    #[test]
    fn cf3d_frag_f03_empty_selection_copies_no_atoms_or_mapping() {
        let graph = topology(vec![rich_atom(0, 1), rich_atom(1, 0)]);
        let original = graph.clone();
        let mut selection = selection(2, &[]);

        let copied = copy_selected_atoms(&graph, &mut selection).unwrap();

        assert!(copied.is_empty());
        assert!(selection.atom_mapping.is_empty());
        assert_eq!(graph, original);
    }
}

#[cfg(test)]
mod cf3d_frag_f04_tests {
    use std::collections::BTreeMap;

    use super::{FragmentSubsetInfo, SelectedBondCopyError, SelectionMask, copy_selected_bonds};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, BondStereo, BondValueError,
        TopologyBlock,
    };
    use cosmolkit_types::Element;

    fn bond(
        id: usize,
        begin: usize,
        end: usize,
        stereo: BondStereo,
        stereo_atoms: Option<[usize; 2]>,
    ) -> Bond {
        let mut spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Double)
            .with_stereo(stereo);
        if let Some([first, second]) = stereo_atoms {
            spec = spec.with_stereo_atoms(AtomId::new(first), AtomId::new(second));
        }
        Bond::from_spec(BondId::new(id), spec)
    }

    fn rich_bond(
        id: usize,
        begin: usize,
        end: usize,
        stereo: BondStereo,
        stereo_atoms: Option<[usize; 2]>,
    ) -> Bond {
        let mut spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Double)
            .with_stereo(stereo);
        if let Some([first, second]) = stereo_atoms {
            spec = spec.with_stereo_atoms(AtomId::new(first), AtomId::new(second));
        }
        let spec = spec
            .with_prop("ordinary", format!("bond-{id}"))
            .expect("fixed ordinary bond property is valid")
            .with_computed_prop("_computed", format!("cache-{id}"))
            .expect("fixed computed bond property is valid");
        let mut result = Bond::from_spec(BondId::new(id), spec);
        result.set_temporary_flags(1_u64 << id);
        result
    }

    fn topology(atom_count: usize, bonds: Vec<Bond>) -> TopologyBlock {
        let atoms = (0..atom_count)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed source-shaped bond topology is valid")
    }

    fn selection(
        bond_count: usize,
        bond_rows: &[usize],
        atom_mapping: &[(usize, usize)],
    ) -> FragmentSubsetInfo {
        let mut selection = FragmentSubsetInfo::default();
        selection.selected_bonds = SelectionMask::default();
        selection.selected_bonds.resize(bond_count);
        for row in bond_rows {
            selection.selected_bonds.set(*row);
        }
        selection.atom_mapping = atom_mapping
            .iter()
            .map(|(old, new)| (AtomId::new(*old), AtomId::new(*new)))
            .collect();
        selection
    }

    fn identity_mapping(atom_count: usize) -> Vec<(usize, usize)> {
        (0..atom_count).map(|index| (index, index)).collect()
    }

    #[test]
    fn cf3d_frag_f04_endpoint_membership_combinations_return_source_error() {
        let source = topology(2, vec![bond(0, 0, 1, BondStereo::None, None)]);

        for membership in 0..4 {
            let mut atom_mapping = Vec::new();
            if membership & 1 != 0 {
                atom_mapping.push((0, 1));
            }
            if membership & 2 != 0 {
                atom_mapping.push((1, 0));
            }
            let mut info = selection(1, &[0], &atom_mapping);
            let copied = copy_selected_bonds(&source, &mut info);

            if membership == 3 {
                let copied = copied.expect("both endpoint mappings are present");
                assert_eq!(copied.len(), 1);
                assert_eq!(copied[0].begin(), AtomId::new(1));
                assert_eq!(copied[0].end(), AtomId::new(0));
                assert_eq!(
                    info.bond_mapping,
                    BTreeMap::from([(BondId::new(0), BondId::new(0))])
                );
            } else {
                let expected_atom = if membership & 1 == 0 { 0 } else { 1 };
                assert_eq!(
                    copied,
                    Err(SelectedBondCopyError::MissingEndpointMapping {
                        bond: BondId::new(0),
                        atom: AtomId::new(expected_atom),
                    })
                );
                assert_eq!(
                    copied.unwrap_err().to_string(),
                    "copyMolSubset: subset bonds contain atoms not contained in subset atoms"
                );
                assert!(info.bond_mapping.is_empty());
            }
        }
    }

    #[test]
    fn cf3d_frag_f04_first_missing_selected_bond_follows_source_row_order() {
        let source = topology(
            4,
            vec![
                bond(0, 0, 1, BondStereo::None, None),
                bond(1, 2, 3, BondStereo::None, None),
            ],
        );
        let mut info = selection(2, &[0, 1], &[(0, 0), (1, 1), (2, 2)]);

        let result = copy_selected_bonds(&source, &mut info);

        assert_eq!(
            result,
            Err(SelectedBondCopyError::MissingEndpointMapping {
                bond: BondId::new(1),
                atom: AtomId::new(3),
            })
        );
        assert_eq!(
            info.bond_mapping,
            BTreeMap::from([(BondId::new(0), BondId::new(0))])
        );
    }

    #[test]
    fn cf3d_frag_f04_every_bond_stereo_tag_and_source_row_is_copied() {
        let tags = [
            BondStereo::None,
            BondStereo::Any,
            BondStereo::Z,
            BondStereo::E,
            BondStereo::Cis,
            BondStereo::Trans,
            BondStereo::AtropCw,
            BondStereo::AtropCcw,
        ];
        let source = topology(
            16,
            tags.iter()
                .enumerate()
                .map(|(row, tag)| rich_bond(row, row * 2, row * 2 + 1, *tag, Some([14, 15])))
                .collect(),
        );
        let all_bonds = (0..tags.len()).collect::<Vec<_>>();
        let mut info = selection(tags.len(), &all_bonds, &identity_mapping(16));

        let copied = copy_selected_bonds(&source, &mut info).expect("all maps are present");

        assert_eq!(copied.len(), tags.len());
        assert_eq!(
            info.bond_mapping,
            (0..tags.len())
                .map(|row| (BondId::new(row), BondId::new(row)))
                .collect()
        );
        for (row, copied_bond) in copied.iter().enumerate() {
            let original = &source.bonds[row];
            assert_eq!(copied_bond.id(), BondId::new(row));
            assert_eq!(copied_bond.begin(), original.begin());
            assert_eq!(copied_bond.end(), original.end());
            assert_eq!(copied_bond.stereo(), tags[row]);
            assert_eq!(
                copied_bond.stereo_atoms(),
                Some([AtomId::new(14), AtomId::new(15)])
            );
            assert_eq!(
                copied_bond.prop("ordinary"),
                Some(&cosmolkit_model::PropertyValue::String(
                    (format!("bond-{row}")).into()
                ))
            );
            assert_eq!(
                copied_bond.prop("_computed"),
                Some(&cosmolkit_model::PropertyValue::String(
                    (format!("cache-{row}")).into()
                ))
            );
            assert!(copied_bond.is_prop_computed("_computed").unwrap());
            assert_eq!(copied_bond.temporary_flags(), 1_u64 << row);
        }
    }

    #[test]
    fn cf3d_frag_f04_missing_stereo_references_clear_without_retagging() {
        for stereo in [BondStereo::Cis, BondStereo::Trans] {
            for retained_refs in 0..4 {
                let source = topology(4, vec![bond(0, 0, 1, stereo, Some([2, 3]))]);
                let mut atom_mapping = vec![(0, 0), (1, 1)];
                if retained_refs & 1 != 0 {
                    atom_mapping.push((2, 2));
                }
                if retained_refs & 2 != 0 {
                    atom_mapping.push((3, 3));
                }
                let mut info = selection(1, &[0], &atom_mapping);

                let copied = copy_selected_bonds(&source, &mut info)
                    .expect("the selected bond endpoints are mapped");
                let output = &copied[0];

                assert_eq!(output.stereo(), stereo);
                if retained_refs == 3 {
                    assert_eq!(
                        output.stereo_atoms(),
                        Some([AtomId::new(2), AtomId::new(3)])
                    );
                    assert_eq!(output.validate(), Ok(()));
                } else {
                    assert_eq!(output.stereo_atoms(), None);
                    assert_eq!(output.validate(), Err(BondValueError::StereoAtomsRequired));
                }
            }
        }

        let source = topology(2, vec![bond(0, 0, 1, BondStereo::Z, None)]);
        let mut info = selection(1, &[0], &[(0, 0), (1, 1)]);
        let copied = copy_selected_bonds(&source, &mut info).expect("endpoints are mapped");
        assert_eq!(copied[0].stereo(), BondStereo::Z);
        assert_eq!(copied[0].stereo_atoms(), None);
        assert_eq!(copied[0].validate(), Ok(()));
    }

    #[test]
    fn cf3d_frag_f04_second_stereo_lookup_remaps_old_new_key_collision() {
        let source = topology(4, vec![bond(0, 1, 3, BondStereo::Cis, Some([1, 3]))]);
        let mut info = selection(1, &[0], &[(1, 0), (3, 1)]);

        let copied = copy_selected_bonds(&source, &mut info).expect("endpoints are mapped");

        assert_eq!(copied[0].begin(), AtomId::new(0));
        assert_eq!(copied[0].end(), AtomId::new(1));
        assert_eq!(copied[0].stereo(), BondStereo::Cis);
        assert_eq!(
            copied[0].stereo_atoms(),
            Some([AtomId::new(0), AtomId::new(0)])
        );
    }

    #[test]
    fn cf3d_frag_f04_gapped_copy_preserves_props_and_flags_through_outer_clear() {
        let source = topology(
            8,
            (0..4)
                .map(|row| {
                    rich_bond(
                        row,
                        row * 2,
                        row * 2 + 1,
                        if row % 2 == 0 {
                            BondStereo::Cis
                        } else {
                            BondStereo::Trans
                        },
                        Some([6, 7]),
                    )
                })
                .collect(),
        );
        let mut info = selection(4, &[1, 3], &identity_mapping(8));

        let mut copied = copy_selected_bonds(&source, &mut info).expect("all maps are present");

        assert_eq!(copied.len(), 2);
        assert_eq!(copied[0].id(), BondId::new(0));
        assert_eq!(copied[1].id(), BondId::new(1));
        assert_eq!(
            info.bond_mapping,
            BTreeMap::from([
                (BondId::new(1), BondId::new(0)),
                (BondId::new(3), BondId::new(1)),
            ])
        );
        for (new_row, source_row) in [(0, 1), (1, 3)] {
            assert_eq!(
                copied[new_row].prop("ordinary"),
                Some(&cosmolkit_model::PropertyValue::String(
                    (format!("bond-{source_row}")).into()
                ))
            );
            assert_eq!(
                copied[new_row].prop("_computed"),
                Some(&cosmolkit_model::PropertyValue::String(
                    (format!("cache-{source_row}")).into()
                ))
            );
            assert!(copied[new_row].is_prop_computed("_computed").unwrap());
            assert_eq!(copied[new_row].temporary_flags(), 1_u64 << source_row);
        }

        for bond in &mut copied {
            bond.clear_computed_props();
        }
        for (new_row, source_row) in [(0, 1), (1, 3)] {
            assert_eq!(
                copied[new_row].prop("ordinary"),
                Some(&cosmolkit_model::PropertyValue::String(
                    (format!("bond-{source_row}")).into()
                ))
            );
            assert_eq!(copied[new_row].prop("_computed"), None);
            assert!(!copied[new_row].is_prop_computed("_computed").unwrap());
            assert_eq!(copied[new_row].temporary_flags(), 1_u64 << source_row);
        }
    }
}

#[cfg(test)]
mod cf3d_frag_f05_tests {
    use std::collections::BTreeMap;

    use super::{
        FragmentSubsetInfo, SelectedSubstanceGroupCopyError, SubstanceGroupAtomList,
        copy_selected_substance_groups, is_selected_substance_group,
    };
    use cosmolkit_model::{
        AtomId, BondId, SGroupAttachPoint, SGroupBondRole, SGroupBracket, SGroupBracketStyle,
        SGroupCState, SGroupConnection, SGroupData, SGroupDisplay, SubstanceGroup,
        SubstanceGroupId, SubstanceGroupKind,
    };

    fn atom(index: usize) -> AtomId {
        AtomId::new(index)
    }

    fn bond(index: usize) -> BondId {
        BondId::new(index)
    }

    fn group(
        id: usize,
        atoms: Vec<AtomId>,
        bonds: Vec<BondId>,
        parent_atoms: Vec<AtomId>,
    ) -> SubstanceGroup {
        SubstanceGroup::new(
            SubstanceGroupId::new(id),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms(atoms)
        .with_bonds(bonds)
        .with_parent_atoms(parent_atoms)
    }

    fn selection(
        atom_count: usize,
        bond_count: usize,
        selected_atoms: &[usize],
        selected_bonds: &[usize],
        atom_mapping: &[(usize, usize)],
        bond_mapping: &[(usize, usize)],
    ) -> FragmentSubsetInfo {
        let mut info = FragmentSubsetInfo::default();
        info.selected_atoms.resize(atom_count);
        for &index in selected_atoms {
            info.selected_atoms.set(index);
        }
        info.selected_bonds.resize(bond_count);
        for &index in selected_bonds {
            info.selected_bonds.set(index);
        }
        info.atom_mapping = atom_mapping
            .iter()
            .map(|&(source, mapped)| (atom(source), atom(mapped)))
            .collect();
        info.bond_mapping = bond_mapping
            .iter()
            .map(|&(source, mapped)| (bond(source), bond(mapped)))
            .collect();
        info
    }

    fn atom_case(case: usize, offset: usize) -> (Vec<AtomId>, Vec<usize>) {
        match case {
            0 => (Vec::new(), Vec::new()),
            1 => (
                vec![atom(offset), atom(offset + 1)],
                vec![offset, offset + 1],
            ),
            2 => (vec![atom(offset), atom(offset + 1)], vec![offset + 1]),
            3 => (vec![atom(offset), atom(offset + 1)], vec![offset]),
            _ => unreachable!("fixed source-selection matrix state"),
        }
    }

    fn bond_case(case: usize) -> (Vec<BondId>, Vec<usize>) {
        match case {
            0 => (Vec::new(), Vec::new()),
            1 => (vec![bond(0), bond(1)], vec![0, 1]),
            2 => (vec![bond(0), bond(1)], vec![1]),
            3 => (vec![bond(0), bond(1)], vec![0]),
            _ => unreachable!("fixed source-selection matrix state"),
        }
    }

    #[test]
    fn cf3d_frag_f05_selection_empty_all_each_omission_and_category_and() {
        for atom_state in 0..4 {
            for bond_state in 0..4 {
                for parent_state in 0..4 {
                    let (atoms, selected_atoms) = atom_case(atom_state, 0);
                    let (parent_atoms, selected_parent_atoms) = atom_case(parent_state, 2);
                    let (bonds, selected_bonds) = bond_case(bond_state);
                    let mut all_selected_atoms = selected_atoms;
                    all_selected_atoms.extend(selected_parent_atoms);
                    let info = selection(4, 2, &all_selected_atoms, &selected_bonds, &[], &[]);
                    let source_group = group(0, atoms, bonds, parent_atoms);
                    let expected = atom_state < 2 && bond_state < 2 && parent_state < 2;
                    assert_eq!(
                        is_selected_substance_group(&source_group, &info),
                        expected,
                        "atom/bond/parent states {atom_state}/{bond_state}/{parent_state}"
                    );
                }
            }
        }
    }

    #[test]
    fn cf3d_frag_f05_copy_maps_member_lists_and_preserves_other_fields_in_order() {
        let display = SGroupDisplay {
            brackets: vec![SGroupBracket::new([
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
            ])],
            field_position: Some([2.5, 3.5]),
            display_tag: Some("source-display".into()),
        };
        let attach_point = SGroupAttachPoint {
            atom: atom(5),
            leaving_atom: Some(atom(0)),
            label: Some("R1".into()),
            order: Some(1),
        };
        let cstate = SGroupCState::new(bond(4), [1.25, -2.5, 3.75]);
        let data = SGroupData {
            field_name: Some("FIELD".into()),
            values: vec!["source-value".into()],
            ..SGroupData::default()
        };
        let mut first = group(0, vec![atom(3), atom(1)], vec![], vec![atom(4)])
            .with_head_crossing_bonds(vec![bond(3), bond(4)])
            .with_crossing_bond_correspondence(vec![bond(4), bond(3)])
            .with_label("polymer")
            .with_connection(SGroupConnection::HeadToTail)
            .with_subtype("source-subtype")
            .with_bracket_style(SGroupBracketStyle::Bracket)
            .with_display(display.clone())
            .with_data(data.clone())
            .with_attach_points(vec![attach_point.clone()])
            .with_cstates(vec![cstate])
            .with_prop("source-key", "source-value")
            .unwrap()
            .with_data_field("field-one")
            .with_data_field("field-two");
        first.push_bond_with_role(bond(2), SGroupBondRole::Contained);
        first.push_bond_with_role(bond(0), SGroupBondRole::Crossing);

        let omitted = group(1, vec![atom(2)], Vec::new(), Vec::new()).with_label("partial");
        let all_empty = group(2, Vec::new(), Vec::new(), Vec::new())
            .with_parent(SubstanceGroupId::new(0))
            .with_label("empty-selected");
        let source_groups = vec![first, omitted, all_empty];
        let source_snapshot = source_groups.clone();
        let info = selection(
            6,
            5,
            &[1, 3, 4],
            &[0, 2],
            &[(1, 0), (3, 1), (4, 2)],
            &[(0, 0), (2, 1)],
        );

        let copied = copy_selected_substance_groups(&source_groups, &info)
            .expect("selected source groups have complete in-range maps");

        assert_eq!(source_groups, source_snapshot);
        assert_eq!(copied.len(), 2);
        assert_eq!(copied[0].id(), SubstanceGroupId::new(0));
        assert_eq!(copied[1].id(), SubstanceGroupId::new(1));
        assert_eq!(copied[0].atoms(), &[atom(1), atom(0)]);
        assert_eq!(copied[0].parent_atoms(), &[atom(2)]);
        assert_eq!(copied[0].bonds(), &[bond(1), bond(0)]);
        assert_eq!(copied[0].bond_role(bond(1)), SGroupBondRole::Contained);
        assert_eq!(copied[0].bond_role(bond(0)), SGroupBondRole::Crossing);
        assert_eq!(copied[0].head_crossing_bonds(), &[bond(3), bond(4)]);
        assert_eq!(
            copied[0].crossing_bond_correspondence(),
            &[bond(4), bond(3)]
        );
        assert_eq!(copied[0].attach_points(), &[attach_point]);
        assert_eq!(copied[0].cstates(), &[cstate]);
        assert_eq!(
            copied[0].bracket_style(),
            Some(&SGroupBracketStyle::Bracket)
        );
        assert_eq!(copied[0].display(), Some(&display));
        assert_eq!(copied[0].data(), Some(&data));
        assert_eq!(
            copied[0].label().map(super::fixed_property_text),
            Some("polymer")
        );
        assert_eq!(copied[0].connection(), Some(&SGroupConnection::HeadToTail));
        assert_eq!(
            copied[0].subtype().map(super::fixed_property_text),
            Some("source-subtype")
        );
        assert_eq!(
            copied[0].props(),
            &BTreeMap::from([(
                "source-key".into(),
                cosmolkit_model::PropertyValue::String("source-value".into())
            )])
        );
        assert_eq!(
            copied[0].data_fields(),
            &["field-one".into(), "field-two".into()]
        );

        assert_eq!(copied[1].atoms(), &[]);
        assert_eq!(copied[1].bonds(), &[]);
        assert_eq!(copied[1].parent_atoms(), &[]);
        assert_eq!(copied[1].parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(
            copied[1].label().map(super::fixed_property_text),
            Some("empty-selected")
        );
    }

    #[test]
    fn cf3d_frag_f05_mapping_errors_keep_category_member_and_source_order() {
        let atom_group = group(7, vec![atom(0)], Vec::new(), Vec::new());
        let selected_atom_without_map = selection(1, 0, &[0], &[], &[], &[]);
        assert_eq!(
            copy_selected_substance_groups(&[atom_group], &selected_atom_without_map),
            Err(SelectedSubstanceGroupCopyError::MissingAtomMapping {
                group: SubstanceGroupId::new(7),
                list: SubstanceGroupAtomList::Atoms,
                atom: atom(0),
            })
        );

        let parent_group = group(8, Vec::new(), Vec::new(), vec![atom(0)]);
        let selected_parent_without_map = selection(1, 0, &[0], &[], &[], &[]);
        assert_eq!(
            copy_selected_substance_groups(&[parent_group], &selected_parent_without_map),
            Err(SelectedSubstanceGroupCopyError::MissingAtomMapping {
                group: SubstanceGroupId::new(8),
                list: SubstanceGroupAtomList::ParentAtoms,
                atom: atom(0),
            })
        );

        let bond_group = group(9, Vec::new(), vec![bond(0)], Vec::new());
        let selected_bond_without_map = selection(0, 1, &[], &[0], &[], &[]);
        assert_eq!(
            copy_selected_substance_groups(&[bond_group], &selected_bond_without_map),
            Err(SelectedSubstanceGroupCopyError::MissingBondMapping {
                group: SubstanceGroupId::new(9),
                bond: bond(0),
            })
        );

        let ordered_group = group(10, vec![atom(0), atom(1)], Vec::new(), Vec::new());
        let missing_later_member = selection(2, 0, &[0, 1], &[], &[(0, 0)], &[]);
        assert_eq!(
            copy_selected_substance_groups(&[ordered_group], &missing_later_member),
            Err(SelectedSubstanceGroupCopyError::MissingAtomMapping {
                group: SubstanceGroupId::new(10),
                list: SubstanceGroupAtomList::Atoms,
                atom: atom(1),
            })
        );

        let range_then_missing = group(11, vec![atom(0), atom(1)], Vec::new(), Vec::new());
        let map_error_precedes_setter_range_error = selection(2, 0, &[0, 1], &[], &[(0, 2)], &[]);
        assert_eq!(
            copy_selected_substance_groups(
                &[range_then_missing],
                &map_error_precedes_setter_range_error
            ),
            Err(SelectedSubstanceGroupCopyError::MissingAtomMapping {
                group: SubstanceGroupId::new(11),
                list: SubstanceGroupAtomList::Atoms,
                atom: atom(1),
            })
        );
    }

    #[test]
    fn cf3d_frag_f05_destination_ranges_follow_three_source_setters() {
        let atom_group = group(12, vec![atom(0)], Vec::new(), Vec::new());
        let atom_range = selection(1, 0, &[0], &[], &[(0, 1)], &[]);
        assert_eq!(
            copy_selected_substance_groups(&[atom_group], &atom_range),
            Err(SelectedSubstanceGroupCopyError::AtomOutOfRange {
                group: SubstanceGroupId::new(12),
                list: SubstanceGroupAtomList::Atoms,
                source_atom: atom(0),
                mapped_atom: atom(1),
                atom_count: 1,
            })
        );

        let parent_group = group(13, vec![atom(0)], Vec::new(), vec![atom(1)]);
        let parent_range = selection(2, 0, &[0, 1], &[], &[(0, 0), (1, 2)], &[]);
        assert_eq!(
            copy_selected_substance_groups(&[parent_group], &parent_range),
            Err(SelectedSubstanceGroupCopyError::AtomOutOfRange {
                group: SubstanceGroupId::new(13),
                list: SubstanceGroupAtomList::ParentAtoms,
                source_atom: atom(1),
                mapped_atom: atom(2),
                atom_count: 2,
            })
        );

        let bond_group = group(14, Vec::new(), vec![bond(0), bond(1)], Vec::new());
        let bond_range = selection(0, 2, &[], &[0, 1], &[], &[(0, 0), (1, 2)]);
        assert_eq!(
            copy_selected_substance_groups(&[bond_group], &bond_range),
            Err(SelectedSubstanceGroupCopyError::BondOutOfRange {
                group: SubstanceGroupId::new(14),
                source_bond: bond(1),
                mapped_bond: bond(2),
                bond_count: 2,
            })
        );
    }
}

#[cfg(test)]
mod cf3d_frag_f10_tests {
    use super::{AtomPathSubsetCopyError, copy_mol_subset_atom_path};
    use crate::{
        PropertyCacheError, SanitizeError, SanitizeParams, SanitizeStage, ValenceError,
        ValencePhase, sanitize_topology,
    };
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer2D, Conformer3D, CoordinateBlock,
        CoordinateDimension, MoleculeProperties, StereoGroup, StereoGroupKind, SubstanceGroup,
        SubstanceGroupId, SubstanceGroupKind, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    fn benzene() -> (TopologyBlock, CoordinateBlock) {
        let atoms = (0..6)
            .map(|index| {
                AtomSpec::new(Element::C)
                    .with_prop("source-row", format!("atom-{index}"))
                    .unwrap()
                    .with_computed_prop("subset-cache", format!("atom-cache-{index}"))
                    .unwrap()
            })
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect();
        let bonds = (0..6)
            .map(|index| {
                BondSpec::new(
                    AtomId::new(index),
                    AtomId::new((index + 1) % 6),
                    if index % 2 == 0 {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    },
                )
                .with_prop("source-row", format!("bond-{index}"))
                .unwrap()
                .with_computed_prop("subset-cache", format!("bond-cache-{index}"))
                .unwrap()
            })
            .enumerate()
            .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
            .collect();
        let substance_group = SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms(vec![AtomId::new(0), AtomId::new(2)])
        .with_bonds(vec![BondId::new(0), BondId::new(2)])
        .with_parent_atoms(vec![AtomId::new(1)])
        .with_label("source-sgroup");
        let stereo_group = StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(5)],
            vec![BondId::new(0)],
        )
        .with_id(17)
        .with_write_id(9);
        let topology =
            TopologyBlock::try_from_parts(atoms, bonds, vec![substance_group], vec![stereo_group])
                .unwrap();
        let coordinates = CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(
                    31,
                    vec![
                        [0.0, 10.0],
                        [1.0, 11.0],
                        [2.0, 12.0],
                        [3.0, 13.0],
                        [4.0, 14.0],
                        [5.0, 15.0],
                    ],
                )
                .with_prop("kind", "source-2d"),
            ],
            conformers_3d: vec![
                Conformer3D::new(
                    45,
                    vec![
                        [0.0, 10.0, 20.0],
                        [1.0, 11.0, 21.0],
                        [2.0, 12.0, 22.0],
                        [3.0, 13.0, 23.0],
                        [4.0, 14.0, 24.0],
                        [5.0, 15.0, 25.0],
                    ],
                    false,
                )
                .with_prop("kind", "source-3d"),
            ],
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            source_conformer_order: None,
        };

        (topology, coordinates)
    }

    fn expected_full_topology(source: &TopologyBlock, sanitize: bool) -> TopologyBlock {
        let mut atoms = source.atoms.clone();
        for atom in &mut atoms {
            atom.clear_computed_props();
        }
        let mut topology = TopologyBlock::try_from_parts(
            atoms,
            source.bonds.clone(),
            source.substance_groups.clone(),
            source.stereo_groups.clone(),
        )
        .unwrap();
        if sanitize {
            topology = sanitize_topology(&topology, &SanitizeParams::default())
                .unwrap()
                .topology;
        }
        for atom in &mut topology.atoms {
            atom.clear_computed_props();
        }
        for bond in &mut topology.bonds {
            bond.clear_computed_props();
        }
        topology
    }

    fn expected_full_coordinates(copy_coordinates: bool) -> CoordinateBlock {
        if !copy_coordinates {
            return CoordinateBlock::default();
        }
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                31,
                vec![
                    [0.0, 10.0],
                    [1.0, 11.0],
                    [2.0, 12.0],
                    [3.0, 13.0],
                    [4.0, 14.0],
                    [5.0, 15.0],
                ],
            )],
            conformers_3d: vec![Conformer3D::new(
                45,
                vec![
                    [0.0, 10.0, 20.0],
                    [1.0, 11.0, 21.0],
                    [2.0, 12.0, 22.0],
                    [3.0, 13.0, 23.0],
                    [4.0, 14.0, 24.0],
                    [5.0, 15.0, 25.0],
                ],
                false,
            )],
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            source_conformer_order: None,
        }
    }

    #[test]
    fn cf3d_frag_f10_full_path_option_product_preserves_source_pipeline() {
        let (source_topology, source_coordinates) = benzene();
        let original_topology = source_topology.clone();
        let original_coordinates = source_coordinates.clone();
        let full_path = [5, 4, 3, 2, 1, 0].map(AtomId::new);

        for sanitize in [false, true] {
            for copy_coordinates in [false, true] {
                let copied = copy_mol_subset_atom_path(
                    &source_topology,
                    &source_coordinates,
                    &full_path,
                    sanitize,
                    copy_coordinates,
                )
                .unwrap();
                assert_eq!(
                    copied.topology,
                    expected_full_topology(&source_topology, sanitize),
                    "sanitize={sanitize}"
                );
                assert_eq!(
                    copied.coordinates,
                    expected_full_coordinates(copy_coordinates),
                    "copy_coordinates={copy_coordinates}"
                );
                assert_eq!(copied.molecule_properties, MoleculeProperties::default());
                assert_eq!(copied.topology.atoms.len(), 6);
                assert_eq!(copied.topology.bonds.len(), 6);
                assert_eq!(
                    copied.topology.substance_groups,
                    source_topology.substance_groups
                );
                assert_eq!(copied.topology.stereo_groups, source_topology.stereo_groups);
                assert_eq!(
                    copied
                        .topology
                        .atoms
                        .iter()
                        .map(|atom| atom.prop("source-row").unwrap().clone())
                        .collect::<Vec<_>>(),
                    ["atom-0", "atom-1", "atom-2", "atom-3", "atom-4", "atom-5",].map(|value| {
                        cosmolkit_model::PropertyValue::String((value.to_owned()).into())
                    })
                );
                assert!(
                    copied
                        .topology
                        .atoms
                        .iter()
                        .all(|atom| atom.prop("subset-cache").is_none())
                );
                assert!(
                    copied
                        .topology
                        .bonds
                        .iter()
                        .all(|bond| bond.prop("subset-cache").is_none())
                );
                if sanitize {
                    assert!(copied.topology.atoms.iter().all(Atom::is_aromatic));
                    assert!(copied.topology.bonds.iter().all(Bond::is_aromatic));
                } else {
                    assert!(copied.topology.atoms.iter().all(|atom| !atom.is_aromatic()));
                    assert!(copied.topology.bonds.iter().all(|bond| !bond.is_aromatic()));
                }
            }
        }

        assert_eq!(source_topology, original_topology);
        assert_eq!(source_coordinates, original_coordinates);
    }

    #[test]
    fn cf3d_frag_f10_empty_and_gapped_paths_match_source_rows() {
        let (source_topology, source_coordinates) = benzene();
        let original_topology = source_topology.clone();
        let original_coordinates = source_coordinates.clone();

        let empty =
            copy_mol_subset_atom_path(&source_topology, &source_coordinates, &[], false, true)
                .unwrap();
        assert_eq!(empty.topology, TopologyBlock::default());
        assert_eq!(empty.molecule_properties, MoleculeProperties::default());
        assert_eq!(
            empty.coordinates,
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(31, Vec::new())],
                conformers_3d: vec![Conformer3D::new(45, Vec::new(), false)],
                source_coordinate_dim: Some(CoordinateDimension::ThreeD),
                source_conformer_order: None,
            }
        );

        let gapped = copy_mol_subset_atom_path(
            &source_topology,
            &source_coordinates,
            &[4, 0, 2, 99].map(AtomId::new),
            true,
            true,
        )
        .unwrap();
        assert_eq!(gapped.topology.atoms.len(), 3);
        assert!(gapped.topology.bonds.is_empty());
        assert!(gapped.topology.substance_groups.is_empty());
        assert!(gapped.topology.stereo_groups.is_empty());
        assert_eq!(
            gapped
                .topology
                .atoms
                .iter()
                .map(|atom| atom.prop("source-row").unwrap().clone())
                .collect::<Vec<_>>(),
            ["atom-0", "atom-2", "atom-4"]
                .map(|value| cosmolkit_model::PropertyValue::String((value.to_owned()).into()))
        );
        assert!(gapped.topology.atoms.iter().all(|atom| !atom.is_aromatic()));
        assert_eq!(
            gapped.coordinates,
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(
                    31,
                    vec![[0.0, 10.0], [2.0, 12.0], [4.0, 14.0]],
                )],
                conformers_3d: vec![Conformer3D::new(
                    45,
                    vec![[0.0, 10.0, 20.0], [2.0, 12.0, 22.0], [4.0, 14.0, 24.0]],
                    false,
                )],
                source_coordinate_dim: Some(CoordinateDimension::ThreeD),
                source_conformer_order: None,
            }
        );
        assert_eq!(source_topology, original_topology);
        assert_eq!(source_coordinates, original_coordinates);
    }

    #[test]
    fn cf3d_frag_f10_returns_first_typed_sanitize_failure() {
        let mut atom_specs = vec![AtomSpec::new(Element::C)];
        atom_specs.extend((0..5).map(|_| AtomSpec::new(Element::H)));
        let atoms = atom_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect();
        let bonds = (1..=5)
            .map(|hydrogen| {
                Bond::from_spec(
                    BondId::new(hydrogen - 1),
                    BondSpec::new(AtomId::new(0), AtomId::new(hydrogen), BondOrder::Single),
                )
            })
            .collect();
        let source_topology =
            TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
        let original_topology = source_topology.clone();

        let error = copy_mol_subset_atom_path(
            &source_topology,
            &CoordinateBlock::default(),
            &[0, 1, 2, 3, 4, 5].map(AtomId::new),
            true,
            false,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            AtomPathSubsetCopyError::Sanitize(SanitizeError::Properties {
                stage: SanitizeStage::Properties,
                source: PropertyCacheError::Valence(ValenceError::InvalidValence {
                    atom,
                    atomic_number: 6,
                    phase: ValencePhase::Explicit,
                    ..
                }),
            }) if atom == AtomId::new(0)
        ));
        assert_eq!(source_topology, original_topology);
    }
}

#[cfg(test)]
mod cf3d_frag_f11_tests {
    use super::{SelectionMask, fragment_has_challenging_features};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, ChiralTag, Element, StereoGroup,
        StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock,
    };
    use cosmolkit_types::BondStereo;

    fn topology(
        chiral_tags: &[ChiralTag],
        bonds: &[(usize, usize, BondStereo)],
        substance_groups: Vec<SubstanceGroup>,
        stereo_groups: Vec<StereoGroup>,
    ) -> TopologyBlock {
        let atoms = chiral_tags
            .iter()
            .copied()
            .enumerate()
            .map(|(index, tag)| {
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C).with_chiral_tag(tag),
                )
            })
            .collect();
        let bonds = bonds
            .iter()
            .copied()
            .enumerate()
            .map(|(index, (begin, end, stereo))| {
                let begin = AtomId::new(begin);
                let end = AtomId::new(end);
                let mut spec = BondSpec::new(begin, end, BondOrder::Double).with_stereo(stereo);
                if matches!(stereo, BondStereo::Cis | BondStereo::Trans) {
                    spec = spec.with_stereo_atoms(begin, end);
                }
                Bond::from_spec(BondId::new(index), spec)
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, substance_groups, stereo_groups)
            .expect("fixed F11 topology must satisfy model invariants")
    }

    fn challenging(
        topology: &TopologyBlock,
        component: &[usize],
        selected_atoms: &[usize],
    ) -> bool {
        let component = component
            .iter()
            .copied()
            .map(AtomId::new)
            .collect::<Vec<_>>();
        let mut atoms_in_fragment = SelectionMask::default();
        atoms_in_fragment.resize(topology.atoms.len());
        for atom in selected_atoms {
            atoms_in_fragment.set(*atom);
        }
        fragment_has_challenging_features(topology, &component, &atoms_in_fragment)
    }

    fn sgroup(atoms: Vec<AtomId>, bonds: Vec<BondId>, parent_atoms: Vec<AtomId>) -> SubstanceGroup {
        SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms(atoms)
        .with_bonds(bonds)
        .with_parent_atoms(parent_atoms)
    }

    #[test]
    fn cf3d_frag_f11_chiral_tag_exclusions_cover_every_modeled_variant() {
        let cases = [
            (ChiralTag::Unspecified, false),
            (ChiralTag::TetrahedralCw, true),
            (ChiralTag::TetrahedralCcw, true),
            (ChiralTag::Other, false),
            (ChiralTag::Tetrahedral, true),
            (ChiralTag::Allene, true),
            (ChiralTag::SquarePlanar, true),
            (ChiralTag::TrigonalBipyramidal, true),
            (ChiralTag::Octahedral, true),
        ];

        for (tag, expected) in cases {
            let topology = topology(&[tag], &[], Vec::new(), Vec::new());
            assert_eq!(challenging(&topology, &[0], &[0]), expected, "{tag:?}");
        }
    }

    #[test]
    fn cf3d_frag_f11_bond_stereo_exclusions_cover_every_modeled_variant() {
        let cases = [
            (BondStereo::None, false),
            (BondStereo::Any, false),
            (BondStereo::Z, true),
            (BondStereo::E, true),
            (BondStereo::Cis, true),
            (BondStereo::Trans, true),
            (BondStereo::AtropCw, true),
            (BondStereo::AtropCcw, true),
        ];

        for (stereo, expected) in cases {
            let topology = topology(
                &[ChiralTag::Unspecified, ChiralTag::Unspecified],
                &[(0, 1, stereo)],
                Vec::new(),
                Vec::new(),
            );
            assert_eq!(
                challenging(&topology, &[0, 1], &[0, 1]),
                expected,
                "{stereo:?}"
            );
        }
    }

    #[test]
    fn cf3d_frag_f11_crossing_stereo_bond_is_ignored() {
        let topology = topology(
            &[ChiralTag::Unspecified, ChiralTag::Unspecified],
            &[(0, 1, BondStereo::Z)],
            Vec::new(),
            Vec::new(),
        );

        assert!(!challenging(&topology, &[0], &[0]));
        assert!(challenging(&topology, &[0, 1], &[0, 1]));
    }

    #[test]
    fn cf3d_frag_f11_substance_group_checks_atoms_and_parents_but_not_bonds() {
        let plain = &[ChiralTag::Unspecified, ChiralTag::Unspecified];
        let cases = [
            (sgroup(vec![AtomId::new(0)], vec![], vec![]), true),
            (sgroup(vec![AtomId::new(1)], vec![], vec![]), false),
            (sgroup(vec![], vec![], vec![AtomId::new(0)]), true),
            (sgroup(vec![], vec![], vec![AtomId::new(1)]), false),
            (sgroup(vec![], vec![], vec![]), false),
        ];

        for (group, expected) in cases {
            let topology = topology(plain, &[], vec![group], Vec::new());
            assert_eq!(challenging(&topology, &[0], &[0]), expected);
        }

        let topology = topology(
            plain,
            &[(0, 1, BondStereo::None)],
            vec![sgroup(vec![], vec![BondId::new(0)], vec![])],
            Vec::new(),
        );
        assert!(!challenging(&topology, &[0, 1], &[0, 1]));
    }

    #[test]
    fn cf3d_frag_f11_stereo_group_atom_and_bond_membership() {
        let plain = &[ChiralTag::Unspecified, ChiralTag::Unspecified];
        let atom_cases = [
            (vec![AtomId::new(0)], true),
            (vec![AtomId::new(1)], false),
            (vec![], false),
        ];
        for (atoms, expected) in atom_cases {
            let group = StereoGroup::new(StereoGroupKind::Absolute, atoms, vec![]);
            let topology = topology(plain, &[], Vec::new(), vec![group]);
            assert_eq!(challenging(&topology, &[0], &[0]), expected);
        }

        let bond_cases = [(vec![0, 1], true), (vec![0], false), (vec![], false)];
        for (selected_atoms, expected) in bond_cases {
            let group = StereoGroup::new(StereoGroupKind::Or, vec![], vec![BondId::new(0)]);
            let topology = topology(plain, &[(0, 1, BondStereo::None)], Vec::new(), vec![group]);
            let component = if selected_atoms.len() == 2 {
                vec![0, 1]
            } else {
                vec![0]
            };
            assert_eq!(
                challenging(&topology, &component, &selected_atoms),
                expected
            );
        }
    }

    #[test]
    fn cf3d_frag_f11_earlier_atom_stereo_short_circuits_later_memberships() {
        let topology = topology(
            &[ChiralTag::TetrahedralCw, ChiralTag::Unspecified],
            &[(0, 1, BondStereo::None)],
            vec![sgroup(vec![AtomId::new(0)], vec![], vec![])],
            vec![StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(0)],
                vec![BondId::new(0)],
            )],
        );

        assert!(challenging(&topology, &[0, 1], &[0, 1]));
    }
}

#[cfg(test)]
mod cf3d_frag_f15_tests {
    use super::{
        FullCopyComponent, SelectionMask, copy_full_molecule_remove_atoms_outside_component,
    };
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, Conformer3D,
        CoordinateBlock, CoordinateDimension, MoleculeProperties, SdfPropertyList,
        SdfPropertyListTarget, StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupId,
        SubstanceGroupKind, TopologyBlock, TopologyMapping,
    };
    use cosmolkit_types::{BondStereo, Element};

    fn input(component_count: usize) -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
        let atom_count = component_count * 2;
        let atoms = (0..atom_count)
            .map(|index| {
                let spec = AtomSpec::new(Element::C)
                    .with_prop("source-row", format!("atom-{index}"))
                    .unwrap()
                    .with_computed_prop("atom-cache", format!("cache-{index}"))
                    .unwrap();
                Atom::from_spec(AtomId::new(index), spec)
            })
            .collect();
        let bonds = (0..component_count)
            .map(|component| {
                let begin = AtomId::new(component * 2);
                let end = AtomId::new(component * 2 + 1);
                let spec = BondSpec::new(begin, end, BondOrder::Double)
                    .with_stereo(BondStereo::Cis)
                    .with_stereo_atoms(begin, end)
                    .with_prop("source-row", format!("bond-{component}"))
                    .unwrap()
                    .with_computed_prop("bond-cache", format!("cache-{component}"))
                    .unwrap();
                Bond::from_spec(BondId::new(component), spec)
            })
            .collect();
        let substance_groups = (0..component_count)
            .map(|component| {
                SubstanceGroup::new(
                    SubstanceGroupId::new(component),
                    SubstanceGroupKind::StructuralRepeatUnit,
                )
                .with_atoms(vec![
                    AtomId::new(component * 2),
                    AtomId::new(component * 2 + 1),
                ])
                .with_bonds(vec![BondId::new(component)])
                .with_label(format!("group-{component}"))
            })
            .collect();
        let stereo_groups = (0..component_count)
            .map(|component| {
                StereoGroup::new(
                    StereoGroupKind::Or,
                    vec![AtomId::new(component * 2)],
                    vec![BondId::new(component)],
                )
                .with_id(100 + component as u32)
                .with_write_id(200 + component as u32)
            })
            .collect();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, substance_groups, stereo_groups)
            .expect("fixed F15 source topology must be valid");

        let conformers_2d = [41, 5]
            .into_iter()
            .map(|id| {
                Conformer2D::new(
                    id,
                    (0..atom_count)
                        .map(|row| [row as f64 + id as f64, row as f64 + 0.25])
                        .collect(),
                )
                .with_prop("source-conformer", format!("2d-{id}"))
            })
            .collect();
        let conformers_3d = [(73, false), (11, true)]
            .into_iter()
            .map(|(id, is_3d)| {
                Conformer3D::new(
                    id,
                    (0..atom_count)
                        .map(|row| [row as f64, row as f64 + 10.0, row as f64 + 20.0])
                        .collect(),
                    is_3d,
                )
                .with_prop("source-conformer", format!("3d-{id}"))
            })
            .collect();
        let coordinates = CoordinateBlock {
            conformers_2d,
            conformers_3d,
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            source_conformer_order: None,
        };
        let properties = MoleculeProperties::default()
            .with_name("source-molecule")
            .with_prop("ordinary", "keep-me")
            .unwrap()
            .with_computed_prop("molecule-cache", "clear-me")
            .unwrap()
            .with_sdf_data_field("FIELD", "source-value")
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Atom,
                "atom-values",
                (0..atom_count)
                    .map(|row| {
                        Some(cosmolkit_model::PropertyValue::String(
                            (format!("atom-value-{row}")).into(),
                        ))
                    })
                    .collect(),
            ))
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Bond,
                "bond-values",
                (0..component_count)
                    .map(|row| {
                        Some(cosmolkit_model::PropertyValue::String(
                            (format!("bond-value-{row}")).into(),
                        ))
                    })
                    .collect(),
            ));

        (topology, coordinates, properties)
    }

    fn mask(atom_count: usize, selected_component: Option<usize>) -> SelectionMask {
        let mut selected = SelectionMask::default();
        selected.resize(atom_count);
        if let Some(component) = selected_component {
            selected.set(component * 2);
            selected.set(component * 2 + 1);
        } else {
            for atom in 0..atom_count {
                selected.set(atom);
            }
        }
        selected
    }

    fn assert_component_copy(
        source: &TopologyBlock,
        source_coordinates: &CoordinateBlock,
        source_properties: &MoleculeProperties,
        copied: &FullCopyComponent,
        component_count: usize,
        selected_component: usize,
    ) {
        let first_atom = selected_component * 2;
        assert_eq!(copied.topology.atoms.len(), 2);
        assert_eq!(copied.topology.bonds.len(), 1);
        assert_eq!(copied.topology.atoms[0].id(), AtomId::new(0));
        assert_eq!(copied.topology.atoms[1].id(), AtomId::new(1));
        assert_eq!(
            copied.topology.atoms[0].prop("source-row"),
            Some(&cosmolkit_model::PropertyValue::String(
                (format!("atom-{first_atom}")).into()
            ))
        );
        assert_eq!(
            copied.topology.atoms[1].prop("source-row"),
            Some(&cosmolkit_model::PropertyValue::String(
                (format!("atom-{}", first_atom + 1)).into()
            ))
        );
        assert!(copied.topology.atoms.iter().all(|atom| {
            atom.prop("atom-cache").is_none()
                && atom
                    .computed_prop_names()
                    .unwrap()
                    .expect("cleared computed StringVector exists")
                    .is_empty()
        }));
        let bond = &copied.topology.bonds[0];
        assert_eq!(bond.id(), BondId::new(0));
        assert_eq!(bond.begin(), AtomId::new(0));
        assert_eq!(bond.end(), AtomId::new(1));
        assert_eq!(bond.stereo(), BondStereo::Cis);
        assert_eq!(bond.stereo_atoms(), Some([AtomId::new(0), AtomId::new(1)]));
        assert_eq!(
            bond.prop("source-row"),
            Some(&cosmolkit_model::PropertyValue::String(
                (format!("bond-{selected_component}")).into()
            ))
        );
        assert!(bond.prop("bond-cache").is_none());
        assert!(
            bond.computed_prop_names()
                .unwrap()
                .expect("cleared computed StringVector exists")
                .is_empty()
        );

        assert_eq!(copied.topology.substance_groups.len(), 1);
        let group = &copied.topology.substance_groups[0];
        assert_eq!(group.id(), SubstanceGroupId::new(0));
        assert_eq!(group.atoms(), &[AtomId::new(0), AtomId::new(1)]);
        assert_eq!(group.bonds(), &[BondId::new(0)]);
        assert_eq!(
            group.label().map(super::fixed_property_text),
            Some(format!("group-{selected_component}")).as_deref()
        );
        assert_eq!(copied.topology.stereo_groups.len(), 1);
        let stereo = &copied.topology.stereo_groups[0];
        assert_eq!(stereo.kind(), StereoGroupKind::Or);
        assert_eq!(stereo.id(), Some(100 + selected_component as u32));
        assert_eq!(stereo.write_id(), 200 + selected_component as u32);
        assert_eq!(stereo.atoms(), &[AtomId::new(0)]);
        assert_eq!(stereo.bonds(), &[BondId::new(0)]);
        assert!(copied.topology.validate().is_ok());

        let expected_atom_old_to_new: Vec<_> = (0..component_count * 2)
            .map(|row| {
                (row == first_atom || row == first_atom + 1).then(|| AtomId::new(row - first_atom))
            })
            .collect();
        assert_eq!(
            copied.mapping.atoms().old_to_new(),
            expected_atom_old_to_new.as_slice()
        );
        assert_eq!(
            copied.mapping.atoms().new_to_old(),
            &[
                Some(AtomId::new(first_atom)),
                Some(AtomId::new(first_atom + 1))
            ]
        );
        let expected_bond_old_to_new: Vec<_> = (0..component_count)
            .map(|row| (row == selected_component).then_some(BondId::new(0)))
            .collect();
        assert_eq!(
            copied.mapping.bonds().old_to_new(),
            expected_bond_old_to_new.as_slice()
        );
        assert_eq!(
            copied.mapping.bonds().new_to_old(),
            &[Some(BondId::new(selected_component))]
        );

        assert_eq!(
            copied.coordinates.source_coordinate_dim,
            source_coordinates.source_coordinate_dim
        );
        assert_eq!(copied.coordinates.conformers_2d.len(), 2);
        for (source_conf, copied_conf) in source_coordinates
            .conformers_2d
            .iter()
            .zip(&copied.coordinates.conformers_2d)
        {
            assert_eq!(copied_conf.id(), source_conf.id());
            assert_eq!(
                copied_conf.coordinates(),
                &[
                    source_conf.coordinates()[first_atom],
                    source_conf.coordinates()[first_atom + 1]
                ]
            );
            assert_eq!(copied_conf.props(), source_conf.props());
        }
        assert_eq!(copied.coordinates.conformers_3d.len(), 2);
        for (source_conf, copied_conf) in source_coordinates
            .conformers_3d
            .iter()
            .zip(&copied.coordinates.conformers_3d)
        {
            assert_eq!(copied_conf.id(), source_conf.id());
            assert_eq!(copied_conf.is_3d(), source_conf.is_3d());
            assert_eq!(
                copied_conf.coordinates(),
                &[
                    source_conf.coordinates()[first_atom],
                    source_conf.coordinates()[first_atom + 1]
                ]
            );
            assert_eq!(copied_conf.props(), source_conf.props());
        }

        assert_eq!(
            copied
                .molecule_properties
                .name()
                .map(super::fixed_property_text),
            Some("source-molecule")
        );
        assert_eq!(
            copied.molecule_properties.prop("ordinary"),
            Some(&cosmolkit_model::PropertyValue::String("keep-me".into()))
        );
        assert!(
            !copied
                .molecule_properties
                .is_prop_computed("molecule-cache")
                .unwrap()
        );
        assert_eq!(copied.molecule_properties.prop("molecule-cache"), None);
        assert_eq!(
            copied.molecule_properties.sdf_data_fields(),
            &[("FIELD".into(), "source-value".into())]
        );
        assert_eq!(copied.molecule_properties.sdf_property_lists().len(), 2);
        let atom_values = copied.molecule_properties.sdf_property_lists()[0].values();
        assert_eq!(
            atom_values,
            &[
                Some(cosmolkit_model::PropertyValue::String(
                    (format!("atom-value-{first_atom}")).into()
                )),
                Some(cosmolkit_model::PropertyValue::String(
                    (format!("atom-value-{}", first_atom + 1)).into()
                ))
            ]
        );
        let bond_values = copied.molecule_properties.sdf_property_lists()[1].values();
        assert_eq!(
            bond_values,
            &[Some(cosmolkit_model::PropertyValue::String(
                (format!("bond-value-{selected_component}")).into()
            ))]
        );
        assert_eq!(source.validate(), Ok(()));
        assert_eq!(
            source_coordinates.validate_for_atom_count(source.atoms.len()),
            Ok(())
        );
        assert_eq!(
            source_properties.name().map(super::fixed_property_text),
            Some("source-molecule")
        );
    }

    #[test]
    fn cf3d_frag_f15_noop_full_copy_keeps_identity_and_computed_state() {
        let (source, coordinates, properties) = input(2);
        let original = (source.clone(), coordinates.clone(), properties.clone());
        let copied = copy_full_molecule_remove_atoms_outside_component(
            &source,
            &coordinates,
            &properties,
            &mask(source.atoms.len(), None),
        )
        .unwrap();

        assert_eq!(copied.topology, source);
        assert_eq!(copied.coordinates, coordinates);
        assert_eq!(copied.molecule_properties, properties);
        assert_eq!(copied.mapping, TopologyMapping::identity(4, 2));
        assert!(copied.topology.atoms.iter().all(|atom| {
            atom.prop("atom-cache").is_some()
                && atom
                    .computed_prop_names()
                    .unwrap()
                    .expect("computed StringVector exists")
                    .contains(&"atom-cache".into())
        }));
        assert!(copied.topology.bonds.iter().all(|bond| {
            bond.prop("bond-cache").is_some()
                && bond
                    .computed_prop_names()
                    .unwrap()
                    .expect("computed StringVector exists")
                    .contains(&"bond-cache".into())
        }));
        assert!(
            copied
                .molecule_properties
                .is_prop_computed("molecule-cache")
                .unwrap()
        );
        assert_eq!((source, coordinates, properties), original);
    }

    #[test]
    fn cf3d_frag_f15_slow_copy_deletes_complement_for_two_three_and_four_components() {
        for component_count in [2, 3, 4] {
            let (source, coordinates, properties) = input(component_count);
            let original = (source.clone(), coordinates.clone(), properties.clone());
            let selected_component = component_count / 2;
            let copied = copy_full_molecule_remove_atoms_outside_component(
                &source,
                &coordinates,
                &properties,
                &mask(source.atoms.len(), Some(selected_component)),
            )
            .unwrap();

            assert_component_copy(
                &source,
                &coordinates,
                &properties,
                &copied,
                component_count,
                selected_component,
            );
            assert_eq!((source, coordinates, properties), original);
        }
    }
}

#[cfg(test)]
mod cf3d_frag_f16_tests {
    use super::{FullCopyComponent, copy_single_full_molecule_component};
    use crate::sanitize::{SanitizeError, SanitizeParams, SanitizeStage, sanitize_topology};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, Conformer3D,
        CoordinateBlock, CoordinateDimension, MoleculeProperties, StereoGroup, StereoGroupKind,
        SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock, TopologyMapping,
    };
    use cosmolkit_types::Element;

    fn coordinates(atom_count: usize) -> CoordinateBlock {
        let conformers_2d = [41, 5]
            .into_iter()
            .map(|id| {
                Conformer2D::new(
                    id,
                    (0..atom_count)
                        .map(|row| [row as f64 + 0.5, row as f64 + id as f64])
                        .collect(),
                )
                .with_prop("source-conformer", format!("2d-{id}"))
            })
            .collect();
        let conformers_3d = [(73, false), (11, true)]
            .into_iter()
            .map(|(id, is_3d)| {
                Conformer3D::new(
                    id,
                    (0..atom_count)
                        .map(|row| [row as f64, row as f64 + 10.0, row as f64 + 20.0])
                        .collect(),
                    is_3d,
                )
                .with_prop("source-conformer", format!("3d-{id}"))
            })
            .collect();
        CoordinateBlock {
            conformers_2d,
            conformers_3d,
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            source_conformer_order: None,
        }
    }

    fn molecule_properties() -> MoleculeProperties {
        MoleculeProperties::default()
            .with_name("source-molecule")
            .with_prop("ordinary-molecule", "retained")
            .unwrap()
            .with_computed_prop("computed-molecule", "retained-cache")
            .unwrap()
            .with_sdf_data_field("SOURCE", "full-copy")
    }

    fn singleton_source() -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
        let atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("ordinary-atom", "singleton")
                .unwrap()
                .with_computed_prop("computed-atom", "singleton-cache")
                .unwrap(),
        );
        let group = SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms(vec![AtomId::new(0)])
        .with_label("singleton-group");
        let stereo_group = StereoGroup::new(StereoGroupKind::Or, vec![AtomId::new(0)], vec![])
            .with_id(3)
            .with_write_id(9);
        let topology =
            TopologyBlock::try_from_parts(vec![atom], vec![], vec![group], vec![stereo_group])
                .expect("fixed zero-bond singleton topology must be structurally valid");
        (topology, coordinates(1), molecule_properties())
    }

    fn connected_overvalent_source() -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
        let atoms = [Element::O, Element::C, Element::C, Element::C]
            .into_iter()
            .enumerate()
            .map(|(index, element)| {
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(element)
                        .with_prop("ordinary-atom", format!("atom-{index}"))
                        .unwrap()
                        .with_computed_prop("computed-atom", format!("atom-cache-{index}"))
                        .unwrap(),
                )
            })
            .collect();
        let bonds = (0..3)
            .map(|index| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(0), AtomId::new(index + 1), BondOrder::Single)
                        .with_prop("ordinary-bond", format!("bond-{index}"))
                        .unwrap()
                        .with_computed_prop("computed-bond", format!("bond-cache-{index}"))
                        .unwrap(),
                )
            })
            .collect();
        let group = SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms((0..4).map(AtomId::new).collect())
        .with_bonds((0..3).map(BondId::new).collect())
        .with_label("connected-group");
        let stereo_group = StereoGroup::new(
            StereoGroupKind::And,
            vec![AtomId::new(0), AtomId::new(1)],
            vec![BondId::new(0)],
        )
        .with_id(17)
        .with_write_id(29);
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![group], vec![stereo_group])
            .expect("fixed connected graph must be structurally valid before sanitation");
        (topology, coordinates(4), molecule_properties())
    }

    fn assert_identity_mapping(copied: &FullCopyComponent, atom_count: usize, bond_count: usize) {
        assert_eq!(
            copied.mapping,
            TopologyMapping::identity(atom_count, bond_count)
        );
        assert_eq!(
            copied.mapping.atoms().old_to_new(),
            &(0..atom_count)
                .map(|row| Some(AtomId::new(row)))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            copied.mapping.atoms().new_to_old(),
            &(0..atom_count)
                .map(|row| Some(AtomId::new(row)))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            copied.mapping.bonds().old_to_new(),
            &(0..bond_count)
                .map(|row| Some(BondId::new(row)))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            copied.mapping.bonds().new_to_old(),
            &(0..bond_count)
                .map(|row| Some(BondId::new(row)))
                .collect::<Vec<_>>()
        );
    }

    fn assert_all_conformers(copied: &FullCopyComponent, atom_count: usize) {
        assert_eq!(copied.coordinates.conformers_2d.len(), 2);
        assert_eq!(
            copied
                .coordinates
                .conformers_2d
                .iter()
                .map(Conformer2D::id)
                .collect::<Vec<_>>(),
            vec![41, 5]
        );
        assert_eq!(copied.coordinates.conformers_3d.len(), 2);
        assert_eq!(
            copied
                .coordinates
                .conformers_3d
                .iter()
                .map(|conformer| (conformer.id(), conformer.is_3d()))
                .collect::<Vec<_>>(),
            vec![(73, false), (11, true)]
        );
        assert!(
            copied
                .coordinates
                .conformers_2d
                .iter()
                .all(|conformer| conformer.coordinates().len() == atom_count)
        );
        assert!(
            copied
                .coordinates
                .conformers_3d
                .iter()
                .all(|conformer| conformer.coordinates().len() == atom_count)
        );
        assert!(copied.coordinates.conformers_2d.iter().all(|conformer| {
            conformer
                .props()
                .get("source-conformer".as_bytes())
                .is_some()
        }));
        assert!(copied.coordinates.conformers_3d.iter().all(|conformer| {
            conformer
                .props()
                .get("source-conformer".as_bytes())
                .is_some()
        }));
        assert_eq!(
            copied.coordinates.source_coordinate_dim,
            Some(CoordinateDimension::ThreeD)
        );
    }

    fn assert_molecule_properties(copied: &FullCopyComponent) {
        assert_eq!(
            copied
                .molecule_properties
                .name()
                .map(super::fixed_property_text),
            Some("source-molecule")
        );
        assert_eq!(
            copied
                .molecule_properties
                .props()
                .get("ordinary-molecule".as_bytes())
                .map(|v| v.as_string().expect("fixed StringTag"))
                .map(super::fixed_property_text),
            Some("retained")
        );
        assert!(
            copied
                .molecule_properties
                .is_prop_computed("computed-molecule")
                .unwrap()
        );
        assert_eq!(
            copied
                .molecule_properties
                .sdf_data_fields()
                .iter()
                .find(|(key, _)| key.as_bytes() == b"SOURCE")
                .map(|(_, value)| super::fixed_property_text(value)),
            Some("full-copy")
        );
    }

    #[test]
    fn cf3d_frag_f16_zero_bond_singleton_full_copy_preserves_state() {
        let (source_topology, source_coordinates, source_properties) = singleton_source();
        let original = (
            source_topology.clone(),
            source_coordinates.clone(),
            source_properties.clone(),
        );

        let copied = copy_single_full_molecule_component(
            &source_topology,
            &source_coordinates,
            &source_properties,
        );

        assert_eq!(copied.topology, source_topology);
        assert_eq!(copied.coordinates, source_coordinates);
        assert_eq!(copied.molecule_properties, source_properties);
        assert_identity_mapping(&copied, 1, 0);
        assert_all_conformers(&copied, 1);
        assert_molecule_properties(&copied);
        assert_eq!(copied.topology.substance_groups.len(), 1);
        assert_eq!(copied.topology.stereo_groups[0].id(), Some(3));
        assert_eq!(copied.topology.stereo_groups[0].write_id(), 9);
        assert_eq!(
            copied.topology.atoms[0].prop("ordinary-atom"),
            Some(&cosmolkit_model::PropertyValue::String(
                ("singleton".to_owned()).into()
            ))
        );
        assert!(
            copied.topology.atoms[0]
                .computed_prop_names()
                .unwrap()
                .expect("computed StringVector exists")
                .contains(&"computed-atom".into())
        );
        assert_eq!(
            (source_topology, source_coordinates, source_properties),
            original
        );
    }

    #[test]
    fn cf3d_frag_f16_connected_copy_defers_source_sanitation_failure() {
        let (source_topology, source_coordinates, source_properties) =
            connected_overvalent_source();
        let original = (
            source_topology.clone(),
            source_coordinates.clone(),
            source_properties.clone(),
        );

        let copied = copy_single_full_molecule_component(
            &source_topology,
            &source_coordinates,
            &source_properties,
        );

        assert_eq!(copied.topology, source_topology);
        assert_eq!(copied.coordinates, source_coordinates);
        assert_eq!(copied.molecule_properties, source_properties);
        assert_identity_mapping(&copied, 4, 3);
        assert_all_conformers(&copied, 4);
        assert_molecule_properties(&copied);
        assert_eq!(copied.topology.substance_groups.len(), 1);
        assert_eq!(copied.topology.stereo_groups[0].id(), Some(17));
        assert_eq!(copied.topology.stereo_groups[0].write_id(), 29);
        assert!(copied.topology.bonds.iter().all(|bond| {
            bond.computed_prop_names()
                .unwrap()
                .expect("computed StringVector exists")
                .contains(&"computed-bond".into())
                && bond.prop("ordinary-bond").is_some()
        }));

        assert!(matches!(
            sanitize_topology(&copied.topology, &SanitizeParams::default()),
            Err(SanitizeError::Properties {
                stage: SanitizeStage::Properties,
                ..
            })
        ));
        assert_eq!(
            (source_topology, source_coordinates, source_properties),
            original
        );
    }
}

#[cfg(test)]
mod cf3d_frag_f17_tests {
    use super::{
        AtomPathSubsetCopyError, FullCopyComponentError, MoleculeFragmentsFailure,
        OrderedFragmentBuildError, build_ordered_fragment_copies, get_molecule_fragments,
    };
    use crate::sanitize::{SanitizeError, SanitizeStage};
    use cosmolkit_model::{
        Atom, AtomId, AtomMapping, AtomSpec, Bond, BondId, BondMapping, BondOrder, BondSpec,
        Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, MoleculeProperties,
        SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock, TopologyMapping,
    };
    use cosmolkit_types::{BondStereo, Element};

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum ComponentShape {
        Singleton,
        PlainDimer,
        ChallengingDimer,
    }

    #[derive(Clone, Debug)]
    struct ExpectedComponent {
        atoms: Vec<AtomId>,
        bonds: Vec<BondId>,
        challenging: bool,
    }

    fn coordinates(atom_count: usize) -> CoordinateBlock {
        CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(
                    41,
                    (0..atom_count)
                        .map(|row| [row as f64 + 0.25, row as f64 + 0.75])
                        .collect(),
                )
                .with_prop("source-conformer", "2d-source"),
            ],
            conformers_3d: vec![
                Conformer3D::new(
                    73,
                    (0..atom_count)
                        .map(|row| [row as f64, row as f64 + 1.0, row as f64 + 2.0])
                        .collect(),
                    true,
                )
                .with_prop("source-conformer", "3d-source"),
            ],
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            source_conformer_order: None,
        }
    }

    fn molecule_properties() -> MoleculeProperties {
        MoleculeProperties::default()
            .with_name("fragment-source")
            .with_prop("ordinary-molecule", "kept-by-full-copy")
            .unwrap()
            .with_computed_prop("computed-molecule", "source-cache")
            .unwrap()
    }

    fn source_for(
        shapes: &[ComponentShape],
    ) -> (
        TopologyBlock,
        CoordinateBlock,
        MoleculeProperties,
        Vec<ExpectedComponent>,
    ) {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut expected = Vec::new();
        for shape in shapes {
            let atom_start = atoms.len();
            let atom_count = usize::from(*shape != ComponentShape::Singleton) + 1;
            let component_atoms = (atom_start..atom_start + atom_count)
                .map(AtomId::new)
                .collect::<Vec<_>>();
            for atom_index in atom_start..atom_start + atom_count {
                let atom = Atom::from_spec(
                    AtomId::new(atom_index),
                    AtomSpec::new(Element::C)
                        .with_prop("source-atom", format!("atom-{atom_index}"))
                        .unwrap()
                        .with_computed_prop("computed-atom", format!("atom-cache-{atom_index}"))
                        .unwrap(),
                );
                atoms.push(atom);
            }
            let mut component_bonds = Vec::new();
            if atom_count == 2 {
                let bond_id = BondId::new(bonds.len());
                let challenging = *shape == ComponentShape::ChallengingDimer;
                let mut spec = BondSpec::new(
                    component_atoms[0],
                    component_atoms[1],
                    if challenging {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    },
                )
                .with_prop("source-bond", format!("bond-{}", bond_id.index()))
                .unwrap()
                .with_computed_prop("computed-bond", format!("bond-cache-{}", bond_id.index()))
                .unwrap();
                if challenging {
                    spec = spec
                        .with_stereo(BondStereo::Cis)
                        .with_stereo_atoms(component_atoms[0], component_atoms[1]);
                }
                bonds.push(Bond::from_spec(bond_id, spec));
                component_bonds.push(bond_id);
            }
            expected.push(ExpectedComponent {
                atoms: component_atoms,
                bonds: component_bonds,
                challenging: *shape == ComponentShape::ChallengingDimer,
            });
        }
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])
            .expect("fixed F17 components must be structurally valid");
        let coordinates = coordinates(topology.atoms.len());
        (topology, coordinates, molecule_properties(), expected)
    }

    fn expected_mapping(
        atom_count: usize,
        bond_count: usize,
        component: &ExpectedComponent,
    ) -> TopologyMapping {
        let mut atom_old_to_new = vec![None; atom_count];
        for (new_row, old_atom) in component.atoms.iter().copied().enumerate() {
            atom_old_to_new[old_atom.index()] = Some(AtomId::new(new_row));
        }
        let mut bond_old_to_new = vec![None; bond_count];
        for (new_row, old_bond) in component.bonds.iter().copied().enumerate() {
            bond_old_to_new[old_bond.index()] = Some(BondId::new(new_row));
        }
        TopologyMapping {
            atoms: AtomMapping {
                old_to_new: atom_old_to_new,
                new_to_old: component.atoms.iter().copied().map(Some).collect(),
            },
            bonds: BondMapping {
                old_to_new: bond_old_to_new,
                new_to_old: component.bonds.iter().copied().map(Some).collect(),
            },
        }
    }

    fn shapes_for(count: usize, pattern: usize) -> Vec<ComponentShape> {
        match pattern {
            0 => vec![ComponentShape::Singleton; count],
            1 => vec![ComponentShape::PlainDimer; count],
            2 => {
                let mut shapes = vec![ComponentShape::PlainDimer; count];
                if let Some(first) = shapes.first_mut() {
                    *first = ComponentShape::ChallengingDimer;
                }
                shapes
            }
            _ => unreachable!("the F17 fixture has exactly three finite patterns"),
        }
    }

    fn assert_successful_case(
        component_count: usize,
        shapes: &[ComponentShape],
        sanitize: bool,
        copy_conformers: bool,
        source: &TopologyBlock,
        source_coordinates: &CoordinateBlock,
        source_properties: &MoleculeProperties,
        expected: &[ExpectedComponent],
    ) {
        let original = (
            source.clone(),
            source_coordinates.clone(),
            source_properties.clone(),
        );
        let copied = build_ordered_fragment_copies(
            source,
            source_coordinates,
            source_properties,
            sanitize,
            copy_conformers,
        )
        .expect("fixed route matrix should build every component");
        assert_eq!(copied.len(), component_count);
        assert_eq!(copied.len(), expected.len());

        for (component_index, (fragment, expected)) in copied.iter().zip(expected).enumerate() {
            let fast_subset = component_count > 1
                && (expected.atoms.len() == 1 || (component_count > 3 && !expected.challenging));
            let coordinates_copied = !fast_subset || copy_conformers;
            let retained_full_properties = component_count == 1 || !fast_subset;
            assert_eq!(fragment.component_atoms, expected.atoms);
            assert_eq!(
                fragment.copy.mapping,
                expected_mapping(source.atoms.len(), source.bonds.len(), expected)
            );
            assert_eq!(fragment.copy.topology.atoms.len(), expected.atoms.len());
            assert_eq!(fragment.copy.topology.bonds.len(), expected.bonds.len());
            for (new_row, old_atom) in expected.atoms.iter().copied().enumerate() {
                let copied_atom = &fragment.copy.topology.atoms[new_row];
                assert_eq!(copied_atom.id(), AtomId::new(new_row));
                assert_eq!(
                    copied_atom.prop("source-atom"),
                    Some(&cosmolkit_model::PropertyValue::String(
                        (format!("atom-{}", old_atom.index())).into()
                    ))
                );
                assert_eq!(
                    copied_atom.is_prop_computed("computed-atom").unwrap(),
                    component_count == 1
                );
            }
            for (new_row, old_bond) in expected.bonds.iter().copied().enumerate() {
                let copied_bond = &fragment.copy.topology.bonds[new_row];
                let source_bond = &source.bonds[old_bond.index()];
                assert_eq!(copied_bond.id(), BondId::new(new_row));
                assert_eq!(copied_bond.stereo(), source_bond.stereo());
                assert_eq!(
                    copied_bond.prop("source-bond"),
                    Some(&cosmolkit_model::PropertyValue::String(
                        (format!("bond-{}", old_bond.index())).into()
                    ))
                );
                assert_eq!(
                    copied_bond.is_prop_computed("computed-bond").unwrap(),
                    component_count == 1
                );
            }

            assert_eq!(
                fragment
                    .copy
                    .molecule_properties
                    .name()
                    .map(super::fixed_property_text),
                retained_full_properties.then_some("fragment-source")
            );
            assert_eq!(
                fragment
                    .copy
                    .molecule_properties
                    .props()
                    .get("ordinary-molecule".as_bytes())
                    .map(|v| v.as_string().expect("fixed StringTag"))
                    .map(super::fixed_property_text),
                retained_full_properties.then_some("kept-by-full-copy")
            );
            assert_eq!(
                fragment
                    .copy
                    .molecule_properties
                    .is_prop_computed("computed-molecule")
                    .unwrap(),
                component_count == 1
            );

            if coordinates_copied {
                let expected_2d = expected
                    .atoms
                    .iter()
                    .map(|atom| source_coordinates.conformers_2d[0].coordinates()[atom.index()])
                    .collect::<Vec<_>>();
                let expected_3d = expected
                    .atoms
                    .iter()
                    .map(|atom| source_coordinates.conformers_3d[0].coordinates()[atom.index()])
                    .collect::<Vec<_>>();
                assert_eq!(fragment.copy.coordinates.conformers_2d[0].id(), 41);
                assert_eq!(
                    fragment.copy.coordinates.conformers_2d[0].coordinates(),
                    expected_2d
                );
                assert_eq!(
                    fragment.copy.coordinates.conformers_2d[0]
                        .props()
                        .get("source-conformer".as_bytes())
                        .map(super::fixed_property_text),
                    (!fast_subset).then_some("2d-source")
                );
                assert_eq!(fragment.copy.coordinates.conformers_3d[0].id(), 73);
                assert!(fragment.copy.coordinates.conformers_3d[0].is_3d());
                assert_eq!(
                    fragment.copy.coordinates.conformers_3d[0].coordinates(),
                    expected_3d
                );
                assert_eq!(
                    fragment.copy.coordinates.conformers_3d[0]
                        .props()
                        .get("source-conformer".as_bytes())
                        .map(super::fixed_property_text),
                    (!fast_subset).then_some("3d-source")
                );
                assert_eq!(
                    fragment.copy.coordinates.source_coordinate_dim,
                    Some(CoordinateDimension::ThreeD)
                );
            } else {
                assert!(fragment.copy.coordinates.conformers_2d.is_empty());
                assert!(fragment.copy.coordinates.conformers_3d.is_empty());
                assert_eq!(fragment.copy.coordinates.source_coordinate_dim, None);
            }
            assert_eq!(fragment.component_atoms, expected.atoms);
            assert_eq!(
                fragment
                    .copy
                    .topology
                    .bonds
                    .iter()
                    .any(|bond| bond.stereo() == BondStereo::Cis),
                shapes[component_index] == ComponentShape::ChallengingDimer
            );
        }
        assert_eq!(
            (
                source.clone(),
                source_coordinates.clone(),
                source_properties.clone()
            ),
            original
        );
    }

    #[test]
    fn cf3d_frag_f17_count_shape_challenge_route_order_and_maps() {
        let (empty_topology, empty_coordinates, empty_properties, empty_expected) = source_for(&[]);
        assert_successful_case(
            0,
            &[],
            false,
            false,
            &empty_topology,
            &empty_coordinates,
            &empty_properties,
            &empty_expected,
        );

        for component_count in 1..=5 {
            for pattern in 0..3 {
                let shapes = shapes_for(component_count, pattern);
                let (topology, coordinates, properties, expected) = source_for(&shapes);
                for sanitize in [false, true] {
                    for copy_conformers in [false, true] {
                        assert_successful_case(
                            component_count,
                            &shapes,
                            sanitize,
                            copy_conformers,
                            &topology,
                            &coordinates,
                            &properties,
                            &expected,
                        );
                    }
                }
            }
        }
    }

    fn source_with_early_fast_sanitize_failure()
    -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let components = [
            (Element::O, 3),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 2),
        ];
        for (element, component_size) in components {
            let start = atoms.len();
            let atom_count = if component_size == 0 {
                1
            } else {
                component_size + 1
            };
            for atom_index in start..start + atom_count {
                atoms.push(Atom::from_spec(
                    AtomId::new(atom_index),
                    AtomSpec::new(if atom_index == 0 { element } else { Element::C }),
                ));
            }
            if element == Element::O {
                for atom_index in start + 1..start + atom_count {
                    let bond_id = BondId::new(bonds.len());
                    bonds.push(Bond::from_spec(
                        bond_id,
                        BondSpec::new(
                            AtomId::new(start),
                            AtomId::new(atom_index),
                            BondOrder::Single,
                        ),
                    ));
                }
            } else if component_size == 2 {
                let bond_id = BondId::new(bonds.len());
                bonds.push(Bond::from_spec(
                    bond_id,
                    BondSpec::new(
                        AtomId::new(start),
                        AtomId::new(start + 1),
                        BondOrder::Double,
                    )
                    .with_stereo(BondStereo::Cis)
                    .with_stereo_atoms(AtomId::new(start), AtomId::new(start + 1)),
                ));
            }
        }
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])
            .expect("fixed failure-order graph must be structurally valid");
        let coordinates = coordinates(topology.atoms.len());
        (topology, coordinates, molecule_properties())
    }

    #[test]
    fn cf3d_frag_f17_early_fast_subset_error_precedes_later_slow_component() {
        let (topology, coordinates, properties) = source_with_early_fast_sanitize_failure();
        let error = build_ordered_fragment_copies(&topology, &coordinates, &properties, true, true)
            .unwrap_err();

        assert!(matches!(
            error,
            OrderedFragmentBuildError::FastSubset {
                component_index: 0,
                source: AtomPathSubsetCopyError::Sanitize(SanitizeError::Properties {
                    stage: SanitizeStage::Properties,
                    ..
                }),
            }
        ));
    }

    #[test]
    fn cf3d_frag_f17_later_slow_error_follows_fast_singletons() {
        let shapes = [
            ComponentShape::Singleton,
            ComponentShape::Singleton,
            ComponentShape::Singleton,
            ComponentShape::Singleton,
            ComponentShape::ChallengingDimer,
        ];
        let (topology, mut coordinates, properties, _) = source_for(&shapes);
        coordinates.conformers_2d[0] =
            Conformer2D::new(41, vec![[0.0, 0.0]; topology.atoms.len() - 1]);

        let error =
            build_ordered_fragment_copies(&topology, &coordinates, &properties, false, false)
                .unwrap_err();

        assert!(matches!(
            error,
            OrderedFragmentBuildError::SlowFullCopy {
                component_index: 4,
                source: FullCopyComponentError::CoordinateValidation(_),
            }
        ));
    }

    #[test]
    fn cf3d_frag_f18_options_return_ordered_components_maps_and_copied_state() {
        for component_count in 0..=5 {
            let pattern_count = if component_count == 0 { 1 } else { 3 };
            for pattern in 0..pattern_count {
                let shapes = shapes_for(component_count, pattern);
                let (topology, coordinates, properties, expected) = source_for(&shapes);
                let original = (topology.clone(), coordinates.clone(), properties.clone());

                for sanitize in [false, true] {
                    for copy_conformers in [false, true] {
                        let fragments = get_molecule_fragments(
                            &topology,
                            &coordinates,
                            &properties,
                            sanitize,
                            copy_conformers,
                        )
                        .expect("fixed F18 option matrix must succeed");
                        assert_eq!(fragments.len(), component_count);

                        for (component_index, (fragment, expected)) in
                            fragments.iter().zip(&expected).enumerate()
                        {
                            let fast_subset = component_count > 1
                                && (expected.atoms.len() == 1
                                    || (component_count > 3 && !expected.challenging));
                            let retained_full_properties = component_count == 1 || !fast_subset;
                            let retained_computed_properties = component_count == 1 && !sanitize;
                            assert_eq!(fragment.component_atoms(), expected.atoms);
                            assert_eq!(
                                fragment.topology_mapping(),
                                &expected_mapping(
                                    topology.atoms.len(),
                                    topology.bonds.len(),
                                    expected,
                                )
                            );
                            assert_eq!(fragment.topology().atoms.len(), expected.atoms.len());
                            assert_eq!(fragment.topology().bonds.len(), expected.bonds.len());
                            for (new_row, old_atom) in expected.atoms.iter().copied().enumerate() {
                                let atom = &fragment.topology().atoms[new_row];
                                assert_eq!(atom.id(), AtomId::new(new_row));
                                assert_eq!(
                                    atom.prop("source-atom"),
                                    Some(&cosmolkit_model::PropertyValue::String(
                                        (format!("atom-{}", old_atom.index())).into()
                                    ))
                                );
                                assert_eq!(
                                    atom.is_prop_computed("computed-atom").unwrap(),
                                    retained_computed_properties
                                );
                            }
                            for (new_row, old_bond) in expected.bonds.iter().copied().enumerate() {
                                let bond = &fragment.topology().bonds[new_row];
                                assert_eq!(bond.id(), BondId::new(new_row));
                                assert_eq!(
                                    bond.prop("source-bond"),
                                    Some(&cosmolkit_model::PropertyValue::String(
                                        (format!("bond-{}", old_bond.index())).into()
                                    ))
                                );
                                assert_eq!(
                                    bond.is_prop_computed("computed-bond").unwrap(),
                                    retained_computed_properties
                                );
                            }

                            assert_eq!(
                                fragment
                                    .molecule_properties()
                                    .name()
                                    .map(super::fixed_property_text),
                                retained_full_properties.then_some("fragment-source")
                            );
                            assert_eq!(
                                fragment.molecule_properties().prop("ordinary-molecule"),
                                retained_full_properties.then_some(
                                    &cosmolkit_model::PropertyValue::String(
                                        "kept-by-full-copy".into()
                                    )
                                )
                            );
                            assert_eq!(
                                fragment
                                    .molecule_properties()
                                    .is_prop_computed("computed-molecule")
                                    .unwrap(),
                                retained_computed_properties
                            );

                            assert_eq!(
                                fragment.coordinates().conformers_2d.len(),
                                usize::from(copy_conformers)
                            );
                            assert_eq!(
                                fragment.coordinates().conformers_3d.len(),
                                usize::from(copy_conformers)
                            );
                            let expected_source_dim = (!fast_subset || copy_conformers)
                                .then_some(CoordinateDimension::ThreeD);
                            assert_eq!(
                                fragment.coordinates().source_coordinate_dim,
                                expected_source_dim
                            );
                            if copy_conformers {
                                let expected_2d = expected
                                    .atoms
                                    .iter()
                                    .map(|atom| {
                                        coordinates.conformers_2d[0].coordinates()[atom.index()]
                                    })
                                    .collect::<Vec<_>>();
                                let expected_3d = expected
                                    .atoms
                                    .iter()
                                    .map(|atom| {
                                        coordinates.conformers_3d[0].coordinates()[atom.index()]
                                    })
                                    .collect::<Vec<_>>();
                                assert_eq!(fragment.coordinates().conformers_2d[0].id(), 41);
                                assert_eq!(
                                    fragment.coordinates().conformers_2d[0].coordinates(),
                                    expected_2d
                                );
                                assert_eq!(fragment.coordinates().conformers_3d[0].id(), 73);
                                assert_eq!(
                                    fragment.coordinates().conformers_3d[0].coordinates(),
                                    expected_3d
                                );
                            }
                            assert_eq!(
                                fragment
                                    .topology()
                                    .bonds
                                    .iter()
                                    .any(|bond| { bond.stereo() == BondStereo::Cis }),
                                shapes[component_index] == ComponentShape::ChallengingDimer
                            );
                        }
                        assert_eq!(
                            (topology.clone(), coordinates.clone(), properties.clone()),
                            original
                        );
                    }
                }
            }
        }
    }

    fn alternating_benzene_source(
        singleton_count: usize,
    ) -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
        let atom_count = 6 + singleton_count;
        let atoms = (0..atom_count)
            .map(|row| {
                Atom::from_spec(
                    AtomId::new(row),
                    AtomSpec::new(Element::C)
                        .with_no_implicit(true)
                        .with_prop("source-atom", format!("atom-{row}"))
                        .unwrap()
                        .with_computed_prop("computed-atom", format!("atom-cache-{row}"))
                        .unwrap(),
                )
            })
            .collect();
        let bonds = (0..6)
            .map(|row| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(
                        AtomId::new(row),
                        AtomId::new((row + 1) % 6),
                        if row % 2 == 0 {
                            BondOrder::Double
                        } else {
                            BondOrder::Single
                        },
                    )
                    .with_prop("source-bond", format!("bond-{row}"))
                    .unwrap()
                    .with_computed_prop("computed-bond", format!("bond-cache-{row}"))
                    .unwrap(),
                )
            })
            .collect();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])
            .expect("fixed alternating benzene topology must be valid");
        (topology, coordinates(atom_count), molecule_properties())
    }

    #[test]
    fn cf3d_frag_f18_final_sanitize_normalizes_full_copy_and_clears_cache() {
        let (topology, coordinates, properties) = alternating_benzene_source(0);
        let original = (topology.clone(), coordinates.clone(), properties.clone());

        let unsanitized =
            get_molecule_fragments(&topology, &coordinates, &properties, false, true).unwrap();
        assert_eq!(unsanitized.len(), 1);
        assert!(
            unsanitized[0]
                .topology()
                .atoms
                .iter()
                .all(|atom| !atom.is_aromatic())
        );
        assert!(
            unsanitized[0]
                .molecule_properties()
                .is_prop_computed("computed-molecule")
                .unwrap()
        );

        let sanitized =
            get_molecule_fragments(&topology, &coordinates, &properties, true, true).unwrap();
        assert_eq!(sanitized.len(), 1);
        assert!(sanitized[0].topology().atoms.iter().all(Atom::is_aromatic));
        assert!(sanitized[0].topology().bonds.iter().all(Bond::is_aromatic));
        assert!(
            !sanitized[0]
                .molecule_properties()
                .is_prop_computed("computed-molecule")
                .unwrap()
        );
        assert!(
            sanitized[0]
                .topology()
                .atoms
                .iter()
                .all(|atom| !atom.is_prop_computed("computed-atom").unwrap())
        );
        assert!(
            sanitized[0]
                .topology()
                .bonds
                .iter()
                .all(|bond| !bond.is_prop_computed("computed-bond").unwrap())
        );
        assert_eq!(
            (topology, coordinates, properties),
            original,
            "the detached source must remain unchanged"
        );
    }

    #[test]
    fn cf3d_frag_f18_outer_sanitize_follows_fast_subset_sanitize() {
        let (topology, coordinates, properties) = alternating_benzene_source(3);
        let original = (topology.clone(), coordinates.clone(), properties.clone());
        let fragments =
            get_molecule_fragments(&topology, &coordinates, &properties, true, false).unwrap();

        assert_eq!(fragments.len(), 4);
        assert_eq!(
            fragments[0].component_atoms(),
            &(0..6).map(AtomId::new).collect::<Vec<_>>()
        );
        assert!(fragments[0].topology().atoms.iter().all(Atom::is_aromatic));
        assert!(fragments[0].topology().bonds.iter().all(Bond::is_aromatic));
        assert!(fragments.iter().all(|fragment| {
            fragment.coordinates().conformers_2d.is_empty()
                && fragment.coordinates().conformers_3d.is_empty()
        }));
        assert_eq!(
            (topology, coordinates, properties),
            original,
            "the detached source must remain unchanged"
        );
    }

    fn earlier_slow_sanitize_failure_later_fast_failure_source()
    -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::O)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(3), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(4), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(5), AtomSpec::new(Element::C)),
            Atom::from_spec(
                AtomId::new(6),
                AtomSpec::new(Element::O).with_explicit_hydrogens(3),
            ),
        ];
        let bonds = (1..4)
            .map(|row| {
                Bond::from_spec(
                    BondId::new(row - 1),
                    BondSpec::new(AtomId::new(0), AtomId::new(row), BondOrder::Single),
                )
            })
            .collect();
        let group = SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms((0..4).map(AtomId::new).collect())
        .with_bonds((0..3).map(BondId::new).collect());
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![group], vec![])
            .expect("fixed competing-sanitize topology must be structurally valid");
        (topology, coordinates(7), molecule_properties())
    }

    #[test]
    fn cf3d_frag_f18_build_errors_precede_outer_sanitize_errors() {
        let (topology, coordinates, properties) =
            earlier_slow_sanitize_failure_later_fast_failure_source();
        let error =
            get_molecule_fragments(&topology, &coordinates, &properties, true, true).unwrap_err();

        assert!(matches!(
            &error.failure,
            MoleculeFragmentsFailure::Build(OrderedFragmentBuildError::FastSubset {
                component_index: 3,
                source: AtomPathSubsetCopyError::Sanitize(SanitizeError::Properties {
                    stage: SanitizeStage::Properties,
                    ..
                }),
            })
        ));
    }

    fn two_distinct_final_sanitize_failures(
        aromatic_first: bool,
    ) -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
        #[derive(Clone, Copy)]
        enum Component {
            OvervalentOxygen,
            AromaticAtomOutsideRing,
        }
        let components = if aromatic_first {
            [
                Component::AromaticAtomOutsideRing,
                Component::OvervalentOxygen,
            ]
        } else {
            [
                Component::OvervalentOxygen,
                Component::AromaticAtomOutsideRing,
            ]
        };
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        for component in components {
            let start = atoms.len();
            match component {
                Component::OvervalentOxygen => {
                    atoms.push(Atom::from_spec(
                        AtomId::new(start),
                        AtomSpec::new(Element::O),
                    ));
                    for row in start + 1..start + 4 {
                        atoms.push(Atom::from_spec(AtomId::new(row), AtomSpec::new(Element::C)));
                        bonds.push(Bond::from_spec(
                            BondId::new(bonds.len()),
                            BondSpec::new(AtomId::new(start), AtomId::new(row), BondOrder::Single),
                        ));
                    }
                }
                Component::AromaticAtomOutsideRing => {
                    atoms.push(Atom::from_spec(
                        AtomId::new(start),
                        AtomSpec::new(Element::C)
                            .with_aromatic(true)
                            .with_no_implicit(true),
                    ));
                    atoms.push(Atom::from_spec(
                        AtomId::new(start + 1),
                        AtomSpec::new(Element::C).with_no_implicit(true),
                    ));
                    bonds.push(Bond::from_spec(
                        BondId::new(bonds.len()),
                        BondSpec::new(
                            AtomId::new(start),
                            AtomId::new(start + 1),
                            BondOrder::Single,
                        ),
                    ));
                }
            }
        }
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])
            .expect("fixed sanitize-failure components must be structurally valid");
        (
            topology.clone(),
            coordinates(topology.atoms.len()),
            molecule_properties(),
        )
    }

    #[test]
    fn cf3d_frag_f18_reports_the_first_of_two_distinct_sanitize_failures() {
        let (topology, coordinates, properties) = two_distinct_final_sanitize_failures(false);
        assert_eq!(topology.atoms.len(), 6);
        let source = (topology.clone(), coordinates.clone(), properties.clone());
        let error =
            get_molecule_fragments(&topology, &coordinates, &properties, true, true).unwrap_err();
        assert_eq!(error.component_index(), Some(0));
        assert!(matches!(
            &error.failure,
            MoleculeFragmentsFailure::FinalSanitize {
                component_index: 0,
                source: SanitizeError::Properties {
                    stage: SanitizeStage::Properties,
                    ..
                }
            }
        ));
        assert_eq!((topology, coordinates, properties), source);

        let (topology, coordinates, properties) = two_distinct_final_sanitize_failures(true);
        let error =
            get_molecule_fragments(&topology, &coordinates, &properties, true, true).unwrap_err();
        assert_eq!(error.component_index(), Some(0));
        assert!(matches!(
            &error.failure,
            MoleculeFragmentsFailure::FinalSanitize {
                component_index: 0,
                source: SanitizeError::Kekulize {
                    stage: SanitizeStage::Kekulize,
                    ..
                }
            }
        ));
    }

    fn ids(indices: &[usize]) -> Vec<AtomId> {
        indices.iter().copied().map(AtomId::new).collect()
    }
    fn bond_ids(indices: &[usize]) -> Vec<BondId> {
        indices.iter().copied().map(BondId::new).collect()
    }

    fn source562_run(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        sanitize: bool,
        copy_coordinates: bool,
        labels: &mut Vec<i32>,
        mapping: &mut Vec<Vec<i32>>,
        metadata: super::FragmentSourceMetadataView<'_>,
    ) -> Result<Vec<super::MoleculeFragment>, super::MoleculeFragmentsError> {
        // These new native fixtures explicitly construct the source conformer
        // vector as [2D(id41), 3D(id73)]. This is fixture input, not production
        // inference for an unknown mixed-dimension source ordering.
        assert_eq!(coordinates.conformers_2d.len(), 1);
        assert_eq!(coordinates.conformers_3d.len(), 1);
        let mut fixture_coordinates = coordinates.clone();
        fixture_coordinates.source_conformer_order =
            Some(vec![CoordinateDimension::TwoD, CoordinateDimension::ThreeD]);
        super::get_molecule_fragments_with_source_outputs(
            topology,
            &super::FragmentCoordinateView::from_coordinate_block(&fixture_coordinates),
            &MoleculeProperties::default(),
            sanitize,
            copy_coordinates,
            Some(labels),
            Some(mapping),
            metadata,
        )
    }

    #[test]
    fn source562_outputs_replace_labels_but_append_maps_in_source_order() {
        let (graph, coords, _, _) =
            source_for(&[ComponentShape::PlainDimer, ComponentShape::Singleton]);
        let mut labels = vec![77, 88];
        let mut mapping = vec![vec![91]];
        let result = source562_run(
            &graph,
            &coords,
            false,
            true,
            &mut labels,
            &mut mapping,
            super::FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap();
        assert_eq!(labels, vec![0, 0, 1]);
        assert_eq!(mapping, vec![vec![91], vec![0, 1], vec![2]]);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn source562_slow_raw_coordinates_preserve_float_bits() {
        let (graph, mut coords, _, _) =
            source_for(&[ComponentShape::PlainDimer, ComponentShape::PlainDimer]);
        let nan = f64::from_bits(0x7ff8_0000_0000_1234);
        coords.conformers_3d[0] = Conformer3D::new(
            73,
            vec![
                [nan, -0.0, f64::INFINITY],
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
            ],
            true,
        );
        let result = source562_run(
            &graph,
            &coords,
            false,
            true,
            &mut vec![],
            &mut vec![],
            super::FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap();
        let row = result[0].coordinates().conformers_3d[0].coordinates()[0];
        assert_eq!(row[0].to_bits(), nan.to_bits());
        assert_eq!(row[1].to_bits(), (-0.0f64).to_bits());
        assert_eq!(row[2].to_bits(), f64::INFINITY.to_bits());
    }

    #[test]
    fn source562_deleted_missing_tail_is_allowed_until_a_retained_row_is_lost() {
        let (graph, mut coords, _, _) =
            source_for(&[ComponentShape::PlainDimer, ComponentShape::PlainDimer]);
        coords.conformers_2d[0] = Conformer2D::new(41, vec![[0.0, 0.0]; 2]);
        let mut labels = vec![77];
        let mut mapping = vec![vec![91]];
        let error = source562_run(
            &graph,
            &coords,
            false,
            false,
            &mut labels,
            &mut mapping,
            super::FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap_err();
        assert_eq!(error.component_index(), Some(1));
        assert_eq!(labels, vec![0, 0, 1, 1]);
        assert_eq!(mapping, vec![vec![91], vec![0, 1]]);
        // copyConformers=false clears copies only after all builders finish;
        // the second slow-copy coordinate invariant still executes first.
        assert!(matches!(
            error.failure,
            MoleculeFragmentsFailure::Build(OrderedFragmentBuildError::SlowFullCopy {
                component_index: 1,
                ..
            })
        ));
    }

    #[test]
    fn source562_fast_sanitize_failure_keeps_only_completed_component_maps() {
        let atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::O).with_explicit_hydrogens(3),
        );
        let graph = TopologyBlock::try_from_parts(
            vec![
                atom,
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        let mut labels = vec![77];
        let mut mapping = vec![vec![91]];
        let error = source562_run(
            &graph,
            &coordinates(2),
            true,
            true,
            &mut labels,
            &mut mapping,
            super::FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap_err();
        assert_eq!(error.component_index(), Some(0));
        assert_eq!(labels, vec![0, 1]);
        assert_eq!(mapping, vec![vec![91]]);
    }

    #[test]
    fn source562_outer_sanitize_failure_follows_all_component_map_appends() {
        let (graph, coords, _) = two_distinct_final_sanitize_failures(false);
        let mut labels = vec![77];
        let mut mapping = vec![vec![91]];
        let error = source562_run(
            &graph,
            &coords,
            true,
            true,
            &mut labels,
            &mut mapping,
            super::FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap_err();
        assert_eq!(error.component_index(), Some(0));
        assert_eq!(mapping.len(), 3);
        assert_eq!(mapping[0], vec![91]);
        assert_eq!(labels.len(), graph.atoms.len());
        assert!(matches!(
            error.failure,
            MoleculeFragmentsFailure::FinalSanitize { .. }
        ));
    }

    #[test]
    fn source562_full_copy_and_slow_delete_transport_actual_bookmarks_and_rings() {
        use std::collections::BTreeMap;
        let (single, single_coords, _, _) = source_for(&[ComponentShape::PlainDimer]);
        let ring = super::source_uninitialized_ring_info();
        let atom_marks = BTreeMap::from([(10, ids(&[1, 0, 1])), (11, vec![])]);
        let bond_marks = BTreeMap::from([(12, bond_ids(&[0, 0])), (13, vec![])]);
        let meta = super::FragmentSourceMetadataView {
            rings: Some(&ring),
            atom_bookmarks: Some(&atom_marks),
            bond_bookmarks: Some(&bond_marks),
        };
        let full = source562_run(
            &single,
            &single_coords,
            false,
            false,
            &mut vec![],
            &mut vec![],
            meta,
        )
        .unwrap();
        assert_eq!(full[0].source_metadata().rings.as_ref(), Some(&ring));
        assert_eq!(
            full[0].source_metadata().atom_bookmarks.as_ref().unwrap(),
            &BTreeMap::from([(10, ids(&[1, 0, 1]))])
        );
        assert_eq!(
            full[0].source_metadata().bond_bookmarks.as_ref().unwrap(),
            &BTreeMap::from([(12, bond_ids(&[0, 0]))])
        );
        assert!(full[0].coordinates().conformers_2d.is_empty());
        let (graph, coords, _, _) =
            source_for(&[ComponentShape::PlainDimer, ComponentShape::PlainDimer]);
        let atoms = BTreeMap::from([(7, ids(&[2, 3]))]);
        let bonds = BTreeMap::from([(8, bond_ids(&[1]))]);
        let fragments = source562_run(
            &graph,
            &coords,
            false,
            true,
            &mut vec![],
            &mut vec![],
            super::FragmentSourceMetadataView {
                rings: Some(&ring),
                atom_bookmarks: Some(&atoms),
                bond_bookmarks: Some(&bonds),
            },
        )
        .unwrap();
        assert_eq!(
            fragments[1]
                .source_metadata()
                .atom_bookmarks
                .as_ref()
                .unwrap()[&7],
            ids(&[0, 1])
        );
        assert_eq!(
            fragments[1]
                .source_metadata()
                .bond_bookmarks
                .as_ref()
                .unwrap()[&8],
            bond_ids(&[0])
        );
        assert!(
            !fragments[1]
                .source_metadata()
                .rings
                .as_ref()
                .unwrap()
                .is_initialized()
        );
    }

    #[test]
    fn source562_slow_delete_keeps_native_one_entry_stereo_vector_and_tag() {
        let atoms = (0..3)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double)
                .with_stereo(BondStereo::Cis)
                .with_stereo_atoms(AtomId::new(0), AtomId::new(2)),
        );
        let graph = TopologyBlock::try_from_parts(atoms, vec![bond], vec![], vec![]).unwrap();
        let result = source562_run(
            &graph,
            &coordinates(3),
            false,
            true,
            &mut vec![],
            &mut vec![],
            super::FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap();
        let bond = &result[0].topology().bonds[0];
        assert_eq!(bond.stereo(), BondStereo::Cis);
        assert_eq!(bond.stereo_atom_references(), &[AtomId::new(0)]);
    }
}

fn source_uninitialized_ring_info() -> crate::RingInfo {
    // RDKit❗✔️: RingInfo() {}
    // RDKit❗✔️: bool df_init{false};
    // Source constructor state, distinct from an unmodeled metadata capability.
    crate::RingInfo::from_persisted_components(
        false,
        crate::RingFindType::OtherOrUnknown,
        0,
        0,
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        None,
        Vec::new(),
        Vec::new(),
    )
    .expect("exact source uninitialized RingInfo constructor state is structurally valid")
}

fn fresh_subset_source_metadata(source: FragmentSourceMetadataView<'_>) -> FragmentSourceMetadata {
    // RDKit❗✔️: auto extracted_mol = std::make_unique<RWMol>();
    // RDKit❗✔️: dp_ringInfo = new RingInfo();
    // A fresh subset owns empty bookmarks and, after clearComputedProps=true,
    // an uninitialized ring value. Only explicitly modeled capabilities are
    // returned; None remains unmodeled rather than a native empty value.
    FragmentSourceMetadata {
        rings: source.rings.map(|_| source_uninitialized_ring_info()),
        atom_bookmarks: source.atom_bookmarks.map(|_| BTreeMap::new()),
        bond_bookmarks: source.bond_bookmarks.map(|_| BTreeMap::new()),
    }
}

fn clone_source_fragment_metadata(
    topology: &TopologyBlock,
    source: FragmentSourceMetadataView<'_>,
) -> Result<FragmentSourceMetadata, FullCopyComponentError> {
    // RDKit❗❌:     // Bookmarks should be copied as well:
    // RDKit❗❌:     for (auto abmI : other.d_atomBookmarks) {
    // RDKit❗❌:       for (const auto *aptr : abmI.second) {
    // RDKit❗❌:         setAtomBookmark(getAtomWithIdx(aptr->getIdx()), abmI.first);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     for (auto bbmI : other.d_bondBookmarks) {
    // RDKit❗❌:       for (const auto *bptr : bbmI.second) {
    // RDKit❗❌:         setBondBookmark(getBondWithIdx(bptr->getIdx()), bbmI.first);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // Native bookmark setters append once per pointer; empty source keys are
    // not copied. Typed IDs retain native duplicate encounter order.
    let atoms = source
        .atom_bookmarks
        .map(|marks| {
            let mut copied = BTreeMap::new();
            for (mark, references) in marks {
                for atom in references {
                    if atom.index() >= topology.atoms.len() {
                        return Err(FullCopyComponentError::TopologyEdit(
                            TopologyEditError::AtomOutOfRange {
                                atom: *atom,
                                atom_count: topology.atoms.len(),
                            },
                        ));
                    }
                    copied.entry(*mark).or_insert_with(Vec::new).push(*atom);
                }
            }
            Ok(copied)
        })
        .transpose()?;
    let bonds = source
        .bond_bookmarks
        .map(|marks| {
            let mut copied = BTreeMap::new();
            for (mark, references) in marks {
                for bond in references {
                    if bond.index() >= topology.bonds.len() {
                        return Err(FullCopyComponentError::TopologyEdit(
                            TopologyEditError::BondOutOfRange {
                                bond: *bond,
                                bond_count: topology.bonds.len(),
                            },
                        ));
                    }
                    copied.entry(*mark).or_insert_with(Vec::new).push(*bond);
                }
            }
            Ok(copied)
        })
        .transpose()?;
    Ok(FragmentSourceMetadata {
        rings: source.rings.cloned(),
        atom_bookmarks: atoms,
        bond_bookmarks: bonds,
    })
}

#[cfg(test)]
mod source566_shared_fragment_tests {
    use super::*;
    use cosmolkit_model::AtomSpec;
    use cosmolkit_types::Element;

    #[test]
    fn source566_shared_fragment_clone_retains_identical_owned_value() {
        let graph = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        let coordinates = CoordinateBlock::default();
        let mut labels = vec![77];
        let mut maps = vec![vec![91]];
        let fragments = get_shared_molecule_fragments_with_source_outputs(
            &graph,
            &FragmentCoordinateView::from_coordinate_block(&coordinates),
            &MoleculeProperties::default(),
            false,
            true,
            Some(&mut labels),
            Some(&mut maps),
            FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap();
        assert_eq!(fragments.len(), 1);
        assert_eq!(labels, vec![0]);
        assert_eq!(maps, vec![vec![91], vec![0]]);
        let shared = std::sync::Arc::clone(&fragments[0]);
        assert!(std::sync::Arc::ptr_eq(&shared, &fragments[0]));
        assert_eq!(std::sync::Arc::strong_count(&shared), 2);
        assert_eq!(shared.topology(), &graph);
        drop(fragments);
        assert_eq!(std::sync::Arc::strong_count(&shared), 1);
        assert_eq!(shared.component_atoms(), &[AtomId::new(0)]);
    }

    #[test]
    fn source566_shared_wrapper_preserves_inner_error_and_partial_outputs() {
        let graph = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(
                    AtomId::new(1),
                    AtomSpec::new(Element::O).with_explicit_hydrogens(3),
                ),
            ],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        let coordinates = CoordinateBlock::default();
        let mut labels = vec![77];
        let mut maps = vec![vec![91]];
        let error = get_shared_molecule_fragments_with_source_outputs(
            &graph,
            &FragmentCoordinateView::from_coordinate_block(&coordinates),
            &MoleculeProperties::default(),
            true,
            true,
            Some(&mut labels),
            Some(&mut maps),
            FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap_err();
        assert_eq!(error.component_index(), Some(1));
        assert_eq!(labels, vec![0, 1]);
        assert_eq!(maps, vec![vec![91], vec![0]]);
        assert!(matches!(
            error.failure,
            MoleculeFragmentsFailure::Build(OrderedFragmentBuildError::FastSubset {
                component_index: 1,
                ..
            })
        ));
    }
}

#[cfg(test)]
mod source570_owned_fragment_tests {
    use super::*;
    use cosmolkit_model::AtomSpec;
    use cosmolkit_types::Element;

    fn single_graph() -> TopologyBlock {
        TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }

    #[test]
    fn source570_owned_output_replaces_old_values_after_success_and_empty_success() {
        let graph = single_graph();
        let coordinates = CoordinateBlock::default();
        let view = FragmentCoordinateView::from_coordinate_block(&coordinates);
        let props = MoleculeProperties::default();
        let mut output = get_molecule_fragments(&graph, &coordinates, &props, false, true).unwrap();
        let mut labels = vec![77];
        let mut maps = vec![vec![91]];
        let count = assign_molecule_fragments_with_source_outputs(
            &TopologyBlock::default(),
            &view,
            &props,
            &mut output,
            false,
            true,
            Some(&mut labels),
            Some(&mut maps),
            FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap();
        assert_eq!(count, 0);
        assert!(output.is_empty());
        assert!(labels.is_empty());
        assert_eq!(maps, vec![vec![91]]);
        let count = assign_molecule_fragments_with_source_outputs(
            &graph,
            &view,
            &props,
            &mut output,
            false,
            true,
            Some(&mut labels),
            Some(&mut maps),
            FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap();
        assert_eq!(count, 1);
        assert_eq!(output.len(), 1);
        assert_eq!(output[0].topology(), &graph);
        assert_eq!(labels, vec![0]);
        assert_eq!(maps, vec![vec![91], vec![0]]);
    }

    #[test]
    fn source570_owned_output_failure_retains_prior_fragments_and_completed_maps() {
        let graph = single_graph();
        let coordinates = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        let mut output = get_molecule_fragments(&graph, &coordinates, &props, false, true).unwrap();
        let before = output.clone();
        let invalid = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(
                    AtomId::new(1),
                    AtomSpec::new(Element::O).with_explicit_hydrogens(3),
                ),
            ],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        let mut labels = vec![77];
        let mut maps = vec![vec![91]];
        let error = assign_molecule_fragments_with_source_outputs(
            &invalid,
            &FragmentCoordinateView::from_coordinate_block(&coordinates),
            &props,
            &mut output,
            true,
            true,
            Some(&mut labels),
            Some(&mut maps),
            FragmentSourceMetadataView::unmodeled(),
        )
        .unwrap_err();
        assert_eq!(error.component_index(), Some(1));
        assert_eq!(output, before);
        assert_eq!(labels, vec![0, 1]);
        assert_eq!(maps, vec![vec![91], vec![0]]);
    }
}
