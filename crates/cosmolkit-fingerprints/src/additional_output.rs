//! Detached storage for the optional RDKit fingerprint metadata outputs.
//!
//! Source attribution: pinned RDKit `GraphMol/Fingerprints/
//! FingerprintGenerator.{cpp,h}` and its `AdditionalOutput` wrapper source;
//! exact source paths, headers and notices are listed in
//! `THIRD_PARTY_NOTICES.md`.

use std::collections::BTreeMap;

use crate::FingerprintError;

/// Optional per-call fingerprint metadata.
///
/// This owner is private to the fingerprint domain. `None` represents a null
/// source pointer; `Some(empty)` represents an allocated empty source
/// container. It intentionally does not implement `Clone`, matching the
/// source aggregate's unique ownership.
#[derive(Debug, PartialEq, Eq)]
pub struct FingerprintAdditionalOutput {
    pub(crate) atom_counts: Option<Vec<u32>>,
    pub(crate) atom_to_bits: Option<Vec<Vec<u64>>>,
    pub(crate) bit_info_map: Option<BTreeMap<u64, Vec<(u32, u32)>>>,
    pub(crate) bit_paths: Option<BTreeMap<u64, Vec<Vec<i32>>>>,
    pub(crate) atoms_per_bit: Option<BTreeMap<u64, Vec<Vec<i32>>>>,
}

impl Default for FingerprintAdditionalOutput {
    fn default() -> Self {
        Self::new()
    }
}

impl FingerprintAdditionalOutput {
    /// Construct source-defined output with all five collection pointers absent.
    pub fn new() -> Self {
        // RDKit✔️✔️: atomToBitsType *atomToBits = nullptr;
        // RDKit✔️✔️: bitInfoMapType *bitInfoMap = nullptr;
        // RDKit✔️✔️: bitPathsType *bitPaths = nullptr;
        // RDKit✔️✔️: atomCountsType *atomCounts = nullptr;
        // RDKit✔️✔️: atomsPerBitType *atomsPerBit = nullptr;
        // One O(1) value initialization preserves all five null states; it
        // performs no container allocation, matching the source defaults.
        Self {
            atom_counts: None,
            atom_to_bits: None,
            bit_info_map: None,
            bit_paths: None,
            atoms_per_bit: None,
        }
    }
    /// Source `allocateAtomCounts`: replace this output with an allocated,
    /// empty count vector.
    pub fn allocate_atom_counts(&mut self) {
        // RDKit✔️🔝: void allocateAtomCounts() {
        // RDKit✔️🔝:   atomCountsHolder.reset(new atomCountsType);
        // RDKit✔️🔝:   atomCounts = atomCountsHolder.get();
        // RDKit✔️🔝: }
        // Rust stores the empty Vec inline in Option, preserving Some(empty)
        // while avoiding the source's separate heap allocation for the
        // container object; no backing buffer is allocated for an empty Vec.
        self.atom_counts = Some(Vec::new());
    }

    /// Source `allocateAtomToBits`: replace this output with an allocated,
    /// empty per-atom bit-vector table.
    pub fn allocate_atom_to_bits(&mut self) {
        // RDKit✔️🔝: void allocateAtomToBits() {
        // RDKit✔️🔝:   atomToBitsHolder.reset(new atomToBitsType);
        // RDKit✔️🔝:   atomToBits = atomToBitsHolder.get();
        // RDKit✔️🔝: }
        // The outer Vec is stored inline in Option; its empty state has no
        // backing allocation, so the source container-object allocation is
        // avoided without changing the observable allocation state.
        self.atom_to_bits = Some(Vec::new());
    }

    /// Source `allocateBitInfoMap`: replace this output with an allocated,
    /// empty bit-information map.
    pub fn allocate_bit_info_map(&mut self) {
        // RDKit✔️🔝: void allocateBitInfoMap() {
        // RDKit✔️🔝:   bitInfoMapHolder.reset(new bitInfoMapType);
        // RDKit✔️🔝:   bitInfoMap = bitInfoMapHolder.get();
        // RDKit✔️🔝: }
        // BTreeMap stores its empty root state inline. Keeping it in Option
        // removes only the separate source container-object allocation.
        self.bit_info_map = Some(BTreeMap::new());
    }

    /// Source `allocateBitPaths`: replace this output with an allocated,
    /// empty bit-path map.
    pub fn allocate_bit_paths(&mut self) {
        // RDKit✔️🔝: void allocateBitPaths() {
        // RDKit✔️🔝:   bitPathsHolder.reset(new bitPathsType);
        // RDKit✔️🔝:   bitPaths = bitPathsHolder.get();
        // RDKit✔️🔝: }
        // BTreeMap stores its empty root state inline. Keeping it in Option
        // removes only the separate source container-object allocation.
        self.bit_paths = Some(BTreeMap::new());
    }

    /// Source `allocateAtomsPerBit`: replace this output with an allocated,
    /// empty atoms-per-bit map.
    pub fn allocate_atoms_per_bit(&mut self) {
        // RDKit✔️🔝: void allocateAtomsPerBit() {
        // RDKit✔️🔝:   atomsPerBitHolder.reset(new atomsPerBitType);
        // RDKit✔️🔝:   atomsPerBit = atomsPerBitHolder.get();
        // RDKit✔️🔝: }
        // BTreeMap stores its empty root state inline. Keeping it in Option
        // removes only the separate source container-object allocation.
        self.atoms_per_bit = Some(BTreeMap::new());
    }

    /// Borrow the atom-count output. `None` preserves a null source pointer;
    /// an allocated but empty output is `Some(&[])`.
    #[must_use]
    pub fn atom_counts(&self) -> Option<&[u32]> {
        // BEGIN RDKIT CPP FUNCTION getAtomCountsHelper
        // RDKit✔️🔝: python::object getAtomCountsHelper(const AdditionalOutput &ao) {
        // RDKit✔️🔝:   if (!ao.atomCounts) {
        // RDKit✔️🔝:     return python::object();
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   python::list res;
        // RDKit✔️🔝:   for (const auto v : *ao.atomCounts) {
        // RDKit✔️🔝:     res.append(v);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   return python::tuple(res);
        // RDKit✔️🔝: }
        // END RDKIT CPP FUNCTION getAtomCountsHelper
        // The Rust owner returns the same optional values by borrow, avoiding
        // the wrapper's tuple copy here; language conversion remains at the
        // binding boundary.
        self.atom_counts.as_deref()
    }

    /// Borrow the per-atom bit output without copying its rows or bit order.
    #[must_use]
    pub fn atom_to_bits(&self) -> Option<&[Vec<u64>]> {
        // BEGIN RDKIT CPP FUNCTION getAtomToBitsHelper
        // RDKit✔️🔝: python::object getAtomToBitsHelper(const AdditionalOutput &ao) {
        // RDKit✔️🔝:   if (!ao.atomToBits) {
        // RDKit✔️🔝:     return python::object();
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   python::list res;
        // RDKit✔️🔝:   for (const auto &lst : *ao.atomToBits) {
        // RDKit✔️🔝:     python::list local;
        // RDKit✔️🔝:     for (const auto v : lst) {
        // RDKit✔️🔝:       local.append(v);
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     res.append(python::tuple(local));
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   return python::tuple(res);
        // RDKit✔️🔝: }
        // END RDKIT CPP FUNCTION getAtomToBitsHelper
        // The Rust owner borrows the original ordered rows instead of building
        // the wrapper's nested tuple copies; binding conversion stays separate.
        self.atom_to_bits.as_deref()
    }

    /// Borrow the bit-path map in source key and path order.
    #[must_use]
    pub fn bit_paths(&self) -> Option<&BTreeMap<u64, Vec<Vec<i32>>>> {
        // BEGIN RDKIT CPP FUNCTION getBitPathsHelper
        // RDKit✔️🔝: python::object getBitPathsHelper(const AdditionalOutput &ao) {
        // RDKit✔️🔝:   if (!ao.bitPaths) {
        // RDKit✔️🔝:     return python::object();
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   python::dict res;
        // RDKit✔️🔝:   for (const auto &pr : *ao.bitPaths) {
        // RDKit✔️🔝:     python::list local;
        // RDKit✔️🔝:     for (const auto &lst : pr.second) {
        // RDKit✔️🔝:       python::list inner;
        // RDKit✔️🔝:       for (const auto v : lst) {
        // RDKit✔️🔝:         inner.append(v);
        // RDKit✔️🔝:       }
        // RDKit✔️🔝:       local.append(python::tuple(inner));
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     res[pr.first] = python::tuple(local);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   return res;
        // RDKit✔️🔝: }
        // END RDKIT CPP FUNCTION getBitPathsHelper
        // Borrowing retains the same key/value contents without materializing
        // the wrapper's nested Python lists and dictionary; binding conversion
        // remains at its own boundary.
        self.bit_paths.as_ref()
    }

    /// Borrow the Morgan bit provenance map without copying centers or order.
    #[must_use]
    pub fn bit_info_map(&self) -> Option<&BTreeMap<u64, Vec<(u32, u32)>>> {
        // BEGIN RDKIT CPP FUNCTION getBitInfoMapHelper
        // RDKit✔️🔝: python::object getBitInfoMapHelper(const AdditionalOutput &ao) {
        // RDKit✔️🔝:   if (!ao.bitInfoMap) {
        // RDKit✔️🔝:     return python::object();
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   python::dict res;
        // RDKit✔️🔝:   for (const auto &pr : *ao.bitInfoMap) {
        // RDKit✔️🔝:     python::list local;
        // RDKit✔️🔝:     for (const auto &v : pr.second) {
        // RDKit✔️🔝:       python::tuple inner = python::make_tuple(v.first, v.second);
        // RDKit✔️🔝:       local.append(inner);
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     res[pr.first] = python::tuple(local);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   return res;
        // RDKit✔️🔝: }
        // END RDKIT CPP FUNCTION getBitInfoMapHelper
        // This borrow preserves each stored ordered center pair; unlike the
        // Python helper it creates no map, tuple, or center copies here.
        self.bit_info_map.as_ref()
    }

    /// Borrow the atoms-per-bit map without copying its ordered atom rows.
    #[must_use]
    pub fn atoms_per_bit(&self) -> Option<&BTreeMap<u64, Vec<Vec<i32>>>> {
        // BEGIN RDKIT CPP FUNCTION getAtomsPerBitHelper
        // RDKit✔️🔝: python::object getAtomsPerBitHelper(const AdditionalOutput &ao) {
        // RDKit✔️🔝:   if (!ao.atomsPerBit) {
        // RDKit✔️🔝:     return python::object();
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   python::dict res;
        // RDKit✔️🔝:   for (const auto &pr : *ao.atomsPerBit) {
        // RDKit✔️🔝:     python::list local;
        // RDKit✔️🔝:     for (const auto &lst : pr.second) {
        // RDKit✔️🔝:       python::list inner;
        // RDKit✔️🔝:       for (const auto v : lst) {
        // RDKit✔️🔝:         inner.append(v);
        // RDKit✔️🔝:       }
        // RDKit✔️🔝:       local.append(python::tuple(inner));
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     res[pr.first] = python::tuple(local);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   return res;
        // RDKit✔️🔝: }
        // END RDKIT CPP FUNCTION getAtomsPerBitHelper
        // Borrowing preserves the same map and list contents without producing
        // wrapper tuples, nested lists, or a temporary dictionary here.
        self.atoms_per_bit.as_ref()
    }

    /// Prepare the temporary metadata container used by the source's
    /// count-simulation caller path and reinitialize this original output.
    /// The caller invokes this only when its FingerprintAdditionalOutput is present and
    /// count simulation is enabled; the returned temporary is reinitialized
    /// later by the generic fingerprint helper.
    pub(crate) fn setup_count_simulation_output(&mut self, num_atoms: usize) -> Self {
        // RDKit✔️✔️: void setupTempAdditionalOutput(RDKit::FingerprintFuncArguments &args,
        // RDKit✔️✔️:                                AdditionalOutput &countSimulationOutput,
        // RDKit✔️✔️:                                size_t numAtoms) {
        // RDKit✔️✔️:   if (args.additionalOutput->atomToBits) {
        // RDKit✔️✔️:     countSimulationOutput.allocateAtomToBits();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (args.additionalOutput->atomCounts) {
        // RDKit✔️✔️:     countSimulationOutput.allocateAtomCounts();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (args.additionalOutput->bitInfoMap) {
        // RDKit✔️✔️:     countSimulationOutput.allocateBitInfoMap();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (args.additionalOutput->bitPaths) {
        // RDKit✔️✔️:     countSimulationOutput.allocateBitPaths();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (args.additionalOutput->atomsPerBit) {
        // RDKit✔️✔️:     countSimulationOutput.allocateAtomsPerBit();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   reinitAdditionalOutput(*args.additionalOutput, numAtoms);
        // RDKit✔️✔️: }
        // The source pointer is the present `self`; preserve its five independent
        // allocation tests and allocate the returned temporary in source order.
        let mut count_simulation_output = Self::default();
        if self.atom_to_bits.is_some() {
            count_simulation_output.allocate_atom_to_bits();
        }
        if self.atom_counts.is_some() {
            count_simulation_output.allocate_atom_counts();
        }
        if self.bit_info_map.is_some() {
            count_simulation_output.allocate_bit_info_map();
        }
        if self.bit_paths.is_some() {
            count_simulation_output.allocate_bit_paths();
        }
        if self.atoms_per_bit.is_some() {
            count_simulation_output.allocate_atoms_per_bit();
        }
        self.reinitialize(num_atoms);

        // Local complexity review: five presence checks and empty-container
        // allocations are O(1) per selected field. Reinitializing the original
        // remains O(A) for per-atom vectors and O(M) for the two maps, exactly
        // the source helper's single call; the generic caller later performs
        // its separately sourced reinitialization on the returned temporary.
        count_simulation_output
    }

    /// Copy one count-simulation base bit's provenance to its threshold bit.
    /// The caller supplies the temporary output as `self` and the original
    /// output as `new_output`, matching RDKit's argument order.
    pub(crate) fn duplicate_bit_to(
        &self,
        new_output: &mut Self,
        orig_bit_id: u64,
        new_bit_id: u64,
    ) -> Result<(), FingerprintError> {
        // RDKit✔️✔️: template <typename OutputType>
        // RDKit✔️✔️: void duplicateAdditionalOutputBit(AdditionalOutput &oldAO,
        // RDKit✔️✔️:                                   AdditionalOutput &newAO, OutputType origBitId,
        // RDKit✔️✔️:                                   OutputType newBitId) {
        // RDKit✔️✔️:   PRECONDITION(!((oldAO.bitInfoMap != nullptr) ^ (newAO.bitInfoMap != nullptr)),
        // RDKit✔️✔️:                "bitInfoMap not allocated");
        if self.bit_info_map.is_some() ^ new_output.bit_info_map.is_some() {
            return Err(FingerprintError::PreconditionViolation {
                what: "bitInfoMap not allocated",
            });
        }

        // RDKit✔️✔️:   PRECONDITION(!((oldAO.atomToBits != nullptr) ^ (newAO.atomToBits != nullptr)),
        // RDKit✔️✔️:                "atomToBits not allocated");
        if self.atom_to_bits.is_some() ^ new_output.atom_to_bits.is_some() {
            return Err(FingerprintError::PreconditionViolation {
                what: "atomToBits not allocated",
            });
        }

        // RDKit✔️✔️:   PRECONDITION(!((oldAO.bitPaths != nullptr) ^ (newAO.bitPaths != nullptr)),
        // RDKit✔️✔️:                "bitPaths not allocated");
        if self.bit_paths.is_some() ^ new_output.bit_paths.is_some() {
            return Err(FingerprintError::PreconditionViolation {
                what: "bitPaths not allocated",
            });
        }

        // RDKit✔️✔️:   // we don't need to do anything with atomCounts

        // RDKit✔️✔️:   if (oldAO.atomToBits) {
        // RDKit✔️✔️:     if (newAO.atomToBits->empty()) {
        // RDKit✔️✔️:       newAO.atomToBits->resize(oldAO.atomToBits->size());
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     for (unsigned int i = 0; i < oldAO.atomToBits->size(); ++i) {
        // RDKit✔️✔️:       const auto &nv = oldAO.atomToBits->at(i);
        // RDKit✔️✔️:       if (std::find(nv.begin(), nv.end(), origBitId) != nv.end()) {
        // RDKit✔️✔️:         newAO.atomToBits->at(i).push_back(newBitId);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        if let Some(old_atom_to_bits) = self.atom_to_bits.as_ref() {
            let new_atom_to_bits = new_output
                .atom_to_bits
                .as_mut()
                .expect("the source allocation check above keeps atomToBits present");
            if new_atom_to_bits.is_empty() {
                new_atom_to_bits.resize_with(old_atom_to_bits.len(), Vec::new);
            }
            for (atom_index, old_bits) in old_atom_to_bits.iter().enumerate() {
                if old_bits.contains(&orig_bit_id) {
                    new_atom_to_bits
                        .get_mut(atom_index)
                        .expect("the source vector::at row access requires a matching atom count")
                        .push(new_bit_id);
                }
            }
        }

        // RDKit✔️✔️:   if (oldAO.bitInfoMap) {
        // RDKit✔️✔️:     const auto v = oldAO.bitInfoMap->find(origBitId);
        // RDKit✔️✔️:     if (v != oldAO.bitInfoMap->end()) {
        // RDKit✔️✔️:       (*newAO.bitInfoMap)[newBitId] = v->second;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        if let Some(old_bit_info_map) = self.bit_info_map.as_ref() {
            if let Some(value) = old_bit_info_map.get(&orig_bit_id) {
                new_output
                    .bit_info_map
                    .as_mut()
                    .expect("the source allocation check above keeps bitInfoMap present")
                    .insert(new_bit_id, value.clone());
            }
        }

        // RDKit✔️✔️:   if (oldAO.bitPaths) {
        // RDKit✔️✔️:     const auto v = oldAO.bitPaths->find(origBitId);
        // RDKit✔️✔️:     if (v != oldAO.bitPaths->end()) {
        // RDKit✔️✔️:       (*newAO.bitPaths)[newBitId] = v->second;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        if let Some(old_bit_paths) = self.bit_paths.as_ref() {
            if let Some(value) = old_bit_paths.get(&orig_bit_id) {
                new_output
                    .bit_paths
                    .as_mut()
                    .expect("the source allocation check above keeps bitPaths present")
                    .insert(new_bit_id, value.clone());
            }
        }

        // RDKit✔️✔️:   if (oldAO.atomsPerBit) {
        // RDKit✔️✔️:     const auto v = oldAO.atomsPerBit->find(origBitId);
        // RDKit✔️✔️:     if (v != oldAO.atomsPerBit->end()) {
        // RDKit✔️✔️:       (*newAO.atomsPerBit)[newBitId] = v->second;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        if let Some(old_atoms_per_bit) = self.atoms_per_bit.as_ref() {
            if let Some(value) = old_atoms_per_bit.get(&orig_bit_id) {
                // The pinned caller mirrors this allocation before dispatch;
                // source has no parity check here and null-dereferences if the
                // original entry exists without its destination container.
                new_output
                    .atoms_per_bit
                    .as_mut()
                    .expect("the source caller must mirror atomsPerBit allocation")
                    .insert(new_bit_id, value.clone());
            }
        }

        // Behavior follows source field order and map assignment semantics.
        // Local complexity review: per call, scan each A atom row once and
        // each row's memberships only until a match, then append at most once
        // for that call; perform one O(log M) lookup/insert per provenance map
        // that contains the source key, cloning only that value. No global
        // atomToBits deduplication or whole-output clone is introduced.
        Ok(())
    }

    /// Reinitialize only the fields handled by RDKit's generic output helper.
    /// `atoms_per_bit` is deliberately retained because that source helper has
    /// no branch for it.
    pub(crate) fn reinitialize(&mut self, num_atoms: usize) {
        // RDKit✔️✔️: void reinitAdditionalOutput(AdditionalOutput &ao, size_t numAtoms) {
        // RDKit✔️✔️:   if (ao.atomCounts) {
        // RDKit✔️✔️:     ao.atomCounts->resize(numAtoms);
        // RDKit✔️✔️:     std::fill(ao.atomCounts->begin(), ao.atomCounts->end(), 0);
        // RDKit✔️✔️:   }
        if let Some(atom_counts) = &mut self.atom_counts {
            atom_counts.resize(num_atoms, 0);
            atom_counts.fill(0);
        }

        // RDKit✔️✔️:   if (ao.atomToBits) {
        // RDKit✔️✔️:     ao.atomToBits->resize(numAtoms);
        // RDKit✔️✔️:     std::fill(ao.atomToBits->begin(), ao.atomToBits->end(),
        // RDKit✔️✔️:               std::vector<std::uint64_t>());
        // RDKit✔️✔️:   }
        if let Some(atom_to_bits) = &mut self.atom_to_bits {
            atom_to_bits.resize_with(num_atoms, Vec::new);
            for bit_list in atom_to_bits {
                // Assignment drops any prior per-atom buffer, matching fill
                // with an empty source vector instead of retaining capacity.
                *bit_list = Vec::new();
            }
        }

        // RDKit✔️✔️:   if (ao.bitInfoMap) {
        // RDKit✔️✔️:     ao.bitInfoMap->clear();
        // RDKit✔️✔️:   }
        if let Some(bit_info_map) = &mut self.bit_info_map {
            bit_info_map.clear();
        }

        // RDKit✔️✔️:   if (ao.bitPaths) {
        // RDKit✔️✔️:     ao.bitPaths->clear();
        // RDKit✔️✔️:   }
        if let Some(bit_paths) = &mut self.bit_paths {
            bit_paths.clear();
        }

        // There is intentionally no atomsPerBit branch in the source helper;
        // preserve its value and allocation state exactly.
        // Local complexity review: the two per-atom vectors take O(A), each
        // ordered-map clear takes O(M), and absent fields are O(1) skips.
        // These costs match the source's vector resize/fill and map clear.
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use crate::FingerprintError;

    use super::FingerprintAdditionalOutput;

    #[test]
    fn fingerprint_additional_output_canonical_constructors_preserve_absent_state() {
        let new: fn() -> FingerprintAdditionalOutput = FingerprintAdditionalOutput::new;
        let default: fn() -> FingerprintAdditionalOutput = FingerprintAdditionalOutput::default;
        let output = new();
        assert_eq!(output, default());
        assert_eq!(output.atom_counts(), None);
        assert_eq!(output.atom_to_bits(), None);
        assert_eq!(output.bit_info_map(), None);
        assert_eq!(output.bit_paths(), None);
        assert_eq!(output.atoms_per_bit(), None);
    }

    #[test]
    fn fingerprint_public_additional_output_allocation_masks_and_borrowed_getters() {
        let mut masks_checked = 0;
        for mask in 0_u8..32 {
            let expected = [
                mask & 0b00001 != 0,
                mask & 0b00010 != 0,
                mask & 0b00100 != 0,
                mask & 0b01000 != 0,
                mask & 0b10000 != 0,
            ];
            let mut output = FingerprintAdditionalOutput::default();
            if expected[0] {
                output.allocate_atom_counts();
            }
            if expected[1] {
                output.allocate_atom_to_bits();
            }
            if expected[2] {
                output.allocate_bit_info_map();
            }
            if expected[3] {
                output.allocate_bit_paths();
            }
            if expected[4] {
                output.allocate_atoms_per_bit();
            }

            assert_eq!(
                output.atom_counts().map(|values| values.len()),
                expected[0].then_some(0),
                "atomCounts None/Some(empty), mask {mask:#07b}"
            );
            assert_eq!(
                output.atom_to_bits().map(|rows| rows.len()),
                expected[1].then_some(0),
                "atomToBits None/Some(empty), mask {mask:#07b}"
            );
            assert_eq!(
                output.bit_info_map().map(BTreeMap::len),
                expected[2].then_some(0),
                "bitInfoMap None/Some(empty), mask {mask:#07b}"
            );
            assert_eq!(
                output.bit_paths().map(BTreeMap::len),
                expected[3].then_some(0),
                "bitPaths None/Some(empty), mask {mask:#07b}"
            );
            assert_eq!(
                output.atoms_per_bit().map(BTreeMap::len),
                expected[4].then_some(0),
                "atomsPerBit None/Some(empty), mask {mask:#07b}"
            );
            masks_checked += 1;
        }
        assert_eq!(masks_checked, 32, "all five-bit allocation masks");

        let mut replaced = FingerprintAdditionalOutput {
            atom_counts: Some(vec![9]),
            atom_to_bits: Some(vec![vec![8, 8]]),
            bit_info_map: Some(BTreeMap::from([(7, vec![(1, 2)])])),
            bit_paths: Some(BTreeMap::from([(6, vec![vec![3, 4]])])),
            atoms_per_bit: Some(BTreeMap::from([(5, vec![vec![0, 2]])])),
        };
        replaced.allocate_atom_counts();
        replaced.allocate_atom_to_bits();
        replaced.allocate_bit_info_map();
        replaced.allocate_bit_paths();
        replaced.allocate_atoms_per_bit();
        assert_eq!(replaced.atom_counts(), Some(&[][..]));
        assert_eq!(replaced.atom_to_bits(), Some(&[][..]));
        assert!(replaced.bit_info_map().is_some_and(BTreeMap::is_empty));
        assert!(replaced.bit_paths().is_some_and(BTreeMap::is_empty));
        assert!(replaced.atoms_per_bit().is_some_and(BTreeMap::is_empty));

        let expected_atom_counts = vec![4, 1, 4];
        let expected_atom_to_bits = vec![vec![9, 9, 3], vec![2], vec![]];
        let expected_bit_info =
            BTreeMap::from([(8, vec![(1, 0), (1, 2), (1, 0)]), (3, vec![(0, 0)])]);
        let expected_bit_paths =
            BTreeMap::from([(11, vec![vec![3, 2, 3], vec![8]]), (2, vec![vec![7, 7]])]);
        let expected_atoms_per_bit =
            BTreeMap::from([(12, vec![vec![0, 2], vec![0, 2]]), (1, vec![vec![4, 5]])]);
        let output = FingerprintAdditionalOutput {
            atom_counts: Some(expected_atom_counts.clone()),
            atom_to_bits: Some(expected_atom_to_bits.clone()),
            bit_info_map: Some(expected_bit_info.clone()),
            bit_paths: Some(expected_bit_paths.clone()),
            atoms_per_bit: Some(expected_atoms_per_bit.clone()),
        };

        assert_eq!(output.atom_counts(), Some(expected_atom_counts.as_slice()));
        assert_eq!(
            output.atom_to_bits(),
            Some(expected_atom_to_bits.as_slice())
        );
        assert_eq!(output.bit_info_map(), Some(&expected_bit_info));
        assert_eq!(output.bit_paths(), Some(&expected_bit_paths));
        assert_eq!(output.atoms_per_bit(), Some(&expected_atoms_per_bit));
        assert!(std::ptr::eq(
            output.atom_counts().unwrap(),
            output.atom_counts.as_ref().unwrap().as_slice()
        ));
        assert!(std::ptr::eq(
            output.atom_to_bits().unwrap(),
            output.atom_to_bits.as_ref().unwrap().as_slice()
        ));
        assert!(std::ptr::eq(
            output.bit_info_map().unwrap(),
            output.bit_info_map.as_ref().unwrap()
        ));
        assert!(std::ptr::eq(
            output.bit_paths().unwrap(),
            output.bit_paths.as_ref().unwrap()
        ));
        assert!(std::ptr::eq(
            output.atoms_per_bit().unwrap(),
            output.atoms_per_bit.as_ref().unwrap()
        ));
        assert_eq!(output.atom_counts(), Some(expected_atom_counts.as_slice()));
        assert_eq!(
            output.atom_to_bits(),
            Some(expected_atom_to_bits.as_slice())
        );
        assert_eq!(output.bit_info_map(), Some(&expected_bit_info));
        assert_eq!(output.bit_paths(), Some(&expected_bit_paths));
        assert_eq!(output.atoms_per_bit(), Some(&expected_atoms_per_bit));
    }

    #[test]
    fn fingerprint_morgan_g03_allocation_masks() {
        const EXPECTED_PRESENCE: [[bool; 5]; 32] = [
            [false, false, false, false, false],
            [true, false, false, false, false],
            [false, true, false, false, false],
            [true, true, false, false, false],
            [false, false, true, false, false],
            [true, false, true, false, false],
            [false, true, true, false, false],
            [true, true, true, false, false],
            [false, false, false, true, false],
            [true, false, false, true, false],
            [false, true, false, true, false],
            [true, true, false, true, false],
            [false, false, true, true, false],
            [true, false, true, true, false],
            [false, true, true, true, false],
            [true, true, true, true, false],
            [false, false, false, false, true],
            [true, false, false, false, true],
            [false, true, false, false, true],
            [true, true, false, false, true],
            [false, false, true, false, true],
            [true, false, true, false, true],
            [false, true, true, false, true],
            [true, true, true, false, true],
            [false, false, false, true, true],
            [true, false, false, true, true],
            [false, true, false, true, true],
            [true, true, false, true, true],
            [false, false, true, true, true],
            [true, false, true, true, true],
            [false, true, true, true, true],
            [true, true, true, true, true],
        ];

        for (mask, expected) in EXPECTED_PRESENCE.iter().enumerate() {
            let mut output = FingerprintAdditionalOutput::default();
            if mask & 0b00001 != 0 {
                output.allocate_atom_counts();
            }
            if mask & 0b00010 != 0 {
                output.allocate_atom_to_bits();
            }
            if mask & 0b00100 != 0 {
                output.allocate_bit_info_map();
            }
            if mask & 0b01000 != 0 {
                output.allocate_bit_paths();
            }
            if mask & 0b10000 != 0 {
                output.allocate_atoms_per_bit();
            }

            assert_eq!(
                output.atom_counts,
                expected[0].then(Vec::new),
                "atomCounts, allocation mask {mask:#07b}"
            );
            assert_eq!(
                output.atom_to_bits,
                expected[1].then(Vec::new),
                "atomToBits, allocation mask {mask:#07b}"
            );
            assert_eq!(
                output.bit_info_map,
                expected[2].then(BTreeMap::new),
                "bitInfoMap, allocation mask {mask:#07b}"
            );
            assert_eq!(
                output.bit_paths,
                expected[3].then(BTreeMap::new),
                "bitPaths, allocation mask {mask:#07b}"
            );
            assert_eq!(
                output.atoms_per_bit,
                expected[4].then(BTreeMap::new),
                "atomsPerBit, allocation mask {mask:#07b}"
            );
        }

        let mut populated = FingerprintAdditionalOutput {
            atom_counts: Some(vec![7]),
            atom_to_bits: Some(vec![vec![11]]),
            bit_info_map: Some(BTreeMap::from([(13, vec![(2, 1)])])),
            bit_paths: Some(BTreeMap::from([(17, vec![vec![3, 4]])])),
            atoms_per_bit: Some(BTreeMap::from([(19, vec![vec![5, 6]])])),
        };
        populated.allocate_atom_counts();
        populated.allocate_atom_to_bits();
        populated.allocate_bit_info_map();
        populated.allocate_bit_paths();
        populated.allocate_atoms_per_bit();
        assert_eq!(
            populated,
            FingerprintAdditionalOutput {
                atom_counts: Some(Vec::new()),
                atom_to_bits: Some(Vec::new()),
                bit_info_map: Some(BTreeMap::new()),
                bit_paths: Some(BTreeMap::new()),
                atoms_per_bit: Some(BTreeMap::new()),
            }
        );
    }

    #[test]
    fn fingerprint_morgan_g04_reinitialize_masks_fresh_and_populated() {
        const EXPECTED_PRESENCE: [[bool; 5]; 32] = [
            [false, false, false, false, false],
            [true, false, false, false, false],
            [false, true, false, false, false],
            [true, true, false, false, false],
            [false, false, true, false, false],
            [true, false, true, false, false],
            [false, true, true, false, false],
            [true, true, true, false, false],
            [false, false, false, true, false],
            [true, false, false, true, false],
            [false, true, false, true, false],
            [true, true, false, true, false],
            [false, false, true, true, false],
            [true, false, true, true, false],
            [false, true, true, true, false],
            [true, true, true, true, false],
            [false, false, false, false, true],
            [true, false, false, false, true],
            [false, true, false, false, true],
            [true, true, false, false, true],
            [false, false, true, false, true],
            [true, false, true, false, true],
            [false, true, true, false, true],
            [true, true, true, false, true],
            [false, false, false, true, true],
            [true, false, false, true, true],
            [false, true, false, true, true],
            [true, true, false, true, true],
            [false, false, true, true, true],
            [true, false, true, true, true],
            [false, true, true, true, true],
            [true, true, true, true, true],
        ];

        for (mask, expected) in EXPECTED_PRESENCE.iter().enumerate() {
            let mut fresh = FingerprintAdditionalOutput::default();
            if mask & 0b00001 != 0 {
                fresh.allocate_atom_counts();
            }
            if mask & 0b00010 != 0 {
                fresh.allocate_atom_to_bits();
            }
            if mask & 0b00100 != 0 {
                fresh.allocate_bit_info_map();
            }
            if mask & 0b01000 != 0 {
                fresh.allocate_bit_paths();
            }
            if mask & 0b10000 != 0 {
                fresh.allocate_atoms_per_bit();
            }

            fresh.reinitialize(3);
            assert_eq!(
                fresh.atom_counts,
                expected[0].then(|| vec![0, 0, 0]),
                "fresh atomCounts, allocation mask {mask:#07b}"
            );
            assert_eq!(
                fresh.atom_to_bits,
                expected[1].then(|| vec![vec![], vec![], vec![]]),
                "fresh atomToBits, allocation mask {mask:#07b}"
            );
            assert_eq!(
                fresh.bit_info_map,
                expected[2].then(BTreeMap::new),
                "fresh bitInfoMap, allocation mask {mask:#07b}"
            );
            assert_eq!(
                fresh.bit_paths,
                expected[3].then(BTreeMap::new),
                "fresh bitPaths, allocation mask {mask:#07b}"
            );
            assert_eq!(
                fresh.atoms_per_bit,
                expected[4].then(BTreeMap::new),
                "fresh atomsPerBit, allocation mask {mask:#07b}"
            );

            let mut populated = FingerprintAdditionalOutput::default();
            if expected[0] {
                populated.allocate_atom_counts();
                populated.atom_counts = Some(vec![7, 8, 9, 10]);
            }
            if expected[1] {
                populated.allocate_atom_to_bits();
                populated.atom_to_bits = Some(vec![vec![11, 12], vec![13], vec![14, 15], vec![16]]);
            }
            if expected[2] {
                populated.allocate_bit_info_map();
                populated.bit_info_map = Some(BTreeMap::from([(17, vec![(2, 1)])]));
            }
            if expected[3] {
                populated.allocate_bit_paths();
                populated.bit_paths = Some(BTreeMap::from([(19, vec![vec![3, 4]])]));
            }
            if expected[4] {
                populated.allocate_atoms_per_bit();
                populated.atoms_per_bit = Some(BTreeMap::from([(23, vec![vec![5, 6]])]));
            }

            populated.reinitialize(3);
            assert_eq!(
                populated.atom_counts,
                expected[0].then(|| vec![0, 0, 0]),
                "populated atomCounts, allocation mask {mask:#07b}"
            );
            assert_eq!(
                populated.atom_to_bits,
                expected[1].then(|| vec![vec![], vec![], vec![]]),
                "populated atomToBits, allocation mask {mask:#07b}"
            );
            assert_eq!(
                populated.bit_info_map,
                expected[2].then(BTreeMap::new),
                "populated bitInfoMap, allocation mask {mask:#07b}"
            );
            assert_eq!(
                populated.bit_paths,
                expected[3].then(BTreeMap::new),
                "populated bitPaths, allocation mask {mask:#07b}"
            );
            assert_eq!(
                populated.atoms_per_bit,
                expected[4].then(|| BTreeMap::from([(23, vec![vec![5, 6]])])),
                "populated atomsPerBit, allocation mask {mask:#07b}"
            );
        }
    }

    #[test]
    fn fingerprint_morgan_g05_count_simulation_output_masks_and_order() {
        const EXPECTED_PRESENCE: [[bool; 5]; 32] = [
            [false, false, false, false, false],
            [true, false, false, false, false],
            [false, true, false, false, false],
            [true, true, false, false, false],
            [false, false, true, false, false],
            [true, false, true, false, false],
            [false, true, true, false, false],
            [true, true, true, false, false],
            [false, false, false, true, false],
            [true, false, false, true, false],
            [false, true, false, true, false],
            [true, true, false, true, false],
            [false, false, true, true, false],
            [true, false, true, true, false],
            [false, true, true, true, false],
            [true, true, true, true, false],
            [false, false, false, false, true],
            [true, false, false, false, true],
            [false, true, false, false, true],
            [true, true, false, false, true],
            [false, false, true, false, true],
            [true, false, true, false, true],
            [false, true, true, false, true],
            [true, true, true, false, true],
            [false, false, false, true, true],
            [true, false, false, true, true],
            [false, true, false, true, true],
            [true, true, false, true, true],
            [false, false, true, true, true],
            [true, false, true, true, true],
            [false, true, true, true, true],
            [true, true, true, true, true],
        ];
        const COUNT_SIMULATION: [bool; 2] = [false, true];

        let mut cases = 0;
        let mut temporary_setups = 0;
        let mut reinitializations = 0;
        for count_simulation in COUNT_SIMULATION {
            for (mask, expected) in EXPECTED_PRESENCE.iter().enumerate() {
                let mut original = FingerprintAdditionalOutput::default();
                if expected[0] {
                    original.allocate_atom_counts();
                    original.atom_counts = Some(vec![7, 8, 9, 10]);
                }
                if expected[1] {
                    original.allocate_atom_to_bits();
                    original.atom_to_bits =
                        Some(vec![vec![11, 12], vec![13], vec![14, 15], vec![16]]);
                }
                if expected[2] {
                    original.allocate_bit_info_map();
                    original.bit_info_map = Some(BTreeMap::from([(17, vec![(2, 1)])]));
                }
                if expected[3] {
                    original.allocate_bit_paths();
                    original.bit_paths = Some(BTreeMap::from([(19, vec![vec![3, 4]])]));
                }
                if expected[4] {
                    original.allocate_atoms_per_bit();
                    original.atoms_per_bit = Some(BTreeMap::from([(23, vec![vec![5, 6]])]));
                }

                let mut temporary = if count_simulation {
                    temporary_setups += 1;
                    reinitializations += 1;
                    Some(original.setup_count_simulation_output(3))
                } else {
                    reinitializations += 1;
                    original.reinitialize(3);
                    None
                };

                // The original output is reset during setup, before the
                // temporary output is handed to the generic helper.
                assert_eq!(
                    original.atom_counts,
                    expected[0].then(|| vec![0, 0, 0]),
                    "original atomCounts, countSimulation={count_simulation}, mask {mask:#07b}"
                );
                assert_eq!(
                    original.atom_to_bits,
                    expected[1].then(|| vec![vec![], vec![], vec![]]),
                    "original atomToBits, countSimulation={count_simulation}, mask {mask:#07b}"
                );
                assert_eq!(
                    original.bit_info_map,
                    expected[2].then(BTreeMap::new),
                    "original bitInfoMap, countSimulation={count_simulation}, mask {mask:#07b}"
                );
                assert_eq!(
                    original.bit_paths,
                    expected[3].then(BTreeMap::new),
                    "original bitPaths, countSimulation={count_simulation}, mask {mask:#07b}"
                );
                assert_eq!(
                    original.atoms_per_bit,
                    expected[4].then(|| BTreeMap::from([(23, vec![vec![5, 6]])])),
                    "original atomsPerBit, countSimulation={count_simulation}, mask {mask:#07b}"
                );
                assert_eq!(temporary.is_some(), count_simulation);

                if let Some(output) = temporary.as_mut() {
                    // setupTempAdditionalOutput only mirrors allocation;
                    // getFingerprintHelper sizes this output afterward.
                    assert_eq!(
                        output.atom_counts,
                        expected[0].then(Vec::new),
                        "allocated temp atomCounts, mask {mask:#07b}"
                    );
                    assert_eq!(
                        output.atom_to_bits,
                        expected[1].then(Vec::new),
                        "allocated temp atomToBits, mask {mask:#07b}"
                    );
                    assert_eq!(
                        output.bit_info_map,
                        expected[2].then(BTreeMap::new),
                        "allocated temp bitInfoMap, mask {mask:#07b}"
                    );
                    assert_eq!(
                        output.bit_paths,
                        expected[3].then(BTreeMap::new),
                        "allocated temp bitPaths, mask {mask:#07b}"
                    );
                    assert_eq!(
                        output.atoms_per_bit,
                        expected[4].then(BTreeMap::new),
                        "allocated temp atomsPerBit, mask {mask:#07b}"
                    );

                    reinitializations += 1;
                    output.reinitialize(3);
                    assert_eq!(
                        output.atom_counts,
                        expected[0].then(|| vec![0, 0, 0]),
                        "sized temp atomCounts, mask {mask:#07b}"
                    );
                    assert_eq!(
                        output.atom_to_bits,
                        expected[1].then(|| vec![vec![], vec![], vec![]]),
                        "sized temp atomToBits, mask {mask:#07b}"
                    );
                    assert_eq!(
                        output.bit_info_map,
                        expected[2].then(BTreeMap::new),
                        "sized temp bitInfoMap, mask {mask:#07b}"
                    );
                    assert_eq!(
                        output.bit_paths,
                        expected[3].then(BTreeMap::new),
                        "sized temp bitPaths, mask {mask:#07b}"
                    );
                    assert_eq!(
                        output.atoms_per_bit,
                        expected[4].then(BTreeMap::new),
                        "sized temp atomsPerBit, mask {mask:#07b}"
                    );
                }
                cases += 1;
            }
        }

        assert_eq!(cases, 64);
        assert_eq!(temporary_setups, 32);
        assert_eq!(reinitializations, 96);
    }

    #[test]
    fn fingerprint_morgan_g06_duplicate_output_masks_collisions_and_centers() {
        const EXPECTED_PRESENCE: [[bool; 5]; 32] = [
            [false, false, false, false, false],
            [true, false, false, false, false],
            [false, true, false, false, false],
            [true, true, false, false, false],
            [false, false, true, false, false],
            [true, false, true, false, false],
            [false, true, true, false, false],
            [true, true, true, false, false],
            [false, false, false, true, false],
            [true, false, false, true, false],
            [false, true, false, true, false],
            [true, true, false, true, false],
            [false, false, true, true, false],
            [true, false, true, true, false],
            [false, true, true, true, false],
            [true, true, true, true, false],
            [false, false, false, false, true],
            [true, false, false, false, true],
            [false, true, false, false, true],
            [true, true, false, false, true],
            [false, false, true, false, true],
            [true, false, true, false, true],
            [false, true, true, false, true],
            [true, true, true, false, true],
            [false, false, false, true, true],
            [true, false, false, true, true],
            [false, true, false, true, true],
            [true, true, false, true, true],
            [false, false, true, true, true],
            [true, false, true, true, true],
            [false, true, true, true, true],
            [true, true, true, true, true],
        ];
        const BIT_IDS: [u64; 2] = [8, 23];
        const DESTINATION_COLLISIONS: [bool; 2] = [false, true];
        const CENTER_CASES: [bool; 2] = [false, true];

        let mut cases = 0;
        for (mask, expected_presence) in EXPECTED_PRESENCE.iter().enumerate() {
            for collision in DESTINATION_COLLISIONS {
                for repeated_centers in CENTER_CASES {
                    let source_atom_bits = if repeated_centers {
                        vec![vec![8, 8], vec![7, 8, 8], vec![8], vec![7]]
                    } else {
                        vec![vec![8, 8], vec![7], vec![7], vec![]]
                    };
                    let source_centers = if repeated_centers {
                        vec![(1, 0), (1, 2), (1, 0)]
                    } else {
                        vec![(1, 0)]
                    };
                    let source_atoms_per_bit = if repeated_centers {
                        vec![vec![1, 4], vec![1, 4]]
                    } else {
                        vec![vec![1, 4]]
                    };

                    let old_output = FingerprintAdditionalOutput {
                        atom_counts: expected_presence[0].then(|| vec![1, 2, 3, 4]),
                        atom_to_bits: expected_presence[1].then_some(source_atom_bits),
                        bit_info_map: expected_presence[2].then(|| {
                            BTreeMap::from([(BIT_IDS[0], source_centers), (91, vec![(9, 1)])])
                        }),
                        bit_paths: expected_presence[3].then(|| {
                            BTreeMap::from([
                                (BIT_IDS[0], vec![vec![1, 2], vec![2, 3]]),
                                (91, vec![vec![9]]),
                            ])
                        }),
                        atoms_per_bit: expected_presence[4]
                            .then(|| BTreeMap::from([(BIT_IDS[0], source_atoms_per_bit)])),
                    };
                    let mut new_output = FingerprintAdditionalOutput {
                        atom_counts: expected_presence[0].then(|| vec![71, 72, 73, 74]),
                        atom_to_bits: expected_presence[1].then(|| {
                            if collision {
                                vec![vec![23], vec![55], vec![23], vec![]]
                            } else {
                                vec![vec![], vec![], vec![], vec![]]
                            }
                        }),
                        bit_info_map: expected_presence[2].then(|| {
                            if collision {
                                BTreeMap::from([(BIT_IDS[1], vec![(90, 9)]), (91, vec![(91, 1)])])
                            } else {
                                BTreeMap::from([(91, vec![(91, 1)])])
                            }
                        }),
                        bit_paths: expected_presence[3].then(|| {
                            if collision {
                                BTreeMap::from([(BIT_IDS[1], vec![vec![90]]), (91, vec![vec![91]])])
                            } else {
                                BTreeMap::from([(91, vec![vec![91]])])
                            }
                        }),
                        atoms_per_bit: expected_presence[4].then(|| {
                            if collision {
                                BTreeMap::from([(BIT_IDS[1], vec![vec![90]]), (91, vec![vec![91]])])
                            } else {
                                BTreeMap::from([(91, vec![vec![91]])])
                            }
                        }),
                    };

                    old_output
                        .duplicate_bit_to(&mut new_output, BIT_IDS[0], BIT_IDS[1])
                        .unwrap();

                    assert_eq!(
                        new_output.atom_counts,
                        expected_presence[0].then(|| vec![71, 72, 73, 74]),
                        "atomCounts remains untouched, mask {mask:#07b}, collision={collision}, repeated={repeated_centers}"
                    );
                    let expected_atom_to_bits = match (collision, repeated_centers) {
                        (false, false) => vec![vec![23], vec![], vec![], vec![]],
                        (true, false) => vec![vec![23, 23], vec![55], vec![23], vec![]],
                        (false, true) => vec![vec![23], vec![23], vec![23], vec![]],
                        (true, true) => {
                            vec![vec![23, 23], vec![55, 23], vec![23, 23], vec![]]
                        }
                    };
                    assert_eq!(
                        new_output.atom_to_bits,
                        expected_presence[1].then_some(expected_atom_to_bits),
                        "atomToBits appends once per matching center without global dedup, mask {mask:#07b}, collision={collision}, repeated={repeated_centers}"
                    );
                    assert_eq!(
                        new_output.bit_info_map,
                        expected_presence[2].then(|| {
                            BTreeMap::from([
                                (
                                    BIT_IDS[1],
                                    if repeated_centers {
                                        vec![(1, 0), (1, 2), (1, 0)]
                                    } else {
                                        vec![(1, 0)]
                                    },
                                ),
                                (91, vec![(91, 1)]),
                            ])
                        }),
                        "bitInfoMap copies the complete ordered source value, mask {mask:#07b}, collision={collision}, repeated={repeated_centers}"
                    );
                    assert_eq!(
                        new_output.bit_paths,
                        expected_presence[3].then(|| {
                            BTreeMap::from([
                                (BIT_IDS[1], vec![vec![1, 2], vec![2, 3]]),
                                (91, vec![vec![91]]),
                            ])
                        }),
                        "bitPaths overwrites the colliding destination value, mask {mask:#07b}, collision={collision}, repeated={repeated_centers}"
                    );
                    assert_eq!(
                        new_output.atoms_per_bit,
                        expected_presence[4].then(|| {
                            BTreeMap::from([
                                (
                                    BIT_IDS[1],
                                    if repeated_centers {
                                        vec![vec![1, 4], vec![1, 4]]
                                    } else {
                                        vec![vec![1, 4]]
                                    },
                                ),
                                (91, vec![vec![91]]),
                            ])
                        }),
                        "atomsPerBit copies only the matching source key, mask {mask:#07b}, collision={collision}, repeated={repeated_centers}"
                    );
                    cases += 1;
                }
            }
        }
        assert_eq!(cases, 128);
    }

    #[test]
    fn fingerprint_morgan_g06_source_preconditions_keep_order_and_type() {
        let cases = [
            (
                FingerprintAdditionalOutput {
                    bit_info_map: Some(BTreeMap::new()),
                    ..FingerprintAdditionalOutput::default()
                },
                "bitInfoMap not allocated",
            ),
            (
                FingerprintAdditionalOutput {
                    atom_to_bits: Some(Vec::new()),
                    ..FingerprintAdditionalOutput::default()
                },
                "atomToBits not allocated",
            ),
            (
                FingerprintAdditionalOutput {
                    bit_paths: Some(BTreeMap::new()),
                    ..FingerprintAdditionalOutput::default()
                },
                "bitPaths not allocated",
            ),
        ];

        for (old_output, expected_message) in cases {
            let mut new_output = FingerprintAdditionalOutput::default();
            assert_eq!(
                old_output.duplicate_bit_to(&mut new_output, 8, 23),
                Err(FingerprintError::PreconditionViolation {
                    what: expected_message,
                })
            );
            assert_eq!(new_output, FingerprintAdditionalOutput::default());
        }
    }

    #[test]
    fn fingerprint_morgan_g06_missing_source_keys_leave_destination_unchanged() {
        let old_output = FingerprintAdditionalOutput {
            atom_counts: Some(vec![1, 2]),
            atom_to_bits: Some(vec![vec![7], vec![4]]),
            bit_info_map: Some(BTreeMap::from([(91, vec![(9, 1)])])),
            bit_paths: Some(BTreeMap::from([(91, vec![vec![9]])])),
            atoms_per_bit: Some(BTreeMap::from([(91, vec![vec![9]])])),
        };
        let mut new_output = FingerprintAdditionalOutput {
            atom_counts: Some(vec![71, 72]),
            atom_to_bits: Some(vec![vec![23], vec![55]]),
            bit_info_map: Some(BTreeMap::from([(23, vec![(90, 9)]), (91, vec![(91, 1)])])),
            bit_paths: Some(BTreeMap::from([(23, vec![vec![90]]), (91, vec![vec![91]])])),
            atoms_per_bit: Some(BTreeMap::from([(23, vec![vec![90]]), (91, vec![vec![91]])])),
        };

        old_output.duplicate_bit_to(&mut new_output, 8, 23).unwrap();

        assert_eq!(
            new_output,
            FingerprintAdditionalOutput {
                atom_counts: Some(vec![71, 72]),
                atom_to_bits: Some(vec![vec![23], vec![55]]),
                bit_info_map: Some(BTreeMap::from([(23, vec![(90, 9)]), (91, vec![(91, 1)]),])),
                bit_paths: Some(BTreeMap::from(
                    [(23, vec![vec![90]]), (91, vec![vec![91]]),]
                )),
                atoms_per_bit: Some(BTreeMap::from(
                    [(23, vec![vec![90]]), (91, vec![vec![91]]),]
                )),
            }
        );
    }

    #[test]
    fn fingerprint_morgan_g06_empty_atom_to_bits_destination_resizes() {
        let old_output = FingerprintAdditionalOutput {
            atom_to_bits: Some(vec![vec![8], vec![5], vec![8, 8]]),
            ..FingerprintAdditionalOutput::default()
        };
        let mut new_output = FingerprintAdditionalOutput {
            atom_to_bits: Some(Vec::new()),
            ..FingerprintAdditionalOutput::default()
        };

        old_output.duplicate_bit_to(&mut new_output, 8, 23).unwrap();

        assert_eq!(
            new_output.atom_to_bits,
            Some(vec![vec![23], vec![], vec![23]])
        );
    }
}
