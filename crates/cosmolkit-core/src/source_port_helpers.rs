#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CountSwapsError {
    SizeMismatch,
    MissingProbeElement,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RingNeighborError {
    MissingRing,
    RingIndexOutOfBounds,
}

pub(crate) fn rdkit_get_iso_map<E>(
    bonds: impl IntoIterator<Item = (crate::AtomId, crate::AtomId)>,
    mut atom_at: impl FnMut(crate::AtomId) -> Result<(u8, Option<u16>), E>,
) -> Result<Vec<(crate::AtomId, Vec<u16>)>, E> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/AddHs.cpp :: getIsoMap
    // Shared source-identical traversal only. Clearing `_isotopicHs` remains in
    // each caller's mutation adapter because COSMolKit does not expose mutable
    // atoms through this read-only helper.
    // RDKit✔️✔️:   std::map<unsigned int, std::vector<unsigned int>> isoMap;
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     auto ba = bond->getBeginAtom();
    // RDKit✔️✔️:     auto ea = bond->getEndAtom();
    // RDKit✔️✔️:     int ha = -1;
    // RDKit✔️✔️:     unsigned int iso;
    // RDKit✔️✔️:     if (ba->getAtomicNum() == 1 && ba->getIsotope() &&
    // RDKit✔️✔️:         ea->getAtomicNum() != 1) {
    // RDKit✔️✔️:       ha = ea->getIdx();
    // RDKit✔️✔️:       iso = ba->getIsotope();
    // RDKit✔️✔️:     } else if (ea->getAtomicNum() == 1 && ea->getIsotope() &&
    // RDKit✔️✔️:                ba->getAtomicNum() != 1) {
    // RDKit✔️✔️:       ha = ba->getIdx();
    // RDKit✔️✔️:       iso = ea->getIsotope();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (ha == -1) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto &v = isoMap[ha];
    // RDKit✔️✔️:     v.push_back(iso);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return isoMap;
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/AddHs.cpp :: getIsoMap
    let mut isotope_map = std::collections::BTreeMap::<crate::AtomId, Vec<u16>>::new();
    for (begin_id, end_id) in bonds {
        let (begin_atomic_number, begin_isotope) = atom_at(begin_id)?;
        let (end_atomic_number, end_isotope) = atom_at(end_id)?;
        let tracked =
            if begin_atomic_number == 1 && begin_isotope.is_some() && end_atomic_number != 1 {
                Some((end_id, begin_isotope.expect("checked above")))
            } else if end_atomic_number == 1 && end_isotope.is_some() && begin_atomic_number != 1 {
                Some((begin_id, end_isotope.expect("checked above")))
            } else {
                None
            };
        if let Some((heavy_atom, isotope)) = tracked {
            isotope_map.entry(heavy_atom).or_default().push(isotope);
        }
    }
    Ok(isotope_map.into_iter().collect())
}

pub(crate) type RingNeighborMap = std::collections::BTreeMap<usize, Vec<usize>>;

pub(crate) fn rdkit_make_ring_neighbor_map<'a, T: PartialEq + 'a>(
    ring_count: usize,
    mut ring_at: impl FnMut(usize) -> &'a [T],
    max_size: usize,
    max_overlap_size: usize,
) -> RingNeighborMap {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Aromaticity.cpp :: RingUtils::makeRingNeighborMap
    // RDKit✔️✔️: void makeRingNeighborMap(const VECT_INT_VECT &brings,
    // RDKit✔️✔️:                          INT_INT_VECT_MAP &neighMap, unsigned int maxSize,
    // RDKit✔️✔️:                          unsigned int maxOverlapSize) {
    // RDKit✔️✔️:   auto nrings = rdcast<int>(brings.size());
    // RDKit✔️✔️:   int i, j;
    // RDKit✔️✔️:   INT_VECT ring1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (i = 0; i < nrings; ++i) {
    // RDKit✔️✔️:     // create an empty INT_VECT at neighMap[i] if it does not yet exist
    // RDKit✔️✔️:     neighMap[i];
    // RDKit✔️✔️:     if (maxSize && brings[i].size() > maxSize) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ring1 = brings[i];
    // RDKit✔️✔️:     for (j = i + 1; j < nrings; ++j) {
    // RDKit✔️✔️:       if (maxSize && brings[j].size() > maxSize) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       INT_VECT inter;
    // RDKit✔️✔️:       Intersect(ring1, brings[j], inter);
    // RDKit✔️✔️:       if (inter.size() > 0 &&
    // RDKit✔️✔️:           (!maxOverlapSize || inter.size() <= maxOverlapSize)) {
    // RDKit✔️✔️:         neighMap[i].push_back(j);
    // RDKit✔️✔️:         neighMap[j].push_back(i);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Aromaticity.cpp :: RingUtils::makeRingNeighborMap
    //
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/RDGeneral/types.cpp :: Intersect
    // RDKit✔️✔️: void Intersect(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res) {
    // RDKit✔️✔️:   res.resize(0);
    // RDKit✔️✔️:   INT_VECT_CI ri;
    // RDKit✔️✔️:   for (ri = r1.begin(); ri != r1.end(); ri++) {
    // RDKit✔️✔️:     if (std::find(r2.begin(), r2.end(), (*ri)) != r2.end()) {
    // RDKit✔️✔️:       res.push_back(*ri);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/RDGeneral/types.cpp :: Intersect
    let mut neighbor_map = RingNeighborMap::new();
    for i in 0..ring_count {
        neighbor_map.entry(i).or_default();
        let ring1 = ring_at(i);
        if max_size != 0 && ring1.len() > max_size {
            continue;
        }
        for j in i + 1..ring_count {
            let ring2 = ring_at(j);
            if max_size != 0 && ring2.len() > max_size {
                continue;
            }
            let overlap = ring1.iter().filter(|item| ring2.contains(item)).count();
            if overlap > 0 && (max_overlap_size == 0 || overlap <= max_overlap_size) {
                neighbor_map.entry(i).or_default().push(j);
                neighbor_map.entry(j).or_default().push(i);
            }
        }
    }
    neighbor_map
}

pub(crate) fn rdkit_pick_fused_rings(
    current: usize,
    neighbor_map: &RingNeighborMap,
    result: &mut Vec<usize>,
    done: &mut [bool],
    depth: usize,
) -> Result<(), RingNeighborError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Aromaticity.cpp :: RingUtils::pickFusedRings
    // RDKit✔️✔️: void pickFusedRings(int curr, const INT_INT_VECT_MAP &neighMap, INT_VECT &res,
    // RDKit✔️✔️:                     boost::dynamic_bitset<> &done, int depth) {
    // RDKit✔️✔️:   auto pos = neighMap.find(curr);
    // RDKit✔️✔️:   PRECONDITION(pos != neighMap.end(), "bad argument");
    // RDKit✔️✔️:   done[curr] = 1;
    // RDKit✔️✔️:   res.push_back(curr);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto &neighs = pos->second;
    // RDKit✔️✔️:   for (int neigh : neighs) {
    // RDKit✔️✔️:     if (!done[neigh]) {
    // RDKit✔️✔️:       pickFusedRings(neigh, neighMap, res, done, depth + 1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Aromaticity.cpp :: RingUtils::pickFusedRings
    let neighbors = neighbor_map
        .get(&current)
        .ok_or(RingNeighborError::MissingRing)?;
    let current_done = done
        .get_mut(current)
        .ok_or(RingNeighborError::RingIndexOutOfBounds)?;
    *current_done = true;
    result.push(current);
    for &neighbor in neighbors {
        let neighbor_done = done
            .get(neighbor)
            .copied()
            .ok_or(RingNeighborError::RingIndexOutOfBounds)?;
        if !neighbor_done {
            rdkit_pick_fused_rings(neighbor, neighbor_map, result, done, depth + 1)?;
        }
    }
    let _ = depth;
    Ok(())
}

pub(crate) fn rdkit_set_ring_angle(hybridization: crate::Hybridization, ring_size: usize) -> f64 {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp :: DGeomHelpers::_setRingAngle
    // RDKit✔️✔️: void _setRingAngle(Atom::HybridizationType aHyb, unsigned int ringSize,
    // RDKit✔️✔️:                    double &angle) {
    // RDKit✔️✔️:   // NOTE: this assumes that all angles in a ring are equal. This is
    // RDKit✔️✔️:   // certainly not always the case, particular in aromatic rings with
    // RDKit✔️✔️:   // heteroatoms
    // RDKit✔️✔️:   // like s1cncc1. This led to GitHub55, which was fixed elsewhere.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if ((aHyb == Atom::SP2 && ringSize <= 8) || (ringSize == 3) ||
    // RDKit✔️✔️:       (ringSize == 4)) {
    // RDKit✔️✔️:     angle = M_PI * (1 - 2.0 / ringSize);
    // RDKit✔️✔️:   } else if (aHyb == Atom::SP3) {
    // RDKit✔️✔️:     if (ringSize == 5) {
    // RDKit✔️✔️:       angle = 104 * M_PI / 180;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       angle = 109.5 * M_PI / 180;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (aHyb == Atom::SP3D) {
    // RDKit✔️✔️:     angle = 105.0 * M_PI / 180;
    // RDKit✔️✔️:   } else if (aHyb == Atom::SP3D2) {
    // RDKit✔️✔️:     angle = 90.0 * M_PI / 180;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     angle = 120 * M_PI / 180;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp :: DGeomHelpers::_setRingAngle
    if (hybridization == crate::Hybridization::Sp2 && ring_size <= 8)
        || ring_size == 3
        || ring_size == 4
    {
        std::f64::consts::PI * (1.0 - 2.0 / ring_size as f64)
    } else if hybridization == crate::Hybridization::Sp3 {
        if ring_size == 5 {
            104.0 * std::f64::consts::PI / 180.0
        } else {
            109.5 * std::f64::consts::PI / 180.0
        }
    } else if hybridization == crate::Hybridization::Sp3d {
        105.0 * std::f64::consts::PI / 180.0
    } else if hybridization == crate::Hybridization::Sp3d2 {
        90.0 * std::f64::consts::PI / 180.0
    } else {
        120.0 * std::f64::consts::PI / 180.0
    }
}

pub(crate) fn count_swaps_to_interconvert<T: Copy + Eq>(
    reference: &[T],
    probe: &[T],
) -> Result<usize, CountSwapsError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/RDGeneral/utils.h :: countSwapsToInterconvert
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: unsigned int countSwapsToInterconvert(const T &ref, T probe) {
    // RDKit✔️✔️:   PRECONDITION(ref.size() == probe.size(), "size mismatch");
    // RDKit✔️✔️:   typename T::const_iterator refIt = ref.begin();
    // RDKit✔️✔️:   typename T::iterator probeIt = probe.begin();
    // RDKit✔️✔️:   typename T::iterator probeIt2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int nSwaps = 0;
    // RDKit✔️✔️:   while (refIt != ref.end()) {
    // RDKit✔️✔️:     if ((*probeIt) != (*refIt)) {
    // RDKit✔️✔️:       bool foundIt = false;
    // RDKit✔️✔️:       probeIt2 = probeIt;
    // RDKit✔️✔️:       while ((*probeIt2) != (*refIt) && probeIt2 != probe.end()) {
    // RDKit✔️✔️:         ++probeIt2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (probeIt2 != probe.end()) {
    // RDKit✔️✔️:         foundIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       CHECK_INVARIANT(foundIt, "could not find probe element");
    // RDKit✔️✔️:
    // RDKit✔️✔️:       std::swap(*probeIt, *probeIt2);
    // RDKit✔️✔️:       nSwaps++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++probeIt;
    // RDKit✔️✔️:     ++refIt;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nSwaps;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/RDGeneral/utils.h :: countSwapsToInterconvert
    if reference.len() != probe.len() {
        return Err(CountSwapsError::SizeMismatch);
    }

    let mut probe = probe.to_vec();
    let mut swaps = 0usize;
    for (index, expected) in reference.iter().copied().enumerate() {
        if probe[index] == expected {
            continue;
        }
        let found_index = probe[index..]
            .iter()
            .position(|value| *value == expected)
            .map(|offset| index + offset)
            .ok_or(CountSwapsError::MissingProbeElement)?;
        probe.swap(index, found_index);
        swaps += 1;
    }
    Ok(swaps)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_get_iso_map_covers_every_source_selection_branch() {
        let atoms = [
            (6, None),
            (1, Some(2)),
            (1, Some(3)),
            (1, None),
            (8, Some(18)),
        ];
        let bonds = [
            (crate::AtomId::new(1), crate::AtomId::new(0)),
            (crate::AtomId::new(0), crate::AtomId::new(2)),
            (crate::AtomId::new(3), crate::AtomId::new(0)),
            (crate::AtomId::new(0), crate::AtomId::new(4)),
            (crate::AtomId::new(1), crate::AtomId::new(2)),
        ];

        let isotope_map =
            rdkit_get_iso_map(bonds, |atom| Ok::<_, ()>(atoms[atom.index()])).unwrap();

        assert_eq!(isotope_map, vec![(crate::AtomId::new(0), vec![2, 3])]);
        assert_eq!(
            rdkit_get_iso_map::<()>([], |_| unreachable!()).unwrap(),
            Vec::<(crate::AtomId, Vec<u16>)>::new()
        );
    }

    #[test]
    fn shared_get_iso_map_preserves_std_map_key_and_bond_value_order() {
        let atoms = [(6, None), (8, None), (1, Some(2)), (1, Some(3))];
        let bonds = [
            (crate::AtomId::new(3), crate::AtomId::new(1)),
            (crate::AtomId::new(2), crate::AtomId::new(0)),
            (crate::AtomId::new(2), crate::AtomId::new(1)),
        ];

        let isotope_map =
            rdkit_get_iso_map(bonds, |atom| Ok::<_, ()>(atoms[atom.index()])).unwrap();

        assert_eq!(
            isotope_map,
            vec![
                (crate::AtomId::new(0), vec![2]),
                (crate::AtomId::new(1), vec![3, 2]),
            ]
        );
    }

    #[test]
    fn shared_get_iso_map_propagates_atom_adapter_failures_without_continuing() {
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        struct MissingAtom(crate::AtomId);

        let missing = crate::AtomId::new(4);
        let error = rdkit_get_iso_map([(crate::AtomId::new(0), missing)], |atom| {
            if atom == missing {
                Err(MissingAtom(atom))
            } else {
                Ok((6, None))
            }
        })
        .unwrap_err();

        assert_eq!(error, MissingAtom(missing));
    }

    #[test]
    fn shared_set_ring_angle_covers_every_source_branch() {
        let pi = std::f64::consts::PI;
        assert_eq!(
            rdkit_set_ring_angle(crate::Hybridization::Sp, 3),
            pi * (1.0 - 2.0 / 3.0)
        );
        assert_eq!(
            rdkit_set_ring_angle(crate::Hybridization::Sp3, 4),
            pi * (1.0 - 2.0 / 4.0)
        );
        assert_eq!(
            rdkit_set_ring_angle(crate::Hybridization::Sp2, 8),
            pi * 0.75
        );
        assert_eq!(
            rdkit_set_ring_angle(crate::Hybridization::Sp2, 9),
            120.0 * pi / 180.0
        );
        assert_eq!(
            rdkit_set_ring_angle(crate::Hybridization::Sp3, 5),
            104.0 * pi / 180.0
        );
        assert_eq!(
            rdkit_set_ring_angle(crate::Hybridization::Sp3, 6),
            109.5 * pi / 180.0
        );
        assert_eq!(
            rdkit_set_ring_angle(crate::Hybridization::Sp3d, 6),
            105.0 * pi / 180.0
        );
        assert_eq!(
            rdkit_set_ring_angle(crate::Hybridization::Sp3d2, 6),
            90.0 * pi / 180.0
        );
        assert_eq!(
            rdkit_set_ring_angle(crate::Hybridization::Unspecified, 6),
            120.0 * pi / 180.0
        );
    }

    #[test]
    fn shared_set_ring_angle_preserves_source_integer_boundaries() {
        assert!(
            rdkit_set_ring_angle(crate::Hybridization::Sp2, 0).is_infinite()
                && rdkit_set_ring_angle(crate::Hybridization::Sp2, 0).is_sign_negative()
        );
        assert_eq!(
            rdkit_set_ring_angle(crate::Hybridization::Sp2, usize::MAX),
            120.0_f64.to_radians()
        );
    }

    #[test]
    fn shared_count_swaps_kernel_covers_source_success_branches() {
        assert_eq!(count_swaps_to_interconvert::<u8>(&[], &[]), Ok(0));
        assert_eq!(count_swaps_to_interconvert(&[1, 2, 3], &[1, 2, 3]), Ok(0));
        assert_eq!(count_swaps_to_interconvert(&[1, 2, 3], &[2, 1, 3]), Ok(1));
        assert_eq!(count_swaps_to_interconvert(&[1, 2, 3], &[2, 3, 1]), Ok(2));
        assert_eq!(count_swaps_to_interconvert(&[1, 1, 2], &[2, 1, 1]), Ok(2));
    }

    #[test]
    fn shared_count_swaps_kernel_reports_source_contract_failures() {
        assert_eq!(
            count_swaps_to_interconvert(&[1, 2], &[1]),
            Err(CountSwapsError::SizeMismatch)
        );
        assert_eq!(
            count_swaps_to_interconvert(&[1, 2], &[1, 3]),
            Err(CountSwapsError::MissingProbeElement)
        );
    }
}
