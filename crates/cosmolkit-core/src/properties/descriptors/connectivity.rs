//! Connectivity descriptor implementation boundary.
//!
//! This module exclusively owns RDKit connectivity deltas, Chi descriptors,
//! Hall-Kier alpha, Kappa descriptors, and Phi. Graph path enumeration and
//! periodic-table data remain shared chemistry infrastructure rather than
//! descriptor-local copies.

use std::sync::Arc;

use super::{DescriptorError, DescriptorResult};
use crate::{Atom, Hybridization, Molecule};

pub(super) fn hk_deltas(molecule: &Molecule, force: bool) -> DescriptorResult<Arc<[f64]>> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::detail::hkDeltas
    // RDKit✔️✔️: void hkDeltas(const ROMol &mol, std::vector<double> &deltas, bool force) {
    // RDKit✔️✔️:   PRECONDITION(deltas.size() >= mol.getNumAtoms(), "bad vector size");
    // RDKit✔️✔️:   if (!force && mol.hasProp(common_properties::_connectivityHKDeltas)) {
    // RDKit✔️✔️:     mol.getProp(common_properties::_connectivityHKDeltas, deltas);
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if !force && let Some(deltas) = molecule.connectivity_hk_deltas_cache() {
        return Ok(deltas);
    }

    // RDKit✔️✔️:   const PeriodicTable *tbl = PeriodicTable::getTable();
    // RDKit✔️✔️:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️✔️:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️✔️:   while (atBegin != atEnd) {
    let mut deltas = Vec::with_capacity(molecule.num_atoms());
    for atom in molecule.atoms() {
        // RDKit✔️✔️:     const Atom *at = mol[*atBegin];
        // RDKit✔️✔️:     unsigned int n = at->getAtomicNum();
        let atomic_number = atom.atomic_number();
        // RDKit✔️✔️:     if (n <= 1) {
        // RDKit✔️✔️:       deltas[at->getIdx()] = 0;
        // RDKit✔️✔️:     } else if (n <= 10) {
        // RDKit✔️✔️:       deltas[at->getIdx()] = tbl->getNouterElecs(n) - at->getTotalNumHs();
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       deltas[at->getIdx()] =
        // RDKit✔️✔️:           double(tbl->getNouterElecs(n) - at->getTotalNumHs()) /
        // RDKit✔️✔️:           (n - tbl->getNouterElecs(n) - 1);
        // RDKit✔️✔️:     }
        let mut delta =
            if atomic_number <= 1 {
                0.0
            } else {
                let outer_electrons = crate::chemistry::valence::periodic_table_outer_electrons(atomic_number)
                    .map_err(|_| DescriptorError::Unsupported {
                        function: "hk_deltas",
                        rdkit_function: "PeriodicTable::getNouterElecs",
                    })?;
                let total_hydrogens = crate::chemistry::valence::total_num_hydrogens(molecule, atom.id(), false)
                    .map_err(|_| DescriptorError::Unsupported {
                        function: "hk_deltas",
                        rdkit_function: "Atom::getTotalNumHs",
                    })?;
                // Both source subtractions undergo C++ usual arithmetic conversion
                // to `unsigned int`; preserve that wrapping before conversion to
                // `double`, including malformed cached-valence edge states.
                let outer_electrons = u32::try_from(outer_electrons).map_err(|_| DescriptorError::Unsupported {
                    function: "hk_deltas",
                    rdkit_function: "PeriodicTable::getNouterElecs",
                })?;
                let valence_delta = outer_electrons.wrapping_sub(total_hydrogens);
                if atomic_number <= 10 {
                    f64::from(valence_delta)
                } else {
                    let denominator = u32::from(atomic_number).wrapping_sub(outer_electrons).wrapping_sub(1);
                    f64::from(valence_delta) / f64::from(denominator)
                }
            };
        // RDKit✔️✔️:     if (deltas[at->getIdx()] != 0.0) {
        // RDKit✔️✔️:       deltas[at->getIdx()] = 1. / sqrt(deltas[at->getIdx()]);
        // RDKit✔️✔️:     }
        if delta != 0.0 {
            delta = 1.0 / delta.sqrt();
        }
        deltas.push(delta);
        // RDKit✔️✔️:     ++atBegin;
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   mol.setProp(common_properties::_connectivityHKDeltas, deltas, true);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::detail::hkDeltas
    let deltas = Arc::<[f64]>::from(deltas);
    molecule.set_connectivity_hk_deltas_cache(Arc::clone(&deltas));
    Ok(deltas)
}

pub(super) fn n_vals(molecule: &Molecule, force: bool) -> DescriptorResult<Arc<[f64]>> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::detail::nVals
    // RDKit✔️✔️: void nVals(const ROMol &mol, std::vector<double> &nVs, bool force) {
    // RDKit✔️✔️:   PRECONDITION(nVs.size() >= mol.getNumAtoms(), "bad vector size");
    // RDKit✔️✔️:   if (!force && mol.hasProp(common_properties::_connectivityNVals)) {
    // RDKit✔️✔️:     mol.getProp(common_properties::_connectivityNVals, nVs);
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if !force && let Some(values) = molecule.connectivity_n_vals_cache() {
        return Ok(values);
    }

    // RDKit✔️✔️:   const PeriodicTable *tbl = PeriodicTable::getTable();
    // RDKit✔️✔️:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️✔️:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️✔️:   while (atBegin != atEnd) {
    let mut values = Vec::with_capacity(molecule.num_atoms());
    for atom in molecule.atoms() {
        // RDKit✔️✔️:     const Atom *at = mol[*atBegin];
        // RDKit✔️✔️:     double v = tbl->getNouterElecs(at->getAtomicNum()) - at->getTotalNumHs();
        let outer_electrons =
            crate::chemistry::valence::periodic_table_outer_electrons(atom.atomic_number()).map_err(|_| {
                DescriptorError::Unsupported {
                    function: "n_vals",
                    rdkit_function: "PeriodicTable::getNouterElecs",
                }
            })?;
        let outer_electrons = u32::try_from(outer_electrons).map_err(|_| DescriptorError::Unsupported {
            function: "n_vals",
            rdkit_function: "PeriodicTable::getNouterElecs",
        })?;
        let total_hydrogens =
            crate::chemistry::valence::total_num_hydrogens(molecule, atom.id(), false).map_err(|_| {
                DescriptorError::Unsupported {
                    function: "n_vals",
                    rdkit_function: "Atom::getTotalNumHs",
                }
            })?;
        // Source C++ converts the signed outer-electron count to `unsigned
        // int` before subtraction from the unsigned total-H count.
        let mut value = f64::from(outer_electrons.wrapping_sub(total_hydrogens));
        // RDKit✔️✔️:     if (v != 0.0) {
        // RDKit✔️✔️:       v = 1. / sqrt(v);
        // RDKit✔️✔️:     }
        if value != 0.0 {
            value = 1.0 / value.sqrt();
        }
        // RDKit✔️✔️:     nVs[at->getIdx()] = v;
        values.push(value);
        // RDKit✔️✔️:     ++atBegin;
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   mol.setProp(common_properties::_connectivityNVals, nVs, true);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::detail::nVals
    let values = Arc::<[f64]>::from(values);
    molecule.set_connectivity_n_vals_cache(Arc::clone(&values));
    Ok(values)
}

pub(super) fn get_alpha(atom: &Atom) -> (f64, bool) {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::detail::getAlpha
    // RDKit✔️✔️: double getAlpha(const Atom &atom, bool &found) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   found = false;
    // RDKit✔️✔️:   switch (atom.getAtomicNum()) {
    let (result, found) = match atom.atomic_number() {
        // RDKit✔️✔️:     case 1:
        // RDKit✔️✔️:       res = 0.0;
        // RDKit✔️✔️:       found = true;
        // RDKit✔️✔️:       break;
        1 => (0.0, true),
        // RDKit✔️✔️:     case 6:
        // RDKit✔️✔️:       switch (atom.getHybridization()) {
        // RDKit✔️✔️:         case Atom::SP:
        // RDKit✔️✔️:           res = -0.22;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         case Atom::SP2:
        // RDKit✔️✔️:           res = -0.13;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         default:
        // RDKit✔️✔️:           res = 0.00;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:       };
        // RDKit✔️✔️:       break;
        6 => match atom.hybridization() {
            Hybridization::Sp => (-0.22, true),
            Hybridization::Sp2 => (-0.13, true),
            _ => (0.00, true),
        },
        // RDKit✔️✔️:     case 7:
        // RDKit✔️✔️:       switch (atom.getHybridization()) {
        // RDKit✔️✔️:         case Atom::SP:
        // RDKit✔️✔️:           res = -0.29;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         case Atom::SP2:
        // RDKit✔️✔️:           res = -0.20;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         default:
        // RDKit✔️✔️:           res = -0.04;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:       };
        // RDKit✔️✔️:       break;
        7 => match atom.hybridization() {
            Hybridization::Sp => (-0.29, true),
            Hybridization::Sp2 => (-0.20, true),
            _ => (-0.04, true),
        },
        // RDKit✔️✔️:     case 8:
        // RDKit✔️✔️:       switch (atom.getHybridization()) {
        // RDKit✔️✔️:         case Atom::SP2:
        // RDKit✔️✔️:           res = -0.20;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         default:
        // RDKit✔️✔️:           res = -0.04;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:       };
        // RDKit✔️✔️:       break;
        8 => match atom.hybridization() {
            Hybridization::Sp2 => (-0.20, true),
            _ => (-0.04, true),
        },
        // RDKit✔️✔️:     case 9:
        // RDKit✔️✔️:       res = -0.07;
        // RDKit✔️✔️:       found = true;
        // RDKit✔️✔️:       break;
        9 => (-0.07, true),
        // RDKit✔️✔️:     case 15:
        // RDKit✔️✔️:       switch (atom.getHybridization()) {
        // RDKit✔️✔️:         case Atom::SP2:
        // RDKit✔️✔️:           res = 0.30;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         default:
        // RDKit✔️✔️:           res = 0.43;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:       };
        // RDKit✔️✔️:       break;
        15 => match atom.hybridization() {
            Hybridization::Sp2 => (0.30, true),
            _ => (0.43, true),
        },
        // RDKit✔️✔️:     case 16:
        // RDKit✔️✔️:       switch (atom.getHybridization()) {
        // RDKit✔️✔️:         case Atom::SP2:
        // RDKit✔️✔️:           res = 0.22;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         default:
        // RDKit✔️✔️:           res = 0.35;
        // RDKit✔️✔️:           found = true;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:       };
        // RDKit✔️✔️:       break;
        16 => match atom.hybridization() {
            Hybridization::Sp2 => (0.22, true),
            _ => (0.35, true),
        },
        // RDKit✔️✔️:     case 17:
        // RDKit✔️✔️:       res = 0.29;
        // RDKit✔️✔️:       found = true;
        // RDKit✔️✔️:       break;
        17 => (0.29, true),
        // RDKit✔️✔️:     case 35:
        // RDKit✔️✔️:       res = 0.48;
        // RDKit✔️✔️:       found = true;
        // RDKit✔️✔️:       break;
        35 => (0.48, true),
        // RDKit✔️✔️:     case 53:
        // RDKit✔️✔️:       res = 0.73;
        // RDKit✔️✔️:       found = true;
        // RDKit✔️✔️:       break;
        53 => (0.73, true),
        // RDKit✔️✔️:     default:
        // RDKit✔️✔️:       break;
        _ => (0.0, false),
        // RDKit✔️✔️:   }
    };
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::detail::getAlpha
    (result, found)
}

pub(super) fn calc_chi_nv(molecule: &Molecule, order: usize, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChiNv
    // RDKit✔️✔️: double calcChiNv(const ROMol &mol, unsigned int n, bool force) {
    // RDKit✔️✔️:   std::vector<double> hkDs(mol.getNumAtoms());
    // RDKit✔️✔️:   detail::hkDeltas(mol, hkDs, force);
    let hk_deltas = hk_deltas(molecule, force)?;
    // RDKit✔️✔️:   PATH_LIST ps = findAllPathsOfLengthN(mol, n + 1, false);
    let target_length = order.checked_add(1).ok_or(DescriptorError::Unsupported {
        function: "calc_chi_nv",
        rdkit_function: "Descriptors::calcChiNv",
    })?;
    let paths =
        crate::chemistry::subgraph::find_all_paths_of_length_n(molecule, target_length, false, false, -1, false)
            .map_err(|_| DescriptorError::Unsupported {
                function: "calc_chi_nv",
                rdkit_function: "findAllPathsOfLengthN",
            })?;
    // RDKit✔️✔️:   double res = 0.0;
    let mut result = 0.0;
    // RDKit✔️✔️:   for (const auto &p : ps) {
    for path in paths {
        // RDKit✔️✔️:     TEST_ASSERT(p.size() == n + 1);
        debug_assert_eq!(path.len(), target_length);
        // RDKit✔️✔️:     double accum = 1.0;
        let mut accumulator = 1.0;
        // RDKit✔️✔️:     for (unsigned int i = 0; i < n; ++i) {
        // RDKit✔️✔️:       accum *= hkDs[p[i]];
        // RDKit✔️✔️:     }
        for &atom_index in &path[..order] {
            accumulator *= hk_deltas[atom_index];
        }
        // RDKit✔️✔️:     // only push on the last element if this isn't a ring; this was github 463:
        // RDKit✔️✔️:     if (p[n] != p[0]) {
        // RDKit✔️✔️:       accum *= hkDs[p[n]];
        // RDKit✔️✔️:     }
        if path[order] != path[0] {
            accumulator *= hk_deltas[path[order]];
        }
        // RDKit✔️✔️:     res += accum;
        result += accumulator;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChiNv
    Ok(result)
}

pub(super) fn calc_chi_nn(molecule: &Molecule, order: usize, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChiNn
    // RDKit✔️✔️: double calcChiNn(const ROMol &mol, unsigned int n, bool force) {
    // RDKit✔️✔️:   std::vector<double> nVs(mol.getNumAtoms());
    // RDKit✔️✔️:   detail::nVals(mol, nVs, force);
    let n_vals = n_vals(molecule, force)?;
    // RDKit✔️✔️:   PATH_LIST ps = findAllPathsOfLengthN(mol, n + 1, false);
    let target_length = order.checked_add(1).ok_or(DescriptorError::Unsupported {
        function: "calc_chi_nn",
        rdkit_function: "Descriptors::calcChiNn",
    })?;
    let paths =
        crate::chemistry::subgraph::find_all_paths_of_length_n(molecule, target_length, false, false, -1, false)
            .map_err(|_| DescriptorError::Unsupported {
                function: "calc_chi_nn",
                rdkit_function: "findAllPathsOfLengthN",
            })?;
    // RDKit✔️✔️:   double res = 0.0;
    let mut result = 0.0;
    // RDKit✔️✔️:   for (const auto &p : ps) {
    for path in paths {
        // RDKit✔️✔️:     TEST_ASSERT(p.size() == n + 1);
        debug_assert_eq!(path.len(), target_length);
        // RDKit✔️✔️:     double accum = 1.0;
        let mut accumulator = 1.0;
        // RDKit✔️✔️:     for (unsigned int i = 0; i < n; ++i) {
        // RDKit✔️✔️:       accum *= nVs[p[i]];
        // RDKit✔️✔️:     }
        for &atom_index in &path[..order] {
            accumulator *= n_vals[atom_index];
        }
        // RDKit✔️✔️:     // only push on the last element if this isn't a ring; this was github 463:
        // RDKit✔️✔️:     if (p[n] != p[0]) {
        // RDKit✔️✔️:       accum *= nVs[p[n]];
        // RDKit✔️✔️:     }
        if path[order] != path[0] {
            accumulator *= n_vals[path[order]];
        }
        // RDKit✔️✔️:     res += accum;
        result += accumulator;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChiNn
    Ok(result)
}

pub(super) fn calc_chi_0v(molecule: &Molecule, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi0v
    // RDKit✔️✔️: double calcChi0v(const ROMol &mol, bool force) {
    // RDKit✔️✔️:   std::vector<double> hkDs(mol.getNumAtoms());
    // RDKit✔️✔️:   detail::hkDeltas(mol, hkDs, force);
    let hk_deltas = hk_deltas(molecule, force)?;
    // RDKit✔️✔️:   return std::accumulate(hkDs.begin(), hkDs.end(), 0.0);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi0v
    Ok(hk_deltas.iter().fold(0.0, |result, delta| result + delta))
}

pub(super) fn calc_chi_1v(molecule: &Molecule, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi1v
    // RDKit✔️✔️: double calcChi1v(const ROMol &mol, bool force) {
    // RDKit✔️✔️:   std::vector<double> hkDs(mol.getNumAtoms());
    // RDKit✔️✔️:   detail::hkDeltas(mol, hkDs, force);
    let hk_deltas = hk_deltas(molecule, force)?;
    // RDKit✔️✔️:   double res = 0.0;
    let mut result = 0.0;
    // RDKit✔️✔️:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️✔️:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit✔️✔️:   while (firstB != lastB) {
    for bond in molecule.bonds() {
        // RDKit✔️✔️:     const Bond *bond = mol[*firstB];
        // RDKit✔️✔️:     res += hkDs[bond->getBeginAtomIdx()] * hkDs[bond->getEndAtomIdx()];
        result += hk_deltas[bond.begin().index()] * hk_deltas[bond.end().index()];
        // RDKit✔️✔️:     ++firstB;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi1v
    Ok(result)
}

pub(super) fn calc_chi_2v(molecule: &Molecule, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi2v
    // RDKit✔️✔️: double calcChi2v(const ROMol &mol, bool force) {
    // RDKit✔️✔️:   return calcChiNv(mol, 2, force);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi2v
    calc_chi_nv(molecule, 2, force)
}

pub(super) fn calc_chi_3v(molecule: &Molecule, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi3v
    // RDKit✔️✔️: double calcChi3v(const ROMol &mol, bool force) {
    // RDKit✔️✔️:   return calcChiNv(mol, 3, force);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi3v
    calc_chi_nv(molecule, 3, force)
}

pub(super) fn calc_chi_4v(molecule: &Molecule, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi4v
    // RDKit✔️✔️: double calcChi4v(const ROMol &mol, bool force) {
    // RDKit✔️✔️:   return calcChiNv(mol, 4, force);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi4v
    calc_chi_nv(molecule, 4, force)
}

pub(super) fn calc_chi_0n(molecule: &Molecule, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi0n
    // RDKit✔️✔️: double calcChi0n(const ROMol &mol, bool force) {
    // RDKit✔️✔️:   std::vector<double> nVs(mol.getNumAtoms());
    // RDKit✔️✔️:   detail::nVals(mol, nVs, force);
    let n_vals = n_vals(molecule, force)?;
    // RDKit✔️✔️:   return std::accumulate(nVs.begin(), nVs.end(), 0.0);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi0n
    Ok(n_vals.iter().fold(0.0, |result, value| result + value))
}

pub(super) fn calc_chi_1n(molecule: &Molecule, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi1n
    // RDKit✔️✔️: double calcChi1n(const ROMol &mol, bool force) {
    // RDKit✔️✔️:   std::vector<double> nVs(mol.getNumAtoms());
    // RDKit✔️✔️:   detail::nVals(mol, nVs, force);
    let n_vals = n_vals(molecule, force)?;
    // RDKit✔️✔️:   double res = 0.0;
    let mut result = 0.0;
    // RDKit✔️✔️:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️✔️:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit✔️✔️:   while (firstB != lastB) {
    for bond in molecule.bonds() {
        // RDKit✔️✔️:     const Bond *bond = mol[*firstB];
        // RDKit✔️✔️:     res += nVs[bond->getBeginAtomIdx()] * nVs[bond->getEndAtomIdx()];
        result += n_vals[bond.begin().index()] * n_vals[bond.end().index()];
        // RDKit✔️✔️:     ++firstB;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi1n
    Ok(result)
}

pub(super) fn calc_chi_2n(molecule: &Molecule, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi2n
    // RDKit✔️✔️: double calcChi2n(const ROMol &mol, bool force) {
    // RDKit✔️✔️:   return calcChiNn(mol, 2, force);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi2n
    calc_chi_nn(molecule, 2, force)
}

pub(super) fn calc_chi_3n(molecule: &Molecule, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi3n
    // RDKit✔️✔️: double calcChi3n(const ROMol &mol, bool force) {
    // RDKit✔️✔️:   return calcChiNn(mol, 3, force);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi3n
    calc_chi_nn(molecule, 3, force)
}

pub(super) fn calc_chi_4n(molecule: &Molecule, force: bool) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi4n
    // RDKit✔️✔️: double calcChi4n(const ROMol &mol, bool force) {
    // RDKit✔️✔️:   return calcChiNn(mol, 4, force);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcChi4n
    calc_chi_nn(molecule, 4, force)
}

pub(super) fn calc_chi_0(molecule: &Molecule) -> f64 {
    // BEGIN RDKIT PYTHON FUNCTION: rdkit.Chem.GraphDescriptors.Chi0
    // RDKit✔️✔️: def Chi0(mol):
    // RDKit✔️✔️:   """ From equations (1),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    // RDKit✔️✔️:
    // RDKit✔️✔️:   """
    // RDKit✔️✔️:   deltas = [x.GetDegree() for x in mol.GetAtoms()]
    // RDKit✔️✔️:   while 0 in deltas:
    // RDKit✔️✔️:     deltas.remove(0)
    // RDKit✔️✔️:   deltas = numpy.array(deltas, float)
    // RDKit✔️✔️:   res = sum(numpy.sqrt(1. / deltas))
    // RDKit✔️✔️:   return res
    // END RDKIT PYTHON FUNCTION: rdkit.Chem.GraphDescriptors.Chi0
    molecule
        .atoms()
        .iter()
        .enumerate()
        .map(|(atom_index, _)| molecule.adjacency().neighbors_of(atom_index).len())
        .filter(|degree| *degree != 0)
        .map(|degree| (1.0 / degree as f64).sqrt())
        .fold(0.0, |result, contribution| result + contribution)
}

pub(super) fn calc_chi_1(molecule: &Molecule) -> f64 {
    // BEGIN RDKIT PYTHON FUNCTION: rdkit.Chem.GraphDescriptors.Chi1
    // RDKit✔️✔️: def Chi1(mol):
    // RDKit✔️✔️:   """ From equations (1),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    // RDKit✔️✔️:
    // RDKit✔️✔️:   """
    // RDKit✔️✔️:   c1s = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    // RDKit✔️✔️:   while 0 in c1s:
    // RDKit✔️✔️:     c1s.remove(0)
    // RDKit✔️✔️:   c1s = numpy.array(c1s, float)
    // RDKit✔️✔️:   res = sum(numpy.sqrt(1. / c1s))
    // RDKit✔️✔️:   return res
    // END RDKIT PYTHON FUNCTION: rdkit.Chem.GraphDescriptors.Chi1
    molecule
        .bonds()
        .iter()
        .map(|bond| {
            molecule.adjacency().neighbors_of(bond.begin().index()).len()
                * molecule.adjacency().neighbors_of(bond.end().index()).len()
        })
        .filter(|degree_product| *degree_product != 0)
        .map(|degree_product| (1.0 / degree_product as f64).sqrt())
        .fold(0.0, |result, contribution| result + contribution)
}

pub(super) fn calc_hall_kier_alpha(molecule: &Molecule, mut atom_contributions: Option<&mut [f64]>) -> f64 {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcHallKierAlpha
    // RDKit✔️✔️: double calcHallKierAlpha(const ROMol &mol, std::vector<double> *atomContribs) {
    // RDKit✔️✔️:   PRECONDITION(!atomContribs || atomContribs->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atomContribs vector");
    assert!(
        atom_contributions
            .as_ref()
            .is_none_or(|contributions| contributions.len() >= molecule.num_atoms()),
        "bad atomContribs vector"
    );
    // RDKit✔️✔️:   const PeriodicTable *tbl = PeriodicTable::getTable();
    // RDKit✔️✔️:   double alphaSum = 0.0;
    let mut alpha_sum = 0.0;
    // RDKit✔️✔️:   double rC = tbl->getRb0(6);
    let carbon_radius = crate::chemistry::valence::rdkit_rb0(6);
    // RDKit✔️✔️:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️✔️:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️✔️:   while (atBegin != atEnd) {
    for atom in molecule.atoms() {
        // RDKit✔️✔️:     const Atom *at = mol[*atBegin];
        // RDKit✔️✔️:     ++atBegin;
        // RDKit✔️✔️:     unsigned int n = at->getAtomicNum();
        let atomic_number = atom.atomic_number();
        // RDKit✔️✔️:     if (!n) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if atomic_number == 0 {
            continue;
        }
        // RDKit✔️✔️:     bool found;
        // RDKit✔️✔️:     double alpha = detail::getAlpha(*(at), found);
        let (mut alpha, found) = get_alpha(atom);
        // RDKit✔️✔️:     if (!found) {
        // RDKit✔️✔️:       double rA = tbl->getRb0(n);
        // RDKit✔️✔️:       alpha = rA / rC - 1.0;
        // RDKit✔️✔️:     }
        if !found {
            let atomic_radius = crate::chemistry::valence::rdkit_rb0(atomic_number);
            alpha = atomic_radius / carbon_radius - 1.0;
        }
        // RDKit✔️✔️:     alphaSum += alpha;
        alpha_sum += alpha;
        // RDKit✔️✔️:     if (atomContribs) {
        // RDKit✔️✔️:       (*atomContribs)[at->getIdx()] = alpha;
        // RDKit✔️✔️:     }
        if let Some(contributions) = atom_contributions.as_deref_mut() {
            contributions[atom.id().index()] = alpha;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return alphaSum;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcHallKierAlpha
    alpha_sum
}

fn kappa_1_helper(path_count: f64, heavy_atom_count: f64, alpha: f64) -> f64 {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::kappa1Helper
    // RDKit✔️✔️: double kappa1Helper(double P1, double A, double alpha) {
    // RDKit✔️✔️:   double denom = P1 + alpha;
    let denominator = path_count + alpha;
    // RDKit✔️✔️:   double kappa = 0.0;
    let mut kappa = 0.0;
    // RDKit✔️✔️:   if (denom) {
    // RDKit✔️✔️:     kappa = (A + alpha) * (A + alpha - 1) * (A + alpha - 1) / (denom * denom);
    // RDKit✔️✔️:   }
    if denominator != 0.0 {
        kappa = (heavy_atom_count + alpha) * (heavy_atom_count + alpha - 1.0) * (heavy_atom_count + alpha - 1.0)
            / (denominator * denominator);
    }
    // RDKit✔️✔️:   return kappa;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::kappa1Helper
    kappa
}

fn kappa_2_helper(path_count: f64, heavy_atom_count: f64, alpha: f64) -> f64 {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::kappa2Helper
    // RDKit✔️✔️: double kappa2Helper(double P2, double A, double alpha) {
    // RDKit✔️✔️:   double denom = (P2 + alpha) * (P2 + alpha);
    let denominator = (path_count + alpha) * (path_count + alpha);
    // RDKit✔️✔️:   double kappa = 0.0;
    let mut kappa = 0.0;
    // RDKit✔️✔️:   if (denom) {
    // RDKit✔️✔️:     kappa = (A + alpha - 1) * (A + alpha - 2) * (A + alpha - 2) / denom;
    // RDKit✔️✔️:   }
    if denominator != 0.0 {
        kappa = (heavy_atom_count + alpha - 1.0) * (heavy_atom_count + alpha - 2.0) * (heavy_atom_count + alpha - 2.0)
            / denominator;
    }
    // RDKit✔️✔️:   return kappa;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::kappa2Helper
    kappa
}

fn kappa_3_helper(path_count: f64, heavy_atom_count: i32, alpha: f64) -> f64 {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::kappa3Helper
    // RDKit✔️✔️: double kappa3Helper(double P3, int A, double alpha) {
    // RDKit✔️✔️:   double denom = (P3 + alpha) * (P3 + alpha);
    let denominator = (path_count + alpha) * (path_count + alpha);
    // RDKit✔️✔️:   double kappa = 0.0;
    let mut kappa = 0.0;
    // RDKit✔️✔️:   if (denom) {
    if denominator != 0.0 {
        let heavy_atom_count_as_double = f64::from(heavy_atom_count);
        // RDKit✔️✔️:     if (A % 2) {
        // RDKit✔️✔️:       kappa = (A + alpha - 1) * (A + alpha - 3) * (A + alpha - 3) / denom;
        if heavy_atom_count % 2 != 0 {
            kappa = (heavy_atom_count_as_double + alpha - 1.0)
                * (heavy_atom_count_as_double + alpha - 3.0)
                * (heavy_atom_count_as_double + alpha - 3.0)
                / denominator;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       kappa = (A + alpha - 2) * (A + alpha - 3) * (A + alpha - 3) / denom;
        } else {
            kappa = (heavy_atom_count_as_double + alpha - 2.0)
                * (heavy_atom_count_as_double + alpha - 3.0)
                * (heavy_atom_count_as_double + alpha - 3.0)
                / denominator;
        }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return kappa;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::kappa3Helper
    kappa
}

pub(super) fn calc_kappa_1(molecule: &Molecule) -> f64 {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcKappa1
    // RDKit✔️✔️: double calcKappa1(const ROMol &mol) {
    // RDKit✔️✔️:   double P1 = mol.getNumBonds();
    let path_count = molecule.num_bonds() as f64;
    // RDKit✔️✔️:   double A = mol.getNumHeavyAtoms();
    let heavy_atom_count = molecule.atoms().iter().filter(|atom| atom.atomic_number() > 1).count() as f64;
    // RDKit✔️✔️:   double alpha = calcHallKierAlpha(mol);
    let alpha = calc_hall_kier_alpha(molecule, None);
    // RDKit✔️✔️:   double kappa = kappa1Helper(P1, A, alpha);
    let kappa = kappa_1_helper(path_count, heavy_atom_count, alpha);
    // RDKit✔️✔️:   return kappa;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcKappa1
    kappa
}

pub(super) fn calc_kappa_2(molecule: &Molecule) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcKappa2
    // RDKit✔️✔️: double calcKappa2(const ROMol &mol) {
    // RDKit✔️✔️:   PATH_LIST ps = findAllPathsOfLengthN(mol, 2);
    let paths =
        crate::chemistry::subgraph::find_all_paths_of_length_n(molecule, 2, true, false, -1, false).map_err(|_| {
            DescriptorError::Unsupported {
                function: "calc_kappa_2",
                rdkit_function: "findAllPathsOfLengthN",
            }
        })?;
    // RDKit✔️✔️:   double P2 = ps.size();
    let path_count = paths.len() as f64;
    // RDKit✔️✔️:   double A = mol.getNumHeavyAtoms();
    let heavy_atom_count = molecule.atoms().iter().filter(|atom| atom.atomic_number() > 1).count() as f64;
    // RDKit✔️✔️:   double alpha = calcHallKierAlpha(mol);
    let alpha = calc_hall_kier_alpha(molecule, None);
    // RDKit✔️✔️:   double kappa = kappa2Helper(P2, A, alpha);
    let kappa = kappa_2_helper(path_count, heavy_atom_count, alpha);
    // RDKit✔️✔️:   return kappa;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcKappa2
    Ok(kappa)
}

pub(super) fn calc_kappa_3(molecule: &Molecule) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcKappa3
    // RDKit✔️✔️: double calcKappa3(const ROMol &mol) {
    // RDKit✔️✔️:   double P3 = findAllPathsOfLengthN(mol, 3).size();
    let path_count = crate::chemistry::subgraph::find_all_paths_of_length_n(molecule, 3, true, false, -1, false)
        .map_err(|_| DescriptorError::Unsupported {
            function: "calc_kappa_3",
            rdkit_function: "findAllPathsOfLengthN",
        })?
        .len() as f64;
    // RDKit✔️✔️:   int A = mol.getNumHeavyAtoms();
    let heavy_atom_count = i32::try_from(molecule.atoms().iter().filter(|atom| atom.atomic_number() > 1).count())
        .map_err(|_| DescriptorError::Unsupported {
            function: "calc_kappa_3",
            rdkit_function: "ROMol::getNumHeavyAtoms int conversion",
        })?;
    // RDKit✔️✔️:   double alpha = calcHallKierAlpha(mol);
    let alpha = calc_hall_kier_alpha(molecule, None);
    // RDKit✔️✔️:   double kappa = kappa3Helper(P3, A, alpha);
    let kappa = kappa_3_helper(path_count, heavy_atom_count, alpha);
    // RDKit✔️✔️:   return kappa;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcKappa3
    Ok(kappa)
}

pub(super) fn calc_phi(molecule: &Molecule) -> DescriptorResult<f64> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcPhi
    // RDKit✔️✔️: double calcPhi(const ROMol &mol) {
    let heavy_atom_count = molecule.atoms().iter().filter(|atom| atom.atomic_number() > 1).count();
    // RDKit✔️✔️:   if (!mol.getNumHeavyAtoms()) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    if heavy_atom_count == 0 {
        return Ok(0.0);
    }
    // RDKit✔️✔️:   auto alpha = calcHallKierAlpha(mol);
    let alpha = calc_hall_kier_alpha(molecule, None);
    // RDKit✔️✔️:   auto P1 = mol.getNumBonds();
    let path_count_1 = molecule.num_bonds() as f64;
    // RDKit✔️✔️:   auto A = mol.getNumHeavyAtoms();
    let heavy_atom_count = heavy_atom_count as f64;
    // RDKit✔️✔️:   auto kappa1 = kappa1Helper(P1, A, alpha);
    let kappa_1 = kappa_1_helper(path_count_1, heavy_atom_count, alpha);
    // RDKit✔️✔️:   auto P2 = findAllPathsOfLengthN(mol, 2).size();
    let path_count_2 = crate::chemistry::subgraph::find_all_paths_of_length_n(molecule, 2, true, false, -1, false)
        .map_err(|_| DescriptorError::Unsupported {
            function: "calc_phi",
            rdkit_function: "findAllPathsOfLengthN",
        })?
        .len() as f64;
    // RDKit✔️✔️:   auto kappa2 = kappa2Helper(P2, A, alpha);
    let kappa_2 = kappa_2_helper(path_count_2, heavy_atom_count, alpha);
    // RDKit✔️✔️:   auto Phi = kappa1 * kappa2 / A;
    let phi = kappa_1 * kappa_2 / heavy_atom_count;
    // RDKit✔️✔️:   return Phi;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcPhi
    Ok(phi)
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::{
        calc_chi_0, calc_chi_0n, calc_chi_0v, calc_chi_1, calc_chi_1n, calc_chi_1v, calc_chi_2n, calc_chi_2v,
        calc_chi_3n, calc_chi_3v, calc_chi_4n, calc_chi_4v, calc_chi_nn, calc_chi_nv, calc_hall_kier_alpha,
        calc_kappa_1, calc_kappa_2, calc_kappa_3, calc_phi, get_alpha, hk_deltas, kappa_1_helper, kappa_2_helper,
        kappa_3_helper, n_vals,
    };
    use crate::{Atom, AtomId, AtomSpec, BondOrder, BondSpec, Element, Hybridization, Molecule, MoleculeBuilder};

    fn atom(element: Element, hybridization: Hybridization) -> Atom {
        Atom::from_spec(AtomId::new(0), AtomSpec::new(element).with_hybridization(hybridization))
    }

    fn assert_f64_bits(actual: f64, expected: f64, context: &str) {
        assert_eq!(
            actual.to_bits(),
            expected.to_bits(),
            "{context}: actual={actual:?} expected={expected:?}"
        );
    }

    fn assert_slice_bits(actual: &[f64], expected: &[f64], context: &str) {
        assert_eq!(actual.len(), expected.len(), "{context}: length");
        for (index, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
            assert_f64_bits(actual, expected, &format!("{context}[{index}]"));
        }
    }

    fn explicit_methane() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        for _ in 0..4 {
            let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
            builder
                .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
                .expect("explicit methane bond");
        }
        builder
            .build()
            .expect("explicit methane graph")
            .sanitize()
            .expect("sanitized explicit methane")
    }

    #[test]
    fn get_alpha_covers_every_source_element_hybridization_and_fallback_branch() {
        let cases = [
            (Element::H, Hybridization::Unspecified, 0.0, true),
            (Element::C, Hybridization::Sp, -0.22, true),
            (Element::C, Hybridization::Sp2, -0.13, true),
            (Element::C, Hybridization::Sp3, 0.00, true),
            (Element::N, Hybridization::Sp, -0.29, true),
            (Element::N, Hybridization::Sp2, -0.20, true),
            (Element::N, Hybridization::Sp3, -0.04, true),
            (Element::O, Hybridization::Sp2, -0.20, true),
            (Element::O, Hybridization::Sp3, -0.04, true),
            (Element::F, Hybridization::Sp3, -0.07, true),
            (Element::P, Hybridization::Sp2, 0.30, true),
            (Element::P, Hybridization::Sp3, 0.43, true),
            (Element::S, Hybridization::Sp2, 0.22, true),
            (Element::S, Hybridization::Sp3, 0.35, true),
            (Element::CL, Hybridization::Sp3, 0.29, true),
            (Element::BR, Hybridization::Sp3, 0.48, true),
            (Element::I, Hybridization::Sp3, 0.73, true),
            (Element::DUMMY, Hybridization::Unspecified, 0.0, false),
            (Element::XE, Hybridization::Sp3, 0.0, false),
        ];

        for (element, hybridization, expected, expected_found) in cases {
            let (actual, found) = get_alpha(&atom(element, hybridization));
            assert_f64_bits(actual, expected, element.symbol());
            assert_eq!(found, expected_found, "{} found state", element.symbol());
        }
    }

    #[test]
    fn connectivity_deltas_match_source_total_h_and_periodic_row_branches() {
        let ethanol = Molecule::from_smiles("CCO").expect("ethanol");
        let inverse_sqrt_2 = 1.0 / 2.0_f64.sqrt();
        let inverse_sqrt_5 = 1.0 / 5.0_f64.sqrt();
        assert_slice_bits(
            &hk_deltas(&ethanol, true).unwrap(),
            &[1.0, inverse_sqrt_2, inverse_sqrt_5],
            "ethanol HK deltas",
        );
        assert_slice_bits(
            &n_vals(&ethanol, true).unwrap(),
            &[1.0, inverse_sqrt_2, inverse_sqrt_5],
            "ethanol nVals",
        );

        let explicit_hydrogens = Molecule::from_smiles("[CH2]").expect("explicit-H carbon");
        assert_slice_bits(
            &hk_deltas(&explicit_hydrogens, true).unwrap(),
            &[inverse_sqrt_2],
            "explicit-H HK delta",
        );
        assert_slice_bits(
            &n_vals(&explicit_hydrogens, true).unwrap(),
            &[inverse_sqrt_2],
            "explicit-H nVal",
        );

        let periodic_rows = Molecule::from_smiles("*.[Li+].[Na+].[K+]").expect("periodic rows");
        assert_slice_bits(
            &hk_deltas(&periodic_rows, true).unwrap(),
            &[0.0, 1.0, 3.0, 17.0_f64.sqrt()],
            "periodic-row HK deltas",
        );
        assert_slice_bits(
            &n_vals(&periodic_rows, true).unwrap(),
            &[0.0, 1.0, 1.0, 1.0],
            "periodic-row nVals",
        );
    }

    #[test]
    fn connectivity_delta_caches_obey_cold_warm_force_and_family_independence() {
        let molecule = Molecule::from_smiles("CCO").expect("cache molecule");
        assert!(molecule.connectivity_hk_deltas_cache().is_none());
        assert!(molecule.connectivity_n_vals_cache().is_none());

        let hk_cold = hk_deltas(&molecule, false).unwrap();
        let hk_warm = hk_deltas(&molecule, false).unwrap();
        assert!(Arc::ptr_eq(&hk_cold, &hk_warm));
        assert!(molecule.connectivity_n_vals_cache().is_none());

        let n_cold = n_vals(&molecule, false).unwrap();
        let n_warm = n_vals(&molecule, false).unwrap();
        assert!(Arc::ptr_eq(&n_cold, &n_warm));
        assert!(Arc::ptr_eq(&hk_cold, &molecule.connectivity_hk_deltas_cache().unwrap()));

        molecule.set_connectivity_hk_deltas_cache(Arc::<[f64]>::from([9.0, 9.0, 9.0]));
        assert_slice_bits(
            &hk_deltas(&molecule, false).unwrap(),
            &[9.0, 9.0, 9.0],
            "warm injected HK cache",
        );
        let hk_forced = hk_deltas(&molecule, true).unwrap();
        assert_slice_bits(&hk_forced, &hk_cold, "forced HK recomputation");
        assert!(!Arc::ptr_eq(&hk_forced, &hk_cold));

        molecule.set_connectivity_n_vals_cache(Arc::<[f64]>::from([8.0, 8.0, 8.0]));
        assert_slice_bits(
            &n_vals(&molecule, false).unwrap(),
            &[8.0, 8.0, 8.0],
            "warm injected nVal cache",
        );
        let n_forced = n_vals(&molecule, true).unwrap();
        assert_slice_bits(&n_forced, &n_cold, "forced nVal recomputation");
        assert!(!Arc::ptr_eq(&n_forced, &n_cold));
    }

    #[test]
    fn chi_descriptors_match_pinned_rdkit_across_path_and_element_branches() {
        #[derive(Clone, Copy)]
        enum Fixture {
            Smiles(&'static str),
            ExplicitMethane,
        }

        const CASES: [(&str, Fixture, [u64; 12]); 5] = [
            (
                "open path",
                Fixture::Smiles("CCCCCC"),
                [
                    0x4013_504f_333f_9de6,
                    0x4007_504f_333f_9de6,
                    0x4013_504f_333f_9de6,
                    0x4007_504f_333f_9de6,
                    0x3ffb_504f_333f_9de4,
                    0x3fee_a09e_667f_3bca,
                    0x3fdf_ffff_ffff_fffd,
                    0x4013_504f_333f_9de6,
                    0x4007_504f_333f_9de6,
                    0x3ffb_504f_333f_9de4,
                    0x3fee_a09e_667f_3bca,
                    0x3fdf_ffff_ffff_fffd,
                ],
            ),
            (
                "closed three-membered ring",
                Fixture::Smiles("C1CC1"),
                [
                    0x4000_f876_ccdf_6cda,
                    0x3ff8_0000_0000_0000,
                    0x4000_f876_ccdf_6cd9,
                    0x3ff7_ffff_ffff_fffe,
                    0x3ff0_f876_ccdf_6cd8,
                    0x3fd6_a09e_667f_3bcb,
                    0x0000_0000_0000_0000,
                    0x4000_f876_ccdf_6cd9,
                    0x3ff7_ffff_ffff_fffe,
                    0x3ff0_f876_ccdf_6cd8,
                    0x3fd6_a09e_667f_3bcb,
                    0x0000_0000_0000_0000,
                ],
            ),
            (
                "disconnected graph",
                Fixture::Smiles("CC.CCC"),
                [
                    0x4012_d413_cccf_e77a,
                    0x4003_504f_333f_9de6,
                    0x4012_d413_cccf_e77a,
                    0x4003_504f_333f_9de6,
                    0x3fe6_a09e_667f_3bcc,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x4012_d413_cccf_e77a,
                    0x4003_504f_333f_9de6,
                    0x3fe6_a09e_667f_3bcc,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                ],
            ),
            (
                "explicit hydrogens",
                Fixture::ExplicitMethane,
                [
                    0x4012_0000_0000_0000,
                    0x4000_0000_0000_0000,
                    0x3fe0_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x4012_0000_0000_0000,
                    0x4000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                ],
            ),
            (
                "heavy elements",
                Fixture::Smiles("[SiH3][GeH3]"),
                [
                    0x4000_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                    0x4020_646e_1721_1cc0,
                    0x402f_2d4a_4563_5640,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x4000_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                    0x0000_0000_0000_0000,
                ],
            ),
        ];

        for (case, fixture, expected) in CASES {
            let molecule = match fixture {
                Fixture::Smiles(smiles) => Molecule::from_smiles(smiles)
                    .unwrap_or_else(|error| panic!("failed to parse {case} fixture {smiles:?}: {error}")),
                Fixture::ExplicitMethane => explicit_methane(),
            };
            let actual = [
                calc_chi_0(&molecule),
                calc_chi_1(&molecule),
                calc_chi_0v(&molecule, true).unwrap(),
                calc_chi_1v(&molecule, true).unwrap(),
                calc_chi_2v(&molecule, true).unwrap(),
                calc_chi_3v(&molecule, true).unwrap(),
                calc_chi_4v(&molecule, true).unwrap(),
                calc_chi_0n(&molecule, true).unwrap(),
                calc_chi_1n(&molecule, true).unwrap(),
                calc_chi_2n(&molecule, true).unwrap(),
                calc_chi_3n(&molecule, true).unwrap(),
                calc_chi_4n(&molecule, true).unwrap(),
            ];
            for (index, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
                assert_eq!(
                    actual.to_bits(),
                    expected,
                    "{case} fixed Chi field {index}: actual={actual:?}"
                );
            }
            for (order, expected_index) in [(2, 4), (3, 5), (4, 6)] {
                let actual = calc_chi_nv(&molecule, order, true).unwrap();
                assert_eq!(
                    actual.to_bits(),
                    expected[expected_index],
                    "{case} generic valence Chi order {order}: actual={actual:?}"
                );
            }
            for (order, expected_index) in [(2, 9), (3, 10), (4, 11)] {
                let actual = calc_chi_nn(&molecule, order, true).unwrap();
                assert_eq!(
                    actual.to_bits(),
                    expected[expected_index],
                    "{case} generic nVal Chi order {order}: actual={actual:?}"
                );
            }
        }
    }

    #[test]
    fn hall_kier_alpha_preserves_atom_order_dummy_skip_and_rb0_fallback() {
        let molecule = Molecule::from_smiles("*.CP(=O)(O)SCl.[Xe]").expect("alpha fixture");
        let mut contributions = vec![99.0; molecule.num_atoms()];
        let actual = calc_hall_kier_alpha(&molecule, Some(&mut contributions));

        assert_eq!(actual.to_bits(), 0x4003_3620_2ecf_b9c8);
        assert_slice_bits(
            &contributions,
            &[99.0, 0.0, 0.43, -0.20, -0.04, 0.35, 0.29, 1.5714285714285712],
            "Hall-Kier atom contributions",
        );
    }

    #[test]
    fn kappa_helpers_cover_zero_denominators_and_kappa3_parity_branches() {
        assert_f64_bits(kappa_1_helper(0.0, 3.0, 0.0), 0.0, "kappa1 zero");
        assert_f64_bits(kappa_2_helper(0.0, 3.0, 0.0), 0.0, "kappa2 zero");
        assert_f64_bits(kappa_3_helper(0.0, 3, 0.0), 0.0, "kappa3 zero");
        assert_f64_bits(kappa_1_helper(2.0, 3.0, 0.0), 3.0, "kappa1");
        assert_f64_bits(kappa_2_helper(1.0, 3.0, 0.0), 2.0, "kappa2");
        assert_f64_bits(kappa_3_helper(1.0, 5, 0.0), 16.0, "kappa3 odd");
        assert_f64_bits(kappa_3_helper(1.0, 6, 0.0), 36.0, "kappa3 even");
    }

    #[test]
    fn kappa_and_phi_match_pinned_rdkit_across_boundary_and_ring_branches() {
        const CASES: [(&str, [u64; 5]); 6] = [
            ("", [0, 0, 0, 0, 0]),
            (
                "CCC",
                [
                    0,
                    0x4008_0000_0000_0000,
                    0x4000_0000_0000_0000,
                    0,
                    0x4000_0000_0000_0000,
                ],
            ),
            (
                "CCCC",
                [
                    0,
                    0x4010_0000_0000_0000,
                    0x4008_0000_0000_0000,
                    0x4000_0000_0000_0000,
                    0x4008_0000_0000_0000,
                ],
            ),
            (
                "C1CC1",
                [
                    0,
                    0x3ff5_5555_5555_5555,
                    0x3fcc_71c7_1c71_c71c,
                    0,
                    0x3fb9_48b0_fcd6_e9e0,
                ],
            ),
            (
                "c1ccccc1",
                [
                    0xbfe8_f5c2_8f5c_28f6,
                    0x400b_4ae5_ac96_d07c,
                    0x3ff9_b13b_4bc5_063b,
                    0x3fe2_a303_c5f0_83cd,
                    0x3fed_3790_5fe9_fff9,
                ],
            ),
            (
                "[Na+].[Cl-]",
                [
                    0x3ff4_a3d7_0a3d_70a4,
                    0x4024_bc52_e1e5_3581,
                    0x4002_51eb_851e_b852,
                    0x3fb0_b08a_703e_3bdb,
                    0x4027_be07_dc40_0b58,
                ],
            ),
        ];

        for (smiles, expected) in CASES {
            let molecule = Molecule::from_smiles(smiles).expect("Kappa fixture");
            let actual = [
                calc_hall_kier_alpha(&molecule, None),
                calc_kappa_1(&molecule),
                calc_kappa_2(&molecule).unwrap(),
                calc_kappa_3(&molecule).unwrap(),
                calc_phi(&molecule).unwrap(),
            ];
            for (index, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
                assert_eq!(
                    actual.to_bits(),
                    expected,
                    "{smiles:?} Hall-Kier/Kappa/Phi field {index}: actual={actual:?}"
                );
            }
        }
    }
}
