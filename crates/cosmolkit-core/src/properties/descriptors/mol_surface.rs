//! Molecular-surface descriptor implementation boundary.
//!
//! This module exclusively owns Labute atom contributions, Labute ASA, shared
//! VSA bin assignment, SlogP-VSA, and SMR-VSA. Crippen atom contributions and
//! periodic-table data remain single shared implementations.

use super::{DescriptorError, DescriptorResult};
use crate::{BondOrder, Molecule, model::molecule::LabuteDescriptorCache};
use std::sync::Arc;

const RDKIT_SLOGP_VSA_BINS: [f64; 11] = [-0.4, -0.2, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6];
const RDKIT_SMR_VSA_BINS: [f64; 9] = [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0];

pub(super) fn assign_contributions_to_bins(
    contributions: &[f64],
    properties: &[f64],
    bins: &[f64],
) -> DescriptorResult<Vec<f64>> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::assignContribsToBins
    // RDKit✔️✔️: void assignContribsToBins(const std::vector<double> &contribs,
    // RDKit✔️✔️:                           const std::vector<double> &binProp,
    // RDKit✔️✔️:                           const std::vector<double> &bins,
    // RDKit✔️✔️:                           std::vector<double> &res) {
    // RDKit✔️✔️:   PRECONDITION(contribs.size() == binProp.size(), "mismatched array sizes");
    if contributions.len() != properties.len() {
        return Err(DescriptorError::MismatchedArraySizes {
            function: "assign_contributions_to_bins",
            contributions: contributions.len(),
            properties: properties.len(),
        });
    }
    // RDKit✔️✔️:   PRECONDITION(res.size() >= bins.size() + 1, "mismatched array sizes");
    let mut result = vec![0.0; bins.len() + 1];
    // RDKit✔️✔️:   for (unsigned int i = 0; i < contribs.size(); ++i) {
    for index in 0..contributions.len() {
        // RDKit✔️✔️:     double cVal = contribs[i];
        let contribution = contributions[index];
        // RDKit✔️✔️:     double bVal = binProp[i];
        let property = properties[index];
        // RDKit✔️✔️:     unsigned int idx =
        // RDKit✔️✔️:         std::upper_bound(bins.begin(), bins.end(), bVal) - bins.begin();
        let mut first = 0_usize;
        let mut count = bins.len();
        while count > 0 {
            let step = count / 2;
            let middle = first + step;
            if !(property < bins[middle]) {
                first = middle + 1;
                count -= step + 1;
            } else {
                count = step;
            }
        }
        // RDKit✔️✔️:     res[idx] += cVal;
        result[first] += contribution;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::assignContribsToBins
    Ok(result)
}

pub(super) fn calc_slogp_vsa(molecule: &Molecule, bins: Option<&[f64]>, force: bool) -> DescriptorResult<Vec<f64>> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcSlogP_VSA
    // RDKit✔️✔️: std::vector<double> calcSlogP_VSA(const ROMol &mol, std::vector<double> *bins,
    // RDKit✔️✔️:                                   bool force) {
    // RDKit✔️✔️:   // FIX: use force value to include caching
    // RDKit✔️✔️:   std::vector<double> lbins;
    // RDKit✔️✔️:   if (!bins) {
    // RDKit✔️✔️:     double blist[11] = {-0.4, -0.2, 0,   0.1, 0.15, 0.2,
    // RDKit✔️✔️:                         0.25, 0.3,  0.4, 0.5, 0.6};
    // RDKit✔️✔️:     lbins.resize(11);
    // RDKit✔️✔️:     std::copy(blist, blist + 11, lbins.begin());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     lbins.resize(bins->size());
    // RDKit✔️✔️:     std::copy(bins->begin(), bins->end(), lbins.begin());
    // RDKit✔️✔️:   }
    let local_bins = bins.unwrap_or(&RDKIT_SLOGP_VSA_BINS).to_vec();
    // RDKit✔️✔️:   std::vector<double> res(lbins.size() + 1, 0);

    // RDKit✔️✔️:   std::vector<double> vsaContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   double tmp;
    // RDKit✔️✔️:   getLabuteAtomContribs(mol, vsaContribs, tmp, true, force);
    let vsa_contributions = get_labute_atom_contribs(molecule, true, force);
    // RDKit✔️✔️:   std::vector<double> logpContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   std::vector<double> mrContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   getCrippenAtomContribs(mol, logpContribs, mrContribs, force);
    let crippen_contributions = super::rdkit_crippen_atom_contribs(molecule, force)?;

    // RDKit✔️✔️:   assignContribsToBins(vsaContribs, logpContribs, lbins, res);
    let result = assign_contributions_to_bins(
        &vsa_contributions.atom_contributions,
        &crippen_contributions.logp,
        &local_bins,
    )?;

    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcSlogP_VSA
    Ok(result)
}

pub(super) fn calc_smr_vsa(molecule: &Molecule, bins: Option<&[f64]>, force: bool) -> DescriptorResult<Vec<f64>> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcSMR_VSA
    // RDKit✔️✔️: std::vector<double> calcSMR_VSA(const ROMol &mol, std::vector<double> *bins,
    // RDKit✔️✔️:                                 bool force) {
    // RDKit✔️✔️:   std::vector<double> lbins;
    // RDKit✔️✔️:   if (!bins) {
    // RDKit✔️✔️:     double blist[9] = {1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0};
    // RDKit✔️✔️:     lbins.resize(9);
    // RDKit✔️✔️:     std::copy(blist, blist + 9, lbins.begin());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     lbins.resize(bins->size());
    // RDKit✔️✔️:     std::copy(bins->begin(), bins->end(), lbins.begin());
    // RDKit✔️✔️:   }
    let local_bins = bins.unwrap_or(&RDKIT_SMR_VSA_BINS).to_vec();
    // RDKit✔️✔️:   std::vector<double> res(lbins.size() + 1, 0);

    // RDKit✔️✔️:   std::vector<double> vsaContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   double tmp;
    // RDKit✔️✔️:   getLabuteAtomContribs(mol, vsaContribs, tmp, true, force);
    let vsa_contributions = get_labute_atom_contribs(molecule, true, force);
    // RDKit✔️✔️:   std::vector<double> logpContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   std::vector<double> mrContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   getCrippenAtomContribs(mol, logpContribs, mrContribs, force);
    let crippen_contributions = super::rdkit_crippen_atom_contribs(molecule, force)?;

    // RDKit✔️✔️:   assignContribsToBins(vsaContribs, mrContribs, lbins, res);
    let result = assign_contributions_to_bins(
        &vsa_contributions.atom_contributions,
        &crippen_contributions.mr,
        &local_bins,
    )?;

    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcSMR_VSA
    Ok(result)
}

pub(super) fn slogp_vsa_bin(molecule: &Molecule, index: usize) -> DescriptorResult<f64> {
    // BEGIN RDKIT PYTHON FUNCTION: rdkit.Chem.MolSurf._InstallDescriptors SlogP projection
    // RDKit✔️✔️:   for i in range(len(logpBins)):
    // RDKit✔️✔️:     fn = lambda x, y=i: SlogP_VSA_(x, force=0)[y]
    // RDKit✔️✔️:   i += 1
    // RDKit✔️✔️:   fn = lambda x, y=i: SlogP_VSA_(x, force=0)[y]
    // END RDKIT PYTHON FUNCTION: rdkit.Chem.MolSurf._InstallDescriptors SlogP projection
    Ok(calc_slogp_vsa(molecule, None, false)?[index])
}

pub(super) fn smr_vsa_bin(molecule: &Molecule, index: usize) -> DescriptorResult<f64> {
    // BEGIN RDKIT PYTHON FUNCTION: rdkit.Chem.MolSurf._InstallDescriptors SMR projection
    // RDKit✔️✔️:   for i in range(len(mrBins)):
    // RDKit✔️✔️:     fn = lambda x, y=i: SMR_VSA_(x, force=0)[y]
    // RDKit✔️✔️:   i += 1
    // RDKit✔️✔️:   fn = lambda x, y=i: SMR_VSA_(x, force=0)[y]
    // END RDKIT PYTHON FUNCTION: rdkit.Chem.MolSurf._InstallDescriptors SMR projection
    Ok(calc_smr_vsa(molecule, None, false)?[index])
}

pub(super) fn get_labute_atom_contribs(
    molecule: &Molecule,
    include_hydrogens: bool,
    force: bool,
) -> LabuteDescriptorCache {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::getLabuteAtomContribs
    // RDKit✔️✔️: double getLabuteAtomContribs(const ROMol &mol, std::vector<double> &Vi,
    // RDKit✔️✔️:                              double &hContrib, bool includeHs, bool force) {
    // RDKit✔️✔️:   TEST_ASSERT(Vi.size() == mol.getNumAtoms());
    // RDKit✔️✔️:   if (!force && mol.hasProp(common_properties::_labuteAtomContribs)) {
    if !force && let Some(cached) = molecule.labute_descriptor_cache() {
        // RDKit✔️✔️:     mol.getProp(common_properties::_labuteAtomContribs, Vi);
        // RDKit✔️✔️:     mol.getProp(common_properties::_labuteAtomHContrib, hContrib);
        // RDKit✔️✔️:     double res;
        // RDKit✔️✔️:     mol.getProp(common_properties::_labuteASA, res);
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        return cached;
    }
    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let atom_count = molecule.atoms().len();
    // RDKit✔️✔️:   std::vector<double> rads(nAtoms);
    let mut radii = vec![0.0; atom_count];
    let mut contributions = vec![0.0; atom_count];
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    for atom_index in 0..atom_count {
        // RDKit✔️✔️:     rads[i] = PeriodicTable::getTable()->getRb0(
        // RDKit✔️✔️:         mol.getAtomWithIdx(i)->getAtomicNum());
        radii[atom_index] = crate::chemistry::valence::rdkit_rb0(molecule.atoms()[atom_index].atomic_number());
        // RDKit✔️✔️:     Vi[i] = 0.0;
        contributions[atom_index] = 0.0;
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   for (ROMol::ConstBondIterator bondIt = mol.beginBonds();
    // RDKit✔️✔️:        bondIt != mol.endBonds(); ++bondIt) {
    for bond in molecule.bonds() {
        // RDKit✔️✔️:     const double bondScaleFacts[4] = {.1, 0, .2, .3};
        const BOND_SCALE_FACTORS: [f64; 4] = [0.1, 0.0, 0.2, 0.3];
        // RDKit✔️✔️:     double Ri = rads[(*bondIt)->getBeginAtomIdx()];
        let begin_radius = radii[bond.begin().index()];
        // RDKit✔️✔️:     double Rj = rads[(*bondIt)->getEndAtomIdx()];
        let end_radius = radii[bond.end().index()];
        // RDKit✔️✔️:     double bij = Ri + Rj;
        let mut bond_distance = begin_radius + end_radius;
        // RDKit✔️✔️:     if (!(*bondIt)->getIsAromatic()) {
        if !bond.is_aromatic() {
            // RDKit✔️✔️:       if ((*bondIt)->getBondType() < 4) {
            // RDKit✔️✔️:         bij -= bondScaleFacts[(*bondIt)->getBondType()];
            // RDKit✔️✔️:       }
            let bond_type = match bond.order() {
                BondOrder::Null | BondOrder::Unspecified => Some(0),
                BondOrder::Single => Some(1),
                BondOrder::Double => Some(2),
                BondOrder::Triple => Some(3),
                _ => None,
            };
            if let Some(bond_type) = bond_type {
                bond_distance -= BOND_SCALE_FACTORS[bond_type];
            }
            // RDKit✔️✔️:     } else {
        } else {
            // RDKit✔️✔️:       bij -= bondScaleFacts[0];
            bond_distance -= BOND_SCALE_FACTORS[0];
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     double dij = std::min(std::max(fabs(Ri - Rj), bij), Ri + Rj);
        let distance = (begin_radius - end_radius)
            .abs()
            .max(bond_distance)
            .min(begin_radius + end_radius);
        // RDKit✔️✔️:     Vi[(*bondIt)->getBeginAtomIdx()] += Rj * Rj - (Ri - dij) * (Ri - dij) / dij;
        contributions[bond.begin().index()] +=
            end_radius * end_radius - (begin_radius - distance) * (begin_radius - distance) / distance;
        // RDKit✔️✔️:     Vi[(*bondIt)->getEndAtomIdx()] += Ri * Ri - (Rj - dij) * (Rj - dij) / dij;
        contributions[bond.end().index()] +=
            begin_radius * begin_radius - (end_radius - distance) * (end_radius - distance) / distance;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   hContrib = 0.0;
    let mut hydrogen_contribution = 0.0;
    // RDKit✔️✔️:   if (includeHs) {
    if include_hydrogens {
        // RDKit✔️✔️:     double Rj = PeriodicTable::getTable()->getRb0(1);
        let hydrogen_radius = crate::chemistry::valence::rdkit_rb0(1);
        // RDKit✔️✔️:     for (unsigned int i = 0; i < nAtoms; ++i) {
        for atom_index in 0..atom_count {
            // RDKit✔️✔️:       double Ri = rads[i];
            let atom_radius = radii[atom_index];
            // RDKit✔️✔️:       double bij = Ri + Rj;
            let bond_distance = atom_radius + hydrogen_radius;
            // RDKit✔️✔️:       double dij = std::min(std::max(fabs(Ri - Rj), bij), Ri + Rj);
            let distance = (atom_radius - hydrogen_radius)
                .abs()
                .max(bond_distance)
                .min(atom_radius + hydrogen_radius);
            // RDKit✔️✔️:       Vi[i] += Rj * Rj - (Ri - dij) * (Ri - dij) / dij;
            contributions[atom_index] +=
                hydrogen_radius * hydrogen_radius - (atom_radius - distance) * (atom_radius - distance) / distance;
            // RDKit✔️✔️:       hContrib += Ri * Ri - (Rj - dij) * (Rj - dij) / dij;
            hydrogen_contribution +=
                atom_radius * atom_radius - (hydrogen_radius - distance) * (hydrogen_radius - distance) / distance;
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   double res = 0.0;
    let mut asa = 0.0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    for atom_index in 0..atom_count {
        // RDKit✔️✔️:     double Ri = rads[i];
        let atom_radius = radii[atom_index];
        // RDKit✔️✔️:     Vi[i] = M_PI * Ri * (4. * Ri - Vi[i]);
        contributions[atom_index] =
            std::f64::consts::PI * atom_radius * (4.0 * atom_radius - contributions[atom_index]);
        // RDKit✔️✔️:     res += Vi[i];
        asa += contributions[atom_index];
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   if (includeHs && fabs(hContrib) > 1e-4) {
    if include_hydrogens && hydrogen_contribution.abs() > 1e-4 {
        // RDKit✔️✔️:     double Rj = PeriodicTable::getTable()->getRb0(1);
        let hydrogen_radius = crate::chemistry::valence::rdkit_rb0(1);
        // RDKit✔️✔️:     hContrib = M_PI * Rj * (4. * Rj - hContrib);
        hydrogen_contribution =
            std::f64::consts::PI * hydrogen_radius * (4.0 * hydrogen_radius - hydrogen_contribution);
        // RDKit✔️✔️:     res += hContrib;
        asa += hydrogen_contribution;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   mol.setProp(common_properties::_labuteAtomContribs, Vi, true);
    // RDKit✔️✔️:   mol.setProp(common_properties::_labuteAtomHContrib, hContrib, true);
    // RDKit✔️✔️:   mol.setProp(common_properties::_labuteASA, res, true);
    let computed = LabuteDescriptorCache {
        atom_contributions: Arc::from(contributions),
        hydrogen_contribution,
        asa,
    };
    molecule.set_labute_descriptor_cache(computed.clone());

    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::getLabuteAtomContribs
    computed
}

pub(super) fn calc_labute_asa(molecule: &Molecule, include_hydrogens: bool, force: bool) -> f64 {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcLabuteASA
    // RDKit✔️✔️: double calcLabuteASA(const ROMol &mol, bool includeHs, bool force) {
    // RDKit✔️✔️:   if (!force && mol.hasProp(common_properties::_labuteASA)) {
    if !force && let Some(cached) = molecule.labute_descriptor_cache() {
        // RDKit✔️✔️:     double res;
        // RDKit✔️✔️:     mol.getProp(common_properties::_labuteASA, res);
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        return cached.asa;
    }
    // RDKit✔️✔️:   std::vector<double> contribs;
    // RDKit✔️✔️:   contribs.resize(mol.getNumAtoms());
    // RDKit✔️✔️:   double hContrib;
    // RDKit✔️✔️:   double res;
    // RDKit✔️✔️:   res = getLabuteAtomContribs(mol, contribs, hContrib, includeHs, force);
    let result = get_labute_atom_contribs(molecule, include_hydrogens, force).asa;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcLabuteASA
    result
}

pub(super) fn labute_helper(molecule: &Molecule, include_hydrogens: bool, force: bool) -> Vec<f64> {
    // BEGIN RDKIT PYTHON FUNCTION: rdkit.Chem.MolSurf._LabuteHelper
    // RDKit✔️✔️: def _LabuteHelper(mol, includeHs=1, force=0):
    // RDKit✔️✔️:   """ *Internal Use Only*
    // RDKit✔️✔️:     helper function for LabuteASA calculation
    // RDKit✔️✔️:     returns an array of atomic contributions to the ASA
    // RDKit✔️✔️:
    // RDKit✔️✔️:   **Note:** Changes here affect the version numbers of all ASA descriptors
    // RDKit✔️✔️:
    // RDKit✔️✔️:   """
    // RDKit✔️✔️:   if not force:
    // RDKit✔️✔️:     try:
    // RDKit✔️✔️:       res = mol._labuteContribs
    // RDKit✔️✔️:     except AttributeError:
    // RDKit✔️✔️:       pass
    // RDKit✔️✔️:     else:
    // RDKit✔️✔️:       if res:
    // RDKit✔️✔️:         return res
    // RDKit✔️✔️:   tpl = rdMolDescriptors._CalcLabuteASAContribs(mol, includeHs)
    let contributions = get_labute_atom_contribs(molecule, include_hydrogens, force);
    // RDKit✔️✔️:   ats, hs = tpl
    // RDKit✔️✔️:   Vi = [hs] + list(ats)
    let mut result = Vec::with_capacity(contributions.atom_contributions.len() + 1);
    result.push(contributions.hydrogen_contribution);
    result.extend_from_slice(&contributions.atom_contributions);
    // RDKit✔️✔️:   mol._labuteContribs = Vi
    // RDKit✔️✔️:   return Vi
    // END RDKIT PYTHON FUNCTION: rdkit.Chem.MolSurf._LabuteHelper
    result
}

#[cfg(test)]
mod tests {
    use super::{
        assign_contributions_to_bins, calc_labute_asa, calc_slogp_vsa, calc_smr_vsa, get_labute_atom_contribs,
        labute_helper, slogp_vsa_bin, smr_vsa_bin,
    };
    use crate::{Molecule, properties::descriptors::DescriptorError};

    fn bits(values: &[f64]) -> Vec<u64> {
        values.iter().map(|value| value.to_bits()).collect()
    }

    #[test]
    fn labute_contributions_match_pinned_rdkit_bits_across_bond_and_hydrogen_branches() {
        const CASES: [(&str, &[u64]); 7] = [
            ("", &[0x0000_0000_0000_0000, 0x0000_0000_0000_0000]),
            (
                "CC",
                &[
                    0x3ff4_1b85_1ccc_a042,
                    0x401b_b1e8_2a1b_1454,
                    0x401b_b1e8_2a1b_1454,
                    0x402e_3558_cdb4_a85c,
                ],
            ),
            (
                "C=C",
                &[
                    0x3ff4_1b85_1ccc_a042,
                    0x401a_50d4_840e_2bfb,
                    0x401a_50d4_840e_2bfb,
                    0x402c_d445_27a7_c003,
                ],
            ),
            (
                "C#N",
                &[
                    0x3ff4_c3cc_4ecd_cb56,
                    0x401a_49de_40b5_ed92,
                    0x4015_0c2d_4cba_cd34,
                    0x402a_437f_5092_16ce,
                ],
            ),
            (
                "c1ccccc1",
                &[
                    0x3ff0_87fd_7763_0132,
                    0x4018_43f5_ba92_4bc9,
                    0x4018_43f5_ba92_4bc9,
                    0x4018_43f5_ba92_4bc9,
                    0x4018_43f5_ba92_4bc9,
                    0x4018_43f5_ba92_4bc9,
                    0x4018_43f5_ba92_4bc9,
                    0x4042_b738_37a8_d0e1,
                ],
            ),
            (
                "[H][H]",
                &[
                    0x3ff7_c1bb_ef62_32e4,
                    0x3ff7_c1bb_ef62_32e4,
                    0x3ff7_c1bb_ef62_32e4,
                    0x4011_d14c_f389_a62b,
                ],
            ),
            (
                "CCO",
                &[
                    0x3ff4_2e34_49f8_93d3,
                    0x401b_b1e8_2a1b_1454,
                    0x401a_6d72_7738_75fb,
                    0x4014_6d15_8473_e02a,
                    0x4033_e5ff_4e11_63db,
                ],
            ),
        ];

        for (smiles, expected) in CASES {
            let molecule = Molecule::from_smiles(smiles).expect("Labute fixture");
            let contribution = get_labute_atom_contribs(&molecule, true, true);
            let mut actual = Vec::with_capacity(contribution.atom_contributions.len() + 2);
            actual.push(contribution.hydrogen_contribution);
            actual.extend_from_slice(&contribution.atom_contributions);
            actual.push(contribution.asa);
            assert_eq!(bits(&actual), expected, "{smiles:?} Labute contributions");
            assert_eq!(
                bits(&labute_helper(&molecule, true, false)),
                expected[..expected.len() - 1],
                "{smiles:?} hydrogen-first helper"
            );
            assert_eq!(
                calc_labute_asa(&molecule, true, false).to_bits(),
                *expected.last().unwrap(),
                "{smiles:?} Labute ASA"
            );
        }
    }

    #[test]
    fn labute_cache_preserves_source_call_order_force_clone_and_invalidation_semantics() {
        const WITHOUT_HYDROGENS: u64 = 0x402b_ca6e_1564_c404;
        const WITH_HYDROGENS: u64 = 0x402e_3558_cdb4_a85c;

        let mut molecule = Molecule::from_smiles("CC").expect("Labute cache fixture");
        assert!(molecule.labute_descriptor_cache().is_none());
        assert_eq!(calc_labute_asa(&molecule, false, false).to_bits(), WITHOUT_HYDROGENS);
        assert_eq!(calc_labute_asa(&molecule, true, false).to_bits(), WITHOUT_HYDROGENS);

        let clone = molecule.clone();
        assert_eq!(calc_labute_asa(&clone, true, true).to_bits(), WITH_HYDROGENS);
        assert_eq!(calc_labute_asa(&molecule, true, false).to_bits(), WITHOUT_HYDROGENS);

        assert_eq!(calc_labute_asa(&molecule, true, true).to_bits(), WITH_HYDROGENS);
        assert_eq!(calc_labute_asa(&molecule, false, false).to_bits(), WITH_HYDROGENS);

        let topology = molecule.topology_block().clone();
        molecule.replace_topology_block(topology);
        assert!(molecule.labute_descriptor_cache().is_none());
        assert_eq!(calc_labute_asa(&molecule, false, false).to_bits(), WITHOUT_HYDROGENS);
    }

    #[test]
    fn shared_vsa_binning_matches_source_boundaries_nonfinite_values_and_order() {
        assert_eq!(
            assign_contributions_to_bins(
                &[1.0, 2.0, 4.0, 8.0, 16.0],
                &[-1.0, 0.0, 1.0, 2.0, 3.0],
                &[0.0, 1.0, 2.0],
            ),
            Ok(vec![1.0, 2.0, 4.0, 24.0])
        );
        assert_eq!(
            assign_contributions_to_bins(&[2.0], &[1.0], &[1.0, 1.0, 1.0]),
            Ok(vec![0.0, 0.0, 0.0, 2.0])
        );
        assert_eq!(
            assign_contributions_to_bins(&[1.0, 2.0, 3.0], &[0.0, 1.0, 2.0], &[]),
            Ok(vec![6.0])
        );
        assert_eq!(
            assign_contributions_to_bins(&[1.0, 2.0, 4.0], &[f64::NEG_INFINITY, f64::INFINITY, f64::NAN], &[0.0],),
            Ok(vec![1.0, 6.0])
        );
        assert_eq!(
            assign_contributions_to_bins(
                &[1.0, 2.0, 4.0, 8.0, 16.0],
                &[-1.0, 0.0, 0.5, 1.0, 2.0],
                &[2.0, 0.0, 1.0],
            ),
            Ok(vec![1.0, 0.0, 6.0, 24.0])
        );

        let ordered = assign_contributions_to_bins(&[1.0e16, 1.0, -1.0e16], &[0.0, 0.0, 0.0], &[]).unwrap();
        assert_eq!(ordered[0].to_bits(), 0.0_f64.to_bits());

        assert_eq!(
            assign_contributions_to_bins(&[1.0, 2.0], &[0.0], &[0.0]),
            Err(DescriptorError::MismatchedArraySizes {
                function: "assign_contributions_to_bins",
                contributions: 2,
                properties: 1,
            })
        );
    }

    #[test]
    fn slogp_and_smr_vsa_match_pinned_rdkit_default_and_custom_vectors() {
        const DEFAULT_CASES: [(&str, &[u64], &[u64]); 3] = [
            ("", &[0; 12], &[0; 10]),
            (
                "CCO",
                &[
                    0,
                    0x4027_6d43_fdd6_2b12,
                    0,
                    0,
                    0x401b_b1e8_2a1b_1454,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                ],
                &[
                    0x4014_6d15_8473_e02a,
                    0,
                    0,
                    0,
                    0x401b_b1e8_2a1b_1454,
                    0x401a_6d72_7738_75fb,
                    0,
                    0,
                    0,
                    0,
                ],
            ),
            (
                "CC(=O)N",
                &[
                    0x4016_ef46_86f2_339f,
                    0x4017_a0f3_b914_a2a9,
                    0x4013_2d9b_27d4_2d79,
                    0,
                    0x401b_b1e8_2a1b_1454,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                ],
                &[
                    0x4013_2d9b_27d4_2d79,
                    0,
                    0,
                    0x4016_ef46_86f2_339f,
                    0x401b_b1e8_2a1b_1454,
                    0,
                    0,
                    0,
                    0,
                    0x4017_a0f3_b914_a2a9,
                ],
            ),
        ];

        for (smiles, expected_slogp, expected_smr) in DEFAULT_CASES {
            let molecule = Molecule::from_smiles(smiles).expect("VSA fixture");
            let slogp = calc_slogp_vsa(&molecule, None, true).unwrap();
            let smr = calc_smr_vsa(&molecule, None, true).unwrap();
            assert_eq!(bits(&slogp), expected_slogp, "{smiles:?} SlogP-VSA");
            assert_eq!(bits(&smr), expected_smr, "{smiles:?} SMR-VSA");
            for (index, expected) in expected_slogp.iter().copied().enumerate() {
                assert_eq!(
                    slogp_vsa_bin(&molecule, index).unwrap().to_bits(),
                    expected,
                    "{smiles:?} SlogP-VSA scalar {}",
                    index + 1
                );
            }
            for (index, expected) in expected_smr.iter().copied().enumerate() {
                assert_eq!(
                    smr_vsa_bin(&molecule, index).unwrap().to_bits(),
                    expected,
                    "{smiles:?} SMR-VSA scalar {}",
                    index + 1
                );
            }
        }

        const CUSTOM_CASES: [(&[f64], &[u64], &[u64]); 3] = [
            (
                &[-0.2, 0.0, 0.2],
                &[0x4027_6d43_fdd6_2b12, 0, 0x401b_b1e8_2a1b_1454, 0],
                &[0, 0, 0, 0x4032_a31c_0971_da9e],
            ),
            (
                &[0.0, 0.0],
                &[0x4027_6d43_fdd6_2b12, 0, 0x401b_b1e8_2a1b_1454],
                &[0, 0, 0x4032_a31c_0971_da9e],
            ),
            (
                &[2.0, 0.0, 1.0],
                &[0x4027_6d43_fdd6_2b12, 0, 0x401b_b1e8_2a1b_1454, 0],
                &[0, 0, 0x4014_6d15_8473_e02a, 0x402b_0fad_50a9_c528],
            ),
        ];
        for (bins, expected_slogp, expected_smr) in CUSTOM_CASES {
            let molecule = Molecule::from_smiles("CCO").expect("custom VSA fixture");
            assert_eq!(
                bits(&calc_slogp_vsa(&molecule, Some(bins), true).unwrap()),
                expected_slogp
            );
            assert_eq!(bits(&calc_smr_vsa(&molecule, Some(bins), true).unwrap()), expected_smr);
        }
    }

    #[test]
    fn vsa_force_controls_shared_labute_cache_without_changing_vector_shape() {
        let molecule = Molecule::from_smiles("CC").expect("VSA force fixture");
        let _ = get_labute_atom_contribs(&molecule, false, false);

        let cached_slogp = calc_slogp_vsa(&molecule, None, false).unwrap();
        assert_eq!(cached_slogp.len(), 12);
        assert_eq!(cached_slogp[4].to_bits(), 0x402b_ca6e_1564_c404);
        let forced_slogp = calc_slogp_vsa(&molecule, None, true).unwrap();
        assert_eq!(forced_slogp.len(), 12);
        assert_eq!(forced_slogp[4].to_bits(), 0x402b_b1e8_2a1b_1454);

        let molecule_smr = Molecule::from_smiles("CC").expect("SMR-VSA force fixture");
        let _ = get_labute_atom_contribs(&molecule_smr, false, false);
        let cached_smr = calc_smr_vsa(&molecule_smr, None, false).unwrap();
        assert_eq!(cached_smr.len(), 10);
        assert_eq!(cached_smr[4].to_bits(), 0x402b_ca6e_1564_c404);
        let forced_smr = calc_smr_vsa(&molecule_smr, None, true).unwrap();
        assert_eq!(forced_smr.len(), 10);
        assert_eq!(forced_smr[4].to_bits(), 0x402b_b1e8_2a1b_1454);
    }
}
