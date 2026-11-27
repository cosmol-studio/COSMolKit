//! Gasteiger charge summaries and charge-classified Labute surface areas.
//! RDKit MiscDescriptors.cpp: Copyright (C) 2026 Glysade LLC and other RDKit contributors.
//! MolSurf.cpp: Copyright (C) 2007-2023 Greg Landrum. See LICENSE for BSD terms.

use cosmolkit_core::{
    GasteigerChargeError, GasteigerChargeParams, ValenceAssignment, gasteiger_charges,
};
use cosmolkit_model::TopologyBlock;

use crate::{DescriptorComputedState, DescriptorError, DescriptorResult, labute_atom_rows};

const PEOE_VSA_CUTOFFS: [f64; 13] = [
    -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
];

#[derive(Clone, Copy)]
struct PartialChargeExtrema {
    minimum: f64,
    maximum: f64,
    minimum_absolute: f64,
    maximum_absolute: f64,
}

impl PartialChargeExtrema {
    fn calculate(
        topology: &TopologyBlock,
        valence: Option<&ValenceAssignment>,
    ) -> DescriptorResult<Self> {
        // RDKit❗✔️: PartialCharges calcPartialCharges(const ROMol &mol) {
        // RDKit❗✔️:   if(mol.getNumAtoms() == 0) {
        // RDKit❗✔️:     return {{}, -500, 500, 500, 500};
        // RDKit❗✔️:   }
        if topology.atoms.is_empty() {
            return Ok(Self {
                minimum: 500.0,
                maximum: -500.0,
                minimum_absolute: 500.0,
                maximum_absolute: 500.0,
            });
        }

        // RDKit❗✔️:   PartialCharges charges { {}, -500, 500, 500, 500 };
        // RDKit❗✔️:   std::vector<double> chgs(mol.getNumAtoms(), 0.0);
        // RDKit❗✔️:   computeGasteigerCharges(mol, chgs);
        // RDKit❗✔️:   for (auto chg : chgs) {
        // RDKit❗✔️:     if(!std::isnan(chg)) { // retain nans
        // RDKit❗✔️:       charges.minPartialCharge = std::min(charges.minPartialCharge, chg);
        // RDKit❗✔️:       charges.maxPartialCharge = std::max(charges.maxPartialCharge, chg);
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       charges.minPartialCharge = chg;
        // RDKit❗✔️:       charges.maxPartialCharge = chg;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        let rows = gasteiger_charges(topology, valence, &GasteigerChargeParams::default())
            .map_err(gasteiger_error)?;
        let mut minimum = 500.0;
        let mut maximum = -500.0;
        for row in &rows {
            let charge = row.atom_charge();
            if charge.is_nan() {
                minimum = charge;
                maximum = charge;
            } else {
                minimum = source_minimum(minimum, charge);
                maximum = source_maximum(maximum, charge);
            }
        }

        // RDKit❗✔️:   charges.charges.swap(chgs);
        // RDKit❗✔️:   double v1 = std::fabs(charges.minPartialCharge);
        // RDKit❗✔️:   double v2 = std::fabs(charges.maxPartialCharge);
        // RDKit❗✔️:   charges.maxAbsPartialCharge = std::max(v1,v2);
        // RDKit❗✔️:   charges.minAbsPartialCharge = std::min(v1,v2);
        // RDKit❗✔️:   return charges;
        let minimum_magnitude = minimum.abs();
        let maximum_magnitude = maximum.abs();
        Ok(Self {
            minimum,
            maximum,
            minimum_absolute: source_minimum(minimum_magnitude, maximum_magnitude),
            maximum_absolute: source_maximum(minimum_magnitude, maximum_magnitude),
        })
    }
}

/// Return the maximum graph-atom charge from the default permissive calculation.
pub fn max_partial_charge(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: double calcMaxPartialCharge(const ROMol &mol) {
    // RDKit❗✔️:   return calcPartialCharges(mol).maxPartialCharge;
    // RDKit❗✔️: }
    Ok(PartialChargeExtrema::calculate(topology, valence)?.maximum)
}

/// Return the minimum graph-atom charge from the default permissive calculation.
pub fn min_partial_charge(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: double calcMinPartialCharge(const ROMol &mol) {
    // RDKit❗✔️:   return calcPartialCharges(mol).minPartialCharge;
    // RDKit❗✔️: }
    Ok(PartialChargeExtrema::calculate(topology, valence)?.minimum)
}

/// Return the larger absolute value of the two signed charge extrema.
pub fn max_abs_partial_charge(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: double calcMaxAbsPartialCharge(const ROMol &mol) {
    // RDKit❗✔️:   return calcPartialCharges(mol).maxAbsPartialCharge;
    // RDKit❗✔️: }
    Ok(PartialChargeExtrema::calculate(topology, valence)?.maximum_absolute)
}

/// Return the smaller absolute value of the two signed charge extrema.
pub fn min_abs_partial_charge(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: double calcMinAbsPartialCharge(const ROMol &mol) {
    // RDKit❗✔️:   return calcPartialCharges(mol).minAbsPartialCharge;
    // RDKit❗✔️: }
    Ok(PartialChargeExtrema::calculate(topology, valence)?.minimum_absolute)
}

/// Calculate the 14 default PEOE-VSA bins in atom order.
///
/// Labute's existing single cache slot is read first, with hydrogen area
/// enabled on a miss and without forcing a prior slot to refresh.
pub fn peoe_vsa(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<Vec<f64>> {
    // RDKit source (MolSurf.cpp, calcPEOE_VSA default path):
    // RDKit❗✔️: std::vector<double> calcPEOE_VSA(const ROMol &mol, std::vector<double> *bins,
    // RDKit❗✔️:                                  bool force) {
    // RDKit❗✔️:   double blist[13] = {-.3, -.25, -.20, -.15, -.10, -.05, 0,
    // RDKit❗✔️:                       .05, .10,  .15,  .20,  .25,  .30};
    // RDKit❗✔️:   std::vector<double> res(lbins.size() + 1, 0);
    // RDKit❗✔️:   std::vector<double> vsaContribs(mol.getNumAtoms());
    // RDKit❗✔️:   getLabuteAtomContribs(mol, vsaContribs, tmp, true, force);
    // RDKit❗✔️:   std::vector<double> chgs(mol.getNumAtoms(), 0.0);
    // RDKit❗✔️:   computeGasteigerCharges(mol, chgs);
    // RDKit❗✔️:   assignContribsToBins(vsaContribs, chgs, lbins, res);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    let areas = labute_atom_rows(topology, state, true, false)?;
    let charges = gasteiger_charges(topology, valence, &GasteigerChargeParams::default())
        .map_err(gasteiger_error)?;

    let atom_charges: Vec<_> = charges.into_iter().map(|row| row.atom_charge()).collect();
    crate::assign_contribs_to_bins(&areas.atom_contributions, &atom_charges, &PEOE_VSA_CUTOFFS)
}

fn gasteiger_error(source: GasteigerChargeError) -> DescriptorError {
    DescriptorError::GasteigerCharge { source }
}

fn source_minimum(left: f64, right: f64) -> f64 {
    if right < left { right } else { left }
}

fn source_maximum(left: f64, right: f64) -> f64 {
    if left < right { right } else { left }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_charge_extrema_use_the_pinned_native_sentinels() {
        let extrema = PartialChargeExtrema::calculate(&TopologyBlock::default(), None).unwrap();
        assert_eq!(extrema.maximum, -500.0);
        assert_eq!(extrema.minimum, 500.0);
        assert_eq!(extrema.maximum_absolute, 500.0);
        assert_eq!(extrema.minimum_absolute, 500.0);
    }
}
