//! Gasteiger charge summaries and charge-classified Labute surface areas.
//! RDKit Descriptors.py charge extrema and MolSurf.cpp surface-area projections.
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
        // RDKit✔️✔️: chgs = rdPartialCharges.ComputeGasteigerCharges(mol)
        // RDKit✔️✔️: minChg = 500.
        // RDKit✔️✔️: maxChg = -500.
        // RDKit✔️✔️: for at in mol.GetAtoms():
        // RDKit✔️✔️:   chg = float(at.GetProp('_GasteigerCharge'))
        // RDKit✔️✔️:   minChg = min(chg, minChg)
        // RDKit✔️✔️:   maxChg = max(chg, maxChg)
        // RDKit✔️✔️: res = (minChg, maxChg)
        let rows = gasteiger_charges(topology, valence, &GasteigerChargeParams::default())
            .map_err(gasteiger_error)?;
        let mut minimum = 500.0;
        let mut maximum = -500.0;
        for row in &rows {
            let charge = row.atom_charge();
            // Python min/max retain their first argument when comparisons
            // are unordered. A later finite atom can therefore clear a NaN.
            minimum = source_minimum(charge, minimum);
            maximum = source_maximum(charge, maximum);
        }

        // RDKit✔️✔️: v1, v2 = _ChargeDescriptors(mol, force)
        // RDKit✔️✔️: return max(abs(v1), abs(v2))
        // RDKit✔️✔️: return min(abs(v1), abs(v2))
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
    // RDKit✔️✔️: _, res = _ChargeDescriptors(mol, force)
    // RDKit✔️✔️: return res
    Ok(PartialChargeExtrema::calculate(topology, valence)?.maximum)
}

/// Return the minimum graph-atom charge from the default permissive calculation.
pub fn min_partial_charge(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
) -> DescriptorResult<f64> {
    // RDKit✔️✔️: res, _ = _ChargeDescriptors(mol, force)
    // RDKit✔️✔️: return res
    Ok(PartialChargeExtrema::calculate(topology, valence)?.minimum)
}

/// Return the larger absolute value of the two signed charge extrema.
pub fn max_abs_partial_charge(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
) -> DescriptorResult<f64> {
    // RDKit✔️✔️: v1, v2 = _ChargeDescriptors(mol, force)
    // RDKit✔️✔️: return max(abs(v1), abs(v2))
    Ok(PartialChargeExtrema::calculate(topology, valence)?.maximum_absolute)
}

/// Return the smaller absolute value of the two signed charge extrema.
pub fn min_abs_partial_charge(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
) -> DescriptorResult<f64> {
    // RDKit✔️✔️: v1, v2 = _ChargeDescriptors(mol, force)
    // RDKit✔️✔️: return min(abs(v1), abs(v2))
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
    fn later_finite_atoms_clear_unordered_charge_extrema() {
        let fixture = crate::original_condition_fixture::Fixture::from_smiles("C[Fe].CC");
        let extrema =
            PartialChargeExtrema::calculate(&fixture.topology, Some(&fixture.valence)).unwrap();
        assert_eq!(
            extrema.minimum.to_bits(),
            (-0.06826238342059164_f64).to_bits()
        );
        assert_eq!(extrema.minimum.to_bits(), extrema.maximum.to_bits());
        assert_eq!(
            extrema.minimum_absolute.to_bits(),
            extrema.minimum.abs().to_bits()
        );
    }

    #[test]
    fn empty_charge_extrema_use_the_python_sentinels() {
        let extrema = PartialChargeExtrema::calculate(&TopologyBlock::default(), None).unwrap();
        assert_eq!(extrema.maximum, -500.0);
        assert_eq!(extrema.minimum, 500.0);
        assert_eq!(extrema.maximum_absolute, 500.0);
        assert_eq!(extrema.minimum_absolute, 500.0);
    }
}
