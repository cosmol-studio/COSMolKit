//! Working-molecule preparation for the eight standard BCUT2D values.
//! RDKit 4ce8a0636f9b8ef28d6c3f4af4f502392c997d6b, GraphMol/Descriptors/BCUT.cpp
//! (Copyright 2020 Brian P. Kelley) and GraphMol/AddHs.cpp (Copyright 2003–2025
//! Greg Landrum and other RDKit contributors); BSD terms are in LICENSE.

mod burden;

use cosmolkit_core::{GasteigerChargeParams, RemoveHsParams};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

use crate::{DescriptorComputedState, DescriptorError, DescriptorInput, DescriptorResult};

fn remove_all_hs_params() -> RemoveHsParams {
    // RDKit❗✔️: void removeAllHs(RWMol &mol, bool sanitize) {
    // RDKit❗✔️:   RemoveHsParameters ps;
    // RDKit❗✔️:   ps.removeDegreeZero = true;
    // RDKit❗✔️:   ps.removeHigherDegrees = true;
    // RDKit❗✔️:   ps.removeOnlyHNeighbors = true;
    // RDKit❗✔️:   ps.removeIsotopes = true;
    // RDKit❗✔️:   ps.removeDummyNeighbors = true;
    // RDKit❗✔️:   ps.removeDefiningBondStereo = true;
    // RDKit❗✔️:   ps.removeWithWedgedBond = true;
    // RDKit❗✔️:   ps.removeWithQuery = true;
    // RDKit❗✔️:   ps.removeNonimplicit = true;
    // RDKit❗✔️:   ps.removeInSGroups = true;
    // RDKit❗✔️:   ps.showWarnings = false;
    // RDKit❗✔️:   ps.removeHydrides = true;
    // RDKit❗✔️:   ps.removeNontetrahedralNeighbors = true;
    // RDKit❗✔️:   removeHs(mol, ps, sanitize);
    // RDKit❗✔️: };
    RemoveHsParams {
        remove_degree_zero: true,
        remove_higher_degrees: true,
        remove_only_h_neighbors: true,
        remove_isotopes: true,
        remove_dummy_neighbors: true,
        remove_defining_bond_stereo: true,
        remove_with_wedged_bond: true,
        remove_with_query: true,
        remove_nonimplicit: true,
        remove_in_sgroups: true,
        show_warnings: false,
        remove_hydrides: true,
        remove_nontetrahedral_neighbors: true,
        ..Default::default()
    }
}

/// Standard BCUT2D on a private hydrogen-removed, sanitized graph.
pub fn bcut2d(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> DescriptorResult<[f64; 8]> {
    // RDKit❗✔️: std::vector<double> BCUT2D(const ROMol &m) {
    // RDKit❗✔️:   if (!m.getNumAtoms()) {
    // RDKit❗✔️:     return std::vector<double>(NUM_BCUTS, 0.0);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::unique_ptr<ROMol> mol(MolOps::removeAllHs(m));
    // RDKit❗✔️:   unsigned int numAtoms = mol->getNumAtoms();
    if topology.atoms.is_empty() {
        return Ok([0.0; 8]);
    }
    // RDKit❗✔️: ROMol *removeAllHs(const ROMol &mol, bool sanitize) {
    // RDKit❗✔️:   auto *res = new RWMol(mol);
    // RDKit❗✔️:   try {
    // RDKit❗✔️:     removeAllHs(*res, sanitize);
    // RDKit❗✔️:   } catch (const MolSanitizeException &) {
    // RDKit❗✔️:     delete res;
    // RDKit❗✔️:     throw;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return static_cast<ROMol *>(res);
    // RDKit❗✔️: }
    let working = cosmolkit_core::remove_hydrogens_with_params(
        topology.clone(),
        coordinates.clone(),
        properties.clone(),
        &remove_all_hs_params(),
    )
    .map_err(|source| DescriptorError::Hydrogens {
        function: "bcut2d",
        source,
    })?;
    let topology = &working.topology;
    if topology.atoms.is_empty() {
        return Ok([0.0; 8]);
    }

    // RDKit❗✔️:   std::vector<double> masses;
    // RDKit❗✔️:   std::vector<double> charges;
    // RDKit❗✔️:   masses.reserve(numAtoms);
    // RDKit❗✔️:   charges.reserve(numAtoms);
    // RDKit❗✔️:
    // RDKit❗✔️:   RDKit::computeGasteigerCharges(*mol, 12, true);
    // RDKit❗✔️:   for (auto &atom : mol->atoms()) {
    // RDKit❗✔️:     masses.push_back(atom->getMass());
    // RDKit❗✔️:     charges.push_back(
    // RDKit❗✔️:         atom->getProp<double>(common_properties::_GasteigerCharge));
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    let mut params = GasteigerChargeParams::default();
    params.set_throw_on_param_failure(true);
    let charges =
        cosmolkit_core::gasteiger_charges(topology, working.final_valence.as_ref(), &params)
            .map_err(|source| DescriptorError::GasteigerCharge { source })?;
    let masses: Vec<_> = topology
        .atoms
        .iter()
        .map(cosmolkit_core::atom_mass)
        .collect();
    let own_charges: Vec<_> = charges.iter().map(|row| row.atom_charge()).collect();

    // RDKit❗✔️:   std::vector<double> slogp(numAtoms, 0.0);
    // RDKit❗✔️:   std::vector<double> cmr(numAtoms, 0.0);
    // RDKit❗✔️:   getCrippenAtomContribs(*mol, slogp, cmr);
    // RDKit❗✔️:   // polarizability? - need model
    // RDKit❗✔️:   // slogp?  sasa?
    let missing = |field| DescriptorError::Unsupported {
        function: "bcut2d",
        detail: format!("sanitized hydrogen removal did not supply {field}"),
    };
    let valence = working
        .final_valence
        .as_ref()
        .ok_or_else(|| missing("valence"))?;
    let rings = working
        .final_rings
        .as_ref()
        .ok_or_else(|| missing("rings"))?;
    let input = DescriptorInput::new(
        topology,
        &working.coordinates,
        &working.properties,
        valence,
        rings,
    );
    let crippen = crate::crippen_contributions(
        &input,
        false,
        &mut DescriptorComputedState::default(),
        None,
        None,
    )?;

    // RDKit❗✔️:   auto burden = make_burden(*mol);
    // RDKit❗✔️:   auto atom_bcut = BCUT2D(burden, masses);
    // RDKit❗✔️:   auto gasteiger = BCUT2D(burden, charges);
    // RDKit❗✔️:   auto logp = BCUT2D(burden, slogp);
    // RDKit❗✔️:   auto mr = BCUT2D(burden, cmr);
    // RDKit❗✔️:   std::vector<double> res = {
    // RDKit❗✔️:       atom_bcut.first, atom_bcut.second, gasteiger.first, gasteiger.second,
    // RDKit❗✔️:       logp.first,      logp.second,      mr.first,        mr.second};
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    burden::bcut2d(
        topology,
        &masses,
        &own_charges,
        &crippen.logp,
        &crippen.molar_refractivity,
    )
    .map_err(|error| match error {
        DescriptorError::UnsupportedBcutBond { bond, order } => working
            .mapping
            .bonds()
            .new_to_old()
            .get(bond.index())
            .copied()
            .flatten()
            .map_or_else(
                || missing("source bond mapping"),
                |bond| DescriptorError::UnsupportedBcutBond { bond, order },
            ),
        other => other,
    })
}
