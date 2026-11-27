//! Gasteiger partial charges from RDKit `4ce8a06`, `Code/GraphMol/PartialCharges`.
//! GasteigerCharges.cpp: Copyright (C) 2003-2011 Greg Landrum and Rational Discovery LLC.
//! GasteigerParams.cpp: Copyright (C) 2003-2006 Rational Discovery LLC. See LICENSE for BSD terms.

use cosmolkit_model::{Hybridization, TopologyBlock};
use thiserror::Error;

use crate::{
    ValenceAssignment, ValenceError, rdkit_element_symbol, total_hydrogen_count_from_cache,
};

const MAX_ITERATIONS: i64 = i32::MAX as i64;
const IONX_HYDROGEN: f64 = 20.02;
const DAMP: f64 = 0.5;
const DAMP_SCALE: f64 = 0.5;

/// Parameters for the pinned RDKit Gasteiger charge iteration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GasteigerChargeParams {
    iterations: u32,
    throw_on_param_failure: bool,
}

impl GasteigerChargeParams {
    /// Construct checked parameters; the source loop accepts a signed 32-bit count.
    pub fn new(
        iterations: i64,
        throw_on_param_failure: bool,
    ) -> Result<Self, GasteigerChargeParamError> {
        Ok(Self {
            iterations: checked_iterations(iterations)?,
            throw_on_param_failure,
        })
    }

    /// Number of electronegativity-equalization iterations.
    #[must_use]
    pub const fn iterations(&self) -> u32 {
        self.iterations
    }

    /// Whether missing atom/mode rows produce an error instead of the source fallback.
    #[must_use]
    pub const fn throw_on_param_failure(&self) -> bool {
        self.throw_on_param_failure
    }

    /// Replace the iteration count only after validating the source's signed-int range.
    pub fn set_iterations(&mut self, iterations: i64) -> Result<(), GasteigerChargeParamError> {
        self.iterations = checked_iterations(iterations)?;
        Ok(())
    }

    /// Set strict parameter lookup behavior.
    pub const fn set_throw_on_param_failure(&mut self, value: bool) {
        self.throw_on_param_failure = value;
    }
}

impl Default for GasteigerChargeParams {
    fn default() -> Self {
        Self {
            iterations: 12,
            throw_on_param_failure: false,
        }
    }
}

fn checked_iterations(iterations: i64) -> Result<u32, GasteigerChargeParamError> {
    if !(0..=MAX_ITERATIONS).contains(&iterations) {
        return Err(GasteigerChargeParamError::InvalidIterations { iterations });
    }

    Ok(iterations as u32)
}

/// Invalid construction or replacement of a Gasteiger iteration count.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum GasteigerChargeParamError {
    #[error("Gasteiger iterations must be in 0..={MAX_ITERATIONS}, got {iterations}")]
    InvalidIterations { iterations: i64 },
}

/// The detached partial charge of one graph atom and its non-vertex hydrogens.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasteigerCharge {
    atom_charge: f64,
    attached_hydrogen_charge: f64,
}

impl GasteigerCharge {
    /// Charge assigned to this graph atom, in elementary-charge units.
    #[must_use]
    pub const fn atom_charge(&self) -> f64 {
        self.atom_charge
    }

    /// Sum of charges on this atom's attached hydrogens that are not graph vertices.
    #[must_use]
    pub const fn attached_hydrogen_charge(&self) -> f64 {
        self.attached_hydrogen_charge
    }
}

/// A Gasteiger calculation could not reproduce the source outcome for its input.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum GasteigerChargeError {
    #[error("Gasteiger parameters are missing for element {element:?}, mode {mode:?}")]
    MissingParameters { element: String, mode: String },
    #[error("Gasteiger attached-hydrogen lookup failed: {source}")]
    Valence {
        #[source]
        source: ValenceError,
    },
    #[error("Gasteiger element-symbol lookup failed for atomic number {atomic_number}: {source}")]
    ElementSymbol {
        atomic_number: u8,
        #[source]
        source: ValenceError,
    },
}

/// Calculate ordered atom charges without modifying the supplied topology.
///
/// The borrowed valence assignment is consulted only when at least one
/// iteration is requested and is never prepared or installed by this function.
pub fn gasteiger_charges(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
    params: &GasteigerChargeParams,
) -> Result<Vec<GasteigerCharge>, GasteigerChargeError> {
    // RDKit❗✔️: void computeGasteigerCharges(const ROMol &mol, std::vector<double> &charges,
    // RDKit❗✔️:                              int nIter, bool throwOnParamFailure) {
    // RDKit❗✔️:   PRECONDITION(charges.size() >= mol.getNumAtoms(), "bad array size");
    // RDKit❗✔️:   const PeriodicTable *table = PeriodicTable::getTable();
    // RDKit❗✔️:   const GasteigerParams *params = GasteigerParams::getParams();
    // RDKit❗✔️:   const auto natms = mol.getNumAtoms();
    // RDKit❗✔️:   // space for parameters for each atom in the molecule
    // RDKit❗✔️:   std::vector<DOUBLE_VECT> atmPs;
    // RDKit❗✔️:   atmPs.reserve(natms);
    // RDKit❗✔️:   std::fill(charges.begin(), charges.end(), 0.0);
    // RDKit❗✔️:   DOUBLE_VECT hChrg(
    // RDKit❗✔️:       natms,
    // RDKit❗✔️:       0.0);  // total charge for the implicit hydrogens on each heavy atom
    // RDKit❗✔️:   DOUBLE_VECT ionX(natms, 0.0);
    // RDKit❗✔️:   DOUBLE_VECT energ(natms, 0.0);
    // RDKit❗✔️:   // deal with the conjugated system - distribute the formal charges on atoms of
    // RDKit❗✔️:   // same type in each conjugated system
    // RDKit❗✔️:   Gasteiger::splitChargeConjugated(mol, charges);
    // RDKit❗✔️:   for (const auto atom : mol.atoms()) {
    // RDKit❗✔️:     const std::string elem = table->getElementSymbol(atom->getAtomicNum());
    // RDKit❗✔️:     std::string mode;
    let atom_count = topology.atoms.len();
    if atom_count == 0 {
        return Ok(Vec::new());
    }

    let mut charges = vec![0.0; atom_count];
    let mut hydrogen_charges = vec![0.0; atom_count];

    split_conjugated_formal_charges(topology, &mut charges);

    let mut atom_parameters = Vec::with_capacity(atom_count);
    let mut ionization_parameters = Vec::with_capacity(atom_count);
    for atom in &topology.atoms {
        let element = rdkit_element_symbol(atom.atomic_number()).map_err(|source| {
            GasteigerChargeError::ElementSymbol {
                atomic_number: atom.atomic_number(),
                source,
            }
        })?;
        let mode = atom_mode(topology, atom.id().index());
        // RDKit❗✔️:     // if we get a unknown mode or element type the following will will throw an
        // RDKit❗✔️:     // exception
        // RDKit❗✔️:     atmPs.push_back(params->getParams(elem, mode, throwOnParamFailure));
        // RDKit❗✔️:     // set ionX parameters
        // RDKit❗✔️:     // if Hydrogen treat differently
        // RDKit❗✔️:     const auto idx = atom->getIdx();
        // RDKit❗✔️:     if (atom->getAtomicNum() == 1) {
        // RDKit❗✔️:       ionX[idx] = IONXH;
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       ionX[idx] = atmPs[idx][0] + atmPs[idx][1] + atmPs[idx][2];
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        let parameter = resolve_parameters(element, mode, params.throw_on_param_failure)?;
        let ionization = if atom.atomic_number() == 1 {
            IONX_HYDROGEN
        } else {
            parameter[0] + parameter[1] + parameter[2]
        };
        atom_parameters.push(parameter);
        ionization_parameters.push(ionization);
    }

    // RDKit❗✔️:   // parameters for hydrogen atoms for case where the hydrogen are not in the
    // RDKit❗✔️:   // graph (implicit hydrogens)
    // RDKit❗✔️:   DOUBLE_VECT hParams = params->getParams("H", "*", throwOnParamFailure);
    let hydrogen_parameters = resolve_parameters("H", "*", params.throw_on_param_failure)?;
    let hydrogen_counts = if params.iterations > 0 {
        topology
            .atoms
            .iter()
            .map(|atom| {
                total_hydrogen_count_from_cache(topology, valence, atom.id(), false)
                    .map_err(|source| GasteigerChargeError::Valence { source })
            })
            .collect::<Result<Vec<_>, _>>()?
    } else {
        vec![0; atom_count]
    };

    let mut energies = vec![0.0; atom_count];
    // RDKit❗✔️:   double damp = DAMP;
    // RDKit❗✔️:   for (auto itx = 0u; itx < static_cast<unsigned int>(nIter); itx++) {
    // RDKit❗✔️:     for (auto aix = 0u; aix < natms; aix++) {
    // RDKit❗✔️:       auto enr = atmPs[aix][0] +
    // RDKit❗✔️:                  charges[aix] * (atmPs[aix][1] + atmPs[aix][2] * charges[aix]);
    // RDKit❗✔️:       energ[aix] = enr;
    // RDKit❗✔️:     }
    let mut damping = DAMP;
    for _ in 0..params.iterations {
        for atom_index in 0..atom_count {
            let parameter = atom_parameters[atom_index];
            energies[atom_index] = parameter[0]
                + charges[atom_index] * (parameter[1] + parameter[2] * charges[atom_index]);
        }

        for atom_index in 0..atom_count {
            // RDKit❗✔️:     for (const auto atom : mol.atoms()) {
            // RDKit❗✔️:       const auto aix = atom->getIdx();
            // RDKit❗✔️:       double dq = 0.0;
            // RDKit❗✔️:       for (const auto nbrAt : mol.atomNeighbors(atom)) {
            // RDKit❗✔️:         const auto nbr = nbrAt->getIdx();
            // RDKit❗✔️:         const auto dx = energ[nbr] - energ[aix];
            // RDKit❗✔️:         int sgn;
            // RDKit❗✔️:         if (dx < 0.0) {
            // RDKit❗✔️:           sgn = 0;
            // RDKit❗✔️:         } else {
            // RDKit❗✔️:           sgn = 1;
            // RDKit❗✔️:         }
            // RDKit❗✔️:         dq += dx / ((sgn * (ionX[aix] - ionX[nbr])) + ionX[nbr]);
            // RDKit❗✔️:       }
            let mut charge_delta = 0.0;
            for neighbor in topology.adjacency.neighbors_of(atom_index) {
                let neighbor_index = neighbor.atom_index;
                let energy_difference = energies[neighbor_index] - energies[atom_index];
                let sign = if energy_difference < 0.0 { 0.0 } else { 1.0 };
                let denominator = (sign
                    * (ionization_parameters[atom_index] - ionization_parameters[neighbor_index]))
                    + ionization_parameters[neighbor_index];
                charge_delta += energy_difference / denominator;
            }

            // RDKit❗✔️:       // now loop over the implicit hydrogens and get their contributions
            // RDKit❗✔️:       // since hydrogens don't connect to anything else, update their charges at
            // RDKit❗✔️:       // the same time
            // RDKit❗✔️:       const auto niHs = atom->getTotalNumHs();
            // RDKit❗✔️:       if (niHs > 0) {
            // RDKit❗✔️:         const auto qHs = hChrg[aix] / niHs;
            // RDKit❗✔️:         const auto enr = hParams[0] + qHs * (hParams[1] + hParams[2] * qHs);
            // RDKit❗✔️:         const auto dx = enr - energ[aix];
            // RDKit❗✔️:         int sgn;
            // RDKit❗✔️:         if (dx < 0.0) {
            // RDKit❗✔️:           sgn = 0;
            // RDKit❗✔️:         } else {
            // RDKit❗✔️:           sgn = 1;
            // RDKit❗✔️:         }
            // RDKit❗✔️:         const auto dqH = dx / ((sgn * (ionX[aix] - IONXH)) + IONXH);
            // RDKit❗✔️:         dq += (niHs * dqH);
            // RDKit❗✔️:         // adjust the charges on the hydrogens simultaneously (possible because
            // RDKit❗✔️:         // each of the hydrogens have no other neighbors)
            // RDKit❗✔️:         hChrg[aix] -= (niHs * dqH * damp);
            // RDKit❗✔️:       }
            let hydrogen_count = hydrogen_counts[atom_index];
            if hydrogen_count > 0 {
                let hydrogen_charge_per_atom =
                    hydrogen_charges[atom_index] / f64::from(hydrogen_count);
                let hydrogen_energy = hydrogen_parameters[0]
                    + hydrogen_charge_per_atom
                        * (hydrogen_parameters[1]
                            + hydrogen_parameters[2] * hydrogen_charge_per_atom);
                let energy_difference = hydrogen_energy - energies[atom_index];
                let sign = if energy_difference < 0.0 { 0.0 } else { 1.0 };
                let hydrogen_charge_delta = energy_difference
                    / ((sign * (ionization_parameters[atom_index] - IONX_HYDROGEN))
                        + IONX_HYDROGEN);
                charge_delta += f64::from(hydrogen_count) * hydrogen_charge_delta;
                hydrogen_charges[atom_index] -=
                    f64::from(hydrogen_count) * hydrogen_charge_delta * damping;
            }

            // RDKit❗✔️:       charges[aix] += (damp * dq);
            // RDKit❗✔️:     }
            charges[atom_index] += damping * charge_delta;
        }

        // RDKit❗✔️:     damp *= DAMP_SCALE;
        // RDKit❗✔️:   }
        damping *= DAMP_SCALE;
    }

    // RDKit❗✔️:   for (auto atom : mol.atoms()) {
    // RDKit❗✔️:     const auto aix = atom->getIdx();
    // RDKit❗✔️:     // set the charge on the heavy atom
    // RDKit❗✔️:     atom->setProp(common_properties::_GasteigerCharge, charges[aix], true);
    // RDKit❗✔️:     // set the implicit hydrogen charges
    // RDKit❗✔️:     atom->setProp(common_properties::_GasteigerHCharge, hChrg[aix], true);
    // RDKit❗✔️:   }
    Ok(charges
        .into_iter()
        .zip(hydrogen_charges)
        .map(|(atom_charge, attached_hydrogen_charge)| GasteigerCharge {
            atom_charge,
            attached_hydrogen_charge,
        })
        .collect())
}

fn atom_mode(topology: &TopologyBlock, atom_index: usize) -> &'static str {
    let atom = &topology.atoms[atom_index];

    // RDKit source (GasteigerCharges.cpp, atom mode dispatch):
    // RDKit❗✔️: switch (atom->getHybridization()) {
    // RDKit❗✔️:   case Atom::SP3:
    // RDKit❗✔️:     mode = "sp3";
    // RDKit❗✔️:     break;
    // RDKit❗✔️:   case Atom::SP2:
    // RDKit❗✔️:     mode = "sp2";
    // RDKit❗✔️:     break;
    // RDKit❗✔️:   case Atom::SP:
    // RDKit❗✔️:     mode = "sp";
    // RDKit❗✔️:     break;
    // RDKit❗✔️:   default:
    // RDKit❗✔️:     if (atom->getAtomicNum() == 1) {
    // RDKit❗✔️:       mode = "*";
    // RDKit❗✔️:     } else if (atom->getAtomicNum() == 16) {
    // RDKit❗✔️:       int no = 0;
    // RDKit❗✔️:       for (const auto nbrAt : mol.atomNeighbors(atom)) {
    // RDKit❗✔️:         if (nbrAt->getAtomicNum() == 8) {
    // RDKit❗✔️:           no++;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (no == 2) {
    // RDKit❗✔️:         mode = "so2";
    // RDKit❗✔️:       } else if (no == 1) {
    // RDKit❗✔️:         mode = "so";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         mode = "sp3";
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️: }
    match atom.hybridization() {
        Hybridization::Sp3 => "sp3",
        Hybridization::Sp2 => "sp2",
        Hybridization::Sp => "sp",
        _ if atom.atomic_number() == 1 => "*",
        _ if atom.atomic_number() == 16 => {
            let oxygen_neighbors = topology
                .adjacency
                .neighbors_of(atom_index)
                .iter()
                .filter(|neighbor| topology.atoms[neighbor.atom_index].atomic_number() == 8)
                .count();
            match oxygen_neighbors {
                2 => "so2",
                1 => "so",
                _ => "sp3",
            }
        }
        _ => "",
    }
}

fn split_conjugated_formal_charges(topology: &TopologyBlock, charges: &mut [f64]) {
    // RDKit source (GasteigerCharges.cpp, splitChargeConjugated):
    // RDKit❗✔️: for (const auto at : mol.atoms()) {
    // RDKit❗✔️:   auto aix = at->getIdx();
    // RDKit❗✔️:   double formal = at->getFormalCharge();
    // RDKit❗✔️:   marker.resize(0);
    // RDKit❗✔️:   if ((fabs(formal) > EPS_DOUBLE) && (fabs(charges[aix]) < EPS_DOUBLE)) {
    // RDKit❗✔️:     marker.push_back(aix);
    // RDKit❗✔️:     for (const auto bnd1 : mol.atomBonds(at)) {
    // RDKit❗✔️:       if (bnd1->getIsConjugated()) {
    // RDKit❗✔️:         auto aax = bnd1->getOtherAtomIdx(aix);
    // RDKit❗✔️:         auto aat = mol.getAtomWithIdx(aax);
    // RDKit❗✔️:         for (const auto bnd2 : mol.atomBonds(aat)) {
    // RDKit❗✔️:           if (bnd1 != bnd2) {
    // RDKit❗✔️:             if (bnd2->getIsConjugated()) {
    // RDKit❗✔️:               auto yax = bnd2->getOtherAtomIdx(aax);
    // RDKit❗✔️:               auto yat = mol.getAtomWithIdx(yax);
    // RDKit❗✔️:               if (at->getAtomicNum() == yat->getAtomicNum()) {
    // RDKit❗✔️:                 formal += yat->getFormalCharge();
    // RDKit❗✔️:                 marker.push_back(yax);
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (const auto mci : marker) {
    // RDKit❗✔️:       charges[mci] = (formal / marker.size());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    for (atom_index, atom) in topology.atoms.iter().enumerate() {
        let mut formal_charge = f64::from(atom.formal_charge());
        if !(formal_charge.abs() > f64::EPSILON && charges[atom_index].abs() < f64::EPSILON) {
            continue;
        }

        let mut marker = vec![atom_index];
        for first_neighbor in topology.adjacency.neighbors_of(atom_index) {
            let first_bond = &topology.bonds[first_neighbor.bond.index()];
            if !first_bond.is_conjugated() {
                continue;
            }

            let middle_index = first_neighbor.atom_index;
            for second_neighbor in topology.adjacency.neighbors_of(middle_index) {
                if first_bond.id() == second_neighbor.bond {
                    continue;
                }

                let second_bond = &topology.bonds[second_neighbor.bond.index()];
                if second_bond.is_conjugated()
                    && atom.atomic_number()
                        == topology.atoms[second_neighbor.atom_index].atomic_number()
                {
                    formal_charge +=
                        f64::from(topology.atoms[second_neighbor.atom_index].formal_charge());
                    marker.push(second_neighbor.atom_index);
                }
            }
        }

        let distributed_charge = formal_charge / marker.len() as f64;
        for marker_index in marker {
            charges[marker_index] = distributed_charge;
        }
    }
}

fn resolve_parameters(
    element: &str,
    mode: &str,
    throw_on_param_failure: bool,
) -> Result<[f64; 3], GasteigerChargeError> {
    // RDKit❗✔️:   DOUBLE_VECT getParams(std::string elem, std::string mode,
    // RDKit❗✔️:                         bool throwOnFailure = false) const {
    // RDKit❗✔️:     std::pair<std::string, std::string> query(elem, mode);
    // RDKit❗✔️:     std::map<std::pair<std::string, std::string>, DOUBLE_VECT>::const_iterator
    // RDKit❗✔️:         iter;
    // RDKit❗✔️:     iter = d_paramMap.find(query);
    // RDKit❗✔️:     if (iter != d_paramMap.end()) {
    // RDKit❗✔️:       return iter->second;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       if (throwOnFailure) {
    // RDKit❗✔️:         std::string message =
    // RDKit❗✔️:             "ERROR: No Gasteiger Partial Charge parameters for Element: ";
    // RDKit❗✔️:         message += elem;
    // RDKit❗✔️:         message += " Mode: ";
    // RDKit❗✔️:         message += mode;
    // RDKit❗✔️:         throw ValueErrorException(message);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         iter =
    // RDKit❗✔️:             d_paramMap.find(std::make_pair(std::string("X"), std::string("*")));
    // RDKit❗✔️:         if (iter != d_paramMap.end()) {
    // RDKit❗✔️:           return iter->second;
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           std::string message =
    // RDKit❗✔️:               "ERROR: Default Gasteiger Partial Charge parameters are missing";
    // RDKit❗✔️:           throw ValueErrorException(message);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    if let Some(parameters) = parameter_row(element, mode) {
        return Ok(parameters);
    }

    if throw_on_param_failure {
        return Err(GasteigerChargeError::MissingParameters {
            element: element.to_owned(),
            mode: mode.to_owned(),
        });
    }

    // RDKit's permissive lookup asks for the explicit X/* table row.
    parameter_row("X", "*").ok_or_else(|| GasteigerChargeError::MissingParameters {
        element: "X".to_owned(),
        mode: "*".to_owned(),
    })
}

fn parameter_row(element: &str, mode: &str) -> Option<[f64; 3]> {
    // RDKit source (GasteigerParams.cpp, defaultParamData):
    // RDKit❗✔️: std::string defaultParamData =
    // RDKit❗✔️:     "H       *      7.17    6.24    -0.56 \n \
    // RDKit❗✔️: C       sp3     7.98    9.18    1.88 \n \
    // RDKit❗✔️: C       sp2     8.79    9.32    1.51 \n \
    // RDKit❗✔️: C       sp      10.39   9.45    0.73 \n \
    // RDKit❗✔️: N       sp3     11.54   10.82   1.36 \n \
    // RDKit❗✔️: N       sp2     12.87   11.15   0.85 \n \
    // RDKit❗✔️: N       sp      15.68   11.7    -0.27 \n \
    // RDKit❗✔️: O       sp3     14.18   12.92   1.39 \n \
    // RDKit❗✔️: O       sp2     17.07   13.79   0.47 \n \
    // RDKit❗✔️: F       sp3     14.66   13.85   2.31 \n \
    // RDKit❗✔️: Cl      sp3     11.00   9.69    1.35 \n \
    // RDKit❗✔️: Br      sp3     10.08   8.47    1.16 \n \
    // RDKit❗✔️: I       sp3     9.9     7.96    0.96 \n \
    // RDKit❗✔️: S       sp3     10.14   9.13    1.38 \n \
    // RDKit❗✔️: S       so      10.14   9.13    1.38 \n \
    // RDKit❗✔️: S       so2     12.00   10.81   1.20 \n \
    // RDKit❗✔️: S       sp2     10.88   9.49    1.33 \n \
    // RDKit❗✔️: P       sp3     8.90    8.24    0.96 \n \
    // RDKit❗✔️: X       *       0.00    0.00    0.00 \n \
    // RDKit source (GasteigerParams.cpp, additionalParamData):
    // RDKit❗✔️: std::string additionalParamData =
    // RDKit❗✔️:     "P       sp2     9.665   8.530   0.735 \n \
    // RDKit❗✔️: Si      sp3     7.300   6.567   0.657 \n \
    // RDKit❗✔️: Si      sp2     7.905   6.748   0.443 \n \
    // RDKit❗✔️: Si      sp      9.065   7.027  -0.002 \n \
    // RDKit❗✔️: B       sp3     5.980   6.820   1.605 \n \
    // RDKit❗✔️: B       sp2     6.420   6.807   1.322 \n \
    // RDKit❗✔️: Be      sp3     3.845   6.755   3.165 \n \
    // RDKit❗✔️: Be      sp2     4.005   6.725   3.035 \n \
    // RDKit❗✔️: Mg      sp2     3.565   5.572   2.197 \n \
    // RDKit❗✔️: Mg      sp3     3.300   5.587   2.447 \n \
    // RDKit❗✔️: Mg      sp      4.040   5.472   1.823 \n \
    // RDKit❗✔️: Al      sp3     5.375   4.953   0.867 \n \
    // RDKit❗✔️: Al      sp2     5.795   5.020   0.695 \n \
    match (element, mode) {
        ("H", "*") => Some([7.17, 6.24, -0.56]),
        ("C", "sp3") => Some([7.98, 9.18, 1.88]),
        ("C", "sp2") => Some([8.79, 9.32, 1.51]),
        ("C", "sp") => Some([10.39, 9.45, 0.73]),
        ("N", "sp3") => Some([11.54, 10.82, 1.36]),
        ("N", "sp2") => Some([12.87, 11.15, 0.85]),
        ("N", "sp") => Some([15.68, 11.7, -0.27]),
        ("O", "sp3") => Some([14.18, 12.92, 1.39]),
        ("O", "sp2") => Some([17.07, 13.79, 0.47]),
        ("F", "sp3") => Some([14.66, 13.85, 2.31]),
        ("Cl", "sp3") => Some([11.00, 9.69, 1.35]),
        ("Br", "sp3") => Some([10.08, 8.47, 1.16]),
        ("I", "sp3") => Some([9.9, 7.96, 0.96]),
        ("S", "sp3") => Some([10.14, 9.13, 1.38]),
        ("S", "so") => Some([10.14, 9.13, 1.38]),
        ("S", "so2") => Some([12.00, 10.81, 1.20]),
        ("S", "sp2") => Some([10.88, 9.49, 1.33]),
        ("P", "sp3") => Some([8.90, 8.24, 0.96]),
        ("X", "*") => Some([0.0, 0.0, 0.0]),
        ("P", "sp2") => Some([9.665, 8.530, 0.735]),
        ("Si", "sp3") => Some([7.300, 6.567, 0.657]),
        ("Si", "sp2") => Some([7.905, 6.748, 0.443]),
        ("Si", "sp") => Some([9.065, 7.027, -0.002]),
        ("B", "sp3") => Some([5.980, 6.820, 1.605]),
        ("B", "sp2") => Some([6.420, 6.807, 1.322]),
        ("Be", "sp3") => Some([3.845, 6.755, 3.165]),
        ("Be", "sp2") => Some([4.005, 6.725, 3.035]),
        ("Mg", "sp2") => Some([3.565, 5.572, 2.197]),
        ("Mg", "sp3") => Some([3.300, 5.587, 2.447]),
        ("Mg", "sp") => Some([4.040, 5.472, 1.823]),
        ("Al", "sp3") => Some([5.375, 4.953, 0.867]),
        ("Al", "sp2") => Some([5.795, 5.020, 0.695]),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec};
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

    fn bond(begin: usize, end: usize, order: BondOrder, conjugated: bool) -> BondSpec {
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order).with_conjugated(conjugated)
    }

    fn valence(explicit: &[i32], implicit: &[i32]) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: explicit.to_vec(),
            implicit_hydrogens: implicit.to_vec(),
        }
    }

    fn charges(
        topology: &TopologyBlock,
        valence: Option<&ValenceAssignment>,
        iterations: i64,
    ) -> Vec<GasteigerCharge> {
        gasteiger_charges(
            topology,
            valence,
            &GasteigerChargeParams::new(iterations, false).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn checked_iteration_setter_preserves_prior_value_on_rejection() {
        let mut params = GasteigerChargeParams::default();
        params.set_iterations(0).unwrap();
        assert_eq!(params.iterations(), 0);
        for invalid in [-1, i64::from(i32::MAX) + 1] {
            assert!(params.set_iterations(invalid).is_err());
            assert_eq!(params.iterations(), 0);
        }
        assert_eq!(GasteigerChargeParams::default().iterations(), 12);
    }

    #[test]
    fn zero_iterations_initializes_formal_charge_without_reading_hydrogens() {
        let graph = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::O),
                AtomSpec::new(Element::O).with_formal_charge(-1),
            ],
            vec![
                bond(0, 1, BondOrder::Single, true),
                bond(1, 2, BondOrder::Double, true),
                bond(1, 3, BondOrder::Single, true),
            ],
        );
        let result = charges(&graph, None, 0);
        assert_eq!(
            result
                .iter()
                .map(GasteigerCharge::atom_charge)
                .collect::<Vec<_>>(),
            [0.0, 0.0, -0.5, -0.5]
        );
        assert!(
            result
                .iter()
                .all(|row| row.attached_hydrogen_charge() == 0.0)
        );
    }

    #[test]
    fn positive_iterations_report_missing_implicit_hydrogen_cache() {
        let graph = topology(vec![AtomSpec::new(Element::C)], Vec::new());
        let error = gasteiger_charges(&graph, None, &GasteigerChargeParams::default()).unwrap_err();
        assert!(
            matches!(error, GasteigerChargeError::Valence { source: ValenceError::ImplicitValenceCacheNotInitialized { atom } } if atom == AtomId::new(0))
        );
    }

    #[test]
    fn strict_unknown_parameter_is_not_affected_by_permissive_fallback() {
        let graph = topology(
            vec![AtomSpec::new(Element::NA).with_formal_charge(1)],
            Vec::new(),
        );
        let default = charges(&graph, None, 0);
        assert_eq!(default[0].atom_charge(), 1.0);
        let strict = GasteigerChargeParams::new(0, true).unwrap();
        assert!(matches!(
            gasteiger_charges(&graph, None, &strict),
            Err(GasteigerChargeError::MissingParameters { element, .. }) if element == "Na"
        ));
    }
}
