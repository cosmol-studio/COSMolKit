use super::angle::calc_angle_force_constant;
use super::bond::{BondMathError, calc_bond_force_constant, calc_bond_rest_length};
use super::nonbonded::{calc_nonbonded_depth, calc_nonbonded_minimum};
use super::params::{
    AtomicParams, ParamCollection, RAD2DEG, UffAngle, UffBond, UffInv, UffTor, UffVdw,
};
use super::torsion::{equation17, is_in_group6};
use super::utils::calc_inversion_coefficients;
use cosmolkit_core::{
    PropertyStringError, ValenceAssignment, ValenceError, periodic_table_outer_electrons,
    property_value_to_string, rdkit_default_valence, rdkit_element_symbol,
};
use cosmolkit_model::{Atom, AtomId, Bond, BondId, Hybridization, TopologyBlock};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum UffTypingInput {
    TotalValence,
    ConjugatedBondPresence,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) enum UffTypingError {
    CoreValence(ValenceError),
    CorePropertyString(PropertyStringError),
    BondMath(BondMathError),
    PreparedStateLength {
        input: UffTypingInput,
        expected: usize,
        actual: usize,
    },
    AtomIdOutOfBounds {
        atom_id: AtomId,
        atom_count: usize,
    },
    AtomIndexOutOfBounds {
        index: usize,
        atom_count: usize,
    },
    BondIdOutOfBounds {
        bond_id: BondId,
        bond_count: usize,
    },
    CenterBondMissingAfterSourceSuccess {
        idx2: usize,
        idx3: usize,
    },
    CentralParameterSlotMissingAfterSourceSuccess {
        atom_index: usize,
    },
}

impl std::fmt::Display for UffTypingError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for UffTypingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::CoreValence(source) => Some(source),
            Self::CorePropertyString(source) => Some(source),
            Self::BondMath(source) => Some(source),
            Self::PreparedStateLength { .. } => None,
            Self::AtomIdOutOfBounds { .. } => None,
            Self::AtomIndexOutOfBounds { .. } => None,
            Self::BondIdOutOfBounds { .. } => None,
            Self::CenterBondMissingAfterSourceSuccess { .. } => None,
            Self::CentralParameterSlotMissingAfterSourceSuccess { .. } => None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum UffAtomStateRef<'a> {
    Cached {
        topology: &'a TopologyBlock,
        assignment: &'a ValenceAssignment,
    },
    SuppliedRows {
        total_valences: &'a [i32],
        conjugated_presence: &'a [bool],
    },
}

#[cfg(test)]
std::thread_local! {
    static UFF_PREP_ATOM_STATE_CONJUGATION_READS: std::cell::Cell<usize> = const {
        std::cell::Cell::new(0)
    };
}

impl<'a> UffAtomStateRef<'a> {
    pub(super) fn cached(
        topology: &'a TopologyBlock,
        assignment: &'a ValenceAssignment,
    ) -> Result<Self, super::builder::UffBuilderError> {
        // Match prepare_parameter_query's old order: validate both cached rows
        // before the complete topology/adjacency validation.
        super::builder::validate_typing_valence_cache(topology, assignment)?;
        topology
            .validate()
            .map_err(super::builder::UffBuilderError::TopologyValidation)?;
        Ok(Self::Cached {
            topology,
            assignment,
        })
    }

    pub(super) fn supplied_rows(
        topology: &TopologyBlock,
        total_valences: &'a [i32],
        conjugated_presence: &'a [bool],
    ) -> Result<Self, UffTypingError> {
        let atom_count = topology.atoms.len();
        if total_valences.len() != atom_count {
            return Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: atom_count,
                actual: total_valences.len(),
            });
        }
        if conjugated_presence.len() != atom_count {
            return Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::ConjugatedBondPresence,
                expected: atom_count,
                actual: conjugated_presence.len(),
            });
        }
        Ok(Self::SuppliedRows {
            total_valences,
            conjugated_presence,
        })
    }

    pub(super) fn total_valence_at(self, atom_index: usize) -> i32 {
        match self {
            Self::Cached {
                topology,
                assignment,
            } => {
                // RDKit✔️✔️: unsigned int Atom::getValence(ValenceType which) const {
                // RDKit✔️✔️:   if (!dp_mol) {
                // RDKit✔️✔️:     return 0;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   PRECONDITION(
                // RDKit✔️✔️:       (which == ValenceType::IMPLICIT || d_explicitValence > -1),
                // RDKit✔️✔️:       "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()");
                // RDKit✔️✔️:   PRECONDITION(
                // RDKit✔️✔️:       (which == ValenceType::EXPLICIT || df_noImplicit ||
                // RDKit✔️✔️:        d_implicitValence > -1),
                // RDKit✔️✔️:       "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()");
                // RDKit✔️✔️:   if (which == ValenceType::EXPLICIT) {
                // RDKit✔️✔️:     return d_explicitValence;
                // RDKit✔️✔️:   } else {
                // RDKit✔️✔️:     return df_noImplicit ? 0 : d_implicitValence;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️: }
                // RDKit✔️✔️: unsigned int Atom::getTotalValence() const {
                // RDKit✔️✔️:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
                // RDKit✔️✔️: }
                // Behavior review — RDKit✔️✔️: cached construction performs the
                // source explicit/implicit precondition and range checks before
                // this read. The source noImplicit flag selects implicit zero;
                // the result remains the same explicit-plus-effective-implicit
                // sum for every admitted cached row.
                // Complexity review — RDKit✔️✔️: two indexed scalar reads and
                // one addition are O(1), matching the source getter pair. This
                // path creates no per-atom projection; caller cache and topology
                // validation costs remain accounted for at their owners.
                let explicit = assignment.explicit_valence[atom_index];
                let implicit = if topology.atoms[atom_index].no_implicit() {
                    0
                } else {
                    assignment.implicit_hydrogens[atom_index]
                };
                explicit + implicit
            }
            Self::SuppliedRows { total_valences, .. } => total_valences[atom_index],
        }
    }

    pub(super) fn conjugated_presence_at(self, atom_index: usize) -> bool {
        #[cfg(test)]
        UFF_PREP_ATOM_STATE_CONJUGATION_READS.with(|reads| reads.set(reads.get() + 1));

        match self {
            Self::Cached { topology, .. } => {
                super::builder::atom_has_conjugated_bond_from_validated_topology(
                    topology, atom_index,
                )
            }
            Self::SuppliedRows {
                conjugated_presence,
                ..
            } => conjugated_presence[atom_index],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum UffTypingDiagnosticKind {
    Warning,
    Error,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct UffTypingDiagnostic {
    pub(super) atom_id: Option<AtomId>,
    pub(super) kind: UffTypingDiagnosticKind,
    pub(super) message_prefix: &'static str,
}

const UNRECOGNIZED_CHARGE_STATE_MESSAGE: &str = "UFFTYPER: Unrecognized charge state for atom: ";
const UNRECOGNIZED_ATOM_TYPE_MESSAGE: &str = "UFFTYPER: Unrecognized atom type: ";
const FORCED_SP3_HYBRIDIZATION_WARNING_MESSAGE: &str =
    "UFFTYPER: Warning: hybridization set to SP3 for atom ";
const FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE: &str =
    "UFFTYPER: Warning: hybridization set to SP for atom ";
const UNRECOGNIZED_HYBRIDIZATION_MESSAGE: &str = "UFFTYPER: Unrecognized hybridization for atom: ";
pub(super) const NEEDS_EXPLICIT_HYDROGENS_WARNING_MESSAGE: &str =
    "Molecule does not have explicit Hs. Consider calling AddHs()";

fn append_fixed_charge_flag(
    atom: &Atom,
    total_valence: i32,
    atom_key: &mut cosmolkit_model::PropertyText,
    tolerate_charge_mismatch: bool,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> bool {
    // RDKit❗✔️: void addAtomChargeFlags(const Atom *atom, std::string &atomKey,
    // RDKit❗✔️:                         bool tolerateChargeMismatch) {
    // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
    // A Rust reference makes the non-null atom precondition unrepresentable.
    // RDKit❗✔️:   int totalValence = atom->getTotalValence();
    // RDKit❗✔️:   int fc = atom->getFormalCharge();
    // RDKit❗✔️:   // FIX: come up with some way of handling metals here
    // RDKit❗✔️:   switch (atom->getAtomicNum()) {
    // RDKit❗✔️:     // Atoms only +1 in default UFF params
    // RDKit❗✔️:     case 29:  // Cu
    // RDKit❗✔️:     case 47:  // Ag
    // RDKit❗✔️:       if (totalValence == 1 || fc == 1 || tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+1";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:             << atom->getIdx() << std::endl;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     // Atoms only +2 in default UFF params
    // RDKit❗✔️:     case 4:   // Be
    // RDKit❗✔️:     case 20:  // Ca
    // RDKit❗✔️:     case 25:  // Mn
    // RDKit❗✔️:     case 26:  // Fe
    // RDKit❗✔️:     case 28:  // Ni
    // RDKit❗✔️:               //  case 30:  // Zn
    // RDKit❗✔️:     case 46:  // Pd
    // RDKit❗✔️:               //  case 48:  // Cd
    // RDKit❗✔️:     case 78:  // Pt
    // RDKit❗✔️:       if (totalValence == 2 || fc == 2 || tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+2";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:             << atom->getIdx() << std::endl;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     // Atoms only +3 in default UFF params
    // RDKit❗✔️:     case 21:   // Sc
    // RDKit❗✔️:     case 24:   // Cr
    // RDKit❗✔️:     case 27:   // Co
    // RDKit❗✔️:                //  case 49:  // In
    // RDKit❗✔️:     case 79:   // Au
    // RDKit❗✔️:     case 89:   // Ac
    // RDKit❗✔️:     case 96:   // Cm
    // RDKit❗✔️:     case 97:   // Bk
    // RDKit❗✔️:     case 98:   // Cf
    // RDKit❗✔️:     case 99:   // Es
    // RDKit❗✔️:     case 100:  // Fm
    // RDKit❗✔️:     case 101:  // Md
    // RDKit❗✔️:     case 102:  // No
    // RDKit❗✔️:     case 103:  // Lr/Lw
    // RDKit❗✔️:       if (totalValence == 3 || fc == 3 || tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+3";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:             << atom->getIdx() << std::endl;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     // Atoms only +4 in default UFF params
    // RDKit❗✔️:     case 2:   // He
    // RDKit❗✔️:     case 18:  // Ar
    // RDKit❗✔️:     case 22:  // Ti
    // RDKit❗✔️:     case 36:  // Kr
    // RDKit❗✔️:     case 54:  // Xe
    // RDKit❗✔️:     case 90:  // Th
    // RDKit❗✔️:     case 91:  // Pa
    // RDKit❗✔️:     case 92:  // U
    // RDKit❗✔️:     case 93:  // Np
    // RDKit❗✔️:     case 94:  // Pu
    // RDKit❗✔️:     case 95:  // Am
    // RDKit❗✔️:       if (totalValence == 4 || fc == 4 || tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+4";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:             << atom->getIdx() << std::endl;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     // Atoms only +5 in default UFF params
    // RDKit❗✔️:     case 23:  // V
    // RDKit❗✔️:     case 41:  // Nb
    // RDKit❗✔️:     case 43:  // Tc
    // RDKit❗✔️:     case 73:  // Ta
    // RDKit❗✔️:       if (totalValence == 5 || fc == 5 || tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+5";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:             << atom->getIdx() << std::endl;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     // Atoms only +6 in default UFF params
    // RDKit❗✔️:     case 42:  // Mo
    // RDKit❗✔️:       if (totalValence == 6 || fc == 6 || tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+6";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:             << atom->getIdx() << std::endl;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let formal_charge = i32::from(atom.formal_charge());
    let (suffix, required_valence) = match atom.atomic_number() {
        29 | 47 => ("+1", 1),
        4 | 20 | 25 | 26 | 28 | 46 | 78 => ("+2", 2),
        21 | 24 | 27 | 79 | 89 | 96..=103 => ("+3", 3),
        2 | 18 | 22 | 36 | 54 | 90..=95 => ("+4", 4),
        23 | 41 | 43 | 73 => ("+5", 5),
        42 => ("+6", 6),
        _ => return false,
    };

    if total_valence == required_valence
        || formal_charge == required_valence
        || tolerate_charge_mismatch
    {
        atom_key.extend_bytes(suffix.as_bytes());
    } else {
        diagnostics.push(UffTypingDiagnostic {
            atom_id: Some(atom.id()),
            kind: UffTypingDiagnosticKind::Error,
            message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
        });
    }
    true
}

fn append_valence_only_charge_flag(
    atom: &Atom,
    total_valence: i32,
    atom_key: &mut cosmolkit_model::PropertyText,
    tolerate_charge_mismatch: bool,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> bool {
    // RDKit❗✔️: case 12:  // Mg
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 2:
    // RDKit❗✔️:       atomKey += "+2";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+2";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 30:  // Zn
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 2:
    // RDKit❗✔️:       atomKey += "+2";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+2";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 31:  // Ga
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 3:
    // RDKit❗✔️:       atomKey += "+3";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+3";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 33:  // As
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 3:
    // RDKit❗✔️:       atomKey += "+3";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+3";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 34:  // Se
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 2:
    // RDKit❗✔️:       atomKey += "+2";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+2";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 48:  // Cd
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 2:
    // RDKit❗✔️:       atomKey += "+2";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+2";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 49:  // In
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 3:
    // RDKit❗✔️:       atomKey += "+3";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+3";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 51:  // Sb
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 3:
    // RDKit❗✔️:       atomKey += "+3";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+3";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 52:  // Te
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 2:
    // RDKit❗✔️:       atomKey += "+2";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+2";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 80:  // Hg
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 2:
    // RDKit❗✔️:       atomKey += "+2";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+2";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 81:  // Tl
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 3:
    // RDKit❗✔️:       atomKey += "+3";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+3";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 82:  // Pb
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 3:
    // RDKit❗✔️:       atomKey += "+3";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+3";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 83:  // Bi
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 3:
    // RDKit❗✔️:       atomKey += "+3";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+3";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 84:  // Po
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 2:
    // RDKit❗✔️:       atomKey += "+2";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+2";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;

    let (required_valence, suffix) = match atom.atomic_number() {
        12 | 30 | 34 | 48 | 52 | 80 | 84 => (2, "+2"),
        31 | 33 | 49 | 51 | 81 | 82 | 83 => (3, "+3"),
        _ => return false,
    };

    if total_valence == required_valence {
        atom_key.extend_bytes(suffix.as_bytes());
    } else {
        if tolerate_charge_mismatch {
            atom_key.extend_bytes(suffix.as_bytes());
        }
        diagnostics.push(UffTypingDiagnostic {
            atom_id: Some(atom.id()),
            kind: UffTypingDiagnosticKind::Error,
            message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
        });
    }
    true
}

fn check_unsuffixed_main_group_charge(
    atom: &Atom,
    total_valence: i32,
    _atom_key: &mut cosmolkit_model::PropertyText,
    _tolerate_charge_mismatch: bool,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> bool {
    // RDKit❗✔️: case 13:  // Al
    // RDKit❗✔️:   if (totalValence != 3) {
    // RDKit❗✔️:     BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:         << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:         << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    // RDKit❗✔️: case 14:  // Si
    // RDKit❗✔️:   if (totalValence != 4) {
    // RDKit❗✔️:     BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:         << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:         << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    let required_valence = match atom.atomic_number() {
        13 => 3,
        14 => 4,
        _ => return false,
    };

    if total_valence != required_valence {
        diagnostics.push(UffTypingDiagnostic {
            atom_id: Some(atom.id()),
            kind: UffTypingDiagnosticKind::Error,
            message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
        });
    }
    true
}

fn append_phosphorus_charge_flag(
    atom: &Atom,
    total_valence: i32,
    atom_key: &mut cosmolkit_model::PropertyText,
    tolerate_charge_mismatch: bool,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> bool {
    // RDKit❗✔️: case 15:  // P
    // RDKit❗✔️:   switch (totalValence) {
    // RDKit❗✔️:     case 3:
    // RDKit❗✔️:       atomKey += "+3";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case 5:
    // RDKit❗✔️:       atomKey += "+5";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       if (tolerateChargeMismatch) {
    // RDKit❗✔️:         atomKey += "+5";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:           << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:           << atom->getIdx() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    let 15 = atom.atomic_number() else {
        return false;
    };

    match total_valence {
        3 => atom_key.extend_bytes(b"+3"),
        5 => atom_key.extend_bytes(b"+5"),
        _ => {
            if tolerate_charge_mismatch {
                atom_key.extend_bytes(b"+5");
            }
            diagnostics.push(UffTypingDiagnostic {
                atom_id: Some(atom.id()),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
            });
        }
    }
    true
}

fn append_sulfur_charge_flag(
    atom: &Atom,
    total_valence: i32,
    atom_key: &mut cosmolkit_model::PropertyText,
    tolerate_charge_mismatch: bool,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> bool {
    // RDKit❗✔️: case 16:  // S
    // RDKit❗✔️:   if (atom->getHybridization() != Atom::SP2) {
    // RDKit❗✔️:     switch (totalValence) {
    // RDKit❗✔️:       case 2:
    // RDKit❗✔️:         atomKey += "+2";
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       case 4:
    // RDKit❗✔️:         atomKey += "+4";
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       case 6:
    // RDKit❗✔️:         atomKey += "+6";
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       default:
    // RDKit❗✔️:         if (tolerateChargeMismatch) {
    // RDKit❗✔️:           atomKey += "+6";
    // RDKit❗✔️:         }
    // RDKit❗✔️:         BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit❗✔️:             << atom->getIdx() << std::endl;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   break;
    if atom.atomic_number() != 16 {
        return false;
    }

    if atom.hybridization() == Hybridization::Sp2 {
        return true;
    }

    match total_valence {
        2 => atom_key.extend_bytes(b"+2"),
        4 => atom_key.extend_bytes(b"+4"),
        6 => atom_key.extend_bytes(b"+6"),
        _ => {
            if tolerate_charge_mismatch {
                atom_key.extend_bytes(b"+6");
            }
            diagnostics.push(UffTypingDiagnostic {
                atom_id: Some(atom.id()),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
            });
        }
    }
    true
}

fn rewrite_rhenium_charge_flag(
    atom: &Atom,
    _total_valence: i32,
    atom_key: &mut cosmolkit_model::PropertyText,
    tolerate_charge_mismatch: bool,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> bool {
    // RDKit✔️✔️:     case 75:  // Re
    // RDKit✔️✔️:       if (tolerateChargeMismatch) {
    // RDKit✔️✔️:         if (atomKey == "Re6") {
    // RDKit✔️✔️:           atomKey = "Re6+5";
    // RDKit✔️✔️:         } else if (atomKey == "Re3") {
    // RDKit✔️✔️:           atomKey = "Re3+7";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       BOOST_LOG(rdErrorLog)
    // RDKit✔️✔️:           << "UFFTYPER: Unrecognized charge state for atom: " << atom->getIdx()
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       break;
    if atom.atomic_number() != 75 {
        return false;
    }

    if tolerate_charge_mismatch {
        if atom_key.as_bytes() == b"Re6" {
            atom_key.clear();
            atom_key.extend_bytes(b"Re6+5");
        } else if atom_key.as_bytes() == b"Re3" {
            atom_key.clear();
            atom_key.extend_bytes(b"Re3+7");
        }
    }

    diagnostics.push(UffTypingDiagnostic {
        atom_id: Some(atom.id()),
        kind: UffTypingDiagnosticKind::Error,
        message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
    });
    true
}

fn append_lanthanide_charge_flag(
    atom: &Atom,
    total_valence: i32,
    atom_key: &mut cosmolkit_model::PropertyText,
    tolerate_charge_mismatch: bool,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> bool {
    // RDKit✔️✔️:   // lanthanides
    // RDKit✔️✔️:   if (atom->getAtomicNum() >= 57 && atom->getAtomicNum() <= 71) {
    // RDKit✔️✔️:     switch (totalValence) {
    // RDKit✔️✔️:       case 6:
    // RDKit✔️✔️:         atomKey += "+3";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         if (tolerateChargeMismatch) {
    // RDKit✔️✔️:           atomKey += "+3";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         BOOST_LOG(rdErrorLog)
    // RDKit✔️✔️:             << "UFFTYPER: Unrecognized charge state for atom: "
    // RDKit✔️✔️:             << atom->getIdx() << std::endl;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if !(atom.atomic_number() >= 57 && atom.atomic_number() <= 71) {
        return false;
    }

    match total_valence {
        6 => atom_key.extend_bytes(b"+3"),
        _ => {
            if tolerate_charge_mismatch {
                atom_key.extend_bytes(b"+3");
            }
            diagnostics.push(UffTypingDiagnostic {
                atom_id: Some(atom.id()),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
            });
        }
    }
    true
}

fn add_atom_charge_flags(
    atom: &Atom,
    total_valence: i32,
    atom_key: &mut cosmolkit_model::PropertyText,
    tolerate_charge_mismatch: bool,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) {
    // RDKit❗✔️: void addAtomChargeFlags(const Atom *atom, std::string &atomKey,
    // RDKit❗✔️:                         bool tolerateChargeMismatch) {
    // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
    // The borrowed Atom reference makes the source non-null precondition
    // structural. `total_valence` is the already prepared source-aligned value.
    // RDKit❗✔️:   int totalValence = atom->getTotalValence();
    // RDKit❗✔️:   int fc = atom->getFormalCharge();
    let _source_formal_charge = i32::from(atom.formal_charge());
    // RDKit❗✔️:   // FIX: come up with some way of handling metals here
    // RDKit❗✔️:   switch (atom->getAtomicNum()) {
    match atom.atomic_number() {
        // RDKit❗✔️:     case 29:  // Cu
        // RDKit❗✔️:     case 47:  // Ag
        29 | 47 => {
            append_fixed_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 4:   // Be
        // RDKit❗✔️:     case 20:  // Ca
        // RDKit❗✔️:     case 25:  // Mn
        // RDKit❗✔️:     case 26:  // Fe
        // RDKit❗✔️:     case 28:  // Ni
        // RDKit❗✔️:     case 46:  // Pd
        // RDKit❗✔️:     case 78:  // Pt
        4 | 20 | 25 | 26 | 28 | 46 | 78 => {
            append_fixed_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 21:   // Sc
        // RDKit❗✔️:     case 24:   // Cr
        // RDKit❗✔️:     case 27:   // Co
        // RDKit❗✔️:     case 79:   // Au
        // RDKit❗✔️:     case 89:   // Ac
        // RDKit❗✔️:     case 96:   // Cm
        // RDKit❗✔️:     case 103:  // Lr/Lw
        21 | 24 | 27 | 79 | 89 | 96..=103 => {
            append_fixed_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 2:   // He
        // RDKit❗✔️:     case 18:  // Ar
        // RDKit❗✔️:     case 22:  // Ti
        // RDKit❗✔️:     case 36:  // Kr
        // RDKit❗✔️:     case 54:  // Xe
        // RDKit❗✔️:     case 90:  // Th
        // RDKit❗✔️:     case 95:  // Am
        2 | 18 | 22 | 36 | 54 | 90..=95 => {
            append_fixed_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 23:  // V
        // RDKit❗✔️:     case 41:  // Nb
        // RDKit❗✔️:     case 43:  // Tc
        // RDKit❗✔️:     case 73:  // Ta
        23 | 41 | 43 | 73 => {
            append_fixed_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 42:  // Mo
        42 => {
            append_fixed_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 12:  // Mg
        12 => {
            append_valence_only_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 13:  // Al
        // RDKit❗✔️:     case 14:  // Si
        13 | 14 => {
            check_unsuffixed_main_group_charge(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 15:  // P
        15 => {
            append_phosphorus_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 16:  // S
        16 => {
            append_sulfur_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 30:  // Zn
        30 => {
            append_valence_only_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 31:  // Ga
        31 => {
            append_valence_only_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 33:  // As
        33 => {
            append_valence_only_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 34:  // Se
        34 => {
            append_valence_only_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 48:  // Cd
        48 => {
            append_valence_only_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 49:  // In
        49 => {
            append_valence_only_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 51:  // Sb
        51 => {
            append_valence_only_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 52:  // Te
        52 => {
            append_valence_only_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 75:  // Re
        75 => {
            rewrite_rhenium_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:     case 80:  // Hg
        // RDKit❗✔️:     case 81:  // Tl
        // RDKit❗✔️:     case 82:  // Pb
        // RDKit❗✔️:     case 83:  // Bi
        // RDKit❗✔️:     case 84:  // Po
        80..=84 => {
            append_valence_only_charge_flag(
                atom,
                total_valence,
                atom_key,
                tolerate_charge_mismatch,
                diagnostics,
            );
        }
        // RDKit❗✔️:   }
        _ => {}
    }
    // RDKit❗✔️:   // lanthanides
    // RDKit❗✔️:   if (atom->getAtomicNum() >= 57 && atom->getAtomicNum() <= 71) {
    append_lanthanide_charge_flag(
        atom,
        total_valence,
        atom_key,
        tolerate_charge_mismatch,
        diagnostics,
    );
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
}

fn atom_label_prefix(
    atom: &Atom,
    hybridization: Hybridization,
    mut atom_has_conjugated_bond: impl FnMut() -> bool,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<cosmolkit_model::PropertyText, UffTypingError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::getAtomLabel prefix/hybridization (AtomTyper.cpp:415-498)
    // RDKit❗✔️: std::string getAtomLabel(const Atom *atom) {
    // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
    // RDKit❗✔️:   int atNum = atom->getAtomicNum();
    // RDKit❗✔️:   std::string atomKey = atom->getSymbol();
    // RDKit❗✔️:   if (atomKey.size() == 1) {
    // RDKit❗✔️:     atomKey += '_';
    // RDKit❗✔️:   }
    // RDKit❗✔️:   PeriodicTable *table = PeriodicTable::getTable();
    // RDKit❗✔️:   // FIX: handle main group/organometallic cases better:
    // RDKit❗✔️:   if (atNum) {
    // RDKit❗✔️:     // do not do hybridization on alkali metals or halogens:
    // RDKit❗✔️:     if (table->getDefaultValence(atNum) == -1 ||
    // RDKit❗✔️:         (table->getNouterElecs(atNum) != 1 &&
    // RDKit❗✔️:          table->getNouterElecs(atNum) != 7)) {
    // RDKit❗✔️:       switch (atom->getAtomicNum()) {
    // RDKit❗✔️:         case 12:
    // RDKit❗✔️:         case 13:
    // RDKit❗✔️:         case 14:
    // RDKit❗✔️:         case 15:
    // RDKit❗✔️:         case 50:
    // RDKit❗✔️:         case 51:
    // RDKit❗✔️:         case 52:
    // RDKit❗✔️:         case 81:
    // RDKit❗✔️:         case 82:
    // RDKit❗✔️:         case 83:
    // RDKit❗✔️:         case 84:
    // RDKit❗✔️:           atomKey += '3';
    // RDKit❗✔️:           if (atom->getHybridization() != Atom::SP3) {
    // RDKit❗✔️:             BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:                 << "UFFTYPER: Warning: hybridization set to SP3 for atom "
    // RDKit❗✔️:                 << atom->getIdx() << std::endl;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case 80:
    // RDKit❗✔️:           atomKey += '1';
    // RDKit❗✔️:           if (atom->getHybridization() != Atom::SP) {
    // RDKit❗✔️:             BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:                 << "UFFTYPER: Warning: hybridization set to SP for atom "
    // RDKit❗✔️:                 << atom->getIdx() << std::endl;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         default:
    // RDKit❗✔️:           switch (atom->getHybridization()) {
    // RDKit❗✔️:             case Atom::S:
    // RDKit❗✔️:               // don't need to do anything here
    // RDKit❗✔️:               break;
    // RDKit❗✔️:             case Atom::SP:
    // RDKit❗✔️:               atomKey += '1';
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             case Atom::SP2:
    // RDKit❗✔️:               if ((atom->getIsAromatic() ||
    // RDKit❗✔️:                    MolOps::atomHasConjugatedBond(atom)) &&
    // RDKit❗✔️:                   (atNum == 6 || atNum == 7 || atNum == 8 || atNum == 16)) {
    // RDKit❗✔️:                 atomKey += 'R';
    // RDKit❗✔️:               } else {
    // RDKit❗✔️:                 atomKey += '2';
    // RDKit❗✔️:               }
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             case Atom::SP3:
    // RDKit❗✔️:               atomKey += '3';
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             case Atom::SP2D:
    // RDKit❗✔️:               atomKey += '4';
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             case Atom::SP3D:
    // RDKit❗✔️:               atomKey += '5';
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             case Atom::SP3D2:
    // RDKit❗✔️:               atomKey += '6';
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             default:
    // RDKit❗✔️:               BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:                   << "UFFTYPER: Unrecognized hybridization for atom: "
    // RDKit❗✔️:                   << atom->getIdx() << std::endl;
    // RDKit❗✔️:           }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // special cases by element type:
    // RDKit❗✔️: }
    //
    // BEGIN RDKIT CPP FUNCTION Atom::getSymbol (Atom.cpp:263-270)
    // RDKit❗✔️: std::string Atom::getSymbol() const {
    // RDKit❗✔️:   std::string res;
    // RDKit❗✔️:   // handle dummies differently:
    // RDKit❗✔️:   if (d_atomicNum != 0 ||
    // RDKit❗✔️:       !getPropIfPresent<std::string>(common_properties::dummyLabel, res)) {
    // RDKit❗✔️:     res = PeriodicTable::getTable()->getElementSymbol(d_atomicNum);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Atom::getSymbol
    // BEGIN RDKIT CPP FUNCTION Dict::getValIfPresent(std::string&) (Dict.h:267-274)
    // RDKit❗✔️: bool getValIfPresent(const std::string_view what, std::string &res) const {
    // RDKit❗✔️:   for (const auto &i : _data) {
    // RDKit❗✔️:     if (i.key == what) {
    // RDKit❗✔️:       rdvalue_tostring(i.val, res);
    // RDKit❗✔️:       return true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Dict::getValIfPresent(std::string&)
    // RDKit::Atom::getSymbol uses the model's string dummyLabel when the
    // atomic number is zero. Its typed scalar spelling is delegated to the
    // core rdvalue_tostring port, which retains the conversion error type.
    // Complexity review — RDKit❗✔️: this performs one property lookup and
    // one value conversion. The source scans its ordered Dict rows linearly;
    // the model's ordered PropertyStore uses one BTreeMap lookup, with no
    // intermediate property clone or additional scan. The converted symbol
    // requires the output byte buffer allocation that the source also performs.
    // Canonical PropertyText retains NUL and opaque dummy-label bytes; the
    // source suffix branches append bytes without any Unicode conversion.
    // Cached state keeps the source's branch-lazy incident-bond scan here;
    // supplied rows use an O(1) indexed conjugation-flag read.
    let atomic_number = atom.atomic_number();
    let atom_symbol = if atomic_number == 0 {
        match atom.prop("dummyLabel") {
            Some(dummy_label) => {
                property_value_to_string(dummy_label).map_err(UffTypingError::CorePropertyString)?
            }
            None => rdkit_element_symbol(atomic_number)
                .map_err(UffTypingError::CoreValence)?
                .into(),
        }
    } else {
        rdkit_element_symbol(atomic_number)
            .map_err(UffTypingError::CoreValence)?
            .into()
    };
    let mut atom_key = atom_symbol;
    if atom_key.len() == 1 {
        atom_key.push_byte(b'_');
    }

    if atomic_number != 0 {
        let default_valence =
            rdkit_default_valence(atomic_number).map_err(UffTypingError::CoreValence)?;
        // Preserve the source short-circuit and its two possible outer-electron
        // getter calls: the second lookup follows only when the first is not 1.
        if default_valence == -1
            || (periodic_table_outer_electrons(atomic_number)
                .map_err(UffTypingError::CoreValence)?
                != 1
                && periodic_table_outer_electrons(atomic_number)
                    .map_err(UffTypingError::CoreValence)?
                    != 7)
        {
            match atomic_number {
                12 | 13 | 14 | 15 | 50 | 51 | 52 | 81 | 82 | 83 | 84 => {
                    atom_key.push_byte(b'3');
                    if hybridization != Hybridization::Sp3 {
                        diagnostics.push(UffTypingDiagnostic {
                            atom_id: Some(atom.id()),
                            kind: UffTypingDiagnosticKind::Warning,
                            message_prefix: FORCED_SP3_HYBRIDIZATION_WARNING_MESSAGE,
                        });
                    }
                }
                80 => {
                    atom_key.push_byte(b'1');
                    if hybridization != Hybridization::Sp {
                        diagnostics.push(UffTypingDiagnostic {
                            atom_id: Some(atom.id()),
                            kind: UffTypingDiagnosticKind::Warning,
                            message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
                        });
                    }
                }
                _ => match hybridization {
                    Hybridization::S => {}
                    Hybridization::Sp => atom_key.push_byte(b'1'),
                    Hybridization::Sp2 => {
                        // This borrowed callable is evaluated at the source
                        // expression; aromatic atoms short-circuit its bond scan.
                        if (atom.is_aromatic() || atom_has_conjugated_bond())
                            && matches!(atomic_number, 6 | 7 | 8 | 16)
                        {
                            atom_key.push_byte(b'R');
                        } else {
                            atom_key.push_byte(b'2');
                        }
                    }
                    Hybridization::Sp3 => atom_key.push_byte(b'3'),
                    Hybridization::Sp2d => atom_key.push_byte(b'4'),
                    Hybridization::Sp3d => atom_key.push_byte(b'5'),
                    Hybridization::Sp3d2 => atom_key.push_byte(b'6'),
                    Hybridization::Unspecified | Hybridization::Other => {
                        diagnostics.push(UffTypingDiagnostic {
                            atom_id: Some(atom.id()),
                            kind: UffTypingDiagnosticKind::Error,
                            message_prefix: UNRECOGNIZED_HYBRIDIZATION_MESSAGE,
                        });
                    }
                },
            }
        }
    }
    Ok(atom_key)
    // END RDKIT CPP FUNCTION UFF::Tools::getAtomLabel prefix/hybridization
}

fn get_atom_label(
    atom: &Atom,
    total_valence: i32,
    hybridization: Hybridization,
    atom_has_conjugated_bond: impl FnMut() -> bool,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<cosmolkit_model::PropertyText, UffTypingError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::getAtomLabel (AtomTyper.cpp:415-503)
    // RDKit❗✔️: std::string getAtomLabel(const Atom *atom) {
    // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
    // RDKit❗✔️:   int atNum = atom->getAtomicNum();
    // RDKit❗✔️:   std::string atomKey = atom->getSymbol();
    // RDKit❗✔️:   if (atomKey.size() == 1) {
    // RDKit❗✔️:     atomKey += '_';
    // RDKit❗✔️:   }
    // RDKit❗✔️:   PeriodicTable *table = PeriodicTable::getTable();
    // RDKit❗✔️:
    // RDKit❗✔️:   // FIX: handle main group/organometallic cases better:
    // RDKit❗✔️:   if (atNum) {
    // RDKit❗✔️:     // do not do hybridization on alkali metals or halogens:
    // RDKit❗✔️:     if (table->getDefaultValence(atNum) == -1 ||
    // RDKit❗✔️:         (table->getNouterElecs(atNum) != 1 &&
    // RDKit❗✔️:          table->getNouterElecs(atNum) != 7)) {
    // RDKit❗✔️:       switch (atom->getAtomicNum()) {
    // RDKit❗✔️:         case 12:
    // RDKit❗✔️:         case 13:
    // RDKit❗✔️:         case 14:
    // RDKit❗✔️:         case 15:
    // RDKit❗✔️:         case 50:
    // RDKit❗✔️:         case 51:
    // RDKit❗✔️:         case 52:
    // RDKit❗✔️:         case 81:
    // RDKit❗✔️:         case 82:
    // RDKit❗✔️:         case 83:
    // RDKit❗✔️:         case 84:
    // RDKit❗✔️:           atomKey += '3';
    // RDKit❗✔️:           if (atom->getHybridization() != Atom::SP3) {
    // RDKit❗✔️:             BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:                 << "UFFTYPER: Warning: hybridization set to SP3 for atom "
    // RDKit❗✔️:                 << atom->getIdx() << std::endl;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case 80:
    // RDKit❗✔️:           atomKey += '1';
    // RDKit❗✔️:           if (atom->getHybridization() != Atom::SP) {
    // RDKit❗✔️:             BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:                 << "UFFTYPER: Warning: hybridization set to SP for atom "
    // RDKit❗✔️:                 << atom->getIdx() << std::endl;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         default:
    // RDKit❗✔️:           switch (atom->getHybridization()) {
    // RDKit❗✔️:             case Atom::S:
    // RDKit❗✔️:               // don't need to do anything here
    // RDKit❗✔️:               break;
    // RDKit❗✔️:             case Atom::SP:
    // RDKit❗✔️:               atomKey += '1';
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             case Atom::SP2:
    // RDKit❗✔️:               if ((atom->getIsAromatic() ||
    // RDKit❗✔️:                    MolOps::atomHasConjugatedBond(atom)) &&
    // RDKit❗✔️:                   (atNum == 6 || atNum == 7 || atNum == 8 || atNum == 16)) {
    // RDKit❗✔️:                 atomKey += 'R';
    // RDKit❗✔️:               } else {
    // RDKit❗✔️:                 atomKey += '2';
    // RDKit❗✔️:               }
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             case Atom::SP3:
    // RDKit❗✔️:               atomKey += '3';
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             case Atom::SP2D:
    // RDKit❗✔️:               atomKey += '4';
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             case Atom::SP3D:
    // RDKit❗✔️:               atomKey += '5';
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             case Atom::SP3D2:
    // RDKit❗✔️:               atomKey += '6';
    // RDKit❗✔️:               break;
    // RDKit❗✔️:
    // RDKit❗✔️:             default:
    // RDKit❗✔️:               BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:                   << "UFFTYPER: Unrecognized hybridization for atom: "
    // RDKit❗✔️:                   << atom->getIdx() << std::endl;
    // RDKit❗✔️:           }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // special cases by element type:
    // RDKit❗✔️:   addAtomChargeFlags(atom, atomKey);
    // RDKit❗✔️:   return atomKey;
    // RDKit❗✔️: }
    //
    // RDKit❗✔️: RDKIT_FORCEFIELDHELPERS_EXPORT void addAtomChargeFlags(
    // RDKit❗✔️:     const Atom *atom, std::string &atomKey, bool tolerateChargeMismatch = true);
    // The source header gives getAtomLabel no tolerance parameter and defaults
    // the nested charge call to true. The borrowed Atom reference preserves
    // the non-null precondition; prepared state is consumed without mutation.
    let mut atom_key =
        atom_label_prefix(atom, hybridization, atom_has_conjugated_bond, diagnostics)?;
    add_atom_charge_flags(atom, total_valence, &mut atom_key, true, diagnostics);
    Ok(atom_key)
    // END RDKIT CPP FUNCTION UFF::Tools::getAtomLabel
}

pub(super) fn get_atom_types<'params>(
    topology: &TopologyBlock,
    total_valences: &[i32],
    atom_has_conjugated_bond: &[bool],
    params: &'params ParamCollection,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<(Vec<Option<&'params AtomicParams>>, bool), UffTypingError> {
    let typing_state =
        UffAtomStateRef::supplied_rows(topology, total_valences, atom_has_conjugated_bond)?;
    get_atom_types_from_state(topology, typing_state, params, diagnostics)
}

pub(super) fn get_atom_types_from_state<'params>(
    supplied_topology: &TopologyBlock,
    typing_state: UffAtomStateRef<'_>,
    params: &'params ParamCollection,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<(Vec<Option<&'params AtomicParams>>, bool), UffTypingError> {
    get_atom_types_with_assignments(supplied_topology, typing_state, params, diagnostics, None)
}

pub(super) fn get_atom_types_with_assignments<'params>(
    supplied_topology: &TopologyBlock,
    typing_state: UffAtomStateRef<'_>,
    params: &'params ParamCollection,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    assignments: Option<(&[Hybridization], &[bool])>,
) -> Result<(Vec<Option<&'params AtomicParams>>, bool), UffTypingError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getAtomTypes (AtomTyper.cpp:507-533)
    // RDKit❗✔️: std::pair<AtomicParamVect, bool> getAtomTypes(const ROMol &mol,
    // RDKit❗✔️:                                               const std::string &) {
    let topology = match typing_state {
        UffAtomStateRef::Cached { topology, .. } => topology,
        UffAtomStateRef::SuppliedRows { .. } => supplied_topology,
    };
    let atom_count = topology.atoms.len();
    // RDKit❗✔️:   bool foundAll = true;
    let mut found_all = true;
    // RDKit❗✔️:   auto params = ParamCollection::getParams();
    // The caller owns the default `get_params("")` Arc and lends this table,
    // keeping returned references valid without rebuilding or self-owning it.
    // RDKit❗✔️:   AtomicParamVect paramVect;
    // RDKit❗✔️:   paramVect.resize(mol.getNumAtoms());
    // This single ordered loop consumes either borrowed cached rows or the
    // supplied-row adapter. Its result vector remains the source N nullable
    // parameter slots; explicit assignments are borrowed and never copied.
    let mut param_vect = Vec::with_capacity(atom_count);
    // RDKit❗✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    for (atom_index, atom) in topology.atoms.iter().enumerate() {
        // RDKit❗✔️:     const Atom *atom = mol.getAtomWithIdx(i);
        // RDKit❗✔️:     // construct the atom key:
        // RDKit❗✔️:     std::string atomKey = Tools::getAtomLabel(atom);
        let atom_key = get_atom_label(
            atom,
            typing_state.total_valence_at(atom_index),
            assignments.map_or_else(
                || atom.hybridization(),
                |(hybridizations, _)| hybridizations[atom_index],
            ),
            || match assignments {
                Some((_, conjugated_bonds)) => topology
                    .adjacency
                    .neighbors_of(atom_index)
                    .iter()
                    .any(|neighbor| conjugated_bonds[neighbor.bond.index()]),
                None => typing_state.conjugated_presence_at(atom_index),
            },
            diagnostics,
        )?;
        // RDKit❗✔️:     // ok, we've got the atom key, now get the parameters:
        // RDKit❗✔️:     const AtomicParams *theParams = (*params)(atomKey);
        let the_params = params.get(&atom_key);
        // RDKit❗✔️:     if (!theParams) {
        if the_params.is_none() {
            // RDKit❗✔️:       foundAll = false;
            found_all = false;
            // RDKit❗✔️:       BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized atom type: " << atomKey
            // RDKit❗✔️:                             << " (" << i << ")" << std::endl;
            diagnostics.push(UffTypingDiagnostic {
                atom_id: Some(atom.id()),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: UNRECOGNIZED_ATOM_TYPE_MESSAGE,
            });
            // Missing parameters are represented by this row's None slot;
            // source continues through later atom indices.
        }
        // RDKit❗✔️:     paramVect[i] = theParams;
        param_vect.push(the_params);
        // RDKit❗✔️:   }
    }
    // RDKit❗✔️:   return std::make_pair(paramVect, foundAll);
    Ok((param_vect, found_all))
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getAtomTypes
}

pub(super) fn uff_has_all_molecule_parameters(
    topology: &TopologyBlock,
    total_valences: &[i32],
    atom_has_conjugated_bond: &[bool],
    params: &ParamCollection,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<bool, UffTypingError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFFHasAllMoleculeParams (rdForceFields.cpp:107-112)
    // RDKit❗✔️: bool UFFHasAllMoleculeParams(const ROMol &mol) {
    // RDKit❗✔️:   UFF::AtomicParamVect types;
    // RDKit❗✔️:   bool foundAll;
    // RDKit❗✔️:   boost::tie(types, foundAll) = UFF::getAtomTypes(mol);
    // RDKit❗✔️:   return foundAll;
    // RDKit❗✔️: }
    // Preserve the same prepared, borrowed typing inputs and diagnostic order
    // as get_atom_types. Its source-sized nullable parameter vector is
    // materialized and discarded just as in the wrapper; typing errors remain
    // typed instead of being converted to a boolean.
    let (_, found_all) = get_atom_types(
        topology,
        total_valences,
        atom_has_conjugated_bond,
        params,
        diagnostics,
    )?;
    Ok(found_all)
    // END RDKIT CPP FUNCTION RDKit::UFFHasAllMoleculeParams
}

fn get_uff_vdw_params(
    topology: &TopologyBlock,
    idx1: usize,
    idx2: usize,
    total_valences: &[i32],
    atom_has_conjugated_bond: &[bool],
    params: &ParamCollection,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<Option<UffVdw>, UffTypingError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFVdWParams (AtomTyper.cpp:713-733)
    // RDKit✔️🔝: bool getUFFVdWParams(const ROMol &mol, unsigned int idx1, unsigned int idx2,
    // RDKit✔️🔝:                      UFFVdW &uffVdWParams) {
    // RDKit✔️🔝:   bool res = true;
    // RDKit✔️🔝:   auto params = ParamCollection::getParams();
    // The caller lends its canonical cached table and prepared scalar state;
    // this keeps the pair query local without reassigning molecule chemistry.
    // RDKit✔️🔝:   unsigned int idx[2] = {idx1, idx2};
    // RDKit✔️🔝:   AtomicParamVect paramVect(2);
    // Two ordered borrowed lookups avoid the source's heap-allocated pointer
    // vector and avoid copying AtomicParams while preserving source ordering.
    // RDKit✔️🔝:   unsigned int i;
    let atom_count = topology.atoms.len();
    if total_valences.len() != atom_count {
        return Err(UffTypingError::PreparedStateLength {
            input: UffTypingInput::TotalValence,
            expected: atom_count,
            actual: total_valences.len(),
        });
    }
    if atom_has_conjugated_bond.len() != atom_count {
        return Err(UffTypingError::PreparedStateLength {
            input: UffTypingInput::ConjugatedBondPresence,
            expected: atom_count,
            actual: atom_has_conjugated_bond.len(),
        });
    }

    // RDKit✔️🔝:   for (i = 0; res && (i < 2); ++i) {
    // RDKit✔️🔝:     const Atom *atom = mol.getAtomWithIdx(idx[i]);
    let atom0 = topology
        .atoms
        .get(idx1)
        .ok_or(UffTypingError::AtomIndexOutOfBounds {
            index: idx1,
            atom_count,
        })?;
    // RDKit✔️🔝:     std::string atomKey = Tools::getAtomLabel(atom);
    let atom_key0 = get_atom_label(
        atom0,
        total_valences[idx1],
        atom0.hybridization(),
        || atom_has_conjugated_bond[idx1],
        diagnostics,
    )?;
    // RDKit✔️🔝:     paramVect[i] = (*params)(atomKey);
    let Some(params0) = params.get(&atom_key0) else {
        // RDKit✔️🔝:     res = paramVect[i] != nullptr;
        return Ok(None);
    };

    // Keeping the second access after the first successful parameter lookup
    // preserves the source loop's short-circuit behavior.
    let atom1 = topology
        .atoms
        .get(idx2)
        .ok_or(UffTypingError::AtomIndexOutOfBounds {
            index: idx2,
            atom_count,
        })?;
    let atom_key1 = get_atom_label(
        atom1,
        total_valences[idx2],
        atom1.hybridization(),
        || atom_has_conjugated_bond[idx2],
        diagnostics,
    )?;
    let Some(params1) = params.get(&atom_key1) else {
        // RDKit✔️🔝:     res = paramVect[i] != nullptr;
        return Ok(None);
    };
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   if (res) {
    // RDKit✔️🔝:     uffVdWParams.x_ij =
    // RDKit✔️🔝:         UFF::Utils::calcNonbondedMinimum(paramVect[0], paramVect[1]);
    let x_ij = calc_nonbonded_minimum(params0, params1);
    // RDKit✔️🔝:     uffVdWParams.D_ij =
    // RDKit✔️🔝:         UFF::Utils::calcNonbondedDepth(paramVect[0], paramVect[1]);
    let d_ij = calc_nonbonded_depth(params0, params1);
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return res;
    Ok(Some(UffVdw { x_ij, d_ij }))
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFVdWParams
}

fn get_uff_bond_stretch_params(
    topology: &TopologyBlock,
    idx1: usize,
    idx2: usize,
    total_valences: &[i32],
    atom_has_conjugated_bond: &[bool],
    params: &ParamCollection,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<Option<UffBond>, UffTypingError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFBondStretchParams (AtomTyper.cpp:535-557)
    // RDKit❗✔️: bool getUFFBondStretchParams(const ROMol &mol, unsigned int idx1,
    // RDKit❗✔️:                              unsigned int idx2, UFFBond &uffBondStretchParams) {
    // RDKit❗✔️:   auto params = ParamCollection::getParams();
    // The caller owns the canonical parameter table and lends it here. This
    // avoids replacing the source flyweight lookup with a per-query table build.
    // RDKit❗✔️:   unsigned int idx[2] = {idx1, idx2};
    let endpoint_indices = [idx1, idx2];
    // RDKit❗✔️:   AtomicParamVect paramVect(2);
    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   const Bond *bond = mol.getBondBetweenAtoms(idx1, idx2);
    let atom_count = topology.atoms.len();
    for atom_index in endpoint_indices {
        if atom_index >= atom_count {
            return Err(UffTypingError::AtomIdOutOfBounds {
                atom_id: AtomId::new(atom_index),
                atom_count,
            });
        }
    }

    if total_valences.len() != atom_count {
        return Err(UffTypingError::PreparedStateLength {
            input: UffTypingInput::TotalValence,
            expected: atom_count,
            actual: total_valences.len(),
        });
    }
    if atom_has_conjugated_bond.len() != atom_count {
        return Err(UffTypingError::PreparedStateLength {
            input: UffTypingInput::ConjugatedBondPresence,
            expected: atom_count,
            actual: atom_has_conjugated_bond.len(),
        });
    }

    let Some(neighbor) = topology
        .adjacency
        .neighbors_of(idx1)
        .iter()
        .find(|neighbor| neighbor.atom_index == idx2)
    else {
        // RDKit❗✔️:   bool res = bond != nullptr;
        // A missing edge is the source's ordinary false result.
        return Ok(None);
    };
    let bond_id = neighbor.bond;
    let Some(bond) = topology.bonds.get(bond_id.index()) else {
        return Err(UffTypingError::BondIdOutOfBounds {
            bond_id,
            bond_count: topology.bonds.len(),
        });
    };

    // RDKit❗✔️:   for (i = 0; res && (i < 2); ++i) {
    // Borrow endpoint parameters directly so this query does not allocate the
    // source's two-pointer vector or copy either AtomicParams value.
    // RDKit❗✔️:     const Atom *atom = mol.getAtomWithIdx(idx[i]);
    let atom0 = &topology.atoms[idx1];
    // RDKit❗✔️:     std::string atomKey = Tools::getAtomLabel(atom);
    let atom_key0 = get_atom_label(
        atom0,
        total_valences[idx1],
        atom0.hybridization(),
        || atom_has_conjugated_bond[idx1],
        diagnostics,
    )?;
    // RDKit❗✔️:     paramVect[i] = (*params)(atomKey);
    let Some(params0) = params.get(&atom_key0) else {
        // RDKit❗✔️:     res = paramVect[i] != nullptr;
        return Ok(None);
    };

    let atom1 = &topology.atoms[idx2];
    let atom_key1 = get_atom_label(
        atom1,
        total_valences[idx2],
        atom1.hybridization(),
        || atom_has_conjugated_bond[idx2],
        diagnostics,
    )?;
    let Some(params1) = params.get(&atom_key1) else {
        // RDKit❗✔️:     res = paramVect[i] != nullptr;
        return Ok(None);
    };
    // RDKit❗✔️:   }

    // RDKit❗✔️:   if (res) {
    // RDKit❗✔️:     double bondOrder = bond->getBondTypeAsDouble();
    let bond_order =
        cosmolkit_core::bond_type_as_double(bond.order()).map_err(UffTypingError::CoreValence)?;
    // RDKit❗✔️:     uffBondStretchParams.r0 =
    // RDKit❗✔️:         UFF::Utils::calcBondRestLength(bondOrder, paramVect[0], paramVect[1]);
    let r0 =
        calc_bond_rest_length(bond_order, params0, params1).map_err(UffTypingError::BondMath)?;
    // RDKit❗✔️:     uffBondStretchParams.kb = UFF::Utils::calcBondForceConstant(
    // RDKit❗✔️:         uffBondStretchParams.r0, paramVect[0], paramVect[1]);
    let kb = calc_bond_force_constant(r0, params0, params1);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    Ok(Some(UffBond { kb, r0 }))
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFBondStretchParams
}

fn get_uff_angle_bend_params(
    topology: &TopologyBlock,
    idx1: usize,
    idx2: usize,
    idx3: usize,
    total_valences: &[i32],
    atom_has_conjugated_bond: &[bool],
    params: &ParamCollection,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<Option<UffAngle>, UffTypingError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFAngleBendParams (AtomTyper.cpp:559-589)
    // RDKit❗✔️: bool getUFFAngleBendParams(const ROMol &mol, unsigned int idx1,
    // RDKit❗✔️:                            unsigned int idx2, unsigned int idx3,
    // RDKit❗✔️:                            UFFAngle &uffAngleBendParams) {
    // RDKit❗✔️:   auto params = ParamCollection::getParams();
    // The caller owns the canonical table and lends it here, avoiding a
    // per-angle table lookup/build. Source query order remains row-local.
    // RDKit❗✔️:   unsigned int idx[3] = {idx1, idx2, idx3};
    let atom_indices = [idx1, idx2, idx3];
    // RDKit❗✔️:   AtomicParamVect paramVect(3);
    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   const Bond *bond[2];
    // Fixed stack slots keep borrowed parameters/bonds in source order without
    // a heap-backed parameter vector or copying AtomicParams records.
    let mut param_vect: [Option<&AtomicParams>; 3] = [None, None, None];
    let mut bonds: [Option<&Bond>; 2] = [None, None];
    // RDKit❗✔️:   bool res = true;
    let atom_count = topology.atoms.len();
    if total_valences.len() != atom_count {
        return Err(UffTypingError::PreparedStateLength {
            input: UffTypingInput::TotalValence,
            expected: atom_count,
            actual: total_valences.len(),
        });
    }
    if atom_has_conjugated_bond.len() != atom_count {
        return Err(UffTypingError::PreparedStateLength {
            input: UffTypingInput::ConjugatedBondPresence,
            expected: atom_count,
            actual: atom_has_conjugated_bond.len(),
        });
    }

    // RDKit❗✔️:   for (i = 0; res && (i < 3); ++i) {
    for i in 0..3 {
        // RDKit❗✔️:     if (i < 2) {
        if i < 2 {
            let atom_index = atom_indices[i];
            let next_atom_index = atom_indices[i + 1];
            // ROMol::getBondBetweenAtoms range-checks both IDs at each query.
            for endpoint_index in [atom_index, next_atom_index] {
                if endpoint_index >= atom_count {
                    return Err(UffTypingError::AtomIdOutOfBounds {
                        atom_id: AtomId::new(endpoint_index),
                        atom_count,
                    });
                }
            }
            // RDKit❗✔️:       bond[i] = mol.getBondBetweenAtoms(idx[i], idx[i + 1]);
            let Some(neighbor) = topology
                .adjacency
                .neighbors_of(atom_index)
                .iter()
                .find(|neighbor| neighbor.atom_index == next_atom_index)
            else {
                // RDKit❗✔️:       res = bond[i] != nullptr;
                return Ok(None);
            };
            let Some(bond) = topology.bonds.get(neighbor.bond.index()) else {
                return Err(UffTypingError::BondIdOutOfBounds {
                    bond_id: neighbor.bond,
                    bond_count: topology.bonds.len(),
                });
            };
            bonds[i] = Some(bond);
            // RDKit❗✔️:       res = bond[i] != nullptr;
        }
        // RDKit❗✔️:     if (res) {
        let atom_index = atom_indices[i];
        let Some(atom) = topology.atoms.get(atom_index) else {
            return Err(UffTypingError::AtomIdOutOfBounds {
                atom_id: AtomId::new(atom_index),
                atom_count,
            });
        };
        // RDKit❗✔️:       const Atom *atom = mol.getAtomWithIdx(idx[i]);
        // RDKit❗✔️:       std::string atomKey = Tools::getAtomLabel(atom);
        let atom_key = get_atom_label(
            atom,
            total_valences[atom_index],
            atom.hybridization(),
            || atom_has_conjugated_bond[atom_index],
            diagnostics,
        )?;
        // RDKit❗✔️:       paramVect[i] = (*params)(atomKey);
        let Some(atom_params) = params.get(&atom_key) else {
            // RDKit❗✔️:       res = paramVect[i] != nullptr;
            return Ok(None);
        };
        param_vect[i] = Some(atom_params);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
    }

    // RDKit❗✔️:   if (res) {
    let [Some(bond12), Some(bond23)] = bonds else {
        return Ok(None);
    };
    let [Some(params0), Some(params1), Some(params2)] = param_vect else {
        return Ok(None);
    };
    // RDKit❗✔️:     double bondOrder12 = bond[0]->getBondTypeAsDouble();
    let bond_order12 =
        cosmolkit_core::bond_type_as_double(bond12.order()).map_err(UffTypingError::CoreValence)?;
    // RDKit❗✔️:     double bondOrder23 = bond[1]->getBondTypeAsDouble();
    let bond_order23 =
        cosmolkit_core::bond_type_as_double(bond23.order()).map_err(UffTypingError::CoreValence)?;
    // RDKit❗✔️:     uffAngleBendParams.theta0 = RAD2DEG * paramVect[1]->theta0;
    let theta0 = RAD2DEG * params1.theta0;
    // RDKit❗✔️:     uffAngleBendParams.ka = UFF::Utils::calcAngleForceConstant(
    // RDKit❗✔️:         paramVect[1]->theta0, bondOrder12, bondOrder23, paramVect[0],
    // RDKit❗✔️:         paramVect[1], paramVect[2]);
    let ka = calc_angle_force_constant(
        params1.theta0,
        bond_order12,
        bond_order23,
        params0,
        params1,
        params2,
    )
    .map_err(UffTypingError::BondMath)?;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFAngleBendParams
    Ok(Some(UffAngle { ka, theta0 }))
}

fn torsion_sp3_amplitude(
    bond_order: f64,
    at_num0: i32,
    at_num1: i32,
    params0: &AtomicParams,
    params1: &AtomicParams,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFTorsionParams SP3/SP3 amplitude (AtomTyper.cpp:625-640)
    // RDKit✔️✔️:     if ((hyb[0] == RDKit::Atom::SP3) && (hyb[1] == RDKit::Atom::SP3)) {
    // RDKit✔️✔️:       // general case:
    // RDKit✔️✔️:       uffTorsionParams.V = sqrt(paramVect[0]->V1 * paramVect[1]->V1);
    let mut amplitude = (params0.v1 * params1.v1).sqrt();
    // RDKit✔️✔️:       // special case for single bonds between group 6 elements:
    // RDKit✔️✔️:       if (((int)(bondOrder * 10) == 10) && UFF::Utils::isInGroup6(atNum[0]) &&
    // RDKit✔️✔️:           UFF::Utils::isInGroup6(atNum[1])) {
    // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::isInGroup6 (TorsionAngle.cpp:37-39)
    // RDKit✔️✔️: bool isInGroup6(int num) {
    // RDKit✔️✔️:   return (num == 8 || num == 16 || num == 34 || num == 52 || num == 84);
    // RDKit✔️✔️: }
    // END RDKIT CPP HELPER ForceFields::UFF::Utils::isInGroup6
    if (bond_order * 10.0) as i32 == 10 && is_in_group6(at_num0) && is_in_group6(at_num1) {
        // RDKit✔️✔️:         double V2 = 6.8;
        // RDKit✔️✔️:         double V3 = 6.8;
        let mut v2: f64 = 6.8;
        let mut v3: f64 = 6.8;
        // RDKit✔️✔️:         if (atNum[0] == 8) {
        // RDKit✔️✔️:           V2 = 2.0;
        // RDKit✔️✔️:         }
        if at_num0 == 8 {
            v2 = 2.0;
        }
        // RDKit✔️✔️:         if (atNum[1] == 8) {
        // RDKit✔️✔️:           V3 = 2.0;
        // RDKit✔️✔️:         }
        if at_num1 == 8 {
            v3 = 2.0;
        }
        // RDKit✔️✔️:         uffTorsionParams.V = sqrt(V2 * V3);
        amplitude = (v2 * v3).sqrt();
    }
    // RDKit✔️✔️:       }
    amplitude
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFTorsionParams SP3/SP3 amplitude
}

fn torsion_mixed_amplitude(
    bond_order: f64,
    at_num0: i32,
    at_num1: i32,
    hyb0: Hybridization,
    hyb1: Hybridization,
    params0: &AtomicParams,
    params1: &AtomicParams,
    has_sp2: bool,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFTorsionParams mixed SP2/SP3 amplitude (AtomTyper.cpp:645-662)
    // RDKit✔️✔️:       // SP2 - SP3,  this is, by default, independent of atom type in UFF:
    // RDKit✔️✔️:       uffTorsionParams.V = 1.0;
    let mut amplitude = 1.0;
    // RDKit✔️✔️:       if ((int)(bondOrder * 10) == 10) {
    if (bond_order * 10.0) as i32 == 10 {
        // RDKit✔️✔️:         // special case between group 6 sp3 and non-group 6 sp2:
        // RDKit✔️✔️:         if (((hyb[0] == RDKit::Atom::SP3) && UFF::Utils::isInGroup6(atNum[0]) &&
        // RDKit✔️✔️:              (!UFF::Utils::isInGroup6(atNum[1]))) ||
        // RDKit✔️✔️:             ((hyb[1] == RDKit::Atom::SP3) && UFF::Utils::isInGroup6(atNum[1]) &&
        // RDKit✔️✔️:              (!UFF::Utils::isInGroup6(atNum[0])))) {
        // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::isInGroup6 (TorsionAngle.cpp:37-40)
        // RDKit✔️✔️: bool isInGroup6(int num) {
        // RDKit✔️✔️:   return (num == 8 || num == 16 || num == 34 || num == 52 || num == 84);
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER ForceFields::UFF::Utils::isInGroup6
        if (hyb0 == Hybridization::Sp3 && is_in_group6(at_num0) && !is_in_group6(at_num1))
            || (hyb1 == Hybridization::Sp3 && is_in_group6(at_num1) && !is_in_group6(at_num0))
        {
            // RDKit✔️✔️:           uffTorsionParams.V =
            // RDKit✔️✔️:               UFF::Utils::equation17(bondOrder, paramVect[0], paramVect[1]);
            // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::equation17 (TorsionAngle.cpp:42-47)
            // RDKit✔️✔️: double equation17(double bondOrder23, const AtomicParams *at2Params,
            // RDKit✔️✔️:                   const AtomicParams *at3Params) {
            // RDKit✔️✔️:   return 5. * sqrt(at2Params->U1 * at3Params->U1) *
            // RDKit✔️✔️:          (1. + 4.18 * log(bondOrder23));
            // RDKit✔️✔️: }
            // END RDKIT CPP HELPER ForceFields::UFF::Utils::equation17
            amplitude = equation17(bond_order, params0, params1);
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:         // special case for sp3 - sp2 - sp2
        // RDKit✔️✔️:         // (i.e. the sp2 has another sp2 neighbor, like propene)
        // RDKit✔️✔️:         else if (hasSP2) {
        // RDKit✔️✔️:           uffTorsionParams.V = 2.0;
        } else if has_sp2 {
            amplitude = 2.0;
        }
        // RDKit✔️✔️:         }
    }
    // RDKit✔️✔️:     }
    amplitude
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFTorsionParams mixed SP2/SP3 amplitude
}

fn get_uff_torsion_params(
    topology: &TopologyBlock,
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    total_valences: &[i32],
    atom_has_conjugated_bond: &[bool],
    params: &ParamCollection,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<Option<UffTor>, UffTypingError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFTorsionParams (AtomTyper.cpp:591-665)
    // RDKit✔️🔝: bool getUFFTorsionParams(const ROMol &mol, unsigned int idx1, unsigned int idx2,
    // RDKit✔️🔝:                          unsigned int idx3, unsigned int idx4,
    // RDKit✔️🔝:                          UFFTor &uffTorsionParams) {
    // RDKit✔️🔝:   auto params = ParamCollection::getParams();
    // The private caller lends the already cached default table for this lookup.
    // RDKit✔️🔝:   unsigned int idx[4] = {idx1, idx2, idx3, idx4};
    let idx = [idx1, idx2, idx3, idx4];
    // RDKit✔️🔝:   AtomicParamVect paramVect(2);
    // This fixed two-slot array preserves source slot order and removes the
    // source vector's heap allocation without changing either lookup.
    let mut param_vect: [Option<&AtomicParams>; 2] = [None, None];
    // RDKit✔️🔝:   unsigned int i;
    // RDKit✔️🔝:   const Bond *bond = mol.getBondBetweenAtoms(idx2, idx3);
    let atom_count = topology.atoms.len();
    if total_valences.len() != atom_count {
        return Err(UffTypingError::PreparedStateLength {
            input: UffTypingInput::TotalValence,
            expected: atom_count,
            actual: total_valences.len(),
        });
    }
    if atom_has_conjugated_bond.len() != atom_count {
        return Err(UffTypingError::PreparedStateLength {
            input: UffTypingInput::ConjugatedBondPresence,
            expected: atom_count,
            actual: atom_has_conjugated_bond.len(),
        });
    }
    for atom_index in [idx2, idx3] {
        if atom_index >= atom_count {
            return Err(UffTypingError::AtomIndexOutOfBounds {
                index: atom_index,
                atom_count,
            });
        }
    }
    // BEGIN RDKIT CPP HELPER RDKit::ROMol::getBondBetweenAtoms (ROMol.cpp:338-350)
    // RDKit✔️🔝: const Bond *ROMol::getBondBetweenAtoms(unsigned int idx1,
    // RDKit✔️🔝:                                        unsigned int idx2) const {
    // RDKit✔️🔝:   URANGE_CHECK(idx1, getNumAtoms());
    // RDKit✔️🔝:   URANGE_CHECK(idx2, getNumAtoms());
    // RDKit✔️🔝:   const Bond *res = nullptr;
    // RDKit✔️🔝:
    // RDKit✔️🔝:   auto [edge, found] = boost::edge(boost::vertex(idx1, d_graph),
    // RDKit✔️🔝:                                    boost::vertex(idx2, d_graph), d_graph);
    // RDKit✔️🔝:   if (found) {
    // RDKit✔️🔝:     res = d_graph[edge];
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    // END RDKIT CPP HELPER RDKit::ROMol::getBondBetweenAtoms
    let center_bond_id = topology
        .adjacency
        .neighbors_of(idx2)
        .iter()
        .find_map(|neighbor| (neighbor.atom_index == idx3).then_some(neighbor.bond));
    // RDKit✔️🔝:   int atNum[2];
    let mut at_num = [0_i32; 2];
    // RDKit✔️🔝:   Atom::HybridizationType hyb[2];
    let mut hyb = [Hybridization::Unspecified; 2];
    // RDKit✔️🔝:   bool res = true;
    let mut res = true;
    // RDKit✔️🔝:   bool hasSP2 = false;
    let mut has_sp2 = false;
    // RDKit✔️🔝:   for (i = 0; res && (i < 4); ++i) {
    let mut i = 0;
    while res && i < 4 {
        // RDKit✔️🔝:     if (i < 3) {
        if i < 3 {
            // RDKit✔️🔝:       res = mol.getBondBetweenAtoms(idx[i], idx[i + 1]) != nullptr;
            for atom_index in [idx[i], idx[i + 1]] {
                if atom_index >= atom_count {
                    return Err(UffTypingError::AtomIndexOutOfBounds {
                        index: atom_index,
                        atom_count,
                    });
                }
            }
            res = topology
                .adjacency
                .neighbors_of(idx[i])
                .iter()
                .any(|neighbor| neighbor.atom_index == idx[i + 1]);
        }
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     const Atom *atom = mol.getAtomWithIdx(idx[i]);
        if idx[i] >= atom_count {
            return Err(UffTypingError::AtomIndexOutOfBounds {
                index: idx[i],
                atom_count,
            });
        }
        let atom = &topology.atoms[idx[i]];
        // RDKit✔️🔝:     if ((i == 1) || (i == 2)) {
        if i == 1 || i == 2 {
            // RDKit✔️🔝:       unsigned int j = i - 1;
            let j = i - 1;
            // RDKit✔️🔝:       atNum[j] = atom->getAtomicNum();
            at_num[j] = i32::from(atom.atomic_number());
            // RDKit✔️🔝:       hyb[j] = atom->getHybridization();
            hyb[j] = atom.hybridization();
            // RDKit✔️🔝:       std::string atomKey = Tools::getAtomLabel(atom);
            let atom_key = get_atom_label(
                atom,
                total_valences[idx[i]],
                atom.hybridization(),
                || atom_has_conjugated_bond[idx[i]],
                diagnostics,
            )?;
            // RDKit✔️🔝:       paramVect[j] = (*params)(atomKey);
            // BEGIN RDKIT CPP HELPER ForceFields::UFF::ParamCollection::operator() (Params.h:132-139)
            // RDKit✔️🔝:   const AtomicParams *operator()(const std::string &symbol) const {
            // RDKit✔️🔝:     std::map<std::string, AtomicParams>::const_iterator res;
            // RDKit✔️🔝:     res = d_params.find(symbol);
            // RDKit✔️🔝:     if (res != d_params.end()) {
            // RDKit✔️🔝:       return &((*res).second);
            // RDKit✔️🔝:     }
            // RDKit✔️🔝:     return nullptr;
            // RDKit✔️🔝:   }
            // END RDKIT CPP HELPER ForceFields::UFF::ParamCollection::operator()
            param_vect[j] = params.get(&atom_key);
            // RDKit✔️🔝:       res = paramVect[j] != nullptr;
            res = param_vect[j].is_some();
            // RDKit✔️🔝:     } else if (atom->getHybridization() == Atom::SP2) {
        } else if atom.hybridization() == Hybridization::Sp2 {
            // RDKit✔️🔝:       hasSP2 = true;
            has_sp2 = true;
            // RDKit✔️🔝:     }
        }
        // RDKit✔️🔝:   }
        i += 1;
    }
    // RDKit✔️🔝:   if (res) {
    if res {
        // RDKit✔️🔝:     res = (((hyb[0] == RDKit::Atom::SP2) || (hyb[0] == RDKit::Atom::SP3)) &&
        // RDKit✔️🔝:            ((hyb[1] == RDKit::Atom::SP2) || (hyb[1] == RDKit::Atom::SP3)));
        res = matches!(hyb[0], Hybridization::Sp2 | Hybridization::Sp3)
            && matches!(hyb[1], Hybridization::Sp2 | Hybridization::Sp3);
    }
    // RDKit✔️🔝:   if (res) {
    if res {
        // RDKit✔️🔝:     double bondOrder = bond->getBondTypeAsDouble();
        // BEGIN RDKIT CPP HELPER RDKit::Bond::getBondTypeAsDouble (Bond.cpp:128-185)
        // RDKit✔️🔝: double Bond::getBondTypeAsDouble() const {
        // RDKit✔️🔝:   double res;
        // RDKit✔️🔝:   switch (getBondType()) {
        // RDKit✔️🔝:     case UNSPECIFIED:
        // RDKit✔️🔝:     case IONIC:
        // RDKit✔️🔝:     case ZERO:
        // RDKit✔️🔝:       res = 0;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case SINGLE:
        // RDKit✔️🔝:       res = 1;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case DOUBLE:
        // RDKit✔️🔝:       res = 2;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case TRIPLE:
        // RDKit✔️🔝:       res = 3;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case QUADRUPLE:
        // RDKit✔️🔝:       res = 4;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case QUINTUPLE:
        // RDKit✔️🔝:       res = 5;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case HEXTUPLE:
        // RDKit✔️🔝:       res = 6;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case ONEANDAHALF:
        // RDKit✔️🔝:       res = 1.5;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case TWOANDAHALF:
        // RDKit✔️🔝:       res = 2.5;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case THREEANDAHALF:
        // RDKit✔️🔝:       res = 3.5;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case FOURANDAHALF:
        // RDKit✔️🔝:       res = 4.5;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case FIVEANDAHALF:
        // RDKit✔️🔝:       res = 5.5;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case AROMATIC:
        // RDKit✔️🔝:       res = 1.5;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     case DATIVEONE:
        // RDKit✔️🔝:       res = 1.0;
        // RDKit✔️🔝:       break;  // FIX: this should probably be different
        // RDKit✔️🔝:     case DATIVE:
        // RDKit✔️🔝:       res = 1.0;
        // RDKit✔️🔝:       break;  // FIX: again probably wrong
        // RDKit✔️🔝:     case HYDROGEN:
        // RDKit✔️🔝:       res = 0.0;
        // RDKit✔️🔝:       break;
        // RDKit✔️🔝:     default:
        // RDKit✔️🔝:       UNDER_CONSTRUCTION("Bad bond type");
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   return res;
        // RDKit✔️🔝: }
        // END RDKIT CPP HELPER RDKit::Bond::getBondTypeAsDouble
        // The core owner preserves every defined bond value and reports its typed
        // bad-bond error at the source's under-construction boundary.
        let center_bond_id = center_bond_id
            .ok_or(UffTypingError::CenterBondMissingAfterSourceSuccess { idx2, idx3 })?;
        let center_bond = topology.bonds.get(center_bond_id.index()).ok_or(
            UffTypingError::BondIdOutOfBounds {
                bond_id: center_bond_id,
                bond_count: topology.bonds.len(),
            },
        )?;
        let bond_order = cosmolkit_core::bond_type_as_double(center_bond.order())
            .map_err(UffTypingError::CoreValence)?;
        let param0 = param_vect[0].ok_or(
            UffTypingError::CentralParameterSlotMissingAfterSourceSuccess { atom_index: idx2 },
        )?;
        let param1 = param_vect[1].ok_or(
            UffTypingError::CentralParameterSlotMissingAfterSourceSuccess { atom_index: idx3 },
        )?;
        // RDKit✔️🔝:     if ((hyb[0] == RDKit::Atom::SP3) && (hyb[1] == RDKit::Atom::SP3)) {
        let v = if hyb[0] == Hybridization::Sp3 && hyb[1] == Hybridization::Sp3 {
            // RDKit✔️🔝:       // general case:
            // RDKit✔️🔝:       uffTorsionParams.V = sqrt(paramVect[0]->V1 * paramVect[1]->V1);
            // RDKit✔️🔝:       // special case for single bonds between group 6 elements:
            // RDKit✔️🔝:       if (((int)(bondOrder * 10) == 10) && UFF::Utils::isInGroup6(atNum[0]) &&
            // RDKit✔️🔝:           UFF::Utils::isInGroup6(atNum[1])) {
            // RDKit✔️🔝:         double V2 = 6.8;
            // RDKit✔️🔝:         double V3 = 6.8;
            // RDKit✔️🔝:         if (atNum[0] == 8) {
            // RDKit✔️🔝:           V2 = 2.0;
            // RDKit✔️🔝:         }
            // RDKit✔️🔝:         if (atNum[1] == 8) {
            // RDKit✔️🔝:           V3 = 2.0;
            // RDKit✔️🔝:         }
            // RDKit✔️🔝:         uffTorsionParams.V = sqrt(V2 * V3);
            // RDKit✔️🔝:       }
            torsion_sp3_amplitude(bond_order, at_num[0], at_num[1], param0, param1)
        // RDKit✔️🔝:     } else if ((hyb[0] == RDKit::Atom::SP2) && (hyb[1] == RDKit::Atom::SP2)) {
        } else if hyb[0] == Hybridization::Sp2 && hyb[1] == Hybridization::Sp2 {
            // RDKit✔️🔝:       uffTorsionParams.V =
            // RDKit✔️🔝:           UFF::Utils::equation17(bondOrder, paramVect[0], paramVect[1]);
            // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::equation17 (TorsionAngle.cpp:42-47)
            // RDKit✔️🔝: double equation17(double bondOrder23, const AtomicParams *at2Params,
            // RDKit✔️🔝:                   const AtomicParams *at3Params) {
            // RDKit✔️🔝:   return 5. * sqrt(at2Params->U1 * at3Params->U1) *
            // RDKit✔️🔝:          (1. + 4.18 * log(bondOrder23));
            // RDKit✔️🔝: }
            // END RDKIT CPP HELPER ForceFields::UFF::Utils::equation17
            equation17(bond_order, param0, param1)
        // RDKit✔️🔝:     } else {
        } else {
            // RDKit✔️🔝:       // SP2 - SP3,  this is, by default, independent of atom type in UFF:
            // RDKit✔️🔝:       uffTorsionParams.V = 1.0;
            // RDKit✔️🔝:       if ((int)(bondOrder * 10) == 10) {
            // RDKit✔️🔝:         // special case between group 6 sp3 and non-group 6 sp2:
            // RDKit✔️🔝:         if (((hyb[0] == RDKit::Atom::SP3) && UFF::Utils::isInGroup6(atNum[0]) &&
            // RDKit✔️🔝:              (!UFF::Utils::isInGroup6(atNum[1]))) ||
            // RDKit✔️🔝:             ((hyb[1] == RDKit::Atom::SP3) && UFF::Utils::isInGroup6(atNum[1]) &&
            // RDKit✔️🔝:              (!UFF::Utils::isInGroup6(atNum[0])))) {
            // RDKit✔️🔝:           uffTorsionParams.V =
            // RDKit✔️🔝:               UFF::Utils::equation17(bondOrder, paramVect[0], paramVect[1]);
            // RDKit✔️🔝:         }
            // RDKit✔️🔝:         // special case for sp3 - sp2 - sp2
            // RDKit✔️🔝:         // (i.e. the sp2 has another sp2 neighbor, like propene)
            // RDKit✔️🔝:         else if (hasSP2) {
            // RDKit✔️🔝:           uffTorsionParams.V = 2.0;
            // RDKit✔️🔝:         }
            torsion_mixed_amplitude(
                bond_order, at_num[0], at_num[1], hyb[0], hyb[1], param0, param1, has_sp2,
            )
        };
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:   }
        return Ok(Some(UffTor { v }));
    }
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFTorsionParams
    Ok(None)
}

fn get_uff_inversion_params(
    topology: &TopologyBlock,
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
) -> Result<Option<UffInv>, UffTypingError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFInversionParams (AtomTyper.cpp:667-711)
    // RDKit❗✔️: bool getUFFInversionParams(const ROMol &mol, unsigned int idx1,
    // RDKit❗✔️:                            unsigned int idx2, unsigned int idx3,
    // RDKit❗✔️:                            unsigned int idx4, UFFInv &uffInversionParams) {
    // RDKit❗✔️:   unsigned int idx[4] = {idx1, idx2, idx3, idx4};
    let idx = [idx1, idx2, idx3, idx4];
    let atom_count = topology.atoms.len();

    // BEGIN RDKIT CPP HELPER RDKit::ROMol::getBondBetweenAtoms (ROMol.cpp:338-350)
    // RDKit❗✔️: const Bond *ROMol::getBondBetweenAtoms(unsigned int idx1,
    // RDKit❗✔️:                                        unsigned int idx2) const {
    // RDKit❗✔️:   URANGE_CHECK(idx1, getNumAtoms());
    // RDKit❗✔️:   URANGE_CHECK(idx2, getNumAtoms());
    // RDKit❗✔️:   const Bond *res = nullptr;
    // RDKit❗✔️:
    // RDKit❗✔️:   auto [edge, found] = boost::edge(boost::vertex(idx1, d_graph),
    // RDKit❗✔️:                                    boost::vertex(idx2, d_graph), d_graph);
    // RDKit❗✔️:   if (found) {
    // RDKit❗✔️:     res = d_graph[edge];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDKit::ROMol::getBondBetweenAtoms
    let has_bond = |atom_a: usize, atom_b: usize| -> Result<bool, UffTypingError> {
        if atom_a >= atom_count {
            return Err(UffTypingError::AtomIndexOutOfBounds {
                index: atom_a,
                atom_count,
            });
        }
        if atom_b >= atom_count {
            return Err(UffTypingError::AtomIndexOutOfBounds {
                index: atom_b,
                atom_count,
            });
        }
        Ok(topology
            .adjacency
            .neighbors_of(atom_a)
            .iter()
            .any(|neighbor| neighbor.atom_index == atom_b))
    };
    // RDKit❗✔️:   bool res = (mol.getBondBetweenAtoms(idx1, idx2) &&
    // RDKit❗✔️:               mol.getBondBetweenAtoms(idx2, idx3) &&
    // RDKit❗✔️:               mol.getBondBetweenAtoms(idx2, idx4));
    // Keep both the source's edge order and its short-circuit range-check
    // order: a missing earlier edge prevents later endpoints from being read.
    let mut res = has_bond(idx[0], idx[1])?;
    if res {
        res = has_bond(idx[1], idx[2])?;
    }
    if res {
        res = has_bond(idx[1], idx[3])?;
    }

    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   // bool isAtom2C = false;
    // RDKit❗✔️:   bool isBoundToSP2O = false;
    let mut is_bound_to_sp2_o = false;
    // RDKit❗✔️:   unsigned int at2AtomicNum = 0;
    let mut at2_atomic_num = 0_i32;
    // RDKit❗✔️:   for (i = 0; res && (i < 4); ++i) {
    let mut i = 0;
    while res && i < 4 {
        // BEGIN RDKIT CPP HELPER RDKit::ROMol::getAtomWithIdx const (ROMol.cpp:206-214)
        // RDKit❗✔️: const Atom *ROMol::getAtomWithIdx(unsigned int idx) const {
        // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
        // RDKit❗✔️:
        // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
        // RDKit❗✔️:   const auto res = d_graph[vd];
        // RDKit❗✔️:
        // RDKit❗✔️:   POSTCONDITION(res, "");
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDKit::ROMol::getAtomWithIdx const
        // The successful source edge checks already range-check every index.
        // RDKit❗✔️:     const Atom *atom = mol.getAtomWithIdx(idx[i]);
        let atom = &topology.atoms[idx[i]];
        // RDKit❗✔️:     if (i == 1) {
        if i == 1 {
            // RDKit❗✔️:       at2AtomicNum = atom->getAtomicNum();
            // RDKit❗✔️: int getAtomicNum() const { return d_atomicNum; }
            at2_atomic_num = i32::from(atom.atomic_number());
            // RDKit❗✔️:       if (res) {
            // RDKit❗✔️:         // if the central atom is not carbon, nitrogen, oxygen,
            // RDKit❗✔️:         // phosphorous, arsenic, antimonium or bismuth, skip it
            // RDKit❗✔️:         res = (!(((at2AtomicNum != 6) && (at2AtomicNum != 7) &&
            // RDKit❗✔️:                   (at2AtomicNum != 8) && (at2AtomicNum != 15) &&
            // RDKit❗✔️:                   (at2AtomicNum != 33) && (at2AtomicNum != 51) &&
            // RDKit❗✔️:                   (at2AtomicNum != 83)) ||
            // RDKit❗✔️:                  (atom->getDegree() != 3)));
            // BEGIN RDKIT CPP HELPERS Atom::getDegree / ROMol::getAtomDegree
            // RDKit❗✔️: unsigned int Atom::getDegree() const {
            // RDKit❗✔️:   return dp_mol ? getOwningMol().getAtomDegree(this) : 0;
            // RDKit❗✔️: }
            // RDKit❗✔️: unsigned int ROMol::getAtomDegree(const Atom *at) const {
            // RDKit❗✔️:   PRECONDITION(at, "no atom");
            // RDKit❗✔️:   PRECONDITION(&at->getOwningMol() == this,
            // RDKit❗✔️:                "atom not associated with this molecule");
            // RDKit❗✔️:   return rdcast<unsigned int>(boost::out_degree(at->getIdx(), d_graph));
            // RDKit❗✔️: };
            // END RDKIT CPP HELPERS Atom::getDegree / ROMol::getAtomDegree
            res = matches!(at2_atomic_num, 6 | 7 | 8 | 15 | 33 | 51 | 83)
                && topology.adjacency.neighbors_of(idx[i]).len() == 3;
            // RDKit❗✔️:       }
            // RDKit❗✔️:       if (res) {
            // RDKit❗✔️:         // if the central atom is carbon, nitrogen or oxygen
            // RDKit❗✔️:         // but hybridization is not sp2, skip it
            // RDKit❗✔️:         res = (!(((at2AtomicNum == 6) || (at2AtomicNum == 7) ||
            // RDKit❗✔️:                   (at2AtomicNum == 8)) &&
            // RDKit❗✔️:                  (atom->getHybridization() != Atom::SP2)));
            // RDKit❗✔️: HybridizationType getHybridization() const {
            // RDKit❗✔️:   return static_cast<HybridizationType>(d_hybrid);
            // RDKit❗✔️: }
            if res {
                res = !(matches!(at2_atomic_num, 6 | 7 | 8)
                    && atom.hybridization() != Hybridization::Sp2);
            }
            // RDKit❗✔️:       }
            // RDKit❗✔️:     } else if ((atom->getAtomicNum() == 8) &&
            // RDKit❗✔️:                (atom->getHybridization() == Atom::SP2)) {
        } else if atom.atomic_number() == 8 && atom.hybridization() == Hybridization::Sp2 {
            // RDKit❗✔️:       isBoundToSP2O = true;
            is_bound_to_sp2_o = true;
            // RDKit❗✔️:     }
        }
        // RDKit❗✔️:   }
        i += 1;
    }
    // RDKit❗✔️:   if (res) {
    if res {
        // RDKit❗✔️:     isBoundToSP2O = (isBoundToSP2O && (at2AtomicNum == 6));
        is_bound_to_sp2_o = is_bound_to_sp2_o && at2_atomic_num == 6;
        // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::calcInversionCoefficientsAndForceConstant (Utils.cpp:41-86)
        // RDKit❗✔️: std::tuple<double, double, double, double>
        // RDKit❗✔️: calcInversionCoefficientsAndForceConstant(int at2AtomicNum, bool isCBoundToO) {
        // RDKit❗✔️:   double res = 0.0;
        // RDKit❗✔️:   double C0 = 0.0;
        // RDKit❗✔️:   double C1 = 0.0;
        // RDKit❗✔️:   double C2 = 0.0;
        // RDKit❗✔️:   // if the central atom is sp2 carbon, nitrogen or oxygen
        // RDKit❗✔️:   if ((at2AtomicNum == 6) || (at2AtomicNum == 7) || (at2AtomicNum == 8)) {
        // RDKit❗✔️:     C0 = 1.0;
        // RDKit❗✔️:     C1 = -1.0;
        // RDKit❗✔️:     C2 = 0.0;
        // RDKit❗✔️:     res = (isCBoundToO ? 50.0 : 6.0);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     // group 5 elements are not clearly explained in the UFF paper
        // RDKit❗✔️:     // the following code was inspired by MCCCS Towhee's ffuff.F
        // RDKit❗✔️:     double w0 = M_PI / 180.0;
        // RDKit❗✔️:     switch (at2AtomicNum) {
        // RDKit❗✔️:       // if the central atom is phosphorous
        // RDKit❗✔️:       case 15:
        // RDKit❗✔️:         w0 *= 84.4339;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       // if the central atom is arsenic
        // RDKit❗✔️:       case 33:
        // RDKit❗✔️:         w0 *= 86.9735;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       // if the central atom is antimonium
        // RDKit❗✔️:       case 51:
        // RDKit❗✔️:         w0 *= 87.7047;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       // if the central atom is bismuth
        // RDKit❗✔️:       case 83:
        // RDKit❗✔️:         w0 *= 90.0;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     C2 = 1.0;
        // RDKit❗✔️:     C1 = -4.0 * cos(w0);
        // RDKit❗✔️:     C0 = -(C1 * cos(w0) + C2 * cos(2.0 * w0));
        // RDKit❗✔️:     res = 22.0 / (C0 + C1 + C2);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   res /= 3.0;
        // RDKit❗✔️:
        // RDKit❗✔️:   return std::make_tuple(res, C0, C1, C2);
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER ForceFields::UFF::Utils::calcInversionCoefficientsAndForceConstant
        let (k, _c0, _c1, _c2) = calc_inversion_coefficients(at2_atomic_num, is_bound_to_sp2_o);
        // RDKit❗✔️:     uffInversionParams.K = std::get<0>(invCoeffForceCon);
        return Ok(Some(UffInv { k }));
    }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFInversionParams
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::super::builder::{self, UffBuilderError};
    use super::{
        AtomicParams, BondMathError, FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
        FORCED_SP3_HYBRIDIZATION_WARNING_MESSAGE, ParamCollection, UNRECOGNIZED_ATOM_TYPE_MESSAGE,
        UNRECOGNIZED_CHARGE_STATE_MESSAGE, UNRECOGNIZED_HYBRIDIZATION_MESSAGE, UffAtomStateRef,
        UffTypingDiagnostic, UffTypingDiagnosticKind, UffTypingError, UffTypingInput,
        add_atom_charge_flags, append_fixed_charge_flag, append_lanthanide_charge_flag,
        append_phosphorus_charge_flag, append_sulfur_charge_flag, append_valence_only_charge_flag,
        atom_label_prefix, check_unsuffixed_main_group_charge, get_atom_label, get_atom_types,
        get_atom_types_from_state, get_uff_angle_bend_params, get_uff_bond_stretch_params,
        get_uff_inversion_params, get_uff_torsion_params, get_uff_vdw_params,
        rewrite_rhenium_charge_flag, torsion_mixed_amplitude, torsion_sp3_amplitude,
        uff_has_all_molecule_parameters,
    };
    use crate::kernel::{
        ForceField, ForceFieldKernelError, cf3d_bld_b05_calc_energy, cf3d_bld_b05_calc_grad,
        cf3d_bld_integration_minimize,
    };
    use crate::uff::{
        angle::AngleBendContrib, bond::BondStretchContrib, inversion::InversionContrib,
        torsion::TorsionAngleContrib,
    };
    use cosmolkit_core::ValenceAssignment;
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element, Hybridization,
        PropertyValue, TopologyBlock,
    };

    #[test]
    fn uff_error_e04_typing_variants_preserve_sources_and_payloads() {
        // These manually constructed variants verify trait dispatch only;
        // source-driven BondMath and Valence failures remain covered by the
        // existing bond-query regressions below.
        let core_child = cosmolkit_core::ValenceError::BadBondType {
            bond: Some(BondId::new(73)),
            order: BondOrder::DativeRight,
        };
        let core_error = UffTypingError::CoreValence(core_child);
        assert_eq!(core_error.to_string(), format!("{core_error:?}"));
        let stored_core_child = match &core_error {
            UffTypingError::CoreValence(child) => child,
            _ => unreachable!(),
        };
        let exposed_core_child = std::error::Error::source(&core_error)
            .expect("CoreValence must expose its stored child")
            .downcast_ref::<cosmolkit_core::ValenceError>()
            .expect("CoreValence source keeps its original type");
        assert!(std::ptr::eq(stored_core_child, exposed_core_child));
        assert_eq!(
            exposed_core_child,
            &cosmolkit_core::ValenceError::BadBondType {
                bond: Some(BondId::new(73)),
                order: BondOrder::DativeRight,
            }
        );

        let bond_child = BondMathError::InvalidBondOrder { bond_order: -2.75 };
        let bond_error = UffTypingError::BondMath(bond_child);
        assert_eq!(bond_error.to_string(), format!("{bond_error:?}"));
        let stored_bond_child = match &bond_error {
            UffTypingError::BondMath(child) => child,
            _ => unreachable!(),
        };
        let exposed_bond_child = std::error::Error::source(&bond_error)
            .expect("BondMath must expose its stored child")
            .downcast_ref::<BondMathError>()
            .expect("BondMath source keeps its original type");
        assert!(std::ptr::eq(stored_bond_child, exposed_bond_child));
        match exposed_bond_child {
            BondMathError::InvalidBondOrder { bond_order } => {
                assert_eq!(bond_order.to_bits(), (-2.75_f64).to_bits());
            }
        }

        let leaves = [
            UffTypingError::PreparedStateLength {
                input: UffTypingInput::ConjugatedBondPresence,
                expected: 37,
                actual: 19,
            },
            UffTypingError::AtomIdOutOfBounds {
                atom_id: AtomId::new(31),
                atom_count: 37,
            },
            UffTypingError::AtomIndexOutOfBounds {
                index: 41,
                atom_count: 43,
            },
            UffTypingError::BondIdOutOfBounds {
                bond_id: BondId::new(47),
                bond_count: 53,
            },
            UffTypingError::CenterBondMissingAfterSourceSuccess { idx2: 59, idx3: 61 },
            UffTypingError::CentralParameterSlotMissingAfterSourceSuccess { atom_index: 67 },
        ];
        for error in leaves {
            assert!(std::error::Error::source(&error).is_none());
            assert_eq!(error.to_string(), format!("{error:?}"));
            match &error {
                UffTypingError::PreparedStateLength {
                    input,
                    expected,
                    actual,
                } => {
                    assert_eq!(*input, UffTypingInput::ConjugatedBondPresence);
                    assert_eq!((*expected, *actual), (37, 19));
                }
                UffTypingError::AtomIdOutOfBounds {
                    atom_id,
                    atom_count,
                } => assert_eq!((*atom_id, *atom_count), (AtomId::new(31), 37)),
                UffTypingError::AtomIndexOutOfBounds { index, atom_count } => {
                    assert_eq!((*index, *atom_count), (41, 43));
                }
                UffTypingError::BondIdOutOfBounds {
                    bond_id,
                    bond_count,
                } => assert_eq!((*bond_id, *bond_count), (BondId::new(47), 53)),
                UffTypingError::CenterBondMissingAfterSourceSuccess { idx2, idx3 } => {
                    assert_eq!((*idx2, *idx3), (59, 61));
                }
                UffTypingError::CentralParameterSlotMissingAfterSourceSuccess { atom_index } => {
                    assert_eq!(*atom_index, 67);
                }
                UffTypingError::CoreValence(_)
                | UffTypingError::CorePropertyString(_)
                | UffTypingError::BondMath(_) => {
                    unreachable!("nested variants were tested above")
                }
            }
        }
    }

    #[test]
    fn uff_prepare_p02_cached_state_borrows_source_rows_and_honors_no_implicit() {
        let carbon = Element::from_atomic_number(6)
            .expect("fixed source carbon atomic number is in the model range");
        let atoms = vec![
            Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(carbon).with_no_implicit(false),
            ),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(carbon).with_no_implicit(true)),
        ];
        let topology = cf3d_bld_integration_topology(atoms, &[(0, 1, BondOrder::Single, true)]);
        let assignment = ValenceAssignment {
            explicit_valence: vec![2, 3],
            implicit_hydrogens: vec![4, 5],
        };

        let state = UffAtomStateRef::cached(&topology, &assignment)
            .expect("aligned source cache and topology rows validate");
        let (borrowed_topology, borrowed_assignment) = match state {
            UffAtomStateRef::Cached {
                topology,
                assignment,
            } => (topology, assignment),
            UffAtomStateRef::SuppliedRows { .. } => {
                panic!("cached construction retains its cache-backed form")
            }
        };
        assert!(std::ptr::eq(borrowed_topology, &topology));
        assert!(std::ptr::eq(borrowed_assignment, &assignment));
        assert_eq!(state.total_valence_at(0), 6);
        assert_eq!(state.total_valence_at(1), 3);
        assert!(state.conjugated_presence_at(0));
        assert!(state.conjugated_presence_at(1));
    }

    #[test]
    fn uff_prepare_p02_supplied_rows_preserve_values_and_shape_errors() {
        let topology = topology_from_atoms(vec![fixed_atom(6, 0, 0), fixed_atom(6, 0, 1)]);
        let total_valences = [17, -3];
        let conjugated_presence = [false, true];

        let state =
            UffAtomStateRef::supplied_rows(&topology, &total_valences, &conjugated_presence)
                .expect("aligned supplied data rows retain their exact values");
        assert!(matches!(
            state,
            UffAtomStateRef::SuppliedRows {
                total_valences: rows,
                conjugated_presence: conjugation,
            } if std::ptr::eq(rows, total_valences.as_slice())
                && std::ptr::eq(conjugation, conjugated_presence.as_slice())
        ));
        assert_eq!(state.total_valence_at(0), 17);
        assert_eq!(state.total_valence_at(1), -3);
        assert!(!state.conjugated_presence_at(0));
        assert!(state.conjugated_presence_at(1));

        assert!(matches!(
            UffAtomStateRef::supplied_rows(&topology, &[9], &[true]),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: 2,
                actual: 1,
            })
        ));
        assert!(matches!(
            UffAtomStateRef::supplied_rows(&topology, &[9, 10], &[true]),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::ConjugatedBondPresence,
                expected: 2,
                actual: 1,
            })
        ));
    }

    #[test]
    fn uff_prepare_p04_empty_cached_and_supplied_typing_keep_source_empty_result() {
        let topology = topology_from_atoms(Vec::new());
        let assignment = ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        };
        let total_valences: [i32; 0] = [];
        let conjugated_presence: [bool; 0] = [];
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let states = [
            UffAtomStateRef::cached(&topology, &assignment)
                .expect("empty source cache rows validate"),
            UffAtomStateRef::supplied_rows(&topology, &total_valences, &conjugated_presence)
                .expect("empty supplied rows validate"),
        ];

        for state in states {
            super::UFF_PREP_ATOM_STATE_CONJUGATION_READS.with(|reads| reads.set(0));
            let mut diagnostics = Vec::new();
            let (slots, found_all) =
                get_atom_types_from_state(&topology, state, &params, &mut diagnostics)
                    .expect("source empty typing returns an empty nullable vector");
            assert!(slots.is_empty());
            assert!(found_all);
            assert!(diagnostics.is_empty());
            super::UFF_PREP_ATOM_STATE_CONJUGATION_READS.with(|reads| assert_eq!(reads.get(), 0));
        }
    }

    #[test]
    fn uff_prepare_p04_cached_and_supplied_typing_keep_slots_order_and_lazy_branches() {
        let atoms = vec![
            label_atom(6, 0, Hybridization::Sp2, false, None),
            label_atom(6, 1, Hybridization::Sp2, false, None),
            label_atom(6, 2, Hybridization::Sp2, false, None),
            label_atom(6, 3, Hybridization::Sp2, false, None),
            label_atom(6, 4, Hybridization::Sp2, false, None),
            label_atom(5, 5, Hybridization::Unspecified, false, None),
            label_atom(6, 6, Hybridization::Sp2, true, None),
        ];
        let topology = cf3d_bld_integration_topology(
            atoms,
            &[
                (1, 2, BondOrder::Single, true),
                (3, 4, BondOrder::Single, false),
            ],
        );
        let assignment = ValenceAssignment {
            explicit_valence: vec![2, 2, 2, 2, 2, 3, 2],
            implicit_hydrogens: vec![0; 7],
        };
        let supplied_total_valences = [2, 2, 2, 2, 2, 3, 2];
        let supplied_conjugation = [false, true, true, false, false, false, false];
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let states = [
            UffAtomStateRef::cached(&topology, &assignment)
                .expect("fixed source cache and graph validate"),
            UffAtomStateRef::supplied_rows(
                &topology,
                &supplied_total_valences,
                &supplied_conjugation,
            )
            .expect("fixed supplied rows align with the source graph"),
        ];
        let expected_symbols: [Option<&str>; 7] = [
            Some("C_2"),
            Some("C_R"),
            Some("C_R"),
            Some("C_2"),
            Some("C_2"),
            None,
            Some("C_R"),
        ];
        assert!(params.get("B_").is_none());

        for state in states {
            super::UFF_PREP_ATOM_STATE_CONJUGATION_READS.with(|reads| reads.set(0));
            let mut diagnostics = Vec::new();
            let (slots, found_all) =
                get_atom_types_from_state(&topology, state, &params, &mut diagnostics)
                    .expect("source typing retains nullable rows for every atom");
            assert_eq!(slots.len(), expected_symbols.len());
            assert!(!found_all);
            for (slot, expected_symbol) in slots.iter().zip(expected_symbols) {
                match expected_symbol {
                    Some(symbol) => {
                        let actual = slot.expect("fixed source parameter slot is present");
                        let expected = params
                            .get(symbol)
                            .expect("literal source atom label has a pinned table entry");
                        assert!(std::ptr::eq(actual, expected));
                    }
                    None => assert!(slot.is_none()),
                }
            }
            assert_eq!(
                diagnostics,
                [
                    UffTypingDiagnostic {
                        atom_id: Some(AtomId::new(5)),
                        kind: UffTypingDiagnosticKind::Error,
                        message_prefix: UNRECOGNIZED_HYBRIDIZATION_MESSAGE,
                    },
                    UffTypingDiagnostic {
                        atom_id: Some(AtomId::new(5)),
                        kind: UffTypingDiagnosticKind::Error,
                        message_prefix: UNRECOGNIZED_ATOM_TYPE_MESSAGE,
                    },
                ],
            );
            super::UFF_PREP_ATOM_STATE_CONJUGATION_READS.with(|reads| assert_eq!(reads.get(), 5));
        }
    }

    const FIXED_CHARGE_GROUPS: &[(&[u8], i32, &str)] = &[
        (&[29, 47], 1, "+1"),
        (&[4, 20, 25, 26, 28, 46, 78], 2, "+2"),
        (
            &[21, 24, 27, 79, 89, 96, 97, 98, 99, 100, 101, 102, 103],
            3,
            "+3",
        ),
        (&[2, 18, 22, 36, 54, 90, 91, 92, 93, 94, 95], 4, "+4"),
        (&[23, 41, 43, 73], 5, "+5"),
        (&[42], 6, "+6"),
    ];

    const VALENCE_ONLY_CHARGE_GROUPS: &[(&[u8], i32, &str)] = &[
        (&[12, 30, 34, 48, 52, 80, 84], 2, "+2"),
        (&[31, 33, 49, 51, 81, 82, 83], 3, "+3"),
    ];

    const UNSUFFIXED_MAIN_GROUP_CHARGES: &[(u8, i32)] = &[(13, 3), (14, 4)];

    const RHENIUM_KEY_CASES: &[(&str, &str)] = &[
        ("Re6", "Re6+5"),
        ("Re3", "Re3+7"),
        ("Re", "Re"),
        ("Re6+5", "Re6+5"),
        ("xRe6", "xRe6"),
        ("Re6x", "Re6x"),
        ("xRe3", "xRe3"),
        ("Re3x", "Re3x"),
        ("custom", "custom"),
    ];

    const REPRESENTATIVE_TOTAL_VALENCES: &[i32] = &[0, 3, 6, 7];
    const REPRESENTATIVE_FORMAL_CHARGES: &[i8] = &[-2, 0, 5];

    const PREDICATE_CASES: [(bool, bool, bool, bool); 8] = [
        (false, false, false, false),
        (false, false, true, true),
        (false, true, false, true),
        (false, true, true, true),
        (true, false, false, true),
        (true, false, true, true),
        (true, true, false, true),
        (true, true, true, true),
    ];

    fn fixed_atom(atomic_number: u8, formal_charge: i8, id: usize) -> Atom {
        let element = Element::from_atomic_number(atomic_number)
            .expect("fixed test atomic number is in the model range");
        Atom::from_spec(
            AtomId::new(id),
            AtomSpec::new(element).with_formal_charge(formal_charge),
        )
    }

    fn label_atom(
        atomic_number: u8,
        id: usize,
        hybridization: Hybridization,
        aromatic: bool,
        dummy_label: Option<&str>,
    ) -> Atom {
        let element = Element::from_atomic_number(atomic_number)
            .expect("fixed label test atomic number is in the model range");
        let spec = AtomSpec::new(element)
            .with_hybridization(hybridization)
            .with_aromatic(aromatic);
        let spec = if let Some(dummy_label) = dummy_label {
            spec.with_prop("dummyLabel", dummy_label)
                .expect("dummy label property has a nonempty key")
        } else {
            spec
        };
        Atom::from_spec(AtomId::new(id), spec)
    }

    fn typed_label_atom(atomic_number: u8, id: usize, dummy_label: PropertyValue) -> Atom {
        let element = Element::from_atomic_number(atomic_number)
            .expect("fixed label test atomic number is in the model range");
        let spec = AtomSpec::new(element)
            .with_hybridization(Hybridization::S)
            .with_prop("dummyLabel", dummy_label)
            .expect("dummy label property has a nonempty key");
        Atom::from_spec(AtomId::new(id), spec)
    }

    fn topology_from_atoms(atoms: Vec<Atom>) -> TopologyBlock {
        TopologyBlock::try_from_parts(atoms, Vec::new(), Vec::new(), Vec::new())
            .expect("fixed atom-only topology is structurally valid")
    }

    fn topology_with_bonds(
        atoms: Vec<Atom>,
        endpoint_orders: &[(usize, usize, BondOrder)],
    ) -> TopologyBlock {
        let bonds = endpoint_orders
            .iter()
            .enumerate()
            .map(|(bond_index, &(begin_atom, end_atom, order))| {
                let spec = BondSpec::new(AtomId::new(begin_atom), AtomId::new(end_atom), order)
                    .with_aromatic(order == BondOrder::Aromatic);
                Bond::from_spec(BondId::new(bond_index), spec)
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed endpoint bond topology is structurally valid")
    }

    fn cf3d_bld_integration_topology(
        atoms: Vec<Atom>,
        edges: &[(usize, usize, BondOrder, bool)],
    ) -> TopologyBlock {
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(bond_index, &(begin, end, order, conjugated))| {
                let mut bond = Bond::from_spec(
                    BondId::new(bond_index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                );
                bond.set_conjugated(conjugated);
                bond
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed integration topology is structurally valid")
    }

    fn cf3d_bld_integration_assignment(
        explicit_valence: &[i32],
        implicit_hydrogens: &[i32],
    ) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: explicit_valence.to_vec(),
            implicit_hydrogens: implicit_hydrogens.to_vec(),
        }
    }

    fn cf3d_bld_integration_prepare<'a>(
        topology: &'a TopologyBlock,
        assignment: &'a ValenceAssignment,
        table: &'a ParamCollection,
        diagnostics: &mut Vec<UffTypingDiagnostic>,
    ) -> (
        builder::PreparedTypingValence<'a>,
        Vec<bool>,
        Vec<Option<&'a AtomicParams>>,
        bool,
        bool,
    ) {
        let prepared = builder::prepare_typing_valence(topology, assignment)
            .expect("fixed source cache rows prepare without chemistry recomputation");
        let implicit_hydrogens = prepared.implicit_hydrogens().collect::<Vec<_>>();
        let needs_hydrogens =
            builder::needs_hydrogens_warning(topology, &implicit_hydrogens, diagnostics)
                .expect("fixed prepared hydrogen rows remain aligned");
        let conjugated = builder::prepare_conjugated_presence(topology)
            .expect("fixed supplied conjugation rows are structurally valid");
        let (slots, found_all) = get_atom_types(
            topology,
            &prepared.total_valences,
            &conjugated,
            table,
            diagnostics,
        )
        .expect("fixed prepared UFF typing inputs remain aligned");
        (prepared, conjugated, slots, found_all, needs_hydrogens)
    }

    fn cf3d_bld_integration_attach<'a>(field: &mut ForceField<'a>, rows: &'a mut [[f64; 3]]) {
        field
            .positions_mut()
            .extend(rows.iter_mut().map(|row| &mut row[..]));
    }

    fn cf3d_bld_integration_add_expected_bond(
        field: &mut ForceField<'_>,
        first: usize,
        second: usize,
        bond_order: f64,
        slots: &[Option<&AtomicParams>],
    ) {
        let contribution = BondStretchContrib::new(
            field.positions(),
            u32::try_from(first).expect("fixed source atom index fits"),
            u32::try_from(second).expect("fixed source atom index fits"),
            bond_order,
            slots[first].expect("fixed source endpoint parameter exists"),
            slots[second].expect("fixed source endpoint parameter exists"),
        )
        .expect("fixed source bond contribution constructs");
        field.add_contribution(Box::new(contribution));
    }

    fn cf3d_bld_integration_add_expected_angle(
        field: &mut ForceField<'_>,
        first: usize,
        center: usize,
        last: usize,
        first_order: f64,
        last_order: f64,
        order: u32,
        slots: &[Option<&AtomicParams>],
    ) {
        let contribution = AngleBendContrib::new(
            field.positions(),
            u32::try_from(first).expect("fixed source atom index fits"),
            u32::try_from(center).expect("fixed source atom index fits"),
            u32::try_from(last).expect("fixed source atom index fits"),
            first_order,
            last_order,
            slots[first].expect("fixed source first parameter exists"),
            slots[center].expect("fixed source center parameter exists"),
            slots[last].expect("fixed source last parameter exists"),
            order,
        )
        .expect("fixed source angle contribution constructs");
        field.add_contribution(Box::new(contribution));
    }

    fn cf3d_bld_integration_add_expected_torsion(
        field: &mut ForceField<'_>,
        topology: &TopologyBlock,
        indices: [usize; 4],
        bond_order: f64,
        terminal_sp2: bool,
        slots: &[Option<&AtomicParams>],
    ) {
        let [first, begin, end, last] = indices;
        let contribution = TorsionAngleContrib::new(
            field.positions(),
            u32::try_from(first).expect("fixed source atom index fits"),
            u32::try_from(begin).expect("fixed source atom index fits"),
            u32::try_from(end).expect("fixed source atom index fits"),
            u32::try_from(last).expect("fixed source atom index fits"),
            bond_order,
            i32::from(topology.atoms[begin].atomic_number()),
            i32::from(topology.atoms[end].atomic_number()),
            topology.atoms[begin].hybridization(),
            topology.atoms[end].hybridization(),
            slots[begin].expect("fixed source begin-center parameter exists"),
            slots[end].expect("fixed source end-center parameter exists"),
            terminal_sp2,
        )
        .expect("fixed source torsion contribution constructs");
        field.add_contribution(Box::new(contribution));
    }

    fn cf3d_bld_integration_add_expected_inversion(
        field: &mut ForceField<'_>,
        indices: [usize; 4],
        atomic_number: i32,
        bound_to_sp2_oxygen: bool,
    ) {
        let [first, center, second, third] = indices;
        let contribution = InversionContrib::new(
            field.positions(),
            u32::try_from(first).expect("fixed source atom index fits"),
            u32::try_from(center).expect("fixed source atom index fits"),
            u32::try_from(second).expect("fixed source atom index fits"),
            u32::try_from(third).expect("fixed source atom index fits"),
            atomic_number,
            bound_to_sp2_oxygen,
        )
        .expect("fixed source inversion contribution constructs");
        field.add_contribution(Box::new(contribution));
    }

    fn cf3d_bld_integration_assert_fields(
        mut actual: ForceField<'_>,
        mut expected: ForceField<'_>,
        coordinates: &[f64],
    ) {
        actual
            .initialize()
            .expect("actual fixed integration field initializes");
        expected
            .initialize()
            .expect("expected fixed source field initializes");
        let actual_energy = cf3d_bld_b05_calc_energy(&mut actual, coordinates)
            .expect("actual fixed integration energy evaluates");
        let expected_energy = cf3d_bld_b05_calc_energy(&mut expected, coordinates)
            .expect("expected fixed source energy evaluates");
        assert_eq!(actual_energy.to_bits(), expected_energy.to_bits());
        let mut actual_gradient = vec![0.0; coordinates.len()];
        let mut expected_gradient = vec![0.0; coordinates.len()];
        cf3d_bld_b05_calc_grad(&mut actual, coordinates, &mut actual_gradient)
            .expect("actual fixed integration gradient evaluates");
        cf3d_bld_b05_calc_grad(&mut expected, coordinates, &mut expected_gradient)
            .expect("expected fixed source gradient evaluates");
        assert_eq!(
            actual_gradient
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            expected_gradient
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
    }

    fn topology_with_bond(atoms: Vec<Atom>, order: BondOrder) -> TopologyBlock {
        topology_with_bonds(atoms, &[(0, 1, order)])
    }

    fn topology_with_angle(
        atoms: Vec<Atom>,
        order12: BondOrder,
        order23: BondOrder,
    ) -> TopologyBlock {
        topology_with_bonds(atoms, &[(0, 1, order12), (1, 2, order23)])
    }

    fn inversion_topology(
        central_atomic_number: u8,
        central_hybridization: Hybridization,
        central_degree: usize,
        terminal_sp2_oxygen_mask: u8,
    ) -> (TopologyBlock, [usize; 4]) {
        let atom_count = match central_degree {
            2 => 3,
            3 => 4,
            4 => 5,
            _ => panic!("fixed T17 fixture only uses degrees 2, 3, and 4"),
        };
        let indices = if central_degree == 2 {
            [0, 1, 2, 0]
        } else {
            [0, 1, 2, 3]
        };
        let atoms = (0..atom_count)
            .map(|atom_index| {
                if atom_index == 1 {
                    label_atom(
                        central_atomic_number,
                        atom_index,
                        central_hybridization,
                        false,
                        None,
                    )
                } else {
                    let oxygen_position = match atom_index {
                        0 => Some(0),
                        2 => Some(1),
                        3 => Some(2),
                        _ => None,
                    };
                    let is_sp2_oxygen = oxygen_position
                        .is_some_and(|position| terminal_sp2_oxygen_mask & (1 << position) != 0);
                    label_atom(
                        if is_sp2_oxygen { 8 } else { 6 },
                        atom_index,
                        if is_sp2_oxygen {
                            Hybridization::Sp2
                        } else {
                            Hybridization::Sp3
                        },
                        false,
                        None,
                    )
                }
            })
            .collect();
        let mut edges = vec![(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)];
        if central_degree >= 3 {
            edges.push((1, 3, BondOrder::Single));
        }
        if central_degree == 4 {
            edges.push((1, 4, BondOrder::Single));
        }
        (topology_with_bonds(atoms, &edges), indices)
    }

    fn expected_unpromoted_inversion_k(atomic_number: u8) -> Option<f64> {
        match atomic_number {
            6 | 7 | 8 => Some(2.0),
            15 => Some(4.496_661_509_237_746),
            33 => Some(4.086_825_183_419_91),
            51 => Some(3.979_001_005_763_598),
            83 => Some(3.666_666_666_666_667_4),
            _ => None,
        }
    }

    fn missing_atom_type_diagnostic(atom_id: AtomId) -> UffTypingDiagnostic {
        UffTypingDiagnostic {
            atom_id: Some(atom_id),
            kind: UffTypingDiagnosticKind::Error,
            message_prefix: UNRECOGNIZED_ATOM_TYPE_MESSAGE,
        }
    }

    #[test]
    fn cf3d_typ_t09_symbol_padding_and_dummy_symbol_precedence() {
        let cases = [
            (6, None, Hybridization::S, "C_"),
            (26, None, Hybridization::S, "Fe"),
            (0, None, Hybridization::Other, "*_"),
            (0, Some("Q"), Hybridization::Other, "Q_"),
            (0, Some("QQ"), Hybridization::Other, "QQ"),
            // std::string::size counts bytes, as does Rust str::len.
            (0, Some("λ"), Hybridization::Other, "λ"),
        ];

        for (offset, (atomic_number, dummy_label, hybridization, expected)) in
            cases.into_iter().enumerate()
        {
            let atom = label_atom(
                atomic_number,
                41_000 + offset,
                hybridization,
                false,
                dummy_label,
            );
            let mut diagnostics = Vec::new();
            assert_eq!(
                (atom_label_prefix(&atom, hybridization, || false, &mut diagnostics)
                    .expect("source-supported element label"))
                .as_bytes(),
                (expected).as_bytes(),
            );
            assert!(diagnostics.is_empty());
        }
    }

    #[test]
    fn uff_sync_symbol_absent_dummy_label_uses_source_fallback() {
        for (offset, (atomic_number, expected)) in [(0, "*_"), (6, "C_")].into_iter().enumerate() {
            let atom = label_atom(
                atomic_number,
                46_000 + offset,
                Hybridization::S,
                false,
                None,
            );
            let before = atom.clone();
            let mut diagnostics = Vec::new();

            assert_eq!(
                (atom_label_prefix(&atom, Hybridization::S, || false, &mut diagnostics)
                    .expect("source fallback returns the element symbol"))
                .as_bytes(),
                (expected).as_bytes(),
            );
            assert_eq!(atom, before);
            assert!(diagnostics.is_empty());
        }
    }

    #[test]
    fn uff_sync_symbol_projects_typed_dummy_labels_only_for_dummies() {
        let cases = [
            (PropertyValue::String("R".into()), "R_"),
            (PropertyValue::Int(0), "0_"),
            (PropertyValue::Int(-17), "-17"),
            (PropertyValue::Double(-0.0), "-0"),
            (PropertyValue::Double(0.1), "0.10000000000000001"),
            (PropertyValue::Bool(false), "0_"),
            (PropertyValue::Bool(true), "1_"),
        ];
        let mut actual_cases = 0;

        for (case_index, (property, expected_dummy_symbol)) in cases.into_iter().enumerate() {
            for atomic_number in [0, 6] {
                let atom_id = 46_100 + case_index * 2 + if atomic_number == 6 { 1 } else { 0 };
                let atom = typed_label_atom(atomic_number, atom_id, property.clone());
                let before = atom.clone();
                let mut diagnostics = Vec::new();
                let expected = if atomic_number == 0 {
                    expected_dummy_symbol
                } else {
                    "C_"
                };

                assert_eq!(
                    (atom_label_prefix(&atom, Hybridization::S, || false, &mut diagnostics)
                        .expect("all four modeled property kinds have source spellings"))
                    .as_bytes(),
                    (expected).as_bytes(),
                    "property={property:?}, atomic_number={atomic_number}",
                );
                assert_eq!(atom, before);
                assert_eq!(atom.prop("dummyLabel"), Some(&property));
                assert!(diagnostics.is_empty());
                actual_cases += 1;
            }
        }

        assert_eq!(actual_cases, 14);
    }

    #[test]
    fn uff_sync_symbol_conversion_error_preserves_stored_source_identity() {
        // This checks the typed source holder only. The current four modeled
        // PropertyValue variants all have source projections, so no actual
        // dummyLabel value can make property_value_to_string return this error.
        let cause = cosmolkit_core::PropertyStringError::UnsupportedKind {
            kind: cosmolkit_model::PropertyValueKind::Bool,
        };
        let error = UffTypingError::CorePropertyString(cause);
        let stored = match &error {
            UffTypingError::CorePropertyString(stored) => stored,
            _ => unreachable!(),
        };
        let exposed = std::error::Error::source(&error)
            .expect("the symbol conversion variant exposes its stored typed cause");
        let typed = exposed
            .downcast_ref::<cosmolkit_core::PropertyStringError>()
            .expect("the conversion cause keeps its concrete type");

        assert!(std::ptr::eq(stored, typed));
        assert_eq!(*typed, cause);
        assert!(std::error::Error::source(typed).is_none());
        assert_eq!(error.to_string(), format!("{error:?}"));
    }

    #[test]
    fn cf3d_typ_t09_forced_hybridization_suffixes_and_warnings() {
        const FORCED_SP3: &[(u8, &str)] = &[
            (12, "Mg"),
            (13, "Al"),
            (14, "Si"),
            (15, "P"),
            (50, "Sn"),
            (51, "Sb"),
            (52, "Te"),
            (81, "Tl"),
            (82, "Pb"),
            (83, "Bi"),
            (84, "Po"),
        ];

        let mut next_id = 42_000;
        for &(atomic_number, symbol) in FORCED_SP3 {
            for (hybridization, expected_warning) in
                [(Hybridization::Sp3, false), (Hybridization::Other, true)]
            {
                let atom = label_atom(atomic_number, next_id, hybridization, false, None);
                let mut diagnostics = Vec::new();
                let expected = format!("{}{}3", symbol, if symbol.len() == 1 { "_" } else { "" });
                assert_eq!(
                    (atom_label_prefix(&atom, hybridization, || false, &mut diagnostics)
                        .expect("source-supported forced label"))
                    .as_bytes(),
                    (expected).as_bytes(),
                    "atomic number {atomic_number}, hybridization {hybridization:?}",
                );
                let expected_diagnostics = if expected_warning {
                    vec![UffTypingDiagnostic {
                        atom_id: Some(atom.id()),
                        kind: UffTypingDiagnosticKind::Warning,
                        message_prefix: FORCED_SP3_HYBRIDIZATION_WARNING_MESSAGE,
                    }]
                } else {
                    Vec::new()
                };
                assert_eq!(diagnostics, expected_diagnostics);
                next_id += 1;
            }
        }

        for (hybridization, expected_warning) in
            [(Hybridization::Sp, false), (Hybridization::Other, true)]
        {
            let atom = label_atom(80, next_id, hybridization, false, None);
            let mut diagnostics = Vec::new();
            assert_eq!(
                (atom_label_prefix(&atom, hybridization, || false, &mut diagnostics)
                    .expect("source-supported Hg label"))
                .as_bytes(),
                ("Hg1").as_bytes(),
            );
            let expected_diagnostics = if expected_warning {
                vec![UffTypingDiagnostic {
                    atom_id: Some(atom.id()),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
                }]
            } else {
                Vec::new()
            };
            assert_eq!(diagnostics, expected_diagnostics);
            next_id += 1;
        }
    }

    #[test]
    fn cf3d_typ_t09_hybridization_aromaticity_and_prepared_conjugation_matrix() {
        const HYBRIDIZATIONS: [Hybridization; 9] = [
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp2,
            Hybridization::Sp3,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Unspecified,
            Hybridization::Other,
        ];
        const ELEMENTS: &[(u8, &str)] = &[(6, "C"), (7, "N"), (8, "O"), (16, "S"), (5, "B")];

        let mut id = 43_000;
        for &(atomic_number, symbol) in ELEMENTS {
            for hybridization in HYBRIDIZATIONS {
                for aromatic in [false, true] {
                    for atom_has_conjugated_bond in [false, true] {
                        let atom = label_atom(atomic_number, id, hybridization, aromatic, None);
                        let (suffix, expected_diagnostic) = match hybridization {
                            Hybridization::S => ("", None),
                            Hybridization::Sp => ("1", None),
                            Hybridization::Sp2 => {
                                if (aromatic || atom_has_conjugated_bond)
                                    && matches!(atomic_number, 6 | 7 | 8 | 16)
                                {
                                    ("R", None)
                                } else {
                                    ("2", None)
                                }
                            }
                            Hybridization::Sp3 => ("3", None),
                            Hybridization::Sp2d => ("4", None),
                            Hybridization::Sp3d => ("5", None),
                            Hybridization::Sp3d2 => ("6", None),
                            Hybridization::Unspecified | Hybridization::Other => {
                                ("", Some(UffTypingDiagnosticKind::Error))
                            }
                        };
                        let expected = format!(
                            "{}{}{}",
                            symbol,
                            if symbol.len() == 1 { "_" } else { "" },
                            suffix,
                        );
                        let mut diagnostics = Vec::new();
                        assert_eq!(
                            (atom_label_prefix(
                                &atom,
                                hybridization,
                                || atom_has_conjugated_bond,
                                &mut diagnostics,
                            )
                            .expect("source-supported element label"))
                            .as_bytes(),
                            (expected).as_bytes(),
                            "atomic number {atomic_number}, hybridization {hybridization:?}, aromatic {aromatic}, conjugated {atom_has_conjugated_bond}",
                        );
                        let expected_diagnostics = expected_diagnostic
                            .map(|kind| {
                                vec![UffTypingDiagnostic {
                                    atom_id: Some(atom.id()),
                                    kind,
                                    message_prefix: UNRECOGNIZED_HYBRIDIZATION_MESSAGE,
                                }]
                            })
                            .unwrap_or_default();
                        assert_eq!(diagnostics, expected_diagnostics);
                        id += 1;
                    }
                }
            }
        }
        assert_eq!(id - 43_000, 180);
    }

    #[test]
    fn uff_prepare_p03_source_routes_keep_literal_labels_and_lazy_reads() {
        const HYBRIDIZATIONS: [Hybridization; 9] = [
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp2,
            Hybridization::Sp3,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Unspecified,
            Hybridization::Other,
        ];
        const DEFAULT_SWITCH_ELEMENTS: &[(u8, &str)] =
            &[(6, "C_"), (7, "N_"), (8, "O_"), (16, "S_"), (5, "B_")];
        const SPECIAL_ELEMENTS: &[(u8, &str, bool)] = &[
            (12, "Mg3", true),
            (13, "Al3", true),
            (14, "Si3", true),
            (15, "P_3", true),
            (50, "Sn3", true),
            (51, "Sb3", true),
            (52, "Te3", true),
            (80, "Hg1", false),
            (81, "Tl3", true),
            (82, "Pb3", true),
            (83, "Bi3", true),
            (84, "Po3", true),
        ];

        let mut id = 47_000;
        let mut default_switch_calls = 0;
        for &(atomic_number, symbol) in DEFAULT_SWITCH_ELEMENTS {
            for hybridization in HYBRIDIZATIONS {
                for aromatic in [false, true] {
                    let atom = label_atom(atomic_number, id, hybridization, aromatic, None);
                    let mut calls = 0;
                    let mut diagnostics = Vec::new();
                    let actual = atom_label_prefix(
                        &atom,
                        hybridization,
                        || {
                            calls += 1;
                            default_switch_calls += 1;
                            false
                        },
                        &mut diagnostics,
                    )
                    .expect("fixed default-switch element has a source label");

                    let suffix = match hybridization {
                        Hybridization::S => "",
                        Hybridization::Sp => "1",
                        Hybridization::Sp2 => {
                            if aromatic && matches!(atomic_number, 6 | 7 | 8 | 16) {
                                "R"
                            } else {
                                "2"
                            }
                        }
                        Hybridization::Sp3 => "3",
                        Hybridization::Sp2d => "4",
                        Hybridization::Sp3d => "5",
                        Hybridization::Sp3d2 => "6",
                        Hybridization::Unspecified | Hybridization::Other => "",
                    };
                    assert_eq!(
                        (actual).as_bytes(),
                        (format!("{symbol}{suffix}")).as_bytes()
                    );
                    let expected_diagnostics = if matches!(
                        hybridization,
                        Hybridization::Unspecified | Hybridization::Other
                    ) {
                        vec![UffTypingDiagnostic {
                            atom_id: Some(atom.id()),
                            kind: UffTypingDiagnosticKind::Error,
                            message_prefix: UNRECOGNIZED_HYBRIDIZATION_MESSAGE,
                        }]
                    } else {
                        Vec::new()
                    };
                    assert_eq!(diagnostics, expected_diagnostics);
                    let expected_calls = if hybridization == Hybridization::Sp2 && !aromatic {
                        1
                    } else {
                        0
                    };
                    assert_eq!(calls, expected_calls);
                    id += 1;
                }
            }
        }
        assert_eq!(default_switch_calls, DEFAULT_SWITCH_ELEMENTS.len());

        let mut special_switch_calls = 0;
        for &(atomic_number, expected_label, forced_sp3) in SPECIAL_ELEMENTS {
            for hybridization in HYBRIDIZATIONS {
                for aromatic in [false, true] {
                    let atom = label_atom(atomic_number, id, hybridization, aromatic, None);
                    let mut diagnostics = Vec::new();
                    let actual = atom_label_prefix(
                        &atom,
                        hybridization,
                        || {
                            special_switch_calls += 1;
                            false
                        },
                        &mut diagnostics,
                    )
                    .expect("fixed source special-element branch has a label");
                    assert_eq!((actual).as_bytes(), (expected_label).as_bytes());
                    let expected_warning = if forced_sp3 {
                        (hybridization != Hybridization::Sp3)
                            .then_some(FORCED_SP3_HYBRIDIZATION_WARNING_MESSAGE)
                    } else {
                        (hybridization != Hybridization::Sp)
                            .then_some(FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE)
                    };
                    let expected_diagnostics = expected_warning
                        .map(|message_prefix| {
                            vec![UffTypingDiagnostic {
                                atom_id: Some(atom.id()),
                                kind: UffTypingDiagnosticKind::Warning,
                                message_prefix,
                            }]
                        })
                        .unwrap_or_default();
                    assert_eq!(diagnostics, expected_diagnostics);
                    id += 1;
                }
            }
        }
        assert_eq!(special_switch_calls, 0);

        for (aromatic, expected, expected_calls) in [(false, "C_2", 1), (true, "C_R", 0)] {
            let atom = label_atom(6, id, Hybridization::Sp2, aromatic, None);
            let mut calls = 0;
            let mut diagnostics = Vec::new();
            assert_eq!(
                (get_atom_label(
                    &atom,
                    2,
                    Hybridization::Sp2,
                    || {
                        calls += 1;
                        false
                    },
                    &mut diagnostics
                )
                .expect("source getAtomLabel forwards its lazy conjugation input"))
                .as_bytes(),
                (expected).as_bytes(),
            );
            assert_eq!(calls, expected_calls);
            assert!(diagnostics.is_empty());
            id += 1;
        }

        assert_eq!(id - 47_000, 308);
    }

    #[test]
    fn cf3d_typ_t09_alkali_and_halogen_gate_classes() {
        const SUPPRESSED_ELEMENTS: &[(u8, &str)] = &[
            (3, "Li"),
            (11, "Na"),
            (19, "K"),
            (37, "Rb"),
            (55, "Cs"),
            (87, "Fr"),
            (9, "F"),
            (17, "Cl"),
            (35, "Br"),
            (53, "I"),
            (85, "At"),
        ];
        const HYBRIDIZATIONS: [Hybridization; 9] = [
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp2,
            Hybridization::Sp3,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Unspecified,
            Hybridization::Other,
        ];

        let mut id = 44_000;
        for &(atomic_number, symbol) in SUPPRESSED_ELEMENTS {
            for hybridization in HYBRIDIZATIONS {
                let atom = label_atom(atomic_number, id, hybridization, false, None);
                let expected = if symbol.len() == 1 {
                    format!("{symbol}_")
                } else {
                    symbol.to_owned()
                };
                let mut diagnostics = Vec::new();
                assert_eq!(
                    (atom_label_prefix(&atom, hybridization, || false, &mut diagnostics)
                        .expect("source-supported alkali/halogen label"))
                    .as_bytes(),
                    (expected).as_bytes(),
                    "atomic number {atomic_number}, hybridization {hybridization:?}",
                );
                assert!(diagnostics.is_empty());
                id += 1;
            }
        }

        // Source tests default valence -1 before outer-electron exceptions, so
        // the sentinel gate still enters the hybridization switch for both a
        // non-halogen transition element and group-7 Ts.
        for (atomic_number, symbol) in [(21, "Sc"), (117, "Ts")] {
            let atom = label_atom(atomic_number, id, Hybridization::Sp, false, None);
            let mut diagnostics = Vec::new();
            assert_eq!(
                (atom_label_prefix(&atom, Hybridization::Sp, || false, &mut diagnostics)
                    .expect("source-supported default-valence sentinel label"))
                .as_bytes(),
                (format!("{symbol}1")).as_bytes(),
            );
            assert!(diagnostics.is_empty());
            id += 1;
        }

        let atom = label_atom(117, id, Hybridization::Unspecified, false, None);
        let mut diagnostics = Vec::new();
        assert_eq!(
            (atom_label_prefix(
                &atom,
                Hybridization::Unspecified,
                || false,
                &mut diagnostics
            )
            .expect("source-supported default-valence sentinel label"))
            .as_bytes(),
            ("Ts").as_bytes(),
        );
        assert_eq!(
            diagnostics,
            &[UffTypingDiagnostic {
                atom_id: Some(atom.id()),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: UNRECOGNIZED_HYBRIDIZATION_MESSAGE,
            }]
        );
    }

    #[test]
    fn cf3d_typ_t10_hybridization_family_labels_and_errors() {
        const CASES: &[(Hybridization, bool, bool, &str, Option<&str>)] = &[
            (Hybridization::S, false, false, "C_", None),
            (Hybridization::Sp, false, false, "C_1", None),
            (Hybridization::Sp2, true, false, "C_R", None),
            (Hybridization::Sp3, false, false, "C_3", None),
            (Hybridization::Sp2d, false, false, "C_4", None),
            (Hybridization::Sp3d, false, false, "C_5", None),
            (Hybridization::Sp3d2, false, false, "C_6", None),
            (
                Hybridization::Unspecified,
                false,
                false,
                "C_",
                Some(UNRECOGNIZED_HYBRIDIZATION_MESSAGE),
            ),
            (
                Hybridization::Other,
                false,
                false,
                "C_",
                Some(UNRECOGNIZED_HYBRIDIZATION_MESSAGE),
            ),
        ];

        for (offset, &(hybridization, aromatic, conjugated, expected, error)) in
            CASES.iter().enumerate()
        {
            let atom = label_atom(6, 45_000 + offset, hybridization, aromatic, None);
            let mut diagnostics = Vec::new();
            assert_eq!(
                (get_atom_label(&atom, 0, hybridization, || conjugated, &mut diagnostics)
                    .expect("carbon has a source-defined label"))
                .as_bytes(),
                (expected).as_bytes(),
                "hybridization {hybridization:?}",
            );
            let expected_diagnostics = error
                .map(|message_prefix| {
                    vec![UffTypingDiagnostic {
                        atom_id: Some(atom.id()),
                        kind: UffTypingDiagnosticKind::Error,
                        message_prefix,
                    }]
                })
                .unwrap_or_default();
            assert_eq!(diagnostics, expected_diagnostics);
        }
    }

    #[test]
    fn cf3d_typ_t10_aromatic_and_prepared_conjugation_label_matrix() {
        const CASES: &[(u8, bool, bool, &str)] = &[
            (6, false, false, "C_2"),
            (6, false, true, "C_R"),
            (6, true, false, "C_R"),
            (6, true, true, "C_R"),
            (7, false, false, "N_2"),
            (7, false, true, "N_R"),
            (7, true, false, "N_R"),
            (7, true, true, "N_R"),
            (8, false, false, "O_2"),
            (8, false, true, "O_R"),
            (8, true, false, "O_R"),
            (8, true, true, "O_R"),
            (16, false, false, "S_2"),
            (16, false, true, "S_R"),
            (16, true, false, "S_R"),
            (16, true, true, "S_R"),
        ];

        for (offset, &(atomic_number, aromatic, conjugated, expected)) in CASES.iter().enumerate() {
            let atom = label_atom(
                atomic_number,
                45_100 + offset,
                Hybridization::Sp2,
                aromatic,
                None,
            );
            let mut diagnostics = Vec::new();
            assert_eq!(
                (get_atom_label(
                    &atom,
                    2,
                    Hybridization::Sp2,
                    || conjugated,
                    &mut diagnostics
                )
                .expect("source-supported aromatic/conjugated label"))
                .as_bytes(),
                (expected).as_bytes(),
                "atomic number {atomic_number}, aromatic {aromatic}, conjugated {conjugated}",
            );
            assert!(diagnostics.is_empty());
        }
    }

    #[test]
    fn cf3d_typ_t10_charge_families_and_special_elements() {
        const CASES: &[(u8, Hybridization, i32, &str, bool)] = &[
            (29, Hybridization::Sp3, 1, "Cu3+1", false),
            (4, Hybridization::Sp3, 2, "Be3+2", false),
            (21, Hybridization::Sp3, 3, "Sc3+3", false),
            (18, Hybridization::Sp3, 4, "Ar3+4", false),
            (23, Hybridization::Sp3, 5, "V_3+5", false),
            (42, Hybridization::Sp3, 6, "Mo3+6", false),
            (12, Hybridization::Sp3, 2, "Mg3+2", false),
            (30, Hybridization::Sp3, 2, "Zn3+2", false),
            (31, Hybridization::Sp3, 3, "Ga3+3", false),
            (34, Hybridization::Sp3, 2, "Se3+2", false),
            (48, Hybridization::Sp3, 2, "Cd3+2", false),
            (49, Hybridization::Sp3, 3, "In3+3", false),
            (51, Hybridization::Sp3, 3, "Sb3+3", false),
            (52, Hybridization::Sp3, 2, "Te3+2", false),
            (80, Hybridization::Sp, 2, "Hg1+2", false),
            (81, Hybridization::Sp3, 3, "Tl3+3", false),
            (82, Hybridization::Sp3, 3, "Pb3+3", false),
            (83, Hybridization::Sp3, 3, "Bi3+3", false),
            (84, Hybridization::Sp3, 2, "Po3+2", false),
            (13, Hybridization::Sp3, 3, "Al3", false),
            (14, Hybridization::Sp3, 4, "Si3", false),
            (15, Hybridization::Sp3, 3, "P_3+3", false),
            (15, Hybridization::Sp3, 5, "P_3+5", false),
            (15, Hybridization::Sp3, 0, "P_3+5", true),
            (16, Hybridization::Sp3, 2, "S_3+2", false),
            (57, Hybridization::Sp3, 6, "La3+3", false),
            (64, Hybridization::Sp3, 6, "Gd3+3", false),
            (71, Hybridization::Sp3, 6, "Lu3+3", false),
            (57, Hybridization::Sp3, 0, "La3+3", true),
            (75, Hybridization::Sp3, 0, "Re3+7", true),
        ];

        for (offset, &(atomic_number, hybridization, total_valence, expected, charge_error)) in
            CASES.iter().enumerate()
        {
            let atom = label_atom(atomic_number, 45_200 + offset, hybridization, false, None);
            let mut diagnostics = Vec::new();
            assert_eq!(
                (get_atom_label(
                    &atom,
                    total_valence,
                    hybridization,
                    || false,
                    &mut diagnostics
                )
                .expect("source-supported charge label"))
                .as_bytes(),
                (expected).as_bytes(),
                "atomic number {atomic_number}, valence {total_valence}",
            );
            let expected_diagnostics = if charge_error {
                vec![UffTypingDiagnostic {
                    atom_id: Some(atom.id()),
                    kind: UffTypingDiagnosticKind::Error,
                    message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
                }]
            } else {
                Vec::new()
            };
            assert_eq!(diagnostics, expected_diagnostics);
        }
    }

    #[test]
    fn cf3d_typ_t10_dummy_and_diagnostic_order() {
        let dummy = label_atom(0, 45_300, Hybridization::Unspecified, false, Some("D"));
        let mut diagnostics = Vec::new();
        assert_eq!(
            (get_atom_label(
                &dummy,
                0,
                Hybridization::Unspecified,
                || false,
                &mut diagnostics,
            )
            .expect("dummyLabel supplies the source symbol"))
            .as_bytes(),
            ("D_").as_bytes(),
        );
        assert!(diagnostics.is_empty());

        let magnesium = label_atom(12, 45_301, Hybridization::Other, false, None);
        diagnostics.clear();
        assert_eq!(
            (get_atom_label(
                &magnesium,
                0,
                Hybridization::Other,
                || false,
                &mut diagnostics
            )
            .expect("source tolerance keeps the charge suffix"))
            .as_bytes(),
            ("Mg3+2").as_bytes(),
        );
        assert_eq!(
            diagnostics,
            [
                UffTypingDiagnostic {
                    atom_id: Some(magnesium.id()),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP3_HYBRIDIZATION_WARNING_MESSAGE,
                },
                UffTypingDiagnostic {
                    atom_id: Some(magnesium.id()),
                    kind: UffTypingDiagnosticKind::Error,
                    message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
                },
            ],
        );

        let sulfur = label_atom(16, 45_302, Hybridization::Unspecified, false, None);
        diagnostics.clear();
        assert_eq!(
            (get_atom_label(
                &sulfur,
                0,
                Hybridization::Unspecified,
                || false,
                &mut diagnostics,
            )
            .expect("non-SP2 sulfur keeps the tolerated source suffix"))
            .as_bytes(),
            ("S_+6").as_bytes(),
        );
        assert_eq!(
            diagnostics,
            [
                UffTypingDiagnostic {
                    atom_id: Some(sulfur.id()),
                    kind: UffTypingDiagnosticKind::Error,
                    message_prefix: UNRECOGNIZED_HYBRIDIZATION_MESSAGE,
                },
                UffTypingDiagnostic {
                    atom_id: Some(sulfur.id()),
                    kind: UffTypingDiagnosticKind::Error,
                    message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
                },
            ],
        );
    }

    #[test]
    fn cf3d_typ_t11_empty_and_all_known_rows() {
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();
        let empty = TopologyBlock::default();
        let (empty_slots, empty_found_all) =
            get_atom_types(&empty, &[], &[], params.as_ref(), &mut diagnostics)
                .expect("empty topology has no missing row");
        assert!(empty_slots.is_empty());
        assert!(empty_found_all);
        assert!(diagnostics.is_empty());

        let topology = topology_from_atoms(vec![
            label_atom(6, 0, Hybridization::Sp3, false, None),
            label_atom(6, 1, Hybridization::Sp3, false, None),
        ]);
        let (slots, found_all) = get_atom_types(
            &topology,
            &[4, 4],
            &[false, false],
            params.as_ref(),
            &mut diagnostics,
        )
        .expect("both carbon labels are present in the default table");
        let carbon_sp3 = params.get("C_3").expect("fixed default C_3 row");
        assert!(found_all);
        assert_eq!(slots.len(), 2);
        assert!(std::ptr::eq(slots[0].expect("first C_3 row"), carbon_sp3));
        assert!(std::ptr::eq(slots[1].expect("second C_3 row"), carbon_sp3));
        assert!(std::ptr::eq(
            slots[0].expect("first C_3 row"),
            slots[1].expect("second C_3 row"),
        ));
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t11_missing_first_middle_last_continue_and_preserve_reference() {
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let topology = topology_from_atoms(vec![
            label_atom(
                0,
                0,
                Hybridization::Unspecified,
                false,
                Some("unknown-first"),
            ),
            label_atom(6, 1, Hybridization::Sp3, false, None),
            label_atom(
                0,
                2,
                Hybridization::Unspecified,
                false,
                Some("unknown-middle"),
            ),
            label_atom(6, 3, Hybridization::Sp3, false, None),
            label_atom(
                0,
                4,
                Hybridization::Unspecified,
                false,
                Some("unknown-last"),
            ),
        ]);
        let mut diagnostics = Vec::new();
        let (slots, found_all) = get_atom_types(
            &topology,
            &[0, 4, 0, 4, 0],
            &[false, false, false, false, false],
            params.as_ref(),
            &mut diagnostics,
        )
        .expect("an absent type produces a null row and does not abort");

        assert!(!found_all);
        assert_eq!(slots.len(), 5);
        assert_eq!(slots[0], None);
        assert_eq!(slots[2], None);
        assert_eq!(slots[4], None);
        let carbon_sp3 = params.get("C_3").expect("fixed default C_3 row");
        assert!(std::ptr::eq(
            slots[1].expect("row after first miss"),
            carbon_sp3
        ));
        assert!(std::ptr::eq(
            slots[3].expect("row after middle miss"),
            carbon_sp3
        ));
        assert!(std::ptr::eq(
            slots[1].expect("first C_3 row"),
            slots[3].expect("second C_3 row"),
        ));
        assert_eq!(
            diagnostics,
            [
                missing_atom_type_diagnostic(AtomId::new(0)),
                missing_atom_type_diagnostic(AtomId::new(2)),
                missing_atom_type_diagnostic(AtomId::new(4)),
            ],
        );
    }

    #[test]
    fn cf3d_typ_t11_prepared_slice_lengths_fail_before_rows() {
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let topology = topology_from_atoms(vec![label_atom(6, 0, Hybridization::Sp3, false, None)]);
        let mut diagnostics = Vec::new();

        assert_eq!(
            get_atom_types(&topology, &[], &[false], params.as_ref(), &mut diagnostics),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: 1,
                actual: 0,
            }),
        );
        assert_eq!(
            get_atom_types(&topology, &[4], &[], params.as_ref(), &mut diagnostics),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::ConjugatedBondPresence,
                expected: 1,
                actual: 0,
            }),
        );
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn uff_worker_w09_returns_source_found_all_and_keeps_diagnostic_order() {
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let empty = TopologyBlock::default();
        let mut diagnostics = Vec::new();
        assert_eq!(
            uff_has_all_molecule_parameters(&empty, &[], &[], params.as_ref(), &mut diagnostics),
            Ok(true)
        );
        assert!(diagnostics.is_empty());

        let known_topology = topology_from_atoms(vec![
            label_atom(6, 0, Hybridization::Sp3, false, None),
            label_atom(6, 1, Hybridization::Sp3, false, None),
        ]);
        assert_eq!(
            uff_has_all_molecule_parameters(
                &known_topology,
                &[4, 4],
                &[false, false],
                params.as_ref(),
                &mut diagnostics,
            ),
            Ok(true)
        );
        assert!(diagnostics.is_empty());

        let missing_topology = topology_from_atoms(vec![
            label_atom(
                0,
                0,
                Hybridization::Unspecified,
                false,
                Some("unknown-first"),
            ),
            label_atom(6, 1, Hybridization::Sp3, false, None),
            label_atom(
                0,
                2,
                Hybridization::Unspecified,
                false,
                Some("unknown-middle"),
            ),
            label_atom(6, 3, Hybridization::Sp3, false, None),
            label_atom(
                0,
                4,
                Hybridization::Unspecified,
                false,
                Some("unknown-last"),
            ),
        ]);
        assert_eq!(
            uff_has_all_molecule_parameters(
                &missing_topology,
                &[0, 4, 0, 4, 0],
                &[false; 5],
                params.as_ref(),
                &mut diagnostics,
            ),
            Ok(false)
        );
        assert_eq!(
            diagnostics,
            [
                missing_atom_type_diagnostic(AtomId::new(0)),
                missing_atom_type_diagnostic(AtomId::new(2)),
                missing_atom_type_diagnostic(AtomId::new(4)),
            ]
        );

        let one_atom = topology_from_atoms(vec![label_atom(6, 0, Hybridization::Sp3, false, None)]);
        let mut error_diagnostics = Vec::new();
        assert_eq!(
            uff_has_all_molecule_parameters(
                &one_atom,
                &[],
                &[false],
                params.as_ref(),
                &mut error_diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: 1,
                actual: 0,
            })
        );
        assert_eq!(
            uff_has_all_molecule_parameters(
                &one_atom,
                &[4],
                &[],
                params.as_ref(),
                &mut error_diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::ConjugatedBondPresence,
                expected: 1,
                actual: 0,
            })
        );
        assert!(error_diagnostics.is_empty());
    }

    fn source_charge_dispatch_expectation(
        atomic_number: u8,
        total_valence: i32,
        formal_charge: i8,
        tolerate_charge_mismatch: bool,
    ) -> (&'static str, bool) {
        if let Some((_, required_valence, suffix)) = FIXED_CHARGE_GROUPS
            .iter()
            .find(|(elements, _, _)| elements.contains(&atomic_number))
        {
            let valid = total_valence == *required_valence
                || i32::from(formal_charge) == *required_valence
                || tolerate_charge_mismatch;
            return (if valid { suffix } else { "" }, !valid);
        }

        if let Some((_, required_valence, suffix)) = VALENCE_ONLY_CHARGE_GROUPS
            .iter()
            .find(|(elements, _, _)| elements.contains(&atomic_number))
        {
            let valid = total_valence == *required_valence;
            return (
                if valid || tolerate_charge_mismatch {
                    suffix
                } else {
                    ""
                },
                !valid,
            );
        }

        if let Some((_, required_valence)) = UNSUFFIXED_MAIN_GROUP_CHARGES
            .iter()
            .find(|(element, _)| *element == atomic_number)
        {
            return ("", total_valence != *required_valence);
        }

        match atomic_number {
            15 => match total_valence {
                3 => ("+3", false),
                5 => ("+5", false),
                _ => (if tolerate_charge_mismatch { "+5" } else { "" }, true),
            },
            16 => match total_valence {
                2 => ("+2", false),
                4 => ("+4", false),
                6 => ("+6", false),
                _ => (if tolerate_charge_mismatch { "+6" } else { "" }, true),
            },
            75 => ("", true),
            57..=71 => {
                if total_valence == 6 {
                    ("+3", false)
                } else {
                    (if tolerate_charge_mismatch { "+3" } else { "" }, true)
                }
            }
            _ => ("", false),
        }
    }

    fn assert_mismatch_diagnostic(diagnostics: &[UffTypingDiagnostic], atom_id: AtomId) {
        assert_eq!(
            diagnostics,
            &[UffTypingDiagnostic {
                atom_id: Some(atom_id),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
            }]
        );
    }

    #[test]
    fn cf3d_typ_t08_all_elements_charge_dispatch_matrix() {
        const FORMAL_CHARGES: [i8; 9] = [-2, 0, 1, 2, 3, 4, 5, 6, 7];

        let mut atom_index = 20_000;
        for element in Element::iter_with_dummy() {
            let atomic_number = element.atomic_number();
            for total_valence in 0..=7 {
                for formal_charge in FORMAL_CHARGES {
                    for tolerate_charge_mismatch in [false, true] {
                        let atom = Atom::from_spec(
                            AtomId::new(atom_index),
                            AtomSpec::new(element)
                                .with_formal_charge(formal_charge)
                                .with_hybridization(Hybridization::Sp3),
                        );
                        let atom_id = atom.id();
                        let (expected_suffix, emits_error) = source_charge_dispatch_expectation(
                            atomic_number,
                            total_valence,
                            formal_charge,
                            tolerate_charge_mismatch,
                        );
                        let mut atom_key = cosmolkit_model::PropertyText::from("matrix-prefix");
                        let mut diagnostics = Vec::new();

                        add_atom_charge_flags(
                            &atom,
                            total_valence,
                            &mut atom_key,
                            tolerate_charge_mismatch,
                            &mut diagnostics,
                        );

                        assert_eq!(
                            (atom_key).as_bytes(),
                            (format!("matrix-prefix{expected_suffix}")).as_bytes(),
                            "atomic number {atomic_number}, valence {total_valence}, formal charge {formal_charge}, tolerance {tolerate_charge_mismatch}",
                        );
                        let expected_diagnostics = if emits_error {
                            vec![UffTypingDiagnostic {
                                atom_id: Some(atom_id),
                                kind: UffTypingDiagnosticKind::Error,
                                message_prefix: UNRECOGNIZED_CHARGE_STATE_MESSAGE,
                            }]
                        } else {
                            Vec::new()
                        };
                        assert_eq!(
                            diagnostics, expected_diagnostics,
                            "atomic number {atomic_number}, valence {total_valence}, formal charge {formal_charge}, tolerance {tolerate_charge_mismatch}",
                        );
                        atom_index += 1;
                    }
                }
            }
        }
        assert_eq!(atom_index, 37_136);
    }

    #[test]
    fn cf3d_typ_t01_all_elements_and_eight_predicate_cases() {
        let mut atom_index = 0;
        for &(elements, required_valence, suffix) in FIXED_CHARGE_GROUPS {
            for &atomic_number in elements {
                for &(valence_matches, formal_charge_matches, tolerate, appends) in &PREDICATE_CASES
                {
                    let formal_charge = if formal_charge_matches {
                        required_valence as i8
                    } else {
                        -1
                    };
                    let atom = fixed_atom(atomic_number, formal_charge, atom_index);
                    let atom_id = atom.id();
                    let total_valence = if valence_matches {
                        required_valence
                    } else {
                        required_valence + 1
                    };
                    let mut atom_key = cosmolkit_model::PropertyText::from("fixed-prefix");
                    let mut diagnostics = Vec::new();

                    assert!(append_fixed_charge_flag(
                        &atom,
                        total_valence,
                        &mut atom_key,
                        tolerate,
                        &mut diagnostics,
                    ));
                    let expected_key = if appends {
                        format!("fixed-prefix{suffix}")
                    } else {
                        String::from("fixed-prefix")
                    };
                    assert_eq!(
                        atom_key.as_bytes(),
                        expected_key.as_bytes(),
                        "atomic number {atomic_number}"
                    );
                    if appends {
                        assert!(diagnostics.is_empty(), "atomic number {atomic_number}");
                    } else {
                        assert_mismatch_diagnostic(&diagnostics, atom_id);
                    }
                    atom_index += 1;
                }
            }
        }
    }

    #[test]
    fn cf3d_typ_t01_positive_nonmatching_formal_charge_logs_without_key_change() {
        let mut atom_index = 10_000;
        for &(elements, required_valence, _) in FIXED_CHARGE_GROUPS {
            for &atomic_number in elements {
                let other_positive_charge = if required_valence == 1 { 2 } else { 1 };
                let atom = fixed_atom(atomic_number, other_positive_charge, atom_index);
                let atom_id = atom.id();
                let mut atom_key = cosmolkit_model::PropertyText::from("kept-prefix");
                let mut diagnostics = Vec::new();

                assert!(append_fixed_charge_flag(
                    &atom,
                    required_valence + 1,
                    &mut atom_key,
                    false,
                    &mut diagnostics,
                ));
                assert_eq!(
                    (atom_key).as_bytes(),
                    ("kept-prefix").as_bytes(),
                    "atomic number {atomic_number}"
                );
                assert_mismatch_diagnostic(&diagnostics, atom_id);
                atom_index += 1;
            }
        }
    }

    #[test]
    fn cf3d_typ_t01_nonmember_is_not_handled_or_mutated() {
        let atom = fixed_atom(1, 6, 20_000);
        for tolerate in [false, true] {
            let mut atom_key = cosmolkit_model::PropertyText::from("kept-prefix");
            let mut diagnostics = Vec::new();
            assert!(!append_fixed_charge_flag(
                &atom,
                0,
                &mut atom_key,
                tolerate,
                &mut diagnostics,
            ));
            assert_eq!((atom_key).as_bytes(), ("kept-prefix").as_bytes());
            assert!(diagnostics.is_empty());
        }
    }

    #[test]
    fn cf3d_typ_t02_all_elements_valence_tolerance_formal_charge_matrix() {
        let mut atom_index = 30_000;
        let mut case_count = 0;

        for &(elements, required_valence, suffix) in VALENCE_ONLY_CHARGE_GROUPS {
            for &atomic_number in elements {
                for valence_matches in [false, true] {
                    for tolerate_charge_mismatch in [false, true] {
                        for formal_charge_matches in [false, true] {
                            let formal_charge = if formal_charge_matches {
                                required_valence as i8
                            } else {
                                -1
                            };
                            let atom = fixed_atom(atomic_number, formal_charge, atom_index);
                            let atom_id = atom.id();
                            let total_valence = if valence_matches {
                                required_valence
                            } else {
                                required_valence + 1
                            };
                            let mut atom_key = cosmolkit_model::PropertyText::from("t02-prefix");
                            let mut diagnostics = Vec::new();

                            assert!(append_valence_only_charge_flag(
                                &atom,
                                total_valence,
                                &mut atom_key,
                                tolerate_charge_mismatch,
                                &mut diagnostics,
                            ));
                            let expected_key = if valence_matches || tolerate_charge_mismatch {
                                format!("t02-prefix{suffix}")
                            } else {
                                String::from("t02-prefix")
                            };
                            assert_eq!(
                                atom_key.as_bytes(),
                                expected_key.as_bytes(),
                                "atomic number {atomic_number}"
                            );
                            if valence_matches {
                                assert!(diagnostics.is_empty(), "atomic number {atomic_number}");
                            } else {
                                assert_mismatch_diagnostic(&diagnostics, atom_id);
                            }

                            atom_index += 1;
                            case_count += 1;
                        }
                    }
                }
            }
        }

        assert_eq!(case_count, 112);
    }

    #[test]
    fn cf3d_typ_t02_nonmember_is_not_handled_or_mutated() {
        let atom = fixed_atom(1, 3, 40_000);
        for tolerate_charge_mismatch in [false, true] {
            let mut atom_key = cosmolkit_model::PropertyText::from("t02-prefix");
            let mut diagnostics = Vec::new();
            assert!(!append_valence_only_charge_flag(
                &atom,
                3,
                &mut atom_key,
                tolerate_charge_mismatch,
                &mut diagnostics,
            ));
            assert_eq!((atom_key).as_bytes(), ("t02-prefix").as_bytes());
            assert!(diagnostics.is_empty());
        }
    }

    #[test]
    fn cf3d_typ_t03_all_elements_valence_tolerance_formal_charge_matrix() {
        let mut atom_index = 50_000;
        let mut case_count = 0;

        for &(atomic_number, required_valence) in UNSUFFIXED_MAIN_GROUP_CHARGES {
            for total_valence in [required_valence - 1, required_valence, required_valence + 1] {
                for tolerate_charge_mismatch in [false, true] {
                    for formal_charge_matches in [false, true] {
                        let formal_charge = if formal_charge_matches {
                            required_valence as i8
                        } else {
                            -1
                        };
                        let atom = fixed_atom(atomic_number, formal_charge, atom_index);
                        let atom_id = atom.id();
                        let mut atom_key = cosmolkit_model::PropertyText::from("t03-prefix");
                        let mut diagnostics = Vec::new();

                        assert!(check_unsuffixed_main_group_charge(
                            &atom,
                            total_valence,
                            &mut atom_key,
                            tolerate_charge_mismatch,
                            &mut diagnostics,
                        ));
                        assert_eq!(
                            (atom_key).as_bytes(),
                            ("t03-prefix").as_bytes(),
                            "atomic number {atomic_number}"
                        );
                        if total_valence == required_valence {
                            assert!(diagnostics.is_empty(), "atomic number {atomic_number}");
                        } else {
                            assert_mismatch_diagnostic(&diagnostics, atom_id);
                        }

                        atom_index += 1;
                        case_count += 1;
                    }
                }
            }
        }

        assert_eq!(case_count, 24);
    }

    #[test]
    fn cf3d_typ_t03_nonmember_is_not_handled_or_mutated() {
        let atom = fixed_atom(1, 0, 60_000);
        for tolerate_charge_mismatch in [false, true] {
            let mut atom_key = cosmolkit_model::PropertyText::from("t03-prefix");
            let mut diagnostics = Vec::new();
            assert!(!check_unsuffixed_main_group_charge(
                &atom,
                3,
                &mut atom_key,
                tolerate_charge_mismatch,
                &mut diagnostics,
            ));
            assert_eq!((atom_key).as_bytes(), ("t03-prefix").as_bytes());
            assert!(diagnostics.is_empty());
        }
    }

    #[test]
    fn cf3d_typ_t04_all_valence_tolerance_formal_charge_combinations() {
        let mut atom_index = 70_000;
        let mut case_count = 0;

        for total_valence in [0, 2, 3, 4, 5, 6] {
            for tolerate_charge_mismatch in [false, true] {
                for formal_charge in [0, 3, 5] {
                    let atom = fixed_atom(15, formal_charge, atom_index);
                    let atom_id = atom.id();
                    let mut atom_key = cosmolkit_model::PropertyText::from("t04-prefix");
                    let mut diagnostics = Vec::new();

                    assert!(append_phosphorus_charge_flag(
                        &atom,
                        total_valence,
                        &mut atom_key,
                        tolerate_charge_mismatch,
                        &mut diagnostics,
                    ));

                    let expected_suffix = match total_valence {
                        3 => "+3",
                        5 => "+5",
                        _ if tolerate_charge_mismatch => "+5",
                        _ => "",
                    };
                    assert_eq!(
                        (atom_key).as_bytes(),
                        (format!("t04-prefix{expected_suffix}")).as_bytes()
                    );
                    if matches!(total_valence, 3 | 5) {
                        assert!(diagnostics.is_empty());
                    } else {
                        assert_mismatch_diagnostic(&diagnostics, atom_id);
                    }

                    atom_index += 1;
                    case_count += 1;
                }
            }
        }

        assert_eq!(case_count, 36);
    }

    #[test]
    fn cf3d_typ_t05_all_hybridization_valence_tolerance_and_formal_charge_values() {
        let hybridizations = [
            Hybridization::Unspecified,
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp2,
            Hybridization::Sp3,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Other,
        ];
        let mut atom_index = 80_000;
        let mut case_count = 0;

        for hybridization in hybridizations {
            for total_valence in [0, 2, 3, 4, 5, 6, 7] {
                for tolerate_charge_mismatch in [false, true] {
                    for formal_charge in i8::MIN..=i8::MAX {
                        let element = Element::from_atomic_number(16)
                            .expect("fixed test atomic number is in the model range");
                        let atom = Atom::from_spec(
                            AtomId::new(atom_index),
                            AtomSpec::new(element)
                                .with_formal_charge(formal_charge)
                                .with_hybridization(hybridization),
                        );
                        let atom_id = atom.id();
                        let mut atom_key = cosmolkit_model::PropertyText::from("t05-prefix");
                        let mut diagnostics = Vec::new();

                        assert!(append_sulfur_charge_flag(
                            &atom,
                            total_valence,
                            &mut atom_key,
                            tolerate_charge_mismatch,
                            &mut diagnostics,
                        ));

                        let expected_suffix = if hybridization == Hybridization::Sp2 {
                            ""
                        } else {
                            match total_valence {
                                2 => "+2",
                                4 => "+4",
                                6 => "+6",
                                _ if tolerate_charge_mismatch => "+6",
                                _ => "",
                            }
                        };
                        assert_eq!(
                            (atom_key).as_bytes(),
                            (format!("t05-prefix{expected_suffix}")).as_bytes(),
                            "hybridization={hybridization:?}, valence={total_valence}, tolerance={tolerate_charge_mismatch}, charge={formal_charge}"
                        );

                        let expects_error = hybridization != Hybridization::Sp2
                            && !matches!(total_valence, 2 | 4 | 6);
                        if expects_error {
                            assert_mismatch_diagnostic(&diagnostics, atom_id);
                        } else {
                            assert!(diagnostics.is_empty());
                        }

                        atom_index += 1;
                        case_count += 1;
                    }
                }
            }
        }

        assert_eq!(case_count, 32_256);
    }

    #[test]
    fn cf3d_typ_t05_non_sulfur_is_not_handled_or_mutated() {
        let atom = fixed_atom(6, 0, 120_000);
        let mut atom_key = cosmolkit_model::PropertyText::from("t05-non-sulfur-prefix");
        let mut diagnostics = Vec::new();

        assert!(!append_sulfur_charge_flag(
            &atom,
            5,
            &mut atom_key,
            true,
            &mut diagnostics,
        ));
        assert_eq!((atom_key).as_bytes(), ("t05-non-sulfur-prefix").as_bytes());
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t06_exact_key_rewrites_and_ignored_charge_inputs() {
        let mut atom_id_value = 0;
        for &(input_key, tolerated_output_key) in RHENIUM_KEY_CASES {
            for tolerate_charge_mismatch in [false, true] {
                for &total_valence in REPRESENTATIVE_TOTAL_VALENCES {
                    for &formal_charge in REPRESENTATIVE_FORMAL_CHARGES {
                        let atom = fixed_atom(75, formal_charge, atom_id_value);
                        let atom_id = atom.id();
                        atom_id_value += 1;
                        let mut atom_key = cosmolkit_model::PropertyText::from(input_key);
                        let mut diagnostics = Vec::new();

                        assert!(rewrite_rhenium_charge_flag(
                            &atom,
                            total_valence,
                            &mut atom_key,
                            tolerate_charge_mismatch,
                            &mut diagnostics,
                        ));

                        let expected_key = if tolerate_charge_mismatch {
                            tolerated_output_key
                        } else {
                            input_key
                        };
                        assert_eq!((atom_key).as_bytes(), (expected_key).as_bytes());
                        assert_mismatch_diagnostic(&diagnostics, atom_id);
                    }
                }
            }
        }
    }

    #[test]
    fn cf3d_typ_t06_non_rhenium_is_not_handled_or_mutated() {
        let atom = fixed_atom(74, 5, 1_000);
        let mut atom_key = cosmolkit_model::PropertyText::from("Re6");
        let mut diagnostics = Vec::new();

        assert!(!rewrite_rhenium_charge_flag(
            &atom,
            6,
            &mut atom_key,
            true,
            &mut diagnostics,
        ));
        assert_eq!((atom_key).as_bytes(), ("Re6").as_bytes());
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t07_all_lanthanides_valence_tolerance_and_boundaries() {
        let mut atom_index = 0;
        let mut case_count = 0;
        for atomic_number in 56..=72 {
            for total_valence in [0, 3, 6, 7] {
                for tolerate_charge_mismatch in [false, true] {
                    let atom = fixed_atom(atomic_number, 0, atom_index);
                    let atom_id = atom.id();
                    atom_index += 1;
                    let mut atom_key = cosmolkit_model::PropertyText::from("t07-prefix");
                    let mut diagnostics = Vec::new();

                    let is_lanthanide = (57..=71).contains(&atomic_number);
                    assert_eq!(
                        append_lanthanide_charge_flag(
                            &atom,
                            total_valence,
                            &mut atom_key,
                            tolerate_charge_mismatch,
                            &mut diagnostics,
                        ),
                        is_lanthanide,
                        "atomic_number={atomic_number}, valence={total_valence}, tolerance={tolerate_charge_mismatch}"
                    );

                    let expected_suffix =
                        if is_lanthanide && (total_valence == 6 || tolerate_charge_mismatch) {
                            "+3"
                        } else {
                            ""
                        };
                    assert_eq!(
                        (atom_key).as_bytes(),
                        (format!("t07-prefix{expected_suffix}")).as_bytes()
                    );

                    if is_lanthanide && total_valence != 6 {
                        assert_mismatch_diagnostic(&diagnostics, atom_id);
                    } else {
                        assert!(diagnostics.is_empty());
                    }
                    case_count += 1;
                }
            }
        }
        assert_eq!(case_count, 136);
    }

    #[test]
    fn cf3d_typ_t12_source_reference_values_reversal_and_pair_scope() {
        let topology = topology_with_bond(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(7, 1, Hybridization::Sp3, false, None),
                label_atom(0, 2, Hybridization::Other, false, Some("Q")),
            ],
            BondOrder::Single,
        );
        let topology_before = topology.clone();
        let total_valences = [4, 3, 0];
        let total_valences_before = total_valences;
        let atom_has_conjugated_bond = [false; 3];
        let conjugation_before = atom_has_conjugated_bond;
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();

        // Pinned RDKit testUFFForceField.cpp fixes Csp3-Nsp3 at r0=1.451071
        // and kb=1057.27; these literals are independent of the Rust helper.
        let forward = get_uff_bond_stretch_params(
            &topology,
            0,
            1,
            &total_valences,
            &atom_has_conjugated_bond,
            &params,
            &mut diagnostics,
        )
        .expect("valid source-shaped bond query")
        .expect("both endpoint types exist");
        assert!((forward.r0 - 1.451_071).abs() < 1.0e-5);
        assert!((forward.kb - 1_057.27).abs() < 1.0e-2);
        assert!(diagnostics.is_empty());

        let reverse = get_uff_bond_stretch_params(
            &topology,
            1,
            0,
            &total_valences,
            &atom_has_conjugated_bond,
            &params,
            &mut diagnostics,
        )
        .expect("reversed valid source-shaped bond query")
        .expect("both endpoint types exist in reverse order");
        assert!((reverse.r0 - 1.451_071).abs() < 1.0e-5);
        assert!((reverse.kb - 1_057.27).abs() < 1.0e-2);
        assert!(diagnostics.is_empty());

        // The disconnected unknown type must not be visited by this local pair
        // query; source getUFFBondStretchParams labels only its two endpoints.
        assert_eq!(topology, topology_before);
        assert_eq!(total_valences, total_valences_before);
        assert_eq!(atom_has_conjugated_bond, conjugation_before);
    }

    #[test]
    fn cf3d_typ_t12_absent_bond_and_invalid_indices_keep_source_boundaries() {
        let atoms = vec![
            label_atom(0, 0, Hybridization::Other, false, Some("Q")),
            label_atom(80, 1, Hybridization::Sp3, false, None),
        ];
        let no_bond = topology_from_atoms(atoms.clone());
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let total_valences = [0, 2];
        let conjugation = [false; 2];
        let mut diagnostics = Vec::new();

        assert_eq!(
            get_uff_bond_stretch_params(
                &no_bond,
                0,
                1,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            )
            .expect("valid endpoints without a bond"),
            None
        );
        assert!(diagnostics.is_empty());

        let with_bond = topology_with_bond(atoms, BondOrder::Single);
        assert_eq!(
            get_uff_bond_stretch_params(
                &with_bond,
                2,
                1,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::AtomIdOutOfBounds {
                atom_id: AtomId::new(2),
                atom_count: 2,
            })
        );
        assert_eq!(
            get_uff_bond_stretch_params(
                &with_bond,
                0,
                2,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::AtomIdOutOfBounds {
                atom_id: AtomId::new(2),
                atom_count: 2,
            })
        );
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t12_missing_endpoint_types_preserve_ordered_short_circuit() {
        let topology = topology_with_bond(
            vec![
                label_atom(0, 0, Hybridization::Other, false, Some("Q")),
                label_atom(80, 1, Hybridization::Sp3, false, None),
            ],
            BondOrder::Single,
        );
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let total_valences = [0, 2];
        let conjugation = [false; 2];
        let mut diagnostics = Vec::new();

        assert_eq!(
            get_uff_bond_stretch_params(
                &topology,
                0,
                1,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            )
            .expect("first endpoint type miss is source false"),
            None
        );
        // The second endpoint would emit this source warning. Its absence
        // proves the first missing table entry short-circuits that label call.
        assert!(diagnostics.is_empty());

        assert_eq!(
            get_uff_bond_stretch_params(
                &topology,
                1,
                0,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            )
            .expect("second endpoint type miss is source false"),
            None
        );
        assert_eq!(
            diagnostics,
            vec![UffTypingDiagnostic {
                atom_id: Some(AtomId::new(1)),
                kind: UffTypingDiagnosticKind::Warning,
                message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
            }]
        );
    }

    #[test]
    fn cf3d_typ_t12_every_modeled_bond_order_and_fixed_numeric_reference() {
        let positive_cases = [
            (BondOrder::Single, 1.514_000_000_000, 699.591_798_712_679),
            (BondOrder::Double, 1.374_216_612_462, 935.527_918_062_998),
            (BondOrder::Triple, 1.292_448_572_528, 1_124.559_714_235_118),
            (
                BondOrder::Quadruple,
                1.234_433_224_924,
                1_290.682_816_763_283,
            ),
            (
                BondOrder::Quintuple,
                1.189_433_025_277,
                1_442.787_452_823_364,
            ),
            (
                BondOrder::Hextuple,
                1.152_665_184_990,
                1_585.304_919_151_708,
            ),
            (
                BondOrder::OneAndHalf,
                1.432_231_960_066,
                826.384_671_555_230,
            ),
            (
                BondOrder::TwoAndHalf,
                1.329_216_412_815,
                1_033.796_957_118_917,
            ),
            (
                BondOrder::ThreeAndHalf,
                1.261_361_806_511,
                1_209.771_377_931_984,
            ),
            (
                BondOrder::FourAndHalf,
                1.210_680_532_595,
                1_368.149_820_204_275,
            ),
            (
                BondOrder::FiveAndHalf,
                1.170_212_316_928,
                1_515.054_795_412_171,
            ),
            (BondOrder::Aromatic, 1.432_231_960_066, 826.384_671_555_230),
            (BondOrder::DativeOne, 1.514_000_000_000, 699.591_798_712_679),
            (BondOrder::Dative, 1.514_000_000_000, 699.591_798_712_679),
        ];
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let total_valences = [4, 4];
        let conjugation = [false; 2];

        // These fixed references use the pinned C_3 row (r1=0.757, Z1=1.912,
        // Xi=5.343) and the verbatim RDKit BondStretch.cpp equations.
        for (order, expected_r0, expected_kb) in positive_cases {
            let topology = topology_with_bond(
                vec![
                    label_atom(6, 0, Hybridization::Sp3, false, None),
                    label_atom(6, 1, Hybridization::Sp3, false, None),
                ],
                order,
            );
            let mut diagnostics = Vec::new();
            let actual = get_uff_bond_stretch_params(
                &topology,
                0,
                1,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            )
            .expect("modeled positive bond order is structurally valid")
            .expect("C_3 endpoint parameters exist");
            assert!((actual.r0 - expected_r0).abs() < 1.0e-9, "order={order:?}");
            assert!((actual.kb - expected_kb).abs() < 1.0e-6, "order={order:?}");
            assert!(diagnostics.is_empty(), "order={order:?}");
        }

        for order in [
            BondOrder::Unspecified,
            BondOrder::Ionic,
            BondOrder::Hydrogen,
            BondOrder::Zero,
        ] {
            let topology = topology_with_bond(
                vec![
                    label_atom(6, 0, Hybridization::Sp3, false, None),
                    label_atom(6, 1, Hybridization::Sp3, false, None),
                ],
                order,
            );
            let mut diagnostics = Vec::new();
            assert!(matches!(
                get_uff_bond_stretch_params(
                    &topology,
                    0,
                    1,
                    &total_valences,
                    &conjugation,
                    &params,
                    &mut diagnostics,
                ),
                Err(UffTypingError::BondMath(BondMathError::InvalidBondOrder {
                    bond_order
                })) if bond_order == 0.0
            ));
        }

        for order in [
            BondOrder::ThreeCenter,
            BondOrder::DativeLeft,
            BondOrder::DativeRight,
            BondOrder::Other,
        ] {
            let topology = topology_with_bond(
                vec![
                    label_atom(6, 0, Hybridization::Sp3, false, None),
                    label_atom(6, 1, Hybridization::Sp3, false, None),
                ],
                order,
            );
            let mut diagnostics = Vec::new();
            assert!(matches!(
                get_uff_bond_stretch_params(
                    &topology,
                    0,
                    1,
                    &total_valences,
                    &conjugation,
                    &params,
                    &mut diagnostics,
                ),
                Err(UffTypingError::CoreValence(
                    cosmolkit_core::ValenceError::BadBondType { order: actual, .. }
                )) if actual == order
            ));
        }
    }

    #[test]
    fn cf3d_typ_t12_prepared_state_length_errors_are_typed() {
        let topology = topology_with_bond(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(6, 1, Hybridization::Sp3, false, None),
            ],
            BondOrder::Single,
        );
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();

        assert!(matches!(
            get_uff_bond_stretch_params(
                &topology,
                0,
                1,
                &[4],
                &[false, false],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: 2,
                actual: 1,
            })
        ));
        assert!(matches!(
            get_uff_bond_stretch_params(
                &topology,
                0,
                1,
                &[4, 4],
                &[false],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::ConjugatedBondPresence,
                expected: 2,
                actual: 1,
            })
        ));
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t13_source_reference_degrees_and_force_constant() {
        let topology = topology_with_angle(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(6, 1, Hybridization::Sp3, false, None),
                label_atom(7, 2, Hybridization::Sp3, false, None),
            ],
            BondOrder::Single,
            BondOrder::Single,
        );
        let topology_before = topology.clone();
        let total_valences = [4, 4, 3];
        let conjugation = [false; 3];
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();

        let angle = get_uff_angle_bend_params(
            &topology,
            0,
            1,
            2,
            &total_valences,
            &conjugation,
            &params,
            &mut diagnostics,
        )
        .expect("valid source-shaped C_3 angle query")
        .expect("all three C_3 parameter rows exist");

        // Pinned testUFFHelpers.cpp:1054-1057 uses atoms 6-7-8 in
        // c1ccccc1CCNN (C_3-C_3-N_3) and fixes the result at
        // round(ka*1000)=303297 and round(theta0*1000)=109470.
        assert_eq!((angle.ka * 1000.0).round() as i64, 303_297);
        assert_eq!((angle.theta0 * 1000.0).round() as i64, 109_470);
        assert!(diagnostics.is_empty());
        assert_eq!(topology, topology_before);
        assert_eq!(total_valences, [4, 4, 3]);
        assert_eq!(conjugation, [false; 3]);
    }

    #[test]
    fn cf3d_typ_t13_asymmetric_params_and_endpoint_reversal() {
        let topology = topology_with_angle(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(7, 1, Hybridization::Sp3, false, None),
                label_atom(8, 2, Hybridization::Sp3, false, None),
            ],
            BondOrder::OneAndHalf,
            BondOrder::Double,
        );
        let total_valences = [4, 3, 2];
        let conjugation = [false; 3];
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();

        let forward = get_uff_angle_bend_params(
            &topology,
            0,
            1,
            2,
            &total_valences,
            &conjugation,
            &params,
            &mut diagnostics,
        )
        .expect("valid asymmetric C_3-N_3-O_3 path")
        .expect("all three parameter rows exist");
        let reverse = get_uff_angle_bend_params(
            &topology,
            2,
            1,
            0,
            &total_valences,
            &conjugation,
            &params,
            &mut diagnostics,
        )
        .expect("valid reversed asymmetric path")
        .expect("all reversed parameter rows exist");

        // Fixed from pinned default_params.tsv C_3/N_3/O_3 r1, theta0, Z1,
        // and Xi rows with BondStretch.cpp and AngleBend.cpp operation order.
        assert!((forward.theta0 - 106.7).abs() < 1.0e-12);
        assert!((forward.ka - 433.523_080_654_281_84).abs() < 1.0e-9);
        assert!((reverse.theta0 - 106.7).abs() < 1.0e-12);
        assert!((reverse.ka - 433.523_080_654_281_84).abs() < 1.0e-9);
        assert!((forward.ka - reverse.ka).abs() < 1.0e-12);
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t13_missing_first_and_second_bonds_short_circuit_labels() {
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let total_valences = [2, 4, 2];
        let conjugation = [false; 3];
        let mut diagnostics = Vec::new();

        let first_bond_missing = topology_from_atoms(vec![
            label_atom(80, 0, Hybridization::Sp3, false, None),
            label_atom(6, 1, Hybridization::Sp3, false, None),
            label_atom(8, 2, Hybridization::Sp3, false, None),
        ]);
        assert_eq!(
            get_uff_angle_bend_params(
                &first_bond_missing,
                0,
                1,
                2,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            )
            .expect("missing first bond is the source false result"),
            None
        );
        // The first edge miss returns before label access, even though atom 0
        // would emit the source SP-hybridization warning and idx3 is invalid.
        assert_eq!(
            get_uff_angle_bend_params(
                &first_bond_missing,
                0,
                1,
                99,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            )
            .expect("later idx3 is not checked after the first edge miss"),
            None
        );
        assert!(diagnostics.is_empty());

        let second_bond_missing = topology_with_bond(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(80, 1, Hybridization::Sp3, false, None),
                label_atom(8, 2, Hybridization::Sp3, false, None),
            ],
            BondOrder::Single,
        );
        assert_eq!(
            get_uff_angle_bend_params(
                &second_bond_missing,
                0,
                1,
                2,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            )
            .expect("missing second bond is the source false result"),
            None
        );
        // Atom 1 would warn if its label lookup preceded the second bond test.
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t13_each_missing_parameter_slot_stops_later_rows() {
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let total_valences = [4, 3, 2];
        let conjugation = [false; 3];

        for missing_slot in 0..3 {
            let atoms = match missing_slot {
                0 => vec![
                    label_atom(0, 0, Hybridization::Other, false, Some("Q")),
                    label_atom(80, 1, Hybridization::Sp3, false, None),
                    label_atom(8, 2, Hybridization::Sp3, false, None),
                ],
                1 => vec![
                    label_atom(6, 0, Hybridization::Sp3, false, None),
                    label_atom(0, 1, Hybridization::Other, false, Some("Q")),
                    label_atom(80, 2, Hybridization::Sp3, false, None),
                ],
                2 => vec![
                    label_atom(6, 0, Hybridization::Sp3, false, None),
                    label_atom(7, 1, Hybridization::Sp3, false, None),
                    label_atom(0, 2, Hybridization::Other, false, Some("Q")),
                ],
                _ => unreachable!("fixed slot range"),
            };
            let topology = topology_with_angle(atoms, BondOrder::Single, BondOrder::Single);
            let mut diagnostics = Vec::new();

            assert_eq!(
                get_uff_angle_bend_params(
                    &topology,
                    0,
                    1,
                    2,
                    &total_valences,
                    &conjugation,
                    &params,
                    &mut diagnostics,
                )
                .expect("a missing table row is source false, not an error"),
                None,
                "missing slot {missing_slot}"
            );
            // Later Hg rows would each emit a warning; an earlier missing
            // parameter must retain source row order and suppress them.
            assert!(diagnostics.is_empty(), "missing slot {missing_slot}");
        }
    }

    #[test]
    fn cf3d_typ_t13_diagnostic_prefix_and_inputs_remain_unchanged() {
        let topology = topology_with_angle(
            vec![
                label_atom(80, 0, Hybridization::Sp3, false, None),
                label_atom(0, 1, Hybridization::Other, false, Some("Q")),
                label_atom(7, 2, Hybridization::Sp3, false, None),
            ],
            BondOrder::Single,
            BondOrder::Single,
        );
        let topology_before = topology.clone();
        let total_valences = [2, 0, 3];
        let total_valences_before = total_valences;
        let conjugation = [false; 3];
        let conjugation_before = conjugation;
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();

        assert_eq!(
            get_uff_angle_bend_params(
                &topology,
                0,
                1,
                2,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            )
            .expect("source nonfatal label diagnostics do not become errors"),
            None
        );
        assert_eq!(
            diagnostics,
            vec![UffTypingDiagnostic {
                atom_id: Some(AtomId::new(0)),
                kind: UffTypingDiagnosticKind::Warning,
                message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
            }]
        );
        assert_eq!(topology, topology_before);
        assert_eq!(total_valences, total_valences_before);
        assert_eq!(conjugation, conjugation_before);
    }

    #[test]
    fn cf3d_typ_t13_order_errors_remain_typed_and_sequential() {
        let atoms = vec![
            label_atom(6, 0, Hybridization::Sp3, false, None),
            label_atom(6, 1, Hybridization::Sp3, false, None),
            label_atom(6, 2, Hybridization::Sp3, false, None),
        ];
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let total_valences = [4; 3];
        let conjugation = [false; 3];

        for (order12, order23, expected_error_order) in [
            (BondOrder::Zero, BondOrder::Single, 0.0),
            (BondOrder::Single, BondOrder::Unspecified, 0.0),
        ] {
            let topology = topology_with_angle(atoms.clone(), order12, order23);
            let mut diagnostics = Vec::new();
            assert!(matches!(
                get_uff_angle_bend_params(
                    &topology,
                    0,
                    1,
                    2,
                    &total_valences,
                    &conjugation,
                    &params,
                    &mut diagnostics,
                ),
                Err(UffTypingError::BondMath(BondMathError::InvalidBondOrder {
                    bond_order
                })) if bond_order == expected_error_order
            ));
            assert!(diagnostics.is_empty());
        }

        for (order12, order23) in [
            (BondOrder::Other, BondOrder::Single),
            (BondOrder::Single, BondOrder::Other),
        ] {
            let topology = topology_with_angle(atoms.clone(), order12, order23);
            let mut diagnostics = Vec::new();
            assert!(matches!(
                get_uff_angle_bend_params(
                    &topology,
                    0,
                    1,
                    2,
                    &total_valences,
                    &conjugation,
                    &params,
                    &mut diagnostics,
                ),
                Err(UffTypingError::CoreValence(
                    cosmolkit_core::ValenceError::BadBondType { order, .. }
                )) if order == if order12 == BondOrder::Other { order12 } else { order23 }
            ));
            assert!(diagnostics.is_empty());
        }
    }

    #[test]
    fn cf3d_typ_t13_prepared_state_lengths_and_late_index_are_typed() {
        let topology = topology_with_angle(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(7, 1, Hybridization::Sp3, false, None),
                label_atom(8, 2, Hybridization::Sp3, false, None),
            ],
            BondOrder::Single,
            BondOrder::Single,
        );
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();

        assert!(matches!(
            get_uff_angle_bend_params(
                &topology,
                0,
                1,
                2,
                &[4, 3],
                &[false; 3],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: 3,
                actual: 2,
            })
        ));
        assert!(matches!(
            get_uff_angle_bend_params(
                &topology,
                0,
                1,
                2,
                &[4, 3, 2],
                &[false; 2],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::ConjugatedBondPresence,
                expected: 3,
                actual: 2,
            })
        ));
        assert_eq!(
            get_uff_angle_bend_params(
                &topology,
                0,
                1,
                3,
                &[4, 3, 2],
                &[false; 3],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::AtomIdOutOfBounds {
                atom_id: AtomId::new(3),
                atom_count: 3,
            })
        );
        assert!(diagnostics.is_empty());
    }

    fn atomic_params_with_v1(v1: f64) -> AtomicParams {
        AtomicParams {
            r1: 0.0,
            theta0: 0.0,
            x1: 0.0,
            d1: 0.0,
            zeta: 0.0,
            z1: 0.0,
            v1,
            u1: 0.0,
            gmp_xi: 0.0,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn atomic_params_with_u1(u1: f64) -> AtomicParams {
        AtomicParams {
            r1: 0.0,
            theta0: 0.0,
            x1: 0.0,
            d1: 0.0,
            zeta: 0.0,
            z1: 0.0,
            v1: 0.0,
            u1,
            gmp_xi: 0.0,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn assert_t14_amplitude(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1.0e-12,
            "expected amplitude {expected}, got {actual}"
        );
    }

    fn assert_t15_amplitude(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1.0e-12,
            "expected amplitude {expected}, got {actual}"
        );
    }

    fn expected_t15_equation17(bond_order: f64) -> f64 {
        // Fixed source equation17 reference with U1 values 4.0 and 9.0.
        5.0 * (4.0_f64 * 9.0).sqrt() * (1.0 + 4.18 * bond_order.ln())
    }

    #[test]
    fn cf3d_typ_t14_general_branch_uses_asymmetric_v1_product() {
        let params0 = atomic_params_with_v1(4.0);
        let params1 = atomic_params_with_v1(9.0);

        for (at_num0, at_num1, first, second) in [
            (6, 7, &params0, &params1),
            (7, 6, &params1, &params0),
            (8, 6, &params0, &params1),
            (6, 8, &params1, &params0),
        ] {
            assert_t14_amplitude(
                torsion_sp3_amplitude(1.0, at_num0, at_num1, first, second),
                6.0,
            );
        }
    }

    #[test]
    fn cf3d_typ_t14_group_six_override_covers_every_ordered_pair() {
        const GROUP_SIX: [i32; 5] = [8, 16, 34, 52, 84];
        let params0 = atomic_params_with_v1(4.0);
        let params1 = atomic_params_with_v1(9.0);

        for at_num0 in GROUP_SIX {
            for at_num1 in GROUP_SIX {
                let v2: f64 = if at_num0 == 8 { 2.0 } else { 6.8 };
                let v3: f64 = if at_num1 == 8 { 2.0 } else { 6.8 };
                assert_t14_amplitude(
                    torsion_sp3_amplitude(1.05, at_num0, at_num1, &params0, &params1),
                    (v2 * v3).sqrt(),
                );
            }
        }
    }

    #[test]
    fn cf3d_typ_t14_cast_boundary_and_non_group_pairs_preserve_general_value() {
        const GROUP_SIX_PAIR: (i32, i32) = (8, 16);
        let params0 = atomic_params_with_v1(4.0);
        let params1 = atomic_params_with_v1(9.0);

        for (bond_order, special_case) in [
            (0.999, false),
            (1.0, true),
            (1.05, true),
            (1.099, true),
            (1.1, false),
        ] {
            let expected = if special_case {
                (2.0_f64 * 6.8).sqrt()
            } else {
                6.0
            };
            assert_t14_amplitude(
                torsion_sp3_amplitude(
                    bond_order,
                    GROUP_SIX_PAIR.0,
                    GROUP_SIX_PAIR.1,
                    &params0,
                    &params1,
                ),
                expected,
            );
        }

        for group_num in [8, 16, 34, 52, 84] {
            for non_group_num in [6, 7, 9, 17] {
                for (at_num0, at_num1) in [(group_num, non_group_num), (non_group_num, group_num)] {
                    assert_t14_amplitude(
                        torsion_sp3_amplitude(1.05, at_num0, at_num1, &params0, &params1),
                        6.0,
                    );
                }
            }
        }

        for (at_num0, at_num1) in [(6, 7), (7, 6)] {
            assert_t14_amplitude(
                torsion_sp3_amplitude(1.05, at_num0, at_num1, &params0, &params1),
                6.0,
            );
        }
    }

    #[test]
    fn cf3d_typ_t15_covers_hybridization_membership_terminal_and_cast_product() {
        const GROUP_SIX: [i32; 5] = [8, 16, 34, 52, 84];
        const NON_GROUP_SIX: [i32; 2] = [6, 7];
        const BOND_ORDERS: [f64; 5] = [0.999, 1.0, 1.05, 1.099, 1.1];
        let params0 = atomic_params_with_u1(4.0);
        let params1 = atomic_params_with_u1(9.0);

        for (hyb0, hyb1) in [
            (Hybridization::Sp3, Hybridization::Sp2),
            (Hybridization::Sp2, Hybridization::Sp3),
        ] {
            for at_num0 in GROUP_SIX.into_iter().chain(NON_GROUP_SIX) {
                for at_num1 in GROUP_SIX.into_iter().chain(NON_GROUP_SIX) {
                    let at_num0_is_group_six = GROUP_SIX.contains(&at_num0);
                    let at_num1_is_group_six = GROUP_SIX.contains(&at_num1);
                    let group_six_sp3_non_group_six_sp2 = if hyb0 == Hybridization::Sp3 {
                        at_num0_is_group_six && !at_num1_is_group_six
                    } else {
                        at_num1_is_group_six && !at_num0_is_group_six
                    };

                    for bond_order in BOND_ORDERS {
                        let source_single_class = (bond_order * 10.0) as i32 == 10;
                        for has_sp2 in [false, true] {
                            let expected = if source_single_class && group_six_sp3_non_group_six_sp2
                            {
                                expected_t15_equation17(bond_order)
                            } else if source_single_class && has_sp2 {
                                2.0
                            } else {
                                1.0
                            };

                            assert_t15_amplitude(
                                torsion_mixed_amplitude(
                                    bond_order, at_num0, at_num1, hyb0, hyb1, &params0, &params1,
                                    has_sp2,
                                ),
                                expected,
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn cf3d_typ_t15_group_six_equation17_precedes_terminal_sp2_override() {
        const EXPECTED_EQUATION17: f64 = 30.0;
        let params0 = atomic_params_with_u1(4.0);
        let params1 = atomic_params_with_u1(9.0);

        for (at_num0, at_num1, hyb0, hyb1) in [
            (8, 6, Hybridization::Sp3, Hybridization::Sp2),
            (6, 8, Hybridization::Sp2, Hybridization::Sp3),
        ] {
            for has_sp2 in [false, true] {
                assert_t15_amplitude(
                    torsion_mixed_amplitude(
                        1.0, at_num0, at_num1, hyb0, hyb1, &params0, &params1, has_sp2,
                    ),
                    EXPECTED_EQUATION17,
                );
            }
        }
    }

    fn t16_path(
        center0: (u8, Hybridization),
        center1: (u8, Hybridization),
        terminal_hybridizations: [Hybridization; 2],
        edges: &[(usize, usize, BondOrder)],
    ) -> TopologyBlock {
        topology_with_bonds(
            vec![
                label_atom(6, 0, terminal_hybridizations[0], false, None),
                label_atom(center0.0, 1, center0.1, false, None),
                label_atom(center1.0, 2, center1.1, false, None),
                label_atom(6, 3, terminal_hybridizations[1], false, None),
            ],
            edges,
        )
    }

    fn t16_query_value(
        topology: &TopologyBlock,
        indices: [usize; 4],
        total_valences: &[i32],
        conjugation: &[bool],
        params: &ParamCollection,
        diagnostics: &mut Vec<UffTypingDiagnostic>,
    ) -> Result<Option<f64>, UffTypingError> {
        get_uff_torsion_params(
            topology,
            indices[0],
            indices[1],
            indices[2],
            indices[3],
            total_valences,
            conjugation,
            params,
            diagnostics,
        )
        .map(|value| value.map(|torsion| torsion.v))
    }

    fn t16_warning(atom_index: usize) -> UffTypingDiagnostic {
        UffTypingDiagnostic {
            atom_id: Some(AtomId::new(atom_index)),
            kind: UffTypingDiagnosticKind::Warning,
            message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
        }
    }

    fn t16_default_params() -> std::sync::Arc<ParamCollection> {
        ParamCollection::get_params("").expect("pinned default UFF parameter table")
    }

    #[test]
    fn cf3d_typ_t16_amplitude_dispatch_and_terminal_sp2_matrix() {
        const EDGES: &[(usize, usize, BondOrder)] = &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
        ];
        const TERMINAL_SP2: [(bool, bool); 4] =
            [(false, false), (true, false), (false, true), (true, true)];
        let params = t16_default_params();

        for (center0, center1, base, terminal_sp2_changes_amplitude) in [
            (
                (6, Hybridization::Sp3),
                (6, Hybridization::Sp3),
                2.119,
                false,
            ),
            (
                (6, Hybridization::Sp2),
                (6, Hybridization::Sp2),
                10.0,
                false,
            ),
            ((6, Hybridization::Sp2), (6, Hybridization::Sp3), 1.0, true),
            ((6, Hybridization::Sp3), (6, Hybridization::Sp2), 1.0, true),
        ] {
            for (first_terminal_sp2, last_terminal_sp2) in TERMINAL_SP2 {
                let topology = t16_path(
                    center0,
                    center1,
                    [
                        if first_terminal_sp2 {
                            Hybridization::Sp2
                        } else {
                            Hybridization::Sp3
                        },
                        if last_terminal_sp2 {
                            Hybridization::Sp2
                        } else {
                            Hybridization::Sp3
                        },
                    ],
                    EDGES,
                );
                let mut diagnostics = Vec::new();
                let value = t16_query_value(
                    &topology,
                    [0, 1, 2, 3],
                    &[0; 4],
                    &[false; 4],
                    &params,
                    &mut diagnostics,
                )
                .expect("valid source-shaped torsion query")
                .expect("both central parameter rows exist");
                let expected = if terminal_sp2_changes_amplitude
                    && (first_terminal_sp2 || last_terminal_sp2)
                {
                    2.0
                } else {
                    base
                };
                assert!(
                    (value - expected).abs() <= 1.0e-12,
                    "centers {center0:?}/{center1:?}, terminals {first_terminal_sp2}/{last_terminal_sp2}: expected {expected}, got {value}"
                );
                assert!(diagnostics.is_empty());
            }
        }
    }

    #[test]
    fn cf3d_typ_t16_group_six_overrides_and_non_single_control() {
        const TERMINAL_SP2: [(bool, bool); 4] =
            [(false, false), (true, false), (false, true), (true, true)];
        let params = t16_default_params();

        for (center0, center1) in [
            ((8, Hybridization::Sp3), (6, Hybridization::Sp2)),
            ((6, Hybridization::Sp2), (8, Hybridization::Sp3)),
        ] {
            for (first_terminal_sp2, last_terminal_sp2) in TERMINAL_SP2 {
                let topology = t16_path(
                    center0,
                    center1,
                    [
                        if first_terminal_sp2 {
                            Hybridization::Sp2
                        } else {
                            Hybridization::Sp3
                        },
                        if last_terminal_sp2 {
                            Hybridization::Sp2
                        } else {
                            Hybridization::Sp3
                        },
                    ],
                    &[
                        (0, 1, BondOrder::Single),
                        (1, 2, BondOrder::Single),
                        (2, 3, BondOrder::Single),
                    ],
                );
                let mut diagnostics = Vec::new();
                let value = t16_query_value(
                    &topology,
                    [0, 1, 2, 3],
                    &[0; 4],
                    &[false; 4],
                    &params,
                    &mut diagnostics,
                )
                .expect("valid group-six mixed torsion")
                .expect("both central parameter rows exist");
                // Default O_3 and C_2 rows both have U1=2.0, so source
                // equation17 on a single bond fixes this to 10.0. It takes
                // precedence over every terminal hasSP2 combination.
                assert!((value - 10.0).abs() <= 1.0e-12);
                assert!(diagnostics.is_empty());
            }
        }

        let oxygen_pair = t16_path(
            (8, Hybridization::Sp3),
            (8, Hybridization::Sp3),
            [Hybridization::Sp3; 2],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
            ],
        );
        let mut diagnostics = Vec::new();
        assert_eq!(
            t16_query_value(
                &oxygen_pair,
                [0, 1, 2, 3],
                &[0; 4],
                &[false; 4],
                &params,
                &mut diagnostics,
            )
            .expect("valid source group-six SP3/SP3 torsion"),
            Some(2.0)
        );
        assert!(diagnostics.is_empty());

        let non_single_group_pair = t16_path(
            (8, Hybridization::Sp3),
            (6, Hybridization::Sp2),
            [Hybridization::Sp2; 2],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
            ],
        );
        assert_eq!(
            t16_query_value(
                &non_single_group_pair,
                [0, 1, 2, 3],
                &[0; 4],
                &[false; 4],
                &params,
                &mut diagnostics,
            )
            .expect("double center bond is a valid query"),
            Some(1.0)
        );
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t16_source_res_assignments_and_missing_edges() {
        let params = t16_default_params();

        let first_edge_missing = t16_path(
            (6, Hybridization::Sp3),
            (6, Hybridization::Sp3),
            [Hybridization::Sp2, Hybridization::Sp3],
            &[(1, 2, BondOrder::Single), (2, 3, BondOrder::Single)],
        );
        let mut diagnostics = Vec::new();
        assert_eq!(
            t16_query_value(
                &first_edge_missing,
                [0, 1, 2, 99],
                &[0; 4],
                &[false; 4],
                &params,
                &mut diagnostics,
            )
            .expect("first-edge false short-circuits before idx4"),
            None
        );
        assert!(diagnostics.is_empty());

        let center_edge_missing = t16_path(
            (6, Hybridization::Sp3),
            (6, Hybridization::Sp3),
            [Hybridization::Sp3; 2],
            &[(0, 1, BondOrder::Single), (2, 3, BondOrder::Single)],
        );
        assert!(matches!(
            t16_query_value(
                &center_edge_missing,
                [0, 1, 2, 3],
                &[0; 4],
                &[false; 4],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::CenterBondMissingAfterSourceSuccess { idx2: 1, idx3: 2 })
        ));
        assert!(diagnostics.is_empty());

        let trailing_edge_missing = t16_path(
            (6, Hybridization::Sp3),
            (6, Hybridization::Sp2),
            [Hybridization::Sp3, Hybridization::Sp2],
            &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)],
        );
        // At i=2, the missing trailing edge sets res=false, then the present
        // central parameter overwrites it; source still returns amplitude 2.
        assert_eq!(
            t16_query_value(
                &trailing_edge_missing,
                [0, 1, 2, 3],
                &[0; 4],
                &[false; 4],
                &params,
                &mut diagnostics,
            )
            .expect("source parameter assignment overwrites trailing edge miss"),
            Some(2.0)
        );
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t16_rejects_every_non_sp2_sp3_central_hybridization() {
        const EDGES: &[(usize, usize, BondOrder)] = &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
        ];
        const REJECTED: [Hybridization; 7] = [
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Unspecified,
            Hybridization::Other,
        ];
        let params = t16_default_params();

        for center_slot in 0..2 {
            for hybridization in REJECTED {
                let center0 = if center_slot == 0 {
                    (80, hybridization)
                } else {
                    (6, Hybridization::Sp3)
                };
                let center1 = if center_slot == 1 {
                    (80, hybridization)
                } else {
                    (6, Hybridization::Sp3)
                };
                let topology = t16_path(center0, center1, [Hybridization::Sp3; 2], EDGES);
                let total_valences = if center_slot == 0 {
                    [0, 2, 4, 0]
                } else {
                    [0, 4, 2, 0]
                };
                let mut diagnostics = Vec::new();
                assert_eq!(
                    t16_query_value(
                        &topology,
                        [0, 1, 2, 3],
                        &total_valences,
                        &[false; 4],
                        &params,
                        &mut diagnostics,
                    )
                    .expect("Hg label has a source parameter row"),
                    None,
                    "center slot {center_slot}, hybridization {hybridization:?}"
                );
                let expected_diagnostics = if hybridization == Hybridization::Sp {
                    Vec::new()
                } else {
                    vec![t16_warning(center_slot + 1)]
                };
                assert_eq!(
                    diagnostics, expected_diagnostics,
                    "center slot {center_slot}, hybridization {hybridization:?}"
                );
            }
        }
    }

    #[test]
    fn cf3d_typ_t16_missing_parameters_diagnostics_and_typed_inputs() {
        const EDGES: &[(usize, usize, BondOrder)] = &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
        ];
        let params = t16_default_params();

        let missing_first_param = t16_path(
            (0, Hybridization::Other),
            (80, Hybridization::Sp2),
            [Hybridization::Sp3; 2],
            EDGES,
        );
        let mut diagnostics = Vec::new();
        assert_eq!(
            t16_query_value(
                &missing_first_param,
                [0, 1, 2, 3],
                &[0, 0, 2, 0],
                &[false; 4],
                &params,
                &mut diagnostics,
            )
            .expect("missing first central key is the source false result"),
            None
        );
        // The dummy Q row has no UFF entry; idx2's Hg label is never visited.
        assert!(diagnostics.is_empty());

        let missing_second_param = t16_path(
            (80, Hybridization::Sp2),
            (0, Hybridization::Other),
            [Hybridization::Sp3; 2],
            EDGES,
        );
        assert_eq!(
            t16_query_value(
                &missing_second_param,
                [0, 1, 2, 3],
                &[0, 2, 0, 0],
                &[false; 4],
                &params,
                &mut diagnostics,
            )
            .expect("missing second central key is the source false result"),
            None
        );
        assert_eq!(diagnostics, vec![t16_warning(1)]);

        let central_hg_pair = t16_path(
            (80, Hybridization::Sp2),
            (80, Hybridization::Sp2),
            [Hybridization::Sp3; 2],
            EDGES,
        );
        let existing_prefix = UffTypingDiagnostic {
            atom_id: Some(AtomId::new(99)),
            kind: UffTypingDiagnosticKind::Error,
            message_prefix: UNRECOGNIZED_ATOM_TYPE_MESSAGE,
        };
        diagnostics = vec![existing_prefix];
        let value = t16_query_value(
            &central_hg_pair,
            [0, 1, 2, 3],
            &[0, 2, 2, 0],
            &[false; 4],
            &params,
            &mut diagnostics,
        )
        .expect("both Hg labels have source parameter rows")
        .expect("both central lookups succeed");
        // Hg1+2 has U1=0.1; equation17 on a single bond yields 0.5.
        assert!((value - 0.5).abs() <= 1.0e-12);
        assert_eq!(
            diagnostics,
            vec![existing_prefix, t16_warning(1), t16_warning(2)]
        );

        diagnostics.clear();
        assert!(matches!(
            t16_query_value(
                &central_hg_pair,
                [0, 1, 2, 3],
                &[0; 3],
                &[false; 4],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: 4,
                actual: 3,
            })
        ));
        assert!(matches!(
            t16_query_value(
                &central_hg_pair,
                [0, 1, 2, 3],
                &[0; 4],
                &[false; 3],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::ConjugatedBondPresence,
                expected: 4,
                actual: 3,
            })
        ));

        let index_topology = t16_path(
            (6, Hybridization::Sp3),
            (6, Hybridization::Sp3),
            [Hybridization::Sp3; 2],
            EDGES,
        );
        for (indices, expected_index) in [
            ([99, 1, 2, 3], 99),
            ([0, 99, 2, 3], 99),
            ([0, 1, 99, 3], 99),
            ([0, 1, 2, 99], 99),
        ] {
            assert_eq!(
                t16_query_value(
                    &index_topology,
                    indices,
                    &[0; 4],
                    &[false; 4],
                    &params,
                    &mut diagnostics,
                ),
                Err(UffTypingError::AtomIndexOutOfBounds {
                    index: expected_index,
                    atom_count: 4,
                })
            );
        }
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t17_central_element_degree_and_hybridization_matrix() {
        for atomic_number in 0_u8..=118 {
            for degree in [2, 3, 4] {
                for hybridization in [Hybridization::Sp2, Hybridization::Sp3] {
                    let (topology, indices) =
                        inversion_topology(atomic_number, hybridization, degree, 0);
                    let actual = get_uff_inversion_params(
                        &topology, indices[0], indices[1], indices[2], indices[3],
                    )
                    .expect("all fixture indices are valid");
                    let source_allows_hybridization =
                        !matches!(atomic_number, 6 | 7 | 8) || hybridization == Hybridization::Sp2;
                    let expected = if degree == 3 && source_allows_hybridization {
                        expected_unpromoted_inversion_k(atomic_number)
                    } else {
                        None
                    };
                    match (expected, actual) {
                        (Some(expected), Some(actual)) => assert!(
                            (actual.k - expected).abs() <= 1.0e-12,
                            "atomic number {atomic_number}, degree {degree}, hybridization {hybridization:?}: K {}, expected {expected}",
                            actual.k
                        ),
                        (None, None) => {}
                        (expected, actual) => panic!(
                            "atomic number {atomic_number}, degree {degree}, hybridization {hybridization:?}: actual {actual:?}, expected {expected:?}"
                        ),
                    }
                }
            }
        }
    }

    #[test]
    fn cf3d_typ_t17_terminal_sp2_oxygen_mask_promotes_only_carbon() {
        for atomic_number in [6, 7, 8, 15, 33, 51, 83] {
            let hybridization = if matches!(atomic_number, 6 | 7 | 8) {
                Hybridization::Sp2
            } else {
                Hybridization::Sp3
            };
            for oxygen_mask in 0_u8..8 {
                let (topology, indices) =
                    inversion_topology(atomic_number, hybridization, 3, oxygen_mask);
                let actual = get_uff_inversion_params(
                    &topology, indices[0], indices[1], indices[2], indices[3],
                )
                .expect("all fixture indices are valid")
                .expect("this center and degree are accepted by the source");
                let expected = if atomic_number == 6 && oxygen_mask != 0 {
                    50.0 / 3.0
                } else {
                    expected_unpromoted_inversion_k(atomic_number)
                        .expect("each fixture center is source-supported")
                };
                assert!(
                    (actual.k - expected).abs() <= 1.0e-12,
                    "atomic number {atomic_number}, oxygen mask {oxygen_mask:#05b}: K {}, expected {expected}",
                    actual.k
                );
            }
        }
    }

    #[test]
    fn cf3d_typ_t17_missing_bonds_and_repeated_indices_return_source_false() {
        let atoms = || {
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(6, 1, Hybridization::Sp2, false, None),
                label_atom(6, 2, Hybridization::Sp3, false, None),
                label_atom(6, 3, Hybridization::Sp3, false, None),
            ]
        };
        let single = BondOrder::Single;
        let missing_each_required_bond = [
            vec![(1, 2, single), (1, 3, single)],
            vec![(0, 1, single), (1, 3, single)],
            vec![(0, 1, single), (1, 2, single)],
        ];
        for edges in missing_each_required_bond {
            let topology = topology_with_bonds(atoms(), &edges);
            assert_eq!(get_uff_inversion_params(&topology, 0, 1, 2, 3), Ok(None));
        }

        let first_edge_missing = topology_with_bonds(atoms(), &[(1, 2, single), (1, 3, single)]);
        assert_eq!(
            get_uff_inversion_params(&first_edge_missing, 0, 1, 2, 99),
            Ok(None),
            "the missing first edge short-circuits before invalid idx4"
        );
        let second_edge_missing = topology_with_bonds(atoms(), &[(0, 1, single), (1, 3, single)]);
        assert_eq!(
            get_uff_inversion_params(&second_edge_missing, 0, 1, 2, 99),
            Ok(None),
            "the missing second edge short-circuits before invalid idx4"
        );

        let (repeated_index_topology, repeated_indices) =
            inversion_topology(6, Hybridization::Sp2, 2, 0);
        assert_eq!(repeated_indices, [0, 1, 2, 0]);
        assert_eq!(
            get_uff_inversion_params(
                &repeated_index_topology,
                repeated_indices[0],
                repeated_indices[1],
                repeated_indices[2],
                repeated_indices[3],
            ),
            Ok(None),
            "repeated terminal indices leave central degree two"
        );
    }

    #[test]
    fn cf3d_typ_t17_invalid_indices_follow_source_edge_access_order() {
        let topology = topology_with_bonds(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(6, 1, Hybridization::Sp2, false, None),
                label_atom(6, 2, Hybridization::Sp3, false, None),
                label_atom(6, 3, Hybridization::Sp3, false, None),
            ],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
            ],
        );
        for (indices, invalid_index) in [
            ([99, 1, 2, 3], 99),
            ([0, 99, 2, 3], 99),
            ([0, 1, 99, 3], 99),
            ([0, 1, 2, 99], 99),
        ] {
            assert_eq!(
                get_uff_inversion_params(&topology, indices[0], indices[1], indices[2], indices[3],),
                Err(UffTypingError::AtomIndexOutOfBounds {
                    index: invalid_index,
                    atom_count: 4,
                })
            );
        }
    }

    #[test]
    fn cf3d_typ_t18_bond_state_pair_shape_and_asymmetric_parameter_values() {
        let atoms = vec![
            label_atom(6, 0, Hybridization::Sp3, false, None),
            label_atom(7, 1, Hybridization::Sp3, false, None),
        ];
        let unbonded = topology_from_atoms(atoms.clone());
        let bonded = topology_with_bond(atoms, BondOrder::Single);
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let total_valences = [4, 3];
        let conjugation = [false; 2];

        // Pinned Params.cpp C_3/N_3 rows: x1=(3.851, 3.66),
        // D1=(0.105, 0.069); expected mixed values are fixed independently.
        let pairs = [
            (0, 1, 3.754_285_551_206_780_5, 0.085_117_565_754_666_65),
            (1, 0, 3.754_285_551_206_780_5, 0.085_117_565_754_666_65),
            (0, 0, 3.851, 0.105),
            (1, 1, 3.66, 0.069),
        ];

        for topology in [&unbonded, &bonded] {
            let topology_before = topology.clone();
            for (idx1, idx2, expected_x, expected_d) in pairs {
                let mut diagnostics = Vec::new();
                let vdw = get_uff_vdw_params(
                    topology,
                    idx1,
                    idx2,
                    &total_valences,
                    &conjugation,
                    &params,
                    &mut diagnostics,
                )
                .expect("source-valid pair lookup")
                .expect("both fixed UFF rows exist");

                assert!((vdw.x_ij - expected_x).abs() <= 1.0e-12);
                assert!((vdw.d_ij - expected_d).abs() <= 1.0e-12);
                assert!(diagnostics.is_empty());
            }
            assert_eq!(*topology, topology_before);
        }
        assert_eq!(total_valences, [4, 3]);
        assert_eq!(conjugation, [false; 2]);
    }

    #[test]
    fn cf3d_typ_t18_missing_first_and_second_types_preserve_short_circuit() {
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let total_valences = [0, 4];
        let conjugation = [false; 2];
        let missing_first_atoms = || {
            vec![
                label_atom(0, 0, Hybridization::Sp3, false, None),
                label_atom(6, 1, Hybridization::Sp3, false, None),
            ]
        };
        let missing_first_unbonded = topology_from_atoms(missing_first_atoms());
        let missing_first_bonded = topology_with_bond(missing_first_atoms(), BondOrder::Single);
        let mut diagnostics = Vec::new();
        for topology in [&missing_first_unbonded, &missing_first_bonded] {
            let topology_before = topology.clone();
            assert_eq!(
                get_uff_vdw_params(
                    topology,
                    0,
                    99,
                    &total_valences,
                    &conjugation,
                    &params,
                    &mut diagnostics,
                ),
                Ok(None),
                "a missing first type short-circuits before the invalid second index"
            );
            assert!(diagnostics.is_empty());
            assert_eq!(*topology, topology_before);
        }

        let missing_second_atoms = || {
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(0, 1, Hybridization::Sp3, false, None),
            ]
        };
        let missing_second_unbonded = topology_from_atoms(missing_second_atoms());
        let missing_second_bonded = topology_with_bond(missing_second_atoms(), BondOrder::Single);
        for topology in [&missing_second_unbonded, &missing_second_bonded] {
            let topology_before = topology.clone();
            diagnostics.clear();
            assert_eq!(
                get_uff_vdw_params(
                    topology,
                    0,
                    1,
                    &total_valences,
                    &conjugation,
                    &params,
                    &mut diagnostics,
                ),
                Ok(None),
                "a missing second type returns the source false result"
            );
            assert!(diagnostics.is_empty());
            assert_eq!(*topology, topology_before);
        }

        diagnostics.clear();
        assert_eq!(
            get_uff_vdw_params(
                &missing_second_unbonded,
                99,
                0,
                &total_valences,
                &conjugation,
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::AtomIndexOutOfBounds {
                index: 99,
                atom_count: 2,
            })
        );
        assert!(diagnostics.is_empty());

        let known_types = topology_from_atoms(vec![
            label_atom(6, 0, Hybridization::Sp3, false, None),
            label_atom(7, 1, Hybridization::Sp3, false, None),
        ]);
        assert_eq!(
            get_uff_vdw_params(
                &known_types,
                0,
                99,
                &[4, 3],
                &[false; 2],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::AtomIndexOutOfBounds {
                index: 99,
                atom_count: 2,
            })
        );
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_t18_label_diagnostics_follow_pair_access_order() {
        let topology = topology_from_atoms(vec![
            label_atom(80, 0, Hybridization::Sp3, false, None),
            label_atom(80, 1, Hybridization::Sp3, false, None),
        ]);
        let topology_before = topology.clone();
        let total_valences = [2, 2];
        let conjugation = [false; 2];
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();

        let vdw = get_uff_vdw_params(
            &topology,
            0,
            1,
            &total_valences,
            &conjugation,
            &params,
            &mut diagnostics,
        )
        .expect("both Hg labels have source parameters")
        .expect("both Hg1+2 parameter rows exist");

        assert_eq!(vdw.x_ij, 2.705);
        assert_eq!(vdw.d_ij, 0.385);
        assert_eq!(
            diagnostics,
            [
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(0)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
                },
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(1)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
                },
            ]
        );
        diagnostics.clear();
        let reversed = get_uff_vdw_params(
            &topology,
            1,
            0,
            &total_valences,
            &conjugation,
            &params,
            &mut diagnostics,
        )
        .expect("reversed source-valid pair lookup")
        .expect("both reversed Hg1+2 parameter rows exist");
        assert_eq!(reversed, vdw);
        assert_eq!(
            diagnostics,
            [
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(1)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
                },
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(0)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
                },
            ]
        );
        assert_eq!(topology, topology_before);
        assert_eq!(total_valences, [2, 2]);
        assert_eq!(conjugation, [false; 2]);
    }

    #[test]
    fn cf3d_typ_t18_prepared_state_lengths_are_typed_errors() {
        let topology = topology_from_atoms(vec![
            label_atom(6, 0, Hybridization::Sp3, false, None),
            label_atom(7, 1, Hybridization::Sp3, false, None),
        ]);
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();

        assert_eq!(
            get_uff_vdw_params(
                &topology,
                0,
                1,
                &[4],
                &[false; 2],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: super::UffTypingInput::TotalValence,
                expected: 2,
                actual: 1,
            })
        );
        assert!(diagnostics.is_empty());
        assert_eq!(
            get_uff_vdw_params(
                &topology,
                0,
                1,
                &[4, 3],
                &[false],
                &params,
                &mut diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: super::UffTypingInput::ConjugatedBondPresence,
                expected: 2,
                actual: 1,
            })
        );
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_typ_integration_local_queries_continue_after_missing_row() {
        let topology = topology_with_bonds(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(6, 1, Hybridization::Sp3, false, None),
                label_atom(7, 2, Hybridization::Sp3, false, None),
                label_atom(6, 3, Hybridization::Sp3, false, None),
                label_atom(6, 4, Hybridization::Sp2, false, None),
                label_atom(80, 5, Hybridization::Sp3, false, None),
                label_atom(81, 6, Hybridization::Sp2, false, None),
                label_atom(0, 7, Hybridization::Sp3, false, None),
            ],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
                (4, 2, BondOrder::Single),
                (4, 3, BondOrder::Single),
            ],
        );
        let topology_before = topology.clone();
        let total_valences = [4, 4, 3, 4, 3, 2, 3, 0];
        let conjugation = [false; 8];
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();

        let (slots, found_all) = get_atom_types(
            &topology,
            &total_valences,
            &conjugation,
            params.as_ref(),
            &mut diagnostics,
        )
        .expect("aligned detached typing state");
        assert!(!found_all, "the final dummy row has no default UFF entry");
        assert_eq!(slots.len(), 8);
        for (index, label) in ["C_3", "C_3", "N_3", "C_3", "C_2", "Hg1+2", "Tl3+3"]
            .into_iter()
            .enumerate()
        {
            assert!(std::ptr::eq(
                slots[index].expect("fixed known label has a parameter row"),
                params.get(label).expect("pinned source parameter row"),
            ));
        }
        assert!(slots[7].is_none());
        assert_eq!(
            diagnostics,
            [
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(5)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
                },
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(6)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP3_HYBRIDIZATION_WARNING_MESSAGE,
                },
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(7)),
                    kind: UffTypingDiagnosticKind::Error,
                    message_prefix: UNRECOGNIZED_ATOM_TYPE_MESSAGE,
                },
            ]
        );

        // Fixed from the pinned RDKit C_3/N_3 pair and testUFFForceField.cpp:
        // C_3-N_3 has r0=1.451071 and kb=1057.27.
        let bond = get_uff_bond_stretch_params(
            &topology,
            1,
            2,
            &total_valences,
            &conjugation,
            params.as_ref(),
            &mut diagnostics,
        )
        .expect("known local edge query")
        .expect("both local endpoint types exist despite the final missing row");
        assert!((bond.r0 - 1.451_071).abs() < 1.0e-5);
        assert!((bond.kb - 1_057.27).abs() < 1.0e-2);

        let angle = get_uff_angle_bend_params(
            &topology,
            0,
            1,
            2,
            &total_valences,
            &conjugation,
            params.as_ref(),
            &mut diagnostics,
        )
        .expect("known local angle query")
        .expect("all three local rows exist");
        // Fixed from pinned testUFFHelpers.cpp:1054-1057 (C_3-C_3-N_3).
        assert_eq!((angle.ka * 1000.0).round() as i64, 303_297);
        assert_eq!((angle.theta0 * 1000.0).round() as i64, 109_470);

        let torsion = get_uff_torsion_params(
            &topology,
            0,
            1,
            2,
            3,
            &total_valences,
            &conjugation,
            params.as_ref(),
            &mut diagnostics,
        )
        .expect("known local torsion path")
        .expect("both central rows exist");
        // Pinned Params.cpp C_3/N_3 V1=(2.119, 0.45); source SP3/SP3 uses sqrt(V1a*V1b).
        assert!((torsion.v - 0.976_498_847_925_587_7).abs() <= 1.0e-12);

        let inversion = get_uff_inversion_params(&topology, 0, 4, 2, 3)
            .expect("valid central carbon star")
            .expect("SP2 degree-three carbon is source-supported");
        // Source C inversion coefficients use K=6, then the helper divides by 3.
        assert_eq!(inversion.k, 2.0);

        let vdw = get_uff_vdw_params(
            &topology,
            5,
            6,
            &total_valences,
            &conjugation,
            params.as_ref(),
            &mut diagnostics,
        )
        .expect("isolated atoms need no bond for a vdw query")
        .expect("Hg1+2 and Tl3+3 are known despite another missing type");
        // Pinned Params.cpp Hg1+2/Tl3+3 rows set x1=(2.705,4.347),
        // D1=(0.385,0.68), and Nonbonded.cpp mixes each pair by sqrt(product).
        assert!((vdw.x_ij - 3.429_086_613_079_349).abs() <= 1.0e-12);
        assert!((vdw.d_ij - 0.511_663_952_218_641_2).abs() <= 1.0e-12);
        assert_eq!(
            diagnostics,
            [
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(5)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
                },
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(6)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP3_HYBRIDIZATION_WARNING_MESSAGE,
                },
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(7)),
                    kind: UffTypingDiagnosticKind::Error,
                    message_prefix: UNRECOGNIZED_ATOM_TYPE_MESSAGE,
                },
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(5)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP_HYBRIDIZATION_WARNING_MESSAGE,
                },
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(6)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: FORCED_SP3_HYBRIDIZATION_WARNING_MESSAGE,
                },
            ]
        );
        assert_eq!(topology, topology_before);
        assert_eq!(total_valences, [4, 4, 3, 4, 3, 2, 3, 0]);
        assert_eq!(conjugation, [false; 8]);
    }

    #[test]
    fn cf3d_typ_integration_prepared_valence_and_conjugation_drive_labels() {
        let topology = topology_from_atoms(vec![
            label_atom(6, 0, Hybridization::Sp2, false, None),
            label_atom(15, 1, Hybridization::Sp3, false, None),
        ]);
        let topology_before = topology.clone();
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut diagnostics = Vec::new();

        let (base_slots, base_found_all) = get_atom_types(
            &topology,
            &[2, 3],
            &[false, false],
            params.as_ref(),
            &mut diagnostics,
        )
        .expect("fixed prepared state");
        assert!(base_found_all);
        assert!(std::ptr::eq(
            base_slots[0].expect("C_2 row"),
            params.get("C_2").expect("pinned C_2 row"),
        ));
        assert!(std::ptr::eq(
            base_slots[1].expect("P_3+3 row"),
            params.get("P_3+3").expect("pinned P_3+3 row"),
        ));
        assert!(diagnostics.is_empty());

        let (conjugated_slots, conjugated_found_all) = get_atom_types(
            &topology,
            &[2, 3],
            &[true, false],
            params.as_ref(),
            &mut diagnostics,
        )
        .expect("only supplied conjugation changed");
        assert!(conjugated_found_all);
        assert!(std::ptr::eq(
            conjugated_slots[0].expect("C_R row"),
            params.get("C_R").expect("pinned C_R row"),
        ));
        assert!(std::ptr::eq(
            conjugated_slots[1].expect("unchanged P_3+3 row"),
            params.get("P_3+3").expect("pinned P_3+3 row"),
        ));
        assert!(diagnostics.is_empty());

        let (valence_changed_slots, valence_changed_found_all) = get_atom_types(
            &topology,
            &[2, 5],
            &[false, false],
            params.as_ref(),
            &mut diagnostics,
        )
        .expect("only supplied total valence changed");
        assert!(valence_changed_found_all);
        assert!(std::ptr::eq(
            valence_changed_slots[0].expect("unchanged C_2 row"),
            params.get("C_2").expect("pinned C_2 row"),
        ));
        assert!(std::ptr::eq(
            valence_changed_slots[1].expect("P_3+5 row"),
            params.get("P_3+5").expect("pinned P_3+5 row"),
        ));
        assert!(diagnostics.is_empty());
        assert_eq!(topology, topology_before);
    }

    #[test]
    fn cf3d_bld_integration_prepares_queries_and_assembles_custom_chain_terms() {
        let topology = cf3d_bld_integration_topology(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(6, 1, Hybridization::Sp3, false, None),
                label_atom(7, 2, Hybridization::Sp3, false, None),
                label_atom(6, 3, Hybridization::Sp3, false, None),
            ],
            &[
                (0, 1, BondOrder::Single, false),
                (1, 2, BondOrder::Single, false),
                (2, 3, BondOrder::Single, false),
            ],
        );
        let assignment = cf3d_bld_integration_assignment(&[4, 4, 3, 4], &[0; 4]);
        let table = ParamCollection::get_params("").expect("pinned default UFF parameter table");
        let mut diagnostics = Vec::new();
        let (prepared, conjugated, slots, found_all, needs_hydrogens) =
            cf3d_bld_integration_prepare(&topology, &assignment, &table, &mut diagnostics);
        assert_eq!(prepared.total_valences, [4, 4, 3, 4]);
        assert_eq!(conjugated, [false; 4]);
        assert!(found_all);
        assert!(!needs_hydrogens);
        assert!(diagnostics.is_empty());
        for (index, label) in ["C_3", "C_3", "N_3", "C_3"].into_iter().enumerate() {
            assert!(std::ptr::eq(
                slots[index].expect("fixed chain atom type exists"),
                table
                    .get(label)
                    .expect("pinned source parameter row exists"),
            ));
        }

        // These fixed local values are pinned by AtomTyper.cpp, Params.cpp,
        // Bond.cpp and testUFFHelpers.cpp. The rows remain borrowed from the
        // one canonical cached table used by the subsequent real builder calls.
        let vdw = get_uff_vdw_params(
            &topology,
            0,
            1,
            &prepared.total_valences,
            &conjugated,
            &table,
            &mut diagnostics,
        )
        .expect("fixed C_3 pair query succeeds")
        .expect("both fixed C_3 rows exist");
        assert!((vdw.x_ij - 3.851).abs() < 1.0e-12);
        assert!((vdw.d_ij - 0.105).abs() < 1.0e-12);

        let bond = get_uff_bond_stretch_params(
            &topology,
            1,
            2,
            &prepared.total_valences,
            &conjugated,
            &table,
            &mut diagnostics,
        )
        .expect("fixed C_3-N_3 edge query succeeds")
        .expect("both fixed endpoint rows exist");
        assert!((bond.r0 - 1.451_071).abs() < 1.0e-5);
        assert!((bond.kb - 1_057.27).abs() < 1.0e-2);

        let angle = get_uff_angle_bend_params(
            &topology,
            0,
            1,
            2,
            &prepared.total_valences,
            &conjugated,
            &table,
            &mut diagnostics,
        )
        .expect("fixed C_3-C_3-N_3 angle query succeeds")
        .expect("all fixed angle rows exist");
        assert_eq!((angle.ka * 1000.0).round() as i64, 303_297);
        assert_eq!((angle.theta0 * 1000.0).round() as i64, 109_470);

        let torsion = get_uff_torsion_params(
            &topology,
            0,
            1,
            2,
            3,
            &prepared.total_valences,
            &conjugated,
            &table,
            &mut diagnostics,
        )
        .expect("fixed C_3-N_3 torsion query succeeds")
        .expect("both fixed central rows exist");
        assert!((torsion.v - 0.976_498_847_925_587_7).abs() <= 1.0e-12);
        assert!(diagnostics.is_empty());

        let rings = cosmolkit_core::fast_find_rings(&topology)
            .expect("fixed acyclic source topology has ring state");
        let points = [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
        ];
        let coordinates = points
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let mut actual_points = points;
        let mut actual = ForceField::new(3);
        cf3d_bld_integration_attach(&mut actual, &mut actual_points);
        builder::add_bonds(&topology, &slots, &mut actual)
            .expect("source-order chain bonds append");
        builder::add_angles(&topology, &slots, &rings, &mut actual)
            .expect("source-order chain angles append");
        builder::add_torsions(
            &topology,
            &rings,
            &assignment,
            &slots,
            "[D2]~[D2]",
            &mut actual,
        )
        .expect("the custom source SMARTS appends the central chain torsion");
        builder::add_inversions(&topology, &slots, &mut actual)
            .expect("acyclic chain has no inversion centers");

        let mut expected_points = points;
        let mut expected = ForceField::new(3);
        cf3d_bld_integration_attach(&mut expected, &mut expected_points);
        cf3d_bld_integration_add_expected_bond(&mut expected, 0, 1, 1.0, &slots);
        cf3d_bld_integration_add_expected_bond(&mut expected, 1, 2, 1.0, &slots);
        cf3d_bld_integration_add_expected_bond(&mut expected, 2, 3, 1.0, &slots);
        cf3d_bld_integration_add_expected_angle(&mut expected, 0, 1, 2, 1.0, 1.0, 0, &slots);
        cf3d_bld_integration_add_expected_angle(&mut expected, 1, 2, 3, 1.0, 1.0, 0, &slots);
        cf3d_bld_integration_add_expected_torsion(
            &mut expected,
            &topology,
            [0, 1, 2, 3],
            1.0,
            false,
            &slots,
        );
        cf3d_bld_integration_assert_fields(actual, expected, &coordinates);
    }

    #[test]
    fn cf3d_bld_integration_ring_angles_and_torsion_triangle_exclusion() {
        let topology = cf3d_bld_integration_topology(
            vec![
                label_atom(6, 0, Hybridization::Sp2, false, None),
                label_atom(6, 1, Hybridization::Sp2, false, None),
                label_atom(6, 2, Hybridization::Sp2, false, None),
            ],
            &[
                (0, 1, BondOrder::Single, false),
                (1, 2, BondOrder::Single, false),
                (2, 0, BondOrder::Single, false),
            ],
        );
        let assignment = cf3d_bld_integration_assignment(&[2; 3], &[0; 3]);
        let table = ParamCollection::get_params("").expect("pinned default UFF parameter table");
        let mut diagnostics = Vec::new();
        let (prepared, conjugated, slots, found_all, needs_hydrogens) =
            cf3d_bld_integration_prepare(&topology, &assignment, &table, &mut diagnostics);
        assert_eq!(prepared.total_valences, [2; 3]);
        assert_eq!(conjugated, [false; 3]);
        assert!(found_all);
        assert!(!needs_hydrogens);
        assert!(diagnostics.is_empty());
        for slot in &slots {
            assert!(std::ptr::eq(
                slot.expect("fixed ring carbon has a UFF row"),
                table.get("C_2").expect("pinned source C_2 row"),
            ));
        }

        let rings = cosmolkit_core::fast_find_rings(&topology)
            .expect("fixed source triangle has ring state");
        for atom in 0..3 {
            assert!(rings.is_atom_in_ring_of_size(AtomId::new(atom), 3));
        }
        let side = 3.0_f64.sqrt() / 2.0;
        let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, side, 0.0]];
        let coordinates = points
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let mut actual_points = points;
        let mut actual = ForceField::new(3);
        cf3d_bld_integration_attach(&mut actual, &mut actual_points);
        builder::add_bonds(&topology, &slots, &mut actual).expect("source-order ring bonds append");
        builder::add_angles(&topology, &slots, &rings, &mut actual)
            .expect("source ring-angle substitutions append");
        builder::add_torsions(
            &topology,
            &rings,
            &assignment,
            &slots,
            "[!$(*#*)&!D1]~[!$(*#*)&!D1]",
            &mut actual,
        )
        .expect("the source triangle produces no four-atom torsions");
        builder::add_inversions(&topology, &slots, &mut actual)
            .expect("degree-two ring atoms are not inversion centers");

        let mut expected_points = points;
        let mut expected = ForceField::new(3);
        cf3d_bld_integration_attach(&mut expected, &mut expected_points);
        cf3d_bld_integration_add_expected_bond(&mut expected, 0, 1, 1.0, &slots);
        cf3d_bld_integration_add_expected_bond(&mut expected, 1, 2, 1.0, &slots);
        cf3d_bld_integration_add_expected_bond(&mut expected, 2, 0, 1.0, &slots);
        cf3d_bld_integration_add_expected_angle(&mut expected, 1, 0, 2, 1.0, 1.0, 35, &slots);
        cf3d_bld_integration_add_expected_angle(&mut expected, 0, 1, 2, 1.0, 1.0, 35, &slots);
        cf3d_bld_integration_add_expected_angle(&mut expected, 1, 2, 0, 1.0, 1.0, 35, &slots);
        cf3d_bld_integration_assert_fields(actual, expected, &coordinates);
    }

    #[test]
    fn cf3d_bld_integration_conjugated_carbonyl_and_inversion_kernel_terms() {
        let topology = cf3d_bld_integration_topology(
            vec![
                label_atom(6, 0, Hybridization::Sp2, false, None),
                label_atom(8, 1, Hybridization::Sp2, false, None),
                label_atom(6, 2, Hybridization::Sp3, false, None),
                label_atom(6, 3, Hybridization::Sp3, false, None),
            ],
            &[
                (0, 1, BondOrder::Double, true),
                (0, 2, BondOrder::Single, false),
                (0, 3, BondOrder::Single, false),
            ],
        );
        let assignment = cf3d_bld_integration_assignment(&[4, 2, 1, 1], &[0; 4]);
        let table = ParamCollection::get_params("").expect("pinned default UFF parameter table");
        let mut diagnostics = Vec::new();
        let (prepared, conjugated, slots, found_all, needs_hydrogens) =
            cf3d_bld_integration_prepare(&topology, &assignment, &table, &mut diagnostics);
        assert_eq!(prepared.total_valences, [4, 2, 1, 1]);
        assert_eq!(conjugated, [true, true, false, false]);
        assert!(found_all);
        assert!(!needs_hydrogens);
        assert!(std::ptr::eq(
            slots[0].expect("conjugated carbonyl carbon type exists"),
            table.get("C_R").expect("pinned source C_R row"),
        ));
        assert!(std::ptr::eq(
            slots[1].expect("conjugated carbonyl oxygen type exists"),
            table.get("O_R").expect("pinned source O_R row"),
        ));
        assert!(diagnostics.is_empty());

        let carbonyl_bond = get_uff_bond_stretch_params(
            &topology,
            0,
            1,
            &prepared.total_valences,
            &conjugated,
            &table,
            &mut diagnostics,
        )
        .expect("fixed C_R=O_R query succeeds")
        .expect("both source parameter rows exist");
        assert!(carbonyl_bond.r0 > 0.0 && carbonyl_bond.kb > 0.0);

        let carbonyl_angle = get_uff_angle_bend_params(
            &topology,
            1,
            0,
            2,
            &prepared.total_valences,
            &conjugated,
            &table,
            &mut diagnostics,
        )
        .expect("fixed O_R-C_R-C_3 angle query succeeds")
        .expect("all source angle parameter rows exist");
        assert!(carbonyl_angle.ka > 0.0);
        assert_eq!((carbonyl_angle.theta0 * 1000.0).round() as i64, 120_000);

        let inversion = get_uff_inversion_params(&topology, 1, 0, 2, 3)
            .expect("fixed carbonyl central star query succeeds")
            .expect("SP2 degree-three carbon is source-supported");
        assert_eq!(inversion.k, 50.0 / 3.0);
        assert!(diagnostics.is_empty());

        let rings = cosmolkit_core::fast_find_rings(&topology)
            .expect("fixed carbonyl topology has ring state");
        let points = [
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let coordinates = points
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let mut actual_points = points;
        let mut actual = ForceField::new(3);
        cf3d_bld_integration_attach(&mut actual, &mut actual_points);
        builder::add_bonds(&topology, &slots, &mut actual)
            .expect("source-order carbonyl bonds append");
        builder::add_angles(&topology, &slots, &rings, &mut actual)
            .expect("source-order carbonyl angles append");
        builder::add_angle_special_cases(&topology, &slots, &mut actual)
            .expect("carbonyl center does not dispatch the TBP special path");
        builder::add_inversions(&topology, &slots, &mut actual)
            .expect("source carbonyl inversion permutations append");

        let mut expected_points = points;
        let mut expected = ForceField::new(3);
        cf3d_bld_integration_attach(&mut expected, &mut expected_points);
        cf3d_bld_integration_add_expected_bond(&mut expected, 0, 1, 2.0, &slots);
        cf3d_bld_integration_add_expected_bond(&mut expected, 0, 2, 1.0, &slots);
        cf3d_bld_integration_add_expected_bond(&mut expected, 0, 3, 1.0, &slots);
        cf3d_bld_integration_add_expected_angle(&mut expected, 1, 0, 2, 2.0, 1.0, 3, &slots);
        cf3d_bld_integration_add_expected_angle(&mut expected, 1, 0, 3, 2.0, 1.0, 3, &slots);
        cf3d_bld_integration_add_expected_angle(&mut expected, 2, 0, 3, 1.0, 1.0, 3, &slots);
        cf3d_bld_integration_add_expected_inversion(&mut expected, [1, 0, 2, 3], 6, true);
        cf3d_bld_integration_add_expected_inversion(&mut expected, [1, 0, 3, 2], 6, true);
        cf3d_bld_integration_add_expected_inversion(&mut expected, [2, 0, 3, 1], 6, true);
        cf3d_bld_integration_assert_fields(actual, expected, &coordinates);
    }

    #[test]
    fn cf3d_bld_integration_tbp_dispatch_keeps_the_source_angle_term_order() {
        let topology = cf3d_bld_integration_topology(
            vec![
                label_atom(15, 0, Hybridization::Sp3d, false, None),
                label_atom(6, 1, Hybridization::Sp3, false, None),
                label_atom(6, 2, Hybridization::Sp3, false, None),
                label_atom(6, 3, Hybridization::Sp3, false, None),
                label_atom(6, 4, Hybridization::Sp3, false, None),
                label_atom(6, 5, Hybridization::Sp3, false, None),
            ],
            &[
                (0, 1, BondOrder::Single, false),
                (0, 2, BondOrder::Single, false),
                (0, 3, BondOrder::Single, false),
                (0, 4, BondOrder::Single, false),
                (0, 5, BondOrder::Single, false),
            ],
        );
        let assignment = cf3d_bld_integration_assignment(&[5, 1, 1, 1, 1, 1], &[0; 6]);
        let table = ParamCollection::get_params("").expect("pinned default UFF parameter table");
        let mut diagnostics = Vec::new();
        let (prepared, conjugated, slots, found_all, needs_hydrogens) =
            cf3d_bld_integration_prepare(&topology, &assignment, &table, &mut diagnostics);
        assert_eq!(prepared.total_valences, [5, 1, 1, 1, 1, 1]);
        assert_eq!(conjugated, [false; 6]);
        assert!(found_all);
        assert!(!needs_hydrogens);
        assert!(std::ptr::eq(
            slots[0].expect("fixed phosphorus type exists"),
            table.get("P_3+5").expect("pinned source P_3+5 row"),
        ));
        for slot in &slots[1..] {
            assert!(std::ptr::eq(
                slot.expect("fixed terminal carbon type exists"),
                table.get("C_3").expect("pinned source C_3 row"),
            ));
        }
        assert_eq!(
            diagnostics,
            [UffTypingDiagnostic {
                atom_id: Some(AtomId::new(0)),
                kind: UffTypingDiagnosticKind::Warning,
                message_prefix: FORCED_SP3_HYBRIDIZATION_WARNING_MESSAGE,
            }]
        );

        let rings =
            cosmolkit_core::fast_find_rings(&topology).expect("fixed TBP star has ring state");
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
        ];
        let coordinates = points
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let mut actual_points = points;
        let mut actual = ForceField::new(3);
        cf3d_bld_integration_attach(&mut actual, &mut actual_points);
        builder::add_bonds(&topology, &slots, &mut actual).expect("source-order TBP bonds append");
        builder::add_angles(&topology, &slots, &rings, &mut actual)
            .expect("ordinary angle pass excludes the SP3D degree-five center");
        builder::add_angle_special_cases(&topology, &slots, &mut actual)
            .expect("source TBP selector appends its axial/equatorial terms");
        builder::add_torsions(
            &topology,
            &rings,
            &assignment,
            &slots,
            "[!$(*#*)&!D1]~[!$(*#*)&!D1]",
            &mut actual,
        )
        .expect("terminal-only TBP star has no torsion matches");
        builder::add_inversions(&topology, &slots, &mut actual)
            .expect("SP3D degree-five center is not an inversion center");

        let mut expected_points = points;
        let mut expected = ForceField::new(3);
        cf3d_bld_integration_attach(&mut expected, &mut expected_points);
        for endpoint in 1..=5 {
            cf3d_bld_integration_add_expected_bond(&mut expected, 0, endpoint, 1.0, &slots);
        }
        for (first, center, last, order) in [
            (1, 0, 2, 2),
            (3, 0, 4, 3),
            (3, 0, 5, 3),
            (4, 0, 5, 3),
            (1, 0, 3, 0),
            (1, 0, 4, 0),
            (1, 0, 5, 0),
            (2, 0, 3, 0),
            (2, 0, 4, 0),
            (2, 0, 5, 0),
        ] {
            cf3d_bld_integration_add_expected_angle(
                &mut expected,
                first,
                center,
                last,
                1.0,
                1.0,
                order,
                &slots,
            );
        }
        cf3d_bld_integration_assert_fields(actual, expected, &coordinates);
    }

    #[test]
    fn cf3d_bld_integration_minimize_gathers_borrowed_coordinate_rows() {
        let topology = cf3d_bld_integration_topology(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(6, 1, Hybridization::Sp3, false, None),
            ],
            &[(0, 1, BondOrder::Single, false)],
        );
        let assignment = cf3d_bld_integration_assignment(&[4, 4], &[0, 0]);
        let table = ParamCollection::get_params("").expect("pinned default UFF parameter table");
        let mut diagnostics = Vec::new();
        let (prepared, conjugated, slots, found_all, needs_hydrogens) =
            cf3d_bld_integration_prepare(&topology, &assignment, &table, &mut diagnostics);
        assert_eq!(prepared.total_valences, [4, 4]);
        assert_eq!(conjugated, [false; 2]);
        assert!(found_all);
        assert!(!needs_hydrogens);
        assert!(diagnostics.is_empty());
        let source_bond = get_uff_bond_stretch_params(
            &topology,
            0,
            1,
            &prepared.total_valences,
            &conjugated,
            &table,
            &mut diagnostics,
        )
        .expect("fixed C_3-C_3 query succeeds")
        .expect("both fixed carbon parameter rows exist");
        assert!((source_bond.r0 - 1.514).abs() <= 1.0e-12);

        let mut rows = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut field = ForceField::new(3);
        cf3d_bld_integration_attach(&mut field, &mut rows);
        builder::add_bonds(&topology, &slots, &mut field)
            .expect("source C_3-C_3 bond contribution appends");
        field
            .initialize()
            .expect("fixed two-point field initializes");
        let initial_coordinates = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let initial_energy = cf3d_bld_b05_calc_energy(&mut field, &initial_coordinates)
            .expect("fixed initial bond energy evaluates");

        let status = cf3d_bld_integration_minimize(&mut field, 200, 1.0e-4, 1.0e-6)
            .expect("successful BFGS status passes through the borrowed kernel context");
        assert_eq!(status, 0);
        let minimized_coordinates = field
            .positions()
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let [first, second] = field.positions() else {
            panic!("fixed two atom source field retains two coordinate rows")
        };
        let distance = (0..3)
            .map(|axis| {
                let delta = second[axis] - first[axis];
                delta * delta
            })
            .sum::<f64>()
            .sqrt();
        assert!((distance - source_bond.r0).abs() < 1.0e-3);
        assert!((distance - 2.0).abs() > 1.0e-2);
        let minimized_energy = cf3d_bld_b05_calc_energy(&mut field, &minimized_coordinates)
            .expect("gathered minimized rows evaluate through the same kernel");
        assert!(minimized_energy < initial_energy);
        let mut gradient = vec![0.0; minimized_coordinates.len()];
        cf3d_bld_b05_calc_grad(&mut field, &minimized_coordinates, &mut gradient)
            .expect("gathered minimized gradient evaluates");
        assert!(
            gradient
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt()
                < 1.0e-3
        );
    }

    #[test]
    fn cf3d_bld_integration_keeps_preparation_builder_and_query_errors_typed() {
        let topology = cf3d_bld_integration_topology(
            vec![
                label_atom(6, 0, Hybridization::Sp3, false, None),
                label_atom(6, 1, Hybridization::Sp3, false, None),
                label_atom(7, 2, Hybridization::Sp3, false, None),
                label_atom(6, 3, Hybridization::Sp3, false, None),
            ],
            &[
                (0, 1, BondOrder::Single, false),
                (1, 2, BondOrder::Single, false),
                (2, 3, BondOrder::Single, false),
            ],
        );
        let bad_assignment = cf3d_bld_integration_assignment(&[4, 4, 3], &[0, 0, 0, 0]);
        assert!(matches!(
            builder::prepare_typing_valence(&topology, &bad_assignment),
            Err(UffBuilderError::ValenceAssignmentLengthMismatch {
                field: builder::PreparedValenceField::Explicit,
                expected: 4,
                actual: 3,
            })
        ));

        let assignment = cf3d_bld_integration_assignment(&[4, 4, 3, 4], &[0; 4]);
        let table = ParamCollection::get_params("").expect("pinned default UFF parameter table");
        let mut diagnostics = Vec::new();
        let (prepared, _, slots, found_all, _) =
            cf3d_bld_integration_prepare(&topology, &assignment, &table, &mut diagnostics);
        assert!(found_all);
        assert!(matches!(
            get_atom_types(
                &topology,
                &prepared.total_valences,
                &[],
                &table,
                &mut diagnostics,
            ),
            Err(UffTypingError::PreparedStateLength {
                input: UffTypingInput::ConjugatedBondPresence,
                expected: 4,
                actual: 0,
            })
        ));

        let mut short_rows = [[0.0; 3]];
        let mut short_field = ForceField::new(3);
        cf3d_bld_integration_attach(&mut short_field, &mut short_rows);
        assert!(matches!(
            builder::add_bonds(&topology, &slots, &mut short_field),
            Err(UffBuilderError::ForceFieldKernel(
                ForceFieldKernelError::BondIndexOutOfRange { .. }
            ))
        ));

        let rings =
            cosmolkit_core::fast_find_rings(&topology).expect("fixed source chain has ring state");
        let mut complete_rows = [[0.0; 3]; 4];
        let mut complete_field = ForceField::new(3);
        cf3d_bld_integration_attach(&mut complete_field, &mut complete_rows);
        assert_eq!(
            builder::add_torsions(
                &topology,
                &rings,
                &assignment,
                &slots[..1],
                "[",
                &mut complete_field,
            ),
            Err(UffBuilderError::ParamsLengthMismatch {
                atoms: 4,
                params: 1,
            })
        );
        assert!(
            builder::add_torsions(
                &topology,
                &rings,
                &assignment,
                &slots,
                "[",
                &mut complete_field,
            )
            .is_err()
        );
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    use cosmolkit_model::{Atom, AtomId, AtomSpec, PropertyValue};
    use cosmolkit_model::{Element, Hybridization};
    // FROZEN UINT CONDITION: TEXT_CONSUMER_forcefields/UFFdummyLabel_0
    #[test]
    fn uint_cell_text_consumer_forcefields_uffdummylabel_0_atom_typer() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(0).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(0_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("0_").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(6).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(0_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("C_").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_forcefields/UFFdummyLabel_1
    #[test]
    fn uint_cell_text_consumer_forcefields_uffdummylabel_1_atom_typer() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(0).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(1_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("1_").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(6).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(1_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("C_").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_forcefields/UFFdummyLabel_2147483646
    #[test]
    fn uint_cell_text_consumer_forcefields_uffdummylabel_2147483646_atom_typer() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(0).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(2147483646_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("2147483646").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(6).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(2147483646_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("C_").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_forcefields/UFFdummyLabel_2147483647
    #[test]
    fn uint_cell_text_consumer_forcefields_uffdummylabel_2147483647_atom_typer() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(0).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(2147483647_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("2147483647").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(6).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(2147483647_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("C_").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_forcefields/UFFdummyLabel_2147483648
    #[test]
    fn uint_cell_text_consumer_forcefields_uffdummylabel_2147483648_atom_typer() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(0).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(2147483648_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("2147483648").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(6).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(2147483648_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("C_").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_forcefields/UFFdummyLabel_4294967295
    #[test]
    fn uint_cell_text_consumer_forcefields_uffdummylabel_4294967295_atom_typer() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(0).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(4294967295_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("4294967295").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::from_atomic_number(6).unwrap())
                .with_hybridization(Hybridization::S)
                .with_prop("dummyLabel", PropertyValue::UInt(4294967295_u32))
                .unwrap(),
        );
        let before = a.clone();
        let mut diagnostics = vec![];
        assert_eq!(
            (atom_label_prefix(&a, Hybridization::S, || false, &mut diagnostics).unwrap())
                .as_bytes(),
            ("C_").as_bytes()
        );
        assert_eq!(diagnostics, vec![]);
        assert_eq!(a, before);
    }
}

#[cfg(test)]
mod canonical_byte_label_regressions {
    use super::*;
    use cosmolkit_model::Element;
    use cosmolkit_model::{AtomSpec, PropertyText, PropertyValue};
    #[test]
    fn opaque_dummy_label_bytes_reach_exact_parameter_key_lookup() {
        let parameters = ParamCollection::get_params("").unwrap();
        for (input, expected) in [
            (b"\xff".as_slice(), b"\xff_".as_slice()),
            (b"\xff\0".as_slice(), b"\xff\0".as_slice()),
            (b"C_3\0".as_slice(), b"C_3\0".as_slice()),
            (b"C_3".as_slice(), b"C_3".as_slice()),
        ] {
            let atom = Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::from_atomic_number(0).unwrap())
                    .with_prop(
                        "dummyLabel",
                        PropertyValue::String(PropertyText::from_bytes(input)),
                    )
                    .unwrap(),
            );
            let before = atom.clone();
            let mut diagnostics = Vec::new();
            let label =
                get_atom_label(&atom, 0, Hybridization::S, || false, &mut diagnostics).unwrap();
            assert_eq!(label.as_bytes(), expected);
            assert_eq!(parameters.get(&label).is_some(), expected == b"C_3");
            assert!(diagnostics.is_empty());
            assert_eq!(atom, before);
        }
    }
}
