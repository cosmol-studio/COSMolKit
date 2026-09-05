use std::ops::{BitOr, BitOrAssign};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SanitizeOps(u32);

impl SanitizeOps {
    pub const NONE: Self = Self(0);
    pub const CLEANUP: Self = Self(1 << 0);
    pub const PROPERTIES: Self = Self(1 << 1);
    pub const SYMMRINGS: Self = Self(1 << 2);
    pub const KEKULIZE: Self = Self(1 << 3);
    pub const FIND_RADICALS: Self = Self(1 << 4);
    pub const SET_AROMATICITY: Self = Self(1 << 5);
    pub const SET_CONJUGATION: Self = Self(1 << 6);
    pub const SET_HYBRIDIZATION: Self = Self(1 << 7);
    pub const CLEANUP_CHIRALITY: Self = Self(1 << 8);
    pub const ADJUST_HYDROGENS: Self = Self(1 << 9);
    pub const CLEANUP_ORGANOMETALLICS: Self = Self(1 << 10);
    pub const CLEANUP_ATROPISOMERS: Self = Self(1 << 11);
    pub const ALL: Self = Self((1 << 12) - 1);

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }
}

impl BitOr for SanitizeOps {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl BitOrAssign for SanitizeOps {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SanitizeStep {
    Cleanup,
    CleanupOrganometallics,
    Properties,
    SymmRings,
    Kekulize,
    FindRadicals,
    SetAromaticity,
    SetConjugation,
    SetHybridization,
    CleanupAtropisomers,
    CleanupChirality,
    AdjustHydrogens,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("sanitize failed at {step:?}: {message}")]
pub struct SanitizeError {
    pub step: SanitizeStep,
    pub message: String,
    pub unsupported: Option<crate::UnsupportedFeatureError>,
}

const DETECT_CHEMISTRY_PROBLEM_STAGES: &[(SanitizeStep, SanitizeOps)] = &[
    (SanitizeStep::Cleanup, SanitizeOps::CLEANUP),
    (
        SanitizeStep::CleanupOrganometallics,
        SanitizeOps::CLEANUP_ORGANOMETALLICS,
    ),
    (SanitizeStep::Properties, SanitizeOps::PROPERTIES),
    (SanitizeStep::SymmRings, SanitizeOps::SYMMRINGS),
    (SanitizeStep::Kekulize, SanitizeOps::KEKULIZE),
    (SanitizeStep::FindRadicals, SanitizeOps::FIND_RADICALS),
    (SanitizeStep::SetAromaticity, SanitizeOps::SET_AROMATICITY),
    (SanitizeStep::SetConjugation, SanitizeOps::SET_CONJUGATION),
    (
        SanitizeStep::SetHybridization,
        SanitizeOps::SET_HYBRIDIZATION,
    ),
    (
        SanitizeStep::CleanupAtropisomers,
        SanitizeOps::CLEANUP_ATROPISOMERS,
    ),
    (
        SanitizeStep::CleanupChirality,
        SanitizeOps::CLEANUP_CHIRALITY,
    ),
    (SanitizeStep::AdjustHydrogens, SanitizeOps::ADJUST_HYDROGENS),
];

#[must_use]
pub fn detect_chemistry_problems(
    molecule: &crate::Molecule,
    ops: SanitizeOps,
) -> Vec<SanitizeError> {
    let mut working = molecule.clone();
    let mut problems = Vec::new();

    for (step, stage_ops) in DETECT_CHEMISTRY_PROBLEM_STAGES {
        if !ops.contains(*stage_ops) {
            continue;
        }
        match working.sanitize_with_ops(*stage_ops) {
            Ok(updated) => working = updated,
            Err(crate::OperationError::Sanitize { source, .. }) => problems.push(source),
            Err(crate::OperationError::UnsupportedFeature { source, .. }) => {
                problems.push(SanitizeError {
                    step: *step,
                    message: source.to_string(),
                    unsupported: Some(source),
                });
            }
            Err(other) => problems.push(SanitizeError {
                step: *step,
                message: other.to_string(),
                unsupported: None,
            }),
        }
    }

    problems
}

#[cfg(test)]
mod tests {
    use super::{SanitizeOps, SanitizeStep, detect_chemistry_problems};
    use crate::{AtomSpec, BondOrder, BondSpec, Element, Molecule, MoleculeBuilder};

    fn aromatic_hypervalent_carbon() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let oxygen_a = builder.add_atom(AtomSpec::new(Element::O));
        let oxygen_b = builder.add_atom(AtomSpec::new(Element::O));
        let oxygen_c = builder.add_atom(AtomSpec::new(Element::O));
        for oxygen in [oxygen_a, oxygen_b, oxygen_c] {
            builder
                .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Double))
                .unwrap();
        }
        builder.build().unwrap()
    }

    #[test]
    fn detect_chemistry_problems_accumulates_property_cache_and_kekulize_failures() {
        let molecule = aromatic_hypervalent_carbon();

        let problems =
            detect_chemistry_problems(&molecule, SanitizeOps::PROPERTIES | SanitizeOps::KEKULIZE);

        assert_eq!(problems.len(), 2);
        assert_eq!(problems[0].step, SanitizeStep::Properties);
        assert!(problems[0].message.contains("greater than permitted"));
        assert_eq!(problems[1].step, SanitizeStep::Kekulize);
        assert!(problems[1].message.contains("aromatic"));
    }

    #[test]
    fn detect_chemistry_problems_does_not_mutate_source_molecule() {
        let molecule = aromatic_hypervalent_carbon();
        let original = molecule.clone();

        let _ =
            detect_chemistry_problems(&molecule, SanitizeOps::PROPERTIES | SanitizeOps::KEKULIZE);

        assert_eq!(molecule, original);
    }
}
