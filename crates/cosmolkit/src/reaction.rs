//! Reaction result transport; chemistry remains in the detached reaction owner.

use crate::Molecule;

/// The source bool and the molecule finalized by the sole operation runtime.
/// Kept private until canonical public registration and integration review.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ReactionApplyResult {
    pub(crate) molecule: Molecule,
    pub(crate) changed: bool,
}

impl From<(Molecule, bool)> for ReactionApplyResult {
    fn from((molecule, changed): (Molecule, bool)) -> Self {
        // Preserve the algorithm's bool, including a true result with no
        // graph difference. Conversion owns the already finalized molecule;
        // it cannot validate, commit, or infer chemistry from graph equality.
        Self { molecule, changed }
    }
}

#[cfg(feature = "cap-reaction")]
pub(crate) fn assemble_product_sets(
    molecules: Vec<Molecule>,
    lengths: Vec<usize>,
) -> Result<Vec<Vec<Molecule>>, crate::OperationError> {
    let mut remaining = molecules.into_iter();
    let mut sets = Vec::with_capacity(lengths.len());
    for length in lengths {
        if length > remaining.len() {
            return Err(crate::OperationError::InvalidAlgorithmResult {
                operation: "reaction_products",
                field: "product set length",
                actual: length,
                expected: remaining.len(),
            });
        }
        sets.push(remaining.by_ref().take(length).collect());
    }
    if remaining.len() != 0 {
        return Err(crate::OperationError::InvalidAlgorithmResult {
            operation: "reaction_products",
            field: "ungrouped products",
            actual: remaining.len(),
            expected: 0,
        });
    }
    Ok(sets)
}
