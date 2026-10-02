//! Descriptor methods on the canonical runtime molecule.

use crate::{Molecule, OperationError};

/// Read-boundary error for public descriptor count queries.
///
/// [`Self::MissingPreparedValence`] means the existing validity-checked
/// runtime accessor returned `None`: the molecule carries no valid prepared
/// valence assignment (for example it was built through a raw/unsanitized
/// constructor). It has no child error. [`Self::Algorithm`] retains the
/// actual owned typed domain cause and borrows it through
/// [`std::error::Error::source`]; it never string-flattens the cause.
#[derive(Debug)]
pub enum DescriptorReadError {
    /// No valid prepared valence assignment is installed in the runtime
    /// cache; these queries never create or install one themselves.
    MissingPreparedValence,
    /// The domain owner returned its typed error; borrowed as the source.
    Algorithm {
        source: cosmolkit_descriptors::DescriptorError,
    },
}
impl std::fmt::Display for DescriptorReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingPreparedValence => write!(
                f,
                "descriptor query requires a prepared valence assignment; the molecule has no valid cached assignment"
            ),
            Self::Algorithm { .. } => write!(f, "descriptor query failed"),
        }
    }
}
impl std::error::Error for DescriptorReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::MissingPreparedValence => None,
            Self::Algorithm { source } => Some(source),
        }
    }
}

/// Private mandatory-borrow gate for the valence-backed descriptor queries.
///
/// Borrows the EXISTING validity-checked cached valence assignment from the
/// runtime molecule's derived-cache block; this helper neither creates nor
/// installs a value and
/// never recomputes valence. `None` from the accessor maps to the typed
/// [`DescriptorReadError::MissingPreparedValence`].
fn required_descriptor_valence<'a>(
    molecule: &'a Molecule,
) -> Result<&'a cosmolkit_core::ValenceAssignment, DescriptorReadError> {
    molecule
        .derived_cache_runtime()
        .valence_assignment()
        .ok_or(DescriptorReadError::MissingPreparedValence)
}

fn descriptor_error(
    operation: &'static str,
    error: cosmolkit_descriptors::DescriptorError,
) -> OperationError {
    OperationError::Algorithm {
        operation,
        detail: error.to_string(),
    }
}

impl Molecule {
    /// Returns the number of heavy atoms (RDKit `CalcNumHeavyAtoms`).
    ///
    /// Topology-only read-only query: the source closure touches no
    /// hydrogen-count, valence or ring state, so no prepared valence is
    /// required or borrowed.
    #[must_use]
    pub fn num_heavy_atoms(&self) -> Result<u32, DescriptorReadError> {
        cosmolkit_descriptors::num_heavy_atoms(self.topology())
            .map_err(|source| DescriptorReadError::Algorithm { source })
    }

    /// Returns the total atom count including implicit and
    /// explicit-property hydrogens (RDKit `CalcNumAtoms`).
    ///
    /// Borrows the MANDATORY cached valence assignment through
    /// [`required_descriptor_valence`]; never recomputes or installs one.
    /// Explicit H atom rows count once as rows (source
    /// `includeNeighbors=false`); `Molecule::num_atoms` keeps its separate
    /// atom-table-length meaning.
    #[must_use]
    pub fn total_atom_count(&self) -> Result<u32, DescriptorReadError> {
        let assignment = required_descriptor_valence(self)?;
        cosmolkit_descriptors::total_atom_count_with_valence(self.topology(), assignment)
            .map_err(|source| DescriptorReadError::Algorithm { source })
    }

    /// Returns the Lipinski hydrogen-bond acceptor count
    /// (RDKit `CalcNumLipinskiHBA`).
    ///
    /// Topology-only DIRECT N/O row count — not the general recursive
    /// SMARTS-based NumHBA; no prepared valence is required or borrowed.
    #[must_use]
    pub fn lipinski_hba(&self) -> Result<u32, DescriptorReadError> {
        cosmolkit_descriptors::lipinski_hba(self.topology())
            .map_err(|source| DescriptorReadError::Algorithm { source })
    }

    /// Returns the Lipinski hydrogen-bond donor hydrogen sum
    /// (RDKit `CalcNumLipinskiHBD`).
    ///
    /// The hydrogen sum on N/O rows (NOT the general SMARTS-based NumHBD
    /// and not the donor-atom count). Borrows the MANDATORY cached
    /// valence assignment; never recomputes or installs one.
    #[must_use]
    pub fn lipinski_hbd(&self) -> Result<u32, DescriptorReadError> {
        let assignment = required_descriptor_valence(self)?;
        cosmolkit_descriptors::lipinski_hbd_with_valence(self.topology(), assignment)
            .map_err(|source| DescriptorReadError::Algorithm { source })
    }

    /// Returns the fraction of CSP3 carbons (RDKit `CalcFractionCSP3`).
    ///
    /// Carbon rows whose SOURCE total degree (explicit neighbors +
    /// implicit/explicit-property hydrogens, `includeNeighbors=false`) is
    /// exactly four — never a hybridization flag; zero carbons yield the
    /// literal 0.0. Borrows the MANDATORY cached valence assignment; no
    /// topology, property or cache write occurs.
    #[must_use]
    pub fn fraction_csp3(&self) -> Result<f64, DescriptorReadError> {
        let assignment = required_descriptor_valence(self)?;
        cosmolkit_descriptors::fraction_csp3_with_valence(self.topology(), assignment)
            .map_err(|source| DescriptorReadError::Algorithm { source })
    }

    /// Returns the RDKit-compatible average molecular weight.
    #[must_use]
    pub fn molecular_weight(&self) -> Result<f64, OperationError> {
        self.molecular_weight_with_options(false)
    }

    /// Returns average molecular weight with an explicit heavy-atom mode.
    #[must_use]
    pub fn molecular_weight_with_options(&self, only_heavy: bool) -> Result<f64, OperationError> {
        cosmolkit_descriptors::molecular_weight_with_valence(
            self.topology(),
            only_heavy,
            self.derived_cache_runtime().valence_assignment(),
        )
        .map_err(|error| descriptor_error("molecular_weight", error))
    }

    /// Returns the RDKit-compatible exact molecular weight.
    #[must_use]
    pub fn exact_molecular_weight(&self) -> Result<f64, OperationError> {
        self.exact_molecular_weight_with_options(false)
    }

    /// Returns exact molecular weight with an explicit heavy-atom mode.
    #[must_use]
    pub fn exact_molecular_weight_with_options(
        &self,
        only_heavy: bool,
    ) -> Result<f64, OperationError> {
        cosmolkit_descriptors::exact_molecular_weight_with_valence(
            self.topology(),
            only_heavy,
            self.derived_cache_runtime().valence_assignment(),
        )
        .map_err(|error| descriptor_error("exact_molecular_weight", error))
    }

    /// Returns the Hill-ordered molecular formula.
    #[must_use]
    pub fn molecular_formula(&self) -> Result<String, OperationError> {
        self.molecular_formula_with_options(false, false)
    }

    /// Returns a molecular formula with isotope formatting controls.
    #[must_use]
    pub fn molecular_formula_with_options(
        &self,
        separate_isotopes: bool,
        abbreviate_h_isotopes: bool,
    ) -> Result<String, OperationError> {
        cosmolkit_descriptors::molecular_formula_with_valence(
            self.topology(),
            separate_isotopes,
            abbreviate_h_isotopes,
            self.derived_cache_runtime().valence_assignment(),
        )
        .map_err(|error| descriptor_error("molecular_formula", error))
    }
}

#[cfg(test)]
mod descriptor_public_error_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondOrder, BondSpec, Element};
    use std::error::Error as _;

    /// Raw unsanitized constructor (builder only, no sanitizing stage):
    /// the derived cache holds NO valid valence assignment.
    fn raw_cco() -> Molecule {
        let mut builder = crate::MoleculeBuilder::new();
        let c0 = builder.add_atom(AtomSpec::new(Element::C));
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let o2 = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(c0, c1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, o2, BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn descriptor_public_error_helper_missing_cache() {
        // The helper's REAL result on a live molecule whose cache was never
        // prepared: the existing validity-checked accessor returns None and
        // the helper maps it to the typed MissingPreparedValence with no
        // child error. This is helper transport on a legitimately raw
        // molecule — a corrupt live molecule is not reachable here.
        let molecule = raw_cco();
        assert!(
            molecule
                .derived_cache_runtime()
                .valence_assignment()
                .is_none(),
            "raw builder molecule carries no valid prepared valence"
        );
        let err = required_descriptor_valence(&molecule).unwrap_err();
        assert!(matches!(err, DescriptorReadError::MissingPreparedValence));
        assert!(err.source().is_none(), "no fabricated cause");
    }

    #[test]
    fn descriptor_public_error_algorithm_source_chain() {
        // A REAL detached adapter failure (malformed explicit rows on the
        // borrowed CCO topology) is retained as the owned typed cause and
        // BORROWED through Error::source — never string-flattened. The
        // public query can never observe this state on a live molecule
        // because the cache validity gate precedes the adapter call; this
        // proof exercises the error transport, not live corruption.
        let molecule = raw_cco();
        let bad = cosmolkit_core::ValenceAssignment {
            explicit_valence: vec![4, 1],
            implicit_hydrogens: vec![2, 1, 1],
        };
        let source =
            cosmolkit_descriptors::total_atom_count_with_valence(molecule.topology(), &bad)
                .unwrap_err();
        let err = DescriptorReadError::Algorithm { source };
        assert_eq!(err.to_string(), "descriptor query failed");
        let borrowed = err
            .source()
            .expect("typed cause is borrowed")
            .downcast_ref::<cosmolkit_descriptors::DescriptorError>()
            .expect("source chain carries the actual domain error type");
        assert_eq!(
            borrowed,
            &cosmolkit_descriptors::DescriptorError::InvalidValenceRows {
                function: "num_atoms",
                field: "explicit_valence",
                actual: 2,
                expected: 3,
            }
        );
    }
}

#[cfg(test)]
mod descriptor_public_storage_tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn descriptor_public_storage_queries_share_immutable_arc_blocks() {
        // Private storage proof: CCO x five queries x two sharing
        // receivers (original + Arc-sharing peer clone) x two repeats =
        // 20 ACTUAL query calls. Before and after EVERY call the FOUR
        // actual MoleculeState Arc blocks (topology, coordinates,
        // properties, derived cache) keep never-refreshed identities and
        // whole immutable values, the original/peer share all four
        // blocks, the cached valence assignment stays valid and the
        // borrowed assignment identity is stable. No public test probe
        // or production observer is involved.
        const CCO_HEAVY: u32 = 3;
        const CCO_TOTAL: u32 = 9;
        const CCO_HBA: u32 = 1;
        const CCO_HBD: u32 = 1;
        const CCO_CSP3_BITS: u64 = 0x3ff0_0000_0000_0000;

        let original = Molecule::from_smiles("CCO").unwrap();
        let peer = original.clone();

        // Sharing: the peer clone shares all FOUR Arc blocks.
        assert!(Arc::ptr_eq(
            &original.topology_arc_runtime(),
            &peer.topology_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &original.coordinates_arc_runtime(),
            &peer.coordinates_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &original.properties_arc_runtime(),
            &peer.properties_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &original.derived_cache_arc_runtime(),
            &peer.derived_cache_arc_runtime()
        ));
        assert_eq!(original.topology(), peer.topology());

        // Never-refreshed identities, captured once before any query.
        let topology_arc = original.topology_arc_runtime();
        let coordinates_arc = original.coordinates_arc_runtime();
        let properties_arc = original.properties_arc_runtime();
        let cache_arc = original.derived_cache_arc_runtime();
        // Once-captured complete baseline values and the original valid
        // assignment pointer, alongside the four Arc baselines above.
        let baseline_topology = original.topology().clone();
        let baseline_coordinates = original.coordinate_block_runtime().clone();
        let baseline_properties = original.properties().clone();
        let baseline_cache = original.derived_cache_runtime().clone();
        let baseline_assignment = original
            .derived_cache_runtime()
            .valence_assignment()
            .expect("baseline valence valid before every query");

        // ONE checkpoint closure: BOTH receivers keep all four
        // never-refreshed Arc identities and complete baseline values,
        // share all four Arc blocks with each other, and carry valid
        // assignments pointer-identical to the original baseline.
        // Invoked immediately BEFORE and AFTER each of the twenty
        // actual query calls.
        let checkpoint = |stage: &str| {
            for (who, molecule) in [("original", &original), ("peer", &peer)] {
                assert!(
                    Arc::ptr_eq(&topology_arc, &molecule.topology_arc_runtime()),
                    "{stage} {who}: topology Arc identity"
                );
                assert!(
                    Arc::ptr_eq(&coordinates_arc, &molecule.coordinates_arc_runtime()),
                    "{stage} {who}: coordinates Arc identity"
                );
                assert!(
                    Arc::ptr_eq(&properties_arc, &molecule.properties_arc_runtime()),
                    "{stage} {who}: properties Arc identity"
                );
                assert!(
                    Arc::ptr_eq(&cache_arc, &molecule.derived_cache_arc_runtime()),
                    "{stage} {who}: derived-cache Arc identity"
                );
                assert_eq!(molecule.topology(), &baseline_topology, "{stage} {who}");
                assert_eq!(
                    molecule.coordinate_block_runtime(),
                    &baseline_coordinates,
                    "{stage} {who}"
                );
                assert_eq!(molecule.properties(), &baseline_properties, "{stage} {who}");
                assert_eq!(
                    molecule.derived_cache_runtime(),
                    &baseline_cache,
                    "{stage} {who}"
                );
                let assignment = molecule
                    .derived_cache_runtime()
                    .valence_assignment()
                    .unwrap_or_else(|| panic!("{stage} {who}: valence must stay valid"));
                assert!(
                    std::ptr::eq(baseline_assignment, assignment),
                    "{stage} {who}: assignment pointer-identical to baseline"
                );
            }
            assert!(
                Arc::ptr_eq(
                    &original.topology_arc_runtime(),
                    &peer.topology_arc_runtime()
                ),
                "{stage}: original/peer topology sharing"
            );
            assert!(
                Arc::ptr_eq(
                    &original.coordinates_arc_runtime(),
                    &peer.coordinates_arc_runtime()
                ),
                "{stage}: original/peer coordinates sharing"
            );
            assert!(
                Arc::ptr_eq(
                    &original.properties_arc_runtime(),
                    &peer.properties_arc_runtime()
                ),
                "{stage}: original/peer properties sharing"
            );
            assert!(
                Arc::ptr_eq(
                    &original.derived_cache_arc_runtime(),
                    &peer.derived_cache_arc_runtime()
                ),
                "{stage}: original/peer derived-cache sharing"
            );
        };

        let run = |name: &str, receiver: &Molecule| -> f64 {
            match name {
                "num_heavy_atoms" => f64::from(receiver.num_heavy_atoms().unwrap()),
                "total_atom_count" => f64::from(receiver.total_atom_count().unwrap()),
                "lipinski_hba" => f64::from(receiver.lipinski_hba().unwrap()),
                "lipinski_hbd" => f64::from(receiver.lipinski_hbd().unwrap()),
                "fraction_csp3" => receiver.fraction_csp3().unwrap(),
                _ => unreachable!("frozen five-query table"),
            }
        };
        let expected = |name: &str| -> f64 {
            match name {
                "num_heavy_atoms" => f64::from(CCO_HEAVY),
                "total_atom_count" => f64::from(CCO_TOTAL),
                "lipinski_hba" => f64::from(CCO_HBA),
                "lipinski_hbd" => f64::from(CCO_HBD),
                "fraction_csp3" => f64::from_bits(CCO_CSP3_BITS),
                _ => unreachable!("frozen five-query table"),
            }
        };

        let mut calls = 0usize;
        for name in [
            "num_heavy_atoms",
            "total_atom_count",
            "lipinski_hba",
            "lipinski_hbd",
            "fraction_csp3",
        ] {
            for (label, receiver) in [("original", &original), ("peer", &peer)] {
                for repeat in 0..2 {
                    checkpoint("before");
                    let topology_before = receiver.topology().clone();
                    let coordinates_before = receiver.coordinate_block_runtime().clone();
                    let properties_before = receiver.properties().clone();
                    let cache_before = receiver.derived_cache_runtime().clone();
                    let assignment_before = receiver
                        .derived_cache_runtime()
                        .valence_assignment()
                        .expect("valence valid before the call");

                    let result = run(name, receiver);
                    calls += 1;
                    checkpoint("after");

                    assert_eq!(
                        result.to_bits(),
                        expected(name).to_bits(),
                        "{name} {label} #{repeat}"
                    );
                    // Never-refreshed identities for all four Arc blocks.
                    assert!(
                        Arc::ptr_eq(&topology_arc, &receiver.topology_arc_runtime()),
                        "{name} {label} #{repeat}: topology Arc never refreshed"
                    );
                    assert!(
                        Arc::ptr_eq(&coordinates_arc, &receiver.coordinates_arc_runtime()),
                        "{name} {label} #{repeat}: coordinates Arc never refreshed"
                    );
                    assert!(
                        Arc::ptr_eq(&properties_arc, &receiver.properties_arc_runtime()),
                        "{name} {label} #{repeat}: properties Arc never refreshed"
                    );
                    assert!(
                        Arc::ptr_eq(&cache_arc, &receiver.derived_cache_arc_runtime()),
                        "{name} {label} #{repeat}: derived-cache Arc never refreshed"
                    );
                    // Whole immutable values unchanged; topology equal.
                    assert_eq!(receiver.topology(), &topology_before, "{name} {label}");
                    assert_eq!(
                        receiver.coordinate_block_runtime(),
                        &coordinates_before,
                        "{name} {label}"
                    );
                    assert_eq!(receiver.properties(), &properties_before, "{name} {label}");
                    assert_eq!(
                        receiver.derived_cache_runtime(),
                        &cache_before,
                        "{name} {label}"
                    );
                    assert_eq!(original.topology(), peer.topology());
                    // Valence stays valid with a stable borrowed identity.
                    let assignment_after = receiver
                        .derived_cache_runtime()
                        .valence_assignment()
                        .expect("valence valid after the call");
                    assert!(
                        std::ptr::eq(assignment_before, assignment_after),
                        "{name} {label} #{repeat}: borrowed assignment identity stable"
                    );
                }
            }
        }
        assert_eq!(calls, 20, "exact 20-call storage census");
    }
}
