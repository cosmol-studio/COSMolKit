//! Descriptor methods on the canonical runtime molecule.
//!
//! # Ring-state queries
//!
//! The eleven `num_*ring*` queries are ring-STATE reads over the molecule's
//! derived-cache ordinary ring assignment: they never calculate, install or
//! upgrade rings, clone the cache, or require a prepared valence. A
//! sanitize-enabled constructor (or `with_assigned_rings`) installs valid
//! rows. `num_rings` requires initialized rows and fails with
//! [`DescriptorReadError::MissingInitializedRings`] on absence/reset; the
//! ten classifiers return the source-defined EMPTY-ROW `0` on absence — a
//! documented empty input, not a swallowed error. These are Experimental
//! commitments; Python/JS projections are declared, not implemented.

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
    /// No valid initialized ordinary ring state is installed in the
    /// runtime cache; the ring-count query never creates, installs or
    /// upgrades one. It has no child error.
    MissingInitializedRings,
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
            Self::MissingInitializedRings => write!(
                f,
                "ring query requires initialized ring state; the molecule has no valid cached ring assignment"
            ),
            Self::Algorithm { .. } => write!(f, "descriptor query failed"),
        }
    }
}
impl std::error::Error for DescriptorReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::MissingPreparedValence | Self::MissingInitializedRings => None,
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

    /// Returns the number of rings (RDKit `CalcNumRings`).
    ///
    /// Ring-state read-only query: requires the molecule's VALID
    /// initialized ordinary ring state and delegates exactly once to the
    /// narrow borrowed-row owner. Legitimate absence/reset maps to the
    /// typed [`DescriptorReadError::MissingInitializedRings`]; an
    /// initialized EMPTY row set is the literal 0. This query never
    /// calculates, installs or upgrades rings, clones no cache and reads
    /// no valence.
    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_rings(&self) -> Result<u32, DescriptorReadError> {
        let rings = self
            .derived_cache_runtime()
            .valid_ring_info()
            .ok_or(DescriptorReadError::MissingInitializedRings)?;
        cosmolkit_descriptors::num_rings_with_ring_info(rings)
            .map_err(|source| DescriptorReadError::Algorithm { source })
    }

    /// Returns the number of heterocycles (RDKit `CalcNumHeterocycles`).
    ///
    /// Borrows the molecule's EXISTING initialized ordinary ring rows into
    /// the narrow owner; it never finds, installs or upgrades rings,
    /// clones no cache and reads no valence. Legitimate absence/reset
    /// returns the source-defined EMPTY-ROW result 0 (the pinned
    /// classifier scans the stored rows only) — an explicitly documented
    /// empty input, not a swallowed error.
    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_heterocycles(&self) -> Result<u32, DescriptorReadError> {
        match self.derived_cache_runtime().valid_ring_info() {
            Some(rings) => {
                cosmolkit_descriptors::num_heterocycles_with_ring_info(self.topology(), rings)
                    .map_err(|source| DescriptorReadError::Algorithm { source })
            }
            None => Ok(0),
        }
    }

    /// Returns the number of aromatic rings (RDKit `CalcNumAromaticRings`).
    ///
    /// Borrows the molecule's EXISTING initialized ordinary ring rows into
    /// the narrow owner; it never finds, installs or upgrades rings, clones
    /// no cache and reads no valence. Legitimate absence/reset returns the
    /// source-defined EMPTY-ROW result 0 — an explicitly documented empty
    /// input, not a swallowed error.
    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_aromatic_rings(&self) -> Result<u32, DescriptorReadError> {
        match self.derived_cache_runtime().valid_ring_info() {
            Some(rings) => {
                cosmolkit_descriptors::num_aromatic_rings_with_ring_info(self.topology(), rings)
                    .map_err(|source| DescriptorReadError::Algorithm { source })
            }
            None => Ok(0),
        }
    }

    /// Returns the number of saturated rings
    /// (RDKit `CalcNumSaturatedRings`).
    ///
    /// Borrows the molecule's EXISTING initialized ordinary ring rows into
    /// the narrow owner; it never finds, installs or upgrades rings, clones
    /// no cache and reads no valence. Legitimate absence/reset returns the
    /// source-defined EMPTY-ROW result 0 — an explicitly documented empty
    /// input, not a swallowed error.
    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_saturated_rings(&self) -> Result<u32, DescriptorReadError> {
        match self.derived_cache_runtime().valid_ring_info() {
            Some(rings) => {
                cosmolkit_descriptors::num_saturated_rings_with_ring_info(self.topology(), rings)
                    .map_err(|source| DescriptorReadError::Algorithm { source })
            }
            None => Ok(0),
        }
    }

    /// Returns the number of aliphatic rings
    /// (RDKit `CalcNumAliphaticRings`).
    ///
    /// Borrows the molecule's EXISTING initialized ordinary ring rows into
    /// the narrow owner; it never finds, installs or upgrades rings, clones
    /// no cache and reads no valence. Legitimate absence/reset returns the
    /// source-defined EMPTY-ROW result 0 — an explicitly documented empty
    /// input, not a swallowed error.
    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_aliphatic_rings(&self) -> Result<u32, DescriptorReadError> {
        match self.derived_cache_runtime().valid_ring_info() {
            Some(rings) => {
                cosmolkit_descriptors::num_aliphatic_rings_with_ring_info(self.topology(), rings)
                    .map_err(|source| DescriptorReadError::Algorithm { source })
            }
            None => Ok(0),
        }
    }

    // The six combined heterocycle/carbocycle queries share one thin
    // borrowed-read shape: delegate to the matching narrow owner with the
    // EXISTING initialized rows; absence/reset is the documented
    // source-defined empty-row 0, never a swallowed error; no query
    // finds, installs or upgrades rings, clones the cache or reads
    // valence.
    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_aromatic_heterocycles(&self) -> Result<u32, DescriptorReadError> {
        match self.derived_cache_runtime().valid_ring_info() {
            Some(rings) => cosmolkit_descriptors::num_aromatic_heterocycles_with_ring_info(
                self.topology(),
                rings,
            )
            .map_err(|source| DescriptorReadError::Algorithm { source }),
            None => Ok(0),
        }
    }

    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_aromatic_carbocycles(&self) -> Result<u32, DescriptorReadError> {
        match self.derived_cache_runtime().valid_ring_info() {
            Some(rings) => cosmolkit_descriptors::num_aromatic_carbocycles_with_ring_info(
                self.topology(),
                rings,
            )
            .map_err(|source| DescriptorReadError::Algorithm { source }),
            None => Ok(0),
        }
    }

    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_aliphatic_heterocycles(&self) -> Result<u32, DescriptorReadError> {
        match self.derived_cache_runtime().valid_ring_info() {
            Some(rings) => cosmolkit_descriptors::num_aliphatic_heterocycles_with_ring_info(
                self.topology(),
                rings,
            )
            .map_err(|source| DescriptorReadError::Algorithm { source }),
            None => Ok(0),
        }
    }

    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_aliphatic_carbocycles(&self) -> Result<u32, DescriptorReadError> {
        match self.derived_cache_runtime().valid_ring_info() {
            Some(rings) => cosmolkit_descriptors::num_aliphatic_carbocycles_with_ring_info(
                self.topology(),
                rings,
            )
            .map_err(|source| DescriptorReadError::Algorithm { source }),
            None => Ok(0),
        }
    }

    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_saturated_heterocycles(&self) -> Result<u32, DescriptorReadError> {
        match self.derived_cache_runtime().valid_ring_info() {
            Some(rings) => cosmolkit_descriptors::num_saturated_heterocycles_with_ring_info(
                self.topology(),
                rings,
            )
            .map_err(|source| DescriptorReadError::Algorithm { source }),
            None => Ok(0),
        }
    }

    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_saturated_carbocycles(&self) -> Result<u32, DescriptorReadError> {
        match self.derived_cache_runtime().valid_ring_info() {
            Some(rings) => cosmolkit_descriptors::num_saturated_carbocycles_with_ring_info(
                self.topology(),
                rings,
            )
            .map_err(|source| DescriptorReadError::Algorithm { source }),
            None => Ok(0),
        }
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

/// Q1 state discriminators on the SAME cube topology: public absence
/// semantics split by query class (num_rings requires initialized rows;
/// num_heterocycles returns the source-defined empty-row 0), plus the
/// real Fast5/Sssr5/Symm6 row sets through the private constructor seam.
#[cfg(all(test, feature = "cap-smiles", feature = "cap-rings"))]
mod ring_live_public_q1_tests {
    use super::*;

    fn raw(input: &str) -> Molecule {
        Molecule::from_smiles_with_params(
            input,
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hydrogens: false,
                ..cosmolkit_smiles::SmilesParseParams::default()
            },
        )
        .unwrap()
    }

    fn supplied(input: &str, rings: Option<cosmolkit_core::RingInfo>) -> Molecule {
        let base = raw(input);
        Molecule::from_smiles_parts_with_derived_state(
            base.topology().clone(),
            base.coordinate_block_runtime().clone(),
            base.properties().clone(),
            None,
            rings,
        )
        .unwrap()
    }

    #[test]
    fn ring_live_public_q1_state_discriminators() {
        // C7 on the cube topology (E-V+1 = 5 SSSR rings, 6 symmetrized).
        // NOTE: the Some(reset/uninitialized) arm is structurally
        // unreachable outside cosmolkit-core (no public constructor leaves
        // initialized == false; RingInfo::reset is pub(crate)); the seam
        // clears it by construction, making it observationally identical
        // to the absent arm proven here. No fake reset API is introduced.
        let topology = raw("C12C3C4C1C5C2C3C45").topology().clone();
        let fast5 = cosmolkit_core::fast_find_rings(&topology).unwrap();
        let sssr5 =
            cosmolkit_core::find_sssr(&topology, &cosmolkit_core::RingSearchParams::default())
                .unwrap();
        let symm6 = cosmolkit_core::symmetrized_sssr(
            &topology,
            &cosmolkit_core::RingSearchParams::default(),
        )
        .unwrap();
        assert_eq!(fast5.atom_rings().len(), 5);
        assert_eq!(sssr5.atom_rings().len(), 5);
        assert_eq!(symm6.atom_rings().len(), 6);
        let mut calls = 0usize;
        for (name, supply, want_rings) in [
            ("absent", None, None),
            (
                "other-empty",
                Some(cosmolkit_core::RingInfo::new(
                    cosmolkit_core::RingFindType::OtherOrUnknown,
                    8,
                    12,
                )),
                Some(0),
            ),
            ("fast5", Some(fast5.clone()), Some(5)),
            ("sssr5", Some(sssr5.clone()), Some(5)),
            ("symm6", Some(symm6.clone()), Some(6)),
        ] {
            let label = name.to_string();
            let molecule = supplied("C12C3C4C1C5C2C3C45", supply);
            calls += 1;
            match want_rings {
                None => {
                    // num_rings: legitimate absence/reset is the typed
                    // MissingInitializedRings with no child error.
                    let error = molecule.num_rings().unwrap_err();
                    assert!(
                        matches!(error, DescriptorReadError::MissingInitializedRings),
                        "{label}: got {error:?}"
                    );
                }
                Some(expected) => {
                    assert_eq!(molecule.num_rings().unwrap(), expected, "{label}");
                }
            }
            calls += 1;
            // num_heterocycles: absence is the source-defined empty-row 0;
            // the all-carbon cube also yields 0 with real rows installed.
            assert_eq!(molecule.num_heterocycles().unwrap(), 0, "{label}");
        }
        assert_eq!(calls, 10, "exact census");
    }
}

/// C9: one default benzene and its peer, each eleven methods x two repeats
/// = 44 real queries. The four Arc blocks and full values/validity/quality/
/// paired rows/memberships are checked before AND after EACH call against
/// never-refreshed baselines; the root acquisition counter proves no query
/// ran a finder.
#[cfg(all(
    test,
    feature = "cap-smiles",
    feature = "cap-rings",
    feature = "cap-descriptors"
))]
mod ring_live_public_storage_tests {
    use super::*;
    use crate::AtomId;
    use crate::BondId;
    use crate::DerivedState;

    #[test]
    fn ring_live_public_storage_four_block_repeat_peer_proof() {
        let molecule = Molecule::from_smiles("c1ccccc1").unwrap();
        let peer = molecule.clone();
        // Never-refreshed baselines.
        let topology_arc = molecule.topology_arc_runtime();
        let coordinates_arc = molecule.coordinates_arc_runtime();
        let properties_arc = molecule.properties_arc_runtime();
        let cache_arc = molecule.derived_cache_arc_runtime();
        let topology_value = molecule.topology().clone();
        let coordinates_value = molecule.coordinate_block_runtime().clone();
        let properties_value = molecule.properties().clone();
        let cache_value = cache_arc.as_ref().clone();

        let queries: [(&str, fn(&Molecule) -> Result<u32, DescriptorReadError>, u32); 11] = [
            ("num_rings", Molecule::num_rings, 1),
            ("num_heterocycles", Molecule::num_heterocycles, 0),
            ("num_aromatic_rings", Molecule::num_aromatic_rings, 1),
            ("num_saturated_rings", Molecule::num_saturated_rings, 0),
            ("num_aliphatic_rings", Molecule::num_aliphatic_rings, 0),
            (
                "num_aromatic_heterocycles",
                Molecule::num_aromatic_heterocycles,
                0,
            ),
            (
                "num_aromatic_carbocycles",
                Molecule::num_aromatic_carbocycles,
                1,
            ),
            (
                "num_aliphatic_heterocycles",
                Molecule::num_aliphatic_heterocycles,
                0,
            ),
            (
                "num_aliphatic_carbocycles",
                Molecule::num_aliphatic_carbocycles,
                0,
            ),
            (
                "num_saturated_heterocycles",
                Molecule::num_saturated_heterocycles,
                0,
            ),
            (
                "num_saturated_carbocycles",
                Molecule::num_saturated_carbocycles,
                0,
            ),
        ];

        let check_state = |target: &Molecule, label: &str| {
            assert!(
                std::sync::Arc::ptr_eq(&target.topology_arc_runtime(), &topology_arc),
                "{label}"
            );
            assert!(
                std::sync::Arc::ptr_eq(&target.coordinates_arc_runtime(), &coordinates_arc),
                "{label}"
            );
            assert!(
                std::sync::Arc::ptr_eq(&target.properties_arc_runtime(), &properties_arc),
                "{label}"
            );
            assert!(
                std::sync::Arc::ptr_eq(&target.derived_cache_arc_runtime(), &cache_arc),
                "{label}"
            );
            assert_eq!(target.topology(), &topology_value, "{label}");
            assert_eq!(
                target.coordinate_block_runtime(),
                &coordinates_value,
                "{label}"
            );
            assert_eq!(target.properties(), &properties_value, "{label}");
            let cache = target.derived_cache_runtime();
            assert_eq!(cache, &cache_value, "{label}");
            assert!(
                cache.valid_states().contains(DerivedState::RINGS),
                "{label}"
            );
            let rings = cache.valid_ring_info().expect("{label}: installed");
            assert!(rings.is_initialized(), "{label}");
            assert_eq!(
                rings.find_type(),
                cosmolkit_core::RingFindType::SymmSssr,
                "{label}"
            );
            assert_eq!(rings.atom_rings().len(), 1, "{label}: rows");
            assert_eq!(rings.bond_rings().len(), 1, "{label}: bond rows");
            let mut atoms_row: Vec<usize> = rings.atom_rings()[0]
                .iter()
                .map(|atom| atom.index())
                .collect();
            atoms_row.sort_unstable();
            let mut bonds_row: Vec<usize> = rings.bond_rings()[0]
                .iter()
                .map(|bond| bond.index())
                .collect();
            bonds_row.sort_unstable();
            assert_eq!(atoms_row, vec![0, 1, 2, 3, 4, 5], "{label}");
            assert_eq!(bonds_row, vec![0, 1, 2, 3, 4, 5], "{label}");
            for index in 0..6usize {
                assert_eq!(
                    rings.atom_members(AtomId::new(index)),
                    &[0],
                    "{label}: member {index}"
                );
                assert_eq!(
                    rings.bond_members(BondId::new(index)),
                    &[0],
                    "{label}: bond member {index}"
                );
            }
        };

        let mut calls = 0usize;
        let acquisitions_before = crate::ops::ring_aromaticity_probe::acquisitions();
        for target in [&molecule, &peer] {
            for (name, query, expected) in queries {
                for repeat in 0..2 {
                    let label = format!("{name}/rep{repeat}");
                    check_state(target, &format!("{label}: before"));
                    assert_eq!(query(target).unwrap(), expected, "{label}");
                    calls += 1;
                    check_state(target, &format!("{label}: after"));
                }
            }
        }
        assert_eq!(calls, 44, "exact census");
        // No query ran a root finder (code shape plus actual counter).
        assert_eq!(
            crate::ops::ring_aromaticity_probe::acquisitions() - acquisitions_before,
            0,
            "no query-local finder"
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
