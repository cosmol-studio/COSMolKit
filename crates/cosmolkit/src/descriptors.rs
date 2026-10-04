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
//!
//! # Topology-only SMARTS counts
//!
//! `num_heteroatoms` runs the pinned `[!#6;!#1]` pattern through the ONE
//! retained-pattern matcher with a narrow topology-only query context: it
//! reads ONLY atomic-number rows (dummy Z0 counts; isotope H does not),
//! never computes valence/ring state, and keeps the source default match
//! parameters — including `maxMatches = 1000`, so more than 1000 matching
//! rows clamp at 1000 exactly like the pinned `SubstructMatch` defaults.
//!
//! # General prepared SMARTS counts
//!
//! `num_hba` is the general recursive-SMARTS acceptor count
//! (RDKit `CalcNumHBA`, pinned `2.0.2` pattern). Unlike the
//! topology-only family above, it is a PREPARED read: it borrows BOTH the
//! existing valid valence assignment and the existing initialized
//! ordinary ring rows (absence fails `MissingPreparedValence` first, then
//! `MissingInitializedRings`; an unsanitized molecule fails before any
//! chemistry runs — a public prepared-state boundary limitation, not a
//! chemical divergence), constructs the borrowed five-block
//! descriptor input, and delegates exactly once to the domain owner. It
//! is NOT the direct Lipinski N/O sum (`lipinski_hba`): thiophene S
//! counts here, acid OH does not, and carbonyl O matches through the
//! H0-v2 branch. Experimental; Python/JS projections are declared, not
//! implemented.
//!
//! `num_hbd` is the general donor count (RDKit `CalcNumHBD`, pinned
//! `2.0.1` pattern). It is a NARROW valence-only prepared read: the fixed
//! pattern reads hydrogen-count and valence predicates but NO ring
//! predicate (ROOT-resolved from the pinned source), so it borrows ONLY
//! the existing valid valence assignment — valid valence with absent or
//! reset rings still SUCCEEDS, unlike `num_hba`. It is NOT the direct
//! Lipinski donor-hydrogen sum (`lipinski_hbd`): sulfur donors match
//! here, `[NH4+]` counts one atom here but four hydrogens there, and
//! isolated S has two H and counts 0. Experimental; Python/JS projections
//! are declared, not implemented.

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

    /// Returns the general SMARTS-based heteroatom count
    /// (RDKit `CalcNumHeteroatoms`).
    ///
    /// Topology-only query: the pinned `[!#6;!#1]` pattern reads ONLY
    /// atomic-number rows, so this thin read delegates once to the
    /// descriptor owner's topology-only path (narrow query context;
    /// default match parameters with maxMatches=1000). No prepared
    /// valence, ring state, cache clone/installation or operation
    /// declaration is involved; dummy rows count and isotope H does not.
    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_heteroatoms(&self) -> Result<u32, DescriptorReadError> {
        cosmolkit_descriptors::num_heteroatoms(self.topology())
            .map_err(|source| DescriptorReadError::Algorithm { source })
    }

    /// Returns the general SMARTS-based hydrogen-bond acceptor count
    /// (RDKit `CalcNumHBA`).
    ///
    /// Mandatory-borrow review: this prepared query borrows BOTH the
    /// EXISTING validity-checked cached valence assignment (through
    /// [`required_descriptor_valence`]) and the EXISTING valid
    /// initialized ordinary [`RingInfo`] rows (through the derived-cache
    /// valid-gated read). It never computes, installs or upgrades either
    /// payload, clones no whole cache state and exposes no fresh-state
    /// wrapper. Absence precedence: a missing prepared valence reports
    /// [`DescriptorReadError::MissingPreparedValence`] FIRST; a valid
    /// valence with absent/reset ordinary rings then reports
    /// [`DescriptorReadError::MissingInitializedRings`]. An unsanitized
    /// molecule therefore fails before any chemistry runs — a public
    /// prepared-state boundary limitation, not a source chemical
    /// divergence.
    ///
    /// Behavior review: this is NOT the direct Lipinski N/O row sum
    /// ([`Molecule::lipinski_hba`]); the pinned recursive SMARTS counts
    /// distinct atom sets (thiophene S here, acid OH not, carbonyl O by
    /// the H0-v2 branch; `c1cccn1C` -> 0 and `c1cccc(=O)n1C` -> 1 are
    /// the pinned RDKit #8997 discriminators). Cost review: borrowing
    /// the payloads does NOT make the call free — the domain owner
    /// validates rows, builds the query-match context (adjacency) and
    /// runs the recursive-pattern matcher pass; no O(1) or
    /// allocation-free claim is made here.
    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_hba(&self) -> Result<u32, DescriptorReadError> {
        // BEGIN RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumHBA
        // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
        // RDKit✔️✔️:   const std::string nm##Version = vers;            \
        // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
        // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
        // RDKit✔️✔️:     return m.get().countMatches(mol);              \
        // RDKit✔️✔️:   }                                                \
        // RDKit✔️✔️:   extern int no_such_variable
        // RDKit✔️✔️: SMARTSCOUNTFUNC(NumHBA,
        // RDKit✔️✔️:                 "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$("
        // RDKit✔️✔️:                 "[N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]",
        // RDKit✔️✔️:                 "2.0.2");
        // END RDKIT CPP MACRO EXPANSION: SMARTSCOUNTFUNC -> calcNumHBA
        let assignment = required_descriptor_valence(self)?;
        let rings = self
            .derived_cache_runtime()
            .valid_ring_info()
            .ok_or(DescriptorReadError::MissingInitializedRings)?;
        let input = cosmolkit_descriptors::DescriptorInput::new(
            self.topology(),
            self.coordinate_block_runtime(),
            self.properties(),
            assignment,
            rings,
        );
        cosmolkit_descriptors::num_hba_prepared(&input)
            .map_err(|source| DescriptorReadError::Algorithm { source })
    }

    /// Returns the general SMARTS-based hydrogen-bond donor count
    /// (RDKit `CalcNumHBD`, pinned `2.0.1` pattern).
    ///
    /// Mandatory-borrow review: this prepared query borrows ONLY the
    /// EXISTING validity-checked cached valence assignment through
    /// [`required_descriptor_valence`] and delegates exactly once to the
    /// narrow domain owner `num_hbd_with_valence`. The fixed pattern
    /// reads hydrogen-count and valence predicates but NO ring predicate
    /// (ROOT-resolved from the pinned source), so this query performs NO
    /// ring read, gate, find, install or upgrade, clones no whole cache
    /// state, exposes no fresh-state wrapper and declares no operation.
    /// A missing prepared valence reports
    /// [`DescriptorReadError::MissingPreparedValence`] with no fabricated
    /// cause; a domain failure maps to
    /// [`DescriptorReadError::Algorithm`] with the borrowed typed source
    /// chain. Valid valence with absent/reset rings SUCCEEDS — ring state
    /// is not an input of this query. An unsanitized molecule fails at the
    /// valence gate before any chemistry runs: a public prepared-state
    /// boundary limitation, not a source chemical divergence.
    ///
    /// Behavior review: this is NOT the direct Lipinski donor-hydrogen sum
    /// ([`Molecule::lipinski_hbd`]); the pinned pattern counts MATCHING
    /// ATOMS (sulfur donors count here, `[NH4+]` counts one atom here but
    /// four attached hydrogens there, isolated S has two H and counts 0).
    /// Cost review: borrowing the assignment does NOT make the call free —
    /// the domain owner validates topology and valence rows, builds the
    /// narrow query-match context (adjacency validation/scratch) and runs
    /// one matcher pass; no O(1) or allocation-free claim is made here.
    #[cfg(feature = "cap-descriptors")]
    #[must_use]
    pub fn num_hbd(&self) -> Result<u32, DescriptorReadError> {
        let assignment = required_descriptor_valence(self)?;
        cosmolkit_descriptors::num_hbd_with_valence(self.topology(), assignment)
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

    #[test]
    fn descriptor_heteroatoms_public_storage_internal_proof() {
        // HETERO-REPAIR frozen INTERNAL storage proof for the
        // topology-only num_heteroatoms query: 3 fixed molecules
        // (CCO, *, [2H]O[2H], literal output 1 each) x 2 receivers
        // (original + Arc-sharing peer clone) x 2 repeats = exactly 12
        // ACTUAL calls. Relocated from the superseded external
        // value-equality proof: Molecule::PartialEq EXCLUDES the derived
        // cache and proves neither Arc identity nor generic
        // NaN/signed-zero equality — this private-module proof checks
        // the four real Arc blocks, whole cloned values, valid_states()
        // and the optional borrowed payloads directly. No production
        // cache clone/access change is involved.
        const LITERALS: [(&str, u32); 3] = [("CCO", 1), ("*", 1), ("[2H]O[2H]", 1)];

        let mut calls = 0usize;
        for (smiles, literal) in LITERALS {
            let original = Molecule::from_smiles(smiles).unwrap();
            let peer = original.clone();

            // Captured ONCE before any query; never refreshed.
            let topology_arc = original.topology_arc_runtime();
            let coordinates_arc = original.coordinates_arc_runtime();
            let properties_arc = original.properties_arc_runtime();
            let cache_arc = original.derived_cache_arc_runtime();
            // Complete independently cloned baselines for whole values.
            let baseline_topology = original.topology().clone();
            let baseline_coordinates = original.coordinate_block_runtime().clone();
            let baseline_properties = original.properties().clone();
            let baseline_cache = original.derived_cache_runtime().clone();
            let baseline_states = baseline_cache.valid_states();
            // Optional borrowed payload identities (None stays None, Some
            // stays pointer-identical); whole-value cache state is
            // already covered by the complete cloned baseline above.
            let baseline_valence: Option<*const cosmolkit_core::ValenceAssignment> = original
                .derived_cache_runtime()
                .valence_assignment()
                .map(|reference| reference as *const _);
            let baseline_ring: Option<*const cosmolkit_core::RingInfo> = original
                .derived_cache_runtime()
                .ring_info()
                .map(|reference| reference as *const _);
            let baseline_valid_ring: Option<*const cosmolkit_core::RingInfo> = original
                .derived_cache_runtime()
                .valid_ring_info()
                .map(|reference| reference as *const _);

            // ONE checkpoint closure: BOTH receivers keep every
            // never-refreshed baseline (four Arc identities, four whole
            // values, valid states, optional payload identities) and
            // share all four Arc blocks with each other. Invoked
            // immediately BEFORE and AFTER each of the twelve actual
            // query calls.
            let checkpoint = |stage: &str| {
                for (who, molecule) in [("original", &original), ("peer", &peer)] {
                    assert!(
                        Arc::ptr_eq(&topology_arc, &molecule.topology_arc_runtime()),
                        "{stage} {who} {smiles}: topology Arc identity"
                    );
                    assert!(
                        Arc::ptr_eq(&coordinates_arc, &molecule.coordinates_arc_runtime()),
                        "{stage} {who} {smiles}: coordinates Arc identity"
                    );
                    assert!(
                        Arc::ptr_eq(&properties_arc, &molecule.properties_arc_runtime()),
                        "{stage} {who} {smiles}: properties Arc identity"
                    );
                    assert!(
                        Arc::ptr_eq(&cache_arc, &molecule.derived_cache_arc_runtime()),
                        "{stage} {who} {smiles}: derived-cache Arc identity"
                    );
                    assert_eq!(
                        molecule.topology(),
                        &baseline_topology,
                        "{stage} {who} {smiles}: whole topology"
                    );
                    assert_eq!(
                        molecule.coordinate_block_runtime(),
                        &baseline_coordinates,
                        "{stage} {who} {smiles}: whole coordinates"
                    );
                    assert_eq!(
                        molecule.properties(),
                        &baseline_properties,
                        "{stage} {who} {smiles}: whole properties"
                    );
                    assert_eq!(
                        molecule.derived_cache_runtime(),
                        &baseline_cache,
                        "{stage} {who} {smiles}: whole derived cache"
                    );
                    assert_eq!(
                        molecule.derived_cache_runtime().valid_states(),
                        baseline_states,
                        "{stage} {who} {smiles}: valid states"
                    );
                    let valence: Option<*const cosmolkit_core::ValenceAssignment> = molecule
                        .derived_cache_runtime()
                        .valence_assignment()
                        .map(|reference| reference as *const _);
                    assert_eq!(
                        valence, baseline_valence,
                        "{stage} {who} {smiles}: valence payload identity (None stays None)"
                    );
                    let ring: Option<*const cosmolkit_core::RingInfo> = molecule
                        .derived_cache_runtime()
                        .ring_info()
                        .map(|reference| reference as *const _);
                    assert_eq!(
                        ring, baseline_ring,
                        "{stage} {who} {smiles}: stored ring payload identity"
                    );
                    let valid_ring: Option<*const cosmolkit_core::RingInfo> = molecule
                        .derived_cache_runtime()
                        .valid_ring_info()
                        .map(|reference| reference as *const _);
                    assert_eq!(
                        valid_ring, baseline_valid_ring,
                        "{stage} {who} {smiles}: valid-ring payload identity"
                    );
                }
                assert!(
                    Arc::ptr_eq(
                        &original.topology_arc_runtime(),
                        &peer.topology_arc_runtime()
                    ),
                    "{stage} {smiles}: original/peer topology sharing"
                );
                assert!(
                    Arc::ptr_eq(
                        &original.coordinates_arc_runtime(),
                        &peer.coordinates_arc_runtime()
                    ),
                    "{stage} {smiles}: original/peer coordinates sharing"
                );
                assert!(
                    Arc::ptr_eq(
                        &original.properties_arc_runtime(),
                        &peer.properties_arc_runtime()
                    ),
                    "{stage} {smiles}: original/peer properties sharing"
                );
                assert!(
                    Arc::ptr_eq(
                        &original.derived_cache_arc_runtime(),
                        &peer.derived_cache_arc_runtime()
                    ),
                    "{stage} {smiles}: original/peer derived-cache sharing"
                );
            };

            for (who, receiver) in [("original", &original), ("peer", &peer)] {
                for repeat in 0..2 {
                    checkpoint("before");
                    let count = receiver.num_heteroatoms().unwrap();
                    calls += 1;
                    checkpoint("after");
                    assert_eq!(count, literal, "{smiles} {who} #{repeat}: literal output");
                }
            }
        }
        assert_eq!(calls, 12, "exact 12-call storage census");
    }

    #[test]
    fn descriptor_hba_storage_prepared_payloads_preserved() {
        // HBA-PUBLIC frozen internal storage proof for the prepared
        // general num_hba query: 3 fixed molecules (empty, CCO, thiophene
        // c1ccsc1; literal outputs 0/1/1) x 2 constructor policies
        // (sanitize=true, remove-H false/true) x 2 receivers (original +
        // Arc-sharing peer clone) x 2 repeats = exactly 24 ACTUAL calls.
        // Once-captured never-refreshed baselines: all four Arc pointers,
        // complete independently cloned topology/property/derived-cache
        // values, FLOAT-BIT coordinate rows (PartialEq alone proves
        // neither NaN nor signed-zero identity), valid_states(), and the
        // optional borrowed valence/stored-ring/valid-ring payload
        // identities (None stays None, Some stays pointer-identical). ONE
        // closure runs immediately BEFORE and AFTER every call over BOTH
        // receivers. No production cache clone/access change.
        const FIXTURES: [(&str, u32); 3] = [("", 0), ("CCO", 1), ("c1ccsc1", 1)];

        // Float-bit coordinate identity: collects every 2D/3D row's f64
        // bits in order. An absent block is the empty vector; a present
        // NaN/signed-zero row keeps its exact bits.
        fn coordinate_bits(block: &cosmolkit_model::CoordinateBlock) -> Vec<u64> {
            let mut bits = Vec::new();
            for conformer in &block.conformers_2d {
                for row in conformer.coordinates() {
                    bits.extend(row.iter().map(|value| value.to_bits()));
                }
            }
            for conformer in &block.conformers_3d {
                for row in conformer.coordinates() {
                    bits.extend(row.iter().map(|value| value.to_bits()));
                }
            }
            bits
        }

        let mut calls = 0usize;
        for (smiles, literal) in FIXTURES {
            for remove_hydrogens in [false, true] {
                let label = format!("{smiles:?}/rh={remove_hydrogens}");
                let original = Molecule::from_smiles_with_params(
                    smiles,
                    &crate::SmilesParseParams {
                        sanitize: true,
                        remove_hydrogens,
                        ..Default::default()
                    },
                )
                .unwrap();
                let peer = original.clone();

                // Valid prepared payloads are asserted BEFORE any query:
                // BOTH borrowed rows really exist on this real fixture.
                assert!(
                    original
                        .derived_cache_runtime()
                        .valence_assignment()
                        .is_some(),
                    "{label}: prepared valence actually valid"
                );
                assert!(
                    original.derived_cache_runtime().valid_ring_info().is_some(),
                    "{label}: ordinary rings actually initialized"
                );

                // Captured ONCE; never refreshed.
                let topology_arc = original.topology_arc_runtime();
                let coordinates_arc = original.coordinates_arc_runtime();
                let properties_arc = original.properties_arc_runtime();
                let cache_arc = original.derived_cache_arc_runtime();
                let baseline_topology = original.topology().clone();
                let baseline_coordinates = original.coordinate_block_runtime().clone();
                let baseline_coordinate_bits = coordinate_bits(&baseline_coordinates);
                let baseline_properties = original.properties().clone();
                let baseline_cache = original.derived_cache_runtime().clone();
                let baseline_states = baseline_cache.valid_states();
                let baseline_valence: Option<*const cosmolkit_core::ValenceAssignment> = original
                    .derived_cache_runtime()
                    .valence_assignment()
                    .map(|reference| reference as *const _);
                let baseline_ring: Option<*const cosmolkit_core::RingInfo> = original
                    .derived_cache_runtime()
                    .ring_info()
                    .map(|reference| reference as *const _);
                let baseline_valid_ring: Option<*const cosmolkit_core::RingInfo> = original
                    .derived_cache_runtime()
                    .valid_ring_info()
                    .map(|reference| reference as *const _);

                let checkpoint = |stage: &str| {
                    for (who, molecule) in [("original", &original), ("peer", &peer)] {
                        assert!(
                            Arc::ptr_eq(&topology_arc, &molecule.topology_arc_runtime()),
                            "{stage} {who} {label}: topology Arc identity"
                        );
                        assert!(
                            Arc::ptr_eq(&coordinates_arc, &molecule.coordinates_arc_runtime()),
                            "{stage} {who} {label}: coordinates Arc identity"
                        );
                        assert!(
                            Arc::ptr_eq(&properties_arc, &molecule.properties_arc_runtime()),
                            "{stage} {who} {label}: properties Arc identity"
                        );
                        assert!(
                            Arc::ptr_eq(&cache_arc, &molecule.derived_cache_arc_runtime()),
                            "{stage} {who} {label}: derived-cache Arc identity"
                        );
                        assert_eq!(
                            molecule.topology(),
                            &baseline_topology,
                            "{stage} {who} {label}: whole topology"
                        );
                        assert_eq!(
                            coordinate_bits(molecule.coordinate_block_runtime()),
                            baseline_coordinate_bits,
                            "{stage} {who} {label}: float-bit coordinates"
                        );
                        // HBA-PRESERVE: the COMPLETE coordinate block is
                        // asserted by whole-value equality against the
                        // never-refreshed baseline clone (conformer IDs,
                        // metadata, dimension bookkeeping and all block
                        // fields — not only the flat f64 bits above). The
                        // bit check stays: it pins the exact float
                        // representation where rows exist; this whole-block
                        // check pins every non-row field. Current fixtures
                        // carry EMPTY coordinate blocks, so today the
                        // nonempty-row portion of the bit proof is
                        // exercised by no case — documented, not claimed.
                        assert_eq!(
                            molecule.coordinate_block_runtime(),
                            &baseline_coordinates,
                            "{stage} {who} {label}: whole coordinate block"
                        );
                        assert_eq!(
                            molecule.properties(),
                            &baseline_properties,
                            "{stage} {who} {label}: whole properties"
                        );
                        assert_eq!(
                            molecule.derived_cache_runtime(),
                            &baseline_cache,
                            "{stage} {who} {label}: whole derived cache"
                        );
                        assert_eq!(
                            molecule.derived_cache_runtime().valid_states(),
                            baseline_states,
                            "{stage} {who} {label}: valid states"
                        );
                        let valence: Option<*const cosmolkit_core::ValenceAssignment> = molecule
                            .derived_cache_runtime()
                            .valence_assignment()
                            .map(|reference| reference as *const _);
                        assert_eq!(
                            valence, baseline_valence,
                            "{stage} {who} {label}: valence payload identity"
                        );
                        let ring: Option<*const cosmolkit_core::RingInfo> = molecule
                            .derived_cache_runtime()
                            .ring_info()
                            .map(|reference| reference as *const _);
                        assert_eq!(
                            ring, baseline_ring,
                            "{stage} {who} {label}: stored ring payload identity"
                        );
                        let valid_ring: Option<*const cosmolkit_core::RingInfo> = molecule
                            .derived_cache_runtime()
                            .valid_ring_info()
                            .map(|reference| reference as *const _);
                        assert_eq!(
                            valid_ring, baseline_valid_ring,
                            "{stage} {who} {label}: valid-ring payload identity"
                        );
                    }
                    assert!(
                        Arc::ptr_eq(
                            &original.topology_arc_runtime(),
                            &peer.topology_arc_runtime()
                        ),
                        "{stage} {label}: original/peer topology sharing"
                    );
                    assert!(
                        Arc::ptr_eq(
                            &original.coordinates_arc_runtime(),
                            &peer.coordinates_arc_runtime()
                        ),
                        "{stage} {label}: original/peer coordinates sharing"
                    );
                    assert!(
                        Arc::ptr_eq(
                            &original.properties_arc_runtime(),
                            &peer.properties_arc_runtime()
                        ),
                        "{stage} {label}: original/peer properties sharing"
                    );
                    assert!(
                        Arc::ptr_eq(
                            &original.derived_cache_arc_runtime(),
                            &peer.derived_cache_arc_runtime()
                        ),
                        "{stage} {label}: original/peer derived-cache sharing"
                    );
                };

                for (who, receiver) in [("original", &original), ("peer", &peer)] {
                    for repeat in 0..2 {
                        checkpoint("before");
                        let count = receiver.num_hba().unwrap();
                        calls += 1;
                        checkpoint("after");
                        assert_eq!(count, literal, "{label} {who} #{repeat}: literal output");
                    }
                }
            }
        }
        assert_eq!(calls, 24, "exact 24-call storage census");
    }

    #[test]
    fn descriptor_hba_missing_initialized_rings_fixture() {
        // Real private fixture: the existing owning constructor seam
        // transports a REAL parsed CCO with its REAL valid valence
        // assignment but an ABSENT ordinary ring carrier (None clears
        // storage and validity). num_hba must then report the typed
        // MissingInitializedRings — a real state transition observed
        // through the public query, not a manually constructed error.
        let base = Molecule::from_smiles("CCO").unwrap();
        let assignment = base
            .derived_cache_runtime()
            .valence_assignment()
            .expect("sanitized base carries a valid assignment")
            .clone();
        let molecule = Molecule::from_smiles_parts_with_derived_state(
            base.topology().clone(),
            base.coordinate_block_runtime().clone(),
            base.properties().clone(),
            Some(assignment),
            None,
        )
        .unwrap();
        // The real reached state: valence valid, ordinary rings absent.
        assert!(
            molecule
                .derived_cache_runtime()
                .valence_assignment()
                .is_some(),
            "fixture keeps the valid valence"
        );
        assert!(
            molecule.derived_cache_runtime().valid_ring_info().is_none(),
            "fixture has no valid ordinary rings"
        );
        let observer = molecule.clone();
        let error = molecule.num_hba().unwrap_err();
        assert!(
            matches!(error, DescriptorReadError::MissingInitializedRings),
            "got {error:?}"
        );
        assert_eq!(molecule, observer, "failed read preserves the input");

        // Absence precedence on a molecule with NEITHER payload: the
        // valence gate reports first.
        let raw_params = crate::SmilesParseParams {
            sanitize: false,
            ..Default::default()
        };
        let raw = Molecule::from_smiles_with_params("CCO", &raw_params).unwrap();
        assert!(matches!(
            raw.num_hba().unwrap_err(),
            DescriptorReadError::MissingPreparedValence
        ));
    }

    #[test]
    fn descriptor_hbd_storage_prepared_payloads_preserved() {
        // HBD-PUBLIC frozen internal storage proof for the narrow
        // valence-only prepared query: 3 fixed molecules (empty, CCO,
        // pyrrole c1cc[nH]c1; literal outputs 0/1/1) x 2 constructor
        // policies (sanitize=true, remove-H false/true) x 2 receivers
        // (original + Arc-sharing peer clone) x 2 repeats = exactly 24
        // ACTUAL calls. Once-captured never-refreshed baselines: all four
        // Arc pointers, complete independently cloned topology/property/
        // derived-cache values, the COMPLETE coordinate block by
        // whole-value equality PLUS its float-bit rows/conformer IDs where
        // present, valid_states(), and the optional borrowed valence/
        // stored-ring/valid-ring payload identities (None stays None, Some
        // stays pointer-identical). These constructor fixtures have EMPTY
        // coordinate blocks — a documented limit, not NaN coverage. ONE
        // checkpoint runs immediately BEFORE and AFTER every call over
        // BOTH receivers. No production cache clone/access change.
        const FIXTURES: [(&str, u32); 3] = [("", 0), ("CCO", 1), ("c1cc[nH]c1", 1)];

        fn coordinate_bits_and_ids(block: &cosmolkit_model::CoordinateBlock) -> Vec<u64> {
            let mut bits = Vec::new();
            for conformer in &block.conformers_2d {
                bits.push(u64::try_from(conformer.id()).unwrap_or(u64::MAX));
                for row in conformer.coordinates() {
                    bits.extend(row.iter().map(|value| value.to_bits()));
                }
            }
            for conformer in &block.conformers_3d {
                bits.push(u64::try_from(conformer.id()).unwrap_or(u64::MAX));
                for row in conformer.coordinates() {
                    bits.extend(row.iter().map(|value| value.to_bits()));
                }
            }
            bits
        }

        let mut calls = 0usize;
        for (smiles, literal) in FIXTURES {
            for remove_hydrogens in [false, true] {
                let label = format!("{smiles:?}/rh={remove_hydrogens}");
                let original = Molecule::from_smiles_with_params(
                    smiles,
                    &crate::SmilesParseParams {
                        sanitize: true,
                        remove_hydrogens,
                        ..Default::default()
                    },
                )
                .unwrap();
                let peer = original.clone();

                // Real prepared payload prerequisite BEFORE any query: the
                // valence assignment exists; ring state may be anything.
                assert!(
                    original
                        .derived_cache_runtime()
                        .valence_assignment()
                        .is_some(),
                    "{label}: prepared valence actually valid"
                );

                let topology_arc = original.topology_arc_runtime();
                let coordinates_arc = original.coordinates_arc_runtime();
                let properties_arc = original.properties_arc_runtime();
                let cache_arc = original.derived_cache_arc_runtime();
                let baseline_topology = original.topology().clone();
                let baseline_coordinates = original.coordinate_block_runtime().clone();
                let baseline_coordinate_bits = coordinate_bits_and_ids(&baseline_coordinates);
                let baseline_properties = original.properties().clone();
                let baseline_cache = original.derived_cache_runtime().clone();
                let baseline_states = baseline_cache.valid_states();
                let baseline_valence: Option<*const cosmolkit_core::ValenceAssignment> = original
                    .derived_cache_runtime()
                    .valence_assignment()
                    .map(|reference| reference as *const _);
                let baseline_ring: Option<*const cosmolkit_core::RingInfo> = original
                    .derived_cache_runtime()
                    .ring_info()
                    .map(|reference| reference as *const _);
                let baseline_valid_ring: Option<*const cosmolkit_core::RingInfo> = original
                    .derived_cache_runtime()
                    .valid_ring_info()
                    .map(|reference| reference as *const _);

                let checkpoint = |stage: &str| {
                    for (who, molecule) in [("original", &original), ("peer", &peer)] {
                        assert!(
                            Arc::ptr_eq(&topology_arc, &molecule.topology_arc_runtime()),
                            "{stage} {who} {label}: topology Arc identity"
                        );
                        assert!(
                            Arc::ptr_eq(&coordinates_arc, &molecule.coordinates_arc_runtime()),
                            "{stage} {who} {label}: coordinates Arc identity"
                        );
                        assert!(
                            Arc::ptr_eq(&properties_arc, &molecule.properties_arc_runtime()),
                            "{stage} {who} {label}: properties Arc identity"
                        );
                        assert!(
                            Arc::ptr_eq(&cache_arc, &molecule.derived_cache_arc_runtime()),
                            "{stage} {who} {label}: derived-cache Arc identity"
                        );
                        assert_eq!(
                            molecule.topology(),
                            &baseline_topology,
                            "{stage} {who} {label}: whole topology"
                        );
                        assert_eq!(
                            molecule.coordinate_block_runtime(),
                            &baseline_coordinates,
                            "{stage} {who} {label}: whole coordinate block"
                        );
                        assert_eq!(
                            coordinate_bits_and_ids(molecule.coordinate_block_runtime()),
                            baseline_coordinate_bits,
                            "{stage} {who} {label}: coordinate component bits/IDs"
                        );
                        assert_eq!(
                            molecule.properties(),
                            &baseline_properties,
                            "{stage} {who} {label}: whole properties"
                        );
                        assert_eq!(
                            molecule.derived_cache_runtime(),
                            &baseline_cache,
                            "{stage} {who} {label}: whole derived cache"
                        );
                        assert_eq!(
                            molecule.derived_cache_runtime().valid_states(),
                            baseline_states,
                            "{stage} {who} {label}: valid states"
                        );
                        let valence: Option<*const cosmolkit_core::ValenceAssignment> = molecule
                            .derived_cache_runtime()
                            .valence_assignment()
                            .map(|reference| reference as *const _);
                        assert_eq!(
                            valence, baseline_valence,
                            "{stage} {who} {label}: valence payload identity/validity"
                        );
                        let ring: Option<*const cosmolkit_core::RingInfo> = molecule
                            .derived_cache_runtime()
                            .ring_info()
                            .map(|reference| reference as *const _);
                        assert_eq!(
                            ring, baseline_ring,
                            "{stage} {who} {label}: ring payload unchanged"
                        );
                        let valid_ring: Option<*const cosmolkit_core::RingInfo> = molecule
                            .derived_cache_runtime()
                            .valid_ring_info()
                            .map(|reference| reference as *const _);
                        assert_eq!(
                            valid_ring, baseline_valid_ring,
                            "{stage} {who} {label}: valid-ring payload unchanged"
                        );
                    }
                    assert!(
                        Arc::ptr_eq(
                            &original.topology_arc_runtime(),
                            &peer.topology_arc_runtime()
                        ),
                        "{stage} {label}: original/peer topology sharing"
                    );
                    assert!(
                        Arc::ptr_eq(
                            &original.coordinates_arc_runtime(),
                            &peer.coordinates_arc_runtime()
                        ),
                        "{stage} {label}: original/peer coordinates sharing"
                    );
                    assert!(
                        Arc::ptr_eq(
                            &original.properties_arc_runtime(),
                            &peer.properties_arc_runtime()
                        ),
                        "{stage} {label}: original/peer properties sharing"
                    );
                    assert!(
                        Arc::ptr_eq(
                            &original.derived_cache_arc_runtime(),
                            &peer.derived_cache_arc_runtime()
                        ),
                        "{stage} {label}: original/peer derived-cache sharing"
                    );
                };

                for (who, receiver) in [("original", &original), ("peer", &peer)] {
                    for repeat in 0..2 {
                        checkpoint("before");
                        let count = receiver.num_hbd().unwrap();
                        calls += 1;
                        checkpoint("after");
                        assert_eq!(count, literal, "{label} {who} #{repeat}: literal output");
                    }
                }
            }
        }
        assert_eq!(calls, 24, "exact 24-call storage census");
    }

    #[test]
    fn descriptor_hbd_ring_absence_and_hba_contrast_fixture() {
        // Real owning constructor seam: final CCO with its REAL valid
        // valence but rings None (cleared storage/validity). The narrow
        // HBD query MUST succeed (num_hbd=1 — the pattern reads no ring
        // predicate), while the SAME molecule's num_hba reports the typed
        // MissingInitializedRings. Snapshots BEFORE and AFTER both calls;
        // no new public debug API.
        let base = Molecule::from_smiles("CCO").unwrap();
        let assignment = base
            .derived_cache_runtime()
            .valence_assignment()
            .expect("sanitized base carries a valid assignment")
            .clone();
        let molecule = Molecule::from_smiles_parts_with_derived_state(
            base.topology().clone(),
            base.coordinate_block_runtime().clone(),
            base.properties().clone(),
            Some(assignment),
            None,
        )
        .unwrap();
        assert!(
            molecule
                .derived_cache_runtime()
                .valence_assignment()
                .is_some(),
            "fixture keeps the valid valence"
        );
        assert!(
            molecule.derived_cache_runtime().valid_ring_info().is_none(),
            "fixture has no valid ordinary rings"
        );
        let snapshot = molecule.clone();
        let donors = molecule.num_hbd().expect("no ring gate on num_hbd");
        assert_eq!(donors, 1, "CCO donor count with absent rings");
        assert_eq!(molecule, snapshot, "num_hbd preserves the input");
        let snapshot_after_hbd = molecule.clone();
        let error = molecule.num_hba().unwrap_err();
        assert!(
            matches!(error, DescriptorReadError::MissingInitializedRings),
            "same molecule num_hba: got {error:?}"
        );
        assert_eq!(
            molecule, snapshot_after_hbd,
            "num_hba failure preserves the input"
        );
    }

    #[test]
    fn descriptor_hbd_malformed_domain_rows_preserve_inputs() {
        // Two malformed valence lengths through the REAL domain narrow
        // call on a live molecule's borrowed rows: whole supplied values
        // are preserved and the typed source chain is borrowed, not
        // flattened.
        let molecule = Molecule::from_smiles("CCO").unwrap();
        let topology = molecule.topology().clone();
        // The validator checks explicit_valence first, so the non-target
        // field carries the valid length in each case.
        for (field, explicit, implicit) in [
            ("explicit_valence", 2usize, 3usize),
            ("implicit_hydrogens", 3, 5),
        ] {
            let malformed = cosmolkit_core::ValenceAssignment {
                explicit_valence: vec![1; explicit],
                implicit_hydrogens: vec![1; implicit],
            };
            let malformed_before = malformed.clone();
            let topology_before = topology.clone();
            let error = cosmolkit_descriptors::num_hbd_with_valence(&topology, &malformed)
                .err()
                .unwrap_or_else(|| panic!("{field} length must be rejected"));
            let cosmolkit_descriptors::DescriptorError::Search {
                source: cosmolkit_descriptors::DescriptorSearchCause::Context(context),
                ..
            } = &error
            else {
                panic!("expected Search/Context, got {error:?}")
            };
            assert!(
                std::error::Error::source(&error).is_some(),
                "{field}: borrowed source chain retained"
            );
            let cosmolkit_search::QueryMatchContextError::ValenceRows {
                field: actual_field,
                expected,
                actual,
            } = context
            else {
                panic!("expected ValenceRows, got {context:?}")
            };
            assert_eq!(*actual_field, field, "exact field");
            assert_eq!(*expected, 3, "expected = CCO atom count");
            assert_eq!(
                *actual,
                if field == "explicit_valence" {
                    explicit
                } else {
                    implicit
                },
                "exact actual length"
            );
            assert_eq!(malformed, malformed_before, "whole valence preserved");
            assert_eq!(topology, topology_before, "whole topology preserved");
        }
    }
}
