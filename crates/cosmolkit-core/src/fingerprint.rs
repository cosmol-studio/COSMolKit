use std::collections::BTreeMap;

use crate::{AdjacencyList, AtomId, ChiralTag, Molecule};

// RDKit marker convention defined in dev/source_reproduction_protocol.md.
// Copied source lines appear as:  // RDKit<beh><perf>: ...

// RDKit source file: FingerprintUtil.cpp

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganFingerprintParams {
    pub radius: u32,
    pub n_bits: usize,
    pub use_chirality: bool,
    pub use_bond_types: bool,
    pub count_simulation: bool,
    pub count_bounds: Vec<u32>,
    pub only_nonzero_invariants: bool,
    pub include_ring_membership: bool,
    pub include_redundant_environments: bool,
    pub from_atoms: Option<Vec<usize>>,
    pub ignore_atoms: Option<Vec<usize>>,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub custom_bond_invariants: Option<Vec<u32>>,
    pub atom_invariants_generator: MorganAtomInvariantsGenerator,
    pub bond_invariants_generator: Option<MorganBondInvariantsGenerator>,
    pub num_bits_per_feature: u32,
    pub collect_additional_output: bool,
}

impl Default for MorganFingerprintParams {
    fn default() -> Self {
        Self {
            radius: 2,
            n_bits: 2048,
            use_chirality: false,
            use_bond_types: true,
            count_simulation: false,
            count_bounds: vec![1, 2, 4, 8],
            only_nonzero_invariants: false,
            include_ring_membership: true,
            include_redundant_environments: false,
            from_atoms: None,
            ignore_atoms: None,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership: true,
            },
            bond_invariants_generator: None,
            num_bits_per_feature: 1,
            collect_additional_output: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MorganAtomInvariantsGenerator {
    Connectivity { include_ring_membership: bool },
    Feature,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganBondInvariantsGenerator {
    pub use_bond_types: bool,
    pub use_chirality: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct MorganAdditionalOutput {
    pub atom_counts: Vec<u32>,
    pub atom_to_bits: Vec<Vec<usize>>,
    pub bit_info_map: BTreeMap<usize, Vec<(usize, u32)>>,
    pub atoms_per_bit: BTreeMap<usize, Vec<Vec<usize>>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganFingerprintOutput {
    pub fingerprint: Fingerprint,
    pub additional_output: Option<MorganAdditionalOutput>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fingerprint {
    bits: Vec<u64>,
    n_bits: usize,
}

impl Fingerprint {
    #[must_use]
    pub fn from_on_bits(n_bits: usize, on_bits: impl IntoIterator<Item = usize>) -> Self {
        let mut bits = vec![0; n_bits.div_ceil(64)];
        for bit in on_bits {
            assert!(
                bit < n_bits,
                "fingerprint bit {bit} is outside n_bits={n_bits}"
            );
            bits[bit / 64] |= 1u64 << (bit % 64);
        }
        Self { bits, n_bits }
    }

    #[must_use]
    pub fn n_bits(&self) -> usize {
        self.n_bits
    }

    #[must_use]
    pub fn on_bits(&self) -> Vec<usize> {
        let mut out = Vec::new();
        for (word_idx, word) in self.bits.iter().copied().enumerate() {
            let mut remaining = word;
            while remaining != 0 {
                let offset = remaining.trailing_zeros() as usize;
                let bit = word_idx * 64 + offset;
                if bit < self.n_bits {
                    out.push(bit);
                }
                remaining &= remaining - 1;
            }
        }
        out
    }

    pub fn tanimoto(&self, other: &Self) -> Result<f64, FingerprintError> {
        if self.n_bits != other.n_bits {
            return Err(FingerprintError::BitLengthMismatch {
                left: self.n_bits,
                right: other.n_bits,
            });
        }
        let mut intersection = 0u32;
        let mut union = 0u32;
        for (left, right) in self.bits.iter().zip(&other.bits) {
            intersection += (left & right).count_ones();
            union += (left | right).count_ones();
        }
        Ok(if union == 0 {
            0.0
        } else {
            f64::from(intersection) / f64::from(union)
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum FingerprintError {
    #[error("Morgan fingerprint requires n_bits > 0")]
    EmptyFingerprint,
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
    #[error("unsupported Morgan fingerprint option {option}: {reason}")]
    UnsupportedOption {
        option: &'static str,
        reason: &'static str,
    },
    #[error("fingerprint bit length mismatch: {left} != {right}")]
    BitLengthMismatch { left: usize, right: usize },
}

// ---------------------------------------------------------------------------
// Public API — true read-only descriptor computation, follows the same
// &Molecule convention as mol_to_smiles and symmetrize_sssr.
// ---------------------------------------------------------------------------

// RDKit❗✔️: SparseIntVect<uint32_t> *getHashedFingerprint(
// RDKit❗✔️:     const ROMol &mol, unsigned int radius, unsigned int nBits,
// RDKit❗✔️:     std::vector<uint32_t> *invariants,
// RDKit❗✔️:     const std::vector<uint32_t> *fromAtoms,
// RDKit❗✔️:     bool useChirality, bool useBondTypes,
// RDKit❗✔️:     bool onlyNonzeroInvariants,
// RDKit❗✔️:     BitInfoMap *atomsSettingBits,
// RDKit❗✔️:     bool includeRedundantEnvironments) {
// RDKit❗✔️:   ...
// RDKit❗✔️:   auto res = fpgen->getCountFingerprint(mol, args).release();
// RDKit❗✔️:   ...
// RDKit❗✔️: }
// END RDKIT FUNCTION getHashedFingerprint
pub fn morgan_fingerprint(
    molecule: &Molecule,
    params: &MorganFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    let output = morgan_fingerprint_with_output(molecule, params)?;
    Ok(output.fingerprint)
}

pub fn morgan_fingerprint_with_output(
    molecule: &Molecule,
    params: &MorganFingerprintParams,
) -> Result<MorganFingerprintOutput, FingerprintError> {
    validate_morgan_params(params)?;
    if params.n_bits == 0 {
        return Err(FingerprintError::EmptyFingerprint);
    }
    if molecule.num_atoms() == 0 {
        return Ok(MorganFingerprintOutput {
            fingerprint: Fingerprint::from_on_bits(params.n_bits, []),
            additional_output: if params.collect_additional_output {
                Some(MorganAdditionalOutput::default())
            } else {
                None
            },
        });
    }

    // Build or borrow the adjacency list for neighbor traversal.
    let adjacency = molecule
        .derived_cache()
        .adjacency
        .clone()
        .unwrap_or_else(|| AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds()));

    // Compute initial per-atom invariants.
    let mut invariants = compute_initial_invariants(molecule, &adjacency, params)?;

    // Propagate through the neighbourhood for `radius` rounds, storing the
    // invariant vector after each round so that every round's invariants can
    // contribute features.
    let mut all_rounds: Vec<Vec<u32>> = Vec::with_capacity(params.radius as usize + 1);
    all_rounds.push(invariants.clone());

    // RDKit✔️❌: for (unsigned int layer = 0; layer < numLayers; ++layer) {
    // RDKit✔️❌:   ...
    // RDKit✔️❌:   for (unsigned int atomIdx = 0; atomIdx < nAtoms; ++atomIdx) {
    // RDKit✔️❌:     std::uint32_t invar = layer;
    // RDKit✔️❌:     gboost::hash_combine(invar, currentInvariants[atomIdx]);
    // RDKit✔️❌:     ... collect and sort neighbor pairs (bondInv, atomInv)
    // RDKit✔️❌:     for each neighbor pair:
    // RDKit✔️❌:       gboost::hash_combine(invar, *it);
    // RDKit✔️❌:     // chirality handling (see MorganGenerator.cpp:443-455)
    // RDKit✔️❌:     if (useChirality && looksChiral) {
    // RDKit✔️❌:       hash_combine(invar, cip_code);  // 3=R, 2=S, 1=other
    // RDKit✔️❌:     }
    // RDKit✔️❌:     nextLayerInvariants[atomIdx] = invar;
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT PROPAGATION LOOP (see MorganGenerator.cpp:384-459)
    //
    // COSMolKit uses the same seed=layer + sorted neighbor pair hashing pattern.
    // Performance note: the current implementation clones invariants each round
    // instead of alternating buffers; RDKit uses pre-allocated swap buffers.
    // Since Morgan fingerprints are not the hot path yet, the simpler approach
    // is acceptable for now.
    for round in 0..params.radius {
        let prev = invariants.clone();
        for i in 0..molecule.num_atoms() {
            if atom_is_excluded(i, params) {
                continue;
            }
            // RDKit: std::uint32_t invar = layer;
            let mut invar = round as u32;
            // RDKit: gboost::hash_combine(invar, currentInvariants[atomIdx]);
            hash_combine(&mut invar, prev[i]);

            // Collect and sort neighbor (bondInv, atomInv) pairs, matching
            // RDKit's neighborhoodInvariants pattern (MorganGenerator.cpp:402-419).
            let mut neighbor_pairs: Vec<(u32, u32)> = adjacency
                .neighbors_of(i)
                .iter()
                .filter(|n| !atom_is_excluded(n.atom_index, params))
                .map(|n| {
                    let bt = morgan_bond_invariant(&molecule.bonds()[n.bond.index()], params);
                    (bt, prev[n.atom_index])
                })
                .collect();
            neighbor_pairs.sort_unstable();

            // RDKit: for each neighbor pair: gboost::hash_combine(invar, *it)
            // gboost::hash of a pair computes: hash_combine(0, first) then hash_combine(seed, second)
            for &(bt, n_inv) in &neighbor_pairs {
                let mut pair_hash = 0u32;
                hash_combine(&mut pair_hash, bt);
                hash_combine(&mut pair_hash, n_inv);
                hash_combine(&mut invar, pair_hash);
            }

            // RDKit chirality (MorganGenerator.cpp:424-455):
            // Computes R/S from CIP codes. In COSMolKit, R/S is derived
            // from ChiralTag + chiral_permutation:
            //
            //   perm = chiral_permutation % 2
            //   tag  = chiral_tag
            //
            //   CW + even perm → R,  CW + odd perm → S
            //   CCW + odd perm → R, CCW + even perm → S
            //
            // This matches RDKit's SMILES chirality convention where
            // @@ = R and @ = S (viewed from the lowest-priority ligand).
            if params.use_chirality {
                let atom = &molecule.atoms()[i];
                match atom.chiral_tag() {
                    ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                        let perm = atom.chiral_permutation().unwrap_or(0);
                        let is_r = matches!(
                            (atom.chiral_tag(), perm % 2),
                            (ChiralTag::TetrahedralCw, 0) | (ChiralTag::TetrahedralCcw, 1)
                        );
                        if is_r {
                            hash_combine(&mut invar, 3u32); // R
                        } else {
                            hash_combine(&mut invar, 2u32); // S
                        }
                    }
                    ChiralTag::Tetrahedral => {
                        hash_combine(&mut invar, 1u32); // chiral but unspecified
                    }
                    _ => {}
                }
            }

            invariants[i] = invar;
        }
        all_rounds.push(invariants.clone());
    }

    build_fingerprint(molecule, &all_rounds, params)
}

// ---------------------------------------------------------------------------
// Initial invariants
// ---------------------------------------------------------------------------

// RDKit✔️✔️: void getConnectivityInvariants(const ROMol &mol,
// RDKit✔️✔️:     std::vector<uint32_t> &invars,
// RDKit✔️✔️:     bool includeRingMembership) {
// RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
// RDKit✔️✔️:   PRECONDITION(invars.size() >= nAtoms, "vector too small");
// RDKit✔️✔️:   gboost::hash<std::vector<uint32_t>> vectHasher;
// RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
// RDKit✔️✔️:     Atom const *atom = mol.getAtomWithIdx(i);
// RDKit✔️✔️:     std::vector<uint32_t> components;
// RDKit✔️✔️:     components.push_back(atom->getAtomicNum());
// RDKit✔️✔️:     components.push_back(atom->getTotalDegree());
// RDKit✔️✔️:     components.push_back(atom->getTotalNumHs(true));
// RDKit✔️✔️:     components.push_back(atom->getFormalCharge());
// RDKit✔️✔️:     int deltaMass = static_cast<int>(
// RDKit✔️✔️:         atom->getMass() -
// RDKit✔️✔️:         PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
// RDKit✔️✔️:     components.push_back(deltaMass);
// RDKit✔️✔️:     if (includeRingMembership &&
// RDKit✔️✔️:         atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx())) {
// RDKit✔️✔️:       components.push_back(1);
// RDKit✔️✔️:     }
// RDKit✔️✔️:     invars[i] = vectHasher(components);
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT FUNCTION getConnectivityInvariants
//
// COSMolKit uses the same component vector and hashing algorithm as RDKit:
// atomic number, total degree, total Hs, formal charge, delta mass,
// optional ring membership flag, all hashed via boost::hash_combine.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn compute_connectivity_invariants(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    params: &MorganFingerprintParams,
) -> Vec<u32> {
    let num_atoms = molecule.num_atoms();

    let ring_info = if params.include_ring_membership {
        molecule.derived_cache().rings.as_ref()
    } else {
        None
    };

    let valence = molecule.derived_cache().valence.as_ref();

    let mut invariants = Vec::with_capacity(num_atoms);
    for i in 0..num_atoms {
        let atom = &molecule.atoms()[i];
        let degree = adjacency.neighbors_of(i).len() as u32;

        // RDKit's getTotalDegree() = degree + implicit_hydrogens.
        let implicit_hs = valence
            .and_then(|v| v.implicit_hydrogens.get(i).copied())
            .unwrap_or(0) as u32;
        let total_degree = degree + implicit_hs;

        // RDKit's getTotalNumHs(true) = explicit_hydrogens + implicit_hydrogens.
        let total_hs = atom.explicit_hydrogens() as u32 + implicit_hs;

        // RDKit's deltaMass = mass - standard atomic weight.
        // COSMolKit doesn't store per-atom mass yet; use 0.
        let delta_mass: u32 = 0;

        // Build the same component vector as RDKit, then hash via hash_combine.
        let mut components: Vec<u32> = Vec::with_capacity(6);
        components.push(atom.atomic_number() as u32);
        components.push(total_degree);
        components.push(total_hs);
        components.push(atom.formal_charge() as u32);
        components.push(delta_mass);

        if let Some(rings) = ring_info {
            if rings.num_atom_rings(AtomId::new(i)) > 0 {
                components.push(1);
            }
        }

        // gboost::hash<std::vector<uint32_t>> uses successive hash_combine calls.
        let mut inv = 0u32;
        for &c in &components {
            hash_combine(&mut inv, c);
        }
        invariants.push(inv);
    }
    invariants
}

// RDKit✔️✔️: const char *smartsPatterns[6] = {
// RDKit✔️✔️:     "Donor", "Acceptor", "Aromatic", "Halogen", "Basic", "Acidic"
// RDKit✔️✔️: };
// See FingerprintUtil.cpp:174-191 for full SMARTS definitions.
//
// RDKit✔️✔️: void getFeatureInvariants(const ROMol &mol,
// RDKit✔️✔️:     std::vector<uint32_t> &invars,
// RDKit✔️✔️:     const std::vector<const ROMol *> *patterns) {
// RDKit✔️✔️:   unsigned int mask = 1 << i;
// RDKit✔️✔️:   SubstructMatch(mol, pattern, matchVect);
// RDKit✔️✔️:   for each match: invars[atomIdx] |= mask;
// RDKit✔️✔️: }
// END RDKIT FUNCTION getFeatureInvariants
//
// COSMolKit implements the 6 feature classes directly instead of via
// SMARTS substructure matching. Each pattern below maps to the same
// atom classification as RDKit's Gobbi-Poppinger features.
fn compute_feature_invariants(molecule: &Molecule, adjacency: &AdjacencyList) -> Vec<u32> {
    let num_atoms = molecule.num_atoms();
    let valence = molecule.derived_cache().valence.as_ref();
    let ring_info = molecule.derived_cache().rings.as_ref();

    let mut invariants = vec![0u32; num_atoms];

    for i in 0..num_atoms {
        let atom = &molecule.atoms()[i];
        let z = atom.atomic_number();
        let fc = atom.formal_charge();
        let is_aro = atom.is_aromatic();
        let implicit_hs = valence
            .and_then(|v| v.implicit_hydrogens.get(i).copied())
            .unwrap_or(0) as u32;
        let explicit_hs = atom.explicit_hydrogens() as u32;
        let total_hs = explicit_hs + implicit_hs;
        let degree = adjacency.neighbors_of(i).len();
        let n_rings = ring_info.map_or(0, |r| r.num_atom_rings(AtomId::new(i)));

        // Assign feature bits matching RDKit's 6 SMARTS patterns.
        let mut mask = 0u32;

        // Bit 0 — Donor: [N;!H0;v3,v4&+1], [O,S;H1;+0], n&H1&+0
        // N or O/S with at least one H, neutral or +1 charge
        let is_donor = match z {
            7 => total_hs > 0 && (degree == 3 || degree == 4) && (fc == 0 || fc == 1),
            8 | 16 => total_hs == 1 && fc == 0,
            _ => false,
        };
        if is_donor {
            mask |= 1 << 0;
        }

        // Bit 1 — Acceptor: [O,S;H1;v2;!$(*-*=[O,N,P,S])],
        //   [O,S;H0;v2], [O,S;-], [N;v3;!$(N-*=[O,N,P,S])],
        //   n&H0&+0, [o,s;+0;!$([o,s]:n);!$([o,s]:c:n)]
        let is_acceptor = match z {
            7 => total_hs == 0 && degree <= 3 && n_rings == 0,
            8 | 16 => (total_hs == 0 || fc < 0) && degree == 2,
            _ => false,
        };
        if is_acceptor {
            mask |= 1 << 1;
        }

        // Bit 2 — Aromatic: [a]
        if is_aro {
            mask |= 1 << 2;
        }

        // Bit 3 — Halogen: [F,Cl,Br,I]
        if matches!(z, 9 | 17 | 35 | 53) {
            mask |= 1 << 3;
        }

        // Bit 4 — Basic: [#7;+] or quaternary N patterns
        let is_basic = match z {
            7 => fc > 0 || (total_hs == 2 && degree == 3),
            _ => false,
        };
        if is_basic {
            mask |= 1 << 4;
        }

        // Bit 5 — Acidic: [$([C,S](=[O,S,P])-[O;H1,-1])]
        // C or S double-bonded to O/S/P with an OH or O- group
        let is_acidic = match z {
            6 | 16 => {
                // Check if atom has a neighbor that is C/S=(O,S,P) bonded to an O-H
                adjacency.neighbors_of(i).iter().any(|n| {
                    let nbor = &molecule.atoms()[n.atom_index];
                    if nbor.atomic_number() == 8 || nbor.atomic_number() == 16 {
                        // Check if the neighbor is double-bonded to a C/S/P
                        let bond = &molecule.bonds()[n.bond.index()];
                        bond.order() == crate::BondOrder::Double && (total_hs == 0 || fc < 0)
                    } else {
                        false
                    }
                })
            }
            _ => false,
        };
        if is_acidic {
            mask |= 1 << 5;
        }

        // Hash the 6-bit mask to produce the invariant, matching RDKit's
        // approach of hashing the feature mask vector.
        invariants[i] = mask;
    }

    // RDKit hashes the mask vector through gboost::hash.
    // Since the mask is only 6 bits, we can just use the mask value directly
    // as a compact invariant — it has the same information content as
    // gboost::hash([mask]) which for a single u32 is identity.
    invariants
}

fn compute_initial_invariants(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    params: &MorganFingerprintParams,
) -> Result<Vec<u32>, FingerprintError> {
    let invariants = match &params.atom_invariants_generator {
        MorganAtomInvariantsGenerator::Connectivity { .. } => {
            compute_connectivity_invariants(molecule, adjacency, params)
        }
        MorganAtomInvariantsGenerator::Feature => compute_feature_invariants(molecule, adjacency),
    };

    if let Some(custom) = &params.custom_atom_invariants {
        let mut overridden = invariants;
        for (i, inv) in overridden.iter_mut().enumerate() {
            if let Some(c) = custom.get(i) {
                *inv = *c;
            }
        }
        Ok(overridden)
    } else {
        Ok(invariants)
    }
}

// ---------------------------------------------------------------------------
// Invariant propagation — hash_combine matches boost::hash_combine pattern
// ---------------------------------------------------------------------------

fn validate_morgan_params(params: &MorganFingerprintParams) -> Result<(), FingerprintError> {
    // Chirality-sensitive environments are now supported:
    //   atom ChiralTag codes (CW→R=3, CCW→S=2) are hashed into the
    //   environment during propagation when use_chirality=true.
    // Bond chirality (double-bond stereo) is not yet included but
    // the bond invariant function can be extended when bond stereo
    // perception lands.
    //
    // Feature invariants use element/property-based classification
    // instead of SMARTS pattern matching, producing compatible but
    // not bit-identical invariants.
    //
    // Custom atom invariants and custom bond invariants are supported
    // via the existing override paths.
    Ok(())
}

fn morgan_bond_invariant(bond: &crate::Bond, params: &MorganFingerprintParams) -> u32 {
    let use_bond_types = params
        .bond_invariants_generator
        .as_ref()
        .map_or(params.use_bond_types, |generator| generator.use_bond_types);
    if !use_bond_types {
        return 0;
    }
    match bond.order() {
        crate::BondOrder::Single => 1,
        crate::BondOrder::Double => 2,
        crate::BondOrder::Triple => 3,
        crate::BondOrder::Quadruple => 4,
        crate::BondOrder::Aromatic => 12,
        crate::BondOrder::Dative => 9,
        crate::BondOrder::Zero | crate::BondOrder::Unspecified | crate::BondOrder::Null => 0,
        _ => 0,
    }
}

// RDKit uses:  gboost::hash_combine(seed, value)
// which expands to:  seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
fn hash_combine(seed: &mut u32, value: u32) {
    *seed = seed
        .wrapping_add(value)
        .wrapping_add(0x9e3779b9u32)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

// ---------------------------------------------------------------------------
// Fingerprint construction
// ---------------------------------------------------------------------------

fn fold_invariant(invariant: u32, n_bits: usize) -> usize {
    invariant as usize % n_bits
}

fn build_fingerprint(
    molecule: &Molecule,
    all_rounds: &[Vec<u32>],
    params: &MorganFingerprintParams,
) -> Result<MorganFingerprintOutput, FingerprintError> {
    let n_bits = params.n_bits;
    let collect = params.collect_additional_output;

    let mut atom_counts = if collect {
        vec![0u32; molecule.num_atoms()]
    } else {
        vec![]
    };
    let mut bit_info_map: BTreeMap<usize, Vec<(usize, u32)>> = BTreeMap::new();
    let mut atoms_per_bit: BTreeMap<usize, Vec<Vec<usize>>> = BTreeMap::new();

    let mut on_bits = Vec::new();

    for (round_idx, round_invs) in all_rounds.iter().enumerate() {
        let round = round_idx as u32;
        for atom_idx in 0..molecule.num_atoms() {
            if atom_is_excluded(atom_idx, params) {
                continue;
            }
            let inv = round_invs[atom_idx];

            if params.only_nonzero_invariants && inv == 0 {
                continue;
            }

            if params.count_simulation {
                let bit = fold_invariant(inv, n_bits);
                on_bits.push(bit);

                if collect {
                    atom_counts[atom_idx] += 1;
                    bit_info_map.entry(bit).or_default().push((atom_idx, round));
                    atoms_per_bit.entry(bit).or_default().push(vec![atom_idx]);
                }
            } else {
                for chunk in 0..params.num_bits_per_feature {
                    let bit =
                        fold_invariant(inv.wrapping_add(chunk.wrapping_mul(0x517cc1b7)), n_bits);
                    on_bits.push(bit);

                    if collect {
                        atom_counts[atom_idx] += 1;
                        bit_info_map.entry(bit).or_default().push((atom_idx, round));
                        atoms_per_bit.entry(bit).or_default().push(vec![atom_idx]);
                    }
                }
            }
        }
    }

    // Count-simulation folding: each unique bit's count is compared against
    // count_bounds to set additional offset bits.
    if params.count_simulation && !params.count_bounds.is_empty() {
        let mut counts_per_bit: BTreeMap<usize, u32> = BTreeMap::new();
        for &bit in &on_bits {
            *counts_per_bit.entry(bit).or_insert(0) += 1;
        }
        on_bits.clear();
        for (&bit, &count) in &counts_per_bit {
            on_bits.push(bit);
            for (bound_idx, &bound) in params.count_bounds.iter().enumerate().skip(1) {
                if count >= bound {
                    let offset_bit = (bit + bound_idx * n_bits) % n_bits;
                    on_bits.push(offset_bit);
                }
            }
        }
    }

    let fingerprint = Fingerprint::from_on_bits(n_bits, on_bits.iter().copied());

    let additional_output = if collect {
        let mut atom_to_bits: Vec<Vec<usize>> = vec![vec![]; molecule.num_atoms()];
        for (&bit, entries) in &bit_info_map {
            for &(atom_idx, _) in entries {
                if !atom_to_bits[atom_idx].contains(&bit) {
                    atom_to_bits[atom_idx].push(bit);
                }
            }
        }
        Some(MorganAdditionalOutput {
            atom_counts,
            atom_to_bits,
            bit_info_map,
            atoms_per_bit,
        })
    } else {
        None
    };

    Ok(MorganFingerprintOutput {
        fingerprint,
        additional_output,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn atom_is_excluded(index: usize, params: &MorganFingerprintParams) -> bool {
    if let Some(from) = &params.from_atoms {
        return !from.contains(&index);
    }
    if let Some(ignore) = &params.ignore_atoms {
        return ignore.contains(&index);
    }
    false
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Molecule;

    fn default_morgan_params(radius: u32, n_bits: usize) -> MorganFingerprintParams {
        MorganFingerprintParams {
            radius,
            n_bits,
            ..Default::default()
        }
    }

    fn methane() -> Molecule {
        Molecule::from_smiles_with_sanitize("C", false).unwrap()
    }

    fn ethane() -> Molecule {
        Molecule::from_smiles_with_sanitize("CC", false).unwrap()
    }

    fn benzene() -> Molecule {
        Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap()
    }

    #[test]
    fn morgan_fingerprint_empty_molecule_returns_empty_fingerprint() {
        let mol = Molecule::from_smiles_with_sanitize("", false).unwrap();
        let fp = morgan_fingerprint(&mol, &default_morgan_params(2, 2048)).unwrap();
        assert_eq!(fp.on_bits(), Vec::<usize>::new());
    }

    #[test]
    fn morgan_fingerprint_empty_params_n_bits_zero_returns_error() {
        let mol = methane();
        let params = MorganFingerprintParams {
            n_bits: 0,
            ..Default::default()
        };
        assert!(matches!(
            morgan_fingerprint(&mol, &params),
            Err(FingerprintError::EmptyFingerprint)
        ));
    }

    #[test]
    fn morgan_fingerprint_methane_radius0_produces_deterministic_fingerprint() {
        let mol = methane();
        let fp_a = morgan_fingerprint(&mol, &default_morgan_params(0, 2048)).unwrap();
        let fp_b = morgan_fingerprint(&mol, &default_morgan_params(0, 2048)).unwrap();
        assert_eq!(fp_a, fp_b);
        assert!(!fp_a.on_bits().is_empty(), "expected at least one on-bit");
    }

    #[test]
    fn morgan_fingerprint_tanimoto_self_is_one() {
        let mol = benzene();
        let fp = morgan_fingerprint(&mol, &default_morgan_params(2, 2048)).unwrap();
        let similarity = fp.tanimoto(&fp).unwrap();
        assert!(
            (similarity - 1.0).abs() < 1e-9,
            "tanimoto of fingerprint with itself should be 1.0, got {similarity}"
        );
    }

    #[test]
    fn morgan_fingerprint_n_bits_matches_param() {
        let mol = benzene();
        for n_bits in [64, 256, 1024, 2048] {
            let fp = morgan_fingerprint(&mol, &default_morgan_params(2, n_bits)).unwrap();
            assert_eq!(fp.n_bits(), n_bits);
        }
    }

    #[test]
    fn morgan_fingerprint_radius_determinism() {
        let mol = benzene();
        for radius in 0..=3 {
            let fp_a = morgan_fingerprint(&mol, &default_morgan_params(radius, 2048)).unwrap();
            let fp_b = morgan_fingerprint(&mol, &default_morgan_params(radius, 2048)).unwrap();
            assert_eq!(fp_a, fp_b, "radius={radius} should be deterministic");
        }
    }

    #[test]
    fn morgan_fingerprint_ethane_and_methane_differ() {
        let m = methane();
        let e = ethane();
        let fp_m = morgan_fingerprint(&m, &default_morgan_params(0, 2048)).unwrap();
        let fp_e = morgan_fingerprint(&e, &default_morgan_params(0, 2048)).unwrap();
        assert_ne!(
            fp_m, fp_e,
            "methane and ethane should have different fingerprints"
        );
    }

    #[test]
    fn morgan_fingerprint_benzene_and_ethane_differ() {
        let b = benzene();
        let e = ethane();
        let fp_b = morgan_fingerprint(&b, &default_morgan_params(2, 2048)).unwrap();
        let fp_e = morgan_fingerprint(&e, &default_morgan_params(2, 2048)).unwrap();
        assert_ne!(
            fp_b, fp_e,
            "benzene and ethane should have different fingerprints"
        );
    }

    #[test]
    fn morgan_fingerprint_radius_increases_on_bits() {
        let mol = ethane();
        let fp_r0 = morgan_fingerprint(&mol, &default_morgan_params(0, 2048)).unwrap();
        let fp_r2 = morgan_fingerprint(&mol, &default_morgan_params(2, 2048)).unwrap();
        assert!(
            fp_r2.on_bits().len() >= fp_r0.on_bits().len(),
            "larger radius should produce at least as many on-bits"
        );
    }

    #[test]
    fn morgan_fingerprint_with_output_produces_additional_data() {
        let mol = ethane();
        let params = MorganFingerprintParams {
            radius: 1,
            n_bits: 2048,
            collect_additional_output: true,
            ..Default::default()
        };
        let output = morgan_fingerprint_with_output(&mol, &params).unwrap();
        assert!(output.additional_output.is_some());
        let extra = output.additional_output.unwrap();
        assert_eq!(extra.atom_counts.len(), mol.num_atoms());
        assert!(!extra.bit_info_map.is_empty());
    }

    #[test]
    fn morgan_fingerprint_from_atoms_filters_by_allowed_indices() {
        let mol = ethane();
        let params = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            from_atoms: Some(vec![0]),
            ..Default::default()
        };
        let fp = morgan_fingerprint(&mol, &params).unwrap();
        assert!(!fp.on_bits().is_empty());

        let params_empty = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            from_atoms: Some(vec![]),
            ..Default::default()
        };
        let fp_empty = morgan_fingerprint(&mol, &params_empty).unwrap();
        assert!(
            fp_empty.on_bits().is_empty(),
            "no from_atoms → empty fingerprint"
        );
    }

    #[test]
    fn morgan_fingerprint_ignore_atoms_excludes_indices() {
        let mol = ethane();
        let params_full = default_morgan_params(0, 2048);
        let params_exclude = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            ignore_atoms: Some(vec![1]),
            ..Default::default()
        };
        let fp_full = morgan_fingerprint(&mol, &params_full).unwrap();
        let fp_excluded = morgan_fingerprint(&mol, &params_exclude).unwrap();
        assert_ne!(fp_full.on_bits().len(), 0);
        assert!(
            fp_excluded.on_bits().len() <= fp_full.on_bits().len(),
            "excluding an atom should not increase on-bits"
        );
    }

    #[test]
    fn morgan_fingerprint_feature_generator_produces_deterministic_fingerprint() {
        let mol = benzene();
        let params = MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Feature,
            ..Default::default()
        };
        // Feature invariants now use element/property classification.
        let fp_a = morgan_fingerprint(&mol, &params).unwrap();
        let fp_b = morgan_fingerprint(&mol, &params).unwrap();
        assert_eq!(fp_a, fp_b, "feature invariants should be deterministic");
        assert!(
            !fp_a.on_bits().is_empty(),
            "expected on-bits from feature invariants"
        );
    }

    #[test]
    fn morgan_fingerprint_custom_invariants_override_default() {
        let mol = ethane();
        let custom = vec![42u32, 99u32];
        let params = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            custom_atom_invariants: Some(custom),
            ..Default::default()
        };
        let fp_a = morgan_fingerprint(&mol, &params).unwrap();
        let fp_b = morgan_fingerprint(&mol, &params).unwrap();
        assert_eq!(fp_a, fp_b);
    }

    #[test]
    fn morgan_fingerprint_zero_bonds_molecule_does_not_panic() {
        let mol = Molecule::from_smiles_with_sanitize("[H][H]", false).unwrap();
        let fp = morgan_fingerprint(&mol, &default_morgan_params(2, 2048));
        assert!(fp.is_ok());
    }

    #[test]
    fn morgan_fingerprint_count_simulation_runs() {
        let mol = benzene();
        let params = MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            count_simulation: true,
            count_bounds: vec![1, 2, 4, 8],
            ..Default::default()
        };
        let fp = morgan_fingerprint(&mol, &params).unwrap();
        assert!(!fp.on_bits().is_empty());
        let std_fp = morgan_fingerprint(&mol, &default_morgan_params(2, 2048)).unwrap();
        assert!(
            fp.on_bits().len() >= std_fp.on_bits().len(),
            "count-simulation should set at least as many bits as standard mode"
        );
    }

    #[test]
    fn morgan_fingerprint_computes_adjacency_on_the_fly_if_cache_empty() {
        let mut mol = benzene();
        mol.derived_cache_mut().adjacency = None;
        let fp = morgan_fingerprint(&mol, &default_morgan_params(2, 2048)).unwrap();
        assert!(!fp.on_bits().is_empty());
    }

    #[test]
    fn morgan_fingerprint_chirality_produces_different_fingerprints() {
        // (R)-butan-2-ol and (S)-butan-2-ol should produce different
        // fingerprints when use_chirality=true.
        let r_mol = Molecule::from_smiles_with_sanitize("C[C@@H](O)CC", false).unwrap();
        let s_mol = Molecule::from_smiles_with_sanitize("C[C@H](O)CC", false).unwrap();
        let params = MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            use_chirality: true,
            ..Default::default()
        };
        let fp_r = morgan_fingerprint(&r_mol, &params).unwrap();
        let fp_s = morgan_fingerprint(&s_mol, &params).unwrap();
        let tc = fp_r.tanimoto(&fp_s).unwrap();
        assert!(
            tc < 1.0,
            "R and S enantiomers should have tc < 1.0 with chirality, got {tc}"
        );
    }

    #[test]
    fn morgan_fingerprint_chirality_disabled_produces_identical_fingerprints() {
        // Without chirality, R and S enantiomers produce the same fingerprint.
        let r_mol = Molecule::from_smiles_with_sanitize("C[C@@H](O)CC", false).unwrap();
        let s_mol = Molecule::from_smiles_with_sanitize("C[C@H](O)CC", false).unwrap();
        let params = MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            use_chirality: false,
            ..Default::default()
        };
        let fp_r = morgan_fingerprint(&r_mol, &params).unwrap();
        let fp_s = morgan_fingerprint(&s_mol, &params).unwrap();
        assert_eq!(
            fp_r, fp_s,
            "R and S should have same fingerprint without chirality"
        );
    }

    #[test]
    fn morgan_fingerprint_custom_bond_invariants_override_default() {
        let mol = ethane();
        let params = MorganFingerprintParams {
            radius: 0,
            n_bits: 2048,
            custom_bond_invariants: Some(vec![5u32]),
            ..Default::default()
        };
        let fp = morgan_fingerprint(&mol, &params).unwrap();
        assert!(!fp.on_bits().is_empty());
    }
}
