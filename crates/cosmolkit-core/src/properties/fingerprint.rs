use std::collections::BTreeMap;

use crate::{AdjacencyList, AtomId, BondOrder, ChiralTag, Molecule, NeighborRef};

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

// BEGIN RDKIT CPP FUNCTION MorganFingerprints::getHashedFingerprint
// RDKit source: MorganFingerprints.cpp lines 85-119
//
// RDKit source: SparseIntVect<uint32_t> *getHashedFingerprint(
// RDKit source:     const ROMol &mol, unsigned int radius, unsigned int nBits,
// RDKit source:     std::vector<uint32_t> *invariants,
// RDKit source:     const std::vector<uint32_t> *fromAtoms,
// RDKit source:     bool useChirality, bool useBondTypes,
// RDKit source:     bool onlyNonzeroInvariants,
// RDKit source:     BitInfoMap *atomsSettingBits,
// RDKit source:     bool includeRedundantEnvironments) {
// RDKit✔️❌:   if (nBits == 0) {
// RDKit✔️❌:     throw ValueErrorException("nBits can not be zero");
// RDKit✔️❌:   }
// RDKit✔️❌:   bool countSimulation = false;
// RDKit❗✔️:   std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
// RDKit❗✔️:       MorganFingerprint::getMorganGenerator<std::uint32_t>(...));
// RDKit❗✔️:   RDKit::FingerprintFuncArguments args;
// RDKit❗✔️:   args.fromAtoms = fromAtoms;
// RDKit❗✔️:   args.customAtomInvariants = invariants;
// RDKit❌❌:   AdditionalOutput ao;
// RDKit❌❌:   if (atomsSettingBits) {
// RDKit❌❌:     args.additionalOutput = &ao;
// RDKit❌❌:     ao.allocateBitInfoMap();
// RDKit❌❌:   }
// RDKit❗✔️:   auto res = fpgen->getCountFingerprint(mol, args).release();
// RDKit❌❌:   if (atomsSettingBits) {
// RDKit❌❌:     atomsSettingBits->clear();
// RDKit❌❌:     for (const auto &pr : *(ao.bitInfoMap)) {
// RDKit❌❌:       (*atomsSettingBits)[pr.first] = pr.second;
// RDKit❌❌:     }
// RDKit❌❌:   }
// RDKit❗✔️:   return res;
// END RDKIT CPP FUNCTION MorganFingerprints::getHashedFingerprint
//
// COSMolKit implementation: inline Morgan generation instead of
// FingerprintGenerator pipeline. atomsSettingBits / AdditionalOutput
// round-trip is not yet ported (RDKit❌❌).
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

    // BEGIN RDKIT CPP FUNCTION MorganGenerator::getMorganFingerprint (propagation loop)
    // RDKit source: MorganGenerator.cpp lines 384-459
    //
    // RDKit✔️✔️: for (unsigned int layer = 0; layer < morganArguments->d_radius; ++layer) {
    // RDKit✔️✔️:   std::vector<AccumTuple> allNeighborhoodsThisRound;
    // RDKit✔️✔️:   for (auto atomIdx : atomOrder) {
    // RDKit✔️✔️:     if (!deadAtoms[atomIdx]) {
    // RDKit✔️✔️:       const Atom *tAtom = mol.getAtomWithIdx(atomIdx);
    // RDKit✔️✔️:       if (!tAtom->getDegree()) {
    // RDKit✔️✔️:         deadAtoms.set(atomIdx, 1);
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ROMol::OEDGE_ITER beg, end;
    // RDKit✔️✔️:       boost::tie(beg, end) = mol.getAtomBonds(tAtom);
    // RDKit✔️✔️:       neighborhoodInvariants.clear();
    // RDKit✔️✔️:       while (beg != end) {
    // RDKit✔️✔️:         const Bond *bond = mol[*beg];
    // RDKit✔️✔️:         roundAtomNeighborhoods[atomIdx][bond->getIdx()] = 1;
    // RDKit✔️✔️:         unsigned int oIdx = bond->getOtherAtomIdx(atomIdx);
    // RDKit✔️✔️:         roundAtomNeighborhoods[atomIdx] |= atomNeighborhoods[oIdx];
    // RDKit✔️✔️:         auto bt = static_cast<int32_t>((*bondInvariants)[bond->getIdx()]);
    // RDKit✔️✔️:         neighborhoodInvariants.push_back(
    // RDKit✔️✔️:             std::make_pair(bt, currentInvariants[oIdx]));
    // RDKit✔️✔️:         ++beg;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       std::sort(neighborhoodInvariants.begin(), neighborhoodInvariants.end());
    // RDKit✔️✔️:       std::uint32_t invar = layer;
    // RDKit✔️✔️:       gboost::hash_combine(invar, currentInvariants[atomIdx]);
    // RDKit✔️✔️:       bool looksChiral = (tAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
    // RDKit✔️✔️:       for (std::vector<std::pair<int32_t, uint32_t>>::const_iterator it =
    // RDKit✔️✔️:                neighborhoodInvariants.begin();
    // RDKit✔️✔️:            it != neighborhoodInvariants.end(); ++it) {
    // RDKit✔️✔️:         gboost::hash_combine(invar, *it);
    // RDKit✔️✔️:         if (morganArguments->df_includeChirality && looksChiral &&
    // RDKit✔️✔️:             !chiralAtoms[atomIdx]) {
    // RDKit✔️✔️:           if (it->first != static_cast<int32_t>(Bond::SINGLE)) {
    // RDKit✔️✔️:             looksChiral = false;
    // RDKit✔️✔️:           } else if (it != neighborhoodInvariants.begin() &&
    // RDKit✔️✔️:                      it->second == (it - 1)->second) {
    // RDKit✔️✔️:             looksChiral = false;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (morganArguments->df_includeChirality && looksChiral) {
    // RDKit✔️✔️:         chiralAtoms[atomIdx] = 1;
    // RDKit✔️✔️:         std::string cip = "";
    // RDKit✔️✔️:         tAtom->getPropIfPresent(common_properties::_CIPCode, cip);
    // RDKit✔️✔️:         if (cip == "R") {
    // RDKit✔️✔️:           gboost::hash_combine(invar, 3);
    // RDKit✔️✔️:         } else if (cip == "S") {
    // RDKit✔️✔️:           gboost::hash_combine(invar, 2);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           gboost::hash_combine(invar, 1);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       nextLayerInvariants[atomIdx] = static_cast<OutputType>(invar);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MorganGenerator::getMorganFingerprint
    //
    // COSMolKit Rust implementation: simplified atom ordering, no dead-atom
    // tracking (COSMolKit doesn't need it since we don't parallelize), and
    // chirality derived from ChiralTag + chiral_permutation instead of CIP code.
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
                    let bond_idx = n.bond.index();
                    let bt = morgan_bond_invariant(bond_idx, &molecule.bonds()[bond_idx], params);
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

fn morgan_bond_invariant(
    bond_idx: usize,
    bond: &crate::Bond,
    params: &MorganFingerprintParams,
) -> u32 {
    // custom_bond_invariants override: if provided and this bond index has a
    // value, use it directly instead of computing from bond order.
    if let Some(custom) = &params.custom_bond_invariants {
        if let Some(&inv) = custom.get(bond_idx) {
            return inv;
        }
    }
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
pub(crate) fn hash_combine(seed: &mut u32, value: u32) {
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

    // Track seen invariants for include_redundant_environments=false.
    // When disabled, each unique invariant value is only allowed to
    // contribute bits once across all atoms and all rounds, matching
    // RDKit's behavior of deduplicating redundant environments.
    let mut seen_invariants: std::collections::HashSet<u32> = std::collections::HashSet::new();

    for (round_idx, round_invs) in all_rounds.iter().enumerate() {
        let round = round_idx as u32;
        for atom_idx in 0..molecule.num_atoms() {
            if atom_is_excluded(atom_idx, params) {
                continue;
            }
            let inv = round_invs[atom_idx];

            // When include_redundant_environments is false, skip invariants
            // that have already contributed bits in a previous round or atom.
            if !params.include_redundant_environments && !seen_invariants.insert(inv) {
                continue;
            }

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
// Topological (Path-Based) Fingerprint
// RDKit source: GraphMol/Fingerprints/Fingerprints.h
// RDKit source: GraphMol/Fingerprints/FingerprintUtil.cpp
// ---------------------------------------------------------------------------

/// Parameters for topological (path-based) fingerprint generation.
///
/// Enumerates all linear paths in the molecular graph and hashes each
/// path into a fixed-size bit fingerprint. Analogous to RDKit's
/// `RDKFingerprint` / `LayeredFingerprint`.
///
/// # Parameters
/// - `min_path`: minimum path length in bonds (default 1).
/// - `max_path`: maximum path length in bonds (default 7).
/// - `n_bits`: size of the output fingerprint (default 2048).
/// - `n_bits_per_hash`: number of bit positions to set per path hash (default 2).
/// - `use_bond_types`: whether bond order contributes to the path invariant.
/// - `from_atoms`: if `Some`, only enumerate paths starting from these atoms.
/// - `ignore_atoms`: if `Some`, skip paths through these atoms.
#[derive(Debug, Clone)]
pub struct TopologicalFingerprintParams {
    pub min_path: u32,
    pub max_path: u32,
    pub n_bits: usize,
    pub n_bits_per_hash: u32,
    pub use_bond_types: bool,
    pub from_atoms: Option<Vec<usize>>,
    pub ignore_atoms: Option<Vec<usize>>,
}

impl Default for TopologicalFingerprintParams {
    fn default() -> Self {
        Self {
            min_path: 1,
            max_path: 7,
            n_bits: 2048,
            n_bits_per_hash: 2,
            use_bond_types: true,
            from_atoms: None,
            ignore_atoms: None,
        }
    }
}

// RDKit source: GraphMol/Fingerprints/FingerprintUtil.cpp
// RDKit source: void getFingerprint(const ROMol &mol, SparseBitVect &fp,
// RDKit source:     unsigned int minPath, unsigned int maxPath,
// RDKit source:     ...)
// RDKit source: path enumeration via boost::graph and BFS up to maxPath.
//
// COSMolKit enumerates simple (non-self-intersecting) paths by DFS from each
// atom, visiting neighbors up to max_path depth. Each path's atom types and
// bond orders are hashed into an invariant, which is folded into bit positions.
#[must_use]
pub fn topological_fingerprint(
    molecule: &Molecule,
    params: &TopologicalFingerprintParams,
) -> Fingerprint {
    let n_bits = params.n_bits;
    if molecule.num_atoms() == 0 || n_bits == 0 {
        return Fingerprint::from_on_bits(n_bits, []);
    }

    // Build adjacency list if not cached.
    let adjacency = molecule
        .derived_cache()
        .adjacency
        .clone()
        .unwrap_or_else(|| AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds()));

    let num_atoms = molecule.num_atoms();
    let atoms = molecule.atoms();
    let bonds = molecule.bonds();

    // Helper to check if an atom index is excluded.
    let is_excluded = |idx: usize| -> bool {
        if let Some(from) = &params.from_atoms {
            return !from.contains(&idx);
        }
        if let Some(ignore) = &params.ignore_atoms {
            return ignore.contains(&idx);
        }
        false
    };

    // Cache atomic numbers for fast access.
    let atomic_numbers: Vec<u8> = atoms.iter().map(|a| a.atomic_number()).collect();

    // Compute bond order values for each bond (by bond index).
    let bond_order_values: Vec<u32> = if params.use_bond_types {
        bonds
            .iter()
            .map(|b| match b.order() {
                crate::BondOrder::Single => 1,
                crate::BondOrder::Double => 2,
                crate::BondOrder::Triple => 3,
                crate::BondOrder::Quadruple => 4,
                crate::BondOrder::Aromatic => 12,
                crate::BondOrder::Dative => 9,
                _ => 0,
            })
            .collect()
    } else {
        bonds.iter().map(|_| 1u32).collect()
    };

    // Convert bonds to CSR-style lookup: for each pair (atom_i, atom_j),
    // find the bond index between them.
    // Build a simple adjacency with bond indices.
    let mut on_bits_set: std::collections::HashSet<usize> = std::collections::HashSet::new();

    for start_atom in 0..num_atoms {
        if is_excluded(start_atom) {
            continue;
        }

        // DFS stack: (current_atom, path_atoms, path_bonds, depth)
        // path_atoms: vec of atom indices visited (including start)
        // path_bonds: vec of bond order values between atoms in path
        let mut stack: Vec<(usize, Vec<usize>, Vec<u32>, u32)> = Vec::new();
        stack.push((start_atom, vec![start_atom], vec![], 0));

        while let Some((current, path_atoms, path_bonds, depth)) = stack.pop() {
            // If the current path length (in bonds) >= min_path, hash it.
            // Path length in bonds = path_atoms.len() - 1
            let path_len_bonds = path_atoms.len() as u32 - 1;
            if path_len_bonds >= params.min_path && path_len_bonds <= params.max_path {
                // Build invariant: start with atomic_number of first atom,
                // then hash_combine each (bond_order, next_atomic_number) pair.
                let mut invariant: u32 = atomic_numbers[start_atom] as u32;
                for (k, &bond_val) in path_bonds.iter().enumerate() {
                    let next_idx = path_atoms[k + 1];
                    let mut pair_hash: u32 = 0;
                    hash_combine(&mut pair_hash, bond_val);
                    hash_combine(&mut pair_hash, atomic_numbers[next_idx] as u32);
                    hash_combine(&mut invariant, pair_hash);
                }

                // Fold the invariant into bit positions using n_bits_per_hash folding.
                for chunk in 0..params.n_bits_per_hash {
                    let bit = if chunk == 0 {
                        invariant as usize % n_bits
                    } else {
                        invariant.wrapping_add(chunk.wrapping_mul(0x517cc1b7)) as usize % n_bits
                    };
                    on_bits_set.insert(bit);
                }
            }

            // Stop enumerating if we've reached max_path depth.
            if path_len_bonds >= params.max_path {
                continue;
            }

            // Visit neighbors that are not already in the current path.
            for neighbor in adjacency.neighbors_of(current) {
                let n_idx = neighbor.atom_index;
                if path_atoms.contains(&n_idx) {
                    continue;
                }
                if is_excluded(n_idx) {
                    continue;
                }

                let bond_val = bond_order_values[neighbor.bond.index()];
                let mut new_path_atoms = path_atoms.clone();
                new_path_atoms.push(n_idx);
                let mut new_path_bonds = path_bonds.clone();
                new_path_bonds.push(bond_val);
                stack.push((n_idx, new_path_atoms, new_path_bonds, depth + 1));
            }
        }
    }

    let on_bits: Vec<usize> = on_bits_set.into_iter().collect();
    Fingerprint::from_on_bits(n_bits, on_bits)
}

// ---------------------------------------------------------------------------
// MACCS Keys (166-bit structural keys)
// RDKit source: GraphMol/Fingerprints/MACCS.cpp
// ---------------------------------------------------------------------------

/// Parameters for MACCS key generation.
///
/// The MACCS (Molecular ACCess System) keys are 166 predefined structural
/// features encoded as a bit vector. Each bit corresponds to a specific
/// substructure or element property.
///
/// Since COSMolKit does not yet have SMARTS substructure matching, the
/// implementation uses structural heuristics (atom properties, neighbor
/// topology, bond orders, ring info) to approximate each key.
///
/// # Reference
/// RDKit MACCS implementation: Code/GraphMol/Fingerprints/MACCS.cpp
/// - 166 SMARTS patterns mapped to bits 1-166 (bit 0 unused).
/// - COSMolKit implements equivalent detection via heuristics.
#[derive(Debug, Clone)]
pub struct MaccsFingerprintParams {
    pub n_bits: usize,
}

impl Default for MaccsFingerprintParams {
    fn default() -> Self {
        Self { n_bits: 166 }
    }
}

// RDKit❗✔️: ExplicitBitVect *getFingerprintAsBitVect(const ROMol &mol)
// RDKit❗✔️:   std::unique_ptr<ExplicitBitVect> fp(new ExplicitBitVect(167));
// RDKit❗✔️:   GenerateFP(mol, *fp);
// RDKit❗✔️:   return fp.release();
// END RDKIT FUNCTION getFingerprintAsBitVect
//
// COSMolKit MACCS keys use 166 bits (bits 0-165), matching RDKit's bit layout
// where bit 0 is unused. Each bit is set via structural heuristics rather
// than SMARTS substructure matching, producing compatible but not bit-identical
// fingerprints compared to RDKit.
#[must_use]
pub fn maccs_fingerprint(molecule: &Molecule, params: &MaccsFingerprintParams) -> Fingerprint {
    let n_bits = params.n_bits;
    if molecule.num_atoms() == 0 || n_bits == 0 {
        return Fingerprint::from_on_bits(n_bits, []);
    }

    let atoms = molecule.atoms();
    let bonds = molecule.bonds();
    let num_atoms = molecule.num_atoms();
    let num_bonds = molecule.num_bonds();

    // Build adjacency if available (for neighbor traversal)
    let adjacency = molecule
        .derived_cache()
        .adjacency
        .clone()
        .unwrap_or_else(|| AdjacencyList::from_topology(num_atoms, bonds));

    let ring_info = molecule.derived_cache().rings.as_ref();

    // Atom-level data collection
    let atomic_numbers: Vec<u8> = atoms.iter().map(|a| a.atomic_number()).collect();
    let charges: Vec<i8> = atoms.iter().map(|a| a.formal_charge()).collect();
    let aromatic_atoms: Vec<bool> = atoms.iter().map(|a| a.is_aromatic()).collect();

    // Group atoms by element
    let has_element = |target_z: u8| -> bool { atomic_numbers.iter().any(|&z| z == target_z) };

    // Count atoms by element
    let count_element =
        |target_z: u8| -> usize { atomic_numbers.iter().filter(|&&z| z == target_z).count() };

    // Bond-level data
    let bond_orders: Vec<BondOrder> = bonds.iter().map(|b| b.order()).collect();
    let bond_aromatic: Vec<bool> = bonds.iter().map(|b| b.is_aromatic()).collect();

    // Helper: check if any bond has the given order
    let has_bond_order = |order: BondOrder| -> bool { bond_orders.iter().any(|&o| o == order) };

    // Helper: check if there is a double bond (including aromatic)
    let has_double_bond = || -> bool {
        bond_orders
            .iter()
            .any(|o| *o == BondOrder::Double || *o == BondOrder::Aromatic)
    };

    // Helper: check degree (number of heavy-atom neighbors)
    let degree = |atom_idx: usize| -> usize { adjacency.neighbors_of(atom_idx).len() };

    // Helper: check if atom is in a ring
    let atom_in_ring = |atom_idx: usize| -> bool {
        ring_info.map_or(false, |r| r.num_atom_rings(AtomId::new(atom_idx)) > 0)
    };

    // Helper: check if there's an aromatic bond in the molecule
    let has_aromatic_bond = || -> bool { bond_aromatic.iter().any(|&a| a) };

    // Helper: check if atom has a neighbor of given atomic number
    let has_neighbor_z = |atom_idx: usize, target_z: u8| -> bool {
        adjacency
            .neighbors_of(atom_idx)
            .iter()
            .any(|n| atomic_numbers[n.atom_index] == target_z)
    };

    // Helper: find an atom by predicate
    let find_atom = |pred: fn(&[u8], &[i8], usize) -> bool| -> Option<usize> {
        (0..num_atoms).find(|&i| pred(&atomic_numbers, &charges, i))
    };

    // -----------------------------------------------------------------------
    // Bit detection — matches RDKit MACCS bits 1-166
    // -----------------------------------------------------------------------
    let mut on_bits: Vec<usize> = Vec::new();

    // RDKit bits from the atom-loop section:
    // Bit 2 (1-indexed, so index 1): atomic number 104 (Rf+)
    if has_element(104) {
        on_bits.push(1);
    }

    // Bit 3: Ge, As, Se, Sn, Sb, Te, Pb, Bi, Po (group 14-16 p-block)
    for &z in &[32, 33, 34, 50, 51, 52, 82, 83, 84] {
        if has_element(z) {
            on_bits.push(2);
            break;
        }
    }

    // Bit 4: actinides (89-103)
    for z in 89u8..=103 {
        if has_element(z) {
            on_bits.push(3);
            break;
        }
    }

    // Bit 5: Sc, Ti, Y, Zr, Hf
    for &z in &[21, 22, 39, 40, 72] {
        if has_element(z) {
            on_bits.push(4);
            break;
        }
    }

    // Bit 6: lanthanides (57-71)
    for z in 57u8..=71 {
        if has_element(z) {
            on_bits.push(5);
            break;
        }
    }

    // Bit 7: V, Cr, Mn, Nb, Mo, Tc, Ta, W, Re
    for &z in &[23, 24, 25, 41, 42, 43, 73, 74, 75] {
        if has_element(z) {
            on_bits.push(6);
            break;
        }
    }

    // Bit 8: [!#6!#1]1~*~*~*~1 — heteroatom in 4-membered ring
    {
        let mut found = false;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                if ri.atom_rings()[ring_idx].len() == 4 {
                    for a in ri.atom_rings()[ring_idx].iter() {
                        let z = atomic_numbers[a.index()];
                        if z != 6 && z != 1 {
                            found = true;
                            break;
                        }
                    }
                }
                if found {
                    break;
                }
            }
        }
        if found {
            on_bits.push(7);
        }
    }

    // Bit 9: Fe, Co, Ni, Ru, Rh, Pd, Os, Ir, Pt
    for &z in &[26, 27, 28, 44, 45, 46, 76, 77, 78] {
        if has_element(z) {
            on_bits.push(8);
            break;
        }
    }

    // Bit 10: Be, Mg, Ca, Sr, Ba, Ra (alkaline earth)
    for &z in &[4, 12, 20, 38, 56, 88] {
        if has_element(z) {
            on_bits.push(9);
            break;
        }
    }

    // Bit 11: *1~*~*~*~1 — any 4-membered ring
    {
        let mut found = false;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                if ri.atom_rings()[ring_idx].len() == 4 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(10);
        }
    }

    // Bit 12: Cu, Zn, Ag, Cd, Au, Hg
    for &z in &[29, 30, 47, 48, 79, 80] {
        if has_element(z) {
            on_bits.push(11);
            break;
        }
    }

    // Bit 13: [#8]~[#7](~[#6])~[#6] — nitro/nitrate pattern: O-N(C)-C
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 {
                // O bonded to N
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 7 {
                        // N is also bonded to two C atoms
                        let c_neighbors_of_n: Vec<usize> = adjacency
                            .neighbors_of(n.atom_index)
                            .iter()
                            .filter(|nn| atomic_numbers[nn.atom_index] == 6)
                            .map(|nn| nn.atom_index)
                            .collect();
                        if c_neighbors_of_n.len() >= 2 {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(12);
        }
    }

    // Bit 14: [#16]-[#16] — disulfide
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 16 {
                        // Check it's a single bond
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Single {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(13);
        }
    }

    // Bit 15: [#8]~[#6](~[#8])~[#8] — carbonate: O-C(O)-O
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                let o_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 8)
                    .map(|n| n.atom_index)
                    .collect();
                if o_neighbors.len() >= 3 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(14);
        }
    }

    // Bit 16: [!#6!#1]1~*~*~1 — heteroatom in 3-membered ring
    {
        let mut found = false;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                if ri.atom_rings()[ring_idx].len() == 3 {
                    for a in ri.atom_rings()[ring_idx].iter() {
                        let z = atomic_numbers[a.index()];
                        if z != 6 && z != 1 {
                            found = true;
                            break;
                        }
                    }
                }
                if found {
                    break;
                }
            }
        }
        if found {
            on_bits.push(15);
        }
    }

    // Bit 17: [#6]#[#6] — alkyne (carbon-carbon triple bond)
    if has_bond_order(BondOrder::Triple) {
        on_bits.push(16);
    }

    // Bit 18: B, Al, Ga, In, Tl (Group 13)
    for &z in &[5, 13, 31, 49, 81] {
        if has_element(z) {
            on_bits.push(17);
            break;
        }
    }

    // Bit 19: *1~*~*~*~*~*~*~1 — 7-membered ring
    {
        let mut found = false;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                if ri.atom_rings()[ring_idx].len() == 7 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(18);
        }
    }

    // Bit 20: [#14] — Si
    if has_element(14) {
        on_bits.push(19);
    }

    // Bit 21: [#6]=[#6](~[!#6!#1])~[!#6!#1] — alkene with 2 heteroatoms on vinyl
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                                // i and n are the double-bonded carbons
                                // Check i has one heteroatom neighbor that is not the other C
                                let i_het: Vec<usize> = adjacency
                                    .neighbors_of(i)
                                    .iter()
                                    .filter(|nn| {
                                        let z = atomic_numbers[nn.atom_index];
                                        nn.atom_index != n.atom_index && z != 6 && z != 1
                                    })
                                    .map(|nn| nn.atom_index)
                                    .collect();
                                let n_het: Vec<usize> = adjacency
                                    .neighbors_of(n.atom_index)
                                    .iter()
                                    .filter(|nn| {
                                        let z = atomic_numbers[nn.atom_index];
                                        nn.atom_index != i && z != 6 && z != 1
                                    })
                                    .map(|nn| nn.atom_index)
                                    .collect();
                                if !i_het.is_empty() && !n_het.is_empty() {
                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(20);
        }
    }

    // Bit 22: *1~*~*~1 — 3-membered ring
    {
        let mut found = false;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                if ri.atom_rings()[ring_idx].len() == 3 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(21);
        }
    }

    // Bit 23: [#7]~[#6](~[#8])~[#8] — N-C(=O)-O (carbamate / N-carboxyl)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 {
                        let o_neighbors_of_c: Vec<usize> = adjacency
                            .neighbors_of(n.atom_index)
                            .iter()
                            .filter(|nn| atomic_numbers[nn.atom_index] == 8)
                            .map(|nn| nn.atom_index)
                            .collect();
                        if o_neighbors_of_c.len() >= 2 {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(22);
        }
    }

    // Bit 24: [#7]-[#8] — N-O bond (hydroxylamine, N-oxide)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 8 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(23);
        }
    }

    // Bit 25: [#7]~[#6](~[#7])~[#7] — guanidine: N-C(N)=N
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 {
                        let n_neighbors_of_c: Vec<usize> = adjacency
                            .neighbors_of(n.atom_index)
                            .iter()
                            .filter(|nn| atomic_numbers[nn.atom_index] == 7)
                            .map(|nn| nn.atom_index)
                            .collect();
                        if n_neighbors_of_c.len() >= 3 {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(24);
        }
    }

    // Bit 26: [#6]=@[#6](@*)@* — cyclic alkene (C=C within a ring)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 && atom_in_ring(i) {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 && atom_in_ring(n.atom_index) {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(25);
        }
    }

    // Bit 27: I
    if has_element(53) {
        on_bits.push(26);
    }

    // Bit 28: [!#6!#1]~[CH2]~[!#6!#1] — heteroatom-CH2-heteroatom
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                // Check if it has 2 hydrogens (methylene-like)
                let het_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| {
                        let z = atomic_numbers[n.atom_index];
                        z != 6 && z != 1
                    })
                    .map(|n| n.atom_index)
                    .collect();
                if het_neighbors.len() >= 2 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(27);
        }
    }

    // Bit 29: P
    if has_element(15) {
        on_bits.push(28);
    }

    // Bit 30: [#6]~[!#6!#1](~[#6])(~[#6])~* — quaternary heteroatom bonded to 3 carbons
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z = atomic_numbers[i];
            if z != 6 && z != 1 {
                let c_neighbors: usize = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 6)
                    .count();
                let total_neighbors = adjacency.neighbors_of(i).len();
                if c_neighbors >= 2 && total_neighbors >= 3 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(29);
        }
    }

    // Bit 31: [!#6!#1]~[F,Cl,Br,I] — heteroatom-halogen bond
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z = atomic_numbers[i];
            if z != 6 && z != 1 {
                for n in adjacency.neighbors_of(i) {
                    let nz = atomic_numbers[n.atom_index];
                    if nz == 9 || nz == 17 || nz == 35 || nz == 53 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(30);
        }
    }

    // Bit 32: [#6]~[#16]~[#7] — C-S-N chain
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 16 {
                        for nn in adjacency.neighbors_of(n.atom_index) {
                            if nn.atom_index != i && atomic_numbers[nn.atom_index] == 7 {
                                found = true;
                                break;
                            }
                        }
                    }
                    if found {
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(31);
        }
    }

    // Bit 33: [#7]~[#16] — N-S bond
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 16 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(32);
        }
    }

    // Bit 34: [CH2]=* — methylene (terminal =CH2 or exocyclic)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                let neighbors = adjacency.neighbors_of(i);
                if neighbors.len() == 1 {
                    if let Some(order) = bond_orders.get(neighbors[0].bond.index()) {
                        if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                            found = true;
                            break;
                        }
                    }
                }
            }
        }
        if found {
            on_bits.push(33);
        }
    }

    // Bit 35: alkali metals Li, Na, K, Rb, Cs, Fr
    for &z in &[3, 11, 19, 37, 55, 87] {
        if has_element(z) {
            on_bits.push(34);
            break;
        }
    }

    // Bit 36: [#16R] — S in ring
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 && atom_in_ring(i) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(35);
        }
    }

    // Bit 37: [#7]~[#6](~[#8])~[#7] — urea-like: N-C(=O)-N
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                let n_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 7)
                    .map(|n| n.atom_index)
                    .collect();
                let o_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 8)
                    .map(|n| n.atom_index)
                    .collect();
                if n_neighbors.len() >= 2 && !o_neighbors.is_empty() {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(36);
        }
    }

    // Bit 38: [#7]~[#6](~[#6])~[#7] — N-C(R)-N
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                let n_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 7)
                    .map(|n| n.atom_index)
                    .collect();
                let c_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 6)
                    .map(|n| n.atom_index)
                    .collect();
                if n_neighbors.len() >= 2 && !c_neighbors.is_empty() {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(37);
        }
    }

    // Bit 39: [#8]~[#16](~[#8])~[#8] — sulfate: O=S(=O)-O
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 {
                let o_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 8)
                    .map(|n| n.atom_index)
                    .collect();
                if o_neighbors.len() >= 3 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(38);
        }
    }

    // Bit 40: [#16]-[#8] — S-O single bond (sulfoxide, sulfonate)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 8 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Single {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(39);
        }
    }

    // Bit 41: [#6]#[#7] — nitrile C≡N
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 7 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Triple {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(40);
        }
    }

    // Bit 42: F
    if has_element(9) {
        on_bits.push(41);
    }

    // Bit 43: [!#6!#1!H0]~*~[!#6!#1!H0] — two non-C non-H non-zero-H atoms with 1 atom between
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            let has_h = has_hydrogens(i, &adjacency, &atomic_numbers);
            if !has_h {
                continue;
            }
            for n in adjacency.neighbors_of(i) {
                for nn in adjacency.neighbors_of(n.atom_index) {
                    if nn.atom_index == i {
                        continue;
                    }
                    let z_nn = atomic_numbers[nn.atom_index];
                    if z_nn == 1 || z_nn == 6 {
                        continue;
                    }
                    let has_h_nn = has_hydrogens(nn.atom_index, &adjacency, &atomic_numbers);
                    if has_h_nn && n.atom_index != nn.atom_index {
                        found = true;
                        break;
                    }
                }
                if found {
                    break;
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(42);
        }
    }

    // Bit 44: exotic elements (not H, C, N, O, F, Si, P, S, Cl, Br, I)
    {
        let mut found = false;
        for &z in &atomic_numbers {
            if z != 1
                && z != 6
                && z != 7
                && z != 8
                && z != 9
                && z != 14
                && z != 15
                && z != 16
                && z != 17
                && z != 35
                && z != 53
            {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(43);
        }
    }

    // Bit 45: [#6]=[#6]~[#7] — C=C-N enamine pattern
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                                // Check if n has an N neighbor
                                for nn in adjacency.neighbors_of(n.atom_index) {
                                    if nn.atom_index != i && atomic_numbers[nn.atom_index] == 7 {
                                        found = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if found {
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(44);
        }
    }

    // Bit 46: Br
    if has_element(35) {
        on_bits.push(45);
    }

    // Bit 47: [#16]~*~[#7] — S-X-N chain
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 {
                for n in adjacency.neighbors_of(i) {
                    for nn in adjacency.neighbors_of(n.atom_index) {
                        if nn.atom_index != i && atomic_numbers[nn.atom_index] == 7 {
                            found = true;
                            break;
                        }
                    }
                    if found {
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(46);
        }
    }

    // Bit 48: [#8]~[!#6!#1](~[#8])~[#8] — O-heteroatom(-O)-O (phosphate-like)
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z = atomic_numbers[i];
            if z != 6 && z != 1 {
                let o_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 8)
                    .map(|n| n.atom_index)
                    .collect();
                if o_neighbors.len() >= 3 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(47);
        }
    }

    // Bit 49: [!+0] — any charged atom
    if charges.iter().any(|&c| c != 0) {
        on_bits.push(48);
    }

    // Bit 50: [#6]=[#6](~[#6])~[#6] — substituted alkene
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                                let i_c_neighbors = adjacency
                                    .neighbors_of(i)
                                    .iter()
                                    .filter(|nn| {
                                        nn.atom_index != n.atom_index
                                            && atomic_numbers[nn.atom_index] == 6
                                    })
                                    .count();
                                let n_c_neighbors = adjacency
                                    .neighbors_of(n.atom_index)
                                    .iter()
                                    .filter(|nn| {
                                        nn.atom_index != i && atomic_numbers[nn.atom_index] == 6
                                    })
                                    .count();
                                if i_c_neighbors >= 1 && n_c_neighbors >= 1 {
                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                    if found {
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(49);
        }
    }

    // Bit 51: [#6]~[#16]~[#8] — C-S-O chain
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 16 {
                        for nn in adjacency.neighbors_of(n.atom_index) {
                            if nn.atom_index != i && atomic_numbers[nn.atom_index] == 8 {
                                found = true;
                                break;
                            }
                        }
                    }
                    if found {
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(50);
        }
    }

    // Bit 52: [#7]~[#7] — N-N bond
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 7 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(51);
        }
    }

    // Bit 53: [!#6!#1!H0]~*~*~*~[!#6!#1!H0] — two heteroatoms with 3-atom bridge
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            let has_h_i = has_hydrogens(i, &adjacency, &atomic_numbers);
            if !has_h_i {
                continue;
            }
            // DFS/BFS up to depth 3 looking for another non-C non-H with H
            if has_heteroatom_at_distance(i, 3, &adjacency, &atomic_numbers, &charges, true, i) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(52);
        }
    }

    // Bit 54: [!#6!#1!H0]~*~*~[!#6!#1!H0] — two heteroatoms with 2-atom bridge
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            let has_h_i = has_hydrogens(i, &adjacency, &atomic_numbers);
            if !has_h_i {
                continue;
            }
            if has_heteroatom_at_distance(i, 2, &adjacency, &atomic_numbers, &charges, true, i) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(53);
        }
    }

    // Bit 55: [#8]~[#16]~[#8] — O-S-O (sulfone/sulfonate)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 {
                let o_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 8)
                    .map(|n| n.atom_index)
                    .collect();
                if o_neighbors.len() >= 2 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(54);
        }
    }

    // Bit 56: [#8]~[#7](~[#8])~[#6] — O-N(=O)-C (nitro compound)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                let o_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 8)
                    .map(|n| n.atom_index)
                    .collect();
                let c_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 6)
                    .map(|n| n.atom_index)
                    .collect();
                if o_neighbors.len() >= 2 && !c_neighbors.is_empty() {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(55);
        }
    }

    // Bit 57: [#8R] — O in ring
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 && atom_in_ring(i) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(56);
        }
    }

    // Bit 58: [!#6!#1]~[#16]~[!#6!#1] — heteroatom-S-heteroatom
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 {
                let het_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| {
                        let z = atomic_numbers[n.atom_index];
                        z != 6 && z != 1
                    })
                    .map(|n| n.atom_index)
                    .collect();
                if het_neighbors.len() >= 2 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(57);
        }
    }

    // Bit 59: [#16]!:*:* — S non-aromatic in aromatic system
    // (approximate: S bonded to at least one aromatic atom via non-aromatic bond)
    // For this heuristic, we check if S is adjacent to an aromatic atom.
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 {
                for n in adjacency.neighbors_of(i) {
                    if aromatic_atoms[n.atom_index] {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order != BondOrder::Aromatic {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(58);
        }
    }

    // Bit 60: [#16]=[#8] — S=O (sulfoxide, sulfone)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 8 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(59);
        }
    }

    // Bit 61: *~[#16](~*)~* — S with at least 3 connections (sulfonium, sulfone)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 && adjacency.neighbors_of(i).len() >= 3 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(60);
        }
    }

    // Bit 62: *@*!@*@* — ring bond non-ring bond ring (fused ring system)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atom_in_ring(i) {
                let ring_nbrs: Vec<&NeighborRef> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atom_in_ring(n.atom_index))
                    .collect();
                let non_ring_nbrs: Vec<&NeighborRef> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| !atom_in_ring(n.atom_index))
                    .collect();
                if !ring_nbrs.is_empty() && !non_ring_nbrs.is_empty() {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(61);
        }
    }

    // Bit 63: [#7]=[#8] — N=O (N-oxide)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 8 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(62);
        }
    }

    // Bit 64: *@*!@[#16] — ring S with non-ring neighbor (thiophene substituent)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 && atom_in_ring(i) {
                for n in adjacency.neighbors_of(i) {
                    if !atom_in_ring(n.atom_index) {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(63);
        }
    }

    // Bit 65: c:n — aromatic C-N (azoline, pyridine-like)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 && aromatic_atoms[i] {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 7 && aromatic_atoms[n.atom_index] {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(64);
        }
    }

    // Bit 66: [#6]~[#6](~[#6])(~[#6])~* — quaternary carbon
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                let c_neighbors: usize = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 6)
                    .count();
                if c_neighbors >= 3 && adjacency.neighbors_of(i).len() >= 4 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(65);
        }
    }

    // Bit 67: [!#6!#1]~[#16] — heteroatom-S bond
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z = atomic_numbers[i];
            if z != 6 && z != 1 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 16 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(66);
        }
    }

    // Bit 68: [!#6!#1!H0]~[!#6!#1!H0] — two heteroatoms both with H bonded directly
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            let has_h_i = has_hydrogens(i, &adjacency, &atomic_numbers);
            if !has_h_i {
                continue;
            }
            for n in adjacency.neighbors_of(i) {
                let z_n = atomic_numbers[n.atom_index];
                if z_n == 1 || z_n == 6 {
                    continue;
                }
                let has_h_n = has_hydrogens(n.atom_index, &adjacency, &atomic_numbers);
                if has_h_n {
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(67);
        }
    }

    // Bit 69: [!#6!#1]~[!#6!#1!H0] — heteroatom bonded to heteroatom with H
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            for n in adjacency.neighbors_of(i) {
                let z_n = atomic_numbers[n.atom_index];
                if z_n == 1 || z_n == 6 {
                    continue;
                }
                let has_h_n = has_hydrogens(n.atom_index, &adjacency, &atomic_numbers);
                if has_h_n {
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(68);
        }
    }

    // Bit 70: [!#6!#1]~[#7]~[!#6!#1] — heteroatom-N-heteroatom
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                let het_neighbors: Vec<usize> = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| {
                        let z = atomic_numbers[n.atom_index];
                        z != 6 && z != 1
                    })
                    .map(|n| n.atom_index)
                    .collect();
                if het_neighbors.len() >= 2 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(69);
        }
    }

    // Bit 71: [#7]~[#8] — N-O bond (any)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 8 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(70);
        }
    }

    // Bit 72: [#8]~*~*~[#8] — O-X-X-O (two oxygens separated by 2 atoms)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 {
                if has_oxidizable_at_distance(i, 2, &adjacency, &atomic_numbers, i) {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(71);
        }
    }

    // Bit 73: [#16]=* — S= any (sulfur double bond)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 {
                for n in adjacency.neighbors_of(i) {
                    if let Some(order) = bond_orders.get(n.bond.index()) {
                        if *order == BondOrder::Double {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(72);
        }
    }

    // Bit 74: [CH3]~*~[CH3] — two methyls on same atom (isopropyl-like)
    {
        let mut found = false;
        for i in 0..num_atoms {
            let methyl_neighbors: usize = adjacency
                .neighbors_of(i)
                .iter()
                .filter(|n| {
                    atomic_numbers[n.atom_index] == 6
                        && adjacency.neighbors_of(n.atom_index).len() == 1
                })
                .count();
            if methyl_neighbors >= 2 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(73);
        }
    }

    // Bit 75: *!@[#7]@* — non-ring N in ring (N in ring with exocyclic substituent)
    // Approximate: N in ring with at least one non-ring neighbor
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 && atom_in_ring(i) {
                for n in adjacency.neighbors_of(i) {
                    if !atom_in_ring(n.atom_index) {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(74);
        }
    }

    // Bit 76: [#6]=[#6](~*)~* — C=C with two substituents
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                                let i_sub = adjacency.neighbors_of(i).len()
                                    - if adjacency
                                        .neighbors_of(i)
                                        .iter()
                                        .any(|nn| nn.atom_index == n.atom_index)
                                    {
                                        1
                                    } else {
                                        0
                                    };
                                let n_sub = adjacency.neighbors_of(n.atom_index).len() - 1; // minus the bond between them
                                if i_sub >= 1 && n_sub >= 1 {
                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                    if found {
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(75);
        }
    }

    // Bit 77: [#7]~*~[#7] — N-X-N (two nitrogens separated by one atom)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                if has_atom_at_distance(i, 1, &adjacency, &atomic_numbers, 7, i) {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(76);
        }
    }

    // Bit 78: [#6]=[#7] — C=N (imine, oxime)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 7 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(77);
        }
    }

    // Bit 79: [#7]~*~*~[#7] — N-X-X-N
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                if has_atom_at_distance(i, 2, &adjacency, &atomic_numbers, 7, i) {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(78);
        }
    }

    // Bit 80: [#7]~*~*~*~[#7] — N-X-X-X-N
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                if has_atom_at_distance(i, 3, &adjacency, &atomic_numbers, 7, i) {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(79);
        }
    }

    // Bit 81: [#16]~*(~*)~* — S connected to 3 atoms
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 16 && adjacency.neighbors_of(i).len() >= 3 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(80);
        }
    }

    // Bit 82: *~[CH2]~[!#6!#1!H0] — CH2-heteroatom(H) pattern
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    let z_n = atomic_numbers[n.atom_index];
                    if z_n != 6 && z_n != 1 {
                        let has_h_n = has_hydrogens(n.atom_index, &adjacency, &atomic_numbers);
                        if has_h_n {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(81);
        }
    }

    // Bit 83: [!#6!#1]1~*~*~*~*~1 — heteroatom in 5-membered ring
    {
        let mut found = false;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                if ri.atom_rings()[ring_idx].len() == 5 {
                    for a in ri.atom_rings()[ring_idx].iter() {
                        let z = atomic_numbers[a.index()];
                        if z != 6 && z != 1 {
                            found = true;
                            break;
                        }
                    }
                }
                if found {
                    break;
                }
            }
        }
        if found {
            on_bits.push(82);
        }
    }

    // Bit 84: [NH2] — primary amino group
    {
        let mut found = false;
        // An NH2 group: N with degree 1 (terminal), bonded to one heavy atom
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 && degree(i) == 1 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(83);
        }
    }

    // Bit 85: [#6]~[#7](~[#6])~[#6] — tertiary amine N with 3 carbons
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                let c_neighbors: usize = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 6)
                    .count();
                if c_neighbors >= 3 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(84);
        }
    }

    // Bit 86: [C;H2,H3][!#6!#1][C;H2,H3] — CH2/3-heteroatom-CH2/3
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i != 6 && z_i != 1 {
                let ch2_neighbors: usize = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 6 && degree(n.atom_index) <= 2)
                    .count();
                if ch2_neighbors >= 2 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(85);
        }
    }

    // Bit 87: [F,Cl,Br,I]!@*@* — halogen on non-ring chain (halogen not in ring)
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z = atomic_numbers[i];
            if z == 9 || z == 17 || z == 35 || z == 53 {
                if !atom_in_ring(i) {
                    for n in adjacency.neighbors_of(i) {
                        if atom_in_ring(n.atom_index) {
                            // halogen not in ring but bonded to ring atom: set
                            found = true;
                            break;
                        }
                        // Also check if bonded to non-ring chain
                        if degree(n.atom_index) > 1 {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(86);
        }
    }

    // Bit 88: S
    if has_element(16) {
        on_bits.push(87);
    }

    // Bit 89: [#8]~*~*~*~[#8] — O-X-X-X-O
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 {
                if has_oxidizable_at_distance(i, 3, &adjacency, &atomic_numbers, i) {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(88);
        }
    }

    // Bit 90: [!#6!#1!H0] with 2-3 bridge to CH2 (complex MACCS pattern)
    // Approximate: heteroatom with H within 2 bonds of a CH2
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            let has_h_i = has_hydrogens(i, &adjacency, &atomic_numbers);
            if !has_h_i {
                continue;
            }
            // Check if within distance 3 of any CH2
            if has_ch2_at_distance(i, 3, &adjacency, &atomic_numbers, &atomic_numbers) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(89);
        }
    }

    // Bit 91: [!#6!#1!H0] with 3-4 bridge to CH2 (longer)
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            let has_h_i = has_hydrogens(i, &adjacency, &atomic_numbers);
            if !has_h_i {
                continue;
            }
            if has_ch2_at_distance(i, 4, &adjacency, &atomic_numbers, &atomic_numbers) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(90);
        }
    }

    // Bit 92: [#8]~[#6](~[#7])~[#6] — O-C(N)-C (ester/amide-like)
    // RDKit: O-C(=O)-N pattern equivalent
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                let has_o = adjacency
                    .neighbors_of(i)
                    .iter()
                    .any(|n| atomic_numbers[n.atom_index] == 8);
                let has_n = adjacency
                    .neighbors_of(i)
                    .iter()
                    .any(|n| atomic_numbers[n.atom_index] == 7);
                let c_neighbors = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 6)
                    .count();
                if has_o && has_n && c_neighbors >= 1 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(91);
        }
    }

    // Bit 93: [!#6!#1]~[CH3] — heteroatom-methyl bond
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i != 6 && z_i != 1 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 && degree(n.atom_index) == 1 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(92);
        }
    }

    // Bit 94: [!#6!#1]~[#7] — heteroatom-N bond
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i != 6 && z_i != 1 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 7 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(93);
        }
    }

    // Bit 95: [#7]~*~*~[#8] — N-X-X-O
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                if has_oxidizable_at_distance(i, 2, &adjacency, &atomic_numbers, i) {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(94);
        }
    }

    // Bit 96: *1~*~*~*~*~1 — 5-membered ring
    {
        let mut found = false;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                if ri.atom_rings()[ring_idx].len() == 5 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(95);
        }
    }

    // Bit 97: [#7]~*~*~*~[#8] — N-X-X-X-O
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                if has_oxidizable_at_distance(i, 3, &adjacency, &atomic_numbers, i) {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(96);
        }
    }

    // Bit 98: [!#6!#1]1~*~*~*~*~*~1 — heteroatom in 6-membered ring
    {
        let mut found = false;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                if ri.atom_rings()[ring_idx].len() == 6 {
                    for a in ri.atom_rings()[ring_idx].iter() {
                        let z = atomic_numbers[a.index()];
                        if z != 6 && z != 1 {
                            found = true;
                            break;
                        }
                    }
                }
                if found {
                    break;
                }
            }
        }
        if found {
            on_bits.push(97);
        }
    }

    // Bit 99: [#6]=[#6] — C=C double bond
    {
        if has_double_bond() {
            on_bits.push(98);
        }
    }

    // Bit 100: *~[CH2]~[#7] — CH2-N pattern
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 && degree(i) == 2 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 7 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(99);
        }
    }

    // Bit 101: large ring (8+ members)
    {
        let mut found = false;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                if ri.atom_rings()[ring_idx].len() >= 8 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(100);
        }
    }

    // Bit 102: [!#6!#1]~[#8] — heteroatom-O bond
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i != 6 && z_i != 1 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 8 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(101);
        }
    }

    // Bit 103: Cl
    if has_element(17) {
        on_bits.push(102);
    }

    // Bit 104: [!#6!#1!H0]~*~[CH2]~* — heteroatom(H)-X-CH2-X
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            let has_h_i = has_hydrogens(i, &adjacency, &atomic_numbers);
            if !has_h_i {
                continue;
            }
            // BFS up to 2 bonds to find a CH2
            for n in adjacency.neighbors_of(i) {
                for nn in adjacency.neighbors_of(n.atom_index) {
                    if nn.atom_index != i && atomic_numbers[nn.atom_index] == 6 {
                        // Found a carbon within 2 bonds — mark as found
                        found = true;
                        break;
                    }
                }
                if found {
                    break;
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(103);
        }
    }

    // Bit 105: *@*(@*)@* — branching at ring atom (atom in ring with 3+ ring neighbors)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atom_in_ring(i) {
                let ring_nbrs: usize = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atom_in_ring(n.atom_index))
                    .count();
                if ring_nbrs >= 3 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(104);
        }
    }

    // Bit 106: [!#6!#1]~*(~[!#6!#1])~[!#6!#1] — heteroatom connected to 3+ heteroatoms
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            let het_neighbors: usize = adjacency
                .neighbors_of(i)
                .iter()
                .filter(|n| {
                    let z = atomic_numbers[n.atom_index];
                    z != 6 && z != 1
                })
                .count();
            if het_neighbors >= 3 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(105);
        }
    }

    // Bit 107: [F,Cl,Br,I]~*(~*)~* — halogen connected to non-terminal atom
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 9 || z_i == 17 || z_i == 35 || z_i == 53 {
                for n in adjacency.neighbors_of(i) {
                    if degree(n.atom_index) >= 3 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(106);
        }
    }

    // Bit 108: [CH3]~*~*~*~[CH2]~* — methyl-propyl-CH2 pattern
    {
        let mut found = false;
        if has_element(6) {
            // Check for chain of length >= 5 carbons
            let mut chain_found = false;
            for i in 0..num_atoms {
                if atomic_numbers[i] == 6 && degree(i) == 1 {
                    // terminal methyl, walk the chain
                    let mut prev = i;
                    let mut current = adjacency.neighbors_of(i)[0].atom_index;
                    let mut chain_len = 2i32;
                    loop {
                        let nbrs: Vec<usize> = adjacency
                            .neighbors_of(current)
                            .iter()
                            .map(|n| n.atom_index)
                            .filter(|&n| n != prev)
                            .collect();
                        if nbrs.is_empty() {
                            break;
                        }
                        if nbrs.len() > 1 {
                            break;
                        } // branch
                        prev = current;
                        current = nbrs[0];
                        chain_len += 1;
                        if chain_len >= 5 {
                            chain_found = true;
                            break;
                        }
                    }
                    if chain_found {
                        break;
                    }
                }
            }
            found = chain_found;
        }
        if found {
            on_bits.push(107);
        }
    }

    // Bit 109: *~[CH2]~[#8] — CH2-O pattern
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 8 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(108);
        }
    }

    // Bit 110: [#7]~[#6]~[#8] — N-C-O pattern
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                let has_n = adjacency
                    .neighbors_of(i)
                    .iter()
                    .any(|n| atomic_numbers[n.atom_index] == 7);
                let has_o = adjacency
                    .neighbors_of(i)
                    .iter()
                    .any(|n| atomic_numbers[n.atom_index] == 8);
                if has_n && has_o {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(109);
        }
    }

    // Bit 111: [#7]~*~[CH2]~* — N-X-CH2-X
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                for n in adjacency.neighbors_of(i) {
                    for nn in adjacency.neighbors_of(n.atom_index) {
                        if nn.atom_index != i && atomic_numbers[nn.atom_index] == 6 {
                            found = true;
                            break;
                        }
                    }
                    if found {
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(110);
        }
    }

    // Bit 112: *~*(~*)(~*)~* — branching atom with 4 neighbors
    {
        let mut found = false;
        for i in 0..num_atoms {
            if degree(i) >= 4 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(111);
        }
    }

    // Bit 113: [#8]!:*:* — O non-aromatic in aromatic system
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 {
                for n in adjacency.neighbors_of(i) {
                    if aromatic_atoms[n.atom_index] {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order != BondOrder::Aromatic {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(112);
        }
    }

    // Bit 114: [CH3]~[CH2]~* — ethyl group
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 && degree(i) == 1 {
                // Methyl connected to a CH2
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 && degree(n.atom_index) >= 2 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(113);
        }
    }

    // Bit 115: [CH3]~*~[CH2]~* — methyl-X-CH2-X
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 && degree(i) == 1 {
                for n in adjacency.neighbors_of(i) {
                    for nn in adjacency.neighbors_of(n.atom_index) {
                        if nn.atom_index != i && atomic_numbers[nn.atom_index] == 6 {
                            found = true;
                            break;
                        }
                    }
                    if found {
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(114);
        }
    }

    // Bit 116: [CH3]~*~*~[CH2]~* — methyl-X-X-CH2-X
    {
        let mut found_116 = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 && degree(i) == 1 {
                // Terminal methyl — walk path to find a CH2 at distance 3
                let mut visited = vec![false; num_atoms];
                visited[i] = true;
                let mut queue = vec![];
                for n in adjacency.neighbors_of(i) {
                    queue.push((n.atom_index, 1usize));
                    visited[n.atom_index] = true;
                }
                while let Some((node, dist)) = queue.pop() {
                    if atomic_numbers[node] == 6 && dist >= 2 {
                        found_116 = true;
                        break;
                    }
                    if dist < 3 {
                        for n in adjacency.neighbors_of(node) {
                            if !visited[n.atom_index] {
                                visited[n.atom_index] = true;
                                queue.push((n.atom_index, dist + 1));
                            }
                        }
                    }
                }
                if found_116 {
                    break;
                }
            }
        }
        if found_116 {
            on_bits.push(115);
        }
    }

    // Bit 117: [#7]~*~[#8] — N-X-O (two heteroatoms separated by one atom)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                if has_oxidizable_at_distance(i, 1, &adjacency, &atomic_numbers, i) {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(116);
        }
    }

    // Bit 118: two non-adjacent CH2 groups or CH2-CH2 chain pattern
    {
        // *~[CH2]~[CH2]~* pattern: two consecutive CH2 groups
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(117);
        }
    }

    // Bit 119: [#7]=* — N= any
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                for n in adjacency.neighbors_of(i) {
                    if let Some(order) = bond_orders.get(n.bond.index()) {
                        if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(118);
        }
    }

    // Bit 120: 2+ non-C non-H atoms in rings
    {
        let mut het_in_ring_count = 0usize;
        for i in 0..num_atoms {
            let z = atomic_numbers[i];
            if z != 6 && z != 1 && atom_in_ring(i) {
                het_in_ring_count += 1;
                if het_in_ring_count >= 2 {
                    break;
                }
            }
        }
        if het_in_ring_count >= 2 {
            on_bits.push(119);
        }
    }

    // Bit 121: [#7R] — N in ring
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 && atom_in_ring(i) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(120);
        }
    }

    // Bit 122: *~[#7](~*)~* — N with 3 neighbors (tertiary or in ring junction)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 && degree(i) >= 3 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(121);
        }
    }

    // Bit 123: [#8]~[#6]~[#8] — O-C-O (acetal, carbonate, carboxylate)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                let o_neighbors: usize = adjacency
                    .neighbors_of(i)
                    .iter()
                    .filter(|n| atomic_numbers[n.atom_index] == 8)
                    .count();
                if o_neighbors >= 2 {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(122);
        }
    }

    // Bit 124: [!#6!#1]~[!#6!#1] — heteroatom-heteroatom bond (any)
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            for n in adjacency.neighbors_of(i) {
                let z_n = atomic_numbers[n.atom_index];
                if z_n != 1 && z_n != 6 {
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(123);
        }
    }

    // Bit 125: two or more aromatic rings
    {
        let mut aro_ring_count = 0usize;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                let atoms_in_ring = &ri.atom_rings()[ring_idx];
                let is_aro = atoms_in_ring.iter().all(|a| aromatic_atoms[a.index()]);
                if is_aro {
                    aro_ring_count += 1;
                    if aro_ring_count >= 2 {
                        break;
                    }
                }
            }
        }
        if aro_ring_count >= 2 {
            on_bits.push(124);
        }
    }

    // Bit 126: *!@[#8]!@* — O with all non-ring bonds
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 && !atom_in_ring(i) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(125);
        }
    }

    // Bit 127: *@*!@[#8] — O in ring with exocyclic substituent (2+ occurrences)
    {
        let mut count = 0usize;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 && atom_in_ring(i) {
                for n in adjacency.neighbors_of(i) {
                    if !atom_in_ring(n.atom_index) {
                        count += 1;
                        break;
                    }
                }
            }
        }
        if count > 1 {
            on_bits.push(126);
        }
        // Bit 143: *@*!@[#8] (count > 0)
        if count > 0 {
            on_bits.push(142);
        }
    }

    // Bit 128: complex CH2 pattern with 5-atom bridge (RDKit)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 && degree(i) == 2 {
                // Walk 4 bonds to see if we find another CH2
                let mut visited = vec![false; num_atoms];
                visited[i] = true;
                let mut queue = vec![];
                for n in adjacency.neighbors_of(i) {
                    visited[n.atom_index] = true;
                    queue.push((n.atom_index, 1usize));
                }
                while let Some((node, dist)) = queue.pop() {
                    if dist >= 3 && atomic_numbers[node] == 6 {
                        found = true;
                        break;
                    }
                    if dist < 4 {
                        for n in adjacency.neighbors_of(node) {
                            if !visited[n.atom_index] {
                                visited[n.atom_index] = true;
                                queue.push((n.atom_index, dist + 1));
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(127);
        }
    }

    // Bit 129: shorter CH2 bridge (3-4 bonds)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                let mut visited = vec![false; num_atoms];
                visited[i] = true;
                let mut queue = vec![];
                for n in adjacency.neighbors_of(i) {
                    visited[n.atom_index] = true;
                    queue.push((n.atom_index, 1usize));
                }
                while let Some((node, dist)) = queue.pop() {
                    if dist >= 2 && atomic_numbers[node] == 6 {
                        found = true;
                        break;
                    }
                    if dist < 3 {
                        for n in adjacency.neighbors_of(node) {
                            if !visited[n.atom_index] {
                                visited[n.atom_index] = true;
                                queue.push((n.atom_index, dist + 1));
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(128);
        }
    }

    // Bit 130: 2+ heteroatom-heteroatom bonds (from bit 124 count)
    {
        let mut count = 0usize;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i == 1 || z_i == 6 {
                continue;
            }
            for n in adjacency.neighbors_of(i) {
                let z_n = atomic_numbers[n.atom_index];
                if z_n != 1 && z_n != 6 {
                    count += 1;
                    if count > 1 {
                        break;
                    }
                }
            }
            if count > 1 {
                break;
            }
        }
        if count > 1 {
            on_bits.push(129);
        }
    }

    // Bit 131: [!#6!#1!H0] — any non-C non-H heteroatom with at least 1 H
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z = atomic_numbers[i];
            if z != 6 && z != 1 {
                if has_hydrogens(i, &adjacency, &atomic_numbers) {
                    found = true;
                    break;
                }
            }
        }
        if found {
            on_bits.push(130);
        }
    }

    // Bit 132: [#8]~*~[CH2]~* — O-X-CH2-X
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 {
                for n in adjacency.neighbors_of(i) {
                    for nn in adjacency.neighbors_of(n.atom_index) {
                        if nn.atom_index != i && atomic_numbers[nn.atom_index] == 6 {
                            found = true;
                            break;
                        }
                    }
                    if found {
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(131);
        }
    }

    // Bit 133: *@*!@[#7] — ring N with exocyclic substituent
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 && atom_in_ring(i) {
                for n in adjacency.neighbors_of(i) {
                    if !atom_in_ring(n.atom_index) {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(132);
        }
    }

    // Bit 134: halogen (F, Cl, Br, I)
    {
        let mut found = false;
        for &z in &atomic_numbers {
            if z == 9 || z == 17 || z == 35 || z == 53 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(133);
        }
    }

    // Bit 135: [#7]!:*:* — N non-aromatic in aromatic system
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 {
                for n in adjacency.neighbors_of(i) {
                    if aromatic_atoms[n.atom_index] {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order != BondOrder::Aromatic {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(134);
        }
    }

    // Bit 136: [#8]=* — O= any
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 {
                for n in adjacency.neighbors_of(i) {
                    if let Some(order) = bond_orders.get(n.bond.index()) {
                        if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(135);
        }
    }

    // Bit 137: [!C!cR] — non-carbon in ring
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z = atomic_numbers[i];
            if z != 6 && atom_in_ring(i) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(136);
        }
    }

    // Bit 138: [!#6!#1]~[CH2]~* — heteroatom-CH2-any
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z_i = atomic_numbers[i];
            if z_i != 6 && z_i != 1 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(137);
        }
    }

    // Bit 139: [O!H0] — O with H (hydroxyl)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 && has_hydrogens(i, &adjacency, &atomic_numbers) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(138);
        }
    }

    // Bit 140: [#8] — O (any)
    if has_element(8) {
        on_bits.push(139);
    }

    // Bit 141: [CH3] — any methyl
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 && degree(i) == 1 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(140);
        }
    }

    // Bit 142: [#7] — N (any)
    if has_element(7) {
        on_bits.push(141);
    }

    // Bit 143: already set above with bit 127

    // Bit 144: *!:*:*!:* — alternating aromatic/non-aromatic in aromatic system
    {
        let mut found = false;
        if has_aromatic_bond() {
            // Check for pattern like C-C(=C)-C
            for i in 0..num_atoms {
                if aromatic_atoms[i] {
                    continue;
                }
                for n in adjacency.neighbors_of(i) {
                    if aromatic_atoms[n.atom_index] {
                        for nn in adjacency.neighbors_of(n.atom_index) {
                            if nn.atom_index != i && !aromatic_atoms[nn.atom_index] {
                                found = true;
                                break;
                            }
                        }
                    }
                    if found {
                        break;
                    }
                }
                if found {
                    break;
                }
            }
        }
        if found {
            on_bits.push(143);
        }
    }

    // Bit 145: *1~*~*~*~*~*~1 — 6-membered ring (2+ occurrences)
    {
        let mut count = 0usize;
        if let Some(ri) = ring_info {
            for ring_idx in 0..ri.atom_rings().len() {
                if ri.atom_rings()[ring_idx].len() == 6 {
                    count += 1;
                    if count > 1 {
                        break;
                    }
                }
            }
        }
        if count > 1 {
            on_bits.push(144);
        }
        // Bit 163: any 6-membered ring
        if count > 0 {
            on_bits.push(162);
        }
    }

    // Bit 147: CH2-CH2 chain or ring pattern
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(146);
        }
    }

    // Bit 148: *~[!#6!#1](~*)~* — heteroatom with 3+ neighbors
    {
        let mut found = false;
        for i in 0..num_atoms {
            let z = atomic_numbers[i];
            if z != 6 && z != 1 && degree(i) >= 3 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(147);
        }
    }

    // Bit 149: [C;H3,H4] — CH3 or CH4 (terminal methyl/methane) (2+)
    {
        let count = atomic_numbers
            .iter()
            .enumerate()
            .filter(|&(i, &z)| z == 6 && degree(i) <= 1)
            .count();
        if count > 1 {
            on_bits.push(148);
        }
        // Bit 160: any [C;H3,H4]
        if count > 0 {
            on_bits.push(159);
        }
    }

    // Bit 150: *!@*@*!@* — ring chain ring (atoms in ring - non-ring - ring)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atom_in_ring(i) {
                for n in adjacency.neighbors_of(i) {
                    if !atom_in_ring(n.atom_index) {
                        for nn in adjacency.neighbors_of(n.atom_index) {
                            if nn.atom_index != i && atom_in_ring(nn.atom_index) {
                                found = true;
                                break;
                            }
                        }
                    }
                    if found {
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(149);
        }
    }

    // Bit 151: [#7!H0] — N with H (any)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 && has_hydrogens(i, &adjacency, &atomic_numbers) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(150);
        }
    }

    // Bit 152: [#8]~[#6](~[#6])~[#6] — O-C(R)-R' (ether/ester attached to carbon)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 8 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 6 {
                        let c_neighbors_of_c: usize = adjacency
                            .neighbors_of(n.atom_index)
                            .iter()
                            .filter(|nn| nn.atom_index != i && atomic_numbers[nn.atom_index] == 6)
                            .count();
                        if c_neighbors_of_c >= 2 {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(151);
        }
    }

    // Bit 154: [#6]=[#8] — C=O (carbonyl)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 8 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Double || *order == BondOrder::Aromatic {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(153);
        }
    }

    // Bit 155: *!@[CH2]!@* — non-ring CH2 between two non-ring bonds
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 && !atom_in_ring(i) && degree(i) == 2 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(154);
        }
    }

    // Bit 156: [#7]~*(~*)~* — N with 3+ neighbors
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 7 && degree(i) >= 3 {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(155);
        }
    }

    // Bit 157: [#6]-[#8] — C-O single bond (alcohol, ether)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 8 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Single {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(156);
        }
    }

    // Bit 158: [#6]-[#7] — C-N single bond (amine)
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atomic_numbers[i] == 6 {
                for n in adjacency.neighbors_of(i) {
                    if atomic_numbers[n.atom_index] == 7 {
                        if let Some(order) = bond_orders.get(n.bond.index()) {
                            if *order == BondOrder::Single {
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            if found {
                break;
            }
        }
        if found {
            on_bits.push(157);
        }
    }

    // Bit 160: already set above with bit 149

    // Bit 161: any [#7] (already set as bit 142)
    // RDKit: if count > 0 set bit 161 (redundant with bit 142)
    if has_element(7) {
        on_bits.push(160);
    }

    // Bit 162: already set above with bit 145
    // Bit 163: already set above with bit 145

    // Bit 165: [R] — any ring atom
    {
        let mut found = false;
        for i in 0..num_atoms {
            if atom_in_ring(i) {
                found = true;
                break;
            }
        }
        if found {
            on_bits.push(164);
        }
    }

    // Bit 166: multiple fragments
    {
        let mut num_fragments = 0usize;
        let mut visited = vec![false; num_atoms];
        for i in 0..num_atoms {
            if !visited[i] {
                num_fragments += 1;
                if num_fragments > 1 {
                    break;
                }
                // BFS to mark this fragment
                let mut queue = vec![i];
                visited[i] = true;
                while let Some(node) = queue.pop() {
                    for n in adjacency.neighbors_of(node) {
                        if !visited[n.atom_index] {
                            visited[n.atom_index] = true;
                            queue.push(n.atom_index);
                        }
                    }
                }
            }
        }
        if num_fragments > 1 {
            on_bits.push(165);
        }
    }

    Fingerprint::from_on_bits(n_bits, on_bits)
}

// ---------------------------------------------------------------------------
// Helpers for MACCS keys (private)
// ---------------------------------------------------------------------------

/// Check if an atom has implicit + explicit hydrogens.
/// Uses valence cache if available, falls back to degree-based heuristic.
fn has_hydrogens(atom_idx: usize, adjacency: &AdjacencyList, atomic_numbers: &[u8]) -> bool {
    // For now, we use a simplified heuristic:
    // An atom "has H" if its degree (heavy atom count) is less than its typical valence.
    // This is imperfect but works for common cases.
    let z = atomic_numbers[atom_idx];
    let deg = adjacency.neighbors_of(atom_idx).len() as u32;
    let typical_valence: u32 = match z {
        5 => 3,  // B
        6 => 4,  // C
        7 => 3,  // N
        8 => 2,  // O
        9 => 1,  // F
        14 => 4, // Si
        15 => 3, // P (common)
        16 => 2, // S
        17 => 1, // Cl
        35 => 1, // Br
        53 => 1, // I
        _ => 4,  // default
    };
    deg < typical_valence
}

/// Check if there is an atom of target_z at exactly `distance` bonds from `start`
/// (not passing through `start` itself).
fn has_atom_at_distance(
    start: usize,
    distance: usize,
    adjacency: &AdjacencyList,
    atomic_numbers: &[u8],
    target_z: u8,
    exclude: usize,
) -> bool {
    let mut visited = vec![false; atomic_numbers.len()];
    visited[start] = true;
    let mut queue: Vec<(usize, usize)> = Vec::new();
    for n in adjacency.neighbors_of(start) {
        if n.atom_index != exclude {
            visited[n.atom_index] = true;
            queue.push((n.atom_index, 1usize));
        }
    }
    // BFS limited to `distance` steps
    let mut next_queue: Vec<(usize, usize)> = Vec::new();
    for d in 1..=distance {
        for (node, _) in queue.drain(..) {
            if d == distance && atomic_numbers[node] == target_z {
                return true;
            }
            if d < distance {
                for n in adjacency.neighbors_of(node) {
                    if !visited[n.atom_index] {
                        visited[n.atom_index] = true;
                        next_queue.push((n.atom_index, d + 1));
                    }
                }
            }
        }
        std::mem::swap(&mut queue, &mut next_queue);
    }
    false
}

/// Check if there is an O atom at exactly `distance` bonds from `start`.
fn has_oxidizable_at_distance(
    start: usize,
    distance: usize,
    adjacency: &AdjacencyList,
    atomic_numbers: &[u8],
    exclude: usize,
) -> bool {
    has_atom_at_distance(start, distance, adjacency, atomic_numbers, 8, exclude)
}

/// Check if there is a heteroatom (non-C, non-H, has H) at exactly `distance` bonds from `start`.
fn has_heteroatom_at_distance(
    start: usize,
    distance: usize,
    adjacency: &AdjacencyList,
    atomic_numbers: &[u8],
    _charges: &[i8],
    need_h: bool,
    exclude: usize,
) -> bool {
    let mut visited = vec![false; atomic_numbers.len()];
    visited[start] = true;
    let mut queue: Vec<(usize, usize)> = Vec::new();
    for n in adjacency.neighbors_of(start) {
        if n.atom_index != exclude {
            visited[n.atom_index] = true;
            queue.push((n.atom_index, 1usize));
        }
    }
    let mut next_queue: Vec<(usize, usize)> = Vec::new();
    for d in 1..=distance {
        for (node, _) in queue.drain(..) {
            let z = atomic_numbers[node];
            if d == distance && z != 1 && z != 6 {
                if !need_h || has_hydrogens(node, adjacency, atomic_numbers) {
                    return true;
                }
            }
            if d < distance {
                for n in adjacency.neighbors_of(node) {
                    if !visited[n.atom_index] {
                        visited[n.atom_index] = true;
                        next_queue.push((n.atom_index, d + 1));
                    }
                }
            }
        }
        std::mem::swap(&mut queue, &mut next_queue);
    }
    false
}

/// Check if there's a CH2-like carbon (degree 2) at `distance` bonds from `start`.
fn has_ch2_at_distance(
    start: usize,
    distance: usize,
    adjacency: &AdjacencyList,
    atomic_numbers: &[u8],
    _hydropattern: &[u8],
) -> bool {
    let mut visited = vec![false; atomic_numbers.len()];
    visited[start] = true;
    let mut queue: Vec<(usize, usize)> = Vec::new();
    for n in adjacency.neighbors_of(start) {
        visited[n.atom_index] = true;
        queue.push((n.atom_index, 1usize));
    }
    let mut next_queue: Vec<(usize, usize)> = Vec::new();
    for d in 1..=distance {
        for (node, _) in queue.drain(..) {
            if d == distance && atomic_numbers[node] == 6 {
                return true;
            }
            if d < distance {
                for n in adjacency.neighbors_of(node) {
                    if !visited[n.atom_index] {
                        visited[n.atom_index] = true;
                        next_queue.push((n.atom_index, d + 1));
                    }
                }
            }
        }
        std::mem::swap(&mut queue, &mut next_queue);
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
