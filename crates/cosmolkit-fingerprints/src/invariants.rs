//! Source-ordered Morgan atom invariant construction over prepared detached
//! topology, valence, and ring values.

use std::sync::OnceLock;

use crate::{FingerprintError, MorganError, hash::hash_range};
use cosmolkit_core::{
    PeriodicTableError, RingInfo, ValenceAssignment, atomic_mass,
    total_hydrogen_count_from_validated,
};
use cosmolkit_model::{BondOrder, BondStereo, TopologyBlock};
use cosmolkit_search::{
    QueryGraph, QueryMatchContext, SearchTarget, SearchTargetAccess, SmartsParseError,
    SmartsParseParams, SubstructMatchParams, parse_smarts,
    try_get_substruct_matches_with_params_and_context,
};

/// Fill the source-sized connectivity-invariant prefix without changing any
/// caller-provided suffix values. The caller supplies a validated topology
/// with its matching prepared valence and ring assignments; this private
/// boundary deliberately avoids repeating whole-topology validation per atom.
pub(crate) fn get_connectivity_invariants(
    topology: &TopologyBlock,
    invars: &mut [u32],
    include_ring_membership: bool,
    valence: &ValenceAssignment,
    rings: &RingInfo,
) -> Result<(), MorganError> {
    // BEGIN RDKIT CPP FUNCTION MorganFingerprints::getConnectivityInvariants
    // RDKit❗🔝: void getConnectivityInvariants(const ROMol &mol, std::vector<uint32_t> &invars,
    // RDKit❗🔝:                                bool includeRingMembership) {
    // RDKit❗🔝:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗🔝:   PRECONDITION(invars.size() >= nAtoms, "vector too small");
    // RDKit❗🔝:   gboost::hash<std::vector<uint32_t>> vectHasher;
    // RDKit❗🔝:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit❗🔝:     Atom const *atom = mol.getAtomWithIdx(i);
    // RDKit❗🔝:     std::vector<uint32_t> components;
    // RDKit❗🔝:     components.push_back(atom->getAtomicNum());
    // RDKit❗🔝:     components.push_back(atom->getTotalDegree());
    // RDKit❗🔝:     components.push_back(atom->getTotalNumHs(true));
    // RDKit❗🔝:     components.push_back(atom->getFormalCharge());
    // RDKit❗🔝:     int deltaMass = static_cast<int>(
    // RDKit❗🔝:         atom->getMass() -
    // RDKit❗🔝:         PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
    // RDKit❗🔝:     components.push_back(deltaMass);
    // RDKit❗🔝:
    // RDKit❗🔝:     if (includeRingMembership &&
    // RDKit❗🔝:         atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx())) {
    // RDKit❗🔝:       components.push_back(1);
    // RDKit❗🔝:     }
    // RDKit❗🔝:     invars[i] = vectHasher(components);
    // RDKit❗🔝:   }
    // RDKit❗🔝: }  // end of getConnectivityInvariants()
    // END RDKIT CPP FUNCTION MorganFingerprints::getConnectivityInvariants

    // BEGIN RDKIT CPP FUNCTION Atom::getTotalDegree / Atom::getDegree / ROMol::getAtomDegree
    // RDKit✔️✔️: unsigned int Atom::getTotalDegree() const {
    // RDKit✔️✔️:   unsigned int res = this->getTotalNumHs(false) + this->getDegree();
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int Atom::getDegree() const {
    // RDKit✔️✔️:   return dp_mol ? getOwningMol().getAtomDegree(this) : 0;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int ROMol::getAtomDegree(const Atom *at) const {
    // RDKit✔️✔️:   PRECONDITION(at, "no atom");
    // RDKit✔️✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit✔️✔️:                "atom not associated with this molecule");
    // RDKit✔️✔️:   return rdcast<unsigned int>(boost::out_degree(at->getIdx(), d_graph));
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION Atom::getTotalDegree / Atom::getDegree / ROMol::getAtomDegree

    // BEGIN RDKIT CPP FUNCTION Atom::getTotalNumHs / Atom::getNumImplicitHs / Atom::getValence
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
    // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int Atom::getNumImplicitHs() const {
    // RDKit✔️✔️:   if (df_noImplicit) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   PRECONDITION(d_implicitValence > -1,
    // RDKit✔️✔️:                "getNumImplicitHs() called without preceding call to "
    // RDKit✔️✔️:                "calcImplicitValence()");
    // RDKit✔️✔️:   return getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
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
    // END RDKIT CPP FUNCTION Atom::getTotalNumHs / Atom::getNumImplicitHs / Atom::getValence

    // BEGIN RDKIT CPP FUNCTION ROMol::atomNeighbors
    // RDKit✔️✔️: CXXAtomIterator<const MolGraph, Atom *const, MolGraph::adjacency_iterator>
    // RDKit✔️✔️: atomNeighbors(Atom const *at) const {
    // RDKit✔️✔️:   auto pr = getAtomNeighbors(at);
    // RDKit✔️✔️:   return {&d_graph, pr.first, pr.second};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ROMol::atomNeighbors

    // BEGIN RDKIT CPP FUNCTION Atom scalar accessors
    // RDKit✔️✔️: int getAtomicNum() const { return d_atomicNum; }
    // RDKit✔️✔️: int getFormalCharge() const { return d_formalCharge; }
    // RDKit✔️✔️: unsigned int getNumExplicitHs() const { return d_numExplicitHs; }
    // RDKit✔️✔️: bool getNoImplicit() const { return df_noImplicit; }
    // RDKit✔️✔️: unsigned int getIsotope() const { return d_isotope; }
    // END RDKIT CPP FUNCTION Atom scalar accessors

    // BEGIN RDKIT CPP FUNCTION Atom::getMass / PeriodicTable mass lookups
    // RDKit✔️✔️: double Atom::getMass() const {
    // RDKit✔️✔️:   if (d_isotope) {
    // RDKit✔️✔️:     double res =
    // RDKit✔️✔️:         PeriodicTable::getTable()->getMassForIsotope(d_atomicNum, d_isotope);
    // RDKit✔️✔️:     if (d_atomicNum != 0 && res == 0.0) {
    // RDKit✔️✔️:       res = d_isotope;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return PeriodicTable::getTable()->getAtomicWeight(d_atomicNum);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double getAtomicWeight(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   double mass = byanum[atomicNumber].Mass();
    // RDKit✔️✔️:   return mass;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double getMassForIsotope(UINT atomicNumber, UINT isotope) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   const std::map<unsigned int, std::pair<double, double>> &m =
    // RDKit✔️✔️:       byanum[atomicNumber].d_isotopeInfoMap;
    // RDKit✔️✔️:   std::map<unsigned int, std::pair<double, double>>::const_iterator item =
    // RDKit✔️✔️:       m.find(isotope);
    // RDKit✔️✔️:   if (item == m.end()) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return item->second.first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getMass / PeriodicTable mass lookups

    // BEGIN RDKIT CPP FUNCTION RingInfo::numAtomRings
    // RDKit✔️✔️: unsigned int RingInfo::numAtomRings(unsigned int idx) const {
    // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (idx < d_atomMembers.size()) {
    // RDKit✔️✔️:     return rdcast<unsigned int>(d_atomMembers[idx].size());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RingInfo::numAtomRings

    // BEGIN RDKIT CPP FUNCTION gboost::hash<vector<uint32_t>> / hash_range / hash_combine
    // RDKit✔️✔️: typedef std::uint32_t hash_result_t;
    // RDKit✔️✔️: std::hash_result_t operator()(T const& val) const { return hash_value(val); }
    // RDKit✔️✔️: std::hash_result_t hash_value(std::vector<T, A> const& v) {
    // RDKit✔️✔️:   return hash_range(v.begin(), v.end());
    // RDKit✔️✔️: }
    // RDKit✔️✔️: inline std::hash_result_t hash_range(It first, It last) {
    // RDKit✔️✔️:   std::hash_result_t seed = 0;
    // RDKit✔️✔️:   for (; first != last; ++first) {
    // RDKit✔️✔️:     hash_combine(seed, *first);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return seed;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: inline void hash_combine(std::hash_result_t& seed, T const& v) {
    // RDKit✔️✔️:   gboost::hash<T> hasher;
    // RDKit✔️✔️:   seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: inline std::hash_result_t hash_value(unsigned int v) {
    // RDKit✔️✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION gboost::hash<vector<uint32_t>> / hash_range / hash_combine

    let num_atoms = topology.atoms.len();
    if invars.len() < num_atoms {
        return Err(FingerprintError::PreconditionViolation {
            what: "vector too small",
        }
        .into());
    }

    for (index, atom) in topology.atoms.iter().enumerate() {
        // Atom::getTotalDegree() adds getTotalNumHs(false) to the owning
        // molecule's graph degree, then getTotalNumHs(true) is evaluated
        // independently for the next source component.
        let hydrogens_for_degree =
            total_hydrogen_count_from_validated(topology, valence, atom.id(), false)?;
        let degree = topology.adjacency.neighbors_of(atom.id().index()).len() as u32;
        let total_degree = hydrogens_for_degree.wrapping_add(degree);
        let total_hydrogens =
            total_hydrogen_count_from_validated(topology, valence, atom.id(), true)?;

        let atomic_weight = atomic_mass(atom.element(), None)?;
        let mass = match atom.isotope() {
            None => atomic_weight,
            Some(isotope) => match atomic_mass(atom.element(), Some(isotope)) {
                Ok(mass) => mass,
                // Atom::getMass() substitutes the mass number only for a
                // non-dummy atom when the source isotope table returned zero.
                Err(PeriodicTableError::UnknownIsotope { .. }) if atom.atomic_number() != 0 => {
                    f64::from(isotope)
                }
                Err(PeriodicTableError::UnknownIsotope { .. }) => 0.0,
                Err(source) => return Err(MorganError::PeriodicTable(source)),
            },
        };
        let delta_mass = (mass - atomic_weight) as i32;

        // A fixed stack array carries the same ordered five/six components
        // while avoiding RDKit's per-atom dynamic-vector allocation.
        let mut components = [0_u32; 6];
        components[0] = u32::from(atom.atomic_number());
        components[1] = total_degree;
        components[2] = total_hydrogens;
        components[3] = i32::from(atom.formal_charge()) as u32;
        components[4] = delta_mass as u32;
        let mut component_count = 5;

        if include_ring_membership {
            if !rings.is_initialized() {
                return Err(FingerprintError::PreconditionViolation {
                    what: "RingInfo not initialized",
                }
                .into());
            }
            if rings.num_atom_rings(atom.id()) != 0 {
                components[5] = 1;
                component_count = 6;
            }
        }

        invars[index] = hash_range(&components[..component_count]);
    }

    // Behavior review: the buffer-size precondition precedes every write;
    // each source atom contributes the exact ordered u32 component sequence,
    // and the optional ring marker is present only for actual ring atoms.
    // The call receives source-prepared valence/rings; chemistry is not
    // recomputed here. The stack array changes allocation strategy only.
    Ok(())
}

/// Private detached analogue of RDKit's default Morgan atom-invariant
/// generator. It consumes the caller's prepared property-cache and ring state.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct MorganAtomInvGenerator {
    include_ring_membership: bool,
}

impl MorganAtomInvGenerator {
    pub(crate) const fn new(include_ring_membership: bool) -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganAtomInvGenerator::MorganAtomInvGenerator
        // RDKit❗🔝: MorganAtomInvGenerator::MorganAtomInvGenerator(const bool includeRingMembership)
        // RDKit❗🔝:     : df_includeRingMembership(includeRingMembership) {}
        // END RDKIT CPP FUNCTION MorganAtomInvGenerator::MorganAtomInvGenerator
        Self {
            include_ring_membership,
        }
    }

    pub(crate) fn get_atom_invariants(
        &self,
        topology: &TopologyBlock,
        valence: &ValenceAssignment,
        rings: &RingInfo,
    ) -> Result<Vec<u32>, MorganError> {
        // BEGIN RDKIT CPP FUNCTION MorganAtomInvGenerator::getAtomInvariants
        // RDKit❗🔝: std::vector<std::uint32_t> *MorganAtomInvGenerator::getAtomInvariants(
        // RDKit❗🔝:     const ROMol &mol) const {
        // RDKit❗🔝:   unsigned int nAtoms = mol.getNumAtoms();
        // RDKit❗🔝:   std::unique_ptr<std::vector<std::uint32_t>> atomInvariants(
        // RDKit❗🔝:       new std::vector<std::uint32_t>(nAtoms));
        // RDKit❗🔝:   getConnectivityInvariants(mol, *atomInvariants, df_includeRingMembership);
        // RDKit❗🔝:   return atomInvariants.release();
        // RDKit❗🔝: }
        // END RDKIT CPP FUNCTION MorganAtomInvGenerator::getAtomInvariants

        let atom_count = topology.atoms.len();
        // The shared helper's fixed stack component array preserves the source
        // sequence while avoiding one heap vector allocation per atom.
        let mut atom_invariants = vec![0; atom_count];
        get_connectivity_invariants(
            topology,
            &mut atom_invariants,
            self.include_ring_membership,
            valence,
            rings,
        )?;
        Ok(atom_invariants)
    }
}

impl Default for MorganAtomInvGenerator {
    fn default() -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganAtomInvGenerator default argument
        // RDKit❗🔝: MorganAtomInvGenerator(const bool includeRingMembership = true);
        // END RDKIT CPP FUNCTION MorganAtomInvGenerator default argument
        Self::new(true)
    }
}

impl Clone for MorganAtomInvGenerator {
    fn clone(&self) -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganAtomInvGenerator::clone
        // RDKit❗🔝: MorganAtomInvGenerator *MorganAtomInvGenerator::clone() const {
        // RDKit❗🔝:   return new MorganAtomInvGenerator(df_includeRingMembership);
        // RDKit❗🔝: }
        // END RDKIT CPP FUNCTION MorganAtomInvGenerator::clone
        // The Rust clone copies one immutable bool by value, avoiding the
        // source's heap-allocated generator while preserving independent state.
        Self::new(self.include_ring_membership)
    }
}

/// Private source-shaped Morgan bond invariants for the pinned legacy stereo
/// perception profile. The modern CIP-labeling branch remains unsupported.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct MorganBondInvGenerator {
    use_bond_types: bool,
    include_chirality: bool,
}

impl MorganBondInvGenerator {
    pub(crate) const fn new(use_bond_types: bool, include_chirality: bool) -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganBondInvGenerator::MorganBondInvGenerator
        // RDKit❗✔️: MorganBondInvGenerator::MorganBondInvGenerator(const bool useBondTypes,
        // RDKit❗✔️:                                                const bool useChirality)
        // RDKit❗✔️:     : df_useBondTypes(useBondTypes), df_useChirality(useChirality) {}
        // END RDKIT CPP FUNCTION MorganBondInvGenerator::MorganBondInvGenerator
        Self {
            use_bond_types,
            include_chirality,
        }
    }

    pub(crate) fn get_bond_invariants(&self, topology: &TopologyBlock) -> Vec<u32> {
        // BEGIN RDKIT CPP FUNCTION MorganBondInvGenerator::getBondInvariants
        // RDKit❗✔️: std::vector<std::uint32_t> *MorganBondInvGenerator::getBondInvariants(
        // RDKit❗✔️:     const ROMol &mol) const {
        // RDKit❗✔️:   std::vector<std::uint32_t> *result =
        // RDKit❗✔️:       new std::vector<std::uint32_t>(mol.getNumBonds());
        // RDKit❗✔️:   for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
        // RDKit❗✔️:     Bond const *bond = mol.getBondWithIdx(i);
        // RDKit❗✔️:     int32_t bondInvariant = 1;
        // RDKit❗✔️:     if (df_useBondTypes) {
        // RDKit❗✔️:       if (!df_useChirality || bond->getBondType() != Bond::DOUBLE ||
        // RDKit❗✔️:           bond->getStereo() == Bond::STEREONONE) {
        // RDKit❗✔️:         bondInvariant = static_cast<int32_t>(bond->getBondType());
        // RDKit❗✔️:       } else {
        // RDKit❗✔️:         auto bondStereo = static_cast<int32_t>(bond->getStereo());
        // RDKit❌❌:         if (!Chirality::getUseLegacyStereoPerception()) {
        // RDKit❌❌:           if (!mol.hasProp(common_properties::_CIPComputed)) {
        // RDKit❌❌:             CIPLabeler::assignCIPLabels(const_cast<ROMol &>(mol));
        // RDKit❌❌:           }
        // RDKit❌❌:           std::string cipCode;
        // RDKit❌❌:           if (bond->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
        // RDKit❌❌:             if (cipCode == "E") {
        // RDKit❌❌:               bondStereo = static_cast<int32_t>(Bond::STEREOE);
        // RDKit❌❌:             } else if (cipCode == "Z") {
        // RDKit❌❌:               bondStereo = static_cast<int32_t>(Bond::STEREOZ);
        // RDKit❌❌:             }
        // RDKit❌❌:           }
        // RDKit❌❌:         }
        // RDKit❗✔️:         const int32_t stereoOffset = 100;
        // RDKit❗✔️:         const int32_t bondTypeOffset = 10;
        // RDKit❗✔️:         bondInvariant =
        // RDKit❗✔️:             stereoOffset +
        // RDKit❗✔️:             bondTypeOffset * static_cast<int32_t>(bond->getBondType()) +
        // RDKit❗✔️:             bondStereo;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     (*result)[bond->getIdx()] = static_cast<int32_t>(bondInvariant);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return result;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION MorganBondInvGenerator::getBondInvariants

        // The pinned reference profile enables legacy stereo perception. Its
        // source uses the stored stereo code directly; modern CIP labeling is
        // deliberately not substituted or claimed by this detached owner.
        let mut result = vec![0; topology.bonds.len()];
        for bond in &topology.bonds {
            let mut bond_invariant = 1_i32;
            if self.use_bond_types {
                let bond_type = bond.order();
                if !self.include_chirality
                    || bond_type != BondOrder::Double
                    || bond.stereo() == BondStereo::None
                {
                    bond_invariant = bond_type.rdkit_code() as i32;
                } else {
                    let bond_stereo = bond.stereo().rdkit_code() as i32;
                    let stereo_offset = 100_i32;
                    let bond_type_offset = 10_i32;
                    bond_invariant = stereo_offset
                        + bond_type_offset * bond_type.rdkit_code() as i32
                        + bond_stereo;
                }
            }
            result[bond.id().index()] = bond_invariant as u32;
        }

        // Behavior review: each bond starts at one; false bond-type policy
        // skips all bond fields, and only typed, non-unspecified double-bond
        // stereo adds the exact 100 + 10*type + stereo packing. The validated
        // topology preserves source bond IDs, which remain the result slots.
        // Complexity review: one exact B-row allocation and one O(B) pass with
        // constant-time enum reads; no topology validation or CIP assignment
        // is repeated here. The function follows the pinned legacy profile.
        result
    }
}

impl Default for MorganBondInvGenerator {
    fn default() -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganBondInvGenerator default arguments
        // RDKit❗✔️: MorganBondInvGenerator(const bool useBondTypes = true,
        // RDKit❗✔️:                         const bool useChirality = false);
        // END RDKIT CPP FUNCTION MorganBondInvGenerator default arguments
        Self::new(true, false)
    }
}

impl Clone for MorganBondInvGenerator {
    fn clone(&self) -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganBondInvGenerator::clone
        // RDKit❗✔️: MorganBondInvGenerator *MorganBondInvGenerator::clone() const {
        // RDKit❗✔️:   return new MorganBondInvGenerator(df_useBondTypes, df_useChirality);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION MorganBondInvGenerator::clone
        Self::new(self.use_bond_types, self.include_chirality)
    }
}

/// Feature-based Morgan atom invariants. `None` retains the six built-in
/// source patterns; `Some(&[])` is an owned, present empty pattern list.
#[derive(Debug)]
pub(crate) struct MorganFeatureAtomInvGenerator {
    patterns: Option<Vec<QueryGraph>>,
}

impl MorganFeatureAtomInvGenerator {
    pub(crate) fn new(patterns: Option<&[QueryGraph]>) -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganFeatureAtomInvGenerator constructor
        // RDKit❗✔️: MorganFeatureAtomInvGenerator::MorganFeatureAtomInvGenerator(
        // RDKit❗✔️:     const std::vector<const ROMol *> *patterns) {
        // RDKit❗✔️:   if (patterns) {
        // RDKit❗✔️:     dp_patterns = new std::vector<const ROMol *>;
        // RDKit❗✔️:     dp_patterns->reserve(patterns->size());
        // RDKit❗✔️:     for (auto pattern : *patterns) {
        // RDKit❗✔️:       dp_patterns->push_back(new ROMol(*pattern));
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION MorganFeatureAtomInvGenerator constructor
        Self {
            patterns: patterns.map(<[QueryGraph]>::to_vec),
        }
    }

    pub(crate) fn get_atom_invariants(
        &self,
        target: &SearchTarget<'_>,
        query_context: &QueryMatchContext<'_>,
    ) -> Result<Vec<u32>, MorganError> {
        // BEGIN RDKIT CPP FUNCTION MorganFeatureAtomInvGenerator::getAtomInvariants
        // RDKit❗🔝: std::vector<std::uint32_t> *MorganFeatureAtomInvGenerator::getAtomInvariants(
        // RDKit❗🔝:     const ROMol &mol) const {
        // RDKit❗🔝:   unsigned int nAtoms = mol.getNumAtoms();
        // RDKit❗🔝:   std::vector<std::uint32_t> *result = new std::vector<std::uint32_t>(nAtoms);
        // RDKit❗🔝:   getFeatureInvariants(mol, *result, dp_patterns);
        // RDKit❗🔝:   return result;
        // RDKit❗🔝: }
        // END RDKIT CPP FUNCTION MorganFeatureAtomInvGenerator::getAtomInvariants

        let mut atom_invariants = vec![0; target.num_atoms()];
        // The retained query graphs are immutable borrows during matching.
        // Reusing them avoids RDKit's per-pattern ROMol copy without changing
        // the query or target state; SMARTS parsing stays outside this path.
        get_feature_invariants_with_patterns(
            target,
            &mut atom_invariants,
            self.patterns.as_deref(),
            query_context,
        )?;
        Ok(atom_invariants)
    }
}

impl Default for MorganFeatureAtomInvGenerator {
    fn default() -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganFeatureAtomInvGenerator default argument
        // RDKit❗✔️: MorganFeatureAtomInvGenerator(
        // RDKit❗✔️:     const std::vector<const ROMol *> *patterns = nullptr);
        // END RDKIT CPP FUNCTION MorganFeatureAtomInvGenerator default argument
        Self::new(None)
    }
}

impl Clone for MorganFeatureAtomInvGenerator {
    fn clone(&self) -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganFeatureAtomInvGenerator::clone
        // RDKit❗✔️: MorganFeatureAtomInvGenerator *MorganFeatureAtomInvGenerator::clone() const {
        // RDKit❗✔️:   return new MorganFeatureAtomInvGenerator(dp_patterns);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION MorganFeatureAtomInvGenerator::clone
        Self {
            patterns: self.patterns.as_ref().map(|patterns| patterns.to_vec()),
        }
    }
}

/// Parse and retain RDKit's six default feature SMARTS for the process
/// lifetime. The returned query values are immutable and borrowed by each
/// target search; custom caller-provided patterns are a separate path.
fn default_feature_queries() -> Result<&'static [QueryGraph; 6], MorganError> {
    // BEGIN RDKIT CPP FUNCTION: FingerprintUtil.cpp default feature patterns and flyweight
    // RDKit❗✔️: const char *smartsPatterns[6] = {
    // RDKit❗✔️:     "[$([N;!H0;v3,v4&+1]),\
    // RDKit❗✔️: $([O,S;H1;+0]),\
    // RDKit❗✔️: n&H1&+0]",                                                  // Donor
    // RDKit❗✔️:     "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
    // RDKit❗✔️: $([O,S;H0;v2]),\
    // RDKit❗✔️: $([O,S;-]),\
    // RDKit❗✔️: $([N;v3;!$(N-*=[O,N,P,S])]),\
    // RDKit❗✔️: n&H0&+0,\
    // RDKit❗✔️: $([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])",                    // Acceptor
    // RDKit❗✔️:     "[a]",                                                  // Aromatic
    // RDKit❗✔️:     "[F,Cl,Br,I]",                                          // Halogen
    // RDKit❗✔️:     "[#7;+,\
    // RDKit❗✔️: $([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
    // RDKit❗✔️: $([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
    // RDKit❗✔️: $([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",  // Basic
    // RDKit❗✔️:     "[$([C,S](=[O,S,P])-[O;H1,-1])]"                        // Acidic
    // RDKit❗✔️: };
    // RDKit❗✔️: std::vector<std::string> defaultFeatureSmarts(smartsPatterns,
    // RDKit❗✔️:                                               smartsPatterns + 6);
    // RDKit❗✔️: typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
    // RDKit❗✔️:                          boost::flyweights::no_tracking>
    // RDKit❗✔️:     pattern_flyweight;
    // RDKit❗✔️: ss_matcher::ss_matcher(const std::string &pattern) {
    // RDKit❗✔️:   RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    // RDKit❗✔️:   TEST_ASSERT(p);
    // RDKit❗✔️:   m_matcher.reset(p);
    // RDKit❗✔️: };
    // END RDKIT CPP FUNCTION: FingerprintUtil.cpp default feature patterns and flyweight

    const DEFAULT_FEATURE_SMARTS: [&str; 6] = [
        "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]",
        "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",
        "[a]",
        "[F,Cl,Br,I]",
        "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",
        "[$([C,S](=[O,S,P])-[O;H1,-1])]",
    ];

    static QUERIES: OnceLock<Result<[QueryGraph; 6], SmartsParseError>> = OnceLock::new();
    let cached = QUERIES.get_or_init(|| {
        let params = SmartsParseParams::default();
        Ok([
            parse_smarts(DEFAULT_FEATURE_SMARTS[0], &params)?,
            parse_smarts(DEFAULT_FEATURE_SMARTS[1], &params)?,
            parse_smarts(DEFAULT_FEATURE_SMARTS[2], &params)?,
            parse_smarts(DEFAULT_FEATURE_SMARTS[3], &params)?,
            parse_smarts(DEFAULT_FEATURE_SMARTS[4], &params)?,
            parse_smarts(DEFAULT_FEATURE_SMARTS[5], &params)?,
        ])
    });

    match cached {
        Ok(queries) => Ok(queries),
        Err(source) => Err(MorganError::SmartsParse(source.clone())),
    }
}

/// Fill source-ordered feature masks using the retained defaults or custom
/// patterns. The target and context borrow the same final prepared state.
pub(crate) fn get_feature_invariants(
    target: &SearchTarget<'_>,
    invars: &mut [u32],
    query_context: &QueryMatchContext<'_>,
) -> Result<(), MorganError> {
    get_feature_invariants_with_patterns(target, invars, None, query_context)
}

/// Apply the pinned nullable-pattern behavior: `None` selects the six defaults
/// and `Some(&[])` is a present empty list that zeroes the output without any
/// searches.
pub(crate) fn get_feature_invariants_with_patterns(
    target: &SearchTarget<'_>,
    invars: &mut [u32],
    patterns: Option<&[QueryGraph]>,
    query_context: &QueryMatchContext<'_>,
) -> Result<(), MorganError> {
    // BEGIN RDKIT CPP FUNCTION MorganFingerprints::getFeatureInvariants
    // RDKit❗❗: void getFeatureInvariants(const ROMol &mol, std::vector<uint32_t> &invars,
    // RDKit❗❗:                           const std::vector<const ROMol *> *patterns) {
    // RDKit❗❗:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗❗:   PRECONDITION(invars.size() >= nAtoms, "vector too small");
    // RDKit❌❗:   auto useLocalPatterns = patterns == nullptr;
    // RDKit❌❗:   std::vector<const ROMol *> featureMatchers;
    // RDKit❌❗:   if (useLocalPatterns) {
    // RDKit❗✔️:     featureMatchers.reserve(defaultFeatureSmarts.size());
    // RDKit❗✔️:     for (const auto &smaIt : defaultFeatureSmarts) {
    // RDKit❗✔️:       const ROMol *matcher = pattern_flyweight(smaIt).get().getMatcher();
    // RDKit❗✔️:       CHECK_INVARIANT(matcher, "bad smarts");
    // RDKit❗✔️:       featureMatchers.push_back(matcher);
    // RDKit❗❗:     }
    // RDKit❌❗:   }
    // RDKit❗❗:   std::fill(invars.begin(), invars.end(), 0);
    // RDKit❗❗:   auto &queries = (useLocalPatterns ? featureMatchers : *patterns);
    // RDKit❗❗:   for (unsigned int i = 0; i < queries.size(); ++i) {
    // RDKit❗❗:     unsigned int mask = 1 << i;
    // RDKit❗❗:     std::vector<MatchVectType> matchVect;
    // RDKit❗❗:     // to maintain thread safety, we have to copy the pattern
    // RDKit❗❗:     // molecules:
    // RDKit❗❗:     SubstructMatch(mol, ROMol(*queries[i], true), matchVect);
    // RDKit❗❗:     for (const auto &mvIt : matchVect) {
    // RDKit❗❗:       for (const auto &mIt : mvIt) {
    // RDKit❗❗:         invars[mIt.second] |= mask;
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // RDKit❗❗: }  // end of getFeatureInvariants()
    // END RDKIT CPP FUNCTION MorganFingerprints::getFeatureInvariants

    let atom_count = target.num_atoms();
    if invars.len() < atom_count {
        return Err(FingerprintError::PreconditionViolation {
            what: "vector too small",
        }
        .into());
    }

    let queries: &[QueryGraph] = match patterns {
        Some(queries) => queries,
        None => default_feature_queries()?.as_slice(),
    };
    invars.fill(0);
    let params = SubstructMatchParams::default();
    for (query_index, query) in queries.iter().enumerate() {
        // C++20 N4861 [expr.shift]/2: https://timsong-cpp.github.io/cppwp/n4861/expr.shift
        // defines `1 << i` modulo 2^N. With the pinned 32-bit `int`, ordinal31
        // is defined as 0x80000000 and assignment to unsigned preserves it.
        // [expr.shift]/1 makes only shift counts >=32 undefined, so retain the
        // existing typed source-UB boundary exactly there.
        if query_index >= i32::BITS as usize {
            return Err(FingerprintError::UndefinedArithmetic {
                site: "FingerprintUtil.cpp::getFeatureInvariants (1 << i)",
            }
            .into());
        }
        let mask = 1_u32 << query_index;
        let matches = try_get_substruct_matches_with_params_and_context(
            target,
            query,
            &params,
            query_context,
        )?;
        for matched in matches {
            for atom_index in matched.atom_mapping {
                invars[atom_index] |= mask;
            }
        }
    }

    // Behavior review: the short-buffer error occurs before pattern selection
    // or mutation. Default SMARTS are acquired before the source-ordered full
    // slice zero-fill; `Some` custom patterns preserve their supplied order,
    // including the empty-list case. Each query uses RDKit's default recursive
    // match policy, and every target index in every full atom mapping receives
    // that pattern's ordinal bit. Matching errors preserve earlier ordered
    // updates. The explicit i32 shift guard rejects only source-undefined masks.
    // Complexity review: the six default queries are parsed once into
    // OnceLock, and caller-provided patterns are borrowed. Search observes the
    // existing query and target graph rows directly and allocates its two
    // query-sized mapping buffers once per invocation for reuse at every goal.
    // One paired row is allocated only after a mapping passes final checks.
    // After enumeration, Search projects each accepted row once into atom and
    // bond maps; resolving each query bond scans the selected target neighbor
    // row. These accepted-output/map allocations and remaining source-profile
    // differences leave total cost unresolved: this path claims neither zero
    // allocation nor whole-Search cost equivalence.
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        MorganAtomInvGenerator, MorganBondInvGenerator, MorganFeatureAtomInvGenerator,
        default_feature_queries, get_connectivity_invariants, get_feature_invariants,
        get_feature_invariants_with_patterns,
    };
    use crate::generator::FingerprintArguments;
    use crate::morgan::{MorganGenerator, MorganParams, get_morgan_generator};
    use crate::{FingerprintError, MorganError};
    use cosmolkit_core::{ValenceParams, assign_valence, fast_find_rings};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, BondStereo, Element,
        TopologyBlock,
    };
    use cosmolkit_search::{
        QueryGraph, QueryMatchContextError, SearchTarget, SmartsParseError, SmartsParseParams,
        SubstructMatchError, SubstructMatchParams, build_prepared_query_match_context,
        parse_smarts,
    };
    use cosmolkit_smiles::{SmilesParseParams, SmilesRecord, parse_smiles};
    use std::error::Error;
    #[test]
    fn fingerprint_morgan_m01_defaults_and_complete_parameter_product() {
        const OPTIONS: [[bool; 6]; 64] = [
            // countSimulation, includeChirality, useBondTypes,
            // includeRingMembership, onlyNonzeroInvariants,
            // includeRedundantEnvironments
            [false, false, false, false, false, false],
            [false, false, false, false, false, true],
            [false, false, false, false, true, false],
            [false, false, false, false, true, true],
            [false, false, false, true, false, false],
            [false, false, false, true, false, true],
            [false, false, false, true, true, false],
            [false, false, false, true, true, true],
            [false, false, true, false, false, false],
            [false, false, true, false, false, true],
            [false, false, true, false, true, false],
            [false, false, true, false, true, true],
            [false, false, true, true, false, false],
            [false, false, true, true, false, true],
            [false, false, true, true, true, false],
            [false, false, true, true, true, true],
            [false, true, false, false, false, false],
            [false, true, false, false, false, true],
            [false, true, false, false, true, false],
            [false, true, false, false, true, true],
            [false, true, false, true, false, false],
            [false, true, false, true, false, true],
            [false, true, false, true, true, false],
            [false, true, false, true, true, true],
            [false, true, true, false, false, false],
            [false, true, true, false, false, true],
            [false, true, true, false, true, false],
            [false, true, true, false, true, true],
            [false, true, true, true, false, false],
            [false, true, true, true, false, true],
            [false, true, true, true, true, false],
            [false, true, true, true, true, true],
            [true, false, false, false, false, false],
            [true, false, false, false, false, true],
            [true, false, false, false, true, false],
            [true, false, false, false, true, true],
            [true, false, false, true, false, false],
            [true, false, false, true, false, true],
            [true, false, false, true, true, false],
            [true, false, false, true, true, true],
            [true, false, true, false, false, false],
            [true, false, true, false, false, true],
            [true, false, true, false, true, false],
            [true, false, true, false, true, true],
            [true, false, true, true, false, false],
            [true, false, true, true, false, true],
            [true, false, true, true, true, false],
            [true, false, true, true, true, true],
            [true, true, false, false, false, false],
            [true, true, false, false, false, true],
            [true, true, false, false, true, false],
            [true, true, false, false, true, true],
            [true, true, false, true, false, false],
            [true, true, false, true, false, true],
            [true, true, false, true, true, false],
            [true, true, false, true, true, true],
            [true, true, true, false, false, false],
            [true, true, true, false, false, true],
            [true, true, true, false, true, false],
            [true, true, true, false, true, true],
            [true, true, true, true, false, false],
            [true, true, true, true, false, true],
            [true, true, true, true, true, false],
            [true, true, true, true, true, true],
        ];

        assert_eq!(
            MorganParams::default(),
            MorganParams {
                radius: 3,
                include_chirality: false,
                use_bond_types: true,
                include_ring_membership: true,
                only_nonzero_invariants: false,
                include_redundant_environments: false,
                fp_size: 2048,
                count_simulation: false,
                count_bounds: vec![1, 2, 4, 8],
                bits_per_feature: 1,
            }
        );

        let mut cases = 0;
        for radius in [0, 1, 2, 3] {
            for [
                count_simulation,
                include_chirality,
                use_bond_types,
                include_ring_membership,
                only_nonzero_invariants,
                include_redundant_environments,
            ] in OPTIONS
            {
                let params = MorganParams {
                    radius,
                    include_chirality,
                    use_bond_types,
                    include_ring_membership,
                    only_nonzero_invariants,
                    include_redundant_environments,
                    fp_size: 2048,
                    count_simulation,
                    count_bounds: vec![1, 2, 4, 8],
                    bits_per_feature: 1,
                };
                let actual = get_morgan_generator(&params)
                    .expect("the source default count bounds are valid");

                let expected = MorganGenerator {
                    radius,
                    only_nonzero_invariants,
                    include_redundant_environments,
                    fingerprint_arguments: FingerprintArguments {
                        count_simulation,
                        include_chirality,
                        count_bounds: vec![1, 2, 4, 8],
                        fp_size: 2048,
                        bits_per_feature: 1,
                    },
                    atom_invariants: MorganAtomInvGenerator {
                        include_ring_membership,
                    },
                    bond_invariants: MorganBondInvGenerator {
                        use_bond_types,
                        include_chirality,
                    },
                };
                let options = [
                    count_simulation,
                    include_chirality,
                    use_bond_types,
                    include_ring_membership,
                    only_nonzero_invariants,
                    include_redundant_environments,
                ];
                assert_eq!(
                    actual, expected,
                    "source Morgan option product at radius {radius}, options {options:?}"
                );
                cases += 1;
            }
        }
        assert_eq!(cases, 256);
    }

    fn evaluate_features(smiles: &str) -> (SmilesRecord, Vec<u32>) {
        let record = parse_smiles(smiles, &SmilesParseParams::default())
            .expect("fixed feature target parses with source defaults");
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed feature target has source-valid valence");
        let rings = fast_find_rings(&record.topology).expect("fixed target has ring information");
        let context = build_prepared_query_match_context(&record.topology, &rings, &valence)
            .expect("prepared context borrows matching final target assignments");
        let target = SearchTarget::new(
            &record.topology,
            &record.coordinates,
            &record.topology.stereo_groups,
            Some(&rings),
            Some(&valence),
        );
        let atom_count = record.topology.atoms.len();
        let mut invariants = vec![u32::MAX; atom_count + 1];

        get_feature_invariants(&target, &mut invariants, &context)
            .expect("fixed default feature searches succeed");
        assert_eq!(invariants[atom_count], 0, "source clears the whole slice");
        invariants.truncate(atom_count);
        (record, invariants)
    }

    fn evaluate_features_with_patterns(
        smiles: &str,
        patterns: Option<&[QueryGraph]>,
    ) -> (SmilesRecord, Vec<u32>) {
        let record = parse_smiles(smiles, &SmilesParseParams::default())
            .expect("fixed feature target parses with source defaults");
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed feature target has source-valid valence");
        let rings = fast_find_rings(&record.topology).expect("fixed target has ring information");
        let context = build_prepared_query_match_context(&record.topology, &rings, &valence)
            .expect("prepared context borrows matching final target assignments");
        let target = SearchTarget::new(
            &record.topology,
            &record.coordinates,
            &record.topology.stereo_groups,
            Some(&rings),
            Some(&valence),
        );
        let atom_count = record.topology.atoms.len();
        let mut invariants = vec![u32::MAX; atom_count + 2];

        get_feature_invariants_with_patterns(&target, &mut invariants, patterns, &context)
            .expect("fixed feature-pattern searches succeed");
        assert_eq!(
            &invariants[atom_count..],
            &[0, 0],
            "source clears full output"
        );
        invariants.truncate(atom_count);
        (record, invariants)
    }

    fn evaluate_with_feature_generator(
        smiles: &str,
        generator: &MorganFeatureAtomInvGenerator,
    ) -> (SmilesRecord, Vec<u32>) {
        let record = parse_smiles(smiles, &SmilesParseParams::default())
            .expect("fixed feature-generator target parses with source defaults");
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed feature-generator target has source-valid valence");
        let rings = fast_find_rings(&record.topology)
            .expect("fixed feature-generator target has ring information");
        let context = build_prepared_query_match_context(&record.topology, &rings, &valence)
            .expect("prepared context borrows matching final target assignments");
        let target = SearchTarget::new(
            &record.topology,
            &record.coordinates,
            &record.topology.stereo_groups,
            Some(&rings),
            Some(&valence),
        );
        let invariants = generator
            .get_atom_invariants(&target, &context)
            .expect("fixed feature-generator search succeeds");
        (record, invariants)
    }

    fn parse_feature_patterns(smarts: &[&str]) -> Vec<QueryGraph> {
        let params = SmartsParseParams::default();
        smarts
            .iter()
            .map(|pattern| parse_smarts(pattern, &params).expect("fixed feature SMARTS parse"))
            .collect()
    }

    fn assert_morgan_cause<T>(error: &MorganError, expected: &T)
    where
        T: Error + std::fmt::Debug + PartialEq + 'static,
    {
        assert_eq!(
            error.source().and_then(|cause| cause.downcast_ref::<T>()),
            Some(expected)
        );
    }

    fn raw_bond_topology(order: BondOrder, stereo: BondStereo) -> TopologyBlock {
        let atoms = (0..4)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let central_bond = BondSpec::new(AtomId::new(1), AtomId::new(2), order)
            .with_stereo(stereo)
            .with_stereo_atoms(AtomId::new(0), AtomId::new(3));
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(BondId::new(1), central_bond),
            Bond::from_spec(
                BondId::new(2),
                BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
            ),
        ];

        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed raw bond rows form a valid detached topology")
    }

    #[test]
    fn fingerprint_morgan_i06_bond_type_stereo_complete_product() {
        const PROFILES: [(bool, bool); 4] =
            [(false, false), (false, true), (true, false), (true, true)];
        const ORDERS: [BondOrder; 5] = [
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Aromatic,
            BondOrder::Dative,
        ];
        const STEREOS: [BondStereo; 5] = [
            BondStereo::None,
            BondStereo::E,
            BondStereo::Z,
            BondStereo::Cis,
            BondStereo::Trans,
        ];
        const EXPECTED_CENTERS: [[[u32; 5]; 5]; 4] = [
            [[1; 5], [1; 5], [1; 5], [1; 5], [1; 5]],
            [[1; 5], [1; 5], [1; 5], [1; 5], [1; 5]],
            [[1; 5], [2; 5], [3; 5], [12; 5], [17; 5]],
            [[1; 5], [2, 123, 122, 124, 125], [3; 5], [12; 5], [17; 5]],
        ];

        let defaulted = MorganBondInvGenerator::default();
        assert!(defaulted.use_bond_types);
        assert!(!defaulted.include_chirality);
        assert_eq!(defaulted.clone(), MorganBondInvGenerator::new(true, false));

        let mut calls = 0;
        for (profile_index, flags) in PROFILES.iter().enumerate() {
            let (use_bond_types, include_chirality) = *flags;
            let generator = MorganBondInvGenerator::new(use_bond_types, include_chirality);
            for (order_index, order) in ORDERS.into_iter().enumerate() {
                for (stereo_index, stereo) in STEREOS.into_iter().enumerate() {
                    let topology = raw_bond_topology(order, stereo);
                    assert_eq!(topology.bonds.len(), 3);
                    let expected = [
                        1,
                        EXPECTED_CENTERS[profile_index][order_index][stereo_index],
                        1,
                    ];
                    assert_eq!(
                        generator.get_bond_invariants(&topology),
                        expected,
                        "useBondTypes={use_bond_types}, includeChirality={include_chirality}, \
                         order={order:?}, stereo={stereo:?}"
                    );
                    calls += 1;
                }
            }
        }
        assert_eq!(calls, 100);
    }

    #[test]
    fn fingerprint_morgan_i02_default_feature_masks_follow_source_order() {
        const CASES: [(&str, &[u32]); 6] = [
            ("[NH4+]", &[17]),
            ("CCO", &[0, 0, 3]),
            ("c1cc[nH]c1", &[4, 4, 4, 5, 4]),
            ("CCl", &[0, 8]),
            ("CN", &[0, 19]),
            ("CC(=O)O", &[0, 32, 2, 1]),
        ];

        assert!(SubstructMatchParams::default().recursion_possible);
        for (smiles, expected) in CASES {
            let (record, actual) = evaluate_features(smiles);
            assert_eq!(record.topology.atoms.len(), expected.len(), "{smiles}");
            assert_eq!(actual, expected, "{smiles}");
        }
    }

    #[test]
    fn fingerprint_morgan_i02_default_queries_are_retained() {
        let first = default_feature_queries().expect("six fixed SMARTS parse");
        let second = default_feature_queries().expect("cached SMARTS remain available");

        assert_eq!(first.len(), 6);
        assert!(std::ptr::eq(first, second));
        for (first_query, second_query) in first.iter().zip(second) {
            assert!(std::ptr::eq(first_query, second_query));
        }
    }

    #[test]
    fn fingerprint_morgan_i02_short_output_preserves_source_precondition_order() {
        let record = parse_smiles("CCO", &SmilesParseParams::default()).expect("fixed SMILES");
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("valid fixed valence");
        let rings = fast_find_rings(&record.topology).expect("fixed ring information");
        let context = build_prepared_query_match_context(&record.topology, &rings, &valence)
            .expect("prepared context");
        let target = SearchTarget::new(
            &record.topology,
            &record.coordinates,
            &record.topology.stereo_groups,
            Some(&rings),
            Some(&valence),
        );
        let mut too_short = [0xfeed_u32; 1];

        let error = get_feature_invariants(&target, &mut too_short, &context)
            .expect_err("source rejects a short output before mutation");
        assert!(matches!(
            error,
            MorganError::Fingerprint(FingerprintError::PreconditionViolation {
                what: "vector too small"
            })
        ));
        assert_eq!(too_short, [0xfeed]);
    }

    #[test]
    fn fingerprint_morgan_i02_errors_retain_typed_causes() {
        let parse_cause = SmartsParseError::UnsupportedFeature("fixed parse cause");
        assert_morgan_cause(&MorganError::from(parse_cause.clone()), &parse_cause);

        let context_cause = QueryMatchContextError::UninitializedRings;
        assert_morgan_cause(&MorganError::from(context_cause.clone()), &context_cause);

        let match_cause = SubstructMatchError::Unsupported {
            branch: "fixed matcher cause",
            rdkit_function: "SubstructMatch(fixed)",
        };
        assert_morgan_cause(&MorganError::from(match_cause.clone()), &match_cause);
    }

    #[test]
    fn fingerprint_morgan_i03_default_empty_custom_order_overlap_absence_and_mapping() {
        let (_, default_features) = evaluate_features_with_patterns("CCO", None);
        assert_eq!(default_features, [0, 0, 3]);

        let empty: [QueryGraph; 0] = [];
        let (_, empty_features) = evaluate_features_with_patterns("CCO", Some(&empty));
        assert_eq!(empty_features, [0, 0, 0]);

        let single = parse_feature_patterns(&["[O]"]);
        let (_, single_features) = evaluate_features_with_patterns("CCO", Some(&single));
        assert_eq!(single_features, [0, 0, 1]);

        let ordered = parse_feature_patterns(&["[O]", "[C]"]);
        let (_, ordered_features) = evaluate_features_with_patterns("CCO", Some(&ordered));
        assert_eq!(ordered_features, [2, 2, 1]);

        let overlapping = parse_feature_patterns(&["[C]", "[O]", "[C]"]);
        let (_, overlapping_features) = evaluate_features_with_patterns("CCO", Some(&overlapping));
        assert_eq!(overlapping_features, [5, 5, 2]);

        let absent = parse_feature_patterns(&["[N]"]);
        let (_, absent_features) = evaluate_features_with_patterns("CCO", Some(&absent));
        assert_eq!(absent_features, [0, 0, 0]);

        let mapped = parse_feature_patterns(&["CO"]);
        let (record, mapped_features) = evaluate_features_with_patterns("CCO", Some(&mapped));
        assert_eq!(record.topology.atoms.len(), 3);
        assert_eq!(mapped_features, [0, 1, 1]);
    }

    #[test]
    fn fingerprint_morgan_i03_signed_mask_boundary_is_typed() {
        let record = parse_smiles("C", &SmilesParseParams::default()).expect("fixed target");
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed target has source-valid valence");
        let rings = fast_find_rings(&record.topology).expect("fixed target has ring information");
        let context = build_prepared_query_match_context(&record.topology, &rings, &valence)
            .expect("prepared context borrows matching final target assignments");
        let target = SearchTarget::new(
            &record.topology,
            &record.coordinates,
            &record.topology.stereo_groups,
            Some(&rings),
            Some(&valence),
        );
        let params = SmartsParseParams::default();
        let mut patterns = (0..32)
            .map(|_| parse_smarts("[C]", &params).expect("fixed carbon SMARTS parses"))
            .collect::<Vec<_>>();
        let mut invariants = [u32::MAX, u32::MAX];

        get_feature_invariants_with_patterns(&target, &mut invariants, Some(&patterns), &context)
            .expect("all 32 shift counts are defined by the pinned C++20 source");
        assert_eq!(invariants, [0xffff_ffff, 0]);

        patterns.push(parse_smarts("[C]", &params).expect("fixed 33rd carbon SMARTS parses"));
        invariants = [u32::MAX, u32::MAX];
        let error = get_feature_invariants_with_patterns(
            &target,
            &mut invariants,
            Some(&patterns),
            &context,
        )
        .expect_err("ordinal 32 is the first undefined shift count");
        assert!(matches!(
            error,
            MorganError::Fingerprint(FingerprintError::UndefinedArithmetic {
                site: "FingerprintUtil.cpp::getFeatureInvariants (1 << i)"
            })
        ));
        assert_eq!(invariants, [0xffff_ffff, 0]);
    }

    #[test]
    fn fingerprint_morgan_mask_fix_private_product_preserves_all_rows_and_inputs() {
        const CASES: [(usize, [u32; 4], bool); 7] = [
            (0, [0, 0, 0, 0], false),
            (1, [1, 1, 0, 0], false),
            (30, [0x3fff_ffff, 0x2000_0000, 0x1fff_ffff, 0], false),
            (31, [0x7fff_ffff, 0x4000_0000, 0x3fff_ffff, 0], false),
            (32, [0xffff_ffff, 0x8000_0000, 0x7fff_ffff, 0], false),
            (33, [0xffff_ffff, 0, 0xffff_ffff, 0], true),
            (34, [0xffff_ffff, 0, 0xffff_ffff, 0], true),
        ];

        let parse_params = SmartsParseParams::default();
        let carbon = parse_smarts("[#6]", &parse_params).expect("fixed carbon SMARTS parses");
        let nitrogen = parse_smarts("[#7]", &parse_params).expect("fixed nitrogen SMARTS parses");
        let mut actual_calls = 0;

        for smiles in ["C", "CCO"] {
            let record = parse_smiles(smiles, &SmilesParseParams::default())
                .expect("fixed mask target parses");
            let valence = assign_valence(&record.topology, &ValenceParams::default())
                .expect("fixed mask target has source-valid valence");
            let rings = fast_find_rings(&record.topology)
                .expect("fixed mask target has paired ring information");
            let context = build_prepared_query_match_context(&record.topology, &rings, &valence)
                .expect("prepared context borrows the target assignments");
            let target = SearchTarget::new(
                &record.topology,
                &record.coordinates,
                &record.topology.stereo_groups,
                Some(&rings),
                Some(&valence),
            );
            let atom_count = record.topology.atoms.len();

            for &(pattern_count, expected_masks, should_error) in &CASES {
                for arrangement in 0..4 {
                    let patterns = (0..pattern_count)
                        .map(|ordinal| {
                            let is_carbon = match arrangement {
                                0 => true,
                                1 => ordinal + 1 == pattern_count,
                                2 => ordinal + 1 != pattern_count,
                                3 => false,
                                _ => unreachable!("four frozen pattern arrangements"),
                            };
                            if is_carbon {
                                carbon.clone()
                            } else {
                                nitrogen.clone()
                            }
                        })
                        .collect::<Vec<_>>();
                    let mut invariants = vec![u32::MAX; atom_count + 1];
                    let topology_before = record.topology.clone();
                    let coordinates_before = record.coordinates.clone();
                    let properties_before = record.properties.clone();
                    let valence_before = valence.clone();
                    let rings_before = rings.clone();

                    let result = get_feature_invariants_with_patterns(
                        &target,
                        &mut invariants,
                        Some(&patterns),
                        &context,
                    );
                    if should_error {
                        let error = result.expect_err("ordinal32 is the first undefined shift");
                        assert!(matches!(
                            error,
                            MorganError::Fingerprint(FingerprintError::UndefinedArithmetic {
                                site: "FingerprintUtil.cpp::getFeatureInvariants (1 << i)"
                            })
                        ));
                    } else {
                        result.expect("all shifts through ordinal31 are defined by C++20");
                    }

                    let carbon_mask = expected_masks[arrangement];
                    let expected = if smiles == "C" {
                        vec![carbon_mask, 0]
                    } else {
                        vec![carbon_mask, carbon_mask, 0, 0]
                    };
                    assert_eq!(
                        invariants, expected,
                        "full output row: target={smiles}, count={pattern_count}, arrangement={arrangement}"
                    );
                    assert_eq!(record.topology, topology_before);
                    assert_eq!(record.coordinates, coordinates_before);
                    assert_eq!(record.properties, properties_before);
                    assert_eq!(valence, valence_before);
                    assert_eq!(rings, rings_before);
                    actual_calls += 1;
                }
            }
        }

        assert_eq!(
            actual_calls, 56,
            "the frozen private product is 56 real calls"
        );
    }

    fn evaluate(
        smiles: &str,
        remove_hydrogens: bool,
        include_ring_membership: bool,
    ) -> (SmilesRecord, Vec<u32>) {
        let mut parse_params = SmilesParseParams::default();
        parse_params.remove_hydrogens = remove_hydrogens;
        let record = parse_smiles(smiles, &parse_params).expect("fixed SMILES parses");
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed topology has source-valid valence");
        let rings = fast_find_rings(&record.topology).expect("fixed topology has ring info");
        let mut invariants = vec![0; record.topology.atoms.len()];

        get_connectivity_invariants(
            &record.topology,
            &mut invariants,
            include_ring_membership,
            &valence,
            &rings,
        )
        .expect("prepared fixed topology produces invariants");

        (record, invariants)
    }

    fn evaluate_with_generator(
        smiles: &str,
        generator: &MorganAtomInvGenerator,
    ) -> (SmilesRecord, Vec<u32>) {
        let record = parse_smiles(smiles, &SmilesParseParams::default())
            .expect("fixed generator target parses");
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed generator target has source-valid valence");
        let rings = fast_find_rings(&record.topology).expect("fixed target has ring information");
        let invariants = generator
            .get_atom_invariants(&record.topology, &valence, &rings)
            .expect("prepared fixed target produces atom invariants");
        (record, invariants)
    }

    #[test]
    fn fingerprint_morgan_i04_ring_membership_full_vectors_and_clone() {
        const CASES: [(&str, bool, [u32; 3]); 4] = [
            ("CCC", false, [2_246_728_737, 2_245_384_272, 2_246_728_737]),
            ("CCC", true, [2_246_728_737, 2_245_384_272, 2_246_728_737]),
            ("C1CC1", false, [2_245_384_272; 3]),
            ("C1CC1", true, [2_968_968_094; 3]),
        ];

        for (smiles, include_ring_membership, expected) in CASES {
            let generator = MorganAtomInvGenerator::new(include_ring_membership);
            let cloned = generator.clone();
            for (label, candidate) in [("original", &generator), ("clone", &cloned)] {
                let (record, invariants) = evaluate_with_generator(smiles, candidate);
                assert_eq!(record.topology.atoms.len(), 3, "{smiles}, {label}");
                assert_eq!(invariants.len(), record.topology.atoms.len());
                assert_eq!(
                    invariants, expected,
                    "{smiles}, include={include_ring_membership}, {label}"
                );
            }
        }

        let defaulted = MorganAtomInvGenerator::default();
        let (record, invariants) = evaluate_with_generator("C1CC1", &defaulted);
        assert_eq!(record.topology.atoms.len(), 3);
        assert_eq!(invariants.len(), 3);
        assert_eq!(invariants, [2_968_968_094; 3]);
    }

    #[test]
    fn fingerprint_morgan_i05_default_custom_ownership_clone_and_reuse() {
        let defaulted = MorganFeatureAtomInvGenerator::default();
        assert!(defaulted.patterns.is_none());
        let default_query_ptr = default_feature_queries()
            .expect("six source default patterns are cached")
            .as_ptr();
        let mut calls = 0;
        for _ in 0..2 {
            let (record, invariants) = evaluate_with_feature_generator("CCO", &defaulted);
            calls += 1;
            assert_eq!(record.topology.atoms.len(), 3);
            assert_eq!(invariants.len(), 3);
            assert_eq!(invariants, [0, 0, 3]);
            assert!(std::ptr::eq(
                default_query_ptr,
                default_feature_queries()
                    .expect("default patterns remain cached across target calls")
                    .as_ptr()
            ));
        }

        let empty_patterns: [QueryGraph; 0] = [];
        let empty = MorganFeatureAtomInvGenerator::new(Some(&empty_patterns));
        assert!(empty.patterns.as_ref().is_some_and(Vec::is_empty));
        let (record, invariants) = evaluate_with_feature_generator("CCO", &empty);
        calls += 1;
        assert_eq!(record.topology.atoms.len(), 3);
        assert_eq!(invariants.len(), 3);
        assert_eq!(invariants, [0, 0, 0]);

        let supplied_patterns = parse_feature_patterns(&["[O]", "[C]"]);
        let supplied_ptr = supplied_patterns.as_ptr();
        let generator = MorganFeatureAtomInvGenerator::new(Some(&supplied_patterns));
        let retained_ptr = generator
            .patterns
            .as_deref()
            .expect("custom patterns remain present")
            .as_ptr();
        assert!(!std::ptr::eq(supplied_ptr, retained_ptr));
        drop(supplied_patterns);

        let cloned = generator.clone();
        let cloned_ptr = cloned
            .patterns
            .as_deref()
            .expect("clone preserves custom patterns")
            .as_ptr();
        assert!(!std::ptr::eq(retained_ptr, cloned_ptr));
        assert_eq!(generator.patterns.as_deref(), cloned.patterns.as_deref());

        for (label, candidate) in [("original", &generator), ("clone", &cloned)] {
            let candidate_ptr = candidate
                .patterns
                .as_deref()
                .expect("custom patterns remain stored for each generator")
                .as_ptr();
            for _ in 0..2 {
                let (record, invariants) = evaluate_with_feature_generator("CCO", candidate);
                calls += 1;
                assert_eq!(record.topology.atoms.len(), 3, "{label}");
                assert_eq!(invariants.len(), 3, "{label}");
                assert_eq!(invariants, [2, 2, 1], "{label}");
                assert!(std::ptr::eq(
                    candidate_ptr,
                    candidate
                        .patterns
                        .as_deref()
                        .expect("target calls retain the same custom pattern storage")
                        .as_ptr()
                ));
            }
        }
        assert_eq!(calls, 7);
    }

    #[test]
    fn fingerprint_morgan_i01_elements_and_nonring_component_omission() {
        let input = "C.N.O.CF.CCl.*";
        const ATOMIC_NUMBERS: [u8; 8] = [6, 7, 8, 6, 9, 6, 17, 0];
        const EXPECTED: [u32; 8] = [
            2_246_733_040,
            847_950_754,
            864_666_390,
            2_246_728_737,
            882_399_112,
            2_246_728_737,
            1_016_841_875,
            2_346_609_317,
        ];

        let (record, without_ring_flag) = evaluate(input, true, false);
        assert_eq!(
            record
                .topology
                .atoms
                .iter()
                .map(|atom| atom.atomic_number())
                .collect::<Vec<_>>(),
            ATOMIC_NUMBERS
        );
        assert_eq!(without_ring_flag, EXPECTED);

        let (_, with_ring_flag) = evaluate(input, true, true);
        assert_eq!(with_ring_flag, EXPECTED);
    }

    #[test]
    fn fingerprint_morgan_i01_isotope_and_charge_components() {
        const CASES: [(&str, u8, u32); 5] = [
            ("[13C]", 6, 2_244_242_024),
            ("[14C]", 6, 2_244_242_025),
            ("[100C]", 6, 2_244_242_115),
            ("[NH4+]", 7, 847_680_145),
            ("[OH-]", 8, 864_922_462),
        ];

        for (smiles, atomic_number, expected) in CASES {
            let (record, invariants) = evaluate(smiles, true, false);
            assert_eq!(record.topology.atoms.len(), 1, "{smiles}");
            assert_eq!(
                record.topology.atoms[0].atomic_number(),
                atomic_number,
                "{smiles}"
            );
            assert_eq!(invariants, [expected], "{smiles}");
        }
    }

    #[test]
    fn fingerprint_morgan_i01_attached_explicit_hydrogen_count() {
        let (record, invariants) = evaluate("C[H]", false, false);

        assert_eq!(record.topology.atoms.len(), 2);
        assert_eq!(record.topology.bonds.len(), 1);
        assert_eq!(record.topology.atoms[0].atomic_number(), 6);
        assert_eq!(record.topology.atoms[1].atomic_number(), 1);
        assert_eq!(invariants, [2_246_733_040, 4_277_593_716]);
    }

    #[test]
    fn fingerprint_morgan_i01_ring_flag_and_ring_membership_product() {
        const CASES: [(&str, bool, [u32; 3]); 4] = [
            ("CCC", false, [2_246_728_737, 2_245_384_272, 2_246_728_737]),
            ("CCC", true, [2_246_728_737, 2_245_384_272, 2_246_728_737]),
            ("C1CC1", false, [2_245_384_272; 3]),
            ("C1CC1", true, [2_968_968_094; 3]),
        ];

        for (smiles, include_ring_membership, expected) in CASES {
            let (record, invariants) = evaluate(smiles, true, include_ring_membership);
            assert_eq!(record.topology.atoms.len(), 3, "{smiles}");
            assert_eq!(
                invariants, expected,
                "{smiles}, include={include_ring_membership}"
            );
        }
    }
}
