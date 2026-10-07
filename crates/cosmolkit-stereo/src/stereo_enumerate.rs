// RDKit marker convention defined in dev/source_reproduction_protocol.md.
//
// Source reproduction protocol: dev/source_reproduction_protocol.md
//
// Ported original d892 implementation, adapted to existing detached records.
// This module reproduces RDKit Python's EnumerateStereoisomers call path with
// one typed potential-stereo owner, one flipper model, and one lazy
// configuration engine.
//
// RDKit source files:
//   EnumerateStereoisomers.cpp  (lines 1-184)
//   EnumerateStereoisomers.h    (lines 1-116)
//   Flippers.cpp                (lines 1-115)
//   Flippers.h                  (lines 1-108)

use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::sync::{Arc, Mutex};

use num_bigint::{BigInt, BigUint};

use cosmolkit_core::{
    PotentialStereoCenter, PotentialStereoError, PotentialStereoInfo, PotentialStereoSpecified,
    PotentialStereoType, ValenceModel,
};
use cosmolkit_model::StereoGroupKind;
use cosmolkit_model::{AtomId, BondId, Conformer3D, CoordinateDimension, PropertyText};
use cosmolkit_smiles::SmilesRecord;
use cosmolkit_types::{BondDirection, BondStereo, ChiralTag};

// ──────────────────────────────────────────────
// Error type
// ──────────────────────────────────────────────

#[derive(Debug, thiserror::Error)]
pub enum EnumerationError {
    #[error(transparent)]
    Valence(#[from] cosmolkit_core::ValenceError),
    #[error(transparent)]
    Rings(#[from] cosmolkit_core::RingFindingError),
    #[error(transparent)]
    Stereo(#[from] cosmolkit_core::LegacyStereoError),
    #[error(transparent)]
    DoubleBond(#[from] cosmolkit_core::DoubleBondStereoError),
    #[error(transparent)]
    SmilesWrite(#[from] cosmolkit_smiles::SmilesParseError),
    #[error(transparent)]
    PotentialStereo(#[from] PotentialStereoError),
    #[error(transparent)]
    Hydrogens(#[from] cosmolkit_core::HydrogenError),
    #[error(transparent)]
    Embedding(#[from] cosmolkit_conformer::GenerationError),
    #[error(transparent)]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error(transparent)]
    BondProperty(#[from] cosmolkit_model::BondValueError),
    #[error(transparent)]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
    #[error(transparent)]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error(transparent)]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error("stereo flipper atom index {atom} is outside the topology")]
    InvalidFlipperAtom { atom: AtomId },
    #[error("stereo flipper bond index {bond} is outside the topology")]
    InvalidFlipperBond { bond: BondId },
    #[error("random bit source failed for {bit_count} centers: {message}")]
    RandomBitsSource { bit_count: usize, message: String },
    #[error("stereoisomer random provider lock poisoned: {0}")]
    RandomSourcePoisoned(#[source] std::sync::PoisonError<()>),
    #[error("successful embedding did not return the required heavy-atom coordinate rows")]
    EmbeddingCoordinates,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FlipperSelectionOptions {
    pub(crate) only_unassigned: bool,
    pub(crate) only_stereo_groups: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum StereoFlipper {
    Atom {
        atom: AtomId,
    },
    Bond {
        bond: BondId,
    },
    StereoGroup {
        original_parities: Vec<(AtomId, ChiralTag)>,
    },
}

impl StereoFlipper {
    pub(crate) fn flip(
        &self,
        molecule: &mut SmilesRecord,
        flag: bool,
    ) -> Result<(), EnumerationError> {
        match self {
            Self::Atom { atom } => {
                // RDKit✔️✔️: class _AtomFlipper(object):
                // RDKit✔️✔️:
                // RDKit✔️✔️:   def __init__(self, atom):
                // RDKit✔️✔️:     self.atom = atom
                // RDKit✔️✔️:
                // RDKit✔️✔️:   def flip(self, flag):
                // RDKit✔️✔️:     if flag:
                // RDKit✔️✔️:       self.atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                // RDKit✔️✔️:     else:
                // RDKit✔️✔️:       self.atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                let atom_state = molecule
                    .topology
                    .atoms
                    .get_mut(atom.index())
                    .ok_or(EnumerationError::InvalidFlipperAtom { atom: *atom })?;
                atom_state.set_chiral_tag(if flag {
                    ChiralTag::TetrahedralCw
                } else {
                    ChiralTag::TetrahedralCcw
                });
            }
            Self::Bond { bond } => {
                // RDKit✔️✔️: class _BondFlipper(object):
                // RDKit✔️✔️:
                // RDKit✔️✔️:   def __init__(self, bond):
                // RDKit✔️✔️:     self.bond = bond
                // RDKit✔️✔️:
                // RDKit✔️✔️:   def flip(self, flag):
                // RDKit✔️✔️:     if flag:
                // RDKit✔️✔️:       self.bond.SetStereo(Chem.BondStereo.STEREOCIS)
                // RDKit✔️✔️:     else:
                // RDKit✔️✔️:       self.bond.SetStereo(Chem.BondStereo.STEREOTRANS)
                let bond_state = molecule
                    .topology
                    .bonds
                    .get_mut(bond.index())
                    .ok_or(EnumerationError::InvalidFlipperBond { bond: *bond })?;
                bond_state.set_stereo(if flag {
                    BondStereo::Cis
                } else {
                    BondStereo::Trans
                })?;
            }
            Self::StereoGroup { original_parities } => {
                // RDKit✔️✔️: class _StereoGroupFlipper(object):
                // RDKit✔️✔️:
                // RDKit✔️✔️:   def __init__(self, group):
                // RDKit✔️✔️:     self._original_parities = [(a, a.GetChiralTag()) for a in group.GetAtoms()]
                // RDKit✔️✔️:
                // RDKit✔️✔️:   def flip(self, flag):
                // RDKit✔️✔️:     if flag:
                // RDKit✔️✔️:       for a, original_parity in self._original_parities:
                // RDKit✔️✔️:         a.SetChiralTag(original_parity)
                // RDKit✔️✔️:     else:
                // RDKit✔️✔️:       for a, original_parity in self._original_parities:
                // RDKit✔️✔️:         if original_parity == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                // RDKit✔️✔️:           a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                // RDKit✔️✔️:         elif original_parity == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                // RDKit✔️✔️:           a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                let topology = &mut molecule.topology;
                for (atom, original_parity) in original_parities {
                    let atom_state = topology
                        .atoms
                        .get_mut(atom.index())
                        .ok_or(EnumerationError::InvalidFlipperAtom { atom: *atom })?;
                    let parity = if flag {
                        *original_parity
                    } else {
                        match original_parity {
                            ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                            ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                            other => *other,
                        }
                    };
                    atom_state.set_chiral_tag(parity);
                }
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct ConfigurationBits {
    value: BigUint,
}

impl ConfigurationBits {
    fn zero() -> Self {
        Self {
            value: BigUint::from(0_u8),
        }
    }

    fn from_biguint(value: BigUint) -> Self {
        Self { value }
    }

    pub(crate) fn bit(&self, index: usize) -> bool {
        self.value.bit(index as u64)
    }

    pub(crate) fn value(&self) -> &BigUint {
        &self.value
    }

    fn embedding_seed(&self) -> i32 {
        // RDKit✔️✔️:       # mask bitflag to fit within C++ int.
        // RDKit✔️✔️:       cid = EmbedMolecule(ntm, randomSeed=(bitflag & 0x7fffffff))
        self.value.iter_u32_digits().next().unwrap_or(0) as i32 & 0x7fff_ffff
    }
}

pub(crate) trait RandomBitsSource {
    fn getrandbits(&mut self, bit_count: usize) -> Result<ConfigurationBits, EnumerationError>;
}

struct CallbackRandomBitsSource<F> {
    callback: F,
}

impl<F> RandomBitsSource for CallbackRandomBitsSource<F>
where
    F: FnMut(usize) -> Result<BigUint, String>,
{
    fn getrandbits(&mut self, bit_count: usize) -> Result<ConfigurationBits, EnumerationError> {
        (self.callback)(bit_count)
            .map(ConfigurationBits::from_biguint)
            .map_err(|message| EnumerationError::RandomBitsSource { bit_count, message })
    }
}

const PYTHON_MT_STATE_SIZE: usize = 624;
const PYTHON_MT_MIDDLE_WORD: usize = 397;
const PYTHON_MT_MATRIX_A: u32 = 0x9908_b0df;
const PYTHON_MT_UPPER_MASK: u32 = 0x8000_0000;
const PYTHON_MT_LOWER_MASK: u32 = 0x7fff_ffff;

#[derive(Clone)]
pub(crate) struct PythonRandom {
    index: usize,
    state: [u32; PYTHON_MT_STATE_SIZE],
}

impl PythonRandom {
    pub(crate) fn from_integer_seed(seed: &BigInt) -> Self {
        // CPython✔️✔️: /* This algorithm relies on the number being unsigned.
        // CPython✔️✔️:  * So: if the arg is a PyLong, use its absolute value.
        // CPython✔️✔️:  * Otherwise use its hash value, cast to unsigned.
        // CPython✔️✔️:  */
        // CPython✔️✔️: if (PyLong_CheckExact(arg)) {
        // CPython✔️✔️:     n = PyNumber_Absolute(arg);
        // CPython✔️✔️: }
        // CPython✔️✔️: /* Now split n into 32-bit chunks, from the right. */
        // CPython✔️✔️: bits = _PyLong_NumBits(n);
        // CPython✔️✔️: /* Figure out how many 32-bit chunks this gives us. */
        // CPython✔️✔️: keyused = bits == 0 ? 1 : (bits - 1) / 32 + 1;
        let mut key = seed.magnitude().to_u32_digits();
        if key.is_empty() {
            key.push(0);
        }
        Self::from_seed_words(&key)
    }

    fn from_seed_words(key: &[u32]) -> Self {
        let mut random = Self {
            index: PYTHON_MT_STATE_SIZE,
            state: [0; PYTHON_MT_STATE_SIZE],
        };
        random.init_by_array(key);
        random
    }

    fn init_genrand(&mut self, seed: u32) {
        // CPython✔️✔️: mt[0]= s;
        // CPython✔️✔️: for (mti=1; mti<N; mti++) {
        // CPython✔️✔️:     mt[mti] =
        // CPython✔️✔️:     (1812433253U * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        // CPython✔️✔️: }
        self.state[0] = seed;
        for index in 1..PYTHON_MT_STATE_SIZE {
            self.state[index] = 1_812_433_253_u32
                .wrapping_mul(self.state[index - 1] ^ (self.state[index - 1] >> 30))
                .wrapping_add(index as u32);
        }
        self.index = PYTHON_MT_STATE_SIZE;
    }

    fn init_by_array(&mut self, key: &[u32]) {
        // CPython✔️✔️: init_genrand(self, 19650218U);
        // CPython✔️✔️: i=1; j=0;
        // CPython✔️✔️: k = (N>key_length ? N : key_length);
        // CPython✔️✔️: for (; k; k--) {
        // CPython✔️✔️:     mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525U))
        // CPython✔️✔️:              + init_key[j] + (uint32_t)j; /* non linear */
        // CPython✔️✔️:     i++; j++;
        // CPython✔️✔️:     if (i>=N) { mt[0] = mt[N-1]; i=1; }
        // CPython✔️✔️:     if (j>=key_length) j=0;
        // CPython✔️✔️: }
        // CPython✔️✔️: for (k=N-1; k; k--) {
        // CPython✔️✔️:     mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941U))
        // CPython✔️✔️:              - (uint32_t)i; /* non linear */
        // CPython✔️✔️:     i++;
        // CPython✔️✔️:     if (i>=N) { mt[0] = mt[N-1]; i=1; }
        // CPython✔️✔️: }
        // CPython✔️✔️: mt[0] = 0x80000000U; /* MSB is 1; assuring non-zero initial array */
        debug_assert!(!key.is_empty());
        self.init_genrand(19_650_218);
        let mut state_index = 1;
        let mut key_index = 0;
        for _ in 0..PYTHON_MT_STATE_SIZE.max(key.len()) {
            self.state[state_index] = (self.state[state_index]
                ^ (self.state[state_index - 1] ^ (self.state[state_index - 1] >> 30))
                    .wrapping_mul(1_664_525))
            .wrapping_add(key[key_index])
            .wrapping_add(key_index as u32);
            state_index += 1;
            key_index += 1;
            if state_index >= PYTHON_MT_STATE_SIZE {
                self.state[0] = self.state[PYTHON_MT_STATE_SIZE - 1];
                state_index = 1;
            }
            if key_index >= key.len() {
                key_index = 0;
            }
        }
        for _ in 0..PYTHON_MT_STATE_SIZE - 1 {
            self.state[state_index] = (self.state[state_index]
                ^ (self.state[state_index - 1] ^ (self.state[state_index - 1] >> 30))
                    .wrapping_mul(1_566_083_941))
            .wrapping_sub(state_index as u32);
            state_index += 1;
            if state_index >= PYTHON_MT_STATE_SIZE {
                self.state[0] = self.state[PYTHON_MT_STATE_SIZE - 1];
                state_index = 1;
            }
        }
        self.state[0] = 0x8000_0000;
    }

    fn genrand_uint32(&mut self) -> u32 {
        // CPython✔️✔️: if (self->index >= N) { /* generate N words at one time */
        // CPython✔️✔️:     int kk;
        // CPython✔️✔️:     for (kk=0;kk<N-M;kk++) {
        // CPython✔️✔️:         y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        // CPython✔️✔️:         mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1U];
        // CPython✔️✔️:     }
        // CPython✔️✔️:     for (;kk<N-1;kk++) {
        // CPython✔️✔️:         y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        // CPython✔️✔️:         mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1U];
        // CPython✔️✔️:     }
        // CPython✔️✔️:     y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        // CPython✔️✔️:     mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1U];
        // CPython✔️✔️:     self->index = 0;
        // CPython✔️✔️: }
        if self.index >= PYTHON_MT_STATE_SIZE {
            for index in 0..PYTHON_MT_STATE_SIZE - PYTHON_MT_MIDDLE_WORD {
                let value = (self.state[index] & PYTHON_MT_UPPER_MASK)
                    | (self.state[index + 1] & PYTHON_MT_LOWER_MASK);
                self.state[index] = self.state[index + PYTHON_MT_MIDDLE_WORD]
                    ^ (value >> 1)
                    ^ if value & 1 == 0 {
                        0
                    } else {
                        PYTHON_MT_MATRIX_A
                    };
            }
            for index in PYTHON_MT_STATE_SIZE - PYTHON_MT_MIDDLE_WORD..PYTHON_MT_STATE_SIZE - 1 {
                let value = (self.state[index] & PYTHON_MT_UPPER_MASK)
                    | (self.state[index + 1] & PYTHON_MT_LOWER_MASK);
                self.state[index] = self.state
                    [index + PYTHON_MT_MIDDLE_WORD - PYTHON_MT_STATE_SIZE]
                    ^ (value >> 1)
                    ^ if value & 1 == 0 {
                        0
                    } else {
                        PYTHON_MT_MATRIX_A
                    };
            }
            let value = (self.state[PYTHON_MT_STATE_SIZE - 1] & PYTHON_MT_UPPER_MASK)
                | (self.state[0] & PYTHON_MT_LOWER_MASK);
            self.state[PYTHON_MT_STATE_SIZE - 1] = self.state[PYTHON_MT_MIDDLE_WORD - 1]
                ^ (value >> 1)
                ^ if value & 1 == 0 {
                    0
                } else {
                    PYTHON_MT_MATRIX_A
                };
            self.index = 0;
        }

        // CPython✔️✔️: y = mt[self->index++];
        // CPython✔️✔️: y ^= (y >> 11);
        // CPython✔️✔️: y ^= (y << 7) & 0x9d2c5680U;
        // CPython✔️✔️: y ^= (y << 15) & 0xefc60000U;
        // CPython✔️✔️: y ^= (y >> 18);
        // CPython✔️✔️: return y;
        let mut value = self.state[self.index];
        self.index += 1;
        value ^= value >> 11;
        value ^= (value << 7) & 0x9d2c_5680;
        value ^= (value << 15) & 0xefc6_0000;
        value ^= value >> 18;
        value
    }
}

impl RandomBitsSource for PythonRandom {
    fn getrandbits(&mut self, bit_count: usize) -> Result<ConfigurationBits, EnumerationError> {
        // CPython✔️✔️: if (k < 0) {
        // CPython✔️✔️:     PyErr_SetString(PyExc_ValueError,
        // CPython✔️✔️:                     "number of bits must be non-negative");
        // CPython✔️✔️:     return NULL;
        // CPython✔️✔️: }
        // The Rust boundary is `usize`, so the source's negative branch is
        // structurally unrepresentable.
        // CPython✔️✔️: if (k == 0)
        // CPython✔️✔️:     return PyLong_FromLong(0);
        if bit_count == 0 {
            return Ok(ConfigurationBits::zero());
        }
        // CPython✔️✔️: if (k <= 32)  /* Fast path */
        // CPython✔️✔️:     return PyLong_FromUnsignedLong(genrand_uint32(self) >> (32 - k));
        if bit_count <= 32 {
            return Ok(ConfigurationBits::from_biguint(BigUint::from(
                self.genrand_uint32() >> (32 - bit_count),
            )));
        }

        // CPython✔️✔️: words = (Py_ssize_t)((k - 1u) / 32u + 1u);
        // CPython✔️✔️: /* Fill-out bits of long integer, by 32-bit words, from least significant
        // CPython✔️✔️:    to most significant. */
        // CPython✔️✔️: for (i = 0; i < words; i++, k -= 32)
        // CPython✔️✔️: {
        // CPython✔️✔️:     r = genrand_uint32(self);
        // CPython✔️✔️:     if (k < 32)
        // CPython✔️✔️:         r >>= (32 - k);  /* Drop least significant bits */
        // CPython✔️✔️:     wordarray[i] = r;
        // CPython✔️✔️: }
        let word_count = (bit_count - 1) / 32 + 1;
        let mut remaining = bit_count;
        let mut value = BigUint::from(0_u8);
        for word_index in 0..word_count {
            let mut word = self.genrand_uint32();
            if remaining < 32 {
                word >>= 32 - remaining;
            }
            value |= BigUint::from(word) << (word_index * 32);
            remaining = remaining.saturating_sub(32);
        }
        Ok(ConfigurationBits::from_biguint(value))
    }
}

fn cpython_tuple_hash(values: &[i64]) -> i64 {
    const XXPRIME_1: u64 = 11_400_714_785_074_694_791;
    const XXPRIME_2: u64 = 14_029_467_366_897_019_727;
    const XXPRIME_5: u64 = 2_870_177_450_012_600_261;

    // CPython✔️✔️: Py_uhash_t acc = _PyHASH_XXPRIME_5;
    // CPython✔️✔️: for (i = 0; i < len; i++) {
    // CPython✔️✔️:     Py_uhash_t lane = PyObject_Hash(item[i]);
    // CPython✔️✔️:     if (lane == (Py_uhash_t)-1) {
    // CPython✔️✔️:         return -1;
    // CPython✔️✔️:     }
    // CPython✔️✔️:     acc += lane * _PyHASH_XXPRIME_2;
    // CPython✔️✔️:     acc = _PyHASH_XXROTATE(acc);
    // CPython✔️✔️:     acc *= _PyHASH_XXPRIME_1;
    // CPython✔️✔️: }
    let mut accumulator = XXPRIME_5;
    for &value in values {
        let lane = value as u64;
        accumulator = accumulator.wrapping_add(lane.wrapping_mul(XXPRIME_2));
        accumulator = accumulator.rotate_left(31);
        accumulator = accumulator.wrapping_mul(XXPRIME_1);
    }
    // CPython✔️✔️: /* Add input length, mangled to keep the historical value of hash(()). */
    // CPython✔️✔️: acc += len ^ (_PyHASH_XXPRIME_5 ^ 3527539UL);
    accumulator = accumulator.wrapping_add((values.len() as u64) ^ (XXPRIME_5 ^ 3_527_539));
    // CPython✔️✔️: if (acc == (Py_uhash_t)-1) {
    // CPython✔️✔️:     return 1546275796;
    // CPython✔️✔️: }
    // CPython✔️✔️: return acc;
    if accumulator == u64::MAX {
        1_546_275_796
    } else {
        accumulator as i64
    }
}

fn default_python_random_seed(molecule: &SmilesRecord) -> BigInt {
    // RDKit✔️✔️: if options.rand is None:
    // RDKit✔️✔️:   # deterministic random seed invariant to input atom order
    // RDKit✔️✔️:   seed = hash(tuple(sorted([(a.GetDegree(), a.GetAtomicNum()) for a in tm.GetAtoms()])))
    // RDKit✔️✔️:   rand = random.Random(seed)
    let adjacency = &molecule.topology.adjacency;
    let mut atom_invariants = molecule
        .topology
        .atoms
        .iter()
        .map(|atom| {
            (
                adjacency.neighbors_of(atom.id().index()).len(),
                atom.atomic_number(),
            )
        })
        .collect::<Vec<_>>();
    atom_invariants.sort_unstable();
    let inner_hashes = atom_invariants
        .into_iter()
        .map(|(degree, atomic_number)| {
            cpython_tuple_hash(&[degree as i64, i64::from(atomic_number)])
        })
        .collect::<Vec<_>>();
    BigInt::from(cpython_tuple_hash(&inner_hashes))
}

fn theoretical_configuration_count(center_count: usize) -> BigUint {
    BigUint::from(1_u8) << center_count
}

pub(crate) struct RangeBitsGenerator {
    next: BigUint,
    end: BigUint,
}

impl RangeBitsGenerator {
    pub(crate) fn new(center_count: usize) -> Self {
        // RDKit✔️✔️: class _RangeBitsGenerator(object):
        // RDKit✔️✔️:
        // RDKit✔️✔️:   def __init__(self, nCenters):
        // RDKit✔️✔️:     self.nCenters = nCenters
        Self {
            next: BigUint::from(0_u8),
            end: theoretical_configuration_count(center_count),
        }
    }
}

impl Iterator for RangeBitsGenerator {
    type Item = ConfigurationBits;

    fn next(&mut self) -> Option<Self::Item> {
        // RDKit✔️✔️:   def __iter__(self):
        // RDKit✔️✔️:     for val in range(2**self.nCenters):
        // RDKit✔️✔️:       yield val
        if self.next >= self.end {
            return None;
        }
        let value = self.next.clone();
        self.next += 1_u8;
        Some(ConfigurationBits::from_biguint(value))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = &self.end - &self.next;
        match usize::try_from(&remaining) {
            Ok(remaining) => (remaining, Some(remaining)),
            Err(_) => (usize::MAX, None),
        }
    }
}

pub(crate) struct UniqueRandomBitsGenerator<R> {
    center_count: usize,
    _max_isomers: usize,
    random: R,
    already_seen: HashSet<ConfigurationBits>,
    theoretical_count: BigUint,
}

impl<R> UniqueRandomBitsGenerator<R> {
    pub(crate) fn new(center_count: usize, max_isomers: usize, random: R) -> Self {
        // RDKit✔️✔️: class _UniqueRandomBitsGenerator(object):
        // RDKit✔️✔️:
        // RDKit✔️✔️:   def __init__(self, nCenters, maxIsomers, rand):
        // RDKit✔️✔️:     self.nCenters = nCenters
        // RDKit✔️✔️:     self.maxIsomers = maxIsomers
        // RDKit✔️✔️:     self.rand = rand
        // RDKit✔️✔️:     self.already_seen = set()
        Self {
            center_count,
            _max_isomers: max_isomers,
            random,
            already_seen: HashSet::new(),
            theoretical_count: theoretical_configuration_count(center_count),
        }
    }
}

impl<R: RandomBitsSource> Iterator for UniqueRandomBitsGenerator<R> {
    type Item = Result<ConfigurationBits, EnumerationError>;

    fn next(&mut self) -> Option<Self::Item> {
        // RDKit✔️✔️:   def __iter__(self):
        // RDKit✔️✔️:     # note: important that this is not 'while True' otherwise it
        // RDKit✔️✔️:     # would be possible to have an infinite loop caused by all
        // RDKit✔️✔️:     # isomers failing the embedding process
        // RDKit✔️✔️:     while len(self.already_seen) < 2**self.nCenters:
        while BigUint::from(self.already_seen.len()) < self.theoretical_count {
            // RDKit✔️✔️:       bits = self.rand.getrandbits(self.nCenters)
            let bits = match self.random.getrandbits(self.center_count) {
                Ok(bits) => bits,
                Err(error) => return Some(Err(error)),
            };
            // RDKit✔️✔️:       if bits in self.already_seen:
            // RDKit✔️✔️:         continue
            if !self.already_seen.insert(bits.clone()) {
                continue;
            }

            // RDKit✔️✔️:       self.already_seen.add(bits)
            // RDKit✔️✔️:       yield bits
            return Some(Ok(bits));
        }
        None
    }
}

pub(crate) fn select_stereo_flippers(
    molecule: &mut SmilesRecord,
    options: FlipperSelectionOptions,
) -> Result<Vec<StereoFlipper>, EnumerationError> {
    // RDKit❗✔️: def _getFlippers(mol, options):
    // RDKit❗✔️:   sinfo = Chem.FindPotentialStereo(mol)
    // RDKit❗✔️:   flippers = []
    // RDKit❗✔️:   if not options.onlyStereoGroups:
    // RDKit❗✔️:     for si in sinfo:
    // RDKit❗✔️:       if options.onlyUnassigned and si.specified not in (Chem.StereoSpecified.Unspecified,
    // RDKit❗✔️:                                                          Chem.StereoSpecified.Unknown):
    // RDKit❗✔️:         continue
    // RDKit❗✔️:       if si.type == Chem.StereoType.Atom_Tetrahedral:
    // RDKit❗✔️:         flippers.append(_AtomFlipper(mol.GetAtomWithIdx(si.centeredOn)))
    // RDKit❗✔️:       elif si.type == Chem.StereoType.Bond_Double:
    // RDKit❗✔️:         bnd = mol.GetBondWithIdx(si.centeredOn)
    // RDKit❗✔️:         if not bnd.GetStereoAtoms():
    // RDKit❗✔️:           if si.controllingAtoms[0] == Chem.Atom.NOATOM or \
    // RDKit❗✔️:             si.controllingAtoms[2] == Chem.Atom.NOATOM:
    // RDKit❗✔️:             continue
    // RDKit❗✔️:           bnd.SetStereoAtoms(si.controllingAtoms[0], si.controllingAtoms[2])
    // RDKit❗✔️:         flippers.append(_BondFlipper(mol.GetBondWithIdx(si.centeredOn)))
    // RDKit❗✔️:       ## FIX: support atropisomers
    // RDKit❗✔️:
    // RDKit❗✔️:   if options.onlyUnassigned:
    // RDKit❗✔️:     # otherwise these will be counted twice
    // RDKit❗✔️:     for group in mol.GetStereoGroups():
    // RDKit❗✔️:       if group.GetGroupType() != Chem.StereoGroupType.STEREO_ABSOLUTE:
    // RDKit❗✔️:         flippers.append(_StereoGroupFlipper(group))
    // RDKit❗✔️:
    // RDKit❗✔️:   return flippers
    let valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &molecule.topology,
        ValenceModel::RdkitLike,
        false,
    )?;
    let rings = cosmolkit_core::symmetrized_sssr(&molecule.topology, &Default::default())?;
    let assignment = cosmolkit_core::potential_stereo(
        &molecule.topology,
        &valence,
        &rings,
        &Default::default(),
    )?;
    if let Some(topology) = assignment.cleaned_topology {
        molecule.topology = topology;
    }
    let stereo_info = assignment.stereo;
    select_stereo_flippers_from_info(molecule, &stereo_info, options)
}

fn select_stereo_flippers_from_info(
    molecule: &mut SmilesRecord,
    stereo_info: &[PotentialStereoInfo],
    options: FlipperSelectionOptions,
) -> Result<Vec<StereoFlipper>, EnumerationError> {
    let mut flippers = Vec::new();

    if !options.only_stereo_groups {
        for info in stereo_info {
            if options.only_unassigned
                && !matches!(
                    info.specified,
                    PotentialStereoSpecified::Unspecified | PotentialStereoSpecified::Unknown
                )
            {
                continue;
            }
            match (info.stereo_type, info.centered_on) {
                (PotentialStereoType::AtomTetrahedral, PotentialStereoCenter::Atom(atom)) => {
                    flippers.push(StereoFlipper::Atom { atom });
                }
                (PotentialStereoType::BondDouble, PotentialStereoCenter::Bond(bond)) => {
                    let bond_state = molecule
                        .topology
                        .bonds
                        .get(bond.index())
                        .ok_or(EnumerationError::InvalidFlipperBond { bond })?;
                    if bond_state.stereo_atoms().is_none() {
                        let Some(Some(begin_controller)) = info.controlling_atoms.first().copied()
                        else {
                            continue;
                        };
                        let Some(Some(end_controller)) = info.controlling_atoms.get(2).copied()
                        else {
                            continue;
                        };
                        let bond_state = molecule
                            .topology
                            .bonds
                            .get_mut(bond.index())
                            .ok_or(EnumerationError::InvalidFlipperBond { bond })?;
                        bond_state.set_stereo_atoms(Some([begin_controller, end_controller]));
                    }
                    flippers.push(StereoFlipper::Bond { bond });
                }
                // The Python source deliberately has no atropisomer flipper branch.
                (PotentialStereoType::BondAtropisomer, PotentialStereoCenter::Bond(_)) => {}
                _ => {}
            }
        }
    }

    if options.only_unassigned {
        for group in &molecule.topology.stereo_groups {
            if group.kind() == StereoGroupKind::Absolute {
                continue;
            }
            let mut original_parities = Vec::with_capacity(group.atoms().len());
            for atom in group.atoms() {
                let atom_state = molecule
                    .topology
                    .atoms
                    .get(atom.index())
                    .ok_or(EnumerationError::InvalidFlipperAtom { atom: *atom })?;
                original_parities.push((*atom, atom_state.chiral_tag()));
            }
            flippers.push(StereoFlipper::StereoGroup { original_parities });
        }
    }

    Ok(flippers)
}

#[derive(Debug)]
pub(crate) struct EnumerationWorkspace {
    molecule: SmilesRecord,
    flippers: Vec<StereoFlipper>,
    project_query: cosmolkit_conformer::PruningQueryProjection,
}

impl EnumerationWorkspace {
    pub(crate) fn prepare(
        source: &SmilesRecord,
        options: FlipperSelectionOptions,
        project_query: cosmolkit_conformer::PruningQueryProjection,
    ) -> Result<Self, EnumerationError> {
        // RDKit❗✔️:   tm = Chem.Mol(m)
        let mut molecule = source.clone();

        // RDKit❗✔️:   for atom in tm.GetAtoms():
        // RDKit❗✔️:     atom.ClearProp("_CIPCode")
        // RDKit❗✔️:   for bond in tm.GetBonds():
        // RDKit❗✔️:     if bond.GetBondDir() == Chem.BondDir.EITHERDOUBLE or bond.GetBondDir() == Chem.BondDir.UNKNOWN:
        // RDKit❗✔️:       bond.SetBondDir(Chem.BondDir.NONE)
        let topology = &mut molecule.topology;
        for atom in &mut topology.atoms {
            atom.clear_prop("_CIPCode")?;
        }
        for bond in &mut topology.bonds {
            if matches!(
                bond.direction(),
                BondDirection::EitherDouble | BondDirection::Unknown
            ) {
                bond.set_direction(BondDirection::None);
            }
        }

        // RDKit❗✔️:   flippers = _getFlippers(tm, options)
        // RDKit❗✔️:   nCenters = len(flippers)
        let flippers = select_stereo_flippers(&mut molecule, options)?;

        // RDKit❗✔️:   if not nCenters:
        // RDKit❗✔️:     yield tm
        // RDKit❗✔️:     return
        // The no-center yield is owned by the later lazy iterator. Preserving
        // an empty flipper list here keeps that branch explicit without
        // applying the chiral flag reserved for enumerated workspaces.
        if !flippers.is_empty() {
            // RDKit❗✔️:   tm.SetProp('_MolFileChiralFlag', '1')
            molecule.properties.set_prop("_MolFileChiralFlag", "1")?;
        }

        Ok(Self {
            molecule,
            flippers,
            project_query,
        })
    }

    pub(crate) fn center_count(&self) -> usize {
        self.flippers.len()
    }

    pub(crate) fn molecule(&self) -> &SmilesRecord {
        &self.molecule
    }

    pub(crate) fn into_molecule(self) -> SmilesRecord {
        self.molecule
    }

    pub(crate) fn apply_configuration(
        &mut self,
        configuration: &ConfigurationBits,
    ) -> Result<(), EnumerationError> {
        // RDKit❗✔️:   for bitflag in bitsource:
        // RDKit❗✔️:     for i in range(nCenters):
        // RDKit❗✔️:       flag = bool(bitflag & (1 << i))
        // RDKit❗✔️:       flippers[i].flip(flag)
        for (index, flipper) in self.flippers.iter().enumerate() {
            flipper.flip(&mut self.molecule, configuration.bit(index))?;
        }
        Ok(())
    }

    pub(crate) fn finalize_configuration(
        &self,
        unique: bool,
        isomers_seen: &mut HashSet<PropertyText>,
    ) -> Result<Option<SmilesRecord>, EnumerationError> {
        // RDKit❗✔️:     # from this point on we no longer need the stereogroups (if any are there), so
        // RDKit❗✔️:     # remove them:
        // RDKit❗✔️:     if tm.GetStereoGroups():
        // RDKit❗✔️:       isomer = Chem.RWMol(tm)
        // RDKit❗✔️:       isomer.SetStereoGroups([])
        // RDKit❗✔️:     else:
        // RDKit❗✔️:       isomer = Chem.Mol(tm)
        let mut isomer = self.molecule.clone();
        if !isomer.topology.stereo_groups.is_empty() {
            isomer.topology.stereo_groups.clear();
        }

        // RDKit❗✔️:     Chem.SetDoubleBondNeighborDirections(isomer)
        let rings = cosmolkit_core::symmetrized_sssr(&isomer.topology, &Default::default())?;
        isomer.topology =
            cosmolkit_core::set_double_bond_neighbor_directions(isomer.topology, &rings, None)?
                .topology;

        // RDKit❗✔️:     isomer.ClearComputedProps(includeRings=False)
        clear_computed_props_preserving_rings(&mut isomer)?;

        // RDKit❗✔️:     Chem.AssignStereochemistry(isomer, cleanIt=True, force=True, flagPossibleStereoCenters=True)
        // `ClearComputedProps()` removed `_StereochemDone`, so the source's
        // `force=True` branch cannot return early. The pinned RDKit reference
        // uses legacy stereo perception, which is the shared implementation
        // called here with the remaining two source arguments set to true.
        let valence = cosmolkit_core::assign_valence_with_options_for_topology(
            &isomer.topology,
            ValenceModel::RdkitLike,
            false,
        )?;
        isomer.topology = cosmolkit_core::assign_legacy_stereochemistry_with_flags(
            isomer.topology,
            &valence,
            &rings,
            true,
            true,
        )?;
        isomer.properties.set_computed_prop(
            "_StereochemDone",
            cosmolkit_model::PropertyValue::Bool(true),
        )?;

        // RDKit❗✔️:     if options.unique:
        // RDKit❗✔️:       cansmi = Chem.MolToSmiles(isomer, isomericSmiles=True)
        // RDKit❗✔️:       if cansmi in isomersSeen:
        // RDKit❗✔️:         continue
        // RDKit❗✔️:
        // RDKit❗✔️:       isomersSeen.add(cansmi)
        if unique {
            let canonical_isomeric_smiles = cosmolkit_smiles::write_smiles(&isomer)?;
            if !isomers_seen.insert(canonical_isomeric_smiles) {
                return Ok(None);
            }
        }

        Ok(Some(isomer))
    }

    pub(crate) fn apply_embedding_filter(
        &self,
        mut isomer: SmilesRecord,
        configuration: &ConfigurationBits,
    ) -> Result<Option<SmilesRecord>, EnumerationError> {
        // RDKit❗✔️:     if options.tryEmbedding:
        // RDKit❗✔️:       ntm = Chem.AddHs(isomer)
        let with_hydrogens = cosmolkit_core::add_hydrogens_with_source_valence(
            isomer.topology.clone(),
            isomer.coordinates.clone(),
            isomer.properties.clone(),
            &Default::default(),
            Some(cosmolkit_core::assign_valence_with_options_for_topology(
                &isomer.topology,
                ValenceModel::RdkitLike,
                false,
            )?),
        )?;

        // RDKit❗✔️:       # mask bitflag to fit within C++ int.
        // RDKit❗✔️:       cid = EmbedMolecule(ntm, randomSeed=(bitflag & 0x7fffffff))
        // The Python call selects rdDistGeom's legacy keyword overload, whose
        // defaults differ from both bare `EmbedParameters` and ETKDGv3.
        // RDKit❗✔️:   python::def(
        // RDKit❗✔️:       "EmbedMolecule", RDKit::EmbedMolecule,
        // RDKit❗✔️:       (python::arg("mol"), python::arg("maxAttempts") = 0,
        // RDKit❗✔️:        python::arg("randomSeed") = -1, python::arg("clearConfs") = true,
        // RDKit❗✔️:        python::arg("useRandomCoords") = false, python::arg("boxSizeMult") = 2.0,
        // RDKit❗✔️:        python::arg("randNegEig") = true, python::arg("numZeroFail") = 1,
        // RDKit❗✔️:        python::arg("coordMap") = python::dict(), python::arg("forceTol") = 1e-3,
        // RDKit❗✔️:        python::arg("ignoreSmoothingFailures") = false,
        // RDKit❗✔️:        python::arg("enforceChirality") = true,
        // RDKit❗✔️:        python::arg("useExpTorsionAnglePrefs") = true,
        // RDKit❗✔️:        python::arg("useBasicKnowledge") = true,
        // RDKit❗✔️:        python::arg("printExpTorsionAngles") = false,
        // RDKit❗✔️:        python::arg("useSmallRingTorsions") = false,
        // RDKit❗✔️:        python::arg("useMacrocycleTorsions") = true,
        // RDKit❗✔️:        python::arg("ETversion") = 2,
        // RDKit❗✔️:        python::arg("useMacrocycle14config") = true),
        // RDKit❗✔️:       docString.c_str());
        let mut params = cosmolkit_conformer::EmbedParams::default();
        params.random_seed = configuration.embedding_seed();
        params.use_exp_torsion_angle_prefs = true;
        params.use_basic_knowledge = true;
        params.use_macrocycle_torsions = true;
        params.use_macrocycle14config = true;
        let generated = cosmolkit_conformer::generate_conformers(
            &with_hydrogens.topology,
            &with_hydrogens.coordinates,
            &with_hydrogens.properties,
            1,
            &mut params,
            self.project_query,
        )?;
        let conformer_id = generated.conf_ids.first().copied().unwrap_or(-1);

        // RDKit❗✔️:       if cid >= 0:
        if conformer_id < 0 {
            // RDKit❗✔️:     if cid >= 0:
            // The caller interprets `None` as the source's failed-embedding
            // branch and advances the finite configuration source.
            return Ok(None);
        }

        // RDKit❗✔️:         conf = Chem.Conformer(isomer.GetNumAtoms())
        // RDKit❗✔️:         for aid in range(isomer.GetNumAtoms()):
        // RDKit❗✔️:           conf.SetAtomPosition(aid, ntm.GetConformer().GetAtomPosition(aid))
        let heavy_atom_count = isomer.topology.atoms.len();
        let embedded_coordinates = generated
            .conformers
            .first()
            .ok_or(EnumerationError::EmbeddingCoordinates)?
            .coordinates();
        if embedded_coordinates.len() < heavy_atom_count {
            return Err(EnumerationError::EmbeddingCoordinates);
        }
        let coordinates = embedded_coordinates[..heavy_atom_count].to_vec();

        // RDKit❗✔️:         isomer.AddConformer(conf)
        // Python's `AddConformer` wrapper defaults `assignId` to false, and a
        // newly constructed RDKit conformer has id zero. Preserve that exact
        // append behavior, including a duplicate id when conformers already
        // exist, inside this private owned workspace.
        let coordinate_block = &mut isomer.coordinates;
        coordinate_block
            .conformers_3d
            .push(Conformer3D::new(0, coordinates, true));
        coordinate_block.record_source_conformer_append(CoordinateDimension::ThreeD)?;

        // RDKit❗✔️:     else:
        // RDKit❗✔️:       cid = 1
        // The non-embedding branch is represented by the caller retaining the
        // finalized isomer directly; this function is called only when the
        // option is true.
        // RDKit❗✔️:     if cid >= 0:
        Ok(Some(isomer))
    }

    pub(crate) fn finalize_configuration_with_embedding(
        &self,
        configuration: &ConfigurationBits,
        unique: bool,
        try_embedding: bool,
        isomers_seen: &mut HashSet<PropertyText>,
    ) -> Result<Option<SmilesRecord>, EnumerationError> {
        let Some(isomer) = self.finalize_configuration(unique, isomers_seen)? else {
            // RDKit❗✔️:       if cansmi in isomersSeen:
            // RDKit❗✔️:         continue
            return Ok(None);
        };

        // RDKit❗✔️:     if options.tryEmbedding:
        if try_embedding {
            self.apply_embedding_filter(isomer, configuration)
        } else {
            // RDKit❗✔️:     else:
            // RDKit❗✔️:       cid = 1
            // RDKit❗✔️:     if cid >= 0:
            Ok(Some(isomer))
        }
    }
}

fn clear_computed_props_preserving_rings(
    molecule: &mut SmilesRecord,
) -> Result<(), EnumerationError> {
    // RDKit❗✔️: void ROMol::clearComputedProps(bool includeRings) const {
    // RDKit❗✔️:   // the SSSR information:
    // RDKit❗✔️:   if (includeRings) {
    // RDKit❗✔️:     this->dp_ringInfo->reset();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   RDProps::clearComputedProps();
    // RDKit❗✔️:   for (auto atom : atoms()) {
    // RDKit❗✔️:     atom->clearComputedProps();
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   for (auto bond : bonds()) {
    // RDKit❗✔️:     bond->clearComputedProps();
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Detached owners hold no runtime cache. The caller retains the ring value
    // across this clear; the runtime owns preservation and invalidation effects.
    molecule.properties.clear_computed_props()?;
    for atom in &mut molecule.topology.atoms {
        atom.clear_computed_props()?;
    }
    for bond in &mut molecule.topology.bonds {
        bond.clear_computed_props()?;
    }
    Ok(())
}

// ──────────────────────────────────────────────
// Python-parity lazy public API
// ──────────────────────────────────────────────

/// Source-defined seed or shared random-bits provider.
///
/// Cloning a seed does not share an advanced generator. Cloning a provider
/// retains its identity and state without locking or drawing any bits.
#[derive(Clone)]
pub struct StereoisomerRandomSource(RandomSourceKind);

#[derive(Clone)]
enum RandomSourceKind {
    IntegerSeed(BigInt),
    RandomBits(
        Arc<Mutex<Box<dyn FnMut(usize) -> Result<BigUint, String> + Send + Sync + 'static>>>,
    ),
}

impl StereoisomerRandomSource {
    #[must_use]
    pub fn from_integer_seed(seed: BigInt) -> Self {
        // RDKit❗✔️:       rand = random.Random(options.rand)
        // CPython integer seeding occurs in each enumeration, not here.
        Self(RandomSourceKind::IntegerSeed(seed))
    }

    #[must_use]
    pub fn from_random_bits(
        callback: Box<dyn FnMut(usize) -> Result<BigUint, String> + Send + Sync + 'static>,
    ) -> Self {
        // RDKit❗✔️:     self.rand = rand
        // This parameter preserves a real provider's identity and mutable state.
        Self(RandomSourceKind::RandomBits(Arc::new(Mutex::new(callback))))
    }
}

impl std::fmt::Debug for StereoisomerRandomSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.0 {
            RandomSourceKind::IntegerSeed(seed) => {
                f.debug_tuple("IntegerSeed").field(seed).finish()
            }
            RandomSourceKind::RandomBits(_) => f.write_str("RandomBits(..)"),
        }
    }
}

impl PartialEq for StereoisomerRandomSource {
    fn eq(&self, other: &Self) -> bool {
        match (&self.0, &other.0) {
            (RandomSourceKind::IntegerSeed(a), RandomSourceKind::IntegerSeed(b)) => a == b,
            (RandomSourceKind::RandomBits(a), RandomSourceKind::RandomBits(b)) => Arc::ptr_eq(a, b),
            _ => false,
        }
    }
}
impl Eq for StereoisomerRandomSource {}
impl Hash for StereoisomerRandomSource {
    fn hash<H: Hasher>(&self, state: &mut H) {
        std::mem::discriminant(&self.0).hash(state);
        match &self.0 {
            RandomSourceKind::IntegerSeed(seed) => seed.hash(state),
            RandomSourceKind::RandomBits(provider) => Arc::as_ptr(provider).hash(state),
        }
    }
}

struct SharedRandomBitsSource(StereoisomerRandomSource);
impl RandomBitsSource for SharedRandomBitsSource {
    fn getrandbits(&mut self, bit_count: usize) -> Result<ConfigurationBits, EnumerationError> {
        let RandomSourceKind::RandomBits(provider) = &self.0.0 else {
            unreachable!("shared source is constructed only for the provider representation");
        };
        let mut callback = provider
            .lock()
            .map_err(|_| EnumerationError::RandomSourcePoisoned(std::sync::PoisonError::new(())))?;
        callback(bit_count)
            .map(ConfigurationBits::from_biguint)
            .map_err(|message| EnumerationError::RandomBitsSource { bit_count, message })
    }
}

/// Options for lazy stereoisomer enumeration.
///
/// The defaults reproduce RDKit Python's `StereoEnumerationOptions`. A
/// `random_source` is consulted only when `max_isomers` selects the random
/// configuration source; `None` uses RDKit's molecule-invariant default seed.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct StereoisomerOptions {
    pub try_embedding: bool,
    pub only_unassigned: bool,
    pub only_stereo_groups: bool,
    pub max_isomers: usize,
    pub random_source: Option<StereoisomerRandomSource>,
    pub unique: bool,
}

impl Default for StereoisomerOptions {
    fn default() -> Self {
        // RDKit✔️✔️:   def __init__(self, tryEmbedding=False, onlyUnassigned=True, maxIsomers=1024, rand=None,
        // RDKit✔️✔️:                unique=True, onlyStereoGroups=False):
        // RDKit✔️✔️:     self.tryEmbedding = tryEmbedding
        // RDKit✔️✔️:     self.onlyUnassigned = onlyUnassigned
        // RDKit✔️✔️:     self.onlyStereoGroups = onlyStereoGroups
        // RDKit✔️✔️:     self.maxIsomers = maxIsomers
        // RDKit✔️✔️:     self.rand = rand
        // RDKit✔️✔️:     self.unique = unique
        Self {
            try_embedding: false,
            only_unassigned: true,
            only_stereo_groups: false,
            max_isomers: 1024,
            random_source: None,
            unique: true,
        }
    }
}

impl StereoisomerOptions {
    #[must_use]
    pub fn new(
        try_embedding: bool,
        only_unassigned: bool,
        max_isomers: usize,
        random_source: Option<StereoisomerRandomSource>,
        unique: bool,
        only_stereo_groups: bool,
    ) -> Self {
        // RDKit❗✔️:   def __init__(self, tryEmbedding=False, onlyUnassigned=True, maxIsomers=1024, rand=None,
        // RDKit❗✔️:                unique=True, onlyStereoGroups=False):
        // RDKit❗✔️:     self.tryEmbedding = tryEmbedding
        // RDKit❗✔️:     self.onlyUnassigned = onlyUnassigned
        // RDKit❗✔️:     self.maxIsomers = maxIsomers
        // RDKit❗✔️:     self.rand = rand
        // RDKit❗✔️:     self.unique = unique
        Self {
            try_embedding,
            only_unassigned,
            only_stereo_groups,
            max_isomers,
            random_source,
            unique,
        }
    }
    #[must_use]
    pub fn try_embedding(&self) -> bool {
        self.try_embedding
    }
    pub fn set_try_embedding(&mut self, value: bool) {
        self.try_embedding = value;
    }
    #[must_use]
    pub fn only_unassigned(&self) -> bool {
        self.only_unassigned
    }
    pub fn set_only_unassigned(&mut self, value: bool) {
        self.only_unassigned = value;
    }
    #[must_use]
    pub fn only_stereo_groups(&self) -> bool {
        self.only_stereo_groups
    }
    pub fn set_only_stereo_groups(&mut self, value: bool) {
        self.only_stereo_groups = value;
    }
    #[must_use]
    pub fn max_isomers(&self) -> usize {
        self.max_isomers
    }
    pub fn set_max_isomers(&mut self, value: usize) {
        self.max_isomers = value;
    }
    #[must_use]
    pub fn random_source(&self) -> Option<&StereoisomerRandomSource> {
        self.random_source.as_ref()
    }
    pub fn set_random_source(&mut self, value: Option<StereoisomerRandomSource>) {
        self.random_source = value;
    }
    #[must_use]
    pub fn unique(&self) -> bool {
        self.unique
    }
    pub fn set_unique(&mut self, value: bool) {
        self.unique = value;
    }
    fn flipper_selection(&self) -> FlipperSelectionOptions {
        FlipperSelectionOptions {
            only_unassigned: self.only_unassigned,
            only_stereo_groups: self.only_stereo_groups,
        }
    }
}

enum ConfigurationSource {
    Exhaustive(RangeBitsGenerator),
    Random(UniqueRandomBitsGenerator<PythonRandom>),
    Callback(Box<dyn Iterator<Item = Result<ConfigurationBits, EnumerationError>> + Send + Sync>),
}

impl Iterator for ConfigurationSource {
    type Item = Result<ConfigurationBits, EnumerationError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Exhaustive(source) => source.next().map(Ok),
            Self::Random(source) => source.next(),
            Self::Callback(source) => source.next(),
        }
    }
}

/// Lazy iterator over source-ordered stereoisomers.
///
/// Construction performs the source-defined preprocessing and candidate
/// discovery. Configuration application, uniqueness, optional embedding, and
/// their errors are deferred until `next()` requests an output.
pub struct StereoisomerIterator {
    workspace: Option<EnumerationWorkspace>,
    configurations: Option<ConfigurationSource>,
    pending_no_center: Option<SmilesRecord>,
    isomers_seen: HashSet<PropertyText>,
    unique: bool,
    try_embedding: bool,
    max_isomers: usize,
    yielded: usize,
    finished: bool,
}

impl StereoisomerIterator {
    fn new(
        molecule: &SmilesRecord,
        options: StereoisomerOptions,
        project_query: cosmolkit_conformer::PruningQueryProjection,
    ) -> Result<Self, EnumerationError> {
        Self::new_with_configuration_source(molecule, options, project_query, None)
    }

    fn new_with_configuration_source(
        molecule: &SmilesRecord,
        options: StereoisomerOptions,
        project_query: cosmolkit_conformer::PruningQueryProjection,
        callback: Option<Box<dyn FnMut(usize) -> Result<BigUint, String> + Send + Sync + 'static>>,
    ) -> Result<Self, EnumerationError> {
        // RDKit✔️✔️:   tm = Chem.Mol(m)
        // RDKit✔️✔️:   for atom in tm.GetAtoms():
        // RDKit✔️✔️:     atom.ClearProp("_CIPCode")
        // RDKit✔️✔️:   for bond in tm.GetBonds():
        // RDKit✔️✔️:     if bond.GetBondDir() == Chem.BondDir.EITHERDOUBLE or bond.GetBondDir() == Chem.BondDir.UNKNOWN:
        // RDKit✔️✔️:       bond.SetBondDir(Chem.BondDir.NONE)
        // RDKit✔️✔️:   flippers = _getFlippers(tm, options)
        // RDKit✔️✔️:   nCenters = len(flippers)
        let workspace =
            EnumerationWorkspace::prepare(molecule, options.flipper_selection(), project_query)?;
        let center_count = workspace.center_count();

        // RDKit✔️✔️:   if not nCenters:
        // RDKit✔️✔️:     yield tm
        // RDKit✔️✔️:     return
        if center_count == 0 {
            return Ok(Self {
                workspace: None,
                configurations: None,
                pending_no_center: Some(workspace.into_molecule()),
                isomers_seen: HashSet::new(),
                unique: options.unique,
                try_embedding: options.try_embedding,
                max_isomers: options.max_isomers,
                yielded: 0,
                finished: false,
            });
        }

        // RDKit✔️✔️:   if (options.maxIsomers == 0 or 2**nCenters <= options.maxIsomers):
        // RDKit✔️✔️:     bitsource = _RangeBitsGenerator(nCenters)
        // RDKit✔️✔️:   else:
        // RDKit✔️✔️:     if options.rand is None:
        // RDKit✔️✔️:       seed = hash(tuple(sorted([(a.GetDegree(), a.GetAtomicNum()) for a in tm.GetAtoms()])))
        // RDKit✔️✔️:       rand = random.Random(seed)
        // RDKit✔️✔️:     else:
        // RDKit✔️✔️:       rand = random.Random(options.rand)
        // RDKit✔️✔️:     bitsource = _UniqueRandomBitsGenerator(nCenters, options.maxIsomers, rand)
        let configuration_count = theoretical_configuration_count(center_count);
        let configurations = if options.max_isomers == 0
            || configuration_count <= BigUint::from(options.max_isomers)
        {
            ConfigurationSource::Exhaustive(RangeBitsGenerator::new(center_count))
        } else if let Some(callback) = callback {
            ConfigurationSource::Callback(Box::new(UniqueRandomBitsGenerator::new(
                center_count,
                options.max_isomers,
                CallbackRandomBitsSource { callback },
            )))
        } else {
            match options.random_source.as_ref() {
                Some(source @ StereoisomerRandomSource(RandomSourceKind::RandomBits(_))) => {
                    ConfigurationSource::Callback(Box::new(UniqueRandomBitsGenerator::new(
                        center_count,
                        options.max_isomers,
                        SharedRandomBitsSource(source.clone()),
                    )))
                }
                Some(StereoisomerRandomSource(RandomSourceKind::IntegerSeed(seed))) => {
                    ConfigurationSource::Random(UniqueRandomBitsGenerator::new(
                        center_count,
                        options.max_isomers,
                        PythonRandom::from_integer_seed(seed),
                    ))
                }
                None => {
                    let seed = default_python_random_seed(workspace.molecule());
                    ConfigurationSource::Random(UniqueRandomBitsGenerator::new(
                        center_count,
                        options.max_isomers,
                        PythonRandom::from_integer_seed(&seed),
                    ))
                }
            }
        };

        Ok(Self {
            workspace: Some(workspace),
            configurations: Some(configurations),
            pending_no_center: None,
            isomers_seen: HashSet::new(),
            unique: options.unique,
            try_embedding: options.try_embedding,
            max_isomers: options.max_isomers,
            yielded: 0,
            finished: false,
        })
    }

    #[must_use]
    pub fn yielded_count(&self) -> usize {
        self.yielded
    }
}

impl Iterator for StereoisomerIterator {
    type Item = Result<SmilesRecord, EnumerationError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        if let Some(isomer) = self.pending_no_center.take() {
            self.finished = true;
            self.yielded = 1;
            return Some(
                validate_candidate(&isomer)
                    .map(|()| isomer)
                    .map_err(EnumerationError::from),
            );
        }

        // RDKit✔️✔️:       if options.maxIsomers != 0 and numIsomers >= options.maxIsomers:
        // RDKit✔️✔️:         break
        if self.max_isomers != 0 && self.yielded >= self.max_isomers {
            self.finished = true;
            return None;
        }

        loop {
            // RDKit✔️✔️:   for bitflag in bitsource:
            let configuration = match self.configurations.as_mut()?.next() {
                Some(Ok(configuration)) => configuration,
                Some(Err(error)) => {
                    self.finished = true;
                    return Some(Err(error));
                }
                None => {
                    self.finished = true;
                    return None;
                }
            };

            let workspace = self
                .workspace
                .as_mut()
                .expect("enumeration workspace exists whenever configurations exist");
            if let Err(error) = workspace.apply_configuration(&configuration) {
                self.finished = true;
                return Some(Err(error));
            }

            let isomer = match workspace.finalize_configuration_with_embedding(
                &configuration,
                self.unique,
                self.try_embedding,
                &mut self.isomers_seen,
            ) {
                Ok(Some(isomer)) => isomer,
                Ok(None) => continue,
                Err(error) => {
                    self.finished = true;
                    return Some(Err(error));
                }
            };

            if let Err(error) = validate_candidate(&isomer) {
                self.finished = true;
                return Some(Err(error.into()));
            }

            // RDKit✔️✔️:     if cid >= 0:
            // RDKit✔️✔️:       yield isomer
            // RDKit✔️✔️:       numIsomers += 1
            self.yielded += 1;
            return Some(Ok(isomer));
        }
    }
}

impl std::iter::FusedIterator for StereoisomerIterator {}

/// Return RDKit Python's upper-bound stereoisomer count for `molecule`.
pub fn stereoisomer_count(
    molecule: &SmilesRecord,
    options: &StereoisomerOptions,
) -> Result<BigUint, EnumerationError> {
    // RDKit✔️✔️: def GetStereoisomerCount(m, options=StereoEnumerationOptions()):
    // RDKit✔️✔️:   tm = Chem.Mol(m)
    // RDKit✔️✔️:   flippers = _getFlippers(tm, options)
    // RDKit✔️✔️:   return 2**len(flippers)
    let mut workspace = molecule.clone();
    let flippers = select_stereo_flippers(&mut workspace, options.flipper_selection())?;
    Ok(theoretical_configuration_count(flippers.len()))
}

/// Create a lazy source-ordered stereoisomer iterator.
pub fn enumerate_stereoisomers(
    molecule: &SmilesRecord,
    options: StereoisomerOptions,
    project_query: cosmolkit_conformer::PruningQueryProjection,
) -> Result<StereoisomerIterator, EnumerationError> {
    StereoisomerIterator::new(molecule, options, project_query)
}

/// Create a lazy stereoisomer iterator using a source-compatible random-bit
/// callback when the option boundary selects bounded random enumeration.
///
/// Configuration uniqueness, finite exhaustion, output uniqueness, embedding,
/// and `max_isomers` accounting remain owned by the canonical iterator. The
/// callback supplies only Python-compatible `getrandbits(n)` values.
pub fn enumerate_stereoisomers_with_random_bits<F>(
    molecule: &SmilesRecord,
    options: StereoisomerOptions,
    callback: F,
    project_query: cosmolkit_conformer::PruningQueryProjection,
) -> Result<StereoisomerIterator, EnumerationError>
where
    F: FnMut(usize) -> Result<BigUint, String> + Send + Sync + 'static,
{
    StereoisomerIterator::new_with_configuration_source(
        molecule,
        options,
        project_query,
        Some(Box::new(callback)),
    )
}

fn validate_candidate(record: &SmilesRecord) -> Result<(), EnumerationError> {
    record.topology.validate()?;
    record
        .coordinates
        .validate_for_atom_count(record.topology.atoms.len())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::VecDeque;
    struct ScriptedRandomBits {
        values: VecDeque<Result<BigUint, String>>,
        requested_bit_counts: Vec<usize>,
    }

    impl ScriptedRandomBits {
        fn values(values: impl IntoIterator<Item = u64>) -> Self {
            Self {
                values: values
                    .into_iter()
                    .map(|value| Ok(BigUint::from(value)))
                    .collect(),
                requested_bit_counts: Vec::new(),
            }
        }

        fn error(message: impl Into<String>) -> Self {
            Self {
                values: VecDeque::from([Err(message.into())]),
                requested_bit_counts: Vec::new(),
            }
        }
    }

    impl RandomBitsSource for ScriptedRandomBits {
        fn getrandbits(&mut self, bit_count: usize) -> Result<ConfigurationBits, EnumerationError> {
            self.requested_bit_counts.push(bit_count);
            match self.values.pop_front() {
                Some(Ok(value)) => Ok(ConfigurationBits::from_biguint(value)),
                Some(Err(message)) => {
                    Err(EnumerationError::RandomBitsSource { bit_count, message })
                }
                None => Err(EnumerationError::RandomBitsSource {
                    bit_count,
                    message: "scripted source exhausted".to_owned(),
                }),
            }
        }
    }

    #[test]
    fn python_range_bits_generator_reproduces_zero_center_and_lsb_first_order() {
        let zero_center_values = RangeBitsGenerator::new(0)
            .map(|bits| bits.value().clone())
            .collect::<Vec<_>>();
        assert_eq!(zero_center_values, vec![BigUint::from(0_u8)]);

        let values = RangeBitsGenerator::new(3).collect::<Vec<_>>();
        assert_eq!(values.len(), 8);
        for (value, bits) in values.iter().enumerate() {
            assert_eq!(bits.value(), &BigUint::from(value));
            assert_eq!(bits.bit(0), value & 0b001 != 0);
            assert_eq!(bits.bit(1), value & 0b010 != 0);
            assert_eq!(bits.bit(2), value & 0b100 != 0);
        }
    }

    #[test]
    fn python_range_bits_generator_preserves_arbitrary_width_counts_and_laziness() {
        assert_eq!(
            theoretical_configuration_count(130),
            BigUint::from(1_u8) << 130
        );

        let prefix = RangeBitsGenerator::new(130)
            .take(4)
            .map(|bits| bits.value().clone())
            .collect::<Vec<_>>();
        assert_eq!(
            prefix,
            [0_u8, 1, 2, 3]
                .into_iter()
                .map(BigUint::from)
                .collect::<Vec<_>>()
        );

        // Constructing a million-center range allocates only its arbitrary-width
        // endpoint and current value; consuming a prefix never materializes 2^N rows.
        let huge_prefix = RangeBitsGenerator::new(1_000_000)
            .take(3)
            .map(|bits| bits.value().clone())
            .collect::<Vec<_>>();
        assert_eq!(
            huge_prefix,
            [0_u8, 1, 2]
                .into_iter()
                .map(BigUint::from)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn python_unique_random_bits_generator_rejects_duplicates_and_exhausts_finitely() {
        let random = ScriptedRandomBits::values([3, 1, 3, 0, 2]);
        let mut generator = UniqueRandomBitsGenerator::new(2, 1, random);
        let values = generator
            .by_ref()
            .map(|result| result.unwrap().value().clone())
            .collect::<Vec<_>>();
        assert_eq!(
            values,
            [3_u8, 1, 0, 2]
                .into_iter()
                .map(BigUint::from)
                .collect::<Vec<_>>()
        );
        assert_eq!(generator.random.requested_bit_counts, vec![2; 5]);
        assert!(generator.next().is_none());

        let mut zero_centers =
            UniqueRandomBitsGenerator::new(0, usize::MAX, ScriptedRandomBits::values([0, 0]));
        assert_eq!(
            zero_centers.next().unwrap().unwrap().value(),
            &BigUint::from(0_u8)
        );
        assert!(zero_centers.next().is_none());
        assert_eq!(zero_centers.random.requested_bit_counts, vec![0]);
    }

    #[test]
    fn python_unique_random_bits_generator_propagates_source_errors_structurally() {
        let mut generator =
            UniqueRandomBitsGenerator::new(130, 1024, ScriptedRandomBits::error("fixture failure"));
        assert!(matches!(generator.next().unwrap().unwrap_err(),
            EnumerationError::RandomBitsSource { bit_count: 130, message } if message == "fixture failure"));
        assert_eq!(generator.random.requested_bit_counts, vec![130]);
    }

    fn python_random_values(seed: &BigInt, widths: &[usize]) -> Vec<BigUint> {
        let mut random = PythonRandom::from_integer_seed(seed);
        widths
            .iter()
            .map(|&width| random.getrandbits(width).unwrap().value().clone())
            .collect()
    }

    fn decimal_biguint(value: &str) -> BigUint {
        value.parse().unwrap()
    }

    #[test]
    fn python_random_matches_cpython_integer_seed_and_getrandbits_sequence_exactly() {
        let widths = [0, 1, 2, 31, 32, 33, 64, 65, 130];
        let values = python_random_values(&BigInt::from(0xdead_beef_u32), &widths);
        let expected = [
            "0",
            "0",
            "1",
            "1008527992",
            "3864962434",
            "6868874809",
            "16491411380257451511",
            "25412851181257523183",
            "668769172828971426222707781417953853951",
        ]
        .map(decimal_biguint);
        assert_eq!(values, expected);
    }

    #[test]
    fn python_random_uses_all_integer_seed_bits_and_absolute_negative_seed_semantics() {
        let seed = (BigInt::from(1_u8) << 100) + BigInt::from(0x1234_5678_u32);
        let widths = [0, 1, 2, 31, 32, 33, 64, 65, 130];
        let positive = python_random_values(&seed, &widths);
        let negative = python_random_values(&(-seed), &widths);
        let expected = [
            "0",
            "1",
            "0",
            "552449039",
            "185676661",
            "7120073926",
            "18340008931984535783",
            "34270592875609026320",
            "1287945476240184999319519200333061618674",
        ]
        .map(decimal_biguint);
        assert_eq!(positive, expected);
        assert_eq!(negative, expected);
    }

    #[test]
    fn python_default_random_seed_is_exact_and_invariant_to_input_atom_order() {
        let cco = cosmolkit_smiles::parse_smiles("CCO", &Default::default()).unwrap();
        let occ = cosmolkit_smiles::parse_smiles("OCC", &Default::default()).unwrap();
        assert_ne!(
            cco.topology
                .atoms
                .iter()
                .map(|atom| atom.atomic_number())
                .collect::<Vec<_>>(),
            occ.topology
                .atoms
                .iter()
                .map(|atom| atom.atomic_number())
                .collect::<Vec<_>>()
        );

        let cco_seed = default_python_random_seed(&cco);
        let occ_seed = default_python_random_seed(&occ);
        assert_eq!(cco_seed, BigInt::from(1_554_427_866_021_819_617_i64));
        assert_eq!(cco_seed, occ_seed);

        let expected = [1_u8, 7, 7, 5, 3, 1, 4, 2, 2, 0, 3, 3]
            .into_iter()
            .map(BigUint::from)
            .collect::<Vec<_>>();
        let first = python_random_values(&cco_seed, &[3; 12]);
        let repeated = python_random_values(&cco_seed, &[3; 12]);
        assert_eq!(first, expected);
        assert_eq!(repeated, expected);
    }

    struct CounterRandom {
        next_value: BigUint,
    }

    impl RandomBitsSource for CounterRandom {
        fn getrandbits(&mut self, bit_count: usize) -> Result<ConfigurationBits, EnumerationError> {
            let value = self.next_value.clone();
            self.next_value += 1_u8;
            if value >= theoretical_configuration_count(bit_count) {
                return Err(EnumerationError::RandomBitsSource {
                    bit_count,
                    message: "counter random source exhausted the declared bit width".to_owned(),
                });
            }
            Ok(ConfigurationBits::from_biguint(value))
        }
    }

    #[test]
    fn python_random_custom_getrandbits_boundary_matches_rdkit_counter_fixture() {
        let random = CounterRandom {
            next_value: BigUint::from(0_u8),
        };
        let values = UniqueRandomBitsGenerator::new(3, 3, random)
            .take(3)
            .map(|result| result.unwrap().value().clone())
            .collect::<Vec<_>>();
        assert_eq!(
            values,
            [0_u8, 1, 2]
                .into_iter()
                .map(BigUint::from)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn python_random_instances_are_parallel_and_repeat_call_isolated() {
        let seed = BigInt::from(-5_473_470_236_788_694_370_i64);
        let expected = [6_u8, 1, 4, 6, 4, 2, 0, 7, 6, 5, 2, 6]
            .into_iter()
            .map(BigUint::from)
            .collect::<Vec<_>>();
        let outputs = std::thread::scope(|scope| {
            (0..8)
                .map(|_| {
                    let seed = &seed;
                    scope.spawn(move || python_random_values(seed, &[3; 12]))
                })
                .collect::<Vec<_>>()
                .into_iter()
                .map(|thread| thread.join().unwrap())
                .collect::<Vec<_>>()
        });
        assert_eq!(outputs, vec![expected; 8]);
    }
    use cosmolkit_model::{
        AtomQueryPredicate, BondQueryPredicate, CoordinateBlock, MoleculeProperties, QueryAtom,
        QueryBond, QueryGraph, QueryGraphError, QueryNode, TopologyBlock,
    };
    fn record_from_smiles(input: &str) -> Result<SmilesRecord, Box<dyn std::error::Error>> {
        let params = cosmolkit_smiles::SmilesParseParams::default();
        let mut record = cosmolkit_smiles::parse_smiles(input, &params)?;
        let prepared = cosmolkit_core::sanitize_topology(&record.topology, &Default::default())?;
        record.topology = prepared.topology;
        let mut valence = prepared.final_valence;
        let mut rings = prepared.final_rings;
        Ok(cosmolkit_smiles::finalize_smiles_stereo(
            record,
            &params,
            &mut valence,
            &mut rings,
        )?)
    }
    fn text(value: PropertyText) -> String {
        std::str::from_utf8(value.as_bytes())
            .expect("source fixture UTF8")
            .to_owned()
    }

    fn project_query(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        properties: &MoleculeProperties,
    ) -> Result<QueryGraph, QueryGraphError> {
        let atoms = topology
            .atoms
            .iter()
            .map(|atom| {
                QueryAtom::from_carrier_parts(
                    atom.clone(),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number())),
                )
            })
            .collect();
        let bonds = topology
            .bonds
            .iter()
            .map(|bond| {
                QueryBond::from_carrier_parts(
                    bond.clone(),
                    QueryNode::predicate(BondQueryPredicate::Order(bond.order())),
                )
            })
            .collect();
        let mut props = properties.props().clone();
        if let Some(name) = properties.name() {
            props.insert(
                "_Name".into(),
                cosmolkit_model::PropertyValue::String(name.clone()),
            );
        }
        let mut query = QueryGraph::from_parts(
            atoms,
            bonds,
            props,
            coordinates.conformers_2d.clone(),
            coordinates.conformers_3d.clone(),
            topology.stereo_groups.clone(),
        )?;
        query.set_source_conformer_order(coordinates.source_conformer_order.clone())?;
        cosmolkit_model::replace_query_substance_groups(
            &mut query,
            topology.substance_groups.clone(),
        )?;
        Ok(query)
    }
    #[test]
    fn python_enumerator_finalization_continues_on_meso_canonical_duplicates_in_source_order() {
        let source = record_from_smiles("FC(Cl)C(Cl)F").unwrap();
        let mut workspace = EnumerationWorkspace::prepare(
            &source,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: false,
            },
            project_query,
        )
        .unwrap();
        assert_eq!(workspace.center_count(), 2);

        let mut seen = HashSet::new();
        let mut outputs = Vec::new();
        let mut emitted_by_configuration = Vec::new();
        for configuration in 0_u8..4 {
            workspace
                .apply_configuration(&ConfigurationBits::from_biguint(BigUint::from(
                    configuration,
                )))
                .unwrap();
            let finalized = workspace.finalize_configuration(true, &mut seen).unwrap();
            emitted_by_configuration.push(finalized.is_some());
            if let Some(isomer) = finalized {
                outputs.push(cosmolkit_smiles::write_smiles(&isomer).map(text).unwrap());
            }
        }

        assert_eq!(emitted_by_configuration, vec![true, true, true, false]);
        assert_eq!(
            outputs,
            vec![
                "F[C@H](Cl)[C@@H](F)Cl",
                "F[C@@H](Cl)[C@@H](F)Cl",
                "F[C@H](Cl)[C@H](F)Cl",
            ]
        );
        assert_eq!(seen.len(), 3);
    }

    #[test]
    fn python_enumerator_preserves_tetrahedral_center_when_double_bond_configurations_change() {
        let source =
            record_from_smiles("Cl.N=C(N)c1ccc(C=CC(=O)NCC(O)CNC(=O)C=Cc2ccc(C(=N)N)cc2)cc1")
                .unwrap();
        let source_before = source.clone();

        let non_unique = enumerate_stereoisomers(
            &source,
            StereoisomerOptions {
                max_isomers: 0,
                unique: false,
                ..Default::default()
            },
            project_query,
        )
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
        assert_eq!(non_unique.len(), 8);
        for isomer in &non_unique[2..6] {
            assert!(
                isomer.topology.atoms[14].prop("_CIPCode").is_some(),
                "enumeration finalization must retain the source-assigned CIP code"
            );
        }
        assert_eq!(
            non_unique
                .iter()
                .map(|isomer| isomer.topology.atoms[14].chiral_tag())
                .collect::<Vec<_>>(),
            vec![
                ChiralTag::Unspecified,
                ChiralTag::Unspecified,
                ChiralTag::TetrahedralCcw,
                ChiralTag::TetrahedralCw,
                ChiralTag::TetrahedralCcw,
                ChiralTag::TetrahedralCw,
                ChiralTag::Unspecified,
                ChiralTag::Unspecified,
            ]
        );

        let unique = enumerate_stereoisomers(
            &source,
            StereoisomerOptions {
                max_isomers: 0,
                unique: true,
                ..Default::default()
            },
            project_query,
        )
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
        assert_eq!(
            unique
                .iter()
                .map(|isomer| { cosmolkit_smiles::write_smiles(isomer).map(text) })
                .collect::<Result<Vec<_>, _>>()
                .unwrap(),
            vec![
                "Cl.N=C(N)c1ccc(/C=C/C(=O)NCC(O)CNC(=O)/C=C/c2ccc(C(=N)N)cc2)cc1",
                "Cl.N=C(N)c1ccc(/C=C\\C(=O)NC[C@H](O)CNC(=O)/C=C/c2ccc(C(=N)N)cc2)cc1",
                "Cl.N=C(N)c1ccc(/C=C\\C(=O)NC[C@@H](O)CNC(=O)/C=C/c2ccc(C(=N)N)cc2)cc1",
                "Cl.N=C(N)c1ccc(/C=C\\C(=O)NCC(O)CNC(=O)/C=C\\c2ccc(C(=N)N)cc2)cc1",
            ]
        );
        assert_eq!(source, source_before);
    }

    #[test]
    fn python_enumerator_try_embedding_masks_configuration_to_signed_cxx_seed() {
        let cases = [
            (0_u64, 0_i32),
            (1, 1),
            (0x7fff_ffff, 0x7fff_ffff),
            (0x8000_0000, 0),
            (0xffff_ffff, 0x7fff_ffff),
            (0x1_0000_0005, 5),
        ];
        for (configuration, expected) in cases {
            assert_eq!(
                ConfigurationBits::from_biguint(BigUint::from(configuration)).embedding_seed(),
                expected,
                "configuration {configuration:#x}"
            );
        }
    }

    #[test]
    fn python_enumerator_try_embedding_matches_rdkit_rejection_map_and_exhausts_finitely() {
        // RDKit UnitTestMol3D.py:
        // testEnumerateStereoisomersMaxIsomersShouldBeReturnedEvenWithTryEmbedding
        // testEnumerateStereoisomersTryEmbeddingShouldNotInfiniteLoopWhenMaxIsomersIsLargerThanActual
        let source = record_from_smiles("BrC=CC1OC(C2)(F)C2(Cl)C1").unwrap();
        let source_before = source.clone();
        let mut workspace = EnumerationWorkspace::prepare(
            &source,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: false,
            },
            project_query,
        )
        .unwrap();
        assert_eq!(workspace.center_count(), 4);

        let mut seen = HashSet::new();
        let mut succeeded = Vec::new();
        let mut rejected = Vec::new();
        let mut attempts = 0_usize;
        for configuration in RangeBitsGenerator::new(workspace.center_count()) {
            attempts += 1;
            workspace.apply_configuration(&configuration).unwrap();
            match workspace
                .finalize_configuration_with_embedding(&configuration, true, true, &mut seen)
                .unwrap()
            {
                Some(isomer) => {
                    assert_eq!(isomer.coordinates.conformers_3d.len(), 1);
                    assert_eq!(
                        isomer.coordinates.conformers_3d[0].coordinates().len(),
                        source.topology.atoms.len()
                    );
                    succeeded.push(configuration.value().clone());
                }
                None => rejected.push(configuration.value().clone()),
            }
        }

        assert_eq!(attempts, 16);
        assert_eq!(
            succeeded,
            [0_u8, 1, 6, 7, 8, 9, 14, 15]
                .into_iter()
                .map(BigUint::from)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            rejected,
            [2_u8, 3, 4, 5, 10, 11, 12, 13]
                .into_iter()
                .map(BigUint::from)
                .collect::<Vec<_>>()
        );
        assert_eq!(succeeded.len(), 8);
        assert_eq!(source, source_before);

        // A max of two counts successful outputs, not attempted or rejected
        // configurations. Configurations 0 and 1 are the first two successes.
        let mut seen = HashSet::new();
        let mut emitted = Vec::new();
        let mut attempts = 0_usize;
        for configuration in RangeBitsGenerator::new(workspace.center_count()) {
            attempts += 1;
            workspace.apply_configuration(&configuration).unwrap();
            if workspace
                .finalize_configuration_with_embedding(&configuration, true, true, &mut seen)
                .unwrap()
                .is_some()
            {
                emitted.push(configuration.value().clone());
                if emitted.len() == 2 {
                    break;
                }
            }
        }
        assert_eq!(attempts, 2);
        assert_eq!(emitted, vec![BigUint::from(0_u8), BigUint::from(1_u8)]);
    }
    use cosmolkit_core::PotentialStereoDescriptor;
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomSpec, Bond, BondOrder, BondSpec, Element, StereoGroup,
        TopologyValidationError,
    };
    fn empty_record() -> SmilesRecord {
        SmilesRecord {
            topology: TopologyBlock::default(),
            coordinates: CoordinateBlock::default(),
            properties: MoleculeProperties::default(),
        }
    }
    fn fixture_add_atom(record: &mut SmilesRecord, spec: AtomSpec) -> AtomId {
        let id = AtomId::new(record.topology.atoms.len());
        record.topology.atoms.push(Atom::from_spec(id, spec));
        record.topology.adjacency =
            AdjacencyList::from_topology(record.topology.atoms.len(), &record.topology.bonds);
        id
    }
    fn fixture_add_bond(
        record: &mut SmilesRecord,
        spec: BondSpec,
    ) -> Result<BondId, cosmolkit_model::AdjacencyError> {
        let id = BondId::new(record.topology.bonds.len());
        record.topology.bonds.push(Bond::from_spec(id, spec));
        record.topology.adjacency =
            AdjacencyList::try_from_topology(record.topology.atoms.len(), &record.topology.bonds)?;
        Ok(id)
    }
    fn fixture_add_stereo_group(
        record: &mut SmilesRecord,
        group: StereoGroup,
    ) -> Result<usize, TopologyValidationError> {
        let id = record.topology.stereo_groups.len();
        record.topology.stereo_groups.push(group);
        Ok(id)
    }
    fn fixture_add_conformer(
        record: &mut SmilesRecord,
        conformer: Conformer3D,
    ) -> Result<(), cosmolkit_model::CoordinateValidationError> {
        record
            .coordinates
            .record_source_conformer_append(CoordinateDimension::ThreeD)?;
        record.coordinates.conformers_3d.push(conformer);
        record
            .coordinates
            .validate_for_atom_count(record.topology.atoms.len())
    }
    fn checked_record(record: SmilesRecord) -> SmilesRecord {
        validate_candidate(&record).unwrap();
        record
    }
    fn stereo_info(
        stereo_type: PotentialStereoType,
        specified: PotentialStereoSpecified,
        center: PotentialStereoCenter,
        controlling_atoms: Vec<Option<AtomId>>,
    ) -> PotentialStereoInfo {
        PotentialStereoInfo {
            stereo_type,
            specified,
            centered_on: center,
            descriptor: PotentialStereoDescriptor::None,
            permutation: 0,
            controlling_atoms,
        }
    }

    fn tetrahedral_center(
        builder: &mut SmilesRecord,
        tag: ChiralTag,
        cip_code: Option<&str>,
    ) -> AtomId {
        let mut center = AtomSpec::new(Element::C)
            .with_no_implicit(true)
            .with_chiral_tag(tag)
            .with_prop("keep_atom", "center")
            .unwrap();
        if let Some(cip_code) = cip_code {
            center = center.with_prop("_CIPCode", cip_code).unwrap();
        }
        let center = fixture_add_atom(builder, center);
        for element in [Element::F, Element::CL, Element::BR, Element::I] {
            let ligand = fixture_add_atom(builder, AtomSpec::new(element));
            fixture_add_bond(builder, BondSpec::new(center, ligand, BondOrder::Single)).unwrap();
        }
        center
    }

    fn grouped_selection_molecule() -> (SmilesRecord, AtomId, AtomId, AtomId, AtomId, BondId) {
        let mut builder = empty_record();
        let assigned = fixture_add_atom(
            &mut builder,
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
        );
        let unassigned = fixture_add_atom(&mut builder, AtomSpec::new(Element::C));
        let begin_controller = fixture_add_atom(
            &mut builder,
            AtomSpec::new(Element::F).with_chiral_tag(ChiralTag::TetrahedralCcw),
        );
        let end_controller = fixture_add_atom(&mut builder, AtomSpec::new(Element::CL));
        let double_bond = fixture_add_bond(
            &mut builder,
            BondSpec::new(assigned, unassigned, BondOrder::Double),
        )
        .unwrap();
        fixture_add_bond(
            &mut builder,
            BondSpec::new(assigned, begin_controller, BondOrder::Single),
        )
        .unwrap();
        fixture_add_bond(
            &mut builder,
            BondSpec::new(unassigned, end_controller, BondOrder::Single),
        )
        .unwrap();
        fixture_add_stereo_group(
            &mut builder,
            StereoGroup::new(StereoGroupKind::Absolute, vec![assigned], Vec::new()),
        )
        .unwrap();
        fixture_add_stereo_group(
            &mut builder,
            StereoGroup::new(StereoGroupKind::Or, vec![assigned, unassigned], Vec::new()),
        )
        .unwrap();
        fixture_add_stereo_group(
            &mut builder,
            StereoGroup::new(StereoGroupKind::And, vec![begin_controller], Vec::new()),
        )
        .unwrap();
        (
            checked_record(builder),
            assigned,
            unassigned,
            begin_controller,
            end_controller,
            double_bond,
        )
    }

    #[test]
    fn python_enumerator_preprocessing_preserves_no_center_source_state_and_owned_blocks() {
        let mut builder = empty_record();
        builder.properties = builder.properties.with_name("preprocessing fixture");
        builder.properties.set_prop("keep_molecule", "yes").unwrap();
        builder
            .properties
            .set_prop("_MolFileChiralFlag", "0")
            .unwrap();
        let atoms = (0..4)
            .map(|index| {
                let spec = if index == 0 {
                    AtomSpec::new(Element::C)
                        .with_prop("_CIPCode", "R")
                        .unwrap()
                        .with_prop("keep_atom", "yes")
                        .unwrap()
                } else {
                    AtomSpec::new(Element::C)
                };
                fixture_add_atom(&mut builder, spec)
            })
            .collect::<Vec<_>>();
        for (index, direction) in [
            BondDirection::EitherDouble,
            BondDirection::Unknown,
            BondDirection::BeginWedge,
        ]
        .into_iter()
        .enumerate()
        {
            fixture_add_bond(
                &mut builder,
                BondSpec::new(atoms[index], atoms[index + 1], BondOrder::Single)
                    .with_direction(direction)
                    .with_prop("keep_bond", index.to_string())
                    .unwrap(),
            )
            .unwrap();
        }
        let coordinates = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [3.0, 1.0, 1.0],
        ];
        fixture_add_conformer(
            &mut builder,
            Conformer3D::new(7, coordinates.clone(), true).with_prop("keep_conf", "yes"),
        )
        .unwrap();
        let source = checked_record(builder);
        let source_before = source.clone();

        let workspace = EnumerationWorkspace::prepare(
            &source,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: false,
            },
            project_query,
        )
        .unwrap();

        assert_eq!(source, source_before);
        assert_eq!(workspace.center_count(), 0);
        assert_eq!(
            workspace.molecule().topology.atoms[0].prop("_CIPCode"),
            None
        );
        assert_eq!(
            workspace.molecule().topology.atoms[0].prop("keep_atom"),
            Some(&cosmolkit_model::PropertyValue::from("yes"))
        );
        assert_eq!(
            workspace
                .molecule()
                .topology
                .bonds
                .iter()
                .map(cosmolkit_model::Bond::direction)
                .collect::<Vec<_>>(),
            vec![
                BondDirection::None,
                BondDirection::None,
                BondDirection::BeginWedge,
            ]
        );
        assert_eq!(
            workspace.molecule().topology.bonds[2].prop("keep_bond"),
            Some(&cosmolkit_model::PropertyValue::from("2"))
        );
        assert_eq!(
            workspace.molecule().properties.prop("keep_molecule"),
            Some(&cosmolkit_model::PropertyValue::from("yes"))
        );
        assert_eq!(
            workspace.molecule().properties.prop("_MolFileChiralFlag"),
            Some(&cosmolkit_model::PropertyValue::from("0"))
        );
        assert_eq!(workspace.molecule().coordinates.conformers_3d.len(), 1);
        assert_eq!(
            workspace.molecule().coordinates.conformers_3d[0].coordinates(),
            coordinates
        );
        assert_eq!(
            workspace.molecule().coordinates.conformers_3d[0]
                .props()
                .get(b"keep_conf".as_slice())
                .map(|value| std::str::from_utf8(value.as_bytes()).unwrap()),
            Some("yes")
        );
    }

    #[test]
    fn python_enumerator_preprocessing_filters_fully_assigned_and_partial_centers_exactly() {
        let mut assigned_builder = empty_record();
        let assigned =
            tetrahedral_center(&mut assigned_builder, ChiralTag::TetrahedralCw, Some("R"));
        let assigned_source = checked_record(assigned_builder);
        let assigned_before = assigned_source.clone();
        let assigned_only = EnumerationWorkspace::prepare(
            &assigned_source,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: false,
            },
            project_query,
        )
        .unwrap();
        assert_eq!(assigned_only.center_count(), 0);
        assert_eq!(
            assigned_only.molecule().topology.atoms[assigned.index()].prop("_CIPCode"),
            None
        );
        assert_eq!(assigned_source, assigned_before);

        let mut partial_builder = empty_record();
        let retained =
            tetrahedral_center(&mut partial_builder, ChiralTag::TetrahedralCw, Some("S"));
        let enumerated =
            tetrahedral_center(&mut partial_builder, ChiralTag::Unspecified, Some("stale"));
        let partial_source = checked_record(partial_builder);
        let partial_before = partial_source.clone();
        let mut partial = EnumerationWorkspace::prepare(
            &partial_source,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: false,
            },
            project_query,
        )
        .unwrap();

        assert_eq!(partial.center_count(), 1);
        assert_eq!(
            partial.flippers,
            vec![StereoFlipper::Atom { atom: enumerated }]
        );
        assert_eq!(
            partial.molecule().properties.prop("_MolFileChiralFlag"),
            Some(&cosmolkit_model::PropertyValue::from("1"))
        );
        partial
            .apply_configuration(&ConfigurationBits::from_biguint(BigUint::from(0_u8)))
            .unwrap();
        assert_eq!(
            partial.molecule().topology.atoms[retained.index()].chiral_tag(),
            ChiralTag::TetrahedralCw
        );
        assert_eq!(
            partial.molecule().topology.atoms[enumerated.index()].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        partial
            .apply_configuration(&ConfigurationBits::from_biguint(BigUint::from(1_u8)))
            .unwrap();
        assert_eq!(
            partial.molecule().topology.atoms[enumerated.index()].chiral_tag(),
            ChiralTag::TetrahedralCw
        );
        assert_eq!(partial_source, partial_before);
    }

    #[test]
    fn python_enumerator_preprocessing_clears_either_double_and_applies_stereoany_bond_bits() {
        let mut builder = empty_record();
        let begin = fixture_add_atom(
            &mut builder,
            AtomSpec::new(Element::C).with_no_implicit(true),
        );
        let end = fixture_add_atom(
            &mut builder,
            AtomSpec::new(Element::C).with_no_implicit(true),
        );
        let begin_first = fixture_add_atom(&mut builder, AtomSpec::new(Element::F));
        let begin_second = fixture_add_atom(&mut builder, AtomSpec::new(Element::CL));
        let end_first = fixture_add_atom(&mut builder, AtomSpec::new(Element::BR));
        let end_second = fixture_add_atom(&mut builder, AtomSpec::new(Element::I));
        let double_bond = fixture_add_bond(
            &mut builder,
            BondSpec::new(begin, end, BondOrder::Double)
                .with_stereo(BondStereo::Any)
                .with_direction(BondDirection::EitherDouble),
        )
        .unwrap();
        fixture_add_bond(
            &mut builder,
            BondSpec::new(begin, begin_first, BondOrder::Single)
                .with_direction(BondDirection::Unknown),
        )
        .unwrap();
        fixture_add_bond(
            &mut builder,
            BondSpec::new(begin, begin_second, BondOrder::Single),
        )
        .unwrap();
        fixture_add_bond(
            &mut builder,
            BondSpec::new(end, end_first, BondOrder::Single),
        )
        .unwrap();
        fixture_add_bond(
            &mut builder,
            BondSpec::new(end, end_second, BondOrder::Single),
        )
        .unwrap();
        let source = checked_record(builder);
        let source_before = source.clone();
        let mut workspace = EnumerationWorkspace::prepare(
            &source,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: false,
            },
            project_query,
        )
        .unwrap();

        assert_eq!(workspace.center_count(), 1);
        assert_eq!(
            workspace.flippers,
            vec![StereoFlipper::Bond { bond: double_bond }]
        );
        assert_eq!(
            workspace.molecule().topology.bonds[double_bond.index()].direction(),
            BondDirection::None
        );
        assert_eq!(
            workspace.molecule().topology.bonds[1].direction(),
            BondDirection::None
        );
        assert_eq!(
            workspace.molecule().topology.bonds[double_bond.index()].stereo(),
            BondStereo::Any
        );
        assert!(
            workspace.molecule().topology.bonds[double_bond.index()]
                .stereo_atoms()
                .is_some()
        );
        workspace
            .apply_configuration(&ConfigurationBits::from_biguint(BigUint::from(0_u8)))
            .unwrap();
        assert_eq!(
            workspace.molecule().topology.bonds[double_bond.index()].stereo(),
            BondStereo::Trans
        );
        workspace
            .apply_configuration(&ConfigurationBits::from_biguint(BigUint::from(1_u8)))
            .unwrap();
        assert_eq!(
            workspace.molecule().topology.bonds[double_bond.index()].stereo(),
            BondStereo::Cis
        );
        assert_eq!(source, source_before);
    }

    #[test]
    fn python_enumerator_preprocessing_selects_enhanced_stereo_groups_once_in_group_order() {
        let mut builder = empty_record();
        let first = tetrahedral_center(&mut builder, ChiralTag::TetrahedralCw, Some("R"));
        let second = tetrahedral_center(&mut builder, ChiralTag::TetrahedralCcw, Some("S"));
        fixture_add_stereo_group(
            &mut builder,
            StereoGroup::new(StereoGroupKind::Absolute, vec![first], Vec::new()),
        )
        .unwrap();
        fixture_add_stereo_group(
            &mut builder,
            StereoGroup::new(StereoGroupKind::Or, vec![first, second], Vec::new()),
        )
        .unwrap();
        let source = checked_record(builder);
        let source_before = source.clone();
        let mut workspace = EnumerationWorkspace::prepare(
            &source,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: true,
            },
            project_query,
        )
        .unwrap();

        assert_eq!(workspace.center_count(), 1);
        assert_eq!(
            workspace.molecule().topology.stereo_groups,
            source.topology.stereo_groups
        );
        assert!(matches!(
            &workspace.flippers[0],
            StereoFlipper::StereoGroup { original_parities }
                if original_parities == &vec![
                    (first, ChiralTag::TetrahedralCw),
                    (second, ChiralTag::TetrahedralCcw),
                ]
        ));
        workspace
            .apply_configuration(&ConfigurationBits::from_biguint(BigUint::from(0_u8)))
            .unwrap();
        assert_eq!(
            workspace.molecule().topology.atoms[first.index()].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(
            workspace.molecule().topology.atoms[second.index()].chiral_tag(),
            ChiralTag::TetrahedralCw
        );
        workspace
            .apply_configuration(&ConfigurationBits::from_biguint(BigUint::from(1_u8)))
            .unwrap();
        assert_eq!(
            workspace.molecule().topology.atoms[first.index()].chiral_tag(),
            ChiralTag::TetrahedralCw
        );
        assert_eq!(
            workspace.molecule().topology.atoms[second.index()].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(source, source_before);
    }

    #[test]
    fn python_enumerator_finalization_clears_computed_state_preserves_rings_and_removes_groups() {
        let mut source = record_from_smiles("FC(Cl)C(Cl)F").unwrap();
        source
            .properties
            .set_prop("keep_molecule", "persistent")
            .unwrap();
        source
            .properties
            .set_computed_prop("drop_molecule", "computed")
            .unwrap();
        source.topology.atoms[1]
            .set_prop("keep_atom", "persistent")
            .unwrap();
        source.topology.atoms[1]
            .set_computed_prop("drop_atom", "computed")
            .unwrap();
        source.topology.bonds[0]
            .set_prop("keep_bond", "persistent")
            .unwrap();
        source.topology.bonds[0]
            .set_computed_prop("drop_bond", "computed")
            .unwrap();
        source.topology.stereo_groups.push(StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(1), AtomId::new(3)],
            Vec::new(),
        ));
        let source_before = source.clone();
        let source_rings =
            cosmolkit_core::symmetrized_sssr(&source.topology, &Default::default()).unwrap();

        let mut workspace = EnumerationWorkspace::prepare(
            &source,
            FlipperSelectionOptions {
                only_unassigned: false,
                only_stereo_groups: false,
            },
            project_query,
        )
        .unwrap();
        workspace
            .apply_configuration(&ConfigurationBits::from_biguint(BigUint::from(0_u8)))
            .unwrap();
        let output = workspace
            .finalize_configuration(false, &mut HashSet::new())
            .unwrap()
            .unwrap();

        assert_eq!(source, source_before);
        assert!(output.topology.stereo_groups.is_empty());
        assert_eq!(
            cosmolkit_core::symmetrized_sssr(&output.topology, &Default::default()).unwrap(),
            source_rings
        );
        assert_eq!(
            output.properties.prop("keep_molecule"),
            Some(&cosmolkit_model::PropertyValue::from("persistent"))
        );
        assert_eq!(output.properties.prop("drop_molecule"), None);
        assert_eq!(
            output.topology.atoms[1].prop("keep_atom"),
            Some(&cosmolkit_model::PropertyValue::from("persistent"))
        );
        assert_eq!(output.topology.atoms[1].prop("drop_atom"), None);
        assert_eq!(
            output.topology.bonds[0].prop("keep_bond"),
            Some(&cosmolkit_model::PropertyValue::from("persistent"))
        );
        assert_eq!(output.topology.bonds[0].prop("drop_bond"), None);
        assert_eq!(
            output.properties.prop("_MolFileChiralFlag"),
            Some(&cosmolkit_model::PropertyValue::from("1"))
        );
        assert_eq!(
            output.properties.prop("_StereochemDone"),
            Some(&cosmolkit_model::PropertyValue::Bool(true))
        );
        assert!(
            output
                .properties
                .is_prop_computed("_StereochemDone")
                .unwrap()
        );
        for center in [1, 3] {
            assert!(output.topology.atoms[center].prop("_CIPCode").is_some());
            assert_eq!(
                output.topology.atoms[center].prop("_ChiralityPossible"),
                Some(&cosmolkit_model::PropertyValue::Int(1))
            );
        }
    }

    #[test]
    fn python_enumerator_finalization_sets_double_bond_directions_and_handles_cumulenes() {
        let mut builder = empty_record();
        let begin = fixture_add_atom(
            &mut builder,
            AtomSpec::new(Element::C).with_no_implicit(true),
        );
        let end = fixture_add_atom(
            &mut builder,
            AtomSpec::new(Element::C).with_no_implicit(true),
        );
        let begin_first = fixture_add_atom(&mut builder, AtomSpec::new(Element::F));
        let begin_second = fixture_add_atom(&mut builder, AtomSpec::new(Element::CL));
        let end_first = fixture_add_atom(&mut builder, AtomSpec::new(Element::BR));
        let end_second = fixture_add_atom(&mut builder, AtomSpec::new(Element::I));
        let double_bond = fixture_add_bond(
            &mut builder,
            BondSpec::new(begin, end, BondOrder::Double)
                .with_stereo(BondStereo::Any)
                .with_stereo_atoms(begin_first, end_first),
        )
        .unwrap();
        for (center, ligand) in [
            (begin, begin_first),
            (begin, begin_second),
            (end, end_first),
            (end, end_second),
        ] {
            fixture_add_bond(
                &mut builder,
                BondSpec::new(center, ligand, BondOrder::Single),
            )
            .unwrap();
        }
        let source = checked_record(builder);
        let mut workspace = EnumerationWorkspace::prepare(
            &source,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: false,
            },
            project_query,
        )
        .unwrap();
        workspace
            .apply_configuration(&ConfigurationBits::from_biguint(BigUint::from(0_u8)))
            .unwrap();
        let trans = workspace
            .finalize_configuration(false, &mut HashSet::new())
            .unwrap()
            .unwrap();
        assert_eq!(
            trans.topology.bonds[double_bond.index()].stereo(),
            BondStereo::E
        );
        assert!(trans.topology.bonds.iter().any(|bond| matches!(
            bond.direction(),
            BondDirection::EndDownRight | BondDirection::EndUpRight
        )));

        let cumulene = record_from_smiles("CC=C=CC").unwrap();
        let cumulene_workspace = EnumerationWorkspace {
            molecule: cumulene,
            flippers: Vec::new(),
            project_query,
        };
        let cumulene_output = cumulene_workspace
            .finalize_configuration(false, &mut HashSet::new())
            .unwrap()
            .unwrap();
        assert_eq!(
            cosmolkit_smiles::write_smiles(&cumulene_output)
                .map(text)
                .unwrap(),
            "CC=C=CC"
        );
        assert_eq!(
            cumulene_output.properties.prop("_StereochemDone"),
            Some(&cosmolkit_model::PropertyValue::Bool(true))
        );
    }

    #[test]
    fn python_enumerator_try_embedding_preserves_existing_conformers_and_copies_heavy_rows() {
        let base = record_from_smiles("CC(F)Cl").unwrap();
        let existing_coordinates = (0..base.topology.atoms.len())
            .map(|index| [index as f64, index as f64 + 0.25, index as f64 + 0.5])
            .collect::<Vec<_>>();
        let mut source = base;
        fixture_add_conformer(
            &mut source,
            Conformer3D::new(0, existing_coordinates.clone(), true),
        )
        .unwrap();
        let source_before = source.clone();
        let mut workspace = EnumerationWorkspace::prepare(
            &source,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: false,
            },
            project_query,
        )
        .unwrap();
        assert_eq!(workspace.center_count(), 1);

        let configuration = ConfigurationBits::from_biguint(BigUint::from(1_u8));
        workspace.apply_configuration(&configuration).unwrap();
        let finalized = workspace
            .finalize_configuration(false, &mut HashSet::new())
            .unwrap()
            .unwrap();
        let embedded_a = workspace
            .apply_embedding_filter(finalized.clone(), &configuration)
            .unwrap()
            .unwrap();
        let embedded_b = workspace
            .apply_embedding_filter(finalized, &configuration)
            .unwrap()
            .unwrap();

        assert_eq!(embedded_a.topology.atoms.len(), source.topology.atoms.len());
        assert_eq!(embedded_a.coordinates.conformers_3d.len(), 2);
        assert_eq!(
            embedded_a.coordinates.conformers_3d[0].coordinates(),
            existing_coordinates
        );
        assert_eq!(embedded_a.coordinates.conformers_3d[0].id(), 0);
        assert_eq!(embedded_a.coordinates.conformers_3d[1].id(), 0);
        assert_eq!(
            embedded_a.coordinates.conformers_3d[1].coordinates().len(),
            source.topology.atoms.len()
        );
        assert!(
            embedded_a.coordinates.conformers_3d[1]
                .coordinates()
                .iter()
                .flatten()
                .all(|coordinate| coordinate.is_finite())
        );
        assert_eq!(
            embedded_a.coordinates.conformers_3d[1].coordinates(),
            embedded_b.coordinates.conformers_3d[1].coordinates(),
            "the configuration-derived explicit seed must be deterministic"
        );
        assert_eq!(source, source_before);
    }

    #[test]
    fn python_flippers_repeat_exact_atom_bond_and_group_setter_directions() {
        let mut builder = empty_record();
        let atom = fixture_add_atom(
            &mut builder,
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
        );
        let other = fixture_add_atom(
            &mut builder,
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw),
        );
        let untouched = fixture_add_atom(
            &mut builder,
            AtomSpec::new(Element::F).with_chiral_tag(ChiralTag::Other),
        );
        let end_controller = fixture_add_atom(&mut builder, AtomSpec::new(Element::CL));
        let bond = fixture_add_bond(
            &mut builder,
            BondSpec::new(atom, other, BondOrder::Double)
                .with_stereo_atoms(untouched, end_controller),
        )
        .unwrap();
        fixture_add_bond(
            &mut builder,
            BondSpec::new(atom, untouched, BondOrder::Single),
        )
        .unwrap();
        fixture_add_bond(
            &mut builder,
            BondSpec::new(other, end_controller, BondOrder::Single),
        )
        .unwrap();
        let mut molecule = checked_record(builder);

        let atom_flipper = StereoFlipper::Atom { atom };
        atom_flipper.flip(&mut molecule, false).unwrap();
        assert_eq!(
            molecule.topology.atoms[atom.index()].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        atom_flipper.flip(&mut molecule, true).unwrap();
        atom_flipper.flip(&mut molecule, true).unwrap();
        assert_eq!(
            molecule.topology.atoms[atom.index()].chiral_tag(),
            ChiralTag::TetrahedralCw
        );

        let bond_flipper = StereoFlipper::Bond { bond };
        bond_flipper.flip(&mut molecule, false).unwrap();
        assert_eq!(
            molecule.topology.bonds[bond.index()].stereo(),
            BondStereo::Trans
        );
        bond_flipper.flip(&mut molecule, true).unwrap();
        bond_flipper.flip(&mut molecule, true).unwrap();
        assert_eq!(
            molecule.topology.bonds[bond.index()].stereo(),
            BondStereo::Cis
        );

        let group_flipper = StereoFlipper::StereoGroup {
            original_parities: vec![
                (atom, ChiralTag::TetrahedralCw),
                (other, ChiralTag::TetrahedralCcw),
                (untouched, ChiralTag::Other),
            ],
        };
        group_flipper.flip(&mut molecule, false).unwrap();
        assert_eq!(
            molecule.topology.atoms[atom.index()].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(
            molecule.topology.atoms[other.index()].chiral_tag(),
            ChiralTag::TetrahedralCw
        );
        assert_eq!(
            molecule.topology.atoms[untouched.index()].chiral_tag(),
            ChiralTag::Other
        );
        group_flipper.flip(&mut molecule, true).unwrap();
        group_flipper.flip(&mut molecule, true).unwrap();
        assert_eq!(
            molecule.topology.atoms[atom.index()].chiral_tag(),
            ChiralTag::TetrahedralCw
        );
        assert_eq!(
            molecule.topology.atoms[other.index()].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(
            molecule.topology.atoms[untouched.index()].chiral_tag(),
            ChiralTag::Other
        );
    }

    #[test]
    fn python_flipper_selection_filters_records_and_groups_in_source_order() {
        let (mut molecule, assigned, unassigned, begin_controller, end_controller, double_bond) =
            grouped_selection_molecule();
        let records = vec![
            stereo_info(
                PotentialStereoType::AtomTetrahedral,
                PotentialStereoSpecified::Specified,
                PotentialStereoCenter::Atom(assigned),
                Vec::new(),
            ),
            stereo_info(
                PotentialStereoType::AtomTetrahedral,
                PotentialStereoSpecified::Unspecified,
                PotentialStereoCenter::Atom(unassigned),
                Vec::new(),
            ),
            stereo_info(
                PotentialStereoType::BondDouble,
                PotentialStereoSpecified::Unknown,
                PotentialStereoCenter::Bond(double_bond),
                vec![Some(begin_controller), None, Some(end_controller), None],
            ),
        ];

        let all = select_stereo_flippers_from_info(
            &mut molecule.clone(),
            &records,
            FlipperSelectionOptions {
                only_unassigned: false,
                only_stereo_groups: false,
            },
        )
        .unwrap();
        assert_eq!(
            all,
            vec![
                StereoFlipper::Atom { atom: assigned },
                StereoFlipper::Atom { atom: unassigned },
                StereoFlipper::Bond { bond: double_bond },
            ]
        );

        let unassigned_only = select_stereo_flippers_from_info(
            &mut molecule,
            &records,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: false,
            },
        )
        .unwrap();
        assert_eq!(unassigned_only.len(), 4);
        assert_eq!(
            &unassigned_only[..2],
            &[
                StereoFlipper::Atom { atom: unassigned },
                StereoFlipper::Bond { bond: double_bond },
            ]
        );
        assert!(matches!(
            &unassigned_only[2],
            StereoFlipper::StereoGroup { original_parities }
                if original_parities == &vec![
                    (assigned, ChiralTag::TetrahedralCw),
                    (unassigned, ChiralTag::Unspecified),
                ]
        ));
        assert!(matches!(
            &unassigned_only[3],
            StereoFlipper::StereoGroup { original_parities }
                if original_parities == &vec![(begin_controller, ChiralTag::TetrahedralCcw)]
        ));

        let groups_only = select_stereo_flippers_from_info(
            &mut molecule,
            &records,
            FlipperSelectionOptions {
                only_unassigned: true,
                only_stereo_groups: true,
            },
        )
        .unwrap();
        assert_eq!(groups_only, unassigned_only[2..]);

        let disabled_groups_only = select_stereo_flippers_from_info(
            &mut molecule,
            &records,
            FlipperSelectionOptions {
                only_unassigned: false,
                only_stereo_groups: true,
            },
        )
        .unwrap();
        assert!(disabled_groups_only.is_empty());
    }

    #[test]
    fn python_bond_flipper_selection_initializes_preserves_and_rejects_controllers_exactly() {
        let (molecule, _, _, begin_controller, end_controller, double_bond) =
            grouped_selection_molecule();
        let complete = stereo_info(
            PotentialStereoType::BondDouble,
            PotentialStereoSpecified::Unspecified,
            PotentialStereoCenter::Bond(double_bond),
            vec![Some(begin_controller), None, Some(end_controller), None],
        );
        let options = FlipperSelectionOptions {
            only_unassigned: false,
            only_stereo_groups: false,
        };

        let mut initialized = molecule.clone();
        let flippers = select_stereo_flippers_from_info(
            &mut initialized,
            std::slice::from_ref(&complete),
            options,
        )
        .unwrap();
        assert_eq!(flippers, vec![StereoFlipper::Bond { bond: double_bond }]);
        assert_eq!(
            initialized.topology.bonds[double_bond.index()].stereo_atoms(),
            Some([begin_controller, end_controller])
        );

        let mut preinitialized = molecule.clone();
        preinitialized.topology.bonds[double_bond.index()]
            .set_stereo_atoms(Some([end_controller, begin_controller]));
        select_stereo_flippers_from_info(
            &mut preinitialized,
            std::slice::from_ref(&complete),
            options,
        )
        .unwrap();
        assert_eq!(
            preinitialized.topology.bonds[double_bond.index()].stereo_atoms(),
            Some([end_controller, begin_controller])
        );

        for missing in [0, 2] {
            let mut controllers = complete.controlling_atoms.to_vec();
            controllers[missing] = None;
            let incomplete = stereo_info(
                PotentialStereoType::BondDouble,
                PotentialStereoSpecified::Unspecified,
                PotentialStereoCenter::Bond(double_bond),
                controllers,
            );
            let mut candidate = molecule.clone();
            let flippers = select_stereo_flippers_from_info(
                &mut candidate,
                std::slice::from_ref(&incomplete),
                options,
            )
            .unwrap();
            assert!(flippers.is_empty());
            assert_eq!(
                candidate.topology.bonds[double_bond.index()].stereo_atoms(),
                None
            );
        }
    }

    #[test]
    fn python_flipper_selection_preserves_represented_atropisomers() {
        let mut builder = empty_record();
        let begin = fixture_add_atom(&mut builder, AtomSpec::new(Element::C));
        let end = fixture_add_atom(&mut builder, AtomSpec::new(Element::C));
        let begin_controller = fixture_add_atom(&mut builder, AtomSpec::new(Element::F));
        let end_controller = fixture_add_atom(&mut builder, AtomSpec::new(Element::CL));
        let bond = fixture_add_bond(
            &mut builder,
            BondSpec::new(begin, end, BondOrder::Single).with_stereo(BondStereo::AtropCw),
        )
        .unwrap();
        fixture_add_bond(
            &mut builder,
            BondSpec::new(begin, begin_controller, BondOrder::Single),
        )
        .unwrap();
        fixture_add_bond(
            &mut builder,
            BondSpec::new(end, end_controller, BondOrder::Single),
        )
        .unwrap();
        let mut molecule = checked_record(builder);
        let info = stereo_info(
            PotentialStereoType::BondAtropisomer,
            PotentialStereoSpecified::Specified,
            PotentialStereoCenter::Bond(bond),
            vec![Some(begin_controller), None, Some(end_controller), None],
        );

        let flippers = select_stereo_flippers_from_info(
            &mut molecule,
            &[info],
            FlipperSelectionOptions {
                only_unassigned: false,
                only_stereo_groups: false,
            },
        )
        .unwrap();
        assert!(flippers.is_empty());
        assert_eq!(
            molecule.topology.bonds[bond.index()].stereo(),
            BondStereo::AtropCw
        );
        assert_eq!(molecule.topology.bonds[bond.index()].stereo_atoms(), None);
    }

    #[test]
    fn source_parameter_identity_and_all_six_fields_do_not_draw_random_bits() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let calls = Arc::new(AtomicUsize::new(0));
        let counter = calls.clone();
        let provider = StereoisomerRandomSource::from_random_bits(Box::new(move |_| {
            counter.fetch_add(1, Ordering::SeqCst);
            Ok(BigUint::from(0_u8))
        }));
        let other =
            StereoisomerRandomSource::from_random_bits(Box::new(|_| Ok(BigUint::from(0_u8))));
        assert_eq!(provider, provider.clone());
        assert_ne!(provider, other);
        assert_eq!(format!("{provider:?}"), "RandomBits(..)");
        let mut options = StereoisomerOptions::default();
        assert!(!options.try_embedding());
        assert!(options.only_unassigned());
        assert!(!options.only_stereo_groups());
        assert_eq!(options.max_isomers(), 1024);
        assert!(options.random_source().is_none());
        assert!(options.unique());
        options.set_try_embedding(true);
        options.set_only_unassigned(false);
        options.set_only_stereo_groups(true);
        options.set_max_isomers(3);
        options.set_random_source(Some(provider.clone()));
        options.set_unique(false);
        let constructed =
            StereoisomerOptions::new(true, false, 3, Some(provider.clone()), false, true);
        assert_eq!(options, constructed);
        assert_eq!(options.random_source(), Some(&provider));
        assert_eq!(calls.load(Ordering::SeqCst), 0);
    }

    #[test]
    fn shared_random_parameter_is_lazy_and_shares_state_across_enumerations() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let calls = Arc::new(AtomicUsize::new(0));
        let counter = calls.clone();
        let source = record_from_smiles("FC(Cl)C(Cl)F").unwrap();
        let options = StereoisomerOptions::new(
            false,
            true,
            1,
            Some(StereoisomerRandomSource::from_random_bits(Box::new(
                move |width| {
                    assert_eq!(width, 2);
                    Ok(BigUint::from(counter.fetch_add(1, Ordering::SeqCst)))
                },
            ))),
            true,
            false,
        );
        let mut first = enumerate_stereoisomers(&source, options.clone(), project_query).unwrap();
        let mut second = enumerate_stereoisomers(&source, options, project_query).unwrap();
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        assert_eq!(
            text(cosmolkit_smiles::write_smiles(&first.next().unwrap().unwrap()).unwrap()),
            "F[C@H](Cl)[C@@H](F)Cl"
        );
        assert_eq!(
            text(cosmolkit_smiles::write_smiles(&second.next().unwrap().unwrap()).unwrap()),
            "F[C@@H](Cl)[C@@H](F)Cl"
        );
        assert!(first.next().is_none());
        assert!(second.next().is_none());
        assert_eq!(calls.load(Ordering::SeqCst), 2);
        // No-center and exhaustive branches never ask the provider for bits.
        let options = StereoisomerOptions::new(
            false,
            true,
            0,
            Some(StereoisomerRandomSource::from_random_bits(Box::new(|_| {
                panic!("unexpected draw")
            }))),
            true,
            false,
        );
        assert_eq!(
            enumerate_stereoisomers(&source, options.clone(), project_query)
                .unwrap()
                .count(),
            3
        );
        let no_center = record_from_smiles("CC").unwrap();
        assert_eq!(
            enumerate_stereoisomers(&no_center, options, project_query)
                .unwrap()
                .count(),
            1
        );
    }

    #[test]
    fn cloned_integer_parameter_reinitializes_each_enumeration() {
        let source = record_from_smiles("FC(Cl)C(Cl)F").unwrap();
        let options = StereoisomerOptions::new(
            false,
            true,
            2,
            Some(StereoisomerRandomSource::from_integer_seed(BigInt::from(
                61_453,
            ))),
            false,
            false,
        );
        let collect = |options| {
            enumerate_stereoisomers(&source, options, project_query)
                .unwrap()
                .map(|row| text(cosmolkit_smiles::write_smiles(&row.unwrap()).unwrap()))
                .collect::<Vec<_>>()
        };
        let first = collect(options.clone());
        assert_eq!(first.len(), 2);
        assert_eq!(first, collect(options));
    }

    #[test]
    fn shared_random_parameter_callback_failure_is_deferred_and_fused() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let calls = Arc::new(AtomicUsize::new(0));
        let counter = calls.clone();
        let source = record_from_smiles("FC(Cl)C(Cl)F").unwrap();
        let before = source.clone();
        let options = StereoisomerOptions::new(
            false,
            true,
            3,
            Some(StereoisomerRandomSource::from_random_bits(Box::new(
                move |width| {
                    assert_eq!(width, 2);
                    if counter.fetch_add(1, Ordering::SeqCst) == 0 {
                        Ok(BigUint::from(0_u8))
                    } else {
                        Err("provider fixture failure".into())
                    }
                },
            ))),
            true,
            false,
        );
        let mut iterator = enumerate_stereoisomers(&source, options, project_query).unwrap();
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        assert!(iterator.next().unwrap().is_ok());
        assert!(
            matches!(iterator.next().unwrap(), Err(EnumerationError::RandomBitsSource { bit_count: 2, message }) if message == "provider fixture failure")
        );
        assert_eq!(iterator.yielded_count(), 1);
        assert!(iterator.next().is_none());
        assert_eq!(calls.load(Ordering::SeqCst), 2);
        assert_eq!(source, before);
    }

    #[test]
    fn shared_random_parameter_poison_keeps_a_typed_error_cause() {
        use std::error::Error;
        let source = record_from_smiles("FC(Cl)C(Cl)F").unwrap();
        let options = StereoisomerOptions::new(
            false,
            true,
            1,
            Some(StereoisomerRandomSource::from_random_bits(Box::new(|_| {
                panic!("provider panic fixture")
            }))),
            true,
            false,
        );
        let mut first = enumerate_stereoisomers(&source, options.clone(), project_query).unwrap();
        assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| first.next())).is_err());
        let mut second = enumerate_stereoisomers(&source, options, project_query).unwrap();
        let error = second.next().unwrap().unwrap_err();
        assert!(matches!(error, EnumerationError::RandomSourcePoisoned(_)));
        assert!(
            error
                .source()
                .unwrap()
                .downcast_ref::<std::sync::PoisonError<()>>()
                .is_some()
        );
        assert_eq!(second.yielded_count(), 0);
        assert!(second.next().is_none());
    }
}
