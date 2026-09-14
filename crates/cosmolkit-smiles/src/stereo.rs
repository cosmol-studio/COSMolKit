use cosmolkit_model::BondId;
use cosmolkit_types::{BondOrder, ChiralTag};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub(crate) enum SmilesStereoOrderError {
    #[error("bond index {bond_index} is out of bounds for bond count {bond_count}")]
    BondIndexOutOfRange {
        bond_index: usize,
        bond_count: usize,
    },
    #[error(
        "probe vector length {probe_count} does not match incident-bond length {incident_count}"
    )]
    ProbeLengthMismatch {
        incident_count: usize,
        probe_count: usize,
    },
    #[error("probe entry at position {probe_index} is missing from the incident-bond order")]
    ProbeMemberMissing { probe_index: usize },
}

pub(crate) fn invert_tetrahedral_tag(tag: ChiralTag) -> ChiralTag {
    match tag {
        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
        ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
        other => other,
    }
}

pub(crate) fn count_swaps_to_interconvert<T: Copy + Eq>(
    reference: &[T],
    mut probe: Vec<T>,
) -> Option<usize> {
    // BEGIN RDKIT CPP FUNCTION RDGeneral::countSwapsToInterconvert
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: unsigned int countSwapsToInterconvert(const T &ref, T probe) {
    // RDKit✔️✔️:   PRECONDITION(ref.size() == probe.size(), "size mismatch");
    // RDKit✔️✔️:   typename T::const_iterator refIt = ref.begin();
    // RDKit✔️✔️:   typename T::iterator probeIt = probe.begin();
    // RDKit✔️✔️:   typename T::iterator probeIt2;
    // RDKit✔️✔️:   unsigned int nSwaps = 0;
    // RDKit✔️✔️:   while (refIt != ref.end()) {
    // RDKit✔️✔️:     if ((*probeIt) != (*refIt)) {
    // RDKit✔️✔️:       bool foundIt = false;
    // RDKit✔️✔️:       probeIt2 = probeIt;
    // RDKit✔️✔️:       while ((*probeIt2) != (*refIt) && probeIt2 != probe.end()) {
    // RDKit✔️✔️:         ++probeIt2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (probeIt2 != probe.end()) {
    // RDKit✔️✔️:         foundIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       CHECK_INVARIANT(foundIt, "could not find probe element");
    // RDKit✔️✔️:       std::swap(*probeIt, *probeIt2);
    // RDKit✔️✔️:       nSwaps++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++probeIt;
    // RDKit✔️✔️:     ++refIt;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nSwaps;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDGeneral::countSwapsToInterconvert
    if reference.len() != probe.len() {
        return None;
    }
    let mut swaps = 0;
    for (index, expected) in reference.iter().enumerate() {
        if probe[index] == *expected {
            continue;
        }
        let found = probe[index..]
            .iter()
            .position(|candidate| candidate == expected)
            .map(|offset| index + offset)?;
        probe.swap(index, found);
        swaps += 1;
    }
    Some(swaps)
}

pub(crate) fn atom_has_fourth_valence(
    explicit_hydrogens: u8,
    implicit_valence_is_one: bool,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION Canon::details::atomHasFourthValence modeled molecule branch
    // RDKit✔️✔️: bool atomHasFourthValence(const Atom *atom) {
    // RDKit✔️✔️:   if (atom->getNumExplicitHs() == 1 ||
    // RDKit✔️✔️:       (!atom->needsUpdatePropertyCache() &&
    // RDKit✔️✔️:        atom->getValence(Atom::ValenceType::IMPLICIT) == 1)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::details::atomHasFourthValence modeled molecule branch
    explicit_hydrogens == 1 || implicit_valence_is_one
}

pub(crate) fn chiral_atom_needs_tag_inversion(
    degree: usize,
    explicit_hydrogens: u8,
    is_atom_first: bool,
    has_fourth_valence: bool,
    num_closures: usize,
    is_unsaturated: bool,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION Canon::chiralAtomNeedsTagInversion
    // RDKit✔️✔️: bool chiralAtomNeedsTagInversion(const RDKit::ROMol &mol,
    // RDKit✔️✔️:                                  const RDKit::Atom *atom, bool isAtomFirst,
    // RDKit✔️✔️:                                  size_t numClosures) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   return atom->getDegree() == 3 &&
    // RDKit✔️✔️:          ((isAtomFirst && atom->getNumExplicitHs() == 1) ||
    // RDKit✔️✔️:           (!details::atomHasFourthValence(atom) && numClosures == 1 &&
    // RDKit✔️✔️:            !details::isUnsaturated(atom, mol)));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::chiralAtomNeedsTagInversion
    degree == 3
        && ((is_atom_first && explicit_hydrogens == 1)
            || (!has_fourth_valence && num_closures == 1 && !is_unsaturated))
}

pub(crate) fn bond_order_as_double(order: BondOrder) -> f64 {
    match order {
        BondOrder::Single
        | BondOrder::Dative
        | BondOrder::DativeOne
        | BondOrder::DativeLeft
        | BondOrder::DativeRight => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Quintuple => 5.0,
        BondOrder::Hextuple => 6.0,
        BondOrder::Aromatic | BondOrder::OneAndHalf => 1.5,
        BondOrder::TwoAndHalf => 2.5,
        BondOrder::ThreeAndHalf => 3.5,
        BondOrder::FourAndHalf => 4.5,
        BondOrder::FiveAndHalf => 5.5,
        BondOrder::Ionic
        | BondOrder::Hydrogen
        | BondOrder::ThreeCenter
        | BondOrder::Other
        | BondOrder::Zero
        | BondOrder::Unspecified => 0.0,
    }
}

// RDKit✔️✔️: constexpr unsigned char swap_squareplanar_table[4][6] = {
// RDKit✔️✔️:     {0, 0, 0, 0, 0, 0},
// RDKit✔️✔️:     {3, 1, 2, 2, 1, 3},  // SP1
// RDKit✔️✔️:     {2, 3, 1, 1, 3, 2},  // SP2
// RDKit✔️✔️:     {1, 2, 3, 3, 2, 1}   // SP3
// RDKit✔️✔️: };
const SWAP_SQUAREPLANAR_TABLE: [[u8; 6]; 4] = [
    [0, 0, 0, 0, 0, 0],
    [3, 1, 2, 2, 1, 3],
    [2, 3, 1, 1, 3, 2],
    [1, 2, 3, 3, 2, 1],
];

// RDKit✔️✔️: constexpr unsigned char swap_trigonalbipyramidal_table[21][10] = {
// RDKit✔️✔️:     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
// RDKit✔️✔️:     {9, 20, 17, 2, 2, 2, 7, 2, 6, 3},        // TB1
// RDKit✔️✔️:     {11, 15, 18, 1, 1, 1, 8, 1, 5, 4},       // TB2
// RDKit✔️✔️:     {10, 19, 4, 18, 4, 8, 4, 5, 4, 1},       // TB3
// RDKit✔️✔️:     {12, 16, 3, 17, 3, 7, 3, 6, 3, 2},       // TB4
// RDKit✔️✔️:     {13, 6, 16, 20, 7, 6, 6, 3, 2, 6},       // TB5
// RDKit✔️✔️:     {14, 5, 19, 15, 8, 5, 5, 4, 1, 5},       // TB6
// RDKit✔️✔️:     {8, 14, 10, 11, 5, 4, 1, 8, 8, 8},       // TB7
// RDKit✔️✔️:     {7, 13, 12, 9, 6, 3, 2, 7, 7, 7},        // TB8
// RDKit✔️✔️:     {1, 11, 11, 8, 15, 18, 11, 11, 14, 10},  // TB9
// RDKit✔️✔️:     {3, 12, 7, 12, 16, 12, 17, 13, 12, 9},   // TB10
// RDKit✔️✔️:     {2, 9, 9, 7, 20, 17, 9, 9, 13, 12},      // TB11
// RDKit✔️✔️:     {4, 10, 8, 10, 19, 10, 18, 14, 10, 11},  // TB12
// RDKit✔️✔️:     {5, 8, 14, 14, 14, 19, 15, 10, 11, 14},  // TB13
// RDKit✔️✔️:     {6, 7, 13, 13, 13, 16, 20, 12, 9, 13},   // TB14
// RDKit✔️✔️:     {20, 2, 20, 6, 9, 20, 13, 17, 20, 16},   // TB15
// RDKit✔️✔️:     {19, 4, 5, 19, 10, 14, 19, 19, 18, 15},  // TB16
// RDKit✔️✔️:     {18, 18, 1, 4, 18, 11, 10, 15, 19, 18},  // TB17
// RDKit✔️✔️:     {17, 17, 2, 3, 17, 9, 12, 20, 16, 17},   // TB18
// RDKit✔️✔️:     {16, 3, 6, 16, 12, 13, 16, 16, 17, 20},  // TB19
// RDKit✔️✔️:     {15, 1, 15, 5, 11, 15, 14, 18, 15, 19}   // TB20
// RDKit✔️✔️: };
const SWAP_TRIGONALBIPYRAMIDAL_TABLE: [[u8; 10]; 21] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 20, 17, 2, 2, 2, 7, 2, 6, 3],
    [11, 15, 18, 1, 1, 1, 8, 1, 5, 4],
    [10, 19, 4, 18, 4, 8, 4, 5, 4, 1],
    [12, 16, 3, 17, 3, 7, 3, 6, 3, 2],
    [13, 6, 16, 20, 7, 6, 6, 3, 2, 6],
    [14, 5, 19, 15, 8, 5, 5, 4, 1, 5],
    [8, 14, 10, 11, 5, 4, 1, 8, 8, 8],
    [7, 13, 12, 9, 6, 3, 2, 7, 7, 7],
    [1, 11, 11, 8, 15, 18, 11, 11, 14, 10],
    [3, 12, 7, 12, 16, 12, 17, 13, 12, 9],
    [2, 9, 9, 7, 20, 17, 9, 9, 13, 12],
    [4, 10, 8, 10, 19, 10, 18, 14, 10, 11],
    [5, 8, 14, 14, 14, 19, 15, 10, 11, 14],
    [6, 7, 13, 13, 13, 16, 20, 12, 9, 13],
    [20, 2, 20, 6, 9, 20, 13, 17, 20, 16],
    [19, 4, 5, 19, 10, 14, 19, 19, 18, 15],
    [18, 18, 1, 4, 18, 11, 10, 15, 19, 18],
    [17, 17, 2, 3, 17, 9, 12, 20, 16, 17],
    [16, 3, 6, 16, 12, 13, 16, 16, 17, 20],
    [15, 1, 15, 5, 11, 15, 14, 18, 15, 19],
];

// RDKit✔️✔️: constexpr unsigned char swap_octahedral_table[31][15] = {
// RDKit✔️✔️:     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
// RDKit✔️✔️:     {17, 16, 30, 21, 2, 14, 2, 10, 25, 8, 2, 22, 4, 7, 3},       // OH1
// RDKit✔️✔️:     {7, 3, 25, 22, 1, 4, 1, 8, 30, 10, 1, 21, 14, 17, 16},       // OH2
// RDKit✔️✔️:     {18, 2, 29, 16, 22, 15, 16, 26, 11, 9, 21, 16, 6, 5, 1},     // OH3
// RDKit✔️✔️:     {15, 18, 19, 28, 14, 2, 8, 14, 27, 14, 10, 24, 1, 6, 5},     // OH4
// RDKit✔️✔️:     {14, 17, 20, 15, 27, 16, 9, 28, 15, 15, 23, 11, 7, 3, 4},    // OH5
// RDKit✔️✔️:     {16, 14, 18, 26, 24, 17, 29, 18, 13, 19, 12, 18, 3, 4, 7},   // OH6
// RDKit✔️✔️:     {2, 15, 17, 23, 25, 18, 30, 12, 17, 20, 17, 13, 5, 1, 6},    // OH7
// RDKit✔️✔️:     {23, 26, 11, 12, 10, 10, 4, 2, 29, 1, 14, 20, 10, 13, 9},    // OH8
// RDKit✔️✔️:     {24, 25, 10, 11, 13, 11, 5, 30, 16, 3, 19, 15, 12, 11, 8},   // OH9
// RDKit✔️✔️:     {20, 29, 9, 13, 8, 8, 14, 1, 26, 2, 4, 23, 8, 12, 11},       // OH10
// RDKit✔️✔️:     {19, 30, 8, 9, 12, 9, 15, 25, 3, 16, 24, 5, 13, 9, 10},      // OH11
// RDKit✔️✔️:     {22, 27, 13, 8, 11, 13, 28, 7, 18, 21, 6, 17, 9, 10, 13},    // OH12
// RDKit✔️✔️:     {21, 28, 12, 10, 9, 12, 27, 17, 6, 22, 18, 7, 11, 8, 12},    // OH13
// RDKit✔️✔️:     {5, 6, 24, 27, 4, 1, 10, 4, 28, 4, 8, 19, 2, 18, 15},        // OH14
// RDKit✔️✔️:     {4, 7, 23, 5, 28, 3, 11, 27, 5, 5, 20, 9, 17, 16, 14},       // OH15
// RDKit✔️✔️:     {6, 1, 26, 3, 21, 5, 3, 29, 9, 11, 22, 3, 18, 15, 2},        // OH16
// RDKit✔️✔️:     {1, 5, 7, 20, 30, 6, 25, 13, 7, 23, 7, 12, 15, 2, 18},       // OH17
// RDKit✔️✔️:     {3, 4, 6, 29, 19, 7, 26, 6, 12, 24, 13, 6, 16, 14, 17},      // OH18
// RDKit✔️✔️:     {11, 24, 4, 30, 18, 25, 23, 24, 22, 6, 9, 14, 21, 24, 20},   // OH19
// RDKit✔️✔️:     {10, 23, 5, 17, 29, 26, 24, 21, 23, 7, 15, 8, 23, 22, 19},   // OH20
// RDKit✔️✔️:     {13, 22, 28, 1, 16, 27, 22, 20, 24, 12, 3, 2, 19, 23, 22},   // OH21
// RDKit✔️✔️:     {12, 21, 27, 2, 3, 28, 21, 23, 19, 13, 16, 1, 24, 20, 21},   // OH22
// RDKit✔️✔️:     {8, 20, 15, 7, 26, 29, 19, 22, 20, 17, 5, 10, 20, 21, 24},   // OH23
// RDKit✔️✔️:     {9, 19, 14, 25, 6, 30, 20, 19, 21, 18, 11, 4, 22, 19, 23},   // OH24
// RDKit✔️✔️:     {30, 9, 2, 24, 7, 19, 17, 11, 1, 29, 30, 28, 27, 30, 26},    // OH25
// RDKit✔️✔️:     {29, 8, 16, 6, 23, 20, 18, 3, 10, 30, 27, 29, 29, 28, 25},   // OH26
// RDKit✔️✔️:     {28, 12, 22, 14, 5, 21, 13, 15, 4, 28, 26, 30, 25, 29, 28},  // OH27
// RDKit✔️✔️:     {27, 13, 21, 4, 15, 22, 12, 5, 14, 27, 29, 25, 30, 26, 27},  // OH28
// RDKit✔️✔️:     {26, 10, 3, 18, 20, 23, 6, 16, 8, 25, 28, 26, 26, 27, 30},   // OH29
// RDKit✔️✔️:     {25, 11, 1, 19, 17, 24, 7, 9, 2, 26, 25, 27, 28, 25, 29}     // OH30
// RDKit✔️✔️: };
const SWAP_OCTAHEDRAL_TABLE: [[u8; 15]; 31] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [17, 16, 30, 21, 2, 14, 2, 10, 25, 8, 2, 22, 4, 7, 3],
    [7, 3, 25, 22, 1, 4, 1, 8, 30, 10, 1, 21, 14, 17, 16],
    [18, 2, 29, 16, 22, 15, 16, 26, 11, 9, 21, 16, 6, 5, 1],
    [15, 18, 19, 28, 14, 2, 8, 14, 27, 14, 10, 24, 1, 6, 5],
    [14, 17, 20, 15, 27, 16, 9, 28, 15, 15, 23, 11, 7, 3, 4],
    [16, 14, 18, 26, 24, 17, 29, 18, 13, 19, 12, 18, 3, 4, 7],
    [2, 15, 17, 23, 25, 18, 30, 12, 17, 20, 17, 13, 5, 1, 6],
    [23, 26, 11, 12, 10, 10, 4, 2, 29, 1, 14, 20, 10, 13, 9],
    [24, 25, 10, 11, 13, 11, 5, 30, 16, 3, 19, 15, 12, 11, 8],
    [20, 29, 9, 13, 8, 8, 14, 1, 26, 2, 4, 23, 8, 12, 11],
    [19, 30, 8, 9, 12, 9, 15, 25, 3, 16, 24, 5, 13, 9, 10],
    [22, 27, 13, 8, 11, 13, 28, 7, 18, 21, 6, 17, 9, 10, 13],
    [21, 28, 12, 10, 9, 12, 27, 17, 6, 22, 18, 7, 11, 8, 12],
    [5, 6, 24, 27, 4, 1, 10, 4, 28, 4, 8, 19, 2, 18, 15],
    [4, 7, 23, 5, 28, 3, 11, 27, 5, 5, 20, 9, 17, 16, 14],
    [6, 1, 26, 3, 21, 5, 3, 29, 9, 11, 22, 3, 18, 15, 2],
    [1, 5, 7, 20, 30, 6, 25, 13, 7, 23, 7, 12, 15, 2, 18],
    [3, 4, 6, 29, 19, 7, 26, 6, 12, 24, 13, 6, 16, 14, 17],
    [11, 24, 4, 30, 18, 25, 23, 24, 22, 6, 9, 14, 21, 24, 20],
    [10, 23, 5, 17, 29, 26, 24, 21, 23, 7, 15, 8, 23, 22, 19],
    [13, 22, 28, 1, 16, 27, 22, 20, 24, 12, 3, 2, 19, 23, 22],
    [12, 21, 27, 2, 3, 28, 21, 23, 19, 13, 16, 1, 24, 20, 21],
    [8, 20, 15, 7, 26, 29, 19, 22, 20, 17, 5, 10, 20, 21, 24],
    [9, 19, 14, 25, 6, 30, 20, 19, 21, 18, 11, 4, 22, 19, 23],
    [30, 9, 2, 24, 7, 19, 17, 11, 1, 29, 30, 28, 27, 30, 26],
    [29, 8, 16, 6, 23, 20, 18, 3, 10, 30, 27, 29, 29, 28, 25],
    [28, 12, 22, 14, 5, 21, 13, 15, 4, 28, 26, 30, 25, 29, 28],
    [27, 13, 21, 4, 15, 22, 12, 5, 14, 27, 29, 25, 30, 26, 27],
    [26, 10, 3, 18, 20, 23, 6, 16, 8, 25, 28, 26, 26, 27, 30],
    [25, 11, 1, 19, 17, 24, 7, 9, 2, 26, 25, 27, 28, 25, 29],
];

fn swap_squareplanar(perm: u32, x: usize, y: usize) -> u32 {
    // RDKit✔️✔️: static unsigned int swap_squareplanar(unsigned int perm, unsigned int x,
    // RDKit✔️✔️:                                       unsigned int y) {
    // RDKit✔️✔️:   constexpr unsigned int offset[3] = {0, 2, 3};
    // RDKit✔️✔️:   unsigned int swapidx;
    // RDKit✔️✔️:   if (x == y) {
    // RDKit✔️✔️:     return perm;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (x < y) {
    // RDKit✔️✔️:     if (y > 3) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     swapidx = offset[x] + (y - 1);
    // RDKit✔️✔️:   } else /* x > y */ {
    // RDKit✔️✔️:     if (x > 3) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     swapidx = offset[y] + (x - 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return perm < 4 ? swap_squareplanar_table[perm][swapidx] : 0;
    // RDKit✔️✔️: }
    swap_nontetrahedral(perm, x, y, &[0, 2, 3], 3, 4, |perm, swap_index| {
        SWAP_SQUAREPLANAR_TABLE[perm][swap_index]
    })
}

fn swap_trigonalbipyramidal(perm: u32, x: usize, y: usize) -> u32 {
    // RDKit✔️✔️: static unsigned int swap_trigonalbipyramidal(unsigned int perm, unsigned int x,
    // RDKit✔️✔️:                                              unsigned int y) {
    // RDKit✔️✔️:   constexpr unsigned int offset[4] = {0, 3, 5, 6};
    // RDKit✔️✔️:   unsigned int swapidx;
    // RDKit✔️✔️:   if (x == y) {
    // RDKit✔️✔️:     return perm;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (x < y) {
    // RDKit✔️✔️:     if (y > 4) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     swapidx = offset[x] + (y - 1);
    // RDKit✔️✔️:   } else /* x > y */ {
    // RDKit✔️✔️:     if (x > 4) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     swapidx = offset[y] + (x - 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return perm < 21 ? swap_trigonalbipyramidal_table[perm][swapidx] : 0;
    // RDKit✔️✔️: }
    swap_nontetrahedral(perm, x, y, &[0, 3, 5, 6], 4, 21, |perm, swap_index| {
        SWAP_TRIGONALBIPYRAMIDAL_TABLE[perm][swap_index]
    })
}

fn swap_octahedral(perm: u32, x: usize, y: usize) -> u32 {
    // RDKit✔️✔️: static unsigned int swap_octahedral(unsigned int perm, unsigned int x,
    // RDKit✔️✔️:                                     unsigned int y) {
    // RDKit✔️✔️:   constexpr unsigned int offset[5] = {0, 4, 7, 9, 10};
    // RDKit✔️✔️:   unsigned int swapidx;
    // RDKit✔️✔️:   if (x == y) {
    // RDKit✔️✔️:     return perm;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (x < y) {
    // RDKit✔️✔️:     if (y > 5) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     swapidx = offset[x] + (y - 1);
    // RDKit✔️✔️:   } else /* x > y */ {
    // RDKit✔️✔️:     if (x > 5) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     swapidx = offset[y] + (x - 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return perm < 31 ? swap_octahedral_table[perm][swapidx] : 0;
    // RDKit✔️✔️: }
    swap_nontetrahedral(perm, x, y, &[0, 4, 7, 9, 10], 5, 31, |perm, swap_index| {
        SWAP_OCTAHEDRAL_TABLE[perm][swap_index]
    })
}

fn swap_nontetrahedral(
    perm: u32,
    x: usize,
    y: usize,
    offsets: &[usize],
    max_index: usize,
    permutation_count: usize,
    table_lookup: impl Fn(usize, usize) -> u8,
) -> u32 {
    if x == y {
        return perm;
    }
    let swap_index = if x < y {
        if y > max_index {
            return 0;
        }
        offsets[x] + y - 1
    } else {
        if x > max_index {
            return 0;
        }
        offsets[y] + x - 1
    };
    if perm as usize >= permutation_count {
        return 0;
    }
    table_lookup(perm as usize, swap_index).into()
}

pub(crate) fn nontetrahedral_max_neighbors(chiral_tag: ChiralTag) -> Option<usize> {
    // RDKit✔️✔️: unsigned int getMaxNbors(const Atom::ChiralType tag) {
    // RDKit✔️✔️:   switch (tag) {
    // RDKit✔️✔️:     case Atom::CHI_TETRAHEDRAL_CW:
    // RDKit✔️✔️:     case Atom::CHI_TETRAHEDRAL_CCW:
    // RDKit✔️✔️:     case Atom::CHI_TETRAHEDRAL:  // fall through
    // RDKit✔️✔️:       return 4;
    // RDKit✔️✔️:     case Atom::CHI_ALLENE:
    // RDKit✔️✔️:       return 2;  // not used other than SMI/SMA parsers?
    // RDKit✔️✔️:     case Atom::CHI_SQUAREPLANAR:
    // RDKit✔️✔️:       return 4;
    // RDKit✔️✔️:     case Atom::CHI_TRIGONALBIPYRAMIDAL:
    // RDKit✔️✔️:       return 5;
    // RDKit✔️✔️:     case Atom::CHI_OCTAHEDRAL:
    // RDKit✔️✔️:       return 6;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "Warning: unexpected chiral tag getMaxNbors(): " << tag
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    match chiral_tag {
        ChiralTag::SquarePlanar => Some(4),
        ChiralTag::TrigonalBipyramidal => Some(5),
        ChiralTag::Octahedral => Some(6),
        _ => None,
    }
}

pub(crate) fn insert_implicit_nontetrahedral_neighbors(
    bonds: &mut Vec<Option<BondId>>,
    chiral_tag: ChiralTag,
    first_atom: bool,
) {
    // RDKit✔️✔️: static void insertImplicitNbors(INT_LIST &bonds, const Atom::ChiralType tag,
    // RDKit✔️✔️:                                 const bool firstAtom) {
    // RDKit✔️✔️:   unsigned int ref_max = Chirality::getMaxNbors(tag);
    // RDKit✔️✔️:   if (bonds.size() < ref_max) {
    // RDKit✔️✔️:     if (firstAtom) {
    // RDKit✔️✔️:       bonds.insert(bonds.begin(), ref_max - bonds.size(), -1);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       bonds.insert(++bonds.begin(), ref_max - bonds.size(), -1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let Some(max_neighbors) = nontetrahedral_max_neighbors(chiral_tag) else {
        return;
    };
    if bonds.len() < max_neighbors {
        let missing = max_neighbors - bonds.len();
        let insert_at = if first_atom { 0 } else { 1.min(bonds.len()) };
        bonds.splice(insert_at..insert_at, std::iter::repeat_n(None, missing));
    }
}

pub(crate) fn nontetrahedral_chiral_permutation(
    mut permutation: u32,
    chiral_tag: ChiralTag,
    bond_count: usize,
    incident_bonds: &[BondId],
    probe: &[Option<BondId>],
    inverse: bool,
) -> Result<u32, SmilesStereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::getChiralPermutation
    // RDKit✔️✔️: unsigned int getChiralPermutation(const Atom *cen, const INT_LIST &probe,
    // RDKit✔️✔️:                                   bool inverse) {
    // RDKit✔️✔️:   PRECONDITION(cen, "bad center pointer");
    // RDKit✔️✔️:   PRECONDITION(cen->hasOwningMol(), "no owning mol");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int perm;
    // RDKit✔️✔️:   if (!cen->getPropIfPresent(common_properties::_chiralPermutation, perm) ||
    // RDKit✔️✔️:       perm <= 0) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   decltype(&swap_octahedral) swap_func = nullptr;
    // RDKit✔️✔️:   switch (cen->getChiralTag()) {
    // RDKit✔️✔️:     case Atom::ChiralType::CHI_OCTAHEDRAL:
    // RDKit✔️✔️:       if (probe.size() > 6) {
    // RDKit✔️✔️:         return 0;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       swap_func = swap_octahedral;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
    // RDKit✔️✔️:       if (probe.size() > 5) {
    // RDKit✔️✔️:         return 0;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       swap_func = swap_trigonalbipyramidal;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Atom::ChiralType::CHI_SQUAREPLANAR:
    // RDKit✔️✔️:       if (probe.size() > 4) {
    // RDKit✔️✔️:         return 0;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       swap_func = swap_squareplanar;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!swap_func) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int nbrIdx = 0;
    // RDKit✔️✔️:   std::vector<int> order(cen->getOwningMol().getNumBonds(), -1);
    // RDKit✔️✔️:   for (const auto bnd : cen->getOwningMol().atomBonds(cen)) {
    // RDKit✔️✔️:     order[bnd->getIdx()] = nbrIdx++;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // nbrPerm maps original index to array position
    // RDKit✔️✔️:   std::vector<unsigned int> nbrPerm(nbrIdx);
    // RDKit✔️✔️:   std::iota(nbrPerm.begin(), nbrPerm.end(), 0);
    // RDKit✔️✔️:   std::vector<unsigned int> probePerm(probe.size());
    // RDKit✔️✔️:   nbrIdx = 0;
    // RDKit✔️✔️:   for (auto v : probe) {
    // RDKit✔️✔️:     probePerm[nbrIdx++] = v < 0 ? -1 : order[v];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Missing (implicit) neighbors are at the end when in storage order
    // RDKit✔️✔️:   if (nbrPerm.size() < nbrIdx) {
    // RDKit✔️✔️:     nbrPerm.insert(nbrPerm.end(), nbrIdx - nbrPerm.size(), -1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CHECK_INVARIANT(nbrPerm.size() == probePerm.size(),
    // RDKit✔️✔️:                   "probe vector size does not match");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (inverse) {
    // RDKit✔️✔️:     std::swap(nbrPerm, probePerm);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   boost::dynamic_bitset<> swapped(probe.size());
    // RDKit✔️✔️:   for (unsigned int i = 0; i < probePerm.size() - 1; ++i) {
    // RDKit✔️✔️:     auto pval = probePerm[i];
    // RDKit✔️✔️:     if (nbrPerm[i] == pval) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto tgt = std::find(nbrPerm.begin() + i, nbrPerm.end(), pval);
    // RDKit✔️✔️:     TEST_ASSERT(tgt != nbrPerm.end());
    // RDKit✔️✔️:     perm = swap_func(perm, i, tgt - nbrPerm.begin());
    // RDKit✔️✔️:     std::swap(*tgt, nbrPerm[i]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return perm;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::getChiralPermutation
    if permutation == 0 {
        return Ok(0);
    }
    let swap: fn(u32, usize, usize) -> u32 = match chiral_tag {
        ChiralTag::SquarePlanar if probe.len() <= 4 => swap_squareplanar,
        ChiralTag::TrigonalBipyramidal if probe.len() <= 5 => swap_trigonalbipyramidal,
        ChiralTag::Octahedral if probe.len() <= 6 => swap_octahedral,
        ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral => {
            return Ok(0);
        }
        _ => return Ok(0),
    };

    let mut order = vec![-1_isize; bond_count];
    for (neighbor_index, bond) in incident_bonds.iter().enumerate() {
        let Some(slot) = order.get_mut(bond.index()) else {
            return Err(SmilesStereoOrderError::BondIndexOutOfRange {
                bond_index: bond.index(),
                bond_count,
            });
        };
        *slot = neighbor_index as isize;
    }
    let mut neighbor_permutation = (0..incident_bonds.len())
        .map(|value| value as isize)
        .collect::<Vec<_>>();
    let mut probe_permutation = Vec::with_capacity(probe.len());
    for bond in probe {
        let position = match bond {
            Some(bond) => {
                *order
                    .get(bond.index())
                    .ok_or(SmilesStereoOrderError::BondIndexOutOfRange {
                        bond_index: bond.index(),
                        bond_count,
                    })?
            }
            None => -1,
        };
        probe_permutation.push(position);
    }
    if neighbor_permutation.len() < probe_permutation.len() {
        neighbor_permutation.extend(std::iter::repeat_n(
            -1,
            probe_permutation.len() - neighbor_permutation.len(),
        ));
    }
    if neighbor_permutation.len() != probe_permutation.len() {
        return Err(SmilesStereoOrderError::ProbeLengthMismatch {
            incident_count: neighbor_permutation.len(),
            probe_count: probe_permutation.len(),
        });
    }
    if inverse {
        std::mem::swap(&mut neighbor_permutation, &mut probe_permutation);
    }
    for index in 0..probe_permutation.len().saturating_sub(1) {
        let probe_value = probe_permutation[index];
        if neighbor_permutation[index] == probe_value {
            continue;
        }
        let target = neighbor_permutation[index..]
            .iter()
            .position(|value| *value == probe_value)
            .map(|offset| index + offset)
            .ok_or(SmilesStereoOrderError::ProbeMemberMissing { probe_index: index })?;
        permutation = swap(permutation, index, target);
        neighbor_permutation.swap(index, target);
    }
    Ok(permutation)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bond(index: usize) -> BondId {
        BondId::new(index)
    }

    #[test]
    fn source_swap_tables_cover_every_permutation_and_transposition() {
        for (permutation, row) in SWAP_SQUAREPLANAR_TABLE.iter().enumerate() {
            let mut swap_index = 0;
            for x in 0..4 {
                for y in x + 1..4 {
                    assert_eq!(
                        swap_squareplanar(permutation as u32, x, y),
                        u32::from(row[swap_index])
                    );
                    assert_eq!(
                        swap_squareplanar(permutation as u32, y, x),
                        u32::from(row[swap_index])
                    );
                    swap_index += 1;
                }
            }
        }
        for (permutation, row) in SWAP_TRIGONALBIPYRAMIDAL_TABLE.iter().enumerate() {
            let mut swap_index = 0;
            for x in 0..5 {
                for y in x + 1..5 {
                    assert_eq!(
                        swap_trigonalbipyramidal(permutation as u32, x, y),
                        u32::from(row[swap_index])
                    );
                    assert_eq!(
                        swap_trigonalbipyramidal(permutation as u32, y, x),
                        u32::from(row[swap_index])
                    );
                    swap_index += 1;
                }
            }
        }
        for (permutation, row) in SWAP_OCTAHEDRAL_TABLE.iter().enumerate() {
            let mut swap_index = 0;
            for x in 0..6 {
                for y in x + 1..6 {
                    assert_eq!(
                        swap_octahedral(permutation as u32, x, y),
                        u32::from(row[swap_index])
                    );
                    assert_eq!(
                        swap_octahedral(permutation as u32, y, x),
                        u32::from(row[swap_index])
                    );
                    swap_index += 1;
                }
            }
        }
    }

    #[test]
    fn source_swap_functions_return_zero_for_invalid_indices_and_permutations() {
        assert_eq!(swap_squareplanar(4, 0, 1), 0);
        assert_eq!(swap_trigonalbipyramidal(21, 0, 1), 0);
        assert_eq!(swap_octahedral(31, 0, 1), 0);
        assert_eq!(swap_squareplanar(1, 0, 4), 0);
        assert_eq!(swap_squareplanar(1, 4, 0), 0);
        assert_eq!(swap_trigonalbipyramidal(1, 0, 5), 0);
        assert_eq!(swap_octahedral(1, 6, 0), 0);

        // RDKit returns the unchanged value before validating the permutation
        // when both ligand indices are identical.
        assert_eq!(swap_squareplanar(99, 2, 2), 99);
        assert_eq!(swap_trigonalbipyramidal(99, 3, 3), 99);
        assert_eq!(swap_octahedral(99, 4, 4), 99);
    }

    #[test]
    fn typed_order_errors_keep_structural_failures_distinct() {
        assert_eq!(
            nontetrahedral_chiral_permutation(
                1,
                ChiralTag::SquarePlanar,
                1,
                &[bond(1)],
                &[Some(bond(1))],
                false,
            ),
            Err(SmilesStereoOrderError::BondIndexOutOfRange {
                bond_index: 1,
                bond_count: 1,
            })
        );
        assert_eq!(
            nontetrahedral_chiral_permutation(
                1,
                ChiralTag::SquarePlanar,
                1,
                &[bond(0)],
                &[Some(bond(1))],
                false,
            ),
            Err(SmilesStereoOrderError::BondIndexOutOfRange {
                bond_index: 1,
                bond_count: 1,
            })
        );
        assert_eq!(
            nontetrahedral_chiral_permutation(
                1,
                ChiralTag::SquarePlanar,
                2,
                &[bond(0), bond(1)],
                &[Some(bond(0))],
                false,
            ),
            Err(SmilesStereoOrderError::ProbeLengthMismatch {
                incident_count: 2,
                probe_count: 1,
            })
        );
        assert_eq!(
            nontetrahedral_chiral_permutation(
                1,
                ChiralTag::SquarePlanar,
                3,
                &[bond(0), bond(1)],
                &[Some(bond(2)), Some(bond(0))],
                false,
            ),
            Err(SmilesStereoOrderError::ProbeMemberMissing { probe_index: 0 })
        );
        // The pinned source examines all but the final position. A duplicated
        // member left in that final position therefore returns the accumulated
        // permutation instead of raising its TEST_ASSERT branch.
        assert_eq!(
            nontetrahedral_chiral_permutation(
                1,
                ChiralTag::SquarePlanar,
                2,
                &[bond(0), bond(1)],
                &[Some(bond(0)), Some(bond(0))],
                false,
            ),
            Ok(1)
        );
    }

    #[test]
    fn nontetrahedral_permutation_preserves_source_zero_and_inverse_rules() {
        let incident = [bond(0), bond(1), bond(2), bond(3)];
        let probe = [Some(bond(2)), Some(bond(0)), Some(bond(3)), Some(bond(1))];
        let forward = nontetrahedral_chiral_permutation(
            1,
            ChiralTag::SquarePlanar,
            4,
            &incident,
            &probe,
            false,
        )
        .unwrap();
        let restored = nontetrahedral_chiral_permutation(
            forward,
            ChiralTag::SquarePlanar,
            4,
            &incident,
            &probe,
            true,
        )
        .unwrap();
        assert_eq!(restored, 1);

        assert_eq!(
            nontetrahedral_chiral_permutation(
                0,
                ChiralTag::SquarePlanar,
                0,
                &[bond(9)],
                &[Some(bond(9))],
                false,
            ),
            Ok(0)
        );
        assert_eq!(
            nontetrahedral_chiral_permutation(
                1,
                ChiralTag::TetrahedralCw,
                4,
                &incident,
                &probe,
                false,
            ),
            Ok(0)
        );
        assert_eq!(
            nontetrahedral_chiral_permutation(
                1,
                ChiralTag::SquarePlanar,
                5,
                &[bond(0), bond(1), bond(2), bond(3), bond(4)],
                &[
                    Some(bond(0)),
                    Some(bond(1)),
                    Some(bond(2)),
                    Some(bond(3)),
                    Some(bond(4)),
                ],
                false,
            ),
            Ok(0)
        );

        let with_implicit = [Some(bond(0)), None, Some(bond(1)), None];
        assert!(
            nontetrahedral_chiral_permutation(
                1,
                ChiralTag::SquarePlanar,
                2,
                &[bond(0), bond(1)],
                &with_implicit,
                false,
            )
            .is_ok()
        );
    }
}
