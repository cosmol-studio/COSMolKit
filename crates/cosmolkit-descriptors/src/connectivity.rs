//! Hall–Kier alpha, Kappa shape indices and Phi over detached values.

use cosmolkit_core::{PathSearchParams, all_paths_of_length, rdkit_rb0};
use cosmolkit_model::{Atom, Element, Hybridization, TopologyBlock};

use crate::{DescriptorError, DescriptorResult};

/// Pinned RDKit Hall–Kier alpha descriptor version.
pub const HALL_KIER_ALPHA_VERSION: &str = "1.2.0";

/// Calculates the scalar Hall–Kier alpha and optionally writes per-atom terms.
///
/// Uses each atom's stored hybridization verbatim. This function performs no
/// sanitization, hybridization assignment, valence calculation or cache work.
/// The atom slice's order defines contribution indices, independently of IDs.
///
/// A supplied sink must contain at least `atoms.len()` entries. Non-wildcard
/// entries are overwritten, wildcard entries and any excess tail are preserved.
/// Consequently, a prefilled sink need not sum to the returned scalar. An
/// undersized sink produces a structured error before any entry is changed.
pub fn hall_kier_alpha(
    atoms: &[Atom],
    mut atom_contribs: Option<&mut [f64]>,
) -> DescriptorResult<f64> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8:
    // Code/GraphMol/Descriptors/ConnectivityDescriptors.cpp:267-295.
    // RDKit✔️✔️: double calcHallKierAlpha(const ROMol &mol, std::vector<double> *atomContribs) {
    // RDKit✔️✔️:   PRECONDITION(!atomContribs || atomContribs->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atomContribs vector");
    // RDKit✔️✔️:   const PeriodicTable *tbl = PeriodicTable::getTable();
    // RDKit✔️✔️:   double alphaSum = 0.0;
    // RDKit✔️✔️:   double rC = tbl->getRb0(6);
    // RDKit✔️✔️:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️✔️:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️✔️:   while (atBegin != atEnd) {
    // RDKit✔️✔️:     const Atom *at = mol[*atBegin];
    // RDKit✔️✔️:     ++atBegin;
    // RDKit✔️✔️:     unsigned int n = at->getAtomicNum();
    // RDKit✔️✔️:     if (!n) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bool found;
    // RDKit✔️✔️:     double alpha = detail::getAlpha(*(at), found);
    // RDKit✔️✔️:     if (!found) {
    // RDKit✔️✔️:       double rA = tbl->getRb0(n);
    // RDKit✔️✔️:       alpha = rA / rC - 1.0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     alphaSum += alpha;
    // RDKit✔️✔️:     if (atomContribs) {
    // RDKit✔️✔️:       (*atomContribs)[at->getIdx()] = alpha;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return alphaSum;
    // RDKit✔️✔️: };
    // RDKit✔️✔️:
    // Input rows are the detached equivalent of ROMol's contiguous atom
    // indices. Element bounds the table index to 0..118. Borrowing a slice
    // exposes no bonds, valence, coordinates or runtime state.
    // Cost: source-shaped O(n) sequential traversal, O(1) workspace, no
    // descriptor allocations/clones; periodic data reuse the core singleton.
    if let Some(sink) = atom_contribs.as_ref() {
        if sink.len() < atoms.len() {
            return Err(DescriptorError::InvalidHallKierContributionRows {
                actual: sink.len(),
                minimum: atoms.len(),
            });
        }
    }
    let mut alpha_sum = 0.0;
    let carbon_radius = rdkit_rb0(6);
    for (index, atom) in atoms.iter().enumerate() {
        let atomic_number = atom.element().atomic_number();
        if atomic_number == 0 {
            continue;
        }
        let alpha = match get_alpha(atom) {
            Some(alpha) => alpha,
            None => rdkit_rb0(atomic_number) / carbon_radius - 1.0,
        };
        alpha_sum += alpha;
        if let Some(sink) = atom_contribs.as_mut() {
            sink[index] = alpha;
        }
    }
    Ok(alpha_sum)
}

fn get_alpha(atom: &Atom) -> Option<f64> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8:
    // Code/GraphMol/Descriptors/ConnectivityDescriptors.cpp:71-167.
    // RDKit✔️✔️: double getAlpha(const Atom &atom, bool &found) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   found = false;
    // RDKit✔️✔️:   switch (atom.getAtomicNum()) {
    // RDKit✔️✔️:     case 1:
    // RDKit✔️✔️:       res = 0.0;
    // RDKit✔️✔️:       found = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 6:
    // RDKit✔️✔️:       switch (atom.getHybridization()) {
    // RDKit✔️✔️:         case Atom::SP:
    // RDKit✔️✔️:           res = -0.22;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Atom::SP2:
    // RDKit✔️✔️:           res = -0.13;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           res = 0.00;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:       };
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 7:
    // RDKit✔️✔️:       switch (atom.getHybridization()) {
    // RDKit✔️✔️:         case Atom::SP:
    // RDKit✔️✔️:           res = -0.29;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Atom::SP2:
    // RDKit✔️✔️:           res = -0.20;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           res = -0.04;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       };
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 8:
    // RDKit✔️✔️:       switch (atom.getHybridization()) {
    // RDKit✔️✔️:         case Atom::SP2:
    // RDKit✔️✔️:           res = -0.20;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           res = -0.04;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       };
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 9:
    // RDKit✔️✔️:       res = -0.07;
    // RDKit✔️✔️:       found = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 15:
    // RDKit✔️✔️:       switch (atom.getHybridization()) {
    // RDKit✔️✔️:         case Atom::SP2:
    // RDKit✔️✔️:           res = 0.30;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           res = 0.43;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       };
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 16:
    // RDKit✔️✔️:       switch (atom.getHybridization()) {
    // RDKit✔️✔️:         case Atom::SP2:
    // RDKit✔️✔️:           res = 0.22;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           res = 0.35;
    // RDKit✔️✔️:           found = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       };
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 17:
    // RDKit✔️✔️:       res = 0.29;
    // RDKit✔️✔️:       found = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 35:
    // RDKit✔️✔️:       res = 0.48;
    // RDKit✔️✔️:       found = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 53:
    // RDKit✔️✔️:       res = 0.73;
    // RDKit✔️✔️:       found = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: }  // namespace detail
    // Option encodes the source's found flag; None requests its radius path.
    // Both implementations perform bounded constant-time switches, without
    // allocation, repeated scans or inferred chemical state.
    Some(match atom.element() {
        Element::H => 0.0,
        Element::C => match atom.hybridization() {
            Hybridization::Sp => -0.22,
            Hybridization::Sp2 => -0.13,
            _ => 0.0,
        },
        Element::N => match atom.hybridization() {
            Hybridization::Sp => -0.29,
            Hybridization::Sp2 => -0.20,
            _ => -0.04,
        },
        Element::O => match atom.hybridization() {
            Hybridization::Sp2 => -0.20,
            _ => -0.04,
        },
        Element::F => -0.07,
        Element::P => match atom.hybridization() {
            Hybridization::Sp2 => 0.30,
            _ => 0.43,
        },
        Element::S => match atom.hybridization() {
            Hybridization::Sp2 => 0.22,
            _ => 0.35,
        },
        Element::CL => 0.29,
        Element::BR => 0.48,
        Element::I => 0.73,
        _ => return None,
    })
}

/// Pinned RDKit Kappa1 descriptor version.
pub const KAPPA_1_VERSION: &str = "1.1.0";
/// Pinned RDKit Kappa2 descriptor version.
pub const KAPPA_2_VERSION: &str = "1.1.0";
/// Pinned RDKit Kappa3 descriptor version.
pub const KAPPA_3_VERSION: &str = "1.1.0";
/// Pinned RDKit Kier Phi descriptor version.
pub const PHI_VERSION: &str = "1.0.0";

/// Calculates Kappa1 from all explicit bonds, heavy atoms and stored alpha.
///
/// Explicit hydrogen bonds enter the bond count; atom hydrogen properties do
/// not. The existing count owner validates topology and no chemistry is
/// prepared or changed. A zero denominator returns positive zero.
pub fn kappa_1(topology: &TopologyBlock) -> DescriptorResult<f64> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8,
    // ConnectivityDescriptors.cpp:327-333.
    // RDKit✔️❌: double calcKappa1(const ROMol &mol) {
    // RDKit✔️❌:   double P1 = mol.getNumBonds();
    // RDKit✔️❌:   double A = mol.getNumHeavyAtoms();
    // RDKit✔️❌:   double alpha = calcHallKierAlpha(mol);
    // RDKit✔️❌:   double kappa = kappa1Helper(P1, A, alpha);
    // RDKit✔️❌:   return kappa;
    // RDKit✔️❌: }
    // The count owner's graph validation adds O(n+m+annotations) and
    // adjacency allocation versus source's already-valid O(n) atom count.
    // No descriptor clone, additional healthy validation or count scan.
    let p1 = explicit_bond_count(topology, "kappa_1")?;
    let a = f64::from(heavy_atom_count(topology, "kappa_1")?);
    let alpha = hall_kier_alpha(&topology.atoms, None)?;
    Ok(kappa_1_helper(p1, a, alpha))
}

/// Calculates Kappa2 using default bond paths of length two.
///
/// Paths exclude H/isotopic H, include wildcards, deduplicate bond sets and
/// include non-shortest paths. All enumeration belongs to the core owner;
/// invalid paths/topology retain their structural error. Stored alpha is used.
pub fn kappa_2(topology: &TopologyBlock) -> DescriptorResult<f64> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8,
    // ConnectivityDescriptors.cpp:334-341.
    // RDKit✔️❌: double calcKappa2(const ROMol &mol) {
    // RDKit✔️❌:   PATH_LIST ps = findAllPathsOfLengthN(mol, 2);
    // RDKit✔️❌:   double P2 = ps.size();
    // RDKit✔️❌:   double A = mol.getNumHeavyAtoms();
    // RDKit✔️❌:   double alpha = calcHallKierAlpha(mol);
    // RDKit✔️❌:   double kappa = kappa2Helper(P2, A, alpha);
    // RDKit✔️❌:   return kappa;
    // RDKit✔️❌: }
    // Existing path/count boundaries each validate their borrowed input;
    // this validation allocation is extra relative to source ROMol inputs.
    // The enumerator is reused, never cloned, replaced or repeated.
    let p2 = path_count(topology, 2, "kappa_2")?;
    let a = f64::from(heavy_atom_count(topology, "kappa_2")?);
    let alpha = hall_kier_alpha(&topology.atoms, None)?;
    Ok(kappa_2_helper(p2, a, alpha))
}

/// Calculates Kappa3 using default length-three bond paths and heavy parity.
///
/// Heavy-atom count follows the source's unsigned-u32 to signed-i32 narrowing.
/// Odd/even selection uses that signed value. No clamp or zero-heavy shortcut
/// is added; negative results and source-defined zero denominators are kept.
pub fn kappa_3(topology: &TopologyBlock) -> DescriptorResult<f64> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8,
    // ConnectivityDescriptors.cpp:342-348.
    // RDKit✔️❌: double calcKappa3(const ROMol &mol) {
    // RDKit✔️❌:   double P3 = findAllPathsOfLengthN(mol, 3).size();
    // RDKit✔️❌:   int A = mol.getNumHeavyAtoms();
    // RDKit✔️❌:   double alpha = calcHallKierAlpha(mol);
    // RDKit✔️❌:   double kappa = kappa3Helper(P3, A, alpha);
    // RDKit✔️❌:   return kappa;
    // RDKit✔️❌: }
    // Core path search plus count validation have the same additional
    // boundary cost as Kappa2; no new path/graph algorithm or clone.
    let p3 = path_count(topology, 3, "kappa_3")?;
    // Pinned source int is 32 bits. Rust as preserves the unsigned-to-signed
    // congruence, rather than saturating, rejecting or changing parity.
    let a = heavy_atom_count(topology, "kappa_3")? as i32;
    let alpha = hall_kier_alpha(&topology.atoms, None)?;
    Ok(kappa_3_helper(p3, a, alpha))
}

/// Calculates Kier Phi, returning positive zero for zero heavy atoms.
///
/// The heavy-zero guard occurs before alpha and path work. Positive-heavy
/// inputs use one alpha, one length-two enumeration and the shared numerical
/// Kappa helpers, preserving multiplication/division order and signed zero.
pub fn phi(topology: &TopologyBlock) -> DescriptorResult<f64> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8,
    // ConnectivityDescriptors.cpp:349-361.
    // RDKit✔️❌: double calcPhi(const ROMol &mol) {
    // RDKit✔️❌:   if (!mol.getNumHeavyAtoms()) {
    // RDKit✔️❌:     return 0.0;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   auto alpha = calcHallKierAlpha(mol);
    // RDKit✔️❌:   auto P1 = mol.getNumBonds();
    // RDKit✔️❌:   auto A = mol.getNumHeavyAtoms();
    // RDKit✔️❌:   auto kappa1 = kappa1Helper(P1, A, alpha);
    // RDKit✔️❌:   auto P2 = findAllPathsOfLengthN(mol, 2).size();
    // RDKit✔️❌:   auto kappa2 = kappa2Helper(P2, A, alpha);
    // RDKit✔️❌:   auto Phi = kappa1 * kappa2 / A;
    // RDKit✔️❌:   return Phi;
    // RDKit✔️❌: }
    // The installed count owner validates the detached boundary before
    // returning A; this allocates adjacency unlike source's valid ROMol.
    // Reuse the first count on immutable input: source's second count would
    // produce the same value, and a second validation/count is unnecessary.
    let heavy = heavy_atom_count(topology, "phi")?;
    if heavy == 0 {
        return Ok(0.0);
    }
    let alpha = hall_kier_alpha(&topology.atoms, None)?;
    let p1 = explicit_bond_count(topology, "phi")?;
    let a = f64::from(heavy);
    let kappa1 = kappa_1_helper(p1, a, alpha);
    let p2 = path_count(topology, 2, "phi")?;
    let kappa2 = kappa_2_helper(p2, a, alpha);
    Ok(kappa1 * kappa2 / a)
}

fn kappa_1_helper(p1: f64, a: f64, alpha: f64) -> f64 {
    // RDKit✔️✔️: double kappa1Helper(double P1, double A, double alpha) {
    // RDKit✔️✔️:   double denom = P1 + alpha;
    // RDKit✔️✔️:   double kappa = 0.0;
    // RDKit✔️✔️:   if (denom) {
    // RDKit✔️✔️:     kappa = (A + alpha) * (A + alpha - 1) * (A + alpha - 1) / (denom * denom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return kappa;
    // RDKit✔️✔️: }
    // Constant scalar work, no allocation; preserve source association.
    let denom = p1 + alpha;
    let mut kappa = 0.0;
    if denom != 0.0 {
        kappa = (a + alpha) * (a + alpha - 1.0) * (a + alpha - 1.0) / (denom * denom);
    }
    kappa
}

fn kappa_2_helper(p2: f64, a: f64, alpha: f64) -> f64 {
    // RDKit✔️✔️: double kappa2Helper(double P2, double A, double alpha) {
    // RDKit✔️✔️:   double denom = (P2 + alpha) * (P2 + alpha);
    // RDKit✔️✔️:   double kappa = 0.0;
    // RDKit✔️✔️:   if (denom) {
    // RDKit✔️✔️:     kappa = (A + alpha - 1) * (A + alpha - 2) * (A + alpha - 2) / denom;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return kappa;
    // RDKit✔️✔️: }
    // Constant scalar work, no allocation; preserve source association.
    let denom = (p2 + alpha) * (p2 + alpha);
    let mut kappa = 0.0;
    if denom != 0.0 {
        kappa = (a + alpha - 1.0) * (a + alpha - 2.0) * (a + alpha - 2.0) / denom;
    }
    kappa
}

fn kappa_3_helper(p3: f64, a: i32, alpha: f64) -> f64 {
    // RDKit✔️✔️: double kappa3Helper(double P3, int A, double alpha) {
    // RDKit✔️✔️:   double denom = (P3 + alpha) * (P3 + alpha);
    // RDKit✔️✔️:   double kappa = 0.0;
    // RDKit✔️✔️:   if (denom) {
    // RDKit✔️✔️:     if (A % 2) {
    // RDKit✔️✔️:       kappa = (A + alpha - 1) * (A + alpha - 3) * (A + alpha - 3) / denom;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       kappa = (A + alpha - 2) * (A + alpha - 3) * (A + alpha - 3) / denom;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return kappa;
    // RDKit✔️✔️: }
    // Constant scalar work, no allocation; signed remainder matches source.
    let denom = (p3 + alpha) * (p3 + alpha);
    let mut kappa = 0.0;
    if denom != 0.0 {
        if a % 2 != 0 {
            kappa = (f64::from(a) + alpha - 1.0)
                * (f64::from(a) + alpha - 3.0)
                * (f64::from(a) + alpha - 3.0)
                / denom;
        } else {
            kappa = (f64::from(a) + alpha - 2.0)
                * (f64::from(a) + alpha - 3.0)
                * (f64::from(a) + alpha - 3.0)
                / denom;
        }
    }
    kappa
}

fn heavy_atom_count(topology: &TopologyBlock, function: &'static str) -> DescriptorResult<u32> {
    // Source heavy-count behavior has exactly one owner in counts.rs.
    // ROMol.cpp:187-196 is delegated, never reimplemented here.
    // RDKit✔️❌: unsigned int ROMol::getNumHeavyAtoms() const {
    // RDKit✔️❌:   unsigned int res = 0;
    // RDKit✔️❌:   for (const auto atom : atoms()) {
    // RDKit✔️❌:     if (atom->getAtomicNum() > 1) {
    // RDKit✔️❌:       ++res;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: };
    // RDKit✔️❌:
    // Healthy input: ONE owner call and its existing O(n+m+annotations)
    // validation/allocation. Only on the owner's known two legacy flattened
    // errors do we recover the structural validation cause or source u32
    // count overflow. Other error variants are propagated without rewriting.
    crate::counts::num_heavy_atoms_kernel(topology).map_err(|error| match error {
        DescriptorError::Unsupported {
            function: "num_heavy_atoms",
            ..
        } => match topology.validate() {
            Err(source) => DescriptorError::InvalidTopology { function, source },
            Ok(()) => DescriptorError::CountOverflow {
                function,
                field: "heavy_atoms",
            },
        },
        other => other,
    })
}

fn explicit_bond_count(topology: &TopologyBlock, function: &'static str) -> DescriptorResult<f64> {
    // RDKit✔️✔️: unsigned int ROMol::getNumBonds(bool onlyHeavy) const {
    // RDKit✔️✔️:   // By default return the bonds that connect only the heavy atoms
    // RDKit✔️✔️:   // hydrogen connecting bonds are ignores
    // RDKit✔️✔️:   auto res = numBonds;
    // RDKit✔️✔️:   if (!onlyHeavy) {
    // RDKit✔️✔️:     // If we need hydrogen connecting bonds add them up
    // RDKit✔️✔️:     for (const auto atom : atoms()) {
    // RDKit✔️✔️:       res += atom->getTotalNumHs();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // Only the source default onlyHeavy=true arm is exposed in this packet.
    // numBonds tracks every explicit added edge (RWMol.cpp:524-526), even
    // H endpoints; it is not a heavy-endpoint filter. Row length is O(1),
    // source u32 width is checked and there is no H-property inference.
    let count =
        u32::try_from(topology.bonds.len()).map_err(|_| DescriptorError::CountOverflow {
            function,
            field: "bonds",
        })?;
    Ok(f64::from(count))
}

fn path_count(
    topology: &TopologyBlock,
    length: usize,
    function: &'static str,
) -> DescriptorResult<f64> {
    // Subgraphs.h:124-126:
    // RDKit✔️✔️: RDKIT_SUBGRAPHS_EXPORT PATH_LIST findAllPathsOfLengthN(
    // RDKit✔️✔️:     const ROMol &mol, unsigned int targetLen, bool useBonds = true,
    // RDKit✔️✔️:     bool useHs = false, int rootedAtAtom = -1, bool onlyShortestPaths = false);
    // Enumeration, ring closure, ordering and bond-set deduplication are the
    // existing core owner's source closure. This adapter adds no enumeration,
    // allocation or clone and releases its unneeded path rows after counting.
    all_paths_of_length(topology, length, &PathSearchParams::default())
        .map(|paths| paths.len() as f64)
        .map_err(|source| DescriptorError::Path { function, source })
}
