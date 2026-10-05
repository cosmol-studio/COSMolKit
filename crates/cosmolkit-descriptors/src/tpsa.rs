//! TPSA descriptor owner (RDKit `MolSurf.cpp` TPSA family).

use cosmolkit_core::ValenceAssignment;
use cosmolkit_model::TopologyBlock;
use cosmolkit_model::{AtomId, BondOrder};

use crate::labute::LabuteCacheSlot;
use crate::{DescriptorError, DescriptorInput, DescriptorResult};

/// One include-sulfur/phosphorus slot of the detached TPSA computed-state
/// cache: the two source property guards (`_tpsa-<bool>` scalar and
/// `_tpsaAtomContribs-<bool>` rows) as INDEPENDENT typed presence, with no
/// string keys and no raw public mutation.
#[derive(Clone, Default, Debug, PartialEq)]
pub(crate) struct TpsaCacheSlot {
    pub(crate) scalar: Option<f64>,
    pub(crate) contributions: Option<Vec<f64>>,
}

/// Detached computed-state carrier for the descriptor domain: default
/// empty, `Clone` preserves bits/row order and independent presence per
/// slot; an explicit [`Self::clear`] empties every slot. This mirrors
/// the source molecule's computed-property store WITHOUT molecule or
/// cache-runtime authority: the caller owns the value.
#[derive(Clone, Default, Debug, PartialEq)]
pub struct DescriptorComputedState {
    tpsa_conventional: TpsaCacheSlot,
    tpsa_include_sand_p: TpsaCacheSlot,
    labute: LabuteCacheSlot,
    crippen: CrippenCacheSlot,
    pub(crate) chi_v_weights: Option<Vec<f64>>,
    pub(crate) chi_n_weights: Option<Vec<f64>>,
}

impl DescriptorComputedState {
    /// Explicit invalidation: empties BOTH TPSA slots (each slot's scalar
    /// and contribution presence) AND all THREE Labute presences (rows,
    /// hydrogens, asa) AND all FOUR Crippen presences (LogP rows, MR
    /// rows, LogP scalar, MR scalar). No hidden invalidation exists.
    pub fn clear(&mut self) {
        self.tpsa_conventional = TpsaCacheSlot::default();
        self.tpsa_include_sand_p = TpsaCacheSlot::default();
        self.labute = LabuteCacheSlot::default();
        self.crippen = CrippenCacheSlot::default();
        self.chi_v_weights = None;
        self.chi_n_weights = None;
    }

    fn slot(&self, include_sulfur_phosphorus: bool) -> &TpsaCacheSlot {
        if include_sulfur_phosphorus {
            &self.tpsa_include_sand_p
        } else {
            &self.tpsa_conventional
        }
    }

    fn slot_mut(&mut self, include_sulfur_phosphorus: bool) -> &mut TpsaCacheSlot {
        if include_sulfur_phosphorus {
            &mut self.tpsa_include_sand_p
        } else {
            &mut self.tpsa_conventional
        }
    }

    // Owning-crate test accessors for seeding/inspecting private slots;
    // not exported outside the crate.
    #[cfg(test)]
    pub(crate) fn slot_for_tests(&self, include_sulfur_phosphorus: bool) -> &TpsaCacheSlot {
        self.slot(include_sulfur_phosphorus)
    }

    #[cfg(test)]
    pub(crate) fn slot_mut_for_tests(
        &mut self,
        include_sulfur_phosphorus: bool,
    ) -> &mut TpsaCacheSlot {
        self.slot_mut(include_sulfur_phosphorus)
    }

    pub(crate) fn labute_slot(&self) -> &LabuteCacheSlot {
        &self.labute
    }

    pub(crate) fn labute_slot_mut(&mut self) -> &mut LabuteCacheSlot {
        &mut self.labute
    }

    pub(crate) fn crippen_slot(&self) -> &CrippenCacheSlot {
        &self.crippen
    }

    pub(crate) fn crippen_slot_mut(&mut self) -> &mut CrippenCacheSlot {
        &mut self.crippen
    }
}

/// The source's FOUR independent Crippen computed properties
/// (Crippen.cpp:83-85 contributions, 120-121 scalars), detached into
/// one private slot with independent presence per field and NO
/// includeHs key (the source caches are keyed only by the property
/// names, never by the includeHs flag).
#[derive(Clone, Default, Debug, PartialEq)]
pub(crate) struct CrippenCacheSlot {
    /// `_crippenLogPContribs` per-atom rows.
    pub logp_rows: Option<Vec<f64>>,
    /// `_crippenMRContribs` per-atom rows.
    pub mr_rows: Option<Vec<f64>>,
    /// `_crippenLogP` scalar total.
    pub logp: Option<f64>,
    /// `_crippenMR` scalar total.
    pub mr: Option<f64>,
}

/// Cached TPSA atom contributions (the source `getTPSAAtomContribs`
/// entrypoint, MolSurf.cpp:103-345).
///
/// Behavior review: tests `!force && contributions-present` FIRST — on a
/// hit it clones the cached rows and reads the cached scalar; a present
/// contribution cache with a MISSING scalar preserves the source's
/// property-read failure as the typed
/// [`DescriptorError::MissingComputedScalar`] (no zero/recompute
/// fallback). Otherwise it borrows the supplied FINAL topology/valence/
/// rings from the [`DescriptorInput`] (with `prepared_valence`'s borrowed
/// row-shape check), runs the ONE existing kernel with zero-initialized
/// rows, and publishes BOTH cache values only after success. Valence
/// failures keep their known typed cause and leave the state unchanged.
///
/// Complexity review: the guard is O(1) presence checks plus the cached
/// `rows.clone()` (the source's `getProp(contribsName, Vi)` cache-vector
/// copy), in addition to the fresh output vector allocated on the cold
/// arm and the kernel's six fact arrays; the cold arm is the
/// already-audited kernel cost plus one Vec allocation (source
/// `contribs.resize`) and the two cache stores (source `setProp`
/// writebacks at MolSurf.cpp:341-343). No reassignment, no property
/// strings, no second kernel.
pub fn tpsa_contributions(
    input: &DescriptorInput<'_>,
    include_sulfur_phosphorus: bool,
    force: bool,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<Vec<f64>> {
    // RDKit source (MolSurf.cpp:108-116):
    //   if (!force && mol.hasProp(contribsName)) {
    //     mol.getProp(contribsName, Vi);
    //     mol.getProp(pname, res);
    //     return res;
    //   }
    // RDKit✔️✔️:   if (!force && mol.hasProp(contribsName)) {
    // RDKit✔️✔️:     mol.getProp(contribsName, Vi);
    // RDKit✔️✔️:     mol.getProp(pname, res);
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    if !force {
        if let Some(cached) = state.slot(include_sulfur_phosphorus).contributions.clone() {
            let scalar = state.slot(include_sulfur_phosphorus).scalar.ok_or(
                DescriptorError::MissingComputedScalar {
                    function: "tpsa_contributions",
                    include_sulfur_phosphorus,
                },
            )?;
            let _ = scalar;
            return Ok(cached);
        }
    }

    let topology = input.topology();
    // Borrowed shape check on the supplied FINAL valence rows.
    crate::prepared_valence(topology, Some(input.valence()), "tpsa_contributions")?;
    let mut rows = vec![0.0f64; topology.atoms.len()];
    let res = tpsa_atom_contribs_kernel(
        topology,
        input.valence(),
        input.ring_info(),
        include_sulfur_phosphorus,
        &mut rows,
    )?;

    // RDKit source (MolSurf.cpp:341-343):
    //   mol.setProp(contribsName, Vi, true);
    //   mol.setProp(pname, res, true);
    //   return res;
    // RDKit✔️✔️:   mol.setProp(contribsName, Vi, true);
    // RDKit✔️✔️:   mol.setProp(pname, res, true);
    // RDKit✔️✔️:   return res;
    let slot = state.slot_mut(include_sulfur_phosphorus);
    slot.contributions = Some(rows.clone());
    slot.scalar = Some(res);
    Ok(rows)
}

/// Cached TPSA scalar dispatch (the source `calcTPSA`, MolSurf.cpp:347-360).
///
/// Behavior review: tests `!force && scalar-present` FIRST — returns that
/// scalar with NO vector clone and NO chemistry recomputation (the
/// source's warm `getProp(pname, res)` read). Otherwise calls the ONE
/// contribution owner and reads its resulting published scalar; the
/// source's cold arm likewise forwards into `getTPSAAtomContribs`.
/// The scalar is never re-derived by summing cached rows.
///
/// Complexity review: the warm arm is a single O(1) typed presence read
/// (source property read) that clones NO vector; the cold arm is the
/// contribution owner's audited cost — including its `rows.clone()`
/// cache-vector copy, fresh output vector and six fact arrays — plus
/// one O(1) scalar read. No second kernel, no summation pass.
pub(crate) fn tpsa_scalar(
    input: &DescriptorInput<'_>,
    include_sulfur_phosphorus: bool,
    force: bool,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit source (MolSurf.cpp:347-360):
    //   double calcTPSA(const ROMol &mol, bool force, bool includeSandP) {
    //     std::string pname =
    //         (boost::format("%s-%s") % common_properties::_tpsa % includeSandP).str();
    //     if (!force && mol.hasProp(pname)) {
    //       double res;
    //       mol.getProp(pname, res);
    //       return res;
    //     }
    //     std::vector<double> contribs;
    //     contribs.resize(mol.getNumAtoms());
    //     double res;
    //     res = getTPSAAtomContribs(mol, contribs, force, includeSandP);
    //     return res;
    //   }
    // RDKit✔️✔️:   double calcTPSA(const ROMol &mol, bool force, bool includeSandP) {
    // RDKit✔️✔️:     std::string pname =
    // RDKit✔️✔️:         (boost::format("%s-%s") % common_properties::_tpsa % includeSandP).str();
    // RDKit✔️✔️:     if (!force && mol.hasProp(pname)) {
    // RDKit✔️✔️:       double res;
    // RDKit✔️✔️:       mol.getProp(pname, res);
    // RDKit✔️✔️:       return res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::vector<double> contribs;
    // RDKit✔️✔️:     contribs.resize(mol.getNumAtoms());
    // RDKit✔️✔️:     double res;
    // RDKit✔️✔️:     res = getTPSAAtomContribs(mol, contribs, force, includeSandP);
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    if !force {
        if let Some(res) = state.slot(include_sulfur_phosphorus).scalar {
            return Ok(res);
        }
    }
    tpsa_contributions(input, include_sulfur_phosphorus, force, state)?;
    state
        .slot(include_sulfur_phosphorus)
        .scalar
        .ok_or(DescriptorError::MissingComputedScalar {
            function: "tpsa",
            include_sulfur_phosphorus,
        })
}

/// The ONE TPSA atom-contribution kernel shared by the entrypoints.
///
/// Behavior review: reproduces `getTPSAAtomContribs` exactly — ONE pass
/// over ALL bonds filling per-atom fact arrays with EXPLICIT-H
/// PRECEDENCE (a bond whose begin OR end atom is hydrogen NEVER reaches
/// the aromatic/order arms: begin==H => nNbrs[end] -= 1, nHs[end] += 1;
/// else end==H symmetric; else aromatic => nArom both; else
/// SINGLE/DOUBLE/TRIPLE increment both endpoints; every other order
/// counts NOTHING). `contribs` must be pre-sized to nAtoms (the real
/// caller zero-fills); skipped atoms keep their entries — matching the
/// zero-filled vector semantics of `calcTPSA`'s sole call path. Per
/// atom of the element filter {N(7), O(8)} (+{P(15), S(16)} iff
/// `include_sand_p`): nHs += `total_hydrogen_count_from_validated`
/// with include_neighbors=FALSE — the source `getTotalNumHs(false)`
/// counts atom-spec explicit Hs PLUS effective implicit Hs and EXCLUDES
/// H-atom neighbors, which the bond pass already tallied (no double
/// counting); nNbrs +=
/// degree from the adjacency row (net nNbrs is the HEAVY degree: the
/// H decrements cancel the H neighbors' degree contribution — the
/// typed `atom->getDegree()` projection); chg = formal_charge; in3Ring
/// = `RingInfo::is_atom_in_ring_of_size(i, 3)`; then the N/O exact
/// tables with the clamped linear fallbacks and the P/S tables (no
/// fallback, tmp starts 0.0). Vi[i] = tmp; res += tmp. This pure
/// kernel owns NO cache: the molecule-level property guards
/// (`_tpsa-<bool>` / `_tpsaAtomContribs-<bool>`, force flag,
/// setProp writebacks) are owned by the current T07 typed entrypoints
/// [`tpsa_contributions`] and [`tpsa_scalar`] on the detached
/// [`DescriptorComputedState`].
///
/// Complexity review: two passes over the input — O(bonds) fact pass +
/// O(atoms) table pass — with ONE fixed set of six O(atoms) i32 arrays
/// (the source's single `std::vector` allocation set); the table
/// branches are straight-line else-if chains preserved as-is; no sort,
/// no map. The 3-ring lookup scans the SUPPLIED indexed ring-membership
/// list per qualifying atom, like the source's
/// `isAtomInRingOfSize(i, 3)`; its summed membership cost comes IN
/// ADDITION to the O(bonds+atoms) passes. `contribs` is caller-owned
/// storage written in place; the returned total is a plain f64, not a
/// heap allocation; no allocation beyond the six arrays and the
/// `Result`.
#[allow(clippy::too_many_lines)]
pub(crate) fn tpsa_atom_contribs_kernel(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    ring_info: &cosmolkit_core::RingInfo,
    include_sand_p: bool,
    contribs: &mut [f64],
) -> DescriptorResult<f64> {
    // RDKit source (MolSurf.cpp:118-150):
    //   unsigned int nAtoms = mol.getNumAtoms();
    //   std::vector<int> nNbrs(nAtoms, 0), nSing(nAtoms, 0), nDoub(nAtoms, 0),
    //       nTrip(nAtoms, 0), nArom(nAtoms, 0), nHs(nAtoms, 0);
    //   for (ROMol::ConstBondIterator bIt = mol.beginBonds(); bIt != mol.endBonds();
    //        ++bIt) {
    //     const Bond *bnd = (*bIt);
    //     if (bnd->getBeginAtom()->getAtomicNum() == 1) {
    //       nNbrs[bnd->getEndAtomIdx()] -= 1;
    //       nHs[bnd->getEndAtomIdx()] += 1;
    //     } else if (bnd->getEndAtom()->getAtomicNum() == 1) {
    //       nNbrs[bnd->getBeginAtomIdx()] -= 1;
    //       nHs[bnd->getBeginAtomIdx()] += 1;
    //     } else if (bnd->getIsAromatic()) {
    //       nArom[bnd->getBeginAtomIdx()] += 1;
    //       nArom[bnd->getEndAtomIdx()] += 1;
    //     } else {
    //       switch (bnd->getBondType()) {
    //         case Bond::SINGLE:
    //           nSing[bnd->getBeginAtomIdx()] += 1;
    //           nSing[bnd->getEndAtomIdx()] += 1;
    //           break;
    //         case Bond::DOUBLE:
    //           nDoub[bnd->getBeginAtomIdx()] += 1;
    //           nDoub[bnd->getEndAtomIdx()] += 1;
    //           break;
    //         case Bond::TRIPLE:
    //           nTrip[bnd->getBeginAtomIdx()] += 1;
    //           nTrip[bnd->getEndAtomIdx()] += 1;
    //           break;
    //         default:
    //           break;
    //       }
    //     }
    //   }
    // RDKit✔️✔️: unsigned int nAtoms = mol.getNumAtoms();
    // RDKit✔️✔️:   std::vector<int> nNbrs(nAtoms, 0), nSing(nAtoms, 0), nDoub(nAtoms, 0),
    // RDKit✔️✔️:       nTrip(nAtoms, 0), nArom(nAtoms, 0), nHs(nAtoms, 0);
    // RDKit✔️✔️:   for (ROMol::ConstBondIterator bIt = mol.beginBonds(); bIt != mol.endBonds();
    // RDKit✔️✔️:        ++bIt) {
    // RDKit✔️✔️:     const Bond *bnd = (*bIt);
    // RDKit✔️✔️:     if (bnd->getBeginAtom()->getAtomicNum() == 1) {
    // RDKit✔️✔️:       nNbrs[bnd->getEndAtomIdx()] -= 1;
    // RDKit✔️✔️:       nHs[bnd->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:     } else if (bnd->getEndAtom()->getAtomicNum() == 1) {
    // RDKit✔️✔️:       nNbrs[bnd->getBeginAtomIdx()] -= 1;
    // RDKit✔️✔️:       nHs[bnd->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:     } else if (bnd->getIsAromatic()) {
    // RDKit✔️✔️:       nArom[bnd->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:       nArom[bnd->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       switch (bnd->getBondType()) {
    // RDKit✔️✔️:         case Bond::SINGLE:
    // RDKit✔️✔️:           nSing[bnd->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:           nSing[bnd->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Bond::DOUBLE:
    // RDKit✔️✔️:           nDoub[bnd->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:           nDoub[bnd->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Bond::TRIPLE:
    // RDKit✔️✔️:           nTrip[bnd->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:           nTrip[bnd->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let n_atoms = topology.atoms.len();
    debug_assert_eq!(
        contribs.len(),
        n_atoms,
        "caller pre-sizes contribs (TEST_ASSERT Vi)"
    );
    let mut n_nbrs = vec![0i32; n_atoms];
    let mut n_sing = vec![0i32; n_atoms];
    let mut n_doub = vec![0i32; n_atoms];
    let mut n_trip = vec![0i32; n_atoms];
    let mut n_arom = vec![0i32; n_atoms];
    let mut n_hs = vec![0i32; n_atoms];
    for bond in &topology.bonds {
        let begin = bond.begin().index();
        let end = bond.end().index();
        if topology.atoms[begin].element().atomic_number() == 1 {
            n_nbrs[end] -= 1;
            n_hs[end] += 1;
        } else if topology.atoms[end].element().atomic_number() == 1 {
            n_nbrs[begin] -= 1;
            n_hs[begin] += 1;
        } else if bond.is_aromatic() {
            n_arom[begin] += 1;
            n_arom[end] += 1;
        } else {
            match bond.order() {
                BondOrder::Single => {
                    n_sing[begin] += 1;
                    n_sing[end] += 1;
                }
                BondOrder::Double => {
                    n_doub[begin] += 1;
                    n_doub[end] += 1;
                }
                BondOrder::Triple => {
                    n_trip[begin] += 1;
                    n_trip[end] += 1;
                }
                _ => {}
            }
        }
    }

    // RDKit source (MolSurf.cpp:152-166):
    //   for (unsigned int i = 0; i < nAtoms; ++i) {
    //     const Atom *atom = mol.getAtomWithIdx(i);
    //     int atNum = atom->getAtomicNum();
    //
    //     if (atNum != 7 && atNum != 8 &&
    //         (!includeSandP || (atNum != 15 && atNum != 16))) {
    //       continue;
    //     }
    //
    //     nHs[i] += atom->getTotalNumHs();
    //     int chg = atom->getFormalCharge();
    //     bool in3Ring = mol.getRingInfo()->isAtomInRingOfSize(i, 3);
    //     nNbrs[i] += atom->getDegree();
    //
    //     double tmp = -1;
    // RDKit✔️✔️: for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️✔️:   const Atom *atom = mol.getAtomWithIdx(i);
    // RDKit✔️✔️:   int atNum = atom->getAtomicNum();
    // RDKit✔️✔️:   if (atNum != 7 && atNum != 8 &&
    // RDKit✔️✔️:       (!includeSandP || (atNum != 15 && atNum != 16))) {
    // RDKit✔️✔️:     continue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   nHs[i] += atom->getTotalNumHs();
    // RDKit✔️✔️:   int chg = atom->getFormalCharge();
    // RDKit✔️✔️:   bool in3Ring = mol.getRingInfo()->isAtomInRingOfSize(i, 3);
    // RDKit✔️✔️:   nNbrs[i] += atom->getDegree();
    // RDKit✔️✔️:   double tmp = -1;
    //
    // The COMPLETE N/O/P/S branch bodies are copied verbatim below,
    // adjacent to their corresponding Rust branches.
    let mut res = 0.0f64;
    for (i, atom) in topology.atoms.iter().enumerate() {
        let at_num = atom.element().atomic_number();
        if at_num != 7 && at_num != 8 && (!include_sand_p || (at_num != 15 && at_num != 16)) {
            continue;
        }
        // include_neighbors=false (getTotalNumHs(false)): atom-spec
        // explicit Hs PLUS effective implicit Hs; H-atom neighbors are
        // excluded here because the bond pass already tallied them.
        n_hs[i] += i32::try_from(
            cosmolkit_core::total_hydrogen_count_from_validated(
                topology,
                valence,
                AtomId::new(i),
                false,
            )
            .map_err(|source| DescriptorError::Valence {
                function: "tpsa_atom_contribs",
                source,
            })?,
        )
        .map_err(|_| DescriptorError::CountOverflow {
            function: "tpsa_atom_contribs",
            field: "hydrogen_count",
        })?;
        let chg = i32::from(atom.formal_charge());
        let in3_ring = ring_info.is_atom_in_ring_of_size(AtomId::new(i), 3);
        // Typed heavy-degree projection: the adjacency row length
        // counts every neighbor; the explicit-H decrements above cancel
        // the H rows, so net nNbrs[i] is the HEAVY degree — exactly
        // atom->getDegree() after the source's H decrements.
        n_nbrs[i] += i32::try_from(topology.adjacency.neighbors_of(i).len()).map_err(|_| {
            DescriptorError::CountOverflow {
                function: "tpsa_atom_contribs",
                field: "degree",
            }
        })?;

        let mut tmp = -1.0f64;
        if at_num == 7 {
            // RDKit source (MolSurf.cpp:168-246) — complete N body:
            // RDKit✔️✔️:     if (atNum == 7) {
            // RDKit✔️✔️:       switch (nNbrs[i]) {
            // RDKit✔️✔️:         case 1:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nTrip[i] == 1) {
            // RDKit✔️✔️:             tmp = 23.79;
            // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 23.85;
            // RDKit✔️✔️:           } else if (nHs[i] == 2 && chg == 0 && nSing[i] == 1) {
            // RDKit✔️✔️:             tmp = 26.02;
            // RDKit✔️✔️:           } else if (nHs[i] == 2 && chg == 1 && nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 25.59;
            // RDKit✔️✔️:           } else if (nHs[i] == 3 && chg == 1 && nSing[i] == 1) {
            // RDKit✔️✔️:             tmp = 27.64;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         case 2:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 1 && nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 12.36;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nTrip[i] == 1 &&
            // RDKit✔️✔️:                      nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 13.60;
            // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nSing[i] == 2 && in3Ring) {
            // RDKit✔️✔️:             tmp = 21.94;
            // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nSing[i] == 2 && !in3Ring) {
            // RDKit✔️✔️:             tmp = 12.03;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 1 && nTrip[i] == 1 &&
            // RDKit✔️✔️:                      nSing[i] == 1) {
            // RDKit✔️✔️:             tmp = 4.36;
            // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 1 && nDoub[i] == 1 &&
            // RDKit✔️✔️:                      nSing[i] == 1) {
            // RDKit✔️✔️:             tmp = 13.97;
            // RDKit✔️✔️:           } else if (nHs[i] == 2 && chg == 1 && nSing[i] == 2) {
            // RDKit✔️✔️:             tmp = 16.61;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nArom[i] == 2) {
            // RDKit✔️✔️:             tmp = 12.89;
            // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nArom[i] == 2) {
            // RDKit✔️✔️:             tmp = 15.79;
            // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 1 && nArom[i] == 2) {
            // RDKit✔️✔️:             tmp = 14.14;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         case 3:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 3 && in3Ring) {
            // RDKit✔️✔️:             tmp = 3.01;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nSing[i] == 3 && !in3Ring) {
            // RDKit✔️✔️:             tmp = 3.24;
            //
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nSing[i] == 1 &&
            // RDKit✔️✔️:                      nDoub[i] == 2) {
            // RDKit✔️✔️:             tmp = 11.68;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 1 && nSing[i] == 2 &&
            // RDKit✔️✔️:                      nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 3.01;
            // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 1 && nSing[i] == 3) {
            // RDKit✔️✔️:             tmp = 4.44;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nArom[i] == 3) {
            // RDKit✔️✔️:             tmp = 4.41;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nSing[i] == 1 &&
            // RDKit✔️✔️:                      nArom[i] == 2) {
            // RDKit✔️✔️:             tmp = 4.93;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nDoub[i] == 1 &&
            // RDKit✔️✔️:                      nArom[i] == 2) {
            // RDKit✔️✔️:             tmp = 8.39;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 1 && nArom[i] == 3) {
            // RDKit✔️✔️:             tmp = 4.10;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 1 && nSing[i] == 1 &&
            // RDKit✔️✔️:                      nArom[i] == 2) {
            // RDKit✔️✔️:             tmp = 3.88;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           if (nHs[i] == 0 && nSing[i] == 4 && chg == 1) {
            // RDKit✔️✔️:             tmp = 0.0;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         default:
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       if (tmp < 0.0) {
            // RDKit✔️✔️:         tmp = 30.5 - nNbrs[i] * 8.2 + nHs[i] * 1.5;
            // RDKit✔️✔️:         if (tmp < 0) {
            // RDKit✔️✔️:           tmp = 0.0;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            match n_nbrs[i] {
                1 => {
                    if n_hs[i] == 0 && chg == 0 && n_trip[i] == 1 {
                        tmp = 23.79;
                    } else if n_hs[i] == 1 && chg == 0 && n_doub[i] == 1 {
                        tmp = 23.85;
                    } else if n_hs[i] == 2 && chg == 0 && n_sing[i] == 1 {
                        tmp = 26.02;
                    } else if n_hs[i] == 2 && chg == 1 && n_doub[i] == 1 {
                        tmp = 25.59;
                    } else if n_hs[i] == 3 && chg == 1 && n_sing[i] == 1 {
                        tmp = 27.64;
                    }
                }
                2 => {
                    if n_hs[i] == 0 && chg == 0 && n_sing[i] == 1 && n_doub[i] == 1 {
                        tmp = 12.36;
                    } else if n_hs[i] == 0 && chg == 0 && n_trip[i] == 1 && n_doub[i] == 1 {
                        tmp = 13.60;
                    } else if n_hs[i] == 1 && chg == 0 && n_sing[i] == 2 && in3_ring {
                        tmp = 21.94;
                    } else if n_hs[i] == 1 && chg == 0 && n_sing[i] == 2 && !in3_ring {
                        tmp = 12.03;
                    } else if n_hs[i] == 0 && chg == 1 && n_trip[i] == 1 && n_sing[i] == 1 {
                        tmp = 4.36;
                    } else if n_hs[i] == 1 && chg == 1 && n_doub[i] == 1 && n_sing[i] == 1 {
                        tmp = 13.97;
                    } else if n_hs[i] == 2 && chg == 1 && n_sing[i] == 2 {
                        tmp = 16.61;
                    } else if n_hs[i] == 0 && chg == 0 && n_arom[i] == 2 {
                        tmp = 12.89;
                    } else if n_hs[i] == 1 && chg == 0 && n_arom[i] == 2 {
                        tmp = 15.79;
                    } else if n_hs[i] == 1 && chg == 1 && n_arom[i] == 2 {
                        tmp = 14.14;
                    }
                }
                3 => {
                    if n_hs[i] == 0 && chg == 0 && n_sing[i] == 3 && in3_ring {
                        tmp = 3.01;
                    } else if n_hs[i] == 0 && chg == 0 && n_sing[i] == 3 && !in3_ring {
                        tmp = 3.24;
                    } else if n_hs[i] == 0 && chg == 0 && n_sing[i] == 1 && n_doub[i] == 2 {
                        tmp = 11.68;
                    } else if n_hs[i] == 0 && chg == 1 && n_sing[i] == 2 && n_doub[i] == 1 {
                        tmp = 3.01;
                    } else if n_hs[i] == 1 && chg == 1 && n_sing[i] == 3 {
                        tmp = 4.44;
                    } else if n_hs[i] == 0 && chg == 0 && n_arom[i] == 3 {
                        tmp = 4.41;
                    } else if n_hs[i] == 0 && chg == 0 && n_sing[i] == 1 && n_arom[i] == 2 {
                        tmp = 4.93;
                    } else if n_hs[i] == 0 && chg == 0 && n_doub[i] == 1 && n_arom[i] == 2 {
                        tmp = 8.39;
                    } else if n_hs[i] == 0 && chg == 1 && n_arom[i] == 3 {
                        tmp = 4.10;
                    } else if n_hs[i] == 0 && chg == 1 && n_sing[i] == 1 && n_arom[i] == 2 {
                        tmp = 3.88;
                    }
                }
                4 => {
                    if n_hs[i] == 0 && n_sing[i] == 4 && chg == 1 {
                        tmp = 0.0;
                    }
                }
                _ => {}
            }
            if tmp < 0.0 {
                tmp = 30.5 - f64::from(n_nbrs[i]) * 8.2 + f64::from(n_hs[i]) * 1.5;
                if tmp < 0.0 {
                    tmp = 0.0;
                }
            }
        } else if at_num == 8 {
            // RDKit source (MolSurf.cpp:252-284) — complete O body:
            // RDKit✔️✔️:     } else if (atNum == 8) {
            // RDKit✔️✔️:       switch (nNbrs[i]) {
            // RDKit✔️✔️:         case 1:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 17.07;
            // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nSing[i] == 1) {
            // RDKit✔️✔️:             tmp = 20.23;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == -1 && nSing[i] == 1) {
            // RDKit✔️✔️:             tmp = 23.06;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         case 2:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 2 && in3Ring) {
            // RDKit✔️✔️:             tmp = 12.53;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nSing[i] == 2 && !in3Ring) {
            // RDKit✔️✔️:             tmp = 9.23;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nArom[i] == 2) {
            // RDKit✔️✔️:             tmp = 13.14;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         default:
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       if (tmp < 0.0) {
            // RDKit✔️✔️:         tmp = 28.5 - nNbrs[i] * 8.6 + nHs[i] * 1.5;
            // RDKit✔️✔️:         if (tmp < 0) {
            // RDKit✔️✔️:           tmp = 0.0;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            match n_nbrs[i] {
                1 => {
                    if n_hs[i] == 0 && chg == 0 && n_doub[i] == 1 {
                        tmp = 17.07;
                    } else if n_hs[i] == 1 && chg == 0 && n_sing[i] == 1 {
                        tmp = 20.23;
                    } else if n_hs[i] == 0 && chg == -1 && n_sing[i] == 1 {
                        tmp = 23.06;
                    }
                }
                2 => {
                    if n_hs[i] == 0 && chg == 0 && n_sing[i] == 2 && in3_ring {
                        tmp = 12.53;
                    } else if n_hs[i] == 0 && chg == 0 && n_sing[i] == 2 && !in3_ring {
                        tmp = 9.23;
                    } else if n_hs[i] == 0 && chg == 0 && n_arom[i] == 2 {
                        tmp = 13.14;
                    }
                }
                _ => {}
            }
            if tmp < 0.0 {
                tmp = 28.5 - f64::from(n_nbrs[i]) * 8.6 + f64::from(n_hs[i]) * 1.5;
                if tmp < 0.0 {
                    tmp = 0.0;
                }
            }
        } else if include_sand_p && at_num == 15 {
            // RDKit source (MolSurf.cpp:285-306) — complete P body:
            // RDKit✔️✔️:     } else if (includeSandP && atNum == 15) {
            // RDKit✔️✔️:       tmp = 0.0;
            // RDKit✔️✔️:       switch (nNbrs[i]) {
            // RDKit✔️✔️:         case 2:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 1 && nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 34.14;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         case 3:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 3) {
            // RDKit✔️✔️:             tmp = 13.59;
            // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nSing[i] == 2 &&
            // RDKit✔️✔️:                      nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 23.47;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 3 && nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 9.81;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         default:
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:       }
            tmp = 0.0;
            match n_nbrs[i] {
                2 => {
                    if n_hs[i] == 0 && chg == 0 && n_sing[i] == 1 && n_doub[i] == 1 {
                        tmp = 34.14;
                    }
                }
                3 => {
                    if n_hs[i] == 0 && chg == 0 && n_sing[i] == 3 {
                        tmp = 13.59;
                    } else if n_hs[i] == 1 && chg == 0 && n_sing[i] == 2 && n_doub[i] == 1 {
                        tmp = 23.47;
                    }
                }
                4 => {
                    if n_hs[i] == 0 && chg == 0 && n_sing[i] == 3 && n_doub[i] == 1 {
                        tmp = 9.81;
                    }
                }
                _ => {}
            }
        } else if include_sand_p && at_num == 16 {
            // RDKit source (MolSurf.cpp:307-337) — complete S body:
            // RDKit✔️✔️:     } else if (includeSandP && atNum == 16) {
            // RDKit✔️✔️:       tmp = 0.0;
            // RDKit✔️✔️:       switch (nNbrs[i]) {
            // RDKit✔️✔️:         case 1:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 32.09;
            // RDKit✔️✔️:           } else if (nHs[i] == 1 && chg == 0 && nSing[i] == 1) {
            // RDKit✔️✔️:             tmp = 38.80;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         case 2:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 2) {
            // RDKit✔️✔️:             tmp = 25.30;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nArom[i] == 2) {
            // RDKit✔️✔️:             tmp = 28.24;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         case 3:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nArom[i] == 2 && nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 21.70;
            // RDKit✔️✔️:           } else if (nHs[i] == 0 && chg == 0 && nSing[i] == 2 &&
            // RDKit✔️✔️:                      nDoub[i] == 1) {
            // RDKit✔️✔️:             tmp = 19.21;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           if (nHs[i] == 0 && chg == 0 && nSing[i] == 2 && nDoub[i] == 2) {
            // RDKit✔️✔️:             tmp = 8.38;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         default:
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     }
            tmp = 0.0;
            match n_nbrs[i] {
                1 => {
                    if n_hs[i] == 0 && chg == 0 && n_doub[i] == 1 {
                        tmp = 32.09;
                    } else if n_hs[i] == 1 && chg == 0 && n_sing[i] == 1 {
                        tmp = 38.80;
                    }
                }
                2 => {
                    if n_hs[i] == 0 && chg == 0 && n_sing[i] == 2 {
                        tmp = 25.30;
                    } else if n_hs[i] == 0 && chg == 0 && n_arom[i] == 2 {
                        tmp = 28.24;
                    }
                }
                3 => {
                    if n_hs[i] == 0 && chg == 0 && n_arom[i] == 2 && n_doub[i] == 1 {
                        tmp = 21.70;
                    } else if n_hs[i] == 0 && chg == 0 && n_sing[i] == 2 && n_doub[i] == 1 {
                        tmp = 19.21;
                    }
                }
                4 => {
                    if n_hs[i] == 0 && chg == 0 && n_sing[i] == 2 && n_doub[i] == 2 {
                        tmp = 8.38;
                    }
                }
                _ => {}
            }
        }
        // RDKit source (MolSurf.cpp:331-333):
        //   Vi[i] = tmp;
        //   res += tmp;
        // RDKit✔️✔️: Vi[i] = tmp;
        // RDKit✔️✔️: res += tmp;
        contribs[i] = tmp;
        res += tmp;
    }
    Ok(res)
}
