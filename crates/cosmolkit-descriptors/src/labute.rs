//! Labute ASA descriptor owner (RDKit `MolSurf.cpp` Labute family).

use crate::{DescriptorComputedState, DescriptorError, DescriptorInput, DescriptorResult};

/// Per-atom Labute ASA contributions plus the implicit-hydrogen term.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LabuteContributions {
    /// Per-atom ASA contributions (source `Vi`).
    pub atoms: Vec<f64>,
    /// Implicit-hydrogen contribution (source `hContrib`).
    pub hydrogens: f64,
}

/// One Labute cache slot: THREE independent typed presences modeling the
/// three source properties written together at MolSurf.cpp:82-84 — the
/// contribution rows (`_labuteAtomContribs`), the hydrogen term
/// (`_labuteAtomHContrib`) and the scalar total (`_labuteASA`).
///
/// The warm guard keys on `rows` alone (the source `hasProp
/// (_labuteAtomContribs)`); the reads then follow the source getProp order
/// rows -> hydrogens -> asa, failing with family-specific errors on partial
/// presence exactly where the source `getProp` would throw.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct LabuteCacheSlot {
    pub(crate) rows: Option<Vec<f64>>,
    pub(crate) hydrogens: Option<f64>,
    pub(crate) asa: Option<f64>,
}

/// Cached Labute atom contributions (the source
/// `getLabuteAtomContribs`, MolSurf.cpp:25-86).
///
/// Behavior review: reproduces the source exactly — warm guard keyed on
/// the rows presence alone (source `hasProp(_labuteAtomContribs)`, NOT
/// includeHs-keyed); the warm arm reads the THREE cached values in source
/// getProp order (`Vi`, then `hContrib`, then `_labuteASA`) with
/// family-specific errors on partial presence; the cold arm
/// reuses the existing `cosmolkit_core::rdkit_rb0` owner for every radius,
/// zeroes `Vi`, applies the bond-overlap pass with the source
/// `bondScaleFacts` indexed BY THE C++ ENUM VALUE under
/// `!aromatic && bondType < 4` (SINGLE=1 -> 0.0, DOUBLE=2 -> 0.2,
/// TRIPLE=3 -> 0.3; aromatic -> facts[0]=0.1; type>=4 no scaling),
/// clamps `dij = min(max(|Ri-Rj|, bij), Ri+Rj)`, accumulates
/// `Vi[begin] += Rj*Rj - (Ri-dij)*(Ri-dij)/dij` /
/// `Vi[end] += Ri*Ri - (Rj-dij)*(Rj-dij)/dij` per bond in iteration
/// order; `include_hydrogens` adds the all-atom pass with
/// `Rj = rdkit_rb0(1)` feeding both `Vi[i]` and `hContrib`; the final
/// pass associates `Vi[i] = PI*Ri*(4*Ri - Vi[i])`, `res += Vi[i]`, and
/// the `|hContrib| > 1e-4` gate transforms/adds `hContrib` only above
/// threshold (below: untransformed, still cached, not added to res).
/// Publication-after-success writes the three slot values.
///
/// Complexity review: the warm arm's presence checks are O(1) but
/// serving the cached rows CLONES the contribution vector — O(V), the
/// source `getProp` cache-vector copy. The cold arm owns the output
/// `Vi` allocation, the `rads` Vec allocation AND the cache-vector copy
/// written into the slot, plus three linear passes (radius O(V), bonds
/// O(E) with O(1) arithmetic, final O(V)) and one extra O(V) pass when
/// `include_hydrogens`. No sort, no map.
pub fn labute_contributions(
    input: &DescriptorInput<'_>,
    include_hydrogens: bool,
    force: bool,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<LabuteContributions> {
    // RDKit source (MolSurf.cpp:26-35):
    //   TEST_ASSERT(Vi.size() == mol.getNumAtoms());
    //   if (!force && mol.hasProp(common_properties::_labuteAtomContribs)) {
    //     mol.getProp(common_properties::_labuteAtomContribs, Vi);
    //     mol.getProp(common_properties::_labuteAtomHContrib, hContrib);
    //     double res;
    //     mol.getProp(common_properties::_labuteASA, res);
    //     return res;
    //   }
    // RDKit✔️✔️:   TEST_ASSERT(Vi.size() == mol.getNumAtoms());
    // RDKit✔️✔️:   if (!force && mol.hasProp(common_properties::_labuteAtomContribs)) {
    // RDKit✔️✔️:     mol.getProp(common_properties::_labuteAtomContribs, Vi);
    // RDKit✔️✔️:     mol.getProp(common_properties::_labuteAtomHContrib, hContrib);
    // RDKit✔️✔️:     double res;
    // RDKit✔️✔️:     mol.getProp(common_properties::_labuteASA, res);
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    if !force {
        if let Some(rows) = state.labute_slot().rows.clone() {
            // Read order per source getProps 30-33: hydrogens BEFORE asa.
            let hydrogens =
                state
                    .labute_slot()
                    .hydrogens
                    .ok_or(DescriptorError::MissingLabuteHydrogens {
                        function: "labute_contributions",
                    })?;
            // The source reads and returns res here; the typed owner
            // requires its presence (same failure point) without
            // re-exposing the scalar through this row API.
            let _asa = state
                .labute_slot()
                .asa
                .ok_or(DescriptorError::MissingLabuteAsa {
                    function: "labute_contributions",
                })?;
            return Ok(LabuteContributions {
                atoms: rows,
                hydrogens,
            });
        }
    }

    let topology = input.topology();
    let n_atoms = topology.atoms.len();
    // RDKit source (MolSurf.cpp:36-42):
    //   unsigned int nAtoms = mol.getNumAtoms();
    //   std::vector<double> rads(nAtoms);
    //   for (unsigned int i = 0; i < nAtoms; ++i) {
    //     rads[i] = PeriodicTable::getTable()->getRb0(
    //         mol.getAtomWithIdx(i)->getAtomicNum());
    //     Vi[i] = 0.0;
    //   }
    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit✔️✔️:   std::vector<double> rads(nAtoms);
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️✔️:     rads[i] = PeriodicTable::getTable()->getRb0(
    // RDKit✔️✔️:         mol.getAtomWithIdx(i)->getAtomicNum());
    // RDKit✔️✔️:     Vi[i] = 0.0;
    // RDKit✔️✔️:   }
    let mut rads = Vec::with_capacity(n_atoms);
    let mut vi = vec![0.0f64; n_atoms];
    for atom in &topology.atoms {
        rads.push(cosmolkit_core::rdkit_rb0(atom.element().atomic_number()));
    }

    // RDKit source (MolSurf.cpp:44-60):
    //   for (ROMol::ConstBondIterator bondIt = mol.beginBonds();
    //        bondIt != mol.endBonds(); ++bondIt) {
    //     const double bondScaleFacts[4] = {.1, 0, .2, .3};
    //     double Ri = rads[(*bondIt)->getBeginAtomIdx()];
    //     double Rj = rads[(*bondIt)->getEndAtomIdx()];
    //     double bij = Ri + Rj;
    //     if (!(*bondIt)->getIsAromatic()) {
    //       if ((*bondIt)->getBondType() < 4) {
    //         bij -= bondScaleFacts[(*bondIt)->getBondType()];
    //       }
    //     } else {
    //       bij -= bondScaleFacts[0];
    //     }
    //     double dij = std::min(std::max(fabs(Ri - Rj), bij), Ri + Rj);
    //     Vi[(*bondIt)->getBeginAtomIdx()] += Rj * Rj - (Ri - dij) * (Ri - dij) / dij;
    //     Vi[(*bondIt)->getEndAtomIdx()] += Ri * Ri - (Rj - dij) * (Rj - dij) / dij;
    //   }
    // RDKit✔️✔️:     const double bondScaleFacts[4] = {.1, 0, .2, .3};
    // RDKit✔️✔️:     double Ri = rads[(*bondIt)->getBeginAtomIdx()];
    // RDKit✔️✔️:     double Rj = rads[(*bondIt)->getEndAtomIdx()];
    // RDKit✔️✔️:     double bij = Ri + Rj;
    // RDKit✔️✔️:     if (!(*bondIt)->getIsAromatic()) {
    // RDKit✔️✔️:       if ((*bondIt)->getBondType() < 4) {
    // RDKit✔️✔️:         bij -= bondScaleFacts[(*bondIt)->getBondType()];
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       bij -= bondScaleFacts[0];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     double dij = std::min(std::max(fabs(Ri - Rj), bij), Ri + Rj);
    // RDKit✔️✔️:     Vi[(*bondIt)->getBeginAtomIdx()] += Rj * Rj - (Ri - dij) * (Ri - dij) / dij;
    // RDKit✔️✔️:     Vi[(*bondIt)->getEndAtomIdx()] += Ri * Ri - (Rj - dij) * (Rj - dij) / dij;
    // Typed enum mapping of `getBondType() < 4` + facts[enumValue]
    // (Bond.h:55-79 / vocabulary.rs:226-250): UNSPECIFIED=0 ->
    // facts[0] = 0.1 (the historical mol_surface.rs:221 branch — NOT
    // ZERO=21, which is neither UNSPECIFIED nor an absent bond);
    // SINGLE=1 -> facts[1] = 0.0; DOUBLE=2 -> facts[2] = 0.2; TRIPLE=3
    // -> facts[3] = 0.3; every OTHER non-aromatic enum value 4..=21
    // subtracts NOTHING; ANY aromatic flag subtracts facts[0] = 0.1
    // regardless of enum.
    for bond in &topology.bonds {
        let ri = rads[bond.begin().index()];
        let rj = rads[bond.end().index()];
        let mut bij = ri + rj;
        use cosmolkit_model::BondOrder;
        if !bond.is_aromatic() {
            match bond.order() {
                BondOrder::Unspecified => bij -= 0.1,
                BondOrder::Single => bij -= 0.0,
                BondOrder::Double => bij -= 0.2,
                BondOrder::Triple => bij -= 0.3,
                _ => {}
            }
        } else {
            bij -= 0.1;
        }
        let dij = (ri - rj).abs().max(bij).min(ri + rj);
        vi[bond.begin().index()] += rj * rj - (ri - dij) * (ri - dij) / dij;
        vi[bond.end().index()] += ri * ri - (rj - dij) * (rj - dij) / dij;
    }

    // RDKit source (MolSurf.cpp:61-71):
    //   hContrib = 0.0;
    //   if (includeHs) {
    //     double Rj = PeriodicTable::getTable()->getRb0(1);
    //     for (unsigned int i = 0; i < nAtoms; ++i) {
    //       double Ri = rads[i];
    //       double bij = Ri + Rj;
    //       double dij = std::min(std::max(fabs(Ri - Rj), bij), Ri + Rj);
    //       Vi[i] += Rj * Rj - (Ri - dij) * (Ri - dij) / dij;
    //       hContrib += Ri * Ri - (Rj - dij) * (Rj - dij) / dij;
    //     }
    //   }
    // RDKit✔️✔️:   hContrib = 0.0;
    // RDKit✔️✔️:   if (includeHs) {
    // RDKit✔️✔️:     double Rj = PeriodicTable::getTable()->getRb0(1);
    // RDKit✔️✔️:     for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️✔️:       double Ri = rads[i];
    // RDKit✔️✔️:       double bij = Ri + Rj;
    // RDKit✔️✔️:       double dij = std::min(std::max(fabs(Ri - Rj), bij), Ri + Rj);
    // RDKit✔️✔️:       Vi[i] += Rj * Rj - (Ri - dij) * (Ri - dij) / dij;
    // RDKit✔️✔️:       hContrib += Ri * Ri - (Rj - dij) * (Rj - dij) / dij;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut h_contrib = 0.0f64;
    if include_hydrogens {
        let rj = cosmolkit_core::rdkit_rb0(1);
        for i in 0..n_atoms {
            let ri = rads[i];
            let bij = ri + rj;
            let dij = (ri - rj).abs().max(bij).min(ri + rj);
            vi[i] += rj * rj - (ri - dij) * (ri - dij) / dij;
            h_contrib += ri * ri - (rj - dij) * (rj - dij) / dij;
        }
    }

    // RDKit source (MolSurf.cpp:72-81):
    //   double res = 0.0;
    //   for (unsigned int i = 0; i < nAtoms; ++i) {
    //     double Ri = rads[i];
    //     Vi[i] = M_PI * Ri * (4. * Ri - Vi[i]);
    //     res += Vi[i];
    //   }
    //   if (includeHs && fabs(hContrib) > 1e-4) {
    //     double Rj = PeriodicTable::getTable()->getRb0(1);
    //     hContrib = M_PI * Rj * (4. * Rj - hContrib);
    //     res += hContrib;
    //   }
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️✔️:     double Ri = rads[i];
    // RDKit✔️✔️:     Vi[i] = M_PI * Ri * (4. * Ri - Vi[i]);
    // RDKit✔️✔️:     res += Vi[i];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (includeHs && fabs(hContrib) > 1e-4) {
    // RDKit✔️✔️:     double Rj = PeriodicTable::getTable()->getRb0(1);
    // RDKit✔️✔️:     hContrib = M_PI * Rj * (4. * Rj - hContrib);
    // RDKit✔️✔️:     res += hContrib;
    // RDKit✔️✔️:   }
    let mut res = 0.0f64;
    for i in 0..n_atoms {
        let ri = rads[i];
        vi[i] = std::f64::consts::PI * ri * (4.0 * ri - vi[i]);
        res += vi[i];
    }
    if include_hydrogens && h_contrib.abs() > 1e-4 {
        let rj = cosmolkit_core::rdkit_rb0(1);
        h_contrib = std::f64::consts::PI * rj * (4.0 * rj - h_contrib);
        res += h_contrib;
    }

    // RDKit source (MolSurf.cpp:82-84):
    //   mol.setProp(common_properties::_labuteAtomContribs, Vi, true);
    //   mol.setProp(common_properties::_labuteAtomHContrib, hContrib, true);
    //   mol.setProp(common_properties::_labuteASA, res, true);
    // RDKit✔️✔️:   mol.setProp(common_properties::_labuteAtomContribs, Vi, true);
    // RDKit✔️✔️:   mol.setProp(common_properties::_labuteAtomHContrib, hContrib, true);
    // RDKit✔️✔️:   mol.setProp(common_properties::_labuteASA, res, true);
    let slot = state.labute_slot_mut();
    slot.rows = Some(vi.clone());
    slot.hydrogens = Some(h_contrib);
    slot.asa = Some(res);
    Ok(LabuteContributions {
        atoms: vi,
        hydrogens: h_contrib,
    })
}

/// Cached Labute ASA total (the source `_labuteASA` scalar on the
/// contribution owner's cache).
///
/// Behavior review: warm scalar-first read (no row clone, no
/// recomputation) whenever `!force` and the scalar is present; the cold
/// arm calls the ONE contribution owner and reads its published scalar.
/// The scalar is never re-derived by summation.
///
/// Complexity review: warm arm is a single O(1) typed read; cold arm is
/// the contribution owner's audited cost plus one O(1) read.
pub fn labute_asa(
    input: &DescriptorInput<'_>,
    include_hydrogens: bool,
    force: bool,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit source (MolSurf.cpp:89-101):
    //   double calcLabuteASA(const ROMol &mol, bool includeHs, bool force) {
    //     if (!force && mol.hasProp(common_properties::_labuteASA)) {
    //       double res;
    //       mol.getProp(common_properties::_labuteASA, res);
    //       return res;
    //     }
    //     std::vector<double> contribs;
    //     contribs.resize(mol.getNumAtoms());
    //     double hContrib;
    //     double res;
    //     res = getLabuteAtomContribs(mol, contribs, hContrib, includeHs, force);
    //     return res;
    //   }
    // RDKit✔️✔️:   if (!force && mol.hasProp(common_properties::_labuteASA)) {
    // RDKit✔️✔️:     double res;
    // RDKit✔️✔️:     mol.getProp(common_properties::_labuteASA, res);
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // Typed mapping: `hasProp(_labuteASA)` is the asa Option presence; the
    // scalar-only warm read touches no row vector, exactly like the source.
    if !force {
        if let Some(res) = state.labute_slot().asa {
            return Ok(res);
        }
    }
    // RDKit✔️✔️:   std::vector<double> contribs;
    // RDKit✔️✔️:   contribs.resize(mol.getNumAtoms());
    // RDKit✔️✔️:   double hContrib;
    // RDKit✔️✔️:   double res;
    // RDKit✔️✔️:   res = getLabuteAtomContribs(mol, contribs, hContrib, includeHs, force);
    // RDKit✔️✔️:   return res;
    // Typed mapping: the source's caller-side out-params (`contribs`
    // resized to nAtoms, `hContrib`) exist only to satisfy the owner's
    // signature and are DISCARDED here; at this boundary the ONE
    // contribution owner returns an owned `LabuteContributions` which
    // this projection likewise discards. `return res` returns the
    // owner's `_labuteASA` — identical to the source's owner-return by
    // publish-after-success — read from the typed slot below. ZERO area
    // arithmetic lives in this projection.
    labute_contributions(input, include_hydrogens, force, state)?;
    state
        .labute_slot()
        .asa
        .ok_or(DescriptorError::MissingLabuteAsa {
            function: "labute_asa",
        })
}
