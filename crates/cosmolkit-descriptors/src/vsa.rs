//! VSA binning owner (RDKit `MolSurf.cpp` assignContribsToBins family).

use crate::{DescriptorComputedState, DescriptorError, DescriptorInput, DescriptorResult};

/// Accumulate per-atom contributions into property bins (the source
/// anonymous-namespace `assignContribsToBins`, MolSurf.cpp:362-377).
///
/// Behavior review: reproduces the source exactly — the two
/// PRECONDITIONs (contribs.len()==bin_prop.len() and a res of at least
/// bins.len()+1; the owned boundary returns a freshly allocated res of
/// EXACTLY bins.len()+1, satisfying the second by construction); the
/// accumulation loop runs in ATOM ORDER with per-atom `cVal`/`bVal` and
/// `res[idx] += cVal` (left-to-right +=, observable in the final bits
/// when several atoms hit one bin). `std::upper_bound` (the index of
/// the FIRST element GREATER than bVal) maps to Rust's
/// `partition_point(|&b| b <= bVal)` — the standard equivalent on a
/// partitioned range: equality lands AFTER the matching edge, a value
/// below all bins lands in bin 0, above all in the final slot. The
/// bins vector's monotonicity is the CALLER's precondition exactly as
/// in the source (upper_bound requires a partitioned range); NO sort
/// and NO validation are added, and a non-partitioned input yields the
/// partition point in both languages (documented, not a heuristic).
///
/// Complexity review: O(n log m) with n = atoms and m = bins — one
/// binary search and one O(1) += per atom, identical to the source;
/// the single allocation is the owned output vector; no map, no
/// per-atom allocation, no repeated scans.
pub fn assign_contribs_to_bins(
    contribs: &[f64],
    bin_prop: &[f64],
    bins: &[f64],
) -> DescriptorResult<Vec<f64>> {
    // RDKit source (MolSurf.cpp:362-377):
    //   void assignContribsToBins(const std::vector<double> &contribs,
    //                             const std::vector<double> &binProp,
    //                             const std::vector<double> &bins,
    //                             std::vector<double> &res) {
    //     PRECONDITION(contribs.size() == binProp.size(), "mismatched array sizes");
    //     PRECONDITION(res.size() >= bins.size() + 1, "mismatched array sizes");
    //     for (unsigned int i = 0; i < contribs.size(); ++i) {
    //       double cVal = contribs[i];
    //       double bVal = binProp[i];
    //       unsigned int idx =
    //           std::upper_bound(bins.begin(), bins.end(), bVal) - bins.begin();
    //       res[idx] += cVal;
    //     }
    //   }
    // RDKit✔️✔️:   PRECONDITION(contribs.size() == binProp.size(), "mismatched array sizes");
    // RDKit✔️✔️:   PRECONDITION(res.size() >= bins.size() + 1, "mismatched array sizes");
    // RDKit✔️✔️:   for (unsigned int i = 0; i < contribs.size(); ++i) {
    // RDKit✔️✔️:     double cVal = contribs[i];
    // RDKit✔️✔️:     double bVal = binProp[i];
    // RDKit✔️✔️:     unsigned int idx = std::upper_bound(bins.begin(), bins.end(), bVal) - bins.begin();
    // RDKit✔️✔️:     res[idx] += cVal;
    // RDKit✔️✔️:   }
    // Typed mapping: the PRECONDITION aborts map to the typed
    // MismatchedBinArrays error; the caller-supplied res maps to the
    // owned bins.len()+1 return value (the second precondition holds by
    // construction); upper_bound maps to partition_point(|&b| b <= bVal).
    if contribs.len() != bin_prop.len() {
        return Err(DescriptorError::MismatchedBinArrays {
            contribs_len: contribs.len(),
            bin_prop_len: bin_prop.len(),
            bins_len: bins.len(),
        });
    }
    let mut res = vec![0.0f64; bins.len() + 1];
    for i in 0..contribs.len() {
        let c_val = contribs[i];
        let b_val = bin_prop[i];
        let idx = bins.partition_point(|&b| b <= b_val);
        res[idx] += c_val;
    }
    Ok(res)
}

/// The frozen default SlogP_VSA bin boundaries (the source `blist[11]`,
/// MolSurf.cpp:382-383, copied verbatim).
pub(crate) const DEFAULT_SLOGP_BINS: [f64; 11] =
    [-0.4, -0.2, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6];

/// SlogP_VSA: Labute VSA contributions accumulated into Crippen-logP
/// bins (the source `calcSlogP_VSA`, MolSurf.cpp:379-404).
///
/// Behavior review: `bins == None` selects the FROZEN default
/// 11-boundary list (12 outputs); `Some` copies the caller's
/// boundaries (m boundaries -> m+1 outputs, empty -> 1). The VSA rows
/// come from the ONE Labute contribution owner with
/// `include_hydrogens = true` exactly as the source passes `true`;
/// the hydrogen lump is the discarded `tmp` scalar. The binning
/// property is the CRIPPEN logP row per atom (the MR rows are
/// discarded exactly as the source discards `mrContribs`), and the
/// accumulation itself is the V01 owner `assign_contribs_to_bins`
/// (equality-to-edge lands in the NEXT bin; sortedness remains the
/// caller's precondition). `force` is forwarded to BOTH owners' cache
/// arms: the Labute slot and the canonical `crippen_contributions`
/// entry, which reads the four-field Crippen cache slot in the
/// CALLER's state (warm hit when both row vectors match the atom
/// count; force=true bypasses to the ONE cold kernel) — the same
/// cache arm as `crippen_totals`.
///
/// Complexity review: two owner calls at their audited costs (Labute
/// O(V+E) cold passes and O(V) warm output-row copying; Crippen
/// per-atom pattern matching on cold and O(V) warm output-row copying)
/// plus the V01 O(n log m) binning. Owned
/// allocations: the lbins copy (m), the owned res (m+1), and the
/// owners' returned contribution vectors; the cold Crippen path also
/// stores cache row copies in the caller's state (the warm hit copies
/// them into the owned outputs without the kernel). No map, no sort,
/// no repeated scans.
pub fn slogp_vsa(
    input: &DescriptorInput<'_>,
    bins: Option<&[f64]>,
    force: bool,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<Vec<f64>> {
    // RDKit source (MolSurf.cpp:379-404):
    //   std::vector<double> calcSlogP_VSA(const ROMol &mol, std::vector<double> *bins,
    //                                     bool force) {
    //     // FIX: use force value to include caching
    //     std::vector<double> lbins;
    //     if (!bins) {
    //       double blist[11] = {-0.4, -0.2, 0,   0.1, 0.15, 0.2,
    //                           0.25, 0.3,  0.4, 0.5, 0.6};
    //       lbins.resize(11);
    //       std::copy(blist, blist + 11, lbins.begin());
    //     } else {
    //       lbins.resize(bins->size());
    //       std::copy(bins->begin(), bins->end(), lbins.begin());
    //     }
    //     std::vector<double> res(lbins.size() + 1, 0);
    //
    //     std::vector<double> vsaContribs(mol.getNumAtoms());
    //     double tmp;
    //     getLabuteAtomContribs(mol, vsaContribs, tmp, true, force);
    //     std::vector<double> logpContribs(mol.getNumAtoms());
    //     std::vector<double> mrContribs(mol.getNumAtoms());
    //     getCrippenAtomContribs(mol, logpContribs, mrContribs, force);
    //
    //     assignContribsToBins(vsaContribs, logpContribs, lbins, res);
    //
    //     return res;
    //   }
    // RDKit✔️✔️:   std::vector<double> calcSlogP_VSA(const ROMol &mol, std::vector<double> *bins,
    // RDKit✔️✔️:                                     bool force) {
    // RDKit✔️✔️:   // FIX: use force value to include caching
    // RDKit✔️✔️:   std::vector<double> lbins;
    // RDKit✔️✔️:   if (!bins) {
    // RDKit✔️✔️:     double blist[11] = {-0.4, -0.2, 0,   0.1, 0.15, 0.2,
    // RDKit✔️✔️:                         0.25, 0.3,  0.4, 0.5, 0.6};
    // RDKit✔️✔️:     lbins.resize(11);
    // RDKit✔️✔️:     std::copy(blist, blist + 11, lbins.begin());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     lbins.resize(bins->size());
    // RDKit✔️✔️:     std::copy(bins->begin(), bins->end(), lbins.begin());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<double> res(lbins.size() + 1, 0);
    // RDKit✔️✔️:   std::vector<double> vsaContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   double tmp;
    // RDKit✔️✔️:   getLabuteAtomContribs(mol, vsaContribs, tmp, true, force);
    // RDKit✔️✔️:   std::vector<double> logpContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   std::vector<double> mrContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   getCrippenAtomContribs(mol, logpContribs, mrContribs, force);
    // RDKit✔️✔️:   assignContribsToBins(vsaContribs, logpContribs, lbins, res);
    // RDKit✔️✔️:   return res;
    // Typed mapping: bins==nullptr -> Option::None selects the frozen
    // DEFAULT_SLOGP_BINS copy; Some copies the caller's slice; the
    // Labute owner is called with include_hydrogens=true (the hydrogen
    // lump is the source's discarded `tmp`); the Crippen logP rows key
    // the bins (MR rows discarded); the accumulation is the V01 owner.
    // The source's `force` also reaches getCrippenAtomContribs; the
    // canonical Crippen entry reads the four-field cache slot in the
    // caller's state (warm hit at exact row lengths; force=true runs
    // the ONE cold kernel) — the same cache arm as crippen_totals.
    let lbins: Vec<f64> = match bins {
        None => DEFAULT_SLOGP_BINS.to_vec(),
        Some(custom) => custom.to_vec(),
    };
    let vsa_contribs = crate::labute::labute_contributions(input, true, force, state)?;
    let logp_contribs = crate::crippen::crippen_contributions(input, force, state, None, None)?;
    assign_contribs_to_bins(&vsa_contribs.atoms, &logp_contribs.logp, &lbins)
}

/// The frozen default SMR_VSA bin boundaries (the source `blist[9]`,
/// MolSurf.cpp:409-410, copied verbatim).
pub(crate) const DEFAULT_SMR_BINS: [f64; 9] = [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0];

/// SMR_VSA: Labute VSA contributions accumulated into Crippen-MR
/// bins (the source `calcSMR_VSA`, MolSurf.cpp:406-429).
///
/// Behavior review: structurally identical to [`slogp_vsa`] with the
/// two source differences — the frozen default is the 9-boundary list
/// (10 outputs) and the binning property is the CRIPPEN MR row per
/// atom (the logP rows are discarded exactly as the source discards
/// `logpContribs`). Same owners, same equality-to-edge placement,
/// same caller sortedness precondition, same force forwarding to BOTH
/// owners' cache arms (Labute slot; the canonical Crippen four-field
/// slot through crippen_contributions).
///
/// Complexity review: identical to [`slogp_vsa`] modulo m=9 — two
/// owner calls at their audited costs, the V01 O(n log m) binning,
/// and the owned lbins/res/contribution copies (the cold Crippen path
/// also stores cache row copies in the caller's state); no map, no
/// sort, no repeated scans.
pub fn smr_vsa(
    input: &DescriptorInput<'_>,
    bins: Option<&[f64]>,
    force: bool,
    state: &mut DescriptorComputedState,
) -> DescriptorResult<Vec<f64>> {
    // RDKit source (MolSurf.cpp:406-429):
    //   std::vector<double> calcSMR_VSA(const ROMol &mol, std::vector<double> *bins,
    //                                   bool force) {
    //     std::vector<double> lbins;
    //     if (!bins) {
    //       double blist[9] = {1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0};
    //       lbins.resize(9);
    //       std::copy(blist, blist + 9, lbins.begin());
    //     } else {
    //       lbins.resize(bins->size());
    //       std::copy(bins->begin(), bins->end(), lbins.begin());
    //     }
    //     std::vector<double> res(lbins.size() + 1, 0);
    //
    //     std::vector<double> vsaContribs(mol.getNumAtoms());
    //     double tmp;
    //     getLabuteAtomContribs(mol, vsaContribs, tmp, true, force);
    //     std::vector<double> logpContribs(mol.getNumAtoms());
    //     std::vector<double> mrContribs(mol.getNumAtoms());
    //     getCrippenAtomContribs(mol, logpContribs, mrContribs, force);
    //
    //     assignContribsToBins(vsaContribs, mrContribs, lbins, res);
    //
    //     return res;
    //   }
    // RDKit✔️✔️:   std::vector<double> calcSMR_VSA(const ROMol &mol, std::vector<double> *bins,
    // RDKit✔️✔️:                                   bool force) {
    // RDKit✔️✔️:   std::vector<double> lbins;
    // RDKit✔️✔️:   if (!bins) {
    // RDKit✔️✔️:     double blist[9] = {1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0};
    // RDKit✔️✔️:     lbins.resize(9);
    // RDKit✔️✔️:     std::copy(blist, blist + 9, lbins.begin());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     lbins.resize(bins->size());
    // RDKit✔️✔️:     std::copy(bins->begin(), bins->end(), lbins.begin());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<double> res(lbins.size() + 1, 0);
    // RDKit✔️✔️:   std::vector<double> vsaContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   double tmp;
    // RDKit✔️✔️:   getLabuteAtomContribs(mol, vsaContribs, tmp, true, force);
    // RDKit✔️✔️:   std::vector<double> logpContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   std::vector<double> mrContribs(mol.getNumAtoms());
    // RDKit✔️✔️:   getCrippenAtomContribs(mol, logpContribs, mrContribs, force);
    // RDKit✔️✔️:   assignContribsToBins(vsaContribs, mrContribs, lbins, res);
    // RDKit✔️✔️:   return res;
    // Typed mapping: identical to slogp_vsa except the frozen 9-boundary
    // default and the binning property — the CRIPPEN MR rows key the
    // bins here (the logP rows are the discarded vector). The Crippen
    // force/cache routing is the same as slogp_vsa.
    let lbins: Vec<f64> = match bins {
        None => DEFAULT_SMR_BINS.to_vec(),
        Some(custom) => custom.to_vec(),
    };
    let vsa_contribs = crate::labute::labute_contributions(input, true, force, state)?;
    let mr_contribs = crate::crippen::crippen_contributions(input, force, state, None, None)?;
    assign_contribs_to_bins(&vsa_contribs.atoms, &mr_contribs.molar_refractivity, &lbins)
}
