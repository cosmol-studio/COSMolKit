//! Chi recomputation arithmetic over borrowed prepared molecule state.
//!
//! Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8 (BSD).
//! These domain entrypoints reproduce the complete source recomputation
//! arithmetic only. Source vector-property cache/force/computed lifecycle
//! belongs to the cached adapters below. The retained recomputation APIs
//! remain available for detached callers without source computed properties.
use crate::{ChiInput, DescriptorError, DescriptorResult};
use cosmolkit_core::{
    GraphPath, PathRepresentation, PathSearchParams, all_paths_of_length, element_info,
    total_hydrogen_count_from_validated,
};
use cosmolkit_model::TopologyBlock;
/// Pinned source Chi0v version.
pub const CHI_0_V_VERSION: &str = "1.2.0";
/// Pinned source Chi1v version.
pub const CHI_1_V_VERSION: &str = "1.2.0";
/// Pinned source Chi2v version.
pub const CHI_2_V_VERSION: &str = "1.2.0";
/// Pinned source Chi3v version.
pub const CHI_3_V_VERSION: &str = "1.2.0";
/// Pinned source Chi4v version.
pub const CHI_4_V_VERSION: &str = "1.2.0";
/// Pinned source ChiNv version.
pub const CHI_N_V_VERSION: &str = "1.2.0";
/// Pinned source Chi0n version.
pub const CHI_0_N_VERSION: &str = "1.2.0";
/// Pinned source Chi1n version.
pub const CHI_1_N_VERSION: &str = "1.2.0";
/// Pinned source Chi2n version.
pub const CHI_2_N_VERSION: &str = "1.2.0";
/// Pinned source Chi3n version.
pub const CHI_3_N_VERSION: &str = "1.2.0";
/// Pinned source Chi4n version.
pub const CHI_4_N_VERSION: &str = "1.2.0";
/// Pinned source ChiNn version.
pub const CHI_N_N_VERSION: &str = "1.2.0";

// Local structural preflight only; supplied valence is ALWAYS borrowed.
// Numeric implicit H is read by the source per-atom branch, never globally.
fn validate_input(input: &ChiInput<'_>, function: &'static str) -> DescriptorResult<()> {
    input
        .topology()
        .validate()
        .map_err(|source| DescriptorError::InvalidTopology { function, source })?;
    crate::prepared_valence(input.topology(), Some(input.valence()), function)?;
    Ok(())
}

// Behavioral scope: recomputation rows at atomic indices, exact source u32
// arithmetic/getter exclusion/IEEE operation order. Cache lines stay ❌❌.
// Complexity: O(V) numeric loop/one Vec, but structural preflight adds
// O(V+E log(E+1)) adjacency rebuild with BTreeMap/Vec allocations, plus
// supplied structural metadata checks, beyond already-valid ROMol.
// No graph/assignment/ring clone, table duplication, chemistry or cache.
fn chi_v_weights_kernel(
    input: &ChiInput<'_>,
    function: &'static str,
) -> DescriptorResult<Vec<f64>> {
    // RDKit❗❌: void hkDeltas(const ROMol &mol, std::vector<double> &deltas, bool force) {
    // RDKit✔️❌:   PRECONDITION(deltas.size() >= mol.getNumAtoms(), "bad vector size");
    // RDKit❌❌:   if (!force && mol.hasProp(common_properties::_connectivityHKDeltas)) {
    // RDKit❌❌:     mol.getProp(common_properties::_connectivityHKDeltas, deltas);
    // RDKit❌❌:     return;
    // RDKit❌❌:   }
    // RDKit✔️❌:   const PeriodicTable *tbl = PeriodicTable::getTable();
    // RDKit✔️❌:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️❌:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️❌:   while (atBegin != atEnd) {
    // RDKit✔️❌:     const Atom *at = mol[*atBegin];
    // RDKit✔️❌:     unsigned int n = at->getAtomicNum();
    // RDKit✔️❌:     if (n <= 1) {
    // RDKit✔️❌:       deltas[at->getIdx()] = 0;
    // RDKit✔️❌:     } else if (n <= 10) {
    // RDKit✔️❌:       deltas[at->getIdx()] = tbl->getNouterElecs(n) - at->getTotalNumHs();
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       deltas[at->getIdx()] =
    // RDKit✔️❌:           double(tbl->getNouterElecs(n) - at->getTotalNumHs()) /
    // RDKit✔️❌:           (n - tbl->getNouterElecs(n) - 1);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (deltas[at->getIdx()] != 0.0) {
    // RDKit✔️❌:       deltas[at->getIdx()] = 1. / sqrt(deltas[at->getIdx()]);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++atBegin;
    // RDKit✔️❌:   }
    // RDKit❌❌:   mol.setProp(common_properties::_connectivityHKDeltas, deltas, true);
    // RDKit✔️❌: }

    // Signature/cache scope above is deliberately not the source cached API.
    // The atom-count output allocation establishes the source sink size.
    validate_input(input, function)?;
    let topology = input.topology();
    let mut rows = vec![0.0; topology.atoms.len()];
    for atom in &topology.atoms {
        let z = u32::from(atom.atomic_number());
        if z <= 1 {
            continue;
        }
        let outer = element_info(atom.element()).outer_electrons as u32;
        let hydrogens =
            total_hydrogen_count_from_validated(topology, input.valence(), atom.id(), false)
                .map_err(|source| DescriptorError::Valence { function, source })?;
        let mut value = f64::from(outer.wrapping_sub(hydrogens));
        if z > 10 {
            value /= f64::from(z.wrapping_sub(outer).wrapping_sub(1));
        }
        if value != 0.0 {
            value = 1.0 / value.sqrt();
        }
        rows[atom.id().index()] = value;
    }
    Ok(rows)
}

// Behavioral scope: every atom including wildcard/H reads the sole table/H
// owners, with unsigned subtraction BEFORE f64 and exact sqrt/reciprocal.
// Complexity: same structural preflight/allocation debt as v, O(V) numeric
// loop/one Vec and no cloning, preparation, second table or state writes.
fn chi_n_weights_kernel(
    input: &ChiInput<'_>,
    function: &'static str,
) -> DescriptorResult<Vec<f64>> {
    // RDKit❗❌: void nVals(const ROMol &mol, std::vector<double> &nVs, bool force) {
    // RDKit✔️❌:   PRECONDITION(nVs.size() >= mol.getNumAtoms(), "bad vector size");
    // RDKit❌❌:   if (!force && mol.hasProp(common_properties::_connectivityNVals)) {
    // RDKit❌❌:     mol.getProp(common_properties::_connectivityNVals, nVs);
    // RDKit❌❌:     return;
    // RDKit❌❌:   }
    // RDKit✔️❌:   const PeriodicTable *tbl = PeriodicTable::getTable();
    // RDKit✔️❌:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️❌:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️❌:   while (atBegin != atEnd) {
    // RDKit✔️❌:     const Atom *at = mol[*atBegin];
    // RDKit✔️❌:     double v = tbl->getNouterElecs(at->getAtomicNum()) - at->getTotalNumHs();
    // RDKit✔️❌:     if (v != 0.0) {
    // RDKit✔️❌:       v = 1. / sqrt(v);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     nVs[at->getIdx()] = v;
    // RDKit✔️❌:     ++atBegin;
    // RDKit✔️❌:   }
    // RDKit❌❌:   mol.setProp(common_properties::_connectivityNVals, nVs, true);
    // RDKit✔️❌: }

    validate_input(input, function)?;
    let topology = input.topology();
    let mut rows = vec![0.0; topology.atoms.len()];
    for atom in &topology.atoms {
        let outer = element_info(atom.element()).outer_electrons as u32;
        let hydrogens =
            total_hydrogen_count_from_validated(topology, input.valence(), atom.id(), false)
                .map_err(|source| DescriptorError::Valence { function, source })?;
        let mut value = f64::from(outer.wrapping_sub(hydrogens));
        if value != 0.0 {
            value = 1.0 / value.sqrt();
        }
        rows[atom.id().index()] = value;
    }
    Ok(rows)
}

// Exact sequential sum, including nonfinite supplied source-weight states.
// O(rows), no allocations or reassociation; comparable source accumulate.
fn chi_zero_from_weights(weights: &[f64]) -> f64 {
    // RDKit✔️✔️:   return std::accumulate(hkDs.begin(), hkDs.end(), 0.0);
    // RDKit✔️✔️:   return std::accumulate(nVs.begin(), nVs.end(), 0.0);
    let mut result = 0.0;
    for &value in weights {
        result += value;
    }
    result
}

// Private callers establish weights[0..numAtoms] before this boundary.
// Preserve stored-edge order/all typed bonds/H bonds and IEEE operations.
// Local preflight adds O(V+E log(E+1)) adjacency rebuild and metadata
// checks with BTreeMap/Vec allocations to the source O(E) numeric loop.
fn chi_one_from_weights(topology: &TopologyBlock, weights: &[f64]) -> DescriptorResult<f64> {
    // RDKit✔️❌:   double res = 0.0;
    // RDKit✔️❌:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️❌:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit✔️❌:   while (firstB != lastB) {
    // RDKit✔️❌:     const Bond *bond = mol[*firstB];
    // RDKit✔️❌:     res += hkDs[bond->getBeginAtomIdx()] * hkDs[bond->getEndAtomIdx()];
    // RDKit✔️❌:     ++firstB;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌:   double res = 0.0;
    // RDKit✔️❌:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️❌:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit✔️❌:   while (firstB != lastB) {
    // RDKit✔️❌:     const Bond *bond = mol[*firstB];
    // RDKit✔️❌:     res += nVs[bond->getBeginAtomIdx()] * nVs[bond->getEndAtomIdx()];
    // RDKit✔️❌:     ++firstB;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;

    topology
        .validate()
        .map_err(|source| DescriptorError::InvalidTopology {
            function: "chi_one_from_weights",
            source,
        })?;
    let mut result = 0.0;
    for bond in &topology.bonds {
        result += weights[bond.begin().index()] * weights[bond.end().index()];
    }
    Ok(result)
}

// The only path engine is core's source-backed Atoms/default search.
// Preserve emitted first-bond-set order, final closure rule, IEEE products
// and u32 order wrap. Weights belong to validated atom rows at all callers.
// Core structural validation and Vec<bool> bond sets add costs relative to
// source ROMol/packed bitsets; no graph/ring clone or new enumeration here.
fn chi_order_from_weights(
    topology: &TopologyBlock,
    weights: &[f64],
    order: u32,
) -> DescriptorResult<f64> {
    // RDKit✔️❌:   PATH_LIST ps = findAllPathsOfLengthN(mol, n + 1, false);
    // RDKit✔️❌:   double res = 0.0;
    // RDKit✔️❌:   for (const auto &p : ps) {
    // RDKit✔️❌:     TEST_ASSERT(p.size() == n + 1);
    // RDKit✔️❌:     double accum = 1.0;
    // RDKit✔️❌:     for (unsigned int i = 0; i < n; ++i) {
    // RDKit✔️❌:       accum *= hkDs[p[i]];
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // only push on the last element if this isn't a ring; this was github 463:
    // RDKit✔️❌:     if (p[n] != p[0]) {
    // RDKit✔️❌:       accum *= hkDs[p[n]];
    // RDKit✔️❌:     }
    // RDKit✔️❌:     res += accum;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌:   PATH_LIST ps = findAllPathsOfLengthN(mol, n + 1, false);
    // RDKit✔️❌:   double res = 0.0;
    // RDKit✔️❌:   for (const auto &p : ps) {
    // RDKit✔️❌:     TEST_ASSERT(p.size() == n + 1);
    // RDKit✔️❌:     double accum = 1.0;
    // RDKit✔️❌:     for (unsigned int i = 0; i < n; ++i) {
    // RDKit✔️❌:       accum *= nVs[p[i]];
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // only push on the last element if this isn't a ring; this was github 463:
    // RDKit✔️❌:     if (p[n] != p[0]) {
    // RDKit✔️❌:       accum *= nVs[p[n]];
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     res += accum;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;

    let expected_rows = order.wrapping_add(1) as usize;
    let params = PathSearchParams {
        representation: PathRepresentation::Atoms,
        ..Default::default()
    };
    let paths = all_paths_of_length(topology, expected_rows, &params).map_err(|source| {
        DescriptorError::Path {
            function: "chi_order_from_weights",
            source,
        }
    })?;
    let mut result = 0.0;
    for path in paths {
        let GraphPath::Atoms(path) = path else {
            return Err(DescriptorError::InvalidConnectivityPath {
                function: "chi_order_from_weights",
                expected_rows,
                actual_rows: None,
            });
        };
        if path.len() != expected_rows || path.is_empty() {
            return Err(DescriptorError::InvalidConnectivityPath {
                function: "chi_order_from_weights",
                expected_rows,
                actual_rows: Some(path.len()),
            });
        }
        let n = order as usize;
        let mut accum = 1.0;
        for index in 0..n {
            accum *= weights[path[index].index()];
        }
        if path[n] != path[0] {
            accum *= weights[path[n].index()];
        }
        result += accum;
    }
    Ok(result)
}

// Keep the outer canonical entry name while retaining the original typed
// cause. This changes no category, value, source ordering or error payload.
fn entry_error(error: DescriptorError, function: &'static str) -> DescriptorError {
    match error {
        DescriptorError::InvalidTopology { source, .. } => {
            DescriptorError::InvalidTopology { function, source }
        }
        DescriptorError::Path { source, .. } => DescriptorError::Path { function, source },
        DescriptorError::InvalidConnectivityPath {
            expected_rows,
            actual_rows,
            ..
        } => DescriptorError::InvalidConnectivityPath {
            function,
            expected_rows,
            actual_rows,
        },
        other => other,
    }
}

fn chi_v_order(input: &ChiInput<'_>, order: u32, function: &'static str) -> DescriptorResult<f64> {
    let weights = chi_v_weights_kernel(input, function)?;
    chi_order_from_weights(input.topology(), &weights, order)
        .map_err(|error| entry_error(error, function))
}
fn chi_n_order(input: &ChiInput<'_>, order: u32, function: &'static str) -> DescriptorResult<f64> {
    let weights = chi_n_weights_kernel(input, function)?;
    chi_order_from_weights(input.topology(), &weights, order)
        .map_err(|error| entry_error(error, function))
}

/// Chi0v arithmetic recomputed from prepared state; no cache or force.
/// Borrows only topology/valence and preserves original typed failures.
/// Shares the private source kernel; extra structural validation costs apply.
pub fn chi_0_v(input: &ChiInput<'_>) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChi0v(const ROMol &mol, bool force) {
    // RDKit✔️❌:   std::vector<double> hkDs(mol.getNumAtoms());
    // RDKit❗❌:   detail::hkDeltas(mol, hkDs, force);
    // RDKit✔️❌:   return std::accumulate(hkDs.begin(), hkDs.end(), 0.0);
    // RDKit✔️❌: };
    let weights = chi_v_weights_kernel(input, "chi_0_v")?;
    Ok(chi_zero_from_weights(&weights))
}

/// Chi1v arithmetic recomputed from prepared state; no cache or force.
/// Borrows only topology/valence and preserves original typed failures.
/// Shares the private source kernel; extra structural validation costs apply.
pub fn chi_1_v(input: &ChiInput<'_>) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChi1v(const ROMol &mol, bool force) {
    // RDKit✔️❌:   std::vector<double> hkDs(mol.getNumAtoms());
    // RDKit❗❌:   detail::hkDeltas(mol, hkDs, force);
    // RDKit✔️❌:
    // RDKit✔️❌:   double res = 0.0;
    // RDKit✔️❌:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️❌:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit✔️❌:   while (firstB != lastB) {
    // RDKit✔️❌:     const Bond *bond = mol[*firstB];
    // RDKit✔️❌:     res += hkDs[bond->getBeginAtomIdx()] * hkDs[bond->getEndAtomIdx()];
    // RDKit✔️❌:     ++firstB;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: };
    let weights = chi_v_weights_kernel(input, "chi_1_v")?;
    chi_one_from_weights(input.topology(), &weights).map_err(|error| entry_error(error, "chi_1_v"))
}

/// Chi2v arithmetic recomputed from prepared state; no cache or force.
/// Borrows only topology/valence and preserves original typed failures.
/// Shares the private source kernel; extra structural validation costs apply.
pub fn chi_2_v(input: &ChiInput<'_>) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChi2v(const ROMol &mol, bool force) {
    // RDKit❗❌:   return calcChiNv(mol, 2, force);
    // RDKit✔️❌: };
    chi_v_order(input, 2, "chi_2_v")
}

/// Chi3v arithmetic recomputed from prepared state; no cache or force.
/// Borrows only topology/valence and preserves original typed failures.
/// Shares the private source kernel; extra structural validation costs apply.
pub fn chi_3_v(input: &ChiInput<'_>) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChi3v(const ROMol &mol, bool force) {
    // RDKit❗❌:   return calcChiNv(mol, 3, force);
    // RDKit✔️❌: };
    chi_v_order(input, 3, "chi_3_v")
}

/// Chi4v arithmetic recomputed from prepared state; no cache or force.
/// Borrows only topology/valence and preserves original typed failures.
/// Shares the private source kernel; extra structural validation costs apply.
pub fn chi_4_v(input: &ChiInput<'_>) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChi4v(const ROMol &mol, bool force) {
    // RDKit❗❌:   return calcChiNv(mol, 4, force);
    // RDKit✔️❌: };
    chi_v_order(input, 4, "chi_4_v")
}

/// General order Chiv recomputation; order0/1 differ from fixed wrappers.
/// Source u32 order wrap and path order are retained, after weight reads.
/// Uses one core enumeration and no reassociated products or input clones.
pub fn chi_n_v(input: &ChiInput<'_>, order: u32) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChiNv(const ROMol &mol, unsigned int n, bool force) {
    // RDKit✔️❌:   std::vector<double> hkDs(mol.getNumAtoms());
    // RDKit❗❌:   detail::hkDeltas(mol, hkDs, force);
    // RDKit✔️❌:   PATH_LIST ps = findAllPathsOfLengthN(mol, n + 1, false);
    // RDKit✔️❌:   double res = 0.0;
    // RDKit✔️❌:   for (const auto &p : ps) {
    // RDKit✔️❌:     TEST_ASSERT(p.size() == n + 1);
    // RDKit✔️❌:     double accum = 1.0;
    // RDKit✔️❌:     for (unsigned int i = 0; i < n; ++i) {
    // RDKit✔️❌:       accum *= hkDs[p[i]];
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // only push on the last element if this isn't a ring; this was github 463:
    // RDKit✔️❌:     if (p[n] != p[0]) {
    // RDKit✔️❌:       accum *= hkDs[p[n]];
    // RDKit✔️❌:     }
    // RDKit✔️❌:     res += accum;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    chi_v_order(input, order, "chi_n_v")
}

/// Chi0n arithmetic recomputed from prepared state; no cache or force.
/// Borrows only topology/valence and preserves original typed failures.
/// Shares the private source kernel; extra structural validation costs apply.
pub fn chi_0_n(input: &ChiInput<'_>) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChi0n(const ROMol &mol, bool force) {
    // RDKit✔️❌:   std::vector<double> nVs(mol.getNumAtoms());
    // RDKit❗❌:   detail::nVals(mol, nVs, force);
    // RDKit✔️❌:   return std::accumulate(nVs.begin(), nVs.end(), 0.0);
    // RDKit✔️❌: };
    let weights = chi_n_weights_kernel(input, "chi_0_n")?;
    Ok(chi_zero_from_weights(&weights))
}

/// Chi1n arithmetic recomputed from prepared state; no cache or force.
/// Borrows only topology/valence and preserves original typed failures.
/// Shares the private source kernel; extra structural validation costs apply.
pub fn chi_1_n(input: &ChiInput<'_>) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChi1n(const ROMol &mol, bool force) {
    // RDKit✔️❌:   std::vector<double> nVs(mol.getNumAtoms());
    // RDKit❗❌:   detail::nVals(mol, nVs, force);
    // RDKit✔️❌:
    // RDKit✔️❌:   double res = 0.0;
    // RDKit✔️❌:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️❌:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit✔️❌:   while (firstB != lastB) {
    // RDKit✔️❌:     const Bond *bond = mol[*firstB];
    // RDKit✔️❌:     res += nVs[bond->getBeginAtomIdx()] * nVs[bond->getEndAtomIdx()];
    // RDKit✔️❌:     ++firstB;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: };
    let weights = chi_n_weights_kernel(input, "chi_1_n")?;
    chi_one_from_weights(input.topology(), &weights).map_err(|error| entry_error(error, "chi_1_n"))
}

/// Chi2n arithmetic recomputed from prepared state; no cache or force.
/// Borrows only topology/valence and preserves original typed failures.
/// Shares the private source kernel; extra structural validation costs apply.
pub fn chi_2_n(input: &ChiInput<'_>) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChi2n(const ROMol &mol, bool force) {
    // RDKit❗❌:   return calcChiNn(mol, 2, force);
    // RDKit✔️❌: };
    chi_n_order(input, 2, "chi_2_n")
}

/// Chi3n arithmetic recomputed from prepared state; no cache or force.
/// Borrows only topology/valence and preserves original typed failures.
/// Shares the private source kernel; extra structural validation costs apply.
pub fn chi_3_n(input: &ChiInput<'_>) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChi3n(const ROMol &mol, bool force) {
    // RDKit❗❌:   return calcChiNn(mol, 3, force);
    // RDKit✔️❌: };
    chi_n_order(input, 3, "chi_3_n")
}

/// Chi4n arithmetic recomputed from prepared state; no cache or force.
/// Borrows only topology/valence and preserves original typed failures.
/// Shares the private source kernel; extra structural validation costs apply.
pub fn chi_4_n(input: &ChiInput<'_>) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChi4n(const ROMol &mol, bool force) {
    // RDKit❗❌:   return calcChiNn(mol, 4, force);
    // RDKit✔️❌: };
    chi_n_order(input, 4, "chi_4_n")
}

/// General order Chin recomputation; order0/1 differ from fixed wrappers.
/// Source u32 order wrap and path order are retained, after weight reads.
/// Uses one core enumeration and no reassociated products or input clones.
pub fn chi_n_n(input: &ChiInput<'_>, order: u32) -> DescriptorResult<f64> {
    // RDKit❗❌: double calcChiNn(const ROMol &mol, unsigned int n, bool force) {
    // RDKit✔️❌:   std::vector<double> nVs(mol.getNumAtoms());
    // RDKit❗❌:   detail::nVals(mol, nVs, force);
    // RDKit✔️❌:   PATH_LIST ps = findAllPathsOfLengthN(mol, n + 1, false);
    // RDKit✔️❌:   double res = 0.0;
    // RDKit✔️❌:   for (const auto &p : ps) {
    // RDKit✔️❌:     TEST_ASSERT(p.size() == n + 1);
    // RDKit✔️❌:     double accum = 1.0;
    // RDKit✔️❌:     for (unsigned int i = 0; i < n; ++i) {
    // RDKit✔️❌:       accum *= nVs[p[i]];
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // only push on the last element if this isn't a ring; this was github 463:
    // RDKit✔️❌:     if (p[n] != p[0]) {
    // RDKit✔️❌:       accum *= nVs[p[n]];
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     res += accum;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    chi_n_order(input, order, "chi_n_n")
}

// Full private unit-test proposal: install only inside the source owner if authorized.
// Fixed supplied-weight arithmetic equals source after vector<double> read.
// These do NOT implement/test a source property-cache adapter or force lifecycle.
#[cfg(test)]
mod chi_weighted_tests {
    use super::{chi_one_from_weights, chi_order_from_weights, chi_zero_from_weights};
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element,
        TopologyBlock,
    };
    fn graph(n: usize, edges: &[(usize, usize)]) -> TopologyBlock {
        let atoms = (0..n)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect::<Vec<_>>();
        let adjacency = AdjacencyList::try_from_topology(n, &bonds).unwrap();
        TopologyBlock {
            atoms,
            bonds,
            adjacency,
            ..Default::default()
        }
    }
    #[test]
    fn chi_weighted_zero_keeps_ieee_and_sequential_order() {
        assert_eq!(
            chi_zero_from_weights(&[9007199254740992.0, 1.0, 1.0]).to_bits(),
            9007199254740992.0f64.to_bits()
        );
        assert_eq!(
            chi_zero_from_weights(&[1.0, 1.0, 9007199254740992.0]).to_bits(),
            9007199254740994.0f64.to_bits()
        );
        assert_eq!(chi_zero_from_weights(&[-0.0]).to_bits(), 0.0f64.to_bits());
        assert_eq!(chi_zero_from_weights(&[-1.0, 4.0]), 3.0);
        assert!(chi_zero_from_weights(&[f64::NAN]).is_nan());
        assert_eq!(chi_zero_from_weights(&[f64::INFINITY]), f64::INFINITY);
        assert_eq!(
            chi_zero_from_weights(&[f64::NEG_INFINITY]),
            f64::NEG_INFINITY
        );
        assert!(chi_zero_from_weights(&[f64::INFINITY, f64::NEG_INFINITY]).is_nan());
    }
    #[test]
    fn chi_weighted_one_keeps_edge_order_and_nonfinite_states() {
        let f = graph(6, &[(0, 1), (2, 3), (4, 5)]);
        let g = graph(6, &[(4, 5), (2, 3), (0, 1)]);
        let w = [9007199254740992.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        assert_eq!(
            chi_one_from_weights(&f, &w).unwrap().to_bits(),
            9007199254740992.0f64.to_bits()
        );
        assert_eq!(
            chi_one_from_weights(&g, &w).unwrap().to_bits(),
            9007199254740994.0f64.to_bits()
        );
        let edge = graph(2, &[(0, 1)]);
        assert!(
            chi_one_from_weights(&edge, &[f64::INFINITY, 0.0])
                .unwrap()
                .is_nan()
        );
        assert_eq!(
            chi_one_from_weights(&edge, &[f64::NEG_INFINITY, 1.0]).unwrap(),
            f64::NEG_INFINITY
        );
        assert_eq!(
            chi_one_from_weights(&edge, &[-0.0, 1.0]).unwrap().to_bits(),
            0.0f64.to_bits()
        );
        assert_eq!(chi_one_from_weights(&edge, &[-1.0, 4.0]).unwrap(), -4.0);
    }
    #[test]
    fn chi_weighted_paths_keep_matrix_order_and_final_closure_rule() {
        let f = graph(6, &[(4, 5), (2, 3), (0, 1)]);
        let before = f.clone();
        assert_eq!(
            chi_order_from_weights(&f, &[9007199254740992.0, 1.0, 1.0, 1.0, 1.0, 1.0], 1)
                .unwrap()
                .to_bits(),
            9007199254740992.0f64.to_bits()
        );
        assert_eq!(f, before);
        let triangle = graph(3, &[(0, 1), (1, 2), (2, 0)]);
        assert_eq!(
            chi_order_from_weights(&triangle, &[2.0, 3.0, 5.0], 2).unwrap(),
            90.0
        );
        assert_eq!(
            chi_order_from_weights(&triangle, &[2.0, 3.0, 5.0], 3).unwrap(),
            30.0
        );
        let tail = graph(4, &[(0, 1), (1, 2), (2, 3), (3, 1)]);
        assert_eq!(
            chi_order_from_weights(&tail, &[1.0, 2.0, 3.0, 4.0], 4).unwrap(),
            48.0
        );
        assert!(
            chi_order_from_weights(&triangle, &[f64::NAN, 3.0, 5.0], 3)
                .unwrap()
                .is_nan()
        );
    }
    #[test]
    fn chi_weighted_generic_zero_and_wrap_skip_all_row_reads() {
        let f = graph(1, &[]);
        assert_eq!(chi_order_from_weights(&f, &[f64::NAN], 0).unwrap(), 1.0);
        assert_eq!(
            chi_order_from_weights(&f, &[f64::NAN], u32::MAX)
                .unwrap()
                .to_bits(),
            0.0f64.to_bits()
        );
    }
}

fn cached_v_weights(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
    function: &'static str,
) -> DescriptorResult<Vec<f64>> {
    // RDKit source (ConnectivityDescriptors.cpp detail::hkDeltas):
    // RDKit❗✔️: if (!force && mol.hasProp(common_properties::_connectivityHKDeltas)) {
    // RDKit❗✔️:   mol.getProp(common_properties::_connectivityHKDeltas, deltas);
    // RDKit❗✔️:   return;
    // RDKit❗✔️: }
    // The source copies cached rows to the caller's output vector before any
    // chemistry work; preserve this O(V) copy and guard ordering exactly.
    if !force {
        if let Some(rows) = &state.chi_v_weights {
            return Ok(rows.clone());
        }
    }
    let rows = chi_v_weights_kernel(input, function)?;
    // RDKit❗✔️: mol.setProp(common_properties::_connectivityHKDeltas, deltas, true);
    // Cold numeric work stays in the one retained kernel, whose additional
    // detached structural-validation cost is documented independently.
    state.chi_v_weights = Some(rows.clone());
    Ok(rows)
}

/// Source cached chi_0_v; reuses the existing weight and arithmetic owners.
pub fn chi_0_v_with_state(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::hkDeltas(mol, hkDs, force);
    let weights = cached_v_weights(input, force, state, "chi_0_v")?;
    Ok(chi_zero_from_weights(&weights))
}

/// Source cached chi_1_v; reuses the existing weight and arithmetic owners.
pub fn chi_1_v_with_state(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::hkDeltas(mol, hkDs, force);
    let weights = cached_v_weights(input, force, state, "chi_1_v")?;
    chi_one_from_weights(input.topology(), &weights).map_err(|error| entry_error(error, "chi_1_v"))
}

/// Source cached chi_2_v; reuses the existing weight and arithmetic owners.
pub fn chi_2_v_with_state(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::hkDeltas(mol, hkDs, force);
    let weights = cached_v_weights(input, force, state, "chi_2_v")?;
    chi_order_from_weights(input.topology(), &weights, 2)
        .map_err(|error| entry_error(error, "chi_2_v"))
}

/// Source cached chi_3_v; reuses the existing weight and arithmetic owners.
pub fn chi_3_v_with_state(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::hkDeltas(mol, hkDs, force);
    let weights = cached_v_weights(input, force, state, "chi_3_v")?;
    chi_order_from_weights(input.topology(), &weights, 3)
        .map_err(|error| entry_error(error, "chi_3_v"))
}

/// Source cached chi_4_v; reuses the existing weight and arithmetic owners.
pub fn chi_4_v_with_state(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::hkDeltas(mol, hkDs, force);
    let weights = cached_v_weights(input, force, state, "chi_4_v")?;
    chi_order_from_weights(input.topology(), &weights, 4)
        .map_err(|error| entry_error(error, "chi_4_v"))
}

/// Source cached chi_n_v; reuses the existing weight and arithmetic owners.
pub fn chi_n_v_with_state(
    input: &ChiInput<'_>,
    order: u32,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::hkDeltas(mol, hkDs, force);
    let weights = cached_v_weights(input, force, state, "chi_n_v")?;
    chi_order_from_weights(input.topology(), &weights, order)
        .map_err(|error| entry_error(error, "chi_n_v"))
}

fn cached_n_weights(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
    function: &'static str,
) -> DescriptorResult<Vec<f64>> {
    // RDKit source (ConnectivityDescriptors.cpp detail::nVals):
    // RDKit❗✔️: if (!force && mol.hasProp(common_properties::_connectivityNVals)) {
    // RDKit❗✔️:   mol.getProp(common_properties::_connectivityNVals, nVs);
    // RDKit❗✔️:   return;
    // RDKit❗✔️: }
    // The source copies cached rows to the caller's output vector before any
    // chemistry work; preserve this O(V) copy and guard ordering exactly.
    if !force {
        if let Some(rows) = &state.chi_n_weights {
            return Ok(rows.clone());
        }
    }
    let rows = chi_n_weights_kernel(input, function)?;
    // RDKit❗✔️: mol.setProp(common_properties::_connectivityNVals, nVs, true);
    // Cold numeric work stays in the one retained kernel, whose additional
    // detached structural-validation cost is documented independently.
    state.chi_n_weights = Some(rows.clone());
    Ok(rows)
}

/// Source cached chi_0_n; reuses the existing weight and arithmetic owners.
pub fn chi_0_n_with_state(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::nVals(mol, nVs, force);
    let weights = cached_n_weights(input, force, state, "chi_0_n")?;
    Ok(chi_zero_from_weights(&weights))
}

/// Source cached chi_1_n; reuses the existing weight and arithmetic owners.
pub fn chi_1_n_with_state(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::nVals(mol, nVs, force);
    let weights = cached_n_weights(input, force, state, "chi_1_n")?;
    chi_one_from_weights(input.topology(), &weights).map_err(|error| entry_error(error, "chi_1_n"))
}

/// Source cached chi_2_n; reuses the existing weight and arithmetic owners.
pub fn chi_2_n_with_state(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::nVals(mol, nVs, force);
    let weights = cached_n_weights(input, force, state, "chi_2_n")?;
    chi_order_from_weights(input.topology(), &weights, 2)
        .map_err(|error| entry_error(error, "chi_2_n"))
}

/// Source cached chi_3_n; reuses the existing weight and arithmetic owners.
pub fn chi_3_n_with_state(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::nVals(mol, nVs, force);
    let weights = cached_n_weights(input, force, state, "chi_3_n")?;
    chi_order_from_weights(input.topology(), &weights, 3)
        .map_err(|error| entry_error(error, "chi_3_n"))
}

/// Source cached chi_4_n; reuses the existing weight and arithmetic owners.
pub fn chi_4_n_with_state(
    input: &ChiInput<'_>,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::nVals(mol, nVs, force);
    let weights = cached_n_weights(input, force, state, "chi_4_n")?;
    chi_order_from_weights(input.topology(), &weights, 4)
        .map_err(|error| entry_error(error, "chi_4_n"))
}

/// Source cached chi_n_n; reuses the existing weight and arithmetic owners.
pub fn chi_n_n_with_state(
    input: &ChiInput<'_>,
    order: u32,
    force: bool,
    state: &mut crate::DescriptorComputedState,
) -> DescriptorResult<f64> {
    // RDKit❗✔️: detail::nVals(mol, nVs, force);
    let weights = cached_n_weights(input, force, state, "chi_n_n")?;
    chi_order_from_weights(input.topology(), &weights, order)
        .map_err(|error| entry_error(error, "chi_n_n"))
}

#[cfg(test)]
mod chi_source_cache_tests {
    use super::*;
    use cosmolkit_core::ValenceAssignment;

    #[test]
    fn chi_cache_warm_force_failure_clear_and_clone() {
        // Source returns supplied computed rows before reading atom valence;
        // force instead evaluates its original kernel. Cached state survives
        // failed cold evaluation; detached clone and clear are independent.
        let topology = TopologyBlock::default();
        let malformed = ValenceAssignment {
            explicit_valence: vec![1],
            implicit_hydrogens: vec![1],
        };
        let input = ChiInput::new(&topology, &malformed);
        let mut original = crate::DescriptorComputedState::default();
        original.chi_v_weights = Some(vec![2.0, 3.0]);
        original.chi_n_weights = Some(vec![4.0, 5.0]);
        let baseline = original.clone();
        assert_eq!(
            chi_0_v_with_state(&input, false, &mut original).unwrap(),
            5.0
        );
        assert_eq!(
            chi_0_n_with_state(&input, false, &mut original).unwrap(),
            9.0
        );
        assert!(chi_0_v_with_state(&input, true, &mut original).is_err());
        assert!(chi_0_n_with_state(&input, true, &mut original).is_err());
        assert_eq!(original, baseline);
        let mut peer = original.clone();
        peer.clear();
        assert!(chi_0_v_with_state(&input, false, &mut peer).is_err());
        assert!(chi_0_n_with_state(&input, false, &mut peer).is_err());
        assert_eq!(original, baseline);
        let valid = ValenceAssignment {
            explicit_valence: vec![],
            implicit_hydrogens: vec![],
        };
        let empty = ChiInput::new(&topology, &valid);
        assert_eq!(
            chi_0_v_with_state(&empty, true, &mut peer)
                .unwrap()
                .to_bits(),
            0_f64.to_bits()
        );
        assert_eq!(peer.chi_v_weights, Some(vec![]));
        assert_eq!(original, baseline);
    }
}
