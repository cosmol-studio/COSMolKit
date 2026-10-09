//! Source-backed Topological Torsion fingerprints over prepared detached values.
use crate::generator::{
    FingerprintArguments, FingerprintEnvironment, FingerprintFuncArguments,
    accumulate_sparse_counts_into, project_count_fingerprint, project_fingerprint,
    project_sparse_fingerprint, with_fingerprint_environment_inputs_and_output,
};
use crate::metadata::{
    FingerprintJsonError, bool_or, common_arguments_from_json, common_arguments_json,
    common_arguments_string, parse_object, u32_or,
};
use crate::{
    AtomPairAtomInvariantsGenerator, AtomPairError, AtomPairPreparedInput, Fingerprint,
    FingerprintAdditionalOutput, FingerprintError, MorganError, SparseBitFingerprint,
    SparseCountFingerprint, SparseCountFingerprint32, topological_torsion_code,
    topological_torsion_hash,
};
use cosmolkit_core::{
    GraphPath, PathError, PathRepresentation, PathSearchParams, all_paths_of_length,
};
use cosmolkit_model::AtomId;
use serde_json::Value;
use std::fmt;

#[derive(Debug)]
pub enum TopologicalTorsionError {
    Fingerprint(FingerprintError),
    AtomInvariants(AtomPairError),
    Path(PathError),
    Preparation(MorganError),
    Json(FingerprintJsonError),
    StatePoisoned,
    ThreadCount(cosmolkit_core::ThreadCountError),
    ThreadSpawn(std::io::Error),
    WorkerPanic,
    WorkerProtocol,
    OutputAllocation(std::collections::TryReserveError),
}
impl fmt::Display for TopologicalTorsionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Fingerprint(e) => e.fmt(f),
            Self::AtomInvariants(e) => e.fmt(f),
            Self::Path(e) => e.fmt(f),
            Self::Preparation(e) => e.fmt(f),
            Self::Json(e) => e.fmt(f),
            Self::ThreadCount(e) => e.fmt(f),
            Self::ThreadSpawn(e) => {
                write!(f, "Topological Torsion bulk worker creation failed: {e}")
            }
            Self::WorkerPanic => f.write_str("Topological Torsion bulk worker panicked"),
            Self::WorkerProtocol => f.write_str(
                "Topological Torsion bulk worker returned an incomplete result sequence",
            ),
            Self::OutputAllocation(e) => {
                write!(f, "Topological Torsion IDs output allocation failed: {e}")
            }
            Self::StatePoisoned => {
                f.write_str("Topological Torsion generator state lock was poisoned")
            }
        }
    }
}
impl std::error::Error for TopologicalTorsionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fingerprint(e) => Some(e),
            Self::AtomInvariants(e) => Some(e),
            Self::Path(e) => Some(e),
            Self::Preparation(e) => Some(e),
            Self::Json(e) => Some(e),
            Self::ThreadCount(e) => Some(e),
            Self::ThreadSpawn(e) => Some(e),
            Self::OutputAllocation(e) => Some(e),
            Self::StatePoisoned | Self::WorkerPanic | Self::WorkerProtocol => None,
        }
    }
}
impl From<FingerprintError> for TopologicalTorsionError {
    fn from(e: FingerprintError) -> Self {
        Self::Fingerprint(e)
    }
}
impl From<AtomPairError> for TopologicalTorsionError {
    fn from(e: AtomPairError) -> Self {
        Self::AtomInvariants(e)
    }
}
impl From<PathError> for TopologicalTorsionError {
    fn from(e: PathError) -> Self {
        Self::Path(e)
    }
}
impl From<MorganError> for TopologicalTorsionError {
    fn from(e: MorganError) -> Self {
        Self::Preparation(e)
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologicalTorsionParams {
    pub torsion_atom_count: u32,
    pub only_shortest_paths: bool,
    pub include_chirality: bool,
    pub count_simulation: bool,
    pub count_bounds: Vec<u32>,
    pub fp_size: u32,
    pub bits_per_feature: u32,
}
impl Default for TopologicalTorsionParams {
    fn default() -> Self {
        Self {
            torsion_atom_count: 4,
            only_shortest_paths: false,
            include_chirality: false,
            count_simulation: true,
            count_bounds: vec![1, 2, 4, 8],
            fp_size: 2048,
            bits_per_feature: 1,
        }
    }
}
impl TopologicalTorsionParams {
    pub fn new(
        include_chirality: bool,
        torsion_atom_count: u32,
        count_simulation: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
    ) -> Result<Self, FingerprintError> {
        let params = Self {
            include_chirality,
            torsion_atom_count,
            count_simulation,
            count_bounds,
            fp_size,
            ..Default::default()
        };
        params.common()?;
        Ok(params)
    }

    #[must_use]
    pub fn info_string(&self) -> String {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 47-51
        // RDKit❗✔️: std::string TopologicalTorsionArguments::infoString() const {
        // RDKit❗✔️:   return "TopologicalTorsionArguments torsionAtomCount=" +
        // RDKit❗✔️:          std::to_string(d_torsionAtomCount) +
        // RDKit❗✔️:          " onlyShortestPaths=" + std::to_string(df_onlyShortestPaths);
        // RDKit❗✔️: };
        format!(
            "TopologicalTorsionArguments torsionAtomCount={} onlyShortestPaths={}",
            self.torsion_atom_count, self.only_shortest_paths as u8
        )
    }

    #[must_use]
    pub fn to_json(&self) -> String {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 52-58
        // RDKit❗✔️: void TopologicalTorsionArguments::toJSON(
        // RDKit❗✔️:     boost::property_tree::ptree &pt) const {
        // RDKit❗✔️:   pt.put("type", "TopologicalTorsionArguments");
        // RDKit❗✔️:   pt.put("torsionAtomCount", d_torsionAtomCount);
        // RDKit❗✔️:   pt.put("onlyShortestPaths", df_onlyShortestPaths);
        // RDKit❗✔️:   FingerprintArguments::toJSON(pt);
        // RDKit❗✔️: }
        let common = common_arguments_json(
            self.count_simulation,
            self.fp_size,
            self.bits_per_feature,
            self.include_chirality,
            &self.count_bounds,
        );
        let common_body = &common[1..common.len().saturating_sub(1)];
        let mut json = format!(
            "{{\"type\":\"TopologicalTorsionArguments\",\"torsionAtomCount\":\"{}\",\"onlyShortestPaths\":\"{}\"",
            self.torsion_atom_count, self.only_shortest_paths
        );
        if !common_body.is_empty() {
            json.push(',');
            json.push_str(common_body);
        }
        json.push('}');
        json
    }

    pub fn from_json(&mut self, json: &str) -> Result<(), FingerprintJsonError> {
        if json.trim().is_empty() {
            return Ok(());
        }
        self.from_json_value(&parse_object(json)?)
    }
    pub(crate) fn common_info_string(&self) -> String {
        common_arguments_string(
            self.count_simulation,
            self.fp_size,
            self.bits_per_feature,
            self.include_chirality,
        )
    }
    pub(crate) fn from_json_value(
        &mut self,
        value: &crate::metadata::SourceNode,
    ) -> Result<(), FingerprintJsonError> {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 59-65
        // RDKit❗✔️: void TopologicalTorsionArguments::fromJSON(
        // RDKit❗✔️:     const boost::property_tree::ptree &pt) {
        // RDKit❗✔️:   d_torsionAtomCount = pt.get<uint32_t>("torsionAtomCount", d_torsionAtomCount);
        // RDKit❗✔️:   df_onlyShortestPaths =
        // RDKit❗✔️:       pt.get<bool>("onlyShortestPaths", df_onlyShortestPaths);
        // RDKit❗✔️:   FingerprintArguments::fromJSON(pt);
        // RDKit❗✔️: }
        self.torsion_atom_count = u32_or(value, "torsionAtomCount", self.torsion_atom_count);
        self.only_shortest_paths = bool_or(value, "onlyShortestPaths", self.only_shortest_paths);
        common_arguments_from_json(
            value,
            &mut self.count_simulation,
            &mut self.fp_size,
            &mut self.bits_per_feature,
            &mut self.include_chirality,
            &mut self.count_bounds,
        )
    }
    fn common(&self) -> Result<FingerprintArguments, FingerprintError> {
        // RDKit❗✔️: TopologicalTorsionArguments::TopologicalTorsionArguments(
        // RDKit❗✔️:     const bool includeChirality, const uint32_t torsionAtomCount,
        // RDKit❗✔️:     const bool countSimulation, const std::vector<std::uint32_t> countBounds,
        // RDKit❗✔️:     const std::uint32_t fpSize)
        // RDKit❗✔️:     : FingerprintArguments(countSimulation, countBounds, fpSize, 1,
        // RDKit❗✔️:                            includeChirality),
        // RDKit❗✔️:       d_torsionAtomCount(torsionAtomCount) {}
        FingerprintArguments::new(
            self.count_simulation,
            self.count_bounds.clone(),
            self.fp_size,
            self.bits_per_feature,
            self.include_chirality,
        )
    }
    fn result_size(&self) -> Result<u64, FingerprintError> {
        // RDKit source: TopologicalTorsionGenerator.cpp lines 33-45
        // RDKit❗✔️: template <typename OutputType>
        // RDKit❗✔️: OutputType TopologicalTorsionEnvGenerator<OutputType>::getResultSize() const {
        // RDKit❗✔️:   OutputType result = 1;
        // RDKit❗✔️:   return (result << ((
        // RDKit❗✔️:               dynamic_cast<const TopologicalTorsionArguments *>(
        // RDKit❗✔️:                   this->dp_fingerprintArguments)
        // RDKit❗✔️:                   ->d_torsionAtomCount *
        // RDKit❗✔️:               (codeSize + (dynamic_cast<const TopologicalTorsionArguments *>(
        // RDKit❗✔️:                                this->dp_fingerprintArguments)
        // RDKit❗✔️:                                    ->df_includeChirality
        // RDKit❗✔️:                                ? numChiralBits
        // RDKit❗✔️:                                : 0)))));
        // RDKit❗✔️: };
        let bits_per_atom = 9 + if self.include_chirality { 2 } else { 0 };
        let shift = self.torsion_atom_count.checked_mul(bits_per_atom).ok_or(
            FingerprintError::InvalidArguments {
                reason: "topological torsion result-size width overflow",
            },
        )?;
        if shift >= u64::BITS {
            return Err(FingerprintError::InvalidArguments {
                reason: "topological torsion result-size shift must be less than 64 bits",
            });
        }
        Ok(1_u64 << shift)
    }
}
#[derive(Debug, Clone, Copy, Default)]
enum TorsionCodeMode {
    #[default]
    Modern,
    LegacyUnfolded,
    // Custom legacy codes use modern modulo, but the compatibility length.
    LegacyCustom,
}
#[derive(Debug, Clone, Copy)]
pub struct TopologicalTorsionCall<'a> {
    pub from_atoms: Option<&'a [u32]>,
    pub ignore_atoms: Option<&'a [u32]>,
    pub custom_atom_invariants: Option<&'a [u32]>,
    pub custom_bond_invariants: Option<&'a [u32]>,
    pub conformer_id: i32,
    pub atom_invariants_generator: Option<AtomPairAtomInvariantsGenerator>,
}
impl Default for TopologicalTorsionCall<'_> {
    fn default() -> Self {
        Self {
            from_atoms: None,
            ignore_atoms: None,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            conformer_id: -1,
            atom_invariants_generator: None,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
struct TorsionEnvironment {
    bit_id: u64,
    path: Vec<AtomId>,
}
fn environments(
    input: &AtomPairPreparedInput<'_>,
    arguments: &TopologicalTorsionParams,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: &[u32],
    hash_results: bool,
    mode: TorsionCodeMode,
) -> Result<Vec<TorsionEnvironment>, TopologicalTorsionError> {
    // RDKit source: TopologicalTorsionGenerator.cpp lines 101-194
    // RDKit❗✔️: template <typename OutputType>
    // RDKit❗✔️: std::vector<AtomEnvironment<OutputType> *>
    // RDKit❗✔️: TopologicalTorsionEnvGenerator<OutputType>::getEnvironments(
    // RDKit❗✔️:     const ROMol &mol, FingerprintArguments *arguments,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *fromAtoms,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *ignoreAtoms,
    // RDKit❗✔️:     const int,                 // confId
    // RDKit❗✔️:     const AdditionalOutput *,  // additionalOutput
    // RDKit❗✔️:     const std::vector<std::uint32_t> *atomInvariants,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *,  // bondInvariants
    // RDKit❗✔️:     const bool hashResults) const {
    // RDKit❗✔️:   auto *topologicalTorsionArguments =
    // RDKit❗✔️:       dynamic_cast<TopologicalTorsionArguments *>(arguments);
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<AtomEnvironment<OutputType> *> result;
    // RDKit❗✔️:
    // RDKit❗✔️:   boost::dynamic_bitset<> *fromAtomsBV = nullptr;
    // RDKit❗✔️:   if (fromAtoms) {
    // RDKit❗✔️:     fromAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
    // RDKit❗✔️:     for (auto fAt : *fromAtoms) {
    // RDKit❗✔️:       fromAtomsBV->set(fAt);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   boost::dynamic_bitset<> *ignoreAtomsBV = nullptr;
    // RDKit❗✔️:   if (ignoreAtoms) {
    // RDKit❗✔️:     ignoreAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
    // RDKit❗✔️:     for (auto fAt : *ignoreAtoms) {
    // RDKit❗✔️:       ignoreAtomsBV->set(fAt);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   boost::dynamic_bitset<> pAtoms(mol.getNumAtoms());
    // RDKit❗✔️:   bool useBonds = false;
    // RDKit❗✔️:   bool useHs = false;
    // RDKit❗✔️:   int rootedAtAtom = -1;
    // RDKit❗✔️:   PATH_LIST paths = findAllPathsOfLengthN(
    // RDKit❗✔️:       mol, topologicalTorsionArguments->d_torsionAtomCount, useBonds, useHs,
    // RDKit❗✔️:       rootedAtAtom, topologicalTorsionArguments->df_onlyShortestPaths);
    // RDKit❗✔️:   for (PATH_LIST::const_iterator pathIt = paths.begin(); pathIt != paths.end();
    // RDKit❗✔️:        ++pathIt) {
    // RDKit❗✔️:     bool keepIt = true;
    // RDKit❗✔️:     if (fromAtomsBV) {
    // RDKit❗✔️:       keepIt = false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     std::vector<std::uint32_t> pathCodes;
    // RDKit❗✔️:     const PATH_TYPE &path = *pathIt;
    // RDKit❗✔️:     if (fromAtomsBV) {
    // RDKit❗✔️:       if (fromAtomsBV->test(static_cast<std::uint32_t>(path.front())) ||
    // RDKit❗✔️:           fromAtomsBV->test(static_cast<std::uint32_t>(path.back()))) {
    // RDKit❗✔️:         keepIt = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (keepIt && ignoreAtomsBV) {
    // RDKit❗✔️:       for (int pElem : path) {
    // RDKit❗✔️:         if (ignoreAtomsBV->test(pElem)) {
    // RDKit❗✔️:           keepIt = false;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (keepIt) {
    // RDKit❗✔️:       pAtoms.reset();
    // RDKit❗✔️:       for (auto pIt = path.begin(); pIt < path.end(); ++pIt) {
    // RDKit❗✔️:         // look for a cycle that doesn't start at the first atom
    // RDKit❗✔️:         // we can't effectively canonicalize these at the moment
    // RDKit❗✔️:         // (was github #811)
    // RDKit❗✔️:         if (pIt != path.begin() && *pIt != *(path.begin()) && pAtoms[*pIt]) {
    // RDKit❗✔️:           pathCodes.clear();
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         pAtoms.set(*pIt);
    // RDKit❗✔️:         unsigned int code = (*atomInvariants)[*pIt] % ((1 << codeSize) - 1) + 1;
    // RDKit❗✔️:         // subtract off the branching number:
    // RDKit❗✔️:         if (pIt != path.begin() && pIt + 1 != path.end()) {
    // RDKit❗✔️:           --code;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         pathCodes.push_back(code);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (pathCodes.size()) {
    // RDKit❗✔️:         OutputType code;
    // RDKit❗✔️:         if (hashResults) {
    // RDKit❗✔️:           code = getTopologicalTorsionHash(pathCodes);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           code = getTopologicalTorsionCode(
    // RDKit❗✔️:               pathCodes, topologicalTorsionArguments->df_includeChirality);
    // RDKit❗✔️:         }
    // RDKit❗✔️:         result.push_back(new TopologicalTorsionAtomEnv<OutputType>(code, path));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   delete fromAtomsBV;
    // RDKit❗✔️:   delete ignoreAtomsBV;
    // RDKit❗✔️:
    // RDKit❗✔️:   return result;
    // RDKit❗✔️: };
    let atom_count = input.topology.atoms.len();
    if from_atoms
        .into_iter()
        .chain(ignore_atoms)
        .flatten()
        .any(|&atom_id| atom_id as usize >= atom_count)
    {
        return Err(FingerprintError::InvalidArguments {
            reason: "atom selection contains an atom index outside the molecule",
        }
        .into());
    }

    let mut from_atoms_bits = from_atoms.map(|_| vec![false; atom_count]);
    if let (Some(bits), Some(atom_ids)) = (from_atoms_bits.as_mut(), from_atoms) {
        for &atom_id in atom_ids {
            bits[atom_id as usize] = true;
        }
    }
    let mut ignore_atoms_bits = ignore_atoms.map(|_| vec![false; atom_count]);
    if let (Some(bits), Some(atom_ids)) = (ignore_atoms_bits.as_mut(), ignore_atoms) {
        for &atom_id in atom_ids {
            bits[atom_id as usize] = true;
        }
    }

    let paths = all_paths_of_length(
        input.topology,
        arguments.torsion_atom_count as usize,
        &PathSearchParams {
            ignore_atoms: None,
            representation: PathRepresentation::Atoms,
            use_hydrogens: false,
            rooted_at_atom: None,
            only_shortest_paths: arguments.only_shortest_paths,
        },
    )?;
    let mut result = Vec::new();
    let mut path_atoms = vec![false; atom_count];
    let code_modulus = (1_u32 << 9) - 1;
    for graph_path in paths {
        let GraphPath::Atoms(path) = graph_path else {
            unreachable!("atom-path selector returned bond path")
        };
        let mut keep = from_atoms_bits.is_none();
        if let Some(bits) = from_atoms_bits.as_ref() {
            let (Some(&first), Some(&last)) = (path.first(), path.last()) else {
                return Err(FingerprintError::InvalidArguments {
                    reason: "enumerated topological torsion path is empty",
                }
                .into());
            };
            keep = bits[first.index()] || bits[last.index()];
        }
        if keep
            && ignore_atoms_bits
                .as_ref()
                .is_some_and(|bits| path.iter().any(|&atom_id| bits[atom_id.index()]))
        {
            keep = false;
        }
        if !keep {
            continue;
        }

        path_atoms.fill(false);
        let mut path_codes = Vec::new();
        for (position, &atom) in path.iter().enumerate() {
            let atom_id = atom.index();
            if position != 0 && atom != path[0] && path_atoms[atom_id] {
                path_codes.clear();
                break;
            }
            path_atoms[atom_id] = true;
            // The modern generator and deprecated unfolded API share this
            // traversal but intentionally use different source formulas.
            // RDKit❗✔️: unsigned int code = (*atomInvariants)[*pIt] % ((1 << codeSize) - 1) + 1;
            // RDKit❗✔️: unsigned int code = atomCodes[*pIt] - 1;
            // Source indexes only paths surviving selections and cycle checks.
            // A short vector on a molecule with no selected torsion is a valid
            // no-op. An actually missing indexed value is kept as a typed
            // bounds failure at this access; no fabricated invariant is used.
            let invariant =
                *atom_invariants
                    .get(atom_id)
                    .ok_or(FingerprintError::InvalidArguments {
                        reason: "bad atom invariants size",
                    })?;
            let mut code = match mode {
                TorsionCodeMode::Modern | TorsionCodeMode::LegacyCustom => {
                    invariant % code_modulus + 1
                }
                TorsionCodeMode::LegacyUnfolded => invariant.wrapping_add(1),
            };
            if position != 0 && position + 1 != path.len() {
                code = code.wrapping_sub(1);
            }
            path_codes.push(code);
        }
        if path_codes.is_empty() {
            continue;
        }
        let bit_id = if hash_results {
            u64::from(topological_torsion_hash(&path_codes)?)
        } else {
            topological_torsion_code(&path_codes, arguments.include_chirality)?
        };
        result.push(TorsionEnvironment { bit_id, path });
    }
    Ok(result)
}
impl FingerprintEnvironment<TopologicalTorsionError> for TorsionEnvironment {
    type Output = u64;
    type State = ();
    fn bit_id(
        &self,
        _args: &FingerprintArguments,
        _atoms: &[u32],
        _bonds: &[u32],
        _output: Option<&mut FingerprintAdditionalOutput>,
        _hashed: bool,
        _fp_size: u64,
    ) -> Result<u64, TopologicalTorsionError> {
        // RDKit❗✔️:   return d_bitId;
        Ok(self.bit_id)
    }
    fn update_output(
        &self,
        output: &mut FingerprintAdditionalOutput,
        bit: u64,
        _state: &mut (),
    ) -> Result<(), TopologicalTorsionError> {
        // RDKit❗✔️:   PRECONDITION(additionalOutput, "bad output pointer");
        // RDKit❗✔️:   if (additionalOutput->atomToBits || additionalOutput->atomCounts) {
        // RDKit❗✔️:     for (auto aid : d_atomPath) {
        // RDKit❗✔️:       if (additionalOutput->atomToBits) {
        // RDKit❗✔️:         additionalOutput->atomToBits->at(aid).push_back(bitId);
        // RDKit❗✔️:       }
        // RDKit❗✔️:       if (additionalOutput->atomCounts) {
        // RDKit❗✔️:         additionalOutput->atomCounts->at(aid)++;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        if output.atom_to_bits.is_some() || output.atom_counts.is_some() {
            for atom in &self.path {
                if let Some(rows) = output.atom_to_bits.as_mut() {
                    rows[atom.index()].push(bit);
                }
                if let Some(rows) = output.atom_counts.as_mut() {
                    rows[atom.index()] = rows[atom.index()].wrapping_add(1);
                }
            }
        }
        // RDKit❗✔️:   if (additionalOutput->bitPaths) {
        // RDKit❗✔️:     (*additionalOutput->bitPaths)[bitId].push_back(d_atomPath);
        // RDKit❗✔️:   }
        if let Some(rows) = output.bit_paths.as_mut() {
            rows.entry(bit)
                .or_default()
                .push(self.path.iter().map(|a| a.index() as i32).collect());
        }
        // RDKit❗✔️:   if (additionalOutput->atomsPerBit) {
        // RDKit❗✔️:     (*additionalOutput->atomsPerBit)[bitId].push_back(d_atomPath);
        // RDKit❗✔️:   }
        if let Some(rows) = output.atoms_per_bit.as_mut() {
            rows.entry(bit)
                .or_default()
                .push(self.path.iter().map(|a| a.index() as i32).collect());
        }
        Ok(())
    }
}
fn count_helper(
    input: &AtomPairPreparedInput<'_>,
    params: &TopologicalTorsionParams,
    common: &FingerprintArguments,
    call: &TopologicalTorsionCall<'_>,
    fp_size: u64,
    output: Option<&mut FingerprintAdditionalOutput>,
    mode: TorsionCodeMode,
) -> Result<SparseCountFingerprint, TopologicalTorsionError> {
    configured_count_helper(input, params, common, call, fp_size, output, mode, None)
}
fn configured_count_helper(
    input: &AtomPairPreparedInput<'_>,
    params: &TopologicalTorsionParams,
    common: &FingerprintArguments,
    call: &TopologicalTorsionCall<'_>,
    fp_size: u64,
    output: Option<&mut FingerprintAdditionalOutput>,
    mode: TorsionCodeMode,
    configured_invariants: Option<Option<AtomPairAtomInvariantsGenerator>>,
) -> Result<SparseCountFingerprint, TopologicalTorsionError> {
    let args = FingerprintFuncArguments {
        from_atoms: call.from_atoms,
        ignore_atoms: call.ignore_atoms,
        custom_atom_invariants: call.custom_atom_invariants,
        custom_bond_invariants: call.custom_bond_invariants,
        conformer_id: call.conformer_id,
    };
    with_fingerprint_environment_inputs_and_output(
        input.topology,
        input.properties,
        input.valence,
        input.rings,
        common,
        &args,
        output,
        || {
            // Restored source nullptr has no default generator. Only an
            // actually selected path may try to index it; represent no
            // generated invariants by an empty vector, never fabricated codes.
            if configured_invariants == Some(None) {
                return Ok(Vec::new());
            }
            Ok(configured_invariants
                .flatten()
                .or(call.atom_invariants_generator)
                .unwrap_or(AtomPairAtomInvariantsGenerator {
                    include_chirality: params.include_chirality,
                    topological_torsion_correction: true,
                })
                .atom_invariants(input)?)
        },
        || Ok(Vec::new()),
        |topology, properties, atoms, bonds, output| {
            let prepared = AtomPairPreparedInput {
                topology,
                properties,
                ..*input
            };
            let envs = environments(
                &prepared,
                params,
                call.from_atoms,
                call.ignore_atoms,
                atoms,
                fp_size != 0,
                mode,
            )?;
            let length = if fp_size != 0 {
                fp_size
            } else {
                let size = params.result_size()?;
                if matches!(mode, TorsionCodeMode::Modern) {
                    size
                } else {
                    // RDKit❗✔️: sz -= 1;
                    // Preserve the source legacy compatibility off-by-one at
                    // construction; no second map allocation or resize copy.
                    size - 1
                }
            };
            let mut result = SparseCountFingerprint::new(length);
            accumulate_sparse_counts_into(
                envs,
                common,
                atoms,
                bonds,
                fp_size,
                output,
                &mut result,
            )?;
            Ok(result)
        },
    )
}
pub fn topological_torsion_sparse_count(
    input: &AtomPairPreparedInput<'_>,
    params: &TopologicalTorsionParams,
    call: &TopologicalTorsionCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<SparseCountFingerprint, TopologicalTorsionError> {
    let common = params.common()?;
    count_helper(
        input,
        params,
        &common,
        call,
        0,
        output,
        TorsionCodeMode::Modern,
    )
}
pub fn topological_torsion_sparse_bits(
    input: &AtomPairPreparedInput<'_>,
    params: &TopologicalTorsionParams,
    call: &TopologicalTorsionCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<SparseBitFingerprint, TopologicalTorsionError> {
    let common = params.common()?;
    project_sparse_fingerprint(
        &common,
        params.result_size()?,
        input.topology.atoms.len(),
        output,
        |fp_size, output| {
            count_helper(
                input,
                params,
                &common,
                call,
                fp_size,
                output,
                TorsionCodeMode::Modern,
            )
        },
    )
}
pub fn topological_torsion_count(
    input: &AtomPairPreparedInput<'_>,
    params: &TopologicalTorsionParams,
    call: &TopologicalTorsionCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<SparseCountFingerprint32, TopologicalTorsionError> {
    let common = params.common()?;
    project_count_fingerprint(
        &common,
        input.topology.atoms.len(),
        output,
        |fp_size, output| {
            count_helper(
                input,
                params,
                &common,
                call,
                fp_size,
                output,
                TorsionCodeMode::Modern,
            )
        },
    )
}
pub fn topological_torsion_bits(
    input: &AtomPairPreparedInput<'_>,
    params: &TopologicalTorsionParams,
    call: &TopologicalTorsionCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<Fingerprint, TopologicalTorsionError> {
    let common = params.common()?;
    project_fingerprint(
        &common,
        input.topology.atoms.len(),
        output,
        |fp_size, output| {
            count_helper(
                input,
                params,
                &common,
                call,
                fp_size,
                output,
                TorsionCodeMode::Modern,
            )
        },
    )
}

/// Source parameters for the deprecated torsion algorithms, with canonical names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegacyTopologicalTorsionParams {
    pub torsion_atom_count: u32,
    pub include_chirality: bool,
    pub fp_size: u32,
    pub bits_per_entry: u32,
}
impl Default for LegacyTopologicalTorsionParams {
    fn default() -> Self {
        Self {
            torsion_atom_count: 4,
            include_chirality: false,
            fp_size: 2048,
            bits_per_entry: 4,
        }
    }
}
fn legacy_arguments(params: &LegacyTopologicalTorsionParams) -> TopologicalTorsionParams {
    TopologicalTorsionParams {
        torsion_atom_count: params.torsion_atom_count,
        include_chirality: params.include_chirality,
        fp_size: params.fp_size,
        ..Default::default()
    }
}
fn legacy_invariant_precondition(
    input: &AtomPairPreparedInput<'_>,
    call: &TopologicalTorsionCall<'_>,
) -> Result<(), TopologicalTorsionError> {
    // RDKit❗✔️:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit❗✔️:                "bad atomInvariants size");
    if call
        .custom_atom_invariants
        .is_some_and(|v| v.len() < input.topology.atoms.len())
    {
        return Err(FingerprintError::PreconditionViolation {
            what: "bad atomInvariants size",
        }
        .into());
    }
    Ok(())
}
/// Legacy unfolded algorithm including its source-compatible vector length.
pub fn legacy_topological_torsion_sparse_count(
    input: &AtomPairPreparedInput<'_>,
    params: &LegacyTopologicalTorsionParams,
    call: &TopologicalTorsionCall<'_>,
) -> Result<SparseCountFingerprint, TopologicalTorsionError> {
    // RDKit source: AtomPairs.cpp lines 159-265
    // A Rust deprecation attribute provides the source warning at API use
    // time; it is compile-time instead of RDKit's runtime log emission.
    // RDKit❗✔️:   RDLog::deprecationWarning("please use TopologicalTorsionGenerator");
    // RDKit❗✔️:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit❗✔️:                "bad atomInvariants size");
    // RDKit❗✔️:   const ROMol *lmol = &mol;
    // RDKit❗✔️:   std::unique_ptr<ROMol> tmol;
    // RDKit❗✔️:   if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:     tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    // RDKit❗✔️:     MolOps::assignStereochemistry(*tmol);
    // RDKit❗✔️:     lmol = tmol.get();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   boost::uint64_t sz = 1;
    // RDKit❗✔️:   sz = (sz << (targetSize *
    // RDKit❗✔️:                (codeSize + (includeChirality ? numChiralBits : 0))));
    // RDKit❗✔️:   // NOTE: this -1 is incorrect but it's needed for backwards compatibility.
    // RDKit❗✔️:   //  hopefully we'll never have a case with a torsion that hits this.
    // RDKit❗✔️:   //
    // RDKit❗✔️:   //  mmm, bug compatible.
    // RDKit❗✔️:   sz -= 1;
    // RDKit❗✔️:   auto *res = new SparseIntVect<boost::int64_t>(sz);
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<std::uint32_t> atomCodes;
    // RDKit❗✔️:   atomCodes.reserve(lmol->getNumAtoms());
    // RDKit❗✔️:   for (ROMol::ConstAtomIterator atomItI = lmol->beginAtoms();
    // RDKit❗✔️:        atomItI != lmol->endAtoms(); ++atomItI) {
    // RDKit❗✔️:     if (!atomInvariants) {
    // RDKit❗✔️:       atomCodes.push_back(getAtomCode(*atomItI, 0, includeChirality));
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // need to add to the atomCode here because we subtract off up to 2 below
    // RDKit❗✔️:       // as part of the branch correction
    // RDKit❗✔️:       atomCodes.push_back(
    // RDKit❗✔️:           (*atomInvariants)[(*atomItI)->getIdx()] % ((1 << codeSize) - 1) + 2);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   boost::dynamic_bitset<> *fromAtomsBV = nullptr;
    // RDKit❗✔️:   if (fromAtoms) {
    // RDKit❗✔️:     fromAtomsBV = new boost::dynamic_bitset<>(lmol->getNumAtoms());
    // RDKit❗✔️:     for (auto fAt : *fromAtoms) {
    // RDKit❗✔️:       fromAtomsBV->set(fAt);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   boost::dynamic_bitset<> *ignoreAtomsBV = nullptr;
    // RDKit❗✔️:   if (ignoreAtoms) {
    // RDKit❗✔️:     ignoreAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
    // RDKit❗✔️:     for (auto fAt : *ignoreAtoms) {
    // RDKit❗✔️:       ignoreAtomsBV->set(fAt);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   boost::dynamic_bitset<> pAtoms(lmol->getNumAtoms());
    // RDKit❗✔️:   PATH_LIST paths = findAllPathsOfLengthN(*lmol, targetSize, false);
    // RDKit❗✔️:   for (PATH_LIST::const_iterator pathIt = paths.begin(); pathIt != paths.end();
    // RDKit❗✔️:        ++pathIt) {
    // RDKit❗✔️:     bool keepIt = true;
    // RDKit❗✔️:     if (fromAtomsBV) {
    // RDKit❗✔️:       keepIt = false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     std::vector<std::uint32_t> pathCodes;
    // RDKit❗✔️:     const PATH_TYPE &path = *pathIt;
    // RDKit❗✔️:     if (fromAtomsBV) {
    // RDKit❗✔️:       if (fromAtomsBV->test(static_cast<std::uint32_t>(path.front())) ||
    // RDKit❗✔️:           fromAtomsBV->test(static_cast<std::uint32_t>(path.back()))) {
    // RDKit❗✔️:         keepIt = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (keepIt && ignoreAtomsBV) {
    // RDKit❗✔️:       for (auto pElem : path) {
    // RDKit❗✔️:         if (ignoreAtomsBV->test(pElem)) {
    // RDKit❗✔️:           keepIt = false;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (keepIt) {
    // RDKit❗✔️:       pAtoms.reset();
    // RDKit❗✔️:       for (auto pIt = path.begin(); pIt < path.end(); ++pIt) {
    // RDKit❗✔️:         // look for a cycle that doesn't start at the first atom
    // RDKit❗✔️:         // we can't effectively canonicalize these at the moment
    // RDKit❗✔️:         // (was github #811)
    // RDKit❗✔️:         if (pIt != path.begin() && *pIt != *(path.begin()) && pAtoms[*pIt]) {
    // RDKit❗✔️:           pathCodes.clear();
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         pAtoms.set(*pIt);
    // RDKit❗✔️:         unsigned int code = atomCodes[*pIt] - 1;
    // RDKit❗✔️:         // subtract off the branching number:
    // RDKit❗✔️:         if (pIt != path.begin() && pIt + 1 != path.end()) {
    // RDKit❗✔️:           --code;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         pathCodes.push_back(code);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (pathCodes.size()) {
    // RDKit❗✔️:         boost::int64_t code =
    // RDKit❗✔️:             getTopologicalTorsionCode(pathCodes, includeChirality);
    // RDKit❗✔️:         updateElement(*res, code);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   delete fromAtomsBV;
    // RDKit❗✔️:   delete ignoreAtomsBV;
    // RDKit❗✔️:
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    legacy_invariant_precondition(input, call)?;
    let args = legacy_arguments(params);
    let common = args.common()?;
    let mode = if call.custom_atom_invariants.is_some() {
        TorsionCodeMode::LegacyCustom
    } else {
        TorsionCodeMode::LegacyUnfolded
    };
    // Default AP invariants already subtract two. The legacy endpoint adds
    // one and internal leaves it unchanged, exactly atomCode-1/-2 without
    // modulo. Custom invariants follow source (inv%511)+2 then -1/-2.
    // Source legacy lmol is the invariant provider as well as the
    // environment input. Reuse the sole conditional preparation owner;
    // modern calls retain their original-input provider unchanged. The
    // no-copy Done-present branch borrows and the missing-Done branch
    // clones only topology/properties once; the inner helper then borrows.
    let prepared = crate::prepared::prepare_morgan_environment(
        input.topology,
        input.properties,
        input.valence,
        input.rings,
        params.include_chirality,
    )?;
    let prepared_input = AtomPairPreparedInput {
        topology: prepared.topology(),
        properties: prepared.properties(),
        ..*input
    };
    count_helper(&prepared_input, &args, &common, call, 0, None, mode)
}
/// Legacy hashed count adapter to the single modern count accumulation owner.
pub fn legacy_topological_torsion_count(
    input: &AtomPairPreparedInput<'_>,
    params: &LegacyTopologicalTorsionParams,
    call: &TopologicalTorsionCall<'_>,
) -> Result<SparseCountFingerprint, TopologicalTorsionError> {
    // RDKit source: AtomPairs.cpp lines 267-297
    // RDKit❗✔️: template <typename T>
    // RDKit❗✔️: void TorsionFpCalc(T *res, const ROMol &mol, unsigned int nBits,
    // RDKit❗✔️:                    unsigned int targetSize,
    // RDKit❗✔️:                    const std::vector<std::uint32_t> *fromAtoms,
    // RDKit❗✔️:                    const std::vector<std::uint32_t> *ignoreAtoms,
    // RDKit❗✔️:                    const std::vector<std::uint32_t> *atomInvariants,
    // RDKit❗✔️:                    bool includeChirality) {
    // RDKit❗✔️:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit❗✔️:                "bad atomInvariants size");
    // RDKit❗✔️:   const ROMol *lmol = &mol;
    // RDKit❗✔️:   std::unique_ptr<ROMol> tmol;
    // RDKit❗✔️:   if (includeChirality && !mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:     tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    // RDKit❗✔️:     MolOps::assignStereochemistry(*tmol);
    // RDKit❗✔️:     lmol = tmol.get();
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpgen{
    // RDKit❗✔️:       RDKit::TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(
    // RDKit❗✔️:           includeChirality, targetSize, nullptr, true, nBits)};
    // RDKit❗✔️:   FingerprintFuncArguments args;
    // RDKit❗✔️:   args.fromAtoms = fromAtoms;
    // RDKit❗✔️:   args.ignoreAtoms = ignoreAtoms;
    // RDKit❗✔️:   args.customAtomInvariants = atomInvariants;
    // RDKit❗✔️:   auto siv = fpgen->getCountFingerprint(*lmol, args);
    // RDKit❗✔️:   for (auto v : siv->getNonzeroElements()) {
    // RDKit❗✔️:     res->setVal(v.first, v.second);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Returning the shared generator's owned count vector preserves every
    // element and its declared size while eliminating the source adapter's
    // second sparse-container allocation and O(k) element copy.
    legacy_invariant_precondition(input, call)?;
    if params.fp_size == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "legacy hashed topological torsion nBits must be greater than zero",
        }
        .into());
    }
    let args = legacy_arguments(params);
    let common = args.common()?;
    // Accumulate directly into the returned source-width 64-bit map: same
    // hash/modulo/count order, avoiding C++ temporary32-bit map and copy.
    // Source legacy lmol is the invariant provider as well as the
    // environment input. Reuse the sole conditional preparation owner;
    // modern calls retain their original-input provider unchanged. The
    // no-copy Done-present branch borrows and the missing-Done branch
    // clones only topology/properties once; the inner helper then borrows.
    let prepared = crate::prepared::prepare_morgan_environment(
        input.topology,
        input.properties,
        input.valence,
        input.rings,
        params.include_chirality,
    )?;
    let prepared_input = AtomPairPreparedInput {
        topology: prepared.topology(),
        properties: prepared.properties(),
        ..*input
    };
    count_helper(
        &prepared_input,
        &args,
        &common,
        call,
        u64::from(params.fp_size),
        None,
        TorsionCodeMode::Modern,
    )
}
/// Source threshold projection for legacy hashed torsion bits.
pub fn legacy_topological_torsion_bits(
    input: &AtomPairPreparedInput<'_>,
    params: &LegacyTopologicalTorsionParams,
    call: &TopologicalTorsionCall<'_>,
) -> Result<Fingerprint, TopologicalTorsionError> {
    // RDKit source: AtomPairs.cpp lines 312-347
    // RDKit❗✔️: ExplicitBitVect *getHashedTopologicalTorsionFingerprintAsBitVect(
    // RDKit❗✔️:     const ROMol &mol, unsigned int nBits, unsigned int targetSize,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *fromAtoms,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *ignoreAtoms,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *atomInvariants,
    // RDKit❗✔️:     unsigned int nBitsPerEntry, bool includeChirality) {
    // RDKit❗✔️:   RDLog::deprecationWarning("please use TopologicalTorsionGenerator");
    // RDKit❗✔️:   PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit❗✔️:                "bad atomInvariants size");
    // RDKit❗✔️:   static int bounds[4] = {1, 2, 4, 8};
    // RDKit❗✔️:   unsigned int blockLength = nBits / nBitsPerEntry;
    // RDKit❗✔️:   auto *sres = new SparseIntVect<boost::int64_t>(blockLength);
    // RDKit❗✔️:   TorsionFpCalc(sres, mol, blockLength, targetSize, fromAtoms, ignoreAtoms,
    // RDKit❗✔️:                 atomInvariants, includeChirality);
    // RDKit❗✔️:   auto *res = new ExplicitBitVect(nBits);
    // RDKit❗✔️:
    // RDKit❗✔️:   if (nBitsPerEntry != 4) {
    // RDKit❗✔️:     for (auto val : sres->getNonzeroElements()) {
    // RDKit❗✔️:       for (unsigned int i = 0; i < nBitsPerEntry; ++i) {
    // RDKit❗✔️:         if (val.second > static_cast<int>(i)) {
    // RDKit❗✔️:           res->setBit(val.first * nBitsPerEntry + i);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     for (auto val : sres->getNonzeroElements()) {
    // RDKit❗✔️:       for (unsigned int i = 0; i < nBitsPerEntry; ++i) {
    // RDKit❗✔️:         if (val.second >= bounds[i]) {
    // RDKit❗✔️:           res->setBit(val.first * nBitsPerEntry + i);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   delete sres;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // Rust's deprecation attribute reports use at compile time. Invalid zero
    // divisors and zero-sized hash blocks are returned as structured errors
    // instead of entering source-undefined arithmetic/generator states.
    legacy_invariant_precondition(input, call)?;
    if params.bits_per_entry == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "legacy topological torsion nBitsPerEntry must be greater than zero",
        }
        .into());
    }
    let block_length = params.fp_size / params.bits_per_entry;
    if block_length == 0 {
        return Err(FingerprintError::InvalidArguments {reason:"legacy topological torsion bit vector requires at least one complete entry block"}.into());
    }
    let counts = legacy_topological_torsion_count(
        input,
        &LegacyTopologicalTorsionParams {
            fp_size: block_length,
            ..params.clone()
        },
        call,
    )?;
    const BOUNDS: [i32; 4] = [1, 2, 4, 8];
    let mut on_bits = Vec::new();
    for (&bit, &count) in counts.nonzero_elements() {
        for i in 0..params.bits_per_entry {
            let set = if params.bits_per_entry == 4 {
                count >= BOUNDS[i as usize]
            } else {
                i64::from(count) > i64::from(i)
            };
            if set {
                on_bits.push((bit * u64::from(params.bits_per_entry) + u64::from(i)) as u32);
            }
        }
    }
    Fingerprint::from_on_bits(params.fp_size, on_bits).map_err(Into::into)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::test_support::TestMolecule;
    #[derive(Debug)]
    struct ObservedEnvironment {
        bit_id: u64,
        atom_path: Vec<usize>,
    }
    fn default_invariants(molecule: &TestMolecule) -> Vec<u32> {
        AtomPairAtomInvariantsGenerator {
            include_chirality: false,
            topological_torsion_correction: true,
        }
        .atom_invariants(&molecule.input())
        .unwrap()
    }
    fn observed_environments(
        molecule: &TestMolecule,
        args: &TopologicalTorsionParams,
        from: Option<&[u32]>,
        ignore: Option<&[u32]>,
        invariants: &[u32],
        hashed: bool,
        mode: TorsionCodeMode,
    ) -> Vec<ObservedEnvironment> {
        super::environments(
            &molecule.input(),
            args,
            from,
            ignore,
            invariants,
            hashed,
            mode,
        )
        .unwrap()
        .into_iter()
        .map(|e| ObservedEnvironment {
            bit_id: e.bit_id,
            atom_path: e.path.into_iter().map(|a| a.index()).collect(),
        })
        .collect()
    }
    fn environments(
        molecule: &TestMolecule,
        args: &TopologicalTorsionParams,
        from: Option<&[u32]>,
        ignore: Option<&[u32]>,
        invariants: &[u32],
        hashed: bool,
    ) -> Vec<ObservedEnvironment> {
        observed_environments(
            molecule,
            args,
            from,
            ignore,
            invariants,
            hashed,
            TorsionCodeMode::Modern,
        )
    }
    fn test_environment(bit_id: u64, rows: Vec<usize>) -> TorsionEnvironment {
        TorsionEnvironment {
            bit_id,
            path: rows.into_iter().map(AtomId::new).collect(),
        }
    }
    #[test]
    fn null_and_empty_from_atom_selections_are_distinct() {
        let molecule = TestMolecule::from_smiles("CCCCC").expect("chain");
        let arguments = TopologicalTorsionParams::default();
        let invariants = default_invariants(&molecule);

        let all = environments(&molecule, &arguments, None, None, &invariants, false);
        let empty = environments(&molecule, &arguments, Some(&[]), None, &invariants, false);

        assert_eq!(all.len(), 2);
        assert!(empty.is_empty());
        assert_eq!(all[0].atom_path, vec![0, 1, 2, 3]);
        assert_eq!(all[1].atom_path, vec![1, 2, 3, 4]);
    }

    #[test]
    fn from_atoms_match_only_path_endpoints() {
        let molecule = TestMolecule::from_smiles("CCCCC").expect("chain");
        let arguments = TopologicalTorsionParams::default();
        let invariants = default_invariants(&molecule);

        let internal_only =
            environments(&molecule, &arguments, Some(&[2]), None, &invariants, false);
        let endpoint_and_internal =
            environments(&molecule, &arguments, Some(&[3]), None, &invariants, false);

        assert!(internal_only.is_empty());
        assert_eq!(endpoint_and_internal.len(), 1);
        assert_eq!(endpoint_and_internal[0].atom_path, vec![0, 1, 2, 3]);
    }

    #[test]
    fn ignored_atoms_remove_paths_when_the_atom_is_internal_or_terminal() {
        let molecule = TestMolecule::from_smiles("CCCCC").expect("chain");
        let arguments = TopologicalTorsionParams::default();
        let invariants = default_invariants(&molecule);

        assert!(
            environments(&molecule, &arguments, None, Some(&[2]), &invariants, false,).is_empty()
        );
        let terminal = environments(&molecule, &arguments, None, Some(&[0]), &invariants, false);
        assert_eq!(terminal.len(), 1);
        assert_eq!(terminal[0].atom_path, vec![1, 2, 3, 4]);
    }

    #[test]
    fn custom_invariants_receive_endpoint_and_internal_source_corrections() {
        let molecule = TestMolecule::from_smiles("CCCCC").expect("chain");
        let arguments = TopologicalTorsionParams::default();
        let invariants = [10, 20, 30, 40, 50];
        let generated = environments(&molecule, &arguments, None, None, &invariants, false);

        assert_eq!(
            generated.iter().map(|env| env.bit_id).collect::<Vec<_>>(),
            vec![
                topological_torsion_code(&[11, 20, 30, 41], false).unwrap(),
                topological_torsion_code(&[21, 30, 40, 51], false).unwrap(),
            ]
        );
    }

    #[test]
    fn legacy_unfolded_atom_codes_bypass_the_modern_base_code_modulus() {
        let molecule = TestMolecule::from_smiles("CCCC").expect("chain");
        let mut arguments = TopologicalTorsionParams::default();
        arguments.include_chirality = true;
        let invariants = [600, 600, 600, 600];

        let modern = observed_environments(
            &molecule,
            &arguments,
            None,
            None,
            &invariants,
            false,
            TorsionCodeMode::Modern,
        );
        let legacy = observed_environments(
            &molecule,
            &arguments,
            None,
            None,
            &invariants,
            false,
            TorsionCodeMode::LegacyUnfolded,
        );

        assert_eq!(modern.len(), 1);
        assert_eq!(legacy.len(), 1);
        assert_eq!(
            modern[0].bit_id,
            topological_torsion_code(&[90, 89, 89, 90], true).unwrap()
        );
        assert_eq!(
            legacy[0].bit_id,
            topological_torsion_code(&[601, 600, 600, 601], true).unwrap()
        );
        assert_ne!(modern[0].bit_id, legacy[0].bit_id);
    }

    #[test]
    fn full_ring_closure_is_kept_but_a_later_illegal_repeat_is_not_enumerated() {
        let molecule = TestMolecule::from_smiles("C1CCC1").expect("four-membered ring");
        let invariants = default_invariants(&molecule);
        let mut arguments = TopologicalTorsionParams::default();
        arguments.torsion_atom_count = 5;

        let closed = environments(&molecule, &arguments, None, None, &invariants, false);
        assert_eq!(closed.len(), 1);
        assert_eq!(closed[0].atom_path, vec![0, 1, 2, 3, 0]);

        arguments.torsion_atom_count = 6;
        assert!(environments(&molecule, &arguments, None, None, &invariants, false).is_empty());
    }

    #[test]
    fn hashed_and_unhashed_ids_preserve_source_collisions() {
        let molecule = TestMolecule::from_smiles("CC(C)C").expect("branched molecule");
        let mut arguments = TopologicalTorsionParams::default();
        arguments.torsion_atom_count = 3;
        let invariants = vec![7; molecule.num_atoms()];

        let packed = environments(&molecule, &arguments, None, None, &invariants, false);
        let hashed = environments(&molecule, &arguments, None, None, &invariants, true);
        assert_eq!(packed.len(), 3);
        assert_eq!(hashed.len(), 3);
        assert!(packed.iter().all(|env| env.bit_id == packed[0].bit_id));
        assert!(hashed.iter().all(|env| env.bit_id == hashed[0].bit_id));
        assert_eq!(
            packed[0].bit_id,
            topological_torsion_code(&[8, 7, 8], false).unwrap()
        );
        assert_eq!(
            hashed[0].bit_id,
            u64::from(topological_torsion_hash(&[8, 7, 8]).unwrap())
        );
    }

    #[test]
    fn atom_environment_updates_every_supported_provenance_allocation() {
        let environment = test_environment(99, vec![0, 2, 4]);
        let mut output = FingerprintAdditionalOutput::default();
        output.allocate_atom_to_bits();
        output.allocate_atom_counts();
        output.allocate_bit_paths();
        output.allocate_atoms_per_bit();
        output.reinitialize(5);

        environment.update_output(&mut output, 7, &mut ()).unwrap();
        assert_eq!(environment.bit_id, 99);
        assert_eq!(
            output.atom_to_bits,
            Some(vec![vec![7], vec![], vec![7], vec![], vec![7]])
        );
        assert_eq!(output.atom_counts, Some(vec![1, 0, 1, 0, 1]));
        assert_eq!(
            output.bit_paths.as_ref().unwrap().get(&7),
            Some(&vec![vec![0, 2, 4]])
        );
        assert_eq!(
            output.atoms_per_bit.as_ref().unwrap().get(&7),
            Some(&vec![vec![0, 2, 4]])
        );
    }

    #[test]
    fn provenance_keeps_duplicate_colliding_paths_and_repeated_atoms() {
        let mut output = FingerprintAdditionalOutput::default();
        output.allocate_atom_to_bits();
        output.allocate_atom_counts();
        output.allocate_bit_paths();
        output.allocate_atoms_per_bit();
        output.reinitialize(3);

        test_environment(4, vec![0, 1, 0])
            .update_output(&mut output, 5, &mut ())
            .unwrap();
        test_environment(4, vec![1, 2])
            .update_output(&mut output, 5, &mut ())
            .unwrap();

        assert_eq!(
            output.atom_to_bits,
            Some(vec![vec![5, 5], vec![5, 5], vec![5]])
        );
        assert_eq!(output.atom_counts, Some(vec![2, 2, 1]));
        assert_eq!(
            output.bit_paths.as_ref().unwrap().get(&5),
            Some(&vec![vec![0, 1, 0], vec![1, 2]])
        );
        assert_eq!(output.atoms_per_bit, output.bit_paths);
    }
}

#[cfg(test)]
mod argument_tests {
    use super::*;
    #[test]
    fn defaults_and_constructor_match_the_pinned_source_contract() {
        let arguments = TopologicalTorsionParams::default();
        let common = &arguments;

        assert_eq!(arguments.torsion_atom_count, 4);
        assert!(!arguments.only_shortest_paths);
        assert!(common.count_simulation);
        assert!(!common.include_chirality);
        assert_eq!(common.count_bounds, vec![1, 2, 4, 8]);
        assert_eq!(common.fp_size, 2048);
        assert_eq!(common.bits_per_feature, 1);

        let constructed = TopologicalTorsionParams::new(true, 5, false, vec![2, 3, 7], 4096)
            .expect("valid custom arguments");
        assert_eq!(constructed.torsion_atom_count, 5);
        assert!(!constructed.only_shortest_paths);
        assert!(constructed.include_chirality);
        assert!(!constructed.count_simulation);
        assert_eq!(constructed.count_bounds, vec![2, 3, 7]);
        assert_eq!(constructed.fp_size, 4096);
        assert_eq!(constructed.bits_per_feature, 1);

        assert!(matches!(
            TopologicalTorsionParams::new(false, 4, true, Vec::new(), 2048),
            Err(FingerprintError::PreconditionViolation { .. })
        ));
    }

    #[test]
    fn mutable_options_and_information_strings_use_the_shared_common_arguments() {
        let mut arguments = TopologicalTorsionParams::default();
        arguments.torsion_atom_count = 3;
        arguments.only_shortest_paths = true;
        arguments.count_simulation = false;
        arguments.include_chirality = true;
        arguments.count_bounds = vec![3, 9];
        arguments.fp_size = 1024;
        arguments.bits_per_feature = 2;

        assert_eq!(
            arguments.info_string(),
            "TopologicalTorsionArguments torsionAtomCount=3 onlyShortestPaths=1"
        );
        assert_eq!(arguments.info_string(), arguments.info_string());
        assert_eq!(
            arguments.common_info_string(),
            "Common arguments : countSimulation=0 fpSize=1024 bitsPerFeature=2 includeChirality=1"
        );
    }

    #[test]
    fn json_roundtrip_updates_derived_and_common_fields() {
        let mut arguments = TopologicalTorsionParams::default();
        arguments.torsion_atom_count = 5;
        arguments.only_shortest_paths = true;
        arguments.count_simulation = false;
        arguments.include_chirality = true;
        arguments.count_bounds = vec![2, 6, 10];
        arguments.fp_size = 4096;
        arguments.bits_per_feature = 3;

        let json = arguments.to_json();
        let value: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
        assert_eq!(value["type"], "TopologicalTorsionArguments");
        assert_eq!(value["torsionAtomCount"], "5");
        assert_eq!(value["onlyShortestPaths"], "true");
        assert_eq!(value["countSimulation"], "false");
        assert_eq!(value["includeChirality"], "true");
        assert_eq!(value["countBounds"], serde_json::json!(["2", "6", "10"]));
        assert_eq!(value["fpSize"], "4096");
        assert_eq!(value["numBitsPerFeature"], "3");
        assert_eq!(arguments.to_json(), json);

        let mut restored = TopologicalTorsionParams::default();
        restored.from_json(&json).expect("roundtrip");
        assert_eq!(restored, arguments);
    }

    #[test]
    fn partial_json_preserves_scalar_defaults_and_clears_missing_count_bounds_like_source() {
        let mut arguments = TopologicalTorsionParams::default();
        arguments
            .from_json(r#"{"torsionAtomCount":5,"onlyShortestPaths":true,"fpSize":1024}"#)
            .expect("partial update");

        assert_eq!(arguments.torsion_atom_count, 5);
        assert!(arguments.only_shortest_paths);
        assert_eq!(arguments.fp_size, 1024);
        assert!(arguments.count_simulation);
        assert!(!arguments.include_chirality);
        assert_eq!(arguments.bits_per_feature, 1);
        assert!(arguments.count_bounds.is_empty());

        arguments
            .from_json(r#"{"countBounds":[1,"4",8],"includeChirality":"true"}"#)
            .expect("property-tree-compatible string values");
        assert_eq!(arguments.count_bounds, vec![1, 4, 8]);
        assert!(arguments.include_chirality);
    }

    #[test]
    fn malformed_json_and_invalid_field_types_return_structured_errors() {
        let mut arguments = TopologicalTorsionParams::default();
        assert!(matches!(
            arguments.from_json("{"),
            Err(FingerprintJsonError::Parse(_))
        ));
        arguments.from_json("[]").unwrap();
        assert_eq!(arguments.torsion_atom_count, 4);
        assert!(arguments.count_bounds.is_empty());
        arguments.from_json(r#"{"torsionAtomCount":-1}"#).unwrap();
        assert_eq!(arguments.torsion_atom_count, u32::MAX);
        arguments
            .from_json(r#"{"onlyShortestPaths":"unknown"}"#)
            .unwrap();
        assert!(!arguments.only_shortest_paths);
        arguments.from_json(r#"{"countBounds":{}}"#).unwrap();
        assert!(arguments.count_bounds.is_empty());
        for json in [
            r#"{"countBounds":[{}]}"#,
            r#"{"countBounds":["invalid"]}"#,
            r#"{"countBounds":[4294967296]}"#,
        ] {
            assert!(matches!(
                arguments.from_json(json),
                Err(FingerprintJsonError::Invalid(_))
            ));
        }
    }

    #[test]
    fn argument_unit_ignores_unknown_type_tags_and_fields_like_the_source_method() {
        let mut arguments = TopologicalTorsionParams::default();
        arguments
        .from_json(
            r#"{"type":"UnknownArguments","unknownOption":17,"torsionAtomCount":2,"countBounds":[1]}"#,
        )
        .expect("the argument unit does not dispatch on type");

        assert_eq!(arguments.torsion_atom_count, 2);
        assert_eq!(arguments.count_bounds, vec![1]);
    }

    #[test]
    fn result_sizes_cover_zero_default_chiral_and_maximum_defined_shifts() {
        let mut arguments = TopologicalTorsionParams::default();

        assert_eq!(arguments.result_size().unwrap(), 1_u64 << 36);
        assert_eq!(arguments.result_size().unwrap(), 1_u64 << 36);

        arguments.torsion_atom_count = 0;
        assert_eq!(arguments.result_size().unwrap(), 1);

        arguments.torsion_atom_count = 7;
        assert_eq!(arguments.result_size().unwrap(), 1_u64 << 63);

        arguments.include_chirality = true;
        arguments.torsion_atom_count = 5;
        assert_eq!(arguments.result_size().unwrap(), 1_u64 << 55);
    }

    #[test]
    fn result_size_rejects_every_source_undefined_width_mapping() {
        let mut arguments = TopologicalTorsionParams::default();

        arguments.torsion_atom_count = 8;
        assert!(matches!(
            arguments.result_size(),
            Err(FingerprintError::InvalidArguments { .. })
        ));

        arguments.include_chirality = true;
        arguments.torsion_atom_count = 6;
        assert!(matches!(
            arguments.result_size(),
            Err(FingerprintError::InvalidArguments { .. })
        ));

        arguments.torsion_atom_count = u32::MAX;
        assert!(matches!(
            arguments.result_size(),
            Err(FingerprintError::InvalidArguments { .. })
        ));
    }
}

#[cfg(test)]
#[path = "topological_torsion_legacy_tests.rs"]
mod legacy_tests;

#[path = "topological_torsion_operator.rs"]
mod operator;
pub use operator::{TopologicalTorsionGenerator, TopologicalTorsionSettings};

impl From<FingerprintJsonError> for TopologicalTorsionError {
    fn from(e: FingerprintJsonError) -> Self {
        Self::Json(e)
    }
}

/// Pinned Python source IDs helper over the sole legacy unfolded owner.
pub fn topological_torsion_ids(
    input: &AtomPairPreparedInput<'_>,
    torsion_atom_count: u32,
) -> Result<Vec<u64>, TopologicalTorsionError> {
    // RDKit❗✔️: def GetTopologicalTorsionFingerprintAsIds(mol, targetSize=4):
    // RDKit❗✔️:   nonZeroElements = GetTopologicalTorsionFingerprint(mol, targetSize).GetNonzeroElements()
    // RDKit❗✔️:   frequencies = sorted(nonZeroElements.items())
    // RDKit❗✔️:   res = []
    // RDKit❗✔️:   for k, v in frequencies:
    // RDKit❗✔️:     res.extend([k] * v)
    // RDKit❗✔️:   return res
    // BTreeMap already has source ascending order; visit once without sorting
    // or constructing the source's temporary [k]*v list. O(k+returned IDs).
    // Negative signed source counts multiply to an empty Python list. Checked
    // reservation propagates allocation failure instead of truncating IDs.
    let counts = legacy_topological_torsion_sparse_count(
        input,
        &LegacyTopologicalTorsionParams {
            torsion_atom_count,
            ..Default::default()
        },
        &TopologicalTorsionCall::default(),
    )?;
    let mut result = Vec::new();
    for (&bit, &count) in counts.nonzero_elements() {
        let repetitions = count.max(0) as usize;
        result
            .try_reserve(repetitions)
            .map_err(TopologicalTorsionError::OutputAllocation)?;
        result.extend(std::iter::repeat_n(bit, repetitions));
    }
    Ok(result)
}
