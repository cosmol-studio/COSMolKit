//! Source-ordered self matches for detached conformer pruning.
use crate::EmbedParams;
use cosmolkit_model::{
    CoordinateBlock, MoleculeProperties, QueryGraph, QueryGraphError, TopologyBlock,
};
use cosmolkit_search::{
    QueryMatchContextError, SearchTarget, SubstructMatchError, SubstructMatchParams,
};
/// Mechanical representation transport supplied by the private facade.
/// A function pointer cannot capture a live molecule, runtime or commit authority.
pub type PruningQueryProjection = fn(
    &TopologyBlock,
    &CoordinateBlock,
    &MoleculeProperties,
) -> Result<QueryGraph, QueryGraphError>;
#[derive(Debug, thiserror::Error)]
pub enum PruningError {
    #[error(transparent)]
    Hydrogens(#[from] cosmolkit_core::HydrogenError),
    #[error(transparent)]
    Projection(#[from] QueryGraphError),
    #[error(transparent)]
    Terminal(#[from] cosmolkit_alignment::AlignmentError),
    #[error(transparent)]
    Match(#[from] SubstructMatchError),
    #[error(transparent)]
    Context(#[from] QueryMatchContextError),
    #[error("{0}")]
    Input(&'static str),
}
#[derive(Debug)]
pub struct PruningMatches {
    pub self_matches: Vec<Vec<usize>>,
    pub hydrogen_warnings: Vec<cosmolkit_core::HydrogenWarning>,
}
#[doc(hidden)]
pub fn pruning_self_matches(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    params: &EmbedParams,
    project: PruningQueryProjection,
) -> Result<PruningMatches, PruningError> {
    // BEGIN RDKIT CPP FUNCTION getMolSelfMatches (Embedder.cpp:1449-1497)
    // RDKit❗❌: std::vector<std::vector<unsigned int>> getMolSelfMatches(
    // RDKit❗❌:     const ROMol &mol, const EmbedParameters &params) {
    // RDKit❗❌:   std::vector<std::vector<unsigned int>> res;
    // RDKit❗❌:   if (params.pruneRmsThresh && params.useSymmetryForPruning) {
    // RDKit❗❌:     RWMol tmol(mol);
    // RDKit❗❌:     MolOps::RemoveHsParameters ps;
    // RDKit❗❌:     bool sanitize = false;
    // RDKit❗❌:     MolOps::removeHs(tmol, ps, sanitize);
    // RDKit❗❌:
    // RDKit❗❌:     std::unique_ptr<RWMol> prbMolSymm;
    // RDKit❗❌:     if (params.symmetrizeConjugatedTerminalGroupsForPruning) {
    // RDKit❗❌:       prbMolSymm.reset(new RWMol(tmol));
    // RDKit❗❌:       MolAlign::details::symmetrizeTerminalAtoms(*prbMolSymm);
    // RDKit❗❌:     }
    // RDKit❗❌:     const auto &prbMolForMatch = prbMolSymm ? *prbMolSymm : tmol;
    // RDKit❗❌:
    // RDKit❗❌:     SubstructMatchParameters sssps;
    // RDKit❗❌:     sssps.maxMatches = 1;
    // RDKit❗❌:     // provides the atom indices in the molecule corresponding
    // RDKit❗❌:     // to the indices in the H-stripped version
    // RDKit❗❌:     auto strippedMatch = SubstructMatch(mol, prbMolForMatch, sssps);
    // RDKit❗❌:     CHECK_INVARIANT(strippedMatch.size() == 1, "expected match not found");
    // RDKit❗❌:
    // RDKit❗❌:     sssps.maxMatches = 1000;
    // RDKit❗❌:     sssps.uniquify = false;
    // RDKit❗❌:     auto heavyAtomMatches = SubstructMatch(tmol, prbMolForMatch, sssps);
    // RDKit❗❌:     for (const auto &match : heavyAtomMatches) {
    // RDKit❗❌:       res.emplace_back(0);
    // RDKit❗❌:       res.back().reserve(match.size());
    // RDKit❗❌:       for (auto midx : match) {
    // RDKit❗❌:         res.back().push_back(strippedMatch[0][midx.second].second);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   } else if (params.onlyHeavyAtomsForRMS) {
    // RDKit❗❌:     res.emplace_back(0);
    // RDKit❗❌:     for (const auto &at : mol.atoms()) {
    // RDKit❗❌:       if (at->getAtomicNum() != 1) {
    // RDKit❗❌:         res.back().push_back(at->getIdx());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     res.emplace_back(0);
    // RDKit❗❌:     res.back().reserve(mol.getNumAtoms());
    // RDKit❗❌:     for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:       res.back().push_back(i);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION getMolSelfMatches (Embedder.cpp:1449-1497)

    let mut result = Vec::new();
    let mut warnings = Vec::new();
    if params.prune_rms_thresh != 0.0 && params.use_symmetry_for_pruning {
        let stripped = cosmolkit_core::remove_hydrogens_with_params(
            topology.clone(),
            coordinates.clone(),
            properties.clone(),
            &cosmolkit_core::RemoveHsParams {
                sanitize: false,
                ..Default::default()
            },
        )?;
        // Source representation projection occurs after RemoveHs, retaining
        // its remapped chemistry and final source-uninitialized ring state.
        // Carrier-derived rows do not turn ordinary Atom/Bond::Match into
        // explicit SMARTS predicates. Only alignment's source terminal rewrite
        // changes bond predicate origins. No mapping shortcut replaces SEARCH.
        let mut query = project(
            &stripped.topology,
            &stripped.coordinates,
            &stripped.properties,
        )?;
        let stripped_target = SearchTarget::new(
            &stripped.topology,
            &stripped.coordinates,
            &stripped.topology.stereo_groups,
            None,
            None,
        );
        if params.symmetrize_conjugated_terminal_groups_for_pruning {
            let terminal_context =
                cosmolkit_search::build_topology_query_match_context(&stripped.topology)?;
            query = cosmolkit_alignment::symmetrize_terminal_query_with_context(
                query,
                &stripped_target,
                &terminal_context,
            )?;
        }
        let original_target =
            SearchTarget::new(topology, coordinates, &topology.stereo_groups, None, None);
        let original_context = cosmolkit_search::build_topology_query_match_context(topology)?;
        let mut matching = SubstructMatchParams {
            max_matches: 1,
            ..Default::default()
        };
        let stripped_match = cosmolkit_search::try_get_substruct_matches_with_params_and_context(
            &original_target,
            &query,
            &matching,
            &original_context,
        )?;
        if stripped_match.len() != 1 {
            return Err(PruningError::Input("expected match not found"));
        }
        matching.max_matches = 1000;
        matching.uniquify = false;
        let stripped_context =
            cosmolkit_search::build_topology_query_match_context(&stripped.topology)?;
        let heavy_matches = cosmolkit_search::try_get_substruct_matches_with_params_and_context(
            &stripped_target,
            &query,
            &matching,
            &stripped_context,
        )?;
        for matched in heavy_matches {
            let mut mapped = Vec::with_capacity(matched.atom_mapping.len());
            for index in matched.atom_mapping {
                mapped.push(
                    *stripped_match[0]
                        .atom_mapping
                        .get(index)
                        .ok_or(PruningError::Input("stripped match atom out of range"))?,
                );
            }
            result.push(mapped);
        }
        warnings = stripped.warnings;
    } else if params.only_heavy_atoms_for_rms {
        result.push(
            topology
                .atoms
                .iter()
                .filter(|atom| atom.atomic_number() != 1)
                .map(|atom| atom.id().index())
                .collect(),
        );
    } else {
        result.push((0..topology.atoms.len()).collect());
    }
    Ok(PruningMatches {
        self_matches: result,
        hydrogen_warnings: warnings,
    })
}
