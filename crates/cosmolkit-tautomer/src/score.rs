//! Source tautomer score terms and signed component accumulation.
use crate::engine::{TautomerRecordView, TautomerRunError, total_hydrogens};
use cosmolkit_core::{RingFindingError, symmetrized_sssr};
use cosmolkit_model::QueryGraph;
use cosmolkit_search::{CompiledQuery, SmartsParseParams, parse_smarts};
use std::sync::OnceLock;

/// One named SMARTS contribution to tautomer canonicalization.
#[derive(Debug, Clone)]
pub struct TautomerScoreTerm {
    name: String,
    smarts: String,
    score: i32,
    matcher: Option<QueryGraph>,
    required_elements: Vec<u8>,
    connectivity_smarts: String,
    connectivity_matcher: Option<QueryGraph>,
}

/// Source-defined components of a tautomer canonicalization score.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TautomerScore {
    ring: i32,
    substructure: i32,
    hetero_hydrogen: i32,
}
impl TautomerScore {
    #[must_use]
    pub const fn ring(self) -> i32 {
        self.ring
    }

    #[must_use]
    pub const fn substructure(self) -> i32 {
        self.substructure
    }

    #[must_use]
    pub const fn hetero_hydrogen(self) -> i32 {
        self.hetero_hydrogen
    }

    #[must_use]
    pub const fn total(self) -> i32 {
        self.ring
            .wrapping_add(self.substructure)
            .wrapping_add(self.hetero_hydrogen)
    }
}

impl TautomerScoreTerm {
    #[must_use]
    pub fn new(name: impl Into<String>, smarts: impl Into<String>, score: i32) -> Self {
        Self::with_prerequisites(name.into(), smarts.into(), score, vec![], String::new())
    }

    fn with_prerequisites(
        name: String,
        smarts: String,
        score: i32,
        required_elements: Vec<u8>,
        connectivity_smarts: String,
    ) -> Self {
        // BEGIN RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::SubstructTerm::metadata_shape
        // RDKit❗✔️: struct RDKIT_MOLSTANDARDIZE_EXPORT SubstructTerm {
        // RDKit❗✔️:   std::string name;
        // RDKit❗✔️:   std::string smarts;
        // RDKit❗✔️:   int score;
        // RDKit❗✔️:   RWMol matcher;  // requires assignment
        // RDKit❗✔️:
        // RDKit❗✔️:   // Pre-screening support: elements that must be present (empty = no filter)
        // RDKit❗✔️:   std::vector<int> requiredElements;
        // RDKit❗✔️:   // Bond-order-agnostic connectivity pattern for pre-screening (may be empty)
        // RDKit❗✔️:   std::string connectivitySmarts;
        // RDKit❗✔️:   RWMol connectivityMatcher;
        // RDKit❗✔️:
        // RDKit❗✔️:   SubstructTerm(std::string aname, std::string asmarts, int ascore,
        // RDKit❗✔️:                 std::vector<int> reqElements = {},
        // RDKit❗✔️:                 std::string connSmarts = "");
        // RDKit❗✔️:   SubstructTerm(const SubstructTerm &rhs) = default;
        // RDKit❗✔️:   SubstructTerm &operator=(const SubstructTerm &rhs) = default;
        // RDKit❗✔️:
        // RDKit❗✔️:   bool operator==(const SubstructTerm &rhs) const {
        // RDKit❗✔️:     return name == rhs.name && smarts == rhs.smarts && score == rhs.score;
        // RDKit❗✔️:   }
        // RDKit❗✔️: };
        // END RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::SubstructTerm::metadata_shape
        // BEGIN RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::SubstructTerm
        // RDKit❗✔️: SubstructTerm::SubstructTerm(std::string aname, std::string asmarts, int ascore,
        // RDKit❗✔️:                              std::vector<int> reqElements,
        // RDKit❗✔️:                              std::string connSmarts)
        // RDKit❗✔️:     : name(std::move(aname)),
        // RDKit❗✔️:       smarts(std::move(asmarts)),
        // RDKit❗✔️:       score(ascore),
        // RDKit❗✔️:       requiredElements(std::move(reqElements)),
        // RDKit❗✔️:       connectivitySmarts(std::move(connSmarts)) {
        // RDKit❗✔️:   std::unique_ptr<ROMol> pattern(SmartsToMol(smarts));
        // RDKit❗✔️:   if (pattern) {
        // RDKit❗✔️:     matcher = std::move(*pattern);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   // Initialize connectivity matcher if provided
        // RDKit❗✔️:   if (!connectivitySmarts.empty()) {
        // RDKit❗✔️:     std::unique_ptr<ROMol> connPattern(SmartsToMol(connectivitySmarts));
        // RDKit❗✔️:     if (connPattern) {
        // RDKit❗✔️:       connectivityMatcher = std::move(*connPattern);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::SubstructTerm
        let matcher = parse_smarts(&smarts, &SmartsParseParams::default()).ok();
        let connectivity_matcher = if connectivity_smarts.is_empty() {
            None
        } else {
            parse_smarts(&connectivity_smarts, &SmartsParseParams::default()).ok()
        };
        Self {
            name,
            smarts,
            score,
            matcher,
            required_elements,
            connectivity_smarts,
            connectivity_matcher,
        }
    }

    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[must_use]
    pub fn smarts(&self) -> &str {
        &self.smarts
    }

    #[must_use]
    pub const fn score(&self) -> i32 {
        self.score
    }

    pub(crate) const fn matcher(&self) -> Option<&QueryGraph> {
        self.matcher.as_ref()
    }
}

impl PartialEq for TautomerScoreTerm {
    fn eq(&self, other: &Self) -> bool {
        // RDKit✔️✔️: bool operator==(const SubstructTerm &rhs) const {
        // RDKit✔️✔️:   return name == rhs.name && smarts == rhs.smarts && score == rhs.score;
        // RDKit✔️✔️: }
        self.name == other.name && self.smarts == other.smarts && self.score == other.score
    }
}

impl Eq for TautomerScoreTerm {}

#[must_use]
pub fn default_tautomer_score_terms() -> &'static [TautomerScoreTerm] {
    // BEGIN RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::getDefaultTautomerScoreSubstructs
    // RDKit❗✔️: const std::vector<SubstructTerm> &getDefaultTautomerScoreSubstructs() {
    // RDKit❗✔️:   // Each term specifies:
    // RDKit❗✔️:   //   - name, SMARTS, score
    // RDKit❗✔️:   //   - requiredElements: atomic numbers that must be present (empty = no filter)
    // RDKit❗✔️:   //   - connectivitySmarts: bond-order-agnostic pattern for pre-screening
    // RDKit❗✔️:   // Since tautomerization only moves H and changes bond orders (never creates/
    // RDKit❗✔️:   // destroys heavy-atom bonds), we can skip patterns whose connectivity
    // RDKit❗✔️:   // prerequisites aren't met by the input molecule.
    // RDKit❗✔️:   static std::vector<SubstructTerm> substructureTerms{
    // RDKit❗✔️:       {"benzoquinone", "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]", 25, {6}, "[#6]1(~[#6]~[#6]~[#6](~[#6]~[#6]~1)~[N,S,O])~[N,S,O]"},
    // RDKit❗✔️:       {"oxim", "[#6]=[N][OH]", 4, {6, 7, 8}, "[#6]~[#7]~[#8]"},
    // RDKit❗✔️:       {"C=O", "[#6]=,:[#8]", 2, {6, 8}, "[#6]~[#8]"},
    // RDKit❗✔️:       {"N=O", "[#7]=,:[#8]", 2, {7, 8}, "[#7]~[#8]"},
    // RDKit❗✔️:       {"P=O", "[#15]=,:[#8]", 2, {15, 8}, "[#15]~[#8]"},
    // RDKit❗✔️:       {"C=hetero", "[C]=[!#1;!#6]", 1, {6}, "[C]~[!#1;!#6]"},
    // RDKit❗✔️:       {"C(=hetero)-hetero", "[C](=[!#1;!#6])[!#1;!#6]", 2, {6}, "[C](~[!#1;!#6])~[!#1;!#6]"},
    // RDKit❗✔️:       {"aromatic C = exocyclic N", "[c]=!@[N]", -1, {6, 7}, "[c]~[N]"},
    // RDKit❗✔️:       {"methyl", "[CX4H3]", 1, {6}, ""},
    // RDKit❗✔️:       {"guanidine terminal=N", "[#7]C(=[NR0])[#7H0]", 1, {6, 7}, "[#7]~[#6]~[#7]"},
    // RDKit❗✔️:       {"guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2, {6, 7}, "[#7]~[#6](~[N])~[#7]"},
    // RDKit❗✔️:       {"aci-nitro", "[#6]=[N+]([O-])[OH]", -4, {6, 7, 8}, "[#6]~[#7](~[#8])~[#8]"}};
    // RDKit❗✔️:   return substructureTerms;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::getDefaultTautomerScoreSubstructs
    static TERMS: OnceLock<Vec<TautomerScoreTerm>> = OnceLock::new();
    TERMS.get_or_init(|| {
        vec![
            TautomerScoreTerm::with_prerequisites(
                "benzoquinone".to_owned(),
                "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]".to_owned(),
                25,
                vec![6],
                "[#6]1(~[#6]~[#6]~[#6](~[#6]~[#6]~1)~[N,S,O])~[N,S,O]".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "oxim".to_owned(),
                "[#6]=[N][OH]".to_owned(),
                4,
                vec![6, 7, 8],
                "[#6]~[#7]~[#8]".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "C=O".to_owned(),
                "[#6]=,:[#8]".to_owned(),
                2,
                vec![6, 8],
                "[#6]~[#8]".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "N=O".to_owned(),
                "[#7]=,:[#8]".to_owned(),
                2,
                vec![7, 8],
                "[#7]~[#8]".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "P=O".to_owned(),
                "[#15]=,:[#8]".to_owned(),
                2,
                vec![15, 8],
                "[#15]~[#8]".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "C=hetero".to_owned(),
                "[C]=[!#1;!#6]".to_owned(),
                1,
                vec![6],
                "[C]~[!#1;!#6]".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "C(=hetero)-hetero".to_owned(),
                "[C](=[!#1;!#6])[!#1;!#6]".to_owned(),
                2,
                vec![6],
                "[C](~[!#1;!#6])~[!#1;!#6]".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "aromatic C = exocyclic N".to_owned(),
                "[c]=!@[N]".to_owned(),
                -1,
                vec![6, 7],
                "[c]~[N]".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "methyl".to_owned(),
                "[CX4H3]".to_owned(),
                1,
                vec![6],
                "".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "guanidine terminal=N".to_owned(),
                "[#7]C(=[NR0])[#7H0]".to_owned(),
                1,
                vec![6, 7],
                "[#7]~[#6]~[#7]".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "guanidine endocyclic=N".to_owned(),
                "[#7;R][#6;R]([N])=[#7;R]".to_owned(),
                2,
                vec![6, 7],
                "[#7]~[#6](~[N])~[#7]".to_owned(),
            ),
            TautomerScoreTerm::with_prerequisites(
                "aci-nitro".to_owned(),
                "[#6]=[N+]([O-])[OH]".to_owned(),
                -4,
                vec![6, 7, 8],
                "[#6]~[#7](~[#8])~[#8]".to_owned(),
            ),
        ]
    })
}

/// Score aromatic rings using RDKit's tautomer-canonicalization weights.
pub fn score_tautomer_rings_(
    molecule: &mut crate::TautomerScoreView<'_>,
) -> Result<i32, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::scoreRings restricted cache loan
    // RDKit❗❌: int scoreRings(const ROMol &mol) {
    // RDKit❗❌:   int score = 0;
    // RDKit❗❌:   auto ringInfo = mol.getRingInfo();
    // RDKit❗❌:   if (!ringInfo->isSymmSssr()) {
    // RDKit❗❌:     MolOps::symmetrizeSSSR(const_cast<ROMol &>(mol));
    // RDKit❗❌:     ringInfo = mol.getRingInfo();
    // RDKit❗❌:   }
    // RDKit❗❌:   boost::dynamic_bitset<> isArom(mol.getNumBonds());
    // RDKit❗❌:   boost::dynamic_bitset<> bothCarbon(mol.getNumBonds());
    // RDKit❗❌:   for (const auto &bnd : mol.bonds()) {
    // RDKit❗❌:     if (bnd->getIsAromatic()) {
    // RDKit❗❌:       isArom.set(bnd->getIdx());
    // RDKit❗❌:       if (bnd->getBeginAtom()->getAtomicNum() == 6 &&
    // RDKit❗❌:           bnd->getEndAtom()->getAtomicNum() == 6) {
    // RDKit❗❌:         bothCarbon.set(bnd->getIdx());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   for (const auto &bring : ringInfo->bondRings()) {
    // RDKit❗❌:     bool allC = true;
    // RDKit❗❌:     bool allAromatic = true;
    // RDKit❗❌:     for (const auto bidx : bring) {
    // RDKit❗❌:       if (!isArom[bidx]) {
    // RDKit❗❌:         allAromatic = false;
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (!bothCarbon[bidx]) {
    // RDKit❗❌:         allC = false;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (allAromatic) {
    // RDKit❗❌:       score += 100;
    // RDKit❗❌:       if (allC) {
    // RDKit❗❌:         score += 150;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return score;
    // RDKit❗❌: };
    // END RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::scoreRings restricted cache loan

    if !molecule.rings.is_symm_sssr() {
        *molecule.rings = symmetrized_sssr(molecule.topology, &Default::default())?;
    }
    score_tautomer_rings_retained_cache(molecule.as_record_view())
}

/// Score every match of the supplied source-ordered tautomer SMARTS terms.
#[must_use]
pub fn score_tautomer_substructures(
    molecule: TautomerRecordView<'_>,
    terms: &[TautomerScoreTerm],
) -> Result<i32, TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::scoreSubstructs
    // RDKit❗❌: int scoreSubstructs(const ROMol &mol,
    // RDKit❗❌:                     const std::vector<SubstructTerm> &substructureTerms) {
    // RDKit❗❌:   int score = 0;
    // RDKit❗❌:   SubstructMatchParameters params;
    // RDKit❗❌:   for (const auto &term : substructureTerms) {
    // RDKit❗❌:     if (!term.matcher.getNumAtoms()) {
    // RDKit❗❌:       BOOST_LOG(rdErrorLog) << " matcher for term " << term.name
    // RDKit❗❌:                             << " is invalid, ignoring it." << std::endl;
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     const auto nMatches = SubstructMatchCount(mol, term.matcher, params);
    // RDKit❗❌:     score += static_cast<int>(nMatches) * term.score;
    // RDKit❗❌:   }
    // RDKit❗❌:   return score;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::scoreSubstructs
    let mut score = 0_i32;
    for term in terms {
        let Some(pattern) = term.matcher().filter(|q| q.num_atoms() > 0) else {
            continue;
        };
        let plan = CompiledQuery::compile(pattern.clone())?;
        let count = crate::engine::query_match_count(molecule, &plan, &Default::default())?;
        score = score.wrapping_add((count as i32).wrapping_mul(term.score));
    }
    Ok(score)
}

/// Penalize hydrogens attached to phosphorus, sulfur, selenium, and tellurium.
pub fn score_tautomer_hetero_hydrogens(
    molecule: TautomerRecordView<'_>,
) -> Result<i32, TautomerRunError> {
    // RDKit✔️✔️: int scoreHeteroHs(const ROMol &mol) {
    // RDKit✔️✔️:   int score = 0;
    let mut score = 0_i32;
    // RDKit✔️✔️:   for (const auto &at : mol.atoms()) {
    for atom in &molecule.topology.atoms {
        // RDKit✔️✔️:     int anum = at->getAtomicNum();
        let atomic_number = atom.atomic_number();
        // RDKit✔️✔️:     if (anum == 15 || anum == 16 || anum == 34 || anum == 52) {
        if matches!(atomic_number, 15 | 16 | 34 | 52) {
            // RDKit✔️✔️:       score -= at->getTotalNumHs();
            let total_hydrogens = total_hydrogens(molecule, atom.id())?;
            score = score.wrapping_sub(total_hydrogens as i32);
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return score;
    Ok(score)
    // RDKit✔️✔️: }
}

/// Compute all source tautomer score components and their signed aggregate.
pub fn score_tautomer_with_terms_(
    molecule: &mut crate::TautomerScoreView<'_>,
    terms: &[TautomerScoreTerm],
) -> Result<TautomerScore, TautomerRunError> {
    // RDKit✔️✔️: inline int scoreTautomer(const ROMol &mol) {
    // RDKit✔️✔️:   return scoreRings(mol) + scoreSubstructs(mol) + scoreHeteroHs(mol);
    // RDKit✔️✔️: }
    let ring = score_tautomer_rings_(molecule)?;
    score_after_ring(molecule.as_record_view(), terms, ring)
}

/// Compute the source tautomer score using the pinned default SMARTS terms.
pub fn score_tautomer_(
    molecule: &mut crate::TautomerScoreView<'_>,
) -> Result<TautomerScore, TautomerRunError> {
    score_tautomer_with_terms_(molecule, default_tautomer_score_terms())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::query_matches;
    #[derive(Clone)]
    struct Fixture {
        record: cosmolkit_smiles::SmilesRecord,
        valence: Option<cosmolkit_core::ValenceAssignment>,
        rings: Option<cosmolkit_core::RingInfo>,
    }
    impl Fixture {
        fn view(&self) -> TautomerRecordView<'_> {
            TautomerRecordView {
                topology: &self.record.topology,
                coordinates: &self.record.coordinates,
                properties: &self.record.properties,
                valence: self.valence.as_ref(),
                rings: self.rings.as_ref(),
            }
        }
    }
    impl Fixture {
        fn score_view(&mut self) -> crate::TautomerScoreView<'_> {
            let rings = self.rings.get_or_insert_with(|| {
                cosmolkit_core::RingInfo::new(
                    cosmolkit_core::RingFindType::OtherOrUnknown,
                    self.record.topology.atoms.len(),
                    self.record.topology.bonds.len(),
                )
            });
            crate::TautomerScoreView::new(
                TautomerRecordView {
                    topology: &self.record.topology,
                    coordinates: &self.record.coordinates,
                    properties: &self.record.properties,
                    valence: self.valence.as_ref(),
                    rings: None,
                },
                rings,
            )
        }
    }

    fn fixture(text: &str) -> Result<Fixture, TautomerRunError> {
        fixture_with_sanitize(text, true)
    }
    fn fixture_with_sanitize(text: &str, sanitize: bool) -> Result<Fixture, TautomerRunError> {
        let record = cosmolkit_smiles::parse_smiles(
            text,
            &cosmolkit_smiles::SmilesParseParams {
                sanitize,
                ..Default::default()
            },
        )?;
        let valence = if sanitize {
            Some(cosmolkit_core::assign_valence(
                &record.topology,
                &cosmolkit_core::ValenceParams {
                    strict: false,
                    ..Default::default()
                },
            )?)
        } else {
            None
        };
        Ok(Fixture {
            record,
            valence,
            rings: None,
        })
    }
    fn query(text: &str) -> QueryGraph {
        parse_smarts(text, &Default::default()).expect("parse query fixture")
    }
    #[test]
    fn score_terms_match_all_twelve_source_rows_in_order() {
        let expected = [
            (
                "benzoquinone",
                "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]",
                25,
            ),
            ("oxim", "[#6]=[N][OH]", 4),
            ("C=O", "[#6]=,:[#8]", 2),
            ("N=O", "[#7]=,:[#8]", 2),
            ("P=O", "[#15]=,:[#8]", 2),
            ("C=hetero", "[C]=[!#1;!#6]", 1),
            ("C(=hetero)-hetero", "[C](=[!#1;!#6])[!#1;!#6]", 2),
            ("aromatic C = exocyclic N", "[c]=!@[N]", -1),
            ("methyl", "[CX4H3]", 1),
            ("guanidine terminal=N", "[#7]C(=[NR0])[#7H0]", 1),
            ("guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2),
            ("aci-nitro", "[#6]=[N+]([O-])[OH]", -4),
        ];
        let terms = default_tautomer_score_terms();

        assert_eq!(terms.len(), expected.len());
        for (index, (term, &(name, smarts, score))) in terms.iter().zip(&expected).enumerate() {
            assert_eq!(term.name(), name, "name at source row {index}");
            assert_eq!(term.smarts(), smarts, "SMARTS at source row {index}");
            assert_eq!(term.score(), score, "score at source row {index}");
            let expected_matcher = query(smarts);
            assert_eq!(
                term.matcher().expect("built-in matcher").atoms(),
                expected_matcher.atoms(),
                "matcher atoms at source row {index}"
            );
            assert_eq!(
                term.matcher().expect("built-in matcher").bonds(),
                expected_matcher.bonds(),
                "matcher bonds at source row {index}"
            );
        }
    }

    #[test]
    fn score_terms_compile_the_default_table_only_once() {
        let first = default_tautomer_score_terms();
        let second = default_tautomer_score_terms();

        assert!(std::ptr::eq(first, second));
        assert!(std::ptr::eq(
            first[0].matcher().expect("matcher").atoms().as_ptr(),
            second[0].matcher().expect("matcher").atoms().as_ptr()
        ));
    }

    #[test]
    fn score_terms_invalid_smarts_keep_metadata_and_an_empty_matcher() {
        let invalid = TautomerScoreTerm::new("invalid", "[", -7);

        assert_eq!(invalid.name(), "invalid");
        assert_eq!(invalid.smarts(), "[");
        assert_eq!(invalid.score(), -7);
        assert!(invalid.matcher().is_none());
    }

    #[test]
    fn score_terms_equality_uses_only_source_metadata() {
        let compiled = TautomerScoreTerm::new("same", "[#6]=[#8]", 2);
        let mut without_matcher = compiled.clone();
        without_matcher.matcher = None;

        assert_eq!(compiled, without_matcher);
        assert_ne!(
            compiled,
            TautomerScoreTerm::new("different", "[#6]=[#8]", 2)
        );
        assert_ne!(compiled, TautomerScoreTerm::new("same", "[#7]=[#8]", 2));
        assert_ne!(compiled, TautomerScoreTerm::new("same", "[#6]=[#8]", -2));
    }

    #[test]
    fn ring_score_is_zero_without_rings() {
        let mut molecule = fixture("CCO").expect("parse acyclic molecule");

        assert_eq!(
            score_tautomer_rings_(&mut molecule.score_view()).expect("score rings"),
            0
        );
    }

    #[test]
    fn ring_score_counts_an_aromatic_heteroring_without_the_all_carbon_bonus() {
        let mut molecule = fixture("c1ccncc1").expect("parse pyridine");

        assert_eq!(
            score_tautomer_rings_(&mut molecule.score_view()).expect("score rings"),
            100
        );
    }

    #[test]
    fn ring_score_adds_both_aromatic_and_all_carbon_weights() {
        let mut molecule = fixture("c1ccccc1").expect("parse benzene");

        assert_eq!(
            score_tautomer_rings_(&mut molecule.score_view()).expect("score rings"),
            250
        );
    }

    #[test]
    fn ring_score_traverses_every_symmetrized_fused_ring() {
        let mut molecule = fixture("c1ccc2ccccc2c1").expect("parse naphthalene");

        assert_eq!(
            symmetrized_sssr(&molecule.record.topology, &Default::default())
                .expect("find fused rings")
                .bond_rings()
                .len(),
            2
        );
        assert_eq!(
            score_tautomer_rings_(&mut molecule.score_view()).expect("score rings"),
            500
        );
    }

    #[test]
    fn ring_score_ignores_nonaromatic_rings() {
        let mut molecule = fixture("C1CCCCC1").expect("parse cyclohexane");

        assert_eq!(
            score_tautomer_rings_(&mut molecule.score_view()).expect("score rings"),
            0
        );
    }

    #[test]
    fn ring_score_matches_with_or_without_a_precomputed_symm_sssr_cache() {
        let mut without_cache =
            fixture_with_sanitize("c1ccncc1", false).expect("parse uncached pyridine");
        assert!(without_cache.rings.is_none());
        let mut with_cache = without_cache.clone();
        with_cache.rings = Some(
            symmetrized_sssr(&with_cache.record.topology, &Default::default())
                .expect("assign symmetric SSSR cache"),
        );
        assert!(
            with_cache
                .rings
                .as_ref()
                .is_some_and(cosmolkit_core::RingInfo::is_symm_sssr)
        );

        assert_eq!(
            score_tautomer_rings_(&mut without_cache.score_view()).expect("score uncached"),
            100
        );
        assert_eq!(
            score_tautomer_rings_(&mut with_cache.score_view()).expect("score cached"),
            100
        );
    }

    #[test]
    fn ring_score_keeps_topology_and_retains_source_symm_cache() {
        let mut molecule =
            fixture_with_sanitize("c1ccccc1", false).expect("parse uncached benzene");
        let topology_before = molecule.record.topology.clone();
        assert!(molecule.rings.is_none());

        assert_eq!(
            score_tautomer_rings_(&mut molecule.score_view()).expect("score rings"),
            250
        );
        assert_eq!(&molecule.record.topology, &topology_before);
        assert!(molecule.rings.as_ref().unwrap().is_symm_sssr());
    }

    #[test]
    fn substructure_score_is_zero_when_no_term_matches() {
        let mut molecule = fixture("CCC").expect("parse propane");
        let terms = [TautomerScoreTerm::new("oxygen", "[#8]", 7)];

        assert_eq!(
            score_tautomer_substructures(molecule.view(), &terms).expect("score substructures"),
            0
        );
    }

    #[test]
    fn substructure_score_counts_one_and_multiple_matches() {
        let ethanol = fixture("CCO").expect("parse ethanol");
        let oxygen = [TautomerScoreTerm::new("oxygen", "[#8]", 7)];
        let carbon = [TautomerScoreTerm::new("carbon", "[#6]", 2)];

        assert_eq!(
            score_tautomer_substructures(ethanol.view(), &oxygen).expect("score substructures"),
            7
        );
        assert_eq!(
            score_tautomer_substructures(ethanol.view(), &carbon).expect("score substructures"),
            4
        );
    }

    #[test]
    fn substructure_score_counts_distinct_overlapping_matches() {
        let propane = fixture("CCC").expect("parse propane");
        let terms = [TautomerScoreTerm::new("carbon bond", "[#6]-[#6]", 3)];

        assert_eq!(
            score_tautomer_substructures(propane.view(), &terms).expect("score substructures"),
            6
        );
    }

    #[test]
    fn substructure_score_evaluates_every_builtin_term_with_shared_match_semantics() {
        let mut molecule =
            fixture("CC(=O)NC(=N)N.c1ccccc1").expect("parse representative scoring molecule");

        for term in default_tautomer_score_terms() {
            let matcher = term.matcher().expect("built-in matcher compiles");
            let expected = (query_matches(
                molecule.view(),
                &CompiledQuery::compile(matcher.clone()).unwrap(),
            )
            .unwrap()
            .len() as i32)
                .wrapping_mul(term.score());
            assert_eq!(
                score_tautomer_substructures(molecule.view(), std::slice::from_ref(term))
                    .expect("score substructures"),
                expected,
                "built-in term {}",
                term.name()
            );
        }
    }

    #[test]
    fn substructure_score_accumulates_custom_positive_and_negative_terms_in_source_order() {
        let mut molecule = fixture("CCO").expect("parse ethanol");
        let terms = [
            TautomerScoreTerm::new("carbon", "[#6]", 3),
            TautomerScoreTerm::new("oxygen", "[#8]", -5),
            TautomerScoreTerm::new("no match", "[#15]", 100),
        ];

        assert_eq!(
            score_tautomer_substructures(molecule.view(), &terms).expect("score substructures"),
            1
        );
        assert_eq!(
            terms
                .iter()
                .map(TautomerScoreTerm::name)
                .collect::<Vec<_>>(),
            ["carbon", "oxygen", "no match"]
        );
    }

    #[test]
    fn substructure_score_skips_invalid_matchers_without_dropping_other_terms() {
        let mut molecule = fixture("CCO").expect("parse ethanol");
        let terms = [
            TautomerScoreTerm::new("invalid", "[", i32::MAX),
            TautomerScoreTerm::new("oxygen", "[#8]", -5),
        ];

        assert_eq!(
            score_tautomer_substructures(molecule.view(), &terms).expect("score substructures"),
            -5
        );
    }

    #[test]
    fn substructure_score_retains_source_int_accumulation_boundaries() {
        let methane = fixture("C").expect("parse methane");
        let terms = [
            TautomerScoreTerm::new("maximum", "[#6]", i32::MAX),
            TautomerScoreTerm::new("one", "[#6]", 1),
        ];

        assert_eq!(
            score_tautomer_substructures(methane.view(), &terms).expect("score substructures"),
            i32::MIN
        );
    }

    #[test]
    fn aggregate_score_penalizes_implicit_hydrogens_on_each_source_element() {
        for (smiles, expected) in [("P", -3), ("S", -2), ("[SeH2]", -2), ("[TeH3]", -3)] {
            let mut molecule = fixture(smiles).expect("parse heteroatom");
            assert_eq!(
                score_tautomer_hetero_hydrogens(molecule.view()).expect("score hetero hydrogens"),
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn aggregate_score_counts_bracket_explicit_but_not_neighbor_isotopic_hydrogen() {
        let explicit = fixture("[PH3]").expect("parse explicit phosphorus Hs");
        let isotopic_neighbor = fixture("[2H]P").expect("parse isotopic hydrogen neighbor");

        assert_eq!(
            score_tautomer_hetero_hydrogens(explicit.view()).unwrap(),
            -3
        );
        assert_eq!(
            score_tautomer_hetero_hydrogens(isotopic_neighbor.view()).unwrap(),
            -2
        );
    }

    #[test]
    fn aggregate_score_counts_hydrogens_on_charged_penalized_atoms() {
        for (smiles, expected) in [("[PH4+]", -4), ("[SH-]", -1)] {
            let mut molecule = fixture(smiles).expect("parse charged heteroatom");
            assert_eq!(
                score_tautomer_hetero_hydrogens(molecule.view()).expect("score charged atom"),
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn aggregate_score_excludes_hydrogen_bearing_elements_outside_the_source_set() {
        for smiles in ["C", "N", "O", "F", "Cl"] {
            let mut molecule = fixture(smiles).expect("parse excluded element");
            assert_eq!(
                score_tautomer_hetero_hydrogens(molecule.view()).expect("score excluded element"),
                0,
                "{smiles}"
            );
        }
    }

    #[test]
    fn aggregate_score_reports_missing_valence_only_when_a_penalized_atom_needs_it() {
        let phosphorus = fixture_with_sanitize("P", false).expect("parse unsanitized phosphorus");
        let oxygen = fixture_with_sanitize("O", false).expect("parse unsanitized oxygen");

        assert!(matches!(score_tautomer_hetero_hydrogens(phosphorus.view()),
            Err(TautomerRunError::Valence(cosmolkit_core::ValenceError::ImplicitValenceCacheNotInitialized {atom})) if atom==cosmolkit_model::AtomId::new(0)));
        assert_eq!(score_tautomer_hetero_hydrogens(oxygen.view()).unwrap(), 0);
    }

    #[test]
    fn aggregate_score_exposes_exact_components_and_signed_sum() {
        let mut benzene = fixture("c1ccccc1").expect("parse benzene");
        let carbon_terms = [TautomerScoreTerm::new("carbon", "[#6]", 2)];
        let score = score_tautomer_with_terms_(&mut benzene.score_view(), &carbon_terms)
            .expect("score benzene");

        assert_eq!(score.ring(), 250);
        assert_eq!(score.substructure(), 12);
        assert_eq!(score.hetero_hydrogen(), 0);
        assert_eq!(score.total(), 262);
    }

    #[test]
    fn aggregate_score_default_delegate_uses_the_single_builtin_term_table() {
        let mut molecule = fixture("CC(=O)N").expect("parse acetamide");
        let delegated = score_tautomer_(&mut molecule.score_view()).expect("score with defaults");
        let explicit =
            score_tautomer_with_terms_(&mut molecule.score_view(), default_tautomer_score_terms())
                .expect("score with explicit defaults");

        assert_eq!(delegated, explicit);
        assert_eq!(
            delegated.total(),
            delegated
                .ring()
                .wrapping_add(delegated.substructure())
                .wrapping_add(delegated.hetero_hydrogen())
        );
    }
}

fn relevant_tautomer_score_indices(
    molecule: TautomerRecordView<'_>,
    terms: &[TautomerScoreTerm],
) -> Result<Vec<usize>, TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION getRelevantSubstructTermIndices
    // RDKit❗❌: std::vector<size_t> getRelevantSubstructTermIndices(
    // RDKit❗❌:     const ROMol &mol, const std::vector<SubstructTerm> &terms) {
    // RDKit❗❌:   // Prepare relevant SubstructTerms for this molecule by filtering in two
    // RDKit❗❌:   // stages:
    // RDKit❗❌:   //   1. Element check: skip terms requiring elements not in the molecule
    // RDKit❗❌:   //   2. Connectivity check: skip terms whose bond-order-agnostic pattern
    // RDKit❗❌:   //      doesn't match (since tautomerization doesn't create/destroy bonds)
    // RDKit❗❌:
    // RDKit❗❌:   std::unordered_set<int> presentElements;
    // RDKit❗❌:   for (const auto atom : mol.atoms()) {
    // RDKit❗❌:     presentElements.insert(atom->getAtomicNum());
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<size_t> relevantIndices;
    // RDKit❗❌:   relevantIndices.reserve(terms.size());
    // RDKit❗❌:
    // RDKit❗❌:   SubstructMatchParameters params;
    // RDKit❗❌:   params.maxMatches = 1;
    // RDKit❗❌:
    // RDKit❗❌:   for (size_t i = 0; i < terms.size(); ++i) {
    // RDKit❗❌:     const auto &term = terms[i];
    // RDKit❗❌:
    // RDKit❗❌:     bool hasAllElements = true;
    // RDKit❗❌:     for (int elem : term.requiredElements) {
    // RDKit❗❌:       if (presentElements.find(elem) == presentElements.end()) {
    // RDKit❗❌:         hasAllElements = false;
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!hasAllElements) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (term.connectivityMatcher.getNumAtoms() > 0) {
    // RDKit❗❌:       auto matches = SubstructMatch(mol, term.connectivityMatcher, params);
    // RDKit❗❌:       if (matches.empty()) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     relevantIndices.push_back(i);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return relevantIndices;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION getRelevantSubstructTermIndices
    // Empty construction and source-order inserts grow with distinct elements,
    // rather than reserving the entire atom iterator's exact size hint.
    let mut present = std::collections::HashSet::new();
    for atom in &molecule.topology.atoms {
        present.insert(atom.atomic_number());
    }
    let mut relevant = Vec::with_capacity(terms.len());
    let params = cosmolkit_search::SubstructMatchParams {
        max_matches: 1,
        ..Default::default()
    };
    for (i, term) in terms.iter().enumerate() {
        if !term.required_elements.iter().all(|e| present.contains(e)) {
            continue;
        }
        if let Some(pattern) = term
            .connectivity_matcher
            .as_ref()
            .filter(|q| q.num_atoms() > 0)
        {
            let compiled = CompiledQuery::compile(pattern.clone())?;
            if crate::engine::query_matches_with_params(molecule, &compiled, &params)?.is_empty() {
                continue;
            }
        }
        relevant.push(i);
    }
    Ok(relevant)
}
fn count_double_or_aromatic_tautomer_bonds(
    molecule: TautomerRecordView<'_>,
    elem1: u8,
    elem2: u8,
) -> u32 {
    // BEGIN RDKIT CPP FUNCTION countDoubleOrAromaticBonds
    // RDKit❗✔️: inline unsigned int countDoubleOrAromaticBonds(const ROMol &mol, int elem1,
    // RDKit❗✔️:                                                int elem2) {
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:     const auto bt = bond->getBondType();
    // RDKit❗✔️:     if (bt != Bond::DOUBLE && bt != Bond::AROMATIC) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     const int a1 = bond->getBeginAtom()->getAtomicNum();
    // RDKit❗✔️:     const int a2 = bond->getEndAtom()->getAtomicNum();
    // RDKit❗✔️:     if ((a1 == elem1 && a2 == elem2) || (a1 == elem2 && a2 == elem1)) {
    // RDKit❗✔️:       ++count;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return count;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION countDoubleOrAromaticBonds
    let mut count = 0_u32;
    for b in &molecule.topology.bonds {
        if !matches!(
            b.order(),
            cosmolkit_types::BondOrder::Double | cosmolkit_types::BondOrder::Aromatic
        ) {
            continue;
        }
        let a = molecule.topology.atoms[b.begin().index()].atomic_number();
        let c = molecule.topology.atoms[b.end().index()].atomic_number();
        if (a == elem1 && c == elem2) || (a == elem2 && c == elem1) {
            count = count.wrapping_add(1);
        }
    }
    count
}
fn count_tautomer_methyls(molecule: TautomerRecordView<'_>) -> Result<u32, TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION countMethyls
    // RDKit❗✔️: inline unsigned int countMethyls(const ROMol &mol) {
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   for (const auto atom : mol.atoms()) {
    // RDKit❗✔️:     if (atom->getAtomicNum() != 6) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // X4 means total degree 4 (including implicit H)
    // RDKit❗✔️:     if (atom->getTotalDegree() != 4) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // H3 means exactly 3 total hydrogens
    // RDKit❗✔️:     if (atom->getTotalNumHs() != 3) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ++count;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return count;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION countMethyls
    let mut count = 0_u32;
    for a in &molecule.topology.atoms {
        if a.atomic_number() != 6 {
            continue;
        }
        let hydrogens = total_hydrogens(molecule, a.id())?;
        if molecule
            .topology
            .adjacency
            .neighbors_of(a.id().index())
            .len() as u32
            + hydrogens
            != 4
        {
            continue;
        }
        if hydrogens != 3 {
            continue;
        }
        count = count.wrapping_add(1);
    }
    Ok(count)
}
fn count_tautomer_carbon_double_hetero(molecule: TautomerRecordView<'_>) -> u32 {
    // BEGIN RDKIT CPP FUNCTION countCarbonDoubleHetero
    // RDKit❗✔️: inline unsigned int countCarbonDoubleHetero(const ROMol &mol) {
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:     if (bond->getBondType() != Bond::DOUBLE) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     const auto *a1 = bond->getBeginAtom();
    // RDKit❗✔️:     const auto *a2 = bond->getEndAtom();
    // RDKit❗✔️:     // [C] is aliphatic carbon (not aromatic)
    // RDKit❗✔️:     // Check C double-bonded to heteroatom (not H, not C)
    // RDKit❗✔️:     auto isAliphaticCarbonToHetero = [](const Atom *c, const Atom *het) {
    // RDKit❗✔️:       return c->getAtomicNum() == 6 && !c->getIsAromatic() &&
    // RDKit❗✔️:              het->getAtomicNum() != 1 && het->getAtomicNum() != 6;
    // RDKit❗✔️:     };
    // RDKit❗✔️:     if (isAliphaticCarbonToHetero(a1, a2) ||
    // RDKit❗✔️:         isAliphaticCarbonToHetero(a2, a1)) {
    // RDKit❗✔️:       ++count;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return count;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION countCarbonDoubleHetero
    let mut count = 0_u32;
    for b in &molecule.topology.bonds {
        if b.order() != cosmolkit_types::BondOrder::Double {
            continue;
        }
        let a = &molecule.topology.atoms[b.begin().index()];
        let c = &molecule.topology.atoms[b.end().index()];
        let matches = |carbon: &cosmolkit_model::Atom, het: &cosmolkit_model::Atom| {
            carbon.atomic_number() == 6
                && !carbon.is_aromatic()
                && !matches!(het.atomic_number(), 1 | 6)
        };
        if matches(a, c) || matches(c, a) {
            count = count.wrapping_add(1);
        }
    }
    count
}
fn count_tautomer_aromatic_exocyclic_n(
    molecule: TautomerRecordView<'_>,
) -> Result<u32, TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION countAromaticCarbonExocyclicN
    // RDKit❗✔️: inline unsigned int countAromaticCarbonExocyclicN(const ROMol &mol) {
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   const auto *ringInfo = mol.getRingInfo();
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:     if (bond->getBondType() != Bond::DOUBLE) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // !@ means bond not in ring
    // RDKit❗✔️:     if (ringInfo->numBondRings(bond->getIdx()) > 0) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     const auto *a1 = bond->getBeginAtom();
    // RDKit❗✔️:     const auto *a2 = bond->getEndAtom();
    // RDKit❗✔️:     // [c] is aromatic carbon, [N] is any nitrogen
    // RDKit❗✔️:     auto isAromaticCarbonToN = [](const Atom *c, const Atom *n) {
    // RDKit❗✔️:       return c->getAtomicNum() == 6 && c->getIsAromatic() &&
    // RDKit❗✔️:              n->getAtomicNum() == 7;
    // RDKit❗✔️:     };
    // RDKit❗✔️:     if (isAromaticCarbonToN(a1, a2) || isAromaticCarbonToN(a2, a1)) {
    // RDKit❗✔️:       ++count;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return count;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION countAromaticCarbonExocyclicN
    let mut count = 0_u32;
    for b in &molecule.topology.bonds {
        if b.order() != cosmolkit_types::BondOrder::Double {
            continue;
        }
        // Source numBondRings requires initialized RingInfo. Normal optimized
        // caller has already scored rings and prepared candidate cache; raw
        // detached unsupported precondition stays a typed error, no fallback.
        let rings =
            molecule
                .rings
                .filter(|r| r.is_initialized())
                .ok_or(RingFindingError::Value {
                    message: "RingInfo not initialized",
                })?;
        if rings.num_bond_rings(b.id()) > 0 {
            continue;
        }
        let a = &molecule.topology.atoms[b.begin().index()];
        let c = &molecule.topology.atoms[b.end().index()];
        let matches = |carbon: &cosmolkit_model::Atom, nitrogen: &cosmolkit_model::Atom| {
            carbon.atomic_number() == 6 && carbon.is_aromatic() && nitrogen.atomic_number() == 7
        };
        if matches(a, c) || matches(c, a) {
            count = count.wrapping_add(1);
        }
    }
    Ok(count)
}
fn score_tautomer_substructures_filtered(
    molecule: TautomerRecordView<'_>,
    terms: &[TautomerScoreTerm],
    indices: &[usize],
) -> Result<i32, TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION scoreSubstructsFiltered
    // RDKit❗❌: int scoreSubstructsFiltered(const ROMol &mol,
    // RDKit❗❌:                             const std::vector<SubstructTerm> &terms,
    // RDKit❗❌:                             const std::vector<size_t> &relevantIndices) {
    // RDKit❗❌:   int score = 0;
    // RDKit❗❌:   SubstructMatchParameters params;
    // RDKit❗❌:
    // RDKit❗❌:   for (const size_t idx : relevantIndices) {
    // RDKit❗❌:     const auto &term = terms[idx];
    // RDKit❗❌:     if (!term.matcher.getNumAtoms()) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     unsigned int nMatches = 0;
    // RDKit❗❌:
    // RDKit❗❌:     // Use specialized matchers for simple patterns, VF2 for complex ones
    // RDKit❗❌:     switch (static_cast<PatternIdx>(idx)) {
    // RDKit❗❌:       case PatternIdx::CarbonylO:
    // RDKit❗❌:         nMatches = countDoubleOrAromaticBonds(mol, 6, 8);  // C=O
    // RDKit❗❌:         break;
    // RDKit❗❌:       case PatternIdx::NO:
    // RDKit❗❌:         nMatches = countDoubleOrAromaticBonds(mol, 7, 8);  // N=O
    // RDKit❗❌:         break;
    // RDKit❗❌:       case PatternIdx::PO:
    // RDKit❗❌:         nMatches = countDoubleOrAromaticBonds(mol, 15, 8);  // P=O
    // RDKit❗❌:         break;
    // RDKit❗❌:       case PatternIdx::CHetero:
    // RDKit❗❌:         nMatches = countCarbonDoubleHetero(mol);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case PatternIdx::AromaticCN:
    // RDKit❗❌:         nMatches = countAromaticCarbonExocyclicN(mol);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case PatternIdx::Methyl:
    // RDKit❗❌:         nMatches = countMethyls(mol);
    // RDKit❗❌:         break;
    // RDKit❗❌:       default:
    // RDKit❗❌:         // Fall back to VF2 for complex patterns (benzoquinone, oxim,
    // RDKit❗❌:         // C(=hetero)-hetero, guanidine, aci-nitro)
    // RDKit❗❌:         nMatches = SubstructMatchCount(mol, term.matcher, params);
    // RDKit❗❌:         break;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     score += static_cast<int>(nMatches) * term.score;
    // RDKit❗❌:   }
    // RDKit❗❌:   return score;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION scoreSubstructsFiltered
    let mut score = 0_i32;
    for &index in indices {
        let term = &terms[index];
        let Some(pattern) = term.matcher.as_ref().filter(|q| q.num_atoms() > 0) else {
            continue;
        };
        let count = match index {
            2 => count_double_or_aromatic_tautomer_bonds(molecule, 6, 8),
            3 => count_double_or_aromatic_tautomer_bonds(molecule, 7, 8),
            4 => count_double_or_aromatic_tautomer_bonds(molecule, 15, 8),
            5 => count_tautomer_carbon_double_hetero(molecule),
            7 => count_tautomer_aromatic_exocyclic_n(molecule)?,
            8 => count_tautomer_methyls(molecule)?,
            _ => crate::engine::query_match_count(
                molecule,
                &CompiledQuery::compile(pattern.clone())?,
                &Default::default(),
            )?,
        };
        score = score.wrapping_add((count as i32).wrapping_mul(term.score));
    }
    Ok(score)
}
pub(crate) struct OptimizedTautomerScorer {
    indices: Vec<usize>,
}
pub(crate) fn make_optimized_tautomer_scorer(
    source: TautomerRecordView<'_>,
) -> Result<OptimizedTautomerScorer, TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::makeOptimizedScorer
    // RDKit❗❌: inline boost::function<int(const ROMol &)> makeOptimizedScorer(
    // RDKit❗❌:     const ROMol &mol) {
    // RDKit❗❌:   auto relevantIndices = getRelevantSubstructTermIndices(mol);
    // RDKit❗❌:   // Capture by value since the indices are small and we want the lambda to
    // RDKit❗❌:   // outlive this function.
    // RDKit❗❌:   return [relevantIndices](const ROMol &taut) {
    // RDKit❗❌:     const auto &terms = getDefaultTautomerScoreSubstructs();
    // RDKit❗❌:     return scoreRings(taut) +
    // RDKit❗❌:            scoreSubstructsFiltered(taut, terms, relevantIndices) +
    // RDKit❗❌:            scoreHeteroHs(taut);
    // RDKit❗❌:   };
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::makeOptimizedScorer
    #[cfg(test)]
    SEARCH03_FACTORY_COUNT.with(|n| n.set(n.get() + 1));
    Ok(OptimizedTautomerScorer {
        indices: relevant_tautomer_score_indices(source, default_tautomer_score_terms())?,
    })
}
impl OptimizedTautomerScorer {
    pub(crate) fn score(
        &self,
        molecule: &mut crate::TautomerScoreView<'_>,
    ) -> Result<i32, TautomerRunError> {
        // BEGIN RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::makeOptimizedScorer::lambda
        // RDKit❗❌:   return [relevantIndices](const ROMol &taut) {
        // RDKit❗❌:     const auto &terms = getDefaultTautomerScoreSubstructs();
        // RDKit❗❌:     return scoreRings(taut) +
        // RDKit❗❌:            scoreSubstructsFiltered(taut, terms, relevantIndices) +
        // RDKit❗❌:            scoreHeteroHs(taut);
        // RDKit❗❌:   };
        // END RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::makeOptimizedScorer::lambda
        let ring = score_tautomer_rings_(molecule)?;
        let subs = score_tautomer_substructures_filtered(
            molecule.as_record_view(),
            default_tautomer_score_terms(),
            &self.indices,
        )?;
        Ok(ring
            .wrapping_add(subs)
            .wrapping_add(score_tautomer_hetero_hydrogens(molecule.as_record_view())?))
    }
}
#[cfg(test)]
thread_local! {static SEARCH03_FACTORY_COUNT:std::cell::Cell<usize>=const {std::cell::Cell::new(0)};}

#[cfg(test)]
mod search03_tests {
    use super::*;
    use cosmolkit_model::CoordinateBlock;
    fn molecule(smiles: &str) -> crate::TautomerRecord {
        crate::engine::stereo_tests::fixture_from_smiles(smiles).unwrap()
    }
    #[test]
    fn search03_default_prerequisites_and_equality_keep_source_order() {
        let terms = default_tautomer_score_terms();
        assert_eq!(terms.len(), 12);
        assert_eq!(terms[1].required_elements, vec![6, 7, 8]);
        assert_eq!(terms[5].connectivity_smarts, "[C]~[!#1;!#6]");
        assert!(terms[8].connectivity_matcher.is_none());
        assert_eq!(
            terms[5],
            TautomerScoreTerm::new(terms[5].name(), terms[5].smarts(), terms[5].score())
        );
        let m = molecule("[He]");
        assert!(
            relevant_tautomer_score_indices(m.view(&CoordinateBlock::default()), terms)
                .unwrap()
                .is_empty()
        );
    }
    #[test]
    fn search03_input_prefilter_retains_aromatic_connectivity_quirk() {
        let input = molecule("Oc1ccccn1");
        let candidate = molecule("O=C1CC=CC=N1");
        let coords = CoordinateBlock::default();
        let terms = default_tautomer_score_terms();
        let indices = relevant_tautomer_score_indices(input.view(&coords), terms).unwrap();
        assert!(!indices.contains(&5));
        assert!(
            relevant_tautomer_score_indices(candidate.view(&coords), terms)
                .unwrap()
                .contains(&5)
        );
        assert_eq!(
            count_tautomer_carbon_double_hetero(candidate.view(&coords)),
            2
        );
        assert_eq!(
            score_tautomer_substructures_filtered(candidate.view(&coords), terms, &[5]).unwrap(),
            2
        );
    }
    #[test]
    fn search03_specialized_predicates_use_attributes_without_invented_gates() {
        let coords = CoordinateBlock::default();
        let mut carbonyl = molecule("C=O");
        assert_eq!(
            count_double_or_aromatic_tautomer_bonds(carbonyl.view(&coords), 6, 8),
            1
        );
        assert_eq!(
            count_tautomer_carbon_double_hetero(carbonyl.view(&coords)),
            1
        );
        carbonyl.topology.atoms[0].set_aromatic(true);
        assert_eq!(
            count_tautomer_carbon_double_hetero(carbonyl.view(&coords)),
            0
        );
        carbonyl.topology.bonds[0].set_order(cosmolkit_types::BondOrder::Aromatic);
        assert_eq!(
            count_double_or_aromatic_tautomer_bonds(carbonyl.view(&coords), 8, 6),
            1
        );
        let mut imine = molecule("C=N");
        imine.topology.atoms[0].set_aromatic(true);
        imine.topology.atoms[1].set_aromatic(true);
        assert_eq!(
            count_tautomer_aromatic_exocyclic_n(imine.view(&coords)).unwrap(),
            1
        );
        let mut ring = molecule("C1=NCCC1");
        ring.topology.atoms[0].set_aromatic(true);
        assert_eq!(
            count_tautomer_aromatic_exocyclic_n(ring.view(&coords)).unwrap(),
            0
        );
        for (smiles, expected) in [("CC", 2), ("C", 0), ("[CH3]", 0), ("[H]C([H])([H])C", 2)] {
            assert_eq!(
                count_tautomer_methyls(molecule(smiles).view(&coords)).unwrap(),
                expected,
                "{smiles}"
            );
        }
        let mut methyl = molecule("CC");
        methyl.topology.atoms[0].set_aromatic(true);
        assert_eq!(count_tautomer_methyls(methyl.view(&coords)).unwrap(), 2);
        let missing = TautomerRecordView {
            rings: None,
            ..imine.view(&coords)
        };
        assert!(matches!(
            count_tautomer_aromatic_exocyclic_n(missing),
            Err(TautomerRunError::Rings(RingFindingError::Value {
                message: "RingInfo not initialized"
            }))
        ));
    }
    #[test]
    fn search03_numeric_index_dispatch_is_uncapped_and_keeps_generic_count_cap() {
        let coords = CoordinateBlock::default();
        let m = molecule(&vec!["CC"; 501].join("."));
        let mut terms = (0..9)
            .map(|_| TautomerScoreTerm::new("custom", "[CX4H3]", 1))
            .collect::<Vec<_>>();
        assert_eq!(
            score_tautomer_substructures(m.view(&coords), &terms[8..9]).unwrap(),
            1000
        );
        assert_eq!(
            score_tautomer_substructures_filtered(m.view(&coords), &terms, &[8]).unwrap(),
            1002
        );
        terms[8] = TautomerScoreTerm::new("invalid", "[", 99);
        assert_eq!(
            score_tautomer_substructures_filtered(m.view(&coords), &terms, &[8]).unwrap(),
            0
        );
        let formaldehyde = molecule("C=O");
        terms[2] = TautomerScoreTerm::new("renamed", "[#6]", -3);
        assert_eq!(
            score_tautomer_substructures_filtered(formaldehyde.view(&coords), &terms, &[2])
                .unwrap(),
            -3
        );
    }
    #[test]
    fn search03_invalid_connectivity_retains_term_but_element_filter_runs_first() {
        let m = molecule("CC");
        let coords = CoordinateBlock::default();
        let mut term = TautomerScoreTerm::with_prerequisites(
            "invalid".to_owned(),
            "C".to_owned(),
            1,
            vec![6],
            "[".to_owned(),
        );
        assert_eq!(
            relevant_tautomer_score_indices(m.view(&coords), std::slice::from_ref(&term)).unwrap(),
            vec![0]
        );
        term.required_elements.push(8);
        assert!(
            relevant_tautomer_score_indices(m.view(&coords), &[term])
                .unwrap()
                .is_empty()
        );
    }
    #[test]
    fn search03_real_canonical_default_constructs_once_custom_never_constructs() {
        let coords = CoordinateBlock::default();
        let mut m = molecule("CC");
        let catalog = crate::TautomerCatalog::current().unwrap();
        SEARCH03_FACTORY_COUNT.with(|n| n.set(0));
        crate::canonicalize_with_catalog(
            m.score_view(&coords),
            &catalog,
            Default::default(),
            None,
            None,
        )
        .unwrap();
        assert_eq!(SEARCH03_FACTORY_COUNT.with(std::cell::Cell::get), 1);
        SEARCH03_FACTORY_COUNT.with(|n| n.set(0));
        let mut calls = 0;
        let mut custom = |_: crate::TautomerScoreView<'_>| {
            calls += 1;
            Ok(0)
        };
        crate::canonicalize_with_catalog(
            m.score_view(&coords),
            &catalog,
            Default::default(),
            None,
            Some(&mut custom),
        )
        .unwrap();
        assert_eq!(SEARCH03_FACTORY_COUNT.with(std::cell::Cell::get), 0);
        assert_eq!(calls, 0);
    }
}

impl OptimizedTautomerScorer {
    pub(crate) fn score_record(
        &self,
        record: &mut crate::TautomerRecord,
        coordinates: &cosmolkit_model::CoordinateBlock,
    ) -> Result<i32, TautomerRunError> {
        self.score(&mut record.score_view(coordinates))
    }
}

fn score_tautomer_rings_retained_cache(
    molecule: TautomerRecordView<'_>,
) -> Result<i32, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::scoreRings restricted cache loan
    // RDKit❗❌: int scoreRings(const ROMol &mol) {
    // RDKit❗❌:   int score = 0;
    // RDKit❗❌:   auto ringInfo = mol.getRingInfo();
    // RDKit❗❌:   if (!ringInfo->isSymmSssr()) {
    // RDKit❗❌:     MolOps::symmetrizeSSSR(const_cast<ROMol &>(mol));
    // RDKit❗❌:     ringInfo = mol.getRingInfo();
    // RDKit❗❌:   }
    // RDKit❗❌:   boost::dynamic_bitset<> isArom(mol.getNumBonds());
    // RDKit❗❌:   boost::dynamic_bitset<> bothCarbon(mol.getNumBonds());
    // RDKit❗❌:   for (const auto &bnd : mol.bonds()) {
    // RDKit❗❌:     if (bnd->getIsAromatic()) {
    // RDKit❗❌:       isArom.set(bnd->getIdx());
    // RDKit❗❌:       if (bnd->getBeginAtom()->getAtomicNum() == 6 &&
    // RDKit❗❌:           bnd->getEndAtom()->getAtomicNum() == 6) {
    // RDKit❗❌:         bothCarbon.set(bnd->getIdx());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   for (const auto &bring : ringInfo->bondRings()) {
    // RDKit❗❌:     bool allC = true;
    // RDKit❗❌:     bool allAromatic = true;
    // RDKit❗❌:     for (const auto bidx : bring) {
    // RDKit❗❌:       if (!isArom[bidx]) {
    // RDKit❗❌:         allAromatic = false;
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (!bothCarbon[bidx]) {
    // RDKit❗❌:         allC = false;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (allAromatic) {
    // RDKit❗❌:       score += 100;
    // RDKit❗❌:       if (allC) {
    // RDKit❗❌:         score += 150;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return score;
    // RDKit❗❌: };
    // END RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::scoreRings restricted cache loan
    let mut score = 0_i32;
    let ring_info =
        molecule
            .rings
            .filter(|rings| rings.is_symm_sssr())
            .ok_or(RingFindingError::Value {
                message: "tautomer score requires retained SymmSSSR cache",
            })?;
    let mut is_aromatic = vec![false; molecule.topology.bonds.len()];
    let mut both_carbon = vec![false; molecule.topology.bonds.len()];
    for bond in &molecule.topology.bonds {
        if bond.is_aromatic() {
            is_aromatic[bond.id().index()] = true;
            if molecule.topology.atoms[bond.begin().index()].atomic_number() == 6
                && molecule.topology.atoms[bond.end().index()].atomic_number() == 6
            {
                both_carbon[bond.id().index()] = true;
            }
        }
    }

    for bond_ring in ring_info.bond_rings() {
        let mut all_carbon = true;
        let mut all_aromatic = true;
        for bond_id in bond_ring {
            if !is_aromatic[bond_id.index()] {
                all_aromatic = false;
                break;
            }
            if !both_carbon[bond_id.index()] {
                all_carbon = false;
            }
        }
        if all_aromatic {
            score += 100;
            if all_carbon {
                score += 150;
            }
        }
    }
    Ok(score)
}

fn score_after_ring(
    molecule: TautomerRecordView<'_>,
    terms: &[TautomerScoreTerm],
    ring: i32,
) -> Result<TautomerScore, TautomerRunError> {
    // RDKit❗❌: inline int scoreTautomer(const ROMol &mol) {
    // RDKit❗❌:   return scoreRings(mol) + scoreSubstructs(mol) + scoreHeteroHs(mol);
    // RDKit❗❌: }
    Ok(TautomerScore {
        ring,
        substructure: score_tautomer_substructures(molecule, terms)?,
        hetero_hydrogen: score_tautomer_hetero_hydrogens(molecule)?,
    })
}
#[doc(hidden)]
pub fn score_tautomer_from_retained_cache(
    molecule: TautomerRecordView<'_>,
    terms: Option<&[TautomerScoreTerm]>,
) -> Result<TautomerScore, TautomerRunError> {
    let ring = score_tautomer_rings_retained_cache(molecule)?;
    score_after_ring(
        molecule,
        terms.unwrap_or_else(|| default_tautomer_score_terms()),
        ring,
    )
}

#[cfg(test)]
mod recovery_search06_cache_tests {
    use super::*;
    #[test]
    fn recovery_search06_ring_loan_cold_fast_warm_actual_cache() {
        let coordinates = cosmolkit_model::CoordinateBlock::default();
        let mut record =
            crate::engine::stereo_tests::fixture_from_smiles("C12C3C4C1C5C2C3C45").unwrap();
        for input in [
            cosmolkit_core::fast_find_rings(&record.topology).unwrap(),
            cosmolkit_core::find_sssr(&record.topology, &Default::default()).unwrap(),
        ] {
            assert!(!input.is_symm_sssr());
            record.rings = input;
            score_tautomer_rings_(&mut record.score_view(&coordinates)).unwrap();
            assert!(record.rings.is_symm_sssr());
            assert_eq!(record.rings.num_rings(), 6);
        }
        let rows = record.rings.bond_rings().as_ptr();
        score_tautomer_rings_(&mut record.score_view(&coordinates)).unwrap();
        assert_eq!(rows, record.rings.bond_rings().as_ptr());
    }
    #[test]
    fn recovery_search06_exclusive_loan_keeps_all_non_cache_blocks() {
        let coordinates = cosmolkit_model::CoordinateBlock::default();
        let mut record = crate::engine::stereo_tests::fixture_from_smiles("c1ccccc1").unwrap();
        record.rings = cosmolkit_core::RingInfo::new(
            cosmolkit_core::RingFindType::OtherOrUnknown,
            record.topology.atoms.len(),
            record.topology.bonds.len(),
        );
        let topology = record.topology.clone();
        let properties = record.properties.clone();
        let valence = record.valence.clone();
        assert_eq!(
            score_tautomer_rings_(&mut record.score_view(&coordinates)).unwrap(),
            250
        );
        assert!(record.rings.is_symm_sssr());
        assert_eq!(record.topology, topology);
        assert_eq!(record.properties, properties);
        assert_eq!(record.valence, valence);
    }
}
