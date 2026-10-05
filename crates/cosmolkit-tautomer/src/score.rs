//! Source tautomer score terms and signed component accumulation.
use crate::engine::{TautomerRecordView, TautomerRunError, query_matches, total_hydrogens};
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
        // RDKit✔️✔️: SubstructTerm::SubstructTerm(std::string aname, std::string asmarts, int ascore)
        // RDKit✔️✔️:     : name(std::move(aname)), smarts(std::move(asmarts)), score(ascore) {
        // RDKit✔️✔️:   std::unique_ptr<ROMol> pattern(SmartsToMol(smarts));
        // RDKit✔️✔️:   if (pattern) {
        // RDKit✔️✔️:     matcher = std::move(*pattern);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        let name = name.into();
        let smarts = smarts.into();
        let matcher = parse_smarts(&smarts, &SmartsParseParams::default()).ok();
        Self {
            name,
            smarts,
            score,
            matcher,
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
    // RDKit✔️✔️: const std::vector<SubstructTerm> &getDefaultTautomerScoreSubstructs() {
    // RDKit✔️✔️:   static std::vector<SubstructTerm> substructureTerms{
    static TERMS: OnceLock<Vec<TautomerScoreTerm>> = OnceLock::new();
    TERMS.get_or_init(|| {
        vec![
            // RDKit✔️✔️:       {"benzoquinone", "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]",
            // RDKit✔️✔️:        25},
            TautomerScoreTerm::new(
                "benzoquinone",
                "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]",
                25,
            ),
            // RDKit✔️✔️:       {"oxim", "[#6]=[N][OH]", 4},
            TautomerScoreTerm::new("oxim", "[#6]=[N][OH]", 4),
            // RDKit✔️✔️:       {"C=O", "[#6]=,:[#8]", 2},
            TautomerScoreTerm::new("C=O", "[#6]=,:[#8]", 2),
            // RDKit✔️✔️:       {"N=O", "[#7]=,:[#8]", 2},
            TautomerScoreTerm::new("N=O", "[#7]=,:[#8]", 2),
            // RDKit✔️✔️:       {"P=O", "[#15]=,:[#8]", 2},
            TautomerScoreTerm::new("P=O", "[#15]=,:[#8]", 2),
            // RDKit✔️✔️:       {"C=hetero", "[C]=[!#1;!#6]", 1},
            TautomerScoreTerm::new("C=hetero", "[C]=[!#1;!#6]", 1),
            // RDKit✔️✔️:       {"C(=hetero)-hetero", "[C](=[!#1;!#6])[!#1;!#6]", 2},
            TautomerScoreTerm::new("C(=hetero)-hetero", "[C](=[!#1;!#6])[!#1;!#6]", 2),
            // RDKit✔️✔️:       {"aromatic C = exocyclic N", "[c]=!@[N]", -1},
            TautomerScoreTerm::new("aromatic C = exocyclic N", "[c]=!@[N]", -1),
            // RDKit✔️✔️:       {"methyl", "[CX4H3]", 1},
            TautomerScoreTerm::new("methyl", "[CX4H3]", 1),
            // RDKit✔️✔️:       {"guanidine terminal=N", "[#7]C(=[NR0])[#7H0]", 1},
            TautomerScoreTerm::new("guanidine terminal=N", "[#7]C(=[NR0])[#7H0]", 1),
            // RDKit✔️✔️:       {"guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2},
            TautomerScoreTerm::new("guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2),
            // RDKit✔️✔️:       {"aci-nitro", "[#6]=[N+]([O-])[OH]", -4}};
            TautomerScoreTerm::new("aci-nitro", "[#6]=[N+]([O-])[OH]", -4),
        ]
    })
    // RDKit✔️✔️:   return substructureTerms;
    // RDKit✔️✔️: }
}

/// Score aromatic rings using RDKit's tautomer-canonicalization weights.
pub fn score_tautomer_rings(molecule: TautomerRecordView<'_>) -> Result<i32, RingFindingError> {
    // RDKit✔️✔️: int scoreRings(const ROMol &mol) {
    // RDKit✔️✔️:   int score = 0;
    let mut score = 0_i32;

    // RDKit✔️✔️:   auto ringInfo = mol.getRingInfo();
    // RDKit✔️✔️:   std::unique_ptr<ROMol> cp;
    // RDKit✔️✔️:   if (!ringInfo->isSymmSssr()) {
    // RDKit✔️✔️:     cp.reset(new ROMol(mol));
    // RDKit✔️✔️:     MolOps::symmetrizeSSSR(*cp);
    // RDKit✔️✔️:     ringInfo = cp->getRingInfo();
    // RDKit✔️✔️:   }
    let computed_rings;
    let ring_info = match molecule.rings {
        Some(rings) if rings.is_symm_sssr() => rings,
        _ => {
            computed_rings = symmetrized_sssr(molecule.topology, &Default::default())?;
            &computed_rings
        }
    };

    // RDKit✔️✔️:   boost::dynamic_bitset<> isArom(mol.getNumBonds());
    // RDKit✔️✔️:   boost::dynamic_bitset<> bothCarbon(mol.getNumBonds());
    let mut is_aromatic = vec![false; molecule.topology.bonds.len()];
    let mut both_carbon = vec![false; molecule.topology.bonds.len()];
    // RDKit✔️✔️:   for (const auto &bnd : mol.bonds()) {
    for bond in &molecule.topology.bonds {
        // RDKit✔️✔️:     if (bnd->getIsAromatic()) {
        if bond.is_aromatic() {
            // RDKit✔️✔️:       isArom.set(bnd->getIdx());
            is_aromatic[bond.id().index()] = true;
            // RDKit✔️✔️:       if (bnd->getBeginAtom()->getAtomicNum() == 6 &&
            // RDKit✔️✔️:           bnd->getEndAtom()->getAtomicNum() == 6) {
            if molecule.topology.atoms[bond.begin().index()].atomic_number() == 6
                && molecule.topology.atoms[bond.end().index()].atomic_number() == 6
            {
                // RDKit✔️✔️:         bothCarbon.set(bnd->getIdx());
                both_carbon[bond.id().index()] = true;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   for (const auto &bring : ringInfo->bondRings()) {
    for bond_ring in ring_info.bond_rings() {
        // RDKit✔️✔️:     bool allC = true;
        // RDKit✔️✔️:     bool allAromatic = true;
        let mut all_carbon = true;
        let mut all_aromatic = true;
        // RDKit✔️✔️:     for (const auto bidx : bring) {
        for bond_id in bond_ring {
            // RDKit✔️✔️:       if (!isArom[bidx]) {
            if !is_aromatic[bond_id.index()] {
                // RDKit✔️✔️:         allAromatic = false;
                all_aromatic = false;
                // RDKit✔️✔️:         break;
                break;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       if (!bothCarbon[bidx]) {
            if !both_carbon[bond_id.index()] {
                // RDKit✔️✔️:         allC = false;
                all_carbon = false;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     if (allAromatic) {
        if all_aromatic {
            // RDKit✔️✔️:       score += 100;
            score += 100;
            // RDKit✔️✔️:       if (allC) {
            if all_carbon {
                // RDKit✔️✔️:         score += 150;
                score += 150;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return score;
    Ok(score)
    // RDKit✔️✔️: };
}

/// Score every match of the supplied source-ordered tautomer SMARTS terms.
#[must_use]
pub fn score_tautomer_substructures(
    molecule: TautomerRecordView<'_>,
    terms: &[TautomerScoreTerm],
) -> Result<i32, TautomerRunError> {
    // RDKit✔️✔️: int scoreSubstructs(const ROMol &mol,
    // RDKit✔️✔️:                     const std::vector<SubstructTerm> &substructureTerms) {
    // RDKit✔️✔️:   int score = 0;
    let mut score = 0_i32;
    // RDKit✔️✔️:   for (const auto &term : substructureTerms) {
    for term in terms {
        // RDKit✔️✔️:     if (!term.matcher.getNumAtoms()) {
        let Some(matcher) = term.matcher() else {
            // RDKit✔️✔️:       BOOST_LOG(rdErrorLog) << " matcher for term " << term.name
            // RDKit✔️✔️:                             << " is invalid, ignoring it." << std::endl;
            // RDKit✔️✔️:       continue;
            continue;
            // RDKit✔️✔️:     }
        };
        // RDKit✔️✔️:     SubstructMatchParameters params;
        // RDKit✔️✔️:     const auto matches = SubstructMatch(mol, term.matcher, params);
        // Valid SMARTS compilation/matching errors remain typed failures.
        // Source invalid/empty SMARTS is the sole skip branch above.
        if matcher.num_atoms() == 0 {
            continue;
        }
        let plan = CompiledQuery::compile(matcher.clone())?;
        let matches = query_matches(molecule, &plan)?;
        // RDKit✔️✔️:     // if (!matches.empty()) {
        // RDKit✔️✔️:     //   std::cerr << " " << matches.size() << " matches to " << term.name
        // RDKit✔️✔️:     //             << std::endl;
        // RDKit✔️✔️:     // }
        // RDKit✔️✔️:     score += static_cast<int>(matches.size()) * term.score;
        // The pinned 32-bit source ABI narrows `size_t` to `int` here; Rust's
        // cast and wrapping arithmetic retain that release-build machine
        // behavior instead of introducing a debug-only overflow branch.
        let match_count = matches.len() as i32;
        score = score.wrapping_add(match_count.wrapping_mul(term.score));
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return score;
    Ok(score)
    // RDKit✔️✔️: }
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
pub fn score_tautomer_with_terms(
    molecule: TautomerRecordView<'_>,
    terms: &[TautomerScoreTerm],
) -> Result<TautomerScore, TautomerRunError> {
    // RDKit✔️✔️: inline int scoreTautomer(const ROMol &mol) {
    // RDKit✔️✔️:   return scoreRings(mol) + scoreSubstructs(mol) + scoreHeteroHs(mol);
    // RDKit✔️✔️: }
    Ok(TautomerScore {
        ring: score_tautomer_rings(molecule)?,
        substructure: score_tautomer_substructures(molecule, terms)?,
        hetero_hydrogen: score_tautomer_hetero_hydrogens(molecule)?,
    })
}

/// Compute the source tautomer score using the pinned default SMARTS terms.
pub fn score_tautomer(molecule: TautomerRecordView<'_>) -> Result<TautomerScore, TautomerRunError> {
    score_tautomer_with_terms(molecule, default_tautomer_score_terms())
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let molecule = fixture("CCO").expect("parse acyclic molecule");

        assert_eq!(
            score_tautomer_rings(molecule.view()).expect("score rings"),
            0
        );
    }

    #[test]
    fn ring_score_counts_an_aromatic_heteroring_without_the_all_carbon_bonus() {
        let molecule = fixture("c1ccncc1").expect("parse pyridine");

        assert_eq!(
            score_tautomer_rings(molecule.view()).expect("score rings"),
            100
        );
    }

    #[test]
    fn ring_score_adds_both_aromatic_and_all_carbon_weights() {
        let molecule = fixture("c1ccccc1").expect("parse benzene");

        assert_eq!(
            score_tautomer_rings(molecule.view()).expect("score rings"),
            250
        );
    }

    #[test]
    fn ring_score_traverses_every_symmetrized_fused_ring() {
        let molecule = fixture("c1ccc2ccccc2c1").expect("parse naphthalene");

        assert_eq!(
            symmetrized_sssr(&molecule.record.topology, &Default::default())
                .expect("find fused rings")
                .bond_rings()
                .len(),
            2
        );
        assert_eq!(
            score_tautomer_rings(molecule.view()).expect("score rings"),
            500
        );
    }

    #[test]
    fn ring_score_ignores_nonaromatic_rings() {
        let molecule = fixture("C1CCCCC1").expect("parse cyclohexane");

        assert_eq!(
            score_tautomer_rings(molecule.view()).expect("score rings"),
            0
        );
    }

    #[test]
    fn ring_score_matches_with_or_without_a_precomputed_symm_sssr_cache() {
        let without_cache =
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
            score_tautomer_rings(without_cache.view()).expect("score uncached"),
            100
        );
        assert_eq!(
            score_tautomer_rings(with_cache.view()).expect("score cached"),
            100
        );
    }

    #[test]
    fn ring_score_does_not_modify_the_source_or_materialize_its_cache() {
        let molecule = fixture_with_sanitize("c1ccccc1", false).expect("parse uncached benzene");
        let topology_before = molecule.record.topology.clone();
        assert!(molecule.rings.is_none());

        assert_eq!(
            score_tautomer_rings(molecule.view()).expect("score rings"),
            250
        );
        assert_eq!(&molecule.record.topology, &topology_before);
        assert!(molecule.rings.is_none());
    }

    #[test]
    fn substructure_score_is_zero_when_no_term_matches() {
        let molecule = fixture("CCC").expect("parse propane");
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
        let molecule =
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
        let molecule = fixture("CCO").expect("parse ethanol");
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
        let molecule = fixture("CCO").expect("parse ethanol");
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
            let molecule = fixture(smiles).expect("parse heteroatom");
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
            let molecule = fixture(smiles).expect("parse charged heteroatom");
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
            let molecule = fixture(smiles).expect("parse excluded element");
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
        let benzene = fixture("c1ccccc1").expect("parse benzene");
        let carbon_terms = [TautomerScoreTerm::new("carbon", "[#6]", 2)];
        let score =
            score_tautomer_with_terms(benzene.view(), &carbon_terms).expect("score benzene");

        assert_eq!(score.ring(), 250);
        assert_eq!(score.substructure(), 12);
        assert_eq!(score.hetero_hydrogen(), 0);
        assert_eq!(score.total(), 262);
    }

    #[test]
    fn aggregate_score_default_delegate_uses_the_single_builtin_term_table() {
        let molecule = fixture("CC(=O)N").expect("parse acetamide");
        let delegated = score_tautomer(molecule.view()).expect("score with defaults");
        let explicit = score_tautomer_with_terms(molecule.view(), default_tautomer_score_terms())
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
