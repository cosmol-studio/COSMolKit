//! Thin read-only canonical substructure queries.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn substruct_match(
        &self,
        query: &ck::QueryGraph,
    ) -> Result<Option<ck::MatchResult>, ck::SubstructMatchError> {
        // COSMolKit❗✔️: self.inner.substruct_match(query)
        self.inner.borrow().substruct_match(query)
    }
    pub fn substruct_matches(
        &self,
        query: &ck::QueryGraph,
    ) -> Result<Vec<ck::MatchResult>, ck::SubstructMatchError> {
        // COSMolKit❗✔️: self.inner.substruct_matches(query)
        self.inner.borrow().substruct_matches(query)
    }
    pub fn has_substruct_match(
        &self,
        query: &ck::QueryGraph,
    ) -> Result<bool, ck::SubstructMatchError> {
        // COSMolKit❗✔️: self.inner.has_substruct_match(query)
        self.inner.borrow().has_substruct_match(query)
    }
    pub fn substruct_matches_with_params(
        &self,
        query: &ck::QueryGraph,
        params: &ck::SubstructMatchParams,
    ) -> Result<Vec<ck::MatchResult>, ck::SubstructMatchError> {
        // COSMolKit❗✔️: self.inner.substruct_matches_with_params(query,params)
        self.inner
            .borrow()
            .substruct_matches_with_params(query, params)
    }
    pub fn substruct_matches_compiled(
        &self,
        query: &ck::CompiledQuery,
    ) -> Result<Vec<ck::MatchResult>, ck::MatchError> {
        // COSMolKit❗✔️: self.inner.substruct_matches_compiled(query)
        self.inner.borrow().substruct_matches_compiled(query)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_five_search_calls_preserve_ordered_maps_full_policies_and_read_only_state() {
        for text in ["CCO", "c1ccccc1", "C1CCC1", "F[C@H](Cl)Br"] {
            let m = Molecule::from_smiles(text).unwrap();
            let owner = ck::Molecule::from_smiles(text).unwrap();
            let before = m.inner.borrow().clone();
            for smarts in ["C", "CO", "[R]", "[$(CO)]", "[#6,#8]", "N", "F[C@H](Cl)Br"] {
                let q = ck::parse_smarts(smarts).unwrap();
                let compiled = ck::compile_query(&q).unwrap();
                assert_eq!(
                    m.substruct_match(&q).map_err(|e| e.to_string()),
                    owner.substruct_match(&q).map_err(|e| e.to_string())
                );
                assert_eq!(
                    m.substruct_matches(&q).map_err(|e| e.to_string()),
                    owner.substruct_matches(&q).map_err(|e| e.to_string())
                );
                assert_eq!(
                    m.has_substruct_match(&q).map_err(|e| e.to_string()),
                    owner.has_substruct_match(&q).map_err(|e| e.to_string())
                );
                assert_eq!(
                    m.substruct_matches_compiled(&compiled)
                        .map_err(|e| e.to_string()),
                    owner
                        .substruct_matches_compiled(&compiled)
                        .map_err(|e| e.to_string())
                );
                let mut policies = Vec::new();
                for max_matches in [0, 1, 2, usize::MAX] {
                    policies.push(ck::SubstructMatchParams {
                        max_matches,
                        ..Default::default()
                    });
                }
                for i in 0..11 {
                    let mut p = ck::SubstructMatchParams::default();
                    match i {
                        0 => p.uniquify = false,
                        1 => p.use_chirality = true,
                        2 => p.use_enhanced_stereo = true,
                        3 => p.specified_stereo_query_matches_unspecified = true,
                        4 => p.use_query_query_matches = true,
                        5 => p.recursion_possible = false,
                        6 => p.aromatic_matches_conjugated = true,
                        7 => p.aromatic_matches_single_or_double = true,
                        8 => p.extra_atom_check_overrides_default_check = true,
                        9 => p.extra_bond_check_overrides_default_check = true,
                        _ => p.use_generic_matchers = true,
                    };
                    policies.push(p);
                }
                for num_threads in [i32::MIN, -1, 0, 1, 2, i32::MAX] {
                    policies.push(ck::SubstructMatchParams {
                        num_threads,
                        ..Default::default()
                    });
                }
                policies.push(ck::SubstructMatchParams {
                    atom_properties: vec!["key".into()],
                    bond_properties: vec!["key".into()],
                    max_recursive_matches: 1,
                    ..Default::default()
                });
                for p in policies {
                    assert_eq!(
                        m.substruct_matches_with_params(&q, &p)
                            .map_err(|e| e.to_string()),
                        owner
                            .substruct_matches_with_params(&q, &p)
                            .map_err(|e| e.to_string())
                    );
                }
            }
            assert_eq!(*m.inner.borrow(), before);
        }
    }
}
