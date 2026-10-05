use cosmolkit_model::{
    Atom, AtomId, AtomSpec, CoordinateBlock, PropertyValue, TopologyBlock, query_substance_groups,
};
use cosmolkit_search::{
    SearchTarget, SmartsParseParams, SubstructMatchParams, parse_smarts,
    try_get_substruct_matches_with_params,
};
use cosmolkit_types::Element;

#[test]
fn pinned_scalar_properties_match_the_source_string_projection() {
    for value in [
        PropertyValue::Int(1),
        PropertyValue::UInt(1),
        PropertyValue::Bool(true),
        PropertyValue::Double(1.0),
    ] {
        for reverse in [false, true] {
            let mut query = parse_smarts("C", &SmartsParseParams::default()).unwrap();
            let text = PropertyValue::String("1".into());
            let (query_value, target_value) = if reverse {
                (text, value.clone())
            } else {
                (value.clone(), text)
            };
            query.atoms_mut()[0].set_prop("probe", query_value).unwrap();
            let atom = Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C)
                    .with_prop("probe", target_value)
                    .unwrap(),
            );
            let topology =
                TopologyBlock::try_from_parts(vec![atom], vec![], vec![], vec![]).unwrap();
            let coordinates = CoordinateBlock::default();
            let target =
                SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None);
            let params = SubstructMatchParams {
                atom_properties: vec!["probe".into()],
                ..Default::default()
            };
            let matched = try_get_substruct_matches_with_params(&target, &query, &params).unwrap();
            assert_eq!(matched.len(), 1, "value={value:?}, reverse={reverse}");
            assert_eq!(matched[0].atom_mapping, vec![0]);
            let mut mismatch = query.clone();
            mismatch.atoms_mut()[0]
                .set_prop("probe", "different")
                .unwrap();
            assert!(
                try_get_substruct_matches_with_params(&target, &mismatch, &params)
                    .unwrap()
                    .is_empty()
            );
            let mut missing = query.clone();
            missing.atoms_mut()[0].clear_prop("probe");
            assert!(
                try_get_substruct_matches_with_params(&target, &missing, &params)
                    .unwrap()
                    .is_empty()
            );
        }
    }
}

#[test]
fn pinned_query_cleanup_preserves_original_sgroup_cx_index() {
    for text in [
        "CC |SgD:0:FIELD:DATA:QUERY:INFO:TAG:|",
        "CC |SgD:0:FIELD:DATA::::|",
    ] {
        let query = parse_smarts(text, &SmartsParseParams::default()).unwrap();
        let groups = query_substance_groups(&query);
        assert_eq!(groups.len(), 1, "{text}");
        assert_eq!(
            groups[0].props().get("_cxsmilesindex").unwrap().as_str(),
            "0",
            "{text}"
        );
        assert!(
            query
                .atoms()
                .iter()
                .all(|a| a.prop("_SmilesStart").is_none() && a.prop("_RingClosures").is_none())
        );
        assert!(
            query
                .bonds()
                .iter()
                .all(|b| b.bond().prop("_cxsmilesBondIdx").is_none())
        );
    }
}
