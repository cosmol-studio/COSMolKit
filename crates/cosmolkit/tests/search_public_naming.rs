#![cfg(feature = "cap-search")]

use cosmolkit::{
    BINDING_CONTRACT, BindingDefault, BindingKind, BindingOwner, QueryGraph, SmartsParseError,
    SmartsParseParams,
};

#[test]
fn flat_functions_and_bound_query_factories_share_the_parser() {
    let _: fn(&str) -> Result<QueryGraph, SmartsParseError> = cosmolkit::parse_smarts;
    let _: fn(&str) -> Result<QueryGraph, SmartsParseError> = cosmolkit::search::from_smarts;
    for (text, atoms, bonds) in [("C(N)O", 3, 2), ("C1CC1", 3, 3), ("N<-C", 2, 1)] {
        let flat = cosmolkit::parse_smarts(text).unwrap();
        let factory = cosmolkit::search::from_smarts(text).unwrap();
        for query in [&flat, &factory] {
            assert_eq!((query.num_atoms(), query.num_bonds()), (atoms, bonds));
        }
        assert_eq!(
            cosmolkit::write_smarts(&flat, &Default::default()).unwrap(),
            cosmolkit::write_smarts(&factory, &Default::default()).unwrap()
        );
    }
    let params = SmartsParseParams {
        merge_hs: true,
        ..Default::default()
    };
    for query in [
        cosmolkit::parse_smarts_with_params("C[H]", &params).unwrap(),
        cosmolkit::search::from_smarts_with_params("C[H]", &params).unwrap(),
    ] {
        assert_eq!((query.num_atoms(), query.num_bonds()), (1, 0));
    }
    assert!(params.merge_hs);
    for text in ["C(", "[#6", "C |notCX| name"] {
        let flat = cosmolkit::parse_smarts(text).unwrap_err();
        let factory = cosmolkit::search::from_smarts(text).unwrap_err();
        assert_eq!(format!("{flat:?}"), format!("{factory:?}"));
    }
}

#[test]
fn registry_records_both_static_factories_and_the_flat_functions() {
    for (id, name, parameter_names) in [
        ("QueryGraph.from_smarts", "from_smarts", &["text"][..]),
        (
            "QueryGraph.from_smarts_with_params",
            "from_smarts_with_params",
            &["text", "params"][..],
        ),
    ] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == id)
            .unwrap();
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.python_name, name);
        assert_eq!(row.feature, "cap-search");
        let callable = row.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Static);
        assert_eq!(callable.receiver, None);
        assert_eq!(
            callable
                .parameters
                .iter()
                .map(|p| p.name)
                .collect::<Vec<_>>(),
            parameter_names
        );
        assert!(
            callable
                .parameters
                .iter()
                .all(|p| p.default == BindingDefault::Required)
        );
    }
    for name in [
        "parse_smarts",
        "parse_smarts_with_params",
        "compile_query",
        "write_smarts",
        "write_cx_smarts",
    ] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == format!("search.{name}"))
            .unwrap();
        assert_eq!(row.owner, BindingOwner::Module);
        assert_eq!(row.python_name, name);
        let path: String = row
            .rust_path
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect();
        assert_eq!(path, format!("crate::{name}"));
    }
}
