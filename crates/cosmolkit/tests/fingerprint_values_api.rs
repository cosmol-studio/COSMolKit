#![cfg(feature = "fingerprints")]

use cosmolkit::{BINDING_CONTRACT, BindingExposure, BindingItem, BindingParity, BindingSupport};

#[test]
fn sparse_count_public_contract_has_exact_registered_surface() {
    let methods = [
        "new",
        "length",
        "value",
        "set_value",
        "nonzero_elements",
        "total_value",
        "fuzzy_and",
        "fuzzy_or",
        "with_added",
        "with_subtracted",
        "with_added_scalar",
        "with_subtracted_scalar",
        "with_multiplied_scalar",
        "with_divided_scalar",
    ];
    for name in ["SparseCountFingerprint", "SparseCountFingerprint32"] {
        let prefix = format!("{name}.");
        let rows: Vec<_> = BINDING_CONTRACT
            .iter()
            .filter(|e| e.semantic_id.starts_with(&prefix))
            .collect();
        assert_eq!(rows.len(), methods.len());
        for method in methods {
            let row = rows
                .iter()
                .find(|e| e.semantic_id == format!("{name}.{method}"))
                .unwrap();
            assert_eq!(row.item, BindingItem::Callable);
            assert_eq!(row.feature, "fingerprints");
            assert_eq!(row.exposure, BindingExposure::Public);
            assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
            assert_eq!(row.parity, BindingParity::RequiredNow);
            assert_eq!(row.python_name, method);
            assert!(row.callable.unwrap().operation_semantic_id.is_none());
        }
    }
}

#[test]
fn runtime_sparse_count_exports_are_owner_values_with_canonical_methods() {
    // Fixed public-boundary regression, not a corpus or oracle runner.
    macro_rules! check {
        ($ty:ty) => {{
            let mut left = <$ty>::new(16);
            let mut right = <$ty>::new(16);
            for (k, v) in [(1, 5), (3, -2), (8, 4)] {
                left.set_value(k, v).unwrap();
            }
            for (k, v) in [(1, 3), (3, -4), (9, 7)] {
                right.set_value(k, v).unwrap();
            }
            let before = (left.clone(), right.clone());
            assert_eq!(left.length(), 16);
            assert_eq!(left.value(2).unwrap(), 0);
            assert_eq!(
                left.fuzzy_and(&right).unwrap().nonzero_elements(),
                &[(1, 3), (3, -4)].into()
            );
            assert_eq!(
                left.fuzzy_or(&right).unwrap().nonzero_elements(),
                &[(1, 5), (3, -2), (8, 4), (9, 7)].into()
            );
            assert_eq!(left.total_value(false).unwrap(), 7);
            assert_eq!(left.total_value(true).unwrap(), 11);
            assert_eq!(left.with_added(&right).unwrap().value(1).unwrap(), 8);
            assert_eq!(left.with_subtracted(&right).unwrap().value(1).unwrap(), 2);
            assert_eq!(left.with_added_scalar(2).unwrap().value(3).unwrap(), 0);
            assert_eq!(left.with_subtracted_scalar(2).unwrap().value(1).unwrap(), 3);
            assert_eq!(
                left.with_multiplied_scalar(2).unwrap().value(1).unwrap(),
                10
            );
            assert_eq!(left.with_divided_scalar(2).unwrap().value(1).unwrap(), 2);
            assert_eq!((left, right), before);
        }};
    }
    check!(cosmolkit::SparseCountFingerprint);
    check!(cosmolkit::SparseCountFingerprint32);
}
