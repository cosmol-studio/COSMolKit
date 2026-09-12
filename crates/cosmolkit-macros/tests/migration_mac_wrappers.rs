#[path = "../src/declaration.rs"]
mod declaration;
#[path = "../src/wrappers.rs"]
mod wrappers;

use declaration::{BioRegistry, MoleculeRegistry};
use proc_macro2::TokenStream;
use wrappers::{expand_bio_wrappers, expand_molecule_wrappers};

fn compact(tokens: &TokenStream) -> String {
    tokens
        .to_string()
        .chars()
        .filter(|character| !character.is_whitespace())
        .collect()
}

fn molecule(source: &str) -> String {
    let registry = syn::parse_str::<MoleculeRegistry>(source)
        .unwrap_or_else(|error| panic!("expected Molecule registry: {error}"));
    compact(
        &expand_molecule_wrappers(&registry)
            .unwrap_or_else(|error| panic!("expected Molecule wrappers: {error}")),
    )
}

fn molecule_error(source: &str) -> String {
    syn::parse_str::<MoleculeRegistry>(source)
        .err()
        .expect("invalid Molecule declaration unexpectedly parsed")
        .to_string()
}

fn bio(source: &str) -> String {
    let registry = syn::parse_str::<BioRegistry>(source)
        .unwrap_or_else(|error| panic!("expected Bio registry: {error}"));
    compact(
        &expand_bio_wrappers(&registry)
            .unwrap_or_else(|error| panic!("expected Bio wrappers: {error}")),
    )
}

fn molecule_operation(name: &str, parameters: &str, extra: &str) -> String {
    format!(
        r#"
        op {name}({parameters}) {{
            method: {name},
            impl_fn: crate::operations::{name}_impl,
            kind: weak,
            access: {{ read: [], write: [] }},
            derived_effects: {{
                recompute: [], preserve: [], invalidate: [], operation_defined: [],
            }},
            cip_state: preserve,
            feature: crate::capabilities::{name}_FEATURE,
            parity: not_applicable,
            invariant_profile: "wrapper",
            {extra}
        }}
        "#
    )
}

#[test]
fn empty_registries_emit_only_empty_runtime_owned_inherent_impls() {
    assert_eq!(molecule(""), "implcrate::Molecule{}");
    assert_eq!(bio(""), "implcrate::BioStructure{}");
}

#[test]
fn single_untyped_value_wrapper_checks_support_forwards_once_and_finishes() {
    let output = molecule(&molecule_operation(
        "inspect",
        "first: usize, second: Option<&'static str>",
        r#"docs: "inspect docs", default_method: inspect_default, default_args: [7, None],"#,
    ));

    for expected in [
        "implcrate::Molecule",
        "#[doc=\"inspectdocs\"]",
        "pubfninspect(&self,first:usize,second:Option<&'staticstr>)->Result<crate::Molecule,crate::ops::OperationError>",
        "crate::UnsupportedFeatureError::from_spec(&crate::capabilities::inspect_FEATURE)",
        "letmutparts=crate::OpParts::new(self,&INSPECT_SPEC)?",
        "crate::operations::inspect_impl(&mutparts,first,second)?",
        "parts.finish()",
        "pubfninspect_default(&self)->Result<crate::Molecule,crate::ops::OperationError>{self.inspect(7,None)}",
    ] {
        assert!(output.contains(expected), "missing {expected} in {output}");
    }
    assert_eq!(output.matches("crate::operations::inspect_impl").count(), 1);
    assert_eq!(output.matches("crate::OpParts::new(").count(), 1);
    assert!(output.find("UnsupportedFeature").unwrap() < output.find("OpParts::new").unwrap());
    assert!(output.find("inspect_impl").unwrap() < output.find("parts.finish").unwrap());
}

#[test]
fn single_typed_value_and_in_place_wrappers_finish_before_exposing_metadata() {
    let output = molecule(&molecule_operation(
        "measure",
        "scale: usize",
        r#"
            result_type: crate::Report,
            inplace: true,
            inplace_method: measure_,
        "#,
    ));

    assert!(output.contains(
        "pubfnmeasure(&self,scale:usize)->Result<(crate::Molecule,crate::Report),crate::ops::OperationError>"
    ));
    assert!(output.contains("letresult=crate::operations::measure_impl(&mutparts,scale)?"));
    assert!(output.contains("letmolecule=parts.finish()?;Ok((molecule,result))"));
    assert!(output.contains(
        "pubfnmeasure_(&mutself,scale:usize)->Result<crate::Report,crate::ops::OperationError>"
    ));
    assert!(output.contains("letmutparts=crate::OpParts::new_in_place(self,&MEASURE_SPEC)?"));
    assert!(output.contains("Err(error)=>{parts.abort_in_place();returnErr(error);}"));
    assert!(output.contains("parts.finish_in_place()?;Ok(result)"));
    assert_eq!(output.matches("crate::operations::measure_impl").count(), 2);
}

#[test]
fn untyped_in_place_and_default_aliases_preserve_abort_finish_and_arguments() {
    let output = molecule(&molecule_operation(
        "normalize",
        "mode: usize, enabled: bool",
        r#"
            default_method: normalize_default,
            default_args: [3, true],
            inplace: true,
            inplace_method: normalize_with_params_,
            default_inplace_method: normalize_,
        "#,
    ));

    assert!(output.contains(
        "ifletErr(error)=crate::operations::normalize_impl(&mutparts,mode,enabled){parts.abort_in_place();returnErr(error);}parts.finish_in_place()"
    ));
    assert!(output.contains("self.normalize(3,true)"));
    assert!(output.contains("self.normalize_with_params_(3,true)"));
    assert_eq!(output.matches("OpParts::new_in_place").count(), 1);
    assert_eq!(output.matches("abort_in_place").count(), 1);
    assert_eq!(output.matches("finish_in_place").count(), 1);
}

#[test]
fn trailing_defaults_keep_leading_operation_parameters_on_both_short_methods() {
    let output = molecule(&molecule_operation(
        "position_with_params",
        "atom: usize, point: [f64; 3], params: &crate::Params",
        r#"
            default_method: position,
            default_args: [&crate::Params::default()],
            inplace: true,
            inplace_method: set_position_with_params_,
            default_inplace_method: set_position_,
        "#,
    ));

    assert!(output.contains(
        "pubfnposition(&self,atom:usize,point:[f64;3])->Result<crate::Molecule,crate::ops::OperationError>"
    ));
    assert!(output.contains("self.position_with_params(atom,point,&crate::Params::default())"));
    assert!(output.contains(
        "pubfnset_position_(&mutself,atom:usize,point:[f64;3])->Result<(),crate::ops::OperationError>"
    ));
    assert!(
        output.contains("self.set_position_with_params_(atom,point,&crate::Params::default())")
    );
}

#[test]
fn plain_multiple_output_finishes_the_unique_transaction_before_return() {
    let output = molecule(&molecule_operation(
        "enumerate",
        "limit: usize",
        "output: multiple,",
    ));

    assert!(output.contains(
        "pubfnenumerate(&self,limit:usize)->Result<Vec<crate::Molecule>,crate::ops::OperationError>"
    ));
    assert!(output.contains("letmutparts=crate::MultiOutputOpParts::new(self,&ENUMERATE_SPEC)?"));
    assert!(output.contains("crate::operations::enumerate_impl(&mutparts,limit)?;parts.finish()"));
    assert_eq!(output.matches("enumerate_impl").count(), 1);
    assert!(!output.contains("&mutself"));
    assert!(!output.contains("assemble"));
}

#[test]
fn typed_multiple_output_assembles_only_after_validated_finish() {
    let output = molecule(&molecule_operation(
        "enumerate_typed",
        "limit: usize",
        r#"
            output: multiple,
            result_type: crate::EnumerationResult,
            assemble_fn: crate::results::assemble_enumeration,
        "#,
    ));

    for expected in [
        "pubfnenumerate_typed(&self,limit:usize)->Result<crate::EnumerationResult,crate::ops::OperationError>",
        "letmetadata=crate::operations::enumerate_typed_impl(&mutparts,limit)?",
        "letmolecules=parts.finish()?",
        "crate::results::assemble_enumeration(molecules,metadata)",
    ] {
        assert!(output.contains(expected), "missing {expected} in {output}");
    }
    assert!(output.find("parts.finish").unwrap() < output.find("assemble_enumeration").unwrap());
    assert_eq!(output.matches("assemble_enumeration").count(), 1);
}

#[test]
fn bio_value_wrapper_checks_support_before_one_typed_transaction() {
    let output = bio(r#"
        #[cfg(all(feature = "bio-a", not(feature = "bio-b")))]
        op remove_waters(selector: &crate::Selector, keep: bool) {
            method: remove_waters,
            impl_fn: crate::bio_operations::remove_waters_impl,
            domain: selection,
            kind: weak,
            edit_kind: local,
            may_mutate: [],
            auto_remap: [],
            must_handle: [],
            needs_update: [],
            feature: crate::bio::features::REMOVE_WATERS,
            parity: not_applicable,
            invariant_profile: "bio-wrapper",
        }
        "#);

    for expected in [
        "implcrate::BioStructure",
        "pubfnremove_waters(&self,selector:&crate::Selector,keep:bool)->Result<crate::BioStructure,crate::bio_ops::BioOperationError>",
        "BIO_REMOVE_WATERS_SPEC.support",
        "let_feature=&crate::bio::features::REMOVE_WATERS",
        "letmutparts=crate::BioOpParts::new(self,&BIO_REMOVE_WATERS_SPEC)",
        "crate::bio_operations::remove_waters_impl(&mutparts,selector,keep)?",
        "parts.finish()",
    ] {
        assert!(output.contains(expected), "missing {expected} in {output}");
    }
    assert!(output.find("SPEC.support").unwrap() < output.find("BioOpParts::new").unwrap());
    assert_eq!(output.matches("remove_waters_impl").count(), 1);
}

#[test]
fn cfg_attributes_and_complete_feature_paths_stay_operation_local() {
    let source = format!(
        "{}{}",
        molecule_operation(
            "first",
            "",
            "docs: \"first\",",
        )
        .replace(
            "        op first",
            "        #[cfg(all(feature = \"first-a\", not(feature = \"first-b\")))]\n        op first",
        ),
        molecule_operation("second", "", "")
            .replace("crate::capabilities::second_FEATURE", "external::SECOND_FEATURE")
    );
    let output = molecule(&source);

    assert_eq!(
        output
            .matches("#[cfg(all(feature=\"first-a\",not(feature=\"first-b\")))]")
            .count(),
        1
    );
    assert!(output.contains("from_spec(&external::SECOND_FEATURE)"));
    assert!(!output.contains("crate::external::SECOND_FEATURE"));
    let first_cfg = output.find("#[cfg(all").unwrap();
    let first_method = output.find("pubfnfirst").unwrap();
    let second_method = output.find("pubfnsecond").unwrap();
    assert!(first_cfg < first_method && first_method < second_method);
    assert!(!output[second_method..].contains("first-a"));
}

#[test]
fn invalid_in_place_multiple_and_assembler_combinations_fail_closed() {
    let multiple_in_place = molecule_operation("many", "", "output: multiple, inplace: true,");
    assert!(
        molecule_error(&multiple_in_place)
            .contains("multiple-output molecule operations cannot generate an in-place wrapper")
    );

    let result_without_assembler = molecule_operation(
        "many_result",
        "",
        "output: multiple, result_type: crate::Report,",
    );
    assert!(molecule_error(&result_without_assembler).contains("requires assemble_fn"));

    let assembler_without_result = molecule_operation(
        "many_assembler",
        "",
        "output: multiple, assemble_fn: crate::assemble,",
    );
    assert!(molecule_error(&assembler_without_result).contains("requires result_type"));

    let single_assembler = molecule_operation(
        "single_assembler",
        "",
        "assemble_fn: crate::assemble, result_type: crate::Report,",
    );
    assert!(molecule_error(&single_assembler).contains("only for multiple-output"));
}

#[test]
fn wrapper_tokens_do_not_duplicate_registry_runtime_or_domain_implementations() {
    let output = molecule(&molecule_operation("inspect", "", ""));
    for forbidden in [
        "pubconstINSPECT_SPEC",
        "MOLECULE_OPS",
        "SUPPORT_MATRIX",
        "OPERATION_INVARIANT_MATRIX",
        "PARITY_MATRIX",
        "structInspectAccess",
        "implcrate::OpParts",
        "fninspect_impl",
        "core_old",
    ] {
        assert!(
            !output.contains(forbidden),
            "unexpected {forbidden} in {output}"
        );
    }
}
