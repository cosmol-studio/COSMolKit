#[path = "../src/binding.rs"]
mod binding;
#[path = "../src/status.rs"]
mod status;

fn declaration(extra: &str) -> String {
    format!(
        r#"pub static API = [
        #[cfg(feature = "cap-smiles")]
        {{semantic_id:"Molecule.from_smiles",item:callable,owner:molecule,
          rust:crate::Molecule::from_smiles,python:"from_smiles",javascript:"fromSmiles",
          feature:"cap-smiles",requires:["cap-batch"],kind:static_,
          parameters:[{{name:input,type:&str,default:required}}],
          output:crate::Molecule,error:crate::Error,state:value_returning,operation:none,
          signature:for<'a> fn(&'a str)->Result<crate::Molecule,crate::Error>,{extra}}}
    ];"#
    )
}

#[test]
fn convenience_inherits_signature_lifetimes_and_all_feature_gates() {
    let output = binding::expand_binding_contract(
        declaration("convenience:mol_from_smiles,").parse().unwrap(),
    )
    .unwrap();
    let tree = syn::parse2::<syn::File>(output.clone()).unwrap();
    let module = tree
        .items
        .iter()
        .find_map(|item| match item {
            syn::Item::Mod(module) if module.ident == "conveniences" => Some(module),
            _ => None,
        })
        .unwrap();
    let syn::Item::Fn(function) = &module.content.as_ref().unwrap().1[0] else {
        panic!("factory function missing")
    };
    assert_eq!(function.sig.ident, "mol_from_smiles");
    assert_eq!(function.sig.generics.lifetimes().count(), 1);
    assert_eq!(function.sig.inputs.len(), 1);
    let compact = output
        .to_string()
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect::<String>();
    assert!(compact.contains("input:&'astr"));
    assert!(compact.contains("crate::Molecule::from_smiles(input)"));
    assert!(compact.contains("javascript_target:Some(\"Molecule.fromSmiles\")"));
    assert!(compact.contains("API_CONVENIENCES"));
    let attributes = &function.attrs;
    assert!(
        quote::quote!(#(#attributes)*)
            .to_string()
            .contains("cap-batch")
    );
    assert!(!compact.contains("clone("));
}

#[test]
fn invalid_duplicate_and_instance_conveniences_are_rejected() {
    for name in ["_hidden", "mutating_", "UpperCase"] {
        let error = binding::expand_binding_contract(
            declaration(&format!("convenience:{name},"))
                .parse()
                .unwrap(),
        )
        .unwrap_err();
        assert!(error.to_string().contains("public snake_case"), "{error}");
    }
    let source = declaration("convenience:mol_from_smiles,");
    let instance = source.replace("kind:static_", "kind:instance");
    assert!(
        binding::expand_binding_contract(instance.parse().unwrap())
            .unwrap_err()
            .to_string()
            .contains("static value-returning factory")
    );
    let duplicate = source.replace(
        "convenience:mol_from_smiles,",
        "convenience:mol_from_smiles,convenience:other,",
    );
    assert!(binding::expand_binding_contract(duplicate.parse().unwrap()).is_err());
}

#[test]
fn python_object_adapter_emits_no_rust_or_javascript_factory() {
    let source = r#"pub static API = [
      {semantic_id:"types.Value",item:type,owner:type_,rust:crate::Value,
       python:"Value",javascript:"Value",feature:"runtime",role:value,
       python_adapters:[{name:from_object,convenience:value_from_object,targets:["Value.new"]}]},
      {semantic_id:"Value.new",item:callable,owner:type_,rust:crate::Value::new,
       python:"new",javascript:"new",feature:"runtime",kind:static_,parameters:[],
       output:crate::Value,error:none,state:value_returning,operation:none,signature:fn()->crate::Value}
    ];"#;
    let output = binding::expand_binding_contract(source.parse().unwrap())
        .unwrap()
        .to_string();
    assert!(output.contains("value_from_object"));
    assert!(output.contains("Value.from_object"));
    assert!(output.contains("javascript_target : None"));
    assert!(!output.contains("fn value_from_object"));
}
