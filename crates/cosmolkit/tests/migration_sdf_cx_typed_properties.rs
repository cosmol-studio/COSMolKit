#![cfg(all(feature = "cap-io", feature = "cap-smiles"))]

use cosmolkit::{CxSmilesFields, CxSmilesWriteParams, Molecule, PropertyValue, SmilesWriteParams};

fn one_atom_sdf(fields: &str) -> String {
    format!(
        concat!(
            "typed CX properties\n  COSMolKit         2D\n\n",
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "M  END\n{}$$$$\n",
        ),
        fields,
    )
}

fn atom_property_params() -> CxSmilesWriteParams {
    CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: false,
            kekule: false,
            ..Default::default()
        },
        fields: CxSmilesFields::ATOM_PROPS,
        ..Default::default()
    }
}

#[test]
fn sdf_typed_atom_properties_flow_to_cx_with_source_order_and_projection() {
    let input = one_atom_sdf(concat!(
        ">  <atom.prop.Text>\nfirst\n\n",
        ">  <atom.iprop.Count>\n+007\n\n",
        ">  <atom.dprop.Real>\n-0\n\n",
        ">  <atom.bprop.Active>\n1\n\n",
        ">  <atom.dprop.Scale>\n1e2\n\n",
        ">  <atom.prop.Text>\na.b\n\n",
        ">  <atom.iprop.atomLabel>\n9\n\n",
        ">  <atom.dprop._private>\n2.5\n\n",
    ));
    let molecule = Molecule::from_sdf(&input).expect("public typed SDF reader");
    let atom = &molecule.atoms()[0];

    assert_eq!(
        atom.prop("Text"),
        Some(&PropertyValue::String("a.b".into()))
    );
    assert_eq!(atom.prop("Count"), Some(&PropertyValue::Int(7)));
    assert_eq!(atom.prop("Real"), Some(&PropertyValue::Double(-0.0)));
    assert_eq!(atom.prop("Active"), Some(&PropertyValue::Bool(true)));
    assert_eq!(atom.prop("Scale"), Some(&PropertyValue::Double(100.0)));
    assert_eq!(atom.prop("atomLabel"), Some(&PropertyValue::Int(9)));
    assert_eq!(atom.prop("_private"), Some(&PropertyValue::Double(2.5)));
    assert_eq!(
        molecule.properties().sdf_data_fields(),
        &[
            ("atom.prop.Text".into(), "first".into()),
            ("atom.iprop.Count".into(), "+007".into()),
            ("atom.dprop.Real".into(), "-0".into()),
            ("atom.bprop.Active".into(), "1".into()),
            ("atom.dprop.Scale".into(), "1e2".into()),
            ("atom.prop.Text".into(), "a.b".into()),
            ("atom.iprop.atomLabel".into(), "9".into()),
            ("atom.dprop._private".into(), "2.5".into()),
        ]
    );

    let topology_before = molecule.topology().clone();
    let properties_before = molecule.properties().clone();
    let coordinates_before = molecule.coordinates_2d().map(<[_]>::to_vec);
    assert_eq!(
        molecule
            .to_cx_smiles_with_params(&atom_property_params())
            .unwrap(),
        cosmolkit::PropertyText::from(
            "C |atomProp:0.Text.a&#46;b:0.Count.7:0.Real.-0:0.Active.1:0.Scale.100|"
        )
    );
    assert_eq!(molecule.topology(), &topology_before);
    assert_eq!(molecule.properties(), &properties_before);
    assert_eq!(molecule.coordinates_2d(), coordinates_before.as_deref());
}
