//! Coordinate configuration transport; installation/generation stays in Rust.
use crate::optional_configuration::Coordinate2DParams;
use crate::selectors::{Coordinate2DInputParams, CoordinateZPolicyInput};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

pub(crate) fn input_2d(
    z_policy: Option<CoordinateZPolicyInput>,
) -> PyResult<Coordinate2DInputParams> {
    Ok(Coordinate2DInputParams {
        inner: crate::parameters::NativeField::Owned(ck::Coordinate2DInputParams {
            z_policy: match z_policy {
                Some(value) => value.into_rust()?,
                None => ck::Coordinate2DInputParams::default().z_policy,
            },
        }),
    })
}

pub(crate) fn layout(
    coordinate_map: Option<BTreeMap<usize, [f64; 2]>>,
    canonical_orientation: Option<bool>,
    clear_existing_2d: Option<bool>,
    flips_per_sample: Option<u32>,
    samples: Option<u32>,
    sample_seed: Option<i32>,
    permute_degree_four: Option<bool>,
    force_rdkit: Option<bool>,
    use_ring_templates: Option<bool>,
    use_path_angle_expansion: Option<bool>,
) -> ck::Coordinate2DParams {
    let defaults = ck::Coordinate2DParams::default();
    ck::Coordinate2DParams {
        coordinate_map: coordinate_map.unwrap_or(defaults.coordinate_map),
        canonical_orientation: canonical_orientation.unwrap_or(defaults.canonical_orientation),
        clear_existing_2d: clear_existing_2d.unwrap_or(defaults.clear_existing_2d),
        flips_per_sample: flips_per_sample.unwrap_or(defaults.flips_per_sample),
        samples: samples.unwrap_or(defaults.samples),
        sample_seed: sample_seed.unwrap_or(defaults.sample_seed),
        permute_degree_four: permute_degree_four.unwrap_or(defaults.permute_degree_four),
        force_rdkit: force_rdkit.unwrap_or(defaults.force_rdkit),
        use_ring_templates: use_ring_templates.unwrap_or(defaults.use_ring_templates),
        use_path_angle_expansion: use_path_angle_expansion
            .unwrap_or(defaults.use_path_angle_expansion),
    }
}

pub(crate) fn with_layout(
    molecule: &crate::native::Molecule,
    py: Python<'_>,
    params: Option<&ck::Coordinate2DParams>,
) -> PyResult<crate::native::Molecule> {
    py.detach(|| match params {
        Some(params) => molecule.inner.with_2d_coordinates_with_params(params),
        None => molecule.inner.with_2d_coordinates(),
    })
    .map(|inner| crate::native::Molecule { inner })
    .map_err(crate::operation_errors::OperationError::from_rust)
}
pub(crate) fn add_layout(
    molecule: &mut crate::native::Molecule,
    py: Python<'_>,
    params: Option<&ck::Coordinate2DParams>,
) -> PyResult<()> {
    py.detach(|| match params {
        Some(params) => molecule.inner.add_2d_coordinates_with_params_(params),
        None => molecule.inner.add_2d_coordinates_(),
    })
    .map_err(crate::operation_errors::OperationError::from_rust)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    use pyo3::ffi::c_str;
    use pyo3::types::{PyDict, PyModule};

    #[test]
    fn coordinate_families_forward_native_bits_enums_results_and_value_inplace_ownership() {
        Python::initialize();
        Python::attach(|py|->PyResult<()> {
            let molecule=ck::Molecule::from_smiles("CCO").unwrap().with_added_3d_conformer(vec![[0.,0.,0.],[1.,0.,0.],[2.,1.,0.]]).unwrap();
            let before=molecule.to_binary().unwrap();
            let rows2=vec![vec![0.,-0.],vec![1.2,0.],vec![2.1,0.7]];
            let rows3=vec![[0.,-0.,0.1],[1.2,0.,0.2],[2.1,0.7,0.3]];
            let options2=ck::Coordinate2DInputParams {z_policy:ck::CoordinateZPolicy::RequireZero};
            let options3=ck::Coordinate3DInputParams {is_3d:false};
            let replace=ck::Replace3DCoordinatesParams {conformer_id:0};
            let expected=PyDict::new(py);
            for (name,default,configured) in [
                ("with_2d_coordinate_block",molecule.with_2d_coordinate_block(rows2.clone()),molecule.with_2d_coordinate_block_with_params(rows2.clone(),&options2)),
                ("with_3d_coordinates",molecule.with_3d_coordinates(rows3.clone()),molecule.with_3d_coordinates_with_params(rows3.clone(),&replace)),
                ("with_added_3d_conformer",molecule.with_added_3d_conformer(rows3.clone()),molecule.with_added_3d_conformer_with_params(rows3.clone(),&options3)),
                ("with_only_3d_conformer",molecule.with_only_3d_conformer(rows3.clone()),molecule.with_only_3d_conformer_with_params(rows3.clone(),&options3)),
            ] {
                expected.set_item(format!("{name}_default"),default.unwrap().to_binary().unwrap())?;
                expected.set_item(name,configured.unwrap().to_binary().unwrap())?;
            }
            let options=ck::Coordinate2DParams {canonical_orientation:false,clear_existing_2d:false,flips_per_sample:1,samples:2,sample_seed:7,permute_degree_four:true,force_rdkit:true,use_ring_templates:false,use_path_angle_expansion:false,coordinate_map:std::collections::BTreeMap::from([(0,[0.,0.])])};
            let layout=molecule.with_2d_coordinates_with_params(&options).unwrap();
            let positions=molecule.with_atom_position_with_params(ck::AtomId::new(1),[2.,3.,4.],&ck::AtomPositionParams {conformer_id:Some(0)}).unwrap();
            let module=PyModule::new(py,"coordinate_family_test")?;
            module.add_class::<Molecule>()?;
            module.add_class::<Coordinate2DInputParams>()?;
            module.add_class::<crate::params::Coordinate3DInputParams>()?;
            module.add_class::<crate::params::Replace3DCoordinatesParams>()?;
            module.add_class::<crate::params::AtomPositionParams>()?;
            module.add_class::<Coordinate2DParams>()?;
            module.add_class::<crate::selectors::CoordinateZPolicy>()?;
            module.add_class::<crate::operation_errors::OperationError>()?;
            let locals=PyDict::new(py);
            locals.set_item("api",module)?;locals.set_item("before",before.clone())?;
            locals.set_item("expected",expected)?;locals.set_item("rows2",rows2)?;locals.set_item("rows3",rows3)?;
            locals.set_item("layout",layout.to_binary().unwrap())?;
            locals.set_item("positions",positions.to_binary().unwrap())?;
            locals.set_item("layout_params",Py::new(py,Coordinate2DParams::from_rust(py, options)?)?)?;
            locals.set_item("shared",Py::new(py,Molecule {inner:molecule.clone()})?)?;
            locals.set_item("peer",Py::new(py,Molecule {inner:molecule.clone()})?)?;
            py.run(c_str!(r#"
def fresh(): return api.Molecule.from_binary(bytes(before))
families = [
    ('with_2d_coordinate_block', 'set_2d_coordinates_', rows2, api.Coordinate2DInputParams(z_policy='require_zero'), {'z_policy':'require_zero'}),
    ('with_3d_coordinates', 'set_3d_coordinates_', rows3, api.Replace3DCoordinatesParams(conformer_id=0), {'conformer_id':0}),
    ('with_added_3d_conformer', 'add_3d_conformer_', rows3, api.Coordinate3DInputParams(is_3d=False), {'is_3d':False}),
    ('with_only_3d_conformer', 'set_only_3d_conformer_', rows3, api.Coordinate3DInputParams(is_3d=False), {'is_3d':False}),
]
for value_name, inplace_name, rows, params, kwargs in families:
    for name in [value_name, inplace_name]:
        inplace = name.endswith('_')
        alias = name[:-1] + '_with_params_' if inplace else name + '_with_params'
        for form in ['default','params','keywords','alias']:
            source = fresh()
            if form == 'default': output = getattr(source,name)(rows)
            elif form == 'params': output = getattr(source,name)(rows,params)
            elif form == 'keywords': output = getattr(source,name)(rows,**kwargs)
            else: output = getattr(source,alias)(rows,params)
            result = source if inplace else output
            assert bytes(result.to_binary()) == bytes(expected[value_name+'_default' if form=='default' else value_name]), (name,form)
            if not inplace: assert bytes(source.to_binary()) == bytes(before)
            elif name in ['add_3d_conformer_','set_only_3d_conformer_']: assert isinstance(output,int)
            else: assert output is None
        try: getattr(fresh(),name)(rows,params,**kwargs)
        except TypeError: pass
        else: raise AssertionError('mixed configuration accepted')
assert bytes(fresh().with_2d_coordinate_block(rows2,z_policy=api.CoordinateZPolicy.RequireZero).to_binary()) == bytes(expected['with_2d_coordinate_block'])
try: fresh().with_2d_coordinate_block(rows2,z_policy='unknown')
except ValueError: pass
else: raise AssertionError('unknown enum string accepted')
kwargs = dict(coordinate_map={0:[0.0,0.0]},canonical_orientation=False,clear_existing_2d=False,
    flips_per_sample=1,samples=2,sample_seed=7,permute_degree_four=True,force_rdkit=True,
    use_ring_templates=False,use_path_angle_expansion=False)
for name in ['with_2d_coordinates','add_2d_coordinates_']:
    inplace=name.endswith('_')
    alias=name[:-1]+'_with_params_' if inplace else name+'_with_params'
    for form in ['params','keywords','alias']:
        source=fresh()
        if form=='params': output=getattr(source,name)(layout_params)
        elif form=='keywords': output=getattr(source,name)(**kwargs)
        else: output=getattr(source,alias)(layout_params)
        result=source if inplace else output
        assert bytes(result.to_binary()) == bytes(layout)
        if inplace: assert output is None
        else: assert bytes(source.to_binary()) == bytes(before)
    source=fresh(); output=getattr(source,name)()
    assert (source if inplace else output).has_2d_coordinates()
for name in ['with_atom_position','set_atom_position_']:
    inplace=name.endswith('_')
    alias=name[:-1]+'_with_params_' if inplace else name+'_with_params'
    for form in ['default','params','keywords','alias']:
        source=fresh(); params=api.AtomPositionParams(conformer_id=0)
        if form=='default': output=getattr(source,name)(1,[2.0,3.0,4.0])
        elif form=='params': output=getattr(source,name)(1,[2.0,3.0,4.0],params)
        elif form=='keywords': output=getattr(source,name)(1,[2.0,3.0,4.0],conformer_id=0)
        else: output=getattr(source,alias)(1,[2.0,3.0,4.0],params)
        assert bytes((source if inplace else output).to_binary()) == bytes(positions)
        if not inplace: assert bytes(source.to_binary()) == bytes(before)
for inplace in [False,True]:
    source=fresh()
    result=source.clear_3d_conformers_() if inplace else source.with_cleared_3d_conformers()
    assert (source if inplace else result).num_3d_conformers()==0
    if not inplace: assert bytes(source.to_binary())==bytes(before)
source=fresh()
try: source.set_3d_coordinates_([[0.0,0.0]])
except api.OperationError as error:
    assert error.coordinate_input is not None and error.cause is error.__cause__
else: raise AssertionError('invalid native coordinates accepted')
assert bytes(source.to_binary()) == bytes(before)
shared.set_atom_position_(1,[2.0,3.0,4.0])
assert bytes(shared.to_binary()) == bytes(positions)
assert bytes(peer.to_binary()) == bytes(before)
"#),Some(&locals),None)?;
            assert_eq!(molecule.to_binary().unwrap(),before);
            Ok(())
        }).unwrap();
    }
}
