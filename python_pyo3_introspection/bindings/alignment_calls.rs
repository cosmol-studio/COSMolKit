//! Coordinate alignment delegates all numerical work and effects to Rust.
use cosmolkit as ck;
pub(crate) fn conformer_keywords(
    atom_indices: Option<Vec<usize>>,
    conformer_ids: Option<Vec<usize>>,
    weights: Option<Vec<f64>>,
    reflect: Option<bool>,
    max_iterations: Option<u32>,
) -> ck::ConformerAlignmentParameters {
    let defaults = ck::ConformerAlignmentParameters::default();
    ck::ConformerAlignmentParameters::new(
        atom_indices,
        conformer_ids,
        weights,
        reflect.unwrap_or(defaults.reflect),
        max_iterations.unwrap_or(defaults.max_iterations),
    )
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    use pyo3::prelude::*;
    use pyo3::types::PyDict;
    #[test]
    fn alignment_value_inplace_families_keep_native_reports_coordinates_defaults_and_atomicity() {
        Python::initialize();
        Python::attach(|py|->PyResult<()>{
let mut builder=ck::Molecule::from_smiles("CCO").unwrap().to_builder();
builder.add_3d_conformer(vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0]]).unwrap();
builder.add_3d_conformer(vec![[5.0,2.0,1.0],[5.0,3.0,1.0],[4.0,3.0,1.0]]).unwrap();
let source=builder.build().unwrap();
let mut builder=ck::Molecule::from_smiles("CCO").unwrap().to_builder();
builder.add_3d_conformer(vec![[2.0,4.0,0.0],[2.0,5.0,0.0],[1.0,5.0,0.0]]).unwrap();
let reference=builder.build().unwrap();let before=source.to_binary().unwrap();let reference_before=reference.to_binary().unwrap();
let projected_reference=Py::new(py,Molecule{inner:reference.clone()})?;
{
let mut config=ck::AlignmentParameters::default();config.max_iterations=17;let params=Py::new(py,crate::collection_params::AlignmentParameters::from_rust(py, config.clone())?)?;
let keywords=PyDict::new(py);keywords.set_item("max_iterations",17)?;
let default=source.with_alignment_to(&reference).unwrap();let configured=source.with_alignment_to_with_params(&reference,&config).unwrap();
for mutates in [false,true]{for form in 0..4{
let actual=Py::new(py,Molecule{inner:source.clone()})?;let name=if form==3{if mutates{"align_to_with_params_"}else{"with_alignment_to_with_params"}}else{if mutates{"align_to_"}else{"with_alignment_to"}};
let output=match form{
0=>actual.bind(py).call_method1(name,(projected_reference.clone_ref(py),))?,
1|3=>actual.bind(py).call_method1(name,(projected_reference.clone_ref(py),params.clone_ref(py),))?,
_=>actual.bind(py).call_method(name,(projected_reference.clone_ref(py),),Some(&keywords))?,
};
let (want_molecule,want_report)=if form==0{&default}else{&configured};
let report=if mutates{assert_eq!(actual.borrow(py).inner.to_binary().unwrap(),want_molecule.to_binary().unwrap());output.extract::<Py<crate::alignment_values::AlignmentResult>>()?}else{let (molecule,report)=output.extract::<(Py<Molecule>,Py<crate::alignment_values::AlignmentResult>)>()?;assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(),want_molecule.to_binary().unwrap());assert_eq!(actual.borrow(py).inner.to_binary().unwrap(),before);report};
assert_eq!(&report.borrow(py).inner,want_report);
assert_eq!(report.borrow(py).inner.rmsd().to_bits(),want_report.rmsd().to_bits());assert_eq!(report.borrow(py).inner.transform().matrix().iter().flatten().map(|v|v.to_bits()).collect::<Vec<_>>(),want_report.transform().matrix().iter().flatten().map(|v|v.to_bits()).collect::<Vec<_>>());
}
let actual=Py::new(py,Molecule{inner:source.clone()})?;let error=actual.bind(py).call_method(if mutates{"align_to_"}else{"with_alignment_to"},(projected_reference.clone_ref(py),params.clone_ref(py),),Some(&keywords)).unwrap_err();assert!(error.is_instance_of::<pyo3::exceptions::PyTypeError>(py));assert_eq!(actual.borrow(py).inner.to_binary().unwrap(),before);
}
}
{
let mut config=ck::ConformerAlignmentParameters::default();config.max_iterations=17;let params=Py::new(py,crate::numerical_params::ConformerAlignmentParameters::from_rust(py, config.clone())?)?;
let keywords=PyDict::new(py);keywords.set_item("max_iterations",17)?;
let default=source.with_aligned_conformers().unwrap();let configured=source.with_aligned_conformers_with_params(&config).unwrap();
for mutates in [false,true]{for form in 0..4{
let actual=Py::new(py,Molecule{inner:source.clone()})?;let name=if form==3{if mutates{"align_conformers_with_params_"}else{"with_aligned_conformers_with_params"}}else{if mutates{"align_conformers_"}else{"with_aligned_conformers"}};
let output=match form{
0=>actual.bind(py).call_method1(name,())?,
1|3=>actual.bind(py).call_method1(name,(params.clone_ref(py),))?,
_=>actual.bind(py).call_method(name,(),Some(&keywords))?,
};
let (want_molecule,want_report)=if form==0{&default}else{&configured};
let report=if mutates{assert_eq!(actual.borrow(py).inner.to_binary().unwrap(),want_molecule.to_binary().unwrap());output.extract::<Py<crate::alignment_values::ConformerAlignmentReport>>()?}else{let (molecule,report)=output.extract::<(Py<Molecule>,Py<crate::alignment_values::ConformerAlignmentReport>)>()?;assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(),want_molecule.to_binary().unwrap());assert_eq!(actual.borrow(py).inner.to_binary().unwrap(),before);report};
assert_eq!(&report.borrow(py).inner,want_report);
assert_eq!(report.borrow(py).inner.rmsds().iter().map(|v|v.to_bits()).collect::<Vec<_>>(),want_report.rmsds().iter().map(|v|v.to_bits()).collect::<Vec<_>>());
}
let actual=Py::new(py,Molecule{inner:source.clone()})?;let error=actual.bind(py).call_method(if mutates{"align_conformers_"}else{"with_aligned_conformers"},(params.clone_ref(py),),Some(&keywords)).unwrap_err();assert!(error.is_instance_of::<pyo3::exceptions::PyTypeError>(py));assert_eq!(actual.borrow(py).inner.to_binary().unwrap(),before);
}
}

let mut config=ck::AlignmentParameters::default();config.probe_conformer_id=i32::MAX;
let mut native=source.clone();let expected=native.align_to_with_params_(&reference,&config).unwrap_err();
let actual=Py::new(py,Molecule{inner:source.clone()})?;
let keywords=PyDict::new(py);keywords.set_item("probe_conformer_id",i32::MAX)?;
let error=actual.bind(py).call_method("align_to_",(projected_reference.clone_ref(py),),Some(&keywords)).unwrap_err();
assert!(error.is_instance_of::<crate::operation_errors::OperationError>(py));assert_eq!(error.value(py).str()?.to_str()?,expected.to_string());assert_eq!(actual.borrow(py).inner.to_binary().unwrap(),before);
assert_eq!(source.to_binary().unwrap(),before);assert_eq!(projected_reference.borrow(py).inner.to_binary().unwrap(),reference_before);
Ok(())}).unwrap();
    }
}
