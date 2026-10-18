//! Native structural-to-molecular conversion, including transparent errors.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::bio_errors::BioStructureError;
use crate::error_framework::{Source, error_projection};
use crate::hydrogen_errors::{HydrogenError, HydrogenSource};
use crate::io_source_errors::ErrorRef;
use crate::model_edit_errors::AtomPropertyError;
use crate::model_errors::{CoordinateValidationError, StereoGroupError, TopologyValidationError};
use crate::operation_errors::OperationError;
use crate::property_numeric_errors::PropertyUIntReadError;
use crate::retained_error::RetainedError;
use crate::sanitize_errors::{SanitizeError, SanitizeSource};
use crate::stereo_errors::{StereoError, StereoSource};
use crate::topology_edit_errors::TopologyEditError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

type PdbPostSource = Source<
    Py<AtomPropertyError>,
    Source<
        Py<TopologyValidationError>,
        Source<
            Py<CoordinateValidationError>,
            Source<Py<ValenceError>, Source<Py<TopologyEditError>, Py<PropertyUIntReadError>>>,
        >,
    >,
>;
error_projection!(PdbPostprocessError,PdbPostprocessErrorKind,ck::PdbPostprocessError,"pdb",{
    AtomProperty(_),Unsupported(_),InvalidState(_),Topology(_),Coordinates(_),Valence(_),TopologyEdit(_),PropertyUInt(_),
},cause=|inner:&ck::PdbPostprocessError| match inner{
    ck::PdbPostprocessError::AtomProperty(value)=>Some(AtomPropertyError::from_rust(value.clone())),
    ck::PdbPostprocessError::Topology(value)=>Some(TopologyValidationError::from_rust(value.clone())),
    ck::PdbPostprocessError::Coordinates(value)=>Some(CoordinateValidationError::from_rust(value.clone())),
    ck::PdbPostprocessError::Valence(value)=>Some(ValenceError::from_rust(value.clone())),
    ck::PdbPostprocessError::TopologyEdit(value)=>Some(TopologyEditError::from_rust(value.clone())),
    ck::PdbPostprocessError::PropertyUInt(value)=>Some(PropertyUIntReadError::from_rust(value.clone())),
    ck::PdbPostprocessError::Unsupported(_)|ck::PdbPostprocessError::InvalidState(_)=>None,
},source=Option<PdbPostSource>,methods={
    #[getter] fn value(&self)->Option<&str>{match &self.inner{
        ck::PdbPostprocessError::Unsupported(value)=>Some(value),ck::PdbPostprocessError::InvalidState(value)=>Some(value),_=>None}}
});

// thiserror(transparent) delegates Error::source to the stored error's source.
// The stored enum payload is separately readable, never invented as a cause.
fn conversion_payload(inner: &ck::BioMoleculeConversionError) -> PyErr {
    match inner {
        ck::BioMoleculeConversionError::Structure(value) => {
            BioStructureError::from_rust(value.clone())
        }
        ck::BioMoleculeConversionError::Topology(value) => {
            TopologyValidationError::from_rust(value.clone())
        }
        ck::BioMoleculeConversionError::Chemistry(value) => {
            PdbPostprocessError::from_rust(value.clone())
        }
        ck::BioMoleculeConversionError::Sanitize(value) => SanitizeError::from_rust(value.clone()),
        ck::BioMoleculeConversionError::Hydrogens(value) => HydrogenError::from_rust(value.clone()),
        ck::BioMoleculeConversionError::Valence(value) => ValenceError::from_rust(value.clone()),
        ck::BioMoleculeConversionError::Stereo(value) => StereoError::from_rust(value.clone()),
    }
}
type ConversionSource = Source<
    PdbPostSource,
    Source<
        SanitizeSource,
        Source<
            HydrogenSource,
            Source<
                StereoSource,
                Source<
                    Py<TopologyValidationError>,
                    Source<Py<StereoGroupError>, Py<TemplateAttachmentOrderError>>,
                >,
            >,
        >,
    >,
>;
type ConversionPayload = Source<
    Py<BioStructureError>,
    Source<
        Py<TopologyValidationError>,
        Source<
            Py<PdbPostprocessError>,
            Source<
                Py<SanitizeError>,
                Source<Py<HydrogenError>, Source<Py<ValenceError>, Py<StereoError>>>,
            >,
        >,
    >,
>;
error_projection!(retained BioMoleculeConversionError,BioMoleculeConversionErrorKind,ck::BioMoleculeConversionError,"bio",{
    Structure(_),Topology(_),Chemistry(_),Sanitize(_),Hydrogens(_),Valence(_),Stereo(_),
},cause=|inner:&RetainedError<ck::BioMoleculeConversionError>|Python::attach(|py|conversion_payload(inner.as_ref()).cause(py)),
source=Option<ConversionSource>,methods={
    #[getter] #[doc=r#"{"role":"error_payload"}"#]
    fn value(&self,py:Python<'_>)->PyResult<ConversionPayload>{
        conversion_payload(self.inner.as_ref()).value(py).extract().map_err(Into::into)
    }
});
error_projection!(retained BioMoleculeError,BioMoleculeErrorKind,ck::BioMoleculeError,"bio",{
    Conversion(_),Construction(_),
},cause=|inner:&RetainedError<ck::BioMoleculeError>|{
    let node=inner.node().source().expect("native BIO molecule error has a source");
    Some(match inner.as_ref(){
        ck::BioMoleculeError::Conversion(_)=>BioMoleculeConversionError::from_ref(node),
        ck::BioMoleculeError::Construction(_)=>OperationError::from_ref(node),
    })
},source=Option<Source<Py<BioMoleculeConversionError>,Py<OperationError>>>);

pub(crate) fn keywords(
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
    flavor: Option<u32>,
    proximity_bonding: Option<bool>,
) -> crate::params::BioMoleculeParams {
    let defaults = ck::BioMoleculeParams::default();
    crate::params::BioMoleculeParams {
        inner: ck::BioMoleculeParams {
            sanitize: sanitize.unwrap_or(defaults.sanitize),
            remove_hs: remove_hs.unwrap_or(defaults.remove_hs),
            flavor: flavor.unwrap_or(defaults.flavor),
            proximity_bonding: proximity_bonding.unwrap_or(defaults.proximity_bonding),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bio_structure::BioStructure;
    use crate::native::Molecule;
    use crate::protein_values::Protein;
    use pyo3::types::PyDict;
    use std::sync::Arc;

    #[test]
    fn bio_molecular_conversion_matches_native_defaults_options_and_transparent_causes() {
        Python::initialize();
        Python::attach(|py|->PyResult<()> {
            let text="ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nEND\n";
            let bio=ck::BioStructure::from_pdb(text).unwrap();
            let protein=bio.protein().unwrap();
            let original=bio.to_mmcif().unwrap();
            let defaults=ck::BioMoleculeParams::default();
            let configured=ck::BioMoleculeParams{sanitize:false,remove_hs:false,flavor:1,proximity_bonding:false};
            let params=Py::new(py,crate::params::BioMoleculeParams{inner:configured})?;
            let kwargs=PyDict::new(py);
            kwargs.set_item("sanitize",false)?;kwargs.set_item("remove_hs",false)?;
            kwargs.set_item("flavor",1)?;kwargs.set_item("proximity_bonding",false)?;
            let python_bio=Py::new(py,BioStructure{inner:bio.clone()})?;
            let python_protein=Py::new(py,Protein{inner:Arc::new(protein.clone())})?;
            let snapshot=|value:&Py<PyAny>|->PyResult<String>{
                if value.bind(py).is_instance_of::<BioStructure>() {
                    Ok(value.extract::<Py<BioStructure>>(py)?.borrow(py).inner.to_mmcif().unwrap())
                } else {
                    Ok(value.extract::<Py<Protein>>(py)?.borrow(py).inner.as_bio_structure().to_mmcif().unwrap())
                }
            };
            for (receiver,expected_default,expected_configured) in [
                (python_bio.into_any(),bio.to_molecule_with_params(&defaults).unwrap(),bio.to_molecule_with_params(&configured).unwrap()),
                (python_protein.into_any(),protein.to_molecule_with_params(&defaults).unwrap(),protein.to_molecule_with_params(&configured).unwrap()),
            ] {
                let before=snapshot(&receiver)?;
                let actual=receiver.bind(py).call_method0("to_molecule")?.extract::<Py<Molecule>>()?;
                assert_eq!(actual.borrow(py).inner.to_binary().unwrap(),expected_default.to_binary().unwrap());
                for actual in [receiver.bind(py).call_method1("to_molecule",(&params,))?,receiver.bind(py).call_method("to_molecule",(),Some(&kwargs))?,receiver.bind(py).call_method1("to_molecule_with_params",(&params,))?] {
                    assert_eq!(actual.extract::<Py<Molecule>>()?.borrow(py).inner.to_binary().unwrap(),expected_configured.to_binary().unwrap());
                }
                assert!(receiver.bind(py).call_method("to_molecule",(&params,),Some(&kwargs)).unwrap_err().is_instance_of::<pyo3::exceptions::PyTypeError>(py));
                assert_eq!(snapshot(&receiver)?,before);
            }
            assert_eq!(bio.to_mmcif().unwrap(),original);
            let values=vec![
                ck::BioMoleculeConversionError::Structure(ck::BioStructureError::RowIndexTooLarge{value:usize::MAX}),
                ck::BioMoleculeConversionError::Topology(ck::TopologyValidationError::AdjacencyMismatch),
                ck::BioMoleculeConversionError::Chemistry(ck::PdbPostprocessError::Unsupported("actual unsupported field")),
                ck::BioMoleculeConversionError::Chemistry(ck::PdbPostprocessError::Topology(ck::TopologyValidationError::AdjacencyMismatch)),
                ck::BioMoleculeConversionError::Sanitize(ck::SanitizeError::InvalidOperations{bits:0x80000000,unknown_bits:0x80000000}),
                ck::BioMoleculeConversionError::Hydrogens(ck::HydrogenError::OnlyOnAtomOutOfRange{atom:ck::AtomId::new(9),atom_count:1}),
                ck::BioMoleculeConversionError::Valence(ck::ValenceError::AtomOutOfRange{atom:ck::AtomId::new(9),atom_count:1}),
                ck::BioMoleculeConversionError::Stereo(ck::StereoError::DuplicateConformerId{id:9}),
            ];
            for native in values {
                let display=native.to_string();let has_source=native.source().is_some();
                let error=BioMoleculeError::from_rust(ck::BioMoleculeError::Conversion(native));
                assert_eq!(error.value(py).str()?.to_str()?,display);
                let conversion=error.value(py).getattr("source")?.unbind();
                assert!(conversion.bind(py).is_instance_of::<BioMoleculeConversionError>());
                assert!(conversion.bind(py).is(error.value(py).getattr("__cause__")?.as_any()));
                assert_eq!(!conversion.bind(py).getattr("source")?.is_none(),has_source);
                let value=conversion.bind(py).getattr("value")?.unbind();
                assert_eq!(value.bind(py).str()?.to_str()?,display);
                assert!(value.bind(py).setattr("variant",py.None()).is_err());
                drop(error);drop(conversion);
                assert_eq!(value.bind(py).str()?.to_str()?,display);
            }
            Ok(())
        }).unwrap();
    }
}
