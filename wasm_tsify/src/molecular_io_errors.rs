//! Public IO roots retain non-Clone native children and their actual sources.
use crate::errors::{enum_error, exception};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "IoError | XyzReadError | XyzWriteError | OperationError | Mol2ReadError | Mol2PostError | MolWriteError | PdbWriteError | SdfError | FromUtf8Error | BioPdbReadError | BioMmcifReadError | BioMoleculeError"
    )]
    #[derive(Clone)]
    pub type MolecularIoCause;
}
fn source(owner: &ReadonlyValue<ck::MolecularIoError>) -> Option<MolecularIoCause> {
    use ck::MolecularIoError as E;
    macro_rules! field {
        ($pattern:pat => $field:expr,$module:ident::$wrapper:ident) => {
            crate::$module::$wrapper::from_owner(owner.child(|e| match e {
                $pattern => $field,
                _ => unreachable!("immutable selected IO variant"),
            }))
        };
    }
    Some(
        match owner.as_ref() {
            E::Io { .. } => field!(E::Io{source,..}=>source,external_errors::IoError),
            E::XyzRead(_) => field!(E::XyzRead(v)=>v,xyz_errors::XyzReadError),
            E::XyzWrite(_) => field!(E::XyzWrite(v)=>v,xyz_errors::XyzWriteError),
            E::Construction(_) => field!(E::Construction(v)=>v,operation_errors::OperationError),
            E::Mol2Read(_) => field!(E::Mol2Read(v)=>v,mol2_errors::Mol2ReadError),
            E::Mol2Post(_) => field!(E::Mol2Post(v)=>v,mol2_errors::Mol2PostError),
            E::MolWrite(_) => field!(E::MolWrite(v)=>v,mol_write_errors::MolWriteError),
            E::PdbWrite(_) => field!(E::PdbWrite(v)=>v,io_write_errors::PdbWriteError),
            E::Sdf(_) => field!(E::Sdf(v)=>v,sdf_errors::SdfError),
            #[cfg(feature = "cap-bio")]
            E::PdbRead(_) => field!(E::PdbRead(v)=>v,pdb_read_errors::BioPdbReadError),
            #[cfg(feature = "cap-bio")]
            E::MmcifRead(_) => field!(E::MmcifRead(v)=>v,bio_read_errors::BioMmcifReadError),
            #[cfg(feature = "cap-bio")]
            E::BioMolecule(_) => {
                field!(E::BioMolecule(v)=>v,bio_conversion_errors::BioMoleculeError)
            }
            E::OutputUtf8 { .. } => {
                field!(E::OutputUtf8{source,..}=>source,external_errors::FromUtf8Error)
            }
            E::NoRecord { .. } | E::Parameter { .. } => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(MolecularIoError,MolecularIoErrorKind,ck::MolecularIoError,{
    #[cfg(feature = "cap-bio")]
PdbRead=>ck::MolecularIoError::PdbRead(_),#[cfg(feature = "cap-bio")]
MmcifRead=>ck::MolecularIoError::MmcifRead(_),#[cfg(feature = "cap-bio")]
BioMolecule=>ck::MolecularIoError::BioMolecule(_),
    Io=>ck::MolecularIoError::Io{..},XyzRead=>ck::MolecularIoError::XyzRead(_),XyzWrite=>ck::MolecularIoError::XyzWrite(_),Construction=>ck::MolecularIoError::Construction(_),
    Mol2Read=>ck::MolecularIoError::Mol2Read(_),Mol2Post=>ck::MolecularIoError::Mol2Post(_),MolWrite=>ck::MolecularIoError::MolWrite(_),PdbWrite=>ck::MolecularIoError::PdbWrite(_),
    Sdf=>ck::MolecularIoError::Sdf(_),NoRecord=>ck::MolecularIoError::NoRecord{..},OutputUtf8=>ck::MolecularIoError::OutputUtf8{..},Parameter=>ck::MolecularIoError::Parameter{..},
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<MolecularIoCause>{source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn path(&self)->Option<String>{match self.inner.as_ref(){ck::MolecularIoError::Io{path,..}=>Some(path.display().to_string()),_=>None}}
    #[wasm_bindgen(getter)] pub fn format(&self)->Option<String>{match self.inner.as_ref(){ck::MolecularIoError::NoRecord{format}|ck::MolecularIoError::OutputUtf8{format,..}=>Some((*format).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn name(&self)->Option<String>{match self.inner.as_ref(){ck::MolecularIoError::Parameter{name,..}=>Some((*name).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn detail(&self)->Option<String>{match self.inner.as_ref(){ck::MolecularIoError::Parameter{detail,..}=>Some((*detail).into()),_=>None}}
},source=MolecularIoCause,cause=source);
