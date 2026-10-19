//! Detached enhanced-stereo groups, using the canonical model's local validation.
use crate::chemical_values::enum_value;
use crate::model_errors::StereoGroupError;
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(StereoGroupKind, { Absolute, Or, And });

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct StereoGroup {
    pub(crate) inner: ck::StereoGroup,
}

#[pymethods]
impl StereoGroup {
    #[staticmethod]
    fn new(kind: StereoGroupKind, atoms: Vec<usize>, bonds: Vec<usize>) -> PyResult<Self> {
        ck::StereoGroup::new(
            kind.into(),
            atoms.into_iter().map(ck::AtomId::new).collect(),
            bonds.into_iter().map(ck::BondId::new).collect(),
        )
        .map(|inner| Self { inner })
        .map_err(StereoGroupError::from_rust)
    }
    fn with_id(&self, id: u32) -> Self {
        Self {
            inner: self.inner.clone().with_id(id),
        }
    }
    fn with_write_id(&self, write_id: u32) -> Self {
        Self {
            inner: self.inner.clone().with_write_id(write_id),
        }
    }
    fn id(&self) -> Option<u32> {
        self.inner.id()
    }
    fn write_id(&self) -> u32 {
        self.inner.write_id()
    }
    fn kind(&self) -> StereoGroupKind {
        self.inner.kind().into()
    }
    fn atoms(&self) -> Vec<usize> {
        self.inner.atoms().iter().map(|id| id.index()).collect()
    }
    fn bonds(&self) -> Vec<usize> {
        self.inner.bonds().iter().map(|id| id.index()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::MoleculeBuilder;
    use crate::model_errors::StereoGroupErrorKind;

    #[test]
    fn enhanced_stereo_values_keep_ids_native_validation_and_detached_builder_reads() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            for kind in [
                StereoGroupKind::Absolute,
                StereoGroupKind::Or,
                StereoGroupKind::And,
            ] {
                let group = Py::new(py, StereoGroup::new(kind, vec![2, 0], vec![1])?)?;
                let changed = group
                    .bind(py)
                    .call_method1("with_id", (0,))?
                    .call_method1("with_write_id", (u32::MAX,))?;
                assert!(group.bind(py).call_method0("id")?.is_none());
                assert_eq!(
                    group.bind(py).call_method0("write_id")?.extract::<u32>()?,
                    0
                );
                assert_eq!(changed.call_method0("id")?.extract::<u32>()?, 0);
                assert_eq!(
                    changed.call_method0("write_id")?.extract::<u32>()?,
                    u32::MAX
                );
                assert!(changed.call_method0("kind")?.extract::<StereoGroupKind>()? == kind);
                assert_eq!(
                    changed.call_method0("atoms")?.extract::<Vec<usize>>()?,
                    vec![2, 0]
                );
                assert_eq!(
                    changed.call_method0("bonds")?.extract::<Vec<usize>>()?,
                    vec![1]
                );
            }
            let error = py
                .get_type::<StereoGroup>()
                .call_method1(
                    "new",
                    (StereoGroupKind::And, vec![0, 0], Vec::<usize>::new()),
                )
                .unwrap_err();
            assert!(error.is_instance_of::<StereoGroupError>(py));
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<StereoGroupErrorKind>()?
                    == StereoGroupErrorKind::DuplicateAtom
            );

            let source = ck::Molecule::from_smiles("C[C@H](F)Cl |o1:1|").unwrap();
            let before = source.to_binary().unwrap();
            let expected = source.to_builder().stereo_groups().to_vec();
            assert_eq!(expected.len(), 1);
            let builder = Py::new(
                py,
                MoleculeBuilder {
                    inner: source.to_builder(),
                },
            )?;
            let values = builder.bind(py).call_method0("stereo_groups")?;
            drop(builder);
            assert_eq!(
                values.get_item(0)?.extract::<PyRef<StereoGroup>>()?.inner,
                expected[0]
            );
            assert_eq!(source.to_binary().unwrap(), before);
            Ok(())
        })
        .unwrap();
    }
}
