//! Actual IO-owned fields, without a binding-only aggregate group switch.
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "fixed_projection._native")]
pub struct BioMmcifWriteParams {
    pub(crate) inner: ck::BioMmcifWriteParams,
}

macro_rules! mmcif_parameters {
    (groups {$($group:ident),+ $(,)?}, flags {$($flag:ident),+ $(,)?}) => {
        paste::paste! {
            #[pymethods]
            impl BioMmcifWriteParams {
                #[new]
                #[pyo3(signature=(*, $($group=None,)+ $($flag=ck::BioMmcifWriteParams::default().$flag,)+ align_pairs=ck::BioMmcifWriteParams::default().align_pairs, align_loops=ck::BioMmcifWriteParams::default().align_loops))]
                fn new($($group: Option<bool>,)+ $($flag: bool,)+ align_pairs: u16, align_loops: u16) -> Self {
                    Self { inner: ck::BioMmcifWriteParams {
                        $($group: $group.unwrap_or_else(|| ck::BioMmcifWriteParams::default().$group),)+
                        $($flag,)+ align_pairs, align_loops,
                    } }
                }
                $(
                    #[getter]
                    fn $group(&self) -> bool { self.inner.$group }
                    #[setter]
                    fn [<set_ $group>](&mut self, value: Option<bool>) {
                        self.inner.$group = value.unwrap_or_else(|| ck::BioMmcifWriteParams::default().$group);
                    }
                )+
                $(
                    #[getter]
                    fn $flag(&self) -> bool { self.inner.$flag }
                    #[setter]
                    fn [<set_ $flag>](&mut self, value: bool) { self.inner.$flag = value; }
                )+
                #[getter]
                fn align_pairs(&self) -> u16 { self.inner.align_pairs }
                #[setter]
                fn set_align_pairs(&mut self, value: u16) { self.inner.align_pairs = value; }
                #[getter]
                fn align_loops(&self) -> u16 { self.inner.align_loops }
                #[setter]
                fn set_align_loops(&mut self, value: u16) { self.inner.align_loops = value; }
                fn __repr__(&self) -> String { format!("{:?}", self.inner) }
            }
        }
    };
}

mmcif_parameters!(
    groups {
        atoms,
        block_name,
        entry,
        database_status,
        author,
        cell,
        symmetry,
        entity,
        entity_poly,
        struct_ref,
        chem_comp,
        exptl,
        diffrn,
        reflns,
        refine,
        title_keywords,
        ncs,
        struct_asym,
        origx,
        struct_conf,
        struct_sheet,
        struct_biol,
        assembly,
        conn,
        cis,
        modres,
        scale,
        atom_type,
        entity_poly_seq,
        tls,
        software,
        group_pdb,
        auth_all,
    },
    flags {
        prefer_pairs,
        compact,
        misuse_hash
    }
);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::batch_fingerprint_values::TopologicalFingerprintParams;
    use pyo3::types::PyDict;

    #[test]
    fn actual_mmcif_defaults_and_topological_ignore_mask_are_forwarded() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let class = py.get_type::<BioMmcifWriteParams>();
            let params = class.call0()?.extract::<Py<BioMmcifWriteParams>>()?;
            assert_eq!(params.borrow(py).inner, ck::BioMmcifWriteParams::default());
            assert!(params.borrow(py).inner.atoms);
            assert!(!params.borrow(py).inner.auth_all);
            let structure = ck::BioStructure::from_pdb("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n").unwrap();
            let baseline = structure.to_mmcif().unwrap();
            params.bind(py).setattr("atoms", false)?;
            let configured = structure.to_mmcif_with_params(&params.borrow(py).inner).unwrap();
            assert_ne!(configured, baseline);
            let kwargs = PyDict::new(py);
            kwargs.set_item("atoms", false)?;
            let constructed = class.call((), Some(&kwargs))?.extract::<Py<BioMmcifWriteParams>>()?;
            assert_eq!(constructed.borrow(py).inner, params.borrow(py).inner);
            assert_eq!(structure.to_mmcif_with_params(&constructed.borrow(py).inner).unwrap(), configured);
            params.bind(py).setattr("atoms", py.None())?;
            assert_eq!(structure.to_mmcif_with_params(&params.borrow(py).inner).unwrap(), baseline);
            params.bind(py).setattr("align_pairs", 65535_u16)?;
            assert!(params.bind(py).setattr("align_pairs", 65536_u32).is_err());
            assert_eq!(params.borrow(py).inner.align_pairs, 65535);

            let class = py.get_type::<TopologicalFingerprintParams>();
            let params = class.call0()?.extract::<Py<TopologicalFingerprintParams>>()?;
            let molecule = ck::Molecule::from_smiles("CCO").unwrap();
            let baseline = molecule.fingerprint_topological_with_params(&params.borrow(py).inner).unwrap();
            assert!(baseline.num_on_bits() > 0);
            params.bind(py).setattr("ignore_atoms", vec![1_u32])?;
            let masked = molecule.fingerprint_topological_with_params(&params.borrow(py).inner).unwrap();
            assert_eq!(masked.num_on_bits(), 0);
            let kwargs = PyDict::new(py);
            kwargs.set_item("ignore_atoms", vec![1_u32])?;
            let constructed = class.call((), Some(&kwargs))?.extract::<Py<TopologicalFingerprintParams>>()?;
            assert_eq!(molecule.fingerprint_topological_with_params(&constructed.borrow(py).inner).unwrap(), masked);
            params.bind(py).setattr("ignore_atoms", Vec::<u32>::new())?;
            assert_eq!(molecule.fingerprint_topological_with_params(&params.borrow(py).inner).unwrap(), baseline);
            params.bind(py).setattr("ignore_atoms", py.None())?;
            assert!(params.borrow(py).inner.ignore_atoms.is_none());
            Ok(())
        }).unwrap();
    }
}
