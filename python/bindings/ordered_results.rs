use crate::native::{Fingerprint, Molecule};
use cosmolkit as ck;
use pyo3::prelude::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn readonly_stereo_queries_preserve_native_order_labels_and_molecule_state() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            for text in ["C[C@H](F)Cl", "CC(F)Cl", "CC"] {
                let native = ck::Molecule::from_smiles(text).unwrap();
                let before = native.to_binary().unwrap();
                let expected_stereo = native.tetrahedral_stereo().unwrap();
                let expected_assigned = native.find_chiral_centers(false).unwrap();
                let expected_all = native.find_chiral_centers(true).unwrap();
                let molecule = Py::new(py, Molecule { inner: native })?;
                let rows = molecule
                    .bind(py)
                    .call_method0("tetrahedral_stereo")?
                    .extract::<Vec<TetrahedralStereo>>()?;
                assert_eq!(
                    rows.into_iter()
                        .map(|value| value.inner)
                        .collect::<Vec<_>>(),
                    expected_stereo
                );
                assert!(
                    molecule
                        .bind(py)
                        .call_method0("perceive_stereochemistry")?
                        .is_none()
                );
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method0("find_chiral_centers")?
                        .extract::<Vec<(usize, String)>>()?,
                    expected_assigned
                );
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method1("find_chiral_centers", (true,))?
                        .extract::<Vec<(usize, String)>>()?,
                    expected_all
                );
                assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
            }
            Ok(())
        })
        .unwrap();
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object, eq)]
    #[derive(Clone, Copy, PartialEq)]
    /// Tetrahedral ligand represented by an explicit atom or an implicit hydrogen.
    pub struct LigandRef {
        inner: ck::LigandRef,
    }
}
#[pymethods]
impl LigandRef {
    #[new]
    #[pyo3(signature = (atom=None))]
    fn py_new(atom: Option<usize>) -> Self {
        match atom {
            Some(value) => Self::atom(value),
            None => Self::IMPLICIT_HYDROGEN,
        }
    }
    #[classattr]
    #[pyo3(name = "ImplicitHydrogen")]
    const IMPLICIT_HYDROGEN: Self = Self {
        inner: ck::LigandRef::ImplicitHydrogen,
    };
    #[staticmethod]
    #[pyo3(name = "Atom")]
    /// Represent a tetrahedral ligand using a zero-based atom identifier.
    fn atom(value: usize) -> Self {
        Self {
            inner: ck::LigandRef::Atom(ck::AtomId::new(value)),
        }
    }
    #[getter]
    /// Ligand atom identifier, or None for an implicit hydrogen.
    fn atom_id(&self) -> Option<usize> {
        match self.inner {
            ck::LigandRef::Atom(value) => Some(value.index()),
            ck::LigandRef::ImplicitHydrogen => None,
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Tetrahedral center and its four ordered ligand references.
    pub struct TetrahedralStereo {
        pub(crate) inner: ck::TetrahedralStereo,
    }
}
#[pymethods]
impl TetrahedralStereo {
    #[new]
    fn new(center: usize, ligands: [LigandRef; 4]) -> Self {
        Self {
            inner: ck::TetrahedralStereo {
                center: ck::AtomId::new(center),
                ligands: ligands.map(|value| value.inner),
            },
        }
    }
    #[getter]
    /// Zero-based atom identifier of the tetrahedral center.
    fn center(&self) -> usize {
        self.inner.center.index()
    }
    #[getter]
    /// Four ordered ligand references defining the tetrahedral configuration.
    fn ligands(&self) -> [LigandRef; 4] {
        self.inner.ligands.map(|inner| LigandRef { inner })
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Atom and bond correspondences returned by native substructure matching.
    pub struct MatchResult {
        pub(crate) inner: ck::MatchResult,
    }
}
#[pymethods]
impl MatchResult {
    /// Query-to-target atom pairs represented by this match.
    fn atom_pairs(&self) -> Vec<(usize, usize)> {
        self.inner.atom_pairs()
    }
    #[new]
    fn new(atom_mapping: Vec<usize>, bond_mapping: Vec<usize>) -> Self {
        Self {
            inner: ck::MatchResult::new(atom_mapping, bond_mapping),
        }
    }
    #[getter]
    /// Target atom indices in query-atom order.
    fn atom_mapping(&self) -> Vec<usize> {
        self.inner.atom_mapping.clone()
    }
    #[getter]
    /// Target bond indices in query-bond order.
    fn bond_mapping(&self) -> Vec<usize> {
        self.inner.bond_mapping.clone()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Layered fingerprint and optional atom contribution counts.
    pub struct LayeredFingerprintResult {
        pub(crate) inner: ck::LayeredFingerprintResult,
    }
}
#[pymethods]
impl LayeredFingerprintResult {
    #[new]
    #[pyo3(signature=(fingerprint, atom_counts))]
    fn new(fingerprint: &Fingerprint, atom_counts: Option<Vec<u32>>) -> Self {
        Self {
            inner: ck::LayeredFingerprintResult {
                fingerprint: fingerprint.inner.clone(),
                atom_counts,
            },
        }
    }
    #[getter]
    /// Generated bit fingerprint.
    fn fingerprint(&self) -> Fingerprint {
        Fingerprint {
            inner: self.inner.fingerprint().clone(),
        }
    }
    #[getter]
    /// Optional per-atom contribution counts requested during fingerprint generation.
    fn atom_counts(&self) -> Option<Vec<u32>> {
        self.inner.atom_counts().map(<[u32]>::to_vec)
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
/// Atom-pair atom code together with the prepared molecule used to compute it.
pub struct AtomPairAtomCodeResult {
        pub(crate) inner: ck::AtomPairAtomCodeResult,
    }
}
#[pymethods]
impl AtomPairAtomCodeResult {
    #[new]
    fn new(code: u32, molecule: &Molecule) -> Self {
        Self {
            inner: ck::AtomPairAtomCodeResult {
                code,
                molecule: molecule.inner.clone(),
            },
        }
    }
    #[getter]
    /// Encoded atom-pair atom invariant.
    fn code(&self) -> u32 {
        self.inner.code
    }
    #[getter]
    /// Prepared molecule associated with this code; returned with copy-on-write value semantics.
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule.clone(),
        }
    }
}
