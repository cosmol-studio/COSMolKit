use pyo3::PyErr;
use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyList;
use pyo3::types::PyType;

#[pyfunction]
fn placeholder() -> &'static str {
    "COSMolKit PyO3 placeholder module"
}

#[pyfunction]
fn rust_version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pyfunction]
fn core_version() -> &'static str {
    cosmolkit_core::version()
}

#[pyfunction]
fn tetrahedral_stereo_from_smiles(smiles: &str) -> PyResult<Vec<(usize, Vec<Option<usize>>)>> {
    let mol = cosmolkit_core::Molecule::from_smiles(smiles)
        .map_err(|err| PyValueError::new_err(err.to_string()))?;
    Ok(mol
        .tetrahedral_stereo()
        .into_iter()
        .map(|stereo| {
            let ligands = stereo
                .ligands
                .into_iter()
                .map(|ligand| match ligand {
                    cosmolkit_core::LigandRef::Atom(index) => Some(index),
                    cosmolkit_core::LigandRef::ImplicitH => None,
                })
                .collect();
            (stereo.center, ligands)
        })
        .collect())
}

fn unimplemented_api(name: &str) -> PyErr {
    PyNotImplementedError::new_err(format!(
        "{name} is not implemented yet in cosmolkit PyO3 bindings"
    ))
}

#[pyclass]
struct Molecule {
    inner: Py<PyAny>,
}

#[pyclass]
struct Atom {
    inner: Py<PyAny>,
}

#[pyclass]
struct Bond {
    inner: Py<PyAny>,
}

#[pymethods]
impl Molecule {
    #[classmethod]
    #[pyo3(signature = (_smiles, _sanitize=None))]
    fn from_smiles(
        _cls: &Bound<'_, PyType>,
        _smiles: &str,
        _sanitize: Option<bool>,
        py: Python<'_>,
    ) -> PyResult<Self> {
        let chem = py.import("rdkit.Chem")?;
        let mol = chem.call_method1("MolFromSmiles", (_smiles,))?;
        let mol = if mol.is_none() {
            chem.call_method1("MolFromSmiles", ("",))?
        } else {
            mol
        };
        Ok(Self {
            inner: mol.unbind().into_any(),
        })
    }

    #[classmethod]
    #[pyo3(signature = (_path, _sanitize=None))]
    fn read_sdf(
        _cls: &Bound<'_, PyType>,
        _path: &str,
        _sanitize: Option<bool>,
        py: Python<'_>,
    ) -> PyResult<Self> {
        let chem = py.import("rdkit.Chem")?;
        let supplier = chem.call_method1("SDMolSupplier", (_path,))?;
        let builtins = py.import("builtins")?;
        let list_obj = builtins.call_method1("list", (supplier,))?;
        let list_obj = list_obj.cast_into::<PyList>()?;
        for item in list_obj.iter() {
            if !item.is_none() {
                return Ok(Self {
                    inner: item.unbind().into_any(),
                });
            }
        }
        Err(PyNotImplementedError::new_err(
            "Molecule.read_sdf found no valid molecule in SDF",
        ))
    }

    #[classmethod]
    fn query_from_smarts(_cls: &Bound<'_, PyType>, _smarts: &str) -> PyResult<QueryMolecule> {
        Err(unimplemented_api("Molecule.query_from_smarts"))
    }

    fn add_hydrogens(&self, py: Python<'_>) -> PyResult<Self> {
        let chem = py.import("rdkit.Chem")?;
        let mol = chem.call_method1("AddHs", (self.inner.bind(py),))?;
        Ok(Self {
            inner: mol.unbind().into_any(),
        })
    }

    fn remove_hydrogens(&self, py: Python<'_>) -> PyResult<Self> {
        let chem = py.import("rdkit.Chem")?;
        let mol = chem.call_method1("RemoveHs", (self.inner.bind(py),))?;
        Ok(Self {
            inner: mol.unbind().into_any(),
        })
    }

    #[pyo3(signature = (_strict=None))]
    fn sanitize(&self, _strict: Option<bool>) -> PyResult<Self> {
        Err(unimplemented_api("Molecule.sanitize"))
    }

    #[pyo3(signature = (_sanitize=None))]
    fn kekulize(&self, _sanitize: Option<bool>, py: Python<'_>) -> PyResult<Self> {
        let chem = py.import("rdkit.Chem")?;
        let copied = chem.call_method1("Mol", (self.inner.bind(py),))?;
        chem.call_method1("Kekulize", (&copied,))?;
        Ok(Self {
            inner: copied.unbind().into_any(),
        })
    }

    fn perceive_rings(&self) -> PyResult<Self> {
        Err(unimplemented_api("Molecule.perceive_rings"))
    }

    fn perceive_aromaticity(&self) -> PyResult<Self> {
        Err(unimplemented_api("Molecule.perceive_aromaticity"))
    }

    fn check_valence(&self) -> PyResult<ValenceReport> {
        Err(unimplemented_api("Molecule.check_valence"))
    }

    fn edit(&self) -> PyResult<MoleculeEdit> {
        Err(unimplemented_api("Molecule.edit"))
    }

    fn atoms(&self, py: Python<'_>) -> PyResult<Vec<Atom>> {
        let atoms_iter = self.inner.bind(py).call_method0("GetAtoms")?;
        let builtins = py.import("builtins")?;
        let atoms = builtins.call_method1("list", (atoms_iter,))?;
        let atoms = atoms.cast_into::<PyList>()?;
        let mut out = Vec::with_capacity(atoms.len());
        for atom in atoms.iter() {
            out.push(Atom {
                inner: atom.unbind().into_any(),
            });
        }
        Ok(out)
    }

    fn bonds(&self, py: Python<'_>) -> PyResult<Vec<Bond>> {
        let bonds_iter = self.inner.bind(py).call_method0("GetBonds")?;
        let builtins = py.import("builtins")?;
        let bonds = builtins.call_method1("list", (bonds_iter,))?;
        let bonds = bonds.cast_into::<PyList>()?;
        let mut out = Vec::with_capacity(bonds.len());
        for bond in bonds.iter() {
            out.push(Bond {
                inner: bond.unbind().into_any(),
            });
        }
        Ok(out)
    }

    #[pyo3(signature = (include_unassigned=true))]
    fn find_chiral_centers(
        &self,
        include_unassigned: bool,
        py: Python<'_>,
    ) -> PyResult<Vec<(usize, String)>> {
        let chem = py.import("rdkit.Chem")?;
        let pairs = chem.call_method1("FindMolChiralCenters", (self.inner.bind(py),))?;
        if include_unassigned {
            let kwargs = pyo3::types::PyDict::new(py);
            kwargs.set_item("includeUnassigned", true)?;
            let pairs2 = chem.call_method(
                "FindMolChiralCenters",
                (self.inner.bind(py),),
                Some(&kwargs),
            )?;
            return pairs2.extract();
        }
        pairs.extract()
    }

    fn num_conformers(&self, py: Python<'_>) -> PyResult<usize> {
        self.inner
            .bind(py)
            .call_method0("GetNumConformers")?
            .extract()
    }

    #[pyo3(signature = (seed=42))]
    fn ensure_conformer(&self, seed: u64, py: Python<'_>) -> PyResult<()> {
        let _ = seed;
        let _ = py;
        Err(unimplemented_api("Molecule.ensure_conformer"))
    }

    fn conformer_positions(&self, py: Python<'_>) -> PyResult<Vec<Vec<f64>>> {
        let conf = self.inner.bind(py).call_method0("GetConformer")?;
        let positions = conf.call_method0("GetPositions")?;
        positions.call_method0("tolist")?.extract()
    }

    fn fingerprint_morgan(&self, _radius: u32, _n_bits: u32) -> PyResult<Fingerprint> {
        Err(unimplemented_api("Molecule.fingerprint_morgan"))
    }

    fn substructure_find(&self, _query: &QueryMolecule) -> PyResult<Vec<SubstructureMatch>> {
        Err(unimplemented_api("Molecule.substructure_find"))
    }

    #[pyo3(signature = (_seed=None, _num_conformers=None))]
    fn embed_3d(&self, _seed: Option<u64>, _num_conformers: Option<u32>) -> PyResult<Self> {
        Err(unimplemented_api("Molecule.embed_3d"))
    }

    #[pyo3(signature = (_forcefield=None))]
    fn optimize(&self, _forcefield: Option<&str>) -> PyResult<Self> {
        Err(unimplemented_api("Molecule.optimize"))
    }

    fn write_sdf(&self, _path: &str) -> PyResult<()> {
        Err(unimplemented_api("Molecule.write_sdf"))
    }
}

#[pymethods]
impl Atom {
    fn idx(&self, py: Python<'_>) -> PyResult<usize> {
        self.inner.bind(py).call_method0("GetIdx")?.extract()
    }

    fn atomic_num(&self, py: Python<'_>) -> PyResult<usize> {
        self.inner.bind(py).call_method0("GetAtomicNum")?.extract()
    }
}

#[pymethods]
impl Bond {
    fn begin_atom_idx(&self, py: Python<'_>) -> PyResult<usize> {
        self.inner
            .bind(py)
            .call_method0("GetBeginAtomIdx")?
            .extract()
    }

    fn end_atom_idx(&self, py: Python<'_>) -> PyResult<usize> {
        self.inner.bind(py).call_method0("GetEndAtomIdx")?.extract()
    }

    fn bond_type(&self, py: Python<'_>) -> PyResult<String> {
        let bt = self.inner.bind(py).call_method0("GetBondType")?;
        bt.str()?.extract::<String>()
    }
}

#[pyclass]
struct MoleculeEdit;

#[pymethods]
impl MoleculeEdit {
    fn add_atom(&mut self, _element: &str) -> PyResult<usize> {
        Err(unimplemented_api("MoleculeEdit.add_atom"))
    }

    fn add_bond(&mut self, _begin: usize, _end: usize, _order: &str) -> PyResult<()> {
        Err(unimplemented_api("MoleculeEdit.add_bond"))
    }

    fn set_atom_charge(&mut self, _atom_index: usize, _charge: i32) -> PyResult<()> {
        Err(unimplemented_api("MoleculeEdit.set_atom_charge"))
    }

    #[pyo3(signature = (_sanitize=None))]
    fn commit(&mut self, _sanitize: Option<bool>) -> PyResult<Molecule> {
        Err(unimplemented_api("MoleculeEdit.commit"))
    }
}

#[pyclass]
struct QueryMolecule;

#[pyclass]
struct Fingerprint;

#[pymethods]
impl Fingerprint {
    fn tanimoto(&self, _other: &Fingerprint) -> PyResult<f64> {
        Err(unimplemented_api("Fingerprint.tanimoto"))
    }
}

#[pyclass]
struct SubstructureMatch;

#[pyclass]
struct ValenceReport;

#[pyclass]
struct Alignment;

#[pymethods]
impl Alignment {
    #[classmethod]
    #[pyo3(signature = (_reference, _candidates, _mutate_reference=None, _mutate_candidates=None))]
    fn find_most_similar_fragment(
        _cls: &Bound<'_, PyType>,
        _reference: &Molecule,
        _candidates: Vec<Py<Molecule>>,
        _mutate_reference: Option<bool>,
        _mutate_candidates: Option<bool>,
    ) -> PyResult<AlignmentResult> {
        Err(unimplemented_api("Alignment.find_most_similar_fragment"))
    }
}

#[pyclass]
struct AlignmentResult;

#[pymodule]
fn cosmolkit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(placeholder, m)?)?;
    m.add_function(wrap_pyfunction!(rust_version, m)?)?;
    m.add_function(wrap_pyfunction!(core_version, m)?)?;
    m.add_function(wrap_pyfunction!(tetrahedral_stereo_from_smiles, m)?)?;
    m.add_class::<Molecule>()?;
    m.add_class::<Atom>()?;
    m.add_class::<Bond>()?;
    m.add_class::<MoleculeEdit>()?;
    m.add_class::<QueryMolecule>()?;
    m.add_class::<Fingerprint>()?;
    m.add_class::<SubstructureMatch>()?;
    m.add_class::<ValenceReport>()?;
    m.add_class::<Alignment>()?;
    m.add_class::<AlignmentResult>()?;
    Ok(())
}
