use pyo3::PyErr;
use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
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
    // Representation contract is defined in tetrahedral_stereo_representation.md.
    let mol = cosmolkit_core::Molecule::from_smiles(smiles)
        .map_err(|err| PyValueError::new_err(err.to_string()))?;
    Ok(to_python_tetrahedral_stereo(&mol))
}

fn to_python_tetrahedral_stereo(mol: &cosmolkit_core::Molecule) -> Vec<(usize, Vec<Option<usize>>)> {
    mol.tetrahedral_stereo()
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
        .collect()
}

fn unimplemented_api(name: &str) -> PyErr {
    PyNotImplementedError::new_err(format!(
        "{name} is not implemented yet in cosmolkit PyO3 bindings"
    ))
}

fn bond_order_name(order: cosmolkit_core::BondOrder) -> &'static str {
    match order {
        cosmolkit_core::BondOrder::Null => "ZERO",
        cosmolkit_core::BondOrder::Single => "SINGLE",
        cosmolkit_core::BondOrder::Double => "DOUBLE",
        cosmolkit_core::BondOrder::Triple => "TRIPLE",
        cosmolkit_core::BondOrder::Quadruple => "QUADRUPLE",
        cosmolkit_core::BondOrder::Aromatic => "AROMATIC",
        cosmolkit_core::BondOrder::Dative => "DATIVE",
    }
}

#[pyclass]
#[derive(Clone)]
struct Molecule {
    inner: cosmolkit_core::Molecule,
}

#[pyclass]
#[derive(Clone)]
struct Atom {
    idx: usize,
    atomic_num: usize,
}

#[pyclass]
#[derive(Clone)]
struct Bond {
    begin_atom_idx: usize,
    end_atom_idx: usize,
    bond_type: String,
}

#[pymethods]
impl Molecule {
    #[classmethod]
    #[pyo3(signature = (smiles, sanitize=None))]
    fn from_smiles(
        _cls: &Bound<'_, PyType>,
        smiles: &str,
        sanitize: Option<bool>,
    ) -> PyResult<Self> {
        let _ = sanitize;
        let mol = cosmolkit_core::Molecule::from_smiles(smiles)
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Self { inner: mol })
    }

    #[classmethod]
    #[pyo3(signature = (path, sanitize=None))]
    fn read_sdf(_cls: &Bound<'_, PyType>, path: &str, sanitize: Option<bool>) -> PyResult<Self> {
        let _ = (path, sanitize);
        Err(unimplemented_api("Molecule.read_sdf"))
    }

    #[classmethod]
    fn query_from_smarts(_cls: &Bound<'_, PyType>, smarts: &str) -> PyResult<QueryMolecule> {
        let _ = smarts;
        Err(unimplemented_api("Molecule.query_from_smarts"))
    }

    fn add_hydrogens(&self) -> PyResult<Self> {
        let mut out = self.inner.clone();
        cosmolkit_core::add_hydrogens_in_place(&mut out)
            .map_err(|err| PyValueError::new_err(format!("add_hydrogens failed: {err:?}")))?;
        Ok(Self { inner: out })
    }

    fn remove_hydrogens(&self) -> PyResult<Self> {
        Err(unimplemented_api("Molecule.remove_hydrogens"))
    }

    #[pyo3(signature = (strict=None))]
    fn sanitize(&self, strict: Option<bool>) -> PyResult<Self> {
        let _ = strict;
        Err(unimplemented_api("Molecule.sanitize"))
    }

    #[pyo3(signature = (sanitize=None))]
    fn kekulize(&self, sanitize: Option<bool>) -> PyResult<Self> {
        let _ = sanitize;
        let mut out = self.inner.clone();
        cosmolkit_core::kekulize::kekulize_in_place(&mut out)
            .map_err(|err| PyValueError::new_err(format!("kekulize failed: {err:?}")))?;
        Ok(Self { inner: out })
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

    fn atoms(&self) -> Vec<Atom> {
        self.inner
            .atoms
            .iter()
            .map(|atom| Atom {
                idx: atom.index,
                atomic_num: atom.atomic_num as usize,
            })
            .collect()
    }

    fn bonds(&self) -> Vec<Bond> {
        self.inner
            .bonds
            .iter()
            .map(|bond| Bond {
                begin_atom_idx: bond.begin_atom,
                end_atom_idx: bond.end_atom,
                bond_type: bond_order_name(bond.order).to_string(),
            })
            .collect()
    }

    #[pyo3(signature = (include_unassigned=true))]
    fn find_chiral_centers(&self, include_unassigned: bool) -> Vec<(usize, String)> {
        self.inner
            .atoms
            .iter()
            .filter_map(|atom| match atom.chiral_tag {
                cosmolkit_core::ChiralTag::Unspecified => {
                    if include_unassigned {
                        Some((atom.index, "?".to_string()))
                    } else {
                        None
                    }
                }
                cosmolkit_core::ChiralTag::TetrahedralCw => {
                    Some((atom.index, "CHI_TETRAHEDRAL_CW".to_string()))
                }
                cosmolkit_core::ChiralTag::TetrahedralCcw => {
                    Some((atom.index, "CHI_TETRAHEDRAL_CCW".to_string()))
                }
            })
            .collect()
    }

    fn tetrahedral_stereo(&self) -> Vec<(usize, Vec<Option<usize>>)> {
        // Representation contract is defined in tetrahedral_stereo_representation.md.
        to_python_tetrahedral_stereo(&self.inner)
    }

    fn num_conformers(&self) -> usize {
        0
    }

    #[pyo3(signature = (seed=42))]
    fn ensure_conformer(&self, seed: u64) -> PyResult<()> {
        let _ = seed;
        Err(unimplemented_api("Molecule.ensure_conformer"))
    }

    fn conformer_positions(&self) -> PyResult<Vec<Vec<f64>>> {
        Err(unimplemented_api("Molecule.conformer_positions"))
    }

    #[pyo3(signature = (radius, n_bits))]
    fn fingerprint_morgan(&self, radius: u32, n_bits: u32) -> PyResult<Fingerprint> {
        let _ = (radius, n_bits);
        Err(unimplemented_api("Molecule.fingerprint_morgan"))
    }

    fn substructure_find(&self, query: &QueryMolecule) -> PyResult<Vec<SubstructureMatch>> {
        let _ = query;
        Err(unimplemented_api("Molecule.substructure_find"))
    }

    #[pyo3(signature = (seed=None, num_conformers=None))]
    fn embed_3d(&self, seed: Option<u64>, num_conformers: Option<u32>) -> PyResult<Self> {
        let _ = (seed, num_conformers);
        Err(unimplemented_api("Molecule.embed_3d"))
    }

    #[pyo3(signature = (forcefield=None))]
    fn optimize(&self, forcefield: Option<&str>) -> PyResult<Self> {
        let _ = forcefield;
        Err(unimplemented_api("Molecule.optimize"))
    }

    fn write_sdf(&self, path: &str) -> PyResult<()> {
        let _ = path;
        Err(unimplemented_api("Molecule.write_sdf"))
    }
}

#[pymethods]
impl Atom {
    fn idx(&self) -> usize {
        self.idx
    }

    fn atomic_num(&self) -> usize {
        self.atomic_num
    }
}

#[pymethods]
impl Bond {
    fn begin_atom_idx(&self) -> usize {
        self.begin_atom_idx
    }

    fn end_atom_idx(&self) -> usize {
        self.end_atom_idx
    }

    fn bond_type(&self) -> String {
        self.bond_type.clone()
    }
}

#[pyclass]
struct MoleculeEdit;

#[pymethods]
impl MoleculeEdit {
    fn add_atom(&mut self, element: &str) -> PyResult<usize> {
        let _ = element;
        Err(unimplemented_api("MoleculeEdit.add_atom"))
    }

    #[pyo3(signature = (begin, end, order))]
    fn add_bond(&mut self, begin: usize, end: usize, order: &str) -> PyResult<()> {
        let _ = (begin, end, order);
        Err(unimplemented_api("MoleculeEdit.add_bond"))
    }

    fn set_atom_charge(&mut self, atom_index: usize, charge: i32) -> PyResult<()> {
        let _ = (atom_index, charge);
        Err(unimplemented_api("MoleculeEdit.set_atom_charge"))
    }

    #[pyo3(signature = (sanitize=None))]
    fn commit(&mut self, sanitize: Option<bool>) -> PyResult<Molecule> {
        let _ = sanitize;
        Err(unimplemented_api("MoleculeEdit.commit"))
    }
}

#[pyclass]
struct QueryMolecule;

#[pyclass]
struct Fingerprint;

#[pymethods]
impl Fingerprint {
    fn tanimoto(&self, other: &Fingerprint) -> PyResult<f64> {
        let _ = other;
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
    #[pyo3(signature = (reference, candidates, mutate_reference=None, mutate_candidates=None))]
    fn find_most_similar_fragment(
        _cls: &Bound<'_, PyType>,
        reference: &Molecule,
        candidates: Vec<Py<Molecule>>,
        mutate_reference: Option<bool>,
        mutate_candidates: Option<bool>,
    ) -> PyResult<AlignmentResult> {
        let _ = (reference, candidates, mutate_reference, mutate_candidates);
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
