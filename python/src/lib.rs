use std::fs::File;
use std::io::{BufReader, Write};
use std::path::Path;

use cosmolkit_core::io::molblock::{self, SdfFormat};
use cosmolkit_core::io::sdf::SdfReader;
use pyo3::PyErr;
use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyType};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::define_stub_info_gatherer;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

#[pyfunction]
fn placeholder() -> &'static str {
    "COSMolKit PyO3 module"
}

#[pyfunction]
fn rust_version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pyfunction]
fn core_version() -> &'static str {
    cosmolkit_core::version()
}

fn to_python_tetrahedral_stereo(
    mol: &cosmolkit_core::Molecule,
) -> Vec<(usize, Vec<Option<usize>>)> {
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
        cosmolkit_core::BondOrder::Null => "UNSPECIFIED",
        cosmolkit_core::BondOrder::Single => "SINGLE",
        cosmolkit_core::BondOrder::Double => "DOUBLE",
        cosmolkit_core::BondOrder::Triple => "TRIPLE",
        cosmolkit_core::BondOrder::Quadruple => "QUADRUPLE",
        cosmolkit_core::BondOrder::Aromatic => "AROMATIC",
        cosmolkit_core::BondOrder::Dative => "DATIVE",
    }
}

fn chiral_tag_name(tag: cosmolkit_core::ChiralTag) -> &'static str {
    match tag {
        cosmolkit_core::ChiralTag::Unspecified => "CHI_UNSPECIFIED",
        cosmolkit_core::ChiralTag::TetrahedralCw => "CHI_TETRAHEDRAL_CW",
        cosmolkit_core::ChiralTag::TetrahedralCcw => "CHI_TETRAHEDRAL_CCW",
    }
}

fn bond_direction_name(direction: cosmolkit_core::BondDirection) -> &'static str {
    match direction {
        cosmolkit_core::BondDirection::None => "NONE",
        cosmolkit_core::BondDirection::EndUpRight => "ENDUPRIGHT",
        cosmolkit_core::BondDirection::EndDownRight => "ENDDOWNRIGHT",
    }
}

fn bond_stereo_name(stereo: cosmolkit_core::BondStereo) -> &'static str {
    match stereo {
        cosmolkit_core::BondStereo::None => "STEREONONE",
        cosmolkit_core::BondStereo::Any => "STEREOANY",
        cosmolkit_core::BondStereo::Cis => "STEREOCIS",
        cosmolkit_core::BondStereo::Trans => "STEREOTRANS",
    }
}

fn parse_sdf_format(format: Option<&str>) -> PyResult<SdfFormat> {
    match format.map(|s| s.to_ascii_lowercase()) {
        None => Ok(SdfFormat::Auto),
        Some(v) if v == "auto" => Ok(SdfFormat::Auto),
        Some(v) if v == "v2000" || v == "v2k" => Ok(SdfFormat::V2000),
        Some(v) if v == "v3000" || v == "v3k" => Ok(SdfFormat::V3000),
        Some(v) => Err(PyValueError::new_err(format!(
            "unsupported SDF format '{v}', expected one of: auto, v2000, v3000"
        ))),
    }
}

fn atomic_number_from_element(element: &str) -> Option<u8> {
    match element {
        "H" => Some(1),
        "B" => Some(5),
        "C" => Some(6),
        "N" => Some(7),
        "O" => Some(8),
        "F" => Some(9),
        "P" => Some(15),
        "S" => Some(16),
        "Cl" => Some(17),
        "Br" => Some(35),
        "I" => Some(53),
        "Na" => Some(11),
        "K" => Some(19),
        "Li" => Some(3),
        "Mg" => Some(12),
        "Ca" => Some(20),
        "Fe" => Some(26),
        "Cu" => Some(29),
        "Zn" => Some(30),
        "Si" => Some(14),
        "Al" => Some(13),
        "*" => Some(0),
        _ => None,
    }
}

fn py_method<'py>(obj: &Bound<'py, PyAny>, method: &str) -> PyResult<Bound<'py, PyAny>> {
    obj.call_method0(method)
        .map_err(|err| PyValueError::new_err(format!("from_rdkit failed calling {method}: {err}")))
}

fn py_method_index<'py>(
    obj: &Bound<'py, PyAny>,
    method: &str,
    index: usize,
) -> PyResult<Bound<'py, PyAny>> {
    obj.call_method1(method, (index,))
        .map_err(|err| PyValueError::new_err(format!("from_rdkit failed calling {method}: {err}")))
}

fn py_method_extract<T>(obj: &Bound<'_, PyAny>, method: &str) -> PyResult<T>
where
    for<'a> T: FromPyObject<'a, 'a>,
{
    py_method(obj, method)?.extract::<T>().map_err(|_| {
        PyValueError::new_err(format!("from_rdkit failed extracting result from {method}"))
    })
}

fn py_method_str(obj: &Bound<'_, PyAny>, method: &str) -> PyResult<String> {
    let value = py_method(obj, method)?;
    Ok(value
        .str()
        .map_err(|err| {
            PyValueError::new_err(format!("from_rdkit failed stringifying {method}: {err}"))
        })?
        .to_string_lossy()
        .into_owned())
}

fn rdkit_chiral_tag_from_name(name: &str) -> PyResult<cosmolkit_core::ChiralTag> {
    match name {
        "CHI_UNSPECIFIED" => Ok(cosmolkit_core::ChiralTag::Unspecified),
        "CHI_TETRAHEDRAL_CW" => Ok(cosmolkit_core::ChiralTag::TetrahedralCw),
        "CHI_TETRAHEDRAL_CCW" => Ok(cosmolkit_core::ChiralTag::TetrahedralCcw),
        other => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported atom chiral tag '{other}'"
        ))),
    }
}

fn rdkit_bond_order_from_name(name: &str) -> PyResult<cosmolkit_core::BondOrder> {
    match name {
        "UNSPECIFIED" | "ZERO" => Ok(cosmolkit_core::BondOrder::Null),
        "SINGLE" => Ok(cosmolkit_core::BondOrder::Single),
        "DOUBLE" => Ok(cosmolkit_core::BondOrder::Double),
        "TRIPLE" => Ok(cosmolkit_core::BondOrder::Triple),
        "QUADRUPLE" => Ok(cosmolkit_core::BondOrder::Quadruple),
        "AROMATIC" => Ok(cosmolkit_core::BondOrder::Aromatic),
        "DATIVE" | "DATIVEL" | "DATIVER" => Ok(cosmolkit_core::BondOrder::Dative),
        other => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported bond type '{other}'"
        ))),
    }
}

fn rdkit_bond_direction_from_name(name: &str) -> PyResult<cosmolkit_core::BondDirection> {
    match name {
        "NONE" => Ok(cosmolkit_core::BondDirection::None),
        "ENDUPRIGHT" => Ok(cosmolkit_core::BondDirection::EndUpRight),
        "ENDDOWNRIGHT" => Ok(cosmolkit_core::BondDirection::EndDownRight),
        other => Err(PyValueError::new_err(format!(
            "from_rdkit unsupported bond direction '{other}'"
        ))),
    }
}

fn clone_for_python_value_api(mol: &cosmolkit_core::Molecule) -> cosmolkit_core::Molecule {
    // Transitional implementation for the immutable-first Python API:
    // public methods return a new Molecule value, but today we still realize
    // that by cloning the Rust molecule eagerly before mutation.
    //
    // Planned follow-up: replace this eager clone with copy-on-write storage
    // in core so unchanged topology/conformer/property blocks stay shared.
    mol.clone()
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct Molecule {
    inner: cosmolkit_core::Molecule,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct Atom {
    idx: usize,
    atomic_num: usize,
    formal_charge: i8,
    chiral_tag: String,
    isotope: Option<u16>,
    is_aromatic: bool,
    explicit_hydrogens: usize,
    no_implicit: bool,
    num_radical_electrons: usize,
    degree: usize,
    explicit_valence: Option<usize>,
    implicit_hydrogens: Option<usize>,
    total_num_hs: Option<usize>,
    total_valence: Option<usize>,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object)]
#[derive(Clone)]
struct Bond {
    idx: usize,
    begin_atom_idx: usize,
    end_atom_idx: usize,
    bond_type: String,
    bond_dir: String,
    stereo: String,
    stereo_atoms: Vec<usize>,
    is_aromatic: bool,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
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
    #[pyo3(signature = (rdmol, sanitize=None))]
    fn from_rdkit(
        _cls: &Bound<'_, PyType>,
        rdmol: &Bound<'_, PyAny>,
        sanitize: Option<bool>,
    ) -> PyResult<Self> {
        rdmol.py().import("rdkit.Chem").map_err(|err| {
            PyValueError::new_err(format!(
                "Molecule.from_rdkit requires rdkit to be installed and importable: {err}"
            ))
        })?;
        let _ = sanitize;
        let atom_count: usize = py_method_extract(rdmol, "GetNumAtoms")?;
        let bond_count: usize = py_method_extract(rdmol, "GetNumBonds")?;
        let mut mol = cosmolkit_core::Molecule::new();

        for idx in 0..atom_count {
            let atom = py_method_index(rdmol, "GetAtomWithIdx", idx)?;
            let atomic_num_raw: usize = py_method_extract(&atom, "GetAtomicNum")?;
            let atomic_num = u8::try_from(atomic_num_raw).map_err(|_| {
                PyValueError::new_err(format!(
                    "from_rdkit atom {idx} atomic number out of u8 range: {atomic_num_raw}"
                ))
            })?;
            let formal_charge_raw: i32 = py_method_extract(&atom, "GetFormalCharge")?;
            let formal_charge = i8::try_from(formal_charge_raw).map_err(|_| {
                PyValueError::new_err(format!(
                    "from_rdkit atom {idx} formal charge out of i8 range: {formal_charge_raw}"
                ))
            })?;
            let explicit_h_raw: usize = py_method_extract(&atom, "GetNumExplicitHs")?;
            let explicit_hydrogens = u8::try_from(explicit_h_raw).map_err(|_| {
                PyValueError::new_err(format!(
                    "from_rdkit atom {idx} explicit H count out of u8 range: {explicit_h_raw}"
                ))
            })?;
            let radical_raw: usize = py_method_extract(&atom, "GetNumRadicalElectrons")?;
            let num_radical_electrons = u8::try_from(radical_raw).map_err(|_| {
                PyValueError::new_err(format!(
                    "from_rdkit atom {idx} radical electron count out of u8 range: {radical_raw}"
                ))
            })?;
            let isotope_raw: u16 = py_method_extract(&atom, "GetIsotope")?;
            let chiral_tag = rdkit_chiral_tag_from_name(&py_method_str(&atom, "GetChiralTag")?)?;

            mol.add_atom(cosmolkit_core::Atom {
                index: 0,
                atomic_num,
                is_aromatic: py_method_extract(&atom, "GetIsAromatic")?,
                formal_charge,
                explicit_hydrogens,
                no_implicit: py_method_extract(&atom, "GetNoImplicit")?,
                num_radical_electrons,
                chiral_tag,
                isotope: (isotope_raw != 0).then_some(isotope_raw),
                atom_map_num: None,
            });
        }

        for idx in 0..bond_count {
            let bond = py_method_index(rdmol, "GetBondWithIdx", idx)?;
            let begin_atom: usize = py_method_extract(&bond, "GetBeginAtomIdx")?;
            let end_atom: usize = py_method_extract(&bond, "GetEndAtomIdx")?;
            if begin_atom >= atom_count || end_atom >= atom_count {
                return Err(PyValueError::new_err(format!(
                    "from_rdkit bond {idx} atom index out of range: {begin_atom}-{end_atom}"
                )));
            }
            let is_aromatic: bool = py_method_extract(&bond, "GetIsAromatic")?;
            let order = rdkit_bond_order_from_name(&py_method_str(&bond, "GetBondType")?)?;
            let direction = rdkit_bond_direction_from_name(&py_method_str(&bond, "GetBondDir")?)?;
            let stereo_atoms: Vec<usize> =
                py_method(&bond, "GetStereoAtoms")?
                    .extract()
                    .map_err(|err| {
                        PyValueError::new_err(format!(
                            "from_rdkit failed extracting result from GetStereoAtoms: {err}"
                        ))
                    })?;

            mol.add_bond(cosmolkit_core::Bond {
                index: 0,
                begin_atom,
                end_atom,
                order,
                is_aromatic,
                direction,
                stereo: cosmolkit_core::BondStereo::None,
                stereo_atoms,
            });
        }

        cosmolkit_core::assign_double_bond_stereo_from_directions(&mut mol);
        mol.rebuild_adjacency();
        Ok(Self { inner: mol })
    }

    #[classmethod]
    #[pyo3(signature = (path, sanitize=None))]
    fn read_sdf(_cls: &Bound<'_, PyType>, path: &str, sanitize: Option<bool>) -> PyResult<Self> {
        let _ = sanitize;
        let file = File::open(path)
            .map_err(|e| PyValueError::new_err(format!("read_sdf open failed: {e}")))?;
        let mut reader = SdfReader::new(BufReader::new(file));
        let Some(record) = reader
            .next_record()
            .map_err(|e| PyValueError::new_err(format!("read_sdf parse failed: {e:?}")))?
        else {
            return Err(PyValueError::new_err("read_sdf found no molecule record"));
        };
        Ok(Self {
            inner: record.molecule,
        })
    }

    fn add_hydrogens(&self) -> PyResult<Self> {
        let mut out = clone_for_python_value_api(&self.inner);
        cosmolkit_core::add_hydrogens_in_place(&mut out)
            .map_err(|err| PyValueError::new_err(format!("add_hydrogens failed: {err:?}")))?;
        Ok(Self { inner: out })
    }

    fn remove_hydrogens(&self) -> PyResult<Self> {
        let mut out = clone_for_python_value_api(&self.inner);
        cosmolkit_core::remove_hydrogens_in_place(&mut out)
            .map_err(|err| PyValueError::new_err(format!("remove_hydrogens failed: {err:?}")))?;
        Ok(Self { inner: out })
    }

    #[pyo3(signature = (sanitize=None))]
    fn kekulize(&self, sanitize: Option<bool>) -> PyResult<Self> {
        let _ = sanitize;
        let mut out = clone_for_python_value_api(&self.inner);
        cosmolkit_core::kekulize::kekulize_in_place(&mut out, false)
            .map_err(|err| PyValueError::new_err(format!("kekulize failed: {err:?}")))?;
        Ok(Self { inner: out })
    }

    fn atoms(&self) -> Vec<Atom> {
        let assignment =
            cosmolkit_core::assign_valence(&self.inner, cosmolkit_core::ValenceModel::RdkitLike)
                .ok();
        let mut degrees = vec![0usize; self.inner.atoms.len()];
        for bond in &self.inner.bonds {
            degrees[bond.begin_atom] += 1;
            degrees[bond.end_atom] += 1;
        }
        self.inner
            .atoms
            .iter()
            .map(|atom| Atom {
                idx: atom.index,
                atomic_num: atom.atomic_num as usize,
                formal_charge: atom.formal_charge,
                chiral_tag: chiral_tag_name(atom.chiral_tag).to_string(),
                isotope: atom.isotope,
                is_aromatic: atom.is_aromatic,
                explicit_hydrogens: atom.explicit_hydrogens as usize,
                no_implicit: atom.no_implicit,
                num_radical_electrons: atom.num_radical_electrons as usize,
                degree: degrees[atom.index],
                explicit_valence: assignment
                    .as_ref()
                    .map(|v| v.explicit_valence[atom.index] as usize),
                implicit_hydrogens: assignment
                    .as_ref()
                    .map(|v| v.implicit_hydrogens[atom.index] as usize),
                total_num_hs: assignment.as_ref().map(|v| {
                    atom.explicit_hydrogens as usize + v.implicit_hydrogens[atom.index] as usize
                }),
                total_valence: assignment.as_ref().map(|v| {
                    v.explicit_valence[atom.index] as usize
                        + v.implicit_hydrogens[atom.index] as usize
                }),
            })
            .collect()
    }

    fn bonds(&self) -> Vec<Bond> {
        self.inner
            .bonds
            .iter()
            .map(|bond| Bond {
                idx: bond.index,
                begin_atom_idx: bond.begin_atom,
                end_atom_idx: bond.end_atom,
                bond_type: bond_order_name(bond.order).to_string(),
                bond_dir: bond_direction_name(bond.direction).to_string(),
                stereo: bond_stereo_name(bond.stereo).to_string(),
                stereo_atoms: bond.stereo_atoms.clone(),
                is_aromatic: bond.is_aromatic,
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
        to_python_tetrahedral_stereo(&self.inner)
    }

    fn compute_2d_coords(&self) -> PyResult<Self> {
        let mut out = clone_for_python_value_api(&self.inner);
        out.compute_2d_coords()
            .map_err(|err| PyValueError::new_err(format!("compute_2d_coords failed: {err}")))?;
        Ok(Self { inner: out })
    }

    fn num_conformers(&self) -> usize {
        if self.inner.coords_2d().is_some() {
            1
        } else {
            0
        }
    }

    fn conformer_positions(&self) -> PyResult<Vec<Vec<f64>>> {
        let Some(coords) = self.inner.coords_2d() else {
            return Err(PyValueError::new_err(
                "no 2D coordinates present; call compute_2d_coords() first",
            ));
        };
        Ok(coords.iter().map(|p| vec![p.x, p.y, 0.0]).collect())
    }

    fn dg_bounds_matrix(&self) -> PyResult<Vec<Vec<f64>>> {
        self.inner.dg_bounds_matrix().map_err(|err| {
            PyValueError::new_err(format!("Molecule.dg_bounds_matrix failed: {err}"))
        })
    }

    #[pyo3(signature = (isomeric_smiles=true))]
    fn to_smiles(&self, isomeric_smiles: bool) -> PyResult<String> {
        self.inner.to_smiles(isomeric_smiles).map_err(|err| {
            PyNotImplementedError::new_err(format!(
                "Molecule.to_smiles(isomeric_smiles={isomeric_smiles}) is not implemented yet: {err}"
            ))
        })
    }

    #[pyo3(signature = (path, format=None))]
    fn write_sdf(&self, path: &str, format: Option<&str>) -> PyResult<()> {
        let fmt = parse_sdf_format(format)?;
        let block = molblock::mol_to_2d_sdf_record(&self.inner, fmt)
            .map_err(|err| PyValueError::new_err(format!("write_sdf failed: {err}")))?;
        let mut f = File::create(path)
            .map_err(|e| PyValueError::new_err(format!("write_sdf create failed: {e}")))?;
        f.write_all(block.as_bytes())
            .map_err(|e| PyValueError::new_err(format!("write_sdf write failed: {e}")))?;
        Ok(())
    }

    #[pyo3(signature = (format=None))]
    fn to_sdf_string(&self, format: Option<&str>) -> PyResult<String> {
        let fmt = parse_sdf_format(format)?;
        molblock::mol_to_2d_sdf_record(&self.inner, fmt)
            .map_err(|err| PyValueError::new_err(format!("to_sdf_string failed: {err}")))
    }

    #[pyo3(signature = (directory, file_name=None, format=None))]
    fn write_sdf_to_directory(
        &self,
        directory: &str,
        file_name: Option<&str>,
        format: Option<&str>,
    ) -> PyResult<String> {
        let dir = Path::new(directory);
        if !dir.exists() {
            return Err(PyValueError::new_err(format!(
                "directory does not exist: {directory}"
            )));
        }
        if !dir.is_dir() {
            return Err(PyValueError::new_err(format!(
                "path is not a directory: {directory}"
            )));
        }
        let name = file_name.unwrap_or("molecule.sdf");
        if name.trim().is_empty() {
            return Err(PyValueError::new_err("file_name cannot be empty"));
        }
        let output = dir.join(name);
        let output_str = output
            .to_str()
            .ok_or_else(|| PyValueError::new_err("output path is not valid UTF-8"))?;
        self.write_sdf(output_str, format)?;
        Ok(output_str.to_string())
    }

    fn edit(&self) -> MoleculeEdit {
        MoleculeEdit {
            // Transitional implementation: the edit context starts from a
            // cloned working molecule today. Once core adopts copy-on-write,
            // this should share storage initially and detach only on mutation.
            working: clone_for_python_value_api(&self.inner),
        }
    }

    #[pyo3(signature = (strict=None))]
    fn sanitize(&self, strict: Option<bool>) -> PyResult<Self> {
        let _ = strict;
        Err(unimplemented_api("Molecule.sanitize"))
    }

    fn perceive_rings(&self) -> PyResult<Self> {
        Err(unimplemented_api("Molecule.perceive_rings"))
    }

    fn perceive_aromaticity(&self) -> PyResult<Self> {
        Err(unimplemented_api("Molecule.perceive_aromaticity"))
    }

    fn __len__(&self) -> usize {
        self.inner.atoms.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "Molecule(num_atoms={}, num_bonds={}, has_2d_coords={})",
            self.inner.atoms.len(),
            self.inner.bonds.len(),
            self.inner.coords_2d().is_some()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Atom {
    fn idx(&self) -> usize {
        self.idx
    }
    fn atomic_num(&self) -> usize {
        self.atomic_num
    }
    fn formal_charge(&self) -> i8 {
        self.formal_charge
    }
    fn chiral_tag(&self) -> String {
        self.chiral_tag.clone()
    }
    fn isotope(&self) -> Option<u16> {
        self.isotope
    }
    fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }
    fn explicit_hydrogens(&self) -> usize {
        self.explicit_hydrogens
    }
    fn no_implicit(&self) -> bool {
        self.no_implicit
    }
    fn num_radical_electrons(&self) -> usize {
        self.num_radical_electrons
    }
    fn degree(&self) -> usize {
        self.degree
    }
    fn explicit_valence(&self) -> Option<usize> {
        self.explicit_valence
    }
    fn implicit_hydrogens(&self) -> Option<usize> {
        self.implicit_hydrogens
    }
    fn total_num_hs(&self) -> Option<usize> {
        self.total_num_hs
    }
    fn total_valence(&self) -> Option<usize> {
        self.total_valence
    }

    fn __repr__(&self) -> String {
        format!(
            "Atom(idx={}, atomic_num={}, formal_charge={}, chiral_tag='{}', isotope={}, is_aromatic={}, degree={})",
            self.idx,
            self.atomic_num,
            self.formal_charge,
            self.chiral_tag,
            self.isotope
                .map(|x| x.to_string())
                .unwrap_or_else(|| "None".to_string()),
            self.is_aromatic,
            self.degree
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Bond {
    fn idx(&self) -> usize {
        self.idx
    }
    fn begin_atom_idx(&self) -> usize {
        self.begin_atom_idx
    }
    fn end_atom_idx(&self) -> usize {
        self.end_atom_idx
    }
    fn bond_type(&self) -> String {
        self.bond_type.clone()
    }
    fn bond_dir(&self) -> String {
        self.bond_dir.clone()
    }
    fn stereo(&self) -> String {
        self.stereo.clone()
    }
    fn stereo_atoms(&self) -> Vec<usize> {
        self.stereo_atoms.clone()
    }
    fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }

    fn __repr__(&self) -> String {
        format!(
            "Bond(idx={}, begin_atom_idx={}, end_atom_idx={}, bond_type='{}', bond_dir='{}', stereo='{}')",
            self.idx,
            self.begin_atom_idx,
            self.end_atom_idx,
            self.bond_type,
            self.bond_dir,
            self.stereo
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
struct MoleculeEdit {
    working: cosmolkit_core::Molecule,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MoleculeEdit {
    fn add_atom(&mut self, element: &str) -> PyResult<usize> {
        let Some(atomic_num) = atomic_number_from_element(element) else {
            return Err(PyValueError::new_err(format!(
                "unsupported element symbol '{element}'"
            )));
        };
        let idx = self.working.add_atom(cosmolkit_core::Atom {
            index: 0,
            atomic_num,
            is_aromatic: false,
            formal_charge: 0,
            explicit_hydrogens: 0,
            no_implicit: false,
            num_radical_electrons: 0,
            chiral_tag: cosmolkit_core::ChiralTag::Unspecified,
            isotope: None,
            atom_map_num: None,
        });
        Ok(idx)
    }

    #[pyo3(signature = (begin, end, order))]
    fn add_bond(&mut self, begin: usize, end: usize, order: &str) -> PyResult<()> {
        if begin >= self.working.atoms.len() || end >= self.working.atoms.len() {
            return Err(PyValueError::new_err("bond atom index out of range"));
        }
        let order = match order.to_ascii_lowercase().as_str() {
            "single" => cosmolkit_core::BondOrder::Single,
            "double" => cosmolkit_core::BondOrder::Double,
            "triple" => cosmolkit_core::BondOrder::Triple,
            "aromatic" => cosmolkit_core::BondOrder::Aromatic,
            "dative" => cosmolkit_core::BondOrder::Dative,
            "unspecified" | "null" => cosmolkit_core::BondOrder::Null,
            _ => {
                return Err(PyValueError::new_err(format!(
                    "unsupported bond order '{order}'"
                )));
            }
        };
        self.working.add_bond(cosmolkit_core::Bond {
            index: 0,
            begin_atom: begin,
            end_atom: end,
            order,
            is_aromatic: matches!(order, cosmolkit_core::BondOrder::Aromatic),
            direction: cosmolkit_core::BondDirection::None,
            stereo: cosmolkit_core::BondStereo::None,
            stereo_atoms: Vec::new(),
        });
        Ok(())
    }

    fn set_atom_charge(&mut self, atom_index: usize, charge: i32) -> PyResult<()> {
        let atom = self
            .working
            .atoms
            .get_mut(atom_index)
            .ok_or_else(|| PyValueError::new_err("atom index out of range"))?;
        let charge =
            i8::try_from(charge).map_err(|_| PyValueError::new_err("charge out of i8 range"))?;
        atom.formal_charge = charge;
        Ok(())
    }

    #[pyo3(signature = (sanitize=None))]
    fn commit(&mut self, sanitize: Option<bool>) -> PyResult<Molecule> {
        let _ = sanitize;
        Ok(Molecule {
            inner: self.working.clone(),
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "MoleculeEdit(num_atoms={}, num_bonds={})",
            self.working.atoms.len(),
            self.working.bonds.len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
struct QueryMolecule;

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
struct Fingerprint;

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Fingerprint {
    fn tanimoto(&self, other: &Fingerprint) -> PyResult<f64> {
        let _ = other;
        Err(unimplemented_api("Fingerprint.tanimoto"))
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
struct SubstructureMatch;

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
struct ValenceReport;

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
struct Alignment;

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
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

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass]
struct AlignmentResult;

#[pymodule]
fn cosmolkit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(placeholder, m)?)?;
    m.add_function(wrap_pyfunction!(rust_version, m)?)?;
    m.add_function(wrap_pyfunction!(core_version, m)?)?;
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

#[cfg(feature = "stubgen")]
define_stub_info_gatherer!(stub_info);
