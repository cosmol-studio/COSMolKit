//! Typed RDKit object transport; chemistry remains canonical Rust operations.
use cosmolkit as ck;
use pyo3::{
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
};

/// A private checked handle, not a public dynamically typed binding parameter.
pub struct RdkitMol {
    handle: Py<PyAny>,
}
impl FromPyObject<'_, '_> for RdkitMol {
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr =
        pyo3::type_hint_identifier!("rdkit.Chem.rdchem", "Mol");
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self> {
        let mol = value.py().import("rdkit.Chem.rdchem")?.getattr("Mol")?;
        if !value.is_instance(&mol)? {
            return Err(PyTypeError::new_err("expected rdkit.Chem.rdchem.Mol"));
        }
        Ok(Self {
            handle: value.to_owned().unbind(),
        })
    }
}

fn method<'py>(object: &Bound<'py, PyAny>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    object.call_method0(name)
}
fn indexed<'py>(
    object: &Bound<'py, PyAny>,
    name: &str,
    index: usize,
) -> PyResult<Bound<'py, PyAny>> {
    object.call_method1(name, (index,))
}
fn extract<T>(object: &Bound<'_, PyAny>, name: &str) -> PyResult<T>
where
    for<'a> T: FromPyObject<'a, 'a>,
{
    method(object, name)?.extract::<T>().map_err(Into::into)
}
fn enum_value<T>(
    object: &Bound<'_, PyAny>,
    name: &str,
    convert: impl FnOnce(i64) -> Option<T>,
) -> PyResult<T> {
    let code = extract(object, name)?;
    convert(code)
        .ok_or_else(|| PyValueError::new_err(format!("from_rdkit unsupported {name} code {code}")))
}
fn coordinate(object: &Bound<'_, PyAny>, axis: &str) -> PyResult<f64> {
    object.getattr(axis)?.extract()
}

#[pyfunction]
#[pyo3(signature=(rdmol, sanitize=None))]
pub fn mol_from_rdkit(
    py: Python<'_>,
    rdmol: RdkitMol,
    sanitize: Option<bool>,
) -> PyResult<crate::native::Molecule> {
    import(py, &rdmol, sanitize)
}

/// Copy the 0.3.0 graph fields and 3D conformers from an RDKit molecule.
///
/// `sanitize=None` prepares valence only, `True` sanitizes, and `False`
/// leaves caches unprepared. 2D conformers are skipped, as in 0.3.0.
/// This is not an RDKit pickle/property/SGroup import.
pub(crate) fn import(
    py: Python<'_>,
    source: &RdkitMol,
    sanitize: Option<bool>,
) -> PyResult<crate::native::Molecule> {
    let rdmol = source.handle.bind(py);
    // COSMolKit 0.3.0 python/src/lib.rs:4971-4973 (tag v0.3.0):
    // let atom_count: usize = py_method_extract(rdmol, "GetNumAtoms")?;
    // let bond_count: usize = py_method_extract(rdmol, "GetNumBonds")?;
    // let mut builder = cosmolkit_core::MoleculeBuilder::new();
    // COSMolKit✔️✔️: Transport uses the current public checked builder.
    // Enum codes use the single canonical RDKit mapping, not a second table.
    // Collect detached rows once and validate the complete graph before build.
    // O(V+E) graph transport avoids the current builder's per-append topology
    // copies/revalidation without bypassing any final construction checks.
    let atom_count = extract::<usize>(rdmol, "GetNumAtoms")?;
    let bond_count = extract::<usize>(rdmol, "GetNumBonds")?;
    let mut atoms = Vec::with_capacity(atom_count);
    let mut bonds = Vec::with_capacity(bond_count);
    for index in 0..atom_count {
        let atom = indexed(rdmol, "GetAtomWithIdx", index)?;
        let number = extract::<u8>(&atom, "GetAtomicNum")?;
        let element = ck::Element::from_atomic_number(number).ok_or_else(|| {
            PyValueError::new_err(format!(
                "from_rdkit atom {index} unsupported atomic number {number}"
            ))
        })?;
        let mut spec = ck::AtomSpec::new(element)
            .with_formal_charge(extract(&atom, "GetFormalCharge")?)
            .with_explicit_hydrogens(extract(&atom, "GetNumExplicitHs")?)
            .with_no_implicit(extract(&atom, "GetNoImplicit")?)
            .with_chiral_tag(enum_value(
                &atom,
                "GetChiralTag",
                ck::ChiralTag::from_rdkit_code,
            )?)
            .with_aromatic(extract(&atom, "GetIsAromatic")?)
            .with_radical_electrons(extract(&atom, "GetNumRadicalElectrons")?)
            .with_hybridization(enum_value(
                &atom,
                "GetHybridization",
                ck::Hybridization::from_rdkit_code,
            )?);
        let isotope = extract::<u16>(&atom, "GetIsotope")?;
        let atom_map = extract::<u32>(&atom, "GetAtomMapNum")?;
        if isotope != 0 {
            spec = spec.with_isotope(isotope);
        }
        if atom_map != 0 {
            spec = spec.with_atom_map(atom_map);
        }
        atoms.push(ck::Atom::from_spec(ck::AtomId::new(index), spec));
    }
    for index in 0..bond_count {
        let bond = indexed(rdmol, "GetBondWithIdx", index)?;
        let begin = extract::<usize>(&bond, "GetBeginAtomIdx")?;
        let end = extract::<usize>(&bond, "GetEndAtomIdx")?;
        if begin >= atom_count || end >= atom_count {
            return Err(PyValueError::new_err(format!(
                "from_rdkit bond {index} atom index out of range: {begin}-{end}"
            )));
        }
        let mut spec = ck::BondSpec::new(
            ck::AtomId::new(begin),
            ck::AtomId::new(end),
            enum_value(&bond, "GetBondType", ck::BondOrder::from_rdkit_code)?,
        )
        .with_aromatic(extract(&bond, "GetIsAromatic")?)
        .with_direction(enum_value(
            &bond,
            "GetBondDir",
            ck::BondDirection::from_rdkit_code,
        )?)
        .with_stereo(enum_value(
            &bond,
            "GetStereo",
            ck::BondStereo::from_rdkit_code,
        )?);
        let stereo_atoms = extract::<Vec<usize>>(&bond, "GetStereoAtoms")?;
        match stereo_atoms.as_slice() {
            [] => {}
            [first, second] => {
                spec = spec.with_stereo_atoms(ck::AtomId::new(*first), ck::AtomId::new(*second));
            }
            _ => {
                return Err(PyValueError::new_err(format!(
                    "from_rdkit bond {index} requires zero or two stereo atoms"
                )));
            }
        }
        bonds.push(ck::Bond::from_spec(ck::BondId::new(index), spec));
    }
    let topology = ck::TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .map_err(|error| {
            crate::operation_errors::OperationError::from_rust(ck::OperationError::InvalidTopology(
                error,
            ))
        })?;
    let mut builder = ck::MoleculeBuilder::from_parts(
        topology,
        ck::CoordinateBlock::default(),
        ck::MoleculeProperties::default(),
    );
    // COSMolKit 0.3.0 python/src/lib.rs:5072-5074:
    // if !py_method_extract::<bool>(&conformer, "Is3D")? {
    //     continue;
    // }
    // COSMolKit✔️✔️: Keep dimension filtering and iteration order. The
    // builder validates row count and finiteness; no geometry is generated.
    for conformer in method(rdmol, "GetConformers")?.try_iter()? {
        let conformer = conformer?;
        if !extract::<bool>(&conformer, "Is3D")? {
            continue;
        }
        let mut positions = Vec::with_capacity(atom_count);
        for atom in 0..atom_count {
            let point = indexed(&conformer, "GetAtomPosition", atom)?;
            positions.push([
                coordinate(&point, "x")?,
                coordinate(&point, "y")?,
                coordinate(&point, "z")?,
            ]);
        }
        builder
            .add_3d_conformer(positions)
            .map_err(|error| crate::operation_errors::OperationError::from_rust(error))?;
    }
    let molecule = builder
        .build()
        .map_err(|error| crate::operation_errors::OperationError::from_rust(error))?;
    // COSMolKit 0.3.0 python/src/lib.rs:5091-5099:
    // let inner = match sanitize {
    //     Some(true) => mol
    //         .sanitize()
    //         .map_err(|err| PyValueError::new_err(err.to_string()))?,
    //     Some(false) => mol,
    //     None => mol
    //         .with_assigned_valence()
    //         .map_err(|err| PyValueError::new_err(err.to_string()))?,
    // };
    // COSMolKit✔️✔️: Chemistry delegates to the same registered operations.
    let inner = match sanitize {
        Some(true) => molecule.sanitize(),
        Some(false) => Ok(molecule),
        None => molecule.with_assigned_valence(),
    }
    .map_err(|error| crate::operation_errors::OperationError::from_rust(error))?;
    Ok(crate::native::Molecule { inner })
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::types::PyDict;
    use std::ffi::CString;

    #[test]
    fn rdkit_import_preserves_fields_dimensions_native_sanitize_and_original_external_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let chem = py.import("rdkit.Chem")?;
            for smiles in [
                "[13CH3:7][C@H](F)Cl",
                "F/C=C/F",
                "c1ccccc1",
                "[NH3+][O-]",
                "[CH2]",
            ] {
                let source = chem.call_method1("MolFromSmiles", (smiles,))?;
                assert!(!source.is_none());
                let before = source.call_method0("ToBinary")?.extract::<Vec<u8>>()?;
                let typed = source.extract::<RdkitMol>()?;
                let raw = import(py, &typed, Some(false))?;
                assert_eq!(
                    raw.inner.num_atoms(),
                    extract::<usize>(&source, "GetNumAtoms")?
                );
                assert_eq!(
                    raw.inner.num_bonds(),
                    extract::<usize>(&source, "GetNumBonds")?
                );
                for (index, atom) in raw.inner.atoms().iter().enumerate() {
                    let rd = indexed(&source, "GetAtomWithIdx", index)?;
                    assert_eq!(atom.atomic_number(), extract::<u8>(&rd, "GetAtomicNum")?);
                    assert_eq!(atom.formal_charge(), extract::<i8>(&rd, "GetFormalCharge")?);
                    assert_eq!(
                        atom.explicit_hydrogens(),
                        extract::<u8>(&rd, "GetNumExplicitHs")?
                    );
                    assert_eq!(atom.no_implicit(), extract::<bool>(&rd, "GetNoImplicit")?);
                    assert_eq!(
                        atom.chiral_tag(),
                        enum_value(&rd, "GetChiralTag", ck::ChiralTag::from_rdkit_code)?
                    );
                    assert_eq!(atom.is_aromatic(), extract::<bool>(&rd, "GetIsAromatic")?);
                    assert_eq!(
                        atom.radical_electrons(),
                        extract::<u8>(&rd, "GetNumRadicalElectrons")?
                    );
                    assert_eq!(
                        atom.hybridization(),
                        enum_value(&rd, "GetHybridization", ck::Hybridization::from_rdkit_code)?
                    );
                    let isotope = extract::<u16>(&rd, "GetIsotope")?;
                    let map = extract::<u32>(&rd, "GetAtomMapNum")?;
                    assert_eq!(atom.isotope(), (isotope != 0).then_some(isotope));
                    assert_eq!(atom.atom_map(), (map != 0).then_some(map));
                }
                for (index, bond) in raw.inner.bonds().iter().enumerate() {
                    let rd = indexed(&source, "GetBondWithIdx", index)?;
                    assert_eq!(
                        bond.begin().index(),
                        extract::<usize>(&rd, "GetBeginAtomIdx")?
                    );
                    assert_eq!(bond.end().index(), extract::<usize>(&rd, "GetEndAtomIdx")?);
                    assert_eq!(
                        bond.order(),
                        enum_value(&rd, "GetBondType", ck::BondOrder::from_rdkit_code)?
                    );
                    assert_eq!(
                        bond.direction(),
                        enum_value(&rd, "GetBondDir", ck::BondDirection::from_rdkit_code)?
                    );
                    assert_eq!(
                        bond.stereo(),
                        enum_value(&rd, "GetStereo", ck::BondStereo::from_rdkit_code)?
                    );
                    assert_eq!(bond.is_aromatic(), extract::<bool>(&rd, "GetIsAromatic")?);
                    let stereo = extract::<Vec<usize>>(&rd, "GetStereoAtoms")?;
                    assert_eq!(
                        bond.stereo_atoms()
                            .map(|ids| ids.map(|id| id.index()).to_vec())
                            .unwrap_or_default(),
                        stereo
                    );
                }
                let default = import(py, &typed, None)?;
                let sanitized = import(py, &typed, Some(true))?;
                let public = py
                    .get_type::<crate::native::Molecule>()
                    .call_method1("from_rdkit", (&source,))?
                    .extract::<Py<crate::native::Molecule>>()?;
                assert_eq!(
                    public.borrow(py).inner.to_binary().unwrap(),
                    default.inner.to_binary().unwrap()
                );
                assert_eq!(
                    default.inner.to_binary().unwrap(),
                    raw.inner
                        .with_assigned_valence()
                        .unwrap()
                        .to_binary()
                        .unwrap()
                );
                assert_eq!(
                    sanitized.inner.to_binary().unwrap(),
                    raw.inner.sanitize().unwrap().to_binary().unwrap()
                );
                assert!(raw.inner.atom_metadata(false).is_err());
                assert!(default.inner.atom_metadata(false).is_ok());
                assert_eq!(
                    source.call_method0("ToBinary")?.extract::<Vec<u8>>()?,
                    before
                );
                let free = mol_from_rdkit(py, source.extract::<RdkitMol>()?, None)?;
                assert_eq!(
                    free.inner.to_binary().unwrap(),
                    default.inner.to_binary().unwrap()
                );
            }
            let scope = PyDict::new(py);
            scope.set_item("Chem", &chem)?;
            py.run(
                &CString::new(
                    r#"
source = Chem.MolFromSmiles("[13CH3:7][C@H](F)Cl")
positions = [
    [[10., 11., 0.], [12., 13., 0.], [14., 15., 0.], [16., 17., 0.]],
    [[-0., .1, -.2], [1., 2., 3.], [4., 5., 6.], [7., 8., 9.]],
    [[.25, .5, .75], [1.25, 1.5, 1.75], [2.25, 2.5, 2.75], [3.25, 3.5, 3.75]],
]
for index, rows in enumerate(positions):
    conf = Chem.Conformer(4)
    conf.Set3D(index != 0)
    for atom, point in enumerate(rows):
        conf.SetAtomPosition(atom, point)
    source.AddConformer(conf, assignId=True)
marker = LookupError("actual RDKit subclass failure")
class Broken(Chem.Mol):
    def GetNumAtoms(self):
        raise marker
broken = Broken(source)
"#,
                )
                .unwrap(),
                Some(&scope),
                Some(&scope),
            )?;
            let source = scope.get_item("source")?.unwrap();
            let before = source.call_method0("ToBinary")?.extract::<Vec<u8>>()?;
            let output = import(py, &source.extract::<RdkitMol>()?, None)?;
            assert_eq!(output.inner.num_3d_conformers(), 2);
            assert!(!output.inner.has_2d_coordinates());
            let expected = scope
                .get_item("positions")?
                .unwrap()
                .extract::<Vec<Vec<[f64; 3]>>>()?;
            for index in 0..2 {
                let actual = output.inner.coordinates_3d(index).unwrap();
                for (a, b) in actual
                    .iter()
                    .flatten()
                    .zip(expected[index + 1].iter().flatten())
                {
                    assert_eq!(a.to_bits(), b.to_bits());
                }
            }
            assert_eq!(
                source.call_method0("ToBinary")?.extract::<Vec<u8>>()?,
                before
            );
            source
                .call_method1("GetAtomWithIdx", (0,))?
                .call_method1("SetIsotope", (12,))?;
            source
                .call_method0("GetConformer")?
                .call_method1("SetAtomPosition", (0, (100., 200., 300.)))?;
            assert_eq!(output.inner.atoms()[0].isotope(), Some(13));
            assert_eq!(
                output.inner.coordinates_3d(0).unwrap()[0][0].to_bits(),
                (-0.0f64).to_bits()
            );
            let broken = scope.get_item("broken")?.unwrap().extract::<RdkitMol>()?;
            let error = import(py, &broken, None)
                .err()
                .expect("actual original external failure");
            assert!(error.value(py).is(scope.get_item("marker")?.unwrap()));
            assert!(
                py.None()
                    .bind(py)
                    .extract::<RdkitMol>()
                    .err()
                    .expect("checked external Mol input rejects None")
                    .is_instance_of::<PyTypeError>(py)
            );
            let conf = source.call_method1("GetConformer", (1,))?;
            conf.call_method1("SetAtomPosition", (0, (f64::NAN, 0., 0.)))?;
            let before = source.call_method0("ToBinary")?.extract::<Vec<u8>>()?;
            let error = import(py, &source.extract::<RdkitMol>()?, None)
                .err()
                .expect("native coordinate validation");
            assert!(error.is_instance(py, <crate::operation_errors::OperationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(
                source.call_method0("ToBinary")?.extract::<Vec<u8>>()?,
                before
            );
            Ok(())
        })
        .unwrap();
    }
}
