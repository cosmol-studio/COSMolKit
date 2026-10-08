//! Python object transport from the 0.3.0 RDKit adapter, not a chemistry owner.
use ::cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};

fn method<'py>(object: &Bound<'py, PyAny>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    object.call_method0(name).map_err(|error| {
        PyValueError::new_err(format!("from_rdkit failed calling {name}: {error}"))
    })
}

fn indexed<'py>(
    object: &Bound<'py, PyAny>,
    name: &str,
    index: usize,
) -> PyResult<Bound<'py, PyAny>> {
    object.call_method1(name, (index,)).map_err(|error| {
        PyValueError::new_err(format!("from_rdkit failed calling {name}: {error}"))
    })
}

fn extract<T>(object: &Bound<'_, PyAny>, name: &str) -> PyResult<T>
where
    for<'a> T: FromPyObject<'a, 'a>,
{
    method(object, name)?.extract::<T>().map_err(|error| {
        let error: PyErr = error.into();
        PyValueError::new_err(format!(
            "from_rdkit failed extracting result from {name}: {error}"
        ))
    })
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
    object
        .getattr(axis)
        .and_then(|value| value.extract())
        .map_err(|error| {
            PyValueError::new_err(format!(
                "from_rdkit failed reading coordinate {axis}: {error}"
            ))
        })
}

/// Copy the 0.3.0 graph fields and 3D conformers from an RDKit molecule.
///
/// `sanitize=None` prepares valence only, `True` sanitizes, and `False`
/// leaves caches unprepared. 2D conformers are skipped, as in 0.3.0.
/// This is not an RDKit pickle/property/SGroup import.
pub(crate) fn from_rdkit(
    rdmol: &Bound<'_, PyAny>,
    sanitize: Option<bool>,
) -> PyResult<crate::drawing_binding::Molecule> {
    let py = rdmol.py();
    py.import("rdkit.Chem").map_err(|error| {
        PyValueError::new_err(format!(
            "Molecule.from_rdkit requires rdkit to be installed and importable: {error}"
        ))
    })?;
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
            crate::drawing_binding::operation_pyerr(py, ck::OperationError::InvalidTopology(error))
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
            .map_err(|error| crate::drawing_binding::operation_pyerr(py, error))?;
    }
    let molecule = builder
        .build()
        .map_err(|error| crate::drawing_binding::operation_pyerr(py, error))?;
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
    .map_err(|error| crate::drawing_binding::operation_pyerr(py, error))?;
    Ok(crate::drawing_binding::Molecule::from_inner(inner))
}
