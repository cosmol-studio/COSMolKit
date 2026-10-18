use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "fixed_projection._native")]
pub struct AvalonFingerprintParams {
    pub(crate) inner: ck::AvalonFingerprintParams,
}
#[pymethods]
impl AvalonFingerprintParams {
    #[new]
    #[pyo3(signature = (*, n_bits=ck::AvalonFingerprintParams::default().n_bits, is_query=ck::AvalonFingerprintParams::default().is_query, bit_flags=ck::AvalonFingerprintParams::default().bit_flags.bits()))]
    fn new(n_bits: u32, is_query: bool, bit_flags: u32) -> Self {
        Self {
            inner: ck::AvalonFingerprintParams {
                n_bits,
                is_query,
                bit_flags: ck::AvalonFingerprintFlags::from_bits_retain(bit_flags),
            },
        }
    }
    #[getter]
    fn n_bits(&self) -> u32 {
        self.inner.n_bits
    }
    #[setter]
    fn set_n_bits(&mut self, value: u32) {
        self.inner.n_bits = value;
    }
    #[getter]
    fn is_query(&self) -> bool {
        self.inner.is_query
    }
    #[setter]
    fn set_is_query(&mut self, value: bool) {
        self.inner.is_query = value;
    }
    #[getter]
    fn bit_flags(&self) -> u32 {
        self.inner.bit_flags.bits()
    }
    #[setter]
    fn set_bit_flags(&mut self, value: u32) {
        self.inner.bit_flags = ck::AvalonFingerprintFlags::from_bits_retain(value);
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct CipLabelOptions {
    pub(crate) inner: ck::CipLabelOptions,
}
impl CipLabelOptions {
    fn make(
        atoms: Option<Vec<usize>>,
        bonds: Option<Vec<usize>>,
        max_recursive_iterations: u32,
    ) -> ck::CipLabelOptions {
        let mut inner =
            ck::CipLabelOptions::default().with_max_recursive_iterations(max_recursive_iterations);
        if let Some(atoms) = atoms {
            inner = inner.with_atoms(atoms.into_iter().map(ck::AtomId::new));
        }
        if let Some(bonds) = bonds {
            inner = inner.with_bonds(bonds.into_iter().map(ck::BondId::new));
        }
        inner
    }
}
#[pymethods]
impl CipLabelOptions {
    #[new]
    #[pyo3(signature = (*, atoms=None, bonds=None, max_recursive_iterations=ck::CipLabelOptions::default().max_recursive_iterations()))]
    fn new(
        atoms: Option<Vec<usize>>,
        bonds: Option<Vec<usize>>,
        max_recursive_iterations: u32,
    ) -> Self {
        Self {
            inner: Self::make(atoms, bonds, max_recursive_iterations),
        }
    }
    #[getter]
    fn atoms(&self) -> Option<Vec<usize>> {
        self.inner
            .atoms()
            .map(|atoms| atoms.iter().map(|id| id.index()).collect())
    }
    #[setter]
    fn set_atoms(&mut self, value: Option<Vec<usize>>) {
        self.inner = Self::make(value, self.bonds(), self.inner.max_recursive_iterations());
    }
    #[getter]
    fn bonds(&self) -> Option<Vec<usize>> {
        self.inner
            .bonds()
            .map(|bonds| bonds.iter().map(|id| id.index()).collect())
    }
    #[setter]
    fn set_bonds(&mut self, value: Option<Vec<usize>>) {
        self.inner = Self::make(self.atoms(), value, self.inner.max_recursive_iterations());
    }
    #[getter]
    fn max_recursive_iterations(&self) -> u32 {
        self.inner.max_recursive_iterations()
    }
    #[setter]
    fn set_max_recursive_iterations(&mut self, value: u32) {
        self.inner = self.inner.clone().with_max_recursive_iterations(value);
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct LegacyTopologicalTorsionParams {
    pub(crate) inner: ck::LegacyTopologicalTorsionParams,
}
#[pymethods]
impl LegacyTopologicalTorsionParams {
    #[new]
    #[pyo3(signature = (*, torsion_atom_count=ck::LegacyTopologicalTorsionParams::default().torsion_atom_count, include_chirality=ck::LegacyTopologicalTorsionParams::default().include_chirality, fp_size=ck::LegacyTopologicalTorsionParams::default().fp_size, bits_per_entry=ck::LegacyTopologicalTorsionParams::default().bits_per_entry, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None))]
    fn construct(
        torsion_atom_count: u32,
        include_chirality: bool,
        fp_size: u32,
        bits_per_entry: u32,
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
    ) -> Self {
        Self {
            inner: ck::LegacyTopologicalTorsionParams::new(
                torsion_atom_count,
                include_chirality,
                fp_size,
                bits_per_entry,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
            ),
        }
    }
    #[staticmethod]
    #[pyo3(signature = (torsion_atom_count=ck::LegacyTopologicalTorsionParams::default().torsion_atom_count, include_chirality=ck::LegacyTopologicalTorsionParams::default().include_chirality, fp_size=ck::LegacyTopologicalTorsionParams::default().fp_size, bits_per_entry=ck::LegacyTopologicalTorsionParams::default().bits_per_entry, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None))]
    fn new(
        torsion_atom_count: u32,
        include_chirality: bool,
        fp_size: u32,
        bits_per_entry: u32,
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
    ) -> Self {
        Self::construct(
            torsion_atom_count,
            include_chirality,
            fp_size,
            bits_per_entry,
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
        )
    }
    #[getter]
    fn torsion_atom_count(&self) -> u32 {
        self.inner.torsion_atom_count
    }
    #[setter]
    fn set_torsion_atom_count(&mut self, value: u32) {
        self.inner.torsion_atom_count = value;
    }
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    #[setter]
    fn set_include_chirality(&mut self, value: bool) {
        self.inner.include_chirality = value;
    }
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    #[setter]
    fn set_fp_size(&mut self, value: u32) {
        self.inner.fp_size = value;
    }
    #[getter]
    fn bits_per_entry(&self) -> u32 {
        self.inner.bits_per_entry
    }
    #[setter]
    fn set_bits_per_entry(&mut self, value: u32) {
        self.inner.bits_per_entry = value;
    }
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    #[setter]
    fn set_from_atoms(&mut self, value: Option<Vec<u32>>) {
        self.inner.from_atoms = value;
    }
    #[getter]
    fn ignore_atoms(&self) -> Option<Vec<u32>> {
        self.inner.ignore_atoms.clone()
    }
    #[setter]
    fn set_ignore_atoms(&mut self, value: Option<Vec<u32>>) {
        self.inner.ignore_atoms = value;
    }
    #[getter]
    fn custom_atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_atom_invariants.clone()
    }
    #[setter]
    fn set_custom_atom_invariants(&mut self, value: Option<Vec<u32>>) {
        self.inner.custom_atom_invariants = value;
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct ConformerAlignmentParameters {
    pub(crate) inner: ck::ConformerAlignmentParameters,
}
#[pymethods]
impl ConformerAlignmentParameters {
    #[new]
    #[pyo3(signature = (*, atom_indices=None, conformer_ids=None, weights=None, reflect=ck::ConformerAlignmentParameters::default().reflect, max_iterations=ck::ConformerAlignmentParameters::default().max_iterations))]
    fn construct(
        atom_indices: Option<Vec<usize>>,
        conformer_ids: Option<Vec<usize>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> Self {
        Self {
            inner: ck::ConformerAlignmentParameters::new(
                atom_indices,
                conformer_ids,
                weights,
                reflect,
                max_iterations,
            ),
        }
    }
    #[staticmethod]
    #[pyo3(signature = (atom_indices=None, conformer_ids=None, weights=None, reflect=ck::ConformerAlignmentParameters::default().reflect, max_iterations=ck::ConformerAlignmentParameters::default().max_iterations))]
    fn new(
        atom_indices: Option<Vec<usize>>,
        conformer_ids: Option<Vec<usize>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> Self {
        Self::construct(
            atom_indices,
            conformer_ids,
            weights,
            reflect,
            max_iterations,
        )
    }
    #[getter]
    fn atom_indices(&self) -> Option<Vec<usize>> {
        self.inner.atom_indices.clone()
    }
    #[setter]
    fn set_atom_indices(&mut self, value: Option<Vec<usize>>) {
        self.inner.atom_indices = value;
    }
    #[getter]
    fn conformer_ids(&self) -> Option<Vec<usize>> {
        self.inner.conformer_ids.clone()
    }
    #[setter]
    fn set_conformer_ids(&mut self, value: Option<Vec<usize>>) {
        self.inner.conformer_ids = value;
    }
    #[getter]
    fn weights(&self) -> Option<Vec<f64>> {
        self.inner.weights.clone()
    }
    #[setter]
    fn set_weights(&mut self, value: Option<Vec<f64>>) {
        self.inner.weights = value;
    }
    #[getter]
    fn reflect(&self) -> bool {
        self.inner.reflect
    }
    #[setter]
    fn set_reflect(&mut self, value: bool) {
        self.inner.reflect = value;
    }
    #[getter]
    fn max_iterations(&self) -> u32 {
        self.inner.max_iterations
    }
    #[setter]
    fn set_max_iterations(&mut self, value: u32) {
        self.inner.max_iterations = value;
    }
}
