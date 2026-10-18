use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "fixed_projection._native")]
pub struct UffForceFieldParams {
    pub(crate) inner: ck::UffForceFieldParams,
}

#[pymethods]
impl UffForceFieldParams {
    #[new]
    #[pyo3(signature = (*, conformer_id=ck::UffForceFieldParams::default().conformer_id(), vdw_threshold=ck::UffForceFieldParams::default().vdw_threshold(), ignore_interfragment_interactions=ck::UffForceFieldParams::default().ignore_interfragment_interactions()))]
    fn new(
        conformer_id: Option<usize>,
        vdw_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        Self {
            inner: ck::UffForceFieldParams::new(
                conformer_id,
                vdw_threshold,
                ignore_interfragment_interactions,
            ),
        }
    }
    #[getter]
    fn conformer_id(&self) -> Option<usize> {
        self.inner.conformer_id()
    }
    #[setter]
    fn set_conformer_id(&mut self, value: Option<usize>) {
        self.inner = ck::UffForceFieldParams::new(
            value,
            self.inner.vdw_threshold(),
            self.inner.ignore_interfragment_interactions(),
        );
    }
    #[getter]
    fn vdw_threshold(&self) -> f64 {
        self.inner.vdw_threshold()
    }
    #[setter]
    fn set_vdw_threshold(&mut self, value: f64) {
        self.inner = ck::UffForceFieldParams::new(
            self.inner.conformer_id(),
            value,
            self.inner.ignore_interfragment_interactions(),
        );
    }
    #[getter]
    fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions()
    }
    #[setter]
    fn set_ignore_interfragment_interactions(&mut self, value: bool) {
        self.inner = ck::UffForceFieldParams::new(
            self.inner.conformer_id(),
            self.inner.vdw_threshold(),
            value,
        );
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct MmffForceFieldParams {
    pub(crate) inner: ck::MmffForceFieldParams,
}

#[pymethods]
impl MmffForceFieldParams {
    #[new]
    #[pyo3(signature = (*, conformer_id=ck::MmffForceFieldParams::default().conformer_id(), mmff_variant=ck::MmffForceFieldParams::default().mmff_variant().to_owned(), non_bonded_threshold=ck::MmffForceFieldParams::default().non_bonded_threshold(), ignore_interfragment_interactions=ck::MmffForceFieldParams::default().ignore_interfragment_interactions()))]
    fn new(
        conformer_id: Option<usize>,
        mmff_variant: String,
        non_bonded_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        Self {
            inner: ck::MmffForceFieldParams::new(
                conformer_id,
                mmff_variant,
                non_bonded_threshold,
                ignore_interfragment_interactions,
            ),
        }
    }
    #[getter]
    fn conformer_id(&self) -> Option<usize> {
        self.inner.conformer_id()
    }
    #[setter]
    fn set_conformer_id(&mut self, value: Option<usize>) {
        self.inner = ck::MmffForceFieldParams::new(
            value,
            self.inner.mmff_variant().to_owned(),
            self.inner.non_bonded_threshold(),
            self.inner.ignore_interfragment_interactions(),
        );
    }
    #[getter]
    fn mmff_variant(&self) -> &str {
        self.inner.mmff_variant()
    }
    #[setter]
    fn set_mmff_variant(&mut self, value: String) {
        self.inner = ck::MmffForceFieldParams::new(
            self.inner.conformer_id(),
            value,
            self.inner.non_bonded_threshold(),
            self.inner.ignore_interfragment_interactions(),
        );
    }
    #[getter]
    fn non_bonded_threshold(&self) -> f64 {
        self.inner.non_bonded_threshold()
    }
    #[setter]
    fn set_non_bonded_threshold(&mut self, value: f64) {
        self.inner = ck::MmffForceFieldParams::new(
            self.inner.conformer_id(),
            self.inner.mmff_variant().to_owned(),
            value,
            self.inner.ignore_interfragment_interactions(),
        );
    }
    #[getter]
    fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions()
    }
    #[setter]
    fn set_ignore_interfragment_interactions(&mut self, value: bool) {
        self.inner = ck::MmffForceFieldParams::new(
            self.inner.conformer_id(),
            self.inner.mmff_variant().to_owned(),
            self.inner.non_bonded_threshold(),
            value,
        );
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct SmilesWriteParams {
    pub(crate) inner: ck::SmilesWriteParams,
}

#[pymethods]
impl SmilesWriteParams {
    #[new]
    #[pyo3(signature = (*, isomeric_smiles=ck::SmilesWriteParams::default().isomeric_smiles, kekule=ck::SmilesWriteParams::default().kekule, canonical=ck::SmilesWriteParams::default().canonical, clean_stereo=ck::SmilesWriteParams::default().clean_stereo, rooted_at_atom=ck::SmilesWriteParams::default().rooted_at_atom.map(ck::AtomId::index), all_bonds_explicit=ck::SmilesWriteParams::default().all_bonds_explicit, all_hydrogens_explicit=ck::SmilesWriteParams::default().all_hydrogens_explicit, include_dative_bonds=ck::SmilesWriteParams::default().include_dative_bonds, ignore_atom_map_numbers=ck::SmilesWriteParams::default().ignore_atom_map_numbers))]
    fn new(
        isomeric_smiles: bool,
        kekule: bool,
        canonical: bool,
        clean_stereo: bool,
        rooted_at_atom: Option<usize>,
        all_bonds_explicit: bool,
        all_hydrogens_explicit: bool,
        include_dative_bonds: bool,
        ignore_atom_map_numbers: bool,
    ) -> Self {
        Self {
            inner: ck::SmilesWriteParams {
                isomeric_smiles,
                kekule,
                canonical,
                clean_stereo,
                rooted_at_atom: rooted_at_atom.map(ck::AtomId::new),
                all_bonds_explicit,
                all_hydrogens_explicit,
                include_dative_bonds,
                ignore_atom_map_numbers,
            },
        }
    }
    #[getter]
    fn isomeric_smiles(&self) -> bool {
        self.inner.isomeric_smiles
    }
    #[setter]
    fn set_isomeric_smiles(&mut self, value: bool) {
        self.inner.isomeric_smiles = value;
    }
    #[getter]
    fn kekule(&self) -> bool {
        self.inner.kekule
    }
    #[setter]
    fn set_kekule(&mut self, value: bool) {
        self.inner.kekule = value;
    }
    #[getter]
    fn canonical(&self) -> bool {
        self.inner.canonical
    }
    #[setter]
    fn set_canonical(&mut self, value: bool) {
        self.inner.canonical = value;
    }
    #[getter]
    fn clean_stereo(&self) -> bool {
        self.inner.clean_stereo
    }
    #[setter]
    fn set_clean_stereo(&mut self, value: bool) {
        self.inner.clean_stereo = value;
    }
    #[getter]
    fn rooted_at_atom(&self) -> Option<usize> {
        self.inner.rooted_at_atom.map(ck::AtomId::index)
    }
    #[setter]
    fn set_rooted_at_atom(&mut self, value: Option<usize>) {
        self.inner.rooted_at_atom = value.map(ck::AtomId::new);
    }
    #[getter]
    fn all_bonds_explicit(&self) -> bool {
        self.inner.all_bonds_explicit
    }
    #[setter]
    fn set_all_bonds_explicit(&mut self, value: bool) {
        self.inner.all_bonds_explicit = value;
    }
    #[getter]
    fn all_hydrogens_explicit(&self) -> bool {
        self.inner.all_hydrogens_explicit
    }
    #[setter]
    fn set_all_hydrogens_explicit(&mut self, value: bool) {
        self.inner.all_hydrogens_explicit = value;
    }
    #[getter]
    fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds
    }
    #[setter]
    fn set_include_dative_bonds(&mut self, value: bool) {
        self.inner.include_dative_bonds = value;
    }
    #[getter]
    fn ignore_atom_map_numbers(&self) -> bool {
        self.inner.ignore_atom_map_numbers
    }
    #[setter]
    fn set_ignore_atom_map_numbers(&mut self, value: bool) {
        self.inner.ignore_atom_map_numbers = value;
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct AddHsParams {
    pub(crate) inner: ck::AddHsParams,
}

#[pymethods]
impl AddHsParams {
    #[new]
    #[pyo3(signature = (*, explicit_only=ck::AddHsParams::default().explicit_only, add_coords=ck::AddHsParams::default().add_coords, add_residue_info=ck::AddHsParams::default().add_residue_info, skip_queries=ck::AddHsParams::default().skip_queries, only_on_atoms=ck::AddHsParams::default().only_on_atoms.map(|atoms| atoms.into_iter().map(ck::AtomId::index).collect::<Vec<_>>())))]
    fn new(
        explicit_only: bool,
        add_coords: bool,
        add_residue_info: bool,
        skip_queries: bool,
        only_on_atoms: Option<Vec<usize>>,
    ) -> Self {
        Self {
            inner: ck::AddHsParams {
                explicit_only,
                add_coords,
                add_residue_info,
                skip_queries,
                only_on_atoms: only_on_atoms
                    .map(|atoms| atoms.into_iter().map(ck::AtomId::new).collect()),
            },
        }
    }
    #[getter]
    fn explicit_only(&self) -> bool {
        self.inner.explicit_only
    }
    #[setter]
    fn set_explicit_only(&mut self, value: bool) {
        self.inner.explicit_only = value;
    }
    #[getter]
    fn add_coords(&self) -> bool {
        self.inner.add_coords
    }
    #[setter]
    fn set_add_coords(&mut self, value: bool) {
        self.inner.add_coords = value;
    }
    #[getter]
    fn add_residue_info(&self) -> bool {
        self.inner.add_residue_info
    }
    #[setter]
    fn set_add_residue_info(&mut self, value: bool) {
        self.inner.add_residue_info = value;
    }
    #[getter]
    fn skip_queries(&self) -> bool {
        self.inner.skip_queries
    }
    #[setter]
    fn set_skip_queries(&mut self, value: bool) {
        self.inner.skip_queries = value;
    }
    #[getter]
    fn only_on_atoms(&self) -> Option<Vec<usize>> {
        self.inner
            .only_on_atoms
            .as_ref()
            .map(|atoms| atoms.iter().map(|id| id.index()).collect())
    }
    #[setter]
    fn set_only_on_atoms(&mut self, value: Option<Vec<usize>>) {
        self.inner.only_on_atoms =
            value.map(|atoms| atoms.into_iter().map(ck::AtomId::new).collect());
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct ForceFieldMinimizeParams {
    pub(crate) inner: ck::ForceFieldMinimizeParams,
}

#[pymethods]
impl ForceFieldMinimizeParams {
    #[new]
    #[pyo3(signature = (*, max_iterations=ck::ForceFieldMinimizeParams::default().max_iterations(), force_tolerance=ck::ForceFieldMinimizeParams::default().force_tolerance(), energy_tolerance=ck::ForceFieldMinimizeParams::default().energy_tolerance()))]
    fn new(max_iterations: u32, force_tolerance: f64, energy_tolerance: f64) -> Self {
        Self {
            inner: ck::ForceFieldMinimizeParams::new(
                max_iterations,
                force_tolerance,
                energy_tolerance,
            ),
        }
    }
    #[getter]
    fn max_iterations(&self) -> u32 {
        self.inner.max_iterations()
    }
    #[setter]
    fn set_max_iterations(&mut self, value: u32) {
        self.inner = ck::ForceFieldMinimizeParams::new(
            value,
            self.inner.force_tolerance(),
            self.inner.energy_tolerance(),
        );
    }
    #[getter]
    fn force_tolerance(&self) -> f64 {
        self.inner.force_tolerance()
    }
    #[setter]
    fn set_force_tolerance(&mut self, value: f64) {
        self.inner = ck::ForceFieldMinimizeParams::new(
            self.inner.max_iterations(),
            value,
            self.inner.energy_tolerance(),
        );
    }
    #[getter]
    fn energy_tolerance(&self) -> f64 {
        self.inner.energy_tolerance()
    }
    #[setter]
    fn set_energy_tolerance(&mut self, value: f64) {
        self.inner = ck::ForceFieldMinimizeParams::new(
            self.inner.max_iterations(),
            self.inner.force_tolerance(),
            value,
        );
    }
}
