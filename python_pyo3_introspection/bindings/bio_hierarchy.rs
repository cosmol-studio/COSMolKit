//! Native owned hierarchy snapshots; no relationship resolution in Python.
use cosmolkit as ck;
use pyo3::prelude::*;

macro_rules! view {
    ($name:ident, $row:ident, $module:ident, {$($methods:tt)*}) => {
        #[pyclass(module = "cosmolkit", frozen, from_py_object)]
        #[derive(Clone)]
        pub struct $name { pub(crate) inner: ck::$name }
        #[pymethods]
        impl $name {
            $($methods)*
            fn __repr__(&self) -> String {
                format!("{}(id={})", stringify!($name), self.inner.id().value())
            }
            #[getter]
            fn id(&self) -> u32 { self.inner.id().value() }
            fn value(&self) -> crate::$module::$row {
                crate::$module::$row { inner: self.inner.value().clone() }
            }
        }
    };
}
view!(BioModel, BioModelRow, bio_values, {
    #[getter]
    fn num_chains(&self) -> usize {
        self.inner.num_chains()
    }
    fn __len__(&self) -> usize {
        self.inner.num_chains()
    }
    #[getter]
    fn chain_span(&self) -> PyResult<crate::bio_values::BioRowSpan> {
        crate::bio_values::BioRowSpan::from_rust(self.inner.chain_span())
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    fn chains(&self) -> Vec<BioChain> {
        self.inner
            .chains()
            .into_iter()
            .map(|inner| BioChain { inner })
            .collect()
    }
    #[getter]
    fn source_model_number(&self) -> Option<i32> {
        self.inner.source_model_number()
    }
});
view!(BioChain, BioChainRow, bio_values, {
    #[getter]
    fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }
    fn __len__(&self) -> usize {
        self.inner.num_residues()
    }
    #[getter]
    fn model_id(&self) -> u32 {
        self.inner.model_id().value()
    }
    #[getter]
    fn entity_id(&self) -> Option<u32> {
        self.inner.entity_id().map(|id| id.value())
    }
    #[getter]
    fn residue_span(&self) -> PyResult<crate::bio_values::BioRowSpan> {
        crate::bio_values::BioRowSpan::from_rust(self.inner.residue_span())
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    fn residues(&self) -> Vec<BioResidue> {
        self.inner
            .residues()
            .into_iter()
            .map(|inner| BioResidue { inner })
            .collect()
    }
    #[getter]
    fn kind(&self) -> crate::bio_values::ChainKind {
        self.inner.kind().into()
    }
    #[getter]
    fn source(&self) -> crate::bio_values::ChainSourceIds {
        crate::bio_values::ChainSourceIds {
            inner: self.inner.source().clone(),
        }
    }
});
view!(BioResidue, BioResidueRow, bio_rows, {
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    fn __len__(&self) -> usize {
        self.inner.num_atoms()
    }
    #[getter]
    fn chain_id(&self) -> u32 {
        self.inner.chain_id().value()
    }
    #[getter]
    fn atom_span(&self) -> PyResult<crate::bio_values::BioRowSpan> {
        crate::bio_values::BioRowSpan::from_rust(self.inner.atom_span())
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    fn residue_info_kind(&self) -> crate::residue_values::ResidueInfoKind {
        self.inner.residue_info_kind().into()
    }
    #[getter]
    fn entity_kind(&self) -> crate::bio_source_values::EntityKind {
        self.inner.entity_kind().into()
    }
    #[getter]
    fn entity_id(&self) -> Option<u32> {
        self.inner.entity_id().map(|id| id.value())
    }
    #[getter]
    fn het_flag(&self) -> Option<u8> {
        self.inner.het_flag()
    }
    #[getter]
    fn source(&self) -> crate::bio_source_values::ResidueSourceIds {
        crate::bio_source_values::ResidueSourceIds {
            inner: self.inner.source().clone(),
        }
    }
    #[getter]
    fn sifts_unp(&self) -> crate::bio_rows::BioSiftsUnpResidue {
        crate::bio_rows::BioSiftsUnpResidue {
            inner: self.inner.sifts_unp(),
        }
    }
    #[getter]
    fn atoms(&self) -> Vec<BioAtom> {
        self.inner
            .atoms()
            .into_iter()
            .map(|inner| BioAtom { inner })
            .collect()
    }
    #[getter]
    fn name(&self) -> String {
        self.inner.name().as_str().to_owned()
    }
    #[getter]
    fn kind(&self) -> crate::bio_source_values::ResidueKind {
        self.inner.kind().into()
    }
});
view!(BioAtom, BioAtomRow, bio_rows, {
    #[getter]
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    #[getter]
    fn residue_id(&self) -> u32 {
        self.inner.residue_id().value()
    }
    #[getter]
    fn isotope_mass_number(&self) -> Option<u16> {
        self.inner.isotope_mass_number()
    }
    #[getter]
    fn altloc(&self) -> Option<crate::bio_source_values::AltLocLabel> {
        self.inner
            .altloc()
            .map(|inner| crate::bio_source_values::AltLocLabel { inner })
    }
    #[getter]
    fn formal_charge(&self) -> i8 {
        self.inner.formal_charge()
    }
    #[getter]
    fn calculation_flag(&self) -> crate::bio_rows::BioCalcFlag {
        self.inner.calculation_flag().into()
    }
    #[getter]
    fn anisou(&self) -> [f64; 6] {
        *self.inner.anisou()
    }
    #[getter]
    fn tls_group_id(&self) -> i16 {
        self.inner.tls_group_id()
    }
    #[getter]
    fn fraction(&self) -> f64 {
        self.inner.fraction()
    }
    #[getter]
    fn source(&self) -> crate::bio_source_values::AtomSourceIds {
        crate::bio_source_values::AtomSourceIds {
            inner: *self.inner.source(),
        }
    }
    #[getter]
    fn position(&self) -> [f64; 3] {
        self.inner.position()
    }
    #[getter]
    fn name(&self) -> crate::bio_values::AtomName {
        crate::bio_values::AtomName {
            inner: self.inner.name(),
        }
    }
    #[getter]
    fn element(&self) -> crate::elements::Element {
        crate::elements::Element {
            inner: self.inner.element(),
        }
    }
    #[getter]
    fn occupancy(&self) -> f64 {
        self.inner.occupancy()
    }
    #[getter]
    fn b_iso(&self) -> f64 {
        self.inner.b_iso()
    }
});
