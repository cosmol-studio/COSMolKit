//! Native owned hierarchy snapshots; no relationship resolution in Python.
use cosmolkit as ck;
use pyo3::prelude::*;

macro_rules! view {
    ($name:ident, $row:ident, $module:ident, {$($methods:tt)*}) => {
        #[doc = concat!("Owned ", stringify!($name), " view retaining its structure version.\n\n",
            "Hierarchy and coordinate reads use that retained structure. Editing another ",
            "COW structure value does not invalidate this object or alter its values. ",
            "value() returns the detached row without its owner context.")]
        #[pyclass(module = "cosmolkit", frozen, from_py_object)]
        #[derive(Clone)]
        pub struct $name { pub(crate) inner: ck::$name }
        #[pymethods]
        impl $name {
            $($methods)*
            fn __repr__(&self) -> String {
                format!("{}(id={})", stringify!($name), self.inner.id().value())
            }
            /// Stable local row identifier in the retained structure, not a source-file serial number.
            #[getter]
            fn id(&self) -> u32 { self.inner.id().value() }
            /// Return a detached row copy without structure ownership or contextual traversal.
            fn value(&self) -> crate::$module::$row {
                crate::$module::$row { inner: self.inner.value().clone() }
            }
        }
    };
}
view!(BioModel, BioModelRow, bio_values, {
    #[getter]
    /// Number of chains in this model.
    fn num_chains(&self) -> usize {
        self.inner.num_chains()
    }
    fn __len__(&self) -> usize {
        self.inner.num_chains()
    }
    #[getter]
    /// Half-open span of this model's chain rows in the retained structure.
    fn chain_span(&self) -> PyResult<crate::bio_values::BioRowSpan> {
        crate::bio_values::BioRowSpan::from_rust(self.inner.chain_span())
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    /// Chains in stored order, each retaining the same structure version.
    fn chains(&self) -> Vec<BioChain> {
        self.inner
            .chains()
            .into_iter()
            .map(|inner| BioChain { inner })
            .collect()
    }
    #[getter]
    /// Model number recorded by the input format, or None if it was not supplied.
    fn source_model_number(&self) -> Option<i32> {
        self.inner.source_model_number()
    }
});
view!(BioChain, BioChainRow, bio_values, {
    #[getter]
    /// Number of residue rows in this chain.
    fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }
    fn __len__(&self) -> usize {
        self.inner.num_residues()
    }
    #[getter]
    /// Local identifier of the model containing this chain.
    fn model_id(&self) -> u32 {
        self.inner.model_id().value()
    }
    #[getter]
    /// Associated entity's local identifier, or None when no entity is assigned.
    /// Associated entity's local row identifier, or None if unassigned.
    fn entity_id(&self) -> Option<u32> {
        self.inner.entity_id().map(|id| id.value())
    }
    #[getter]
    /// Half-open span of this chain's residue rows in the retained structure.
    fn residue_span(&self) -> PyResult<crate::bio_values::BioRowSpan> {
        crate::bio_values::BioRowSpan::from_rust(self.inner.residue_span())
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    /// Residues in stored order, each retaining the same structure version.
    fn residues(&self) -> Vec<BioResidue> {
        self.inner
            .residues()
            .into_iter()
            .map(|inner| BioResidue { inner })
            .collect()
    }
    #[getter]
    /// Stored chain classification; this read does not infer a new classification.
    fn kind(&self) -> crate::bio_values::ChainKind {
        self.inner.kind().into()
    }
    #[getter]
    /// Source-format chain identifiers, distinct from local row identifiers.
    fn source(&self) -> crate::bio_values::ChainSourceIds {
        crate::bio_values::ChainSourceIds {
            inner: self.inner.source().clone(),
        }
    }
});
view!(BioResidue, BioResidueRow, bio_rows, {
    #[getter]
    /// Number of atom rows in this residue, including retained alternate locations.
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    fn __len__(&self) -> usize {
        self.inner.num_atoms()
    }
    #[getter]
    /// Local identifier of the chain containing this residue.
    fn chain_id(&self) -> u32 {
        self.inner.chain_id().value()
    }
    #[getter]
    /// Half-open span of this residue's atom rows in the retained structure.
    fn atom_span(&self) -> PyResult<crate::bio_values::BioRowSpan> {
        crate::bio_values::BioRowSpan::from_rust(self.inner.atom_span())
            .map_err(crate::bio_errors::BioStructureError::from_rust)
    }
    #[getter]
    /// Chemical residue classification from the canonical residue information table.
    fn residue_info_kind(&self) -> crate::residue_values::ResidueInfoKind {
        self.inner.residue_info_kind().into()
    }
    #[getter]
    /// Stored entity category associated with this residue.
    fn entity_kind(&self) -> crate::bio_source_values::EntityKind {
        self.inner.entity_kind().into()
    }
    #[getter]
    /// Local entity identifier associated with this residue, or None when unassigned.
    fn entity_id(&self) -> Option<u32> {
        self.inner.entity_id().map(|id| id.value())
    }
    #[getter]
    /// Source ATOM/HETATM classification marker, or None when unavailable.
    fn het_flag(&self) -> Option<u8> {
        self.inner.het_flag()
    }
    #[getter]
    /// Source-file identifiers retained separately from local row indices.
    fn source(&self) -> crate::bio_source_values::ResidueSourceIds {
        crate::bio_source_values::ResidueSourceIds {
            inner: self.inner.source().clone(),
        }
    }
    #[getter]
    /// Stored SIFTS UniProt mapping for this residue.
    fn sifts_unp(&self) -> crate::bio_rows::BioSiftsUnpResidue {
        crate::bio_rows::BioSiftsUnpResidue {
            inner: self.inner.sifts_unp(),
        }
    }
    #[getter]
    /// Atom views in stored residue order, each retaining the same structure version.
    fn atoms(&self) -> Vec<BioAtom> {
        self.inner
            .atoms()
            .into_iter()
            .map(|inner| BioAtom { inner })
            .collect()
    }
    #[getter]
    /// Source atom or residue name retained by this hierarchy object.
    fn name(&self) -> String {
        self.inner.name().as_str().to_owned()
    }
    #[getter]
    /// Stored residue classification.
    fn kind(&self) -> crate::bio_source_values::ResidueKind {
        self.inner.kind().into()
    }
});
view!(BioAtom, BioAtomRow, bio_rows, {
    #[getter]
    /// Atomic number of this atom's element.
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    #[getter]
    /// Local identifier of the residue containing this atom.
    fn residue_id(&self) -> u32 {
        self.inner.residue_id().value()
    }
    #[getter]
    /// Explicit isotope mass number, or None if unspecified.
    fn isotope_mass_number(&self) -> Option<u16> {
        self.inner.isotope_mass_number()
    }
    #[getter]
    /// Alternate-location label, or None for an unlabeled location.
    fn altloc(&self) -> Option<crate::bio_source_values::AltLocLabel> {
        self.inner
            .altloc()
            .map(|inner| crate::bio_source_values::AltLocLabel { inner })
    }
    #[getter]
    /// Stored signed formal charge of the atom.
    fn formal_charge(&self) -> i8 {
        self.inner.formal_charge()
    }
    #[getter]
    /// Source calculation status of this atom.
    fn calculation_flag(&self) -> crate::bio_rows::BioCalcFlag {
        self.inner.calculation_flag().into()
    }
    #[getter]
    /// Anisotropic displacement components in native order: U11, U22, U33, U12, U13, U23.
    fn anisou(&self) -> [f64; 6] {
        *self.inner.anisou()
    }
    #[getter]
    /// Stored translation/libration/screw group identifier.
    fn tls_group_id(&self) -> i16 {
        self.inner.tls_group_id()
    }
    #[getter]
    /// Stored atom fraction value from the source structure.
    fn fraction(&self) -> f64 {
        self.inner.fraction()
    }
    #[getter]
    /// Source-file identifiers retained separately from local row indices.
    fn source(&self) -> crate::bio_source_values::AtomSourceIds {
        crate::bio_source_values::AtomSourceIds {
            inner: *self.inner.source(),
        }
    }
    #[getter]
    /// Cartesian position (x, y, z) from the retained structure's coordinate block.
    fn position(&self) -> [f64; 3] {
        self.inner.position()
    }
    #[getter]
    /// Source atom or residue name retained by this hierarchy object.
    fn name(&self) -> crate::bio_values::AtomName {
        crate::bio_values::AtomName {
            inner: self.inner.name(),
        }
    }
    #[getter]
    /// Element of this atom.
    fn element(&self) -> crate::elements::Element {
        crate::elements::Element {
            inner: self.inner.element(),
        }
    }
    #[getter]
    /// Stored site occupancy.
    fn occupancy(&self) -> f64 {
        self.inner.occupancy()
    }
    #[getter]
    /// Stored isotropic B factor.
    fn b_iso(&self) -> f64 {
        self.inner.b_iso()
    }
});
