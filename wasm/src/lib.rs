//! Binding-facing API for compiling COSMolKit to WebAssembly.
//!
//! This crate deliberately contains no chemistry implementation.  It owns the
//! language-boundary shape only and keeps the canonical [`cosmolkit::Molecule`]
//! behind a private field.  Alef can therefore generate wasm-bindgen bindings
//! from this small, stable surface without exposing operation internals or
//! platform-specific Rust types.

#[cfg(feature = "cap-alignment")]
mod alignment;
#[cfg(feature = "cap-aromaticity")]
mod aromaticity;
#[cfg(feature = "cap-batch")]
mod batch;
#[cfg(feature = "cap-batch")]
pub use batch::{BatchRecord, MoleculeBatch};
#[cfg(feature = "cap-inchi")]
mod inchi;
#[cfg(feature = "cap-reaction")]
mod reaction;
#[cfg(feature = "cap-smiles")]
mod smiles;
#[cfg(feature = "cap-reaction")]
pub use reaction::{ReactionApplyResult, reaction_run};

pub use cosmolkit as rust;
/// Complete Rust facade re-export.
///
/// The binding crate mirrors the canonical Rust facade at its root. This
/// keeps every existing `cosmolkit` feature available to Rust consumers while
/// the `Molecule` type below provides the ABI-safe projection used by Alef and
/// wasm-bindgen.
pub use cosmolkit::*;

/// Returns the COSMolKit version through the binding-facing crate.
#[must_use]
pub fn version() -> &'static str {
    cosmolkit::version()
}

// Decode counted chemical text only at the existing language projection.
// Invalid UTF-8 is a binding error; bytes are never replaced or truncated.
fn language_text(text: &cosmolkit::PropertyText) -> Result<String, String> {
    std::str::from_utf8(text.as_bytes())
        .map(str::to_owned)
        .map_err(|error| format!("binding text encoding: {error}"))
}

/// Immutable identity of an element from the canonical public facade.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Element {
    inner: cosmolkit::Element,
}

impl Element {
    /// Returns the source element for this atomic number, or no value.
    pub fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        cosmolkit::Element::from_atomic_number(atomic_number).map(|inner| Self { inner })
    }

    /// Returns the source element for this symbol, including source aliases.
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        cosmolkit::Element::from_symbol(symbol).map(|inner| Self { inner })
    }

    /// Reads the copied value identity without consuming the binding object.
    pub fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }

    /// Reads the canonical source symbol.
    pub fn symbol(&self) -> String {
        self.inner.symbol().to_owned()
    }
}

/// Immutable projection of the public periodic-table result's eight fields.
#[derive(Clone, Copy, Debug)]
pub struct ElementInfo {
    inner: cosmolkit::ElementInfo,
}

impl ElementInfo {
    pub fn element(&self) -> Element {
        Element {
            inner: self.inner.element,
        }
    }

    pub fn symbol(&self) -> String {
        self.inner.symbol.to_owned()
    }

    pub fn atomic_number(&self) -> u8 {
        self.inner.atomic_number
    }

    pub fn period(&self) -> u8 {
        self.inner.period
    }

    pub fn outer_electrons(&self) -> i32 {
        self.inner.outer_electrons
    }

    pub fn valences(&self) -> Vec<i32> {
        self.inner.valences.to_vec()
    }

    pub fn rb0(&self) -> f64 {
        self.inner.rb0
    }

    pub fn atomic_weight(&self) -> f64 {
        self.inner.atomic_weight
    }
}

/// Looks up table metadata through the public facade with source defaults.
#[cfg(feature = "cap-valence")]
pub fn element_info(element: &Element) -> ElementInfo {
    ElementInfo {
        inner: cosmolkit::element_info(element.inner),
    }
}

#[derive(Clone)]
pub struct Molecule {
    inner: std::cell::RefCell<cosmolkit::Molecule>,
}

impl Molecule {
    /// Creates an empty molecule.
    pub fn new() -> Self {
        Self {
            inner: cosmolkit::Molecule::new().into(),
        }
    }

    /// Parses SMILES while explicitly selecting the source sanitization path.
    ///
    /// All other parser options keep their pinned source defaults.
    #[cfg(feature = "cap-smiles")]
    pub fn from_smiles_with_sanitize(
        smiles: &str,
        sanitize: bool,
    ) -> Result<Self, cosmolkit::SmilesError> {
        let params = cosmolkit::SmilesParseParams {
            sanitize,
            ..cosmolkit::SmilesParseParams::default()
        };
        cosmolkit::Molecule::from_smiles_with_params(smiles, &params).map(|inner| Self {
            inner: inner.into(),
        })
    }

    /// Number of atoms in the molecular graph.
    pub fn num_atoms(&self) -> u32 {
        self.inner.borrow().num_atoms() as u32
    }

    /// Number of bonds in the molecular graph.
    pub fn num_bonds(&self) -> u32 {
        self.inner.borrow().num_bonds() as u32
    }

    /// Returns the molecule name, or an empty string when no name is stored.
    pub fn name_or_empty(&self) -> Result<String, String> {
        self.inner
            .borrow()
            .properties()
            .name()
            .map(language_text)
            .transpose()
            .map(|text| text.unwrap_or_default())
    }

    /// Returns a molecule with its name replaced.
    pub fn with_name(&self, name: &str) -> Result<Self, String> {
        self.inner
            .borrow()
            .to_builder()
            .with_name(name.to_owned())
            .build()
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with a string property replaced or inserted.
    pub fn with_property(&self, key: &str, value: &str) -> Result<Self, String> {
        self.inner
            .borrow()
            .to_builder()
            .with_property(key.to_owned(), value.to_owned())
            .and_then(|builder| builder.build())
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with an SDF data field appended.
    pub fn with_sdf_data_field(&self, key: &str, value: &str) -> Result<Self, String> {
        self.inner
            .borrow()
            .to_builder()
            .with_sdf_data_field(key.to_owned(), value.to_owned())
            .build()
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(|error| error.to_string())
    }

    /// Returns a String-tag property, or an empty string when it is absent.
    /// Wrong value kinds and invalid UTF-8 propagate as binding errors.
    pub fn property_or_empty(&self, key: &str) -> Result<String, String> {
        self.inner
            .borrow()
            .property(key)
            .map(|value| {
                value
                    .as_string()
                    .map_err(|error| error.to_string())
                    .and_then(language_text)
            })
            .transpose()
            .map(|text| text.unwrap_or_default())
    }

    /// Returns user/computed property names in stable source insertion order.
    pub fn property_keys(&self) -> Result<Vec<String>, String> {
        self.inner
            .borrow()
            .properties()
            .ordered_props()
            .map(|(key, _)| language_text(key))
            .collect()
    }

    /// Returns SDF data-field names in their source order, including duplicates.
    pub fn sdf_data_field_names(&self) -> Result<Vec<String>, String> {
        self.inner
            .borrow()
            .properties()
            .sdf_data_fields()
            .iter()
            .map(|(name, _)| language_text(name))
            .collect()
    }

    /// Returns the first SDF data-field value with `name`, or an empty string.
    pub fn sdf_data_field_or_empty(&self, name: &str) -> Result<String, String> {
        self.inner
            .borrow()
            .properties()
            .sdf_data_fields()
            .iter()
            .find(|(field_name, _)| field_name.as_bytes() == name.as_bytes())
            .map(|(_, value)| language_text(value))
            .transpose()
            .map(|text| text.unwrap_or_default())
    }

    /// Returns atomic numbers in molecule atom order.
    pub fn atomic_numbers(&self) -> Vec<u8> {
        self.inner
            .borrow()
            .atoms()
            .iter()
            .map(|atom| atom.atomic_number())
            .collect()
    }

    /// Returns the first 2D conformer as a flattened `[x0, y0, ...]` array.
    pub fn coordinates_2d(&self) -> Vec<f64> {
        self.inner
            .borrow()
            .coordinates_2d()
            .map(|coordinates| {
                coordinates
                    .iter()
                    .flat_map(|point| point.iter().copied())
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Returns the number of 3D conformers.
    pub fn num_conformers_3d(&self) -> u32 {
        self.inner.borrow().conformers_3d().len() as u32
    }

    /// Returns a 3D conformer as a flattened `[x0, y0, z0, ...]` array.
    pub fn coordinates_3d(&self, conformer_id: u32) -> Vec<f64> {
        self.inner
            .borrow()
            .conformers_3d()
            .get(conformer_id as usize)
            .map(|conformer| {
                conformer
                    .coordinates()
                    .iter()
                    .flat_map(|point| point.iter().copied())
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Returns a molecule with the default valence cache assigned.
    #[cfg(feature = "cap-valence")]
    pub fn with_assigned_valence(&self) -> Result<Self, cosmolkit::OperationError> {
        self.inner
            .borrow()
            .with_assigned_valence()
            .map(|inner| Self {
                inner: inner.into(),
            })
    }

    /// Returns an independent value prepared with explicit valence policies.
    #[cfg(feature = "cap-valence")]
    pub fn with_assigned_valence_with_params(
        &self,
        params: &cosmolkit::ValenceParams,
    ) -> Result<Self, cosmolkit::OperationError> {
        self.inner
            .borrow()
            .with_assigned_valence_with_params(params)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }

    /// Prepares this value through the canonical in-place operation.
    #[cfg(feature = "cap-valence")]
    pub fn assign_valence_(&self) -> Result<(), cosmolkit::OperationError> {
        self.inner.borrow_mut().assign_valence_()
    }

    /// Prepares this value through the canonical explicit-policy operation.
    #[cfg(feature = "cap-valence")]
    pub fn assign_valence_with_params_(
        &self,
        params: &cosmolkit::ValenceParams,
    ) -> Result<(), cosmolkit::OperationError> {
        self.inner.borrow_mut().assign_valence_with_params_(params)
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(all(test, feature = "full"))]
mod tests {
    use super::Molecule;

    #[test]
    fn valence_projection_keeps_typed_failures_and_preparation_state() {
        let raw = |text| {
            Molecule::from_smiles_with_params(
                text,
                &cosmolkit::SmilesParseParams {
                    sanitize: false,
                    remove_hydrogens: false,
                    ..Default::default()
                },
            )
            .unwrap()
        };
        let value = raw("CCO");
        assert!(value.atom_pair_fingerprint().is_err());
        let prepared = value.with_assigned_valence().unwrap();
        assert!(prepared.atom_pair_fingerprint().is_ok());
        assert!(value.atom_pair_fingerprint().is_err());
        value.assign_valence_().unwrap();
        assert!(value.atom_pair_fingerprint().is_ok());

        let invalid = raw("C(F)(F)(F)(F)F");
        let before = invalid.to_smiles().unwrap();
        assert!(matches!(
            invalid.with_assigned_valence(),
            Err(cosmolkit::OperationError::Valence(_))
        ));
        assert!(matches!(
            invalid.assign_valence_(),
            Err(cosmolkit::OperationError::Valence(_))
        ));
        assert_eq!(invalid.to_smiles().unwrap(), before);
        assert!(invalid.atom_pair_fingerprint().is_err());
        let params = cosmolkit::ValenceParams {
            strict: false,
            ..Default::default()
        };
        let relaxed = invalid.with_assigned_valence_with_params(&params).unwrap();
        assert!(relaxed.atom_pair_fingerprint().is_ok());
        assert!(invalid.atom_pair_fingerprint().is_err());
        invalid.assign_valence_with_params_(&params).unwrap();
        assert!(invalid.atom_pair_fingerprint().is_ok());
    }

    const ETHANOL_SDF: &str = "\
ethanol
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
";

    #[test]
    fn wrapper_forwards_construction_and_counts() {
        let molecule = Molecule::from_smiles("CCO").expect("CCO parses");
        assert_eq!(molecule.num_atoms(), 3);
        assert_eq!(molecule.num_bonds(), 2);
        assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
        assert_eq!(
            Molecule::from_smiles_with_sanitize("CCO", false)
                .expect("unsanitized parse")
                .num_atoms(),
            molecule.num_atoms()
        );
        let from_sdf = Molecule::from_sdf(ETHANOL_SDF).expect("SDF record parses");
        assert_eq!(from_sdf.num_atoms(), 3);
        assert_eq!(from_sdf.num_bonds(), 2);
    }

    #[test]
    fn wrapper_propagates_runtime_errors() {
        assert!(Molecule::from_smiles("[").is_err());
        assert!(Molecule::from_sdf("not an SDF record").is_err());
    }

    #[test]
    fn wrapper_forwards_descriptor_defaults() {
        let molecule = Molecule::from_smiles("CCO").expect("CCO parses");
        let weight = molecule.molecular_weight().expect("average weight");
        assert!((weight - 46.069).abs() < 1e-3, "got {weight}");
        let exact = molecule.exact_molecular_weight().expect("exact weight");
        assert!((exact - 46.041_864_812).abs() < 1e-9, "got {exact}");
        assert_eq!(molecule.molecular_formula().expect("formula"), "C2H6O");
    }

    #[test]
    fn wrapper_round_trips_properties_through_builder() {
        let molecule = Molecule::from_smiles("CCO").expect("CCO parses");
        let named = molecule.with_name("ethanol").expect("name replacement");
        assert_eq!(named.name_or_empty().unwrap(), "ethanol");
        let propertied = named
            .with_property("source", "binding-test")
            .expect("property replacement");
        assert_eq!(
            propertied.property_or_empty("source").unwrap(),
            "binding-test"
        );
        assert_eq!(propertied.property_or_empty("absent").unwrap(), "");
        assert!(
            propertied
                .property_keys()
                .unwrap()
                .contains(&"source".to_owned())
        );
        let with_field = propertied
            .with_sdf_data_field("ID", "ethanol-1")
            .expect("SDF data field append");
        assert_eq!(
            with_field.sdf_data_field_names().unwrap(),
            vec!["ID".to_owned()]
        );
        assert_eq!(
            with_field.sdf_data_field_or_empty("ID").unwrap(),
            "ethanol-1"
        );
        assert_eq!(with_field.sdf_data_field_or_empty("absent").unwrap(), "");
        assert_eq!(named.num_atoms(), molecule.num_atoms());
    }

    #[test]
    fn wrapper_forwards_operations_and_coordinate_accessors() {
        let molecule = Molecule::from_smiles("CCO").expect("CCO parses");
        let with_hydrogens = molecule.with_hydrogens().expect("hydrogen addition");
        assert_eq!(with_hydrogens.num_atoms(), 9);
        assert_eq!(
            with_hydrogens
                .without_hydrogens()
                .expect("hydrogen removal")
                .num_atoms(),
            3
        );
        assert_eq!(molecule.num_conformers_3d(), 0);
        assert!(molecule.coordinates_2d().is_empty());
        let drawn = molecule.with_2d_coordinates().expect("2D layout");
        assert_eq!(drawn.coordinates_2d().len(), 6);
        assert_eq!(molecule.num_atoms(), 3);
        let benzene = Molecule::from_smiles("c1ccccc1").expect("benzene parses");
        benzene
            .with_kekulized_bonds()
            .expect("kekulization defaults");
        benzene
            .with_assigned_aromaticity()
            .expect("aromaticity defaults");
        molecule.sanitize().expect("sanitize defaults");
        molecule.with_assigned_valence().expect("valence defaults");
        molecule.with_assigned_rings().expect("ring defaults");
        molecule
            .with_assigned_ring_families()
            .expect("ring family defaults");
        molecule.with_assigned_radicals().expect("radical defaults");
    }
    #[test]
    fn language_projection_preserves_nul_and_source_order() {
        use cosmolkit::{MoleculeBuilder, MoleculeProperties};
        let props = MoleculeProperties::default()
            .with_name("na\0mé")
            .with_prop("z", "v\0é")
            .unwrap()
            .with_prop("a", "second")
            .unwrap()
            .with_sdf_data_field("ID\0é", "first\0é")
            .with_sdf_data_field("ID\0é", "second");
        let molecule = Molecule {
            inner: MoleculeBuilder::new()
                .with_properties(props)
                .build()
                .unwrap()
                .into(),
        };
        assert_eq!(molecule.name_or_empty().unwrap(), "na\0mé");
        assert_eq!(molecule.property_or_empty("z").unwrap(), "v\0é");
        assert_eq!(molecule.property_keys().unwrap(), ["z", "a"]);
        assert_eq!(molecule.sdf_data_field_names().unwrap(), ["ID\0é", "ID\0é"]);
        assert_eq!(
            molecule.sdf_data_field_or_empty("ID\0é").unwrap(),
            "first\0é"
        );
        assert_eq!(molecule.sdf_data_field_or_empty("ID").unwrap(), "");
    }

    #[test]
    fn language_projection_rejects_opaque_bytes_and_wrong_property_kinds() {
        use cosmolkit::{MoleculeBuilder, MoleculeProperties, PropertyText, PropertyValue};
        let opaque = PropertyText::from_bytes(b"a\xff\0");
        let props = MoleculeProperties::default()
            .with_name(opaque.clone())
            .with_prop("text", opaque.clone())
            .unwrap()
            .with_prop(opaque.clone(), "value")
            .unwrap()
            .with_prop("native", PropertyValue::Int(1))
            .unwrap()
            .with_sdf_data_field("ID", opaque.clone())
            .with_sdf_data_field(opaque.clone(), "value");
        let inner = MoleculeBuilder::new()
            .with_properties(props)
            .build()
            .unwrap();
        let before = inner.clone();
        let molecule = Molecule {
            inner: inner.into(),
        };
        assert!(
            molecule
                .name_or_empty()
                .unwrap_err()
                .starts_with("binding text encoding:")
        );
        assert!(molecule.property_or_empty("text").is_err());
        assert!(
            molecule
                .property_or_empty("native")
                .unwrap_err()
                .contains("Int")
        );
        assert!(molecule.property_keys().is_err());
        assert!(molecule.sdf_data_field_names().is_err());
        assert!(molecule.sdf_data_field_or_empty("ID").is_err());
        assert_eq!(molecule.property_or_empty("absent").unwrap(), "");
        assert_eq!(*molecule.inner.borrow(), before);
    }
}

#[cfg(all(test, feature = "full"))]
mod bio_read_tests;

#[cfg(all(test, feature = "full"))]
mod bio_hierarchy_tests;

#[cfg(all(test, feature = "full"))]
mod bio_parts_tests;

#[cfg(all(test, feature = "full"))]
mod bio_selection_tests;

#[cfg(all(test, feature = "full"))]
mod bio_writers_tests;

#[cfg(all(test, feature = "full"))]
mod bio_protein_tests;

#[cfg(feature = "cap-conformer")]
mod conformer;
#[cfg(feature = "cap-conformer")]
pub use conformer::{EmbedMoleculeResult, EmbedMultipleConfsResult};

#[cfg(all(test, feature = "cap-conformer"))]
mod conformer_tests;

#[cfg(feature = "cap-depict")]
mod depict;

#[cfg(feature = "cap-descriptors")]
mod descriptors;

#[cfg(feature = "cap-fingerprints")]
mod path_codes;
#[cfg(feature = "cap-fingerprints")]
pub use path_codes::AtomPairAtomCodeResult;

#[cfg(feature = "cap-fingerprints")]
mod atom_pair;

#[cfg(feature = "cap-fingerprints")]
mod morgan;
#[cfg(feature = "cap-fingerprints")]
pub use morgan::{
    morgan_generator_counts, morgan_generator_fingerprints, morgan_generator_sparse_counts,
    morgan_generator_sparse_fingerprints,
};

#[cfg(feature = "cap-fingerprints")]
mod topological_torsion;
#[cfg(feature = "cap-fingerprints")]
pub use topological_torsion::{
    topological_torsion_generator_counts, topological_torsion_generator_fingerprints,
    topological_torsion_generator_sparse_counts, topological_torsion_generator_sparse_fingerprints,
};

#[cfg(feature = "cap-fingerprints")]
mod maccs;

#[cfg(feature = "cap-fingerprints")]
mod layered_pattern;

#[cfg(feature = "cap-fingerprints")]
mod topological;

#[cfg(feature = "cap-forcefields")]
mod forcefield_properties;

#[cfg(feature = "cap-forcefields")]
mod uff;
#[cfg(feature = "cap-forcefields")]
pub use uff::{UffConformerOptimizationResult, UffOptimizationResult};

#[cfg(feature = "cap-forcefields")]
mod mmff;
#[cfg(feature = "cap-forcefields")]
pub use mmff::{MmffOptimizeMoleculeConfsResult, MmffOptimizeMoleculeResult};

#[cfg(feature = "cap-hashing")]
mod hashing;

#[cfg(feature = "cap-hydrogens")]
mod hydrogens;

#[cfg(all(feature = "cap-bio", feature = "cap-io"))]
mod bio_conversion;
#[cfg(all(feature = "cap-bio", feature = "cap-io"))]
pub use bio_conversion::{
    bio_structure_to_molecule, bio_structure_to_molecule_with_params, protein_to_molecule,
    protein_to_molecule_with_params,
};

#[cfg(feature = "cap-io")]
mod sdf_reading;
#[cfg(feature = "cap-io")]
pub use sdf_reading::{SdfGraph, SdfRecord};

#[cfg(feature = "cap-io")]
mod property_strings;

#[cfg(feature = "cap-io")]
mod xyz_mol2;

#[cfg(feature = "cap-io")]
mod mol_sdf_writing;

#[cfg(feature = "cap-io")]
mod sdf_datasets;
#[cfg(feature = "cap-io")]
pub use sdf_datasets::{SdfDataset, SdfDatasetIterator, SdfReader, SdfRecordStream};

#[cfg(all(feature = "cap-io", feature = "cap-batch"))]
mod sdf_batches;
#[cfg(all(feature = "cap-io", feature = "cap-batch"))]
pub use sdf_batches::{SdfBatchIterator, SdfReaderBatchIterator};

#[cfg(feature = "cap-kekulize")]
mod kekulize;

#[cfg(feature = "cap-matrices")]
mod matrices;

#[cfg(feature = "cap-radicals")]
mod radicals;

#[cfg(feature = "cap-rings")]
mod rings;

#[cfg(feature = "cap-sanitize")]
mod sanitize;

#[cfg(feature = "cap-search")]
#[path = "search.rs"]
mod search_projection;
