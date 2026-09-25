//! Binding-facing API for compiling COSMolKit to WebAssembly.
//!
//! This crate deliberately contains no chemistry implementation.  It owns the
//! language-boundary shape only and keeps the canonical [`cosmolkit::Molecule`]
//! behind a private field.  Alef can therefore generate wasm-bindgen bindings
//! from this small, stable surface without exposing operation internals or
//! platform-specific Rust types.

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

#[derive(Clone)]
pub struct Molecule {
    inner: cosmolkit::Molecule,
}

impl Molecule {
    /// Creates an empty molecule.
    pub fn new() -> Self {
        Self {
            inner: cosmolkit::Molecule::new(),
        }
    }

    /// Parses a SMILES string into a binding-owned molecule value.
    pub fn from_smiles(smiles: &str) -> Result<Self, String> {
        cosmolkit::Molecule::from_smiles(smiles)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Parses SMILES while explicitly selecting the source sanitization path.
    ///
    /// All other parser options keep their pinned source defaults.
    pub fn from_smiles_with_sanitize(smiles: &str, sanitize: bool) -> Result<Self, String> {
        let params = cosmolkit::SmilesParseParams {
            sanitize,
            ..cosmolkit::SmilesParseParams::default()
        };
        cosmolkit::Molecule::from_smiles_with_params(smiles, &params)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Reads the first SDF record using the default source policy.
    pub fn from_sdf(text: &str) -> Result<Self, String> {
        cosmolkit::Molecule::from_sdf(text)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Number of atoms in the molecular graph.
    pub fn num_atoms(&self) -> u32 {
        self.inner.num_atoms() as u32
    }

    /// Number of bonds in the molecular graph.
    pub fn num_bonds(&self) -> u32 {
        self.inner.num_bonds() as u32
    }

    /// Returns the molecule name, or an empty string when no name is stored.
    pub fn name_or_empty(&self) -> String {
        self.inner
            .properties()
            .name()
            .unwrap_or_default()
            .to_owned()
    }

    /// Returns a molecule with its name replaced.
    pub fn with_name(&self, name: &str) -> Result<Self, String> {
        self.inner
            .to_builder()
            .with_name(name.to_owned())
            .build()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with a string property replaced or inserted.
    pub fn with_property(&self, key: &str, value: &str) -> Result<Self, String> {
        self.inner
            .to_builder()
            .with_property(key.to_owned(), value.to_owned())
            .and_then(|builder| builder.build())
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with an SDF data field appended.
    pub fn with_sdf_data_field(&self, key: &str, value: &str) -> Result<Self, String> {
        self.inner
            .to_builder()
            .with_sdf_data_field(key.to_owned(), value.to_owned())
            .build()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a property value, or an empty string when it is absent.
    pub fn property_or_empty(&self, key: &str) -> String {
        self.inner.property(key).unwrap_or_default().to_owned()
    }

    /// Returns user/computed property names in stable source order.
    pub fn property_keys(&self) -> Vec<String> {
        self.inner.properties().props().keys().cloned().collect()
    }

    /// Returns SDF data-field names in their source order.
    pub fn sdf_data_field_names(&self) -> Vec<String> {
        self.inner
            .properties()
            .sdf_data_fields()
            .iter()
            .map(|(name, _)| name.clone())
            .collect()
    }

    /// Returns the first SDF data-field value with `name`, or an empty string.
    pub fn sdf_data_field_or_empty(&self, name: &str) -> String {
        self.inner
            .properties()
            .sdf_data_fields()
            .iter()
            .find_map(|(field_name, value)| (field_name == name).then_some(value.as_str()))
            .unwrap_or_default()
            .to_owned()
    }

    /// Returns atomic numbers in molecule atom order.
    pub fn atomic_numbers(&self) -> Vec<u8> {
        self.inner
            .atoms()
            .iter()
            .map(|atom| atom.atomic_number())
            .collect()
    }

    /// Returns the first 2D conformer as a flattened `[x0, y0, ...]` array.
    pub fn coordinates_2d(&self) -> Vec<f64> {
        self.inner
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
        self.inner.conformers_3d().len() as u32
    }

    /// Returns a 3D conformer as a flattened `[x0, y0, z0, ...]` array.
    pub fn coordinates_3d(&self, conformer_id: u32) -> Vec<f64> {
        self.inner
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

    /// Calculates the average molecular weight with pinned source defaults.
    pub fn molecular_weight(&self) -> Result<f64, String> {
        self.inner
            .molecular_weight()
            .map_err(|error| error.to_string())
    }

    /// Calculates the exact molecular weight with pinned source defaults.
    pub fn exact_molecular_weight(&self) -> Result<f64, String> {
        self.inner
            .exact_molecular_weight()
            .map_err(|error| error.to_string())
    }

    /// Returns the Hill-ordered molecular formula with source defaults.
    pub fn molecular_formula(&self) -> Result<String, String> {
        self.inner
            .molecular_formula()
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with explicit hydrogens added through the operation contract.
    pub fn with_hydrogens(&self) -> Result<Self, String> {
        self.inner
            .with_hydrogens()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with explicit hydrogens removed through the operation contract.
    pub fn without_hydrogens(&self) -> Result<Self, String> {
        self.inner
            .without_hydrogens()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with aromatic bonds kekulized under source defaults.
    pub fn with_kekulized_bonds(&self) -> Result<Self, String> {
        self.inner
            .with_kekulized_bonds()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a sanitized molecule using the default operation pipeline.
    pub fn sanitize(&self) -> Result<Self, String> {
        self.inner
            .sanitize()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with the default valence cache assigned.
    pub fn with_assigned_valence(&self) -> Result<Self, String> {
        self.inner
            .with_assigned_valence()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with the default ring cache assigned.
    pub fn with_assigned_rings(&self) -> Result<Self, String> {
        self.inner
            .with_assigned_rings()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with ring-family data assigned.
    pub fn with_assigned_ring_families(&self) -> Result<Self, String> {
        self.inner
            .with_assigned_ring_families()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with default aromaticity assigned.
    pub fn with_assigned_aromaticity(&self) -> Result<Self, String> {
        self.inner
            .with_assigned_aromaticity()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns a molecule with default radical assignments applied.
    pub fn with_assigned_radicals(&self) -> Result<Self, String> {
        self.inner
            .with_assigned_radicals()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Computes and returns a molecule with direct 2D coordinates.
    pub fn with_2d_coordinates(&self) -> Result<Self, String> {
        self.inner
            .with_2d_coordinates()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
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
        assert_eq!(named.name_or_empty(), "ethanol");
        let propertied = named
            .with_property("source", "binding-test")
            .expect("property replacement");
        assert_eq!(propertied.property_or_empty("source"), "binding-test");
        assert_eq!(propertied.property_or_empty("absent"), "");
        assert!(propertied.property_keys().contains(&"source".to_owned()));
        let with_field = propertied
            .with_sdf_data_field("ID", "ethanol-1")
            .expect("SDF data field append");
        assert_eq!(with_field.sdf_data_field_names(), vec!["ID".to_owned()]);
        assert_eq!(with_field.sdf_data_field_or_empty("ID"), "ethanol-1");
        assert_eq!(with_field.sdf_data_field_or_empty("absent"), "");
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
}
