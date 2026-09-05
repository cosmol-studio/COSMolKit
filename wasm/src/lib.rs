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

fn fingerprint_bits(fingerprint: cosmolkit::Fingerprint) -> Vec<u32> {
    fingerprint
        .on_bits()
        .into_iter()
        .map(|bit| bit as u32)
        .collect()
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
    pub fn from_smiles_with_sanitize(smiles: &str, sanitize: bool) -> Result<Self, String> {
        cosmolkit::Molecule::from_smiles_with_sanitize(smiles, sanitize)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Parses a MOL/SDF block into a binding-owned molecule value.
    pub fn from_mol_block(block: &str) -> Result<Self, String> {
        cosmolkit::Molecule::from_mol_block(block)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Parses an XYZ block into a binding-owned molecule value.
    pub fn from_xyz_block(block: &str) -> Result<Self, String> {
        cosmolkit::Molecule::from_xyz_block(block)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Parses a PDB block using the default RDKit-compatible structure policy.
    pub fn from_pdb_block(block: &str) -> Result<Self, String> {
        cosmolkit::Molecule::from_pdb_block(block)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Parses an mmCIF block using the default RDKit-compatible structure policy.
    pub fn from_mmcif_block(block: &str) -> Result<Self, String> {
        cosmolkit::Molecule::from_mmcif_block(block)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Parses a protein one-letter sequence into a molecule.
    pub fn from_protein_sequence(sequence: &str) -> Result<Self, String> {
        cosmolkit::mol_from_sequence(sequence)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Parses a DNA or RNA one-letter sequence into a molecule.
    pub fn from_nucleic_sequence(sequence: &str) -> Result<Self, String> {
        cosmolkit::mol_from_sequence_with_type(sequence, false)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Parses a Tripos MOL2 block with the default parser policy.
    pub fn from_mol2_block(block: &str) -> Result<Option<Self>, String> {
        cosmolkit::mol_from_mol2_block_like_rdkit(block, cosmolkit::Mol2ReadParams::default())
            .map(|record| {
                record.map(|record| Self {
                    inner: record.molecule,
                })
            })
            .map_err(|error| error.to_string())
    }

    /// Parses a MOL2 block and turns the source no-graph result into an
    /// explicit binding error instead of silently returning `None`.
    pub fn from_mol2_or_error(block: &str) -> Result<Self, String> {
        Self::from_mol2_block(block)?
            .ok_or_else(|| "MOL2 input did not produce a molecule".to_owned())
    }

    /// Parses an InChI and preserves the source API's no-graph result.
    pub fn from_inchi(
        inchi: &str,
        sanitize: bool,
        remove_hs: bool,
    ) -> Result<Option<Self>, String> {
        cosmolkit::mol_from_inchi(inchi.as_bytes(), sanitize, remove_hs)
            .map(|output| output.molecule.map(|inner| Self { inner }))
            .map_err(|error| error.to_string())
    }

    /// Parses an InChI and makes the source no-graph result explicit.
    pub fn from_inchi_or_error(
        inchi: &str,
        sanitize: bool,
        remove_hs: bool,
    ) -> Result<Self, String> {
        Self::from_inchi(inchi, sanitize, remove_hs)?
            .ok_or_else(|| "InChI input did not produce a molecule".to_owned())
    }

    /// Writes an isomeric SMILES string without mutating this molecule.
    pub fn to_smiles(&self) -> Result<String, String> {
        self.inner
            .to_smiles(true)
            .map_err(|error| error.to_string())
    }

    /// Writes SMILES with an explicit isomeric-stereo switch.
    pub fn to_smiles_with_isomeric(&self, isomeric: bool) -> Result<String, String> {
        self.inner
            .to_smiles(isomeric)
            .map_err(|error| error.to_string())
    }

    /// Writes the molecule as a V2000 MOL block.
    pub fn to_mol_block_v2000(&self) -> Result<String, String> {
        cosmolkit::io::molblock::mol_to_v2000_block(&self.inner).map_err(|error| error.to_string())
    }

    /// Writes the molecule as a V3000 MOL block.
    pub fn to_mol_block_v3000(&self) -> Result<String, String> {
        cosmolkit::io::molblock::mol_to_v3000_block(&self.inner).map_err(|error| error.to_string())
    }

    /// Writes a V2000 SDF record.
    pub fn to_sdf_record_v2000(&self) -> Result<String, String> {
        cosmolkit::io::molblock::mol_to_2d_sdf_record(
            &self.inner,
            cosmolkit::io::molblock::SdfFormat::V2000,
        )
        .map_err(|error| error.to_string())
    }

    /// Writes a V3000 SDF record.
    pub fn to_sdf_record_v3000(&self) -> Result<String, String> {
        cosmolkit::io::molblock::mol_to_2d_sdf_record(
            &self.inner,
            cosmolkit::io::molblock::SdfFormat::V3000,
        )
        .map_err(|error| error.to_string())
    }

    /// Writes the molecule as SMARTS using default writer parameters.
    pub fn to_smarts(&self) -> Result<String, String> {
        cosmolkit::mol_to_smarts(&self.inner, &cosmolkit::SmilesWriteParams::default())
            .map_err(|error| error.to_string())
    }

    /// Writes the molecule as CXSMARTS using default writer parameters.
    pub fn to_cx_smarts(&self) -> Result<String, String> {
        cosmolkit::mol_to_cx_smarts(&self.inner, &cosmolkit::SmilesWriteParams::default())
            .map_err(|error| error.to_string())
    }

    /// Writes the source-backed InChI string without mutating this molecule.
    pub fn to_inchi(&self) -> Result<String, String> {
        cosmolkit::mol_to_inchi(&self.inner, None)
            .map(|output| String::from_utf8_lossy(&output.inchi).into_owned())
            .map_err(|error| error.to_string())
    }

    /// Writes the source-backed InChIKey without mutating this molecule.
    pub fn to_inchi_key(&self) -> Result<String, String> {
        cosmolkit::mol_to_inchi_key(&self.inner, None)
            .map(|output| String::from_utf8_lossy(&output.key).into_owned())
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
    pub fn with_name(&self, name: &str) -> Self {
        Self {
            inner: self.inner.with_name(name),
        }
    }

    /// Returns a molecule with a string property replaced or inserted.
    pub fn with_property(&self, key: &str, value: &str) -> Self {
        Self {
            inner: self.inner.with_prop(key, value),
        }
    }

    /// Returns a molecule with an SDF data field appended.
    pub fn with_sdf_data_field(&self, key: &str, value: &str) -> Self {
        Self {
            inner: self.inner.with_sdf_data_field(key, value),
        }
    }

    /// Returns a property value, or an empty string when it is absent.
    pub fn property_or_empty(&self, key: &str) -> String {
        self.inner.prop(key).unwrap_or_default().to_owned()
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

    /// Returns the source coordinate dimension as `"2D"`, `"3D"`, or `""`.
    pub fn source_coordinate_dimension_or_empty(&self) -> String {
        match self.inner.source_coordinate_dim() {
            Some(cosmolkit::CoordinateDimension::TwoD) => "2D".to_owned(),
            Some(cosmolkit::CoordinateDimension::ThreeD) => "3D".to_owned(),
            None => String::new(),
        }
    }

    /// Returns atomic numbers in molecule atom order.
    pub fn atomic_numbers(&self) -> Vec<u8> {
        self.inner.atomic_numbers()
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

    /// Renders this molecule as SVG at the requested pixel dimensions.
    pub fn to_svg(&self, width: u32, height: u32) -> Result<String, String> {
        self.inner
            .to_svg(width, height)
            .map_err(|error| error.to_string())
    }

    /// Renders this molecule as PNG bytes at the requested pixel dimensions.
    pub fn to_png(&self, width: u32, height: u32) -> Result<Vec<u8>, String> {
        self.inner
            .to_png(width, height)
            .map_err(|error| error.to_string())
    }

    /// Calculates RDKit-compatible average molecular weight including hydrogen.
    pub fn molecular_weight(&self) -> Result<f64, String> {
        cosmolkit::calc_mol_wt(&self.inner, false).map_err(|error| error.to_string())
    }

    /// Calculates RDKit-compatible exact molecular weight including hydrogen.
    pub fn exact_molecular_weight(&self) -> Result<f64, String> {
        cosmolkit::calc_exact_mol_wt(&self.inner, false).map_err(|error| error.to_string())
    }

    /// Returns the Crippen logP value.
    pub fn crippen_log_p(&self) -> Result<f64, String> {
        cosmolkit::calc_crippen_descriptors(&self.inner, false, false)
            .map(|values| values.logp)
            .map_err(|error| error.to_string())
    }

    /// Returns the Crippen molar refractivity value.
    pub fn crippen_molar_refractivity(&self) -> Result<f64, String> {
        cosmolkit::calc_crippen_descriptors(&self.inner, false, false)
            .map(|values| values.molar_refractivity)
            .map_err(|error| error.to_string())
    }

    /// Calculates topological polar surface area.
    pub fn tpsa(&self) -> Result<f64, String> {
        cosmolkit::calc_tpsa(&self.inner, false, false).map_err(|error| error.to_string())
    }

    /// Counts Lipinski hydrogen-bond acceptors.
    pub fn h_bond_acceptors(&self) -> Result<u32, String> {
        cosmolkit::calc_num_hba(&self.inner).map_err(|error| error.to_string())
    }

    /// Counts Lipinski hydrogen-bond donors.
    pub fn h_bond_donors(&self) -> Result<u32, String> {
        cosmolkit::calc_num_hbd(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the Hill-formula molecular formula.
    pub fn formula(&self) -> Result<String, String> {
        cosmolkit::calc_mol_formula(&self.inner, false, false).map_err(|error| error.to_string())
    }

    /// Returns the fraction of sp3-hybridized carbon atoms.
    pub fn fraction_csp3(&self) -> Result<f64, String> {
        cosmolkit::calc_fraction_csp3(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the number of heavy atoms.
    pub fn num_heavy_atoms(&self) -> Result<u32, String> {
        cosmolkit::calc_num_heavy_atoms(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the number of rings.
    pub fn num_rings(&self) -> Result<u32, String> {
        cosmolkit::calc_num_rings(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the number of aromatic rings.
    pub fn num_aromatic_rings(&self) -> Result<u32, String> {
        cosmolkit::calc_num_aromatic_rings(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the number of rotatable bonds under the default source policy.
    pub fn num_rotatable_bonds(&self) -> Result<u32, String> {
        cosmolkit::calc_num_rotatable_bonds(
            &self.inner,
            cosmolkit::NumRotatableBondsOptions::Default,
        )
        .map_err(|error| error.to_string())
    }

    /// Returns the Chi connectivity descriptor of the requested order.
    pub fn chi(&self, order: u32, path_type: bool, force: bool) -> Result<f64, String> {
        let order = order as usize;
        if path_type {
            cosmolkit::calc_chi_nv(&self.inner, order, force)
        } else {
            cosmolkit::calc_chi_nn(&self.inner, order, force)
        }
        .map_err(|error| error.to_string())
    }

    /// Returns the zero-order Chi connectivity descriptor.
    pub fn chi_0(&self) -> f64 {
        cosmolkit::calc_chi_0(&self.inner)
    }

    /// Returns the first-order Chi connectivity descriptor.
    pub fn chi_1(&self) -> f64 {
        cosmolkit::calc_chi_1(&self.inner)
    }

    /// Returns the Hall-Kier Phi shape descriptor.
    pub fn phi(&self) -> Result<f64, String> {
        cosmolkit::calc_phi(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the number of spiro atoms.
    pub fn num_spiro_atoms(&self) -> Result<u32, String> {
        cosmolkit::calc_num_spiro_atoms(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the number of bridgehead atoms.
    pub fn num_bridgehead_atoms(&self) -> Result<u32, String> {
        cosmolkit::calc_num_bridgehead_atoms(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns Labute's approximate surface area.
    pub fn labute_asa(&self, include_hydrogens: bool, force: bool) -> f64 {
        cosmolkit::calc_labute_asa(&self.inner, include_hydrogens, force)
    }

    /// Returns the number of atom stereo centers.
    pub fn num_atom_stereo_centers(&self) -> Result<u32, String> {
        cosmolkit::calc_num_atom_stereo_centers(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the number of unspecified atom stereo centers.
    pub fn num_unspecified_atom_stereo_centers(&self) -> Result<u32, String> {
        cosmolkit::calc_num_unspecified_atom_stereo_centers(&self.inner)
            .map_err(|error| error.to_string())
    }

    /// Returns the QED score.
    pub fn qed(&self) -> Result<f64, String> {
        cosmolkit::calc_qed(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the Hall-Kier alpha descriptor.
    pub fn hall_kier_alpha(&self) -> f64 {
        cosmolkit::calc_hall_kier_alpha(&self.inner)
    }

    /// Returns the first kappa shape descriptor.
    pub fn kappa_1(&self) -> f64 {
        cosmolkit::calc_kappa_1(&self.inner)
    }

    /// Returns the second kappa shape descriptor.
    pub fn kappa_2(&self) -> Result<f64, String> {
        cosmolkit::calc_kappa_2(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the third kappa shape descriptor.
    pub fn kappa_3(&self) -> Result<f64, String> {
        cosmolkit::calc_kappa_3(&self.inner).map_err(|error| error.to_string())
    }

    /// Returns the set bits of the default 2,048-bit Pattern fingerprint.
    pub fn pattern_fingerprint(&self) -> Result<Vec<u32>, String> {
        let fingerprint = self
            .inner
            .pattern_fingerprint(&cosmolkit::PatternFingerprintParams::default())
            .map_err(|error| error.to_string())?;
        Ok(fingerprint
            .on_bits()
            .into_iter()
            .map(|bit| bit as u32)
            .collect())
    }

    /// Returns the default Morgan fingerprint set bits.
    pub fn morgan_fingerprint(&self) -> Result<Vec<u32>, String> {
        self.inner
            .morgan_fingerprint(&cosmolkit::MorganFingerprintParams::default())
            .map(fingerprint_bits)
            .map_err(|error| error.to_string())
    }

    /// Returns the default AtomPair fingerprint set bits.
    pub fn atom_pair_fingerprint(&self) -> Result<Vec<u32>, String> {
        self.inner
            .atom_pair_fingerprint(&cosmolkit::AtomPairFingerprintParams::default())
            .map(fingerprint_bits)
            .map_err(|error| error.to_string())
    }

    /// Returns the default Layered fingerprint set bits.
    pub fn layered_fingerprint(&self) -> Result<Vec<u32>, String> {
        self.inner
            .layered_fingerprint(&cosmolkit::LayeredFingerprintParams::default())
            .map(fingerprint_bits)
            .map_err(|error| error.to_string())
    }

    /// Returns the default topological fingerprint set bits.
    pub fn topological_fingerprint(&self) -> Result<Vec<u32>, String> {
        self.inner
            .topological_fingerprint(&cosmolkit::TopologicalFingerprintParams::default())
            .map(fingerprint_bits)
            .map_err(|error| error.to_string())
    }

    /// Returns the default MACCS fingerprint set bits.
    pub fn maccs_fingerprint(&self) -> Result<Vec<u32>, String> {
        self.inner
            .maccs_fingerprint(&cosmolkit::MaccsFingerprintParams::default())
            .map(fingerprint_bits)
            .map_err(|error| error.to_string())
    }

    /// Returns the default Avalon fingerprint set bits.
    pub fn avalon_fingerprint(&self) -> Result<Vec<u32>, String> {
        self.inner
            .avalon_fingerprint(&cosmolkit::AvalonFingerprintParams::default())
            .map(fingerprint_bits)
            .map_err(|error| error.to_string())
    }

    /// Returns the default Topological Torsion fingerprint set bits.
    pub fn topological_torsion_fingerprint(&self) -> Result<Vec<u32>, String> {
        let fingerprint = cosmolkit::topological_torsion_fingerprint(
            &self.inner,
            &cosmolkit::TopologicalTorsionFingerprintParams::default(),
        )
        .map_err(|error| error.to_string())?;
        Ok(fingerprint_bits(fingerprint))
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

    /// Returns a molecule with aromatic bonds kekulized.
    pub fn with_kekulized_bonds(&self) -> Result<Self, String> {
        self.inner
            .with_kekulized_bonds(false)
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

    /// Embeds one ETKDG-v3 conformer and returns the resulting molecule.
    pub fn with_3d_conformer(&self) -> Result<Self, String> {
        self.inner
            .with_3d_conformer()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Perceives stereochemistry using the canonical source operation.
    pub fn perceive_stereochemistry(&self) -> Result<(), String> {
        self.inner
            .perceive_stereochemistry()
            .map_err(|error| error.to_string())
    }

    /// Returns an upper bound for the number of source-defined stereoisomers.
    pub fn stereoisomer_count_with_options(
        &self,
        only_unassigned: bool,
        only_stereo_groups: bool,
        max_isomers: u32,
        unique: bool,
    ) -> Result<String, String> {
        let mut options = cosmolkit::StereoisomerOptions::default();
        options.only_unassigned = only_unassigned;
        options.only_stereo_groups = only_stereo_groups;
        options.max_isomers = max_isomers as usize;
        options.unique = unique;
        self.inner
            .stereoisomer_count(&options)
            .map(|count| count.to_string())
            .map_err(|error| error.to_string())
    }

    /// Returns the canonical molecule hash.
    pub fn hash(&self) -> Result<u64, String> {
        self.inner.hash().map_err(|error| error.to_string())
    }

    /// Serializes this molecule using the canonical binary representation.
    pub fn to_binary(&self) -> Result<Vec<u8>, String> {
        cosmolkit::mol_to_binary(&self.inner).map_err(|error| error.to_string())
    }

    /// Restores a molecule from the canonical binary representation.
    pub fn from_binary(data: &[u8]) -> Result<Self, String> {
        cosmolkit::mol_from_binary(data)
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Writes a PDB block using the selected conformer and writer flavor.
    pub fn to_pdb_block(&self, conformer_id: i32, flavor: u32) -> String {
        self.inner.to_pdb_block(conformer_id, flavor)
    }

    /// Returns the distance-geometry bounds matrix in row-major order.
    pub fn dg_bounds_matrix(&self) -> Result<Vec<f64>, String> {
        self.inner
            .dg_bounds_matrix()
            .map(|rows| rows.into_iter().flatten().collect())
            .map_err(|error| error.to_string())
    }

    /// Returns all disconnected fragments as independent molecule values.
    pub fn fragments(&self) -> Result<Vec<Self>, String> {
        self.inner
            .fragments()
            .map(|molecules| molecules.into_iter().map(|inner| Self { inner }).collect())
            .map_err(|error| error.to_string())
    }

    /// Returns the largest disconnected fragment.
    pub fn largest_fragment(&self) -> Result<Self, String> {
        self.inner
            .largest_fragment()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns the Murcko scaffold.
    pub fn murcko_scaffold(&self) -> Result<Self, String> {
        self.inner
            .murcko_scaffold()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Returns the net scaffold.
    pub fn net_scaffold(&self) -> Result<Self, String> {
        self.inner
            .net_scaffold()
            .map(|inner| Self { inner })
            .map_err(|error| error.to_string())
    }

    /// Checks whether a SMARTS query matches this molecule.
    pub fn has_substruct_match(&self, smarts: &str) -> Result<bool, String> {
        let query = cosmolkit::parse_smarts(smarts, &cosmolkit::SmartsParseParams::default())
            .map_err(|error| error.to_string())?;
        Ok(cosmolkit::has_substruct_match(&self.inner, &query))
    }

    /// Returns the number of stereoisomers under RDKit-compatible defaults.
    pub fn stereoisomer_count(&self) -> Result<String, String> {
        self.inner
            .stereoisomer_count(&cosmolkit::StereoisomerOptions::default())
            .map(|count| count.to_string())
            .map_err(|error| error.to_string())
    }
}
