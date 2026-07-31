//! COSMolKit `Molecule` bridge for the toolkit-neutral InChI source port.
//!
//! The reusable engine remains in `cosmolkit-inchi`. This module only maps the
//! core graph model to that boundary and implements its RDKit toolkit callbacks
//! with existing source-backed core chemistry and registered operations.
//!
//! The public scope is limited to four scalar APIs. Source-defined behavior is
//! checked exactly against pinned RDKit 2026.03.1 and official InChI v1.07.5.
//! The official C `NormalizeAndCompare` initial-buffer allocation-failure path
//! is undefined; Rust returns a deterministic structured allocation error for
//! that path instead of claiming an exact C result. MolBlock, SDF/V3000, IXA,
//! AuxInfo conversion, incremental INCHIGEN, version-query, and extended-polymer
//! entry points are outside this public boundary.

pub use cosmolkit_inchi::{
    InchiDiagnostic, InchiDiagnosticLevel, InchiError, InchiErrorKind, InchiReturnValues,
    InchiToInchiKeyOutput, MolToInchiKeyOutput, MolToInchiOutput,
};

use cosmolkit_inchi::{
    InchiAtom, InchiBond, InchiBondDirection, InchiBondStereo, InchiBondType, InchiChiralTag,
    InchiMolecule, InchiToMolToolkit, InchiToolkitError, MolToInchiToolkit,
};

use crate::{
    AtomId, AtomSpec, BondDirection, BondOrder, BondSpec, BondStereo, ChiralTag,
    CoordinateDimension, Element, Molecule, MoleculeBuilder, TopologyTrust, ValenceAssignment,
    ValenceModel, read_parts::MoleculeReadParts,
};

/// Exact-parity contract for one public scalar InChI API.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InchiApiParitySpec {
    pub public_api: &'static str,
    pub rust_api: &'static str,
    pub official_roots: &'static [&'static str],
    pub rdkit_version: &'static str,
    pub official_inchi_version: &'static str,
    pub source_defined_behavior: &'static str,
    pub undefined_behavior: &'static str,
}

/// Structured parity matrix for the complete public InChI surface.
pub const INCHI_API_PARITY_MATRIX: &[InchiApiParitySpec] = &[
    InchiApiParitySpec {
        public_api: "Chem.MolToInchi",
        rust_api: "mol_to_inchi",
        official_roots: &["GetINCHI", "FreeINCHI"],
        rdkit_version: "2026.03.1",
        official_inchi_version: "1.07.5",
        source_defined_behavior: "exact field and branch parity in the audited scalar boundary",
        undefined_behavior: "not reached by the audited generation boundary",
    },
    InchiApiParitySpec {
        public_api: "Chem.MolToInchiKey",
        rust_api: "mol_to_inchi_key",
        official_roots: &["GetINCHI", "FreeINCHI", "GetINCHIKeyFromINCHI"],
        rdkit_version: "2026.03.1",
        official_inchi_version: "1.07.5",
        source_defined_behavior: "exact field and branch parity in the audited scalar boundary",
        undefined_behavior: "not reached by the audited generation and key boundary",
    },
    InchiApiParitySpec {
        public_api: "InchiToInchiKey",
        rust_api: "inchi_to_inchi_key",
        official_roots: &["GetINCHIKeyFromINCHI"],
        rdkit_version: "2026.03.1",
        official_inchi_version: "1.07.5",
        source_defined_behavior: "exact field and branch parity in the audited scalar boundary",
        undefined_behavior: "none identified in the audited key boundary",
    },
    InchiApiParitySpec {
        public_api: "Chem.MolFromInchi",
        rust_api: "mol_from_inchi",
        official_roots: &["GetStructFromINCHI", "FreeStructFromINCHI"],
        rdkit_version: "2026.03.1",
        official_inchi_version: "1.07.5",
        source_defined_behavior: "exact field and branch parity in the audited scalar boundary",
        undefined_behavior: "NormalizeAndCompare initial-buffer allocation failure returns structured allocation_failed instead of executing strcpy with a null source",
    },
];

/// Result of parsing an InChI into a COSMolKit molecule.
#[derive(Clone, Debug, PartialEq)]
pub struct MolFromInchiOutput {
    pub molecule: Option<Molecule>,
    pub return_values: InchiReturnValues,
    pub diagnostics: Vec<InchiDiagnostic>,
}

fn bridge_error(
    operation: &'static str,
    kind: InchiErrorKind,
    detail: impl Into<String>,
) -> InchiError {
    InchiError {
        operation,
        kind,
        detail: detail.into(),
    }
}

fn toolkit_error(kind: &'static str, detail: impl Into<String>) -> InchiToolkitError {
    InchiToolkitError {
        kind,
        message: detail.into(),
    }
}

fn adapter_chiral_tag(tag: ChiralTag) -> InchiChiralTag {
    match tag {
        ChiralTag::Unspecified => InchiChiralTag::Unspecified,
        ChiralTag::TetrahedralCw => InchiChiralTag::TetrahedralCw,
        ChiralTag::TetrahedralCcw => InchiChiralTag::TetrahedralCcw,
        ChiralTag::Other => InchiChiralTag::Other,
        ChiralTag::Tetrahedral => InchiChiralTag::Tetrahedral,
        ChiralTag::Allene => InchiChiralTag::Allene,
        ChiralTag::SquarePlanar => InchiChiralTag::SquarePlanar,
        ChiralTag::TrigonalBipyramidal => InchiChiralTag::TrigonalBipyramidal,
        ChiralTag::Octahedral => InchiChiralTag::Octahedral,
    }
}

fn core_chiral_tag(tag: InchiChiralTag) -> ChiralTag {
    match tag {
        InchiChiralTag::Unspecified => ChiralTag::Unspecified,
        InchiChiralTag::TetrahedralCw => ChiralTag::TetrahedralCw,
        InchiChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCcw,
        InchiChiralTag::Other => ChiralTag::Other,
        InchiChiralTag::Tetrahedral => ChiralTag::Tetrahedral,
        InchiChiralTag::Allene => ChiralTag::Allene,
        InchiChiralTag::SquarePlanar => ChiralTag::SquarePlanar,
        InchiChiralTag::TrigonalBipyramidal => ChiralTag::TrigonalBipyramidal,
        InchiChiralTag::Octahedral => ChiralTag::Octahedral,
    }
}

fn adapter_bond_type(order: BondOrder) -> InchiBondType {
    match order {
        BondOrder::Null | BondOrder::Unspecified => InchiBondType::Unspecified,
        BondOrder::Single => InchiBondType::Single,
        BondOrder::Double => InchiBondType::Double,
        BondOrder::Triple => InchiBondType::Triple,
        BondOrder::Quadruple => InchiBondType::Quadruple,
        BondOrder::Quintuple => InchiBondType::Quintuple,
        BondOrder::Hextuple => InchiBondType::Hextuple,
        BondOrder::OneAndHalf => InchiBondType::OneAndAHalf,
        BondOrder::TwoAndHalf => InchiBondType::TwoAndAHalf,
        BondOrder::ThreeAndHalf => InchiBondType::ThreeAndAHalf,
        BondOrder::FourAndHalf => InchiBondType::FourAndAHalf,
        BondOrder::FiveAndHalf => InchiBondType::FiveAndAHalf,
        BondOrder::Aromatic => InchiBondType::Aromatic,
        BondOrder::Ionic => InchiBondType::Ionic,
        BondOrder::DativeOne => InchiBondType::DativeOne,
        BondOrder::Dative => InchiBondType::Dative,
        BondOrder::DativeLeft => InchiBondType::DativeL,
        BondOrder::DativeRight => InchiBondType::DativeR,
        BondOrder::Hydrogen => InchiBondType::Hydrogen,
        BondOrder::ThreeCenter => InchiBondType::ThreeCenter,
        BondOrder::Other => InchiBondType::Other,
        BondOrder::Zero => InchiBondType::Zero,
    }
}

fn core_bond_order(bond_type: InchiBondType) -> BondOrder {
    match bond_type {
        InchiBondType::Unspecified => BondOrder::Unspecified,
        InchiBondType::Single => BondOrder::Single,
        InchiBondType::Double => BondOrder::Double,
        InchiBondType::Triple => BondOrder::Triple,
        InchiBondType::Quadruple => BondOrder::Quadruple,
        InchiBondType::Quintuple => BondOrder::Quintuple,
        InchiBondType::Hextuple => BondOrder::Hextuple,
        InchiBondType::OneAndAHalf => BondOrder::OneAndHalf,
        InchiBondType::TwoAndAHalf => BondOrder::TwoAndHalf,
        InchiBondType::ThreeAndAHalf => BondOrder::ThreeAndHalf,
        InchiBondType::FourAndAHalf => BondOrder::FourAndHalf,
        InchiBondType::FiveAndAHalf => BondOrder::FiveAndHalf,
        InchiBondType::Aromatic => BondOrder::Aromatic,
        InchiBondType::Ionic => BondOrder::Ionic,
        InchiBondType::Hydrogen => BondOrder::Hydrogen,
        InchiBondType::ThreeCenter => BondOrder::ThreeCenter,
        InchiBondType::DativeOne => BondOrder::DativeOne,
        InchiBondType::Dative => BondOrder::Dative,
        InchiBondType::DativeL => BondOrder::DativeLeft,
        InchiBondType::DativeR => BondOrder::DativeRight,
        InchiBondType::Other => BondOrder::Other,
        InchiBondType::Zero => BondOrder::Zero,
    }
}

fn adapter_bond_direction(direction: BondDirection, unknown_stereo: bool) -> InchiBondDirection {
    if unknown_stereo && direction == BondDirection::None {
        return InchiBondDirection::Unknown;
    }
    match direction {
        BondDirection::None => InchiBondDirection::None,
        BondDirection::BeginWedge => InchiBondDirection::BeginWedge,
        BondDirection::BeginDash => InchiBondDirection::BeginDash,
        BondDirection::EndDownRight => InchiBondDirection::EndDownRight,
        BondDirection::EndUpRight => InchiBondDirection::EndUpRight,
        BondDirection::EitherDouble => InchiBondDirection::EitherDouble,
        BondDirection::Unknown => InchiBondDirection::Unknown,
    }
}

fn core_bond_direction(direction: InchiBondDirection) -> BondDirection {
    match direction {
        InchiBondDirection::None => BondDirection::None,
        InchiBondDirection::BeginWedge => BondDirection::BeginWedge,
        InchiBondDirection::BeginDash => BondDirection::BeginDash,
        InchiBondDirection::EndDownRight => BondDirection::EndDownRight,
        InchiBondDirection::EndUpRight => BondDirection::EndUpRight,
        InchiBondDirection::EitherDouble => BondDirection::EitherDouble,
        InchiBondDirection::Unknown => BondDirection::Unknown,
    }
}

fn adapter_bond_stereo(stereo: BondStereo) -> InchiBondStereo {
    match stereo {
        BondStereo::None => InchiBondStereo::None,
        BondStereo::Any => InchiBondStereo::Any,
        BondStereo::Z => InchiBondStereo::Z,
        BondStereo::E => InchiBondStereo::E,
        BondStereo::Cis => InchiBondStereo::Cis,
        BondStereo::Trans => InchiBondStereo::Trans,
        BondStereo::AtropCw => InchiBondStereo::AtropCw,
        BondStereo::AtropCcw => InchiBondStereo::AtropCcw,
    }
}

fn core_bond_stereo(stereo: InchiBondStereo) -> BondStereo {
    match stereo {
        InchiBondStereo::None => BondStereo::None,
        InchiBondStereo::Any => BondStereo::Any,
        InchiBondStereo::Z => BondStereo::Z,
        InchiBondStereo::E => BondStereo::E,
        InchiBondStereo::Cis => BondStereo::Cis,
        InchiBondStereo::Trans => BondStereo::Trans,
        InchiBondStereo::AtropCw => BondStereo::AtropCw,
        InchiBondStereo::AtropCcw => BondStereo::AtropCcw,
    }
}

fn validate_supported_source(read: MoleculeReadParts<'_>) -> Result<(), String> {
    if read.capabilities().topology_trust() != TopologyTrust::TrustedGraph {
        return Err("InChI requires trusted molecular graph topology".to_owned());
    }
    if !read.substance_groups().is_empty() {
        return Err("InChI bridge does not model substance groups".to_owned());
    }
    if !read.stereo_groups().is_empty() {
        return Err("InChI bridge does not model enhanced stereo groups".to_owned());
    }
    if read.atoms().iter().any(|atom| atom.query().is_some()) {
        return Err("InChI bridge does not model query atoms".to_owned());
    }
    if read.bonds().iter().any(|bond| bond.query().is_some()) {
        return Err("InChI bridge does not model query bonds".to_owned());
    }
    if read.atoms().iter().any(|atom| atom.unknown_stereo()) {
        return Err("InChI bridge cannot preserve atom-level unknown stereo state".to_owned());
    }
    if read.atoms().len() > u32::MAX as usize || read.bonds().len() > u32::MAX as usize {
        return Err("InChI adapter graph index exceeds the source u32 range".to_owned());
    }
    let coordinates = read.coordinates();
    if !coordinates.conformers_2d.is_empty() && !coordinates.conformers_3d.is_empty() {
        return Err(
            "InChI bridge cannot recover source conformer ordering from mixed 2D and 3D stores"
                .to_owned(),
        );
    }
    Ok(())
}

fn adapter_from_read_parts(
    read: MoleculeReadParts<'_>,
) -> Result<(InchiMolecule, Vec<CoordinateDimension>), String> {
    let atoms = read
        .atoms()
        .iter()
        .map(|atom| InchiAtom {
            atomic_number: i32::from(atom.atomic_number()),
            formal_charge: i32::from(atom.formal_charge()),
            num_explicit_hydrogens: u32::from(atom.explicit_hydrogens()),
            is_aromatic: atom.is_aromatic(),
            isotope: u32::from(atom.isotope().unwrap_or(0)),
            num_radical_electrons: u32::from(atom.radical_electrons()),
            no_implicit: atom.no_implicit(),
            chiral_tag: adapter_chiral_tag(atom.chiral_tag()),
            cip_rank: atom
                .prop("_CIPRank")
                .and_then(|rank| rank.parse::<u32>().ok()),
        })
        .collect::<Vec<_>>();

    let mut bonds = Vec::with_capacity(read.num_bonds());
    for bond in read.bonds() {
        let begin = u32::try_from(bond.begin().index())
            .map_err(|_| "bond begin atom index exceeds u32".to_owned())?;
        let end = u32::try_from(bond.end().index())
            .map_err(|_| "bond end atom index exceeds u32".to_owned())?;
        let mut converted = InchiBond::new(begin, end, adapter_bond_type(bond.order()));
        converted.direction = adapter_bond_direction(bond.direction(), bond.unknown_stereo());
        converted.is_aromatic = bond.is_aromatic();
        converted.stereo = adapter_bond_stereo(bond.stereo());
        converted.stereo_atoms = bond
            .stereo_atoms()
            .map(|atoms| {
                atoms
                    .into_iter()
                    .map(|atom| {
                        u32::try_from(atom.index())
                            .map_err(|_| "bond stereo atom index exceeds u32".to_owned())
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .transpose()?
            .unwrap_or_default();
        bonds.push(converted);
    }

    let coordinates = read.coordinates();
    let mut conformers = Vec::new();
    let mut coordinate_dimensions = Vec::new();
    for conformer in &coordinates.conformers_2d {
        conformers.push(
            conformer
                .coordinates()
                .iter()
                .map(|coordinate| [coordinate[0], coordinate[1], 0.0])
                .collect(),
        );
        coordinate_dimensions.push(CoordinateDimension::TwoD);
    }
    for conformer in &coordinates.conformers_3d {
        conformers.push(conformer.coordinates().to_vec());
        coordinate_dimensions.push(if conformer.is_3d() {
            CoordinateDimension::ThreeD
        } else {
            CoordinateDimension::TwoD
        });
    }

    let molecule = InchiMolecule::try_from_graph(atoms, bonds, conformers).map_err(|error| {
        format!(
            "adapter graph bond {} references atom {} outside {} atoms",
            error.bond_index, error.atom_index, error.atom_count
        )
    })?;
    Ok((molecule, coordinate_dimensions))
}

fn core_from_adapter(
    molecule: &InchiMolecule,
    coordinate_dimensions: &[CoordinateDimension],
) -> Result<Molecule, String> {
    if molecule.conformers().len() != coordinate_dimensions.len() {
        return Err("adapter conformer dimension metadata is misaligned".to_owned());
    }

    let mut builder = MoleculeBuilder::new();
    for (index, atom) in molecule.atoms().iter().enumerate() {
        let atomic_number = u8::try_from(atom.atomic_number)
            .ok()
            .filter(|number| *number <= 118)
            .ok_or_else(|| {
                format!(
                    "atom {index} has unsupported atomic number {}",
                    atom.atomic_number
                )
            })?;
        let formal_charge = i8::try_from(atom.formal_charge)
            .map_err(|_| format!("atom {index} formal charge exceeds i8"))?;
        let explicit_hydrogens = u8::try_from(atom.num_explicit_hydrogens)
            .map_err(|_| format!("atom {index} explicit hydrogen count exceeds u8"))?;
        let isotope =
            u16::try_from(atom.isotope).map_err(|_| format!("atom {index} isotope exceeds u16"))?;
        let radical_electrons = u8::try_from(atom.num_radical_electrons)
            .map_err(|_| format!("atom {index} radical electron count exceeds u8"))?;
        let element = Element::from_atomic_number(atomic_number)
            .ok_or_else(|| format!("atom {index} atomic number is unavailable"))?;
        let mut spec = AtomSpec::new(element)
            .with_formal_charge(formal_charge)
            .with_explicit_hydrogens(explicit_hydrogens)
            .with_aromatic(atom.is_aromatic)
            .with_no_implicit(atom.no_implicit)
            .with_radical_electrons(radical_electrons)
            .with_chiral_tag(core_chiral_tag(atom.chiral_tag));
        if isotope != 0 {
            spec = spec.with_isotope(isotope);
        }
        if let Some(rank) = atom.cip_rank {
            spec = spec.with_prop("_CIPRank", rank.to_string());
        }
        builder.add_atom(spec);
    }

    for (index, bond) in molecule.bonds().iter().enumerate() {
        let begin = usize::try_from(bond.begin_atom_index())
            .map_err(|_| format!("bond {index} begin atom index exceeds usize"))?;
        let end = usize::try_from(bond.end_atom_index())
            .map_err(|_| format!("bond {index} end atom index exceeds usize"))?;
        let mut spec = BondSpec::new(
            AtomId::new(begin),
            AtomId::new(end),
            core_bond_order(bond.bond_type),
        )
        .with_aromatic(bond.is_aromatic)
        .with_direction(core_bond_direction(bond.direction))
        .with_stereo(core_bond_stereo(bond.stereo))
        .with_unknown_stereo(bond.direction == InchiBondDirection::Unknown);
        match bond.stereo_atoms.as_slice() {
            [] => {}
            [left, right] => {
                let left = usize::try_from(*left)
                    .map_err(|_| format!("bond {index} left stereo atom exceeds usize"))?;
                let right = usize::try_from(*right)
                    .map_err(|_| format!("bond {index} right stereo atom exceeds usize"))?;
                spec = spec.with_stereo_atoms(AtomId::new(left), AtomId::new(right));
            }
            _ => {
                return Err(format!(
                    "bond {index} has {} stereo atoms instead of zero or two",
                    bond.stereo_atoms.len()
                ));
            }
        }
        builder
            .add_bond(spec)
            .map_err(|error| format!("bond {index} cannot be built: {error}"))?;
    }

    for (conformer, dimension) in molecule
        .conformers()
        .iter()
        .zip(coordinate_dimensions.iter().copied())
    {
        match dimension {
            CoordinateDimension::TwoD => builder
                .add_2d_conformer(
                    conformer
                        .iter()
                        .map(|coordinate| [coordinate[0], coordinate[1]])
                        .collect(),
                )
                .map_err(|error| format!("2D conformer cannot be built: {error}"))?,
            CoordinateDimension::ThreeD => builder
                .add_3d_conformer(conformer.clone())
                .map_err(|error| format!("3D conformer cannot be built: {error}"))?,
        };
    }

    builder
        .build()
        .map_err(|error| format!("adapter molecule cannot be built: {error}"))
}

struct CoreInchiToolkit {
    source_needs_property_cache_update: bool,
    coordinate_dimensions: Vec<CoordinateDimension>,
}

impl CoreInchiToolkit {
    fn for_generation(read: MoleculeReadParts<'_>, dimensions: Vec<CoordinateDimension>) -> Self {
        Self {
            source_needs_property_cache_update: read.derived_cache().valence.is_none(),
            coordinate_dimensions: dimensions,
        }
    }

    fn for_structure() -> Self {
        Self {
            source_needs_property_cache_update: true,
            coordinate_dimensions: Vec::new(),
        }
    }

    fn molecule(&self, adapter: &InchiMolecule) -> Result<Molecule, InchiToolkitError> {
        core_from_adapter(adapter, &self.coordinate_dimensions)
            .map_err(|detail| toolkit_error("COSMolKit molecule conversion", detail))
    }

    fn replace_adapter(
        &mut self,
        adapter: &mut InchiMolecule,
        molecule: &Molecule,
    ) -> Result<(), InchiToolkitError> {
        let read = MoleculeReadParts::from_molecule(molecule);
        let (replacement, dimensions) = adapter_from_read_parts(read)
            .map_err(|detail| toolkit_error("COSMolKit adapter conversion", detail))?;
        adapter
            .replace_graph(
                replacement.atoms().to_vec(),
                replacement.bonds().to_vec(),
                replacement.conformers().to_vec(),
            )
            .map_err(|error| {
                toolkit_error(
                    "COSMolKit adapter graph",
                    format!(
                        "bond {} references atom {} outside {} atoms",
                        error.bond_index, error.atom_index, error.atom_count
                    ),
                )
            })?;
        self.coordinate_dimensions = dimensions;
        Ok(())
    }

    fn valence(&self, adapter: &InchiMolecule) -> Result<ValenceAssignment, InchiToolkitError> {
        let molecule = self.molecule(adapter)?;
        MoleculeReadParts::from_molecule(&molecule)
            .assign_valence_with_options(ValenceModel::RdkitLike, false)
            .map_err(|error| toolkit_error("COSMolKit valence error", error.to_string()))
    }

    fn atom_index(
        &self,
        molecule: &InchiMolecule,
        atom_index: u32,
    ) -> Result<usize, InchiToolkitError> {
        let index = usize::try_from(atom_index).map_err(|_| {
            toolkit_error(
                "Range Error",
                format!("atom index {atom_index} exceeds usize"),
            )
        })?;
        if index >= molecule.atoms().len() {
            return Err(toolkit_error(
                "Range Error",
                format!(
                    "atom index {atom_index} is outside {} atoms",
                    molecule.atoms().len()
                ),
            ));
        }
        Ok(index)
    }
}

impl MolToInchiToolkit for CoreInchiToolkit {
    fn needs_update_property_cache(
        &mut self,
        _molecule: &InchiMolecule,
    ) -> Result<bool, InchiToolkitError> {
        Ok(self.source_needs_property_cache_update)
    }

    fn update_property_cache(
        &mut self,
        molecule: &mut InchiMolecule,
        strict: bool,
    ) -> Result<(), InchiToolkitError> {
        let core = self.molecule(molecule)?;
        let _updated = if strict {
            core.with_assigned_valence()
        } else {
            MoleculeReadParts::from_molecule(&core)
                .assign_valence_with_options(ValenceModel::RdkitLike, false)
                .map(|_| core)
                .map_err(|source| crate::OperationError::Valence {
                    operation: &crate::ASSIGNED_VALENCE_SPEC,
                    source,
                })
        }
        .map_err(|error| toolkit_error("COSMolKit property cache error", error.to_string()))?;
        self.source_needs_property_cache_update = false;
        Ok(())
    }

    fn kekulize(
        &mut self,
        molecule: &mut InchiMolecule,
        clear_aromatic_flags: bool,
    ) -> Result<(), InchiToolkitError> {
        let core = self.molecule(molecule)?;
        let updated = core
            .with_kekulized_bonds(clear_aromatic_flags)
            .map_err(|error| toolkit_error("COSMolKit kekulize error", error.to_string()))?;
        self.replace_adapter(molecule, &updated)
    }

    fn element_symbol(&mut self, atomic_number: i32) -> Result<Vec<u8>, InchiToolkitError> {
        let atomic_number = u8::try_from(atomic_number)
            .map_err(|_| toolkit_error("Pre-condition Violation", "Atomic number not found"))?;
        crate::rdkit_element_symbol(atomic_number)
            .map(|symbol| symbol.as_bytes().to_vec())
            .map_err(|error| toolkit_error("Pre-condition Violation", error.to_string()))
    }

    fn atomic_weight(&mut self, atomic_number: i32) -> Result<f64, InchiToolkitError> {
        let atomic_number = u8::try_from(atomic_number)
            .ok()
            .filter(|number| *number <= 118)
            .ok_or_else(|| toolkit_error("Pre-condition Violation", "Atomic number not found"))?;
        Ok(crate::valence::rdkit_atomic_mass(atomic_number, None))
    }

    fn total_num_hydrogens(
        &mut self,
        molecule: &InchiMolecule,
        atom_index: u32,
    ) -> Result<u32, InchiToolkitError> {
        let index = self.atom_index(molecule, atom_index)?;
        let implicit = self.valence(molecule)?.implicit_hydrogens[index].max(0) as u32;
        molecule.atoms()[index]
            .num_explicit_hydrogens
            .checked_add(implicit)
            .ok_or_else(|| {
                toolkit_error(
                    "COSMolKit integer range error",
                    "total hydrogen count overflows u32",
                )
            })
    }

    fn calc_implicit_valence(
        &mut self,
        molecule: &mut InchiMolecule,
        atom_index: u32,
    ) -> Result<i32, InchiToolkitError> {
        let index = self.atom_index(molecule, atom_index)?;
        Ok(self.valence(molecule)?.implicit_hydrogens[index])
    }

    fn total_degree(
        &mut self,
        molecule: &InchiMolecule,
        atom_index: u32,
    ) -> Result<u32, InchiToolkitError> {
        let index = self.atom_index(molecule, atom_index)?;
        let core = self.molecule(molecule)?;
        let read = MoleculeReadParts::from_molecule(&core);
        let graph_degree =
            u32::try_from(read.adjacency().neighbors_of(index).len()).map_err(|_| {
                toolkit_error("COSMolKit integer range error", "atom degree exceeds u32")
            })?;
        graph_degree
            .checked_add(self.total_num_hydrogens(molecule, atom_index)?)
            .ok_or_else(|| {
                toolkit_error(
                    "COSMolKit integer range error",
                    "total degree overflows u32",
                )
            })
    }
}

impl InchiToMolToolkit for CoreInchiToolkit {
    fn atomic_number(&mut self, element: &[u8]) -> Result<i32, InchiToolkitError> {
        let element = std::str::from_utf8(element).map_err(|_| {
            toolkit_error(
                "Post-condition Violation",
                "element symbol is not valid UTF-8",
            )
        })?;
        crate::rdkit_atomic_number_from_symbol(element)
            .map(i32::from)
            .ok_or_else(|| {
                toolkit_error(
                    "Post-condition Violation",
                    format!("Element '{element}' not found"),
                )
            })
    }

    fn average_atomic_weight(&mut self, atomic_number: i32) -> Result<f64, InchiToolkitError> {
        <Self as MolToInchiToolkit>::atomic_weight(self, atomic_number)
    }

    fn update_property_cache(
        &mut self,
        molecule: &mut InchiMolecule,
        strict: bool,
    ) -> Result<(), InchiToolkitError> {
        <Self as MolToInchiToolkit>::update_property_cache(self, molecule, strict)
    }

    fn assign_atom_cip_ranks(
        &mut self,
        molecule: &mut InchiMolecule,
    ) -> Result<Vec<u32>, InchiToolkitError> {
        let core = self.molecule(molecule)?;
        MoleculeReadParts::from_molecule(&core)
            .rank_mol_atoms()
            .map_err(|error| toolkit_error("COSMolKit canonical rank error", error.to_string()))?
            .into_iter()
            .map(|rank| {
                u32::try_from(rank).map_err(|_| {
                    toolkit_error("COSMolKit integer range error", "CIP rank exceeds u32")
                })
            })
            .collect()
    }

    fn remove_hydrogens(&mut self, molecule: &mut InchiMolecule) -> Result<(), InchiToolkitError> {
        let core = self.molecule(molecule)?;
        let updated = core.without_hydrogens().map_err(|error| {
            toolkit_error("COSMolKit remove hydrogens error", error.to_string())
        })?;
        self.replace_adapter(molecule, &updated)
    }

    fn sanitize_molecule(&mut self, molecule: &mut InchiMolecule) -> Result<(), InchiToolkitError> {
        let core = self.molecule(molecule)?;
        let updated = core
            .sanitize()
            .map_err(|error| toolkit_error("COSMolKit sanitize error", error.to_string()))?;
        self.replace_adapter(molecule, &updated)
    }

    fn assign_stereochemistry(
        &mut self,
        molecule: &mut InchiMolecule,
        clean_it: bool,
        force: bool,
    ) -> Result<(), InchiToolkitError> {
        if !clean_it || !force {
            return Err(toolkit_error(
                "Unsupported Feature",
                "the audited InchiToMol path requires cleanIt=true and force=true",
            ));
        }
        let core = self.molecule(molecule)?;
        core.perceive_stereochemistry()
            .map_err(|error| toolkit_error("COSMolKit stereochemistry error", error.to_string()))
    }
}

/// Generates an InChI from a COSMolKit molecule without mutating it.
pub fn mol_to_inchi(
    molecule: &Molecule,
    options: Option<&[u8]>,
) -> Result<MolToInchiOutput, InchiError> {
    let read = MoleculeReadParts::from_molecule(molecule);
    validate_supported_source(read)
        .map_err(|detail| bridge_error("mol_to_inchi", InchiErrorKind::UnsupportedState, detail))?;
    let (adapter, dimensions) = adapter_from_read_parts(read)
        .map_err(|detail| bridge_error("mol_to_inchi", InchiErrorKind::InvalidInput, detail))?;
    let mut toolkit = CoreInchiToolkit::for_generation(read, dimensions);
    cosmolkit_inchi::mol_to_inchi(&mut toolkit, &adapter, options)
}

/// Generates an InChIKey from a COSMolKit molecule without mutating it.
pub fn mol_to_inchi_key(
    molecule: &Molecule,
    options: Option<&[u8]>,
) -> Result<MolToInchiKeyOutput, InchiError> {
    let read = MoleculeReadParts::from_molecule(molecule);
    validate_supported_source(read).map_err(|detail| {
        bridge_error("mol_to_inchi_key", InchiErrorKind::UnsupportedState, detail)
    })?;
    let (adapter, dimensions) = adapter_from_read_parts(read)
        .map_err(|detail| bridge_error("mol_to_inchi_key", InchiErrorKind::InvalidInput, detail))?;
    let mut toolkit = CoreInchiToolkit::for_generation(read, dimensions);
    cosmolkit_inchi::mol_to_inchi_key(&mut toolkit, &adapter, options)
}

/// Generates an InChIKey directly from an InChI byte string.
pub fn inchi_to_inchi_key(inchi: &[u8]) -> Result<InchiToInchiKeyOutput, InchiError> {
    cosmolkit_inchi::inchi_to_inchi_key(inchi)
}

/// Parses an InChI into a COSMolKit molecule through the source-backed bridge.
pub fn mol_from_inchi(
    inchi: &[u8],
    sanitize: bool,
    remove_hs: bool,
) -> Result<MolFromInchiOutput, InchiError> {
    let mut toolkit = CoreInchiToolkit::for_structure();
    let output = cosmolkit_inchi::mol_from_inchi(&mut toolkit, inchi, sanitize, remove_hs)?;
    let molecule = output
        .molecule
        .as_ref()
        .map(|molecule| core_from_adapter(molecule, &toolkit.coordinate_dimensions))
        .transpose()
        .map_err(|detail| {
            bridge_error(
                "mol_from_inchi",
                InchiErrorKind::InvalidSourceOutput,
                detail,
            )
        })?;
    Ok(MolFromInchiOutput {
        molecule,
        return_values: output.return_values,
        diagnostics: output.diagnostics,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        AtomQueryPredicate, BondId, QueryNode, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
    };

    #[test]
    fn inchi_api_parity_matrix_covers_exactly_the_four_public_scalar_apis() {
        assert_eq!(crate::INCHI_FEATURE.name, "inchi.scalar");
        assert!(matches!(
            crate::INCHI_FEATURE.status,
            crate::SupportStatus::SupportedWithRdkitParity {
                rdkit_version: "2026.03.1"
            }
        ));
        assert!(
            crate::PUBLIC_FEATURES
                .iter()
                .any(|feature| feature.name == crate::INCHI_FEATURE.name)
        );
        assert_eq!(INCHI_API_PARITY_MATRIX.len(), 4);
        assert_eq!(
            INCHI_API_PARITY_MATRIX
                .iter()
                .map(|entry| entry.public_api)
                .collect::<Vec<_>>(),
            [
                "Chem.MolToInchi",
                "Chem.MolToInchiKey",
                "InchiToInchiKey",
                "Chem.MolFromInchi",
            ]
        );
        assert!(INCHI_API_PARITY_MATRIX.iter().all(|entry| {
            entry.rdkit_version == "2026.03.1"
                && entry.official_inchi_version == "1.07.5"
                && entry.source_defined_behavior.starts_with("exact")
        }));
        assert!(
            INCHI_API_PARITY_MATRIX[3]
                .undefined_behavior
                .contains("structured allocation_failed")
        );
    }

    fn methane() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        builder.build().unwrap()
    }

    #[test]
    fn inchi_core_bridge_public_scalar_apis_match_exact_source_defined_methane_results() {
        let molecule = methane();

        let inchi = mol_to_inchi(&molecule, None).unwrap();
        assert_eq!(inchi.inchi, b"InChI=1S/CH4/h1H4");
        assert_eq!(inchi.return_values.return_code, 0);
        assert!(inchi.diagnostics.is_empty());

        let molecule_key = mol_to_inchi_key(&molecule, None).unwrap();
        assert_eq!(molecule_key.key, b"VNWKTOKETHGBQD-UHFFFAOYSA-N");
        assert!(molecule_key.diagnostics.is_empty());

        let direct_key = inchi_to_inchi_key(&inchi.inchi).unwrap();
        assert_eq!(direct_key.key, molecule_key.key);
        assert!(direct_key.diagnostics.is_empty());

        let parsed = mol_from_inchi(&inchi.inchi, false, false).unwrap();
        assert_eq!(parsed.return_values.return_code, 0);
        let parsed = parsed.molecule.unwrap();
        assert_eq!(parsed.num_atoms(), 1);
        assert_eq!(parsed.num_bonds(), 0);
        assert_eq!(parsed.atoms()[0].atomic_number(), 6);
        assert_eq!(parsed.atoms()[0].explicit_hydrogens(), 4);
    }

    #[test]
    fn inchi_core_bridge_parsed_graph_preserves_isotope_charge_bond_and_stereo_fields() {
        let output =
            mol_from_inchi(b"InChI=1S/CHBrClF/c2-1(3)4/t1-/m0/s1/i1+1", false, false).unwrap();
        assert_eq!(output.return_values.return_code, 1);
        let molecule = output.molecule.unwrap();

        assert_eq!(molecule.num_atoms(), 4);
        assert_eq!(molecule.num_bonds(), 3);
        assert_eq!(molecule.atoms()[0].atomic_number(), 6);
        assert_eq!(molecule.atoms()[0].isotope(), Some(13));
        assert_eq!(molecule.atoms()[0].formal_charge(), 0);
        assert!(matches!(
            molecule.atoms()[0].chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ));
        for bond in molecule.bonds() {
            assert!(bond.begin().index() < molecule.num_atoms());
            assert!(bond.end().index() < molecule.num_atoms());
            assert_eq!(bond.order(), BondOrder::Single);
        }
        let adjacency_edges = (0..molecule.num_atoms())
            .map(|atom| molecule.topology_block().adjacency.neighbors_of(atom).len())
            .sum::<usize>();
        assert_eq!(adjacency_edges, molecule.num_bonds() * 2);
    }

    #[test]
    fn inchi_core_bridge_enum_and_coordinate_conversion_is_lossless_in_modeled_state() {
        let orders = [
            BondOrder::Null,
            BondOrder::Unspecified,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Quadruple,
            BondOrder::Quintuple,
            BondOrder::Hextuple,
            BondOrder::OneAndHalf,
            BondOrder::TwoAndHalf,
            BondOrder::ThreeAndHalf,
            BondOrder::FourAndHalf,
            BondOrder::FiveAndHalf,
            BondOrder::Aromatic,
            BondOrder::Ionic,
            BondOrder::Dative,
            BondOrder::DativeOne,
            BondOrder::DativeLeft,
            BondOrder::DativeRight,
            BondOrder::Hydrogen,
            BondOrder::ThreeCenter,
            BondOrder::Other,
            BondOrder::Zero,
        ];
        for order in orders {
            let round_trip = core_bond_order(adapter_bond_type(order));
            if order == BondOrder::Null {
                assert_eq!(round_trip, BondOrder::Unspecified);
            } else {
                assert_eq!(round_trip, order);
            }
        }

        let chiral_tags = [
            ChiralTag::Unspecified,
            ChiralTag::TetrahedralCw,
            ChiralTag::TetrahedralCcw,
            ChiralTag::Other,
            ChiralTag::Tetrahedral,
            ChiralTag::Allene,
            ChiralTag::SquarePlanar,
            ChiralTag::TrigonalBipyramidal,
            ChiralTag::Octahedral,
        ];
        for tag in chiral_tags {
            assert_eq!(core_chiral_tag(adapter_chiral_tag(tag)), tag);
        }

        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_isotope(13)
                .with_formal_charge(-1)
                .with_radical_electrons(1)
                .with_chiral_tag(ChiralTag::TetrahedralCw),
        );
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(
                BondSpec::new(carbon, oxygen, BondOrder::Double)
                    .with_direction(BondDirection::EitherDouble)
                    .with_stereo(BondStereo::Any),
            )
            .unwrap();
        builder
            .add_3d_conformer(vec![[1.0, -0.0, 2.0], [-3.0, 4.0, 5.0]])
            .unwrap();
        let molecule = builder.build().unwrap();
        let read = MoleculeReadParts::from_molecule(&molecule);
        let (adapter, dimensions) = adapter_from_read_parts(read).unwrap();
        let round_trip = core_from_adapter(&adapter, &dimensions).unwrap();

        assert_eq!(round_trip.atoms(), molecule.atoms());
        assert_eq!(round_trip.bonds(), molecule.bonds());
        assert_eq!(round_trip.conformers_3d(), molecule.conformers_3d());
        assert_eq!(
            round_trip.source_coordinate_dim(),
            Some(CoordinateDimension::ThreeD)
        );
    }

    #[test]
    fn inchi_core_bridge_registered_remove_hydrogens_remaps_coordinates_and_graph() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, -0.0], [1.0, 0.0]])
            .unwrap();
        let molecule = builder.build().unwrap();
        let read = MoleculeReadParts::from_molecule(&molecule);
        let (mut adapter, dimensions) = adapter_from_read_parts(read).unwrap();
        let mut toolkit = CoreInchiToolkit::for_generation(read, dimensions);

        InchiToMolToolkit::remove_hydrogens(&mut toolkit, &mut adapter).unwrap();

        assert_eq!(adapter.atoms().len(), 1);
        assert!(adapter.bonds().is_empty());
        assert_eq!(adapter.conformers(), &[vec![[0.0, -0.0, 0.0]]]);
        let rebuilt = core_from_adapter(&adapter, &toolkit.coordinate_dimensions).unwrap();
        assert_eq!(rebuilt.num_atoms(), 1);
        assert_eq!(rebuilt.coordinates_2d(), Some(&[[0.0, -0.0]][..]));
        assert!(rebuilt.derived_cache().valence.is_none());
        let recomputed = rebuilt.with_assigned_valence().unwrap();
        assert!(recomputed.derived_cache().valence.is_some());
    }

    #[test]
    fn inchi_core_bridge_generation_preserves_source_cow_properties_and_cache() {
        let source = methane()
            .with_name("source methane")
            .with_prop("record", "kept")
            .with_assigned_valence()
            .unwrap();
        let before = source.clone();
        assert!(source.derived_cache().valence.is_some());

        let generated = mol_to_inchi(&source, Some(b"-AuxNone")).unwrap();

        assert_eq!(generated.inchi, b"InChI=1S/CH4/h1H4");
        assert_eq!(source, before);
        assert_eq!(source.properties().name(), Some("source methane"));
        assert_eq!(source.prop("record"), Some("kept"));
        assert!(source.derived_cache().valence.is_some());
    }

    #[test]
    fn inchi_core_bridge_rejects_unmodeled_state_with_structured_unsupported_errors() {
        let unknown = Molecule::new();
        let error = mol_to_inchi(&unknown, None).unwrap_err();
        assert_eq!(error.kind, InchiErrorKind::UnsupportedState);
        assert!(error.detail.contains("trusted"));

        let mut query_builder = MoleculeBuilder::new();
        query_builder.add_atom(
            AtomSpec::new(Element::C).with_query(QueryNode::predicate(AtomQueryPredicate::Any)),
        );
        let query = query_builder.build().unwrap();
        let error = mol_to_inchi(&query, None).unwrap_err();
        assert_eq!(error.kind, InchiErrorKind::UnsupportedState);
        assert!(error.detail.contains("query atoms"));

        let mut sgroup_builder = MoleculeBuilder::new();
        sgroup_builder.add_atom(AtomSpec::new(Element::C));
        sgroup_builder
            .add_substance_group(SubstanceGroup::new(
                SubstanceGroupId::new(0),
                SubstanceGroupKind::Data,
            ))
            .unwrap();
        let sgroup = sgroup_builder.build().unwrap();
        let error = mol_to_inchi(&sgroup, None).unwrap_err();
        assert_eq!(error.kind, InchiErrorKind::UnsupportedState);
        assert!(error.detail.contains("substance groups"));

        let mut mixed_builder = MoleculeBuilder::new();
        mixed_builder.add_atom(AtomSpec::new(Element::C));
        mixed_builder.set_2d_coordinates(vec![[0.0, 0.0]]).unwrap();
        mixed_builder
            .add_3d_conformer(vec![[0.0, 0.0, 0.0]])
            .unwrap();
        let mixed = mixed_builder.build().unwrap();
        let error = mol_to_inchi(&mixed, None).unwrap_err();
        assert_eq!(error.kind, InchiErrorKind::UnsupportedState);
        assert!(error.detail.contains("mixed 2D and 3D"));
    }

    #[test]
    fn inchi_core_bridge_invalid_adapter_output_is_rejected_without_partial_molecule() {
        let invalid_charge = InchiMolecule::try_from_graph(
            vec![InchiAtom {
                atomic_number: 6,
                formal_charge: i32::from(i8::MAX) + 1,
                ..InchiAtom::default()
            }],
            Vec::new(),
            Vec::new(),
        )
        .unwrap();
        let error = core_from_adapter(&invalid_charge, &[]).unwrap_err();
        assert!(error.contains("formal charge exceeds i8"));

        let invalid_stereo = InchiMolecule::try_from_graph(
            vec![InchiAtom::default(), InchiAtom::default()],
            vec![{
                let mut bond = InchiBond::new(0, 1, InchiBondType::Double);
                bond.stereo_atoms = vec![0];
                bond
            }],
            Vec::new(),
        )
        .unwrap();
        let error = core_from_adapter(&invalid_stereo, &[]).unwrap_err();
        assert!(error.contains("instead of zero or two"));
    }

    #[test]
    fn inchi_core_bridge_preserves_structured_allocation_error_category() {
        let error = bridge_error(
            "mol_to_inchi",
            InchiErrorKind::AllocationFailed,
            "AllocationFailed",
        );
        assert_eq!(error.operation, "mol_to_inchi");
        assert_eq!(error.kind, InchiErrorKind::AllocationFailed);
        assert_eq!(error.detail, "AllocationFailed");
    }

    #[test]
    fn inchi_core_bridge_bond_index_types_remain_valid_after_adapter_round_trip() {
        let mut builder = MoleculeBuilder::new();
        let atoms = (0..4)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        for pair in atoms.windows(2) {
            builder
                .add_bond(BondSpec::new(pair[0], pair[1], BondOrder::Single))
                .unwrap();
        }
        let molecule = builder.build().unwrap();
        let (adapter, dimensions) =
            adapter_from_read_parts(MoleculeReadParts::from_molecule(&molecule)).unwrap();
        let round_trip = core_from_adapter(&adapter, &dimensions).unwrap();

        for (index, bond) in round_trip.bonds().iter().enumerate() {
            assert_eq!(bond.id(), BondId::new(index));
            assert!(bond.begin().index() < round_trip.num_atoms());
            assert!(bond.end().index() < round_trip.num_atoms());
        }
    }
}
