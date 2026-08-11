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
    structure_cache: Option<Molecule>,
}

impl CoreInchiToolkit {
    fn for_generation(read: MoleculeReadParts<'_>, dimensions: Vec<CoordinateDimension>) -> Self {
        Self {
            source_needs_property_cache_update: read.derived_cache().valence.is_none(),
            coordinate_dimensions: dimensions,
            structure_cache: None,
        }
    }

    fn for_structure() -> Self {
        Self {
            source_needs_property_cache_update: true,
            coordinate_dimensions: Vec::new(),
            structure_cache: None,
        }
    }

    fn molecule(&self, adapter: &InchiMolecule) -> Result<Molecule, InchiToolkitError> {
        core_from_adapter(adapter, &self.coordinate_dimensions)
            .map_err(|detail| toolkit_error("COSMolKit molecule conversion", detail))
    }

    fn cache_structure(&mut self, adapter: &InchiMolecule) -> Result<(), InchiToolkitError> {
        self.structure_cache = Some(self.molecule(adapter)?);
        Ok(())
    }

    fn take_cached_or_build(
        &mut self,
        adapter: &InchiMolecule,
    ) -> Result<Molecule, InchiToolkitError> {
        self.structure_cache
            .take()
            .map_or_else(|| self.molecule(adapter), Ok)
    }

    fn finish_structure(&mut self, adapter: &InchiMolecule) -> Result<Molecule, InchiToolkitError> {
        self.take_cached_or_build(adapter)
    }

    fn replace_adapter(
        &mut self,
        adapter: &mut InchiMolecule,
        molecule: Molecule,
    ) -> Result<(), InchiToolkitError> {
        // BEGIN RDKIT C++ SOURCE: External/INCHI-API/inchi.cpp:1651-1668
        // RDKit✔️❌:     cleanUp(*m);
        // RDKit✔️✔️:     if (sanitize) {
        // RDKit✔️✔️:       if (removeHs) {
        // RDKit✔️✔️:         MolOps::removeHs(*m);
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         MolOps::sanitizeMol(*m);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     MolOps::assignStereochemistry(*m, true, true);
        // The source mutates one RWMol through this sequence. Moving the owned
        // Molecule into the cache preserves that identity across the adapter
        // boundary and avoids cloning the completed graph. The marker remains
        // performance-negative because the adapter materialization itself has
        // no RDKit counterpart.
        let read = MoleculeReadParts::from_molecule(&molecule);
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
        self.structure_cache = Some(molecule);
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
        let mut core = self.take_cached_or_build(molecule)?;
        core.assign_valence_strict_(strict)
            .map_err(|error| toolkit_error("COSMolKit property cache error", error.to_string()))?;
        self.structure_cache = Some(core);
        self.source_needs_property_cache_update = false;
        Ok(())
    }

    fn kekulize(
        &mut self,
        molecule: &mut InchiMolecule,
        clear_aromatic_flags: bool,
    ) -> Result<(), InchiToolkitError> {
        let mut core = self.take_cached_or_build(molecule)?;
        core.kekulize_(clear_aromatic_flags)
            .map_err(|error| toolkit_error("COSMolKit kekulize error", error.to_string()))?;
        self.replace_adapter(molecule, core)
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
        let core = self.take_cached_or_build(molecule)?;
        let ranks = crate::stereo::assign_atom_cip_ranks(&core)
            .map_err(|error| toolkit_error("COSMolKit CIP rank error", error.to_string()))?;
        self.structure_cache = Some(core);
        Ok(ranks)
    }

    fn remove_hydrogens(&mut self, molecule: &mut InchiMolecule) -> Result<(), InchiToolkitError> {
        let mut core = self.take_cached_or_build(molecule)?;
        core.remove_hydrogens_with_sanitize_(true)
            .map_err(|error| match error {
                crate::OperationError::Sanitize { .. } => {
                    toolkit_error("MolSanitizeException", error.to_string())
                }
                _ => toolkit_error("COSMolKit remove hydrogens error", error.to_string()),
            })?;
        self.replace_adapter(molecule, core)
    }

    fn sanitize_molecule(&mut self, molecule: &mut InchiMolecule) -> Result<(), InchiToolkitError> {
        let mut core = self.take_cached_or_build(molecule)?;
        core.sanitize_().map_err(|error| match error {
            crate::OperationError::Sanitize { .. } => {
                toolkit_error("MolSanitizeException", error.to_string())
            }
            _ => toolkit_error("COSMolKit sanitize error", error.to_string()),
        })?;
        self.replace_adapter(molecule, core)
    }

    fn synchronize_after_cleanup(
        &mut self,
        molecule: &InchiMolecule,
    ) -> Result<(), InchiToolkitError> {
        // RDKit's `cleanUp`, `removeHs`/`sanitizeMol`, and
        // `assignStereochemistry` all mutate the same RWMol. The adapter cleanup
        // runs in the source crate, so this explicit boundary moves its exact
        // graph state into the corresponding native molecule once.
        self.cache_structure(molecule)
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
        let core = self.take_cached_or_build(molecule)?;
        core.perceive_stereochemistry()
            .map_err(|error| toolkit_error("COSMolKit stereochemistry error", error.to_string()))?;
        self.structure_cache = Some(core);
        Ok(())
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
        .map(|molecule| toolkit.finish_structure(molecule))
        .transpose()
        .map_err(|error| {
            bridge_error(
                "mol_from_inchi",
                InchiErrorKind::InvalidSourceOutput,
                error.message,
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
    fn inchi_generation_matches_chematic_issue_11_pointer_regressions() {
        let cases = [
            (
                440,
                "C[C@H](NC(=O)[C@@H](CC(N)=O)N(O)C(=O)[C@H](CO)NC(=O)CNC(=O)C1CCN=C(C(CO)NC(=O)[C@@H]2CCNc3c(NC(=O)CC(O)C(N)=O)cc4cc(O)c(O)cc4[n+]32)N1)C(=O)NCC(=O)N[C@H](C)C(=O)NCC(=O)N[C@H]1CCCN(O)C1=O",
                "InChI=1S/C48H67N17O20/c1-20(42(77)54-16-37(74)58-24-4-3-9-63(84)47(24)82)56-36(73)15-53-43(78)21(2)57-46(81)30(13-34(49)71)65(85)48(83)27(19-67)60-38(75)17-55-44(79)23-5-7-51-40(61-23)26(18-66)62-45(80)28-6-8-52-41-25(59-35(72)14-33(70)39(50)76)10-22-11-31(68)32(69)12-29(22)64(28)41/h10-12,20-21,23-24,26-28,30,33,66-67,70,84-85H,3-9,13-19H2,1-2H3,(H16,49,50,51,52,53,54,55,56,57,58,59,60,61,62,68,69,71,72,73,74,75,76,77,78,79,80,81)/p+1/t20-,21+,23?,24+,26?,27+,28+,30-,33?/m1/s1",
                "PZXWNCGKYBGCFQ-MLNIRCADSA-O",
            ),
            (
                904,
                "CC[C@H](C)[C@H]1NC(=O)c2cc3ccc(O)cc3[n+](C)c21",
                "InChI=1S/C16H18N2O2/c1-4-9(2)14-15-12(16(20)17-14)7-10-5-6-11(19)8-13(10)18(15)3/h5-9,14H,4H2,1-3H3,(H,17,20)/p+1/t9-,14+/m0/s1",
                "GTNBOKIWQUJZFH-LKFCYVNXSA-O",
            ),
            (
                2040,
                "COC(=O)c1cc2c[n+](C)c3c(N)c(Cl)c(=N)c([nH]1)c23",
                "InChI=1S/C13H11ClN4O2/c1-18-4-5-3-6(13(19)20-2)17-11-7(5)12(18)10(16)8(14)9(11)15/h3-4H,1-2H3,(H3,15,16)/p+1",
                "WVQCUPZVEGKBBJ-UHFFFAOYSA-O",
            ),
            (
                2556,
                "CC1=[N+](C)[C@H]2[C@@H](OC1=C(N)C(=O)O)O[C@H](CO)[C@H](O)[C@@H]2O",
                "InChI=1S/C12H18N2O7/c1-4-10(6(13)11(18)19)21-12-7(14(4)2)9(17)8(16)5(3-15)20-12/h5,7-9,12-13,15-17H,3H2,1-2H3,(H,18,19)/p+1/t5-,7-,8+,9-,12-/m1/s1",
                "WSWRZDSZMARMQJ-GXIXWSSHSA-O",
            ),
            (
                3248,
                "CC[n+]1c(-c2ccccc2)c2cc(N)ccc2c2ccc(N)cc21.[Br-]",
                "InChI=1S/C21H19N3.BrH/c1-2-24-20-13-16(23)9-11-18(20)17-10-8-15(22)12-19(17)21(24)14-6-4-3-5-7-14;/h3-13,23H,2,22H2,1H3;1H",
                "ZMMJGEGLRURXTF-UHFFFAOYSA-N",
            ),
            (
                3620,
                "Cc1cc[n+]2c(NC(=O)c3ccccc3)c(-c3cccs3)[nH]c2c1",
                "InChI=1S/C19H15N3OS/c1-13-9-10-22-16(12-13)20-17(15-8-5-11-24-15)18(22)21-19(23)14-6-3-2-4-7-14/h2-12H,1H3,(H,21,23)/p+1",
                "OVFLFOQKWYSFNY-UHFFFAOYSA-O",
            ),
            (
                3811,
                "Cc1cc2c(s1)Nc1ccccc1NC2=[N+]1CCN(C)CC1",
                "InChI=1S/C17H20N4S/c1-12-11-13-16(21-9-7-20(2)8-10-21)18-14-5-3-4-6-15(14)19-17(13)22-12/h3-6,11H,7-10H2,1-2H3,(H,18,19)/p+1",
                "HZFBJGGXPNEAPR-UHFFFAOYSA-O",
            ),
            (
                4944,
                "Cc1cc[n+]2c(-c3ccnc(N)n3)c(C)[nH]c2c1",
                "InChI=1S/C13H13N5/c1-8-4-6-18-11(7-8)16-9(2)12(18)10-3-5-15-13(14)17-10/h3-7H,1-2H3,(H2,14,15,17)/p+1",
                "OQODWFCIKXTEHB-UHFFFAOYSA-O",
            ),
        ];

        for (row, smiles, expected_inchi, expected_key) in cases {
            let molecule = Molecule::from_smiles(smiles)
                .unwrap_or_else(|error| panic!("row {row} SMILES: {error}"));
            let output = mol_to_inchi(&molecule, None)
                .unwrap_or_else(|error| panic!("row {row} InChI generation: {error}"));
            assert_eq!(output.inchi, expected_inchi.as_bytes(), "row {row} InChI");
            assert_eq!(output.return_values.return_code, 1, "row {row} return code");

            let key = inchi_to_inchi_key(&output.inchi)
                .unwrap_or_else(|error| panic!("row {row} InChIKey generation: {error}"));
            assert_eq!(key.key, expected_key.as_bytes(), "row {row} InChIKey");
        }
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
    fn inchi_core_bridge_uses_legacy_cip_ranks_for_double_bond_stereo_neighbors() {
        let output = mol_from_inchi(
            b"InChI=1S/C19H21NO3S2/c1-14(12-15-8-4-2-5-9-15)13-16-18(23)20(19(24)25-16)11-7-3-6-10-17(21)22/h2,4-5,8-9,12-13H,3,6-7,10-11H2,1H3,(H,21,22)/b14-12+,16-13-",
            true,
            true,
        )
        .expect("source-defined double-bond stereo InChI must parse");
        let molecule = output
            .molecule
            .expect("successful parse must return a graph");
        let bond = molecule
            .bonds()
            .iter()
            .find(|bond| {
                let endpoints = [bond.begin().index(), bond.end().index()];
                endpoints == [12, 15] || endpoints == [15, 12]
            })
            .expect("audited double bond must be present");

        assert_eq!(bond.stereo(), BondStereo::Z);
        assert_eq!(
            bond.stereo_atoms(),
            Some([AtomId::new(13), AtomId::new(24)])
        );
        assert_eq!(
            molecule.to_smiles(true).unwrap(),
            "CC(/C=C1\\SC(=S)N(CCCCCC(=O)O)C1=O)=C\\c1ccccc1"
        );
    }

    #[test]
    fn inchi_core_bridge_classifies_rdkit_mol_sanitize_exception() {
        let error = mol_from_inchi(
            b"InChI=1S/C8H16O6S2/c9-5-8(14-16(11,12)13)7(10)6-15-3-1-2-4-15/h7-10H,1-6H2/t7-,8+/m0/s1",
            true,
            true,
        )
        .expect_err("RDKit source path must reject the invalid sanitized valence");

        assert_eq!(error.kind, InchiErrorKind::SanitizeFailed);
    }

    #[test]
    fn inchi_core_bridge_accepts_zero_z_v2000_2d_input_without_duplicate_conformer() {
        let input = concat!(
            "zero-z\n",
            "  PubChem          2D\n",
            "\n",
            "  2  1  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "  1  2  1  0\n",
            "M  END\n",
            "$$$$\n",
        );
        let record = crate::io::sdf::read_sdf_from_str(input).unwrap();

        assert_eq!(record.molecule.conformers_2d().len(), 1);
        assert!(record.molecule.conformers_3d().is_empty());
        let output = mol_to_inchi(&record.molecule, None).unwrap();
        assert_eq!(output.inchi, b"InChI=1S/C2H6/c1-2/h1-2H3");
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
