use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
};

use cosmolkit_types::{ChiralTag, Element, Hybridization};

/// Stable atom-table index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomId(usize);

impl AtomId {
    #[must_use]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

impl fmt::Display for AtomId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}", self.0)
    }
}

/// An atom-local string property could not be stored.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AtomPropertyError {
    #[error("atom property key cannot be empty")]
    EmptyKey,
}

/// One ordered template-attachment entry carried by an atom.
///
/// `target` is a canonical COSMolKit atom-table id. Source row numbers and
/// bookmarks must be resolved before this value is constructed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TemplateAttachment {
    target: AtomId,
    label: String,
}

impl TemplateAttachment {
    #[must_use]
    pub fn new(target: AtomId, label: impl Into<String>) -> Self {
        Self {
            target,
            label: label.into(),
        }
    }

    #[must_use]
    pub const fn target(&self) -> AtomId {
        self.target
    }

    #[must_use]
    pub fn label(&self) -> &str {
        &self.label
    }
}

/// Ordered template attachment state associated with one carrier atom.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TemplateAttachmentOrder {
    entries: Vec<TemplateAttachment>,
}

/// A template attachment order violates its local structural invariants.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum TemplateAttachmentOrderError {
    #[error("template attachment order must contain at least one entry")]
    Empty,
    #[error(
        "template attachment entry {duplicate_position} repeats target {target} from entry {first_position}"
    )]
    DuplicateTarget {
        first_position: usize,
        duplicate_position: usize,
        target: AtomId,
    },
    #[error(
        "template attachment entry {duplicate_position} repeats label {label:?} from entry {first_position}"
    )]
    DuplicateLabel {
        first_position: usize,
        duplicate_position: usize,
        label: String,
    },
    #[error(
        "template attachment entry {position} references atom {target}, out of range for {atom_count} atoms"
    )]
    TargetOutOfRange {
        position: usize,
        target: AtomId,
        atom_count: usize,
    },
    #[error(
        "template attachment entry {position} references atom {target}, outside a mapping of length {mapping_len}"
    )]
    MappingTargetOutOfRange {
        position: usize,
        target: AtomId,
        mapping_len: usize,
    },
    #[error("template attachment entry {position} loses referenced atom {target}")]
    TargetRemoved { position: usize, target: AtomId },
}

impl TemplateAttachmentOrder {
    pub fn new(entries: Vec<TemplateAttachment>) -> Result<Self, TemplateAttachmentOrderError> {
        if entries.is_empty() {
            return Err(TemplateAttachmentOrderError::Empty);
        }
        let mut targets = BTreeMap::new();
        let mut labels = BTreeMap::new();
        for (position, entry) in entries.iter().enumerate() {
            if let Some(first_position) = targets.insert(entry.target(), position) {
                return Err(TemplateAttachmentOrderError::DuplicateTarget {
                    first_position,
                    duplicate_position: position,
                    target: entry.target(),
                });
            }
            if let Some(first_position) = labels.insert(entry.label(), position) {
                return Err(TemplateAttachmentOrderError::DuplicateLabel {
                    first_position,
                    duplicate_position: position,
                    label: entry.label().to_owned(),
                });
            }
        }
        Ok(Self { entries })
    }

    #[must_use]
    pub fn entries(&self) -> &[TemplateAttachment] {
        &self.entries
    }

    pub fn validate_for_atom_count(
        &self,
        atom_count: usize,
    ) -> Result<(), TemplateAttachmentOrderError> {
        for (position, entry) in self.entries.iter().enumerate() {
            if entry.target().index() >= atom_count {
                return Err(TemplateAttachmentOrderError::TargetOutOfRange {
                    position,
                    target: entry.target(),
                    atom_count,
                });
            }
        }
        Ok(())
    }

    /// Remap every referenced atom through one validated old-to-new index map.
    #[doc(hidden)]
    pub fn remapped(
        &self,
        old_to_new: &[Option<AtomId>],
    ) -> Result<Self, TemplateAttachmentOrderError> {
        let mut entries = Vec::with_capacity(self.entries.len());
        for (position, entry) in self.entries.iter().enumerate() {
            let Some(mapped) = old_to_new.get(entry.target().index()) else {
                return Err(TemplateAttachmentOrderError::MappingTargetOutOfRange {
                    position,
                    target: entry.target(),
                    mapping_len: old_to_new.len(),
                });
            };
            let Some(target) = *mapped else {
                return Err(TemplateAttachmentOrderError::TargetRemoved {
                    position,
                    target: entry.target(),
                });
            };
            entries.push(TemplateAttachment::new(target, entry.label()));
        }
        Self::new(entries)
    }
}

/// Typed atom-level PDB residue metadata.
///
/// This models the RDKit `AtomPDBResidueInfo` subset needed by hydrogen
/// addition. It is atom state, not generic string props, so topology operations
/// remap it with the atom row.
#[derive(Debug, Clone)]
pub struct AtomPdbResidueInfo {
    atom_name: String,
    serial_number: i32,
    alt_loc: String,
    residue_name: String,
    residue_number: i32,
    chain_id: String,
    insertion_code: String,
    occupancy: f64,
    temp_factor: f64,
    is_hetero_atom: bool,
    secondary_structure: u32,
    segment_number: u32,
    monomer_class: String,
}

impl PartialEq for AtomPdbResidueInfo {
    fn eq(&self, other: &Self) -> bool {
        self.atom_name == other.atom_name
            && self.serial_number == other.serial_number
            && self.alt_loc == other.alt_loc
            && self.residue_name == other.residue_name
            && self.residue_number == other.residue_number
            && self.chain_id == other.chain_id
            && self.insertion_code == other.insertion_code
            && self.occupancy.to_bits() == other.occupancy.to_bits()
            && self.temp_factor.to_bits() == other.temp_factor.to_bits()
            && self.is_hetero_atom == other.is_hetero_atom
            && self.secondary_structure == other.secondary_structure
            && self.segment_number == other.segment_number
            && self.monomer_class == other.monomer_class
    }
}

impl Eq for AtomPdbResidueInfo {}

impl Default for AtomPdbResidueInfo {
    fn default() -> Self {
        // RDKit✔️✔️: AtomPDBResidueInfo() : AtomMonomerInfo(PDBRESIDUE) {}
        Self::new("", 0, "", 0, "", false)
    }
}

impl AtomPdbResidueInfo {
    #[must_use]
    pub fn new(
        atom_name: impl Into<String>,
        serial_number: i32,
        residue_name: impl Into<String>,
        residue_number: i32,
        chain_id: impl Into<String>,
        is_hetero_atom: bool,
    ) -> Self {
        // BEGIN RDKIT CPP FUNCTION AtomPDBResidueInfo::AtomPDBResidueInfo
        // RDKit✔️✔️: AtomPDBResidueInfo(const std::string &atomName, int serialNumber = 0,
        // RDKit✔️✔️:                    std::string altLoc = "", std::string residueName = "",
        // RDKit✔️✔️:                    int residueNumber = 0, std::string chainId = "",
        // RDKit✔️✔️:                    std::string insertionCode = "", double occupancy = 1.0,
        // RDKit✔️✔️:                    double tempFactor = 0.0, bool isHeteroAtom = false,
        // RDKit✔️✔️:                    unsigned int secondaryStructure = 0,
        // RDKit✔️✔️:                    unsigned int segmentNumber = 0,
        // RDKit✔️✔️:                    std::string monomerClass = "")
        // RDKit✔️✔️:     : AtomMonomerInfo(PDBRESIDUE, atomName, residueName, residueNumber, chainId,
        // RDKit✔️✔️:                       monomerClass),
        // RDKit✔️✔️:       d_serialNumber(serialNumber),
        // RDKit✔️✔️:       d_altLoc(std::move(altLoc)),
        // RDKit✔️✔️:       d_insertionCode(std::move(insertionCode)),
        // RDKit✔️✔️:       d_occupancy(occupancy),
        // RDKit✔️✔️:       d_tempFactor(tempFactor),
        // RDKit✔️✔️:       df_heteroAtom(isHeteroAtom),
        // RDKit✔️✔️:       d_secondaryStructure(secondaryStructure),
        // RDKit✔️✔️:       d_segmentNumber(segmentNumber) {}
        // END RDKIT CPP FUNCTION AtomPDBResidueInfo::AtomPDBResidueInfo
        Self {
            atom_name: atom_name.into(),
            serial_number,
            alt_loc: String::new(),
            residue_name: residue_name.into(),
            residue_number,
            chain_id: chain_id.into(),
            insertion_code: String::new(),
            occupancy: 1.0,
            temp_factor: 0.0,
            is_hetero_atom,
            secondary_structure: 0,
            segment_number: 0,
            monomer_class: String::new(),
        }
    }

    #[must_use]
    pub fn with_alt_loc(mut self, alt_loc: impl Into<String>) -> Self {
        // RDKit✔️✔️: void setAltLoc(const std::string &val) { d_altLoc = val; }
        self.alt_loc = alt_loc.into();
        self
    }

    #[must_use]
    pub fn with_insertion_code(mut self, insertion_code: impl Into<String>) -> Self {
        // RDKit✔️✔️: void setInsertionCode(const std::string &val) { d_insertionCode = val; }
        self.insertion_code = insertion_code.into();
        self
    }

    #[must_use]
    pub const fn with_occupancy(mut self, occupancy: f64) -> Self {
        // RDKit✔️✔️: void setOccupancy(double val) { d_occupancy = val; }
        self.occupancy = occupancy;
        self
    }

    #[must_use]
    pub const fn with_temp_factor(mut self, temp_factor: f64) -> Self {
        // RDKit✔️✔️: void setTempFactor(double val) { d_tempFactor = val; }
        self.temp_factor = temp_factor;
        self
    }

    #[must_use]
    pub const fn with_secondary_structure(mut self, secondary_structure: u32) -> Self {
        // RDKit✔️✔️: void setSecondaryStructure(unsigned int val) { d_secondaryStructure = val; }
        self.secondary_structure = secondary_structure;
        self
    }

    #[must_use]
    pub const fn with_segment_number(mut self, segment_number: u32) -> Self {
        // RDKit✔️✔️: void setSegmentNumber(unsigned int val) { d_segmentNumber = val; }
        self.segment_number = segment_number;
        self
    }

    #[must_use]
    pub fn with_monomer_class(mut self, monomer_class: impl Into<String>) -> Self {
        // RDKit✔️✔️: void setMonomerClass(const std::string &val) { d_monomerClass = val; }
        self.monomer_class = monomer_class.into();
        self
    }

    #[must_use]
    pub fn atom_name(&self) -> &str {
        // RDKit✔️✔️: const std::string &getName() const { return d_name; }
        &self.atom_name
    }

    #[must_use]
    pub const fn serial_number(&self) -> i32 {
        // RDKit✔️✔️: int getSerialNumber() const { return d_serialNumber; }
        self.serial_number
    }

    #[must_use]
    pub fn alt_loc(&self) -> &str {
        // RDKit✔️✔️: const std::string &getAltLoc() const { return d_altLoc; }
        &self.alt_loc
    }

    #[must_use]
    pub fn residue_name(&self) -> &str {
        // RDKit✔️✔️: const std::string &getResidueName() const { return d_residueName; }
        &self.residue_name
    }

    #[must_use]
    pub const fn residue_number(&self) -> i32 {
        // RDKit✔️✔️: int getResidueNumber() const { return d_residueNumber; }
        self.residue_number
    }

    #[must_use]
    pub fn chain_id(&self) -> &str {
        // RDKit✔️✔️: const std::string &getChainId() const { return d_chainId; }
        &self.chain_id
    }

    #[must_use]
    pub fn insertion_code(&self) -> &str {
        // RDKit✔️✔️: const std::string &getInsertionCode() const { return d_insertionCode; }
        &self.insertion_code
    }

    #[must_use]
    pub const fn occupancy(&self) -> f64 {
        // RDKit✔️✔️: double getOccupancy() const { return d_occupancy; }
        self.occupancy
    }

    #[must_use]
    pub const fn temp_factor(&self) -> f64 {
        // RDKit✔️✔️: double getTempFactor() const { return d_tempFactor; }
        self.temp_factor
    }

    #[must_use]
    pub const fn is_hetero_atom(&self) -> bool {
        // RDKit✔️✔️: bool getIsHeteroAtom() const { return df_heteroAtom; }
        self.is_hetero_atom
    }

    #[must_use]
    pub const fn secondary_structure(&self) -> u32 {
        // RDKit✔️✔️: unsigned int getSecondaryStructure() const { return d_secondaryStructure; }
        self.secondary_structure
    }

    #[must_use]
    pub const fn segment_number(&self) -> u32 {
        // RDKit✔️✔️: unsigned int getSegmentNumber() const { return d_segmentNumber; }
        self.segment_number
    }

    #[must_use]
    pub fn monomer_class(&self) -> &str {
        // RDKit✔️✔️: const std::string &getMonomerClass() const { return d_monomerClass; }
        &self.monomer_class
    }
}

/// Atom construction payload.
///
/// `AtomSpec` is deliberately separate from `Atom`: callers provide facts, and
/// builders assign indices. Future agents must not add an `index` field here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomSpec {
    element: Element,
    formal_charge: i8,
    explicit_hydrogens: u8,
    chiral_tag: ChiralTag,
    chiral_permutation: Option<u32>,
    unknown_stereo: bool,
    mol_parity: Option<i32>,
    mol_inversion_flag: Option<i32>,
    implicit_hydrogen: bool,
    tracked_isotopic_hydrogens: Vec<u16>,
    is_aromatic: bool,
    isotope: Option<u16>,
    atom_map: Option<u32>,
    no_implicit: bool,
    radical_electrons: u8,
    hybridization: Hybridization,
    props: BTreeMap<String, String>,
    computed_props: BTreeSet<String>,
    pdb_residue_info: Option<AtomPdbResidueInfo>,
    template_attachment_order: Option<TemplateAttachmentOrder>,
}

impl AtomSpec {
    #[must_use]
    pub const fn new(element: Element) -> Self {
        Self {
            element,
            formal_charge: 0,
            explicit_hydrogens: 0,
            chiral_tag: ChiralTag::Unspecified,
            chiral_permutation: None,
            unknown_stereo: false,
            mol_parity: None,
            mol_inversion_flag: None,
            implicit_hydrogen: false,
            tracked_isotopic_hydrogens: Vec::new(),
            is_aromatic: false,
            isotope: None,
            atom_map: None,
            no_implicit: false,
            radical_electrons: 0,
            hybridization: Hybridization::Unspecified,
            props: BTreeMap::new(),
            computed_props: BTreeSet::new(),
            pdb_residue_info: None,
            template_attachment_order: None,
        }
    }

    #[must_use]
    pub const fn with_element(mut self, element: Element) -> Self {
        self.element = element;
        self
    }

    #[must_use]
    pub const fn with_formal_charge(mut self, formal_charge: i8) -> Self {
        self.formal_charge = formal_charge;
        self
    }

    #[must_use]
    pub const fn with_explicit_hydrogens(mut self, explicit_hydrogens: u8) -> Self {
        self.explicit_hydrogens = explicit_hydrogens;
        self
    }

    #[must_use]
    pub const fn with_chiral_tag(mut self, chiral_tag: ChiralTag) -> Self {
        self.chiral_tag = chiral_tag;
        self
    }

    #[must_use]
    pub const fn with_chiral_permutation(mut self, chiral_permutation: u32) -> Self {
        self.chiral_permutation = Some(chiral_permutation);
        self
    }

    #[must_use]
    pub const fn without_chiral_permutation(mut self) -> Self {
        self.chiral_permutation = None;
        self
    }

    #[must_use]
    pub const fn with_unknown_stereo(mut self, unknown_stereo: bool) -> Self {
        self.unknown_stereo = unknown_stereo;
        self
    }

    #[must_use]
    pub const fn with_mol_parity(mut self, mol_parity: i32) -> Self {
        self.mol_parity = Some(mol_parity);
        self
    }

    #[must_use]
    pub const fn without_mol_parity(mut self) -> Self {
        self.mol_parity = None;
        self
    }

    #[must_use]
    pub const fn with_mol_inversion_flag(mut self, mol_inversion_flag: i32) -> Self {
        self.mol_inversion_flag = Some(mol_inversion_flag);
        self
    }

    #[must_use]
    pub const fn without_mol_inversion_flag(mut self) -> Self {
        self.mol_inversion_flag = None;
        self
    }

    #[must_use]
    pub const fn with_implicit_hydrogen(mut self, implicit_hydrogen: bool) -> Self {
        self.implicit_hydrogen = implicit_hydrogen;
        self
    }

    #[must_use]
    pub fn with_tracked_isotopic_hydrogens(mut self, isotopes: Vec<u16>) -> Self {
        self.tracked_isotopic_hydrogens = isotopes;
        self
    }

    #[must_use]
    pub fn without_tracked_isotopic_hydrogens(mut self) -> Self {
        self.tracked_isotopic_hydrogens.clear();
        self
    }

    #[must_use]
    pub const fn with_aromatic(mut self, is_aromatic: bool) -> Self {
        self.is_aromatic = is_aromatic;
        self
    }

    #[must_use]
    pub const fn with_isotope(mut self, isotope: u16) -> Self {
        self.isotope = Some(isotope);
        self
    }

    #[must_use]
    pub const fn without_isotope(mut self) -> Self {
        self.isotope = None;
        self
    }

    #[must_use]
    pub const fn with_atom_map(mut self, atom_map: u32) -> Self {
        self.atom_map = Some(atom_map);
        self
    }

    #[must_use]
    pub const fn without_atom_map(mut self) -> Self {
        self.atom_map = None;
        self
    }

    #[must_use]
    pub const fn with_no_implicit(mut self, no_implicit: bool) -> Self {
        self.no_implicit = no_implicit;
        self
    }

    #[must_use]
    pub const fn with_radical_electrons(mut self, radical_electrons: u8) -> Self {
        self.radical_electrons = radical_electrons;
        self
    }

    #[must_use]
    pub const fn with_hybridization(mut self, hybridization: Hybridization) -> Self {
        self.hybridization = hybridization;
        self
    }

    #[must_use]
    pub fn with_prop(
        mut self,
        key: impl Into<String>,
        value: impl Into<String>,
    ) -> Result<Self, AtomPropertyError> {
        let key = key.into();
        validate_property_key(&key)?;
        // RDKit✔️✔️: d_props.setVal(key, val);
        self.props.insert(key, value.into());
        Ok(self)
    }

    #[must_use]
    pub fn with_computed_prop(
        mut self,
        key: impl Into<String>,
        value: impl Into<String>,
    ) -> Result<Self, AtomPropertyError> {
        let key = key.into();
        validate_property_key(&key)?;
        // RDKit✔️✔️: if (computed) {
        // RDKit✔️✔️:   STR_VECT compLst;
        // RDKit✔️✔️:   getPropIfPresent(RDKit::detail::computedPropName, compLst);
        // RDKit✔️✔️:   if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
        // RDKit✔️✔️:     compLst.emplace_back(key);
        // RDKit✔️✔️:     d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: d_props.setVal(key, val);
        self.props.insert(key.clone(), value.into());
        self.computed_props.insert(key);
        Ok(self)
    }

    #[must_use]
    pub fn with_pdb_residue_info(mut self, info: AtomPdbResidueInfo) -> Self {
        self.pdb_residue_info = Some(info);
        self
    }

    #[must_use]
    pub fn without_pdb_residue_info(mut self) -> Self {
        self.pdb_residue_info = None;
        self
    }

    #[must_use]
    pub fn with_template_attachment_order(mut self, order: TemplateAttachmentOrder) -> Self {
        self.template_attachment_order = Some(order);
        self
    }

    #[must_use]
    pub fn without_template_attachment_order(mut self) -> Self {
        self.template_attachment_order = None;
        self
    }

    #[must_use]
    pub const fn element(&self) -> Element {
        self.element
    }

    #[must_use]
    pub const fn formal_charge(&self) -> i8 {
        self.formal_charge
    }

    #[must_use]
    pub const fn explicit_hydrogens(&self) -> u8 {
        self.explicit_hydrogens
    }

    #[must_use]
    pub const fn chiral_tag(&self) -> ChiralTag {
        self.chiral_tag
    }

    #[must_use]
    pub const fn chiral_permutation(&self) -> Option<u32> {
        self.chiral_permutation
    }

    #[must_use]
    pub const fn unknown_stereo(&self) -> bool {
        self.unknown_stereo
    }

    #[must_use]
    pub const fn mol_parity(&self) -> Option<i32> {
        self.mol_parity
    }

    #[must_use]
    pub const fn mol_inversion_flag(&self) -> Option<i32> {
        self.mol_inversion_flag
    }

    #[must_use]
    pub const fn implicit_hydrogen(&self) -> bool {
        self.implicit_hydrogen
    }

    #[must_use]
    pub fn tracked_isotopic_hydrogens(&self) -> &[u16] {
        &self.tracked_isotopic_hydrogens
    }

    #[must_use]
    pub const fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }

    #[must_use]
    pub const fn isotope(&self) -> Option<u16> {
        self.isotope
    }

    #[must_use]
    pub const fn atom_map(&self) -> Option<u32> {
        self.atom_map
    }

    #[must_use]
    pub const fn no_implicit(&self) -> bool {
        self.no_implicit
    }

    #[must_use]
    pub const fn radical_electrons(&self) -> u8 {
        self.radical_electrons
    }

    #[must_use]
    pub const fn hybridization(&self) -> Hybridization {
        self.hybridization
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    #[must_use]
    pub fn is_prop_computed(&self, key: &str) -> bool {
        self.computed_props.contains(key)
    }

    #[must_use]
    pub fn computed_prop_names(&self) -> &BTreeSet<String> {
        &self.computed_props
    }

    #[must_use]
    pub const fn pdb_residue_info(&self) -> Option<&AtomPdbResidueInfo> {
        self.pdb_residue_info.as_ref()
    }

    #[must_use]
    pub const fn template_attachment_order(&self) -> Option<&TemplateAttachmentOrder> {
        self.template_attachment_order.as_ref()
    }
}

fn validate_property_key(key: &str) -> Result<(), AtomPropertyError> {
    // BEGIN RDKIT CPP FUNCTION RDProps::setProp empty-key precondition
    // RDKit✔️✔️: if(key.empty()) {
    // RDKit✔️✔️:   throw ValueErrorException("Cannot set property with empty key");
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDProps::setProp empty-key precondition
    if key.is_empty() {
        Err(AtomPropertyError::EmptyKey)
    } else {
        Ok(())
    }
}

/// Immutable atom record owned by `Molecule`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Atom {
    id: AtomId,
    element: Element,
    formal_charge: i8,
    explicit_hydrogens: u8,
    chiral_tag: ChiralTag,
    chiral_permutation: Option<u32>,
    unknown_stereo: bool,
    mol_parity: Option<i32>,
    mol_inversion_flag: Option<i32>,
    implicit_hydrogen: bool,
    tracked_isotopic_hydrogens: Vec<u16>,
    is_aromatic: bool,
    isotope: Option<u16>,
    atom_map: Option<u32>,
    no_implicit: bool,
    radical_electrons: u8,
    hybridization: Hybridization,
    props: BTreeMap<String, String>,
    computed_props: BTreeSet<String>,
    pdb_residue_info: Option<AtomPdbResidueInfo>,
    template_attachment_order: Option<TemplateAttachmentOrder>,
}

impl Atom {
    pub fn from_spec(id: AtomId, spec: AtomSpec) -> Self {
        Self {
            id,
            element: spec.element,
            formal_charge: spec.formal_charge,
            explicit_hydrogens: spec.explicit_hydrogens,
            chiral_tag: spec.chiral_tag,
            chiral_permutation: spec.chiral_permutation,
            unknown_stereo: spec.unknown_stereo,
            mol_parity: spec.mol_parity,
            mol_inversion_flag: spec.mol_inversion_flag,
            implicit_hydrogen: spec.implicit_hydrogen,
            tracked_isotopic_hydrogens: spec.tracked_isotopic_hydrogens,
            is_aromatic: spec.is_aromatic,
            isotope: spec.isotope,
            atom_map: spec.atom_map,
            no_implicit: spec.no_implicit,
            radical_electrons: spec.radical_electrons,
            hybridization: spec.hybridization,
            props: spec.props,
            computed_props: spec.computed_props,
            pdb_residue_info: spec.pdb_residue_info,
            template_attachment_order: spec.template_attachment_order,
        }
    }

    #[doc(hidden)]
    pub fn with_id(mut self, id: AtomId) -> Self {
        self.id = id;
        self
    }

    #[doc(hidden)]
    pub fn set_element(&mut self, element: Element) {
        // RDKit✔️✔️: void setAtomicNum(int newNum) { d_atomicNum = newNum; }
        self.element = element;
    }

    #[must_use]
    pub const fn id(&self) -> AtomId {
        // RDKit✔️✔️: unsigned int getIdx() const { return d_index; }
        self.id
    }

    #[must_use]
    pub const fn element(&self) -> Element {
        self.element
    }

    #[must_use]
    pub const fn atomic_number(&self) -> u8 {
        // RDKit✔️✔️: int getAtomicNum() const { return d_atomicNum; }
        self.element.atomic_number()
    }

    #[must_use]
    pub const fn formal_charge(&self) -> i8 {
        // RDKit✔️✔️: int getFormalCharge() const { return d_formalCharge; }
        self.formal_charge
    }

    #[must_use]
    pub const fn explicit_hydrogens(&self) -> u8 {
        // RDKit✔️✔️: unsigned int getNumExplicitHs() const { return d_numExplicitHs; }
        self.explicit_hydrogens
    }

    #[must_use]
    pub const fn chiral_tag(&self) -> ChiralTag {
        // RDKit✔️✔️: ChiralType getChiralTag() const {
        // RDKit✔️✔️:   return static_cast<ChiralType>(d_chiralTag);
        // RDKit✔️✔️: }
        self.chiral_tag
    }

    #[must_use]
    pub const fn chiral_permutation(&self) -> Option<u32> {
        self.chiral_permutation
    }

    #[must_use]
    pub const fn unknown_stereo(&self) -> bool {
        self.unknown_stereo
    }

    #[must_use]
    pub const fn mol_parity(&self) -> Option<i32> {
        self.mol_parity
    }

    #[must_use]
    pub const fn mol_inversion_flag(&self) -> Option<i32> {
        self.mol_inversion_flag
    }

    #[must_use]
    pub const fn implicit_hydrogen(&self) -> bool {
        self.implicit_hydrogen
    }

    #[must_use]
    pub fn tracked_isotopic_hydrogens(&self) -> &[u16] {
        &self.tracked_isotopic_hydrogens
    }

    #[must_use]
    pub const fn is_aromatic(&self) -> bool {
        // RDKit✔️✔️: bool getIsAromatic() const { return df_isAromatic; }
        self.is_aromatic
    }

    #[must_use]
    pub const fn isotope(&self) -> Option<u16> {
        // RDKit✔️✔️: unsigned int getIsotope() const { return d_isotope; }
        // Source zero is projected as absence in the detached value.
        self.isotope
    }

    #[must_use]
    pub const fn atom_map(&self) -> Option<u32> {
        self.atom_map
    }

    #[must_use]
    pub const fn no_implicit(&self) -> bool {
        // RDKit✔️✔️: bool getNoImplicit() const { return df_noImplicit; }
        self.no_implicit
    }

    #[must_use]
    pub const fn radical_electrons(&self) -> u8 {
        // RDKit✔️✔️: unsigned int getNumRadicalElectrons() const { return d_numRadicalElectrons; }
        self.radical_electrons
    }

    #[must_use]
    pub const fn hybridization(&self) -> Hybridization {
        // RDKit✔️✔️: HybridizationType getHybridization() const {
        // RDKit✔️✔️:   return static_cast<HybridizationType>(d_hybrid);
        // RDKit✔️✔️: }
        self.hybridization
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    /// Returns whether a property is registered as computed state.
    #[must_use]
    pub fn is_prop_computed(&self, key: &str) -> bool {
        self.computed_props.contains(key)
    }

    #[must_use]
    pub fn computed_prop_names(&self) -> &BTreeSet<String> {
        &self.computed_props
    }

    /// Returns the modern CIP descriptor persisted on this atom, if present.
    pub fn cip_descriptor(
        &self,
    ) -> Result<Option<crate::CipDescriptor>, crate::CipDescriptorError> {
        crate::cip::descriptor_from_property(self.prop("_CIPCode"))
    }

    #[must_use]
    pub const fn pdb_residue_info(&self) -> Option<&AtomPdbResidueInfo> {
        self.pdb_residue_info.as_ref()
    }

    #[must_use]
    pub const fn template_attachment_order(&self) -> Option<&TemplateAttachmentOrder> {
        self.template_attachment_order.as_ref()
    }

    /// Apply the shared detached atom-index remap to typed attachment state.
    #[doc(hidden)]
    pub fn remap_template_attachment_order(
        &mut self,
        old_to_new: &[Option<AtomId>],
    ) -> Result<(), TemplateAttachmentOrderError> {
        let Some(order) = &self.template_attachment_order else {
            return Ok(());
        };
        let remapped = order.remapped(old_to_new)?;
        self.template_attachment_order = Some(remapped);
        Ok(())
    }

    #[doc(hidden)]
    pub fn set_chiral_tag(&mut self, chiral_tag: ChiralTag) {
        // RDKit✔️✔️: void setChiralTag(ChiralType what) { d_chiralTag = what; }
        self.chiral_tag = chiral_tag;
    }

    #[doc(hidden)]
    pub fn set_chiral_permutation(&mut self, chiral_permutation: Option<u32>) {
        self.chiral_permutation = chiral_permutation;
    }

    #[doc(hidden)]
    pub fn set_unknown_stereo(&mut self, unknown_stereo: bool) {
        self.unknown_stereo = unknown_stereo;
    }

    #[doc(hidden)]
    pub fn set_mol_parity(&mut self, mol_parity: Option<i32>) {
        self.mol_parity = mol_parity;
    }

    #[doc(hidden)]
    pub fn set_mol_inversion_flag(&mut self, mol_inversion_flag: Option<i32>) {
        self.mol_inversion_flag = mol_inversion_flag;
    }

    #[doc(hidden)]
    pub fn set_implicit_hydrogen(&mut self, implicit_hydrogen: bool) {
        self.implicit_hydrogen = implicit_hydrogen;
    }

    #[doc(hidden)]
    pub fn set_tracked_isotopic_hydrogens(&mut self, isotopes: Vec<u16>) {
        self.tracked_isotopic_hydrogens = isotopes;
    }

    #[doc(hidden)]
    pub fn set_aromatic(&mut self, is_aromatic: bool) {
        // RDKit✔️✔️: void setIsAromatic(bool what) { df_isAromatic = what; }
        self.is_aromatic = is_aromatic;
    }

    #[doc(hidden)]
    pub fn set_formal_charge(&mut self, formal_charge: i8) {
        // RDKit✔️✔️: void setFormalCharge(int what) { d_formalCharge = what; }
        self.formal_charge = formal_charge;
    }

    #[doc(hidden)]
    pub fn set_explicit_hydrogens(&mut self, explicit_hydrogens: u8) {
        // RDKit✔️✔️: void setNumExplicitHs(unsigned int what) { d_numExplicitHs = what; }
        self.explicit_hydrogens = explicit_hydrogens;
    }

    #[doc(hidden)]
    pub fn set_isotope(&mut self, isotope: Option<u16>) {
        // RDKit✔️✔️: void Atom::setIsotope(unsigned int what) { d_isotope = what; }
        // Source zero is projected as `None` in the detached value.
        self.isotope = isotope;
    }

    #[doc(hidden)]
    pub fn set_atom_map(&mut self, atom_map: Option<u32>) {
        self.atom_map = atom_map;
    }

    #[doc(hidden)]
    pub fn set_no_implicit(&mut self, no_implicit: bool) {
        // RDKit✔️✔️: void setNoImplicit(bool what) { df_noImplicit = what; }
        self.no_implicit = no_implicit;
    }

    #[doc(hidden)]
    pub fn set_radical_electrons(&mut self, radical_electrons: u8) {
        // RDKit✔️✔️: void setNumRadicalElectrons(unsigned int num) { d_numRadicalElectrons = num; }
        self.radical_electrons = radical_electrons;
    }

    #[doc(hidden)]
    pub fn set_prop(
        &mut self,
        key: impl Into<String>,
        value: impl Into<String>,
    ) -> Result<(), AtomPropertyError> {
        let key = key.into();
        validate_property_key(&key)?;
        // RDKit✔️✔️: d_props.setVal(key, val);
        // A non-computed write does not remove an existing computed marker.
        self.props.insert(key, value.into());
        Ok(())
    }

    #[doc(hidden)]
    pub fn set_computed_prop(
        &mut self,
        key: impl Into<String>,
        value: impl Into<String>,
    ) -> Result<(), AtomPropertyError> {
        // RDKit✔️🔝: if (computed) {
        // RDKit✔️🔝:   STR_VECT compLst;
        // RDKit✔️🔝:   getPropIfPresent(RDKit::detail::computedPropName, compLst);
        // RDKit✔️🔝:   if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
        // RDKit✔️🔝:     compLst.emplace_back(key);
        // RDKit✔️🔝:     d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // RDKit✔️🔝: d_props.setVal(key, val);
        // The ordered set preserves membership semantics while replacing the
        // source vector's linear duplicate scan with logarithmic insertion.
        let key = key.into();
        validate_property_key(&key)?;
        self.props.insert(key.clone(), value.into());
        self.computed_props.insert(key);
        Ok(())
    }

    #[doc(hidden)]
    pub fn clear_prop(&mut self, key: &str) {
        // RDKit✔️🔝: auto svi = std::find(compLst.begin(), compLst.end(), key);
        // RDKit✔️🔝: if (svi != compLst.end()) {
        // RDKit✔️🔝:   compLst.erase(svi);
        // RDKit✔️🔝:   d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️🔝: }
        // RDKit✔️🔝: d_props.clearVal(key);
        // BTreeSet removal preserves the source transition with logarithmic
        // lookup instead of the source vector's linear search and erase.
        self.props.remove(key);
        self.computed_props.remove(key);
    }

    #[doc(hidden)]
    pub fn clear_computed_props(&mut self) {
        // RDKit✔️🔝: for (const auto &key : compLst) {
        // RDKit✔️🔝:   d_props.clearVal(key);
        // RDKit✔️🔝: }
        // Moving the set avoids the source vector copy while preserving exact
        // membership-based clearing, including non-computed properties that
        // happen to use a conventional computed-property name.
        for key in std::mem::take(&mut self.computed_props) {
            self.props.remove(&key);
        }
    }

    #[doc(hidden)]
    pub fn set_hybridization(&mut self, hybridization: Hybridization) {
        // RDKit✔️✔️: void setHybridization(HybridizationType what) { d_hybrid = what; }
        self.hybridization = hybridization;
    }

    #[doc(hidden)]
    pub fn set_pdb_residue_info(&mut self, info: Option<AtomPdbResidueInfo>) {
        self.pdb_residue_info = info;
    }
}
