//! Source-backed detached atom-code implementation.
//! Pinned RDKit351f8f378f8ad6bbd517980c38896e66bf907af8.
//! Reuse foundational numPi, property string and modern CIP owners.
use cosmolkit_core::{
    PropertyStringError, ValenceError, num_pi_electrons_for_topology, property_value_to_string,
};
use cosmolkit_model::{
    AtomId, ChiralTag, MoleculeProperties, TopologyBlock, TopologyValidationError,
};
use cosmolkit_stereo::{CipLabelOptions, CipLabelerError, assign_cip_labels};
use std::{borrow::Cow, fmt};

/// Exact source table, including its implicit final zero.
pub(crate) const ATOM_NUMBER_TYPES: [u32; 16] =
    [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53, 0];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomCodeOptions {
    pub branch_subtract: u32,
    pub include_chirality: bool,
    /// Already resolved source stereo configuration; no global environment IO.
    pub use_legacy_stereo_perception: bool,
}
impl Default for AtomCodeOptions {
    fn default() -> Self {
        Self {
            branch_subtract: 0,
            include_chirality: false,
            use_legacy_stereo_perception: true,
        }
    }
}

/// Validated canonical detached blocks, without runtime/cache ownership.
/// Validate once during construction; returned input can be reused for another
/// selected atom without whole-state cloning or per-atom global validation.
#[derive(Debug)]
pub struct AtomCodeInput<'a> {
    topology: Cow<'a, TopologyBlock>,
    properties: Cow<'a, MoleculeProperties>,
}
impl AtomCodeInput<'static> {
    pub fn new(
        topology: TopologyBlock,
        properties: MoleculeProperties,
    ) -> Result<Self, TopologyValidationError> {
        topology.validate()?;
        Ok(Self {
            topology: Cow::Owned(topology),
            properties: Cow::Owned(properties),
        })
    }
}
impl<'a> AtomCodeInput<'a> {
    /// Scoped COW input from a declared runtime capability, with local validation.
    pub fn from_cow(
        topology: Cow<'a, TopologyBlock>,
        properties: Cow<'a, MoleculeProperties>,
    ) -> Result<Self, TopologyValidationError> {
        topology.validate()?;
        Ok(Self {
            topology,
            properties,
        })
    }
    /// Borrow already detached state; modern CIP creates a private copy only on its source guard.
    pub fn from_ref(
        topology: &'a TopologyBlock,
        properties: &'a MoleculeProperties,
    ) -> Result<Self, TopologyValidationError> {
        topology.validate()?;
        Ok(Self {
            topology: Cow::Borrowed(topology),
            properties: Cow::Borrowed(properties),
        })
    }
    pub fn topology(&self) -> &TopologyBlock {
        &self.topology
    }
    pub fn properties(&self) -> &MoleculeProperties {
        &self.properties
    }
    pub fn into_parts(self) -> (TopologyBlock, MoleculeProperties) {
        (self.topology.into_owned(), self.properties.into_owned())
    }
}
#[derive(Debug)]
pub struct AtomCodeAssignment<'a> {
    code: u32,
    input: AtomCodeInput<'a>,
}
impl<'a> AtomCodeAssignment<'a> {
    /// Return the changed owner blocks only when source assignment created them.
    /// A borrowed result explicitly records no state effect, preserving COW.
    pub fn into_optional_owned_parts(self) -> (u32, Option<(TopologyBlock, MoleculeProperties)>) {
        let pair = match (self.input.topology, self.input.properties) {
            (Cow::Borrowed(_), Cow::Borrowed(_)) => None,
            (topology, properties) => Some((topology.into_owned(), properties.into_owned())),
        };
        (self.code, pair)
    }
    pub fn code(&self) -> u32 {
        self.code
    }
    pub fn topology(&self) -> &TopologyBlock {
        self.input.topology()
    }
    pub fn properties(&self) -> &MoleculeProperties {
        self.input.properties()
    }
    pub fn into_input(self) -> AtomCodeInput<'a> {
        self.input
    }
    pub fn into_parts(self) -> (u32, TopologyBlock, MoleculeProperties) {
        let (topology, properties) = self.input.into_parts();
        (self.code, topology, properties)
    }
}
#[derive(Debug, Clone, PartialEq)]
pub enum AtomCodeError {
    Precondition { what: &'static str },
    Postcondition { what: &'static str, code: u32 },
    Valence(ValenceError),
    Cip(CipLabelerError),
    PropertyString(PropertyStringError),
}
impl fmt::Display for AtomCodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Precondition { what } => f.write_str(what),
            Self::Postcondition { what, .. } => f.write_str(what),
            Self::Valence(source) => write!(f, "{source}"),
            Self::Cip(source) => write!(f, "{source}"),
            Self::PropertyString(source) => write!(f, "{source}"),
        }
    }
}
impl std::error::Error for AtomCodeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Valence(source) => Some(source),
            Self::Cip(source) => Some(source),
            Self::PropertyString(source) => Some(source),
            Self::Precondition { .. } | Self::Postcondition { .. } => None,
        }
    }
}
/// Source getAtomCode, with all owning-molecule CIP effects returned explicitly.
/// `explicit_valence` is the selected stored i8 cache, not a recomputation.
/// Null atom representation fails before chemistry; unsupported independent
/// source value kinds and modern CIP failures propagate from their sole owners.
pub fn atom_code<'a>(
    mut input: AtomCodeInput<'a>,
    atom_id: Option<AtomId>,
    explicit_valence: Option<i8>,
    options: &AtomCodeOptions,
) -> Result<AtomCodeAssignment<'a>, AtomCodeError> {
    // RDKit❌❌: const unsigned int numTypeBits = 4;
    // RDKit❌❌: const unsigned int atomNumberTypes[1 << numTypeBits] = {
    // RDKit❌❌:     5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53};
    // RDKit❌❌: const unsigned int numPiBits = 2;
    // RDKit❌❌: const unsigned int maxNumPi = (1 << numPiBits) - 1;
    // RDKit❌❌: const unsigned int numBranchBits = 3;
    // RDKit❌❌: const unsigned int maxNumBranches = (1 << numBranchBits) - 1;
    // RDKit❌❌: const unsigned int numChiralBits = 2;
    // RDKit❌❌: const unsigned int codeSize = numTypeBits + numPiBits + numBranchBits;
    // RDKit❌❌: const unsigned int numPathBits = 5;
    // RDKit❌❌: const unsigned int maxPathLen = (1 << numPathBits) - 1;
    // RDKit❌❌: const unsigned int numAtomPairFingerprintBits =
    // RDKit❌❌:     numPathBits + 2 * codeSize;  // note that this is only accurate if chirality
    // RDKit❌❌: std::uint32_t getAtomCode(const Atom *atom, unsigned int branchSubtract,
    // RDKit❌❌:                           bool includeChirality) {
    // RDKit❌❌:   PRECONDITION(atom, "no atom");
    // RDKit❌❌:   std::uint32_t code;
    // RDKit❌❌:
    // RDKit❌❌:   unsigned int numBranches = 0;
    // RDKit❌❌:   if (atom->getDegree() > branchSubtract) {
    // RDKit❌❌:     numBranches = atom->getDegree() - branchSubtract;
    // RDKit❌❌:   }
    // RDKit❌❌:
    // RDKit❌❌:   code = numBranches % maxNumBranches;
    // RDKit❌❌:   unsigned int nPi = numPiElectrons(*atom) % maxNumPi;
    // RDKit❌❌:   code |= nPi << numBranchBits;
    // RDKit❌❌:
    // RDKit❌❌:   unsigned int typeIdx = 0;
    // RDKit❌❌:   unsigned int nTypes = 1 << numTypeBits;
    // RDKit❌❌:   while (typeIdx < nTypes) {
    // RDKit❌❌:     if (atomNumberTypes[typeIdx] ==
    // RDKit❌❌:         static_cast<unsigned int>(atom->getAtomicNum())) {
    // RDKit❌❌:       break;
    // RDKit❌❌:     } else if (atomNumberTypes[typeIdx] >
    // RDKit❌❌:                static_cast<unsigned int>(atom->getAtomicNum())) {
    // RDKit❌❌:       typeIdx = nTypes;
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit❌❌:     ++typeIdx;
    // RDKit❌❌:   }
    // RDKit❌❌:   if (typeIdx == nTypes) {
    // RDKit❌❌:     --typeIdx;
    // RDKit❌❌:   }
    // RDKit❌❌:   code |= typeIdx << (numBranchBits + numPiBits);
    // RDKit❌❌:   if (includeChirality) {
    // RDKit❌❌:     // if we aren't using legacy stereo, we need to compute the CIP codes
    // RDKit❌❌:     if (!Chirality::getUseLegacyStereoPerception() &&
    // RDKit❌❌:         atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❌❌:         !atom->getOwningMol().hasProp(common_properties::_CIPComputed)) {
    // RDKit❌❌:       CIPLabeler::assignCIPLabels(atom->getOwningMol());
    // RDKit❌❌:     }
    // RDKit❌❌:     std::string cipCode;
    // RDKit❌❌:     if (atom->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
    // RDKit❌❌:       std::uint32_t offset = numBranchBits + numPiBits + numTypeBits;
    // RDKit❌❌:       if (cipCode == "R") {
    // RDKit❌❌:         code |= 1 << offset;
    // RDKit❌❌:       } else if (cipCode == "S") {
    // RDKit❌❌:         code |= 2 << offset;
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit❌❌:   POSTCONDITION(code < static_cast<std::uint32_t>(
    // RDKit❌❌:                            1 << (codeSize + (includeChirality ? 2 : 0))),
    // RDKit❌❌:                 "code exceeds number of bits");
    // RDKit❌❌:   return code;
    // RDKit❌❌: };
    // RDKit❌❌: unsigned int Atom::getDegree() const {
    // RDKit❌❌:   return dp_mol ? getOwningMol().getAtomDegree(this) : 0;
    // RDKit❌❌: }
    // RDKit❌❌:   template <typename T>
    // RDKit❌❌:   bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit❌❌:     return d_props.getValIfPresent(key, res);
    // RDKit❌❌:   }
    // RDKit❌❌:
    // RDKit❌❌:   //! \overload
    // RDKit❌❌:   bool hasProp(const std::string_view key) const { return d_props.hasVal(key); }
    // RDKit❌❌:
    // RDKit❌❌:   bool getValIfPresent(const std::string_view what, std::string &res) const {
    // RDKit❌❌:     for (const auto &i : _data) {
    // RDKit❌❌:       if (i.key == what) {
    // RDKit❌❌:         rdvalue_tostring(i.val, res);
    // RDKit❌❌:         return true;
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:     return false;
    // RDKit❌❌:   }
    // RDKit❌❌:
    // RDKit❌❌: inline bool rdvalue_tostring(RDValue_cast_t val, std::string &res) {
    // RDKit❌❌:   switch (val.getTag()) {
    // RDKit❌❌:     case RDTypeTag::StringTag:
    // RDKit❌❌:       res = rdvalue_cast<std::string>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     case RDTypeTag::IntTag:
    // RDKit❌❌:       res = boost::lexical_cast<std::string>(rdvalue_cast<int>(val));
    // RDKit❌❌:       break;
    // RDKit❌❌:     case RDTypeTag::DoubleTag: {
    // RDKit❌❌:       Utils::LocaleSwitcher ls;  // for lexical cast...
    // RDKit❌❌:       res = boost::lexical_cast<std::string>(rdvalue_cast<double>(val));
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit❌❌:     case RDTypeTag::UnsignedIntTag:
    // RDKit❌❌:       res = boost::lexical_cast<std::string>(rdvalue_cast<unsigned int>(val));
    // RDKit❌❌:       break;
    // RDKit❌❌: #ifdef RDVALUE_HASBOOL
    // RDKit❌❌:     case RDTypeTag::BoolTag:
    // RDKit❌❌:       res = boost::lexical_cast<std::string>(rdvalue_cast<bool>(val));
    // RDKit❌❌:       break;
    // RDKit❌❌: #endif
    // RDKit❌❌:     case RDTypeTag::FloatTag: {
    // RDKit❌❌:       Utils::LocaleSwitcher ls;  // for lexical cast...
    // RDKit❌❌:       res = boost::lexical_cast<std::string>(rdvalue_cast<float>(val));
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit❌❌:     case RDTypeTag::VecDoubleTag: {
    // RDKit❌❌:       // vectToString uses std::imbue for locale
    // RDKit❌❌:       res = vectToString<double>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit❌❌:     case RDTypeTag::VecFloatTag: {
    // RDKit❌❌:       // vectToString uses std::imbue for locale
    // RDKit❌❌:       res = vectToString<float>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit❌❌:     case RDTypeTag::VecIntTag:
    // RDKit❌❌:       res = vectToString<int>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     case RDTypeTag::VecUnsignedIntTag:
    // RDKit❌❌:       res = vectToString<unsigned int>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     case RDTypeTag::VecStringTag:
    // RDKit❌❌:       res = vectToString<std::string>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     case RDTypeTag::AnyTag: {
    // RDKit❌❌:       Utils::LocaleSwitcher ls;  // for lexical cast...
    // RDKit❌❌:       try {
    // RDKit❌❌:         res = std::any_cast<std::string>(rdvalue_cast<std::any &>(val));
    // RDKit❌❌:       } catch (const std::bad_any_cast &) {
    // RDKit❌❌:         auto &rdtype = rdvalue_cast<std::any &>(val).type();
    // RDKit❌❌:         if (rdtype == typeid(long)) {
    // RDKit❌❌:           res = boost::lexical_cast<std::string>(
    // RDKit❌❌:               std::any_cast<long>(rdvalue_cast<std::any &>(val)));
    // RDKit❌❌:         } else if (rdtype == typeid(int64_t)) {
    // RDKit❌❌:           res = boost::lexical_cast<std::string>(
    // RDKit❌❌:               std::any_cast<int64_t>(rdvalue_cast<std::any &>(val)));
    // RDKit❌❌:         } else if (rdtype == typeid(uint64_t)) {
    // RDKit❌❌:           res = boost::lexical_cast<std::string>(
    // RDKit❌❌:               std::any_cast<uint64_t>(rdvalue_cast<std::any &>(val)));
    // RDKit❌❌:         } else if (rdtype == typeid(unsigned long)) {
    // RDKit❌❌:           res = boost::lexical_cast<std::string>(
    // RDKit❌❌:               std::any_cast<unsigned long>(rdvalue_cast<std::any &>(val)));
    // RDKit❌❌:         } else {
    // RDKit❌❌:           throw;
    // RDKit❌❌:           return false;
    // RDKit❌❌:         }
    // RDKit❌❌:       }
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit❌❌:     default:
    // RDKit❌❌:       res = "";
    // RDKit❌❌:   }
    // RDKit❌❌:   return true;
    // RDKit❌❌: }
    let atom_id = atom_id.ok_or(AtomCodeError::Precondition { what: "no atom" })?;
    let atom = input.topology.atoms.get(atom_id.index()).ok_or_else(|| {
        AtomCodeError::Valence(ValenceError::AtomOutOfRange {
            atom: atom_id,
            atom_count: input.topology.atoms.len(),
        })
    })?;
    let degree = input.topology.adjacency.neighbors_of(atom_id.index()).len() as u32;
    let branches = if degree > options.branch_subtract {
        degree - options.branch_subtract
    } else {
        0
    };
    let mut code = branches % 7;
    let n_pi = num_pi_electrons_for_topology(&input.topology, atom_id, explicit_valence)
        .map_err(AtomCodeError::Valence)?
        % 3;
    code |= n_pi << 3;
    // The C++ array has15 initializers and an implicit final zero. Preserve
    // the exact unsigned comparison order; no periodictable heuristic.

    let atomic_number = u32::from(atom.atomic_number());
    let mut type_index = 0usize;
    while type_index < 16 {
        if ATOM_NUMBER_TYPES[type_index] == atomic_number {
            break;
        }
        if ATOM_NUMBER_TYPES[type_index] > atomic_number {
            type_index = 16;
            break;
        }
        type_index += 1;
    }
    if type_index == 16 {
        type_index -= 1;
    }
    code |= (type_index as u32) << 5;
    let needs_cip = options.include_chirality
        && !options.use_legacy_stereo_perception
        && atom.chiral_tag() != ChiralTag::Unspecified
        && input.properties.prop("_CIPComputed").is_none();
    if needs_cip {
        let assignment = assign_cip_labels(
            input.topology.into_owned(),
            input.properties.into_owned(),
            &CipLabelOptions::default(),
        )
        .map_err(AtomCodeError::Cip)?;
        let (topology, properties) = assignment.into_parts();
        input = AtomCodeInput {
            topology: Cow::Owned(topology),
            properties: Cow::Owned(properties),
        };
    }
    if options.include_chirality {
        if let Some(value) = input.topology.atoms[atom_id.index()].prop("_CIPCode") {
            let cip_code =
                property_value_to_string(value).map_err(AtomCodeError::PropertyString)?;
            if cip_code.as_bytes() == b"R" {
                code |= 1 << 9;
            } else if cip_code.as_bytes() == b"S" {
                code |= 2 << 9;
            }
        }
    }
    let bits = 9 + if options.include_chirality { 2 } else { 0 };
    if code >= (1u32 << bits) {
        return Err(AtomCodeError::Postcondition {
            what: "code exceeds number of bits",
            code,
        });
    }
    Ok(AtomCodeAssignment { code, input })
}
