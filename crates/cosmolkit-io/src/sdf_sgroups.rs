//! Typed V3000 SGroup and enhanced-stereo collection lowering.

use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_model::{
    AtomId, BondId, PropertyText, PropertyValue, SGroupAttachPoint, SGroupBondRole, SGroupBracket,
    SGroupBracketStyle, SGroupCState, SGroupConnection, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock,
};

use crate::sdf::{
    SdfReadError, SdfWriteError, get_v3000_line, parse_rdkit_atof, parse_rdkit_double,
    parse_rdkit_int, parse_rdkit_unsigned, rdkit_substr,
};

fn sgroup_kind_from_rdkit_type(value: &str) -> SubstanceGroupKind {
    match value {
        "DAT" => SubstanceGroupKind::Data,
        "SUP" => SubstanceGroupKind::Superatom,
        "MUL" => SubstanceGroupKind::MultipleGroup,
        "SRU" => SubstanceGroupKind::StructuralRepeatUnit,
        "MON" => SubstanceGroupKind::Monomer,
        "COP" => SubstanceGroupKind::Copolymer,
        "CRO" => SubstanceGroupKind::Crosslink,
        "GRA" => SubstanceGroupKind::Graft,
        "MOD" => SubstanceGroupKind::Modification,
        "MER" => SubstanceGroupKind::Mer,
        "ANY" => SubstanceGroupKind::AnyPolymer,
        "COM" => SubstanceGroupKind::MixtureComponent,
        "MIX" => SubstanceGroupKind::Mixture,
        "FOR" => SubstanceGroupKind::Formulation,
        other => SubstanceGroupKind::Generic(other.into()),
    }
}

fn is_valid_rdkit_sgroup_type(value: &str) -> bool {
    matches!(
        value,
        "SRU"
            | "MON"
            | "COP"
            | "CRO"
            | "GRA"
            | "MOD"
            | "MER"
            | "ANY"
            | "COM"
            | "MIX"
            | "FOR"
            | "SUP"
            | "MUL"
            | "DAT"
            | "GEN"
    )
}

fn is_valid_rdkit_sgroup_subtype(value: &str) -> bool {
    // BEGIN RDKIT CPP FUNCTION SubstanceGroupChecks::isValidSubType
    // RDKit✔️✔️: const std::vector<std::string> sGroupSubtypes = {"ALT", "RAN", "BLO"};
    // RDKit✔️✔️: bool SubstanceGroupChecks::isValidSubType(const std::string &type) {
    // RDKit✔️✔️:   return std::find(SubstanceGroupChecks::sGroupSubtypes.begin(),
    // RDKit✔️✔️:                    SubstanceGroupChecks::sGroupSubtypes.end(),
    // RDKit✔️✔️:                    type) != SubstanceGroupChecks::sGroupSubtypes.end();
    // RDKit✔️✔️: }
    matches!(value, "ALT" | "RAN" | "BLO")
    // The fixed three-value match is behaviorally identical to the source
    // vector search and has constant time with no allocation.
    // END RDKIT CPP FUNCTION
}

fn sgroup_connection_from_rdkit(value: &str) -> Option<SGroupConnection> {
    // BEGIN RDKIT CPP FUNCTION SubstanceGroupChecks::isValidConnectType
    // RDKit✔️✔️: const std::vector<std::string> sGroupConnectTypes = {"HH", "HT", "EU"};
    // RDKit✔️✔️: bool SubstanceGroupChecks::isValidConnectType(const std::string &type) {
    // RDKit✔️✔️:   return std::find(SubstanceGroupChecks::sGroupConnectTypes.begin(),
    // RDKit✔️✔️:                    SubstanceGroupChecks::sGroupConnectTypes.end(),
    // RDKit✔️✔️:                    type) != SubstanceGroupChecks::sGroupConnectTypes.end();
    // RDKit✔️✔️: }
    match value {
        "HH" => Some(SGroupConnection::HeadToHead),
        "HT" => Some(SGroupConnection::HeadToTail),
        "EU" => Some(SGroupConnection::Either),
        _ => None,
    }
    // The source validates one of three exact strings and stores that string.
    // The canonical enum stores the same three states without a second raw
    // property. Matching is constant time and allocation-free.
    // END RDKIT CPP FUNCTION
}

fn is_valid_rdkit_sgroup_class(value: &str) -> bool {
    // BEGIN RDKIT CPP FUNCTION SubstanceGroupChecks::isValidClass
    // RDKit✔️✔️: const std::vector<std::string> sGroupClasses = {
    // RDKit✔️✔️:     "AA",        "dAA",    "DNA",     "RNA",      "SUGAR",    "BASE",
    // RDKit✔️✔️:     "PHOSPHATE", "LINKER", "CHEM",    "LGRP",     "MODAA",    "MODdAA",
    // RDKit✔️✔️:     "MODDNA",    "MODRNA", "XLINKAA", "XLINKdAA", "XLINKDNA", "XLINKRNA",
    // RDKit✔️✔️: };
    // RDKit✔️✔️: bool SubstanceGroupChecks::isValidClass(const std::string &sgroupClass) {
    // RDKit✔️✔️:   return std::find(SubstanceGroupChecks::sGroupClasses.begin(),
    // RDKit✔️✔️:                    SubstanceGroupChecks::sGroupClasses.end(),
    // RDKit✔️✔️:                    sgroupClass) != SubstanceGroupChecks::sGroupClasses.end();
    // RDKit✔️✔️: }
    matches!(
        value,
        "AA" | "dAA"
            | "DNA"
            | "RNA"
            | "SUGAR"
            | "BASE"
            | "PHOSPHATE"
            | "LINKER"
            | "CHEM"
            | "LGRP"
            | "MODAA"
            | "MODdAA"
            | "MODDNA"
            | "MODRNA"
            | "XLINKAA"
            | "XLINKdAA"
            | "XLINKDNA"
            | "XLINKRNA"
    )
    // The fixed match has the exact case-sensitive source vocabulary, does
    // not allocate, and is no worse than the source's bounded linear search.
    // END RDKIT CPP FUNCTION
}

fn split_sgroup_line(
    line: &str,
    line_number: usize,
) -> Result<(&str, &str, &str, &str), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock
    // RDKit❗✔️: std::string ParseV3000SGroupsBlock(std::istream *inStream, unsigned int &line,
    // RDKit❗✔️:                                    unsigned int nSgroups, RWMol *mol,
    // RDKit❗✔️:                                    bool strictParsing) {
    // RDKit❗✔️:   PRECONDITION(inStream, "no stream");
    // RDKit❗✔️:   PRECONDITION(mol, "no molecule");
    // RDKit❗✔️:   unsigned int defaultLineNum = 0;
    // RDKit❗✔️:   std::string defaultString;
    // RDKit❗✔️:
    // RDKit❗✔️:   // SGroups may be written in unsorted ID order, according to spec, so we will
    // RDKit❗✔️:   // temporarily store them in a map before adding them to the mol
    // RDKit❗✔️:   IDX_TO_SGROUP_MAP sGroupMap;
    // RDKit❗✔️:
    // RDKit❗✔️:   std::unordered_map<std::string, std::stringstream> defaultLabels;
    // RDKit❗✔️:
    // RDKit❗✔️:   auto tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit❗✔️:
    // RDKit❗✔️:   // Store defaults
    // RDKit❗✔️:   if (tempStr.substr(0, 7) == "DEFAULT" && tempStr.length() > 8) {
    // RDKit❗✔️:     defaultString = tempStr.substr(7);
    // RDKit❗✔️:     defaultLineNum = line;
    // RDKit❗✔️:     boost::trim_right(defaultString);
    // RDKit❗✔️:     tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit❗✔️:     boost::trim_right(tempStr);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   for (unsigned int si = 0; si < nSgroups; ++si) {
    // RDKit❗✔️:     unsigned int sequenceId;
    // RDKit❗✔️:     unsigned int externalId;
    // RDKit❗✔️:     std::string type;
    // RDKit❗✔️:
    // RDKit❗✔️:     std::stringstream lineStream(tempStr);
    // RDKit❗✔️:     lineStream >> sequenceId;
    // RDKit❗✔️:     lineStream >> type;
    // RDKit❗✔️:     lineStream >> externalId;
    // RDKit❗✔️:
    // RDKit❗✔️:     std::set<std::string> parsedLabels;
    // RDKit❗✔️:     if (strictParsing && !SubstanceGroupChecks::isValidType(type)) {
    // RDKit❗✔️:       std::ostringstream errout;
    // RDKit❗✔️:       errout << "Unsupported SGroup type '" << type << "' on line " << line;
    // RDKit❗✔️:       throw MolFileUnhandledFeatureException(errout.str());
    // RDKit❗✔️:     } else if (!strictParsing &&
    // RDKit❗✔️:                nSgroups == std::numeric_limits<unsigned int>::max() &&
    // RDKit❗✔️:                lineStream.fail()) {
    // RDKit❗✔️:       // something went wrong and we didn't know how many SGroups to expect, and
    // RDKit❗✔️:       // now we have seen something that doesn't look like an SGroup start.
    // RDKit❗✔️:       // So we assume we're done.
    // RDKit❗✔️:       nSgroups = 0;
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     SubstanceGroup sgroup(mol, type);
    // RDKit❗✔️:     STR_VECT dataFields;
    // RDKit❗✔️:
    // RDKit❗✔️:     sgroup.setProp<unsigned int>("index", sequenceId);
    // RDKit❗✔️:     if (externalId > 0) {
    // RDKit❗✔️:       if (!SubstanceGroupChecks::isSubstanceGroupIdFree(*mol, externalId)) {
    // RDKit❗✔️:         std::ostringstream errout;
    // RDKit❗✔️:         errout << "Existing SGroup ID '" << externalId
    // RDKit❗✔️:                << "' assigned to a second SGroup on line " << line;
    // RDKit❗✔️:         if (strictParsing) {
    // RDKit❗✔️:           throw FileParseException(errout.str());
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:           sgroup.setIsValid(false);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       sgroup.setProp<unsigned int>("ID", externalId);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     while (sgroup.getIsValid() && !lineStream.eof() && !lineStream.fail()) {
    // RDKit❗✔️:       char spacer;
    // RDKit❗✔️:       std::string label;
    // RDKit❗✔️:
    // RDKit❗✔️:       lineStream.get(spacer);
    // RDKit❗✔️:       if (lineStream.gcount() == 0) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       } else if (spacer != ' ') {
    // RDKit❗✔️:         std::ostringstream errout;
    // RDKit❗✔️:         errout << "Found character '" << spacer
    // RDKit❗✔️:                << "' when expecting a separator (space) on line " << line;
    // RDKit❗✔️:         if (strictParsing) {
    // RDKit❗✔️:           throw FileParseException(errout.str());
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:           sgroup.setIsValid(false);
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       std::getline(lineStream, label, '=');
    // RDKit❗✔️:       if (label.empty()) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       ParseV3000ParseLabel(label, lineStream, dataFields, line, sgroup,
    // RDKit❗✔️:                            nSgroups, mol, strictParsing);
    // RDKit❗✔️:       parsedLabels.insert(label);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // Process defaults
    // RDKit❗✔️:     lineStream.clear();
    // RDKit❗✔️:     lineStream.str(defaultString);
    // RDKit❗✔️:     while (sgroup.getIsValid() && !lineStream.eof() && !lineStream.fail()) {
    // RDKit❗✔️:       char spacer;
    // RDKit❗✔️:       std::string label;
    // RDKit❗✔️:
    // RDKit❗✔️:       lineStream.get(spacer);
    // RDKit❗✔️:       if (lineStream.gcount() == 0) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       } else if (spacer != ' ') {
    // RDKit❗✔️:         std::ostringstream errout;
    // RDKit❗✔️:         errout << "Found character '" << spacer
    // RDKit❗✔️:                << "' when expecting a separator (space) in DEFAULTS on line "
    // RDKit❗✔️:                << defaultLineNum;
    // RDKit❗✔️:         if (strictParsing) {
    // RDKit❗✔️:           throw FileParseException(errout.str());
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:           sgroup.setIsValid(false);
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       std::getline(lineStream, label, '=');
    // RDKit❗✔️:       if (label.empty()) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (std::find(parsedLabels.begin(), parsedLabels.end(), label) ==
    // RDKit❗✔️:           parsedLabels.end()) {
    // RDKit❗✔️:         ParseV3000ParseLabel(label, lineStream, dataFields, defaultLineNum,
    // RDKit❗✔️:                              sgroup, nSgroups, mol, strictParsing);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         spacer = lineStream.peek();
    // RDKit❗✔️:         if (spacer == ' ') {
    // RDKit❗✔️:           std::ostringstream errout;
    // RDKit❗✔️:           errout << "Found unexpected whitespace at DEFAULT label " << label;
    // RDKit❗✔️:           if (strictParsing) {
    // RDKit❗✔️:             throw FileParseException(errout.str());
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:             sgroup.setIsValid(false);
    // RDKit❗✔️:             continue;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         } else if (spacer == '(') {
    // RDKit❗✔️:           std::getline(lineStream, label, ')');
    // RDKit❗✔️:           lineStream.get(spacer);
    // RDKit❗✔️:         } else if (spacer == '"') {
    // RDKit❗✔️:           lineStream.get(spacer);
    // RDKit❗✔️:           std::getline(lineStream, label, '"');
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           std::getline(lineStream, label, ' ');
    // RDKit❗✔️:           lineStream.putback(' ');
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     sgroup.setProp("DATAFIELDS", dataFields);
    // RDKit❗✔️:     sGroupMap.emplace(sequenceId, sgroup);
    // RDKit❗✔️:
    // RDKit❗✔️:     tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit❗✔️:     boost::trim_right(tempStr);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (sGroupMap.size() != nSgroups) {
    // RDKit❗✔️:     std::ostringstream errout;
    // RDKit❗✔️:     errout << "Found " << sGroupMap.size() << " SGroups when " << nSgroups
    // RDKit❗✔️:            << " were expected." << std::endl;
    // RDKit❗✔️:     if (strictParsing) {
    // RDKit❗✔️:       throw FileParseException(errout.str());
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // SGroups successfully parsed, now add them to the molecule
    // RDKit❗✔️:   for (const auto &sg : sGroupMap) {
    // RDKit❗✔️:     if (sg.second.getIsValid()) {
    // RDKit❗✔️:       addSubstanceGroup(*mol, sg.second);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << "SGroup " << sg.first
    // RDKit❗✔️:                               << " is invalid and will be ignored" << std::endl;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return tempStr;
    // RDKit❗✔️: }
    // RDKit❗✔️:
    // END RDKIT CPP FUNCTION
    // Three source formatted extractions consume only a numeric prefix for
    // unsigned fields; type consumes a whitespace-delimited string. Keep the
    // untouched suffix for the existing literal-space label parser. Reuse the
    // sole cursor's sign/overflow handling; malformed/overflow source state
    // remains a structural error, never a fabricated destination value.
    // Cost: linear scans of the three consumed fields, borrowed slices, no
    // line clone, new buffering or chemistry work. No all-input marker upgrade.
    let mut cursor = SGroupLineCursor::new(line);
    cursor.skip_c_locale_whitespace();
    let sequence_start = cursor.position;
    cursor.read_unsigned();
    let sequence = &line[sequence_start..cursor.position];
    parse_sgroup_header_unsigned(sequence, "sequence ID", line_number)?;
    cursor.skip_c_locale_whitespace();
    let kind_start = cursor.position;
    while cursor
        .peek()
        .is_some_and(|byte| !byte.is_ascii_whitespace())
    {
        cursor.position += 1;
    }
    let kind = &line[kind_start..cursor.position];
    cursor.skip_c_locale_whitespace();
    let external_start = cursor.position;
    cursor.read_unsigned();
    let external_id = &line[external_start..cursor.position];
    if kind.is_empty() || external_id.is_empty() {
        return Err(SdfReadError::Parse(format!(
            "SGroup line too short: '{line}' on line {line_number}"
        )));
    }
    parse_sgroup_header_unsigned(external_id, "external ID", line_number)?;
    Ok((sequence, kind, external_id, &line[cursor.position..]))
}

fn parse_sgroup_header_unsigned(
    text: &str,
    field: &str,
    line_number: usize,
) -> Result<u32, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock (formatted header fields)
    // RDKit❗✔️:     unsigned int sequenceId;
    // RDKit❗✔️:     unsigned int externalId;
    // RDKit❗✔️:     std::string type;
    // RDKit❗✔️:
    // RDKit❗✔️:     std::stringstream lineStream(tempStr);
    // RDKit❗✔️:     lineStream >> sequenceId;
    // RDKit❗✔️:     lineStream >> type;
    // RDKit❗✔️:     lineStream >> externalId;
    let (negative, digits) = match text.as_bytes().first() {
        Some(b'+') => (false, &text[1..]),
        Some(b'-') => (true, &text[1..]),
        _ => (false, text),
    };
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(SdfReadError::Parse(format!(
            "Cannot convert '{text}' to unsigned {field} on line {line_number}"
        )));
    }
    let magnitude = digits.parse::<u32>().map_err(|_| {
        SdfReadError::Parse(format!(
            "Unsigned {field} '{text}' is out of range on line {line_number}"
        ))
    })?;
    // Formatted unsigned extraction accepts a sign. For a representable
    // magnitude, a leading minus is assigned modulo the destination width;
    // fixed RDKit 2026.03.1 therefore observes `-1` as `UINT_MAX`.
    let value = if negative {
        0_u32.wrapping_sub(magnitude)
    } else {
        magnitude
    };
    // Behavioral review: complete signed decimal tokens match the pinned
    // formatted `unsigned int` assignment, including `+` and negative wrap.
    // The source leaves these destinations uninitialized on no conversion and
    // subsequently reads them, and overflow also poisons later extraction;
    // those C++ undefined/uninitialized cases are rejected structurally here
    // rather than assigned a fabricated zero. This boundary is not marked as
    // all-input source equivalence.
    // Complexity review: one linear validation scan and one bounded decimal
    // conversion match the source formatted extraction asymptotically, with
    // no allocation beyond an error message on failure.
    Ok(value)
    // END RDKIT CPP FUNCTION
}

#[derive(Clone, Copy)]
struct SGroupLineCursor<'a> {
    text: &'a str,
    position: usize,
    failed: bool,
}

impl<'a> SGroupLineCursor<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            text,
            position: 0,
            failed: false,
        }
    }

    fn get(&mut self) -> Option<u8> {
        if self.failed {
            return None;
        }
        let byte = self.text.as_bytes().get(self.position).copied()?;
        self.position += 1;
        Some(byte)
    }

    fn peek(&self) -> Option<u8> {
        (!self.failed)
            .then(|| self.text.as_bytes().get(self.position).copied())
            .flatten()
    }

    fn skip_c_locale_whitespace(&mut self) {
        while self
            .peek()
            .is_some_and(|byte| matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c))
        {
            self.position += 1;
        }
    }

    fn read_unsigned(&mut self) -> Option<u32> {
        if self.failed {
            return None;
        }
        self.skip_c_locale_whitespace();
        let negative = self.peek() == Some(b'-');
        if matches!(self.peek(), Some(b'+') | Some(b'-')) {
            self.position += 1;
        }
        let start = self.position;
        let mut value = 0_u32;
        let mut overflow = false;
        while let Some(byte) = self.peek().filter(u8::is_ascii_digit) {
            overflow |= value
                .checked_mul(10)
                .and_then(|current| current.checked_add(u32::from(byte - b'0')))
                .map_or(true, |next| {
                    value = next;
                    false
                });
            self.position += 1;
        }
        if self.position == start {
            self.failed = true;
            return Some(0);
        }
        if overflow {
            self.failed = true;
            return Some(u32::MAX);
        }
        Some(if negative {
            value.wrapping_neg()
        } else {
            value
        })
    }

    fn read_double(&mut self) -> Option<f64> {
        // BEGIN LIBSTDC++ CPP FUNCTION num_get::_M_extract_float (C locale)
        // libstdc++✔️✔️: else if ((__c == __lit[__num_base::_S_ie]
        // libstdc++✔️✔️:           || __c == __lit[__num_base::_S_iE])
        // libstdc++✔️✔️:          && !__found_sci && __found_mantissa)
        // libstdc++✔️✔️:   {
        // libstdc++✔️✔️:     // Scientific notation.
        // libstdc++✔️✔️:     __xtrc += 'e';
        // libstdc++✔️✔️:     __found_sci = true;
        // libstdc++✔️✔️:
        // libstdc++✔️✔️:     // Remove optional plus or minus sign, if they exist.
        // libstdc++✔️✔️:     if (++__beg != __end)
        // libstdc++✔️✔️:       {
        // libstdc++✔️✔️:         __c = *__beg;
        // libstdc++✔️✔️:         const bool __plus = __c == __lit[__num_base::_S_iplus];
        // libstdc++✔️✔️:         if (__plus || __c == __lit[__num_base::_S_iminus])
        // libstdc++✔️✔️:           __xtrc += __plus ? '+' : '-';
        // libstdc++✔️✔️:         else
        // libstdc++✔️✔️:           continue;
        // libstdc++✔️✔️:       }
        // libstdc++✔️✔️:     else
        // libstdc++✔️✔️:       {
        // libstdc++✔️✔️:         __testeof = true;
        // libstdc++✔️✔️:         break;
        // libstdc++✔️✔️:       }
        // libstdc++✔️✔️:   }
        // END LIBSTDC++ CPP FUNCTION num_get::_M_extract_float (C locale)
        // BEGIN LIBSTDC++ CPP FUNCTION num_get::do_get(double&)
        // libstdc++✔️✔️: string __xtrc;
        // libstdc++✔️✔️: __xtrc.reserve(32);
        // libstdc++✔️✔️: __beg = _M_extract_float(__beg, __end, __io, __err, __xtrc);
        // libstdc++✔️✔️: std::__convert_to_v(__xtrc.c_str(), __v, __err, _S_get_c_locale());
        // libstdc++✔️✔️: if (__beg == __end)
        // libstdc++✔️✔️:   __err |= ios_base::eofbit;
        // libstdc++✔️✔️: return __beg;
        // END LIBSTDC++ CPP FUNCTION num_get::do_get(double&)
        // BEGIN LIBSTDC++ CPP FUNCTION __convert_to_v(double&)
        // libstdc++✔️✔️: char* __sanity;
        // libstdc++✔️✔️: __v = __strtod_l(__s, &__sanity, __cloc);
        // libstdc++✔️✔️:
        // libstdc++✔️✔️: // _GLIBCXX_RESOLVE_LIB_DEFECTS
        // libstdc++✔️✔️: // 23. Num_get overflow result.
        // libstdc++✔️✔️: if (__sanity == __s || *__sanity != '\0')
        // libstdc++✔️✔️:   {
        // libstdc++✔️✔️:     __v = 0.0;
        // libstdc++✔️✔️:     __err = ios_base::failbit;
        // libstdc++✔️✔️:   }
        // libstdc++✔️✔️: else if (__v == numeric_limits<double>::infinity())
        // libstdc++✔️✔️:   {
        // libstdc++✔️✔️:     __v = numeric_limits<double>::max();
        // libstdc++✔️✔️:     __err = ios_base::failbit;
        // libstdc++✔️✔️:   }
        // libstdc++✔️✔️: else if (__v == -numeric_limits<double>::infinity())
        // libstdc++✔️✔️:   {
        // libstdc++✔️✔️:     __v = -numeric_limits<double>::max();
        // libstdc++✔️✔️:     __err = ios_base::failbit;
        // libstdc++✔️✔️:   }
        // END LIBSTDC++ CPP FUNCTION __convert_to_v(double&)
        // Behavioral review: the scanner consumes an exponent marker and its
        // optional sign even when no exponent digit follows. The fixed GCC
        // 15.2 libstdc++ conversion then assigns zero and failbit; decimal
        // overflow assigns the signed finite maximum and failbit. Once failed,
        // the cursor models the stream sentry by refusing subsequent reads.
        // The existing raw-atof owner remains responsible for supported finite
        // decimal conversion. Exceptionally long significands, subnormal
        // boundaries, hexadecimal floats, and NaN payloads retain the approved
        // non-equivalence boundary and are not claimed here.
        // Complexity review: this is one allocation-free linear byte scan plus
        // the existing linear decimal conversion, matching formatted extraction
        // asymptotically without repeated scans or new buffering.
        if self.failed {
            return None;
        }
        self.skip_c_locale_whitespace();
        let start = self.position;
        if matches!(self.peek(), Some(b'+') | Some(b'-')) {
            self.position += 1;
        }
        let integer_start = self.position;
        while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
            self.position += 1;
        }
        let mut has_digit = self.position != integer_start;
        if self.peek() == Some(b'.') {
            self.position += 1;
            let fraction_start = self.position;
            while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                self.position += 1;
            }
            has_digit |= self.position != fraction_start;
        }
        if !has_digit {
            self.failed = true;
            return Some(0.0);
        }
        if matches!(self.peek(), Some(b'e') | Some(b'E')) {
            self.position += 1;
            if matches!(self.peek(), Some(b'+') | Some(b'-')) {
                self.position += 1;
            }
            let exponent_digits = self.position;
            while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                self.position += 1;
            }
            if self.position == exponent_digits {
                self.failed = true;
                return Some(0.0);
            }
        }
        let value = parse_rdkit_atof(&self.text[start..self.position]);
        if value.is_infinite() {
            self.failed = true;
            Some(value.signum() * f64::MAX)
        } else {
            Some(value)
        }
    }

    fn read_string_token(&mut self) -> Option<&'a str> {
        if self.failed {
            return None;
        }
        self.skip_c_locale_whitespace();
        let start = self.position;
        while self
            .peek()
            .is_some_and(|byte| !matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c))
        {
            self.position += 1;
        }
        if self.position == start {
            self.failed = true;
            None
        } else {
            Some(&self.text[start..self.position])
        }
    }

    fn read_label(&mut self) -> Option<&'a str> {
        if self.failed || self.position == self.text.len() {
            return None;
        }
        let start = self.position;
        while let Some(byte) = self.get() {
            if byte == b'=' {
                return Some(&self.text[start..self.position - 1]);
            }
        }
        (start != self.position).then_some(&self.text[start..self.position])
    }
}

fn parse_array<T: Copy>(
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    max_count: Option<usize>,
    strict_parsing: bool,
    mut parse: impl FnMut(&mut SGroupLineCursor<'_>) -> Option<T>,
) -> Result<Vec<T>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000Array
    // RDKit✔️✔️:   auto paren = stream.get();  // discard parentheses
    // RDKit✔️✔️:   if (paren != '(') {
    // RDKit✔️✔️:   BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "WARNING: first character of V3000 array is not '('" << std::endl;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:   unsigned int count = 0;
    // RDKit✔️✔️:   stream >> count;
    // RDKit✔️✔️:   std::vector<T> values;
    // RDKit✔️✔️:   if (maxV >= 0 && count > static_cast<unsigned int>(maxV)) {
    // RDKit✔️✔️:     SGroupWarnOrThrow(strictParsing, "invalid count value");
    // RDKit✔️✔️:   return values;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:   values.reserve(count);
    // RDKit✔️✔️:     T value;
    // RDKit✔️✔️:   for (unsigned i = 0; i < count; ++i) {
    // RDKit✔️✔️:     stream >> value;
    // RDKit✔️✔️:     values.push_back(value);
    // RDKit✔️✔️: }
    // RDKit✔️✔️:   paren = stream.get();  // discard parentheses
    // RDKit✔️✔️:   if (paren != ')') {
    // RDKit✔️✔️:   BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "WARNING: final character of V3000 array is not ')'" << std::endl;
    // RDKit✔️✔️: }
    // Missing parentheses are warnings only. The cursor deliberately consumes
    // the same bytes as the source stream. Formatted extraction assigns zero
    // to an arithmetic destination when no characters can be converted and
    // sets failbit; later sentry failures leave that zero in place. Overflow
    // assigns the corresponding saturated boundary and also sets failbit.
    let _opening_parenthesis = cursor.get();
    let count = cursor.read_unsigned().unwrap_or(0) as usize;
    if max_count.is_some_and(|maximum| count > maximum) {
        return if strict_parsing {
            Err(SdfReadError::Parse("invalid count value".to_owned()))
        } else {
            Ok(Vec::new())
        };
    }
    let mut values = Vec::with_capacity(count);
    let mut previous = None;
    for _ in 0..count {
        let value = parse(cursor).or(previous).ok_or_else(|| {
            SdfReadError::Parse(format!(
                "V3000 array has no initialized value on line {line_number}"
            ))
        })?;
        values.push(value);
        previous = Some(value);
    }
    let _closing_parenthesis = cursor.get();
    Ok(values)
    // END RDKIT CPP FUNCTION
}

fn parse_u32_array(
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    max_count: Option<usize>,
    strict_parsing: bool,
) -> Result<Vec<u32>, SdfReadError> {
    parse_array(cursor, line_number, max_count, strict_parsing, |cursor| {
        cursor.read_unsigned()
    })
}

fn parse_f64_array(
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    max_count: Option<usize>,
    strict_parsing: bool,
) -> Result<Vec<f64>, SdfReadError> {
    parse_array(cursor, line_number, max_count, strict_parsing, |cursor| {
        cursor.read_double()
    })
}

fn parse_string_property(cursor: &mut SGroupLineCursor<'_>) -> String {
    // BEGIN RDKIT CPP FUNCTION ParseV3000StringPropLabel
    // RDKit✔️✔️: std::string ParseV3000StringPropLabel(std::stringstream &stream) {
    // RDKit✔️✔️:   std::string strValue;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto nextChar = stream.peek();
    // RDKit✔️✔️:   if (nextChar == ' ') {
    // RDKit✔️✔️:     // empty value, we peeked at the next field's separator
    // RDKit✔️✔️:     return strValue;
    // RDKit✔️✔️:   } else if (nextChar == '"') {
    // RDKit✔️✔️:     // skip the opening quote:
    // RDKit✔️✔️:     stream.get();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // this is a bit gross because it's legal to include a \" in a value,
    // RDKit✔️✔️:     // but the way that's done is by doubling it. So
    // RDKit✔️✔️:     // FIELDINFO=""""
    // RDKit✔️✔️:     // should assign the value \" to FIELDINFO
    // RDKit✔️✔️:     char chr;
    // RDKit✔️✔️:     while (stream.get(chr)) {
    // RDKit✔️✔️:       if (chr == '"') {
    // RDKit✔️✔️:         nextChar = stream.peek();
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // if the next element in the stream is a \" then we have a quoted \".
    // RDKit✔️✔️:         // Otherwise we're done
    // RDKit✔️✔️:         if (nextChar != '"') {
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           // skip the second \"
    // RDKit✔️✔️:           stream.get();
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       strValue += chr;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (nextChar == '\'') {
    // RDKit✔️✔️:     std::getline(stream, strValue, '\'');
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     stream >> strValue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   boost::trim_right(strValue);
    // RDKit✔️✔️:   return strValue;
    // RDKit✔️✔️: }
    let mut value = Vec::new();
    match cursor.peek() {
        Some(b' ') => {}
        Some(b'"') => {
            cursor.get();
            while let Some(byte) = cursor.get() {
                if byte == b'"' {
                    if cursor.peek() != Some(b'"') {
                        break;
                    }
                    cursor.get();
                }
                value.push(byte);
            }
        }
        Some(b'\'') => {
            // `std::getline(stream, value, '\'')` starts on the opening quote,
            // so the source returns an empty value and consumes that quote.
            cursor.get();
        }
        Some(_) => {
            cursor.skip_c_locale_whitespace();
            while let Some(byte) = cursor.peek() {
                if matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c) {
                    break;
                }
                value.push(byte);
                cursor.position += 1;
            }
        }
        None => {}
    }
    while value.last().is_some_and(|byte| byte.is_ascii_whitespace()) {
        value.pop();
    }
    String::from_utf8(value)
        .expect("SGroup cursor removes only ASCII delimiters from valid UTF-8 input")
    // END RDKIT CPP FUNCTION
}

fn parse_cstate(
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    group: &mut SubstanceGroup,
    bonds: &BTreeMap<u32, BondId>,
    bond_endpoints: &dyn Fn(BondId) -> Option<(AtomId, AtomId)>,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000CStateLabel
    // RDKit✔️✔️:   stream.get();  // discard parentheses
    // RDKit✔️✔️:   unsigned int count;
    // RDKit✔️✔️:   unsigned int bondMark;
    // RDKit✔️✔️:   stream >> count >> bondMark;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string type = sgroup.getProp<std::string>("TYPE");
    // RDKit✔️✔️:   if ((type != "SUP" && count != 1) || (type == "SUP" && count != 4)) {
    // RDKit✔️✔️:   std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Unexpected number of fields for CSTATE field on line " << line;
    // RDKit✔️✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
    // RDKit✔️✔️:     sgroup.setIsValid(false);
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:   Bond *bond = mol->getUniqueBondWithBookmark(bondMark);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     RDGeom::Point3D vector;
    // RDKit✔️✔️:   if (type == "SUP") {
    // RDKit✔️✔️:     stream >> vector.x >> vector.y >> vector.z;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:     sgroup.addCState(bond->getIdx(), vector);
    // RDKit✔️✔️:   } catch (const std::exception &e) {
    // RDKit✔️✔️:       SGroupWarnOrThrow<>(strictParsing, e.what());
    // RDKit✔️✔️:     sgroup.setIsValid(false);
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   stream.get();  // discard final parentheses
    // END RDKIT CPP FUNCTION
    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::addCState
    // RDKit✔️✔️: void SubstanceGroup::addCState(unsigned int bondIdx,
    // RDKit✔️✔️:                                const RDGeom::Point3D &vector) {
    // RDKit✔️✔️:   PRECONDITION(dp_mol, "bad mol");
    // RDKit✔️✔️:   PRECONDITION(!d_bonds.empty(), "no bonds");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (getBondType(bondIdx) != SubstanceGroup::BondType::XBOND) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Bond with index " << bondIdx
    // RDKit✔️✔️:            << " is not an XBOND for current SubstanceGroup";
    // RDKit✔️✔️:     throw SubstanceGroupException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   d_cstates.push_back({bondIdx, vector});
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::getBondType
    // RDKit✔️✔️: SubstanceGroup::BondType SubstanceGroup::getBondType(
    // RDKit✔️✔️:     unsigned int bondIdx) const {
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       std::find(d_bonds.begin(), d_bonds.end(), bondIdx) != d_bonds.end(),
    // RDKit✔️✔️:       "bond is not part of the SubstanceGroup")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto bond = dp_mol->getBondWithIdx(bondIdx);
    // RDKit✔️✔️:   bool begin_atom_in_sgroup =
    // RDKit✔️✔️:       std::find(d_atoms.begin(), d_atoms.end(), bond->getBeginAtomIdx()) !=
    // RDKit✔️✔️:       d_atoms.end();
    // RDKit✔️✔️:   bool end_atom_in_sgroup = std::find(d_atoms.begin(), d_atoms.end(),
    // RDKit✔️✔️:                                       bond->getEndAtomIdx()) != d_atoms.end();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (begin_atom_in_sgroup && end_atom_in_sgroup) {
    // RDKit✔️✔️:     return SubstanceGroup::BondType::CBOND;
    // RDKit✔️✔️:   } else if (begin_atom_in_sgroup || end_atom_in_sgroup) {
    // RDKit✔️✔️:     return SubstanceGroup::BondType::XBOND;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Neither beginning nor ending atoms of bond " << bond->getIdx()
    // RDKit✔️✔️:            << " is in this SubstanceGroup.";
    // RDKit✔️✔️:     throw SubstanceGroupException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION

    let _opening_parenthesis = cursor.get();
    let count = cursor.read_unsigned().ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Unexpected number of fields for CSTATE field on line {line_number}"
        ))
    })?;
    let expected = if group.kind() == &SubstanceGroupKind::Superatom {
        4
    } else {
        1
    };
    if count != expected {
        return Err(SdfReadError::Parse(format!(
            "Unexpected number of fields for CSTATE field on line {line_number}"
        )));
    }
    let bookmark = cursor.read_unsigned().ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Unexpected number of fields for CSTATE field on line {line_number}"
        ))
    })?;
    let bond = *bonds.get(&bookmark).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup bond index {bookmark} out of range on line {line_number}"
        ))
    })?;
    let mut vector = [0.0; 3];
    if group.kind() == &SubstanceGroupKind::Superatom {
        // Point3D is zero-initialized. Formatted extraction stops after its
        // first failed component, leaving that and all later components zero.
        for component in &mut vector {
            let Some(parsed) = cursor.read_double() else {
                break;
            };
            *component = parsed;
        }
    }
    let (begin, end) = bond_endpoints(bond).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup bond index {bookmark} out of range on line {line_number}"
        ))
    })?;
    let begin_is_member = group.atoms().contains(&begin);
    let end_is_member = group.atoms().contains(&end);
    if !group.bonds().contains(&bond) || begin_is_member == end_is_member {
        return Err(SdfReadError::Parse(format!(
            "Bond with index {} is not an XBOND for current SubstanceGroup",
            bond.index()
        )));
    }
    group.push_cstate(SGroupCState { bond, vector });
    let _closing_parenthesis = cursor.get();
    // Behavior review: the parser retains all XYZ components for SUP groups,
    // uses Point3D's zero vector for every other group kind, resolves the
    // externally unique V3000 bookmark to the canonical BondId, and admits
    // only the group's typed crossing-bond membership. Missing/invalid SUP
    // components reproduce formatted-stream fail-state propagation rather
    // than being parsed independently. A missing count/bookmark is rejected
    // structurally because the source reads an uninitialized destination in
    // that malformed case. Complexity review: this is one bounded scalar
    // scan plus one ordered-map lookup and one linear membership lookup,
    // matching the source helper chain without a second representation,
    // whole-table clone, or repeated graph scan.
    Ok(())
}

fn parse_sap(
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    group: &mut SubstanceGroup,
    atoms: &BTreeMap<u32, AtomId>,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SAPLabel
    // RDKit✔️✔️:   stream.get();  // discard parentheses
    // RDKit✔️✔️:   unsigned int count = 0;
    // RDKit✔️✔️:   unsigned int aIdxMark = 0;
    // RDKit✔️✔️:   std::string lvIdxStr;  // In V3000 this may be a string
    // RDKit✔️✔️:   std::string sapIdStr;
    // RDKit✔️✔️:   stream >> count >> aIdxMark >> lvIdxStr >> sapIdStr;
    // RDKit✔️✔️:   sapIdStr.pop_back();
    // RDKit✔️✔️:     unsigned int aIdx = mol->getAtomWithBookmark(aIdxMark)->getIdx();
    // RDKit✔️✔️:     int lvIdx = -1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   boost::to_upper(lvIdxStr);
    // RDKit✔️✔️:   if (lvIdxStr == "AIDX") {
    // RDKit✔️✔️:     lvIdx = aIdx;
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:     unsigned int lvIdxTmp = FileParserUtils::toInt(lvIdxStr);
    // RDKit✔️✔️:     if (lvIdxTmp > 0) {
    // RDKit✔️✔️:       lvIdx = mol->getAtomWithBookmark(lvIdxTmp)->getIdx();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:     sgroup.addAttachPoint(aIdx, lvIdx, sapIdStr);
    // RDKit✔️✔️:   } catch (const std::exception &e) {
    // RDKit✔️✔️:       SGroupWarnOrThrow<>(strictParsing, e.what());
    // RDKit✔️✔️:     sgroup.setIsValid(false);
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    let _discarded_opening = cursor.get();
    let _count = cursor.read_unsigned().ok_or_else(|| {
        SdfReadError::Parse(format!("SGroup SAP count missing on line {line_number}"))
    })?;
    let atom_bookmark = cursor.read_unsigned().ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup attach atom index missing on line {line_number}"
        ))
    })?;
    let leaving_text = cursor.read_string_token().ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup leaving atom index missing on line {line_number}"
        ))
    })?;
    let label_with_final_byte = cursor.read_string_token().ok_or_else(|| {
        SdfReadError::Parse(format!("SGroup SAP label missing on line {line_number}"))
    })?;
    let label_end = label_with_final_byte.len().checked_sub(1).ok_or_else(|| {
        SdfReadError::Parse(format!("SGroup SAP label missing on line {line_number}"))
    })?;
    let label = label_with_final_byte.get(..label_end).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup SAP label has invalid UTF-8 byte truncation on line {line_number}"
        ))
    })?;
    let atom = *atoms.get(&atom_bookmark).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup attach atom index {atom_bookmark} out of range on line {line_number}"
        ))
    })?;
    let leaving_atom = if leaving_text.eq_ignore_ascii_case("AIDX") {
        Some(atom)
    } else {
        let bookmark = parse_rdkit_int(leaving_text).map_err(|()| {
            SdfReadError::Parse(format!(
                "Cannot convert '{leaving_text}' to int on line {line_number}"
            ))
        })? as u32;
        if bookmark == 0 {
            None
        } else {
            Some(*atoms.get(&bookmark).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup leaving atom index {bookmark} out of range on line {line_number}"
                ))
            })?)
        }
    };
    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::addAttachPoint
    // RDKit✔️✔️: void SubstanceGroup::addAttachPoint(unsigned int aIdx, int lvIdx,
    // RDKit✔️✔️:                                     const std::string &idStr) {
    // RDKit✔️✔️:   d_saps.push_back({aIdx, lvIdx, idStr});
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SubstanceGroup::addAttachPoint
    group.push_attach_point(SGroupAttachPoint {
        atom,
        leaving_atom,
        label: Some(label.into()),
        order: None,
    });
    // Behavior review: the count is source-read but deliberately not checked;
    // atom bookmarks use formatted unsigned extraction, AIDX is uppercased by
    // the source, and numeric leaving bookmarks use the distinct toInt
    // contract (so leading '+' is no-conversion zero). Both nonzero bookmarks
    // must resolve in the molecule, but addAttachPoint imposes no ATOMS-member
    // precondition. The final byte of the label token is removed regardless of
    // whether it is ')'; invalid UTF-8 byte truncation is a structured Rust
    // boundary error instead of constructing an invalid string.
    // Complexity review: this advances the existing shared cursor once, does
    // two ordered bookmark lookups at most, and appends one typed value. It is
    // linear in the four token lengths with no token vector or parallel model,
    // matching the source asymptotic and allocation shape.
    Ok(())
    // END RDKIT CPP FUNCTION
}

fn parse_label(
    label: &str,
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    group: &mut SubstanceGroup,
    parents: &mut BTreeMap<u32, u32>,
    atoms: &BTreeMap<u32, AtomId>,
    bonds: &BTreeMap<u32, BondId>,
    bond_endpoints: &dyn Fn(BondId) -> Option<(AtomId, AtomId)>,
    strict_parsing: bool,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000ParseLabel
    // RDKit✔️✔️:     if (label == "XBHEAD" || label == "XBCORR") {
    // RDKit✔️✔️:       std::vector<unsigned int> bvect = ParseV3000Array<unsigned int>(
    // RDKit✔️✔️:           lineStream, mol->getNumBonds(), strictParsing);
    // RDKit✔️✔️:       std::transform(bvect.begin(), bvect.end(), bvect.begin(),
    // RDKit✔️✔️:                      [](unsigned int v) -> unsigned int { return v - 1; });
    // RDKit✔️✔️:       sgroup.setProp(label, bvect);
    // RDKit✔️✔️:     } else if (label == "ATOMS") {
    // RDKit❗✔️:       for (auto atomIdx : ParseV3000Array<unsigned int>(
    // RDKit❗✔️:                lineStream, mol->getNumAtoms(), strictParsing)) {
    // RDKit❗✔️:         sgroup.addAtomWithBookmark(atomIdx);
    // RDKit❗✔️:   }
    // RDKit❗✔️:     } else if (label == "PATOMS") {
    // RDKit❗✔️:       for (auto patomIdx : ParseV3000Array<unsigned int>(
    // RDKit❗✔️:                lineStream, mol->getNumAtoms(), strictParsing)) {
    // RDKit❗✔️:         sgroup.addParentAtomWithBookmark(patomIdx);
    // RDKit❗✔️:   }
    // RDKit❗✔️:     } else if (label == "CBONDS" || label == "XBONDS") {
    // RDKit❗✔️:       for (auto bondIdx : ParseV3000Array<unsigned int>(
    // RDKit❗✔️:                lineStream, mol->getNumBonds(), strictParsing)) {
    // RDKit❗✔️:         sgroup.addBondWithBookmark(bondIdx);
    // RDKit❗✔️:   }
    // RDKit✔️✔️:     } else if (label == "BRKXYZ") {
    // RDKit✔️✔️:       auto coords = ParseV3000Array<double>(lineStream, 9, strictParsing);
    // RDKit✔️✔️:       if (coords.size() != 9) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:         errout << "Unexpected number of coordinates for BRKXYZ on line "
    // RDKit✔️✔️:            << line;
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   SubstanceGroup::Bracket bracket;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < 3; ++i) {
    // RDKit✔️✔️:         bracket[i] = RDGeom::Point3D(*(coords.begin() + (3 * i)),
    // RDKit✔️✔️:                                      *(coords.begin() + (3 * i) + 1),
    // RDKit✔️✔️:                                      *(coords.begin() + (3 * i) + 2));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:     sgroup.addBracket(bracket);
    // RDKit✔️✔️:     } else if (label == "CSTATE") {
    // RDKit❗✔️:       ParseV3000CStateLabel(mol, sgroup, lineStream, line, strictParsing);
    // RDKit❗✔️:     } else if (label == "SAP") {
    // RDKit❗✔️:       ParseV3000SAPLabel(mol, sgroup, lineStream, strictParsing);
    // RDKit❗✔️:     } else if (label == "PARENT") {
    // RDKit❗✔️:       lineStream >> parentIdx;
    // RDKit❗✔️:       sgroup.setProp<unsigned int>("PARENT", parentIdx);
    // RDKit❗✔️:     } else if (label == "COMPNO") {
    // RDKit❗✔️:       lineStream >> compno;
    // RDKit❗✔️:     if (compno > 256u) {
    // RDKit❗✔️:     std::ostringstream errout;
    // RDKit❗✔️:       errout << "SGroup SNC value over 256: '" << compno << "' on line "
    // RDKit❗✔️:            << line;
    // RDKit❗✔️:     throw FileParseException(errout.str());
    // RDKit❗✔️:   }
    // RDKit❗✔️:       sgroup.setProp<unsigned int>("COMPNO", compno);
    // RDKit❗✔️:     } else if (label == "FIELDDATA") {
    // RDKit❗✔️:       auto strValue = ParseV3000StringPropLabel(lineStream);
    // RDKit❗✔️:   if (strictParsing) {
    // RDKit❗✔️:         strValue = strValue.substr(0, 200);
    // RDKit❗✔️:   }
    // RDKit❗✔️:       dataFields.push_back(strValue);
    // RDKit❗✔️: } else {
    // RDKit❗✔️:       auto strValue = ParseV3000StringPropLabel(lineStream);
    // RDKit❗✔️:       sgroup.setProp(label, strValue);
    // RDKit❗✔️: }
    match label {
        "XBHEAD" | "XBCORR" => {
            // Unlike XBONDS/CBONDS, these values are one-based bond-table
            // positions, not V3000 bond bookmarks. The canonical model stores
            // the corresponding zero-based BondId so remapping and writers use
            // the same typed reference. RDKit's unsigned `v - 1` can retain an
            // invalid wrapped/out-of-range property value; COSMolKit rejects
            // that value at the detached-model boundary because persistent
            // typed references must satisfy local topology invariants.
            for source_position in
                parse_u32_array(cursor, line_number, Some(bonds.len()), strict_parsing)?
            {
                let index = source_position.checked_sub(1).ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "SGroup {label} bond-row position 0 is invalid on line {line_number}"
                    ))
                })? as usize;
                if index >= bonds.len() {
                    return Err(SdfReadError::Parse(format!(
                        "SGroup {label} bond-row position {source_position} out of range on line {line_number}"
                    )));
                }
                let bond = BondId::new(index);
                if label == "XBHEAD" {
                    group.push_head_crossing_bond(bond);
                } else {
                    group.push_crossing_bond_correspondence(bond);
                }
            }
        }
        "ATOMS" | "PATOMS" => {
            for bookmark in parse_u32_array(cursor, line_number, Some(atoms.len()), strict_parsing)?
            {
                let atom = *atoms.get(&bookmark).ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "SGroup atom index {bookmark} out of range on line {line_number}"
                    ))
                })?;
                if label == "ATOMS" {
                    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::addAtomWithBookmark
                    // RDKit✔️✔️: void SubstanceGroup::addAtomWithBookmark(int mark) {
                    // RDKit✔️✔️:   PRECONDITION(dp_mol, "bad mol");
                    // RDKit✔️✔️:   Atom *atom = dp_mol->getUniqueAtomWithBookmark(mark);
                    // RDKit✔️✔️:   PRECONDITION(atom, "atom not found");
                    // RDKit✔️✔️:   d_atoms.push_back(atom->getIdx());
                    // RDKit✔️✔️: }
                    group.push_atom(atom);
                    // END RDKIT CPP FUNCTION
                } else {
                    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::addParentAtomWithBookmark
                    // RDKit✔️✔️: void SubstanceGroup::addParentAtomWithBookmark(int mark) {
                    // RDKit✔️✔️:   PRECONDITION(dp_mol, "bad mol");
                    // RDKit✔️✔️:
                    // RDKit✔️✔️:   Atom *atom = dp_mol->getUniqueAtomWithBookmark(mark);
                    // RDKit✔️✔️:   unsigned int idx = atom->getIdx();
                    // RDKit✔️✔️:   if (std::find(d_atoms.begin(), d_atoms.end(), idx) == d_atoms.end()) {
                    // RDKit✔️✔️:     std::ostringstream errout;
                    // RDKit✔️✔️:     errout << "Atom with bookmark " << mark
                    // RDKit✔️✔️:            << " is not a member of current SubstanceGroup ";
                    // RDKit✔️✔️:     throw SubstanceGroupException(errout.str());
                    // RDKit✔️✔️:   }
                    // RDKit✔️✔️:
                    // RDKit✔️✔️:   d_patoms.push_back(idx);
                    // RDKit✔️✔️: }
                    if !group.atoms().contains(&atom) {
                        return Err(SdfReadError::Parse(format!(
                            "Atom with bookmark {bookmark} is not a member of current SubstanceGroup on line {line_number}"
                        )));
                    }
                    group.push_parent_atom(atom);
                    // END RDKIT CPP FUNCTION
                }
            }
        }
        "CBONDS" | "XBONDS" => {
            let role = if label == "CBONDS" {
                SGroupBondRole::Contained
            } else {
                SGroupBondRole::Crossing
            };
            for bookmark in parse_u32_array(cursor, line_number, Some(bonds.len()), strict_parsing)?
            {
                let bond = *bonds.get(&bookmark).ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "SGroup bond index {bookmark} out of range on line {line_number}"
                    ))
                })?;
                // BEGIN RDKIT CPP FUNCTION SubstanceGroup::addBondWithBookmark
                // RDKit✔️✔️: void SubstanceGroup::addBondWithBookmark(int mark) {
                // RDKit✔️✔️:   PRECONDITION(dp_mol, "bad mol");
                // RDKit✔️✔️:   Bond *bond = dp_mol->getUniqueBondWithBookmark(mark);
                // RDKit✔️✔️:   d_bonds.push_back(bond->getIdx());
                // RDKit✔️✔️: }
                group.push_bond_with_role(bond, role);
                // END RDKIT CPP FUNCTION
            }
        }
        "BRKXYZ" => {
            let coordinates = parse_f64_array(cursor, line_number, Some(9), strict_parsing)?;
            if coordinates.len() != 9 {
                return Err(SdfReadError::Parse(format!(
                    "Unexpected number of coordinates for BRKXYZ on line {line_number}"
                )));
            }
            group.display_mut().brackets.push(SGroupBracket {
                points: [
                    [coordinates[0], coordinates[1], coordinates[2]],
                    [coordinates[3], coordinates[4], coordinates[5]],
                    [coordinates[6], coordinates[7], coordinates[8]],
                ],
            });
            // The fixed nine-value array is copied once into the canonical
            // three-point display value. This remains O(1), allocates only the
            // source-equivalent appended bracket slot, and never enters the
            // molecule coordinate/conformer block.
        }
        _ => {
            parse_non_array_label(
                label,
                cursor,
                line_number,
                group,
                parents,
                atoms,
                bonds,
                bond_endpoints,
                strict_parsing,
            )?;
        }
    }
    Ok(())
    // END RDKIT CPP FUNCTION
}

#[allow(clippy::too_many_arguments)]
fn parse_non_array_label(
    label: &str,
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    group: &mut SubstanceGroup,
    parents: &mut BTreeMap<u32, u32>,
    atoms: &BTreeMap<u32, AtomId>,
    bonds: &BTreeMap<u32, BondId>,
    bond_endpoints: &dyn Fn(BondId) -> Option<(AtomId, AtomId)>,
    strict_parsing: bool,
) -> Result<(), SdfReadError> {
    match label {
        "CSTATE" => {
            // The source passes the shared lineStream into the helper. In
            // particular, a formatted-extraction failbit inside CSTATE stops
            // the outer label loop; isolating the value in a second cursor
            // would incorrectly allow later labels to be consumed.
            parse_cstate(cursor, line_number, group, bonds, bond_endpoints)?;
        }
        "SAP" => {
            parse_sap(cursor, line_number, group, atoms)?;
        }
        "PARENT" => {
            // BEGIN RDKIT CPP FUNCTION ParseV3000ParseLabel (PARENT)
            // RDKit✔️✔️:     } else if (label == "PARENT") {
            // RDKit✔️✔️:       // Store relationship until all SGroups have been read
            // RDKit✔️✔️:       unsigned int parentIdx;
            // RDKit✔️✔️:       if (lineStream.eof()) {
            // RDKit✔️✔️:     std::ostringstream errout;
            // RDKit✔️✔️:         errout << "PARENT label not found on line " << line;
            // RDKit✔️✔️:     throw FileParseException(errout.str());
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:       lineStream >> parentIdx;
            // RDKit✔️✔️:       if (lineStream.fail()) {
            // RDKit✔️✔️:     std::ostringstream errout;
            // RDKit✔️✔️:         errout << "Invalid PARENT label found on line " << line;
            // RDKit✔️✔️:     throw FileParseException(errout.str());
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:       sgroup.setProp<unsigned int>("PARENT", parentIdx);
            if cursor.peek().is_none() {
                return Err(SdfReadError::Parse(format!(
                    "PARENT label not found on line {line_number}"
                )));
            }
            let parent = cursor.read_unsigned().ok_or_else(|| {
                SdfReadError::Parse(format!("Invalid PARENT label found on line {line_number}"))
            })?;
            if cursor.failed {
                return Err(SdfReadError::Parse(format!(
                    "Invalid PARENT label found on line {line_number}"
                )));
            }
            group.set_prop("PARENT", PropertyValue::UInt(parent))?;
            if let Some(sequence) = group.rdkit_sequence_id() {
                parents.insert(sequence, parent);
            }
            // Behavioral review: the shared cursor performs the source's
            // formatted unsigned extraction, including C-locale whitespace,
            // both signs, prefix consumption and fail state. The temporary
            // sequence relation is resolved only after all surviving rows are
            // known, so no raw or dangling PARENT property enters the typed
            // canonical model. Undefined/uninitialized overflow state remains
            // a structured error, matching the documented identity boundary.
            // Complexity review: one bounded decimal scan plus one ordered-map
            // insertion is O(token length + log S), with no duplicate model or
            // additional line scan; this is comparable to formatted extraction
            // followed by source property insertion.
            // END RDKIT CPP FUNCTION
        }
        "COMPNO" => {
            // BEGIN RDKIT CPP FUNCTION ParseV3000ParseLabel (COMPNO)
            // RDKit✔️✔️:     } else if (label == "COMPNO") {
            // RDKit✔️✔️:   unsigned int compno;
            // RDKit✔️✔️:       lineStream >> compno;
            // RDKit✔️✔️:     if (compno > 256u) {
            // RDKit✔️✔️:     std::ostringstream errout;
            // RDKit✔️✔️:       errout << "SGroup SNC value over 256: '" << compno << "' on line "
            // RDKit✔️✔️:            << line;
            // RDKit✔️✔️:     throw FileParseException(errout.str());
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:       sgroup.setProp<unsigned int>("COMPNO", compno);
            let number = cursor.read_unsigned().ok_or_else(|| {
                SdfReadError::Parse(format!("Invalid COMPNO label found on line {line_number}"))
            })?;
            if cursor.failed {
                return Err(SdfReadError::Parse(format!(
                    "Invalid COMPNO label found on line {line_number}"
                )));
            }
            if number > 256 {
                return Err(SdfReadError::Parse(format!(
                    "SGroup SNC value over 256: '{number}' on line {line_number}"
                )));
            }
            group.set_prop("COMPNO", PropertyValue::UInt(number))?;
            group.set_component_number(number);
            // Formatted unsigned extraction advances the shared label stream;
            // the bound check and one typed assignment are constant time and
            // allocate no intermediate token.
            // END RDKIT CPP FUNCTION
        }
        "FIELDDATA" => {
            // BEGIN RDKIT CPP FUNCTION ParseV3000ParseLabel (DAT fields)
            // RDKit❗✔️:     } else if (label == "FIELDDATA") {
            // RDKit❗✔️:       auto strValue = ParseV3000StringPropLabel(lineStream);
            // RDKit❗✔️:   if (strictParsing) {
            // RDKit❗✔️:         strValue = strValue.substr(0, 200);
            // RDKit❗✔️:   }
            // RDKit❗✔️:       dataFields.push_back(strValue);
            let parsed = parse_string_property(cursor);
            let parsed = PropertyText::from_bytes(
                &parsed.as_bytes()[..if strict_parsing {
                    parsed.len().min(200)
                } else {
                    parsed.len()
                }],
            );
            group.data_mut().values.push(parsed.into());
            // The canonical DAT value vector preserves append order. Strict
            // mode performs the source's 200-byte prefix operation; an invalid
            // UTF-8 cut is a structured Rust text-boundary error, so behavior
            // is intentionally partial at that representation boundary.
            // Non-strict mode moves the complete parsed string without a copy.
            // END RDKIT CPP FUNCTION
        }
        "SUBTYPE" => {
            let parsed = parse_string_property(cursor);
            if !is_valid_rdkit_sgroup_subtype(&parsed) {
                return Err(SdfReadError::Parse(format!(
                    "Unsupported SGroup subtype '{parsed}' on line {line_number}"
                )));
            }
            group.set_prop("SUBTYPE", parsed.clone())?;
            group.set_subtype(parsed);
        }
        "CONNECT" => {
            let parsed = parse_string_property(cursor);
            let connection = sgroup_connection_from_rdkit(&parsed).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "Unsupported SGroup connection type '{parsed}' on line {line_number}"
                ))
            })?;
            group.set_prop("CONNECT", parsed)?;
            group.set_connection(connection);
        }
        "CLASS" => {
            let parsed = parse_string_property(cursor);
            if !is_valid_rdkit_sgroup_class(&parsed) {
                return Err(SdfReadError::Parse(format!(
                    "Unsupported SGroup template class '{parsed}' on line {line_number}"
                )));
            }
            group.set_prop("CLASS", parsed.clone())?;
            group.set_class(parsed);
        }
        "LABEL" => {
            let value = parse_string_property(cursor);
            group.set_prop("LABEL", value.clone())?;
            group.set_label(value);
        }
        // BEGIN RDKIT CPP FUNCTION ParseV3000ParseLabel (string-property tail)
        // RDKit✔️✔️: } else {
        // RDKit✔️✔️:       // Parse string props
        // RDKit✔️✔️:       auto strValue = ParseV3000StringPropLabel(lineStream);
        // RDKit✔️✔️:
        // RDKit✔️✔️:       if (label == "SUBTYPE" &&
        // RDKit✔️✔️:           !SubstanceGroupChecks::isValidSubType(strValue)) {
        // RDKit✔️✔️:     std::ostringstream errout;
        // RDKit✔️✔️:         errout << "Unsupported SGroup subtype '" << strValue << "' on line "
        // RDKit✔️✔️:            << line;
        // RDKit✔️✔️:     throw FileParseException(errout.str());
        // RDKit✔️✔️:       } else if (label == "CONNECT" &&
        // RDKit✔️✔️:                  !SubstanceGroupChecks::isValidConnectType(strValue)) {
        // RDKit✔️✔️:     std::ostringstream errout;
        // RDKit✔️✔️:         errout << "Unsupported SGroup connection type '" << strValue
        // RDKit✔️✔️:              << "' on line " << line;
        // RDKit✔️✔️:     throw FileParseException(errout.str());
        // RDKit✔️✔️:       } else if (label == "CLASS" &&
        // RDKit✔️✔️:                  !SubstanceGroupChecks::isValidClass(strValue)) {
        // RDKit✔️✔️:     std::ostringstream errout;
        // RDKit✔️✔️:         errout << "Unsupported SGroup template class '" << strValue
        // RDKit✔️✔️:              << "' on line " << line;
        // RDKit✔️✔️:     throw FileParseException(errout.str());
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:       // NATREPLACE is not validated nor used
        // RDKit✔️✔️:
        // RDKit✔️✔️:       sgroup.setProp(label, strValue);
        // RDKit✔️✔️: }
        "FIELDNAME" => {
            let value = parse_string_property(cursor);
            group.set_prop("FIELDNAME", value.clone())?;
            group.data_mut().field_name = Some(value.into());
        }
        "FIELDTYPE" => {
            let value = parse_string_property(cursor);
            group.set_prop("FIELDTYPE", value.clone())?;
            group.data_mut().field_type = Some(value.into());
        }
        "FIELDINFO" => {
            let value = parse_string_property(cursor);
            group.set_prop("FIELDINFO", value.clone())?;
            group.data_mut().field_info = Some(value.into());
        }
        "FIELDDISP" => {
            let value = parse_string_property(cursor);
            group.set_prop("FIELDDISP", value.clone())?;
            group.data_mut().field_display = Some(value.into());
        }
        "QUERYTYPE" => {
            let value = parse_string_property(cursor);
            group.set_prop("QUERYTYPE", value.clone())?;
            group.data_mut().query_type = Some(value.into());
        }
        "QUERYOP" => {
            let value = parse_string_property(cursor);
            group.set_prop("QUERYOP", value.clone())?;
            group.data_mut().query_op = Some(value.into());
        }
        // The six DAT metadata labels use the one shared source string parser
        // and overwrite their single canonical typed slots in O(value length),
        // matching the source property assignment without a parallel raw prop.
        // END RDKIT CPP FUNCTION
        "ESTATE" => {
            let value = parse_string_property(cursor);
            group.set_prop("ESTATE", value.clone())?;
            group.set_expansion_state(value);
        }
        "BRKTYP" => {
            let parsed = parse_string_property(cursor);
            let style = match parsed.as_str() {
                "BRACKET" => SGroupBracketStyle::Bracket,
                "PAREN" => SGroupBracketStyle::Parenthesis,
                "" => SGroupBracketStyle::None,
                other => SGroupBracketStyle::Unknown(other.into()),
            };
            group.set_prop("BRKTYP", parsed)?;
            group.set_bracket_style(style);
        }
        other => group.set_prop(other, parse_string_property(cursor))?,
    }
    Ok(())
}

fn skip_overridden_default_value(
    cursor: &mut SGroupLineCursor<'_>,
    label: &str,
    line_number: usize,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock (overridden default)
    // RDKit✔️✔️:         spacer = lineStream.peek();
    // RDKit✔️✔️:         if (spacer == ' ') {
    // RDKit✔️✔️:   std::ostringstream errout;
    // RDKit✔️✔️:           errout << "Found unexpected whitespace at DEFAULT label " << label;
    // RDKit✔️✔️:   if (strictParsing) {
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️✔️:     sgroup.setIsValid(false);
    // RDKit✔️✔️:     continue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:         } else if (spacer == '(') {
    // RDKit✔️✔️:           std::getline(lineStream, label, ')');
    // RDKit✔️✔️:       lineStream.get(spacer);
    // RDKit✔️✔️:         } else if (spacer == '"') {
    // RDKit✔️✔️:       lineStream.get(spacer);
    // RDKit✔️✔️:           std::getline(lineStream, label, '"');
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:           std::getline(lineStream, label, ' ');
    // RDKit✔️✔️:           lineStream.putback(' ');
    // RDKit✔️✔️: }
    match cursor.peek() {
        Some(b' ') => {
            return Err(SdfReadError::Parse(format!(
                "Found unexpected whitespace at DEFAULT label {label} on line {line_number}"
            )));
        }
        Some(b'(') => {
            while cursor.get().is_some_and(|byte| byte != b')') {}
            let _ = cursor.get();
        }
        Some(b'"') => {
            let _ = cursor.get();
            while cursor.get().is_some_and(|byte| byte != b'"') {}
        }
        Some(_) => {
            while cursor.peek().is_some_and(|byte| byte != b' ') {
                cursor.position += 1;
            }
        }
        None => {}
    }
    // Behavioral review: this preserves the source's distinct skip grammar
    // instead of parsing the overridden value. A parenthesized value consumes
    // one byte after `)`, a quoted value stops at the first quote without
    // doubled-quote handling, and a scalar leaves its separating space for
    // the outer label loop. A literal-space empty value invalidates the group.
    // Complexity review: each skipped value is scanned once in place with no
    // allocation, matching the source's linear getline/stream operations.
    Ok(())
    // END RDKIT CPP FUNCTION
}

fn apply_labels(
    text: &str,
    line_number: usize,
    group: &mut SubstanceGroup,
    parents: &mut BTreeMap<u32, u32>,
    atoms: &BTreeMap<u32, AtomId>,
    bonds: &BTreeMap<u32, BondId>,
    bond_endpoints: &dyn Fn(BondId) -> Option<(AtomId, AtomId)>,
    seen: &mut BTreeSet<String>,
    defaults_only: bool,
    strict_parsing: bool,
) -> Result<bool, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock (label loop)
    // RDKit❗✔️:     while (sgroup.getIsValid() && !lineStream.eof() && !lineStream.fail()) {
    // RDKit❗✔️:       lineStream.get(spacer);
    // RDKit❗✔️:       if (lineStream.gcount() == 0) {
    // RDKit❗✔️:     continue;
    // RDKit❗✔️:       } else if (spacer != ' ') {
    // RDKit❗✔️:     if (strictParsing) {
    // RDKit❗✔️:       throw FileParseException(errout.str());
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:     sgroup.setIsValid(false);
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:       ParseV3000ParseLabel(label, lineStream, dataFields, line, sgroup,
    // RDKit❗✔️:                            nSgroups, mol, strictParsing);
    // RDKit❗✔️: }
    // RDKit❗✔️:   } catch (const std::exception &e) {
    // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, e.what());
    // RDKit❗✔️:     sgroup.setIsValid(false);
    // RDKit❗✔️:   return;
    // RDKit❗✔️: }
    let mut cursor = SGroupLineCursor::new(text);
    while cursor.position < cursor.text.len() && !cursor.failed {
        let Some(spacer) = cursor.get() else {
            continue;
        };
        if spacer != b' ' {
            let error = SdfReadError::Parse(format!(
                "Found character when expecting a separator (space) on line {line_number}"
            ));
            return if strict_parsing {
                Err(error)
            } else {
                Ok(false)
            };
        }
        let Some(label) = cursor.read_label() else {
            continue;
        };
        if label.is_empty() {
            continue;
        }
        if defaults_only && seen.contains(label) {
            if let Err(error) = skip_overridden_default_value(&mut cursor, label, line_number) {
                return if strict_parsing {
                    Err(error)
                } else {
                    Ok(false)
                };
            }
        } else {
            if let Err(error) = parse_label(
                label,
                &mut cursor,
                line_number,
                group,
                parents,
                atoms,
                bonds,
                bond_endpoints,
                strict_parsing,
            ) {
                if strict_parsing {
                    return Err(error);
                }
                if matches!(&error, SdfReadError::Parse(message) if message == "invalid count value")
                {
                    seen.insert(label.to_owned());
                    continue;
                }
                return Ok(false);
            }
        }
        // The source's `parsedLabels` records only the row's explicit labels;
        // defaults do not add to it, so repeated non-overridden defaults are
        // all parsed in encounter order and later assignments may overwrite.
        if !defaults_only {
            seen.insert(label.to_owned());
        }
    }
    Ok(true)
    // END RDKIT CPP FUNCTION
}

pub(super) fn parse_v3000_sgroup_block(
    lines: &crate::sdf::MolBlockLines<'_>,
    cursor: &mut usize,
    expected_count: usize,
    atoms: &BTreeMap<u32, AtomId>,
    bonds: &BTreeMap<u32, BondId>,
    bond_endpoints: &dyn Fn(BondId) -> Option<(AtomId, AtomId)>,
    strict_parsing: bool,
) -> Result<Vec<SubstanceGroup>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock
    // RDKit❗✔️:   // SGroups may be written in unsorted ID order, according to spec, so we will
    // RDKit❗✔️:   // temporarily store them in a map before adding them to the mol
    // RDKit❗✔️:   IDX_TO_SGROUP_MAP sGroupMap;
    // RDKit❗✔️:   auto tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit❗✔️:   if (tempStr.substr(0, 7) == "DEFAULT" && tempStr.length() > 8) {
    // RDKit❗✔️:     defaultString = tempStr.substr(7);
    // RDKit❗✔️:     tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit❗✔️: }
    // RDKit❗✔️:   for (unsigned int si = 0; si < nSgroups; ++si) {
    // RDKit❗✔️:     lineStream >> sequenceId;
    // RDKit❗✔️:     lineStream >> type;
    // RDKit❗✔️:     lineStream >> externalId;
    // RDKit❗✔️:     if (strictParsing && !SubstanceGroupChecks::isValidType(type)) {
    // RDKit❗✔️:       throw MolFileUnhandledFeatureException(errout.str());
    // RDKit❗✔️:   }
    // RDKit❗✔️:     SubstanceGroup sgroup(mol, type);
    // RDKit❗✔️:       sgroup.setProp<unsigned int>("index", sequenceId);
    // RDKit❗✔️:     if (externalId > 0) {
    // RDKit❗✔️:       sgroup.setProp<unsigned int>("ID", externalId);
    // RDKit❗✔️:   }
    // RDKit❗✔️:       ParseV3000ParseLabel(label, lineStream, dataFields, line, sgroup,
    // RDKit❗✔️:                            nSgroups, mol, strictParsing);
    // RDKit❗✔️:     // Process defaults
    // RDKit❗✔️:       if (std::find(parsedLabels.begin(), parsedLabels.end(), label) ==
    // RDKit❗✔️:           parsedLabels.end()) {
    // RDKit❗✔️:         ParseV3000ParseLabel(label, lineStream, dataFields, defaultLineNum,
    // RDKit❗✔️:                              sgroup, nSgroups, mol, strictParsing);
    // RDKit❗✔️:   }
    // RDKit❗✔️:       sGroupMap.emplace(sequenceId, sgroup);
    // RDKit❗✔️:     tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit❗✔️: }
    // RDKit❗✔️:   if (sGroupMap.size() != nSgroups) {
    // RDKit❗✔️:   std::ostringstream errout;
    // RDKit❗✔️:     errout << "Found " << sGroupMap.size() << " SGroups when " << nSgroups
    // RDKit❗✔️:            << " were expected." << std::endl;
    // RDKit❗✔️:   if (strictParsing) {
    // RDKit❗✔️:     throw FileParseException(errout.str());
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️:   for (const auto &sg : sGroupMap) {
    // RDKit❗✔️:     if (sg.second.getIsValid()) {
    // RDKit❗✔️:       addSubstanceGroup(*mol, sg.second);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << "SGroup " << sg.first
    // RDKit❗✔️:                               << " is invalid and will be ignored" << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let (mut current, mut line_number) = get_v3000_line(lines, cursor)?;
    let mut defaults = String::new();
    let mut default_line = line_number;
    if current.starts_with("DEFAULT") && current.len() > 8 {
        defaults = current[7..].trim_end().to_owned();
        default_line = line_number;
        (current, line_number) = get_v3000_line(lines, cursor)?;
    }

    let mut groups = BTreeMap::<u32, SubstanceGroup>::new();
    let mut invalid_sequences = BTreeSet::new();
    let mut parents = BTreeMap::<u32, u32>::new();
    for _ in 0..expected_count {
        if current.to_ascii_uppercase().starts_with("END SGROUP") {
            break;
        }
        let (sequence_text, kind_text, external_text, labels) =
            split_sgroup_line(current.trim_end(), line_number)?;
        let sequence = parse_sgroup_header_unsigned(sequence_text, "sequence ID", line_number)?;
        if strict_parsing && !is_valid_rdkit_sgroup_type(kind_text) {
            return Err(SdfReadError::Parse(format!(
                "Unsupported SGroup type '{kind_text}' on line {line_number}"
            )));
        }
        let external_id = parse_sgroup_header_unsigned(external_text, "external ID", line_number)?;
        let mut group = SubstanceGroup::new(
            SubstanceGroupId::new(groups.len()),
            sgroup_kind_from_rdkit_type(kind_text),
        );
        group.set_rdkit_sequence_id(sequence);
        group.set_prop("TYPE", kind_text)?;
        group.set_prop("index", PropertyValue::UInt(sequence))?;
        if external_id != 0 {
            group.set_prop("ID", PropertyValue::UInt(external_id))?;
            group.set_external_id(external_id);
        }
        let mut seen = BTreeSet::new();
        let mut candidate_parents = BTreeMap::new();
        let valid = apply_labels(
            &labels,
            line_number,
            &mut group,
            &mut candidate_parents,
            atoms,
            bonds,
            bond_endpoints,
            &mut seen,
            false,
            strict_parsing,
        )?;
        let valid = if valid {
            apply_labels(
                &defaults,
                default_line,
                &mut group,
                &mut candidate_parents,
                atoms,
                bonds,
                bond_endpoints,
                &mut seen,
                true,
                strict_parsing,
            )?
        } else {
            false
        };
        group.set_prop(
            "DATAFIELDS",
            PropertyValue::StringVector(
                group
                    .data()
                    .map_or_else(Vec::new, |data| data.values.clone()),
            ),
        )?;
        if let std::collections::btree_map::Entry::Vacant(entry) = groups.entry(sequence) {
            entry.insert(group);
            parents.extend(candidate_parents);
            if !valid {
                invalid_sequences.insert(sequence);
            }
        }
        (current, line_number) = get_v3000_line(lines, cursor)?;
    }
    if !current.to_ascii_uppercase().starts_with("END SGROUP") {
        if strict_parsing {
            return Err(SdfReadError::Parse(format!(
                "END SGROUP line not found on line {line_number}"
            )));
        }
        *cursor = cursor.saturating_sub(1);
    }
    if groups.len() != expected_count && strict_parsing {
        return Err(SdfReadError::Parse(format!(
            "Found {} SGroups when {expected_count} were expected.",
            groups.len()
        )));
    }

    groups.retain(|sequence, _| !invalid_sequences.contains(sequence));

    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock (installation)
    // RDKit❗✔️:   // SGroups successfully parsed, now add them to the molecule
    // RDKit❗✔️:   for (const auto &sg : sGroupMap) {
    // RDKit❗✔️:     if (sg.second.getIsValid()) {
    // RDKit❗✔️:       addSubstanceGroup(*mol, sg.second);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << "SGroup " << sg.first
    // RDKit❗✔️:                               << " is invalid and will be ignored" << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit stores PARENT as an unchecked unsigned property. The canonical
    // COSMolKit model instead requires a typed in-range SubstanceGroupId. In
    // strict mode an unresolved source sequence is therefore a structured
    // model-boundary error. In non-strict mode the complete child is removed,
    // recursively, so removing an invalid parent can never leave a dangling
    // or partially retained hierarchy. Self-parent and multi-row cycles remain
    // representable because all referenced source rows survive and neither the
    // pinned source nor canonical local validation forbids those relations.
    let children_by_parent = parents.iter().fold(
        BTreeMap::<u32, Vec<u32>>::new(),
        |mut children, (sequence, parent_sequence)| {
            children
                .entry(*parent_sequence)
                .or_default()
                .push(*sequence);
            children
        },
    );
    let mut unresolved = parents
        .iter()
        .filter_map(|(sequence, parent_sequence)| {
            (groups.contains_key(sequence) && !groups.contains_key(parent_sequence))
                .then_some((*sequence, *parent_sequence))
        })
        .collect::<Vec<_>>();
    if strict_parsing && let Some((sequence, parent_sequence)) = unresolved.first().copied() {
        return Err(SdfReadError::Parse(format!(
            "SGroup {sequence} references missing parent SGroup {parent_sequence}"
        )));
    }
    while let Some((sequence, _)) = unresolved.pop() {
        if groups.remove(&sequence).is_none() {
            continue;
        }
        if let Some(children) = children_by_parent.get(&sequence) {
            unresolved.extend(
                children
                    .iter()
                    .filter_map(|child| groups.contains_key(child).then_some((*child, sequence))),
            );
        }
    }

    let ids = groups
        .keys()
        .enumerate()
        .map(|(index, sequence)| (*sequence, SubstanceGroupId::new(index)))
        .collect::<BTreeMap<_, _>>();
    for (sequence, id) in &ids {
        groups
            .get_mut(sequence)
            .expect("key came from group map")
            .set_id(*id);
    }
    for (sequence, parent_sequence) in parents {
        let Some(group) = groups.get_mut(&sequence) else {
            continue;
        };
        let parent = ids
            .get(&parent_sequence)
            .expect("unresolved parent groups were removed before compact ID assignment");
        group.set_parent(*parent);
    }
    // Behavioral review: sorted surviving source rows receive compact IDs once
    // and every retained PARENT relation is installed against that same map;
    // strict failure returns no partial record and non-strict cascading removal
    // retains no child whose parent was discarded. This is an explicit typed
    // model boundary beyond RDKit's unchecked property storage.
    // Complexity review: the reverse parent index and removal work are linear
    // apart from ordered-map/set lookups, O(S log S), with each row removed at
    // most once. This preserves the source installation scale without repeated
    // hierarchy scans or cloning group state.
    // END RDKIT CPP FUNCTION
    Ok(groups.into_values().collect())
    // END RDKIT CPP FUNCTION
}

fn parse_stereo_collection_line(
    line: &str,
    line_number: usize,
    atom_count: usize,
) -> Result<Option<StereoGroup>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION parseEnhancedStereo (line recognition)
    // RDKit✔️✔️:   const regex stereo_label(
    // RDKit✔️✔️:       R"regex(MDLV30/STE(...)([0-9]*) +ATOMS=\(([0-9]+) +(.*)\) *)regex");
    // Recognition and the matched tag/index branches are reproduced below;
    // the caller owns duplicate-ABS strictness and collection installation.
    // Marker review: behavior verified by the closed collection regressions
    // (skip/recognition, id rules, row positions, strictness); the manual
    // single-pass scan replaces the source's per-call `std::regex`
    // construction and match without changing any accept/reject outcome.
    // `regex_match` requires the whole line to match: the tag is exactly three
    // characters, the optional group id is digits only, one or more spaces
    // separate it from `ATOMS=(`, the count is `[0-9]+` followed by spaces,
    // and only spaces may follow the closing parenthesis. The block caller
    // reproduces the source's asymmetric normalization: its first logical
    // payload is uppercased, while subsequent payloads reach this matcher
    // with their original case.
    // Non-matching lines are unrecognized collection types and are skipped,
    // not parsed and not errors.
    let Some(rest) = line.strip_prefix("MDLV30/STE") else {
        return Ok(None);
    };
    if rest.len() < 3 {
        return Ok(None);
    }
    // The source regex consumes exactly three bytes before ASCII digits/spaces.
    // A split inside UTF-8 cannot satisfy that suffix; it is an unrecognized
    // collection line, not a Rust string slicing panic.
    let Some(tag) = rest.get(..3) else {
        return Ok(None);
    };
    let after_tag = &rest[3..];
    let id_digits = after_tag.bytes().take_while(u8::is_ascii_digit).count();
    let after_digits = &after_tag[id_digits..];
    let after_id_separator = after_digits.trim_start_matches(' ');
    if after_id_separator.len() == after_digits.len() {
        return Ok(None);
    }
    let Some(after_atoms) = after_id_separator.strip_prefix("ATOMS=(") else {
        return Ok(None);
    };
    let count_digits = after_atoms.bytes().take_while(u8::is_ascii_digit).count();
    if count_digits == 0 {
        return Ok(None);
    }
    let count_text = &after_atoms[..count_digits];
    let after_count = &after_atoms[count_digits..];
    let atoms_region = after_count.trim_start_matches(' ');
    if atoms_region.len() == after_count.len() {
        return Ok(None);
    }
    // `(.*)\)` is greedy: the group ends at the last ')' that is followed
    // only by the regex's trailing spaces.
    let Some(close) = atoms_region.rfind(')') else {
        return Ok(None);
    };
    if !atoms_region[close + 1..].bytes().all(|byte| byte == b' ') {
        return Ok(None);
    }
    let atoms_text = &atoms_region[..close];
    // BEGIN RDKIT CPP FUNCTION parseEnhancedStereo (tag dispatch)
    // RDKit✔️✔️:       if (match[1] == "ABS") {
    // RDKit✔️✔️:         grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
    // RDKit✔️✔️:       } else if (match[1] == "REL") {
    // RDKit✔️✔️:         grouptype = RDKit::StereoGroupType::STEREO_OR;
    // RDKit✔️✔️:         groupid = FileParserUtils::toUnsigned(match[2], true);
    // RDKit✔️✔️:       } else if (match[1] == "RAC") {
    // RDKit✔️✔️:         grouptype = RDKit::StereoGroupType::STEREO_AND;
    // RDKit✔️✔️:         groupid = FileParserUtils::toUnsigned(match[2], true);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         errout << "Unrecognized stereogroup type : '" << tempStr << "' on line"
    // RDKit✔️✔️:                << line;
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       }
    let kind = match tag {
        "ABS" => StereoGroupKind::Absolute,
        "REL" => StereoGroupKind::Or,
        "RAC" => StereoGroupKind::And,
        _ => {
            return Err(SdfReadError::Parse(format!(
                "Unrecognized stereogroup type : '{line}' on line{line_number}"
            )));
        }
    };
    // The source constructs every group with its `groupid`: ABS keeps the
    // initialized zero, while REL/RAC run their digits through
    // `FileParserUtils::toUnsigned` — which resolves empty or overflowing
    // digit strings to zero because the `std::from_chars` error code is
    // ignored (the regex already guarantees digits, so no cast throw exists).
    let group_id = if tag == "ABS" {
        0
    } else {
        parse_rdkit_unsigned(&after_tag[..id_digits]).map_err(|()| {
            SdfReadError::Parse(format!(
                "Cannot convert stereo group id on line {line_number}"
            ))
        })?
    };
    // RDKit✔️✔️:       const unsigned int count = FileParserUtils::toUnsigned(match[3], true);
    // RDKit✔️✔️:       std::vector<Atom *> atoms;
    // RDKit✔️✔️:       std::stringstream ss(match[4]);
    // RDKit✔️✔️:       unsigned int index;
    // RDKit✔️✔️:       for (size_t i = 0; i < count; ++i) {
    // RDKit✔️✔️:         ss >> index;
    // RDKit✔️✔️:         // atoms are 1 indexed in molfiles
    // RDKit✔️✔️:         atoms.push_back(mol->getAtomWithIdx(index - 1));
    // RDKit✔️✔️:       }
    // COLLECTION atom values are 1-based main CTAB atom row positions
    // (`getAtomWithIdx(index - 1)`), not V3000 atom bookmarks: a nonsequential
    // bookmark table cannot change collection membership. `ss >> index`
    // is decimal formatted extraction, not whitespace tokenization: it consumes
    // a numeric prefix, so `1+2` supplies two values and `1x` supplies one.
    // Its sentry skips C-locale whitespace; EOF leaves the previous index
    // untouched. No digits produces zero, unsigned negation wraps, and overflow
    // produces UINT_MAX/failbit. Zero/out-of-range indices fail immediately.
    // EOF before the first index would use uninitialized C++ state: reject it
    // structurally rather than fabricate an atom. Later EOF repeats the last
    // initialized index, including when the declared count exceeds the text.
    let count = parse_rdkit_unsigned(count_text).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{count_text}' to unsigned int on line {line_number}"
        ))
    })? as usize;
    let bytes = atoms_text.as_bytes();
    let mut position_cursor = 0;
    let mut previous_position = None;
    let mut extraction_failed = false;
    let group_atoms = (0..count)
        .map(|_| {
            let out_of_range = |value: String| {
                SdfReadError::Parse(format!(
                    "Stereo group atom index {value} out of range on line {line_number}"
                ))
            };
            while bytes
                .get(position_cursor)
                .is_some_and(|byte| matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c))
            {
                position_cursor += 1;
            }
            let position = if extraction_failed || position_cursor == bytes.len() {
                previous_position.ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "Stereo group has no initialized atom index on line {line_number}"
                    ))
                })?
            } else {
                let negative = bytes[position_cursor] == b'-';
                if matches!(bytes[position_cursor], b'+' | b'-') {
                    position_cursor += 1;
                }
                let digits_start = position_cursor;
                let mut value = Some(0_u32);
                while let Some(byte) = bytes.get(position_cursor).filter(|b| b.is_ascii_digit()) {
                    value =
                        value.and_then(|v| v.checked_mul(10)?.checked_add(u32::from(byte - b'0')));
                    position_cursor += 1;
                }
                extraction_failed = value.is_none() || position_cursor == digits_start;
                let value = match value {
                    None => u32::MAX,
                    Some(value) if negative => value.wrapping_neg(),
                    Some(value) => value,
                };
                previous_position = Some(value);
                value
            };
            let row = usize::try_from(
                position
                    .checked_sub(1)
                    .ok_or_else(|| out_of_range(position.to_string()))?,
            )
            .map_err(|_| out_of_range(position.to_string()))?;
            if row >= atom_count {
                return Err(out_of_range(position.to_string()));
            }
            Ok(AtomId::new(row))
        })
        .collect::<Result<Vec<_>, _>>()?;
    // RDKit✔️✔️:       groups.emplace_back(grouptype, std::move(atoms), std::move(newBonds),
    // RDKit✔️✔️:                           groupid);
    // This constructor argument is the source/read ID; the independent
    // write ID remains zero until an explicit source forwarding operation.
    Ok(Some(
        StereoGroup::new(kind, group_atoms, Vec::new()).with_id(group_id),
    ))
}

pub(super) fn parse_v3000_collection_block(
    lines: &crate::sdf::MolBlockLines<'_>,
    cursor: &mut usize,
    atom_count: usize,
    strict_parsing: bool,
) -> Result<Vec<StereoGroup>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION parseEnhancedStereo
    // RDKit❗✔️: std::string parseEnhancedStereo(std::istream *inStream, unsigned int &line,
    // RDKit❗✔️:                                 RWMol *mol, bool strictParsing) {
    // RDKit❗✔️:   // Lines like (absolute, relative, racemic):
    // RDKit❗✔️:   // M  V30 MDLV30/STEABS ATOMS=(2 2 3)
    // RDKit❗✔️:   // M  V30 MDLV30/STEREL1 ATOMS=(1 12)
    // RDKit❗✔️:   // M  V30 MDLV30/STERAC1 ATOMS=(1 12)
    // RDKit❗✔️:   const regex stereo_label(
    // RDKit❗✔️:       R"regex(MDLV30/STE(...)([0-9]*) +ATOMS=\(([0-9]+) +(.*)\) *)regex");
    // RDKit❗✔️:
    // RDKit❗✔️:   smatch match;
    // RDKit❗✔️:   std::vector<StereoGroup> groups;
    // RDKit❗✔️:
    // RDKit❗✔️:   // Read the collection until the end
    // RDKit❗✔️:   auto tempStr = getV3000Line(inStream, line);
    // RDKit❗✔️:   boost::to_upper(tempStr);
    // RDKit❗✔️:   unsigned abs_group_seen = 0;
    // RDKit❗✔️:   while (!startsWith(tempStr, "END", 3)) {
    // RDKit❗✔️:     // If this line in the collection is part of a stereo group
    // RDKit❗✔️:     if (regex_match(tempStr, match, stereo_label)) {
    // RDKit❗✔️:       StereoGroupType grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
    // RDKit❗✔️:       unsigned groupid = 0;
    // RDKit❗✔️:
    // RDKit❗✔️:       if (match[1] == "ABS") {
    // RDKit❗✔️:         grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
    // RDKit❗✔️:         // Warn only one per mol about multiple ABS groups
    // RDKit❗✔️:         if (abs_group_seen == 1) {
    // RDKit❗✔️:           std::ostringstream errout;
    // RDKit❗✔️:           errout << "Seen a second ABS stereo group on line " << line
    // RDKit❗✔️:                  << std::endl;
    // RDKit❗✔️:           if (strictParsing) {
    // RDKit❗✔️:             throw FileParseException(errout.str());
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:         ++abs_group_seen;
    // RDKit❗✔️:       } else if (match[1] == "REL") {
    // RDKit❗✔️:         grouptype = RDKit::StereoGroupType::STEREO_OR;
    // RDKit❗✔️:         groupid = FileParserUtils::toUnsigned(match[2], true);
    // RDKit❗✔️:       } else if (match[1] == "RAC") {
    // RDKit❗✔️:         grouptype = RDKit::StereoGroupType::STEREO_AND;
    // RDKit❗✔️:         groupid = FileParserUtils::toUnsigned(match[2], true);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         std::ostringstream errout;
    // RDKit❗✔️:         errout << "Unrecognized stereogroup type : '" << tempStr << "' on line"
    // RDKit❗✔️:                << line;
    // RDKit❗✔️:         throw FileParseException(errout.str());
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       const unsigned int count = FileParserUtils::toUnsigned(match[3], true);
    // RDKit❗✔️:       std::vector<Atom *> atoms;
    // RDKit❗✔️:       std::stringstream ss(match[4]);
    // RDKit❗✔️:       unsigned int index;
    // RDKit❗✔️:       for (size_t i = 0; i < count; ++i) {
    // RDKit❗✔️:         ss >> index;
    // RDKit❗✔️:         // atoms are 1 indexed in molfiles
    // RDKit❗✔️:         atoms.push_back(mol->getAtomWithIdx(index - 1));
    // RDKit❗✔️:       }
    // RDKit❗✔️:       std::vector<Bond *> newBonds;
    // RDKit❗✔️:       groups.emplace_back(grouptype, std::move(atoms), std::move(newBonds),
    // RDKit❗✔️:                           groupid);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // skip collection types we don't know how to read. Only one documented
    // RDKit❗✔️:       // is MDLV30/HILITE
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << "Skipping unrecognized collection type at "
    // RDKit❗✔️:                                  "line "
    // RDKit❗✔️:                               << line << ": " << tempStr << std::endl;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     tempStr = getV3000Line(inStream, line);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (!groups.empty()) {
    // RDKit❗✔️:     mol->setStereoGroups(std::move(groups));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   tempStr = getV3000Line(inStream, line);
    // RDKit❗✔️:   return tempStr;
    // RDKit❗✔️: }
    // The recognition/index extraction is delegated to the private helper;
    // installation is owned by read_v3000_record_detached. The full source
    // remains here; this correction does not certify every parser branch.
    // Behavior review: get_v3000_line preserves case. Exactly the first
    // assembled payload is uppercased before the loop; later payloads and the
    // END-prefix check retain their source case, matching the placement of
    // boost::to_upper above. A new call resets this first-line rule. The
    // remaining partial marker is solely the documented first-extraction
    // undefined-state boundary in parse_stereo_collection_line.
    // Complexity review: each logical line is read and scanned once. Only the
    // first payload allocates an uppercase copy, reducing rather than
    // increasing the source-corresponding linear normalization work; group
    // accumulation remains linear in total recognized membership.
    let mut groups = Vec::new();
    let mut absolute_count = 0_usize;
    let (first_line, mut line_number) = get_v3000_line(lines, cursor)?;
    let mut line = first_line.to_ascii_uppercase();
    while !line.starts_with("END") {
        if let Some(group) = parse_stereo_collection_line(&line, line_number, atom_count)? {
            if group.kind() == StereoGroupKind::Absolute {
                absolute_count += 1;
                if absolute_count > 1 && strict_parsing {
                    return Err(SdfReadError::Parse(format!(
                        "Seen a second ABS stereo group on line {line_number}\n"
                    )));
                }
            }
            groups.push(group);
        }
        (line, line_number) = get_v3000_line(lines, cursor)?;
    }
    Ok(groups)
    // END RDKIT CPP FUNCTION
}

#[derive(Debug)]
pub(super) struct V2000SgroupState {
    strict_parsing: bool,
    groups: BTreeMap<u32, SubstanceGroup>,
    invalid_sequences: BTreeSet<u32>,
    parent_by_sequence: BTreeMap<u32, u32>,
    scd_counter: u32,
    last_data_group: u32,
    current_data_field: String,
    pending_attach_points: Vec<(u32, AtomId, Option<String>)>,
}

impl V2000SgroupState {
    pub(super) fn new(strict_parsing: bool) -> Self {
        Self {
            strict_parsing,
            groups: BTreeMap::new(),
            invalid_sequences: BTreeSet::new(),
            parent_by_sequence: BTreeMap::new(),
            scd_counter: 0,
            last_data_group: 0,
            current_data_field: String::new(),
            pending_attach_points: Vec::new(),
        }
    }
}

impl Default for V2000SgroupState {
    fn default() -> Self {
        Self::new(true)
    }
}

fn parse_v2000_int_field(
    line: &str,
    line_number: usize,
    position: &mut usize,
    is_counter: bool,
) -> Result<u32, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseSGroupIntField
    // RDKit❗✔️:   ++pos;  // Account for separation space
    // RDKit❗✔️:   size_t len = 3 - isFieldCounter;  // field counters are smaller
    // RDKit❗✔️:     fieldValue = FileParserUtils::toInt(text.substr(pos, len));
    // RDKit❗✔️:   pos += len;
    // RDKit❗✔️:   return fieldValue;
    *position += 1;
    let length = 3 - usize::from(is_counter);
    if *position >= line.len() {
        return Err(SdfReadError::Parse(format!(
            "SGroup line too short: '{line}' on line {line_number}"
        )));
    }
    let text = rdkit_substr(line, *position, length);
    let value = parse_rdkit_int(text).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{text}' to int on line {line_number}"
        ))
    })? as u32;
    *position += length;
    Ok(value)
    // END RDKIT CPP FUNCTION
}

fn parse_v2000_double_field(
    line: &str,
    line_number: usize,
    position: &mut usize,
) -> Result<f64, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseSGroupDoubleField
    // RDKit❗✔️:   size_t len = 10;
    // RDKit❗✔️:     fieldValue = FileParserUtils::toDouble(text.substr(pos, len));
    // RDKit❗✔️:   pos += len;
    // RDKit❗✔️:   return fieldValue;
    let length = 10;
    if *position >= line.len() {
        return Err(SdfReadError::Parse(format!(
            "SGroup line too short: '{line}' on line {line_number}"
        )));
    }
    let text = rdkit_substr(line, *position, length);
    let value = parse_rdkit_double(text).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{text}' to double on line {line_number}"
        ))
    })?;
    *position += length;
    Ok(value)
    // END RDKIT CPP FUNCTION
}

impl V2000SgroupState {
    fn group_mut(&mut self, sequence: u32) -> Option<&mut SubstanceGroup> {
        // BEGIN RDKIT CPP FUNCTION FindSgIdx
        // RDKit❗✔️:   auto sgIt = sGroupMap.find(sgIdx);
        // RDKit❗✔️:   if (sgIt == sGroupMap.end()) {
        // RDKit❗✔️:   return nullptr;
        // RDKit❗✔️: }
        // RDKit❗✔️:   return &sgIt->second;
        self.groups.get_mut(&sequence)
        // END RDKIT CPP FUNCTION
    }

    fn group_mut_if_present(&mut self, sequence: u32) -> Option<&mut SubstanceGroup> {
        self.group_mut(sequence)
    }

    fn recoverable<T>(&self, result: Result<T, SdfReadError>) -> Result<Option<T>, SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupIntField (policy overload)
        // RDKit❗✔️:     try {
        // RDKit❗✔️:     res = ParseSGroupIntField(text, line, pos, isFieldCounter);
        // RDKit❗✔️:   } catch (const std::exception &e) {
        // RDKit❗✔️:   if (strictParsing) {
        // RDKit❗✔️:     throw;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     ok = false;
        // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << e.what() << std::endl;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        match result {
            Ok(value) => Ok(Some(value)),
            Err(error) if self.strict_parsing => Err(error),
            Err(_) => Ok(None),
        }
        // END RDKIT CPP FUNCTION
    }

    fn invalidate_or_throw(
        &mut self,
        sequence: u32,
        error: SdfReadError,
    ) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION SGroupWarnOrThrow
        // RDKit❗✔️:     if (strictParsing) {
        // RDKit❗✔️:     throw Exc(msg);
        // RDKit❗✔️: } else {
        // RDKit❗✔️:   BOOST_LOG(rdWarningLog) << msg << std::endl;
        // RDKit❗✔️: }
        if self.strict_parsing {
            Err(error)
        } else {
            self.invalid_sequences.insert(sequence);
            Ok(())
        }
        // END RDKIT CPP FUNCTION
    }

    pub(super) fn parse_line(
        &mut self,
        line: &str,
        line_number: usize,
        atom_count: usize,
        bond_endpoints: &[(AtomId, AtomId)],
    ) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseMolBlockProperties (SGroup dispatch)
        // RDKit❗✔️:     } else if (lineBeg == "M  STY") {
        // RDKit❗✔️:       ParseSGroupV2000STYLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SST") {
        // RDKit❗✔️:       ParseSGroupV2000SSTLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SLB") {
        // RDKit❗✔️:       ParseSGroupV2000SLBLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SCN") {
        // RDKit❗✔️:       ParseSGroupV2000SCNLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SDS") {
        // RDKit❗✔️:       ParseSGroupV2000SDSLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SAL" || lineBeg == "M  SBL" ||
        // RDKit❗✔️:                lineBeg == "M  SPA") {
        // RDKit❗✔️:       ParseSGroupV2000VectorDataLine(sGroupMap, mol, tempStr, line,
        // RDKit❗✔️:                                      strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SMT") {
        // RDKit❗✔️:       ParseSGroupV2000SMTLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SDI") {
        // RDKit❗✔️:       ParseSGroupV2000SDILine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SBV") {
        // RDKit❗✔️:       ParseSGroupV2000SBVLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SDT") {
        // RDKit❗✔️:       ParseSGroupV2000SDTLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SDD") {
        // RDKit❗✔️:       ParseSGroupV2000SDDLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SCD" || lineBeg == "M  SED") {
        // RDKit❗✔️:       ParseSGroupV2000SCDSEDLine(sGroupMap, dataFieldsMap, mol, tempStr, line,
        // RDKit❗✔️:                                  strictParsing, SCDcounter, lastDataSGroup,
        // RDKit❗✔️:                                  currentDataField);
        // RDKit❗✔️:     } else if (lineBeg == "M  SPL") {
        // RDKit❗✔️:       ParseSGroupV2000SPLLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SNC") {
        // RDKit❗✔️:       ParseSGroupV2000SNCLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SAP") {
        // RDKit❗✔️:       ParseSGroupV2000SAPLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SCL") {
        // RDKit❗✔️:       ParseSGroupV2000SCLLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️:     } else if (lineBeg == "M  SBT") {
        // RDKit❗✔️:       ParseSGroupV2000SBTLine(sGroupMap, mol, tempStr, line, strictParsing);
        match rdkit_substr(line, 0, 6) {
            "M  STY" => self.parse_sty(line, line_number),
            "M  SST" => self.parse_sst(line, line_number),
            "M  SLB" => self.parse_slb(line, line_number),
            "M  SCN" => self.parse_scn(line, line_number),
            "M  SDS" => self.parse_sds(line, line_number),
            "M  SAL" | "M  SBL" | "M  SPA" => {
                self.parse_vector(line, line_number, atom_count, bond_endpoints.len())
            }
            "M  SMT" => self.parse_smt(line, line_number),
            "M  SDI" => self.parse_sdi(line, line_number),
            "M  SBV" => self.parse_sbv(line, line_number, bond_endpoints),
            "M  SDT" => self.parse_sdt(line, line_number),
            "M  SDD" => self.parse_sdd(line, line_number),
            "M  SCD" | "M  SED" => self.parse_scd_sed(line, line_number),
            "M  SPL" => self.parse_spl(line, line_number),
            "M  SNC" => self.parse_snc(line, line_number),
            "M  SAP" => self.parse_sap(line, line_number, atom_count),
            "M  SCL" => self.parse_scl(line, line_number),
            "M  SBT" => self.parse_sbt(line, line_number),
            prefix => Err(SdfReadError::Parse(format!(
                "Unknown SGroup line '{prefix}' on line {line_number}"
            ))),
        }
        // END RDKIT CPP FUNCTION
    }

    fn parse_sty(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000STYLine
        // RDKit❗✔️:   unsigned int pos = 6;
        // RDKit❗✔️:   unsigned int nent =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️:     unsigned int sequenceId =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     std::string typ = text.substr(pos + 1, 3);
        // RDKit❗✔️:     if (SubstanceGroupChecks::isValidType(typ)) {
        // RDKit❗✔️:       auto sgroup = SubstanceGroup(mol, typ);
        // RDKit❗✔️:       sgroup.setProp<unsigned int>("index", sequenceId);
        // RDKit❗✔️:       sGroupMap.emplace(sequenceId, sgroup);
        // RDKit❗✔️: }
        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup STY line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            let kind_text = rdkit_substr(line, position + 1, 3);
            if !is_valid_rdkit_sgroup_type(kind_text) {
                if self.strict_parsing {
                    return Err(SdfReadError::Parse(format!(
                        "S group {kind_text} on line {line_number}"
                    )));
                }
                position += 4;
                continue;
            }
            let mut group = SubstanceGroup::new(
                SubstanceGroupId::new(self.groups.len()),
                sgroup_kind_from_rdkit_type(kind_text),
            );
            group.set_rdkit_sequence_id(sequence);
            group.set_prop("TYPE", kind_text)?;
            group.set_prop("index", PropertyValue::UInt(sequence))?;
            self.groups.entry(sequence).or_insert(group);
            position += 4;
        }
        Ok(())
        // END RDKIT CPP FUNCTION
    }

    fn parse_sst(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SSTLine
        // RDKit❗✔️: void ParseSGroupV2000SSTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit❗✔️:                              const std::string &text, unsigned int &line,
        // RDKit❗✔️:                              bool strictParsing) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(text.substr(0, 6) == "M  SST", "bad SST line");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 6;
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int nent =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   for (unsigned int ie = 0; ie < nent; ++ie) {
        // RDKit❗✔️:     if (text.size() < pos + 8) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "SGroup SST line too short: '" << text << "' on line " << line;
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     unsigned int sgIdx =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:     if (!sgroup) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     };
        // RDKit❗✔️:
        // RDKit❗✔️:     std::string subType = text.substr(++pos, 3);
        // RDKit❗✔️:
        // RDKit❗✔️:     if (!SubstanceGroupChecks::isValidSubType(subType)) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "Unsupported SGroup subtype '" << subType << "' on line "
        // RDKit❗✔️:              << line;
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     sgroup->setProp("SUBTYPE", subType);
        // RDKit❗✔️:     pos += 3;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.

        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SST line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            let subtype = rdkit_substr(line, position + 1, 3);
            if !matches!(subtype, "ALT" | "RAN" | "BLO") {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "Unsupported SGroup subtype '{subtype}' on line {line_number}"
                    )),
                );
            }
            {
                let group = self
                    .group_mut_if_present(sequence)
                    .expect("presence checked before parsing SGroup subtype");
                group.set_prop("SUBTYPE", subtype)?;
                group.set_subtype(subtype);
            }
            position += 4;
        }
        Ok(())
    }

    fn parse_slb(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SLBLine
        // RDKit❗✔️: void ParseSGroupV2000SLBLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit❗✔️:                              const std::string &text, unsigned int line,
        // RDKit❗✔️:                              bool strictParsing) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(text.substr(0, 6) == "M  SLB", "bad SLB line");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 6;
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int nent =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   for (unsigned int ie = 0; ie < nent; ++ie) {
        // RDKit❗✔️:     if (text.size() < pos + 8) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "SGroup SLB line too short: '" << text << "' on line " << line;
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     unsigned int sgIdx =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:     if (!sgroup) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     unsigned int id = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (id != 0 && !SubstanceGroupChecks::isSubstanceGroupIdFree(*mol, id)) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "SGroup ID '" << id
        // RDKit❗✔️:              << "' is assigned to more than one SGroup, on line " << line;
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     sgroup->setProp<unsigned int>("ID", id);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.
        // This MolBlock-only owner starts with no committed SGroups; source
        // isSubstanceGroupIdFree scans only mol groups, not sGroupMap staging.

        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SLB line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            let Some(id) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            {
                let group = self
                    .group_mut_if_present(sequence)
                    .expect("presence checked before parsing SGroup label ID");
                group.set_prop("ID", PropertyValue::UInt(id))?;
                group.set_external_id(id);
            }
        }
        Ok(())
    }

    fn parse_scn(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // RDKit✔️✔️: void ParseSGroupV2000SCNLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit✔️✔️:                              const std::string &text, unsigned int line,
        // RDKit✔️✔️:                              bool strictParsing) {
        // RDKit✔️✔️:   PRECONDITION(mol, "bad mol");
        // RDKit✔️✔️:   PRECONDITION(text.substr(0, 6) == "M  SCN", "bad SCN line");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   unsigned int pos = 6;
        // RDKit✔️✔️:   bool ok;
        // RDKit✔️✔️:   unsigned int nent =
        // RDKit✔️✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit✔️✔️:   if (!ok) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   for (unsigned int ie = 0; ie < nent; ++ie) {
        // RDKit✔️✔️:     if (text.size() < pos + 7) {
        // RDKit✔️✔️:       std::ostringstream errout;
        // RDKit✔️✔️:       errout << "SGroup SCN line too short: '" << text << "' on line " << line;
        // RDKit✔️✔️:       errout << "\n needed: " << pos + 7 << " found: " << text.size();
        // RDKit✔️✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit✔️✔️:       return;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     unsigned int sgIdx =
        // RDKit✔️✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit✔️✔️:     if (!ok) {
        // RDKit✔️✔️:       return;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit✔️✔️:     if (!sgroup) {
        // RDKit✔️✔️:       return;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     std::string connect = text.substr(++pos, 2);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     if (!SubstanceGroupChecks::isValidConnectType(connect)) {
        // RDKit✔️✔️:       std::ostringstream errout;
        // RDKit✔️✔️:       errout << "Unsupported SGroup connection type '" << connect
        // RDKit✔️✔️:              << "' on line " << line;
        // RDKit✔️✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit✔️✔️:       sgroup->setIsValid(false);
        // RDKit✔️✔️:       return;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     sgroup->setProp("CONNECT", connect);
        // RDKit✔️✔️:     pos += 3;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // Source cursor advances one separator plus the final three bytes
        // for every pair; splitting multiple entries adds no allocations.
        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 7 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SCN line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            position += 1;
            let text = rdkit_substr(line, position, 2);
            let Some(connection) = sgroup_connection_from_rdkit(text) else {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "Unsupported SGroup connection type '{text}' on line {line_number}"
                    )),
                );
            };
            {
                let group = self
                    .group_mut_if_present(sequence)
                    .expect("presence checked before parsing SGroup connection");
                group.set_prop("CONNECT", text)?;
                group.set_connection(connection);
            }
            position += 3;
        }
        Ok(())
    }

    fn parse_sds(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SDSLine
        // RDKit❗✔️: void ParseSGroupV2000SDSLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit❗✔️:                              const std::string &text, unsigned int line,
        // RDKit❗✔️:                              bool strictParsing) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(text.substr(0, 10) == "M  SDS EXP", "bad SDS line");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 10;
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int nent =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   for (unsigned int ie = 0; ie < nent; ++ie) {
        // RDKit❗✔️:     if (text.size() < pos + 4) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "SGroup SDS line too short: '" << text << "' on line " << line;
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     unsigned int sgIdx =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:     if (!sgroup) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     sgroup->setProp("ESTATE", "E");
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.

        if !line.starts_with("M  SDS EXP") {
            return Err(SdfReadError::Parse(format!(
                "bad SDS line on line {line_number}"
            )));
        }
        let mut position = 10;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 4 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SDS line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            let Some(group) = self.group_mut_if_present(sequence) else {
                return Ok(());
            };
            group.set_expansion_state("E");
        }
        Ok(())
    }

    fn parse_vector(
        &mut self,
        line: &str,
        line_number: usize,
        atom_count: usize,
        bond_count: usize,
    ) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000VectorDataLine
        // RDKit❗✔️:   if (typ == "SAL") {
        // RDKit❗✔️:     sGroupAddIndexedElement = &SubstanceGroup::addAtomWithBookmark;
        // RDKit❗✔️:   } else if (typ == "SBL") {
        // RDKit❗✔️:     sGroupAddIndexedElement = &SubstanceGroup::addBondWithBookmark;
        // RDKit❗✔️:   } else if (typ == "SPA") {
        // RDKit❗✔️:     sGroupAddIndexedElement = &SubstanceGroup::addParentAtomWithBookmark;
        // RDKit❗✔️: }
        // RDKit❗✔️:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   unsigned int nent =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️:       if (!ok) {
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️:     if (text.size() < pos + 4) {
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️:     unsigned int nbr = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:       if (!ok) {
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️:       (sgroup->*sGroupAddIndexedElement)(nbr);
        let kind = rdkit_substr(line, 3, 3);
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            self.invalid_sequences.insert(sequence);
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 4 {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "SGroup line too short: '{line}' on line {line_number}"
                    )),
                );
            }
            let Some(bookmark) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            let Some(index) = bookmark.checked_sub(1).map(|index| index as usize) else {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "SGroup {kind} index 0 out of range on line {line_number}"
                    )),
                );
            };
            if (kind == "SBL" && index >= bond_count) || (kind != "SBL" && index >= atom_count) {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "SGroup {kind} index {bookmark} out of range on line {line_number}"
                    )),
                );
            }
            let group = self
                .group_mut_if_present(sequence)
                .expect("presence checked before parsing vector entries");
            match kind {
                "SAL" => group.push_atom(AtomId::new(index)),
                "SBL" => group.push_bond(BondId::new(index)),
                "SPA" => {
                    if !group.atoms().contains(&AtomId::new(index)) {
                        return self.invalidate_or_throw(sequence, SdfReadError::Parse(format!(
                            "Atom with bookmark {bookmark} is not a member of current SubstanceGroup on line {line_number}"
                        )));
                    }
                    group.push_parent_atom(AtomId::new(index));
                }
                _ => unreachable!("dispatch limits V2000 vector tags"),
            }
        }
        Ok(())
        // END RDKIT CPP FUNCTION
    }

    fn parse_smt(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SMTLine
        // RDKit❗✔️: void ParseSGroupV2000SMTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit❗✔️:                              const std::string &text, unsigned int &line,
        // RDKit❗✔️:                              bool strictParsing) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(text.substr(0, 6) == "M  SMT", "bad SMT line");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 6;
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:   if (!sgroup) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   ++pos;
        // RDKit❗✔️:
        // RDKit❗✔️:   if (pos >= text.length()) {
        // RDKit❗✔️:     std::ostringstream errout;
        // RDKit❗✔️:     errout << "SGroup line too short: '" << text << "' on line " << line;
        // RDKit❗✔️:     SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   std::string label = text.substr(pos, text.length() - pos);
        // RDKit❗✔️:
        // RDKit❗✔️:   if (sgroup->getProp<std::string>("TYPE") ==
        // RDKit❗✔️:       "MUL") {  // Case of multiple groups
        // RDKit❗✔️:     sgroup->setProp("MULT", label);
        // RDKit❗✔️:
        // RDKit❗✔️:   } else {  // Case of abbreviation groups, but we might not have seen a SCL
        // RDKit❗✔️:             // line yet
        // RDKit❗✔️:     sgroup->setProp("LABEL", label);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.

        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        position += 1;
        if position >= line.len() {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "SGroup line too short: '{line}' on line {line_number}"
                )),
            );
        }
        let label = &line[position..];
        let group = self
            .group_mut_if_present(sequence)
            .expect("presence checked before parsing SGroup label");
        if group.kind() == &SubstanceGroupKind::MultipleGroup {
            group.set_prop("MULT", label)?;
        } else {
            group.set_prop("LABEL", label)?;
            group.set_label(label);
        }
        Ok(())
    }

    fn parse_sdi(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SDILine
        // RDKit❗✔️: void ParseSGroupV2000SDILine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit❗✔️:                              const std::string &text, unsigned int line,
        // RDKit❗✔️:                              bool strictParsing) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(text.substr(0, 6) == "M  SDI", "bad SDI line");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 6;
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:   if (!sgroup) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int nCoords =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (nCoords != 4) {
        // RDKit❗✔️:     std::ostringstream errout;
        // RDKit❗✔️:     errout << "Unexpected number of coordinates for SDI on line " << line;
        // RDKit❗✔️:     SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   SubstanceGroup::Bracket bracket;
        // RDKit❗✔️:   for (unsigned int i = 0; i < 2; ++i) {
        // RDKit❗✔️:     double x = ParseSGroupDoubleField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     double y = ParseSGroupDoubleField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     double z = 0.;
        // RDKit❗✔️:     bracket[i] = RDGeom::Point3D(x, y, z);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   bracket[2] = RDGeom::Point3D(0., 0., 0.);
        // RDKit❗✔️:   try {
        // RDKit❗✔️:     sgroup->addBracket(bracket);
        // RDKit❗✔️:   } catch (const std::exception &e) {
        // RDKit❗✔️:     SGroupWarnOrThrow<>(strictParsing, e.what());
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.

        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            self.invalid_sequences.insert(sequence);
            return Ok(());
        };
        if count != 4 {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "Unexpected number of coordinates for SDI on line {line_number}"
                )),
            );
        }
        let mut coordinate = [0.0; 4];
        for value in &mut coordinate {
            let Some(parsed) =
                self.recoverable(parse_v2000_double_field(line, line_number, &mut position))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            *value = parsed;
        }
        self.group_mut_if_present(sequence)
            .expect("presence checked before parsing SGroup bracket")
            .display_mut()
            .brackets
            .push(SGroupBracket {
                points: [
                    [coordinate[0], coordinate[1], 0.0],
                    [coordinate[2], coordinate[3], 0.0],
                    [0.0, 0.0, 0.0],
                ],
            });
        Ok(())
    }

    fn parse_sbv(
        &mut self,
        line: &str,
        line_number: usize,
        bond_endpoints: &[(AtomId, AtomId)],
    ) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SBVLine
        // RDKit❗✔️:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:       if (!ok) {
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️:   SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:   if (!sgroup) {
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️:   unsigned int bondMark =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:       if (!ok) {
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️:   if (sgroup->getProp<std::string>("TYPE") == "SUP") {
        // RDKit❗✔️:     vector.x = ParseSGroupDoubleField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:     vector.y = ParseSGroupDoubleField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:     try {
        // RDKit❗✔️:     sgroup->addCState(bond->getIdx(), vector);
        // RDKit❗✔️:   } catch (const std::exception &e) {
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, e.what());
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        let Some(group) = self.group_mut_if_present(sequence) else {
            return Ok(());
        };
        let is_superatom = group.kind() == &SubstanceGroupKind::Superatom;
        let Some(bookmark) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            self.invalid_sequences.insert(sequence);
            return Ok(());
        };
        let Some(index) = bookmark.checked_sub(1).map(|index| index as usize) else {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "SGroup bond index 0 out of range on line {line_number}"
                )),
            );
        };
        let Some(&(begin, end)) = bond_endpoints.get(index) else {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "SGroup bond index {bookmark} out of range on line {line_number}"
                )),
            );
        };
        let vector = if is_superatom {
            let Some(x) =
                self.recoverable(parse_v2000_double_field(line, line_number, &mut position))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            let Some(y) =
                self.recoverable(parse_v2000_double_field(line, line_number, &mut position))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            [x, y, 0.0]
        } else {
            [0.0, 0.0, 0.0]
        };
        let group = self
            .group_mut_if_present(sequence)
            .expect("presence checked before parsing CState");
        let bond = BondId::new(index);
        let bond_is_listed = group.bonds().contains(&bond);
        let begin_is_member = group.atoms().contains(&begin);
        let end_is_member = group.atoms().contains(&end);
        if !bond_is_listed || begin_is_member == end_is_member {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "Bond with index {index} is not an XBOND for current SubstanceGroup"
                )),
            );
        }
        group.push_cstate(SGroupCState { bond, vector });
        Ok(())
        // END RDKIT CPP FUNCTION
    }

    fn parse_sdt(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // RDKit✔️✔️: void ParseSGroupV2000SDTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit✔️✔️:                              const std::string &text, unsigned int line,
        // RDKit✔️✔️:                              bool strictParsing) {
        // RDKit✔️✔️:   PRECONDITION(mol, "bad mol");
        // RDKit✔️✔️:   PRECONDITION(text.substr(0, 6) == "M  SDT", "bad SDT line");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   unsigned int pos = 6;
        // RDKit✔️✔️:   bool ok;
        // RDKit✔️✔️:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit✔️✔️:   if (!ok) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit✔️✔️:   if (!sgroup) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::string fieldName;
        // RDKit✔️✔️:   std::string fieldType;
        // RDKit✔️✔️:   std::string fieldInfo;
        // RDKit✔️✔️:   std::string queryType;
        // RDKit✔️✔️:   std::string queryOp;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   try {
        // RDKit✔️✔️:     fieldName = text.substr(++pos, 30);
        // RDKit✔️✔️:     boost::trim_right(fieldName);
        // RDKit✔️✔️:     pos += 30;
        // RDKit✔️✔️:     fieldType = text.substr(pos, 2);
        // RDKit✔️✔️:     boost::trim_right(fieldType);
        // RDKit✔️✔️:     pos += 2;
        // RDKit✔️✔️:     fieldInfo = text.substr(pos, 20);
        // RDKit✔️✔️:     boost::trim_right(fieldInfo);
        // RDKit✔️✔️:     pos += 20;
        // RDKit✔️✔️:     queryType = text.substr(pos, 2);
        // RDKit✔️✔️:     boost::trim_right(queryType);
        // RDKit✔️✔️:     pos += 2;
        // RDKit✔️✔️:     queryOp = text.substr(pos, text.length() - pos);
        // RDKit✔️✔️:     boost::trim_right(queryOp);
        // RDKit✔️✔️:   } catch (const std::out_of_range &) {
        // RDKit✔️✔️:     // all kinds of wild things out there... this insulates us from them without
        // RDKit✔️✔️:     // making the code super complicated
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // only add entries for the remaining properties if they aren't blank
        // RDKit✔️✔️:   if (!fieldName.empty()) {
        // RDKit✔️✔️:     sgroup->setProp("FIELDNAME", fieldName);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!fieldType.empty()) {
        // RDKit✔️✔️:     sgroup->setProp("FIELDTYPE", fieldType);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!fieldInfo.empty()) {
        // RDKit✔️✔️:     sgroup->setProp("FIELDINFO", fieldInfo);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!queryType.empty()) {
        // RDKit✔️✔️:     sgroup->setProp("QUERYTYPE", queryType);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!queryOp.empty()) {
        // RDKit✔️✔️:     sgroup->setProp("QUERYOP", queryOp);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        position += 1;
        let field_name = rdkit_substr(line, position, 30).trim_end().to_owned();
        position += 30;
        let field_type = rdkit_substr(line, position, 2).trim_end().to_owned();
        position += 2;
        let field_info = rdkit_substr(line, position, 20).trim_end().to_owned();
        position += 20;
        let query_type = rdkit_substr(line, position, 2).trim_end().to_owned();
        position += 2;
        let query_op = rdkit_substr(line, position, line.len().saturating_sub(position))
            .trim_end()
            .to_owned();
        let group = self
            .group_mut_if_present(sequence)
            .expect("presence checked before parsing SGroup data header");
        if !field_name.is_empty() {
            group.set_prop("FIELDNAME", field_name.clone())?;
            group.data_mut().field_name = Some(field_name.into());
        }
        if !field_type.is_empty() {
            group.set_prop("FIELDTYPE", field_type.clone())?;
            group.data_mut().field_type = Some(field_type.into());
        }
        if !field_info.is_empty() {
            group.set_prop("FIELDINFO", field_info.clone())?;
            group.data_mut().field_info = Some(field_info.into());
        }
        if !query_type.is_empty() {
            group.set_prop("QUERYTYPE", query_type.clone())?;
            group.data_mut().query_type = Some(query_type.into());
        }
        if !query_op.is_empty() {
            group.set_prop("QUERYOP", query_op.clone())?;
            group.data_mut().query_op = Some(query_op.into());
        }

        Ok(())
    }

    fn parse_sdd(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SDDLine
        // RDKit❗✔️: void ParseSGroupV2000SDDLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit❗✔️:                              const std::string &text, unsigned int line,
        // RDKit❗✔️:                              bool strictParsing) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(text.substr(0, 6) == "M  SDD", "bad SDD line");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 6;
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:   if (!sgroup) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // Store the rest of the line as is.
        // RDKit❗✔️:   ++pos;
        // RDKit❗✔️:   if (pos < text.length()) {
        // RDKit❗✔️:     sgroup->setProp("FIELDDISP", text.substr(pos, text.length() - pos));
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.

        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        position += 1;
        if position < line.len() {
            let group = self
                .group_mut_if_present(sequence)
                .expect("presence checked before parsing SGroup display data");
            group.set_prop("FIELDDISP", &line[position..])?;
            group.data_mut().field_display = Some(line[position..].into());
        }
        Ok(())
    }

    fn parse_scd_sed(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SCDSEDLine
        // RDKit❗✔️: void ParseSGroupV2000SCDSEDLine(IDX_TO_SGROUP_MAP &sGroupMap,
        // RDKit❗✔️:                                 IDX_TO_STR_VECT_MAP &dataFieldsMap, RWMol *mol,
        // RDKit❗✔️:                                 const std::string &text, unsigned int line,
        // RDKit❗✔️:                                 bool strictParsing, unsigned int &counter,
        // RDKit❗✔️:                                 unsigned int &lastDataSGroup,
        // RDKit❗✔️:                                 std::ostringstream &currentDataField) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 3;
        // RDKit❗✔️:   std::string type = text.substr(pos, 3);
        // RDKit❗✔️:   pos += 3;
        // RDKit❗✔️:
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:   if (!sgroup) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   if (lastDataSGroup != 0 && lastDataSGroup != sgIdx) {
        // RDKit❗✔️:     std::ostringstream errout;
        // RDKit❗✔️:     errout << "Found a Data Field not matching the SGroup of the last Data "
        // RDKit❗✔️:               "Field at line "
        // RDKit❗✔️:            << line;
        // RDKit❗✔️:     SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   } else if (lastDataSGroup == 0 && type == "SCD") {
        // RDKit❗✔️:     lastDataSGroup = sgIdx;
        // RDKit❗✔️:   } else if (type == "SED") {
        // RDKit❗✔️:     lastDataSGroup = 0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // have we already seen an SDT line?
        // RDKit❗✔️:   if (!sgroup->hasProp("FIELDNAME")) {
        // RDKit❗✔️:     // one can read the docs and draw the conclusion that this is mandatory,
        // RDKit❗✔️:     // but it's also possible to interpret them the other way, and we know
        // RDKit❗✔️:     // that there are CTABs out there with empty fieldnames in SDT lines,
        // RDKit❗✔️:     // so let's just issue a warning and accept it.
        // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
        // RDKit❗✔️:         << "Found a SCD/SED line with missing/empty SDT specification at line "
        // RDKit❗✔️:         << line << std::endl;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   if (strictParsing) {
        // RDKit❗✔️:     if (type == "SCD" && counter > 2) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "Found too many consecutive SCD lines, (#" << (counter + 1)
        // RDKit❗✔️:              << " at line " << line << ") for SGroup " << sgIdx;
        // RDKit❗✔️:       throw FileParseException(errout.str());
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   if (pos + 1 < text.length()) {
        // RDKit❗✔️:     currentDataField << text.substr(++pos, 69);
        // RDKit❗✔️:
        // RDKit❗✔️:     if (type == "SED") {
        // RDKit❗✔️:       std::string trimmedData = boost::trim_right_copy(currentDataField.str());
        // RDKit❗✔️:       dataFieldsMap[sgIdx].push_back(trimmedData.substr(0, 200));
        // RDKit❗✔️:       currentDataField.str("");
        // RDKit❗✔️:       counter = 0;
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       ++counter;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.

        let kind = rdkit_substr(line, 3, 3);
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        if self.last_data_group != 0 && self.last_data_group != sequence {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "Found a Data Field not matching the SGroup of the last Data Field at line {line_number}"
                )),
            );
        } else if self.last_data_group == 0 && kind == "SCD" {
            self.last_data_group = sequence;
        } else if kind == "SED" {
            self.last_data_group = 0;
        }
        if self.strict_parsing && kind == "SCD" && self.scd_counter > 2 {
            return Err(SdfReadError::Parse(format!(
                "Found too many consecutive SCD lines, (#{} at line {line_number}) for SGroup {sequence}",
                self.scd_counter + 1
            )));
        }
        if position + 1 < line.len() {
            self.current_data_field
                .push_str(rdkit_substr(line, position + 1, 69));
            if kind == "SED" {
                let value = rdkit_substr(self.current_data_field.trim_end(), 0, 200).to_owned();
                self.group_mut_if_present(sequence)
                    .expect("presence checked before parsing SGroup data value")
                    .data_mut()
                    .values
                    .push(value.into());
                self.current_data_field.clear();
                self.scd_counter = 0;
            } else {
                self.scd_counter += 1;
            }
        }
        Ok(())
    }

    fn parse_spl(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SPLLine
        // RDKit❗✔️: void ParseSGroupV2000SPLLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit❗✔️:                              const std::string &text, unsigned int line,
        // RDKit❗✔️:                              bool strictParsing) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(text.substr(0, 6) == "M  SPL", "bad SPL line");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 6;
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int nent =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   for (unsigned int ie = 0; ie < nent; ++ie) {
        // RDKit❗✔️:     if (text.size() < pos + 8) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "SGroup SPL line too short: '" << text << "' on line " << line;
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     unsigned int sgIdx =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:     if (!sgroup) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     unsigned int parentIdx = ParseSGroupIntField(text, line, pos);
        // RDKit❗✔️:
        // RDKit❗✔️:     sgroup->setProp<unsigned int>("PARENT", parentIdx);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.

        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SPL line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            // RDKit intentionally uses the throwing overload for PARENT even
            // when the surrounding parse is non-strict.
            let parent = parse_v2000_int_field(line, line_number, &mut position, false)?;
            self.group_mut_if_present(sequence)
                .expect("presence checked")
                .set_prop("PARENT", PropertyValue::UInt(parent))?;
            self.parent_by_sequence.insert(sequence, parent);
        }
        Ok(())
    }

    fn parse_snc(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SNCLine
        // RDKit❗✔️: void ParseSGroupV2000SNCLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit❗✔️:                              const std::string &text, unsigned int line,
        // RDKit❗✔️:                              bool strictParsing) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(text.substr(0, 6) == "M  SNC", "bad SNC line");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 6;
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int nent =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   for (unsigned int ie = 0; ie < nent; ++ie) {
        // RDKit❗✔️:     if (text.size() < pos + 8) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "SGroup SNC line too short: '" << text << "' on line " << line;
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     unsigned int sgIdx =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:     if (!sgroup) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     unsigned int compno =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (compno > 256u) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "SGroup SNC value over 256: '" << compno << "' on line "
        // RDKit❗✔️:              << line;
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     sgroup->setProp<unsigned int>("COMPNO", compno);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.

        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SNC line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            let Some(component) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            if component > 256 {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "SGroup SNC value over 256: '{component}' on line {line_number}"
                    )),
                );
            }
            {
                let group = self
                    .group_mut_if_present(sequence)
                    .expect("presence checked before parsing SGroup component number");
                group.set_prop("COMPNO", PropertyValue::UInt(component))?;
                group.set_component_number(component);
            }
        }
        Ok(())
    }

    fn parse_sap(
        &mut self,
        line: &str,
        line_number: usize,
        atom_count: usize,
    ) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SAPLine
        // RDKit❗✔️:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:       if (!ok) {
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️:   SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:   if (!sgroup) {
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️:   unsigned int nent =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️:       if (!ok) {
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️:     int lvIdx = -1;
        // RDKit❗✔️:     if (text.size() < pos + 11) {
        // RDKit❗✔️:   if (strictParsing) {
        // RDKit❗✔️:     throw FileParseException(errout.str());
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     if (text.size() < pos + 4) {
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:         lvIdx = mol->getNumAtoms();
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:     std::string id = "  ";
        // RDKit❗✔️:     unsigned int aIdxMark =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (lvIdx == -1) {
        // RDKit❗✔️:       unsigned int lvIdxMark =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:       if (lvIdxMark != 0) {
        // RDKit❗✔️:         lvIdx = mol->getAtomWithBookmark(lvIdxMark)->getIdx();
        // RDKit❗✔️:   }
        // RDKit❗✔️:       if (text.size() >= pos + 3) {
        // RDKit❗✔️:         id = text.substr(pos + 1, 2);
        // RDKit❗✔️:     pos += 3;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️:       sgroup->addAttachPoint(aIdx, lvIdx, id);
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            self.invalid_sequences.insert(sequence);
            return Ok(());
        };
        for _ in 0..count {
            let missing_leaving_atom = line.len() < position + 11;
            if missing_leaving_atom && self.strict_parsing {
                return Err(SdfReadError::Parse(format!(
                    "SGroup SAP line too short: '{line}' on line {line_number}"
                )));
            }
            if missing_leaving_atom && line.len() < position + 4 {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            }
            let Some(atom_bookmark) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            let atom_index = atom_bookmark.checked_sub(1).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup attach atom index 0 out of range on line {line_number}"
                ))
            })? as usize;
            if atom_index >= atom_count {
                return Err(SdfReadError::Parse(format!(
                    "SGroup attach atom index out of range on line {line_number}"
                )));
            }
            let atom = AtomId::new(atom_index);
            if missing_leaving_atom {
                self.pending_attach_points
                    .push((sequence, atom, Some("  ".to_owned())));
                continue;
            }
            let Some(leaving_bookmark) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            if leaving_bookmark as usize > atom_count {
                return Err(SdfReadError::Parse(format!(
                    "SGroup attach atom index out of range on line {line_number}"
                )));
            }
            let label = if line.len() >= position + 3 {
                let label = rdkit_substr(line, position + 1, 2).to_owned();
                position += 3;
                Some(label)
            } else {
                Some("  ".to_owned())
            };
            self.group_mut_if_present(sequence)
                .expect("presence checked before parsing attachment points")
                .push_attach_point(SGroupAttachPoint {
                    atom,
                    leaving_atom: (leaving_bookmark != 0)
                        .then(|| AtomId::new(leaving_bookmark as usize - 1)),
                    label: label.map(Into::into),
                    order: None,
                });
        }
        Ok(())
        // END RDKIT CPP FUNCTION
    }

    fn parse_scl(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SCLLine
        // RDKit❗✔️: void ParseSGroupV2000SCLLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit❗✔️:                              const std::string &text, unsigned int line,
        // RDKit❗✔️:                              bool strictParsing) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(text.substr(0, 6) == "M  SCL", "bad SCL line");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 6;
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:   if (!sgroup) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (pos + 1 >= text.length()) {
        // RDKit❗✔️:     std::ostringstream errout;
        // RDKit❗✔️:     errout << "SGroup SCL line too short: '" << text << "' on line " << line;
        // RDKit❗✔️:     SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   ++pos;
        // RDKit❗✔️:   sgroup->setProp("CLASS", text.substr(pos, text.length() - pos));
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.

        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        position += 1;
        if position >= line.len() {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "SGroup SCL line too short: '{line}' on line {line_number}"
                )),
            );
        }
        {
            let group = self
                .group_mut_if_present(sequence)
                .expect("presence checked before parsing SGroup class");
            group.set_prop("CLASS", &line[position..])?;
            group.set_class(&line[position..]);
        }
        Ok(())
    }

    fn parse_sbt(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SBTLine
        // RDKit❗✔️: void ParseSGroupV2000SBTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
        // RDKit❗✔️:                              const std::string &text, unsigned int line,
        // RDKit❗✔️:                              bool strictParsing) {
        // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
        // RDKit❗✔️:   PRECONDITION(text.substr(0, 6) == "M  SBT", "bad SBT line");
        // RDKit❗✔️:
        // RDKit❗✔️:   unsigned int pos = 6;
        // RDKit❗✔️:   bool ok;
        // RDKit❗✔️:   unsigned int nent =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   for (unsigned int ie = 0; ie < nent; ++ie) {
        // RDKit❗✔️:     if (text.size() < pos + 8) {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "SGroup SBT line too short: '" << text << "' on line " << line;
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     unsigned int sgIdx =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️:     if (!sgroup) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     unsigned int bracketType =
        // RDKit❗✔️:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:     if (!ok) {
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     if (bracketType == 0) {
        // RDKit❗✔️:       sgroup->setProp("BRKTYP", "BRACKET");
        // RDKit❗✔️:     } else if (bracketType == 1) {
        // RDKit❗✔️:       sgroup->setProp("BRKTYP", "PAREN");
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       std::ostringstream errout;
        // RDKit❗✔️:       errout << "Invalid SBT value '" << bracketType << "' on line " << line;
        // RDKit❗✔️:       SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION
        // Behavior review: retain source fixed fields, early exits and strict/
        // non-strict invalidation in detached staging; diagnostics/warning text
        // are not claimed as complete reproduction.
        // Cost review: each entry performs bounded fixed-field work and O(log S)
        // BTreeMap lookup like source std::map; owned label/data text remains
        // linear in copied bytes, with no topology or predicate clone.

        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SBT line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            let Some(bracket_type) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            let style = match bracket_type {
                0 => SGroupBracketStyle::Bracket,
                1 => SGroupBracketStyle::Parenthesis,
                _ => {
                    return self.invalidate_or_throw(
                        sequence,
                        SdfReadError::Parse(format!(
                            "Invalid SBT value '{bracket_type}' on line {line_number}"
                        )),
                    );
                }
            };
            let group = self
                .group_mut_if_present(sequence)
                .expect("presence checked before parsing SGroup bracket type");
            group.set_prop(
                "BRKTYP",
                if bracket_type == 0 {
                    "BRACKET"
                } else {
                    "PAREN"
                },
            )?;
            group.set_bracket_style(style);
        }
        Ok(())
    }

    pub(super) fn finish(
        mut self,
        bond_endpoints: &[(AtomId, AtomId)],
    ) -> Result<Vec<SubstanceGroup>, SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseMolBlockProperties (SGroup finalization)
        // RDKit❗✔️:     for (auto &sgroup : sGroupMap) {
        // RDKit❗✔️:       if (sgroup.second.getIsValid()) {
        // RDKit❗✔️:         sgroup.second.setProp("DATAFIELDS", dataFieldsMap[sgroup.first]);
        // RDKit❗✔️:         sgroup.second.setIsValid(checkAttachmentPointsAreValid(mol, sgroup));
        // RDKit❗✔️:   }
        // RDKit❗✔️:       if (sgroup.second.getIsValid()) {
        // RDKit❗✔️:         addSubstanceGroup(*mol, sgroup.second);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     if (strictParsing) {
        // RDKit❗✔️:       throw FileParseException(errout.str());
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       BOOST_LOG(rdWarningLog)
        // RDKit❗✔️:               << errout.str() << " and will be ignored" << std::endl;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // BEGIN RDKIT CPP FUNCTION checkAttachmentPointsAreValid
        // RDKit❗✔️:   for (auto &attachPoint : attachPoints) {
        // RDKit❗✔️:     if (attachPoint.lvIdx == nAtoms) {
        // RDKit❗✔️:       const std::vector<unsigned int> &bonds = sgroup.second.getBonds();
        // RDKit❗✔️:       if (bonds.size() == 1) {
        // RDKit❗✔️:         const auto bond = mol->getBondWithIdx(bonds.front());
        // RDKit❗✔️:         if (bond->getBeginAtomIdx() == attachPoint.aIdx ||
        // RDKit❗✔️:             bond->getEndAtomIdx() == attachPoint.aIdx) {
        // RDKit❗✔️:           attachPoint.lvIdx = bond->getOtherAtomIdx(attachPoint.aIdx);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:     if (attachPoint.lvIdx == nAtoms) {
        // RDKit❗✔️:     res = false;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        for (sequence, group) in &mut self.groups {
            if !self.invalid_sequences.contains(sequence) {
                // RDKit❗❌:         sgroup.second.setProp("DATAFIELDS", dataFieldsMap[sgroup.first]);
                // Actual tag/order is retained in the sole PropertyStore; the
                // existing explicit SGroupData field requires one extra O(D)
                // payload copy over source, retained as performance debt.
                group.set_prop(
                    "DATAFIELDS",
                    PropertyValue::StringVector(
                        group
                            .data()
                            .map_or_else(Vec::new, |data| data.values.clone()),
                    ),
                )?;
            }
        }
        for (sequence, atom, label) in std::mem::take(&mut self.pending_attach_points) {
            if self.invalid_sequences.contains(&sequence) {
                continue;
            }
            let leaving_atom = self.groups.get(&sequence).and_then(|group| {
                let [bond] = group.bonds() else {
                    return None;
                };
                let &(begin, end) = bond_endpoints.get(bond.index())?;
                if begin == atom {
                    Some(end)
                } else if end == atom {
                    Some(begin)
                } else {
                    None
                }
            });
            if let Some(leaving_atom) = leaving_atom {
                self.groups
                    .get_mut(&sequence)
                    .expect("pending attachment point belongs to a parsed group")
                    .push_attach_point(SGroupAttachPoint {
                        atom,
                        leaving_atom: Some(leaving_atom),
                        label: label.map(Into::into),
                        order: None,
                    });
            } else {
                self.invalid_sequences.insert(sequence);
            }
        }
        // END RDKIT CPP FUNCTION
        if let Some(sequence) = self.invalid_sequences.iter().next().copied()
            && self.strict_parsing
        {
            return Err(SdfReadError::Parse(format!("SGroup {sequence} is invalid")));
        }
        self.groups
            .retain(|sequence, _| !self.invalid_sequences.contains(sequence));
        let ids = self
            .groups
            .keys()
            .enumerate()
            .map(|(index, sequence)| (*sequence, SubstanceGroupId::new(index)))
            .collect::<BTreeMap<_, _>>();
        for (sequence, id) in &ids {
            self.groups
                .get_mut(sequence)
                .expect("key came from group map")
                .set_id(*id);
        }
        for (sequence, parent_sequence) in self.parent_by_sequence {
            let parent = *ids.get(&parent_sequence).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup {parent_sequence} referenced as parent not found."
                ))
            })?;
            self.groups
                .get_mut(&sequence)
                .ok_or_else(|| SdfReadError::Parse(format!("SGroup {sequence} missing")))?
                .set_parent(parent);
        }
        Ok(self.groups.into_values().collect())
        // END RDKIT CPP FUNCTION
    }
}

fn sgroup_string_value(
    group: &SubstanceGroup,
    key: &str,
    detached: Option<&[u8]>,
) -> Result<Option<PropertyText>, SdfWriteError> {
    // RDKit❗✔️:   template <typename T>
    // RDKit❗✔️:   bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit❗✔️:     return d_props.getValIfPresent(key, res);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   bool getValIfPresent(const std::string_view what, std::string &res) const {
    // RDKit❗✔️:     for (const auto &i : _data) {
    // RDKit❗✔️:       if (i.key == what) {
    // RDKit❗✔️:         rdvalue_tostring(i.val, res);
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // Actual counted property bytes take precedence. The explicit detached
    // field is the same modeled value when no source property was supplied.
    // O(log P + bytes), one canonical conversion; no decoding or retagging.
    match group.props().get(key.as_bytes()) {
        Some(value) => Ok(Some(cosmolkit_core::property_value_to_string(value)?)),
        None => Ok(detached.map(PropertyText::from_bytes)),
    }
}
fn sgroup_uint_value(
    group: &SubstanceGroup,
    key: &str,
    detached: Option<u32>,
) -> Result<Option<u32>, SdfWriteError> {
    // RDKit❗✔️:   template <typename T>
    // RDKit❗✔️:   bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit❗✔️:     return d_props.getValIfPresent(key, res);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   template <typename T>
    // RDKit❗✔️:   bool getValIfPresent(const std::string_view what, T &res) const {
    // RDKit❗✔️:     for (const auto &data : _data) {
    // RDKit❗✔️:       if (data.key == what) {
    // RDKit❗✔️:         res = from_rdvalue<T>(data.val);
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // Preserve from_rdvalue<unsigned int> conversion and failure order.
    // O(log P + bytes) for text, O(log P) for typed scalar, no buffering.
    match group.props().get(key.as_bytes()) {
        Some(value) => Ok(Some(cosmolkit_core::property_value_to_uint(value)?)),
        None => Ok(detached),
    }
}
fn sgroup_parent_value(group: &SubstanceGroup) -> Result<Option<u32>, SdfWriteError> {
    // RDKit❗✔️: std::string FormatV3000ParentBlock(const SubstanceGroup &sgroup) {
    // RDKit❗✔️:   std::ostringstream ret;
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int parentIdx = -1;
    // RDKit❗✔️:   if (sgroup.getPropIfPresent("PARENT", parentIdx)) {
    // RDKit❗✔️:     ret << " PARENT=" << parentIdx;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return ret.str();
    // RDKit❗✔️: }
    // Convert the explicit identifier only if the actual source property is
    // absent; a present property retains its own conversion/error ordering.
    if let Some(value) = group.props().get(b"PARENT".as_slice()) {
        return Ok(Some(cosmolkit_core::property_value_to_uint(value)?));
    }
    group
        .parent()
        .map(|parent| {
            parent
                .index()
                .checked_add(1)
                .and_then(|value| u32::try_from(value).ok())
                .ok_or_else(|| {
                    SdfWriteError::SubstanceGroup(
                        "parent sequence exceeds source unsigned int".into(),
                    )
                })
        })
        .transpose()
}
fn sgroup_data_values(group: &SubstanceGroup) -> Result<&[PropertyText], SdfWriteError> {
    // RDKit❗✔️:   template <typename T>
    // RDKit❗✔️:   bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit❗✔️:     return d_props.getValIfPresent(key, res);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   template <typename T>
    // RDKit❗✔️:   bool getValIfPresent(const std::string_view what, T &res) const {
    // RDKit❗✔️:     for (const auto &data : _data) {
    // RDKit❗✔️:       if (data.key == what) {
    // RDKit❗✔️:         res = from_rdvalue<T>(data.val);
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️: template <>
    // RDKit❗✔️: inline std::vector<std::string> rdvalue_cast<std::vector<std::string>>(
    // RDKit❗✔️:     RDValue_cast_t v) {
    // RDKit❗✔️:   if (rdvalue_is<std::vector<std::string>>(v)) {
    // RDKit❗✔️:     return *v.ptrCast<std::vector<std::string>>();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   throw std::bad_any_cast();
    // RDKit❗✔️: }
    // Wrong-kind errors remain visible; the absent explicit detached value
    // is borrowed. O(log P), no second property map or payload copy.
    match group.props().get(b"DATAFIELDS".as_slice()) {
        Some(value) => Ok(value.as_string_vector()?),
        None => Ok(group
            .data()
            .map_or(group.data_fields(), |data| data.values.as_slice())),
    }
}

fn rdkit_sgroup_type(group: &SubstanceGroup) -> Result<PropertyText, SdfWriteError> {
    if let Some(value) = group.props().get(b"TYPE".as_slice()) {
        return cosmolkit_core::property_value_to_string(value).map_err(SdfWriteError::Property);
    }
    // Behavior: absent typed property uses the existing explicit detached kind;
    // present source property uses the canonical byte-preserving converter.
    // Complexity: one indexed lookup plus source conversion, no second map.
    let kind: &[u8] = match group.kind() {
        SubstanceGroupKind::Data => b"DAT",
        SubstanceGroupKind::Superatom => b"SUP",
        SubstanceGroupKind::MultipleGroup => b"MUL",
        SubstanceGroupKind::StructuralRepeatUnit => b"SRU",
        SubstanceGroupKind::Monomer => b"MON",
        SubstanceGroupKind::Copolymer => b"COP",
        SubstanceGroupKind::Crosslink => b"CRO",
        SubstanceGroupKind::Graft => b"GRA",
        SubstanceGroupKind::Modification => b"MOD",
        SubstanceGroupKind::Mer => b"MER",
        SubstanceGroupKind::AnyPolymer => b"ANY",
        SubstanceGroupKind::MixtureComponent => b"COM",
        SubstanceGroupKind::Mixture => b"MIX",
        SubstanceGroupKind::Formulation => b"FOR",
        SubstanceGroupKind::Generic(value) => value.as_bytes(),
    };
    Ok(PropertyText::from(kind))
}

fn v3000_index_block<T>(
    name: &str,
    values: impl IntoIterator<Item = T>,
    index: impl Fn(T) -> usize,
) -> PropertyText {
    // BEGIN RDKIT CPP FUNCTION BuildV3000IdxVectorDataBlock
    // RDKit✔️✔️:   size_t size = dataVectorEnd - dataVectorBegin;
    // RDKit✔️✔️:   if (size) {
    // RDKit✔️✔️:     ret << ' ' << key << "=(" << size;
    // RDKit✔️✔️:     for (auto itr = dataVectorBegin; itr < dataVectorEnd; ++itr) {
    // RDKit✔️✔️:       ret << ' ' << 1 + *itr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:     ret << ')';
    // RDKit✔️✔️: }
    let indices = values.into_iter().map(index).collect::<Vec<_>>();
    if indices.is_empty() {
        return PropertyText::new();
    }
    let mut output = format!(" {name}=({}", indices.len());
    for value in indices {
        output.push_str(&format!(" {}", value + 1));
    }
    output.push(')');
    output.into()
    // END RDKIT CPP FUNCTION
}

fn v3000_string_block(name: &str, value: Option<&[u8]>) -> PropertyText {
    // RDKit❗✔️: std::string FormatV3000StringPropertyBlock(const std::string &prop,
    // RDKit❗✔️:                                            const SubstanceGroup &sgroup) {
    // RDKit❗✔️:   std::ostringstream ret;
    // RDKit❗✔️:
    // RDKit❗✔️:   std::string propValue;
    // RDKit❗✔️:   if (sgroup.getPropIfPresent(prop, propValue)) {
    // RDKit❗✔️:     if (!propValue.empty()) {
    // RDKit❗✔️:       ret << ' ' << prop << '=';
    // RDKit❗✔️:       // CTAB spec says: "Strings that contain blank spaces or start with left
    // RDKit❗✔️:       // parenthesis or double quote, must be surrounded by double quotes A
    // RDKit❗✔️:       // double quote can be entered literally by doubling it."
    // RDKit❗✔️:       // However, BIOVIA Draw 2020 doesn't correctly parse values like
    // RDKit❗✔️:       // foo"" or foo(bar) but does fine with "foo""" and "foo(bar)"
    // RDKit❗✔️:       // and both BIOVIA Draw and Marvin Sketch happily ignore the theoretically
    // RDKit❗✔️:       // extra quotes.
    // RDKit❗✔️:       bool needsQuotes = propValue.find(' ') != std::string::npos ||
    // RDKit❗✔️:                          propValue.find('"') != std::string::npos ||
    // RDKit❗✔️:                          propValue.find('(') != std::string::npos;
    // RDKit❗✔️:       if (needsQuotes) {
    // RDKit❗✔️:         ret << "\"";
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       for (auto chr : propValue) {
    // RDKit❗✔️:         ret << chr;
    // RDKit❗✔️:         // double quotes need to be doubled on output:
    // RDKit❗✔️:         if (chr == '"') {
    // RDKit❗✔️:           ret << chr;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       if (needsQuotes) {
    // RDKit❗✔️:         ret << "\"";
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return ret.str();
    // RDKit❗✔️: }

    let Some(value) = value.filter(|value| !value.is_empty()) else {
        return PropertyText::new();
    };
    // BEGIN RDKIT CPP FUNCTION FormatV3000StringPropertyBlock
    // RDKit✔️✔️:       bool needsQuotes = propValue.find(' ') != std::string::npos ||
    // RDKit✔️✔️:                          propValue.find('"') != std::string::npos ||
    // RDKit✔️✔️:                          propValue.find('(') != std::string::npos;
    // RDKit✔️✔️:       if (needsQuotes) {
    // RDKit✔️✔️:         ret << "\"";
    // RDKit✔️✔️: }
    // RDKit✔️✔️:       for (auto chr : propValue) {
    // RDKit✔️✔️:         ret << chr;
    // RDKit✔️✔️:         if (chr == '"') {
    // RDKit✔️✔️:         ret << chr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️:       if (needsQuotes) {
    // RDKit✔️✔️:         ret << "\"";
    // RDKit✔️✔️: }
    // Behavior: source byte scanning and quote doubling, without text decoding.
    // Complexity: O(bytes), one output buffer with source-required expansion.
    let needs_quotes = value.iter().any(|byte| matches!(byte, b' ' | b'"' | b'('));
    let mut output = PropertyText::from(format!(" {name}="));
    if needs_quotes {
        output.push_byte(b'"');
    }
    for &byte in value {
        output.push_byte(byte);
        if byte == b'"' {
            output.push_byte(byte);
        }
    }
    if needs_quotes {
        output.push_byte(b'"');
    }
    output
    // END RDKIT CPP FUNCTION
}

fn add_v3000_block(block: &[u8], current: &mut PropertyText, output: &mut PropertyText) {
    // RDKit❗✔️: void addBlockToSGroupString(std::string block, std::string &currentLine,
    // RDKit❗✔️:                             std::ostringstream &os) {
    // RDKit❗✔️:   if (block.empty()) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (currentLine.length() + block.length() < 78) {
    // RDKit❗✔️:     currentLine += block;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     os << currentLine << " -\n";
    // RDKit❗✔️:     unsigned int length = block.size();
    // RDKit❗✔️:     unsigned int start = 0;
    // RDKit❗✔️:     while (length - start >= 73) {
    // RDKit❗✔️:       os << "M  V30";
    // RDKit❗✔️:       if (start) {
    // RDKit❗✔️:         os << ' ';
    // RDKit❗✔️:       }
    // RDKit❗✔️:       os << block.substr(start, 72);
    // RDKit❗✔️:       start += 72;
    // RDKit❗✔️:       if (start < length) {
    // RDKit❗✔️:         // need to write more, so add another "-"
    // RDKit❗✔️:         os << "-\n";
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (start < length) {
    // RDKit❗✔️:       currentLine =
    // RDKit❗✔️:           "M  V30" + std::string(start ? " " : "") + block.substr(start, 73);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: }  // namespace

    // BEGIN RDKIT CPP FUNCTION addBlockToSGroupString
    // RDKit❗✔️:   if (block.empty()) {
    // RDKit❗✔️:   return;
    // RDKit❗✔️: }
    // RDKit❗✔️:   if (currentLine.length() + block.length() < 78) {
    // RDKit❗✔️:     currentLine += block;
    // RDKit❗✔️: } else {
    // RDKit❗✔️:     os << currentLine << " -\n";
    // RDKit❗✔️:     unsigned int length = block.size();
    // RDKit❗✔️:   unsigned int start = 0;
    // RDKit❗✔️:     while (length - start >= 73) {
    // RDKit❗✔️:       os << "M  V30";
    // RDKit❗✔️:       if (start) {
    // RDKit❗✔️:         os << ' ';
    // RDKit❗✔️:     }
    // RDKit❗✔️:       os << block.substr(start, 72);
    // RDKit❗✔️:       start += 72;
    // RDKit❗✔️:     if (start < length) {
    // RDKit❗✔️:         // need to write more, so add another "-"
    // RDKit❗✔️:         os << "-\n";
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Behavior: counted substrings split at source 72-byte offsets, including
    // non-UTF8 and embedded NUL. Generated syntax stays ASCII.
    // Complexity: O(bytes), append-only output and bounded slices, as source.
    if block.is_empty() {
        return;
    }
    if current.len() + block.len() < 78 {
        current.extend_bytes(block);
        return;
    }
    output.extend_bytes(current.as_bytes());
    output.extend_bytes(b" -\n");
    let mut start = 0;
    while block.len() - start >= 73 {
        output.extend_bytes(b"M  V30");
        if start != 0 {
            output.push_byte(b' ');
        }
        output.extend_bytes(&block[start..start + 72]);
        start += 72;
        if start < block.len() {
            output.extend_bytes(b"-\n");
        }
    }
    if start < block.len() {
        current.clear();
        current.extend_bytes(b"M  V30");
        if start != 0 {
            current.push_byte(b' ');
        }
        current.extend_bytes(&block[start..]);
    }
    // END RDKIT CPP FUNCTION
}

fn connection_text(value: &SGroupConnection) -> &[u8] {
    match value {
        SGroupConnection::HeadToHead => b"HH",
        SGroupConnection::HeadToTail => b"HT",
        SGroupConnection::Either => b"EU",
        SGroupConnection::Unknown(value) => value.as_bytes(),
    }
}

fn bracket_style_text(value: &SGroupBracketStyle) -> &[u8] {
    match value {
        SGroupBracketStyle::Bracket => b"BRACKET",
        SGroupBracketStyle::Parenthesis => b"PAREN",
        SGroupBracketStyle::None => b"",
        SGroupBracketStyle::Unknown(value) => value.as_bytes(),
    }
}

fn assigned_stereo_group_ids(groups: &[StereoGroup]) -> Vec<Option<u32>> {
    // BEGIN RDKIT CPP FUNCTION storeIdsInUse/assignMissingIds/assignStereoGroupIds
    // RDKit❗❗: void storeIdsInUse(boost::dynamic_bitset<> &ids, StereoGroup &sg) {
    // RDKit❗❗:   const auto groupId = sg.getWriteId();
    // RDKit❗❗:   if (groupId == 0) {
    // RDKit❗❗:     return;
    // RDKit❗❗:   } else if (groupId >= ids.size()) {
    // RDKit❗❗:     ids.resize(groupId + 1);
    // RDKit❗❗:   }
    // RDKit❗❗:   if (ids[groupId]) {
    // RDKit❗❗:     BOOST_LOG(rdWarningLog)
    // RDKit❗❗:         << "StereoGroup ID " << groupId
    // RDKit❗❗:         << " is used by more than one group, and will be reassined"
    // RDKit❗❗:         << std::endl;
    // RDKit❗❗:     sg.setWriteId(0);
    // RDKit❗❗:   } else {
    // RDKit❗❗:     ids[groupId] = true;
    // RDKit❗❗:   }
    // RDKit❗❗: }
    // RDKit❗❗: void assignMissingIds(const boost::dynamic_bitset<> &ids, unsigned &nextId,
    // RDKit❗❗:                       StereoGroup &sg) {
    // RDKit❗❗:   if (sg.getWriteId() == 0) {
    // RDKit❗❗:     ++nextId;
    // RDKit❗❗:     while (nextId < ids.size() && ids[nextId]) {
    // RDKit❗❗:       ++nextId;
    // RDKit❗❗:     }
    // RDKit❗❗:     sg.setWriteId(nextId);
    // RDKit❗❗:   }
    // RDKit❗❗: }
    // RDKit❗❗: void assignStereoGroupIds(std::vector<StereoGroup> &groups) {
    // RDKit❗❗:   if (groups.empty()) {
    // RDKit❗❗:     return;
    // RDKit❗❗:   }
    // RDKit❗❗:   boost::dynamic_bitset<> andIds;
    // RDKit❗❗:   boost::dynamic_bitset<> orIds;
    // RDKit❗❗:   for (auto &sg : groups) {
    // RDKit❗❗:     if (sg.getGroupType() == StereoGroupType::STEREO_AND) {
    // RDKit❗❗:       storeIdsInUse(andIds, sg);
    // RDKit❗❗:     } else if (sg.getGroupType() == StereoGroupType::STEREO_OR) {
    // RDKit❗❗:       storeIdsInUse(orIds, sg);
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // RDKit❗❗:   unsigned andId = 0;
    // RDKit❗❗:   unsigned orId = 0;
    // RDKit❗❗:   for (auto &sg : groups) {
    // RDKit❗❗:     if (sg.getGroupType() == StereoGroupType::STEREO_AND) {
    // RDKit❗❗:       assignMissingIds(andIds, andId, sg);
    // RDKit❗❗:     } else if (sg.getGroupType() == StereoGroupType::STEREO_OR) {
    // RDKit❗❗:       assignMissingIds(orIds, orId, sg);
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // RDKit❗❗: }
    // Behavior remains partial because the pinned duplicate-ID warning is not
    // represented by this private value-only helper.
    // Complexity review: the source bitsets use O(maximum ID) bits with O(1)
    // lookup; BTreeSet uses O(group count) entries with O(log group count)
    // lookup. It avoids sparse-ID-sized allocation but is slower for dense
    // inputs, so the performance axis remains unresolved without measurements.
    let mut assigned = groups
        .iter()
        .map(|group| match group.kind() {
            StereoGroupKind::Absolute => None,
            StereoGroupKind::Or | StereoGroupKind::And => {
                let id = group.write_id();
                (id != 0).then_some(id)
            }
        })
        .collect::<Vec<_>>();
    for kind in [StereoGroupKind::Or, StereoGroupKind::And] {
        let mut used = BTreeSet::new();
        for (group, id) in groups.iter().zip(&mut assigned) {
            if group.kind() == kind
                && let Some(value) = *id
                && !used.insert(value)
            {
                *id = None;
            }
        }
        let mut next = 0_u32;
        for (group, id) in groups.iter().zip(&mut assigned) {
            if group.kind() == kind && id.is_none() {
                next += 1;
                while used.contains(&next) {
                    next += 1;
                }
                *id = Some(next);
            }
        }
    }
    assigned
}

fn sgroup_bond_membership(
    group: &SubstanceGroup,
    bond: &cosmolkit_model::Bond,
) -> Result<SGroupBondRole, SdfWriteError> {
    // RDKit✔️✔️: SubstanceGroup::BondType SubstanceGroup::getBondType(
    // RDKit✔️✔️:     unsigned int bondIdx) const {
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       std::find(d_bonds.begin(), d_bonds.end(), bondIdx) != d_bonds.end(),
    // RDKit✔️✔️:       "bond is not part of the SubstanceGroup")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto bond = dp_mol->getBondWithIdx(bondIdx);
    // RDKit✔️✔️:   bool begin_atom_in_sgroup =
    // RDKit✔️✔️:       std::find(d_atoms.begin(), d_atoms.end(), bond->getBeginAtomIdx()) !=
    // RDKit✔️✔️:       d_atoms.end();
    // RDKit✔️✔️:   bool end_atom_in_sgroup = std::find(d_atoms.begin(), d_atoms.end(),
    // RDKit✔️✔️:                                       bond->getEndAtomIdx()) != d_atoms.end();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (begin_atom_in_sgroup && end_atom_in_sgroup) {
    // RDKit✔️✔️:     return SubstanceGroup::BondType::CBOND;
    // RDKit✔️✔️:   } else if (begin_atom_in_sgroup || end_atom_in_sgroup) {
    // RDKit✔️✔️:     return SubstanceGroup::BondType::XBOND;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Neither beginning nor ending atoms of bond " << bond->getIdx()
    // RDKit✔️✔️:            << " is in this SubstanceGroup.";
    // RDKit✔️✔️:     throw SubstanceGroupException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior review: callers visit the group's own validated bond rows;
    // endpoint membership defines the live role regardless of input label.
    // Neither endpoint in the group is a source error, not unsupported input.
    // Complexity review: two membership scans match the pinned vector scans.
    match (
        group.atoms().contains(&bond.begin()),
        group.atoms().contains(&bond.end()),
    ) {
        (true, true) => Ok(SGroupBondRole::Contained),
        (true, false) | (false, true) => Ok(SGroupBondRole::Crossing),
        (false, false) => Err(SdfWriteError::SGroupBondMembership {
            group: group.id(),
            bond: bond.id(),
        }),
    }
}

pub(super) fn write_v3000_sgroup(
    sequence: usize,
    group: &SubstanceGroup,
    bonds: &[cosmolkit_model::Bond],
) -> Result<PropertyText, SdfWriteError> {
    // BEGIN RDKIT CPP FUNCTION GetV3000MolFileSGroupLines
    // RDKit❗✔️:   std::string currLine = (boost::format("M  V30 %d %s %d") % idx %
    // RDKit❗✔️:                           sgroup.getProp<std::string>("TYPE") % id)
    // RDKit❗✔️:                              .str();
    // RDKit❗✔️:   addBlockToSGroupString(
    // RDKit❗✔️:       BuildV3000IdxVectorDataBlock("ATOMS", sgroup.getAtoms()), currLine, os);
    // RDKit✔️✔️: std::string BuildV3000BondsBlock(const SubstanceGroup &sgroup) {
    // RDKit✔️✔️:   std::ostringstream ret;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto isXBond = [&sgroup](unsigned int bondIdx) {
    // RDKit✔️✔️:     return SubstanceGroup::BondType::XBOND == sgroup.getBondType(bondIdx);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto bonds = sgroup.getBonds();
    // RDKit✔️✔️:   auto first_cbond = std::stable_partition(bonds.begin(), bonds.end(), isXBond);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   ret << BuildV3000IdxVectorDataBlock("XBONDS", bonds.begin(), first_cbond);
    // RDKit✔️✔️:   ret << BuildV3000IdxVectorDataBlock("CBONDS", first_cbond, bonds.end());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (sgroup.hasProp("XBHEAD")) {
    // RDKit✔️✔️:     auto v = sgroup.getProp<std::vector<unsigned int>>("XBHEAD");
    // RDKit✔️✔️:     ret << BuildV3000IdxVectorDataBlock("XBHEAD", v.begin(), v.end());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (sgroup.hasProp("XBCORR")) {
    // RDKit✔️✔️:     auto v = sgroup.getProp<std::vector<unsigned int>>("XBCORR");
    // RDKit✔️✔️:     ret << BuildV3000IdxVectorDataBlock("XBCORR", v.begin(), v.end());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return ret.str();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: addBlockToSGroupString(BuildV3000BondsBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(
    // RDKit❗✔️:     BuildV3000IdxVectorDataBlock("PATOMS", sgroup.getParentAtoms()), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("SUBTYPE", sgroup),
    // RDKit❗✔️:                        currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("CONNECT", sgroup),
    // RDKit❗✔️:                        currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000ParentBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000CompNoBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("LABEL", sgroup),
    // RDKit❗✔️:                        currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000BracketBlock(sgroup.getBrackets()),
    // RDKit❗✔️:                        currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000CStateBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000FieldDataBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000AttachPointBlock(sgroup.getAttachPoints()),
    // RDKit❗✔️:                        currLine, os);
    let mut output = PropertyText::new();
    let mut current = PropertyText::from(format!("M  V30 {sequence} "));
    current.extend_bytes(rdkit_sgroup_type(group)?.as_bytes());
    current.extend_bytes(
        format!(
            " {}",
            sgroup_uint_value(group, "ID", group.external_id())?.unwrap_or(0)
        )
        .as_bytes(),
    );
    let mut blocks = vec![v3000_index_block(
        "ATOMS",
        group.atoms().iter().copied(),
        AtomId::index,
    )];
    // Stable classification preserves each category's input order, using the
    // current topology even if group atom membership changed after parsing.
    let mut crossing = Vec::new();
    let mut contained = Vec::new();
    for id in group.bonds() {
        let bond = &bonds[id.index()];
        match sgroup_bond_membership(group, bond)? {
            SGroupBondRole::Crossing => crossing.push(*id),
            SGroupBondRole::Contained => contained.push(*id),
        }
    }
    blocks.push(v3000_index_block(
        "XBONDS",
        crossing.into_iter(),
        BondId::index,
    ));
    blocks.push(v3000_index_block(
        "CBONDS",
        contained.into_iter(),
        BondId::index,
    ));
    blocks.push(v3000_index_block(
        "XBHEAD",
        group.head_crossing_bonds().iter().copied(),
        BondId::index,
    ));
    blocks.push(v3000_index_block(
        "XBCORR",
        group.crossing_bond_correspondence().iter().copied(),
        BondId::index,
    ));
    blocks.push(v3000_index_block(
        "PATOMS",
        group.parent_atoms().iter().copied(),
        AtomId::index,
    ));
    blocks.push(v3000_string_block(
        "SUBTYPE",
        sgroup_string_value(
            group,
            "SUBTYPE",
            group.subtype().map(PropertyText::as_bytes),
        )?
        .as_ref()
        .map(PropertyText::as_bytes),
    ));
    blocks.push(v3000_string_block(
        "MULT",
        group
            .props()
            .get(b"MULT".as_slice())
            .map(cosmolkit_core::property_value_to_string)
            .transpose()
            .map_err(SdfWriteError::Property)?
            .as_ref()
            .map(PropertyText::as_bytes),
    ));
    blocks.push(v3000_string_block(
        "CONNECT",
        sgroup_string_value(group, "CONNECT", group.connection().map(connection_text))?
            .as_ref()
            .map(PropertyText::as_bytes),
    ));
    if let Some(parent) = sgroup_parent_value(group)? {
        blocks.push(format!(" PARENT={parent}").into());
    }
    if let Some(component) = sgroup_uint_value(group, "COMPNO", group.component_number())? {
        blocks.push(format!(" COMPNO={component}").into());
    }
    blocks.push(v3000_string_block(
        "LABEL",
        sgroup_string_value(group, "LABEL", group.label().map(PropertyText::as_bytes))?
            .as_ref()
            .map(PropertyText::as_bytes),
    ));
    if let Some(display) = group.display() {
        for bracket in &display.brackets {
            blocks.push(
                format!(
                    " BRKXYZ=(9 {:.4} {:.4} 0 {:.4} {:.4} 0 0 0 0)",
                    bracket.points[0][0],
                    bracket.points[0][1],
                    bracket.points[1][0],
                    bracket.points[1][1]
                )
                .into(),
            );
        }
    }
    blocks.push(v3000_string_block(
        "ESTATE",
        sgroup_string_value(
            group,
            "ESTATE",
            group.expansion_state().map(PropertyText::as_bytes),
        )?
        .as_ref()
        .map(PropertyText::as_bytes),
    ));
    for cstate in group.cstates() {
        blocks.push(
            (if group.kind() == &SubstanceGroupKind::Superatom {
                format!(
                    " CSTATE=(4 {} {:.4} {:.4} 0)",
                    cstate.bond.index() + 1,
                    cstate.vector[0],
                    cstate.vector[1]
                )
            } else {
                format!(" CSTATE=(1 {})", cstate.bond.index() + 1)
            })
            .into(),
        );
    }
    let data = group.data();
    blocks.push(v3000_string_block(
        "FIELDNAME",
        sgroup_string_value(
            group,
            "FIELDNAME",
            data.and_then(|data| data.field_name.as_ref().map(PropertyText::as_bytes)),
        )?
        .as_ref()
        .map(PropertyText::as_bytes),
    ));
    blocks.push(v3000_string_block(
        "FIELDINFO",
        sgroup_string_value(
            group,
            "FIELDINFO",
            data.and_then(|data| data.field_info.as_ref().map(PropertyText::as_bytes)),
        )?
        .as_ref()
        .map(PropertyText::as_bytes),
    ));
    blocks.push(v3000_string_block(
        "FIELDDISP",
        sgroup_string_value(
            group,
            "FIELDDISP",
            data.and_then(|data| data.field_display.as_ref().map(PropertyText::as_bytes)),
        )?
        .as_ref()
        .map(PropertyText::as_bytes),
    ));
    blocks.push(v3000_string_block(
        "QUERYTYPE",
        sgroup_string_value(
            group,
            "QUERYTYPE",
            data.and_then(|data| data.query_type.as_ref().map(PropertyText::as_bytes)),
        )?
        .as_ref()
        .map(PropertyText::as_bytes),
    ));
    blocks.push(v3000_string_block(
        "QUERYOP",
        sgroup_string_value(
            group,
            "QUERYOP",
            data.and_then(|data| data.query_op.as_ref().map(PropertyText::as_bytes)),
        )?
        .as_ref()
        .map(PropertyText::as_bytes),
    ));
    // RDKit✔️✔️: std::string FormatV3000FieldDataBlock(const SubstanceGroup &sgroup) {
    // RDKit✔️✔️:   std::ostringstream ret;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   STR_VECT dataFields;
    // RDKit✔️✔️:   if (sgroup.getPropIfPresent("DATAFIELDS", dataFields)) {
    // RDKit✔️✔️:     for (const auto &data : dataFields) {
    // RDKit✔️✔️:       ret << " FIELDDATA=\"" << data << "\"";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return ret.str();
    // RDKit✔️✔️: }
    // Copy every payload byte once, including quotes; the source does not
    // escape DATAFIELDS. Linear append cost, with no decoded payload copy.
    for value in sgroup_data_values(group)? {
        let mut block = PropertyText::from(" FIELDDATA=\"");
        for &byte in value.as_bytes() {
            block.push_byte(byte);
        }
        block.push_byte(b'"');
        blocks.push(block);
    }
    blocks.push(v3000_string_block(
        "CLASS",
        sgroup_string_value(group, "CLASS", group.class().map(PropertyText::as_bytes))?
            .as_ref()
            .map(PropertyText::as_bytes),
    ));
    for point in group.attach_points() {
        let leaving = if point.leaving_atom == Some(point.atom) {
            "aidx".to_owned()
        } else {
            point
                .leaving_atom
                .map_or_else(|| "0".to_owned(), |atom| (atom.index() + 1).to_string())
        };
        let mut block =
            PropertyText::from(format!(" SAP=(3 {} {leaving} ", point.atom.index() + 1));
        block.extend_bytes(
            point
                .label
                .as_ref()
                .map(PropertyText::as_bytes)
                .unwrap_or_default(),
        );
        block.push_byte(b')');
        blocks.push(block);
    }
    blocks.push(v3000_string_block(
        "BRKTYP",
        sgroup_string_value(
            group,
            "BRKTYP",
            group.bracket_style().map(bracket_style_text),
        )?
        .as_ref()
        .map(PropertyText::as_bytes),
    ));
    for block in blocks {
        add_v3000_block(block.as_bytes(), &mut current, &mut output);
    }
    if !current.is_empty() {
        output.extend_bytes((&current).as_bytes());
        output.push_byte(b'\n');
    } else if output.as_bytes().ends_with(b" -\n") {
        output = PropertyText::from(&output.as_bytes()[..output.len() - 3]);
        output.push_byte(b'\n');
    }
    Ok(output)
    // END RDKIT CPP FUNCTION
}

pub(super) fn write_v3000_typed_blocks(
    topology: &TopologyBlock,
) -> Result<PropertyText, SdfWriteError> {
    let mut output = PropertyText::new();
    if !topology.substance_groups.is_empty() {
        // BEGIN RDKIT CPP FUNCTION getV3000CTAB (SGroup block)
        // RDKit✔️✔️: if (nSGroups > 0) {
        // RDKit✔️✔️:   res += "M  V30 BEGIN SGROUP\n";
        // RDKit✔️✔️:   unsigned int idx = 0;
        // RDKit✔️✔️:   for (const auto &sgroup : sgroups) {
        // RDKit✔️✔️:     res += GetV3000MolFileSGroupLines(++idx, sgroup);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   res += "M  V30 END SGROUP\n";
        // RDKit✔️✔️: }
        output.extend_bytes(("M  V30 BEGIN SGROUP\n").as_bytes());
        for (index, group) in topology.substance_groups.iter().enumerate() {
            output
                .extend_bytes((write_v3000_sgroup(index + 1, group, &topology.bonds)?).as_bytes());
        }
        output.extend_bytes(("M  V30 END SGROUP\n").as_bytes());
        // END RDKIT CPP FUNCTION
    }
    let atoms = topology
        .stereo_groups
        .iter()
        .map(|group| group.atoms().to_vec())
        .collect::<Vec<_>>();
    output.extend_bytes(write_v3000_collection_rows(&topology.stereo_groups, &atoms).as_bytes());

    Ok(output)
}

fn v2000_int(value: usize) -> String {
    format!(" {value:>3}")
}

fn v2000_count(value: usize) -> String {
    format!(" {value:>2}")
}

fn v2000_string_field(value: &[u8], width: usize, pad: bool, separator: bool) -> PropertyText {
    // RDKit❗✔️: inline std::string FormatV2000StringField(const std::string &value,
    // RDKit❗✔️:                                           unsigned int fieldSize, bool pad,
    // RDKit❗✔️:                                           bool addSeparator) {
    // RDKit❗✔️:   std::ostringstream os;
    // RDKit❗✔️:   if (addSeparator) {
    // RDKit❗✔️:     os << ' ';
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (value.size() >= fieldSize) {
    // RDKit❗✔️:     os << value.substr(0, fieldSize);
    // RDKit❗✔️:   } else if (pad) {
    // RDKit❗✔️:     os << std::setw(fieldSize) << std::left << value;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     os << value;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return os.str();
    // RDKit❗✔️: }

    // RDKit❗✔️: inline std::string FormatV2000StringField(const std::string &value,
    // RDKit❗✔️:                                           unsigned int fieldSize, bool pad,
    // RDKit❗✔️:                                           bool addSeparator) {
    // RDKit❗✔️:   std::ostringstream os;
    // RDKit❗✔️:   if (addSeparator) {
    // RDKit❗✔️:     os << ' ';
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (value.size() >= fieldSize) {
    // RDKit❗✔️:     os << value.substr(0, fieldSize);
    // RDKit❗✔️:   } else if (pad) {
    // RDKit❗✔️:     os << std::setw(fieldSize) << std::left << value;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     os << value;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return os.str();
    // RDKit❗✔️: }
    // RDKit❗✔️:
    // RDKit❗✔️: inline std::string FormatV3000DoubleField(double value) {
    // Behavior: byte-counted prefix/truncation and ASCII padding, no decoding.
    // Complexity: one output buffer, O(min(bytes,width)+padding), as source.
    let mut output = PropertyText::new();
    if separator {
        output.push_byte(b' ');
    }
    output.extend_bytes(&value[..value.len().min(width)]);
    if pad {
        for _ in value.len()..width {
            output.push_byte(b' ');
        }
    }
    output
}

fn append_v2000_pairs(
    output: &mut PropertyText,
    code: &str,
    pairs: &[(usize, PropertyText)],
    per_line: usize,
) {
    for chunk in pairs.chunks(per_line) {
        output.extend_bytes((&format!("M  {code}{}", v2000_count(chunk.len()))).as_bytes());
        for (index, value) in chunk {
            output.extend_bytes((&v2000_int(*index)).as_bytes());
            output.push_byte(b' ');
            output.extend_bytes((value).as_bytes());
        }
        output.push_byte(b'\n');
    }
}

fn append_v2000_indices(
    output: &mut PropertyText,
    code: &str,
    group: usize,
    indices: impl IntoIterator<Item = usize>,
) {
    let values = indices.into_iter().collect::<Vec<_>>();
    for chunk in values.chunks(15) {
        output.extend_bytes(
            (&format!("M  {code}{}{}", v2000_int(group), v2000_count(chunk.len()))).as_bytes(),
        );
        for index in chunk {
            output.extend_bytes((&v2000_int(index + 1)).as_bytes());
        }
        output.push_byte(b'\n');
    }
}

pub(super) fn write_v2000_sgroups(topology: &TopologyBlock) -> Result<PropertyText, SdfWriteError> {
    if topology.substance_groups.is_empty() {
        return Ok(PropertyText::new());
    }
    // BEGIN RDKIT CPP FUNCTION GetMolFileSGroupInfo
    // RDKit❗✔️:   ret << BuildV2000STYLines(mol);
    // RDKit❗✔️:   ret << BuildV2000SLBLines(mol);
    // RDKit❗✔️:   ret << BuildV2000StringPropLines(8, mol, "SUBTYPE", "SST", 3);
    // RDKit❗✔️:   ret << BuildV2000StringPropLines(8, mol, "CONNECT", "SCN", 3);
    // RDKit❗✔️:   ret << BuildV2000SDSLines(mol);
    // RDKit❗✔️:   ret << BuildV2000SPLLines(mol);
    // RDKit❗✔️:   ret << BuildV2000SNCLines(mol);
    // RDKit❗✔️:   ret << BuildV2000SBTLines(mol);
    // RDKit❗✔️:   for (const auto &sgroup : getSubstanceGroups(mol)) {
    // RDKit❗✔️:     ret << BuildV2000IdxVectorDataLines(15, idx, "SAL", sgroup.getAtoms());
    // RDKit❗✔️:     ret << BuildV2000IdxVectorDataLines(15, idx, "SPA",
    // RDKit❗✔️:                                         sgroup.getParentAtoms());
    // RDKit❗✔️:     ret << BuildV2000IdxVectorDataLines(15, idx, "SBL", sgroup.getBonds());
    // RDKit❗✔️:     ret << BuildV2000SDILine(idx, sgroup);
    // RDKit❗✔️:     ret << BuildV2000SMTLine(idx, sgroup);
    // RDKit❗✔️:     ret << BuildV2000SBVLine(idx, sgroup);
    // RDKit❗✔️:     ret << BuildV2000SDTLine(idx, sgroup);
    // RDKit❗✔️:     ret << BuildV2000SDDLine(idx, sgroup);
    // RDKit❗✔️:     ret << BuildV2000SCDSEDLines(idx, sgroup);
    // RDKit❗✔️:     ret << BuildV2000SAPLines(idx, sgroup);
    // RDKit❗✔️:     ret << BuildV2000SCLLine(idx, sgroup);
    // RDKit❗✔️: }
    let mut output = PropertyText::new();
    for chunk in topology.substance_groups.chunks(8).enumerate() {
        output.extend_bytes(format!("M  STY{}", v2000_count(chunk.1.len())).as_bytes());
        for (offset, group) in chunk.1.iter().enumerate() {
            output.extend_bytes(v2000_int(chunk.0 * 8 + offset + 1).as_bytes());
            output.push_byte(b' ');
            let kind = rdkit_sgroup_type(group)?;
            output.extend_bytes(&kind.as_bytes()[..kind.len().min(3)]);
            for _ in kind.len()..3 {
                output.push_byte(b' ');
            }
        }
        output.push_byte(b'\n');
    }
    // RDKit❗✔️: std::string BuildV2000SLBLines(const ROMol &mol) {
    // RDKit❗✔️:   std::ostringstream ret;
    // RDKit❗✔️:   std::ostringstream temp;
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   const auto &sgroups = getSubstanceGroups(mol);
    // RDKit❗✔️:   for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    // RDKit❗✔️:     unsigned int id;
    // RDKit❗✔️:     // Write value if assigned, else 0
    // RDKit❗✔️:     if (sg->getPropIfPresent("ID", id)) {
    // RDKit❗✔️:       temp << FormatV2000IntField(1 + (sg - sgroups.begin()))
    // RDKit❗✔️:            << FormatV2000IntField(id);
    // RDKit❗✔️:       if (++count == 8) {
    // RDKit❗✔️:         ret << "M  SLB" << FormatV2000NumEntriesField(8) << temp.str() << "\n";
    // RDKit❗✔️:         temp.str("");
    // RDKit❗✔️:         count = 0;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (count) {
    // RDKit❗✔️:     ret << "M  SLB" << FormatV2000NumEntriesField(count) << temp.str() << "\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return ret.str();
    // RDKit❗✔️: }
    // RDKit❗✔️: std::string BuildV2000StringPropLines(const unsigned int entriesPerLine,
    // RDKit❗✔️:                                       const ROMol &mol,
    // RDKit❗✔️:                                       const std::string &propName,
    // RDKit❗✔️:                                       const std::string &propCode,
    // RDKit❗✔️:                                       const unsigned int fieldWitdh) {
    // RDKit❗✔️:   std::ostringstream ret;
    // RDKit❗✔️:   std::ostringstream temp;
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   const auto &sgroups = getSubstanceGroups(mol);
    // RDKit❗✔️:   for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    // RDKit❗✔️:     std::string propValue;
    // RDKit❗✔️:     // Write field only if defined
    // RDKit❗✔️:     if (sg->getPropIfPresent(propName, propValue)) {
    // RDKit❗✔️:       temp << FormatV2000IntField(1 + (sg - sgroups.begin()))
    // RDKit❗✔️:            << FormatV2000StringField(propValue, fieldWitdh, true, true);
    // RDKit❗✔️:       if (++count == entriesPerLine) {
    // RDKit❗✔️:         ret << "M  " << propCode << FormatV2000NumEntriesField(entriesPerLine)
    // RDKit❗✔️:             << temp.str() << "\n";
    // RDKit❗✔️:         temp.str("");
    // RDKit❗✔️:         count = 0;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (count) {
    // RDKit❗✔️:     ret << "M  " << propCode << FormatV2000NumEntriesField(count) << temp.str()
    // RDKit❗✔️:         << "\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return ret.str();
    // RDKit❗✔️: }
    // RDKit❗✔️: std::string BuildV2000SDSLines(const ROMol &mol) {
    // RDKit❗✔️:   std::ostringstream ret;
    // RDKit❗✔️:   std::ostringstream temp;
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   const auto &sgroups = getSubstanceGroups(mol);
    // RDKit❗✔️:   for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    // RDKit❗✔️:     // Write field only if defined
    // RDKit❗✔️:     std::string eState;
    // RDKit❗✔️:     if (sg->getPropIfPresent("ESTATE", eState) && eState == "E") {
    // RDKit❗✔️:       temp << FormatV2000IntField(1 + (sg - sgroups.begin()));
    // RDKit❗✔️:       if (++count == 15) {
    // RDKit❗✔️:         ret << "M  SDS EXP" << FormatV2000NumEntriesField(15) << temp.str()
    // RDKit❗✔️:             << "\n";
    // RDKit❗✔️:         temp.str("");
    // RDKit❗✔️:         count = 0;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (count) {
    // RDKit❗✔️:     ret << "M  SDS EXP" << FormatV2000NumEntriesField(count) << temp.str()
    // RDKit❗✔️:         << "\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return ret.str();
    // RDKit❗✔️: }
    // RDKit❗✔️: std::string BuildV2000SPLLines(const ROMol &mol) {
    // RDKit❗✔️:   std::ostringstream ret;
    // RDKit❗✔️:   std::ostringstream temp;
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   const auto &sgroups = getSubstanceGroups(mol);
    // RDKit❗✔️:   for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    // RDKit❗✔️:     // Write field only if a parent is defined
    // RDKit❗✔️:     unsigned int parentIdx = -1;
    // RDKit❗✔️:     if (sg->getPropIfPresent("PARENT", parentIdx)) {
    // RDKit❗✔️:       temp << FormatV2000IntField(1 + (sg - sgroups.begin()))
    // RDKit❗✔️:            << FormatV2000IntField(parentIdx);
    // RDKit❗✔️:       if (++count == 8) {
    // RDKit❗✔️:         ret << "M  SPL" << FormatV2000NumEntriesField(8) << temp.str() << "\n";
    // RDKit❗✔️:         temp.str("");
    // RDKit❗✔️:         count = 0;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (count) {
    // RDKit❗✔️:     ret << "M  SPL" << FormatV2000NumEntriesField(count) << temp.str() << "\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return ret.str();
    // RDKit❗✔️: }
    // RDKit❗✔️: std::string BuildV2000SNCLines(const ROMol &mol) {
    // RDKit❗✔️:   std::ostringstream ret;
    // RDKit❗✔️:   std::ostringstream temp;
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   const auto &sgroups = getSubstanceGroups(mol);
    // RDKit❗✔️:   for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    // RDKit❗✔️:     unsigned int compno;
    // RDKit❗✔️:     // Write field only if compno is set
    // RDKit❗✔️:     if (sg->getPropIfPresent("COMPNO", compno)) {
    // RDKit❗✔️:       temp << FormatV2000IntField(1 + (sg - sgroups.begin()))
    // RDKit❗✔️:            << FormatV2000IntField(compno);
    // RDKit❗✔️:       if (++count == 8) {
    // RDKit❗✔️:         ret << "M  SNC" << FormatV2000NumEntriesField(8) << temp.str() << "\n";
    // RDKit❗✔️:         temp.str("");
    // RDKit❗✔️:         count = 0;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (count) {
    // RDKit❗✔️:     ret << "M  SNC" << FormatV2000NumEntriesField(count) << temp.str() << "\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return ret.str();
    // RDKit❗✔️: }
    // RDKit❗✔️: std::string BuildV2000SBTLines(const ROMol &mol) {
    // RDKit❗✔️:   std::ostringstream ret;
    // RDKit❗✔️:   std::ostringstream temp;
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int count = 0;
    // RDKit❗✔️:   const auto &sgroups = getSubstanceGroups(mol);
    // RDKit❗✔️:   for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    // RDKit❗✔️:     std::string bracketType;
    // RDKit❗✔️:     if (sg->getPropIfPresent("BRKTYP", bracketType)) {
    // RDKit❗✔️:       unsigned int idx = 1 + (sg - sgroups.begin());
    // RDKit❗✔️:       if (bracketType == "BRACKET") {
    // RDKit❗✔️:         temp << FormatV2000IntField(idx) << FormatV2000IntField(0);
    // RDKit❗✔️:       } else if (bracketType == "PAREN") {
    // RDKit❗✔️:         temp << FormatV2000IntField(idx) << FormatV2000IntField(1);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         std::ostringstream errout;
    // RDKit❗✔️:         errout << "Invalid BRKTYP value '" << bracketType << "' for SGroup "
    // RDKit❗✔️:                << idx;
    // RDKit❗✔️:         throw SubstanceGroupException(errout.str());
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (++count == 8) {
    // RDKit❗✔️:         ret << "M  SBT" << FormatV2000NumEntriesField(8) << temp.str() << "\n";
    // RDKit❗✔️:         temp.str("");
    // RDKit❗✔️:         count = 0;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (count) {
    // RDKit❗✔️:     ret << "M  SBT" << FormatV2000NumEntriesField(count) << temp.str() << "\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return ret.str();
    // RDKit❗✔️: }
    let mut slb = Vec::new();
    for (index, group) in topology.substance_groups.iter().enumerate() {
        if let Some(id) = sgroup_uint_value(group, "ID", group.external_id())? {
            // Source FormatV2000IntField(int) receives unsigned int by the
            // pinned platform's narrowing conversion, preserving its bits.
            slb.push((index + 1, PropertyText::from(format!("{:>3}", id as i32))));
        }
    }
    append_v2000_pairs(&mut output, "SLB", &slb, 8);
    let mut subtype = Vec::new();
    for (index, group) in topology.substance_groups.iter().enumerate() {
        if let Some(value) = sgroup_string_value(
            group,
            "SUBTYPE",
            group.subtype().map(PropertyText::as_bytes),
        )? {
            subtype.push((
                index + 1,
                v2000_string_field(value.as_bytes(), 3, true, false),
            ));
        }
    }
    append_v2000_pairs(&mut output, "SST", &subtype, 8);
    let mut connections = Vec::new();
    for (index, group) in topology.substance_groups.iter().enumerate() {
        if let Some(value) =
            sgroup_string_value(group, "CONNECT", group.connection().map(connection_text))?
        {
            connections.push((
                index + 1,
                v2000_string_field(value.as_bytes(), 3, true, false),
            ));
        }
    }
    append_v2000_pairs(&mut output, "SCN", &connections, 8);
    let mut expanded = Vec::new();
    for (index, group) in topology.substance_groups.iter().enumerate() {
        if sgroup_string_value(
            group,
            "ESTATE",
            group.expansion_state().map(PropertyText::as_bytes),
        )?
        .is_some_and(|value| value.as_bytes() == b"E")
        {
            expanded.push(index + 1);
        }
    }
    for chunk in expanded.chunks(15) {
        output.extend_bytes(format!("M  SDS EXP{}", v2000_count(chunk.len())).as_bytes());
        for index in chunk {
            output.extend_bytes(v2000_int(*index).as_bytes());
        }
        output.push_byte(b'\n');
    }
    let mut parents = Vec::new();
    for (index, group) in topology.substance_groups.iter().enumerate() {
        if let Some(parent) = sgroup_parent_value(group)? {
            parents.push((
                index + 1,
                PropertyText::from(format!("{:>3}", parent as i32)),
            ));
        }
    }
    append_v2000_pairs(&mut output, "SPL", &parents, 8);
    let mut components = Vec::new();
    for (index, group) in topology.substance_groups.iter().enumerate() {
        if let Some(value) = sgroup_uint_value(group, "COMPNO", group.component_number())? {
            components.push((
                index + 1,
                PropertyText::from(format!("{:>3}", value as i32)),
            ));
        }
    }
    append_v2000_pairs(&mut output, "SNC", &components, 8);
    let mut brackets = Vec::new();
    for (index, group) in topology.substance_groups.iter().enumerate() {
        let Some(style) = sgroup_string_value(
            group,
            "BRKTYP",
            group.bracket_style().map(bracket_style_text),
        )?
        else {
            continue;
        };
        let value = match style.as_bytes() {
            b"BRACKET" => 0,
            b"PAREN" => 1,
            _ => {
                return Err(SdfWriteError::SubstanceGroup(format!(
                    "Invalid BRKTYP value {:?} for SGroup {}",
                    style,
                    index + 1
                )));
            }
        };
        brackets.push((index + 1, PropertyText::from(format!("{value:>3}"))));
    }
    append_v2000_pairs(&mut output, "SBT", &brackets, 8);

    for (zero_index, group) in topology.substance_groups.iter().enumerate() {
        let index = zero_index + 1;
        append_v2000_indices(
            &mut output,
            "SAL",
            index,
            group.atoms().iter().map(|atom| atom.index()),
        );
        append_v2000_indices(
            &mut output,
            "SPA",
            index,
            group.parent_atoms().iter().map(|atom| atom.index()),
        );
        append_v2000_indices(
            &mut output,
            "SBL",
            index,
            group.bonds().iter().map(|bond| bond.index()),
        );
        if let Some(display) = group.display() {
            for bracket in &display.brackets {
                output.extend_bytes(
                    (&format!(
                        "M  SDI{}{}{:>10.4}{:>10.4}{:>10.4}{:>10.4}\n",
                        v2000_int(index),
                        v2000_count(4),
                        bracket.points[0][0],
                        bracket.points[0][1],
                        bracket.points[1][0],
                        bracket.points[1][1]
                    ))
                        .as_bytes(),
                );
            }
        }
        // RDKit❗✔️: std::string BuildV2000SMTLine(const int idx, const SubstanceGroup &sgroup) {
        // RDKit❗✔️:   std::ostringstream ret;
        // RDKit❗✔️:
        // RDKit❗✔️:   std::string smtValue;
        // RDKit❗✔️:   if ((sgroup.getProp<std::string>("TYPE") == "MUL" &&
        // RDKit❗✔️:        sgroup.getPropIfPresent("MULT", smtValue)) ||
        // RDKit❗✔️:       sgroup.getPropIfPresent("LABEL", smtValue)) {
        // RDKit❗✔️:     ret << "M  SMT" << FormatV2000IntField(idx)
        // RDKit❗✔️:         << FormatV2000StringField(smtValue, 69, false, true) << "\n";
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return ret.str();
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: std::string BuildV2000SDILine(const int idx, const SubstanceGroup &sgroup) {
        let multiple = if rdkit_sgroup_type(group)?.as_bytes() == b"MUL" {
            group
                .props()
                .get(b"MULT".as_slice())
                .map(cosmolkit_core::property_value_to_string)
                .transpose()
                .map_err(SdfWriteError::Property)?
        } else {
            None
        };
        let label = match multiple {
            Some(value) => Some(value),
            None => sgroup_string_value(group, "LABEL", group.label().map(PropertyText::as_bytes))?,
        };
        if let Some(label) = label {
            let label = label.as_bytes();
            output.extend_bytes(format!("M  SMT{} ", v2000_int(index)).as_bytes());
            output.extend_bytes(&label[..label.len().min(69)]);
            output.push_byte(b'\n');
        }
        for cstate in group.cstates() {
            output.extend_bytes(
                (&format!(
                    "M  SBV{}{}",
                    v2000_int(index),
                    v2000_int(cstate.bond.index() + 1)
                ))
                    .as_bytes(),
            );
            if group.kind() == &SubstanceGroupKind::Superatom {
                output.extend_bytes(
                    (&format!("{:>10.4}{:>10.4}", cstate.vector[0], cstate.vector[1])).as_bytes(),
                );
            }
            output.push_byte(b'\n');
        }
        let data = group.data();
        let field_name = sgroup_string_value(
            group,
            "FIELDNAME",
            data.and_then(|data| data.field_name.as_ref().map(PropertyText::as_bytes)),
        )?;
        // BEGIN RDKIT CPP FUNCTION BuildV2000SDTLine
        // RDKit❗✔️: std::string BuildV2000SDTLine(const int idx, const SubstanceGroup &sgroup) {
        // RDKit❗✔️:   std::ostringstream ret;
        // RDKit❗✔️:
        // RDKit❗✔️:   std::string sdtValue;
        // RDKit❗✔️:   if (sgroup.getPropIfPresent("FIELDNAME", sdtValue)) {
        // RDKit❗✔️:     ret << "M  SDT" << FormatV2000IntField(idx);
        // RDKit❗✔️:     ret << FormatV2000StringField(sdtValue, 30, true, true);
        // RDKit❗✔️:
        // RDKit❗✔️:     if (sgroup.getPropIfPresent("FIELDTYPE", sdtValue)) {
        // RDKit❗✔️:       ret << FormatV2000StringField(sdtValue, 2, true, false);
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       ret << " T";
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     if (sgroup.getPropIfPresent("FIELDINFO", sdtValue)) {
        // RDKit❗✔️:       ret << FormatV2000StringField(sdtValue, 20, true, false);
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     if (sgroup.getPropIfPresent("QUERYTYPE", sdtValue)) {
        // RDKit❗✔️:       ret << FormatV2000StringField(sdtValue, 2, true, false);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (sgroup.getPropIfPresent("QUERYOP", sdtValue)) {
        // RDKit❗✔️:       ret << FormatV2000StringField(sdtValue, 15, true, false);
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     ret << "\n";
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return ret.str();
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION
        if let Some(field_name) = field_name {
            output.extend_bytes(format!("M  SDT{}", v2000_int(index)).as_bytes());
            output
                .extend_bytes(v2000_string_field(field_name.as_bytes(), 30, true, true).as_bytes());
            let field_type = sgroup_string_value(
                group,
                "FIELDTYPE",
                data.and_then(|d| d.field_type.as_ref().map(PropertyText::as_bytes)),
            )?;
            if let Some(field_type) = field_type {
                output.extend_bytes(
                    v2000_string_field(field_type.as_bytes(), 2, true, false).as_bytes(),
                );
            } else {
                output.extend_bytes(b" T");
            }
            for (key, field, width) in [
                ("FIELDINFO", data.and_then(|d| d.field_info.as_ref()), 20),
                ("QUERYTYPE", data.and_then(|d| d.query_type.as_ref()), 2),
                ("QUERYOP", data.and_then(|d| d.query_op.as_ref()), 15),
            ] {
                if let Some(value) =
                    sgroup_string_value(group, key, field.map(PropertyText::as_bytes))?
                {
                    output.extend_bytes(
                        v2000_string_field(value.as_bytes(), width, true, false).as_bytes(),
                    );
                }
            }
            output.push_byte(b'\n');
        }
        // Preserve source field presence independently from empty bytes, and
        // source QUERYOP truncation/padding at fifteen counted bytes. Four
        // optional lookups and bounded byte formatting match the source cost.
        if let Some(display) = sgroup_string_value(
            group,
            "FIELDDISP",
            data.and_then(|d| d.field_display.as_ref().map(PropertyText::as_bytes)),
        )? {
            output.extend_bytes(format!("M  SDD{} ", v2000_int(index)).as_bytes());
            output.extend_bytes(display.as_bytes());
            output.push_byte(b'\n');
        }
        for value in sgroup_data_values(group)? {
            if value.len() > 200 {
                return Err(SdfWriteError::SubstanceGroup(format!(
                    "data field in SGroup {index} is longer than 200 bytes"
                )));
            }
            if value.is_empty() {
                output.extend_bytes(format!("M  SED{} \n", v2000_int(index)).as_bytes());
            } else {
                let count = value.len().div_ceil(69);
                for (i, chunk) in value.as_bytes().chunks(69).enumerate() {
                    let code = if i + 1 == count { "SED" } else { "SCD" };
                    output.extend_bytes(format!("M  {code}{} ", v2000_int(index)).as_bytes());
                    output.extend_bytes(chunk);
                    output.push_byte(b'\n');
                }
            }
        }
        for chunk in group.attach_points().chunks(6) {
            output.extend_bytes(
                (&format!("M  SAP{}{}", v2000_int(index), v2000_count(chunk.len()))).as_bytes(),
            );
            for point in chunk {
                output.extend_bytes((&v2000_int(point.atom.index() + 1)).as_bytes());
                output.extend_bytes(
                    (&v2000_int(point.leaving_atom.map_or(0, |atom| atom.index() + 1))).as_bytes(),
                );
                output.push_byte(b' ');
                output.extend_bytes(
                    v2000_string_field(
                        point
                            .label
                            .as_ref()
                            .map(PropertyText::as_bytes)
                            .unwrap_or_default(),
                        2,
                        true,
                        false,
                    )
                    .as_bytes(),
                );
            }
            output.push_byte(b'\n');
        }
        if let Some(class) =
            sgroup_string_value(group, "CLASS", group.class().map(PropertyText::as_bytes))?
        {
            output.extend_bytes(format!("M  SCL{} ", v2000_int(index)).as_bytes());
            output.extend_bytes(class.as_bytes());
            output.push_byte(b'\n');
        }
    }
    Ok(output)
    // END RDKIT CPP FUNCTION
}

#[cfg(test)]
mod cf3d_sgids_io_4_tests {
    use super::{parse_stereo_collection_line, parse_v3000_collection_block};
    use crate::sdf::SdfReadError;
    use cosmolkit_model::{AtomId, StereoGroupKind};

    #[test]
    fn cf3d_sgids_io_4_collection_numbers_are_read_ids_only() {
        let lines = [
            "M  V30 MDLV30/HILITE ATOMS=(2 1 2)",
            "M  V30 MDLV30/STEABS ATOMS=(2 3 1)",
            "M  V30 MDLV30/STEREL17 ATOMS=(2 2 4)",
            "M  V30 MDLV30/STERAC23 ATOMS=(1 3)",
            "M  V30 END COLLECTION",
        ];
        let mut cursor = 0;
        let progress = crate::sdf::MolBlockReadProgress::default();
        let lines = crate::sdf::MolBlockLines::new(&lines, &progress);
        let groups = parse_v3000_collection_block(&lines, &mut cursor, 4, true).unwrap();

        assert_eq!(cursor, lines.len());
        assert_eq!(groups.len(), 3);
        assert_eq!(groups[0].kind(), StereoGroupKind::Absolute);
        assert_eq!(groups[0].id(), Some(0));
        assert_eq!(groups[0].write_id(), 0);
        assert_eq!(groups[0].atoms(), &[AtomId::new(2), AtomId::new(0)]);

        assert_eq!(groups[1].kind(), StereoGroupKind::Or);
        assert_eq!(groups[1].id(), Some(17));
        assert_eq!(groups[1].write_id(), 0);
        assert_eq!(groups[1].atoms(), &[AtomId::new(1), AtomId::new(3)]);

        assert_eq!(groups[2].kind(), StereoGroupKind::And);
        assert_eq!(groups[2].id(), Some(23));
        assert_eq!(groups[2].write_id(), 0);
        assert_eq!(groups[2].atoms(), &[AtomId::new(2)]);
    }

    #[test]
    fn cf3d_sgids_io_4_invalid_members_keep_source_error() {
        let error = parse_stereo_collection_line("MDLV30/STEREL17 ATOMS=(1 5)", 8, 4).unwrap_err();
        assert!(matches!(
            error,
            SdfReadError::Parse(message)
                if message == "Stereo group atom index 5 out of range on line 8"
        ));
    }
}

#[cfg(test)]
mod cf3d_sgids_io_5_tests {
    use super::{assigned_stereo_group_ids, write_v3000_typed_blocks};
    use cosmolkit_model::{Atom, AtomId, AtomSpec, StereoGroup, StereoGroupKind, TopologyBlock};
    use cosmolkit_types::Element;

    fn group(kind: StereoGroupKind, atom: usize, read_id: u32, write_id: u32) -> StereoGroup {
        StereoGroup::new(kind, vec![AtomId::new(atom)], vec![])
            .with_id(read_id)
            .with_write_id(write_id)
    }

    #[test]
    fn cf3d_sgids_io_5_writer_uses_write_ids_and_keeps_groups_unchanged() {
        let stereo_groups = vec![
            group(StereoGroupKind::Or, 0, 17, 0),
            group(StereoGroupKind::Or, 1, 17, 9),
            group(StereoGroupKind::Or, 2, 77, 9),
            group(StereoGroupKind::Or, 3, 4, 4),
            group(StereoGroupKind::And, 4, 17, 0),
            group(StereoGroupKind::And, 5, 17, 4),
            group(StereoGroupKind::Absolute, 6, 17, 0),
        ];
        let atoms = (0..7)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let topology = TopologyBlock::try_from_parts(atoms, vec![], vec![], stereo_groups)
            .expect("fixed writer topology is valid");
        let original_groups = topology.stereo_groups.clone();

        assert_eq!(
            assigned_stereo_group_ids(&topology.stereo_groups),
            vec![Some(1), Some(9), Some(2), Some(4), Some(1), Some(4), None]
        );
        assert_eq!(
            write_v3000_typed_blocks(&topology)
                .map(super::fixture_writer_text)
                .unwrap(),
            "M  V30 BEGIN COLLECTION\n\
             M  V30 MDLV30/STEREL1 ATOMS=(1 1)\n\
             M  V30 MDLV30/STEREL9 ATOMS=(1 2)\n\
             M  V30 MDLV30/STEREL2 ATOMS=(1 3)\n\
             M  V30 MDLV30/STEREL4 ATOMS=(1 4)\n\
             M  V30 MDLV30/STERAC1 ATOMS=(1 5)\n\
             M  V30 MDLV30/STERAC4 ATOMS=(1 6)\n\
             M  V30 MDLV30/STEABS ATOMS=(1 7)\n\
             M  V30 END COLLECTION\n"
        );
        assert_eq!(topology.stereo_groups, original_groups);
    }
}

pub(super) fn write_v3000_collection_rows(
    groups: &[StereoGroup],
    atom_ids: &[Vec<cosmolkit_model::AtomId>],
) -> PropertyText {
    // RDKit❗❌: void appendEnhancedStereoGroups(
    // RDKit❗❌:     std::string &res, const RWMol &tmol,
    // RDKit❗❌:     std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> &wedgeBonds) {
    // RDKit❗❌:   if (!tmol.getStereoGroups().empty()) {
    // RDKit❗❌:     auto stereo_groups = tmol.getStereoGroups();
    // RDKit❗❌:     assignStereoGroupIds(stereo_groups);
    // RDKit❗❌:     res += "M  V30 BEGIN COLLECTION\n";
    // RDKit❗❌:     std::string tmp;
    // RDKit❗❌:     tmp.reserve(80);
    // RDKit❗❌:     for (auto &&group : stereo_groups) {
    // RDKit❗❌:       tmp += "M  V30 MDLV30/";
    // RDKit❗❌:       switch (group.getGroupType()) {
    // RDKit❗❌:         case RDKit::StereoGroupType::STEREO_ABSOLUTE:
    // RDKit❗❌:           tmp += "STEABS";
    // RDKit❗❌:           break;
    // RDKit❗❌:         case RDKit::StereoGroupType::STEREO_OR:
    // RDKit❗❌:           tmp += "STEREL";
    // RDKit❗❌:           tmp += std::to_string(group.getWriteId());
    // RDKit❗❌:           break;
    // RDKit❗❌:         case RDKit::StereoGroupType::STEREO_AND:
    // RDKit❗❌:           tmp += "STERAC";
    // RDKit❗❌:           tmp += std::to_string(group.getWriteId());
    // RDKit❗❌:           break;
    // RDKit❗❌:       }
    // RDKit❗❌:       tmp += " ATOMS=(";
    // RDKit❗❌:
    // RDKit❗❌:       std::vector<unsigned int> atomIds;
    // RDKit❗❌:       Atropisomers::getAllAtomIdsForStereoGroup(tmol, group, atomIds,
    // RDKit❗❌:                                                 wedgeBonds);
    // RDKit❗❌:
    // RDKit❗❌:       tmp += std::to_string(atomIds.size());
    // RDKit❗❌:       for (auto &&atom : atomIds) {
    // RDKit❗❌:         tmp += ' ';
    // RDKit❗❌:         // atoms are 1 indexed in molfiles
    // RDKit❗❌:         auto idxStr = std::to_string(atom + 1);
    // RDKit❗❌:         if (tmp.size() + idxStr.size() >= 78) {
    // RDKit❗❌:           res += tmp + "-\n";
    // RDKit❗❌:           tmp = "M  V30 ";
    // RDKit❗❌:         }
    // RDKit❗❌:         tmp += idxStr;
    // RDKit❗❌:       }
    // RDKit❗❌:       res += tmp + ")\n";
    // RDKit❗❌:       tmp.clear();
    // RDKit❗❌:     }
    // RDKit❗❌:     res += tmp + "M  V30 END COLLECTION\n";
    // RDKit❗❌:   }
    // RDKit❗❌: }
    let mut output = PropertyText::new();
    if !groups.is_empty() {
        // BEGIN RDKIT CPP FUNCTION appendEnhancedStereoGroups
        // RDKit❗✔️: auto stereo_groups = tmol.getStereoGroups();
        // RDKit❗✔️: assignStereoGroupIds(stereo_groups);
        // RDKit❗✔️: res += "M  V30 BEGIN COLLECTION\n";
        // RDKit❗✔️: std::string tmp;
        // RDKit❗✔️: tmp.reserve(80);
        // RDKit❗✔️: for (auto &&group : stereo_groups) {
        // RDKit❗✔️:   tmp += "M  V30 MDLV30/";
        // RDKit❗✔️:   switch (group.getGroupType()) {
        // RDKit❗✔️:     case RDKit::StereoGroupType::STEREO_ABSOLUTE:
        // RDKit❗✔️:       tmp += "STEABS";
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case RDKit::StereoGroupType::STEREO_OR:
        // RDKit❗✔️:       tmp += "STEREL";
        // RDKit❗✔️:       tmp += std::to_string(group.getWriteId());
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case RDKit::StereoGroupType::STEREO_AND:
        // RDKit❗✔️:       tmp += "STERAC";
        // RDKit❗✔️:       tmp += std::to_string(group.getWriteId());
        // RDKit❗✔️:       break;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   tmp += " ATOMS=(";
        // RDKit❗✔️:   tmp += std::to_string(atomIds.size());
        // RDKit❗✔️:   for (auto &&atom : atomIds) {
        // RDKit❗✔️:     tmp += ' ';
        // RDKit❗✔️:     // atoms are 1 indexed in molfiles
        // RDKit❗✔️:     auto idxStr = std::to_string(atom + 1);
        // RDKit❗✔️:     if (tmp.size() + idxStr.size() >= 78) {
        // RDKit❗✔️:       res += tmp + "-\n";
        // RDKit❗✔️:       tmp = "M  V30 ";
        // RDKit❗✔️:     }
        // RDKit❗✔️:     tmp += idxStr;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   res += tmp + ")\n";
        // RDKit❗✔️:   tmp.clear();
        // RDKit❗✔️: }
        // RDKit❗✔️: res += tmp + "M  V30 END COLLECTION\n";
        output.extend_bytes(b"M  V30 BEGIN COLLECTION\n");
        let group_ids = assigned_stereo_group_ids(&groups);
        for ((group, assigned_id), atoms) in groups.iter().zip(group_ids).zip(atom_ids) {
            let label = match group.kind() {
                StereoGroupKind::Absolute => "STEABS".to_owned(),
                StereoGroupKind::Or => format!("STEREL{}", assigned_id.expect("OR ID assigned")),
                StereoGroupKind::And => {
                    format!("STERAC{}", assigned_id.expect("AND ID assigned"))
                }
            };
            let mut current = format!("M  V30 MDLV30/{label} ATOMS=({}", atoms.len());
            for atom in atoms {
                current.push(' ');
                let index = (atom.index() + 1).to_string();
                if current.len() + index.len() >= 78 {
                    output.extend_bytes(current.as_bytes());
                    output.extend_bytes(b"-\n");
                    current.clear();
                    current.push_str("M  V30 ");
                }
                current.push_str(&index);
            }
            output.extend_bytes(current.as_bytes());
            output.extend_bytes(b")\n");
        }
        output.extend_bytes(b"M  V30 END COLLECTION\n");
        // END RDKIT CPP FUNCTION
    }
    output
}

#[cfg(test)]
fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(value.as_bytes())
        .expect("original text fixture must retain exact UTF-8 bytes")
}
#[cfg(test)]
fn fixture_writer_text(value: cosmolkit_model::PropertyText) -> String {
    String::from_utf8(value.into_bytes())
        .expect("original writer fixture must retain exact UTF-8 bytes")
}

#[cfg(test)]
mod native_sdt_field_tests {
    use super::*;

    #[test]
    fn v3000_field_data_retains_quotes_verbatim_like_rdkit() {
        let group = SubstanceGroup::new(
            cosmolkit_model::SubstanceGroupId::new(0),
            SubstanceGroupKind::Data,
        )
        .with_data(cosmolkit_model::SGroupData {
            values: vec!["a\"b".into(), "".into(), "a\"\"b".into()],
            ..Default::default()
        });
        let output = write_v3000_sgroup(1, &group, &[]).unwrap();
        assert_eq!(
            output.as_bytes(),
            b"M  V30 1 DAT 0 FIELDDATA=\"a\"b\" FIELDDATA=\"\" FIELDDATA=\"a\"\"b\"\n"
        );
    }

    fn sdt(data: cosmolkit_model::SGroupData) -> Vec<u8> {
        let mut topology = TopologyBlock::default();
        topology.substance_groups.push(
            SubstanceGroup::new(
                cosmolkit_model::SubstanceGroupId::new(0),
                SubstanceGroupKind::Data,
            )
            .with_data(data),
        );
        let output = write_v2000_sgroups(&topology).unwrap();
        output
            .as_bytes()
            .split(|b| *b == b'\n')
            .find(|l| l.starts_with(b"M  SDT"))
            .unwrap()
            .to_vec()
    }

    #[test]
    fn native_sdt_omits_absent_fields_and_distinguishes_empty_field_type() {
        for field_type in [None, Some(PropertyText::new())] {
            let actual = sdt(cosmolkit_model::SGroupData {
                field_name: Some("FIELD".into()),
                field_type: field_type.clone(),
                ..Default::default()
            });
            let mut expected = b"M  SDT   1 FIELD                         ".to_vec();
            expected.extend_from_slice(if field_type.is_none() { b" T" } else { b"  " });
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn native_sdt_queryop_preserves_fifteen_byte_padding_and_truncation() {
        for (input, expected_tail) in [
            (b"OP".as_slice(), b"OP             ".as_slice()),
            (
                b"abcdefghijklmnopq".as_slice(),
                b"abcdefghijklmno".as_slice(),
            ),
            (
                &[0xff, 0x00][..],
                &[
                    0xff, 0x00, b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ', b' ',
                    b' ', b' ',
                ][..],
            ),
        ] {
            let actual = sdt(cosmolkit_model::SGroupData {
                field_name: Some("FIELD".into()),
                query_op: Some(PropertyText::from(input.to_vec())),
                ..Default::default()
            });
            let mut expected = b"M  SDT   1 FIELD                          T".to_vec();
            expected.extend_from_slice(expected_tail);
            assert_eq!(actual, expected);
        }
    }
}
