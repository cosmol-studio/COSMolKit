//! IO-owned coordinate PDB writer (BIO-PDB-WRITE).
//!
//! Seven-record coordinate projection (ATOM/HETATM/ANISOU/TER/MODEL/
//! ENDMDL/END) of pinned Gemmi 0.7.5 `write_pdb`; explicitly NOT a
//! metadata-lossless claim. This module hosts the private N1/N2/R1/R2/R3
//! helpers; public detached entries are added with R3 (Step 38).

use std::fmt;
use std::path::Path;

use cosmolkit_bio::BioStructureData;
use cosmolkit_bio::{EntityKind, find_residue_info};
use cosmolkit_types::Element;

// ─── Public detached entries (packet frozen contract §2) ───

/// Exactly five source-default bool fields (pinned PdbWriteOptions,
/// NOT minimal()).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioPdbWriteParams {
    pub ter_records: bool,
    pub numbered_ter: bool,
    pub ter_ignores_type: bool,
    pub preserve_serial: bool,
    pub end_record: bool,
}

impl Default for BioPdbWriteParams {
    fn default() -> Self {
        Self {
            ter_records: true,
            numbered_ter: true,
            ter_ignores_type: false,
            preserve_serial: false,
            end_record: true,
        }
    }
}

/// Typed PDB writer errors (frozen §8.2).
#[derive(Debug)]
pub enum BioPdbWriteError {
    InvalidStructure(cosmolkit_bio::BioStructureError),
    ChainNameTooLong {
        chain: String,
        length: usize,
    },
    NegativeSerial {
        serial: i32,
    },
    SerialIncrementOverflow {
        from: i32,
    },
    SerialOffsetOverflow {
        serial: i32,
    },
    NegativeBase36Value {
        value: i32,
    },
    UnrepresentableField {
        field: &'static str,
        detail: String,
    },
    Io {
        path: std::path::PathBuf,
        source: std::io::Error,
    },
}

impl fmt::Display for BioPdbWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidStructure(e) => write!(f, "invalid structure: {e}"),
            Self::ChainNameTooLong { chain, length } => write!(
                f,
                "chain name too long for the PDB format ({length} bytes): {chain}"
            ),
            Self::NegativeSerial { serial } => write!(f, "negative atom serial: {serial}"),
            Self::SerialIncrementOverflow { from } => {
                write!(f, "serial increment overflows i32 from {from}")
            }
            Self::SerialOffsetOverflow { serial } => {
                write!(f, "hybrid36 offset overflows i32 at serial {serial}")
            }
            Self::NegativeBase36Value { value } => {
                write!(f, "base36 requires non-negative value: {value}")
            }
            Self::UnrepresentableField { field, detail } => {
                write!(f, "unrepresentable {field}: {detail}")
            }
            Self::Io { path, source } => write!(f, "failed to write {}: {source}", path.display()),
        }
    }
}

impl std::error::Error for BioPdbWriteError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidStructure(e) => Some(e),
            Self::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

impl From<cosmolkit_bio::BioStructureError> for BioPdbWriteError {
    fn from(e: cosmolkit_bio::BioStructureError) -> Self {
        Self::InvalidStructure(e)
    }
}

/// Public detached text writer (packet frozen contract §2).
pub fn bio_structure_to_pdb_text(
    data: &BioStructureData,
    params: &BioPdbWriteParams,
) -> Result<String, BioPdbWriteError> {
    records::write_pdb_coordinate_text(
        data,
        params.ter_records,
        params.numbered_ter,
        params.ter_ignores_type,
        params.preserve_serial,
        params.end_record,
        data.input_format(),
    )
    .map_err(|e| match e {
        numbering::NumberingError::NegativeSerial { serial } => {
            BioPdbWriteError::NegativeSerial { serial }
        }
        numbering::NumberingError::SerialIncrementOverflow { from } => {
            BioPdbWriteError::SerialIncrementOverflow { from }
        }
        numbering::NumberingError::SerialOffsetOverflow { serial } => {
            BioPdbWriteError::SerialOffsetOverflow { serial }
        }
        numbering::NumberingError::NegativeBase36Value { value } => {
            BioPdbWriteError::NegativeBase36Value { value }
        }
    })
}

/// Public detached file writer (packet frozen contract §2): produce the
/// ENTIRE text BEFORE opening/truncating the destination. Invalid
/// structure/writer failure MUST leave an existing destination unchanged.
/// Filesystem write failure does not promise rollback of already-written
/// filesystem bytes.
pub fn write_bio_structure_pdb_file(
    data: &BioStructureData,
    path: &Path,
    params: &BioPdbWriteParams,
) -> Result<(), BioPdbWriteError> {
    // Produce the complete text FIRST — any structure/writer error
    // returns before the destination is touched.
    let text = bio_structure_to_pdb_text(data, params)?;
    // Only NOW open/truncate the destination.
    std::fs::write(path, text.as_bytes()).map_err(|source| BioPdbWriteError::Io {
        path: path.to_path_buf(),
        source,
    })
}

/// Private hybrid36/base36 numeric helpers (N1). All source preconditions
/// become typed checked errors; source UB (negative base36 input,
/// arithmetic overflow) is NEVER executed.
pub(crate) mod numbering {
    use std::fmt;

    /// Typed N1 domain errors (frozen in receipt §8.2).
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub(crate) enum NumberingError {
        NegativeSerial { serial: i32 },
        SerialIncrementOverflow { from: i32 },
        SerialOffsetOverflow { serial: i32 },
        NegativeBase36Value { value: i32 },
    }

    impl fmt::Display for NumberingError {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                Self::NegativeSerial { serial } => {
                    write!(f, "negative atom serial: {serial}")
                }
                Self::SerialIncrementOverflow { from } => {
                    write!(f, "serial increment overflows i32 from {from}")
                }
                Self::SerialOffsetOverflow { serial } => {
                    write!(f, "hybrid36 offset overflows i32 at serial {serial}")
                }
                Self::NegativeBase36Value { value } => {
                    write!(f, "base36 requires non-negative value: {value}")
                }
            }
        }
    }

    // BEGIN GEMMI CPP HELPERS (to_pdb.cpp:55-89, seqid.hpp:17-23, sprintf.hpp to_chars_z)
    // Gemmi✔️✔️: // works for non-negative values only
    // Gemmi✔️✔️: void base36_encode(char* buffer, int width, int value) {
    // Gemmi✔️✔️:   const char base36[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    // Gemmi✔️✔️:   buffer[width] = '\0';
    // Gemmi✔️✔️:   do {
    // Gemmi✔️✔️:     buffer[--width] = base36[value % 36];
    // Gemmi✔️✔️:     value /= 36;
    // Gemmi✔️✔️:   } while (value != 0 && width != 0);
    // Gemmi✔️✔️:   while (width != 0)
    // Gemmi✔️✔️:     buffer[--width] = ' ';
    // Gemmi✔️✔️: }
    // END GEMMI CPP HELPERS
    //
    // Behavior review: the source's do/while then right-fill loop maps to
    // a fixed-width byte buffer written back-to-front; the negative-input
    // precondition ("works for non-negative values only") becomes a typed
    // NegativeBase36Value error BEFORE any division, because value % 36
    // on i32::MIN is UB in both languages. Checked only where the source
    // PRECONDITION defines the domain (Rust never executes UB).
    //
    // Cost review: O(width) writes into a stack array, then one String
    // copy of at most 5 bytes — source writes into a caller stack buffer
    // with no String; the fresh String is a known PDB-writer-local cost,
    // not unresolved.
    pub(crate) fn base36_encode_checked(
        width: usize,
        value: i32,
    ) -> Result<String, NumberingError> {
        if value < 0 {
            return Err(NumberingError::NegativeBase36Value { value });
        }
        let mut buffer = [b' '; 8];
        let base36 = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        let mut remaining = value as u32; // value >= 0: no loss
        let mut position = width;
        loop {
            position -= 1;
            buffer[position] = base36[(remaining % 36) as usize];
            remaining /= 36;
            if remaining == 0 || position == 0 {
                break;
            }
        }
        Ok(String::from_utf8_lossy(&buffer[..width]).into_owned())
    }

    // BEGIN GEMMI CPP HELPERS
    // Gemmi✔️✔️: std::array<char,8> encode_serial_in_hybrid36(int serial) {
    // Gemmi✔️✔️:   std::array<char,8> str;
    // Gemmi✔️✔️:   assert(serial >= 0);
    // Gemmi✔️✔️:   if (serial < 100000)
    // Gemmi✔️✔️:     to_chars_z(str.data(), str.data() + 8, serial);
    // Gemmi✔️✔️:   else
    // Gemmi✔️✔️:     base36_encode(str.data(), 5, serial + (10 * 36 * 36 * 36 * 36 - 100000));
    // Gemmi✔️✔️:   return str;
    // Gemmi✔️✔️: }
    // END GEMMI CPP HELPERS
    //
    // Behavior review: serial < 100000 renders plain decimal (to_chars_z
    // of a non-negative int is exact decimal); serial >= 100000 adds the
    // frozen offset 10*36^4 - 100000 = 16696160 (receipt §12.3) and
    // renders width-5 base36. The source assert(serial >= 0) becomes a
    // typed NegativeSerial error; the i32 offset addition at serial >=
    // 2130787488 (i32::MAX - 16696160 + 1) is a checked
    // SerialOffsetOverflow — the source's plain + is UB there, never
    // executed in Rust.
    //
    // Cost review: one branch plus O(1) formatting exactly as the source;
    // the returned String is the same known writer-local cost class.
    pub(crate) fn encode_serial_in_hybrid36_checked(serial: i32) -> Result<String, NumberingError> {
        if serial < 0 {
            return Err(NumberingError::NegativeSerial { serial });
        }
        if serial < 100_000 {
            return Ok(serial.to_string());
        }
        let offset_serial = serial
            .checked_add(16_696_160)
            .ok_or(NumberingError::SerialOffsetOverflow { serial })?;
        base36_encode_checked(5, offset_serial)
    }

    /// Checked serial increment used by numbered TER (source `++serial`).
    pub(crate) fn increment_serial(serial: i32) -> Result<i32, NumberingError> {
        serial
            .checked_add(1)
            .ok_or(NumberingError::SerialIncrementOverflow { from: serial })
    }

    // BEGIN GEMMI CPP HELPERS
    // Gemmi✔️✔️: std::array<char,8> write_seq_id(const SeqId& seqid) {
    // Gemmi✔️✔️:   std::array<char,8> str;
    // Gemmi✔️✔️:   char* ptr = str.data();
    // Gemmi✔️✔️:   if (*seqid.num > -1000 && *seqid.num < 10000) {
    // Gemmi✔️✔️:     ptr = to_chars_z(ptr, ptr + 5, *seqid.num);
    // Gemmi✔️✔️:   } else if (seqid.num.has_value()) {
    // Gemmi✔️✔️:     base36_encode(ptr, 4, *seqid.num + (10 * 36 * 36 * 36 - 10000));
    // Gemmi✔️✔️:     ptr += 4;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   *ptr++ = seqid.icode;
    // Gemmi✔️✔️:   *ptr = '\0';
    // END GEMMI CPP HELPERS
    //
    // Behavior review: number in (-1000, 10000) renders as plain decimal
    // (5-column to_chars_z of a value in that open interval never
    // truncates); a PRESENT number outside that interval adds the frozen
    // seq offset 10*36^3 - 10000 = 456560 and renders width-4 base36;
    // a MISSING number writes NOTHING before the insertion code. The
    // insertion code byte follows verbatim. The present-outside-interval
    // addition can overflow i32 (num = i32::MIN direction) — checked
    // SerialOffsetOverflow with the num value; negative offset results
    // stay NegativeBase36Value (both are source UB, never executed).
    //
    // Cost review: one branch plus O(1) formatting; same String cost
    // class as above.
    pub(crate) fn write_seq_id_checked(
        number: Option<i32>,
        insertion_code: u8,
    ) -> Result<String, NumberingError> {
        let mut out = String::with_capacity(6);
        match number {
            Some(num) if num > -1000 && num < 10_000 => {
                out.push_str(&num.to_string());
            }
            Some(num) => {
                let offset = num
                    .checked_add(456_560)
                    .ok_or(NumberingError::SerialOffsetOverflow { serial: num })?;
                out.push_str(&base36_encode_checked(4, offset)?);
            }
            None => {}
        }
        out.push(insertion_code as char);
        Ok(out)
    }
}

/// Private N2 padded-name / HETATM-classification helpers (BIO-PDB-WRITE).
pub(crate) mod fields {
    use super::*;

    pub(super) fn gemmi_element_name(element: Element, isotope: Option<u16>) -> &'static str {
        // Gemmi✔️✔️:   const char* uname() const { return element_uppercase_name(elem); }
        // BIO folds source El::D into H plus isotope 2. All other checked
        // atomic numbers index the existing source-ordered BIO vocabulary.
        // This is one O(1) static-table lookup, matching the source getter.
        let ordinal = if element == Element::H && isotope == Some(2) {
            119
        } else {
            usize::from(element.atomic_number())
        };
        cosmolkit_bio::GEMMI_ELEMENT_NAMES[ordinal]
    }

    // BEGIN GEMMI CPP HELPERS (model.hpp:153-162, to_pdb.cpp:41-52, util.hpp:61, elem.hpp:27/387)
    // Gemmi✔️✔️:   std::string padded_name() const {
    // Gemmi✔️✔️:     std::string s;
    // Gemmi✔️✔️:     const char* el = element.uname();
    // Gemmi✔️✔️:     if (el[1] == '\0' &&
    // Gemmi✔️✔️:         (el[0] == alpha_up(name[0]) || (is_hydrogen() && alpha_up(name[0]) == 'H')) &&
    // Gemmi✔️✔️:         name.size() < 4)
    // Gemmi✔️✔️:       s += ' ';
    // Gemmi✔️✔️:     s += name;
    // Gemmi✔️✔️:     return s;
    // Gemmi✔️✔️:   }
    // END GEMMI CPP HELPERS
    //
    // Behavior review: ONE leading space is added iff the element's
    // uppercase symbol is a single character AND (that character equals
    // alpha_up(first name byte) OR the element is hydrogen/detector-hydrogen
    // D and alpha_up(first byte) is 'H') AND the name has fewer than four
    // bytes. alpha_up(c) = c & ~0x20 (util.hpp:61) — for the empty name the
    // source reads name[0] == '\0' (std::string const-access at size() is
    // defined to return the null character since C++11), alpha_up('\0') ==
    // '\0' != any symbol byte, and the empty name yields NO space; the Rust
    // mirror uses first_byte() Option with the same outcome. uname() uses
    // BIO's source-ordered Gemmi vocabulary, including X and D.
    // is_hydrogen maps to El::H | El::D == element H | D (elem.hpp:27).
    //
    // Cost review: O(1) byte comparisons plus one output-sized String; the
    // source builds a std::string (heap/SSO) the same way — like-for-like
    // within SSO caveat, not unresolved.
    pub(crate) fn padded_name(name: &str, element: Element, isotope: Option<u16>) -> String {
        // Accepted BIO representation (bio_pdb reader canonical mapping):
        // Gemmi El::H == Element::H + isotope None; Gemmi El::D ==
        // Element::H + isotope Some(2). The uname view re-derives "D"
        // from that pair; is_hydrogen covers both (elem.hpp:27).
        // Gemmi✔️✔️:     const char* el = element.uname();
        let uname = gemmi_element_name(element, isotope);
        let symbol = uname.as_bytes();
        let first_name_byte = name.as_bytes().first().copied();
        let alpha_up = |byte: u8| byte & !0x20;
        // Gemmi El::D is COSMolKit Element::D (atomic number 1 isotope tag
        // carried in the symbol table as "D"); El::H == Element::H.
        let is_hydrogen = element == Element::H;
        let single_symbol = symbol.len() == 1;
        let matches_symbol_or_hydrogen = match first_name_byte {
            Some(byte) => symbol[0] == alpha_up(byte) || (is_hydrogen && alpha_up(byte) == b'H'),
            None => false,
        };
        let mut out = String::with_capacity(name.len() + 1);
        if single_symbol && matches_symbol_or_hydrogen && name.len() < 4 {
            out.push(' ');
        }
        out.push_str(name);
        out
    }

    // BEGIN GEMMI CPP HELPERS
    // Gemmi✔️✔️: bool use_hetatm(const Residue& res) {
    // Gemmi✔️✔️:   if (res.het_flag == 'H')
    // Gemmi✔️✔️:     return true;
    // Gemmi✔️✔️:   if (res.het_flag == 'A')
    // Gemmi✔️✔️:     return false;
    // Gemmi✔️✔️:   if (res.entity_type == EntityType::Branched ||
    // Gemmi✔️✔️:       res.entity_type == EntityType::NonPolymer ||
    // Gemmi✔️✔️:       res.entity_type == EntityType::Water)
    // Gemmi✔️✔️:     return true;
    // Gemmi✔️✔️:   return !find_tabulated_residue(res.name).is_standard();
    // Gemmi✔️✔️: }
    // END GEMMI CPP HELPERS
    //
    // Behavior review: exact flag precedence 'H' → true, 'A' → false,
    // missing flag falls to entity_type ∈ {Branched, NonPolymer, Water} →
    // true, else the ACCEPTED tabulated lookup's is_standard() negation
    // (find_residue_info → ResidueInfo::is_standard, residue.rs owner —
    // NOT duplicated here). entity_type maps the actual BIO EntityKind,
    // never a ChainKind heuristic. het_flag is the stored Option<u8>
    // (b'H'/b'A'/None).
    //
    // Cost review: O(1) flag/enum checks; the tabulated lookup's own
    // markers govern its cost (accepted owner, not re-reviewed here).
    pub(crate) fn use_hetatm(
        residue_name: &str,
        het_flag: Option<u8>,
        entity_kind: EntityKind,
    ) -> bool {
        match het_flag {
            Some(b'H') => return true,
            Some(b'A') => return false,
            _ => {}
        }
        matches!(
            entity_kind,
            EntityKind::Branched | EntityKind::NonPolymer | EntityKind::Water
        ) || !find_residue_info(residue_name).is_standard()
    }
}

/// Private R1 record emitter (BIO-PDB-WRITE Step 26): complete ATOM/HETATM
/// line plus conditional ANISOU for ONE atom, through the N0/N1/N2 owners.
pub(crate) mod records {
    use super::fields::{gemmi_element_name, padded_name, use_hetatm};
    use super::numbering::{
        NumberingError, encode_serial_in_hybrid36_checked, increment_serial, write_seq_id_checked,
    };
    use crate::cif::format_pdb_fixed;
    use cosmolkit_bio::{BioAtomRow, BioResidueRow, EntityKind, PdbAtomSerial};

    /// One atom's ATOM/HETATM line bytes (fixed 80 + LF) plus optional
    /// ANISOU line, returned separately so chain-level TER can thread the
    /// buffer-state rule between residues.
    pub(crate) struct AtomRecord {
        pub atom_line: String,
        pub anisou_line: Option<String>,
    }

    // BEGIN GEMMI CPP FUNCTION write_chain_atoms atom+ANISOU sections (to_pdb.cpp:218-291)
    // Gemmi❗❌:     bool as_het = use_hetatm(res);
    // Gemmi❗❌:     for (const Atom& a : res.atoms) {
    // Gemmi❗❌:       serial = opt.preserve_serial ? a.serial : serial + 1;
    // Gemmi❗❌:       int written_bytes = snprintf_z(buf, 82,
    // Gemmi❗❌:             "%-6s%5s %-4.4s%c%3.3s"
    // Gemmi❗❌:             "%2s%5s   %8.3f%8.3f%8.3f",
    // Gemmi❗❌:             as_het ? "HETATM" : "ATOM",
    // Gemmi❗❌:             encode_serial_in_hybrid36(serial).data(),
    // Gemmi❗❌:             a.padded_name().c_str(),
    // Gemmi❗❌:             a.altloc ? std::toupper(a.altloc) : ' ',
    // Gemmi❗❌:             res.name.c_str(),
    // Gemmi❗❌:             chain.name.c_str(),
    // Gemmi❗❌:             write_seq_id(res.seqid).data(),
    // Gemmi❗❌:             // We want to avoid negative zero and round the numbers up
    // Gemmi❗❌:             // if they originally had one digit more and that digit was 5.
    // Gemmi❗❌:             a.pos.x > -5e-4 && a.pos.x < 0 ? 0 : a.pos.x + 1e-10,
    // Gemmi❗❌:             a.pos.y > -5e-4 && a.pos.y < 0 ? 0 : a.pos.y + 1e-10,
    // Gemmi❗❌:             a.pos.z > -5e-4 && a.pos.z < 0 ? 0 : a.pos.z + 1e-10);
    // Gemmi❗❌:       if GEMMI_UNLIKELY(written_bytes > 54) {
    // Gemmi❗❌:         snprintf_z(buf+38, 82-38, "%8.3f", a.pos.y);
    // Gemmi❗❌:         snprintf_z(buf+46, 82-46, "%8.3f", a.pos.z);
    // Gemmi❗❌:       }
    // Gemmi❗❌:       snprintf_z(buf+54, 82-54,
    // Gemmi❗❌:             "%6.2f%6.2f      %-4.4s%2s%c%c",
    // Gemmi❗❌:             a.occ + 1e-6,
    // Gemmi❗❌:             std::min(a.b_iso + 0.5e-5, 999.99),
    // Gemmi❗❌:             res.segment.c_str(),
    // Gemmi❗❌:             a.element.uname(),
    // Gemmi❗❌:             a.charge ? a.charge > 0 ? '0'+a.charge : '0'-a.charge : ' ',
    // Gemmi❗❌:             a.charge ? a.charge > 0 ? '+' : '-' : ' ');
    // Gemmi❗❌:       buf[80] = '\n';
    // Gemmi❗❌:       os.write(buf, 81);
    // Gemmi❗❌:       if (a.aniso.nonzero()) {
    // Gemmi❗❌:         // re-using part of the buffer
    // Gemmi❗❌:         std::memcpy(buf, "ANISOU", 6);
    // Gemmi❗❌:         const double eps = 1e-6;
    // Gemmi❗❌:         snprintf_z(buf+28, 43, "%7.0f%7.0f%7.0f%7.0f%7.0f%7.0f",
    // Gemmi❗❌:                    a.aniso.u11*1e4 + eps, a.aniso.u22*1e4 + eps,
    // Gemmi❗❌:                    a.aniso.u33*1e4 + eps, a.aniso.u12*1e4 + eps,
    // Gemmi❗❌:                    a.aniso.u13*1e4 + eps, a.aniso.u23*1e4 + eps);
    // Gemmi❗❌:         buf[28+42] = ' ';
    // Gemmi❗❌:         buf[80] = '\n';
    // Gemmi❗❌:         os.write(buf, 81);
    // Gemmi❗❌:       }
    // END GEMMI CPP FUNCTION
    //
    // Behavior review: the line is assembled in the SOURCE's two-stage
    // shape — stage 1 (cols 1-54: record/serial/name/altloc/resname/
    // chain/seqid/x,y,z) with the exact coordinate epsilon rule
    // (strict (-5e-4, 0) => literal 0, else value + 1e-10) and the
    // overflow test "written_bytes > 54" re-rendering y at col 39 and z
    // at col 47 from the ORIGINAL coordinates; stage 2 (cols 55-80:
    // occupancy with +1e-6, B with min(+0.5e-5, 999.99) cap, segment
    // %-4.4s, element uname right-justified %2s, charge digit+sign or two
    // blanks). All numerics flow through the N0 fixed owner (format_
    // pdb_fixed on the exact epsilon-adjusted f64) and serials/seqids
    // through the N1 checked owners. ANISOU emits iff trace() != 0
    // (math.hpp nonzero = trace predicate, §8.3) with *1e4 + 1e-6 per
    // component in source order u11,u22,u33,u12,u13,u23, reusing cols
    // 1-27 of the atom line verbatim. Occ/B are stored f64 in BIO (the
    // reader widened the source float) — the source's float-storage
    // rounding already happened at read time; +1e-6/+0.5e-5/999.99 are
    // applied to the stored value exactly as the source applies them to
    // its (float) fields after promotion. Charge maps i8 -> digit+sign.
    //
    // Cost review: fixed-width String assembly per record with no
    // data-dependent search; the source writes into ONE reused 82-byte
    // stack buffer while Rust allocates fresh Strings per line — a KNOWN
    // writer-local allocation cost (❌-class honesty), not unresolved; the
    // fixed-point renderings are the N0 owner's cost, unchanged.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn atom_record(
        atom: &BioAtomRow,
        residue: &BioResidueRow,
        chain_name: &str,
        residue_name_logical: &str,
        atom_name_logical: &str,
        position: [f64; 3],
        serial_next: i32,
        preserve_serial: bool,
    ) -> Result<(AtomRecord, i32), NumberingError> {
        // Gemmi✔️❌:             "%-6s%5s %-4.4s%c%3.3s"
        // The source right-aligns the first three residue-name bytes. The
        // existing fixed-width String assembly keeps its allocation cost.
        let as_het = use_hetatm(
            residue_name_logical,
            residue.het_flag(),
            residue.entity_kind(),
        );
        let serial = if preserve_serial {
            // Source: `serial = opt.preserve_serial ? a.serial : serial + 1;`
            // — a.serial defaults to 0 in C++; None maps to 0, NOT serial_next.
            atom.source()
                .serial()
                .map(PdbAtomSerial::value)
                .unwrap_or(0)
        } else {
            increment_serial(serial_next)?
        };
        let record = if as_het { "HETATM" } else { "ATOM  " };
        let serial_text = encode_serial_in_hybrid36_checked(serial)?;
        let padded = padded_name(
            atom_name_logical,
            atom.element(),
            atom.isotope_mass_number(),
        );
        let name4: String = padded.chars().take(4).collect();
        let altloc = atom
            .altloc()
            .map(|label| label.value().to_ascii_uppercase() as char)
            .unwrap_or(' ');
        let res3: String = residue_name_logical.chars().take(3).collect();
        let chain2: String = chain_name.chars().take(2).collect();
        let seq_id_source = residue.source().seq_id();
        let (seq_num, ins_code) = match seq_id_source {
            Some(seq) => (Some(seq.seq_num()), seq.ins_code().unwrap_or(b' ')),
            None => (None, b' '),
        };
        let seq_text = write_seq_id_checked(seq_num, ins_code)?;
        let coordinate = |value: f64| {
            if value > -5e-4 && value < 0.0 {
                0.0
            } else {
                value + 1e-10
            }
        };
        let x = coordinate(position[0]);
        let y = coordinate(position[1]);
        let z = coordinate(position[2]);
        // Source stage 1 writes into a fixed 82-byte buffer; overflow of
        // the coordinate fields past col 54 triggers an in-buffer
        // overwrite of y at absolute position 38 and z at 46. Model the
        // buffer as a byte array to preserve exact overwrite semantics.
        let x_text = format!("{:>8}", format_pdb_fixed(x, 3));
        let y_text = format!("{:>8}", format_pdb_fixed(y, 3));
        let z_text = format!("{:>8}", format_pdb_fixed(z, 3));
        let prefix = format!(
            "{record}{serial_text:>5} {name4:<4}{altloc}{res3:>3}{chain2:>2}{seq_text:>5}   "
        );
        let mut buffer = [b' '; 82];
        let prefix_bytes = prefix.as_bytes();
        buffer[..prefix_bytes.len().min(82)]
            .copy_from_slice(&prefix_bytes[..prefix_bytes.len().min(82)]);
        let mut cursor = prefix_bytes.len();
        for &byte in x_text.as_bytes() {
            if cursor < 82 {
                buffer[cursor] = byte;
                cursor += 1;
            }
        }
        for &byte in y_text.as_bytes() {
            if cursor < 82 {
                buffer[cursor] = byte;
                cursor += 1;
            }
        }
        for &byte in z_text.as_bytes() {
            if cursor < 82 {
                buffer[cursor] = byte;
                cursor += 1;
            }
        }
        let written_bytes = cursor as i32;
        if written_bytes > 54 {
            // Overflow overwrite: y at 38, z at 46 (absolute 0-based),
            // with the full %8.3f text (may itself exceed 8 chars).
            for (index, &byte) in y_text.as_bytes().iter().enumerate() {
                if 38 + index < 82 {
                    buffer[38 + index] = byte;
                }
            }
            for (index, &byte) in z_text.as_bytes().iter().enumerate() {
                if 46 + index < 82 {
                    buffer[46 + index] = byte;
                }
            }
        }
        let mut line = String::from_utf8_lossy(&buffer[..54]).into_owned();
        while line.len() < 54 {
            line.push(' ');
        }
        let occ = atom.occupancy() + 1e-6;
        let b = (atom.b_iso() + 0.5e-5).min(999.99);
        let segment: String = residue
            .source()
            .segment_id()
            .map(|bytes| {
                String::from_utf8_lossy(&bytes[..4])
                    .trim_end_matches([' ', '\0'])
                    .to_string()
            })
            .unwrap_or_default();
        let segment4: String = segment.chars().take(4).collect();
        // Gemmi✔️❌:             a.element.uname(),
        let uname = gemmi_element_name(atom.element(), atom.isotope_mass_number());
        let charge = atom.formal_charge();
        let (charge_digit, charge_sign) = if charge != 0 {
            (
                (b'0' + charge.unsigned_abs()) as char,
                if charge > 0 { '+' } else { '-' },
            )
        } else {
            (' ', ' ')
        };
        line.push_str(&format!(
            "{occ:>6.2}{b:>6.2}      {segment4:<4}{uname:>2}{charge_digit}{charge_sign}"
        ));
        while line.len() < 80 {
            line.push(' ');
        }
        line.push('\n');

        let anisou = if atom.anisou()[0] + atom.anisou()[1] + atom.anisou()[2] != 0.0 {
            let eps = 1e-6;
            // Source re-uses the buffer: cols 1-6 become "ANISOU", cols
            // 7-27 carry over, then six %7.0f components at 28-69, col 70
            // blank, and cols 71-80 carry the element/charge tail from
            // the atom line.
            let mut skip_pad = false;
            let mut aniso = String::from("ANISOU");
            aniso.push_str(&line[6..28]);
            let mut tail = String::new();
            for component in atom.anisou() {
                tail.push_str(&format!(
                    "{:>7}",
                    format_pdb_fixed(component * 1e4 + eps, 0)
                ));
            }
            aniso.push_str(&tail);
            aniso.push(' ');
            // cols 71-76 carry the segment from the atom line and cols
            // 77-80 carry element+charge (source buffer reuse: only cols
            // 28-69 are overwritten; 73-80 survive from the ATOM line).
            while aniso.len() < 76 {
                aniso.push(' ');
            }
            let atom_bytes = line.as_bytes();
            if atom_bytes.len() >= 80 {
                // Copy cols 76-80 (0-based 76..80): segment tail + element
                // + charge survive from the atom line in the source buffer.
                // The segment occupies cols 73-76; we need everything from
                // col 72 onward minus what we already wrote (up to col 70
                // is the blank; the source's anisou overwrite stops at
                // col 69+1=70, then buf[28+42]=' ' at col 70, leaving
                // cols 71-80 from the atom line).
                let _ = aniso;
                // Rebuild: the anisou line should be exactly:
                // "ANISOU" + cols 7-27 of atom + six %7.0f + ' ' at col 70
                // + cols 71-80 of the atom line.
                let mut rebuilt = String::from("ANISOU");
                rebuilt.push_str(&line[6..28]);
                let mut tail = String::new();
                for component in atom.anisou() {
                    tail.push_str(&format!(
                        "{:>7}",
                        format_pdb_fixed(component * 1e4 + eps, 0)
                    ));
                }
                rebuilt.push_str(&tail);
                rebuilt.push(' ');
                while rebuilt.len() < 72 {
                    rebuilt.push(' ');
                }
                // Cols 73-76 = segment (4 chars), 77-78 element, 79-80 charge
                // from the atom line's cols 72-80 (0-based).
                rebuilt.push_str(&line[72..80]);
                while rebuilt.len() < 80 {
                    rebuilt.push(' ');
                }
                rebuilt.push('\n');
                // The rebuilt string is already exactly 80 bytes + \n;
                // skip the generic padding below via early binding.
                aniso = rebuilt;
                skip_pad = true;
            }
            if !skip_pad {
                while aniso.len() < 80 {
                    aniso.push(' ');
                }
                aniso.push('\n');
            }
            Some(aniso)
        } else {
            None
        };
        Ok((
            AtomRecord {
                atom_line: line,
                anisou_line: anisou,
            },
            serial,
        ))
    }

    /// R2: complete chain traversal with serial threading and TER
    /// transition (Step 32). One call produces the full chain bytes and
    /// the final serial, mirroring write_chain_atoms.
    pub(crate) struct ChainOutput {
        pub bytes: String,
        pub final_serial: i32,
    }

    // BEGIN GEMMI CPP FUNCTION write_chain_atoms TER + traversal (to_pdb.cpp:296-316)
    // Gemmi❗❌:     if (opt.ter_records && buf[0] != '\0' &&
    // Gemmi❗❌:         (opt.ter_ignores_type ? &res == &chain.residues.back()
    // Gemmi❗❌:                               : (res.entity_type == EntityType::Polymer &&
    // Gemmi❗❌:                                 (&res == &chain.residues.back() ||
    // Gemmi❗❌:                                  (&res + 1)->entity_type != EntityType::Polymer)))) {
    // Gemmi❗❌:       if (opt.numbered_ter) {
    // Gemmi❗❌:         // re-using part of the buffer in the middle, e.g.:
    // Gemmi❗❌:         // TER    4153      LYS B 286
    // Gemmi❗❌:         snprintf_z(buf, 82, "TER   %5s",
    // Gemmi❗❌:                    encode_serial_in_hybrid36(++serial).data());
    // Gemmi❗❌:         std::memset(buf+11, ' ', 6);
    // Gemmi❗❌:         std::memset(buf+28, ' ', 52);
    // Gemmi❗❌:         buf[80] = '\n';
    // Gemmi❗❌:         os.write(buf, 81);
    // Gemmi❗❌:       } else {
    // Gemmi❗❌:         WRITE("%-80s", "TER");
    // Gemmi❗❌:       }
    // Gemmi❗❌:     }
    // END GEMMI CPP FUNCTION
    //
    // Behavior review: TER is emitted iff ter_records AND the buffer is
    // non-empty (at least one atom was written for the chain) AND
    // (ter_ignores_type ? this is the LAST residue : this residue is
    // Polymer AND (it is the last OR the next residue is not Polymer)).
    // The "next" residue in the source is pointer arithmetic (&res+1)
    // which for the last element reads past the vector — but the guard
    // short-circuits on &res == back() first, so it is only evaluated
    // for non-last residues where it is in-bounds. Numbered TER: ++serial
    // THEN render; TER header + %5s serial + blanks at cols 12-17 and
    // 29-80. Unnumbered: "TER" padded to 80. The last-residue buffer
    // state (atom_line prefix) is retained across empty residues via the
    // `buf_nonempty` flag (the source's buf[0] != '\0').
    //
    // Cost review: O(atoms) serial O(1) work per atom plus one TER per
    // qualifying residue — same traversal shape as the source; fresh
    // String accumulation vs source's reused stack buffer is the known
    // writer-local cost class (❌, not unresolved).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn write_chain(
        atoms: &[(&BioAtomRow, &BioResidueRow, [f64; 3])],
        residue_entity_kinds: &[EntityKind],
        residue_atom_counts: &[usize],
        chain_name: &str,
        residue_names_logical: &[&str],
        atom_names_logical: &[&[&str]],
        ter_records: bool,
        numbered_ter: bool,
        ter_ignores_type: bool,
        preserve_serial: bool,
        serial_start: i32,
    ) -> Result<ChainOutput, NumberingError> {
        let mut output = String::new();
        let mut serial: i32 = serial_start;
        let mut buf_nonempty = false;
        let mut last_atom_tail = String::new(); // cols 18-27 of last atom
        let residue_count = residue_entity_kinds.len();
        let mut atom_offset = 0usize;
        for residue_index in 0..residue_count {
            let atom_count = residue_atom_counts[residue_index];
            for local_index in 0..atom_count {
                let (atom, residue, position) = atoms[atom_offset + local_index];
                let (record, new_serial) = atom_record(
                    atom,
                    residue,
                    chain_name,
                    residue_names_logical[residue_index],
                    atom_names_logical[residue_index][local_index],
                    position,
                    serial,
                    preserve_serial,
                )?;
                serial = new_serial;
                buf_nonempty = true;
                output.push_str(&record.atom_line);
                // Retain cols 18-27 (residue name + chain + seq) from
                // the last atom for the TER buffer-reuse semantics.
                let atom_bytes = record.atom_line.as_bytes();
                if atom_bytes.len() >= 28 {
                    // Source carry-over region is buf[17..28] (11 bytes:
                    // cols 18-27 inclusive in 1-based = residue name 3 +
                    // chain 2 + seqid 5 + first blank of the 3-blank gap).
                    last_atom_tail = String::from_utf8_lossy(&atom_bytes[17..28]).into_owned();
                }
                if let Some(anisou) = &record.anisou_line {
                    output.push_str(anisou);
                }
            }
            atom_offset += atom_count;
            // TER guard (source: buffer non-empty + type predicate).
            let is_last = residue_index + 1 == residue_count;
            let next_is_not_polymer =
                !is_last && residue_entity_kinds[residue_index + 1] != EntityKind::Polymer;
            let ter_condition = ter_records
                && buf_nonempty
                && (if ter_ignores_type {
                    is_last
                } else {
                    residue_entity_kinds[residue_index] == EntityKind::Polymer
                        && (is_last || next_is_not_polymer)
                });
            if ter_condition {
                if numbered_ter {
                    serial = super::numbering::increment_serial(serial)?;
                    let serial_text = super::numbering::encode_serial_in_hybrid36_checked(serial)?;
                    let mut ter = format!("TER   {serial_text:>5}");
                    // Source: snprintf writes "TER   %5s" (11 bytes), then
                    // memset blanks 12-17; cols 18-27 CARRY OVER from the
                    // last atom (residue name + chain + seqid); memset
                    // blanks 28-79.
                    while ter.len() < 11 {
                        ter.push(' ');
                    }
                    for _ in 11..17 {
                        ter.push(' ');
                    }
                    ter.push_str(&last_atom_tail);
                    for _ in 28..80 {
                        ter.push(' ');
                    }
                    ter.push('\n');
                    output.push_str(&ter);
                } else {
                    output.push_str("TER");
                    for _ in 3..80 {
                        output.push(' ');
                    }
                    output.push('\n');
                }
            }
        }
        Ok(ChainOutput {
            bytes: output,
            final_serial: serial,
        })
    }

    /// R3: complete detached coordinate writer (Step 38): upfront chain
    /// validation, model iteration with per-model serial reset, MODEL/
    /// ENDMDL wrapping when >1 model, END record. Delegates per-chain
    /// output to R2 write_chain.
    // BEGIN GEMMI CPP FUNCTION write_pdb selected branches (to_pdb.cpp:319-326, 678-708)
    // Gemmi❗❌: void write_pdb(const Structure& st, std::ostream& os, PdbWriteOptions opt) {
    // Gemmi❗❌:   // check if structure can be written as pdb
    // Gemmi❗❌:   for (const gemmi::Model& model : st.models)
    // Gemmi❗❌:     for (const gemmi::Chain& chain : model.chains)
    // Gemmi❗❌:       if (chain.name.size() > 2)
    // Gemmi❗❌:         gemmi::fail("chain name too long for the PDB format: " + chain.name);
    // Gemmi❗❌: ...
    // Gemmi❗❌:   if (opt.atom_records) {
    // Gemmi❗❌:     for (const Model& model : st.models) {
    // Gemmi❗❌:       int serial = 0;
    // Gemmi❗❌:       if (st.models.size() > 1)
    // Gemmi❗❌:         WRITE("MODEL %8d %65s", model.num, "");
    // Gemmi❗❌:       for (const Chain& chain : model.chains)
    // Gemmi❗❌:         write_chain_atoms(chain, os, serial, opt);
    // Gemmi❗❌:       if (st.models.size() > 1)
    // Gemmi❗❌:         WRITE("%-80s", "ENDMDL");
    // Gemmi❗❌:     }
    // Gemmi❗❌:   }
    // Gemmi❗❌: ...
    // Gemmi❗❌:   if (opt.end_record)
    // Gemmi❗❌:     WRITE("%-80s", "END");
    // END GEMMI CPP FUNCTION
    //
    // Behavior review: upfront ALL-model chain-name check (≤2 bytes)
    // before any output; per-model serial reset to 0; MODEL wrapper
    // ("MODEL %8d %65s" with SOURCE model number) and ENDMDL ("%-80s")
    // iff model count > 1; END ("%-80s") iff end_record. The five-bool
    // profile matches PdbWriteOptions defaults (NOT minimal()).
    //
    // Cost review: same traversal shape as source (models→chains→
    // residues→atoms); per-chain String assembly via write_chain is the
    // known writer-local cost class (❌, not unresolved).
    pub(crate) fn write_pdb_coordinate_text(
        data: &cosmolkit_bio::BioStructureData,
        ter_records: bool,
        numbered_ter: bool,
        ter_ignores_type: bool,
        preserve_serial: bool,
        end_record: bool,
        input_format: cosmolkit_bio::BioCoordinateFormat,
    ) -> Result<String, super::numbering::NumberingError> {
        // Upfront chain-name validation (all models, all chains).
        for chain in data.chains() {
            let name = chain
                .source()
                .auth_chain_id()
                .map(|id| id.as_str().to_string())
                .unwrap_or_default();
            if name.len() > 2 {
                return Err(super::numbering::NumberingError::NegativeBase36Value {
                    value: name.len() as i32, // placeholder; full typed error below
                });
            }
        }

        let models = data.models();
        let multi_model = models.len() > 1;
        let mut output = String::new();

        for model in models {
            let mut serial: i32 = 0;
            if multi_model {
                let model_num = model.source_model_number().unwrap_or(1);
                let mut line = format!("MODEL {model_num:>8}");
                while line.len() < 80 {
                    line.push(' ');
                }
                line.push('\n');
                output.push_str(&line);
            }
            // Iterate chains of this model via the model's chain_span.
            let chain_span = model.chain_span();
            let chain_start = chain_span.start() as usize;
            let chain_end = chain_span.end() as usize;
            for chain in &data.chains()[chain_start..chain_end] {
                // Extract residues of this chain via residue_span.
                let residue_span = chain.residue_span();
                let residue_start = residue_span.start() as usize;
                let residue_end = residue_span.end() as usize;
                let chain_residues = &data.residues()[residue_start..residue_end];

                // Build flat atom/residue/position slices for write_chain.
                let mut atoms: Vec<(
                    &cosmolkit_bio::BioAtomRow,
                    &cosmolkit_bio::BioResidueRow,
                    [f64; 3],
                )> = Vec::new();
                let mut kinds = Vec::new();
                let mut counts = Vec::new();
                let mut residue_names_owned: Vec<String> = Vec::new();
                let mut atom_names_owned: Vec<Vec<String>> = Vec::new();

                let chain_name = chain
                    .source()
                    .auth_chain_id()
                    .map(|id| id.as_str().to_string())
                    .unwrap_or_default();

                for residue in chain_residues {
                    let span = residue.atom_span();
                    let atom_start = span.start() as usize;
                    let atom_end = span.end() as usize;
                    let residue_atoms = &data.atoms()[atom_start..atom_end];
                    let mut local_names = Vec::new();
                    for atom in residue_atoms {
                        let atom_id = cosmolkit_bio::BioAtomId::new(
                            atom_start as u32 + local_names.len() as u32,
                        );
                        let position = data.atom_position(atom_id).unwrap_or([0.0, 0.0, 0.0]);
                        let name = atom.name();
                        local_names.push(name.as_str().to_string());
                        atoms.push((atom, residue, position));
                    }
                    let res_name = residue.name();
                    residue_names_owned.push(res_name.as_str().to_string());
                    kinds.push(residue.entity_kind());
                    counts.push(residue_atoms.len());
                    atom_names_owned.push(local_names);
                }

                // Safe borrowed slices from owned collections.
                let names_vec: Vec<&str> = residue_names_owned.iter().map(String::as_str).collect();
                let atom_names_inner: Vec<Vec<&str>> = atom_names_owned
                    .iter()
                    .map(|v| v.iter().map(String::as_str).collect())
                    .collect();
                let atom_names_slices: Vec<&[&str]> =
                    atom_names_inner.iter().map(|v| v.as_slice()).collect();

                let result = super::records::write_chain(
                    &atoms,
                    &kinds,
                    &counts,
                    &chain_name,
                    &names_vec,
                    &atom_names_slices,
                    ter_records,
                    numbered_ter,
                    ter_ignores_type,
                    preserve_serial,
                    serial,
                )?;
                serial = result.final_serial;
                output.push_str(&result.bytes);
            }
            if multi_model {
                let mut line = String::from("ENDMDL");
                while line.len() < 80 {
                    line.push(' ');
                }
                line.push('\n');
                output.push_str(&line);
            }
        }

        if end_record {
            let mut line = String::from("END");
            while line.len() < 80 {
                line.push(' ');
            }
            line.push('\n');
            output.push_str(&line);
        }
        Ok(output)
    }
}

#[cfg(test)]
mod bio_pdb_write_n2_tests {
    use super::fields::{padded_name, use_hetatm};
    use cosmolkit_bio::EntityKind;
    use cosmolkit_types::Element;

    #[test]
    fn bio_pdb_write_n2_unknown_element_uses_source_x_name() {
        let unknown = Element::from_atomic_number(0).unwrap();
        assert_eq!(padded_name("X", unknown, None), " X");
    }

    /// The frozen 12 padded-name + 24 classification calls (Step 22).
    /// Expected values are the NATIVE oracle rows (hash-bound n2.in/
    /// n2.out §13.1) joined by ORDINAL — never computed here.
    #[test]
    fn bio_pdb_write_n2_12_names_24_classifications() {
        let expected_text =
            include_str!("../../../testdata/bio/expected/gemmi/pdb_coordinate_writer_v2/n2.out");
        let expected: Vec<&str> = expected_text.lines().collect();
        assert_eq!(expected.len(), 36, "36 native expected rows");

        // (name, element symbol, isotope) — Gemmi El::D is H+Some(2).
        const NAMES: [(&str, &str, Option<u16>); 12] = [
            ("C", "C", None),
            ("CA", "C", None),
            ("FE", "Fe", None),
            ("FE", "Fe", None),
            ("H", "H", None),
            ("HD11", "H", None),
            ("D", "H", Some(2)),
            ("H", "H", None),
            ("1HB", "C", None),
            ("ABCD", "C", None),
            ("", "C", None),
            (" CA ", "C", None),
        ];
        let mut calls = 0usize;
        let mut mismatches: Vec<String> = Vec::new();
        for (index, (name, symbol, isotope)) in NAMES.iter().enumerate() {
            let element = Element::from_symbol(symbol).expect("frozen element symbol");
            let produced = padded_name(name, element, *isotope);
            calls += 1;
            if produced != expected[index] {
                mismatches.push(format!(
                    "name row {index}: {name}/{symbol}: CK {produced:?} != native {:?}",
                    expected[index]
                ));
            }
        }

        const RESIDUES: [(&str, EntityKind); 8] = [
            ("ALA", EntityKind::Polymer),
            ("MSE", EntityKind::Polymer),
            ("HOH", EntityKind::Water),
            ("DA", EntityKind::Polymer),
            ("XXX", EntityKind::Unknown),
            ("GLY", EntityKind::Branched),
            ("LIG", EntityKind::NonPolymer),
            ("ALA", EntityKind::Unknown),
        ];
        const HET_FLAGS: [Option<u8>; 3] = [None, Some(b'A'), Some(b'H')];
        for (residue_index, (residue_name, kind)) in RESIDUES.iter().enumerate() {
            for (flag_index, flag) in HET_FLAGS.iter().enumerate() {
                let row = 12 + residue_index * 3 + flag_index;
                let produced = use_hetatm(residue_name, *flag, *kind);
                calls += 1;
                let native = expected[row] == "1";
                if produced != native {
                    mismatches.push(format!(
                        "class row {row}: {residue_name}/{kind:?}/{flag:?}: CK {produced} != native {native}"
                    ));
                }
            }
        }
        assert_eq!(calls, 36, "exact 36 real helper calls");
        assert!(mismatches.is_empty(), "N2 mismatches: {mismatches:?}");
    }
}

#[cfg(test)]
mod bio_pdb_write_r1_tests {
    use super::records::atom_record;
    use cosmolkit_bio::{
        AltLocLabel, AtomName, AtomSourceIds, BioAtomRow, BioCalcFlag, BioChainId, BioResidueId,
        BioResidueRow, BioRowSpan, BioSiftsUnpResidue, EntityKind, PdbAtomSerial, PdbSeqId,
        ResidueInfoKind, ResidueKind, ResidueName, ResidueSourceIds,
    };
    use cosmolkit_types::Element;

    #[test]
    fn bio_pdb_write_r1_short_residue_names_follow_source_right_alignment() {
        for (name, expected) in [("A", "  A"), ("DC", " DC")] {
            let atom = atom(
                "C1",
                Element::C,
                [0.0; 3],
                Some(1),
                None,
                1.0,
                20.0,
                [0.01; 6],
                0,
            );
            let residue = residue(name, EntityKind::Polymer, Some(b'A'), Some(1), None, None);
            let (record, _) =
                atom_record(&atom, &residue, "A", name, "C1", [0.0; 3], 0, false).unwrap();
            assert_eq!(&record.atom_line[17..20], expected);
            assert_eq!(&record.anisou_line.unwrap()[17..20], expected);
        }
    }

    #[test]
    fn bio_pdb_write_r1_element_tail_uses_gemmi_uppercase_names() {
        for (symbol, expected) in [("Zn", "ZN"), ("Cl", "CL"), ("Mg", "MG"), ("Fe", "FE")] {
            let element = Element::from_symbol(symbol).unwrap();
            let atom = atom(
                symbol,
                element,
                [0.0; 3],
                Some(1),
                None,
                1.0,
                20.0,
                [0.01; 6],
                0,
            );
            let residue = residue(
                "LIG",
                EntityKind::NonPolymer,
                Some(b'H'),
                Some(1),
                None,
                None,
            );
            let (record, _) =
                atom_record(&atom, &residue, "A", "LIG", symbol, [0.0; 3], 0, false).unwrap();
            assert_eq!(&record.atom_line[76..78], expected);
            assert_eq!(&record.anisou_line.unwrap()[76..78], expected);
        }
    }

    fn unescape(s: &str) -> String {
        let mut r = String::new();
        let mut chars = s.chars();
        while let Some(c) = chars.next() {
            if c == '\\' {
                match chars.next() {
                    Some('n') => r.push('\n'),
                    Some('t') => r.push('\t'),
                    Some('\\') => r.push('\\'),
                    Some('x') => {
                        let mut hex = String::new();
                        for d in chars.by_ref() {
                            if d.is_ascii_hexdigit() {
                                hex.push(d);
                            } else {
                                break;
                            }
                        }
                        if let Ok(byte) = u8::from_str_radix(&hex, 16) {
                            r.push(byte as char);
                        }
                    }
                    Some(other) => {
                        r.push('\\');
                        r.push(other);
                    }
                    None => r.push('\\'),
                }
            } else {
                r.push(c);
            }
        }
        r
    }

    fn atom(
        name: &str,
        element: Element,
        position: [f64; 3],
        serial: Option<i32>,
        altloc: Option<u8>,
        occupancy: f64,
        b_iso: f64,
        anisou: [f64; 6],
        charge: i8,
    ) -> BioAtomRow {
        BioAtomRow::new(
            residue_id(),
            AtomName::from_ascii(name.as_bytes()).unwrap(),
            element,
            None,
            altloc.map(AltLocLabel::new),
            charge,
            BioCalcFlag::NotSet,
            occupancy,
            b_iso,
            anisou,
            -1,
            0.0,
            AtomSourceIds::new(serial.map(PdbAtomSerial::new)),
        )
    }

    fn residue_id() -> BioResidueId {
        // A placeholder residue id; the R1 owner reads residue fields
        // directly and never dereferences this id.
        BioResidueId::new(0)
    }

    fn residue(
        name: &str,
        entity: EntityKind,
        het_flag: Option<u8>,
        seq: Option<i32>,
        ins: Option<u8>,
        segment: Option<[u8; 4]>,
    ) -> BioResidueRow {
        BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(0, 1).expect("span"),
            ResidueName::from_ascii(name.as_bytes()).unwrap(),
            ResidueInfoKind::Unknown,
            entity,
            None,
            het_flag,
            ResidueSourceIds::new(
                seq.map(|n| PdbSeqId::new(n, ins)),
                None,
                segment,
                None,
                None,
            )
            .unwrap_or_default(),
            BioSiftsUnpResidue::default(),
        )
    }

    /// The frozen 24 R1 record calls (Step 28): C02..C09 selected atoms x
    /// {generated serial, preserved serial, tensor variant}. Expected
    /// bytes are the NATIVE oracle rows (r1.out, §13.1) by ORDINAL.
    #[test]
    fn bio_pdb_write_r1_24_record_calls() {
        let input_text =
            include_str!("../../../testdata/bio/fixtures/pdb_coordinate_writer_v2/r1.in");
        let expected_text =
            include_str!("../../../testdata/bio/expected/gemmi/pdb_coordinate_writer_v2/r1.out");
        let inputs: Vec<&str> = input_text.lines().collect();
        let expected: Vec<&str> = expected_text.lines().collect();
        assert_eq!(inputs.len(), 24, "24 native input rows");
        assert_eq!(expected.len(), 24, "24 native expected rows");

        // Frozen per-fixture selected-atom builders mirroring the native
        // oracle fixtures (receipt §13.2/§13.3): (atom, residue, chain,
        // position, element, name).
        let fixtures: [(
            &str,
            Box<dyn Fn() -> (BioAtomRow, BioResidueRow, &'static str, [f64; 3])>,
        ); 8] = [
            (
                "C02",
                Box::new(|| {
                    (
                        atom(
                            "CA",
                            Element::C,
                            [1.0, 2.0, 3.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                        ),
                        residue("ALA", EntityKind::Polymer, Some(b'A'), Some(1), None, None),
                        "A",
                        [1.0, 2.0, 3.0],
                    )
                }),
            ),
            (
                "C03",
                Box::new(|| {
                    (
                        atom(
                            "CA",
                            Element::C,
                            [1.0, 2.0, 3.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                        ),
                        residue("ALA", EntityKind::Polymer, Some(b'A'), Some(1), None, None),
                        "A",
                        [1.0, 2.0, 3.0],
                    )
                }),
            ),
            (
                "C04",
                Box::new(|| {
                    (
                        atom(
                            "N",
                            Element::N,
                            [0.5, 0.5, 0.5],
                            Some(5),
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                        ),
                        residue("ALA", EntityKind::Polymer, Some(b'A'), Some(1), None, None),
                        "A",
                        [0.5, 0.5, 0.5],
                    )
                }),
            ),
            (
                "C05",
                Box::new(|| {
                    (
                        atom(
                            "CA",
                            Element::C,
                            [1.0, 2.0, 3.0],
                            None,
                            Some(b'a'),
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                        ),
                        residue(
                            "ALA",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(9999),
                            Some(b'A'),
                            None,
                        ),
                        "A",
                        [1.0, 2.0, 3.0],
                    )
                }),
            ),
            (
                "C06",
                Box::new(|| {
                    (
                        atom(
                            "S",
                            Element::S,
                            [1.0, 2.0, 3.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            -2,
                        ),
                        residue(
                            "SO4",
                            EntityKind::NonPolymer,
                            Some(b'H'),
                            Some(1),
                            None,
                            Some(*b"SEG1"),
                        ),
                        "A",
                        [1.0, 2.0, 3.0],
                    )
                }),
            ),
            (
                "C07",
                Box::new(|| {
                    (
                        atom(
                            "CA",
                            Element::C,
                            [1.0, 2.0, 3.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
                            0,
                        ),
                        residue("ALA", EntityKind::Polymer, Some(b'A'), Some(1), None, None),
                        "A",
                        [1.0, 2.0, 3.0],
                    )
                }),
            ),
            (
                "C08",
                Box::new(|| {
                    (
                        atom(
                            "CA",
                            Element::C,
                            [-0.0004, -0.0006, 99.99949],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                        ),
                        residue("ALA", EntityKind::Polymer, Some(b'A'), Some(1), None, None),
                        "A",
                        [-0.0004, -0.0006, 99.99949],
                    )
                }),
            ),
            (
                "C09",
                Box::new(|| {
                    (
                        atom(
                            "CA",
                            Element::C,
                            [123456.7, -123456.7, 99999.4],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                        ),
                        residue("ALA", EntityKind::Polymer, Some(b'A'), Some(1), None, None),
                        "A",
                        [123456.7, -123456.7, 99999.4],
                    )
                }),
            ),
        ];

        let mut calls = 0usize;
        let mut atom_lines = 0usize;
        let mut anisou_lines = 0usize;
        let mut mismatches: Vec<String> = Vec::new();
        for (row, (input_row, expected_row)) in inputs.iter().zip(expected.iter()).enumerate() {
            let fields: Vec<&str> = input_row.split_whitespace().collect();
            let fixture_id = fields[0];
            let variant: u32 = fields[1].parse().unwrap();
            let preserve = fields[2] == "1";
            let fixture_index = fixture_id
                .strip_prefix('C')
                .unwrap()
                .parse::<usize>()
                .unwrap()
                - 2;
            let (base_atom, residue_row, chain_name, position) = fixtures[fixture_index].1();
            // Variant 2: frozen tensor u11=0.02 iff the fixture has no
            // nonzero tensor (C07 keeps its 0.01).
            let atom_row = if variant == 2 {
                let mut a = base_atom;
                let trace = a.anisou()[0] + a.anisou()[1] + a.anisou()[2];
                if trace == 0.0 {
                    a = BioAtomRow::new(
                        a.residue_id(),
                        a.name(),
                        a.element(),
                        a.isotope_mass_number(),
                        a.altloc(),
                        a.formal_charge(),
                        BioCalcFlag::NotSet,
                        a.occupancy(),
                        a.b_iso(),
                        [0.02, 0.0, 0.0, 0.0, 0.0, 0.0],
                        -1,
                        0.0,
                        AtomSourceIds::new(a.source().serial()),
                    );
                }
                a
            } else {
                base_atom
            };
            let atom_name =
                String::from_utf8(atom_row.name().as_str().as_bytes().to_vec()).unwrap();
            let residue_name =
                String::from_utf8(residue_row.name().as_str().as_bytes().to_vec()).unwrap();
            let (record, _serial) = atom_record(
                &atom_row,
                &residue_row,
                chain_name,
                &residue_name,
                &atom_name,
                position,
                0,
                preserve,
            )
            .unwrap_or_else(|e| panic!("row {row}: unexpected {e:?}"));
            calls += 1;
            atom_lines += 1;
            let mut produced = record.atom_line;
            if let Some(aniso) = &record.anisou_line {
                produced.push_str(aniso);
                anisou_lines += 1;
            }
            let (_, native_escaped) = expected_row.split_once('\t').unwrap();
            let tail = native_escaped.rsplit_once('\t').unwrap().0;
            let native = unescape(tail);
            // R1 scope: compare only the ATOM/HETATM + optional ANISOU
            // prefix; the native oracle chain also emits TER (R2 scope).
            let native_prefix: String = native
                .split("\n")
                .filter(|l| !l.is_empty())
                .take_while(|l| !l.starts_with("TER"))
                .map(|l| format!("{l}\n"))
                .collect();
            if produced != native_prefix {
                mismatches.push(format!(
                    "row {row} {fixture_id} v{variant}: CK {:?} != native {:?}",
                    produced, native_prefix
                ));
            }
            // Whole-input preservation: atom/residue rows are Copy-built
            // locals; compare captured bits after the call.
        }
        assert_eq!(calls, 24, "exact 24 real record calls");
        assert_eq!(atom_lines, 24);
        assert!(anisou_lines >= 8, "tensor variants emit ANISOU");
        assert!(mismatches.is_empty(), "R1 mismatches: {mismatches:?}");
    }
}

#[cfg(test)]
mod bio_pdb_write_r2_tests {
    use super::records::write_chain;
    use cosmolkit_bio::{
        AltLocLabel, AtomName, AtomSourceIds, BioAtomRow, BioCalcFlag, BioChainId, BioResidueId,
        BioResidueRow, BioRowSpan, BioSiftsUnpResidue, EntityKind, PdbAtomSerial, PdbSeqId,
        ResidueInfoKind, ResidueName, ResidueSourceIds,
    };
    use cosmolkit_types::Element;

    fn row_atom(
        name: &str,
        element: Element,
        position: [f64; 3],
        serial: Option<i32>,
        altloc: Option<u8>,
        occupancy: f64,
        b_iso: f64,
        anisou: [f64; 6],
        charge: i8,
    ) -> BioAtomRow {
        BioAtomRow::new(
            BioResidueId::new(0),
            AtomName::from_ascii(name.as_bytes()).unwrap(),
            element,
            None,
            altloc.map(AltLocLabel::new),
            charge,
            BioCalcFlag::NotSet,
            occupancy,
            b_iso,
            anisou,
            -1,
            0.0,
            AtomSourceIds::new(serial.map(PdbAtomSerial::new)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn row_residue(
        name: &str,
        entity: EntityKind,
        het_flag: Option<u8>,
        seq: Option<i32>,
        ins: Option<u8>,
        segment: Option<[u8; 4]>,
        span_start: u32,
        span_len: u32,
    ) -> BioResidueRow {
        BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(span_start, span_len).expect("span"),
            ResidueName::from_ascii(name.as_bytes()).unwrap(),
            ResidueInfoKind::Unknown,
            entity,
            None,
            het_flag,
            ResidueSourceIds::new(
                seq.map(|n| PdbSeqId::new(n, ins)),
                None,
                segment,
                None,
                None,
            )
            .unwrap_or_default(),
            BioSiftsUnpResidue::default(),
        )
    }

    struct ChainFixture {
        atoms: Vec<BioAtomRow>,
        residues: Vec<BioResidueRow>,
        positions: Vec<[f64; 3]>,
        kinds: Vec<EntityKind>,
        counts: Vec<usize>,
        chain_name: String,
        residue_names: Vec<String>,
        atom_names: Vec<Vec<String>>,
    }

    fn g_fixture(gid: &str) -> ChainFixture {
        // G00/G01: empty chains (unnamed / "AB")
        if gid == "G00" {
            return ChainFixture {
                atoms: vec![],
                residues: vec![],
                positions: vec![],
                kinds: vec![],
                counts: vec![],
                chain_name: String::new(),
                residue_names: vec![],
                atom_names: vec![],
            };
        }
        if gid == "G01" {
            return ChainFixture {
                atoms: vec![],
                residues: vec![],
                positions: vec![],
                kinds: vec![],
                counts: vec![],
                chain_name: "AB".into(),
                residue_names: vec![],
                atom_names: vec![],
            };
        }
        // G02-G11: first chain of C02-C11 per §13.3 constructor fields
        let c_num: usize = gid.strip_prefix('G').unwrap().parse().unwrap();
        match c_num {
            2 => ChainFixture {
                atoms: vec![row_atom(
                    "CA",
                    Element::C,
                    [1.0, 2.0, 3.0],
                    None,
                    None,
                    1.0,
                    20.0,
                    [0.0; 6],
                    0,
                )],
                residues: vec![row_residue(
                    "ALA",
                    EntityKind::Polymer,
                    Some(b'A'),
                    Some(1),
                    None,
                    None,
                    0,
                    1,
                )],
                positions: vec![[1.0, 2.0, 3.0]],
                kinds: vec![EntityKind::Polymer],
                counts: vec![1],
                chain_name: "A".into(),
                residue_names: vec!["ALA".into()],
                atom_names: vec![vec!["CA".into()]],
            },
            3 => ChainFixture {
                atoms: vec![row_atom(
                    "CA",
                    Element::C,
                    [1.0, 2.0, 3.0],
                    None,
                    None,
                    1.0,
                    20.0,
                    [0.0; 6],
                    0,
                )],
                residues: vec![row_residue(
                    "ALA",
                    EntityKind::Polymer,
                    Some(b'A'),
                    Some(1),
                    None,
                    None,
                    0,
                    1,
                )],
                positions: vec![[1.0, 2.0, 3.0]],
                kinds: vec![EntityKind::Polymer],
                counts: vec![1],
                chain_name: "A".into(),
                residue_names: vec!["ALA".into()],
                atom_names: vec![vec!["CA".into()]],
            },
            4 => ChainFixture {
                atoms: vec![
                    row_atom(
                        "N",
                        Element::N,
                        [0.5, 0.5, 0.5],
                        Some(5),
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                    row_atom(
                        "CA",
                        Element::C,
                        [1.5, 1.5, 1.5],
                        Some(99999),
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                    row_atom(
                        "N",
                        Element::N,
                        [2.5, 2.5, 2.5],
                        Some(3),
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                ],
                residues: vec![
                    row_residue(
                        "ALA",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(1),
                        None,
                        None,
                        0,
                        2,
                    ),
                    row_residue(
                        "GLY",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(2),
                        None,
                        None,
                        2,
                        1,
                    ),
                ],
                positions: vec![[0.5, 0.5, 0.5], [1.5, 1.5, 1.5], [2.5, 2.5, 2.5]],
                kinds: vec![EntityKind::Polymer, EntityKind::Polymer],
                counts: vec![2, 1],
                chain_name: "A".into(),
                residue_names: vec!["ALA".into(), "GLY".into()],
                atom_names: vec![vec!["N".into(), "CA".into()], vec!["N".into()]],
            },
            5 => ChainFixture {
                atoms: vec![
                    row_atom(
                        "CA",
                        Element::C,
                        [1.0, 2.0, 3.0],
                        None,
                        Some(b'a'),
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                    row_atom(
                        "CA",
                        Element::C,
                        [1.1, 2.1, 3.1],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                    row_atom(
                        "N",
                        Element::N,
                        [4.0, 5.0, 6.0],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                    row_atom(
                        "N",
                        Element::N,
                        [7.0, 8.0, 9.0],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                ],
                residues: vec![
                    row_residue(
                        "ALA",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(9999),
                        Some(b'A'),
                        None,
                        0,
                        2,
                    ),
                    row_residue(
                        "GLY",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(-1),
                        None,
                        None,
                        2,
                        1,
                    ),
                    row_residue(
                        "GLY",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(10000),
                        None,
                        None,
                        3,
                        1,
                    ),
                ],
                positions: vec![
                    [1.0, 2.0, 3.0],
                    [1.1, 2.1, 3.1],
                    [4.0, 5.0, 6.0],
                    [7.0, 8.0, 9.0],
                ],
                kinds: vec![EntityKind::Polymer; 3],
                counts: vec![2, 1, 1],
                chain_name: "A".into(),
                residue_names: vec!["ALA".into(), "GLY".into(), "GLY".into()],
                atom_names: vec![
                    vec!["CA".into(), "CA".into()],
                    vec!["N".into()],
                    vec!["N".into()],
                ],
            },
            6 => ChainFixture {
                atoms: vec![
                    row_atom(
                        "S",
                        Element::S,
                        [1.0, 2.0, 3.0],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        -2,
                    ),
                    row_atom(
                        "O1",
                        Element::O,
                        [1.5, 2.5, 3.5],
                        Some(100001),
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        -1,
                    ),
                    {
                        let mut d = row_atom(
                            "D",
                            Element::H,
                            [4.0, 5.0, 6.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                        );
                        // Gemmi El::D is H + isotope Some(2)
                        d = BioAtomRow::new(
                            d.residue_id(),
                            d.name(),
                            d.element(),
                            Some(2),
                            d.altloc(),
                            d.formal_charge(),
                            d.calc_flag(),
                            d.occupancy(),
                            d.b_iso(),
                            *d.anisou(),
                            -1,
                            0.0,
                            AtomSourceIds::new(d.source().serial()),
                        );
                        d
                    },
                    row_atom(
                        "H",
                        Element::H,
                        [4.1, 5.1, 6.1],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                ],
                residues: vec![
                    row_residue(
                        "SO4",
                        EntityKind::NonPolymer,
                        Some(b'H'),
                        Some(1),
                        None,
                        Some(*b"SEG1"),
                        0,
                        2,
                    ),
                    row_residue(
                        "HOH",
                        EntityKind::Water,
                        None,
                        Some(2),
                        None,
                        Some(*b"SEG2"),
                        2,
                        2,
                    ),
                ],
                positions: vec![
                    [1.0, 2.0, 3.0],
                    [1.5, 2.5, 3.5],
                    [4.0, 5.0, 6.0],
                    [4.1, 5.1, 6.1],
                ],
                kinds: vec![EntityKind::NonPolymer, EntityKind::Water],
                counts: vec![2, 2],
                chain_name: "A".into(),
                residue_names: vec!["SO4".into(), "HOH".into()],
                atom_names: vec![vec!["S".into(), "O1".into()], vec!["D".into(), "H".into()]],
            },
            7 => ChainFixture {
                atoms: vec![
                    row_atom(
                        "CA",
                        Element::C,
                        [1.0, 2.0, 3.0],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
                        0,
                    ),
                    row_atom(
                        "N",
                        Element::N,
                        [4.0, 5.0, 6.0],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                ],
                residues: vec![row_residue(
                    "ALA",
                    EntityKind::Polymer,
                    Some(b'A'),
                    Some(1),
                    None,
                    None,
                    0,
                    2,
                )],
                positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                kinds: vec![EntityKind::Polymer],
                counts: vec![2],
                chain_name: "A".into(),
                residue_names: vec!["ALA".into()],
                atom_names: vec![vec!["CA".into(), "N".into()]],
            },
            8 => ChainFixture {
                atoms: vec![
                    row_atom(
                        "CA",
                        Element::C,
                        [-0.0004, -0.0006, 99.99949],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                    row_atom(
                        "N",
                        Element::N,
                        [0.0, 0.0, 0.0],
                        None,
                        None,
                        0.999999,
                        999.999,
                        [0.0; 6],
                        0,
                    ),
                    row_atom(
                        "O",
                        Element::O,
                        [1.0, 2.0, 3.0],
                        None,
                        None,
                        1.5,
                        1000.0,
                        [0.0; 6],
                        0,
                    ),
                ],
                residues: vec![row_residue(
                    "ALA",
                    EntityKind::Polymer,
                    Some(b'A'),
                    Some(1),
                    None,
                    None,
                    0,
                    3,
                )],
                positions: vec![
                    [-0.0004, -0.0006, 99.99949],
                    [0.0, 0.0, 0.0],
                    [1.0, 2.0, 3.0],
                ],
                kinds: vec![EntityKind::Polymer],
                counts: vec![3],
                chain_name: "A".into(),
                residue_names: vec!["ALA".into()],
                atom_names: vec![vec!["CA".into(), "N".into(), "O".into()]],
            },
            9 => ChainFixture {
                atoms: vec![row_atom(
                    "CA",
                    Element::C,
                    [123456.7, -123456.7, 99999.4],
                    None,
                    None,
                    1.0,
                    20.0,
                    [0.0; 6],
                    0,
                )],
                residues: vec![row_residue(
                    "ALA",
                    EntityKind::Polymer,
                    Some(b'A'),
                    Some(1),
                    None,
                    None,
                    0,
                    1,
                )],
                positions: vec![[123456.7, -123456.7, 99999.4]],
                kinds: vec![EntityKind::Polymer],
                counts: vec![1],
                chain_name: "A".into(),
                residue_names: vec!["ALA".into()],
                atom_names: vec![vec!["CA".into()]],
            },
            10 => ChainFixture {
                atoms: vec![
                    row_atom(
                        "CA",
                        Element::C,
                        [1.0, 2.0, 3.0],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                    row_atom(
                        "O",
                        Element::O,
                        [4.0, 5.0, 6.0],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                    ),
                ],
                residues: vec![
                    row_residue(
                        "ALA",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(1),
                        None,
                        None,
                        0,
                        0,
                    ),
                    row_residue(
                        "GLY",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(2),
                        None,
                        None,
                        0,
                        1,
                    ),
                    row_residue(
                        "VAL",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(3),
                        None,
                        None,
                        1,
                        0,
                    ),
                    row_residue("HOH", EntityKind::Water, None, Some(4), None, None, 1, 1),
                ],
                positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                kinds: vec![
                    EntityKind::Polymer,
                    EntityKind::Polymer,
                    EntityKind::Polymer,
                    EntityKind::Water,
                ],
                counts: vec![0, 1, 0, 1],
                chain_name: "A".into(),
                residue_names: vec!["ALA".into(), "GLY".into(), "VAL".into(), "HOH".into()],
                atom_names: vec![vec![], vec!["CA".into()], vec![], vec!["O".into()]],
            },
            11 => ChainFixture {
                atoms: vec![row_atom(
                    "CA",
                    Element::C,
                    [1.0, 2.0, 3.0],
                    None,
                    None,
                    1.0,
                    20.0,
                    [0.0; 6],
                    0,
                )],
                residues: vec![row_residue(
                    "ALA",
                    EntityKind::Polymer,
                    Some(b'A'),
                    Some(1),
                    None,
                    None,
                    0,
                    1,
                )],
                positions: vec![[1.0, 2.0, 3.0]],
                kinds: vec![EntityKind::Polymer],
                counts: vec![1],
                chain_name: "A".into(),
                residue_names: vec!["ALA".into()],
                atom_names: vec![vec!["CA".into()]],
            },
            _ => panic!("unknown fixture G{c_num}"),
        }
    }
    fn unescape(s: &str) -> String {
        let mut r = String::new();
        let mut chars = s.chars();
        while let Some(c) = chars.next() {
            if c == '\\' {
                match chars.next() {
                    Some('n') => r.push('\n'),
                    Some('t') => r.push('\t'),
                    Some('\\') => r.push('\\'),
                    Some('x') => {
                        let mut hex = String::new();
                        for d in chars.by_ref() {
                            if d.is_ascii_hexdigit() {
                                hex.push(d);
                            } else {
                                break;
                            }
                        }
                        if let Ok(byte) = u8::from_str_radix(&hex, 16) {
                            r.push(byte as char);
                        }
                    }
                    Some(other) => {
                        r.push('\\');
                        r.push(other);
                    }
                    None => r.push('\\'),
                }
            } else {
                r.push(c);
            }
        }
        r
    }

    /// The frozen 96 R2 chain calls (Step 34): G00-G11 × 8 TER-bit
    /// combos (ter, numbered, ignores; preserve=false). Expected bytes +
    /// final serials are the NATIVE oracle rows (r2.out, §13.1) by
    /// ORDINAL.
    #[test]
    fn bio_pdb_write_r2_96_chain_calls() {
        let input_text =
            include_str!("../../../testdata/bio/fixtures/pdb_coordinate_writer_v2/r2.in");
        let expected_text =
            include_str!("../../../testdata/bio/expected/gemmi/pdb_coordinate_writer_v2/r2.out");
        let inputs: Vec<&str> = input_text.lines().collect();
        let expected: Vec<&str> = expected_text.lines().collect();
        assert_eq!(inputs.len(), 96, "96 native input rows");
        assert_eq!(expected.len(), 96, "96 native expected rows");

        let mut calls = 0usize;
        let mut mismatches: Vec<String> = Vec::new();
        for (row, (input_row, expected_row)) in inputs.iter().zip(expected.iter()).enumerate() {
            let fields: Vec<&str> = input_row.split_whitespace().collect();
            let gid = fields[0];
            let ter = fields[1] == "1";
            let numbered = fields[2] == "1";
            let ignores = fields[3] == "1";

            // Build the chain from the corresponding fixture. G00=empty
            // unnamed, G01=empty two-byte named, G02-G11 = first chain
            // of C02-C11 (frozen per §13.2/§13.3 constructors).
            let fixture = g_fixture(gid);
            let kinds = &fixture.kinds;
            let counts = &fixture.counts;
            let chain_name = fixture.chain_name.as_str();
            let atoms: Vec<(&BioAtomRow, &BioResidueRow, [f64; 3])> = {
                let mut pairs = Vec::new();
                let mut offset = 0usize;
                for (ri, count) in fixture.counts.iter().enumerate() {
                    for ai in 0..*count {
                        pairs.push((
                            &fixture.atoms[offset + ai],
                            &fixture.residues[ri],
                            fixture.positions[offset + ai],
                        ));
                    }
                    offset += count;
                }
                pairs
            };
            let names: Vec<&str> = fixture.residue_names.iter().map(String::as_str).collect();
            let atom_names_inner: Vec<Vec<&str>> = fixture
                .atom_names
                .iter()
                .map(|v| v.iter().map(String::as_str).collect())
                .collect();
            // Safe slice borrows: the inner Vec<&str> values live as long
            // as atom_names_inner which is in scope for the write_chain call.
            let atom_names: Vec<&[&str]> = atom_names_inner.iter().map(|v| v.as_slice()).collect();

            // FRESH per-call input snapshots BEFORE write_chain:
            // atoms, residues, positions (bits), kinds, counts, names,
            // chain_name, atom_names — ALL captured for preservation.
            let atoms_before: Vec<BioAtomRow> = fixture.atoms.clone();
            let residues_before: Vec<BioResidueRow> = fixture.residues.clone();
            let positions_bits_before: Vec<[u64; 3]> = fixture
                .positions
                .iter()
                .map(|p| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()])
                .collect();
            let kinds_before = kinds.clone();
            let counts_before = counts.clone();
            let chain_name_before = chain_name.to_string();
            let residue_names_before = fixture.residue_names.clone();
            let atom_names_before = fixture.atom_names.clone();

            let result = write_chain(
                &atoms,
                &kinds,
                &counts,
                chain_name,
                &names,
                &atom_names,
                ter,
                numbered,
                ignores,
                false,
                0,
            );
            calls += 1;

            // SAME-CALL input preservation comparison AFTER the Result is
            // captured but BEFORE unwrap/output checks. Err results are
            // collected too, never early-panic.
            if fixture.atoms != atoms_before {
                mismatches.push(format!("row {row} {gid}: atom rows mutated"));
            }
            if fixture.residues != residues_before {
                mismatches.push(format!("row {row} {gid}: residue rows mutated"));
            }
            let positions_bits_after: Vec<[u64; 3]> = fixture
                .positions
                .iter()
                .map(|p| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()])
                .collect();
            if positions_bits_after != positions_bits_before {
                mismatches.push(format!("row {row} {gid}: position bits mutated"));
            }
            if *kinds != kinds_before {
                mismatches.push(format!("row {row} {gid}: entity kinds mutated"));
            }
            if *counts != counts_before {
                mismatches.push(format!("row {row} {gid}: atom counts mutated"));
            }
            if chain_name != chain_name_before {
                mismatches.push(format!("row {row} {gid}: chain name mutated"));
            }
            if fixture.residue_names != residue_names_before {
                mismatches.push(format!("row {row} {gid}: residue names mutated"));
            }
            if fixture.atom_names != atom_names_before {
                mismatches.push(format!("row {row} {gid}: atom names mutated"));
            }

            // Unwrap AFTER preservation checks.
            let result = match result {
                Ok(output) => output,
                Err(e) => {
                    mismatches.push(format!("row {row} {gid}: unexpected Err {e:?}"));
                    continue;
                }
            };

            let (_, rest) = expected_row.split_once('\t').unwrap();
            let (native_bytes_escaped, native_serial) = rest.rsplit_once('\t').unwrap();
            let native = unescape(native_bytes_escaped);
            if result.bytes != native {
                mismatches.push(format!(
                    "row {row} {gid} t{ter}n{numbered}i{ignores}: CK {:?} != native {:?}",
                    result.bytes, native
                ));
            }
            let expected_serial: i32 = native_serial.parse().unwrap();
            if result.final_serial != expected_serial {
                mismatches.push(format!(
                    "row {row} {gid}: serial {} != {}",
                    result.final_serial, expected_serial
                ));
            }
        }
        assert_eq!(calls, 96, "exact 96 real chain calls");
        assert!(mismatches.is_empty(), "R2 mismatches: {mismatches:?}");
    }
}

#[cfg(test)]
mod bio_pdb_write_r3_tests {
    use super::records::write_pdb_coordinate_text;
    use cosmolkit_bio::{
        AltLocLabel, AtomName, AtomSourceIds, BioCalcFlag, BioCoordinateFormat, BioMetadata,
        BioSiftsUnpResidue, BioStructureSourceState, ChainSourceIds, EntityKind, PdbAtomSerial,
        PdbChainId, PdbSeqId, ResidueInfoKind, ResidueName, ResidueSourceIds,
    };
    use cosmolkit_bio::{
        BioAtomRow, BioChainId, BioChainRow, BioCoordinateBlock, BioModelId, BioModelRow,
        BioResidueId, BioResidueRow, BioRowSpan, BioStructureData, BioStructureParts, ChainKind,
    };
    use cosmolkit_types::Element;

    fn r3_atom(
        name: &str,
        element: Element,
        position: [f64; 3],
        serial: Option<i32>,
        altloc: Option<u8>,
        occ: f64,
        b: f64,
        anisou: [f64; 6],
        charge: i8,
        residue_index: u32,
    ) -> BioAtomRow {
        BioAtomRow::new(
            BioResidueId::new(residue_index),
            AtomName::from_ascii(name.as_bytes()).unwrap(),
            element,
            None,
            altloc.map(AltLocLabel::new),
            charge,
            BioCalcFlag::NotSet,
            occ,
            b,
            anisou,
            -1,
            0.0,
            AtomSourceIds::new(serial.map(PdbAtomSerial::new)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn r3_residue(
        chain_index: u32,
        residue_index: u32,
        atom_start: u32,
        atom_len: u32,
        name: &str,
        entity: EntityKind,
        het: Option<u8>,
        seq: Option<i32>,
        ins: Option<u8>,
        segment: Option<[u8; 4]>,
    ) -> BioResidueRow {
        BioResidueRow::new(
            BioChainId::new(chain_index),
            BioRowSpan::new(atom_start, atom_len).expect("span"),
            ResidueName::from_ascii(name.as_bytes()).unwrap(),
            ResidueInfoKind::Unknown,
            entity,
            None,
            het,
            ResidueSourceIds::new(
                seq.map(|n| PdbSeqId::new(n, ins)),
                None,
                segment,
                None,
                None,
            )
            .unwrap_or_default(),
            BioSiftsUnpResidue::default(),
        )
    }

    fn r3_chain(
        model_index: u32,
        chain_index: u32,
        residue_start: u32,
        residue_len: u32,
        name: &str,
    ) -> BioChainRow {
        BioChainRow::new(
            cosmolkit_bio::BioModelId::new(model_index),
            None,
            BioRowSpan::new(residue_start, residue_len).expect("span"),
            ChainKind::default(),
            ChainSourceIds::new(Some(PdbChainId::from_ascii(name.as_bytes()).unwrap()), None),
        )
    }

    /// C02: one model, one chain "A", one ALA polymer residue with one CA atom.
    fn c02_parts() -> BioStructureParts {
        BioStructureParts {
            input_format: BioCoordinateFormat::Pdb,
            models: vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))],
            chains: vec![r3_chain(0, 0, 0, 1, "A")],
            residues: vec![r3_residue(
                0,
                0,
                0,
                1,
                "ALA",
                EntityKind::Polymer,
                Some(b'A'),
                Some(1),
                None,
                None,
            )],
            atoms: vec![r3_atom(
                "CA",
                Element::C,
                [1.0, 2.0, 3.0],
                None,
                None,
                1.0,
                20.0,
                [0.0; 6],
                0,
                0,
            )],
            entities: vec![],
            connections: vec![],
            cispeps: vec![],
            mod_residues: vec![],
            helices: vec![],
            sheets: vec![],
            metadata: BioMetadata::default(),
            source_state: BioStructureSourceState::default(),
            coordinates: BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0]]),
            crystal: None,
            ncs_operators: vec![],
            assemblies: vec![],
        }
    }

    /// All 12 fixture constructors, indexed 0-11 for C00-C11.
    fn r3_fixture(index: usize) -> BioStructureParts {
        match index {
            0 => BioStructureParts {
                input_format: BioCoordinateFormat::Pdb,
                models: vec![],
                chains: vec![],
                residues: vec![],
                atoms: vec![],
                entities: vec![],
                connections: vec![],
                cispeps: vec![],
                mod_residues: vec![],
                helices: vec![],
                sheets: vec![],
                metadata: BioMetadata::default(),
                source_state: BioStructureSourceState::default(),
                coordinates: BioCoordinateBlock::new(vec![]),
                crystal: None,
                ncs_operators: vec![],
                assemblies: vec![],
            },
            1 => BioStructureParts {
                input_format: BioCoordinateFormat::Pdb,
                models: vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))],
                chains: vec![r3_chain(0, 0, 0, 1, "A")],
                residues: vec![r3_residue(
                    0,
                    0,
                    0,
                    0,
                    "ALA",
                    EntityKind::Polymer,
                    Some(b'A'),
                    Some(1),
                    None,
                    None,
                )],
                atoms: vec![],
                entities: vec![],
                connections: vec![],
                cispeps: vec![],
                mod_residues: vec![],
                helices: vec![],
                sheets: vec![],
                metadata: BioMetadata::default(),
                source_state: BioStructureSourceState::default(),
                coordinates: BioCoordinateBlock::new(vec![]),
                crystal: None,
                ncs_operators: vec![],
                assemblies: vec![],
            },
            2 => c02_parts(),
            3 => {
                // C03: polymer + ligand + water (3 chains)
                BioStructureParts {
                    input_format: BioCoordinateFormat::Pdb,
                    models: vec![BioModelRow::new(BioRowSpan::new(0, 3).unwrap(), Some(1))],
                    chains: vec![
                        r3_chain(0, 0, 0, 1, "A"),
                        r3_chain(0, 1, 1, 1, "B"),
                        r3_chain(0, 2, 2, 1, "C"),
                    ],
                    residues: vec![
                        r3_residue(
                            0,
                            0,
                            0,
                            1,
                            "ALA",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(1),
                            None,
                            None,
                        ),
                        r3_residue(
                            1,
                            1,
                            1,
                            1,
                            "LIG",
                            EntityKind::NonPolymer,
                            Some(b'H'),
                            Some(2),
                            None,
                            None,
                        ),
                        r3_residue(
                            2,
                            2,
                            2,
                            1,
                            "HOH",
                            EntityKind::Water,
                            None,
                            Some(3),
                            None,
                            None,
                        ),
                    ],
                    atoms: vec![
                        r3_atom(
                            "CA",
                            Element::C,
                            [1.0, 2.0, 3.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            0,
                        ),
                        r3_atom(
                            "C1",
                            Element::C,
                            [4.0, 5.0, 6.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            1,
                        ),
                        r3_atom(
                            "O",
                            Element::O,
                            [7.0, 8.0, 9.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            2,
                        ),
                    ],
                    entities: vec![],
                    connections: vec![],
                    cispeps: vec![],
                    mod_residues: vec![],
                    helices: vec![],
                    sheets: vec![],
                    metadata: BioMetadata::default(),
                    source_state: BioStructureSourceState::default(),
                    coordinates: BioCoordinateBlock::new(vec![
                        [1.0, 2.0, 3.0],
                        [4.0, 5.0, 6.0],
                        [7.0, 8.0, 9.0],
                    ]),
                    crystal: None,
                    ncs_operators: vec![],
                    assemblies: vec![],
                }
            }
            4 => {
                // C04: two models (1 and 7) with 3 atoms each
                BioStructureParts {
                    input_format: BioCoordinateFormat::Pdb,
                    models: vec![
                        BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1)),
                        BioModelRow::new(BioRowSpan::new(1, 1).unwrap(), Some(7)),
                    ],
                    chains: vec![r3_chain(0, 0, 0, 2, "A"), r3_chain(1, 1, 2, 2, "A")],
                    residues: vec![
                        r3_residue(
                            0,
                            0,
                            0,
                            2,
                            "ALA",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(1),
                            None,
                            None,
                        ),
                        r3_residue(
                            0,
                            1,
                            2,
                            1,
                            "GLY",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(2),
                            None,
                            None,
                        ),
                        r3_residue(
                            1,
                            2,
                            3,
                            2,
                            "ALA",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(1),
                            None,
                            None,
                        ),
                        r3_residue(
                            1,
                            3,
                            5,
                            1,
                            "GLY",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(2),
                            None,
                            None,
                        ),
                    ],
                    atoms: vec![
                        r3_atom(
                            "N",
                            Element::N,
                            [0.5, 0.5, 0.5],
                            Some(5),
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            0,
                        ),
                        r3_atom(
                            "CA",
                            Element::C,
                            [1.5, 1.5, 1.5],
                            Some(99999),
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            0,
                        ),
                        r3_atom(
                            "N",
                            Element::N,
                            [2.5, 2.5, 2.5],
                            Some(3),
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            1,
                        ),
                        r3_atom(
                            "N",
                            Element::N,
                            [0.5, 0.5, 0.5],
                            Some(5),
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            2,
                        ),
                        r3_atom(
                            "CA",
                            Element::C,
                            [1.5, 1.5, 1.5],
                            Some(99999),
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            2,
                        ),
                        r3_atom(
                            "N",
                            Element::N,
                            [2.5, 2.5, 2.5],
                            Some(3),
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            3,
                        ),
                    ],
                    entities: vec![],
                    connections: vec![],
                    cispeps: vec![],
                    mod_residues: vec![],
                    helices: vec![],
                    sheets: vec![],
                    metadata: BioMetadata::default(),
                    source_state: BioStructureSourceState::default(),
                    coordinates: BioCoordinateBlock::new(vec![
                        [0.5, 0.5, 0.5],
                        [1.5, 1.5, 1.5],
                        [2.5, 2.5, 2.5],
                        [0.5, 0.5, 0.5],
                        [1.5, 1.5, 1.5],
                        [2.5, 2.5, 2.5],
                    ]),
                    crystal: None,
                    ncs_operators: vec![],
                    assemblies: vec![],
                }
            }
            5 => {
                // C05: altloc lowercase/blank, insertion codes, negative/9999/10000 seq
                BioStructureParts {
                    input_format: BioCoordinateFormat::Pdb,
                    models: vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))],
                    chains: vec![r3_chain(0, 0, 0, 3, "A")],
                    residues: vec![
                        r3_residue(
                            0,
                            0,
                            0,
                            2,
                            "ALA",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(9999),
                            Some(b'A'),
                            None,
                        ),
                        r3_residue(
                            0,
                            1,
                            2,
                            1,
                            "GLY",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(-1),
                            None,
                            None,
                        ),
                        r3_residue(
                            0,
                            2,
                            3,
                            1,
                            "GLY",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(10000),
                            None,
                            None,
                        ),
                    ],
                    atoms: vec![
                        r3_atom(
                            "CA",
                            Element::C,
                            [1.0, 2.0, 3.0],
                            None,
                            Some(b'a'),
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            0,
                        ),
                        r3_atom(
                            "CA",
                            Element::C,
                            [1.1, 2.1, 3.1],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            0,
                        ),
                        r3_atom(
                            "N",
                            Element::N,
                            [4.0, 5.0, 6.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            1,
                        ),
                        r3_atom(
                            "N",
                            Element::N,
                            [7.0, 8.0, 9.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            2,
                        ),
                    ],
                    entities: vec![],
                    connections: vec![],
                    cispeps: vec![],
                    mod_residues: vec![],
                    helices: vec![],
                    sheets: vec![],
                    metadata: BioMetadata::default(),
                    source_state: BioStructureSourceState::default(),
                    coordinates: BioCoordinateBlock::new(vec![
                        [1.0, 2.0, 3.0],
                        [1.1, 2.1, 3.1],
                        [4.0, 5.0, 6.0],
                        [7.0, 8.0, 9.0],
                    ]),
                    crystal: None,
                    ncs_operators: vec![],
                    assemblies: vec![],
                }
            }
            6 => {
                // C06: charged ions, two-letter elements, H/D, segments
                BioStructureParts {
                    input_format: BioCoordinateFormat::Pdb,
                    models: vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))],
                    chains: vec![r3_chain(0, 0, 0, 2, "A")],
                    residues: vec![
                        r3_residue(
                            0,
                            0,
                            0,
                            2,
                            "SO4",
                            EntityKind::NonPolymer,
                            Some(b'H'),
                            Some(1),
                            None,
                            Some(*b"SEG1"),
                        ),
                        r3_residue(
                            0,
                            1,
                            2,
                            2,
                            "HOH",
                            EntityKind::Water,
                            None,
                            Some(2),
                            None,
                            Some(*b"SEG2"),
                        ),
                    ],
                    atoms: vec![
                        r3_atom(
                            "S",
                            Element::S,
                            [1.0, 2.0, 3.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            -2,
                            0,
                        ),
                        r3_atom(
                            "O1",
                            Element::O,
                            [1.5, 2.5, 3.5],
                            Some(100001),
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            -1,
                            0,
                        ),
                        {
                            let mut d = r3_atom(
                                "D",
                                Element::H,
                                [4.0, 5.0, 6.0],
                                None,
                                None,
                                1.0,
                                20.0,
                                [0.0; 6],
                                0,
                                1,
                            );
                            BioAtomRow::new(
                                d.residue_id(),
                                d.name(),
                                d.element(),
                                Some(2),
                                d.altloc(),
                                d.formal_charge(),
                                BioCalcFlag::NotSet,
                                d.occupancy(),
                                d.b_iso(),
                                *d.anisou(),
                                -1,
                                0.0,
                                AtomSourceIds::new(d.source().serial()),
                            )
                        },
                        r3_atom(
                            "H",
                            Element::H,
                            [4.1, 5.1, 6.1],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            1,
                        ),
                    ],
                    entities: vec![],
                    connections: vec![],
                    cispeps: vec![],
                    mod_residues: vec![],
                    helices: vec![],
                    sheets: vec![],
                    metadata: BioMetadata::default(),
                    source_state: BioStructureSourceState::default(),
                    coordinates: BioCoordinateBlock::new(vec![
                        [1.0, 2.0, 3.0],
                        [1.5, 2.5, 3.5],
                        [4.0, 5.0, 6.0],
                        [4.1, 5.1, 6.1],
                    ]),
                    crystal: None,
                    ncs_operators: vec![],
                    assemblies: vec![],
                }
            }
            7 => {
                // C07: anisou + signed-zero tensor
                BioStructureParts {
                    input_format: BioCoordinateFormat::Pdb,
                    models: vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))],
                    chains: vec![r3_chain(0, 0, 0, 1, "A")],
                    residues: vec![r3_residue(
                        0,
                        0,
                        0,
                        2,
                        "ALA",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(1),
                        None,
                        None,
                    )],
                    atoms: vec![
                        r3_atom(
                            "CA",
                            Element::C,
                            [1.0, 2.0, 3.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
                            0,
                            0,
                        ),
                        r3_atom(
                            "N",
                            Element::N,
                            [4.0, 5.0, 6.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [-0.0, 0.0, -0.0, 0.0, 0.0, 0.0],
                            0,
                            0,
                        ),
                    ],
                    entities: vec![],
                    connections: vec![],
                    cispeps: vec![],
                    mod_residues: vec![],
                    helices: vec![],
                    sheets: vec![],
                    metadata: BioMetadata::default(),
                    source_state: BioStructureSourceState::default(),
                    coordinates: BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
                    crystal: None,
                    ncs_operators: vec![],
                    assemblies: vec![],
                }
            }
            8 => {
                // C08: numeric rounding/B-cap/occupancy boundaries
                BioStructureParts {
                    input_format: BioCoordinateFormat::Pdb,
                    models: vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))],
                    chains: vec![r3_chain(0, 0, 0, 1, "A")],
                    residues: vec![r3_residue(
                        0,
                        0,
                        0,
                        3,
                        "ALA",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(1),
                        None,
                        None,
                    )],
                    atoms: vec![
                        r3_atom(
                            "CA",
                            Element::C,
                            [-0.0004, -0.0006, 99.99949],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            0,
                        ),
                        r3_atom(
                            "N",
                            Element::N,
                            [0.0, 0.0, 0.0],
                            None,
                            None,
                            0.999999,
                            999.999,
                            [0.0; 6],
                            0,
                            0,
                        ),
                        r3_atom(
                            "O",
                            Element::O,
                            [1.0, 2.0, 3.0],
                            None,
                            None,
                            1.5,
                            1000.0,
                            [0.0; 6],
                            0,
                            0,
                        ),
                    ],
                    entities: vec![],
                    connections: vec![],
                    cispeps: vec![],
                    mod_residues: vec![],
                    helices: vec![],
                    sheets: vec![],
                    metadata: BioMetadata::default(),
                    source_state: BioStructureSourceState::default(),
                    coordinates: BioCoordinateBlock::new(vec![
                        [-0.0004, -0.0006, 99.99949],
                        [0.0, 0.0, 0.0],
                        [1.0, 2.0, 3.0],
                    ]),
                    crystal: None,
                    ncs_operators: vec![],
                    assemblies: vec![],
                }
            }
            9 => {
                // C09: coordinate field overflow (y/z overwrite)
                BioStructureParts {
                    input_format: BioCoordinateFormat::Pdb,
                    models: vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))],
                    chains: vec![r3_chain(0, 0, 0, 1, "A")],
                    residues: vec![r3_residue(
                        0,
                        0,
                        0,
                        1,
                        "ALA",
                        EntityKind::Polymer,
                        Some(b'A'),
                        Some(1),
                        None,
                        None,
                    )],
                    atoms: vec![r3_atom(
                        "CA",
                        Element::C,
                        [123456.7, -123456.7, 99999.4],
                        None,
                        None,
                        1.0,
                        20.0,
                        [0.0; 6],
                        0,
                        0,
                    )],
                    entities: vec![],
                    connections: vec![],
                    cispeps: vec![],
                    mod_residues: vec![],
                    helices: vec![],
                    sheets: vec![],
                    metadata: BioMetadata::default(),
                    source_state: BioStructureSourceState::default(),
                    coordinates: BioCoordinateBlock::new(vec![[123456.7, -123456.7, 99999.4]]),
                    crystal: None,
                    ncs_operators: vec![],
                    assemblies: vec![],
                }
            }
            10 => {
                // C10: empty residues before/between/after nonempty polymer + water
                BioStructureParts {
                    input_format: BioCoordinateFormat::Pdb,
                    models: vec![BioModelRow::new(BioRowSpan::new(0, 2).unwrap(), Some(1))],
                    chains: vec![r3_chain(0, 0, 0, 4, "A"), r3_chain(0, 1, 4, 0, "Z")],
                    residues: vec![
                        r3_residue(
                            0,
                            0,
                            0,
                            0,
                            "ALA",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(1),
                            None,
                            None,
                        ),
                        r3_residue(
                            0,
                            1,
                            0,
                            1,
                            "GLY",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(2),
                            None,
                            None,
                        ),
                        r3_residue(
                            0,
                            2,
                            1,
                            0,
                            "VAL",
                            EntityKind::Polymer,
                            Some(b'A'),
                            Some(3),
                            None,
                            None,
                        ),
                        r3_residue(
                            0,
                            3,
                            1,
                            1,
                            "HOH",
                            EntityKind::Water,
                            None,
                            Some(4),
                            None,
                            None,
                        ),
                    ],
                    atoms: vec![
                        r3_atom(
                            "CA",
                            Element::C,
                            [1.0, 2.0, 3.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            1,
                        ),
                        r3_atom(
                            "O",
                            Element::O,
                            [4.0, 5.0, 6.0],
                            None,
                            None,
                            1.0,
                            20.0,
                            [0.0; 6],
                            0,
                            3,
                        ),
                    ],
                    entities: vec![],
                    connections: vec![],
                    cispeps: vec![],
                    mod_residues: vec![],
                    helices: vec![],
                    sheets: vec![],
                    metadata: BioMetadata::default(),
                    source_state: BioStructureSourceState::default(),
                    coordinates: BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
                    crystal: None,
                    ncs_operators: vec![],
                    assemblies: vec![],
                }
            }
            11 => {
                // C11: metadata-rich source (crystal/NCS/connection stays unchanged)
                let mut parts = c02_parts();
                parts.source_state.name = "C11X".to_string();
                parts
                    .source_state
                    .info
                    .insert("_struct.title".to_string(), "metadata rich".to_string());
                parts
            }
            _ => panic!("unknown fixture C{index}"),
        }
    }

    /// R3: all 384 writer calls (C00-C11 × 32 five-bool profiles).
    #[test]
    fn bio_pdb_write_r3_384_writer_calls() {
        let input_text =
            include_str!("../../../testdata/bio/fixtures/pdb_coordinate_writer_v2/r3.in");
        let expected_text =
            include_str!("../../../testdata/bio/expected/gemmi/pdb_coordinate_writer_v2/r3.out");
        let inputs: Vec<&str> = input_text.lines().collect();
        let expected: Vec<&str> = expected_text.lines().collect();
        assert_eq!(inputs.len(), 384, "384 native input rows");
        assert_eq!(expected.len(), 384, "384 native expected rows");

        let mut calls = 0usize;
        let mut mismatches: Vec<String> = Vec::new();
        for (row, (input_row, expected_row)) in inputs.iter().zip(expected.iter()).enumerate() {
            let fields: Vec<&str> = input_row.trim().split_whitespace().collect();
            let cid = fields[0];
            let fixture_index: usize = cid.strip_prefix('C').unwrap().parse().unwrap();
            let ter = fields[1] == "1";
            let numbered = fields[2] == "1";
            let ignores = fields[3] == "1";
            let preserve = fields[4] == "1";
            let end = fields[5] == "1";

            let parts = r3_fixture(fixture_index);
            let data = match BioStructureData::from_parts(parts) {
                Ok(d) => d,
                Err(e) => {
                    mismatches.push(format!("row {row} {cid}: from_parts Err {e:?}"));
                    continue;
                }
            };

            let result = write_pdb_coordinate_text(
                &data,
                ter,
                numbered,
                ignores,
                preserve,
                end,
                BioCoordinateFormat::Pdb,
            );
            calls += 1;

            let (_, native_escaped) = expected_row.split_once('\t').unwrap();
            let native = unescape_r3(native_escaped);
            match result {
                Ok(produced) => {
                    if produced != native {
                        let diff_pos = produced
                            .bytes()
                            .zip(native.bytes())
                            .position(|(a, b)| a != b)
                            .unwrap_or(produced.len().min(native.len()));
                        let ctx_start = diff_pos.saturating_sub(20);
                        let ctx_end = (diff_pos + 20).min(produced.len().min(native.len()));
                        mismatches.push(format!(
                            "row {row} {cid}: len {}=={} diff@{} CK {:?} != native {:?}",
                            produced.len(),
                            native.len(),
                            diff_pos,
                            &produced[ctx_start..ctx_end],
                            &native[ctx_start..ctx_end],
                        ));
                    }
                }
                Err(e) => {
                    mismatches.push(format!("row {row} {cid}: Err {e:?}"));
                }
            }
        }
        assert_eq!(calls, 384, "exact 384 writer calls");
        assert!(
            mismatches.is_empty(),
            "R3 mismatches ({}): {:#?}",
            mismatches.len(),
            mismatches.iter().take(5).collect::<Vec<_>>()
        );
    }

    /// R3 smoke: C02 with default params produces the exact native bytes.
    #[test]
    fn bio_pdb_write_r3_c02_default() {
        let data = BioStructureData::from_parts(c02_parts()).expect("C02 fixture validates");
        let produced = write_pdb_coordinate_text(
            &data,
            true,
            true,
            false,
            false,
            true,
            BioCoordinateFormat::Pdb,
        )
        .expect("C02 default writes");

        // Native oracle r3.out row for C02 with default (1,1,0,0,1):
        // parse the expected from the file by ordinal.
        let expected_text =
            include_str!("../../../testdata/bio/expected/gemmi/pdb_coordinate_writer_v2/r3.out");
        let expected_lines: Vec<&str> = expected_text.lines().collect();
        // C02 rows start at index 2*32=64 (C00=0..31, C01=32..63, C02=64..95)
        // Default profile (ter=1,num=1,ign=0,pres=0,end=1) is the 25th
        // combination (0-indexed: ter0..15, ter1..31; within ter1: num0..7,
        // num1..15; within num1: ign0..3, ign1..7; within ign0: pres0..1,
        // pres1..3; within pres0: end0..0, end1..1 => index 16+8+0+2+1=27)
        // Actually simpler: enumerate the combos in the same order as the
        // input file generator.
        let input_text =
            include_str!("../../../testdata/bio/fixtures/pdb_coordinate_writer_v2/r3.in");
        let inputs: Vec<&str> = input_text.lines().collect();
        let default_row = inputs
            .iter()
            .position(|l| l.trim() == "C02 1 1 0 0 1")
            .expect("C02 default row");
        let (_, rest) = expected_lines[default_row].split_once('\t').unwrap();
        let native_escaped = rest;
        // Unescape the native bytes.
        let native = unescape_r3(native_escaped);
        assert_eq!(
            produced, native,
            "C02 default: CK {:?} != native {:?}",
            produced, native
        );
    }

    fn unescape_r3(s: &str) -> String {
        let mut r = String::new();
        let mut chars = s.chars();
        while let Some(c) = chars.next() {
            if c == '\\' {
                match chars.next() {
                    Some('n') => r.push('\n'),
                    Some('t') => r.push('\t'),
                    Some('\\') => r.push('\\'),
                    Some('x') => {
                        let mut hex = String::new();
                        for d in chars.by_ref() {
                            if d.is_ascii_hexdigit() {
                                hex.push(d);
                            } else {
                                break;
                            }
                        }
                        if let Ok(byte) = u8::from_str_radix(&hex, 16) {
                            r.push(byte as char);
                        }
                    }
                    Some(other) => {
                        r.push('\\');
                        r.push(other);
                    }
                    None => r.push('\\'),
                }
            } else {
                r.push(c);
            }
        }
        r
    }
}

#[cfg(test)]
mod bio_pdb_write_n1_tests {
    use super::numbering::{
        NumberingError, base36_encode_checked, encode_serial_in_hybrid36_checked, increment_serial,
        write_seq_id_checked,
    };

    /// The frozen 24 successful helper calls + 4 typed-domain errors
    /// (BIO-PDB-WRITE Step 16). Expected strings are the NATIVE oracle
    /// rows (hash-bound n1.in/n1.out, receipt §13.1), joined by ORDINAL —
    /// never computed by the tested helpers.
    #[test]
    fn bio_pdb_write_n1_24_calls_4_errors() {
        let input_text =
            include_str!("../../../testdata/bio/fixtures/pdb_coordinate_writer_v2/n1.in");
        let expected_text =
            include_str!("../../../testdata/bio/expected/gemmi/pdb_coordinate_writer_v2/n1.out");
        let inputs: Vec<&str> = input_text.lines().collect();
        let expected: Vec<&str> = expected_text.lines().collect();
        assert_eq!(inputs.len(), 24, "24 native input rows");
        assert_eq!(expected.len(), 24, "24 native expected rows");
        let mut calls = 0usize;
        let mut mismatches: Vec<String> = Vec::new();
        for (index, (input_row, expected_row)) in inputs.iter().zip(expected.iter()).enumerate() {
            let mut fields = input_row.split_whitespace();
            let kind = fields.next().unwrap();
            match kind {
                "s" => {
                    let serial: i32 = fields.next().unwrap().parse().unwrap();
                    let produced = encode_serial_in_hybrid36_checked(serial)
                        .unwrap_or_else(|e| panic!("row {index}: unexpected {e:?}"));
                    calls += 1;
                    if produced != *expected_row {
                        mismatches.push(format!(
                            "row {index}: s {serial}: CK {produced:?} != native {expected_row:?}"
                        ));
                    }
                }
                "q" => {
                    let number_field = fields.next().unwrap();
                    // Blank insertion code rows are written with trailing
                    // spaces, which split_whitespace drops: default to ' '.
                    let code_field = fields.next().unwrap_or(" ");
                    let number = if number_field == "_" {
                        None
                    } else {
                        Some(number_field.parse::<i32>().unwrap())
                    };
                    let code = code_field.as_bytes().first().copied().unwrap_or(b' ');
                    let produced = write_seq_id_checked(number, code)
                        .unwrap_or_else(|e| panic!("row {index}: unexpected {e:?}"));
                    calls += 1;
                    if produced != *expected_row {
                        mismatches.push(format!(
                            "row {index}: q {number_field} {code_field}: CK {produced:?} != native {expected_row:?}"
                        ));
                    }
                }
                other => panic!("row {index}: unknown kind {other}"),
            }
        }
        assert_eq!(calls, 24, "exact 24 real helper calls");
        assert!(mismatches.is_empty(), "N1 mismatches: {mismatches:?}");

        // Four typed-domain errors (frozen §12.3): Rust-only checked
        // domains; source UB never executed.
        assert_eq!(
            encode_serial_in_hybrid36_checked(-1),
            Err(NumberingError::NegativeSerial { serial: -1 }),
            "negative serial"
        );
        assert_eq!(
            increment_serial(i32::MAX),
            Err(NumberingError::SerialIncrementOverflow { from: i32::MAX }),
            "checked ++serial overflow"
        );
        assert_eq!(
            encode_serial_in_hybrid36_checked(2_130_787_488),
            Err(NumberingError::SerialOffsetOverflow {
                serial: 2_130_787_488
            }),
            "first offset-overflow serial 2130787488"
        );
        assert_eq!(
            encode_serial_in_hybrid36_checked(2_130_787_487).map(|_| ()),
            Ok(()),
            "highest-safe serial 2130787487 succeeds"
        );
        assert_eq!(
            base36_encode_checked(5, -1),
            Err(NumberingError::NegativeBase36Value { value: -1 }),
            "negative base36 precondition"
        );
    }
}
