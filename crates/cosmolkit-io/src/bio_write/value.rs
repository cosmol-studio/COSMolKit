//! Numeric text profiles selected by the pinned mmCIF atom-site writer.
//!
//! Coordinates are doubles and use Gemmi `to_str(double)` (`%.9g`);
//! occupancy, B-factor, deuterium fraction and anisotropic tensors are floats
//! in the source model and use Gemmi `to_str(float)` (`%.6g`). The underlying
//! conversions are the existing `format_cif_f64`/`format_cif_f32` machinery;
//! this module only binds the writer's field selection to those owners.

use crate::cif::{format_cif_f32, format_cif_f64, quote_cif_value};
use cosmolkit_bio::{AtomName, BioCoordinateFormat, PdbSeqId, ResidueName, ResidueSourceIds};

/// Gemmi `to_str(double)` for atom-site coordinate columns (`Cartn_x/y/z`).
pub(crate) fn coordinate_text(value: f64) -> String {
    // Gemmi✔️✔️: inline std::string to_str(double d) {
    // Gemmi✔️✔️:   char buf[24];
    // Gemmi✔️✔️:   int len = sprintf_z(buf, "%.9g", d);
    // Gemmi✔️✔️:   return std::string(buf, len > 0 ? len : 0);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: vv.emplace_back(to_str(atom.pos.x));
    // Gemmi✔️✔️: vv.emplace_back(to_str(atom.pos.y));
    // Gemmi✔️✔️: vv.emplace_back(to_str(atom.pos.z));
    // Behavior: source Position is double; nine significant digits.
    // Complexity: delegated to the existing format_cif_f64 owner.
    format_cif_f64(value)
}

/// Gemmi `to_str(float)` for occupancy, B-factor, fraction and anisou cells.
pub(crate) fn float_field_text(value: f64) -> String {
    // Gemmi✔️✔️: inline std::string to_str(float d) {
    // Gemmi✔️✔️:   char buf[16];
    // Gemmi✔️✔️:   int len = sprintf_z(buf, "%.6g", d);
    // Gemmi✔️✔️:   return std::string(buf, len > 0 ? len : 0);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: vv.emplace_back(to_str(atom.occ));
    // Gemmi✔️✔️: vv.emplace_back(to_str(atom.b_iso));
    // Gemmi✔️✔️: vv.emplace_back(to_str(atom.fraction));
    // Gemmi✔️✔️: aniso_val.emplace_back(to_str(atom->aniso.u11));
    // Behavior: the source stores these fields as float, so the BIO f64
    // storage is narrowed through f32 exactly once before the six-significant
    // -digit conversion; this preserves source precision semantics.
    // Complexity: delegated to the existing format_cif_f32 owner.
    format_cif_f32(value as f32)
}

/// Gemmi `number_or_dot`: NaN or absent becomes `.`.
pub(crate) fn number_or_dot(value: Option<f64>) -> String {
    // Gemmi✔️✔️: inline std::string number_or_dot(double d) {
    // Gemmi✔️✔️:   return std::isnan(d) ? "." : to_str(d);
    // Gemmi✔️✔️: }
    // Behavior: the BIO Option encodes the source NaN-absence; a present
    // NaN also maps to "." exactly as the source double does.
    // Complexity: one branch plus the v01 double profile.
    value
        .filter(|value| !value.is_nan())
        .map_or_else(|| ".".to_string(), coordinate_text)
}

/// Gemmi `number_or_qmark`: NaN or absent becomes `?`.
pub(crate) fn number_or_qmark(value: Option<f64>) -> String {
    // Gemmi✔️✔️: inline std::string number_or_qmark(double d) {
    // Gemmi✔️✔️:   return std::isnan(d) ? "?" : to_str(d);
    // Gemmi✔️✔️: }
    value
        .filter(|value| !value.is_nan())
        .map_or_else(|| "?".to_string(), coordinate_text)
}

/// Gemmi `int_or_dot`: the source `-1` N/A sentinel (or absence) becomes `.`.
pub(crate) fn int_or_dot(value: Option<i32>) -> String {
    // Gemmi✔️✔️: // for use with non-negative Metadata fields that use -1 for N/A
    // Gemmi✔️✔️: inline std::string int_or_dot(int n) {
    // Gemmi✔️✔️:   return n == -1 ? "." : std::to_string(n);
    // Gemmi✔️✔️: }
    // Behavior: BIO represents "no value" as None; a stored -1 is the
    // source's own N/A sentinel and maps to "." exactly like the source.
    if matches!(value, None | Some(-1)) {
        ".".to_string()
    } else {
        value.expect("checked above").to_string()
    }
}

/// Gemmi `int_or_qmark`: the source `-1` N/A sentinel (or absence) becomes `?`.
pub(crate) fn int_or_qmark(value: Option<i32>) -> String {
    // Gemmi✔️✔️: inline std::string int_or_qmark(int n) {
    // Gemmi✔️✔️:   return n == -1 ? "?" : std::to_string(n);
    // Gemmi✔️✔️: }
    if matches!(value, None | Some(-1)) {
        "?".to_string()
    } else {
        value.expect("checked above").to_string()
    }
}

/// Gemmi `string_or_dot`: an empty string becomes `.`; others are CIF-quoted.
pub(crate) fn string_or_dot(value: &str) -> String {
    // Gemmi✔️✔️: inline std::string string_or_dot(const std::string& s) {
    // Gemmi✔️✔️:   return s.empty() ? "." : cif::quote(s);
    // Gemmi✔️✔️: }
    if value.is_empty() {
        ".".to_string()
    } else {
        quote_cif_value(value.to_string())
    }
}

/// Gemmi `string_or_qmark`: an empty string becomes `?`; others are CIF-quoted.
pub(crate) fn string_or_qmark(value: &str) -> String {
    // Gemmi✔️✔️: inline std::string string_or_qmark(const std::string& s) {
    // Gemmi✔️✔️:   return s.empty() ? "?" : cif::quote(s);
    // Gemmi✔️✔️: }
    if value.is_empty() {
        "?".to_string()
    } else {
        quote_cif_value(value.to_string())
    }
}

/// Gemmi `qchain`: quote the author chain id; unlike string_or_* an empty
/// name becomes `''`, never `.` or `?`.
pub(crate) fn qchain(name: &str) -> String {
    // Gemmi✔️✔️: // Quote chain name or entity id if necessary. It is necessary
    // Gemmi✔️✔️: // only if the chain name is missing, which was OK in the past.
    // Gemmi✔️✔️: // Here we use '' rather than . or ?.
    // Gemmi✔️✔️: inline std::string qchain(const std::string& s) {
    // Gemmi✔️✔️:   return cif::quote(s);
    // Gemmi✔️✔️: }
    // Behavior: plain cif::quote; the empty string becomes '' because
    // quote() always wraps empty input in quotes.
    // Complexity: delegated to quote_cif_value.
    quote_cif_value(name.to_string())
}

/// Gemmi `subchain_or_dot`: the label subchain id or `.` when absent/empty.
pub(crate) fn subchain_or_dot(residue: &ResidueSourceIds) -> String {
    // Gemmi✔️✔️: inline std::string subchain_or_dot(const Residue& res) {
    // Gemmi✔️✔️:   return res.subchain.empty() ? "." : cif::quote(res.subchain);
    // Gemmi✔️✔️: }
    // Behavior: BIO Option<String> encodes absence; a present empty string
    // still maps to "." exactly like the source.
    // Complexity: one branch plus quoting.
    match residue.subchain_id() {
        Some(subchain) if !subchain.is_empty() => quote_cif_value(subchain.to_string()),
        _ => ".".to_string(),
    }
}

/// Gemmi `pdbx_icode`: the insertion code character or `?` when unset.
pub(crate) fn pdbx_icode(seq_id: Option<PdbSeqId>) -> String {
    // Gemmi✔️✔️: inline std::string pdbx_icode(const SeqId& seqid) {
    // Gemmi✔️✔️:   return std::string(1, seqid.has_icode() ? seqid.icode : '?');
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: char has_icode() const { return icode != ' '; }  // seqid.hpp:79
    // Behavior: Gemmi treats a space icode as unset; BIO encodes unset as
    // None, and a stored b' ' byte also maps to '?' exactly like has_icode.
    // Complexity: O(1).
    match seq_id.and_then(PdbSeqId::ins_code) {
        Some(code) if code != b' ' => String::from(code as char),
        _ => "?".to_string(),
    }
}

/// Atom-name text for the writer, keyed on authoritative document provenance.
///
/// `input_format` decides: only `Pdb` documents carry raw four-column
/// names, which Gemmi's reader trims via `read_string`; every CIF-family
/// format carries decoded logical names that must be preserved verbatim
/// (quoted mmCIF `' CA '` keeps its spaces).
pub(crate) fn atom_name_text(name: AtomName, input_format: BioCoordinateFormat) -> String {
    // BIO-ROWS R01: delegated to the canonical BIO owner
    // (`cosmolkit_bio::atom_name_logical_view`, which carries the pinned
    // `read_string` anchors); this wrapper only materializes the String
    // the writer loop appends to.
    cosmolkit_bio::atom_name_logical_view(&name, input_format).to_string()
}

/// Residue-name text for the writer, keyed on the same provenance rule.
pub(crate) fn residue_name_text(name: &ResidueName, input_format: BioCoordinateFormat) -> &str {
    // BIO-ROWS R02: delegated to the canonical BIO owner
    // (`cosmolkit_bio::residue_name_logical_view`).
    cosmolkit_bio::residue_name_logical_view(name, input_format)
}

#[cfg(test)]
mod tests {
    use super::{
        atom_name_text, coordinate_text, float_field_text, int_or_dot, int_or_qmark, number_or_dot,
        number_or_qmark, pdbx_icode, qchain, residue_name_text, string_or_dot, string_or_qmark,
        subchain_or_dot,
    };
    use cosmolkit_bio::{AtomName, BioCoordinateFormat, PdbSeqId, ResidueName, ResidueSourceIds};

    #[test]
    fn bio_pdbscope_namefix_reader_provenance_routes_pdb_trim_and_cif_verbatim() {
        use crate::bio_mmcif::read_mmcif_bio_structure;
        use crate::bio_pdb::{BioPdbReadParams, read_pdb_bio_structure};
        use cosmolkit_bio::BioCoordinateFormat;

        // Real PDB reader path: raw four-column names reach the writer with
        // Pdb provenance, so the Gemmi read_string trim applies (pdb.cpp:33-47,
        // 943) — padded/unpadded/leading/trailing/all-space/maximum-width.
        let pdb_text = concat!(
            "ATOM      1  CA  ALA A   1      11.111  2.222   3.333  1.00 20.00           C  ",
            "\nATOM      2 HD11 ALA A   1      11.111  2.222   3.333  1.00 20.00           H  ",
            "\nHETATM    3 O    HOH B   2      11.111  2.222   3.333  1.00 20.00           O  ",
            "\nEND     \n",
        );
        let pdb =
            read_pdb_bio_structure(pdb_text, "namefix.pdb", &BioPdbReadParams::default()).unwrap();
        assert_eq!(pdb.input_format(), BioCoordinateFormat::Pdb);
        let atoms: Vec<_> = pdb.atoms().to_vec();
        assert_eq!(
            atom_name_text(atoms[0].name(), pdb.input_format()),
            "CA",
            "raw ' CA ' columns trim per read_string"
        );
        assert_eq!(atom_name_text(atoms[1].name(), pdb.input_format()), "HD11");
        assert_eq!(atom_name_text(atoms[2].name(), pdb.input_format()), "O");
        let residues: Vec<_> = pdb.residues().to_vec();
        assert_eq!(
            residue_name_text(&residues[0].name(), pdb.input_format()),
            "ALA"
        );
        assert_eq!(
            residue_name_text(&residues[1].name(), pdb.input_format()),
            "HOH"
        );
        // Same stored raw bytes through a non-Pdb discriminator would be
        // preserved verbatim; the discriminator, not the byte shape, routes.
        assert_eq!(
            atom_name_text(atoms[0].name(), BioCoordinateFormat::Mmcif),
            " CA "
        );

        // Real mmCIF reader path: quoted logical names ' CA ' / ' A ' keep
        // their spaces (mmcif.cpp:248 row.str decoding), verbatim at the
        // writer for Mmcif provenance; the PDB-style trim must not run.
        let cif_text = concat!(
            "data_namefix\n",
            "loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n",
            "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n",
            "_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n",
            "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n",
            "_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
            "ATOM 1 N ' CA ' . ' A ' A 1 ? ? 1.0 2.0 3.0 1.0 20.0\n",
        );
        let cif = read_mmcif_bio_structure(cif_text, "namefix.cif").unwrap();
        assert_eq!(cif.input_format(), BioCoordinateFormat::Mmcif);
        let cif_atoms: Vec<_> = cif.atoms().to_vec();
        assert_eq!(cif_atoms[0].name().as_bytes(), b" CA ");
        assert_eq!(
            atom_name_text(cif_atoms[0].name(), cif.input_format()),
            " CA ",
            "quoted mmCIF logical name keeps its spaces"
        );
        let cif_residues: Vec<_> = cif.residues().to_vec();
        assert_eq!(cif_residues[0].name().as_bytes(), b" A ");
        assert_eq!(
            residue_name_text(&cif_residues[0].name(), cif.input_format()),
            " A "
        );
        // The writer quoting is applied downstream by cif::quote on exactly
        // this preserved text; the helper boundary here returns the text.
    }

    #[test]
    fn bio_pdbscope_v03_names_icode_subchain_and_chain_quoting() {
        // Gemmi read_string (pdb.cpp:33-47) trims raw PDB four-column names
        // only for Pdb provenance; CIF-family decoded logical names keep
        // their spaces verbatim (mmcif.cpp:248 row.str, json.cpp:31).
        assert_eq!(
            atom_name_text(
                AtomName::from_ascii(b" CA ").unwrap(),
                BioCoordinateFormat::Pdb
            ),
            "CA"
        );
        assert_eq!(
            atom_name_text(
                AtomName::from_ascii(b" CA ").unwrap(),
                BioCoordinateFormat::Mmcif
            ),
            " CA "
        );
        assert_eq!(
            atom_name_text(
                AtomName::from_ascii(b"CA").unwrap(),
                BioCoordinateFormat::Mmjson
            ),
            "CA"
        );
        assert_eq!(
            atom_name_text(
                AtomName::from_ascii(b"HD11").unwrap(),
                BioCoordinateFormat::Pdb
            ),
            "HD11"
        );
        assert_eq!(
            atom_name_text(
                AtomName::from_ascii(b"    ").unwrap(),
                BioCoordinateFormat::Pdb
            ),
            ""
        );
        assert_eq!(
            atom_name_text(
                AtomName::from_ascii(b"    ").unwrap(),
                BioCoordinateFormat::ChemComp
            ),
            "    "
        );
        assert_eq!(
            residue_name_text(
                &ResidueName::from_ascii(b"ALA").unwrap(),
                BioCoordinateFormat::Pdb
            ),
            "ALA"
        );
        assert_eq!(
            residue_name_text(
                &ResidueName::from_ascii(b" A ").unwrap(),
                BioCoordinateFormat::Mmcif
            ),
            " A "
        );
        assert_eq!(
            residue_name_text(
                &ResidueName::from_ascii(b" ").unwrap(),
                BioCoordinateFormat::Pdb
            ),
            ""
        );
        assert_eq!(pdbx_icode(None), "?");
        assert_eq!(pdbx_icode(Some(PdbSeqId::new(7, None))), "?");
        assert_eq!(pdbx_icode(Some(PdbSeqId::new(7, Some(b' ')))), "?");
        assert_eq!(pdbx_icode(Some(PdbSeqId::new(7, Some(b'A')))), "A");
        assert_eq!(subchain_or_dot(&ResidueSourceIds::default()), ".");
        assert_eq!(
            subchain_or_dot(
                &ResidueSourceIds::new(None, None, None, Some(String::new()), None).unwrap()
            ),
            "."
        );
        assert_eq!(
            subchain_or_dot(
                &ResidueSourceIds::new(None, None, None, Some("ABC".to_string()), None).unwrap()
            ),
            "ABC"
        );
        assert_eq!(
            subchain_or_dot(
                &ResidueSourceIds::new(None, None, None, Some("two words".to_string()), None)
                    .unwrap()
            ),
            "'two words'"
        );
        assert_eq!(qchain("A"), "A");
        assert_eq!(qchain(""), "''");
        assert_eq!(qchain("two words"), "'two words'");
        // Invalid-text typed error remains the structured channel even though
        // current AtomName/ResidueName/PdbChainId construction makes it
        // unreachable for these fields (ASCII enforced at construction).
        let error = crate::bio_write::BioMmcifWriteError::InvalidText { field: "atom name" };
        assert_eq!(error.to_string(), "invalid UTF-8 in atom name");
    }

    #[test]
    fn bio_pdbscope_v02_sentinels_distinguish_absence_nan_negone_and_empty() {
        assert_eq!(number_or_dot(None), ".");
        assert_eq!(number_or_qmark(None), "?");
        assert_eq!(number_or_dot(Some(f64::NAN)), ".");
        assert_eq!(number_or_qmark(Some(f64::NAN)), "?");
        assert_eq!(number_or_dot(Some(0.0)), "0");
        assert_eq!(number_or_qmark(Some(-0.0)), "-0");
        assert_eq!(number_or_dot(Some(1.5)), "1.5");
        assert_eq!(number_or_qmark(Some(123456789.0)), "123456789");
        assert_eq!(int_or_dot(None), ".");
        assert_eq!(int_or_qmark(None), "?");
        assert_eq!(int_or_dot(Some(-1)), ".");
        assert_eq!(int_or_qmark(Some(-1)), "?");
        assert_eq!(int_or_dot(Some(0)), "0");
        assert_eq!(int_or_qmark(Some(42)), "42");
        assert_eq!(string_or_dot(""), ".");
        assert_eq!(string_or_qmark(""), "?");
        assert_eq!(string_or_dot("ALA"), "ALA");
        assert_eq!(string_or_qmark("two words"), "'two words'");
        assert_eq!(string_or_dot("."), "'.'");
        assert_eq!(string_or_qmark("?"), "'?'");
        assert_eq!(string_or_dot("O'Brien"), "\"O'Brien\"");
    }

    #[test]
    fn bio_pdbscope_v01_coordinate_profile_is_nine_significant_digits() {
        assert_eq!(coordinate_text(0.0), "0");
        assert_eq!(coordinate_text(-0.0), "-0");
        assert_eq!(coordinate_text(1.0), "1");
        assert_eq!(coordinate_text(0.1), "0.1");
        assert_eq!(coordinate_text(-12.75), "-12.75");
        assert_eq!(coordinate_text(0.0001), "0.0001");
        assert_eq!(coordinate_text(0.00001), "1e-05");
        assert_eq!(coordinate_text(1234567890.0), "1.23456789e+09");
        assert_eq!(coordinate_text(123456789.0), "123456789");
        assert_eq!(coordinate_text(0.9999999996), "1");
        assert_eq!(coordinate_text(0.123456789), "0.123456789");
    }

    #[test]
    fn bio_pdbscope_v01_float_fields_use_the_f32_six_digit_profile() {
        assert_eq!(float_field_text(1.0), "1");
        assert_eq!(float_field_text(-0.0), "-0");
        assert_eq!(float_field_text(0.5), "0.5");
        assert_eq!(float_field_text(20.0), "20");
        assert_eq!(float_field_text(3.14159265), "3.14159");
        assert_eq!(float_field_text(0.001), "0.001");
        assert_eq!(float_field_text(0.00001), "1e-05");
        assert_eq!(float_field_text(1234567.0), "1.23457e+06");
        assert_eq!(float_field_text(999999.5), "1e+06");
    }
}
