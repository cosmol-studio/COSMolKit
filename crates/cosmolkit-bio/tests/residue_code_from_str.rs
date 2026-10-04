//! BIO-PARSE1-36: frozen regressions for the restored legacy permissive
//! `FromStr` text parser on `ResidueCode`. Expected values are pasted from
//! the ROOT-frozen ANNEX table / literal contract — never computed by
//! invoking the methods under test beyond the parse itself.

use cosmolkit_bio::{ResidueCode, ResidueCodeParseError};
use std::str::FromStr;

/// Private copy of the exact ROOT-frozen ANNEX368 (variant, ordinal, wire)
/// from residue_code_serde.rs — immutable expected data for this file only.
/// The first tuple field is the Rust variant name, the second the frozen
/// discriminant, the third the canonical wire string.
const ANNEX: [(&str, u16, &str); 368] = [
    ("ALA", 0, "ALA"),
    ("ARG", 1, "ARG"),
    ("ASN", 2, "ASN"),
    ("ABA", 3, "ABA"),
    ("ASP", 4, "ASP"),
    ("ASX", 5, "ASX"),
    ("CYS", 6, "CYS"),
    ("CSH", 7, "CSH"),
    ("GLN", 8, "GLN"),
    ("GLU", 9, "GLU"),
    ("GLX", 10, "GLX"),
    ("GLY", 11, "GLY"),
    ("HIS", 12, "HIS"),
    ("ILE", 13, "ILE"),
    ("LEU", 14, "LEU"),
    ("LYS", 15, "LYS"),
    ("MET", 16, "MET"),
    ("MSE", 17, "MSE"),
    ("ORN", 18, "ORN"),
    ("PHE", 19, "PHE"),
    ("PRO", 20, "PRO"),
    ("SER", 21, "SER"),
    ("THR", 22, "THR"),
    ("TRP", 23, "TRP"),
    ("TYR", 24, "TYR"),
    ("UNK", 25, "UNK"),
    ("VAL", 26, "VAL"),
    ("SEC", 27, "SEC"),
    ("PYL", 28, "PYL"),
    ("SEP", 29, "SEP"),
    ("TPO", 30, "TPO"),
    ("PCA", 31, "PCA"),
    ("CSO", 32, "CSO"),
    ("PTR", 33, "PTR"),
    ("KCX", 34, "KCX"),
    ("CSD", 35, "CSD"),
    ("LLP", 36, "LLP"),
    ("CME", 37, "CME"),
    ("MLY", 38, "MLY"),
    ("DAL", 39, "DAL"),
    ("TYS", 40, "TYS"),
    ("OCS", 41, "OCS"),
    ("M3L", 42, "M3L"),
    ("FME", 43, "FME"),
    ("ALY", 44, "ALY"),
    ("HYP", 45, "HYP"),
    ("CAS", 46, "CAS"),
    ("CRO", 47, "CRO"),
    ("CSX", 48, "CSX"),
    ("DPR", 49, "DPR"),
    ("DGL", 50, "DGL"),
    ("DVA", 51, "DVA"),
    ("CSS", 52, "CSS"),
    ("DPN", 53, "DPN"),
    ("DSN", 54, "DSN"),
    ("DLE", 55, "DLE"),
    ("HIC", 56, "HIC"),
    ("NLE", 57, "NLE"),
    ("MVA", 58, "MVA"),
    ("MLZ", 59, "MLZ"),
    ("CR2", 60, "CR2"),
    ("SAR", 61, "SAR"),
    ("DAR", 62, "DAR"),
    ("DLY", 63, "DLY"),
    ("YCM", 64, "YCM"),
    ("NRQ", 65, "NRQ"),
    ("CGU", 66, "CGU"),
    ("R0TD", 67, "0TD"),
    ("MLE", 68, "MLE"),
    ("DAS", 69, "DAS"),
    ("DTR", 70, "DTR"),
    ("CXM", 71, "CXM"),
    ("TPQ", 72, "TPQ"),
    ("DCY", 73, "DCY"),
    ("DSG", 74, "DSG"),
    ("DTY", 75, "DTY"),
    ("DHI", 76, "DHI"),
    ("MEN", 77, "MEN"),
    ("DTH", 78, "DTH"),
    ("SAC", 79, "SAC"),
    ("DGN", 80, "DGN"),
    ("AIB", 81, "AIB"),
    ("SMC", 82, "SMC"),
    ("IAS", 83, "IAS"),
    ("CIR", 84, "CIR"),
    ("BMT", 85, "BMT"),
    ("DIL", 86, "DIL"),
    ("FGA", 87, "FGA"),
    ("PHI", 88, "PHI"),
    ("CRQ", 89, "CRQ"),
    ("SME", 90, "SME"),
    ("GHP", 91, "GHP"),
    ("MHO", 92, "MHO"),
    ("NEP", 93, "NEP"),
    ("TRQ", 94, "TRQ"),
    ("TOX", 95, "TOX"),
    ("ALC", 96, "ALC"),
    ("R3FG", 97, "3FG"),
    ("SCH", 98, "SCH"),
    ("MDO", 99, "MDO"),
    ("MAA", 100, "MAA"),
    ("GYS", 101, "GYS"),
    ("MK8", 102, "MK8"),
    ("CR8", 103, "CR8"),
    ("KPI", 104, "KPI"),
    ("SCY", 105, "SCY"),
    ("DHA", 106, "DHA"),
    ("OMY", 107, "OMY"),
    ("CAF", 108, "CAF"),
    ("R0AF", 109, "0AF"),
    ("SNN", 110, "SNN"),
    ("MHS", 111, "MHS"),
    ("MLU", 112, "MLU"),
    ("SNC", 113, "SNC"),
    ("PHD", 114, "PHD"),
    ("B3E", 115, "B3E"),
    ("MEA", 116, "MEA"),
    ("MED", 117, "MED"),
    ("OAS", 118, "OAS"),
    ("GL3", 119, "GL3"),
    ("FVA", 120, "FVA"),
    ("PHL", 121, "PHL"),
    ("CRF", 122, "CRF"),
    ("OMZ", 123, "OMZ"),
    ("BFD", 124, "BFD"),
    ("MEQ", 125, "MEQ"),
    ("DAB", 126, "DAB"),
    ("AGM", 127, "AGM"),
    ("PSU", 128, "PSU"),
    ("R5MU", 129, "5MU"),
    ("R7MG", 130, "7MG"),
    ("OMG", 131, "OMG"),
    ("UR3", 132, "UR3"),
    ("OMC", 133, "OMC"),
    ("R2MG", 134, "2MG"),
    ("H2U", 135, "H2U"),
    ("R4SU", 136, "4SU"),
    ("OMU", 137, "OMU"),
    ("R4OC", 138, "4OC"),
    ("MA6", 139, "MA6"),
    ("M2G", 140, "M2G"),
    ("R1MA", 141, "1MA"),
    ("R6MZ", 142, "6MZ"),
    ("CCC", 143, "CCC"),
    ("R2MA", 144, "2MA"),
    ("R1MG", 145, "1MG"),
    ("R5BU", 146, "5BU"),
    ("MIA", 147, "MIA"),
    ("DOC", 148, "DOC"),
    ("R8OG", 149, "8OG"),
    ("R5CM", 150, "5CM"),
    ("R3DR", 151, "3DR"),
    ("BRU", 152, "BRU"),
    ("CBR", 153, "CBR"),
    ("HOH", 154, "HOH"),
    ("DOD", 155, "DOD"),
    ("HEM", 156, "HEM"),
    ("SO4", 157, "SO4"),
    ("GOL", 158, "GOL"),
    ("EDO", 159, "EDO"),
    ("NAG", 160, "NAG"),
    ("PO4", 161, "PO4"),
    ("ACT", 162, "ACT"),
    ("PEG", 163, "PEG"),
    ("MAN", 164, "MAN"),
    ("FAD", 165, "FAD"),
    ("BMA", 166, "BMA"),
    ("ADP", 167, "ADP"),
    ("DMS", 168, "DMS"),
    ("ACE", 169, "ACE"),
    ("NH2", 170, "NH2"),
    ("MPD", 171, "MPD"),
    ("MES", 172, "MES"),
    ("NAD", 173, "NAD"),
    ("NAP", 174, "NAP"),
    ("TRS", 175, "TRS"),
    ("ATP", 176, "ATP"),
    ("PG4", 177, "PG4"),
    ("GDP", 178, "GDP"),
    ("FUC", 179, "FUC"),
    ("FMT", 180, "FMT"),
    ("GAL", 181, "GAL"),
    ("PGE", 182, "PGE"),
    ("FMN", 183, "FMN"),
    ("PLP", 184, "PLP"),
    ("EPE", 185, "EPE"),
    ("SF4", 186, "SF4"),
    ("BME", 187, "BME"),
    ("CIT", 188, "CIT"),
    ("BE7", 189, "BE7"),
    ("MRD", 190, "MRD"),
    ("MHA", 191, "MHA"),
    ("BU3", 192, "BU3"),
    ("PGO", 193, "PGO"),
    ("BU2", 194, "BU2"),
    ("PDO", 195, "PDO"),
    ("BU1", 196, "BU1"),
    ("PG6", 197, "PG6"),
    ("R1BO", 198, "1BO"),
    ("PE7", 199, "PE7"),
    ("PG5", 200, "PG5"),
    ("TFP", 201, "TFP"),
    ("DHD", 202, "DHD"),
    ("PEU", 203, "PEU"),
    ("TAU", 204, "TAU"),
    ("SBT", 205, "SBT"),
    ("SAL", 206, "SAL"),
    ("IOH", 207, "IOH"),
    ("IPA", 208, "IPA"),
    ("PIG", 209, "PIG"),
    ("B3P", 210, "B3P"),
    ("BTB", 211, "BTB"),
    ("NHE", 212, "NHE"),
    ("C8E", 213, "C8E"),
    ("OTE", 214, "OTE"),
    ("PE4", 215, "PE4"),
    ("XPE", 216, "XPE"),
    ("PE8", 217, "PE8"),
    ("P33", 218, "P33"),
    ("N8E", 219, "N8E"),
    ("R2OS", 220, "2OS"),
    ("R1PS", 221, "1PS"),
    ("CPS", 222, "CPS"),
    ("DMX", 223, "DMX"),
    ("MPO", 224, "MPO"),
    ("GCD", 225, "GCD"),
    ("DXG", 226, "DXG"),
    ("CM5", 227, "CM5"),
    ("ACA", 228, "ACA"),
    ("ACN", 229, "ACN"),
    ("CCN", 230, "CCN"),
    ("GLC", 231, "GLC"),
    ("DR6", 232, "DR6"),
    ("NH4", 233, "NH4"),
    ("AZI", 234, "AZI"),
    ("BNG", 235, "BNG"),
    ("BOG", 236, "BOG"),
    ("BGC", 237, "BGC"),
    ("BCN", 238, "BCN"),
    ("BRO", 239, "BRO"),
    ("CAC", 240, "CAC"),
    ("CBX", 241, "CBX"),
    ("ACY", 242, "ACY"),
    ("CBM", 243, "CBM"),
    ("CLO", 244, "CLO"),
    ("R3CO", 245, "3CO"),
    ("NCO", 246, "NCO"),
    ("CU1", 247, "CU1"),
    ("CYN", 248, "CYN"),
    ("MA4", 249, "MA4"),
    ("TAR", 250, "TAR"),
    ("GLO", 251, "GLO"),
    ("MTL", 252, "MTL"),
    ("SOR", 253, "SOR"),
    ("DMU", 254, "DMU"),
    ("DDQ", 255, "DDQ"),
    ("DMF", 256, "DMF"),
    ("DIO", 257, "DIO"),
    ("DOX", 258, "DOX"),
    ("R12P", 259, "12P"),
    ("SDS", 260, "SDS"),
    ("LMT", 261, "LMT"),
    ("EOH", 262, "EOH"),
    ("EEE", 263, "EEE"),
    ("EGL", 264, "EGL"),
    ("FLO", 265, "FLO"),
    ("TRT", 266, "TRT"),
    ("FCY", 267, "FCY"),
    ("FRU", 268, "FRU"),
    ("GBL", 269, "GBL"),
    ("GPX", 270, "GPX"),
    ("HTO", 271, "HTO"),
    ("HTG", 272, "HTG"),
    ("B7G", 273, "B7G"),
    ("C10", 274, "C10"),
    ("R16D", 275, "16D"),
    ("HEZ", 276, "HEZ"),
    ("IOD", 277, "IOD"),
    ("IDO", 278, "IDO"),
    ("ICI", 279, "ICI"),
    ("ICT", 280, "ICT"),
    ("TLA", 281, "TLA"),
    ("LAT", 282, "LAT"),
    ("LBT", 283, "LBT"),
    ("LDA", 284, "LDA"),
    ("MN3", 285, "MN3"),
    ("MRY", 286, "MRY"),
    ("MOH", 287, "MOH"),
    ("BEQ", 288, "BEQ"),
    ("C15", 289, "C15"),
    ("MG8", 290, "MG8"),
    ("POL", 291, "POL"),
    ("NO3", 292, "NO3"),
    ("JEF", 293, "JEF"),
    ("P4C", 294, "P4C"),
    ("CE1", 295, "CE1"),
    ("DIA", 296, "DIA"),
    ("CXE", 297, "CXE"),
    ("IPH", 298, "IPH"),
    ("PIN", 299, "PIN"),
    ("R15P", 300, "15P"),
    ("CRY", 301, "CRY"),
    ("PGR", 302, "PGR"),
    ("PGQ", 303, "PGQ"),
    ("SPD", 304, "SPD"),
    ("SPK", 305, "SPK"),
    ("SPM", 306, "SPM"),
    ("SUC", 307, "SUC"),
    ("TBU", 308, "TBU"),
    ("TMA", 309, "TMA"),
    ("TEP", 310, "TEP"),
    ("SCN", 311, "SCN"),
    ("TRE", 312, "TRE"),
    ("ETF", 313, "ETF"),
    ("R144", 314, "144"),
    ("UMQ", 315, "UMQ"),
    ("URE", 316, "URE"),
    ("YT3", 317, "YT3"),
    ("ZN2", 318, "ZN2"),
    ("FE2", 319, "FE2"),
    ("R3NI", 320, "3NI"),
    ("SIA", 321, "SIA"),
    ("XYP", 322, "XYP"),
    ("A2G", 323, "A2G"),
    ("GLA", 324, "GLA"),
    ("NDG", 325, "NDG"),
    ("NGA", 326, "NGA"),
    ("A", 327, "A"),
    ("C", 328, "C"),
    ("G", 329, "G"),
    ("I", 330, "I"),
    ("U", 331, "U"),
    ("N", 332, "N"),
    ("F", 333, "F"),
    ("K", 334, "K"),
    ("DA", 335, "DA"),
    ("DC", 336, "DC"),
    ("DG", 337, "DG"),
    ("DI", 338, "DI"),
    ("DT", 339, "DT"),
    ("DU", 340, "DU"),
    ("DN", 341, "DN"),
    ("AG", 342, "AG"),
    ("AL", 343, "AL"),
    ("BA", 344, "BA"),
    ("BR", 345, "BR"),
    ("CA", 346, "CA"),
    ("CD", 347, "CD"),
    ("CL", 348, "CL"),
    ("CM", 349, "CM"),
    ("CN", 350, "CN"),
    ("CO", 351, "CO"),
    ("CS", 352, "CS"),
    ("CU", 353, "CU"),
    ("FE", 354, "FE"),
    ("HG", 355, "HG"),
    ("LI", 356, "LI"),
    ("MG", 357, "MG"),
    ("MN", 358, "MN"),
    ("NA", 359, "NA"),
    ("NI", 360, "NI"),
    ("NO", 361, "NO"),
    ("PB", 362, "PB"),
    ("RB", 363, "RB"),
    ("SR", 364, "SR"),
    ("Y1", 365, "Y1"),
    ("ZN", 366, "ZN"),
    ("UNKNOWN", 367, "UNKNOWN"),
];

/// Alternating-case spelling by BYTE index: even index lower, odd index
/// upper, ASCII digits unchanged (input construction only — the expected
/// value stays the frozen ANNEX ordinal).
fn alternating(name: &str) -> String {
    name.bytes()
        .enumerate()
        .map(|(index, byte)| {
            if byte.is_ascii_digit() {
                byte
            } else if index % 2 == 0 {
                byte.to_ascii_lowercase()
            } else {
                byte.to_ascii_uppercase()
            }
        })
        .collect::<Vec<u8>>()
        .into_iter()
        .map(|byte| byte as char)
        .collect()
}

/// Frozen case-sensitive length-2 rows (name, ordinal) from the ROOT
/// source-backed correction — exact literals, never derived from a
/// lookup/parser/Serde output. The length-2 lookup branch is an EXACT
/// match (uppercase D/+ first byte with uppercase second byte, then exact
/// string match), so only the original-uppercase spellings of these rows
/// parse; their lowercase/alternating spellings are rejected.
const CASE_SENSITIVE_LEN2: [(&str, u16); 32] = [
    ("DA", 335),
    ("DC", 336),
    ("DG", 337),
    ("DI", 338),
    ("DT", 339),
    ("DU", 340),
    ("DN", 341),
    ("AG", 342),
    ("AL", 343),
    ("BA", 344),
    ("BR", 345),
    ("CA", 346),
    ("CD", 347),
    ("CL", 348),
    ("CM", 349),
    ("CN", 350),
    ("CO", 351),
    ("CS", 352),
    ("CU", 353),
    ("FE", 354),
    ("HG", 355),
    ("LI", 356),
    ("MG", 357),
    ("MN", 358),
    ("NA", 359),
    ("NI", 360),
    ("NO", 361),
    ("PB", 362),
    ("RB", 363),
    ("SR", 364),
    ("Y1", 365),
    ("ZN", 366),
];

/// Frozen 2220-call MIXED product (source-backed correction): 367 canonical
/// rows x3 spellings x2 repeats (original ASCII, ASCII lowercase,
/// alternating by byte index) = 2202, plus UNKNOWN/unknown/UnKnOwN x2 = 6,
/// plus WAT/H2O/TRY upper+lower each x2 = 12 aliases. The 32 length-2
/// rows' original-uppercase spellings SUCCEED with their literal ANNEX
/// ordinals; their lowercase+alternating spellings are REJECTED (32x2x2 =
/// 128 Err). Every other original canonical spelling succeeds with its
/// literal ordinal. Total 2220 = 2092 Ok + 128 Err. Every Err call checks
/// exact input bytes, literal Display and source=None; every Ok checks the
/// literal ordinal; fresh byte baseline before each parse, compared after
/// the Result BEFORE branching; census increments only after invocation.
#[test]
fn residue_code_from_str_all_rows_mixed_exact_2220() {
    // Collected discrepancy record: every input/outcome/ordinal/error-field
    // mismatch is recorded, never a per-case panic, so ALL 2220 actual
    // calls run before the final censuses and the empty-discrepancy
    // assertion.
    let mut discrepancies: Vec<String> = Vec::new();
    let mut calls = 0usize;
    let mut ok_calls = 0usize;
    let mut err_calls = 0usize;
    for (variant, ordinal, name) in ANNEX {
        if name == "UNKNOWN" {
            continue;
        }
        let case_sensitive = CASE_SENSITIVE_LEN2.iter().any(|(row, _)| *row == name);
        let spellings = [name.to_string(), name.to_lowercase(), alternating(name)];
        for (spelling_index, spelling) in spellings.iter().enumerate() {
            for repeat in 0..2 {
                let label = format!("{variant}#{spelling}#{repeat}");
                let baseline = spelling.as_bytes().to_vec();
                let result = ResidueCode::from_str(&spelling);
                calls += 1;
                if spelling.as_bytes() != baseline.as_slice() {
                    discrepancies.push(format!("{label}: input bytes changed"));
                }
                // The original spelling (index 0) always succeeds; for the
                // 32 case-sensitive length-2 rows the lowercase (1) and
                // alternating (2) spellings are source-rejected Err.
                let expect_ok = !case_sensitive || spelling_index == 0;
                match result {
                    Ok(parsed) => {
                        ok_calls += 1;
                        if !expect_ok {
                            discrepancies.push(format!(
                                "{label}: unexpected Ok({}) for rejected spelling",
                                parsed as u16
                            ));
                        }
                        if parsed as u16 != ordinal {
                            discrepancies.push(format!(
                                "{label}: Ok ordinal {} != literal {ordinal}",
                                parsed as u16
                            ));
                        }
                    }
                    Err(error) => {
                        err_calls += 1;
                        if expect_ok {
                            discrepancies
                                .push(format!("{label}: unexpected Err for accepted spelling"));
                        }
                        if error.input() != spelling.as_str() {
                            discrepancies.push(format!("{label}: Err input bytes differ"));
                        }
                        let expected_display = format!("unknown residue code name '{spelling}'");
                        if error.to_string() != expected_display {
                            discrepancies.push(format!("{label}: Err Display differs"));
                        }
                        let dynamic: &dyn std::error::Error = &error;
                        if dynamic.source().is_some() {
                            discrepancies.push(format!("{label}: Err source not None"));
                        }
                    }
                }
            }
        }
    }

    for spelling in ["UNKNOWN", "unknown", "UnKnOwN"] {
        for repeat in 0..2 {
            let label = format!("sentinel#{spelling}#{repeat}");
            let baseline = spelling.as_bytes().to_vec();
            let result = ResidueCode::from_str(spelling);
            calls += 1;
            if spelling.as_bytes() != baseline.as_slice() {
                discrepancies.push(format!("{label}: input bytes changed"));
            }
            match result {
                Ok(parsed) => {
                    ok_calls += 1;
                    if parsed as u16 != 367 {
                        discrepancies.push(format!(
                            "{label}: sentinel ordinal {} != 367",
                            parsed as u16
                        ));
                    }
                }
                Err(_) => {
                    err_calls += 1;
                    discrepancies.push(format!("{label}: sentinel unexpected Err"))
                }
            }
        }
    }

    for (alias, ordinal) in [
        ("WAT", 154u16),
        ("wat", 154),
        ("H2O", 154),
        ("h2o", 154),
        ("TRY", 23),
        ("try", 23),
    ] {
        for repeat in 0..2 {
            let label = format!("alias#{alias}#{repeat}");
            let baseline = alias.as_bytes().to_vec();
            let result = ResidueCode::from_str(alias);
            calls += 1;
            if alias.as_bytes() != baseline.as_slice() {
                discrepancies.push(format!("{label}: input bytes changed"));
            }
            match result {
                Ok(parsed) => {
                    ok_calls += 1;
                    if parsed as u16 != ordinal {
                        discrepancies.push(format!(
                            "{label}: alias ordinal {} != literal {ordinal}",
                            parsed as u16
                        ));
                    }
                }
                Err(_) => {
                    err_calls += 1;
                    discrepancies.push(format!("{label}: alias unexpected Err"))
                }
            }
        }
    }

    assert_eq!(calls, 2220, "exact 2220-call census");
    assert_eq!(ok_calls, 2092, "exact 2092 Ok census");
    assert_eq!(err_calls, 128, "exact 128 Err census");
    assert!(
        discrepancies.is_empty(),
        "collected discrepancies after all 2220 calls: {discrepancies:?}"
    );
}

/// Frozen 36-call MIXED product (source-backed correction): the exact same
/// 18 strings x2. A -> Ok(327) and AL -> Ok(343) are 4 source-recognized
/// Ok controls with literal ordinals; the other 16 strings (32 calls) are
/// true errors retaining the EXACT original bytes via input(), the Display
/// literal `unknown residue code name '{input}'`, std::error::Error source
/// None, Clone equality, and the error's OWNED input payload surviving the
/// drop of the ORIGINAL caller-owned input String. Each call constructs an
/// OWNED input String from the SAME literal, snapshots its bytes, borrows
/// it into FromStr, compares bytes after the Result, and only THEN drops
/// the ORIGINAL String before inspecting the returned error/clone fields
/// against the unchanged literal. Loops collect discrepancies rather than
/// early-panicking until the final 36/4/32 censuses.
#[test]
fn residue_code_from_str_error_and_recognized_controls_exact_36() {
    const INVALID: [&str; 18] = [
        "",
        " ",
        "ALA ",
        " ALA",
        "\tALA",
        "UNKNOWN ",
        " UNK",
        "????",
        "A",
        "AL",
        "ALAA",
        "\0",
        "ALA\0",
        "é",
        "ＡＬＡ",
        "😀",
        "unknowns",
        "NO_SUCH_RESIDUE",
    ];

    let mut discrepancies: Vec<String> = Vec::new();
    let mut calls = 0usize;
    let mut ok_calls = 0usize;
    let mut err_calls = 0usize;
    for input in INVALID {
        for repeat in 0..2 {
            let label = format!("{input:?}#{repeat}");
            // OWNED original input String from the SAME literal; the error
            // must retain its own payload after this String is dropped.
            let original_input = input.to_string();
            let baseline = original_input.as_bytes().to_vec();
            let result = ResidueCode::from_str(&original_input);
            calls += 1;
            if original_input.as_bytes() != baseline.as_slice() {
                discrepancies.push(format!("{label}: input bytes changed"));
            }
            // Source-recognized controls: the exact length-1/length-2
            // tabulated rows A (327) and AL (343) parse successfully.
            let recognized = match input {
                "A" => Some(327u16),
                "AL" => Some(343u16),
                _ => None,
            };
            // Drop the ORIGINAL caller-owned input String BEFORE error
            // field inspection; every later check compares against the
            // unchanged literal, not against the dropped storage.
            drop(original_input);
            match (recognized, result) {
                (Some(expected), Ok(parsed)) => {
                    ok_calls += 1;
                    if parsed as u16 != expected {
                        discrepancies.push(format!(
                            "{label}: control ordinal {} != literal {expected}",
                            parsed as u16
                        ));
                    }
                }
                (None, Err(error)) => {
                    let error: ResidueCodeParseError = error;
                    err_calls += 1;
                    if error.input() != input {
                        discrepancies.push(format!("{label}: exact original bytes differ"));
                    }
                    let expected_display = format!("unknown residue code name '{input}'");
                    if error.to_string() != expected_display {
                        discrepancies.push(format!("{label}: literal Display differs"));
                    }
                    let dynamic: &dyn std::error::Error = &error;
                    if dynamic.source().is_some() {
                        discrepancies.push(format!("{label}: source not None"));
                    }
                    let clone = error.clone();
                    if clone != error {
                        discrepancies.push(format!("{label}: Clone inequality"));
                    }
                    // The error's OWNED payload still equals the literal
                    // after the ORIGINAL String was dropped above.
                    if error.input() != input {
                        discrepancies.push(format!(
                            "{label}: owned payload changed after original drop"
                        ));
                    }
                }
                (Some(_), Err(_)) => {
                    discrepancies.push(format!("{label}: unexpected Err for control"));
                    err_calls += 1;
                }
                (None, Ok(_)) => {
                    discrepancies.push(format!("{label}: unexpected Ok for invalid input"));
                    ok_calls += 1;
                }
            }
        }
    }
    assert_eq!(calls, 36, "exact 36-call census");
    assert_eq!(ok_calls, 4, "exact 4 Ok control census");
    assert_eq!(err_calls, 32, "exact 32 Err census");
    assert!(
        discrepancies.is_empty(),
        "collected discrepancies after all 36 calls: {discrepancies:?}"
    );
}

/// Frozen parser/Serde separation: six strings, each exercised through BOTH
/// the permissive FromStr and the accepted STRICT serde_json wire with
/// separate 6+6 counters. Wire expectations are frozen independently of the
/// parser: ALA both ALA; ala parse ALA / Serde Err; H2O parse HOH / Serde
/// Err; unknown parse UNKNOWN / Serde Err; UNKNOWN both UNKNOWN; UNK both
/// UNK.
#[test]
fn residue_code_from_str_and_serde_wire_are_separate_contracts() {
    let mut parse_calls = 0usize;
    let mut serde_calls = 0usize;

    let table: [(&str, Option<u16>, Option<u16>); 6] = [
        ("ALA", Some(0), Some(0)),
        ("ala", Some(0), None),
        ("H2O", Some(154), None),
        ("unknown", Some(367), None),
        ("UNKNOWN", Some(367), Some(367)),
        ("UNK", Some(25), Some(25)),
    ];
    for (input, from_str_expected, serde_expected) in table {
        let label = format!("separation#{input}");
        let baseline = input.as_bytes().to_vec();
        let parsed =
            ResidueCode::from_str(input).unwrap_or_else(|error| panic!("{label} parse: {error}"));
        assert_eq!(
            input.as_bytes(),
            baseline.as_slice(),
            "{label}: parse input unchanged"
        );
        assert_eq!(
            parsed as u16,
            from_str_expected.expect("{label}: parser expectation"),
            "{label}: frozen FromStr ordinal"
        );
        parse_calls += 1;

        let json = format!("\"{input}\"");
        let json_baseline = json.as_bytes().to_vec();
        let wire: Result<ResidueCode, _> = serde_json::from_str(&json);
        assert_eq!(
            json.as_bytes(),
            json_baseline.as_slice(),
            "{label}: serde input unchanged"
        );
        match (serde_expected, wire) {
            (Some(expected), Ok(code)) => {
                assert_eq!(code as u16, expected, "{label}: frozen Serde ordinal");
            }
            (None, Err(_)) => {}
            (Some(_), Err(error)) => panic!("{label}: unexpected Serde error: {error}"),
            (None, Ok(code)) => panic!("{label}: unexpected Serde success: {code:?}"),
        }
        serde_calls += 1;
    }

    assert_eq!(parse_calls, 6, "exact 6-call FromStr census");
    assert_eq!(serde_calls, 6, "exact 6-call Serde census");
}
