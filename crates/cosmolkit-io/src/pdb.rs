//! Source-backed PDB coordinate and explicit-topology IO over detached values.
//!
//! The reader owns fixed-column parsing, alternate-location filtering, explicit
//! `CONECT` topology, and multi-model coordinates. Proximity bonding,
//! sanitization, hydrogen removal, and 3D stereochemistry assignment are
//! chemistry operations and intentionally remain outside this format layer.

use std::collections::{BTreeMap, HashMap};

use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomPdbResidueInfo, AtomSpec, Bond, BondId, BondSpec, Conformer2D,
    Conformer3D, CoordinateBlock, CoordinateDimension, MoleculeProperties, TopologyBlock,
};
use cosmolkit_types::{BondOrder, Element};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PdbReadError {
    #[error("PDB contains no atom records after source-defined filtering")]
    Empty,
    #[error("invalid PDB {field} on line {line}: {value}")]
    Field {
        field: &'static str,
        line: usize,
        value: String,
    },
    #[error("detached PDB reader does not support {feature}")]
    Unsupported { feature: &'static str },
    #[error("invalid detached PDB topology: {0}")]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error("invalid detached PDB coordinates: {0}")]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error("PDB chemistry postprocessing failed: {0}")]
    Postprocess(#[from] crate::PdbPostprocessError),
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PdbWriteError {
    #[error("PDB conformer id {id} was not found")]
    ConformerNotFound { id: usize },
    #[error("detached PDB writer does not support {feature}")]
    Unsupported { feature: &'static str },
    #[error("detached PDB topology is invalid: {0}")]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error("detached PDB coordinates are invalid: {0}")]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
}

/// Controls the source-defined PDB block-reader flavor flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PdbReadParams {
    /// Bit 1 keeps alternate conformations and pseudo/dummy records. Bit 8
    /// requests standard-residue bond-order correction, which is not yet
    /// available in the detached format layer.
    pub flavor: u32,
}

/// Controls PDB conformer selection and writer flavor flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PdbWriteParams {
    /// Select a conformer by its detached id. `None` writes all conformers when
    /// there is more than one, matching RDKit's negative `confId` behavior.
    pub conformer_id: Option<usize>,
    /// Bit 2 suppresses single/aromatic `CONECT` entries, bit 4 writes both
    /// directions, bit 8 suppresses bond-order multiplicity, bit 16 adds
    /// `MASTER`, and bit 32 adds `TER`.
    pub flavor: u32,
}

fn bytes_field(line: &str, start: usize, end: usize) -> &str {
    line.as_bytes()
        .get(start..end)
        .and_then(|bytes| std::str::from_utf8(bytes).ok())
        .unwrap_or("")
}

fn byte_at(line: &str, index: usize) -> u8 {
    line.as_bytes().get(index).copied().unwrap_or(0)
}

fn parse_field<T: std::str::FromStr>(
    line: &str,
    start: usize,
    end: usize,
    name: &'static str,
    line_number: usize,
) -> Result<T, PdbReadError> {
    let value = bytes_field(line, start, end);
    value.trim().parse().map_err(|_| PdbReadError::Field {
        field: name,
        line: line_number,
        value: value.to_owned(),
    })
}

fn pdb_atom_from_symbol(symbol: &str) -> Option<(Element, Option<u16>)> {
    // RDKit source: PDBParser.cpp lines 36-49
    // RDKit✔️✔️: Atom *PDBAtomFromSymbol(const char *symb) {
    // RDKit✔️✔️:   PRECONDITION(symb, "bad char ptr");
    // RDKit✔️✔️:   if (symb[0] == 'D' && !symb[1]) {
    // RDKit✔️✔️:     auto *result = new Atom(1);
    // RDKit✔️✔️:     result->setIsotope(2);
    // RDKit✔️✔️:     return result;
    // RDKit✔️✔️:   } else if (symb[0] == 'T' && !symb[1]) {
    // RDKit✔️✔️:     auto *result = new Atom(1);
    // RDKit✔️✔️:     result->setIsotope(3);
    // RDKit✔️✔️:     return result;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   int elemno = PeriodicTable::getTable()->getAtomicNumber(symb);
    // RDKit✔️✔️:   return elemno > 0 ? new Atom(elemno) : (Atom *)nullptr;
    // RDKit✔️✔️: }
    match symbol {
        "D" => Some((Element::H, Some(2))),
        "T" => Some((Element::H, Some(3))),
        _ => Element::from_symbol(symbol).map(|element| (element, None)),
    }
}

fn include_atom_record(line: &str, flavor: u32) -> bool {
    // RDKit source: PDBParser.cpp lines 57-78
    // RDKit✔️✔️:   if (len < 16) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if ((flavor & 1) == 0) {
    // RDKit✔️✔️:     // Ignore alternate locations of atoms.
    // RDKit✔️✔️:     if (len >= 17 && ptr[16] != ' ' && ptr[16] != 'A' && ptr[16] != '1') {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // Ignore XPLOR pseudo atoms
    // RDKit✔️✔️:     if (len >= 54 && !memcmp(ptr + 30, "9999.0009999.0009999.000", 24)) {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // Ignore NMR pseudo atoms
    // RDKit✔️✔️:     if (ptr[12] == ' ' && ptr[13] == 'Q') {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // Ignore PDB dummy residues
    // RDKit✔️✔️:     if (len >= 20 && !memcmp(ptr + 18, "DUM", 3)) {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if line.len() < 16 {
        return false;
    }
    if flavor & 1 == 0 {
        let alt_loc = byte_at(line, 16);
        if line.len() >= 17 && !matches!(alt_loc, b' ' | b'A' | b'1') {
            return false;
        }
        if line.len() >= 54 && bytes_field(line, 30, 54) == "9999.0009999.0009999.000" {
            return false;
        }
        if byte_at(line, 12) == b' ' && byte_at(line, 13) == b'Q' {
            return false;
        }
        if line.len() >= 21 && bytes_field(line, 18, 21) == "DUM" {
            return false;
        }
    }
    true
}

fn element_from_record(
    line: &str,
    serial: i32,
    line_number: usize,
) -> Result<(Element, Option<u16>), PdbReadError> {
    // RDKit source: PDBParser.cpp lines 93-160
    // RDKit✔️✔️:   // Attempt #1:  Atomic Symbol in columns 76 and 77
    // RDKit✔️✔️:   if (len >= 78) {
    // RDKit✔️✔️:     if (ptr[76] >= 'A' && ptr[76] <= 'Z') {
    // RDKit✔️✔️:       symb[0] = ptr[76];
    // RDKit✔️✔️:       if (ptr[77] >= 'A' && ptr[77] <= 'Z') {
    // RDKit✔️✔️:         symb[1] = ptr[77] + 32;  // tolower
    // RDKit✔️✔️:         symb[2] = '\0';
    // RDKit✔️✔️:       } else if (ptr[77] >= 'a' && ptr[77] <= 'z') {
    // RDKit✔️✔️:         symb[1] = ptr[77];
    // RDKit✔️✔️:         symb[2] = '\0';
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         symb[1] = '\0';
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ptr[76] == ' ' && ptr[77] >= 'A' && ptr[77] <= 'Z') {
    // RDKit✔️✔️:       symb[0] = ptr[77];
    // RDKit✔️✔️:       symb[1] = '\0';
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       symb[0] = '\0';
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (len == 77) {
    // RDKit✔️✔️:     if (ptr[76] >= 'A' && ptr[76] <= 'Z') {
    // RDKit✔️✔️:       symb[0] = ptr[76];
    // RDKit✔️✔️:       symb[1] = '\0';
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       symb[0] = '\0';
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     symb[0] = '\0';
    // RDKit✔️✔️:   }
    let mut symbol = [0_u8; 2];
    let mut symbol_len = 0_usize;
    if line.len() >= 78 {
        let first = byte_at(line, 76);
        let second = byte_at(line, 77);
        if first.is_ascii_uppercase() {
            symbol[0] = first;
            symbol_len = 1;
            if second.is_ascii_uppercase() {
                symbol[1] = second.to_ascii_lowercase();
                symbol_len = 2;
            } else if second.is_ascii_lowercase() {
                symbol[1] = second;
                symbol_len = 2;
            }
        } else if first == b' ' && second.is_ascii_uppercase() {
            symbol[0] = second;
            symbol_len = 1;
        }
    } else if line.len() == 77 && byte_at(line, 76).is_ascii_uppercase() {
        symbol[0] = byte_at(line, 76);
        symbol_len = 1;
    }
    let symbol_text = std::str::from_utf8(&symbol[..symbol_len])
        .expect("PDB element symbol was assembled from ASCII bytes");
    if let Some(atom) = pdb_atom_from_symbol(symbol_text) {
        return Ok(atom);
    }

    // RDKit✔️✔️:   if (!atom) {
    // RDKit✔️✔️:     // Attempt #2: Atomic Symbol from PDB atom name
    // RDKit✔️✔️:     if (ptr[13] >= 'A' && ptr[13] <= 'Z') {
    // RDKit✔️✔️:       if (ptr[12] == ' ') {
    // RDKit✔️✔️:         symb[0] = ptr[13];
    // RDKit✔️✔️:         if (ptr[14] >= 'a' && ptr[14] <= 'z') {
    // RDKit✔️✔️:           symb[1] = ptr[14];
    // RDKit✔️✔️:           symb[2] = '\0';
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           symb[1] = '\0';
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else if (ptr[12] >= 'A' && ptr[12] <= 'Z') {
    // RDKit✔️✔️:         symb[0] = ptr[12];
    // RDKit✔️✔️:         symb[1] = ptr[13] + 32;  // tolower
    // RDKit✔️✔️:         symb[2] = '\0';
    // RDKit✔️✔️:         if (ptr[12] == 'H' && ptr[0] == 'A') {
    // RDKit✔️✔️:           // No He, Hf, Hg, Ho or Hs in ATOM records
    // RDKit✔️✔️:           symb[0] = 'H';
    // RDKit✔️✔️:           symb[1] = '\0';
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else if (ptr[12] >= '0' && ptr[12] <= '9') {
    // RDKit✔️✔️:         symb[0] = ptr[13];
    // RDKit✔️✔️:         symb[1] = '\0';
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         symb[0] = '\0';
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       symb[0] = '\0';
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (symb[0]) {
    // RDKit✔️✔️:       atom = PDBAtomFromSymbol(symb);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    symbol_len = 0;
    let first = byte_at(line, 12);
    let second = byte_at(line, 13);
    let third = byte_at(line, 14);
    if second.is_ascii_uppercase() {
        if first == b' ' {
            symbol[0] = second;
            symbol_len = 1;
            if third.is_ascii_lowercase() {
                symbol[1] = third;
                symbol_len = 2;
            }
        } else if first.is_ascii_uppercase() {
            symbol[0] = first;
            symbol[1] = second.to_ascii_lowercase();
            symbol_len = 2;
            if first == b'H' && line.starts_with("ATOM  ") {
                symbol_len = 1;
            }
        } else if first.is_ascii_digit() {
            symbol[0] = second;
            symbol_len = 1;
        }
    }
    let symbol_text = std::str::from_utf8(&symbol[..symbol_len])
        .expect("PDB element symbol was assembled from ASCII bytes");
    pdb_atom_from_symbol(symbol_text).ok_or_else(|| PdbReadError::Field {
        field: "element",
        line: line_number,
        value: format!("atom #{serial}"),
    })
}

fn formal_charge_from_record(line: &str) -> i8 {
    // RDKit source: PDBParser.cpp lines 202-238
    // RDKit✔️✔️:   if (len >= 79) {
    // RDKit✔️✔️:     int charge = 0;
    // RDKit✔️✔️:     if (ptr[78] >= '1' && ptr[78] <= '9') {
    // RDKit✔️✔️:       if (ptr[79] == '-') {
    // RDKit✔️✔️:         charge = -(ptr[78] - '0');
    // RDKit✔️✔️:       } else if (ptr[79] == '+' || ptr[79] == ' ' || !ptr[79]) {
    // RDKit✔️✔️:         charge = ptr[78] - '0';
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ptr[78] == '+') {
    // RDKit✔️✔️:       if (ptr[79] >= '1' && ptr[79] <= '9') {
    // RDKit✔️✔️:         charge = ptr[79] - '0';
    // RDKit✔️✔️:       } else if (ptr[79] == '+') {
    // RDKit✔️✔️:         charge = 2;
    // RDKit✔️✔️:       } else if (ptr[79] != '0') {
    // RDKit✔️✔️:         charge = 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ptr[78] == '-') {
    // RDKit✔️✔️:       if (ptr[79] >= '1' && ptr[79] <= '9') {
    // RDKit✔️✔️:         charge = ptr[79] - '0';
    // RDKit✔️✔️:       } else if (ptr[79] == '-') {
    // RDKit✔️✔️:         charge = -2;
    // RDKit✔️✔️:       } else if (ptr[79] != '0') {
    // RDKit✔️✔️:         charge = -1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ptr[78] == ' ') {
    // RDKit✔️✔️:       if (ptr[79] >= '1' && ptr[79] <= '9') {
    // RDKit✔️✔️:         charge = ptr[79] - '0';
    // RDKit✔️✔️:       } else if (ptr[79] == '+') {
    // RDKit✔️✔️:         charge = 1;
    // RDKit✔️✔️:       } else if (ptr[79] == '-') {
    // RDKit✔️✔️:         charge = -1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (charge != 0) {
    // RDKit✔️✔️:       atom->setFormalCharge(charge);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if line.len() < 79 {
        return 0;
    }
    let first = byte_at(line, 78);
    let second = byte_at(line, 79);
    if (b'1'..=b'9').contains(&first) {
        if second == b'-' {
            -((first - b'0') as i8)
        } else if matches!(second, b'+' | b' ' | 0) {
            (first - b'0') as i8
        } else {
            0
        }
    } else if first == b'+' {
        if (b'1'..=b'9').contains(&second) {
            (second - b'0') as i8
        } else if second == b'+' {
            2
        } else if second != b'0' {
            1
        } else {
            0
        }
    } else if first == b'-' {
        if (b'1'..=b'9').contains(&second) {
            (second - b'0') as i8
        } else if second == b'-' {
            -2
        } else if second != b'0' {
            -1
        } else {
            0
        }
    } else if first == b' ' {
        if (b'1'..=b'9').contains(&second) {
            (second - b'0') as i8
        } else if second == b'+' {
            1
        } else if second == b'-' {
            -1
        } else {
            0
        }
    } else {
        0
    }
}

fn coordinates_from_record(
    line: &str,
    line_number: usize,
) -> Result<Option<[f64; 3]>, PdbReadError> {
    // RDKit source: PDBParser.cpp lines 171-185
    // RDKit✔️✔️:   if (len >= 38) {
    // RDKit✔️✔️:     RDGeom::Point3D pos;
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:       pos.x = FileParserUtils::toDouble(std::string(ptr + 30, 8));
    // RDKit✔️✔️:       if (len >= 46) {
    // RDKit✔️✔️:         pos.y = FileParserUtils::toDouble(std::string(ptr + 38, 8));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (len >= 54) {
    // RDKit✔️✔️:         pos.z = FileParserUtils::toDouble(std::string(ptr + 46, 8));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:       std::ostringstream errout;
    // RDKit✔️✔️:       errout << "Problem with coordinates for PDB atom #" << serialno;
    // RDKit✔️✔️:       throw FileParseException(errout.str());
    // RDKit✔️✔️:     }
    if line.len() < 38 {
        return Ok(None);
    }
    let x = parse_field(line, 30, 38, "x coordinate", line_number)?;
    let y = if line.len() >= 46 {
        parse_field(line, 38, 46, "y coordinate", line_number)?
    } else {
        0.0
    };
    let z = if line.len() >= 54 {
        parse_field(line, 46, 54, "z coordinate", line_number)?
    } else {
        0.0
    };
    Ok(Some([x, y, z]))
}

fn atom_from_record(
    line: &str,
    line_number: usize,
    id: AtomId,
) -> Result<(Atom, i32, Option<[f64; 3]>), PdbReadError> {
    let serial = parse_field(line, 6, 11, "serial", line_number)?;
    let (element, isotope) = element_from_record(line, serial, line_number)?;
    let coordinates = coordinates_from_record(line, line_number)?;
    let residue_number = if line.len() >= 26 {
        parse_field(line, 22, 26, "residue number", line_number)?
    } else {
        1
    };
    let occupancy = if line.len() >= 60 {
        parse_field(line, 54, 60, "occupancy", line_number)?
    } else {
        1.0
    };
    let temp_factor = if line.len() >= 66 {
        parse_field(line, 60, 66, "temperature factor", line_number)?
    } else {
        0.0
    };

    // RDKit source: PDBParser.cpp lines 240-307
    // RDKit✔️✔️:   tmp = std::string(ptr + 12, 4);
    // RDKit✔️✔️:   AtomPDBResidueInfo *info = new AtomPDBResidueInfo(tmp, serialno);
    // RDKit✔️✔️:   atom->setMonomerInfo(info);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (len >= 20) {
    // RDKit✔️✔️:     tmp = std::string(ptr + 17, 3);
    // RDKit✔️✔️:     // boost::trim(tmp);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     tmp = "UNL";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   info->setResidueName(tmp);
    // RDKit✔️✔️:   if (ptr[0] == 'H') {
    // RDKit✔️✔️:     info->setIsHeteroAtom(true);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (len >= 17) {
    // RDKit✔️✔️:     tmp = std::string(ptr + 16, 1);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     tmp = " ";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   info->setAltLoc(tmp);
    // RDKit✔️✔️:   if (len >= 22) {
    // RDKit✔️✔️:     tmp = std::string(ptr + 21, 1);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     tmp = " ";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   info->setChainId(tmp);
    // RDKit✔️✔️:   if (len >= 27) {
    // RDKit✔️✔️:     tmp = std::string(ptr + 26, 1);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     tmp = " ";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   info->setInsertionCode(tmp);
    let residue = AtomPdbResidueInfo::new(
        bytes_field(line, 12, 16),
        serial,
        if line.len() >= 20 {
            bytes_field(line, 17, 20)
        } else {
            "UNL"
        },
        residue_number,
        if line.len() >= 22 {
            bytes_field(line, 21, 22)
        } else {
            " "
        },
        line.starts_with("HETATM"),
    )
    .with_alt_loc(if line.len() >= 17 {
        bytes_field(line, 16, 17)
    } else {
        " "
    })
    .with_insertion_code(if line.len() >= 27 {
        bytes_field(line, 26, 27)
    } else {
        " "
    })
    .with_occupancy(occupancy)
    .with_temp_factor(temp_factor);
    let mut spec = AtomSpec::new(element)
        .with_formal_charge(formal_charge_from_record(line))
        .with_pdb_residue_info(residue);
    if let Some(isotope) = isotope {
        spec = spec.with_isotope(isotope);
    }
    Ok((Atom::from_spec(id, spec), serial, coordinates))
}

fn title_from_record(line: &str, max_len: usize) -> Option<String> {
    // RDKit source: PDBParser.cpp lines 405-420
    // RDKit✔️✔️: void PDBTitleLine(RWMol *mol, const char *ptr, unsigned int len) {
    // RDKit✔️✔️:   PRECONDITION(mol, "bad mol");
    // RDKit✔️✔️:   PRECONDITION(ptr, "bad char ptr");
    // RDKit✔️✔️:   std::string title;
    // RDKit✔️✔️:   while (ptr[len - 1] == ' ') {
    // RDKit✔️✔️:     len--;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (ptr[len - 1] == ';') {
    // RDKit✔️✔️:     len--;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (len > 21 && !strncmp(ptr + 10, " MOLECULE: ", 11)) {
    // RDKit✔️✔️:     title = std::string(ptr + 21, len - 21);
    // RDKit✔️✔️:   } else if (len > 10) {
    // RDKit✔️✔️:     title = std::string(ptr + 10, len - 10);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!title.empty()) {
    // RDKit✔️✔️:     mol->setProp(common_properties::_Name, title);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let bytes = line.as_bytes();
    let mut len = bytes.len().min(max_len);
    while len > 0 && bytes[len - 1] == b' ' {
        len -= 1;
    }
    if len > 0 && bytes[len - 1] == b';' {
        len -= 1;
    }
    let start = if len > 21 && bytes.get(10..21) == Some(b" MOLECULE: ") {
        21
    } else if len > 10 {
        10
    } else {
        return None;
    };
    std::str::from_utf8(&bytes[start..len])
        .ok()
        .filter(|title| !title.is_empty())
        .map(str::to_owned)
}

pub(super) fn same_pdb_residue(left: &AtomPdbResidueInfo, right: &AtomPdbResidueInfo) -> bool {
    // RDKit source: ProximityBonds.cpp `SamePDBResidue`
    // RDKit✔️✔️: bool SamePDBResidue(AtomPDBResidueInfo *p, AtomPDBResidueInfo *q) {
    // RDKit✔️✔️:   return p->getResidueNumber() == q->getResidueNumber() &&
    // RDKit✔️✔️:          p->getResidueName() == q->getResidueName() &&
    // RDKit✔️✔️:          p->getChainId() == q->getChainId() &&
    // RDKit✔️✔️:          p->getInsertionCode() == q->getInsertionCode();
    // RDKit✔️✔️: }
    left.residue_number() == right.residue_number()
        && left.residue_name() == right.residue_name()
        && left.chain_id() == right.chain_id()
        && left.insertion_code() == right.insertion_code()
}

fn is_blacklisted_atom(atomic_number: u8) -> bool {
    // RDKit source: ProximityBonds.cpp `IsBlacklistedAtom`
    // RDKit✔️✔️: static bool IsBlacklistedAtom(Atom *atom) {
    // RDKit✔️✔️:   // blacklist metals, noble gasses and halogens
    // RDKit✔️✔️:   int elem = atom->getAtomicNum();
    // RDKit✔️✔️:   // make an inverse query (non-metals and metaloids)
    // RDKit✔️✔️:   return !((5 <= elem && elem <= 8) || (14 <= elem && elem <= 16) ||
    // RDKit✔️✔️:            (32 <= elem && elem <= 34) || (51 <= elem && elem <= 52));
    // RDKit✔️✔️: }
    !((5..=8).contains(&atomic_number)
        || (14..=16).contains(&atomic_number)
        || (32..=34).contains(&atomic_number)
        || (51..=52).contains(&atomic_number))
}

pub(super) fn is_blacklisted_pair(atoms: &[Atom], begin: AtomId, end: AtomId) -> bool {
    // RDKit source: ProximityBonds.cpp `IsBlacklistedPair`
    // RDKit✔️✔️:   if (!SamePDBResidue(beg_info, end_info)) {
    // RDKit✔️✔️:     if (IsBlacklistedAtom(beg_atom) || IsBlacklistedAtom(end_atom)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // Dont make bonds to waters
    // RDKit✔️✔️:     if (beg_info->getResidueName() == "HOH" ||
    // RDKit✔️✔️:         end_info->getResidueName() == "HOH") {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    let Some(begin_atom) = atoms.get(begin.index()) else {
        return false;
    };
    let Some(end_atom) = atoms.get(end.index()) else {
        return false;
    };
    let (Some(begin_info), Some(end_info)) =
        (begin_atom.pdb_residue_info(), end_atom.pdb_residue_info())
    else {
        return false;
    };
    !same_pdb_residue(begin_info, end_info)
        && (is_blacklisted_atom(begin_atom.atomic_number())
            || is_blacklisted_atom(end_atom.atomic_number())
            || begin_info.residue_name() == "HOH"
            || end_info.residue_name() == "HOH")
}

fn apply_conect_target(
    atoms: &[Atom],
    serial_to_index: &HashMap<i32, AtomId>,
    bonds: &mut Vec<Bond>,
    bond_by_atoms: &mut HashMap<(usize, usize), usize>,
    seen_by_bond: &mut HashMap<usize, u8>,
    source: i32,
    target: i32,
) -> bool {
    let (Some(&begin), Some(&end)) = (serial_to_index.get(&source), serial_to_index.get(&target))
    else {
        return true;
    };
    if source == target {
        return true;
    }
    let key = (
        begin.index().min(end.index()),
        begin.index().max(end.index()),
    );
    if let Some(&bond_index) = bond_by_atoms.get(&key) {
        if bonds[bond_index].order() == BondOrder::Zero {
            return false;
        }
        // RDKit source: PDBParser.cpp lines 349-382
        // RDKit✔️✔️:         // Here we use a single byte bitmap to count duplicates
        // RDKit✔️✔️:         // Low nibble counts src < dst, high nibble for src > dst
        // RDKit✔️✔️:         int seen = bmap[bond];
        // RDKit✔️✔️:         if (src < dst) {
        // RDKit✔️✔️:           if ((seen & 0x0f) == 0x01) {
        // RDKit✔️✔️:             bmap[bond] = seen | 0x02;
        // RDKit✔️✔️:             if ((seen & 0x20) == 0) {
        // RDKit✔️✔️:               bond->setBondType(Bond::DOUBLE);
        // RDKit✔️✔️:             }
        // RDKit✔️✔️:           } else if ((seen & 0x0f) == 0x03) {
        // RDKit✔️✔️:             bmap[bond] = seen | 0x04;
        // RDKit✔️✔️:             if ((seen & 0x40) == 0) {
        // RDKit✔️✔️:               bond->setBondType(Bond::TRIPLE);
        // RDKit✔️✔️:             }
        // RDKit✔️✔️:           } else if ((seen & 0x0f) == 0x07) {
        // RDKit✔️✔️:             bmap[bond] = seen | 0x08;
        // RDKit✔️✔️:             if ((seen & 0x80) == 0) {
        // RDKit✔️✔️:               bond->setBondType(Bond::QUADRUPLE);
        // RDKit✔️✔️:             }
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         } else /* src < dst */ {
        // RDKit✔️✔️:           if ((seen & 0xf0) == 0x10) {
        // RDKit✔️✔️:             bmap[bond] = seen | 0x20;
        // RDKit✔️✔️:             if ((seen & 0x02) == 0) {
        // RDKit✔️✔️:               bond->setBondType(Bond::DOUBLE);
        // RDKit✔️✔️:             }
        // RDKit✔️✔️:           } else if ((seen & 0xf0) == 0x30) {
        // RDKit✔️✔️:             bmap[bond] = seen | 0x40;
        // RDKit✔️✔️:             if ((seen & 0x04) == 0) {
        // RDKit✔️✔️:               bond->setBondType(Bond::TRIPLE);
        // RDKit✔️✔️:             }
        // RDKit✔️✔️:           } else if ((seen & 0xf0) == 0x70) {
        // RDKit✔️✔️:             bmap[bond] = seen | 0x80;
        // RDKit✔️✔️:             if ((seen & 0x08) == 0) {
        // RDKit✔️✔️:               bond->setBondType(Bond::QUADRUPLE);
        // RDKit✔️✔️:             }
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        let seen = seen_by_bond.get(&bond_index).copied().unwrap_or(0);
        let (next_seen, next_order) = if source < target {
            match seen & 0x0f {
                0x01 => (seen | 0x02, (seen & 0x20 == 0).then_some(BondOrder::Double)),
                0x03 => (seen | 0x04, (seen & 0x40 == 0).then_some(BondOrder::Triple)),
                0x07 => (
                    seen | 0x08,
                    (seen & 0x80 == 0).then_some(BondOrder::Quadruple),
                ),
                _ => (seen, None),
            }
        } else {
            match seen & 0xf0 {
                0x10 => (seen | 0x20, (seen & 0x02 == 0).then_some(BondOrder::Double)),
                0x30 => (seen | 0x40, (seen & 0x04 == 0).then_some(BondOrder::Triple)),
                0x70 => (
                    seen | 0x80,
                    (seen & 0x08 == 0).then_some(BondOrder::Quadruple),
                ),
                _ => (seen, None),
            }
        };
        seen_by_bond.insert(bond_index, next_seen);
        if let Some(order) = next_order {
            bonds[bond_index].set_order(order);
        }
        return true;
    }

    // RDKit✔️✔️:       } else if (!bond) {
    // RDKit✔️✔️:         // Bonds in PDB file are explicit
    // RDKit✔️✔️:         // if they are not sanitize friendly, set their order to zero
    // RDKit✔️✔️:         if (IsBlacklistedPair(amap[src], amap[dst])) {
    // RDKit✔️✔️:           bond = new Bond(Bond::ZERO);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           bond = new Bond(Bond::SINGLE);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         mol->addBond(bond, true);
    // RDKit✔️✔️:         bmap[bond] = (src < dst) ? 0x01 : 0x10;
    let order = if is_blacklisted_pair(atoms, begin, end) {
        BondOrder::Zero
    } else {
        BondOrder::Single
    };
    let bond_index = bonds.len();
    bonds.push(Bond::from_spec(
        BondId::new(bond_index),
        BondSpec::new(begin, end, order),
    ));
    bond_by_atoms.insert(key, bond_index);
    seen_by_bond.insert(bond_index, if source < target { 0x01 } else { 0x10 });
    true
}

fn apply_conect_record(
    line: &str,
    line_number: usize,
    atoms: &[Atom],
    serial_to_index: &HashMap<i32, AtomId>,
    bonds: &mut Vec<Bond>,
    bond_by_atoms: &mut HashMap<(usize, usize), usize>,
    seen_by_bond: &mut HashMap<usize, u8>,
) -> Result<(), PdbReadError> {
    // RDKit source: PDBParser.cpp lines 315-344
    // RDKit✔️✔️:   if (len < 16) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string tmp(ptr + 6, 5);
    // RDKit✔️✔️:   bool fail = false;
    // RDKit✔️✔️:   int src, dst;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     src = FileParserUtils::toInt(tmp);
    // RDKit✔️✔️:     if (amap.find(src) == amap.end()) {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:     fail = true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!fail) {
    // RDKit✔️✔️:     if (len > 41) {
    // RDKit✔️✔️:       len = 41;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int pos = 11; pos + 5 <= len; pos += 5) {
    // RDKit✔️✔️:       if (!memcmp(ptr + pos, "     ", 5)) {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    if line.len() < 16 {
        return Ok(());
    }
    let source: i32 = parse_field(line, 6, 11, "CONECT source", line_number)?;
    if !serial_to_index.contains_key(&source) {
        return Ok(());
    }
    let len = line.len().min(41);
    let mut position = 11;
    while position + 5 <= len {
        let value = bytes_field(line, position, position + 5);
        if value == "     " {
            break;
        }
        let target = value.trim().parse().map_err(|_| PdbReadError::Field {
            field: "CONECT target",
            line: line_number,
            value: value.to_owned(),
        })?;
        if !apply_conect_target(
            atoms,
            serial_to_index,
            bonds,
            bond_by_atoms,
            seen_by_bond,
            source,
            target,
        ) {
            break;
        }
        position += 5;
    }
    Ok(())
}

#[derive(Debug)]
struct WorkingConformer {
    coordinates: Vec<[f64; 3]>,
    is_3d: bool,
}

impl WorkingConformer {
    fn new(atom_count: usize) -> Self {
        Self {
            coordinates: vec![[0.0, 0.0, 0.0]; atom_count],
            is_3d: false,
        }
    }

    fn set(&mut self, atom: usize, point: [f64; 3]) {
        self.coordinates[atom] = point;
        self.is_3d |= point[2] != 0.0;
    }
}

/// Read source-aligned PDB fixed-column state and explicit `CONECT` topology.
pub fn read_pdb_detached_with_params(
    block: &str,
    params: PdbReadParams,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), PdbReadError> {
    let mut atoms = Vec::new();
    let mut serial_to_index = HashMap::new();
    let mut bonds = Vec::new();
    let mut bond_by_atoms = HashMap::new();
    let mut seen_by_bond = HashMap::new();
    let mut conformers = Vec::<WorkingConformer>::new();
    let mut primary_conformer = None::<usize>;
    let mut multi_conformer = false;
    let mut conformer_atom_index = 0_usize;
    let mut active_multi_conformer = None::<usize>;
    let mut title = None;

    // RDKit source: PDBParser.cpp `parsePdbBlock` record dispatch
    // RDKit✔️✔️:     // ATOM records
    // RDKit✔️✔️:     if (str[0] == 'A' && str[1] == 'T' && str[2] == 'O' && str[3] == 'M' &&
    // RDKit✔️✔️:         str[4] == ' ' && str[5] == ' ') {
    // RDKit✔️✔️:       if (!multi_conformer) {
    // RDKit✔️✔️:         if (!mol) {
    // RDKit✔️✔️:           mol.reset(new RWMol());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         PDBAtomLine(mol.get(), str, len, flavor, amap);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         PDBConformerLine(mol.get(), str, len, conf, conformer_atmidx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // HETATM records
    // RDKit✔️✔️:     } else if (str[0] == 'H' && str[1] == 'E' && str[2] == 'T' &&
    // RDKit✔️✔️:                str[3] == 'A' && str[4] == 'T' && str[5] == 'M') {
    // RDKit✔️✔️:       if (!multi_conformer) {
    // RDKit✔️✔️:         if (!mol) {
    // RDKit✔️✔️:           mol.reset(new RWMol());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         PDBAtomLine(mol.get(), str, len, flavor, amap);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         PDBConformerLine(mol.get(), str, len, conf, conformer_atmidx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // CONECT records
    for (line_index, line) in block.lines().enumerate() {
        let line_number = line_index + 1;
        if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
            if !multi_conformer {
                if !include_atom_record(line, params.flavor) {
                    continue;
                }
                let id = AtomId::new(atoms.len());
                let (atom, serial, point) = atom_from_record(line, line_number, id)?;
                atoms.push(atom);
                serial_to_index.insert(serial, id);
                if let Some(point) = point {
                    let conformer = *primary_conformer.get_or_insert_with(|| {
                        conformers.push(WorkingConformer::new(atoms.len()));
                        conformers.len() - 1
                    });
                    if conformers[conformer].coordinates.len() < atoms.len() {
                        conformers[conformer]
                            .coordinates
                            .resize(atoms.len(), [0.0; 3]);
                    }
                    conformers[conformer].set(id.index(), point);
                } else if let Some(conformer) = primary_conformer {
                    conformers[conformer].coordinates.push([0.0; 3]);
                }
            } else if let Some(point) = coordinates_from_record(line, line_number)? {
                if conformer_atom_index >= atoms.len() {
                    continue;
                }
                let conformer = *active_multi_conformer.get_or_insert_with(|| {
                    conformers.push(WorkingConformer::new(atoms.len()));
                    conformers.len() - 1
                });
                conformers[conformer].set(conformer_atom_index, point);
                conformer_atom_index += 1;
            }
        } else if line.starts_with("CONECT") {
            if !multi_conformer {
                apply_conect_record(
                    line,
                    line_number,
                    &atoms,
                    &serial_to_index,
                    &mut bonds,
                    &mut bond_by_atoms,
                    &mut seen_by_bond,
                )?;
            }
        } else if line.starts_with("COMPND") {
            if line.len() > 10
                && (byte_at(line, 9) == b' ' || line.as_bytes().get(9..21) == Some(b"2 MOLECULE: "))
                && let Some(value) = title_from_record(line, line.len())
            {
                title = Some(value);
            }
        } else if line.starts_with("HEADER") {
            if let Some(value) = title_from_record(line, 50) {
                title = Some(value);
            }
        } else if line.starts_with("ENDMDL") {
            if atoms.is_empty() {
                break;
            }
            multi_conformer = true;
            conformer_atom_index = 0;
            active_multi_conformer = None;
        }
    }

    if atoms.is_empty() {
        return Err(PdbReadError::Empty);
    }
    let mut topology = TopologyBlock {
        adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
        atoms,
        bonds,
        ..TopologyBlock::default()
    };
    topology.validate()?;
    let coordinate_block = CoordinateBlock {
        conformers_2d: Vec::new(),
        conformers_3d: conformers
            .into_iter()
            .enumerate()
            .map(|(id, conformer)| Conformer3D::new(id, conformer.coordinates, conformer.is_3d))
            .collect(),
        source_coordinate_dim: Some(CoordinateDimension::ThreeD),
    };
    coordinate_block.validate_for_atom_count(topology.atoms.len())?;
    crate::postprocess_pdb_detached(
        &mut topology,
        &coordinate_block,
        crate::PdbPostprocessParams {
            proximity_bonding: false,
            flavor: params.flavor,
        },
    )?;
    let properties = title.map_or_else(MoleculeProperties::default, |value| {
        MoleculeProperties::default().with_name(value)
    });
    Ok((topology, coordinate_block, properties))
}

/// Read PDB fixed-column state and explicit `CONECT` topology with flavor 0.
pub fn read_pdb_detached(
    block: &str,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), PdbReadError> {
    read_pdb_detached_with_params(block, PdbReadParams::default())
}

fn default_atom_number(atomic_number: u8, counts: &mut BTreeMap<u8, u32>) -> String {
    // RDKit source: PDBWriter.cpp lines 141-165
    // RDKit✔️✔️: std::string GetDefaultAtomNumber(const Atom *atom,
    // RDKit✔️✔️:                                  std::map<unsigned int, unsigned int> &elem) {
    // RDKit✔️✔️:   std::string ret = "  ";
    // RDKit✔️✔️:   unsigned int atno = atom->getAtomicNum();
    // RDKit✔️✔️:   if (elem.find(atno) == elem.end()) {
    // RDKit✔️✔️:     elem[atno] = 1;
    // RDKit✔️✔️:     ret[0] = '1';
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     unsigned int tmp = elem[atno] + 1;
    // RDKit✔️✔️:     elem[atno] = tmp;
    // RDKit✔️✔️:     if (tmp < 10) {
    // RDKit✔️✔️:       ret[0] = tmp + '0';
    // RDKit✔️✔️:     } else if (tmp < 100) {
    // RDKit✔️✔️:       ret[0] = (tmp / 10) + '0';
    // RDKit✔️✔️:       ret[1] = (tmp % 10) + '0';
    // RDKit✔️✔️:     } else if (tmp < 360) {
    // RDKit✔️✔️:       ret[0] = ((tmp - 100) / 10) + 'A';
    // RDKit✔️✔️:       ret[1] = ((tmp - 100) % 10) + '0';
    // RDKit✔️✔️:     } else if (tmp < 1036) {
    // RDKit✔️✔️:       ret[0] = ((tmp - 360) / 26) + 'A';
    // RDKit✔️✔️:       ret[1] = ((tmp - 360) % 26) + 'A';
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return ret;
    // RDKit✔️✔️: }
    let count = counts.entry(atomic_number).or_default();
    *count += 1;
    let value = *count;
    let mut result = *b"  ";
    if value < 10 {
        result[0] = b'0' + value as u8;
    } else if value < 100 {
        result[0] = b'0' + (value / 10) as u8;
        result[1] = b'0' + (value % 10) as u8;
    } else if value < 360 {
        result[0] = b'A' + ((value - 100) / 10) as u8;
        result[1] = b'0' + ((value - 100) % 10) as u8;
    } else if value < 1036 {
        result[0] = b'A' + ((value - 360) / 26) as u8;
        result[1] = b'A' + ((value - 360) % 26) as u8;
    }
    String::from_utf8(result.to_vec()).expect("default PDB atom name is ASCII")
}

#[derive(Clone, Copy)]
enum SelectedConformer<'a> {
    TwoD(&'a Conformer2D),
    ThreeD(&'a Conformer3D),
}

impl SelectedConformer<'_> {
    fn point(self, atom: usize) -> [f64; 3] {
        match self {
            Self::TwoD(conformer) => {
                let point = conformer.coordinates()[atom];
                [point[0], point[1], 0.0]
            }
            Self::ThreeD(conformer) => conformer.coordinates()[atom],
        }
    }
}

fn pdb_atom_line(
    atom: &Atom,
    point: Option<[f64; 3]>,
    element_counts: &mut BTreeMap<u8, u32>,
) -> String {
    // RDKit source: PDBWriter.cpp lines 50-68
    // RDKit✔️✔️:   std::string symb = atom->getSymbol();
    // RDKit✔️✔️:   char at1, at2;
    // RDKit✔️✔️:   switch (symb.length()) {
    // RDKit✔️✔️:     case 0:
    // RDKit✔️✔️:       at1 = ' ';
    // RDKit✔️✔️:       at2 = 'X';
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 1:
    // RDKit✔️✔️:       at1 = ' ';
    // RDKit✔️✔️:       at2 = symb[0];
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       at1 = symb[0];
    // RDKit✔️✔️:       at2 = symb[1];
    // RDKit✔️✔️:       if (at2 >= 'a' && at2 <= 'z') {
    // RDKit✔️✔️:         at2 -= 32;  // toupper
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    let symbol = atom.element().symbol().as_bytes();
    let (element_first, element_second) = match symbol {
        [] => (b' ', b'X'),
        [only] => (b' ', *only),
        [first, second, ..] => (*first, second.to_ascii_uppercase()),
    };
    let mut line = String::with_capacity(80);
    if let Some(info) = atom.pdb_residue_info() {
        // RDKit source: PDBWriter.cpp lines 70-100
        // RDKit✔️✔️:   if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
        // RDKit✔️✔️:     ss << (info->getIsHeteroAtom() ? "HETATM" : "ATOM  ");
        // RDKit✔️✔️:     ss << std::setw(5) << atom->getIdx() + 1;
        // RDKit✔️✔️:     ss << ' ';
        // RDKit✔️✔️:     const std::string &name = info->getName();
        // RDKit✔️✔️:     if (name.empty()) {
        // RDKit✔️✔️:       std::string atnum = GetDefaultAtomNumber(atom, elem);
        // RDKit✔️✔️:       ss << at1 << at2 << atnum;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       ss << std::setw(4) << name.substr(0, 4);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     const char *ptr = info->getAltLoc().c_str();
        // RDKit✔️✔️:     if (*ptr == '\0') {
        // RDKit✔️✔️:       ptr = " ";
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ss << *ptr;
        // RDKit✔️✔️:     ss << std::setw(3) << info->getResidueName().substr(0, 3);
        // RDKit✔️✔️:     ss << ' ';
        // RDKit✔️✔️:     ptr = info->getChainId().c_str();
        // RDKit✔️✔️:     if (*ptr == '\0') {
        // RDKit✔️✔️:       ptr = " ";
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ss << *ptr;
        // RDKit✔️✔️:     ss << std::setw(4) << info->getResidueNumber();
        // RDKit✔️✔️:     ptr = info->getInsertionCode().c_str();
        // RDKit✔️✔️:     if (*ptr == '\0') {
        // RDKit✔️✔️:       ptr = " ";
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ss << *ptr;
        // RDKit✔️✔️:     ss << "   ";
        line.push_str(if info.is_hetero_atom() {
            "HETATM"
        } else {
            "ATOM  "
        });
        line.push_str(&format!("{:>5} ", atom.id().index() + 1));
        if info.atom_name().is_empty() {
            line.push(char::from(element_first));
            line.push(char::from(element_second));
            line.push_str(&default_atom_number(atom.atomic_number(), element_counts));
        } else {
            let name = info.atom_name().chars().take(4).collect::<String>();
            line.push_str(&format!("{name:>4}"));
        }
        line.push(info.alt_loc().chars().next().unwrap_or(' '));
        let residue = info.residue_name().chars().take(3).collect::<String>();
        line.push_str(&format!("{residue:>3} "));
        line.push(info.chain_id().chars().next().unwrap_or(' '));
        line.push_str(&format!("{:>4}", info.residue_number()));
        line.push(info.insertion_code().chars().next().unwrap_or(' '));
        line.push_str("   ");
    } else {
        // RDKit source: PDBWriter.cpp lines 101-109
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     info = (AtomPDBResidueInfo *)nullptr;
        // RDKit✔️✔️:     std::string atnum = GetDefaultAtomNumber(atom, elem);
        // RDKit✔️✔️:     ss << "HETATM";
        // RDKit✔️✔️:     ss << std::setw(5) << atom->getIdx() + 1;
        // RDKit✔️✔️:     ss << ' ';
        // RDKit✔️✔️:     ss << at1 << at2 << atnum;
        // RDKit✔️✔️:     ss << " UNL     1    ";
        // RDKit✔️✔️:   }
        line.push_str(&format!("HETATM{:>5} ", atom.id().index() + 1));
        line.push(char::from(element_first));
        line.push(char::from(element_second));
        line.push_str(&default_atom_number(atom.atomic_number(), element_counts));
        line.push_str(" UNL     1    ");
    }

    // RDKit source: PDBWriter.cpp lines 111-124
    // RDKit✔️✔️:   if (conf) {
    // RDKit✔️✔️:     const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    // RDKit✔️✔️:     ss << boost::format("%8.3f%8.3f%8.3f") % pos.x % pos.y % pos.z;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     ss << "   0.000   0.000   0.000";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (info) {
    // RDKit✔️✔️:     ss << boost::format("%6.2f%6.2f") % info->getOccupancy() %
    // RDKit✔️✔️:               info->getTempFactor();
    // RDKit✔️✔️:     ss << "          ";
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     ss << "  1.00  0.00          ";
    // RDKit✔️✔️:   }
    let point = point.unwrap_or([0.0; 3]);
    line.push_str(&format!(
        "{:>8.3}{:>8.3}{:>8.3}",
        point[0], point[1], point[2]
    ));
    if let Some(info) = atom.pdb_residue_info() {
        line.push_str(&format!(
            "{:>6.2}{:>6.2}",
            info.occupancy(),
            info.temp_factor()
        ));
        line.push_str("          ");
    } else {
        line.push_str("  1.00  0.00          ");
    }
    line.push(char::from(element_first));
    line.push(char::from(element_second));
    let charge = atom.formal_charge();
    if (1..10).contains(&charge) {
        line.push(char::from(b'0' + charge as u8));
        line.push('+');
    } else if (-9..0).contains(&charge) {
        line.push(char::from(b'0' + (-charge) as u8));
        line.push('-');
    } else {
        line.push_str("  ");
    }
    line
}

fn pdb_bond_lines(
    atom: AtomId,
    bonds: &[Bond],
    all: bool,
    both: bool,
    multiplicity: bool,
    conect_count: &mut usize,
) -> String {
    // RDKit source: PDBWriter.cpp lines 167-225
    // RDKit✔️✔️: std::string GetPDBBondLines(const Atom *atom, bool all, bool both, bool mult,
    // RDKit✔️✔️:                             unsigned int &conect_count) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   unsigned int src = atom->getIdx() + 1;
    // RDKit✔️✔️:   std::vector<unsigned int> v;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   ROMol *mol = &atom->getOwningMol();
    // RDKit✔️✔️:   for (ROMol::OBOND_ITER_PAIR bondIt = mol->getAtomBonds(atom);
    // RDKit✔️✔️:        bondIt.first != bondIt.second; ++bondIt.first) {
    // RDKit✔️✔️:     Bond *bptr = (*mol)[*bondIt.first];
    // RDKit✔️✔️:     Atom *nptr = bptr->getOtherAtom(atom);
    // RDKit✔️✔️:     unsigned int dst = nptr->getIdx() + 1;
    // RDKit✔️✔️:     if (dst < src && !both) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     Bond::BondType btype = Bond::SINGLE;
    // RDKit✔️✔️:     if (mult) {
    // RDKit✔️✔️:       btype = bptr->getBondType();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     switch (btype) {
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:       case Bond::SINGLE:
    // RDKit✔️✔️:       case Bond::AROMATIC:
    // RDKit✔️✔️:         if (all) {
    // RDKit✔️✔️:           v.push_back(dst);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Bond::QUADRUPLE:
    // RDKit✔️✔️:         v.push_back(dst);
    // RDKit✔️✔️:         /* FALLTHRU */
    // RDKit✔️✔️:       case Bond::TRIPLE:
    // RDKit✔️✔️:         v.push_back(dst);
    // RDKit✔️✔️:         /* FALLTHRU */
    // RDKit✔️✔️:       case Bond::DOUBLE:
    // RDKit✔️✔️:         v.push_back(dst);
    // RDKit✔️✔️:         v.push_back(dst);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let source = atom.index() + 1;
    let mut destinations = Vec::new();
    for bond in bonds {
        let destination = if bond.begin() == atom {
            bond.end().index() + 1
        } else if bond.end() == atom {
            bond.begin().index() + 1
        } else {
            continue;
        };
        if destination < source && !both {
            continue;
        }
        let order = if multiplicity {
            bond.order()
        } else {
            BondOrder::Single
        };
        let copies = match order {
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            _ if all => 1,
            _ => 0,
        };
        destinations.extend(std::iter::repeat_n(destination, copies));
    }
    destinations.sort_unstable();
    let mut output = String::new();
    for chunk in destinations.chunks(4) {
        output.push_str(&format!("CONECT{source:>5}"));
        *conect_count += 1;
        for destination in chunk {
            output.push_str(&format!("{destination:>5}"));
        }
        output.push('\n');
    }
    output
}

fn pdb_body(
    topology: &TopologyBlock,
    conformer: Option<SelectedConformer<'_>>,
    flavor: u32,
    atom_count: &mut usize,
    ter_count: &mut usize,
    conect_count: &mut usize,
) -> String {
    let mut output = String::new();
    let mut last = String::new();
    let mut element_counts = BTreeMap::new();
    for atom in &topology.atoms {
        last = pdb_atom_line(
            atom,
            conformer.map(|selected| selected.point(atom.id().index())),
            &mut element_counts,
        );
        output.push_str(&last);
        output.push('\n');
        *atom_count += 1;
    }

    // RDKit source: PDBWriter.cpp lines 242-263
    // RDKit✔️✔️:   if (ter_count == 0 && atm_count && (flavor & 32)) {
    // RDKit✔️✔️:     std::stringstream ss;
    // RDKit✔️✔️:     ss << "TER   ";
    // RDKit✔️✔️:     ss << std::setw(5) << atm_count + 1;
    // RDKit✔️✔️:     if (last.length() >= 27) {
    // RDKit✔️✔️:       ss << "      ";
    // RDKit✔️✔️:       ss << last.substr(17, 10);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ss << '\n';
    // RDKit✔️✔️:     res += ss.str();
    // RDKit✔️✔️:     ter_count = 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool all = (flavor & 2) == 0;
    // RDKit✔️✔️:   bool both = (flavor & 4) != 0;
    // RDKit✔️✔️:   bool mult = (flavor & 8) == 0;
    if *ter_count == 0 && *atom_count > 0 && flavor & 32 != 0 {
        output.push_str(&format!("TER   {:>5}", *atom_count + 1));
        if last.len() >= 27 {
            output.push_str("      ");
            output.push_str(&last[17..27]);
        }
        output.push('\n');
        *ter_count = 1;
    }
    let all = flavor & 2 == 0;
    let both = flavor & 4 != 0;
    let multiplicity = flavor & 8 == 0;
    if all || multiplicity {
        for atom in &topology.atoms {
            output.push_str(&pdb_bond_lines(
                atom.id(),
                &topology.bonds,
                all,
                both,
                multiplicity,
                conect_count,
            ));
        }
    }
    output
}

fn selected_conformers(
    coordinates: &CoordinateBlock,
    id: Option<usize>,
) -> Result<Vec<SelectedConformer<'_>>, PdbWriteError> {
    if let Some(id) = id {
        return coordinates
            .conformers_3d
            .iter()
            .find(|conformer| conformer.id() == id)
            .map(SelectedConformer::ThreeD)
            .or_else(|| {
                coordinates
                    .conformers_2d
                    .iter()
                    .find(|conformer| conformer.id() == id)
                    .map(SelectedConformer::TwoD)
            })
            .map(|conformer| vec![conformer])
            .ok_or(PdbWriteError::ConformerNotFound { id });
    }
    if !coordinates.conformers_2d.is_empty() && !coordinates.conformers_3d.is_empty() {
        return Err(PdbWriteError::Unsupported {
            feature: "implicit ordering of mixed 2D and 3D conformer stores",
        });
    }
    Ok(coordinates
        .conformers_3d
        .iter()
        .map(SelectedConformer::ThreeD)
        .chain(
            coordinates
                .conformers_2d
                .iter()
                .map(SelectedConformer::TwoD),
        )
        .collect())
}

/// Write a source-aligned PDB block with explicit conformer/flavor controls.
pub fn write_pdb_detached_with_params(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    params: PdbWriteParams,
) -> Result<String, PdbWriteError> {
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    if topology
        .bonds
        .iter()
        .any(|bond| bond.order() == BondOrder::Aromatic || bond.is_aromatic())
    {
        // RDKit source: PDBWriter.cpp lines 267-270
        // RDKit❌❌: std::string MolToPDBBlock(const ROMol &imol, int confId, unsigned int flavor) {
        // RDKit❌❌:   RWMol rwmol(imol);
        // RDKit❌❌:   MolOps::Kekulize(rwmol);
        // RDKit❌❌:   Utils::LocaleSwitcher ls;
        return Err(PdbWriteError::Unsupported {
            feature: "RDKit Kekulize preprocessing for aromatic PDB output",
        });
    }
    let conformers = selected_conformers(coordinates, params.conformer_id)?;
    let mut output = String::new();
    if let Some(name) = properties.name().filter(|name| !name.is_empty()) {
        // RDKit✔️✔️:   if (rwmol.getPropIfPresent(common_properties::_Name, name)) {
        // RDKit✔️✔️:     if (!name.empty()) {
        // RDKit✔️✔️:       res += "COMPND    ";
        // RDKit✔️✔️:       res += name;
        // RDKit✔️✔️:       res += '\n';
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        output.push_str("COMPND    ");
        output.push_str(name);
        output.push('\n');
    }
    let mut atom_count = 0;
    let mut ter_count = 0;
    let mut conect_count = 0;
    if params.conformer_id.is_none() && conformers.len() > 1 {
        // RDKit source: PDBWriter.cpp lines 286-299
        // RDKit✔️✔️:   if (confId < 0 && rwmol.getNumConformers() > 1) {
        // RDKit✔️✔️:     int count = rwmol.getNumConformers();
        // RDKit✔️✔️:     for (confId = 0; confId < count; confId++) {
        // RDKit✔️✔️:       conf = &(rwmol.getConformer(confId));
        // RDKit✔️✔️:       std::stringstream ss;
        // RDKit✔️✔️:       ss << "MODEL     ";
        // RDKit✔️✔️:       ss << std::setw(4) << (confId + 1);
        // RDKit✔️✔️:       ss << "\n";
        // RDKit✔️✔️:       res += ss.str();
        // RDKit✔️✔️:       res +=
        // RDKit✔️✔️:           MolToPDBBody(rwmol, conf, flavor, atm_count, ter_count, conect_count);
        // RDKit✔️✔️:       res += "ENDMDL\n";
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        for (index, conformer) in conformers.iter().copied().enumerate() {
            output.push_str(&format!("MODEL     {:>4}\n", index + 1));
            output.push_str(&pdb_body(
                topology,
                Some(conformer),
                params.flavor,
                &mut atom_count,
                &mut ter_count,
                &mut conect_count,
            ));
            output.push_str("ENDMDL\n");
        }
    } else {
        output.push_str(&pdb_body(
            topology,
            conformers.first().copied(),
            params.flavor,
            &mut atom_count,
            &mut ter_count,
            &mut conect_count,
        ));
    }
    if params.flavor & 16 != 0 {
        // RDKit source: PDBWriter.cpp lines 310-317
        // RDKit✔️✔️:   if (flavor & 16) {
        // RDKit✔️✔️:     std::stringstream ss;
        // RDKit✔️✔️:     ss << "MASTER        0    0    0    0    0    0    0    0";
        // RDKit✔️✔️:     ss << std::setw(5) << atm_count;
        // RDKit✔️✔️:     ss << std::setw(5) << ter_count;
        // RDKit✔️✔️:     ss << std::setw(5) << conect_count;
        // RDKit✔️✔️:     ss << "    0\n";
        // RDKit✔️✔️:     res += ss.str();
        // RDKit✔️✔️:   }
        output.push_str(&format!(
            "MASTER        0    0    0    0    0    0    0    0{atom_count:>5}{ter_count:>5}{conect_count:>5}    0\n"
        ));
    }
    output.push_str("END\n");
    Ok(output)
}

/// Write a PDB block with RDKit's negative-conformer-id and flavor-0 defaults.
pub fn write_pdb_detached(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> Result<String, PdbWriteError> {
    write_pdb_detached_with_params(topology, coordinates, properties, PdbWriteParams::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pdb_atom(
        record: &str,
        serial: i32,
        name: &str,
        alt_loc: char,
        residue: &str,
        chain: char,
        residue_number: i32,
        point: [f64; 3],
        occupancy: f64,
        temp_factor: f64,
        element: &str,
        charge: &str,
    ) -> String {
        format!(
            "{record}{serial:>5} {name:<4}{alt_loc}{residue:>3} {chain}{residue_number:>4}    {x:>8.3}{y:>8.3}{z:>8.3}{occupancy:>6.2}{temp_factor:>6.2}          {element:>2}{charge:>2}",
            x = point[0],
            y = point[1],
            z = point[2],
        )
    }

    #[test]
    fn reader_preserves_fixed_column_atom_state_and_rdkit_element_rules() {
        let oxygen = pdb_atom(
            "ATOM  ",
            19,
            " O  ",
            'A',
            "HOH",
            'B',
            7,
            [1.0, 2.0, 0.0],
            0.5,
            12.25,
            "O",
            "1-",
        );
        let deuterium = pdb_atom(
            "HETATM",
            22,
            " D  ",
            ' ',
            "LIG",
            'C',
            8,
            [0.0, 0.0, -0.0],
            1.0,
            0.0,
            "D",
            "",
        );
        let input = format!(
            "HEADER    FIXED STATE                                      \n{oxygen}\n{deuterium}\nEND\n"
        );
        let (topology, coordinates, properties) = read_pdb_detached(&input).expect("read PDB");
        assert_eq!(topology.atoms.len(), 2);
        assert_eq!(topology.atoms[0].formal_charge(), -1);
        assert_eq!(topology.atoms[1].element(), Element::H);
        assert_eq!(topology.atoms[1].isotope(), Some(2));
        let info = topology.atoms[0].pdb_residue_info().expect("residue info");
        assert_eq!(info.atom_name(), " O  ");
        assert_eq!(info.alt_loc(), "A");
        assert_eq!(info.occupancy(), 0.5);
        assert_eq!(info.temp_factor(), 12.25);
        assert_eq!(properties.name(), Some("FIXED STATE"));
        assert!(!coordinates.conformers_3d[0].is_3d());
        assert!(coordinates.conformers_3d[0].coordinates()[1][2].is_sign_negative());
    }

    #[test]
    fn reader_filters_default_altloc_and_pseudo_records_but_flavor_one_keeps_them() {
        let normal = pdb_atom(
            "ATOM  ", 1, " C  ", ' ', "GLY", 'A', 1, [0.0; 3], 1.0, 0.0, "C", "",
        );
        let alternate = pdb_atom(
            "ATOM  ",
            2,
            " C  ",
            'B',
            "GLY",
            'A',
            1,
            [1.0, 0.0, 0.0],
            1.0,
            0.0,
            "C",
            "",
        );
        let nmr = pdb_atom(
            "ATOM  ",
            3,
            " Q1 ",
            ' ',
            "GLY",
            'A',
            1,
            [2.0, 0.0, 0.0],
            1.0,
            0.0,
            "C",
            "",
        );
        let input = format!("{normal}\n{alternate}\n{nmr}\n");
        assert_eq!(
            read_pdb_detached(&input).expect("filtered").0.atoms.len(),
            1
        );
        assert_eq!(
            read_pdb_detached_with_params(&input, PdbReadParams { flavor: 1 })
                .expect("unfiltered")
                .0
                .atoms
                .len(),
            3
        );
    }

    #[test]
    fn reader_decodes_directional_conect_multiplicity_and_blacklisted_zero_bonds() {
        let carbon = pdb_atom(
            "HETATM", 10, " C1 ", ' ', "LIG", 'A', 1, [0.0; 3], 1.0, 0.0, "C", "",
        );
        let oxygen = pdb_atom(
            "HETATM",
            20,
            " O1 ",
            ' ',
            "LIG",
            'A',
            1,
            [1.0, 0.0, 0.0],
            1.0,
            0.0,
            "O",
            "",
        );
        let sodium = pdb_atom(
            "HETATM",
            30,
            "NA  ",
            ' ',
            "NA ",
            'B',
            2,
            [2.0, 0.0, 0.0],
            1.0,
            0.0,
            "NA",
            "",
        );
        let input = format!(
            "{carbon}\n{oxygen}\n{sodium}\nCONECT   10   20   20   20\nCONECT   20   10\nCONECT   30   20\n"
        );
        let topology = read_pdb_detached(&input).expect("PDB").0;
        assert_eq!(topology.bonds.len(), 2);
        assert_eq!(topology.bonds[0].order(), BondOrder::Triple);
        assert_eq!(topology.bonds[1].order(), BondOrder::Zero);
    }

    #[test]
    fn reader_collects_following_models_as_ordered_conformers() {
        let first_a = pdb_atom(
            "ATOM  ",
            1,
            " C  ",
            ' ',
            "GLY",
            'A',
            1,
            [0.0, 0.0, 0.0],
            1.0,
            0.0,
            "C",
            "",
        );
        let first_b = pdb_atom(
            "ATOM  ",
            2,
            " O  ",
            ' ',
            "GLY",
            'A',
            1,
            [1.0, 0.0, 0.0],
            1.0,
            0.0,
            "O",
            "",
        );
        let second_a = pdb_atom(
            "ATOM  ",
            1,
            " C  ",
            ' ',
            "GLY",
            'A',
            1,
            [0.0, 1.0, 0.0],
            1.0,
            0.0,
            "C",
            "",
        );
        let second_b = pdb_atom(
            "ATOM  ",
            2,
            " O  ",
            ' ',
            "GLY",
            'A',
            1,
            [1.0, 1.0, 2.0],
            1.0,
            0.0,
            "O",
            "",
        );
        let input = format!(
            "MODEL        1\n{first_a}\n{first_b}\nENDMDL\nMODEL        2\n{second_a}\n{second_b}\nENDMDL\n"
        );
        let coordinates = read_pdb_detached(&input).expect("models").1;
        assert_eq!(coordinates.conformers_3d.len(), 2);
        assert_eq!(
            coordinates.conformers_3d[1].coordinates()[1],
            [1.0, 1.0, 2.0]
        );
        assert!(coordinates.conformers_3d[1].is_3d());
    }

    #[test]
    fn reader_flavor_eight_applies_standard_residue_bond_orders() {
        let carbon = pdb_atom(
            "ATOM  ", 1, " C  ", ' ', "ALA", 'A', 1, [0.0; 3], 1.0, 0.0, "C", "",
        );
        let oxygen = pdb_atom(
            "ATOM  ",
            2,
            " O  ",
            ' ',
            "ALA",
            'A',
            1,
            [1.2, 0.0, 0.0],
            1.0,
            0.0,
            "O",
            "",
        );
        let input = format!("{carbon}\n{oxygen}\nCONECT    1    2\n");
        let topology = read_pdb_detached_with_params(&input, PdbReadParams { flavor: 8 })
            .expect("flavor-8 PDB")
            .0;
        assert_eq!(topology.bonds.len(), 1);
        assert_eq!(topology.bonds[0].order(), BondOrder::Double);
    }

    #[test]
    fn detached_postprocess_adds_proximity_bonds_and_ignores_hydrogen_contacts() {
        let carbon = pdb_atom(
            "HETATM", 1, " C1 ", ' ', "LIG", 'A', 1, [0.0; 3], 1.0, 0.0, "C", "",
        );
        let oxygen = pdb_atom(
            "HETATM",
            2,
            " O1 ",
            ' ',
            "LIG",
            'A',
            1,
            [1.2, 0.0, 0.0],
            1.0,
            0.0,
            "O",
            "",
        );
        let hydrogen_a = pdb_atom(
            "HETATM",
            3,
            " H1 ",
            ' ',
            "LIG",
            'A',
            1,
            [10.0, 0.0, 0.0],
            1.0,
            0.0,
            "H",
            "",
        );
        let hydrogen_b = pdb_atom(
            "HETATM",
            4,
            " H2 ",
            ' ',
            "LIG",
            'A',
            1,
            [10.7, 0.0, 0.0],
            1.0,
            0.0,
            "H",
            "",
        );
        let input = format!("{carbon}\n{oxygen}\n{hydrogen_a}\n{hydrogen_b}\n");
        let (mut topology, coordinates, _) = read_pdb_detached(&input).expect("PDB syntax");
        crate::postprocess_pdb_detached(
            &mut topology,
            &coordinates,
            crate::PdbPostprocessParams {
                proximity_bonding: true,
                flavor: 0,
            },
        )
        .expect("proximity bonding");
        assert_eq!(topology.bonds.len(), 1);
        assert_eq!(topology.bonds[0].begin(), AtomId::new(1));
        assert_eq!(topology.bonds[0].end(), AtomId::new(0));
        assert_eq!(topology.bonds[0].order(), BondOrder::Single);
    }

    #[test]
    fn detached_postprocess_matches_rdkit_multivalent_h_and_neutral_n_cleanup() {
        let hydrogen = pdb_atom(
            "HETATM", 1, " H1 ", ' ', "LIG", 'A', 1, [0.0; 3], 1.0, 0.0, "H", "",
        );
        let carbon = pdb_atom(
            "HETATM",
            2,
            " C1 ",
            ' ',
            "LIG",
            'A',
            1,
            [1.0, 0.0, 0.0],
            1.0,
            0.0,
            "C",
            "",
        );
        let oxygen = pdb_atom(
            "HETATM",
            3,
            " O1 ",
            ' ',
            "LIG",
            'A',
            1,
            [-1.1, 0.0, 0.0],
            1.0,
            0.0,
            "O",
            "",
        );
        let (mut topology, coordinates, _) =
            read_pdb_detached(&format!("{hydrogen}\n{carbon}\n{oxygen}\n")).expect("PDB");
        crate::postprocess_pdb_detached(
            &mut topology,
            &coordinates,
            crate::PdbPostprocessParams {
                proximity_bonding: true,
                flavor: 0,
            },
        )
        .expect("multivalent-H cleanup");
        assert_eq!(topology.bonds.len(), 1);
        assert_eq!(topology.bonds[0].begin(), AtomId::new(1));
        assert_eq!(topology.bonds[0].end(), AtomId::new(0));

        let nitrogen = pdb_atom(
            "HETATM", 1, " N1 ", ' ', "LIG", 'A', 1, [0.0; 3], 1.0, 0.0, "N", "",
        );
        let neighbors = [
            pdb_atom(
                "HETATM",
                2,
                " C1 ",
                ' ',
                "LIG",
                'A',
                1,
                [1.4, 0.0, 0.0],
                1.0,
                0.0,
                "C",
                "",
            ),
            pdb_atom(
                "HETATM",
                3,
                " C2 ",
                ' ',
                "LIG",
                'A',
                1,
                [-1.4, 0.0, 0.0],
                1.0,
                0.0,
                "C",
                "",
            ),
            pdb_atom(
                "HETATM",
                4,
                " C3 ",
                ' ',
                "LIG",
                'A',
                1,
                [0.0, 1.4, 0.0],
                1.0,
                0.0,
                "C",
                "",
            ),
            pdb_atom(
                "HETATM",
                5,
                " C4 ",
                ' ',
                "LIG",
                'A',
                1,
                [0.0, -1.4, 0.0],
                1.0,
                0.0,
                "C",
                "",
            ),
        ];
        let input = format!(
            "{nitrogen}\n{}\n{}\n{}\n{}\nCONECT    1    2    3    4    5\n",
            neighbors[0], neighbors[1], neighbors[2], neighbors[3]
        );
        let topology = read_pdb_detached(&input)
            .expect("four-coordinate nitrogen")
            .0;
        assert_eq!(topology.atoms[0].formal_charge(), 1);
    }

    #[test]
    fn writer_matches_rdkit_atom_numbering_metadata_charge_and_conect_shape() {
        let (mut topology, coordinates, properties) = read_pdb_detached(&format!(
            "{}\n{}\nCONECT   41   99   99\n",
            pdb_atom(
                "ATOM  ",
                41,
                " N  ",
                ' ',
                "GLY",
                'A',
                4,
                [1.25, -2.5, 0.0],
                0.75,
                11.5,
                "N",
                "1+",
            ),
            pdb_atom(
                "HETATM",
                99,
                " O1 ",
                ' ',
                "LIG",
                'B',
                8,
                [0.0, 0.0, 3.0],
                1.0,
                2.0,
                "O",
                "",
            ),
        ))
        .expect("seed PDB");
        topology.bonds[0].set_order(BondOrder::Double);
        let output = write_pdb_detached(&topology, &coordinates, &properties).expect("write PDB");
        assert_eq!(
            output,
            "ATOM      1  N   GLY A   4       1.250  -2.500   0.000  0.75 11.50           N1+\n\
             HETATM    2  O1  LIG B   8       0.000   0.000   3.000  1.00  2.00           O  \n\
             CONECT    1    2    2\n\
             END\n"
        );
    }

    #[test]
    fn writer_handles_multiple_conformers_and_flavor_records() {
        let atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
        let topology = TopologyBlock {
            atoms: vec![atom],
            adjacency: AdjacencyList::from_topology(1, &[]),
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(4, vec![[1.0, 2.0, 3.0]], true),
                Conformer3D::new(9, vec![[4.0, 5.0, 6.0]], true),
            ],
            ..CoordinateBlock::default()
        };
        let output = write_pdb_detached_with_params(
            &topology,
            &coordinates,
            &MoleculeProperties::default().with_name("models"),
            PdbWriteParams {
                conformer_id: None,
                flavor: 16 | 32,
            },
        )
        .expect("write models");
        assert!(output.contains("COMPND    models\nMODEL        1\n"));
        assert!(output.contains("ENDMDL\nMODEL        2\n"));
        assert_eq!(
            output
                .lines()
                .filter(|line| line.starts_with("TER"))
                .count(),
            1
        );
        assert!(output.contains("MASTER"));
        let selected = write_pdb_detached_with_params(
            &topology,
            &coordinates,
            &MoleculeProperties::default(),
            PdbWriteParams {
                conformer_id: Some(9),
                flavor: 0,
            },
        )
        .expect("select conformer by id");
        assert!(selected.contains("   4.000   5.000   6.000"));
    }

    #[test]
    fn writer_reproduces_rdkit_conect_flavor_interactions() {
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
        )];
        let topology = TopologyBlock {
            atoms,
            adjacency: AdjacencyList::from_topology(2, &bonds),
            bonds,
            ..TopologyBlock::default()
        };
        let write = |flavor| {
            write_pdb_detached_with_params(
                &topology,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
                PdbWriteParams {
                    conformer_id: None,
                    flavor,
                },
            )
            .expect("write flavored PDB")
        };
        assert!(write(2).contains("CONECT    1    2    2\n"));
        assert!(write(4).contains("CONECT    2    1    1\n"));
        assert!(write(8).contains("CONECT    1    2\n"));
        assert!(!write(10).contains("CONECT"));
    }

    #[test]
    fn unsupported_source_dependent_paths_fail_closed() {
        let atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
        let topology = TopologyBlock {
            atoms: vec![atom],
            adjacency: AdjacencyList::from_topology(1, &[]),
            ..TopologyBlock::default()
        };
        let missing = write_pdb_detached_with_params(
            &topology,
            &CoordinateBlock::default(),
            &MoleculeProperties::default(),
            PdbWriteParams {
                conformer_id: Some(8),
                flavor: 0,
            },
        )
        .expect_err("missing conformer");
        assert_eq!(missing, PdbWriteError::ConformerNotFound { id: 8 });
    }
}
