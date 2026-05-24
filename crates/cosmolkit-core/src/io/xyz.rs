//! RDKit-aligned XYZ reader.
//!
//! XYZ files contain atom identities and Cartesian coordinates only. This
//! reader intentionally does not infer connectivity or bond orders.

use std::num::ParseFloatError;

use crate::{AtomSpec, Element, Molecule, MoleculeBuildError, MoleculeBuilder};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum XyzReadError {
    #[error("empty XYZ block")]
    EmptyBlock,
    #[error(
        "unable to recognize the number of atoms: cannot convert '{value}' to unsigned int on line 0"
    )]
    AtomCount { value: String },
    #[error("EOF hit while reading atoms")]
    UnexpectedEof,
    #[error("missing coordinates on line {line}")]
    MissingCoordinates { line: usize },
    #[error("cannot convert '{value}' to double on line {line}")]
    Coordinate {
        value: String,
        line: usize,
        #[source]
        source: ParseFloatError,
    },
    #[error("{message}")]
    AtomSymbol { message: String },
    #[error(transparent)]
    Build(#[from] MoleculeBuildError),
}

fn parse_float_error() -> ParseFloatError {
    "x".parse::<f64>().unwrap_err()
}

fn rdkit_to_unsigned(value: &str) -> Result<usize, XyzReadError> {
    // BEGIN RDKIT CPP FUNCTION FileParserUtils::toUnsigned(input, true)
    // RDKit✔️✔️: const char *txt = input.data();
    // RDKit✔️✔️: for (size_t i = 0u; i < input.size() && *txt != '\x00'; ++i) {
    // RDKit✔️✔️:   if ((*txt >= '0' && *txt <= '9') || (acceptSpaces && *txt == ' ') ||
    // RDKit✔️✔️:       *txt == '+') {
    // RDKit✔️✔️:     ++txt;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     throw boost::bad_lexical_cast();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let checked = value.split_once('\0').map_or(value, |(prefix, _)| prefix);
    if !checked
        .bytes()
        .all(|byte| byte.is_ascii_digit() || byte == b' ' || byte == b'+')
    {
        return Err(XyzReadError::AtomCount {
            value: value.to_string(),
        });
    }

    // RDKit✔️✔️: txt = input.data();
    // RDKit✔️✔️: unsigned int sz = input.size();
    // RDKit✔️✔️: if (acceptSpaces) {
    // RDKit✔️✔️:   while (*txt == ' ') {
    // RDKit✔️✔️:     ++txt;
    // RDKit✔️✔️:     --sz;
    // RDKit✔️✔️:     if (sz < 1U) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let trimmed_start = checked.trim_start_matches(' ');
    if checked.len() != 0 && trimmed_start.is_empty() {
        return Ok(0);
    }

    // RDKit✔️✔️: unsigned int res = 0;
    // RDKit✔️✔️: std::from_chars(txt, txt + sz, res);
    let digit_prefix = trimmed_start.bytes().take_while(u8::is_ascii_digit).count();
    if digit_prefix == 0 {
        return Ok(0);
    }
    trimmed_start[..digit_prefix]
        .parse::<usize>()
        .map_err(|_| XyzReadError::AtomCount {
            value: value.to_string(),
        })
    // END RDKIT CPP FUNCTION
}

fn rdkit_to_double_no_spaces(value: &str, line: usize) -> Result<f64, XyzReadError> {
    // BEGIN RDKIT CPP FUNCTION FileParserUtils::toDouble(input, false)
    // RDKit✔️✔️: const char *txt = input.data();
    // RDKit✔️✔️: for (size_t i = 0u; i < input.size() && *txt != '\x00'; ++i) {
    // RDKit✔️✔️:   if ((*txt >= '0' && *txt <= '9') || (acceptSpaces && *txt == ' ') ||
    // RDKit✔️✔️:       *txt == '+' || *txt == '-' || *txt == ',' || *txt == '.') {
    // RDKit✔️✔️:     ++txt;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     throw boost::bad_lexical_cast();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let checked = value.split_once('\0').map_or(value, |(prefix, _)| prefix);
    if !checked.bytes().all(|byte| {
        byte.is_ascii_digit() || byte == b'+' || byte == b'-' || byte == b',' || byte == b'.'
    }) {
        return Err(XyzReadError::Coordinate {
            value: value.to_string(),
            line,
            source: value.parse::<f64>().err().unwrap_or_else(parse_float_error),
        });
    }

    // RDKit✔️✔️: double res = atof(input.data());
    let numeric_prefix = checked
        .split_once(',')
        .map_or(checked, |(prefix, _)| prefix);
    numeric_prefix
        .parse::<f64>()
        .map_err(|source| XyzReadError::Coordinate {
            value: value.to_string(),
            line,
            source,
        })
    // END RDKIT CPP FUNCTION
}

#[must_use]
fn normalize_xyz_symbol(raw: &str) -> String {
    let mut chars = raw.chars().collect::<Vec<_>>();
    if chars.len() == 2 && chars[1].is_ascii_uppercase() {
        chars[1] = chars[1].to_ascii_lowercase();
    }
    chars.into_iter().collect()
}

fn element_from_rdkit_symbol(symbol: &str) -> Result<Element, XyzReadError> {
    for atomic_number in 0..=118u8 {
        if crate::chemistry::valence::rdkit_element_symbol(atomic_number).ok() == Some(symbol) {
            return Element::from_atomic_number(atomic_number).ok_or_else(|| {
                XyzReadError::AtomSymbol {
                    message: format!("Bad atomic number for element symbol {symbol}"),
                }
            });
        }
    }
    Err(XyzReadError::AtomSymbol {
        message: format!("Element '{symbol}' not found"),
    })
}

fn parse_xyz_atom_line(atom_line: &str, line: usize) -> Result<(Element, [f64; 3]), XyzReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseXYZFileAtomLine
    // RDKit✔️✔️: std::string whitespace{" \t"};
    // RDKit✔️✔️: size_t delims[8];
    // RDKit✔️✔️: size_t prev = 0;
    // RDKit✔️✔️: for (unsigned int i = 0; i < 7; i++) {
    // RDKit✔️✔️:   if (i % 2 == 0) {
    // RDKit✔️✔️:     delims[i] = atomLine.find_first_not_of(whitespace, prev);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     delims[i] = atomLine.find_first_of(whitespace, prev);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (delims[i] == std::string::npos) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Missing coordinates on line " << line << std::endl;
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   prev = delims[i];
    // RDKit✔️✔️: }
    // RDKit✔️✔️: delims[7] = atomLine.find_last_not_of(whitespace) + 1;
    let mut delims = [0usize; 8];
    let mut prev = 0usize;
    for (i, delim) in delims.iter_mut().take(7).enumerate() {
        *delim = if i % 2 == 0 {
            atom_line[prev..]
                .find(|ch: char| ch != ' ' && ch != '\t')
                .map(|idx| prev + idx)
        } else {
            atom_line[prev..].find([' ', '\t']).map(|idx| prev + idx)
        }
        .ok_or(XyzReadError::MissingCoordinates { line })?;
        prev = *delim;
    }
    delims[7] = atom_line
        .rfind(|ch: char| ch != ' ' && ch != '\t')
        .map_or(0, |idx| idx + 1);

    // RDKit✔️✔️: pos.x = FileParserUtils::toDouble(
    // RDKit✔️✔️:     atomLine.substr(delims[2], delims[3] - delims[2]), false);
    // RDKit✔️✔️: pos.y = FileParserUtils::toDouble(
    // RDKit✔️✔️:     atomLine.substr(delims[4], delims[5] - delims[4]), false);
    // RDKit✔️✔️: pos.z = FileParserUtils::toDouble(
    // RDKit✔️✔️:     atomLine.substr(delims[6], delims[7] - delims[6]), false);
    let pos = [
        rdkit_to_double_no_spaces(&atom_line[delims[2]..delims[3]], line)?,
        rdkit_to_double_no_spaces(&atom_line[delims[4]..delims[5]], line)?,
        rdkit_to_double_no_spaces(&atom_line[delims[6]..delims[7]], line)?,
    ];

    // RDKit✔️✔️: std::string symb{atomLine.substr(delims[0], delims[1] - delims[0])};
    // RDKit✔️✔️: if (symb.size() == 2 && symb[1] >= 'A' && symb[1] <= 'Z') {
    // RDKit✔️✔️:   symb[1] = static_cast<char>(tolower(symb[1]));
    // RDKit✔️✔️: }
    // RDKit✔️✔️: atom = new Atom(PeriodicTable::getTable()->getAtomicNumber(symb));
    let symbol = normalize_xyz_symbol(&atom_line[delims[0]..delims[1]]);
    let element = element_from_rdkit_symbol(&symbol)?;
    Ok((element, pos))
    // END RDKIT CPP FUNCTION
}

fn parse_extra_line(extra_line: &str) -> Result<(), XyzReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseExtraLine
    // RDKit✔️✔️: std::string whitespace{" \t"};
    // RDKit✔️✔️: if (extraLine.find_first_not_of(whitespace) != std::string::npos) {
    // RDKit✔️✔️:   std::ostringstream errout;
    // RDKit✔️✔️:   errout << "More lines than expected" << std::endl;
    // RDKit✔️✔️:   throw FileParseException(errout.str());
    // RDKit✔️✔️: }
    if extra_line.trim_matches([' ', '\t']).is_empty() {
        Ok(())
    } else {
        Err(XyzReadError::AtomSymbol {
            message: "More lines than expected".to_string(),
        })
    }
    // END RDKIT CPP FUNCTION
}

pub fn read_xyz_from_str(block: &str) -> Result<Molecule, XyzReadError> {
    // BEGIN RDKIT CPP FUNCTION MolFromXYZDataStream / MolFromXYZBlock
    // RDKit✔️✔️: std::unique_ptr<RWMol> MolFromXYZBlock(const std::string &xyzBlock) {
    // RDKit✔️✔️:   std::istringstream xyz(xyzBlock);
    // RDKit✔️✔️:   xyz.peek();
    // RDKit✔️✔️:   if (!xyz.eof()) {
    // RDKit✔️✔️:     return MolFromXYZDataStream(xyz);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   }
    if block.is_empty() {
        return Err(XyzReadError::EmptyBlock);
    }

    let mut lines = block.lines();
    // RDKit✔️✔️: std::string num{getLine(inStream)};
    // RDKit✔️✔️: numAtoms = FileParserUtils::toUnsigned(num);
    let atom_count_line = lines.next().ok_or(XyzReadError::EmptyBlock)?;
    let num_atoms = rdkit_to_unsigned(atom_count_line)?;

    // RDKit✔️✔️: std::string comment{getLine(inStream)};
    let comment = lines.next().unwrap_or_default();

    // RDKit✔️✔️: auto mol = std::make_unique<RWMol>();
    // RDKit✔️✔️: if (numAtoms) {
    // RDKit✔️✔️:   Conformer *conf = new Conformer(numAtoms);
    // RDKit✔️✔️:   if (!comment.empty()) {
    // RDKit✔️✔️:     mol->setProp("_FileComments", comment);
    // RDKit✔️✔️:   }
    let mut builder = MoleculeBuilder::new();
    let mut coords = Vec::with_capacity(num_atoms);
    if num_atoms > 0 && !comment.is_empty() {
        builder = builder.with_property("_FileComments", comment);
    }

    // RDKit✔️✔️: for (unsigned int i = 0; i < numAtoms; i++) {
    // RDKit✔️✔️:   if (inStream.eof()) {
    // RDKit✔️✔️:     throw FileParseException("EOF hit while reading atoms");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDGeom::Point3D pos;
    // RDKit✔️✔️:   std::string atomLine{getLine(inStream)};
    // RDKit✔️✔️:   Atom *atom = ParseXYZFileAtomLine(atomLine, pos, i + 2);
    // RDKit✔️✔️:   unsigned int idx = mol->addAtom(atom, false, true);
    // RDKit✔️✔️:   conf->setAtomPos(idx, pos);
    // RDKit✔️✔️: }
    for i in 0..num_atoms {
        let atom_line = lines.next().ok_or(XyzReadError::UnexpectedEof)?;
        let (element, coord) = parse_xyz_atom_line(atom_line, i + 2)?;
        builder.add_atom(AtomSpec::new(element));
        coords.push(coord);
    }
    if num_atoms > 0 {
        builder.add_3d_conformer(coords)?;
    }

    // RDKit✔️✔️: while (!inStream.eof()) {
    // RDKit✔️✔️:   std::string extraLine{getLine(inStream)};
    // RDKit✔️✔️:   ParseExtraLine(extraLine);
    // RDKit✔️✔️: }
    for extra_line in lines {
        parse_extra_line(extra_line)?;
    }

    builder.build().map_err(Into::into)
    // END RDKIT CPP FUNCTION
}

#[cfg(test)]
mod tests {
    use super::read_xyz_from_str;

    #[test]
    fn xyz_reader_parses_atoms_coordinates_comment_and_no_bonds_like_rdkit() {
        let mol = read_xyz_from_str(
            "5\nmethane\nC 0.0 0.0 0.0\nH -0.635 -0.635 0.635\nH -0.635 0.635 -0.635\nH 0.635 -0.635 -0.635\nH 0.635 0.635 0.635\n",
        )
        .expect("xyz parse");

        assert_eq!(mol.num_atoms(), 5);
        assert_eq!(mol.num_bonds(), 0);
        assert_eq!(mol.atomic_numbers(), vec![6, 1, 1, 1, 1]);
        assert_eq!(mol.properties().prop("_FileComments"), Some("methane"));
        assert_eq!(mol.conformers_3d().len(), 1);
        assert_eq!(mol.conformers_3d()[0].coords()[1], [-0.635, -0.635, 0.635]);
    }

    #[test]
    fn xyz_reader_accepts_uppercase_second_symbol_letter_like_rdkit() {
        let mol = read_xyz_from_str("1\n\nCL 1 2 3\n").expect("xyz parse");

        assert_eq!(mol.atomic_numbers(), vec![17]);
        assert_eq!(mol.conformers_3d()[0].coords()[0], [1.0, 2.0, 3.0]);
    }

    #[test]
    fn xyz_reader_rejects_extra_atom_line_fields_like_rdkit() {
        let err = read_xyz_from_str("1\n\nC 1 2 3 extra\n").expect_err("xyz parse should fail");

        assert!(err.to_string().contains("cannot convert '3 extra'"));
    }

    #[test]
    fn xyz_reader_rejects_scientific_notation_like_rdkit_to_double_false() {
        let err = read_xyz_from_str("1\n\nC 1e0 2 3\n").expect_err("xyz parse should fail");

        assert!(err.to_string().contains("cannot convert '1e0'"));
    }

    #[test]
    fn xyz_reader_atom_count_accepts_only_rdkit_unsigned_space_rules() {
        let mol = read_xyz_from_str("   \ncomment\n").expect("space-only count is zero");
        assert_eq!(mol.num_atoms(), 0);

        let mol = read_xyz_from_str("1 2\n\nC 1 2 3\n").expect("from_chars reads prefix");
        assert_eq!(mol.num_atoms(), 1);

        let err = read_xyz_from_str("\t1\n\nC 1 2 3\n").expect_err("tab is not accepted");
        assert!(
            err.to_string()
                .contains("unable to recognize the number of atoms")
        );
    }
}
