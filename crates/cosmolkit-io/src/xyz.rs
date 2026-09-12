//! RDKit-aligned XYZ reader over detached model values.

use std::num::ParseFloatError;

use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Conformer3D, CoordinateBlock, MoleculeProperties,
    TopologyBlock,
};
use cosmolkit_types::Element;

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
    #[error("invalid detached XYZ topology: {0}")]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error("invalid detached XYZ coordinates: {0}")]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error("invalid molecule property: {0}")]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum XyzWriteError {
    #[error("detached XYZ topology is invalid: {0}")]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error("detached XYZ coordinates are invalid: {0}")]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error("XYZ conformer id {id} was not found")]
    ConformerNotFound { id: usize },
}

/// Controls RDKit-aligned XYZ conformer selection and coordinate precision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct XyzWriteParams {
    /// Select a conformer by id. `None` uses the first available conformer.
    pub conformer_id: Option<usize>,
    /// Number of digits after the decimal point.
    pub precision: u32,
}

impl Default for XyzWriteParams {
    fn default() -> Self {
        Self {
            conformer_id: None,
            precision: 6,
        }
    }
}

fn parse_float_error() -> ParseFloatError {
    "x".parse::<f64>().unwrap_err()
}

fn rdkit_to_unsigned(value: &str) -> Result<usize, XyzReadError> {
    // BEGIN RDKIT CPP FUNCTION FileParserUtils::toUnsigned
    // RDKit✔️✔️:   for (size_t i = 0u; i < input.size() && *txt != '\x00'; ++i) {
    // RDKit✔️✔️:     if ((*txt >= '0' && *txt <= '9') || (acceptSpaces && *txt == ' ') ||
    // RDKit✔️✔️:         *txt == '+') {
    // RDKit✔️✔️:       ++txt;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       throw boost::bad_lexical_cast();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   std::from_chars(txt, txt + sz, res);
    // RDKit✔️✔️:   return res;
    let checked = value.split_once('\0').map_or(value, |(prefix, _)| prefix);
    if !checked
        .bytes()
        .all(|byte| byte.is_ascii_digit() || byte == b' ' || byte == b'+')
    {
        return Err(XyzReadError::AtomCount {
            value: value.to_string(),
        });
    }
    let trimmed = checked.trim_start_matches(' ');
    if !checked.is_empty() && trimmed.is_empty() {
        return Ok(0);
    }
    let count = trimmed.bytes().take_while(u8::is_ascii_digit).count();
    if count == 0 {
        return Ok(0);
    }
    let result = Ok(trimmed[..count].parse::<u32>().unwrap_or(0) as usize);
    // END RDKIT CPP FUNCTION
    result
}

fn rdkit_to_double_no_spaces(value: &str, line: usize) -> Result<f64, XyzReadError> {
    // BEGIN RDKIT CPP FUNCTION FileParserUtils::toDouble
    // RDKit✔️✔️:   for (size_t i = 0u; i < input.size() && *txt != '\x00'; ++i) {
    // RDKit✔️✔️:     if ((*txt >= '0' && *txt <= '9') || (acceptSpaces && *txt == ' ') ||
    // RDKit✔️✔️:         *txt == '+' || *txt == '-' || *txt == ',' || *txt == '.') {
    // RDKit✔️✔️:       ++txt;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       throw boost::bad_lexical_cast();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double res = atof(input.data());
    // RDKit✔️✔️:   return res;
    let checked = value.split_once('\0').map_or(value, |(prefix, _)| prefix);
    if !checked
        .bytes()
        .all(|byte| byte.is_ascii_digit() || matches!(byte, b'+' | b'-' | b',' | b'.'))
    {
        return Err(XyzReadError::Coordinate {
            value: value.to_string(),
            line,
            source: value.parse::<f64>().err().unwrap_or_else(parse_float_error),
        });
    }
    let result = checked
        .split_once(',')
        .map_or(checked, |(prefix, _)| prefix)
        .parse::<f64>()
        .map_err(|source| XyzReadError::Coordinate {
            value: value.to_string(),
            line,
            source,
        });
    // END RDKIT CPP FUNCTION
    result
}

fn normalize_symbol(raw: &str) -> String {
    let mut chars = raw.chars().collect::<Vec<_>>();
    if chars.len() == 2 && chars[1].is_ascii_uppercase() {
        chars[1] = chars[1].to_ascii_lowercase();
    }
    chars.into_iter().collect()
}

fn parse_atom_line(line_text: &str, line: usize) -> Result<(Element, [f64; 3]), XyzReadError> {
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
            line_text[prev..]
                .find(|ch: char| ch != ' ' && ch != '\t')
                .map(|index| prev + index)
        } else {
            line_text[prev..]
                .find([' ', '\t'])
                .map(|index| prev + index)
        }
        .ok_or(XyzReadError::MissingCoordinates { line })?;
        prev = *delim;
    }
    delims[7] = line_text
        .rfind(|ch: char| ch != ' ' && ch != '\t')
        .map_or(0, |index| index + 1);
    let symbol = normalize_symbol(&line_text[delims[0]..delims[1]]);
    let element = Element::from_symbol(&symbol).ok_or_else(|| XyzReadError::AtomSymbol {
        message: format!("Element '{symbol}' not found"),
    })?;
    let coord = [
        rdkit_to_double_no_spaces(&line_text[delims[2]..delims[3]], line)?,
        rdkit_to_double_no_spaces(&line_text[delims[4]..delims[5]], line)?,
        rdkit_to_double_no_spaces(&line_text[delims[6]..delims[7]], line)?,
    ];
    Ok((element, coord))
}

/// Read an XYZ block into detached topology, coordinates, and properties.
pub fn read_xyz_detached(
    block: &str,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), XyzReadError> {
    // BEGIN RDKIT CPP FUNCTION MolFromXYZDataStream
    // RDKit✔️✔️:   unsigned int numAtoms = 0;
    // RDKit✔️✔️:   std::string num{getLine(inStream)};
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     numAtoms = FileParserUtils::toUnsigned(num);
    // RDKit✔️✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Unable to recognize the number of atoms: cannot convert '" << num
    // RDKit✔️✔️:            << "' to unsigned int on line 0" << std::endl;
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::string comment{getLine(inStream)};
    // RDKit✔️✔️:   auto mol = std::make_unique<RWMol>();
    // RDKit✔️✔️:   if (numAtoms) {
    // RDKit✔️✔️:     Conformer *conf = new Conformer(numAtoms);
    // RDKit✔️✔️:     if (!comment.empty()) {
    // RDKit✔️✔️:       mol->setProp("_FileComments", comment);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int i = 0; i < numAtoms; i++) {
    // RDKit✔️✔️:       if (inStream.eof()) {
    // RDKit✔️✔️:         throw FileParseException("EOF hit while reading atoms");
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       RDGeom::Point3D pos;
    // RDKit✔️✔️:       std::string atomLine{getLine(inStream)};
    // RDKit✔️✔️:       Atom *atom = ParseXYZFileAtomLine(atomLine, pos, i + 2);
    // RDKit✔️✔️:       unsigned int idx = mol->addAtom(atom, false, true);
    // RDKit✔️✔️:       conf->setAtomPos(idx, pos);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     mol->addConformer(conf);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   while (!inStream.eof()) {
    // RDKit✔️✔️:     std::string extraLine{getLine(inStream)};
    // RDKit✔️✔️:     ParseExtraLine(extraLine);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return mol;
    if block.is_empty() {
        return Err(XyzReadError::EmptyBlock);
    }
    let mut lines = block.lines();
    let count_line = lines.next().ok_or(XyzReadError::EmptyBlock)?;
    let atom_count = rdkit_to_unsigned(count_line)?;
    let comment = lines.next().unwrap_or_default();
    let mut atoms = Vec::with_capacity(atom_count.min(4096));
    let mut coords = Vec::with_capacity(atom_count.min(4096));
    for index in 0..atom_count {
        let line_text = lines.next().ok_or(XyzReadError::UnexpectedEof)?;
        let (element, coord) = parse_atom_line(line_text, index + 2)?;
        atoms.push(Atom::from_spec(AtomId::new(index), AtomSpec::new(element)));
        coords.push(coord);
    }
    for extra in lines {
        if !extra.trim_matches([' ', '\t']).is_empty() {
            return Err(XyzReadError::AtomSymbol {
                message: "More lines than expected".into(),
            });
        }
    }
    let bonds = Vec::new();
    let topology = TopologyBlock {
        adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
        atoms,
        bonds,
        substance_groups: Vec::new(),
        stereo_groups: Vec::new(),
    };
    topology.validate()?;
    let mut properties = MoleculeProperties::default();
    if atom_count > 0 && !comment.is_empty() {
        properties = properties.with_prop("_FileComments", comment)?;
    }
    let coordinates = CoordinateBlock {
        conformers_2d: Vec::new(),
        conformers_3d: if atom_count > 0 {
            vec![Conformer3D::new(0, coords, true)]
        } else {
            Vec::new()
        },
        source_coordinate_dim: (atom_count > 0)
            .then_some(cosmolkit_model::CoordinateDimension::ThreeD),
    };
    coordinates.validate_for_atom_count(atom_count)?;
    let result = Ok((topology, coordinates, properties));
    // END RDKIT CPP FUNCTION
    result
}

/// Write detached model values as an XYZ block using RDKit's defaults.
pub fn write_xyz_detached(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> Result<String, XyzWriteError> {
    write_xyz_detached_with_params(topology, coordinates, properties, XyzWriteParams::default())
}

/// Write detached model values as an XYZ block with explicit writer options.
pub fn write_xyz_detached_with_params(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    params: XyzWriteParams,
) -> Result<String, XyzWriteError> {
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;

    // BEGIN RDKIT CPP FUNCTION MolToXYZBlock
    // RDKit✔️✔️: if (!mol.getNumConformers()) {
    // RDKit✔️✔️:   BOOST_LOG(rdErrorLog)
    // RDKit✔️✔️:       << "Cannot write molecules with no conformers to XYZ block\n";
    // RDKit✔️✔️:   return "";
    // RDKit✔️✔️: }
    if coordinates.conformers_3d.is_empty() && coordinates.conformers_2d.is_empty() {
        return Ok(String::new());
    }

    // RDKit✔️✔️: const auto &conf = mol.getConformer(confId);
    // RDKit✔️✔️: const unsigned int nAtoms = mol.getNumAtoms();
    enum Points<'a> {
        TwoD(&'a [[f64; 2]]),
        ThreeD(&'a [[f64; 3]]),
    }
    let points = match params.conformer_id {
        Some(id) => coordinates
            .conformers_3d
            .iter()
            .find(|conformer| conformer.id() == id)
            .map(|conformer| Points::ThreeD(conformer.coordinates()))
            .or_else(|| {
                coordinates
                    .conformers_2d
                    .iter()
                    .find(|conformer| conformer.id() == id)
                    .map(|conformer| Points::TwoD(conformer.coordinates()))
            })
            .ok_or(XyzWriteError::ConformerNotFound { id })?,
        None => coordinates
            .conformers_3d
            .first()
            .map(|conformer| Points::ThreeD(conformer.coordinates()))
            .or_else(|| {
                coordinates
                    .conformers_2d
                    .first()
                    .map(|conformer| Points::TwoD(conformer.coordinates()))
            })
            .expect("the empty conformer case returned above"),
    };

    // RDKit✔️✔️: std::stringstream ss;
    // RDKit✔️✔️: ss << nAtoms << '\n';
    // RDKit✔️✔️: unsigned fieldWidth = 5 + precision;
    // RDKit✔️✔️: std::stringstream formatString;
    // RDKit✔️✔️: formatString << "%-3s %" << fieldWidth << "." << precision << "f %"
    // RDKit✔️✔️:              << fieldWidth << "." << precision << "f %" << fieldWidth << "."
    // RDKit✔️✔️:              << precision << "f\n";
    let precision = params.precision as usize;
    let field_width = precision.saturating_add(5);
    let mut output = format!("{}\n", topology.atoms.len());

    // RDKit✔️✔️: std::string name;
    // RDKit✔️✔️: if (mol.getPropIfPresent(common_properties::_Name, name)) {
    // RDKit✔️✔️:   ss << name.substr(0, name.find_first_of('\n'));
    // RDKit✔️✔️: }
    // RDKit✔️✔️: ss << '\n';
    if let Some(name) = properties.name() {
        output.push_str(name.split_once('\n').map_or(name, |(first, _)| first));
    }
    output.push('\n');

    // RDKit✔️✔️: for (unsigned int i = 0; i < nAtoms; i++) {
    // RDKit✔️✔️:   const auto &symbol = mol.getAtomWithIdx(i)->getSymbol();
    // RDKit✔️✔️:   const auto &pos = conf.getAtomPos(i);
    // RDKit✔️✔️:   ss << boost::format{formatString.str()} % symbol % pos.x % pos.y % pos.z;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return ss.str();
    for (index, atom) in topology.atoms.iter().enumerate() {
        let point = match points {
            Points::TwoD(values) => [values[index][0], values[index][1], 0.0],
            Points::ThreeD(values) => values[index],
        };
        output.push_str(&format!(
            "{:<3} {:>width$.precision$} {:>width$.precision$} {:>width$.precision$}\n",
            atom.element().symbol(),
            point[0],
            point[1],
            point[2],
            width = field_width,
            precision = precision,
        ));
    }
    // END RDKIT CPP FUNCTION
    Ok(output)
}

#[cfg(test)]
mod tests {
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomId, AtomSpec, Conformer2D, Conformer3D, CoordinateBlock,
        MoleculeProperties, TopologyBlock,
    };
    use cosmolkit_types::Element;

    use super::{
        XyzWriteParams, read_xyz_detached, write_xyz_detached, write_xyz_detached_with_params,
    };

    #[test]
    fn parses_detached_xyz_without_bonds() {
        let (topology, coordinates, properties) =
            read_xyz_detached("2\nwater\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\n").expect("xyz parse");
        assert_eq!(topology.atoms.len(), 2);
        assert!(topology.bonds.is_empty());
        assert_eq!(
            coordinates.conformers_3d[0].coordinates()[1],
            [0.0, 0.0, 1.0]
        );
        assert_eq!(properties.prop("_FileComments"), Some("water"));
    }

    #[test]
    fn preserves_rdkit_xyz_edge_validation() {
        let (topology, coordinates, _) =
            read_xyz_detached("1\n\nCL 1 2 3\n").expect("uppercase symbol");
        assert_eq!(topology.atoms[0].element().atomic_number(), 17);
        assert_eq!(
            coordinates.conformers_3d[0].coordinates()[0],
            [1.0, 2.0, 3.0]
        );

        let err = read_xyz_detached("1\n\nC 1 2 3 extra\n").expect_err("extra field");
        assert!(err.to_string().contains("cannot convert '3 extra'"));

        let err = read_xyz_detached("1\n\nC 1e0 2 3\n").expect_err("scientific notation");
        assert!(err.to_string().contains("cannot convert '1e0'"));
    }

    #[test]
    fn preserves_rdkit_xyz_count_rules_and_truncation() {
        let (topology, coordinates, _) = read_xyz_detached("   \ncomment\n").expect("zero count");
        assert!(topology.atoms.is_empty());
        assert!(coordinates.conformers_3d.is_empty());

        let (topology, _, _) = read_xyz_detached("1 2\n\nC 1 2 3\n").expect("prefix count");
        assert_eq!(topology.atoms.len(), 1);

        let err = read_xyz_detached("\t1\n\nC 1 2 3\n").expect_err("tab count");
        assert!(
            err.to_string()
                .contains("unable to recognize the number of atoms")
        );

        let err = read_xyz_detached("4294967295\ncomment\n").expect_err("truncated count");
        assert_eq!(err, super::XyzReadError::UnexpectedEof);
    }

    #[test]
    fn detached_xyz_writer_round_trips_coordinates_and_comment() {
        let input = "2\nwater\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\n";
        let (topology, coordinates, properties) = read_xyz_detached(input).expect("read");
        let output = write_xyz_detached(&topology, &coordinates, &properties).expect("write");
        let (roundtrip, coords, props) = read_xyz_detached(&output).expect("roundtrip");
        assert_eq!(roundtrip.atoms.len(), 2);
        assert_eq!(coords.conformers_3d[0].coordinates()[1], [0.0, 0.0, 1.0]);
        assert_eq!(props.prop("_FileComments"), None);
    }

    fn methane_parts() -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
        let elements = [Element::C, Element::H, Element::H, Element::H, Element::H];
        let atoms = elements
            .into_iter()
            .enumerate()
            .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
            .collect::<Vec<_>>();
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &[]),
            atoms,
            ..TopologyBlock::default()
        };
        let coordinates = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [-0.635, -0.635, 0.635],
                    [-0.635, 0.635, -0.635],
                    [0.635, -0.635, -0.635],
                    [0.635, 0.635, 0.635],
                ],
                true,
            )],
            ..CoordinateBlock::default()
        };
        let properties = MoleculeProperties::default()
            .with_name("methane\nthis part should not be output")
            .with_prop("_FileComments", "ignored comment")
            .expect("the internal file-comments property key is non-empty");
        (topology, coordinates, properties)
    }

    #[test]
    fn writer_matches_rdkit_default_format_and_name_rule() {
        let (topology, coordinates, properties) = methane_parts();
        let output = write_xyz_detached(&topology, &coordinates, &properties).expect("write");
        assert_eq!(
            output,
            "5\nmethane\nC      0.000000    0.000000    0.000000\nH     -0.635000   -0.635000    0.635000\nH     -0.635000    0.635000   -0.635000\nH      0.635000   -0.635000   -0.635000\nH      0.635000    0.635000    0.635000\n"
        );
    }

    #[test]
    fn writer_selects_conformer_and_precision_like_rdkit() {
        let (topology, mut coordinates, properties) = methane_parts();
        coordinates
            .conformers_3d
            .push(Conformer3D::new(7, vec![[1.25, -2.5, 3.75]; 5], true));
        let output = write_xyz_detached_with_params(
            &topology,
            &coordinates,
            &properties,
            XyzWriteParams {
                conformer_id: Some(7),
                precision: 2,
            },
        )
        .expect("selected write");
        assert!(output.contains("C      1.25   -2.50    3.75\n"));

        let error = write_xyz_detached_with_params(
            &topology,
            &coordinates,
            &properties,
            XyzWriteParams {
                conformer_id: Some(8),
                precision: 6,
            },
        )
        .expect_err("missing conformer");
        assert!(error.to_string().contains("conformer id 8"));
    }

    #[test]
    fn writer_returns_empty_block_without_conformers_and_supports_2d() {
        let (topology, _, properties) = methane_parts();
        assert_eq!(
            write_xyz_detached(&topology, &CoordinateBlock::default(), &properties)
                .expect("empty output"),
            ""
        );

        let coordinates = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(3, vec![[1.0, 2.0]; 5])],
            ..CoordinateBlock::default()
        };
        let output = write_xyz_detached(&topology, &coordinates, &properties).expect("2D write");
        assert!(output.contains("C      1.000000    2.000000    0.000000\n"));
    }

    #[test]
    fn reader_reproduces_unsigned_overflow_as_zero() {
        let (topology, coordinates, properties) =
            read_xyz_detached("4294967296\ncomment\n").expect("overflow count");
        assert!(topology.atoms.is_empty());
        assert!(coordinates.conformers_2d.is_empty());
        assert!(coordinates.conformers_3d.is_empty());
        assert_eq!(properties.prop("_FileComments"), None);
    }
}
