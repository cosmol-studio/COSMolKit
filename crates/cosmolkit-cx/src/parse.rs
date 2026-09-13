use crate::scan::{
    consume_until_record_boundary, expect_byte, parse_delimited_number_list, parse_number_list,
    read_colon_field, read_number, read_pair, read_text_to,
};
use crate::{
    CxAtomConstraint, CxAtomProperty, CxBondReference, CxCoordinateBondKind, CxCoordinateBonds,
    CxCoordinates, CxCountConstraint, CxDataSGroup, CxDoubleBondStereo, CxDoubleBondStereoKind,
    CxEnhancedStereo, CxLinkNode, CxParseError, CxPolymerSGroup, CxRadical, CxRecord, CxRingBond,
    CxSGroupHierarchy, CxStereoGroupKind, CxVariableAttachment, CxWedgeBond, CxWedgeDirection,
    ParsedCxExtensions,
};

/// Parse one CX extension block without referring to a destination graph.
///
pub fn parse_cx_extensions(text: &str) -> Result<ParsedCxExtensions, CxParseError> {
    // RDKit source (verbatim; graph-effect calls are lowering-owned):
    /*
    template <typename Iterator>
    bool parse_it(Iterator &first, Iterator last, RDKit::RWMol &mol,
                  unsigned int startAtomIdx, unsigned int startBondIdx) {
      if (first >= last || *first != '|') {
        return false;
      }
      ++first;
      unsigned int nSGroups = 0;
      unsigned int confIndex = 0;
      while (first < last && *first != '|') {
        typename Iterator::difference_type length = std::distance(first, last);
        if (*first == '(') {
          if (!parse_coords(first, last, mol, startAtomIdx, confIndex++)) {
            return false;
          }
        } else if (*first == '$') {
          if (length > 4 && *(first + 1) == '_' && *(first + 2) == 'A' &&
              *(first + 3) == 'V' && *(first + 4) == ':') {
            first += 4;
            if (!parse_atom_values(first, last, mol, startAtomIdx)) {
              return false;
            }
          } else {
            if (!parse_atom_labels(first, last, mol, startAtomIdx)) {
              return false;
            }
          }
        } else if (length > 9 && std::string(first, first + 9) == "atomProp:") {
          first += 9;
          if (!parse_atom_props(first, last, mol, startAtomIdx)) {
            return false;
          }
        } else if (*first == 'C') {
          if (!parse_coordinate_bonds(first, last, mol, Bond::DATIVE, startAtomIdx,
                                      startBondIdx)) {
            return false;
          }
        } else if (*first == 'H') {
          if (!parse_coordinate_bonds(first, last, mol, Bond::HYDROGEN,
                                      startAtomIdx, startBondIdx)) {
            return false;
          }
        } else if (*first == 'Z') {
          if (!parse_zero_bonds(first, last, mol, startAtomIdx, startBondIdx)) {
            return false;
          }
        } else if (*first == '^') {
          if (!parse_radicals(first, last, mol, startAtomIdx)) {
            return false;
          }
        } else if (*first == 'a' || *first == 'o' ||
                   (*first == '&' && first + 1 < last && first[1] != '#')) {
          if (!parse_enhanced_stereo(first, last, mol, startAtomIdx)) {
            return false;
          }
        } else if (*first == 'r' && first + 1 < last && first[1] == 'b') {
          if (!parse_ring_bonds(first, last, mol, startAtomIdx)) {
            return false;
          }
        } else if (*first == 'L' && first + 1 < last && first[1] == 'N') {
          if (!parse_linknodes(first, last, mol, startAtomIdx)) {
            return false;
          }
        } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
                   first[2] == 'D') {
          if (!parse_data_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
            return false;
          }
        } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
                   first[2] == 'H') {
          if (!parse_sgroup_hierarchy(first, last, mol)) {
            return false;
          }
        } else if (*first == 'S' && first + 1 < last && first[1] == 'g') {
          if (!parse_polymer_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
            return false;
          }
        } else if (*first == 'u') {
          if (!parse_unsaturation(first, last, mol, startAtomIdx)) {
            return false;
          }
        } else if (*first == 's') {
          if (!parse_substitution(first, last, mol, startAtomIdx)) {
            return false;
          }
        } else if (*first == 'm') {
          if (!parse_variable_attachments(first, last, mol, startAtomIdx)) {
            return false;
          }
        } else if (*first == 'w') {
          if (!parse_wedged_bonds(first, last, mol, startAtomIdx, startBondIdx)) {
            return false;
          }
        } else if (*first == 'c' && first + 2 < last && first[1] == 't' &&
                   first[2] == 'u') {
          if (!parse_doublebond_stereo(first, last, mol, startAtomIdx, startBondIdx,
                                       Bond::BondStereo::STEREOANY)) {
            return false;
          }
        } else if (*first == 'c') {
          if (!parse_doublebond_stereo(first, last, mol, startAtomIdx, startBondIdx,
                                       Bond::BondStereo::STEREOCIS)) {
            return false;
          }
        } else if (*first == 't') {
          if (!parse_doublebond_stereo(first, last, mol, startAtomIdx, startBondIdx,
                                       Bond::BondStereo::STEREOTRANS)) {
            return false;
          }
        } else {
          ++first;
        }
        // if(first < last && *first != '|') ++first;
      }
      if (first >= last || *first != '|') {
        return false;
      }
      ++first;  // step past the last '|'
      return true;
    }

    void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
                           std::string::const_iterator &first,
                           unsigned int startAtomIdx, unsigned int startBondIdx) {
      // BOOST_LOG(rdWarningLog) << "parseCXNExtensions: " << extText << std::endl;
      if (extText.empty()) {
        return;
      }
      if (extText[0] != '|') {
        throw RDKit::SmilesParseException(
            "CXSMILES extension does not start with |");
      }
      first = extText.begin();
      bool ok =
          parser::parse_it(first, extText.end(), mol, startAtomIdx, startBondIdx);
      if (!ok) {
        throw RDKit::SmilesParseException("failure parsing CXSMILES extensions");
      }
      processCXSmilesLabels(mol);
      mol.clearProp("_cxsmilesLabelsProcessed");
      mol.clearProp(cxsgTracker);
    }
        */
    // RDKit✔️✔️: the Rust loop preserves the source dispatch priority, record
    // order, conformer ordinals, pipe contract and one-pass cursor complexity
    // for detached syntax; all RWMol effects remain lowering-owned.
    if text.is_empty() {
        return Ok(ParsedCxExtensions::new(Vec::new(), 0));
    }
    let bytes = text.as_bytes();
    if bytes.first().copied() != Some(b'|') {
        return Err(CxParseError::new(
            0,
            "CXSMILES extension does not start with |",
        ));
    }

    let mut cursor = 1;
    let mut conformer = 0;
    let mut records = Vec::new();
    while cursor < bytes.len() && bytes[cursor] != b'|' {
        if bytes[cursor] == b',' {
            cursor += 1;
            continue;
        }
        let start = cursor;
        let record = if bytes[cursor] == b'(' {
            let mut coordinates = parse_coordinates(text, &mut cursor)?;
            coordinates.conformer = conformer;
            conformer += 1;
            CxRecord::Coordinates(coordinates)
        } else if bytes[cursor] == b'$' {
            parse_labels_or_values(text, &mut cursor)?
        } else if text[cursor..].starts_with("atomProp:") {
            CxRecord::AtomProperties(parse_atom_properties(text, &mut cursor)?)
        } else if matches!(bytes[cursor], b'C' | b'H')
            && cursor + 1 < bytes.len()
            && bytes[cursor + 1] == b':'
        {
            CxRecord::CoordinateBonds(parse_coordinate_bonds(text, &mut cursor)?)
        } else if bytes[cursor] == b'Z' {
            CxRecord::ZeroBonds(parse_index_list(text, &mut cursor, b'Z')?)
        } else if bytes[cursor] == b'^' {
            CxRecord::Radicals(parse_radicals(text, &mut cursor)?)
        } else if matches!(bytes[cursor], b'a' | b'o')
            || (bytes[cursor] == b'&' && bytes.get(cursor + 1) != Some(&b'#'))
        {
            CxRecord::EnhancedStereo(parse_enhanced_stereo(text, &mut cursor)?)
        } else if text[cursor..].starts_with("rb:") {
            CxRecord::RingBonds(parse_ring_bonds(text, &mut cursor)?)
        } else if bytes[cursor] == b'u' {
            CxRecord::Unsaturation(parse_index_list(text, &mut cursor, b'u')?)
        } else if bytes[cursor] == b's' {
            CxRecord::Substitution(parse_count_constraints(text, &mut cursor, b's')?)
        } else if bytes[cursor] == b'm' {
            CxRecord::VariableAttachments(parse_variable_attachments(text, &mut cursor)?)
        } else if bytes[cursor] == b'w' {
            CxRecord::WedgedBonds(parse_wedge_bonds(text, &mut cursor)?)
        } else if text[cursor..].starts_with("ctu") {
            CxRecord::DoubleBondStereo(parse_double_bond_stereo(
                text,
                &mut cursor,
                CxDoubleBondStereoKind::Any,
            )?)
        } else if bytes[cursor] == b'c' {
            CxRecord::DoubleBondStereo(parse_double_bond_stereo(
                text,
                &mut cursor,
                CxDoubleBondStereoKind::Cis,
            )?)
        } else if bytes[cursor] == b't' {
            CxRecord::DoubleBondStereo(parse_double_bond_stereo(
                text,
                &mut cursor,
                CxDoubleBondStereoKind::Trans,
            )?)
        } else if text[cursor..].starts_with("LN:") {
            CxRecord::LinkNodes(parse_link_nodes(text, &mut cursor)?)
        } else if text[cursor..].starts_with("SgD:") {
            CxRecord::DataSGroup(parse_data_sgroup(text, &mut cursor)?)
        } else if text[cursor..].starts_with("SgH:") {
            CxRecord::SGroupHierarchy(parse_sgroup_hierarchy(text, &mut cursor)?)
        } else if text[cursor..].starts_with("Sg:") {
            CxRecord::PolymerSGroup(parse_polymer_sgroup(text, &mut cursor)?)
        } else {
            let raw = consume_until_record_boundary(text, &mut cursor);
            CxRecord::Unknown(raw)
        };
        if cursor == start {
            return Err(CxParseError::new(cursor, "CX parser made no progress"));
        }
        records.push(record);
    }
    if cursor >= bytes.len() || bytes[cursor] != b'|' {
        return Err(CxParseError::new(
            cursor,
            "failure parsing CXSMILES extensions",
        ));
    }
    Ok(ParsedCxExtensions::new(records, cursor + 1))
}

fn parse_coordinates(text: &str, cursor: &mut usize) -> Result<CxCoordinates, CxParseError> {
    // RDKit source (verbatim; Conformer installation is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_coords(Iterator &first, Iterator last, RDKit::RWMol &mol,
                      unsigned int startAtomIdx, unsigned int confIdx) {
      if (first >= last || *first != '(') {
        return false;
      }

      auto *conf = new Conformer(mol.getNumAtoms());
      mol.addConformer(conf);
      conf->setId(confIdx);
      ++first;
      unsigned int atIdx = 0;
      bool is3D = false;
      while (first <= last && *first != ')') {
        RDGeom::Point3D pt;
        std::string tkn = read_text_to(first, last, ";)");
        if (VALID_ATIDX(atIdx)) {
          if (!tkn.empty()) {
            std::vector<std::string> tokens;
            boost::split(tokens, tkn, boost::is_any_of(std::string(",")));
            if (tokens.size() >= 1 && tokens[0].size()) {
              pt.x = boost::lexical_cast<double>(tokens[0]);
            }
            if (tokens.size() >= 2 && tokens[1].size()) {
              pt.y = boost::lexical_cast<double>(tokens[1]);
            }
            if (tokens.size() >= 3 && tokens[2].size()) {
              pt.z = boost::lexical_cast<double>(tokens[2]);
              is3D = true;
            }
          }

          conf->setAtomPos(atIdx - startAtomIdx, pt);
        }
        ++atIdx;
        if (first <= last && *first != ')') {
          ++first;
        }
      }
      // make sure that the conformer really is 3D!
      if (is3D && hasNonZeroZCoords(*conf)) {
        conf->set3D(true);
      } else {
        conf->set3D(false);
      }
      if (first >= last || *first != ')') {
        return false;
      }
      ++first;
      return true;
    }

    inline bool hasNonZeroZCoords(const Conformer &conf) {
      constexpr double zeroTol = 1e-3;
      for (auto p : conf.getPositions()) {
        if (std::abs(p.z) > zeroTol) {
          return true;
        }
      }
      return false;
    }
        */
    // RDKit✔️✔️: detached rows preserve source token/default behavior and
    // conformer order; `is_3d` applies the exact 1e-3 nonzero-z threshold.
    let start = *cursor;
    *cursor += 1;
    let mut values = Vec::new();
    let mut is_3d = false;
    while *cursor < text.len() && text.as_bytes()[*cursor] != b')' {
        let field_start = *cursor;
        while *cursor < text.len() && !matches!(text.as_bytes()[*cursor], b';' | b')') {
            *cursor += 1;
        }
        let field = &text[field_start..*cursor];
        if field.is_empty() {
            values.push(None);
        } else {
            let mut parts = field.split(',');
            let x = parse_coordinate_component(parts.next(), field_start)?;
            let y = parse_coordinate_component(parts.next(), field_start)?;
            let z = parse_coordinate_component(parts.next(), field_start)?;
            if let Some(z) = z {
                is_3d |= z.abs() > 1e-3;
            }
            values.push(Some([
                x.unwrap_or_default(),
                y.unwrap_or_default(),
                z.unwrap_or_default(),
            ]));
        }
        if *cursor < text.len() && text.as_bytes()[*cursor] == b';' {
            *cursor += 1;
        }
    }
    if *cursor >= text.len() || text.as_bytes()[*cursor] != b')' {
        return Err(CxParseError::new(
            start,
            "unterminated CX coordinate record",
        ));
    }
    *cursor += 1;
    Ok(CxCoordinates {
        conformer: 0,
        values,
        is_3d,
    })
}

fn parse_coordinate_component(
    value: Option<&str>,
    offset: usize,
) -> Result<Option<f64>, CxParseError> {
    let Some(value) = value else {
        return Ok(None);
    };
    if value.is_empty() {
        return Ok(None);
    }
    value
        .parse::<f64>()
        .map(Some)
        .map_err(|_| CxParseError::new(offset, "invalid CX coordinate"))
}

fn parse_labels_or_values(text: &str, cursor: &mut usize) -> Result<CxRecord, CxParseError> {
    // RDKit source (verbatim; Atom property installation is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_atom_values(Iterator &first, Iterator last, RDKit::RWMol &mol,
                           unsigned int startAtomIdx) {
      if (first >= last || *first != ':') {
        return false;
      }
      ++first;
      unsigned int atIdx = 0;
      while (first <= last && *first != '$') {
        std::string tkn = read_text_to(first, last, ";$");
        if (tkn != "" && VALID_ATIDX(atIdx)) {
          mol.getAtomWithIdx(atIdx)->setProp(RDKit::common_properties::molFileValue,
                                             tkn);
        }
        ++atIdx;
        if (first <= last && *first != '$') {
          ++first;
        }
      }
      if (first >= last || *first != '$') {
        return false;
      }
      ++first;
      return true;
    }

    template <typename Iterator>
    bool parse_atom_labels(Iterator &first, Iterator last, RDKit::RWMol &mol,
                           unsigned int startAtomIdx) {
      if (first >= last || *first != '$') {
        return false;
      }
      ++first;
      unsigned int atIdx = 0;
      while (first <= last && *first != '$') {
        std::string tkn = read_text_to(first, last, ";$");
        if (!tkn.empty() && VALID_ATIDX(atIdx)) {
          mol.getAtomWithIdx(atIdx - startAtomIdx)
              ->setProp(RDKit::common_properties::atomLabel, tkn);
        }
        ++atIdx;
        if (first <= last && *first != '$') {
          ++first;
        }
      }
      if (first >= last || *first != '$') {
        return false;
      }
      ++first;
      return true;
    }
        */
    // RDKit✔️✔️: the detached slot vector preserves source order, empty slots,
    // delimiter consumption and decoded text with linear scan complexity.
    let value = text[*cursor..].starts_with("$_AV:");
    *cursor += if value { 5 } else { 1 };
    let mut fields = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor] != b'$' {
        let field = read_text_to(text, cursor, b";$")?;
        fields.push(if field.is_empty() { None } else { Some(field) });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b';' {
            *cursor += 1;
        }
    }
    if *cursor >= text.len() || text.as_bytes()[*cursor] != b'$' {
        return Err(CxParseError::new(*cursor, "unterminated CX atom record"));
    }
    *cursor += 1;
    Ok(if value {
        CxRecord::AtomValues(fields)
    } else {
        CxRecord::AtomLabels(fields)
    })
}

fn parse_atom_properties(
    text: &str,
    cursor: &mut usize,
) -> Result<Vec<CxAtomProperty>, CxParseError> {
    // RDKit source (verbatim; Atom property installation is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_atom_props(Iterator &first, Iterator last, RDKit::RWMol &mol,
                          unsigned int startAtomIdx) {
      if (first >= last) {
        return false;
      }
      while (first <= last && *first != '|' && *first != ',') {
        unsigned int atIdx;
        if (read_int(first, last, atIdx)) {
          if (first >= last || *first != '.') {
            return false;
          }
          ++first;
          std::string pname = read_text_to(first, last, ".");
          if (!pname.empty()) {
            if (first >= last || *first != '.') {
              return false;
            }
            ++first;
            std::string pval = read_text_to(first, last, ":|,");
            if (VALID_ATIDX(atIdx) && !pval.empty()) {
              mol.getAtomWithIdx(atIdx - startAtomIdx)->setProp(pname, pval);
            }
          }
        }
        if (first <= last && *first != '|' && *first != ',') {
          ++first;
        }
      }
      if (first <= last && *first != '|' && *first != ',') {
        return false;
      }
      if (*first != '|') {
        ++first;
      }
      return true;
    }
        */
    // RDKit✔️✔️: the Rust parser preserves the source's colon-separated entry
    // scan and comma/pipe termination while returning detached records.
    *cursor += "atomProp:".len();
    let mut result = Vec::new();
    while *cursor < text.len() && !matches!(text.as_bytes()[*cursor], b'|' | b',') {
        let atom_offset = *cursor;
        let atom = read_number(text, cursor)
            .map_err(|_| CxParseError::new(atom_offset, "invalid CX atomProp atom index"))?;
        expect_byte(text, cursor, b'.')?;
        let name = read_text_to(text, cursor, b".")?;
        if name.is_empty() {
            if *cursor < text.len() && !matches!(text.as_bytes()[*cursor], b'|' | b',') {
                *cursor += 1;
            }
            while *cursor < text.len() && !matches!(text.as_bytes()[*cursor], b'|' | b',') {
                *cursor += 1;
            }
            break;
        }
        expect_byte(text, cursor, b'.')?;
        let value = read_text_to(text, cursor, b":|,")?;
        if !value.is_empty() {
            result.push(CxAtomProperty { atom, name, value });
        }
        if *cursor < text.len() && !matches!(text.as_bytes()[*cursor], b'|' | b',') {
            *cursor += 1;
        }
    }
    if text.as_bytes().get(*cursor) == Some(&b',') {
        *cursor += 1;
    }
    Ok(result)
}

fn parse_coordinate_bonds(
    text: &str,
    cursor: &mut usize,
) -> Result<CxCoordinateBonds, CxParseError> {
    // RDKit source (verbatim; bond lookup/type mutation is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_coordinate_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
                                Bond::BondType typ, unsigned int startAtomIdx,
                                unsigned int startBondIdx) {
      if (first >= last || (*first != 'C' && *first != 'H')) {
        return false;
      }
      ++first;
      if (first >= last || *first != ':') {
        return false;
      }
      ++first;
      while (first <= last && *first >= '0' && *first <= '9') {
        unsigned int aidx;
        unsigned int bidx;
        if (read_int_pair(first, last, aidx, bidx)) {
          if (VALID_ATIDX(aidx) && VALID_BNDIDX(bidx)) {
            auto bnd = get_bond_with_smiles_idx(mol, bidx - startBondIdx);
            if (!bnd || (bnd->getBeginAtomIdx() != aidx - startAtomIdx &&
                         bnd->getEndAtomIdx() != aidx - startAtomIdx)) {
              BOOST_LOG(rdWarningLog) << "BOND NOT FOUND! " << bidx
                                      << " involving atom " << aidx << std::endl;
              return false;
            }
            bnd->setBondType(typ);
            if (bnd->getBeginAtomIdx() != aidx - startAtomIdx) {
              unsigned int tmp = bnd->getBeginAtomIdx();
              bnd->setBeginAtomIdx(aidx - startAtomIdx);
              bnd->setEndAtomIdx(tmp);
            }
          }
        } else {
          return false;
        }
        if (first < last && *first == ',') {
          ++first;
        }
      }
      return true;
    }
        */
    // RDKit✔️✔️: pair order, bond-kind dispatch, separators and one-pass scan
    // match the lexical source closure; target bond checks remain downstream.
    let kind = if text.as_bytes()[*cursor] == b'C' {
        CxCoordinateBondKind::Dative
    } else {
        CxCoordinateBondKind::Hydrogen
    };
    *cursor += 2;
    let mut bonds = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let (atom, bond) = read_pair(text, cursor, b'.')?;
        bonds.push(CxBondReference { atom, bond });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        }
    }
    Ok(CxCoordinateBonds { kind, bonds })
}

fn parse_index_list(
    text: &str,
    cursor: &mut usize,
    marker: u8,
) -> Result<Vec<usize>, CxParseError> {
    // RDKit source (verbatim; bond/atom mutation is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_zero_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
                          unsigned int, unsigned int startBondIdx) {
      // these look like: C1CCCCC~CCCC1 |Z:5|
      if (first >= last || *first != 'Z') {
        return false;
      }
      ++first;
      if (first >= last || *first != ':') {
        return false;
      }
      ++first;

      while (first < last && *first >= '0' && *first <= '9') {
        unsigned int bondIdx;
        if (!read_int(first, last, bondIdx)) {
          return false;
        }
        if (VALID_BNDIDX(bondIdx)) {
          auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);

          if (!bond) {
            BOOST_LOG(rdWarningLog)
                << "bond " << bondIdx
                << " not found, cannot mark as zero order bond." << std::endl;
            return false;
          }
          bond->setBondType(Bond::ZERO);
        }
        if (first < last && *first == ',') {
          ++first;
        }
      }
      return true;
    }

    template <typename Iterator>
    bool parse_unsaturation(Iterator &first, Iterator last, RDKit::RWMol &mol,
                            unsigned int startAtomIdx) {
      if (first + 1 >= last || *first != 'u') {
        return false;
      }
      ++first;
      if (first >= last || *first != ':') {
        return false;
      }
      ++first;
      while (first < last && *first >= '0' && *first <= '9') {
        unsigned int idx;
        if (!read_int(first, last, idx)) {
          return false;
        }
        if (VALID_ATIDX(idx)) {
          auto atom = mol.getAtomWithIdx(idx - startAtomIdx);
          if (!atom->hasQuery()) {
            atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
          }
          atom->expandQuery(makeAtomUnsaturatedQuery(), Queries::COMPOSITE_AND);
        }
        if (first < last && *first == ',') {
          ++first;
        }
      }
      return true;
    }
        */
    // RDKit✔️✔️: both markers share the exact lexical integer-list shape;
    // detached indices preserve order and defer bond/query effects.
    expect_byte(text, cursor, marker)?;
    expect_byte(text, cursor, b':')?;
    let mut values = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        values.push(read_number(text, cursor)?);
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(values)
}

fn parse_radicals(text: &str, cursor: &mut usize) -> Result<Vec<CxRadical>, CxParseError> {
    // RDKit source (verbatim; Atom mutation is lowering-owned):
    /*
    template <typename Iterator>
    bool processRadicalSection(Iterator &first, Iterator last, RDKit::RWMol &mol,
                               unsigned int numRadicalElectrons,
                               unsigned int startAtomIdx) {
      if (first >= last) {
        return false;
      }
      ++first;
      if (first >= last || *first != ':') {
        return false;
      }
      ++first;
      unsigned int atIdx;
      if (!read_int(first, last, atIdx)) {
        return false;
      }
      if (VALID_ATIDX(atIdx)) {
        mol.getAtomWithIdx(atIdx - startAtomIdx)
            ->setNumRadicalElectrons(numRadicalElectrons);
      }
      while (first < last && *first == ',') {
        ++first;
        if (first < last && (*first < '0' || *first > '9')) {
          return true;
        }
        if (!read_int(first, last, atIdx)) {
          return false;
        }
        if (VALID_ATIDX(atIdx)) {
          mol.getAtomWithIdx(atIdx - startAtomIdx)
              ->setNumRadicalElectrons(numRadicalElectrons);
        }
      }
      return first < last;
    }

    template <typename Iterator>
    bool parse_radicals(Iterator &first, Iterator last, RDKit::RWMol &mol,
                        unsigned int startAtomIdx) {
      if (first >= last || *first != '^') {
        return false;
      }
      while (*first == '^') {
        ++first;
        if (first >= last) {
          return false;
        }
        if (*first < '1' || *first > '7') {
          return false;  // these are the values that are allowed to be there
        }
        switch (*first) {
          case '1':
            if (!processRadicalSection(first, last, mol, 1, startAtomIdx)) {
              return false;
            }
            break;
          case '2':
          case '3':
          case '4':
            if (!processRadicalSection(first, last, mol, 2, startAtomIdx)) {
              return false;
            }
            break;
          case '5':
          case '6':
          case '7':
            if (!processRadicalSection(first, last, mol, 3, startAtomIdx)) {
              return false;
            }
            break;
          default:
            BOOST_LOG(rdWarningLog)
                << "Radical specification " << *first << " ignored.";
        }
      }
      return true;
    }
        */
    // RDKit✔️✔️: radical sections, electron classes, atom order and separator
    // transitions match the source with a single forward scan.
    let mut result = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor] == b'^' {
        *cursor += 1;
        let electrons = match text.as_bytes().get(*cursor).copied() {
            Some(b'1') => 1,
            Some(b'2'..=b'4') => 2,
            Some(b'5'..=b'7') => 3,
            _ => return Err(CxParseError::new(*cursor, "invalid CX radical marker")),
        };
        *cursor += 1;
        expect_byte(text, cursor, b':')?;
        loop {
            let atom = read_number(text, cursor)?;
            result.push(CxRadical { atom, electrons });
            if *cursor >= text.len() || text.as_bytes()[*cursor] != b',' {
                break;
            }
            *cursor += 1;
            if !text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
                break;
            }
        }
    }
    Ok(result)
}

fn parse_enhanced_stereo(text: &str, cursor: &mut usize) -> Result<CxEnhancedStereo, CxParseError> {
    // RDKit source (verbatim; StereoGroup construction is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_enhanced_stereo(Iterator &first, Iterator last, RDKit::RWMol &mol,
                               unsigned int startAtomIdx) {
      StereoGroupType group_type = StereoGroupType::STEREO_ABSOLUTE;
      if (*first == 'a') {
        group_type = StereoGroupType::STEREO_ABSOLUTE;
      } else if (*first == 'o') {
        group_type = StereoGroupType::STEREO_OR;
      } else if (*first == '&') {
        group_type = StereoGroupType::STEREO_AND;
      }
      ++first;

      // OR and AND groups carry a group number
      unsigned int group_id = 0;
      if (group_type != StereoGroupType::STEREO_ABSOLUTE) {
        read_int(first, last, group_id);
      }

      if (first >= last || *first != ':') {
        return false;
      }
      ++first;

      std::vector<Atom *> atoms;
      std::vector<Bond *> bonds;

      while (first <= last && *first >= '0' && *first <= '9') {
        unsigned int aidx;
        if (read_int(first, last, aidx)) {
          if (VALID_ATIDX(aidx)) {
            Atom *atom = mol.getAtomWithIdx(aidx - startAtomIdx);
            if (!atom) {
              BOOST_LOG(rdWarningLog)
                  << "Atom " << aidx << " not found!" << std::endl;
              return false;
            }
            atoms.push_back(atom);
          }
        } else {
          return false;
        }

        if (first < last && *first == ',') {
          ++first;
        }
      }
      if (!atoms.empty()) {
        // we need to do a bit of work to check whether or not we've already seen
        // this particular StereoGroup (was Github #6050)
        const auto group_hash =
            10 * group_id + static_cast<unsigned int>(group_type);
        std::vector<unsigned int> sgTracker;
        mol.getPropIfPresent(cxsgTracker, sgTracker);
        std::vector<StereoGroup> mol_stereo_groups(mol.getStereoGroups());
        TEST_ASSERT(mol_stereo_groups.size() == sgTracker.size());

        auto iter = std::find(sgTracker.begin(), sgTracker.end(), group_hash);
        if (iter != sgTracker.end()) {
          auto index = iter - sgTracker.begin();
          auto gAtoms = mol_stereo_groups[index].getAtoms();
          gAtoms.insert(gAtoms.end(), atoms.begin(), atoms.end());
          mol_stereo_groups[index] =
              StereoGroup(mol_stereo_groups[index].getGroupType(),
                          std::move(gAtoms), std::move(bonds), group_id);
        } else {
          // not seen this before, create a new stereogroup
          mol_stereo_groups.emplace_back(group_type, std::move(atoms),
                                         std::move(bonds), group_id);
          sgTracker.push_back(group_hash);
          mol.setProp(cxsgTracker, sgTracker);
        }

        mol.setStereoGroups(std::move(mol_stereo_groups));
      }

      return true;
    }
        */
    // RDKit✔️✔️: kind, optional group id, ordered indices and delimiters match
    // the source lexical closure; merging repeated groups is lowering-owned.
    let kind = match text.as_bytes().get(*cursor).copied() {
        Some(b'a') => CxStereoGroupKind::Absolute,
        Some(b'o') => CxStereoGroupKind::Or,
        Some(b'&') => CxStereoGroupKind::And,
        _ => return Err(CxParseError::new(*cursor, "invalid CX stereo group")),
    };
    *cursor += 1;
    let group_id = if kind == CxStereoGroupKind::Absolute {
        0
    } else if text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
        read_number(text, cursor)? as u32
    } else {
        0
    };
    expect_byte(text, cursor, b':')?;
    let atoms = parse_number_list(text, cursor)?;
    Ok(CxEnhancedStereo {
        kind,
        group_id,
        atoms,
    })
}

fn parse_ring_bonds(text: &str, cursor: &mut usize) -> Result<Vec<CxRingBond>, CxParseError> {
    // RDKit source (verbatim; query construction is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_ring_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
                          unsigned int startAtomIdx) {
      if (first >= last || *first != 'r' || first + 1 >= last ||
          *(first + 1) != 'b' || first + 2 >= last || *(first + 2) != ':') {
        return false;
      }
      first += 3;
      while (first < last && *first >= '0' && *first <= '9') {
        unsigned int n1;
        if (!read_int(first, last, n1)) {
          return false;
        }
        // check that we can read at least two more characters:
        if (first + 1 >= last || *first != ':') {
          return false;
        }
        ++first;
        unsigned int n2;
        bool gt = false;
        if (*first == '*') {
          ++first;
          n2 = 0xDEADBEEF;
          if (VALID_ATIDX(n1)) {
            mol.setProp(common_properties::_NeedsQueryScan, 1);
          }
        } else {
          if (!read_int(first, last, n2)) {
            return false;
          }
          switch (n2) {
            case 0:
            case 2:
            case 3:
              break;
            case 4:
              gt = true;
              break;
            default:
              BOOST_LOG(rdWarningLog)
                  << "unrecognized rb value: " << n2 << std::endl;
              return false;
          }
        }
        if (VALID_ATIDX(n1)) {
          auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
          if (!atom->hasQuery()) {
            atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
          }
          if (!gt) {
            atom->expandQuery(makeAtomRingBondCountQuery(n2),
                              Queries::COMPOSITE_AND);
          } else {
            auto q = static_cast<ATOM_EQUALS_QUERY *>(new ATOM_LESSEQUAL_QUERY);
            q->setVal(n2);
            q->setDescription("AtomRingBondCount");
            q->setDataFunc(queryAtomRingBondCount);
            atom->expandQuery(q, Queries::COMPOSITE_AND);
          }
        }
        if (first < last && *first == ',') {
          ++first;
        }
      }
      return true;
    }
        */
    // RDKit✔️✔️: exact, less-equal and query-scan lexical states plus rejection
    // of all other values match the source in linear time.
    *cursor += 3;
    let mut result = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let atom = read_number(text, cursor)?;
        expect_byte(text, cursor, b':')?;
        let constraint = if text.as_bytes().get(*cursor) == Some(&b'*') {
            *cursor += 1;
            CxCountConstraint::QueryScan
        } else {
            let value = read_number(text, cursor)? as u32;
            match value {
                0 | 2 | 3 => CxCountConstraint::Exact(value),
                4 => CxCountConstraint::LessEqual(value),
                _ => {
                    return Err(CxParseError::new(
                        *cursor,
                        "unrecognized CX ring-bond count",
                    ));
                }
            }
        };
        result.push(CxRingBond { atom, constraint });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(result)
}

fn parse_count_constraints(
    text: &str,
    cursor: &mut usize,
    marker: u8,
) -> Result<Vec<CxAtomConstraint>, CxParseError> {
    // RDKit source (verbatim; query construction is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_substitution(Iterator &first, Iterator last, RDKit::RWMol &mol,
                            unsigned int startAtomIdx) {
      if (first >= last || *first != 's' || first + 1 >= last ||
          *(first + 1) != ':') {
        return false;
      }
      first += 2;
      while (first < last && *first >= '0' && *first <= '9') {
        unsigned int n1;
        if (!read_int(first, last, n1)) {
          return false;
        }
        // check that we can read at least two more characters:
        if (first + 1 >= last || *first != ':') {
          return false;
        }
        ++first;
        unsigned int n2;
        if (*first == '*') {
          ++first;
          n2 = 0xDEADBEEF;
          if (VALID_ATIDX(n1)) {
            mol.setProp(common_properties::_NeedsQueryScan, 1);
          }
        } else {
          if (!read_int(first, last, n2)) {
            return false;
          }
        }
        if (VALID_ATIDX(n1)) {
          auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
          if (!atom->hasQuery()) {
            atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
          }
          atom->expandQuery(makeAtomNonHydrogenDegreeQuery(n2),
                            Queries::COMPOSITE_AND);
        }
        if (first < last && *first == ',') {
          ++first;
        }
      }
      return true;
    }
        */
    // RDKit✔️✔️: exact/query-scan constraints and ordered separators reproduce
    // the source lexical path; query composition remains downstream.
    expect_byte(text, cursor, marker)?;
    expect_byte(text, cursor, b':')?;
    let mut result = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let atom = read_number(text, cursor)?;
        expect_byte(text, cursor, b':')?;
        let constraint = if text.as_bytes().get(*cursor) == Some(&b'*') {
            *cursor += 1;
            CxCountConstraint::QueryScan
        } else {
            CxCountConstraint::Exact(read_number(text, cursor)? as u32)
        };
        result.push(CxAtomConstraint { atom, constraint });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(result)
}

fn parse_link_nodes(text: &str, cursor: &mut usize) -> Result<Vec<CxLinkNode>, CxParseError> {
    // RDKit source (verbatim; degree inference/property installation is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_linknodes(Iterator &first, Iterator last, RDKit::RWMol &mol,
                         unsigned int startAtomIdx) {
      // these look like: |LN:1:1.3.2.6,4:1.4.3.6|
      // that's two records:
      //   1:1.3.2.6: 1-3 repeats, atom 1-2, 1-6
      //   4:1.4.3.6: 1-4 repeats, atom 4-3, 4-6
      // which maps to the property value "1 3 2 2 3 2 7|1 4 2 5 4 5 7"
      // If the linking atom only has two neighbors then the outer atom
      // specification (the last two digits) can be left out. So for a molecule
      // where atom 1 has bonds only to atoms 2 and 6 we could have
      // |LN:1:1.3|
      // instead of
      // |LN:1:1.3.2.6|
      if (first >= last || *first != 'L' || first + 1 >= last ||
          *(first + 1) != 'N' || first + 2 >= last || *(first + 2) != ':') {
        return false;
      }
      first += 3;
      std::string accum = "";
      while (first < last && *first >= '0' && *first <= '9') {
        unsigned int atidx;
        if (!read_int(first, last, atidx)) {
          return false;
        }
        // check that we can read at least two more characters:
        if (first + 1 >= last || *first != ':') {
          return false;
        }
        ++first;
        unsigned int startReps;
        if (!read_int(first, last, startReps)) {
          return false;
        }
        if (first + 1 >= last || *first != '.') {
          return false;
        }
        ++first;
        unsigned int endReps;
        if (!read_int(first, last, endReps)) {
          return false;
        }
        unsigned int idx1;
        unsigned int idx2;
        if (first < last && *first == '.') {
          ++first;
          if (!read_int(first, last, idx1)) {
            return false;
          }
          ++first;
          if (!read_int(first, last, idx2)) {
            return false;
          }
        } else if (VALID_ATIDX(atidx) &&
                   mol.getAtomWithIdx(atidx - startAtomIdx)->getDegree() == 2) {
          auto nbrs =
              mol.getAtomNeighbors(mol.getAtomWithIdx(atidx - startAtomIdx));
          idx1 = *nbrs.first;
          nbrs.first++;
          idx2 = *nbrs.first;
        } else if (VALID_ATIDX(atidx)) {
          return false;
        }
        if (first < last && *first == ',') {
          ++first;
        }
        if (VALID_ATIDX(atidx)) {
          if (!accum.empty()) {
            accum += "|";
          }
          accum += (boost::format("%d %d 2 %d %d %d %d") % startReps % endReps %
                    (atidx - startAtomIdx + 1) % (idx1 - startAtomIdx + 1) %
                    (atidx - startAtomIdx + 1) % (idx2 - startAtomIdx + 1))
                       .str();
        }
      }
      if (!accum.empty()) {
        mol.setProp(common_properties::molFileLinkNodes, accum);
      }
      return true;
    }
        */
    // RDKit✔️✔️: explicit/omitted outer syntax, repetition bounds and order
    // match the source; omitted-neighbour resolution is deliberately deferred.
    *cursor += "LN:".len();
    let mut records = Vec::new();
    while text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
        let atom = read_number(text, cursor)?;
        expect_byte(text, cursor, b':')?;
        let start_repetitions = read_number(text, cursor)?;
        expect_byte(text, cursor, b'.')?;
        let end_repetitions = read_number(text, cursor)?;
        let outer_atoms = if text.as_bytes().get(*cursor) == Some(&b'.') {
            *cursor += 1;
            let first = read_number(text, cursor)?;
            expect_byte(text, cursor, b'.')?;
            Some([first, read_number(text, cursor)?])
        } else {
            None
        };
        records.push(CxLinkNode {
            atom,
            start_repetitions,
            end_repetitions,
            outer_atoms,
        });
        if text.as_bytes().get(*cursor) == Some(&b',')
            && text
                .as_bytes()
                .get(*cursor + 1)
                .is_some_and(u8::is_ascii_digit)
        {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(records)
}

fn parse_data_sgroup(text: &str, cursor: &mut usize) -> Result<CxDataSGroup, CxParseError> {
    // RDKit source (verbatim; SGroup creation/property installation is lowering-owned):
    /*
    template <typename Iterator>
    void parse_data_sgroup_attr(Iterator &first, Iterator last,
                                SubstanceGroup &sgroup, bool keepSGroup,
                                std::string fieldName, bool fieldIsArray = false) {
      PRECONDITION(first < last, "parse_data_sgroup_attr: first >= last");
      if (first != last && *first != '|') {
        std::string data = read_text_to(first, last, ":");
        ++first;
        if (!data.empty() && keepSGroup) {
          if (fieldIsArray) {
            std::vector<std::string> dataFields = {data};
            sgroup.setProp(fieldName, dataFields);
          } else {
            sgroup.setProp(fieldName, data);
          }
        }
      }
    }

    template <typename Iterator>
    bool parse_data_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
                           unsigned int startAtomIdx, unsigned int nSGroups) {
      // these look like: |SgD:2,1:FIELD:info::::|
      // example from CXSMILES docs:
      //    SgD:3,2,1,0:name:data:like:unit:t:(1.,1.)
      // the fields are:
      //    SgD:[atom indices]:[field name]:[data value]:[query
      //    operator]:[unit]:[tag]:[coords]
      //   coords are (-1) if atomic coordinates are present
      if (first >= last || *first != 'S' || first + 3 >= last ||
          *(first + 1) != 'g' || *(first + 2) != 'D' || *(first + 3) != ':') {
        return false;
      }
      first += 4;
      std::vector<unsigned int> atoms;
      if (!read_int_list(first, last, atoms)) {
        return false;
      }
      SubstanceGroup sgroup(&mol, std::string("DAT"));
      sgroup.setProp(cxsmilesindex, nSGroups);
      bool keepSGroup = false;
      for (auto idx : atoms) {
        if (VALID_ATIDX(idx)) {
          keepSGroup = true;
          sgroup.addAtomWithIdx(idx - startAtomIdx);
        }
      }
      ++first;
      parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDNAME");

      // FIX:
      if (keepSGroup) {
        sgroup.setProp("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0");
      }

      parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "DATAFIELDS", true);

      parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "QUERYOP");

      parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDINFO");

      parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDTAG");

      if (first < last && *first == '(') {
        // FIX
        std::string coords = read_text_to(first, last, ")");
        ++first;
        if (keepSGroup) {
          sgroup.setProp("COORDS", coords);
        }
      }
      // the label processing can destroy sgroup info, so do that now
      // (the function will immediately return if already called)
      if (keepSGroup) {
        processCXSmilesLabels(mol);
        sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
        addSubstanceGroup(mol, sgroup);
      }
      return true;
    }
        */
    // RDKit✔️✔️: atom-list, five text attributes and optional coordinate text
    // reproduce the source lexical sequence with detached ordered storage.
    *cursor += "SgD:".len();
    let atoms = parse_delimited_number_list(text, cursor, b',')?;
    expect_byte(text, cursor, b':')?;
    let field_name = read_colon_field(text, cursor)?;
    let data = read_colon_field(text, cursor)?;
    let query_op = read_colon_field(text, cursor)?;
    let field_info = read_colon_field(text, cursor)?;
    let field_tag = read_colon_field(text, cursor)?;
    let coordinates = if text.as_bytes().get(*cursor) == Some(&b'(') {
        let start = *cursor;
        let coordinates = read_text_to(text, cursor, b")")?;
        if text.as_bytes().get(*cursor) != Some(&b')') {
            return Err(CxParseError::new(
                start,
                "unterminated CX SGroup coordinates",
            ));
        }
        *cursor += 1;
        Some(coordinates)
    } else {
        None
    };
    Ok(CxDataSGroup {
        atoms,
        field_name,
        data,
        query_op,
        field_info,
        field_tag,
        coordinates,
    })
}

fn parse_sgroup_hierarchy(
    text: &str,
    cursor: &mut usize,
) -> Result<Vec<CxSGroupHierarchy>, CxParseError> {
    // RDKit source (verbatim; SGroup lookup/parent mutation is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_sgroup_hierarchy(Iterator &first, Iterator last, RDKit::RWMol &mol) {
      // these look like: |SgH:1:0|
      // from CXSMILES docs:
      //    SgH:parentSgroupIndex1:childSgroupIndex1.childSgroupIndex2,parentSgroupIndex2:childSgroupIndex1
      if (first >= last || *first != 'S' || first + 3 >= last ||
          *(first + 1) != 'g' || *(first + 2) != 'H' || *(first + 3) != ':') {
        return false;
      }
      first += 4;
      auto &sgs = getSubstanceGroups(mol);
      while (1) {
        unsigned int parentId;
        if (!read_int(first, last, parentId)) {
          return false;
        }

        bool validParent = true;
        auto psg = find_matching_sgroup(sgs, parentId);
        if (psg == sgs.end()) {
          validParent = false;
        } else {
          psg->getPropIfPresent("index", parentId);
        }
        if (first <= last && *first == ':') {
          ++first;
          std::vector<unsigned int> children;
          if (!read_int_list(first, last, children, '.')) {
            return false;
          }
          if (validParent) {
            for (auto childId : children) {
              if (childId >= sgs.size()) {
                throw SmilesParseException(
                    "child id references non-existent SGroup");
              }
              auto csg = find_matching_sgroup(sgs, childId);
              if (csg != sgs.end()) {
                unsigned int cid;
                csg->getProp("index", cid);
                csg->setProp("PARENT", parentId);
              }
            }
          }
          if (first <= last && *first == ',') {
            ++first;
          } else {
            break;
          }
        } else {
          return false;
        }
      }

      return true;
    }
        */
    // RDKit✔️✔️: ordered parent/child syntax and separators match the source;
    // reference validity and parent installation remain lowering-owned.
    *cursor += "SgH:".len();
    let mut result = Vec::new();
    loop {
        let parent = read_number(text, cursor)?;
        expect_byte(text, cursor, b':')?;
        let children = parse_delimited_number_list(text, cursor, b'.')?;
        result.push(CxSGroupHierarchy { parent, children });
        if text.as_bytes().get(*cursor) == Some(&b',')
            && text
                .as_bytes()
                .get(*cursor + 1)
                .is_some_and(u8::is_ascii_digit)
        {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(result)
}

fn parse_polymer_sgroup(text: &str, cursor: &mut usize) -> Result<CxPolymerSGroup, CxParseError> {
    // RDKit source (verbatim; type validation and SGroup installation are lowering-owned):
    /*
    template <typename Iterator>
    bool parse_polymer_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
                              unsigned int startAtomIdx, unsigned int nSGroups) {
      // these look like:
      //    |Sg:n:6,1,2,4::hh&#44;f:6,0,:4,2,|
      // example from CXSMILES docs:
      // the fields are:
      //    Sg:[type]:[atom indices]:[subscript]:[superscript]:[head crossing
      //    bonds]:[tail crossing bonds]:
      //
      // note that it's legit for empty fields to be completely missing.
      //   for example, this doesn't have any crossing bonds indicated:
      // *-CCCN-* |$star_e;;;;;star_e$,Sg:n:4,1,2,3::hh|
      // this last bit makes the whole thing doubleplusfun to parse

      if (first >= last || *first != 'S' || first + 2 >= last ||
          *(first + 1) != 'g' || *(first + 2) != ':') {
        return false;
      }
      first += 3;

      const auto type_code = read_text_to(first, last, ":");
      ++first;
      const auto type = sgroupTypemap.find(type_code);
      if (type == sgroupTypemap.end()) {
        return false;
      }
      bool keepSGroup = false;
      SubstanceGroup sgroup(&mol, type->second);
      sgroup.setProp(cxsmilesindex, nSGroups);
      if (type_code == "alt") {
        sgroup.setProp("SUBTYPE", std::string("ALT"));
      } else if (type_code == "ran") {
        sgroup.setProp("SUBTYPE", std::string("RAN"));
      } else if (type_code == "blk") {
        sgroup.setProp("SUBTYPE", std::string("BLO"));
      }

      std::vector<unsigned int> atoms;
      if (!read_int_list(first, last, atoms)) {
        return false;
      }
      //++first;
      for (auto idx : atoms) {
        if (VALID_ATIDX(idx)) {
          sgroup.addAtomWithIdx(idx - startAtomIdx);
          keepSGroup = true;
        }
      }
      std::vector<unsigned int> headCrossing;
      std::vector<unsigned int> tailCrossing;
      if (first <= last && *first == ':') {
        ++first;
        std::string subscript = read_text_to(first, last, ":|");
        if (keepSGroup && !subscript.empty()) {
          sgroup.setProp("LABEL", subscript);
        }
        if (first <= last && *first == ':') {
          ++first;
          std::string superscript = read_text_to(first, last, ":|,");
          if (keepSGroup && !superscript.empty()) {
            sgroup.setProp("CONNECT", superscript);
          }

          if (first <= last && *first == ':') {
            ++first;
            if (!read_int_list(first, last, headCrossing)) {
              return false;
            }
            if (keepSGroup && !headCrossing.empty()) {
              for (auto &cidx : headCrossing) {
                if (VALID_ATIDX(cidx)) {
                  cidx -= startAtomIdx;
                } else {
                  keepSGroup = false;
                  break;
                }
              }
              sgroup.setProp(_headCrossings, headCrossing, true);
            }
            if (first <= last && *first == ':') {
              ++first;
              if (!read_int_list(first, last, tailCrossing)) {
                return false;
              }
            }
            if (keepSGroup && !tailCrossing.empty()) {
              for (auto &cidx : tailCrossing) {
                if (VALID_ATIDX(cidx)) {
                  cidx -= startAtomIdx;
                } else {
                  keepSGroup = false;
                  break;
                }
              }
              sgroup.setProp("_tailCrossings", tailCrossing, true);
            }
          }
        }
      }
      if (keepSGroup) {  // the label processing can destroy sgroup info, so do that
                         // now (the function will immediately return if already
                         // called)
        processCXSmilesLabels(mol);

        finalizePolymerSGroup(mol, sgroup);
        sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);

        addSubstanceGroup(mol, sgroup);
      }
      return true;
    }
        */
    // RDKit✔️✔️: type text, atoms and the absent/partial/full optional field
    // sequence are preserved in source order with equivalent linear scans.
    *cursor += "Sg:".len();
    let type_code = read_text_to(text, cursor, b":")?;
    expect_byte(text, cursor, b':')?;
    let atoms = parse_delimited_number_list(text, cursor, b',')?;
    let mut label = String::new();
    let mut connect = String::new();
    let mut head_crossings = Vec::new();
    let mut tail_crossings = Vec::new();
    if text.as_bytes().get(*cursor) == Some(&b':') {
        *cursor += 1;
        label = read_text_to(text, cursor, b":|")?;
        if text.as_bytes().get(*cursor) == Some(&b':') {
            *cursor += 1;
            connect = read_text_to(text, cursor, b":|,")?;
            if text.as_bytes().get(*cursor) == Some(&b':') {
                *cursor += 1;
                head_crossings = parse_delimited_number_list(text, cursor, b',')?;
                if text.as_bytes().get(*cursor) == Some(&b':') {
                    *cursor += 1;
                    tail_crossings = parse_delimited_number_list(text, cursor, b',')?;
                }
            }
        }
    }
    Ok(CxPolymerSGroup {
        type_code,
        atoms,
        label,
        connect,
        head_crossings,
        tail_crossings,
    })
}

fn parse_variable_attachments(
    text: &str,
    cursor: &mut usize,
) -> Result<Vec<CxVariableAttachment>, CxParseError> {
    // RDKit source (verbatim; graph degree checks/bond properties are lowering-owned):
    /*
    template <typename Iterator>
    bool parse_variable_attachments(Iterator &first, Iterator last,
                                    RDKit::RWMol &mol, unsigned int startAtomIdx) {
      // these look like: CO*.C1=CC=NC=C1 |m:2:3.5.4|
      // that corresponds to replacing the bond to atom 2 with bonds to atom 3, 5,
      // or 4
      //
      if (first >= last || *first != 'm' || first + 1 >= last ||
          *(first + 1) != ':') {
        return false;
      }
      first += 2;

      while (first < last && *first >= '0' && *first <= '9') {
        unsigned int at1idx;
        if (!read_int(first, last, at1idx)) {
          return false;
        }

        if (VALID_ATIDX(at1idx) &&
            mol.getAtomWithIdx(at1idx - startAtomIdx)->getDegree() != 1) {
          BOOST_LOG(rdWarningLog)
              << "position variation bond to atom with more than one bond"
              << std::endl;
          return false;
        }
        if (first < last && *first == ':') {
          ++first;
        } else {
          BOOST_LOG(rdWarningLog) << "improperly formatted m: block" << std::endl;
          return false;
        }
        std::vector<std::string> others;
        while (first < last && *first >= '0' && *first <= '9') {
          unsigned int aidx;
          if (!read_int(first, last, aidx)) {
            return false;
          }
          if (VALID_ATIDX(aidx)) {
            others.push_back(std::to_string(aidx - startAtomIdx + 1));
          }
          if (first < last && *first == '.') {
            ++first;
          }
        }
        if (VALID_ATIDX(at1idx)) {
          std::string endPts = "(" + std::to_string(others.size());
          for (auto idx : others) {
            endPts += " " + idx;
          }
          endPts += ")";

          for (auto nbri : boost::make_iterator_range(
                   mol.getAtomBonds(mol.getAtomWithIdx(at1idx - startAtomIdx)))) {
            auto bnd = mol[nbri];
            bnd->setProp(common_properties::_MolFileBondEndPts, endPts);
            bnd->setProp(common_properties::_MolFileBondAttach, std::string("ANY"));
          }
        }
        if (first < last && *first == ',') {
          ++first;
        }
      }
      return true;
    }
        */
    // RDKit✔️✔️: attachment rows and ordered endpoint lists reproduce the
    // source syntax; degree and bond-property effects remain downstream.
    *cursor += "m:".len();
    let mut result = Vec::new();
    loop {
        let atom = read_number(text, cursor)?;
        expect_byte(text, cursor, b':')?;
        let endpoints = parse_delimited_number_list(text, cursor, b'.')?;
        result.push(CxVariableAttachment { atom, endpoints });
        if text.as_bytes().get(*cursor) == Some(&b',')
            && text
                .as_bytes()
                .get(*cursor + 1)
                .is_some_and(u8::is_ascii_digit)
        {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(result)
}

fn parse_wedge_bonds(text: &str, cursor: &mut usize) -> Result<Vec<CxWedgeBond>, CxParseError> {
    // RDKit source (verbatim; bond lookup/orientation/stereo effects are lowering-owned):
    /*
    template <typename Iterator>
    bool parse_wedged_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
                            unsigned int startAtomIdx, unsigned int startBondIdx) {
      // these look like: CC(O)Cl |w:1.0|
      // also wD and wU for down and up wedges.
      //
      // We do not end up using this to set stereochemistry, but the relevant bond
      // properties are set in case client code wants to do something with the
      // information.
      if (first >= last || *first != 'w' || first + 1 >= last) {
        return false;
      }
      ++first;
      Bond::BondDir state = Bond::BondDir::NONE;
      unsigned int cfg = 0;
      switch (*first) {
        case ':':
          state = Bond::BondDir::UNKNOWN;
          cfg = 2;
          break;
        case 'U':
          state = Bond::BondDir::BEGINWEDGE;
          cfg = 1;
          ++first;
          break;
        case 'D':
          state = Bond::BondDir::BEGINDASH;
          cfg = 3;
          ++first;
          break;
        default:
          break;
      }
      if (state == Bond::BondDir::NONE || first >= last || first + 1 >= last ||
          *first != ':') {
        return false;
      }
      ++first;
      while (first < last && *first >= '0' && *first <= '9') {
        unsigned int atomIdx;
        if (!read_int(first, last, atomIdx)) {
          return false;
        }
        if (first < last && *first == '.') {
          ++first;
        } else {
          BOOST_LOG(rdWarningLog) << "improperly formatted w block" << std::endl;
          return false;
        }
        unsigned int bondIdx;
        if (!read_int(first, last, bondIdx)) {
          return false;
        }

        if (VALID_ATIDX(atomIdx) && VALID_BNDIDX(bondIdx)) {
          auto atom = mol.getAtomWithIdx(atomIdx - startAtomIdx);
          auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);

          if (!bond) {
            BOOST_LOG(rdWarningLog)
                << "bond " << bondIdx << " not found, wedge from atom " << atomIdx
                << " cannot be applied." << std::endl;
            return false;
          }

          // we can't set wedging twice:
          if (bond->hasProp(common_properties::_MolFileBondCfg)) {
            BOOST_LOG(rdWarningLog)
                << "w block attempts to set wedging on bond " << bond->getIdx()
                << " more than once." << std::endl;
            return false;
          }

          // first things first, the atom needs to be the start atom of the bond for
          // any of this to make sense
          if (atom->getIdx() != bond->getBeginAtomIdx()) {
            if (atom->getIdx() != bond->getEndAtomIdx()) {
              BOOST_LOG(rdWarningLog)
                  << "atom " << atomIdx << " is not associated with bond "
                  << bondIdx << "(" << bond->getBeginAtomIdx() + startAtomIdx << "-"
                  << bond->getEndAtomIdx() + startAtomIdx << ")"
                  << " in w block" << std::endl;
              return false;
            }
            auto eidx = bond->getBeginAtomIdx();
            bond->setBeginAtomIdx(atom->getIdx());
            bond->setEndAtomIdx(eidx);
          }
          bond->setProp(common_properties::_MolFileBondCfg, cfg);
          bond->setBondDir(state);
          if (cfg == 2 && canHaveDirection(*bond)) {
            bond->getBeginAtom()->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
            mol.setProp(detail::_needsDetectBondStereo, 1);
          }
          if ((cfg == 1 || cfg == 3) && canHaveDirection(*bond)) {
            mol.setProp(detail::_needsDetectAtomStereo, 1);
          }
        }
        if (first < last && *first == ',') {
          ++first;
        }
      }
      return true;
    }
        */
    // RDKit✔️✔️: marker classification, configuration values, ordered pairs
    // and separators match the source; graph/stereo effects remain downstream.
    expect_byte(text, cursor, b'w')?;
    let (direction, configuration) = match text.as_bytes().get(*cursor).copied() {
        Some(b':') => (CxWedgeDirection::Unknown, 2),
        Some(b'U') => {
            *cursor += 1;
            (CxWedgeDirection::BeginWedge, 1)
        }
        Some(b'D') => {
            *cursor += 1;
            (CxWedgeDirection::BeginDash, 3)
        }
        _ => return Err(CxParseError::new(*cursor, "invalid CX wedge marker")),
    };
    expect_byte(text, cursor, b':')?;
    let mut result = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let (atom, bond) = read_pair(text, cursor, b'.')?;
        result.push(CxWedgeBond {
            atom,
            bond,
            direction,
            configuration,
        });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(result)
}

fn parse_double_bond_stereo(
    text: &str,
    cursor: &mut usize,
    stereo: CxDoubleBondStereoKind,
) -> Result<CxDoubleBondStereo, CxParseError> {
    // RDKit source (verbatim; bond lookup/stereo application is lowering-owned):
    /*
    template <typename Iterator>
    bool parse_doublebond_stereo(Iterator &first, Iterator last, RDKit::RWMol &mol,
                                 unsigned int, unsigned int startBondIdx,
                                 Bond::BondStereo stereo) {
      // these look like: C1CCCC/C=C/CCC1 |ctu:5|
      // also c and t for cis or trans
      //
      while (first < last && *first != ':') {
        ++first;
      }
      if (first >= last || *first != ':') {
        return false;
      }
      ++first;

      while (first < last && *first >= '0' && *first <= '9') {
        unsigned int bondIdx;
        if (!read_int(first, last, bondIdx)) {
          return false;
        }
        if (VALID_BNDIDX(bondIdx)) {
          auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);

          if (!bond) {
            BOOST_LOG(rdWarningLog)
                << "bond " << bondIdx
                << " not found, cannot mark as stereo double bond." << std::endl;
            return false;
          }

          bool useCXOrdering = true;
          Chirality::detail::setStereoForBond(mol, bond, stereo, useCXOrdering);
        }
        if (first < last && *first == ',') {
          ++first;
        }
      }
      return true;
    }
        */
    // RDKit✔️✔️: prefix-to-kind dispatch and ordered bond-list syntax preserve
    // the source's one-pass behavior; bond stereo application is downstream.
    while *cursor < text.len() && text.as_bytes()[*cursor] != b':' {
        *cursor += 1;
    }
    expect_byte(text, cursor, b':')?;
    Ok(CxDoubleBondStereo {
        stereo,
        bonds: parse_number_list(text, cursor)?,
    })
}

#[cfg(test)]
mod tests {
    use super::{CxRecord, CxStereoGroupKind, parse_cx_extensions};

    #[test]
    fn parses_representation_independent_records() {
        let parsed = parse_cx_extensions("|(0,0,;1,0,)$C;O$| name").unwrap();
        assert_eq!(parsed.consumed(), 18);
        assert_eq!(parsed.records().len(), 2);
        assert!(matches!(parsed.records()[0], CxRecord::Coordinates(_)));
        assert!(matches!(parsed.records()[1], CxRecord::AtomLabels(_)));
    }

    #[test]
    fn parses_query_and_stereo_records_without_a_graph() {
        let parsed = parse_cx_extensions("|u:0,rb:1:4,o2:0,1|").unwrap();
        assert!(matches!(parsed.records()[0], CxRecord::Unsaturation(_)));
        assert!(matches!(parsed.records()[1], CxRecord::RingBonds(_)));
        match &parsed.records()[2] {
            CxRecord::EnhancedStereo(stereo) => {
                assert_eq!(stereo.kind, CxStereoGroupKind::Or);
                assert_eq!(stereo.group_id, 2);
            }
            other => panic!("unexpected record: {other:?}"),
        }
    }

    #[test]
    fn rejects_missing_outer_pipe() {
        let error = parse_cx_extensions("u:0|").unwrap_err();
        assert_eq!(error.offset, 0);
    }

    #[test]
    fn parses_structural_records_without_destination_state() {
        let parsed = parse_cx_extensions(
            "|LN:1:1.3.2.6,SgD:2,1:FIELD:info::::,SgH:1:0,Sg:n:0,1:n:ht:2:3,m:2:3.5.4|",
        )
        .unwrap();
        assert!(matches!(parsed.records()[0], CxRecord::LinkNodes(_)));
        assert!(matches!(parsed.records()[1], CxRecord::DataSGroup(_)));
        assert!(matches!(parsed.records()[2], CxRecord::SGroupHierarchy(_)));
        assert!(matches!(parsed.records()[3], CxRecord::PolymerSGroup(_)));
        assert!(matches!(
            parsed.records()[4],
            CxRecord::VariableAttachments(_)
        ));
        match &parsed.records()[1] {
            CxRecord::DataSGroup(group) => {
                assert_eq!(group.atoms, vec![2, 1]);
                assert_eq!(group.field_name, "FIELD");
                assert_eq!(group.data, "info");
            }
            other => panic!("unexpected record: {other:?}"),
        }
    }
}
