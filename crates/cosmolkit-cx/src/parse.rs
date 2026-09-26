use crate::scan::{expect_byte, parse_delimited_number_list, read_number, read_pair, read_text_to};
use crate::{
    CxAtomConstraint, CxAtomProperty, CxBondReference, CxCoordinateBondKind, CxCoordinateBonds,
    CxCoordinates, CxCountConstraint, CxDataSGroup, CxDoubleBondStereo, CxDoubleBondStereoKind,
    CxEnhancedStereo, CxLinkNode, CxParseError, CxParseProgress, CxPolymerSGroup,
    CxProgressCheckpoint, CxProgressPhase, CxRadical, CxRecord, CxRingBond, CxSGroupHierarchy,
    CxStereoGroupKind, CxVariableAttachment, CxWedgeBond, CxWedgeDirection, ParsedCxExtensions,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CxDispatch {
    Coordinates,
    AtomLabelsOrValues,
    AtomProperties,
    CoordinateBonds,
    ZeroBonds,
    Radicals,
    EnhancedStereo,
    RingBonds,
    LinkNodes,
    DataSGroup,
    SGroupHierarchy,
    PolymerSGroup,
    Unsaturation,
    Substitution,
    VariableAttachments,
    WedgedBonds,
    DoubleBondAny,
    DoubleBondCis,
    DoubleBondTrans,
    Unknown,
}

fn classify_cx_dispatch(bytes: &[u8], cursor: usize) -> CxDispatch {
    // RDKit source (verbatim; source-order dispatch recognition):
    /*
    if (*first == '(') {
    } else if (*first == '$') {
      if (length > 4 && *(first + 1) == '_' && *(first + 2) == 'A' &&
          *(first + 3) == 'V' && *(first + 4) == ':') {
      } else {
      }
    } else if (length > 9 && std::string(first, first + 9) == "atomProp:") {
    } else if (*first == 'C') {
    } else if (*first == 'H') {
    } else if (*first == 'Z') {
    } else if (*first == '^') {
    } else if (*first == 'a' || *first == 'o' ||
               (*first == '&' && first + 1 < last && first[1] != '#')) {
    } else if (*first == 'r' && first + 1 < last && first[1] == 'b') {
    } else if (*first == 'L' && first + 1 < last && first[1] == 'N') {
    } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
               first[2] == 'D') {
    } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
               first[2] == 'H') {
    } else if (*first == 'S' && first + 1 < last && first[1] == 'g') {
    } else if (*first == 'u') {
    } else if (*first == 's') {
    } else if (*first == 'm') {
    } else if (*first == 'w') {
    } else if (*first == 'c' && first + 2 < last && first[1] == 't' &&
               first[2] == 'u') {
    } else if (*first == 'c') {
    } else if (*first == 't') {
    } else {
      ++first;
    }
        */
    let Some(&first) = bytes.get(cursor) else {
        return CxDispatch::Unknown;
    };
    let remaining = bytes.len() - cursor;
    match first {
        b'(' => CxDispatch::Coordinates,
        b'$' if remaining > 4 && bytes[cursor + 1..cursor + 5] == *b"_AV:" => {
            CxDispatch::AtomLabelsOrValues
        }
        b'$' => CxDispatch::AtomLabelsOrValues,
        b'a' if remaining > 9 && bytes[cursor..cursor + 9] == *b"atomProp:" => {
            CxDispatch::AtomProperties
        }
        b'a' | b'o' => CxDispatch::EnhancedStereo,
        b'&' if remaining > 1 && bytes[cursor + 1] != b'#' => CxDispatch::EnhancedStereo,
        b'r' if remaining > 1 && bytes[cursor + 1] == b'b' => CxDispatch::RingBonds,
        b'L' if remaining > 1 && bytes[cursor + 1] == b'N' => CxDispatch::LinkNodes,
        b'S' if remaining > 2 && bytes[cursor + 1..cursor + 3] == *b"gD" => CxDispatch::DataSGroup,
        b'S' if remaining > 2 && bytes[cursor + 1..cursor + 3] == *b"gH" => {
            CxDispatch::SGroupHierarchy
        }
        b'S' if remaining > 1 && bytes[cursor + 1] == b'g' => CxDispatch::PolymerSGroup,
        b'c' if remaining > 2 && bytes[cursor + 1..cursor + 3] == *b"tu" => {
            CxDispatch::DoubleBondAny
        }
        b'c' => CxDispatch::DoubleBondCis,
        b't' => CxDispatch::DoubleBondTrans,
        b'C' | b'H' => CxDispatch::CoordinateBonds,
        b'Z' => CxDispatch::ZeroBonds,
        b'^' => CxDispatch::Radicals,
        b'u' => CxDispatch::Unsaturation,
        b's' => CxDispatch::Substitution,
        b'm' => CxDispatch::VariableAttachments,
        b'w' => CxDispatch::WedgedBonds,
        _ => CxDispatch::Unknown,
    }
    // RDKit✔️✔️: first-byte and bounded lookahead predicates follow source
    // priority; byte-slice matching keeps non-ASCII input safe in linear time.
}

/// Parse one CX extension block without referring to a destination graph.
pub fn parse_cx_extensions(text: &str) -> Result<ParsedCxExtensions, CxParseError> {
    // RDKit source (verbatim; this is the strict API adapter):
    /*
    void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
                           std::string::const_iterator &first,
                           unsigned int startAtomIdx, unsigned int startBondIdx) {
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
    let progress = parse_cx_extensions_progress(text);
    let (records, _, consumed, complete, error) = progress.into_parts();
    if !complete {
        return Err(error.unwrap_or_else(|| {
            CxParseError::new(consumed, "failure parsing CXSMILES extensions")
        }));
    }
    Ok(ParsedCxExtensions::new(records, consumed))
}

/// Parse one CX block while retaining its source cursor and committed record
/// checkpoints, including when a later helper fails.
pub fn parse_cx_extensions_progress(text: &str) -> CxParseProgress {
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

        */
    // RDKit❗✔️: the Rust loop now follows source dispatch priority, record
    // order, conformer ordinals, bytewise unknown advancement, pipe contract
    // and linear cursor complexity; helper-specific partial commits are added
    // by the following source-ordered family steps.
    if text.is_empty() {
        return CxParseProgress::from_parts(Vec::new(), Vec::new(), 0, true, None);
    }
    let bytes = text.as_bytes();
    if bytes.first().copied() != Some(b'|') {
        return CxParseProgress::from_parts(
            Vec::new(),
            Vec::new(),
            0,
            false,
            Some(CxParseError::new(
                0,
                "CXSMILES extension does not start with |",
            )),
        );
    }

    let mut cursor = 1;
    let mut conformer = 0;
    let mut records = Vec::new();
    let mut checkpoints = Vec::new();
    while cursor < bytes.len() && bytes[cursor] != b'|' {
        let start = cursor;
        let dispatch = classify_cx_dispatch(bytes, cursor);
        if dispatch == CxDispatch::Unknown {
            while cursor < bytes.len()
                && bytes[cursor] != b'|'
                && classify_cx_dispatch(bytes, cursor) == CxDispatch::Unknown
            {
                // RDKit's default branch increments its iterator by one byte.
                cursor += 1;
            }
            let raw = text
                .get(start..cursor)
                .expect("an unknown CX span ends at a UTF-8 boundary")
                .to_owned();
            records.push(CxRecord::Unknown(raw));
            continue;
        }
        if dispatch == CxDispatch::Coordinates {
            let conformer_index = conformer;
            conformer += 1;
            if let Err(error) = parse_coordinates_progress(
                text,
                &mut cursor,
                conformer_index,
                &mut records,
                &mut checkpoints,
            ) {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::AtomLabelsOrValues {
            if let Err(error) =
                parse_labels_or_values_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::AtomProperties {
            if let Err(error) =
                parse_atom_properties_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::CoordinateBonds {
            if let Err(error) =
                parse_coordinate_bonds_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::ZeroBonds {
            if let Err(error) =
                parse_zero_bonds_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::Radicals {
            if let Err(error) =
                parse_radicals_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::EnhancedStereo {
            if let Err(error) =
                parse_enhanced_stereo_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::Unsaturation {
            if let Err(error) =
                parse_unsaturation_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::RingBonds {
            if let Err(error) =
                parse_ring_bonds_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::Substitution {
            if let Err(error) =
                parse_substitution_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::LinkNodes {
            if let Err(error) =
                parse_link_nodes_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::DataSGroup {
            if let Err(error) =
                parse_data_sgroup_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::SGroupHierarchy {
            if let Err(error) =
                parse_sgroup_hierarchy_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::PolymerSGroup {
            if let Err(error) =
                parse_polymer_sgroup_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::VariableAttachments {
            if let Err(error) = parse_variable_attachments_progress(
                text,
                &mut cursor,
                &mut records,
                &mut checkpoints,
            ) {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if dispatch == CxDispatch::WedgedBonds {
            if let Err(error) =
                parse_wedge_bonds_progress(text, &mut cursor, &mut records, &mut checkpoints)
            {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        if matches!(
            dispatch,
            CxDispatch::DoubleBondAny | CxDispatch::DoubleBondCis | CxDispatch::DoubleBondTrans
        ) {
            let stereo = match dispatch {
                CxDispatch::DoubleBondAny => CxDoubleBondStereoKind::Any,
                CxDispatch::DoubleBondCis => CxDoubleBondStereoKind::Cis,
                CxDispatch::DoubleBondTrans => CxDoubleBondStereoKind::Trans,
                _ => unreachable!("only double-bond stereo dispatch reaches this branch"),
            };
            if let Err(error) = parse_double_bond_stereo_progress(
                text,
                &mut cursor,
                &mut records,
                &mut checkpoints,
                stereo,
            ) {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
            continue;
        }
        let record = (|| -> Result<CxRecord, CxParseError> {
            match dispatch {
                CxDispatch::Coordinates => unreachable!("coordinates are handled above"),
                CxDispatch::AtomLabelsOrValues => {
                    unreachable!("atom labels and values are handled above")
                }
                CxDispatch::AtomProperties => unreachable!("atom properties are handled above"),
                CxDispatch::CoordinateBonds => {
                    unreachable!("coordinate bonds are handled above")
                }
                CxDispatch::ZeroBonds => unreachable!("zero bonds are handled above"),
                CxDispatch::Radicals => unreachable!("radicals are handled above"),
                CxDispatch::EnhancedStereo => unreachable!("enhanced stereo is handled above"),
                CxDispatch::Unsaturation => unreachable!("unsaturation is handled above"),
                CxDispatch::RingBonds => unreachable!("ring bonds are handled above"),
                CxDispatch::Substitution => unreachable!("substitution is handled above"),
                CxDispatch::LinkNodes => unreachable!("link nodes are handled above"),
                CxDispatch::DataSGroup => unreachable!("data SGroups are handled above"),
                CxDispatch::SGroupHierarchy => {
                    unreachable!("SGroup hierarchy records are handled above")
                }
                CxDispatch::PolymerSGroup => {
                    unreachable!("polymer SGroups are handled above")
                }
                CxDispatch::VariableAttachments => {
                    unreachable!("variable attachments are handled above")
                }
                CxDispatch::WedgedBonds => {
                    unreachable!("wedged bonds are handled above")
                }
                CxDispatch::DoubleBondAny
                | CxDispatch::DoubleBondCis
                | CxDispatch::DoubleBondTrans => {
                    unreachable!("double-bond stereo is handled above")
                }
                CxDispatch::Unknown => unreachable!("unknown dispatch is consumed above"),
            }
        })();
        let record = match record {
            Ok(record) => record,
            Err(error) => {
                return CxParseProgress::from_parts(
                    records,
                    checkpoints,
                    cursor,
                    false,
                    Some(error),
                );
            }
        };
        if cursor == start {
            return CxParseProgress::from_parts(
                records,
                checkpoints,
                cursor,
                false,
                Some(CxParseError::new(cursor, "CX parser made no progress")),
            );
        }
        records.push(record);
        checkpoints.push(CxProgressCheckpoint {
            record_index: records.len() - 1,
            item_index: None,
            cursor,
            phase: CxProgressPhase::Complete,
        });
    }
    if cursor >= bytes.len() || bytes[cursor] != b'|' {
        return CxParseProgress::from_parts(
            records,
            checkpoints,
            cursor,
            false,
            Some(CxParseError::new(
                cursor,
                "failure parsing CXSMILES extensions",
            )),
        );
    }
    CxParseProgress::from_parts(records, checkpoints, cursor + 1, true, None)
}

fn parse_coordinates_progress(
    text: &str,
    cursor: &mut usize,
    conformer: usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; destination mutations are emitted as progress):
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
    // RDKit❗✔️: a conformer begins before scanning; each completed coordinate
    // row commits before its separator is consumed. The final dimensionality
    // update precedes the closing-parenthesis check. The byte scan and row
    // component visits are linear, with no prefix reparse or whole-input copy.
    let start = *cursor;
    if *cursor >= text.len() || text.as_bytes()[*cursor] != b'(' {
        return Err(CxParseError::new(*cursor, "invalid CX coordinate record"));
    }
    let record_index = records.len();
    records.push(CxRecord::Coordinates(CxCoordinates {
        conformer,
        values: Vec::new(),
        // RDKit Conformer starts with df_is3D=true; parse_coords only resets
        // it after the row loop, so a component conversion failure keeps it.
        is_3d: true,
    }));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });
    *cursor += 1;
    let mut has_z_component = false;
    while *cursor < text.len() && text.as_bytes()[*cursor] != b')' {
        let field_start = *cursor;
        while *cursor < text.len() && !matches!(text.as_bytes()[*cursor], b';' | b')') {
            *cursor += 1;
        }
        let field = &text[field_start..*cursor];
        let (value, has_z) = parse_coordinate_row(field, field_start)?;
        has_z_component |= has_z;
        let item_index = match &mut records[record_index] {
            CxRecord::Coordinates(coordinates) => {
                let item_index = coordinates.values.len();
                coordinates.values.push(value);
                item_index
            }
            _ => unreachable!("coordinate progress record changed kind"),
        };
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b';' {
            *cursor += 1;
        }
    }
    if let CxRecord::Coordinates(coordinates) = &mut records[record_index] {
        coordinates.is_3d = has_z_component
            && coordinates
                .values
                .iter()
                .flatten()
                .any(|point| point[2].abs() > 1e-3);
    }
    if *cursor >= text.len() || text.as_bytes()[*cursor] != b')' {
        return Err(CxParseError::new(
            start,
            "unterminated CX coordinate record",
        ));
    }
    *cursor += 1;
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_coordinate_row(
    field: &str,
    offset: usize,
) -> Result<(Option<[f64; 3]>, bool), CxParseError> {
    if field.is_empty() {
        return Ok((None, false));
    }
    let mut parts = field.split(',');
    let x = parse_coordinate_component(parts.next(), offset)?;
    let y = parse_coordinate_component(parts.next(), offset)?;
    let z = parse_coordinate_component(parts.next(), offset)?;
    let has_z = z.is_some();
    Ok((
        Some([
            x.unwrap_or_default(),
            y.unwrap_or_default(),
            z.unwrap_or_default(),
        ]),
        has_z,
    ))
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

fn parse_labels_or_values_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; each property write is an item checkpoint):
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
    // RDKit❗✔️: each nonempty decoded slot is retained and emitted at its
    // source write cursor; empty slots advance the atom index without a write.
    // The parser makes one forward pass with no prefix reparse. Progress adds
    // an owned record vector and one checkpoint per nonempty slot.
    let value = text
        .get(*cursor..)
        .is_some_and(|remaining| remaining.starts_with("$_AV:"));
    if value {
        *cursor += 4;
        expect_byte(text, cursor, b':')?;
    } else {
        expect_byte(text, cursor, b'$')?;
    }
    let record_index = records.len();
    records.push(if value {
        CxRecord::AtomValues(Vec::new())
    } else {
        CxRecord::AtomLabels(Vec::new())
    });
    let mut atom_index = 0;
    while *cursor < text.len() && text.as_bytes()[*cursor] != b'$' {
        let field = read_text_to(text, cursor, b";$")?;
        let written = !field.is_empty();
        match &mut records[record_index] {
            CxRecord::AtomValues(fields) | CxRecord::AtomLabels(fields) => {
                fields.push(if written { Some(field) } else { None });
            }
            _ => unreachable!("atom progress record changed kind"),
        }
        if written {
            checkpoints.push(CxProgressCheckpoint {
                record_index,
                item_index: Some(atom_index),
                cursor: *cursor,
                phase: CxProgressPhase::Item,
            });
        }
        atom_index += 1;
        if *cursor < text.len() && text.as_bytes()[*cursor] == b';' {
            *cursor += 1;
        }
    }
    if *cursor >= text.len() || text.as_bytes()[*cursor] != b'$' {
        return Err(CxParseError::new(*cursor, "unterminated CX atom record"));
    }
    *cursor += 1;
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_atom_properties_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; each property write is an item checkpoint):
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
    // RDKit❗✔️: each complete nonempty property is retained at its source
    // write cursor; invalid nonnumeric bytes follow read_int's false branch
    // and one-byte source advance. The scan is linear with no prefix reparse.
    // Progress adds owned properties and one checkpoint per property write.
    *cursor += "atomProp:".len();
    if *cursor >= text.len() {
        return Err(CxParseError::new(*cursor, "empty CX atomProp record"));
    }
    let record_index = records.len();
    records.push(CxRecord::AtomProperties(Vec::new()));
    while *cursor < text.len() && !matches!(text.as_bytes()[*cursor], b'|' | b',') {
        let atom_offset = *cursor;
        let atom = match read_number(text, cursor) {
            Ok(atom) => Some(atom),
            Err(_) if *cursor == atom_offset => None,
            Err(error) => return Err(error),
        };
        if let Some(atom) = atom {
            expect_byte(text, cursor, b'.')?;
            let name = read_text_to(text, cursor, b".")?;
            if !name.is_empty() {
                expect_byte(text, cursor, b'.')?;
                let value = read_text_to(text, cursor, b":|,")?;
                if !value.is_empty() {
                    let property_index = match &mut records[record_index] {
                        CxRecord::AtomProperties(properties) => {
                            let property_index = properties.len();
                            properties.push(CxAtomProperty { atom, name, value });
                            property_index
                        }
                        _ => unreachable!("atomProp progress record changed kind"),
                    };
                    checkpoints.push(CxProgressCheckpoint {
                        record_index,
                        item_index: Some(property_index),
                        cursor: *cursor,
                        phase: CxProgressPhase::Item,
                    });
                }
            }
        }
        if *cursor < text.len() && !matches!(text.as_bytes()[*cursor], b'|' | b',') {
            *cursor += 1;
        }
    }
    match text.as_bytes().get(*cursor) {
        Some(&b'|') => {}
        Some(&b',') => *cursor += 1,
        _ => {
            return Err(CxParseError::new(
                *cursor,
                "unterminated CX atomProp record",
            ));
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_coordinate_bonds_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; the downstream query consumer applies each pair):
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
    // RDKit❗❌: pair scanning stays linear, while graph-neutral progress stores
    // each source pair and its checkpoint for ordered downstream mutation.
    let kind = match text.as_bytes().get(*cursor).copied() {
        Some(b'C') => CxCoordinateBondKind::Dative,
        Some(b'H') => CxCoordinateBondKind::Hydrogen,
        _ => return Err(CxParseError::new(*cursor, "invalid CX coordinate bond")),
    };
    *cursor += 1;
    expect_byte(text, cursor, b':')?;
    let record_index = records.len();
    records.push(CxRecord::CoordinateBonds(CxCoordinateBonds {
        kind,
        bonds: Vec::new(),
    }));
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let (atom, bond) = read_pair(text, cursor, b'.')?;
        let item_index = match &mut records[record_index] {
            CxRecord::CoordinateBonds(annotation) => {
                let item_index = annotation.bonds.len();
                annotation.bonds.push(CxBondReference { atom, bond });
                item_index
            }
            _ => unreachable!("coordinate bond progress record changed kind"),
        };
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_zero_bonds_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; each valid bond-type write is an item checkpoint):
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
    */
    // RDKit❗❌: the scanner is linear; item retention and ordered checkpoints
    // add per-index storage so the consumer can preserve partial writes.
    expect_byte(text, cursor, b'Z')?;
    expect_byte(text, cursor, b':')?;
    let record_index = records.len();
    records.push(CxRecord::ZeroBonds(Vec::new()));
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let bond = read_number(text, cursor)?;
        let item_index = match &mut records[record_index] {
            CxRecord::ZeroBonds(indices) => {
                let item_index = indices.len();
                indices.push(bond);
                item_index
            }
            _ => unreachable!("zero bond progress record changed kind"),
        };
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_unsaturation_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; atom mutation is Search-lowering-owned):
    /*
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
    // RDKit✔️❌: retain each parsed index and source cursor so Search can apply
    // its query effect before a later index fails; checkpoints add linear storage.
    if *cursor + 1 >= text.len() {
        return Err(CxParseError::new(*cursor, "unterminated CX index list"));
    }
    expect_byte(text, cursor, b'u')?;
    expect_byte(text, cursor, b':')?;
    let record_index = records.len();
    records.push(CxRecord::Unsaturation(Vec::new()));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let atom = read_number(text, cursor)?;
        let item_index = match &mut records[record_index] {
            CxRecord::Unsaturation(indices) => {
                let item_index = indices.len();
                indices.push(atom);
                item_index
            }
            _ => unreachable!("unsaturation progress record changed kind"),
        };
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_radicals_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; progress captures each source atom assignment):
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
    // RDKit❗❌: radical sections scan once in source order; detached radical
    // values and item checkpoints add per-assignment storage for partial writes.
    let record_index = records.len();
    records.push(CxRecord::Radicals(Vec::new()));
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

        let atom = read_number(text, cursor)?;
        let item_index = match &mut records[record_index] {
            CxRecord::Radicals(radicals) => {
                let item_index = radicals.len();
                radicals.push(CxRadical { atom, electrons });
                item_index
            }
            _ => unreachable!("radical progress record changed kind"),
        };
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });

        while *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
            if *cursor < text.len() && !text.as_bytes()[*cursor].is_ascii_digit() {
                break;
            }
            let atom = read_number(text, cursor)?;
            let item_index = match &mut records[record_index] {
                CxRecord::Radicals(radicals) => {
                    let item_index = radicals.len();
                    radicals.push(CxRadical { atom, electrons });
                    item_index
                }
                _ => unreachable!("radical progress record changed kind"),
            };
            checkpoints.push(CxProgressCheckpoint {
                record_index,
                item_index: Some(item_index),
                cursor: *cursor,
                phase: CxProgressPhase::Item,
            });
        }

        if *cursor >= text.len() {
            return Err(CxParseError::new(*cursor, "unterminated CX radical record"));
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_enhanced_stereo_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; Search lowering commits only after this helper):
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
    // RDKit❗❌: the CX parser preserves each ordered member on failure; one
    // progress checkpoint per parsed index adds linear storage so Search can
    // defer the group effect until the pinned helper's completion point.
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
    let record_index = records.len();
    records.push(CxRecord::EnhancedStereo(CxEnhancedStereo {
        kind,
        group_id,
        atoms: Vec::new(),
    }));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });
    while text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
        let atom = read_number(text, cursor)?;
        let item_index = match &mut records[record_index] {
            CxRecord::EnhancedStereo(stereo) => {
                let item_index = stereo.atoms.len();
                stereo.atoms.push(atom);
                item_index
            }
            _ => unreachable!("enhanced-stereo progress record changed kind"),
        };
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_ring_bonds_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; query construction is Search-lowering-owned):
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
    // RDKit✔️❌: retain complete ordered pairs on later parse failure; each pair
    // checkpoint adds linear storage before Search applies its query effect.
    if text.as_bytes().get(*cursor..*cursor + 3) != Some(b"rb:") {
        return Err(CxParseError::new(*cursor, "invalid CX ring-bond record"));
    }
    *cursor += 3;
    let record_index = records.len();
    records.push(CxRecord::RingBonds(Vec::new()));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let atom = read_number(text, cursor)?;
        if (*cursor).saturating_add(1) >= text.len() {
            return Err(CxParseError::new(
                *cursor,
                "expected ':' before CX ring-bond count",
            ));
        }
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
        let item_index = match &mut records[record_index] {
            CxRecord::RingBonds(constraints) => {
                let item_index = constraints.len();
                constraints.push(CxRingBond { atom, constraint });
                item_index
            }
            _ => unreachable!("ring-bond progress record changed kind"),
        };
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_substitution_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; query construction is Search-lowering-owned):
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
    // RDKit✔️❌: retain only complete atom/value pairs before failure; each pair
    // checkpoint adds linear storage before Search applies its query effect.
    if text.as_bytes().get(*cursor..*cursor + 2) != Some(b"s:") {
        return Err(CxParseError::new(*cursor, "invalid CX substitution record"));
    }
    *cursor += 1;
    expect_byte(text, cursor, b':')?;
    let record_index = records.len();
    records.push(CxRecord::Substitution(Vec::new()));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let atom = read_number(text, cursor)?;
        if (*cursor).saturating_add(1) >= text.len() {
            return Err(CxParseError::new(
                *cursor,
                "expected ':' before CX substitution count",
            ));
        }
        expect_byte(text, cursor, b':')?;
        let constraint = if text.as_bytes().get(*cursor) == Some(&b'*') {
            *cursor += 1;
            CxCountConstraint::QueryScan
        } else {
            CxCountConstraint::Exact(read_number(text, cursor)? as u32)
        };
        let item_index = match &mut records[record_index] {
            CxRecord::Substitution(constraints) => {
                let item_index = constraints.len();
                constraints.push(CxAtomConstraint { atom, constraint });
                item_index
            }
            _ => unreachable!("substitution progress record changed kind"),
        };
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_link_nodes_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
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
    // RDKit❗❌: retain each parsed item and source cursor, but defer the
    // whole-helper destination effect until Complete; per-item detached
    // records and checkpoints add allocation beyond the source accumulator.
    if text.as_bytes().get(*cursor..*cursor + 3) != Some(b"LN:") {
        return Err(CxParseError::new(*cursor, "invalid CX link-node record"));
    }
    *cursor += 3;
    let record_index = records.len();
    records.push(CxRecord::LinkNodes(Vec::new()));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });

    let bytes = text.as_bytes();
    while bytes.get(*cursor).is_some_and(u8::is_ascii_digit) {
        let atom = read_number(text, cursor)?;
        if cursor.saturating_add(1) >= bytes.len() || bytes.get(*cursor) != Some(&b':') {
            return Err(CxParseError::new(*cursor, "expected ':' in CX link node"));
        }
        *cursor += 1;
        let start_repetitions = read_number(text, cursor)?;
        if cursor.saturating_add(1) >= bytes.len() || bytes.get(*cursor) != Some(&b'.') {
            return Err(CxParseError::new(*cursor, "expected '.' in CX link node"));
        }
        *cursor += 1;
        let end_repetitions = read_number(text, cursor)?;
        let outer_atoms = if text.as_bytes().get(*cursor) == Some(&b'.') {
            *cursor += 1;
            let first = read_number(text, cursor)?;
            // RDKit increments past one separator byte here without checking
            // its value before reading the second explicit outer atom.
            if *cursor < bytes.len() {
                *cursor += 1;
            }
            Some([first, read_number(text, cursor)?])
        } else {
            None
        };

        if bytes.get(*cursor) == Some(&b',') {
            *cursor += 1;
        }

        let item_index = match &mut records[record_index] {
            CxRecord::LinkNodes(nodes) => {
                let item_index = nodes.len();
                nodes.push(CxLinkNode {
                    atom,
                    start_repetitions,
                    end_repetitions,
                    outer_atoms,
                });
                item_index
            }
            _ => unreachable!("link-node progress record changed kind"),
        };
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_data_sgroup_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
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
    // RDKit❗✔️: retain each source field mutation and cursor on failure while
    // deferring the destination SGroup effect until this helper completes.
    if text.as_bytes().get(*cursor..*cursor + 4) != Some(b"SgD:") {
        return Err(CxParseError::new(*cursor, "invalid CX data SGroup record"));
    }
    *cursor += 4;
    let atoms = parse_delimited_number_list(text, cursor, b',')?;
    let record_index = records.len();
    records.push(CxRecord::DataSGroup(CxDataSGroup {
        atoms,
        field_name: String::new(),
        data: String::new(),
        query_op: String::new(),
        field_info: String::new(),
        field_tag: String::new(),
        coordinates: None,
    }));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });

    // RDKit unconditionally advances one byte after read_int_list rather than
    // validating the separator. Keep that cursor mutation source-shaped.
    if *cursor < text.len() {
        *cursor += 1;
    }
    for (item_index, field) in [
        "FIELDNAME",
        "DATAFIELDS",
        "QUERYOP",
        "FIELDINFO",
        "FIELDTAG",
    ]
    .into_iter()
    .enumerate()
    {
        let value = read_data_sgroup_attr(text, cursor)?;
        if let Some(value) = value {
            let CxRecord::DataSGroup(group) = &mut records[record_index] else {
                unreachable!("data SGroup progress record changed kind");
            };
            match field {
                "FIELDNAME" => group.field_name = value,
                "DATAFIELDS" => group.data = value,
                "QUERYOP" => group.query_op = value,
                "FIELDINFO" => group.field_info = value,
                "FIELDTAG" => group.field_tag = value,
                _ => unreachable!("data SGroup field order changed"),
            }
        }
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
    }

    if text.as_bytes().get(*cursor) == Some(&b'(') {
        let coordinates = read_text_to(text, cursor, b")")?;
        // The source increments after read_text_to without validating that a
        // closing parenthesis was found. Cap the end-iterator case safely.
        if *cursor < text.len() {
            *cursor += 1;
        }
        let CxRecord::DataSGroup(group) = &mut records[record_index] else {
            unreachable!("data SGroup progress record changed kind");
        };
        group.coordinates = Some(coordinates);
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(5),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
    }

    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn read_data_sgroup_attr(text: &str, cursor: &mut usize) -> Result<Option<String>, CxParseError> {
    // RDKit source (verbatim; delimiter and advancement from
    // parse_data_sgroup_attr):
    /*
      if (first != last && *first != '|') {
        std::string data = read_text_to(first, last, ":");
        ++first;
    */
    if *cursor >= text.len() || text.as_bytes().get(*cursor) == Some(&b'|') {
        return Ok(None);
    }
    let value = read_text_to(text, cursor, b":")?;
    if *cursor < text.len() {
        *cursor += 1;
    }
    Ok(Some(value))
}

fn parse_sgroup_hierarchy_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
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
    // RDKit❗✔️: retain incomplete parent relationships and parsed child values;
    // emit child checkpoints only after the entire source child list parses,
    // when the pinned helper begins its mutation loop.
    if text.as_bytes().get(*cursor..*cursor + 4) != Some(b"SgH:") {
        return Err(CxParseError::new(
            *cursor,
            "invalid CX SGroup hierarchy record",
        ));
    }
    *cursor += 4;
    let record_index = records.len();
    records.push(CxRecord::SGroupHierarchy(Vec::new()));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });
    let mut next_item_index = 0;
    loop {
        let parent = read_number(text, cursor)?;
        let hierarchy_index = match &mut records[record_index] {
            CxRecord::SGroupHierarchy(hierarchies) => {
                let hierarchy_index = hierarchies.len();
                hierarchies.push(CxSGroupHierarchy {
                    parent,
                    children: Vec::new(),
                });
                hierarchy_index
            }
            _ => unreachable!("SGroup hierarchy progress record changed kind"),
        };
        expect_byte(text, cursor, b':')?;
        let first_item_index = next_item_index;
        loop {
            if text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
                let child = read_number(text, cursor)?;
                let CxRecord::SGroupHierarchy(hierarchies) = &mut records[record_index] else {
                    unreachable!("SGroup hierarchy progress record changed kind");
                };
                hierarchies[hierarchy_index].children.push(child);
                next_item_index += 1;
            }
            if text.as_bytes().get(*cursor) != Some(&b'.') {
                break;
            }
            *cursor += 1;
        }

        for item_index in first_item_index..next_item_index {
            checkpoints.push(CxProgressCheckpoint {
                record_index,
                item_index: Some(item_index),
                cursor: *cursor,
                phase: CxProgressPhase::Item,
            });
        }

        if text.as_bytes().get(*cursor) == Some(&b',') {
            *cursor += 1;
        } else {
            break;
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

#[derive(Clone, Copy)]
enum PolymerSGroupIndexField {
    Atoms,
    HeadCrossings,
    TailCrossings,
}

fn parse_polymer_index_list_progress(
    text: &str,
    cursor: &mut usize,
    separator: u8,
    field: PolymerSGroupIndexField,
    record_index: usize,
    next_item_index: &mut usize,
    records: &mut [CxRecord],
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; parser::read_int and parser::read_int_list):
    // RDKit✔️🔝: template <typename Iterator>
    // RDKit✔️🔝: bool read_int(Iterator &first, Iterator last, unsigned int &res) {
    // RDKit✔️🔝:   std::string num = "";
    // RDKit✔️🔝:   while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️🔝:     num += *first;
    // RDKit✔️🔝:     ++first;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   if (num.empty()) {
    // RDKit✔️🔝:     return false;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   res = boost::lexical_cast<unsigned int>(num);
    // RDKit✔️🔝:   return true;
    // RDKit✔️🔝: }
    // RDKit✔️🔝: template <typename Iterator>
    // RDKit✔️🔝: bool read_int_list(Iterator &first, Iterator last,
    // RDKit✔️🔝:                    std::vector<unsigned int> &res, char sep = ',') {
    // RDKit✔️🔝:   while (1) {
    // RDKit✔️🔝:     std::string num = "";
    // RDKit✔️🔝:     while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️🔝:       num += *first;
    // RDKit✔️🔝:       ++first;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     if (!num.empty()) {
    // RDKit✔️🔝:       res.push_back(boost::lexical_cast<unsigned int>(num));
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     if (first >= last || *first != sep) {
    // RDKit✔️🔝:       break;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     ++first;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return true;
    // RDKit✔️🔝: }
    // Each parsed number is retained in its partial CX record and receives an
    // ordered checkpoint; the source list accepts empty slots and trailing separators.
    loop {
        if text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
            let value = read_number(text, cursor)?;
            let CxRecord::PolymerSGroup(polymer) = &mut records[record_index] else {
                unreachable!("polymer progress record changed kind");
            };
            match field {
                PolymerSGroupIndexField::Atoms => polymer.atoms.push(value),
                PolymerSGroupIndexField::HeadCrossings => polymer.head_crossings.push(value),
                PolymerSGroupIndexField::TailCrossings => polymer.tail_crossings.push(value),
            }
            checkpoints.push(CxProgressCheckpoint {
                record_index,
                item_index: Some(*next_item_index),
                cursor: *cursor,
                phase: CxProgressPhase::Item,
            });
            *next_item_index += 1;
        }
        if text.as_bytes().get(*cursor) != Some(&separator) {
            break;
        }
        *cursor += 1;
    }
    Ok(())
}

fn parse_polymer_sgroup_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; type map and polymer helper from CXSmilesOps.cpp):
    // RDKit✔️✔️: const std::map<std::string, std::string> sgroupTypemap = {
    // RDKit✔️✔️:     {"n", "SRU"},   {"mon", "MON"}, {"mer", "MER"}, {"co", "COP"},
    // RDKit✔️✔️:     {"xl", "CRO"},  {"mod", "MOD"}, {"mix", "MIX"}, {"f", "FOR"},
    // RDKit✔️✔️:     {"any", "ANY"}, {"gen", "GEN"}, {"c", "COM"},   {"grf", "GRA"},
    // RDKit✔️✔️:     {"alt", "COP"}, {"ran", "COP"}, {"blk", "COP"}};
    // RDKit❗❌: template <typename Iterator>
    // RDKit❗❌: bool parse_polymer_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit❗❌:                           unsigned int startAtomIdx, unsigned int nSGroups) {
    // RDKit❗❌:   if (first >= last || *first != 'S' || first + 2 >= last ||
    // RDKit❗❌:       *(first + 1) != 'g' || *(first + 2) != ':') {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:   first += 3;
    // RDKit❗❌:   const auto type_code = read_text_to(first, last, ":");
    // RDKit❗❌:   ++first;
    // RDKit❗❌:   const auto type = sgroupTypemap.find(type_code);
    // RDKit❗❌:   if (type == sgroupTypemap.end()) {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<unsigned int> atoms;
    // RDKit❗❌:   if (!read_int_list(first, last, atoms)) {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (first <= last && *first == ':') {
    // RDKit❗❌:     ++first;
    // RDKit❗❌:     std::string subscript = read_text_to(first, last, ":|");
    // RDKit❗❌:     if (keepSGroup && !subscript.empty()) {
    // RDKit❗❌:       sgroup.setProp("LABEL", subscript);
    // RDKit❗❌:     }
    // RDKit❗❌:     if (first <= last && *first == ':') {
    // RDKit❗❌:       ++first;
    // RDKit❗❌:       std::string superscript = read_text_to(first, last, ":|,");
    // RDKit❗❌:       if (keepSGroup && !superscript.empty()) {
    // RDKit❗❌:         sgroup.setProp("CONNECT", superscript);
    // RDKit❗❌:       }
    // RDKit❗❌:       if (first <= last && *first == ':') {
    // RDKit❗❌:         ++first;
    // RDKit❗❌:         if (!read_int_list(first, last, headCrossing)) {
    // RDKit❗❌:           return false;
    // RDKit❗❌:         }
    // RDKit❗❌:         if (first <= last && *first == ':') {
    // RDKit❗❌:           ++first;
    // RDKit❗❌:           if (!read_int_list(first, last, tailCrossing)) {
    // RDKit❗❌:             return false;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (keepSGroup) {
    // RDKit❗❌:     processCXSmilesLabels(mol);
    // RDKit❗❌:     finalizePolymerSGroup(mol, sgroup);
    // RDKit❗❌:     sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
    // RDKit❗❌:     addSubstanceGroup(mol, sgroup);
    // RDKit❗❌:   }
    // RDKit❗❌:   return true;
    // RDKit❗❌: }
    // The progress record retains the source's partial local SGroup state; the
    // Search owner installs it only at Complete, matching the source helper's
    // single destination commit point.
    if text.as_bytes().get(*cursor..*cursor + 3) != Some(b"Sg:") {
        return Err(CxParseError::new(
            *cursor,
            "invalid CX polymer SGroup record",
        ));
    }
    *cursor += 3;
    let type_code = read_text_to(text, cursor, b":")?;
    if *cursor < text.len() {
        *cursor += 1;
    }
    if !matches!(
        type_code.as_str(),
        "n" | "mon"
            | "mer"
            | "co"
            | "xl"
            | "mod"
            | "mix"
            | "f"
            | "any"
            | "gen"
            | "c"
            | "grf"
            | "alt"
            | "ran"
            | "blk"
    ) {
        return Err(CxParseError::new(*cursor, "unknown CX polymer SGroup type"));
    }

    let record_index = records.len();
    records.push(CxRecord::PolymerSGroup(CxPolymerSGroup {
        type_code,
        atoms: Vec::new(),
        label: String::new(),
        connect: String::new(),
        head_crossings: Vec::new(),
        tail_crossings: Vec::new(),
    }));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });
    let mut next_item_index = 0;
    parse_polymer_index_list_progress(
        text,
        cursor,
        b',',
        PolymerSGroupIndexField::Atoms,
        record_index,
        &mut next_item_index,
        records,
        checkpoints,
    )?;

    if text.as_bytes().get(*cursor) == Some(&b':') {
        *cursor += 1;
        let label = read_text_to(text, cursor, b":|")?;
        if !label.is_empty() {
            let CxRecord::PolymerSGroup(polymer) = &mut records[record_index] else {
                unreachable!("polymer progress record changed kind");
            };
            polymer.label = label;
            checkpoints.push(CxProgressCheckpoint {
                record_index,
                item_index: Some(next_item_index),
                cursor: *cursor,
                phase: CxProgressPhase::Item,
            });
            next_item_index += 1;
        }
        if text.as_bytes().get(*cursor) == Some(&b':') {
            *cursor += 1;
            let connect = read_text_to(text, cursor, b":|,")?;
            if !connect.is_empty() {
                let CxRecord::PolymerSGroup(polymer) = &mut records[record_index] else {
                    unreachable!("polymer progress record changed kind");
                };
                polymer.connect = connect;
                checkpoints.push(CxProgressCheckpoint {
                    record_index,
                    item_index: Some(next_item_index),
                    cursor: *cursor,
                    phase: CxProgressPhase::Item,
                });
                next_item_index += 1;
            }
            if text.as_bytes().get(*cursor) == Some(&b':') {
                *cursor += 1;
                parse_polymer_index_list_progress(
                    text,
                    cursor,
                    b',',
                    PolymerSGroupIndexField::HeadCrossings,
                    record_index,
                    &mut next_item_index,
                    records,
                    checkpoints,
                )?;
                if text.as_bytes().get(*cursor) == Some(&b':') {
                    *cursor += 1;
                    parse_polymer_index_list_progress(
                        text,
                        cursor,
                        b',',
                        PolymerSGroupIndexField::TailCrossings,
                        record_index,
                        &mut next_item_index,
                        records,
                        checkpoints,
                    )?;
                }
            }
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_variable_attachments_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
    // RDKit source (verbatim; CXSmilesOps.cpp::parse_variable_attachments):
    // RDKit❗❌: bool parse_variable_attachments(Iterator &first, Iterator last,
    // RDKit❗❌:                                 RDKit::RWMol &mol, unsigned int startAtomIdx) {
    // RDKit❗❌:   // these look like: CO*.C1=CC=NC=C1 |m:2:3.5.4|
    // RDKit❗❌:   // that corresponds to replacing the bond to atom 2 with bonds to atom 3, 5,
    // RDKit❗❌:   // or 4
    // RDKit❗❌:   //
    // RDKit❗❌:   if (first >= last || *first != 'm' || first + 1 >= last ||
    // RDKit❗❌:       *(first + 1) != ':') {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:   first += 2;
    // RDKit❗❌:
    // RDKit❗❌:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit❗❌:     unsigned int at1idx;
    // RDKit❗❌:     if (!read_int(first, last, at1idx)) {
    // RDKit❗❌:       return false;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (VALID_ATIDX(at1idx) &&
    // RDKit❗❌:         mol.getAtomWithIdx(at1idx - startAtomIdx)->getDegree() != 1) {
    // RDKit❗❌:       BOOST_LOG(rdWarningLog)
    // RDKit❗❌:           << "position variation bond to atom with more than one bond"
    // RDKit❗❌:           << std::endl;
    // RDKit❗❌:       return false;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (first < last && *first == ':') {
    // RDKit❗❌:       ++first;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       BOOST_LOG(rdWarningLog) << "improperly formatted m: block" << std::endl;
    // RDKit❗❌:       return false;
    // RDKit❗❌:     }
    // RDKit❗❌:     std::vector<std::string> others;
    // RDKit❗❌:     while (first < last && *first >= '0' && *first <= '9') {
    // RDKit❗❌:       unsigned int aidx;
    // RDKit❗❌:       if (!read_int(first, last, aidx)) {
    // RDKit❗❌:         return false;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (VALID_ATIDX(aidx)) {
    // RDKit❗❌:         others.push_back(std::to_string(aidx - startAtomIdx + 1));
    // RDKit❗❌:       }
    // RDKit❗❌:       if (first < last && *first == '.') {
    // RDKit❗❌:         ++first;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (VALID_ATIDX(at1idx)) {
    // RDKit❗❌:       std::string endPts = "(" + std::to_string(others.size());
    // RDKit❗❌:       for (auto idx : others) {
    // RDKit❗❌:         endPts += " " + idx;
    // RDKit❗❌:       }
    // RDKit❗❌:       endPts += ")";
    // RDKit❗❌:
    // RDKit❗❌:       for (auto nbri : boost::make_iterator_range(
    // RDKit❗❌:                mol.getAtomBonds(mol.getAtomWithIdx(at1idx - startAtomIdx)))) {
    // RDKit❗❌:         auto bnd = mol[nbri];
    // RDKit❗❌:         bnd->setProp(common_properties::_MolFileBondEndPts, endPts);
    // RDKit❗❌:         bnd->setProp(common_properties::_MolFileBondAttach, std::string("ANY"));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (first < last && *first == ',') {
    // RDKit❗❌:       ++first;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return true;
    // RDKit❗❌: }
    // RDKit❗❌: Primary atom checkpoints precede colon/degree-sensitive source
    // behavior; completed-row checkpoints precede source comma consumption.
    if text.as_bytes().get(*cursor..*cursor + 2) != Some(b"m:") {
        return Err(CxParseError::new(
            *cursor,
            "invalid CX variable-attachment record",
        ));
    }
    *cursor += 2;
    let record_index = records.len();
    records.push(CxRecord::VariableAttachments(Vec::new()));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });
    let mut next_item_index = 0;
    while text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
        let atom = read_number(text, cursor)?;
        let attachment_index = match &mut records[record_index] {
            CxRecord::VariableAttachments(attachments) => {
                let attachment_index = attachments.len();
                attachments.push(CxVariableAttachment {
                    atom,
                    endpoints: Vec::new(),
                });
                attachment_index
            }
            _ => unreachable!("variable-attachment progress record changed kind"),
        };
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(next_item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        next_item_index += 1;

        if text.as_bytes().get(*cursor) == Some(&b':') {
            *cursor += 1;
        } else {
            return Err(CxParseError::new(*cursor, "improperly formatted m: block"));
        }
        while text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
            let endpoint = read_number(text, cursor)?;
            let CxRecord::VariableAttachments(attachments) = &mut records[record_index] else {
                unreachable!("variable-attachment progress record changed kind");
            };
            attachments[attachment_index].endpoints.push(endpoint);
            if text.as_bytes().get(*cursor) == Some(&b'.') {
                *cursor += 1;
            }
        }
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(next_item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        next_item_index += 1;
        if text.as_bytes().get(*cursor) == Some(&b',') {
            *cursor += 1;
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_wedge_bonds_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
) -> Result<(), CxParseError> {
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
    // RDKit❗❌: parser checkpoints add one record and one row marker while
    // retaining the pinned cursor and per-pair mutation boundary.
    let bytes = text.as_bytes();
    let last = bytes
        .iter()
        .enumerate()
        .skip(*cursor)
        .find_map(|(index, byte)| (*byte == b'|').then_some(index))
        .unwrap_or(bytes.len());
    if *cursor >= last || bytes.get(*cursor) != Some(&b'w') || *cursor + 1 >= last {
        return Err(CxParseError::new(*cursor, "invalid CX wedge record"));
    }
    *cursor += 1;
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
    if *cursor >= last || *cursor + 1 >= last || bytes.get(*cursor) != Some(&b':') {
        return Err(CxParseError::new(*cursor, "invalid CX wedge record"));
    }
    *cursor += 1;
    let record_index = records.len();
    records.push(CxRecord::WedgedBonds(Vec::new()));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });
    let mut item_index = 0;
    while *cursor < last && bytes[*cursor].is_ascii_digit() {
        let (atom, bond) = read_pair(text, cursor, b'.')?;
        let CxRecord::WedgedBonds(wedges) = &mut records[record_index] else {
            unreachable!("wedge progress record changed kind");
        };
        wedges.push(CxWedgeBond {
            atom,
            bond,
            direction,
            configuration,
        });
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        item_index += 1;
        if *cursor < last && bytes[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

fn parse_double_bond_stereo_progress(
    text: &str,
    cursor: &mut usize,
    records: &mut Vec<CxRecord>,
    checkpoints: &mut Vec<CxProgressCheckpoint>,
    stereo: CxDoubleBondStereoKind,
) -> Result<(), CxParseError> {
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
    // RDKit❗❌: item checkpoints are emitted at each parsed source bond index,
    // before source comma consumption, so earlier target effects survive errors.
    let bytes = text.as_bytes();
    let last = bytes
        .iter()
        .enumerate()
        .skip(*cursor)
        .find_map(|(index, byte)| (*byte == b'|').then_some(index))
        .unwrap_or(bytes.len());
    while *cursor < last && bytes[*cursor] != b':' {
        *cursor += 1;
    }
    if *cursor >= last {
        return Err(CxParseError::new(
            *cursor,
            "invalid CX double-bond stereo record",
        ));
    }
    *cursor += 1;
    let record_index = records.len();
    records.push(CxRecord::DoubleBondStereo(CxDoubleBondStereo {
        stereo,
        bonds: Vec::new(),
    }));
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Begin,
    });
    let mut item_index = 0;
    while *cursor < last && bytes[*cursor].is_ascii_digit() {
        let bond = read_number(text, cursor)?;
        let CxRecord::DoubleBondStereo(parsed) = &mut records[record_index] else {
            unreachable!("double-bond stereo progress record changed kind");
        };
        parsed.bonds.push(bond);
        checkpoints.push(CxProgressCheckpoint {
            record_index,
            item_index: Some(item_index),
            cursor: *cursor,
            phase: CxProgressPhase::Item,
        });
        item_index += 1;
        if *cursor < last && bytes[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    checkpoints.push(CxProgressCheckpoint {
        record_index,
        item_index: None,
        cursor: *cursor,
        phase: CxProgressPhase::Complete,
    });
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        CxAtomConstraint, CxAtomProperty, CxBondReference, CxCoordinateBondKind, CxCountConstraint,
        CxDispatch, CxDoubleBondStereoKind, CxLinkNode, CxProgressPhase, CxRecord, CxRingBond,
        CxSGroupHierarchy, CxStereoGroupKind, CxVariableAttachment, CxWedgeDirection,
        classify_cx_dispatch, parse_cx_extensions, parse_cx_extensions_progress,
    };

    #[test]
    fn cx_progress_dispatch_uses_source_order_and_length_guards() {
        for (text, expected) in [
            (b"(".as_slice(), CxDispatch::Coordinates),
            (b"$", CxDispatch::AtomLabelsOrValues),
            (b"$_AV:", CxDispatch::AtomLabelsOrValues),
            (b"atomProp:|", CxDispatch::AtomProperties),
            (b"atomProp:", CxDispatch::EnhancedStereo),
            (b"C", CxDispatch::CoordinateBonds),
            (b"H", CxDispatch::CoordinateBonds),
            (b"Z", CxDispatch::ZeroBonds),
            (b"^", CxDispatch::Radicals),
            (b"a", CxDispatch::EnhancedStereo),
            (b"o", CxDispatch::EnhancedStereo),
            (b"&x", CxDispatch::EnhancedStereo),
            (b"&#", CxDispatch::Unknown),
            (b"rb", CxDispatch::RingBonds),
            (b"LN", CxDispatch::LinkNodes),
            (b"SgD", CxDispatch::DataSGroup),
            (b"SgH", CxDispatch::SGroupHierarchy),
            (b"Sg", CxDispatch::PolymerSGroup),
            (b"u", CxDispatch::Unsaturation),
            (b"s", CxDispatch::Substitution),
            (b"m", CxDispatch::VariableAttachments),
            (b"w", CxDispatch::WedgedBonds),
            (b"ctu", CxDispatch::DoubleBondAny),
            (b"c", CxDispatch::DoubleBondCis),
            (b"t", CxDispatch::DoubleBondTrans),
            ("☃".as_bytes(), CxDispatch::Unknown),
        ] {
            assert_eq!(classify_cx_dispatch(text, 0), expected, "{text:?}");
        }

        let progress = parse_cx_extensions_progress("|x$label$|");
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), "|x$label$|".len());
        assert!(matches!(progress.records()[0], CxRecord::Unknown(_)));
        assert!(matches!(progress.records()[1], CxRecord::AtomLabels(_)));
    }

    #[test]
    fn cx_progress_dispatch_reports_each_branch_failure_cursor() {
        for (text, expected_cursor) in [
            ("|(", 2),
            ("|$", 2),
            ("|atomProp:0x.foo|", 11),
            ("|C?", 2),
            ("|H?", 2),
            ("|Z?", 2),
            ("|^x", 2),
            ("|a?", 2),
            ("|o?", 2),
            ("|&?", 2),
            ("|rb?", 1),
            ("|LN?", 1),
            ("|SgD?", 1),
            ("|SgH?", 1),
            ("|Sg?", 1),
            ("|u?", 2),
            ("|s?", 1),
            ("|m?", 1),
            ("|w?", 2),
            ("|ctu?", 5),
            ("|c?", 3),
            ("|t?", 3),
            ("|?", 2),
        ] {
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert_eq!(progress.consumed(), expected_cursor, "{text:?}");
            assert!(progress.error().is_some(), "{text:?}");
        }
    }

    #[test]
    fn cx_progress_dispatch_keeps_prior_effects_and_utf8_cursor_distinct_from_error() {
        let text = "|☃,$label$ C?";
        let progress = parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), text.find('?').expect("failure marker"));
        assert_eq!(progress.records().len(), 3);
        assert!(matches!(&progress.records()[0], CxRecord::Unknown(raw) if raw == "☃,"));
        assert!(matches!(progress.records()[1], CxRecord::AtomLabels(_)));
        assert!(matches!(&progress.records()[2], CxRecord::Unknown(raw) if raw == " "));
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[0].record_index, 1);
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Item);
        assert_eq!(progress.checkpoints()[1].record_index, 1);
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Complete);

        let malformed_coordinate = parse_cx_extensions_progress("|(no)|");
        assert_eq!(malformed_coordinate.consumed(), 4);
        assert_eq!(
            malformed_coordinate.error().map(|error| error.offset),
            Some(2)
        );
        assert_ne!(
            malformed_coordinate.consumed(),
            malformed_coordinate
                .error()
                .expect("coordinate error")
                .offset
        );
    }

    #[test]
    fn cx_progress_coordinates_commits_rows_at_source_cursors() {
        let text = "|(1,,2;;3,4,,9;)|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        let CxRecord::Coordinates(coordinates) = &progress.records()[0] else {
            panic!("coordinate record expected");
        };
        assert_eq!(
            coordinates.values,
            vec![Some([1.0, 0.0, 2.0]), None, Some([3.0, 4.0, 0.0])]
        );
        assert!(coordinates.is_3d);

        let separators = text
            .match_indices(';')
            .map(|(cursor, _)| cursor)
            .collect::<Vec<_>>();
        let checkpoints = progress.checkpoints();
        assert_eq!(checkpoints.len(), 5);
        assert_eq!(checkpoints[0].phase, CxProgressPhase::Begin);
        assert_eq!(checkpoints[0].cursor, 1);
        for (row, separator) in separators.into_iter().enumerate() {
            assert_eq!(checkpoints[row + 1].phase, CxProgressPhase::Item);
            assert_eq!(checkpoints[row + 1].item_index, Some(row));
            assert_eq!(checkpoints[row + 1].cursor, separator);
        }
        assert_eq!(checkpoints[4].phase, CxProgressPhase::Complete);
        assert_eq!(checkpoints[4].cursor, text.find(')').unwrap() + 1);
    }

    #[test]
    fn cx_progress_coordinates_reports_each_component_failure_and_keeps_prior_rows() {
        for text in ["|(bad,0,0)|", "|(1,bad,0)|", "|(1,0,bad)|"] {
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text}");
            assert_eq!(progress.consumed(), text.find(')').unwrap(), "{text}");
            assert_eq!(
                progress.error().map(|error| error.offset),
                Some(2),
                "{text}"
            );
            let CxRecord::Coordinates(coordinates) = &progress.records()[0] else {
                panic!("coordinate record expected");
            };
            assert!(coordinates.values.is_empty(), "{text}");
            assert!(coordinates.is_3d, "conversion failure keeps source default");
            assert_eq!(progress.checkpoints().len(), 1, "{text}");
            assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Begin);
        }

        let text = "|(0,0,0.0011;1,bad,2)|";
        let progress = parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), text.find(')').unwrap());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(text.find("1,bad").unwrap())
        );
        let CxRecord::Coordinates(coordinates) = &progress.records()[0] else {
            panic!("coordinate record expected");
        };
        assert_eq!(coordinates.values, vec![Some([0.0, 0.0, 0.0011])]);
        assert!(coordinates.is_3d);
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Item);
        assert_eq!(
            progress.checkpoints()[1].cursor,
            text.find(';').expect("row separator")
        );
    }

    #[test]
    fn cx_progress_coordinates_finalizes_z_before_a_missing_close() {
        for (text, expected_is_3d) in [("|(0,0,0.001)|", false), ("|(0,0,0.0011)|", true)] {
            let progress = parse_cx_extensions_progress(text);
            assert!(progress.is_complete(), "{text}");
            let CxRecord::Coordinates(coordinates) = &progress.records()[0] else {
                panic!("coordinate record expected");
            };
            assert_eq!(coordinates.is_3d, expected_is_3d, "{text}");
        }

        let text = "|(0,0,0.0011";
        let progress = parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        let CxRecord::Coordinates(coordinates) = &progress.records()[0] else {
            panic!("coordinate record expected");
        };
        assert!(coordinates.is_3d);
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[1].cursor, text.len());
    }

    #[test]
    fn cx_progress_coordinates_retains_earlier_record_and_prior_row_effects() {
        let text = "|$label$(0,0,1;bad,0,0)|";
        let progress = parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), text.find(')').unwrap());
        assert!(matches!(progress.records()[0], CxRecord::AtomLabels(_)));
        let CxRecord::Coordinates(coordinates) = &progress.records()[1] else {
            panic!("partial coordinate record expected");
        };
        assert_eq!(coordinates.values, vec![Some([0.0, 0.0, 1.0])]);
        assert_eq!(progress.checkpoints().len(), 4);
        assert_eq!(progress.checkpoints()[0].record_index, 0);
        assert_eq!(progress.checkpoints()[0].item_index, Some(0));
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Item);
        assert_eq!(progress.checkpoints()[1].record_index, 0);
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Complete);
        assert_eq!(progress.checkpoints()[2].record_index, 1);
        assert_eq!(progress.checkpoints()[2].phase, CxProgressPhase::Begin);
        assert_eq!(progress.checkpoints()[3].record_index, 1);
        assert_eq!(progress.checkpoints()[3].item_index, Some(0));
        assert_eq!(progress.checkpoints()[3].phase, CxProgressPhase::Item);
    }

    #[test]
    fn cx_progress_labels_and_values_keep_empty_slots_and_decode_entities() {
        let label_text = "|$a&#321;b;;line&#59;break$|";
        let labels = parse_cx_extensions_progress(label_text);
        assert!(labels.is_complete());
        let CxRecord::AtomLabels(values) = &labels.records()[0] else {
            panic!("atom labels expected");
        };
        assert_eq!(
            values,
            &vec![Some("aAb".to_owned()), None, Some("line;break".to_owned())]
        );
        assert_eq!(labels.checkpoints().len(), 3);
        assert_eq!(labels.checkpoints()[0].phase, CxProgressPhase::Item);
        assert_eq!(labels.checkpoints()[0].item_index, Some(0));
        assert_eq!(
            labels.checkpoints()[0].cursor,
            label_text.find(";;").unwrap()
        );
        assert_eq!(labels.checkpoints()[1].phase, CxProgressPhase::Item);
        assert_eq!(labels.checkpoints()[1].item_index, Some(2));
        assert_eq!(labels.checkpoints()[2].phase, CxProgressPhase::Complete);

        let value_text = "|$_AV:value0;;value2$|";
        let parsed_values = parse_cx_extensions_progress(value_text);
        assert!(parsed_values.is_complete());
        let CxRecord::AtomValues(values) = &parsed_values.records()[0] else {
            panic!("atom values expected");
        };
        assert_eq!(
            values,
            &vec![Some("value0".to_owned()), None, Some("value2".to_owned())]
        );
        assert_eq!(parsed_values.checkpoints()[0].item_index, Some(0));
        assert_eq!(parsed_values.checkpoints()[1].item_index, Some(2));

        let empty = parse_cx_extensions_progress("|$$|");
        assert!(empty.is_complete());
        assert!(matches!(&empty.records()[0], CxRecord::AtomLabels(values) if values.is_empty()));
        assert_eq!(empty.checkpoints().len(), 1);
        assert_eq!(empty.checkpoints()[0].phase, CxProgressPhase::Complete);
    }

    #[test]
    fn cx_progress_labels_and_values_retain_slots_before_missing_closers() {
        for (text, expected_value) in [("|$first;second", false), ("|$_AV:first;second", true)] {
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text}");
            assert_eq!(progress.consumed(), text.len(), "{text}");
            assert_eq!(progress.error().map(|error| error.offset), Some(text.len()));
            let values = match (&progress.records()[0], expected_value) {
                (CxRecord::AtomLabels(values), false) | (CxRecord::AtomValues(values), true) => {
                    values
                }
                other => panic!("unexpected CX slot record: {other:?}"),
            };
            assert_eq!(
                values,
                &vec![Some("first".to_owned()), Some("second".to_owned())]
            );
            assert_eq!(progress.checkpoints().len(), 2);
            assert_eq!(progress.checkpoints()[0].item_index, Some(0));
            assert_eq!(
                progress.checkpoints()[0].cursor,
                text.find(';').expect("first slot separator")
            );
            assert_eq!(progress.checkpoints()[1].item_index, Some(1));
            assert_eq!(progress.checkpoints()[1].cursor, text.len());
        }
    }

    #[test]
    fn cx_progress_labels_keeps_prior_items_when_entity_decoding_fails() {
        for (text, entity) in [
            ("|$first;bad&#oops$|", "&#oops"),
            ("|$_AV:first;bad&#2147483648;$|", "&#2147483648;"),
        ] {
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text}");
            assert_eq!(progress.consumed(), text.find(entity).unwrap(), "{text}");
            assert_eq!(
                progress.error().map(|error| error.offset),
                Some(text.find(entity).unwrap())
            );
            let values = match &progress.records()[0] {
                CxRecord::AtomLabels(values) | CxRecord::AtomValues(values) => values,
                other => panic!("unexpected CX slot record: {other:?}"),
            };
            assert_eq!(values, &vec![Some("first".to_owned())]);
            assert_eq!(progress.checkpoints().len(), 1);
            assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Item);
            assert_eq!(progress.checkpoints()[0].item_index, Some(0));
        }
    }

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
            "|LN:1:1.3.2.6,SgD:2,1:FIELD:info::::,Sg:n:0,1:n:ht:2:3,m:2:3.5.4,SgH:1:0|",
        )
        .unwrap();
        let typed_records = parsed
            .records()
            .iter()
            .filter(|record| !matches!(record, CxRecord::Unknown(_)))
            .collect::<Vec<_>>();
        assert_eq!(typed_records.len(), 5);
        assert!(matches!(typed_records[0], CxRecord::LinkNodes(_)));
        assert!(matches!(typed_records[1], CxRecord::DataSGroup(_)));
        assert!(matches!(typed_records[2], CxRecord::PolymerSGroup(_)));
        assert!(matches!(typed_records[3], CxRecord::VariableAttachments(_)));
        assert!(matches!(typed_records[4], CxRecord::SGroupHierarchy(_)));
        match typed_records[1] {
            CxRecord::DataSGroup(group) => {
                assert_eq!(group.atoms, vec![2, 1]);
                assert_eq!(group.field_name, "FIELD");
                assert_eq!(group.data, "info");
            }
            other => panic!("unexpected record: {other:?}"),
        }
    }

    #[test]
    fn cx_progress_properties_commits_items_at_source_cursors() {
        let text = "|atomProp:0.first.one:1.second.two|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());

        let CxRecord::AtomProperties(properties) = &progress.records()[0] else {
            panic!("expected atomProp record");
        };
        assert_eq!(
            properties,
            &vec![
                CxAtomProperty {
                    atom: 0,
                    name: "first".to_owned(),
                    value: "one".to_owned(),
                },
                CxAtomProperty {
                    atom: 1,
                    name: "second".to_owned(),
                    value: "two".to_owned(),
                },
            ]
        );
        assert_eq!(progress.checkpoints().len(), 3);
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Item);
        assert_eq!(progress.checkpoints()[0].item_index, Some(0));
        assert_eq!(
            progress.checkpoints()[0].cursor,
            text.find(":1").expect("first property separator")
        );
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Item);
        assert_eq!(progress.checkpoints()[1].item_index, Some(1));
        let close = text.rfind('|').expect("outer closing pipe");
        assert_eq!(progress.checkpoints()[1].cursor, close);
        assert_eq!(progress.checkpoints()[2].phase, CxProgressPhase::Complete);
        assert_eq!(progress.checkpoints()[2].item_index, None);
        assert_eq!(progress.checkpoints()[2].cursor, close);
    }

    #[test]
    fn cx_progress_properties_preserves_source_skip_and_empty_field_rules() {
        let text = "|atomProp:0..ignored:1.empty.:x2.good.kept|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());

        let CxRecord::AtomProperties(properties) = &progress.records()[0] else {
            panic!("expected atomProp record");
        };
        assert_eq!(
            properties,
            &[CxAtomProperty {
                atom: 2,
                name: "good".to_owned(),
                value: "kept".to_owned(),
            }]
        );
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Item);
        assert_eq!(progress.checkpoints()[0].item_index, Some(0));
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Complete);
    }

    #[test]
    fn cx_progress_properties_reports_source_failures_and_retains_prior_items() {
        let missing_first_dot = "|atomProp:0x.foo|";
        let progress = parse_cx_extensions_progress(missing_first_dot);
        let failure = missing_first_dot.find('x').expect("bad first separator");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));

        let missing_second_dot = "|atomProp:0.name|";
        let progress = parse_cx_extensions_progress(missing_second_dot);
        let failure = missing_second_dot.len();
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));

        let malformed_value = "|atomProp:0.first.kept:1.second.&#oops|";
        let progress = parse_cx_extensions_progress(malformed_value);
        let failure = malformed_value.find("&#oops").expect("bad entity");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));
        let CxRecord::AtomProperties(properties) = &progress.records()[0] else {
            panic!("expected partial atomProp record");
        };
        assert_eq!(properties.len(), 1);
        assert_eq!(properties[0].name, "first");
        assert_eq!(properties[0].value, "kept");
        assert_eq!(progress.checkpoints().len(), 1);
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Item);
        assert_eq!(
            progress.checkpoints()[0].cursor,
            malformed_value.find(":1").unwrap()
        );

        let overflow = "|atomProp:0.first.kept:4294967296.second.bad|";
        let progress = parse_cx_extensions_progress(overflow);
        let number_start = overflow.find("4294967296").expect("overflowing source int");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), number_start + 10);
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(number_start)
        );
        let CxRecord::AtomProperties(properties) = &progress.records()[0] else {
            panic!("expected partial atomProp record");
        };
        assert_eq!(properties.len(), 1);
        assert_eq!(properties[0].name, "first");
        assert_eq!(progress.checkpoints().len(), 1);

        let missing_close = "|atomProp:0.first.kept";
        let progress = parse_cx_extensions_progress(missing_close);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), missing_close.len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(missing_close.len())
        );
        let CxRecord::AtomProperties(properties) = &progress.records()[0] else {
            panic!("expected partial atomProp record");
        };
        assert_eq!(properties.len(), 1);
        assert_eq!(properties[0].value, "kept");
        assert_eq!(progress.checkpoints().len(), 1);
        assert_eq!(progress.checkpoints()[0].cursor, missing_close.len());
    }

    #[test]
    fn cx_progress_bonds_emits_source_order_items_and_completions() {
        let coordinate = "|C:1.0,0.1|";
        let progress = parse_cx_extensions_progress(coordinate);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), coordinate.len());
        let CxRecord::CoordinateBonds(annotation) = &progress.records()[0] else {
            panic!("expected coordinate-bond record");
        };
        assert_eq!(annotation.kind, CxCoordinateBondKind::Dative);
        assert_eq!(
            annotation.bonds,
            vec![
                CxBondReference { atom: 1, bond: 0 },
                CxBondReference { atom: 0, bond: 1 },
            ]
        );
        assert_eq!(progress.checkpoints().len(), 3);
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Item);
        assert_eq!(progress.checkpoints()[0].item_index, Some(0));
        assert_eq!(
            progress.checkpoints()[0].cursor,
            coordinate.find(',').expect("pair separator")
        );
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Item);
        assert_eq!(progress.checkpoints()[1].item_index, Some(1));
        let close = coordinate.rfind('|').expect("closing pipe");
        assert_eq!(progress.checkpoints()[1].cursor, close);
        assert_eq!(progress.checkpoints()[2].phase, CxProgressPhase::Complete);
        assert_eq!(progress.checkpoints()[2].cursor, close);

        let hydrogen = parse_cx_extensions_progress("|H:0.0|");
        assert!(hydrogen.is_complete());
        assert!(matches!(
            &hydrogen.records()[0],
            CxRecord::CoordinateBonds(annotation)
                if annotation.kind == CxCoordinateBondKind::Hydrogen
                    && annotation.bonds == [CxBondReference { atom: 0, bond: 0 }]
        ));

        let zero = parse_cx_extensions_progress("|Z:0,2|");
        assert!(zero.is_complete());
        assert!(matches!(&zero.records()[0], CxRecord::ZeroBonds(indices) if indices == &[0, 2]));
        assert_eq!(zero.checkpoints().len(), 3);
        assert_eq!(zero.checkpoints()[0].phase, CxProgressPhase::Item);
        assert_eq!(
            zero.checkpoints()[0].cursor,
            "|Z:0,2|".find(',').expect("zero-bond separator")
        );
        assert_eq!(zero.checkpoints()[1].phase, CxProgressPhase::Item);
        assert_eq!(zero.checkpoints()[1].cursor, "|Z:0,2|".rfind('|').unwrap());
        assert_eq!(zero.checkpoints()[2].phase, CxProgressPhase::Complete);
        assert_eq!(zero.checkpoints()[2].cursor, "|Z:0,2|".rfind('|').unwrap());
    }

    #[test]
    fn cx_progress_bonds_retains_items_and_cursors_on_pair_and_integer_failures() {
        let missing_coordinate_colon = "|C-1.0|";
        let progress = parse_cx_extensions_progress(missing_coordinate_colon);
        let failure = missing_coordinate_colon
            .find('-')
            .expect("missing coordinate colon");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));

        let missing_zero_colon = "|Zx|";
        let progress = parse_cx_extensions_progress(missing_zero_colon);
        let failure = missing_zero_colon
            .find('x')
            .expect("missing zero-bond colon");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));

        let missing_first_separator = "|C:1x.0|";
        let progress = parse_cx_extensions_progress(missing_first_separator);
        let failure = missing_first_separator
            .find('x')
            .expect("invalid pair separator");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));
        assert_eq!(progress.checkpoints().len(), 0);

        let missing_second_integer = "|H:1.|";
        let progress = parse_cx_extensions_progress(missing_second_integer);
        let failure = missing_second_integer.rfind('|').expect("closing pipe");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));
        assert_eq!(progress.checkpoints().len(), 0);

        let partial_pair = "|C:1.0,0x.1|";
        let progress = parse_cx_extensions_progress(partial_pair);
        let failure = partial_pair.find('x').expect("invalid later pair");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));
        let CxRecord::CoordinateBonds(annotation) = &progress.records()[0] else {
            panic!("expected partial coordinate-bond record");
        };
        assert_eq!(annotation.bonds, [CxBondReference { atom: 1, bond: 0 }]);
        assert_eq!(progress.checkpoints().len(), 1);
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Item);
        assert_eq!(
            progress.checkpoints()[0].cursor,
            partial_pair.find(',').expect("first pair separator")
        );

        let coordinate_overflow = "|C:1.0,4294967296.1|";
        let progress = parse_cx_extensions_progress(coordinate_overflow);
        let number_start = coordinate_overflow
            .find("4294967296")
            .expect("overflowing coordinate atom index");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), number_start + 10);
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(number_start)
        );
        let CxRecord::CoordinateBonds(annotation) = &progress.records()[0] else {
            panic!("expected partial coordinate-bond record");
        };
        assert_eq!(annotation.bonds, [CxBondReference { atom: 1, bond: 0 }]);
        assert_eq!(progress.checkpoints().len(), 1);

        let zero_overflow = "|Z:0,4294967296|";
        let progress = parse_cx_extensions_progress(zero_overflow);
        let number_start = zero_overflow
            .find("4294967296")
            .expect("overflowing zero-bond index");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), number_start + 10);
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(number_start)
        );
        assert!(matches!(&progress.records()[0], CxRecord::ZeroBonds(indices) if indices == &[0]));
        assert_eq!(progress.checkpoints().len(), 1);
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Item);
        assert_eq!(
            progress.checkpoints()[0].cursor,
            zero_overflow.find(',').unwrap()
        );
    }

    #[test]
    fn cx_progress_radicals_dispatches_chained_sections_at_source_cursors() {
        let text = "|^1:0,1^2:2^3:3^4:4^5:5^6:6^7:7|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        let CxRecord::Radicals(radicals) = &progress.records()[0] else {
            panic!("radical record expected");
        };
        assert_eq!(
            radicals
                .iter()
                .map(|radical| (radical.atom, radical.electrons))
                .collect::<Vec<_>>(),
            [
                (0, 1),
                (1, 1),
                (2, 2),
                (3, 2),
                (4, 2),
                (5, 3),
                (6, 3),
                (7, 3)
            ]
        );

        let close = text.rfind('|').expect("closing pipe");
        let expected_cursors = [
            text.find(',').expect("first atom separator"),
            text.find("^2").expect("second radical section"),
            text.find("^3").expect("third radical section"),
            text.find("^4").expect("fourth radical section"),
            text.find("^5").expect("fifth radical section"),
            text.find("^6").expect("sixth radical section"),
            text.find("^7").expect("seventh radical section"),
            close,
        ];
        assert_eq!(progress.checkpoints().len(), expected_cursors.len() + 1);
        for (index, (checkpoint, expected_cursor)) in progress.checkpoints()[..8]
            .iter()
            .zip(expected_cursors)
            .enumerate()
        {
            assert_eq!(checkpoint.record_index, 0);
            assert_eq!(checkpoint.item_index, Some(index));
            assert_eq!(checkpoint.phase, CxProgressPhase::Item);
            assert_eq!(checkpoint.cursor, expected_cursor);
        }
        let complete = progress
            .checkpoints()
            .last()
            .expect("completion checkpoint");
        assert_eq!(complete.phase, CxProgressPhase::Complete);
        assert_eq!(complete.item_index, None);
        assert_eq!(complete.cursor, close);

        let trailing_nondigit = "|^1:0,|";
        let progress = parse_cx_extensions_progress(trailing_nondigit);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), trailing_nondigit.len());
        assert_eq!(
            progress.checkpoints()[0].cursor,
            trailing_nondigit.find(',').expect("trailing comma")
        );
        assert_eq!(
            progress.checkpoints()[1].cursor,
            trailing_nondigit.rfind('|').expect("closing pipe")
        );
        assert!(matches!(
            &progress.records()[0],
            CxRecord::Radicals(radicals)
                if radicals.iter().map(|radical| (radical.atom, radical.electrons)).collect::<Vec<_>>() == [(0, 1)]
        ));
    }

    #[test]
    fn cx_progress_radicals_preserves_failure_cursors_and_prior_items() {
        for (text, expected_cursor) in [
            ("|^", 2),
            ("|^0:0|", 2),
            ("|^8:0|", 2),
            ("|^1x0|", 3),
            ("|^1:|", 4),
        ] {
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert_eq!(progress.consumed(), expected_cursor, "{text:?}");
            assert_eq!(
                progress.error().map(|error| error.offset),
                Some(expected_cursor),
                "{text:?}"
            );
            assert!(progress.checkpoints().is_empty(), "{text:?}");
            assert!(matches!(
                progress.records().first(),
                Some(CxRecord::Radicals(radicals)) if radicals.is_empty()
            ));
        }

        let first_overflow = "|^1:4294967296|";
        let progress = parse_cx_extensions_progress(first_overflow);
        let number_start = first_overflow
            .find("4294967296")
            .expect("overflowing atom index");
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), number_start + "4294967296".len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(number_start)
        );
        assert!(progress.checkpoints().is_empty());
        assert!(matches!(
            progress.records().first(),
            Some(CxRecord::Radicals(radicals)) if radicals.is_empty()
        ));

        for text in ["|^1:0", "|^1:0,", "|^1:0,4294967296|"] {
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert!(progress.error().is_some(), "{text:?}");
            assert!(matches!(
                progress.records().first(),
                Some(CxRecord::Radicals(radicals))
                    if radicals.first().is_some_and(|radical| radical.atom == 0 && radical.electrons == 1)
            ));
            assert_eq!(progress.checkpoints().len(), 1, "{text:?}");
            assert_eq!(progress.checkpoints()[0].item_index, Some(0));
            assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Item);
        }

        let missing_later_index = "|^1:0,";
        let progress = parse_cx_extensions_progress(missing_later_index);
        assert_eq!(progress.consumed(), missing_later_index.len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(missing_later_index.len())
        );
        assert_eq!(
            progress.checkpoints()[0].cursor,
            missing_later_index
                .find(',')
                .expect("later-index separator")
        );

        let later_overflow = "|^1:0,4294967296|";
        let progress = parse_cx_extensions_progress(later_overflow);
        let number_start = later_overflow.find("4294967296").expect("later overflow");
        assert_eq!(progress.consumed(), number_start + "4294967296".len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(number_start)
        );
        assert_eq!(
            progress.checkpoints()[0].cursor,
            later_overflow.find(',').expect("first-item separator")
        );
    }

    #[test]
    fn cx_progress_stereo_records_ordered_groups_at_source_cursors() {
        let text = "|a:0,1,o2:2,3|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        assert_eq!(
            progress.records(),
            &[
                CxRecord::EnhancedStereo(super::CxEnhancedStereo {
                    kind: CxStereoGroupKind::Absolute,
                    group_id: 0,
                    atoms: vec![0, 1],
                }),
                CxRecord::EnhancedStereo(super::CxEnhancedStereo {
                    kind: CxStereoGroupKind::Or,
                    group_id: 2,
                    atoms: vec![2, 3],
                }),
            ]
        );

        let first_atom = text.find("0,1").expect("first group members");
        let second_group = text.find("o2").expect("second group");
        let second_atom = text.find("2,3").expect("second group members");
        let expected = [
            (0, None, CxProgressPhase::Begin, text.find(':').unwrap() + 1),
            (0, Some(0), CxProgressPhase::Item, first_atom + 1),
            (
                0,
                Some(1),
                CxProgressPhase::Item,
                text.find("1,o2").unwrap() + 1,
            ),
            (0, None, CxProgressPhase::Complete, second_group),
            (
                1,
                None,
                CxProgressPhase::Begin,
                text.find("2:").unwrap() + 2,
            ),
            (1, Some(0), CxProgressPhase::Item, second_atom + 1),
            (
                1,
                Some(1),
                CxProgressPhase::Item,
                text.find("3|").unwrap() + 1,
            ),
            (1, None, CxProgressPhase::Complete, text.rfind('|').unwrap()),
        ];
        assert_eq!(progress.checkpoints().len(), expected.len());
        for (checkpoint, (record_index, item_index, phase, cursor)) in
            progress.checkpoints().iter().zip(expected)
        {
            assert_eq!(checkpoint.record_index, record_index);
            assert_eq!(checkpoint.item_index, item_index);
            assert_eq!(checkpoint.phase, phase);
            assert_eq!(checkpoint.cursor, cursor);
        }
    }

    #[test]
    fn cx_progress_stereo_failures_keep_exact_cursor_and_partial_group() {
        for (text, expected_cursor, expected_error_offset) in [
            ("|a!|", 2, 2),
            ("|o12x:0|", 4, 4),
            ("|o4294967296:0|", 12, 2),
        ] {
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert_eq!(progress.consumed(), expected_cursor, "{text:?}");
            assert_eq!(
                progress.error().map(|error| error.offset),
                Some(expected_error_offset),
                "{text:?}"
            );
            assert!(progress.records().is_empty(), "{text:?}");
            assert!(progress.checkpoints().is_empty(), "{text:?}");
        }

        let text = "|a:0,4294967296|";
        let number_start = text.find("4294967296").expect("overflowing atom index");
        let progress = parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), number_start + "4294967296".len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(number_start)
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::EnhancedStereo(stereo)]
                if stereo.kind == CxStereoGroupKind::Absolute
                    && stereo.group_id == 0
                    && stereo.atoms == [0]
        ));
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Begin);
        assert_eq!(progress.checkpoints()[1].item_index, Some(0));
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Item);
        assert!(
            progress
                .checkpoints()
                .iter()
                .all(|checkpoint| checkpoint.phase != CxProgressPhase::Complete)
        );
    }

    #[test]
    fn cx_progress_linknodes_emits_source_order_items_and_cursors() {
        let text = "|LN:0:1.3.1x2,1:2.4.2.0|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        assert_eq!(
            progress.records(),
            &[CxRecord::LinkNodes(vec![
                CxLinkNode {
                    atom: 0,
                    start_repetitions: 1,
                    end_repetitions: 3,
                    outer_atoms: Some([1, 2]),
                },
                CxLinkNode {
                    atom: 1,
                    start_repetitions: 2,
                    end_repetitions: 4,
                    outer_atoms: Some([2, 0]),
                },
            ])]
        );

        let checkpoints = progress.checkpoints();
        assert_eq!(checkpoints.len(), 4);
        assert_eq!(checkpoints[0].phase, CxProgressPhase::Begin);
        assert_eq!(checkpoints[0].cursor, text.find("LN:").expect("prefix") + 3);
        assert_eq!(checkpoints[1].item_index, Some(0));
        assert_eq!(checkpoints[1].phase, CxProgressPhase::Item);
        assert_eq!(
            checkpoints[1].cursor,
            text.find(',').expect("first link-node separator") + 1
        );
        assert_eq!(checkpoints[2].item_index, Some(1));
        assert_eq!(checkpoints[2].phase, CxProgressPhase::Item);
        assert_eq!(
            checkpoints[2].cursor,
            text.rfind('|').expect("closing pipe")
        );
        assert_eq!(checkpoints[3].phase, CxProgressPhase::Complete);
        assert_eq!(checkpoints[3].cursor, checkpoints[2].cursor);
    }

    #[test]
    fn cx_progress_linknodes_failures_keep_exact_cursor_and_partial_items() {
        for (text, expected_cursor, expected_error_offset) in [
            ("|LN?|", 1, 1),
            ("|LN:0|", 5, 5),
            ("|LN:0x|", 5, 5),
            ("|LN:0:.2.3|", 6, 6),
            ("|LN:0:1x2.3|", 7, 7),
            ("|LN:0:1.|", 8, 8),
            ("|LN:0:1.2.|", 10, 10),
            ("|LN:0:1.2.3.|", 12, 12),
            ("|LN:4294967296:1.2|", 14, 4),
        ] {
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert_eq!(progress.consumed(), expected_cursor, "{text:?}");
            assert_eq!(
                progress.error().map(|error| error.offset),
                Some(expected_error_offset),
                "{text:?}"
            );
            if text.starts_with("|LN:") {
                assert!(matches!(
                    progress.records(),
                    [CxRecord::LinkNodes(nodes)] if nodes.is_empty()
                ));
                assert_eq!(progress.checkpoints().len(), 1, "{text:?}");
                assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Begin);
                assert_eq!(progress.checkpoints()[0].cursor, 4, "{text:?}");
            } else {
                assert!(progress.records().is_empty());
                assert!(progress.checkpoints().is_empty());
            }
        }

        let text = "|LN:0:1.2.3.4,4294967296:1.2.3.4|";
        let number_start = text.find("4294967296").expect("later center overflow");
        let progress = parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), number_start + "4294967296".len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(number_start)
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::LinkNodes(nodes)]
                if nodes == &[CxLinkNode {
                    atom: 0,
                    start_repetitions: 1,
                    end_repetitions: 2,
                    outer_atoms: Some([3, 4]),
                }]
        ));
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Begin);
        assert_eq!(progress.checkpoints()[1].item_index, Some(0));
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Item);
    }

    #[test]
    fn cx_progress_constraints_emits_source_order_items_and_cursors() {
        let text = "|u:0,1,rb:0:3,1:4,2:*,s:0:2,1:*|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        assert_eq!(
            progress.records(),
            &[
                CxRecord::Unsaturation(vec![0, 1]),
                CxRecord::RingBonds(vec![
                    CxRingBond {
                        atom: 0,
                        constraint: CxCountConstraint::Exact(3),
                    },
                    CxRingBond {
                        atom: 1,
                        constraint: CxCountConstraint::LessEqual(4),
                    },
                    CxRingBond {
                        atom: 2,
                        constraint: CxCountConstraint::QueryScan,
                    },
                ]),
                CxRecord::Substitution(vec![
                    CxAtomConstraint {
                        atom: 0,
                        constraint: CxCountConstraint::Exact(2),
                    },
                    CxAtomConstraint {
                        atom: 1,
                        constraint: CxCountConstraint::QueryScan,
                    },
                ]),
            ]
        );

        let checkpoints = progress.checkpoints();
        assert_eq!(checkpoints.len(), 13);
        assert_eq!(
            checkpoints
                .iter()
                .map(|checkpoint| (
                    checkpoint.record_index,
                    checkpoint.item_index,
                    checkpoint.phase
                ))
                .collect::<Vec<_>>(),
            [
                (0, None, CxProgressPhase::Begin),
                (0, Some(0), CxProgressPhase::Item),
                (0, Some(1), CxProgressPhase::Item),
                (0, None, CxProgressPhase::Complete),
                (1, None, CxProgressPhase::Begin),
                (1, Some(0), CxProgressPhase::Item),
                (1, Some(1), CxProgressPhase::Item),
                (1, Some(2), CxProgressPhase::Item),
                (1, None, CxProgressPhase::Complete),
                (2, None, CxProgressPhase::Begin),
                (2, Some(0), CxProgressPhase::Item),
                (2, Some(1), CxProgressPhase::Item),
                (2, None, CxProgressPhase::Complete),
            ]
        );
        assert_eq!(checkpoints[0].cursor, text.find("0,1").unwrap());
        assert_eq!(checkpoints[1].cursor, text.find("0,1").unwrap() + 1);
        assert_eq!(checkpoints[2].cursor, text.find("1,rb").unwrap() + 1);
        assert_eq!(checkpoints[3].cursor, text.find("rb:").unwrap());
        assert_eq!(checkpoints[4].cursor, text.find("rb:").unwrap() + 3);
        assert_eq!(checkpoints[5].cursor, text.find("0:3").unwrap() + 3);
        assert_eq!(checkpoints[6].cursor, text.find("1:4").unwrap() + 3);
        assert_eq!(checkpoints[7].cursor, text.find("2:*").unwrap() + 3);
        assert_eq!(checkpoints[8].cursor, text.find("s:0:2").unwrap());
        assert_eq!(checkpoints[9].cursor, text.find("s:").unwrap() + 2);
        assert_eq!(checkpoints[10].cursor, text.find("0:2").unwrap() + 3);
        assert_eq!(checkpoints[11].cursor, text.find("1:*").unwrap() + 3);
        assert_eq!(checkpoints[12].cursor, text.rfind('|').unwrap());
    }

    #[test]
    fn cx_progress_constraints_integer_overflow_keeps_consumed_cursor() {
        for text in ["|rb:4294967296:3|", "|rb:0:4294967296|"] {
            let number_start = text.find("4294967296").expect("overflowing ring integer");
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert_eq!(
                progress.consumed(),
                number_start + "4294967296".len(),
                "{text:?}"
            );
            assert_eq!(
                progress.error().map(|error| error.offset),
                Some(number_start),
                "{text:?}"
            );
            assert!(matches!(
                progress.records(),
                [CxRecord::RingBonds(constraints)] if constraints.is_empty()
            ));
            assert_eq!(progress.checkpoints().len(), 1, "{text:?}");
            assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Begin);
        }

        for text in ["|s:4294967296:3|", "|s:0:4294967296|"] {
            let number_start = text
                .find("4294967296")
                .expect("overflowing substitution integer");
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert_eq!(
                progress.consumed(),
                number_start + "4294967296".len(),
                "{text:?}"
            );
            assert_eq!(
                progress.error().map(|error| error.offset),
                Some(number_start),
                "{text:?}"
            );
            assert!(matches!(
                progress.records(),
                [CxRecord::Substitution(constraints)] if constraints.is_empty()
            ));
            assert_eq!(progress.checkpoints().len(), 1, "{text:?}");
            assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Begin);
        }
    }

    #[test]
    fn cx_progress_constraints_fail_at_source_cursors_and_keep_prior_items() {
        let missing_unsaturation_colon = "|ux|";
        let progress = parse_cx_extensions_progress(missing_unsaturation_colon);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), 2);
        assert_eq!(progress.error().map(|error| error.offset), Some(2));
        assert!(progress.records().is_empty());
        assert!(progress.checkpoints().is_empty());

        let unsaturation_overflow = "|u:0,4294967296|";
        let number_start = unsaturation_overflow
            .find("4294967296")
            .expect("overflowing unsaturation index");
        let progress = parse_cx_extensions_progress(unsaturation_overflow);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), number_start + "4294967296".len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(number_start)
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::Unsaturation(indices)] if indices.as_slice() == [0]
        ));
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Item);

        for (text, expected_error, prior_constraints) in [
            ("|rb:0x3|", 5, 0),
            ("|rb:0:", 5, 0),
            ("|rb:0:|", 6, 0),
            ("|rb:0:1|", 7, 0),
            ("|rb:0:3,1:1|", 11, 1),
        ] {
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert_eq!(progress.consumed(), expected_error, "{text:?}");
            assert_eq!(
                progress.error().map(|error| error.offset),
                Some(expected_error),
                "{text:?}"
            );
            assert!(matches!(
                progress.records().first(),
                Some(CxRecord::RingBonds(constraints))
                    if constraints.len() == prior_constraints
            ));
            assert_eq!(
                progress.checkpoints().len(),
                1 + prior_constraints,
                "{text:?}"
            );
            assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Begin);
            assert!(
                progress
                    .checkpoints()
                    .iter()
                    .all(|checkpoint| checkpoint.phase != CxProgressPhase::Complete)
            );
        }

        for (text, expected_error, prior_constraints) in [
            ("|s:0:x|", 5, 0),
            ("|s:0:", 4, 0),
            ("|s:0:|", 5, 0),
            ("|s:0:2,1:x|", 9, 1),
        ] {
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert_eq!(progress.consumed(), expected_error, "{text:?}");
            assert_eq!(
                progress.error().map(|error| error.offset),
                Some(expected_error),
                "{text:?}"
            );
            assert!(matches!(
                progress.records().first(),
                Some(CxRecord::Substitution(constraints))
                    if constraints.len() == prior_constraints
            ));
            assert_eq!(
                progress.checkpoints().len(),
                1 + prior_constraints,
                "{text:?}"
            );
            assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Begin);
            assert!(
                progress
                    .checkpoints()
                    .iter()
                    .all(|checkpoint| checkpoint.phase != CxProgressPhase::Complete)
            );
        }
    }

    #[test]
    fn cx_progress_sgroups_tracks_field_sequence_and_source_cursors() {
        let text = "|SgD:2,0:NAME:value,with,comma&#58;embedded:op:unit:tag:(1,2)|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        let parsed = parse_cx_extensions(text).expect("complete data SGroup");
        assert_eq!(parsed.records(), progress.records());

        let [CxRecord::DataSGroup(group)] = progress.records() else {
            panic!("expected one data SGroup");
        };
        assert_eq!(group.atoms, vec![2, 0]);
        assert_eq!(group.field_name, "NAME");
        assert_eq!(group.data, "value,with,comma:embedded");
        assert_eq!(group.query_op, "op");
        assert_eq!(group.field_info, "unit");
        assert_eq!(group.field_tag, "tag");
        assert_eq!(group.coordinates.as_deref(), Some("(1,2"));

        let expected = [
            (None, CxProgressPhase::Begin, text.find(":NAME").unwrap()),
            (
                Some(0),
                CxProgressPhase::Item,
                text.find(":value").unwrap() + 1,
            ),
            (
                Some(1),
                CxProgressPhase::Item,
                text.find(":op").unwrap() + 1,
            ),
            (
                Some(2),
                CxProgressPhase::Item,
                text.find(":unit").unwrap() + 1,
            ),
            (
                Some(3),
                CxProgressPhase::Item,
                text.find(":tag").unwrap() + 1,
            ),
            (Some(4), CxProgressPhase::Item, text.find("(1,2)").unwrap()),
            (Some(5), CxProgressPhase::Item, text.find(")|").unwrap() + 1),
            (
                None,
                CxProgressPhase::Complete,
                text.find(")|").unwrap() + 1,
            ),
        ];
        assert_eq!(progress.checkpoints().len(), expected.len());
        for (record_index, (checkpoint, (item_index, phase, cursor))) in
            progress.checkpoints().iter().zip(expected).enumerate()
        {
            assert_eq!(checkpoint.record_index, 0, "checkpoint {record_index}");
            assert_eq!(
                checkpoint.item_index, item_index,
                "checkpoint {record_index}"
            );
            assert_eq!(checkpoint.phase, phase, "checkpoint {record_index}");
            assert_eq!(checkpoint.cursor, cursor, "checkpoint {record_index}");
        }
    }

    #[test]
    fn cx_progress_sgroups_failures_keep_exact_cursor_and_partial_fields() {
        let invalid_prefix = "|SgD?0:F:D:Q:I:T:|";
        let progress = parse_cx_extensions_progress(invalid_prefix);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), 1);
        assert_eq!(progress.error().map(|error| error.offset), Some(1));
        assert!(progress.records().is_empty());
        assert!(progress.checkpoints().is_empty());

        let overflow = "|SgD:4294967296:F:D:Q:I:T:|";
        let overflow_start = overflow.find("4294967296").unwrap();
        let progress = parse_cx_extensions_progress(overflow);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), overflow_start + "4294967296".len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(overflow_start)
        );
        assert!(progress.records().is_empty());
        assert!(progress.checkpoints().is_empty());

        for (text, committed_fields) in [
            ("|SgD:0:&#oops:D:Q:I:T:|", 0),
            ("|SgD:0:F:&#oops:Q:I:T:|", 1),
            ("|SgD:0:F:D:&#oops:I:T:|", 2),
            ("|SgD:0:F:D:Q:&#oops:T:|", 3),
            ("|SgD:0:F:D:Q:I:&#oops:|", 4),
            ("|SgD:0:F:D:Q:I:T:(&#oops)|", 5),
        ] {
            let failure = text.find("&#oops").unwrap();
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert_eq!(progress.consumed(), failure, "{text:?}");
            assert_eq!(progress.error().map(|error| error.offset), Some(failure));
            let [CxRecord::DataSGroup(group)] = progress.records() else {
                panic!("expected partial data SGroup for {text:?}");
            };
            assert_eq!(
                group.field_name.is_empty(),
                committed_fields == 0,
                "{text:?}"
            );
            assert_eq!(group.data.is_empty(), committed_fields < 2, "{text:?}");
            assert_eq!(group.query_op.is_empty(), committed_fields < 3, "{text:?}");
            assert_eq!(
                group.field_info.is_empty(),
                committed_fields < 4,
                "{text:?}"
            );
            assert_eq!(group.field_tag.is_empty(), committed_fields < 5, "{text:?}");
            assert!(group.coordinates.is_none(), "{text:?}");
            assert_eq!(progress.checkpoints().len(), 1 + committed_fields);
            assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Begin);
            assert!(
                progress
                    .checkpoints()
                    .iter()
                    .all(|checkpoint| checkpoint.phase != CxProgressPhase::Complete)
            );
        }

        let after_complete = "|SgD:0:F:one::::SgD:0:G:bad&#oops:Q:I:T:|";
        let failure = after_complete.find("&#oops").unwrap();
        let progress = parse_cx_extensions_progress(after_complete);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.records().len(), 2);
        assert!(progress.checkpoints().iter().any(|checkpoint| {
            checkpoint.record_index == 0 && checkpoint.phase == CxProgressPhase::Complete
        }));
        assert!(progress.checkpoints().iter().all(|checkpoint| {
            checkpoint.record_index != 1 || checkpoint.phase != CxProgressPhase::Complete
        }));
    }

    #[test]
    fn cx_progress_sgroups_complete_helper_before_outer_pipe_failure() {
        let text = "|SgD:0:F:D:Q:I:T:(raw|";
        let progress = parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        assert_eq!(progress.error().map(|error| error.offset), Some(text.len()));
        let [CxRecord::DataSGroup(group)] = progress.records() else {
            panic!("expected completed data SGroup syntax");
        };
        assert_eq!(group.coordinates.as_deref(), Some("(raw|"));
        assert!(progress.checkpoints().iter().any(|checkpoint| {
            checkpoint.phase == CxProgressPhase::Complete && checkpoint.cursor == text.len()
        }));
    }

    #[test]
    fn cx_progress_polymer_tracks_optional_fields_and_dispatch_types() {
        let text = "|Sg:n:2,0,2:repeat:hh&#44;f:1,0,:2,|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        let parsed = parse_cx_extensions(text).expect("complete polymer SGroup");
        assert_eq!(parsed.records(), progress.records());
        let [CxRecord::PolymerSGroup(polymer)] = progress.records() else {
            panic!("expected one polymer SGroup");
        };
        assert_eq!(polymer.type_code, "n");
        assert_eq!(polymer.atoms, [2, 0, 2]);
        assert_eq!(polymer.label, "repeat");
        assert_eq!(polymer.connect, "hh,f");
        assert_eq!(polymer.head_crossings, [1, 0]);
        assert_eq!(polymer.tail_crossings, [2]);

        let atom_start = text.find("2,0,2").expect("atom list");
        let head_start = text.find("1,0,").expect("head list");
        let tail_start = text.rfind(":2,").expect("tail list") + 1;
        let expected = [
            (None, CxProgressPhase::Begin, atom_start),
            (Some(0), CxProgressPhase::Item, atom_start + 1),
            (Some(1), CxProgressPhase::Item, atom_start + 3),
            (Some(2), CxProgressPhase::Item, atom_start + 5),
            (
                Some(3),
                CxProgressPhase::Item,
                text.find(":hh").expect("label delimiter"),
            ),
            (
                Some(4),
                CxProgressPhase::Item,
                text.find(":1,0,").expect("connect delimiter"),
            ),
            (Some(5), CxProgressPhase::Item, head_start + 1),
            (Some(6), CxProgressPhase::Item, head_start + 3),
            (Some(7), CxProgressPhase::Item, tail_start + 1),
            (None, CxProgressPhase::Complete, text.len() - 1),
        ];
        assert_eq!(progress.checkpoints().len(), expected.len());
        for (checkpoint, (item_index, phase, cursor)) in progress.checkpoints().iter().zip(expected)
        {
            assert_eq!(checkpoint.record_index, 0);
            assert_eq!(checkpoint.item_index, item_index);
            assert_eq!(checkpoint.phase, phase);
            assert_eq!(checkpoint.cursor, cursor);
        }

        for type_code in [
            "n", "mon", "mer", "co", "xl", "mod", "mix", "f", "any", "gen", "c", "grf", "alt",
            "ran", "blk",
        ] {
            let input = format!("|Sg:{type_code}:0|");
            let progress = parse_cx_extensions_progress(&input);
            assert!(progress.is_complete(), "{input:?}");
            assert!(matches!(
                progress.records(),
                [CxRecord::PolymerSGroup(polymer)] if polymer.type_code == type_code
            ));
            assert_eq!(
                parse_cx_extensions(&input)
                    .expect("supported source polymer type")
                    .records(),
                progress.records(),
                "{input:?}"
            );
        }
    }

    #[test]
    fn cx_progress_polymer_failures_keep_partial_fields_and_cursors() {
        let invalid_prefix = "|Sgx:n:0|";
        let progress = parse_cx_extensions_progress(invalid_prefix);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), 1);
        assert_eq!(progress.error().map(|error| error.offset), Some(1));
        assert!(progress.records().is_empty());
        assert!(progress.checkpoints().is_empty());

        let unknown_type = "|Sg:unknown:0|";
        let progress = parse_cx_extensions_progress(unknown_type);
        let unknown_end = unknown_type.find(":0").expect("type delimiter") + 1;
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), unknown_end);
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(unknown_end)
        );
        assert!(progress.records().is_empty());
        assert!(progress.checkpoints().is_empty());

        let overflow_after_item = "|Sg:n:0,4294967296|";
        let overflow_start = overflow_after_item.find("4294967296").unwrap();
        let progress = parse_cx_extensions_progress(overflow_after_item);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), overflow_start + "4294967296".len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(overflow_start)
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::PolymerSGroup(polymer)] if polymer.atoms == [0]
        ));
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[1].item_index, Some(0));
        assert!(
            progress
                .checkpoints()
                .iter()
                .all(|checkpoint| checkpoint.phase != CxProgressPhase::Complete)
        );

        for (text, bad_number, expected_atoms, expected_head, expected_tail, committed_items) in [
            ("|Sg:n:4294967296|", "4294967296", vec![], vec![], vec![], 0),
            (
                "|Sg:n:0:lab:eu:4294967296|",
                "4294967296",
                vec![0],
                vec![],
                vec![],
                3,
            ),
            (
                "|Sg:n:0:lab:eu:1:4294967296|",
                "4294967296",
                vec![0],
                vec![1],
                vec![],
                4,
            ),
        ] {
            let error_offset = text.find(bad_number).unwrap();
            let progress = parse_cx_extensions_progress(text);
            assert!(!progress.is_complete(), "{text:?}");
            assert_eq!(
                progress.consumed(),
                error_offset + bad_number.len(),
                "{text:?}"
            );
            assert_eq!(
                progress.error().map(|error| error.offset),
                Some(error_offset)
            );
            let [CxRecord::PolymerSGroup(polymer)] = progress.records() else {
                panic!("expected partial polymer SGroup for {text:?}");
            };
            assert_eq!(polymer.atoms, expected_atoms, "{text:?}");
            assert_eq!(polymer.head_crossings, expected_head, "{text:?}");
            assert_eq!(polymer.tail_crossings, expected_tail, "{text:?}");
            assert_eq!(
                progress.checkpoints().len(),
                1 + committed_items,
                "{text:?}"
            );
            assert!(
                progress
                    .checkpoints()
                    .iter()
                    .all(|checkpoint| checkpoint.phase != CxProgressPhase::Complete)
            );
        }

        let malformed_entity = "|Sg:n:0:bad&#oops:eu|";
        let entity_start = malformed_entity.find("&#oops").unwrap();
        let progress = parse_cx_extensions_progress(malformed_entity);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), entity_start);
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(entity_start)
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::PolymerSGroup(polymer)] if polymer.atoms == [0] && polymer.label.is_empty()
        ));
        assert_eq!(progress.checkpoints().len(), 2);
    }

    #[test]
    fn cx_progress_polymer_keeps_completed_helper_before_later_errors() {
        let after_complete = "|Sg:n:0,Sg:unknown:1|";
        let progress = parse_cx_extensions_progress(after_complete);
        let error_offset = after_complete.find(":1").unwrap() + 1;
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), error_offset);
        assert_eq!(progress.records().len(), 1);
        assert!(matches!(
            progress.records(),
            [CxRecord::PolymerSGroup(polymer)] if polymer.atoms == [0]
        ));
        assert!(progress.checkpoints().iter().any(|checkpoint| {
            checkpoint.record_index == 0 && checkpoint.phase == CxProgressPhase::Complete
        }));

        let missing_outer_pipe = "|Sg:n:0:repeat:hh";
        let progress = parse_cx_extensions_progress(missing_outer_pipe);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), missing_outer_pipe.len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(missing_outer_pipe.len())
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::PolymerSGroup(polymer)]
                if polymer.atoms == [0] && polymer.label == "repeat" && polymer.connect == "hh"
        ));
        assert!(progress.checkpoints().iter().any(|checkpoint| {
            checkpoint.phase == CxProgressPhase::Complete
                && checkpoint.cursor == missing_outer_pipe.len()
        }));
    }

    #[test]
    fn cx_progress_attachments_tracks_nested_lists_at_source_cursors() {
        let text = "|m:2:3.5.4,6:7.|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        let parsed = parse_cx_extensions(text).expect("complete variable attachments");
        assert_eq!(parsed.records(), progress.records());
        assert!(matches!(
            progress.records(),
            [CxRecord::VariableAttachments(attachments)]
                if attachments == &[
                    CxVariableAttachment { atom: 2, endpoints: vec![3, 5, 4] },
                    CxVariableAttachment { atom: 6, endpoints: vec![7] },
                ]
        ));

        let row0 = text.find("2:").expect("first attachment atom");
        let row1 = text.find("6:").expect("second attachment atom");
        let row0_end = text.find(",6:").expect("first row separator");
        let closing_pipe = text.rfind('|').expect("outer closing pipe");
        let expected = [
            (None, CxProgressPhase::Begin, text.find("2:").unwrap()),
            (Some(0), CxProgressPhase::Item, row0 + 1),
            (Some(1), CxProgressPhase::Item, row0_end),
            (Some(2), CxProgressPhase::Item, row1 + 1),
            (Some(3), CxProgressPhase::Item, closing_pipe),
            (None, CxProgressPhase::Complete, closing_pipe),
        ];
        assert_eq!(progress.checkpoints().len(), expected.len());
        for (checkpoint, (item_index, phase, cursor)) in progress.checkpoints().iter().zip(expected)
        {
            assert_eq!(checkpoint.record_index, 0);
            assert_eq!(checkpoint.item_index, item_index);
            assert_eq!(checkpoint.phase, phase);
            assert_eq!(checkpoint.cursor, cursor);
        }

        for input in ["|m:|", "|m:2:3,|"] {
            let progress = parse_cx_extensions_progress(input);
            assert!(progress.is_complete(), "{input:?}");
            assert_eq!(progress.consumed(), input.len(), "{input:?}");
            assert!(matches!(
                progress.records(),
                [CxRecord::VariableAttachments(attachments)]
                    if (input == "|m:|" && attachments.is_empty())
                        || (input == "|m:2:3,|"
                            && attachments == &[CxVariableAttachment {
                                atom: 2,
                                endpoints: vec![3],
                            }])
            ));
        }
    }

    #[test]
    fn cx_progress_attachments_keeps_partial_rows_and_prior_checkpoints() {
        let invalid_prefix = "|m?";
        let progress = parse_cx_extensions_progress(invalid_prefix);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), 1);
        assert_eq!(progress.error().map(|error| error.offset), Some(1));
        assert!(progress.records().is_empty());
        assert!(progress.checkpoints().is_empty());

        let primary_overflow = "|m:4294967296|";
        let overflow_start = primary_overflow.find("4294967296").unwrap();
        let progress = parse_cx_extensions_progress(primary_overflow);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), overflow_start + "4294967296".len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(overflow_start)
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::VariableAttachments(attachments)] if attachments.is_empty()
        ));
        assert_eq!(progress.checkpoints().len(), 1);
        assert_eq!(progress.checkpoints()[0].phase, CxProgressPhase::Begin);

        let missing_colon = "|m:2x|";
        let error_offset = missing_colon.find('x').unwrap();
        let progress = parse_cx_extensions_progress(missing_colon);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), error_offset);
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(error_offset)
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::VariableAttachments(attachments)]
                if attachments == &[CxVariableAttachment { atom: 2, endpoints: vec![] }]
        ));
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[1].item_index, Some(0));
        assert_eq!(progress.checkpoints()[1].cursor, error_offset);

        let endpoint_overflow = "|m:2:3.4294967296|";
        let overflow_start = endpoint_overflow.find("4294967296").unwrap();
        let progress = parse_cx_extensions_progress(endpoint_overflow);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), overflow_start + "4294967296".len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(overflow_start)
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::VariableAttachments(attachments)]
                if attachments == &[CxVariableAttachment { atom: 2, endpoints: vec![3] }]
        ));
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[1].item_index, Some(0));
        assert!(
            progress
                .checkpoints()
                .iter()
                .all(|checkpoint| checkpoint.phase != CxProgressPhase::Complete)
        );

        let later_row_failure = "|m:2:3,4x:5|";
        let error_offset = later_row_failure.find('x').unwrap();
        let progress = parse_cx_extensions_progress(later_row_failure);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), error_offset);
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(error_offset)
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::VariableAttachments(attachments)]
                if attachments == &[
                    CxVariableAttachment { atom: 2, endpoints: vec![3] },
                    CxVariableAttachment { atom: 4, endpoints: vec![] },
                ]
        ));
        assert_eq!(progress.checkpoints().len(), 4);
        assert_eq!(progress.checkpoints()[2].item_index, Some(1));
        assert_eq!(progress.checkpoints()[2].cursor, error_offset - 2);
        assert_eq!(progress.checkpoints()[3].item_index, Some(2));
        assert_eq!(progress.checkpoints()[3].cursor, error_offset);
        assert!(
            progress
                .checkpoints()
                .iter()
                .all(|checkpoint| checkpoint.phase != CxProgressPhase::Complete)
        );
    }

    #[test]
    fn cx_progress_hierarchy_tracks_order_and_source_mutation_cursors() {
        let text = "|SgH:4:1.3,7:0.2.|";
        let progress = parse_cx_extensions_progress(text);
        let final_list_cursor = text.rfind('|').expect("closing pipe");
        let second_parent_cursor = text.find(",7").expect("second parent separator");

        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        assert_eq!(
            progress.records(),
            &[CxRecord::SGroupHierarchy(vec![
                CxSGroupHierarchy {
                    parent: 4,
                    children: vec![1, 3],
                },
                CxSGroupHierarchy {
                    parent: 7,
                    children: vec![0, 2],
                },
            ])]
        );
        let checkpoints = progress.checkpoints();
        assert_eq!(checkpoints.len(), 6);
        assert_eq!(checkpoints[0].phase, CxProgressPhase::Begin);
        assert_eq!(checkpoints[0].cursor, text.find('4').expect("first parent"));
        for (checkpoint, (item_index, cursor)) in checkpoints[1..5].iter().zip([
            (0, second_parent_cursor),
            (1, second_parent_cursor),
            (2, final_list_cursor),
            (3, final_list_cursor),
        ]) {
            assert_eq!(checkpoint.phase, CxProgressPhase::Item);
            assert_eq!(checkpoint.item_index, Some(item_index));
            assert_eq!(checkpoint.cursor, cursor);
        }
        assert_eq!(checkpoints[5].phase, CxProgressPhase::Complete);
        assert_eq!(checkpoints[5].cursor, final_list_cursor);
    }

    #[test]
    fn cx_progress_hierarchy_failures_keep_incomplete_entries_and_prior_items() {
        let invalid_prefix = "|SgH?1:0|";
        let invalid_prefix_progress = parse_cx_extensions_progress(invalid_prefix);
        assert!(!invalid_prefix_progress.is_complete());
        assert_eq!(invalid_prefix_progress.consumed(), 1);
        assert_eq!(
            invalid_prefix_progress.error().map(|error| error.offset),
            Some(1)
        );
        assert!(invalid_prefix_progress.records().is_empty());
        assert!(invalid_prefix_progress.checkpoints().is_empty());

        let missing_parent = "|SgH::0|";
        let missing_parent_progress = parse_cx_extensions_progress(missing_parent);
        assert!(!missing_parent_progress.is_complete());
        assert_eq!(missing_parent_progress.consumed(), 5);
        assert_eq!(
            missing_parent_progress.error().map(|error| error.offset),
            Some(5)
        );
        assert!(matches!(
            missing_parent_progress.records(),
            [CxRecord::SGroupHierarchy(hierarchies)] if hierarchies.is_empty()
        ));
        assert_eq!(missing_parent_progress.checkpoints().len(), 1);
        assert_eq!(
            missing_parent_progress.checkpoints()[0].phase,
            CxProgressPhase::Begin
        );

        let parent_overflow = "|SgH:4294967296:0|";
        let parent_overflow_start = parent_overflow
            .find("4294967296")
            .expect("overflowing parent index");
        let parent_overflow_progress = parse_cx_extensions_progress(parent_overflow);
        assert!(!parent_overflow_progress.is_complete());
        assert_eq!(
            parent_overflow_progress.consumed(),
            parent_overflow_start + "4294967296".len()
        );
        assert_eq!(
            parent_overflow_progress.error().map(|error| error.offset),
            Some(parent_overflow_start)
        );
        assert!(matches!(
            parent_overflow_progress.records(),
            [CxRecord::SGroupHierarchy(hierarchies)] if hierarchies.is_empty()
        ));
        assert_eq!(parent_overflow_progress.checkpoints().len(), 1);

        let missing_colon = "|SgH:1x:0|";
        let missing_colon_error = missing_colon.find('x').expect("bad parent separator");
        let missing_colon_progress = parse_cx_extensions_progress(missing_colon);
        assert!(!missing_colon_progress.is_complete());
        assert_eq!(missing_colon_progress.consumed(), missing_colon_error);
        assert_eq!(
            missing_colon_progress.error().map(|error| error.offset),
            Some(missing_colon_error)
        );
        assert!(matches!(
            missing_colon_progress.records(),
            [CxRecord::SGroupHierarchy(hierarchies)]
                if hierarchies == &[CxSGroupHierarchy { parent: 1, children: Vec::new() }]
        ));
        assert_eq!(missing_colon_progress.checkpoints().len(), 1);

        let child_overflow = "|SgH:1:0.4294967296|";
        let child_overflow_start = child_overflow
            .find("4294967296")
            .expect("overflowing child index");
        let child_overflow_progress = parse_cx_extensions_progress(child_overflow);
        assert!(!child_overflow_progress.is_complete());
        assert_eq!(
            child_overflow_progress.consumed(),
            child_overflow_start + "4294967296".len()
        );
        assert_eq!(
            child_overflow_progress.error().map(|error| error.offset),
            Some(child_overflow_start)
        );
        assert!(matches!(
            child_overflow_progress.records(),
            [CxRecord::SGroupHierarchy(hierarchies)]
                if hierarchies == &[CxSGroupHierarchy { parent: 1, children: vec![0] }]
        ));
        assert_eq!(child_overflow_progress.checkpoints().len(), 1);

        let later_parent_error = "|SgH:1:0,2x:3|";
        let later_error_offset = later_parent_error.find('x').expect("later parent error");
        let later_parent_progress = parse_cx_extensions_progress(later_parent_error);
        assert!(!later_parent_progress.is_complete());
        assert_eq!(later_parent_progress.consumed(), later_error_offset);
        assert_eq!(
            later_parent_progress.error().map(|error| error.offset),
            Some(later_error_offset)
        );
        assert!(matches!(
            later_parent_progress.records(),
            [CxRecord::SGroupHierarchy(hierarchies)] if hierarchies == &[
                CxSGroupHierarchy { parent: 1, children: vec![0] },
                CxSGroupHierarchy { parent: 2, children: Vec::new() },
            ]
        ));
        assert_eq!(later_parent_progress.checkpoints().len(), 2);
        assert_eq!(
            later_parent_progress.checkpoints()[1].phase,
            CxProgressPhase::Item
        );

        let trailing_comma = "|SgH:1:0,|";
        let comma_error = trailing_comma.rfind('|').expect("closing pipe");
        let trailing_comma_progress = parse_cx_extensions_progress(trailing_comma);
        assert!(!trailing_comma_progress.is_complete());
        assert_eq!(trailing_comma_progress.consumed(), comma_error);
        assert_eq!(
            trailing_comma_progress.error().map(|error| error.offset),
            Some(comma_error)
        );
        assert_eq!(trailing_comma_progress.checkpoints().len(), 2);
        assert_eq!(trailing_comma_progress.checkpoints()[1].item_index, Some(0));

        let empty_children = parse_cx_extensions_progress("|SgH:1:|");
        assert!(empty_children.is_complete());
        assert!(matches!(
            empty_children.records(),
            [CxRecord::SGroupHierarchy(hierarchies)]
                if hierarchies == &[CxSGroupHierarchy { parent: 1, children: Vec::new() }]
        ));
        assert_eq!(empty_children.checkpoints().len(), 2);
        assert_eq!(
            empty_children.checkpoints()[1].phase,
            CxProgressPhase::Complete
        );
    }

    #[test]
    fn cx_progress_directions_dispatches_all_kinds_at_source_cursors() {
        let text = "|wU:1.0,0.1wD:2.2w:0.0ctu:1,0c:2t:3|";
        let progress = parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        let parsed = parse_cx_extensions(text).expect("complete direction records");
        assert_eq!(parsed.records(), progress.records());

        assert!(matches!(
            &progress.records()[0],
            CxRecord::WedgedBonds(wedges)
                if wedges.len() == 2
                    && wedges[0].atom == 1
                    && wedges[0].bond == 0
                    && wedges[0].direction == CxWedgeDirection::BeginWedge
                    && wedges[0].configuration == 1
                    && wedges[1].atom == 0
                    && wedges[1].bond == 1
        ));
        assert!(matches!(
            &progress.records()[1],
            CxRecord::WedgedBonds(wedges)
                if wedges == &[crate::CxWedgeBond {
                    atom: 2,
                    bond: 2,
                    direction: CxWedgeDirection::BeginDash,
                    configuration: 3,
                }]
        ));
        assert!(matches!(
            &progress.records()[2],
            CxRecord::WedgedBonds(wedges)
                if wedges == &[crate::CxWedgeBond {
                    atom: 0,
                    bond: 0,
                    direction: CxWedgeDirection::Unknown,
                    configuration: 2,
                }]
        ));
        assert!(matches!(
            &progress.records()[3],
            CxRecord::DoubleBondStereo(stereo)
                if stereo.stereo == CxDoubleBondStereoKind::Any
                    && stereo.bonds == [1, 0]
        ));
        assert!(matches!(
            &progress.records()[4],
            CxRecord::DoubleBondStereo(stereo)
                if stereo.stereo == CxDoubleBondStereoKind::Cis
                    && stereo.bonds == [2]
        ));
        assert!(matches!(
            &progress.records()[5],
            CxRecord::DoubleBondStereo(stereo)
                if stereo.stereo == CxDoubleBondStereoKind::Trans
                    && stereo.bonds == [3]
        ));

        let first_item = text.find(",0.1").expect("second first-record pair");
        let any_item = text.find(",0c:").expect("second any-stereo bond");
        let second_record = text.find("wD:").expect("second wedge record");
        let third_record = text.find("w:").expect("unknown wedge record");
        let fourth_record = text.find("ctu:").expect("any stereo record");
        let fifth_record = text.find("c:2").expect("cis stereo record");
        let sixth_record = text.find("t:3").expect("trans stereo record");
        let close = text.rfind('|').expect("closing pipe");
        let expected = [
            (0, None, CxProgressPhase::Begin, text.find("1.0").unwrap()),
            (0, Some(0), CxProgressPhase::Item, first_item),
            (0, Some(1), CxProgressPhase::Item, second_record),
            (0, None, CxProgressPhase::Complete, second_record),
            (1, None, CxProgressPhase::Begin, text.find("2.2").unwrap()),
            (1, Some(0), CxProgressPhase::Item, third_record),
            (1, None, CxProgressPhase::Complete, third_record),
            (2, None, CxProgressPhase::Begin, third_record + 2),
            (2, Some(0), CxProgressPhase::Item, fourth_record),
            (2, None, CxProgressPhase::Complete, fourth_record),
            (3, None, CxProgressPhase::Begin, fourth_record + 4),
            (3, Some(0), CxProgressPhase::Item, any_item),
            (3, Some(1), CxProgressPhase::Item, fifth_record),
            (3, None, CxProgressPhase::Complete, fifth_record),
            (4, None, CxProgressPhase::Begin, fifth_record + 2),
            (4, Some(0), CxProgressPhase::Item, sixth_record),
            (4, None, CxProgressPhase::Complete, sixth_record),
            (5, None, CxProgressPhase::Begin, sixth_record + 2),
            (5, Some(0), CxProgressPhase::Item, close),
            (5, None, CxProgressPhase::Complete, close),
        ];
        assert_eq!(progress.checkpoints().len(), expected.len());
        for (checkpoint, (record_index, item_index, phase, cursor)) in
            progress.checkpoints().iter().zip(expected)
        {
            assert_eq!(checkpoint.record_index, record_index);
            assert_eq!(checkpoint.item_index, item_index);
            assert_eq!(checkpoint.phase, phase);
            assert_eq!(checkpoint.cursor, cursor);
        }
    }

    #[test]
    fn cx_progress_directions_retains_completed_items_and_failure_cursors() {
        let invalid_marker = "|wX:0.0|";
        let failure = invalid_marker.find('X').unwrap();
        let progress = parse_cx_extensions_progress(invalid_marker);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));
        assert!(progress.records().is_empty());
        assert!(progress.checkpoints().is_empty());

        let short_unknown = "|w:|";
        let failure = short_unknown.find(':').unwrap();
        let progress = parse_cx_extensions_progress(short_unknown);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));
        assert!(progress.records().is_empty());

        let malformed_pair = "|wU:0.0,1x.1|";
        let failure = malformed_pair.find('x').unwrap();
        let progress = parse_cx_extensions_progress(malformed_pair);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));
        assert!(matches!(
            &progress.records()[0],
            CxRecord::WedgedBonds(wedges)
                if wedges.len() == 1 && wedges[0].atom == 0 && wedges[0].bond == 0
        ));
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Item);
        assert_eq!(
            progress.checkpoints()[1].cursor,
            malformed_pair.find(',').unwrap()
        );

        let later_overflow = "|wU:0.0,1.4294967296|";
        let overflow = later_overflow.find("4294967296").unwrap();
        let progress = parse_cx_extensions_progress(later_overflow);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), overflow + "4294967296".len());
        assert_eq!(progress.error().map(|error| error.offset), Some(overflow));
        assert_eq!(progress.checkpoints().len(), 2);
        assert!(matches!(
            &progress.records()[0],
            CxRecord::WedgedBonds(wedges) if wedges.len() == 1
        ));

        let missing_stereo_colon = "|ctu|";
        let close = missing_stereo_colon.rfind('|').unwrap();
        let progress = parse_cx_extensions_progress(missing_stereo_colon);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), close);
        assert_eq!(progress.error().map(|error| error.offset), Some(close));
        assert!(progress.records().is_empty());
        assert!(progress.checkpoints().is_empty());

        let stereo_overflow = "|ctu:1,4294967296|";
        let overflow = stereo_overflow.find("4294967296").unwrap();
        let progress = parse_cx_extensions_progress(stereo_overflow);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), overflow + "4294967296".len());
        assert_eq!(progress.error().map(|error| error.offset), Some(overflow));
        assert!(matches!(
            &progress.records()[0],
            CxRecord::DoubleBondStereo(stereo)
                if stereo.stereo == CxDoubleBondStereoKind::Any
                    && stereo.bonds == [1]
        ));
        assert_eq!(progress.checkpoints().len(), 2);
        assert_eq!(progress.checkpoints()[1].item_index, Some(0));

        let missing_outer_pipe = "|wU:0.0";
        let progress = parse_cx_extensions_progress(missing_outer_pipe);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), missing_outer_pipe.len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(missing_outer_pipe.len())
        );
        assert_eq!(progress.checkpoints().len(), 3);
        assert_eq!(progress.checkpoints()[1].phase, CxProgressPhase::Item);
        assert_eq!(progress.checkpoints()[2].phase, CxProgressPhase::Complete);
        assert_eq!(progress.checkpoints()[2].cursor, missing_outer_pipe.len());
    }
}
