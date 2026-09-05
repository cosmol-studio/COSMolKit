use super::*;

pub(super) fn handle_cx_part_and_name(
    state: &mut SmilesBuildState,
    params: &SmilesParseParams,
    cx_part: &str,
    name: &mut String,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION handleCXPartAndName
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: void handleCXPartAndName(RWMol *res, const T &params, const std::string &cxPart,
    // RDKit✔️✔️:                          std::string &name) {
    // RDKit✔️✔️:   if (!res || cxPart.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::string::const_iterator pos = cxPart.cbegin();
    // RDKit✔️✔️:   bool cxfailed = false;
    // RDKit✔️✔️:   if (params.allowCXSMILES) {
    // RDKit✔️✔️:     if (*pos == '|') {
    // RDKit✔️✔️:       try {
    // RDKit✔️✔️:         SmilesParseOps::parseCXExtensions(*res, cxPart, pos);
    // RDKit✔️✔️:       } catch (...) {
    // RDKit✔️✔️:         cxfailed = true;
    // RDKit✔️✔️:         if (params.strictCXSMILES) {
    // RDKit✔️✔️:           throw;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       res->setProp("_CXSMILES_Data", std::string(cxPart.cbegin(), pos));
    // RDKit✔️✔️:     } else if (params.strictCXSMILES && !params.parseName &&
    // RDKit✔️✔️:                pos != cxPart.cend()) {
    // RDKit✔️✔️:       throw RDKit::SmilesParseException(
    // RDKit✔️✔️:           "CXSMILES extension does not start with | and parseName=false");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!cxfailed && params.parseName && pos != cxPart.end()) {
    // RDKit✔️✔️:     std::string nmpart(pos, cxPart.cend());
    // RDKit✔️✔️:     name = boost::trim_copy(nmpart);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION handleCXPartAndName
    if cx_part.is_empty() {
        return Ok(());
    }

    if params.allow_cxsmiles && cx_part.starts_with('|') {
        match parse_cx_extensions(state, cx_part) {
            Ok(consumed) => {
                state.set_property("_CXSMILES_Data", &cx_part[..consumed]);
                if params.parse_name && consumed < cx_part.len() {
                    let name_part = cx_part[consumed..].trim();
                    if !name_part.is_empty() {
                        *name = name_part.to_string();
                    }
                }
                return Ok(());
            }
            Err(error) => {
                if params.strict_cxsmiles {
                    return Err(error);
                }
                state.set_property("_CXSMILES_Data", "");
                return Ok(());
            }
        }
    }

    if params.allow_cxsmiles && params.strict_cxsmiles && !params.parse_name {
        return Err(SmilesParseError::ParseError(
            "CXSMILES extension does not start with | and parseName=false".to_string(),
        ));
    }

    if params.parse_name {
        *name = cx_part.trim().to_string();
        return Ok(());
    }

    Ok(())
}

pub(super) fn parse_cx_extensions(state: &mut SmilesBuildState, ext_text: &str) -> Result<usize, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parseCXExtensions / parser::parse_it
    // RDKit✔️✔️: void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
    // RDKit✔️✔️:                        std::string::const_iterator &first,
    // RDKit✔️✔️:                        unsigned int startAtomIdx, unsigned int startBondIdx) {
    // RDKit✔️✔️:   if (extText.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (extText[0] != '|') {
    // RDKit✔️✔️:     throw RDKit::SmilesParseException(
    // RDKit✔️✔️:         "CXSMILES extension does not start with |");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first = extText.begin();
    // RDKit✔️✔️:   bool ok =
    // RDKit✔️✔️:       parser::parse_it(first, extText.end(), mol, startAtomIdx, startBondIdx);
    // RDKit✔️✔️:   if (!ok) {
    // RDKit✔️✔️:     throw RDKit::SmilesParseException("failure parsing CXSMILES extensions");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   processCXSmilesLabels(mol);
    // RDKit✔️✔️:   mol.clearProp("_cxsmilesLabelsProcessed");
    // RDKit✔️✔️:   mol.clearProp(cxsgTracker);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: bool parse_it(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:               unsigned int startAtomIdx, unsigned int startBondIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != '|') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   while (first < last && *first != '|') {
    // RDKit✔️✔️:     typename Iterator::difference_type length = std::distance(first, last);
    // RDKit✔️✔️:     if (*first == '(') {
    // RDKit✔️✔️:       if (!parse_coords(first, last, mol, startAtomIdx, confIndex++)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == '$') {
    // RDKit✔️✔️:       if (length > 4 && *(first + 1) == '_' && *(first + 2) == 'A' &&
    // RDKit✔️✔️:           *(first + 3) == 'V' && *(first + 4) == ':') {
    // RDKit✔️✔️:         first += 4;
    // RDKit✔️✔️:         if (!parse_atom_values(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (!parse_atom_labels(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (length > 9 && std::string(first, first + 9) == "atomProp:") {
    // RDKit✔️✔️:       first += 9;
    // RDKit✔️✔️:       if (!parse_atom_props(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'C') {
    // RDKit✔️✔️:       if (!parse_coordinate_bonds(first, last, mol, Bond::DATIVE, startAtomIdx,
    // RDKit✔️✔️:                                   startBondIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'H') {
    // RDKit✔️✔️:       if (!parse_coordinate_bonds(first, last, mol, Bond::HYDROGEN,
    // RDKit✔️✔️:                                   startAtomIdx, startBondIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'Z') {
    // RDKit✔️✔️:       if (!parse_zero_bonds(first, last, mol, startAtomIdx, startBondIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == '^') {
    // RDKit✔️✔️:       if (!parse_radicals(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'a' || *first == 'o' ||
    // RDKit✔️✔️:                (*first == '&' && first + 1 < last && first[1] != '#')) {
    // RDKit✔️✔️:       if (!parse_enhanced_stereo(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'r' && first + 1 < last && first[1] == 'b') {
    // RDKit✔️✔️:       if (!parse_ring_bonds(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'L' && first + 1 < last && first[1] == 'N') {
    // RDKit✔️✔️:       if (!parse_linknodes(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
    // RDKit✔️✔️:                first[2] == 'D') {
    // RDKit✔️✔️:       if (!parse_data_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
    // RDKit✔️✔️:                first[2] == 'H') {
    // RDKit✔️✔️:       if (!parse_sgroup_hierarchy(first, last, mol)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'S' && first + 1 < last && first[1] == 'g') {
    // RDKit✔️✔️:       if (!parse_polymer_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'u') {
    // RDKit✔️✔️:       if (!parse_unsaturation(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 's') {
    // RDKit✔️✔️:       if (!parse_substitution(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'm') {
    // RDKit✔️✔️:       if (!parse_variable_attachments(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'w') {
    // RDKit✔️✔️:       if (!parse_wedged_bonds(first, last, mol, startAtomIdx, startBondIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'c' && first + 2 < last && first[1] == 't' &&
    // RDKit✔️✔️:                first[2] == 'u') {
    // RDKit✔️✔️:       if (!parse_doublebond_stereo(first, last, mol, startAtomIdx, startBondIdx,
    // RDKit✔️✔️:                                    Bond::BondStereo::STEREOANY)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'c') {
    // RDKit✔️✔️:       if (!parse_doublebond_stereo(first, last, mol, startAtomIdx, startBondIdx,
    // RDKit✔️✔️:                                    Bond::BondStereo::STEREOCIS)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 't') {
    // RDKit✔️✔️:       if (!parse_doublebond_stereo(first, last, mol, startAtomIdx, startBondIdx,
    // RDKit✔️✔️:                                    Bond::BondStereo::STEREOTRANS)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // if(first < last && *first != '|') ++first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != '|') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parseCXExtensions / parser::parse_it
    if ext_text.is_empty() {
        return Ok(0);
    }
    let bytes = ext_text.as_bytes();
    if bytes.first().copied() != Some(b'|') {
        return Err(SmilesParseError::ParseError(
            "CXSMILES extension does not start with |".to_string(),
        ));
    }
    let mut pos = 1;
    let mut conformer_idx = 0_usize;
    let mut sgroup_idx = 0_usize;
    while pos < bytes.len() && bytes[pos] != b'|' {
        if bytes[pos] == b'(' {
            parse_cx_coords(state, ext_text, &mut pos, conformer_idx)?;
            conformer_idx += 1;
        } else if bytes[pos] == b'$' {
            if ext_text[pos..].starts_with("$_AV:") {
                pos += 4;
                parse_cx_atom_values(state, ext_text, &mut pos)?;
            } else {
                parse_cx_atom_labels(state, ext_text, &mut pos)?;
            }
        } else if ext_text[pos..].starts_with("atomProp:") {
            pos += "atomProp:".len();
            parse_cx_atom_props(state, ext_text, &mut pos)?;
        } else if bytes[pos] == b'C' {
            parse_cx_coordinate_bonds(state, ext_text, &mut pos, BondOrder::Dative)?;
        } else if bytes[pos] == b'H' {
            parse_cx_coordinate_bonds(state, ext_text, &mut pos, BondOrder::Hydrogen)?;
        } else if bytes[pos] == b'Z' {
            parse_cx_zero_bonds(state, ext_text, &mut pos)?;
        } else if bytes[pos] == b'^' {
            parse_cx_radicals(state, ext_text, &mut pos)?;
        } else if matches!(bytes[pos], b'a' | b'o')
            || (bytes[pos] == b'&' && pos + 1 < bytes.len() && bytes[pos + 1] != b'#')
        {
            parse_cx_enhanced_stereo(state, ext_text, &mut pos)?;
        } else if ext_text[pos..].starts_with("rb:") {
            parse_cx_ring_bonds(state, ext_text, &mut pos)?;
        } else if ext_text[pos..].starts_with("LN:") {
            parse_cx_linknodes(state, ext_text, &mut pos)?;
        } else if ext_text[pos..].starts_with("SgD:") {
            parse_cx_data_sgroup(state, ext_text, &mut pos, sgroup_idx)?;
            sgroup_idx += 1;
        } else if ext_text[pos..].starts_with("SgH:") {
            // BEGIN RDKIT CPP FUNCTION parse_sgroup_hierarchy (called from parser::parse_it)
            // RDKit✔️✔️:     } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
            // RDKit✔️✔️:                first[2] == 'H') {
            // RDKit✔️✔️:       if (!parse_sgroup_hierarchy(first, last, mol)) {
            // RDKit✔️✔️:         return false;
            // RDKit✔️✔️:       }
            // END RDKIT CPP FUNCTION call-site in parser::parse_it
            parse_cx_sgroup_hierarchy(state, ext_text, &mut pos)?;
        } else if ext_text[pos..].starts_with("Sg:") {
            // BEGIN RDKIT CPP FUNCTION parse_polymer_sgroup (called from parser::parse_it)
            // RDKit✔️✔️:     } else if (*first == 'S' && first + 1 < last && first[1] == 'g') {
            // RDKit✔️✔️:       if (!parse_polymer_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
            // RDKit✔️✔️:         return false;
            // RDKit✔️✔️:       }
            // END RDKIT CPP FUNCTION call-site in parser::parse_it
            parse_cx_polymer_sgroup(state, ext_text, &mut pos, sgroup_idx)?;
            sgroup_idx += 1;
        } else if bytes[pos] == b'u' {
            parse_cx_unsaturation(state, ext_text, &mut pos)?;
        } else if bytes[pos] == b's' {
            parse_cx_substitution(state, ext_text, &mut pos)?;
        } else if bytes[pos] == b'm' {
            parse_cx_variable_attachments(state, ext_text, &mut pos)?;
        } else if bytes[pos] == b'w' {
            parse_cx_wedged_bonds(state, ext_text, &mut pos)?;
        } else if ext_text[pos..].starts_with("ctu") {
            parse_cx_doublebond_stereo(state, ext_text, &mut pos, BondStereo::Any)?;
        } else if bytes[pos] == b'c' {
            parse_cx_doublebond_stereo(state, ext_text, &mut pos, BondStereo::Cis)?;
        } else if bytes[pos] == b't' {
            parse_cx_doublebond_stereo(state, ext_text, &mut pos, BondStereo::Trans)?;
        } else if bytes[pos] == b',' {
            pos += 1;
        } else {
            // Skip unrecognized CXSMILES extension characters (matching RDKit
            // behavior which silently skips unknown extensions rather than
            // failing). Advance past the single-byte extension token.
            pos += 1;
        }
    }
    if pos >= bytes.len() || bytes[pos] != b'|' {
        return Err(SmilesParseError::ParseError(
            "failure parsing CXSMILES extensions".to_string(),
        ));
    }
    Ok(pos + 1)
}

pub(super) fn parse_cx_coords(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    conformer_idx: usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_coords / hasNonZeroZCoords
    // RDKit✔️✔️: bool parse_coords(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                   unsigned int startAtomIdx, unsigned int confIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != '(') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto *conf = new Conformer(mol.getNumAtoms());
    // RDKit✔️✔️:   mol.addConformer(conf);
    // RDKit✔️✔️:   conf->setId(confIdx);
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   unsigned int atIdx = 0;
    // RDKit✔️✔️:   bool is3D = false;
    // RDKit✔️✔️:   while (first <= last && *first != ')') {
    // RDKit✔️✔️:     RDGeom::Point3D pt;
    // RDKit✔️✔️:     std::string tkn = read_text_to(first, last, ";)");
    // RDKit✔️✔️:     if (VALID_ATIDX(atIdx)) {
    // RDKit✔️✔️:       if (!tkn.empty()) {
    // RDKit✔️✔️:         std::vector<std::string> tokens;
    // RDKit✔️✔️:         boost::split(tokens, tkn, boost::is_any_of(std::string(",")));
    // RDKit✔️✔️:         if (tokens.size() >= 1 && tokens[0].size()) {
    // RDKit✔️✔️:           pt.x = boost::lexical_cast<double>(tokens[0]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (tokens.size() >= 2 && tokens[1].size()) {
    // RDKit✔️✔️:           pt.y = boost::lexical_cast<double>(tokens[1]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (tokens.size() >= 3 && tokens[2].size()) {
    // RDKit✔️✔️:           pt.z = boost::lexical_cast<double>(tokens[2]);
    // RDKit✔️✔️:           is3D = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       conf->setAtomPos(atIdx - startAtomIdx, pt);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atIdx;
    // RDKit✔️✔️:     if (first <= last && *first != ')') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // make sure that the conformer really is 3D!
    // RDKit✔️✔️:   if (is3D && hasNonZeroZCoords(*conf)) {
    // RDKit✔️✔️:     conf->set3D(true);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     conf->set3D(false);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != ')') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: inline bool hasNonZeroZCoords(const Conformer &conf) {
    // RDKit✔️✔️:   constexpr double zeroTol = 1e-3;
    // RDKit✔️✔️:   for (auto p : conf.getPositions()) {
    // RDKit✔️✔️:     if (std::abs(p.z) > zeroTol) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_coords / hasNonZeroZCoords
    expect_byte(ext_text, *pos, b'(')?;
    *pos += 1;
    let atom_count = state.builder.atoms().len();
    let mut coords = vec![[0.0_f64, 0.0_f64, 0.0_f64]; atom_count];
    let mut atom_idx = 0_usize;
    let mut saw_z = false;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b')' {
        let token = read_cx_text_to(ext_text, pos, b";)")?;
        if atom_idx < atom_count && !token.is_empty() {
            let mut pieces = token.split(',');
            if let Some(x) = pieces.next()
                && !x.is_empty()
            {
                coords[atom_idx][0] = x.parse::<f64>().map_err(|_| cx_parse_failure())?;
            }
            if let Some(y) = pieces.next()
                && !y.is_empty()
            {
                coords[atom_idx][1] = y.parse::<f64>().map_err(|_| cx_parse_failure())?;
            }
            if let Some(z) = pieces.next()
                && !z.is_empty()
            {
                coords[atom_idx][2] = z.parse::<f64>().map_err(|_| cx_parse_failure())?;
                saw_z = true;
            }
        }
        atom_idx += 1;
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b')' {
            *pos += 1;
        }
    }
    let is_3d = saw_z && coords.iter().any(|coord| coord[2].abs() > 1e-3);
    state
        .builder
        .add_conformer(Conformer3D::new(conformer_idx, coords, is_3d))
        .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
    expect_byte(ext_text, *pos, b')')?;
    *pos += 1;
    Ok(())
}

pub(super) fn parse_cx_atom_labels(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_atom_labels
    // RDKit✔️✔️: bool parse_atom_labels(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                        unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != '$') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   unsigned int atIdx = 0;
    // RDKit✔️✔️:   while (first <= last && *first != '$') {
    // RDKit✔️✔️:     std::string tkn = read_text_to(first, last, ";$");
    // RDKit✔️✔️:     if (!tkn.empty() && VALID_ATIDX(atIdx)) {
    // RDKit✔️✔️:       mol.getAtomWithIdx(atIdx - startAtomIdx)
    // RDKit✔️✔️:           ->setProp(RDKit::common_properties::atomLabel, tkn);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atIdx;
    // RDKit✔️✔️:     if (first <= last && *first != '$') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != '$') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_atom_labels
    expect_byte(ext_text, *pos, b'$')?;
    *pos += 1;
    let mut atom_idx = 0_usize;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'$' {
        let token = read_cx_text_to(ext_text, pos, b";$")?;
        if !token.is_empty()
            && let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx))
        {
            atom.set_prop("atomLabel", token);
        }
        atom_idx += 1;
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'$' {
            *pos += 1;
        }
    }
    expect_byte(ext_text, *pos, b'$')?;
    *pos += 1;
    Ok(())
}

pub(super) fn parse_cx_atom_values(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_atom_values
    // RDKit✔️✔️: bool parse_atom_values(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                        unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   unsigned int atIdx = 0;
    // RDKit✔️✔️:   while (first <= last && *first != '$') {
    // RDKit✔️✔️:     std::string tkn = read_text_to(first, last, ";$");
    // RDKit✔️✔️:     if (tkn != "" && VALID_ATIDX(atIdx)) {
    // RDKit✔️✔️:       mol.getAtomWithIdx(atIdx)->setProp(RDKit::common_properties::molFileValue,
    // RDKit✔️✔️:                                          tkn);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atIdx;
    // RDKit✔️✔️:     if (first <= last && *first != '$') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != '$') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_atom_values
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    let mut atom_idx = 0_usize;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'$' {
        let token = read_cx_text_to(ext_text, pos, b";$")?;
        if !token.is_empty()
            && let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx))
        {
            atom.set_prop("molFileValue", token);
        }
        atom_idx += 1;
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'$' {
            *pos += 1;
        }
    }
    expect_byte(ext_text, *pos, b'$')?;
    *pos += 1;
    Ok(())
}

pub(super) fn parse_cx_atom_props(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_atom_props
    // RDKit✔️✔️: bool parse_atom_props(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                       unsigned int startAtomIdx) {
    // RDKit✔️✔️:   while (first <= last && *first != '|' && *first != ',') {
    // RDKit✔️✔️:     unsigned int atIdx;
    // RDKit✔️✔️:     if (read_int(first, last, atIdx)) {
    // RDKit✔️✔️:       if (first >= last || *first != '.') {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       std::string pname = read_text_to(first, last, ".");
    // RDKit✔️✔️:       if (!pname.empty()) {
    // RDKit✔️✔️:         if (first >= last || *first != '.') {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         ++first;
    // RDKit✔️✔️:         std::string pval = read_text_to(first, last, ":|,");
    // RDKit✔️✔️:         if (VALID_ATIDX(atIdx) && !pval.empty()) {
    // RDKit✔️✔️:           mol.getAtomWithIdx(atIdx - startAtomIdx)->setProp(pname, pval);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first <= last && *first != '|' && *first != ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (*first != '|') {
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_atom_props
    while *pos < ext_text.len() && !matches!(ext_text.as_bytes()[*pos], b'|' | b',') {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b'.')?;
        *pos += 1;
        let prop_name = read_cx_text_to(ext_text, pos, b".")?;
        expect_byte(ext_text, *pos, b'.')?;
        *pos += 1;
        let prop_value = read_cx_text_to(ext_text, pos, b":|,")?;
        if !prop_name.is_empty()
            && !prop_value.is_empty()
            && let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx))
        {
            atom.set_prop(prop_name, prop_value);
        }
        if *pos < ext_text.len() && !matches!(ext_text.as_bytes()[*pos], b'|' | b',') {
            *pos += 1;
        }
    }
    if *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'|' {
        *pos += 1;
    }
    Ok(())
}

pub(super) fn cx_bond_with_smiles_idx(state: &SmilesBuildState, idx: usize) -> Option<BondId> {
    // BEGIN RDKIT CPP FUNCTION get_bond_with_smiles_idx
    // RDKit✔️✔️: Bond *get_bond_with_smiles_idx(const ROMol &mol, unsigned idx) {
    // RDKit✔️✔️:   for (auto bnd : mol.bonds()) {
    // RDKit✔️✔️:     unsigned int smilesIdx;
    // RDKit✔️✔️:     if (bnd->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
    // RDKit✔️✔️:         smilesIdx == idx) {
    // RDKit✔️✔️:       return bnd;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_bond_with_smiles_idx
    state
        .builder
        .bonds()
        .iter()
        .find(|bond| {
            bond.prop(CXSMILES_BOND_IDX_PROP)
                .and_then(|value| value.parse::<usize>().ok())
                == Some(idx)
        })
        .map(|bond| bond.id())
}

pub(super) fn parse_cx_coordinate_bonds(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    order: BondOrder,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_coordinate_bonds
    // RDKit✔️✔️: bool parse_coordinate_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                             Bond::BondType typ, unsigned int startAtomIdx,
    // RDKit✔️✔️:                             unsigned int startBondIdx) {
    // RDKit✔️✔️:   if (first >= last || (*first != 'C' && *first != 'H')) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int aidx;
    // RDKit✔️✔️:     unsigned int bidx;
    // RDKit✔️✔️:     if (read_int_pair(first, last, aidx, bidx)) {
    // RDKit✔️✔️:       if (VALID_ATIDX(aidx) && VALID_BNDIDX(bidx)) {
    // RDKit✔️✔️:         auto bnd = get_bond_with_smiles_idx(mol, bidx - startBondIdx);
    // RDKit✔️✔️:         if (!bnd || (bnd->getBeginAtomIdx() != aidx - startAtomIdx &&
    // RDKit✔️✔️:                      bnd->getEndAtomIdx() != aidx - startAtomIdx)) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog) << "BOND NOT FOUND! " << bidx
    // RDKit✔️✔️:                                   << " involving atom " << aidx << std::endl;
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         bnd->setBondType(typ);
    // RDKit✔️✔️:         if (bnd->getBeginAtomIdx() != aidx - startAtomIdx) {
    // RDKit✔️✔️:           unsigned int tmp = bnd->getBeginAtomIdx();
    // RDKit✔️✔️:           bnd->setBeginAtomIdx(aidx - startAtomIdx);
    // RDKit✔️✔️:           bnd->setEndAtomIdx(tmp);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_coordinate_bonds
    if *pos >= ext_text.len() || !matches!(ext_text.as_bytes()[*pos], b'C' | b'H') {
        return Err(cx_parse_failure());
    }
    *pos += 1;
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b'.')?;
        *pos += 1;
        let bond_idx = read_cx_usize(ext_text, pos)?;
        if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
            let bond_id = cx_bond_with_smiles_idx(state, bond_idx).ok_or_else(cx_parse_failure)?;
            let (begin, end) = state
                .builder
                .bond(bond_id)
                .map(|bond| (bond.begin(), bond.end()))
                .ok_or_else(cx_parse_failure)?;
            let atom_id = AtomId::new(atom_idx);
            if begin != atom_id && end != atom_id {
                return Err(cx_parse_failure());
            }
            let bond = state.builder.bond_mut(bond_id).ok_or_else(cx_parse_failure)?;
            bond.set_order(order);
            if begin != atom_id {
                bond.set_endpoints(atom_id, begin);
            }
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

pub(super) fn parse_cx_zero_bonds(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_zero_bonds
    // RDKit✔️✔️: bool parse_zero_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                       unsigned int, unsigned int startBondIdx) {
    // RDKit✔️✔️:   // these look like: C1CCCCC~CCCC1 |Z:5|
    // RDKit✔️✔️:   if (first >= last || *first != 'Z') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int bondIdx;
    // RDKit✔️✔️:     if (!read_int(first, last, bondIdx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_BNDIDX(bondIdx)) {
    // RDKit✔️✔️:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!bond) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "bond " << bondIdx
    // RDKit✔️✔️:             << " not found, cannot mark as zero order bond." << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       bond->setBondType(Bond::ZERO);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_zero_bonds
    expect_byte(ext_text, *pos, b'Z')?;
    *pos += 1;
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let bond_idx = read_cx_usize(ext_text, pos)?;
        let bond_id = cx_bond_with_smiles_idx(state, bond_idx).ok_or_else(cx_parse_failure)?;
        state
            .builder
            .bond_mut(bond_id)
            .ok_or_else(cx_parse_failure)?
            .set_order(BondOrder::Zero);
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

pub(super) fn parse_cx_enhanced_stereo(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_enhanced_stereo
    // RDKit✔️✔️: bool parse_enhanced_stereo(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                            unsigned int startAtomIdx) {
    // RDKit✔️✔️:   StereoGroupType group_type = StereoGroupType::STEREO_ABSOLUTE;
    // RDKit✔️✔️:   if (*first == 'a') {
    // RDKit✔️✔️:     group_type = StereoGroupType::STEREO_ABSOLUTE;
    // RDKit✔️✔️:   } else if (*first == 'o') {
    // RDKit✔️✔️:     group_type = StereoGroupType::STEREO_OR;
    // RDKit✔️✔️:   } else if (*first == '&') {
    // RDKit✔️✔️:     group_type = StereoGroupType::STEREO_AND;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // OR and AND groups carry a group number
    // RDKit✔️✔️:   unsigned int group_id = 0;
    // RDKit✔️✔️:   if (group_type != StereoGroupType::STEREO_ABSOLUTE) {
    // RDKit✔️✔️:     read_int(first, last, group_id);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Atom *> atoms;
    // RDKit✔️✔️:   std::vector<Bond *> bonds;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int aidx;
    // RDKit✔️✔️:     if (read_int(first, last, aidx)) {
    // RDKit✔️✔️:       if (VALID_ATIDX(aidx)) {
    // RDKit✔️✔️:         Atom *atom = mol.getAtomWithIdx(aidx - startAtomIdx);
    // RDKit✔️✔️:         if (!atom) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "Atom " << aidx << " not found!" << std::endl;
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atoms.push_back(atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!atoms.empty()) {
    // RDKit✔️✔️:     // we need to do a bit of work to check whether or not we've already seen
    // RDKit✔️✔️:     // this particular StereoGroup (was Github #6050)
    // RDKit✔️✔️:     const auto group_hash =
    // RDKit✔️✔️:         10 * group_id + static_cast<unsigned int>(group_type);
    // RDKit✔️✔️:     std::vector<unsigned int> sgTracker;
    // RDKit✔️✔️:     mol.getPropIfPresent(cxsgTracker, sgTracker);
    // RDKit✔️✔️:     std::vector<StereoGroup> mol_stereo_groups(mol.getStereoGroups());
    // RDKit✔️✔️:     TEST_ASSERT(mol_stereo_groups.size() == sgTracker.size());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto iter = std::find(sgTracker.begin(), sgTracker.end(), group_hash);
    // RDKit✔️✔️:     if (iter != sgTracker.end()) {
    // RDKit✔️✔️:       auto index = iter - sgTracker.begin();
    // RDKit✔️✔️:       auto gAtoms = mol_stereo_groups[index].getAtoms();
    // RDKit✔️✔️:       gAtoms.insert(gAtoms.end(), atoms.begin(), atoms.end());
    // RDKit✔️✔️:       mol_stereo_groups[index] =
    // RDKit✔️✔️:           StereoGroup(mol_stereo_groups[index].getGroupType(),
    // RDKit✔️✔️:                       std::move(gAtoms), std::move(bonds), group_id);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // not seen this before, create a new stereogroup
    // RDKit✔️✔️:       mol_stereo_groups.emplace_back(group_type, std::move(atoms),
    // RDKit✔️✔️:                                      std::move(bonds), group_id);
    // RDKit✔️✔️:       sgTracker.push_back(group_hash);
    // RDKit✔️✔️:       mol.setProp(cxsgTracker, sgTracker);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     mol.setStereoGroups(std::move(mol_stereo_groups));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_enhanced_stereo
    let kind = match ext_text.as_bytes().get(*pos).copied() {
        Some(b'a') => StereoGroupKind::Absolute,
        Some(b'o') => StereoGroupKind::Or,
        Some(b'&') => StereoGroupKind::And,
        _ => return Err(cx_parse_failure()),
    };
    *pos += 1;
    let group_id = if kind == StereoGroupKind::Absolute {
        0
    } else if *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        read_cx_usize(ext_text, pos)? as u32
    } else {
        0
    };
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    let mut atoms = Vec::new();
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        let atom_id = AtomId::new(atom_idx);
        if state.builder.atom_mut(atom_id).is_some() {
            atoms.push(atom_id);
        } else {
            return Err(cx_parse_failure());
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    if atoms.is_empty() {
        return Ok(());
    }
    let kind_code = match kind {
        StereoGroupKind::Absolute => 1,
        StereoGroupKind::Or => 2,
        StereoGroupKind::And => 3,
    };
    if let Some(index) = state.cx_stereo_group_tracker.get(&(kind_code, group_id)).copied() {
        if let Some(group) = state.builder.stereo_groups_mut().get_mut(index) {
            for atom in atoms {
                group.push_atom(atom);
            }
        }
    } else {
        let index = state.builder.stereo_groups_mut().len();
        state
            .builder
            .add_stereo_group(StereoGroup::new(kind, atoms, Vec::new()).with_id(group_id))
            .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
        state.cx_stereo_group_tracker.insert((kind_code, group_id), index);
    }
    Ok(())
}

pub(super) fn expand_cx_atom_query(state: &mut SmilesBuildState, atom_idx: usize, predicate: AtomQueryPredicate) {
    if let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx)) {
        let next = QueryNode::predicate(predicate);
        let combined = match atom.query().cloned() {
            Some(QueryNode::And(mut children)) => {
                children.push(next);
                QueryNode::And(children)
            }
            Some(existing) => QueryNode::and(vec![existing, next]),
            None => next,
        };
        atom.set_query(Some(combined));
    }
}

pub(super) fn parse_cx_unsaturation(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_unsaturation
    // RDKit✔️✔️: bool parse_unsaturation(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                         unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first + 1 >= last || *first != 'u') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int idx;
    // RDKit✔️✔️:     if (!read_int(first, last, idx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(idx)) {
    // RDKit✔️✔️:       auto atom = mol.getAtomWithIdx(idx - startAtomIdx);
    // RDKit✔️✔️:       if (!atom->hasQuery()) {
    // RDKit✔️✔️:         atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       atom->expandQuery(makeAtomUnsaturatedQuery(), Queries::COMPOSITE_AND);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_unsaturation
    expect_byte(ext_text, *pos, b'u')?;
    *pos += 1;
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expand_cx_atom_query(state, atom_idx, AtomQueryPredicate::IsUnsaturated);
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

pub(super) fn parse_cx_ring_bonds(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_ring_bonds
    // RDKit✔️✔️: bool parse_ring_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                       unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != 'r' || first + 1 >= last ||
    // RDKit✔️✔️:       *(first + 1) != 'b' || first + 2 >= last || *(first + 2) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 3;
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int n1;
    // RDKit✔️✔️:     if (!read_int(first, last, n1)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // check that we can read at least two more characters:
    // RDKit✔️✔️:     if (first + 1 >= last || *first != ':') {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     unsigned int n2;
    // RDKit✔️✔️:     bool gt = false;
    // RDKit✔️✔️:     if (*first == '*') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       n2 = 0xDEADBEEF;
    // RDKit✔️✔️:       if (VALID_ATIDX(n1)) {
    // RDKit✔️✔️:         mol.setProp(common_properties::_NeedsQueryScan, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if (!read_int(first, last, n2)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       switch (n2) {
    // RDKit✔️✔️:         case 0:
    // RDKit✔️✔️:         case 2:
    // RDKit✔️✔️:         case 3:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 4:
    // RDKit✔️✔️:           gt = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "unrecognized rb value: " << n2 << std::endl;
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(n1)) {
    // RDKit✔️✔️:       auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
    // RDKit✔️✔️:       if (!atom->hasQuery()) {
    // RDKit✔️✔️:         atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!gt) {
    // RDKit✔️✔️:         atom->expandQuery(makeAtomRingBondCountQuery(n2),
    // RDKit✔️✔️:                           Queries::COMPOSITE_AND);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         auto q = static_cast<ATOM_EQUALS_QUERY *>(new ATOM_LESSEQUAL_QUERY);
    // RDKit✔️✔️:         q->setVal(n2);
    // RDKit✔️✔️:         q->setDescription("AtomRingBondCount");
    // RDKit✔️✔️:         q->setDataFunc(queryAtomRingBondCount);
    // RDKit✔️✔️:         atom->expandQuery(q, Queries::COMPOSITE_AND);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_ring_bonds
    if !ext_text[*pos..].starts_with("rb:") {
        return Err(cx_parse_failure());
    }
    *pos += 3;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b':')?;
        *pos += 1;
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'*' {
            *pos += 1;
            if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
                state.set_property("_NeedsQueryScan", "1");
                expand_cx_atom_query(
                    state,
                    atom_idx,
                    AtomQueryPredicate::RingBondCount(crate::search::query::QUERY_SCAN_MAGIC_VALUE),
                );
            }
        } else {
            let ring_bonds = read_cx_usize(ext_text, pos)?;
            let predicate = match ring_bonds {
                0 | 2 | 3 => AtomQueryPredicate::RingBondCount(ring_bonds as u32),
                4 => AtomQueryPredicate::RingBondCountLessEqual(4),
                _ => return Err(cx_parse_failure()),
            };
            expand_cx_atom_query(state, atom_idx, predicate);
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

pub(super) fn atom_neighbors(state: &SmilesBuildState, atom: AtomId) -> Vec<AtomId> {
    state
        .builder
        .bonds()
        .iter()
        .filter_map(|bond| {
            if bond.begin() == atom {
                Some(bond.end())
            } else if bond.end() == atom {
                Some(bond.begin())
            } else {
                None
            }
        })
        .collect()
}

pub(super) fn parse_cx_linknodes(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_linknodes
    // RDKit✔️✔️: bool parse_linknodes(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                      unsigned int startAtomIdx) {
    // RDKit✔️✔️:   // these look like: |LN:1:1.3.2.6,4:1.4.3.6|
    // RDKit✔️✔️:   // that's two records:
    // RDKit✔️✔️:   //   1:1.3.2.6: 1-3 repeats, atom 1-2, 1-6
    // RDKit✔️✔️:   //   4:1.4.3.6: 1-4 repeats, atom 4-3, 4-6
    // RDKit✔️✔️:   // which maps to the property value "1 3 2 2 3 2 7|1 4 2 5 4 5 7"
    // RDKit✔️✔️:   // If the linking atom only has two neighbors then the outer atom
    // RDKit✔️✔️:   // specification (the last two digits) can be left out. So for a molecule
    // RDKit✔️✔️:   // where atom 1 has bonds only to atoms 2 and 6 we could have
    // RDKit✔️✔️:   // |LN:1:1.3|
    // RDKit✔️✔️:   // instead of
    // RDKit✔️✔️:   // |LN:1:1.3.2.6|
    // RDKit✔️✔️:   if (first >= last || *first != 'L' || first + 1 >= last ||
    // RDKit✔️✔️:       *(first + 1) != 'N' || first + 2 >= last || *(first + 2) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 3;
    // RDKit✔️✔️:   std::string accum = "";
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int atidx;
    // RDKit✔️✔️:     if (!read_int(first, last, atidx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // check that we can read at least two more characters:
    // RDKit✔️✔️:     if (first + 1 >= last || *first != ':') {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     unsigned int startReps;
    // RDKit✔️✔️:     if (!read_int(first, last, startReps)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first + 1 >= last || *first != '.') {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     unsigned int endReps;
    // RDKit✔️✔️:     if (!read_int(first, last, endReps)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned int idx1;
    // RDKit✔️✔️:     unsigned int idx2;
    // RDKit✔️✔️:     if (first < last && *first == '.') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       if (!read_int(first, last, idx1)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       if (!read_int(first, last, idx2)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (VALID_ATIDX(atidx) &&
    // RDKit✔️✔️:                mol.getAtomWithIdx(atidx - startAtomIdx)->getDegree() == 2) {
    // RDKit✔️✔️:       auto nbrs =
    // RDKit✔️✔️:           mol.getAtomNeighbors(mol.getAtomWithIdx(atidx - startAtomIdx));
    // RDKit✔️✔️:       idx1 = *nbrs.first;
    // RDKit✔️✔️:       nbrs.first++;
    // RDKit✔️✔️:       idx2 = *nbrs.first;
    // RDKit✔️✔️:     } else if (VALID_ATIDX(atidx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(atidx)) {
    // RDKit✔️✔️:       if (!accum.empty()) {
    // RDKit✔️✔️:         accum += "|";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       accum += (boost::format("%d %d 2 %d %d %d %d") % startReps % endReps %
    // RDKit✔️✔️:                 (atidx - startAtomIdx + 1) % (idx1 - startAtomIdx + 1) %
    // RDKit✔️✔️:                 (atidx - startAtomIdx + 1) % (idx2 - startAtomIdx + 1))
    // RDKit✔️✔️:                    .str();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!accum.empty()) {
    // RDKit✔️✔️:     mol.setProp(common_properties::molFileLinkNodes, accum);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_linknodes
    if !ext_text[*pos..].starts_with("LN:") {
        return Err(cx_parse_failure());
    }
    *pos += 3;
    let mut records = Vec::new();
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b':')?;
        *pos += 1;
        let start_reps = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b'.')?;
        *pos += 1;
        let end_reps = read_cx_usize(ext_text, pos)?;
        let atom_id = AtomId::new(atom_idx);
        let (idx1, idx2) = if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'.' {
            *pos += 1;
            let idx1 = read_cx_usize(ext_text, pos)?;
            expect_byte(ext_text, *pos, b'.')?;
            *pos += 1;
            let idx2 = read_cx_usize(ext_text, pos)?;
            (idx1, idx2)
        } else if state.builder.atom_mut(atom_id).is_some() {
            let neighbors = atom_neighbors(state, atom_id);
            if neighbors.len() != 2 {
                return Err(cx_parse_failure());
            }
            (neighbors[0].index(), neighbors[1].index())
        } else {
            (0, 0)
        };
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
        if state.builder.atom_mut(atom_id).is_some() {
            records.push(format!(
                "{} {} 2 {} {} {} {}",
                start_reps,
                end_reps,
                atom_idx + 1,
                idx1 + 1,
                atom_idx + 1,
                idx2 + 1
            ));
        }
    }
    if !records.is_empty() {
        state.set_property("_MolFileLinkNodes", &records.join("|"));
    }
    Ok(())
}

pub(super) fn parse_cx_data_sgroup_attr(
    sgroup: &mut SubstanceGroup,
    ext_text: &str,
    pos: &mut usize,
    field_name: &str,
    field_is_array: bool,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_data_sgroup_attr
    // RDKit✔️✔️: void parse_data_sgroup_attr(Iterator &first, Iterator last,
    // RDKit✔️✔️:                             SubstanceGroup &sgroup, bool keepSGroup,
    // RDKit✔️✔️:                             std::string fieldName, bool fieldIsArray = false) {
    // RDKit✔️✔️:   PRECONDITION(first < last, "parse_data_sgroup_attr: first >= last");
    // RDKit✔️✔️:   if (first != last && *first != '|') {
    // RDKit✔️✔️:     std::string data = read_text_to(first, last, ":");
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     if (!data.empty() && keepSGroup) {
    // RDKit✔️✔️:       if (fieldIsArray) {
    // RDKit✔️✔️:         std::vector<std::string> dataFields = {data};
    // RDKit✔️✔️:         sgroup.setProp(fieldName, dataFields);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         sgroup.setProp(fieldName, data);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_data_sgroup_attr
    if *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'|' {
        let data = read_cx_text_to(ext_text, pos, b":")?;
        if *pos < ext_text.len() {
            *pos += 1;
        }
        if !data.is_empty() {
            if field_is_array {
                sgroup.set_prop(field_name, data.clone());
                sgroup.push_data_field(data);
            } else {
                sgroup.set_prop(field_name, data);
            }
        }
    }
    Ok(())
}

pub(super) fn parse_cx_data_sgroup(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    sgroup_idx: usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_data_sgroup
    // RDKit✔️✔️: bool parse_data_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                        unsigned int startAtomIdx, unsigned int nSGroups) {
    // RDKit✔️✔️:   // these look like: |SgD:2,1:FIELD:info::::|
    // RDKit✔️✔️:   // example from CXSMILES docs:
    // RDKit✔️✔️:   //    SgD:3,2,1,0:name:data:like:unit:t:(1.,1.)
    // RDKit✔️✔️:   // the fields are:
    // RDKit✔️✔️:   //    SgD:[atom indices]:[field name]:[data value]:[query
    // RDKit✔️✔️:   //    operator]:[unit]:[tag]:[coords]
    // RDKit✔️✔️:   //   coords are (-1) if atomic coordinates are present
    // RDKit✔️✔️:   if (first >= last || *first != 'S' || first + 3 >= last ||
    // RDKit✔️✔️:       *(first + 1) != 'g' || *(first + 2) != 'D' || *(first + 3) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 4;
    // RDKit✔️✔️:   std::vector<unsigned int> atoms;
    // RDKit✔️✔️:   if (!read_int_list(first, last, atoms)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   SubstanceGroup sgroup(&mol, std::string("DAT"));
    // RDKit✔️✔️:   sgroup.setProp(cxsmilesindex, nSGroups);
    // RDKit✔️✔️:   bool keepSGroup = false;
    // RDKit✔️✔️:   for (auto idx : atoms) {
    // RDKit✔️✔️:     if (VALID_ATIDX(idx)) {
    // RDKit✔️✔️:       keepSGroup = true;
    // RDKit✔️✔️:       sgroup.addAtomWithIdx(idx - startAtomIdx);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDNAME");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // FIX:
    // RDKit✔️✔️:   if (keepSGroup) {
    // RDKit✔️✔️:     sgroup.setProp("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "DATAFIELDS", true);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "QUERYOP");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDINFO");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDTAG");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (first < last && *first == '(') {
    // RDKit✔️✔️:     // FIX
    // RDKit✔️✔️:     std::string coords = read_text_to(first, last, ")");
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     if (keepSGroup) {
    // RDKit✔️✔️:       sgroup.setProp("COORDS", coords);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // the label processing can destroy sgroup info, so do that now
    // RDKit✔️✔️:   // (the function will immediately return if already called)
    // RDKit✔️✔️:   if (keepSGroup) {
    // RDKit✔️✔️:     processCXSmilesLabels(mol);
    // RDKit✔️✔️:     sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
    // RDKit✔️✔️:     addSubstanceGroup(mol, sgroup);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_data_sgroup
    if !ext_text[*pos..].starts_with("SgD:") {
        return Err(cx_parse_failure());
    }
    *pos += 4;
    let mut atoms = Vec::new();
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
            atoms.push(AtomId::new(atom_idx));
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        } else {
            break;
        }
    }
    let mut sgroup = SubstanceGroup::new(SubstanceGroupId::new(sgroup_idx), SubstanceGroupKind::Data);
    sgroup.set_prop("_cxsmilesindex", sgroup_idx.to_string());
    let keep_sgroup = !atoms.is_empty();
    for atom in atoms {
        sgroup.push_atom(atom);
    }
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    parse_cx_data_sgroup_attr(&mut sgroup, ext_text, pos, "FIELDNAME", false)?;
    if keep_sgroup {
        sgroup.set_prop("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0");
    }
    parse_cx_data_sgroup_attr(&mut sgroup, ext_text, pos, "DATAFIELDS", true)?;
    parse_cx_data_sgroup_attr(&mut sgroup, ext_text, pos, "QUERYOP", false)?;
    parse_cx_data_sgroup_attr(&mut sgroup, ext_text, pos, "FIELDINFO", false)?;
    parse_cx_data_sgroup_attr(&mut sgroup, ext_text, pos, "FIELDTAG", false)?;
    if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'(' {
        let coords = read_cx_text_to(ext_text, pos, b")")?;
        if *pos < ext_text.len() {
            *pos += 1;
        }
        if keep_sgroup {
            sgroup.set_prop("COORDS", coords);
        }
    }
    if keep_sgroup {
        sgroup.set_prop("index", (state.builder.substance_groups_len() + 1).to_string());
        state
            .builder
            .add_substance_group(sgroup)
            .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
    }
    Ok(())
}

pub(super) fn parse_cx_sgroup_hierarchy(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_sgroup_hierarchy (CXSmilesOps.cpp)
    // RDKit✔️✔️: bool parse_sgroup_hierarchy(Iterator &first, Iterator last, RDKit::RWMol &mol) {
    // RDKit✔️✔️:   // these look like: |SgH:1:0|
    // RDKit✔️✔️:   // from CXSMILES docs:
    // RDKit✔️✔️:   //    SgH:parentSgroupIndex1:childSgroupIndex1.childSgroupIndex2,parentSgroupIndex2:childSgroupIndex1
    // RDKit✔️✔️:   if (first >= last || *first != 'S' || first + 3 >= last ||
    // RDKit✔️✔️:       *(first + 1) != 'g' || *(first + 2) != 'H' || *(first + 3) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 4;
    // RDKit✔️✔️:   auto &sgs = getSubstanceGroups(mol);
    // RDKit✔️✔️:   while (true) {
    // RDKit✔️✔️:     unsigned int parentId;
    // RDKit✔️✔️:     if (!read_int(first, last, parentId)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bool validParent = true;
    // RDKit✔️✔️:     auto psg = find_matching_sgroup(sgs, parentId);
    // RDKit✔️✔️:     if (psg == sgs.end()) {
    // RDKit✔️✔️:       validParent = false;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       psg->getPropIfPresent(\"index\", parentId);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first <= last && *first == ':') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       std::vector<unsigned int> children;
    // RDKit✔️✔️:       if (!read_int_list(first, last, children, '.')) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (validParent) {
    // RDKit✔️✔️:         for (auto childId : children) {
    // RDKit✔️✔️:           if (childId >= sgs.size()) {
    // RDKit✔️✔️:             throw SmilesParseException("child id references non-existent SGroup");
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           auto csg = find_matching_sgroup(sgs, childId);
    // RDKit✔️✔️:           if (csg != sgs.end()) {
    // RDKit✔️✔️:             unsigned int cid;
    // RDKit✔️✔️:             csg->getProp("index", cid);
    // RDKit✔️✔️:             csg->setProp("PARENT", parentId);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (first <= last && *first == ',') {
    // RDKit✔️✔️:         ++first;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_sgroup_hierarchy
    if !ext_text[*pos..].starts_with("SgH:") {
        return Err(cx_parse_failure());
    }
    *pos += 4;
    // Build a mapping from _cxsmilesindex -> (position, "index" value)
    let cx_index_to_parent: Vec<(usize, usize)> = state
        .builder
        .substance_groups()
        .iter()
        .enumerate()
        .filter_map(|(pos, sg)| {
            sg.props().get("_cxsmilesindex").and_then(|cx_idx| {
                let idx_val = sg
                    .props()
                    .get("index")
                    .and_then(|v| v.parse::<usize>().ok())
                    .unwrap_or(pos);
                cx_idx.parse::<usize>().ok().map(|parsed| (parsed, idx_val))
            })
        })
        .collect();
    loop {
        let parent_cx_idx = read_cx_usize(ext_text, pos)?;
        let parent_actual_idx = cx_index_to_parent
            .iter()
            .find(|(cx, _)| *cx == parent_cx_idx)
            .map(|(_, idx)| *idx);
        expect_byte(ext_text, *pos, b':')?;
        *pos += 1;
        // Read child indices separated by '.'
        let mut children = Vec::new();
        while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
            let child_cx_idx = read_cx_usize(ext_text, pos)?;
            children.push(child_cx_idx);
            if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'.' {
                *pos += 1;
            } else {
                break;
            }
        }
        // Set PARENT prop on each child
        if let Some(actual_parent_idx) = parent_actual_idx {
            for child_cx_idx in &children {
                if *child_cx_idx >= state.builder.substance_groups_len() {
                    return Err(SmilesParseError::ParseError(
                        "child id references non-existent SGroup".to_string(),
                    ));
                }
                if let Some((child_pos, _)) = cx_index_to_parent.iter().find(|(cx, _)| *cx == *child_cx_idx) {
                    if let Some(child_sg) = state.builder.substance_group_mut(*child_pos) {
                        child_sg.set_prop("PARENT", actual_parent_idx.to_string());
                    }
                }
            }
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        } else {
            break;
        }
    }
    Ok(())
}

pub(super) fn parse_cx_polymer_sgroup(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    sgroup_idx: usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_polymer_sgroup (CXSmilesOps.cpp)
    // RDKit✔️✔️: bool parse_polymer_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                           unsigned int startAtomIdx, unsigned int nSGroups) {
    // RDKit✔️✔️:   // format: |Sg:type:atomIndices:subscript:superscript:headCrossings:tailCrossings:|
    // RDKit✔️✔️:   // type codes: n->SRU, mon->MON, mer->MER, co/xl/mod/mix/f/any/gen/c/grf/alt/ran/blk
    // RDKit✔️✔️:   if (first >= last || *first != 'S' || first + 2 >= last ||
    // RDKit✔️✔️:       *(first + 1) != 'g' || *(first + 2) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 3;
    // RDKit✔️✔️:   const auto type_code = read_text_to(first, last, ":");
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   const auto type = sgroupTypemap.find(type_code);
    // RDKit✔️✔️:   if (type == sgroupTypemap.end()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool keepSGroup = false;
    // RDKit✔️✔️:   SubstanceGroup sgroup(&mol, type->second);
    // RDKit✔️✔️:   sgroup.setProp(cxsmilesindex, nSGroups);
    // RDKit✔️✔️:   if (type_code == "alt") {
    // RDKit✔️✔️:     sgroup.setProp("SUBTYPE", "ALT");
    // RDKit✔️✔️:   } else if (type_code == "ran") {
    // RDKit✔️✔️:     sgroup.setProp("SUBTYPE", "RAN");
    // RDKit✔️✔️:   } else if (type_code == "blk") {
    // RDKit✔️✔️:     sgroup.setProp("SUBTYPE", "BLO");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<unsigned int> atoms;
    // RDKit✔️✔️:   if (!read_int_list(first, last, atoms)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto idx : atoms) {
    // RDKit✔️✔️:     if (VALID_ATIDX(idx)) {
    // RDKit✔️✔️:       sgroup.addAtomWithIdx(idx - startAtomIdx);
    // RDKit✔️✔️:       keepSGroup = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<unsigned int> headCrossing;
    // RDKit✔️✔️:   std::vector<unsigned int> tailCrossing;
    // RDKit✔️✔️:   if (first <= last && *first == ':') {
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     std::string subscript = read_text_to(first, last, ":|");
    // RDKit✔️✔️:     if (keepSGroup && !subscript.empty()) {
    // RDKit✔️✔️:       sgroup.setProp("LABEL", subscript);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first <= last && *first == ':') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       std::string superscript = read_text_to(first, last, ":|,");
    // RDKit✔️✔️:       if (keepSGroup && !superscript.empty()) {
    // RDKit✔️✔️:         sgroup.setProp("CONNECT", superscript);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // ... headCrossings and tailCrossings follow, omitted for brevity
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (keepSGroup) {
    // RDKit✔️✔️:     processCXSmilesLabels(mol);
    // RDKit✔️✔️:     finalizePolymerSGroup(mol, sgroup);
    // RDKit✔️✔️:     sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
    // RDKit✔️✔️:     addSubstanceGroup(mol, sgroup);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_polymer_sgroup
    if !ext_text[*pos..].starts_with("Sg:") {
        return Err(cx_parse_failure());
    }
    *pos += 3;
    let type_code = read_cx_text_to(ext_text, pos, b":")?;
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    let Some(kind) = cx_sgroup_type_to_kind(&type_code) else {
        return Err(cx_parse_failure());
    };
    let mut keep_sgroup = false;
    let mut sgroup = SubstanceGroup::new(SubstanceGroupId::new(sgroup_idx), kind);
    sgroup.set_prop("_cxsmilesindex", sgroup_idx.to_string());
    // Set subtype for alt/ran/blk
    match type_code.as_str() {
        "alt" => sgroup.set_prop("SUBTYPE", "ALT"),
        "ran" => sgroup.set_prop("SUBTYPE", "RAN"),
        "blk" => sgroup.set_prop("SUBTYPE", "BLO"),
        _ => {}
    }
    // Read atom indices
    let mut atoms = Vec::new();
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
            atoms.push(AtomId::new(atom_idx));
            keep_sgroup = true;
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        } else {
            break;
        }
    }
    for atom in &atoms {
        sgroup.push_atom(*atom);
    }
    // Read subscript, superscript, headCrossings, tailCrossings
    if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b':' {
        *pos += 1;
        let subscript = read_cx_text_to(ext_text, pos, b":|")?;
        if keep_sgroup && !subscript.is_empty() {
            sgroup.set_prop("LABEL", &subscript);
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b':' {
            *pos += 1;
            let superscript = read_cx_text_to(ext_text, pos, b":|,")?;
            if keep_sgroup && !superscript.is_empty() {
                sgroup.set_prop("CONNECT", &superscript);
            }
            // Head crossing bonds
            if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b':' {
                *pos += 1;
                let mut head_crossing = Vec::new();
                while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
                    let cidx = read_cx_usize(ext_text, pos)?;
                    head_crossing.push(cidx);
                    if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
                        *pos += 1;
                    } else {
                        break;
                    }
                }
                if keep_sgroup && !head_crossing.is_empty() {
                    sgroup.set_prop(
                        "_headCrossings",
                        head_crossing
                            .iter()
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                            .join(","),
                    );
                }
                // Tail crossing bonds
                if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b':' {
                    *pos += 1;
                    let mut tail_crossing = Vec::new();
                    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
                        let cidx = read_cx_usize(ext_text, pos)?;
                        tail_crossing.push(cidx);
                        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
                            *pos += 1;
                        } else {
                            break;
                        }
                    }
                    if keep_sgroup && !tail_crossing.is_empty() {
                        sgroup.set_prop(
                            "_tailCrossings",
                            tail_crossing
                                .iter()
                                .map(|v| v.to_string())
                                .collect::<Vec<_>>()
                                .join(","),
                        );
                    }
                }
            }
        }
    }
    if keep_sgroup {
        sgroup.set_prop("index", (state.builder.substance_groups_len() + 1).to_string());
        state
            .builder
            .add_substance_group(sgroup)
            .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
    }
    Ok(())
}

pub(super) fn cx_sgroup_type_to_kind(type_code: &str) -> Option<SubstanceGroupKind> {
    // RDKit✔️✔️: sgroupTypemap: n->SRU, mon->MON, mer->MER, co->COP, xl->CRO,
    // RDKit✔️✔️: mod->MOD, mix->MIX, f->FOR, any->ANY, gen->GEN, c->COM,
    // RDKit✔️✔️: grf->GRA, alt/ran/blk->COP (with subtype)
    // END RDKIT CPP map lookup
    match type_code {
        "n" => Some(SubstanceGroupKind::StructuralRepeatUnit),
        "mon" => Some(SubstanceGroupKind::Monomer),
        "mer" => Some(SubstanceGroupKind::Mer),
        "co" => Some(SubstanceGroupKind::Copolymer),
        "xl" => Some(SubstanceGroupKind::Crosslink),
        "mod" => Some(SubstanceGroupKind::Modification),
        "mix" => Some(SubstanceGroupKind::MixtureComponent),
        "f" => Some(SubstanceGroupKind::Formulation),
        "any" => Some(SubstanceGroupKind::AnyPolymer),
        "gen" => Some(SubstanceGroupKind::Generic("GEN".to_string())),
        "c" => Some(SubstanceGroupKind::Generic("COM".to_string())),
        "grf" => Some(SubstanceGroupKind::Graft),
        "alt" | "ran" | "blk" => Some(SubstanceGroupKind::Copolymer),
        _ => None,
    }
}

pub(super) fn parse_cx_substitution(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_substitution
    // RDKit✔️✔️: bool parse_substitution(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                         unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != 's' || first + 1 >= last ||
    // RDKit✔️✔️:       *(first + 1) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 2;
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int n1;
    // RDKit✔️✔️:     if (!read_int(first, last, n1)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // check that we can read at least two more characters:
    // RDKit✔️✔️:     if (first + 1 >= last || *first != ':') {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     unsigned int n2;
    // RDKit✔️✔️:     if (*first == '*') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       n2 = 0xDEADBEEF;
    // RDKit✔️✔️:       if (VALID_ATIDX(n1)) {
    // RDKit✔️✔️:         mol.setProp(common_properties::_NeedsQueryScan, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if (!read_int(first, last, n2)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(n1)) {
    // RDKit✔️✔️:       auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
    // RDKit✔️✔️:       if (!atom->hasQuery()) {
    // RDKit✔️✔️:         atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       atom->expandQuery(makeAtomNonHydrogenDegreeQuery(n2),
    // RDKit✔️✔️:                         Queries::COMPOSITE_AND);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_substitution
    if !ext_text[*pos..].starts_with("s:") {
        return Err(cx_parse_failure());
    }
    *pos += 2;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b':')?;
        *pos += 1;
        let non_hydrogen_degree = if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'*' {
            *pos += 1;
            if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
                state.set_property("_NeedsQueryScan", "1");
            }
            crate::search::query::QUERY_SCAN_MAGIC_VALUE
        } else {
            read_cx_usize(ext_text, pos)? as u32
        };
        expand_cx_atom_query(
            state,
            atom_idx,
            AtomQueryPredicate::NonHydrogenDegree(non_hydrogen_degree),
        );
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

pub(super) fn bond_can_have_direction(order: BondOrder) -> bool {
    // BEGIN RDKIT CPP FUNCTION canHaveDirection
    // RDKit✔️✔️: inline bool canHaveDirection(const Bond &bond) {
    // RDKit✔️✔️:   auto bondType = bond.getBondType();
    // RDKit✔️✔️:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION canHaveDirection
    matches!(order, BondOrder::Single | BondOrder::Aromatic)
}

pub(super) fn parse_cx_wedged_bonds(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_wedged_bonds
    // RDKit✔️✔️: bool parse_wedged_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                         unsigned int startAtomIdx, unsigned int startBondIdx) {
    // RDKit✔️✔️:   // these look like: CC(O)Cl |w:1.0|
    // RDKit✔️✔️:   // also wD and wU for down and up wedges.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // We do not end up using this to set stereochemistry, but the relevant bond
    // RDKit✔️✔️:   // properties are set in case client code wants to do something with the
    // RDKit✔️✔️:   // information.
    // RDKit✔️✔️:   if (first >= last || *first != 'w' || first + 1 >= last) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   Bond::BondDir state = Bond::BondDir::NONE;
    // RDKit✔️✔️:   unsigned int cfg = 0;
    // RDKit✔️✔️:   switch (*first) {
    // RDKit✔️✔️:     case ':':
    // RDKit✔️✔️:       state = Bond::BondDir::UNKNOWN;
    // RDKit✔️✔️:       cfg = 2;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 'U':
    // RDKit✔️✔️:       state = Bond::BondDir::BEGINWEDGE;
    // RDKit✔️✔️:       cfg = 1;
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 'D':
    // RDKit✔️✔️:       state = Bond::BondDir::BEGINDASH;
    // RDKit✔️✔️:       cfg = 3;
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (state == Bond::BondDir::NONE || first >= last || first + 1 >= last ||
    // RDKit✔️✔️:       *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int atomIdx;
    // RDKit✔️✔️:     if (!read_int(first, last, atomIdx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == '.') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog) << "improperly formatted w block" << std::endl;
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned int bondIdx;
    // RDKit✔️✔️:     if (!read_int(first, last, bondIdx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (VALID_ATIDX(atomIdx) && VALID_BNDIDX(bondIdx)) {
    // RDKit✔️✔️:       auto atom = mol.getAtomWithIdx(atomIdx - startAtomIdx);
    // RDKit✔️✔️:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!bond) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "bond " << bondIdx << " not found, wedge from atom " << atomIdx
    // RDKit✔️✔️:             << " cannot be applied." << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // we can't set wedging twice:
    // RDKit✔️✔️:       if (bond->hasProp(common_properties::_MolFileBondCfg)) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "w block attempts to set wedging on bond " << bond->getIdx()
    // RDKit✔️✔️:             << " more than once." << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // first things first, the atom needs to be the start atom of the bond for
    // RDKit✔️✔️:       // any of this to make sense
    // RDKit✔️✔️:       if (atom->getIdx() != bond->getBeginAtomIdx()) {
    // RDKit✔️✔️:         if (atom->getIdx() != bond->getEndAtomIdx()) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "atom " << atomIdx << " is not associated with bond "
    // RDKit✔️✔️:               << bondIdx << "(" << bond->getBeginAtomIdx() + startAtomIdx << "-"
    // RDKit✔️✔️:               << bond->getEndAtomIdx() + startAtomIdx << ")"
    // RDKit✔️✔️:               << " in w block" << std::endl;
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         auto eidx = bond->getBeginAtomIdx();
    // RDKit✔️✔️:         bond->setBeginAtomIdx(atom->getIdx());
    // RDKit✔️✔️:         bond->setEndAtomIdx(eidx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       bond->setProp(common_properties::_MolFileBondCfg, cfg);
    // RDKit✔️✔️:       bond->setBondDir(state);
    // RDKit✔️✔️:       if (cfg == 2 && canHaveDirection(*bond)) {
    // RDKit✔️✔️:         bond->getBeginAtom()->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit✔️✔️:         mol.setProp(detail::_needsDetectBondStereo, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if ((cfg == 1 || cfg == 3) && canHaveDirection(*bond)) {
    // RDKit✔️✔️:         mol.setProp(detail::_needsDetectAtomStereo, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_wedged_bonds
    expect_byte(ext_text, *pos, b'w')?;
    *pos += 1;
    let (direction, cfg) = match ext_text.as_bytes().get(*pos).copied() {
        Some(b':') => (BondDirection::Unknown, "2"),
        Some(b'U') => {
            *pos += 1;
            (BondDirection::BeginWedge, "1")
        }
        Some(b'D') => {
            *pos += 1;
            (BondDirection::BeginDash, "3")
        }
        _ => return Err(cx_parse_failure()),
    };
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b'.')?;
        *pos += 1;
        let bond_idx = read_cx_usize(ext_text, pos)?;
        if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
            let bond_id = cx_bond_with_smiles_idx(state, bond_idx).ok_or_else(cx_parse_failure)?;
            let (begin, end, order, has_cfg) = state
                .builder
                .bond(bond_id)
                .map(|bond| {
                    (
                        bond.begin(),
                        bond.end(),
                        bond.order(),
                        bond.prop("_MolFileBondCfg").is_some(),
                    )
                })
                .ok_or_else(cx_parse_failure)?;
            if has_cfg {
                return Err(cx_parse_failure());
            }
            let atom_id = AtomId::new(atom_idx);
            if begin != atom_id && end != atom_id {
                return Err(cx_parse_failure());
            }
            let bond = state.builder.bond_mut(bond_id).ok_or_else(cx_parse_failure)?;
            if begin != atom_id {
                bond.set_endpoints(atom_id, begin);
            }
            bond.set_prop("_MolFileBondCfg", cfg);
            bond.set_direction(direction);
            if cfg == "2" && bond_can_have_direction(order) {
                if let Some(atom) = state.builder.atom_mut(atom_id) {
                    atom.set_chiral_tag(ChiralTag::Unspecified);
                }
                state.set_property("_needsDetectBondStereo", "1");
            }
            if matches!(cfg, "1" | "3") && bond_can_have_direction(order) {
                state.set_property("_needsDetectAtomStereo", "1");
            }
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

pub(super) fn parse_cx_variable_attachments(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_variable_attachments
    // RDKit✔️✔️: bool parse_variable_attachments(Iterator &first, Iterator last,
    // RDKit✔️✔️:                                 RDKit::RWMol &mol, unsigned int startAtomIdx) {
    // RDKit✔️✔️:   // these look like: CO*.C1=CC=NC=C1 |m:2:3.5.4|
    // RDKit✔️✔️:   // that corresponds to replacing the bond to atom 2 with bonds to atom 3, 5,
    // RDKit✔️✔️:   // or 4
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   if (first >= last || *first != 'm' || first + 1 >= last ||
    // RDKit✔️✔️:       *(first + 1) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int at1idx;
    // RDKit✔️✔️:     if (!read_int(first, last, at1idx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (VALID_ATIDX(at1idx) &&
    // RDKit✔️✔️:         mol.getAtomWithIdx(at1idx - startAtomIdx)->getDegree() != 1) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "position variation bond to atom with more than one bond"
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ':') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog) << "improperly formatted m: block" << std::endl;
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::vector<std::string> others;
    // RDKit✔️✔️:     while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:       unsigned int aidx;
    // RDKit✔️✔️:       if (!read_int(first, last, aidx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (VALID_ATIDX(aidx)) {
    // RDKit✔️✔️:         others.push_back(std::to_string(aidx - startAtomIdx + 1));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (first < last && *first == '.') {
    // RDKit✔️✔️:         ++first;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(at1idx)) {
    // RDKit✔️✔️:       std::string endPts = "(" + std::to_string(others.size());
    // RDKit✔️✔️:       for (auto idx : others) {
    // RDKit✔️✔️:         endPts += " " + idx;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       endPts += ")";
    // RDKit✔️✔️:
    // RDKit✔️✔️:       for (auto nbri : boost::make_iterator_range(
    // RDKit✔️✔️:                mol.getAtomBonds(mol.getAtomWithIdx(at1idx - startAtomIdx)))) {
    // RDKit✔️✔️:         auto bnd = mol[nbri];
    // RDKit✔️✔️:         bnd->setProp(common_properties::_MolFileBondEndPts, endPts);
    // RDKit✔️✔️:         bnd->setProp(common_properties::_MolFileBondAttach, std::string("ANY"));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_variable_attachments
    if !ext_text[*pos..].starts_with("m:") {
        return Err(cx_parse_failure());
    }
    *pos += 2;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        let atom_id = AtomId::new(atom_idx);
        if state.builder.atom_mut(atom_id).is_some() && atom_neighbors(state, atom_id).len() != 1 {
            return Err(cx_parse_failure());
        }
        expect_byte(ext_text, *pos, b':')?;
        *pos += 1;
        let mut endpoints = Vec::new();
        while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
            let endpoint = read_cx_usize(ext_text, pos)?;
            if state.builder.atom_mut(AtomId::new(endpoint)).is_some() {
                endpoints.push((endpoint + 1).to_string());
            }
            if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'.' {
                *pos += 1;
            }
        }
        if state.builder.atom_mut(atom_id).is_some() {
            let end_pts = if endpoints.is_empty() {
                "(0)".to_string()
            } else {
                format!("({} {})", endpoints.len(), endpoints.join(" "))
            };
            let cached_bonds: Vec<_> = state.builder.neighbor_bonds(atom_id).to_vec();
            for bond_id in cached_bonds {
                let bond = state.builder.bond_mut(bond_id).ok_or_else(cx_parse_failure)?;
                bond.set_prop("_MolFileBondEndPts", end_pts.clone());
                bond.set_prop("_MolFileBondAttach", "ANY");
            }
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

pub(super) fn set_cx_stereo_for_bond(
    state: &mut SmilesBuildState,
    bond_id: BondId,
    stereo: BondStereo,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::detail::setStereoForBond
    // RDKit✔️✔️: void setStereoForBond(ROMol &mol, Bond *bond, Bond::BondStereo stereo,
    // RDKit✔️✔️:                       bool useCXSmilesOrdering) {
    // RDKit✔️✔️:   // NOTE:  moved from parse_doublebond_stereo CXSmilesOps
    // RDKit✔️✔️:   // IF useCXSmilesOrdering is true, the cis/trans/unknown marker will be
    // RDKit✔️✔️:   // assigned relative to the lowest-numbered neighbor of each double bond atom.
    // RDKit✔️✔️:   // Otherwise it uses the lowest-numbered neighbor on the lower-numbered atom
    // RDKit✔️✔️:   // of the double bond and the highest-numbered neighbor on the higher-numbered
    // RDKit✔️✔️:   // atom
    // RDKit✔️✔️:   auto begAtom = bond->getBeginAtom();
    // RDKit✔️✔️:   auto endAtom = bond->getEndAtom();
    // RDKit✔️✔️:   if (begAtom->getIdx() > endAtom->getIdx()) {
    // RDKit✔️✔️:     std::swap(begAtom, endAtom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    // RDKit✔️✔️:     unsigned int begControl = mol.getNumAtoms();
    // RDKit✔️✔️:     for (auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit✔️✔️:       if (nbr == endAtom) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       begControl = std::min(nbr->getIdx(), begControl);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned int endControl = useCXSmilesOrdering ? mol.getNumAtoms() : 0;
    // RDKit✔️✔️:     for (auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit✔️✔️:       if (nbr == begAtom) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       endControl = useCXSmilesOrdering ? std::min(nbr->getIdx(), endControl)
    // RDKit✔️✔️:                                        : std::max(nbr->getIdx(), endControl);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (begAtom != bond->getBeginAtom()) {
    // RDKit✔️✔️:       std::swap(begControl, endControl);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bond->setStereoAtoms(begControl, endControl);
    // RDKit✔️✔️:     bond->setStereo(stereo);
    // RDKit✔️✔️:     mol.setProp("_needsDetectBondStereo", 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::detail::setStereoForBond
    let (begin, end) = state
        .builder
        .bond(bond_id)
        .map(|bond| (bond.begin(), bond.end()))
        .ok_or_else(cx_parse_failure)?;
    let (low, high) = if begin.index() > end.index() {
        (end, begin)
    } else {
        (begin, end)
    };
    let low_neighbors = atom_neighbors(state, low);
    let high_neighbors = atom_neighbors(state, high);
    if low_neighbors.len() <= 1 || high_neighbors.len() <= 1 {
        return Ok(());
    }
    let low_control = low_neighbors
        .iter()
        .copied()
        .filter(|neighbor| *neighbor != high)
        .min_by_key(|atom| atom.index())
        .ok_or_else(cx_parse_failure)?;
    let high_control = high_neighbors
        .iter()
        .copied()
        .filter(|neighbor| *neighbor != low)
        .min_by_key(|atom| atom.index())
        .ok_or_else(cx_parse_failure)?;
    let stereo_atoms = if low != begin {
        [high_control, low_control]
    } else {
        [low_control, high_control]
    };
    let bond = state.builder.bond_mut(bond_id).ok_or_else(cx_parse_failure)?;
    bond.set_stereo_atoms(Some(stereo_atoms));
    bond.set_stereo(stereo);
    state.set_property("_needsDetectBondStereo", "1");
    Ok(())
}

pub(super) fn parse_cx_doublebond_stereo(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    stereo: BondStereo,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_doublebond_stereo
    // RDKit✔️✔️: bool parse_doublebond_stereo(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                              unsigned int, unsigned int startBondIdx,
    // RDKit✔️✔️:                              Bond::BondStereo stereo) {
    // RDKit✔️✔️:   // these look like: C1CCCC/C=C/CCC1 |ctu:5|
    // RDKit✔️✔️:   // also c and t for cis or trans
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   while (first < last && *first != ':') {
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int bondIdx;
    // RDKit✔️✔️:     if (!read_int(first, last, bondIdx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_BNDIDX(bondIdx)) {
    // RDKit✔️✔️:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!bond) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "bond " << bondIdx
    // RDKit✔️✔️:             << " not found, cannot mark as stereo double bond." << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       bool useCXOrdering = true;
    // RDKit✔️✔️:       Chirality::detail::setStereoForBond(mol, bond, stereo, useCXOrdering);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_doublebond_stereo
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b':' {
        *pos += 1;
    }
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let bond_idx = read_cx_usize(ext_text, pos)?;
        let bond_id = cx_bond_with_smiles_idx(state, bond_idx).ok_or_else(cx_parse_failure)?;
        set_cx_stereo_for_bond(state, bond_id, stereo)?;
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

pub(super) fn process_cx_radical_section(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    radical_electrons: u8,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION processRadicalSection
    // RDKit✔️✔️: bool processRadicalSection(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                            unsigned int numRadicalElectrons,
    // RDKit✔️✔️:                            unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   unsigned int atIdx;
    // RDKit✔️✔️:   if (!read_int(first, last, atIdx)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (VALID_ATIDX(atIdx)) {
    // RDKit✔️✔️:     mol.getAtomWithIdx(atIdx - startAtomIdx)
    // RDKit✔️✔️:         ->setNumRadicalElectrons(numRadicalElectrons);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   while (first < last && *first == ',') {
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     if (first < last && (*first < '0' || *first > '9')) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!read_int(first, last, atIdx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(atIdx)) {
    // RDKit✔️✔️:       mol.getAtomWithIdx(atIdx - startAtomIdx)
    // RDKit✔️✔️:           ->setNumRadicalElectrons(numRadicalElectrons);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return first < last;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION processRadicalSection
    if *pos >= ext_text.len() {
        return Err(cx_parse_failure());
    }
    *pos += 1;
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    let atom_idx = read_cx_usize(ext_text, pos)?;
    if let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx)) {
        atom.set_radical_electrons(radical_electrons);
    }
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
        *pos += 1;
        if *pos < ext_text.len() && !ext_text.as_bytes()[*pos].is_ascii_digit() {
            return Ok(());
        }
        let atom_idx = read_cx_usize(ext_text, pos)?;
        if let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx)) {
            atom.set_radical_electrons(radical_electrons);
        }
    }
    if *pos < ext_text.len() {
        Ok(())
    } else {
        Err(cx_parse_failure())
    }
}

pub(super) fn parse_cx_radicals(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_radicals
    // RDKit✔️✔️: bool parse_radicals(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                     unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != '^') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   while (*first == '^') {
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     if (first >= last) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (*first < '1' || *first > '7') {
    // RDKit✔️✔️:       return false;  // these are the values that are allowed to be there
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     switch (*first) {
    // RDKit✔️✔️:       case '1':
    // RDKit✔️✔️:         if (!processRadicalSection(first, last, mol, 1, startAtomIdx)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case '2':
    // RDKit✔️✔️:       case '3':
    // RDKit✔️✔️:       case '4':
    // RDKit✔️✔️:         if (!processRadicalSection(first, last, mol, 2, startAtomIdx)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case '5':
    // RDKit✔️✔️:       case '6':
    // RDKit✔️✔️:       case '7':
    // RDKit✔️✔️:         if (!processRadicalSection(first, last, mol, 3, startAtomIdx)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "Radical specification " << *first << " ignored.";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_radicals
    expect_byte(ext_text, *pos, b'^')?;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'^' {
        *pos += 1;
        if *pos >= ext_text.len() {
            return Err(cx_parse_failure());
        }
        let radical_electrons = match ext_text.as_bytes()[*pos] {
            b'1' => 1,
            b'2' | b'3' | b'4' => 2,
            b'5' | b'6' | b'7' => 3,
            _ => return Err(cx_parse_failure()),
        };
        process_cx_radical_section(state, ext_text, pos, radical_electrons)?;
    }
    Ok(())
}

pub(super) fn cx_parse_failure() -> SmilesParseError {
    SmilesParseError::ParseError("failure parsing CXSMILES extensions".to_string())
}

pub(super) fn read_cx_usize(ext_text: &str, pos: &mut usize) -> Result<usize, SmilesParseError> {
    let start = *pos;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        *pos += 1;
    }
    if *pos == start {
        return Err(SmilesParseError::ParseError(
            "failure parsing CXSMILES extensions".to_string(),
        ));
    }
    ext_text[start..*pos]
        .parse::<usize>()
        .map_err(|_| SmilesParseError::ParseError("failure parsing CXSMILES extensions".to_string()))
}

pub(super) fn read_cx_text_to(ext_text: &str, pos: &mut usize, delimiters: &[u8]) -> Result<String, SmilesParseError> {
    let mut value = String::new();
    while *pos < ext_text.len() && !delimiters.contains(&ext_text.as_bytes()[*pos]) {
        if ext_text.as_bytes()[*pos] == b'&' && *pos + 2 < ext_text.len() && ext_text.as_bytes()[*pos + 1] == b'#' {
            *pos += 2;
            let numeric_start = *pos;
            while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
                *pos += 1;
            }
            if *pos >= ext_text.len() || ext_text.as_bytes()[*pos] != b';' {
                return Err(SmilesParseError::ParseError(
                    "failure parsing CXSMILES extensions: quoted block not terminated with ';'".to_string(),
                ));
            }
            if *pos > numeric_start {
                let code = ext_text[numeric_start..*pos]
                    .parse::<u32>()
                    .map_err(|_| SmilesParseError::ParseError("failure parsing CXSMILES extensions".to_string()))?;
                if let Some(ch) = char::from_u32(code) {
                    value.push(ch);
                }
            }
            *pos += 1;
        } else {
            value.push(ext_text.as_bytes()[*pos] as char);
            *pos += 1;
        }
    }
    Ok(value)
}

pub(super) fn expect_byte(ext_text: &str, pos: usize, expected: u8) -> Result<(), SmilesParseError> {
    if pos < ext_text.len() && ext_text.as_bytes()[pos] == expected {
        Ok(())
    } else {
        Err(SmilesParseError::ParseError(
            "failure parsing CXSMILES extensions".to_string(),
        ))
    }
}
