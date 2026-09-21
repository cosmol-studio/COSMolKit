use cosmolkit_cx::{
    CxCoordinateBondKind, CxDoubleBondStereoKind, CxRecord, CxStereoGroupKind, CxWedgeDirection,
    ParsedCxExtensions,
};
use cosmolkit_model::{
    AdjacencyList, AtomId, BondDirection, BondId, BondOrder, BondStereo, ChiralTag, Conformer2D,
    Conformer3D, CoordinateDimension, SGroupConnection, SGroupData, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock,
};

use crate::{CXSMILES_BOND_IDX_PROP, SmilesParseError, SmilesRecord};

fn cx_failure() -> SmilesParseError {
    SmilesParseError::Cx("failure parsing CXSMILES extensions".to_owned())
}

fn model_failure(error: impl std::fmt::Display) -> SmilesParseError {
    SmilesParseError::Model(error.to_string())
}

fn bond_with_smiles_index(topology: &TopologyBlock, index: usize) -> Option<BondId> {
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
    topology
        .bonds
        .iter()
        .find(|bond| {
            bond.prop(CXSMILES_BOND_IDX_PROP)
                .and_then(|value| value.parse::<usize>().ok())
                == Some(index)
        })
        .map(|bond| bond.id())
}

fn atom_neighbors(topology: &TopologyBlock, atom: AtomId) -> &[cosmolkit_model::NeighborRef] {
    topology.adjacency.neighbors_of(atom.index())
}

fn can_have_direction(order: BondOrder) -> bool {
    // BEGIN RDKIT CPP FUNCTION canHaveDirection
    // RDKit✔️✔️: inline bool canHaveDirection(const Bond &bond) {
    // RDKit✔️✔️:   auto bondType = bond.getBondType();
    // RDKit✔️✔️:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION canHaveDirection
    matches!(order, BondOrder::Single | BondOrder::Aromatic)
}

fn set_double_bond_stereo(
    record: &mut SmilesRecord,
    bond_id: BondId,
    stereo: BondStereo,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::detail::setStereoForBond
    // RDKit✔️✔️: auto begAtom = bond->getBeginAtom();
    // RDKit✔️✔️: auto endAtom = bond->getEndAtom();
    // RDKit✔️✔️: if (begAtom->getIdx() > endAtom->getIdx()) {
    // RDKit✔️✔️:   std::swap(begAtom, endAtom);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    // RDKit✔️✔️:   unsigned int begControl = mol.getNumAtoms();
    // RDKit✔️✔️:   for (auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit✔️✔️:     if (nbr == endAtom) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     begControl = std::min(nbr->getIdx(), begControl);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int endControl = useCXSmilesOrdering ? mol.getNumAtoms() : 0;
    // RDKit✔️✔️:   for (auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit✔️✔️:     if (nbr == begAtom) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     endControl = useCXSmilesOrdering ? std::min(nbr->getIdx(), endControl)
    // RDKit✔️✔️:                                      : std::max(nbr->getIdx(), endControl);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (begAtom != bond->getBeginAtom()) {
    // RDKit✔️✔️:     std::swap(begControl, endControl);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bond->setStereoAtoms(begControl, endControl);
    // RDKit✔️✔️:   bond->setStereo(stereo);
    // RDKit✔️✔️:   mol.setProp("_needsDetectBondStereo", 1);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::detail::setStereoForBond
    let (begin, end) = record
        .topology
        .bonds
        .get(bond_id.index())
        .map(|bond| (bond.begin(), bond.end()))
        .ok_or_else(cx_failure)?;
    let (low, high) = if begin.index() > end.index() {
        (end, begin)
    } else {
        (begin, end)
    };
    let low_neighbors = atom_neighbors(&record.topology, low);
    let high_neighbors = atom_neighbors(&record.topology, high);
    if low_neighbors.len() <= 1 || high_neighbors.len() <= 1 {
        return Ok(());
    }
    let low_control = low_neighbors
        .iter()
        .map(|neighbor| AtomId::new(neighbor.atom_index))
        .filter(|neighbor| *neighbor != high)
        .min_by_key(|atom| atom.index())
        .ok_or_else(cx_failure)?;
    let high_control = high_neighbors
        .iter()
        .map(|neighbor| AtomId::new(neighbor.atom_index))
        .filter(|neighbor| *neighbor != low)
        .min_by_key(|atom| atom.index())
        .ok_or_else(cx_failure)?;
    let stereo_atoms = if low != begin {
        [high_control, low_control]
    } else {
        [low_control, high_control]
    };
    let bond = record
        .topology
        .bonds
        .get_mut(bond_id.index())
        .ok_or_else(cx_failure)?;
    bond.set_stereo_atoms(Some(stereo_atoms));
    bond.set_stereo(stereo).map_err(model_failure)?;
    record
        .properties
        .set_prop("_needsDetectBondStereo", "1")
        .map_err(model_failure)?;
    Ok(())
}

fn sgroup_kind(type_code: &str) -> Option<SubstanceGroupKind> {
    // BEGIN RDKIT CPP VALUE sgroupTypemap
    // RDKit✔️🔝: const std::map<std::string, std::string> sgroupTypemap = {
    // RDKit✔️🔝:     {"n", "SRU"},   {"mon", "MON"}, {"mer", "MER"}, {"co", "COP"},
    // RDKit✔️🔝:     {"xl", "CRO"},  {"mod", "MOD"}, {"mix", "MIX"}, {"f", "FOR"},
    // RDKit✔️🔝:     {"any", "ANY"}, {"gen", "GEN"}, {"c", "COM"},   {"grf", "GRA"},
    // RDKit✔️🔝:     {"alt", "COP"}, {"ran", "COP"}, {"blk", "COP"}};
    // END RDKIT CPP VALUE sgroupTypemap
    // The fixed string match preserves the 15-entry source mapping without
    // RDKit's O(log n) tree lookup or runtime map storage.
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
        "gen" => Some(SubstanceGroupKind::Generic("GEN".to_owned())),
        "c" => Some(SubstanceGroupKind::Generic("COM".to_owned())),
        "grf" => Some(SubstanceGroupKind::Graft),
        "alt" | "ran" | "blk" => Some(SubstanceGroupKind::Copolymer),
        _ => None,
    }
}

fn unsupported_query_label(label: &str) -> bool {
    // BEGIN RDKIT CPP FUNCTION processCXSmilesLabels (query labels)
    // RDKit❌❌: if (symb == "star_e") { addquery(makeAtomNullQuery(), symb, mol, atom->getIdx()); }
    // RDKit❌❌: else if (symb == "Q_e") { addquery(makeQAtomQuery(), symb, mol, atom->getIdx()); }
    // RDKit❌❌: else if (symb == "QH_p") { addquery(makeQHAtomQuery(), symb, mol, atom->getIdx()); }
    // RDKit❌❌: else if (symb == "AH_p") { addquery(makeAHAtomQuery(), symb, mol, atom->getIdx()); }
    // RDKit❌❌: else if (symb == "X_p") { addquery(makeXAtomQuery(), symb, mol, atom->getIdx()); }
    // RDKit❌❌: else if (symb == "XH_p") { addquery(makeXHAtomQuery(), symb, mol, atom->getIdx()); }
    // RDKit❌❌: else if (symb == "M_p") { addquery(makeMAtomQuery(), symb, mol, atom->getIdx()); }
    // RDKit❌❌: else if (symb == "MH_p") { addquery(makeMHAtomQuery(), symb, mol, atom->getIdx()); }
    // END RDKIT CPP FUNCTION processCXSmilesLabels (query labels)
    // Concrete Atom values do not contain a QueryGraph. The caller returns a
    // structured unsupported error instead of preserving a misleading label.
    matches!(
        label,
        "star_e" | "Q_e" | "QH_p" | "AH_p" | "X_p" | "XH_p" | "M_p" | "MH_p"
    )
}

fn normalize_source_coordinate_dimension(record: &mut SmilesRecord) {
    record.coordinates.source_coordinate_dim = if record.coordinates.conformers_3d.is_empty() {
        (!record.coordinates.conformers_2d.is_empty()).then_some(CoordinateDimension::TwoD)
    } else {
        Some(CoordinateDimension::ThreeD)
    };
}

fn polymer_crossing_bond(
    topology: &TopologyBlock,
    source_index: usize,
) -> Result<Option<BondId>, SmilesParseError> {
    // RDKit validates these source values with VALID_ATIDX before later using
    // them as bond indices. Preserve the source filter, then turn the unsafe
    // downstream bond access into a structured detached-model failure.
    if source_index >= topology.atoms.len() {
        return Ok(None);
    }
    if source_index >= topology.bonds.len() {
        return Err(cx_failure());
    }
    Ok(Some(BondId::new(source_index)))
}

fn infer_unmarked_polymer_crossings(
    topology: &TopologyBlock,
    atoms: &[AtomId],
) -> Result<(Vec<BondId>, Vec<BondId>), SmilesParseError> {
    if atoms.is_empty() {
        return Err(cx_failure());
    }
    let mut head = Vec::new();
    let mut tail = Vec::new();
    let first = atoms[0];
    for neighbor in atom_neighbors(topology, first) {
        let neighbor_atom = AtomId::new(neighbor.atom_index);
        if atoms.contains(&neighbor_atom) {
            continue;
        }
        if atoms.len() > 1 || head.is_empty() {
            head.push(neighbor.bond);
        } else if tail.is_empty() {
            tail.push(neighbor.bond);
        }
    }
    if atoms.len() > 1 {
        let last = *atoms.last().expect("polymer atoms are non-empty");
        for neighbor in atom_neighbors(topology, last) {
            if !atoms.contains(&AtomId::new(neighbor.atom_index)) {
                tail.push(neighbor.bond);
            }
        }
    }
    Ok((head, tail))
}

fn normalize_polymer_connection(value: &str) -> (SGroupConnection, String, bool) {
    let flipped = value.contains(",f");
    let without_flip = value.replace(",f", "");
    match without_flip.as_str() {
        "hh" | "HH" => (SGroupConnection::HeadToHead, "HH".to_owned(), flipped),
        "ht" | "HT" => (SGroupConnection::HeadToTail, "HT".to_owned(), flipped),
        "eu" | "EU" | "" => (SGroupConnection::Either, "EU".to_owned(), flipped),
        _ => (SGroupConnection::Either, "EU".to_owned(), flipped),
    }
}

fn finalize_polymer_sgroup(
    topology: &TopologyBlock,
    group: &mut SubstanceGroup,
    source_connect: &str,
    source_head: &[usize],
    source_tail: &[usize],
) -> Result<bool, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION finalizePolymerSGroup
    // RDKit✔️✔️: if (connect.find(",f") != std::string::npos) { isFlipped = true; }
    // RDKit✔️✔️: if (connect == "hh") connect = "HH";
    // RDKit✔️✔️: else if (connect == "ht") connect = "HT";
    // RDKit✔️✔️: else if (connect == "eu") connect = "EU";
    // RDKit✔️✔️: else connect = "EU";
    // RDKit✔️✔️: if (headCrossings.empty() && tailCrossings.empty()) {
    // RDKit✔️✔️:   setupUnmarkedPolymerSGroup(mol, sgroup, headCrossings, tailCrossings);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: for (auto &bondIdx : headCrossings) sgroup.addBondWithIdx(bondIdx);
    // RDKit✔️✔️: sgroup.setProp("XBHEAD", headCrossings);
    // RDKit✔️✔️: for (auto &bondIdx : tailCrossings) sgroup.addBondWithIdx(bondIdx);
    // RDKit✔️✔️: for (unsigned int i = 0;
    // RDKit✔️✔️:      i < std::min(headCrossings.size(), tailCrossings.size()); ++i) {
    // RDKit✔️✔️:   unsigned tailIdx = isFlipped
    // RDKit✔️✔️:       ? tailCrossings[tailCrossings.size() - i - 1] : tailCrossings[i];
    // RDKit✔️✔️:   xbcorr.push_back(headCrossings[i]); xbcorr.push_back(tailIdx);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION finalizePolymerSGroup
    let (connection, connect, flipped) = normalize_polymer_connection(source_connect);
    group.set_connection(connection);
    group.set_prop("CONNECT", connect);

    let mut head = Vec::new();
    let mut tail = Vec::new();
    let mut valid = true;
    for index in source_head {
        match polymer_crossing_bond(topology, *index)? {
            Some(bond) => head.push(bond),
            None => valid = false,
        }
    }
    for index in source_tail {
        match polymer_crossing_bond(topology, *index)? {
            Some(bond) => tail.push(bond),
            None => valid = false,
        }
    }
    if !valid {
        return Ok(false);
    }
    if head.is_empty() && tail.is_empty() {
        (head, tail) = infer_unmarked_polymer_crossings(topology, group.atoms())?;
    }
    if head.is_empty() && tail.is_empty() {
        return Ok(true);
    }
    for bond in head.iter().chain(&tail) {
        group.push_bond(*bond);
    }
    for bond in &head {
        group.push_head_crossing_bond(*bond);
    }
    let mut xbcorr = Vec::with_capacity(head.len().min(tail.len()) * 2);
    for index in 0..head.len().min(tail.len()) {
        xbcorr.push(head[index]);
        xbcorr.push(if flipped {
            tail[tail.len() - index - 1]
        } else {
            tail[index]
        });
    }
    for bond in xbcorr {
        group.push_crossing_bond_correspondence(bond);
    }
    Ok(true)
}

/// Apply representation-independent CX records to detached concrete model
/// blocks. Query-only records are rejected before any destination mutation.
pub(crate) fn apply_cx_to_smiles_record(
    record: &mut SmilesRecord,
    parsed: &ParsedCxExtensions,
) -> Result<(), SmilesParseError> {
    // RDKit✔️❌: parseCXExtensions mutates the destination while parsing.
    // COSMolKit stages one detached clone so non-strict recovery cannot expose
    // a partially installed CX record. This preserves behavior on success but
    // adds one O(molecule-size) clone to the source's record-linear work.
    let mut staged = record.clone();
    apply_cx_to_smiles_record_in_place(&mut staged, parsed)?;
    *record = staged;
    Ok(())
}

fn apply_cx_to_smiles_record_in_place(
    record: &mut SmilesRecord,
    parsed: &ParsedCxExtensions,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parser::parse_it dispatch
    // RDKit✔️✔️: } else if (*first == 'C') {
    // RDKit✔️✔️:   if (!parse_coordinate_bonds(first, last, mol, Bond::DATIVE, startAtomIdx,
    // RDKit✔️✔️:                               startBondIdx)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (*first == 'H') {
    // RDKit✔️✔️:   if (!parse_coordinate_bonds(first, last, mol, Bond::HYDROGEN,
    // RDKit✔️✔️:                               startAtomIdx, startBondIdx)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (*first == 'Z') {
    // RDKit✔️✔️:   if (!parse_zero_bonds(first, last, mol, startAtomIdx, startBondIdx)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (*first == '^') {
    // RDKit✔️✔️:   if (!parse_radicals(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (*first == 'a' || *first == 'o' ||
    // RDKit✔️✔️:            (*first == '&' && first + 1 < last && first[1] != '#')) {
    // RDKit✔️✔️:   if (!parse_enhanced_stereo(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (*first == 'L' && first + 1 < last && first[1] == 'N') {
    // RDKit✔️✔️:   if (!parse_linknodes(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
    // RDKit✔️✔️:            first[2] == 'D') {
    // RDKit✔️✔️:   if (!parse_data_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
    // RDKit✔️✔️:            first[2] == 'H') {
    // RDKit✔️✔️:   if (!parse_sgroup_hierarchy(first, last, mol)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (*first == 'S' && first + 1 < last && first[1] == 'g') {
    // RDKit✔️✔️:   if (!parse_polymer_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (*first == 'm') {
    // RDKit✔️✔️:   if (!parse_variable_attachments(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (*first == 'w') {
    // RDKit✔️✔️:   if (!parse_wedged_bonds(first, last, mol, startAtomIdx, startBondIdx)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION parser::parse_it dispatch
    for item in parsed.records() {
        let reason = match item {
            CxRecord::Unsaturation(_) => Some(
                "CX atom-query records require a QueryGraph and are not representable in a concrete Molecule",
            ),
            CxRecord::RingBonds(_) => Some(
                "CX ring-bond query records require a QueryGraph and are not representable in a concrete Molecule",
            ),
            CxRecord::Substitution(_) => Some(
                "CX substitution query records require a QueryGraph and are not representable in a concrete Molecule",
            ),
            CxRecord::AtomLabels(values)
                if values
                    .iter()
                    .flatten()
                    .any(|label| unsupported_query_label(label)) =>
            {
                Some(
                    "CX special query-atom labels require a QueryGraph and are not representable in a concrete Molecule",
                )
            }
            _ => None,
        };
        if let Some(reason) = reason {
            return Err(SmilesParseError::UnsupportedCx(reason));
        }
    }

    let atom_count = record.topology.atoms.len();
    let mut sgroup_index = 0;
    for item in parsed.records() {
        match item {
            CxRecord::Coordinates(coordinates) => {
                // RDKit✔️✔️: auto *conf = new Conformer(mol.getNumAtoms());
                // RDKit✔️✔️: mol.addConformer(conf);
                // RDKit✔️✔️: conf->setId(confIdx);
                // RDKit✔️✔️: conf->setAtomPos(atIdx - startAtomIdx, pt);
                // RDKit✔️✔️: if (is3D && hasNonZeroZCoords(*conf)) conf->set3D(true);
                // RDKit✔️✔️: else conf->set3D(false);
                let mut points = vec![[0.0; 3]; atom_count];
                for (index, value) in coordinates.values.iter().enumerate() {
                    if index >= atom_count {
                        break;
                    }
                    if let Some(value) = value {
                        points[index] = *value;
                    }
                }
                if coordinates.is_3d {
                    record.coordinates.conformers_3d.push(Conformer3D::new(
                        coordinates.conformer,
                        points,
                        true,
                    ));
                } else {
                    record.coordinates.conformers_2d.push(Conformer2D::new(
                        coordinates.conformer,
                        points
                            .into_iter()
                            .map(|point| [point[0], point[1]])
                            .collect(),
                    ));
                }
            }
            CxRecord::AtomLabels(values) => {
                // RDKit✔️✔️: mol.getAtomWithIdx(atIdx - startAtomIdx)
                // RDKit✔️✔️:     ->setProp(RDKit::common_properties::atomLabel, tkn);
                for (index, value) in values.iter().enumerate() {
                    if let Some(value) = value
                        && let Some(atom) = record.topology.atoms.get_mut(index)
                    {
                        atom.set_prop("atomLabel", value).map_err(model_failure)?;
                        if matches!(value.as_str(), "Pol_p" | "Mod_p") {
                            // RDKit✔️✔️: atom->setProp(common_properties::dummyLabel,
                            // RDKit✔️✔️:               symb.substr(0, symb.size() - 2));
                            // RDKit✔️✔️: atom->clearProp(common_properties::atomLabel);
                            atom.set_prop("dummyLabel", &value[..value.len() - 2])
                                .map_err(model_failure)?;
                            atom.clear_prop("atomLabel");
                        }
                    }
                }
            }
            CxRecord::AtomValues(values) => {
                // RDKit✔️✔️: mol.getAtomWithIdx(atIdx)->setProp(
                // RDKit✔️✔️:     RDKit::common_properties::molFileValue, tkn);
                for (index, value) in values.iter().enumerate() {
                    if let Some(value) = value
                        && let Some(atom) = record.topology.atoms.get_mut(index)
                    {
                        atom.set_prop("molFileValue", value)
                            .map_err(model_failure)?;
                    }
                }
            }
            CxRecord::AtomProperties(properties) => {
                // RDKit✔️✔️: mol.getAtomWithIdx(atIdx - startAtomIdx)->setProp(pname, pval);
                for property in properties {
                    if let Some(atom) = record.topology.atoms.get_mut(property.atom) {
                        atom.set_prop(property.name.clone(), property.value.clone())
                            .map_err(model_failure)?;
                    }
                }
            }
            CxRecord::CoordinateBonds(annotation) => {
                // RDKit✔️✔️: bnd->setBondType(typ);
                // RDKit✔️✔️: if (bnd->getBeginAtomIdx() != aidx - startAtomIdx) {
                // RDKit✔️✔️:   unsigned int tmp = bnd->getBeginAtomIdx();
                // RDKit✔️✔️:   bnd->setBeginAtomIdx(aidx - startAtomIdx);
                // RDKit✔️✔️:   bnd->setEndAtomIdx(tmp);
                // RDKit✔️✔️: }
                let order = match annotation.kind {
                    CxCoordinateBondKind::Dative => BondOrder::Dative,
                    CxCoordinateBondKind::Hydrogen => BondOrder::Hydrogen,
                };
                for reference in &annotation.bonds {
                    let atom = AtomId::new(reference.atom);
                    if atom.index() >= atom_count || reference.bond >= record.topology.bonds.len() {
                        continue;
                    }
                    let bond_id = bond_with_smiles_index(&record.topology, reference.bond)
                        .ok_or_else(cx_failure)?;
                    let (begin, end) = record
                        .topology
                        .bonds
                        .get(bond_id.index())
                        .map(|bond| (bond.begin(), bond.end()))
                        .ok_or_else(cx_failure)?;
                    if begin != atom && end != atom {
                        return Err(cx_failure());
                    }
                    let bond = record
                        .topology
                        .bonds
                        .get_mut(bond_id.index())
                        .ok_or_else(cx_failure)?;
                    bond.set_order(order);
                    if begin != atom {
                        bond.set_endpoints(atom, begin);
                    }
                }
                record.topology.adjacency = AdjacencyList::from_topology(
                    record.topology.atoms.len(),
                    &record.topology.bonds,
                );
            }
            CxRecord::ZeroBonds(indices) => {
                // RDKit✔️✔️: bond->setBondType(Bond::ZERO);
                for index in indices {
                    if *index >= record.topology.bonds.len() {
                        continue;
                    }
                    let bond_id =
                        bond_with_smiles_index(&record.topology, *index).ok_or_else(cx_failure)?;
                    record.topology.bonds[bond_id.index()].set_order(BondOrder::Zero);
                }
            }
            CxRecord::EnhancedStereo(stereo) => {
                // RDKit✔️✔️: if (iter != sgTracker.end()) {
                // RDKit✔️✔️:   auto gAtoms = mol_stereo_groups[index].getAtoms();
                // RDKit✔️✔️:   gAtoms.insert(gAtoms.end(), atoms.begin(), atoms.end());
                // RDKit✔️✔️: } else {
                // RDKit✔️✔️:   mol_stereo_groups.emplace_back(group_type, std::move(atoms),
                // RDKit✔️✔️:                                  std::move(bonds), group_id);
                // RDKit✔️✔️: }
                let kind = match stereo.kind {
                    CxStereoGroupKind::Absolute => StereoGroupKind::Absolute,
                    CxStereoGroupKind::Or => StereoGroupKind::Or,
                    CxStereoGroupKind::And => StereoGroupKind::And,
                };
                let atoms = stereo
                    .atoms
                    .iter()
                    .filter(|index| **index < atom_count)
                    .map(|index| AtomId::new(*index))
                    .collect::<Vec<_>>();
                if atoms.is_empty() {
                    continue;
                }
                if let Some(group) = record
                    .topology
                    .stereo_groups
                    .iter_mut()
                    .find(|group| group.kind() == kind && group.id() == Some(stereo.group_id))
                {
                    for atom in atoms {
                        group.push_atom(atom);
                    }
                } else {
                    record
                        .topology
                        .stereo_groups
                        .push(StereoGroup::new(kind, atoms, Vec::new()).with_id(stereo.group_id));
                }
            }
            CxRecord::WedgedBonds(wedges) => {
                // RDKit✔️✔️: bond->setProp(common_properties::_MolFileBondCfg, cfg);
                // RDKit✔️✔️: bond->setBondDir(state);
                // RDKit✔️✔️: if (cfg == 2 && canHaveDirection(*bond)) {
                // RDKit✔️✔️:   bond->getBeginAtom()->setChiralTag(Atom::CHI_UNSPECIFIED);
                // RDKit✔️✔️:   mol.setProp(detail::_needsDetectBondStereo, 1);
                // RDKit✔️✔️: }
                // RDKit✔️✔️: if ((cfg == 1 || cfg == 3) && canHaveDirection(*bond)) {
                // RDKit✔️✔️:   mol.setProp(detail::_needsDetectAtomStereo, 1);
                // RDKit✔️✔️: }
                for wedge in wedges {
                    if wedge.atom >= atom_count || wedge.bond >= record.topology.bonds.len() {
                        continue;
                    }
                    let bond_id = bond_with_smiles_index(&record.topology, wedge.bond)
                        .ok_or_else(cx_failure)?;
                    let (begin, end, order, has_cfg) = record
                        .topology
                        .bonds
                        .get(bond_id.index())
                        .map(|bond| {
                            (
                                bond.begin(),
                                bond.end(),
                                bond.order(),
                                bond.prop("_MolFileBondCfg").is_some(),
                            )
                        })
                        .ok_or_else(cx_failure)?;
                    if has_cfg {
                        return Err(cx_failure());
                    }
                    let atom = AtomId::new(wedge.atom);
                    if begin != atom && end != atom {
                        return Err(cx_failure());
                    }
                    let (cfg, direction) = match wedge.direction {
                        CxWedgeDirection::Unknown => ("2", BondDirection::Unknown),
                        CxWedgeDirection::BeginWedge => ("1", BondDirection::BeginWedge),
                        CxWedgeDirection::BeginDash => ("3", BondDirection::BeginDash),
                    };
                    let bond = &mut record.topology.bonds[bond_id.index()];
                    if begin != atom {
                        bond.set_endpoints(atom, begin);
                    }
                    bond.set_prop("_MolFileBondCfg", cfg)
                        .map_err(model_failure)?;
                    bond.set_direction(direction);
                    if cfg == "2" && can_have_direction(order) {
                        record.topology.atoms[atom.index()].set_chiral_tag(ChiralTag::Unspecified);
                        record
                            .properties
                            .set_prop("_needsDetectBondStereo", "1")
                            .map_err(model_failure)?;
                    }
                    if matches!(cfg, "1" | "3") && can_have_direction(order) {
                        record
                            .properties
                            .set_prop("_needsDetectAtomStereo", "1")
                            .map_err(model_failure)?;
                    }
                }
                record.topology.adjacency = AdjacencyList::from_topology(
                    record.topology.atoms.len(),
                    &record.topology.bonds,
                );
            }
            CxRecord::DoubleBondStereo(stereo) => {
                let value = match stereo.stereo {
                    CxDoubleBondStereoKind::Any => BondStereo::Any,
                    CxDoubleBondStereoKind::Cis => BondStereo::Cis,
                    CxDoubleBondStereoKind::Trans => BondStereo::Trans,
                };
                for index in &stereo.bonds {
                    if *index >= record.topology.bonds.len() {
                        continue;
                    }
                    let bond_id =
                        bond_with_smiles_index(&record.topology, *index).ok_or_else(cx_failure)?;
                    set_double_bond_stereo(record, bond_id, value)?;
                }
            }
            CxRecord::Radicals(radicals) => {
                // RDKit✔️✔️: mol.getAtomWithIdx(atIdx - startAtomIdx)
                // RDKit✔️✔️:     ->setNumRadicalElectrons(numRadicalElectrons);
                for radical in radicals {
                    if let Some(atom) = record.topology.atoms.get_mut(radical.atom) {
                        atom.set_radical_electrons(radical.electrons);
                    }
                }
            }
            CxRecord::LinkNodes(link_nodes) => {
                // RDKit✔️✔️: accum += (boost::format("%d %d 2 %d %d %d %d") %
                // RDKit✔️✔️:     startReps % endReps % (atidx - startAtomIdx + 1) %
                // RDKit✔️✔️:     (idx1 - startAtomIdx + 1) % (atidx - startAtomIdx + 1) %
                // RDKit✔️✔️:     (idx2 - startAtomIdx + 1)).str();
                // RDKit✔️✔️: mol.setProp(common_properties::molFileLinkNodes, accum);
                let mut lowered = Vec::new();
                for link in link_nodes {
                    let atom = AtomId::new(link.atom);
                    if atom.index() >= atom_count {
                        continue;
                    }
                    let outer_atoms = if let Some(outer_atoms) = link.outer_atoms {
                        outer_atoms
                    } else {
                        let neighbors = atom_neighbors(&record.topology, atom);
                        if neighbors.len() != 2 {
                            return Err(cx_failure());
                        }
                        [neighbors[0].atom_index, neighbors[1].atom_index]
                    };
                    lowered.push(format!(
                        "{} {} 2 {} {} {} {}",
                        link.start_repetitions,
                        link.end_repetitions,
                        link.atom + 1,
                        outer_atoms[0] + 1,
                        link.atom + 1,
                        outer_atoms[1] + 1
                    ));
                }
                if !lowered.is_empty() {
                    record
                        .properties
                        .set_prop("_MolFileLinkNodes", lowered.join("|"))
                        .map_err(model_failure)?;
                }
            }
            CxRecord::DataSGroup(data) => {
                // RDKit✔️✔️: SubstanceGroup sgroup(&mol, std::string("DAT"));
                // RDKit✔️✔️: sgroup.setProp(cxsmilesindex, nSGroups);
                // RDKit✔️✔️: sgroup.addAtomWithIdx(idx - startAtomIdx);
                // RDKit✔️✔️: sgroup.setProp("FIELDDISP",
                // RDKit✔️✔️:     "    0.0000    0.0000    DR    ALL  0       0");
                // RDKit✔️✔️: addSubstanceGroup(mol, sgroup);
                let atoms = data
                    .atoms
                    .iter()
                    .filter(|index| **index < atom_count)
                    .map(|index| AtomId::new(*index))
                    .collect::<Vec<_>>();
                if !atoms.is_empty() {
                    let group_id = SubstanceGroupId::new(record.topology.substance_groups.len());
                    let typed_data = SGroupData {
                        field_name: (!data.field_name.is_empty()).then(|| data.field_name.clone()),
                        field_info: (!data.field_info.is_empty()).then(|| data.field_info.clone()),
                        field_display: Some(
                            "    0.0000    0.0000    DR    ALL  0       0".to_owned(),
                        ),
                        query_op: (!data.query_op.is_empty()).then(|| data.query_op.clone()),
                        values: (!data.data.is_empty())
                            .then(|| vec![data.data.clone()])
                            .unwrap_or_default(),
                        ..SGroupData::default()
                    };
                    let mut group = SubstanceGroup::new(group_id, SubstanceGroupKind::Data)
                        .with_atoms(atoms)
                        .with_data(typed_data);
                    group.set_prop("_cxsmilesindex", sgroup_index.to_string());
                    group.set_prop(
                        "index",
                        (record.topology.substance_groups.len() + 1).to_string(),
                    );
                    group.set_prop("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0");
                    if !data.field_name.is_empty() {
                        group.set_prop("FIELDNAME", data.field_name.clone());
                    }
                    if !data.data.is_empty() {
                        group.set_prop("DATAFIELDS", data.data.clone());
                        group.push_data_field(data.data.clone());
                    }
                    if !data.query_op.is_empty() {
                        group.set_prop("QUERYOP", data.query_op.clone());
                    }
                    if !data.field_info.is_empty() {
                        group.set_prop("FIELDINFO", data.field_info.clone());
                    }
                    if !data.field_tag.is_empty() {
                        group.set_prop("FIELDTAG", data.field_tag.clone());
                    }
                    if let Some(coordinates) = &data.coordinates {
                        group.set_prop("COORDS", coordinates.clone());
                    }
                    record.topology.substance_groups.push(group);
                }
                sgroup_index += 1;
            }
            CxRecord::SGroupHierarchy(relationships) => {
                // RDKit✔️✔️: auto psg = find_matching_sgroup(sgs, parentId);
                // RDKit✔️✔️: psg->getPropIfPresent("index", parentId);
                // RDKit✔️✔️: if (childId >= sgs.size()) {
                // RDKit✔️✔️:   throw SmilesParseException(
                // RDKit✔️✔️:       "child id references non-existent SGroup");
                // RDKit✔️✔️: }
                // RDKit✔️✔️: csg->setProp("PARENT", parentId);
                let cx_indices = record
                    .topology
                    .substance_groups
                    .iter()
                    .enumerate()
                    .filter_map(|(position, group)| {
                        let cx_index =
                            group.props().get("_cxsmilesindex")?.parse::<usize>().ok()?;
                        let index = group
                            .props()
                            .get("index")
                            .and_then(|value| value.parse::<usize>().ok())
                            .unwrap_or(position);
                        Some((cx_index, position, index, group.id()))
                    })
                    .collect::<Vec<_>>();
                for relationship in relationships {
                    let Some((_, _, parent_index, parent_id)) = cx_indices
                        .iter()
                        .find(|(cx_index, _, _, _)| *cx_index == relationship.parent)
                    else {
                        continue;
                    };
                    for child in &relationship.children {
                        if *child >= record.topology.substance_groups.len() {
                            return Err(SmilesParseError::Cx(
                                "child id references non-existent SGroup".to_owned(),
                            ));
                        }
                        if let Some((_, child_position, _, _)) = cx_indices
                            .iter()
                            .find(|(cx_index, _, _, _)| cx_index == child)
                        {
                            let child_group =
                                &mut record.topology.substance_groups[*child_position];
                            child_group.set_parent(*parent_id);
                            child_group.set_prop("PARENT", parent_index.to_string());
                        }
                    }
                }
            }
            CxRecord::PolymerSGroup(polymer) => {
                let kind = sgroup_kind(&polymer.type_code).ok_or_else(cx_failure)?;
                let atoms = polymer
                    .atoms
                    .iter()
                    .filter(|index| **index < atom_count)
                    .map(|index| AtomId::new(*index))
                    .collect::<Vec<_>>();
                if !atoms.is_empty() {
                    let mut group = SubstanceGroup::new(
                        SubstanceGroupId::new(record.topology.substance_groups.len()),
                        kind,
                    )
                    .with_atoms(atoms);
                    group.set_prop("_cxsmilesindex", sgroup_index.to_string());
                    group.set_prop(
                        "index",
                        (record.topology.substance_groups.len() + 1).to_string(),
                    );
                    match polymer.type_code.as_str() {
                        "alt" => {
                            group.set_prop("SUBTYPE", "ALT");
                            group.set_subtype("ALT");
                        }
                        "ran" => {
                            group.set_prop("SUBTYPE", "RAN");
                            group.set_subtype("RAN");
                        }
                        "blk" => {
                            group.set_prop("SUBTYPE", "BLO");
                            group.set_subtype("BLO");
                        }
                        _ => {}
                    }
                    if !polymer.label.is_empty() {
                        group.set_prop("LABEL", polymer.label.clone());
                        group.set_label(polymer.label.clone());
                    }
                    let keep_group = finalize_polymer_sgroup(
                        &record.topology,
                        &mut group,
                        &polymer.connect,
                        &polymer.head_crossings,
                        &polymer.tail_crossings,
                    )?;
                    if keep_group {
                        record.topology.substance_groups.push(group);
                    }
                }
                sgroup_index += 1;
            }
            CxRecord::VariableAttachments(attachments) => {
                // RDKit✔️✔️: bnd->setProp(common_properties::_MolFileBondEndPts, endPts);
                // RDKit✔️✔️: bnd->setProp(common_properties::_MolFileBondAttach,
                // RDKit✔️✔️:              std::string("ANY"));
                for attachment in attachments {
                    let atom = AtomId::new(attachment.atom);
                    if atom.index() >= atom_count {
                        continue;
                    }
                    if atom_neighbors(&record.topology, atom).len() != 1 {
                        return Err(cx_failure());
                    }
                    let endpoints = attachment
                        .endpoints
                        .iter()
                        .filter(|index| **index < atom_count)
                        .map(|index| (index + 1).to_string())
                        .collect::<Vec<_>>();
                    let value = if endpoints.is_empty() {
                        "(0)".to_owned()
                    } else {
                        format!("({} {})", endpoints.len(), endpoints.join(" "))
                    };
                    let bond_ids = record
                        .topology
                        .bonds
                        .iter()
                        .filter(|bond| bond.begin() == atom || bond.end() == atom)
                        .map(|bond| bond.id())
                        .collect::<Vec<_>>();
                    for bond_id in bond_ids {
                        let bond = &mut record.topology.bonds[bond_id.index()];
                        bond.set_prop("_MolFileBondEndPts", value.clone())
                            .map_err(model_failure)?;
                        bond.set_prop("_MolFileBondAttach", "ANY")
                            .map_err(model_failure)?;
                    }
                }
            }
            CxRecord::Unsaturation(_) | CxRecord::RingBonds(_) | CxRecord::Substitution(_) => {
                unreachable!("query records were preflighted")
            }
            CxRecord::Unknown(_) => {}
        }
    }
    if record.properties.prop("_needsDetectAtomStereo").is_some() {
        // BEGIN RDKIT CPP FUNCTION SmilesParse.cpp CX wedge post-processing
        // RDKit✔️✔️: if (res->hasProp(SmilesParseOps::detail::_needsDetectAtomStereo)) {
        // RDKit✔️✔️:   res->clearProp(SmilesParseOps::detail::_needsDetectAtomStereo);
        // RDKit✔️✔️:   if (conf) {
        // RDKit✔️✔️:     MolOps::assignChiralTypesFromBondDirs(*res, conf->getId());
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION SmilesParse.cpp CX wedge post-processing
        record.properties.clear_prop("_needsDetectAtomStereo");
        if let Some(conformer) = record.coordinates.conformers_2d.first() {
            let conformer = Conformer3D::new(
                conformer.id(),
                conformer
                    .coordinates()
                    .iter()
                    .map(|point| [point[0], point[1], 0.0])
                    .collect(),
                false,
            );
            cosmolkit_core::assign_chiral_types_from_bond_dirs(
                &mut record.topology,
                &conformer,
                false,
            )
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        }
    }
    record.topology.adjacency =
        AdjacencyList::from_topology(record.topology.atoms.len(), &record.topology.bonds);
    normalize_source_coordinate_dimension(record);
    record
        .coordinates
        .validate_for_atom_count(record.topology.atoms.len())
        .map_err(model_failure)?;
    record
        .topology
        .validate()
        .map_err(|error| SmilesParseError::Model(error.to_string()))?;
    Ok(())
}
