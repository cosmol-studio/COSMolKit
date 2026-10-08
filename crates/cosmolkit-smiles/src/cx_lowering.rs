use cosmolkit_cx::{
    CxCoordinateBondKind, CxDoubleBondStereoKind, CxRecord, CxStereoGroupKind, CxWedgeDirection,
    ParsedCxExtensions,
};
use cosmolkit_model::{
    AdjacencyList, AtomId, BondDirection, BondId, BondOrder, BondStereo, ChiralTag, Conformer2D,
    Conformer3D, CoordinateDimension, PropertyValue, SGroupConnection, SGroupData, StereoGroup,
    StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock,
};

use crate::{CXSMILES_BOND_IDX_PROP, SmilesParseError, SmilesRecord};

fn cx_failure() -> SmilesParseError {
    SmilesParseError::Cx("failure parsing CXSMILES extensions".to_owned())
}

fn model_failure(error: impl std::fmt::Display) -> SmilesParseError {
    SmilesParseError::Model(error.to_string())
}

fn warn_cx_coordinate_bond_not_found(reference: cosmolkit_cx::CxBondReference) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog) << "BOND NOT FOUND! " << bidx
    // RDKit✔️✔️:                                   << " involving atom " << aidx << std::endl;
    // Reproduce this source payload/default output plus newline/flush, with
    // non-throwing I/O like the source ostream default. Independent RDLog
    // enable flags, alternate sinks and prefix formatting are unmodeled.
    use std::io::Write;
    let mut output = std::io::stderr().lock();
    let _ = writeln!(
        output,
        "BOND NOT FOUND! {} involving atom {}",
        reference.bond, reference.atom
    )
    .and_then(|_| output.flush());
}

fn warn_cx_zero_bond_not_found(index: usize) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "bond " << bondIdx
    // RDKit✔️✔️:             << " not found, cannot mark as zero order bond." << std::endl;
    // Source payload/default stderr, decimal native index, newline and flush;
    // default ostream I/O does not throw. Independent RDLog enable flags,
    // alternate sinks and prefixes remain unmodeled, not inferred here.
    use std::io::Write;
    let mut output = std::io::stderr().lock();
    let _ = writeln!(
        output,
        "bond {} not found, cannot mark as zero order bond.",
        index
    )
    .and_then(|_| output.flush());
}

fn bond_with_smiles_index(
    topology: &TopologyBlock,
    index: usize,
) -> Result<Option<BondId>, SmilesParseError> {
    // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:488-497
    // RDKit❗✔️: template <>
    // RDKit❗✔️: inline unsigned int rdvalue_cast<unsigned int>(RDValue_cast_t v) {
    // RDKit❗✔️:   if (rdvalue_is<unsigned int>(v)) {
    // RDKit❗✔️:     return v.value.u;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (rdvalue_is<int>(v)) {
    // RDKit❗✔️:     return boost::numeric_cast<unsigned int>(v.value.i);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   throw std::bad_any_cast();
    // RDKit❗✔️: }
    // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:488-497

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
    // Read only visited rows; a first match prevents all later property reads.
    for bond in &topology.bonds {
        let Some(value) = bond.prop(CXSMILES_BOND_IDX_PROP) else {
            continue;
        };
        // RDKit❗✔️: template <class T>
        // RDKit❗✔️: typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
        // RDKit❗✔️:     RDValue_cast_t arg) {
        // RDKit❗✔️:   T res;
        // RDKit❗✔️:   if (arg.getTag() == RDTypeTag::StringTag) {
        // RDKit❗✔️:     Utils::LocaleSwitcher ls;
        // RDKit❗✔️:     try {
        // RDKit❗✔️:       res = rdvalue_cast<T>(arg);
        // RDKit❗✔️:     } catch (const std::bad_any_cast &exc) {
        // RDKit❗✔️:       try {
        // RDKit❗✔️: 	std::string val = rdvalue_cast<std::string>(arg);
        // RDKit❗✔️: 	// trim only the right characters, this mimics how SD values
        // RDKit❗✔️: 	//  work on read, they will be trimmed by the MolFile parser
        // RDKit❗✔️: 	boost::trim_right(val);
        // RDKit❗✔️:         res = boost::lexical_cast<T>(val);
        // RDKit❗✔️:       } catch (...) {
        // RDKit❗✔️:         throw exc;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     res = rdvalue_cast<T>(arg);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        let source_index = cosmolkit_core::property_value_to_uint(value)
            .map_err(SmilesParseError::WriterNumeric)? as usize;
        if source_index == index {
            return Ok(Some(bond.id()));
        }
    }
    Ok(None)
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
    // BEGIN COMPLETE SF201 concrete canonical selector caller
    // RDKit❗✔️: void setStereoForBond(ROMol &mol, Bond *bond, Bond::BondStereo stereo,
    // RDKit❗✔️:                       bool useCXSmilesOrdering) {
    // RDKit❗✔️:   // NOTE:  moved from parse_doublebond_stereo CXSmilesOps
    // RDKit❗✔️:   // IF useCXSmilesOrdering is true, the cis/trans/unknown marker will be
    // RDKit❗✔️:   // assigned relative to the lowest-numbered neighbor of each double bond atom.
    // RDKit❗✔️:   // Otherwise it uses the lowest-numbered neighbor on the lower-numbered atom
    // RDKit❗✔️:   // of the double bond and the highest-numbered neighbor on the higher-numbered
    // RDKit❗✔️:   // atom
    // RDKit❗✔️:   auto begAtom = bond->getBeginAtom();
    // RDKit❗✔️:   auto endAtom = bond->getEndAtom();
    // RDKit❗✔️:   if (begAtom->getIdx() > endAtom->getIdx()) {
    // RDKit❗✔️:     std::swap(begAtom, endAtom);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    // RDKit❗✔️:     unsigned int begControl = mol.getNumAtoms();
    // RDKit❗✔️:     for (auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit❗✔️:       if (nbr == endAtom) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       begControl = std::min(nbr->getIdx(), begControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     unsigned int endControl = useCXSmilesOrdering ? mol.getNumAtoms() : 0;
    // RDKit❗✔️:     for (auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit❗✔️:       if (nbr == begAtom) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       endControl = useCXSmilesOrdering ? std::min(nbr->getIdx(), endControl)
    // RDKit❗✔️:                                        : std::max(nbr->getIdx(), endControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (begAtom != bond->getBeginAtom()) {
    // RDKit❗✔️:       std::swap(begControl, endControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     bond->setStereoAtoms(begControl, endControl);
    // RDKit❗✔️:     bond->setStereo(stereo);
    // RDKit❗✔️:     mol.setProp("_needsDetectBondStereo", 1);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END COMPLETE SF201 concrete canonical selector caller
    // CORE owns the existing neighbor/control algorithm. This adapter writes
    // the returned references, stereo and source Int flag in source order.
    let bond = record
        .topology
        .bonds
        .get(bond_id.index())
        .ok_or_else(cx_failure)?;
    let references = cosmolkit_core::double_bond_stereo_reference_atoms(
        bond_id,
        record.topology.atoms.len(),
        bond.begin(),
        bond.end(),
        true,
        |atom| atom_neighbors(&record.topology, atom).len(),
        |atom| {
            atom_neighbors(&record.topology, atom)
                .iter()
                .map(|neighbor| AtomId::new(neighbor.atom_index))
        },
    )
    .map_err(model_failure)?;
    let Some(references) = references else {
        return Ok(());
    };
    let bond = &mut record.topology.bonds[bond_id.index()];
    bond.set_stereo_atoms(Some(references));
    bond.set_stereo(stereo).map_err(model_failure)?;
    record
        .properties
        .set_prop("_needsDetectBondStereo", 1_i32)?;
    Ok(())
}

fn sgroup_kind(type_code: &[u8]) -> Option<(SubstanceGroupKind, &'static str)> {
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
        b"n" => Some((SubstanceGroupKind::StructuralRepeatUnit, "SRU")),
        b"mon" => Some((SubstanceGroupKind::Monomer, "MON")),
        b"mer" => Some((SubstanceGroupKind::Mer, "MER")),
        b"co" => Some((SubstanceGroupKind::Copolymer, "COP")),
        b"xl" => Some((SubstanceGroupKind::Crosslink, "CRO")),
        b"mod" => Some((SubstanceGroupKind::Modification, "MOD")),
        b"mix" => Some((SubstanceGroupKind::MixtureComponent, "MIX")),
        b"f" => Some((SubstanceGroupKind::Formulation, "FOR")),
        b"any" => Some((SubstanceGroupKind::AnyPolymer, "ANY")),
        b"gen" => Some((SubstanceGroupKind::Generic("GEN".into()), "GEN")),
        b"c" => Some((SubstanceGroupKind::Generic("COM".into()), "COM")),
        b"grf" => Some((SubstanceGroupKind::Graft, "GRA")),
        b"alt" | b"ran" | b"blk" => Some((SubstanceGroupKind::Copolymer, "COP")),
        _ => None,
    }
}

fn unsupported_query_label(label: &[u8]) -> bool {
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
        b"star_e" | b"Q_e" | b"QH_p" | b"AH_p" | b"X_p" | b"XH_p" | b"M_p" | b"MH_p"
    )
}

fn normalize_source_coordinate_dimension(record: &mut SmilesRecord) {
    record.coordinates.source_coordinate_dim = if record.coordinates.conformers_3d.is_empty() {
        (!record.coordinates.conformers_2d.is_empty()).then_some(CoordinateDimension::TwoD)
    } else {
        Some(CoordinateDimension::ThreeD)
    };
}

fn polymer_crossing_bond(topology: &TopologyBlock, source_index: usize) -> Option<BondId> {
    // RDKit✔️✔️: #define VALID_ATIDX(_atidx_) \
    // RDKit✔️✔️:   ((_atidx_) >= startAtomIdx && (_atidx_) < startAtomIdx + mol.getNumAtoms())
    // Behavior: canonical local startAtomIdx=0. Test only the source atom
    // window; source bond bounds belong to finalizer addBondWithIdx later.
    // Complexity: O(1), no allocation or graph access beyond atom count.
    (source_index < topology.atoms.len()).then(|| BondId::new(source_index))
}

fn finalize_polymer_sgroup(
    topology: &TopologyBlock,
    group: &mut SubstanceGroup,
    source_connect: &[u8],
    source_head: &[usize],
    source_tail: &[usize],
) -> Result<bool, SmilesParseError> {
    // BEGIN COMPLETE PINNED SF198 crossing-to-finalizer boundary
    // RDKit❗❌: bool parse_polymer_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit❗❌:                           unsigned int startAtomIdx, unsigned int nSGroups) {
    // RDKit❗❌:   // these look like:
    // RDKit❗❌:   //    |Sg:n:6,1,2,4::hh&#44;f:6,0,:4,2,|
    // RDKit❗❌:   // example from CXSMILES docs:
    // RDKit❗❌:   // the fields are:
    // RDKit❗❌:   //    Sg:[type]:[atom indices]:[subscript]:[superscript]:[head crossing
    // RDKit❗❌:   //    bonds]:[tail crossing bonds]:
    // RDKit❗❌:   //
    // RDKit❗❌:   // note that it's legit for empty fields to be completely missing.
    // RDKit❗❌:   //   for example, this doesn't have any crossing bonds indicated:
    // RDKit❗❌:   // *-CCCN-* |$star_e;;;;;star_e$,Sg:n:4,1,2,3::hh|
    // RDKit❗❌:   // this last bit makes the whole thing doubleplusfun to parse
    // RDKit❗❌:
    // RDKit❗❌:   if (first >= last || *first != 'S' || first + 2 >= last ||
    // RDKit❗❌:       *(first + 1) != 'g' || *(first + 2) != ':') {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:   first += 3;
    // RDKit❗❌:
    // RDKit❗❌:   const auto type_code = read_text_to(first, last, ":");
    // RDKit❗❌:   ++first;
    // RDKit❗❌:   const auto type = sgroupTypemap.find(type_code);
    // RDKit❗❌:   if (type == sgroupTypemap.end()) {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:   bool keepSGroup = false;
    // RDKit❗❌:   SubstanceGroup sgroup(&mol, type->second);
    // RDKit❗❌:   sgroup.setProp(cxsmilesindex, nSGroups);
    // RDKit❗❌:   if (type_code == "alt") {
    // RDKit❗❌:     sgroup.setProp("SUBTYPE", std::string("ALT"));
    // RDKit❗❌:   } else if (type_code == "ran") {
    // RDKit❗❌:     sgroup.setProp("SUBTYPE", std::string("RAN"));
    // RDKit❗❌:   } else if (type_code == "blk") {
    // RDKit❗❌:     sgroup.setProp("SUBTYPE", std::string("BLO"));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<unsigned int> atoms;
    // RDKit❗❌:   if (!read_int_list(first, last, atoms)) {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:   //++first;
    // RDKit❗❌:   for (auto idx : atoms) {
    // RDKit❗❌:     if (VALID_ATIDX(idx)) {
    // RDKit❗❌:       sgroup.addAtomWithIdx(idx - startAtomIdx);
    // RDKit❗❌:       keepSGroup = true;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<unsigned int> headCrossing;
    // RDKit❗❌:   std::vector<unsigned int> tailCrossing;
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
    // RDKit❗❌:
    // RDKit❗❌:       if (first <= last && *first == ':') {
    // RDKit❗❌:         ++first;
    // RDKit❗❌:         if (!read_int_list(first, last, headCrossing)) {
    // RDKit❗❌:           return false;
    // RDKit❗❌:         }
    // RDKit❗❌:         if (keepSGroup && !headCrossing.empty()) {
    // RDKit❗❌:           for (auto &cidx : headCrossing) {
    // RDKit❗❌:             if (VALID_ATIDX(cidx)) {
    // RDKit❗❌:               cidx -= startAtomIdx;
    // RDKit❗❌:             } else {
    // RDKit❗❌:               keepSGroup = false;
    // RDKit❗❌:               break;
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:           sgroup.setProp(_headCrossings, headCrossing, true);
    // RDKit❗❌:         }
    // RDKit❗❌:         if (first <= last && *first == ':') {
    // RDKit❗❌:           ++first;
    // RDKit❗❌:           if (!read_int_list(first, last, tailCrossing)) {
    // RDKit❗❌:             return false;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         if (keepSGroup && !tailCrossing.empty()) {
    // RDKit❗❌:           for (auto &cidx : tailCrossing) {
    // RDKit❗❌:             if (VALID_ATIDX(cidx)) {
    // RDKit❗❌:               cidx -= startAtomIdx;
    // RDKit❗❌:             } else {
    // RDKit❗❌:               keepSGroup = false;
    // RDKit❗❌:               break;
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:           sgroup.setProp("_tailCrossings", tailCrossing, true);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (keepSGroup) {  // the label processing can destroy sgroup info, so do that
    // RDKit❗❌:                      // now (the function will immediately return if already
    // RDKit❗❌:                      // called)
    // RDKit❗❌:     processCXSmilesLabels(mol);
    // RDKit❗❌:
    // RDKit❗❌:     finalizePolymerSGroup(mol, sgroup);
    // RDKit❗❌:     sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
    // RDKit❗❌:
    // RDKit❗❌:     addSubstanceGroup(mol, sgroup);
    // RDKit❗❌:   }
    // RDKit❗❌:   return true;
    // RDKit❗❌: }
    // END COMPLETE PINNED SF198 crossing-to-finalizer boundary
    let mut head = Vec::new();
    let mut tail = Vec::new();
    let mut valid = true;
    for index in source_head {
        match polymer_crossing_bond(topology, *index) {
            Some(bond) => head.push(bond),
            None => valid = false,
        }
    }
    for index in source_tail {
        match polymer_crossing_bond(topology, *index) {
            Some(bond) => tail.push(bond),
            None => valid = false,
        }
    }
    if !valid {
        return Ok(false);
    }
    // BEGIN COMPLETE RDProps::setProp computed crossing bookkeeping
    // RDKit✔️❌: void setProp(const std::string_view key, T val, bool computed = false) const {
    // RDKit✔️❌:     if(key.empty()) {
    // RDKit✔️❌:       throw ValueErrorException("Cannot set property with empty key");
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (computed) {
    // RDKit✔️❌:       STR_VECT compLst;
    // RDKit✔️❌:       getPropIfPresent(RDKit::detail::computedPropName, compLst);
    // RDKit✔️❌:       if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
    // RDKit✔️❌:         compLst.emplace_back(key);
    // RDKit✔️❌:         d_props.setVal(RDKit::detail::computedPropName, compLst);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     d_props.setVal(key, val);
    // RDKit✔️❌:   }
    // END COMPLETE RDProps::setProp computed crossing bookkeeping
    // Source computed crossing payloads travel through the existing typed
    // head/tail arguments. The modeled StringVector bookkeeping remains an
    // actual ordered property record; CORE clearProp removes these names and
    // leaves the source empty reserved entry, including its insertion position.
    // No opaque unsigned-vector tag or parallel computed classification store.
    let mut computed_crossings = Vec::<cosmolkit_model::PropertyText>::new();
    if !head.is_empty() {
        computed_crossings.push("_headCrossings".into());
        group.set_prop("__computedProps", computed_crossings.clone())?;
    }
    if !tail.is_empty() {
        computed_crossings.push("_tailCrossings".into());
        group.set_prop("__computedProps", computed_crossings)?;
    }
    // CX parsing stores CONNECT only when the source superscript is nonempty.
    // The core helper owns source normalization, inferred crossings, and the
    // ordered typed SGroup updates shared with the SMARTS lowerer.
    cosmolkit_core::finalize_polymer_sgroup(
        group,
        (!source_connect.is_empty()).then_some(source_connect),
        &head,
        &tail,
        topology.atoms.len(),
        topology.bonds.len(),
        |atom| {
            topology
                .adjacency
                .neighbors_of(atom.index())
                .iter()
                .map(|neighbor| (AtomId::new(neighbor.atom_index), neighbor.bond))
        },
    )
    .map_err(model_failure)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(input: &str) -> crate::SmilesRecord {
        crate::parse_smiles(input, &Default::default())
            .unwrap_or_else(|error| panic!("failed to parse {input:?}: {error}"))
    }

    #[test]
    fn query_sgroups_smiles_polymer_accepted_crossings_keep_existing_output() {
        let record = parse("CCCCC |Sg:n:1,2,3:repeat:ht:0,0,3:3,3,0:|");
        let group = &record.topology.substance_groups[0];
        assert_eq!(
            group.bonds(),
            &[
                BondId::new(0),
                BondId::new(0),
                BondId::new(3),
                BondId::new(3),
                BondId::new(3),
                BondId::new(0),
            ]
        );
        assert_eq!(
            group.head_crossing_bonds(),
            &[BondId::new(0), BondId::new(0), BondId::new(3)]
        );
        assert_eq!(
            group.crossing_bond_correspondence(),
            &[
                BondId::new(0),
                BondId::new(3),
                BondId::new(0),
                BondId::new(3),
                BondId::new(3),
                BondId::new(0),
            ]
        );
        assert_eq!(
            group.props().get(b"CONNECT".as_slice()),
            Some(&cosmolkit_model::PropertyValue::String("HT".into()))
        );

        let output = crate::write_cx_smiles(&record).expect("accepted SGroup remains writable");
        assert!(
            output
                .as_bytes()
                .windows(b"Sg:n:1,2,3:repeat:ht:0,0,3:3,3,0:".len())
                .any(|bytes| bytes == b"Sg:n:1,2,3:repeat:ht:0,0,3:3,3,0:"),
            "{output:?}"
        );
    }

    #[test]
    fn query_sgroups_smiles_polymer_uppercase_connect_uses_source_fallback() {
        let record = parse("CCCCC |Sg:n:1,2,3:repeat:HH:0,0,3:2,3,1:|");
        let group = &record.topology.substance_groups[0];
        assert_eq!(
            group.connection(),
            Some(&cosmolkit_model::SGroupConnection::Either)
        );
        assert_eq!(
            group.props().get(b"CONNECT".as_slice()),
            Some(&cosmolkit_model::PropertyValue::String("EU".into()))
        );
        assert_eq!(
            group.crossing_bond_correspondence(),
            &[
                BondId::new(0),
                BondId::new(2),
                BondId::new(0),
                BondId::new(3),
                BondId::new(3),
                BondId::new(1),
            ]
        );
    }

    #[test]
    fn query_sgroups_smiles_polymer_repeated_flip_markers_reverse_tail_pairs() {
        let record = parse("CCCCC |Sg:n:1,2,3:repeat:hh&#44;f&#44;f:0,0,3:2,3,1:|");
        let group = &record.topology.substance_groups[0];
        assert_eq!(
            group.connection(),
            Some(&cosmolkit_model::SGroupConnection::HeadToHead)
        );
        assert_eq!(
            group.props().get(b"CONNECT".as_slice()),
            Some(&cosmolkit_model::PropertyValue::String("HH".into()))
        );
        assert_eq!(
            group.crossing_bond_correspondence(),
            &[
                BondId::new(0),
                BondId::new(1),
                BondId::new(0),
                BondId::new(3),
                BondId::new(3),
                BondId::new(2),
            ]
        );
    }

    #[test]
    fn query_sgroups_smiles_polymer_keeps_source_crossing_index_filter() {
        let dropped = parse("CCCCC |Sg:n:1,2,3::ht:9:|");
        assert!(dropped.topology.substance_groups.is_empty());

        let bond_out_of_range =
            crate::parse_smiles("CCCCC |Sg:n:1,2,3::ht:4:|", &Default::default())
                .expect_err("source-valid atom index is still checked as a bond index");
        assert_eq!(
            bond_out_of_range,
            crate::SmilesParseError::Model(
                "polymer SGroup bond 4 is outside the topology with 4 bonds".into()
            )
        );
    }
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
                    .any(|label| unsupported_query_label(label.as_bytes())) =>
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
    let mut sgroup_index = 0_u32;
    let mut stereo_hashes = Vec::new();
    for item in parsed.records() {
        match item {
            CxRecord::Coordinates(coordinates) => {
                // BEGIN COMPLETE PINNED SF188
                // RDKit✔️❌: bool parse_coords(Iterator &first, Iterator last, RDKit::RWMol &mol,
                // RDKit✔️❌:                   unsigned int startAtomIdx, unsigned int confIdx) {
                // RDKit✔️❌:   if (first >= last || *first != '(') {
                // RDKit✔️❌:     return false;
                // RDKit✔️❌:   }
                // RDKit✔️❌:
                // RDKit✔️❌:   auto *conf = new Conformer(mol.getNumAtoms());
                // RDKit✔️❌:   mol.addConformer(conf);
                // RDKit✔️❌:   conf->setId(confIdx);
                // RDKit✔️❌:   ++first;
                // RDKit✔️❌:   unsigned int atIdx = 0;
                // RDKit✔️❌:   bool is3D = false;
                // RDKit✔️❌:   while (first <= last && *first != ')') {
                // RDKit✔️❌:     RDGeom::Point3D pt;
                // RDKit✔️❌:     std::string tkn = read_text_to(first, last, ";)");
                // RDKit✔️❌:     if (VALID_ATIDX(atIdx)) {
                // RDKit✔️❌:       if (!tkn.empty()) {
                // RDKit✔️❌:         std::vector<std::string> tokens;
                // RDKit✔️❌:         boost::split(tokens, tkn, boost::is_any_of(std::string(",")));
                // RDKit✔️❌:         if (tokens.size() >= 1 && tokens[0].size()) {
                // RDKit✔️❌:           pt.x = boost::lexical_cast<double>(tokens[0]);
                // RDKit✔️❌:         }
                // RDKit✔️❌:         if (tokens.size() >= 2 && tokens[1].size()) {
                // RDKit✔️❌:           pt.y = boost::lexical_cast<double>(tokens[1]);
                // RDKit✔️❌:         }
                // RDKit✔️❌:         if (tokens.size() >= 3 && tokens[2].size()) {
                // RDKit✔️❌:           pt.z = boost::lexical_cast<double>(tokens[2]);
                // RDKit✔️❌:           is3D = true;
                // RDKit✔️❌:         }
                // RDKit✔️❌:       }
                // RDKit✔️❌:
                // RDKit✔️❌:       conf->setAtomPos(atIdx - startAtomIdx, pt);
                // RDKit✔️❌:     }
                // RDKit✔️❌:     ++atIdx;
                // RDKit✔️❌:     if (first <= last && *first != ')') {
                // RDKit✔️❌:       ++first;
                // RDKit✔️❌:     }
                // RDKit✔️❌:   }
                // RDKit✔️❌:   // make sure that the conformer really is 3D!
                // RDKit✔️❌:   if (is3D && hasNonZeroZCoords(*conf)) {
                // RDKit✔️❌:     conf->set3D(true);
                // RDKit✔️❌:   } else {
                // RDKit✔️❌:     conf->set3D(false);
                // RDKit✔️❌:   }
                // RDKit✔️❌:   if (first >= last || *first != ')') {
                // RDKit✔️❌:     return false;
                // RDKit✔️❌:   }
                // RDKit✔️❌:   ++first;
                // RDKit✔️❌:   return true;
                // RDKit✔️❌: }
                // END COMPLETE PINNED SF188
                let source_atom_count =
                    u32::try_from(atom_count).map_err(|_| SmilesParseError::Syntax {
                        offset: 0,
                        message: "CX atom count exceeds source unsigned32 domain".to_owned(),
                    })?;
                let mut points = vec![[0.0; 3]; atom_count];
                for (slot, value) in coordinates.values.iter().enumerate() {
                    let index = slot as u32;
                    if index < source_atom_count {
                        points[index as usize] = value.unwrap_or([0.0; 3]);
                    }
                }
                // Source Conformer always owns Point3D; set3D(false) does not
                // erase z, change numerical bits or move it to projected XY.
                record
                    .coordinates
                    .record_source_conformer_append(cosmolkit_model::CoordinateDimension::ThreeD)
                    .map_err(SmilesParseError::Coordinates)?;
                // This complete-record lowerer also accepts syntax-only records;
                // derive the final flag from the actual mapped destination rows.
                let is_3d = coordinates.is_3d && points.iter().any(|point| point[2].abs() > 1e-3);
                record.coordinates.conformers_3d.push(Conformer3D::new(
                    coordinates.conformer,
                    points,
                    is_3d,
                ));
            }
            CxRecord::AtomLabels(values) => {
                // BEGIN COMPLETE PINNED SF187 GRAPH WRITE
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
                // END COMPLETE PINNED SF187 GRAPH WRITE
                // RDKit✔️✔️: mol.getAtomWithIdx(atIdx - startAtomIdx)
                // RDKit✔️✔️:     ->setProp(RDKit::common_properties::atomLabel, tkn);
                // RDKit✔️✔️: unsigned int atIdx = 0;
                // RDKit✔️✔️: ++atIdx;
                // Transport slot remains usize; source atom ordinal is u32.
                for (slot, value) in values.iter().enumerate() {
                    let index = (slot as u32) as usize;
                    if let Some(value) = value.as_ref().filter(|value| !value.is_empty())
                        && let Some(atom) = record.topology.atoms.get_mut(index)
                    {
                        atom.set_prop("atomLabel", value)?;
                        if matches!(value.as_bytes(), b"Pol_p" | b"Mod_p") {
                            // RDKit✔️✔️: atom->setProp(common_properties::dummyLabel,
                            // RDKit✔️✔️:               symb.substr(0, symb.size() - 2));
                            // RDKit✔️✔️: atom->clearProp(common_properties::atomLabel);
                            atom.set_prop(
                                "dummyLabel",
                                cosmolkit_model::PropertyText::from_bytes(
                                    &value.as_bytes()[..value.len() - 2],
                                ),
                            )?;
                            atom.clear_prop("atomLabel")?;
                        }
                    }
                }
            }
            CxRecord::AtomValues(values) => {
                // BEGIN COMPLETE PINNED SF185 GRAPH WRITE
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
                // END COMPLETE PINNED SF185 GRAPH WRITE
                // RDKit✔️✔️: mol.getAtomWithIdx(atIdx)->setProp(
                // RDKit✔️✔️:     RDKit::common_properties::molFileValue, tkn);
                // RDKit✔️✔️: unsigned int atIdx = 0;
                // RDKit✔️✔️: ++atIdx;
                // Record-slot ordinals are transport positions; chemistry uses
                // the source 32-bit unsigned atom index, including wraparound.
                for (slot, value) in values.iter().enumerate() {
                    let index = (slot as u32) as usize;
                    if let Some(value) = value.as_ref().filter(|value| !value.is_empty())
                        && let Some(atom) = record.topology.atoms.get_mut(index)
                    {
                        atom.set_prop("molFileValue", value)?;
                    }
                }
            }
            CxRecord::AtomProperties(properties) => {
                // RDKit✔️✔️: mol.getAtomWithIdx(atIdx - startAtomIdx)->setProp(pname, pval);
                for property in properties {
                    // RDKit✔️✔️: if (!pname.empty()) {
                    // RDKit✔️✔️:   if (VALID_ATIDX(atIdx) && !pval.empty()) {
                    if property.name.is_empty() || property.value.is_empty() {
                        continue;
                    }
                    if let Some(atom) = record.topology.atoms.get_mut(property.atom) {
                        atom.set_prop(property.name.clone(), property.value.clone())?;
                    }
                }
            }
            CxRecord::CoordinateBonds(annotation) => {
                // BEGIN COMPLETE PINNED SF189
                // RDKit✔️❌: bool parse_coordinate_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
                // RDKit✔️❌:                             Bond::BondType typ, unsigned int startAtomIdx,
                // RDKit✔️❌:                             unsigned int startBondIdx) {
                // RDKit✔️❌:   if (first >= last || (*first != 'C' && *first != 'H')) {
                // RDKit✔️❌:     return false;
                // RDKit✔️❌:   }
                // RDKit✔️❌:   ++first;
                // RDKit✔️❌:   if (first >= last || *first != ':') {
                // RDKit✔️❌:     return false;
                // RDKit✔️❌:   }
                // RDKit✔️❌:   ++first;
                // RDKit✔️❌:   while (first <= last && *first >= '0' && *first <= '9') {
                // RDKit✔️❌:     unsigned int aidx;
                // RDKit✔️❌:     unsigned int bidx;
                // RDKit✔️❌:     if (read_int_pair(first, last, aidx, bidx)) {
                // RDKit✔️❌:       if (VALID_ATIDX(aidx) && VALID_BNDIDX(bidx)) {
                // RDKit✔️❌:         auto bnd = get_bond_with_smiles_idx(mol, bidx - startBondIdx);
                // RDKit✔️❌:         if (!bnd || (bnd->getBeginAtomIdx() != aidx - startAtomIdx &&
                // RDKit✔️❌:                      bnd->getEndAtomIdx() != aidx - startAtomIdx)) {
                // RDKit✔️❌:           BOOST_LOG(rdWarningLog) << "BOND NOT FOUND! " << bidx
                // RDKit✔️❌:                                   << " involving atom " << aidx << std::endl;
                // RDKit✔️❌:           return false;
                // RDKit✔️❌:         }
                // RDKit✔️❌:         bnd->setBondType(typ);
                // RDKit✔️❌:         if (bnd->getBeginAtomIdx() != aidx - startAtomIdx) {
                // RDKit✔️❌:           unsigned int tmp = bnd->getBeginAtomIdx();
                // RDKit✔️❌:           bnd->setBeginAtomIdx(aidx - startAtomIdx);
                // RDKit✔️❌:           bnd->setEndAtomIdx(tmp);
                // RDKit✔️❌:         }
                // RDKit✔️❌:       }
                // RDKit✔️❌:     } else {
                // RDKit✔️❌:       return false;
                // RDKit✔️❌:     }
                // RDKit✔️❌:     if (first < last && *first == ',') {
                // RDKit✔️❌:       ++first;
                // RDKit✔️❌:     }
                // RDKit✔️❌:   }
                // RDKit✔️❌:   return true;
                // RDKit✔️❌: }
                // END COMPLETE PINNED SF189
                // BEGIN COMPLETE REACHED Code/GraphMol/Bond.h void setBondType(BondType bT)
                // RDKit✔️✔️:   void setBondType(BondType bT) { d_bondType = bT; }
                // END COMPLETE REACHED Code/GraphMol/Bond.h void setBondType(BondType bT)
                // BEGIN COMPLETE REACHED Code/GraphMol/Bond.h unsigned int getBeginAtomIdx() const
                // RDKit✔️✔️:   unsigned int getBeginAtomIdx() const { return d_beginAtomIdx; }
                // END COMPLETE REACHED Code/GraphMol/Bond.h unsigned int getBeginAtomIdx() const
                // BEGIN COMPLETE REACHED Code/GraphMol/Bond.h unsigned int getEndAtomIdx() const
                // RDKit✔️✔️:   unsigned int getEndAtomIdx() const { return d_endAtomIdx; }
                // END COMPLETE REACHED Code/GraphMol/Bond.h unsigned int getEndAtomIdx() const
                // BEGIN COMPLETE REACHED Code/GraphMol/Bond.cpp void Bond::setBeginAtomIdx(unsigned int what)
                // RDKit✔️✔️: void Bond::setBeginAtomIdx(unsigned int what) {
                // RDKit✔️✔️:   if (dp_mol) {
                // RDKit✔️✔️:     URANGE_CHECK(what, getOwningMol().getNumAtoms());
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   d_beginAtomIdx = what;
                // RDKit✔️✔️: }
                // END COMPLETE REACHED Code/GraphMol/Bond.cpp void Bond::setBeginAtomIdx(unsigned int what)
                // BEGIN COMPLETE REACHED Code/GraphMol/Bond.cpp void Bond::setEndAtomIdx(unsigned int what)
                // RDKit✔️✔️: void Bond::setEndAtomIdx(unsigned int what) {
                // RDKit✔️✔️:   if (dp_mol) {
                // RDKit✔️✔️:     URANGE_CHECK(what, getOwningMol().getNumAtoms());
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   d_endAtomIdx = what;
                // RDKit✔️✔️: }
                // END COMPLETE REACHED Code/GraphMol/Bond.cpp void Bond::setEndAtomIdx(unsigned int what)
                // Endpoints are already known valid after the range/membership
                // checks; both native setters pass URANGE_CHECK. The detached
                // setter preserves their assignments without another failure.
                let order = match annotation.kind {
                    CxCoordinateBondKind::Dative => BondOrder::Dative,
                    CxCoordinateBondKind::Hydrogen => BondOrder::Hydrogen,
                };
                for reference in &annotation.bonds {
                    let atom = AtomId::new(reference.atom);
                    if atom.index() >= atom_count || reference.bond >= record.topology.bonds.len() {
                        continue;
                    }
                    let bond_id = match bond_with_smiles_index(&record.topology, reference.bond)? {
                        Some(bond_id) => bond_id,
                        None => {
                            warn_cx_coordinate_bond_not_found(*reference);
                            return Err(cx_failure());
                        }
                    };
                    let (begin, end) = record
                        .topology
                        .bonds
                        .get(bond_id.index())
                        .map(|bond| (bond.begin(), bond.end()))
                        .ok_or_else(cx_failure)?;
                    if begin != atom && end != atom {
                        warn_cx_coordinate_bond_not_found(*reference);
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
                // BEGIN COMPLETE PINNED SF190
                // RDKit✔️❌: bool parse_zero_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
                // RDKit✔️❌:                       unsigned int, unsigned int startBondIdx) {
                // RDKit✔️❌:   // these look like: C1CCCCC~CCCC1 |Z:5|
                // RDKit✔️❌:   if (first >= last || *first != 'Z') {
                // RDKit✔️❌:     return false;
                // RDKit✔️❌:   }
                // RDKit✔️❌:   ++first;
                // RDKit✔️❌:   if (first >= last || *first != ':') {
                // RDKit✔️❌:     return false;
                // RDKit✔️❌:   }
                // RDKit✔️❌:   ++first;
                // RDKit✔️❌:
                // RDKit✔️❌:   while (first < last && *first >= '0' && *first <= '9') {
                // RDKit✔️❌:     unsigned int bondIdx;
                // RDKit✔️❌:     if (!read_int(first, last, bondIdx)) {
                // RDKit✔️❌:       return false;
                // RDKit✔️❌:     }
                // RDKit✔️❌:     if (VALID_BNDIDX(bondIdx)) {
                // RDKit✔️❌:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
                // RDKit✔️❌:
                // RDKit✔️❌:       if (!bond) {
                // RDKit✔️❌:         BOOST_LOG(rdWarningLog)
                // RDKit✔️❌:             << "bond " << bondIdx
                // RDKit✔️❌:             << " not found, cannot mark as zero order bond." << std::endl;
                // RDKit✔️❌:         return false;
                // RDKit✔️❌:       }
                // RDKit✔️❌:       bond->setBondType(Bond::ZERO);
                // RDKit✔️❌:     }
                // RDKit✔️❌:     if (first < last && *first == ',') {
                // RDKit✔️❌:       ++first;
                // RDKit✔️❌:     }
                // RDKit✔️❌:   }
                // RDKit✔️❌:   return true;
                // RDKit✔️❌: }
                // END COMPLETE PINNED SF190
                for index in indices {
                    if *index >= record.topology.bonds.len() {
                        continue;
                    }
                    let bond_id = match bond_with_smiles_index(&record.topology, *index)? {
                        Some(bond_id) => bond_id,
                        None => {
                            warn_cx_zero_bond_not_found(*index);
                            return Err(cx_failure());
                        }
                    };
                    record.topology.bonds[bond_id.index()].set_order(BondOrder::Zero);
                }
            }
            CxRecord::EnhancedStereo(stereo) => {
                // BEGIN COMPLETE PINNED SF205 concrete caller
                // RDKit❗✔️: bool parse_enhanced_stereo(Iterator &first, Iterator last, RDKit::RWMol &mol,
                // RDKit❗✔️:                            unsigned int startAtomIdx) {
                // RDKit❗✔️:   StereoGroupType group_type = StereoGroupType::STEREO_ABSOLUTE;
                // RDKit❗✔️:   if (*first == 'a') {
                // RDKit❗✔️:     group_type = StereoGroupType::STEREO_ABSOLUTE;
                // RDKit❗✔️:   } else if (*first == 'o') {
                // RDKit❗✔️:     group_type = StereoGroupType::STEREO_OR;
                // RDKit❗✔️:   } else if (*first == '&') {
                // RDKit❗✔️:     group_type = StereoGroupType::STEREO_AND;
                // RDKit❗✔️:   }
                // RDKit❗✔️:   ++first;
                // RDKit❗✔️:
                // RDKit❗✔️:   // OR and AND groups carry a group number
                // RDKit❗✔️:   unsigned int group_id = 0;
                // RDKit❗✔️:   if (group_type != StereoGroupType::STEREO_ABSOLUTE) {
                // RDKit❗✔️:     read_int(first, last, group_id);
                // RDKit❗✔️:   }
                // RDKit❗✔️:
                // RDKit❗✔️:   if (first >= last || *first != ':') {
                // RDKit❗✔️:     return false;
                // RDKit❗✔️:   }
                // RDKit❗✔️:   ++first;
                // RDKit❗✔️:
                // RDKit❗✔️:   std::vector<Atom *> atoms;
                // RDKit❗✔️:   std::vector<Bond *> bonds;
                // RDKit❗✔️:
                // RDKit❗✔️:   while (first <= last && *first >= '0' && *first <= '9') {
                // RDKit❗✔️:     unsigned int aidx;
                // RDKit❗✔️:     if (read_int(first, last, aidx)) {
                // RDKit❗✔️:       if (VALID_ATIDX(aidx)) {
                // RDKit❗✔️:         Atom *atom = mol.getAtomWithIdx(aidx - startAtomIdx);
                // RDKit❗✔️:         if (!atom) {
                // RDKit❗✔️:           BOOST_LOG(rdWarningLog)
                // RDKit❗✔️:               << "Atom " << aidx << " not found!" << std::endl;
                // RDKit❗✔️:           return false;
                // RDKit❗✔️:         }
                // RDKit❗✔️:         atoms.push_back(atom);
                // RDKit❗✔️:       }
                // RDKit❗✔️:     } else {
                // RDKit❗✔️:       return false;
                // RDKit❗✔️:     }
                // RDKit❗✔️:
                // RDKit❗✔️:     if (first < last && *first == ',') {
                // RDKit❗✔️:       ++first;
                // RDKit❗✔️:     }
                // RDKit❗✔️:   }
                // RDKit❗✔️:   if (!atoms.empty()) {
                // RDKit❗✔️:     // we need to do a bit of work to check whether or not we've already seen
                // RDKit❗✔️:     // this particular StereoGroup (was Github #6050)
                // RDKit❗✔️:     const auto group_hash =
                // RDKit❗✔️:         10 * group_id + static_cast<unsigned int>(group_type);
                // RDKit❗✔️:     std::vector<unsigned int> sgTracker;
                // RDKit❗✔️:     mol.getPropIfPresent(cxsgTracker, sgTracker);
                // RDKit❗✔️:     std::vector<StereoGroup> mol_stereo_groups(mol.getStereoGroups());
                // RDKit❗✔️:     TEST_ASSERT(mol_stereo_groups.size() == sgTracker.size());
                // RDKit❗✔️:
                // RDKit❗✔️:     auto iter = std::find(sgTracker.begin(), sgTracker.end(), group_hash);
                // RDKit❗✔️:     if (iter != sgTracker.end()) {
                // RDKit❗✔️:       auto index = iter - sgTracker.begin();
                // RDKit❗✔️:       auto gAtoms = mol_stereo_groups[index].getAtoms();
                // RDKit❗✔️:       gAtoms.insert(gAtoms.end(), atoms.begin(), atoms.end());
                // RDKit❗✔️:       mol_stereo_groups[index] =
                // RDKit❗✔️:           StereoGroup(mol_stereo_groups[index].getGroupType(),
                // RDKit❗✔️:                       std::move(gAtoms), std::move(bonds), group_id);
                // RDKit❗✔️:     } else {
                // RDKit❗✔️:       // not seen this before, create a new stereogroup
                // RDKit❗✔️:       mol_stereo_groups.emplace_back(group_type, std::move(atoms),
                // RDKit❗✔️:                                      std::move(bonds), group_id);
                // RDKit❗✔️:       sgTracker.push_back(group_hash);
                // RDKit❗✔️:       mol.setProp(cxsgTracker, sgTracker);
                // RDKit❗✔️:     }
                // RDKit❗✔️:
                // RDKit❗✔️:     mol.setStereoGroups(std::move(mol_stereo_groups));
                // RDKit❗✔️:   }
                // RDKit❗✔️:
                // RDKit❗✔️:   return true;
                // RDKit❗✔️: }
                // END COMPLETE PINNED SF205 concrete caller
                let kind = match stereo.kind {
                    CxStereoGroupKind::Absolute => StereoGroupKind::Absolute,
                    CxStereoGroupKind::Or => StereoGroupKind::Or,
                    CxStereoGroupKind::And => StereoGroupKind::And,
                };
                let atoms = stereo
                    .atoms
                    .iter()
                    .filter(|&&index| index < atom_count)
                    .map(|&index| AtomId::new(index))
                    .collect();
                let groups = cosmolkit_core::parser_stereo_order::merge_cx_stereo_group(
                    &record.topology.stereo_groups,
                    &mut stereo_hashes,
                    record.properties.prop(b"_sgTracker"),
                    kind,
                    stereo.group_id,
                    atoms,
                )
                .map_err(model_failure)?;
                if let Some(groups) = groups {
                    record.topology.stereo_groups = groups;
                }
            }
            CxRecord::WedgedBonds(wedges) => {
                // BEGIN COMPLETE PINNED SF200 concrete caller
                // RDKit❗❌: bool parse_wedged_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
                // RDKit❗❌:                         unsigned int startAtomIdx, unsigned int startBondIdx) {
                // RDKit❗❌:   // these look like: CC(O)Cl |w:1.0|
                // RDKit❗❌:   // also wD and wU for down and up wedges.
                // RDKit❗❌:   //
                // RDKit❗❌:   // We do not end up using this to set stereochemistry, but the relevant bond
                // RDKit❗❌:   // properties are set in case client code wants to do something with the
                // RDKit❗❌:   // information.
                // RDKit❗❌:   if (first >= last || *first != 'w' || first + 1 >= last) {
                // RDKit❗❌:     return false;
                // RDKit❗❌:   }
                // RDKit❗❌:   ++first;
                // RDKit❗❌:   Bond::BondDir state = Bond::BondDir::NONE;
                // RDKit❗❌:   unsigned int cfg = 0;
                // RDKit❗❌:   switch (*first) {
                // RDKit❗❌:     case ':':
                // RDKit❗❌:       state = Bond::BondDir::UNKNOWN;
                // RDKit❗❌:       cfg = 2;
                // RDKit❗❌:       break;
                // RDKit❗❌:     case 'U':
                // RDKit❗❌:       state = Bond::BondDir::BEGINWEDGE;
                // RDKit❗❌:       cfg = 1;
                // RDKit❗❌:       ++first;
                // RDKit❗❌:       break;
                // RDKit❗❌:     case 'D':
                // RDKit❗❌:       state = Bond::BondDir::BEGINDASH;
                // RDKit❗❌:       cfg = 3;
                // RDKit❗❌:       ++first;
                // RDKit❗❌:       break;
                // RDKit❗❌:     default:
                // RDKit❗❌:       break;
                // RDKit❗❌:   }
                // RDKit❗❌:   if (state == Bond::BondDir::NONE || first >= last || first + 1 >= last ||
                // RDKit❗❌:       *first != ':') {
                // RDKit❗❌:     return false;
                // RDKit❗❌:   }
                // RDKit❗❌:   ++first;
                // RDKit❗❌:   while (first < last && *first >= '0' && *first <= '9') {
                // RDKit❗❌:     unsigned int atomIdx;
                // RDKit❗❌:     if (!read_int(first, last, atomIdx)) {
                // RDKit❗❌:       return false;
                // RDKit❗❌:     }
                // RDKit❗❌:     if (first < last && *first == '.') {
                // RDKit❗❌:       ++first;
                // RDKit❗❌:     } else {
                // RDKit❗❌:       BOOST_LOG(rdWarningLog) << "improperly formatted w block" << std::endl;
                // RDKit❗❌:       return false;
                // RDKit❗❌:     }
                // RDKit❗❌:     unsigned int bondIdx;
                // RDKit❗❌:     if (!read_int(first, last, bondIdx)) {
                // RDKit❗❌:       return false;
                // RDKit❗❌:     }
                // RDKit❗❌:
                // RDKit❗❌:     if (VALID_ATIDX(atomIdx) && VALID_BNDIDX(bondIdx)) {
                // RDKit❗❌:       auto atom = mol.getAtomWithIdx(atomIdx - startAtomIdx);
                // RDKit❗❌:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
                // RDKit❗❌:
                // RDKit❗❌:       if (!bond) {
                // RDKit❗❌:         BOOST_LOG(rdWarningLog)
                // RDKit❗❌:             << "bond " << bondIdx << " not found, wedge from atom " << atomIdx
                // RDKit❗❌:             << " cannot be applied." << std::endl;
                // RDKit❗❌:         return false;
                // RDKit❗❌:       }
                // RDKit❗❌:
                // RDKit❗❌:       // we can't set wedging twice:
                // RDKit❗❌:       if (bond->hasProp(common_properties::_MolFileBondCfg)) {
                // RDKit❗❌:         BOOST_LOG(rdWarningLog)
                // RDKit❗❌:             << "w block attempts to set wedging on bond " << bond->getIdx()
                // RDKit❗❌:             << " more than once." << std::endl;
                // RDKit❗❌:         return false;
                // RDKit❗❌:       }
                // RDKit❗❌:
                // RDKit❗❌:       // first things first, the atom needs to be the start atom of the bond for
                // RDKit❗❌:       // any of this to make sense
                // RDKit❗❌:       if (atom->getIdx() != bond->getBeginAtomIdx()) {
                // RDKit❗❌:         if (atom->getIdx() != bond->getEndAtomIdx()) {
                // RDKit❗❌:           BOOST_LOG(rdWarningLog)
                // RDKit❗❌:               << "atom " << atomIdx << " is not associated with bond "
                // RDKit❗❌:               << bondIdx << "(" << bond->getBeginAtomIdx() + startAtomIdx << "-"
                // RDKit❗❌:               << bond->getEndAtomIdx() + startAtomIdx << ")"
                // RDKit❗❌:               << " in w block" << std::endl;
                // RDKit❗❌:           return false;
                // RDKit❗❌:         }
                // RDKit❗❌:         auto eidx = bond->getBeginAtomIdx();
                // RDKit❗❌:         bond->setBeginAtomIdx(atom->getIdx());
                // RDKit❗❌:         bond->setEndAtomIdx(eidx);
                // RDKit❗❌:       }
                // RDKit❗❌:       bond->setProp(common_properties::_MolFileBondCfg, cfg);
                // RDKit❗❌:       bond->setBondDir(state);
                // RDKit❗❌:       if (cfg == 2 && canHaveDirection(*bond)) {
                // RDKit❗❌:         bond->getBeginAtom()->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
                // RDKit❗❌:         mol.setProp(detail::_needsDetectBondStereo, 1);
                // RDKit❗❌:       }
                // RDKit❗❌:       if ((cfg == 1 || cfg == 3) && canHaveDirection(*bond)) {
                // RDKit❗❌:         mol.setProp(detail::_needsDetectAtomStereo, 1);
                // RDKit❗❌:       }
                // RDKit❗❌:     }
                // RDKit❗❌:     if (first < last && *first == ',') {
                // RDKit❗❌:       ++first;
                // RDKit❗❌:     }
                // RDKit❗❌:   }
                // RDKit❗❌:   return true;
                // RDKit❗❌: }
                // END COMPLETE PINNED SF200 concrete caller
                // BEGIN COMPLETE SF200 concrete reached helpers
                // RDKit❗❌: Bond *get_bond_with_smiles_idx(const ROMol &mol, unsigned idx) {
                // RDKit❗❌:   for (auto bnd : mol.bonds()) {
                // RDKit❗❌:     unsigned int smilesIdx;
                // RDKit❗❌:     if (bnd->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
                // RDKit❗❌:         smilesIdx == idx) {
                // RDKit❗❌:       return bnd;
                // RDKit❗❌:     }
                // RDKit❗❌:   }
                // RDKit❗❌:   return nullptr;
                // RDKit❗❌: }
                // RDKit❗❌:
                // RDKit❗❌: inline bool canHaveDirection(const Bond &bond) {
                // RDKit❗❌:   auto bondType = bond.getBondType();
                // RDKit❗❌:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
                // RDKit❗❌: }
                // RDKit❗❌:
                // RDKit❗❌: void Bond::setBeginAtomIdx(unsigned int what) {
                // RDKit❗❌:   if (dp_mol) {
                // RDKit❗❌:     URANGE_CHECK(what, getOwningMol().getNumAtoms());
                // RDKit❗❌:   }
                // RDKit❗❌:   d_beginAtomIdx = what;
                // RDKit❗❌: }
                // RDKit❗❌:
                // RDKit❗❌: void Bond::setEndAtomIdx(unsigned int what) {
                // RDKit❗❌:   if (dp_mol) {
                // RDKit❗❌:     URANGE_CHECK(what, getOwningMol().getNumAtoms());
                // RDKit❗❌:   }
                // RDKit❗❌:   d_endAtomIdx = what;
                // RDKit❗❌: }
                // END COMPLETE SF200 concrete reached helpers
                for wedge in wedges {
                    if wedge.atom >= atom_count || wedge.bond >= record.topology.bonds.len() {
                        continue;
                    }
                    let bond_id = match bond_with_smiles_index(&record.topology, wedge.bond)? {
                        Some(id) => id,
                        None => {
                            eprintln!(
                                "bond {} not found, wedge from atom {} cannot be applied.",
                                wedge.bond, wedge.atom
                            );
                            return Err(cx_failure());
                        }
                    };
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
                        eprintln!(
                            "w block attempts to set wedging on bond {} more than once.",
                            bond_id.index()
                        );
                        return Err(cx_failure());
                    }
                    let atom = AtomId::new(wedge.atom);
                    if begin != atom && end != atom {
                        eprintln!(
                            "atom {} is not associated with bond {}({}-{}) in w block",
                            wedge.atom,
                            wedge.bond,
                            begin.index(),
                            end.index()
                        );
                        return Err(cx_failure());
                    }
                    let (cfg, direction) = match wedge.direction {
                        CxWedgeDirection::Unknown => (2_u32, BondDirection::Unknown),
                        CxWedgeDirection::BeginWedge => (1_u32, BondDirection::BeginWedge),
                        CxWedgeDirection::BeginDash => (3_u32, BondDirection::BeginDash),
                    };
                    let bond = &mut record.topology.bonds[bond_id.index()];
                    if begin != atom {
                        bond.set_endpoints(atom, begin);
                    }
                    bond.set_prop("_MolFileBondCfg", cfg)?;
                    bond.set_direction(direction);
                    if cfg == 2 && can_have_direction(order) {
                        record.topology.atoms[atom.index()].set_chiral_tag(ChiralTag::Unspecified);
                        record
                            .properties
                            .set_prop("_needsDetectBondStereo", 1_i32)?;
                    }
                    if matches!(cfg, 1 | 3) && can_have_direction(order) {
                        record
                            .properties
                            .set_prop("_needsDetectAtomStereo", 1_i32)?;
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
                    let bond_id = match bond_with_smiles_index(&record.topology, *index)? {
                        Some(id) => id,
                        None => {
                            eprintln!(
                                "bond {} not found, cannot mark as stereo double bond.",
                                index
                            );
                            return Err(cx_failure());
                        }
                    };
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
                // BEGIN COMPLETE PINNED SF193
                // RDKit✔️❌: bool parse_linknodes(Iterator &first, Iterator last, RDKit::RWMol &mol,
                // RDKit✔️❌:                      unsigned int startAtomIdx) {
                // RDKit✔️❌:   // these look like: |LN:1:1.3.2.6,4:1.4.3.6|
                // RDKit✔️❌:   // that's two records:
                // RDKit✔️❌:   //   1:1.3.2.6: 1-3 repeats, atom 1-2, 1-6
                // RDKit✔️❌:   //   4:1.4.3.6: 1-4 repeats, atom 4-3, 4-6
                // RDKit✔️❌:   // which maps to the property value "1 3 2 2 3 2 7|1 4 2 5 4 5 7"
                // RDKit✔️❌:   // If the linking atom only has two neighbors then the outer atom
                // RDKit✔️❌:   // specification (the last two digits) can be left out. So for a molecule
                // RDKit✔️❌:   // where atom 1 has bonds only to atoms 2 and 6 we could have
                // RDKit✔️❌:   // |LN:1:1.3|
                // RDKit✔️❌:   // instead of
                // RDKit✔️❌:   // |LN:1:1.3.2.6|
                // RDKit✔️❌:   if (first >= last || *first != 'L' || first + 1 >= last ||
                // RDKit✔️❌:       *(first + 1) != 'N' || first + 2 >= last || *(first + 2) != ':') {
                // RDKit✔️❌:     return false;
                // RDKit✔️❌:   }
                // RDKit✔️❌:   first += 3;
                // RDKit✔️❌:   std::string accum = "";
                // RDKit✔️❌:   while (first < last && *first >= '0' && *first <= '9') {
                // RDKit✔️❌:     unsigned int atidx;
                // RDKit✔️❌:     if (!read_int(first, last, atidx)) {
                // RDKit✔️❌:       return false;
                // RDKit✔️❌:     }
                // RDKit✔️❌:     // check that we can read at least two more characters:
                // RDKit✔️❌:     if (first + 1 >= last || *first != ':') {
                // RDKit✔️❌:       return false;
                // RDKit✔️❌:     }
                // RDKit✔️❌:     ++first;
                // RDKit✔️❌:     unsigned int startReps;
                // RDKit✔️❌:     if (!read_int(first, last, startReps)) {
                // RDKit✔️❌:       return false;
                // RDKit✔️❌:     }
                // RDKit✔️❌:     if (first + 1 >= last || *first != '.') {
                // RDKit✔️❌:       return false;
                // RDKit✔️❌:     }
                // RDKit✔️❌:     ++first;
                // RDKit✔️❌:     unsigned int endReps;
                // RDKit✔️❌:     if (!read_int(first, last, endReps)) {
                // RDKit✔️❌:       return false;
                // RDKit✔️❌:     }
                // RDKit✔️❌:     unsigned int idx1;
                // RDKit✔️❌:     unsigned int idx2;
                // RDKit✔️❌:     if (first < last && *first == '.') {
                // RDKit✔️❌:       ++first;
                // RDKit✔️❌:       if (!read_int(first, last, idx1)) {
                // RDKit✔️❌:         return false;
                // RDKit✔️❌:       }
                // RDKit✔️❌:       ++first;
                // RDKit✔️❌:       if (!read_int(first, last, idx2)) {
                // RDKit✔️❌:         return false;
                // RDKit✔️❌:       }
                // RDKit✔️❌:     } else if (VALID_ATIDX(atidx) &&
                // RDKit✔️❌:                mol.getAtomWithIdx(atidx - startAtomIdx)->getDegree() == 2) {
                // RDKit✔️❌:       auto nbrs =
                // RDKit✔️❌:           mol.getAtomNeighbors(mol.getAtomWithIdx(atidx - startAtomIdx));
                // RDKit✔️❌:       idx1 = *nbrs.first;
                // RDKit✔️❌:       nbrs.first++;
                // RDKit✔️❌:       idx2 = *nbrs.first;
                // RDKit✔️❌:     } else if (VALID_ATIDX(atidx)) {
                // RDKit✔️❌:       return false;
                // RDKit✔️❌:     }
                // RDKit✔️❌:     if (first < last && *first == ',') {
                // RDKit✔️❌:       ++first;
                // RDKit✔️❌:     }
                // RDKit✔️❌:     if (VALID_ATIDX(atidx)) {
                // RDKit✔️❌:       if (!accum.empty()) {
                // RDKit✔️❌:         accum += "|";
                // RDKit✔️❌:       }
                // RDKit✔️❌:       accum += (boost::format("%d %d 2 %d %d %d %d") % startReps % endReps %
                // RDKit✔️❌:                 (atidx - startAtomIdx + 1) % (idx1 - startAtomIdx + 1) %
                // RDKit✔️❌:                 (atidx - startAtomIdx + 1) % (idx2 - startAtomIdx + 1))
                // RDKit✔️❌:                    .str();
                // RDKit✔️❌:     }
                // RDKit✔️❌:   }
                // RDKit✔️❌:   if (!accum.empty()) {
                // RDKit✔️❌:     mol.setProp(common_properties::molFileLinkNodes, accum);
                // RDKit✔️❌:   }
                // RDKit✔️❌:   return true;
                // RDKit✔️❌: }
                // END COMPLETE PINNED SF193
                // BEGIN COMPLETE source common_properties::molFileLinkNodes
                // RDKit✔️✔️: inline constexpr std::string_view molFileLinkNodes = "_molLinkNodes";
                // END COMPLETE source common_properties::molFileLinkNodes
                // BEGIN COMPLETE Boost1.81::parse_printf_directive
                // Boost✔️✔️:     bool parse_printf_directive(Iter & start, const Iter& last,
                // Boost✔️✔️:                                 detail::format_item<Ch, Tr, Alloc> * fpar,
                // Boost✔️✔️:                                 const Facet& fac,
                // Boost✔️✔️:                                 std::size_t offset, unsigned char exceptions)
                // Boost✔️✔️:     {
                // Boost✔️✔️:         typedef typename basic_format<Ch, Tr, Alloc>::format_item_t format_item_t;
                // Boost✔️✔️:
                // Boost✔️✔️:         fpar->argN_ = format_item_t::argN_no_posit;  // if no positional-directive
                // Boost✔️✔️:         bool precision_set = false;
                // Boost✔️✔️:         bool in_brackets=false;
                // Boost✔️✔️:         Iter start0 = start;
                // Boost✔️✔️:         std::size_t fstring_size = last-start0+offset;
                // Boost✔️✔️:         char mssiz = 0;
                // Boost✔️✔️:
                // Boost✔️✔️:         if(start>= last) { // empty directive : this is a trailing %
                // Boost✔️✔️:                 maybe_throw_exception(exceptions, start-start0 + offset, fstring_size);
                // Boost✔️✔️:                 return false;
                // Boost✔️✔️:         }
                // Boost✔️✔️:
                // Boost✔️✔️:         if(*start== const_or_not(fac).widen( '|')) {
                // Boost✔️✔️:             in_brackets=true;
                // Boost✔️✔️:             if( ++start >= last ) {
                // Boost✔️✔️:                 maybe_throw_exception(exceptions, start-start0 + offset, fstring_size);
                // Boost✔️✔️:                 return false;
                // Boost✔️✔️:             }
                // Boost✔️✔️:         }
                // Boost✔️✔️:
                // Boost✔️✔️:         // the flag '0' would be picked as a digit for argument order, but here it's a flag :
                // Boost✔️✔️:         if(*start== const_or_not(fac).widen( '0'))
                // Boost✔️✔️:             goto parse_flags;
                // Boost✔️✔️:
                // Boost✔️✔️:         // handle argument order (%2$d)  or possibly width specification: %2d
                // Boost✔️✔️:         if(wrap_isdigit(fac, *start)) {
                // Boost✔️✔️:             int n;
                // Boost✔️✔️:             start = str2int(start, last, n, fac);
                // Boost✔️✔️:             if( start >= last ) {
                // Boost✔️✔️:                 maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
                // Boost✔️✔️:                 return false;
                // Boost✔️✔️:             }
                // Boost✔️✔️:
                // Boost✔️✔️:             // %N% case : this is already the end of the directive
                // Boost✔️✔️:             if( *start ==  const_or_not(fac).widen( '%') ) {
                // Boost✔️✔️:                 fpar->argN_ = n-1;
                // Boost✔️✔️:                 ++start;
                // Boost✔️✔️:                 if( in_brackets)
                // Boost✔️✔️:                     maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
                // Boost✔️✔️:                 return true;
                // Boost✔️✔️:             }
                // Boost✔️✔️:
                // Boost✔️✔️:             if ( *start== const_or_not(fac).widen( '$') ) {
                // Boost✔️✔️:                 fpar->argN_ = n-1;
                // Boost✔️✔️:                 ++start;
                // Boost✔️✔️:             }
                // Boost✔️✔️:             else {
                // Boost✔️✔️:                 // non-positional directive
                // Boost✔️✔️:                 fpar->fmtstate_.width_ = n;
                // Boost✔️✔️:                 fpar->argN_  = format_item_t::argN_no_posit;
                // Boost✔️✔️:                 goto parse_precision;
                // Boost✔️✔️:             }
                // Boost✔️✔️:         }
                // Boost✔️✔️:
                // Boost✔️✔️:       parse_flags:
                // Boost✔️✔️:         // handle flags
                // Boost✔️✔️:         while (start != last) { // as long as char is one of + - = _ # 0 or ' '
                // Boost✔️✔️:             switch ( wrap_narrow(fac, *start, 0)) {
                // Boost✔️✔️:                 case '\'':
                // Boost✔️✔️:                     break; // no effect yet. (painful to implement)
                // Boost✔️✔️:                 case '-':
                // Boost✔️✔️:                     fpar->fmtstate_.flags_ |= std::ios_base::left;
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case '=':
                // Boost✔️✔️:                     fpar->pad_scheme_ |= format_item_t::centered;
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case '_':
                // Boost✔️✔️:                     fpar->fmtstate_.flags_ |= std::ios_base::internal;
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case ' ':
                // Boost✔️✔️:                     fpar->pad_scheme_ |= format_item_t::spacepad;
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case '+':
                // Boost✔️✔️:                     fpar->fmtstate_.flags_ |= std::ios_base::showpos;
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case '0':
                // Boost✔️✔️:                     fpar->pad_scheme_ |= format_item_t::zeropad;
                // Boost✔️✔️:                     // need to know alignment before really setting flags,
                // Boost✔️✔️:                     // so just add 'zeropad' flag for now, it will be processed later.
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case '#':
                // Boost✔️✔️:                     fpar->fmtstate_.flags_ |= std::ios_base::showpoint | std::ios_base::showbase;
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 default:
                // Boost✔️✔️:                     goto parse_width;
                // Boost✔️✔️:             }
                // Boost✔️✔️:             ++start;
                // Boost✔️✔️:         } // loop on flag.
                // Boost✔️✔️:
                // Boost✔️✔️:         if( start>=last) {
                // Boost✔️✔️:             maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
                // Boost✔️✔️:             return true;
                // Boost✔️✔️:         }
                // Boost✔️✔️:
                // Boost✔️✔️:       // first skip 'asterisk fields' : * or num (length)
                // Boost✔️✔️:       parse_width:
                // Boost✔️✔️:         if(*start == const_or_not(fac).widen( '*') )
                // Boost✔️✔️:             ++start;
                // Boost✔️✔️:         else if(start!=last && wrap_isdigit(fac, *start))
                // Boost✔️✔️:             start = str2int(start, last, fpar->fmtstate_.width_, fac);
                // Boost✔️✔️:
                // Boost✔️✔️:       parse_precision:
                // Boost✔️✔️:         if( start>= last) {
                // Boost✔️✔️:             maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
                // Boost✔️✔️:             return true;
                // Boost✔️✔️:         }
                // Boost✔️✔️:         // handle precision spec
                // Boost✔️✔️:         if (*start== const_or_not(fac).widen( '.')) {
                // Boost✔️✔️:             ++start;
                // Boost✔️✔️:             if(start != last && *start == const_or_not(fac).widen( '*') )
                // Boost✔️✔️:                 ++start;
                // Boost✔️✔️:             else if(start != last && wrap_isdigit(fac, *start)) {
                // Boost✔️✔️:                 start = str2int(start, last, fpar->fmtstate_.precision_, fac);
                // Boost✔️✔️:                 precision_set = true;
                // Boost✔️✔️:             }
                // Boost✔️✔️:             else
                // Boost✔️✔️:                 fpar->fmtstate_.precision_ =0;
                // Boost✔️✔️:         }
                // Boost✔️✔️:
                // Boost✔️✔️:       // argument type modifiers
                // Boost✔️✔️:         while (start != last) {
                // Boost✔️✔️:             switch (wrap_narrow(fac, *start, 0)) {
                // Boost✔️✔️:                 case 'h':
                // Boost✔️✔️:                 case 'l':
                // Boost✔️✔️:                 case 'j':
                // Boost✔️✔️:                 case 'z':
                // Boost✔️✔️:                 case 'L':
                // Boost✔️✔️:                     // boost::format ignores argument type modifiers as it relies on
                // Boost✔️✔️:                     // the type of the argument fed into it by operator %
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:
                // Boost✔️✔️:                 // Note that the ptrdiff_t argument type 't' from C++11 is not honored
                // Boost✔️✔️:                 // because it was already in use as the tabulation specifier in boost::format
                // Boost✔️✔️:                 // case 't':
                // Boost✔️✔️:
                // Boost✔️✔️:                 // Microsoft extensions:
                // Boost✔️✔️:                 // https://msdn.microsoft.com/en-us/library/tcxf1dw6.aspx
                // Boost✔️✔️:
                // Boost✔️✔️:                 case 'w':
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case 'I':
                // Boost✔️✔️:                     mssiz = 'I';
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case '3':
                // Boost✔️✔️:                     if (mssiz != 'I') {
                // Boost✔️✔️:                         maybe_throw_exception(exceptions, start - start0 + offset, fstring_size);
                // Boost✔️✔️:                         return true;
                // Boost✔️✔️:                     }
                // Boost✔️✔️:                     mssiz = '3';
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case '2':
                // Boost✔️✔️:                     if (mssiz != '3') {
                // Boost✔️✔️:                         maybe_throw_exception(exceptions, start - start0 + offset, fstring_size);
                // Boost✔️✔️:                         return true;
                // Boost✔️✔️:                     }
                // Boost✔️✔️:                     mssiz = 0x00;
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case '6':
                // Boost✔️✔️:                     if (mssiz != 'I') {
                // Boost✔️✔️:                         maybe_throw_exception(exceptions, start - start0 + offset, fstring_size);
                // Boost✔️✔️:                         return true;
                // Boost✔️✔️:                     }
                // Boost✔️✔️:                     mssiz = '6';
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 case '4':
                // Boost✔️✔️:                     if (mssiz != '6') {
                // Boost✔️✔️:                         maybe_throw_exception(exceptions, start - start0 + offset, fstring_size);
                // Boost✔️✔️:                         return true;
                // Boost✔️✔️:                     }
                // Boost✔️✔️:                     mssiz = 0x00;
                // Boost✔️✔️:                     break;
                // Boost✔️✔️:                 default:
                // Boost✔️✔️:                     if (mssiz && mssiz == 'I') {
                // Boost✔️✔️:                         mssiz = 0;
                // Boost✔️✔️:                     }
                // Boost✔️✔️:                     goto parse_conversion_specification;
                // Boost✔️✔️:             }
                // Boost✔️✔️:             ++start;
                // Boost✔️✔️:         } // loop on argument type modifiers to pick up 'hh', 'll', and the more complex microsoft ones
                // Boost✔️✔️:
                // Boost✔️✔️:       parse_conversion_specification:
                // Boost✔️✔️:         if (start >= last || mssiz) {
                // Boost✔️✔️:             maybe_throw_exception(exceptions, start - start0 + offset, fstring_size);
                // Boost✔️✔️:             return true;
                // Boost✔️✔️:         }
                // Boost✔️✔️:
                // Boost✔️✔️:         if( in_brackets && *start== const_or_not(fac).widen( '|') ) {
                // Boost✔️✔️:             ++start;
                // Boost✔️✔️:             return true;
                // Boost✔️✔️:         }
                // Boost✔️✔️:
                // Boost✔️✔️:         // The default flags are "dec" and "skipws"
                // Boost✔️✔️:         // so if changing the base, need to unset basefield first
                // Boost✔️✔️:
                // Boost✔️✔️:         switch (wrap_narrow(fac, *start, 0))
                // Boost✔️✔️:         {
                // Boost✔️✔️:             // Boolean
                // Boost✔️✔️:             case 'b':
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::boolalpha;
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:
                // Boost✔️✔️:             // Decimal
                // Boost✔️✔️:             case 'u':
                // Boost✔️✔️:             case 'd':
                // Boost✔️✔️:             case 'i':
                // Boost✔️✔️:                 // Defaults are sufficient
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:
                // Boost✔️✔️:             // Hex
                // Boost✔️✔️:             case 'X':
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
                // Boost✔️✔️:                 BOOST_FALLTHROUGH;
                // Boost✔️✔️:             case 'x':
                // Boost✔️✔️:             case 'p': // pointer => set hex.
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::hex;
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:
                // Boost✔️✔️:             // Octal
                // Boost✔️✔️:             case 'o':
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::oct;
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:
                // Boost✔️✔️:             // Floating
                // Boost✔️✔️:             case 'A':
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
                // Boost✔️✔️:                 BOOST_FALLTHROUGH;
                // Boost✔️✔️:             case 'a':
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::fixed;
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::scientific;
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:             case 'E':
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
                // Boost✔️✔️:                 BOOST_FALLTHROUGH;
                // Boost✔️✔️:             case 'e':
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::scientific;
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:             case 'F':
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
                // Boost✔️✔️:                 BOOST_FALLTHROUGH;
                // Boost✔️✔️:             case 'f':
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::fixed;
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:             case 'G':
                // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
                // Boost✔️✔️:                 BOOST_FALLTHROUGH;
                // Boost✔️✔️:             case 'g':
                // Boost✔️✔️:                 // default flags are correct here
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:
                // Boost✔️✔️:             // Tabulation (a boost::format extension)
                // Boost✔️✔️:             case 'T':
                // Boost✔️✔️:                 ++start;
                // Boost✔️✔️:                 if( start >= last) {
                // Boost✔️✔️:                     maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
                // Boost✔️✔️:                     return false;
                // Boost✔️✔️:                 } else {
                // Boost✔️✔️:                     fpar->fmtstate_.fill_ = *start;
                // Boost✔️✔️:                 }
                // Boost✔️✔️:                 fpar->pad_scheme_ |= format_item_t::tabulation;
                // Boost✔️✔️:                 fpar->argN_ = format_item_t::argN_tabulation;
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:             case 't':
                // Boost✔️✔️:                 fpar->fmtstate_.fill_ = const_or_not(fac).widen( ' ');
                // Boost✔️✔️:                 fpar->pad_scheme_ |= format_item_t::tabulation;
                // Boost✔️✔️:                 fpar->argN_ = format_item_t::argN_tabulation;
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:
                // Boost✔️✔️:             // Character
                // Boost✔️✔️:             case 'C':
                // Boost✔️✔️:             case 'c':
                // Boost✔️✔️:                 fpar->truncate_ = 1;
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:
                // Boost✔️✔️:             // String
                // Boost✔️✔️:             case 'S':
                // Boost✔️✔️:             case 's':
                // Boost✔️✔️:                 if(precision_set) // handle truncation manually, with own parameter.
                // Boost✔️✔️:                     fpar->truncate_ = fpar->fmtstate_.precision_;
                // Boost✔️✔️:                 fpar->fmtstate_.precision_ = 6; // default stream precision.
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:
                // Boost✔️✔️:             // %n is insecure and ignored by boost::format
                // Boost✔️✔️:             case 'n' :
                // Boost✔️✔️:                 fpar->argN_ = format_item_t::argN_ignored;
                // Boost✔️✔️:                 break;
                // Boost✔️✔️:
                // Boost✔️✔️:             default:
                // Boost✔️✔️:                 maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
                // Boost✔️✔️:         }
                // Boost✔️✔️:         ++start;
                // Boost✔️✔️:
                // Boost✔️✔️:         if( in_brackets ) {
                // Boost✔️✔️:             if( start != last && *start== const_or_not(fac).widen( '|') ) {
                // Boost✔️✔️:                 ++start;
                // Boost✔️✔️:                 return true;
                // Boost✔️✔️:             }
                // Boost✔️✔️:             else  maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
                // Boost✔️✔️:         }
                // Boost✔️✔️:         return true;
                // Boost✔️✔️:     }
                // END COMPLETE Boost1.81::parse_printf_directive
                // BEGIN COMPLETE Boost1.81::put
                // Boost✔️✔️:     void put( T x,
                // Boost✔️✔️:               const format_item<Ch, Tr, Alloc>& specs,
                // Boost✔️✔️:               typename basic_format<Ch, Tr, Alloc>::string_type& res,
                // Boost✔️✔️:               typename basic_format<Ch, Tr, Alloc>::internal_streambuf_t & buf,
                // Boost✔️✔️:               io::detail::locale_t *loc_p = NULL)
                // Boost✔️✔️:     {
                // Boost✔️✔️: #ifdef BOOST_MSVC
                // Boost✔️✔️:        // If std::min<unsigned> or std::max<unsigned> are already instantiated
                // Boost✔️✔️:        // at this point then we get a blizzard of warning messages when we call
                // Boost✔️✔️:        // those templates with std::size_t as arguments.  Weird and very annoyning...
                // Boost✔️✔️: #pragma warning(push)
                // Boost✔️✔️: #pragma warning(disable:4267)
                // Boost✔️✔️: #endif
                // Boost✔️✔️:         // does the actual conversion of x, with given params, into a string
                // Boost✔️✔️:         // using the supplied stringbuf.
                // Boost✔️✔️:
                // Boost✔️✔️:         typedef typename basic_format<Ch, Tr, Alloc>::string_type   string_type;
                // Boost✔️✔️:         typedef typename basic_format<Ch, Tr, Alloc>::format_item_t format_item_t;
                // Boost✔️✔️:         typedef typename string_type::size_type size_type;
                // Boost✔️✔️:
                // Boost✔️✔️:         basic_oaltstringstream<Ch, Tr, Alloc>  oss( &buf);
                // Boost✔️✔️:
                // Boost✔️✔️: #if !defined(BOOST_NO_STD_LOCALE)
                // Boost✔️✔️:         if(loc_p != NULL)
                // Boost✔️✔️:             oss.imbue(*loc_p);
                // Boost✔️✔️: #endif
                // Boost✔️✔️:
                // Boost✔️✔️:         specs.fmtstate_.apply_on(oss, loc_p);
                // Boost✔️✔️:
                // Boost✔️✔️:         // the stream format state can be modified by manipulators in the argument :
                // Boost✔️✔️:         put_head( oss, x );
                // Boost✔️✔️:         // in case x is a group, apply the manip part of it,
                // Boost✔️✔️:         // in order to find width
                // Boost✔️✔️:
                // Boost✔️✔️:         const std::ios_base::fmtflags fl=oss.flags();
                // Boost✔️✔️:         const bool internal = (fl & std::ios_base::internal) != 0;
                // Boost✔️✔️:         const std::streamsize w = oss.width();
                // Boost✔️✔️:         const bool two_stepped_padding= internal && (w!=0);
                // Boost✔️✔️:
                // Boost✔️✔️:         res.resize(0);
                // Boost✔️✔️:         if(! two_stepped_padding) {
                // Boost✔️✔️:             if(w>0) // handle padding via mk_str, not natively in stream
                // Boost✔️✔️:                 oss.width(0);
                // Boost✔️✔️:             put_last( oss, x);
                // Boost✔️✔️:             const Ch * res_beg = buf.pbase();
                // Boost✔️✔️:             Ch prefix_space = 0;
                // Boost✔️✔️:             if(specs.pad_scheme_ & format_item_t::spacepad)
                // Boost✔️✔️:                 if(buf.pcount()== 0 ||
                // Boost✔️✔️:                    (res_beg[0] !=oss.widen('+') && res_beg[0] !=oss.widen('-')  ))
                // Boost✔️✔️:                     prefix_space = oss.widen(' ');
                // Boost✔️✔️:             size_type res_size = (std::min)(
                // Boost✔️✔️:                 (static_cast<size_type>((specs.truncate_ & (std::numeric_limits<size_type>::max)())) - !!prefix_space),
                // Boost✔️✔️:                 buf.pcount() );
                // Boost✔️✔️:             mk_str(res, res_beg, res_size, w, oss.fill(), fl,
                // Boost✔️✔️:                    prefix_space, (specs.pad_scheme_ & format_item_t::centered) !=0 );
                // Boost✔️✔️:         }
                // Boost✔️✔️:         else  { // 2-stepped padding
                // Boost✔️✔️:             // internal can be implied by zeropad, or user-set.
                // Boost✔️✔️:             // left, right, and centered alignment overrule internal,
                // Boost✔️✔️:             // but spacepad or truncate might be mixed with internal (using manipulator)
                // Boost✔️✔️:             put_last( oss, x); // may pad
                // Boost✔️✔️:             const Ch * res_beg = buf.pbase();
                // Boost✔️✔️:             size_type res_size = buf.pcount();
                // Boost✔️✔️:             bool prefix_space=false;
                // Boost✔️✔️:             if(specs.pad_scheme_ & format_item_t::spacepad)
                // Boost✔️✔️:                 if(buf.pcount()== 0 ||
                // Boost✔️✔️:                    (res_beg[0] !=oss.widen('+') && res_beg[0] !=oss.widen('-')  ))
                // Boost✔️✔️:                     prefix_space = true;
                // Boost✔️✔️:             if(res_size == static_cast<size_type>(w) && w<=specs.truncate_ && !prefix_space) {
                // Boost✔️✔️:                 // okay, only one thing was printed and padded, so res is fine
                // Boost✔️✔️:                 res.assign(res_beg, res_size);
                // Boost✔️✔️:             }
                // Boost✔️✔️:             else { //   length w exceeded
                // Boost✔️✔️:                 // either it was multi-output with first output padding up all width..
                // Boost✔️✔️:                 // either it was one big arg and we are fine.
                // Boost✔️✔️:                 // Note that res_size<w is possible  (in case of bad user-defined formatting)
                // Boost✔️✔️:                 res.assign(res_beg, res_size);
                // Boost✔️✔️:                 res_beg=NULL;  // invalidate pointers.
                // Boost✔️✔️:
                // Boost✔️✔️:                 // make a new stream, to start re-formatting from scratch :
                // Boost✔️✔️:                 buf.clear_buffer();
                // Boost✔️✔️:                 basic_oaltstringstream<Ch, Tr, Alloc>  oss2( &buf);
                // Boost✔️✔️:                 specs.fmtstate_.apply_on(oss2, loc_p);
                // Boost✔️✔️:                 put_head( oss2, x );
                // Boost✔️✔️:
                // Boost✔️✔️:                 oss2.width(0);
                // Boost✔️✔️:                 if(prefix_space)
                // Boost✔️✔️:                     oss2 << ' ';
                // Boost✔️✔️:                 put_last(oss2, x );
                // Boost✔️✔️:                 if(buf.pcount()==0 && specs.pad_scheme_ & format_item_t::spacepad) {
                // Boost✔️✔️:                     prefix_space =true;
                // Boost✔️✔️:                     oss2 << ' ';
                // Boost✔️✔️:                 }
                // Boost✔️✔️:                 // we now have the minimal-length output
                // Boost✔️✔️:                 const Ch * tmp_beg = buf.pbase();
                // Boost✔️✔️:                 size_type tmp_size = (std::min)(
                // Boost✔️✔️:                     (static_cast<size_type>(specs.truncate_ & (std::numeric_limits<size_type>::max)())),
                // Boost✔️✔️:                     buf.pcount());
                // Boost✔️✔️:
                // Boost✔️✔️:                 if(static_cast<size_type>(w) <= tmp_size) {
                // Boost✔️✔️:                     // minimal length is already >= w, so no padding (cool!)
                // Boost✔️✔️:                         res.assign(tmp_beg, tmp_size);
                // Boost✔️✔️:                 }
                // Boost✔️✔️:                 else { // hum..  we need to pad (multi_output, or spacepad present)
                // Boost✔️✔️:                     //find where we should pad
                // Boost✔️✔️:                     size_type sz = (std::min)(res_size + (prefix_space ? 1 : 0), tmp_size);
                // Boost✔️✔️:                     size_type i = prefix_space;
                // Boost✔️✔️:                     for(; i<sz && tmp_beg[i] == res[i - (prefix_space ? 1 : 0)]; ++i) {}
                // Boost✔️✔️:                     if(i>=tmp_size) i=prefix_space;
                // Boost✔️✔️:                     res.assign(tmp_beg, i);
                // Boost✔️✔️:                                         std::streamsize d = w - static_cast<std::streamsize>(tmp_size);
                // Boost✔️✔️:                                         BOOST_ASSERT(d>0);
                // Boost✔️✔️:                     res.append(static_cast<size_type>( d ), oss2.fill());
                // Boost✔️✔️:                     res.append(tmp_beg+i, tmp_size-i);
                // Boost✔️✔️:                     BOOST_ASSERT(i+(tmp_size-i)+(std::max)(d,(std::streamsize)0)
                // Boost✔️✔️:                                  == static_cast<size_type>(w));
                // Boost✔️✔️:                     BOOST_ASSERT(res.size() == static_cast<size_type>(w));
                // Boost✔️✔️:                 }
                // Boost✔️✔️:             }
                // Boost✔️✔️:         }
                // Boost✔️✔️:         buf.clear_buffer();
                // Boost✔️✔️: #ifdef BOOST_MSVC
                // Boost✔️✔️: #pragma warning(pop)
                // Boost✔️✔️: #endif
                // Boost✔️✔️:     }
                // END COMPLETE Boost1.81::put
                // BEGIN COMPLETE Boost1.81::put_last<unsigned>
                // Boost✔️✔️:     void put_last( BOOST_IO_STD basic_ostream<Ch, Tr> & os, const T& x ) {
                // Boost✔️✔️:         os << x ;
                // Boost✔️✔️:     }
                // END COMPLETE Boost1.81::put_last<unsigned>
                // Same local window and source insertion-order neighbors;
                // explicit outer indices are not checked against graph size.
                // The native property key is _molLinkNodes, not a display name.
                // Formatting wraps each u32 +1 before decimal emission.
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
                    let source_uint = |value: usize| u32::try_from(value).map_err(|_| cx_failure());
                    let center_one = source_uint(link.atom)?.wrapping_add(1);
                    let outer_one = source_uint(outer_atoms[0])?.wrapping_add(1);
                    let outer_two = source_uint(outer_atoms[1])?.wrapping_add(1);
                    let start_repetitions = source_uint(link.start_repetitions)?;
                    let end_repetitions = source_uint(link.end_repetitions)?;
                    lowered.push(format!(
                        "{start_repetitions} {end_repetitions} 2 {center_one} {outer_one} {center_one} {outer_two}"
                    ));
                }
                if !lowered.is_empty() {
                    record
                        .properties
                        .set_prop("_molLinkNodes", lowered.join("|"))?;
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
                        field_display: Some("    0.0000    0.0000    DR    ALL  0       0".into()),
                        query_op: (!data.query_op.is_empty()).then(|| data.query_op.clone()),
                        values: (!data.data.is_empty())
                            .then(|| vec![data.data.clone()])
                            .unwrap_or_default(),
                        ..SGroupData::default()
                    };
                    let mut group = SubstanceGroup::new(group_id, SubstanceGroupKind::Data)
                        .with_atoms(atoms)
                        .with_data(typed_data);
                    // BEGIN COMPLETE PINNED SF194 property-field writes
                    // RDKit✔️❌: void parse_data_sgroup_attr(Iterator &first, Iterator last,
                    // RDKit✔️❌:                             SubstanceGroup &sgroup, bool keepSGroup,
                    // RDKit✔️❌:                             std::string fieldName, bool fieldIsArray = false) {
                    // RDKit✔️❌:   PRECONDITION(first < last, "parse_data_sgroup_attr: first >= last");
                    // RDKit✔️❌:   if (first != last && *first != '|') {
                    // RDKit✔️❌:     std::string data = read_text_to(first, last, ":");
                    // RDKit✔️❌:     ++first;
                    // RDKit✔️❌:     if (!data.empty() && keepSGroup) {
                    // RDKit✔️❌:       if (fieldIsArray) {
                    // RDKit✔️❌:         std::vector<std::string> dataFields = {data};
                    // RDKit✔️❌:         sgroup.setProp(fieldName, dataFields);
                    // RDKit✔️❌:       } else {
                    // RDKit✔️❌:         sgroup.setProp(fieldName, data);
                    // RDKit✔️❌:       }
                    // RDKit✔️❌:     }
                    // RDKit✔️❌:   }
                    // RDKit✔️❌: }
                    // END COMPLETE PINNED SF194 property-field writes
                    // BEGIN COMPLETE SubstanceGroup::SubstanceGroup(TYPE)
                    // RDKit✔️❌: SubstanceGroup::SubstanceGroup(ROMol *owning_mol, const std::string &type)
                    // RDKit✔️❌:     : RDProps(), dp_mol(owning_mol) {
                    // RDKit✔️❌:   PRECONDITION(owning_mol, "supplied owning molecule is bad");
                    // RDKit✔️❌:
                    // RDKit✔️❌:   // TYPE is required to be set , as other properties will depend on it.
                    // RDKit✔️❌:   setProp<std::string>("TYPE", type);
                    // RDKit✔️❌: }
                    // END COMPLETE SubstanceGroup::SubstanceGroup(TYPE)
                    // Source property insertion order is retained by the sole store: TYPE,
                    // sequence index, FIELDNAME when nonempty, FIELDDISP, remaining fields,
                    // optional COORDS, then dense source index at helper completion.
                    group.set_prop("TYPE", "DAT")?;
                    group.set_prop("_cxsmilesindex", sgroup_index)?;
                    if !data.field_name.is_empty() {
                        group.set_prop("FIELDNAME", data.field_name.clone())?;
                    }
                    group.set_prop("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0")?;
                    if !data.data.is_empty() {
                        group.set_prop("DATAFIELDS", vec![data.data.clone()])?;
                        group.push_data_field(data.data.clone());
                    }
                    if !data.query_op.is_empty() {
                        group.set_prop("QUERYOP", data.query_op.clone())?;
                    }
                    if !data.field_info.is_empty() {
                        group.set_prop("FIELDINFO", data.field_info.clone())?;
                    }
                    if !data.field_tag.is_empty() {
                        group.set_prop("FIELDTAG", data.field_tag.clone())?;
                    }
                    if let Some(coordinates) = &data.coordinates {
                        group.set_prop("COORDS", coordinates.clone())?;
                    }
                    group.set_prop(
                        "index",
                        (record.topology.substance_groups.len() as u32).wrapping_add(1),
                    )?;
                    record.topology.substance_groups.push(group);
                }
                sgroup_index = sgroup_index.wrapping_add(1);
            }
            CxRecord::SGroupHierarchy(relationships) => {
                // RDKit❗❌: bool parse_sgroup_hierarchy(Iterator &first, Iterator last, RDKit::RWMol &mol) {
                // RDKit❗❌:   // these look like: |SgH:1:0|
                // RDKit❗❌:   // from CXSMILES docs:
                // RDKit❗❌:   //    SgH:parentSgroupIndex1:childSgroupIndex1.childSgroupIndex2,parentSgroupIndex2:childSgroupIndex1
                // RDKit❗❌:   if (first >= last || *first != 'S' || first + 3 >= last ||
                // RDKit❗❌:       *(first + 1) != 'g' || *(first + 2) != 'H' || *(first + 3) != ':') {
                // RDKit❗❌:     return false;
                // RDKit❗❌:   }
                // RDKit❗❌:   first += 4;
                // RDKit❗❌:   auto &sgs = getSubstanceGroups(mol);
                // RDKit❗❌:   while (1) {
                // RDKit❗❌:     unsigned int parentId;
                // RDKit❗❌:     if (!read_int(first, last, parentId)) {
                // RDKit❗❌:       return false;
                // RDKit❗❌:     }
                // RDKit❗❌:
                // RDKit❗❌:     bool validParent = true;
                // RDKit❗❌:     auto psg = find_matching_sgroup(sgs, parentId);
                // RDKit❗❌:     if (psg == sgs.end()) {
                // RDKit❗❌:       validParent = false;
                // RDKit❗❌:     } else {
                // RDKit❗❌:       psg->getPropIfPresent("index", parentId);
                // RDKit❗❌:     }
                // RDKit❗❌:     if (first <= last && *first == ':') {
                // RDKit❗❌:       ++first;
                // RDKit❗❌:       std::vector<unsigned int> children;
                // RDKit❗❌:       if (!read_int_list(first, last, children, '.')) {
                // RDKit❗❌:         return false;
                // RDKit❗❌:       }
                // RDKit❗❌:       if (validParent) {
                // RDKit❗❌:         for (auto childId : children) {
                // RDKit❗❌:           if (childId >= sgs.size()) {
                // RDKit❗❌:             throw SmilesParseException(
                // RDKit❗❌:                 "child id references non-existent SGroup");
                // RDKit❗❌:           }
                // RDKit❗❌:           auto csg = find_matching_sgroup(sgs, childId);
                // RDKit❗❌:           if (csg != sgs.end()) {
                // RDKit❗❌:             unsigned int cid;
                // RDKit❗❌:             csg->getProp("index", cid);
                // RDKit❗❌:             csg->setProp("PARENT", parentId);
                // RDKit❗❌:           }
                // RDKit❗❌:         }
                // RDKit❗❌:       }
                // RDKit❗❌:       if (first <= last && *first == ',') {
                // RDKit❗❌:         ++first;
                // RDKit❗❌:       } else {
                // RDKit❗❌:         break;
                // RDKit❗❌:       }
                // RDKit❗❌:     } else {
                // RDKit❗❌:       return false;
                // RDKit❗❌:     }
                // RDKit❗❌:   }
                // RDKit❗❌:
                // RDKit❗❌:   return true;
                // RDKit❗❌: }
                // RDKit❗❌: std::vector<RDKit::SubstanceGroup>::iterator find_matching_sgroup(
                // RDKit❗❌:     std::vector<RDKit::SubstanceGroup> &sgs, unsigned int targetId) {
                // RDKit❗❌:   return std::find_if(sgs.begin(), sgs.end(), [targetId](const auto &sg) {
                // RDKit❗❌:     unsigned int pval;
                // RDKit❗❌:     if (sg.getPropIfPresent(cxsmilesindex, pval)) {
                // RDKit❗❌:       if (pval == targetId) {
                // RDKit❗❌:         return true;
                // RDKit❗❌:       }
                // RDKit❗❌:     }
                // RDKit❗❌:     return false;
                // RDKit❗❌:   });
                // RDKit❗❌: }
                // Visit exactly the source parent/child lookups. Building an
                // eager index would cast unvisited later properties before a
                // first match or before the source child bound/missing-index
                // failure. The existing CORE conversion owns each reached cast.
                // Full parser/progress behavior remains under its own source
                // pair; this detached typed consumer preserves lookup order.
                // Cost: source-shaped linear group scans, logarithmic property
                // access in the sole MODEL store, no temporary index map.
                for relationship in relationships {
                    let mut parent_index =
                        u32::try_from(relationship.parent).map_err(model_failure)?;
                    let mut parent = None;
                    for group in &record.topology.substance_groups {
                        let Some(value) = group.props().get(b"_cxsmilesindex".as_slice()) else {
                            continue;
                        };
                        if cosmolkit_core::property_value_to_uint(value)
                            .map_err(SmilesParseError::WriterNumeric)?
                            == parent_index
                        {
                            if let Some(value) = group.props().get(b"index".as_slice()) {
                                parent_index = cosmolkit_core::property_value_to_uint(value)
                                    .map_err(SmilesParseError::WriterNumeric)?;
                            }
                            parent = Some(group.id());
                            break;
                        }
                    }
                    let Some(parent_id) = parent else {
                        continue;
                    };
                    for &child in &relationship.children {
                        if child >= record.topology.substance_groups.len() {
                            return Err(SmilesParseError::Cx(
                                "child id references non-existent SGroup".to_owned(),
                            ));
                        }
                        let child = u32::try_from(child).map_err(model_failure)?;
                        for group in &mut record.topology.substance_groups {
                            let Some(value) = group.props().get(b"_cxsmilesindex".as_slice())
                            else {
                                continue;
                            };
                            if cosmolkit_core::property_value_to_uint(value)
                                .map_err(SmilesParseError::WriterNumeric)?
                                != child
                            {
                                continue;
                            }
                            let value =
                                group.props().get(b"index".as_slice()).ok_or_else(|| {
                                    SmilesParseError::Cx(
                                        "SGroup child is missing its source index property"
                                            .to_owned(),
                                    )
                                })?;
                            cosmolkit_core::property_value_to_uint(value)
                                .map_err(SmilesParseError::WriterNumeric)?;
                            group.set_prop("PARENT", parent_index)?;
                            group.set_parent(parent_id);
                            break;
                        }
                    }
                }
            }
            CxRecord::PolymerSGroup(polymer) => {
                let (kind, source_type) =
                    sgroup_kind(polymer.type_code.as_bytes()).ok_or_else(cx_failure)?;
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
                    // Source SubstanceGroup constructor writes TYPE before the
                    // sequence/subtype fields, using the pinned typemap value.
                    // RDKit✔️❌: SubstanceGroup::SubstanceGroup(ROMol *owning_mol, const std::string &type)
                    // RDKit✔️❌:     : RDProps(), dp_mol(owning_mol) {
                    // RDKit✔️❌:   PRECONDITION(owning_mol, "supplied owning molecule is bad");
                    // RDKit✔️❌:
                    // RDKit✔️❌:   // TYPE is required to be set , as other properties will depend on it.
                    // RDKit✔️❌:   setProp<std::string>("TYPE", type);
                    // RDKit✔️❌: }
                    // Ordinary TYPE assignment uses the source String tag; the sole
                    // store preserves order, with its known extra key allocation cost.
                    group.set_prop("TYPE", source_type)?;
                    group.set_prop("_cxsmilesindex", sgroup_index)?;
                    match polymer.type_code.as_bytes() {
                        b"alt" => {
                            group.set_prop("SUBTYPE", "ALT")?;
                            group.set_subtype("ALT");
                        }
                        b"ran" => {
                            group.set_prop("SUBTYPE", "RAN")?;
                            group.set_subtype("RAN");
                        }
                        b"blk" => {
                            group.set_prop("SUBTYPE", "BLO")?;
                            group.set_subtype("BLO");
                        }
                        _ => {}
                    }
                    if !polymer.label.is_empty() {
                        group.set_prop("LABEL", polymer.label.clone())?;
                        group.set_label(polymer.label.clone());
                    }
                    if !polymer.connect.is_empty() {
                        group.set_prop("CONNECT", polymer.connect.clone())?;
                    }
                    let keep_group = finalize_polymer_sgroup(
                        &record.topology,
                        &mut group,
                        polymer.connect.as_bytes(),
                        &polymer.head_crossings,
                        &polymer.tail_crossings,
                    )?;
                    if keep_group {
                        // RDKit❗❌: bool parse_polymer_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
                        // RDKit❗❌:                           unsigned int startAtomIdx, unsigned int nSGroups) {
                        // RDKit❗❌:   // these look like:
                        // RDKit❗❌:   //    |Sg:n:6,1,2,4::hh&#44;f:6,0,:4,2,|
                        // RDKit❗❌:   // example from CXSMILES docs:
                        // RDKit❗❌:   // the fields are:
                        // RDKit❗❌:   //    Sg:[type]:[atom indices]:[subscript]:[superscript]:[head crossing
                        // RDKit❗❌:   //    bonds]:[tail crossing bonds]:
                        // RDKit❗❌:   //
                        // RDKit❗❌:   // note that it's legit for empty fields to be completely missing.
                        // RDKit❗❌:   //   for example, this doesn't have any crossing bonds indicated:
                        // RDKit❗❌:   // *-CCCN-* |$star_e;;;;;star_e$,Sg:n:4,1,2,3::hh|
                        // RDKit❗❌:   // this last bit makes the whole thing doubleplusfun to parse
                        // RDKit❗❌:
                        // RDKit❗❌:   if (first >= last || *first != 'S' || first + 2 >= last ||
                        // RDKit❗❌:       *(first + 1) != 'g' || *(first + 2) != ':') {
                        // RDKit❗❌:     return false;
                        // RDKit❗❌:   }
                        // RDKit❗❌:   first += 3;
                        // RDKit❗❌:
                        // RDKit❗❌:   const auto type_code = read_text_to(first, last, ":");
                        // RDKit❗❌:   ++first;
                        // RDKit❗❌:   const auto type = sgroupTypemap.find(type_code);
                        // RDKit❗❌:   if (type == sgroupTypemap.end()) {
                        // RDKit❗❌:     return false;
                        // RDKit❗❌:   }
                        // RDKit❗❌:   bool keepSGroup = false;
                        // RDKit❗❌:   SubstanceGroup sgroup(&mol, type->second);
                        // RDKit❗❌:   sgroup.setProp(cxsmilesindex, nSGroups);
                        // RDKit❗❌:   if (type_code == "alt") {
                        // RDKit❗❌:     sgroup.setProp("SUBTYPE", std::string("ALT"));
                        // RDKit❗❌:   } else if (type_code == "ran") {
                        // RDKit❗❌:     sgroup.setProp("SUBTYPE", std::string("RAN"));
                        // RDKit❗❌:   } else if (type_code == "blk") {
                        // RDKit❗❌:     sgroup.setProp("SUBTYPE", std::string("BLO"));
                        // RDKit❗❌:   }
                        // RDKit❗❌:
                        // RDKit❗❌:   std::vector<unsigned int> atoms;
                        // RDKit❗❌:   if (!read_int_list(first, last, atoms)) {
                        // RDKit❗❌:     return false;
                        // RDKit❗❌:   }
                        // RDKit❗❌:   //++first;
                        // RDKit❗❌:   for (auto idx : atoms) {
                        // RDKit❗❌:     if (VALID_ATIDX(idx)) {
                        // RDKit❗❌:       sgroup.addAtomWithIdx(idx - startAtomIdx);
                        // RDKit❗❌:       keepSGroup = true;
                        // RDKit❗❌:     }
                        // RDKit❗❌:   }
                        // RDKit❗❌:   std::vector<unsigned int> headCrossing;
                        // RDKit❗❌:   std::vector<unsigned int> tailCrossing;
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
                        // RDKit❗❌:
                        // RDKit❗❌:       if (first <= last && *first == ':') {
                        // RDKit❗❌:         ++first;
                        // RDKit❗❌:         if (!read_int_list(first, last, headCrossing)) {
                        // RDKit❗❌:           return false;
                        // RDKit❗❌:         }
                        // RDKit❗❌:         if (keepSGroup && !headCrossing.empty()) {
                        // RDKit❗❌:           for (auto &cidx : headCrossing) {
                        // RDKit❗❌:             if (VALID_ATIDX(cidx)) {
                        // RDKit❗❌:               cidx -= startAtomIdx;
                        // RDKit❗❌:             } else {
                        // RDKit❗❌:               keepSGroup = false;
                        // RDKit❗❌:               break;
                        // RDKit❗❌:             }
                        // RDKit❗❌:           }
                        // RDKit❗❌:           sgroup.setProp(_headCrossings, headCrossing, true);
                        // RDKit❗❌:         }
                        // RDKit❗❌:         if (first <= last && *first == ':') {
                        // RDKit❗❌:           ++first;
                        // RDKit❗❌:           if (!read_int_list(first, last, tailCrossing)) {
                        // RDKit❗❌:             return false;
                        // RDKit❗❌:           }
                        // RDKit❗❌:         }
                        // RDKit❗❌:         if (keepSGroup && !tailCrossing.empty()) {
                        // RDKit❗❌:           for (auto &cidx : tailCrossing) {
                        // RDKit❗❌:             if (VALID_ATIDX(cidx)) {
                        // RDKit❗❌:               cidx -= startAtomIdx;
                        // RDKit❗❌:             } else {
                        // RDKit❗❌:               keepSGroup = false;
                        // RDKit❗❌:               break;
                        // RDKit❗❌:             }
                        // RDKit❗❌:           }
                        // RDKit❗❌:           sgroup.setProp("_tailCrossings", tailCrossing, true);
                        // RDKit❗❌:         }
                        // RDKit❗❌:       }
                        // RDKit❗❌:     }
                        // RDKit❗❌:   }
                        // RDKit❗❌:   if (keepSGroup) {  // the label processing can destroy sgroup info, so do that
                        // RDKit❗❌:                      // now (the function will immediately return if already
                        // RDKit❗❌:                      // called)
                        // RDKit❗❌:     processCXSmilesLabels(mol);
                        // RDKit❗❌:
                        // RDKit❗❌:     finalizePolymerSGroup(mol, sgroup);
                        // RDKit❗❌:     sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
                        // RDKit❗❌:
                        // RDKit❗❌:     addSubstanceGroup(mol, sgroup);
                        // RDKit❗❌:   }
                        // RDKit❗❌:   return true;
                        // RDKit❗❌: }
                        // Native explicit unsigned construction happens after
                        // finalization. A failed/skipped local group never writes
                        // index; conversion keeps the u32 tag and wrapping sum.
                        group.set_prop(
                            "index",
                            (record.topology.substance_groups.len() as u32).wrapping_add(1),
                        )?;
                        record.topology.substance_groups.push(group);
                    }
                }
                sgroup_index = sgroup_index.wrapping_add(1);
            }
            CxRecord::VariableAttachments(attachments) => {
                // BEGIN COMPLETE PINNED SF199 concrete caller
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
                // END COMPLETE PINNED SF199 concrete caller
                for attachment in attachments {
                    let atom = AtomId::new(attachment.atom);
                    if atom.index() >= atom_count {
                        continue;
                    }
                    if atom_neighbors(&record.topology, atom).len() != 1 {
                        eprintln!("position variation bond to atom with more than one bond");
                        return Err(cx_failure());
                    }
                    let endpoints = attachment
                        .endpoints
                        .iter()
                        .filter(|index| **index < atom_count)
                        .map(|index| (*index as u32).wrapping_add(1).to_string())
                        .collect::<Vec<_>>();
                    let value = if endpoints.is_empty() {
                        "(0)".to_owned()
                    } else {
                        format!("({} {})", endpoints.len(), endpoints.join(" "))
                    };
                    let bond_ids = atom_neighbors(&record.topology, atom)
                        .iter()
                        .map(|neighbor| neighbor.bond)
                        .collect::<Vec<_>>();
                    for bond_id in bond_ids {
                        let bond = &mut record.topology.bonds[bond_id.index()];
                        bond.set_prop("_MolFileBondEndPts", value.clone())?;
                        bond.set_prop("_MolFileBondAttach", "ANY")?;
                    }
                }
            }
            CxRecord::Unsaturation(_) | CxRecord::RingBonds(_) | CxRecord::Substitution(_) => {
                unreachable!("query records were preflighted")
            }
            CxRecord::Unknown(_) => {}
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

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    fn graph(props: Vec<cosmolkit_model::PropertyValue>) -> cosmolkit_model::TopologyBlock {
        let atoms = (0..props.len() + 1)
            .map(|i| {
                cosmolkit_model::Atom::from_spec(
                    cosmolkit_model::AtomId::new(i),
                    cosmolkit_model::AtomSpec::new(cosmolkit_types::Element::C),
                )
            })
            .collect();
        let bonds = props
            .into_iter()
            .enumerate()
            .map(|(i, v)| {
                cosmolkit_model::Bond::from_spec(
                    cosmolkit_model::BondId::new(i),
                    cosmolkit_model::BondSpec::new(
                        cosmolkit_model::AtomId::new(i),
                        cosmolkit_model::AtomId::new(i + 1),
                        cosmolkit_types::BondOrder::Single,
                    )
                    .with_prop("_cxsmilesBondIdx", v)
                    .unwrap(),
                )
            })
            .collect();
        cosmolkit_model::TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }

    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXlower_0
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxlower_0_cx_lowering() {
        let g = graph(vec![cosmolkit_model::PropertyValue::UInt(0_u32)]);
        let before = g.clone();
        assert_eq!(
            bond_with_smiles_index(&g, 0_usize),
            Ok(Some(BondId::new(0)))
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXlower_1
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxlower_1_cx_lowering() {
        let g = graph(vec![cosmolkit_model::PropertyValue::UInt(1_u32)]);
        let before = g.clone();
        assert_eq!(
            bond_with_smiles_index(&g, 1_usize),
            Ok(Some(BondId::new(0)))
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXlower_2147483646
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxlower_2147483646_cx_lowering() {
        let g = graph(vec![cosmolkit_model::PropertyValue::UInt(2147483646_u32)]);
        let before = g.clone();
        assert_eq!(
            bond_with_smiles_index(&g, 2147483646_usize),
            Ok(Some(BondId::new(0)))
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXlower_2147483647
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxlower_2147483647_cx_lowering() {
        let g = graph(vec![cosmolkit_model::PropertyValue::UInt(2147483647_u32)]);
        let before = g.clone();
        assert_eq!(
            bond_with_smiles_index(&g, 2147483647_usize),
            Ok(Some(BondId::new(0)))
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXlower_2147483648
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxlower_2147483648_cx_lowering() {
        let g = graph(vec![cosmolkit_model::PropertyValue::UInt(2147483648_u32)]);
        let before = g.clone();
        assert_eq!(
            bond_with_smiles_index(&g, 2147483648_usize),
            Ok(Some(BondId::new(0)))
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXlower_4294967295
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxlower_4294967295_cx_lowering() {
        let g = graph(vec![cosmolkit_model::PropertyValue::UInt(4294967295_u32)]);
        let before = g.clone();
        assert_eq!(
            bond_with_smiles_index(&g, 4294967295_usize),
            Ok(Some(BondId::new(0)))
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: CX_FIRST_MATCH
    #[test]
    fn uint_cell_cx_first_match_cx_lowering() {
        let g = graph(vec![
            cosmolkit_model::PropertyValue::UInt(0),
            cosmolkit_model::PropertyValue::IntVector(vec![]),
        ]);
        let before = g.clone();
        assert_eq!(bond_with_smiles_index(&g, 0), Ok(Some(BondId::new(0))));
        assert_eq!(g, before);
    }
}

#[cfg(test)]
fn fixed_property_text(value: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(value.as_bytes()).expect("original fixed fixture text is UTF8")
}
