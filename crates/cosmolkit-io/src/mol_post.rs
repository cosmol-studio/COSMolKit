//! Molfile-specific detached postprocessing.

use cosmolkit_core::{
    AtropisomerConformer, RemoveHsParams, RingSearchParams, SanitizeOperations, SanitizeParams,
    StructureTagParams, ValenceModel, assign_chiral_tags_from_structure,
    assign_chiral_types_from_bond_dirs, assign_legacy_stereochemistry_with_query_state,
    assign_valence_for_topology, calculate_explicit_valence_for_topology,
    clear_single_bond_directions, detect_atropisomer_chirality, expand_attachment_points,
    remove_hydrogens_with_query_state, sanitize_topology_with_query_state,
    set_double_bond_neighbor_directions, symmetrized_sssr,
};
use cosmolkit_model::{
    AdjacencyList, AtomId, AtomQueryPredicate, BondQueryPredicate, Conformer3D, CoordinateBlock,
    PropertyText, PropertyValue, QueryAtom, QueryAtomConversionError, QueryBond, QueryGraph,
    QueryNode, QueryStateRef, RecursiveStructureQuery, SubstanceGroup, SubstanceGroupId,
    TopologyBlock, TopologyMapping, query_substance_groups, remap_query_rows,
    replace_query_substance_groups,
};
use cosmolkit_types::BondOrder;

use crate::sdf::{
    MolBlockRecord, QueryMolBlockRecord, parse_rdkit_int, parse_rdkit_unsigned,
    query_from_concrete_atom_value,
};

/// Source options applied after Molfile syntax parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolPostParams {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub expand_attachment_points: bool,
}

impl Default for MolPostParams {
    fn default() -> Self {
        // BEGIN RDKIT CPP FUNCTION MolFileParserParams
        // RDKit✔️✔️:   bool sanitize = true;      /**< sanitize the molecule after building it */
        // RDKit✔️✔️:   bool removeHs = true;      /**< remove Hs after constructing the molecule */
        // RDKit✔️✔️:   bool expandAttachmentPoints =
        // RDKit✔️✔️:       false; /**< toggle conversion of attachment points into dummy atoms */
        // END RDKIT CPP FUNCTION
        Self {
            sanitize: true,
            remove_hs: true,
            expand_attachment_points: false,
        }
    }
}

/// Typed reasons carried inside the existing Molfile postprocessing categories.
/// This value grants no runtime or commit authority.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MolProcessingError {
    #[error("{0}")]
    Valence(#[from] cosmolkit_core::ValenceError),
    #[error("{0}")]
    PropertyInt(#[from] cosmolkit_core::PropertyIntReadError),
    #[error("{0}")]
    PropertyString(#[from] cosmolkit_core::PropertyStringError),
    #[error("{0}")]
    Rings(#[from] cosmolkit_core::RingFindingError),
    #[error("{0}")]
    Stereo(#[from] cosmolkit_core::StereoError),
    #[error("{0}")]
    DoubleBondStereo(#[from] cosmolkit_core::DoubleBondStereoError),
    #[error("{0}")]
    BondDirectionStereo(#[from] cosmolkit_core::BondDirectionStereoError),
    #[error("{0}")]
    Atropisomer(#[from] cosmolkit_core::AtropisomerError),
    #[error("{0}")]
    Sanitize(#[from] cosmolkit_core::SanitizeError),
    #[error("{0}")]
    Hydrogen(#[from] cosmolkit_core::HydrogenError),
    #[error("{0}")]
    LegacyStereo(#[from] cosmolkit_core::LegacyStereoError),
    #[error("{0}")]
    AttachmentExpansion(#[from] cosmolkit_core::AttachmentExpansionError),
    #[error("{0}")]
    QueryAtom(#[from] cosmolkit_model::QueryAtomConversionError),
    #[error("{0}")]
    QueryGraph(#[from] cosmolkit_model::QueryGraphError),
    #[error("{0}")]
    QueryState(#[from] cosmolkit_model::QueryStateError),
    #[error("{0}")]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error("{0}")]
    TopologyEdit(#[from] cosmolkit_model::TopologyEditError),
    #[error("{0}")]
    PropertyValue(#[from] cosmolkit_model::PropertyValueError),
    #[error("{0}")]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error("{0}")]
    BondValue(#[from] cosmolkit_model::BondValueError),
    #[error("{0}")]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error("{0}")]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MolPostError {
    #[error(
        "atom {atom} property {property} unsigned value {value} causes positive_overflow converting UInt to signed int"
    )]
    UnsignedPropertyOverflow {
        atom: AtomId,
        property: &'static str,
        value: u32,
    },
    #[error("atom {atom} property {property} has invalid kind {kind:?}")]
    InvalidPropertyKind {
        atom: AtomId,
        property: &'static str,
        kind: cosmolkit_model::PropertyValueKind,
    },
    #[error("invalid attachment value on atom {atom}: {value:?}")]
    AttachmentValue { atom: AtomId, value: PropertyText },
    #[error(
        "bad_any_cast reading SGroup {group:?} property {property} kind {kind:?} as StringVector"
    )]
    SGroupPropertyKind {
        group: SubstanceGroupId,
        property: &'static str,
        kind: cosmolkit_model::PropertyValueKind,
    },
    #[error("attachment-point expansion failed: {0}")]
    AttachmentExpansion(#[source] MolProcessingError),
    #[error("Molfile postprocessing property is outside the detached model: {0}")]
    Representation(&'static str),
    #[error("invalid numeric data {value:?} in DAT SGroup {field}")]
    DataFieldNumber {
        field: &'static str,
        value: PropertyText,
    },
    #[error(transparent)]
    QueryAtomConversion(#[from] QueryAtomConversionError),
    #[error("Molfile postprocessing failed: {0}")]
    Processing(#[source] MolProcessingError),
}

fn parse_int_property(value: &PropertyValue) -> Result<i32, ()> {
    // The canonical CORE converter owns RDValue arithmetic reads, including
    // byte strings, C-locale right trimming, signed bounds and wrong kinds.
    cosmolkit_core::property_value_to_int(value).map_err(|_| ())
}

fn source_int_property_or_zero(
    atom: &cosmolkit_model::Atom,
    key: &'static str,
) -> Result<i32, MolPostError> {
    // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:441-450
    // RDKit❗✔️: template <>
    // RDKit❗✔️: inline int rdvalue_cast<int>(RDValue_cast_t v) {
    // RDKit❗✔️:   if (rdvalue_is<int>(v)) {
    // RDKit❗✔️:     return v.value.i;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (rdvalue_is<unsigned int>(v)) {
    // RDKit❗✔️:     return boost::numeric_cast<int>(v.value.u);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   throw std::bad_any_cast();
    // RDKit❗✔️: }
    // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:441-450

    // RDKit❗✔️:     if (atom->getPropIfPresent(common_properties::molSubstCount, ival) &&
    // Vector cast errors are lazy, source-key-specific and cannot become zero.
    // Other scalar conditions are retained; no allocation or full-map preflight.
    match atom.prop(key) {
        None => Ok(0),
        Some(value) => cosmolkit_core::property_value_to_int(value).map_err(|error| match error {
            cosmolkit_core::PropertyIntReadError::UnsignedOverflow { value } => {
                MolPostError::UnsignedPropertyOverflow {
                    atom: atom.id(),
                    property: key,
                    value,
                }
            }
            cosmolkit_core::PropertyIntReadError::InvalidKind { kind } => {
                MolPostError::InvalidPropertyKind {
                    atom: atom.id(),
                    property: key,
                    kind,
                }
            }
            error => MolPostError::Processing(MolProcessingError::PropertyInt(error)),
        }),
    }
}

fn property_diagnostic(value: &PropertyValue) -> PropertyText {
    match value {
        PropertyValue::String(value) => value.clone(),
        PropertyValue::UInt(value) => value.to_string().into(),
        _ => format!("<{:?}>", value.kind()).into(),
    }
}

fn data_values(group: &SubstanceGroup) -> Result<&[PropertyText], MolPostError> {
    // RDKit❗✔️:   template <typename T>
    // RDKit❗✔️:   bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit❗✔️:     return d_props.getValIfPresent(key, res);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   template <typename T>
    // RDKit❗✔️:   bool getValIfPresent(const std::string_view what, T &res) const {
    // RDKit❗✔️:     for (const auto &data : _data) {
    // RDKit❗✔️:       if (data.key == what) {
    // RDKit❗✔️:         res = from_rdvalue<T>(data.val);
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️: template <>
    // RDKit❗✔️: inline std::vector<std::string> rdvalue_cast<std::vector<std::string>>(
    // RDKit❗✔️:     RDValue_cast_t v) {
    // RDKit❗✔️:   if (rdvalue_is<std::vector<std::string>>(v)) {
    // RDKit❗✔️:     return *v.ptrCast<std::vector<std::string>>();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   throw std::bad_any_cast();
    // RDKit❗✔️: }
    // Actual present DATAFIELDS is authoritative: only tag StringVector has
    // the source vector cast. Absence selects the existing detached projection,
    // never a fallback for a reached wrong-kind value. Borrowing avoids the
    // source vector copy; lookup and payload traversal stay linear.
    match group.props().get(b"DATAFIELDS".as_slice()) {
        Some(PropertyValue::StringVector(values)) => Ok(values),
        Some(value) => Err(MolPostError::SGroupPropertyKind {
            group: group.id(),
            property: "DATAFIELDS",
            kind: value.kind(),
        }),
        None => Ok(group
            .data()
            .map_or(group.data_fields(), |data| data.values.as_slice())),
    }
}
fn source_sgroup_string(
    group: &SubstanceGroup,
    key: &str,
    projection: Option<&PropertyText>,
) -> Result<Option<PropertyText>, MolPostError> {
    // Scalar source reads use the single CORE getProp<string> conversion.
    // Only absent source records select existing explicit detached fields;
    // neither actual tags nor raw bytes are reconstructed from field names.
    group
        .props()
        .get(key.as_bytes())
        .map(cosmolkit_core::property_value_to_string)
        .transpose()
        .map(|value| value.or_else(|| projection.cloned()))
        .map_err(|error| MolPostError::Processing(MolProcessingError::PropertyString(error)))
}
fn source_sgroup_is_data(group: &SubstanceGroup) -> Result<bool, MolPostError> {
    let source = source_sgroup_string(group, "TYPE", None)?;
    Ok(source.as_ref().map_or(
        matches!(group.kind(), cosmolkit_model::SubstanceGroupKind::Data),
        |kind| kind.as_bytes() == b"DAT",
    ))
}

fn retain_substance_groups(
    groups: Vec<SubstanceGroup>,
    remove: &[bool],
) -> Result<Vec<SubstanceGroup>, MolPostError> {
    let map_len = groups
        .iter()
        .map(|group| group.id().index())
        .max()
        .map_or(0, |maximum| maximum + 1);
    let mut old_to_new = vec![None; map_len];
    let mut retained =
        Vec::with_capacity(groups.len() - remove.iter().filter(|flag| **flag).count());
    for (position, group) in groups.into_iter().enumerate() {
        if !remove[position] {
            old_to_new[group.id().index()] = Some(SubstanceGroupId::new(retained.len()));
            retained.push(group);
        }
    }
    for group in &mut retained {
        let new_id = old_to_new[group.id().index()].ok_or(MolPostError::Representation(
            "retained SGroup id was not mapped",
        ))?;
        group.set_id(new_id);
        if let Some(old_parent) = group.parent() {
            let new_parent = old_to_new
                .get(old_parent.index())
                .and_then(|mapped| *mapped)
                .ok_or(MolPostError::Representation(
                    "retained SGroup references a consumed parent",
                ))?;
            group.set_parent(new_parent);
        }
    }
    Ok(retained)
}

fn process_mrv_coordinate_bond(
    atoms: &mut [cosmolkit_model::Atom],
    bonds: &mut [cosmolkit_model::Bond],
    query_bonds: Option<&mut [QueryBond]>,
    group: &SubstanceGroup,
) -> Result<(), MolPostError> {
    // RDKit✔️❌: void processMrvCoordinateBond(RWMol &mol, const SubstanceGroup &sg) {
    // RDKit✔️❌:   std::vector<std::string> dataFields;
    // RDKit✔️❌:   if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    // RDKit✔️❌:     if (dataFields.empty()) {
    // RDKit✔️❌:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:           << "ignoring MRV_COORDINATE_BOND_TYPE SGroup without data fields."
    // RDKit✔️❌:           << std::endl;
    // RDKit✔️❌:       return;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     auto coordinate_bond_idx =
    // RDKit✔️❌:         FileParserUtils::toUnsigned(dataFields[0], true) - 1;
    // RDKit✔️❌:
    // RDKit✔️❌:     if (dataFields.size() > 1) {
    // RDKit✔️❌:       BOOST_LOG(rdWarningLog) << "ignoring extra data fields in "
    // RDKit✔️❌:                                  "MRV_COORDINATE_BOND_TYPE SGroup for bond "
    // RDKit✔️❌:                               << coordinate_bond_idx << '.' << std::endl;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     Bond *old_bond = nullptr;
    // RDKit✔️❌:     try {
    // RDKit✔️❌:       old_bond = mol.getBondWithIdx(coordinate_bond_idx);
    // RDKit✔️❌:     } catch (const Invar::Invariant &) {
    // RDKit✔️❌:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:           << "molecule does not contain a bond matching the "
    // RDKit✔️❌:              "MRV_COORDINATE_BOND_TYPE SGroup for bond "
    // RDKit✔️❌:           << coordinate_bond_idx << ", ignoring." << std::endl;
    // RDKit✔️❌:       return;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (!old_bond || old_bond->getBondType() != Bond::BondType::UNSPECIFIED) {
    // RDKit✔️❌:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:           << "MRV_COORDINATE_BOND_TYPE SGroup with value "
    // RDKit✔️❌:           << coordinate_bond_idx
    // RDKit✔️❌:           << " does not reference a query bond, ignoring." << std::endl;
    // RDKit✔️❌:       return;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     Bond new_bond(Bond::BondType::DATIVE);
    // RDKit✔️❌:     auto preserveProps = true;
    // RDKit✔️❌:     auto keepSGroups = true;
    // RDKit✔️❌:     mol.replaceBond(coordinate_bond_idx, &new_bond, preserveProps, keepSGroups);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // RDKit✔️❌: void RWMol::replaceBond(unsigned int idx, Bond *bond_pin, bool preserveProps,
    // RDKit✔️❌:                         bool keepSGroups) {
    // RDKit✔️❌:   PRECONDITION(bond_pin, "bad bond passed to replaceBond");
    // RDKit✔️❌:   URANGE_CHECK(idx, getNumBonds());
    // RDKit✔️❌:   auto bIter = getEdges();
    // RDKit✔️❌:   for (unsigned int i = 0; i < idx; i++) {
    // RDKit✔️❌:     ++bIter.first;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   const auto *obond = d_graph[*(bIter.first)];
    // RDKit✔️❌:   auto *bond_p = bond_pin->copy();
    // RDKit✔️❌:   bond_p->setOwningMol(this);
    // RDKit✔️❌:   bond_p->setIdx(idx);
    // RDKit✔️❌:   bond_p->setBeginAtomIdx(obond->getBeginAtomIdx());
    // RDKit✔️❌:   bond_p->setEndAtomIdx(obond->getEndAtomIdx());
    // RDKit✔️❌:
    // RDKit✔️❌:   // Update explicit Hs, if set, on both ends. This was github #7128
    // RDKit✔️❌:   auto orderDifference =
    // RDKit✔️❌:       bond_p->getBondTypeAsDouble() - obond->getBondTypeAsDouble();
    // RDKit✔️❌:   if (orderDifference > 0) {
    // RDKit✔️❌:     for (auto atom : {bond_p->getBeginAtom(), bond_p->getEndAtom()}) {
    // RDKit✔️❌:       if (auto explicit_hs = atom->getNumExplicitHs(); explicit_hs > 0) {
    // RDKit✔️❌:         auto new_hs = static_cast<int>(explicit_hs - orderDifference);
    // RDKit✔️❌:         atom->setNumExplicitHs(std::max(new_hs, 0));
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (preserveProps) {
    // RDKit✔️❌:     const bool replaceExistingData = false;
    // RDKit✔️❌:     bond_p->updateProps(*d_graph[*(bIter.first)], replaceExistingData);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   const auto orig_p = d_graph[*(bIter.first)];
    // RDKit✔️❌:   delete orig_p;
    // RDKit✔️❌:   d_graph[*(bIter.first)] = bond_p;
    // RDKit✔️❌:
    // RDKit✔️❌:   if (!keepSGroups) {
    // RDKit✔️❌:     removeSubstanceGroupsReferencingBond(*this, idx);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // handle bookmarks
    // RDKit✔️❌:   for (auto &ab : d_bondBookmarks) {
    // RDKit✔️❌:     for (auto &elem : ab.second) {
    // RDKit✔️❌:       if (elem == orig_p) {
    // RDKit✔️❌:         elem = bond_p;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: };
    // RDKit✔️❌: double Bond::getBondTypeAsDouble() const {
    // RDKit✔️❌:   double res;
    // RDKit✔️❌:   switch (getBondType()) {
    // RDKit✔️❌:     case UNSPECIFIED:
    // RDKit✔️❌:     case IONIC:
    // RDKit✔️❌:     case ZERO:
    // RDKit✔️❌:       res = 0;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case SINGLE:
    // RDKit✔️❌:       res = 1;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case DOUBLE:
    // RDKit✔️❌:       res = 2;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case TRIPLE:
    // RDKit✔️❌:       res = 3;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case QUADRUPLE:
    // RDKit✔️❌:       res = 4;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case QUINTUPLE:
    // RDKit✔️❌:       res = 5;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case HEXTUPLE:
    // RDKit✔️❌:       res = 6;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case ONEANDAHALF:
    // RDKit✔️❌:       res = 1.5;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case TWOANDAHALF:
    // RDKit✔️❌:       res = 2.5;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case THREEANDAHALF:
    // RDKit✔️❌:       res = 3.5;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case FOURANDAHALF:
    // RDKit✔️❌:       res = 4.5;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case FIVEANDAHALF:
    // RDKit✔️❌:       res = 5.5;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case AROMATIC:
    // RDKit✔️❌:       res = 1.5;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case DATIVEONE:
    // RDKit✔️❌:       res = 1.0;
    // RDKit✔️❌:       break;  // FIX: this should probably be different
    // RDKit✔️❌:     case DATIVE:
    // RDKit✔️❌:       res = 1.0;
    // RDKit✔️❌:       break;  // FIX: again probably wrong
    // RDKit✔️❌:     case HYDROGEN:
    // RDKit✔️❌:       res = 0.0;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     default:
    // RDKit✔️❌:       UNDER_CONSTRUCTION("Bad bond type");
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // Behavior review: a new ordinary DATIVE Bond replaces only an
    // UNSPECIFIED target; row ID/endpoints and ordered typed properties survive,
    // while direction, stereo, flags and query identity start at defaults.
    // Source order difference is 1 - 0, so positive explicit Hs at both ends
    // decrease by one, including the donor. Other orders/out-of-range rows are
    // source-defined no-ops; malformed numeric data propagates before lookup.
    // Complexity review: one numeric scan and indexed lookup; copying only the
    // replaced bond's properties follows preserveProps. Ordered tree insertion
    // and the extra owned query carrier copy exceed the single source Bond;
    // the performance axis records that cost. No entire group payload clone.
    let Some(value) = data_values(group)?.first() else {
        return Ok(());
    };
    let row = parse_rdkit_unsigned(value)
        .map_err(|()| MolPostError::DataFieldNumber {
            field: "MRV_COORDINATE_BOND_TYPE",
            value: value.clone(),
        })?
        .wrapping_sub(1) as usize;
    let Some(old) = bonds
        .get(row)
        .filter(|bond| bond.order() == BondOrder::Unspecified)
    else {
        return Ok(());
    };
    let mut replacement = cosmolkit_model::Bond::from_spec(
        old.id(),
        cosmolkit_model::BondSpec::new(old.begin(), old.end(), BondOrder::Dative),
    );
    replacement
        .replace_property_records(
            cosmolkit_model::ordered_bond_properties(old)
                .map(|(key, value)| (key.clone(), value.clone())),
        )
        .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    for id in [old.begin(), old.end()] {
        let atom = &mut atoms[id.index()];
        atom.set_explicit_hydrogens(atom.explicit_hydrogens().saturating_sub(1));
    }
    if let Some(query_bonds) = query_bonds {
        query_bonds[row] = QueryBond::from_carrier_parts(
            replacement.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Dative)),
        );
    }
    bonds[row] = replacement;
    Ok(())
}

fn process_mrv_implicit_h(
    atoms: &mut [cosmolkit_model::Atom],
    bonds: &[cosmolkit_model::Bond],
    adjacency: &AdjacencyList,
    group: &SubstanceGroup,
) -> Result<(), MolPostError> {
    // RDKit✔️✔️: void processMrvImplicitH(RWMol &mol, const SubstanceGroup &sg) {
    // RDKit✔️✔️:   std::vector<std::string> dataFields;
    // RDKit✔️✔️:   if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    // RDKit✔️✔️:     for (const auto &df : dataFields) {
    // RDKit✔️✔️:       if (df.substr(0, 6) == "IMPL_H") {
    // RDKit✔️✔️:         auto val = FileParserUtils::toInt(df.substr(6));
    // RDKit✔️✔️:         for (auto atIdx : sg.getAtoms()) {
    // RDKit✔️✔️:           if (atIdx < mol.getNumAtoms()) {
    // RDKit✔️✔️:             // if the atom has aromatic bonds to it, then set the explicit
    // RDKit✔️✔️:             // value, otherwise skip it.
    // RDKit✔️✔️:             auto atom = mol.getAtomWithIdx(atIdx);
    // RDKit✔️✔️:             bool hasAromaticBonds = false;
    // RDKit✔️✔️:             for (auto bndI :
    // RDKit✔️✔️:                  boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit✔️✔️:               auto bnd = (mol)[bndI];
    // RDKit✔️✔️:               if (bnd->getIsAromatic() ||
    // RDKit✔️✔️:                   bnd->getBondType() == Bond::AROMATIC) {
    // RDKit✔️✔️:                 hasAromaticBonds = true;
    // RDKit✔️✔️:                 break;
    // RDKit✔️✔️:               }
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             if (hasAromaticBonds) {
    // RDKit✔️✔️:               atom->setNumExplicitHs(val);
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:                   << "MRV_IMPLICIT_H SGroup on atom without aromatic "
    // RDKit✔️✔️:                      "bonds, "
    // RDKit✔️✔️:                   << atIdx << ", ignored." << std::endl;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:                 << "bad atom index, " << atIdx
    // RDKit✔️✔️:                 << ", found in MRV_IMPLICIT_H SGroup. Ignoring it."
    // RDKit✔️✔️:                 << std::endl;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior review: convert each IMPL_H suffix before visiting atoms;
    // lexical errors propagate, while only aromatic neighbors trigger a write.
    // Checked model width errors occur only at the actual source setter.
    // Complexity review: one data scan plus direct per-atom adjacency traversal;
    // groups/data and bonds are borrowed, with no repeated graph scan.
    for value in data_values(group)? {
        let Some(raw) = value.as_bytes().strip_prefix(b"IMPL_H") else {
            continue;
        };
        let count = parse_rdkit_int(raw).map_err(|()| MolPostError::DataFieldNumber {
            field: "MRV_IMPLICIT_H",
            value: raw.to_vec().into(),
        })?;
        for id in group.atoms() {
            let Some(atom) = atoms.get_mut(id.index()) else {
                continue;
            };
            if adjacency.neighbors_of(id.index()).iter().any(|neighbor| {
                let bond = &bonds[neighbor.bond.index()];
                bond.is_aromatic() || bond.order() == BondOrder::Aromatic
            }) {
                atom.set_explicit_hydrogens(u8::try_from(count).map_err(|_| {
                    MolPostError::Representation("MRV_IMPLICIT_H count outside u8")
                })?);
            }
        }
    }
    Ok(())
}

fn process_zbo(bonds: &mut [cosmolkit_model::Bond], group: &SubstanceGroup) {
    // RDKit✔️✔️: void processZBO(RWMol &mol, const SubstanceGroup &sg) {
    // RDKit✔️✔️:   for (auto bidx : sg.getBonds()) {
    // RDKit✔️✔️:     auto bond = mol.getBondWithIdx(bidx);
    // RDKit✔️✔️:     bond->setBondType(Bond::BondType::ZERO);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior/complexity review: validated SGroup membership indexes the
    // existing rows once; only the bond order changes, without replacement.
    for id in group.bonds() {
        bonds[id.index()].set_order(BondOrder::Zero);
    }
}

fn process_zch(
    atoms: &mut [cosmolkit_model::Atom],
    group: &SubstanceGroup,
) -> Result<(), MolPostError> {
    // RDKit✔️✔️: void processZCH(RWMol &mol, const SubstanceGroup &sg) {
    // RDKit✔️✔️:   RDUNUSED_PARAM(mol);
    // RDKit✔️✔️:   std::vector<std::string> dataFields;
    // RDKit✔️✔️:   if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    // RDKit✔️✔️:     if (dataFields.empty()) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "ignoring ZCHG SGroup without data fields." << std::endl;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto &df : dataFields) {
    // RDKit✔️✔️:       std::string trimmed = boost::trim_copy(df);
    // RDKit✔️✔️:       std::vector<std::string> splitLine;
    // RDKit✔️✔️:       boost::split(splitLine, trimmed, boost::is_any_of(";"),
    // RDKit✔️✔️:                    boost::token_compress_off);
    // RDKit✔️✔️:       const auto &aids = sg.getAtoms();
    // RDKit✔️✔️:       if (splitLine.size() < aids.size()) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "DATAFIELDS in ZCH SGroup is shorter than the number of atoms in the SGroup. Ignoring it."
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (auto i = 0u; i < aids.size(); ++i) {
    // RDKit✔️✔️:         auto aid = aids[i];
    // RDKit✔️✔️:         auto atom = mol.getAtomWithIdx(aid);
    // RDKit✔️✔️:         auto val = 0;
    // RDKit✔️✔️:         if (!splitLine[i].empty()) {
    // RDKit✔️✔️:           val = FileParserUtils::toInt(splitLine[i]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atom->setFormalCharge(val);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior review: short rows are ignored; empty tokens mean zero.
    // Invalid nonempty tokens propagate the source numeric conversion error.
    // Complexity review: splitting is linear and no group/data is cloned.
    for value in data_values(group)? {
        let values = value
            .as_bytes()
            .trim_ascii()
            .split(|byte| *byte == b';')
            .collect::<Vec<_>>();
        if values.len() < group.atoms().len() {
            continue;
        }
        for (id, text) in group.atoms().iter().zip(values) {
            let parsed = if text.is_empty() {
                0
            } else {
                parse_rdkit_int(text).map_err(|()| MolPostError::DataFieldNumber {
                    field: "ZCH",
                    value: text.to_vec().into(),
                })?
            };
            atoms[id.index()].set_formal_charge(
                i8::try_from(parsed)
                    .map_err(|_| MolPostError::Representation("ZCH charge outside i8"))?,
            );
        }
    }
    Ok(())
}

fn process_hyd(
    atoms: &mut [cosmolkit_model::Atom],
    group: &SubstanceGroup,
) -> Result<(), MolPostError> {
    // RDKit✔️✔️: void processHYD(RWMol &mol, const SubstanceGroup &sg) {
    // RDKit✔️✔️:   std::vector<std::string> dataFields;
    // RDKit✔️✔️:   if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    // RDKit✔️✔️:     if (dataFields.empty()) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "ignoring HYD SGroup without data fields." << std::endl;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto &df : dataFields) {
    // RDKit✔️✔️:       std::string trimmed = boost::trim_copy(df);
    // RDKit✔️✔️:       std::vector<std::string> splitLine;
    // RDKit✔️✔️:       boost::split(splitLine, trimmed, boost::is_any_of(";"),
    // RDKit✔️✔️:                    boost::token_compress_off);
    // RDKit✔️✔️:       const auto &aids = sg.getAtoms();
    // RDKit✔️✔️:       if (splitLine.size() < aids.size()) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "DATAFIELDS in HYD SGroup is shorter than the number of atoms in the SGroup. Ignoring it."
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (auto i = 0u; i < aids.size(); ++i) {
    // RDKit✔️✔️:         auto aid = aids[i];
    // RDKit✔️✔️:         auto atom = mol.getAtomWithIdx(aid);
    // RDKit✔️✔️:         auto val = 0;
    // RDKit✔️✔️:         if (!splitLine[i].empty()) {
    // RDKit✔️✔️:           val = FileParserUtils::toInt(splitLine[i]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atom->setProp("_ZBO_H", true);
    // RDKit✔️✔️:         atom->setNumExplicitHs(val);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior review: preserves the source short-row/empty-token rules and
    // typed numeric error, then sets the established _ZBO_H marker and H count.
    // Complexity review: one split/application pass; all group payload borrowed.
    for value in data_values(group)? {
        let values = value
            .as_bytes()
            .trim_ascii()
            .split(|byte| *byte == b';')
            .collect::<Vec<_>>();
        if values.len() < group.atoms().len() {
            continue;
        }
        for (id, text) in group.atoms().iter().zip(values) {
            let parsed = if text.is_empty() {
                0
            } else {
                parse_rdkit_int(text).map_err(|()| MolPostError::DataFieldNumber {
                    field: "HYD",
                    value: text.to_vec().into(),
                })?
            };
            let atom = &mut atoms[id.index()];
            atom.set_prop("_ZBO_H", PropertyValue::Bool(true))
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            atom.set_explicit_hydrogens(
                u8::try_from(parsed)
                    .map_err(|_| MolPostError::Representation("HYD count outside u8"))?,
            );
        }
    }
    Ok(())
}

fn process_groups_on_topology(
    topology: &mut TopologyBlock,
    mut query_bonds: Option<&mut [QueryBond]>,
) -> Result<(), MolPostError> {
    // RDKit✔️✔️: void processSGroups(RWMol *mol) {
    // RDKit✔️✔️:   std::vector<unsigned int> sgsToRemove;
    // RDKit✔️✔️:   unsigned int sgIdx = 0;
    // RDKit✔️✔️:   for (auto &sg : getSubstanceGroups(*mol)) {
    // RDKit✔️✔️:     if (sg.getProp<std::string>("TYPE") == "DAT") {
    // RDKit✔️✔️:       std::string field;
    // RDKit✔️✔️:       if (sg.getPropIfPresent("FIELDNAME", field)) {
    // RDKit✔️✔️:         if (field == "MRV_COORDINATE_BOND_TYPE") {
    // RDKit✔️✔️:           // V2000 support for coordinate bonds
    // RDKit✔️✔️:           processMrvCoordinateBond(*mol, sg);
    // RDKit✔️✔️:           sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         } else if (field == "MRV_IMPLICIT_H") {
    // RDKit✔️✔️:           // CXN extension to specify implicit Hs, used for aromatic rings
    // RDKit✔️✔️:           processMrvImplicitH(*mol, sg);
    // RDKit✔️✔️:           sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         } else if (field == "ZBO") {
    // RDKit✔️✔️:           // RDKit extension for zero-order bonds
    // RDKit✔️✔️:           processZBO(*mol, sg);
    // RDKit✔️✔️:           sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         } else if (field == "ZCH") {
    // RDKit✔️✔️:           // RDKit extension for charge on atoms involved in zero-order bonds
    // RDKit✔️✔️:           processZCH(*mol, sg);
    // RDKit✔️✔️:           sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         } else if (field == "HYD") {
    // RDKit✔️✔️:           // RDKit extension for hydrogen-count on atoms involved in
    // RDKit✔️✔️:           // zero-order bonds
    // RDKit✔️✔️:           processHYD(*mol, sg);
    // RDKit✔️✔️:           sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (sg.getPropIfPresent("QUERYTYPE", field) &&
    // RDKit✔️✔️:           (field == "SMARTSQ" || field == "SQ")) {
    // RDKit✔️✔️:         processSMARTSQ(*mol, sg);
    // RDKit✔️✔️:         sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++sgIdx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // now remove the S groups we processed, we saved indices so do this in
    // RDKit✔️✔️:   // backwards
    // RDKit✔️✔️:   auto &sgs = getSubstanceGroups(*mol);
    // RDKit✔️✔️:   for (auto it = sgsToRemove.rbegin(); it != sgsToRemove.rend(); ++it) {
    // RDKit✔️✔️:     sgs.erase(sgs.begin() + *it);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior review: process recognized DAT groups in input order; SMARTSQ
    // remains on the query carrier until its query action. Consume recognized
    // groups only after successful processing, preserving other groups/parents.
    // Complexity review: split field borrows avoid cloning any group/data.
    // Membership-local work and one ordered retain pass replace source erases.
    let TopologyBlock {
        atoms,
        bonds,
        adjacency,
        substance_groups,
        ..
    } = topology;
    let mut remove = vec![false; substance_groups.len()];
    for (index, group) in substance_groups.iter().enumerate() {
        if !source_sgroup_is_data(group)? {
            continue;
        }
        let Some(field) = source_sgroup_string(
            group,
            "FIELDNAME",
            group.data().and_then(|data| data.field_name.as_ref()),
        )?
        else {
            continue;
        };
        let field = field.as_bytes();
        match field {
            b"MRV_COORDINATE_BOND_TYPE" => {
                process_mrv_coordinate_bond(atoms, bonds, query_bonds.as_deref_mut(), group)?
            }
            b"MRV_IMPLICIT_H" => process_mrv_implicit_h(atoms, bonds, adjacency, group)?,
            b"ZBO" => process_zbo(bonds, group),
            b"ZCH" => process_zch(atoms, group)?,
            b"HYD" => process_hyd(atoms, group)?,
            _ => continue,
        }
        remove[index] = true;
    }
    *substance_groups = retain_substance_groups(std::mem::take(substance_groups), &remove)?;
    Ok(())
}

fn process_atom_properties(
    topology: &mut TopologyBlock,
    mut query_atoms: Option<&mut [QueryAtom]>,
) -> Result<(), MolPostError> {
    // RDKit✔️✔️: void ProcessMolProps(RWMol *mol) {
    // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
    // RDKit✔️✔️:   // we have to loop the ugly way because we may need to actually replace an
    // RDKit✔️✔️:   // atom
    // RDKit✔️✔️:   for (unsigned int aidx = 0; aidx < mol->getNumAtoms(); ++aidx) {
    // RDKit✔️✔️:     auto atom = mol->getAtomWithIdx(aidx);
    // RDKit✔️✔️:     int ival = 0;
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molSubstCount, ival) &&
    // RDKit✔️✔️:         ival != 0) {
    // RDKit✔️✔️:       if (!atom->hasQuery()) {
    // RDKit✔️✔️:         atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       bool gtQuery = false;
    // RDKit✔️✔️:       if (ival == -1) {
    // RDKit✔️✔️:         ival = 0;
    // RDKit✔️✔️:       } else if (ival == -2) {
    // RDKit✔️✔️:         // as drawn
    // RDKit✔️✔️:         ival = atom->getDegree();
    // RDKit✔️✔️:       } else if (ival >= 6) {
    // RDKit✔️✔️:         // 6 or more
    // RDKit✔️✔️:         gtQuery = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!gtQuery) {
    // RDKit✔️✔️:         atom->expandQuery(makeAtomExplicitDegreeQuery(ival));
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         // create a temp query the normal way so that we can be sure to get
    // RDKit✔️✔️:         // the description right
    // RDKit✔️✔️:         std::unique_ptr<ATOM_EQUALS_QUERY> tmp{
    // RDKit✔️✔️:             makeAtomExplicitDegreeQuery(ival)};
    // RDKit✔️✔️:         atom->expandQuery(makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(
    // RDKit✔️✔️:             ival, tmp->getDataFunc(),
    // RDKit✔️✔️:             std::string("less_") + tmp->getDescription()));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molTotValence, ival) &&
    // RDKit✔️✔️:         ival != 0 && !atom->hasProp("_ZBO_H")) {
    // RDKit✔️✔️:       atom->setNoImplicit(true);
    // RDKit✔️✔️:       if (ival == 15     // V2000
    // RDKit✔️✔️:           || ival == -1  // v3000
    // RDKit✔️✔️:       ) {
    // RDKit✔️✔️:         atom->setNumExplicitHs(0);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (static_cast<int>(atom->getValence(Atom::ValenceType::EXPLICIT)) >
    // RDKit✔️✔️:             ival) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "atom " << atom->getIdx() << " has specified valence (" << ival
    // RDKit✔️✔️:               << ") smaller than the drawn valence "
    // RDKit✔️✔️:               << atom->getValence(Atom::ValenceType::EXPLICIT) << "."
    // RDKit✔️✔️:               << std::endl;
    // RDKit✔️✔️:           atom->setNumExplicitHs(0);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           atom->setNumExplicitHs(ival -
    // RDKit✔️✔️:                                  atom->getValence(Atom::ValenceType::EXPLICIT));
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     atom->clearProp(common_properties::molTotValence);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   processSGroups(mol);
    // RDKit✔️✔️: }
    // Behavior review: the optional query slice is present exactly after the
    // Molfile owner has promoted a record containing SUBST/SMARTSQ/query
    // syntax. SUBST expands the existing predicate in source order and marks
    // the carrier as a real source QueryAtom for downstream RemoveHs. The
    // total-valence branch observes the same current carrier topology.
    // Complexity review: move the existing tree through the sole source AND
    // helper; no deep predicate clone. The temporary empty vector allocates
    // nothing and cannot escape. Per-atom valence traversal remains source-shaped.
    for index in 0..topology.atoms.len() {
        let substitution = source_int_property_or_zero(&topology.atoms[index], "molSubstCount")?;
        if substitution != 0 {
            let atoms = query_atoms
                .as_deref_mut()
                .ok_or(MolPostError::Representation(
                    "molSubstCount requires query record promotion",
                ))?;
            let degree = if substitution == -1 {
                0
            } else if substitution == -2 {
                i64::try_from(topology.adjacency.neighbors_of(index).len())
                    .map_err(|_| MolPostError::Representation("molSubstCount degree outside i64"))?
            } else {
                i64::from(substitution)
            };
            let predicate = if substitution >= 6 {
                let degree = u8::try_from(degree).map_err(|_| {
                    MolPostError::Representation("molSubstCount range target outside u8")
                })?;
                AtomQueryPredicate::ExplicitDegreeLessEqual(degree)
            } else {
                let degree = i32::try_from(degree).map_err(|_| {
                    MolPostError::Representation("molSubstCount query target outside i32")
                })?;
                AtomQueryPredicate::ExplicitDegree(degree)
            };
            let current =
                std::mem::replace(atoms[index].predicate_mut(), QueryNode::and(Vec::new()));
            atoms[index].set_predicate(crate::sdf::expand_molfile_atom_query(
                Some(current),
                QueryNode::predicate(predicate),
            ));
            topology.atoms[index]
                .set_prop("_MolFileAtomQuery", "1")
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            atoms[index]
                .set_prop("_MolFileAtomQuery", "1")
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        }
        let value = source_int_property_or_zero(&topology.atoms[index], "molTotValence")?;
        if value != 0 && topology.atoms[index].prop("_ZBO_H").is_none() {
            let explicit =
                calculate_explicit_valence_for_topology(topology, AtomId::new(index), false, false)
                    .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            let hydrogens = if value == 15 || value == -1 || explicit > value {
                0
            } else {
                u8::try_from(value - explicit)
                    .map_err(|_| MolPostError::Representation("molTotValence H count outside u8"))?
            };
            topology.atoms[index].set_no_implicit(true);
            topology.atoms[index].set_explicit_hydrogens(hydrogens);
        }
        topology.atoms[index]
            .clear_prop("molTotValence")
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    }
    Ok(())
}

fn apply_stereo_and_sanitize(
    mut topology: TopologyBlock,
    mut coordinates: CoordinateBlock,
    mut properties: cosmolkit_model::MoleculeProperties,
    chirality_possible: bool,
    params: MolPostParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<
    (
        TopologyBlock,
        CoordinateBlock,
        cosmolkit_model::MoleculeProperties,
        TopologyMapping,
        Option<(Vec<QueryAtom>, Vec<QueryBond>)>,
        MolPostDerivedState,
    ),
    MolPostError,
> {
    // BEGIN RDKIT CPP FUNCTION finishMolProcessing applicable stereo/sanitize closure
    // RDKit✔️✔️:   // update the chirality and stereo-chemistry
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // NOTE: we detect the stereochemistry before sanitizing/removing
    // RDKit✔️✔️:   // hydrogens because the removal of H atoms may actually remove
    // RDKit✔️✔️:   // the wedged bond from the molecule.  This wipes out the only
    // RDKit✔️✔️:   // sign that chirality ever existed and makes us sad... so first
    // RDKit✔️✔️:   // perceive chirality, then remove the Hs and sanitize.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   const Conformer &conf = res->getConformer();
    // RDKit✔️✔️:   if (chiralityPossible || conf.is3D()) {
    // RDKit✔️✔️:     if (!conf.is3D()) {
    // RDKit✔️✔️:       bool replaceExistingTags = true;
    // RDKit✔️✔️:       MolOps::assignChiralTypesFromBondDirs(*res, conf.getId(),
    // RDKit✔️✔️:                                             replaceExistingTags);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       res->updatePropertyCache(false);
    // RDKit✔️✔️:       MolOps::assignChiralTypesFrom3D(*res, conf.getId(), true);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Atropisomers::detectAtropisomerChirality(*res, &conf);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now that atom stereochem has been perceived, the wedging
    // RDKit✔️✔️:   // information is no longer needed, so we clear
    // RDKit✔️✔️:   // single bond dir flags:
    // RDKit✔️✔️:   MolOps::clearSingleBondDirFlags(*res);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (params.sanitize) {
    // RDKit✔️✔️:     if (params.removeHs) {
    // RDKit✔️✔️:       // Bond stereo detection must happen before H removal, or
    // RDKit✔️✔️:       // else we might be removing stereogenic H atoms in double
    // RDKit✔️✔️:       // bonds (e.g. imines). But before we run stereo detection,
    // RDKit✔️✔️:       // we need to run mol cleanup so don't have trouble with
    // RDKit✔️✔️:       // e.g. nitro groups. Sadly, this a;; means we will find
    // RDKit✔️✔️:       // run both cleanup and ring finding twice (a fast find
    // RDKit✔️✔️:       // rings in bond stereo detection, and another in
    // RDKit✔️✔️:       // sanitization's SSSR symmetrization).
    // RDKit✔️✔️:       unsigned int failedOp = 0;
    // RDKit✔️✔️:       MolOps::sanitizeMol(*res, failedOp, MolOps::SANITIZE_CLEANUP);
    // RDKit✔️✔️:       MolOps::detectBondStereochemistry(*res);
    // RDKit✔️✔️:       MolOps::removeHs(*res);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       MolOps::sanitizeMol(*res);
    // RDKit✔️✔️:       MolOps::detectBondStereochemistry(*res);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     MolOps::assignStereochemistry(*res, true, true, true);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     MolOps::detectBondStereochemistry(*res);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION finishMolProcessing applicable stereo/sanitize closure
    // Behavior review: source RWMol dynamic QueryAtom/QueryBond identity remains
    // available to every called owner. The detached adaptation carries the
    // validated typed query rows through sanitize, hydrogen removal and legacy
    // stereo, and remaps them with the same authoritative topology mapping.
    // Coordinate-driven double-bond directions retain source ordering. The
    // preceding bookmark, attachment, explicit-valence and ProcessMolProps
    // statements have distinct owners or the intentionally open attachment
    // gate and are not claimed by this helper.
    // Complexity review: each owner call is linear or owner-defined; this
    // orchestration introduces no repeated whole-graph clone beyond the owned
    // source-equivalent transform results.
    let mut final_state = MolPostDerivedState::default();
    let original_atom_count = topology.atoms.len();
    let original_bond_count = topology.bonds.len();
    let mut mapping = TopologyMapping::identity(original_atom_count, original_bond_count);
    let first_3d = coordinates.conformers_3d.first();
    if let Some(conformer) = first_3d {
        if conformer.is_3d() {
            let valence = assign_valence_for_topology(&topology, ValenceModel::RdkitLike)
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            topology = assign_chiral_tags_from_structure(
                &topology,
                &coordinates,
                &valence,
                &StructureTagParams::default(),
            )
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?
            .topology;
        } else if chirality_possible {
            assign_chiral_types_from_bond_dirs(&mut topology, conformer, true)
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        }
        let assignment =
            detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::ThreeD(conformer)))
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        apply_atropisomer_assignment(&mut topology, assignment)?;
    } else if let Some(conformer) = coordinates.conformers_2d.first() {
        if chirality_possible {
            let pseudo = Conformer3D::new(
                conformer.id(),
                conformer
                    .coordinates()
                    .iter()
                    .map(|xy| [xy[0], xy[1], 0.0])
                    .collect(),
                false,
            );
            assign_chiral_types_from_bond_dirs(&mut topology, &pseudo, true)
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        }
        let assignment =
            detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::TwoD(conformer)))
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        apply_atropisomer_assignment(&mut topology, assignment)?;
    }
    topology = clear_single_bond_directions(topology, false)
        .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    if params.sanitize {
        if params.remove_hs {
            topology = sanitize_topology_with_query_state(
                &topology,
                &SanitizeParams {
                    operations: SanitizeOperations::CLEANUP,
                },
                query_state,
            )
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?
            .topology;
            topology = detect_double_bond_stereochemistry(topology, &coordinates, &mut properties)?;
            let removed = remove_hydrogens_with_query_state(
                topology,
                coordinates,
                properties,
                &RemoveHsParams::default(),
                query_state,
            )
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            topology = removed.topology;
            coordinates = removed.coordinates;
            properties = removed.properties;
            mapping = removed.mapping;
        } else {
            topology = sanitize_topology_with_query_state(
                &topology,
                &SanitizeParams::default(),
                query_state,
            )
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?
            .topology;
            topology = detect_double_bond_stereochemistry(topology, &coordinates, &mut properties)?;
        }
        let remapped_query_rows = query_state
            .map(|state| remap_query_rows(state, &topology, &mapping))
            .transpose()
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        let final_query_state = remapped_query_rows
            .as_ref()
            .map(|(atoms, bonds)| QueryStateRef::try_for_topology(atoms, bonds, &topology))
            .transpose()
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        let valence = assign_valence_for_topology(&topology, ValenceModel::RdkitLike)
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        let rings = symmetrized_sssr(&topology, &RingSearchParams::default())
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        topology = assign_legacy_stereochemistry_with_query_state(
            topology,
            &valence,
            &rings,
            final_query_state,
        )
        .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        // RDKit❗✔️:     MolOps::assignStereochemistry(*res, true, true, true);
        // ROOT io44-move-sanitized-ring-carrier: retain the exact final owner
        // values already used above. Assignment moves are O(1); no finder,
        // intermediate-topology row, property inference or live installation.
        final_state = MolPostDerivedState {
            valence: Some(valence),
            rings: Some(rings),
        };
        // RDKit✔️✔️:   mol.setProp(common_properties::_StereochemDone, 1, true);
        // Behavior review: molecule-level computed properties are carried by
        // the detached `MoleculeProperties` block, so this is the implementing
        // location for the source wrapper's final property write.
        // Complexity review: one ordered-map property update matches the
        // source computed-property write and adds no graph traversal.
        properties
            .set_computed_prop("_StereochemDone", 1_i32)
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    } else {
        topology = detect_double_bond_stereochemistry(topology, &coordinates, &mut properties)?;
    }
    let query_rows = query_state
        .map(|state| remap_query_rows(state, &topology, &mapping))
        .transpose()
        .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    Ok((
        topology,
        coordinates,
        properties,
        mapping,
        query_rows,
        final_state,
    ))
}

fn apply_atropisomer_assignment(
    topology: &mut TopologyBlock,
    assignment: cosmolkit_core::AtropisomerAssignment,
) -> Result<(), MolPostError> {
    // RDKit✔️✔️:       bondToTry->getBeginAtom()->updatePropertyCache(false);
    // RDKit✔️✔️:       bondToTry->getEndAtom()->updatePropertyCache(false);
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // RDKit✔️✔️:     MolOps::setConjugation(mol);
    // RDKit✔️✔️:     MolOps::setHybridization(mol);
    // The core owner returns the source's ordered writes instead of mutating
    // the borrowed topology. Apply every effect, not only the final stereo.
    // Linear assignment over returned rows; no recomputation or graph copy.
    for (id, facts) in assignment.atom_valence_updates {
        topology.atoms[id.index()].set_source_valence_facts(facts);
    }
    if let Some(conjugated) = assignment.conjugated_bonds {
        for (bond, value) in topology.bonds.iter_mut().zip(conjugated) {
            bond.set_conjugated(value);
        }
    }
    if let Some(hybridization) = assignment.hybridization {
        for (atom, value) in topology.atoms.iter_mut().zip(hybridization.values) {
            atom.set_hybridization(value);
        }
    }
    for update in assignment.bond_updates {
        topology.bonds[update.bond.index()]
            .set_stereo(update.stereo)
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    }
    Ok(())
}

pub(super) fn detect_double_bond_stereochemistry(
    topology: TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &mut cosmolkit_model::MoleculeProperties,
) -> Result<TopologyBlock, MolPostError> {
    // BEGIN RDKIT CPP FUNCTION detectBondStereochemistry
    // RDKit✔️❌: void detectBondStereochemistry(ROMol &mol, int confId) {
    // RDKit✔️❌:   if (!mol.getNumConformers()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️❌:   setDoubleBondNeighborDirections(mol, &conf);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION detectBondStereochemistry
    // Behavior review: the first stored conformer is used regardless of its
    // independent is_3d flag; XY is lifted with exact positive-zero Z only for
    // the detached geometry kernel. This prepares directions without assigning
    // final bond stereo, matching the source phase boundary.
    // Complexity review: empty coordinates return before ring perception.
    // With 2D coordinates, lifting all n points allocates O(n) extra geometry
    // unlike source conformer access; retain a negative cost marker for this
    // detached adaptation. No topology clone is introduced here.
    if coordinates.conformers_3d.is_empty() && coordinates.conformers_2d.is_empty() {
        return Ok(topology);
    }
    let rings = symmetrized_sssr(&topology, &RingSearchParams::default())
        .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    if let Some(conformer) = coordinates.conformers_3d.first() {
        let update = set_double_bond_neighbor_directions(topology, &rings, Some(conformer))
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        if update.needs_detect_bond_stereo {
            // RDKit❗✔️:     mol.setProp("_needsDetectBondStereo", 1);
            properties
                .set_prop("_needsDetectBondStereo", 1_i32)
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        }
        return Ok(update.topology);
    }
    let Some(conformer) = coordinates.conformers_2d.first() else {
        return Ok(topology);
    };
    let lifted = Conformer3D::new(
        conformer.id(),
        conformer
            .coordinates()
            .iter()
            .map(|xy| [xy[0], xy[1], 0.0])
            .collect(),
        false,
    );
    let update = set_double_bond_neighbor_directions(topology, &rings, Some(&lifted))
        .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    if update.needs_detect_bond_stereo {
        // RDKit❗✔️:     mol.setProp("_needsDetectBondStereo", 1);
        properties
            .set_prop("_needsDetectBondStereo", 1_i32)
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    }
    Ok(update.topology)
}

fn concrete_to_query(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: cosmolkit_model::MoleculeProperties,
) -> Result<QueryMolBlockRecord, MolPostError> {
    let substance_groups = topology.substance_groups;
    let atoms = topology
        .atoms
        .into_iter()
        .map(|atom| {
            let predicate = query_from_concrete_atom_value(&atom);
            QueryAtom::from_carrier_parts(atom, predicate)
        })
        .collect();
    let bonds = topology
        .bonds
        .into_iter()
        .map(|bond| {
            let predicate = if bond.order() == BondOrder::Unspecified {
                QueryNode::predicate(BondQueryPredicate::Any)
            } else {
                QueryNode::predicate(BondQueryPredicate::Order(bond.order()))
            };
            QueryBond::from_carrier_parts(bond, predicate)
        })
        .collect();
    let mut props = properties.props().clone();
    if let Some(name) = properties.name() {
        props.insert("_Name".into(), PropertyValue::String(name.clone()));
    }
    let mut query = QueryGraph::from_parts(
        atoms,
        bonds,
        props,
        coordinates.conformers_2d,
        coordinates.conformers_3d,
        topology.stereo_groups,
    )
    .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    query
        .set_source_conformer_order(coordinates.source_conformer_order)
        .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    replace_query_substance_groups(&mut query, substance_groups)
        .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    Ok(QueryMolBlockRecord {
        query,
        properties,
        source_coordinate_dim: coordinates.source_coordinate_dim,
    })
}

fn is_source_query_atom(atom: &QueryAtom) -> bool {
    !atom.predicate_is_carrier_derived()
}

fn is_source_query_bond(bond: &QueryBond) -> bool {
    !bond.predicate_is_carrier_derived()
}

fn synchronize_query_atom(source: &mut QueryAtom, carrier: cosmolkit_model::Atom) {
    // RDKit❗✔️: Atom *replaceAtomWithQueryAtom(RWMol *mol, Atom *atom) {
    // RDKit❗✔️:   PRECONDITION(mol, "bad molecule");
    // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
    // RDKit❗✔️:   if (atom->hasQuery()) {
    // RDKit❗✔️:     return atom;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   QueryAtom qa(*atom);
    // RDKit❗✔️:   unsigned int idx = atom->getIdx();
    // RDKit❗✔️:
    // RDKit❗✔️:   if (atom->hasProp(common_properties::_hasMassQuery)) {
    // RDKit❗✔️:     qa.expandQuery(makeAtomMassQuery(static_cast<int>(atom->getMass())));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   mol->replaceAtom(idx, &qa);
    // RDKit❗✔️:   return mol->getAtomWithIdx(idx);
    // RDKit❗✔️: }
    // Behavior review: typed construction provenance distinguishes an
    // explicit detached query from Molfile IO's uniform wrapper around an
    // ordinary atom. Only the latter receives the constructor snapshot of the
    // final carrier; optional Molfile properties are not provenance.
    // Complexity review: one enum test; explicit carriers move the existing
    // predicate and the processed Atom into the new QueryAtom without cloning
    // the query tree. Synthesized carriers build one bounded predicate leaf.
    if is_source_query_atom(source) {
        let predicate = std::mem::replace(source.predicate_mut(), QueryNode::and(Vec::new()));
        *source = QueryAtom::from_parts(carrier, predicate);
    } else {
        let predicate = query_from_concrete_atom_value(&carrier);
        *source = QueryAtom::from_carrier_parts(carrier, predicate);
    }
}

fn synchronize_query_bond(mut source: QueryBond, carrier: cosmolkit_model::Bond) -> QueryBond {
    // BEGIN RDKIT CPP FUNCTION QueryBond::QueryBond(const Bond &)
    // RDKit✔️✔️:   explicit QueryBond(const Bond &other)
    // RDKit✔️✔️:       : Bond(other), dp_query(makeBondOrderEqualsQuery(other.getBondType())) {}
    // END RDKIT CPP FUNCTION QueryBond::QueryBond(const Bond &)
    // Behavior review: explicit detached query types retain their predicate;
    // a Molfile-owned uniform wrapper for an ordinary source Bond is
    // reconstructed from the final processed carrier, exactly where RDKit
    // still owns an ordinary Bond. Optional properties are not provenance.
    // Complexity review: classification is one enum test and reconstruction
    // allocates one leaf predicate.
    if is_source_query_bond(&source) {
        *source.bond_mut() = carrier;
        source
    } else {
        let predicate = if carrier.order() == BondOrder::Unspecified {
            QueryNode::predicate(BondQueryPredicate::Any)
        } else {
            QueryNode::predicate(BondQueryPredicate::Order(carrier.order()))
        };
        QueryBond::from_carrier_parts(carrier, predicate)
    }
}

fn record_requires_query(record: &MolBlockRecord) -> Result<bool, MolPostError> {
    match record {
        MolBlockRecord::Concrete { topology, .. } => {
            for atom in &topology.atoms {
                if source_int_property_or_zero(atom, "molSubstCount")? != 0 {
                    return Ok(true);
                }
            }
            for group in &topology.substance_groups {
                if !source_sgroup_is_data(group)? {
                    continue;
                }
                let query_type = source_sgroup_string(
                    group,
                    "QUERYTYPE",
                    group.data().and_then(|data| data.query_type.as_ref()),
                )?;
                if matches!(
                    query_type.as_ref().map(PropertyText::as_bytes),
                    Some(b"SMARTSQ" | b"SQ")
                ) {
                    return Ok(true);
                }
            }
            Ok(false)
        }
        MolBlockRecord::Query(_) => Ok(false),
    }
}

fn promote_record_to_query(record: &mut MolBlockRecord) -> Result<(), MolPostError> {
    if record_requires_query(record)? {
        let old = std::mem::replace(
            record,
            MolBlockRecord::Concrete {
                topology: TopologyBlock::default(),
                coordinates: CoordinateBlock::default(),
                properties: cosmolkit_model::MoleculeProperties::default(),
            },
        );
        let MolBlockRecord::Concrete {
            topology,
            coordinates,
            properties,
        } = old
        else {
            unreachable!()
        };
        *record = MolBlockRecord::Query(concrete_to_query(topology, coordinates, properties)?);
    }
    Ok(())
}

fn process_smarts_groups(record: &mut MolBlockRecord) -> Result<(), MolPostError> {
    // RDKit✔️❌: void processSMARTSQ(RWMol &mol, const SubstanceGroup &sg) {
    // RDKit✔️❌:   std::string field;
    // RDKit✔️❌:   if (sg.getPropIfPresent("QUERYOP", field) && field != "=") {
    // RDKit✔️❌:     BOOST_LOG(rdWarningLog) << "unrecognized QUERYOP '" << field
    // RDKit✔️❌:                             << "' for SMARTSQ. Query ignored." << std::endl;
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   std::vector<std::string> dataFields;
    // RDKit✔️❌:   if (!sg.getPropIfPresent("DATAFIELDS", dataFields) || dataFields.empty()) {
    // RDKit✔️❌:     BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:         << "empty FIELDDATA for SMARTSQ. Query ignored." << std::endl;
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (dataFields.size() > 1) {
    // RDKit✔️❌:     BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:         << "multiple FIELDDATA values for SMARTSQ. Taking the first."
    // RDKit✔️❌:         << std::endl;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   const std::string &sma = dataFields[0];
    // RDKit✔️❌:   if (sma.empty()) {
    // RDKit✔️❌:     BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:         << "Skipping empty SMARTS value for SMARTSQ." << std::endl;
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto aidx : sg.getAtoms()) {
    // RDKit✔️❌:     auto at = mol.getAtomWithIdx(aidx);
    // RDKit✔️❌:
    // RDKit✔️❌:     std::unique_ptr<RWMol> m;
    // RDKit✔️❌:     try {
    // RDKit✔️❌:       m.reset(SmartsToMol(sma));
    // RDKit✔️❌:     } catch (...) {
    // RDKit✔️❌:       // Is this ever used?
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (!m || !m->getNumAtoms()) {
    // RDKit✔️❌:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:           << "SMARTS for SMARTSQ '" << sma
    // RDKit✔️❌:           << "' could not be parsed or has no atoms. Ignoring it." << std::endl;
    // RDKit✔️❌:       return;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (!at->hasQuery()) {
    // RDKit✔️❌:       QueryAtom qAt(*at);
    // RDKit✔️❌:       int oidx = at->getIdx();
    // RDKit✔️❌:       mol.replaceAtom(oidx, &qAt);
    // RDKit✔️❌:       at = mol.getAtomWithIdx(oidx);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     QueryAtom::QUERYATOM_QUERY *query = nullptr;
    // RDKit✔️❌:     if (m->getNumAtoms() == 1) {
    // RDKit✔️❌:       query = m->getAtomWithIdx(0)->getQuery()->copy();
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       query = new RecursiveStructureQuery(m.release());
    // RDKit✔️❌:     }
    // RDKit✔️❌:     at->setQuery(query);
    // RDKit✔️❌:     at->setProp(common_properties::MRV_SMA, sma);
    // RDKit✔️❌:     at->setProp(common_properties::_MolFileAtomQuery, 1);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // Behavior review: checked source atom lookup, per-atom parsing and
    // source-only parser failure fallback precede query replacement, then
    // MRV_SMA bytes and integer _MolFileAtomQuery are written in that order.
    // Complexity review: the borrowed QueryGraph API requires a group/data
    // snapshot before mutable atom writes; this deep clone exceeds RDKit's
    // reference iteration. Each actual atom invokes the canonical parser in
    // source order; no atom means no parser call. The snapshot cost remains.

    let MolBlockRecord::Query(query_record) = record else {
        return Ok(());
    };
    let mut groups = query_substance_groups(&query_record.query).to_vec();
    let mut remove = vec![false; groups.len()];
    for (index, group) in groups.iter().enumerate() {
        if !source_sgroup_is_data(group)? {
            continue;
        }
        let query_type = source_sgroup_string(
            group,
            "QUERYTYPE",
            group.data().and_then(|data| data.query_type.as_ref()),
        )?;
        if !matches!(
            query_type.as_ref().map(PropertyText::as_bytes),
            Some(b"SMARTSQ" | b"SQ")
        ) {
            continue;
        }
        remove[index] = true;
        let query_op = source_sgroup_string(
            group,
            "QUERYOP",
            group.data().and_then(|data| data.query_op.as_ref()),
        )?;
        if query_op.as_ref().is_some_and(|op| op.as_bytes() != b"=") {
            continue;
        }
        let Some(smarts) = data_values(group)?
            .first()
            .filter(|value| !value.is_empty())
        else {
            continue;
        };
        for atom_id in group.atoms() {
            // Source getAtomWithIdx precedes parsing; invalid members fail here.
            if query_record.query.atom(atom_id.index()).is_none() {
                return Err(MolPostError::Representation(
                    "SMARTSQ atom index out of range",
                ));
            }
            // The source catches parser failures only, then returns from this
            // SGroup. Empty membership never calls the canonical parser.
            let Ok(parsed) = cosmolkit_search::parse_smarts(
                smarts,
                &cosmolkit_search::SmartsParseParams::default(),
            ) else {
                break;
            };
            if parsed.num_atoms() == 0 {
                break;
            }
            let predicate = if parsed.num_atoms() == 1 {
                parsed.atoms()[0].predicate().clone()
            } else {
                QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(
                    RecursiveStructureQuery::from_query_graph(parsed, 0)
                        .with_source_smarts(smarts.clone()),
                ))
            };
            let atom = query_record.query.atom_mut(atom_id.index()).ok_or(
                MolPostError::Representation("SMARTSQ atom index out of range"),
            )?;
            atom.set_predicate(predicate);
            atom.set_prop("MRV SMA", smarts.clone())
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            atom.set_prop("_MolFileAtomQuery", PropertyValue::Int(1))
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        }
    }
    groups = retain_substance_groups(groups, &remove)?;
    replace_query_substance_groups(&mut query_record.query, groups)
        .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
    Ok(())
}

fn attachment_values(topology: &TopologyBlock) -> Result<Vec<Option<i32>>, MolPostError> {
    // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:441-450
    // RDKit❗✔️: template <>
    // RDKit❗✔️: inline int rdvalue_cast<int>(RDValue_cast_t v) {
    // RDKit❗✔️:   if (rdvalue_is<int>(v)) {
    // RDKit❗✔️:     return v.value.i;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (rdvalue_is<unsigned int>(v)) {
    // RDKit❗✔️:     return boost::numeric_cast<int>(v.value.u);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   throw std::bad_any_cast();
    // RDKit❗✔️: }
    // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:441-450

    topology
        .atoms
        .iter()
        .map(|atom| {
            atom.prop("molAttachPoint")
                .map(|value| {
                    if let PropertyValue::UInt(number) = value {
                        return i32::try_from(*number).map_err(|_| {
                            MolPostError::UnsignedPropertyOverflow {
                                atom: atom.id(),
                                property: "molAttachPoint",
                                value: *number,
                            }
                        });
                    }
                    parse_int_property(value).map_err(|()| MolPostError::AttachmentValue {
                        atom: atom.id(),
                        value: property_diagnostic(value),
                    })
                })
                .transpose()
        })
        .collect()
}

fn expand_record_attachment_points(record: MolBlockRecord) -> Result<MolBlockRecord, MolPostError> {
    // BEGIN RDKIT CPP FUNCTION finishMolProcessing
    // RDKit❗❌:   res->clearAllAtomBookmarks();
    // RDKit❗❌:   res->clearAllBondBookmarks();
    // RDKit❗❌:
    // RDKit❗❌:   if (params.expandAttachmentPoints) {
    // RDKit❗❌:     MolOps::expandAttachmentPoints(*res);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // calculate explicit valence on each atom:
    // RDKit❗❌:   for (auto atom : res->atoms()) {
    // RDKit❗❌:     atom->calcExplicitValence(false);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION finishMolProcessing
    // Behavior review: this dispatcher runs before the existing property,
    // stereo and sanitize pipeline. The source logger is projected to stderr;
    // CK-COORD-001 intentionally isolates mixed-conformer degree-one direction.
    // Complexity review: conversion to detached query rows and validated block
    // reconstruction add allocations compared with RDKit's mutable RWMol.
    match record {
        MolBlockRecord::Concrete {
            topology,
            coordinates,
            properties,
        } => {
            let values = attachment_values(&topology)?;
            if !values
                .iter()
                .flatten()
                .any(|value| matches!(value, 1 | 2 | -1))
            {
                let result =
                    expand_attachment_points(topology, coordinates, None, &values, true, true)
                        .map_err(|error| {
                            MolPostError::AttachmentExpansion(MolProcessingError::from(error))
                        })?;
                for warning in result.warnings {
                    eprintln!(
                        "Invalid value for molAttachPoint: {} on atom {}. Not expanding this atttachment point.",
                        warning.value,
                        warning.atom.index()
                    );
                }
                return Ok(MolBlockRecord::Concrete {
                    topology: result.topology,
                    coordinates: result.coordinates,
                    properties,
                });
            }
            let query = concrete_to_query(topology, coordinates, properties)?;
            expand_record_attachment_points(MolBlockRecord::Query(query))
        }
        MolBlockRecord::Query(query_record) => {
            let old_atoms = query_record.query.atoms().to_vec();
            let old_bonds = query_record.query.bonds().to_vec();
            let topology = TopologyBlock::try_from_parts(
                old_atoms
                    .iter()
                    .map(QueryAtom::try_to_atom)
                    .collect::<Result<Vec<_>, _>>()?,
                old_bonds.iter().map(|bond| bond.bond().clone()).collect(),
                query_substance_groups(&query_record.query).to_vec(),
                query_record.query.stereo_groups().to_vec(),
            )
            .map_err(|error| MolPostError::AttachmentExpansion(MolProcessingError::from(error)))?;
            let values = attachment_values(&topology)?;
            let state = QueryStateRef::try_for_topology(&old_atoms, &old_bonds, &topology)
                .map_err(|error| {
                    MolPostError::AttachmentExpansion(MolProcessingError::from(error))
                })?;
            let result = expand_attachment_points(
                topology,
                query_record
                    .query
                    .coordinate_block(query_record.source_coordinate_dim),
                Some(state),
                &values,
                true,
                true,
            )
            .map_err(|error| MolPostError::AttachmentExpansion(MolProcessingError::from(error)))?;
            for warning in result.warnings {
                eprintln!(
                    "Invalid value for molAttachPoint: {} on atom {}. Not expanding this atttachment point.",
                    warning.value,
                    warning.atom.index()
                );
            }
            let (atoms, bonds) = result.query_rows.ok_or(MolPostError::Representation(
                "attachment query rows missing after expansion",
            ))?;
            let substance_groups = result.topology.substance_groups;
            let stereo_groups = result.topology.stereo_groups;
            let mut query = QueryGraph::from_parts(
                atoms,
                bonds,
                query_record.query.props().clone(),
                result.coordinates.conformers_2d,
                result.coordinates.conformers_3d,
                stereo_groups,
            )
            .map_err(|error| MolPostError::AttachmentExpansion(MolProcessingError::from(error)))?;
            query
                .set_source_conformer_order(result.coordinates.source_conformer_order)
                .map_err(|error| {
                    MolPostError::AttachmentExpansion(MolProcessingError::from(error))
                })?;
            replace_query_substance_groups(&mut query, substance_groups).map_err(|error| {
                MolPostError::AttachmentExpansion(MolProcessingError::from(error))
            })?;
            Ok(MolBlockRecord::Query(QueryMolBlockRecord {
                query,
                properties: query_record.properties,
                source_coordinate_dim: result.coordinates.source_coordinate_dim,
            }))
        }
    }
}

fn calculate_record_explicit_valence(record: &mut MolBlockRecord) -> Result<(), MolPostError> {
    // BEGIN RDKIT CPP FUNCTION finishMolProcessing explicit-valence prepass
    // RDKit✔️❌:   // calculate explicit valence on each atom:
    // RDKit✔️❌:   for (auto atom : res->atoms()) {
    // RDKit✔️❌:     atom->calcExplicitValence(false);
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION
    // Behavior review: the core valence owner checks each current carrier row
    // in source order with strict=false and writes the signed-byte result
    // back, preserving the independently stored implicit valence.
    // RDKit✔️✔️: int Atom::calcExplicitValence(bool strict) {
    // RDKit✔️✔️:   bool checkIt = false;
    // RDKit✔️✔️:   d_explicitValence = calculateExplicitValence(*this, strict, checkIt);
    // RDKit✔️✔️:   return d_explicitValence;
    // RDKit✔️✔️: }
    // Complexity review: Concrete borrows its topology, while Query must
    // materialize a validated topology from its carrier rows. The latter adds
    // a full O(V+E) clone/allocation not present in RWMol's in-place pass.
    let mut query_topology = match record {
        MolBlockRecord::Concrete { .. } => None,
        MolBlockRecord::Query(query) => Some(
            TopologyBlock::try_from_parts(
                query
                    .query
                    .atoms()
                    .iter()
                    .map(QueryAtom::try_to_atom)
                    .collect::<Result<Vec<_>, _>>()?,
                query
                    .query
                    .bonds()
                    .iter()
                    .map(|row| row.bond().clone())
                    .collect(),
                query_substance_groups(&query.query).to_vec(),
                query.query.stereo_groups().to_vec(),
            )
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?,
        ),
    };
    let topology = match record {
        MolBlockRecord::Concrete { topology, .. } => topology,
        MolBlockRecord::Query(_) => query_topology
            .as_mut()
            .expect("query topology was materialized"),
    };
    for index in 0..topology.atoms.len() {
        let explicit = calculate_explicit_valence_for_topology(
            topology,
            topology.atoms[index].id(),
            false,
            false,
        )
        .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
        let mut facts = topology.atoms[index].source_valence_facts();
        facts.explicit_valence = explicit as i8;
        topology.atoms[index].set_source_valence_facts(facts);
    }
    if let MolBlockRecord::Query(query) = record {
        for (atom, carrier) in query.query.atoms_mut().iter_mut().zip(
            &query_topology
                .expect("query topology was materialized")
                .atoms,
        ) {
            atom.set_source_valence_facts(carrier.source_valence_facts());
        }
    }
    Ok(())
}

/// Detached final chemistry assignments for the final postprocessed topology.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MolPostDerivedState {
    pub valence: Option<cosmolkit_core::ValenceAssignment>,
    pub rings: Option<cosmolkit_core::RingInfo>,
}

/// Apply the ordered Molfile postprocessing closure to a detached record.
/// `final_state` receives moved final assignments; no live cache is installed.
pub fn finish_mol_block_record(
    mut record: MolBlockRecord,
    chirality_possible: bool,
    params: MolPostParams,
    final_state: Option<&mut MolPostDerivedState>,
) -> Result<MolBlockRecord, MolPostError> {
    if params.expand_attachment_points {
        record = expand_record_attachment_points(record)?;
    }
    calculate_record_explicit_valence(&mut record)?;
    promote_record_to_query(&mut record)?;
    match record {
        MolBlockRecord::Concrete {
            mut topology,
            coordinates,
            properties,
        } => {
            process_atom_properties(&mut topology, None)?;
            process_groups_on_topology(&mut topology, None)?;
            let (topology, coordinates, properties, _, _, state) = apply_stereo_and_sanitize(
                topology,
                coordinates,
                properties,
                chirality_possible,
                params,
                None,
            )?;
            if let Some(output) = final_state {
                *output = state;
            }
            Ok(MolBlockRecord::Concrete {
                topology,
                coordinates,
                properties,
            })
        }
        MolBlockRecord::Query(mut query_record) => {
            let mut topology = TopologyBlock {
                atoms: query_record
                    .query
                    .atoms()
                    .iter()
                    .map(QueryAtom::try_to_atom)
                    .collect::<Result<Vec<_>, _>>()?,
                bonds: query_record
                    .query
                    .bonds()
                    .iter()
                    .map(|bond| bond.bond().clone())
                    .collect(),
                adjacency: AdjacencyList::from_topology(
                    query_record.query.num_atoms(),
                    &query_record
                        .query
                        .bonds()
                        .iter()
                        .map(|bond| bond.bond().clone())
                        .collect::<Vec<_>>(),
                ),
                substance_groups: query_substance_groups(&query_record.query).to_vec(),
                stereo_groups: query_record.query.stereo_groups().to_vec(),
            };
            process_atom_properties(&mut topology, Some(query_record.query.atoms_mut()))?;
            process_groups_on_topology(&mut topology, Some(query_record.query.bonds_mut()))?;
            for (query_atom, atom) in query_record
                .query
                .atoms_mut()
                .iter_mut()
                .zip(&topology.atoms)
            {
                synchronize_query_atom(query_atom, atom.clone());
            }
            for (query_bond, bond) in query_record
                .query
                .bonds_mut()
                .iter_mut()
                .zip(&topology.bonds)
            {
                *query_bond.bond_mut() = bond.clone();
            }
            replace_query_substance_groups(
                &mut query_record.query,
                topology.substance_groups.clone(),
            )
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            let mut wrapped = MolBlockRecord::Query(query_record);
            process_smarts_groups(&mut wrapped)?;
            let MolBlockRecord::Query(mut query_record) = wrapped else {
                unreachable!()
            };
            for (index, query_atom) in query_record.query.atoms().iter().enumerate() {
                topology.atoms[index] = query_atom.try_to_atom()?;
            }
            topology.substance_groups = query_substance_groups(&query_record.query).to_vec();

            let coordinates = query_record
                .query
                .coordinate_block(query_record.source_coordinate_dim);
            let mut old_query = query_record.query;
            let query_props = old_query.props().clone();
            let query_state =
                QueryStateRef::try_for_topology(old_query.atoms(), old_query.bonds(), &topology)
                    .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            let (topology, coordinates, properties, mapping, query_rows, state) =
                apply_stereo_and_sanitize(
                    topology,
                    coordinates,
                    query_record.properties,
                    chirality_possible,
                    params,
                    Some(query_state),
                )?;
            let (mut query_atoms, mut query_bonds) = query_rows.ok_or(
                MolPostError::Representation("query state missing after mol-post finalization"),
            )?;
            // RDKit✔️✔️:   ProcessMolProps(res);
            // RDKit✔️✔️:       MolOps::removeHs(*res);
            // finishMolProcessing mutates the owned RWMol; it does not copy
            // surviving QueryAtom queries. RecursiveStructureQuery::copy is
            // observably different: quickCopy clears its nested properties
            // and conformers. Move the original predicates onto the mapped
            // carriers, retaining the authoritative row mapping and origins.
            // O(V+E) pointer/vector swaps; no recursive predicate allocation.
            for (query, old) in query_atoms.iter_mut().zip(mapping.atoms().new_to_old()) {
                if let Some(old) = old
                    && is_source_query_atom(query)
                {
                    std::mem::swap(
                        query.predicate_mut(),
                        old_query.atoms_mut()[old.index()].predicate_mut(),
                    );
                }
            }
            for (query, old) in query_bonds.iter_mut().zip(mapping.bonds().new_to_old()) {
                if let Some(old) = old
                    && is_source_query_bond(query)
                {
                    std::mem::swap(
                        query.predicate_mut(),
                        old_query.bonds_mut()[old.index()].predicate_mut(),
                    );
                }
            }
            let query_atoms = query_atoms
                .into_iter()
                .zip(&topology.atoms)
                .map(|(mut query, carrier)| {
                    synchronize_query_atom(&mut query, carrier.clone());
                    query
                })
                .collect();
            let query_bonds = query_bonds
                .into_iter()
                .zip(&topology.bonds)
                .map(|(query, carrier)| synchronize_query_bond(query, carrier.clone()))
                .collect();
            let source_coordinate_dim = coordinates.source_coordinate_dim;
            let mut query = QueryGraph::from_parts(
                query_atoms,
                query_bonds,
                query_props,
                coordinates.conformers_2d,
                coordinates.conformers_3d,
                topology.stereo_groups,
            )
            .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            query
                .set_source_conformer_order(coordinates.source_conformer_order)
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            replace_query_substance_groups(&mut query, topology.substance_groups)
                .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
            query_record.query = query;
            query_record.properties = properties;
            query_record.source_coordinate_dim = source_coordinate_dim;
            if query_record.query.prop("_NeedsQueryScan").is_some() {
                query_record
                    .query
                    .clear_prop("_NeedsQueryScan")
                    .map_err(|error| MolPostError::Processing(MolProcessingError::from(error)))?;
                cosmolkit_search::complete_mol_queries(
                    &mut query_record.query,
                    cosmolkit_search::QUERY_SCAN_MAGIC_VALUE,
                );
            }
            if let Some(output) = final_state {
                *output = state;
            }
            Ok(MolBlockRecord::Query(query_record))
        }
    }
}

#[cfg(test)]
mod uint_post_proposed_tests {
    use super::*;
    #[test]
    fn proposed_uint_signed_post_keys_and_attachment_overflow_category() {
        for (value, expected) in [
            (0_u32, Some(0)),
            (1, Some(1)),
            (2147483646, Some(2147483646)),
            (2147483647, Some(2147483647)),
            (2147483648, None),
            (4294967295, None),
        ] {
            for key in ["molSubstCount", "molTotValence", "molAttachPoint"] {
                let atom = cosmolkit_model::Atom::from_spec(
                    AtomId::new(0),
                    cosmolkit_model::AtomSpec::new(cosmolkit_types::Element::C)
                        .with_prop(key, PropertyValue::UInt(value))
                        .unwrap(),
                );
                let wanted = expected.ok_or(MolPostError::UnsignedPropertyOverflow {
                    atom: AtomId::new(0),
                    property: key,
                    value,
                });
                assert_eq!(source_int_property_or_zero(&atom, key), wanted);
                if key == "molAttachPoint" {
                    let graph =
                        TopologyBlock::try_from_parts(vec![atom], vec![], vec![], vec![]).unwrap();
                    assert_eq!(
                        attachment_values(&graph),
                        expected.map(|x| vec![Some(x)]).ok_or(
                            MolPostError::UnsignedPropertyOverflow {
                                atom: AtomId::new(0),
                                property: key,
                                value
                            }
                        )
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec};
    use cosmolkit_types::Element;

    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molSubstCount_0
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molsubstcount_0() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molSubstCount", PropertyValue::UInt(0_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(source_int_property_or_zero(&a, "molSubstCount"), Ok(0));
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molTotValence_0
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_moltotvalence_0() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molTotValence", PropertyValue::UInt(0_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(source_int_property_or_zero(&a, "molTotValence"), Ok(0));
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molAttachPoint_0
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molattachpoint_0() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molAttachPoint", PropertyValue::UInt(0_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(source_int_property_or_zero(&a, "molAttachPoint"), Ok(0));
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molSubstCount_1
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molsubstcount_1() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molSubstCount", PropertyValue::UInt(1_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(source_int_property_or_zero(&a, "molSubstCount"), Ok(1));
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molTotValence_1
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_moltotvalence_1() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molTotValence", PropertyValue::UInt(1_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(source_int_property_or_zero(&a, "molTotValence"), Ok(1));
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molAttachPoint_1
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molattachpoint_1() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molAttachPoint", PropertyValue::UInt(1_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(source_int_property_or_zero(&a, "molAttachPoint"), Ok(1));
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molSubstCount_2147483646
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molsubstcount_2147483646() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molSubstCount", PropertyValue::UInt(2147483646_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molSubstCount"),
            Ok(2147483646)
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molTotValence_2147483646
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_moltotvalence_2147483646() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molTotValence", PropertyValue::UInt(2147483646_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molTotValence"),
            Ok(2147483646)
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molAttachPoint_2147483646
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molattachpoint_2147483646() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molAttachPoint", PropertyValue::UInt(2147483646_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molAttachPoint"),
            Ok(2147483646)
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molSubstCount_2147483647
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molsubstcount_2147483647() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molSubstCount", PropertyValue::UInt(2147483647_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molSubstCount"),
            Ok(2147483647)
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molTotValence_2147483647
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_moltotvalence_2147483647() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molTotValence", PropertyValue::UInt(2147483647_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molTotValence"),
            Ok(2147483647)
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molAttachPoint_2147483647
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molattachpoint_2147483647() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molAttachPoint", PropertyValue::UInt(2147483647_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molAttachPoint"),
            Ok(2147483647)
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molSubstCount_2147483648
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molsubstcount_2147483648() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molSubstCount", PropertyValue::UInt(2147483648_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molSubstCount"),
            Err(MolPostError::UnsignedPropertyOverflow {
                atom: AtomId::new(0),
                property: "molSubstCount",
                value: 2147483648_u32
            })
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molTotValence_2147483648
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_moltotvalence_2147483648() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molTotValence", PropertyValue::UInt(2147483648_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molTotValence"),
            Err(MolPostError::UnsignedPropertyOverflow {
                atom: AtomId::new(0),
                property: "molTotValence",
                value: 2147483648_u32
            })
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molAttachPoint_2147483648
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molattachpoint_2147483648() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molAttachPoint", PropertyValue::UInt(2147483648_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molAttachPoint"),
            Err(MolPostError::UnsignedPropertyOverflow {
                atom: AtomId::new(0),
                property: "molAttachPoint",
                value: 2147483648_u32
            })
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molSubstCount_4294967295
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molsubstcount_4294967295() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molSubstCount", PropertyValue::UInt(4294967295_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molSubstCount"),
            Err(MolPostError::UnsignedPropertyOverflow {
                atom: AtomId::new(0),
                property: "molSubstCount",
                value: 4294967295_u32
            })
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molTotValence_4294967295
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_moltotvalence_4294967295() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molTotValence", PropertyValue::UInt(4294967295_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molTotValence"),
            Err(MolPostError::UnsignedPropertyOverflow {
                atom: AtomId::new(0),
                property: "molTotValence",
                value: 4294967295_u32
            })
        );
        assert_eq!(a, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_io/mol_post_molAttachPoint_4294967295
    #[test]
    fn uint_cell_signed_consumer_io_mol_post_molattachpoint_4294967295() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molAttachPoint", PropertyValue::UInt(4294967295_u32))
                .unwrap(),
        );
        let before = a.clone();
        assert_eq!(
            source_int_property_or_zero(&a, "molAttachPoint"),
            Err(MolPostError::UnsignedPropertyOverflow {
                atom: AtomId::new(0),
                property: "molAttachPoint",
                value: 4294967295_u32
            })
        );
        assert_eq!(a, before);
    }
}
