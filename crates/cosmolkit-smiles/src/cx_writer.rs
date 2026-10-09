#[cfg(test)]
mod native_source_conformer_tests {
    use super::*;
    use cosmolkit_model::{Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension};

    #[test]
    fn source_default_selects_first_of_multiple_homogeneous_sets_without_id_sorting() {
        let two = CoordinateBlock {
            conformers_2d: vec![
                Conformer2D::new(91, vec![[9.0, 1.0]]),
                Conformer2D::new(0, vec![[0.0, 0.0]]),
            ],
            ..Default::default()
        };
        let before = two.clone();
        match source_conformer_by_id(&two, -1).unwrap() {
            CoordinateSource::TwoD(conf) => assert!(std::ptr::eq(conf, &two.conformers_2d[0])),
            CoordinateSource::ThreeD(_) => panic!("wrong stored conformer"),
        }
        assert_eq!(two, before);
        let three = CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(91, vec![[9.0, 1.0, 2.0]], false),
                Conformer3D::new(0, vec![[0.0, 0.0, 0.0]], true),
            ],
            ..Default::default()
        };
        let before = three.clone();
        match source_conformer_by_id(&three, -1).unwrap() {
            CoordinateSource::ThreeD(conf) => assert!(std::ptr::eq(conf, &three.conformers_3d[0])),
            CoordinateSource::TwoD(_) => panic!("wrong stored conformer"),
        }
        assert_eq!(three, before);
    }

    #[test]
    fn source_default_consumes_actual_cross_dimension_order_and_preserves_explicit_choices() {
        for first in [CoordinateDimension::ThreeD, CoordinateDimension::TwoD] {
            let order = if first == CoordinateDimension::ThreeD {
                vec![
                    CoordinateDimension::ThreeD,
                    CoordinateDimension::TwoD,
                    CoordinateDimension::ThreeD,
                ]
            } else {
                vec![
                    CoordinateDimension::TwoD,
                    CoordinateDimension::ThreeD,
                    CoordinateDimension::ThreeD,
                ]
            };
            let coords = CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(1, vec![[1.0, 2.0]])],
                conformers_3d: vec![
                    Conformer3D::new(91, vec![[9.0, 1.0, 2.0]], false),
                    Conformer3D::new(0, vec![[0.0, 0.0, 3.0]], true),
                ],
                source_conformer_order: Some(order),
                ..Default::default()
            };
            let before = coords.clone();
            match (first, source_conformer_by_id(&coords, -1).unwrap()) {
                (CoordinateDimension::ThreeD, CoordinateSource::ThreeD(conf)) => {
                    assert!(std::ptr::eq(conf, &coords.conformers_3d[0]))
                }
                (CoordinateDimension::TwoD, CoordinateSource::TwoD(conf)) => {
                    assert!(std::ptr::eq(conf, &coords.conformers_2d[0]))
                }
                _ => panic!("source insertion order not retained"),
            }
            match select_cx_coordinates(&coords, CxCoordinateSelection::ThreeD { id: 0 })
                .unwrap()
                .unwrap()
            {
                CoordinateSource::ThreeD(conf) => {
                    assert!(std::ptr::eq(conf, &coords.conformers_3d[1]))
                }
                _ => panic!("explicit dimension/ID changed"),
            }
            match select_cx_coordinates(&coords, CxCoordinateSelection::TwoD { id: 1 })
                .unwrap()
                .unwrap()
            {
                CoordinateSource::TwoD(conf) => {
                    assert!(std::ptr::eq(conf, &coords.conformers_2d[0]))
                }
                _ => panic!("explicit dimension/ID changed"),
            }
            assert_eq!(coords, before);
        }
    }
}

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Write as _,
};

use cosmolkit_core::stereo_graph::StereoGraphAccess;
use cosmolkit_core::{
    AtropisomerConformer, CrossedBondContext, RingInfo, RingSearchParams, ValenceAssignment,
    ValenceModel, WedgeAssignments, WedgeInfo, assign_valence_with_options_for_topology,
    atropisomer_carriers_for_bonds, fast_find_rings_from_parts, find_sssr,
    get_all_atom_ids_for_stereo_groups, get_molfile_bond_stereo_info,
    pick_bonds_to_wedge_with_ring_info, property_value_to_string, wedge_bonds_from_atropisomers,
};
use cosmolkit_model::{
    AtomId, Bond, BondId, Conformer2D, Conformer3D, PropertyText, PropertyValue, SGroupConnection,
    StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupKind, ordered_atom_properties,
    set_stereo_group_write_id, stereo_group_write_id,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo};

use crate::{SmilesParseError, SmilesRecord, SmilesWriteParams, writer::write_smiles_for_cx};

fn source_string_property(
    value: Option<&PropertyValue>,
    _name: &str,
) -> Result<Option<cosmolkit_model::PropertyText>, SmilesParseError> {
    // RDKit✔️🔝: bool getValIfPresent(const std::string_view what, std::string &res) const {
    // RDKit✔️🔝:     for (const auto &i : _data) {
    // RDKit✔️🔝:       if (i.key == what) {
    // RDKit✔️🔝:         rdvalue_tostring(i.val, res);
    // RDKit✔️🔝:         return true;
    // RDKit✔️🔝:       }
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     return false;
    // RDKit✔️🔝:   }
    // Behavior: MODEL's canonical lookup supplies absence separately from a
    // present tagged value. Absence returns false/None without formatting.
    // Every present value reaches the single independently compared SF379
    // rdvalue_tostring owner; empty String is Some(empty), never absent. No
    // as_string/UTF8 checks, second formatter, silent error-to-None or fallback.
    // Complexity: source lookup is replaced by the canonical O(log P) tree
    // lookup at callers; one owning converted byte buffer for a present value.
    value
        .map(property_value_to_string)
        .transpose()
        .map_err(SmilesParseError::WriterProperty)
}

fn source_unsigned_property(value: &PropertyValue, name: &str) -> Result<u32, SmilesParseError> {
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
    let _ = name;
    cosmolkit_core::property_value_to_uint(value).map_err(SmilesParseError::WriterNumeric)
}

/// CXSMILES extension families selected for writing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxSmilesFields(u32);

impl CxSmilesFields {
    pub const NONE: Self = Self(0);
    pub const ATOM_LABELS: Self = Self(1 << 0);
    pub const MOLFILE_VALUES: Self = Self(1 << 1);
    pub const COORDS: Self = Self(1 << 2);
    pub const RADICALS: Self = Self(1 << 3);
    pub const ATOM_PROPS: Self = Self(1 << 4);
    pub const LINKNODES: Self = Self(1 << 5);
    pub const ENHANCED_STEREO: Self = Self(1 << 6);
    pub const SGROUPS: Self = Self(1 << 7);
    pub const POLYMER: Self = Self(1 << 8);
    pub const BOND_CFG: Self = Self(1 << 9);
    pub const BOND_ATROPISOMER: Self = Self(1 << 10);
    pub const COORDINATE_BONDS: Self = Self(1 << 11);
    pub const HYDROGEN_BONDS: Self = Self(1 << 12);
    pub const ZERO_BONDS: Self = Self(1 << 13);
    pub const ALL: Self = Self(0x7fff_ffff);
    pub const ALL_BUT_COORDS: Self = Self(Self::ALL.0 ^ Self::COORDS.0);

    #[must_use]
    pub const fn bits(self) -> u32 {
        self.0
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }
}

impl std::ops::BitOr for CxSmilesFields {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

/// Options for detached CXSMILES writing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxSmilesWriteParams {
    pub smiles: SmilesWriteParams,
    pub fields: CxSmilesFields,
    pub coordinate_selection: CxCoordinateSelection,
}

/// Selects which stored coordinate set a CXSMILES export uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CxCoordinateSelection {
    /// Use the first inserted conformer, and omit coordinates when none exist.
    Auto,
    /// Select a stored 2D layout by its dimension-scoped ID.
    TwoD { id: usize },
    /// Select a stored 3D conformer by its dimension-scoped ID.
    ThreeD { id: usize },
}

impl Default for CxSmilesWriteParams {
    fn default() -> Self {
        Self {
            smiles: SmilesWriteParams::default(),
            fields: CxSmilesFields::ALL,
            coordinate_selection: CxCoordinateSelection::Auto,
        }
    }
}

/// Writes canonical CXSMILES with every modeled extension family enabled.
pub fn write_cx_smiles<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
) -> Result<PropertyText, SmilesParseError> {
    let record = record.into();
    write_cx_smiles_with_params(record, &CxSmilesWriteParams::default())
}

/// Writes CXSMILES from detached values with explicit traversal and field policy.
pub fn write_cx_smiles_with_params<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    params: &CxSmilesWriteParams,
) -> Result<PropertyText, SmilesParseError> {
    let record = record.into();
    record
        .coordinates
        .validate_for_atom_count(record.topology.atoms.len())
        .map_err(|error| SmilesParseError::Model(error.to_string()))?;

    // BEGIN RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles copy and kekulization order
    // RDKit❗✔️:   RWMol trwmol(romol);
    // RDKit❗✔️:   SmilesWriteParams params = paramsInput;
    // RDKit❗✔️:   if (params.doKekule) {
    // RDKit❗✔️:     MolOps::Kekulize(trwmol);
    // RDKit❗✔️:     params.doKekule = false;
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles copy and kekulization order
    // The CX wrapper kekulizes its private full input before writer stereo
    // preparation. Keep the caller's detached record unchanged and prevent
    // the later fragment writer from repeating that phase after preparation.
    let mut writer_params = params.smiles;
    let mut prepared_record = record.to_owned_record();
    if writer_params.kekule {
        prepared_record.topology = cosmolkit_core::kekulize(
            &prepared_record.topology,
            &cosmolkit_core::KekulizeParams::default(),
        )
        .map_err(SmilesParseError::WriterKekulize)?
        .topology;
        writer_params.kekule = false;
    }
    // RWMol's private copy above also retains RingInfo. Borrow that source
    // state alongside the prepared values instead of dropping it at the
    // owned-record adapter; no additional state copy is needed.
    let output = write_smiles_for_cx(
        crate::SmilesRecordView {
            topology: &prepared_record.topology,
            coordinates: &prepared_record.coordinates,
            properties: &prepared_record.properties,
            rings: record.rings,
        },
        &writer_params,
    )?;

    // BEGIN RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles empty SMILES return
    // RDKit❗✔️:   if (res.empty()) {
    // RDKit❗✔️:     return res;
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles empty SMILES return
    if output.text.is_empty() {
        return Ok(output.text);
    }

    // BEGIN RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles bond direction restoration
    // RDKit❌❌:   if (restoreBondDirs == RestoreBondDirOptionTrue) {
    // RDKit❌❌:     RDKit::Chirality::reapplyMolBlockWedging(trwmol);
    // RDKit❌❌:   } else if (restoreBondDirs == RestoreBondDirOptionClear) {
    // RDKit❗✔️:     for (auto bond : trwmol.bonds()) {
    // RDKit❗✔️:       if (!canHaveDirection(*bond)) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       using RDKit::common_properties::_MolFileBondCfg;
    // RDKit❗✔️:       if (auto cfg = 0u;
    // RDKit❗✔️:           bond->getPropIfPresent<unsigned int>(_MolFileBondCfg, cfg) &&
    // RDKit❗✔️:           cfg == 2) {
    // RDKit❗✔️:         bond->setBondDir(Bond::BondDir::UNKNOWN);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         if (bond->getBondDir() != Bond::BondDir::NONE) {
    // RDKit❗✔️:           bond->setBondDir(Bond::BondDir::NONE);
    // RDKit❗✔️:         }
    // RDKit❗✔️:         bond->clearProp(_MolFileBondCfg);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles bond direction restoration
    // This API carries RDKit's default Clear choice. Its explicit True choice
    // requires the absent shared reapplyMolBlockWedging owner and is not
    // replaced here by a writer-local chemistry implementation.
    for bond in &mut prepared_record.topology.bonds {
        if !can_have_direction(bond) {
            continue;
        }
        let cfg = bond
            .prop("_MolFileBondCfg")
            .map(|value| source_unsigned_property(value, "_MolFileBondCfg"))
            .transpose()?;
        if cfg == Some(2) {
            bond.set_direction(BondDirection::Unknown);
        } else {
            if bond.direction() != BondDirection::None {
                bond.set_direction(BondDirection::None);
            }
            bond.clear_prop("_MolFileBondCfg")?;
        }
    }

    // BEGIN RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles nonisomeric CX field mask
    // RDKit❗✔️:   if (!params.doIsomericSmiles) {
    // RDKit❗✔️:     flags &= ~(SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO |
    // RDKit❗✔️:                  SmilesWrite::CXSmilesFields::CX_BOND_CFG);
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles nonisomeric CX field mask
    let mut fields = params.fields;
    if !writer_params.isomeric_smiles {
        fields.0 &= !(CxSmilesFields::ENHANCED_STEREO.0 | CxSmilesFields::BOND_CFG.0);
    }

    // The source's outer cleanStereo cache update occurs after the field mask
    // and before assignment/extension writing. The resulting detached valence
    // assignment is consumed by the remaining stage on this same prepared
    // record; extensions still use the maps emitted by the base write.
    let post_base_valence =
        prepare_cx_post_base_valence(&prepared_record, writer_params.clean_stereo)?;
    if let Some(valence) = post_base_valence.as_ref() {
        apply_cx_post_base_stereochemistry(&mut prepared_record, valence)?;
    }

    let extension = write_cx_extensions(
        &prepared_record,
        fields,
        &output.atom_order,
        &output.bond_order,
        post_base_valence,
        params.coordinate_selection,
    )?;
    if extension.is_empty() {
        Ok(output.text)
    } else if output.text.is_empty() {
        Ok(extension)
    } else {
        {
            let mut text = output.text;
            text.push_byte(b' ');
            text.extend_bytes(extension.as_bytes());
            Ok(text)
        }
    }
}

fn prepare_cx_post_base_valence<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    clean_stereo: bool,
) -> Result<Option<cosmolkit_core::ValenceAssignment>, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles cleanStereo outer cache update
    // RDKit❗❌:   if (params.cleanStereo) {
    // RDKit❗❌:     if (trwmol.needsUpdatePropertyCache()) {
    // RDKit❗❌:       trwmol.updatePropertyCache(false);
    // RDKit❗❌:     }
    // END RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles cleanStereo outer cache update
    if !clean_stereo {
        return Ok(None);
    }

    // Behavior review: SmilesRecord has no source valence-cache values or
    // validity bit, so every supported nonempty record reaching this stage
    // maps the source predicate to needs-update. Keep the non-strict owner
    // assignment available for the exact later cleanStereo call; a property
    // string or computed-property marker does not imply cache validity.
    // Complexity review: the detached owner validates topology and allocates
    // row-aligned valence vectors before its atom scan, while ROMol updates
    // atom caches in place and then scans bonds. This preserves the required
    // detached assignment shape but has extra allocation and validation work.
    let valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &record.topology,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?;
    Ok(Some(valence))
}

fn apply_cx_post_base_stereochemistry(
    prepared_record: &mut SmilesRecord,
    valence: &cosmolkit_core::ValenceAssignment,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles cleanStereo assignment and cleanup
    // RDKit❗❌:   if (params.cleanStereo) {
    // RDKit❗❌:     MolOps::assignStereochemistry(trwmol, true);
    // RDKit❗❌:     Chirality::cleanupStereoGroups(trwmol);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION SmilesWrite.cpp::MolToCXSmiles cleanStereo assignment and cleanup
    // BEGIN RDKIT CPP FUNCTION Chirality.cpp::assignStereochemistry done-property effects
    // RDKit❗❌:   if (!force && mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   mol.setProp(common_properties::_StereochemDone, 1, true);
    // END RDKIT CPP FUNCTION Chirality.cpp::assignStereochemistry done-property effects
    // Behavior review: property presence is the force=false skip guard. An
    // absent marker runs the fixed legacy owner with clean=true and
    // possible=false, clears the legacy pending marker, and stores the
    // computed done marker. Explicit group cleanup remains outside that guard.
    // Complexity review: ring and valence assignments use row-aligned
    // detached storage while RDKit updates its molecule in place. This moves
    // the existing topology through the owner without another record clone,
    // but retains the additional linear storage cost.
    if prepared_record.properties.prop("_StereochemDone").is_none() {
        // BEGIN RDKIT CPP FUNCTION Chirality.cpp::legacyStereoPerception molecule property effect
        // RDKit❗❌: void legacyStereoPerception(ROMol &mol, bool cleanIt,
        // RDKit❗❌:                             bool flagPossibleStereoCenters) {
        // RDKit❗❌:   mol.clearProp("_needsDetectBondStereo");
        // END RDKIT CPP FUNCTION Chirality.cpp::legacyStereoPerception molecule property effect
        prepared_record
            .properties
            .clear_prop("_needsDetectBondStereo")?;
        let rings = fast_find_rings_from_parts(
            prepared_record.topology.atoms.len(),
            &prepared_record.topology.bonds,
            &prepared_record.topology.adjacency,
        )
        .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        let topology = std::mem::take(&mut prepared_record.topology);
        prepared_record.topology = cosmolkit_core::assign_legacy_stereochemistry_with_flags(
            topology, valence, &rings, true, false,
        )
        .map_err(|error| match error {
            cosmolkit_core::LegacyStereoError::StereoGroup(cause) => {
                SmilesParseError::StereoGroup(cause)
            }
            other => SmilesParseError::WriterStereo(other.to_string()),
        })?;
        prepared_record
            .properties
            .set_computed_prop("_StereochemDone", 1_i32)
            .map_err(|error| SmilesParseError::Model(error.to_string()))?;
    }
    cosmolkit_core::cleanup_stereo_groups(&mut prepared_record.topology)?;
    Ok(())
}

/// Append counted CX bytes using the sole source delimiter rule.
#[doc(hidden)]
pub fn append_cx_extension_source(addition: impl AsRef<[u8]>, output: &mut PropertyText) {
    // RDKit❗✔️: void appendToCXExtension(const std::string &addition, std::string &base) {
    // RDKit❗✔️:   if (!addition.empty()) {
    // RDKit❗✔️:     if (base.size() > 1) {
    // RDKit❗✔️:       base += ",";
    // RDKit❗✔️:     }
    // RDKit❗✔️:     base += addition;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let addition = addition.as_ref();
    if !addition.is_empty() {
        if output.len() > 1 {
            output.push_byte(b',');
        }
        output.extend_bytes(addition);
    }
    // O(addition bytes), one growing destination, no scan/clone/UTF8 parsing;
    // Rust lacks Native SSO but no material asymptotic/buffering difference.
}

fn append_extension(addition: impl AsRef<[u8]>, output: &mut PropertyText) {
    append_cx_extension_source(addition, output);
}

pub(crate) fn write_cx_extensions<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    fields: CxSmilesFields,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    post_base_valence: Option<ValenceAssignment>,
    coordinate_selection: CxCoordinateSelection,
) -> Result<PropertyText, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::getCXExtensions field order
    // RDKit✔️✔️: if ((flags & SmilesWrite::CXSmilesFields::CX_COORDS) &&
    // RDKit✔️✔️:     mol.getNumConformers()) {
    // RDKit✔️✔️:   res += "(" + get_coords_block(mol, atomOrder) + ")";
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if ((flags & SmilesWrite::CXSmilesFields::CX_ATOM_LABELS) && needLabels) {
    // RDKit✔️✔️:   auto lbls = get_atomlabel_block(mol, atomOrder);
    // RDKit✔️✔️:   if (!lbls.empty()) {
    // RDKit✔️✔️:     if (res.size() > 1) { res += ","; }
    // RDKit✔️✔️:     res += "$" + lbls + "$";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit❗✔️: if ((flags & SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES) && needValues) {
    // RDKit❗✔️:   if (res.size() > 1) { res += ","; }
    // RDKit❗✔️:   res += "$_AV:" +
    // RDKit❗✔️:          get_value_block(mol, atomOrder, common_properties::molFileValue) + "$";
    // RDKit❗✔️: }
    // RDKit✔️✔️: auto radblock = get_radical_block(mol, atomOrder);
    // RDKit✔️✔️: if ((flags & SmilesWrite::CXSmilesFields::CX_RADICALS) && radblock.size()) {
    // RDKit✔️✔️:   if (res.size() > 1) { res += ","; }
    // RDKit✔️✔️:   res += radblock;
    // RDKit✔️✔️:   if (res.back() == ',') { res.erase(res.size() - 1); }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_ATOM_PROPS) {
    // RDKit✔️✔️:   const auto atomblock = get_atom_props_block(mol, atomOrder);
    // RDKit✔️✔️:   appendToCXExtension(atomblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_BOND_CFG) {
    // RDKit✔️✔️:   bool includeCoords = flags & SmilesWrite::CX_COORDS && mol.getNumConformers();
    // RDKit✔️✔️:   const auto cfgblock = get_bond_config_block(mol, atomOrder, bondOrder,
    // RDKit✔️✔️:                                               includeCoords, wedgeBonds);
    // RDKit✔️✔️:   appendToCXExtension(cfgblock, res);
    // RDKit✔️✔️:   const auto cistransblock = get_ringbond_cistrans_block(mol, atomOrder, bondOrder);
    // RDKit✔️✔️:   appendToCXExtension(cistransblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_COORDINATE_BONDS) {
    // RDKit✔️✔️:   const auto block = get_coord_or_hydrogen_bonds_block(
    // RDKit✔️✔️:       mol, Bond::BondType::DATIVE, "C", atomOrder, bondOrder);
    // RDKit✔️✔️:   appendToCXExtension(block, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_HYDROGEN_BONDS) {
    // RDKit✔️✔️:   const auto block = get_coord_or_hydrogen_bonds_block(
    // RDKit✔️✔️:       mol, Bond::BondType::HYDROGEN, "H", atomOrder, bondOrder);
    // RDKit✔️✔️:   appendToCXExtension(block, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_ZERO_BONDS) {
    // RDKit✔️✔️:   const auto block = get_zerobonds_block(mol, atomOrder, bondOrder);
    // RDKit✔️✔️:   appendToCXExtension(block, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_LINKNODES) {
    // RDKit✔️✔️:   const auto linknodeblock = get_linknodes_block(mol, atomOrder);
    // RDKit✔️✔️:   appendToCXExtension(linknodeblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO) {
    // RDKit✔️✔️:   const auto stereoblock = get_enhanced_stereo_block(mol, atomOrder, wedgeBonds);
    // RDKit✔️✔️:   appendToCXExtension(stereoblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_SGROUPS) {
    // RDKit✔️✔️:   const auto sgroupdatablock = get_sgroup_data_block(mol, atomOrder);
    // RDKit✔️✔️:   appendToCXExtension(sgroupdatablock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_POLYMER) {
    // RDKit✔️✔️:   const auto sgrouppolyblock = get_sgroup_polymer_block(mol, atomOrder, bondOrder);
    // RDKit✔️✔️:   appendToCXExtension(sgrouppolyblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & (SmilesWrite::CXSmilesFields::CX_SGROUPS |
    // RDKit✔️✔️:              SmilesWrite::CXSmilesFields::CX_POLYMER)) {
    // RDKit✔️✔️:   const auto sgrouphierarchyblock = get_sgroup_hierarchy_block(mol);
    // RDKit✔️✔️:   appendToCXExtension(sgrouphierarchyblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (res.size() > 1) { res += "|"; } else { res = ""; }
    // END RDKIT CPP FUNCTION SmilesWrite::getCXExtensions field order
    let mut result = PropertyText::from("|");

    // BEGIN RDKIT CPP FUNCTION getCXExtensions field-data preflight
    // RDKit✔️✔️:   bool needLabels = false;
    // RDKit✔️✔️:   bool needValues = false;
    // RDKit✔️✔️:   for (auto idx : atomOrder) {
    // RDKit✔️✔️:     const auto at = mol.getAtomWithIdx(idx);
    // RDKit✔️✔️:     if (at->hasProp(common_properties::atomLabel) ||
    // RDKit✔️✔️:         at->hasProp(common_properties::_QueryAtomGenericLabel) ||
    // RDKit✔️✔️:         at->hasProp(common_properties::dummyLabel) ||
    // RDKit✔️✔️:         at->hasProp(common_properties::_fromAttachPoint)) {
    // RDKit✔️✔️:       needLabels = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (at->hasProp(common_properties::molFileValue)) {
    // RDKit✔️✔️:       needValues = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION getCXExtensions field-data preflight
    // Behavior review: field-presence checks use the base writer's atom map
    // and gate only the matching optional block; selected block order and
    // empty-block separator behavior remain source-shaped.
    // Complexity review: one O(A) source-shaped preflight avoids two
    // atom-order-sized temporary String vectors when labels/values are absent
    // under the default all-fields mask; active emitters retain their order.
    let mut need_labels = false;
    let mut need_values = false;
    for atom_id in atom_order {
        let atom = &record.topology.atoms[atom_id.index()];
        if atom.prop("atomLabel").is_some()
            || atom.prop("_QueryAtomGenericLabel").is_some()
            || atom.prop("dummyLabel").is_some()
            || atom.prop("_fromAttchpt").is_some()
        {
            need_labels = true;
        }
        if atom.prop("molFileValue").is_some() {
            need_values = true;
        }
    }

    let selected_source = if fields.contains(CxSmilesFields::COORDS) {
        coordinate_source(record, coordinate_selection)?
    } else {
        None
    };
    if let Some(source) = selected_source {
        let coords = write_cx_coordinates_from_source(source, atom_order)?;
        result.push_byte(b'(');
        result.extend_bytes((&coords).as_ref());
        result.push_byte(b')');
    }
    if fields.contains(CxSmilesFields::ATOM_LABELS) && need_labels {
        let labels = write_atom_labels(record, atom_order)?;
        if !labels.is_empty() {
            let mut framed = PropertyText::from("$");
            framed.extend_bytes(labels.as_bytes());
            framed.push_byte(b'$');
            append_extension(framed, &mut result);
        }
    }
    if fields.contains(CxSmilesFields::MOLFILE_VALUES) && need_values {
        let values = write_atom_values(record, atom_order)?;
        let mut framed = PropertyText::from("$_AV:");
        framed.extend_bytes(values.as_bytes());
        framed.push_byte(b'$');
        append_extension(framed, &mut result);
    }
    if fields.contains(CxSmilesFields::RADICALS) {
        append_extension(write_radicals(record, atom_order)?, &mut result);
    }
    if fields.contains(CxSmilesFields::ATOM_PROPS) {
        append_extension(write_atom_properties(record, atom_order)?, &mut result);
    }

    // BEGIN RDKIT CPP FUNCTION SmilesWrite::getCXExtensions conformer selection
    // RDKit❗✔️:   const Conformer *conf = nullptr;
    // RDKit❗✔️:   if (mol.getNumConformers() && (flags & SmilesWrite::CXSmilesFields::CX_COORDS)) {
    // RDKit❗✔️:     conf = &mol.getConformer();
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION SmilesWrite::getCXExtensions conformer selection
    // Coordinate output and wedge geometry use the same selected record. The
    // source does not pass stored conformers to wedge helpers when CX_COORDS
    // is disabled.
    let conformer = selected_source.map(|source| match source {
        CoordinateSource::ThreeD(value) => AtropisomerConformer::ThreeD(value),
        CoordinateSource::TwoD(value) => AtropisomerConformer::TwoD(value),
    });
    let coords_included = conformer.is_some();
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::getCXExtensions shared wedge map
    // RDKit❗❌:   std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>> wedgeBonds;
    // RDKit❗❌:   if (flags & SmilesWrite::CXSmilesFields::CX_BOND_CFG) {
    // RDKit❗❌:     wedgeBonds = Chirality::pickBondsToWedge(mol, nullptr, conf);
    // RDKit❗❌:   } else if (flags & SmilesWrite::CXSmilesFields::CX_BOND_ATROPISOMER) {
    // RDKit❗❌:     Atropisomers::wedgeBondsFromAtropisomers(mol, conf, wedgeBonds);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION SmilesWrite::getCXExtensions shared wedge map
    // Behavior review: CX_BOND_CFG takes precedence; the selected detached
    // map remains alive for the later enhanced-stereo source phase.
    // Complexity review: the map and SSSR are produced once and shared by
    // extension consumers; the core helper returns its promoted ring state.
    let wedge_assignments;
    if fields.contains(CxSmilesFields::BOND_CFG) {
        let (assignments, ring_info) =
            pick_bonds_to_wedge_with_ring_info(&record.topology, conformer)
                .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        wedge_assignments = assignments;

        let wedge_valence = if coords_included {
            Some(match post_base_valence {
                Some(valence) => valence,
                None => assign_valence_with_options_for_topology(
                    &record.topology,
                    ValenceModel::RdkitLike,
                    false,
                )
                .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?,
            })
        } else {
            None
        };
        let crossed_bonds = if let Some(valence) = wedge_valence.as_ref() {
            Some(
                CrossedBondContext::new(&record.topology, valence, &ring_info, true)
                    .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?,
            )
        } else {
            None
        };
        append_extension(
            write_bond_config(
                record,
                atom_order,
                bond_order,
                coords_included,
                false,
                &wedge_assignments,
                crossed_bonds.as_ref(),
                conformer,
            )?,
            &mut result,
        );
        append_extension(
            write_ring_bond_stereo(record, atom_order, bond_order, &ring_info),
            &mut result,
        );
    } else if fields.contains(CxSmilesFields::BOND_ATROPISOMER) {
        let ring_info = find_sssr(&record.topology, &RingSearchParams::default())
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        let atropisomer_assignments = wedge_bonds_from_atropisomers(
            &record.topology,
            &ring_info,
            conformer,
            &BTreeSet::new(),
        )
        .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        wedge_assignments =
            WedgeAssignments::from_atropisomer_wedge_assignment(atropisomer_assignments);
        append_extension(
            write_bond_config(
                record,
                atom_order,
                bond_order,
                coords_included,
                true,
                &wedge_assignments,
                None,
                conformer,
            )?,
            &mut result,
        );
    } else {
        wedge_assignments = WedgeAssignments::default();
    }
    if fields.contains(CxSmilesFields::COORDINATE_BONDS) {
        append_extension(
            write_coordinate_bonds(record, atom_order, bond_order, BondOrder::Dative, "C")?,
            &mut result,
        );
    }
    if fields.contains(CxSmilesFields::HYDROGEN_BONDS) {
        append_extension(
            write_coordinate_bonds(record, atom_order, bond_order, BondOrder::Hydrogen, "H")?,
            &mut result,
        );
    }
    if fields.contains(CxSmilesFields::ZERO_BONDS) {
        append_extension(write_zero_bonds(record, bond_order)?, &mut result);
    }
    if fields.contains(CxSmilesFields::LINKNODES) {
        append_extension(write_link_nodes(record, atom_order)?, &mut result);
    }
    if fields.contains(CxSmilesFields::ENHANCED_STEREO) {
        append_extension(
            write_enhanced_stereo(record, atom_order, &wedge_assignments)?,
            &mut result,
        );
    }
    if fields.contains(CxSmilesFields::SGROUPS) {
        append_extension(write_data_sgroups(record, atom_order)?, &mut result);
    }
    if fields.contains(CxSmilesFields::POLYMER) {
        append_extension(
            write_polymer_sgroups(record, atom_order, bond_order)?,
            &mut result,
        );
    }
    if fields.contains(CxSmilesFields::SGROUPS) || fields.contains(CxSmilesFields::POLYMER) {
        append_extension(
            write_sgroup_hierarchy(
                record,
                fields.contains(CxSmilesFields::SGROUPS),
                fields.contains(CxSmilesFields::POLYMER),
            )?,
            &mut result,
        );
    }

    if result.len() == 1 {
        Ok(PropertyText::new())
    } else {
        result.push_byte(b'|');
        Ok(result)
    }
}

fn atom_positions(atom_order: &[AtomId], atom_count: usize) -> Vec<Option<usize>> {
    let mut positions = vec![None; atom_count];
    for (position, atom) in atom_order.iter().copied().enumerate() {
        positions[atom.index()] = Some(position);
    }
    positions
}

#[doc(hidden)]
pub fn zero_small_cx_coordinate(value: f64) -> f64 {
    // BEGIN RDKIT CPP FUNCTION CXSmilesOps.cpp::zero_small_vals
    // RDKit✔️✔️: double zero_small_vals(double val) {
    // RDKit✔️✔️:   if (fabs(val) < 1e-4) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return val;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION CXSmilesOps.cpp::zero_small_vals
    // Behavior review: the strict comparison preserves both signed values at
    // exactly 1e-4 and returns positive zero below it, on either side.
    // Complexity review: one abs, comparison and branch; constant time, no
    // allocation, equivalent to the source fabs and conditional.
    if value.abs() < 1e-4 { 0.0 } else { value }
}

// BEGIN RDKIT CPP FUNCTION CXSmilesOps.cpp::get_coords_block general conversion
// RDKit❗✔️:     res += boost::str(boost::format("%g,%g,") % zero_small_vals(pt.x) %
// RDKit❗✔️:                       zero_small_vals(pt.y));
// RDKit❗✔️:       auto zc = boost::str(boost::format("%g") % zero_small_vals(pt.z));
// END RDKIT CPP FUNCTION CXSmilesOps.cpp::get_coords_block general conversion
// Behavior review: `%g` defaults to six significant digits, switches to
// scientific notation when the rounded exponent is below -4 or at least 6,
// trims trailing zeroes, and preserves the sign of zero.
// Complexity review: one formatted String allocation and constant work over
// six digits; exponent parsing and in-place conversion keep this comparable
// to the pinned single `%g` formatting result.
fn format_general(value: f64) -> String {
    // RDKit❗❌:     res += boost::str(boost::format("%g,%g,") % zero_small_vals(pt.x) %
    // RDKit❗❌:                       zero_small_vals(pt.y));
    // RDKit❗❌:       auto zc = boost::str(boost::format("%g") % zero_small_vals(pt.z));
    if !value.is_finite() {
        // Pinned lower-case %g and the canonical CORE Boost conversion share
        // inf/-inf/nan/-nan spelling and sign. Double is an infallible modeled
        // source tag; do not use Rust Display's uppercase, unsigned NaN.
        let text = property_value_to_string(&PropertyValue::Double(value))
            .expect("source Double-to-string conversion is infallible");
        return String::from_utf8(text.into_bytes())
            .expect("canonical numeric string conversion emits ASCII");
    }
    if value == 0.0 {
        return if value.is_sign_negative() {
            "-0".to_owned()
        } else {
            "0".to_owned()
        };
    }

    let mut text = format!("{value:.5e}");
    let exponent_at = text
        .find('e')
        .expect("finite scientific formatting includes an exponent");
    let exponent = text[exponent_at + 1..]
        .parse::<i32>()
        .expect("scientific exponent is a signed decimal integer");

    if !(-4..6).contains(&exponent) {
        text.truncate(exponent_at);
        while text.ends_with('0') {
            text.pop();
        }
        if text.ends_with('.') {
            text.pop();
        }
        std::fmt::Write::write_fmt(&mut text, format_args!("e{exponent:+03}"))
            .expect("writing a formatted exponent to String cannot fail");
        return text;
    }

    let negative = text.starts_with('-');
    let mut digits = [0_u8; 6];
    let mut digit_count = 0;
    for byte in text[..exponent_at].bytes() {
        if byte.is_ascii_digit() {
            digits[digit_count] = byte;
            digit_count += 1;
        }
    }
    debug_assert_eq!(digit_count, digits.len());

    text.clear();
    if negative {
        text.push('-');
    }
    let decimal_position = exponent + 1;
    if decimal_position <= 0 {
        text.push_str("0.");
        for _ in 0..-decimal_position {
            text.push('0');
        }
        for digit in digits {
            text.push(char::from(digit));
        }
    } else {
        let decimal_position =
            usize::try_from(decimal_position).expect("a non-negative decimal position fits usize");
        if decimal_position >= digits.len() {
            for digit in digits {
                text.push(char::from(digit));
            }
            for _ in digits.len()..decimal_position {
                text.push('0');
            }
        } else {
            for digit in &digits[..decimal_position] {
                text.push(char::from(*digit));
            }
            text.push('.');
            for digit in &digits[decimal_position..] {
                text.push(char::from(*digit));
            }
        }
    }
    if text.contains('.') {
        while text.ends_with('0') {
            text.pop();
        }
        if text.ends_with('.') {
            text.pop();
        }
    }
    text
}

/// The canonical SOURCE coordinate number conversion, shared by query output.
#[doc(hidden)]
pub fn format_cx_coordinate(value: f64) -> String {
    format_general(zero_small_cx_coordinate(value))
}

#[derive(Clone, Copy)]
#[doc(hidden)]
pub enum CoordinateSource<'a> {
    ThreeD(&'a Conformer3D),
    TwoD(&'a Conformer2D),
}

fn coordinate_source<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    selection: CxCoordinateSelection,
) -> Result<Option<CoordinateSource<'record>>, SmilesParseError> {
    select_cx_coordinates(record.into().coordinates, selection)
}

#[doc(hidden)]
pub fn select_cx_coordinates<'record>(
    coordinates: &'record cosmolkit_model::CoordinateBlock,
    selection: CxCoordinateSelection,
) -> Result<Option<CoordinateSource<'record>>, SmilesParseError> {
    select_cx_coordinates_from_sets(
        &coordinates.conformers_2d,
        &coordinates.conformers_3d,
        selection,
    )
}

/// Borrow canonical split conformer storage using the same sole selector.
#[doc(hidden)]
pub fn select_cx_coordinates_from_sets<'record>(
    conformers_2d: &'record [Conformer2D],
    conformers_3d: &'record [Conformer3D],
    selection: CxCoordinateSelection,
) -> Result<Option<CoordinateSource<'record>>, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION CXSmilesOps.cpp::get_coords_block conformer selection
    // RDKit❗✔️:   const auto &conf = mol.getConformer();
    // END RDKIT CPP FUNCTION CXSmilesOps.cpp::get_coords_block conformer selection
    // BEGIN RDKIT CPP FUNCTION ROMol.cpp::getConformer(int id) const
    // RDKit❗✔️: const Conformer &ROMol::getConformer(int id) const {
    // RDKit❗✔️:   // make sure we have more than one conformation
    // RDKit❗✔️:   if (d_confs.size() == 0) {
    // RDKit❗✔️:     throw ConformerException("No conformations available on the molecule");
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (id < 0) {
    // RDKit❗✔️:     return *(d_confs.front());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto cid = (unsigned int)id;
    // RDKit❗✔️:   for (auto conf : d_confs) {
    // RDKit❗✔️:     if (conf->getId() == cid) {
    // RDKit❗✔️:       return *conf;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // we did not find a conformation with the specified ID
    // RDKit❗✔️:   std::string mesg = "Can't find conformation with ID: ";
    // RDKit❗✔️:   mesg += id;
    // RDKit❗✔️:   throw ConformerException(mesg);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION ROMol.cpp::getConformer(int id) const
    // CK-COORD-002 (approved): Auto selects no set for zero, the unique set
    // for one, and a typed ambiguity for multiple. Explicit selection matches
    // dimension plus the stored ID; it never infers from vector position,
    // smallest ID, source dimension, or Conformer3D::is_3d(). This is the
    // deliberate default-selection difference from RDKit's first insertion.
    // Behavior review: source's implicit first conformer is intentionally
    // replaced only at selection; formatting and geometry still consume this
    // exact borrowed set.
    // Complexity review: Auto uses collection lengths and a direct singleton
    // access (O(1)); explicit lookup scans only the selected dimension (O(n)),
    // matching ROMol's explicit-ID scan and without allocating a new set.
    match selection {
        CxCoordinateSelection::Auto => {
            let two_d_count = conformers_2d.len();
            let three_d_count = conformers_3d.len();
            match (two_d_count, three_d_count) {
                (0, 0) => Ok(None),
                (1, 0) => Ok(Some(CoordinateSource::TwoD(&conformers_2d[0]))),
                (0, 1) => Ok(Some(CoordinateSource::ThreeD(&conformers_3d[0]))),
                (two_d_count, three_d_count) => {
                    Err(SmilesParseError::AmbiguousCoordinateSelection {
                        two_d_count,
                        three_d_count,
                    })
                }
            }
        }
        CxCoordinateSelection::TwoD { id } => {
            let conformer = conformers_2d
                .iter()
                .find(|conformer| conformer.id() == id)
                .ok_or(SmilesParseError::MissingCoordinateSelection { selection })?;
            Ok(Some(CoordinateSource::TwoD(conformer)))
        }
        CxCoordinateSelection::ThreeD { id } => {
            let conformer = conformers_3d
                .iter()
                .find(|conformer| conformer.id() == id)
                .ok_or(SmilesParseError::MissingCoordinateSelection { selection })?;
            Ok(Some(CoordinateSource::ThreeD(conformer)))
        }
    }
}

#[derive(Clone, Copy)]
enum SourceConformerIndex {
    TwoD(usize),
    ThreeD(usize),
}

fn source_conformer_index_by_id(
    coordinates: &cosmolkit_model::CoordinateBlock,
    id: i32,
) -> Result<SourceConformerIndex, SmilesParseError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION ROMol.cpp::getConformer(int id) const
    // RDKit✔️✔️: const Conformer &ROMol::getConformer(int id) const {
    // RDKit✔️✔️:   // make sure we have more than one conformation
    // RDKit✔️✔️:   if (d_confs.size() == 0) {
    // RDKit✔️✔️:     throw ConformerException("No conformations available on the molecule");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (id < 0) {
    // RDKit✔️✔️:     return *(d_confs.front());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto cid = (unsigned int)id;
    // RDKit✔️✔️:   for (auto conf : d_confs) {
    // RDKit✔️✔️:     if (conf->getId() == cid) {
    // RDKit✔️✔️:       return *conf;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // we did not find a conformation with the specified ID
    // RDKit✔️✔️:   std::string mesg = "Can't find conformation with ID: ";
    // RDKit✔️✔️:   mesg += id;
    // RDKit✔️✔️:   throw ConformerException(mesg);
    // RDKit✔️✔️: }
    // END RDKIT COMPLETE CPP FUNCTION ROMol.cpp::getConformer(int id) const
    // MODEL supplies the genuine source front; indexes identify that same
    // existing row for both immutable and mutable borrowed projections below.
    // O(1) default / O(C) ID scan, constant counters and no success allocation.
    // No pointer reconstruction, guessed ordering or second lookup algorithm.
    // Native mesg += id appends a narrowed char: decimal error context below
    // is project-structured exception translation, not source message bytes.
    if coordinates.conformers_2d.is_empty() && coordinates.conformers_3d.is_empty() {
        return Err(SmilesParseError::Model(
            "No conformations available on the molecule".into(),
        ));
    }
    if id < 0 {
        let first = coordinates
            .first_source_conformer()
            .map_err(SmilesParseError::Coordinates)?
            .ok_or(SmilesParseError::Coordinates(
                cosmolkit_model::CoordinateValidationError::MissingSourceConformerOrder,
            ))?;
        return Ok(match first {
            cosmolkit_model::CoordinateSourceConformer::TwoD(_) => SourceConformerIndex::TwoD(0),
            cosmolkit_model::CoordinateSourceConformer::ThreeD(_) => {
                SourceConformerIndex::ThreeD(0)
            }
        });
    }
    let cid = id as u32 as usize;
    if let Some(order) = coordinates.source_conformer_order.as_deref() {
        let mut two_d_index = 0;
        let mut three_d_index = 0;
        for dimension in order {
            let (stored_id, index) =
                match dimension {
                    cosmolkit_model::CoordinateDimension::TwoD => {
                        let row = coordinates
                            .conformers_2d
                            .get(two_d_index)
                            .ok_or(SmilesParseError::Coordinates(
                            cosmolkit_model::CoordinateValidationError::MissingSourceConformerOrder,
                        ))?;
                        let index = SourceConformerIndex::TwoD(two_d_index);
                        two_d_index += 1;
                        (row.id(), index)
                    }
                    cosmolkit_model::CoordinateDimension::ThreeD => {
                        let row = coordinates
                            .conformers_3d
                            .get(three_d_index)
                            .ok_or(SmilesParseError::Coordinates(
                            cosmolkit_model::CoordinateValidationError::MissingSourceConformerOrder,
                        ))?;
                        let index = SourceConformerIndex::ThreeD(three_d_index);
                        three_d_index += 1;
                        (row.id(), index)
                    }
                };
            if stored_id == cid {
                return Ok(index);
            }
        }
        if two_d_index != coordinates.conformers_2d.len()
            || three_d_index != coordinates.conformers_3d.len()
        {
            return Err(SmilesParseError::Coordinates(
                cosmolkit_model::CoordinateValidationError::MissingSourceConformerOrder,
            ));
        }
    } else if coordinates.conformers_2d.is_empty() {
        if let Some(index) = coordinates
            .conformers_3d
            .iter()
            .position(|row| row.id() == cid)
        {
            return Ok(SourceConformerIndex::ThreeD(index));
        }
    } else if coordinates.conformers_3d.is_empty() {
        if let Some(index) = coordinates
            .conformers_2d
            .iter()
            .position(|row| row.id() == cid)
        {
            return Ok(SourceConformerIndex::TwoD(index));
        }
    } else {
        return Err(SmilesParseError::Coordinates(
            cosmolkit_model::CoordinateValidationError::MissingSourceConformerOrder,
        ));
    }
    Err(SmilesParseError::Model(format!(
        "Can't find conformation with ID: {id}"
    )))
}

pub(super) fn source_conformer_by_id(
    coordinates: &cosmolkit_model::CoordinateBlock,
    id: i32,
) -> Result<CoordinateSource<'_>, SmilesParseError> {
    Ok(match source_conformer_index_by_id(coordinates, id)? {
        SourceConformerIndex::TwoD(index) => {
            CoordinateSource::TwoD(&coordinates.conformers_2d[index])
        }
        SourceConformerIndex::ThreeD(index) => {
            CoordinateSource::ThreeD(&coordinates.conformers_3d[index])
        }
    })
}

enum CoordinateSourceMut<'a> {
    TwoD(&'a mut Conformer2D),
    ThreeD(&'a mut Conformer3D),
}

fn source_conformer_by_id_mut(
    coordinates: &mut cosmolkit_model::CoordinateBlock,
    id: i32,
) -> Result<CoordinateSourceMut<'_>, SmilesParseError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION ROMol.cpp::getConformer(int id)
    // RDKit✔️✔️: Conformer &ROMol::getConformer(int id) {
    // RDKit✔️✔️:   return const_cast<Conformer &>(
    // RDKit✔️✔️:       static_cast<const ROMol *>(this)->getConformer(id));
    // RDKit✔️✔️: }
    // END RDKIT COMPLETE CPP FUNCTION ROMol.cpp::getConformer(int id)
    // Both overloads delegate to the sole compared const lookup body. The
    // returned mutable reference aliases the same original detached row;
    // no clone, Cow, unrestricted runtime access, or duplicated selection.
    // Lookup O(1)/O(C) then O(1) borrowed projection; no new allocation/scan.
    Ok(match source_conformer_index_by_id(coordinates, id)? {
        SourceConformerIndex::TwoD(index) => {
            CoordinateSourceMut::TwoD(&mut coordinates.conformers_2d[index])
        }
        SourceConformerIndex::ThreeD(index) => {
            CoordinateSourceMut::ThreeD(&mut coordinates.conformers_3d[index])
        }
    })
}

/// Write the canonical CX coordinate body from an actual borrowed conformer.
#[doc(hidden)]
pub fn write_cx_coordinates_from_source(
    source: CoordinateSource<'_>,
    atom_order: &[AtomId],
) -> Result<String, SmilesParseError> {
    // RDKit❗❌: std::string get_coords_block(const ROMol &mol,
    // RDKit❗❌:                              const std::vector<unsigned int> &atomOrder) {
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:   const auto &conf = mol.getConformer();
    // RDKit❗❌:   bool first = true;
    // RDKit❗❌:   for (auto idx : atomOrder) {
    // RDKit❗❌:     const auto &pt = conf.getAtomPos(idx);
    // RDKit❗❌:     if (!first) {
    // RDKit❗❌:       res += ";";
    // RDKit❗❌:     } else {
    // RDKit❗❌:       first = false;
    // RDKit❗❌:     }
    // RDKit❗❌:     res += boost::str(boost::format("%g,%g,") % zero_small_vals(pt.x) %
    // RDKit❗❌:                       zero_small_vals(pt.y));
    // RDKit❗❌:     if (conf.is3D()) {
    // RDKit❗❌:       auto zc = boost::str(boost::format("%g") % zero_small_vals(pt.z));
    // RDKit❗❌:       if (zc != "0") {
    // RDKit❗❌:         res += zc;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // RDKit❗✔️: const RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) const {
    // RDKit❗✔️:   if (dp_mol) {
    // RDKit❗✔️:     PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   URANGE_CHECK(atomId, d_positions.size());
    // RDKit❗✔️:   return d_positions.at(atomId);
    // RDKit❗✔️: }
    // Source getConformer/front selection is performed by the caller's sole
    // MODEL selector. This input is that existing borrowed row, never guessed
    // from an ID, dimension preference, geometry or import-provenance hint.
    let mut output = String::new();
    for (position, atom) in atom_order.iter().enumerate() {
        let (x, y, z) = match source {
            CoordinateSource::ThreeD(conformer) => {
                let point = conformer.coordinates().get(atom.index()).ok_or(
                    SmilesParseError::CxCoordinateAtomOutOfRange {
                        atom: *atom,
                        atom_count: conformer.coordinates().len(),
                    },
                )?;
                (
                    point[0],
                    point[1],
                    if conformer.is_3d() {
                        Some(point[2])
                    } else {
                        None
                    },
                )
            }
            CoordinateSource::TwoD(conformer) => {
                let point = conformer.coordinates().get(atom.index()).ok_or(
                    SmilesParseError::CxCoordinateAtomOutOfRange {
                        atom: *atom,
                        atom_count: conformer.coordinates().len(),
                    },
                )?;
                (point[0], point[1], None)
            }
        };
        if position != 0 {
            output.push(';');
        }
        output.push_str(&format_general(zero_small_cx_coordinate(x)));
        output.push(',');
        output.push_str(&format_general(zero_small_cx_coordinate(y)));
        output.push(',');
        if let Some(z) = z {
            let z = format_general(zero_small_cx_coordinate(z));
            if z != "0" {
                output.push_str(&z);
            }
        }
    }
    // Cost ❌: the shared scalar formatter owns heap strings and reparses its
    // rounded exponent/digits, and output has no SSO. One output buffer avoids
    // the old per-point String vector/join and keeps source member traversal.
    Ok(output)
}

fn write_atom_labels<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
) -> Result<PropertyText, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION get_atomlabel_block
    // RDKit❗🔝: std::string res = "";
    // RDKit❗🔝: for (auto idx : atomOrder) {
    // RDKit❗🔝:   if (idx != atomOrder.front()) {
    // RDKit❗🔝:     res += ";";
    // RDKit❗🔝:   }
    // RDKit❗🔝:   std::string lbl;
    // RDKit❗🔝:   int val;
    // RDKit❗🔝:   const auto atom = mol.getAtomWithIdx(idx);
    // RDKit❗🔝:   if (atom->getPropIfPresent(common_properties::_QueryAtomGenericLabel,
    // RDKit❗🔝:                               lbl)) {
    // RDKit❗🔝:     res += quote_string(lbl + "_p");
    // RDKit❗🔝:   } else if (!atom->getAtomicNum() &&
    // RDKit❗🔝:              atom->getPropIfPresent(common_properties::dummyLabel,
    // RDKit❗🔝:                                     lbl) &&
    // RDKit❗🔝:              std::find(SmilesParseOps::pseudoatoms.begin(),
    // RDKit❗🔝:                        SmilesParseOps::pseudoatoms.end(), lbl) !=
    // RDKit❗🔝:                  SmilesParseOps::pseudoatoms.end()) {
    // RDKit❗🔝:     res += quote_string(lbl + "_p");
    // RDKit❗🔝:   } else if (!atom->getAtomicNum() &&
    // RDKit❗🔝:              atom->getPropIfPresent(common_properties::_fromAttachPoint,
    // RDKit❗🔝:                                         val) &&
    // RDKit❗🔝:              (val == 1 || val == 2)) {
    // RDKit❗🔝:     res += quote_string("_AP" + std::to_string(val));
    // RDKit❗🔝:   } else if (atom->getPropIfPresent(common_properties::atomLabel,
    // RDKit❗🔝:                                         lbl)) {
    // RDKit❗🔝:     res += quote_string(lbl);
    // RDKit❗🔝:   }
    // RDKit❗🔝: }
    // RDKit❗🔝: // if we didn't find anything return an empty string
    // RDKit❗🔝: if (std::find_if_not(res.begin(), res.end(),
    // RDKit❗🔝:                      [](const auto c) { return c == ';'; }) ==
    // RDKit❗🔝:     res.end()) {
    // RDKit❗🔝:   res.clear();
    // RDKit❗🔝: }
    // RDKit❗🔝: return res;
    // END RDKIT CPP FUNCTION get_atomlabel_block
    // BEGIN RDKIT CPP FUNCTION quote_string
    // RDKit❗🔝: std::string quote_string(const std::string &txt) {
    // RDKit❗🔝:   // FIX
    // RDKit❗🔝:   return txt;
    // RDKit❗🔝: }
    // END RDKIT CPP FUNCTION quote_string
    // BEGIN RDKIT CPP FUNCTION SmilesParseOps::pseudoatoms
    // RDKit❗🔝: constexpr std::array<std::string_view, 2> pseudoatoms{"Pol", "Mod"};
    // END RDKIT CPP FUNCTION SmilesParseOps::pseudoatoms
    // Behavior review: branch order, identity quoting, output-map separators
    // and semicolon-only clearing follow the pinned source; keep behavior
    // unresolved until the S40 boundary tests run.
    // Complexity review: one pass appends into one String, then scans its
    // bytes once, matching source O(atoms + output bytes); this avoids the
    // old Vec<String> and per-label formatting allocations and is expected
    // to allocate less on label-heavy molecules.
    const PSEUDOATOMS: [&[u8]; 2] = [b"Pol", b"Mod"];
    let Some(first_atom) = atom_order.first().copied() else {
        return Ok(PropertyText::new());
    };

    let mut labels = PropertyText::new();
    for atom_id in atom_order {
        if *atom_id != first_atom {
            labels.push_byte(b';');
        }
        let atom = &record.topology.atoms[atom_id.index()];
        if let Some(label) = atom.prop("_QueryAtomGenericLabel") {
            let label = source_string_property(Some(label), "_QueryAtomGenericLabel")?
                .expect("present canonical property conversion preserves presence");
            labels.extend_bytes((&label).as_ref());
            labels.extend_bytes(("_p").as_ref());
        } else if atom.atomic_number() == 0 {
            let dummy_label = source_string_property(atom.prop("dummyLabel"), "dummyLabel")?;
            if dummy_label
                .as_ref()
                .is_some_and(|label| PSEUDOATOMS.contains(&label.as_bytes()))
            {
                labels.extend_bytes((&dummy_label.expect("checked above")).as_ref());
                labels.extend_bytes(("_p").as_ref());
            } else if let Some(value) = atom.prop("_fromAttchpt") {
                let value = source_unsigned_property(value, "_fromAttchpt")?;
                if matches!(value, 1 | 2) {
                    labels.extend_bytes(("_AP").as_ref());
                    labels.extend_bytes((&value.to_string()).as_ref());
                } else if let Some(label) = atom.prop("atomLabel") {
                    labels.extend_bytes(
                        (&source_string_property(Some(label), "atomLabel")?
                            .expect("present canonical property conversion preserves presence"))
                            .as_ref(),
                    );
                }
            } else if let Some(label) = atom.prop("atomLabel") {
                labels.extend_bytes(
                    (&source_string_property(Some(label), "atomLabel")?
                        .expect("present canonical property conversion preserves presence"))
                        .as_ref(),
                );
            }
        } else if let Some(label) = atom.prop("atomLabel") {
            labels.extend_bytes(
                (&source_string_property(Some(label), "atomLabel")?
                    .expect("present canonical property conversion preserves presence"))
                    .as_ref(),
            );
        }
    }

    if labels.as_bytes().iter().all(|byte| *byte == b';') {
        Ok(PropertyText::new())
    } else {
        Ok(labels)
    }
}

fn write_atom_values<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
) -> Result<PropertyText, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION get_value_block
    // RDKit❗✔️: std::string res = "";
    // RDKit❗✔️: bool first = true;
    // RDKit❗✔️: for (auto idx : atomOrder) {
    // RDKit❗✔️:   if (!first) {
    // RDKit❗✔️:     res += ";";
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     first = false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::string lbl;
    // RDKit❗✔️:   if (mol.getAtomWithIdx(idx)->getPropIfPresent(prop, lbl)) {
    // RDKit❗✔️:     res += quote_string(lbl);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: return res;
    // END RDKIT CPP FUNCTION get_value_block
    // BEGIN RDKIT CPP FUNCTION quote_string
    // RDKit❗✔️: std::string quote_string(const std::string &txt) {
    // RDKit❗✔️:   // FIX
    // RDKit❗✔️:   return txt;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION quote_string
    // Behavior review: each base output-map slot gets exactly one separator
    // after the first and present values are appended literally, including
    // empty values and punctuation; do not collapse empty slot sequences.
    // Complexity review: one atom-order pass and one output String are
    // source-shaped O(atoms + emitted bytes) with no per-slot String vector.
    let mut values = PropertyText::new();
    let mut first = true;
    for atom_id in atom_order {
        if !first {
            values.push_byte(b';');
        } else {
            first = false;
        }
        if let Some(value) = record.topology.atoms[atom_id.index()].prop("molFileValue") {
            values.extend_bytes(
                (&source_string_property(Some(value), "molFileValue")?
                    .expect("present canonical property conversion preserves presence"))
                    .as_ref(),
            );
        }
    }
    Ok(values)
}

fn write_radicals<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
) -> Result<String, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION get_radical_block
    // RDKit❗🔝: std::string get_radical_block(const ROMol &mol,
    // RDKit❗🔝:                             const std::vector<unsigned int> &atomOrder) {
    // RDKit❗🔝:   std::string res = "";
    // RDKit❗🔝:   std::map<unsigned int, std::vector<unsigned int>> rads;
    // RDKit❗🔝:   for (unsigned int i = 0; i < atomOrder.size(); ++i) {
    // RDKit❗🔝:     auto idx = atomOrder[i];
    // RDKit❗🔝:     auto nrad = mol.getAtomWithIdx(idx)->getNumRadicalElectrons();
    // RDKit❗🔝:     if (nrad) {
    // RDKit❗🔝:       rads[nrad].push_back(i);
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:   if (rads.size()) {
    // RDKit❗🔝:     for (const auto &pr : rads) {
    // RDKit❗🔝:       switch (pr.first) {
    // RDKit❗🔝:         case 1:
    // RDKit❗🔝:           res += "^1:";
    // RDKit❗🔝:           break;
    // RDKit❗🔝:         case 2:
    // RDKit❗🔝:           res += "^2:";
    // RDKit❗🔝:           break;
    // RDKit❗🔝:         case 3:
    // RDKit❗🔝:           res += "^5:";
    // RDKit❗🔝:           break;
    // RDKit❗🔝:         default:
    // RDKit❗🔝:           BOOST_LOG(rdWarningLog) << "unsupported number of radical electrons "
    // RDKit❗🔝:                                   << pr.first << std::endl;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       for (auto aidx : pr.second) {
    // RDKit❗🔝:         res += boost::str(boost::format("%d,") % aidx);
    // RDKit❗🔝:       }
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return res;
    // RDKit❗🔝: }
    // END RDKIT CPP FUNCTION get_radical_block
    // The domain writer preserves the source's serialized default branch
    // (unprefixed indices for counts outside 1..=3); RDKit's Boost warning
    // logger has no corresponding detached-writer channel here. Formatting
    // indices directly into one String avoids the source's per-index temporary
    // boost::format String while keeping the same ordered groups and bytes.
    let mut radicals = BTreeMap::<u8, Vec<usize>>::new();
    for (position, atom) in atom_order.iter().copied().enumerate() {
        let count = record.topology.atoms[atom.index()].radical_electrons();
        if count > 0 {
            radicals.entry(count).or_default().push(position);
        }
    }
    let mut block = String::new();
    for (count, atoms) in radicals {
        match count {
            1 => block.push_str("^1:"),
            2 => block.push_str("^2:"),
            3 => block.push_str("^5:"),
            _ => {} // Source warns, then appends positions without a prefix.
        }
        for position in atoms {
            write!(&mut block, "{position},")
                .expect("formatting an atom output index into String cannot fail");
        }
    }
    if block.ends_with(',') {
        block.pop();
    }
    Ok(block)
}

#[doc(hidden)]
pub fn quote_cx_atom_property(text: &[u8]) -> PropertyText {
    // BEGIN RDKIT CPP FUNCTION quote_atomprop_string
    // RDKit❗❌: std::string quote_atomprop_string(const std::string &txt) {
    // RDKit❗❌:   // at a bare minimum, . needs to be escaped
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:   for (auto c : txt) {
    // RDKit❗❌:     if (c == '.') {
    // RDKit❗❌:       res += "&#46;";
    // RDKit❗❌:     } else {
    // RDKit❗❌:       res += c;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION quote_atomprop_string
    // Behavior review: Rust replaces only ASCII periods with the source token;
    // empty text and every other counted byte, including NUL/non-UTF8, remain unchanged.
    // Complexity review: one linear string replacement creates one output
    // buffer, matching the source's linear scan and output construction.
    // Cost ❌: short Native std::string results use SSO; the counted byte
    // vector allocates even for a single period. The scan/buffer stay linear.
    let mut result = PropertyText::new();
    for &byte in text {
        if byte == b'.' {
            result.extend_bytes(b"&#46;");
        } else {
            result.push_byte(byte);
        }
    }
    result
}

trait SourceCxAtomProperties {
    fn source_atomic_number(&self) -> u8;
    fn source_prop(&self, key: &[u8]) -> Option<&PropertyValue>;
    fn source_prop_required(
        &self,
        key: &PropertyText,
    ) -> Result<&PropertyValue, cosmolkit_model::MissingPropertyError>;
    fn source_property_records(
        &self,
    ) -> Result<
        impl Iterator<Item = (&PropertyText, &PropertyValue)> + '_,
        cosmolkit_model::AtomPropertyError,
    >;
}
impl SourceCxAtomProperties for cosmolkit_model::Atom {
    fn source_atomic_number(&self) -> u8 {
        self.atomic_number()
    }
    fn source_prop(&self, key: &[u8]) -> Option<&PropertyValue> {
        self.prop(key)
    }
    fn source_prop_required(
        &self,
        key: &PropertyText,
    ) -> Result<&PropertyValue, cosmolkit_model::MissingPropertyError> {
        self.prop_required(key)
    }
    fn source_property_records(
        &self,
    ) -> Result<
        impl Iterator<Item = (&PropertyText, &PropertyValue)> + '_,
        cosmolkit_model::AtomPropertyError,
    > {
        self.property_records(false, false)
    }
}
impl SourceCxAtomProperties for cosmolkit_model::QueryAtom {
    fn source_atomic_number(&self) -> u8 {
        self.atomic_number()
    }
    fn source_prop(&self, key: &[u8]) -> Option<&PropertyValue> {
        self.prop(key)
    }
    fn source_prop_required(
        &self,
        key: &PropertyText,
    ) -> Result<&PropertyValue, cosmolkit_model::MissingPropertyError> {
        self.prop_required(key)
    }
    fn source_property_records(
        &self,
    ) -> Result<
        impl Iterator<Item = (&PropertyText, &PropertyValue)> + '_,
        cosmolkit_model::AtomPropertyError,
    > {
        self.property_records(false, false)
    }
}

fn write_atom_properties<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
) -> Result<PropertyText, SmilesParseError> {
    write_cx_atom_properties_impl(&record.into().topology.atoms, atom_order)
}

/// Serialize actual borrowed query-carrier properties through the sole CX body.
#[doc(hidden)]
pub fn write_query_cx_atom_properties_source(
    query: &cosmolkit_model::QueryGraph,
    atom_order: &[AtomId],
) -> Result<PropertyText, SmilesParseError> {
    write_cx_atom_properties_impl(query.atoms(), atom_order)
}

fn write_cx_atom_properties_impl<A: SourceCxAtomProperties>(
    atoms: &[A],
    atom_order: &[AtomId],
) -> Result<PropertyText, SmilesParseError> {
    // RDKit❗❌: std::string get_atom_props_block(const ROMol &mol,
    // RDKit❗❌:                                  const std::vector<unsigned int> &atomOrder) {
    // RDKit❗❌:   constexpr std::array<std::string_view, 7> skip = {
    // RDKit❗❌:       common_properties::atomLabel,       common_properties::molFileValue,
    // RDKit❗❌:       common_properties::molParity,       common_properties::molAtomMapNumber,
    // RDKit❗❌:       common_properties::molStereoCare,   common_properties::molRxnExactChange,
    // RDKit❗❌:       common_properties::molInversionFlag};
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:   unsigned int which = 0;
    // RDKit❗❌:   for (auto idx : atomOrder) {
    // RDKit❗❌:     const auto atom = mol.getAtomWithIdx(idx);
    // RDKit❗❌:     bool isAttachmentPoint = !atom->getAtomicNum() &&
    // RDKit❗❌:                              atom->hasProp(common_properties::_fromAttachPoint);
    // RDKit❗❌:     bool includePrivate = false, includeComputed = false;
    // RDKit❗❌:     for (const auto &pn : atom->getPropList(includePrivate, includeComputed)) {
    // RDKit❗❌:       if (std::find(skip.begin(), skip.end(), pn) == skip.end()) {
    // RDKit❗❌:         std::string pv = atom->getProp<std::string>(pn);
    // RDKit❗❌:         if (pn == "dummyLabel" &&
    // RDKit❗❌:             (isAttachmentPoint || pv == "*" ||
    // RDKit❗❌:              std::find(SmilesParseOps::pseudoatoms.begin(),
    // RDKit❗❌:                        SmilesParseOps::pseudoatoms.end(),
    // RDKit❗❌:                        pv) != SmilesParseOps::pseudoatoms.end())) {
    // RDKit❗❌:           // it's a pseudoatom or attachment point, skip it
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         if (res.empty()) {
    // RDKit❗❌:           res += "atomProp";
    // RDKit❗❌:         }
    // RDKit❗❌:         res +=
    // RDKit❗❌:             boost::str(boost::format(":%d.%s.%s") % which %
    // RDKit❗❌:                        quote_atomprop_string(pn) % quote_atomprop_string(pv));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     ++which;
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // RDKit❗✔️: inline constexpr std::string_view molRxnExactChange = "molRxnExachg";
    const SKIP: [&[u8]; 7] = [
        b"atomLabel",
        b"molFileValue",
        b"molParity",
        b"molAtomMapNumber",
        b"molStereoCare",
        b"molRxnExachg",
        b"molInversionFlag",
    ];
    const PSEUDOATOMS: [&[u8]; 2] = [b"Pol", b"Mod"];
    let mut result = PropertyText::new();
    let mut which = 0u32;
    for atom_id in atom_order {
        let atom =
            atoms
                .get(atom_id.index())
                .ok_or(SmilesParseError::CxAtomPropertyAtomOutOfRange {
                    atom: *atom_id,
                    atom_count: atoms.len(),
                })?;
        let attachment =
            atom.source_atomic_number() == 0 && atom.source_prop(b"_fromAttchpt").is_some();
        for (name, _) in atom.source_property_records().map_err(|source| {
            SmilesParseError::WriterPropertyList {
                atom: *atom_id,
                source,
            }
        })? {
            if SKIP.contains(&name.as_bytes()) {
                continue;
            }
            let value =
                cosmolkit_core::required_property_value_to_string(atom.source_prop_required(name))
                    .map_err(|source| SmilesParseError::WriterRequiredProperty {
                        atom: *atom_id,
                        source,
                    })?;
            if name.as_bytes() == b"dummyLabel"
                && (attachment
                    || value.as_bytes() == b"*"
                    || PSEUDOATOMS.contains(&value.as_bytes()))
            {
                continue;
            }
            if result.is_empty() {
                result.extend_bytes(b"atomProp");
            }
            write!(&mut result, ":{which}.").expect("counted byte formatting is infallible");
            result.extend_bytes(quote_cx_atom_property(name.as_bytes()).as_bytes());
            result.push_byte(b'.');
            result.extend_bytes(quote_cx_atom_property(value.as_bytes()).as_bytes());
        }
        which = which.wrapping_add(1);
    }
    // Cost ❌: canonical scalar/quoted text owns heap buffers without Native
    // SSO. One borrowed property iterator and one output buffer replace the
    // old QUERY entry vector/join; no Atom/QueryAtom/AST/property payload clone.
    // The closed private trait reads actual canonical facts only; one body
    // preserves all source branches for ordinary and query carriers.
    Ok(result)
}

fn can_have_direction(bond: &Bond) -> bool {
    matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic)
}

fn normalized_wedge_direction(bond: &Bond) -> BondDirection {
    match bond.direction() {
        BondDirection::BeginDash | BondDirection::BeginWedge | BondDirection::Unknown => {
            bond.direction()
        }
        _ => BondDirection::None,
    }
}

fn molfile_cfg_bond_direction(cfg: Option<u32>) -> BondDirection {
    // BEGIN RDKIT CPP FUNCTION CXSmilesOps::get_bond_config_block _MolFileBondCfg switch
    // RDKit❗❌:         switch (cfg) {
    // RDKit❗❌:           case 1:
    // RDKit❗❌:             bd = Bond::BondDir::BEGINWEDGE;
    // RDKit❗❌:             break;
    // RDKit❗❌:           case 2:
    // RDKit❗❌:             bd = Bond::BondDir::UNKNOWN;
    // RDKit❗❌:             break;
    // RDKit❗❌:           case 3:
    // RDKit❗❌:             bd = Bond::BondDir::BEGINDASH;
    // RDKit❗❌:             break;
    // RDKit❗❌:
    // RDKit❗❌:           default:
    // RDKit❗❌:             bd = Bond::BondDir::NONE;
    // RDKit❗❌:         }
    // END RDKIT CPP FUNCTION CXSmilesOps::get_bond_config_block _MolFileBondCfg switch
    // Behavior review: this property maps 3 to BEGINDASH and 4 to NONE.
    // Complexity review: one integer match is constant time with no allocation.
    match cfg {
        Some(1) => BondDirection::BeginWedge,
        Some(2) => BondDirection::Unknown,
        Some(3) => BondDirection::BeginDash,
        _ => BondDirection::None,
    }
}

fn cx_writer_direction_from_molfile_code(dir_code: i32) -> BondDirection {
    // BEGIN RDKIT CPP FUNCTION CXSmilesOps::get_bond_config_block GetMolFile code switch
    // RDKit❗❌:         switch (dirCode) {
    // RDKit❗❌:           case 1:
    // RDKit❗❌:             bd = Bond::BondDir::BEGINWEDGE;
    // RDKit❗❌:             break;
    // RDKit❗❌:           case 3:
    // RDKit❗❌:             bd = Bond::BondDir::UNKNOWN;
    // RDKit❗❌:             break;
    // RDKit❗❌:           case 6:
    // RDKit❗❌:             bd = Bond::BondDir::BEGINDASH;
    // RDKit❗❌:             break;
    // RDKit❗❌:           default:
    // RDKit❗❌:             bd = Bond::BondDir::NONE;
    // RDKit❗❌:         }
    // END RDKIT CPP FUNCTION CXSmilesOps::get_bond_config_block GetMolFile code switch
    // Behavior review: GetMolFile code 3 maps to UNKNOWN; code 4 uses NONE.
    // Complexity review: one integer match is constant time with no allocation.
    match dir_code {
        1 => BondDirection::BeginWedge,
        3 => BondDirection::Unknown,
        6 => BondDirection::BeginDash,
        _ => BondDirection::None,
    }
}

fn cx_atom_position(atom_order: &[AtomId], atom: AtomId) -> usize {
    atom_order
        .iter()
        .position(|candidate| *candidate == atom)
        .unwrap_or(atom_order.len())
}

fn other_bond_atom(bond: &Bond, atom: AtomId) -> AtomId {
    if bond.begin() == atom {
        bond.end()
    } else {
        bond.begin()
    }
}

fn flip_atropisomer_wedge_for_output_order<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
    axial_bond: BondId,
    wedge_start_atom: AtomId,
    direction: BondDirection,
    carriers: &[cosmolkit_core::AtropisomerCarrierEnd; 2],
) -> BondDirection {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION CXSmilesOps::get_bond_config_block atropisomer reordering
    // RDKit❗❌:               unsigned int firstReorderedIdx =
    // RDKit❗❌:                   std::find(atomOrder.begin(), atomOrder.end(),
    // RDKit❗❌:                             bondNbr->getBeginAtom()->getIdx()) -
    // RDKit❗❌:                   atomOrder.begin();
    // RDKit❗❌:               unsigned int secondReorderedIdx =
    // RDKit❗❌:                   std::find(atomOrder.begin(), atomOrder.end(),
    // RDKit❗❌:                             bondNbr->getEndAtom()->getIdx()) -
    // RDKit❗❌:                   atomOrder.begin();
    // RDKit❗❌:               if (firstReorderedIdx > secondReorderedIdx) {
    // RDKit❗❌:                 ++swaps;
    // RDKit❗❌:               }
    // RDKit❗❌:               for (unsigned int bondAtomIndex = 0; bondAtomIndex < 2;
    // RDKit❗❌:                    ++bondAtomIndex) {
    // RDKit❗❌:                 if (atomAndBondVecs[bondAtomIndex].first == firstAtom) {
    // RDKit❗❌:                   continue;  // swapped atoms on the side where the wedge bond
    // RDKit❗❌:                              // is does NOT change the wedge bond
    // RDKit❗❌:                 }
    // RDKit❗❌:                 if (atomAndBondVecs[bondAtomIndex].second.size() == 2) {
    // RDKit❗❌:                   unsigned int firstOtherAtomIdx =
    // RDKit❗❌:                       atomAndBondVecs[bondAtomIndex]
    // RDKit❗❌:                           .second[0]
    // RDKit❗❌:                           ->getOtherAtom(atomAndBondVecs[bondAtomIndex].first)
    // RDKit❗❌:                           ->getIdx();
    // RDKit❗❌:                   unsigned int secondOtherAtomIdx =
    // RDKit❗❌:                       atomAndBondVecs[bondAtomIndex]
    // RDKit❗❌:                           .second[1]
    // RDKit❗❌:                           ->getOtherAtom(atomAndBondVecs[bondAtomIndex].first)
    // RDKit❗❌:                           ->getIdx();
    // RDKit❗❌:                   unsigned int firstReorderedAtomIdx =
    // RDKit❗❌:                       std::find(atomOrder.begin(), atomOrder.end(),
    // RDKit❗❌:                                 firstOtherAtomIdx) -
    // RDKit❗❌:                       atomOrder.begin();
    // RDKit❗❌:                   unsigned int secondReorderedAtomIdx =
    // RDKit❗❌:                       std::find(atomOrder.begin(), atomOrder.end(),
    // RDKit❗❌:                                 secondOtherAtomIdx) -
    // RDKit❗❌:                       atomOrder.begin();
    // RDKit❗❌:                   if (firstReorderedAtomIdx > secondReorderedAtomIdx) {
    // RDKit❗❌:                     ++swaps;
    // RDKit❗❌:                   }
    // RDKit❗❌:                 }
    // RDKit❗❌:               }
    // RDKit❗❌:               if (swaps % 2) {
    // RDKit❗❌:                 bd = (bd == Bond::BondDir::BEGINWEDGE)
    // RDKit❗❌:                          ? Bond::BondDir::BEGINDASH
    // RDKit❗❌:                          : Bond::BondDir::BEGINWEDGE;
    // RDKit❗❌:               }
    // END RDKIT CPP FUNCTION CXSmilesOps::get_bond_config_block atropisomer reordering
    // Behavior is still awaiting the fixed CX atropisomer output-order cases.
    // Complexity review: the source performs four `std::find` operations for
    // the axial/carrier atoms; the detached path performs the same lookups and
    // obtains carrier ends through one batch validation instead of per bond.
    let axial = &record.topology.bonds[axial_bond.index()];
    let mut swaps = usize::from(
        cx_atom_position(atom_order, axial.begin()) > cx_atom_position(atom_order, axial.end()),
    );
    for carrier_end in carriers {
        if carrier_end.focus() == wedge_start_atom {
            continue;
        }
        let carrier_bonds = carrier_end.carrier_bonds();
        if carrier_bonds.len() == 2 {
            let focus = carrier_end.focus();
            let first = other_bond_atom(&record.topology.bonds[carrier_bonds[0].index()], focus);
            let second = other_bond_atom(&record.topology.bonds[carrier_bonds[1].index()], focus);
            swaps += usize::from(
                cx_atom_position(atom_order, first) > cx_atom_position(atom_order, second),
            );
        }
    }
    if swaps % 2 == 0 {
        direction
    } else {
        match direction {
            BondDirection::BeginWedge => BondDirection::BeginDash,
            BondDirection::BeginDash => BondDirection::BeginWedge,
            value => value,
        }
    }
}

fn write_bond_config<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    coords_included: bool,
    atropisomer_only: bool,
    wedge_assignments: &WedgeAssignments,
    crossed_bonds: Option<&CrossedBondContext<'_>>,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<String, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION get_bond_config_block core emission
    // RDKit❗❌: for (unsigned int i = 0; i < bondOrder.size(); ++i) {
    // RDKit❗❌:   const auto bond = mol.getBondWithIdx(bondOrder[i]);
    // RDKit❗❌:   unsigned int wedgeStartAtomIdx = bond->getBeginAtomIdx();
    // RDKit❗❌:   if (!canHaveDirection(*bond)) { continue; }
    // RDKit❗❌:   Bond::BondDir bd = bond->getBondDir();
    // RDKit❗❌:   switch (bd) {
    // RDKit❗❌:     case Bond::BondDir::BEGINDASH:
    // RDKit❗❌:     case Bond::BondDir::BEGINWEDGE:
    // RDKit❗❌:     case Bond::BondDir::UNKNOWN:
    // RDKit❗❌:       break;
    // RDKit❗❌:     default:
    // RDKit❗❌:       bd = Bond::BondDir::NONE;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!atropisomerOnly && bd == Bond::BondDir::NONE &&
    // RDKit❗❌:       bond->getPropIfPresent(common_properties::_MolFileBondCfg, cfg)) {
    // RDKit❗❌:     switch (cfg) {
    // RDKit❗❌:       case 1:
    // RDKit❗❌:         bd = Bond::BondDir::BEGINWEDGE;
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 2:
    // RDKit❗❌:         bd = Bond::BondDir::UNKNOWN;
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 3:
    // RDKit❗❌:         bd = Bond::BondDir::BEGINDASH;
    // RDKit❗❌:         break;
    // RDKit❗❌:       default:
    // RDKit❗❌:         bd = Bond::BondDir::NONE;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (bd == Bond::BondDir::UNKNOWN) {
    // RDKit❗❌:     wType = "w";
    // RDKit❗❌:   }
    // RDKit❗❌:   else if (coordsIncluded || isAnAtropisomer) {
    // RDKit❗❌:     if (bd == Bond::BondDir::BEGINWEDGE) {
    // RDKit❗❌:       wType = "wU";
    // RDKit❗❌:     } else if (bd == Bond::BondDir::BEGINDASH) {
    // RDKit❗❌:       wType = "wD";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   wParts[wType].push_back(format("%d.%d", begAtomOrder, i));
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION get_bond_config_block core emission
    // Behavior review: current direction precedes _MolFileBondCfg and the
    // source GetMolFile switch; atrop-only emission keeps its own early gate.
    // Complexity review: output atom mapping is direct-indexed. When there
    // are axial bonds and no coordinates, one validated batch carrier lookup
    // adds an O(V+E) pass before the source-order bond loop.
    let positions = atom_positions(atom_order, record.topology.atoms.len());
    let atropisomer_carriers = if coords_included {
        BTreeMap::new()
    } else {
        let axial_bonds = record
            .topology
            .bonds
            .iter()
            .filter(|bond| matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw))
            .map(Bond::id)
            .collect::<Vec<_>>();
        if axial_bonds.is_empty() {
            BTreeMap::new()
        } else {
            atropisomer_carriers_for_bonds(&record.topology, &axial_bonds)
                .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?
                .into_iter()
                .collect::<BTreeMap<_, _>>()
        }
    };
    let mut parts = BTreeMap::<&'static str, Vec<String>>::new();
    for (bond_position, bond_id) in bond_order.iter().copied().enumerate() {
        let bond = &record.topology.bonds[bond_id.index()];
        if !can_have_direction(bond) {
            continue;
        }
        let mut wedge_start_atom = bond.begin();
        let mut direction = normalized_wedge_direction(bond);

        // BEGIN RDKIT CPP FUNCTION Atropisomers::WedgeBondFromAtropisomerOneBondNoConf state handoff
        // RDKit❗❌:     auto bestBond = atomAndBondVecs[bestBondEnd].second[bestBondNumber];
        // RDKit❗❌:     if (bestBond->getBeginAtom() != atomAndBondVecs[bestBondEnd].first) {
        // RDKit❗❌:       bestBond->setEndAtom(bestBond->getBeginAtom());
        // RDKit❗❌:       bestBond->setBeginAtom(atomAndBondVecs[bestBondEnd].first);
        // RDKit❗❌:     }
        // RDKit❗❌:     bestBond->setBondDir(bestBondDir);
        // END RDKIT CPP FUNCTION Atropisomers::WedgeBondFromAtropisomerOneBondNoConf state handoff
        // The detached owner returns this mutation as a typed map update. Project
        // it into the serialization-local state before source gates inspect the
        // direction or begin atom; the borrowed record remains unchanged.
        // Complexity review: one BTreeMap lookup costs O(log A) per directional
        // bond, versus the source's direct read after its in-place O(1) mutation.
        if let Some(WedgeInfo::Atropisomer { update }) = wedge_assignments.get(bond_id) {
            wedge_start_atom = update.begin;
            direction = update.direction;
        }

        // BEGIN RDKIT CPP FUNCTION CXSmilesOps::get_bond_config_block atrop-only gate
        // RDKit❗❌:   if (atropisomerOnly && bd == Bond::BondDir::NONE) {
        // RDKit❗❌:     continue;
        // RDKit❗❌:   }
        // END RDKIT CPP FUNCTION CXSmilesOps::get_bond_config_block atrop-only gate
        if atropisomer_only && direction == BondDirection::None {
            continue;
        }

        // BEGIN RDKIT CPP FUNCTION CXSmilesOps::get_bond_config_block atropisomer detection
        // RDKit❗❌:   const Atom *firstAtom = bond->getBeginAtom();
        // RDKit❗❌:   if (bd == Bond::BondDir::BEGINDASH || bd == Bond::BondDir::BEGINWEDGE) {
        // RDKit❗❌:     for (auto bondNbr : mol.atomBonds(firstAtom)) {
        // RDKit❗❌:       if (bondNbr->getIdx() == bond->getIdx()) {
        // RDKit❗❌:         continue;  // a bond is not its own neighbor
        // RDKit❗❌:       }
        // RDKit❗❌:       if (bondNbr->getStereo() == Bond::BondStereo::STEREOATROPCW ||
        // RDKit❗❌:           bondNbr->getStereo() == Bond::BondStereo::STEREOATROPCCW) {
        // RDKit❗❌:         isAnAtropisomer = true;
        // RDKit❗❌:         break;
        // RDKit❗❌:       }
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // END RDKIT CPP FUNCTION CXSmilesOps::get_bond_config_block atropisomer detection
        // The source tests the current begin-wedge/dash before reading config
        // or consulting the shared assignment map. Complexity review: inspect
        // the begin atom's incident bonds once and stop at the first axial bond.
        let mut atropisomer_bond = None;
        if matches!(
            direction,
            BondDirection::BeginWedge | BondDirection::BeginDash
        ) {
            for neighbor in record
                .topology
                .adjacency
                .neighbors_of(wedge_start_atom.index())
            {
                if neighbor.bond == bond_id {
                    continue;
                }
                if matches!(
                    record.topology.bonds[neighbor.bond.index()].stereo(),
                    BondStereo::AtropCw | BondStereo::AtropCcw
                ) {
                    atropisomer_bond = Some(neighbor.bond);
                    break;
                }
            }
        }
        let is_atropisomer = atropisomer_bond.is_some();

        if !coords_included
            && matches!(
                direction,
                BondDirection::BeginWedge | BondDirection::BeginDash
            )
            && let Some(axial_bond) = atropisomer_bond
        {
            let Some(Some(carriers)) = atropisomer_carriers.get(&axial_bond) else {
                return Err(SmilesParseError::WriterStereo(
                    "Internal error - should not occur".to_owned(),
                ));
            };
            direction = flip_atropisomer_wedge_for_output_order(
                record,
                atom_order,
                axial_bond,
                wedge_start_atom,
                direction,
                carriers,
            );
        }

        if atropisomer_only && !is_atropisomer {
            continue;
        }

        if !atropisomer_only && direction == BondDirection::None {
            // This property uses its own switch: value 3 means BEGINDASH,
            // while value 4 and every other default branch mean NONE.
            let cfg = bond
                .prop("_MolFileBondCfg")
                .map(|value| source_unsigned_property(value, "_MolFileBondCfg"))
                .transpose()?;
            direction = molfile_cfg_bond_direction(cfg);
        }

        if !atropisomer_only && direction == BondDirection::None && coords_included {
            // BEGIN RDKIT CPP FUNCTION Chirality::GetMolFileBondStereoInfo CX dir-code projection
            // RDKit❗❌:         Chirality::GetMolFileBondStereoInfo(
            // RDKit❗❌:             bond, wedgeBonds, &mol.getConformer(), dirCode, reverse);
            // RDKit❗❌:         switch (dirCode) {
            // RDKit❗❌:           case 1: bd = Bond::BondDir::BEGINWEDGE; break;
            // RDKit❗❌:           case 3: bd = Bond::BondDir::UNKNOWN; break;
            // RDKit❗❌:           case 6: bd = Bond::BondDir::BEGINDASH; break;
            // RDKit❗❌:           default: bd = Bond::BondDir::NONE;
            // RDKit❗❌:         }
            // RDKit❗❌:         if (reverse) {
            // RDKit❗❌:           wedgeStartAtomIdx = bond->getEndAtomIdx();
            // RDKit❗❌:         }
            // END RDKIT CPP FUNCTION Chirality::GetMolFileBondStereoInfo CX dir-code projection
            let crossed_bonds = crossed_bonds.ok_or_else(|| {
                SmilesParseError::WriterStereo(
                    "missing crossed-bond context for coordinate CX bond config".to_owned(),
                )
            })?;
            let info =
                get_molfile_bond_stereo_info(crossed_bonds, wedge_assignments, bond_id, conformer)
                    .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
            // Keep this writer's source switch distinct from `_MolFileBondCfg`:
            // helper code 3 maps to UNKNOWN, code 4 takes the default NONE arm.
            direction = cx_writer_direction_from_molfile_code(info.direction_code);
            if info.reverse {
                wedge_start_atom = bond.end();
            }
        }

        let kind = match direction {
            BondDirection::Unknown => Some("w"),
            BondDirection::BeginWedge if coords_included || is_atropisomer => Some("wU"),
            BondDirection::BeginDash if coords_included || is_atropisomer => Some("wD"),
            _ => None,
        };
        let (Some(kind), Some(begin_position)) = (kind, positions[wedge_start_atom.index()]) else {
            continue;
        };
        parts
            .entry(kind)
            .or_default()
            .push(format!("{begin_position}.{bond_position}"));
    }
    Ok(parts
        .into_iter()
        .map(|(kind, entries)| format!("{kind}:{}", entries.join(",")))
        .collect::<Vec<_>>()
        .join(","))
}

fn write_coordinate_bonds<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    order: BondOrder,
    symbol: &str,
) -> Result<PropertyText, SmilesParseError> {
    write_cx_coord_or_hydrogen_bonds_source(
        cosmolkit_core::stereo_graph::BondRows::Concrete(&record.into().topology.bonds),
        atom_order,
        bond_order,
        order,
        symbol.as_bytes(),
    )
}

/// Source CX typed-bond emitter over the actual borrowed canonical bond rows.
#[doc(hidden)]
pub fn write_cx_coord_or_hydrogen_bonds_source(
    bonds: cosmolkit_core::stereo_graph::BondRows<'_>,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    order: BondOrder,
    symbol: &[u8],
) -> Result<PropertyText, SmilesParseError> {
    // RDKit❗🔝: std::string get_coord_or_hydrogen_bonds_block(
    // RDKit❗🔝:     const ROMol &mol, Bond::BondType bondType, std::string symbol,
    // RDKit❗🔝:     const std::vector<unsigned int> &atomOrder,
    // RDKit❗🔝:     const std::vector<unsigned int> &bondOrder) {
    // RDKit❗🔝:   std::string res = "";
    // RDKit❗🔝:   for (unsigned int i = 0; i < bondOrder.size(); ++i) {
    // RDKit❗🔝:     auto idx = bondOrder[i];
    // RDKit❗🔝:     const auto bond = mol.getBondWithIdx(idx);
    // RDKit❗🔝:     if (bond->getBondType() != bondType) {
    // RDKit❗🔝:       continue;
    // RDKit❗🔝:     }
    // RDKit❗🔝:     auto begAtomOrder =
    // RDKit❗🔝:         std::find(atomOrder.begin(), atomOrder.end(), bond->getBeginAtomIdx()) -
    // RDKit❗🔝:         atomOrder.begin();
    // RDKit❗🔝:     if (!res.empty()) {
    // RDKit❗🔝:       res += ",";
    // RDKit❗🔝:     } else {
    // RDKit❗🔝:       res = symbol + ":";
    // RDKit❗🔝:     }
    // RDKit❗🔝:     res += boost::str(boost::format("%d.%d") % begAtomOrder % i);
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return res;
    // RDKit❗🔝: }
    let mut result = PropertyText::new();
    for (i, &id) in bond_order.iter().enumerate() {
        let bond = bonds
            .get(id.index())
            .ok_or(SmilesParseError::CxTypedBondOutOfRange {
                bond: id,
                bond_count: bonds.len(),
            })?;
        if bond.order() != order {
            continue;
        }
        let begin = atom_order
            .iter()
            .position(|&atom| atom == bond.begin())
            .unwrap_or(atom_order.len());
        if !result.is_empty() {
            result.push_byte(b',');
        } else {
            result.extend_bytes(symbol);
            result.push_byte(b':');
        }
        write!(&mut result, "{begin}.{i}").expect("counted byte formatting is infallible");
    }
    // Cost 🔝: one growing output buffer, identical literal first std::find
    // and source loop order. Direct decimal writing removes each Native
    // boost::format construction/parse, dynamic format-item/argument storage
    // and materialized item string; borrowed closed bond rows add no graph
    // copy, atom conversion, query loss, map or per-item output buffer.
    Ok(result)
}

fn write_zero_bonds<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    bond_order: &[BondId],
) -> Result<PropertyText, SmilesParseError> {
    write_cx_zero_bonds_source(
        cosmolkit_core::stereo_graph::BondRows::Concrete(&record.into().topology.bonds),
        bond_order,
    )
}

/// Source zero-order CX emission over actual ordinary or query bond rows.
#[doc(hidden)]
pub fn write_cx_zero_bonds_source(
    bonds: cosmolkit_core::stereo_graph::BondRows<'_>,
    bond_order: &[BondId],
) -> Result<PropertyText, SmilesParseError> {
    // RDKit❗🔝: std::string get_zerobonds_block(const ROMol &mol,
    // RDKit❗🔝:                                 const std::vector<unsigned int> &,
    // RDKit❗🔝:                                 const std::vector<unsigned int> &bondOrder) {
    // RDKit❗🔝:   std::string res = "";
    // RDKit❗🔝:   for (unsigned int i = 0; i < bondOrder.size(); ++i) {
    // RDKit❗🔝:     auto idx = bondOrder[i];
    // RDKit❗🔝:     const auto bond = mol.getBondWithIdx(idx);
    // RDKit❗🔝:     if (bond->getBondType() != Bond::BondType::ZERO) {
    // RDKit❗🔝:       continue;
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (!res.empty()) {
    // RDKit❗🔝:       res += ",";
    // RDKit❗🔝:     } else {
    // RDKit❗🔝:       res = "Z:";
    // RDKit❗🔝:     }
    // RDKit❗🔝:     res += boost::str(boost::format("%d") % i);
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return res;
    // RDKit❗🔝: }
    let mut result = PropertyText::new();
    for (i, &id) in bond_order.iter().enumerate() {
        let bond = bonds
            .get(id.index())
            .ok_or(SmilesParseError::CxZeroBondOutOfRange {
                bond: id,
                bond_count: bonds.len(),
            })?;
        if bond.order() != BondOrder::Zero {
            continue;
        }
        if !result.is_empty() {
            result.push_byte(b',');
        } else {
            result.extend_bytes(b"Z:");
        }
        write!(&mut result, "{i}").expect("counted byte formatting is infallible");
    }
    // Cost 🔝: identical source-order loop, O(B) work and output storage;
    // direct decimal emission eliminates per-item Boost format construction,
    // parsed format/argument buffers and item strings, plus previous Rust
    // Vec<String>/join. Borrowed actual bond rows preserve query predicates.
    Ok(result)
}

fn other_atom(bond: &Bond, atom: AtomId) -> Option<AtomId> {
    if bond.begin() == atom {
        Some(bond.end())
    } else if bond.end() == atom {
        Some(bond.begin())
    } else {
        None
    }
}

fn write_ring_bond_stereo<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    rings: &RingInfo,
) -> String {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION get_ringbond_cistrans_block
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto rinfo = mol.getRingInfo();
    // RDKit✔️✔️:   std::string c = "", t = "", ctu = "";
    // RDKit✔️✔️:   for (unsigned int i = 0; i < bondOrder.size(); ++i) {
    // RDKit✔️✔️:     auto idx = bondOrder[i];
    // RDKit✔️✔️: if (!rinfo->numBondRings(idx) ||
    // RDKit✔️✔️:     rinfo->minBondRingSize(idx) < Chirality::minRingSizeForDoubleBondStereo) continue;
    // RDKit✔️✔️: if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit✔️✔️:     bond->getBondType() != Bond::BondType::AROMATIC) {
    // RDKit✔️✔️:   continue;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (bstereo == Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:   if (ctu.empty()) { ctu += "ctu:"; } else { ctu += ","; }
    // RDKit✔️✔️:   ctu += label;
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   Atom *begAtom = bond->getBeginAtom();
    // RDKit✔️✔️:   Atom *endAtom = bond->getEndAtom();
    // RDKit✔️✔️:   bool needSwap = false;
    // RDKit✔️✔️:   if (begAtom->getDegree() > 2) {
    // RDKit✔️✔️:     unsigned int o1 = atomOrder[bond->getStereoAtoms()[0]];
    // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit✔️✔️:       if (nbr == endAtom ||
    // RDKit✔️✔️:           nbr->getIdx() ==
    // RDKit✔️✔️:               static_cast<unsigned>(bond->getStereoAtoms()[0])) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (atomOrder[nbr->getIdx() < o1]) {
    // RDKit✔️✔️:         needSwap = !needSwap;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (endAtom->getDegree() > 2) {
    // RDKit✔️✔️:     unsigned int o1 = atomOrder[bond->getStereoAtoms()[1]];
    // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit✔️✔️:       if (nbr == begAtom ||
    // RDKit✔️✔️:           nbr->getIdx() ==
    // RDKit✔️✔️:               static_cast<unsigned>(bond->getStereoAtoms()[1])) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (atomOrder[nbr->getIdx() < o1]) {
    // RDKit✔️✔️:         needSwap = !needSwap;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bstereo == Bond::BondStereo::STEREOCIS || needSwap) {
    // RDKit✔️✔️:     if (c.empty()) { c += "c:"; } else { c += ","; }
    // RDKit✔️✔️:     c += label;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if (t.empty()) { t += "t:"; } else { t += ","; }
    // RDKit✔️✔️:     t += label;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return c + t + ctu;
    // END RDKIT CPP FUNCTION get_ringbond_cistrans_block
    // Behavior review: the controlling-neighbor condition intentionally
    // preserves the pinned C++ boolean subscript `atomOrder[nbrIdx < o1]`;
    // it is not an output-position comparison. Category blocks are likewise
    // concatenated without an added separator, exactly as `c + t + ctu`.
    // Complexity review: one bond-order pass plus adjacency scans matches the
    // source O(B + incident-neighbor visits) shape; all order lookups are O(1)
    // and the three category buffers are allocated once.
    const MIN_RING_SIZE: usize = 8;
    let mut cis = Vec::new();
    let mut trans = Vec::new();
    let mut unknown = Vec::new();
    for (bond_position, bond_id) in bond_order.iter().copied().enumerate() {
        if rings.num_bond_rings(bond_id) == 0 || rings.min_bond_ring_size(bond_id) < MIN_RING_SIZE {
            continue;
        }
        let bond = &record.topology.bonds[bond_id.index()];
        if !matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic)
            || !matches!(
                bond.stereo(),
                BondStereo::Any | BondStereo::Cis | BondStereo::Trans
            )
        {
            continue;
        }
        if bond.stereo() == BondStereo::Any {
            unknown.push(bond_position.to_string());
            continue;
        }
        let Some([begin_ref, end_ref]) = bond.stereo_atoms() else {
            continue;
        };
        let mut swap = false;
        for (center, opposite, reference) in [
            (bond.begin(), bond.end(), begin_ref),
            (bond.end(), bond.begin(), end_ref),
        ] {
            let neighbors = record.topology.adjacency.neighbors_of(center.index());
            if neighbors.len() <= 2 {
                continue;
            }
            let reference_order_value = atom_order[reference.index()].index();
            for neighbor in neighbors {
                let neighbor_atom =
                    other_atom(&record.topology.bonds[neighbor.bond.index()], center)
                        .expect("adjacency endpoint");
                if neighbor_atom != opposite
                    && neighbor_atom != reference
                    && atom_order[usize::from(neighbor_atom.index() < reference_order_value)]
                        .index()
                        != 0
                {
                    swap = !swap;
                }
            }
        }
        if bond.stereo() == BondStereo::Cis || swap {
            cis.push(bond_position.to_string());
        } else {
            trans.push(bond_position.to_string());
        }
    }
    let mut result = String::new();
    if !cis.is_empty() {
        result.push_str(&format!("c:{}", cis.join(",")));
    }
    if !trans.is_empty() {
        result.push_str(&format!("t:{}", trans.join(",")));
    }
    if !unknown.is_empty() {
        result.push_str(&format!("ctu:{}", unknown.join(",")));
    }
    result
}

fn stereo_kind_order(kind: StereoGroupKind) -> u8 {
    match kind {
        StereoGroupKind::Absolute => 0,
        StereoGroupKind::Or => 1,
        StereoGroupKind::And => 2,
    }
}

#[doc(hidden)]
pub fn assign_stereo_group_ids(groups: &mut [(StereoGroup, Vec<usize>)]) {
    // BEGIN RDKIT CPP FUNCTION StereoGroup.cpp::storeIdsInUse
    // RDKit✔️✔️: void storeIdsInUse(boost::dynamic_bitset<> &ids, StereoGroup &sg) {
    // RDKit✔️✔️:   const auto groupId = sg.getWriteId();
    // RDKit✔️✔️:   if (groupId == 0) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   } else if (groupId >= ids.size()) {
    // RDKit✔️✔️:     ids.resize(groupId + 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (ids[groupId]) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "StereoGroup ID " << groupId
    // RDKit✔️✔️:         << " is used by more than one group, and will be reassined"
    // RDKit✔️✔️:         << std::endl;
    // RDKit✔️✔️:     sg.setWriteId(0);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     ids[groupId] = true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION StereoGroup.cpp::storeIdsInUse
    // BEGIN RDKIT CPP FUNCTION StereoGroup.cpp::assignMissingIds
    // RDKit✔️✔️: void assignMissingIds(const boost::dynamic_bitset<> &ids, unsigned &nextId,
    // RDKit✔️✔️:                       StereoGroup &sg) {
    // RDKit✔️✔️:   if (sg.getWriteId() == 0) {
    // RDKit✔️✔️:     ++nextId;
    // RDKit✔️✔️:     while (nextId < ids.size() && ids[nextId]) {
    // RDKit✔️✔️:       ++nextId;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     sg.setWriteId(nextId);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION StereoGroup.cpp::assignMissingIds
    // BEGIN RDKIT CPP FUNCTION StereoGroup.cpp::assignStereoGroupIds
    // RDKit✔️✔️: void assignStereoGroupIds(std::vector<StereoGroup> &groups) {
    // RDKit✔️✔️:   if (groups.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   boost::dynamic_bitset<> andIds;
    // RDKit✔️✔️:   boost::dynamic_bitset<> orIds;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto &sg : groups) {
    // RDKit✔️✔️:     if (sg.getGroupType() == StereoGroupType::STEREO_AND) {
    // RDKit✔️✔️:       storeIdsInUse(andIds, sg);
    // RDKit✔️✔️:     } else if (sg.getGroupType() == StereoGroupType::STEREO_OR) {
    // RDKit✔️✔️:       storeIdsInUse(orIds, sg);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned andId = 0;
    // RDKit✔️✔️:   unsigned orId = 0;
    // RDKit✔️✔️:   for (auto &sg : groups) {
    // RDKit✔️✔️:     if (sg.getGroupType() == StereoGroupType::STEREO_AND) {
    // RDKit✔️✔️:       assignMissingIds(andIds, andId, sg);
    // RDKit✔️✔️:     } else if (sg.getGroupType() == StereoGroupType::STEREO_OR) {
    // RDKit✔️✔️:       assignMissingIds(orIds, orId, sg);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION StereoGroup.cpp::assignStereoGroupIds
    // Behavior review: use independent packed ID bitmaps, retain the first
    // unique nonzero write ID in each type, zero later duplicates, then fill
    // missing IDs in a separate source-order pass. Read IDs are never read.
    // The duplicate-ID warning is emitted in the same source branch.
    // Complexity review: two O(n) group passes, direct bit lookups, and
    // dynamic bit storage proportional to the largest explicit ID, matching
    // the source's dynamic_bitset rather than a sparse-map substitute.
    let mut and_ids = Vec::<bool>::new();
    let mut or_ids = Vec::<bool>::new();

    for (group, _) in groups.iter_mut() {
        let used_ids = match group.kind() {
            StereoGroupKind::And => &mut and_ids,
            StereoGroupKind::Or => &mut or_ids,
            StereoGroupKind::Absolute => continue,
        };
        let write_id = stereo_group_write_id(group);
        if write_id == 0 {
            continue;
        }
        let index = write_id as usize;
        if index >= used_ids.len() {
            used_ids.resize(index + 1, false);
        }
        if used_ids[index] {
            eprintln!(
                "StereoGroup ID {write_id} is used by more than one group, and will be reassined"
            );
            set_stereo_group_write_id(group, 0);
        } else {
            used_ids[index] = true;
        }
    }

    let mut next_and = 0_u32;
    let mut next_or = 0_u32;
    for (group, _) in groups.iter_mut() {
        if group.kind() == StereoGroupKind::Absolute || stereo_group_write_id(group) != 0 {
            continue;
        }
        let (used_ids, next_id) = match group.kind() {
            StereoGroupKind::And => (&and_ids, &mut next_and),
            StereoGroupKind::Or => (&or_ids, &mut next_or),
            StereoGroupKind::Absolute => continue,
        };
        *next_id = next_id.wrapping_add(1);
        while (*next_id as usize) < used_ids.len() && used_ids[*next_id as usize] {
            *next_id = next_id.wrapping_add(1);
        }
        set_stereo_group_write_id(group, *next_id);
    }
}

fn write_enhanced_stereo<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
    wedge_bonds: &WedgeAssignments,
) -> Result<String, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION getSortedMappedIndexes
    // RDKit✔️❌: std::vector<unsigned> getSortedMappedIndexes(
    // RDKit✔️❌:     const std::vector<unsigned int> &atomIds,
    // RDKit✔️❌:     const std::vector<unsigned> &revOrder) {
    // RDKit✔️❌:   std::vector<unsigned> res;
    // RDKit✔️❌:   res.reserve(atomIds.size());
    // RDKit✔️❌:   for (auto atomId : atomIds) {
    // RDKit✔️❌:     res.push_back(revOrder[atomId]);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   std::sort(res.begin(), res.end());
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION getSortedMappedIndexes
    // BEGIN RDKIT CPP FUNCTION getSortedStereoGroupsAndIndices
    // RDKit✔️❌:   auto &groups = mol.getStereoGroups();
    // RDKit✔️❌:   std::vector<StGrpIdxPair> sortingGroups;
    // RDKit✔️❌:   sortingGroups.reserve(groups.size());
    // RDKit✔️❌:   for (const auto &sg : groups) {
    // RDKit✔️❌:     std::vector<unsigned int> atomIds;
    // RDKit✔️❌:     Atropisomers::getAllAtomIdsForStereoGroup(mol, sg, atomIds, wedgeBonds);
    // RDKit✔️❌:     const auto newAtomIndexes = getSortedMappedIndexes(atomIds, revOrder);
    // RDKit✔️❌:     if (!newAtomIndexes.empty()) {
    // RDKit✔️❌:       sortingGroups.emplace_back(sg, newAtomIndexes);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   std::sort(sortingGroups.begin(), sortingGroups.end(),
    // RDKit✔️❌:             [](const StGrpIdxPair &a, const StGrpIdxPair &b) {
    // RDKit✔️❌:               const auto &[sgA, idxsA] = a;
    // RDKit✔️❌:               const auto &[sgB, idxsB] = b;
    // RDKit✔️❌:               if (sgA.getGroupType() == sgB.getGroupType()) {
    // RDKit✔️❌:                 if (sgA.getWriteId() == sgB.getWriteId()) {
    // RDKit✔️❌:                   return idxsA < idxsB;
    // RDKit✔️❌:                 }
    // RDKit✔️❌:                 return sgA.getWriteId() < sgB.getWriteId();
    // RDKit✔️❌:               }
    // RDKit✔️❌:               return sgA.getGroupType() < sgB.getGroupType();
    // RDKit✔️❌:             });
    // END RDKIT CPP FUNCTION getSortedStereoGroupsAndIndices
    // BEGIN RDKIT CPP FUNCTION get_enhanced_stereo_block
    // RDKit✔️❌:   if (mol.getStereoGroups().empty()) {
    // RDKit✔️❌:     return "";
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (unsigned i = 0; i < atomOrder.size(); ++i) {
    // RDKit✔️❌:     revOrder[atomOrder[i]] = i;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   auto [groups, groupsAtoms] =
    // RDKit✔️❌:       getSortedStereoGroupsAndIndices(mol, revOrder, wedgeBonds);
    // RDKit✔️❌:   assignStereoGroupIds(groups);
    // RDKit✔️❌:   switch (sgItr->getGroupType()) {
    // RDKit✔️❌:     case StereoGroupType::STEREO_ABSOLUTE:
    // RDKit✔️❌:       res << "a:";
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case StereoGroupType::STEREO_OR:
    // RDKit✔️❌:       res << "o" << sgItr->getWriteId() << ":";
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     case StereoGroupType::STEREO_AND:
    // RDKit✔️❌:       res << "&" << sgItr->getWriteId() << ":";
    // RDKit✔️❌:       break;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (const auto &aid : *grpAtomsItr) {
    // RDKit✔️❌:     res << aid << ",";
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION get_enhanced_stereo_block
    // Behavior review: the checked batch returns source-ordered rows, with
    // direct members followed by qualifying group-bond endpoints. This loop
    // zips those rows to source group order, maps and sorts each row, omits
    // empty mapped groups, then keeps the source type/write-ID/index sort,
    // assignment and text stages unchanged.
    // Complexity review: unlike the pinned pointer-valid source, the detached
    // boundary validates O(V+E) topology once for the whole batch. It removes
    // the former O(G*(V+E)) repeated validation. The returned batch retains
    // O(total collected IDs) raw rows until mapping, in addition to the source
    // output rows; this known peak-memory cost keeps the complexity marker ❌.
    if record.topology.stereo_groups.is_empty() {
        return Ok(String::new());
    }
    let positions = atom_positions(atom_order, record.topology.atoms.len());
    let mut groups = Vec::with_capacity(record.topology.stereo_groups.len());
    let atom_ids_by_group = get_all_atom_ids_for_stereo_groups(
        &record.topology,
        &record.topology.stereo_groups,
        wedge_bonds,
    )
    .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
    for (group, atom_ids) in record.topology.stereo_groups.iter().zip(atom_ids_by_group) {
        let mut atoms = Vec::with_capacity(atom_ids.len());
        for atom in atom_ids {
            if let Some(position) = positions[atom.index()] {
                atoms.push(position);
            }
        }
        atoms.sort_unstable();
        if !atoms.is_empty() {
            groups.push((group.clone(), atoms));
        }
    }
    groups.sort_by(|(left_group, left_atoms), (right_group, right_atoms)| {
        stereo_kind_order(left_group.kind())
            .cmp(&stereo_kind_order(right_group.kind()))
            .then_with(|| {
                stereo_group_write_id(left_group).cmp(&stereo_group_write_id(right_group))
            })
            .then_with(|| left_atoms.cmp(right_atoms))
    });
    assign_stereo_group_ids(&mut groups);
    Ok(groups
        .into_iter()
        .map(|(group, atoms)| {
            let prefix = match group.kind() {
                StereoGroupKind::Absolute => "a".to_owned(),
                StereoGroupKind::Or => format!("o{}", stereo_group_write_id(&group)),
                StereoGroupKind::And => format!("&{}", stereo_group_write_id(&group)),
            };
            format!(
                "{prefix}:{}",
                atoms
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect::<Vec<_>>()
        .join(","))
}

#[derive(Debug, PartialEq, Eq)]
struct SourceLinkNode {
    min_repetitions: u32,
    max_repetitions: u32,
    bond_atoms: [(u32, u32); 2],
}

#[derive(Debug, PartialEq, Eq)]
enum SourceLinkIssueKind {
    UnsignedConversion,
    TooFewValues,
    BadCounts,
    UnderConstruction,
    DifferentCenters,
    MissingMappedAtom(u32),
    MissingBond,
}

#[derive(Debug, PartialEq, Eq)]
struct SourceLinkIssue {
    kind: SourceLinkIssueKind,
    raw_node: PropertyText,
}

fn source_link_issue(
    strict: bool,
    warnings: &mut Vec<SourceLinkIssue>,
    kind: SourceLinkIssueKind,
    raw_node: &[u8],
) -> Result<(), SmilesParseError> {
    // RDKit❗✔️: #define UNDER_CONSTRUCTION(fn)                                           \
    // RDKit❗✔️:   Invar::Invariant inv("Incomplete Code",                                \
    // RDKit❗✔️:                        "This routine is still under development", fn,    \
    // RDKit❗✔️:                        __FILE__, __LINE__);                              \
    // RDKit❗✔️:   BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n" << std::endl; \
    // RDKit❗✔️:   throw inv;
    // Preserve raw diagnostic bytes in the detached warning. The existing CX
    // error display translates native ValueError/Incomplete Code structurally;
    // this source-defined exception is never labeled an unsupported input.
    let message = source_link_diagnostic(&kind, raw_node, strict);
    if strict {
        return Err(SmilesParseError::Cx(
            String::from_utf8_lossy(message.as_bytes()).into_owned(),
        ));
    }
    warnings.push(SourceLinkIssue {
        kind,
        raw_node: PropertyText::from_bytes(raw_node),
    });
    Ok(())
}

fn source_link_diagnostic(kind: &SourceLinkIssueKind, raw: &[u8], strict: bool) -> PropertyText {
    // Each literal is the corresponding pinned LinkNode.h stream insertion;
    // raw node bytes are never interpreted as UTF-8 for a warning.
    let (prefix, suffix) = match kind {
        SourceLinkIssueKind::UnsignedConversion => {
            ("Cannot convert values in LINKNODE '", "' to unsigned ints")
        }
        SourceLinkIssueKind::TooFewValues => ("not enough values in LINKNODE '", "'"),
        SourceLinkIssueKind::BadCounts => ("bad counts in LINKNODE '", "'"),
        SourceLinkIssueKind::DifferentCenters => {
            ("bonds don't start at the same atom for LINKNODE '", "'")
        }
        SourceLinkIssueKind::MissingBond => ("bond not found between atoms in LINKNODE '", "'"),
        SourceLinkIssueKind::UnderConstruction => {
            return PropertyText::from(if strict {
                "Incomplete Code: only link nodes with 2 bonds are currently supported"
            } else {
                "only link nodes with 2 bonds are currently supported"
            });
        }
        SourceLinkIssueKind::MissingMappedAtom(index) => {
            let mut out = PropertyText::new();
            write!(
                &mut out,
                "atom index {index} cannot be found in molecule for LINKNODE '"
            )
            .expect("byte formatting");
            out.extend_bytes(raw);
            out.push_byte(b'\'');
            return out;
        }
    };
    let mut out = PropertyText::from(prefix);
    out.extend_bytes(raw);
    out.extend_bytes(suffix.as_bytes());
    out
}

/// Deliver one counted-byte Native warning, preserving the stream newline.
#[doc(hidden)]
pub fn emit_cx_link_node_warning_source(message: &[u8]) {
    use std::io::Write as _;
    let mut sink = std::io::stderr().lock();
    // Native log stream failure does not change chemistry control flow.
    let _ = sink.write_all(message);
    let _ = sink.write_all(b"\n");
}

fn source_link_tokens(input: &[u8], delimiter: u8) -> impl Iterator<Item = &[u8]> {
    // Boost❗🔝:     explicit
    // Boost❗🔝:     char_separator(const Char* dropped_delims,
    // Boost❗🔝:                    const Char* kept_delims = 0,
    // Boost❗🔝:                    empty_token_policy empty_tokens = drop_empty_tokens)
    // Boost❗🔝:       : m_dropped_delims(dropped_delims),
    // Boost❗🔝:         m_use_ispunct(false),
    // Boost❗🔝:         m_use_isspace(false),
    // Boost❗🔝:         m_empty_tokens(empty_tokens),
    // Boost❗🔝:         m_output_done(false)
    // Boost❗🔝:     {
    // Boost❗🔝:       // Borland workaround
    // Boost❗🔝:       if (kept_delims)
    // Boost❗🔝:         m_kept_delims = kept_delims;
    // Boost❗🔝:     }
    // Boost❗🔝:
    // Boost❗🔝:     template <typename InputIterator, typename Token>
    // Boost❗🔝:     bool operator()(InputIterator& next, InputIterator end, Token& tok)
    // Boost❗🔝:     {
    // Boost❗🔝:       typedef tokenizer_detail::assign_or_plus_equal<
    // Boost❗🔝:         BOOST_DEDUCED_TYPENAME tokenizer_detail::get_iterator_category<
    // Boost❗🔝:           InputIterator
    // Boost❗🔝:         >::iterator_category
    // Boost❗🔝:       > assigner;
    // Boost❗🔝:
    // Boost❗🔝:       assigner::clear(tok);
    // Boost❗🔝:
    // Boost❗🔝:       // skip past all dropped_delims
    // Boost❗🔝:       if (m_empty_tokens == drop_empty_tokens)
    // Boost❗🔝:         for (; next != end  && is_dropped(*next); ++next)
    // Boost❗🔝:           { }
    // Boost❗🔝:
    // Boost❗🔝:       InputIterator start(next);
    // Boost❗🔝:
    // Boost❗🔝:       if (m_empty_tokens == drop_empty_tokens) {
    // Boost❗🔝:
    // Boost❗🔝:         if (next == end)
    // Boost❗🔝:           return false;
    // Boost❗🔝:
    // Boost❗🔝:
    // Boost❗🔝:         // if we are on a kept_delims move past it and stop
    // Boost❗🔝:         if (is_kept(*next)) {
    // Boost❗🔝:           assigner::plus_equal(tok,*next);
    // Boost❗🔝:           ++next;
    // Boost❗🔝:         } else
    // Boost❗🔝:           // append all the non delim characters
    // Boost❗🔝:           for (; next != end && !is_dropped(*next) && !is_kept(*next); ++next)
    // Boost❗🔝:             assigner::plus_equal(tok,*next);
    // Boost❗🔝:       }
    // Boost❗🔝:       else { // m_empty_tokens == keep_empty_tokens
    // Boost❗🔝:
    // Boost❗🔝:         // Handle empty token at the end
    // Boost❗🔝:         if (next == end)
    // Boost❗🔝:         {
    // Boost❗🔝:           if (m_output_done == false)
    // Boost❗🔝:           {
    // Boost❗🔝:             m_output_done = true;
    // Boost❗🔝:             assigner::assign(start,next,tok);
    // Boost❗🔝:             return true;
    // Boost❗🔝:           }
    // Boost❗🔝:           else
    // Boost❗🔝:             return false;
    // Boost❗🔝:         }
    // Boost❗🔝:
    // Boost❗🔝:         if (is_kept(*next)) {
    // Boost❗🔝:           if (m_output_done == false)
    // Boost❗🔝:             m_output_done = true;
    // Boost❗🔝:           else {
    // Boost❗🔝:             assigner::plus_equal(tok,*next);
    // Boost❗🔝:             ++next;
    // Boost❗🔝:             m_output_done = false;
    // Boost❗🔝:           }
    // Boost❗🔝:         }
    // Boost❗🔝:         else if (m_output_done == false && is_dropped(*next)) {
    // Boost❗🔝:           m_output_done = true;
    // Boost❗🔝:         }
    // Boost❗🔝:         else {
    // Boost❗🔝:           if (is_dropped(*next))
    // Boost❗🔝:             start=++next;
    // Boost❗🔝:           for (; next != end && !is_dropped(*next) && !is_kept(*next); ++next)
    // Boost❗🔝:             assigner::plus_equal(tok,*next);
    // Boost❗🔝:           m_output_done = true;
    // Boost❗🔝:         }
    // Boost❗🔝:       }
    // Boost❗🔝:       assigner::assign(start,next,tok);
    // Boost❗🔝:       return true;
    // Boost❗🔝:     }
    // Boost❗🔝:     bool is_kept(Char E) const
    // Boost❗🔝:     {
    // Boost❗🔝:       if (m_kept_delims.length())
    // Boost❗🔝:         return m_kept_delims.find(E) != string_type::npos;
    // Boost❗🔝:       else if (m_use_ispunct) {
    // Boost❗🔝:         return Traits::ispunct(E) != 0;
    // Boost❗🔝:       } else
    // Boost❗🔝:         return false;
    // Boost❗🔝:     }
    // Boost❗🔝:     bool is_dropped(Char E) const
    // Boost❗🔝:     {
    // Boost❗🔝:       if (m_dropped_delims.length())
    // Boost❗🔝:         return m_dropped_delims.find(E) != string_type::npos;
    // Boost❗🔝:       else if (m_use_isspace) {
    // Boost❗🔝:         return Traits::isspace(E) != 0;
    // Boost❗🔝:       } else
    // Boost❗🔝:         return false;
    // Boost❗🔝:     }
    // Explicit char_separator("|") / char_separator(" ") only: no kept
    // delimiters, no locale whitespace/punctuation and drop_empty_tokens.
    // Borrowing slices replaces native per-token copies without changing bytes.
    input
        .split(move |byte| *byte == delimiter)
        .filter(|token| !token.is_empty())
}

fn source_link_unsigned(token: &[u8]) -> Option<u32> {
    // lexical_cast does not trim whitespace. The sole CORE from_rdvalue uint
    // owner has the full Boost unsigned converter anchor; prevent its explicit
    // RDValue trim_right adaptation before delegating the same decimal/sign loop.
    if token
        .iter()
        .any(|byte| matches!(*byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c))
    {
        return None;
    }
    cosmolkit_core::property_value_to_uint(&PropertyValue::String(PropertyText::from_bytes(token)))
        .ok()
}

fn source_link_has_bond<G: StereoGraphAccess>(
    topology: &G,
    first: u32,
    second: u32,
) -> Result<bool, SmilesParseError> {
    // RDKit❗✔️: const Bond *ROMol::getBondBetweenAtoms(unsigned int idx1,
    // RDKit❗✔️:                                        unsigned int idx2) const {
    // RDKit❗✔️:   URANGE_CHECK(idx1, getNumAtoms());
    // RDKit❗✔️:   URANGE_CHECK(idx2, getNumAtoms());
    // RDKit❗✔️:   const Bond *res = nullptr;
    // RDKit❗✔️:
    // RDKit❗✔️:   auto [edge, found] = boost::edge(boost::vertex(idx1, d_graph),
    // RDKit❗✔️:                                    boost::vertex(idx2, d_graph), d_graph);
    // RDKit❗✔️:   if (found) {
    // RDKit❗✔️:     res = d_graph[edge];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // Native range checks precede edge lookup and remain errors in non-strict
    // mode. Lookup follows the first endpoint's adjacency, including self edges.
    let atom_count = topology.atoms().len();
    for atom in [first, second] {
        if atom as usize >= atom_count {
            return Err(SmilesParseError::Cx(format!(
                "link-node atom index {atom} is out of range for {atom_count} atoms"
            )));
        }
    }
    Ok(topology
        .adjacency()
        .neighbors_of(first as usize)
        .iter()
        .any(|neighbor| neighbor.atom_index == second as usize))
}

#[cfg(test)]
fn get_mol_link_nodes_source(
    topology: &cosmolkit_model::TopologyBlock,
    properties: &cosmolkit_model::MoleculeProperties,
    strict: bool,
    atom_index_map: Option<&BTreeMap<u32, AtomId>>,
    warnings: &mut Vec<SourceLinkIssue>,
) -> Result<Vec<SourceLinkNode>, SmilesParseError> {
    get_mol_link_nodes_impl(
        topology,
        properties.prop("_molLinkNodes"),
        strict,
        atom_index_map,
        warnings,
        &mut |_| {},
    )
}

fn get_mol_link_nodes_impl<G: StereoGraphAccess>(
    topology: &G,
    property: Option<&PropertyValue>,
    strict: bool,
    atom_index_map: Option<&BTreeMap<u32, AtomId>>,
    warnings: &mut Vec<SourceLinkIssue>,
    warning_sink: &mut dyn FnMut(&[u8]),
) -> Result<Vec<SourceLinkNode>, SmilesParseError> {
    // RDKit❗✔️: inline std::vector<LinkNode> getMolLinkNodes(
    // RDKit❗✔️:     const ROMol &mol, bool strict = true,
    // RDKit❗✔️:     const std::map<unsigned, Atom *> *atomIdxMap = nullptr) {
    // RDKit❗✔️:   std::vector<LinkNode> res;
    // RDKit❗✔️:   std::string pval;
    // RDKit❗✔️:   if (!mol.getPropIfPresent(common_properties::molFileLinkNodes, pval)) {
    // RDKit❗✔️:     return res;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::vector<int> mapping;
    // RDKit❗✔️:
    // RDKit❗✔️:   boost::char_separator<char> pipesep("|");
    // RDKit❗✔️:   boost::char_separator<char> spacesep(" ");
    // RDKit❗✔️:   for (auto linknodetext : tokenizer(pval, pipesep)) {
    // RDKit❗✔️:     LinkNode node;
    // RDKit❗✔️:     tokenizer tokens(linknodetext, spacesep);
    // RDKit❗✔️:     std::vector<unsigned int> data;
    // RDKit❗✔️:     try {
    // RDKit❗✔️:       std::transform(tokens.begin(), tokens.end(), std::back_inserter(data),
    // RDKit❗✔️:                      [](const std::string &token) -> unsigned int {
    // RDKit❗✔️:                        return boost::lexical_cast<unsigned int>(token);
    // RDKit❗✔️:                      });
    // RDKit❗✔️:     } catch (boost::bad_lexical_cast &) {
    // RDKit❗✔️:       std::ostringstream errout;
    // RDKit❗✔️:       errout << "Cannot convert values in LINKNODE '" << linknodetext
    // RDKit❗✔️:              << "' to unsigned ints";
    // RDKit❗✔️:       if (strict) {
    // RDKit❗✔️:         throw ValueErrorException(errout.str());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // the second test here is for the atom-pairs defining the bonds
    // RDKit❗✔️:     // data[2] contains the number of bonds
    // RDKit❗✔️:     if (data.size() < 5 || data.size() < 3 + 2 * data[2]) {
    // RDKit❗✔️:       std::ostringstream errout;
    // RDKit❗✔️:       errout << "not enough values in LINKNODE '" << linknodetext << "'";
    // RDKit❗✔️:       if (strict) {
    // RDKit❗✔️:         throw ValueErrorException(errout.str());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     node.minRep = data[0];
    // RDKit❗✔️:     node.maxRep = data[1];
    // RDKit❗✔️:     if (node.minRep == 0 || node.maxRep < node.minRep) {
    // RDKit❗✔️:       std::ostringstream errout;
    // RDKit❗✔️:       errout << "bad counts in LINKNODE '" << linknodetext << "'";
    // RDKit❗✔️:       if (strict) {
    // RDKit❗✔️:         throw ValueErrorException(errout.str());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     node.nBonds = data[2];
    // RDKit❗✔️:     if (node.nBonds != 2) {
    // RDKit❗✔️:       if (strict) {
    // RDKit❗✔️:         UNDER_CONSTRUCTION(
    // RDKit❗✔️:             "only link nodes with 2 bonds are currently supported");
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "only link nodes with 2 bonds are currently supported"
    // RDKit❗✔️:             << std::endl;
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // both bonds must start from the same atom:
    // RDKit❗✔️:     if (data[3] != data[5]) {
    // RDKit❗✔️:       std::ostringstream errout;
    // RDKit❗✔️:       errout << "bonds don't start at the same atom for LINKNODE '"
    // RDKit❗✔️:              << linknodetext << "'";
    // RDKit❗✔️:       if (strict) {
    // RDKit❗✔️:         throw ValueErrorException(errout.str());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     if (atomIdxMap) {
    // RDKit❗✔️:       // map the indices back to the original atom numbers
    // RDKit❗✔️:       for (unsigned int i = 3; i <= 6; ++i) {
    // RDKit❗✔️:         const auto aidx = atomIdxMap->find(data[i] - 1);
    // RDKit❗✔️:         if (aidx == atomIdxMap->end()) {
    // RDKit❗✔️:           std::ostringstream errout;
    // RDKit❗✔️:           errout << "atom index " << data[i]
    // RDKit❗✔️:                  << " cannot be found in molecule for LINKNODE '"
    // RDKit❗✔️:                  << linknodetext << "'";
    // RDKit❗✔️:           if (strict) {
    // RDKit❗✔️:             throw ValueErrorException(errout.str());
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:             continue;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           data[i] = aidx->second->getIdx();
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       for (unsigned int i = 3; i <= 6; ++i) {
    // RDKit❗✔️:         --data[i];
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     node.bondAtoms.push_back(std::make_pair(data[3], data[4]));
    // RDKit❗✔️:     node.bondAtoms.push_back(std::make_pair(data[5], data[6]));
    // RDKit❗✔️:     if (!mol.getBondBetweenAtoms(data[4], data[3]) ||
    // RDKit❗✔️:         !mol.getBondBetweenAtoms(data[6], data[5])) {
    // RDKit❗✔️:       std::ostringstream errout;
    // RDKit❗✔️:       errout << "bond not found between atoms in LINKNODE '" << linknodetext
    // RDKit❗✔️:              << "'";
    // RDKit❗✔️:       if (strict) {
    // RDKit❗✔️:         throw ValueErrorException(errout.str());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     res.push_back(std::move(node));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // Behavior: all native strict/non-strict branches and both index-map paths.
    // Missing mapped indices in non-strict mode continue only the inner loop,
    // retaining that original one-based value; validation then uses real bonds.
    // Complexity: linear byte/token work, one u32 vector per node and one result
    // vector; four logarithmic map lookups and at most two degree-bound edge
    // searches. No topology/property clone or synthetic mapping/cached facts.
    let mut result = Vec::new();
    let Some(raw) = property else {
        return Ok(result);
    };
    let raw = property_value_to_string(raw).map_err(SmilesParseError::WriterProperty)?;
    let mut issue = |kind, item: &[u8]| -> Result<(), SmilesParseError> {
        source_link_issue(strict, warnings, kind, item)?;
        let warning = warnings.last().expect("non-strict issue records warning");
        warning_sink(
            source_link_diagnostic(&warning.kind, warning.raw_node.as_bytes(), false).as_bytes(),
        );
        Ok(())
    };
    for item in source_link_tokens(raw.as_bytes(), b'|') {
        let Some(mut data) = source_link_tokens(item, b' ')
            .map(source_link_unsigned)
            .collect::<Option<Vec<_>>>()
        else {
            issue(SourceLinkIssueKind::UnsignedConversion, item)?;
            continue;
        };
        if data.len() < 5 || data.len() < 3u32.wrapping_add(2u32.wrapping_mul(data[2])) as usize {
            issue(SourceLinkIssueKind::TooFewValues, item)?;
            continue;
        }
        let min_repetitions = data[0];
        let max_repetitions = data[1];
        if min_repetitions == 0 || max_repetitions < min_repetitions {
            issue(SourceLinkIssueKind::BadCounts, item)?;
            continue;
        }
        if data[2] != 2 {
            issue(SourceLinkIssueKind::UnderConstruction, item)?;
            continue;
        }
        if data[3] != data[5] {
            issue(SourceLinkIssueKind::DifferentCenters, item)?;
            continue;
        }
        if let Some(map) = atom_index_map {
            for value in &mut data[3..=6] {
                if let Some(atom) = map.get(&value.wrapping_sub(1)) {
                    *value = atom.index() as u32;
                } else {
                    issue(SourceLinkIssueKind::MissingMappedAtom(*value), item)?;
                    continue;
                }
            }
        } else {
            for value in &mut data[3..=6] {
                *value = value.wrapping_sub(1);
            }
        }
        let bond_atoms = [(data[3], data[4]), (data[5], data[6])];
        if !source_link_has_bond(topology, data[4], data[3])?
            || !source_link_has_bond(topology, data[6], data[5])?
        {
            issue(SourceLinkIssueKind::MissingBond, item)?;
            continue;
        }
        result.push(SourceLinkNode {
            min_repetitions,
            max_repetitions,
            bond_atoms,
        });
    }
    Ok(result)
}

fn write_link_nodes<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
) -> Result<String, SmilesParseError> {
    let record = record.into();
    write_cx_link_nodes_impl(
        record.topology,
        record.properties.prop("_molLinkNodes"),
        atom_order,
        &mut emit_cx_link_node_warning_source,
    )
}

/// Borrow actual query carriers and the actual counted source property.
#[doc(hidden)]
pub fn write_query_cx_link_nodes_source(
    query: &cosmolkit_model::QueryGraph,
    atom_order: &[AtomId],
    warning_sink: &mut dyn FnMut(&[u8]),
) -> Result<String, SmilesParseError> {
    write_cx_link_nodes_impl(query, query.prop("_molLinkNodes"), atom_order, warning_sink)
}

fn write_cx_link_nodes_impl<G: StereoGraphAccess>(
    graph: &G,
    property: Option<&PropertyValue>,
    atom_order: &[AtomId],
    warning_sink: &mut dyn FnMut(&[u8]),
) -> Result<String, SmilesParseError> {
    // RDKit❗❌: std::string get_linknodes_block(const ROMol &mol,
    // RDKit❗❌:                                 const std::vector<unsigned int> &atomOrder) {
    // RDKit❗❌:   bool strict = false;
    // RDKit❗❌:   auto linkNodes = MolEnumerator::utils::getMolLinkNodes(mol, strict);
    // RDKit❗❌:   if (linkNodes.empty()) {
    // RDKit❗❌:     return "";
    // RDKit❗❌:   }
    // RDKit❗❌:   // we need a map from original atom idx to output idx:
    // RDKit❗❌:   std::vector<unsigned int> revOrder(mol.getNumAtoms());
    // RDKit❗❌:   for (unsigned i = 0; i < atomOrder.size(); ++i) {
    // RDKit❗❌:     revOrder[atomOrder[i]] = i;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::stringstream res;
    // RDKit❗❌:   res << "LN:";
    // RDKit❗❌:   for (const auto &ln : linkNodes) {
    // RDKit❗❌:     unsigned int atomIdx = atomOrder[ln.bondAtoms[0].first];
    // RDKit❗❌:     res << atomIdx << ":" << ln.minRep << "." << ln.maxRep;
    // RDKit❗❌:     if (mol.getAtomWithIdx(ln.bondAtoms[0].first)->getDegree() > 2) {
    // RDKit❗❌:       // include the outer atom indices
    // RDKit❗❌:       res << "." << atomOrder[ln.bondAtoms[0].second] << "."
    // RDKit❗❌:           << atomOrder[ln.bondAtoms[1].second];
    // RDKit❗❌:     }
    // RDKit❗❌:     res << ",";
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::string resStr = res.str();
    // RDKit❗❌:   if (!resStr.empty() && resStr.back() == ',') {
    // RDKit❗❌:     resStr.pop_back();
    // RDKit❗❌:   }
    // RDKit❗❌:   return resStr;
    // RDKit❗❌: }
    let mut warnings = Vec::new();
    let nodes = get_mol_link_nodes_impl(graph, property, false, None, &mut warnings, warning_sink)?;
    if nodes.is_empty() {
        return Ok(String::new());
    }
    // Preserve the source zero-initialization and last-wins writes even
    // though the pinned source does not subsequently read this vector.
    let mut reverse_order = vec![0u32; graph.atoms().len()];
    for (i, &id) in atom_order.iter().enumerate() {
        let count = reverse_order.len();
        *reverse_order
            .get_mut(id.index())
            .ok_or(SmilesParseError::CxLinkAtomOutOfRange {
                atom: id.index(),
                atom_count: count,
            })? = i as u32;
    }
    let output_at = |index: usize| {
        atom_order.get(index).map(|a| a.index() as u32).ok_or(
            SmilesParseError::CxLinkOrderOutOfRange {
                index,
                count: atom_order.len(),
            },
        )
    };
    let mut result = String::from("LN:");
    for node in nodes {
        let center = node.bond_atoms[0].0 as usize;
        write!(
            &mut result,
            "{}:{}.{}",
            output_at(center)?,
            node.min_repetitions,
            node.max_repetitions
        )
        .expect("string formatting");
        if graph.adjacency().neighbors_of(center).len() > 2 {
            write!(
                &mut result,
                ".{}.{}",
                output_at(node.bond_atoms[0].1 as usize)?,
                output_at(node.bond_atoms[1].1 as usize)?
            )
            .expect("string formatting");
        }
        result.push(',');
    }
    if result.as_bytes().last() == Some(&b',') {
        result.pop();
    }
    // Cost ❌: the existing detached parser collects warning records as well
    // as emitting source diagnostics. Source tokens/result/revOrder loops,
    // edge checks and forward indexes retained. One output buffer replaces
    // previous per-node strings/Vec/join; no query/Atom/property-map clones.
    Ok(result)
}

fn is_data_sgroup(group: &SubstanceGroup) -> Result<bool, SmilesParseError> {
    // BEGIN COMPLETE RDProps::getPropIfPresent
    // RDKit✔️❌: bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit✔️❌:     return d_props.getValIfPresent(key, res);
    // RDKit✔️❌:   }
    // END COMPLETE RDProps::getPropIfPresent

    // A present source property is read through the sole CORE conversion;
    // wrong supported types/conversions propagate. The typed detached kind
    // supplies TYPE only when no explicit source property was supplied.
    if let Some(value) = group.props().get(b"TYPE".as_slice()) {
        return Ok(cosmolkit_core::property_value_to_string(value)
            .map_err(SmilesParseError::WriterProperty)?
            .as_bytes()
            == b"DAT");
    }
    Ok(matches!(group.kind(), SubstanceGroupKind::Data))
}

fn data_sgroup_value(group: &SubstanceGroup, key: &str) -> Result<PropertyText, SmilesParseError> {
    // BEGIN COMPLETE RDProps::getPropIfPresent
    // RDKit✔️❌: bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit✔️❌:     return d_props.getValIfPresent(key, res);
    // RDKit✔️❌:   }
    // END COMPLETE RDProps::getPropIfPresent

    // Scalar source getProp<string> uses canonical RDValue conversion. This
    // helper is never used to reinterpret DATAFIELDS as a scalar string.
    if let Some(value) = group.props().get(key.as_bytes()) {
        return cosmolkit_core::property_value_to_string(value)
            .map_err(SmilesParseError::WriterProperty);
    }
    Ok(group
        .data()
        .and_then(|data| match key {
            "FIELDNAME" => data.field_name.as_ref(),
            "QUERYOP" => data.query_op.as_ref(),
            "FIELDINFO" => data.field_info.as_ref(),
            _ => None,
        })
        .cloned()
        .unwrap_or_default())
}

fn data_sgroup_values(group: &SubstanceGroup) -> Result<&[PropertyText], SmilesParseError> {
    // BEGIN COMPLETE RDProps::getPropIfPresent
    // RDKit✔️❌: bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit✔️❌:     return d_props.getValIfPresent(key, res);
    // RDKit✔️❌:   }
    // END COMPLETE RDProps::getPropIfPresent

    // Source DATAFIELDS is vector<string>: exact tag, bytes and element order.
    // A present scalar is a native cast failure, never split or inferred.
    if let Some(value) = group.props().get(b"DATAFIELDS".as_slice()) {
        return value
            .as_string_vector()
            .map_err(SmilesParseError::WriterPropertyKind);
    }
    Ok(group
        .data()
        .map(|data| data.values.as_slice())
        .filter(|values| !values.is_empty())
        .unwrap_or_else(|| group.data_fields()))
}

fn write_data_sgroups<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
) -> Result<PropertyText, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION get_sgroup_data_block
    // RDKit✔️✔️: std::string get_sgroup_data_block(const ROMol &mol,
    // RDKit✔️✔️:                                   const std::vector<unsigned int> &atomOrder) {
    // RDKit✔️✔️:   const auto &sgs = getSubstanceGroups(mol);
    // RDKit✔️✔️:   if (sgs.empty()) {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int sgroupOutputIndex = 0;
    // RDKit✔️✔️:   mol.getPropIfPresent("_cxsmilesOutputIndex", sgroupOutputIndex);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::stringstream res;
    // RDKit✔️✔️:   // we need a map from original atom idx to output idx:
    // RDKit✔️✔️:   std::vector<unsigned int> revOrder(mol.getNumAtoms());
    // RDKit✔️✔️:   for (unsigned i = 0; i < atomOrder.size(); ++i) {
    // RDKit✔️✔️:     revOrder[atomOrder[i]] = i;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto &sg : sgs) {
    // RDKit✔️✔️:     if (sg.hasProp("TYPE") && sg.getProp<std::string>("TYPE") == "DAT") {
    // RDKit✔️✔️:       sg.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);
    // RDKit✔️✔️:       ++sgroupOutputIndex;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       res << "SgD:";
    // RDKit✔️✔️:       // we don't attempt to canonicalize the atom order because the user
    // RDKit✔️✔️:       // may ascribe some significance to the ordering of the atoms
    // RDKit✔️✔️:       for (const auto oaid : sg.getAtoms()) {
    // RDKit✔️✔️:         res << revOrder[oaid] << ",";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // remove the extra ",":
    // RDKit✔️✔️:       res.seekp(-1, res.cur);
    // RDKit✔️✔️:       res << ":";
    // RDKit✔️✔️:       std::string prop;
    // RDKit✔️✔️:       if (sg.getPropIfPresent("FIELDNAME", prop) && !prop.empty()) {
    // RDKit✔️✔️:         res << prop;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       res << ":";
    // RDKit✔️✔️:       std::vector<std::string> vprop;
    // RDKit✔️✔️:       if (sg.getPropIfPresent("DATAFIELDS", vprop) && !vprop.empty()) {
    // RDKit✔️✔️:         for (const auto &pv : vprop) {
    // RDKit✔️✔️:           res << pv << ",";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // remove the extra ",":
    // RDKit✔️✔️:         res.seekp(-1, res.cur);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       res << ":";
    // RDKit✔️✔️:       if (sg.getPropIfPresent("QUERYOP", prop) && !prop.empty()) {
    // RDKit✔️✔️:         res << prop;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       res << ":";
    // RDKit✔️✔️:       if (sg.getPropIfPresent("FIELDINFO", prop) && !prop.empty()) {
    // RDKit✔️✔️:         res << prop;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       res << ":";
    // RDKit✔️✔️:       if (sg.getPropIfPresent("FIELDTAG", prop) && !prop.empty()) {
    // RDKit✔️✔️:         res << prop;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       res << ":";
    // RDKit✔️✔️:       // FIX: do something about the coordinates
    // RDKit✔️✔️:       res << ",";  // only add a comma if we wrote something
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string resStr = res.str();
    // RDKit✔️✔️:   if (!resStr.empty() && resStr.back() == ',') {
    // RDKit✔️✔️:     resStr.pop_back();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   mol.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return resStr;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_sgroup_data_block
    // Behavior review: canonical typed SGroupData fields project the same
    // RDKit property types; legacy raw properties remain a fallback for CX
    // lowered records. Field and value bytes are appended literally because
    // this source function never calls quote_string. The empty-member case
    // deliberately retains the source seekp result, which removes the usual
    // atom/field delimiter. Output-index bookkeeping is computed without
    // mutation by write_sgroup_hierarchy on the same immutable prepared input.
    // Complexity review: one reverse-order vector and one output buffer match
    // the source's O(atom count + total members + output bytes) time and
    // O(atom count + output bytes) storage. Typed values are borrowed and
    // written once; no group, field, or full output clone is introduced.
    let positions = atom_positions(atom_order, record.topology.atoms.len());
    let mut output = PropertyText::new();
    for group in &record.topology.substance_groups {
        if !is_data_sgroup(group)? {
            continue;
        }
        if !output.is_empty() {
            output.push_byte(b',');
        }
        output.extend_bytes(("SgD:").as_ref());
        for (member_index, atom) in group.atoms().iter().enumerate() {
            if member_index != 0 {
                output.push_byte(b',');
            }
            // The pinned source reverse vector is zero-initialized, so an
            // SGroup atom omitted from a selected fragment maps to row zero.
            let position = positions[atom.index()].unwrap_or_default();
            write!(&mut output, "{position}").expect("writing to String cannot fail");
        }
        if !group.atoms().is_empty() {
            output.push_byte(b':');
        }
        output.extend_bytes((data_sgroup_value(group, "FIELDNAME")?).as_ref());
        output.push_byte(b':');

        let values = data_sgroup_values(group)?;
        {
            for (value_index, value) in values.iter().enumerate() {
                if value_index != 0 {
                    output.push_byte(b',');
                }
                output.extend_bytes((value).as_ref());
            }
        }
        output.push_byte(b':');
        output.extend_bytes((data_sgroup_value(group, "QUERYOP")?).as_ref());
        output.push_byte(b':');
        output.extend_bytes((data_sgroup_value(group, "FIELDINFO")?).as_ref());
        output.push_byte(b':');
        output.extend_bytes((data_sgroup_value(group, "FIELDTAG")?).as_ref());
        output.push_byte(b':');
    }
    Ok(output)
}

fn polymer_type(group: &SubstanceGroup) -> Result<Option<&'static str>, SmilesParseError> {
    // BEGIN RDKIT CPP VALUE sgroupTypemap
    // RDKit✔️🔝: const std::map<std::string, std::string> sgroupTypemap = {
    // RDKit✔️🔝:     {"n", "SRU"},   {"mon", "MON"}, {"mer", "MER"}, {"co", "COP"},
    // RDKit✔️🔝:     {"xl", "CRO"},  {"mod", "MOD"}, {"mix", "MIX"}, {"f", "FOR"},
    // RDKit✔️🔝:     {"any", "ANY"}, {"gen", "GEN"}, {"c", "COM"},   {"grf", "GRA"},
    // RDKit✔️🔝:     {"alt", "COP"}, {"ran", "COP"}, {"blk", "COP"}};
    // END RDKIT CPP VALUE sgroupTypemap
    // BEGIN RDKIT CPP FUNCTION get_sgroup_polymer_block (type selection)
    // RDKit✔️🔝: std::map<std::string, std::string> reverseTypemap;
    // RDKit✔️🔝: for (const auto &pr : SmilesParseOps::sgroupTypemap) {
    // RDKit✔️🔝:   if (reverseTypemap.find(pr.second) == reverseTypemap.end()) {
    // RDKit✔️🔝:     reverseTypemap[pr.second] = pr.first;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝: }
    // RDKit✔️🔝: std::string subtype;
    // RDKit✔️🔝: if (typ == "COP" && sg.getPropIfPresent("SUBTYPE", subtype)) {
    // RDKit✔️🔝:   if (subtype == "ALT") {
    // RDKit✔️🔝:     res << "alt";
    // RDKit✔️🔝:   } else if (subtype == "RAN") {
    // RDKit✔️🔝:     res << "ran";
    // RDKit✔️🔝:   } else if (subtype == "BLO") {
    // RDKit✔️🔝:     res << "blk";
    // RDKit✔️🔝:   } else {
    // RDKit✔️🔝:     res << reverseTypemap["COP"];
    // RDKit✔️🔝:   }
    // RDKit✔️🔝: } else {
    // RDKit✔️🔝:   res << reverseTypemap[typ];
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION get_sgroup_polymer_block (type selection)
    // A fixed enum match removes both source tree maps while preserving their
    // lexicographically first reverse spelling, including COP -> "alt".
    Ok(match group.kind() {
        SubstanceGroupKind::StructuralRepeatUnit => Some("n"),
        SubstanceGroupKind::Monomer => Some("mon"),
        SubstanceGroupKind::Mer => Some("mer"),
        SubstanceGroupKind::Copolymer => {
            let raw_subtype = if group.subtype().is_none() {
                group
                    .props()
                    .get(b"SUBTYPE".as_slice())
                    .map(cosmolkit_core::property_value_to_string)
                    .transpose()
                    .map_err(SmilesParseError::WriterProperty)?
            } else {
                None
            };
            match group
                .subtype()
                .map(PropertyText::as_bytes)
                .or_else(|| raw_subtype.as_ref().map(PropertyText::as_bytes))
            {
                Some(b"ALT") => Some("alt"),
                Some(b"RAN") => Some("ran"),
                Some(b"BLO") => Some("blk"),
                _ => Some("alt"),
            }
        }
        SubstanceGroupKind::Crosslink => Some("xl"),
        SubstanceGroupKind::Modification => Some("mod"),
        SubstanceGroupKind::Mixture => Some("mix"),
        SubstanceGroupKind::MixtureComponent => Some("c"),
        SubstanceGroupKind::Formulation => Some("f"),
        SubstanceGroupKind::AnyPolymer => Some("any"),
        SubstanceGroupKind::Graft => Some("grf"),
        SubstanceGroupKind::Generic(value) if value.as_bytes() == b"GEN" => Some("gen"),
        SubstanceGroupKind::Generic(value) if value.as_bytes() == b"COM" => Some("c"),
        _ => None,
    })
}

fn connection_text(group: &SubstanceGroup) -> Result<PropertyText, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION get_sgroup_polymer_block (connectivity)
    // RDKit✔️✔️: std::string connect;
    // RDKit✔️✔️: if (sg.getPropIfPresent("CONNECT", connect)) {
    // RDKit✔️✔️:   boost::algorithm::to_lower(connect);
    // RDKit✔️✔️:   res << connect;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_sgroup_polymer_block (connectivity)
    if let Some(connection) = group.connection() {
        return Ok(match connection {
            SGroupConnection::HeadToHead => PropertyText::from("hh"),
            SGroupConnection::HeadToTail => PropertyText::from("ht"),
            SGroupConnection::Either => PropertyText::from("eu"),
            SGroupConnection::Unknown(value) => {
                PropertyText::from(value.as_bytes().to_ascii_lowercase())
            }
        });
    }
    let raw = group
        .props()
        .get(b"CONNECT".as_slice())
        .map(cosmolkit_core::property_value_to_string)
        .transpose()
        .map_err(SmilesParseError::WriterProperty)?;
    Ok(raw
        .map(|value| PropertyText::from(value.as_bytes().to_ascii_lowercase()))
        .unwrap_or_default())
}

fn write_polymer_sgroups<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> Result<PropertyText, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION get_sgroup_polymer_block (field order)
    // RDKit✔️❌: for (const auto &sg : sgs) {
    // RDKit✔️❌:   std::string typ;
    // RDKit✔️❌:   if (sg.getPropIfPresent("TYPE", typ) &&
    // RDKit✔️❌:       reverseTypemap.find(typ) != reverseTypemap.end()) {
    // RDKit✔️❌:     sg.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);
    // RDKit✔️❌:     ++sgroupOutputIndex;
    // RDKit✔️❌:
    // RDKit✔️❌:     res << "Sg:";
    // RDKit✔️❌:     res << ":";
    // RDKit✔️❌:     for (const auto oaid : sg.getAtoms()) {
    // RDKit✔️❌:       res << revAtomOrder[oaid] << ",";
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // remove the extra ",":
    // RDKit✔️❌:     res.seekp(-1, res.cur);
    // RDKit✔️❌:     res << ":";
    // RDKit✔️❌:     std::string label;
    // RDKit✔️❌:     if (sg.getPropIfPresent("LABEL", label)) {
    // RDKit✔️❌:       res << label;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     res << ":";
    // RDKit✔️❌:     std::string connect;
    // RDKit✔️❌:     if (sg.getPropIfPresent("CONNECT", connect)) {
    // RDKit✔️❌:       boost::algorithm::to_lower(connect);
    // RDKit✔️❌:       res << connect;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     res << ":";
    // END RDKIT CPP FUNCTION get_sgroup_polymer_block (field order)
    // The Rust path has source-equivalent order but currently allocates
    // per-group atom/crossing vectors and a formatted block before joining.
    let atom_positions = atom_positions(atom_order, record.topology.atoms.len());
    let mut blocks = Vec::new();
    for group in &record.topology.substance_groups {
        let Some(kind) = polymer_type(group)? else {
            continue;
        };
        let atoms = group
            .atoms()
            .iter()
            .map(|atom| atom_positions[atom.index()].unwrap_or_default())
            .map(|position| position.to_string())
            .collect::<Vec<_>>();
        if atoms.is_empty() {
            continue;
        }
        // BEGIN RDKIT CPP FUNCTION get_sgroup_polymer_block (crossing bonds)
        // RDKit✔️❌: std::vector<unsigned int> headCrossings;
        // RDKit✔️❌: if (sg.getPropIfPresent("XBHEAD", headCrossings) &&
        // RDKit✔️❌:     headCrossings.size() > 1) {
        // RDKit✔️❌:   for (auto v : headCrossings) {
        // RDKit✔️❌:     res << bondOrder[v] << ",";
        // RDKit✔️❌:   }
        // RDKit✔️❌:   // remove the extra ",":
        // RDKit✔️❌:   res.seekp(-1, res.cur);
        // RDKit✔️❌: }
        // RDKit✔️❌: res << ":";
        // RDKit✔️❌: std::vector<unsigned int> tailCrossings;
        // RDKit✔️❌: if (sg.getPropIfPresent("XBCORR", tailCrossings) &&
        // RDKit✔️❌:     tailCrossings.size() > 2) {
        // RDKit✔️❌:   for (unsigned int i = 1; i < tailCrossings.size(); i += 2) {
        // RDKit✔️❌:     res << bondOrder[tailCrossings[i]] << ",";
        // RDKit✔️❌:   }
        // RDKit✔️❌:   // remove the extra ",":
        // RDKit✔️❌:   res.seekp(-1, res.cur);
        // RDKit✔️❌: }
        // RDKit✔️❌: res << ":";
        // END RDKIT CPP FUNCTION get_sgroup_polymer_block (crossing bonds)
        // The direct source lookup is intentionally not an inverse mapping.
        // Rust retains extra per-field collection allocations before joining.
        let crossing_position = |bond: BondId| {
            bond_order
                .get(bond.index())
                .copied()
                .map(|ordered_bond| ordered_bond.index())
                .ok_or_else(|| {
                    SmilesParseError::Model(format!(
                        "SGroup {} crossing bond {} cannot index the CX bond order",
                        group.id().index(),
                        bond.index()
                    ))
                })
        };
        let head_text = if group.head_crossing_bonds().len() > 1 {
            group
                .head_crossing_bonds()
                .iter()
                .copied()
                .map(crossing_position)
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .map(|position| position.to_string())
                .collect::<Vec<_>>()
                .join(",")
        } else {
            String::new()
        };
        let tail_text = if group.crossing_bond_correspondence().len() > 2 {
            group
                .crossing_bond_correspondence()
                .iter()
                .skip(1)
                .step_by(2)
                .copied()
                .map(crossing_position)
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .map(|position| position.to_string())
                .collect::<Vec<_>>()
                .join(",")
        } else {
            String::new()
        };
        let mut block = PropertyText::new();
        write!(&mut block, "Sg:{kind}:{}:", atoms.join(","))
            .expect("byte output formatting cannot fail");
        let raw_label = if group.label().is_none() {
            group
                .props()
                .get(b"LABEL".as_slice())
                .map(cosmolkit_core::property_value_to_string)
                .transpose()
                .map_err(SmilesParseError::WriterProperty)?
        } else {
            None
        };
        if let Some(label) = group.label().or(raw_label.as_ref()) {
            block.extend_bytes(label.as_bytes());
        }
        block.push_byte(b':');
        block.extend_bytes(connection_text(group)?.as_bytes());
        write!(&mut block, ":{head_text}:{tail_text}:")
            .expect("byte output formatting cannot fail");
        blocks.push(block);
    }
    let mut output = PropertyText::new();
    for (index, block) in blocks.into_iter().enumerate() {
        if index != 0 {
            output.push_byte(b',');
        }
        output.extend_bytes(block.as_bytes());
    }
    Ok(output)
}

fn sgroup_index(group: &SubstanceGroup) -> Result<usize, SmilesParseError> {
    // RDKit❗❌: unsigned int sgidx = sg.getIndexInMol();
    // RDKit❗❌: sg.getPropIfPresent("index", sgidx);
    // String-only detached SGroup metadata requires one temporary tagged copy;
    // conversion still uses the canonical from_rdvalue owner and propagates errors.
    group
        .props()
        .get(b"index".as_slice())
        .map(|value| cosmolkit_core::property_value_to_uint(value))
        .transpose()
        .map_err(SmilesParseError::WriterNumeric)
        .map(|value| {
            value
                .map(|value| value as usize)
                .unwrap_or_else(|| group.id().index())
        })
}

fn write_sgroup_hierarchy<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    include_data: bool,
    include_polymer: bool,
) -> Result<String, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION get_sgroup_hierarchy_block
    // RDKit✔️✔️: std::string get_sgroup_hierarchy_block(const ROMol &mol) {
    // RDKit✔️✔️:   const auto &sgs = getSubstanceGroups(mol);
    // RDKit✔️✔️:   if (sgs.empty()) {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::stringstream res;
    // RDKit✔️✔️:   // we need a map from sgroup index to output index;
    // RDKit✔️✔️:   std::map<unsigned int, unsigned int> sgroupOrder;
    // RDKit✔️✔️:   bool parentPresent = false;
    // RDKit✔️✔️:   for (const auto &sg : sgs) {
    // RDKit✔️✔️:     if (sg.hasProp("_cxsmilesOutputIndex")) {
    // RDKit✔️✔️:       unsigned int sgidx = sg.getIndexInMol();
    // RDKit✔️✔️:       sg.getPropIfPresent("index", sgidx);
    // RDKit✔️✔️:       sgroupOrder[sgidx] = sg.getProp<unsigned int>("_cxsmilesOutputIndex");
    // RDKit✔️✔️:       sg.clearProp("_cxsmilesOutputIndex");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (sg.hasProp("PARENT")) {
    // RDKit✔️✔️:       parentPresent = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (parentPresent) {
    // RDKit✔️✔️:     // now loop over them and add the information
    // RDKit✔️✔️:     std::map<unsigned int, std::vector<unsigned int>> accum;
    // RDKit✔️✔️:     for (const auto &sg : sgs) {
    // RDKit✔️✔️:       unsigned pidx;
    // RDKit✔️✔️:       if (sg.getPropIfPresent("PARENT", pidx) &&
    // RDKit✔️✔️:           sgroupOrder.find(pidx) != sgroupOrder.end()) {
    // RDKit✔️✔️:         unsigned int sgidx = sg.getIndexInMol();
    // RDKit✔️✔️:         sg.getPropIfPresent("index", sgidx);
    // RDKit✔️✔️:         if (sgroupOrder.find(sgidx) != sgroupOrder.end()) {
    // RDKit✔️✔️:           accum[sgroupOrder[pidx]].push_back(sgroupOrder[sgidx]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!accum.empty()) {
    // RDKit✔️✔️:       res << "SgH:";
    // RDKit✔️✔️:       for (const auto &pr : accum) {
    // RDKit✔️✔️:         res << pr.first << ":";
    // RDKit✔️✔️:         for (auto v : pr.second) {
    // RDKit✔️✔️:           res << v << ".";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // remove the extra ".":
    // RDKit✔️✔️:         res.seekp(-1, res.cur);
    // RDKit✔️✔️:         res << ",";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::string resStr = res.str();
    // RDKit✔️✔️:     while (!resStr.empty() && resStr.back() == ',') {
    // RDKit✔️✔️:       resStr.pop_back();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return resStr;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_sgroup_hierarchy_block
    let mut output_indices = BTreeMap::new();
    let mut next = 0;
    if include_data {
        for group in &record.topology.substance_groups {
            if is_data_sgroup(group)? {
                output_indices.insert(sgroup_index(group)?, next);
                next += 1;
            }
        }
    }
    if include_polymer {
        for group in &record.topology.substance_groups {
            if polymer_type(group)?.is_some() {
                output_indices.insert(sgroup_index(group)?, next);
                next += 1;
            }
        }
    }
    let mut hierarchy = BTreeMap::<usize, Vec<usize>>::new();
    for group in &record.topology.substance_groups {
        let Some(child) = output_indices.get(&sgroup_index(group)?).copied() else {
            continue;
        };
        // RDKit resolves the parent through the `PARENT` property, which
        // holds the parent's source `index`, and looks it up in the map keyed
        // by that same `index` space. The typed `parent` stores the parent's
        // `SubstanceGroupId`, so translate it through the parent row's own
        // `index` rather than using the row id directly. Topology validation
        // fixes each SubstanceGroupId to its row position, so direct indexing
        // avoids an additional tree map. Fall back to the preserved `PARENT`
        // property when only that is available.
        let parent_key = if let Some(parent_id) = group.parent() {
            record
                .topology
                .substance_groups
                .get(parent_id.index())
                .map(sgroup_index)
                .transpose()?
        } else {
            group
                .props()
                .get(b"PARENT".as_slice())
                .map(|value| cosmolkit_core::property_value_to_uint(value))
                .transpose()
                .map_err(SmilesParseError::WriterNumeric)?
                .map(|value| value as usize)
        };
        let Some(parent) = parent_key.and_then(|parent| output_indices.get(&parent).copied())
        else {
            continue;
        };
        hierarchy.entry(parent).or_default().push(child);
    }
    if hierarchy.is_empty() {
        Ok(String::new())
    } else {
        // Complexity review: like the source, this uses one ordered map for
        // source-to-output indices, one ordered parent accumulator, and one
        // output buffer. Validated typed parent IDs provide O(1) row lookup;
        // each selected group and emitted child is visited once apart from the
        // O(log n) ordered-map operations.
        let mut output = String::from("SgH:");
        for (parent_position, (parent, children)) in hierarchy.into_iter().enumerate() {
            if parent_position != 0 {
                output.push(',');
            }
            write!(&mut output, "{parent}:").expect("writing to String cannot fail");
            for (child_position, child) in children.into_iter().enumerate() {
                if child_position != 0 {
                    output.push('.');
                }
                write!(&mut output, "{child}").expect("writing to String cannot fail");
            }
        }
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SmilesParseParams, parse_smiles};
    use cosmolkit_model::{
        Atom, AtomSpec, BondSpec, CoordinateBlock, MoleculeProperties, TopologyBlock,
    };
    use cosmolkit_types::Element;

    fn parse(input: &str) -> SmilesRecord {
        parse_smiles(input, &SmilesParseParams::default()).expect("parse CXSMILES")
    }

    fn write_noncanonical(input: &str) -> String {
        String::from_utf8(
            write_cx_smiles_with_params(
                &parse(input),
                &CxSmilesWriteParams {
                    smiles: SmilesWriteParams {
                        canonical: false,
                        ..Default::default()
                    },
                    ..Default::default()
                },
            )
            .expect("write CXSMILES")
            .into_bytes(),
        )
        .expect("original fixture CX text is UTF-8")
    }

    fn record_with_coordinate_sets(
        conformers_2d: Vec<Conformer2D>,
        conformers_3d: Vec<Conformer3D>,
    ) -> SmilesRecord {
        let mut record = parse("CC");
        record.coordinates = CoordinateBlock {
            conformers_2d,
            conformers_3d,
            ..Default::default()
        };
        record
    }

    #[test]
    fn cx_coordinate_selection_auto_returns_none_or_the_only_stored_set() {
        let empty = record_with_coordinate_sets(Vec::new(), Vec::new());
        let empty_before = empty.clone();
        assert!(matches!(
            coordinate_source(&empty, CxCoordinateSelection::Auto),
            Ok(None)
        ));
        assert_eq!(empty, empty_before);

        let only_2d = record_with_coordinate_sets(
            vec![Conformer2D::new(17, vec![[1.0, 2.0], [3.0, 4.0]])],
            Vec::new(),
        );
        let only_2d_before = only_2d.clone();
        match coordinate_source(&only_2d, CxCoordinateSelection::Auto)
            .unwrap()
            .expect("one stored 2D set is selected")
        {
            CoordinateSource::TwoD(conformer) => assert_eq!(conformer.id(), 17),
            CoordinateSource::ThreeD(_) => panic!("Auto selected the wrong dimension"),
        }
        assert_eq!(only_2d, only_2d_before);

        let only_3d = record_with_coordinate_sets(
            Vec::new(),
            vec![Conformer3D::new(
                23,
                vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                true,
            )],
        );
        let only_3d_before = only_3d.clone();
        match coordinate_source(&only_3d, CxCoordinateSelection::Auto)
            .unwrap()
            .expect("one stored 3D set is selected")
        {
            CoordinateSource::ThreeD(conformer) => assert_eq!(conformer.id(), 23),
            CoordinateSource::TwoD(_) => panic!("Auto selected the wrong dimension"),
        }
        assert_eq!(only_3d, only_3d_before);
    }

    #[test]
    fn cx_coordinate_selection_auto_rejects_multiple_sets_and_reversed_mixed_order() {
        let multiple_2d = record_with_coordinate_sets(
            vec![
                Conformer2D::new(4, vec![[0.0, 0.0], [1.0, 0.0]]),
                Conformer2D::new(9, vec![[2.0, 0.0], [3.0, 0.0]]),
            ],
            Vec::new(),
        );
        let multiple_3d = record_with_coordinate_sets(
            Vec::new(),
            vec![
                Conformer3D::new(4, vec![[0.0, 0.0, 1.0], [1.0, 0.0, 2.0]], true),
                Conformer3D::new(9, vec![[2.0, 0.0, 3.0], [3.0, 0.0, 4.0]], true),
            ],
        );

        for (record, expected_2d, expected_3d) in [(multiple_2d, 2, 0), (multiple_3d, 0, 2)] {
            let before = record.clone();
            assert_eq!(record.coordinates.conformers_2d.len(), expected_2d);
            assert_eq!(record.coordinates.conformers_3d.len(), expected_3d);
            assert!(matches!(
                coordinate_source(&record, CxCoordinateSelection::Auto),
                Err(SmilesParseError::AmbiguousCoordinateSelection {
                    two_d_count,
                    three_d_count,
                }) if two_d_count == expected_2d && three_d_count == expected_3d
            ));
            assert_eq!(record, before);
        }

        for input in ["CC |(0,0;1,0)(0,0,1;1,0,2)|", "CC |(0,0,1;1,0,2)(0,0;1,0)|"] {
            let record = parse(input);
            let before = record.clone();
            // Native parse_coords appends Point3D-backed conformers and stores
            // is3D separately; both source rows remain in append order.
            assert_eq!(record.coordinates.conformers_2d.len(), 0);
            assert_eq!(record.coordinates.conformers_3d.len(), 2);
            assert!(matches!(
                coordinate_source(&record, CxCoordinateSelection::Auto),
                Err(SmilesParseError::AmbiguousCoordinateSelection {
                    two_d_count: 0,
                    three_d_count: 2,
                })
            ));
            assert_eq!(record, before, "selection must preserve {input}");
        }
    }

    #[test]
    fn cx_coordinate_selection_explicit_ids_are_dimension_scoped_and_not_positions() {
        let shared_id = record_with_coordinate_sets(
            vec![Conformer2D::new(5, vec![[10.0, 20.0], [30.0, 40.0]])],
            vec![Conformer3D::new(
                5,
                vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                true,
            )],
        );
        let shared_id_before = shared_id.clone();
        assert!(matches!(
            coordinate_source(&shared_id, CxCoordinateSelection::Auto),
            Err(SmilesParseError::AmbiguousCoordinateSelection {
                two_d_count: 1,
                three_d_count: 1,
            })
        ));
        match coordinate_source(&shared_id, CxCoordinateSelection::TwoD { id: 5 })
            .unwrap()
            .expect("explicit 2D ID exists")
        {
            CoordinateSource::TwoD(conformer) => {
                assert_eq!(conformer.coordinates()[0], [10.0, 20.0]);
            }
            CoordinateSource::ThreeD(_) => panic!("explicit 2D selected the wrong dimension"),
        }
        match coordinate_source(&shared_id, CxCoordinateSelection::ThreeD { id: 5 })
            .unwrap()
            .expect("explicit 3D ID exists")
        {
            CoordinateSource::ThreeD(conformer) => {
                assert_eq!(conformer.coordinates()[0], [1.0, 2.0, 3.0]);
            }
            CoordinateSource::TwoD(_) => panic!("explicit 3D selected the wrong dimension"),
        }
        assert_eq!(shared_id, shared_id_before);

        let noncontiguous_ids = record_with_coordinate_sets(
            Vec::new(),
            vec![
                Conformer3D::new(7, vec![[7.0, 0.0, 0.0], [7.0, 1.0, 0.0]], true),
                Conformer3D::new(21, vec![[21.0, 0.0, 0.0], [21.0, 1.0, 0.0]], true),
            ],
        );
        let noncontiguous_before = noncontiguous_ids.clone();
        assert!(matches!(
            coordinate_source(&noncontiguous_ids, CxCoordinateSelection::ThreeD { id: 1 }),
            Err(SmilesParseError::MissingCoordinateSelection {
                selection: CxCoordinateSelection::ThreeD { id: 1 }
            })
        ));
        match coordinate_source(&noncontiguous_ids, CxCoordinateSelection::ThreeD { id: 21 })
            .unwrap()
            .expect("stored noncontiguous ID resolves by ID")
        {
            CoordinateSource::ThreeD(conformer) => {
                assert_eq!(conformer.coordinates()[0], [21.0, 0.0, 0.0]);
            }
            CoordinateSource::TwoD(_) => panic!("explicit 3D selected the wrong dimension"),
        }
        assert_eq!(noncontiguous_ids, noncontiguous_before);

        let empty = record_with_coordinate_sets(Vec::new(), Vec::new());
        let empty_before = empty.clone();
        assert!(matches!(
            coordinate_source(&empty, CxCoordinateSelection::TwoD { id: 99 }),
            Err(SmilesParseError::MissingCoordinateSelection {
                selection: CxCoordinateSelection::TwoD { id: 99 }
            })
        ));
        let other_dimension_only = record_with_coordinate_sets(
            Vec::new(),
            vec![Conformer3D::new(
                3,
                vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                true,
            )],
        );
        let other_dimension_before = other_dimension_only.clone();
        assert!(matches!(
            coordinate_source(&other_dimension_only, CxCoordinateSelection::TwoD { id: 3 }),
            Err(SmilesParseError::MissingCoordinateSelection {
                selection: CxCoordinateSelection::TwoD { id: 3 }
            })
        ));
        assert_eq!(empty, empty_before);
        assert_eq!(other_dimension_only, other_dimension_before);
    }

    fn coordinate_consumer_record() -> SmilesRecord {
        let mut record = parse("C[C@H](O)Cl");
        record.coordinates = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                101,
                vec![
                    [3.9163, 5.4767],
                    [3.9163, 3.9367],
                    [2.5826, 3.1667],
                    [5.25, 3.1667],
                ],
            )],
            conformers_3d: vec![Conformer3D::new(
                202,
                vec![
                    [-3.9163, 5.4767, 0.0],
                    [-3.9163, 3.9367, 1.0],
                    [-2.5826, 3.1667, 2.0],
                    [-5.25, 3.1667, 3.0],
                ],
                true,
            )],
            ..Default::default()
        };
        record
    }

    fn coordinate_consumer_params(
        selection: CxCoordinateSelection,
        fields: CxSmilesFields,
    ) -> CxSmilesWriteParams {
        CxSmilesWriteParams {
            smiles: SmilesWriteParams {
                canonical: true,
                isomeric_smiles: true,
                clean_stereo: true,
                ..Default::default()
            },
            fields,
            coordinate_selection: selection,
        }
    }

    #[test]
    fn cx_coordinate_consumer_explicit_selection_drives_text_and_wedge_geometry() {
        // Pinned RDKit 2026.03.1 MolToCXSmiles outputs with one selected
        // conformer, canonical/isomeric/cleanStereo true, and COORDS|BOND_CFG.
        // The two geometries intentionally produce opposite wedge directions.
        let record = coordinate_consumer_record();
        let before = record.clone();
        let fields = CxSmilesFields::COORDS | CxSmilesFields::BOND_CFG;

        assert_eq!(
            write_cx_smiles_with_params(
                &record,
                &coordinate_consumer_params(CxCoordinateSelection::TwoD { id: 101 }, fields),
            )
            .unwrap(),
            "C[C@H](O)Cl |(3.9163,5.4767,;3.9163,3.9367,;2.5826,3.1667,;5.25,3.1667,),wU:1.0|"
                .into()
        );
        assert_eq!(write_cx_smiles_with_params(
                &record,
                &coordinate_consumer_params(CxCoordinateSelection::ThreeD { id: 202 }, fields),
            )
            .unwrap(),
            "C[C@H](O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,1;-2.5826,3.1667,2;-5.25,3.1667,3),wD:1.0|".into()
        );
        assert_eq!(
            record, before,
            "coordinate consumers must preserve the input"
        );
    }

    #[test]
    fn cx_coordinate_consumer_auto_ambiguity_is_structured_and_preserves_input() {
        let record = coordinate_consumer_record();
        let before = record.clone();
        let error = write_cx_smiles_with_params(
            &record,
            &coordinate_consumer_params(
                CxCoordinateSelection::Auto,
                CxSmilesFields::COORDS | CxSmilesFields::BOND_CFG,
            ),
        )
        .expect_err("multiple stored sets require explicit selection under CK-COORD-002");
        assert!(matches!(
            error,
            SmilesParseError::AmbiguousCoordinateSelection {
                two_d_count: 1,
                three_d_count: 1,
            }
        ));
        assert_eq!(record, before, "failed selection must preserve the input");
    }

    #[test]
    fn cx_coordinate_consumer_disabled_coords_bypasses_missing_selection_and_geometry() {
        let record = coordinate_consumer_record();
        let before = record.clone();
        assert_eq!(
            write_cx_smiles_with_params(
                &record,
                &coordinate_consumer_params(
                    CxCoordinateSelection::TwoD { id: 999 },
                    CxSmilesFields::BOND_CFG,
                ),
            )
            .unwrap(),
            "C[C@H](O)Cl".into()
        );
        assert_eq!(
            record, before,
            "coordinate-disabled export must preserve input"
        );
    }

    fn stereo_group_with_write_id(
        kind: StereoGroupKind,
        read_id: u32,
        write_id: u32,
    ) -> (StereoGroup, Vec<usize>) {
        let mut group = StereoGroup::new(kind, Vec::new(), Vec::new())
            .expect("valid distinct stereo members")
            .with_id(read_id);
        set_stereo_group_write_id(&mut group, write_id);
        (group, Vec::new())
    }

    fn stereo_group_with_members_and_write_id(
        kind: StereoGroupKind,
        read_id: u32,
        write_id: u32,
        atoms: Vec<AtomId>,
        bonds: Vec<BondId>,
    ) -> StereoGroup {
        let mut group = StereoGroup::new(kind, atoms, bonds)
            .expect("valid distinct stereo members")
            .with_id(read_id);
        set_stereo_group_write_id(&mut group, write_id);
        group
    }

    fn write_extensions_for_test(
        record: &SmilesRecord,
        fields: CxSmilesFields,
    ) -> Result<String, SmilesParseError> {
        // These owner tests target the extension phase with explicit identity
        // output maps. The detached base writer does not accept atrop bond
        // stereo, while the pinned CX extension phase consumes precomputed maps.
        let atom_order = (0..record.topology.atoms.len())
            .map(AtomId::new)
            .collect::<Vec<_>>();
        let bond_order = (0..record.topology.bonds.len())
            .map(BondId::new)
            .collect::<Vec<_>>();
        write_extensions_with_orders_for_test(record, fields, &atom_order, &bond_order)
    }

    fn write_extensions_with_orders_for_test(
        record: &SmilesRecord,
        fields: CxSmilesFields,
        atom_order: &[AtomId],
        bond_order: &[BondId],
    ) -> Result<String, SmilesParseError> {
        write_cx_extensions(
            record,
            fields,
            atom_order,
            bond_order,
            None,
            CxCoordinateSelection::Auto,
        )
        .map(|text| {
            String::from_utf8(text.into_bytes()).expect("original fixture CX extension is UTF-8")
        })
    }

    fn atrop_record(stereo: BondStereo, reverse_carrier_bonds: bool) -> SmilesRecord {
        let atom_specs = (0..4)
            .map(|_| {
                AtomSpec::new(Element::from_atomic_number(6).expect("carbon"))
                    .with_no_implicit(true)
            })
            .collect::<Vec<_>>();
        let atoms = atom_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect();
        let first_carrier = if reverse_carrier_bonds {
            (2, 0)
        } else {
            (0, 2)
        };
        let second_carrier = if reverse_carrier_bonds {
            (3, 1)
        } else {
            (1, 3)
        };
        let bonds = [
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single).with_stereo(stereo),
            BondSpec::new(
                AtomId::new(first_carrier.0),
                AtomId::new(first_carrier.1),
                BondOrder::Single,
            ),
            BondSpec::new(
                AtomId::new(second_carrier.0),
                AtomId::new(second_carrier.1),
                BondOrder::Single,
            ),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
        .collect();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("valid atrop writer fixture");
        SmilesRecord {
            topology,
            coordinates: CoordinateBlock::default(),
            properties: MoleculeProperties::default(),
        }
    }

    #[test]
    fn stereo_write_id_assignment_keeps_unique_namespaces_and_fills_duplicates_and_holes() {
        let mut groups = vec![
            stereo_group_with_write_id(StereoGroupKind::Or, 7, 5),
            stereo_group_with_write_id(StereoGroupKind::And, 2, 5),
            stereo_group_with_write_id(StereoGroupKind::Or, 11, 5),
            stereo_group_with_write_id(StereoGroupKind::Or, 13, 1),
            stereo_group_with_write_id(StereoGroupKind::Or, 17, 3),
            stereo_group_with_write_id(StereoGroupKind::Or, 19, 0),
            stereo_group_with_write_id(StereoGroupKind::And, 23, 5),
            stereo_group_with_write_id(StereoGroupKind::And, 29, 0),
            stereo_group_with_write_id(StereoGroupKind::Absolute, 31, 77),
        ];
        let read_ids = groups
            .iter()
            .map(|(group, _)| group.id())
            .collect::<Vec<_>>();

        assign_stereo_group_ids(&mut groups);

        assert_eq!(
            groups
                .iter()
                .map(|(group, _)| stereo_group_write_id(group))
                .collect::<Vec<_>>(),
            vec![5, 5, 2, 1, 3, 4, 1, 2, 77]
        );
        assert_eq!(
            groups
                .iter()
                .map(|(group, _)| group.id())
                .collect::<Vec<_>>(),
            read_ids,
            "assignment must never forward or overwrite read IDs"
        );
    }

    #[test]
    fn stereo_write_id_assignment_accepts_an_empty_group_list() {
        let mut groups = Vec::<(StereoGroup, Vec<usize>)>::new();
        assign_stereo_group_ids(&mut groups);
        assert!(groups.is_empty());
    }

    #[test]
    fn writes_coordinates_labels_values_radicals_and_atom_properties_in_source_order() {
        let input = "*C |$foo;$,$_AV:x;y$,atomProp:0.a&#46;b.c&#46;d,^1:1,(1,2,3;3,4,0)|";
        assert_eq!(
            write_noncanonical(input),
            "*[CH2] |(1,2,3;3,4,),$foo;$,$_AV:x;y$,^1:1,atomProp:0.a&#46;b.c&#46;d|"
        );
    }

    #[test]
    fn canonical_output_indices_follow_final_atom_and_bond_order() {
        for (input, expected) in [
            ("OC |$oxygen;carbon$|", "CO |$carbon;oxygen$|"),
            ("OC |$_AV:o;c$|", "CO |$_AV:c;o$|"),
            ("OC |atomProp:0.x.y|", "CO |atomProp:1.x.y|"),
            ("OC |^1:0|", "C[O] |^1:1|"),
            ("OC |w:0.0|", "CO |w:1.0|"),
            ("CC |C:0.0|", "C[CH4] |C:1.0|"),
            ("CC |H:0.0|", "CC |H:0.0|"),
            ("CC |Z:0|", "C~C |Z:0|"),
        ] {
            assert_eq!(
                write_cx_smiles(&parse(input)).unwrap(),
                expected.into(),
                "{input}"
            );
        }
    }

    #[test]
    fn writes_enhanced_stereo_with_source_sorting_and_fresh_write_ids() {
        let input = "F[C@H](Cl)Br.O[C@@H](N)I |&7:1,o2:5|";
        assert_eq!(
            write_cx_smiles(&parse(input)).unwrap(),
            "F[C@H](Cl)Br.N[C@H](O)I |o1:5,&1:1|".into()
        );
        assert_eq!(
            write_noncanonical(input),
            "F[C@H](Cl)Br.O[C@@H](N)I |o1:5,&1:1|"
        );
    }

    #[test]
    fn enhanced_stereo_default_write_ids_do_not_forward_read_ids() {
        let mut record = parse("CC");
        record.topology.stereo_groups = vec![
            stereo_group_with_members_and_write_id(
                StereoGroupKind::And,
                2,
                0,
                vec![AtomId::new(1)],
                Vec::new(),
            ),
            stereo_group_with_members_and_write_id(
                StereoGroupKind::And,
                7,
                0,
                vec![AtomId::new(0)],
                Vec::new(),
            ),
        ];
        let before = record.clone();

        assert_eq!(
            write_extensions_for_test(&record, CxSmilesFields::ENHANCED_STEREO).unwrap(),
            "|&1:0,&2:1|"
        );
        assert_eq!(record, before, "extension writing must not mutate input");
    }

    #[test]
    fn enhanced_stereo_sorts_explicit_write_ids_before_mapped_indices() {
        let mut record = parse("CC");
        record.topology.stereo_groups = vec![
            stereo_group_with_members_and_write_id(
                StereoGroupKind::And,
                7,
                7,
                vec![AtomId::new(0)],
                Vec::new(),
            ),
            stereo_group_with_members_and_write_id(
                StereoGroupKind::And,
                2,
                2,
                vec![AtomId::new(1)],
                Vec::new(),
            ),
        ];

        assert_eq!(
            write_extensions_for_test(&record, CxSmilesFields::ENHANCED_STEREO).unwrap(),
            "|&2:1,&7:0|"
        );
    }

    #[test]
    fn enhanced_stereo_duplicate_write_ids_break_ties_by_mapped_indices() {
        let mut record = parse("CC");
        record.topology.stereo_groups = vec![
            stereo_group_with_members_and_write_id(
                StereoGroupKind::And,
                11,
                5,
                vec![AtomId::new(0)],
                Vec::new(),
            ),
            stereo_group_with_members_and_write_id(
                StereoGroupKind::And,
                13,
                5,
                vec![AtomId::new(1)],
                Vec::new(),
            ),
        ];
        let atom_order = [AtomId::new(1), AtomId::new(0)];
        let bond_order = [BondId::new(0)];

        assert_eq!(
            write_extensions_with_orders_for_test(
                &record,
                CxSmilesFields::ENHANCED_STEREO,
                &atom_order,
                &bond_order,
            )
            .unwrap(),
            "|&5:0,&1:1|"
        );
    }

    #[test]
    fn enhanced_stereo_atrop_bond_group_uses_shared_wedge_assignments() {
        let mut record = atrop_record(BondStereo::AtropCw, true);
        record.topology.stereo_groups = vec![stereo_group_with_members_and_write_id(
            StereoGroupKind::Or,
            7,
            0,
            Vec::new(),
            vec![BondId::new(0)],
        )];
        let before = record.clone();

        assert_eq!(
            write_extensions_for_test(
                &record,
                CxSmilesFields::BOND_CFG | CxSmilesFields::ENHANCED_STEREO,
            )
            .unwrap(),
            "|wU:1.2,o1:1|"
        );
        assert_eq!(record, before, "extension writing must not mutate input");
    }

    #[test]
    fn enhanced_stereo_mixed_group_maps_members_and_atrop_endpoints() {
        let mut record = atrop_record(BondStereo::AtropCw, true);
        record.topology.stereo_groups = vec![stereo_group_with_members_and_write_id(
            StereoGroupKind::And,
            7,
            0,
            vec![AtomId::new(3)],
            vec![BondId::new(0)],
        )];
        let before = record.clone();

        assert_eq!(
            write_extensions_for_test(
                &record,
                CxSmilesFields::BOND_CFG | CxSmilesFields::ENHANCED_STEREO,
            )
            .unwrap(),
            "|wU:1.2,&1:1,3|"
        );
        assert_eq!(record, before, "extension writing must not mutate input");
    }

    #[test]
    fn enhanced_stereo_omits_a_group_with_no_mapped_atoms() {
        let mut record = parse("CC");
        record.topology.stereo_groups = vec![stereo_group_with_members_and_write_id(
            StereoGroupKind::Or,
            7,
            4,
            Vec::new(),
            vec![BondId::new(0)],
        )];

        assert_eq!(
            write_extensions_for_test(&record, CxSmilesFields::ENHANCED_STEREO).unwrap(),
            ""
        );
    }

    #[test]
    fn cx_wedges_with_2d_coordinates_assign_and_write_tetrahedral_stereo() {
        for (input, expected) in [
            (
                "CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wU:1.0|",
                "C[C@@H](O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wU:1.0|",
            ),
            (
                "CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wD:1.0|",
                "C[C@H](O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wD:1.0|",
            ),
        ] {
            assert_eq!(
                write_cx_smiles(&parse(input)).unwrap(),
                expected.into(),
                "{input}"
            );
        }
        assert_eq!(
            write_cx_smiles(&parse("CC |(0,0,;1,0,),wU:0.0|")).unwrap(),
            "CC |(0,0,;1,0,)|".into()
        );
    }

    #[test]
    fn cx_wedge_flags_require_an_enabled_present_coordinate_conformer() {
        let with_coordinates = parse("CC |(0,0,;1,0,),wU:0.0|");
        let with_coordinates_before = with_coordinates.clone();
        assert_eq!(
            write_extensions_for_test(
                &with_coordinates,
                CxSmilesFields::COORDS | CxSmilesFields::BOND_CFG,
            )
            .expect("coordinates enable wedge output"),
            "|(0,0,;1,0,),wU:0.0|"
        );
        assert_eq!(
            write_extensions_for_test(&with_coordinates, CxSmilesFields::BOND_CFG)
                .expect("stored coordinates are not passed without CX_COORDS"),
            ""
        );
        assert_eq!(with_coordinates, with_coordinates_before);

        let without_coordinates = parse("CC |wU:0.0|");
        let without_coordinates_before = without_coordinates.clone();
        assert_eq!(
            write_extensions_for_test(
                &without_coordinates,
                CxSmilesFields::COORDS | CxSmilesFields::BOND_CFG,
            )
            .expect("CX_COORDS without a conformer does not enable wedge output"),
            ""
        );
        assert_eq!(without_coordinates, without_coordinates_before);
    }

    #[test]
    fn cx_explicit_bond_direction_precedes_molfile_bond_config() {
        let mut record = parse("CC |(0,0,;1,0,),wD:0.0|");
        record.topology.bonds[0]
            .set_prop("_MolFileBondCfg", "1")
            .expect("valid molfile bond configuration property");
        let before = record.clone();

        assert_eq!(
            write_extensions_for_test(&record, CxSmilesFields::COORDS | CxSmilesFields::BOND_CFG,)
                .expect("explicit direction is kept before config fallback"),
            "|(0,0,;1,0,),wD:0.0|"
        );
        assert_eq!(record, before, "extension writing must not mutate input");
    }

    #[test]
    fn cx_no_coordinate_atrop_assignments_write_reoriented_carrier_bonds() {
        for (stereo, expected, fields) in [
            (BondStereo::AtropCw, "|wU:1.2|", CxSmilesFields::BOND_CFG),
            (
                BondStereo::AtropCcw,
                "|wU:0.1|",
                CxSmilesFields::BOND_ATROPISOMER,
            ),
        ] {
            let record = atrop_record(stereo, true);
            let before = record.clone();
            assert_eq!(
                write_extensions_for_test(&record, fields)
                    .expect("no-coordinate atrop assignment is serializable"),
                expected,
                "stereo={stereo:?}, fields={}",
                fields.bits()
            );
            assert_eq!(record, before, "extension writing must not mutate input");
        }

        let both_flags = atrop_record(BondStereo::AtropCw, true);
        let before = both_flags.clone();
        assert_eq!(
            write_extensions_for_test(
                &both_flags,
                CxSmilesFields::BOND_CFG | CxSmilesFields::BOND_ATROPISOMER,
            )
            .expect("CX_BOND_CFG takes precedence over CX_BOND_ATROPISOMER"),
            "|wU:1.2|"
        );
        assert_eq!(
            both_flags, before,
            "extension writing must not mutate input"
        );
    }

    #[test]
    fn writes_data_polymer_and_hierarchy_records() {
        for (input, expected) in [
            ("OCC |SgD:0,2:foo:bar::::|", "CCO |SgD:2,0:foo:bar::::|"),
            ("OCC |Sg:n:0,1:lab:ht:::|", "CCO |Sg:n:2,1:lab:ht:::|"),
            (
                "OCC |SgD:0:a:b::::,Sg:n:1,2::ht:::,SgH:1:0|",
                "CCO |SgD:2:a:b::::,Sg:n:1,0::ht:::,SgH:1:0|",
            ),
        ] {
            assert_eq!(
                write_cx_smiles(&parse(input)).unwrap(),
                expected.into(),
                "{input}"
            );
        }
    }

    #[test]
    fn cx_sgroup_typed_crossings_follow_explicit_bond_traversal_order() {
        let record = parse("CCCCC |Sg:n:1,2,3:repeat:ht:0,0,3:3,3,0:|");
        let atom_order = (0..record.topology.atoms.len())
            .map(AtomId::new)
            .collect::<Vec<_>>();
        let bond_order = [
            BondId::new(3),
            BondId::new(2),
            BondId::new(1),
            BondId::new(0),
        ];
        assert_eq!(
            write_polymer_sgroups(&record, &atom_order, &bond_order)
                .expect("typed crossings map through traversal order"),
            "Sg:n:1,2,3:repeat:ht:3,3,0:0,0,3:".into()
        );
    }

    #[test]
    fn review_hierarchy_resolves_sparse_source_indices_and_selected_outputs() {
        use cosmolkit_model::SubstanceGroupId;
        let mut record = parse("C");
        let mut child = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data);
        child.set_prop("index", "42");
        child.set_parent(SubstanceGroupId::new(2));
        let mut sibling = SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data);
        sibling.set_prop("index", "7");
        // Preserve the existing raw-property path when no typed parent exists.
        sibling.set_prop("PARENT", "900");
        let mut parent = SubstanceGroup::new(
            SubstanceGroupId::new(2),
            SubstanceGroupKind::StructuralRepeatUnit,
        );
        parent.set_prop("index", "900");
        record.topology.substance_groups = vec![child, sibling, parent];
        assert_eq!(
            write_sgroup_hierarchy(&record, true, true).unwrap(),
            "SgH:2:0.1"
        );
        assert_eq!(write_sgroup_hierarchy(&record, true, false).unwrap(), "");
        assert_eq!(write_sgroup_hierarchy(&record, false, true).unwrap(), "");
        assert_eq!(write_sgroup_hierarchy(&record, false, false).unwrap(), "");
    }

    #[test]
    fn review_hierarchy_many_typed_children_keep_row_order() {
        use cosmolkit_model::SubstanceGroupId;
        let mut record = parse("C");
        let count = 512;
        let parent_id = SubstanceGroupId::new(count - 1);
        record.topology.substance_groups = (0..count)
            .map(|row| {
                let mut group =
                    SubstanceGroup::new(SubstanceGroupId::new(row), SubstanceGroupKind::Data);
                group.set_prop("index", (1000 + row * 7).to_string());
                if row + 1 != count {
                    group.set_parent(parent_id);
                }
                group
            })
            .collect();
        let children = (0..count - 1)
            .map(|row| row.to_string())
            .collect::<Vec<_>>()
            .join(".");
        assert_eq!(
            write_sgroup_hierarchy(&record, true, false).unwrap(),
            format!("SgH:{}:{children}", count - 1)
        );
    }

    #[test]
    fn writes_large_ring_cis_trans_and_unknown_blocks() {
        // The input remains raw CX state at this writer boundary. Pinned
        // RDKit 2026.03.1 with sanitize=false/removeHs=false, legacy stereo,
        // canonical/isomeric/cleanStereo=true, CX_ALL and default Clear runs
        // the post-base clean stage before extension emission: explicit ring
        // cis/trans annotations are removed when that stage resolves them.
        for (input, expected) in [
            ("C1CCCC/C=C/CCC1 |t:5|", "C1=C/CCCCCCCC/1"),
            ("C1CCCCC=CCCC1 |c:5|", "C1=CCCCCCCCC1"),
            ("C1=CCCCCCCCC1 |ctu:0|", "C1=CCCCCCCCC1 |ctu:0|"),
        ] {
            assert_eq!(
                write_cx_smiles(&parse(input)).unwrap(),
                expected.into(),
                "{input}"
            );
        }
    }

    #[test]
    fn writes_link_nodes_with_rdkit_canonical_index_behavior() {
        for (input, expected) in [
            ("OC1CCC(F)C1 |LN:1:1.3.2.6|", "OC1CCC(F)C1 |LN:1:1.3.2.6|"),
            (
                "FC1CCC(O)C1 |LN:1:1.3.2.6,4:1.4.3.6|",
                "OC1CCC(F)C1 |LN:4:1.3.3.6,1:1.4.2.6|",
            ),
            ("C1OCCC1C |LN:0:1.5,1:1.3|", "CC1CCOC1 |LN:5:1.5,4:1.3|"),
        ] {
            assert_eq!(
                write_cx_smiles(&parse(input)).unwrap(),
                expected.into(),
                "{input}"
            );
        }
    }

    #[test]
    fn field_selection_can_omit_coordinates_without_losing_other_records() {
        let record = parse("CC |(0,0,;1,0,),$_AV:left;right$|");
        let result = write_cx_smiles_with_params(
            &record,
            &CxSmilesWriteParams {
                fields: CxSmilesFields::ALL_BUT_COORDS,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(result, "CC |$_AV:left;right$|".into());
    }

    #[test]
    fn cx_post_base_cache_preparation_skips_when_clean_stereo_is_disabled() {
        let record = parse("CC");
        let before = record.clone();

        assert!(
            prepare_cx_post_base_valence(&record, false)
                .expect("cleanStereo=false skips the outer cache update")
                .is_none()
        );
        assert_eq!(record, before, "cache preparation must borrow its record");
    }

    #[test]
    fn cx_post_base_cache_preparation_runs_for_present_done_property() {
        // The pinned wrapper performs its outer needsUpdatePropertyCache check
        // before assignStereochemistry's force=false _StereochemDone presence
        // guard. The detached input has no cache fields; a zero-like marker,
        // whether registered computed or ordinary, cannot stand in for them.
        for computed in [false, true] {
            let mut record = parse("CC");
            if computed {
                record
                    .properties
                    .set_computed_prop("_StereochemDone", "0")
                    .expect("set computed done marker");
            } else {
                record
                    .properties
                    .set_prop("_StereochemDone", "0")
                    .expect("set ordinary done marker");
            }
            let before = record.clone();

            let assignment = prepare_cx_post_base_valence(&record, true)
                .expect("missing detached cache is recomputed non-strictly")
                .expect("cleanStereo=true prepares valence");
            assert_eq!(assignment.explicit_valence, [1, 1]);
            assert_eq!(assignment.implicit_hydrogens, [3, 3]);
            assert_eq!(record, before, "cache preparation must not mutate input");
        }
    }

    #[test]
    fn cx_post_base_cache_preparation_uses_update_property_cache_false_value() {
        // Pinned Atom::updatePropertyCache(false) calls
        // calculateExplicitValence(false), whose checkIt=false path skips the
        // `strict || checkIt` validation gate and keeps explicit valence 5.
        let parser = SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        };
        let record = parse_smiles("C(F)(F)(F)(F)F", &parser).expect("parse raw pentavalent C");

        let assignment = prepare_cx_post_base_valence(&record, true)
            .expect("non-strict cache preparation keeps calculated overvalence")
            .expect("cleanStereo=true prepares valence");
        assert_eq!(assignment.explicit_valence[0], 5);
    }

    #[test]
    fn cx_post_base_assignment_clears_pending_and_sets_computed_done_marker() {
        let parser = SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        };
        let input = parse_smiles("C1CCCCC=CCCC1 |c:5|", &parser).expect("parse raw ring CXSMILES");
        assert_eq!(
            input.properties.prop("_needsDetectBondStereo"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
        assert_eq!(input.properties.prop("_StereochemDone"), None);
        let input_before = input.clone();
        let mut prepared_record = input.clone();

        let valence = prepare_cx_post_base_valence(&prepared_record, true)
            .expect("outer cache preparation succeeds")
            .expect("cleanStereo=true prepares valence");
        apply_cx_post_base_stereochemistry(&mut prepared_record, &valence)
            .expect("absent done marker runs legacy assignment");

        assert_eq!(
            prepared_record.properties.prop("_needsDetectBondStereo"),
            None
        );
        assert_eq!(
            prepared_record.properties.prop("_StereochemDone"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
        assert!(
            prepared_record
                .properties
                .is_prop_computed("_StereochemDone")
                .unwrap()
        );
        assert!(
            !prepared_record
                .properties
                .is_prop_computed("_needsDetectBondStereo")
                .unwrap()
        );
        assert_eq!(
            input, input_before,
            "post-base work must preserve the caller"
        );
    }

    #[test]
    fn cx_post_base_zero_done_marker_skips_assignment_but_cleans_groups() {
        let mut input = parse("F[C@](Cl)(Br)I.F[C@](Cl)(Br)I");
        input.topology.atoms[6].set_chiral_tag(cosmolkit_types::ChiralTag::Unspecified);
        input.topology.stereo_groups = vec![
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(1), AtomId::new(6)],
                Vec::new(),
            )
            .expect("valid distinct stereo members")
            .with_id(7),
        ];
        input
            .properties
            .set_prop("_needsDetectBondStereo", "1")
            .expect("set pending source marker");
        input
            .properties
            .set_prop("_StereochemDone", "0")
            .expect("set present zero-like source marker");
        let input_before = input.clone();
        let original_first_tag = input.topology.atoms[1].chiral_tag();
        let mut prepared_record = input.clone();

        let valence = prepare_cx_post_base_valence(&prepared_record, true)
            .expect("outer cache preparation precedes the guard")
            .expect("cleanStereo=true prepares valence");
        apply_cx_post_base_stereochemistry(&mut prepared_record, &valence)
            .expect("zero-like present marker skips assignment and reaches cleanup");

        assert_eq!(
            prepared_record.properties.prop("_StereochemDone"),
            Some(&cosmolkit_model::PropertyValue::String("0".into()))
        );
        assert!(
            !prepared_record
                .properties
                .is_prop_computed("_StereochemDone")
                .unwrap()
        );
        assert_eq!(
            prepared_record.properties.prop("_needsDetectBondStereo"),
            Some(&cosmolkit_model::PropertyValue::String("1".into())),
            "the legacy perception clear occurs only when assignment runs"
        );
        assert_eq!(
            prepared_record.topology.atoms[1].chiral_tag(),
            original_first_tag,
            "the force=false presence guard skips atom assignment"
        );
        assert_eq!(
            prepared_record.topology.atoms[6].chiral_tag(),
            cosmolkit_types::ChiralTag::Unspecified
        );
        assert_eq!(prepared_record.topology.stereo_groups.len(), 1);
        let group = &prepared_record.topology.stereo_groups[0];
        assert_eq!(group.kind(), StereoGroupKind::And);
        assert_eq!(group.id(), Some(7));
        assert_eq!(group.atoms(), &[AtomId::new(1)]);
        assert!(group.bonds().is_empty());
        assert_eq!(
            input, input_before,
            "post-base work must preserve the caller"
        );
    }

    #[test]
    fn cx_post_base_cache_preparation_propagates_owner_errors() {
        // This malformed detached topology is not a valid source ROMol; this
        // owner-boundary test only verifies that the existing typed valence
        // validation error is not replaced with a guessed assignment.
        let mut record = parse("CC");
        record.topology.atoms.swap(0, 1);
        let before = record.clone();

        let error = prepare_cx_post_base_valence(&record, true)
            .expect_err("invalid detached row identities must propagate");
        assert!(matches!(error, SmilesParseError::WriterValence(_)));
        assert_eq!(
            record, before,
            "failed cache preparation must not mutate input"
        );
    }

    #[test]
    fn cx_post_base_empty_output_skips_the_post_base_stage() {
        let record = SmilesRecord {
            topology: cosmolkit_model::TopologyBlock::default(),
            coordinates: cosmolkit_model::CoordinateBlock::default(),
            properties: cosmolkit_model::MoleculeProperties::default(),
        };
        let before = record.clone();
        let params = CxSmilesWriteParams {
            smiles: SmilesWriteParams {
                clean_stereo: true,
                ..Default::default()
            },
            fields: CxSmilesFields::ALL,
            ..Default::default()
        };

        assert_eq!(
            write_cx_smiles_with_params(&record, &params).unwrap(),
            "".into()
        );
        assert_eq!(record, before, "empty CX writing must preserve its caller");
    }

    #[test]
    fn coordinate_numbers_use_c_percent_g_shape() {
        for (value, expected) in [
            (0.0, "0"),
            (0.0001, "0.0001"),
            (0.00123456789, "0.00123457"),
            (1.23456789, "1.23457"),
            (12_345.6789, "12345.7"),
            (123_456.789, "123457"),
            (1_234_567.89, "1.23457e+06"),
            (-1_234_567.89, "-1.23457e+06"),
            (999_999.9, "1e+06"),
        ] {
            assert_eq!(format_general(value), expected, "{value:?}");
        }
    }
}

#[cfg(test)]
mod uint_cx_proposed_tests {
    use super::*;
    #[test]
    fn proposed_uint_cx_unsigned_and_string_full_width() {
        for (value, text) in [
            (0_u32, "0"),
            (1, "1"),
            (2147483646, "2147483646"),
            (2147483647, "2147483647"),
            (2147483648, "2147483648"),
            (4294967295, "4294967295"),
        ] {
            assert_eq!(
                source_unsigned_property(&PropertyValue::UInt(value), "_MolFileBondCfg"),
                Ok(value)
            );
            assert_eq!(
                source_string_property(Some(&PropertyValue::UInt(value)), "atomLabel")
                    .unwrap()
                    .expect("original present property cell")
                    .as_bytes(),
                text.as_bytes()
            );
        }
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    use cosmolkit_model::PropertyValue;
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

    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXwriter_0
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxwriter_0_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(0_u32)]);
        let value = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_unsigned_property(value, "_cxsmilesBondIdx"),
            Ok(0_u32)
        );
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_smiles/CXtext_0
    #[test]
    fn uint_cell_text_consumer_smiles_cxtext_0_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(0_u32)]);
        let v = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_string_property(Some(v), "_cxsmilesBondIdx")
                .unwrap()
                .expect("original present property cell")
                .as_bytes(),
            b"0".as_slice()
        );
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXwriter_1
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxwriter_1_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(1_u32)]);
        let value = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_unsigned_property(value, "_cxsmilesBondIdx"),
            Ok(1_u32)
        );
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_smiles/CXtext_1
    #[test]
    fn uint_cell_text_consumer_smiles_cxtext_1_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(1_u32)]);
        let v = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_string_property(Some(v), "_cxsmilesBondIdx")
                .unwrap()
                .expect("original present property cell")
                .as_bytes(),
            b"1".as_slice()
        );
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXwriter_2147483646
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxwriter_2147483646_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(2147483646_u32)]);
        let value = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_unsigned_property(value, "_cxsmilesBondIdx"),
            Ok(2147483646_u32)
        );
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_smiles/CXtext_2147483646
    #[test]
    fn uint_cell_text_consumer_smiles_cxtext_2147483646_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(2147483646_u32)]);
        let v = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_string_property(Some(v), "_cxsmilesBondIdx")
                .unwrap()
                .expect("original present property cell")
                .as_bytes(),
            b"2147483646".as_slice()
        );
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXwriter_2147483647
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxwriter_2147483647_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(2147483647_u32)]);
        let value = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_unsigned_property(value, "_cxsmilesBondIdx"),
            Ok(2147483647_u32)
        );
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_smiles/CXtext_2147483647
    #[test]
    fn uint_cell_text_consumer_smiles_cxtext_2147483647_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(2147483647_u32)]);
        let v = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_string_property(Some(v), "_cxsmilesBondIdx")
                .unwrap()
                .expect("original present property cell")
                .as_bytes(),
            b"2147483647".as_slice()
        );
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXwriter_2147483648
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxwriter_2147483648_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(2147483648_u32)]);
        let value = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_unsigned_property(value, "_cxsmilesBondIdx"),
            Ok(2147483648_u32)
        );
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_smiles/CXtext_2147483648
    #[test]
    fn uint_cell_text_consumer_smiles_cxtext_2147483648_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(2147483648_u32)]);
        let v = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_string_property(Some(v), "_cxsmilesBondIdx")
                .unwrap()
                .expect("original present property cell")
                .as_bytes(),
            b"2147483648".as_slice()
        );
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_smiles/CXwriter_4294967295
    #[test]
    fn uint_cell_unsigned_consumer_smiles_cxwriter_4294967295_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(4294967295_u32)]);
        let value = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_unsigned_property(value, "_cxsmilesBondIdx"),
            Ok(4294967295_u32)
        );
    }
    // FROZEN UINT CONDITION: TEXT_CONSUMER_smiles/CXtext_4294967295
    #[test]
    fn uint_cell_text_consumer_smiles_cxtext_4294967295_cx_writer() {
        let g = graph(vec![PropertyValue::UInt(4294967295_u32)]);
        let v = g.bonds[0].prop("_cxsmilesBondIdx").unwrap();
        assert_eq!(
            source_string_property(Some(v), "_cxsmilesBondIdx")
                .unwrap()
                .expect("original present property cell")
                .as_bytes(),
            b"4294967295".as_slice()
        );
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_0
    #[test]
    fn uint_cell_unsigned_cfg_0_cx_writer() {
        let mut g = graph(vec![PropertyValue::UInt(0)]);
        g.bonds[0]
            .set_prop("_MolFileBondCfg", PropertyValue::UInt(0_u32))
            .unwrap();
        let before = g.clone();
        let cfg = source_unsigned_property(
            g.bonds[0].prop("_MolFileBondCfg").unwrap(),
            "_MolFileBondCfg",
        )
        .unwrap();
        assert_eq!(cfg, 0_u32);
        assert_eq!(molfile_cfg_bond_direction(Some(cfg)), BondDirection::None);
        assert_eq!(g, before);
        let record = crate::SmilesRecord {
            topology: g.clone(),
            coordinates: Default::default(),
            properties: Default::default(),
        };
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .unwrap();
        let rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let context = CrossedBondContext::new(&record.topology, &valence, &rings, true).unwrap();
        assert_eq!(
            write_bond_config(
                &record,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                false,
                &WedgeAssignments::default(),
                Some(&context),
                None
            ),
            Ok("".into())
        );
        assert_eq!(record.topology, g);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_1
    #[test]
    fn uint_cell_unsigned_cfg_1_cx_writer() {
        let mut g = graph(vec![PropertyValue::UInt(0)]);
        g.bonds[0]
            .set_prop("_MolFileBondCfg", PropertyValue::UInt(1_u32))
            .unwrap();
        let before = g.clone();
        let cfg = source_unsigned_property(
            g.bonds[0].prop("_MolFileBondCfg").unwrap(),
            "_MolFileBondCfg",
        )
        .unwrap();
        assert_eq!(cfg, 1_u32);
        assert_eq!(
            molfile_cfg_bond_direction(Some(cfg)),
            BondDirection::BeginWedge
        );
        assert_eq!(g, before);
        let record = crate::SmilesRecord {
            topology: g.clone(),
            coordinates: Default::default(),
            properties: Default::default(),
        };
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .unwrap();
        let rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let context = CrossedBondContext::new(&record.topology, &valence, &rings, true).unwrap();
        assert_eq!(
            write_bond_config(
                &record,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                false,
                &WedgeAssignments::default(),
                Some(&context),
                None
            ),
            Ok("wU:0.0".into())
        );
        assert_eq!(record.topology, g);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_2
    #[test]
    fn uint_cell_unsigned_cfg_2_cx_writer() {
        let mut g = graph(vec![PropertyValue::UInt(0)]);
        g.bonds[0]
            .set_prop("_MolFileBondCfg", PropertyValue::UInt(2_u32))
            .unwrap();
        let before = g.clone();
        let cfg = source_unsigned_property(
            g.bonds[0].prop("_MolFileBondCfg").unwrap(),
            "_MolFileBondCfg",
        )
        .unwrap();
        assert_eq!(cfg, 2_u32);
        assert_eq!(
            molfile_cfg_bond_direction(Some(cfg)),
            BondDirection::Unknown
        );
        assert_eq!(g, before);
        let record = crate::SmilesRecord {
            topology: g.clone(),
            coordinates: Default::default(),
            properties: Default::default(),
        };
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .unwrap();
        let rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let context = CrossedBondContext::new(&record.topology, &valence, &rings, true).unwrap();
        assert_eq!(
            write_bond_config(
                &record,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                false,
                &WedgeAssignments::default(),
                Some(&context),
                None
            ),
            Ok("w:0.0".into())
        );
        assert_eq!(record.topology, g);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_3
    #[test]
    fn uint_cell_unsigned_cfg_3_cx_writer() {
        let mut g = graph(vec![PropertyValue::UInt(0)]);
        g.bonds[0]
            .set_prop("_MolFileBondCfg", PropertyValue::UInt(3_u32))
            .unwrap();
        let before = g.clone();
        let cfg = source_unsigned_property(
            g.bonds[0].prop("_MolFileBondCfg").unwrap(),
            "_MolFileBondCfg",
        )
        .unwrap();
        assert_eq!(cfg, 3_u32);
        assert_eq!(
            molfile_cfg_bond_direction(Some(cfg)),
            BondDirection::BeginDash
        );
        assert_eq!(g, before);
        let record = crate::SmilesRecord {
            topology: g.clone(),
            coordinates: Default::default(),
            properties: Default::default(),
        };
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .unwrap();
        let rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let context = CrossedBondContext::new(&record.topology, &valence, &rings, true).unwrap();
        assert_eq!(
            write_bond_config(
                &record,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                false,
                &WedgeAssignments::default(),
                Some(&context),
                None
            ),
            Ok("wD:0.0".into())
        );
        assert_eq!(record.topology, g);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_4
    #[test]
    fn uint_cell_unsigned_cfg_4_cx_writer() {
        let mut g = graph(vec![PropertyValue::UInt(0)]);
        g.bonds[0]
            .set_prop("_MolFileBondCfg", PropertyValue::UInt(4_u32))
            .unwrap();
        let before = g.clone();
        let cfg = source_unsigned_property(
            g.bonds[0].prop("_MolFileBondCfg").unwrap(),
            "_MolFileBondCfg",
        )
        .unwrap();
        assert_eq!(cfg, 4_u32);
        assert_eq!(molfile_cfg_bond_direction(Some(cfg)), BondDirection::None);
        assert_eq!(g, before);
        let record = crate::SmilesRecord {
            topology: g.clone(),
            coordinates: Default::default(),
            properties: Default::default(),
        };
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .unwrap();
        let rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let context = CrossedBondContext::new(&record.topology, &valence, &rings, true).unwrap();
        assert_eq!(
            write_bond_config(
                &record,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                false,
                &WedgeAssignments::default(),
                Some(&context),
                None
            ),
            Ok("".into())
        );
        assert_eq!(record.topology, g);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_255
    #[test]
    fn uint_cell_unsigned_cfg_255_cx_writer() {
        let mut g = graph(vec![PropertyValue::UInt(0)]);
        g.bonds[0]
            .set_prop("_MolFileBondCfg", PropertyValue::UInt(255_u32))
            .unwrap();
        let before = g.clone();
        let cfg = source_unsigned_property(
            g.bonds[0].prop("_MolFileBondCfg").unwrap(),
            "_MolFileBondCfg",
        )
        .unwrap();
        assert_eq!(cfg, 255_u32);
        assert_eq!(molfile_cfg_bond_direction(Some(cfg)), BondDirection::None);
        assert_eq!(g, before);
        let record = crate::SmilesRecord {
            topology: g.clone(),
            coordinates: Default::default(),
            properties: Default::default(),
        };
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .unwrap();
        let rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let context = CrossedBondContext::new(&record.topology, &valence, &rings, true).unwrap();
        assert_eq!(
            write_bond_config(
                &record,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                false,
                &WedgeAssignments::default(),
                Some(&context),
                None
            ),
            Ok("".into())
        );
        assert_eq!(record.topology, g);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_256
    #[test]
    fn uint_cell_unsigned_cfg_256_cx_writer() {
        let mut g = graph(vec![PropertyValue::UInt(0)]);
        g.bonds[0]
            .set_prop("_MolFileBondCfg", PropertyValue::UInt(256_u32))
            .unwrap();
        let before = g.clone();
        let cfg = source_unsigned_property(
            g.bonds[0].prop("_MolFileBondCfg").unwrap(),
            "_MolFileBondCfg",
        )
        .unwrap();
        assert_eq!(cfg, 256_u32);
        assert_eq!(molfile_cfg_bond_direction(Some(cfg)), BondDirection::None);
        assert_eq!(g, before);
        let record = crate::SmilesRecord {
            topology: g.clone(),
            coordinates: Default::default(),
            properties: Default::default(),
        };
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .unwrap();
        let rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let context = CrossedBondContext::new(&record.topology, &valence, &rings, true).unwrap();
        assert_eq!(
            write_bond_config(
                &record,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                false,
                &WedgeAssignments::default(),
                Some(&context),
                None
            ),
            Ok("".into())
        );
        assert_eq!(record.topology, g);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_2147483647
    #[test]
    fn uint_cell_unsigned_cfg_2147483647_cx_writer() {
        let mut g = graph(vec![PropertyValue::UInt(0)]);
        g.bonds[0]
            .set_prop("_MolFileBondCfg", PropertyValue::UInt(2147483647_u32))
            .unwrap();
        let before = g.clone();
        let cfg = source_unsigned_property(
            g.bonds[0].prop("_MolFileBondCfg").unwrap(),
            "_MolFileBondCfg",
        )
        .unwrap();
        assert_eq!(cfg, 2147483647_u32);
        assert_eq!(molfile_cfg_bond_direction(Some(cfg)), BondDirection::None);
        assert_eq!(g, before);
        let record = crate::SmilesRecord {
            topology: g.clone(),
            coordinates: Default::default(),
            properties: Default::default(),
        };
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .unwrap();
        let rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let context = CrossedBondContext::new(&record.topology, &valence, &rings, true).unwrap();
        assert_eq!(
            write_bond_config(
                &record,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                false,
                &WedgeAssignments::default(),
                Some(&context),
                None
            ),
            Ok("".into())
        );
        assert_eq!(record.topology, g);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_2147483648
    #[test]
    fn uint_cell_unsigned_cfg_2147483648_cx_writer() {
        let mut g = graph(vec![PropertyValue::UInt(0)]);
        g.bonds[0]
            .set_prop("_MolFileBondCfg", PropertyValue::UInt(2147483648_u32))
            .unwrap();
        let before = g.clone();
        let cfg = source_unsigned_property(
            g.bonds[0].prop("_MolFileBondCfg").unwrap(),
            "_MolFileBondCfg",
        )
        .unwrap();
        assert_eq!(cfg, 2147483648_u32);
        assert_eq!(molfile_cfg_bond_direction(Some(cfg)), BondDirection::None);
        assert_eq!(g, before);
        let record = crate::SmilesRecord {
            topology: g.clone(),
            coordinates: Default::default(),
            properties: Default::default(),
        };
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .unwrap();
        let rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let context = CrossedBondContext::new(&record.topology, &valence, &rings, true).unwrap();
        assert_eq!(
            write_bond_config(
                &record,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                false,
                &WedgeAssignments::default(),
                Some(&context),
                None
            ),
            Ok("".into())
        );
        assert_eq!(record.topology, g);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_4294967295
    #[test]
    fn uint_cell_unsigned_cfg_4294967295_cx_writer() {
        let mut g = graph(vec![PropertyValue::UInt(0)]);
        g.bonds[0]
            .set_prop("_MolFileBondCfg", PropertyValue::UInt(4294967295_u32))
            .unwrap();
        let before = g.clone();
        let cfg = source_unsigned_property(
            g.bonds[0].prop("_MolFileBondCfg").unwrap(),
            "_MolFileBondCfg",
        )
        .unwrap();
        assert_eq!(cfg, 4294967295_u32);
        assert_eq!(molfile_cfg_bond_direction(Some(cfg)), BondDirection::None);
        assert_eq!(g, before);
        let record = crate::SmilesRecord {
            topology: g.clone(),
            coordinates: Default::default(),
            properties: Default::default(),
        };
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .unwrap();
        let rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let context = CrossedBondContext::new(&record.topology, &valence, &rings, true).unwrap();
        assert_eq!(
            write_bond_config(
                &record,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                false,
                &WedgeAssignments::default(),
                Some(&context),
                None
            ),
            Ok("".into())
        );
        assert_eq!(record.topology, g);
    }
}

#[cfg(test)]
fn fixed_property_text(value: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(value.as_bytes()).expect("original fixed fixture text is UTF8")
}

#[cfg(test)]
mod source_getter_tests {
    use super::*;
    use cosmolkit_model::{CoordinateBlock, CoordinateDimension, CoordinateValidationError};

    #[test]
    fn source_getter_handles_negative_ids_ordered_duplicate_ids_and_errors() {
        for first in [CoordinateDimension::TwoD, CoordinateDimension::ThreeD] {
            let order = if first == CoordinateDimension::TwoD {
                vec![
                    CoordinateDimension::TwoD,
                    CoordinateDimension::ThreeD,
                    CoordinateDimension::ThreeD,
                ]
            } else {
                vec![
                    CoordinateDimension::ThreeD,
                    CoordinateDimension::TwoD,
                    CoordinateDimension::ThreeD,
                ]
            };
            let rows = CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(19, vec![[1.0, 2.0]])],
                conformers_3d: vec![
                    Conformer3D::new(19, vec![[3.0, 4.0, 5.0]], false),
                    Conformer3D::new(0, vec![[6.0, 7.0, 8.0]], true),
                ],
                source_conformer_order: Some(order),
                ..Default::default()
            };
            let before = rows.clone();
            for id in [-1, i32::MIN, 19] {
                match (first, source_conformer_by_id(&rows, id).unwrap()) {
                    (CoordinateDimension::TwoD, CoordinateSource::TwoD(row)) => {
                        assert!(std::ptr::eq(row, &rows.conformers_2d[0]))
                    }
                    (CoordinateDimension::ThreeD, CoordinateSource::ThreeD(row)) => {
                        assert!(std::ptr::eq(row, &rows.conformers_3d[0]))
                    }
                    _ => panic!("not the first inserted matching conformer"),
                }
            }
            match source_conformer_by_id(&rows, 0).unwrap() {
                CoordinateSource::ThreeD(row) => assert!(std::ptr::eq(row, &rows.conformers_3d[1])),
                _ => panic!("wrong explicit source ID"),
            }
            for id in [1, i32::MAX] {
                assert!(
                    matches!(source_conformer_by_id(&rows, id), Err(SmilesParseError::Model(message)) if message == format!("Can't find conformation with ID: {id}"))
                );
            }
            let mut missing_order = rows.clone();
            missing_order.source_conformer_order = None;
            for id in [-1, 19, 0] {
                assert!(matches!(
                    source_conformer_by_id(&missing_order, id),
                    Err(SmilesParseError::Coordinates(
                        CoordinateValidationError::MissingSourceConformerOrder
                    ))
                ));
            }
            assert_eq!(rows, before);
        }
        let empty = CoordinateBlock::default();
        for id in [-1, i32::MIN, 0, i32::MAX] {
            assert!(
                matches!(source_conformer_by_id(&empty, id), Err(SmilesParseError::Model(message)) if message == "No conformations available on the molecule")
            );
        }
        for rows in [
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(19, vec![[1.0, 2.0]])],
                ..Default::default()
            },
            CoordinateBlock {
                conformers_3d: vec![Conformer3D::new(19, vec![[3.0, 4.0, 5.0]], false)],
                ..Default::default()
            },
        ] {
            assert!(source_conformer_by_id(&rows, 19).is_ok());
            assert!(matches!(
                source_conformer_by_id(&rows, 0),
                Err(SmilesParseError::Model(_))
            ));
        }
    }
}

#[cfg(test)]
mod source_mutable_getter_tests {
    use super::*;
    use cosmolkit_model::{CoordinateBlock, CoordinateDimension};

    #[test]
    fn source_mutable_getter_aliases_const_first_matching_row_without_other_changes() {
        for first in [CoordinateDimension::TwoD, CoordinateDimension::ThreeD] {
            for id in [-1, i32::MIN, 19, 0] {
                let order = if first == CoordinateDimension::TwoD {
                    vec![
                        CoordinateDimension::TwoD,
                        CoordinateDimension::ThreeD,
                        CoordinateDimension::ThreeD,
                    ]
                } else {
                    vec![
                        CoordinateDimension::ThreeD,
                        CoordinateDimension::TwoD,
                        CoordinateDimension::ThreeD,
                    ]
                };
                let mut rows = CoordinateBlock {
                    conformers_2d: vec![Conformer2D::new(19, vec![[1.0, 2.0]])],
                    conformers_3d: vec![
                        Conformer3D::new(19, vec![[3.0, 4.0, 5.0]], false),
                        Conformer3D::new(0, vec![[6.0, 7.0, 8.0]], true),
                    ],
                    source_conformer_order: Some(order),
                    ..Default::default()
                };
                let before = rows.clone();
                let address = match source_conformer_by_id(&rows, id).unwrap() {
                    CoordinateSource::TwoD(row) => row as *const Conformer2D as usize,
                    CoordinateSource::ThreeD(row) => row as *const Conformer3D as usize,
                };
                match source_conformer_by_id_mut(&mut rows, id).unwrap() {
                    CoordinateSourceMut::TwoD(row) => {
                        assert_eq!(row as *const Conformer2D as usize, address);
                        row.coordinates_mut()[0][0] = 42.0;
                    }
                    CoordinateSourceMut::ThreeD(row) => {
                        assert_eq!(row as *const Conformer3D as usize, address);
                        row.coordinates_mut()[0][0] = 42.0;
                    }
                }
                let mut expected = before;
                if id == 0 {
                    expected.conformers_3d[1].coordinates_mut()[0][0] = 42.0;
                } else if first == CoordinateDimension::TwoD {
                    expected.conformers_2d[0].coordinates_mut()[0][0] = 42.0;
                } else {
                    expected.conformers_3d[0].coordinates_mut()[0][0] = 42.0;
                }
                assert_eq!(rows, expected);
                let before_failure = rows.clone();
                assert!(matches!(
                    source_conformer_by_id_mut(&mut rows, 1),
                    Err(SmilesParseError::Model(_))
                ));
                assert_eq!(rows, before_failure);
            }
        }
        let mut empty = CoordinateBlock::default();
        assert!(
            matches!(source_conformer_by_id_mut(&mut empty, -1), Err(SmilesParseError::Model(message)) if message == "No conformations available on the molecule")
        );
        assert_eq!(empty, CoordinateBlock::default());
    }
}

#[cfg(test)]
mod complete_link_node_source_tests {
    use super::*;
    use crate::{SmilesParseParams, parse_smiles};

    fn record(bytes: &[u8]) -> SmilesRecord {
        let mut record = parse_smiles("CC(C)C", &SmilesParseParams::default()).unwrap();
        record
            .properties
            .set_prop(
                "_molLinkNodes",
                PropertyValue::String(PropertyText::from_bytes(bytes)),
            )
            .unwrap();
        record
    }
    fn read(
        record: &SmilesRecord,
        strict: bool,
        map: Option<&BTreeMap<u32, AtomId>>,
        warnings: &mut Vec<SourceLinkIssue>,
    ) -> Result<Vec<SourceLinkNode>, SmilesParseError> {
        get_mol_link_nodes_source(&record.topology, &record.properties, strict, map, warnings)
    }

    #[test]
    fn absent_empty_and_numeric_property_use_the_real_string_owner() {
        let mut r = parse_smiles("CC", &SmilesParseParams::default()).unwrap();
        let mut warnings = Vec::new();
        assert!(read(&r, true, None, &mut warnings).unwrap().is_empty());
        r.properties.set_prop("_molLinkNodes", "").unwrap();
        assert!(read(&r, true, None, &mut warnings).unwrap().is_empty());
        r.properties
            .set_prop("_molLinkNodes", PropertyValue::UInt(123))
            .unwrap();
        assert!(read(&r, false, None, &mut warnings).unwrap().is_empty());
        assert_eq!(warnings[0].kind, SourceLinkIssueKind::TooFewValues);
        assert_eq!(warnings[0].raw_node.as_bytes(), b"123");
    }

    #[test]
    fn every_native_rejection_has_strict_error_and_nonstrict_outer_continue() {
        for (bad, kind, expected_error) in [
            (
                "x 3 2 2 1 2 3",
                SourceLinkIssueKind::UnsignedConversion,
                "Cannot convert values",
            ),
            (
                "1 3 2 2 1",
                SourceLinkIssueKind::TooFewValues,
                "not enough values",
            ),
            (
                "0 3 2 2 1 2 3",
                SourceLinkIssueKind::BadCounts,
                "bad counts",
            ),
            (
                "4 3 2 2 1 2 3",
                SourceLinkIssueKind::BadCounts,
                "bad counts",
            ),
            (
                "1 3 1 2 1",
                SourceLinkIssueKind::UnderConstruction,
                "Incomplete Code",
            ),
            (
                "1 3 2 2 1 1 3",
                SourceLinkIssueKind::DifferentCenters,
                "bonds don't start",
            ),
            (
                "1 3 2 2 2 2 3",
                SourceLinkIssueKind::MissingBond,
                "bond not found",
            ),
        ] {
            let r = record(format!("{bad}|1 3 2 2 1 2 3").as_bytes());
            let mut warnings = Vec::new();
            assert!(
                read(&r, true, None, &mut warnings)
                    .unwrap_err()
                    .to_string()
                    .contains(expected_error),
                "{bad}"
            );
            assert!(warnings.is_empty());
            let nodes = read(&r, false, None, &mut warnings).unwrap();
            assert_eq!(nodes.len(), 1, "{bad}");
            assert_eq!(warnings.len(), 1);
            assert_eq!(warnings[0].kind, kind);
            assert_eq!(warnings[0].raw_node.as_bytes(), bad.as_bytes());
            assert_eq!(nodes[0].bond_atoms, [(1, 0), (1, 2)]);
        }
    }

    #[test]
    fn separators_drop_only_explicit_empty_tokens_and_keep_raw_bad_bytes() {
        let r = record(b"||  +1   3  2  2 1 2 3  ||1 3 2 2 1 2 4||");
        let mut warnings = Vec::new();
        let nodes = read(&r, true, None, &mut warnings).unwrap();
        assert_eq!(nodes.len(), 2);
        assert_eq!(nodes[1].bond_atoms[1], (1, 3));
        assert!(warnings.is_empty());
        for suffix in [b'\t', b'\n', b'\r', 0x0b, 0x0c, 0, 0xff] {
            let mut raw = b"1 3 2 2 1 2 3".to_vec();
            raw.push(suffix);
            let r = record(&raw);
            warnings.clear();
            assert!(read(&r, false, None, &mut warnings).unwrap().is_empty());
            assert_eq!(warnings[0].kind, SourceLinkIssueKind::UnsignedConversion);
            assert_eq!(warnings[0].raw_node.as_bytes(), raw);
        }
    }

    #[test]
    fn uint_sign_boundaries_and_required_size_wrap_precede_bond_count_guard() {
        let mut warnings = Vec::new();
        let r = record(b"-1 -1 2 2 1 2 3");
        let nodes = read(&r, true, None, &mut warnings).unwrap();
        assert_eq!(nodes[0].min_repetitions, u32::MAX);
        assert_eq!(nodes[0].max_repetitions, u32::MAX);
        for (raw, kind) in [
            (
                b"1 3 2147483649 2 1".as_slice(),
                SourceLinkIssueKind::UnderConstruction,
            ),
            (
                b"1 3 2147483647 2 1".as_slice(),
                SourceLinkIssueKind::UnderConstruction,
            ),
            (
                b"1 3 2147483646 2 1".as_slice(),
                SourceLinkIssueKind::TooFewValues,
            ),
            (
                b"1 4294967296 2 2 1 2 3".as_slice(),
                SourceLinkIssueKind::UnsignedConversion,
            ),
        ] {
            warnings.clear();
            assert!(
                read(&record(raw), false, None, &mut warnings)
                    .unwrap()
                    .is_empty()
            );
            assert_eq!(warnings[0].kind, kind);
        }
    }

    #[test]
    fn supplied_map_uses_actual_atom_ids_and_checks_centers_before_mapping() {
        let r = record(b"1 3 2 11 21 11 31");
        let map = BTreeMap::from([
            (10, AtomId::new(1)),
            (20, AtomId::new(0)),
            (30, AtomId::new(3)),
        ]);
        let mut warnings = Vec::new();
        assert_eq!(
            read(&r, true, Some(&map), &mut warnings).unwrap()[0].bond_atoms,
            [(1, 0), (1, 3)]
        );
        let r = record(b"1 3 2 11 21 12 31");
        assert!(
            read(&r, false, Some(&map), &mut warnings)
                .unwrap()
                .is_empty()
        );
        assert_eq!(warnings.len(), 1);
        assert_eq!(warnings[0].kind, SourceLinkIssueKind::DifferentCenters);
    }

    #[test]
    fn missing_map_entries_continue_inner_loop_and_can_still_accept_node() {
        let r = record(b"1 3 2 1 21 1 31");
        let map = BTreeMap::from([(20, AtomId::new(0)), (30, AtomId::new(3))]);
        let mut warnings = Vec::new();
        assert_eq!(
            read(&r, false, Some(&map), &mut warnings).unwrap()[0].bond_atoms,
            [(1, 0), (1, 3)]
        );
        assert_eq!(
            warnings.iter().map(|w| &w.kind).collect::<Vec<_>>(),
            vec![
                &SourceLinkIssueKind::MissingMappedAtom(1),
                &SourceLinkIssueKind::MissingMappedAtom(1)
            ]
        );
        warnings.clear();
        assert!(
            read(&r, true, Some(&map), &mut warnings)
                .unwrap_err()
                .to_string()
                .contains("atom index 1 cannot be found")
        );
        assert!(warnings.is_empty());
        let r = record(b"1 3 2 0 21 0 31");
        let map = BTreeMap::from([
            (u32::MAX, AtomId::new(1)),
            (20, AtomId::new(0)),
            (30, AtomId::new(3)),
        ]);
        assert_eq!(
            read(&r, true, Some(&map), &mut warnings).unwrap()[0].bond_atoms[0],
            (1, 0)
        );
    }

    #[test]
    fn bond_lookup_range_errors_ignore_strict_and_first_missing_bond_short_circuits() {
        for strict in [false, true] {
            let mut warnings = Vec::new();
            assert!(
                read(&record(b"1 3 2 2 0 2 3"), strict, None, &mut warnings)
                    .unwrap_err()
                    .to_string()
                    .contains("4294967295")
            );
            assert!(warnings.is_empty());
            assert!(
                read(&record(b"1 3 2 0 1 0 3"), strict, None, &mut warnings)
                    .unwrap_err()
                    .to_string()
                    .contains("4294967295")
            );
            assert!(warnings.is_empty());
        }
        let mut warnings = Vec::new();
        assert!(
            read(&record(b"1 3 2 2 2 2 0"), false, None, &mut warnings)
                .unwrap()
                .is_empty()
        );
        assert_eq!(warnings[0].kind, SourceLinkIssueKind::MissingBond);
        warnings.clear();
        assert!(read(&record(b"1 3 2 2 1 2 0"), false, None, &mut warnings).is_err());
        assert!(warnings.is_empty());
    }

    #[test]
    fn strict_late_failure_and_writer_preserve_source_order_and_atomorder_indexing() {
        let r = record(b"1 3 2 2 1 2 3|bad|2 4 2 2 3 2 4");
        let mut warnings = Vec::new();
        assert!(read(&r, true, None, &mut warnings).is_err());
        assert!(warnings.is_empty());
        let nodes = read(&r, false, None, &mut warnings).unwrap();
        assert_eq!(
            nodes.iter().map(|n| n.min_repetitions).collect::<Vec<_>>(),
            vec![1, 2]
        );
        let order = [
            AtomId::new(3),
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
        ];
        assert_eq!(
            write_link_nodes(&r, &order).unwrap(),
            "LN:0:1.3.3.1,0:2.4.1.2"
        );
    }
}

#[cfg(test)]
mod cx_coordinate_shared_source_tests {
    use super::*;
    #[test]
    fn nonfinite_coordinate_strings_reuse_source_case_and_sign() {
        for (value, expected) in [
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            (f64::from_bits(0x7ff8000000001234), "nan"),
            (f64::from_bits(0xfff8000000001234), "-nan"),
            (f64::from_bits(0x7ff0000000001234), "nan"),
        ] {
            assert_eq!(format_cx_coordinate(value), expected);
        }
    }
    #[test]
    fn two_d_positions_retain_duplicates_and_empty_selection() {
        let row = Conformer2D::new(9, vec![[1.0, 2.0], [3.0, 4.0]]);
        let source = CoordinateSource::TwoD(&row);
        assert_eq!(
            write_cx_coordinates_from_source(
                source,
                &[AtomId::new(1), AtomId::new(1), AtomId::new(0)]
            )
            .unwrap(),
            "3,4,;3,4,;1,2,"
        );
        assert_eq!(write_cx_coordinates_from_source(source, &[]).unwrap(), "");
    }
    #[test]
    fn is3d_flag_controls_z_independently_of_three_column_storage() {
        for (flag, expected) in [(false, "1,2,"), (true, "1,2,-nan")] {
            let row = Conformer3D::new(
                3,
                vec![[1.0, 2.0, f64::from_bits(0xfff8000000001234)]],
                flag,
            );
            assert_eq!(
                write_cx_coordinates_from_source(CoordinateSource::ThreeD(&row), &[AtomId::new(0)])
                    .unwrap(),
                expected
            );
        }
    }
    #[test]
    fn invalid_position_has_structural_source_range_error() {
        let row = Conformer2D::new(0, vec![[1.0, 2.0]]);
        assert!(matches!(
            write_cx_coordinates_from_source(CoordinateSource::TwoD(&row), &[AtomId::new(1)]),
            Err(SmilesParseError::CxCoordinateAtomOutOfRange { atom_count: 1, .. })
        ));
    }
    #[test]
    fn source_zero_threshold_six_digit_format_and_z_string_guard_are_shared() {
        let row = Conformer3D::new(
            0,
            vec![[999999.9, 1e-4, -0.0], [-1e6, 0.0000999, 1e-4]],
            true,
        );
        assert_eq!(
            write_cx_coordinates_from_source(
                CoordinateSource::ThreeD(&row),
                &[AtomId::new(0), AtomId::new(1)]
            )
            .unwrap(),
            "1e+06,0.0001,;-1e+06,0,0.0001"
        );
    }
}

#[cfg(test)]
mod zero_small_source_tests {
    use super::*;
    #[test]
    fn signed_zeros_and_subnormals_return_positive_zero() {
        for value in [0.0, -0.0, f64::from_bits(1), -f64::from_bits(1)] {
            assert_eq!(zero_small_cx_coordinate(value).to_bits(), 0.0f64.to_bits());
        }
    }
    #[test]
    fn threshold_is_strict_on_both_sides_and_signs() {
        let edge = 1e-4f64;
        let below = f64::from_bits(edge.to_bits() - 1);
        let above = f64::from_bits(edge.to_bits() + 1);
        for value in [below, -below] {
            assert_eq!(zero_small_cx_coordinate(value).to_bits(), 0.0f64.to_bits());
        }
        for value in [edge, -edge, above, -above] {
            assert_eq!(zero_small_cx_coordinate(value).to_bits(), value.to_bits());
        }
    }
    #[test]
    fn ordinary_and_max_finite_values_keep_original_bits() {
        for value in [1.0, -1.0, 42.25, -42.25, f64::MAX, -f64::MAX] {
            assert_eq!(zero_small_cx_coordinate(value).to_bits(), value.to_bits());
        }
    }
    #[test]
    fn infinities_and_nan_payloads_return_original_input_bits() {
        for value in [
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8000000001234),
            f64::from_bits(0xfff8000000001234),
            f64::from_bits(0x7ff0000000001234),
        ] {
            assert_eq!(zero_small_cx_coordinate(value).to_bits(), value.to_bits());
        }
    }
}

#[cfg(test)]
mod quote_cx_atom_property_source_tests {
    use super::*;
    #[test]
    fn empty_source_string_stays_empty() {
        assert!(quote_cx_atom_property(b"").is_empty());
    }
    #[test]
    fn only_periods_are_escaped_and_existing_entities_are_unchanged() {
        assert_eq!(
            quote_cx_atom_property(b"a..b&#46;:c,$|\\"),
            PropertyText::from("a&#46;&#46;b&#46;:c,$|\\")
        );
    }
    #[test]
    fn every_other_counted_byte_including_nul_and_nonutf8_is_preserved() {
        let input = (0..=255).collect::<Vec<u8>>();
        let mut expected = input[..46].to_vec();
        expected.extend_from_slice(b"&#46;");
        expected.extend_from_slice(&input[47..]);
        assert_eq!(
            quote_cx_atom_property(&input).as_bytes(),
            expected.as_slice()
        );
    }
}

#[cfg(test)]
mod atom_properties_source_tests {
    use super::*;
    use cosmolkit_model::{
        Atom, AtomQueryPredicate, AtomSpec, QueryAtom, QueryAtomIdentity, QueryNode,
    };
    use cosmolkit_types::Element;
    fn atom(i: usize, n: u8) -> QueryAtom {
        QueryAtom::from_identity_parts(
            AtomId::new(i),
            QueryAtomIdentity::from_atomic_number(n),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        )
    }
    fn write(atoms: &[QueryAtom]) -> PropertyText {
        write_cx_atom_properties_impl(
            atoms,
            &(0..atoms.len()).map(AtomId::new).collect::<Vec<_>>(),
        )
        .unwrap()
    }
    #[test]
    fn empty_order_and_no_properties_have_no_prefix() {
        assert!(
            write_cx_atom_properties_impl::<QueryAtom>(&[], &[])
                .unwrap()
                .is_empty()
        );
        assert!(write(&[atom(0, 6)]).is_empty());
    }
    #[test]
    fn source_insertion_order_and_overwrite_position_are_retained() {
        let mut a = atom(0, 6);
        a.set_prop("z", "old").unwrap();
        a.set_prop("a", "first").unwrap();
        a.set_prop("z", "new").unwrap();
        assert_eq!(
            write(&[a]),
            PropertyText::from("atomProp:0.z.new:0.a.first")
        );
    }
    #[test]
    fn exact_skip_private_and_computed_flags_filter_before_output() {
        let mut a = atom(0, 6);
        for key in [
            "atomLabel",
            "molFileValue",
            "molParity",
            "molAtomMapNumber",
            "molStereoCare",
            "molRxnExachg",
            "molInversionFlag",
            "_private",
        ] {
            a.set_prop(key, "ignored").unwrap();
        }
        a.set_computed_prop("derived", "ignored").unwrap();
        a.set_prop("visible", "value").unwrap();
        assert_eq!(write(&[a]), PropertyText::from("atomProp:0.visible.value"));
    }
    #[test]
    fn dummy_special_values_skip_on_any_actual_atomic_identity_but_empty_does_not() {
        for n in [0, 6, 119] {
            for value in ["*", "Pol", "Mod"] {
                let mut a = atom(0, n);
                a.set_prop("dummyLabel", value).unwrap();
                assert!(write(&[a]).is_empty());
            }
        }
        let mut a = atom(0, 0);
        a.set_prop("dummyLabel", "").unwrap();
        assert_eq!(write(&[a]), PropertyText::from("atomProp:0.dummyLabel."));
    }
    #[test]
    fn attachment_uses_property_presence_without_integer_conversion() {
        for (n, expected) in [(0, ""), (6, "atomProp:0.dummyLabel.custom")] {
            let mut a = atom(0, n);
            a.set_prop("_fromAttchpt", "bad-int").unwrap();
            a.set_prop("dummyLabel", "custom").unwrap();
            assert_eq!(write(&[a]), PropertyText::from(expected));
        }
    }
    #[test]
    fn counter_advances_for_empty_rows_and_repeated_atom_ids() {
        let mut a = atom(0, 6);
        a.set_prop("key", "v").unwrap();
        let atoms = [a, atom(1, 6)];
        assert_eq!(
            write_cx_atom_properties_impl(
                &atoms,
                &[AtomId::new(0), AtomId::new(1), AtomId::new(0)]
            )
            .unwrap(),
            PropertyText::from("atomProp:0.key.v:2.key.v")
        );
    }
    #[test]
    fn counted_binary_property_keys_and_values_escape_only_periods() {
        let mut a = atom(0, 6);
        a.set_prop(
            PropertyText::from(vec![255, 0, b'.']),
            PropertyValue::String(PropertyText::from(vec![128, b'.', 0])),
        )
        .unwrap();
        let mut expected = PropertyText::from("atomProp:0.");
        expected.extend_bytes(&[255, 0]);
        expected.extend_bytes(b"&#46;.");
        expected.push_byte(128);
        expected.extend_bytes(b"&#46;");
        expected.push_byte(0);
        assert_eq!(write(&[a]), expected);
    }
    #[test]
    fn all_seven_modeled_property_tags_use_source_string_projection() {
        let mut a = atom(0, 6);
        for (k, v) in [
            ("s", PropertyValue::String("x".into())),
            ("i", PropertyValue::Int(-7)),
            ("u", PropertyValue::UInt(u32::MAX)),
            ("d", PropertyValue::Double(-0.0)),
            ("b", PropertyValue::Bool(true)),
            ("iv", PropertyValue::IntVector(vec![1, -2])),
            (
                "sv",
                PropertyValue::StringVector(vec!["a".into(), "b".into()]),
            ),
        ] {
            a.set_prop(k, v).unwrap();
        }
        assert_eq!(
            write(&[a]),
            PropertyText::from(
                "atomProp:0.s.x:0.i.-7:0.u.4294967295:0.d.-0:0.b.1:0.iv.[1,-2]:0.sv.[a,b]"
            )
        );
    }
    #[test]
    fn malformed_computed_list_is_not_replaced_by_empty_filter() {
        let mut a = atom(0, 6);
        a.set_prop("__computedProps", "bad-tag").unwrap();
        a.set_prop("molAtomMapNumber", "skipped").unwrap();
        assert!(matches!(
            write_cx_atom_properties_impl(&[a], &[AtomId::new(0)]),
            Err(SmilesParseError::WriterPropertyList { .. })
        ));
    }
    #[test]
    fn invalid_atom_order_returns_source_structural_range_error() {
        assert!(matches!(
            write_cx_atom_properties_impl::<QueryAtom>(&[], &[AtomId::new(0)]),
            Err(SmilesParseError::CxAtomPropertyAtomOutOfRange { atom_count: 0, .. })
        ));
    }
    #[test]
    fn ordinary_atom_consumer_uses_the_same_single_body() {
        let mut a = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
        a.set_prop("a.b", "c.d").unwrap();
        assert_eq!(
            write_cx_atom_properties_impl(&[a], &[AtomId::new(0)]).unwrap(),
            PropertyText::from("atomProp:0.a&#46;b.c&#46;d")
        );
    }
}

#[cfg(test)]
mod typed_bonds_source_tests {
    use super::*;
    use cosmolkit_core::stereo_graph::BondRows;
    use cosmolkit_model::{BondSpec, QueryBond};
    fn bonds(orders: &[BondOrder]) -> Vec<Bond> {
        orders
            .iter()
            .enumerate()
            .map(|(i, &o)| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(i), AtomId::new(i + 1), o),
                )
            })
            .collect()
    }
    fn write(
        bonds: &[Bond],
        ao: &[usize],
        bo: &[usize],
        order: BondOrder,
        symbol: &[u8],
    ) -> PropertyText {
        write_cx_coord_or_hydrogen_bonds_source(
            BondRows::Concrete(bonds),
            &ao.iter().map(|&i| AtomId::new(i)).collect::<Vec<_>>(),
            &bo.iter().map(|&i| BondId::new(i)).collect::<Vec<_>>(),
            order,
            symbol,
        )
        .unwrap()
    }
    #[test]
    fn empty_and_unselected_order_never_emit_bare_prefix() {
        assert!(write(&[], &[], &[], BondOrder::Dative, b"C").is_empty());
        assert!(
            write(
                &bonds(&[BondOrder::Single]),
                &[0, 1],
                &[0],
                BondOrder::Dative,
                b"C"
            )
            .is_empty()
        );
    }
    #[test]
    fn bond_type_comparison_is_exact_without_dative_one_folding() {
        let b = bonds(&[
            BondOrder::DativeOne,
            BondOrder::Dative,
            BondOrder::Hydrogen,
            BondOrder::Single,
        ]);
        assert_eq!(
            write(&b, &[0, 1, 2, 3, 4], &[3, 0, 2, 1], BondOrder::Dative, b"C"),
            PropertyText::from("C:1.3")
        );
        assert_eq!(
            write(
                &b,
                &[0, 1, 2, 3, 4],
                &[3, 0, 2, 1],
                BondOrder::Hydrogen,
                b"H"
            ),
            PropertyText::from("H:2.2")
        );
        assert_eq!(
            write(
                &b,
                &[0, 1, 2, 3, 4],
                &[3, 0, 2, 1],
                BondOrder::DativeOne,
                b"X"
            ),
            PropertyText::from("X:0.1")
        );
    }
    #[test]
    fn first_duplicate_atom_position_matches_std_find() {
        assert_eq!(
            write(
                &bonds(&[BondOrder::Dative]),
                &[1, 0, 0],
                &[0],
                BondOrder::Dative,
                b"C"
            ),
            PropertyText::from("C:1.0")
        );
    }
    #[test]
    fn missing_start_atom_is_end_distance_even_for_empty_order() {
        let b = bonds(&[BondOrder::Hydrogen]);
        assert_eq!(
            write(&b, &[1], &[0], BondOrder::Hydrogen, b"H"),
            PropertyText::from("H:1.0")
        );
        assert_eq!(
            write(&b, &[], &[0], BondOrder::Hydrogen, b"H"),
            PropertyText::from("H:0.0")
        );
    }
    #[test]
    fn repeated_bonds_keep_output_positions_including_filtered_rows() {
        assert_eq!(
            write(
                &bonds(&[BondOrder::Dative, BondOrder::Single]),
                &[0, 1, 2],
                &[1, 0, 0],
                BondOrder::Dative,
                b"C"
            ),
            PropertyText::from("C:0.1,0.2")
        );
    }
    #[test]
    fn actual_oriented_begin_is_not_canonicalized_to_lower_endpoint() {
        let b = [Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Dative),
        )];
        assert_eq!(
            write(&b, &[0, 1], &[0], BondOrder::Dative, b"C"),
            PropertyText::from("C:1.0")
        );
    }
    #[test]
    fn counted_binary_symbol_is_preserved_once_without_escaping() {
        let b = bonds(&[BondOrder::Dative, BondOrder::Dative]);
        let mut expected = PropertyText::from(vec![255, 0, b'.', b':']);
        expected.extend_bytes(b"0.0,1.1");
        assert_eq!(
            write(&b, &[0, 1, 2], &[0, 1], BondOrder::Dative, &[255, 0, b'.']),
            expected
        );
        assert_eq!(
            write(&b, &[0, 1, 2], &[0], BondOrder::Dative, b""),
            PropertyText::from(":0.0")
        );
    }
    #[test]
    fn ordinary_and_query_carriers_share_one_body_and_propagate_range_error() {
        let b = bonds(&[BondOrder::Dative]);
        let q = [QueryBond::from_parts(
            b[0].clone(),
            cosmolkit_model::QueryNode::predicate(cosmolkit_model::BondQueryPredicate::Any),
        )];
        let ao = [AtomId::new(0), AtomId::new(1)];
        let bo = [BondId::new(0)];
        assert_eq!(
            write_cx_coord_or_hydrogen_bonds_source(
                BondRows::Concrete(&b),
                &ao,
                &bo,
                BondOrder::Dative,
                b"C"
            )
            .unwrap(),
            write_cx_coord_or_hydrogen_bonds_source(
                BondRows::Query(&q),
                &ao,
                &bo,
                BondOrder::Dative,
                b"C"
            )
            .unwrap()
        );
        assert!(matches!(
            write_cx_coord_or_hydrogen_bonds_source(
                BondRows::Query(&q),
                &ao,
                &[BondId::new(1)],
                BondOrder::Dative,
                b"C"
            ),
            Err(SmilesParseError::CxTypedBondOutOfRange { bond_count: 1, .. })
        ));
    }
}

#[cfg(test)]
mod zero_bonds_owner_source_tests {
    use super::*;
    use cosmolkit_core::stereo_graph::BondRows;
    use cosmolkit_model::{BondQueryPredicate, BondSpec, QueryBond, QueryNode};
    fn bonds(orders: &[BondOrder]) -> Vec<Bond> {
        orders
            .iter()
            .enumerate()
            .map(|(i, &o)| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(i), AtomId::new(i + 1), o),
                )
            })
            .collect()
    }
    fn write(bonds: &[Bond], bo: &[usize]) -> PropertyText {
        write_cx_zero_bonds_source(
            BondRows::Concrete(bonds),
            &bo.iter().map(|&i| BondId::new(i)).collect::<Vec<_>>(),
        )
        .unwrap()
    }
    #[test]
    fn empty_and_unselected_bonds_have_no_prefix() {
        assert!(write(&[], &[]).is_empty());
        assert!(write(&bonds(&[BondOrder::Single]), &[0]).is_empty());
    }
    #[test]
    fn repeated_selected_bonds_keep_positions_after_filtered_rows() {
        assert_eq!(
            write(
                &bonds(&[BondOrder::Zero, BondOrder::Zero, BondOrder::Single]),
                &[2, 1, 0, 1]
            ),
            PropertyText::from("Z:1,2,3")
        );
    }
    #[test]
    fn all_zero_bonds_emit_one_prefix_and_no_trailing_comma() {
        assert_eq!(
            write(
                &bonds(&[BondOrder::Zero, BondOrder::Zero, BondOrder::Zero]),
                &[0, 1, 2]
            ),
            PropertyText::from("Z:0,1,2")
        );
    }
    #[test]
    fn other_modeled_bond_types_are_not_zero_aliases() {
        assert!(
            write(
                &bonds(&[
                    BondOrder::Single,
                    BondOrder::Dative,
                    BondOrder::Hydrogen,
                    BondOrder::DativeOne,
                    BondOrder::Aromatic
                ]),
                &[0, 1, 2, 3, 4]
            )
            .is_empty()
        );
    }
    #[test]
    fn actual_query_carrier_order_wins_without_reading_or_dropping_predicate() {
        let b = bonds(&[BondOrder::Zero]);
        let q = [QueryBond::from_parts(
            b[0].clone(),
            QueryNode::predicate(BondQueryPredicate::Any),
        )];
        let before = q.clone();
        assert_eq!(
            write_cx_zero_bonds_source(BondRows::Query(&q), &[BondId::new(0)]).unwrap(),
            PropertyText::from("Z:0")
        );
        assert_eq!(q, before);
    }
    #[test]
    fn invalid_source_bond_is_error_before_type_filter() {
        assert!(matches!(
            write_cx_zero_bonds_source(
                BondRows::Concrete(&bonds(&[BondOrder::Single])),
                &[BondId::new(1)]
            ),
            Err(SmilesParseError::CxZeroBondOutOfRange { bond_count: 1, .. })
        ));
    }
}

#[cfg(test)]
mod link_nodes_owner_source_tests {
    use super::*;
    use cosmolkit_model::{
        AtomQueryPredicate, BondSpec, QueryAtom, QueryAtomIdentity, QueryBond, QueryGraph,
        QueryNode,
    };
    fn graph(raw: Option<&[u8]>, branched: bool) -> QueryGraph {
        let atoms = (0..4)
            .map(|i| {
                QueryAtom::from_identity_parts(
                    AtomId::new(i),
                    QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                )
            })
            .collect();
        let edges = if branched {
            vec![(1, 0), (1, 2), (1, 3)]
        } else {
            vec![(1, 0), (1, 2)]
        };
        QueryGraph::from_parts(
            atoms,
            edges
                .into_iter()
                .enumerate()
                .map(|(i, (a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            raw.map(|v| {
                (
                    PropertyText::from("_molLinkNodes"),
                    PropertyValue::String(PropertyText::from_bytes(v)),
                )
            }),
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn run(q: &QueryGraph, ao: &[usize]) -> (Result<String, SmilesParseError>, Vec<Vec<u8>>) {
        let mut warnings = vec![];
        let output = write_query_cx_link_nodes_source(
            q,
            &ao.iter().copied().map(AtomId::new).collect::<Vec<_>>(),
            &mut |v| warnings.push(v.to_vec()),
        );
        (output, warnings)
    }
    #[test]
    fn absent_or_empty_nodes_return_before_invalid_reverse_order() {
        for raw in [None, Some(&b"|||"[..])] {
            let (out, w) = run(&graph(raw, false), &[99]);
            assert_eq!(out.unwrap(), "");
            assert!(w.is_empty());
        }
    }
    #[test]
    fn degree_two_uses_forward_order_and_does_not_read_outer_order() {
        let q = graph(Some(b"1 3 2 2 1 2 3"), false);
        assert_eq!(run(&q, &[3, 0, 2, 1]).0.unwrap(), "LN:0:1.3");
        assert_eq!(run(&q, &[0, 1]).0.unwrap(), "LN:1:1.3");
    }
    #[test]
    fn degree_three_includes_source_outer_pair_in_property_order() {
        let q = graph(Some(b"1 3 2 2 3 2 1"), true);
        assert_eq!(run(&q, &[3, 0, 2, 1]).0.unwrap(), "LN:0:1.3.2.3");
    }
    #[test]
    fn repeated_nodes_and_empty_delimiter_tokens_preserve_output_order() {
        let q = graph(Some(b"|1 3 2 2 1 2 3||2 4 2 2 1 2 3|"), false);
        assert_eq!(run(&q, &[0, 1, 2, 3]).0.unwrap(), "LN:1:1.3,1:2.4");
    }
    #[test]
    fn all_nonstrict_rejections_warn_then_continue_to_valid_node() {
        let cases: [(&[u8], &[u8]); 7] = [
            (
                b"x",
                b"Cannot convert values in LINKNODE 'x' to unsigned ints",
            ),
            (b"1 3 2 2 1", b"not enough values in LINKNODE '1 3 2 2 1'"),
            (b"0 3 2 2 1 2 3", b"bad counts in LINKNODE '0 3 2 2 1 2 3'"),
            (b"4 3 2 2 1 2 3", b"bad counts in LINKNODE '4 3 2 2 1 2 3'"),
            (
                b"1 3 3 2 1 2 3 2 4",
                b"only link nodes with 2 bonds are currently supported",
            ),
            (
                b"1 3 2 2 1 3 2",
                b"bonds don't start at the same atom for LINKNODE '1 3 2 2 1 3 2'",
            ),
            (
                b"1 3 2 2 1 2 4",
                b"bond not found between atoms in LINKNODE '1 3 2 2 1 2 4'",
            ),
        ];
        for (raw, expected) in cases {
            let mut v = raw.to_vec();
            v.extend_from_slice(b"|1 3 2 2 1 2 3");
            let (out, w) = run(&graph(Some(&v), false), &[0, 1, 2, 3]);
            assert_eq!(out.unwrap(), "LN:1:1.3");
            assert_eq!(w, vec![expected.to_vec()]);
        }
    }
    #[test]
    fn binary_warning_bytes_are_not_utf8_normalized() {
        let (out, w) = run(&graph(Some(b"\xff\0"), false), &[]);
        assert_eq!(out.unwrap(), "");
        assert_eq!(
            w,
            vec![b"Cannot convert values in LINKNODE '\xff\0' to unsigned ints".to_vec()]
        );
    }
    #[test]
    fn warning_prefix_is_emitted_before_later_hard_range_error() {
        let (out, w) = run(&graph(Some(b"x|1 3 2 2 0 2 3"), false), &[0, 1, 2, 3]);
        assert!(matches!(out,Err(SmilesParseError::Cx(ref s))if s.contains("4294967295")));
        assert_eq!(w.len(), 1);
        assert_eq!(
            w[0],
            b"Cannot convert values in LINKNODE 'x' to unsigned ints"
        );
    }
    #[test]
    fn reverse_order_writes_are_checked_only_after_nonempty_nodes() {
        let (out, w) = run(&graph(Some(b"1 3 2 2 1 2 3"), false), &[99]);
        assert!(matches!(
            out,
            Err(SmilesParseError::CxLinkAtomOutOfRange {
                atom: 99,
                atom_count: 4
            })
        ));
        assert!(w.is_empty());
    }
    #[test]
    fn missing_forward_order_is_error_at_reached_source_read() {
        for (ao, index) in [(vec![0], 1), (vec![0, 1], 2)] {
            let (out, _) = run(&graph(Some(b"1 3 2 2 1 2 3"), true), &ao);
            assert!(
                matches!(out,Err(SmilesParseError::CxLinkOrderOutOfRange{index:i,..})if i==index)
            );
        }
    }
    #[test]
    fn scalar_property_uses_canonical_string_conversion() {
        let mut q = graph(None, false);
        q.set_prop("_molLinkNodes", PropertyValue::UInt(123))
            .unwrap();
        let (out, w) = run(&q, &[]);
        assert_eq!(out.unwrap(), "");
        assert_eq!(w, vec![b"not enough values in LINKNODE '123'".to_vec()]);
    }
    #[test]
    fn actual_query_identity_predicates_properties_and_cache_are_immutable() {
        let q = graph(Some(b"1 3 2 2 1 2 3"), true);
        let before = q.clone();
        assert_eq!(run(&q, &[0, 1, 2, 3]).0.unwrap(), "LN:1:1.3.0.2");
        assert_eq!(q, before);
    }
}
