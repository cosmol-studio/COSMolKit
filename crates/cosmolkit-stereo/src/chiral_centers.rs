//! Read-only modern RDKit chiral-center query over detached values.

use std::borrow::Cow;

use cosmolkit_core::{
    PotentialStereoCenter, PotentialStereoDescriptor, PotentialStereoParams,
    PotentialStereoSpecified, PotentialStereoType, RingInfo, RingSearchParams, ValenceAssignment,
    ValenceModel, ValenceParams, assign_valence, potential_stereo, symmetrized_sssr,
    symmetrized_sssr_with_existing,
};
use cosmolkit_model::{MoleculeProperties, TopologyBlock};

use crate::{
    CipLabelOptions, StereoReadError, cip_labels::assign_cip_labels_cow_with_prepared_state,
};

/// Finds tetrahedral stereocenters using the modern perception and CIP owners.
///
/// This projects pinned RDKit `FindMolChiralCenters` with `includeCIP=True`
/// and `useLegacyImplementation=False`. Candidate order and lowercase CIP
/// labels are preserved. Neither the input topology nor its properties change.
pub fn find_chiral_centers(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    include_unassigned: bool,
) -> Result<Vec<(usize, String)>, StereoReadError> {
    find_chiral_centers_with_prepared_state(topology, properties, include_unassigned, None, None)
}

/// Read-only query reusing explicitly supplied valid state for the same topology.
#[doc(hidden)]
pub fn find_chiral_centers_with_prepared_state(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    include_unassigned: bool,
    source_valence: Option<&ValenceAssignment>,
    source_rings: Option<&RingInfo>,
) -> Result<Vec<(usize, String)>, StereoReadError> {
    // Source: RDKit 2026.03.1, 351f8f378f8ad6bbd517980c38896e66bf907af8,
    // rdkit/Chem/__init__.py, FindMolChiralCenters modern branch.
    // BEGIN RDKIT PYTHON modern candidate perception and CIP selection
    // RDKit✔️❌:       centers = []
    // RDKit✔️❌:       itms = FindPotentialStereo(mol)
    // RDKit✔️❌:       if includeCIP:
    // RDKit✔️❌:         atomsToLabel = []
    // RDKit✔️❌:         bondsToLabel = []
    // RDKit✔️❌:         for si in itms:
    // RDKit✔️❌:           if si.type == StereoType.Atom_Tetrahedral:
    // RDKit✔️❌:             atomsToLabel.append(si.centeredOn)
    // RDKit✔️❌:           elif si.type == StereoType.Bond_Double:
    // RDKit✔️❌:             bondsToLabel.append(si.centeredOn)
    // RDKit✔️❌:         AssignCIPLabels(mol, atomsToLabel=atomsToLabel, bondsToLabel=bondsToLabel)
    // END RDKIT PYTHON modern candidate perception and CIP selection
    // FindPotentialStereo prepares a non-strict property cache and symmetric
    // rings. Reuse their detached owners instead of using raw chiral tags as
    // potential-center evidence. All source writes stay in local COW values.
    // Perception retains its necessary detached working graph. Valid caller
    // valence/SymmSSSR are borrowed, not recalculated or cloned; the same actual
    // rows are passed into the modern CIP molecule shared by all configurations.
    let (valence, rings) = prepare_chiral_state(topology, source_valence, source_rings)?;
    let perception = potential_stereo(
        topology,
        &valence,
        &rings,
        &PotentialStereoParams::default(),
    )?;
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    for info in &perception.stereo {
        match (info.stereo_type, info.centered_on) {
            (PotentialStereoType::AtomTetrahedral, PotentialStereoCenter::Atom(atom)) => {
                atoms.push(atom)
            }
            (PotentialStereoType::BondDouble, PotentialStereoCenter::Bond(bond)) => {
                bonds.push(bond)
            }
            _ => {}
        }
    }
    let options = CipLabelOptions::default()
        .with_atoms(atoms)
        .with_bonds(bonds);
    let (labelled, _properties, _completion) = assign_cip_labels_cow_with_prepared_state(
        Cow::Borrowed(topology),
        Cow::Borrowed(properties),
        &options,
        Some(&valence),
        Some(&rings),
    )?;

    // BEGIN RDKIT PYTHON modern center filtering and source-defined label fallback
    // RDKit✔️✔️:       for si in itms:
    // RDKit✔️✔️:         if si.type == StereoType.Atom_Tetrahedral and (includeUnassigned or si.specified
    // RDKit✔️✔️:                                                        == StereoSpecified.Specified):
    // RDKit✔️✔️:           idx = si.centeredOn
    // RDKit✔️✔️:           atm = mol.GetAtomWithIdx(idx)
    // RDKit✔️✔️:           if includeCIP and atm.HasProp("_CIPCode"):
    // RDKit✔️✔️:             code = atm.GetProp("_CIPCode")
    // RDKit✔️✔️:           else:
    // RDKit✔️✔️:             if si.specified:
    // RDKit✔️✔️:               code = str(si.descriptor)
    // RDKit✔️✔️:             else:
    // RDKit✔️✔️:               code = '?'
    // RDKit✔️✔️:               atm.SetIntProp('_ChiralityPossible', 1)
    // RDKit✔️✔️:           centers.append((idx, code))
    // END RDKIT PYTHON modern center filtering and source-defined label fallback
    // This query returns only the detached rows, so the upstream's temporary
    // _ChiralityPossible property is not installed on the read-only input.
    // Unknown (enum value 2) is truthy upstream: without CIP it returns NoValue,
    // not '?'. Only Specified rows survive when include_unassigned is false.
    let mut centers = Vec::new();
    for info in perception.stereo {
        if info.stereo_type != PotentialStereoType::AtomTetrahedral
            || (!include_unassigned && info.specified != PotentialStereoSpecified::Specified)
        {
            continue;
        }
        let PotentialStereoCenter::Atom(atom) = info.centered_on else {
            unreachable!("tetrahedral perception uses atom identifiers")
        };
        let code = if let Some(value) = labelled.atoms[atom.index()].prop("_CIPCode") {
            // GetProp uses canonical lexical conversion of the newly assigned
            // selected label. Candidate CIP clears stale selected labels first.
            let text = cosmolkit_core::property_value_to_string(value)?;
            String::from_utf8(text.into_bytes())
                .map_err(|source| StereoReadError::CipLabelEncoding { atom, source })?
        } else if info.specified == PotentialStereoSpecified::Unspecified {
            "?".to_owned()
        } else {
            match info.descriptor {
                PotentialStereoDescriptor::None => "NoValue",
                PotentialStereoDescriptor::TetrahedralClockwise => "Tet_CW",
                PotentialStereoDescriptor::TetrahedralCounterclockwise => "Tet_CCW",
                PotentialStereoDescriptor::BondCis => "Bond_Cis",
                PotentialStereoDescriptor::BondTrans => "Bond_Trans",
                PotentialStereoDescriptor::BondAtropCw => "Bond_AtropCW",
                PotentialStereoDescriptor::BondAtropCcw => "Bond_AtropCCW",
            }
            .to_owned()
        };
        centers.push((atom.index(), code));
    }
    Ok(centers)
}

fn prepare_chiral_state<'a>(
    topology: &'a TopologyBlock,
    source_valence: Option<&'a ValenceAssignment>,
    source_rings: Option<&'a RingInfo>,
) -> Result<(Cow<'a, ValenceAssignment>, Cow<'a, RingInfo>), StereoReadError> {
    // BEGIN RDKIT CPP FUNCTION FindStereo.cpp::findPotentialStereo cache guards
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit✔️✔️:     MolOps::symmetrizeSSSR(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (mol.needsUpdatePropertyCache()) {
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION FindStereo.cpp::findPotentialStereo cache guards
    // Preserve the source ring-then-property-cache preparation order. A
    // valid-gated cache is supplied by the facade, not inferred from atom
    // contents. Missing state
    // goes through the original owners. When an insufficient ring cache already
    // owns RDL families, a private working cache preserves that source state.
    let rings = match source_rings {
        Some(rings) if rings.is_symm_sssr() => Cow::Borrowed(rings),
        Some(rings) if rings.are_ring_families_initialized() => {
            let mut prepared = rings.clone();
            symmetrized_sssr_with_existing(
                topology,
                &mut prepared,
                false,
                &RingSearchParams::default(),
            )?;
            Cow::Owned(prepared)
        }
        _ => Cow::Owned(symmetrized_sssr(topology, &RingSearchParams::default())?),
    };
    let valence = match source_valence {
        Some(valence) => Cow::Borrowed(valence),
        None => Cow::Owned(assign_valence(
            topology,
            &ValenceParams {
                model: ValenceModel::RdkitLike,
                strict: false,
            },
        )?),
    };
    Ok((valence, rings))
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec};
    use cosmolkit_types::{BondDirection, BondOrder, ChiralTag, Element};

    #[test]
    fn chiral_queries_borrow_prepared_state_and_rebuild_only_missing_quality() {
        let topology = tetrahedron(ChiralTag::TetrahedralCcw, BondDirection::None);
        let valence = assign_valence(
            &topology,
            &ValenceParams {
                model: ValenceModel::RdkitLike,
                strict: false,
            },
        )
        .unwrap();
        let rings = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
        let before = (topology.clone(), valence.clone(), rings.clone());
        let (borrowed_valence, borrowed_rings) =
            prepare_chiral_state(&topology, Some(&valence), Some(&rings)).unwrap();
        assert!(matches!(borrowed_valence, Cow::Borrowed(value) if std::ptr::eq(value, &valence)));
        assert!(matches!(borrowed_rings, Cow::Borrowed(value) if std::ptr::eq(value, &rings)));
        let properties = MoleculeProperties::default();
        for _ in 0..2 {
            assert_eq!(
                find_chiral_centers_with_prepared_state(
                    &topology,
                    &properties,
                    false,
                    Some(&valence),
                    Some(&rings),
                )
                .unwrap(),
                [(0, "S".to_owned())],
            );
        }
        let fast = cosmolkit_core::fast_find_rings(&topology).unwrap();
        let (_, prepared_rings) =
            prepare_chiral_state(&topology, Some(&valence), Some(&fast)).unwrap();
        assert!(matches!(prepared_rings, Cow::Owned(_)));
        assert!(prepared_rings.is_symm_sssr());
        assert_eq!((topology, valence, rings), before);
    }

    #[test]
    fn chiral_query_keeps_nonstrict_property_cache_for_later_cip_reads() {
        // Pinned RDKit 2026.09.1: FindMolChiralCenters on unsanitized
        // [C@](F)(Cl)(Br)C(F)(F)(F)F returns [(0, 'R')]. FindPotentialStereo
        // prepares non-strict valence (atom 4 has explicit valence 5); Digraph
        // then reads that actual atom cache instead of sanitizing all atoms.
        let atoms = [
            Element::C,
            Element::F,
            Element::CL,
            Element::BR,
            Element::C,
            Element::F,
            Element::F,
            Element::F,
            Element::F,
        ]
        .into_iter()
        .enumerate()
        .map(|(index, element)| {
            let spec = AtomSpec::new(element);
            Atom::from_spec(
                AtomId::new(index),
                if index == 0 {
                    spec.with_chiral_tag(ChiralTag::TetrahedralCcw)
                } else {
                    spec
                },
            )
        })
        .collect();
        let bonds = [
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 4),
            (4, 5),
            (4, 6),
            (4, 7),
            (4, 8),
        ]
        .into_iter()
        .enumerate()
        .map(|(index, (begin, end))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            )
        })
        .collect();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        let before = topology.clone();
        assert_eq!(
            find_chiral_centers(&topology, &MoleculeProperties::default(), false).unwrap(),
            [(0, "R".to_owned())],
        );
        assert_eq!(topology, before);
    }

    fn tetrahedron(tag: ChiralTag, direction: BondDirection) -> TopologyBlock {
        let atoms = [Element::C, Element::F, Element::CL, Element::BR, Element::I]
            .into_iter()
            .enumerate()
            .map(|(index, element)| {
                let spec = AtomSpec::new(element);
                Atom::from_spec(
                    AtomId::new(index),
                    if index == 0 {
                        spec.with_chiral_tag(tag)
                    } else {
                        spec
                    },
                )
            })
            .collect();
        let bonds = (1..5)
            .map(|index| {
                Bond::from_spec(
                    BondId::new(index - 1),
                    BondSpec::new(AtomId::new(0), AtomId::new(index), BondOrder::Single)
                        .with_direction(if index == 1 {
                            direction
                        } else {
                            BondDirection::None
                        }),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    #[test]
    fn chiral_centers_exclude_ordinary_atoms_and_empty_graphs() {
        let empty =
            TopologyBlock::try_from_parts(Vec::new(), Vec::new(), Vec::new(), Vec::new()).unwrap();
        let methane = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap();
        for topology in [empty, methane] {
            for include in [false, true] {
                assert_eq!(
                    find_chiral_centers(&topology, &MoleculeProperties::default(), include)
                        .unwrap(),
                    Vec::<(usize, String)>::new()
                );
            }
        }
    }

    #[test]
    fn chiral_centers_filter_potential_centers_before_cip_labels() {
        // Pinned RDKit 2026.03.1 modern FindMolChiralCenters on
        // [C@](F)(Cl)(Br)I / [C@@](F)(Cl)(Br)I / C(F)(Cl)(Br)I.
        for (tag, label) in [
            (ChiralTag::TetrahedralCcw, "S"),
            (ChiralTag::TetrahedralCw, "R"),
            (ChiralTag::Unspecified, "?"),
        ] {
            let topology = tetrahedron(tag, BondDirection::None);
            let before = topology.clone();
            let properties = MoleculeProperties::default().with_name("source");
            let properties_before = properties.clone();
            for include in [false, true] {
                let expected = if include || tag != ChiralTag::Unspecified {
                    vec![(0, label.to_owned())]
                } else {
                    Vec::new()
                };
                assert_eq!(
                    find_chiral_centers(&topology, &properties, include).unwrap(),
                    expected
                );
                assert_eq!(topology, before);
                assert_eq!(properties, properties_before);
            }
            assert!(topology.atoms[0].cip_descriptor().unwrap().is_none());
        }
    }

    #[test]
    fn chiral_centers_unknown_stereo_retains_source_no_value_fallback() {
        // RDKit's Unknown enum is truthy; its absent-CIP descriptor is NoValue.
        let topology = tetrahedron(ChiralTag::Unspecified, BondDirection::Unknown);
        assert_eq!(
            find_chiral_centers(&topology, &MoleculeProperties::default(), true).unwrap(),
            [(0, "NoValue".to_owned())]
        );
        assert!(
            find_chiral_centers(&topology, &MoleculeProperties::default(), false)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn chiral_centers_clear_stale_selected_property_strings_before_labeling() {
        use cosmolkit_model::PropertyValue;
        // Candidate CIPLabeler.cpp::clearSelectedCIPProperties clears every
        // selected potential center, including unassigned non-configurations.
        for value in [
            PropertyValue::String("foo".into()),
            PropertyValue::String("".into()),
            PropertyValue::Bool(true),
            PropertyValue::Int(17),
        ] {
            let mut topology = tetrahedron(ChiralTag::Unspecified, BondDirection::None);
            topology.atoms[0].set_prop("_CIPCode", value).unwrap();
            let before = topology.clone();
            assert_eq!(
                find_chiral_centers(&topology, &MoleculeProperties::default(), true).unwrap(),
                [(0, "?".to_owned())]
            );
            assert!(
                find_chiral_centers(&topology, &MoleculeProperties::default(), false)
                    .unwrap()
                    .is_empty()
            );
            assert_eq!(topology, before);
        }
    }

    #[test]
    fn chiral_centers_clear_invalid_stale_label_without_mutating_input() {
        use cosmolkit_model::PropertyValue;
        let mut topology = tetrahedron(ChiralTag::Unspecified, BondDirection::None);
        topology.atoms[0]
            .set_prop("_CIPCode", PropertyValue::String(vec![0xff].into()))
            .unwrap();
        let before = topology.clone();
        // No string decoding occurs: upstream clears this selected stale
        // property before configuration discovery. The input stays read-only.
        assert_eq!(
            find_chiral_centers(&topology, &MoleculeProperties::default(), true).unwrap(),
            [(0, "?".to_owned())]
        );
        assert_eq!(topology, before);
    }

    #[test]
    fn chiral_centers_errors_do_not_become_empty_results() {
        use std::error::Error as _;
        let topology = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::ThreeCenter),
            )],
            Vec::new(),
            Vec::new(),
        )
        .unwrap();
        let before = topology.clone();
        let error =
            find_chiral_centers(&topology, &MoleculeProperties::default(), true).unwrap_err();
        assert!(matches!(error, StereoReadError::Valence(_)), "{error:?}");
        assert!(error.source().is_some());
        assert_eq!(topology, before);
    }
}
