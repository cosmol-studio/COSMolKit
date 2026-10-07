//! Source CX vector insertion over canonical query graphs, retaining every AST.
use crate::{SmartsWriteError, SmartsWriteOutput};
use cosmolkit_model::{
    AtomId, BondId, Conformer3D, PropertyValue, QueryGraph, SubstanceGroupId, insert_stereo_groups,
    query_substance_groups, replace_query_substance_groups,
};
use cosmolkit_smiles::{CoordinateSource, CxCoordinateSelection, select_cx_coordinates_from_sets};
use std::collections::BTreeMap;

/// The same query graph type and writer evidence, with no second AST, live
/// molecule, cache authority or compatibility topology.
#[doc(hidden)]
pub struct QueryCxComposition {
    pub graph: QueryGraph,
    pub atom_order: Vec<AtomId>,
    pub bond_order: Vec<BondId>,
}

fn check_cx_features(query: &QueryGraph) {
    // RDKit❗✔️: void checkCXFeatures(const ROMol &mol) {
    // RDKit❗✔️:   std::string lns;
    // RDKit❗✔️:   if (mol.getPropIfPresent(common_properties::molFileLinkNodes, lns)) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "CX Extensions: mol has link nodes which are not currently supported"
    // RDKit❗✔️:         << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const auto &sgs = getSubstanceGroups(mol);
    // RDKit❗✔️:   auto parent_check =
    // RDKit❗✔️:       std::any_of(sgs.cbegin(), sgs.cend(), [&](const SubstanceGroup &sg) {
    // RDKit❗✔️:         if (sg.hasProp("PARENT")) {
    // RDKit❗✔️:           return true;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         return false;
    // RDKit❗✔️:       });
    // RDKit❗✔️:   if (parent_check) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "CX Extensions: Substance group hierarchy is not always preserved."
    // RDKit❗✔️:         << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    if query.prop("_molLinkNodes").is_some() {
        eprintln!("CX Extensions: mol has link nodes which are not currently supported");
    }
    if query_substance_groups(query)
        .iter()
        .any(|group| group.parent().is_some() || group.props().contains_key(b"PARENT".as_slice()))
    {
        eprintln!("CX Extensions: Substance group hierarchy is not always preserved.");
    }
}

#[doc(hidden)]
pub fn compose_query_cx_templates(
    templates: &[(&QueryGraph, &SmartsWriteOutput)],
    selections: &[CxCoordinateSelection],
    coordinates_enabled: bool,
) -> Result<QueryCxComposition, SmartsWriteError> {
    // RDKit❗✔️: std::string getCXExtensions(const std::vector<ROMol *> &mols,
    // RDKit❗✔️:                             std::uint32_t flags) {
    // RDKit❗✔️:   for (const auto &mol : mols) {
    // RDKit❗✔️:     checkCXFeatures(*mol);
    // RDKit❗✔️:     if (!mol->hasProp(RDKit::common_properties::_smilesAtomOutputOrder) ||
    // RDKit❗✔️:         !mol->hasProp(RDKit::common_properties::_smilesBondOutputOrder)) {
    // RDKit❗✔️:       throw ValueErrorException(
    // RDKit❗✔️:           "Input molecule does not have the required "
    // RDKit❗✔️:           "smiles ordering properties set");
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   RDKit::RWMol rwmol;
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<unsigned int> atomOrdering;
    // RDKit❗✔️:   std::vector<unsigned int> bondOrdering;
    // RDKit❗✔️:
    // RDKit❗✔️:   for (const auto &mol : mols) {
    // RDKit❗✔️:     const auto at_count = rwmol.getNumAtoms();
    // RDKit❗✔️:     const auto bond_count = rwmol.getNumBonds();
    // RDKit❗✔️:
    // RDKit❗✔️:     std::vector<unsigned int> prevAtomOrdering;
    // RDKit❗✔️:     std::vector<unsigned int> prevBondOrdering;
    // RDKit❗✔️:
    // RDKit❗✔️:     rwmol.insertMol(*mol);
    // RDKit❗✔️:
    // RDKit❗✔️:     mol->getProp(RDKit::common_properties::_smilesAtomOutputOrder,
    // RDKit❗✔️:                  prevAtomOrdering);
    // RDKit❗✔️:     mol->getProp(RDKit::common_properties::_smilesBondOutputOrder,
    // RDKit❗✔️:                  prevBondOrdering);
    // RDKit❗✔️:     for (auto i : prevAtomOrdering) {
    // RDKit❗✔️:       atomOrdering.push_back(i + at_count);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (auto i : prevBondOrdering) {
    // RDKit❗✔️:       bondOrdering.push_back(i + bond_count);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   rwmol.setProp(RDKit::common_properties::_smilesAtomOutputOrder, atomOrdering,
    // RDKit❗✔️:                 true);
    // RDKit❗✔️:   rwmol.setProp(RDKit::common_properties::_smilesBondOutputOrder, bondOrdering,
    // RDKit❗✔️:                 true);
    // RDKit❗✔️:
    // RDKit❗✔️:   return getCXExtensions(rwmol, flags);
    // RDKit❗✔️: }
    if selections.len() != templates.len() {
        return Err(SmartsWriteError::CxCoordinateSelectionArity {
            expected: templates.len(),
            actual: selections.len(),
        });
    }
    // Full source preflight before any insertions, and original template order.
    for (template, (query, output)) in templates.iter().enumerate() {
        check_cx_features(query);
        // MolToSmarts returns before recording either output order on zero
        // atoms. Fresh canonical writer evidence must not fabricate presence.
        if !output.source_orders_written
            && (query.prop("_smilesAtomOutputOrder").is_none()
                || query.prop("_smilesBondOutputOrder").is_none())
        {
            return Err(SmartsWriteError::CxMissingOutputOrder { template });
        }
    }
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut atom_order = Vec::new();
    let mut bond_order = Vec::new();
    let mut stereo_groups = Vec::new();
    let mut substance_groups = Vec::new();
    let mut conformer: Option<(usize, bool, Vec<[f64; 3]>)> = None;
    for (template, ((query, output), selection)) in templates.iter().zip(selections).enumerate() {
        let atom_offset = atoms.len();
        let bond_offset = bonds.len();
        let offset_u32 = u32::try_from(atom_offset).map_err(|_| SmartsWriteError::CxRowCount {
            kind: "atom",
            count: atom_offset,
        })?;
        let final_atoms =
            atom_offset
                .checked_add(query.num_atoms())
                .ok_or(SmartsWriteError::CxRowCount {
                    kind: "atom",
                    count: usize::MAX,
                })?;
        let final_bonds =
            bond_offset
                .checked_add(query.num_bonds())
                .ok_or(SmartsWriteError::CxRowCount {
                    kind: "bond",
                    count: usize::MAX,
                })?;
        if final_atoms > u32::MAX as usize {
            return Err(SmartsWriteError::CxRowCount {
                kind: "atom",
                count: final_atoms,
            });
        }
        if final_bonds > u32::MAX as usize {
            return Err(SmartsWriteError::CxRowCount {
                kind: "bond",
                count: final_bonds,
            });
        }
        insert_rows(
            query,
            template,
            atom_offset,
            bond_offset,
            offset_u32,
            &mut atoms,
            &mut bonds,
        )?;
        // SOURCE addAtom resizes all existing conformers with zero rows before
        // insertMol copies the incoming positions. D4 selects one per template.
        if let Some((_, _, points)) = &mut conformer {
            points.resize(final_atoms, [0.0; 3]);
        }
        // D4 selectors are consumed only when CX_COORDS is enabled. Source
        // insertion coordinates are unobservable to its null-conf branches.
        let selected = if coordinates_enabled {
            select_cx_coordinates_from_sets(
                query.conformers_2d(),
                query.conformers_3d(),
                *selection,
            )
            .map_err(|source| SmartsWriteError::CxCoordinateSelection { template, source })?
        } else {
            None
        };
        if let Some(selected) = selected {
            let (id, is_3d) = match selected {
                CoordinateSource::ThreeD(conf) => (conf.id(), conf.is_3d()),
                CoordinateSource::TwoD(conf) => (conf.id(), false),
            };
            let (_, _, points) =
                conformer.get_or_insert_with(|| (id, is_3d, vec![[0.0; 3]; final_atoms]));
            // SOURCE retains the first conformer's is3D; later insertions do
            // not promote it. A false flag still preserves raw Z coordinates.
            for row in 0..query.num_atoms() {
                points[atom_offset + row] = match selected {
                    CoordinateSource::ThreeD(conf) => conf.coordinates()[row],
                    CoordinateSource::TwoD(conf) => {
                        let p = conf.coordinates()[row];
                        [p[0], p[1], 0.0]
                    }
                };
            }
        }
        stereo_groups = insert_stereo_groups(
            &stereo_groups,
            query.stereo_groups(),
            atom_offset,
            bond_offset,
        );
        for group in query_substance_groups(query) {
            let copied = group.clone().with_inserted_offsets(
                SubstanceGroupId::new(substance_groups.len()),
                atom_offset,
                bond_offset,
            );
            substance_groups.push(copied);
        }
        // SOURCE retrieves preexisting order properties after insertMol. The
        // QueryGraph molecule property carrier stores strings; on the source
        // zero-atom early return these must raise the corresponding type error
        // rather than fabricate presence or an empty UInt vector.
        if !output.source_orders_written {
            return Err(SmartsWriteError::CxOutputOrderPropertyType {
                template,
                property: "_smilesAtomOutputOrder",
            });
        }
        for row in &output.atom_order {
            if row.index() >= query.num_atoms() {
                return Err(SmartsWriteError::CxOutputOrder {
                    template,
                    kind: "atom",
                    row: row.index(),
                    count: query.num_atoms(),
                });
            }
            atom_order.push(AtomId::new(atom_offset + row.index()));
        }
        for row in &output.bond_order {
            if row.index() >= query.num_bonds() {
                return Err(SmartsWriteError::CxOutputOrder {
                    template,
                    kind: "bond",
                    row: row.index(),
                    count: query.num_bonds(),
                });
            }
            bond_order.push(BondId::new(bond_offset + row.index()));
        }
    }
    // insertMol never copies molecule-level properties; output-order evidence
    // is explicit. Link-node warnings do not invent an aggregated link record.
    let conformers = conformer
        .into_iter()
        .map(|(id, is_3d, points)| Conformer3D::new(id, points, is_3d))
        .collect();
    let mut graph = QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::new(),
        Vec::new(),
        conformers,
        stereo_groups,
    )?;
    replace_query_substance_groups(&mut graph, substance_groups)?;
    Ok(QueryCxComposition {
        graph,
        atom_order,
        bond_order,
    })
}

fn insert_rows(
    query: &QueryGraph,
    template: usize,
    atom_offset: usize,
    bond_offset: usize,
    offset_u32: u32,
    atoms: &mut Vec<cosmolkit_model::QueryAtom>,
    bonds: &mut Vec<cosmolkit_model::QueryBond>,
) -> Result<(), SmartsWriteError> {
    // RDKit❗✔️: void RWMol::insertMol(const ROMol &other) {
    // RDKit❗✔️:   auto origNumAtoms = getNumAtoms();
    // RDKit❗✔️:   auto origNumBonds = getNumBonds();
    // RDKit❗✔️:   for (const auto oatom : other.atoms()) {
    // RDKit❗✔️:     Atom *newAt = oatom->copy();
    // RDKit❗✔️:     const bool updateLabel = false;
    // RDKit❗✔️:     const bool takeOwnership = true;
    // RDKit❗✔️:     addAtom(newAt, updateLabel, takeOwnership);
    // RDKit❗✔️:     // take care of atom-numbering-dependent properties:
    // RDKit❗✔️:     if (INT_VECT nAtoms;
    // RDKit❗✔️:         newAt->getPropIfPresent(common_properties::_ringStereoAtoms, nAtoms)) {
    // RDKit❗✔️:       for (auto &val : nAtoms) {
    // RDKit❗✔️:         if (val < 0) {
    // RDKit❗✔️:           val = -1 * (-val + origNumAtoms);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           val += origNumAtoms;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       newAt->setProp(common_properties::_ringStereoAtoms, nAtoms, true);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (unsigned int val;
    // RDKit❗✔️:         oatom->getPropIfPresent(common_properties::_ringStereoOtherAtom, val)) {
    // RDKit❗✔️:       newAt->setProp(common_properties::_ringStereoOtherAtom,
    // RDKit❗✔️:                      val + origNumAtoms, true);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   for (const auto obond : other.bonds()) {
    // RDKit❗✔️:     Bond *bond_p = obond->copy();
    // RDKit❗✔️:     unsigned int idx1, idx2;
    // RDKit❗✔️:     idx1 = bond_p->getBeginAtomIdx() + origNumAtoms;
    // RDKit❗✔️:     idx2 = bond_p->getEndAtomIdx() + origNumAtoms;
    // RDKit❗✔️:     bond_p->setOwningMol(this);
    // RDKit❗✔️:     bond_p->setBeginAtomIdx(idx1);
    // RDKit❗✔️:     bond_p->setEndAtomIdx(idx2);
    // RDKit❗✔️:     for (auto &v : bond_p->getStereoAtoms()) {
    // RDKit❗✔️:       v += origNumAtoms;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     const bool takeOwnership = true;
    // RDKit❗✔️:     addBond(bond_p, takeOwnership);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // add atom to any conformers as well, if we have any
    // RDKit❗✔️:   if (other.getNumConformers() && !getNumConformers()) {
    // RDKit❗✔️:     for (const auto &oconf : other.d_confs) {
    // RDKit❗✔️:       auto *nconf = new Conformer(getNumAtoms());
    // RDKit❗✔️:       nconf->set3D(oconf->is3D());
    // RDKit❗✔️:       nconf->setId(oconf->getId());
    // RDKit❗✔️:       for (unsigned int i = 0; i < oconf->getNumAtoms(); ++i) {
    // RDKit❗✔️:         nconf->setAtomPos(i + origNumAtoms, oconf->getAtomPos(i));
    // RDKit❗✔️:       }
    // RDKit❗✔️:       const bool assignId = false;
    // RDKit❗✔️:       addConformer(nconf, assignId);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else if (getNumConformers()) {
    // RDKit❗✔️:     if (other.getNumConformers() == getNumConformers()) {
    // RDKit❗✔️:       ConformerIterator cfi;
    // RDKit❗✔️:       ConstConformerIterator ocfi;
    // RDKit❗✔️:       for (cfi = beginConformers(), ocfi = other.beginConformers();
    // RDKit❗✔️:            cfi != endConformers(); ++cfi, ++ocfi) {
    // RDKit❗✔️:         for (unsigned int i = 0; i < (*ocfi)->getNumAtoms(); ++i) {
    // RDKit❗✔️:           (*cfi)->setAtomPos(i + origNumAtoms, (*ocfi)->getAtomPos(i));
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // add stereo groups
    // RDKit❗✔️:   insertStereoGroups(*this, other, origNumAtoms, origNumBonds);
    // RDKit❗✔️:   // add substance groups
    // RDKit❗✔️:   insertSubstanceGroups(*this, other, origNumAtoms, origNumBonds);
    // RDKit❗✔️: }
    // One copy of the source query carrier/predicate, as RWMol::insertMol.
    // Explicit signed/unsigned wrap reproduces its ring-reference arithmetic.
    for source in query.atoms() {
        let mut atom = source
            .clone()
            .with_id(AtomId::new(atom_offset + source.id().index()));
        if let Some(value) = atom.prop("_ringStereoAtoms") {
            let PropertyValue::IntVector(values) = value else {
                return Err(SmartsWriteError::CxAtomPropertyKind {
                    template,
                    atom: source.id(),
                    property: "_ringStereoAtoms",
                    kind: value.kind(),
                });
            };
            let values: Vec<i32> = values
                .iter()
                .map(|value| {
                    if *value < 0 {
                        value
                            .wrapping_neg()
                            .wrapping_add(offset_u32 as i32)
                            .wrapping_neg()
                    } else {
                        value.wrapping_add(offset_u32 as i32)
                    }
                })
                .collect();
            atom.set_computed_prop("_ringStereoAtoms", PropertyValue::IntVector(values))
                .expect("nonempty source property key");
        }
        if let Some(value) = source.prop("_ringStereoOtherAtom") {
            let value = cosmolkit_core::property_value_to_uint(value).map_err(|source_error| {
                SmartsWriteError::CxAtomPropertyUInt {
                    template,
                    atom: source.id(),
                    property: "_ringStereoOtherAtom",
                    source: source_error,
                }
            })?;
            atom.set_computed_prop("_ringStereoOtherAtom", value.wrapping_add(offset_u32))
                .expect("nonempty source property key");
        }
        atoms.push(atom);
    }
    for source in query.bonds() {
        let carrier = source.bond();
        let references = carrier
            .stereo_atoms()
            .map(|pair| pair.map(|row| AtomId::new(row.index() + atom_offset)));
        let bond = source.clone().remapped(
            BondId::new(bond_offset + source.id().index()),
            AtomId::new(carrier.begin().index() + atom_offset),
            AtomId::new(carrier.end().index() + atom_offset),
            references,
        );
        bonds.push(bond);
    }
    Ok(())
}
