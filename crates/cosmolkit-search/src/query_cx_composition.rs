//! Source CX vector insertion over canonical query graphs, retaining every AST.
use crate::{SmartsWriteError, SmartsWriteOutput};
use cosmolkit_model::{
    AtomId, BondId, Conformer3D, PropertyValue, QueryAtom, QueryBond, QueryGraph, StereoGroup,
    SubstanceGroup, SubstanceGroupId, insert_stereo_groups, query_substance_groups,
    replace_query_substance_groups,
};
use cosmolkit_smiles::{CoordinateSource, CxSmilesFields};
use std::collections::BTreeMap;

/// The same query graph type and writer evidence, with no second AST, live
/// molecule, cache authority or compatibility topology.
#[doc(hidden)]
pub struct QueryCxComposition {
    pub graph: QueryGraph,
    pub extension: cosmolkit_model::PropertyText,
    pub source_orders_written: bool,
    pub atom_order: Vec<AtomId>,
    pub bond_order: Vec<BondId>,
}

// Private detached output assembly, using the existing canonical graph rows.
// No second AST, adjacency/cache, live molecule or commit authority.
#[derive(Default)]
struct QueryAssembly {
    atoms: Vec<QueryAtom>,
    bonds: Vec<QueryBond>,
    conformers: Vec<Conformer3D>,
    stereo_groups: Vec<StereoGroup>,
    substance_groups: Vec<SubstanceGroup>,
}

fn check_cx_features(
    query: &QueryGraph,
    warning: &mut dyn FnMut(&'static str),
) -> Result<(), SmartsWriteError> {
    // RDKit❗❌: void checkCXFeatures(const ROMol &mol) {
    // RDKit❗❌:   std::string lns;
    // RDKit❗❌:   if (mol.getPropIfPresent(common_properties::molFileLinkNodes, lns)) {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "CX Extensions: mol has link nodes which are not currently supported"
    // RDKit❗❌:         << std::endl;
    // RDKit❗❌:   }
    // RDKit❗❌:   const auto &sgs = getSubstanceGroups(mol);
    // RDKit❗❌:   auto parent_check =
    // RDKit❗❌:       std::any_of(sgs.cbegin(), sgs.cend(), [&](const SubstanceGroup &sg) {
    // RDKit❗❌:         if (sg.hasProp("PARENT")) {
    // RDKit❗❌:           return true;
    // RDKit❗❌:         }
    // RDKit❗❌:         return false;
    // RDKit❗❌:       });
    // RDKit❗❌:   if (parent_check) {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "CX Extensions: Substance group hierarchy is not always preserved."
    // RDKit❗❌:         << std::endl;
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Execute the source string read even though only property presence is
    // consumed. Delegate every modeled tag conversion to its sole CORE owner.
    // Any conversion error precedes both the link-node and hierarchy warnings.
    if let Some(value) = query.prop("_molLinkNodes") {
        let _lns = cosmolkit_core::property_value_to_string(value)?;
        warning("CX Extensions: mol has link nodes which are not currently supported");
    }
    // Typed parent is the model projection of the same source PARENT fact;
    // raw property presence also counts, without reading a numeric parent.
    let parent_check = query_substance_groups(query)
        .iter()
        .any(|group| group.parent().is_some() || group.props().contains_key(b"PARENT".as_slice()));
    if parent_check {
        warning("CX Extensions: Substance group hierarchy is not always preserved.");
    }
    // Private warning sink corresponds to source rdWarningLog; production
    // emits each warning at once, and tests observe the same ordered effects.
    // Cost ❌: PropertyText lacks Native short-string storage; property lookup
    // uses modeled tree-backed stores. Groups remain borrowed and short-circuit.
    Ok(())
}

#[doc(hidden)]
pub fn compose_query_cx_templates(
    templates: &[(&QueryGraph, &SmartsWriteOutput)],
    flags: CxSmilesFields,
) -> Result<QueryCxComposition, SmartsWriteError> {
    // RDKit❗❌: std::string getCXExtensions(const std::vector<ROMol *> &mols,
    // RDKit❗❌:                             std::uint32_t flags) {
    // RDKit❗❌:   for (const auto &mol : mols) {
    // RDKit❗❌:     checkCXFeatures(*mol);
    // RDKit❗❌:     if (!mol->hasProp(RDKit::common_properties::_smilesAtomOutputOrder) ||
    // RDKit❗❌:         !mol->hasProp(RDKit::common_properties::_smilesBondOutputOrder)) {
    // RDKit❗❌:       throw ValueErrorException(
    // RDKit❗❌:           "Input molecule does not have the required "
    // RDKit❗❌:           "smiles ordering properties set");
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   RDKit::RWMol rwmol;
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<unsigned int> atomOrdering;
    // RDKit❗❌:   std::vector<unsigned int> bondOrdering;
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto &mol : mols) {
    // RDKit❗❌:     const auto at_count = rwmol.getNumAtoms();
    // RDKit❗❌:     const auto bond_count = rwmol.getNumBonds();
    // RDKit❗❌:
    // RDKit❗❌:     std::vector<unsigned int> prevAtomOrdering;
    // RDKit❗❌:     std::vector<unsigned int> prevBondOrdering;
    // RDKit❗❌:
    // RDKit❗❌:     rwmol.insertMol(*mol);
    // RDKit❗❌:
    // RDKit❗❌:     mol->getProp(RDKit::common_properties::_smilesAtomOutputOrder,
    // RDKit❗❌:                  prevAtomOrdering);
    // RDKit❗❌:     mol->getProp(RDKit::common_properties::_smilesBondOutputOrder,
    // RDKit❗❌:                  prevBondOrdering);
    // RDKit❗❌:     for (auto i : prevAtomOrdering) {
    // RDKit❗❌:       atomOrdering.push_back(i + at_count);
    // RDKit❗❌:     }
    // RDKit❗❌:     for (auto i : prevBondOrdering) {
    // RDKit❗❌:       bondOrdering.push_back(i + bond_count);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   rwmol.setProp(RDKit::common_properties::_smilesAtomOutputOrder, atomOrdering,
    // RDKit❗❌:                 true);
    // RDKit❗❌:   rwmol.setProp(RDKit::common_properties::_smilesBondOutputOrder, bondOrdering,
    // RDKit❗❌:                 true);
    // RDKit❗❌:
    // RDKit❗❌:   return getCXExtensions(rwmol, flags);
    // RDKit❗❌: }
    // Complete source preflight occurs before any insertion/conformer read.
    // Source writer outputs are the explicit detached typed getter values;
    // false does not fabricate presence for a zero-atom writer early return.
    for (template, (query, output)) in templates.iter().enumerate() {
        check_cx_features(query, &mut |message| eprintln!("{message}"))?;
        if !output.source_orders_written
            && (query.prop("_smilesAtomOutputOrder").is_none()
                || query.prop("_smilesBondOutputOrder").is_none())
        {
            return Err(SmartsWriteError::CxMissingOutputOrder { template });
        }
    }
    let mut assembly = QueryAssembly::default();
    let mut atom_order = Vec::new();
    let mut bond_order = Vec::new();
    for (template, (query, output)) in templates.iter().enumerate() {
        let atom_offset = assembly.atoms.len() as u32;
        let bond_offset = assembly.bonds.len() as u32;
        insert_query_mol(query, template, &mut assembly, || {
            // Native insertMol iterates ALL actual conformers irrespective of
            // CX_COORDS. Source iterator errors stay at this coordinate stage,
            // after source atom/bond/property insertion prefixes.
            query
                .source_conformers()
                .map(|rows| {
                    rows.into_iter()
                        .map(|row| match row {
                            cosmolkit_model::CoordinateSourceConformer::TwoD(c) => {
                                CoordinateSource::TwoD(c)
                            }
                            cosmolkit_model::CoordinateSourceConformer::ThreeD(c) => {
                                CoordinateSource::ThreeD(c)
                            }
                        })
                        .collect()
                })
                .map_err(SmartsWriteError::CxCoordinateSource)
        })?;
        if !output.source_orders_written {
            // Generic modeled RDValue tags cannot masquerade as UInt vectors.
            // The first source read is atom order; insertion already happened.
            return Err(SmartsWriteError::CxOutputOrderPropertyType {
                template,
                property: "_smilesAtomOutputOrder",
            });
        }
        for row in &output.atom_order {
            let row = u32::try_from(row.index()).map_err(|_| SmartsWriteError::CxRowCount {
                kind: "atom output order",
                count: row.index(),
            })?;
            atom_order.push(AtomId::new(row.wrapping_add(atom_offset) as usize));
        }
        for row in &output.bond_order {
            let row = u32::try_from(row.index()).map_err(|_| SmartsWriteError::CxRowCount {
                kind: "bond output order",
                count: row.index(),
            })?;
            bond_order.push(BondId::new(row.wrapping_add(bond_offset) as usize));
        }
    }
    let mut graph = QueryGraph::from_parts(
        assembly.atoms,
        assembly.bonds,
        [],
        vec![],
        vec![],
        assembly.stereo_groups,
    )?;
    // Native addConformer checks row count, never finiteness or duplicate IDs.
    // The existing canonical source append boundary records exact order.
    for conformer in assembly.conformers {
        graph.add_conformer_3d(conformer)?;
    }
    replace_query_substance_groups(&mut graph, assembly.substance_groups)?;
    // The two UInt-vector values are retained in the typed output fields,
    // rather than inserting wrong-tag IntVector/String surrogates. Actual
    // computed-name dictionary effects use the canonical MODEL setter owner.
    let mut props = graph.source_molecule_properties();
    props.register_transient_computed_name("_smilesAtomOutputOrder")?;
    props.register_transient_computed_name("_smilesBondOutputOrder")?;
    graph.replace_source_molecule_properties(&props);
    let extension = crate::smarts_write::write_query_cx_extensions(
        &mut graph,
        &atom_order,
        &bond_order,
        flags,
    )?;
    // Known extra copied scalar/property/output views retain second❌. Source
    // insertion/order/flag loops run once with no atom-order permutation guess,
    // per-template row filtering, map normalization or conformer selection.
    Ok(QueryCxComposition {
        graph,
        extension,
        atom_order,
        bond_order,
        source_orders_written: true,
    })
}

fn insert_query_mol<'a>(
    query: &'a QueryGraph,
    template: usize,
    assembly: &mut QueryAssembly,
    conformer_reader: impl FnOnce() -> Result<Vec<CoordinateSource<'a>>, SmartsWriteError>,
) -> Result<(), SmartsWriteError> {
    // BEGIN RDKIT CPP FUNCTION RWMol::insertMol
    // RDKit❗❌: void RWMol::insertMol(const ROMol &other) {
    // RDKit❗❌:   auto origNumAtoms = getNumAtoms();
    // RDKit❗❌:   auto origNumBonds = getNumBonds();
    // RDKit❗❌:   for (const auto oatom : other.atoms()) {
    // RDKit❗❌:     Atom *newAt = oatom->copy();
    // RDKit❗❌:     const bool updateLabel = false;
    // RDKit❗❌:     const bool takeOwnership = true;
    // RDKit❗❌:     addAtom(newAt, updateLabel, takeOwnership);
    // RDKit❗❌:     // take care of atom-numbering-dependent properties:
    // RDKit❗❌:     if (INT_VECT nAtoms;
    // RDKit❗❌:         newAt->getPropIfPresent(common_properties::_ringStereoAtoms, nAtoms)) {
    // RDKit❗❌:       for (auto &val : nAtoms) {
    // RDKit❗❌:         if (val < 0) {
    // RDKit❗❌:           val = -1 * (-val + origNumAtoms);
    // RDKit❗❌:         } else {
    // RDKit❗❌:           val += origNumAtoms;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       newAt->setProp(common_properties::_ringStereoAtoms, nAtoms, true);
    // RDKit❗❌:     }
    // RDKit❗❌:     if (unsigned int val;
    // RDKit❗❌:         oatom->getPropIfPresent(common_properties::_ringStereoOtherAtom, val)) {
    // RDKit❗❌:       newAt->setProp(common_properties::_ringStereoOtherAtom,
    // RDKit❗❌:                      val + origNumAtoms, true);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto obond : other.bonds()) {
    // RDKit❗❌:     Bond *bond_p = obond->copy();
    // RDKit❗❌:     unsigned int idx1, idx2;
    // RDKit❗❌:     idx1 = bond_p->getBeginAtomIdx() + origNumAtoms;
    // RDKit❗❌:     idx2 = bond_p->getEndAtomIdx() + origNumAtoms;
    // RDKit❗❌:     bond_p->setOwningMol(this);
    // RDKit❗❌:     bond_p->setBeginAtomIdx(idx1);
    // RDKit❗❌:     bond_p->setEndAtomIdx(idx2);
    // RDKit❗❌:     for (auto &v : bond_p->getStereoAtoms()) {
    // RDKit❗❌:       v += origNumAtoms;
    // RDKit❗❌:     }
    // RDKit❗❌:     const bool takeOwnership = true;
    // RDKit❗❌:     addBond(bond_p, takeOwnership);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // add atom to any conformers as well, if we have any
    // RDKit❗❌:   if (other.getNumConformers() && !getNumConformers()) {
    // RDKit❗❌:     for (const auto &oconf : other.d_confs) {
    // RDKit❗❌:       auto *nconf = new Conformer(getNumAtoms());
    // RDKit❗❌:       nconf->set3D(oconf->is3D());
    // RDKit❗❌:       nconf->setId(oconf->getId());
    // RDKit❗❌:       for (unsigned int i = 0; i < oconf->getNumAtoms(); ++i) {
    // RDKit❗❌:         nconf->setAtomPos(i + origNumAtoms, oconf->getAtomPos(i));
    // RDKit❗❌:       }
    // RDKit❗❌:       const bool assignId = false;
    // RDKit❗❌:       addConformer(nconf, assignId);
    // RDKit❗❌:     }
    // RDKit❗❌:   } else if (getNumConformers()) {
    // RDKit❗❌:     if (other.getNumConformers() == getNumConformers()) {
    // RDKit❗❌:       ConformerIterator cfi;
    // RDKit❗❌:       ConstConformerIterator ocfi;
    // RDKit❗❌:       for (cfi = beginConformers(), ocfi = other.beginConformers();
    // RDKit❗❌:            cfi != endConformers(); ++cfi, ++ocfi) {
    // RDKit❗❌:         for (unsigned int i = 0; i < (*ocfi)->getNumAtoms(); ++i) {
    // RDKit❗❌:           (*cfi)->setAtomPos(i + origNumAtoms, (*ocfi)->getAtomPos(i));
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // add stereo groups
    // RDKit❗❌:   insertStereoGroups(*this, other, origNumAtoms, origNumBonds);
    // RDKit❗❌:   // add substance groups
    // RDKit❗❌:   insertSubstanceGroups(*this, other, origNumAtoms, origNumBonds);
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RWMol::insertMol
    // The complete insertion has one owner. Atom append and existing-conformer
    // zero writes precede the ring-property reads/setters, including errors.
    // Computed setter errors propagate structurally, not as an expect panic.
    // Source three conformer-count branches use physical positions, never IDs,
    // and only newly created conformers copy incoming ID/is3D (not properties).
    // The actual source caller supplies every conformer in physical append
    // order; no D4 projection or CX flag changes this insertion behavior.
    // Boundary gaps: native pointer ownership and unsigned32 index widths,
    // signed-negation UB at INT_MIN, model two-reference stereo cardinality,
    // final detached group/reference validation versus native setter timing.
    // Cost: linear rows, reached properties, coordinates and groups; coordinate
    // view Vec allocation and detached group-role/index reconstruction add
    // source-absent costs, so overall complexity marker remains a known loss.
    let atom_offset = assembly.atoms.len();
    let bond_offset = assembly.bonds.len();
    let offset_u32 = u32::try_from(atom_offset).map_err(|_| SmartsWriteError::CxRowCount {
        kind: "atom",
        count: atom_offset,
    })?;
    for (kind, old, added) in [
        ("atom", atom_offset, query.num_atoms()),
        ("bond", bond_offset, query.num_bonds()),
    ] {
        let count = old.checked_add(added).ok_or(SmartsWriteError::CxRowCount {
            kind,
            count: usize::MAX,
        })?;
        if count > u32::MAX as usize {
            return Err(SmartsWriteError::CxRowCount { kind, count });
        }
    }
    for source in query.atoms() {
        let atom_id = AtomId::new(atom_offset + source.id().index());
        assembly.atoms.push(source.clone().with_id(atom_id));
        for conformer in &mut assembly.conformers {
            conformer
                .source_set_atom_position(atom_id.index(), [0.0; 3])
                .map_err(|source| SmartsWriteError::CxCoordinateStorage { template, source })?;
        }
        let atom = assembly
            .atoms
            .last_mut()
            .expect("just appended source atom");
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
                .map_err(|source_error| SmartsWriteError::CxAtomPropertyWrite {
                    template,
                    atom: source.id(),
                    property: "_ringStereoAtoms",
                    source: source_error,
                })?;
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
                .map_err(|source_error| SmartsWriteError::CxAtomPropertyWrite {
                    template,
                    atom: source.id(),
                    property: "_ringStereoOtherAtom",
                    source: source_error,
                })?;
        }
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
        assembly.bonds.push(bond);
    }

    let incoming = conformer_reader()?;
    if !incoming.is_empty() && assembly.conformers.is_empty() {
        for source in &incoming {
            let (id, is_3d) = match source {
                CoordinateSource::TwoD(conf) => (conf.id(), false),
                CoordinateSource::ThreeD(conf) => (conf.id(), conf.is_3d()),
            };
            let mut conformer = Conformer3D::new(id, vec![[0.0; 3]; assembly.atoms.len()], is_3d);
            copy_source_conformer_positions(&mut conformer, *source, atom_offset, template)?;
            if conformer.coordinates().len() != assembly.atoms.len() {
                return Err(SmartsWriteError::CxCoordinateStorage {
                    template,
                    source: cosmolkit_model::CoordinateValidationError::RowCount {
                        dimension: "3D",
                        conformer: conformer.id(),
                        rows: conformer.coordinates().len(),
                        atom_count: assembly.atoms.len(),
                    },
                });
            }
            assembly.conformers.push(conformer);
        }
    } else if !assembly.conformers.is_empty() && incoming.len() == assembly.conformers.len() {
        for (target, source) in assembly.conformers.iter_mut().zip(&incoming) {
            copy_source_conformer_positions(target, *source, atom_offset, template)?;
        }
    }
    assembly.stereo_groups = insert_stereo_groups(
        &assembly.stereo_groups,
        query.stereo_groups(),
        atom_offset,
        bond_offset,
    )
    .map_err(SmartsWriteError::StereoGroup)?;
    for group in query_substance_groups(query) {
        let copy = group.clone().with_inserted_offsets(
            SubstanceGroupId::new(assembly.substance_groups.len()),
            atom_offset,
            bond_offset,
        );
        assembly.substance_groups.push(copy);
    }
    Ok(())
}

fn copy_source_conformer_positions(
    target: &mut Conformer3D,
    source: CoordinateSource<'_>,
    offset: usize,
    template: usize,
) -> Result<(), SmartsWriteError> {
    // RDKit❗✔️:       for (unsigned int i = 0; i < oconf->getNumAtoms(); ++i) {
    // RDKit❗✔️:         nconf->setAtomPos(i + origNumAtoms, oconf->getAtomPos(i));
    // RDKit❗✔️:       }
    // Both native coordinate branches call the same get/set pair. Borrowed
    // raw coordinate rows are copied once, no flags/ID/property conversion.
    // TwoD is the existing explicit planar projection with source Z=0.
    match source {
        CoordinateSource::ThreeD(conformer) => {
            for (row, &point) in conformer.coordinates().iter().enumerate() {
                target
                    .source_set_atom_position(offset + row, point)
                    .map_err(|source| SmartsWriteError::CxCoordinateStorage { template, source })?;
            }
        }
        CoordinateSource::TwoD(conformer) => {
            for (row, &[x, y]) in conformer.coordinates().iter().enumerate() {
                target
                    .source_set_atom_position(offset + row, [x, y, 0.0])
                    .map_err(|source| SmartsWriteError::CxCoordinateStorage { template, source })?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod source_insert_mol_complete_tests {
    use super::*;
    use cosmolkit_model::{
        AtomSpec, Bond, BondQueryPredicate, BondSpec, Conformer2D, QueryNode, StereoGroupKind,
        SubstanceGroupKind,
    };
    use cosmolkit_types::{BondOrder, BondStereo, Element};

    fn graph(count: usize) -> QueryGraph {
        QueryGraph::from_parts(
            (0..count)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            vec![],
            BTreeMap::new(),
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn initialized() -> QueryAssembly {
        let mut assembly = QueryAssembly::default();
        insert_query_mol(&graph(1), 0, &mut assembly, || Ok(Vec::new())).unwrap();
        assembly
    }

    #[test]
    fn source_new_conformers_copy_all_ids_flags_and_positions_but_not_input_properties() {
        let mut assembly = initialized();
        let query = graph(1);
        let first = Conformer3D::new(7, vec![[1.0, 2.0, 3.0]], false).with_prop("input", "ignored");
        let second = Conformer3D::new(7, vec![[4.0, 5.0, 6.0]], true).with_prop("input", "ignored");
        insert_query_mol(&query, 1, &mut assembly, || {
            Ok(vec![
                CoordinateSource::ThreeD(&first),
                CoordinateSource::ThreeD(&second),
            ])
        })
        .unwrap();
        assert_eq!(assembly.conformers.len(), 2);
        assert_eq!(assembly.conformers[0].id(), 7);
        assert_eq!(assembly.conformers[1].id(), 7);
        assert!(!assembly.conformers[0].is_3d());
        assert!(assembly.conformers[1].is_3d());
        assert_eq!(
            assembly.conformers[0].coordinates(),
            [[0.0; 3], [1.0, 2.0, 3.0]]
        );
        assert_eq!(
            assembly.conformers[1].coordinates(),
            [[0.0; 3], [4.0, 5.0, 6.0]]
        );
        assert!(assembly.conformers[0].props().is_empty());
        assert!(assembly.conformers[1].props().is_empty());
        assert_eq!(first.coordinates(), [[1.0, 2.0, 3.0]]);
        assert_eq!(first.props().len(), 1);
    }

    #[test]
    fn source_equal_counts_pair_by_append_position_and_keep_existing_ids_flags_properties() {
        let mut assembly = initialized();
        let query = graph(1);
        assembly.conformers = vec![
            Conformer3D::new(9, vec![[1.0, 2.0, 3.0]], false).with_prop("target", "first"),
            Conformer3D::new(8, vec![[4.0, 5.0, 6.0]], true).with_prop("target", "second"),
        ];
        let first =
            Conformer3D::new(8, vec![[7.0, 8.0, 9.0]], true).with_prop("input", "not copied");
        let second = Conformer2D::new(9, vec![[10.0, 11.0]]);
        let before: Vec<_> = assembly
            .conformers
            .iter()
            .map(|conf| (conf.id(), conf.is_3d(), conf.props().clone()))
            .collect();
        insert_query_mol(&query, 1, &mut assembly, || {
            Ok(vec![
                CoordinateSource::ThreeD(&first),
                CoordinateSource::TwoD(&second),
            ])
        })
        .unwrap();
        assert_eq!(
            assembly.conformers[0].coordinates(),
            [[1.0, 2.0, 3.0], [7.0, 8.0, 9.0]]
        );
        assert_eq!(
            assembly.conformers[1].coordinates(),
            [[4.0, 5.0, 6.0], [10.0, 11.0, 0.0]]
        );
        for (conf, (id, is_3d, props)) in assembly.conformers.iter().zip(before) {
            assert_eq!(conf.id(), id);
            assert_eq!(conf.is_3d(), is_3d);
            assert_eq!(conf.props(), &props);
        }
    }

    #[test]
    fn source_different_or_zero_counts_keep_zero_extensions_without_choosing_an_input_frame() {
        for incoming_count in [0, 1, 3] {
            let mut assembly = initialized();
            assembly.conformers = vec![
                Conformer3D::new(7, vec![[1.0, 2.0, 3.0]], false),
                Conformer3D::new(9, vec![[4.0, 5.0, 6.0]], true),
            ];
            let query = graph(1);
            let incoming = Conformer3D::new(99, vec![[100.0, 101.0, 102.0]], true);
            insert_query_mol(&query, 1, &mut assembly, || {
                Ok(vec![CoordinateSource::ThreeD(&incoming); incoming_count])
            })
            .unwrap();
            assert_eq!(assembly.conformers.len(), 2);
            assert_eq!(
                assembly.conformers[0].coordinates(),
                [[1.0, 2.0, 3.0], [0.0; 3]]
            );
            assert_eq!(
                assembly.conformers[1].coordinates(),
                [[4.0, 5.0, 6.0], [0.0; 3]]
            );
        }
        let mut assembly = initialized();
        insert_query_mol(&graph(1), 1, &mut assembly, || Ok(Vec::new())).unwrap();
        assert!(assembly.conformers.is_empty());
    }

    #[test]
    fn source_ring_reference_writes_wrap_unsigned_and_shift_stereo_and_group_rows() {
        let mut query = graph(4);
        query
            .atom_mut(0)
            .unwrap()
            .set_prop(
                "_ringStereoAtoms",
                PropertyValue::IntVector(vec![-3, -1, 0, 2]),
            )
            .unwrap();
        query
            .atom_mut(0)
            .unwrap()
            .set_prop("_ringStereoOtherAtom", PropertyValue::UInt(u32::MAX))
            .unwrap();
        let mut bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double)
                .with_stereo(BondStereo::E)
                .with_stereo_atoms(AtomId::new(0), AtomId::new(3)),
        );
        bond.set_prop("raw", PropertyValue::UInt(u32::MAX)).unwrap();
        let mut query = QueryGraph::from_parts(
            query.atoms().to_vec(),
            vec![QueryBond::from_carrier_parts(
                bond,
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
            )],
            BTreeMap::new(),
            vec![],
            vec![],
            vec![
                StereoGroup::new(
                    StereoGroupKind::And,
                    vec![AtomId::new(1)],
                    vec![BondId::new(0)],
                )
                .expect("valid distinct stereo members")
                .with_id(7)
                .with_write_id(9),
            ],
        )
        .unwrap();
        replace_query_substance_groups(
            &mut query,
            vec![
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                    .with_atoms(vec![AtomId::new(0), AtomId::new(3)]),
            ],
        )
        .unwrap();
        let before = query.clone();
        let mut assembly = QueryAssembly::default();
        insert_query_mol(&graph(2), 0, &mut assembly, || Ok(Vec::new())).unwrap();
        insert_query_mol(&query, 1, &mut assembly, || Ok(Vec::new())).unwrap();
        assert_eq!(
            assembly.atoms[2].prop("_ringStereoAtoms"),
            Some(&PropertyValue::IntVector(vec![-5, -3, 2, 4]))
        );
        assert_eq!(
            assembly.atoms[2].prop("_ringStereoOtherAtom"),
            Some(&PropertyValue::UInt(1))
        );
        assert!(
            assembly.atoms[2]
                .is_prop_computed("_ringStereoAtoms")
                .unwrap()
        );
        assert!(
            assembly.atoms[2]
                .is_prop_computed("_ringStereoOtherAtom")
                .unwrap()
        );
        assert_eq!(assembly.bonds[0].begin(), AtomId::new(3));
        assert_eq!(assembly.bonds[0].end(), AtomId::new(4));
        assert_eq!(
            assembly.bonds[0].bond().stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(5)])
        );
        assert_eq!(
            assembly.bonds[0].bond().prop("raw"),
            Some(&PropertyValue::UInt(u32::MAX))
        );
        assert_eq!(assembly.stereo_groups[0].atoms(), [AtomId::new(3)]);
        assert_eq!(assembly.stereo_groups[0].id(), Some(7));
        assert_eq!(assembly.stereo_groups[0].write_id(), 0);
        assert_eq!(
            assembly.substance_groups[0].atoms(),
            [AtomId::new(2), AtomId::new(5)]
        );
        assert_eq!(query, before);
    }

    #[test]
    fn source_property_failure_retains_appended_rows_and_prior_zero_writes_before_coordinate_reader()
     {
        let mut query = graph(2);
        query
            .atom_mut(1)
            .unwrap()
            .set_prop("_ringStereoAtoms", PropertyValue::UInt(3))
            .unwrap();
        let mut assembly = initialized();
        assembly
            .conformers
            .push(Conformer3D::new(7, vec![[1.0, 2.0, 3.0]], true));
        let result = insert_query_mol(&query, 1, &mut assembly, || {
            panic!("native row property error precedes coordinate stage")
        });
        assert!(
            matches!(result,Err(SmartsWriteError::CxAtomPropertyKind {atom,property:"_ringStereoAtoms",..}) if atom==AtomId::new(1))
        );
        assert_eq!(assembly.atoms.len(), 3);
        assert_eq!(
            assembly.conformers[0].coordinates(),
            [[1.0, 2.0, 3.0], [0.0; 3], [0.0; 3]]
        );
        assert!(assembly.bonds.is_empty());
        assert!(assembly.stereo_groups.is_empty());
        assert!(assembly.substance_groups.is_empty());
        let mut query = graph(1);
        query
            .atom_mut(0)
            .unwrap()
            .set_prop("_ringStereoAtoms", PropertyValue::IntVector(vec![0]))
            .unwrap();
        query
            .atom_mut(0)
            .unwrap()
            .set_prop("__computedProps", PropertyValue::UInt(1))
            .unwrap();
        let mut assembly = initialized();
        let result = insert_query_mol(&query, 1, &mut assembly, || {
            panic!("computed-setter conversion error must propagate before coordinate stage")
        });
        assert!(matches!(
            result,
            Err(SmartsWriteError::CxAtomPropertyWrite {
                property: "_ringStereoAtoms",
                ..
            })
        ));
        assert_eq!(assembly.atoms.len(), 2);
        assert_eq!(
            assembly.atoms[1].prop("_ringStereoAtoms"),
            Some(&PropertyValue::IntVector(vec![0]))
        );
    }
}

#[cfg(test)]
mod check_cx_features_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyText, QueryAtom, SubstanceGroupKind};
    use cosmolkit_types::Element;
    const LINK: &str = "CX Extensions: mol has link nodes which are not currently supported";
    const PARENT: &str = "CX Extensions: Substance group hierarchy is not always preserved.";
    fn graph() -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn group(i: usize) -> SubstanceGroup {
        SubstanceGroup::new(SubstanceGroupId::new(i), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(0)])
    }
    fn warnings(q: &QueryGraph) -> Vec<&'static str> {
        let mut output = vec![];
        check_cx_features(q, &mut |s| output.push(s)).unwrap();
        output
    }
    #[test]
    fn missing_properties_emit_no_warning_and_keep_graph_unchanged() {
        let q = graph();
        let before = q.clone();
        assert!(warnings(&q).is_empty());
        assert_eq!(q, before);
    }
    #[test]
    fn empty_present_link_node_string_still_warns() {
        let mut q = graph();
        q.set_prop("_molLinkNodes", PropertyValue::String("".into()))
            .unwrap();
        assert_eq!(warnings(&q), vec![LINK]);
    }
    #[test]
    fn every_modeled_link_node_tag_uses_source_string_projection() {
        for value in [
            PropertyValue::String("x".into()),
            PropertyValue::Int(-3),
            PropertyValue::UInt(u32::MAX),
            PropertyValue::Double(1.25),
            PropertyValue::Bool(false),
            PropertyValue::IntVector(vec![1, -2]),
            PropertyValue::StringVector(vec!["".into(), "x".into()]),
        ] {
            let mut q = graph();
            q.set_prop("_molLinkNodes", value).unwrap();
            assert_eq!(warnings(&q), vec![LINK]);
        }
    }
    #[test]
    fn parent_warning_tests_presence_without_numeric_conversion() {
        let mut q = graph();
        let g = group(0)
            .with_prop("PARENT", PropertyValue::String("not-a-number".into()))
            .unwrap();
        replace_query_substance_groups(&mut q, vec![g]).unwrap();
        assert_eq!(warnings(&q), vec![PARENT]);
    }
    #[test]
    fn typed_parent_projection_is_the_same_source_presence_fact() {
        let mut q = graph();
        replace_query_substance_groups(
            &mut q,
            vec![group(0), group(1).with_parent(SubstanceGroupId::new(0))],
        )
        .unwrap();
        assert_eq!(warnings(&q), vec![PARENT]);
    }
    #[test]
    fn several_parent_properties_emit_one_warning() {
        let mut q = graph();
        replace_query_substance_groups(
            &mut q,
            vec![
                group(0)
                    .with_prop("PARENT", PropertyValue::Bool(false))
                    .unwrap(),
                group(1).with_prop("PARENT", PropertyValue::Int(0)).unwrap(),
            ],
        )
        .unwrap();
        assert_eq!(warnings(&q), vec![PARENT]);
    }
    #[test]
    fn link_warning_precedes_hierarchy_and_preserves_binary_payload() {
        let mut q = graph();
        let mut text = PropertyText::new();
        text.extend_bytes(&[255, 0, b'.']);
        q.set_prop("_molLinkNodes", PropertyValue::String(text))
            .unwrap();
        replace_query_substance_groups(
            &mut q,
            vec![
                group(0)
                    .with_prop("PARENT", PropertyValue::String("".into()))
                    .unwrap(),
            ],
        )
        .unwrap();
        let before = q.clone();
        assert_eq!(warnings(&q), vec![LINK, PARENT]);
        assert_eq!(q, before);
    }
}
#[cfg(test)]
mod vector_cx_composition_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, Element, PropertyText};
    use cosmolkit_smiles::CxSmilesFields as F;
    use cosmolkit_types::BondOrder;
    fn graph(n: usize) -> QueryGraph {
        QueryGraph::from_parts(
            (0..n)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn output(q: &QueryGraph) -> SmartsWriteOutput {
        SmartsWriteOutput {
            text: PropertyText::default(),
            atom_order: (0..q.num_atoms()).map(AtomId::new).collect(),
            bond_order: (0..q.num_bonds()).map(BondId::new).collect(),
            source_orders_written: true,
        }
    }
    fn run(q: &QueryGraph) -> QueryCxComposition {
        compose_query_cx_templates(&[(q, &output(q))], F::NONE).unwrap()
    }
    #[test]
    fn empty_templates_record_both_computed_names_and_reach_final_writer() {
        let result = compose_query_cx_templates(&[], F::NONE).unwrap();
        assert!(result.source_orders_written);
        assert!(result.extension.is_empty());
        assert_eq!(result.graph.num_atoms(), 0);
        assert_eq!(
            result.graph.computed_prop_names().unwrap().unwrap(),
            [
                PropertyText::from("_smilesAtomOutputOrder"),
                PropertyText::from("_smilesBondOutputOrder")
            ]
        );
    }
    #[test]
    fn forward_offsets_keep_repetitions_and_source_graphs_unchanged() {
        let a = graph(2);
        let b = graph(3);
        let before = (a.clone(), b.clone());
        let mut ao = output(&a);
        ao.atom_order = vec![AtomId::new(1), AtomId::new(1), AtomId::new(0)];
        let mut bo = output(&b);
        bo.atom_order = vec![AtomId::new(2), AtomId::new(0)];
        let r = compose_query_cx_templates(&[(&a, &ao), (&b, &bo)], F::NONE).unwrap();
        assert_eq!(r.atom_order, [1, 1, 0, 4, 2].map(AtomId::new));
        assert_eq!((a, b), before);
    }
    #[test]
    fn all_template_presence_checks_precede_first_insertion_property_read() {
        let mut a = graph(1);
        a.atom_mut(0)
            .unwrap()
            .set_prop("_ringStereoAtoms", PropertyValue::UInt(3))
            .unwrap();
        let b = graph(0);
        let mut missing = output(&b);
        missing.source_orders_written = false;
        assert!(matches!(
            compose_query_cx_templates(&[(&a, &output(&a)), (&b, &missing)], F::NONE),
            Err(SmartsWriteError::CxMissingOutputOrder { template: 1 })
        ));
    }
    #[test]
    fn insertion_errors_precede_order_value_cast_after_presence_preflight() {
        let mut q = graph(1);
        q.set_prop("_smilesAtomOutputOrder", PropertyValue::IntVector(vec![0]))
            .unwrap();
        q.set_prop("_smilesBondOutputOrder", PropertyValue::IntVector(vec![]))
            .unwrap();
        let mut o = output(&q);
        o.source_orders_written = false;
        assert!(matches!(
            compose_query_cx_templates(&[(&q, &o)], F::NONE),
            Err(SmartsWriteError::CxOutputOrderPropertyType {
                property: "_smilesAtomOutputOrder",
                ..
            })
        ));
        q.atom_mut(0)
            .unwrap()
            .set_prop("_ringStereoAtoms", PropertyValue::UInt(3))
            .unwrap();
        assert!(matches!(
            compose_query_cx_templates(&[(&q, &o)], F::NONE),
            Err(SmartsWriteError::CxAtomPropertyKind { .. })
        ));
    }
    #[test]
    fn flags_off_still_insert_all_mixed_conformers_and_duplicate_ids() {
        let mut q = graph(1);
        q.add_conformer_3d(Conformer3D::new(7, vec![[1., 2., 3.]], false))
            .unwrap();
        q = q.with_2d_coordinate_block(vec![[4., 5.]]).unwrap();
        q.add_conformer_3d(Conformer3D::new(7, vec![[6., 7., 8.]], true))
            .unwrap();
        let r = run(&q);
        let c = r.graph.conformers_3d();
        assert_eq!(c.len(), 3);
        assert_eq!(c.iter().map(|x| x.id()).collect::<Vec<_>>(), [7, 0, 7]);
        assert_eq!(c[1].coordinates(), [[4., 5., 0.]]);
        assert!(!c[0].is_3d());
        assert!(!c[1].is_3d());
        assert!(c[2].is_3d());
        assert!(r.extension.is_empty());
    }
    #[test]
    fn equal_counts_pair_by_position_keep_first_ids_and_flags() {
        let (mut a, mut b) = (graph(1), graph(1));
        for (id, x, flag) in [(8, 1., false), (3, 2., true)] {
            a.add_conformer_3d(Conformer3D::new(id, vec![[x, 0., 0.]], flag))
                .unwrap();
        }
        for (id, x, flag) in [(3, 4., true), (8, 5., false)] {
            b.add_conformer_3d(Conformer3D::new(id, vec![[x, 0., 0.]], flag))
                .unwrap();
        }
        let r =
            compose_query_cx_templates(&[(&a, &output(&a)), (&b, &output(&b))], F::NONE).unwrap();
        let c = r.graph.conformers_3d();
        assert_eq!(c[0].id(), 8);
        assert!(!c[0].is_3d());
        assert_eq!(c[0].coordinates(), [[1., 0., 0.], [4., 0., 0.]]);
        assert_eq!(c[1].id(), 3);
        assert_eq!(c[1].coordinates(), [[2., 0., 0.], [5., 0., 0.]]);
    }
    #[test]
    fn unequal_counts_only_extend_existing_frames_with_zero_positions() {
        let (mut a, mut b) = (graph(1), graph(1));
        for id in [8, 3] {
            a.add_conformer_3d(Conformer3D::new(id, vec![[1., 2., 3.]], true))
                .unwrap();
        }
        b.add_conformer_3d(Conformer3D::new(9, vec![[4., 5., 6.]], true))
            .unwrap();
        let r =
            compose_query_cx_templates(&[(&a, &output(&a)), (&b, &output(&b))], F::NONE).unwrap();
        assert_eq!(r.graph.conformers_3d().len(), 2);
        for c in r.graph.conformers_3d() {
            assert_eq!(c.coordinates(), [[1., 2., 3.], [0.; 3]]);
        }
    }
    #[test]
    fn native_add_conformer_keeps_nonfinite_and_duplicate_id_values() {
        let mut q = graph(1);
        let bits = 0x7ff8_0000_0000_0042;
        for _ in 0..2 {
            q.add_conformer_3d(Conformer3D::new(
                3,
                vec![[f64::from_bits(bits), f64::INFINITY, f64::NEG_INFINITY]],
                false,
            ))
            .unwrap();
        }
        let r = run(&q);
        for c in r.graph.conformers_3d() {
            assert_eq!(c.id(), 3);
            assert_eq!(c.coordinates()[0][0].to_bits(), bits);
            assert_eq!(c.coordinates()[0][1], f64::INFINITY);
        }
        assert_eq!(q.conformers_3d()[0].coordinates()[0][0].to_bits(), bits);
    }
    #[test]
    fn unsigned_offset_wraps_before_final_global_access() {
        let a = graph(1);
        let b = graph(1);
        let mut o = output(&b);
        o.atom_order = vec![AtomId::new(u32::MAX as usize)];
        let r = compose_query_cx_templates(&[(&a, &output(&a)), (&b, &o)], F::NONE).unwrap();
        assert_eq!(r.atom_order, [AtomId::new(0), AtomId::new(0)]);
    }
    #[test]
    fn local_out_of_range_rows_may_resolve_to_actual_other_template_atoms() {
        let a = graph(1);
        let b = graph(1);
        let mut o = output(&a);
        o.atom_order = vec![AtomId::new(1)];
        let r = compose_query_cx_templates(&[(&a, &o), (&b, &output(&b))], F::NONE).unwrap();
        assert_eq!(r.atom_order, [AtomId::new(1), AtomId::new(1)]);
    }
    #[test]
    fn final_flags_decide_if_invalid_bond_order_is_actually_read() {
        let q = graph(1);
        let mut o = output(&q);
        o.bond_order = vec![BondId::new(99)];
        assert!(compose_query_cx_templates(&[(&q, &o)], F::NONE).is_ok());
        assert!(compose_query_cx_templates(&[(&q, &o)], F::ZERO_BONDS).is_err());
    }
    #[test]
    fn final_single_writer_emits_atom_label_and_acquires_actual_ring_cache() {
        let mut q = QueryGraph::from_parts(
            graph(8).atoms().to_vec(),
            (0..8)
                .map(|i| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(
                            AtomId::new(i),
                            AtomId::new((i + 1) % 8),
                            if i == 0 {
                                BondOrder::Double
                            } else {
                                BondOrder::Single
                            },
                        ),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        q.atom_mut(2)
            .unwrap()
            .set_prop("atomLabel", "source")
            .unwrap();
        let before = q.clone();
        let r =
            compose_query_cx_templates(&[(&q, &output(&q))], F::ATOM_LABELS | F::BOND_CFG).unwrap();
        assert!(String::from_utf8_lossy(r.extension.as_bytes()).contains("source"));
        assert!(r.graph.source_ring_info().initialized);
        assert_eq!(q, before);
    }
}
