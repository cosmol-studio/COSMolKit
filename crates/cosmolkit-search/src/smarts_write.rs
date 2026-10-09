//! RDKit SMARTS serialization for canonical [`QueryGraph`] values.

use cosmolkit_model::{
    AdjacencyList, AtomId, AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, Bond,
    BondId, BondQueryPredicate, PropertyText, QueryAtom, QueryBond, QueryGraph, QueryNode,
    RecursiveStructureQuery, SGroupConnection, StereoGroup, StereoGroupKind, SubstanceGroup,
    SubstanceGroupKind, query_substance_groups,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Element, Hybridization};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

#[cfg(not(feature = "smiles-integration"))]
#[path = "smarts_write/minimal_traversal.rs"]
mod minimal_traversal;

/// Options for QueryGraph-native SMARTS serialization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmartsWriteParams {
    /// Emit atom-map labels stored on query atoms.
    pub include_atom_maps: bool,
    /// Emit atom and bond stereochemistry.
    pub isomeric_smiles: bool,
    /// Preserve directional dative bond tokens.
    pub include_dative_bonds: bool,
    /// Start graph traversal at this atom index when present.
    pub rooted_at_atom: Option<usize>,
}

impl Default for SmartsWriteParams {
    fn default() -> Self {
        Self {
            include_atom_maps: true,
            isomeric_smiles: true,
            include_dative_bonds: true,
            rooted_at_atom: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct QueryBoolFeatures(u32);

impl QueryBoolFeatures {
    const HAS_AND: Self = Self(0x1);
    const HAS_LOW_AND: Self = Self(0x2);
    const HAS_OR: Self = Self(0x4);
    const HAS_RECURSION: Self = Self(0x8);

    const fn contains(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }

    fn insert(&mut self, other: Self) {
        self.0 |= other.0;
    }
}

impl std::ops::BitOrAssign for QueryBoolFeatures {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SmartsWriteError {
    #[error("{0}")]
    StereoGroup(#[from] cosmolkit_model::StereoGroupError),

    #[cfg(feature = "smiles-integration")]
    #[error("source SMARTS traversal preparation failed: {0}")]
    Traversal(#[from] cosmolkit_smiles::SmartsTraversalError),
    #[error("source CX coordinate selection failed: {0}")]
    CxCoordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error("CX source RingInfo: {0}")]
    CxRingInfo(#[source] cosmolkit_core::RingFindingError),
    #[error("CX ring atomOrder index {index} is outside {count} items")]
    CxRingAtomOrderIndex { index: usize, count: usize },
    #[error("CX ring bond {bond} has no source stereo atoms for side {side}")]
    CxRingStereoReferenceMissing { bond: BondId, side: usize },

    #[error("CX bond wedge: {0}")]
    CxWedge(#[source] cosmolkit_core::WedgeError),
    #[error("Internal error - should not occur (atrop bond {bond})")]
    CxBondConfigAtropMissingCarriers { bond: BondId },

    #[cfg(feature = "smiles-integration")]
    #[error("query CX atom-property output: {0}")]
    CxAtomPropertyOutput(#[source] cosmolkit_smiles::SmilesParseError),

    #[error("No conformations available on the molecule")]
    CxMissingConformer,
    #[error("query CX source conformer order: {0}")]
    CxCoordinateSource(#[source] cosmolkit_model::CoordinateValidationError),
    #[cfg(feature = "smiles-integration")]
    #[error("query CX coordinate output: {0}")]
    CxCoordinateOutput(#[source] cosmolkit_smiles::SmilesParseError),

    #[error("CX source bond {bond} is outside order of length {bond_count}")]
    CxSourceBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("CX molecule property {property}: {source}")]
    CxMoleculePropertyUInt {
        property: &'static str,
        #[source]
        source: cosmolkit_core::PropertyUIntReadError,
    },
    #[error("SGroup {group:?} property {property}: bad_any_cast from {actual:?} to {expected}")]
    CxSgroupVectorCast {
        group: cosmolkit_model::SubstanceGroupId,
        property: &'static str,
        actual: cosmolkit_model::PropertyValueKind,
        expected: &'static str,
    },

    #[error("SGroup {group:?} property {property} write: {source}")]
    CxSgroupPropertyWrite {
        group: cosmolkit_model::SubstanceGroupId,
        property: &'static str,
        #[source]
        source: cosmolkit_model::MoleculePropertyError,
    },

    #[error("CX source stereo-group collection: {0}")]
    CxStereoGroup(#[from] cosmolkit_core::AtropisomerError),

    #[error("CX source atom {atom} is outside reverse order of length {atom_count}")]
    CxSourceAtomOutOfRange { atom: AtomId, atom_count: usize },

    #[error("SMARTS source molecule property write: {0}")]
    MoleculePropertyWrite(#[from] cosmolkit_model::MoleculePropertyError),
    #[error("SMARTS source atom count {count} exceeds unsigned int")]
    SourceAtomCount { count: usize },

    #[error("atom {atom} source property write: {source}")]
    AtomPropertyWrite {
        atom: AtomId,
        #[source]
        source: cosmolkit_model::AtomPropertyError,
    },

    #[cfg(feature = "smiles-integration")]
    #[error("SMARTS canonical traversal: {0}")]
    CanonicalTraversal(#[from] cosmolkit_smiles::SmilesParseError),

    #[error("Can't write smarts for this query bond type: {description}")]
    UnwritableBondQuery { description: &'static str },
    #[error("atomToLeftIdx {index} exceeds source signed int")]
    SourceAtomToLeftIndex { index: usize },
    #[error("bond begin index {index} exceeds source unsigned int")]
    SourceBondBeginIndex { index: usize },
    #[error("atom {atom} property molAtomMapNumber: {source}")]
    AtomMapInt {
        atom: AtomId,
        #[source]
        source: cosmolkit_core::PropertyIntReadError,
    },
    #[error("Atomic number not found: {atomic_number}")]
    AtomTypeAtomicNumber { atomic_number: i32 },
    #[error("source abs(INT_MIN) is undefined for formal-charge query {value}")]
    ChargeMagnitudeOverflow { value: i32 },
    #[error("CX property value has the wrong source type: {0}")]
    PropertyValue(#[from] cosmolkit_model::PropertyValueError),
    #[error("CX required property read failed at atom {atom:?}: {source}")]
    CxRequiredProperty {
        atom: AtomId,
        #[source]
        source: cosmolkit_core::RequiredPropertyStringError,
    },
    #[error("property listing failed on atom {atom}: {source}")]
    CxPropertyList {
        atom: AtomId,
        #[source]
        source: cosmolkit_model::AtomPropertyError,
    },
    #[error("received {actual} coordinate selectors for {expected} CX templates")]
    CxCoordinateSelectionArity { expected: usize, actual: usize },
    #[error("template {template} does not have the source SMARTS output-order properties")]
    CxMissingOutputOrder { template: usize },
    #[error(
        "template {template} property {property}: bad_any_cast reading string as source unsigned output-order vector"
    )]
    CxOutputOrderPropertyType {
        template: usize,
        property: &'static str,
    },
    #[error("template {template}, atom {atom}, property {property}: invalid kind {kind:?}")]
    CxAtomPropertyKind {
        template: usize,
        atom: AtomId,
        property: &'static str,
        kind: cosmolkit_model::PropertyValueKind,
    },
    #[error("template {template}, atom {atom}, property {property}: {source}")]
    CxAtomPropertyUInt {
        template: usize,
        atom: AtomId,
        property: &'static str,
        #[source]
        source: cosmolkit_core::PropertyUIntReadError,
    },
    #[error("template {template}, atom {atom}, property {property} write: {source}")]
    CxAtomPropertyWrite {
        template: usize,
        atom: AtomId,
        property: &'static str,
        #[source]
        source: cosmolkit_model::AtomPropertyError,
    },
    #[error("template {template} source conformer update: {source}")]
    CxCoordinateStorage {
        template: usize,
        #[source]
        source: cosmolkit_model::CoordinateValidationError,
    },
    #[cfg(feature = "smiles-integration")]
    #[error("template {template} coordinate selection: {source}")]
    CxCoordinateSelection {
        template: usize,
        #[source]
        source: cosmolkit_smiles::SmilesParseError,
    },
    #[error("query CX composition: {0}")]
    CxComposition(#[from] cosmolkit_model::QueryGraphError),
    #[error("template {template} {kind} output-order row {row} exceeds {count} rows")]
    CxOutputOrder {
        template: usize,
        kind: &'static str,
        row: usize,
        count: usize,
    },
    #[error("CX {kind} row count {count} exceeds source unsigned int")]
    CxRowCount { kind: &'static str, count: usize },
    #[error("bond {bond} property {property}: {source}")]
    CxBondPropertyUInt {
        bond: BondId,
        property: &'static str,
        #[source]
        source: cosmolkit_core::PropertyUIntReadError,
    },
    #[error("bond {bond} property {property} has invalid kind {kind:?}")]
    InvalidPropertyKind {
        bond: BondId,
        property: &'static str,
        kind: cosmolkit_model::PropertyValueKind,
    },
    #[error("atom {atom} property {property}: {source}")]
    CxAtomPropertyInt {
        atom: AtomId,
        property: &'static str,
        #[source]
        source: cosmolkit_core::PropertyIntReadError,
    },
    #[error("SGroup {group:?} property {property}: {source}")]
    CxSgroupPropertyUInt {
        group: cosmolkit_model::SubstanceGroupId,
        property: &'static str,
        #[source]
        source: cosmolkit_core::PropertyUIntReadError,
    },
    #[error("SMARTS property string conversion failed: {0}")]
    Property(#[from] cosmolkit_core::PropertyStringError),
    #[error("SMARTS writer periodic-table lookup failed: {0}")]
    Valence(#[from] cosmolkit_core::ValenceError),
    #[error("query graph is invalid: {0}")]
    InvalidGraph(String),
    #[error("SMARTS writer query-graph traversal is not available for this graph: {detail}")]
    QueryGraphTraversalUnsupported { detail: &'static str },
    #[error("This is a non-smartable query - OR above and below AND in the binary tree")]
    OrAboveAndBelowAnd,
    #[error("Don't know how to combine using {description}")]
    UnknownCombination { description: String },
    #[error("recursive SMARTS query has no query molecule")]
    MissingRecursiveQueryMolecule,
    #[error("Can't write smarts for this bond dir type: {direction:?}")]
    SourceBondDirection { direction: BondDirection },
    #[error("Can't write smarts for this query bond type: {predicate:?}")]
    UnsupportedBondQuery { predicate: BondQueryPredicate },
    #[error("Can't write smarts for this query atom type: {predicate:?}")]
    UnsupportedAtomQuery { predicate: AtomQueryPredicate },
    #[error("SMARTS {kind} composite query requires at least two children")]
    CompositeChildCount { kind: &'static str },
    #[error("SMARTS writer does not support XOR query composites")]
    XorComposite,
    #[error("CXSMARTS extensions for QueryGraph are not yet supported: {detail}")]
    QueryGraphCxExtensionsUnsupported { detail: &'static str },
    #[error("rooted atom index {atom} is out of range")]
    RootedAtomOutOfRange { atom: usize },
    #[error("SMARTS fragment requires at least one atom")]
    EmptyAtomSelection,
    #[error("an explicit SMARTS fragment bond selection cannot be empty")]
    EmptyBondSelection,
    #[error("SMARTS fragment atom index {atom} is out of range")]
    FragmentAtomOutOfRange { atom: usize },
    #[error("SMARTS fragment bond index {bond} is out of range")]
    FragmentBondOutOfRange { bond: usize },
    #[error("atom {atom} is not an endpoint of bond {bond}")]
    BondAtomNotEndpoint { bond: usize, atom: usize },
}

/// Serialize concrete detached rows through the existing SMARTS traversal.
#[doc(hidden)]
#[cfg(feature = "smiles-integration")]
pub fn topology_to_smarts(
    topology: &cosmolkit_model::TopologyBlock,
    coordinates: &cosmolkit_model::CoordinateBlock,
    properties: &cosmolkit_model::MoleculeProperties,
    params: &SmartsWriteParams,
    include_cx: bool,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit✔️✔️:   const unsigned int nAtoms = mol.getNumAtoms();
    // RDKit✔️✔️:   if (!nAtoms) {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    if topology.atoms.is_empty() {
        return Ok(PropertyText::new());
    }
    // RDKit❗❌:   ROMol mol(inmol);
    // RDKit❗❌:   for (auto &atom : mol.atoms()) {
    // RDKit❗❌:     atom->updatePropertyCache(false);
    // RDKit❗❌:   }
    if let Some(atom) = params.rooted_at_atom {
        if atom >= topology.atoms.len() {
            return Err(SmartsWriteError::RootedAtomOutOfRange { atom });
        }
    }
    // RDKit✔️✔️:   ROMol mol(inmol);
    // molToSmarts owns the sole Canon traversal. Running the old preparation
    // traversal here as well inverts tetrahedral tags twice and carries its
    // _TraversalRingClosureBond scratch into the real emission traversal.
    // Move the detached carrier into the existing writer without that pass.
    let prepared = topology.clone();
    let mut effective = *params;
    if include_cx {
        effective.include_dative_bonds = false;
    }
    let prepared_graph = concrete_smarts_graph(prepared, None, properties)?;
    let output = query_graph_to_smarts_output(&prepared_graph, &effective).map_err(|error| {
        // Preserve the public concrete-writer property-error contract while
        // retaining the exact checked signed-int conversion and its payload.
        match error {
            SmartsWriteError::AtomMapInt { atom, source } => SmartsWriteError::CxAtomPropertyInt {
                atom,
                property: "molAtomMapNumber",
                source,
            },
            error => error,
        }
    })?;
    let mut text = output.text;
    if include_cx && !text.is_empty() {
        // RDKit✔️✔️:   auto res = MolToSmarts(mol, ps);
        // RDKit✔️✔️:     auto cxext = SmilesWrite::getCXExtensions(mol);
        // Native extensions read the ORIGINAL molecule, not the temporary
        // traversal's inverted tags/directions or computed CIP properties.
        // Reusing QueryGraph CX serialization needs a second detached carrier
        // here: extra O(V+E) copy cost, explicitly retained as ❌ above.
        drop(prepared_graph);
        let mut graph = concrete_smarts_graph(topology.clone(), Some(coordinates), properties)?;
        let extension = write_query_cx_extensions(
            &mut graph,
            &output.atom_order,
            &output.bond_order,
            cosmolkit_smiles::CxSmilesFields::ALL,
        )?;
        if !extension.is_empty() {
            text.push_byte(b' ');
            text.extend_bytes(extension.as_bytes());
        }
    }
    Ok(text)
}

#[cfg(feature = "smiles-integration")]
fn concrete_smarts_graph(
    topology: cosmolkit_model::TopologyBlock,
    coordinates: Option<&cosmolkit_model::CoordinateBlock>,
    properties: &cosmolkit_model::MoleculeProperties,
) -> Result<QueryGraph, SmartsWriteError> {
    // Rust-only transport into the existing uniform graph carrier; origin,
    // ordered properties and groups remain exact. No predicates are exposed
    // as queries for ordinary rows. Atoms/bonds move rather than clone again.
    let atoms = topology
        .atoms
        .into_iter()
        .map(|atom| {
            let predicate =
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number()));
            QueryAtom::from_carrier_parts(atom, predicate)
        })
        .collect();
    let bonds = topology
        .bonds
        .into_iter()
        .map(|bond| {
            let predicate = QueryNode::predicate(BondQueryPredicate::Order(bond.order()));
            QueryBond::from_carrier_parts(bond, predicate)
        })
        .collect();
    let mut graph = QueryGraph::from_parts(
        atoms,
        bonds,
        properties
            .ordered_props()
            .map(|(key, value)| (key.clone(), value.clone())),
        coordinates.map_or_else(Vec::new, |coordinates| coordinates.conformers_2d.clone()),
        coordinates.map_or_else(Vec::new, |coordinates| coordinates.conformers_3d.clone()),
        topology.stereo_groups,
    )?;
    if let Some(coordinates) = coordinates {
        graph.set_source_conformer_order(coordinates.source_conformer_order.clone())?;
    }
    cosmolkit_model::replace_query_substance_groups(&mut graph, topology.substance_groups)?;
    Ok(graph)
}

/// Serialize an independent query graph.
///
/// Uses the shared detached source Canon traversal with explicit query state.
/// No live `Molecule` or runtime storage is constructed.
pub fn query_graph_to_smarts(
    query: &QueryGraph,
    params: &SmartsWriteParams,
) -> Result<PropertyText, SmartsWriteError> {
    query_graph_to_smarts_result(query, params).map(|result| result.smarts)
}

fn query_graph_to_smarts_result(
    query: &QueryGraph,
    params: &SmartsWriteParams,
) -> Result<SmartsWriteResult, SmartsWriteError> {
    #[cfg(feature = "smiles-integration")]
    {
        mol_to_smarts_wrapper_source(query, params)
    }
    #[cfg(not(feature = "smiles-integration"))]
    {
        minimal_traversal::write(query, params, None, None, false)
    }
}

/// Serialize with the existing source traversal and retain output-order
/// evidence for reaction-wide CX composition.
#[doc(hidden)]
pub fn query_graph_to_smarts_output(
    query: &QueryGraph,
    params: &SmartsWriteParams,
) -> Result<SmartsWriteOutput, SmartsWriteError> {
    let result = query_graph_to_smarts_result(query, params)?;
    Ok(SmartsWriteOutput {
        text: result.smarts,
        atom_order: result.atom_ordering,
        bond_order: result.bond_ordering,
        source_orders_written: result.source_orders_written,
    })
}

/// Detached writer evidence; no runtime authority or duplicated graph.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmartsWriteOutput {
    pub text: PropertyText,
    pub atom_order: Vec<AtomId>,
    pub bond_order: Vec<BondId>,
    /// Explicit proof of reaching the SOURCE order-property recording branch.
    pub source_orders_written: bool,
}

#[cfg(feature = "smiles-integration")]
pub fn query_graph_to_cx_smarts(
    query: &QueryGraph,
    params: &SmartsWriteParams,
) -> Result<PropertyText, SmartsWriteError> {
    // The existing immutable detached value API stages one actual query.
    // Native const ROMol mutates properties, bonds and ring cache. Those exact
    // effects belong to the sole mutable source kernel below, never runtime.
    // This API adapter's extra graph clone is an explicit known cost.
    let mut staged = query.clone();
    mol_to_cx_smarts_source(&mut staged, params)
}

pub fn query_graph_fragment_to_smarts(
    query: &QueryGraph,
    params: &SmartsWriteParams,
    atoms: &[AtomId],
    bonds: Option<&[BondId]>,
) -> Result<PropertyText, SmartsWriteError> {
    #[cfg(feature = "smiles-integration")]
    {
        mol_fragment_to_smarts_source(query, params, atoms, bonds).map(|output| output.smarts)
    }
    #[cfg(not(feature = "smiles-integration"))]
    {
        minimal_traversal::write(query, params, Some(atoms), bonds, false)
            .map(|output| output.smarts)
    }
}

#[cfg(feature = "smiles-integration")]
pub fn query_graph_fragment_to_cx_smarts(
    query: &QueryGraph,
    params: &SmartsWriteParams,
    atoms: &[AtomId],
    bonds: Option<&[BondId]>,
) -> Result<PropertyText, SmartsWriteError> {
    // Thin immutable detached value adapter; Native's same-graph property,
    // bond/cache mutation effects are implemented in the source owner below.
    // One explicit extra clone, without a second traversal or CX algorithm.
    let mut staged = query.clone();
    mol_fragment_to_cx_smarts_source(&mut staged, params, atoms, bonds)
}

#[cfg(feature = "smiles-integration")]
fn append_query_cx_extension(addition: impl AsRef<[u8]>, output: &mut PropertyText) {
    // RDKit❗✔️: void appendToCXExtension(const std::string &addition, std::string &base) {
    // RDKit❗✔️:   if (!addition.empty()) {
    // RDKit❗✔️:     if (base.size() > 1) {
    // RDKit❗✔️:       base += ",";
    // RDKit❗✔️:     }
    // RDKit❗✔️:     base += addition;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    cosmolkit_smiles::append_cx_extension_source(addition, output);
}

#[cfg(feature = "smiles-integration")]
fn query_cx_atom_positions(atom_order: &[AtomId], atom_count: usize) -> Vec<Option<usize>> {
    let mut positions = vec![None; atom_count];
    for (position, atom) in atom_order.iter().copied().enumerate() {
        positions[atom.index()] = Some(position);
    }
    positions
}

#[cfg(feature = "smiles-integration")]
fn query_cx_source_reverse_atom_order(atom_order: &[AtomId], atom_count: usize) -> Vec<usize> {
    // RDKit value-initializes the reverse vector. For fragment CXSMARTS this
    // intentionally leaves every unselected source atom mapped to output 0.
    let mut positions = vec![0; atom_count];
    for (position, atom) in atom_order.iter().copied().enumerate() {
        positions[atom.index()] = position;
    }
    positions
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_coordinates(
    query: &QueryGraph,
    atom_order: &[AtomId],
) -> Result<String, SmartsWriteError> {
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
    let source = query
        .first_source_conformer()
        .map_err(SmartsWriteError::CxCoordinateSource)?
        .ok_or(SmartsWriteError::CxMissingConformer)?;
    let source = match source {
        cosmolkit_model::CoordinateSourceConformer::TwoD(row) => {
            cosmolkit_smiles::CoordinateSource::TwoD(row)
        }
        cosmolkit_model::CoordinateSourceConformer::ThreeD(row) => {
            cosmolkit_smiles::CoordinateSource::ThreeD(row)
        }
    };
    // MODEL alone owns source-front selection; SMILES alone owns CX number and
    // coordinate serialization. Both consume the original borrowed conformer.
    // Empty storage is Native getConformer failure, even for empty atomOrder;
    // the enclosing Native getNumConformers guard owns the legitimate skip.
    // Cost ❌: canonical number formatting/output retain the documented owned
    // scalar buffers and no SSO. This adapter does no allocation or clone.
    cosmolkit_smiles::write_cx_coordinates_from_source(source, atom_order)
        .map_err(SmartsWriteError::CxCoordinateOutput)
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_atom_labels(
    query: &QueryGraph,
    atom_order: &[AtomId],
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string get_atomlabel_block(const ROMol &mol,
    // RDKit❗❌:                                 const std::vector<unsigned int> &atomOrder) {
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:   for (auto idx : atomOrder) {
    // RDKit❗❌:     if (idx != atomOrder.front()) {
    // RDKit❗❌:       res += ";";
    // RDKit❗❌:     }
    // RDKit❗❌:     std::string lbl;
    // RDKit❗❌:     int val;
    // RDKit❗❌:     const auto atom = mol.getAtomWithIdx(idx);
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::_QueryAtomGenericLabel,
    // RDKit❗❌:                                lbl)) {
    // RDKit❗❌:       res += quote_string(lbl + "_p");
    // RDKit❗❌:     } else if (!atom->getAtomicNum() &&
    // RDKit❗❌:                atom->getPropIfPresent(common_properties::dummyLabel, lbl) &&
    // RDKit❗❌:                std::find(SmilesParseOps::pseudoatoms.begin(),
    // RDKit❗❌:                          SmilesParseOps::pseudoatoms.end(),
    // RDKit❗❌:                          lbl) != SmilesParseOps::pseudoatoms.end()) {
    // RDKit❗❌:       res += quote_string(lbl + "_p");
    // RDKit❗❌:     } else if (!atom->getAtomicNum() &&
    // RDKit❗❌:                atom->getPropIfPresent(common_properties::_fromAttachPoint,
    // RDKit❗❌:                                       val) &&
    // RDKit❗❌:                (val == 1 || val == 2)) {
    // RDKit❗❌:       res += quote_string("_AP" + std::to_string(val));
    // RDKit❗❌:     } else if (atom->getPropIfPresent(common_properties::atomLabel, lbl)) {
    // RDKit❗❌:       res += quote_string(lbl);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   // if we didn't find anything return an empty string
    // RDKit❗❌:   if (std::find_if_not(res.begin(), res.end(),
    // RDKit❗❌:                        [](const auto c) { return c == ';'; }) == res.end()) {
    // RDKit❗❌:     res.clear();
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // RDKit❗✔️: constexpr std::array<std::string_view, 2> pseudoatoms{"Pol", "Mod"};
    // RDKit❗✔️: inline constexpr std::string_view _fromAttachPoint = "_fromAttchpt";
    const PSEUDOATOMS: [&[u8]; 2] = [b"Pol", b"Mod"];
    let Some(first) = atom_order.first() else {
        return Ok(PropertyText::new());
    };
    let mut result = PropertyText::new();
    for id in atom_order {
        // Native compares the atom identity to front(), not loop position.
        if id != first {
            result.push_byte(b';');
        }
        let atom =
            query
                .atoms()
                .get(id.index())
                .ok_or(SmartsWriteError::CxSourceAtomOutOfRange {
                    atom: *id,
                    atom_count: query.num_atoms(),
                })?;
        if let Some(value) = atom.prop("_QueryAtomGenericLabel") {
            let mut label = cosmolkit_core::property_value_to_string(value)?;
            label.extend_bytes(b"_p");
            result.extend_bytes(quote_query_cx_string(label.as_bytes()).as_bytes());
            continue;
        }
        if atom.atomic_number() == 0 {
            if let Some(value) = atom.prop("dummyLabel") {
                let mut label = cosmolkit_core::property_value_to_string(value)?;
                if PSEUDOATOMS.contains(&label.as_bytes()) {
                    label.extend_bytes(b"_p");
                    result.extend_bytes(quote_query_cx_string(label.as_bytes()).as_bytes());
                    continue;
                }
            }
            if let Some(value) = atom.prop("_fromAttchpt") {
                let value = cosmolkit_core::property_value_to_int(value).map_err(|source| {
                    SmartsWriteError::CxAtomPropertyInt {
                        atom: *id,
                        property: "_fromAttchpt",
                        source,
                    }
                })?;
                if matches!(value, 1 | 2) {
                    let label = if value == 1 { b"_AP1" } else { b"_AP2" };
                    result.extend_bytes(quote_query_cx_string(label).as_bytes());
                    continue;
                }
            }
        }
        if let Some(value) = atom.prop("atomLabel") {
            let label = cosmolkit_core::property_value_to_string(value)?;
            result.extend_bytes(quote_query_cx_string(label.as_bytes()).as_bytes());
        }
    }
    if result.as_bytes().iter().all(|c| *c == b';') {
        result.clear();
    }
    // Cost ❌: counted strings have no Native SSO. Source branch order and
    // one linear output buffer are retained; no identity/predicate conversion,
    // property scans or rewritten AST, and quote_string is the sole copy owner.
    Ok(result)
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_atom_values(
    query: &QueryGraph,
    atom_order: &[AtomId],
    property: &[u8],
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string get_value_block(const ROMol &mol,
    // RDKit❗❌:                             const std::vector<unsigned int> &atomOrder,
    // RDKit❗❌:                             const std::string_view &prop) {
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:   bool first = true;
    // RDKit❗❌:   for (auto idx : atomOrder) {
    // RDKit❗❌:     if (!first) {
    // RDKit❗❌:       res += ";";
    // RDKit❗❌:     } else {
    // RDKit❗❌:       first = false;
    // RDKit❗❌:     }
    // RDKit❗❌:     std::string lbl;
    // RDKit❗❌:     if (mol.getAtomWithIdx(idx)->getPropIfPresent(prop, lbl)) {
    // RDKit❗❌:       res += quote_string(lbl);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    let mut result = PropertyText::new();
    for (position, id) in atom_order.iter().enumerate() {
        if position != 0 {
            result.push_byte(b';');
        }
        let atom =
            query
                .atoms()
                .get(id.index())
                .ok_or(SmartsWriteError::CxSourceAtomOutOfRange {
                    atom: *id,
                    atom_count: query.num_atoms(),
                })?;
        if let Some(value) = atom.prop(property) {
            let label = cosmolkit_core::property_value_to_string(value)?;
            result.extend_bytes(quote_query_cx_string(label.as_bytes()).as_bytes());
        }
    }
    // Source returns all semicolons, unlike get_atomlabel_block. The enclosing
    // presence preflight owns whether $_AV is emitted, not this helper.
    // Cost ❌: no Native SSO; source scalar/quote copy and one linear output
    // retained without per-position vector or UTF-8 conversion of key/value.
    Ok(result)
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_radicals(
    query: &QueryGraph,
    atom_order: &[AtomId],
    warn: &mut dyn FnMut(u32),
) -> Result<String, SmartsWriteError> {
    // RDKit❗🔝: std::string get_radical_block(const ROMol &mol,
    // RDKit❗🔝:                               const std::vector<unsigned int> &atomOrder) {
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
    // RDKit❗✔️: unsigned int getNumRadicalElectrons() const { return d_numRadicalElectrons; }
    let mut radicals = BTreeMap::<u32, Vec<usize>>::new();
    for (position, id) in atom_order.iter().enumerate() {
        let atom =
            query
                .atoms()
                .get(id.index())
                .ok_or(SmartsWriteError::CxSourceAtomOutOfRange {
                    atom: *id,
                    atom_count: query.num_atoms(),
                })?;
        let count = u32::from(atom.radical_electrons());
        if count != 0 {
            radicals.entry(count).or_default().push(position);
        }
    }
    let mut output = String::new();
    for (count, positions) in radicals {
        match count {
            1 => output.push_str("^1:"),
            2 => output.push_str("^2:"),
            3 => output.push_str("^5:"),
            _ => warn(count),
        }
        // Source warnings never suppress the group's positions. The trailing
        // comma belongs to this helper; only the enclosing writer removes one.
        for position in positions {
            write!(&mut output, "{position},").expect("String formatting is infallible");
        }
    }
    // Cost 🔝: BTreeMap preserves Native ordered grouping and vector traversal.
    // Direct ASCII decimal appends preserve fixed "%d," under pinned C locale
    // and avoid constructing/parsing a boost::format, its dynamic item vector,
    // temporary stream/string and subsequent string append for every index.
    // Native format_implementation.hpp constructor calls parse(s), and
    // parsing.hpp allocates items_.resize(num_items); this is a concrete
    // per-index allocation/parser removal, not a benchmark-based claim.
    Ok(output)
}

#[cfg(feature = "smiles-integration")]
fn quote_query_cx_string(text: &[u8]) -> PropertyText {
    // RDKit❗❌: std::string quote_string(const std::string &txt) {
    // RDKit❗❌:   // FIX
    // RDKit❗❌:   return txt;
    // RDKit❗❌: }
    // Native deliberately returns an owned, unescaped copy, including NULs.
    // Cost ❌: the modeled counted byte vector has no C++ short-string storage.
    PropertyText::from_bytes(text)
}

#[cfg(feature = "smiles-integration")]
fn quote_query_cx_atom_property(text: &[u8]) -> PropertyText {
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
    // Reuse the existing SMILES byte serializer; no second escaping algorithm.
    // The used SEARCH adapter adds no allocation/copy; canonical output has
    // the known short-string allocation cost vs Native SSO (❌).
    cosmolkit_smiles::quote_cx_atom_property(text)
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_atom_properties(
    query: &QueryGraph,
    atom_order: &[AtomId],
) -> Result<PropertyText, SmartsWriteError> {
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
    // SMILES owns the one borrowed property serializer. The actual query
    // carrier identity/ordered properties are read directly, without coercing
    // to an ordinary Atom, dropping the predicate, or cloning detached state.
    // Cost ❌ remains the canonical scalar/quote output's missing SSO; this
    // used adapter adds no loops, temporary vectors, allocations or copies.
    cosmolkit_smiles::write_query_cx_atom_properties_source(query, atom_order).map_err(|error| {
        match error {
            cosmolkit_smiles::SmilesParseError::WriterPropertyList { atom, source } => {
                SmartsWriteError::CxPropertyList { atom, source }
            }
            cosmolkit_smiles::SmilesParseError::WriterRequiredProperty { atom, source } => {
                SmartsWriteError::CxRequiredProperty { atom, source }
            }
            cosmolkit_smiles::SmilesParseError::CxAtomPropertyAtomOutOfRange {
                atom,
                atom_count,
            } => SmartsWriteError::CxSourceAtomOutOfRange { atom, atom_count },
            other => SmartsWriteError::CxAtomPropertyOutput(other),
        }
    })
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_bond_config(
    query: &QueryGraph,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    coordinates_included: bool,
    wedge_bonds: &cosmolkit_core::WedgeAssignments,
    atropisomer_only: bool,
) -> Result<String, SmartsWriteError> {
    // RDKit❗❌: std::string get_bond_config_block(
    // RDKit❗❌:     const ROMol &mol, const std::vector<unsigned int> &atomOrder,
    // RDKit❗❌:     const std::vector<unsigned int> &bondOrder, bool coordsIncluded,
    // RDKit❗❌:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>> &wedgeBonds,
    // RDKit❗❌:     bool atropisomerOnly = false) {
    // RDKit❗❌:   std::map<std::string, std ::vector<std::string>> wParts;
    // RDKit❗❌:   for (unsigned int i = 0; i < bondOrder.size(); ++i) {
    // RDKit❗❌:     auto idx = bondOrder[i];
    // RDKit❗❌:     const auto bond = mol.getBondWithIdx(idx);
    // RDKit❗❌:     unsigned int wedgeStartAtomIdx = bond->getBeginAtomIdx();
    // RDKit❗❌:
    // RDKit❗❌:     if (!canHaveDirection(*bond)) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     // when figuring out what to output for the bond, favor the wedge state:
    // RDKit❗❌:     Bond::BondDir bd = bond->getBondDir();
    // RDKit❗❌:     switch (bd) {
    // RDKit❗❌:       case Bond::BondDir::BEGINDASH:
    // RDKit❗❌:       case Bond::BondDir::BEGINWEDGE:
    // RDKit❗❌:       case Bond::BondDir::UNKNOWN:
    // RDKit❗❌:         break;
    // RDKit❗❌:       default:
    // RDKit❗❌:         bd = Bond::BondDir::NONE;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (atropisomerOnly && bd == Bond::BondDir::NONE) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // see if this one is an atropisomer
    // RDKit❗❌:
    // RDKit❗❌:     bool isAnAtropisomer = false;
    // RDKit❗❌:
    // RDKit❗❌:     const Atom *firstAtom = bond->getBeginAtom();
    // RDKit❗❌:     if (bd == Bond::BondDir::BEGINDASH || bd == Bond::BondDir::BEGINWEDGE) {
    // RDKit❗❌:       for (auto bondNbr : mol.atomBonds(firstAtom)) {
    // RDKit❗❌:         if (bondNbr->getIdx() == bond->getIdx()) {
    // RDKit❗❌:           continue;  // a bond is not its own neighbor
    // RDKit❗❌:         }
    // RDKit❗❌:         if (bondNbr->getStereo() == Bond::BondStereo::STEREOATROPCW ||
    // RDKit❗❌:             bondNbr->getStereo() == Bond::BondStereo::STEREOATROPCCW) {
    // RDKit❗❌:           isAnAtropisomer = true;
    // RDKit❗❌:
    // RDKit❗❌:           // if it is for an atropisomer and there are no coords, check to see
    // RDKit❗❌:           // if the wedge needs to be flipped based on the smiles reordering
    // RDKit❗❌:           if (!coordsIncluded && isAnAtropisomer) {
    // RDKit❗❌:             Atropisomers::AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit❗❌:             if (!Atropisomers::getAtropisomerAtomsAndBonds(
    // RDKit❗❌:                     bondNbr, atomAndBondVecs, mol)) {
    // RDKit❗❌:               throw ValueErrorException("Internal error - should not occur");
    // RDKit❗❌:               // should not happen
    // RDKit❗❌:             } else {
    // RDKit❗❌:               unsigned int swaps = 0;
    // RDKit❗❌:
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
    // RDKit❗❌:
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
    // RDKit❗❌:
    // RDKit❗❌:                   unsigned int firstReorderedAtomIdx =
    // RDKit❗❌:                       std::find(atomOrder.begin(), atomOrder.end(),
    // RDKit❗❌:                                 firstOtherAtomIdx) -
    // RDKit❗❌:                       atomOrder.begin();
    // RDKit❗❌:                   unsigned int secondReorderedAtomIdx =
    // RDKit❗❌:                       std::find(atomOrder.begin(), atomOrder.end(),
    // RDKit❗❌:                                 secondOtherAtomIdx) -
    // RDKit❗❌:                       atomOrder.begin();
    // RDKit❗❌:
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
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:
    // RDKit❗❌:           break;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (atropisomerOnly) {
    // RDKit❗❌:       // one of the bonds on the beginning atom of this bond must be an
    // RDKit❗❌:       // atropisomer
    // RDKit❗❌:
    // RDKit❗❌:       if (!isAnAtropisomer) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {  //  atropisomeronly is FALSE - check for a wedging caused by
    // RDKit❗❌:               //  chiral atom
    // RDKit❗❌:       unsigned int cfg = 0;
    // RDKit❗❌:       if (bd == Bond::BondDir::NONE &&
    // RDKit❗❌:           bond->getPropIfPresent(common_properties::_MolFileBondCfg, cfg)) {
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
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       if (bd == Bond::BondDir::NONE && coordsIncluded) {
    // RDKit❗❌:         int dirCode;
    // RDKit❗❌:         bool reverse;
    // RDKit❗❌:         Chirality::GetMolFileBondStereoInfo(
    // RDKit❗❌:             bond, wedgeBonds, &mol.getConformer(), dirCode, reverse);
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
    // RDKit❗❌:         if (reverse) {
    // RDKit❗❌:           wedgeStartAtomIdx = bond->getEndAtomIdx();
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     auto begAtomOrder =
    // RDKit❗❌:         std::find(atomOrder.begin(), atomOrder.end(), wedgeStartAtomIdx) -
    // RDKit❗❌:         atomOrder.begin();
    // RDKit❗❌:
    // RDKit❗❌:     std::string wType = "";
    // RDKit❗❌:     if (bd == Bond::BondDir::UNKNOWN) {
    // RDKit❗❌:       wType = "w";
    // RDKit❗❌:     } else if (coordsIncluded || isAnAtropisomer) {
    // RDKit❗❌:       // we only do wedgeUp and wedgeDown if coordinates are being output
    // RDKit❗❌:       // or its an atropisomer
    // RDKit❗❌:       if (bd == Bond::BondDir::BEGINWEDGE) {
    // RDKit❗❌:         wType = "wU";
    // RDKit❗❌:       } else if (bd == Bond::BondDir::BEGINDASH) {
    // RDKit❗❌:         wType = "wD";
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (wType != "") {
    // RDKit❗❌:       if (wParts.find(wType) == wParts.end()) {
    // RDKit❗❌:         wParts[wType] = std::vector<std::string>();
    // RDKit❗❌:       }
    // RDKit❗❌:       wParts[wType].push_back(
    // RDKit❗❌:           boost::str(boost::format("%d.%d") % begAtomOrder % i));
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:
    // RDKit❗❌:   for (auto wPart : wParts) {
    // RDKit❗❌:     if (res != "") {
    // RDKit❗❌:       res += ",";
    // RDKit❗❌:     }
    // RDKit❗❌:     res += wPart.first + ":" + boost::algorithm::join(wPart.second, ",");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // Actual source-mutated atrop carrier endpoints/directions are retained in
    // the typed shared assignment. Chiral map entries do not mutate them and
    // are consulted only by the coordinate GetMolFile child below.
    let source_bond = |id: BondId| {
        query
            .bonds()
            .get(id.index())
            .map(cosmolkit_model::QueryBond::bond)
            .ok_or(SmartsWriteError::CxSourceBondOutOfRange {
                bond: id,
                bond_count: query.num_bonds(),
            })
    };
    let current_endpoints = |bond: &cosmolkit_model::Bond| match wedge_bonds.get(bond.id()) {
        Some(cosmolkit_core::WedgeInfo::Atropisomer { update }) => (update.begin, update.end),
        _ => (bond.begin(), bond.end()),
    };
    // Native std::find chooses the first occurrence and returns end on absent;
    // its unsigned local casts retain uint32 bits, including end-distance.
    let output_position = |id: AtomId| {
        atom_order
            .iter()
            .position(|&item| item == id)
            .unwrap_or(atom_order.len()) as u32
    };
    let mut parts = BTreeMap::<&'static str, Vec<String>>::new();
    for (i, &id) in bond_order.iter().enumerate() {
        let bond = source_bond(id)?;
        let (begin, end) = current_endpoints(bond);
        let mut wedge_start = begin;
        if !matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic) {
            continue;
        }
        let actual_direction = match wedge_bonds.get(id) {
            Some(cosmolkit_core::WedgeInfo::Atropisomer { update }) => update.direction,
            _ => bond.direction(),
        };
        let mut direction = match actual_direction {
            BondDirection::BeginDash | BondDirection::BeginWedge | BondDirection::Unknown => {
                actual_direction
            }
            _ => BondDirection::None,
        };
        if atropisomer_only && direction == BondDirection::None {
            continue;
        }
        let mut is_atropisomer = false;
        if matches!(
            direction,
            BondDirection::BeginDash | BondDirection::BeginWedge
        ) {
            for &(_, neighbor_id) in &query.adjacency()[begin.index()] {
                let neighbor_id = BondId::new(neighbor_id);
                if neighbor_id == id {
                    continue;
                }
                let neighbor = source_bond(neighbor_id)?;
                if matches!(
                    neighbor.stereo(),
                    BondStereo::AtropCw | BondStereo::AtropCcw
                ) {
                    is_atropisomer = true;
                    if !coordinates_included {
                        let carriers =
                            cosmolkit_core::query_atropisomer_carriers_source(query, neighbor_id)?
                                .ok_or(SmartsWriteError::CxBondConfigAtropMissingCarriers {
                                    bond: neighbor_id,
                                })?;
                        let (axial_begin, axial_end) = current_endpoints(neighbor);
                        let mut swaps =
                            u32::from(output_position(axial_begin) > output_position(axial_end));
                        for side in &carriers {
                            if side.focus() == begin {
                                continue;
                            }
                            if side.carrier_bonds().len() == 2 {
                                let other = |id: BondId| -> Result<AtomId, SmartsWriteError> {
                                    let carrier = source_bond(id)?;
                                    let (a, b) = current_endpoints(carrier);
                                    if a == side.focus() {
                                        Ok(b)
                                    } else if b == side.focus() {
                                        Ok(a)
                                    } else {
                                        Err(SmartsWriteError::CxWedge(
                                            cosmolkit_core::WedgeError::CenterNotIncident {
                                                bond: id,
                                                center: side.focus(),
                                            },
                                        ))
                                    }
                                };
                                if output_position(other(side.carrier_bonds()[0])?)
                                    > output_position(other(side.carrier_bonds()[1])?)
                                {
                                    swaps = swaps.wrapping_add(1);
                                }
                            }
                        }
                        if swaps % 2 != 0 {
                            direction = if direction == BondDirection::BeginWedge {
                                BondDirection::BeginDash
                            } else {
                                BondDirection::BeginWedge
                            };
                        }
                    }
                    break;
                }
            }
        }
        if atropisomer_only {
            if !is_atropisomer {
                continue;
            }
        } else {
            if direction == BondDirection::None {
                let cfg = bond
                    .prop("_MolFileBondCfg")
                    .map(cosmolkit_core::property_value_to_uint)
                    .transpose()
                    .map_err(|source| SmartsWriteError::CxBondPropertyUInt {
                        bond: bond.id(),
                        property: "_MolFileBondCfg",
                        source,
                    })?;
                direction = match cfg {
                    Some(1) => BondDirection::BeginWedge,
                    Some(2) => BondDirection::Unknown,
                    Some(3) => BondDirection::BeginDash,
                    _ => BondDirection::None,
                };
            }
            if direction == BondDirection::None && coordinates_included {
                let source = query
                    .first_source_conformer()
                    .map_err(SmartsWriteError::CxCoordinateSource)?
                    .ok_or(SmartsWriteError::CxMissingConformer)?;
                let conformer = match source {
                    cosmolkit_model::CoordinateSourceConformer::TwoD(c) => {
                        cosmolkit_core::AtropisomerConformer::TwoD(c)
                    }
                    cosmolkit_model::CoordinateSourceConformer::ThreeD(c) => {
                        cosmolkit_core::AtropisomerConformer::ThreeD(c)
                    }
                };
                let info = cosmolkit_core::get_query_directional_bond_stereo_info_source(
                    query,
                    wedge_bonds,
                    id,
                    Some(conformer),
                )
                .map_err(SmartsWriteError::CxWedge)?;
                direction = match info.direction_code {
                    1 => BondDirection::BeginWedge,
                    3 => BondDirection::Unknown,
                    6 => BondDirection::BeginDash,
                    _ => BondDirection::None,
                };
                if info.reverse {
                    wedge_start = end;
                }
            }
        }
        let position = output_position(wedge_start);
        let kind = match direction {
            BondDirection::Unknown => Some("w"),
            BondDirection::BeginWedge if coordinates_included || is_atropisomer => Some("wU"),
            BondDirection::BeginDash if coordinates_included || is_atropisomer => Some("wD"),
            _ => None,
        };
        if let Some(kind) = kind {
            parts
                .entry(kind)
                .or_default()
                .push(format!("{position}.{i}"));
        }
    }
    // Same lexical map grouping, per-entry buffers and ordered joins. Cost ❌:
    // Rust item strings lack Native SSO; actual projected mutation lookup adds
    // O(log W), while source reads the already-mutated bond directly. Reached
    // CORE kernels borrow QueryGraph; no eager whole-graph work or coercion.
    Ok(parts
        .into_iter()
        .map(|(kind, entries)| format!("{kind}:{}", entries.join(",")))
        .collect::<Vec<_>>()
        .join(","))
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_typed_bonds(
    query: &QueryGraph,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    order: BondOrder,
    symbol: &str,
) -> Result<PropertyText, SmartsWriteError> {
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
    // One SMILES-owned source emitter reads the actual QueryBond's Bond
    // facts without AST coercion or a second loop/serialization algorithm.
    cosmolkit_smiles::write_cx_coord_or_hydrogen_bonds_source(
        cosmolkit_core::stereo_graph::BondRows::Query(query.bonds()),
        atom_order,
        bond_order,
        order,
        symbol.as_bytes(),
    )
    .map_err(|error| match error {
        cosmolkit_smiles::SmilesParseError::CxTypedBondOutOfRange { bond, bond_count } => {
            SmartsWriteError::CxSourceBondOutOfRange { bond, bond_count }
        }
        other => SmartsWriteError::CanonicalTraversal(other),
    })
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_zero_bonds(
    query: &QueryGraph,
    bond_order: &[BondId],
) -> Result<PropertyText, SmartsWriteError> {
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
    // The used SMILES owner performs the sole actual bond-order loop. Cost
    // 🔝 inherits its direct output instead of Boost item-format machinery;
    // this adapter adds no loops, AST changes, carrier clones or buffers.
    cosmolkit_smiles::write_cx_zero_bonds_source(
        cosmolkit_core::stereo_graph::BondRows::Query(query.bonds()),
        bond_order,
    )
    .map_err(|error| match error {
        cosmolkit_smiles::SmilesParseError::CxZeroBondOutOfRange { bond, bond_count } => {
            SmartsWriteError::CxSourceBondOutOfRange { bond, bond_count }
        }
        other => SmartsWriteError::CanonicalTraversal(other),
    })
}

#[cfg(feature = "smiles-integration")]
fn query_cx_stereo_kind_order(kind: StereoGroupKind) -> u8 {
    match kind {
        StereoGroupKind::Absolute => 0,
        StereoGroupKind::Or => 1,
        StereoGroupKind::And => 2,
    }
}

#[cfg(feature = "smiles-integration")]
fn get_sorted_mapped_indexes(
    atom_ids: &[AtomId],
    rev_order: &[usize],
) -> Result<Vec<usize>, SmartsWriteError> {
    // RDKit❗✔️: std::vector<unsigned> getSortedMappedIndexes(
    // RDKit❗✔️:     const std::vector<unsigned int> &atomIds,
    // RDKit❗✔️:     const std::vector<unsigned> &revOrder) {
    // RDKit❗✔️:   std::vector<unsigned> res;
    // RDKit❗✔️:   res.reserve(atomIds.size());
    // RDKit❗✔️:   for (auto atomId : atomIds) {
    // RDKit❗✔️:     res.push_back(revOrder[atomId]);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::sort(res.begin(), res.end());
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // Map in incoming source order, then sort the mapped indexes. Duplicate
    // members and default-zero reverse entries are retained exactly.
    let mut result = Vec::with_capacity(atom_ids.len());
    for atom in atom_ids {
        let index =
            rev_order
                .get(atom.index())
                .ok_or(SmartsWriteError::CxSourceAtomOutOfRange {
                    atom: *atom,
                    atom_count: rev_order.len(),
                })?;
        result.push(*index);
    }
    result.sort_unstable();
    // Cost: one reserved O(n) result, O(n log n) integer sort, no set, query
    // clone, sparse filtering or extra reverse-vector reconstruction. Checked
    // indexing translates invalid source accesses structurally, never defaults.
    Ok(result)
}

#[cfg(feature = "smiles-integration")]
fn get_sorted_stereo_groups_and_indices(
    query: &QueryGraph,
    rev_order: &[usize],
    wedge_bonds: &cosmolkit_core::WedgeAssignments,
) -> Result<Vec<(StereoGroup, Vec<usize>)>, SmartsWriteError> {
    // RDKit❗🔝: getSortedStereoGroupsAndIndices(
    // RDKit❗🔝:     const ROMol &mol, const std::vector<unsigned int> &revOrder,
    // RDKit❗🔝:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗🔝:         &wedgeBonds) {
    // RDKit❗🔝:   using StGrpIdxPair = std::pair<StereoGroup, std::vector<unsigned>>;
    // RDKit❗🔝:
    // RDKit❗🔝:   auto &groups = mol.getStereoGroups();
    // RDKit❗🔝:
    // RDKit❗🔝:   std::vector<StGrpIdxPair> sortingGroups;
    // RDKit❗🔝:   sortingGroups.reserve(groups.size());
    // RDKit❗🔝:
    // RDKit❗🔝:   for (const auto &sg : groups) {
    // RDKit❗🔝:     std::vector<unsigned int> atomIds;
    // RDKit❗🔝:     Atropisomers::getAllAtomIdsForStereoGroup(mol, sg, atomIds, wedgeBonds);
    // RDKit❗🔝:     const auto newAtomIndexes = getSortedMappedIndexes(atomIds, revOrder);
    // RDKit❗🔝:     if (!newAtomIndexes.empty()) {
    // RDKit❗🔝:       sortingGroups.emplace_back(sg, newAtomIndexes);
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:
    // RDKit❗🔝:   // sort by 1) StereoGroup type; 2) StereoGroup id; 3) atom indexes
    // RDKit❗🔝:   std::sort(sortingGroups.begin(), sortingGroups.end(),
    // RDKit❗🔝:             [](const StGrpIdxPair &a, const StGrpIdxPair &b) {
    // RDKit❗🔝:               const auto &[sgA, idxsA] = a;
    // RDKit❗🔝:               const auto &[sgB, idxsB] = b;
    // RDKit❗🔝:               if (sgA.getGroupType() == sgB.getGroupType()) {
    // RDKit❗🔝:                 if (sgA.getWriteId() == sgB.getWriteId()) {
    // RDKit❗🔝:                   return idxsA < idxsB;
    // RDKit❗🔝:                 }
    // RDKit❗🔝:                 return sgA.getWriteId() < sgB.getWriteId();
    // RDKit❗🔝:               }
    // RDKit❗🔝:               return sgA.getGroupType() < sgB.getGroupType();
    // RDKit❗🔝:             });
    // RDKit❗🔝:
    // RDKit❗🔝:   std::vector<StereoGroup> sgs;
    // RDKit❗🔝:   std::vector<std::vector<unsigned>> sgAtomIdxs;
    // RDKit❗🔝:   sgs.reserve(sortingGroups.size());
    // RDKit❗🔝:   sgAtomIdxs.reserve(sortingGroups.size());
    // RDKit❗🔝:
    // RDKit❗🔝:   for (auto &&p : sortingGroups) {
    // RDKit❗🔝:     sgs.push_back(std::move(p.first));
    // RDKit❗🔝:     sgAtomIdxs.push_back(std::move(p.second));
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return {std::move(sgs), std::move(sgAtomIdxs)};
    // RDKit❗🔝: }
    let mut sorting_groups = Vec::with_capacity(query.stereo_groups().len());
    for group in query.stereo_groups() {
        let mut atom_ids = Vec::new();
        cosmolkit_core::collect_query_stereo_group_atom_ids_source(
            query,
            group,
            &mut atom_ids,
            wedge_bonds,
        )?;
        let new_atom_indexes = get_sorted_mapped_indexes(&atom_ids, rev_order)?;
        if !new_atom_indexes.is_empty() {
            sorting_groups.push((group.clone(), new_atom_indexes));
        }
    }
    sorting_groups.sort_unstable_by(|(a, indexes_a), (b, indexes_b)| {
        query_cx_stereo_kind_order(a.kind())
            .cmp(&query_cx_stereo_kind_order(b.kind()))
            .then_with(|| a.write_id().cmp(&b.write_id()))
            .then_with(|| indexes_a.cmp(indexes_b))
    });
    // Cost 🔝: the private paired return carries the exact same associated
    // source group/index rows, avoiding Native's two final split-vector
    // allocations and traversal. Moving mapped indexes also avoids Native's
    // copy of const newAtomIndexes into the sorting pair. Group members still
    // copy once, and the source comparator/unstable O(n log n) sort is retained.
    // Equal-key library sort order remains an explicit final source review
    // item; neither read IDs nor other source-absent fields break ties here.
    Ok(sorting_groups)
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_enhanced_stereo(
    query: &QueryGraph,
    atom_order: &[AtomId],
    wedge_bonds: &cosmolkit_core::WedgeAssignments,
) -> Result<String, SmartsWriteError> {
    // RDKit❗❌: std::string get_enhanced_stereo_block(
    // RDKit❗❌:     const ROMol &mol, const std::vector<unsigned int> &atomOrder,
    // RDKit❗❌:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗❌:         &wedgeBonds) {
    // RDKit❗❌:   if (mol.getStereoGroups().empty()) {
    // RDKit❗❌:     return "";
    // RDKit❗❌:   }
    // RDKit❗❌:   std::stringstream res;
    // RDKit❗❌:   // we need a map from original atom idx to output idx:
    // RDKit❗❌:   std::vector<unsigned int> revOrder(mol.getNumAtoms());
    // RDKit❗❌:   for (unsigned i = 0; i < atomOrder.size(); ++i) {
    // RDKit❗❌:     revOrder[atomOrder[i]] = i;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   auto [groups, groupsAtoms] =
    // RDKit❗❌:       getSortedStereoGroupsAndIndices(mol, revOrder, wedgeBonds);
    // RDKit❗❌:
    // RDKit❗❌:   assignStereoGroupIds(groups);
    // RDKit❗❌:
    // RDKit❗❌:   auto grpAtomsItr = groupsAtoms.begin();
    // RDKit❗❌:   for (auto sgItr = groups.begin(); sgItr != groups.end();
    // RDKit❗❌:        ++sgItr, ++grpAtomsItr) {
    // RDKit❗❌:     switch (sgItr->getGroupType()) {
    // RDKit❗❌:       case StereoGroupType::STEREO_ABSOLUTE:
    // RDKit❗❌:         res << "a:";
    // RDKit❗❌:         break;
    // RDKit❗❌:       case StereoGroupType::STEREO_OR:
    // RDKit❗❌:         res << "o" << sgItr->getWriteId() << ":";
    // RDKit❗❌:         break;
    // RDKit❗❌:       case StereoGroupType::STEREO_AND:
    // RDKit❗❌:         res << "&" << sgItr->getWriteId() << ":";
    // RDKit❗❌:         break;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     for (const auto &aid : *grpAtomsItr) {
    // RDKit❗❌:       res << aid << ",";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::string resStr = res.str();
    // RDKit❗❌:   if (!resStr.empty() && resStr.back() == ',') {
    // RDKit❗❌:     resStr.pop_back();
    // RDKit❗❌:   }
    // RDKit❗❌:   return resStr;
    // RDKit❗❌: }
    if query.stereo_groups().is_empty() {
        return Ok(String::new());
    }
    let mut reverse_order = vec![0; query.num_atoms()];
    for (index, atom) in atom_order.iter().enumerate() {
        *reverse_order
            .get_mut(atom.index())
            .ok_or(SmartsWriteError::CxSourceAtomOutOfRange {
                atom: *atom,
                atom_count: query.num_atoms(),
            })? = index;
    }
    let mut groups = get_sorted_stereo_groups_and_indices(query, &reverse_order, wedge_bonds)?;
    cosmolkit_smiles::assign_stereo_group_ids(&mut groups);
    let mut result = String::new();
    for (group, atom_indexes) in groups {
        match group.kind() {
            StereoGroupKind::Absolute => result.push_str("a:"),
            StereoGroupKind::Or => {
                write!(&mut result, "o{}:", group.write_id()).expect("String writing is infallible")
            }
            StereoGroupKind::And => {
                write!(&mut result, "&{}:", group.write_id()).expect("String writing is infallible")
            }
        }
        for atom in atom_indexes {
            write!(&mut result, "{atom},").expect("String writing is infallible");
        }
    }
    if result.ends_with(',') {
        result.pop();
    }
    // Every generated byte is source ASCII framing or unsigned numeric text.
    // Input group IDs/membership stay unchanged; only sorted local copies get
    // missing/duplicate write IDs. One trailing comma is removed once.
    // Cost ❌: Rust String lacks Native SSO on tiny outputs. On longer output,
    // direct one-buffer emission avoids Native res.str()'s final full copy and
    // the former Vec-of-per-number strings; all source scans/sorts are retained.
    Ok(result)
}

#[cfg(feature = "smiles-integration")]
fn query_cx_other_atom(bond: &Bond, atom: AtomId) -> Option<AtomId> {
    if bond.begin() == atom {
        Some(bond.end())
    } else if bond.end() == atom {
        Some(bond.begin())
    } else {
        None
    }
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_ring_bond_stereo(
    query: &QueryGraph,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> Result<String, SmartsWriteError> {
    // RDKit❗❌: std::string get_ringbond_cistrans_block(
    // RDKit❗❌:     const ROMol &mol, const std::vector<unsigned int> &atomOrder,
    // RDKit❗❌:     const std::vector<unsigned int> &bondOrder) {
    // RDKit❗❌:   if (!mol.getRingInfo()->isInitialized()) {
    // RDKit❗❌:     return "";
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   const auto rinfo = mol.getRingInfo();
    // RDKit❗❌:   std::string c = "", t = "", ctu = "";
    // RDKit❗❌:   for (unsigned int i = 0; i < bondOrder.size(); ++i) {
    // RDKit❗❌:     auto idx = bondOrder[i];
    // RDKit❗❌:     if (!rinfo->numBondRings(idx) ||
    // RDKit❗❌:         rinfo->minBondRingSize(idx) <
    // RDKit❗❌:             Chirality::minRingSizeForDoubleBondStereo) {
    // RDKit❗❌:       // we only do ring bonds of a minimum size
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     const auto bond = mol.getBondWithIdx(idx);
    // RDKit❗❌:     if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit❗❌:         bond->getBondType() != Bond::BondType::AROMATIC) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     Bond::BondStereo bstereo = bond->getStereo();
    // RDKit❗❌:     if (bstereo != Bond::BondStereo::STEREOANY &&
    // RDKit❗❌:         bstereo != Bond::BondStereo::STEREOCIS &&
    // RDKit❗❌:         bstereo != Bond::BondStereo::STEREOTRANS) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     auto label = std::to_string(i);
    // RDKit❗❌:
    // RDKit❗❌:     if (bstereo == Bond::BondStereo::STEREOANY) {
    // RDKit❗❌:       // this one's easy because we don't care about the atom order.
    // RDKit❗❌:       if (ctu.empty()) {
    // RDKit❗❌:         ctu += "ctu:";
    // RDKit❗❌:       } else {
    // RDKit❗❌:         ctu += ",";
    // RDKit❗❌:       }
    // RDKit❗❌:       ctu += label;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       Atom *begAtom = bond->getBeginAtom();
    // RDKit❗❌:       Atom *endAtom = bond->getEndAtom();
    // RDKit❗❌:       bool needSwap = false;
    // RDKit❗❌:       if (begAtom->getDegree() > 2) {
    // RDKit❗❌:         unsigned int o1 = atomOrder[bond->getStereoAtoms()[0]];
    // RDKit❗❌:         for (const auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit❗❌:           if (nbr == endAtom ||
    // RDKit❗❌:               nbr->getIdx() ==
    // RDKit❗❌:                   static_cast<unsigned>(bond->getStereoAtoms()[0])) {
    // RDKit❗❌:             continue;
    // RDKit❗❌:           }
    // RDKit❗❌:           if (atomOrder[nbr->getIdx() < o1]) {
    // RDKit❗❌:             // this neighbor came first, we need to swap:
    // RDKit❗❌:             needSwap = !needSwap;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (endAtom->getDegree() > 2) {
    // RDKit❗❌:         unsigned int o1 = atomOrder[bond->getStereoAtoms()[1]];
    // RDKit❗❌:         for (const auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit❗❌:           if (nbr == begAtom ||
    // RDKit❗❌:               nbr->getIdx() ==
    // RDKit❗❌:                   static_cast<unsigned>(bond->getStereoAtoms()[1])) {
    // RDKit❗❌:             continue;
    // RDKit❗❌:           }
    // RDKit❗❌:           if (atomOrder[nbr->getIdx() < o1]) {
    // RDKit❗❌:             // this neighbor came first, we need to swap:
    // RDKit❗❌:             needSwap = !needSwap;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (bstereo == Bond::BondStereo::STEREOCIS || needSwap) {
    // RDKit❗❌:         if (c.empty()) {
    // RDKit❗❌:           c += "c:";
    // RDKit❗❌:         } else {
    // RDKit❗❌:           c += ",";
    // RDKit❗❌:         }
    // RDKit❗❌:         c += label;
    // RDKit❗❌:       } else {
    // RDKit❗❌:         if (t.empty()) {
    // RDKit❗❌:           t += "t:";
    // RDKit❗❌:         } else {
    // RDKit❗❌:           t += ",";
    // RDKit❗❌:         }
    // RDKit❗❌:         t += label;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return c + t + ctu;
    // RDKit❗❌: }
    // RDKit❗✔️: constexpr unsigned int minRingSizeForDoubleBondStereo = 8;
    // The Native initialized guard precedes every order/cache access. No ring
    // finding, copied Bond rows, rebuilt adjacency or guessed fresh cache.
    if !query.source_ring_info().initialized {
        return Ok(String::new());
    }
    let rings = cosmolkit_core::RingInfo::from_source_snapshot(query.source_ring_info())
        .map_err(SmartsWriteError::CxRingInfo)?;
    let mut cis = String::new();
    let mut trans = String::new();
    let mut unknown = String::new();
    let order_at = |index: usize| {
        atom_order
            .get(index)
            .copied()
            .ok_or(SmartsWriteError::CxRingAtomOrderIndex {
                index,
                count: atom_order.len(),
            })
    };
    for (i, &id) in bond_order.iter().enumerate() {
        if rings.num_bond_rings(id) == 0 || rings.min_bond_ring_size(id) < 8 {
            continue;
        }
        let bond = query
            .bonds()
            .get(id.index())
            .ok_or(SmartsWriteError::CxSourceBondOutOfRange {
                bond: id,
                bond_count: query.num_bonds(),
            })?
            .bond();
        if !matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic) {
            continue;
        }
        let stereo = bond.stereo();
        if !matches!(
            stereo,
            BondStereo::Any | BondStereo::Cis | BondStereo::Trans
        ) {
            continue;
        }
        let (target, prefix) = if stereo == BondStereo::Any {
            (&mut unknown, "ctu:")
        } else {
            let mut need_swap = false;
            for (side, center, opposite) in
                [(0, bond.begin(), bond.end()), (1, bond.end(), bond.begin())]
            {
                let neighbors = &query.adjacency()[center.index()];
                if neighbors.len() > 2 {
                    let references = bond
                        .stereo_atoms()
                        .ok_or(SmartsWriteError::CxRingStereoReferenceMissing { bond: id, side })?;
                    let reference = references[side];
                    let o1 = order_at(reference.index())?.index();
                    for &(neighbor, _) in neighbors {
                        if neighbor == opposite.index() || neighbor == reference.index() {
                            continue;
                        }
                        // Literal pinned source bool subscript; do NOT replace
                        // this with reverse-order lookup or compare positions.
                        if order_at(usize::from(neighbor < o1))?.index() != 0 {
                            need_swap = !need_swap;
                        }
                    }
                }
            }
            if stereo == BondStereo::Cis || need_swap {
                (&mut cis, "c:")
            } else {
                (&mut trans, "t:")
            }
        };
        if target.is_empty() {
            target.push_str(prefix);
        } else {
            target.push(',');
        }
        write!(target, "{i}").expect("string formatting is infallible");
    }
    // Native returns c+t+ctu with no commas between the three groups.
    cis.push_str(&trans);
    cis.push_str(&unknown);
    // Cost ❌: actual cache import deep-copies represented vectors once versus
    // the source borrowed pointer; String/decimal output lacks Native SSO.
    // Source bond/member/neighbor/control order otherwise stays literal, with
    // three growing buffers and no cloned graph, cache rebuild or repair.
    Ok(cis)
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_link_nodes(
    query: &QueryGraph,
    atom_order: &[AtomId],
) -> Result<String, SmartsWriteError> {
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
    // The canonical complete parser and emitter own both carrier paths.
    // No query-to-Atom conversion, inferred bonds, partial parser or skipped
    // malformed order. Source non-strict warnings precede later errors.
    cosmolkit_smiles::write_query_cx_link_nodes_source(
        query,
        atom_order,
        &mut cosmolkit_smiles::emit_cx_link_node_warning_source,
    )
    .map_err(SmartsWriteError::CanonicalTraversal)
}

#[cfg(feature = "smiles-integration")]
fn query_cx_is_data_sgroup(group: &SubstanceGroup) -> Result<bool, SmartsWriteError> {
    // TYPE is the actual raw property or canonical constructor-kind fact.
    // Generic("DAT") therefore follows the same source branch as Data.
    Ok(query_cx_source_sgroup_type(group)?.as_bytes() == b"DAT")
}

#[cfg(feature = "smiles-integration")]
fn query_cx_sgroup_value(
    group: &SubstanceGroup,
    key: &str,
) -> Result<PropertyText, SmartsWriteError> {
    // BEGIN COMPLETE RDProps::getPropIfPresent
    // RDKit✔️❌: bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit✔️❌:     return d_props.getValIfPresent(key, res);
    // RDKit✔️❌:   }
    // END COMPLETE RDProps::getPropIfPresent

    // Scalar source getProp<string> uses canonical RDValue conversion. This
    // helper is never used to reinterpret DATAFIELDS as a scalar string.
    if let Some(value) = group.props().get(key.as_bytes()) {
        return cosmolkit_core::property_value_to_string(value).map_err(SmartsWriteError::from);
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

#[cfg(feature = "smiles-integration")]
fn query_cx_data_sgroup_values(
    group: &SubstanceGroup,
) -> Result<&[PropertyText], SmartsWriteError> {
    // BEGIN COMPLETE RDProps::getPropIfPresent
    // RDKit✔️❌: bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit✔️❌:     return d_props.getValIfPresent(key, res);
    // RDKit✔️❌:   }
    // END COMPLETE RDProps::getPropIfPresent

    // RDKit❗❌: template <>
    // RDKit❗❌: inline std::vector<std::string> rdvalue_cast<std::vector<std::string>>(
    // RDKit❗❌:     RDValue_cast_t v) {
    // RDKit❗❌:   if (rdvalue_is<std::vector<std::string>>(v)) {
    // RDKit❗❌:     return *v.ptrCast<std::vector<std::string>>();
    // RDKit❗❌:   }
    // RDKit❗❌:   throw std::bad_any_cast();
    // RDKit❗❌: }
    // Source DATAFIELDS is vector<string>: exact tag, bytes and element order.
    // A present scalar is a native cast failure, never split or inferred.
    if let Some(value) = group.props().get(b"DATAFIELDS".as_slice()) {
        return value
            .as_string_vector()
            .map_err(SmartsWriteError::PropertyValue);
    }
    Ok(group
        .data()
        .map(|data| data.values.as_slice())
        .unwrap_or_else(|| group.data_fields()))
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_data_sgroups(
    query: &QueryGraph,
    state: &mut QueryCxSgroupState,
    atom_order: &[AtomId],
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string get_sgroup_data_block(const ROMol &mol,
    // RDKit❗❌:                                   const std::vector<unsigned int> &atomOrder) {
    // RDKit❗❌:   const auto &sgs = getSubstanceGroups(mol);
    // RDKit❗❌:   if (sgs.empty()) {
    // RDKit❗❌:     return "";
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int sgroupOutputIndex = 0;
    // RDKit❗❌:   mol.getPropIfPresent("_cxsmilesOutputIndex", sgroupOutputIndex);
    // RDKit❗❌:
    // RDKit❗❌:   std::stringstream res;
    // RDKit❗❌:   // we need a map from original atom idx to output idx:
    // RDKit❗❌:   std::vector<unsigned int> revOrder(mol.getNumAtoms());
    // RDKit❗❌:   for (unsigned i = 0; i < atomOrder.size(); ++i) {
    // RDKit❗❌:     revOrder[atomOrder[i]] = i;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto &sg : sgs) {
    // RDKit❗❌:     if (sg.hasProp("TYPE") && sg.getProp<std::string>("TYPE") == "DAT") {
    // RDKit❗❌:       sg.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);
    // RDKit❗❌:       ++sgroupOutputIndex;
    // RDKit❗❌:
    // RDKit❗❌:       res << "SgD:";
    // RDKit❗❌:       // we don't attempt to canonicalize the atom order because the user
    // RDKit❗❌:       // may ascribe some significance to the ordering of the atoms
    // RDKit❗❌:       for (const auto oaid : sg.getAtoms()) {
    // RDKit❗❌:         res << revOrder[oaid] << ",";
    // RDKit❗❌:       }
    // RDKit❗❌:       // remove the extra ",":
    // RDKit❗❌:       res.seekp(-1, res.cur);
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       std::string prop;
    // RDKit❗❌:       if (sg.getPropIfPresent("FIELDNAME", prop) && !prop.empty()) {
    // RDKit❗❌:         res << prop;
    // RDKit❗❌:       }
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       std::vector<std::string> vprop;
    // RDKit❗❌:       if (sg.getPropIfPresent("DATAFIELDS", vprop) && !vprop.empty()) {
    // RDKit❗❌:         for (const auto &pv : vprop) {
    // RDKit❗❌:           res << pv << ",";
    // RDKit❗❌:         }
    // RDKit❗❌:         // remove the extra ",":
    // RDKit❗❌:         res.seekp(-1, res.cur);
    // RDKit❗❌:       }
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       if (sg.getPropIfPresent("QUERYOP", prop) && !prop.empty()) {
    // RDKit❗❌:         res << prop;
    // RDKit❗❌:       }
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       if (sg.getPropIfPresent("FIELDINFO", prop) && !prop.empty()) {
    // RDKit❗❌:         res << prop;
    // RDKit❗❌:       }
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       if (sg.getPropIfPresent("FIELDTAG", prop) && !prop.empty()) {
    // RDKit❗❌:         res << prop;
    // RDKit❗❌:       }
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       // FIX: do something about the coordinates
    // RDKit❗❌:       res << ",";  // only add a comma if we wrote something
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::string resStr = res.str();
    // RDKit❗❌:   if (!resStr.empty() && resStr.back() == ',') {
    // RDKit❗❌:     resStr.pop_back();
    // RDKit❗❌:   }
    // RDKit❗❌:   mol.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);
    // RDKit❗❌:
    // RDKit❗❌:   return resStr;
    // RDKit❗❌: }
    if state.groups.is_empty() {
        return Ok(PropertyText::new());
    }
    let mut output_index = match state.properties.prop("_cxsmilesOutputIndex") {
        Some(value) => cosmolkit_core::property_value_to_uint(value).map_err(|source| {
            SmartsWriteError::CxMoleculePropertyUInt {
                property: "_cxsmilesOutputIndex",
                source,
            }
        })?,
        None => 0,
    };
    let mut reverse_order = vec![0; query.num_atoms()];
    for (i, atom) in atom_order.iter().enumerate() {
        *reverse_order
            .get_mut(atom.index())
            .ok_or(SmartsWriteError::CxSourceAtomOutOfRange {
                atom: *atom,
                atom_count: query.num_atoms(),
            })? = i;
    }
    let mut result = PropertyText::new();
    for group in &mut state.groups {
        if !query_cx_is_data_sgroup(group)? {
            continue;
        }
        group
            .set_prop(
                "_cxsmilesOutputIndex",
                cosmolkit_model::PropertyValue::UInt(output_index),
            )
            .map_err(|source| SmartsWriteError::CxSgroupPropertyWrite {
                group: group.id(),
                property: "_cxsmilesOutputIndex",
                source,
            })?;
        output_index = output_index.wrapping_add(1);
        if !result.is_empty() {
            result.push_byte(b',');
        }
        result.extend_bytes(b"SgD:");
        // Member order is source user metadata, not canonicalized or deduped.
        for (i, atom) in group.atoms().iter().enumerate() {
            if i != 0 {
                result.push_byte(b',');
            }
            let index = reverse_order.get(atom.index()).ok_or(
                SmartsWriteError::CxSourceAtomOutOfRange {
                    atom: *atom,
                    atom_count: reverse_order.len(),
                },
            )?;
            write!(&mut result, "{index}").expect("byte formatting is infallible");
        }
        // Empty members make source seekp(-1) overwrite the prefix colon.
        if !group.atoms().is_empty() {
            result.push_byte(b':');
        }
        result.extend_bytes(query_cx_sgroup_value(group, "FIELDNAME")?.as_bytes());
        result.push_byte(b':');
        let values = query_cx_data_sgroup_values(group)?;
        for (i, value) in values.iter().enumerate() {
            if i != 0 {
                result.push_byte(b',');
            }
            result.extend_bytes(value.as_bytes());
        }
        result.push_byte(b':');
        for key in ["QUERYOP", "FIELDINFO", "FIELDTAG"] {
            result.extend_bytes(query_cx_sgroup_value(group, key)?.as_bytes());
            result.push_byte(b':');
        }
    }
    // This source molecule write is after the complete loop only. A later
    // bad DATAFIELDS cast retains prior/current group writes but not this write.
    state.properties.set_prop(
        "_cxsmilesOutputIndex",
        cosmolkit_model::PropertyValue::UInt(output_index),
    )?;
    // Cost ❌: detached canonical group/property copies and no short-string
    // storage. Source reverse scan and one counted byte buffer are retained;
    // numeric formatting adds no per-atom owned String, vectors stay borrowed.
    Ok(result)
}

#[cfg(feature = "smiles-integration")]
fn query_cx_source_sgroup_type(group: &SubstanceGroup) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: SubstanceGroup::SubstanceGroup(ROMol *owning_mol, const std::string &type)
    // RDKit❗❌:     : RDProps(), dp_mol(owning_mol) {
    // RDKit❗❌:   PRECONDITION(owning_mol, "supplied owning molecule is bad");
    // RDKit❗❌:
    // RDKit❗❌:   // TYPE is required to be set , as other properties will depend on it.
    // RDKit❗❌:   setProp<std::string>("TYPE", type);
    // RDKit❗❌: }
    // Canonical kind is the detached constructor TYPE fact. An actual raw
    // TYPE overrides that projected fact and keeps source string conversion.
    if let Some(value) = group.props().get(b"TYPE".as_slice()) {
        return Ok(cosmolkit_core::property_value_to_string(value)?);
    }
    Ok(match group.kind() {
        SubstanceGroupKind::Data => "DAT".into(),
        SubstanceGroupKind::Superatom => "SUP".into(),
        SubstanceGroupKind::MultipleGroup => "MUL".into(),
        SubstanceGroupKind::StructuralRepeatUnit => "SRU".into(),
        SubstanceGroupKind::Monomer => "MON".into(),
        SubstanceGroupKind::Copolymer => "COP".into(),
        SubstanceGroupKind::Crosslink => "CRO".into(),
        SubstanceGroupKind::Graft => "GRA".into(),
        SubstanceGroupKind::Modification => "MOD".into(),
        SubstanceGroupKind::Mer => "MER".into(),
        SubstanceGroupKind::AnyPolymer => "ANY".into(),
        SubstanceGroupKind::MixtureComponent => "COM".into(),
        SubstanceGroupKind::Mixture => "MIX".into(),
        SubstanceGroupKind::Formulation => "FOR".into(),
        SubstanceGroupKind::Generic(value) => value.clone(),
    })
}

#[cfg(feature = "smiles-integration")]
fn query_cx_connection_text(group: &SubstanceGroup) -> Result<PropertyText, SmartsWriteError> {
    // Boost❗❌:         template<typename WritableRangeT>
    // Boost❗❌:         inline void to_lower(
    // Boost❗❌:             WritableRangeT& Input,
    // Boost❗❌:             const std::locale& Loc=std::locale())
    // Boost❗❌:         {
    // Boost❗❌:             ::boost::algorithm::detail::transform_range(
    // Boost❗❌:                 ::boost::as_literal(Input),
    // Boost❗❌:                 ::boost::algorithm::detail::to_lowerF<
    // Boost❗❌:                     typename range_value<WritableRangeT>::type >(Loc));
    // Boost❗❌:         }
    // Boost❗❌:             template<typename CharT>
    // Boost❗❌:             struct to_lowerF
    // Boost❗❌:             {
    // Boost❗❌:                 typedef CharT argument_type;
    // Boost❗❌:                 typedef CharT result_type;
    // Boost❗❌:                 // Constructor
    // Boost❗❌:                 to_lowerF( const std::locale& Loc ) : m_Loc( &Loc ) {}
    // Boost❗❌:
    // Boost❗❌:                 // Operation
    // Boost❗❌:                 CharT operator ()( CharT Ch ) const
    // Boost❗❌:                 {
    // Boost❗❌:                     #if defined(BOOST_BORLANDC) && (BOOST_BORLANDC >= 0x560) && (BOOST_BORLANDC <= 0x564) && !defined(_USE_OLD_RW_STL)
    // Boost❗❌:                         return std::tolower( static_cast<typename boost::make_unsigned <CharT>::type> ( Ch ));
    // Boost❗❌:                     #else
    // Boost❗❌:                         return std::tolower<CharT>( Ch, *m_Loc );
    // Boost❗❌:                     #endif
    // Boost❗❌:                 }
    // Boost❗❌:             private:
    // Boost❗❌:                 const std::locale* m_Loc;
    // Boost❗❌:             };
    // Boost❗❌:             template<typename RangeT, typename FunctorT>
    // Boost❗❌:             void transform_range(
    // Boost❗❌:                 const RangeT& Input,
    // Boost❗❌:                 FunctorT Functor)
    // Boost❗❌:             {
    // Boost❗❌:                 std::transform(
    // Boost❗❌:                     ::boost::begin(Input),
    // Boost❗❌:                     ::boost::end(Input),
    // Boost❗❌:                     ::boost::begin(Input),
    // Boost❗❌:                     Functor);
    // Boost❗❌:             }
    let text = if let Some(value) = group.props().get(b"CONNECT".as_slice()) {
        cosmolkit_core::property_value_to_string(value)?
    } else {
        match group.connection() {
            Some(SGroupConnection::HeadToHead) => "HH".into(),
            Some(SGroupConnection::HeadToTail) => "HT".into(),
            Some(SGroupConnection::Either) => "EU".into(),
            Some(SGroupConnection::Unknown(value)) => value.clone(),
            None => PropertyText::new(),
        }
    };
    // Source boost::algorithm::to_lower in the pinned C-locale profile.
    Ok(PropertyText::from(text.as_bytes().to_ascii_lowercase()))
}

#[cfg(feature = "smiles-integration")]
fn query_cx_sgroup_crossings<'a>(
    group: &'a SubstanceGroup,
    property: &'static str,
    typed: &'a [BondId],
) -> Result<&'a [BondId], SmartsWriteError> {
    // RDKit❗✔️: template <>
    // RDKit❗✔️: inline std::vector<unsigned int> rdvalue_cast<std::vector<unsigned int>>(
    // RDKit❗✔️:     RDValue_cast_t v) {
    // RDKit❗✔️:   if (rdvalue_is<std::vector<unsigned int>>(v)) {
    // RDKit❗✔️:     return *v.ptrCast<std::vector<unsigned int>>();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   throw std::bad_any_cast();
    // RDKit❗✔️: }
    // Source requires vector<unsigned int>, represented by the actual typed
    // BondId field. Every modeled generic RDValue tag is a different type,
    // including IntVector: a present raw tag must fail, never be reinterpreted.
    if let Some(value) = group.props().get(property.as_bytes()) {
        return Err(SmartsWriteError::CxSgroupVectorCast {
            group: group.id(),
            property,
            actual: value.kind(),
            expected: "vector<unsigned int>",
        });
    }
    Ok(typed)
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_polymer_sgroups(
    query: &QueryGraph,
    state: &mut QueryCxSgroupState,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string get_sgroup_polymer_block(
    // RDKit❗❌:     const ROMol &mol, const std::vector<unsigned int> &atomOrder,
    // RDKit❗❌:     const std::vector<unsigned int> &bondOrder) {
    // RDKit❗❌:   const auto &sgs = getSubstanceGroups(mol);
    // RDKit❗❌:   if (sgs.empty()) {
    // RDKit❗❌:     return "";
    // RDKit❗❌:   }
    // RDKit❗❌:   unsigned int sgroupOutputIndex = 0;
    // RDKit❗❌:   mol.getPropIfPresent("_cxsmilesOutputIndex", sgroupOutputIndex);
    // RDKit❗❌:   std::stringstream res;
    // RDKit❗❌:   // we need a map from original atom idx to output idx:
    // RDKit❗❌:   std::vector<unsigned int> revAtomOrder(mol.getNumAtoms());
    // RDKit❗❌:   for (unsigned i = 0; i < atomOrder.size(); ++i) {
    // RDKit❗❌:     revAtomOrder[atomOrder[i]] = i;
    // RDKit❗❌:   }
    // RDKit❗❌:   // we need a map from original bond idx to output idx:
    // RDKit❗❌:   std::vector<unsigned int> revBondOrder(mol.getNumBonds());
    // RDKit❗❌:   for (unsigned i = 0; i < bondOrder.size(); ++i) {
    // RDKit❗❌:     revBondOrder[bondOrder[i]] = i;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::map<std::string, std::string> reverseTypemap;
    // RDKit❗❌:   for (const auto &pr : SmilesParseOps::sgroupTypemap) {
    // RDKit❗❌:     if (reverseTypemap.find(pr.second) == reverseTypemap.end()) {
    // RDKit❗❌:       reverseTypemap[pr.second] = pr.first;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto &sg : sgs) {
    // RDKit❗❌:     std::string typ;
    // RDKit❗❌:     if (sg.getPropIfPresent("TYPE", typ) &&
    // RDKit❗❌:         reverseTypemap.find(typ) != reverseTypemap.end()) {
    // RDKit❗❌:       sg.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);
    // RDKit❗❌:       ++sgroupOutputIndex;
    // RDKit❗❌:
    // RDKit❗❌:       res << "Sg:";
    // RDKit❗❌:       std::string subtype;
    // RDKit❗❌:       if (typ == "COP" && sg.getPropIfPresent("SUBTYPE", subtype)) {
    // RDKit❗❌:         if (subtype == "ALT") {
    // RDKit❗❌:           res << "alt";
    // RDKit❗❌:         } else if (subtype == "RAN") {
    // RDKit❗❌:           res << "ran";
    // RDKit❗❌:         } else if (subtype == "BLO") {
    // RDKit❗❌:           res << "blk";
    // RDKit❗❌:         } else {
    // RDKit❗❌:           res << reverseTypemap["COP"];
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res << reverseTypemap[typ];
    // RDKit❗❌:       }
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       for (const auto oaid : sg.getAtoms()) {
    // RDKit❗❌:         res << revAtomOrder[oaid] << ",";
    // RDKit❗❌:       }
    // RDKit❗❌:       // remove the extra ",":
    // RDKit❗❌:       res.seekp(-1, res.cur);
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       std::string label;
    // RDKit❗❌:       if (sg.getPropIfPresent("LABEL", label)) {
    // RDKit❗❌:         res << label;
    // RDKit❗❌:       }
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       std::string connect;
    // RDKit❗❌:       if (sg.getPropIfPresent("CONNECT", connect)) {
    // RDKit❗❌:         boost::algorithm::to_lower(connect);
    // RDKit❗❌:         res << connect;
    // RDKit❗❌:       }
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       std::vector<unsigned int> headCrossings;
    // RDKit❗❌:       if (sg.getPropIfPresent("XBHEAD", headCrossings) &&
    // RDKit❗❌:           headCrossings.size() > 1) {
    // RDKit❗❌:         for (auto v : headCrossings) {
    // RDKit❗❌:           res << bondOrder[v] << ",";
    // RDKit❗❌:         }
    // RDKit❗❌:         // remove the extra ",":
    // RDKit❗❌:         res.seekp(-1, res.cur);
    // RDKit❗❌:       }
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       std::vector<unsigned int> tailCrossings;
    // RDKit❗❌:       if (sg.getPropIfPresent("XBCORR", tailCrossings) &&
    // RDKit❗❌:           tailCrossings.size() > 2) {
    // RDKit❗❌:         for (unsigned int i = 1; i < tailCrossings.size(); i += 2) {
    // RDKit❗❌:           res << bondOrder[tailCrossings[i]] << ",";
    // RDKit❗❌:         }
    // RDKit❗❌:         // remove the extra ",":
    // RDKit❗❌:         res.seekp(-1, res.cur);
    // RDKit❗❌:       }
    // RDKit❗❌:       res << ":";
    // RDKit❗❌:       res << ",";  // only add a comma if we wrote something
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::string resStr = res.str();
    // RDKit❗❌:   while (!resStr.empty() && resStr.back() == ',') {
    // RDKit❗❌:     resStr.pop_back();
    // RDKit❗❌:   }
    // RDKit❗❌:   mol.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);
    // RDKit❗❌:
    // RDKit❗❌:   return resStr;
    // RDKit❗❌: }
    // RDKit❗✔️: const std::map<std::string, std::string> sgroupTypemap = {
    // RDKit❗✔️:     {"n", "SRU"},   {"mon", "MON"}, {"mer", "MER"}, {"co", "COP"},
    // RDKit❗✔️:     {"xl", "CRO"},  {"mod", "MOD"}, {"mix", "MIX"}, {"f", "FOR"},
    // RDKit❗✔️:     {"any", "ANY"}, {"gen", "GEN"}, {"c", "COM"},   {"grf", "GRA"},
    // RDKit❗✔️:     {"alt", "COP"}, {"ran", "COP"}, {"blk", "COP"}};
    if state.groups.is_empty() {
        return Ok(PropertyText::new());
    }
    let mut output_index = match state.properties.prop("_cxsmilesOutputIndex") {
        Some(value) => cosmolkit_core::property_value_to_uint(value).map_err(|source| {
            SmartsWriteError::CxMoleculePropertyUInt {
                property: "_cxsmilesOutputIndex",
                source,
            }
        })?,
        None => 0,
    };
    let mut reverse_atoms = vec![0; query.num_atoms()];
    for (i, atom) in atom_order.iter().enumerate() {
        *reverse_atoms
            .get_mut(atom.index())
            .ok_or(SmartsWriteError::CxSourceAtomOutOfRange {
                atom: *atom,
                atom_count: query.num_atoms(),
            })? = i;
    }
    let mut reverse_bonds = vec![0; query.num_bonds()];
    for (i, bond) in bond_order.iter().enumerate() {
        *reverse_bonds
            .get_mut(bond.index())
            .ok_or(SmartsWriteError::CxSourceBondOutOfRange {
                bond: *bond,
                bond_count: query.num_bonds(),
            })? = i;
    }
    // Exact lexical traversal of the pinned static std::map above. The first
    // reverse entry wins (COP -> alt), including source duplicate values.
    const TYPEMAP: [(&str, &str); 15] = [
        ("alt", "COP"),
        ("any", "ANY"),
        ("blk", "COP"),
        ("c", "COM"),
        ("co", "COP"),
        ("f", "FOR"),
        ("gen", "GEN"),
        ("grf", "GRA"),
        ("mer", "MER"),
        ("mix", "MIX"),
        ("mod", "MOD"),
        ("mon", "MON"),
        ("n", "SRU"),
        ("ran", "COP"),
        ("xl", "CRO"),
    ];
    let mut reverse_typemap = BTreeMap::new();
    for (code, typ) in TYPEMAP {
        reverse_typemap.entry(typ.as_bytes()).or_insert(code);
    }
    let mut result = PropertyText::new();
    for group in &mut state.groups {
        let typ = query_cx_source_sgroup_type(group)?;
        let Some(default_code) = reverse_typemap.get(typ.as_bytes()) else {
            continue;
        };
        group
            .set_prop(
                "_cxsmilesOutputIndex",
                cosmolkit_model::PropertyValue::UInt(output_index),
            )
            .map_err(|source| SmartsWriteError::CxSgroupPropertyWrite {
                group: group.id(),
                property: "_cxsmilesOutputIndex",
                source,
            })?;
        output_index = output_index.wrapping_add(1);
        if !result.is_empty() {
            result.push_byte(b',');
        }
        result.extend_bytes(b"Sg:");
        let raw_subtype;
        let code = if typ.as_bytes() == b"COP" {
            let subtype = if let Some(value) = group.props().get(b"SUBTYPE".as_slice()) {
                raw_subtype = cosmolkit_core::property_value_to_string(value)?;
                Some(raw_subtype.as_bytes())
            } else {
                group.subtype().map(PropertyText::as_bytes)
            };
            match subtype {
                Some(b"ALT") => "alt",
                Some(b"RAN") => "ran",
                Some(b"BLO") => "blk",
                _ => default_code,
            }
        } else {
            default_code
        };
        result.extend_bytes(code.as_bytes());
        result.push_byte(b':');
        for (i, atom) in group.atoms().iter().enumerate() {
            if i != 0 {
                result.push_byte(b',');
            }
            let mapped = reverse_atoms.get(atom.index()).ok_or(
                SmartsWriteError::CxSourceAtomOutOfRange {
                    atom: *atom,
                    atom_count: reverse_atoms.len(),
                },
            )?;
            write!(&mut result, "{mapped}").expect("byte formatting is infallible");
        }
        // With empty members Native seekp(-1) overwrites the prefix colon.
        if !group.atoms().is_empty() {
            result.push_byte(b':');
        }
        let raw_label;
        let label = if let Some(value) = group.props().get(b"LABEL".as_slice()) {
            raw_label = cosmolkit_core::property_value_to_string(value)?;
            Some(&raw_label)
        } else {
            group.label()
        };
        if let Some(label) = label {
            result.extend_bytes(label.as_bytes());
        }
        result.push_byte(b':');
        result.extend_bytes(query_cx_connection_text(group)?.as_bytes());
        result.push_byte(b':');
        let head = query_cx_sgroup_crossings(group, "XBHEAD", group.head_crossing_bonds())?;
        if head.len() > 1 {
            for (i, crossing) in head.iter().enumerate() {
                if i != 0 {
                    result.push_byte(b',');
                }
                let bond = bond_order.get(crossing.index()).ok_or(
                    SmartsWriteError::CxSourceBondOutOfRange {
                        bond: *crossing,
                        bond_count: bond_order.len(),
                    },
                )?;
                write!(&mut result, "{}", bond.index()).expect("byte formatting is infallible");
            }
        }
        result.push_byte(b':');
        let tail =
            query_cx_sgroup_crossings(group, "XBCORR", group.crossing_bond_correspondence())?;
        if tail.len() > 2 {
            for (i, crossing) in tail.iter().skip(1).step_by(2).enumerate() {
                if i != 0 {
                    result.push_byte(b',');
                }
                let bond = bond_order.get(crossing.index()).ok_or(
                    SmartsWriteError::CxSourceBondOutOfRange {
                        bond: *crossing,
                        bond_count: bond_order.len(),
                    },
                )?;
                write!(&mut result, "{}", bond.index()).expect("byte formatting is infallible");
            }
        }
        result.push_byte(b':');
    }
    state.properties.set_prop(
        "_cxsmilesOutputIndex",
        cosmolkit_model::PropertyValue::UInt(output_index),
    )?;
    // Cost ❌: detached group/property scratch copies and PropertyText lack of
    // short string storage. One counted output buffer, borrowed crossing rows,
    // source-shaped reverse vectors/maps; no per-field string-vector buffering.
    Ok(result)
}

#[cfg(feature = "smiles-integration")]
fn query_cx_sgroup_index(group: &SubstanceGroup) -> Result<u32, SmartsWriteError> {
    match group.props().get(b"index".as_slice()) {
        Some(value) => cosmolkit_core::property_value_to_uint(value).map_err(|source| {
            SmartsWriteError::CxSgroupPropertyUInt {
                group: group.id(),
                property: "index",
                source,
            }
        }),
        None => Ok(group.id().index() as u32),
    }
}

// Transient writer-owned copies of canonical detached values. This is the
// source mutable-property scratch state for immutable public query writers,
// not a second graph/AST or persistent molecule representation. The data and
// polymer writers share these exact group rows and molecule property records.
#[cfg(feature = "smiles-integration")]
struct QueryCxSgroupState {
    groups: Vec<SubstanceGroup>,
    properties: cosmolkit_model::MoleculeProperties,
}

#[cfg(feature = "smiles-integration")]
impl QueryCxSgroupState {
    fn from_query(query: &QueryGraph) -> Self {
        Self {
            groups: query_substance_groups(query).to_vec(),
            properties: query.source_molecule_properties(),
        }
    }
}

#[cfg(feature = "smiles-integration")]
fn query_cx_sgroup_uint(
    group: &SubstanceGroup,
    property: &'static str,
    value: &cosmolkit_model::PropertyValue,
) -> Result<u32, SmartsWriteError> {
    cosmolkit_core::property_value_to_uint(value).map_err(|source| {
        SmartsWriteError::CxSgroupPropertyUInt {
            group: group.id(),
            property,
            source,
        }
    })
}

#[cfg(feature = "smiles-integration")]
fn write_query_cx_sgroup_hierarchy(
    state: &mut QueryCxSgroupState,
) -> Result<String, SmartsWriteError> {
    // RDKit❗❌: std::string get_sgroup_hierarchy_block(const ROMol &mol) {
    // RDKit❗❌:   const auto &sgs = getSubstanceGroups(mol);
    // RDKit❗❌:   if (sgs.empty()) {
    // RDKit❗❌:     return "";
    // RDKit❗❌:   }
    // RDKit❗❌:   std::stringstream res;
    // RDKit❗❌:   // we need a map from sgroup index to output index;
    // RDKit❗❌:   std::map<unsigned int, unsigned int> sgroupOrder;
    // RDKit❗❌:   bool parentPresent = false;
    // RDKit❗❌:   for (const auto &sg : sgs) {
    // RDKit❗❌:     if (sg.hasProp("_cxsmilesOutputIndex")) {
    // RDKit❗❌:       unsigned int sgidx = sg.getIndexInMol();
    // RDKit❗❌:       sg.getPropIfPresent("index", sgidx);
    // RDKit❗❌:       sgroupOrder[sgidx] = sg.getProp<unsigned int>("_cxsmilesOutputIndex");
    // RDKit❗❌:       sg.clearProp("_cxsmilesOutputIndex");
    // RDKit❗❌:     }
    // RDKit❗❌:     if (sg.hasProp("PARENT")) {
    // RDKit❗❌:       parentPresent = true;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (parentPresent) {
    // RDKit❗❌:     // now loop over them and add the information
    // RDKit❗❌:     std::map<unsigned int, std::vector<unsigned int>> accum;
    // RDKit❗❌:     for (const auto &sg : sgs) {
    // RDKit❗❌:       unsigned pidx;
    // RDKit❗❌:       if (sg.getPropIfPresent("PARENT", pidx) &&
    // RDKit❗❌:           sgroupOrder.find(pidx) != sgroupOrder.end()) {
    // RDKit❗❌:         unsigned int sgidx = sg.getIndexInMol();
    // RDKit❗❌:         sg.getPropIfPresent("index", sgidx);
    // RDKit❗❌:         if (sgroupOrder.find(sgidx) != sgroupOrder.end()) {
    // RDKit❗❌:           accum[sgroupOrder[pidx]].push_back(sgroupOrder[sgidx]);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!accum.empty()) {
    // RDKit❗❌:       res << "SgH:";
    // RDKit❗❌:       for (const auto &pr : accum) {
    // RDKit❗❌:         res << pr.first << ":";
    // RDKit❗❌:         for (auto v : pr.second) {
    // RDKit❗❌:           res << v << ".";
    // RDKit❗❌:         }
    // RDKit❗❌:         // remove the extra ".":
    // RDKit❗❌:         res.seekp(-1, res.cur);
    // RDKit❗❌:         res << ",";
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     std::string resStr = res.str();
    // RDKit❗❌:     while (!resStr.empty() && resStr.back() == ',') {
    // RDKit❗❌:       resStr.pop_back();
    // RDKit❗❌:     }
    // RDKit❗❌:     return resStr;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     return "";
    // RDKit❗❌:   }
    // RDKit❗❌: }
    let groups = &mut state.groups;
    if groups.is_empty() {
        return Ok(String::new());
    }
    let mut sgroup_order = BTreeMap::new();
    let mut parent_present = false;
    for group in groups.iter_mut() {
        if let Some(value) = group.props().get(b"_cxsmilesOutputIndex".as_slice()) {
            let source_index = query_cx_sgroup_index(group)?;
            let output_index = query_cx_sgroup_uint(group, "_cxsmilesOutputIndex", value)?;
            sgroup_order.insert(source_index, output_index);
            group.clear_prop("_cxsmilesOutputIndex").map_err(|source| {
                SmartsWriteError::CxSgroupPropertyWrite {
                    group: group.id(),
                    property: "_cxsmilesOutputIndex",
                    source,
                }
            })?;
        }
        parent_present |=
            group.parent().is_some() || group.props().contains_key(b"PARENT".as_slice());
    }
    if !parent_present {
        return Ok(String::new());
    }
    let mut accum = BTreeMap::<u32, Vec<u32>>::new();
    for group in groups.iter() {
        // Source reads PARENT before looking up the parent and only then
        // reads the child index. A nonemitted child cannot hide a bad PARENT.
        let parent = if let Some(value) = group.props().get(b"PARENT".as_slice()) {
            Some(query_cx_sgroup_uint(group, "PARENT", value)?)
        } else if let Some(parent) = group.parent() {
            let parent_group = groups.get(parent.index()).ok_or_else(|| {
                SmartsWriteError::InvalidGraph(format!(
                    "SGroup {} parent {} is outside {} groups",
                    group.id().index(),
                    parent.index(),
                    groups.len()
                ))
            })?;
            // Existing typed relation projects to its source external index.
            Some(query_cx_sgroup_index(parent_group)?)
        } else {
            None
        };
        if let Some(parent_output) = parent.and_then(|index| sgroup_order.get(&index).copied()) {
            let child_source_index = query_cx_sgroup_index(group)?;
            if let Some(child_output) = sgroup_order.get(&child_source_index).copied() {
                accum.entry(parent_output).or_default().push(child_output);
            }
        }
    }
    let mut result = String::new();
    if !accum.is_empty() {
        result.push_str("SgH:");
        for (parent, children) in accum {
            write!(&mut result, "{parent}:").expect("String writing is infallible");
            for child in children {
                write!(&mut result, "{child}.").expect("String writing is infallible");
            }
            // Every accum vector was created by a successful push, so source
            // seekp(-1) always overwrites the final dot with this comma.
            result.pop();
            result.push(',');
        }
    }
    while result.ends_with(',') {
        result.pop();
    }
    // Cost ❌: immutable query projection owns group/property copies, and
    // short String output lacks Native SSO. Map ordering/scans are source
    // shaped; canonical physical group IDs avoid Native linear pointer search.
    Ok(result)
}

#[cfg(feature = "smiles-integration")]
pub(crate) fn write_query_cx_extensions(
    query: &mut QueryGraph,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    fields: cosmolkit_smiles::CxSmilesFields,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string getCXExtensions(const ROMol &mol, std::uint32_t flags) {
    // RDKit❗❌:   std::string res = "|";
    // RDKit❗❌:   const std::vector<unsigned int> &atomOrder =
    // RDKit❗❌:       mol.getProp<std::vector<unsigned int>>(
    // RDKit❗❌:           common_properties::_smilesAtomOutputOrder);
    // RDKit❗❌:   const std::vector<unsigned int> &bondOrder =
    // RDKit❗❌:       mol.getProp<std::vector<unsigned int>>(
    // RDKit❗❌:           common_properties::_smilesBondOutputOrder);
    // RDKit❗❌:
    // RDKit❗❌:   bool needLabels = false;
    // RDKit❗❌:   bool needValues = false;
    // RDKit❗❌:   for (auto idx : atomOrder) {
    // RDKit❗❌:     const auto at = mol.getAtomWithIdx(idx);
    // RDKit❗❌:     if (at->hasProp(common_properties::atomLabel) ||
    // RDKit❗❌:         at->hasProp(common_properties::_QueryAtomGenericLabel) ||
    // RDKit❗❌:         at->hasProp(common_properties::dummyLabel) ||
    // RDKit❗❌:         at->hasProp(common_properties::_fromAttachPoint)) {
    // RDKit❗❌:       needLabels = true;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (at->hasProp(common_properties::molFileValue)) {
    // RDKit❗❌:       needValues = true;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if ((flags & SmilesWrite::CXSmilesFields::CX_COORDS) &&
    // RDKit❗❌:       mol.getNumConformers()) {
    // RDKit❗❌:     res += "(" + get_coords_block(mol, atomOrder) + ")";
    // RDKit❗❌:   }
    // RDKit❗❌:   if ((flags & SmilesWrite::CXSmilesFields::CX_ATOM_LABELS) && needLabels) {
    // RDKit❗❌:     auto lbls = get_atomlabel_block(mol, atomOrder);
    // RDKit❗❌:     if (!lbls.empty()) {
    // RDKit❗❌:       if (res.size() > 1) {
    // RDKit❗❌:         res += ",";
    // RDKit❗❌:       }
    // RDKit❗❌:       res += "$" + lbls + "$";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if ((flags & SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES) && needValues) {
    // RDKit❗❌:     if (res.size() > 1) {
    // RDKit❗❌:       res += ",";
    // RDKit❗❌:     }
    // RDKit❗❌:     res += "$_AV:" +
    // RDKit❗❌:            get_value_block(mol, atomOrder, common_properties::molFileValue) +
    // RDKit❗❌:            "$";
    // RDKit❗❌:   }
    // RDKit❗❌:   auto radblock = get_radical_block(mol, atomOrder);
    // RDKit❗❌:   if ((flags & SmilesWrite::CXSmilesFields::CX_RADICALS) && radblock.size()) {
    // RDKit❗❌:     if (res.size() > 1) {
    // RDKit❗❌:       res += ",";
    // RDKit❗❌:     }
    // RDKit❗❌:     res += radblock;
    // RDKit❗❌:     if (res.back() == ',') {
    // RDKit❗❌:       res.erase(res.size() - 1);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (flags & SmilesWrite::CXSmilesFields::CX_ATOM_PROPS) {
    // RDKit❗❌:     const auto atomblock = get_atom_props_block(mol, atomOrder);
    // RDKit❗❌:     appendToCXExtension(atomblock, res);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   const Conformer *conf = nullptr;
    // RDKit❗❌:   if (mol.getNumConformers() && (flags & SmilesWrite::CX_COORDS)) {
    // RDKit❗❌:     conf = &mol.getConformer();
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>> wedgeBonds;
    // RDKit❗❌:   if (flags & SmilesWrite::CXSmilesFields::CX_BOND_CFG) {
    // RDKit❗❌:     wedgeBonds = Chirality::pickBondsToWedge(mol, nullptr, conf);
    // RDKit❗❌:
    // RDKit❗❌:     bool includeCoords = flags & SmilesWrite::CXSmilesFields::CX_COORDS &&
    // RDKit❗❌:                          mol.getNumConformers();
    // RDKit❗❌:     const auto cfgblock = get_bond_config_block(mol, atomOrder, bondOrder,
    // RDKit❗❌:                                                 includeCoords, wedgeBonds);
    // RDKit❗❌:     appendToCXExtension(cfgblock, res);
    // RDKit❗❌:     const auto cistransblock =
    // RDKit❗❌:         get_ringbond_cistrans_block(mol, atomOrder, bondOrder);
    // RDKit❗❌:     appendToCXExtension(cistransblock, res);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // do the CX_BOND_ATROPISOMER only if CX_BOND_CFG s not done.  CX_BOND_CFG
    // RDKit❗❌:   // includes the atropisomer wedging
    // RDKit❗❌:   else if (flags & SmilesWrite::CXSmilesFields::CX_BOND_ATROPISOMER) {
    // RDKit❗❌:     Atropisomers::wedgeBondsFromAtropisomers(mol, conf, wedgeBonds);
    // RDKit❗❌:     const auto cfgblock = get_bond_config_block(
    // RDKit❗❌:         mol, atomOrder, bondOrder, conf != nullptr, wedgeBonds, true);
    // RDKit❗❌:     appendToCXExtension(cfgblock, res);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (flags & SmilesWrite::CXSmilesFields::CX_COORDINATE_BONDS) {
    // RDKit❗❌:     const auto block = get_coord_or_hydrogen_bonds_block(
    // RDKit❗❌:         mol, Bond::BondType::DATIVE, "C", atomOrder, bondOrder);
    // RDKit❗❌:     appendToCXExtension(block, res);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (flags & SmilesWrite::CXSmilesFields::CX_HYDROGEN_BONDS) {
    // RDKit❗❌:     const auto block = get_coord_or_hydrogen_bonds_block(
    // RDKit❗❌:         mol, Bond::BondType::HYDROGEN, "H", atomOrder, bondOrder);
    // RDKit❗❌:     appendToCXExtension(block, res);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (flags & SmilesWrite::CXSmilesFields::CX_ZERO_BONDS) {
    // RDKit❗❌:     const auto block = get_zerobonds_block(mol, atomOrder, bondOrder);
    // RDKit❗❌:     appendToCXExtension(block, res);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (flags & SmilesWrite::CXSmilesFields::CX_LINKNODES) {
    // RDKit❗❌:     const auto linknodeblock = get_linknodes_block(mol, atomOrder);
    // RDKit❗❌:     appendToCXExtension(linknodeblock, res);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (flags & SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO) {
    // RDKit❗❌:     const auto stereoblock =
    // RDKit❗❌:         get_enhanced_stereo_block(mol, atomOrder, wedgeBonds);
    // RDKit❗❌:     appendToCXExtension(stereoblock, res);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (flags & SmilesWrite::CXSmilesFields::CX_SGROUPS) {
    // RDKit❗❌:     const auto sgroupdatablock = get_sgroup_data_block(mol, atomOrder);
    // RDKit❗❌:     appendToCXExtension(sgroupdatablock, res);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (flags & SmilesWrite::CXSmilesFields::CX_POLYMER) {
    // RDKit❗❌:     const auto sgrouppolyblock =
    // RDKit❗❌:         get_sgroup_polymer_block(mol, atomOrder, bondOrder);
    // RDKit❗❌:     appendToCXExtension(sgrouppolyblock, res);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (flags & (SmilesWrite::CXSmilesFields::CX_SGROUPS |
    // RDKit❗❌:                SmilesWrite::CXSmilesFields::CX_POLYMER)) {
    // RDKit❗❌:     const auto sgrouphierarchyblock = get_sgroup_hierarchy_block(mol);
    // RDKit❗❌:     appendToCXExtension(sgrouphierarchyblock, res);
    // RDKit❗❌:   }
    // RDKit❗❌:   mol.clearProp("_cxsmilesOutputIndex");
    // RDKit❗❌:   if (res.size() > 1) {
    // RDKit❗❌:     res += "|";
    // RDKit❗❌:   } else {
    // RDKit❗❌:     res = "";
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    use cosmolkit_smiles::CxSmilesFields as F;
    let mut need_labels = false;
    let mut need_values = false;
    for id in atom_order {
        let atom =
            query
                .atoms()
                .get(id.index())
                .ok_or(SmartsWriteError::CxSourceAtomOutOfRange {
                    atom: *id,
                    atom_count: query.num_atoms(),
                })?;
        need_labels |= [
            "atomLabel",
            "_QueryAtomGenericLabel",
            "dummyLabel",
            "_fromAttchpt",
        ]
        .iter()
        .any(|key| atom.prop(key).is_some());
        need_values |= atom.prop("molFileValue").is_some();
    }
    let mut result = PropertyText::from("|");
    let coordinates = if fields.contains(F::COORDS)
        && (!query.conformers_2d().is_empty() || !query.conformers_3d().is_empty())
    {
        Some(write_query_cx_coordinates(query, atom_order)?)
    } else {
        None
    };
    if let Some(coordinates) = &coordinates {
        result.push_byte(b'(');
        result.extend_bytes((coordinates).as_ref());
        result.push_byte(b')');
    }
    if fields.contains(F::ATOM_LABELS) && need_labels {
        let labels = write_query_cx_atom_labels(query, atom_order)?;
        if !labels.is_empty() {
            let mut framed = PropertyText::from("$");
            framed.extend_bytes(labels.as_bytes());
            framed.push_byte(b'$');
            append_query_cx_extension(framed, &mut result);
        }
    }
    if fields.contains(F::MOLFILE_VALUES) && need_values {
        append_query_cx_extension(
            {
                let mut framed = PropertyText::from("$_AV:");
                framed.extend_bytes(
                    write_query_cx_atom_values(query, atom_order, b"molFileValue")?.as_bytes(),
                );
                framed.push_byte(b'$');
                framed
            },
            &mut result,
        );
    }
    // SOURCE evaluates radical grouping/logging even with CX_RADICALS off.
    let radicals = write_query_cx_radicals(query, atom_order, &mut |count| {
        eprintln!("unsupported number of radical electrons {count}");
    })?;
    if fields.contains(F::RADICALS) && !radicals.is_empty() {
        append_query_cx_extension(radicals, &mut result);
        if result.as_bytes().last() == Some(&b',') {
            // Native erases exactly one trailing comma after appending radblock.
            // Moving the owned vector out and back retains its allocation.
            let mut bytes = result.into_bytes();
            bytes.pop();
            result = PropertyText::from(bytes);
        }
    }
    if fields.contains(F::ATOM_PROPS) {
        append_query_cx_extension(
            write_query_cx_atom_properties(query, atom_order)?,
            &mut result,
        );
    }
    // Source selects front() here even if neither bond-generation flag is on.
    // Only the selected actual conformer is copied to permit detached graph
    // mutation; no graph/AST/atom coercion or coordinate preference is used.
    enum OwnedConformer {
        TwoD(cosmolkit_model::Conformer2D),
        ThreeD(cosmolkit_model::Conformer3D),
    }
    let conformer = if (query
        .conformers_2d()
        .len()
        .wrapping_add(query.conformers_3d().len()) as u32)
        != 0
        && fields.contains(F::COORDS)
    {
        Some(
            match query
                .first_source_conformer()
                .map_err(SmartsWriteError::CxCoordinateSource)?
                .ok_or(SmartsWriteError::CxMissingConformer)?
            {
                cosmolkit_model::CoordinateSourceConformer::TwoD(c) => {
                    OwnedConformer::TwoD(c.clone())
                }
                cosmolkit_model::CoordinateSourceConformer::ThreeD(c) => {
                    OwnedConformer::ThreeD(c.clone())
                }
            },
        )
    } else {
        None
    };
    let conf = conformer.as_ref().map(|c| match c {
        OwnedConformer::TwoD(c) => cosmolkit_core::AtropisomerConformer::TwoD(c),
        OwnedConformer::ThreeD(c) => cosmolkit_core::AtropisomerConformer::ThreeD(c),
    });
    let mut wedge_bonds = cosmolkit_core::WedgeAssignments::default();
    if fields.contains(F::BOND_CFG) || fields.contains(F::BOND_ATROPISOMER) {
        let mut rings = cosmolkit_core::RingInfo::from_source_snapshot(query.source_ring_info())
            .map_err(SmartsWriteError::CxRingInfo)?;
        let mut properties = query.source_molecule_properties();
        let generated = if fields.contains(F::BOND_CFG) {
            cosmolkit_core::pick_query_bonds_to_wedge_source(
                query,
                conf,
                &mut rings,
                Some(&mut properties),
            )
            .map(|w| {
                wedge_bonds = w;
            })
            .map_err(SmartsWriteError::CxWedge)
        } else {
            cosmolkit_core::wedge_query_bonds_from_atropisomers_source(
                query,
                &mut rings,
                Some(&mut properties),
                conf,
                &mut wedge_bonds,
            )
            .map_err(SmartsWriteError::CxStereoGroup)
        };
        // Source cache/dictionary and every actual Bond write survive a later
        // source error; this is not an atomic runtime commit boundary.
        query.replace_source_ring_info(rings.into_source_snapshot());
        query.replace_source_molecule_properties(&properties);
        generated?;
    }
    if fields.contains(F::BOND_CFG) {
        append_query_cx_extension(
            write_query_cx_bond_config(
                query,
                atom_order,
                bond_order,
                conf.is_some(),
                &wedge_bonds,
                false,
            )?,
            &mut result,
        );
        append_query_cx_extension(
            write_query_cx_ring_bond_stereo(query, atom_order, bond_order)?,
            &mut result,
        );
    } else if fields.contains(F::BOND_ATROPISOMER) {
        append_query_cx_extension(
            write_query_cx_bond_config(
                query,
                atom_order,
                bond_order,
                conf.is_some(),
                &wedge_bonds,
                true,
            )?,
            &mut result,
        );
    }
    if fields.contains(F::COORDINATE_BONDS) {
        append_query_cx_extension(
            write_query_cx_typed_bonds(query, atom_order, bond_order, BondOrder::Dative, "C")?,
            &mut result,
        );
    }
    if fields.contains(F::HYDROGEN_BONDS) {
        append_query_cx_extension(
            write_query_cx_typed_bonds(query, atom_order, bond_order, BondOrder::Hydrogen, "H")?,
            &mut result,
        );
    }
    if fields.contains(F::ZERO_BONDS) {
        append_query_cx_extension(write_query_cx_zero_bonds(query, bond_order)?, &mut result);
    }
    if fields.contains(F::LINKNODES) {
        append_query_cx_extension(write_query_cx_link_nodes(query, atom_order)?, &mut result);
    }
    if fields.contains(F::ENHANCED_STEREO) {
        append_query_cx_extension(
            write_query_cx_enhanced_stereo(query, atom_order, &wedge_bonds)?,
            &mut result,
        );
    }
    if fields.contains(F::SGROUPS) || fields.contains(F::POLYMER) {
        let mut sgroup_state = QueryCxSgroupState::from_query(query);
        let stage = (|| -> Result<(), SmartsWriteError> {
            if fields.contains(F::SGROUPS) {
                append_query_cx_extension(
                    write_query_cx_data_sgroups(query, &mut sgroup_state, atom_order)?,
                    &mut result,
                );
            }
            if fields.contains(F::POLYMER) {
                append_query_cx_extension(
                    write_query_cx_polymer_sgroups(
                        query,
                        &mut sgroup_state,
                        atom_order,
                        bond_order,
                    )?,
                    &mut result,
                );
            }
            append_query_cx_extension(
                write_query_cx_sgroup_hierarchy(&mut sgroup_state)?,
                &mut result,
            );
            Ok(())
        })();
        query.replace_source_molecule_properties(&sgroup_state.properties);
        cosmolkit_model::replace_query_substance_groups(query, sgroup_state.groups)
            .map_err(|e| SmartsWriteError::InvalidGraph(e.to_string()))?;
        stage?;
    }
    query.clear_prop("_cxsmilesOutputIndex")?;
    if result.len() == 1 {
        Ok(PropertyText::new())
    } else {
        result.push_byte(b'|');
        Ok(result)
    }
}

#[doc(hidden)]
pub fn query_atom_to_smarts(
    atom: &QueryAtom,
    params: &SmartsWriteParams,
) -> Result<PropertyText, SmartsWriteError> {
    query_atom_to_smarts_with_state(atom, params, &mut false, params.isomeric_smiles)
}

fn query_atom_to_smarts_with_state(
    atom: &QueryAtom,
    params: &SmartsWriteParams,
    stereo_written: &mut bool,
    owner_iso: bool,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string GetAtomSmarts(const Atom *atom, const SmilesWriteParams &params) {
    // RDKit❗❌:   PRECONDITION(atom, "bad atom");
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:   bool needParen = false;
    // RDKit❗❌:
    // RDKit❗❌:   // BOOST_LOG(rdInfoLog)<<"Atom: " <<qatom->getIdx()<<std::endl;
    // RDKit❗❌:   if (!atom->hasQuery()) {
    // RDKit❗❌:     res = getNonQueryAtomSmarts(atom);
    // RDKit❗❌:     // BOOST_LOG(rdInfoLog)<<"\tno query:" <<res;
    // RDKit❗❌:     return res;
    // RDKit❗❌:   }
    // RDKit❗❌:   const auto query = atom->getQuery();
    // RDKit❗❌:   PRECONDITION(query, "atom has no query");
    // RDKit❗❌:   unsigned int queryFeatures = 0;
    // RDKit❗❌:   std::string descrip = query->getDescription();
    // RDKit❗❌:   if (descrip.empty()) {
    // RDKit❗❌:     // we have simple atom - just generate the smiles and return
    // RDKit❗❌:     res = SmilesWrite::GetAtomSmiles(atom);
    // RDKit❗❌:     return res;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     if ((descrip == "AtomOr") || (descrip == "AtomAnd")) {
    // RDKit❗❌:       const QueryAtom *qatom = dynamic_cast<const QueryAtom *>(atom);
    // RDKit❗❌:       PRECONDITION(qatom, "could not convert atom to query atom");
    // RDKit❗❌:       // we have a composite query
    // RDKit❗❌:       needParen = true;
    // RDKit❗❌:       res = _recurseGetSmarts(qatom, query, query->getNegation(), queryFeatures,
    // RDKit❗❌:                               params);
    // RDKit❗❌:       if (res.length() == 1) {  // single atom symbol we don't need parens
    // RDKit❗❌:         needParen = false;
    // RDKit❗❌:       }
    // RDKit❗❌:     } else if (descrip == "RecursiveStructure") {
    // RDKit❗❌:       // it's a bare recursive structure query:
    // RDKit❗❌:       res = getRecursiveStructureQuerySmarts(query, params);
    // RDKit❗❌:       needParen = true;
    // RDKit❗❌:     } else {  // we have a simple smarts
    // RDKit❗❌:       const QueryAtom *qatom = dynamic_cast<const QueryAtom *>(atom);
    // RDKit❗❌:       PRECONDITION(qatom, "could not convert atom to query atom");
    // RDKit❗❌:       res = getAtomSmartsSimple(qatom, query, needParen, true, params);
    // RDKit❗❌:       if (query->getNegation()) {
    // RDKit❗❌:         res = "!" + res;
    // RDKit❗❌:         needParen = true;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     std::string mapNum;
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
    // RDKit❗❌:       needParen = true;
    // RDKit❗❌:       res += ":" + mapNum;
    // RDKit❗❌:     }
    // RDKit❗❌:     std::string symbol;
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::smilesSymbol, symbol)) {
    // RDKit❗❌:       needParen = true;
    // RDKit❗❌:       if (!res.empty()) {
    // RDKit❗❌:         res = symbol + ";" + res;
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res = symbol;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (needParen) {
    // RDKit❗❌:       res = "[" + res + "]";
    // RDKit❗❌:     }
    // RDKit❗❌:     return res;
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // CarrierDerived is the canonical source !hasQuery fact, not a query.
    // Return through the sole ordinary writer before reading the synthetic
    // uniform-storage predicate. Empty Native query description is not a
    // representable explicit QueryNode and is never fabricated as a fallback.
    // Not wrappers project a query's own negation bit (including cancellation).
    // Existing include_atom_maps=false is a projection option hiding that
    // property; enabled maps follow source generic string conversion exactly.
    // Cost ❌: preserve source tree order and prefix-copy sequence, but leaf
    // formatter and Vec short-output allocations retain known SSO differences.
    if atom.predicate_is_carrier_derived() {
        return non_query_atom_to_smarts(atom, owner_iso, stereo_written);
    }
    let mut features = QueryBoolFeatures::default();
    let mut needs_brackets = false;
    let (query, negated) = atom_query_without_not(atom.predicate(), false);
    let mut result = match query {
        QueryNode::And(_) | QueryNode::Or(_) => {
            let result = recurse_get_smarts_with_owner(
                atom,
                query,
                negated,
                &mut features,
                params,
                stereo_written,
                &mut query_graph_to_smarts,
                owner_iso,
            )?;
            needs_brackets = result.len() != 1;
            result
        }
        QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => {
            needs_brackets = true;
            get_recursive_structure_query_smarts(recursive, negated, params, query_graph_to_smarts)?
        }
        QueryNode::Predicate(predicate) => {
            let mut result = get_atom_smarts_simple(
                atom,
                predicate,
                &mut needs_brackets,
                true,
                owner_iso,
                stereo_written,
            )?;
            if negated {
                result.insert_byte(0, b'!');
                needs_brackets = true;
            }
            result
        }
        QueryNode::Xor(_) => return Err(SmartsWriteError::XorComposite),
        QueryNode::Not(_) => unreachable!("all source negation flags projected above"),
    };
    if params.include_atom_maps {
        let typed = atom.atom_map().map(cosmolkit_model::PropertyValue::UInt);
        if let Some(map) = typed.as_ref().or_else(|| atom.prop("molAtomMapNumber")) {
            let map = cosmolkit_core::property_value_to_string(map)?;
            needs_brackets = true;
            result.push_byte(b':');
            result.extend_bytes(map.as_bytes());
        }
    }
    if let Some(symbol) = atom.prop("smilesSymbol") {
        let mut symbol = cosmolkit_core::property_value_to_string(symbol)?;
        needs_brackets = true;
        if !result.is_empty() {
            symbol.push_byte(b';');
            symbol.extend_bytes(result.as_bytes());
        }
        result = symbol;
    }
    if needs_brackets {
        let mut framed = PropertyText::with_capacity(result.len() + 2);
        framed.push_byte(b'[');
        framed.extend_bytes(result.as_bytes());
        framed.push_byte(b']');
        result = framed;
    }
    Ok(result)
}

#[doc(hidden)]
pub fn query_bond_to_smarts(
    bond: &QueryBond,
    params: &SmartsWriteParams,
    atom_to_left_idx: Option<usize>,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string GetBondSmarts(const Bond *bond, const SmilesWriteParams &params,
    // RDKit❗❌:                           int atomToLeftIdx) {
    // RDKit❗❌:   PRECONDITION(bond, "bad bond");
    // RDKit❗❌:
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:
    // RDKit❗❌:   // BOOST_LOG(rdInfoLog) << "bond: " << bond->getIdx() << std::endl;
    // RDKit❗❌:   ;
    // RDKit❗❌:   // it is possible that we are regular single bond and we don't need to write
    // RDKit❗❌:   // anything
    // RDKit❗❌:   if (!bond->hasQuery()) {
    // RDKit❗❌:     res = getNonQueryBondSmarts(bond, atomToLeftIdx, params);
    // RDKit❗❌:     // BOOST_LOG(rdInfoLog) << "\tno query:" << res << std::endl;
    // RDKit❗❌:     return res;
    // RDKit❗❌:   }
    // RDKit❗❌:   // describeQuery(bond->getQuery());
    // RDKit❗❌:   auto qbond = dynamic_cast<const QueryBond *>(bond);
    // RDKit❗❌:   if (!qbond && ((bond->getBondType() == Bond::SINGLE) ||
    // RDKit❗❌:                  (bond->getBondType() == Bond::AROMATIC))) {
    // RDKit❗❌:     BOOST_LOG(rdInfoLog) << "\tbasic:" << res << std::endl;
    // RDKit❗❌:     return res;
    // RDKit❗❌:   }
    // RDKit❗❌:   CHECK_INVARIANT(qbond, "could not convert bond to QueryBond");
    // RDKit❗❌:
    // RDKit❗❌:   const auto query = qbond->getQuery();
    // RDKit❗❌:   CHECK_INVARIANT(query, "bond has no query");
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int queryFeatures = 0;
    // RDKit❗❌:   auto descrip = query->getDescription();
    // RDKit❗❌:   if ((descrip == "BondAnd") || (descrip == "BondOr")) {
    // RDKit❗❌:     // composite query
    // RDKit❗❌:     res = _recurseBondSmarts(bond, query, query->getNegation(), atomToLeftIdx,
    // RDKit❗❌:                              queryFeatures, params);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // simple query
    // RDKit❗❌:     if (query->getNegation()) {
    // RDKit❗❌:       res = "!";
    // RDKit❗❌:     }
    // RDKit❗❌:     res += getBondSmartsSimple(bond, query, atomToLeftIdx, params);
    // RDKit❗❌:   }
    // RDKit❗❌:   // BOOST_LOG(rdInfoLog) << "\t  query:" << descrip << " " << res << std::endl;
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // Canonical origin is the source hasQuery test; a carrier-derived helper
    // predicate must never be rendered. The explicit QueryBond projection
    // guarantees the source dynamic QueryBond class/non-null query invariant;
    // arbitrary third-party Bond subclasses are not represented by this type.
    // Preserve source own flag, composite binary recursion, and simple prefix
    // before appending the leaf text (including an empty leaf's literal !).
    // Cost ❌: ordinary kernel is O(1); explicit tree rendering retains source
    // cumulative prefix copies and known leaf/Vec short-string allocations
    // versus Native SSO. Borrow all inputs, with no query/carrier cloning.
    if bond.predicate_is_carrier_derived() {
        return Ok(non_query_bond_to_smarts(bond.bond(), atom_to_left_idx, params)?.into());
    }
    let (query, negated) = bond_query_without_not(bond.predicate(), false);
    let mut features = QueryBoolFeatures::default();
    match query {
        QueryNode::And(_) | QueryNode::Or(_) => recurse_bond_smarts(
            bond.bond(),
            query,
            negated,
            atom_to_left_idx,
            &mut features,
            params,
        ),
        QueryNode::Predicate(predicate) => {
            let leaf = get_bond_smarts_simple(bond.bond(), predicate, atom_to_left_idx, params)?;
            let mut result = PropertyText::new();
            if negated {
                result.push_byte(b'!');
            }
            result.extend_bytes(leaf.as_bytes());
            Ok(result)
        }
        QueryNode::Xor(_) => Err(SmartsWriteError::XorComposite),
        QueryNode::Not(_) => unreachable!("source negation wrappers projected above"),
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
struct SmartsWriteResult {
    smarts: PropertyText,
    atom_ordering: Vec<AtomId>,
    bond_ordering: Vec<BondId>,
    source_orders_written: bool,
    // Actual computed-dictionary effects of typed source order setters.
    // No IntVector/String surrogate for the independent UInt-vector tag.
    source_properties: Option<cosmolkit_model::MoleculeProperties>,
}

fn combine_child_smarts(
    child1: PropertyText,
    features1: QueryBoolFeatures,
    child2: PropertyText,
    features2: QueryBoolFeatures,
    description: &str,
    features: &mut QueryBoolFeatures,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❗: std::string _combineChildSmarts(std::string cs1, unsigned int features1,
    // RDKit❗❗:                                 std::string cs2, unsigned int features2,
    // RDKit❗❗:                                 std::string descrip, unsigned int &features) {
    // RDKit❗❗:   std::string res = "";
    // RDKit❗❗:   if ((descrip.find("Or") > 0) && (descrip.find("Or") < descrip.length())) {
    // RDKit❗❗:     // if either of child smarts already have a "," and ";" we can't have one
    // RDKit❗❗:     // more OR here
    // RDKit❗❗:     if ((features1 & static_cast<unsigned int>(QueryBoolFeatures::HAS_LOWAND) &&
    // RDKit❗❗:          features1 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR)) ||
    // RDKit❗❗:         (features2 & static_cast<unsigned int>(QueryBoolFeatures::HAS_LOWAND) &&
    // RDKit❗❗:          features2 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR))) {
    // RDKit❗❗:       throw ValueErrorException(
    // RDKit❗❗:           "This is a non-smartable query - OR above and below AND in the "
    // RDKit❗❗:           "binary tree");
    // RDKit❗❗:     }
    // RDKit❗❗:     res += cs1;
    // RDKit❗❗:     if (!(cs1.empty() || cs2.empty())) {
    // RDKit❗❗:       res += ",";
    // RDKit❗❗:     }
    // RDKit❗❗:     res += cs2;
    // RDKit❗❗:     features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_OR);
    // RDKit❗❗:   } else if ((descrip.find("And") > 0) &&
    // RDKit❗❗:              (descrip.find("And") < descrip.length())) {
    // RDKit❗❗:     std::string symb;
    // RDKit❗❗:     if (features1 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR) ||
    // RDKit❗❗:         features2 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR)) {
    // RDKit❗❗:       symb = ";";
    // RDKit❗❗:       features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_LOWAND);
    // RDKit❗❗:     } else {
    // RDKit❗❗:       symb = "&";
    // RDKit❗❗:       features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_AND);
    // RDKit❗❗:     }
    // RDKit❗❗:     res += cs1;
    // RDKit❗❗:     if (!(cs1.empty() || cs2.empty())) {
    // RDKit❗❗:       res += symb;
    // RDKit❗❗:     }
    // RDKit❗❗:     res += cs2;
    // RDKit❗❗:   } else {
    // RDKit❗❗:     std::stringstream err;
    // RDKit❗❗:     err << "Don't know how to combine using " << descrip;
    // RDKit❗❗:     throw ValueErrorException(err.str());
    // RDKit❗❗:   }
    // RDKit❗❗:   features |= features1;
    // RDKit❗❗:   features |= features2;
    // RDKit❗❗:
    // RDKit❗❗:   return res;
    // RDKit❗❗: }
    // Source first-occurrence searches require a positive byte offset. OR
    // takes precedence even if And occurs earlier. Error paths leave features
    // unchanged; only successful combinations OR in both child words.
    // Cost review: O(description bytes + child bytes), no tree scan or clone.
    // A single Vec allocation/copy can save long-string append growth, but
    // source SSO avoids heap allocation for short/empty results. Mixed regimes
    // leave aggregate cost unresolved; child buffers are copied, not merged.
    let has_or = description.find("Or").is_some_and(|position| position > 0);
    let has_and = description.find("And").is_some_and(|position| position > 0);
    let separator = if has_or {
        if (features1.contains(QueryBoolFeatures::HAS_LOW_AND)
            && features1.contains(QueryBoolFeatures::HAS_OR))
            || (features2.contains(QueryBoolFeatures::HAS_LOW_AND)
                && features2.contains(QueryBoolFeatures::HAS_OR))
        {
            return Err(SmartsWriteError::OrAboveAndBelowAnd);
        }
        features.insert(QueryBoolFeatures::HAS_OR);
        ","
    } else if has_and {
        if features1.contains(QueryBoolFeatures::HAS_OR)
            || features2.contains(QueryBoolFeatures::HAS_OR)
        {
            features.insert(QueryBoolFeatures::HAS_LOW_AND);
            ";"
        } else {
            features.insert(QueryBoolFeatures::HAS_AND);
            "&"
        }
    } else {
        return Err(SmartsWriteError::UnknownCombination {
            description: description.to_owned(),
        });
    };

    let mut result = PropertyText::with_capacity(child1.len() + child2.len() + 1);
    result.extend_bytes((&child1).as_ref());
    if !child1.is_empty() && !child2.is_empty() {
        result.extend_bytes((separator).as_ref());
    }
    result.extend_bytes((&child2).as_ref());
    *features |= features1;
    *features |= features2;
    Ok(result)
}

fn describe_query<T>(query: &QueryNode<T>, leader: String) {
    // RDKit❗✔️: void describeQuery(const T *query, std::string leader = "\t") {
    // RDKit❗✔️:   // BOOST_LOG(rdInfoLog) << leader << query->getDescription() << std::endl;
    // RDKit❗✔️:   typename T::CHILD_VECT_CI iter;
    // RDKit❗✔️:   for (iter = query->beginChildren(); iter != query->endChildren(); ++iter) {
    // RDKit❗✔️:     describeQuery(iter->get(), leader + "\t");
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️:   void setNegation(bool what) { this->df_negate = what; }
    // RDKit❗✔️:   bool getNegation() const { return this->df_negate; }
    // RDKit❗✔️:   CHILD_VECT_CI beginChildren() const { return this->d_children.begin(); }
    // RDKit❗✔️:   CHILD_VECT_CI endChildren() const { return this->d_children.end(); }
    // Native logging is commented out: do not inspect, stringify, evaluate,
    // clone or mutate a predicate. Visit actual child queries in their order.
    // Model Not is a negation wrapper around the same source query object,
    // not an extra source child edge: no extra leader tab/allocation for Not.
    // Every actual child edge allocates one extended leader as Native does;
    // no child sorting, recursive-molecule traversal or secondary AST.
    let mut node = query;
    while let QueryNode::Not(child) = node {
        node = child;
    }
    match node {
        QueryNode::Predicate(_) => {}
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
            for child in children {
                describe_query(child, format!("{leader}\t"));
            }
        }
        QueryNode::Not(_) => unreachable!("negation wrappers consumed above"),
    }
}

fn range_prefix(data_function: AtomRangeDataFunction) -> &'static str {
    match data_function {
        AtomRangeDataFunction::ExplicitDegree => "D",
        AtomRangeDataFunction::NonHydrogenDegree => "d",
        AtomRangeDataFunction::TotalDegree => "X",
        AtomRangeDataFunction::TotalValence => "v",
        AtomRangeDataFunction::NumAtomRings => "R",
        AtomRangeDataFunction::NumHeteroatomNeighbors => "z",
        AtomRangeDataFunction::NumAliphaticHeteroatomNeighbors => "Z",
        AtomRangeDataFunction::MinRingSize => "r",
        AtomRangeDataFunction::RingBondCount => "x",
        AtomRangeDataFunction::ImplicitHydrogenCount => "h",
        AtomRangeDataFunction::FormalCharge => "+",
        AtomRangeDataFunction::NegativeFormalCharge => "-",
        AtomRangeDataFunction::AtomRingSize { .. } => "k",
    }
}

fn get_atom_smarts_simple(
    atom: &QueryAtom,
    query: &AtomQueryPredicate,
    need_paren: &mut bool,
    check_for_symbol: bool,
    do_isomeric_smarts: bool,
    stereo_written: &mut bool,
) -> Result<PropertyText, SmartsWriteError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/SmilesParse/SmartsWrite.cpp :: getAtomSmartsSimple
    // RDKit❗❌: std::string getAtomSmartsSimple(const QueryAtom *qatom,
    // RDKit❗❌:                                 const Atom::QUERYATOM_QUERY *query,
    // RDKit❗❌:                                 bool &needParen, bool checkForSymbol,
    // RDKit❗❌:                                 const SmilesWriteParams &) {
    // RDKit❗❌:   PRECONDITION(query, "bad query");
    // RDKit❗❌:
    // RDKit❗❌:   auto *equery = dynamic_cast<const ATOM_EQUALS_QUERY *>(query);
    // RDKit❗❌:
    // RDKit❗❌:   std::string descrip = query->getDescription();
    // RDKit❗❌:   bool hasVal = false;
    // RDKit❗❌:   enum class Modifiers : std::uint8_t {
    // RDKit❗❌:     NONE,
    // RDKit❗❌:     RANGE,
    // RDKit❗❌:     LESS,
    // RDKit❗❌:     GREATER
    // RDKit❗❌:   };
    // RDKit❗❌:   Modifiers mods = Modifiers::NONE;
    // RDKit❗❌:   if (boost::starts_with(descrip, "range_")) {
    // RDKit❗❌:     mods = Modifiers::RANGE;
    // RDKit❗❌:     descrip = descrip.substr(6);
    // RDKit❗❌:   } else if (boost::starts_with(descrip, "less_")) {
    // RDKit❗❌:     mods = Modifiers::LESS;
    // RDKit❗❌:     descrip = descrip.substr(5);
    // RDKit❗❌:   } else if (boost::starts_with(descrip, "greater_")) {
    // RDKit❗❌:     mods = Modifiers::GREATER;
    // RDKit❗❌:     descrip = descrip.substr(8);
    // RDKit❗❌:   }
    // RDKit❗❌:   std::stringstream res;
    // RDKit❗❌:   if (descrip == "AtomImplicitHCount") {
    // RDKit❗❌:     res << "h";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomHasImplicitH") {
    // RDKit❗❌:     res << "h";
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomTotalValence") {
    // RDKit❗❌:     res << "v";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomAtomicNum") {
    // RDKit❗❌:     if (!qatom->hasProp(common_properties::smilesSymbol)) {
    // RDKit❗❌:       res << "#";
    // RDKit❗❌:       hasVal = true;
    // RDKit❗❌:       needParen = true;
    // RDKit❗❌:     }
    // RDKit❗❌:   } else if (descrip == "AtomExplicitDegree") {
    // RDKit❗❌:     res << "D";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomNonHydrogenDegree") {
    // RDKit❗❌:     res << "d";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomTotalDegree") {
    // RDKit❗❌:     res << "X";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomHasRingBond") {
    // RDKit❗❌:     res << "x";
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomHCount") {
    // RDKit❗❌:     res << "H";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomIsAliphatic") {
    // RDKit❗❌:     res << "A";
    // RDKit❗❌:     needParen = false;
    // RDKit❗❌:   } else if (descrip == "AtomIsAromatic") {
    // RDKit❗❌:     res << "a";
    // RDKit❗❌:     needParen = false;
    // RDKit❗❌:   } else if (descrip == "AtomNull") {
    // RDKit❗❌:     res << "*";
    // RDKit❗❌:     needParen = false;
    // RDKit❗❌:   } else if (descrip == "AtomInRing") {
    // RDKit❗❌:     res << "R";
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomMinRingSize") {
    // RDKit❗❌:     res << "r";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomRingSize") {
    // RDKit❗❌:     res << "k";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomInNRings") {
    // RDKit❗❌:     res << "R";
    // RDKit❗❌:     if (mods == Modifiers::NONE && equery && equery->getVal() >= 0) {
    // RDKit❗❌:       hasVal = true;
    // RDKit❗❌:     }
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomHasHeteroatomNeighbors") {
    // RDKit❗❌:     res << "z";
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomNumHeteroatomNeighbors") {
    // RDKit❗❌:     res << "z";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomHasAliphaticHeteroatomNeighbors") {
    // RDKit❗❌:     res << "Z";
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomNumAliphaticHeteroatomNeighbors") {
    // RDKit❗❌:     res << "Z";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomFormalCharge") {
    // RDKit❗❌:     int val = equery ? equery->getVal() : 0;
    // RDKit❗❌:
    // RDKit❗❌:     if (val < 0) {
    // RDKit❗❌:       res << "-";
    // RDKit❗❌:     } else {
    // RDKit❗❌:       res << "+";
    // RDKit❗❌:     }
    // RDKit❗❌:     if (mods == Modifiers::NONE && abs(val) != 1) {
    // RDKit❗❌:       res << abs(val);
    // RDKit❗❌:     }
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomNegativeFormalCharge") {
    // RDKit❗❌:     int val = equery ? equery->getVal() : 0;
    // RDKit❗❌:     if (val < 0) {
    // RDKit❗❌:       res << "+";
    // RDKit❗❌:     } else {
    // RDKit❗❌:       res << "-";
    // RDKit❗❌:     }
    // RDKit❗❌:     if (mods == Modifiers::NONE && abs(val) != 1) {
    // RDKit❗❌:       res << abs(val);
    // RDKit❗❌:     }
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomHybridization" && equery) {
    // RDKit❗❌:     res << "^";
    // RDKit❗❌:     switch (equery->getVal()) {
    // RDKit❗❌:       case Atom::S:
    // RDKit❗❌:         res << "0";
    // RDKit❗❌:         break;
    // RDKit❗❌:       case Atom::SP:
    // RDKit❗❌:         res << "1";
    // RDKit❗❌:         break;
    // RDKit❗❌:       case Atom::SP2:
    // RDKit❗❌:         res << "2";
    // RDKit❗❌:         break;
    // RDKit❗❌:       case Atom::SP3:
    // RDKit❗❌:         res << "3";
    // RDKit❗❌:         break;
    // RDKit❗❌:       case Atom::SP3D:
    // RDKit❗❌:         res << "4";
    // RDKit❗❌:         break;
    // RDKit❗❌:       case Atom::SP3D2:
    // RDKit❗❌:         res << "5";
    // RDKit❗❌:         break;
    // RDKit❗❌:     }
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomMass" && equery) {
    // RDKit❗❌:     res << equery->getVal() / massIntegerConversionFactor << "*";
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomIsotope" && equery) {
    // RDKit❗❌:     res << equery->getVal() << "*";
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomRingBondCount") {
    // RDKit❗❌:     res << "x";
    // RDKit❗❌:     hasVal = true;
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomUnsaturated") {
    // RDKit❗❌:     res << "$(*=,:,#*)";
    // RDKit❗❌:     needParen = true;
    // RDKit❗❌:   } else if (descrip == "AtomType" && equery) {
    // RDKit❗❌:     int atNum;
    // RDKit❗❌:     bool isAromatic;
    // RDKit❗❌:     parseAtomType(equery->getVal(), atNum, isAromatic);
    // RDKit❗❌:     if (!checkForSymbol || !qatom->hasProp(common_properties::smilesSymbol)) {
    // RDKit❗❌:       std::string symbol = PeriodicTable::getTable()->getElementSymbol(atNum);
    // RDKit❗❌:       if (isAromatic) {
    // RDKit❗❌:         symbol[0] += ('a' - 'A');
    // RDKit❗❌:       }
    // RDKit❗❌:       res << symbol;
    // RDKit❗❌:
    // RDKit❗❌:       if (!SmilesWrite::inOrganicSubset(atNum)) {
    // RDKit❗❌:         needParen = true;
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       if (isAromatic) {
    // RDKit❗❌:         res << "a";
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res << "A";
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "Cannot write SMARTS for query type : " << descrip
    // RDKit❗❌:         << ". Ignoring it." << std::endl;
    // RDKit❗❌:     res << "*";
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (mods != Modifiers::NONE) {
    // RDKit❗❌:     res << "{";
    // RDKit❗❌:     const ATOM_RANGE_QUERY *rquery = nullptr;
    // RDKit❗❌:     switch (mods) {
    // RDKit❗❌:       case Modifiers::LESS:
    // RDKit❗❌:         res << equery->getVal() << "-";
    // RDKit❗❌:         break;
    // RDKit❗❌:       case Modifiers::RANGE:
    // RDKit❗❌:         rquery = dynamic_cast<const ATOM_RANGE_QUERY *>(query);
    // RDKit❗❌:         CHECK_INVARIANT(rquery, "query could not be converted to range query");
    // RDKit❗❌:         res << ((const ATOM_RANGE_QUERY *)query)->getLower() << "-"
    // RDKit❗❌:             << ((const ATOM_RANGE_QUERY *)query)->getUpper();
    // RDKit❗❌:         break;
    // RDKit❗❌:       case Modifiers::GREATER:
    // RDKit❗❌:         res << "-" << equery->getVal();
    // RDKit❗❌:         break;
    // RDKit❗❌:       default:
    // RDKit❗❌:         break;
    // RDKit❗❌:     }
    // RDKit❗❌:     res << "}";
    // RDKit❗❌:   } else if (hasVal) {
    // RDKit❗❌:     res << equery->getVal();
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // handle atomic stereochemistry
    // RDKit❗❌:   if (qatom->hasOwningMol() &&
    // RDKit❗❌:       qatom->getOwningMol().hasProp(common_properties::_doIsoSmiles)) {
    // RDKit❗❌:     if (qatom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗❌:         !qatom->hasProp(_qatomHasStereoSet) &&
    // RDKit❗❌:         !qatom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❗❌:       qatom->setProp(_qatomHasStereoSet, 1);
    // RDKit❗❌:       switch (qatom->getChiralTag()) {
    // RDKit❗❌:         case Atom::CHI_TETRAHEDRAL_CW:
    // RDKit❗❌:           res << "@@";
    // RDKit❗❌:           needParen = true;
    // RDKit❗❌:           break;
    // RDKit❗❌:         case Atom::CHI_TETRAHEDRAL_CCW:
    // RDKit❗❌:           res << "@";
    // RDKit❗❌:           needParen = true;
    // RDKit❗❌:           break;
    // RDKit❗❌:         default:
    // RDKit❗❌:           break;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return res.str();
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // needParen is an in/out parameter. Atomic-number symbol suppression,
    // organic AtomType and source warning fallback do not overwrite it.
    // Detached caller booleans project owner _doIsoSmiles presence and the
    // temporary stereo-written property; raw atom marker presence still wins.
    // Cost ❌: typed formatting uses temporary strings versus a single source
    // stream. No graph traversal, matching heuristic or reordering is added.
    match query {
        AtomQueryPredicate::Any
        | AtomQueryPredicate::IsAromatic(_)
        | AtomQueryPredicate::AtomicNumber(_)
        | AtomQueryPredicate::AtomType { .. } => {}
        AtomQueryPredicate::RecursiveSmarts(_)
        | AtomQueryPredicate::AtomicNumberIn(_)
        | AtomQueryPredicate::AtomicNumberNotIn(_)
        | AtomQueryPredicate::NumRadicalElectrons(_)
        | AtomQueryPredicate::HasChiralTag
        | AtomQueryPredicate::MissingChiralTag
        | AtomQueryPredicate::ImplicitValence(_)
        | AtomQueryPredicate::ExplicitValence(_)
        | AtomQueryPredicate::HeavyAtomDegree(_)
        | AtomQueryPredicate::IsBridgehead
        | AtomQueryPredicate::HasProperty(_)
        | AtomQueryPredicate::PropertyValue { .. }
        | AtomQueryPredicate::RGroupLabel(_)
        | AtomQueryPredicate::MolFileAlias(_)
        | AtomQueryPredicate::ChiralTagMatch(_)
        | AtomQueryPredicate::ChiralPermutationMatch(_)
        | AtomQueryPredicate::UnsupportedFeature(_) => {}
        _ => *need_paren = true,
    }
    let result = match query {
        AtomQueryPredicate::Any => {
            *need_paren = false;
            "*".to_owned()
        }
        AtomQueryPredicate::AtomicNumber(value) => {
            if atom.prop("smilesSymbol").is_some() {
                String::new()
            } else {
                *need_paren = true;
                format!("#{value}")
            }
        }
        AtomQueryPredicate::AtomType {
            atomic_number,
            aromatic,
        } => {
            let (atomic_number, aromatic) = crate::query_behavior::parse_atom_type(
                crate::query_behavior::make_atom_type(i32::from(*atomic_number), *aromatic),
            );
            if check_for_symbol && atom.prop("smilesSymbol").is_some() {
                if aromatic { "a".into() } else { "A".into() }
            } else {
                let element = u8::try_from(atomic_number)
                    .ok()
                    .and_then(Element::from_atomic_number)
                    .ok_or(SmartsWriteError::AtomTypeAtomicNumber { atomic_number })?;
                let mut symbol = element.symbol().as_bytes().to_vec();
                if aromatic {
                    symbol[0] = symbol[0].wrapping_add(b'a' - b'A');
                }
                if !cosmolkit_core::is_rdkit_organic_subset(element.atomic_number()) {
                    *need_paren = true;
                }
                String::from_utf8(symbol).expect(
                    "canonical element symbols and source first-byte case arithmetic are ASCII",
                )
            }
        }
        AtomQueryPredicate::ImplicitHydrogenCount(value) => format!("h{value}"),
        AtomQueryPredicate::HasImplicitHydrogen => "h".into(),
        AtomQueryPredicate::TotalValence(value) => format!("v{value}"),
        AtomQueryPredicate::ExplicitDegree(value) => format!("D{value}"),
        AtomQueryPredicate::NonHydrogenDegree(value) => format!("d{}", *value as i32),
        AtomQueryPredicate::TotalDegree(value) => format!("X{value}"),
        AtomQueryPredicate::HasRingBond => "x".into(),
        AtomQueryPredicate::HydrogenCount(value) => format!("H{value}"),
        AtomQueryPredicate::IsAromatic(false) => {
            *need_paren = false;
            "A".into()
        }
        AtomQueryPredicate::IsAromatic(true) => {
            *need_paren = false;
            "a".into()
        }
        AtomQueryPredicate::InRing => "R".into(),
        AtomQueryPredicate::SmallestRingSize(value) => format!("r{value}"),
        AtomQueryPredicate::InRingOfSize(value) => format!("k{value}"),
        AtomQueryPredicate::NumAtomRings(value) if *value >= 0 => format!("R{value}"),
        AtomQueryPredicate::NumAtomRings(_) => "R".into(),
        AtomQueryPredicate::HasHeteroatomNeighbors => "z".into(),
        AtomQueryPredicate::NumHeteroatomNeighbors(value) => format!("z{value}"),
        AtomQueryPredicate::HasAliphaticHeteroatomNeighbors => "Z".into(),
        AtomQueryPredicate::NumAliphaticHeteroatomNeighbors(value) => format!("Z{value}"),
        AtomQueryPredicate::FormalCharge(value) => {
            if *value == i32::MIN {
                return Err(SmartsWriteError::ChargeMagnitudeOverflow { value: *value });
            }
            let sign = if *value < 0 { '-' } else { '+' };
            if value.unsigned_abs() == 1 {
                sign.to_string()
            } else {
                format!("{sign}{}", value.unsigned_abs())
            }
        }
        AtomQueryPredicate::NegativeFormalCharge(value) => {
            if *value == i32::MIN {
                return Err(SmartsWriteError::ChargeMagnitudeOverflow { value: *value });
            }
            let sign = if *value < 0 { '+' } else { '-' };
            if value.unsigned_abs() == 1 {
                sign.to_string()
            } else {
                format!("{sign}{}", value.unsigned_abs())
            }
        }
        AtomQueryPredicate::HybridizationMatch(value) => format!(
            "^{}",
            match value {
                Hybridization::S => "0",
                Hybridization::Sp => "1",
                Hybridization::Sp2 => "2",
                Hybridization::Sp3 => "3",
                Hybridization::Sp3d => "4",
                Hybridization::Sp3d2 => "5",
                Hybridization::Unspecified | Hybridization::Sp2d | Hybridization::Other => "",
            }
        ),
        AtomQueryPredicate::Mass(value) => format!("{value}*"),
        AtomQueryPredicate::Isotope(value) => format!("{value}*"),
        AtomQueryPredicate::RingBondCount(value) => format!("x{value}"),
        AtomQueryPredicate::IsUnsaturated => "$(*=,:,#*)".into(),
        AtomQueryPredicate::Range(range) => {
            let (bounds, data_function) = range.writer_parts();
            // Less/Greater derive from EqualityQuery; RangeQuery does not.
            // Thus charge prefix uses the actual comparison threshold for the
            // former and the source zero-valued null-equery fallback for range.
            let equery_val = match bounds {
                AtomRangeBounds::LessEqual(v) | AtomRangeBounds::GreaterEqual(v) => v,
                AtomRangeBounds::Inclusive { .. } => 0,
            };
            let prefix = match data_function {
                AtomRangeDataFunction::FormalCharge => {
                    if equery_val < 0 {
                        "-"
                    } else {
                        "+"
                    }
                }
                AtomRangeDataFunction::NegativeFormalCharge => {
                    if equery_val < 0 {
                        "+"
                    } else {
                        "-"
                    }
                }
                _ => range_prefix(data_function),
            };
            let bounds = match bounds {
                AtomRangeBounds::LessEqual(value) => format!("{value}-"),
                AtomRangeBounds::GreaterEqual(value) => format!("-{value}"),
                AtomRangeBounds::Inclusive { lower, upper, .. } => format!("{lower}-{upper}"),
            };
            format!("{prefix}{{{bounds}}}")
        }
        AtomQueryPredicate::ExplicitDegreeLessEqual(value) => format!("D{{{value}-}}"),
        AtomQueryPredicate::NonHydrogenDegreeLessEqual(value) => format!("d{{{}-}}", *value as i32),
        AtomQueryPredicate::NonHydrogenDegreeGreaterEqual(value) => {
            format!("d{{-{}}}", *value as i32)
        }
        AtomQueryPredicate::TotalDegreeLessEqual(value) => format!("X{{{value}-}}"),
        AtomQueryPredicate::TotalDegreeGreaterEqual(value) => format!("X{{-{value}}}"),
        AtomQueryPredicate::TotalValenceLessEqual(value) => format!("v{{{value}-}}"),
        AtomQueryPredicate::TotalValenceGreaterEqual(value) => format!("v{{-{value}}}"),
        AtomQueryPredicate::RingBondCountLessEqual(value) => format!("x{{{value}-}}"),
        AtomQueryPredicate::ImplicitHydrogenCountLessEqual(value) => format!("h{{{value}-}}"),
        AtomQueryPredicate::InRingOfSizeLessEqual(value) => format!("k{{{value}-}}"),
        AtomQueryPredicate::InRingOfSizeGreaterEqual(value) => format!("k{{-{value}}}"),
        AtomQueryPredicate::SmallestRingSizeLessEqual(value) => format!("r{{{value}-}}"),
        AtomQueryPredicate::SmallestRingSizeGreaterEqual(value) => format!("r{{-{value}}}"),
        AtomQueryPredicate::DegreeLessEqual(value) => format!("D{{{value}-}}"),
        AtomQueryPredicate::DegreeGreaterEqual(value) => format!("D{{-{value}}}"),
        AtomQueryPredicate::RecursiveSmarts(_)
        | AtomQueryPredicate::NumRadicalElectrons(_)
        | AtomQueryPredicate::HasChiralTag
        | AtomQueryPredicate::MissingChiralTag
        | AtomQueryPredicate::ImplicitValence(_)
        | AtomQueryPredicate::ExplicitValence(_)
        | AtomQueryPredicate::HeavyAtomDegree(_)
        | AtomQueryPredicate::IsBridgehead
        | AtomQueryPredicate::HasProperty(_)
        | AtomQueryPredicate::PropertyValue { .. } => {
            let description = match query {
                AtomQueryPredicate::RecursiveSmarts(_) => "RecursiveStructure",
                AtomQueryPredicate::NumRadicalElectrons(_) => "AtomNumRadicalElectrons",
                AtomQueryPredicate::HasChiralTag => "AtomHasChiralTag",
                AtomQueryPredicate::MissingChiralTag => "AtomMissingChiralTag",
                AtomQueryPredicate::ImplicitValence(_) => "AtomImplicitValence",
                AtomQueryPredicate::ExplicitValence(_) => "AtomExplicitValence",
                AtomQueryPredicate::HeavyAtomDegree(_) => "AtomHeavyAtomDegree",
                AtomQueryPredicate::IsBridgehead => "AtomIsBridgehead",
                AtomQueryPredicate::HasProperty(_) => "HasProp",
                AtomQueryPredicate::PropertyValue { .. } => "HasPropWithValue",
                _ => unreachable!(),
            };
            // This is the literal source-defined fallback, not an inferred
            // replacement for an unmodeled predicate's matching behavior.
            eprintln!("Cannot write SMARTS for query type : {description}. Ignoring it.");
            "*".into()
        }
        AtomQueryPredicate::AtomicNumberIn(_)
        | AtomQueryPredicate::AtomicNumberNotIn(_)
        | AtomQueryPredicate::RGroupLabel(_)
        | AtomQueryPredicate::MolFileAlias(_)
        | AtomQueryPredicate::ChiralTagMatch(_)
        | AtomQueryPredicate::ChiralPermutationMatch(_)
        | AtomQueryPredicate::UnsupportedFeature(_) => {
            // These independent projections do not carry a Native leaf query
            // class/description. Do not invent a description or value from it.
            return Err(SmartsWriteError::UnsupportedAtomQuery {
                predicate: query.clone(),
            });
        }
    };
    let mut result = PropertyText::from(result);

    if do_isomeric_smarts
        && atom.chiral_tag() != ChiralTag::Unspecified
        && !*stereo_written
        && atom.prop("_qatomHasStereoSet").is_none()
        && atom.prop("_brokenChirality").is_none()
    {
        // Native sets the temporary property before switching, including
        // non-tetrahedral tags which emit no stereo token.
        *stereo_written = true;
        match atom.chiral_tag() {
            ChiralTag::TetrahedralCw => {
                result.extend_bytes(("@@").as_ref());
                *need_paren = true;
                *stereo_written = true;
            }
            ChiralTag::TetrahedralCcw => {
                result.push_byte(b'@');
                *need_paren = true;
                *stereo_written = true;
            }
            _ => {}
        }
    }
    Ok(result)
}

fn get_recursive_structure_query_smarts<F>(
    query: &RecursiveStructureQuery,
    negated: bool,
    params: &SmartsWriteParams,
    write_molecule: F,
) -> Result<PropertyText, SmartsWriteError>
where
    F: FnOnce(&crate::QueryGraph, &SmartsWriteParams) -> Result<PropertyText, SmartsWriteError>,
{
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/SmilesParse/SmartsWrite.cpp :: getRecursiveStructureQuerySmarts
    // RDKit❗❗: std::string getRecursiveStructureQuerySmarts(
    // RDKit❗❗:     const QueryAtom::QUERYATOM_QUERY *query, const SmilesWriteParams &params) {
    // RDKit❗❗:   PRECONDITION(query, "bad query");
    // RDKit❗❗:   PRECONDITION(query->getDescription() == "RecursiveStructure", "bad query");
    // RDKit❗❗:   const auto *rquery = dynamic_cast<const RecursiveStructureQuery *>(query);
    // RDKit❗❗:   PRECONDITION(rquery, "could not convert query to RecursiveStructureQuery");
    // RDKit❗❗:   auto *qmol = const_cast<ROMol *>(rquery->getQueryMol());
    // RDKit❗❗:   std::string res = MolToSmarts(*qmol, params);
    // RDKit❗❗:   res = "$(" + res + ")";
    // RDKit❗❗:   if (rquery->getNegation()) {
    // RDKit❗❗:     res = "!" + res;
    // RDKit❗❗:   }
    // RDKit❗❗:   return res;
    // RDKit❗❗: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Reuse the complete canonical writer once, before allocating output.
    // Typed RecursiveStructureQuery enforces source class/description identity;
    // absent inner graph is a structural error at the native null dereference,
    // never a cached-text fallback. Caller negation is the enclosing query
    // node's source flag. The immutable detached projection does not expose
    // Native owner mutation/const_cast scratch state.
    // Cost ❗: both implementations copy O(n) bytes with one writer call.
    // libstdc++ strings keep up to 15 chars in local storage; this Vec-backed
    // wrapper always allocates. For longer queries, one exact-capacity buffer
    // can reduce source concatenation allocations/copies. Neither regime
    // establishes an aggregate win. No query clone or graph/set/serial scan.
    let query_molecule = query
        .query_graph()
        .ok_or(SmartsWriteError::MissingRecursiveQueryMolecule)?;
    let inner = write_molecule(query_molecule, params)?;
    let mut result = PropertyText::with_capacity(inner.len() + if negated { 4 } else { 3 });
    if negated {
        result.push_byte(b'!');
    }
    result.extend_bytes(("$(").as_ref());
    result.extend_bytes((&inner).as_ref());
    result.push_byte(b')');
    Ok(result)
}

fn get_basic_bond_repr(
    bond_order: BondOrder,
    direction: BondDirection,
    reverse_dative: bool,
    params: &SmartsWriteParams,
) -> &'static str {
    // RDKit❗✔️: std::string getBasicBondRepr(Bond::BondType typ, Bond::BondDir dir,
    // RDKit❗✔️:                              bool reverseDative,
    // RDKit❗✔️:                              const SmilesWriteParams &params) {
    // RDKit❗✔️:   std::string res;
    // RDKit❗✔️:   switch (typ) {
    // RDKit❗✔️:     case Bond::SINGLE:
    // RDKit❗✔️:       res = "-";
    // RDKit❗✔️:       if (params.doIsomericSmiles) {
    // RDKit❗✔️:         if (dir == Bond::ENDDOWNRIGHT) {
    // RDKit❗✔️:           res = "\\";
    // RDKit❗✔️:         } else if (dir == Bond::ENDUPRIGHT) {
    // RDKit❗✔️:           res = "/";
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::DOUBLE:
    // RDKit❗✔️:       res = "=";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::TRIPLE:
    // RDKit❗✔️:       res = "#";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::QUADRUPLE:
    // RDKit❗✔️:       res = "$";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::AROMATIC:
    // RDKit❗✔️:       res = ":";
    // RDKit❗✔️:       if (params.doIsomericSmiles) {
    // RDKit❗✔️:         if (dir == Bond::ENDDOWNRIGHT) {
    // RDKit❗✔️:           res = "\\";
    // RDKit❗✔️:         } else if (dir == Bond::ENDUPRIGHT) {
    // RDKit❗✔️:           res = "/";
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::DATIVE:
    // RDKit❗✔️:       if (params.includeDativeBonds) {
    // RDKit❗✔️:         if (reverseDative) {
    // RDKit❗✔️:           res = "<-";
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           res = "->";
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         res = "-";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::ZERO:
    // RDKit❗✔️:       res = "~";  // Actually means "any", but we use ~ for unknown bond types
    // RDKit❗✔️:                   // in SMILES,
    // RDKit❗✔️:       break;      // and this will match a ZOB.
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       res = "";
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // Reuse the source-exact constant switch and directional/dative guards.
    // Native results have at most two bytes and fit std::string SSO. Return
    // borrowed literal bytes here; the sole caller constructs canonical owned
    // text at its output boundary. This removes Rust's avoidable intermediate
    // heap String without changing any token, byte order or source branch.
    // Cost ✔️: one bounded switch, O(1) and no heap allocation in either
    // Native's short result or this kernel; no scan, clone or dynamic lookup.
    match bond_order {
        BondOrder::Single | BondOrder::Aromatic if params.isomeric_smiles => match direction {
            BondDirection::EndDownRight => "\\",
            BondDirection::EndUpRight => "/",
            _ if bond_order == BondOrder::Single => "-",
            _ => ":",
        },
        BondOrder::Single => "-",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Quadruple => "$",
        BondOrder::Aromatic => ":",
        BondOrder::Dative if params.include_dative_bonds && reverse_dative => "<-",
        BondOrder::Dative if params.include_dative_bonds => "->",
        BondOrder::Dative => "-",
        BondOrder::Zero => "~",
        _ => "",
    }
}

fn source_reverse_dative(
    bond: &Bond,
    atom_to_left_idx: Option<usize>,
) -> Result<bool, SmartsWriteError> {
    // RDKit❗✔️:     bool reverseDative =
    // RDKit❗✔️:         (atomToLeftIdx >= 0 &&
    // RDKit❗✔️:          bond->getBeginAtomIdx() != static_cast<unsigned int>(atomToLeftIdx));
    // Sole detached transport of the identical source expression reached by
    // simple-query and non-query writers. Absent/negative left short-circuits
    // before the unsigned begin read. Source-width errors stay structural.
    // O(1), no allocation, clone, traversal or property lookup.
    if let Some(atom_idx) = atom_to_left_idx {
        let left = i32::try_from(atom_idx)
            .map_err(|_| SmartsWriteError::SourceAtomToLeftIndex { index: atom_idx })?;
        let begin = u32::try_from(bond.begin().index()).map_err(|_| {
            SmartsWriteError::SourceBondBeginIndex {
                index: bond.begin().index(),
            }
        })?;
        Ok(begin != left as u32)
    } else {
        Ok(false)
    }
}

fn get_bond_smarts_simple(
    bond: &Bond,
    query: &BondQueryPredicate,
    atom_to_left_idx: Option<usize>,
    params: &SmartsWriteParams,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string getBondSmartsSimple(const Bond *bond,
    // RDKit❗❌:                                 const QueryBond::QUERYBOND_QUERY *bquery,
    // RDKit❗❌:                                 int atomToLeftIdx,
    // RDKit❗❌:                                 const SmilesWriteParams &params) {
    // RDKit❗❌:   PRECONDITION(bond, "bad bond");
    // RDKit❗❌:   PRECONDITION(bquery, "bad query");
    // RDKit❗❌:
    // RDKit❗❌:   auto *equery = dynamic_cast<const BOND_EQUALS_QUERY *>(bquery);
    // RDKit❗❌:
    // RDKit❗❌:   std::string descrip = bquery->getDescription();
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:   if (descrip == "BondNull") {
    // RDKit❗❌:     res += "~";
    // RDKit❗❌:   } else if (descrip == "BondInRing") {
    // RDKit❗❌:     res += "@";
    // RDKit❗❌:   } else if (descrip == "SingleOrAromaticBond") {
    // RDKit❗❌:     auto dir = bond->getBondDir();
    // RDKit❗❌:     switch (dir) {
    // RDKit❗❌:       case Bond::ENDDOWNRIGHT: {
    // RDKit❗❌:         res += "\\";
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:       case Bond::ENDUPRIGHT: {
    // RDKit❗❌:         res += "/";
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:       default:
    // RDKit❗❌:         break;
    // RDKit❗❌:     }
    // RDKit❗❌:   } else if (descrip == "SingleOrDoubleBond") {
    // RDKit❗❌:     res += "-,=";
    // RDKit❗❌:   } else if (descrip == "DoubleOrAromaticBond") {
    // RDKit❗❌:     res += "=,:";
    // RDKit❗❌:   } else if (descrip == "SingleOrDoubleOrAromaticBond") {
    // RDKit❗❌:     res += "-,=,:";
    // RDKit❗❌:   } else if (descrip == "BondDir" && equery) {
    // RDKit❗❌:     int val = equery->getVal();
    // RDKit❗❌:     if (val == static_cast<int>(Bond::ENDDOWNRIGHT)) {
    // RDKit❗❌:       res += "\\";
    // RDKit❗❌:     } else if (val == static_cast<int>(Bond::ENDUPRIGHT)) {
    // RDKit❗❌:       res += "/";
    // RDKit❗❌:     } else {
    // RDKit❗❌:       throw "Can't write smarts for this bond dir type";
    // RDKit❗❌:     }
    // RDKit❗❌:   } else if (descrip == "BondOrder" && equery) {
    // RDKit❗❌:     bool reverseDative =
    // RDKit❗❌:         (atomToLeftIdx >= 0 &&
    // RDKit❗❌:          bond->getBeginAtomIdx() != static_cast<unsigned int>(atomToLeftIdx));
    // RDKit❗❌:     res += getBasicBondRepr(static_cast<Bond::BondType>(equery->getVal()),
    // RDKit❗❌:                             bond->getBondDir(), reverseDative, params);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     std::stringstream msg;
    // RDKit❗❌:     msg << "Can't write smarts for this query bond type: " << descrip;
    // RDKit❗❌:     throw msg.str().c_str();
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // Preserve source query-class dispatch: the four canonical OrderIn shapes
    // identify source factories, not arbitrary set membership. No sorting,
    // duplicate removal, heuristic tokens or fallbacks for other identities.
    // SingleOrAromaticBond intentionally ignores doIsomericSmiles. BondDir
    // failures and known source query descriptions are source errors, distinct
    // from independent predicates lacking a modeled Native class identity.
    // Cost ❌: bounded comparisons/dispatch for successful variants, but owned
    // Vec text allocates for short outputs that fit Native std::string SSO.
    // Unknown independent list errors clone their actual list payload only.
    match query {
        BondQueryPredicate::Any => Ok("~".into()),
        BondQueryPredicate::IsInRing(_) => Ok("@".into()),
        BondQueryPredicate::OrderIn(orders)
            if orders.as_slice() == [BondOrder::Single, BondOrder::Aromatic] =>
        {
            Ok(match bond.direction() {
                BondDirection::EndDownRight => "\\".into(),
                BondDirection::EndUpRight => "/".into(),
                _ => PropertyText::new(),
            })
        }
        BondQueryPredicate::OrderIn(orders)
            if orders.as_slice() == [BondOrder::Single, BondOrder::Double] =>
        {
            Ok("-,=".into())
        }
        BondQueryPredicate::OrderIn(orders)
            if orders.as_slice() == [BondOrder::Double, BondOrder::Aromatic] =>
        {
            Ok("=,:".into())
        }
        BondQueryPredicate::OrderIn(orders)
            if orders.as_slice() == [BondOrder::Single, BondOrder::Double, BondOrder::Aromatic] =>
        {
            Ok("-,=,:".into())
        }
        BondQueryPredicate::Direction(BondDirection::EndDownRight) => Ok("\\".into()),
        BondQueryPredicate::Direction(BondDirection::EndUpRight) => Ok("/".into()),
        BondQueryPredicate::Direction(direction) => Err(SmartsWriteError::SourceBondDirection {
            direction: *direction,
        }),
        BondQueryPredicate::Order(order) => {
            let reverse_dative = source_reverse_dative(bond, atom_to_left_idx)?;
            Ok(get_basic_bond_repr(*order, bond.direction(), reverse_dative, params).into())
        }
        BondQueryPredicate::HasStereo
        | BondQueryPredicate::NumRingBonds(_)
        | BondQueryPredicate::InRingOfSize(_)
        | BondQueryPredicate::MinRingSize(_)
        | BondQueryPredicate::HasProperty(_)
        | BondQueryPredicate::PropertyValue { .. } => {
            let description = match query {
                BondQueryPredicate::HasStereo => "BondStereo",
                BondQueryPredicate::NumRingBonds(_) => "BondInNRings",
                BondQueryPredicate::InRingOfSize(_) => "BondRingSize",
                BondQueryPredicate::MinRingSize(_) => "BondMinRingSize",
                BondQueryPredicate::HasProperty(_) => "HasProp",
                BondQueryPredicate::PropertyValue { .. } => "HasPropWithValue",
                _ => unreachable!(),
            };
            Err(SmartsWriteError::UnwritableBondQuery { description })
        }
        _ => Err(SmartsWriteError::UnsupportedBondQuery {
            predicate: query.clone(),
        }),
    }
}

fn atom_query_without_not<'a>(
    mut node: &'a QueryNode<AtomQueryPredicate>,
    mut negate: bool,
) -> (&'a QueryNode<AtomQueryPredicate>, bool) {
    while let QueryNode::Not(child) = node {
        negate = !negate;
        node = child;
    }
    (node, negate)
}

fn recurse_get_smarts<F>(
    atom: &QueryAtom,
    node: &QueryNode<AtomQueryPredicate>,
    negate: bool,
    features: &mut QueryBoolFeatures,
    params: &SmartsWriteParams,
    stereo_written: &mut bool,
    write_molecule: &mut F,
) -> Result<PropertyText, SmartsWriteError>
where
    F: FnMut(&crate::QueryGraph, &SmartsWriteParams) -> Result<PropertyText, SmartsWriteError>,
{
    recurse_get_smarts_with_owner(
        atom,
        node,
        negate,
        features,
        params,
        stereo_written,
        write_molecule,
        params.isomeric_smiles,
    )
}

fn recurse_get_smarts_with_owner<F>(
    atom: &QueryAtom,
    node: &QueryNode<AtomQueryPredicate>,
    negate: bool,
    features: &mut QueryBoolFeatures,
    params: &SmartsWriteParams,
    stereo_written: &mut bool,
    write_molecule: &mut F,
    owner_iso: bool,
) -> Result<PropertyText, SmartsWriteError>
where
    F: FnMut(&crate::QueryGraph, &SmartsWriteParams) -> Result<PropertyText, SmartsWriteError>,
{
    // RDKit❗❌: std::string _recurseGetSmarts(const QueryAtom *qatom,
    // RDKit❗❌:                               const QueryAtom::QUERYATOM_QUERY *node,
    // RDKit❗❌:                               bool negate, unsigned int &features,
    // RDKit❗❌:                               const SmilesWriteParams &params) {
    // RDKit❗❌:   PRECONDITION(node, "bad node");
    // RDKit❗❌:   // the algorithm goes something like this
    // RDKit❗❌:   // - recursively get the smarts for the child queries
    // RDKit❗❌:   // - combine the child smarts using the following rules:
    // RDKit❗❌:   //      - if we are currently at an OR query, combine the subqueries with a
    // RDKit❗❌:   //      ",",
    // RDKit❗❌:   //        but only if neither of child smarts do not contain "," and ";"
    // RDKit❗❌:   //        This situation leads to a no smartable situation and throw an
    // RDKit❗❌:   //        error
    // RDKit❗❌:   //      - if we are currently at an and query, combine the child smarts with
    // RDKit❗❌:   //      "&"
    // RDKit❗❌:   //        if neither of the child smarts contain a "," - otherwise combine
    // RDKit❗❌:   //        them
    // RDKit❗❌:   //        the child smarts with a ";"
    // RDKit❗❌:   //
    // RDKit❗❌:   // There is an additional complication with composite nodes that carry a
    // RDKit❗❌:   // negation - in this
    // RDKit❗❌:   // case we will propagate the negation to the child nodes using the
    // RDKit❗❌:   // following rules
    // RDKit❗❌:   //   NOT (a AND b) = ( NOT (a)) AND ( NOT (b))
    // RDKit❗❌:   //   NOT (a OR b) = ( NOT (a)) OR ( NOT (b))
    // RDKit❗❌:
    // RDKit❗❌:   auto descrip = node->getDescription();
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int child1Features = 0;
    // RDKit❗❌:   unsigned int child2Features = 0;
    // RDKit❗❌:   auto chi = node->beginChildren();
    // RDKit❗❌:   auto child1 = chi->get();
    // RDKit❗❌:   auto dsc1 = child1->getDescription();
    // RDKit❗❌:
    // RDKit❗❌:   ++chi;
    // RDKit❗❌:   CHECK_INVARIANT(chi != node->endChildren(),
    // RDKit❗❌:                   "Not enough children on the query");
    // RDKit❗❌:
    // RDKit❗❌:   bool needParen;
    // RDKit❗❌:   std::string csmarts1;
    // RDKit❗❌:   // deal with the first child
    // RDKit❗❌:   if (dsc1 == "RecursiveStructure") {
    // RDKit❗❌:     csmarts1 = getRecursiveStructureQuerySmarts(child1, params);
    // RDKit❗❌:     features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_RECURSION);
    // RDKit❗❌:   } else if ((dsc1 != "AtomOr") && (dsc1 != "AtomAnd")) {
    // RDKit❗❌:     // child 1 is a simple node, but we only check for the smilesSymbol
    // RDKit❗❌:     //  if descrip=="AtomAnd"
    // RDKit❗❌:     csmarts1 = getAtomSmartsSimple(qatom, child1, needParen,
    // RDKit❗❌:                                    descrip == "AtomAnd", params);
    // RDKit❗❌:     bool nneg = (negate) ^ (child1->getNegation());
    // RDKit❗❌:     if (nneg) {
    // RDKit❗❌:       csmarts1 = "!" + csmarts1;
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // child 1 is composite node - recurse
    // RDKit❗❌:     bool nneg = (negate) ^ (child1->getNegation());
    // RDKit❗❌:     csmarts1 = _recurseGetSmarts(qatom, child1, nneg, child1Features, params);
    // RDKit❗❌:   }
    // RDKit❗❌:   // ok if we have a negation and we have an OR , we have to change to
    // RDKit❗❌:   // an AND since we propagated the negation
    // RDKit❗❌:   // i.e NOT (A OR B) = (NOT (A)) AND (NOT(B))
    // RDKit❗❌:   if (negate) {
    // RDKit❗❌:     if (descrip == "AtomOr") {
    // RDKit❗❌:       descrip = "AtomAnd";
    // RDKit❗❌:     } else if (descrip == "AtomAnd") {
    // RDKit❗❌:       descrip = "AtomOr";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   auto res = csmarts1;
    // RDKit❗❌:   while (chi != node->endChildren()) {
    // RDKit❗❌:     auto child2 = chi->get();
    // RDKit❗❌:     ++chi;
    // RDKit❗❌:
    // RDKit❗❌:     auto dsc2 = child2->getDescription();
    // RDKit❗❌:     std::string csmarts2;
    // RDKit❗❌:
    // RDKit❗❌:     // deal with the next child
    // RDKit❗❌:     if (dsc2 == "RecursiveStructure") {
    // RDKit❗❌:       csmarts2 = getRecursiveStructureQuerySmarts(child2, params);
    // RDKit❗❌:       features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_RECURSION);
    // RDKit❗❌:     } else if ((dsc2 != "AtomOr") && (dsc2 != "AtomAnd")) {
    // RDKit❗❌:       // child 2 is a simple node
    // RDKit❗❌:       csmarts2 = getAtomSmartsSimple(qatom, child2, needParen, false, params);
    // RDKit❗❌:       bool nneg = (negate) ^ (child2->getNegation());
    // RDKit❗❌:       if (nneg) {
    // RDKit❗❌:         csmarts2 = "!" + csmarts2;
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       bool nneg = (negate) ^ (child2->getNegation());
    // RDKit❗❌:       csmarts2 = _recurseGetSmarts(qatom, child2, nneg, child2Features, params);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     res = _combineChildSmarts(res, child1Features, csmarts2, child2Features,
    // RDKit❗❌:                               descrip, features);
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // Not nodes project source per-node negation flags; no predicate/child
    // reordering or boolean simplification. First-child symbol permission is
    // taken before the source De Morgan description swap.
    // Direct recursive children use ONLY their own negation, and set the
    // caller's feature word only AFTER their writer succeeds. Composite child
    // words stay local until combine succeeds. Source child2Features persists
    // across later siblings; child1Features never becomes aggregate features.
    // Cost ❌: source prefix concatenation also copies accumulated output,
    // potentially quadratic for many siblings. Vec-backed short strings and
    // leaf formatting add allocation costs versus source SSO. Each represented
    // node/Not wrapper is visited without cloning graphs or sorting children.
    let (node, negate) = atom_query_without_not(node, negate);
    let (children, original_description) = match node {
        QueryNode::And(children) => (children, "AtomAnd"),
        QueryNode::Or(children) => (children, "AtomOr"),
        QueryNode::Xor(_) => return Err(SmartsWriteError::XorComposite),
        QueryNode::Predicate(_) | QueryNode::Not(_) => {
            return Err(SmartsWriteError::CompositeChildCount { kind: "atom" });
        }
    };
    if children.len() < 2 {
        return Err(SmartsWriteError::CompositeChildCount { kind: "atom" });
    }
    let render_child = |child: &QueryNode<AtomQueryPredicate>,
                        check_for_symbol: bool,
                        child_features: &mut QueryBoolFeatures,
                        stereo_written: &mut bool,
                        write_molecule: &mut F|
     -> Result<(PropertyText, bool), SmartsWriteError> {
        let (child, own_negation) = atom_query_without_not(child, false);
        match child {
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => Ok((
                get_recursive_structure_query_smarts(
                    recursive,
                    own_negation,
                    params,
                    write_molecule,
                )?,
                true,
            )),
            QueryNode::Predicate(predicate) => {
                let mut need_paren = false;
                let mut result = get_atom_smarts_simple(
                    atom,
                    predicate,
                    &mut need_paren,
                    check_for_symbol,
                    owner_iso,
                    stereo_written,
                )?;
                if negate ^ own_negation {
                    result.insert_byte(0, b'!');
                }
                Ok((result, false))
            }
            QueryNode::And(_) | QueryNode::Or(_) | QueryNode::Xor(_) | QueryNode::Not(_) => Ok((
                recurse_get_smarts_with_owner(
                    atom,
                    child,
                    negate ^ own_negation,
                    child_features,
                    params,
                    stereo_written,
                    write_molecule,
                    owner_iso,
                )?,
                false,
            )),
        }
    };
    let mut first_features = QueryBoolFeatures::default();
    let mut second_features = QueryBoolFeatures::default();
    let (mut result, first_recursive) = render_child(
        &children[0],
        original_description == "AtomAnd",
        &mut first_features,
        stereo_written,
        write_molecule,
    )?;
    if first_recursive {
        features.insert(QueryBoolFeatures::HAS_RECURSION);
    }
    let description = if negate {
        if original_description == "AtomOr" {
            "AtomAnd"
        } else {
            "AtomOr"
        }
    } else {
        original_description
    };
    for child in &children[1..] {
        let (child_smarts, recursive) = render_child(
            child,
            false,
            &mut second_features,
            stereo_written,
            write_molecule,
        )?;
        if recursive {
            features.insert(QueryBoolFeatures::HAS_RECURSION);
        }
        result = combine_child_smarts(
            result,
            first_features,
            child_smarts,
            second_features,
            description,
            features,
        )?;
    }
    Ok(result)
}

fn bond_query_without_not<'a>(
    mut node: &'a QueryNode<BondQueryPredicate>,
    mut negate: bool,
) -> (&'a QueryNode<BondQueryPredicate>, bool) {
    while let QueryNode::Not(child) = node {
        negate = !negate;
        node = child;
    }
    (node, negate)
}

fn recurse_bond_smarts(
    bond: &Bond,
    node: &QueryNode<BondQueryPredicate>,
    negate: bool,
    atom_to_left_idx: Option<usize>,
    features: &mut QueryBoolFeatures,
    params: &SmartsWriteParams,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string _recurseBondSmarts(const Bond *bond,
    // RDKit❗❌:                                const QueryBond::QUERYBOND_QUERY *node,
    // RDKit❗❌:                                bool negate, int atomToLeftIdx,
    // RDKit❗❌:                                unsigned int &features,
    // RDKit❗❌:                                const SmilesWriteParams &params) {
    // RDKit❗❌:   // the algorithm goes something like this
    // RDKit❗❌:   // - recursively get the smarts for the child query bonds
    // RDKit❗❌:   // - combine the child smarts using the following rules:
    // RDKit❗❌:   //      - if we are currently at an OR query, combine the subqueries with a
    // RDKit❗❌:   //      ",",
    // RDKit❗❌:   //        but only if neither of child smarts do not contain "," and ";"
    // RDKit❗❌:   //        This situation leads to a no smartable situation and throw an
    // RDKit❗❌:   //        error
    // RDKit❗❌:   //      - if we are currently at an and query, combine the child smarts with
    // RDKit❗❌:   //      "&"
    // RDKit❗❌:   //        if neither of the child smarts contain a "," - otherwise combine
    // RDKit❗❌:   //        them
    // RDKit❗❌:   //        the child smarts with a ";"
    // RDKit❗❌:   //
    // RDKit❗❌:   // There is an additional complication with composite nodes that carry a
    // RDKit❗❌:   // negation - in this
    // RDKit❗❌:   // case we will propagate the negation to the child nodes using the
    // RDKit❗❌:   // following rules
    // RDKit❗❌:   //   NOT (a AND b) = ( NOT (a)) AND ( NOT (b))
    // RDKit❗❌:   //   NOT (a OR b) = ( NOT (a)) OR ( NOT (b))
    // RDKit❗❌:   PRECONDITION(bond, "bad bond");
    // RDKit❗❌:   PRECONDITION(node, "bad node");
    // RDKit❗❌:   std::string descrip = node->getDescription();
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:
    // RDKit❗❌:   const QueryBond::QUERYBOND_QUERY *child1;
    // RDKit❗❌:   const QueryBond::QUERYBOND_QUERY *child2;
    // RDKit❗❌:   unsigned int child1Features = 0;
    // RDKit❗❌:   unsigned int child2Features = 0;
    // RDKit❗❌:   QueryBond::QUERYBOND_QUERY::CHILD_VECT_CI chi;
    // RDKit❗❌:
    // RDKit❗❌:   chi = node->beginChildren();
    // RDKit❗❌:   child1 = chi->get();
    // RDKit❗❌:   chi++;
    // RDKit❗❌:   child2 = chi->get();
    // RDKit❗❌:   chi++;
    // RDKit❗❌:   // OK we should be at the end of vector by now - since we can have only two
    // RDKit❗❌:   // children,
    // RDKit❗❌:   // well - at least in this case
    // RDKit❗❌:   CHECK_INVARIANT(chi == node->endChildren(), "Too many children on the query");
    // RDKit❗❌:
    // RDKit❗❌:   std::string dsc1, dsc2;
    // RDKit❗❌:   dsc1 = child1->getDescription();
    // RDKit❗❌:   dsc2 = child2->getDescription();
    // RDKit❗❌:   std::string csmarts1, csmarts2;
    // RDKit❗❌:
    // RDKit❗❌:   if ((dsc1 != "BondOr") && (dsc1 != "BondAnd")) {
    // RDKit❗❌:     // child1 is  simple node get the smarts directly
    // RDKit❗❌:     const auto *tchild = static_cast<const BOND_EQUALS_QUERY *>(child1);
    // RDKit❗❌:     csmarts1 = getBondSmartsSimple(bond, tchild, atomToLeftIdx, params);
    // RDKit❗❌:     bool nneg = (negate) ^ (tchild->getNegation());
    // RDKit❗❌:     if (nneg) {
    // RDKit❗❌:       csmarts1 = "!" + csmarts1;
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // child1 is a composite node recurse further
    // RDKit❗❌:     bool nneg = (negate) ^ (child1->getNegation());
    // RDKit❗❌:     csmarts1 = _recurseBondSmarts(bond, child1, nneg, atomToLeftIdx,
    // RDKit❗❌:                                   child1Features, params);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // now deal with the second child
    // RDKit❗❌:   if ((dsc2 != "BondOr") && (dsc2 != "BondAnd")) {
    // RDKit❗❌:     // child 2 is a simple node
    // RDKit❗❌:     const auto *tchild = static_cast<const BOND_EQUALS_QUERY *>(child2);
    // RDKit❗❌:     csmarts2 = getBondSmartsSimple(bond, tchild, atomToLeftIdx, params);
    // RDKit❗❌:     bool nneg = (negate) ^ (tchild->getNegation());
    // RDKit❗❌:     if (nneg) {
    // RDKit❗❌:       csmarts2 = "!" + csmarts2;
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // child two is a composite node - recurse
    // RDKit❗❌:     bool nneg = (negate) ^ (child2->getNegation());
    // RDKit❗❌:     csmarts1 = _recurseBondSmarts(bond, child2, nneg, atomToLeftIdx,
    // RDKit❗❌:                                   child2Features, params);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // ok if we have a negation and we have to change the underlying logic,
    // RDKit❗❌:   // since we propagated the negation i.e NOT (A OR B) = (NOT (A)) AND
    // RDKit❗❌:   // (NOT(B))
    // RDKit❗❌:   if (negate) {
    // RDKit❗❌:     if (descrip == "BondOr") {
    // RDKit❗❌:       descrip = "BondAnd";
    // RDKit❗❌:     } else if (descrip == "BondAnd") {
    // RDKit❗❌:       descrip = "BondOr";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   res += _combineChildSmarts(csmarts1, child1Features, csmarts2, child2Features,
    // RDKit❗❌:                              descrip, features);
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // Reuse source-exact binary recursion and its observable asymmetric second
    // composite assignment: overwrite child1 text, leave child2 text empty,
    // retain both child feature words. Do not correct or simplify the source.
    // Child writers/error effects run in physical order before connective
    // negation swap and combination. Not wrappers project source node flags.
    // Cost ❌: one visit per query/Not node, O(height) stack, with cumulative
    // copied child text proportional to output at each depth (not just O(n)).
    // Source does these prefix copies too, but Vec short strings/leaf writer
    // formatting retain known allocation costs versus Native SSO.
    let (node, negate) = bond_query_without_not(node, negate);
    let (children, mut description) = match node {
        QueryNode::And(children) => (children, "BondAnd"),
        QueryNode::Or(children) => (children, "BondOr"),
        QueryNode::Xor(_) => return Err(SmartsWriteError::XorComposite),
        QueryNode::Predicate(_) | QueryNode::Not(_) => {
            return Err(SmartsWriteError::CompositeChildCount { kind: "bond" });
        }
    };
    if children.len() != 2 {
        return Err(SmartsWriteError::CompositeChildCount { kind: "bond" });
    }

    let render_child = |child: &QueryNode<BondQueryPredicate>,
                        child_features: &mut QueryBoolFeatures|
     -> Result<PropertyText, SmartsWriteError> {
        let (child, child_negate) = bond_query_without_not(child, negate);
        match child {
            QueryNode::Predicate(predicate) => {
                let mut result = get_bond_smarts_simple(bond, predicate, atom_to_left_idx, params)?;
                if child_negate {
                    result.insert_byte(0, b'!');
                }
                Ok(result)
            }
            QueryNode::And(_) | QueryNode::Or(_) | QueryNode::Xor(_) | QueryNode::Not(_) => {
                recurse_bond_smarts(
                    bond,
                    child,
                    child_negate,
                    atom_to_left_idx,
                    child_features,
                    params,
                )
            }
        }
    };

    let mut child1_features = QueryBoolFeatures::default();
    let mut child2_features = QueryBoolFeatures::default();
    let mut child1_smarts = render_child(&children[0], &mut child1_features)?;
    let child2_base = bond_query_without_not(&children[1], false).0;
    let child2_smarts = if matches!(child2_base, QueryNode::And(_) | QueryNode::Or(_)) {
        // Preserve the pinned source assignment to csmarts1 in this branch.
        child1_smarts = render_child(&children[1], &mut child2_features)?;
        PropertyText::new()
    } else {
        render_child(&children[1], &mut child2_features)?
    };
    if negate {
        description = if description == "BondOr" {
            "BondAnd"
        } else {
            "BondOr"
        };
    }
    combine_child_smarts(
        child1_smarts,
        child1_features,
        child2_smarts,
        child2_features,
        description,
        features,
    )
}

/// Serialize a complete detached query graph, including cycles and recursive
/// SMARTS predicates.
pub fn write_smarts(
    graph: &QueryGraph,
    params: &SmartsWriteParams,
) -> Result<PropertyText, SmartsWriteError> {
    graph
        .validate()
        .map_err(|error| SmartsWriteError::InvalidGraph(error.to_string()))?;
    query_graph_to_smarts(graph, params)
}

struct NonQueryAtomSmartsValue<'a> {
    id: AtomId,
    atomic_number: u8,
    isotope: Option<u16>,
    explicit_hydrogens: u8,
    formal_charge: i8,
    chiral_tag: ChiralTag,
    atom_map: Option<u32>,
    props: &'a BTreeMap<PropertyText, cosmolkit_model::PropertyValue>,
}
impl<'a> NonQueryAtomSmartsValue<'a> {
    fn prop(&self, key: &str) -> Option<&'a cosmolkit_model::PropertyValue> {
        self.props.get(key.as_bytes())
    }
}
impl<'a> From<&'a cosmolkit_model::Atom> for NonQueryAtomSmartsValue<'a> {
    fn from(atom: &'a cosmolkit_model::Atom) -> Self {
        Self {
            id: atom.id(),
            atomic_number: atom.atomic_number(),
            isotope: atom.isotope(),
            explicit_hydrogens: atom.explicit_hydrogens(),
            formal_charge: atom.formal_charge(),
            chiral_tag: atom.chiral_tag(),
            atom_map: atom.atom_map(),
            props: atom.props(),
        }
    }
}
impl<'a> From<&'a QueryAtom> for NonQueryAtomSmartsValue<'a> {
    fn from(atom: &'a QueryAtom) -> Self {
        Self {
            id: atom.id(),
            atomic_number: atom.atomic_number(),
            isotope: atom.isotope(),
            explicit_hydrogens: atom.explicit_hydrogens(),
            formal_charge: atom.formal_charge(),
            chiral_tag: atom.chiral_tag(),
            atom_map: atom.atom_map(),
            props: atom.props(),
        }
    }
}
// The view borrows canonical properties and copies scalar getters only.
// No clone, query conversion, chemistry dispatch or ownership change occurs.

fn non_query_atom_to_smarts<'a>(
    atom: impl Into<NonQueryAtomSmartsValue<'a>>,
    owner_has_do_iso_smiles: bool,
    stereo_written: &mut bool,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❗: std::string getNonQueryAtomSmarts(const Atom *atom) {
    // RDKit❗❗:   PRECONDITION(atom, "bad atom");
    // RDKit❗❗:   PRECONDITION(!atom->hasQuery(), "atom should not have query");
    // RDKit❗❗:   std::stringstream res;
    // RDKit❗❗:   res << "[";
    // RDKit❗❗:
    // RDKit❗❗:   int isotope = atom->getIsotope();
    // RDKit❗❗:   if (isotope) {
    // RDKit❗❗:     res << isotope;
    // RDKit❗❗:   }
    // RDKit❗❗:
    // RDKit❗❗:   std::string symbol;
    // RDKit❗❗:   if (atom->getPropIfPresent(common_properties::smilesSymbol, symbol)) {
    // RDKit❗❗:     res << symbol;
    // RDKit❗❗:   } else if (SmilesWrite::inOrganicSubset(atom->getAtomicNum())) {
    // RDKit❗❗:     res << "#" << atom->getAtomicNum();
    // RDKit❗❗:   } else {
    // RDKit❗❗:     res << atom->getSymbol();
    // RDKit❗❗:   }
    // RDKit❗❗:   bool addedChirality = false;
    // RDKit❗❗:   if (atom->hasOwningMol() &&
    // RDKit❗❗:       atom->getOwningMol().hasProp(common_properties::_doIsoSmiles)) {
    // RDKit❗❗:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗❗:         !atom->hasProp(_qatomHasStereoSet) &&
    // RDKit❗❗:         !atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❗❗:       atom->setProp(_qatomHasStereoSet, 1);
    // RDKit❗❗:       switch (atom->getChiralTag()) {
    // RDKit❗❗:         case Atom::CHI_TETRAHEDRAL_CW:
    // RDKit❗❗:           res << "@@";
    // RDKit❗❗:           addedChirality = true;
    // RDKit❗❗:           break;
    // RDKit❗❗:         case Atom::CHI_TETRAHEDRAL_CCW:
    // RDKit❗❗:           res << "@";
    // RDKit❗❗:           addedChirality = true;
    // RDKit❗❗:           break;
    // RDKit❗❗:         default:
    // RDKit❗❗:           break;
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // RDKit❗❗:
    // RDKit❗❗:   if (addedChirality && atom->getNumExplicitHs() == 1) {
    // RDKit❗❗:     // FIX: this isn't really correct in many cases, but
    // RDKit❗❗:     //   fixing it requires opening a fairly large construction site on the
    // RDKit❗❗:     //   SMARTS handling side. We'll do this later.
    // RDKit❗❗:     res << "H";
    // RDKit❗❗:   }
    // RDKit❗❗:   auto chg = atom->getFormalCharge();
    // RDKit❗❗:   if (chg) {
    // RDKit❗❗:     if (chg == -1) {
    // RDKit❗❗:       res << "-";
    // RDKit❗❗:     } else if (chg == 1) {
    // RDKit❗❗:       res << "+";
    // RDKit❗❗:     } else if (chg < 0) {
    // RDKit❗❗:       res << atom->getFormalCharge();
    // RDKit❗❗:     } else {
    // RDKit❗❗:       res << "+" << atom->getFormalCharge();
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // RDKit❗❗:   int mapNum;
    // RDKit❗❗:   if (atom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
    // RDKit❗❗:     res << ":";
    // RDKit❗❗:     res << mapNum;
    // RDKit❗❗:   }
    // RDKit❗❗:   res << "]";
    // RDKit❗❗:   return res.str();
    // RDKit❗❗: }
    // Caller established source !hasQuery. Shared ordinary Atom/carrier
    // getters borrow the same properties without try_to_atom cloning.
    // The explicit owner-presence input
    // projects source hasOwningMol && owner.hasProp(_doIsoSmiles), not a parsed
    // property's boolean value. stereo_written projects the mutable scratch
    // marker; raw marker presence still prevents a second emission.
    // Stream directly into canonical bytes in source order, retaining its
    // narrow explicit-H stereo rule, even where source comments call it FIX.
    // Cost ❗: linear bytes and logarithmic MODEL property lookup vs Native
    // dictionary scan. Vec and stringstream/SSO have different short-buffer
    // allocation/growth regimes; no graph scans/clones or extra decimal buffer.
    let atom = atom.into();
    let mut result = PropertyText::from("[");
    if let Some(isotope) = atom.isotope.filter(|v| *v != 0) {
        write!(result, "{isotope}").expect("byte formatter");
    }
    if let Some(symbol) = atom.prop("smilesSymbol") {
        result.extend_bytes(cosmolkit_core::property_value_to_string(symbol)?.as_bytes());
    } else if cosmolkit_core::is_rdkit_organic_subset(atom.atomic_number) {
        write!(result, "#{}", atom.atomic_number).expect("byte formatter");
    } else {
        // RDKit❗❗: std::string Atom::getSymbol() const {
        // RDKit❗❗:   std::string res;
        // RDKit❗❗:   // handle dummies differently:
        // RDKit❗❗:   if (d_atomicNum != 0 ||
        // RDKit❗❗:       !getPropIfPresent<std::string>(common_properties::dummyLabel, res)) {
        // RDKit❗❗:     res = PeriodicTable::getTable()->getElementSymbol(d_atomicNum);
        // RDKit❗❗:   }
        // RDKit❗❗:   return res;
        // RDKit❗❗: }
        if atom.atomic_number == 0
            && let Some(label) = atom.prop("dummyLabel")
        {
            result.extend_bytes(cosmolkit_core::property_value_to_string(label)?.as_bytes());
        } else {
            result
                .extend_bytes(cosmolkit_core::rdkit_element_symbol(atom.atomic_number)?.as_bytes());
        }
    }
    let mut added_chirality = false;
    if owner_has_do_iso_smiles
        && atom.chiral_tag != ChiralTag::Unspecified
        && !*stereo_written
        && atom.prop("_qatomHasStereoSet").is_none()
        && atom.prop("_brokenChirality").is_none()
    {
        *stereo_written = true;
        match atom.chiral_tag {
            ChiralTag::TetrahedralCw => {
                result.extend_bytes(b"@@");
                added_chirality = true;
            }
            ChiralTag::TetrahedralCcw => {
                result.push_byte(b'@');
                added_chirality = true;
            }
            _ => {}
        }
    }
    if added_chirality && atom.explicit_hydrogens == 1 {
        result.push_byte(b'H');
    }
    let charge = atom.formal_charge;
    match charge {
        0 => {}
        -1 => result.push_byte(b'-'),
        1 => result.push_byte(b'+'),
        i8::MIN..=-2 => write!(result, "{charge}").expect("byte formatter"),
        2..=i8::MAX => write!(result, "+{charge}").expect("byte formatter"),
    }
    // Canonical typed slot is authoritative; otherwise the actual raw source
    // property is read as int (unlike query atom writer's string conversion).
    let typed = atom.atom_map.map(cosmolkit_model::PropertyValue::UInt);
    if let Some(value) = typed.as_ref().or_else(|| atom.prop("molAtomMapNumber")) {
        let value = cosmolkit_core::property_value_to_int(value).map_err(|source| {
            SmartsWriteError::AtomMapInt {
                atom: atom.id,
                source,
            }
        })?;
        write!(result, ":{value}").expect("byte formatter");
    }
    result.push_byte(b']');
    Ok(result)
}

/// Serialize one detached query atom.
pub fn atom_to_smarts(
    graph: &QueryGraph,
    atom_id: AtomId,
    params: &SmartsWriteParams,
) -> Result<PropertyText, SmartsWriteError> {
    graph
        .validate()
        .map_err(|error| SmartsWriteError::InvalidGraph(error.to_string()))?;
    let atom =
        graph
            .atom(atom_id.index())
            .ok_or_else(|| SmartsWriteError::FragmentAtomOutOfRange {
                atom: atom_id.index(),
            })?;
    query_atom_to_smarts(atom, params)
}

fn non_query_bond_to_smarts(
    bond: &Bond,
    atom_to_left_idx: Option<usize>,
    params: &SmartsWriteParams,
) -> Result<&'static str, SmartsWriteError> {
    // RDKit❗✔️: std::string getNonQueryBondSmarts(const Bond *qbond, int atomToLeftIdx,
    // RDKit❗✔️:                                   const SmilesWriteParams &params) {
    // RDKit❗✔️:   PRECONDITION(qbond, "bad bond");
    // RDKit❗✔️:   std::string res;
    // RDKit❗✔️:
    // RDKit❗✔️:   if (qbond->getIsAromatic()) {
    // RDKit❗✔️:     res = ":";
    // RDKit❗✔️:     if (params.doIsomericSmiles) {
    // RDKit❗✔️:       if (qbond->getBondDir() == Bond::ENDDOWNRIGHT) {
    // RDKit❗✔️:         res = "\\";
    // RDKit❗✔️:       } else if (qbond->getBondDir() == Bond::ENDUPRIGHT) {
    // RDKit❗✔️:         res = "/";
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     bool reverseDative =
    // RDKit❗✔️:         (atomToLeftIdx >= 0 &&
    // RDKit❗✔️:          qbond->getBeginAtomIdx() != static_cast<unsigned int>(atomToLeftIdx));
    // RDKit❗✔️:     res = getBasicBondRepr(qbond->getBondType(), qbond->getBondDir(),
    // RDKit❗✔️:                            reverseDative, params);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // Aromatic carrier flag wins over the stored order before any left/index
    // read. Reuse the canonical basic-token aromatic switch and the single
    // source endpoint expression; never infer aromaticity from the order.
    // This function intentionally does not read a query payload: the Native
    // helper has no !hasQuery precondition and callers own query dispatch.
    // Cost ✔️: one scalar branch and bounded shared switch, O(1). Borrowed
    // result has no heap allocation, like Native's at-most-two-byte SSO text;
    // the caller acquires canonical owned output, no intermediate heap String.
    if bond.is_aromatic() {
        Ok(get_basic_bond_repr(
            BondOrder::Aromatic,
            bond.direction(),
            false,
            params,
        ))
    } else {
        let reverse = source_reverse_dative(bond, atom_to_left_idx)?;
        Ok(get_basic_bond_repr(
            bond.order(),
            bond.direction(),
            reverse,
            params,
        ))
    }
}

/// Serialize one detached query bond.
pub fn bond_to_smarts(
    graph: &QueryGraph,
    bond_id: BondId,
) -> Result<PropertyText, SmartsWriteError> {
    graph
        .validate()
        .map_err(|error| SmartsWriteError::InvalidGraph(error.to_string()))?;
    let bond =
        graph
            .bond(bond_id.index())
            .ok_or_else(|| SmartsWriteError::FragmentBondOutOfRange {
                bond: bond_id.index(),
            })?;
    query_bond_to_smarts(bond, &SmartsWriteParams::default(), None)
}

#[cfg(all(test, feature = "smiles-integration"))]
mod uint_complete_source_condition_cells {
    use super::*;
    fn query_graph(props: Vec<cosmolkit_model::PropertyValue>) -> cosmolkit_model::QueryGraph {
        let atoms: Vec<cosmolkit_model::QueryAtom> = (0..props.len() + 1)
            .map(|i| {
                cosmolkit_model::QueryAtom::new(
                    cosmolkit_model::AtomId::new(i),
                    cosmolkit_model::AtomSpec::new(cosmolkit_types::Element::C),
                )
            })
            .collect();
        let bonds = props
            .into_iter()
            .enumerate()
            .map(|(i, v)| {
                cosmolkit_model::QueryBond::new(
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
        let atom_count = atoms.len();
        cosmolkit_model::QueryGraph::from_parts(
            atoms,
            bonds,
            Vec::<(
                cosmolkit_model::PropertyText,
                cosmolkit_model::PropertyValue,
            )>::new(),
            vec![cosmolkit_model::Conformer2D::new(
                0,
                vec![[0.0, 0.0]; atom_count],
            )],
            vec![],
            vec![],
        )
        .unwrap()
    }

    // FROZEN UINT CONDITION: UNSIGNED_CFG_0
    #[test]
    fn uint_cell_unsigned_cfg_0_smarts_write() {
        let mut q = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop(
                "_MolFileBondCfg",
                cosmolkit_model::PropertyValue::UInt(0_u32),
            )
            .unwrap();
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                &cosmolkit_core::WedgeAssignments::default(),
                false
            ),
            Ok("".into())
        );
        assert_eq!(q, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_1
    #[test]
    fn uint_cell_unsigned_cfg_1_smarts_write() {
        let mut q = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop(
                "_MolFileBondCfg",
                cosmolkit_model::PropertyValue::UInt(1_u32),
            )
            .unwrap();
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                &cosmolkit_core::WedgeAssignments::default(),
                false
            ),
            Ok("wU:0.0".into())
        );
        assert_eq!(q, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_2
    #[test]
    fn uint_cell_unsigned_cfg_2_smarts_write() {
        let mut q = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop(
                "_MolFileBondCfg",
                cosmolkit_model::PropertyValue::UInt(2_u32),
            )
            .unwrap();
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                &cosmolkit_core::WedgeAssignments::default(),
                false
            ),
            Ok("w:0.0".into())
        );
        assert_eq!(q, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_3
    #[test]
    fn uint_cell_unsigned_cfg_3_smarts_write() {
        let mut q = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop(
                "_MolFileBondCfg",
                cosmolkit_model::PropertyValue::UInt(3_u32),
            )
            .unwrap();
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                &cosmolkit_core::WedgeAssignments::default(),
                false
            ),
            Ok("wD:0.0".into())
        );
        assert_eq!(q, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_4
    #[test]
    fn uint_cell_unsigned_cfg_4_smarts_write() {
        let mut q = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop(
                "_MolFileBondCfg",
                cosmolkit_model::PropertyValue::UInt(4_u32),
            )
            .unwrap();
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                &cosmolkit_core::WedgeAssignments::default(),
                false
            ),
            Ok("".into())
        );
        assert_eq!(q, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_255
    #[test]
    fn uint_cell_unsigned_cfg_255_smarts_write() {
        let mut q = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop(
                "_MolFileBondCfg",
                cosmolkit_model::PropertyValue::UInt(255_u32),
            )
            .unwrap();
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                &cosmolkit_core::WedgeAssignments::default(),
                false
            ),
            Ok("".into())
        );
        assert_eq!(q, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_256
    #[test]
    fn uint_cell_unsigned_cfg_256_smarts_write() {
        let mut q = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop(
                "_MolFileBondCfg",
                cosmolkit_model::PropertyValue::UInt(256_u32),
            )
            .unwrap();
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                &cosmolkit_core::WedgeAssignments::default(),
                false
            ),
            Ok("".into())
        );
        assert_eq!(q, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_2147483647
    #[test]
    fn uint_cell_unsigned_cfg_2147483647_smarts_write() {
        let mut q = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop(
                "_MolFileBondCfg",
                cosmolkit_model::PropertyValue::UInt(2147483647_u32),
            )
            .unwrap();
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                &cosmolkit_core::WedgeAssignments::default(),
                false
            ),
            Ok("".into())
        );
        assert_eq!(q, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_2147483648
    #[test]
    fn uint_cell_unsigned_cfg_2147483648_smarts_write() {
        let mut q = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop(
                "_MolFileBondCfg",
                cosmolkit_model::PropertyValue::UInt(2147483648_u32),
            )
            .unwrap();
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                &cosmolkit_core::WedgeAssignments::default(),
                false
            ),
            Ok("".into())
        );
        assert_eq!(q, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CFG_4294967295
    #[test]
    fn uint_cell_unsigned_cfg_4294967295_smarts_write() {
        let mut q = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop(
                "_MolFileBondCfg",
                cosmolkit_model::PropertyValue::UInt(4294967295_u32),
            )
            .unwrap();
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &[AtomId::new(0), AtomId::new(1)],
                &[BondId::new(0)],
                true,
                &cosmolkit_core::WedgeAssignments::default(),
                false
            ),
            Ok("".into())
        );
        assert_eq!(q, before);
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod smarts_source_state_tests {
    use super::*;
    use cosmolkit_model::QueryAtomIdentity;

    fn graph(identity: u8, predicate: AtomQueryPredicate) -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::from_identity_parts(
                AtomId::new(0),
                QueryAtomIdentity::from_atomic_number(identity),
                QueryNode::predicate(predicate),
            )],
            vec![],
            Vec::<(
                cosmolkit_model::PropertyText,
                cosmolkit_model::PropertyValue,
            )>::new(),
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }

    #[test]
    fn periodic_table_precondition_uses_carrier_and_precedes_graph_emission() {
        // Independently constructed carrier/query states distinguish the
        // property's table key from numeric values in the predicate tree.
        let invalid = graph(u8::MAX, AtomQueryPredicate::Any);
        let expected =
            SmartsWriteError::Valence(cosmolkit_core::ValenceError::PeriodicTableLookup {
                atomic_number: u8::MAX,
                field: "valences",
            });
        assert_eq!(
            mol_to_smarts_source(
                &invalid,
                &Default::default(),
                vec![cosmolkit_smiles::AtomColor::White],
                &[true],
                None
            )
            .map(|_| ()),
            Err(expected.clone())
        );
        assert_eq!(
            query_atom_to_smarts(invalid.atom(0).unwrap(), &Default::default()).unwrap(),
            "*".into()
        );
        let before = invalid.clone();
        assert_eq!(
            query_graph_to_smarts(&invalid, &Default::default()),
            Err(expected.clone())
        );
        #[cfg(feature = "smiles-integration")]
        assert_eq!(
            query_graph_to_cx_smarts(&invalid, &Default::default()),
            Err(expected.clone())
        );
        assert_eq!(invalid, before);
        let valid = graph(0, AtomQueryPredicate::AtomicNumber(u8::MAX));
        assert!(
            mol_to_smarts_source(
                &valid,
                &Default::default(),
                vec![cosmolkit_smiles::AtomColor::White],
                &[true],
                None
            )
            .map(|_| ())
            .is_ok()
        );
        assert!(query_graph_to_smarts(&valid, &Default::default()).is_ok());

        let inner = RecursiveStructureQuery::from_query_graph(invalid, 0);
        let recursive = graph(0, AtomQueryPredicate::RecursiveSmarts(inner));
        assert_eq!(
            query_graph_to_smarts(&recursive, &Default::default()),
            Err(expected)
        );
    }
}

#[cfg(test)]
mod complete_atom_smarts_simple_source_tests {
    use super::*;
    use cosmolkit_model::{AtomRangeQuery, AtomSpec, PropertyValue};
    fn atom() -> QueryAtom {
        QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))
    }
    fn write(
        a: &QueryAtom,
        q: &AtomQueryPredicate,
        incoming: bool,
        check: bool,
        iso: bool,
        written: bool,
    ) -> Result<(PropertyText, bool, bool), SmartsWriteError> {
        let mut paren = incoming;
        let mut written = written;
        let value = get_atom_smarts_simple(a, q, &mut paren, check, iso, &mut written)?;
        Ok((value, paren, written))
    }
    fn text(q: AtomQueryPredicate, expected: &str, paren: bool) {
        for incoming in [false, true] {
            let (v, p, w) = write(&atom(), &q, incoming, true, false, false).unwrap();
            assert_eq!(v, PropertyText::from(expected));
            assert_eq!(p, paren);
            assert!(!w);
        }
    }
    #[test]
    fn scalar_descriptions_have_exact_tokens_values_and_bracket_assignments() {
        for (q, token, paren) in [
            (AtomQueryPredicate::Any, "*", false),
            (AtomQueryPredicate::IsAromatic(false), "A", false),
            (AtomQueryPredicate::IsAromatic(true), "a", false),
            (AtomQueryPredicate::ImplicitHydrogenCount(-3), "h-3", true),
            (AtomQueryPredicate::HasImplicitHydrogen, "h", true),
            (AtomQueryPredicate::TotalValence(0), "v0", true),
            (AtomQueryPredicate::ExplicitDegree(-1), "D-1", true),
            (AtomQueryPredicate::TotalDegree(5), "X5", true),
            (AtomQueryPredicate::HasRingBond, "x", true),
            (AtomQueryPredicate::HydrogenCount(2), "H2", true),
            (AtomQueryPredicate::InRing, "R", true),
            (AtomQueryPredicate::SmallestRingSize(-2), "r-2", true),
            (AtomQueryPredicate::InRingOfSize(6), "k6", true),
            (AtomQueryPredicate::NumAtomRings(-1), "R", true),
            (AtomQueryPredicate::NumAtomRings(0), "R0", true),
            (AtomQueryPredicate::NumAtomRings(2), "R2", true),
            (AtomQueryPredicate::HasHeteroatomNeighbors, "z", true),
            (AtomQueryPredicate::NumHeteroatomNeighbors(2), "z2", true),
            (
                AtomQueryPredicate::HasAliphaticHeteroatomNeighbors,
                "Z",
                true,
            ),
            (
                AtomQueryPredicate::NumAliphaticHeteroatomNeighbors(3),
                "Z3",
                true,
            ),
            (AtomQueryPredicate::RingBondCount(2), "x2", true),
            (AtomQueryPredicate::IsUnsaturated, "$(*=,:,#*)", true),
        ] {
            text(q, token, paren);
        }
    }
    #[test]
    fn atomic_number_symbol_suppression_is_presence_only_and_preserves_incoming_brackets() {
        for value in [0, 6, 255] {
            text(
                AtomQueryPredicate::AtomicNumber(value),
                &format!("#{value}"),
                true,
            );
            for prop in [
                PropertyValue::Bool(false),
                PropertyValue::UInt(0),
                PropertyValue::String("different symbol".into()),
            ] {
                let mut a = atom();
                a.set_prop("smilesSymbol", prop).unwrap();
                for incoming in [false, true] {
                    let (v, p, _) = write(
                        &a,
                        &AtomQueryPredicate::AtomicNumber(value),
                        incoming,
                        false,
                        false,
                        false,
                    )
                    .unwrap();
                    assert!(v.is_empty());
                    assert_eq!(p, incoming);
                }
            }
        }
    }
    #[test]
    fn organic_atom_types_leave_bracket_state_unchanged_and_nonorganic_force_it() {
        for (n, aromatic, token, organic) in [
            (0, false, "*", true),
            (5, false, "B", true),
            (6, true, "c", true),
            (7, false, "N", true),
            (17, false, "Cl", true),
            (35, true, "br", true),
            (14, true, "si", false),
            (34, true, "se", false),
            (118, true, "og", false),
            (1, false, "H", false),
        ] {
            for incoming in [false, true] {
                let (v, p, _) = write(
                    &atom(),
                    &AtomQueryPredicate::AtomType {
                        atomic_number: n,
                        aromatic,
                    },
                    incoming,
                    true,
                    false,
                    false,
                )
                .unwrap();
                assert_eq!(v, PropertyText::from(token));
                assert_eq!(p, if organic { incoming } else { true });
            }
        }
    }
    #[test]
    fn source_atom_type_encoding_threshold_and_periodic_table_errors_are_not_unsupported() {
        for incoming in [false, true] {
            for (n, aromatic, decoded) in [(0, true, 1000), (119, false, 119), (255, true, 255)] {
                let mut paren = incoming;
                let mut written = false;
                assert!(
                    matches!(get_atom_smarts_simple(&atom(),&AtomQueryPredicate::AtomType{atomic_number:n,aromatic},&mut paren,true,false,&mut written),Err(SmartsWriteError::AtomTypeAtomicNumber{atomic_number})if atomic_number==decoded)
                );
                assert_eq!(paren, incoming);
                assert!(!written);
            }
        }
        let mut a = atom();
        a.set_prop("smilesSymbol", PropertyValue::Bool(false))
            .unwrap();
        for (n, aromatic, expected) in [(0, true, "A"), (255, true, "a"), (255, false, "A")] {
            for incoming in [false, true] {
                let (v, p, _) = write(
                    &a,
                    &AtomQueryPredicate::AtomType {
                        atomic_number: n,
                        aromatic,
                    },
                    incoming,
                    true,
                    false,
                    false,
                )
                .unwrap();
                assert_eq!(v, PropertyText::from(expected));
                assert_eq!(p, incoming);
                assert!(
                    write(
                        &a,
                        &AtomQueryPredicate::AtomType {
                            atomic_number: n,
                            aromatic
                        },
                        incoming,
                        false,
                        false,
                        false
                    )
                    .is_err()
                );
            }
        }
    }
    #[test]
    fn unsigned_degree_storage_projects_source_signed_int_values_at_writer_boundary() {
        for (value, expected) in [
            (0, "d0"),
            (2147483647, "d2147483647"),
            (2147483648, "d-2147483648"),
            (u32::MAX, "d-1"),
        ] {
            text(AtomQueryPredicate::NonHydrogenDegree(value), expected, true);
        }
        text(
            AtomQueryPredicate::NonHydrogenDegreeLessEqual(u32::MAX),
            "d{-1-}",
            true,
        );
        text(
            AtomQueryPredicate::NonHydrogenDegreeGreaterEqual(2147483648),
            "d{--2147483648}",
            true,
        );
    }
    #[test]
    fn formal_charge_sign_unit_magnitude_and_source_undefined_abs_boundary_are_explicit() {
        for (q, expected) in [
            (AtomQueryPredicate::FormalCharge(-5), "-5"),
            (AtomQueryPredicate::FormalCharge(-1), "-"),
            (AtomQueryPredicate::FormalCharge(0), "+0"),
            (AtomQueryPredicate::FormalCharge(1), "+"),
            (AtomQueryPredicate::FormalCharge(2), "+2"),
            (AtomQueryPredicate::NegativeFormalCharge(-2), "+2"),
            (AtomQueryPredicate::NegativeFormalCharge(-1), "+"),
            (AtomQueryPredicate::NegativeFormalCharge(0), "-0"),
            (AtomQueryPredicate::NegativeFormalCharge(1), "-"),
            (AtomQueryPredicate::NegativeFormalCharge(5), "-5"),
        ] {
            text(q, expected, true);
        }
        for q in [
            AtomQueryPredicate::FormalCharge(i32::MIN),
            AtomQueryPredicate::NegativeFormalCharge(i32::MIN),
        ] {
            assert!(matches!(
                write(&atom(), &q, false, true, false, false),
                Err(SmartsWriteError::ChargeMagnitudeOverflow { value: i32::MIN })
            ));
        }
    }
    #[test]
    fn ranges_preserve_signed_bounds_openness_independence_and_equality_base_charge_signs() {
        for (function, prefix) in [
            (AtomRangeDataFunction::ExplicitDegree, "D"),
            (AtomRangeDataFunction::NonHydrogenDegree, "d"),
            (AtomRangeDataFunction::TotalDegree, "X"),
            (AtomRangeDataFunction::TotalValence, "v"),
            (AtomRangeDataFunction::NumAtomRings, "R"),
            (AtomRangeDataFunction::NumHeteroatomNeighbors, "z"),
            (AtomRangeDataFunction::NumAliphaticHeteroatomNeighbors, "Z"),
            (AtomRangeDataFunction::MinRingSize, "r"),
            (AtomRangeDataFunction::RingBondCount, "x"),
            (AtomRangeDataFunction::ImplicitHydrogenCount, "h"),
            (
                AtomRangeDataFunction::AtomRingSize {
                    lower: 0,
                    upper: 99,
                    lower_open: true,
                    upper_open: false,
                },
                "k",
            ),
        ] {
            for (bounds, suffix) in [
                (AtomRangeBounds::LessEqual(-3), "{-3-}"),
                (AtomRangeBounds::GreaterEqual(-3), "{--3}"),
                (
                    AtomRangeBounds::Inclusive {
                        lower: 8,
                        upper: -2,
                        lower_open: true,
                        upper_open: false,
                    },
                    "{8--2}",
                ),
            ] {
                text(
                    AtomQueryPredicate::Range(AtomRangeQuery::new(bounds, function)),
                    &format!("{prefix}{suffix}"),
                    true,
                );
            }
        }
        for (function, bounds, expected) in [
            (
                AtomRangeDataFunction::FormalCharge,
                AtomRangeBounds::LessEqual(-3),
                "-{-3-}",
            ),
            (
                AtomRangeDataFunction::NegativeFormalCharge,
                AtomRangeBounds::GreaterEqual(-3),
                "+{--3}",
            ),
            (
                AtomRangeDataFunction::FormalCharge,
                AtomRangeBounds::LessEqual(3),
                "+{3-}",
            ),
            (
                AtomRangeDataFunction::NegativeFormalCharge,
                AtomRangeBounds::GreaterEqual(3),
                "-{-3}",
            ),
            (
                AtomRangeDataFunction::FormalCharge,
                AtomRangeBounds::Inclusive {
                    lower: -5,
                    upper: -2,
                    lower_open: false,
                    upper_open: false,
                },
                "+{-5--2}",
            ),
            (
                AtomRangeDataFunction::NegativeFormalCharge,
                AtomRangeBounds::Inclusive {
                    lower: -5,
                    upper: -2,
                    lower_open: true,
                    upper_open: true,
                },
                "-{-5--2}",
            ),
        ] {
            text(
                AtomQueryPredicate::Range(AtomRangeQuery::new(bounds, function)),
                expected,
                true,
            );
        }
    }
    #[test]
    fn dedicated_comparison_variants_retain_literal_source_less_greater_tokens() {
        for (q, expected) in [
            (AtomQueryPredicate::ExplicitDegreeLessEqual(2), "D{2-}"),
            (AtomQueryPredicate::TotalDegreeLessEqual(2), "X{2-}"),
            (AtomQueryPredicate::TotalDegreeGreaterEqual(2), "X{-2}"),
            (AtomQueryPredicate::TotalValenceLessEqual(2), "v{2-}"),
            (AtomQueryPredicate::TotalValenceGreaterEqual(2), "v{-2}"),
            (AtomQueryPredicate::RingBondCountLessEqual(2), "x{2-}"),
            (
                AtomQueryPredicate::ImplicitHydrogenCountLessEqual(2),
                "h{2-}",
            ),
            (AtomQueryPredicate::InRingOfSizeLessEqual(2), "k{2-}"),
            (AtomQueryPredicate::InRingOfSizeGreaterEqual(2), "k{-2}"),
            (AtomQueryPredicate::SmallestRingSizeLessEqual(2), "r{2-}"),
            (AtomQueryPredicate::SmallestRingSizeGreaterEqual(2), "r{-2}"),
            (AtomQueryPredicate::DegreeLessEqual(2), "D{2-}"),
            (AtomQueryPredicate::DegreeGreaterEqual(2), "D{-2}"),
        ] {
            text(q, expected, true);
        }
    }
    #[test]
    fn hybridization_mass_and_isotope_follow_source_scalar_switch_and_mass_projection() {
        for (tag, expected) in [
            (Hybridization::S, "^0"),
            (Hybridization::Sp, "^1"),
            (Hybridization::Sp2, "^2"),
            (Hybridization::Sp3, "^3"),
            (Hybridization::Sp3d, "^4"),
            (Hybridization::Sp3d2, "^5"),
            (Hybridization::Unspecified, "^"),
            (Hybridization::Sp2d, "^"),
            (Hybridization::Other, "^"),
        ] {
            text(AtomQueryPredicate::HybridizationMatch(tag), expected, true);
        }
        text(AtomQueryPredicate::Mass(0), "0*", true);
        text(AtomQueryPredicate::Mass(65535), "65535*", true);
        text(AtomQueryPredicate::Isotope(-1), "-1*", true);
    }
    #[test]
    fn known_native_unknown_leaf_descriptions_use_literal_source_warning_fallback_and_keep_bracket_state()
     {
        for q in [
            AtomQueryPredicate::NumRadicalElectrons(2),
            AtomQueryPredicate::HasChiralTag,
            AtomQueryPredicate::MissingChiralTag,
            AtomQueryPredicate::ImplicitValence(3),
            AtomQueryPredicate::ExplicitValence(3),
            AtomQueryPredicate::HeavyAtomDegree(2),
            AtomQueryPredicate::IsBridgehead,
            AtomQueryPredicate::HasProperty("name".into()),
            AtomQueryPredicate::PropertyValue {
                name: "name".into(),
                value: "value".into(),
            },
        ] {
            for incoming in [false, true] {
                let (v, p, _) = write(&atom(), &q, incoming, true, false, false).unwrap();
                assert_eq!(v, PropertyText::from("*"));
                assert_eq!(p, incoming);
            }
        }
        for q in [
            AtomQueryPredicate::AtomicNumberIn(vec![6, 7]),
            AtomQueryPredicate::MolFileAlias("x".into()),
            AtomQueryPredicate::UnsupportedFeature("unmodeled identity"),
        ] {
            assert!(matches!(
                write(&atom(), &q, false, true, false, false),
                Err(SmartsWriteError::UnsupportedAtomQuery { .. })
            ));
        }
    }
    #[test]
    fn tetrahedral_emission_uses_source_property_presence_and_does_not_convert_marker_values() {
        for (tag, token) in [
            (ChiralTag::TetrahedralCw, "*@@"),
            (ChiralTag::TetrahedralCcw, "*@"),
        ] {
            for iso in [false, true] {
                for initial_written in [false, true] {
                    for guard in 0..3 {
                        let mut a = atom();
                        a.set_chiral_tag(tag);
                        if guard == 1 {
                            a.set_prop("_qatomHasStereoSet", PropertyValue::Bool(false))
                                .unwrap();
                        }
                        if guard == 2 {
                            a.set_prop(
                                "_brokenChirality",
                                PropertyValue::String("not converted".into()),
                            )
                            .unwrap();
                        }
                        let (v, p, w) = write(
                            &a,
                            &AtomQueryPredicate::Any,
                            false,
                            true,
                            iso,
                            initial_written,
                        )
                        .unwrap();
                        let emits = iso && !initial_written && guard == 0;
                        assert_eq!(v, PropertyText::from(if emits { token } else { "*" }));
                        assert_eq!(p, emits);
                        assert_eq!(w, initial_written || emits);
                    }
                }
            }
        }
    }
    #[test]
    fn non_tetrahedral_tags_mark_stereo_before_switch_even_when_they_emit_no_token() {
        for tag in [
            ChiralTag::Other,
            ChiralTag::Tetrahedral,
            ChiralTag::Allene,
            ChiralTag::SquarePlanar,
            ChiralTag::TrigonalBipyramidal,
            ChiralTag::Octahedral,
        ] {
            let mut a = atom();
            a.set_chiral_tag(tag);
            let (v, p, w) = write(&a, &AtomQueryPredicate::Any, false, true, true, false).unwrap();
            assert_eq!(v, PropertyText::from("*"));
            assert!(!p && w);
            a.set_chiral_tag(ChiralTag::TetrahedralCw);
            let (v, p, w) = write(&a, &AtomQueryPredicate::Any, false, true, true, w).unwrap();
            assert_eq!(v, PropertyText::from("*"));
            assert!(!p && w);
        }
        let (v, _, w) = write(&atom(), &AtomQueryPredicate::Any, false, true, true, false).unwrap();
        assert_eq!(v, PropertyText::from("*"));
        assert!(!w);
    }
}

#[cfg(test)]
mod complete_recursive_structure_smarts_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue};
    use std::cell::Cell;
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
    #[test]
    fn missing_inner_graph_errors_before_writer_even_if_cached_source_text_exists() {
        let q = RecursiveStructureQuery::new().with_source_smarts("C");
        for negated in [false, true] {
            assert!(matches!(
                get_recursive_structure_query_smarts(
                    &q,
                    negated,
                    &Default::default(),
                    |_, _| panic!("must not read cached text or call writer")
                ),
                Err(SmartsWriteError::MissingRecursiveQueryMolecule)
            ));
        }
    }
    #[test]
    fn writer_borrows_exact_query_and_all_parameter_fields_once_before_wrapping() {
        let q = RecursiveStructureQuery::from_query_graph(graph(), 9);
        let params = SmartsWriteParams {
            include_atom_maps: false,
            isomeric_smiles: false,
            include_dative_bonds: false,
            rooted_at_atom: Some(77),
        };
        let calls = Cell::new(0);
        let out = get_recursive_structure_query_smarts(&q, false, &params, |g, p| {
            calls.set(calls.get() + 1);
            assert!(std::ptr::eq(g, q.query_graph().unwrap()));
            assert!(std::ptr::eq(p, &params));
            assert_eq!(*p, params);
            Ok("[C:7]".into())
        })
        .unwrap();
        assert_eq!(calls.get(), 1);
        assert_eq!(out, PropertyText::from("$([C:7])"));
    }
    #[test]
    fn negation_is_exactly_one_prefix_even_if_nested_smarts_has_its_own_negation() {
        let q = RecursiveStructureQuery::from_query_graph(graph(), 0);
        for (negated, expected) in [(false, "$(!$([N]))"), (true, "!$(!$([N]))")] {
            let out =
                get_recursive_structure_query_smarts(&q, negated, &Default::default(), |_, _| {
                    Ok("!$([N])".into())
                })
                .unwrap();
            assert_eq!(out, PropertyText::from(expected));
        }
    }
    #[test]
    fn writer_error_propagates_unchanged_with_one_call_and_no_cached_text_fallback() {
        let q =
            RecursiveStructureQuery::from_query_graph(graph(), 0).with_source_smarts("fallback");
        for negated in [false, true] {
            let calls = Cell::new(0);
            let out =
                get_recursive_structure_query_smarts(&q, negated, &Default::default(), |_, _| {
                    calls.set(calls.get() + 1);
                    Err(SmartsWriteError::UnknownCombination {
                        description: "exact error".into(),
                    })
                });
            assert_eq!(calls.get(), 1);
            assert!(
                matches!(out,Err(SmartsWriteError::UnknownCombination{description})if description=="exact error")
            );
        }
    }
    #[test]
    fn empty_graph_is_distinct_from_null_and_still_calls_canonical_writer() {
        let q = RecursiveStructureQuery::from_query_graph(
            QueryGraph::from_parts(vec![], vec![], [], vec![], vec![], vec![]).unwrap(),
            0,
        );
        let calls = Cell::new(0);
        let out = get_recursive_structure_query_smarts(&q, true, &Default::default(), |g, p| {
            calls.set(calls.get() + 1);
            assert!(g.atoms().is_empty());
            query_graph_to_smarts(g, p)
        })
        .unwrap();
        assert_eq!(calls.get(), 1);
        assert_eq!(out, PropertyText::from("!$()"));
    }
    #[test]
    fn source_std_string_wrapping_preserves_nul_and_non_utf8_bytes() {
        let q = RecursiveStructureQuery::from_query_graph(graph(), 0);
        let inner = PropertyText::from(vec![0xff, 0, b'!', b')']);
        let out = get_recursive_structure_query_smarts(&q, true, &Default::default(), |_, _| {
            Ok(inner.clone())
        })
        .unwrap();
        assert_eq!(out.as_ref(), &[b'!', b'$', b'(', 0xff, 0, b'!', b')', b')']);
    }
    #[test]
    fn serial_membership_and_source_provenance_do_not_replace_actual_query_graph() {
        let mut q = RecursiveStructureQuery::from_query_graph(graph(), u32::MAX)
            .with_source_smarts("not actual graph");
        q.insert_atom_index(-1);
        q.insert_atom_index(i32::MAX);
        let before = format!("{q:?}");
        let out = get_recursive_structure_query_smarts(
            &q,
            false,
            &Default::default(),
            query_graph_to_smarts,
        )
        .unwrap();
        assert_eq!(out, PropertyText::from("$([#6])"));
        assert_eq!(format!("{q:?}"), before);
    }
    #[test]
    fn represented_graph_properties_and_all_query_inputs_are_immutable() {
        let mut g = graph();
        g.atoms_mut()[0]
            .set_prop("molAtomMapNumber", PropertyValue::Int(8))
            .unwrap();
        let q = RecursiveStructureQuery::from_query_graph(g, 7);
        let before = format!("{q:?}");
        let out = get_recursive_structure_query_smarts(
            &q,
            true,
            &Default::default(),
            query_graph_to_smarts,
        )
        .unwrap();
        assert_eq!(out, PropertyText::from("!$([#6:8])"));
        assert_eq!(format!("{q:?}"), before);
    }
}

#[cfg(test)]
mod complete_combine_child_smarts_source_tests {
    use super::*;

    #[test]
    fn source_feature_words_and_connectives_cover_all_flag_pairs() {
        for a in 0..16u32 {
            for b in 0..16u32 {
                for or in [false, true] {
                    let mut flags = QueryBoolFeatures(0x8000_0010);
                    let out = combine_child_smarts(
                        "C".into(),
                        QueryBoolFeatures(a | 0x100),
                        "N".into(),
                        QueryBoolFeatures(b | 0x200),
                        if or { "AtomOr" } else { "AtomAnd" },
                        &mut flags,
                    );
                    if or && (a & 6 == 6 || b & 6 == 6) {
                        assert_eq!(out, Err(SmartsWriteError::OrAboveAndBelowAnd));
                        assert_eq!(flags.0, 0x8000_0010);
                    } else {
                        let (sep, added) = if or {
                            (",", 4)
                        } else if (a | b) & 4 != 0 {
                            (";", 2)
                        } else {
                            ("&", 1)
                        };
                        assert_eq!(out.unwrap(), PropertyText::from(format!("C{sep}N")));
                        assert_eq!(flags.0, 0x8000_0010 | 0x300 | a | b | added);
                    }
                }
            }
        }
    }

    #[test]
    fn empty_children_suppress_separator_but_keep_feature_effects() {
        for (a, b, want) in [("", "", ""), ("", "N", "N"), ("C", "", "C")] {
            for (desc, flag) in [("AtomOr", 4), ("AtomAnd", 1)] {
                let mut f = QueryBoolFeatures(8);
                assert_eq!(
                    combine_child_smarts(
                        a.into(),
                        QueryBoolFeatures(0),
                        b.into(),
                        QueryBoolFeatures(0),
                        desc,
                        &mut f
                    )
                    .unwrap(),
                    PropertyText::from(want)
                );
                assert_eq!(f.0, 8 | flag);
            }
        }
    }

    #[test]
    fn exact_description_first_occurrence_and_or_priority() {
        for (desc, sep) in [
            ("xOr", ","),
            ("xAnd", "&"),
            ("xAndOr", ","),
            ("OrxAnd", "&"),
            ("ÅOr", ","),
            ("x\0And", "&"),
        ] {
            let mut f = QueryBoolFeatures::default();
            assert_eq!(
                combine_child_smarts(
                    "C".into(),
                    QueryBoolFeatures(0),
                    "N".into(),
                    QueryBoolFeatures(0),
                    desc,
                    &mut f
                )
                .unwrap(),
                PropertyText::from(format!("C{sep}N"))
            );
        }
        for desc in [
            "", "Or", "And", "OrxOr", "AndxAnd", "atomor", "xor", "AtomXor",
        ] {
            let mut f = QueryBoolFeatures(0xffff_ffff);
            assert_eq!(
                combine_child_smarts(
                    "".into(),
                    QueryBoolFeatures(0),
                    "".into(),
                    QueryBoolFeatures(0),
                    desc,
                    &mut f
                ),
                Err(SmartsWriteError::UnknownCombination {
                    description: desc.into()
                })
            );
            assert_eq!(f.0, 0xffff_ffff);
        }
    }

    #[test]
    fn nonsmartable_or_rejects_before_empty_child_and_feature_updates() {
        for (a, b) in [(6, 0), (0, 6), (6, 6)] {
            let mut f = QueryBoolFeatures(0x1000);
            assert_eq!(
                combine_child_smarts(
                    "".into(),
                    QueryBoolFeatures(a),
                    "".into(),
                    QueryBoolFeatures(b),
                    "AtomOr",
                    &mut f
                ),
                Err(SmartsWriteError::OrAboveAndBelowAnd)
            );
            assert_eq!(f.0, 0x1000);
        }
        let mut f = QueryBoolFeatures(0);
        assert_eq!(
            combine_child_smarts(
                "C".into(),
                QueryBoolFeatures(6),
                "N".into(),
                QueryBoolFeatures(0),
                "AtomAnd",
                &mut f
            )
            .unwrap(),
            PropertyText::from("C;N")
        );
        assert_eq!(f.0, 6);
    }

    #[test]
    fn output_bytes_are_preserved_without_parsing_or_deduplication() {
        let a = PropertyText::from(vec![0xff, 0, b',', b';']);
        let b = PropertyText::from(vec![0xfe, b'!']);
        let mut f = QueryBoolFeatures::default();
        let out = combine_child_smarts(
            a,
            QueryBoolFeatures(0),
            b,
            QueryBoolFeatures(0),
            "AtomOr",
            &mut f,
        )
        .unwrap();
        assert_eq!(out.as_bytes(), &[0xff, 0, b',', b';', b',', 0xfe, b'!']);
        assert_eq!(f.0, 4);
    }

    #[test]
    fn incoming_flags_do_not_choose_separator_or_cause_rejection() {
        for (desc, want) in [("AtomOr", "C,C"), ("AtomAnd", "C&C")] {
            let mut f = QueryBoolFeatures(u32::MAX);
            assert_eq!(
                combine_child_smarts(
                    "C".into(),
                    QueryBoolFeatures(0),
                    "C".into(),
                    QueryBoolFeatures(0),
                    desc,
                    &mut f
                )
                .unwrap(),
                PropertyText::from(want)
            );
            assert_eq!(f.0, u32::MAX);
        }
    }
}

#[cfg(test)]
mod complete_recurse_get_smarts_source_tests {
    use super::*;
    use cosmolkit_model::AtomSpec;
    fn atom() -> QueryAtom {
        QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))
    }
    fn leaf(n: u8) -> QueryNode<AtomQueryPredicate> {
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(n))
    }
    fn not(n: QueryNode<AtomQueryPredicate>) -> QueryNode<AtomQueryPredicate> {
        QueryNode::Not(Box::new(n))
    }
    fn rec() -> QueryNode<AtomQueryPredicate> {
        QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(
            RecursiveStructureQuery::from_query_graph(
                QueryGraph::from_parts(vec![atom()], vec![], [], vec![], vec![], vec![]).unwrap(),
                9,
            ),
        ))
    }
    fn bad() -> QueryNode<AtomQueryPredicate> {
        QueryNode::Predicate(AtomQueryPredicate::UnsupportedFeature(
            "explicit test identity",
        ))
    }
    fn render(
        a: &QueryAtom,
        n: &QueryNode<AtomQueryPredicate>,
        neg: bool,
        flags: &mut QueryBoolFeatures,
    ) -> Result<PropertyText, SmartsWriteError> {
        recurse_get_smarts(
            a,
            n,
            neg,
            flags,
            &Default::default(),
            &mut false,
            &mut |_, _| Ok("R".into()),
        )
    }

    #[test]
    fn de_morgan_switches_connective_and_xors_each_simple_child() {
        for or in [false, true] {
            for neg in [false, true] {
                let n = if or {
                    QueryNode::Or(vec![leaf(6), leaf(7)])
                } else {
                    QueryNode::And(vec![leaf(6), leaf(7)])
                };
                let mut f = QueryBoolFeatures(0x1000);
                let sep = if or ^ neg { "," } else { "&" };
                let prefix = if neg { "!" } else { "" };
                assert_eq!(
                    render(&atom(), &n, neg, &mut f).unwrap(),
                    PropertyText::from(format!("{prefix}#6{sep}{prefix}#7"))
                );
                assert_eq!(f.0, 0x1000 | if or ^ neg { 4 } else { 1 });
            }
        }
    }
    #[test]
    fn child_negation_flags_and_nested_not_projection_preserve_order() {
        let n = QueryNode::And(vec![not(leaf(6)), not(not(leaf(7))), leaf(8)]);
        assert_eq!(
            render(&atom(), &n, false, &mut QueryBoolFeatures(0)).unwrap(),
            PropertyText::from("!#6&#7&#8")
        );
        assert_eq!(
            render(&atom(), &n, true, &mut QueryBoolFeatures(0)).unwrap(),
            PropertyText::from("#6,!#7,!#8")
        );
        let n = QueryNode::And(vec![not(QueryNode::Or(vec![leaf(6), leaf(7)])), leaf(8)]);
        assert_eq!(
            render(&atom(), &n, false, &mut QueryBoolFeatures(0)).unwrap(),
            PropertyText::from("!#6&!#7&#8")
        );
    }
    #[test]
    fn only_first_child_uses_original_and_symbol_permission() {
        let mut a = atom();
        a.set_prop("smilesSymbol", "X").unwrap();
        for (or, neg, want) in [
            (false, false, "A&N"),
            (false, true, "!A,!N"),
            (true, false, "C,N"),
            (true, true, "!C&!N"),
        ] {
            let n = if or {
                QueryNode::Or(vec![
                    QueryNode::Predicate(AtomQueryPredicate::AtomType {
                        atomic_number: 6,
                        aromatic: false,
                    }),
                    QueryNode::Predicate(AtomQueryPredicate::AtomType {
                        atomic_number: 7,
                        aromatic: false,
                    }),
                ])
            } else {
                QueryNode::And(vec![
                    QueryNode::Predicate(AtomQueryPredicate::AtomType {
                        atomic_number: 6,
                        aromatic: false,
                    }),
                    QueryNode::Predicate(AtomQueryPredicate::AtomType {
                        atomic_number: 7,
                        aromatic: false,
                    }),
                ])
            };
            assert_eq!(
                render(&a, &n, neg, &mut QueryBoolFeatures(0)).unwrap(),
                PropertyText::from(want)
            );
        }
    }
    #[test]
    fn recursive_children_ignore_parent_negation_but_keep_their_own_flag() {
        for first in [false, true] {
            for own in [false, true] {
                for parent in [false, true] {
                    let r = if own { not(rec()) } else { rec() };
                    let nodes = if first {
                        vec![r, leaf(6)]
                    } else {
                        vec![leaf(6), r]
                    };
                    let n = QueryNode::And(nodes);
                    let mut f = QueryBoolFeatures(0);
                    let recursive = if own { "!$(R)" } else { "$(R)" };
                    let scalar = if parent { "!#6" } else { "#6" };
                    let sep = if parent { "," } else { "&" };
                    let want = if first {
                        format!("{recursive}{sep}{scalar}")
                    } else {
                        format!("{scalar}{sep}{recursive}")
                    };
                    assert_eq!(
                        render(&atom(), &n, parent, &mut f).unwrap(),
                        PropertyText::from(want)
                    );
                    assert_eq!(f.0, 8 | if parent { 4 } else { 1 });
                }
            }
        }
    }
    #[test]
    fn source_second_child_feature_word_persists_across_all_later_siblings() {
        let n = QueryNode::And(vec![
            leaf(6),
            QueryNode::Or(vec![leaf(7), leaf(8)]),
            leaf(9),
            QueryNode::And(vec![leaf(15), leaf(16)]),
        ]);
        let mut f = QueryBoolFeatures(0);
        assert_eq!(
            render(&atom(), &n, false, &mut f).unwrap(),
            PropertyText::from("#6;#7,#8;#9;#15&#16")
        );
        assert_eq!(f.0, 7);
    }
    #[test]
    fn first_child_word_stays_fixed_instead_of_becoming_aggregate_output_features() {
        let n = QueryNode::And(vec![leaf(6), QueryNode::Or(vec![leaf(7), leaf(8)]), rec()]);
        let mut f = QueryBoolFeatures(0);
        assert_eq!(
            render(&atom(), &n, false, &mut f).unwrap(),
            PropertyText::from("#6;#7,#8;$(R)")
        );
        assert_eq!(f.0, 14);
        // A first composite's own feature word does remain in every combine.
        let n = QueryNode::And(vec![
            QueryNode::Or(vec![leaf(6), leaf(7)]),
            leaf(8),
            leaf(9),
        ]);
        assert_eq!(
            render(&atom(), &n, false, &mut QueryBoolFeatures(0)).unwrap(),
            PropertyText::from("#6,#7;#8;#9")
        );
    }
    #[test]
    fn first_composite_features_are_not_published_before_a_failed_later_leaf() {
        let n = QueryNode::And(vec![QueryNode::Or(vec![leaf(6), leaf(7)]), bad()]);
        let mut f = QueryBoolFeatures(0x100);
        assert!(matches!(
            render(&atom(), &n, false, &mut f),
            Err(SmartsWriteError::UnsupportedAtomQuery { .. })
        ));
        assert_eq!(f.0, 0x100);
    }
    #[test]
    fn recursion_feature_effect_occurs_after_writer_success_and_before_later_error() {
        let n = QueryNode::And(vec![rec(), bad()]);
        let mut f = QueryBoolFeatures(0x100);
        assert!(render(&atom(), &n, false, &mut f).is_err());
        assert_eq!(f.0, 0x108);
        for n in [
            QueryNode::And(vec![rec(), leaf(6)]),
            QueryNode::And(vec![leaf(6), rec()]),
        ] {
            let mut f = QueryBoolFeatures(0x100);
            let out = recurse_get_smarts(
                &atom(),
                &n,
                false,
                &mut f,
                &Default::default(),
                &mut false,
                &mut |_, _| Err(SmartsWriteError::OrAboveAndBelowAnd),
            );
            assert_eq!(out, Err(SmartsWriteError::OrAboveAndBelowAnd));
            assert_eq!(f.0, 0x100);
        }
        let n = QueryNode::Or(vec![
            rec(),
            QueryNode::And(vec![QueryNode::Or(vec![leaf(6), leaf(7)]), leaf(8)]),
        ]);
        let mut f = QueryBoolFeatures(0);
        assert_eq!(
            render(&atom(), &n, false, &mut f),
            Err(SmartsWriteError::OrAboveAndBelowAnd)
        );
        assert_eq!(f.0, 8);
    }
    #[test]
    fn malformed_composites_are_structural_errors_before_any_child_writer() {
        for n in [QueryNode::And(vec![]), QueryNode::Or(vec![rec()]), leaf(6)] {
            let mut f = QueryBoolFeatures(0x100);
            assert_eq!(
                recurse_get_smarts(
                    &atom(),
                    &n,
                    false,
                    &mut f,
                    &Default::default(),
                    &mut false,
                    &mut |_, _| panic!("unreached")
                ),
                Err(SmartsWriteError::CompositeChildCount { kind: "atom" })
            );
            assert_eq!(f.0, 0x100);
        }
    }
    #[test]
    fn recursive_callbacks_follow_child_order_forward_options_and_keep_graphs_immutable() {
        let n = QueryNode::And(vec![rec(), rec()]);
        let before = format!("{n:?}");
        let a = atom();
        let before_a = format!("{a:?}");
        let p = SmartsWriteParams {
            rooted_at_atom: Some(8),
            include_atom_maps: false,
            include_dative_bonds: false,
            isomeric_smiles: false,
        };
        let mut calls = 0;
        let out = recurse_get_smarts(
            &a,
            &n,
            false,
            &mut QueryBoolFeatures(0),
            &p,
            &mut false,
            &mut |g, params| {
                assert_eq!(g.atoms().len(), 1);
                assert!(std::ptr::eq(params, &p));
                calls += 1;
                Ok(PropertyText::from(vec![0xff, calls]))
            },
        )
        .unwrap();
        assert_eq!(calls, 2);
        assert_eq!(
            out.as_bytes(),
            &[b'$', b'(', 0xff, 1, b')', b'&', b'$', b'(', 0xff, 2, b')']
        );
        assert_eq!(format!("{n:?}"), before);
        assert_eq!(format!("{a:?}"), before_a);
    }
    #[test]
    fn source_stereo_mark_is_shared_across_leaves_and_survives_later_error() {
        let mut a = atom();
        a.set_chiral_tag(ChiralTag::TetrahedralCw);
        let n = QueryNode::And(vec![
            QueryNode::Predicate(AtomQueryPredicate::Any),
            QueryNode::Predicate(AtomQueryPredicate::Any),
        ]);
        let mut written = false;
        assert_eq!(
            recurse_get_smarts(
                &a,
                &n,
                false,
                &mut QueryBoolFeatures(0),
                &Default::default(),
                &mut written,
                &mut |_, _| panic!("no recursion")
            )
            .unwrap(),
            PropertyText::from("*@@&*")
        );
        assert!(written);
        let n = QueryNode::And(vec![QueryNode::Predicate(AtomQueryPredicate::Any), bad()]);
        let mut written = false;
        assert!(
            recurse_get_smarts(
                &a,
                &n,
                false,
                &mut QueryBoolFeatures(0),
                &Default::default(),
                &mut written,
                &mut |_, _| panic!("no recursion")
            )
            .is_err()
        );
        assert!(written);
    }
}

#[cfg(test)]
mod complete_non_query_atom_smarts_source_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, PropertyValue};
    fn atom(spec: AtomSpec) -> Atom {
        Atom::from_spec(AtomId::new(4), spec)
    }
    fn write(a: &Atom, iso: bool) -> Result<PropertyText, SmartsWriteError> {
        non_query_atom_to_smarts(a, iso, &mut false)
    }
    #[test]
    fn every_element_uses_source_atomic_number_or_element_symbol_in_brackets() {
        for z in 0..=118 {
            let e = Element::from_atomic_number(z).unwrap();
            let want = if [0, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53].contains(&z) {
                format!("[#{z}]")
            } else {
                format!("[{}]", e.symbol())
            };
            for aromatic in [false, true] {
                assert_eq!(
                    write(&atom(AtomSpec::new(e).with_aromatic(aromatic)), false).unwrap(),
                    PropertyText::from(&want)
                );
            }
        }
    }
    #[test]
    fn isotope_precedes_symbol_and_zero_is_not_written() {
        for (isotope, want) in [(0, "[#6]"), (13, "[13#6]"), (u16::MAX, "[65535#6]")] {
            assert_eq!(
                write(
                    &atom(AtomSpec::new(Element::C).with_isotope(isotope)),
                    false
                )
                .unwrap(),
                PropertyText::from(want)
            );
        }
    }
    #[test]
    fn symbol_override_preserves_raw_bytes_and_dummy_label_is_unreached() {
        let mut a = atom(AtomSpec::new(Element::DUMMY));
        a.set_prop("dummyLabel", PropertyValue::IntVector(vec![8]))
            .unwrap();
        assert_eq!(write(&a, false).unwrap(), PropertyText::from("[#0]"));
        a.set_prop(
            "smilesSymbol",
            PropertyValue::String(PropertyText::from(vec![0xff, 0, b'X'])),
        )
        .unwrap();
        assert_eq!(
            write(&a, false).unwrap().as_bytes(),
            &[b'[', 0xff, 0, b'X', b']']
        );
        a.set_prop("smilesSymbol", PropertyValue::Int(7)).unwrap();
        assert_eq!(write(&a, false).unwrap(), PropertyText::from("[7]"));
        a.set_prop("smilesSymbol", PropertyValue::String(PropertyText::new()))
            .unwrap();
        assert_eq!(write(&a, false).unwrap(), PropertyText::from("[]"));
    }
    #[test]
    fn source_stereo_and_explicit_hydrogen_rule_for_all_tags() {
        for (tag, token, marked) in [
            (ChiralTag::Unspecified, "", false),
            (ChiralTag::TetrahedralCw, "@@", true),
            (ChiralTag::TetrahedralCcw, "@", true),
            (ChiralTag::Other, "", true),
            (ChiralTag::Tetrahedral, "", true),
            (ChiralTag::Allene, "", true),
            (ChiralTag::SquarePlanar, "", true),
            (ChiralTag::TrigonalBipyramidal, "", true),
            (ChiralTag::Octahedral, "", true),
        ] {
            for h in [0, 1, 2, u8::MAX] {
                for iso in [false, true] {
                    let a = atom(
                        AtomSpec::new(Element::C)
                            .with_chiral_tag(tag)
                            .with_explicit_hydrogens(h),
                    );
                    let mut written = false;
                    let t = if iso { token } else { "" };
                    let hs = if !t.is_empty() && h == 1 { "H" } else { "" };
                    assert_eq!(
                        non_query_atom_to_smarts(&a, iso, &mut written).unwrap(),
                        PropertyText::from(format!("[#6{t}{hs}]"))
                    );
                    assert_eq!(written, iso && marked);
                }
            }
        }
    }
    #[test]
    fn raw_presence_guards_and_existing_scratch_prevent_stereo_and_h() {
        for key in ["_qatomHasStereoSet", "_brokenChirality"] {
            for value in [PropertyValue::Bool(false), PropertyValue::String("".into())] {
                let mut a = atom(
                    AtomSpec::new(Element::C)
                        .with_chiral_tag(ChiralTag::TetrahedralCw)
                        .with_explicit_hydrogens(1),
                );
                a.set_prop(key, value).unwrap();
                let mut written = false;
                assert_eq!(
                    non_query_atom_to_smarts(&a, true, &mut written).unwrap(),
                    PropertyText::from("[#6]")
                );
                assert!(!written);
            }
        }
        let a = atom(
            AtomSpec::new(Element::C)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_explicit_hydrogens(1),
        );
        let mut written = true;
        assert_eq!(
            non_query_atom_to_smarts(&a, true, &mut written).unwrap(),
            PropertyText::from("[#6]")
        );
        assert!(written);
    }
    #[test]
    fn all_signed_byte_charges_follow_literal_source_formatting() {
        for charge in i8::MIN..=i8::MAX {
            let text = match charge {
                0 => String::new(),
                -1 => "-".into(),
                1 => "+".into(),
                v if v < 0 => v.to_string(),
                v => format!("+{v}"),
            };
            assert_eq!(
                write(
                    &atom(AtomSpec::new(Element::C).with_formal_charge(charge)),
                    false
                )
                .unwrap(),
                PropertyText::from(format!("[#6{text}]"))
            );
        }
    }
    #[test]
    fn map_property_is_read_as_signed_int_and_typed_slot_has_precedence() {
        for (value, want) in [
            (PropertyValue::Int(i32::MIN), "[#6:-2147483648]"),
            (PropertyValue::UInt(i32::MAX as u32), "[#6:2147483647]"),
            (PropertyValue::String("-03 \t".into()), "[#6:-3]"),
            (PropertyValue::Int(0), "[#6:0]"),
        ] {
            let mut a = atom(AtomSpec::new(Element::C));
            a.set_prop("molAtomMapNumber", value).unwrap();
            assert_eq!(write(&a, false).unwrap(), PropertyText::from(want));
        }
        let mut a = atom(AtomSpec::new(Element::C).with_atom_map(9));
        a.set_prop("molAtomMapNumber", PropertyValue::Bool(false))
            .unwrap();
        assert_eq!(write(&a, false).unwrap(), PropertyText::from("[#6:9]"));
        a.set_atom_map(Some(u32::MAX));
        assert!(matches!(
            write(&a, false),
            Err(SmartsWriteError::AtomMapInt {
                source: cosmolkit_core::PropertyIntReadError::UnsignedOverflow { .. },
                ..
            })
        ));
    }
    #[test]
    fn map_conversion_failure_retains_prior_stereo_marker_but_inputs_remain_immutable() {
        for value in [
            PropertyValue::Double(8.0),
            PropertyValue::Bool(false),
            PropertyValue::String("bad".into()),
        ] {
            let mut a = atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw));
            a.set_prop("molAtomMapNumber", value).unwrap();
            let before = format!("{a:?}");
            let mut written = false;
            assert!(
                matches!(non_query_atom_to_smarts(&a,true,&mut written),Err(SmartsWriteError::AtomMapInt{atom,..})if atom==AtomId::new(4))
            );
            assert!(written);
            assert_eq!(format!("{a:?}"), before);
        }
    }
    #[test]
    fn complete_output_order_is_isotope_symbol_stereo_h_charge_map() {
        let mut a = atom(
            AtomSpec::new(Element::C)
                .with_isotope(13)
                .with_chiral_tag(ChiralTag::TetrahedralCcw)
                .with_explicit_hydrogens(1)
                .with_formal_charge(2)
                .with_atom_map(17),
        );
        a.set_prop("smilesSymbol", "Q").unwrap();
        let before = format!("{a:?}");
        let mut written = false;
        assert_eq!(
            non_query_atom_to_smarts(&a, true, &mut written).unwrap(),
            PropertyText::from("[13Q@H+2:17]")
        );
        assert!(written);
        assert_eq!(format!("{a:?}"), before);
    }
}

#[cfg(test)]
mod complete_query_atom_smarts_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue};
    fn atom(node: QueryNode<AtomQueryPredicate>) -> QueryAtom {
        let mut a = QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C));
        a.set_predicate(node);
        a
    }
    fn leaf(p: AtomQueryPredicate) -> QueryNode<AtomQueryPredicate> {
        QueryNode::Predicate(p)
    }
    fn num(n: u8) -> QueryNode<AtomQueryPredicate> {
        leaf(AtomQueryPredicate::AtomicNumber(n))
    }
    fn not(q: QueryNode<AtomQueryPredicate>) -> QueryNode<AtomQueryPredicate> {
        QueryNode::Not(Box::new(q))
    }
    fn rec() -> QueryNode<AtomQueryPredicate> {
        leaf(AtomQueryPredicate::RecursiveSmarts(
            RecursiveStructureQuery::from_query_graph(
                QueryGraph::from_parts(vec![], vec![], [], vec![], vec![], vec![]).unwrap(),
                0,
            ),
        ))
    }
    fn write(a: &QueryAtom) -> Result<PropertyText, SmartsWriteError> {
        query_atom_to_smarts(a, &Default::default())
    }
    #[test]
    fn leaf_source_brackets_are_determined_by_actual_leaf_writer() {
        for (p, want) in [
            (AtomQueryPredicate::Any, "*"),
            (AtomQueryPredicate::IsAromatic(true), "a"),
            (AtomQueryPredicate::IsAromatic(false), "A"),
            (AtomQueryPredicate::AtomicNumber(6), "[#6]"),
            (
                AtomQueryPredicate::AtomType {
                    atomic_number: 6,
                    aromatic: false,
                },
                "C",
            ),
            (
                AtomQueryPredicate::AtomType {
                    atomic_number: 26,
                    aromatic: false,
                },
                "[Fe]",
            ),
            (AtomQueryPredicate::FormalCharge(0), "[+0]"),
        ] {
            assert_eq!(write(&atom(leaf(p))).unwrap(), PropertyText::from(want));
        }
    }
    #[test]
    fn net_leaf_negation_controls_prefix_and_forced_brackets() {
        assert_eq!(
            write(&atom(not(leaf(AtomQueryPredicate::Any)))).unwrap(),
            PropertyText::from("[!*]")
        );
        assert_eq!(
            write(&atom(not(not(leaf(AtomQueryPredicate::Any))))).unwrap(),
            PropertyText::from("*")
        );
        assert_eq!(
            write(&atom(not(not(not(num(6)))))).unwrap(),
            PropertyText::from("[!#6]")
        );
    }
    #[test]
    fn composite_dispatch_uses_own_flag_and_source_child_order() {
        assert_eq!(
            write(&atom(QueryNode::Or(vec![num(6), num(7)]))).unwrap(),
            PropertyText::from("[#6,#7]")
        );
        assert_eq!(
            write(&atom(not(QueryNode::Or(vec![num(6), num(7)])))).unwrap(),
            PropertyText::from("[!#6&!#7]")
        );
        assert_eq!(
            write(&atom(QueryNode::And(vec![num(6), not(num(7))]))).unwrap(),
            PropertyText::from("[#6&!#7]")
        );
    }
    #[test]
    fn bare_recursive_query_wraps_actual_empty_molecule_and_own_negation() {
        assert_eq!(write(&atom(rec())).unwrap(), PropertyText::from("[$()]"));
        assert_eq!(
            write(&atom(not(rec()))).unwrap(),
            PropertyText::from("[!$()]")
        );
        assert_eq!(
            write(&atom(not(not(rec())))).unwrap(),
            PropertyText::from("[$()]")
        );
    }
    #[test]
    fn missing_recursive_molecule_does_not_use_cached_text_or_ordinary_atom() {
        let q = leaf(AtomQueryPredicate::RecursiveSmarts(
            RecursiveStructureQuery::new().with_source_smarts("C"),
        ));
        let mut a = atom(q);
        a.set_prop("molAtomMapNumber", 8).unwrap();
        a.set_prop("smilesSymbol", "X").unwrap();
        assert_eq!(
            write(&a),
            Err(SmartsWriteError::MissingRecursiveQueryMolecule)
        );
    }
    #[test]
    fn raw_map_string_conversion_covers_all_scalar_source_kinds() {
        for (value, want) in [
            (PropertyValue::Int(8), "[#6:8]"),
            (PropertyValue::Int(-8), "[#6:-8]"),
            (PropertyValue::UInt(u32::MAX), "[#6:4294967295]"),
            (PropertyValue::Bool(false), "[#6:0]"),
            (PropertyValue::Double(8.0), "[#6:8]"),
            (PropertyValue::String(" 008 \t".into()), "[#6: 008 \t]"),
            (PropertyValue::String("".into()), "[#6:]"),
        ] {
            let mut a = atom(num(6));
            a.set_prop("molAtomMapNumber", value).unwrap();
            assert_eq!(write(&a).unwrap(), PropertyText::from(want));
        }
    }
    #[test]
    fn map_and_symbol_raw_bytes_are_retained_in_source_output_order() {
        let mut a = atom(leaf(AtomQueryPredicate::Any));
        a.set_prop(
            "molAtomMapNumber",
            PropertyValue::String(PropertyText::from(vec![0xff, 0])),
        )
        .unwrap();
        a.set_prop(
            "smilesSymbol",
            PropertyValue::String(PropertyText::from(vec![0xfe, b'X'])),
        )
        .unwrap();
        assert_eq!(
            write(&a).unwrap().as_bytes(),
            &[b'[', 0xfe, b'X', b';', b'*', b':', 0xff, 0, b']']
        );
    }
    #[test]
    fn authoritative_typed_map_precedes_raw_map_without_signed_conversion() {
        let mut a = atom(num(6));
        a.set_atom_map(Some(u32::MAX));
        a.set_prop("molAtomMapNumber", PropertyValue::String("bad".into()))
            .unwrap();
        assert_eq!(write(&a).unwrap(), PropertyText::from("[#6:4294967295]"));
        a.set_atom_map(None);
        assert_eq!(write(&a).unwrap(), PropertyText::from("[#6:bad]"));
    }
    #[test]
    fn existing_projection_can_hide_both_raw_and_typed_map_properties() {
        let p = SmartsWriteParams {
            include_atom_maps: false,
            ..Default::default()
        };
        let mut a = atom(num(6));
        a.set_prop("molAtomMapNumber", 8).unwrap();
        assert_eq!(
            query_atom_to_smarts(&a, &p).unwrap(),
            PropertyText::from("[#6]")
        );
        a.set_atom_map(Some(9));
        assert_eq!(
            query_atom_to_smarts(&a, &p).unwrap(),
            PropertyText::from("[#6]")
        );
    }
    #[test]
    fn source_symbol_prefix_handles_empty_leaf_result_before_and_after_map() {
        let mut a = atom(num(6));
        a.set_prop("smilesSymbol", "X").unwrap();
        assert_eq!(write(&a).unwrap(), PropertyText::from("[X]"));
        a.set_prop("molAtomMapNumber", 8).unwrap();
        assert_eq!(write(&a).unwrap(), PropertyText::from("[X;:8]"));
        a.set_predicate(not(num(6)));
        assert_eq!(write(&a).unwrap(), PropertyText::from("[X;!:8]"));
        a.set_predicate(leaf(AtomQueryPredicate::AtomType {
            atomic_number: 6,
            aromatic: false,
        }));
        assert_eq!(write(&a).unwrap(), PropertyText::from("[X;A:8]"));
    }
    #[test]
    fn symbol_property_uses_source_string_conversion_and_empty_string_presence() {
        let mut a = atom(num(6));
        a.set_prop("smilesSymbol", PropertyValue::Int(7)).unwrap();
        assert_eq!(write(&a).unwrap(), PropertyText::from("[7]"));
        a.set_prop("smilesSymbol", "").unwrap();
        assert_eq!(write(&a).unwrap(), PropertyText::from("[]"));
        a.set_predicate(leaf(AtomQueryPredicate::Any));
        assert_eq!(write(&a).unwrap(), PropertyText::from("[;*]"));
    }
    #[test]
    fn source_defined_unknown_leaf_fallback_remains_distinct_from_unmodeled_identity() {
        assert_eq!(
            write(&atom(leaf(AtomQueryPredicate::NumRadicalElectrons(2)))).unwrap(),
            PropertyText::from("*")
        );
        assert_eq!(
            write(&atom(not(leaf(AtomQueryPredicate::NumRadicalElectrons(2))))).unwrap(),
            PropertyText::from("[!*]")
        );
        assert!(matches!(
            write(&atom(leaf(AtomQueryPredicate::UnsupportedFeature(
                "unmodeled identity"
            )))),
            Err(SmartsWriteError::UnsupportedAtomQuery { .. })
        ));
    }
    #[test]
    fn independent_xor_and_malformed_composite_remain_explicit_errors() {
        assert_eq!(
            write(&atom(QueryNode::Xor(vec![num(6), num(7)]))),
            Err(SmartsWriteError::XorComposite)
        );
        assert_eq!(
            write(&atom(QueryNode::And(vec![num(6)]))),
            Err(SmartsWriteError::CompositeChildCount { kind: "atom" })
        );
    }
    #[test]
    fn query_leaf_stereo_is_written_once_and_does_not_copy_ordinary_h_rule() {
        let mut a = atom(leaf(AtomQueryPredicate::AtomType {
            atomic_number: 6,
            aromatic: false,
        }));
        a.set_chiral_tag(ChiralTag::TetrahedralCw);
        a.set_explicit_hydrogens(1);
        assert_eq!(write(&a).unwrap(), PropertyText::from("[C@@]"));
        a.set_predicate(QueryNode::And(vec![
            leaf(AtomQueryPredicate::Any),
            leaf(AtomQueryPredicate::Any),
        ]));
        assert_eq!(write(&a).unwrap(), PropertyText::from("[*@@&*]"));
        a.set_prop("_qatomHasStereoSet", false).unwrap();
        assert_eq!(write(&a).unwrap(), PropertyText::from("[*&*]"));
    }
    #[test]
    fn successful_and_failed_writes_leave_all_represented_atom_state_immutable() {
        let mut a = atom(not(QueryNode::Or(vec![num(6), num(7)])));
        a.set_atom_map(Some(17));
        a.set_prop("other", PropertyValue::IntVector(vec![3, 2, 1]))
            .unwrap();
        let before = format!("{a:?}");
        assert_eq!(write(&a).unwrap(), PropertyText::from("[!#6&!#7:17]"));
        assert_eq!(format!("{a:?}"), before);
        a.set_predicate(leaf(AtomQueryPredicate::UnsupportedFeature("unmodeled")));
        let before = format!("{a:?}");
        assert!(write(&a).is_err());
        assert_eq!(format!("{a:?}"), before);
    }
}

#[cfg(test)]
mod complete_basic_bond_repr_source_tests {
    use super::*;
    #[test]
    fn pinned_bond_type_table_including_source_empty_default() {
        let tokens = [
            "", "-", "=", "#", "$", "", "", "", "", "", "", "", ":", "", "", "", "", "->", "", "",
            "", "~",
        ];
        for (code, token) in tokens.into_iter().enumerate() {
            let order = BondOrder::from_rdkit_code(code as i64).unwrap();
            assert_eq!(
                get_basic_bond_repr(order, BondDirection::None, false, &Default::default()),
                token
            );
        }
    }
    #[test]
    fn only_single_and_aromatic_types_use_source_end_direction_tokens() {
        for code in 0..=21 {
            for dir in 0..=6 {
                for iso in [false, true] {
                    let order = BondOrder::from_rdkit_code(code).unwrap();
                    let direction = BondDirection::from_rdkit_code(dir).unwrap();
                    let p = SmartsWriteParams {
                        isomeric_smiles: iso,
                        ..Default::default()
                    };
                    let ordinary = [
                        "", "-", "=", "#", "$", "", "", "", "", "", "", "", ":", "", "", "", "",
                        "->", "", "", "", "~",
                    ][code as usize];
                    let want = if iso && [1, 12].contains(&code) && dir == 3 {
                        "\\"
                    } else if iso && [1, 12].contains(&code) && dir == 4 {
                        "/"
                    } else {
                        ordinary
                    };
                    assert_eq!(
                        get_basic_bond_repr(order, direction, false, &p),
                        want,
                        "type{code} dir{dir} iso{iso}"
                    );
                }
            }
        }
    }
    #[test]
    fn dative_projection_and_left_right_orientation_follow_all_guards() {
        for include in [false, true] {
            for reverse in [false, true] {
                for iso in [false, true] {
                    let p = SmartsWriteParams {
                        include_dative_bonds: include,
                        isomeric_smiles: iso,
                        ..Default::default()
                    };
                    for dir in 0..=6 {
                        assert_eq!(
                            get_basic_bond_repr(
                                BondOrder::Dative,
                                BondDirection::from_rdkit_code(dir).unwrap(),
                                reverse,
                                &p
                            ),
                            if !include {
                                "-"
                            } else if reverse {
                                "<-"
                            } else {
                                "->"
                            }
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn irrelevant_options_and_reverse_flag_do_not_change_other_source_types() {
        for code in 0..=21 {
            if code == 17 {
                continue;
            }
            let order = BondOrder::from_rdkit_code(code).unwrap();
            let p = SmartsWriteParams {
                include_atom_maps: false,
                include_dative_bonds: false,
                rooted_at_atom: Some(usize::MAX),
                ..Default::default()
            };
            let token = get_basic_bond_repr(order, BondDirection::None, false, &Default::default());
            assert_eq!(
                get_basic_bond_repr(order, BondDirection::None, true, &p),
                token
            );
            // Ownership is acquired only at the caller boundary, with no text conversion.
            assert_eq!(PropertyText::from(token).as_bytes(), token.as_bytes());
        }
    }
}

#[cfg(test)]
mod complete_bond_smarts_simple_source_tests {
    use super::*;
    use cosmolkit_model::{BondSpec, PropertyValue};
    fn bond(dir: BondDirection) -> Bond {
        Bond::from_spec(
            BondId::new(2),
            BondSpec::new(AtomId::new(3), AtomId::new(4), BondOrder::Triple).with_direction(dir),
        )
    }
    fn write(
        b: &Bond,
        q: &BondQueryPredicate,
        left: Option<usize>,
    ) -> Result<PropertyText, SmartsWriteError> {
        get_bond_smarts_simple(b, q, left, &Default::default())
    }
    #[test]
    fn literal_null_ring_and_canonical_multi_order_queries_follow_source_tokens() {
        for (q, want) in [
            (BondQueryPredicate::Any, "~"),
            (BondQueryPredicate::IsInRing(true), "@"),
            (BondQueryPredicate::IsInRing(false), "@"),
            (
                BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Double]),
                "-,=",
            ),
            (
                BondQueryPredicate::OrderIn(vec![BondOrder::Double, BondOrder::Aromatic]),
                "=,:",
            ),
            (
                BondQueryPredicate::OrderIn(vec![
                    BondOrder::Single,
                    BondOrder::Double,
                    BondOrder::Aromatic,
                ]),
                "-,=,:",
            ),
        ] {
            assert_eq!(
                write(&bond(BondDirection::None), &q, Some(usize::MAX)).unwrap(),
                PropertyText::from(want)
            );
        }
    }
    #[test]
    fn single_or_aromatic_reads_carrier_direction_without_isomeric_guard() {
        let q = BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Aromatic]);
        for code in 0..=6 {
            for iso in [false, true] {
                for dative in [false, true] {
                    let b = bond(BondDirection::from_rdkit_code(code).unwrap());
                    let p = SmartsWriteParams {
                        isomeric_smiles: iso,
                        include_dative_bonds: dative,
                        ..Default::default()
                    };
                    assert_eq!(
                        get_bond_smarts_simple(&b, &q, Some(usize::MAX), &p).unwrap(),
                        PropertyText::from(match code {
                            3 => "\\",
                            4 => "/",
                            _ => "",
                        })
                    );
                }
            }
        }
    }
    #[test]
    fn direction_query_uses_its_value_and_source_rejects_all_other_directions() {
        for code in 0..=6 {
            let direction = BondDirection::from_rdkit_code(code).unwrap();
            let q = BondQueryPredicate::Direction(direction);
            let p = SmartsWriteParams {
                isomeric_smiles: false,
                ..Default::default()
            };
            let out =
                get_bond_smarts_simple(&bond(BondDirection::EndUpRight), &q, Some(usize::MAX), &p);
            match code {
                3 => assert_eq!(out.unwrap(), PropertyText::from("\\")),
                4 => assert_eq!(out.unwrap(), PropertyText::from("/")),
                _ => assert_eq!(
                    out,
                    Err(SmartsWriteError::SourceBondDirection { direction })
                ),
            }
        }
    }
    #[test]
    fn order_query_uses_query_value_and_source_endpoint_orientation() {
        let b = bond(BondDirection::EndDownRight);
        for (order, left, want) in [
            (BondOrder::Single, None, "\\"),
            (BondOrder::Aromatic, None, "\\"),
            (BondOrder::Double, Some(9), "="),
            (BondOrder::Dative, None, "->"),
            (BondOrder::Dative, Some(3), "->"),
            (BondOrder::Dative, Some(4), "<-"),
            (BondOrder::Unspecified, None, ""),
            (BondOrder::Zero, Some(9), "~"),
        ] {
            assert_eq!(
                write(&b, &BondQueryPredicate::Order(order), left).unwrap(),
                PropertyText::from(want)
            );
        }
        let p = SmartsWriteParams {
            isomeric_smiles: false,
            include_dative_bonds: false,
            ..Default::default()
        };
        assert_eq!(
            get_bond_smarts_simple(&b, &BondQueryPredicate::Order(BondOrder::Single), None, &p)
                .unwrap(),
            PropertyText::from("-")
        );
        assert_eq!(
            get_bond_smarts_simple(
                &b,
                &BondQueryPredicate::Order(BondOrder::Dative),
                Some(4),
                &p
            )
            .unwrap(),
            PropertyText::from("-")
        );
    }
    #[test]
    fn known_source_unwritable_queries_preserve_native_description_errors() {
        for (q, description) in [
            (BondQueryPredicate::HasStereo, "BondStereo"),
            (BondQueryPredicate::NumRingBonds(-1), "BondInNRings"),
            (BondQueryPredicate::InRingOfSize(6), "BondRingSize"),
            (BondQueryPredicate::MinRingSize(3), "BondMinRingSize"),
            (BondQueryPredicate::HasProperty("x".into()), "HasProp"),
            (
                BondQueryPredicate::PropertyValue {
                    name: "x".into(),
                    value: "v".into(),
                },
                "HasPropWithValue",
            ),
        ] {
            let e = write(&bond(BondDirection::None), &q, Some(usize::MAX)).unwrap_err();
            assert_eq!(e, SmartsWriteError::UnwritableBondQuery { description });
            assert_eq!(
                e.to_string(),
                format!("Can't write smarts for this query bond type: {description}")
            );
        }
    }
    #[test]
    fn arbitrary_order_vectors_are_not_guessed_sorted_or_deduplicated() {
        for orders in [
            vec![],
            vec![BondOrder::Single],
            vec![BondOrder::Aromatic, BondOrder::Single],
            vec![BondOrder::Single, BondOrder::Single, BondOrder::Aromatic],
            vec![BondOrder::Single, BondOrder::Triple],
        ] {
            let q = BondQueryPredicate::OrderIn(orders);
            assert_eq!(
                write(&bond(BondDirection::None), &q, None),
                Err(SmartsWriteError::UnsupportedBondQuery {
                    predicate: q.clone()
                })
            );
        }
    }
    #[test]
    fn source_width_checks_follow_actual_order_branch_and_left_short_circuit() {
        let b = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(usize::MAX), AtomId::new(0), BondOrder::Single),
        );
        assert_eq!(
            write(&b, &BondQueryPredicate::Any, Some(usize::MAX)).unwrap(),
            PropertyText::from("~")
        );
        let q = BondQueryPredicate::Order(BondOrder::Dative);
        assert_eq!(write(&b, &q, None).unwrap(), PropertyText::from("->"));
        assert_eq!(
            write(&b, &q, Some(usize::MAX)),
            Err(SmartsWriteError::SourceAtomToLeftIndex { index: usize::MAX })
        );
        assert_eq!(
            write(&b, &q, Some(0)),
            Err(SmartsWriteError::SourceBondBeginIndex { index: usize::MAX })
        );
        assert_eq!(
            write(&bond(BondDirection::None), &q, Some(i32::MAX as usize)).unwrap(),
            PropertyText::from("<-")
        );
        assert!(matches!(
            write(
                &bond(BondDirection::None),
                &BondQueryPredicate::Order(BondOrder::Double),
                Some(i32::MAX as usize + 1)
            ),
            Err(SmartsWriteError::SourceAtomToLeftIndex { .. })
        ));
    }
    #[test]
    fn represented_inputs_and_unreached_properties_are_not_mutated_or_interpreted() {
        let mut b = bond(BondDirection::EndDownRight);
        b.set_prop(
            "unused",
            PropertyValue::String(PropertyText::from(vec![0xff, 0])),
        )
        .unwrap();
        let before = format!("{b:?}");
        let q = BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Aromatic]);
        let before_q = q.clone();
        assert_eq!(write(&b, &q, None).unwrap(), PropertyText::from("\\"));
        assert_eq!(format!("{b:?}"), before);
        assert_eq!(q, before_q);
        assert!(write(&b, &BondQueryPredicate::HasStereo, None).is_err());
        assert_eq!(format!("{b:?}"), before);
    }
}

#[cfg(test)]
mod complete_recurse_bond_smarts_source_tests {
    use super::*;
    use cosmolkit_model::BondSpec;
    fn bond() -> Bond {
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(3), AtomId::new(4), BondOrder::Single),
        )
    }
    fn order(o: BondOrder) -> QueryNode<BondQueryPredicate> {
        QueryNode::Predicate(BondQueryPredicate::Order(o))
    }
    fn not(q: QueryNode<BondQueryPredicate>) -> QueryNode<BondQueryPredicate> {
        QueryNode::Not(Box::new(q))
    }
    fn render(
        q: &QueryNode<BondQueryPredicate>,
        neg: bool,
        f: &mut QueryBoolFeatures,
    ) -> Result<PropertyText, SmartsWriteError> {
        recurse_bond_smarts(&bond(), q, neg, None, f, &Default::default())
    }
    #[test]
    fn source_binary_connectives_and_parent_negation_preserve_feature_words() {
        for or in [false, true] {
            for neg in [false, true] {
                let q = if or {
                    QueryNode::Or(vec![order(BondOrder::Single), order(BondOrder::Double)])
                } else {
                    QueryNode::And(vec![order(BondOrder::Single), order(BondOrder::Double)])
                };
                let mut f = QueryBoolFeatures(0x8000_0000);
                let sep = if or ^ neg { "," } else { "&" };
                let prefix = if neg { "!" } else { "" };
                assert_eq!(
                    render(&q, neg, &mut f).unwrap(),
                    PropertyText::from(format!("{prefix}-{sep}{prefix}="))
                );
                assert_eq!(f.0, 0x8000_0000 | if or ^ neg { 4 } else { 1 });
            }
        }
    }
    #[test]
    fn source_own_negation_and_not_chain_projection_do_not_reorder_children() {
        let q = QueryNode::And(vec![
            not(order(BondOrder::Single)),
            not(not(order(BondOrder::Double))),
        ]);
        assert_eq!(
            render(&q, false, &mut QueryBoolFeatures(0)).unwrap(),
            PropertyText::from("!-&=")
        );
        assert_eq!(
            render(&q, true, &mut QueryBoolFeatures(0)).unwrap(),
            PropertyText::from("-,!=")
        );
        let q = not(QueryNode::Or(vec![
            order(BondOrder::Single),
            order(BondOrder::Double),
        ]));
        assert_eq!(
            render(&q, false, &mut QueryBoolFeatures(0)).unwrap(),
            PropertyText::from("!-&!=")
        );
    }
    #[test]
    fn second_composite_overwrites_first_text_exactly_as_pinned_source() {
        let q = QueryNode::And(vec![
            order(BondOrder::Single),
            QueryNode::Or(vec![order(BondOrder::Double), order(BondOrder::Triple)]),
        ]);
        let mut f = QueryBoolFeatures(0);
        assert_eq!(
            render(&q, false, &mut f).unwrap(),
            PropertyText::from("=,#")
        );
        assert_eq!(f.0, 6);
        assert_eq!(
            render(&q, true, &mut QueryBoolFeatures(0)).unwrap(),
            PropertyText::from("!=&!#")
        );
        let q = QueryNode::And(vec![
            QueryNode::Or(vec![order(BondOrder::Double), order(BondOrder::Triple)]),
            order(BondOrder::Single),
        ]);
        assert_eq!(
            render(&q, false, &mut QueryBoolFeatures(0)).unwrap(),
            PropertyText::from("=,#;-")
        );
    }
    #[test]
    fn both_composites_keep_feature_words_even_when_first_text_is_overwritten() {
        let q = QueryNode::And(vec![
            QueryNode::Or(vec![order(BondOrder::Single), order(BondOrder::Double)]),
            QueryNode::And(vec![order(BondOrder::Triple), order(BondOrder::Zero)]),
        ]);
        let mut f = QueryBoolFeatures(0);
        assert_eq!(
            render(&q, false, &mut f).unwrap(),
            PropertyText::from("#&~")
        );
        assert_eq!(f.0, 7);
        let q = QueryNode::Or(vec![
            QueryNode::And(vec![order(BondOrder::Single), order(BondOrder::Double)]),
            QueryNode::Or(vec![order(BondOrder::Double), order(BondOrder::Triple)]),
        ]);
        let mut f = QueryBoolFeatures(0);
        assert_eq!(
            render(&q, false, &mut f).unwrap(),
            PropertyText::from("=,#")
        );
        assert_eq!(f.0, 5);
    }
    #[test]
    fn child_error_order_and_failed_combine_leave_outer_features_unchanged() {
        let bad = QueryNode::Predicate(BondQueryPredicate::HasStereo);
        let wrong_dir = QueryNode::Predicate(BondQueryPredicate::Direction(BondDirection::None));
        let q = QueryNode::And(vec![bad, wrong_dir]);
        let mut f = QueryBoolFeatures(0x100);
        assert_eq!(
            render(&q, false, &mut f),
            Err(SmartsWriteError::UnwritableBondQuery {
                description: "BondStereo"
            })
        );
        assert_eq!(f.0, 0x100);
        let q = QueryNode::And(vec![
            QueryNode::Or(vec![order(BondOrder::Single), order(BondOrder::Double)]),
            QueryNode::Predicate(BondQueryPredicate::HasStereo),
        ]);
        assert!(render(&q, false, &mut f).is_err());
        assert_eq!(f.0, 0x100);
        let q = QueryNode::Or(vec![
            QueryNode::And(vec![
                QueryNode::Or(vec![order(BondOrder::Single), order(BondOrder::Double)]),
                order(BondOrder::Triple),
            ]),
            order(BondOrder::Zero),
        ]);
        assert_eq!(
            render(&q, false, &mut f),
            Err(SmartsWriteError::OrAboveAndBelowAnd)
        );
        assert_eq!(f.0, 0x100);
    }
    #[test]
    fn atom_to_left_and_all_output_options_reach_actual_child_writers() {
        let q = QueryNode::Or(vec![
            order(BondOrder::Dative),
            QueryNode::Predicate(BondQueryPredicate::Direction(BondDirection::EndUpRight)),
        ]);
        let b = bond();
        assert_eq!(
            recurse_bond_smarts(
                &b,
                &q,
                false,
                Some(4),
                &mut QueryBoolFeatures(0),
                &Default::default()
            )
            .unwrap(),
            PropertyText::from("<-,/")
        );
        let p = SmartsWriteParams {
            include_dative_bonds: false,
            isomeric_smiles: false,
            ..Default::default()
        };
        assert_eq!(
            recurse_bond_smarts(&b, &q, false, Some(4), &mut QueryBoolFeatures(0), &p).unwrap(),
            PropertyText::from("-,/")
        );
        assert_eq!(
            recurse_bond_smarts(&b, &q, true, Some(4), &mut QueryBoolFeatures(0), &p).unwrap(),
            PropertyText::from("!-&!/")
        );
    }
    #[test]
    fn exact_binary_child_count_is_checked_before_any_child_dispatch() {
        for q in [
            QueryNode::And(vec![]),
            QueryNode::Or(vec![order(BondOrder::Single)]),
            QueryNode::And(vec![
                QueryNode::Predicate(BondQueryPredicate::HasStereo),
                order(BondOrder::Double),
                order(BondOrder::Triple),
            ]),
            order(BondOrder::Single),
        ] {
            let mut f = QueryBoolFeatures(0x100);
            assert_eq!(
                recurse_bond_smarts(
                    &bond(),
                    &q,
                    false,
                    Some(usize::MAX),
                    &mut f,
                    &Default::default()
                ),
                Err(SmartsWriteError::CompositeChildCount { kind: "bond" })
            );
            assert_eq!(f.0, 0x100);
        }
    }
    #[test]
    fn tree_and_bond_inputs_are_immutable_including_error_paths_and_empty_tokens() {
        let b = bond();
        let before_b = format!("{b:?}");
        let q = QueryNode::Or(vec![
            order(BondOrder::Unspecified),
            QueryNode::Predicate(BondQueryPredicate::Any),
        ]);
        let before_q = format!("{q:?}");
        let mut f = QueryBoolFeatures(0);
        assert_eq!(
            recurse_bond_smarts(&b, &q, false, None, &mut f, &Default::default()).unwrap(),
            PropertyText::from("~")
        );
        assert_eq!(f.0, 4);
        assert_eq!(format!("{b:?}"), before_b);
        assert_eq!(format!("{q:?}"), before_q);
        let q = QueryNode::And(vec![
            order(BondOrder::Single),
            QueryNode::Predicate(BondQueryPredicate::HasStereo),
        ]);
        let before_q = format!("{q:?}");
        assert!(recurse_bond_smarts(&b, &q, false, None, &mut f, &Default::default()).is_err());
        assert_eq!(format!("{b:?}"), before_b);
        assert_eq!(format!("{q:?}"), before_q);
    }
}

#[cfg(test)]
mod complete_non_query_bond_smarts_source_tests {
    use super::*;
    use cosmolkit_model::{BondSpec, PropertyValue};
    fn bond(order: BondOrder, aromatic: bool, dir: BondDirection) -> Bond {
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(3), AtomId::new(4), order)
                .with_aromatic(aromatic)
                .with_direction(dir),
        )
    }
    #[test]
    fn aromatic_flag_overrides_all_stored_orders_before_endpoint_reads() {
        let tokens = [
            "", "-", "=", "#", "$", "", "", "", "", "", "", "", ":", "", "", "", "", "->", "", "",
            "", "~",
        ];
        for code in 0..=21 {
            let o = BondOrder::from_rdkit_code(code).unwrap();
            assert_eq!(
                non_query_bond_to_smarts(
                    &bond(o, false, BondDirection::None),
                    None,
                    &Default::default()
                )
                .unwrap(),
                tokens[code as usize]
            );
            assert_eq!(
                non_query_bond_to_smarts(
                    &bond(o, true, BondDirection::None),
                    Some(usize::MAX),
                    &Default::default()
                )
                .unwrap(),
                ":"
            );
        }
    }
    #[test]
    fn source_aromatic_directions_require_isomeric_output_regardless_of_stored_order() {
        for code in 0..=21 {
            for d in 0..=6 {
                for iso in [false, true] {
                    let b = bond(
                        BondOrder::from_rdkit_code(code).unwrap(),
                        true,
                        BondDirection::from_rdkit_code(d).unwrap(),
                    );
                    let p = SmartsWriteParams {
                        isomeric_smiles: iso,
                        include_dative_bonds: false,
                        ..Default::default()
                    };
                    let want = if iso && d == 3 {
                        "\\"
                    } else if iso && d == 4 {
                        "/"
                    } else {
                        ":"
                    };
                    assert_eq!(
                        non_query_bond_to_smarts(&b, Some(usize::MAX), &p).unwrap(),
                        want
                    );
                }
            }
        }
    }
    #[test]
    fn nonaromatic_source_dative_orientation_and_projection_are_shared() {
        let b = bond(BondOrder::Dative, false, BondDirection::EndDownRight);
        for include in [false, true] {
            for (left, reversed) in [
                (None, false),
                (Some(3), false),
                (Some(4), true),
                (Some(i32::MAX as usize), true),
            ] {
                let p = SmartsWriteParams {
                    include_dative_bonds: include,
                    ..Default::default()
                };
                assert_eq!(
                    non_query_bond_to_smarts(&b, left, &p).unwrap(),
                    if !include {
                        "-"
                    } else if reversed {
                        "<-"
                    } else {
                        "->"
                    }
                );
            }
        }
        let p = SmartsWriteParams {
            isomeric_smiles: false,
            ..Default::default()
        };
        assert_eq!(
            non_query_bond_to_smarts(
                &bond(BondOrder::Single, false, BondDirection::EndUpRight),
                None,
                &p
            )
            .unwrap(),
            "-"
        );
    }
    #[test]
    fn width_errors_and_missing_left_follow_actual_source_short_circuits() {
        let b = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(usize::MAX), AtomId::new(0), BondOrder::Double),
        );
        assert_eq!(
            non_query_bond_to_smarts(&b, None, &Default::default()).unwrap(),
            "="
        );
        assert_eq!(
            non_query_bond_to_smarts(&b, Some(usize::MAX), &Default::default()),
            Err(SmartsWriteError::SourceAtomToLeftIndex { index: usize::MAX })
        );
        assert_eq!(
            non_query_bond_to_smarts(&b, Some(0), &Default::default()),
            Err(SmartsWriteError::SourceBondBeginIndex { index: usize::MAX })
        );
        let b = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(usize::MAX), AtomId::new(0), BondOrder::Double)
                .with_aromatic(true),
        );
        assert_eq!(
            non_query_bond_to_smarts(&b, Some(usize::MAX), &Default::default()).unwrap(),
            ":"
        );
    }
    #[test]
    fn helper_does_not_dispatch_query_or_owner_properties_and_never_mutates_input() {
        let mut b = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(3), AtomId::new(4), BondOrder::Triple)
                .with_query(QueryNode::Predicate(BondQueryPredicate::Any)),
        );
        b.set_prop("_doIsoSmiles", PropertyValue::Bool(true))
            .unwrap();
        let before = format!("{b:?}");
        assert_eq!(
            non_query_bond_to_smarts(&b, None, &Default::default()).unwrap(),
            "#"
        );
        assert_eq!(format!("{b:?}"), before);
        assert!(non_query_bond_to_smarts(&b, Some(usize::MAX), &Default::default()).is_err());
        assert_eq!(format!("{b:?}"), before);
    }
}

#[cfg(test)]
mod complete_query_bond_smarts_source_tests {
    use super::*;
    use cosmolkit_model::BondSpec;
    fn carrier(order: BondOrder, aromatic: bool, dir: BondDirection) -> Bond {
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(3), AtomId::new(4), order)
                .with_aromatic(aromatic)
                .with_direction(dir),
        )
    }
    fn leaf(o: BondOrder) -> QueryNode<BondQueryPredicate> {
        QueryNode::Predicate(BondQueryPredicate::Order(o))
    }
    fn not(q: QueryNode<BondQueryPredicate>) -> QueryNode<BondQueryPredicate> {
        QueryNode::Not(Box::new(q))
    }
    fn explicit(q: QueryNode<BondQueryPredicate>) -> QueryBond {
        QueryBond::from_parts(carrier(BondOrder::Single, false, BondDirection::None), q)
    }
    fn write(b: &QueryBond) -> Result<PropertyText, SmartsWriteError> {
        query_bond_to_smarts(b, &Default::default(), None)
    }
    #[test]
    fn explicit_order_dispatch_uses_query_value_for_every_source_type() {
        let tokens = [
            "", "-", "=", "#", "$", "", "", "", "", "", "", "", ":", "", "", "", "", "->", "", "",
            "", "~",
        ];
        for code in 0..=21 {
            assert_eq!(
                write(&explicit(leaf(BondOrder::from_rdkit_code(code).unwrap()))).unwrap(),
                PropertyText::from(tokens[code as usize])
            );
        }
        assert_eq!(
            write(&QueryBond::new(
                BondId::new(0),
                BondSpec::new(AtomId::new(3), AtomId::new(4), BondOrder::Unspecified)
            ))
            .unwrap(),
            PropertyText::from("~")
        );
    }
    #[test]
    fn ordinary_origin_returns_before_reading_the_uniform_storage_predicate() {
        let b = QueryBond::from_carrier_parts(
            carrier(BondOrder::Double, true, BondDirection::None),
            QueryNode::Predicate(BondQueryPredicate::HasStereo),
        );
        assert!(b.predicate_is_carrier_derived());
        assert_eq!(write(&b).unwrap(), PropertyText::from(":"));
        let b = QueryBond::from_parts(
            carrier(BondOrder::Double, true, BondDirection::None),
            QueryNode::Predicate(BondQueryPredicate::HasStereo),
        );
        assert_eq!(
            write(&b),
            Err(SmartsWriteError::UnwritableBondQuery {
                description: "BondStereo"
            })
        );
    }
    #[test]
    fn own_negation_cancels_and_keeps_literal_prefix_for_an_empty_leaf() {
        let q = QueryNode::Predicate(BondQueryPredicate::OrderIn(vec![
            BondOrder::Single,
            BondOrder::Aromatic,
        ]));
        assert_eq!(
            write(&explicit(not(q.clone()))).unwrap(),
            PropertyText::from("!")
        );
        assert_eq!(write(&explicit(not(not(q)))).unwrap(), PropertyText::new());
        assert_eq!(
            write(&explicit(not(leaf(BondOrder::Double)))).unwrap(),
            PropertyText::from("!=")
        );
    }
    #[test]
    fn source_composite_flag_and_second_composite_overwrite_reach_public_dispatch() {
        let q = QueryNode::And(vec![
            leaf(BondOrder::Single),
            QueryNode::Or(vec![leaf(BondOrder::Double), leaf(BondOrder::Triple)]),
        ]);
        assert_eq!(
            write(&explicit(q.clone())).unwrap(),
            PropertyText::from("=,#")
        );
        assert_eq!(
            write(&explicit(not(q))).unwrap(),
            PropertyText::from("!=&!#")
        );
    }
    #[test]
    fn actual_source_params_are_forwarded_for_simple_and_ordinary_branches() {
        let b = QueryBond::from_parts(
            carrier(BondOrder::Single, false, BondDirection::EndDownRight),
            leaf(BondOrder::Single),
        );
        let p = SmartsWriteParams {
            isomeric_smiles: false,
            include_dative_bonds: false,
            ..Default::default()
        };
        assert_eq!(
            query_bond_to_smarts(&b, &p, None).unwrap(),
            PropertyText::from("-")
        );
        let b = explicit(leaf(BondOrder::Dative));
        assert_eq!(
            query_bond_to_smarts(&b, &Default::default(), Some(4)).unwrap(),
            PropertyText::from("<-")
        );
        assert_eq!(
            query_bond_to_smarts(&b, &p, Some(4)).unwrap(),
            PropertyText::from("-")
        );
        let b = QueryBond::from_carrier_parts(
            carrier(BondOrder::Single, true, BondDirection::EndDownRight),
            leaf(BondOrder::Double),
        );
        assert_eq!(
            query_bond_to_smarts(&b, &p, None).unwrap(),
            PropertyText::from(":")
        );
    }
    #[test]
    fn known_source_writer_errors_are_not_hidden_by_carrier_or_generic_fallback() {
        assert_eq!(
            write(&explicit(QueryNode::Predicate(
                BondQueryPredicate::HasStereo
            ))),
            Err(SmartsWriteError::UnwritableBondQuery {
                description: "BondStereo"
            })
        );
        assert_eq!(
            write(&explicit(QueryNode::Predicate(
                BondQueryPredicate::Direction(BondDirection::None)
            ))),
            Err(SmartsWriteError::SourceBondDirection {
                direction: BondDirection::None
            })
        );
        assert_eq!(
            write(&explicit(QueryNode::Xor(vec![
                leaf(BondOrder::Single),
                leaf(BondOrder::Double)
            ]))),
            Err(SmartsWriteError::XorComposite)
        );
    }
    #[test]
    fn origin_and_source_aromatic_branch_determine_whether_left_index_is_read() {
        let b = QueryBond::from_carrier_parts(
            carrier(BondOrder::Double, true, BondDirection::None),
            leaf(BondOrder::Dative),
        );
        assert_eq!(
            query_bond_to_smarts(&b, &Default::default(), Some(usize::MAX)).unwrap(),
            PropertyText::from(":")
        );
        let b = QueryBond::from_parts(
            carrier(BondOrder::Double, true, BondDirection::None),
            leaf(BondOrder::Double),
        );
        assert!(matches!(
            query_bond_to_smarts(&b, &Default::default(), Some(usize::MAX)),
            Err(SmartsWriteError::SourceAtomToLeftIndex { .. })
        ));
        let b = explicit(QueryNode::Predicate(BondQueryPredicate::Any));
        assert_eq!(
            query_bond_to_smarts(&b, &Default::default(), Some(usize::MAX)).unwrap(),
            PropertyText::from("~")
        );
    }
    #[test]
    fn malformed_source_composites_fail_before_any_leaf_or_endpoint_read() {
        for nodes in [
            vec![],
            vec![leaf(BondOrder::Single)],
            vec![
                leaf(BondOrder::Single),
                leaf(BondOrder::Double),
                leaf(BondOrder::Triple),
            ],
        ] {
            assert_eq!(
                query_bond_to_smarts(
                    &explicit(QueryNode::And(nodes)),
                    &Default::default(),
                    Some(usize::MAX)
                ),
                Err(SmartsWriteError::CompositeChildCount { kind: "bond" })
            );
        }
    }
    #[test]
    fn canonical_carrier_constructor_preserves_an_actual_incoming_query() {
        let b = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(3), AtomId::new(4), BondOrder::Triple)
                .with_query(QueryNode::Predicate(BondQueryPredicate::Any)),
        );
        let b =
            QueryBond::from_carrier_parts(b, QueryNode::Predicate(BondQueryPredicate::HasStereo));
        assert!(!b.predicate_is_carrier_derived());
        assert_eq!(write(&b).unwrap(), PropertyText::from("~"));
    }
    #[test]
    fn successful_and_failed_dispatch_do_not_mutate_origin_payload_or_carrier() {
        let b = explicit(not(QueryNode::Or(vec![
            leaf(BondOrder::Single),
            leaf(BondOrder::Double),
        ])));
        let before = format!("{b:?}");
        assert_eq!(write(&b).unwrap(), PropertyText::from("!-&!="));
        assert_eq!(format!("{b:?}"), before);
        let b = explicit(QueryNode::Predicate(BondQueryPredicate::HasStereo));
        let before = format!("{b:?}");
        assert!(write(&b).is_err());
        assert_eq!(format!("{b:?}"), before);
    }
}

#[cfg(test)]
mod ordinary_atom_dispatch_source_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, PropertyValue};
    #[test]
    fn carrier_derived_atom_uses_actual_attributes_and_ignores_synthetic_query() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_isotope(13)
                .with_chiral_tag(ChiralTag::TetrahedralCcw)
                .with_explicit_hydrogens(1)
                .with_formal_charge(2)
                .with_atom_map(17),
        );
        let a = QueryAtom::from_carrier_parts(
            a,
            QueryNode::Predicate(AtomQueryPredicate::UnsupportedFeature(
                "not an actual query",
            )),
        );
        assert_eq!(
            query_atom_to_smarts(&a, &Default::default()).unwrap(),
            PropertyText::from("[13#6@H+2:17]")
        );
    }
    #[test]
    fn ordinary_signed_map_read_differs_from_explicit_query_string_read() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("molAtomMapNumber", "008 \t")
                .unwrap(),
        );
        let a = QueryAtom::from_carrier_parts(a, QueryNode::Predicate(AtomQueryPredicate::Any));
        assert_eq!(
            query_atom_to_smarts(&a, &Default::default()).unwrap(),
            PropertyText::from("[#6:8]")
        );
        let mut a = a;
        a.set_predicate(QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6)));
        assert_eq!(
            query_atom_to_smarts(&a, &Default::default()).unwrap(),
            PropertyText::from("[#6:008 \t]")
        );
    }
    #[test]
    fn ordinary_early_return_skips_fake_recursion_and_retains_source_map_behavior() {
        let a = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C).with_atom_map(9));
        let a = QueryAtom::from_carrier_parts(
            a,
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(
                RecursiveStructureQuery::new(),
            )),
        );
        let p = SmartsWriteParams {
            include_atom_maps: false,
            ..Default::default()
        };
        // Native ordinary early return does not reach query-only map policy.
        assert_eq!(
            query_atom_to_smarts(&a, &p).unwrap(),
            PropertyText::from("[#6:9]")
        );
    }
    #[test]
    fn borrowed_ordinary_transport_retains_raw_bytes_properties_and_error_immutability() {
        let a = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop(
                    "smilesSymbol",
                    PropertyValue::String(PropertyText::from(vec![0xff, 0])),
                )
                .unwrap()
                .with_prop("molAtomMapNumber", PropertyValue::Int(7))
                .unwrap()
                .with_prop(
                    "__computedProps",
                    PropertyValue::StringVector(vec!["other".into()]),
                )
                .unwrap(),
        );
        let mut a = QueryAtom::from_carrier_parts(a, QueryNode::Predicate(AtomQueryPredicate::Any));
        let before = format!("{a:?}");
        assert_eq!(
            query_atom_to_smarts(&a, &Default::default())
                .unwrap()
                .as_bytes(),
            &[b'[', 0xff, 0, b':', b'7', b']']
        );
        assert_eq!(format!("{a:?}"), before);
        a.set_prop("molAtomMapNumber", PropertyValue::Bool(false))
            .unwrap();
        let before = format!("{a:?}");
        assert!(matches!(
            query_atom_to_smarts(&a, &Default::default()),
            Err(SmartsWriteError::AtomMapInt { .. })
        ));
        assert_eq!(format!("{a:?}"), before);
    }
}

#[allow(clippy::too_many_arguments)]
#[cfg(feature = "smiles-integration")]
fn fragment_smarts_construct(
    query: &mut QueryGraph,
    rings: &mut cosmolkit_core::RingInfo,
    atom: usize,
    colors: &mut [cosmolkit_smiles::AtomColor],
    ranks: &[u32],
    params: &SmartsWriteParams,
    atom_ordering: &mut Vec<AtomId>,
    bond_ordering: &mut Vec<BondId>,
    atoms_in_play: &[bool],
    bonds_in_play: Option<&[bool]>,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string FragmentSmartsConstruct(
    // RDKit❗❌:     ROMol &mol, unsigned int atomIdx, std::vector<Canon::AtomColors> &colors,
    // RDKit❗❌:     UINT_VECT &ranks, const SmilesWriteParams &params,
    // RDKit❗❌:     std::vector<unsigned int> &atomOrdering,
    // RDKit❗❌:     std::vector<unsigned int> &bondOrdering,
    // RDKit❗❌:     const boost::dynamic_bitset<> &atomsInPlay,
    // RDKit❗❌:     const boost::dynamic_bitset<> *bondsInPlay) {
    // RDKit❗❌:   // this is dirty trick get around the fact that canonicalizeFragment
    // RDKit❗❌:   // thinks we already called findSSSR - to do some atom ranking
    // RDKit❗❌:   // but for smarts we are going to ignore that part. We will artificially
    // RDKit❗❌:   // set the "SSSR" property to an empty property
    // RDKit❗❌:
    // RDKit❗❌:   mol.getRingInfo()->reset();
    // RDKit❗❌:   mol.getRingInfo()->initialize(FIND_RING_TYPE_SYMM_SSSR);
    // RDKit❗❌:   for (auto &atom : mol.atoms()) {
    // RDKit❗❌:     atom->updatePropertyCache(false);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // For Smarts, we avoid reordering of chiral atoms in canonicalizeFragment.
    // RDKit❗❌:   bool doRandom = false;
    // RDKit❗❌:   bool doChiralInversions = true;
    // RDKit❗❌:   Canon::MolStack molStack;
    // RDKit❗❌:   molStack.reserve(mol.getNumAtoms() + mol.getNumBonds());
    // RDKit❗❌:   Canon::canonicalizeFragment(
    // RDKit❗❌:       mol, atomIdx, colors, ranks, molStack, &atomsInPlay, bondsInPlay, nullptr,
    // RDKit❗❌:       params.doIsomericSmiles, doRandom, doChiralInversions);
    // RDKit❗❌:
    // RDKit❗❌:   // now clear the "SSSR" property
    // RDKit❗❌:   mol.getRingInfo()->reset();
    // RDKit❗❌:
    // RDKit❗❌:   std::stringstream res;
    // RDKit❗❌:   for (auto &msCI : molStack) {
    // RDKit❗❌:     switch (msCI.type) {
    // RDKit❗❌:       case Canon::MOL_STACK_ATOM: {
    // RDKit❗❌:         auto *atm = msCI.obj.atom;
    // RDKit❗❌:         res << SmartsWrite::GetAtomSmarts(atm, params);
    // RDKit❗❌:         atomOrdering.push_back(atm->getIdx());
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:       case Canon::MOL_STACK_BOND: {
    // RDKit❗❌:         auto *bnd = msCI.obj.bond;
    // RDKit❗❌:         res << SmartsWrite::GetBondSmarts(bnd, params, msCI.number);
    // RDKit❗❌:         bondOrdering.push_back(bnd->getIdx());
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:       case Canon::MOL_STACK_RING: {
    // RDKit❗❌:         if (msCI.number < 10) {
    // RDKit❗❌:           res << msCI.number;
    // RDKit❗❌:         } else {
    // RDKit❗❌:           res << "%" << msCI.number;
    // RDKit❗❌:         }
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:       case Canon::MOL_STACK_BRANCH_OPEN: {
    // RDKit❗❌:         res << "(";
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:       case Canon::MOL_STACK_BRANCH_CLOSE: {
    // RDKit❗❌:         res << ")";
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:       default:
    // RDKit❗❌:         break;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res.str();
    // RDKit❗❌: }
    // Source ring reset/empty SymmSSSR, all-carrier cache updates, sole full
    // query-aware Canon pipeline, then reset BEFORE emission. No custom DFS,
    // pre-emission rank heuristic, inferred query bond type, or fallback.
    // Cost ❌: known carrier/property transport copies in shared Canon adapter,
    // Vec-backed short strings and packed-bitset projection costs are explicit.
    rings.reset();
    rings.initialize(cosmolkit_core::RingFindType::SymmSssr);
    let bonds = query
        .bonds()
        .iter()
        .map(|b| b.bond().clone())
        .collect::<Vec<_>>();
    let adjacency = cosmolkit_model::AdjacencyList::from_topology(query.num_atoms(), &bonds);
    let complex_bonds = query
        .bonds()
        .iter()
        .map(cosmolkit_core::query_bond_has_complex_type_query)
        .collect::<Vec<_>>();
    for i in 0..query.num_atoms() {
        cosmolkit_core::update_query_atom_property_cache_source(
            query.atoms_mut(),
            &bonds,
            &adjacency,
            AtomId::new(i),
            false,
            Some(&complex_bonds),
        )?;
    }
    let mut stack = Vec::with_capacity(query.num_atoms() + query.num_bonds());
    cosmolkit_smiles::canonicalize_query_fragment_source(
        query,
        rings,
        atom,
        colors,
        ranks,
        &mut stack,
        Some(atoms_in_play),
        bonds_in_play,
        params.isomeric_smiles,
    )?;
    rings.reset();
    let mut result = PropertyText::new();
    for item in stack {
        match item {
            cosmolkit_smiles::MolStackElem::Atom(i) => {
                let atom = query
                    .atom(i)
                    .ok_or(SmartsWriteError::FragmentAtomOutOfRange { atom: i })?;
                let mut stereo_written = false;
                let written = query_atom_to_smarts_with_state(
                    atom,
                    params,
                    &mut stereo_written,
                    query.prop("_doIsoSmiles").is_some(),
                );
                // Native sets the temporary property inside the leaf writer,
                // before any later map/symbol conversion may fail. Return that
                // already completed effect even when the writer returns Err.
                if stereo_written {
                    query
                        .atom_mut(i)
                        .unwrap()
                        .set_prop("_qatomHasStereoSet", cosmolkit_model::PropertyValue::Int(1))
                        .map_err(|source| SmartsWriteError::AtomPropertyWrite {
                            atom: AtomId::new(i),
                            source,
                        })?;
                }
                result.extend_bytes(written?.as_bytes());
                atom_ordering.push(AtomId::new(i));
            }
            cosmolkit_smiles::MolStackElem::Bond { bond, atom_to_left } => {
                let value = query
                    .bond(bond.index())
                    .ok_or(SmartsWriteError::FragmentBondOutOfRange { bond: bond.index() })?;
                result.extend_bytes(
                    query_bond_to_smarts(value, params, Some(atom_to_left))?.as_bytes(),
                );
                bond_ordering.push(bond);
            }
            cosmolkit_smiles::MolStackElem::Ring(n) => {
                if n >= 10 {
                    result.push_byte(b'%');
                }
                result.extend_bytes(n.to_string().as_bytes());
            }
            cosmolkit_smiles::MolStackElem::BranchOpen(_) => result.push_byte(b'('),
            cosmolkit_smiles::MolStackElem::BranchClose(_) => result.push_byte(b')'),
        }
    }
    Ok(result)
}

#[cfg(all(test, feature = "smiles-integration"))]
mod fragment_smarts_construct_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, QueryAtomIdentity, SourceAtomValenceFacts};
    use cosmolkit_smiles::AtomColor;
    use cosmolkit_types::Element;
    fn graph(nums: &[u8], edges: &[(usize, usize)]) -> QueryGraph {
        let atoms = nums
            .iter()
            .enumerate()
            .map(|(i, &n)| {
                QueryAtom::from_identity_parts(
                    AtomId::new(i),
                    QueryAtomIdentity::from_atomic_number(n),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(n)),
                )
            })
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                QueryBond::new(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect();
        QueryGraph::from_parts(atoms, bonds, [], vec![], vec![], vec![]).unwrap()
    }
    fn rings() -> cosmolkit_core::RingInfo {
        cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 0, 0)
    }
    fn run(
        q: &mut QueryGraph,
        start: usize,
        mask: &[bool],
        bonds: Option<&[bool]>,
    ) -> Result<(PropertyText, Vec<AtomId>, Vec<BondId>), SmartsWriteError> {
        let mut colors = vec![AtomColor::White; q.num_atoms()];
        for (i, &included) in mask.iter().enumerate() {
            if !included {
                colors[i] = AtomColor::Black;
            }
        }
        let ranks = (0..q.num_atoms()).map(|i| i as u32).collect::<Vec<_>>();
        let (mut a, mut b) = (vec![], vec![]);
        fragment_smarts_construct(
            q,
            &mut rings(),
            start,
            &mut colors,
            &ranks,
            &SmartsWriteParams::default(),
            &mut a,
            &mut b,
            mask,
            bonds,
        )
        .map(|s| (s, a, b))
    }
    #[test]
    fn chain_emits_source_stack_including_left_endpoints_and_exact_orders() {
        let mut q = graph(&[6, 8, 7], &[(0, 1), (1, 2)]);
        let (s, a, b) = run(&mut q, 0, &[true; 3], None).unwrap();
        assert_eq!(s, PropertyText::from("[#6]-[#8]-[#7]"));
        assert_eq!(a, vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]);
        assert_eq!(b, vec![BondId::new(0), BondId::new(1)]);
        assert_eq!(
            q.atom(0).unwrap().prop("_TraversalStartPoint"),
            Some(&cosmolkit_model::PropertyValue::Bool(true))
        );
        assert!(q.prop("_StereochemDone").is_some());
    }
    #[test]
    fn ranks_drive_branch_order_through_shared_canon_not_custom_neighbor_sort() {
        let mut q = graph(&[6, 8, 7, 9], &[(0, 1), (0, 2), (0, 3)]);
        let n = q.num_atoms();
        let (mut a, mut b) = (vec![], vec![]);
        let s = fragment_smarts_construct(
            &mut q,
            &mut rings(),
            0,
            &mut vec![AtomColor::White; n],
            &[0, 3, 2, 1],
            &SmartsWriteParams::default(),
            &mut a,
            &mut b,
            &[true; 4],
            None,
        )
        .unwrap();
        assert_eq!(
            a,
            vec![
                AtomId::new(0),
                AtomId::new(3),
                AtomId::new(2),
                AtomId::new(1)
            ]
        );
        assert_eq!(s, PropertyText::from("[#6](-[#9])(-[#7])-[#8]"));
        assert_eq!(b, vec![BondId::new(2), BondId::new(1), BondId::new(0)]);
    }
    #[test]
    fn ring_stack_records_only_emitted_bonds_and_resets_fake_sssr() {
        let mut q = graph(&[6; 3], &[(0, 1), (1, 2), (2, 0)]);
        let mut r = rings();
        let (mut a, mut b) = (vec![], vec![]);
        let s = fragment_smarts_construct(
            &mut q,
            &mut r,
            0,
            &mut [AtomColor::White; 3],
            &[0, 1, 2],
            &SmartsWriteParams::default(),
            &mut a,
            &mut b,
            &[true; 3],
            None,
        )
        .unwrap();
        assert_eq!(s, PropertyText::from("[#6]1-[#6]-[#6]-1"));
        assert_eq!(a.len(), 3);
        assert_eq!(b.len(), 3);
        assert!(!r.is_initialized());
        assert_eq!(r.num_rings(), 0);
    }
    #[test]
    fn property_cache_updates_excluded_atoms_before_any_traversal() {
        let mut q = graph(&[6, 6], &[]);
        let (s, a, b) = run(&mut q, 0, &[true, false], None).unwrap();
        assert_eq!(s, PropertyText::from("[#6]"));
        assert_eq!(a, vec![AtomId::new(0)]);
        assert!(b.is_empty());
        assert_eq!(
            q.atom(1).unwrap().source_valence_facts(),
            SourceAtomValenceFacts {
                explicit_valence: 0,
                implicit_valence: 4
            }
        );
        assert!(q.atom(1).unwrap().prop("_TraversalStartPoint").is_none());
    }
    #[test]
    fn later_invalid_identity_retains_earlier_cache_and_initialized_fake_ring_state() {
        let mut q = graph(&[6, 255], &[]);
        let mut r = rings();
        let (mut a, mut b) = (vec![AtomId::new(99)], vec![]);
        let mut colors = [AtomColor::White; 2];
        assert!(matches!(
            fragment_smarts_construct(
                &mut q,
                &mut r,
                0,
                &mut colors,
                &[0, 1],
                &SmartsWriteParams::default(),
                &mut a,
                &mut b,
                &[true, false],
                None
            ),
            Err(SmartsWriteError::Valence(_))
        ));
        assert_eq!(
            q.atom(0).unwrap().source_valence_facts(),
            SourceAtomValenceFacts {
                explicit_valence: 0,
                implicit_valence: 4
            }
        );
        assert_eq!(
            q.atom(1).unwrap().source_valence_facts(),
            SourceAtomValenceFacts::UNINITIALIZED
        );
        assert!(r.is_symm_sssr());
        assert_eq!(a, vec![AtomId::new(99)]);
        assert_eq!(colors, [AtomColor::White; 2]);
    }
    #[test]
    fn canon_precondition_failure_occurs_after_cache_updates_and_before_reset() {
        let mut q = graph(&[6], &[]);
        let mut r = rings();
        let (mut a, mut b) = (vec![], vec![]);
        assert!(matches!(
            fragment_smarts_construct(
                &mut q,
                &mut r,
                0,
                &mut [],
                &[0],
                &SmartsWriteParams::default(),
                &mut a,
                &mut b,
                &[true],
                None
            ),
            Err(SmartsWriteError::CanonicalTraversal(_))
        ));
        assert_eq!(
            q.atom(0).unwrap().source_valence_facts().implicit_valence,
            4
        );
        assert!(r.is_symm_sssr());
        assert!(a.is_empty());
    }
    #[test]
    fn atom_writer_failure_happens_after_ring_reset_and_keeps_only_emitted_order_prefix() {
        let mut q = graph(&[6, 8], &[(0, 1)]);
        q.atom_mut(1).unwrap().set_predicate(QueryNode::predicate(
            AtomQueryPredicate::HydrogenCount(i32::MIN),
        ));
        // This leaf writes a valid source signed count; use independent XOR error
        // to ensure failure in GetAtomSmarts rather than preparation/classification.
        q.atom_mut(1).unwrap().set_predicate(QueryNode::xor(vec![
            QueryNode::predicate(AtomQueryPredicate::Any),
            QueryNode::predicate(AtomQueryPredicate::Any),
        ]));
        let mut r = rings();
        let (mut a, mut b) = (vec![], vec![]);
        let mut colors = [AtomColor::White; 2];
        assert!(matches!(
            fragment_smarts_construct(
                &mut q,
                &mut r,
                0,
                &mut colors,
                &[0, 1],
                &SmartsWriteParams::default(),
                &mut a,
                &mut b,
                &[true; 2],
                None
            ),
            Err(SmartsWriteError::XorComposite)
        ));
        assert!(!r.is_initialized());
        assert_eq!(a, vec![AtomId::new(0)]);
        assert_eq!(b, vec![BondId::new(0)]);
        assert_eq!(colors, [AtomColor::Black; 2]);
    }
    #[test]
    fn bond_writer_failure_does_not_append_the_failing_bond_order() {
        let mut q = graph(&[6, 8], &[(0, 1)]);
        q.bonds_mut()[0].set_predicate(QueryNode::predicate(BondQueryPredicate::HasStereo));
        let mut r = rings();
        let (mut a, mut b) = (vec![], vec![]);
        assert!(matches!(
            fragment_smarts_construct(
                &mut q,
                &mut r,
                0,
                &mut [AtomColor::White; 2],
                &[0, 1],
                &SmartsWriteParams::default(),
                &mut a,
                &mut b,
                &[true; 2],
                None
            ),
            Err(SmartsWriteError::UnwritableBondQuery {
                description: "BondStereo"
            })
        ));
        assert_eq!(a, vec![AtomId::new(0)]);
        assert!(b.is_empty());
        assert!(!r.is_initialized());
    }
    #[test]
    fn stale_direction_is_cleared_by_canon_before_single_bond_rendering() {
        let mut q = graph(&[6, 8], &[(0, 1)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_direction(BondDirection::EndUpRight);
        let (s, _, _) = run(&mut q, 0, &[true; 2], None).unwrap();
        assert_eq!(s, PropertyText::from("[#6]-[#8]"));
        assert_eq!(q.bond(0).unwrap().bond().direction(), BondDirection::None);
    }
    #[test]
    fn complex_query_type_controls_cache_without_copying_query_into_bond_carrier() {
        let mut q = graph(&[6, 6], &[(0, 1)]);
        q.bonds_mut()[0].set_predicate(QueryNode::or(vec![
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        ]));
        let predicate = q.bond(0).unwrap().predicate().clone();
        let _ = run(&mut q, 0, &[true; 2], None).unwrap();
        assert_eq!(
            q.atom(0).unwrap().source_valence_facts(),
            SourceAtomValenceFacts {
                explicit_valence: 1,
                implicit_valence: 0
            }
        );
        assert!(q.bond(0).unwrap().bond().query().is_none());
        assert_eq!(q.bond(0).unwrap().predicate(), &predicate);
    }
    #[test]
    fn ordinary_carrier_origin_and_raw_member_effects_survive_canon_projection() {
        let atom = cosmolkit_model::Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
        let qatom =
            QueryAtom::from_carrier_parts(atom, QueryNode::predicate(AtomQueryPredicate::Any));
        let mut q =
            QueryGraph::from_parts(vec![qatom], vec![], [], vec![], vec![], vec![]).unwrap();
        q.set_prop(
            "_StereochemDone",
            cosmolkit_model::PropertyValue::Bool(false),
        )
        .unwrap();
        let (s, _, _) = run(&mut q, 0, &[true], None).unwrap();
        assert_eq!(s, PropertyText::from("[#6]"));
        assert!(q.atom(0).unwrap().predicate_is_carrier_derived());
        assert_eq!(
            q.prop("_StereochemDone"),
            Some(&cosmolkit_model::PropertyValue::Bool(false))
        );
    }
    #[test]
    fn excluded_bond_breaks_chirality_in_shared_canonical_pipeline() {
        let mut q = graph(&[6, 9, 17, 35], &[(0, 1), (0, 2), (0, 3)]);
        q.atom_mut(0)
            .unwrap()
            .set_chiral_tag(ChiralTag::TetrahedralCw);
        q.atom_mut(0).unwrap().set_explicit_hydrogens(1);
        q.set_prop(
            "_StereochemDone",
            cosmolkit_model::PropertyValue::Bool(false),
        )
        .unwrap();
        let _ = run(
            &mut q,
            0,
            &[true, true, true, false],
            Some(&[true, true, false]),
        )
        .unwrap();
        assert_eq!(
            q.atom(0).unwrap().prop("_brokenChirality"),
            Some(&cosmolkit_model::PropertyValue::Bool(true))
        );
        assert_eq!(q.atom(0).unwrap().chiral_tag(), ChiralTag::TetrahedralCw);
    }
    #[test]
    fn query_h_count_changes_source_ring_inversion_without_changing_carrier_hydrogens() {
        fn prepared(n: i32) -> QueryGraph {
            let mut q = graph(&[16, 6, 6, 9], &[(0, 1), (1, 2), (2, 0), (0, 3)]);
            let center = q.atom_mut(0).unwrap();
            center.set_chiral_tag(ChiralTag::TetrahedralCw);
            center.set_no_implicit(true);
            center.set_formal_charge(1);
            center.set_predicate(QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(16)),
                QueryNode::predicate(AtomQueryPredicate::HydrogenCount(n)),
            ]));
            q.set_prop(
                "_StereochemDone",
                cosmolkit_model::PropertyValue::Bool(false),
            )
            .unwrap();
            q
        }
        let (mut h0, mut h1) = (prepared(0), prepared(1));
        let before = h1.atom(0).unwrap().predicate().clone();
        let _ = run(&mut h0, 0, &[true; 4], None).unwrap();
        let _ = run(&mut h1, 0, &[true; 4], None).unwrap();
        assert_ne!(
            h0.atom(0).unwrap().chiral_tag(),
            h1.atom(0).unwrap().chiral_tag()
        );
        assert_eq!(h1.atom(0).unwrap().explicit_hydrogens(), 0);
        assert_eq!(
            h1.atom(0).unwrap().source_valence_facts().implicit_valence,
            0
        );
        assert_eq!(h1.atom(0).unwrap().predicate(), &before);
    }
    #[test]
    fn source_atom_stereo_marker_returns_to_working_graph_and_suppresses_second_emission() {
        let mut q = graph(&[6], &[]);
        q.atom_mut(0)
            .unwrap()
            .set_chiral_tag(ChiralTag::TetrahedralCw);
        q.set_prop(
            "_StereochemDone",
            cosmolkit_model::PropertyValue::Bool(false),
        )
        .unwrap();
        q.set_prop("_doIsoSmiles", cosmolkit_model::PropertyValue::Int(1))
            .unwrap();
        let first = run(&mut q, 0, &[true], None).unwrap().0;
        assert_eq!(first, PropertyText::from("[#6@@]"));
        assert_eq!(
            q.atom(0).unwrap().prop("_qatomHasStereoSet"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
        let second = run(&mut q, 0, &[true], None).unwrap().0;
        assert_eq!(second, PropertyText::from("[#6]"));
    }
    #[test]
    fn failed_map_conversion_retains_source_stereo_marker_without_atom_order_append() {
        let atom = cosmolkit_model::Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
        let mut atom =
            QueryAtom::from_carrier_parts(atom, QueryNode::predicate(AtomQueryPredicate::Any));
        atom.set_chiral_tag(ChiralTag::TetrahedralCw);
        atom.set_prop(
            "molAtomMapNumber",
            cosmolkit_model::PropertyValue::Double(0.5),
        )
        .unwrap();
        let mut q = QueryGraph::from_parts(vec![atom], vec![], [], vec![], vec![], vec![]).unwrap();
        q.set_prop(
            "_StereochemDone",
            cosmolkit_model::PropertyValue::Bool(false),
        )
        .unwrap();
        q.set_prop("_doIsoSmiles", cosmolkit_model::PropertyValue::Int(1))
            .unwrap();
        let (mut a, mut b) = (vec![], vec![]);
        let mut r = rings();
        assert!(matches!(
            fragment_smarts_construct(
                &mut q,
                &mut r,
                0,
                &mut [AtomColor::White],
                &[0],
                &SmartsWriteParams::default(),
                &mut a,
                &mut b,
                &[true],
                None
            ),
            Err(SmartsWriteError::AtomMapInt { .. })
        ));
        assert_eq!(
            q.atom(0).unwrap().prop("_qatomHasStereoSet"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
        assert!(a.is_empty());
        assert!(!r.is_initialized());
    }

    #[test]
    fn source_atom_owner_flag_is_presence_only_and_independent_of_params() {
        let mut q = graph(&[6], &[]);
        q.atom_mut(0)
            .unwrap()
            .set_chiral_tag(ChiralTag::TetrahedralCw);
        q.set_prop(
            "_StereochemDone",
            cosmolkit_model::PropertyValue::Bool(false),
        )
        .unwrap();
        assert_eq!(
            run(&mut q, 0, &[true], None).unwrap().0,
            PropertyText::from("[#6]")
        );
        assert!(q.atom(0).unwrap().prop("_qatomHasStereoSet").is_none());
        q.set_prop("_doIsoSmiles", cosmolkit_model::PropertyValue::Bool(false))
            .unwrap();
        let (mut a, mut b) = (vec![], vec![]);
        let mut params = SmartsWriteParams::default();
        params.isomeric_smiles = false;
        assert_eq!(
            fragment_smarts_construct(
                &mut q,
                &mut rings(),
                0,
                &mut [AtomColor::White],
                &[0],
                &params,
                &mut a,
                &mut b,
                &[true],
                None
            )
            .unwrap(),
            PropertyText::from("[#6@@]")
        );
    }
}

#[cfg(feature = "smiles-integration")]
fn mol_to_smarts_source(
    query: &QueryGraph,
    params: &SmartsWriteParams,
    mut colors: Vec<cosmolkit_smiles::AtomColor>,
    atoms_in_play: &[bool],
    bonds_in_play: Option<&[bool]>,
) -> Result<SmartsWriteResult, SmartsWriteError> {
    // RDKit❗❌: std::string molToSmarts(const ROMol &inmol, const SmilesWriteParams &params,
    // RDKit❗❌:                         std::vector<AtomColors> &&colors,
    // RDKit❗❌:                         const boost::dynamic_bitset<> &atomsInPlay,
    // RDKit❗❌:                         const boost::dynamic_bitset<> *bondsInPlay) {
    // RDKit❗❌:   PRECONDITION(params.rootedAtAtom < static_cast<int>(inmol.getNumAtoms()),
    // RDKit❗❌:                "bad atom index");
    // RDKit❗❌:   ROMol mol(inmol);
    // RDKit❗❌:   const unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗❌:   UINT_VECT ranks;
    // RDKit❗❌:   ranks.reserve(nAtoms);
    // RDKit❗❌:   // For smiles writing we would be canonicalizing but we will not do that
    // RDKit❗❌:   // here. We will simply use the atom indices as the rank
    // RDKit❗❌:   for (const auto &atom : mol.atoms()) {
    // RDKit❗❌:     ranks.push_back(atom->getIdx());
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (params.doIsomericSmiles) {
    // RDKit❗❌:     mol.setProp(common_properties::_doIsoSmiles, 1);
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<unsigned int> atomOrdering;
    // RDKit❗❌:   std::vector<unsigned int> bondOrdering;
    // RDKit❗❌:
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:   auto colorIt = std::find(colors.begin(), colors.end(), Canon::WHITE_NODE);
    // RDKit❗❌:   while (colorIt != colors.end()) {
    // RDKit❗❌:     unsigned int nextAtomIdx = 0;
    // RDKit❗❌:     std::string subSmi;
    // RDKit❗❌:
    // RDKit❗❌:     if (params.rootedAtAtom > -1 &&
    // RDKit❗❌:         colors[params.rootedAtAtom] == Canon::WHITE_NODE) {
    // RDKit❗❌:       nextAtomIdx = params.rootedAtAtom;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // Try to find a non-chiral atom we have not processed yet.
    // RDKit❗❌:       // If we can't find non-chiral atom, use the chiral atom with
    // RDKit❗❌:       // the lowest rank (we are guaranteed to find an unprocessed atom).
    // RDKit❗❌:       unsigned nextRank = nAtoms + 1;
    // RDKit❗❌:       for (auto atom : mol.atoms()) {
    // RDKit❗❌:         if (colors[atom->getIdx()] == Canon::WHITE_NODE) {
    // RDKit❗❌:           if (atom->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW &&
    // RDKit❗❌:               atom->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW) {
    // RDKit❗❌:             nextAtomIdx = atom->getIdx();
    // RDKit❗❌:             break;
    // RDKit❗❌:           }
    // RDKit❗❌:           if (ranks[atom->getIdx()] < nextRank) {
    // RDKit❗❌:             nextRank = ranks[atom->getIdx()];
    // RDKit❗❌:             nextAtomIdx = atom->getIdx();
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     subSmi = FragmentSmartsConstruct(mol, nextAtomIdx, colors, ranks, params,
    // RDKit❗❌:                                      atomOrdering, bondOrdering, atomsInPlay,
    // RDKit❗❌:                                      bondsInPlay);
    // RDKit❗❌:     res += subSmi;
    // RDKit❗❌:
    // RDKit❗❌:     colorIt = std::find(colors.begin(), colors.end(), Canon::WHITE_NODE);
    // RDKit❗❌:     if (colorIt != colors.end()) {
    // RDKit❗❌:       res += ".";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   inmol.setProp(common_properties::_smilesAtomOutputOrder, atomOrdering, true);
    // RDKit❗❌:   inmol.setProp(common_properties::_smilesBondOutputOrder, bondOrdering, true);
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // Same source ROMol working copy, physical atom-index ranks and repeated
    // first-white/root/non-tetrahedral scan. No global chiral sort or neighbor
    // ordering approximation. Each component uses the complete source fragment.
    // The immutable detached API returns the source output-order effect proof;
    // its authorized caller records properties, rather than granting runtime
    // commit authority or pretending the borrowed input was mutated here.
    // Cost ❌: source O(V) scans/components and O(V+E+query/property bytes) copy;
    // known packed bitmap/Vec and Fragment common-field transport costs remain.
    if params
        .rooted_at_atom
        .is_some_and(|i| i >= query.num_atoms())
    {
        return Err(SmartsWriteError::RootedAtomOutOfRange {
            atom: params.rooted_at_atom.unwrap(),
        });
    }
    let mut mol = query.clone();
    let n = u32::try_from(mol.num_atoms()).map_err(|_| SmartsWriteError::SourceAtomCount {
        count: mol.num_atoms(),
    })?;
    let ranks = mol
        .atoms()
        .iter()
        .map(|a| {
            u32::try_from(a.id().index()).map_err(|_| SmartsWriteError::SourceAtomCount {
                count: mol.num_atoms(),
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    if params.isomeric_smiles {
        mol.set_prop("_doIsoSmiles", cosmolkit_model::PropertyValue::Int(1))?;
    }
    let mut result = SmartsWriteResult::default();
    let mut rings =
        cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 0, 0);
    while colors.contains(&cosmolkit_smiles::AtomColor::White) {
        let mut next = 0;
        if params
            .rooted_at_atom
            .is_some_and(|i| colors[i] == cosmolkit_smiles::AtomColor::White)
        {
            next = params.rooted_at_atom.unwrap();
        } else {
            let mut rank = n.wrapping_add(1);
            for atom in mol.atoms() {
                let i = atom.id().index();
                if colors[i] != cosmolkit_smiles::AtomColor::White {
                    continue;
                }
                if !matches!(
                    atom.chiral_tag(),
                    ChiralTag::TetrahedralCcw | ChiralTag::TetrahedralCw
                ) {
                    next = i;
                    break;
                }
                if ranks[i] < rank {
                    rank = ranks[i];
                    next = i;
                }
            }
        }
        let fragment = fragment_smarts_construct(
            &mut mol,
            &mut rings,
            next,
            &mut colors,
            &ranks,
            params,
            &mut result.atom_ordering,
            &mut result.bond_ordering,
            atoms_in_play,
            bonds_in_play,
        )?;
        result.smarts.extend_bytes(fragment.as_bytes());
        if colors.contains(&cosmolkit_smiles::AtomColor::White) {
            result.smarts.push_byte(b'.');
        }
    }
    // At the exact source inmol.setProp point, retain computed bookkeeping
    // and propagate its real conversion errors. The typed vectors above are
    // the source getter values; generic UInt-vector RDValue stays unmodeled.
    // This preserves the immutable detached adapter while mutable source
    // callers install these exact effects before invoking getCXExtensions.
    // Known extra cost: detached molecule-property copy versus Native's
    // direct computed-list reads/writes, inherited by the enclosing owner.
    let mut source_properties = query.source_molecule_properties();
    source_properties.register_transient_computed_name("_smilesAtomOutputOrder")?;
    source_properties.register_transient_computed_name("_smilesBondOutputOrder")?;
    result.source_properties = Some(source_properties);
    result.source_orders_written = true;
    Ok(result)
}

#[cfg(all(test, feature = "smiles-integration"))]
mod mol_to_smarts_source_tests {
    use super::*;
    use cosmolkit_model::{BondSpec, PropertyValue, QueryAtomIdentity, SourceAtomValenceFacts};
    use cosmolkit_smiles::AtomColor;
    fn graph(nums: &[u8], edges: &[(usize, usize)]) -> QueryGraph {
        QueryGraph::from_parts(
            nums.iter()
                .enumerate()
                .map(|(i, &n)| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(n),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(n)),
                    )
                })
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn write(q: &QueryGraph, p: &SmartsWriteParams) -> Result<SmartsWriteResult, SmartsWriteError> {
        mol_to_smarts_source(
            q,
            p,
            vec![AtomColor::White; q.num_atoms()],
            &vec![true; q.num_atoms()],
            None,
        )
    }
    fn preserve_stereo(q: &mut QueryGraph) {
        q.set_prop("_StereochemDone", PropertyValue::Bool(false))
            .unwrap();
    }
    #[test]
    fn physical_non_chiral_scan_precedes_lower_rank_tetrahedral_component() {
        let mut q = graph(&[6, 8, 7], &[]);
        q.atom_mut(0)
            .unwrap()
            .set_chiral_tag(ChiralTag::TetrahedralCw);
        preserve_stereo(&mut q);
        let before = q.clone();
        let out = write(&q, &Default::default()).unwrap();
        assert_eq!(
            out.atom_ordering,
            vec![AtomId::new(1), AtomId::new(2), AtomId::new(0)]
        );
        assert_eq!(out.smarts, PropertyText::from("[#8].[#7].[#6@@]"));
        assert!(out.source_orders_written);
        assert_eq!(q, before);
    }
    #[test]
    fn all_tetrahedral_remaining_atoms_use_lowest_actual_index_rank() {
        let mut q = graph(&[6, 6], &[]);
        for atom in q.atoms_mut() {
            atom.set_chiral_tag(ChiralTag::TetrahedralCw);
        }
        preserve_stereo(&mut q);
        let out = write(&q, &Default::default()).unwrap();
        assert_eq!(out.atom_ordering, vec![AtomId::new(0), AtomId::new(1)]);
        assert_eq!(out.smarts, PropertyText::from("[#6@@].[#6@@]"));
    }
    #[test]
    fn non_tetrahedral_chiral_tag_is_a_source_non_tetrahedral_start_candidate() {
        let mut q = graph(&[6, 15], &[]);
        q.atom_mut(0)
            .unwrap()
            .set_chiral_tag(ChiralTag::TetrahedralCw);
        q.atom_mut(1)
            .unwrap()
            .set_chiral_tag(ChiralTag::SquarePlanar);
        preserve_stereo(&mut q);
        let out = write(&q, &Default::default()).unwrap();
        assert_eq!(out.atom_ordering, vec![AtomId::new(1), AtomId::new(0)]);
    }
    #[test]
    fn rooted_component_is_first_then_source_physical_scan_resumes() {
        let q = graph(&[6, 8, 7], &[]);
        let mut p = SmartsWriteParams::default();
        p.rooted_at_atom = Some(2);
        let out = write(&q, &p).unwrap();
        assert_eq!(out.smarts, PropertyText::from("[#7].[#6].[#8]"));
        assert_eq!(
            out.atom_ordering,
            vec![AtomId::new(2), AtomId::new(0), AtomId::new(1)]
        );
    }
    #[test]
    fn already_black_root_is_ignored_and_selected_component_has_no_trailing_dot() {
        let q = graph(&[6, 8, 7], &[]);
        let mut p = SmartsWriteParams::default();
        p.rooted_at_atom = Some(0);
        let out = mol_to_smarts_source(
            &q,
            &p,
            vec![AtomColor::Black, AtomColor::White, AtomColor::White],
            &[false, true, true],
            Some(&[]),
        )
        .unwrap();
        assert_eq!(out.smarts, PropertyText::from("[#8].[#7]"));
        assert_eq!(out.atom_ordering, vec![AtomId::new(1), AtomId::new(2)]);
        assert!(out.bond_ordering.is_empty());
    }
    #[test]
    fn no_white_atoms_skip_fragment_cache_but_record_empty_source_orders() {
        let q = graph(&[255], &[]);
        let out = mol_to_smarts_source(
            &q,
            &Default::default(),
            vec![AtomColor::Black],
            &[false],
            None,
        )
        .unwrap();
        assert!(out.smarts.is_empty());
        assert!(out.source_orders_written);
        assert!(out.atom_ordering.is_empty());
        assert_eq!(
            q.atom(0).unwrap().source_valence_facts(),
            SourceAtomValenceFacts::UNINITIALIZED
        );
    }
    #[test]
    fn rooted_index_precondition_is_checked_even_when_no_atom_is_white() {
        let q = graph(&[255], &[]);
        let mut p = SmartsWriteParams::default();
        p.rooted_at_atom = Some(1);
        assert!(matches!(
            mol_to_smarts_source(&q, &p, vec![AtomColor::Black], &[false], None),
            Err(SmartsWriteError::RootedAtomOutOfRange { atom: 1 })
        ));
    }
    #[test]
    fn false_isomeric_option_preserves_existing_owner_presence_without_mutating_input() {
        let mut q = graph(&[6], &[]);
        q.atom_mut(0)
            .unwrap()
            .set_chiral_tag(ChiralTag::TetrahedralCw);
        preserve_stereo(&mut q);
        q.set_prop("_doIsoSmiles", PropertyValue::Bool(false))
            .unwrap();
        let before = q.clone();
        let mut p = SmartsWriteParams::default();
        p.isomeric_smiles = false;
        assert_eq!(write(&q, &p).unwrap().smarts, PropertyText::from("[#6@@]"));
        assert_eq!(q, before);
        q.clear_prop("_doIsoSmiles").unwrap();
        assert_eq!(write(&q, &p).unwrap().smarts, PropertyText::from("[#6]"));
    }
    #[test]
    fn public_output_route_uses_real_canonical_ring_stack_and_preserves_input() {
        let q = graph(&[6; 3], &[(0, 1), (1, 2), (2, 0)]);
        let before = q.clone();
        let out = query_graph_to_smarts_output(&q, &Default::default()).unwrap();
        assert_eq!(out.text, PropertyText::from("[#6]1-[#6]-[#6]-1"));
        assert_eq!(out.atom_order.len(), 3);
        assert_eq!(out.bond_order.len(), 3);
        assert!(out.source_orders_written);
        assert_eq!(q, before);
    }
}

#[cfg(feature = "smiles-integration")]
fn mol_to_smarts_wrapper_source(
    query: &QueryGraph,
    params: &SmartsWriteParams,
) -> Result<SmartsWriteResult, SmartsWriteError> {
    // RDKit❗❌: std::string MolToSmarts(const ROMol &mol, const SmilesWriteParams &ps) {
    // RDKit❗❌:   const unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗❌:   if (!nAtoms) {
    // RDKit❗❌:     return "";
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<AtomColors> colors(nAtoms, Canon::WHITE_NODE);
    // RDKit❗❌:   boost::dynamic_bitset<> atomsInPlay(nAtoms);
    // RDKit❗❌:   atomsInPlay.set();  // all atoms are in play
    // RDKit❗❌:   return molToSmarts(mol, ps, std::move(colors), atomsInPlay, nullptr);
    // RDKit❗❌: }
    // Empty source returns before rooted-index checks, graph copy, properties
    // or output-order recording. Nonempty colors/atom mask include every atom;
    // bondsInPlay is the actual absent source pointer. Reuse one molToSmarts.
    // Known ❌: byte mask versus packed bitset and inherited writer costs;
    // no extra BTreeSet selection/ring renumbering or duplicate traversal.
    let n = query.num_atoms();
    if n == 0 {
        return Ok(SmartsWriteResult::default());
    }
    mol_to_smarts_source(
        query,
        params,
        vec![cosmolkit_smiles::AtomColor::White; n],
        &vec![true; n],
        None,
    )
}

#[cfg(all(test, feature = "smiles-integration"))]
mod mol_to_smarts_wrapper_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue};
    fn q(n: usize) -> QueryGraph {
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
    #[test]
    fn empty_early_return_skips_root_validation_copy_properties_and_orders() {
        let mut graph = q(0);
        graph
            .set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        let before = graph.clone();
        let mut p = SmartsWriteParams::default();
        p.rooted_at_atom = Some(999);
        let out = mol_to_smarts_wrapper_source(&graph, &p).unwrap();
        assert!(out.smarts.is_empty());
        assert!(!out.source_orders_written);
        assert!(out.atom_ordering.is_empty());
        assert_eq!(graph, before);
    }
    #[test]
    fn nonempty_does_not_skip_root_precondition() {
        let mut p = SmartsWriteParams::default();
        p.rooted_at_atom = Some(1);
        assert!(matches!(
            mol_to_smarts_wrapper_source(&q(1), &p),
            Err(SmartsWriteError::RootedAtomOutOfRange { atom: 1 })
        ));
    }
    #[test]
    fn all_components_and_atoms_are_in_play_with_source_physical_order() {
        let graph = q(3);
        let out = mol_to_smarts_wrapper_source(&graph, &Default::default()).unwrap();
        assert_eq!(out.smarts, PropertyText::from("[#6].[#6].[#6]"));
        assert_eq!(
            out.atom_ordering,
            vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]
        );
        assert!(out.bond_ordering.is_empty());
        assert!(out.source_orders_written);
    }
    #[test]
    fn public_text_and_output_apis_share_source_empty_and_nonempty_order_semantics() {
        for n in [0, 1, 3] {
            let graph = q(n);
            let out = query_graph_to_smarts_output(&graph, &Default::default()).unwrap();
            assert_eq!(
                out.text,
                query_graph_to_smarts(&graph, &Default::default()).unwrap()
            );
            assert_eq!(out.source_orders_written, n != 0);
            assert_eq!(out.atom_order.len(), n);
        }
    }
    #[test]
    fn flags_and_root_are_forwarded_without_wrapper_overrides() {
        let mut graph = q(2);
        graph.atom_mut(0).unwrap().set_atom_map(Some(7));
        graph.atom_mut(1).unwrap().set_atom_map(Some(8));
        let mut p = SmartsWriteParams::default();
        p.rooted_at_atom = Some(1);
        p.include_atom_maps = false;
        assert_eq!(
            mol_to_smarts_wrapper_source(&graph, &p).unwrap().smarts,
            PropertyText::from("[#6].[#6]")
        );
        p.include_atom_maps = true;
        assert_eq!(
            mol_to_smarts_wrapper_source(&graph, &p).unwrap().smarts,
            PropertyText::from("[#6:8].[#6:7]")
        );
    }
    #[test]
    fn source_cache_and_stereo_working_effects_never_modify_borrowed_input() {
        let mut graph = q(1);
        graph
            .atom_mut(0)
            .unwrap()
            .set_chiral_tag(ChiralTag::TetrahedralCw);
        graph
            .set_prop("_StereochemDone", PropertyValue::Bool(false))
            .unwrap();
        let before = graph.clone();
        let one = mol_to_smarts_wrapper_source(&graph, &Default::default()).unwrap();
        let two = mol_to_smarts_wrapper_source(&graph, &Default::default()).unwrap();
        assert_eq!(one.smarts, PropertyText::from("[#6@@]"));
        assert_eq!(one.smarts, two.smarts);
        assert_eq!(graph, before);
        assert!(graph.atom(0).unwrap().prop("_qatomHasStereoSet").is_none());
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod get_sorted_mapped_indexes_source_tests {
    use super::*;
    fn ids(v: &[usize]) -> Vec<AtomId> {
        v.iter().copied().map(AtomId::new).collect()
    }
    #[test]
    fn empty_input_stays_empty_without_reading_reverse_vector() {
        assert_eq!(
            get_sorted_mapped_indexes(&[], &[]).unwrap(),
            Vec::<usize>::new()
        );
        assert!(get_sorted_mapped_indexes(&[], &[9, 8]).unwrap().is_empty());
    }
    #[test]
    fn reverse_permutation_is_applied_before_sorting() {
        assert_eq!(
            get_sorted_mapped_indexes(&ids(&[2, 0, 1]), &[2, 0, 1]).unwrap(),
            vec![0, 1, 2]
        );
    }
    #[test]
    fn repeated_input_members_remain_repeated() {
        assert_eq!(
            get_sorted_mapped_indexes(&ids(&[1, 1, 0]), &[9, 3]).unwrap(),
            vec![3, 3, 9]
        );
    }
    #[test]
    fn sparse_default_zero_mappings_are_not_filtered() {
        assert_eq!(
            get_sorted_mapped_indexes(&ids(&[0, 1, 2]), &[0, 0, 4]).unwrap(),
            vec![0, 0, 4]
        );
    }
    #[test]
    fn mapped_value_order_does_not_assume_input_atom_index_order() {
        let input = ids(&[3, 0]);
        let reverse = vec![9, 100, 100, 2];
        let before = (input.clone(), reverse.clone());
        assert_eq!(
            get_sorted_mapped_indexes(&input, &reverse).unwrap(),
            vec![2, 9]
        );
        assert_eq!((input, reverse), before);
    }
    #[test]
    fn invalid_reverse_access_is_a_structural_error_not_zero() {
        assert!(
            matches!(get_sorted_mapped_indexes(&ids(&[0,2,1]),&[7,9]),Err(SmartsWriteError::CxSourceAtomOutOfRange{atom,atom_count:2})if atom==AtomId::new(2))
        );
        assert!(matches!(
            get_sorted_mapped_indexes(&ids(&[0]), &[]),
            Err(SmartsWriteError::CxSourceAtomOutOfRange { atom_count: 0, .. })
        ));
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod get_sorted_stereo_groups_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn group(
        kind: StereoGroupKind,
        atoms: &[usize],
        bonds: &[usize],
        read: u32,
        write: u32,
    ) -> StereoGroup {
        StereoGroup::new(
            kind,
            atoms.iter().copied().map(AtomId::new).collect(),
            bonds.iter().copied().map(BondId::new).collect(),
        )
        .expect("valid distinct stereo members")
        .with_id(read)
        .with_write_id(write)
    }
    fn graph(groups: Vec<StereoGroup>, direction: BondDirection) -> QueryGraph {
        QueryGraph::from_parts(
            (0..3)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            vec![
                QueryBond::new(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
                ),
                QueryBond::new(
                    BondId::new(1),
                    BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single)
                        .with_direction(direction),
                ),
            ],
            [],
            vec![],
            vec![],
            groups,
        )
        .unwrap()
    }
    fn sorted(q: &QueryGraph, rev: &[usize]) -> Vec<(StereoGroup, Vec<usize>)> {
        get_sorted_stereo_groups_and_indices(q, rev, &Default::default()).unwrap()
    }
    #[test]
    fn empty_collected_groups_are_filtered_after_source_collection() {
        let q = graph(
            vec![
                group(StereoGroupKind::Or, &[], &[], 7, 4),
                group(StereoGroupKind::And, &[], &[0], 8, 5),
            ],
            BondDirection::None,
        );
        assert!(sorted(&q, &[0, 1, 2]).is_empty());
    }
    #[test]
    fn source_group_type_precedes_write_id_and_members() {
        let q = graph(
            vec![
                group(StereoGroupKind::And, &[0], &[], 1, 1),
                group(StereoGroupKind::Or, &[2], &[], 2, 99),
                group(StereoGroupKind::Absolute, &[1], &[], 3, 200),
            ],
            BondDirection::None,
        );
        assert_eq!(
            sorted(&q, &[0, 1, 2])
                .iter()
                .map(|(g, _)| g.kind())
                .collect::<Vec<_>>(),
            vec![
                StereoGroupKind::Absolute,
                StereoGroupKind::Or,
                StereoGroupKind::And
            ]
        );
    }
    #[test]
    fn same_kind_orders_by_write_id_before_mapped_members() {
        let q = graph(
            vec![
                group(StereoGroupKind::Or, &[0], &[], 1, 9),
                group(StereoGroupKind::Or, &[2], &[], 2, 2),
            ],
            BondDirection::None,
        );
        assert_eq!(
            sorted(&q, &[0, 1, 2])
                .iter()
                .map(|(g, _)| g.write_id())
                .collect::<Vec<_>>(),
            vec![2, 9]
        );
    }
    #[test]
    fn equal_write_ids_use_mapped_members_without_read_id_tiebreak() {
        let q = graph(
            vec![
                group(StereoGroupKind::Or, &[0], &[], 1, 7),
                group(StereoGroupKind::Or, &[2], &[], 99, 7),
            ],
            BondDirection::None,
        );
        let out = sorted(&q, &[5, 0, 2]);
        assert_eq!(out[0].0.id(), Some(99));
        assert_eq!(out[0].1, vec![2]);
        assert_eq!(out[1].1, vec![5]);
    }
    #[test]
    fn duplicate_members_and_sparse_zero_reverse_entries_survive() {
        let q = graph(
            vec![group(StereoGroupKind::And, &[2, 0, 1], &[], 4, 0)],
            BondDirection::None,
        );
        assert_eq!(sorted(&q, &[0, 0, 4])[0].1, vec![0, 0, 4]);
    }
    #[test]
    fn bond_only_group_uses_actual_neighbor_wedge_direction() {
        let q = graph(
            vec![group(StereoGroupKind::And, &[], &[0], 4, 3)],
            BondDirection::BeginDash,
        );
        let out = sorted(&q, &[2, 0, 1]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].1, vec![0]);
        assert_eq!(out[0].0.bonds(), &[BondId::new(0)]);
    }
    #[test]
    fn actual_atrop_info_reaches_canonical_group_collection() {
        let q = graph(
            vec![group(StereoGroupKind::Or, &[], &[0], 4, 3)],
            BondDirection::None,
        );
        let update = cosmolkit_core::AtropisomerWedgeUpdate {
            bond: BondId::new(1),
            begin: AtomId::new(1),
            end: AtomId::new(2),
            direction: BondDirection::None,
            atropisomer_bond: BondId::new(0),
        };
        let map = cosmolkit_core::WedgeAssignments::from_atropisomer_wedge_assignment(
            cosmolkit_core::AtropisomerWedgeAssignment {
                bond_updates: vec![update],
                source_map_writes: vec![BondId::new(1)],
                diagnostics: vec![],
            },
        );
        assert_eq!(
            get_sorted_stereo_groups_and_indices(&q, &[0, 2, 1], &map).unwrap()[0].1,
            vec![2]
        );
    }
    #[test]
    fn sorting_preserves_group_read_write_ids_and_input_query() {
        let q = graph(
            vec![
                group(StereoGroupKind::Or, &[2], &[], 9, 0),
                group(StereoGroupKind::Or, &[0], &[], 8, 0),
            ],
            BondDirection::None,
        );
        let before = q.clone();
        let out = sorted(&q, &[0, 1, 2]);
        assert_eq!(
            out.iter()
                .map(|(g, _)| (g.id(), g.write_id()))
                .collect::<Vec<_>>(),
            vec![(Some(8), 0), (Some(9), 0)]
        );
        assert_eq!(q, before);
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod enhanced_stereo_block_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn group(
        kind: StereoGroupKind,
        atoms: &[usize],
        bonds: &[usize],
        read: u32,
        write: u32,
    ) -> StereoGroup {
        StereoGroup::new(
            kind,
            atoms.iter().copied().map(AtomId::new).collect(),
            bonds.iter().copied().map(BondId::new).collect(),
        )
        .expect("valid distinct stereo members")
        .with_id(read)
        .with_write_id(write)
    }
    fn graph(groups: Vec<StereoGroup>) -> QueryGraph {
        QueryGraph::from_parts(
            (0..3)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            vec![
                QueryBond::new(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
                ),
                QueryBond::new(
                    BondId::new(1),
                    BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
                ),
            ],
            [],
            vec![],
            vec![],
            groups,
        )
        .unwrap()
    }
    fn write(q: &QueryGraph, order: &[usize]) -> String {
        write_query_cx_enhanced_stereo(
            q,
            &order.iter().copied().map(AtomId::new).collect::<Vec<_>>(),
            &Default::default(),
        )
        .unwrap()
    }
    #[test]
    fn no_groups_return_before_reverse_order_index_access() {
        assert!(write(&graph(vec![]), &[999]).is_empty());
    }
    #[test]
    fn present_groups_make_invalid_reverse_access_structural() {
        let q = graph(vec![group(StereoGroupKind::Or, &[], &[], 0, 0)]);
        assert!(matches!(
            write_query_cx_enhanced_stereo(&q, &[AtomId::new(999)], &Default::default()),
            Err(SmartsWriteError::CxSourceAtomOutOfRange { atom_count: 3, .. })
        ));
    }
    #[test]
    fn source_atom_output_order_controls_emitted_group_members() {
        let q = graph(vec![group(StereoGroupKind::Absolute, &[0, 2], &[], 99, 8)]);
        assert_eq!(write(&q, &[1, 2, 0]), "a:1,2");
    }
    #[test]
    fn fragment_unselected_members_keep_native_zero_mapping() {
        let q = graph(vec![group(StereoGroupKind::Absolute, &[0, 2], &[], 99, 8)]);
        assert_eq!(write(&q, &[1]), "a:0,0");
    }
    #[test]
    fn or_and_missing_ids_use_separate_source_sequences() {
        let q = graph(vec![
            group(StereoGroupKind::And, &[1], &[], 5, 0),
            group(StereoGroupKind::Or, &[0], &[], 17, 0),
        ]);
        assert_eq!(write(&q, &[0, 1, 2]), "o1:0,&1:1");
    }
    #[test]
    fn explicit_write_id_is_preserved_independently_of_read_id() {
        let q = graph(vec![group(StereoGroupKind::Or, &[0], &[], 99, 9)]);
        assert_eq!(write(&q, &[0, 1, 2]), "o9:0");
    }
    #[test]
    fn duplicate_write_ids_are_reassigned_after_sort_without_resorting() {
        let q = graph(vec![
            group(StereoGroupKind::Or, &[1], &[], 1, 5),
            group(StereoGroupKind::Or, &[0], &[], 2, 5),
        ]);
        assert_eq!(write(&q, &[0, 1, 2]), "o5:0,o1:1");
    }
    #[test]
    fn actual_atrop_info_participates_in_bond_only_group_emission() {
        let q = graph(vec![group(StereoGroupKind::Or, &[], &[0], 4, 3)]);
        let update = cosmolkit_core::AtropisomerWedgeUpdate {
            bond: BondId::new(1),
            begin: AtomId::new(1),
            end: AtomId::new(2),
            direction: BondDirection::None,
            atropisomer_bond: BondId::new(0),
        };
        let map = cosmolkit_core::WedgeAssignments::from_atropisomer_wedge_assignment(
            cosmolkit_core::AtropisomerWedgeAssignment {
                bond_updates: vec![update],
                source_map_writes: vec![BondId::new(1)],
                diagnostics: vec![],
            },
        );
        assert_eq!(
            write_query_cx_enhanced_stereo(
                &q,
                &[AtomId::new(0), AtomId::new(2), AtomId::new(1)],
                &map
            )
            .unwrap(),
            "o3:2"
        );
    }
    #[test]
    fn nonempty_group_collection_with_no_marked_members_emits_empty() {
        let q = graph(vec![group(StereoGroupKind::And, &[], &[0], 8, 5)]);
        assert!(write(&q, &[0, 1, 2]).is_empty());
    }
    #[test]
    fn duplicate_members_and_final_comma_follow_source_and_input_stays_unchanged() {
        let q = graph(vec![group(StereoGroupKind::Or, &[0, 1], &[], 88, 0)]);
        let before = q.clone();
        assert_eq!(write(&q, &[0, 1, 2]), "o1:0,1");
        assert_eq!(q, before);
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod sgroup_hierarchy_source_tests {
    use super::*;
    use cosmolkit_model::{PropertyValue, SubstanceGroupId};
    fn group(i: usize) -> SubstanceGroup {
        SubstanceGroup::new(SubstanceGroupId::new(i), SubstanceGroupKind::Data)
    }
    fn prop(mut g: SubstanceGroup, k: &str, v: PropertyValue) -> SubstanceGroup {
        g.set_prop(k, v).unwrap();
        g
    }
    fn emitted(i: usize, index: u32, output: u32) -> SubstanceGroup {
        prop(
            prop(group(i), "index", PropertyValue::UInt(index)),
            "_cxsmilesOutputIndex",
            PropertyValue::UInt(output),
        )
    }
    fn parent(g: SubstanceGroup, index: u32) -> SubstanceGroup {
        prop(g, "PARENT", PropertyValue::UInt(index))
    }
    fn state(groups: Vec<SubstanceGroup>) -> QueryCxSgroupState {
        QueryCxSgroupState {
            groups,
            properties: Default::default(),
        }
    }
    fn cleared(g: &SubstanceGroup) -> bool {
        !g.props().contains_key(b"_cxsmilesOutputIndex".as_slice())
    }
    #[test]
    fn empty_groups_return_without_touching_molecule_properties() {
        let mut s = state(vec![]);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::String("bad".into()))
            .unwrap();
        let before = s.properties.clone();
        assert!(write_query_cx_sgroup_hierarchy(&mut s).unwrap().is_empty());
        assert_eq!(s.properties, before);
    }
    #[test]
    fn no_parent_still_clears_actual_output_indexes() {
        let mut s = state(vec![emitted(0, 10, 7), emitted(1, 20, 3)]);
        assert!(write_query_cx_sgroup_hierarchy(&mut s).unwrap().is_empty());
        assert!(s.groups.iter().all(cleared));
    }
    #[test]
    fn actual_source_and_output_indexes_determine_hierarchy() {
        let mut s = state(vec![emitted(0, 10, 7), parent(emitted(1, 20, 3), 10)]);
        assert_eq!(write_query_cx_sgroup_hierarchy(&mut s).unwrap(), "SgH:7:3");
        assert!(s.groups.iter().all(cleared));
    }
    #[test]
    fn group_types_never_fabricate_missing_output_index_properties() {
        let mut s = state(vec![group(0), parent(group(1), 0)]);
        assert!(write_query_cx_sgroup_hierarchy(&mut s).unwrap().is_empty());
    }
    #[test]
    fn parent_output_keys_sort_but_children_keep_physical_order() {
        let mut s = state(vec![
            emitted(0, 10, 7),
            emitted(1, 20, 2),
            parent(emitted(2, 30, 5), 10),
            parent(emitted(3, 40, 4), 20),
            parent(emitted(4, 50, 3), 10),
        ]);
        assert_eq!(
            write_query_cx_sgroup_hierarchy(&mut s).unwrap(),
            "SgH:2:4,7:5.3"
        );
    }
    #[test]
    fn duplicate_source_indexes_overwrite_map_last() {
        let mut s = state(vec![
            emitted(0, 10, 7),
            emitted(1, 10, 2),
            parent(emitted(2, 30, 5), 10),
        ]);
        assert_eq!(write_query_cx_sgroup_hierarchy(&mut s).unwrap(), "SgH:2:5");
    }
    #[test]
    fn default_source_index_is_the_physical_group_index() {
        let mut s = state(vec![
            prop(group(0), "_cxsmilesOutputIndex", PropertyValue::UInt(8)),
            parent(
                prop(group(1), "_cxsmilesOutputIndex", PropertyValue::UInt(9)),
                0,
            ),
        ]);
        assert_eq!(write_query_cx_sgroup_hierarchy(&mut s).unwrap(), "SgH:8:9");
    }
    #[test]
    fn bad_parent_on_nonemitted_group_is_not_skipped() {
        let mut s = state(vec![
            emitted(0, 10, 7),
            prop(group(1), "PARENT", PropertyValue::String("bad".into())),
        ]);
        assert!(matches!(
            write_query_cx_sgroup_hierarchy(&mut s),
            Err(SmartsWriteError::CxSgroupPropertyUInt {
                property: "PARENT",
                ..
            })
        ));
        assert!(cleared(&s.groups[0]));
    }
    #[test]
    fn unmatched_parent_short_circuits_bad_child_index() {
        let g = parent(
            prop(group(1), "index", PropertyValue::String("bad".into())),
            999,
        );
        let mut s = state(vec![emitted(0, 10, 7), g]);
        assert!(write_query_cx_sgroup_hierarchy(&mut s).unwrap().is_empty());
    }
    #[test]
    fn matched_parent_reads_bad_child_index_even_if_child_not_emitted() {
        let g = parent(
            prop(group(1), "index", PropertyValue::String("bad".into())),
            10,
        );
        let mut s = state(vec![emitted(0, 10, 7), g]);
        assert!(matches!(
            write_query_cx_sgroup_hierarchy(&mut s),
            Err(SmartsWriteError::CxSgroupPropertyUInt {
                property: "index",
                ..
            })
        ));
        assert!(cleared(&s.groups[0]));
    }
    #[test]
    fn no_parent_does_not_read_unemitted_bad_index() {
        let mut s = state(vec![prop(
            group(0),
            "index",
            PropertyValue::String("bad".into()),
        )]);
        assert!(write_query_cx_sgroup_hierarchy(&mut s).unwrap().is_empty());
    }
    #[test]
    fn source_index_conversion_precedes_output_conversion_and_clear() {
        let bad = prop(
            prop(group(1), "index", PropertyValue::String("bad-index".into())),
            "_cxsmilesOutputIndex",
            PropertyValue::String("bad-output".into()),
        );
        let mut s = state(vec![emitted(0, 10, 7), bad]);
        assert!(matches!(
            write_query_cx_sgroup_hierarchy(&mut s),
            Err(SmartsWriteError::CxSgroupPropertyUInt {
                property: "index",
                ..
            })
        ));
        assert!(cleared(&s.groups[0]));
        assert!(!cleared(&s.groups[1]));
    }
    #[test]
    fn output_conversion_failure_keeps_current_and_later_group_properties() {
        let mut s = state(vec![
            emitted(0, 10, 7),
            prop(
                group(1),
                "_cxsmilesOutputIndex",
                PropertyValue::String("bad".into()),
            ),
            emitted(2, 30, 2),
        ]);
        assert!(matches!(
            write_query_cx_sgroup_hierarchy(&mut s),
            Err(SmartsWriteError::CxSgroupPropertyUInt {
                property: "_cxsmilesOutputIndex",
                ..
            })
        ));
        assert!(cleared(&s.groups[0]));
        assert!(!cleared(&s.groups[1]));
        assert!(!cleared(&s.groups[2]));
    }
    #[test]
    fn computed_list_failure_preserves_prior_clears_and_current_output_property() {
        let mut s = state(vec![
            emitted(0, 10, 7),
            prop(emitted(1, 20, 3), "__computedProps", PropertyValue::Int(9)),
            emitted(2, 30, 2),
        ]);
        assert!(matches!(
            write_query_cx_sgroup_hierarchy(&mut s),
            Err(SmartsWriteError::CxSgroupPropertyWrite {
                property: "_cxsmilesOutputIndex",
                source: cosmolkit_model::MoleculePropertyError::ComputedListKind(_),
                ..
            })
        ));
        assert!(cleared(&s.groups[0]));
        assert!(!cleared(&s.groups[1]));
        assert!(!cleared(&s.groups[2]));
    }
    #[test]
    fn typed_parent_relation_projects_to_existing_external_index() {
        let mut s = state(vec![
            emitted(0, 11, 7),
            emitted(1, 22, 3).with_parent(SubstanceGroupId::new(0)),
        ]);
        assert_eq!(write_query_cx_sgroup_hierarchy(&mut s).unwrap(), "SgH:7:3");
        assert_eq!(s.groups[1].parent(), Some(SubstanceGroupId::new(0)));
    }
    #[test]
    fn clearing_output_property_removes_only_first_computed_membership() {
        let names = vec![
            "x".into(),
            "_cxsmilesOutputIndex".into(),
            "_cxsmilesOutputIndex".into(),
        ];
        let mut s = state(vec![prop(
            emitted(0, 11, 7),
            "__computedProps",
            PropertyValue::StringVector(names),
        )]);
        assert!(write_query_cx_sgroup_hierarchy(&mut s).unwrap().is_empty());
        assert_eq!(
            s.groups[0].props().get(b"__computedProps".as_slice()),
            Some(&PropertyValue::StringVector(vec![
                "x".into(),
                "_cxsmilesOutputIndex".into()
            ]))
        );
    }
    #[test]
    fn absent_output_index_does_not_read_bad_computed_list() {
        let mut s = state(vec![prop(
            group(0),
            "__computedProps",
            PropertyValue::Int(9),
        )]);
        assert!(write_query_cx_sgroup_hierarchy(&mut s).unwrap().is_empty());
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod polymer_sgroups_source_tests {
    use super::*;
    use cosmolkit_model::{
        AtomSpec, BondSpec, PropertyValue, SubstanceGroupId, replace_query_substance_groups,
    };
    use cosmolkit_types::Element;
    fn group(i: usize, kind: SubstanceGroupKind) -> SubstanceGroup {
        SubstanceGroup::new(SubstanceGroupId::new(i), kind).with_atoms(vec![AtomId::new(0)])
    }
    fn prop(mut g: SubstanceGroup, k: &str, v: PropertyValue) -> SubstanceGroup {
        g.set_prop(k, v).unwrap();
        g
    }
    fn graph(groups: Vec<SubstanceGroup>) -> QueryGraph {
        let mut q = QueryGraph::from_parts(
            (0..4)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            (0..3)
                .map(|i| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(i), AtomId::new(i + 1), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        replace_query_substance_groups(&mut q, groups).unwrap();
        q
    }
    fn write(q: &QueryGraph, s: &mut QueryCxSgroupState) -> PropertyText {
        write_query_cx_polymer_sgroups(
            q,
            s,
            &(0..4).map(AtomId::new).collect::<Vec<_>>(),
            &(0..3).map(BondId::new).collect::<Vec<_>>(),
        )
        .unwrap()
    }
    fn counter(s: &QueryCxSgroupState) -> Option<&PropertyValue> {
        s.properties.prop("_cxsmilesOutputIndex")
    }
    #[test]
    fn empty_groups_return_before_bad_counter_and_reverse_order_reads() {
        let q = graph(vec![]);
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::String("bad".into()))
            .unwrap();
        assert!(
            write_query_cx_polymer_sgroups(&q, &mut s, &[AtomId::new(999)], &[BondId::new(999)])
                .unwrap()
                .is_empty()
        );
        assert_eq!(counter(&s), Some(&PropertyValue::String("bad".into())));
    }
    #[test]
    fn molecule_counter_conversion_precedes_reverse_access() {
        let q = graph(vec![group(0, SubstanceGroupKind::StructuralRepeatUnit)]);
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::String("bad".into()))
            .unwrap();
        assert!(matches!(
            write_query_cx_polymer_sgroups(&q, &mut s, &[AtomId::new(999)], &[]),
            Err(SmartsWriteError::CxMoleculePropertyUInt { .. })
        ));
        assert!(
            !s.groups[0]
                .props()
                .contains_key(b"_cxsmilesOutputIndex".as_slice())
        );
    }
    #[test]
    fn raw_unknown_type_overrides_kind_and_skips_later_wrong_vector_tag() {
        let g = prop(
            prop(
                group(0, SubstanceGroupKind::StructuralRepeatUnit),
                "TYPE",
                PropertyValue::String("unknown".into()),
            ),
            "XBHEAD",
            PropertyValue::IntVector(vec![]),
        );
        let q = graph(vec![g]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert!(write(&q, &mut s).is_empty());
        assert_eq!(counter(&s), Some(&PropertyValue::UInt(0)));
    }
    #[test]
    fn raw_type_and_subtype_override_constructor_kind() {
        let g = prop(
            prop(
                group(0, SubstanceGroupKind::Data),
                "TYPE",
                PropertyValue::String("COP".into()),
            ),
            "SUBTYPE",
            PropertyValue::String("RAN".into()),
        );
        let q = graph(vec![g]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert_eq!(write(&q, &mut s), PropertyText::from("Sg:ran:0:::::"));
    }
    #[test]
    fn lexical_reverse_typemap_chooses_alt_for_unknown_copolymer_subtype() {
        let q = graph(vec![
            group(0, SubstanceGroupKind::Copolymer).with_subtype("unknown"),
        ]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert_eq!(write(&q, &mut s), PropertyText::from("Sg:alt:0:::::"));
    }
    #[test]
    fn raw_label_and_connect_take_precedence_and_keep_counted_binary_bytes() {
        let mut g = group(0, SubstanceGroupKind::StructuralRepeatUnit);
        g.set_label("typed");
        g.set_connection(SGroupConnection::HeadToTail);
        let mut label = PropertyText::new();
        label.extend_bytes(&[255, 0, b'.']);
        g = prop(
            prop(g, "LABEL", PropertyValue::String(label.clone())),
            "CONNECT",
            PropertyValue::String("HH".into()),
        );
        let q = graph(vec![g]);
        let mut s = QueryCxSgroupState::from_query(&q);
        let mut expected = PropertyText::from("Sg:n:0:");
        expected.extend_bytes(label.as_bytes());
        expected.extend_bytes(b":hh:::");
        assert_eq!(write(&q, &mut s), expected);
    }
    #[test]
    fn crossing_bonds_use_forward_order_and_odd_tail_positions() {
        let g = group(0, SubstanceGroupKind::StructuralRepeatUnit)
            .with_head_crossing_bonds(vec![BondId::new(0), BondId::new(1)])
            .with_crossing_bond_correspondence(vec![
                BondId::new(2),
                BondId::new(0),
                BondId::new(1),
                BondId::new(2),
            ]);
        let q = graph(vec![g]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert_eq!(
            write_query_cx_polymer_sgroups(
                &q,
                &mut s,
                &[AtomId::new(0)],
                &[BondId::new(2), BondId::new(0), BondId::new(1)]
            )
            .unwrap(),
            PropertyText::from("Sg:n:0:::2,0:2,1:")
        );
    }
    #[test]
    fn wrong_generic_vector_tag_fails_after_group_write_without_final_counter_commit() {
        let q = graph(vec![prop(
            group(0, SubstanceGroupKind::StructuralRepeatUnit),
            "XBHEAD",
            PropertyValue::IntVector(vec![]),
        )]);
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(7))
            .unwrap();
        assert!(matches!(
            write_query_cx_polymer_sgroups(&q, &mut s, &[AtomId::new(0)], &[]),
            Err(SmartsWriteError::CxSgroupVectorCast {
                property: "XBHEAD",
                ..
            })
        ));
        assert_eq!(
            s.groups[0].props().get(b"_cxsmilesOutputIndex".as_slice()),
            Some(&PropertyValue::UInt(7))
        );
        assert_eq!(counter(&s), Some(&PropertyValue::UInt(7)));
    }
    #[test]
    fn source_unsigned_output_counter_wraps_after_group_write() {
        let q = graph(vec![group(0, SubstanceGroupKind::StructuralRepeatUnit)]);
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(u32::MAX))
            .unwrap();
        write(&q, &mut s);
        assert_eq!(
            s.groups[0].props().get(b"_cxsmilesOutputIndex".as_slice()),
            Some(&PropertyValue::UInt(u32::MAX))
        );
        assert_eq!(counter(&s), Some(&PropertyValue::UInt(0)));
    }
    #[test]
    fn unused_native_reverse_bond_vector_is_still_checked_before_group_writes() {
        let q = graph(vec![group(0, SubstanceGroupKind::StructuralRepeatUnit)]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert!(matches!(
            write_query_cx_polymer_sgroups(&q, &mut s, &[], &[BondId::new(999)]),
            Err(SmartsWriteError::CxSourceBondOutOfRange { .. })
        ));
        assert!(
            !s.groups[0]
                .props()
                .contains_key(b"_cxsmilesOutputIndex".as_slice())
        );
    }
    #[test]
    fn empty_atom_members_follow_source_seek_overwrite() {
        let g = SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        );
        let q = graph(vec![g]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert_eq!(write(&q, &mut s), PropertyText::from("Sg:n:::::"));
    }
    #[test]
    fn hierarchy_consumes_actual_polymer_writes_and_input_query_stays_unchanged() {
        let q = graph(vec![
            group(0, SubstanceGroupKind::StructuralRepeatUnit),
            group(1, SubstanceGroupKind::StructuralRepeatUnit)
                .with_parent(SubstanceGroupId::new(0)),
        ]);
        let before = q.clone();
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(9))
            .unwrap();
        write(&q, &mut s);
        assert_eq!(write_query_cx_sgroup_hierarchy(&mut s).unwrap(), "SgH:9:10");
        assert!(
            s.groups
                .iter()
                .all(|g| !g.props().contains_key(b"_cxsmilesOutputIndex".as_slice()))
        );
        assert_eq!(counter(&s), Some(&PropertyValue::UInt(11)));
        assert_eq!(q, before);
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod data_sgroups_source_tests {
    use super::*;
    use cosmolkit_model::{
        AtomSpec, BondSpec, PropertyValue, SubstanceGroupId, replace_query_substance_groups,
    };
    use cosmolkit_types::Element;
    fn group(i: usize, kind: SubstanceGroupKind) -> SubstanceGroup {
        SubstanceGroup::new(SubstanceGroupId::new(i), kind).with_atoms(vec![AtomId::new(0)])
    }
    fn prop(mut g: SubstanceGroup, k: &str, v: PropertyValue) -> SubstanceGroup {
        g.set_prop(k, v).unwrap();
        g
    }
    fn graph(groups: Vec<SubstanceGroup>) -> QueryGraph {
        let mut q = QueryGraph::from_parts(
            (0..4)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            (0..3)
                .map(|i| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(i), AtomId::new(i + 1), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        replace_query_substance_groups(&mut q, groups).unwrap();
        q
    }
    fn write(q: &QueryGraph, s: &mut QueryCxSgroupState) -> PropertyText {
        write_query_cx_data_sgroups(q, s, &(0..4).map(AtomId::new).collect::<Vec<_>>()).unwrap()
    }
    fn counter(s: &QueryCxSgroupState) -> Option<&PropertyValue> {
        s.properties.prop("_cxsmilesOutputIndex")
    }
    #[test]
    fn empty_groups_return_before_bad_counter_and_reverse_access() {
        let q = graph(vec![]);
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::String("bad".into()))
            .unwrap();
        assert!(
            write_query_cx_data_sgroups(&q, &mut s, &[AtomId::new(999)])
                .unwrap()
                .is_empty()
        );
        assert_eq!(counter(&s), Some(&PropertyValue::String("bad".into())));
    }
    #[test]
    fn bad_counter_precedes_invalid_reverse_order_and_group_write() {
        let q = graph(vec![group(0, SubstanceGroupKind::Data)]);
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::String("bad".into()))
            .unwrap();
        assert!(matches!(
            write_query_cx_data_sgroups(&q, &mut s, &[AtomId::new(999)]),
            Err(SmartsWriteError::CxMoleculePropertyUInt { .. })
        ));
        assert!(
            !s.groups[0]
                .props()
                .contains_key(b"_cxsmilesOutputIndex".as_slice())
        );
    }
    #[test]
    fn raw_type_overrides_constructor_and_skips_wrong_datafields_tag() {
        let g = prop(
            prop(
                group(0, SubstanceGroupKind::Data),
                "TYPE",
                PropertyValue::String("SRU".into()),
            ),
            "DATAFIELDS",
            PropertyValue::String("wrong".into()),
        );
        let q = graph(vec![g]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert!(write(&q, &mut s).is_empty());
        assert_eq!(counter(&s), Some(&PropertyValue::UInt(0)));
    }
    #[test]
    fn generic_dat_is_actual_source_type() {
        let q = graph(vec![group(0, SubstanceGroupKind::Generic("DAT".into()))]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert_eq!(write(&q, &mut s), PropertyText::from("SgD:0::::::"));
        assert_eq!(counter(&s), Some(&PropertyValue::UInt(1)));
    }
    #[test]
    fn user_member_order_and_duplicates_are_not_canonicalized() {
        let q = graph(vec![group(0, SubstanceGroupKind::Data).with_atoms(vec![
            AtomId::new(2),
            AtomId::new(0),
            AtomId::new(2),
        ])]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert_eq!(
            write_query_cx_data_sgroups(&q, &mut s, &[AtomId::new(2), AtomId::new(0)]).unwrap(),
            PropertyText::from("SgD:0,1,0::::::")
        );
    }
    #[test]
    fn all_fields_keep_counted_binary_strings_and_datafields_vector_order() {
        let binary = PropertyText::from(vec![255, 0, b'.']);
        let mut g = group(0, SubstanceGroupKind::Data);
        for key in ["FIELDNAME", "QUERYOP", "FIELDINFO", "FIELDTAG"] {
            g = prop(g, key, PropertyValue::String(binary.clone()));
        }
        g = prop(
            g,
            "DATAFIELDS",
            PropertyValue::StringVector(vec![binary.clone(), "tail".into()]),
        );
        let q = graph(vec![g]);
        let mut s = QueryCxSgroupState::from_query(&q);
        let mut expected = PropertyText::from("SgD:0:");
        expected.extend_bytes(binary.as_bytes());
        expected.push_byte(b':');
        expected.extend_bytes(binary.as_bytes());
        expected.extend_bytes(b",tail:");
        for _ in 0..3 {
            expected.extend_bytes(binary.as_bytes());
            expected.push_byte(b':');
        }
        assert_eq!(write(&q, &mut s), expected);
    }
    #[test]
    fn empty_members_reproduce_seek_over_prefix_colon() {
        let q = graph(vec![SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::Data,
        )]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert_eq!(write(&q, &mut s), PropertyText::from("SgD::::::"));
    }
    #[test]
    fn empty_vector_and_empty_vector_elements_keep_distinct_framing() {
        for (values, expected) in [
            (vec![], "SgD:0::::::"),
            (
                vec![PropertyText::new(), PropertyText::new()],
                "SgD:0::,::::",
            ),
        ] {
            let q = graph(vec![prop(
                group(0, SubstanceGroupKind::Data),
                "DATAFIELDS",
                PropertyValue::StringVector(values),
            )]);
            let mut s = QueryCxSgroupState::from_query(&q);
            assert_eq!(write(&q, &mut s), PropertyText::from(expected));
        }
    }
    #[test]
    fn wrong_datafields_tag_keeps_group_write_without_final_counter_commit() {
        let q = graph(vec![prop(
            group(0, SubstanceGroupKind::Data),
            "DATAFIELDS",
            PropertyValue::String("wrong".into()),
        )]);
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(7))
            .unwrap();
        assert!(matches!(
            write_query_cx_data_sgroups(&q, &mut s, &[AtomId::new(0)]),
            Err(SmartsWriteError::PropertyValue(_))
        ));
        assert_eq!(
            s.groups[0].props().get(b"_cxsmilesOutputIndex".as_slice()),
            Some(&PropertyValue::UInt(7))
        );
        assert_eq!(counter(&s), Some(&PropertyValue::UInt(7)));
    }
    #[test]
    fn unsigned_counter_wraps_after_actual_group_write() {
        let q = graph(vec![group(0, SubstanceGroupKind::Data)]);
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(u32::MAX))
            .unwrap();
        write(&q, &mut s);
        assert_eq!(
            s.groups[0].props().get(b"_cxsmilesOutputIndex".as_slice()),
            Some(&PropertyValue::UInt(u32::MAX))
        );
        assert_eq!(counter(&s), Some(&PropertyValue::UInt(0)));
    }
    #[test]
    fn non_data_groups_keep_existing_counter() {
        let q = graph(vec![group(0, SubstanceGroupKind::StructuralRepeatUnit)]);
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(44))
            .unwrap();
        assert!(write(&q, &mut s).is_empty());
        assert_eq!(counter(&s), Some(&PropertyValue::UInt(44)));
    }
    #[test]
    fn invalid_reverse_order_fails_before_group_counter_write() {
        let q = graph(vec![group(0, SubstanceGroupKind::Data)]);
        let mut s = QueryCxSgroupState::from_query(&q);
        assert!(matches!(
            write_query_cx_data_sgroups(&q, &mut s, &[AtomId::new(999)]),
            Err(SmartsWriteError::CxSourceAtomOutOfRange { .. })
        ));
        assert!(
            !s.groups[0]
                .props()
                .contains_key(b"_cxsmilesOutputIndex".as_slice())
        );
    }
    #[test]
    fn data_then_polymer_share_actual_counter_and_hierarchy_indexes() {
        let q = graph(vec![
            group(0, SubstanceGroupKind::Data),
            group(1, SubstanceGroupKind::StructuralRepeatUnit)
                .with_parent(SubstanceGroupId::new(0)),
        ]);
        let mut s = QueryCxSgroupState::from_query(&q);
        s.properties
            .set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(9))
            .unwrap();
        write(&q, &mut s);
        write_query_cx_polymer_sgroups(&q, &mut s, &[AtomId::new(0)], &[]).unwrap();
        assert_eq!(counter(&s), Some(&PropertyValue::UInt(11)));
        assert_eq!(write_query_cx_sgroup_hierarchy(&mut s).unwrap(), "SgH:9:10");
    }
    #[test]
    fn enclosing_writer_emits_data_polymer_hierarchy_without_mutating_query() {
        let mut q = graph(vec![
            group(0, SubstanceGroupKind::Data),
            group(1, SubstanceGroupKind::StructuralRepeatUnit)
                .with_parent(SubstanceGroupId::new(0)),
        ]);
        let before = q.clone();
        let result = write_query_cx_extensions(
            &mut q,
            &(0..4).map(AtomId::new).collect::<Vec<_>>(),
            &(0..3).map(BondId::new).collect::<Vec<_>>(),
            cosmolkit_smiles::CxSmilesFields::SGROUPS | cosmolkit_smiles::CxSmilesFields::POLYMER,
        )
        .unwrap();
        assert_eq!(
            result,
            PropertyText::from("|SgD:0::::::,Sg:n:0:::::,SgH:0:1|")
        );
        assert_eq!(q, before);
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod data_sgroup_empty_presence_source_tests {
    use super::*;
    use cosmolkit_model::{
        AtomSpec, PropertyValue, SGroupData, SubstanceGroupId, replace_query_substance_groups,
    };
    use cosmolkit_types::Element;
    #[test]
    fn present_empty_vector_does_not_select_another_nonempty_value() {
        for raw in [false, true] {
            let mut group = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(0)])
                .with_data_field("other")
                .with_data(SGroupData::default());
            if raw {
                group
                    .set_prop("DATAFIELDS", PropertyValue::StringVector(vec![]))
                    .unwrap();
            }
            let mut q = QueryGraph::from_parts(
                vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
                vec![],
                [],
                vec![],
                vec![],
                vec![],
            )
            .unwrap();
            replace_query_substance_groups(&mut q, vec![group]).unwrap();
            let mut state = QueryCxSgroupState::from_query(&q);
            assert_eq!(
                write_query_cx_data_sgroups(&q, &mut state, &[AtomId::new(0)]).unwrap(),
                PropertyText::from("SgD:0::::::")
            );
        }
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod quote_string_source_tests {
    use super::*;
    #[test]
    fn empty_string_remains_empty() {
        assert!(quote_query_cx_string(b"").is_empty());
    }
    #[test]
    fn source_fix_does_not_escape_cx_delimiters() {
        assert_eq!(
            quote_query_cx_string(b"a.b;c:d,$|\\\"<>&"),
            PropertyText::from_bytes(b"a.b;c:d,$|\\\"<>&")
        );
    }
    #[test]
    fn native_counted_copy_preserves_every_byte() {
        let bytes = (0..=255).collect::<Vec<u8>>();
        assert_eq!(quote_query_cx_string(&bytes).as_bytes(), bytes.as_slice());
    }
    #[test]
    fn returned_value_is_independent_of_source_storage() {
        let mut bytes = vec![b'x'; 128];
        bytes[1] = 0;
        let output = quote_query_cx_string(&bytes);
        bytes.fill(b'y');
        assert_eq!(output.as_bytes()[1], 0);
        assert_eq!(output.as_bytes()[127], b'x');
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod atom_labels_source_tests {
    use super::*;
    use cosmolkit_model::{PropertyValue, QueryAtomIdentity};
    fn atom(i: usize, number: u8, props: Vec<(&str, PropertyValue)>) -> QueryAtom {
        let mut a = QueryAtom::from_identity_parts(
            AtomId::new(i),
            QueryAtomIdentity::from_atomic_number(number),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        );
        for (k, v) in props {
            a.set_prop(k, v).unwrap();
        }
        a
    }
    fn graph(atoms: Vec<QueryAtom>) -> QueryGraph {
        QueryGraph::from_parts(atoms, vec![], [], vec![], vec![], vec![]).unwrap()
    }
    fn text(t: &str) -> PropertyValue {
        PropertyValue::String(t.into())
    }
    fn write(q: &QueryGraph) -> PropertyText {
        write_query_cx_atom_labels(q, &(0..q.num_atoms()).map(AtomId::new).collect::<Vec<_>>())
            .unwrap()
    }
    #[test]
    fn empty_order_returns_empty() {
        let q = graph(vec![atom(0, 0, vec![("_fromAttchpt", text("bad"))])]);
        assert!(write_query_cx_atom_labels(&q, &[]).unwrap().is_empty());
    }
    #[test]
    fn generic_label_presence_wins_even_when_empty_and_keeps_binary() {
        let binary = PropertyText::from(vec![255, 0, b';']);
        let q = graph(vec![
            atom(
                0,
                6,
                vec![
                    ("_QueryAtomGenericLabel", text("")),
                    ("atomLabel", text("other")),
                ],
            ),
            atom(
                1,
                0,
                vec![
                    (
                        "_QueryAtomGenericLabel",
                        PropertyValue::String(binary.clone()),
                    ),
                    ("_fromAttchpt", text("bad")),
                ],
            ),
        ]);
        let before = q.clone();
        let mut expected = PropertyText::from("_p;");
        expected.extend_bytes(binary.as_bytes());
        expected.extend_bytes(b"_p");
        assert_eq!(write(&q), expected);
        assert_eq!(q, before);
    }
    #[test]
    fn only_source_pseudoatoms_and_zero_actual_identity_get_suffix() {
        let q = graph(vec![
            atom(0, 0, vec![("dummyLabel", text("Pol"))]),
            atom(1, 0, vec![("dummyLabel", text("Mod"))]),
            atom(
                2,
                6,
                vec![("dummyLabel", text("Pol")), ("atomLabel", text("C"))],
            ),
            atom(
                3,
                0,
                vec![("dummyLabel", text("Pol_p")), ("atomLabel", text("other"))],
            ),
        ]);
        assert_eq!(write(&q), PropertyText::from("Pol_p;Mod_p;C;other"));
    }
    #[test]
    fn unknown_dummy_label_falls_through_to_source_attach_points() {
        for n in [1, 2] {
            let q = graph(vec![atom(
                0,
                0,
                vec![
                    ("dummyLabel", text("other")),
                    ("_fromAttchpt", PropertyValue::Int(n)),
                    ("atomLabel", text("ordinary")),
                ],
            )]);
            assert_eq!(write(&q), PropertyText::from(format!("_AP{n}")));
        }
    }
    #[test]
    fn other_attach_point_values_reach_ordinary_label() {
        for n in [-1, 0, 3] {
            let q = graph(vec![atom(
                0,
                0,
                vec![
                    ("_fromAttchpt", PropertyValue::Int(n)),
                    ("atomLabel", text("ordinary")),
                ],
            )]);
            assert_eq!(write(&q), PropertyText::from("ordinary"));
        }
    }
    #[test]
    fn nonzero_identity_skips_bad_attach_point_conversion() {
        let q = graph(vec![atom(
            0,
            119,
            vec![
                ("_fromAttchpt", text("bad")),
                ("atomLabel", text("ordinary")),
            ],
        )]);
        assert_eq!(write(&q), PropertyText::from("ordinary"));
    }
    #[test]
    fn zero_identity_bad_attach_point_errors_before_ordinary_label() {
        let q = graph(vec![atom(
            0,
            0,
            vec![
                ("_fromAttchpt", text("bad")),
                ("atomLabel", text("ordinary")),
            ],
        )]);
        assert!(matches!(
            write_query_cx_atom_labels(&q, &[AtomId::new(0)]),
            Err(SmartsWriteError::CxAtomPropertyInt {
                property: "_fromAttchpt",
                ..
            })
        ));
    }
    #[test]
    fn recognized_pseudoatom_short_circuits_bad_attach_point() {
        let q = graph(vec![atom(
            0,
            0,
            vec![("dummyLabel", text("Mod")), ("_fromAttchpt", text("bad"))],
        )]);
        assert_eq!(write(&q), PropertyText::from("Mod_p"));
    }
    #[test]
    fn duplicate_first_atom_uses_source_identity_delimiter_comparison() {
        let q = graph(vec![
            atom(0, 6, vec![("atomLabel", text("A"))]),
            atom(1, 6, vec![("atomLabel", text("B"))]),
            atom(2, 6, vec![("atomLabel", text("C"))]),
        ]);
        assert_eq!(
            write_query_cx_atom_labels(
                &q,
                &[
                    AtomId::new(0),
                    AtomId::new(1),
                    AtomId::new(0),
                    AtomId::new(2)
                ]
            )
            .unwrap(),
            PropertyText::from("A;BA;C")
        );
    }
    #[test]
    fn semicolons_only_clear_but_empty_positions_and_delimiters_are_preserved() {
        let q = graph(vec![
            atom(0, 6, vec![("atomLabel", text(";;"))]),
            atom(1, 6, vec![]),
        ]);
        assert!(write(&q).is_empty());
        let q = graph(vec![
            atom(0, 6, vec![]),
            atom(1, 6, vec![("atomLabel", text("a.b:c,$|"))]),
            atom(2, 6, vec![]),
        ]);
        assert_eq!(write(&q), PropertyText::from(";a.b:c,$|;"));
    }
    #[test]
    fn source_atom_access_is_structurally_checked() {
        let q = graph(vec![]);
        assert!(matches!(
            write_query_cx_atom_labels(&q, &[AtomId::new(0)]),
            Err(SmartsWriteError::CxSourceAtomOutOfRange { atom_count: 0, .. })
        ));
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod atom_values_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue};
    use cosmolkit_types::Element;
    fn atom(i: usize, props: Vec<(&[u8], PropertyValue)>) -> QueryAtom {
        let mut a = QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C));
        for (k, v) in props {
            a.set_prop(PropertyText::from_bytes(k), v).unwrap();
        }
        a
    }
    fn text(t: &str) -> PropertyValue {
        PropertyValue::String(t.into())
    }
    fn graph(atoms: Vec<QueryAtom>) -> QueryGraph {
        QueryGraph::from_parts(atoms, vec![], [], vec![], vec![], vec![]).unwrap()
    }
    fn write(q: &QueryGraph, key: &[u8]) -> PropertyText {
        write_query_cx_atom_values(
            q,
            &(0..q.num_atoms()).map(AtomId::new).collect::<Vec<_>>(),
            key,
        )
        .unwrap()
    }
    #[test]
    fn empty_order_does_not_access_atoms_or_properties() {
        let q = graph(vec![]);
        assert!(write(&q, b"key").is_empty());
    }
    #[test]
    fn supplied_property_is_not_hardcoded_to_molfilevalue() {
        let q = graph(vec![atom(
            0,
            vec![
                (b"custom", text("custom")),
                (b"molFileValue", text("other")),
            ],
        )]);
        assert_eq!(write(&q, b"custom"), PropertyText::from("custom"));
        assert_eq!(write(&q, b"molFileValue"), PropertyText::from("other"));
    }
    #[test]
    fn counted_binary_key_and_value_remain_independent() {
        let bytes = vec![255, 0, b';', b'.', b'$'];
        let q = graph(vec![atom(
            0,
            vec![(
                &[255, 0, b'k'],
                PropertyValue::String(PropertyText::from(bytes.clone())),
            )],
        )]);
        let before = q.clone();
        assert_eq!(write(&q, &[255, 0, b'k']).as_bytes(), bytes.as_slice());
        assert_eq!(q, before);
    }
    #[test]
    fn absent_and_empty_values_preserve_all_position_separators() {
        let q = graph(vec![
            atom(0, vec![]),
            atom(1, vec![(b"key", text(""))]),
            atom(2, vec![]),
        ]);
        assert_eq!(write(&q, b"key"), PropertyText::from(";;"));
    }
    #[test]
    fn repeated_first_atom_is_separated_by_position() {
        let q = graph(vec![
            atom(0, vec![(b"key", text("A"))]),
            atom(1, vec![(b"key", text("B"))]),
        ]);
        assert_eq!(
            write_query_cx_atom_values(
                &q,
                &[AtomId::new(0), AtomId::new(1), AtomId::new(0)],
                b"key"
            )
            .unwrap(),
            PropertyText::from("A;B;A")
        );
    }
    #[test]
    fn scalar_reads_reuse_canonical_source_conversion_and_unescaped_quote() {
        let q = graph(vec![
            atom(0, vec![(b"key", PropertyValue::Int(-7))]),
            atom(1, vec![(b"key", PropertyValue::UInt(u32::MAX))]),
            atom(2, vec![(b"key", text("a.b:c,$|"))]),
        ]);
        assert_eq!(
            write(&q, b"key"),
            PropertyText::from("-7;4294967295;a.b:c,$|")
        );
    }
    #[test]
    fn invalid_atom_order_is_a_structured_source_access_failure() {
        let q = graph(vec![atom(0, vec![])]);
        assert!(matches!(
            write_query_cx_atom_values(&q, &[AtomId::new(1)], b"key"),
            Err(SmartsWriteError::CxSourceAtomOutOfRange { atom_count: 1, .. })
        ));
    }
    #[test]
    fn enclosing_presence_preflight_emits_empty_values_only_when_present() {
        for present in [false, true] {
            let mut q = graph(vec![
                atom(
                    0,
                    if present {
                        vec![(b"molFileValue", text(""))]
                    } else {
                        vec![]
                    },
                ),
                atom(1, vec![]),
            ]);
            let output = write_query_cx_extensions(
                &mut q,
                &[AtomId::new(0), AtomId::new(1)],
                &[],
                cosmolkit_smiles::CxSmilesFields::MOLFILE_VALUES,
            )
            .unwrap();
            assert_eq!(
                output,
                PropertyText::from(if present { "|$_AV:;$|" } else { "" })
            );
        }
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod radicals_source_tests {
    use super::*;
    use cosmolkit_model::AtomSpec;
    use cosmolkit_types::Element;
    fn graph(counts: &[u8]) -> QueryGraph {
        QueryGraph::from_parts(
            counts
                .iter()
                .enumerate()
                .map(|(i, c)| {
                    QueryAtom::new(
                        AtomId::new(i),
                        AtomSpec::new(Element::C).with_radical_electrons(*c),
                    )
                })
                .collect(),
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn write(q: &QueryGraph, order: &[usize]) -> (String, Vec<u32>) {
        let mut warnings = vec![];
        let text = write_query_cx_radicals(
            q,
            &order.iter().copied().map(AtomId::new).collect::<Vec<_>>(),
            &mut |n| warnings.push(n),
        )
        .unwrap();
        (text, warnings)
    }
    #[test]
    fn zero_radicals_and_empty_order_produce_no_output_or_warnings() {
        let q = graph(&[0, 0]);
        assert_eq!(write(&q, &[0, 1]), (String::new(), vec![]));
        assert_eq!(write(&q, &[]), (String::new(), vec![]));
    }
    #[test]
    fn ordered_counts_group_output_positions_and_keep_trailing_comma() {
        let q = graph(&[3, 1, 2, 1, 0]);
        assert_eq!(
            write(&q, &[2, 0, 3, 1, 2]),
            ("^1:2,3,^2:0,4,^5:1,".into(), vec![])
        );
    }
    #[test]
    fn unsupported_counts_warn_once_per_group_and_keep_positions_without_header() {
        let q = graph(&[4, 1, 255, 4]);
        assert_eq!(
            write(&q, &[0, 1, 2, 3]),
            ("^1:1,0,3,2,".into(), vec![4, 255])
        );
    }
    #[test]
    fn repeated_atom_ids_keep_distinct_output_positions() {
        let q = graph(&[2]);
        assert_eq!(write(&q, &[0, 0]), ("^2:0,1,".into(), vec![]));
    }
    #[test]
    fn collection_error_occurs_before_warning_iteration() {
        let q = graph(&[4]);
        let mut warnings = vec![];
        assert!(matches!(
            write_query_cx_radicals(&q, &[AtomId::new(0), AtomId::new(99)], &mut |n| warnings
                .push(n)),
            Err(SmartsWriteError::CxSourceAtomOutOfRange { .. })
        ));
        assert!(warnings.is_empty());
    }
    #[test]
    fn actual_atom_facts_are_used_and_query_remains_unchanged() {
        let mut q = graph(&[3]);
        q.atoms_mut()[0].set_predicate(QueryNode::predicate(
            AtomQueryPredicate::NumRadicalElectrons(0),
        ));
        let before = q.clone();
        assert_eq!(write(&q, &[0]), ("^5:0,".into(), vec![]));
        assert_eq!(q, before);
    }
    #[test]
    fn enclosing_writer_removes_one_source_trailing_comma() {
        let mut q = graph(&[1, 2, 3]);
        assert_eq!(
            write_query_cx_extensions(
                &mut q,
                &[AtomId::new(0), AtomId::new(1), AtomId::new(2)],
                &[],
                cosmolkit_smiles::CxSmilesFields::RADICALS
            )
            .unwrap(),
            PropertyText::from("|^1:0,^2:1,^5:2|")
        );
    }
    #[test]
    fn enclosing_disabled_radicals_keeps_no_output() {
        let mut q = graph(&[4]);
        assert!(
            write_query_cx_extensions(
                &mut q,
                &[AtomId::new(0)],
                &[],
                cosmolkit_smiles::CxSmilesFields::NONE
            )
            .unwrap()
            .is_empty()
        );
    }
    #[test]
    fn source_invalid_atom_access_is_structural_even_with_flags_off() {
        let mut q = graph(&[]);
        assert!(matches!(
            write_query_cx_extensions(
                &mut q,
                &[AtomId::new(0)],
                &[],
                cosmolkit_smiles::CxSmilesFields::NONE
            ),
            Err(SmartsWriteError::CxSourceAtomOutOfRange { .. })
        ));
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod coordinates_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, Conformer2D, Conformer3D, CoordinateDimension};
    use cosmolkit_types::Element;
    fn graph(
        n: usize,
        two: Vec<Conformer2D>,
        three: Vec<Conformer3D>,
        order: Option<Vec<CoordinateDimension>>,
    ) -> QueryGraph {
        let mut q = QueryGraph::from_parts(
            (0..n)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            vec![],
            [],
            two,
            three,
            vec![],
        )
        .unwrap();
        q.set_source_conformer_order(order).unwrap();
        q
    }
    fn write(q: &QueryGraph) -> String {
        write_query_cx_coordinates(q, &(0..q.num_atoms()).map(AtomId::new).collect::<Vec<_>>())
            .unwrap()
    }
    #[test]
    fn missing_conformer_errors_before_even_empty_or_invalid_atom_order() {
        let q = graph(0, vec![], vec![], None);
        for ids in [vec![], vec![AtomId::new(99)]] {
            assert!(matches!(
                write_query_cx_coordinates(&q, &ids),
                Err(SmartsWriteError::CxMissingConformer)
            ));
        }
    }
    #[test]
    fn actual_conformer_with_empty_order_returns_empty_string() {
        let q = graph(1, vec![Conformer2D::new(0, vec![[1.0, 2.0]])], vec![], None);
        assert_eq!(write_query_cx_coordinates(&q, &[]).unwrap(), "");
    }
    #[test]
    fn two_d_coordinates_share_strict_zero_threshold_and_general_format() {
        let q = graph(
            2,
            vec![Conformer2D::new(0, vec![[1.25, 2.0], [-9e-5, 1e-4]])],
            vec![],
            None,
        );
        assert_eq!(write(&q), "1.25,2,;0,0.0001,");
    }
    #[test]
    fn three_column_storage_does_not_override_actual_is3d_flag() {
        for (flag, expected) in [(false, "1,2,"), (true, "1,2,7")] {
            let q = graph(
                1,
                vec![],
                vec![Conformer3D::new(0, vec![[1.0, 2.0, 7.0]], flag)],
                None,
            );
            assert_eq!(write(&q), expected);
        }
    }
    #[test]
    fn z_is_omitted_by_formatted_zero_string_only() {
        let q = graph(
            2,
            vec![],
            vec![Conformer3D::new(
                0,
                vec![[1.0, 2.0, -0.0], [3.0, 4.0, 1e-4]],
                true,
            )],
            None,
        );
        assert_eq!(write(&q), "1,2,;3,4,0.0001");
    }
    #[test]
    fn mixed_dimensions_follow_actual_front_not_id_or_three_d_preference() {
        for (order, expected) in [
            (
                vec![CoordinateDimension::TwoD, CoordinateDimension::ThreeD],
                "9,8,",
            ),
            (
                vec![CoordinateDimension::ThreeD, CoordinateDimension::TwoD],
                "1,2,3",
            ),
        ] {
            let q = graph(
                1,
                vec![Conformer2D::new(900, vec![[9.0, 8.0]])],
                vec![Conformer3D::new(1, vec![[1.0, 2.0, 3.0]], true)],
                Some(order),
            );
            assert_eq!(write(&q), expected);
        }
    }
    #[test]
    fn mixed_missing_source_order_remains_a_structural_error() {
        let q = graph(
            1,
            vec![Conformer2D::new(900, vec![[9.0, 8.0]])],
            vec![Conformer3D::new(1, vec![[1.0, 2.0, 3.0]], true)],
            None,
        );
        assert!(matches!(
            write_query_cx_coordinates(&q, &[AtomId::new(0)]),
            Err(SmartsWriteError::CxCoordinateSource(
                cosmolkit_model::CoordinateValidationError::MissingSourceConformerOrder
            ))
        ));
    }
    #[test]
    fn output_order_repeats_positions_without_mutating_input() {
        let q = graph(
            2,
            vec![Conformer2D::new(99, vec![[1.0, 2.0], [3.0, 4.0]])],
            vec![],
            None,
        );
        let before = q.clone();
        assert_eq!(
            write_query_cx_coordinates(&q, &[AtomId::new(1), AtomId::new(1), AtomId::new(0)])
                .unwrap(),
            "3,4,;3,4,;1,2,"
        );
        assert_eq!(q, before);
    }
    #[test]
    fn nonfinite_values_keep_source_case_sign_and_original_storage_bits() {
        let point = [
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0xfff8000000001234),
        ];
        let q = graph(
            1,
            vec![],
            vec![Conformer3D::new(0, vec![point], true)],
            None,
        );
        assert_eq!(write(&q), "inf,-inf,-nan");
        assert_eq!(
            q.conformers_3d()[0].coordinates()[0].map(f64::to_bits),
            point.map(f64::to_bits)
        );
    }
    #[test]
    fn enclosing_count_guard_owns_no_conformer_skip_and_empty_order_parentheses() {
        let mut q = graph(1, vec![], vec![], None);
        assert!(
            write_query_cx_extensions(
                &mut q,
                &[AtomId::new(0)],
                &[],
                cosmolkit_smiles::CxSmilesFields::COORDS
            )
            .unwrap()
            .is_empty()
        );
        let mut q = graph(1, vec![Conformer2D::new(0, vec![[1.0, 2.0]])], vec![], None);
        assert_eq!(
            write_query_cx_extensions(&mut q, &[], &[], cosmolkit_smiles::CxSmilesFields::COORDS)
                .unwrap(),
            PropertyText::from("|()|")
        );
        assert_eq!(
            write_query_cx_extensions(
                &mut q,
                &[AtomId::new(0)],
                &[],
                cosmolkit_smiles::CxSmilesFields::COORDS
            )
            .unwrap(),
            PropertyText::from("|(1,2,)|")
        );
    }
    #[test]
    fn coordinate_range_errors_propagate_from_canonical_owner() {
        let q = graph(1, vec![Conformer2D::new(0, vec![[1.0, 2.0]])], vec![], None);
        assert!(matches!(
            write_query_cx_coordinates(&q, &[AtomId::new(1)]),
            Err(SmartsWriteError::CxCoordinateOutput(
                cosmolkit_smiles::SmilesParseError::CxCoordinateAtomOutOfRange {
                    atom_count: 1,
                    ..
                }
            ))
        ));
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod quote_atomprop_consumer_source_tests {
    use super::*;
    use cosmolkit_model::AtomSpec;
    use cosmolkit_types::Element;
    #[test]
    fn actual_atom_properties_consumer_quotes_names_and_values_through_canonical_owner() {
        let mut atom = QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C));
        atom.set_prop("a.b", "c.d").unwrap();
        let query = QueryGraph::from_parts(vec![atom], vec![], [], vec![], vec![], vec![]).unwrap();
        let before = query.clone();
        assert_eq!(
            write_query_cx_atom_properties(&query, &[AtomId::new(0)]).unwrap(),
            PropertyText::from("atomProp:0.a&#46;b.c&#46;d")
        );
        assert_eq!(query, before);
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod atom_properties_consumer_source_tests {
    use super::*;
    use cosmolkit_model::AtomSpec;
    use cosmolkit_smiles::CxSmilesFields;
    use cosmolkit_types::Element;
    #[test]
    fn enclosing_cx_writer_uses_canonical_properties_and_preserves_structural_errors() {
        let mut atom = QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C));
        atom.set_prop("a.b", "c.d").unwrap();
        let mut query =
            QueryGraph::from_parts(vec![atom], vec![], [], vec![], vec![], vec![]).unwrap();
        let before = query.clone();
        assert_eq!(
            write_query_cx_extensions(
                &mut query,
                &[AtomId::new(0)],
                &[],
                CxSmilesFields::ATOM_PROPS
            )
            .unwrap(),
            PropertyText::from("|atomProp:0.a&#46;b.c&#46;d|")
        );
        assert_eq!(query, before);
        assert!(matches!(
            write_query_cx_atom_properties(&query, &[AtomId::new(1)]),
            Err(SmartsWriteError::CxSourceAtomOutOfRange { atom_count: 1, .. })
        ));
        let mut atom = QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C));
        atom.set_prop("__computedProps", "bad-tag").unwrap();
        let mut query =
            QueryGraph::from_parts(vec![atom], vec![], [], vec![], vec![], vec![]).unwrap();
        assert!(matches!(
            write_query_cx_atom_properties(&query, &[AtomId::new(0)]),
            Err(SmartsWriteError::CxPropertyList { .. })
        ));
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod bond_config_source_tests {
    use super::*;
    use cosmolkit_core::{AtropisomerWedgeAssignment, AtropisomerWedgeUpdate, WedgeAssignments};
    use cosmolkit_model::{
        AtomQueryPredicate, BondSpec, Conformer2D, QueryAtomIdentity, QueryBond, QueryNode,
    };
    fn graph(n: usize, edges: &[(usize, usize)], coords: bool) -> QueryGraph {
        QueryGraph::from_parts(
            (0..n)
                .map(|i| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            if coords {
                vec![Conformer2D::new(0, vec![[0.0, 0.0]; n])]
            } else {
                vec![]
            },
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn write(
        q: &QueryGraph,
        ao: &[usize],
        bo: &[usize],
        coords: bool,
        only: bool,
    ) -> Result<String, SmartsWriteError> {
        write_query_cx_bond_config(
            q,
            &ao.iter().map(|&i| AtomId::new(i)).collect::<Vec<_>>(),
            &bo.iter().map(|&i| BondId::new(i)).collect::<Vec<_>>(),
            coords,
            &WedgeAssignments::default(),
            only,
        )
    }
    fn axial() -> QueryGraph {
        let mut q = graph(6, &[(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)], false);
        q.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::AtropCw)
            .unwrap();
        q.bonds_mut()[1]
            .bond_mut()
            .set_direction(BondDirection::BeginWedge);
        q
    }
    #[test]
    fn empty_order_does_not_eagerly_read_missing_conformer() {
        assert_eq!(
            write(&graph(0, &[], false), &[], &[], true, false).unwrap(),
            ""
        );
    }
    #[test]
    fn explicit_direction_precedes_bad_cfg_and_unreached_conformer() {
        let mut q = graph(2, &[(0, 1)], false);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop("_MolFileBondCfg", "bad-uint")
            .unwrap();
        for (d, expected) in [
            (BondDirection::BeginWedge, "wU:0.0"),
            (BondDirection::BeginDash, "wD:0.0"),
            (BondDirection::Unknown, "w:0.0"),
        ] {
            q.bonds_mut()[0].bond_mut().set_direction(d);
            assert_eq!(write(&q, &[0, 1], &[0], true, false).unwrap(), expected);
        }
    }
    #[test]
    fn source_uint_cfg_requires_conformer_only_after_none_and_preserves_bad_cast() {
        let mut q = graph(2, &[(0, 1)], false);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop("_MolFileBondCfg", 0u32)
            .unwrap();
        assert!(matches!(
            write(&q, &[0, 1], &[0], true, false),
            Err(SmartsWriteError::CxMissingConformer)
        ));
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop("_MolFileBondCfg", -1i32)
            .unwrap();
        assert!(matches!(
            write(&q, &[0, 1], &[0], false, false),
            Err(SmartsWriteError::CxBondPropertyUInt { .. })
        ));
    }
    #[test]
    fn can_have_direction_gate_precedes_cfg_read_and_keeps_aromatic() {
        let mut q = graph(2, &[(0, 1)], false);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop("_MolFileBondCfg", "bad")
            .unwrap();
        q.bonds_mut()[0].bond_mut().set_order(BondOrder::Double);
        assert_eq!(write(&q, &[0, 1], &[0], true, false).unwrap(), "");
        q.bonds_mut()[0].bond_mut().set_order(BondOrder::Aromatic);
        q.bonds_mut()[0]
            .bond_mut()
            .set_direction(BondDirection::Unknown);
        assert_eq!(write(&q, &[0, 1], &[0], false, false).unwrap(), "w:0.0");
    }
    #[test]
    fn first_duplicate_and_missing_atom_use_native_find_end_distance() {
        let mut q = graph(2, &[(0, 1)], false);
        q.bonds_mut()[0]
            .bond_mut()
            .set_direction(BondDirection::Unknown);
        assert_eq!(
            write(&q, &[1, 0, 0], &[0, 0], false, false).unwrap(),
            "w:1.0,1.1"
        );
        assert_eq!(write(&q, &[1], &[0], false, false).unwrap(), "w:1.0");
    }
    #[test]
    fn lexical_group_order_and_output_bond_positions_are_source_defined() {
        let mut q = graph(4, &[(0, 1), (1, 2), (2, 3)], false);
        for (i, d) in [
            BondDirection::BeginWedge,
            BondDirection::Unknown,
            BondDirection::BeginDash,
        ]
        .into_iter()
        .enumerate()
        {
            q.bonds_mut()[i].bond_mut().set_direction(d);
        }
        assert_eq!(
            write(&q, &[0, 1, 2, 3], &[0, 2, 1, 0], true, false).unwrap(),
            "w:1.2,wD:2.1,wU:0.0,0.3"
        );
    }
    #[test]
    fn atrop_only_early_and_final_gates_skip_cfg_and_non_atrop_unknown() {
        let mut q = graph(2, &[(0, 1)], false);
        q.bonds_mut()[0]
            .bond_mut()
            .set_prop("_MolFileBondCfg", "bad")
            .unwrap();
        assert_eq!(write(&q, &[0, 1], &[0], true, true).unwrap(), "");
        q.bonds_mut()[0]
            .bond_mut()
            .set_direction(BondDirection::Unknown);
        assert_eq!(write(&q, &[0, 1], &[0], false, true).unwrap(), "");
    }
    #[test]
    fn missing_atrop_carrier_errors_only_when_coordinate_free_flip_is_reached() {
        let mut q = graph(3, &[(0, 1), (0, 2)], false);
        q.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::AtropCw)
            .unwrap();
        q.bonds_mut()[1]
            .bond_mut()
            .set_direction(BondDirection::BeginWedge);
        assert!(matches!(
            write(&q, &[0, 1, 2], &[1], false, false),
            Err(SmartsWriteError::CxBondConfigAtropMissingCarriers { .. })
        ));
        assert_eq!(write(&q, &[0, 1, 2], &[1], true, false).unwrap(), "wU:0.0");
    }
    #[test]
    fn atrop_flip_uses_axial_order_and_only_opposite_side_carrier_order() {
        let q = axial();
        let before = q.clone();
        for (order, expected) in [
            ([0, 1, 2, 3, 4, 5], "wU:0.0"),
            ([1, 0, 2, 3, 4, 5], "wD:1.0"),
            ([0, 1, 2, 3, 5, 4], "wD:0.0"),
            ([0, 1, 3, 2, 4, 5], "wU:0.0"),
            ([1, 0, 2, 3, 5, 4], "wU:1.0"),
        ] {
            assert_eq!(write(&q, &order, &[1], false, true).unwrap(), expected);
        }
        assert_eq!(q, before);
    }
    #[test]
    fn first_axial_neighbor_wins_without_probing_later_invalid_carriers() {
        let mut q = graph(7, &[(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (0, 6)], false);
        for i in [0, 5] {
            q.bonds_mut()[i]
                .bond_mut()
                .set_stereo(BondStereo::AtropCw)
                .unwrap();
        }
        q.bonds_mut()[1]
            .bond_mut()
            .set_direction(BondDirection::BeginWedge);
        assert_eq!(
            write(&q, &[0, 1, 2, 3, 4, 5, 6], &[1], false, true).unwrap(),
            "wU:0.0"
        );
    }
    #[test]
    fn cfg_wedge_does_not_repeat_the_earlier_atrop_detection() {
        let mut q = axial();
        q.bonds_mut()[1]
            .bond_mut()
            .set_direction(BondDirection::None);
        q.bonds_mut()[1]
            .bond_mut()
            .set_prop("_MolFileBondCfg", 1u32)
            .unwrap();
        assert_eq!(
            write(&q, &[0, 1, 2, 3, 4, 5], &[1], false, false).unwrap(),
            ""
        );
        assert_eq!(
            write(&q, &[0, 1, 2, 3, 4, 5], &[1], true, false).unwrap(),
            "wU:0.0"
        );
    }
    #[test]
    fn coordinate_molfile_child_uses_literal_code_three_switch() {
        let mut q = graph(2, &[(0, 1)], true);
        q.bonds_mut()[0]
            .bond_mut()
            .set_direction(BondDirection::EitherDouble);
        assert_eq!(write(&q, &[0, 1], &[0], true, false).unwrap(), "w:0.0");
        q.bonds_mut()[0]
            .bond_mut()
            .set_direction(BondDirection::None);
        assert_eq!(write(&q, &[0, 1], &[0], true, false).unwrap(), "");
        assert!(matches!(
            write(&q, &[0, 1], &[1], true, false),
            Err(SmartsWriteError::CxSourceBondOutOfRange { bond_count: 1, .. })
        ));
    }
    #[test]
    fn shared_atrop_map_retains_source_mutated_endpoint_and_direction_without_graph_clone() {
        let mut q = graph(6, &[(0, 1), (2, 0), (0, 3), (1, 4), (1, 5)], false);
        q.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::AtropCw)
            .unwrap();
        let w = WedgeAssignments::from_atropisomer_wedge_assignment(AtropisomerWedgeAssignment {
            source_map_writes: vec![BondId::new(1)],
            bond_updates: vec![AtropisomerWedgeUpdate {
                bond: BondId::new(1),
                begin: AtomId::new(0),
                end: AtomId::new(2),
                direction: BondDirection::BeginDash,
                atropisomer_bond: BondId::new(0),
            }],
            diagnostics: vec![],
        });
        let before = q.clone();
        assert_eq!(
            write_query_cx_bond_config(
                &q,
                &(0..6).map(AtomId::new).collect::<Vec<_>>(),
                &[BondId::new(1)],
                false,
                &w,
                true
            )
            .unwrap(),
            "wD:0.0"
        );
        assert_eq!(q, before);
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod typed_bonds_consumer_source_tests {
    use super::*;
    use cosmolkit_model::{AtomQueryPredicate, BondSpec, QueryAtomIdentity, QueryBond, QueryNode};
    use cosmolkit_smiles::CxSmilesFields as F;
    fn graph() -> QueryGraph {
        QueryGraph::from_parts(
            (0..4)
                .map(|i| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect(),
            [BondOrder::DativeOne, BondOrder::Dative, BondOrder::Hydrogen]
                .into_iter()
                .enumerate()
                .map(|(i, o)| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(i), AtomId::new(i + 1), o),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn actual_enclosing_writer_emits_only_exact_coordinate_and_hydrogen_types() {
        let mut q = graph();
        let before = q.clone();
        let ao = (0..4).map(AtomId::new).collect::<Vec<_>>();
        assert_eq!(
            write_query_cx_extensions(
                &mut q,
                &ao,
                &[BondId::new(2), BondId::new(0), BondId::new(1)],
                F::COORDINATE_BONDS | F::HYDROGEN_BONDS
            )
            .unwrap(),
            PropertyText::from("|C:1.2,H:2.0|")
        );
        assert_eq!(q, before);
    }
    #[test]
    fn query_adapter_preserves_first_find_end_and_structured_range_error() {
        let q = graph();
        assert_eq!(
            write_query_cx_typed_bonds(
                &q,
                &[AtomId::new(1), AtomId::new(1)],
                &[BondId::new(1)],
                BondOrder::Dative,
                "C"
            )
            .unwrap(),
            PropertyText::from("C:0.0")
        );
        assert_eq!(
            write_query_cx_typed_bonds(&q, &[], &[BondId::new(1)], BondOrder::Dative, "C").unwrap(),
            PropertyText::from("C:0.0")
        );
        assert!(matches!(
            write_query_cx_typed_bonds(&q, &[], &[BondId::new(3)], BondOrder::Dative, "C"),
            Err(SmartsWriteError::CxSourceBondOutOfRange { bond_count: 3, .. })
        ));
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod zero_bonds_consumer_source_tests {
    use super::*;
    use cosmolkit_model::{AtomQueryPredicate, BondSpec, QueryAtomIdentity, QueryBond, QueryNode};
    fn graph() -> QueryGraph {
        QueryGraph::from_parts(
            (0..4)
                .map(|i| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect(),
            [BondOrder::Zero, BondOrder::Single, BondOrder::Zero]
                .into_iter()
                .enumerate()
                .map(|(i, o)| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(i), AtomId::new(i + 1), o),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn enclosing_cx_writer_uses_source_output_positions_and_shared_owner() {
        let mut q = graph();
        let before = q.clone();
        assert_eq!(
            write_query_cx_extensions(
                &mut q,
                &(0..4).map(AtomId::new).collect::<Vec<_>>(),
                &[BondId::new(1), BondId::new(2), BondId::new(0)],
                cosmolkit_smiles::CxSmilesFields::ZERO_BONDS
            )
            .unwrap(),
            PropertyText::from("|Z:1,2|")
        );
        assert_eq!(q, before);
    }
    #[test]
    fn query_adapter_preserves_range_error_instead_of_empty_output() {
        assert!(matches!(
            write_query_cx_zero_bonds(&graph(), &[BondId::new(3)]),
            Err(SmartsWriteError::CxSourceBondOutOfRange { bond_count: 3, .. })
        ));
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod ring_bond_stereo_source_tests {
    use super::*;
    use cosmolkit_core::{RingFindType, RingInfo};
    use cosmolkit_model::{AtomQueryPredicate, BondSpec, QueryAtomIdentity, QueryBond, QueryNode};
    fn graph(n: usize, edges: &[(usize, usize)]) -> QueryGraph {
        QueryGraph::from_parts(
            (0..n)
                .map(|i| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn cycle(n: usize, extra: &[(usize, usize)]) -> QueryGraph {
        let mut edges = (0..n).map(|i| (i, (i + 1) % n)).collect::<Vec<_>>();
        edges.extend_from_slice(extra);
        let count = edges.iter().map(|&(a, b)| a.max(b)).max().unwrap() + 1;
        graph(count, &edges)
    }
    fn cache(q: &mut QueryGraph, atom_rows: &[Vec<usize>], bond_rows: &[Vec<usize>]) {
        let mut r = RingInfo::new(RingFindType::Sssr, q.num_atoms(), q.num_bonds());
        for (a, b) in atom_rows.iter().zip(bond_rows) {
            r.add_ring(a, b).unwrap();
        }
        q.replace_source_ring_info(r.into_source_snapshot());
    }
    fn activate(q: &mut QueryGraph, n: usize) {
        cache(q, &[(0..n).collect()], &[(0..n).collect()]);
    }
    fn tag(q: &mut QueryGraph, id: usize, stereo: BondStereo, refs: [usize; 2]) {
        let b = q.bonds_mut()[id].bond_mut();
        b.set_order(BondOrder::Double);
        b.set_stereo_atoms(Some(refs.map(AtomId::new)));
        b.set_stereo(stereo).unwrap();
    }
    fn write(q: &QueryGraph, ao: &[usize], bo: &[usize]) -> Result<String, SmartsWriteError> {
        write_query_cx_ring_bond_stereo(
            q,
            &ao.iter().map(|&i| AtomId::new(i)).collect::<Vec<_>>(),
            &bo.iter().map(|&i| BondId::new(i)).collect::<Vec<_>>(),
        )
    }
    #[test]
    fn uninitialized_source_ring_cache_returns_before_orders_and_never_finds_cycle() {
        let mut q = cycle(8, &[]);
        tag(&mut q, 0, BondStereo::Any, [7, 2]);
        assert_eq!(write(&q, &[], &[99]).unwrap(), "");
        assert_eq!(write(&q, &(0..8).collect::<Vec<_>>(), &[0]).unwrap(), "");
        let mut source = q.source_ring_info().clone();
        source.find_type = 17;
        q.replace_source_ring_info(source);
        assert_eq!(write(&q, &[], &[0]).unwrap(), "");
    }
    #[test]
    fn initialized_empty_actual_cache_is_not_rebuilt_from_cyclic_graph() {
        let mut q = cycle(8, &[]);
        tag(&mut q, 0, BondStereo::Any, [7, 2]);
        cache(&mut q, &[], &[]);
        assert!(q.source_ring_info().initialized);
        assert_eq!(write(&q, &[], &[0]).unwrap(), "");
    }
    #[test]
    fn native_minimum_ring_size_gate_distinguishes_seven_and_eight() {
        for (n, expected) in [(7, ""), (8, "ctu:0")] {
            let mut q = cycle(n, &[]);
            tag(&mut q, 0, BondStereo::Any, [n - 1, 2]);
            activate(&mut q, n);
            assert_eq!(write(&q, &[], &[0]).unwrap(), expected);
        }
    }
    #[test]
    fn cache_membership_gate_precedes_get_bond_and_preserves_reached_range_error() {
        let mut q = cycle(8, &[]);
        activate(&mut q, 8);
        assert_eq!(write(&q, &[], &[99]).unwrap(), "");
        cache(&mut q, &[(0..8).collect()], &[vec![99; 8]]);
        assert!(matches!(
            write(&q, &[], &[99]),
            Err(SmartsWriteError::CxSourceBondOutOfRange { bond_count: 8, .. })
        ));
    }
    #[test]
    fn only_double_or_aromatic_with_any_cis_trans_tags_are_included() {
        let mut q = cycle(8, &[]);
        activate(&mut q, 8);
        tag(&mut q, 0, BondStereo::Any, [7, 2]);
        for (o, expected) in [
            (BondOrder::Single, ""),
            (BondOrder::Double, "ctu:0"),
            (BondOrder::Aromatic, "ctu:0"),
        ] {
            q.bonds_mut()[0].bond_mut().set_order(o);
            assert_eq!(write(&q, &[], &[0]).unwrap(), expected);
        }
        q.bonds_mut()[0].bond_mut().set_order(BondOrder::Double);
        for s in [BondStereo::None, BondStereo::E, BondStereo::Z] {
            q.bonds_mut()[0].bond_mut().set_stereo(s).unwrap();
            assert_eq!(write(&q, &[], &[0]).unwrap(), "");
        }
    }
    #[test]
    fn output_positions_include_unselected_rows_and_repeat_bond_ids() {
        let mut q = cycle(8, &[]);
        tag(&mut q, 0, BondStereo::Any, [7, 2]);
        activate(&mut q, 8);
        assert_eq!(write(&q, &[], &[1, 0, 0]).unwrap(), "ctu:1,2");
    }
    #[test]
    fn three_group_return_is_literal_concatenation_without_intergroup_commas() {
        let edges = (0..3)
            .flat_map(|r| (0..8).map(move |i| (r * 8 + i, r * 8 + (i + 1) % 8)))
            .collect::<Vec<_>>();
        let mut q = graph(24, &edges);
        cache(
            &mut q,
            &[(0..8).collect(), (8..16).collect(), (16..24).collect()],
            &[(0..8).collect(), (8..16).collect(), (16..24).collect()],
        );
        tag(&mut q, 0, BondStereo::Cis, [7, 2]);
        tag(&mut q, 8, BondStereo::Trans, [15, 10]);
        tag(&mut q, 16, BondStereo::Any, [23, 18]);
        assert_eq!(
            write(&q, &(0..24).collect::<Vec<_>>(), &[0, 8, 16]).unwrap(),
            "c:0t:1ctu:2"
        );
    }
    #[test]
    fn stereo_any_never_reads_atom_order_or_stereo_references() {
        let mut q = cycle(8, &[(0, 8)]);
        tag(&mut q, 0, BondStereo::Any, [7, 2]);
        q.bonds_mut()[0].bond_mut().set_stereo_atoms(None);
        activate(&mut q, 8);
        assert_eq!(write(&q, &[], &[0]).unwrap(), "ctu:0");
    }
    #[test]
    fn false_boolean_subscript_uses_atom_order_zero_and_cis_is_not_xor_flipped() {
        let mut q = cycle(8, &[(0, 8)]);
        tag(&mut q, 0, BondStereo::Trans, [7, 2]);
        activate(&mut q, 8);
        assert_eq!(write(&q, &(0..9).collect::<Vec<_>>(), &[0]).unwrap(), "t:0");
        let ao = [1, 0, 2, 3, 4, 5, 6, 7, 8];
        assert_eq!(write(&q, &ao, &[0]).unwrap(), "c:0");
        q.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::Cis)
            .unwrap();
        assert_eq!(write(&q, &ao, &[0]).unwrap(), "c:0");
    }
    #[test]
    fn true_boolean_subscript_uses_atom_order_one_as_integer_truth() {
        let edges = (1..=8)
            .map(|i| (i, if i == 8 { 1 } else { i + 1 }))
            .chain([(1, 0)])
            .collect::<Vec<_>>();
        let mut q = graph(9, &edges);
        cache(&mut q, &[(1..=8).collect()], &[(0..8).collect()]);
        tag(&mut q, 0, BondStereo::Trans, [8, 3]);
        assert_eq!(write(&q, &(0..9).collect::<Vec<_>>(), &[0]).unwrap(), "c:0");
        assert_eq!(
            write(&q, &[1, 0, 2, 3, 4, 5, 6, 7, 8], &[0]).unwrap(),
            "t:0"
        );
    }
    #[test]
    fn end_atom_branch_and_both_endpoint_toggles_follow_native_parity() {
        for (extra, ao, expected) in [
            (vec![(1, 8)], vec![1, 0, 2, 3, 4, 5, 6, 7, 8], "c:0"),
            (
                vec![(1, 8), (0, 9)],
                vec![1, 0, 2, 3, 4, 5, 6, 7, 8, 9],
                "t:0",
            ),
        ] {
            let mut q = cycle(8, &extra);
            tag(&mut q, 0, BondStereo::Trans, [7, 2]);
            activate(&mut q, 8);
            assert_eq!(write(&q, &ao, &[0]).unwrap(), expected);
        }
    }
    #[test]
    fn reached_missing_reference_and_atom_order_errors_are_structural_and_immutable() {
        let mut q = cycle(8, &[(0, 8)]);
        tag(&mut q, 0, BondStereo::Trans, [7, 2]);
        activate(&mut q, 8);
        let before = q.clone();
        assert!(matches!(
            write(&q, &[], &[0]),
            Err(SmartsWriteError::CxRingAtomOrderIndex { index: 7, count: 0 })
        ));
        assert_eq!(q, before);
        q.bonds_mut()[0].bond_mut().set_stereo_atoms(None);
        let before = q.clone();
        assert!(matches!(
            write(&q, &(0..9).collect::<Vec<_>>(), &[0]),
            Err(SmartsWriteError::CxRingStereoReferenceMissing { side: 0, .. })
        ));
        assert_eq!(q, before);
    }

    #[test]
    fn enclosing_cx_writer_promotes_rings_before_ring_stereo_output() {
        let mut q = cycle(8, &[]);
        tag(&mut q, 0, BondStereo::Any, [7, 2]);
        let ao = (0..8).map(AtomId::new).collect::<Vec<_>>();
        let bo = [BondId::new(0)];
        let before = q.clone();
        assert!(
            write_query_cx_extensions(&mut q, &ao, &bo, cosmolkit_smiles::CxSmilesFields::NONE)
                .unwrap()
                .is_empty()
        );
        assert_eq!(q, before);
        // getCXExtensions(CX_BOND_CFG) calls pickBondsToWedge first. Its
        // atrop stage calls findSSSR even without an atrop bond; only the
        // standalone get_ringbond_cistrans_block skips an uninitialized cache.
        assert_eq!(
            write_query_cx_extensions(&mut q, &ao, &bo, cosmolkit_smiles::CxSmilesFields::BOND_CFG)
                .unwrap(),
            PropertyText::from("|ctu:0|")
        );
        let source = q.source_ring_info();
        assert!(source.initialized);
        assert_eq!(source.atom_rings.len(), 1);
        assert_eq!(source.bond_rings.len(), 1);
        assert_eq!(
            source.atom_rings[0]
                .iter()
                .copied()
                .collect::<BTreeSet<_>>(),
            (0..8).map(AtomId::new).collect()
        );
        assert_eq!(
            source.bond_rings[0]
                .iter()
                .copied()
                .collect::<BTreeSet<_>>(),
            (0..8).map(BondId::new).collect()
        );
        let mut expected = before;
        expected.replace_source_ring_info(source.clone());
        assert_eq!(q, expected, "no non-cache state is changed");
        activate(&mut q, 8);
        let before = q.clone();
        assert_eq!(
            write_query_cx_extensions(&mut q, &ao, &bo, cosmolkit_smiles::CxSmilesFields::BOND_CFG)
                .unwrap(),
            PropertyText::from("|ctu:0|")
        );
        assert_eq!(q, before);
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod link_nodes_consumer_source_tests {
    use super::*;
    use cosmolkit_model::{
        AtomQueryPredicate, BondSpec, PropertyValue, QueryAtomIdentity, QueryBond, QueryNode,
    };
    fn graph(raw: &[u8]) -> QueryGraph {
        QueryGraph::from_parts(
            (0..4)
                .map(|i| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect(),
            [(1, 0), (1, 2), (1, 3)]
                .into_iter()
                .enumerate()
                .map(|(i, (a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [(
                PropertyText::from("_molLinkNodes"),
                PropertyValue::String(PropertyText::from_bytes(raw)),
            )],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn enclosing_cx_uses_full_source_link_parser_and_forward_indexes() {
        let mut q = graph(b"1 3 2 2 1 2 3");
        let before = q.clone();
        assert_eq!(
            write_query_cx_extensions(
                &mut q,
                &[3, 0, 2, 1].map(AtomId::new),
                &[0, 1, 2].map(BondId::new),
                cosmolkit_smiles::CxSmilesFields::LINKNODES
            )
            .unwrap(),
            PropertyText::from("|LN:0:1.3.3.2|")
        );
        assert_eq!(q, before);
    }
    #[test]
    fn invalid_counts_and_three_bonds_are_rejected_before_valid_node() {
        let q = graph(b"0 3 2 2 1 2 3|1 3 3 2 1 2 3 2 4|1 3 2 2 1 2 3");
        assert_eq!(
            write_query_cx_link_nodes(&q, &[0, 1, 2, 3].map(AtomId::new)).unwrap(),
            "LN:1:1.3.0.2"
        );
    }
    #[test]
    fn reached_source_range_errors_propagate_through_search() {
        assert!(matches!(
            write_query_cx_link_nodes(&graph(b"1 3 2 2 0 2 3"), &[]),
            Err(SmartsWriteError::CanonicalTraversal(
                cosmolkit_smiles::SmilesParseError::Cx(_)
            ))
        ));
        assert!(matches!(
            write_query_cx_link_nodes(&graph(b"1 3 2 2 1 2 3"), &[AtomId::new(99)]),
            Err(SmartsWriteError::CanonicalTraversal(
                cosmolkit_smiles::SmilesParseError::CxLinkAtomOutOfRange { .. }
            ))
        ));
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod append_cx_source_tests {
    use super::*;
    #[test]
    fn empty_addition_never_changes_any_base() {
        for b in [b"".as_slice(), b"|", b"||", b"abc", b"\xff\0"] {
            let mut out = PropertyText::from_bytes(b);
            append_query_cx_extension(b"", &mut out);
            assert_eq!(out.as_bytes(), b);
        }
    }
    #[test]
    fn separator_depends_only_on_base_byte_length() {
        for (base, expected) in [
            (b"".as_slice(), b"x".as_slice()),
            (b"|", b"|x"),
            (b"a", b"ax"),
            (b"||", b"||,x"),
            (b"a,", b"a,,x"),
        ] {
            let mut out = PropertyText::from_bytes(base);
            append_query_cx_extension(b"x", &mut out);
            assert_eq!(out.as_bytes(), expected);
        }
    }
    #[test]
    fn binary_input_remains_counted_and_unescaped() {
        let mut out = PropertyText::from_bytes(b"|\0");
        append_query_cx_extension(b"\xff\0,|", &mut out);
        assert_eq!(out.as_bytes(), b"|\0,\xff\0,|");
    }
    #[test]
    fn successive_appends_share_ordinary_owner_and_keep_literal_punctuation() {
        let mut a = PropertyText::from("|");
        let mut b = a.clone();
        for v in [b"c:1".as_slice(), b"", b",", b"t:2"] {
            append_query_cx_extension(v, &mut a);
            cosmolkit_smiles::append_cx_extension_source(v, &mut b);
        }
        assert_eq!(a, b);
        assert_eq!(a.as_bytes(), b"|c:1,,,t:2");
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod complete_cx_source_tests {
    use super::*;
    use cosmolkit_model::{
        AtomQueryPredicate, BondSpec, Conformer2D, PropertyValue, QueryAtomIdentity, QueryBond,
        QueryNode, SourceAtomValenceFacts, SubstanceGroupId, replace_query_substance_groups,
    };
    use cosmolkit_smiles::CxSmilesFields as F;
    use cosmolkit_types::ChiralTag;
    fn graph(n: usize, edges: &[(usize, usize)]) -> QueryGraph {
        QueryGraph::from_parts(
            (0..n)
                .map(|i| {
                    let mut a = QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    );
                    a.set_source_valence_facts(SourceAtomValenceFacts {
                        explicit_valence: 0,
                        implicit_valence: 0,
                    });
                    a
                })
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn output(q: &mut QueryGraph, flags: F) -> Result<PropertyText, SmartsWriteError> {
        let ao = (0..q.num_atoms()).map(AtomId::new).collect::<Vec<_>>();
        let bo = (0..q.num_bonds()).map(BondId::new).collect::<Vec<_>>();
        write_query_cx_extensions(q, &ao, &bo, flags)
    }
    fn cycle() -> QueryGraph {
        let mut q = graph(8, &(0..8).map(|i| (i, (i + 1) % 8)).collect::<Vec<_>>());
        q.bonds_mut()[0].bond_mut().set_order(BondOrder::Double);
        q.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::Any)
            .unwrap();
        q
    }
    fn atrop() -> QueryGraph {
        let mut q = graph(6, &[(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::AtropCw)
            .unwrap();
        q
    }
    #[test]
    fn flags_off_clear_actual_molecule_counter_without_acquiring_rings() {
        let mut q = graph(1, &[]);
        q.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(9))
            .unwrap();
        q.set_prop("keep", PropertyValue::Int(7)).unwrap();
        assert!(output(&mut q, F::NONE).unwrap().is_empty());
        assert!(q.prop("_cxsmilesOutputIndex").is_none());
        assert_eq!(q.prop("keep"), Some(&PropertyValue::Int(7)));
        assert!(!q.source_ring_info().initialized);
    }
    #[test]
    fn atom_preflight_error_precedes_cleanup_and_cache_acquisition() {
        let mut q = graph(1, &[]);
        q.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(9))
            .unwrap();
        assert!(matches!(
            write_query_cx_extensions(&mut q, &[AtomId::new(1)], &[], F::BOND_CFG),
            Err(SmartsWriteError::CxSourceAtomOutOfRange { .. })
        ));
        assert!(q.prop("_cxsmilesOutputIndex").is_some());
        assert!(!q.source_ring_info().initialized);
    }
    #[test]
    fn bond_cfg_initializes_actual_sssr_before_ring_extension() {
        let mut q = cycle();
        assert!(!q.source_ring_info().initialized);
        assert_eq!(
            output(&mut q, F::BOND_CFG).unwrap(),
            PropertyText::from("|ctu:0|")
        );
        assert!(q.source_ring_info().initialized);
        assert_eq!(q.source_ring_info().find_type, 1);
        assert_eq!(q.source_ring_info().bond_rings.len(), 1);
    }
    #[test]
    fn atrop_only_acquires_sssr_but_does_not_emit_ring_cis_trans() {
        let mut q = cycle();
        assert!(output(&mut q, F::BOND_ATROPISOMER).unwrap().is_empty());
        assert!(q.source_ring_info().initialized);
        assert_eq!(q.source_ring_info().find_type, 1);
    }
    #[test]
    fn initialized_empty_sssr_is_trusted_without_rebuilding_cycle() {
        let mut q = cycle();
        let r = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Sssr, 0, 0)
            .into_source_snapshot();
        q.replace_source_ring_info(r.clone());
        assert!(output(&mut q, F::BOND_CFG).unwrap().is_empty());
        assert_eq!(q.source_ring_info(), &r);
    }
    #[test]
    fn unknown_cache_transport_is_not_read_with_unrelated_flags() {
        let mut q = graph(1, &[]);
        let mut r = q.source_ring_info().clone();
        r.find_type = 77;
        q.replace_source_ring_info(r);
        assert!(output(&mut q, F::NONE).unwrap().is_empty());
        assert!(matches!(
            output(&mut q, F::BOND_CFG),
            Err(SmartsWriteError::CxRingInfo(
                cosmolkit_core::RingFindingError::SourceRingFindType { tag: 77 }
            ))
        ));
    }
    #[test]
    fn shared_atrop_map_reaches_enhanced_stereo_and_actual_bond_writes() {
        let mut q = atrop();
        cosmolkit_model::replace_query_stereo_groups(
            &mut q,
            vec![
                StereoGroup::new(StereoGroupKind::Absolute, vec![], vec![BondId::new(0)])
                    .expect("valid distinct stereo members"),
            ],
        )
        .unwrap();
        let atoms = q.atoms().to_vec();
        let text = output(&mut q, F::BOND_ATROPISOMER | F::ENHANCED_STEREO).unwrap();
        assert!(
            text.as_bytes()
                .windows(3)
                .any(|w| w == b"wU:" || w == b"wD:")
        );
        assert!(text.as_bytes().windows(2).any(|w| w == b"a:"));
        assert!(q.bonds().iter().any(|b| matches!(
            b.bond().direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash
        )));
        assert_eq!(q.atoms(), atoms);
    }
    #[test]
    fn existing_wedge_direction_write_without_new_map_entry_is_retained() {
        let mut q = atrop();
        q.bonds_mut()[1]
            .bond_mut()
            .set_direction(BondDirection::BeginWedge);
        output(&mut q, F::BOND_ATROPISOMER).unwrap();
        assert_eq!(q.bonds()[1].bond().direction(), BondDirection::BeginDash);
    }
    #[test]
    fn later_output_order_error_preserves_earlier_bond_and_cache_writes() {
        let mut q = atrop();
        q.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(7))
            .unwrap();
        let ao = (0..6).map(AtomId::new).collect::<Vec<_>>();
        assert!(matches!(
            write_query_cx_extensions(&mut q, &ao, &[BondId::new(99)], F::BOND_ATROPISOMER),
            Err(SmartsWriteError::CxSourceBondOutOfRange { .. })
        ));
        assert!(q.source_ring_info().initialized);
        assert!(q.bonds().iter().any(|b| matches!(
            b.bond().direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash
        )));
        assert_eq!(
            q.prop("_cxsmilesOutputIndex"),
            Some(&PropertyValue::UInt(7))
        );
    }
    #[test]
    fn cache_getter_error_retains_ring_acquisition_and_skips_final_cleanup() {
        let mut q = atrop();
        q.atoms_mut()[0].set_source_valence_facts(SourceAtomValenceFacts::UNINITIALIZED);
        q.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(7))
            .unwrap();
        assert!(matches!(
            output(&mut q, F::BOND_CFG),
            Err(SmartsWriteError::CxWedge(_))
        ));
        assert!(q.source_ring_info().initialized);
        assert!(q.prop("_cxsmilesOutputIndex").is_some());
    }
    #[test]
    fn bond_cfg_takes_precedence_over_atrop_only_and_uses_selected_source_coords() {
        let mut q = QueryGraph::from_parts(
            (0..4)
                .map(|i| {
                    let mut a = QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    );
                    if i == 1 {
                        a.set_chiral_tag(ChiralTag::TetrahedralCcw);
                    }
                    a
                })
                .collect(),
            [(0, 1), (1, 2), (1, 3)]
                .into_iter()
                .enumerate()
                .map(|(i, (a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![Conformer2D::new(
                19,
                vec![[1.0, 0.0], [0.0, 0.0], [-0.5, 0.866], [-0.5, -0.866]],
            )],
            vec![],
            vec![],
        )
        .unwrap();
        let mut only = q.clone();
        let flags = F::COORDS | F::BOND_ATROPISOMER;
        let only_text = output(&mut only, flags).unwrap();
        let text = output(&mut q, flags | F::BOND_CFG).unwrap();
        assert!(
            !only_text
                .as_bytes()
                .windows(2)
                .any(|w| w == b"wU" || w == b"wD")
        );
        assert!(text.as_bytes().windows(2).any(|w| w == b"wU" || w == b"wD"));
    }
    #[test]
    fn sgroup_counter_and_group_cleanup_are_installed_on_actual_query() {
        let mut q = graph(1, &[]);
        let mut g = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data);
        g = g.with_atoms(vec![AtomId::new(0)]);
        replace_query_substance_groups(&mut q, vec![g]).unwrap();
        q.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(9))
            .unwrap();
        assert!(
            output(&mut q, F::SGROUPS)
                .unwrap()
                .as_bytes()
                .starts_with(b"|SgD:0:")
        );
        assert!(q.prop("_cxsmilesOutputIndex").is_none());
        assert!(
            query_substance_groups(&q)[0]
                .props()
                .get(b"_cxsmilesOutputIndex".as_slice())
                .is_none()
        );
    }
    #[test]
    fn sgroup_hierarchy_error_retains_source_counter_and_group_prefix_writes() {
        let mut q = graph(1, &[]);
        let mut g = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data);
        g = g.with_atoms(vec![AtomId::new(0)]);
        g.set_prop("index", PropertyValue::String(PropertyText::from("bad")))
            .unwrap();
        replace_query_substance_groups(&mut q, vec![g]).unwrap();
        assert!(output(&mut q, F::SGROUPS).is_err());
        assert_eq!(
            q.prop("_cxsmilesOutputIndex"),
            Some(&PropertyValue::UInt(1))
        );
        assert_eq!(
            query_substance_groups(&q)[0]
                .props()
                .get(b"_cxsmilesOutputIndex".as_slice()),
            Some(&PropertyValue::UInt(0))
        );
    }
}
#[cfg(test)]
mod describe_query_source_tests {
    use super::*;
    #[test]
    fn disabled_logging_keeps_opaque_nonclone_payloads_borrowed_and_unmodified() {
        use std::{cell::Cell, rc::Rc};
        struct Opaque {
            drops: Rc<Cell<usize>>,
            value: u32,
        }
        impl Drop for Opaque {
            fn drop(&mut self) {
                self.drops.set(self.drops.get() + 1);
            }
        }
        let drops = Rc::new(Cell::new(0));
        let leaf = |value| {
            QueryNode::Predicate(Opaque {
                drops: drops.clone(),
                value,
            })
        };
        let nodes = [
            leaf(7),
            QueryNode::And(vec![leaf(1), QueryNode::Not(Box::new(leaf(2)))]),
            QueryNode::Or(vec![leaf(3), leaf(4)]),
            QueryNode::Xor(vec![leaf(5), leaf(6)]),
            QueryNode::Not(Box::new(QueryNode::Not(Box::new(leaf(8))))),
            QueryNode::And(vec![]),
        ];
        for node in &nodes {
            describe_query(node, "\0\t".to_owned());
        }
        assert_eq!(drops.get(), 0);
        match &nodes[0] {
            QueryNode::Predicate(p) => assert_eq!(p.value, 7),
            _ => panic!(),
        }
        drop(nodes);
        assert_eq!(drops.get(), 8);
    }
}

#[cfg(feature = "smiles-integration")]
fn mol_fragment_to_smarts_source(
    query: &QueryGraph,
    params: &SmartsWriteParams,
    atoms: &[AtomId],
    bonds: Option<&[BondId]>,
) -> Result<SmartsWriteResult, SmartsWriteError> {
    // RDKit❗❌: std::string MolFragmentToSmarts(const ROMol &mol,
    // RDKit❗❌:                                 const SmilesWriteParams &params,
    // RDKit❗❌:                                 const std::vector<int> &atomsToUse,
    // RDKit❗❌:                                 const std::vector<int> *bondsToUse) {
    // RDKit❗❌:   PRECONDITION(!atomsToUse.empty(), "no atoms provided");
    // RDKit❗❌:   PRECONDITION(!bondsToUse || !bondsToUse->empty(), "no bonds provided");
    // RDKit❗❌:
    // RDKit❗❌:   auto nAtoms = mol.getNumAtoms();
    // RDKit❗❌:   if (!nAtoms) {
    // RDKit❗❌:     return "";
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::unique_ptr<boost::dynamic_bitset<>> bondsInPlay(nullptr);
    // RDKit❗❌:   if (bondsToUse != nullptr) {
    // RDKit❗❌:     bondsInPlay.reset(new boost::dynamic_bitset<>(mol.getNumBonds(), 0));
    // RDKit❗❌:     for (auto bidx : *bondsToUse) {
    // RDKit❗❌:       bondsInPlay->set(bidx);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // Mark all atoms except the ones in atomIndices as already processed.
    // RDKit❗❌:   // white: unprocessed
    // RDKit❗❌:   // grey: partial
    // RDKit❗❌:   // black: complete
    // RDKit❗❌:   boost::dynamic_bitset<> atomsInPlay(nAtoms);
    // RDKit❗❌:   std::vector<AtomColors> colors(nAtoms, Canon::BLACK_NODE);
    // RDKit❗❌:   for (const auto &idx : atomsToUse) {
    // RDKit❗❌:     colors[idx] = Canon::WHITE_NODE;
    // RDKit❗❌:     atomsInPlay.set(idx);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   SmilesWriteParams ps(params);
    // RDKit❗❌:   ps.rootedAtAtom = -1;
    // RDKit❗❌:   return molToSmarts(mol, ps, std::move(colors), atomsInPlay,
    // RDKit❗❌:                      bondsInPlay.get());
    // RDKit❗❌: }
    // Preserve source precondition/empty-return/read order. Bit writes use
    // checked structural error translation, never sorting or filtering rows.
    if atoms.is_empty() {
        return Err(SmartsWriteError::EmptyAtomSelection);
    }
    if bonds.is_some_and(<[BondId]>::is_empty) {
        return Err(SmartsWriteError::EmptyBondSelection);
    }
    let n = query.num_atoms();
    if n == 0 {
        return Ok(SmartsWriteResult::default());
    }
    let bonds_in_play = if let Some(rows) = bonds {
        let mut mask = vec![false; query.num_bonds()];
        for row in rows {
            let slot = mask
                .get_mut(row.index())
                .ok_or(SmartsWriteError::FragmentBondOutOfRange { bond: row.index() })?;
            *slot = true;
        }
        Some(mask)
    } else {
        None
    };
    let mut atoms_in_play = vec![false; n];
    let mut colors = vec![cosmolkit_smiles::AtomColor::Black; n];
    for row in atoms {
        let color = colors
            .get_mut(row.index())
            .ok_or(SmartsWriteError::FragmentAtomOutOfRange { atom: row.index() })?;
        *color = cosmolkit_smiles::AtomColor::White;
        atoms_in_play[row.index()] = true;
    }
    let mut ps = *params;
    ps.rooted_at_atom = None;
    // None remains Native nullptr, not an eagerly synthesized induced-bond
    // mask. The sole complete molToSmarts owns all chemistry and source order.
    // Cost ❌: byte masks versus packed Native bitsets; no BTreeSet rebuild,
    // local ranking/duplicate traversal or cloned query in this wrapper.
    mol_to_smarts_source(query, &ps, colors, &atoms_in_play, bonds_in_play.as_deref())
}
#[cfg(all(test, feature = "smiles-integration"))]
mod mol_fragment_to_smarts_wrapper_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, PropertyValue};
    fn q(n: usize, edges: &[(usize, usize)]) -> QueryGraph {
        QueryGraph::from_parts(
            (0..n)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn run(
        q: &QueryGraph,
        a: &[usize],
        b: Option<&[usize]>,
    ) -> Result<SmartsWriteResult, SmartsWriteError> {
        let a: Vec<_> = a.iter().copied().map(AtomId::new).collect();
        let b = b.map(|b| b.iter().copied().map(BondId::new).collect::<Vec<_>>());
        mol_fragment_to_smarts_source(q, &Default::default(), &a, b.as_deref())
    }
    #[test]
    fn source_preconditions_precede_empty_graph_return() {
        let g = q(0, &[]);
        assert!(matches!(
            run(&g, &[], Some(&[])),
            Err(SmartsWriteError::EmptyAtomSelection)
        ));
        assert!(matches!(
            run(&g, &[99], Some(&[])),
            Err(SmartsWriteError::EmptyBondSelection)
        ));
    }
    #[test]
    fn empty_graph_returns_before_index_or_root_access() {
        let g = q(0, &[]);
        let o = mol_fragment_to_smarts_source(
            &g,
            &SmartsWriteParams {
                rooted_at_atom: Some(usize::MAX),
                ..Default::default()
            },
            &[AtomId::new(99)],
            Some(&[BondId::new(99)]),
        )
        .unwrap();
        assert!(o.smarts.is_empty());
        assert!(!o.source_orders_written);
    }
    #[test]
    fn bond_mask_access_precedes_atom_mask_access() {
        assert!(matches!(
            run(&q(1, &[]), &[99], Some(&[7])),
            Err(SmartsWriteError::FragmentBondOutOfRange { bond: 7 })
        ));
    }
    #[test]
    fn invalid_bond_rows_are_read_in_input_order_without_set_sort() {
        assert!(matches!(
            run(&q(1, &[]), &[0], Some(&[99, 7])),
            Err(SmartsWriteError::FragmentBondOutOfRange { bond: 99 })
        ));
    }
    #[test]
    fn invalid_atom_rows_are_read_in_input_order() {
        assert!(matches!(
            run(&q(1, &[]), &[99, 7], None),
            Err(SmartsWriteError::FragmentAtomOutOfRange { atom: 99 })
        ));
    }
    #[test]
    fn duplicate_reversed_atom_rows_become_source_bitmap_not_output_order() {
        let o = run(&q(3, &[]), &[2, 0, 2], None).unwrap();
        assert_eq!(o.smarts, PropertyText::from("[#6].[#6]"));
        assert_eq!(o.atom_ordering, [AtomId::new(0), AtomId::new(2)]);
        assert!(o.source_orders_written);
    }
    #[test]
    fn source_fragment_resets_out_of_range_root_before_canonical_call() {
        let g = q(2, &[(0, 1)]);
        let o = mol_fragment_to_smarts_source(
            &g,
            &SmartsWriteParams {
                rooted_at_atom: Some(usize::MAX),
                ..Default::default()
            },
            &[AtomId::new(1)],
            None,
        )
        .unwrap();
        assert_eq!(o.smarts, PropertyText::from("[#6]"));
        assert_eq!(o.atom_ordering, [AtomId::new(1)]);
    }
    #[test]
    fn absent_bond_pointer_and_explicit_mask_keep_distinct_source_inputs() {
        let g = q(3, &[(0, 1), (1, 2)]);
        let all = run(&g, &[0, 1, 2], None).unwrap();
        assert_eq!(all.smarts, PropertyText::from("[#6]-[#6]-[#6]"));
        assert_eq!(all.bond_ordering, [BondId::new(0), BondId::new(1)]);
        let cut = run(&g, &[0, 1, 2], Some(&[1])).unwrap();
        assert_eq!(cut.smarts, PropertyText::from("[#6].[#6]-[#6]"));
        assert_eq!(cut.bond_ordering, [BondId::new(1)]);
    }
    #[test]
    fn repeated_bond_rows_are_one_true_bit_without_duplicate_emission() {
        let o = run(&q(2, &[(0, 1)]), &[0, 1], Some(&[0, 0])).unwrap();
        assert_eq!(o.smarts, PropertyText::from("[#6]-[#6]"));
        assert_eq!(o.bond_ordering, [BondId::new(0)]);
    }
    #[test]
    fn unselected_query_bond_is_not_emitted_or_evaluated() {
        let mut g = q(2, &[(0, 1)]);
        g.bonds_mut()[0].set_predicate(QueryNode::predicate(BondQueryPredicate::HasStereo));
        let o = run(&g, &[0], None).unwrap();
        assert_eq!(o.smarts, PropertyText::from("[#6]"));
        assert!(o.bond_ordering.is_empty());
    }
    #[test]
    fn ring_fragment_reuses_canonical_ring_stack_and_preserves_input_cache() {
        let g = q(3, &[(0, 1), (1, 2), (2, 0)]);
        let before = g.clone();
        let o = run(&g, &[2, 1, 0], None).unwrap();
        assert_eq!(o.smarts, PropertyText::from("[#6]1-[#6]-[#6]-1"));
        assert_eq!(o.atom_ordering, [0, 1, 2].map(AtomId::new));
        assert_eq!(g, before);
        assert!(!g.source_ring_info().initialized);
    }
    #[test]
    fn detached_public_adapter_reuses_complete_source_and_keeps_input_values() {
        let mut g = q(2, &[(0, 1)]);
        g.set_prop("_Name", PropertyValue::String("kept".into()))
            .unwrap();
        g.add_conformer_3d(cosmolkit_model::Conformer3D::new(
            7,
            vec![[1., 2., 3.], [4., 5., 6.]],
            false,
        ))
        .unwrap();
        let before = g.clone();
        assert_eq!(
            query_graph_fragment_to_smarts(&g, &Default::default(), &[AtomId::new(1)], None)
                .unwrap(),
            PropertyText::from("[#6]")
        );
        assert_eq!(g, before);
    }
}

#[cfg(feature = "smiles-integration")]
fn mol_to_cx_smarts_source(
    query: &mut QueryGraph,
    params: &SmartsWriteParams,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string MolToCXSmarts(const ROMol &mol, const SmilesWriteParams &params) {
    // RDKit❗❌:   SmilesWriteParams ps(params);
    // RDKit❗❌:   ps.includeDativeBonds = false;
    // RDKit❗❌:   auto res = MolToSmarts(mol, ps);
    // RDKit❗❌:   if (!res.empty()) {
    // RDKit❗❌:     auto cxext = SmilesWrite::getCXExtensions(mol);
    // RDKit❗❌:     if (!cxext.empty()) {
    // RDKit❗❌:       res += " " + cxext;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    let mut ps = *params;
    ps.include_dative_bonds = false;
    let mut result = mol_to_smarts_wrapper_source(query, &ps)?;
    if let Some(properties) = result.source_properties.take() {
        query.replace_source_molecule_properties(&properties);
    }
    if !result.smarts.is_empty() {
        // Use the SAME original source graph after actual order-property
        // recording. Never substitute molToSmarts' private working ring cache.
        let cxext = write_query_cx_extensions(
            query,
            &result.atom_ordering,
            &result.bond_ordering,
            cosmolkit_smiles::CxSmilesFields::ALL,
        )?;
        if !cxext.is_empty() {
            result.smarts.push_byte(b' ');
            result.smarts.extend_bytes(cxext.as_bytes());
        }
    }
    // One source option copy, one canonical SMARTS call, conditional CX call
    // and counted-byte append. No Graph clone in this mutable source body;
    // inherited detached property/cache/writer costs remain known❌.
    Ok(result.smarts)
}
#[cfg(all(test, feature = "smiles-integration"))]
mod mol_to_cx_smarts_wrapper_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, MoleculePropertyError, PropertyValue};
    fn q(n: usize) -> QueryGraph {
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
    fn dative() -> QueryGraph {
        QueryGraph::from_parts(
            q(2).atoms().to_vec(),
            vec![QueryBond::new(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Dative),
            )],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn empty_smarts_skips_cx_setters_cache_and_cleanup() {
        let mut g = q(0);
        g.set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        g.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(7))
            .unwrap();
        let before = g.clone();
        assert!(
            mol_to_cx_smarts_source(
                &mut g,
                &SmartsWriteParams {
                    rooted_at_atom: Some(99),
                    ..Default::default()
                }
            )
            .unwrap()
            .is_empty()
        );
        assert_eq!(g, before);
    }
    #[test]
    fn whole_source_forces_dative_false_but_retains_coordinate_bond_extension() {
        let mut g = dative();
        let p = SmartsWriteParams::default();
        assert!(p.include_dative_bonds);
        assert_eq!(
            mol_to_cx_smarts_source(&mut g, &p).unwrap(),
            PropertyText::from("[#6]-[#6] |C:0.0|")
        );
        assert!(p.include_dative_bonds);
        assert!(g.source_ring_info().initialized);
    }
    #[test]
    fn root_order_is_shared_with_cx_coordinate_bond_writer() {
        let mut g = dative();
        let p = SmartsWriteParams {
            rooted_at_atom: Some(1),
            ..Default::default()
        };
        assert_eq!(
            mol_to_cx_smarts_source(&mut g, &p).unwrap(),
            PropertyText::from("[#6]-[#6] |C:1.0|")
        );
    }
    #[test]
    fn order_setter_failure_occurs_before_any_cx_ring_or_cleanup_effect() {
        let mut g = q(1);
        g.set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        g.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(7))
            .unwrap();
        // Canon.cpp first perceives stereo if this property is absent. Its
        // computed setter would fail earlier, before output-order setters.
        // Native tests presence, so false still bypasses that earlier stage.
        g.set_prop("_StereochemDone", PropertyValue::Bool(false))
            .unwrap();
        let before = g.clone();
        assert!(matches!(
            mol_to_cx_smarts_source(&mut g, &Default::default()),
            Err(SmartsWriteError::MoleculePropertyWrite(
                MoleculePropertyError::ComputedListKind(_)
            ))
        ));
        assert_eq!(g, before);
    }
    #[test]
    fn computed_names_keep_existing_duplicates_and_native_order() {
        let mut g = q(1);
        g.set_prop(
            "__computedProps",
            PropertyValue::StringVector(vec![
                "kept".into(),
                "_smilesAtomOutputOrder".into(),
                "_smilesAtomOutputOrder".into(),
            ]),
        )
        .unwrap();
        assert_eq!(
            mol_to_cx_smarts_source(&mut g, &Default::default()).unwrap(),
            PropertyText::from("[#6]")
        );
        assert_eq!(
            g.computed_prop_names().unwrap().unwrap(),
            [
                "kept",
                "_smilesAtomOutputOrder",
                "_smilesAtomOutputOrder",
                "_smilesBondOutputOrder"
            ]
            .map(PropertyText::from)
        );
        assert!(g.source_ring_info().initialized);
    }
    #[test]
    fn source_error_after_order_recording_preserves_order_names_and_ring_prefix() {
        let mut g = q(4);
        g.set_prop("_molLinkNodes", "1 3 2 2 0 2 3").unwrap();
        g.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(7))
            .unwrap();
        assert!(mol_to_cx_smarts_source(&mut g, &Default::default()).is_err());
        assert_eq!(
            g.computed_prop_names().unwrap().unwrap(),
            ["_smilesAtomOutputOrder", "_smilesBondOutputOrder"].map(PropertyText::from)
        );
        assert!(g.source_ring_info().initialized);
        assert_eq!(
            g.prop("_cxsmilesOutputIndex"),
            Some(&PropertyValue::UInt(7))
        );
    }
    #[test]
    fn immutable_public_value_adapter_keeps_original_cx_properties_and_cache() {
        let mut g = q(1);
        g.atom_mut(0)
            .unwrap()
            .set_prop("atomLabel", "label")
            .unwrap();
        g.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(7))
            .unwrap();
        let before = g.clone();
        assert_eq!(
            query_graph_to_cx_smarts(&g, &Default::default()).unwrap(),
            PropertyText::from("[#6] |$label$|")
        );
        assert_eq!(g, before);
        assert!(!g.source_ring_info().initialized);
    }
    #[test]
    fn rooted_index_error_precedes_source_computed_property_read() {
        let mut g = q(1);
        g.set_prop("__computedProps", PropertyValue::UInt(4))
            .unwrap();
        assert!(matches!(
            mol_to_cx_smarts_source(
                &mut g,
                &SmartsWriteParams {
                    rooted_at_atom: Some(99),
                    ..Default::default()
                }
            ),
            Err(SmartsWriteError::RootedAtomOutOfRange { atom: 99 })
        ));
        assert!(!g.source_ring_info().initialized);
    }
}

#[cfg(all(test, feature = "smiles-integration"))]
mod cx_smarts_earlier_stereo_computed_failure_source_tests {
    use super::*;
    #[test]
    fn missing_stereo_done_fails_its_computed_setter_before_output_order_recording() {
        let mut g = QueryGraph::from_parts(
            vec![QueryAtom::new(
                AtomId::new(0),
                cosmolkit_model::AtomSpec::new(Element::C),
            )],
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        g.set_prop(
            "__computedProps",
            cosmolkit_model::PropertyValue::Bool(false),
        )
        .unwrap();
        let before = g.clone();
        assert!(matches!(
            mol_to_cx_smarts_source(&mut g, &Default::default()),
            Err(SmartsWriteError::CanonicalTraversal(
                cosmolkit_smiles::SmilesParseError::MoleculeProperty(
                    cosmolkit_model::MoleculePropertyError::ComputedListKind(_)
                )
            ))
        ));
        assert_eq!(g, before);
        assert!(!g.source_ring_info().initialized);
    }
}

#[cfg(feature = "smiles-integration")]
fn mol_fragment_to_cx_smarts_source(
    query: &mut QueryGraph,
    params: &SmartsWriteParams,
    atoms: &[AtomId],
    bonds: Option<&[BondId]>,
) -> Result<PropertyText, SmartsWriteError> {
    // RDKit❗❌: std::string MolFragmentToCXSmarts(const ROMol &mol,
    // RDKit❗❌:                                   const SmilesWriteParams &params,
    // RDKit❗❌:                                   const std::vector<int> &atomsToUse,
    // RDKit❗❌:                                   const std::vector<int> *bondsToUse) {
    // RDKit❗❌:   auto res = MolFragmentToSmarts(mol, params, atomsToUse, bondsToUse);
    // RDKit❗❌:   if (!res.empty()) {
    // RDKit❗❌:     auto cxext = SmilesWrite::getCXExtensions(mol);
    // RDKit❗❌:     if (!cxext.empty()) {
    // RDKit❗❌:       res += " " + cxext;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // Fragment wrapper preserves the caller's dative option. The canonical
    // FragmentToSmarts owner alone clears root and creates source bitmasks.
    let mut result = mol_fragment_to_smarts_source(query, params, atoms, bonds)?;
    if let Some(properties) = result.source_properties.take() {
        query.replace_source_molecule_properties(&properties);
    }
    if !result.smarts.is_empty() {
        let cxext = write_query_cx_extensions(
            query,
            &result.atom_ordering,
            &result.bond_ordering,
            cosmolkit_smiles::CxSmilesFields::ALL,
        )?;
        if !cxext.is_empty() {
            result.smarts.push_byte(b' ');
            result.smarts.extend_bytes(cxext.as_bytes());
        }
    }
    // Same source graph and actual emitted orders; excluded rows retain the
    // source writer's reverse-order zero initialization. No inverse guesses,
    // subset reconstruction or working-ring-cache substitution.
    // O(1) wrapper plus linear append; inherited detached property/cache and
    // byte-mask costs remain known❌. No graph clone in this source kernel.
    Ok(result.smarts)
}
#[cfg(all(test, feature = "smiles-integration"))]
mod mol_fragment_to_cx_smarts_wrapper_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, PropertyValue};
    fn q(n: usize) -> QueryGraph {
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
    #[test]
    fn empty_selection_error_precedes_all_cx_effects() {
        let mut g = q(1);
        g.set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        let before = g.clone();
        assert!(matches!(
            mol_fragment_to_cx_smarts_source(&mut g, &Default::default(), &[], None),
            Err(SmartsWriteError::EmptyAtomSelection)
        ));
        assert_eq!(g, before);
    }
    #[test]
    fn empty_graph_nonempty_selection_returns_before_cx_or_row_checks() {
        let mut g = q(0);
        g.set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        let before = g.clone();
        assert!(
            mol_fragment_to_cx_smarts_source(
                &mut g,
                &Default::default(),
                &[AtomId::new(99)],
                Some(&[BondId::new(99)])
            )
            .unwrap()
            .is_empty()
        );
        assert_eq!(g, before);
    }
    #[test]
    fn fragment_preserves_both_explicit_dative_options() {
        for include in [false, true] {
            let mut g = QueryGraph::from_parts(
                q(2).atoms().to_vec(),
                vec![QueryBond::new(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Dative),
                )],
                [],
                vec![],
                vec![],
                vec![],
            )
            .unwrap();
            let p = SmartsWriteParams {
                include_dative_bonds: include,
                ..Default::default()
            };
            let out = mol_fragment_to_cx_smarts_source(
                &mut g,
                &p,
                &[AtomId::new(0), AtomId::new(1)],
                None,
            )
            .unwrap();
            assert_eq!(
                out,
                PropertyText::from(if include {
                    "[#6]->[#6] |C:0.0|"
                } else {
                    "[#6]-[#6] |C:0.0|"
                })
            );
        }
    }
    #[test]
    fn fragment_ignores_root_and_cx_uses_actual_selected_atom_order() {
        let mut g = q(2);
        g.atom_mut(0)
            .unwrap()
            .set_prop("atomLabel", "first")
            .unwrap();
        g.atom_mut(1)
            .unwrap()
            .set_prop("atomLabel", "second")
            .unwrap();
        let p = SmartsWriteParams {
            rooted_at_atom: Some(usize::MAX),
            ..Default::default()
        };
        assert_eq!(
            mol_fragment_to_cx_smarts_source(&mut g, &p, &[AtomId::new(1)], None).unwrap(),
            PropertyText::from("[#6] |$second$|")
        );
        assert!(g.source_ring_info().initialized);
    }
    #[test]
    fn empty_extension_adds_no_space_but_retains_recorded_names_and_cache() {
        let mut g = q(1);
        assert_eq!(
            mol_fragment_to_cx_smarts_source(&mut g, &Default::default(), &[AtomId::new(0)], None)
                .unwrap(),
            PropertyText::from("[#6]")
        );
        assert!(g.source_ring_info().initialized);
        assert_eq!(
            g.computed_prop_names().unwrap().unwrap(),
            ["_smilesAtomOutputOrder", "_smilesBondOutputOrder"].map(PropertyText::from)
        );
    }
    #[test]
    fn later_cx_failure_retains_actual_order_metadata_and_ring_prefix() {
        let mut g = q(4);
        g.set_prop("_molLinkNodes", "1 3 2 2 0 2 3").unwrap();
        g.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(7))
            .unwrap();
        assert!(
            mol_fragment_to_cx_smarts_source(&mut g, &Default::default(), &[AtomId::new(1)], None)
                .is_err()
        );
        assert!(g.source_ring_info().initialized);
        assert_eq!(
            g.computed_prop_names().unwrap().unwrap(),
            ["_smilesAtomOutputOrder", "_smilesBondOutputOrder"].map(PropertyText::from)
        );
        assert_eq!(
            g.prop("_cxsmilesOutputIndex"),
            Some(&PropertyValue::UInt(7))
        );
    }
    #[test]
    fn fragment_coordinates_use_source_front_and_selected_coordinate_row() {
        let mut g = q(2);
        g.add_conformer_3d(cosmolkit_model::Conformer3D::new(
            9,
            vec![[1., 2., 3.], [4., 5., 6.]],
            true,
        ))
        .unwrap();
        g.add_conformer_3d(cosmolkit_model::Conformer3D::new(
            1,
            vec![[7., 8., 9.], [10., 11., 12.]],
            true,
        ))
        .unwrap();
        assert_eq!(
            mol_fragment_to_cx_smarts_source(&mut g, &Default::default(), &[AtomId::new(1)], None)
                .unwrap(),
            PropertyText::from("[#6] |(4,5,6)|")
        );
        assert_eq!(g.conformers_3d()[0].id(), 9);
        assert_eq!(g.conformers_3d().len(), 2);
    }
    #[test]
    fn public_immutable_adapter_keeps_graph_properties_and_coordinate_state() {
        let mut g = q(2);
        g.atom_mut(1)
            .unwrap()
            .set_prop("atomLabel", "selected")
            .unwrap();
        g.set_prop("_cxsmilesOutputIndex", PropertyValue::UInt(7))
            .unwrap();
        let before = g.clone();
        assert_eq!(
            query_graph_fragment_to_cx_smarts(&g, &Default::default(), &[AtomId::new(1)], None)
                .unwrap(),
            PropertyText::from("[#6] |$selected$|")
        );
        assert_eq!(g, before);
        assert!(!g.source_ring_info().initialized);
    }
}
