//! Public SDF record value. Parsing and chemistry finalization belong to
//! `cosmolkit-io` and its detached algorithm dependencies.

use std::fmt;

use crate::{
    CoordinateDimension, Molecule, MoleculeProperties, OperationError, QueryGraph, SubstanceGroup,
};

pub use cosmolkit_io::SdfCoordinateMode;

/// Public policy for reading one SDF record into its final Rust graph value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SdfReadParams {
    pub sanitize: bool,
    pub remove_hydrogens: bool,
    pub strict_parsing: bool,
    pub expand_attachment_points: bool,
    pub process_property_lists: bool,
    pub coordinate_mode: SdfCoordinateMode,
}

impl Default for SdfReadParams {
    fn default() -> Self {
        Self {
            sanitize: true,
            remove_hydrogens: true,
            strict_parsing: true,
            expand_attachment_points: false,
            process_property_lists: true,
            coordinate_mode: SdfCoordinateMode::Preserve,
        }
    }
}

/// The one existing chemistry value carried by a finalized SDF record.
#[derive(Clone, Debug)]
pub enum SdfGraph {
    /// A runtime-validated concrete molecule.
    Molecule(Molecule),
    /// The canonical detached query graph, without lossy concrete lowering.
    Query(QueryGraph),
}

impl SdfGraph {
    const fn kind(&self) -> &'static str {
        match self {
            Self::Molecule(_) => "molecule",
            Self::Query(_) => "query_graph",
        }
    }
}

/// Failure of the public SDF record and molecule-reading boundary.
#[derive(Debug)]
pub enum SdfError {
    /// Detached record framing, parsing, or property-list processing failed.
    Read(cosmolkit_io::SdfReadError),
    /// Source-ordered chemistry finalization failed.
    Post(cosmolkit_io::MolPostError),
    /// Validated live molecule construction failed.
    Construction(OperationError),
    /// A concrete-only reader encountered a finalized query record.
    QueryRecord,
    /// An accessor requested the other tagged graph payload.
    WrongGraphKind {
        expected: &'static str,
        actual: &'static str,
    },
}

impl fmt::Display for SdfError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Read(error) => write!(formatter, "SDF reading failed: {error}"),
            Self::Post(error) => write!(formatter, "SDF finalization failed: {error}"),
            Self::Construction(error) => write!(formatter, "SDF construction failed: {error}"),
            Self::QueryRecord => formatter
                .write_str("query-bearing SDF record cannot be represented as a concrete molecule"),
            Self::WrongGraphKind { expected, actual } => {
                write!(formatter, "expected SDF {expected} payload, found {actual}")
            }
        }
    }
}

impl std::error::Error for SdfError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Read(error) => Some(error),
            Self::Post(error) => Some(error),
            Self::Construction(error) => Some(error),
            Self::QueryRecord | Self::WrongGraphKind { .. } => None,
        }
    }
}

impl From<cosmolkit_io::SdfReadError> for SdfError {
    fn from(error: cosmolkit_io::SdfReadError) -> Self {
        match error {
            cosmolkit_io::SdfReadError::QueryRecord => Self::QueryRecord,
            other => Self::Read(other),
        }
    }
}

impl From<cosmolkit_io::MolPostError> for SdfError {
    fn from(error: cosmolkit_io::MolPostError) -> Self {
        Self::Post(error)
    }
}

impl From<OperationError> for SdfError {
    fn from(error: OperationError) -> Self {
        Self::Construction(error)
    }
}

/// One finalized SDF record with its ordered raw fields and typed state.
///
/// The record does not define a third molecule model. Its graph is exactly
/// one `Molecule` or `QueryGraph` value; accessors borrow without mutation.
#[derive(Clone, Debug)]
pub struct SdfRecord {
    graph: SdfGraph,
    data_fields: Vec<(String, String)>,
    properties: MoleculeProperties,
    substance_groups: Vec<SubstanceGroup>,
    source_coordinate_dim: Option<CoordinateDimension>,
}

impl SdfRecord {
    /// Read and finalize the first SDF record using the default source policy.
    pub fn from_sdf(text: &str) -> Result<Self, SdfError> {
        Self::from_sdf_with_params(text, &SdfReadParams::default())
    }

    /// Read and finalize the first SDF record without lowering query graphs.
    pub fn from_sdf_with_params(text: &str, params: &SdfReadParams) -> Result<Self, SdfError> {
        // One-record framing, parsing and property lists stay in the IO owner;
        // finalization uses that owner's parser-retained chirality bit.
        let parsed = cosmolkit_io::read_sdf_graph_record_detached_with_params(
            text,
            cosmolkit_io::SdfDataReadParams {
                strict_parsing: params.strict_parsing,
                process_property_lists: params.process_property_lists,
                coordinate_mode: params.coordinate_mode,
            },
        )?;
        let finalized = parsed.finish_mol_post(cosmolkit_io::MolPostParams {
            sanitize: params.sanitize,
            remove_hs: params.remove_hydrogens,
            expand_attachment_points: params.expand_attachment_points,
        })?;
        let data_fields = finalized.data_fields;
        let mol_block = finalized.mol_block;
        match mol_block {
            cosmolkit_io::MolBlockRecord::Concrete {
                topology,
                coordinates,
                properties,
            } => {
                let substance_groups = topology.substance_groups.clone();
                let source_coordinate_dim = coordinates.source_coordinate_dim;
                let graph = SdfGraph::Molecule(Molecule::from_validated_parts(
                    topology,
                    coordinates,
                    properties.clone(),
                )?);
                Ok(Self::from_finalized_graph(
                    graph,
                    data_fields,
                    properties,
                    substance_groups,
                    source_coordinate_dim,
                ))
            }
            cosmolkit_io::MolBlockRecord::Query(record) => Ok(Self::from_finalized_graph(
                SdfGraph::Query(record.query),
                data_fields,
                record.properties,
                record.substance_groups,
                record.source_coordinate_dim,
            )),
        }
    }

    /// Receives a graph only after detached finalization and, for a concrete
    /// payload, private runtime construction validation have succeeded.
    pub(crate) fn from_finalized_graph(
        graph: SdfGraph,
        data_fields: Vec<(String, String)>,
        properties: MoleculeProperties,
        substance_groups: Vec<SubstanceGroup>,
        source_coordinate_dim: Option<CoordinateDimension>,
    ) -> Self {
        Self {
            graph,
            data_fields,
            properties,
            substance_groups,
            source_coordinate_dim,
        }
    }

    #[must_use]
    pub const fn graph(&self) -> &SdfGraph {
        &self.graph
    }

    pub fn molecule(&self) -> Result<&Molecule, SdfError> {
        match &self.graph {
            SdfGraph::Molecule(molecule) => Ok(molecule),
            other => Err(SdfError::WrongGraphKind {
                expected: "molecule",
                actual: other.kind(),
            }),
        }
    }

    pub fn query_graph(&self) -> Result<&QueryGraph, SdfError> {
        match &self.graph {
            SdfGraph::Query(query) => Ok(query),
            other => Err(SdfError::WrongGraphKind {
                expected: "query_graph",
                actual: other.kind(),
            }),
        }
    }

    #[must_use]
    pub fn data_fields(&self) -> &[(String, String)] {
        &self.data_fields
    }

    #[must_use]
    pub const fn properties(&self) -> &MoleculeProperties {
        &self.properties
    }

    #[must_use]
    pub fn substance_groups(&self) -> &[SubstanceGroup] {
        &self.substance_groups
    }

    #[must_use]
    pub const fn source_coordinate_dim(&self) -> Option<CoordinateDimension> {
        self.source_coordinate_dim
    }
}

impl Molecule {
    /// Read the first SDF record, rejecting a finalized query-bearing graph.
    pub fn from_sdf(text: &str) -> Result<Self, SdfError> {
        Self::from_sdf_with_params(text, &SdfReadParams::default())
    }

    /// Read one SDF record into the validated concrete live-molecule boundary.
    pub fn from_sdf_with_params(text: &str, params: &SdfReadParams) -> Result<Self, SdfError> {
        match SdfRecord::from_sdf_with_params(text, params)?.graph {
            SdfGraph::Molecule(molecule) => Ok(molecule),
            SdfGraph::Query(_) => Err(SdfError::QueryRecord),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use crate::{CoordinateBlock, SubstanceGroupId, SubstanceGroupKind, TopologyBlock};

    use super::*;

    #[test]
    fn sdf_record_result_boundary_concrete_preserves_ordered_typed_state() {
        let fields = vec![
            ("ID".into(), "first".into()),
            ("ID".into(), "second".into()),
        ];
        let properties = MoleculeProperties::default()
            .with_name("record")
            .with_sdf_data_field("ID", "first")
            .with_sdf_data_field("ID", "second");
        let molecule = Molecule::from_validated_parts(
            TopologyBlock::default(),
            CoordinateBlock::default(),
            properties.clone(),
        )
        .expect("private validated construction");
        let group = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_label("typed");
        let record = SdfRecord::from_finalized_graph(
            SdfGraph::Molecule(molecule),
            fields.clone(),
            properties.clone(),
            vec![group.clone()],
            Some(CoordinateDimension::TwoD),
        );

        assert!(matches!(record.graph(), SdfGraph::Molecule(_)));
        assert!(record.molecule().is_ok());
        assert!(matches!(
            record.query_graph(),
            Err(SdfError::WrongGraphKind {
                expected: "query_graph",
                actual: "molecule"
            })
        ));
        assert_eq!(record.data_fields(), fields);
        assert_eq!(record.properties(), &properties);
        assert_eq!(record.substance_groups(), &[group]);
        assert_eq!(
            record.source_coordinate_dim(),
            Some(CoordinateDimension::TwoD)
        );
        // Every accessor is an immutable borrow; inspecting a clone cannot
        // mutate the original record's graph or ordered metadata.
        let clone = record.clone();
        assert_eq!(clone.data_fields(), record.data_fields());
        assert_eq!(clone.properties(), record.properties());
        assert!(record.molecule().is_ok());
    }

    #[test]
    fn sdf_record_result_boundary_query_preserves_graph_and_wrong_kind_error() {
        let query = QueryGraph::from_parts(vec![], vec![], BTreeMap::new(), vec![], vec![], vec![])
            .expect("empty detached query graph is valid");
        let record = SdfRecord::from_finalized_graph(
            SdfGraph::Query(query.clone()),
            vec![("Q".into(), "one".into())],
            MoleculeProperties::default(),
            vec![],
            Some(CoordinateDimension::ThreeD),
        );

        assert!(matches!(record.graph(), SdfGraph::Query(_)));
        assert_eq!(record.query_graph().expect("query payload"), &query);
        assert!(matches!(
            record.molecule(),
            Err(SdfError::WrongGraphKind {
                expected: "molecule",
                actual: "query_graph"
            })
        ));
        assert_eq!(record.data_fields(), &[("Q".into(), "one".into())]);
        assert_eq!(
            record.source_coordinate_dim(),
            Some(CoordinateDimension::ThreeD)
        );
    }
}
