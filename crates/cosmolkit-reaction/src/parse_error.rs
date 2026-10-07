use crate::ReactionRole;
use cosmolkit_model::{AtomId, QueryGraphError};

#[derive(Debug, thiserror::Error)]
pub enum ReactionParseError {
    #[error("parser cleanup at {role:?} template {template}: {source}")]
    ParserCleanup {
        role: ReactionRole,
        template: usize,
        #[source]
        source: cosmolkit_search::SmartsParseError,
    },
    #[error("reaction template property read failed: {source}")]
    TemplateProperty {
        #[source]
        source: crate::ReactionValidationError,
    },
    #[error("a reaction requires at least two > characters (found {count})")]
    Separators { count: usize },
    #[error("multi-step reactions not supported (found {count} separators)")]
    MultiStep { count: usize },
    #[error("reaction component byte slice {start}..{end} is invalid")]
    ComponentBounds { start: usize, end: usize },
    #[error("reaction atom or bond offset exceeds source unsigned range")]
    OffsetOverflow,
    #[error("problems constructing {role:?} template {template} from SMARTS {text:?}: {source}")]
    Smarts {
        role: ReactionRole,
        template: usize,
        text: String,
        #[source]
        source: cosmolkit_search::SmartsParseError,
    },
    #[error("problems constructing {role:?} template {template} from SMILES {text:?}: {source}")]
    Smiles {
        role: ReactionRole,
        template: usize,
        text: String,
        #[source]
        source: cosmolkit_smiles::SmilesParseError,
    },
    #[error("invalid {role:?} template {template}: {source}")]
    Model {
        role: ReactionRole,
        template: usize,
        #[source]
        source: QueryGraphError,
    },
    #[error("failed fragmenting agent template: {source}")]
    AgentFragments {
        #[source]
        source: cosmolkit_search::SmartsParseError,
    },
    #[error(
        "CX parsing for {role:?} template {template} at atom {start_atom}, bond {start_bond}: {source}"
    )]
    CxParse {
        role: ReactionRole,
        template: usize,
        start_atom: usize,
        start_bond: usize,
        #[source]
        source: cosmolkit_cx::CxParseError,
    },
    #[error(
        "CX lowering for {role:?} template {template} at atom {start_atom}, bond {start_bond}: {source}"
    )]
    CxLowering {
        role: ReactionRole,
        template: usize,
        start_atom: usize,
        start_bond: usize,
        #[source]
        source: cosmolkit_search::CxQueryLoweringError,
    },
    #[error("stereo order at product {product}, atom {atom}: {source}")]
    StereoOrder {
        product: usize,
        atom: AtomId,
        #[source]
        source: cosmolkit_core::StereoOrderError,
    },
}
