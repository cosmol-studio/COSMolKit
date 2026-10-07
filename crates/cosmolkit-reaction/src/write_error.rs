use crate::ReactionRole;
#[derive(Debug, thiserror::Error)]
pub enum ReactionWriteError {
    #[error("{role:?} template {template} SMARTS write: {source}")]
    Template {
        role: ReactionRole,
        template: usize,
        #[source]
        source: cosmolkit_search::SmartsWriteError,
    },
    #[error("{role:?} template {template} connected components: {source}")]
    Connectivity {
        role: ReactionRole,
        template: usize,
        #[source]
        source: cosmolkit_core::PathError,
    },
    #[error("received {actual} coordinate selectors for {expected} reaction templates")]
    CoordinateSelectionArity { expected: usize, actual: usize },
}
