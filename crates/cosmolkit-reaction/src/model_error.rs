use cosmolkit_model::QueryGraphError;

/// Source template role, retaining source vector order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionRole {
    Reactant,
    Product,
    Agent,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ReactionModelError {
    #[error("invalid {role:?} template {template}: {source}")]
    Template {
        role: ReactionRole,
        template: usize,
        #[source]
        source: QueryGraphError,
    },
    #[error("{role:?} template index {index} is outside {count} templates")]
    TemplateIndex {
        role: ReactionRole,
        index: usize,
        count: usize,
    },
}
