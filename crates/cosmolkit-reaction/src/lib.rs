//! Detached reaction values and source-backed execution.
//!
//! The facade retains the sole live Molecule and runtime commit authority.

mod local_params;
mod management;
mod model_error;
mod product;
mod reaction;

pub use local_params::ReactionTemplateRemovalParams;
pub use management::ReactionTemplateRemoval;
#[doc(hidden)]
pub use management::{without_agents, without_unmapped_products, without_unmapped_reactants};
pub use model_error::{ReactionModelError, ReactionRole};
#[doc(hidden)]
pub use product::{ReactionInput, ReactionProduct, ReactionRowOrigin};
pub use reaction::Reaction;

mod parse;
mod parse_error;
mod template_stereo;
pub use parse::{ReactionParseParams, parse_smirks, parse_smirks_with_params};
pub use parse_error::ReactionParseError;

mod validation;
mod validation_error;
pub use validation::{
    ReactionValidationIssue, ReactionValidationIssueKind, ReactionValidationParams,
    ReactionValidationReport, ReactionValidationSeverity,
};
#[doc(hidden)]
pub use validation::{initialize_reaction, validate_reaction};
pub use validation_error::{ReactionInitializationError, ReactionValidationError};

mod matching;
mod run_error;
pub use run_error::ReactionRunError;

mod materialize;
mod product_error;
#[doc(hidden)]
pub use product_error::ReactionProductError;

mod coordinate_selection;
pub use coordinate_selection::ReactionCoordinateSelection;
mod product_stereo;

mod run_params;
pub use run_params::{ReactionRunParams, ReactionSingleRunParams};
mod runner;
#[doc(hidden)]
pub use runner::{run_reactant, run_reactants};

mod apply_params;
pub use apply_params::ReactionApplyParams;
mod apply_error;
pub use apply_error::ReactionApplyError;
mod apply;
#[doc(hidden)]
pub use apply::{ReactionApplyChanges, apply_reaction};

mod write_params;
pub use write_params::ReactionWriteParams;
mod write_error;
pub use write_error::ReactionWriteError;
mod write;
