// Exercise the private DRAW-templates owner without exposing a mutation API.
// The source module is compiled as this test target; no second production
// implementation or public template registry is introduced.
// Template-match regressions also exercise the existing private fragment owner.
#[path = "../src/embedded_frag.rs"]
mod embedded_frag;
#[path = "../src/geometry.rs"]
mod geometry;
#[path = "../src/nontetrahedral.rs"]
mod nontetrahedral;
#[path = "../src/templates.rs"]
mod templates;

use cosmolkit_depict::{Compute2DCoordinatesParams, compute_2d_coordinates};
