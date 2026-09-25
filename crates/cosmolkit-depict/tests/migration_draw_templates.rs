// Exercise the private DRAW-templates owner without exposing a mutation API.
// The source module is compiled as this test target; no second production
// implementation or public template registry is introduced.
#[path = "../src/templates.rs"]
mod templates;
