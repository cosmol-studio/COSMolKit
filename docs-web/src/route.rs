use dioxus::prelude::*;

use crate::{
    component::Navbar,
    page::{
        Api, Batch, Benchmarks, Confseq, Descriptors, Fingerprints, Genindex, Home, Installation,
        Io, JavaScript, Molecule, Protein, PyModindex, Python, Quickstart, SearchPage, Validation,
    },
};

#[derive(Debug, Clone, Routable, PartialEq)]
#[rustfmt::skip]
pub enum Route {
    #[layout(Navbar)]
    #[route("/")]
    Home {},
    #[route("/python")]
    Python {},
    #[route("/javascript")]
    JavaScript {},

    #[route("/benchmarks")]
    Benchmarks {},
    #[route("/validation")]
    Validation {},
    #[route("/api")]
    Api {},
    #[route("/installation")]
    Installation {},
    #[route("/quickstart")]
    Quickstart {},
    #[route("/confseq")]
    Confseq {},
    #[route("/molecule")]
    Molecule {},
    #[route("/batch")]
    Batch {},
    #[route("/fingerprints")]
    Fingerprints {},
    #[route("/descriptors")]
    Descriptors {},
    #[route("/protein")]
    Protein {},
    #[route("/io")]
    Io {},

    #[route("/search")]
    SearchPage {},
    #[route("/genindex")]
    Genindex {},
    #[route("/py-modindex")]
    PyModindex {},
}

#[cfg(all(feature = "ssg", not(target_arch = "wasm32")))]
#[server(endpoint = "static_routes", output = server_fn::codec::Json)]
async fn static_routes() -> Result<Vec<String>, ServerFnError> {
    Ok(Route::static_routes()
        .iter()
        .map(ToString::to_string)
        .collect())
}
