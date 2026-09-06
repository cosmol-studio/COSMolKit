use dioxus::prelude::*;

use crate::{
    component::Navbar,
    page::{
        Api, Batch, Benchmarks, Confseq, Descriptors, Fingerprints, Genindex, Home, HomeHtml,
        Installation, Io, JavaScript, JavaScriptHtml, Molecule, Protein, PyModindex, Python,
        Quickstart, SearchPage, Validation,
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
    #[route("/javascript.html")]
    JavaScriptHtml {},
    #[route("/benchmarks")]
    Benchmarks {},
    #[route("/validation")]
    Validation {},
    #[route("/api.html")]
    Api {},
    #[route("/installation.html")]
    Installation {},
    #[route("/quickstart.html")]
    Quickstart {},
    #[route("/confseq.html")]
    Confseq {},
    #[route("/molecule.html")]
    Molecule {},
    #[route("/batch.html")]
    Batch {},
    #[route("/fingerprints.html")]
    Fingerprints {},
    #[route("/descriptors.html")]
    Descriptors {},
    #[route("/protein.html")]
    Protein {},
    #[route("/io.html")]
    Io {},
    #[route("/index.html")]
    HomeHtml {},
    #[route("/search.html")]
    SearchPage {},
    #[route("/genindex.html")]
    Genindex {},
    #[route("/py-modindex.html")]
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
