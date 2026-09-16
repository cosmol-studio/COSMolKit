use crate::component::Navbar;
use dioxus::prelude::*;
include!(concat!(env!("OUT_DIR"), "/routes.rs"));

#[cfg(all(feature = "ssg", not(target_arch = "wasm32")))]
#[server(endpoint = "static_routes", output = server_fn::codec::Json)]
async fn static_routes() -> Result<Vec<String>, ServerFnError> {
    Ok(Route::static_routes()
        .iter()
        .map(|route| route.path().to_string())
        .collect())
}
