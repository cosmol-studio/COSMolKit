mod component;
mod page;
mod route;

use dioxus::prelude::*;
use route::Route;

#[cfg_attr(all(target_arch = "wasm32", feature = "ssg"), allow(dead_code))]
const FAVICON: Asset = asset!("/assets/logo.svg");

#[cfg(all(feature = "ssg", not(target_arch = "wasm32")))]
fn main() {
    dioxus::LaunchBuilder::new()
        .with_cfg(server_only! {
            dioxus::server::ServeConfig::builder()
                .incremental(
                    dioxus::server::IncrementalRendererConfig::new()
                        .static_dir(
                            std::env::current_exe()
                                .expect("current executable path")
                                .parent()
                                .expect("server executable directory")
                                .join("public"),
                        )
                        // SSG output is a deployable artifact; never reuse a
                        // route from an older build when source or head data changed.
                        .clear_cache(true),
                )
        })
        .launch(App);
}

#[cfg(any(not(feature = "ssg"), target_arch = "wasm32"))]
fn main() {
    dioxus::launch(App);
}

#[component]
fn App() -> Element {
    rsx! {
        span { id: "wasm-ready", hidden: true }
        InitialLoader {}
        document::Link { rel: "icon", href: FAVICON }
        document::Link { rel: "stylesheet", href: asset!("/assets/main.css") }
        document::Link { rel: "stylesheet", href: asset!("/assets/tailwind.css") }
        document::Link { rel: "stylesheet", href: asset!("/assets/docs.css") }
        Router::<Route> {}
    }
}

#[component]
fn InitialLoader() -> Element {
    #[allow(unused_mut)]
    let mut visible = use_signal(|| cfg!(target_arch = "wasm32"));

    #[cfg(target_arch = "wasm32")]
    use_effect(move || visible.set(false));

    if !visible() {
        return rsx! {};
    }

    rsx! {
        div {
            id: "initial-loader",
            role: "status",
            aria_live: "polite",
            aria_label: "Loading COSMolKit documentation",
            div { class: "initial-loader-content",
                p { class: "initial-loader-brand", span { "COSMolKit" } " Docs" }
                p { class: "initial-loader-status", "Loading documentation" }
                span { class: "initial-loader-progress", aria_hidden: "true" }
            }
        }
    }
}
