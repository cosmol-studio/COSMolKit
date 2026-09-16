use dioxus::prelude::*;

use crate::{
    component::{MDI_OPEN_IN_NEW, MdiIcon},
    route::Route,
};

#[component]
pub fn Navbar() -> Element {
    let current = use_route::<Route>();
    let is_javascript = current.metadata().binding == "javascript";

    rsx! {
        header {
            class: "docs-masthead",
            nav {
                class: "docs-masthead-inner",
                aria_label: "Main navigation",
                Link {
                    class: "docs-brand",
                    to: Route::Home {},
                    "COSMolKit"
                    span { class: "docs-brand-badge", "DOCS" }
                }
                div { class: "docs-language-switch", role: "group", aria_label: "Documentation language",
                    for (binding, label) in [("python", "Python"), ("javascript", "JavaScript")] {
                        {
                            let counterpart = current.counterpart(binding);
                            let missing_topic = counterpart.is_none();
                            let target = counterpart.unwrap_or_else(|| if binding == "python" { Route::Python {} } else { Route::JavaScript {} });
                            rsx! { a {
                                class: if current.metadata().binding == binding { "docs-language-option is-active" } else { "docs-language-option" },
                                href: target.path(),
                                aria_current: if current.metadata().binding == binding { Some("page") } else { None },
                                title: if missing_topic { Some("Open this language's documentation overview; this topic is not yet available") } else { None },
                                "{label}"
                            } }
                        }
                    }
                }
                div {
                    class: "docs-top-links",
                    if !is_javascript {
                        Link { to: Route::SearchPage { q: String::new(), fragment: String::new() }, "Search" }
                    }
                    Link { to: if is_javascript { Route::JavaScript {} } else { Route::Python {} }, "Guides" }
                    if !is_javascript {
                        Link { to: Route::Api { fragment: String::new() }, "API reference" }
                    }
                    Link { class: "docs-top-secondary", to: Route::Validation {}, "Validation" }
                    a {
                        class: "docs-top-secondary",
                        href: "https://github.com/cosmol-studio/COSMolKit",
                        target: "_blank",
                        rel: "noreferrer",
                        "GitHub"
                        MdiIcon { size: 13, path: MDI_OPEN_IN_NEW }
                    }
                }
            }
        }
        SuspenseBoundary {
            fallback: |_| rsx! {
                main { class: "grid min-h-[calc(100vh-74px)] place-items-center bg-[#071426] px-6 pt-[74px] font-sans text-[#e8edf5]", span { class: "text-sm font-semibold text-[#91a1b5]", "Loading documentation" } }
            },
            Outlet::<Route> {}
        }
        footer {
            class: "border-t border-white/8 bg-[#081321] px-6 py-5 font-sans text-xs text-[#718299] max-[640px]:px-3.5",
            div {
                class: "mx-auto flex w-full max-w-6xl flex-wrap items-center gap-x-5 gap-y-2",
                span { class: "font-semibold text-[#9caabd]", "COSMolKit documentation" }
                nav { class: "cosmolkit-project-links flex flex-wrap gap-x-5 gap-y-2", aria_label: "Project and Rust crates",
                    for (label, href) in [
                        ("GitHub", "https://github.com/cosmol-studio/COSMolKit"),
                        ("Python package", "https://pypi.org/project/cosmolkit/"),
                        ("cosmolkit", "https://crates.io/crates/cosmolkit"),
                        ("cosmolkit-core", "https://crates.io/crates/cosmolkit-core"),
                        ("cosmolkit-inchi", "https://crates.io/crates/cosmolkit-inchi"),
                        ("cosmolkit-ringdecomposer", "https://crates.io/crates/cosmolkit-ringdecomposer"),
                        ("Rust API", "https://docs.rs/cosmolkit/latest/cosmolkit/"),
                        ("Documentation", "https://kit.cosmol.org/"),
                        ("Web tools", "https://tools.cosmol.org/"),
                    ] {
                        a { class: "text-[#7ab5ff] no-underline hover:text-white", href, target: "_blank", rel: "noreferrer", "{label}" }
                    }
                }
            }
        }
    }
}
