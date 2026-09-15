use dioxus::prelude::*;

use crate::{
    component::{MDI_OPEN_IN_NEW, MdiIcon},
    route::Route,
};

#[component]
pub fn Navbar() -> Element {
    let is_javascript = matches!(use_route::<Route>(), Route::JavaScript {});

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
                    Link {
                        class: if is_javascript { "docs-language-option" } else { "docs-language-option is-active" },
                        to: Route::Python {},
                        aria_current: if is_javascript { None } else { Some("true") },
                        "Python"
                    }
                    Link {
                        class: if is_javascript { "docs-language-option is-active" } else { "docs-language-option" },
                        to: Route::JavaScript {},
                        aria_current: if is_javascript { Some("true") } else { None },
                        title: "JavaScript / WebAssembly documentation — not yet available",
                        "JavaScript"
                    }
                }
                div {
                    class: "docs-top-links",
                    Link { to: if is_javascript { Route::JavaScript {} } else { Route::Python {} }, "Guides" }
                    if !is_javascript {
                        Link { to: Route::Api {}, "API reference" }
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
                a { class: "text-[#7ab5ff] no-underline hover:text-white", href: "https://github.com/cosmol-studio/COSMolKit", target: "_blank", rel: "noreferrer", "GitHub" }
                a { class: "text-[#7ab5ff] no-underline hover:text-white", href: "https://pypi.org/project/cosmolkit/", target: "_blank", rel: "noreferrer", "Python package" }
                a { class: "text-[#7ab5ff] no-underline hover:text-white", href: "https://crates.io/crates/cosmolkit", target: "_blank", rel: "noreferrer", "crates.io" }
                a { class: "text-[#7ab5ff] no-underline hover:text-white", href: "https://tools.cosmol.org/", target: "_blank", rel: "noreferrer", "Web tools" }
            }
        }
    }
}
