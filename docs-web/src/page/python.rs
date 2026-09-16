use dioxus::prelude::*;

use crate::{component::Seo, route::Route};

include!(concat!(env!("OUT_DIR"), "/sphinx_docs.rs"));
include!(concat!(env!("OUT_DIR"), "/search_asset.rs"));

#[component]
pub fn Python() -> Element {
    rsx! {
        Seo { title: "Python API Documentation | COSMolKit", description: "Explore COSMolKit Python guides for molecules, fingerprints, descriptors, conformers, structural biology, and file IO, with a complete API reference." }
        div { class: "min-h-screen uu-backdrop m-0 pt-[74px]",
            main { role: "main", class: "mx-auto w-full max-w-6xl px-6 py-10 font-sans text-[#e8edf5] max-[640px]:px-3.5 max-[640px]:py-7",
                div { class: "border-b border-white/10 pb-7",
                    span { class: "text-xs font-bold tracking-[0.08em] text-[#4b96ff]", "PYTHON API" }
                    h1 { class: "mb-3 mt-2 text-[28px] leading-[1.35] font-bold text-white", "Python documentation" }
                    p { class: "m-0 max-w-[760px] text-[15px] leading-6 text-[#9caabd]", "Explore Python guides, examples, and the complete COSMolKit API reference." }
                    Link { class: "docs-text-action", to: Route::SearchPage { q: String::new(), fragment: String::new() }, "Search Python documentation" }
                }
                section { class: "mt-8 grid grid-cols-2 gap-4 max-[760px]:grid-cols-1", aria_label: "Python documentation sections",
                    for route in Route::navigation("python").into_iter().filter(|r| !r.metadata().summary.is_empty()) {
                        PythonCard { title: route.metadata().label, summary: route.metadata().summary, to: route }
                    }
                }
            }
        }
    }
}

#[component]
fn PythonCard(to: Route, title: &'static str, summary: &'static str) -> Element {
    rsx! {
        Link { to,
            class: "group flex min-h-[148px] flex-col rounded-lg border border-[#28415f] bg-[#0b1727] p-5 no-underline transition-colors hover:border-[#438ee9] hover:bg-[#0d1b2d]",
            h2 { class: "m-0 text-lg font-bold text-white", "{title}" }
            p { class: "mt-3 text-[13px] leading-5 text-[#9caabd]", "{summary}" }
            div { class: "mt-auto flex items-center justify-between border-t border-white/8 pt-4 text-xs font-semibold text-[#7ab5ff]",
                span { "Open section" }
                span { class: "text-base transition-transform group-hover:translate-x-1", ">" }
            }
        }
    }
}

#[component]
fn DocumentationNavigation(current: &'static str) -> Element {
    let pages = Route::navigation("python");
    rsx! {
        nav { class: "docs-navigation", aria_label: "Python documentation",
            Link { class: "docs-sidebar-brand", to: Route::Python {},
                span { class: "docs-eyebrow", "COSMOLKIT" }
                span { "Python documentation" }
            }
            p { class: "docs-nav-label", "USER GUIDE & REFERENCE" }
            for route in pages {
                Link {
                    to: route.clone(),
                    class: if current == route.metadata().docname { "docs-nav-link is-active" } else { "docs-nav-link" },
                    "{route.metadata().label}"
                }
            }
        }
    }
}

macro_rules! sphinx_page {
    ($name:ident, $constant:ident, $source:literal, $title:literal $(, $parameter:ident)*) => {
        #[component]
        pub fn $name($($parameter: String),*) -> Element {
            $(let _ = &$parameter;)*
            #[cfg(all(target_arch = "wasm32", feature = "web"))]
            super::anchor::use_fragment_scroll();
            let description = sphinx_metadata($source);
            rsx! {
                document::Style { "{STYLESHEET}" }
                Seo { title: $title.to_string(), description,  }
                div { class: "min-h-screen uu-backdrop m-0 pt-[74px]",
                    div { class: "docs-layout",
                        a { class: "docs-skip-link", href: "#docs-article", "Skip to content" }
                        aside { class: "docs-sidebar",
                            DocumentationNavigation { current: $source }
                        }
                        details { class: "docs-mobile-navigation",
                            summary { "Browse Python documentation" }
                            DocumentationNavigation { current: $source }
                        }
                        main { role: "main", class: "docs-main", id: "docs-article", tabindex: "-1",
                            div { class: "docs-article-toolbar",
                                Link { to: Route::Python {}, "Python documentation" }
                                if !matches!($source, "search" | "genindex" | "py-modindex") {
                                    a { href: concat!("https://github.com/cosmol-studio/COSMolKit/blob/main/python/docs/source/", $source, ".rst"), target: "_blank", rel: "noreferrer", "View source ↗" }
                                }
                            }
                            if !sphinx_toc($source).is_empty() {
                                details { class: "docs-mobile-toc",
                                    summary { "On this page" }
                                    nav { class: "docs-toc-tree", aria_label: "On this page", dangerous_inner_html: sphinx_toc($source) }
                                }
                            }
                            article { class: "sphinx-content docs-article",
                                if $source == "search" {
                                    h1 { "Search documentation" }
                                    form { id: "docs-search-form", action: Route::SearchPage { q: String::new(), fragment: String::new() }.path(), method: "get", role: "search",
                                        label { r#for: "docs-query", "Search guides and API reference" }
                                        input { id: "docs-query", name: "q", r#type: "search", placeholder: "Type to search: Molecule, fingerprint, from_smiles…", autocomplete: "off", maxlength: "512" }
                                        button { r#type: "reset", "Clear" }
                                    }
                                    p { id: "docs-search-status", role: "status", aria_live: "polite", "Loading search… If this message remains, reload the page to retry." }
                                    noscript { p { "Enable JavaScript to search, or browse the documentation using the navigation links." } }
                                }
                                div { dangerous_inner_html: $constant }
                            }
                        }
                        if !sphinx_toc($source).is_empty() {
                            aside { class: "docs-toc",
                                p { class: "docs-nav-label", "ON THIS PAGE" }
                                nav { class: "docs-toc-tree", aria_label: "On this page", dangerous_inner_html: sphinx_toc($source) }
                            }
                        }
                    }
                }
                if $source == "search" {
                    SearchScripts {}
                }
            }
        }
    };
}

#[component]
fn SearchScripts() -> Element {
    // wasm-bindgen generates the bindings. This small loader runs only when the
    // search component mounts; module imports and init reuse the loaded engine.
    let loader = format!(
        r#"try {{
            const {{ default: init, mount_search }} = await import("{SEARCH_BINDINGS}");
            await init({{ module_or_path: "{SEARCH_WASM}" }});
            mount_search();
        }} catch (error) {{
            const status = document.getElementById("docs-search-status");
            if (status) status.textContent = "Search could not load. Reload the page to retry.";
            console.error(error);
        }}"#
    );
    rsx! {
        document::Script {
            id: "docs-search-loader", r#type: "module",
            "data-search-bindings": SEARCH_BINDINGS,
            "data-search-wasm": SEARCH_WASM,
            "{loader}"
        }
    }
}

include!(concat!(env!("OUT_DIR"), "/sphinx_pages.rs"));
