use dioxus::prelude::*;

use crate::{component::Seo, route::Route};

const WEBSITE_JSON_LD: &str = r#"{
  "@context": "https://schema.org",
  "@type": "WebSite",
  "name": "COSMolKit",
  "alternateName": "COSMolKit Documentation",
  "url": "https://kit.cosmol.org/"
}"#;

#[component]
pub fn Home() -> Element {
    rsx! {
        Seo {
            title: "COSMolKit documentation — Guides and API reference",
            description: "Learn COSMolKit with installation instructions, Python guides, API reference, and source-backed validation evidence.",

        }
        document::Script { r#type: "application/ld+json", "{WEBSITE_JSON_LD}" }
        div { class: "docs-home",
            main { role: "main", class: "docs-home-inner",
                section { class: "docs-intro", aria_label: "Documentation overview",
                    div {
                        p { class: "home-kicker", "THE COSMOLKIT HANDBOOK" }
                        h1 { "Documentation" }
                        p { class: "docs-intro-description", "A practical guide to the Rust-native chemistry toolkit. Learn the Python API, build molecular workflows, and find the details you need." }
                        div { class: "docs-intro-actions",
                            Link { class: "docs-primary-action", to: Route::Quickstart { fragment: String::new() }, "Start the quick guide", span { aria_hidden: "true", "→" } }
                            Link { class: "docs-text-action", to: Route::Api { fragment: String::new() }, "Browse API reference", span { aria_hidden: "true", "→" } }
                        }
                    }
                    aside { class: "docs-reading-path", aria_label: "Getting started",
                        p { class: "home-kicker", "NEW TO COSMOLKIT?" }
                        StartLink { number: "01", title: "Install the package", detail: "Set up your Python environment", to: Route::Installation { fragment: String::new() } }
                        StartLink { number: "02", title: "Create your first molecule", detail: "Follow the quick start", to: Route::Quickstart { fragment: String::new() } }
                        StartLink { number: "03", title: "Understand molecule values", detail: "Learn the core data model", to: Route::Molecule { fragment: String::new() } }
                    }
                }
                section { class: "docs-guide-section", aria_label: "Explore the documentation",
                    div { class: "docs-section-heading",
                        div { p { class: "home-kicker", "GUIDES & REFERENCE" } h2 { "Find your next step" } }
                        Link { class: "docs-text-action", to: Route::Python {}, "All Python documentation →" }
                    }
                    div { class: "docs-topic-grid",
                        TopicGroup { number: "01", title: "Molecules & data", description: "Build, inspect, and exchange molecular structures.",
                            TopicLink { title: "Molecule values", to: Route::Molecule { fragment: String::new() } }
                            TopicLink { title: "File IO and arrays", to: Route::Io { fragment: String::new() } }
                            TopicLink { title: "Batch workflows", to: Route::Batch { fragment: String::new() } }
                        }
                        TopicGroup { number: "02", title: "Chemical workflows", description: "Explore representations, properties, and 3D structures.",
                            TopicLink { title: "Fingerprints", to: Route::Fingerprints { fragment: String::new() } }
                            TopicLink { title: "Molecular descriptors", to: Route::Descriptors { fragment: String::new() } }
                            TopicLink { title: "Conformers · ConfSeq", to: Route::Confseq { fragment: String::new() } }
                            TopicLink { title: "Protein structures", to: Route::Protein { fragment: String::new() } }
                        }
                        TopicGroup { number: "03", title: "Look up the details", description: "Check signatures, locate symbols, and review evidence.",
                            TopicLink { title: "Search documentation", to: Route::SearchPage { q: String::new(), fragment: String::new() } }
                            TopicLink { title: "Python API reference", to: Route::Api { fragment: String::new() } }
                            TopicLink { title: "General index", to: Route::Genindex { fragment: String::new() } }
                            TopicLink { title: "Validation evidence", to: Route::Validation {} }
                        }
                    }
                    div { class: "docs-availability",
                        span { class: "docs-status-label", "NOT YET AVAILABLE" }
                        Link { to: Route::JavaScript {}, "JavaScript / WebAssembly documentation" }
                        span { aria_hidden: "true", "·" }
                        Link { to: Route::Benchmarks {}, "Benchmark reports" }
                    }
                }
                section { class: "docs-related", aria_label: "Related COSMolKit websites",
                    div { class: "docs-section-heading",
                        div { p { class: "home-kicker", "BEYOND THE DOCUMENTATION" } h2 { "Explore the wider project" } }
                        p { class: "docs-related-note", "Companion resources on tools.cosmol.org" }
                    }
                    div { class: "docs-related-grid",
                        a { class: "docs-external-card docs-external-tools", href: "https://tools.cosmol.org/tools", target: "_blank", rel: "noreferrer",
                            div { class: "docs-external-meta", span { "INTERACTIVE WORKSPACE" } span { "EXTERNAL ↗" } }
                            h3 { "COSMolKit Web tools" }
                            p { "Convert molecular formats, render structures, and explore chemistry directly in your browser. A separate workspace for hands-on tasks." }
                            div { class: "docs-external-destination", span { "Open Web tools" } span { "tools.cosmol.org/tools ↗" } }
                        }
                        a { class: "docs-external-card docs-external-blog", href: "https://tools.cosmol.org/blog", target: "_blank", rel: "noreferrer",
                            div { class: "docs-external-meta", span { "ARTICLES & PERSPECTIVES" } span { "EXTERNAL ↗" } }
                            h3 { "From the blog" }
                            p { "Read beyond the reference: articles about cheminformatics, the toolkit, and the ideas behind the project." }
                            div { class: "docs-external-destination", span { "Read the blog" } span { "tools.cosmol.org/blog ↗" } }
                        }
                    }
                    p { class: "docs-external-disclosure", "These links open the companion website in a new tab. Guides and API reference stay here." }
                }
            }
        }
    }
}

#[component]
fn StartLink(
    number: &'static str,
    title: &'static str,
    detail: &'static str,
    to: Route,
) -> Element {
    rsx! {
        Link { class: "docs-start-link", to,
            span { class: "docs-step-number", "{number}" }
            div { strong { "{title}" } span { "{detail}" } }
            span { class: "docs-link-arrow", aria_hidden: "true", "→" }
        }
    }
}

#[component]
fn TopicGroup(
    number: &'static str,
    title: &'static str,
    description: &'static str,
    children: Element,
) -> Element {
    rsx! {
        div { class: "docs-topic-group",
            span { class: "docs-topic-number", "{number}" }
            h3 { "{title}" }
            p { "{description}" }
            div { class: "docs-topic-links", {children} }
        }
    }
}

#[component]
fn TopicLink(title: &'static str, to: Route) -> Element {
    rsx! {
        Link { class: "docs-topic-link", to, "{title}", span { aria_hidden: "true", "→" } }
    }
}
