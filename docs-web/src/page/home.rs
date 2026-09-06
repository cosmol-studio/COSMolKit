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
            title: "COSMolKit — Rust-native cheminformatics documentation",
            description: "COSMolKit is a pure-Rust cheminformatics and structural biology toolkit with Python bindings, source-backed RDKit parity, and documented validation evidence.",
            canonical: "https://kit.cosmol.org/",
        }
        document::Script { r#type: "application/ld+json", "{WEBSITE_JSON_LD}" }
        div { class: "min-h-screen uu-backdrop m-0 pt-[74px]",
            main { class: "mx-auto w-full max-w-6xl px-6 py-5 font-sans text-[#e8edf5] max-[640px]:px-3.5",
                section { class: "home-hero relative mx-auto mb-4 flex max-w-5xl items-center gap-8 border-b border-white/10 pb-6 max-[800px]:flex-col max-[800px]:items-start max-[800px]:gap-5",
                    div {
                        class: "home-hero-copy z-10 flex flex-col items-start justify-center",
                        span { class: "text-xs font-bold tracking-[0.08em] text-[#4b96ff]", "PURE-RUST CHEMINFORMATICS" }
                        h1 { class: "mb-4 mt-3 text-[28px] leading-[1.35] font-bold text-white", "COSMolKit documentation" }
                        p { class: "mb-4 max-w-[500px] leading-relaxed text-white opacity-75", "A source-backed Rust toolkit for molecular graphs, structural biology, fingerprints, descriptors, stereochemistry, conformers, and high-throughput workflows." }
                        div { class: "flex w-full justify-start gap-4",
                            Link { class: "cursor-pointer rounded bg-[#3082ff] px-4 py-2 text-white", to: Route::Python {}, "Open Python docs" }
                            Link { class: "rounded border border-white px-4 py-2 text-white", to: Route::Validation {}, "Read validation" }
                        }
                    }
                    div { class: "home-hero-art flex items-center justify-center overflow-visible",
                        img { class: "home-logo opacity-90", src: asset!("/assets/logo.svg"), alt: "COSMolKit molecular logo" }
                    }
                }
                section { class: "mt-8 grid grid-cols-2 gap-4 max-[760px]:grid-cols-1", aria_label: "Documentation sections",
                    DocCard { eyebrow: "PYTHON API", title: "Python documentation", summary: "Install the package, learn the Molecule value API, use batch workflows, fingerprints, descriptors, protein structures, and file IO.", to: Route::Python {} }
                    DocCard { eyebrow: "JAVASCRIPT", title: "JavaScript documentation", summary: "A reserved documentation surface for a future JavaScript and WebAssembly API.", to: Route::JavaScriptHtml {} }
                    DocCard { eyebrow: "VALIDATION", title: "Validation evidence", summary: "Read the maintained parity boundary, ChEMBL 37 corpus scope, exact comparison policy, and reproducibility notes.", to: Route::Validation {} }
                    DocCard { eyebrow: "BENCHMARKS", title: "Benchmark reports", summary: "A reserved home for reproducible performance measurements across core chemistry and browser workflows.", to: Route::Benchmarks {} }
                }
            }
        }
    }
}

#[component]
pub fn HomeHtml() -> Element {
    Home()
}

#[component]
fn DocCard(
    eyebrow: &'static str,
    title: &'static str,
    summary: &'static str,
    to: Route,
) -> Element {
    rsx! {
        Link { to: to, class: "group flex min-h-[185px] flex-col rounded-lg border border-[#28415f] bg-[#0b1727] p-5 no-underline transition-colors hover:border-[#438ee9] hover:bg-[#0d1b2d]",
            span { class: "text-[10px] font-bold tracking-[0.08em] text-[#72bddb]", "{eyebrow}" }
            h2 { class: "mb-2 mt-3 text-xl font-bold text-white", "{title}" }
            p { class: "m-0 text-[13px] leading-5 text-[#9caabd]", "{summary}" }
            div { class: "mt-auto flex items-center justify-between border-t border-white/8 pt-4 text-xs font-semibold text-[#7ab5ff]", span { "Open section" }, span { class: "text-base transition-transform group-hover:translate-x-1", ">" } }
        }
    }
}
