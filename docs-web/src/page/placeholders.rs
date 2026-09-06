use dioxus::prelude::*;

use crate::{component::Seo, route::Route};

#[component]
pub fn JavaScript() -> Element {
    rsx! { Placeholder { eyebrow: "JAVASCRIPT API", title: "JavaScript documentation is reserved", body: "This route is reserved for the future JavaScript and WebAssembly API. The current browser tools remain available at tools.cosmol.org.", canonical: "https://kit.cosmol.org/javascript.html" } }
}

#[component]
pub fn JavaScriptHtml() -> Element {
    JavaScript()
}

#[component]
pub fn Benchmarks() -> Element {
    rsx! { Placeholder { eyebrow: "BENCHMARKS", title: "Benchmark reports are reserved", body: "Reproducible benchmark reports will be published here as browser and native measurements are promoted into the maintained documentation surface.", canonical: "https://kit.cosmol.org/benchmarks" } }
}

#[component]
pub fn Validation() -> Element {
    rsx! {
        Seo { title: "COSMolKit Validation | ChEMBL 37 and RDKit parity", description: "COSMolKit validation scope, exact parity comparisons, ChEMBL 37 corpus coverage, and reproducible evidence.", canonical: "https://kit.cosmol.org/validation" }
        div { class: "min-h-screen uu-backdrop m-0 pt-[74px]",
            main { class: "mx-auto w-full max-w-4xl px-6 py-10 font-sans text-[#e8edf5] max-[640px]:px-3.5 max-[640px]:py-7",
                span { class: "text-xs font-bold tracking-[0.08em] text-[#4b96ff]", "VALIDATION" }
                h1 { class: "mb-4 mt-3 text-[28px] leading-[1.35] font-bold text-white", "Source-backed parity evidence" }
                p { class: "max-w-[760px] text-[16px] leading-7 text-[#aebacd]", "The maintained validation record covers exact source comparisons across chemistry, molecular state, fingerprints, descriptors, stereochemistry, coordinates, file IO, batch execution, and concurrency." }
                div { class: "mt-8 grid grid-cols-3 gap-4 max-[760px]:grid-cols-1",
                    Stat { value: "ChEMBL 37", label: "Primary corpus" }
                    Stat { value: "RDKit 2026.03.1", label: "Pinned reference" }
                    Stat { value: "Exact + bounded", label: "Comparison policy" }
                }
                section { class: "mt-10 border-t border-white/10 pt-8",
                    h2 { class: "text-xl font-bold text-white", "Read the maintained record" }
                    p { class: "mt-3 text-[14px] leading-6 text-[#9caabd]", "The full validation scope, evidence tables, source boundaries, and reproduction commands are maintained in the repository's VALIDATION.md and dev/gap_reports documents." }
                    a { class: "mt-5 inline-flex rounded-md bg-[#3082ff] px-4 py-2.5 text-sm font-bold text-white no-underline hover:bg-[#438ee9]", href: "https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md", target: "_blank", rel: "noreferrer", "Open VALIDATION.md" }
                }
            }
        }
    }
}

#[component]
fn Placeholder(
    eyebrow: &'static str,
    title: &'static str,
    body: &'static str,
    canonical: &'static str,
) -> Element {
    rsx! {
        Seo { title: "{title} | COSMolKit", description: "{body}", canonical: "{canonical}" }
        div { class: "min-h-screen uu-backdrop m-0 pt-[74px]",
            main { class: "mx-auto w-full max-w-4xl px-6 py-10 font-sans text-[#e8edf5] max-[640px]:px-3.5 max-[640px]:py-7",
                Link { class: "text-[13px] font-semibold text-[#7ab5ff] no-underline hover:text-[#b4d6ff]", to: Route::Home {}, "Back to documentation" }
                span { class: "mt-8 block text-xs font-bold tracking-[0.08em] text-[#4b96ff]", "{eyebrow}" }
                h1 { class: "mb-4 mt-3 text-[28px] leading-[1.35] font-bold text-white", "{title}" }
                p { class: "max-w-[680px] text-[16px] leading-7 text-[#aebacd]", "{body}" }
            }
        }
    }
}

#[component]
fn Stat(value: &'static str, label: &'static str) -> Element {
    rsx! { div { class: "rounded-md border border-[#28415f] bg-[#0b1727] p-4", strong { class: "block text-sm font-bold text-white", "{value}" }, span { class: "mt-1 block text-xs text-[#8495aa]", "{label}" } } }
}
