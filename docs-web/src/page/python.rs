use dioxus::prelude::*;

use crate::{component::Seo, route::Route};

include!(concat!(env!("OUT_DIR"), "/sphinx_docs.rs"));

#[component]
pub fn Python() -> Element {
    rsx! {
        Seo { title: "Python API Documentation | COSMolKit", description: "COSMolKit Python documentation generated from the project's Sphinx sources.", canonical: "https://kit.cosmol.org/python" }
        div { class: "min-h-screen uu-backdrop m-0 pt-[74px]",
            main { class: "mx-auto w-full max-w-6xl px-6 py-10 font-sans text-[#e8edf5] max-[640px]:px-3.5 max-[640px]:py-7",
                div { class: "border-b border-white/10 pb-7",
                    span { class: "text-xs font-bold tracking-[0.08em] text-[#4b96ff]", "PYTHON API" }
                    h1 { class: "mb-3 mt-2 text-[28px] leading-[1.35] font-bold text-white", "Python documentation" }
                    p { class: "m-0 max-w-[760px] text-[15px] leading-6 text-[#9caabd]", "The reference pages below are compiled from python/docs with Sphinx and presented inside the COSMolKit documentation shell." }
                }
                section { class: "mt-8 grid grid-cols-2 gap-4 max-[760px]:grid-cols-1", aria_label: "Python documentation sections",
                    PythonCard { to: Route::Installation {}, title: "Installation", summary: "Install COSMolKit and verify the Python package." }
                    PythonCard { to: Route::Quickstart {}, title: "Quick Start", summary: "Build molecules and run the first core workflows." }
                    PythonCard { to: Route::Confseq {}, title: "ConfSeq", summary: "Generate and inspect reproducible conformers." }
                    PythonCard { to: Route::Molecule {}, title: "Molecule Values", summary: "Work with molecular graphs, stereochemistry, and coordinates." }
                    PythonCard { to: Route::Batch {}, title: "Batch Workflows", summary: "Process ordered molecule collections with scalar parity." }
                    PythonCard { to: Route::Fingerprints {}, title: "Fingerprints", summary: "Compute source-backed molecular fingerprint vectors." }
                    PythonCard { to: Route::Descriptors {}, title: "Molecular Descriptors", summary: "Calculate descriptors and composition properties." }
                    PythonCard { to: Route::Protein {}, title: "Protein Structures", summary: "Read and inspect structural-biology values." }
                    PythonCard { to: Route::Io {}, title: "File IO and Arrays", summary: "Use molecular file readers, writers, and NumPy arrays." }
                    PythonCard { to: Route::Api {}, title: "API Reference", summary: "Browse the complete generated Python API reference." }
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

macro_rules! sphinx_page {
    ($name:ident, $constant:ident, $source:literal, $title:literal, $canonical:literal) => {
        #[component]
        pub fn $name() -> Element {
            rsx! {
                document::Style { "{STYLESHEET}" }
                Seo { title: $title.to_string(), description: format!("{title} — COSMolKit Python documentation", title = $title), canonical: $canonical.to_string() }
                div { class: "min-h-screen uu-backdrop m-0 pt-[74px]",
                    main { class: "mx-auto w-full max-w-6xl px-6 py-10 font-sans text-[#e8edf5] max-[640px]:px-3.5 max-[640px]:py-7",
                        div { class: "mb-7 flex items-center justify-between gap-4 border-b border-white/10 pb-5 max-[640px]:items-start max-[640px]:flex-col",
                            Link { class: "text-[13px] font-semibold text-[#7ab5ff] no-underline hover:text-[#b4d6ff]", to: Route::Python {}, "Python documentation" }
                            a { class: "text-[12px] text-[#718299] no-underline hover:text-white", href: concat!("https://github.com/cosmol-studio/COSMolKit/blob/main/python/docs/source/", $source, ".rst"), target: "_blank", rel: "noreferrer", "View source" }
                        }
                        div { class: "sphinx-content w-full max-w-6xl rounded-lg border border-[#28415f] bg-[#0b1727] p-7 max-[640px]:p-5", dangerous_inner_html: $constant, }
                    }
                }
            }
        }
    };
}

sphinx_page!(
    Installation,
    INSTALLATION,
    "installation",
    "Installation | COSMolKit",
    "https://kit.cosmol.org/installation.html"
);
sphinx_page!(
    Quickstart,
    QUICKSTART,
    "quickstart",
    "Quick Start | COSMolKit",
    "https://kit.cosmol.org/quickstart.html"
);
sphinx_page!(
    Confseq,
    CONFSEQ,
    "confseq",
    "ConfSeq | COSMolKit",
    "https://kit.cosmol.org/confseq.html"
);
sphinx_page!(
    Molecule,
    MOLECULE,
    "molecule",
    "Molecule Values | COSMolKit",
    "https://kit.cosmol.org/molecule.html"
);
sphinx_page!(
    Batch,
    BATCH,
    "batch",
    "Batch Workflows | COSMolKit",
    "https://kit.cosmol.org/batch.html"
);
sphinx_page!(
    Fingerprints,
    FINGERPRINTS,
    "fingerprints",
    "Fingerprints | COSMolKit",
    "https://kit.cosmol.org/fingerprints.html"
);
sphinx_page!(
    Descriptors,
    DESCRIPTORS,
    "descriptors",
    "Molecular Descriptors | COSMolKit",
    "https://kit.cosmol.org/descriptors.html"
);
sphinx_page!(
    Protein,
    PROTEIN,
    "protein",
    "Protein Structures | COSMolKit",
    "https://kit.cosmol.org/protein.html"
);
sphinx_page!(
    Io,
    IO,
    "io",
    "File IO and Arrays | COSMolKit",
    "https://kit.cosmol.org/io.html"
);
sphinx_page!(
    Api,
    API,
    "api",
    "API Reference | COSMolKit",
    "https://kit.cosmol.org/api.html"
);
sphinx_page!(
    SearchPage,
    SEARCH,
    "search",
    "Search | COSMolKit",
    "https://kit.cosmol.org/search.html"
);
sphinx_page!(
    Genindex,
    GENINDEX,
    "genindex",
    "General Index | COSMolKit",
    "https://kit.cosmol.org/genindex.html"
);
sphinx_page!(
    PyModindex,
    PY_MODINDEX,
    "py-modindex",
    "Python Module Index | COSMolKit",
    "https://kit.cosmol.org/py-modindex.html"
);
