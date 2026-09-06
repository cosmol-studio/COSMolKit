use dioxus::prelude::*;

use crate::{
    component::{MDI_OPEN_IN_NEW, MdiIcon},
    route::Route,
};

#[component]
pub fn Navbar() -> Element {
    rsx! {
        header {
            class: "fixed top-0 z-20 w-full px-6 pt-5 font-sans max-[800px]:px-3.5 max-[800px]:pt-3.5",
            nav {
                class: "mx-auto flex min-h-[54px] w-full max-w-6xl items-center justify-between gap-6 rounded-lg border border-white/8 bg-[#0c1828]/90 px-4 shadow-[0_10px_32px_rgba(0,0,0,0.12)] backdrop-blur-md max-[800px]:min-h-[58px] max-[800px]:gap-3.5 max-[800px]:px-[13px]",
                Link {
                    class: "flex min-w-0 items-baseline gap-[5px] whitespace-nowrap text-[17px] font-extrabold text-white no-underline max-[800px]:block max-[800px]:text-[15px] max-[800px]:leading-[1.2]",
                    to: Route::Home {},
                    span { class: "text-[#4b96ff] max-[800px]:block", "COSMolkit" }
                }
                div {
                    class: "flex items-center gap-6 max-[800px]:gap-3.5 max-[480px]:gap-[11px]",
                    a { class: "uu-nav-link whitespace-nowrap text-sm font-semibold text-[#9baabd] no-underline hover:text-white focus-visible:text-white max-[800px]:text-xs", href: "https://tools.cosmol.org/tools", "Tools" }
                    a { class: "uu-nav-link whitespace-nowrap text-sm font-semibold text-[#9baabd] no-underline hover:text-white focus-visible:text-white max-[800px]:text-xs", href: "https://tools.cosmol.org/blog", "Blog" }
                    div { class: "group relative",
                        button { r#type: "button", class: "uu-nav-button inline-flex cursor-pointer items-center gap-[5px] whitespace-nowrap bg-transparent text-sm font-semibold text-[#9baabd] outline-none hover:text-white focus-visible:text-white max-[800px]:text-xs", aria_label: "Open COSMol documentation", "Docs", MdiIcon { size: 15, path: MDI_OPEN_IN_NEW } }
                        div { class: "invisible absolute top-full left-0 w-[300px] translate-y-1 pt-3 opacity-0 transition-all group-hover:visible group-hover:translate-y-0 group-hover:opacity-100 group-focus-within:visible group-focus-within:translate-y-0 group-focus-within:opacity-100 max-[560px]:fixed max-[560px]:top-[67px] max-[560px]:right-3.5 max-[560px]:left-auto max-[560px]:w-[min(300px,calc(100vw-28px))]",
                            div { class: "rounded-md border border-[#2a3b52] bg-[#0c1828] p-2 shadow-[0_16px_38px_rgba(0,0,0,0.32)]",
                                span { class: "block px-2.5 pb-2 pt-1 text-[10px] font-bold tracking-[0.08em] text-[#718299]", "COSMOL DOCUMENTATION" }
                                a { class: "block rounded-[5px] px-2.5 py-2.5 text-left no-underline hover:bg-[#13253b]", href: "https://kit.cosmol.org/", target: "_blank", rel: "noreferrer", span { class: "block text-xs font-bold text-white", "COSMolKit" }, span { class: "mt-0.5 block text-[11px] text-[#8495aa]", "Python API and toolkit documentation" } }
                                a { class: "block rounded-[5px] px-2.5 py-2.5 text-left no-underline hover:bg-[#13253b]", href: "/javascript.html", span { class: "block text-xs font-bold text-white", "JavaScript API" }, span { class: "mt-0.5 block text-[11px] text-[#8495aa]", "Reserved WebAssembly documentation" } }
                            }
                        }
                    }
                    a { class: "uu-nav-link whitespace-nowrap text-sm font-semibold text-[#9baabd] no-underline hover:text-white focus-visible:text-white max-[800px]:text-xs max-[640px]:hidden", href: "https://tools.cosmol.org/ecosystem", "Ecosystem" }
                    div { class: "group relative",
                        button { r#type: "button", class: "uu-nav-button inline-flex cursor-pointer items-center gap-[5px] whitespace-nowrap bg-transparent text-sm font-semibold text-[#9baabd] outline-none hover:text-white focus-visible:text-white max-[800px]:text-xs", aria_label: "Open COSMol GitHub projects", "Github", MdiIcon { size: 15, path: MDI_OPEN_IN_NEW } }
                        div { class: "invisible absolute top-full right-0 w-[300px] translate-y-1 pt-3 opacity-0 transition-all group-hover:visible group-hover:translate-y-0 group-hover:opacity-100 group-focus-within:visible group-focus-within:translate-y-0 group-focus-within:opacity-100 max-[560px]:fixed max-[560px]:top-[67px] max-[560px]:right-3.5 max-[560px]:w-[min(300px,calc(100vw-28px))]",
                            div { class: "rounded-md border border-[#2a3b52] bg-[#0c1828] p-2 shadow-[0_16px_38px_rgba(0,0,0,0.32)]",
                                span { class: "block px-2.5 pb-2 pt-1 text-[10px] font-bold tracking-[0.08em] text-[#718299]", "COSMOL GITHUB" }
                                a { class: "block rounded-[5px] px-2.5 py-2.5 text-left no-underline hover:bg-[#13253b]", href: "https://github.com/cosmol-studio/COSMolKit", target: "_blank", rel: "noreferrer", span { class: "block text-xs font-bold text-white", "COSMolKit" }, span { class: "mt-0.5 block text-[11px] text-[#8495aa]", "Core chemistry toolkit" } }
                            }
                        }
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
                span { class: "font-semibold text-[#9caabd]", "Rust cheminformatics powered by COSMolKit" }
                a { class: "text-[#7ab5ff] no-underline hover:text-white", href: "https://github.com/cosmol-studio/COSMolKit", target: "_blank", rel: "noreferrer", "GitHub" }
                a { class: "text-[#7ab5ff] no-underline hover:text-white", href: "https://pypi.org/project/cosmolkit/", target: "_blank", rel: "noreferrer", "Python package" }
                a { class: "text-[#7ab5ff] no-underline hover:text-white", href: "https://crates.io/crates/cosmolkit", target: "_blank", rel: "noreferrer", "crates.io" }
                a { class: "text-[#7ab5ff] no-underline hover:text-white", href: "https://tools.cosmol.org/", target: "_blank", rel: "noreferrer", "Web tools" }
            }
        }
    }
}
