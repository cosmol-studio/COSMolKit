"""Generate Rust route declarations and Sphinx wrappers from routes.toml."""

from pathlib import Path
import sys

from extract_sphinx_metadata import rust_literal as lit
from route_contract import BASE_URL, SOCIAL_IMAGE_URL, PAGES, counterpart
from build_search_bundle import build_search_bundle


def variant(page):
    fields = page['query'] + ': String::new(), ' if page.get('query') else ''
    if page.get('source') == 'sphinx':
        fields += 'fragment: String::new()'
    return f'Self::{page["component"]} {{{fields}}}'


def generate(out):
    lines = ["// Generated from routes.toml; do not edit.", "use crate::page::*;",
             f'pub const SOCIAL_IMAGE_URL: &str = {lit(SOCIAL_IMAGE_URL)};',
             "#[derive(Debug, Clone, Routable, PartialEq)]", "pub enum Route {", "#[layout(Navbar)]"]
    for p in PAGES:
        suffix = "?:" + p["query"] if p.get("query") else ""
        fields = p["query"] + ": String, " if p.get("query") else ""
        if p.get("source") == "sphinx":
            suffix += "#:fragment"
            fields += "fragment: String"
        lines += [f'#[route({lit(p["path"] + suffix)})]', p["component"] + " {" + fields + "},"]
    lines += ["}", "pub struct RouteMetadata { pub binding: &'static str, pub label: &'static str, pub summary: &'static str, pub docname: &'static str, pub indexable: bool }",
              "impl Route { pub fn metadata(&self) -> RouteMetadata { match self {"]
    for p in PAGES:
        fields = {k: lit(p.get(k, "")) for k in ("binding", "label", "summary", "docname")}
        fields.update(indexable=str(p["indexable"]).lower())
        lines.append(f'Self::{p["component"]} {{..}} => RouteMetadata {{' + ",".join(f"{k}: {v}" for k,v in fields.items()) + "},")
    lines += ["}}",
              'pub fn canonical(&self) -> String { format!("{}{}", ' + lit(BASE_URL.rstrip('/')) + ', self.path()) }',
              'pub fn navigation(binding: &str) -> Vec<Self> { vec![' + ",".join(variant(p) for p in sorted(PAGES,key=lambda p:p.get('order',-1)) if 'order' in p and p['status']=='ready') + '].into_iter().filter(|route| route.metadata().binding == binding).collect() }', "}"]
    pairing = ['impl Route { pub fn counterpart(&self, binding: &str) -> Option<Self> { match (self, binding) {']
    for page in PAGES:
        for binding in ("python", "javascript"):
            target = counterpart(PAGES, page, binding)
            if target:
                pairing.append(f'(Self::{page["component"]} {{..}}, {lit(binding)}) => Some({variant(target)}),')
    pairing += ['_ => None, } } }']
    lines.extend(pairing)
    lines.append("impl Route { pub fn path(&self) -> &'static str { match self {" + "".join(f'Self::{p["component"]} {{..}} => {lit(p["path"])},' for p in PAGES) + "} } }")
    out.mkdir(parents=True, exist_ok=True)
    build_search_bundle(out)
    (out / "routes.rs").write_text("\n".join(lines), encoding="utf-8")
    wrappers = []
    exports = []
    for p in PAGES:
        if p.get("source") == "sphinx":
            wrappers.append(f'sphinx_page!({p["component"]}, {p["docname"].upper().replace("-", "_").replace("/", "_")}, {lit(p["docname"])}, {lit(p["title"])}{", " + p["query"] if p.get("query") else ""}, fragment);')
            exports.append(p['component'])
    (out / "sphinx_pages.rs").write_text("\n".join(wrappers), encoding="utf-8")
    (out / "sphinx_exports.rs").write_text('pub use python::{Python,' + ','.join(exports) + '};', encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) == 2:
        generate(Path(sys.argv[1]))
    # Small line protocol consumed by build.rs. Validated paths cannot contain tabs.
    for page in PAGES:
        print("\t".join((page.get("docname", ""), page["path"], page.get("source", "native"))))
