use std::{
    env, fs,
    path::{Path, PathBuf},
    process::Command,
};

struct Page {
    docname: String,
    route: String,
    source: String,
}

fn contract() -> &'static Vec<Page> {
    static CONTRACT: std::sync::OnceLock<Vec<Page>> = std::sync::OnceLock::new();
    CONTRACT.get_or_init(|| {
        let manifest = env::var_os("CARGO_MANIFEST_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|| {
                if Path::new("docs-web/routes.toml").is_file() {
                    PathBuf::from("docs-web")
                } else {
                    PathBuf::from(".")
                }
            });
        let manifest = manifest.canonicalize().expect("documentation workspace");
        let root = manifest.parent().expect("repository root");
        let mut command = Command::new(docs_python(root));
        command.arg(manifest.join("scripts/generate_routes.py"));
        if let Some(out) = env::var_os("OUT_DIR") {
            command.arg(out);
        }
        let output = command
            .output()
            .expect("generate documentation route contract");
        assert!(
            output.status.success(),
            "route generation failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8(output.stdout)
            .expect("UTF-8 contract")
            .lines()
            .map(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                assert_eq!(fields.len(), 3, "invalid route generator record");
                Page {
                    docname: fields[0].into(),
                    route: fields[1].into(),
                    source: fields[2].into(),
                }
            })
            .collect()
    })
}

fn pages() -> Vec<&'static str> {
    contract()
        .iter()
        .filter(|p| p.source == "sphinx")
        .map(|p| p.docname.as_str())
        .collect()
}

fn main() {
    println!("cargo:rerun-if-changed=../python/docs/source");
    println!("cargo:rerun-if-changed=../python/docs/build/html");
    println!("cargo:rerun-if-changed=scripts/extract_sphinx_metadata.py");
    println!("cargo:rerun-if-changed=scripts/html_metadata.py");
    println!("cargo:rerun-if-changed=routes.toml");
    println!("cargo:rerun-if-changed=scripts/route_contract.py");
    println!("cargo:rerun-if-changed=scripts/generate_routes.py");
    println!("cargo:rerun-if-changed=scripts/build_search_bundle.py");
    println!("cargo:rerun-if-changed=scripts/generate_search_index.py");
    println!("cargo:rerun-if-changed=src/search");
    println!("cargo:rerun-if-env-changed=COSMOLKIT_WASM_BINDGEN");
    println!("cargo:rerun-if-env-changed=COSMOLKIT_DOCS_PYTHON");

    let manifest_dir =
        PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").expect("manifest directory"));
    let repository_root = manifest_dir.parent().expect("repository root");
    let docs_source = repository_root.join("python/docs/source");
    let docs_output = repository_root.join("python/docs/build/html");

    if !docs_output.join("index.html").is_file() {
        build_sphinx(&repository_root, &docs_source, &docs_output);
    }

    let out_dir = PathBuf::from(env::var_os("OUT_DIR").expect("OUT_DIR"));
    // The library target builds only the search engine. It must not recursively
    // build the website or generate browser bindings for itself.
    if env::var_os("CARGO_FEATURE_SEARCH_ENGINE").is_some() {
        let status = Command::new(docs_python(repository_root))
            .arg(manifest_dir.join("scripts/generate_search_index.py"))
            .arg(&docs_output)
            .arg(out_dir.join("search_index.json"))
            .status()
            .expect("generate search records");
        assert!(status.success(), "search index generation failed");
        return;
    }
    let mut generated = generate_document_module(&docs_output);
    let metadata = Command::new(docs_python(repository_root))
        .arg(manifest_dir.join("scripts/extract_sphinx_metadata.py"))
        .arg(&docs_output)
        .args(pages())
        .output()
        .expect("extract Sphinx page metadata");
    assert!(
        metadata.status.success(),
        "Sphinx metadata extraction failed: {}",
        String::from_utf8_lossy(&metadata.stderr)
    );
    generated.push_str(&String::from_utf8(metadata.stdout).expect("UTF-8 Sphinx metadata"));
    fs::write(out_dir.join("sphinx_docs.rs"), generated).expect("write generated Sphinx module");
}

fn docs_python(repository_root: &Path) -> PathBuf {
    env::var_os("COSMOLKIT_DOCS_PYTHON")
        .map(PathBuf::from)
        .filter(|path| path.is_file())
        .or_else(|| {
            let path = repository_root.join(if cfg!(windows) {
                ".venv/Scripts/python.exe"
            } else {
                ".venv/bin/python"
            });
            path.is_file().then_some(path)
        })
        .unwrap_or_else(|| PathBuf::from("python3"))
}

fn build_sphinx(repository_root: &Path, source: &Path, output: &Path) {
    let python = docs_python(repository_root);

    let status = Command::new(&python)
        .args(["-m", "sphinx", "-W", "--keep-going", "-E", "-b", "html"])
        .arg(source)
        .arg(output)
        .current_dir(repository_root)
        .status()
        .unwrap_or_else(|error| panic!("failed to run Sphinx with {}: {error}", python.display()));
    assert!(
        status.success(),
        "Sphinx documentation build failed with status {status}"
    );
}

fn generate_document_module(output: &Path) -> String {
    let mut module = String::from("// Generated by docs-web/build.rs; do not edit.\n");
    let mut toc =
        String::from("pub fn sphinx_toc(page: &str) -> &'static str {\n    match page {\n");
    for page in pages() {
        let path = output.join(format!("{page}.html"));
        let html = fs::read_to_string(&path).unwrap_or_else(|error| {
            panic!("read compiled Sphinx page {}: {error}", path.display())
        });
        let body = extract_article(&html, &path);
        module.push_str(&format!(
            "pub const {}: &str = {};\n",
            constant_name(page),
            rust_literal(&rewrite_page_links(body))
        ));
        toc.push_str(&format!(
            "        {page:?} => {},\n",
            rust_literal(&rewrite_page_links(extract_toc(&html, &path)))
        ));
    }
    toc.push_str("        _ => \"\",\n    }\n}\n");
    module.push_str(&toc);

    let mut stylesheet = String::new();
    for relative in ["_static/pygments.css", "_static/cosmolkit.css"] {
        let path = output.join(relative);
        if let Ok(css) = fs::read_to_string(&path) {
            stylesheet.push_str(&css);
            stylesheet.push('\n');
        }
    }
    module.push_str(&format!(
        "pub const STYLESHEET: &str = {};\n",
        rust_literal(&stylesheet)
    ));
    module
}

// Only rewrite known page destinations, never arbitrary .html resources.
fn canonical_page_link(href: &str) -> Option<String> {
    let suffix_start = href.find(['?', '#']).unwrap_or(href.len());
    let (path, suffix) = href.split_at(suffix_start);
    let path = path
        .strip_prefix("./")
        .or_else(|| path.strip_prefix('/'))
        .unwrap_or(path);
    let page = path.strip_suffix(".html").unwrap_or(path);
    if page.is_empty() {
        return None;
    }
    let route = contract()
        .iter()
        .find(|entry| entry.docname == page || entry.route.trim_start_matches('/') == page)?;
    Some(format!("{}{suffix}", route.route))
}

fn rewrite_anchor(tag: &str) -> String {
    let bytes = tag.as_bytes();
    let mut cursor = 1;
    while cursor < bytes.len() && bytes[cursor].is_ascii_alphabetic() {
        cursor += 1;
    }
    let mut href = None;
    let mut download = false;
    while cursor < bytes.len() {
        while cursor < bytes.len() && bytes[cursor].is_ascii_whitespace() {
            cursor += 1;
        }
        if cursor == bytes.len() || matches!(bytes[cursor], b'>' | b'/') {
            break;
        }
        let start = cursor;
        while cursor < bytes.len()
            && !bytes[cursor].is_ascii_whitespace()
            && !matches!(bytes[cursor], b'=' | b'>' | b'/')
        {
            cursor += 1;
        }
        let name = &tag[start..cursor];
        download |= name.eq_ignore_ascii_case("download");
        while cursor < bytes.len() && bytes[cursor].is_ascii_whitespace() {
            cursor += 1;
        }
        if bytes.get(cursor) != Some(&b'=') {
            continue;
        }
        cursor += 1;
        while cursor < bytes.len() && bytes[cursor].is_ascii_whitespace() {
            cursor += 1;
        }
        let quote = bytes
            .get(cursor)
            .copied()
            .filter(|b| matches!(b, b'\'' | b'"'));
        if quote.is_some() {
            cursor += 1;
        }
        let value_start = cursor;
        while cursor < bytes.len()
            && match quote {
                Some(quote) => bytes[cursor] != quote,
                None => !bytes[cursor].is_ascii_whitespace() && bytes[cursor] != b'>',
            }
        {
            cursor += 1;
        }
        if name.eq_ignore_ascii_case("href") {
            href = Some(value_start..cursor);
        }
        if quote.is_some() && cursor < bytes.len() {
            cursor += 1;
        }
    }
    let mut result = tag.to_owned();
    if !download {
        if let Some(range) = href {
            if let Some(canonical) = canonical_page_link(&tag[range.clone()]) {
                result.replace_range(range, &canonical);
            }
        }
    }
    result
}

fn rewrite_page_links(html: &str) -> String {
    let mut result = String::with_capacity(html.len());
    let lowercase = html.to_ascii_lowercase();
    let mut cursor = 0;
    while let Some(offset) = html[cursor..].find('<') {
        let start = cursor + offset;
        result.push_str(&html[cursor..start]);
        if html[start..].starts_with("<!--") {
            let end = html[start + 4..]
                .find("-->")
                .map_or(html.len(), |n| start + 4 + n + 3);
            result.push_str(&html[start..end]);
            cursor = end;
            continue;
        }
        let mut quote = None;
        let end = html.as_bytes()[start + 1..].iter().position(|&byte| {
            if let Some(current) = quote {
                if byte == current {
                    quote = None;
                }
            } else if matches!(byte, b'\'' | b'"') {
                quote = Some(byte);
            } else if byte == b'>' {
                return true;
            }
            false
        });
        let Some(end) = end.map(|n| start + 1 + n + 1) else {
            cursor = start;
            break;
        };
        let tag = &html[start..end];
        let name = lowercase[start + 1..end - 1]
            .split(|c: char| c.is_ascii_whitespace() || c == '/')
            .next()
            .unwrap_or("");
        if matches!(name, "pre" | "code" | "script" | "style" | "textarea") {
            // These regions may contain literal markup rather than navigation.
            let close = lowercase[end..]
                .find(&format!("</{name}"))
                .map_or(html.len(), |n| end + n);
            result.push_str(&html[start..close]);
            cursor = close;
        } else {
            if name == "a" {
                result.push_str(&rewrite_anchor(tag));
            } else {
                result.push_str(tag);
            }
            cursor = end;
        }
    }
    result.push_str(&html[cursor..]);
    result
}

fn constant_name(page: &str) -> String {
    page.replace(['-', '/'], "_").to_ascii_uppercase()
}

fn extract_article<'a>(html: &'a str, path: &Path) -> &'a str {
    let start_marker = "<article role=\"main\" id=\"furo-main-content\">";
    let start = html
        .find(start_marker)
        .unwrap_or_else(|| panic!("compiled Sphinx page {} has no article", path.display()));
    let content_start = start + start_marker.len();
    let end = html[content_start..]
        .find("</article>")
        .map(|offset| content_start + offset)
        .unwrap_or_else(|| panic!("compiled Sphinx page {} has no article end", path.display()));
    &html[content_start..end]
}

fn extract_toc<'a>(html: &'a str, path: &Path) -> &'a str {
    // Furo's local TOC is outside the article. Keep its markup and IDs verbatim,
    // including nested API entries, rather than deriving a second heading tree.
    let marker = "<div class=\"toc-tree\"";
    let Some(start) = html.find(marker) else {
        return "";
    };
    let malformed = || -> ! {
        panic!(
            "compiled Sphinx page {} has malformed toc-tree container",
            path.display()
        )
    };
    let opening = &html[start + marker.len()..];
    if !opening.starts_with('>') {
        malformed();
    }
    let content_start = start + marker.len() + 1;
    let mut cursor = content_start;
    let mut depth = 1usize;
    while let Some(offset) = html[cursor..].find('<') {
        let tag_start = cursor + offset;
        if html[tag_start..].starts_with("<!--") {
            let Some(end) = html[tag_start + 4..].find("-->") else {
                malformed();
            };
            cursor = tag_start + 4 + end + 3;
            continue;
        }
        let Some(end) = html[tag_start..].find('>') else {
            malformed();
        };
        let tag = &html[tag_start + 1..tag_start + end];
        if tag.contains('<') {
            malformed();
        }
        let name = tag.split_ascii_whitespace().next().unwrap_or("");
        match name {
            "div" => depth += 1,
            "/div" => {
                depth -= 1;
                if depth == 0 {
                    return &html[content_start..tag_start];
                }
            }
            "/aside" | "/body" | "/html" => malformed(),
            _ => {}
        }
        cursor = tag_start + end + 1;
    }
    malformed()
}

fn rust_literal(value: &str) -> String {
    let mut hashes = 1usize;
    while value.contains(&format!("\"{}", "#".repeat(hashes))) {
        hashes += 1;
    }
    let marker = "#".repeat(hashes);
    format!("r{marker}\"{value}\"{marker}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rewrites_known_pages_and_preserves_suffixes() {
        for page in pages().into_iter().chain(["index", "javascript"]) {
            let route = if page == "index" {
                "python".to_string()
            } else if page == "javascript" {
                "javascript".to_string()
            } else {
                format!("python/{page}")
            };
            for prefix in ["", "./", "/"] {
                for extension in [".html", ""] {
                    let input = format!(
                        "<a href=\"{prefix}{page}{extension}?q=a&amp;b=2#section\">{page}.html</a>"
                    );
                    let expected =
                        format!("<a href=\"/{route}?q=a&amp;b=2#section\">{page}.html</a>");
                    assert_eq!(rewrite_page_links(&input), expected);
                }
            }
        }
    }

    #[test]
    fn preserves_non_page_links_and_downloads() {
        for href in [
            "https://example.com/api.html",
            "http://example.com/api.html",
            "//example.com/api.html",
            "mailto:api.html",
            "javascript:api.html",
            "#api.html",
            "?next=api.html",
            "_static/api.html",
            "_downloads/api.html",
            "api.png",
            "unknown.html",
            "../api.html",
            "/other/api.html",
        ] {
            let html = format!("<a href=\"{href}\">link</a>");
            assert_eq!(rewrite_page_links(&html), html);
        }
        for html in [
            "<a href='api.html' download>API</a>",
            "<a DOWNLOAD='api.html' href='api.html'>API</a>",
            "<img src='api.html'><link href='api.html'>",
            "<a data-href='api.html' title=\"href='api.html'\">API</a>",
        ] {
            assert_eq!(rewrite_page_links(html), html);
        }
    }

    #[test]
    fn preserves_code_comments_and_raw_text() {
        for html in [
            "<pre><a href='api.html'>api.html</a></pre>",
            "<code>&lt;a href=\"api.html\"&gt;</code>",
            "<script>const link = '<a href=\"api.html\">';</script>",
            "<style>/* <a href='api.html'> */</style>",
            "<textarea><a href='api.html'>API</a></textarea>",
            "<!-- <a href='api.html'> -->",
            "Text href=\"api.html\" and api.html: 分子",
        ] {
            let input = format!("{html}<a href='api.html'>API</a>");
            assert_eq!(
                rewrite_page_links(&input),
                format!("{html}<a href='/python/api'>API</a>")
            );
        }
    }

    #[test]
    fn handles_anchor_attribute_syntax_without_touching_other_attributes() {
        let html = "<A title='x > y' HREF = 'api.html#atom' data-href='api.html'>分子</A><a href=quickstart.html>Start</a>";
        assert_eq!(
            rewrite_page_links(html),
            "<A title='x > y' HREF = '/python/api#atom' data-href='api.html'>分子</A><a href=/python/quickstart>Start</a>"
        );
    }

    #[test]
    fn preserves_nested_api_toc_verbatim() {
        let toc = "\n<ul><li><a href=\"#cosmolkit.Molecule\">Molecule</a><ul><li><a href=\"#cosmolkit.Molecule.add_hydrogens_\"><code>add_hydrogens_()</code></a></li></ul></li></ul>\n";
        let html =
            format!("<div class=\"toc-tree-container\"><div class=\"toc-tree\">{toc}</div></div>");
        assert_eq!(extract_toc(&html, Path::new("api.html")), toc);
    }

    #[test]
    fn balances_nested_divs_and_skips_comments() {
        let toc =
            "<!-- </div> --><div class=\"nested\"><a href=\"#section\">Section</a></div>after";
        let html = format!("<div class=\"toc-tree\">{toc}</div><div>outside</div>");
        assert_eq!(extract_toc(&html, Path::new("guide.html")), toc);
    }

    #[test]
    fn absent_or_empty_toc_is_empty() {
        for html in [
            "<article>Search</article>",
            "<div class=\"toc-tree-container\"></div>",
            "<div class=\"toc-tree\"></div>",
        ] {
            assert_eq!(extract_toc(html, Path::new("search.html")), "");
        }
    }

    #[test]
    #[should_panic(expected = "compiled Sphinx page broken.html has malformed toc-tree container")]
    fn rejects_unterminated_opening() {
        extract_toc("<div class=\"toc-tree\"", Path::new("broken.html"));
    }

    #[test]
    #[should_panic(expected = "compiled Sphinx page broken.html has malformed toc-tree container")]
    fn rejects_missing_close() {
        extract_toc(
            "<div class=\"toc-tree\"><ul></ul>",
            Path::new("broken.html"),
        );
    }

    #[test]
    #[should_panic(expected = "compiled Sphinx page broken.html has malformed toc-tree container")]
    fn rejects_unclosed_nested_div() {
        extract_toc(
            "<aside><div class=\"toc-tree\"><div>nested</div></aside></body>",
            Path::new("broken.html"),
        );
    }

    #[test]
    fn generates_page_mapping_with_empty_fallback() {
        let directory = env::temp_dir().join(format!(
            "cosmolkit-toc-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&directory).unwrap();
        let toc = "<ul><li><a href=\"api.html#cosmolkit.Atom\">Atom</a></li></ul>";
        let article = "<h1>Document</h1><a href=\"index.html?from=api#top\">Home</a>";
        for page in pages() {
            let mut html =
                format!("<article role=\"main\" id=\"furo-main-content\">{article}</article>");
            if page == "api" {
                html.push_str(&format!("<div class=\"toc-tree\">{toc}</div>"));
            }
            fs::write(directory.join(format!("{page}.html")), html).unwrap();
        }
        let module = generate_document_module(&directory);
        fs::remove_dir_all(&directory).unwrap();
        assert!(module.contains("pub fn sphinx_toc(page: &str) -> &'static str {"));
        for page in pages() {
            let expected = if page == "api" {
                "<ul><li><a href=\"/python/api#cosmolkit.Atom\">Atom</a></li></ul>"
            } else {
                ""
            };
            assert!(module.contains(&format!("{page:?} => {},", rust_literal(expected))));
        }
        assert!(module.contains("_ => \"\","));
        assert!(module.contains(&format!(
            "pub const API: &str = {};",
            rust_literal("<h1>Document</h1><a href=\"/python?from=api#top\">Home</a>")
        )));
    }

    #[test]
    fn raw_literal_preserves_fragment_quotes() {
        assert_eq!(rust_literal("href=\"#api\""), "r##\"href=\"#api\"\"##");
    }
}
