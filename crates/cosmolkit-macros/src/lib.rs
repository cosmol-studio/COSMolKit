use proc_macro::{Literal, TokenStream, TokenTree};
use std::path::PathBuf;

#[proc_macro]
pub fn rdkit_uff_params(input: TokenStream) -> TokenStream {
    let path = manifest_relative_path(input);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()));
    let mut out = String::from("&[");
    for (line_no, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(
            cols.len() == 3,
            "{}:{} expected 3 tab-separated columns",
            path.display(),
            line_no + 1
        );
        let key = escape_str(cols[0]);
        let r1: f64 = cols[1]
            .parse()
            .unwrap_or_else(|err| panic!("{}:{} invalid r1: {err}", path.display(), line_no + 1));
        let xi: f64 = cols[2]
            .parse()
            .unwrap_or_else(|err| panic!("{}:{} invalid Xi: {err}", path.display(), line_no + 1));
        out.push_str(&format!("(\"{key}\",{r1:?}_f64,{xi:?}_f64),"));
    }
    out.push(']');
    out.parse().expect("generated UFF params tokens")
}

#[proc_macro]
pub fn rdkit_periodic_rvdw(input: TokenStream) -> TokenStream {
    let path = manifest_relative_path(input);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()));
    let mut out = String::from("&[");
    for (line_no, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(
            cols.len() == 2,
            "{}:{} expected 2 tab-separated columns",
            path.display(),
            line_no + 1
        );
        let atomic_num: u8 = cols[0].parse().unwrap_or_else(|err| {
            panic!(
                "{}:{} invalid atomic number: {err}",
                path.display(),
                line_no + 1
            )
        });
        let rvdw: f64 = cols[1]
            .parse()
            .unwrap_or_else(|err| panic!("{}:{} invalid rvdw: {err}", path.display(), line_no + 1));
        out.push_str(&format!("({atomic_num}_u8,{rvdw:?}_f64),"));
    }
    out.push(']');
    out.parse().expect("generated periodic rvdw tokens")
}

fn manifest_relative_path(input: TokenStream) -> PathBuf {
    let Some(TokenTree::Literal(lit)) = input.into_iter().next() else {
        panic!("expected a single string literal path");
    };
    let relative = parse_string_literal(lit);
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR must be set");
    PathBuf::from(manifest_dir).join(relative)
}

fn parse_string_literal(lit: Literal) -> String {
    let source = lit.to_string();
    source
        .strip_prefix('"')
        .and_then(|value| value.strip_suffix('"'))
        .expect("expected a plain string literal path")
        .to_string()
}

fn escape_str(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace('\t', "\\t")
}
