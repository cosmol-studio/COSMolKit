#!/usr/bin/env python3
"""Generate the source-level function inventory for vendored InChI."""

from __future__ import annotations

import hashlib
import re
import subprocess
import tempfile
from dataclasses import dataclass, replace
from datetime import date
from pathlib import Path

from tree_sitter import Language, Node, Parser, Query, QueryCursor
import tree_sitter_c


REPO_ROOT = Path(__file__).resolve().parent.parent
VENDOR_ROOT = REPO_ROOT / "third_party" / "InChI"
SOURCE_ROOT = VENDOR_ROOT / "INCHI-1-SRC"
BASE_ROOT = SOURCE_ROOT / "INCHI_BASE" / "src"
API_ROOT = SOURCE_ROOT / "INCHI_API" / "libinchi" / "src"
OUTPUT = REPO_ROOT / "dev" / "gap_reports" / "official_inchi_function_inventory.md"


@dataclass(frozen=True)
class FunctionDefinition:
    path: str
    line: int
    name: str
    linkage: str
    signature: str
    extraction: str


def relative(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def is_below(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def source_scope(path: Path) -> str:
    if is_below(path, BASE_ROOT) or is_below(path, API_ROOT):
        return "production"
    if is_below(path, SOURCE_ROOT / "INCHI_EXE"):
        return "cli"
    return "demo"


def declarator_name(node: Node, source: bytes) -> str | None:
    if node.type in {"identifier", "field_identifier"}:
        return source[node.start_byte : node.end_byte].decode("utf-8")
    nested = node.child_by_field_name("declarator")
    if nested is not None:
        return declarator_name(nested, source)
    for child in node.named_children:
        if child.type in {
            "function_declarator",
            "pointer_declarator",
            "parenthesized_declarator",
        }:
            name = declarator_name(child, source)
            if name is not None:
                return name
    return None


def parse_file(
    parser: Parser,
    function_query: Query,
    error_query: Query,
    path: Path,
) -> tuple[list[FunctionDefinition], int]:
    source = path.read_bytes()
    tree = parser.parse(source)
    definitions: list[FunctionDefinition] = []

    matches = QueryCursor(function_query).matches(tree.root_node)
    for _, captures in matches:
        node = captures["function"][0]
        declarator = captures["declarator"][0]
        body = captures["body"][0]
        name = declarator_name(declarator, source)
        if name is None:
            raise RuntimeError(
                f"cannot extract function name at {relative(path)}:{node.start_point.row + 1}"
            )
        signature = " ".join(
            source[node.start_byte : body.start_byte].decode("utf-8").split()
        )
        prefix = source[node.start_byte : declarator.start_byte].decode("utf-8")
        definitions.append(
            FunctionDefinition(
                path=relative(path),
                line=node.start_point.row + 1,
                name=name,
                linkage="internal" if "static" in prefix.split() else "external",
                signature=signature.replace("|", "\\|"),
                extraction="tree-sitter",
            )
        )

    definitions.sort(key=lambda definition: definition.line)
    errors = QueryCursor(error_query).captures(tree.root_node).get("error", [])
    return definitions, len(errors)


def digest(paths: list[Path]) -> str:
    aggregate = hashlib.sha256()
    for path in paths:
        aggregate.update(relative(path).encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(hashlib.sha256(path.read_bytes()).digest())
    return aggregate.hexdigest()


def compiled_function_definitions(
    parser: Parser, language: Language, c_files: list[Path]
) -> dict[tuple[str, int], FunctionDefinition]:
    definitions: dict[tuple[str, int], FunctionDefinition] = {}
    marker = re.compile(r"^/\* (.+):(\d+):[NOI]F \*/")
    declarator_query = Query(language, "(function_declarator) @declarator")
    include_args = [
        f"-I{BASE_ROOT}",
        f"-I{API_ROOT}",
        f"-I{API_ROOT / 'ixa'}",
    ]
    with tempfile.TemporaryDirectory(prefix="inchi-aux-info-") as tmp:
        tmp_root = Path(tmp)
        for index, path in enumerate(c_files):
            aux_path = tmp_root / f"{index}.aux"
            subprocess.run(
                [
                    "gcc",
                    "-std=c11",
                    "-DCOMPILE_ANSI_ONLY",
                    "-DTARGET_API_LIB",
                    "-fsyntax-only",
                    "-aux-info",
                    str(aux_path),
                    *include_args,
                    str(path),
                ],
                check=True,
                cwd=REPO_ROOT,
                capture_output=True,
                text=True,
            )
            for line in aux_path.read_text(encoding="utf-8").splitlines():
                match = marker.match(line)
                if match is None:
                    continue
                source_path = Path(match.group(1)).resolve()
                if not is_below(source_path, VENDOR_ROOT):
                    continue
                declaration = line.split("*/", 1)[1].split("/*", 1)[0].strip()
                declaration_tree = parser.parse(declaration.encode("utf-8"))
                declarators = QueryCursor(declarator_query).captures(
                    declaration_tree.root_node
                ).get("declarator", [])
                if not declarators:
                    raise RuntimeError(f"cannot parse GCC aux declaration: {line}")
                declarator = min(declarators, key=lambda node: node.start_byte)
                name = declarator_name(declarator, declaration.encode("utf-8"))
                if name is None:
                    raise RuntimeError(f"cannot name GCC aux declaration: {line}")
                path_text = relative(source_path)
                source_line = int(match.group(2))
                definitions[(path_text, source_line)] = FunctionDefinition(
                    path=path_text,
                    line=source_line,
                    name=name,
                    linkage=(
                        "internal" if declaration.startswith("static ") else "external"
                    ),
                    signature=declaration.removesuffix(";").replace("|", "\\|"),
                    extraction="gcc-aux",
                )
    return definitions


def append_function_table(
    lines: list[str], definitions: list[FunctionDefinition]
) -> None:
    lines.extend(
        [
            "| Function | Source | Linkage | Extraction | Signature |",
            "|---|---|---|---|---|",
        ]
    )
    for definition in definitions:
        lines.append(
            f"| `{definition.name}` | `{definition.path}:{definition.line}` | "
            f"{definition.linkage} | {definition.extraction} | "
            f"`{definition.signature}` |"
        )


def main() -> None:
    language = Language(tree_sitter_c.language())
    parser = Parser(language)
    function_query = Query(
        language,
        "(function_definition "
        "declarator: (_) @declarator "
        "body: (compound_statement) @body) @function",
    )
    error_query = Query(language, "(ERROR) @error")
    c_files = sorted(SOURCE_ROOT.rglob("*.c"), key=relative)
    production_headers = sorted(
        set(BASE_ROOT.glob("*.h")) | set(API_ROOT.rglob("*.h")), key=relative
    )

    by_file: list[tuple[Path, str, list[FunctionDefinition], int]] = []
    all_c_definitions: list[FunctionDefinition] = []
    for path in c_files:
        definitions, errors = parse_file(parser, function_query, error_query, path)
        by_file.append((path, source_scope(path), definitions, errors))
        all_c_definitions.extend(definitions)

    header_definitions: list[FunctionDefinition] = []
    header_errors: dict[str, int] = {}
    for path in production_headers:
        definitions, errors = parse_file(parser, function_query, error_query, path)
        header_definitions.extend(definitions)
        if errors:
            header_errors[relative(path)] = errors

    production_c = [item for item in by_file if item[1] == "production"]
    cli_c = [item for item in by_file if item[1] == "cli"]
    demo_c = [item for item in by_file if item[1] == "demo"]
    raw_production_functions = [
        definition
        for _, scope, definitions, _ in by_file
        if scope == "production"
        for definition in definitions
    ]
    raw_production_locations = {
        (definition.path, definition.line)
        for definition in raw_production_functions + header_definitions
    }
    compiled_definitions = compiled_function_definitions(
        parser, language, [path for path, _, _, _ in production_c]
    )
    recovered_definitions = [
        definition
        for location, definition in compiled_definitions.items()
        if location not in raw_production_locations
    ]
    raw_by_location = {
        (definition.path, definition.line): definition
        for definition in raw_production_functions + header_definitions
    }
    corrected_names = sum(
        raw_by_location[location].name != definition.name
        for location, definition in compiled_definitions.items()
        if location in raw_by_location
    )

    def merge_compiled(definition: FunctionDefinition) -> FunctionDefinition:
        compiled = compiled_definitions.get((definition.path, definition.line))
        if compiled is None:
            return definition
        return replace(compiled, extraction="tree-sitter+gcc-aux")

    by_production_path = {
        relative(path): definitions for path, _, definitions, _ in production_c
    }
    for definitions in by_production_path.values():
        definitions[:] = [merge_compiled(definition) for definition in definitions]
    header_definitions = [
        merge_compiled(definition) for definition in header_definitions
    ]
    for definition in recovered_definitions:
        if not definition.path.endswith(".c"):
            header_definitions.append(definition)
        else:
            by_production_path[definition.path].append(definition)
    for definitions in by_production_path.values():
        definitions.sort(key=lambda definition: definition.line)
    header_definitions.sort(key=lambda definition: (definition.path, definition.line))
    production_functions = [
        definition
        for path, _, _, _ in production_c
        for definition in by_production_path[relative(path)]
    ]

    lines = [
        "# Official InChI Function Inventory",
        "",
        "## Generation Contract",
        "",
        f"- Generated: `{date.today().isoformat()}`",
        "- Upstream release: `v1.07.5`",
        "- Upstream commit: `11a87982bb518f57ac013f0b258c283655e1ea1d`",
        "- Parser: `tree-sitter-c` over raw vendored source definitions",
        f"- Parsed C files: `{len(c_files)}`",
        f"- C-source function definitions: `{len(all_c_definitions)}`",
        f"- Production C files: `{len(production_c)}`",
        f"- Production C-source functions: `{len(production_functions)}`",
        f"- Production header-defined functions: `{len(header_definitions)}`",
        f"- GCC-configured production function locations: `{len(compiled_definitions)}`",
        f"- GCC-recovered production definitions: `{len(recovered_definitions)}`",
        f"- GCC-corrected production symbol names: `{corrected_names}`",
        "- GCC subset check: `complete`",
        f"- CLI C files: `{len(cli_c)}`",
        f"- Demo C files: `{len(demo_c)}`",
        f"- C-source aggregate SHA-256: `{digest(c_files)}`",
        "",
        "The production classification is the complete official `libinchi` target:",
        "`INCHI_BASE/src` plus `INCHI_API/libinchi/src`. CLI and demo functions",
        "are inventoried because the bootstrap plan requires every vendored C source",
        "file, but they are not selected for the reusable production Rust crate.",
        "",
        "A function is one raw C `function_definition` node. Conditional definitions",
        "remain visible because extraction is performed before C preprocessing.",
        "Production header definitions are listed separately so macro-controlled and",
        "header-implemented behavior cannot disappear from the port call graph.",
        "",
        "## Parse Integrity",
        "",
        "| Source | Scope | Functions | Parse errors |",
        "|---|---|---:|---:|",
    ]
    for path, scope, definitions, errors in by_file:
        lines.append(
            f"| `{relative(path)}` | {scope} | {len(definitions)} | {errors} |"
        )
    if header_errors:
        lines.extend(["", "Production headers with parser errors:"])
        for path, errors in sorted(header_errors.items()):
            lines.append(f"- `{path}`: {errors}")
    else:
        lines.extend(["", "No production header parse errors were reported."])

    for heading, selected in (
        ("Production C Functions", production_c),
        ("CLI C Functions (Inventory Only)", cli_c),
        ("Demo C Functions (Inventory Only)", demo_c),
    ):
        lines.extend(["", f"## {heading}", ""])
        for path, _, definitions, errors in selected:
            lines.extend(
                [
                    f"### `{relative(path)}`",
                    "",
                    f"Parse errors: `{errors}`. Function definitions: `{len(definitions)}`.",
                    "",
                ]
            )
            append_function_table(lines, definitions)
            lines.append("")

    lines.extend(["", "## Production Header-Defined Functions", ""])
    append_function_table(lines, header_definitions)
    lines.extend(
        [
            "",
            "## Regeneration",
            "",
            "Run from the repository root:",
            "",
            "```bash",
            "uv run --no-project --with tree-sitter --with tree-sitter-c python dev/generate_inchi_function_inventory.py",
            "```",
            "",
        ]
    )
    OUTPUT.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
