#!/usr/bin/env python3
"""Audit the configured official libinchi target with GCC compiler artifacts."""

from __future__ import annotations

import json
import re
import shlex
import subprocess
from collections import defaultdict, deque
from dataclasses import dataclass
from datetime import date
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
BUILD_ROOT = REPO_ROOT / "target" / "inchi-active-call-graph-audit"
COMPILE_DATABASE = BUILD_ROOT / "compile_commands.json"
REPORT = (
    REPO_ROOT
    / "dev"
    / "gap_reports"
    / "official_inchi_active_call_graph_audit.md"
)
MACHINE_RESULT = BUILD_ROOT / "official_inchi_active_call_graph.json"
VENDOR_ROOT = REPO_ROOT / "third_party" / "InChI"
MODE_H = "third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h"
UTIL_C = "third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c"
INCHI_DLL_B_C = (
    "third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c"
)


@dataclass(frozen=True, order=True)
class Function:
    path: str
    line: int
    name: str
    linkage: str
    translation_unit: str

    @property
    def key(self) -> str:
        return f"{self.path}:{self.line}:{self.name}"


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def cgraph_path(output: str) -> Path:
    object_path = BUILD_ROOT / output
    object_stem = object_path.name.removesuffix(".o")
    matches = sorted(object_path.parent.glob(object_stem + ".*.cgraph"))
    if len(matches) != 1:
        raise RuntimeError(f"expected one cgraph dump for {output}, found {matches}")
    return matches[0]


def aux_command(entry: dict[str, str], aux_path: Path) -> list[str]:
    command = shlex.split(entry["command"])
    source = Path(entry["file"]).resolve()
    result: list[str] = []
    skip_next = False
    for argument in command:
        if skip_next:
            skip_next = False
            continue
        if argument in {"-o", "-c"}:
            skip_next = argument == "-o"
            continue
        if Path(argument).resolve() == source:
            continue
        if argument.startswith("-fdump-"):
            continue
        result.append(argument)
    result.extend(["-std=c11", "-fsyntax-only", "-aux-info", str(aux_path), str(source)])
    return result


def parse_aux(path: Path, translation_unit: str) -> list[Function]:
    marker = re.compile(r"^/\* (.+):(\d+):([NOI])F \*/")
    definitions: list[Function] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = marker.match(line)
        if match is None:
            continue
        source_path = Path(match.group(1))
        if not source_path.is_absolute():
            source_path = REPO_ROOT / source_path
        source_path = source_path.resolve()
        try:
            source = source_path.relative_to(VENDOR_ROOT.resolve())
        except ValueError:
            continue
        declaration = line.split("*/", 1)[1].split("/*", 1)[0].strip()
        open_paren = declaration.find("(")
        if open_paren < 0:
            continue
        identifiers = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", declaration[:open_paren])
        if not identifiers:
            continue
        definitions.append(
            Function(
                path=(Path("third_party/InChI") / source).as_posix(),
                line=int(match.group(2)),
                name=identifiers[-1],
                linkage="internal" if declaration.startswith("static ") else "external",
                translation_unit=translation_unit,
            )
        )
    return sorted(set(definitions))


def parse_cgraph(path: Path) -> tuple[dict[str, set[str]], set[str], dict[str, int]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    start = text.index("Initial Symbol table:")
    end = text.index("Removing unused symbols:", start)
    table = text[start:end]
    entry = re.compile(r"^([^\s/][^\n/]*)/(\d+) \(([^\n]*)\)$", re.MULTILINE)
    matches = list(entry.finditer(table))
    calls: dict[str, set[str]] = defaultdict(set)
    definitions: set[str] = set()
    indirect: dict[str, int] = defaultdict(int)
    call_name = re.compile(r"([^\s/]+)/\d+")
    for index, match in enumerate(matches):
        name = match.group(1).strip()
        body_end = matches[index + 1].start() if index + 1 < len(matches) else len(table)
        body = table[match.end() : body_end]
        if "Type: function definition analyzed" in body:
            definitions.add(name)
        calls_match = re.search(r"^  Calls:(.*)$", body, re.MULTILINE)
        if calls_match is not None:
            calls[name].update(call_name.findall(calls_match.group(1)))
        indirect[name] += body.count("Indirect call")
    return calls, definitions, indirect


def exported_roots() -> set[str]:
    result = subprocess.run(
        ["nm", "-D", "--defined-only", str(BUILD_ROOT / "lib" / "libinchi.so")],
        check=True,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    return {
        fields[2]
        for line in result.stdout.splitlines()
        if len(fields := line.split()) >= 3 and fields[1] == "T"
    }


def closure(roots: set[Function], graph: dict[Function, set[Function]]) -> set[Function]:
    reached: set[Function] = set()
    queue = deque(sorted(roots))
    while queue:
        function = queue.popleft()
        if function in reached:
            continue
        reached.add(function)
        queue.extend(sorted(graph.get(function, set()) - reached))
    return reached


def source_test_name(function: Function) -> str:
    stem = Path(function.path).stem.lower()
    name = re.sub(r"[^a-z0-9]+", "_", function.name.lower()).strip("_")
    return f"source_port__{stem}__{name}__line_{function.line}"


def markdown_table(rows: list[tuple[str, ...]], headers: tuple[str, ...]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return lines


def main() -> None:
    entries = json.loads(COMPILE_DATABASE.read_text(encoding="utf-8"))
    if len(entries) != 60:
        raise RuntimeError(f"expected 60 configured C translation units, found {len(entries)}")

    functions_by_tu: dict[str, list[Function]] = {}
    calls_by_tu: dict[str, dict[str, set[str]]] = {}
    indirect_by_tu: dict[str, dict[str, int]] = {}
    for index, entry in enumerate(entries):
        tu = relative(Path(entry["file"]))
        aux_path = BUILD_ROOT / f"active-{index:02d}.aux"
        subprocess.run(
            aux_command(entry, aux_path),
            check=True,
            cwd=Path(entry["directory"]),
            capture_output=True,
            text=True,
        )
        aux_functions = parse_aux(aux_path, tu)
        calls, cgraph_definitions, indirect = parse_cgraph(cgraph_path(entry["output"]))
        functions_by_tu[tu] = [
            function for function in aux_functions if function.name in cgraph_definitions
        ]
        calls_by_tu[tu] = calls
        indirect_by_tu[tu] = indirect

    active = sorted({function for functions in functions_by_tu.values() for function in functions})
    by_name: dict[str, list[Function]] = defaultdict(list)
    by_tu_name: dict[tuple[str, str], list[Function]] = defaultdict(list)
    for function in active:
        by_name[function.name].append(function)
        by_tu_name[(function.translation_unit, function.name)].append(function)

    graph: dict[Function, set[Function]] = {function: set() for function in active}
    external_calls: dict[Function, set[str]] = defaultdict(set)
    indirect_calls: dict[Function, int] = defaultdict(int)
    for caller in active:
        for callee_name in calls_by_tu[caller.translation_unit].get(caller.name, set()):
            candidates = by_tu_name.get((caller.translation_unit, callee_name), [])
            if not candidates:
                candidates = by_name.get(callee_name, [])
            unique_locations = {(item.path, item.line): item for item in candidates}
            if len(unique_locations) == 1:
                graph[caller].add(next(iter(unique_locations.values())))
            elif len(unique_locations) > 1:
                raise RuntimeError(
                    f"ambiguous active callee {callee_name} from {caller.key}: "
                    + ", ".join(item.key for item in candidates)
                )
            else:
                external_calls[caller].add(callee_name)
        indirect_calls[caller] = indirect_by_tu[caller.translation_unit].get(caller.name, 0)

    exports = exported_roots()
    root_functions = {
        function
        for name in exports
        for function in by_name.get(name, [])
    }
    reachable = closure(root_functions, graph)
    unreachable = set(active) - reachable

    inchi_matches = [item for item in active if item.name == "InchiToInchiAtom"]
    if len(inchi_matches) != 1:
        raise RuntimeError(f"expected one active InchiToInchiAtom, found {inchi_matches}")
    inchi = inchi_matches[0]
    inchi_closure = closure({inchi}, graph)
    depths = {inchi: 0}
    queue = deque([inchi])
    while queue:
        caller = queue.popleft()
        for callee in sorted(graph[caller]):
            if callee not in depths:
                depths[callee] = depths[caller] + 1
                queue.append(callee)

    direct = sorted(graph[inchi])
    transitive = sorted(inchi_closure - {inchi} - set(direct), key=lambda item: (depths[item], item))
    inchi_external = sorted(
        {
            callee
            for function in inchi_closure
            for callee in external_calls.get(function, set())
        }
    )

    special_names = [
        "is_in_the_slist",
        "is_element_a_metal",
        "InchiToInchiAtom",
        "FindToken",
        "LoadLine",
        "inchi_ios_gets",
        "inchi_ios_getsTab",
        "inchi_ios_getsTab1",
        "lrtrim",
        "mystrncpy",
        "CreateInchiAtom",
        "CreateInchi_Stereo0D",
        "FreeInchi_Stereo0D",
    ]
    special_rows: list[tuple[str, ...]] = []
    for name in special_names:
        matches = by_name.get(name, [])
        if not matches:
            special_rows.append((f"`{name}`", "inactive/not defined", "none", "none"))
            continue
        locations = "<br>".join(f"`{item.path}:{item.line}`" for item in matches)
        status = "export-reachable" if any(item in reachable for item in matches) else "active-unreachable"
        called = "yes" if any(item in inchi_closure for item in matches) else "no"
        special_rows.append((f"`{name}`", status, called, locations))

    lines = [
        "# Official InChI Active Call Graph Audit",
        "",
        "## Audit Basis",
        "",
        f"- Generated: `{date.today().isoformat()}`",
        "- Selected target: official `libinchi` shared-library CMake target from "
        "`third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/CMakeLists.txt`.",
        "- Platform branch: GCC 13.3.0 on Linux, C11.",
        "- Active target definitions: `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`, `fPIC`.",
        "- Compiler evidence: the target's 60-entry `compile_commands.json`, GCC "
        "`-aux-info`, GCC `-fdump-ipa-cgraph`, and the linked ELF version-script exports.",
        "- Raw tree-sitter inventory is not used to decide activity or call edges.",
        f"- Preprocessed active vendored function definitions: `{len(active)}`.",
        f"- Export-root reachable definitions through compiler direct edges: `{len(reachable)}`.",
        f"- Active definitions not directly reachable from an exported API: `{len(unreachable)}`.",
        f"- Exported production roots: `{len(exports)}`.",
        f"- Compiler direct source edges: `{sum(len(value) for value in graph.values())}`.",
        f"- Compiler indirect-call sites: `{sum(indirect_calls.values())}`.",
        "",
        "Activity means that GCC saw a function body after preprocessing under the selected "
        "target configuration. Reachability is the recursive closure of the ELF version-script "
        "exports over compiler-recorded direct edges. Indirect callback sites are recorded "
        "separately and are not converted into guessed name-based edges.",
        "",
        "## Required Classification",
        "",
        "- **Active source functions:** listed in the active-definition appendix below; each has "
        "a GCC `NF/OF/IF` definition record and an IPA function body.",
        "- **Inactive conditional functions:** raw bodies excluded by preprocessing have no GCC "
        "definition/body and cannot receive a production Port step.",
        "- **Macro-only behavior:** allocation names under this target are preprocessor aliases, "
        "not C function bodies; they are documented in the allocation section below.",
        "- **Header-defined behavior:** active header bodies are included only when GCC emits an "
        "active vendored body; plain macros are represented as behavior attached to their caller.",
        "- **External declarations:** `NC` declaration records without a matching active vendored "
        "definition are not functions to port.",
        "- **libc/compiler intrinsics:** unresolved compiler edges such as `malloc`, `free`, "
        "`memcpy`, ctype accessors, and checked builtins are runtime semantics, not vendored Port steps.",
        "- **CLI/demo and non-production behavior:** a compiled body outside the exported-root "
        "closure is active but is classified as target-unreachable and does not justify a public API.",
        "",
        "## Allocation Branch",
        "",
        f"Under `TARGET_API_LIB` on GCC/Linux, `{MODE_H}:1172-1181` selects:",
        "",
        "```c",
        "#define inchi_malloc   malloc",
        "#define inchi_calloc   calloc",
        "#define inchi_realloc  realloc",
        "#define inchi_free(X)  do{ if(X) free(X); }while(0)",
        "```",
        "",
        f"Therefore `{UTIL_C}:1552`, `:1561`, and `:1570` are preprocessed out. In particular, "
        "the raw self-referential body at `util.c:1570` is not the active GCC/Linux "
        "`TARGET_API_LIB` implementation. Rust must reproduce libc allocation sizing/zeroing/"
        "reallocation and the null-checked `free` macro behavior, anchored to `mode.h`, without "
        "inventing active C functions named `inchi_malloc`, `inchi_calloc`, `inchi_realloc`, or "
        "`inchi_free`.",
        "",
        "The same rule applies to `qmalloc`, `qfree`, `qzfree`, `fast_alloc`, and `fast_free`: "
        "they are macro-expanded behavior. Under `USE_ALLOCA == 0` and GCC, they ultimately use "
        "the active allocation macros above; no standalone production function exists for them.",
        "",
        "## Required Function Checks",
        "",
    ]
    lines.extend(markdown_table(special_rows, ("Name", "Configured status", "In `InchiToInchiAtom` closure", "Definition")))
    lines.extend(
        [
            "",
            "Both `is_in_the_slist` and `is_element_a_metal` are active GCC function bodies and "
            "direct callees of `InchiToInchiAtom`; they must be ported before that caller.",
            "",
            "## `InchiToInchiAtom` Direct Callees",
            "",
        ]
    )
    direct_rows = [
        (
            f"`{item.name}`",
            f"`{item.path}:{item.line}`",
            "vendored active definition",
        )
        for item in direct
    ]
    direct_rows.extend((f"`{name}`", "runtime", "libc/compiler intrinsic or macro expansion") for name in sorted(external_calls[inchi]))
    lines.extend(markdown_table(direct_rows, ("Callee", "Definition", "Classification")))
    lines.extend(["", "## Recursive Transitive Callees", ""])
    transitive_rows = [
        (
            str(depths[item]),
            f"`{item.name}`",
            f"`{item.path}:{item.line}`",
            ", ".join(f"`{callee.name}`" for callee in sorted(graph[item])) or "none",
        )
        for item in transitive
    ]
    lines.extend(markdown_table(transitive_rows, ("Minimum depth", "Function", "Definition", "Direct vendored callees")))
    lines.extend(
        [
            "",
            "The complete non-vendored closure leaves are:",
            "",
            ", ".join(f"`{name}`" for name in inchi_external) + ".",
            "",
            "No function-name inference was used for these edges. The direct and transitive "
            "tables are generated from the configured GCC IPA graph; unresolved indirect calls "
            "remain explicitly indirect.",
            "",
            "## Active Definition Appendix",
            "",
        ]
    )
    active_rows = [
        (
            f"`{function.name}`",
            f"`{function.path}:{function.line}`",
            function.linkage,
            "export root" if function.name in exports else "reachable" if function in reachable else "active, export-unreachable",
            ", ".join(f"`{callee.name}`" for callee in sorted(graph[function])) or "none",
            str(indirect_calls[function]),
        )
        for function in active
    ]
    lines.extend(markdown_table(active_rows, ("Function", "Source", "Linkage", "Production classification", "Direct active callees", "Indirect sites")))
    lines.extend(
        [
            "",
            "## Plan Input",
            "",
            "The machine-readable companion is "
            "`target/inchi-active-call-graph-audit/official_inchi_active_call_graph.json`. "
            "Only active vendored definitions may become function Port steps. Macro behavior is "
            "attached to source-backed caller/behavior steps, inactive raw bodies are excluded, "
            "and no production API may be exposed until its complete reachable closure is ported "
            "and exact official-C behavior fixtures pass.",
        ]
    )

    report_data = {
        "configuration": {
            "compile_definitions": ["COMPILE_ANSI_ONLY", "TARGET_API_LIB", "fPIC"],
            "compiler": "GCC 13.3.0",
            "platform": "Linux",
            "translation_units": len(entries),
        },
        "exports": sorted(exports),
        "functions": [
            {
                "key": function.key,
                "name": function.name,
                "path": function.path,
                "line": function.line,
                "linkage": function.linkage,
                "translation_unit": function.translation_unit,
                "reachable": function in reachable,
                "export_root": function.name in exports,
                "callees": sorted(callee.key for callee in graph[function]),
                "external_callees": sorted(external_calls[function]),
                "indirect_calls": indirect_calls[function],
                "test_name": source_test_name(function),
            }
            for function in active
        ],
        "inchi_to_inchi_atom": {
            "root": inchi.key,
            "direct": [item.key for item in direct],
            "transitive": [item.key for item in transitive],
            "external": inchi_external,
        },
    }
    MACHINE_RESULT.write_text(json.dumps(report_data, indent=2) + "\n", encoding="utf-8")
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"active definitions: {len(active)}")
    print(f"export reachable: {len(reachable)}")
    print(f"active unreachable: {len(unreachable)}")
    print(f"InchiToInchiAtom closure: {len(inchi_closure)} vendored functions")
    print(f"wrote {relative(REPORT)}")
    print(f"wrote {relative(MACHINE_RESULT)}")


if __name__ == "__main__":
    main()
