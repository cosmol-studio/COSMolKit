#!/usr/bin/env python3
"""Regenerate the remaining InChI plan from the RDKit core API closure."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
BUILD_ROOT = REPO_ROOT / "target" / "inchi-active-call-graph-audit"
AUDIT_JSON = BUILD_ROOT / "official_inchi_active_call_graph.json"
COMPILE_DATABASE = BUILD_ROOT / "compile_commands.json"
PLAN_PATH = REPO_ROOT / "dev" / "rdkit_inchi_full_port_plan.md"
REPORT_PATH = (
    REPO_ROOT / "dev" / "gap_reports" / "rdkit_inchi_core_scope_audit.md"
)
PREFIX_END_STEP = 4720
ROOT_NAMES = (
    "GetINCHI",
    "FreeINCHI",
    "GetStructFromINCHI",
    "FreeStructFromINCHI",
    "GetINCHIKeyFromINCHI",
)
READ_ACTION = (
    "Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` "
    "to reload and follow the required execution standard, source reproduction "
    "rules, artifact requirements, no-git rule, and completion criteria for the "
    "next task."
)
CORE_READ_ACTION = (
    "Read `dev/policy_invariants.md`, `dev/source_reproduction_protocol.md`, "
    "and `dev/README.md` to reload and follow the required execution standard, "
    "source reproduction rules, core operation rules, artifact requirements, "
    "no-git rule, and completion criteria for the next task."
)


@dataclass(frozen=True, order=True)
class Function:
    key: str
    name: str
    path: str
    line: int
    linkage: str
    translation_unit: str
    test_name: str

    @property
    def location(self) -> tuple[str, int]:
        return self.path, self.line


@dataclass(frozen=True)
class AdapterFunction:
    name: str
    path: str
    line: int
    responsibility: str

    @property
    def test_name(self) -> str:
        stem = Path(self.path).stem.lower()
        normalized = normalize(self.name)
        return f"source_port__{stem}__{normalized}__line_{self.line}"


ADAPTER_FUNCTIONS = (
    AdapterFunction(
        "assignBondDirs",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        91,
        "double-bond direction constraint propagation used by `InchiToMol`",
    ),
    AdapterFunction(
        "findAlternatingBonds",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        195,
        "recursive alternating-bond path search used by cleanup helpers",
    ),
    AdapterFunction(
        "getNumDoubleBondedNegativelyChargedNeighboringSi",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        306,
        "silicon-neighbor count used by valence cleanup",
    ),
    AdapterFunction(
        "_Valence4NCleanUp1",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        324,
        "post-parse nitrogen cleanup",
    ),
    AdapterFunction(
        "_Valence4NCleanUp2",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        376,
        "post-parse nitrogen cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUp1",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        393,
        "post-parse nitrogen cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUp2",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        417,
        "post-parse nitrogen cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUp3",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        443,
        "post-parse nitrogen cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUp4",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        477,
        "post-parse nitrogen/silicon cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUp5",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        544,
        "post-parse nitrogen/heteroatom cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUp6",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        625,
        "post-parse nitrogen ring cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUp7",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        679,
        "post-parse nitrogen/oxygen ring cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUp8",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        750,
        "post-parse nitrogen ring cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUp9",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        804,
        "post-parse nitrogen ring cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUpA",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        855,
        "post-parse tin-substitution cleanup",
    ),
    AdapterFunction(
        "_Valence5NCleanUpB",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        916,
        "post-parse carbon-path cleanup",
    ),
    AdapterFunction(
        "_Valence7SCleanUp1",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        941,
        "post-parse sulfur cleanup",
    ),
    AdapterFunction(
        "_Valence7SCleanUp2",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        989,
        "post-parse sulfur/nitrogen cleanup",
    ),
    AdapterFunction(
        "_Valence7SCleanUp3",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        1021,
        "post-parse sulfur/nitrogen cleanup",
    ),
    AdapterFunction(
        "_Valence8SCleanUp1",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        1044,
        "post-parse sulfur cleanup",
    ),
    AdapterFunction(
        "_Valence8ClCleanUp1",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        1079,
        "post-parse chlorine cleanup",
    ),
    AdapterFunction(
        "_Valence5ClCleanUp1",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        1114,
        "post-parse chlorine cleanup",
    ),
    AdapterFunction(
        "_Valence3ClCleanUp1",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        1134,
        "post-parse chlorine/sulfur cleanup",
    ),
    AdapterFunction(
        "cleanUp",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        1151,
        "ordered `InchiToMol` post-parse cleanup dispatcher",
    ),
    AdapterFunction(
        "InchiToMol",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        1254,
        "`Chem.MolFromInchi` graph, isotope, hydrogen, bond, stereo, and diagnostic conversion",
    ),
    AdapterFunction(
        "fixOptionSymbol",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        1674,
        "GCC/Linux option-prefix normalization used by `MolToInchi`",
    ),
    AdapterFunction(
        "rCleanUp",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        1694,
        "pre-generation reverse cleanup used by `MolToInchi`",
    ),
    AdapterFunction(
        "MolToInchi",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        1747,
        "`Chem.MolToInchi` graph-to-`inchi_Input` conversion and diagnostics",
    ),
    AdapterFunction(
        "InchiToInchiKey",
        "third_party/rdkit/External/INCHI-API/inchi.cpp",
        2145,
        "`InchiToInchiKey` return-code mapping",
    ),
    AdapterFunction(
        "MolToInchiKey",
        "third_party/rdkit/External/INCHI-API/inchi.h",
        107,
        "`Chem.MolToInchiKey` composition",
    ),
)


def normalize(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def cgraph_path(output: str) -> Path:
    object_path = BUILD_ROOT / output
    stem = object_path.name.removesuffix(".o")
    matches = sorted(object_path.parent.glob(stem + ".*.cgraph"))
    if len(matches) != 1:
        raise RuntimeError(f"expected one cgraph dump for {output}, found {matches}")
    return matches[0]


def load_graph() -> tuple[
    dict[str, Function],
    dict[Function, set[Function]],
    dict[Function, set[Function]],
]:
    data = json.loads(AUDIT_JSON.read_text(encoding="utf-8"))
    functions = {
        item["key"]: Function(
            key=item["key"],
            name=item["name"],
            path=item["path"],
            line=item["line"],
            linkage=item["linkage"],
            translation_unit=item["translation_unit"],
            test_name=item["test_name"],
        )
        for item in data["functions"]
    }
    graph = {
        functions[item["key"]]: {functions[key] for key in item["callees"]}
        for item in data["functions"]
    }
    by_name: dict[str, list[Function]] = defaultdict(list)
    for function in functions.values():
        by_name[function.name].append(function)
    ambiguous = {
        name: matches
        for name, matches in by_name.items()
        if len({item.location for item in matches}) != 1
    }

    entries = json.loads(COMPILE_DATABASE.read_text(encoding="utf-8"))
    entry_pattern = re.compile(
        r"^([^\s/][^\n/]*)/(\d+) \(([^\n]*)\)$", re.MULTILINE
    )
    address_pattern = re.compile(r"([^\s/]+)/\d+ \(addr\)")
    address_names: dict[tuple[str, str], set[str]] = defaultdict(set)
    for entry in entries:
        translation_unit = (
            Path(entry["file"]).resolve().relative_to(REPO_ROOT).as_posix()
        )
        text = cgraph_path(entry["output"]).read_text(
            encoding="utf-8", errors="replace"
        )
        start = text.index("Initial Symbol table:")
        end = text.index("Removing unused symbols:", start)
        table = text[start:end]
        matches = list(entry_pattern.finditer(table))
        for index, match in enumerate(matches):
            body_end = (
                matches[index + 1].start()
                if index + 1 < len(matches)
                else len(table)
            )
            body = table[match.end() : body_end]
            references = re.search(r"^  References:(.*)$", body, re.MULTILINE)
            if references is not None:
                address_names[(translation_unit, match.group(1).strip())].update(
                    address_pattern.findall(references.group(1))
                )

    callback_edges: dict[Function, set[Function]] = defaultdict(set)
    selected_names = set(ROOT_NAMES)
    roots = {
        function
        for function in functions.values()
        if function.name in selected_names
    }
    if {item.name for item in roots} != selected_names or len(roots) != len(
        selected_names
    ):
        raise RuntimeError(f"selected root resolution failed: {sorted(roots)}")

    while True:
        selected = closure(roots, graph)
        changed = False
        for caller in selected:
            names = address_names.get(
                (caller.translation_unit, caller.name), set()
            )
            for name in names:
                if name in ambiguous:
                    raise RuntimeError(
                        f"ambiguous address-taken callback {name} from {caller.key}"
                    )
                matches = by_name.get(name, [])
                if not matches:
                    continue
                target = matches[0]
                if target not in graph[caller]:
                    graph[caller].add(target)
                    callback_edges[caller].add(target)
                    changed = True
        if not changed:
            break
    return functions, graph, callback_edges


def closure(
    roots: set[Function], graph: dict[Function, set[Function]]
) -> set[Function]:
    reached: set[Function] = set()
    queue = deque(sorted(roots))
    while queue:
        function = queue.popleft()
        if function in reached:
            continue
        reached.add(function)
        queue.extend(sorted(graph[function] - reached))
    return reached


def completed_locations(plan: str) -> set[tuple[str, int]]:
    step_pattern = re.compile(
        r"^Step (\d+) \[x\]: (?P<action>.+)$", re.MULTILINE
    )
    source_pattern = re.compile(r"`(third_party/[^`]+):(\d+)`")
    result: set[tuple[str, int]] = set()
    for match in step_pattern.finditer(plan):
        if int(match.group(1)) > PREFIX_END_STEP:
            continue
        action = match.group("action")
        if not action.startswith("Port "):
            continue
        result.update(
            (path, int(line))
            for path, line in source_pattern.findall(action)
        )
    return result


def components(
    functions: set[Function], graph: dict[Function, set[Function]]
) -> list[list[Function]]:
    next_index = 0
    indices: dict[Function, int] = {}
    lowlinks: dict[Function, int] = {}
    stack: list[Function] = []
    on_stack: set[Function] = set()
    result: list[list[Function]] = []

    def visit(node: Function) -> None:
        nonlocal next_index
        indices[node] = next_index
        lowlinks[node] = next_index
        next_index += 1
        stack.append(node)
        on_stack.add(node)
        for callee in sorted(graph[node] & functions):
            if callee not in indices:
                visit(callee)
                lowlinks[node] = min(lowlinks[node], lowlinks[callee])
            elif callee in on_stack:
                lowlinks[node] = min(lowlinks[node], indices[callee])
        if lowlinks[node] != indices[node]:
            return
        component: list[Function] = []
        while True:
            member = stack.pop()
            on_stack.remove(member)
            component.append(member)
            if member == node:
                break
        result.append(sorted(component))

    for function in sorted(functions):
        if function not in indices:
            visit(function)
    return result


def callee_first_components(
    functions: set[Function], graph: dict[Function, set[Function]]
) -> list[list[Function]]:
    grouped = components(functions, graph)
    component_index = {
        function: index
        for index, component in enumerate(grouped)
        for function in component
    }
    dependencies: dict[int, set[int]] = defaultdict(set)
    for caller in functions:
        caller_index = component_index[caller]
        for callee in graph[caller] & functions:
            callee_index = component_index[callee]
            if caller_index != callee_index:
                dependencies[caller_index].add(callee_index)
    visited: set[int] = set()
    ordered: list[int] = []

    def visit(index: int) -> None:
        if index in visited:
            return
        visited.add(index)
        for dependency in sorted(dependencies[index]):
            visit(dependency)
        ordered.append(index)

    for index in range(len(grouped)):
        visit(index)
    return [grouped[index] for index in ordered]


def add_step(
    lines: list[str], number: int, checked: bool, action: str
) -> int:
    if "\n" in action:
        raise RuntimeError(f"step action must be one line: {action!r}")
    lines.append(f"Step {number} [{'x' if checked else ' '}]: {action}")
    return number + 1


def add_read(
    lines: list[str],
    number: int,
    *,
    checked: bool = False,
    core: bool = False,
) -> int:
    return add_step(
        lines,
        number,
        checked,
        CORE_READ_ACTION if core else READ_ACTION,
    )


def unit_names(functions: list[Function]) -> tuple[str, str, str]:
    if len(functions) == 1:
        function = functions[0]
        label = f"`{function.name}`"
        source = f"`{function.path}:{function.line}`"
        test_name = function.test_name
    else:
        label = " and ".join(f"`{item.name}`" for item in functions)
        source = "; ".join(
            f"`{item.path}:{item.line}`" for item in functions
        )
        joined = "_".join(normalize(item.name) for item in functions)
        lines = "_".join(str(item.line) for item in functions)
        test_name = f"source_port__scc__{joined}__lines_{lines}"
    return label, source, test_name


def add_official_unit(
    lines: list[str], number: int, functions: list[Function]
) -> int:
    label, source, test_name = unit_names(functions)
    oracle_name = (
        "official_c_oracle__"
        + "_".join(normalize(item.name) for item in functions)
        + "__exact"
    )
    number = add_read(lines, number)
    if len(functions) == 1:
        port_action = (
            f"Port {label} from {source} into `cosmolkit-inchi` with its "
            f"complete verbatim source frame, complete Rust behavior, active "
            f"macro behavior, source-defined integer, allocation, ownership, "
            f"output-mutation, error, and cleanup semantics, plus focused test "
            f"`{test_name}` covering every active branch; retain an `INCHI❗` "
            "first-axis marker until the official C oracle passes exactly, and "
            "do not use a heuristic, fallback, placeholder, partial branch, "
            "companion implementation, production FFI, or external command."
        )
    else:
        port_action = (
            f"Port the mutual-recursion SCC {label} from {source} as one "
            f"indivisible source-backed unit in `cosmolkit-inchi`, including "
            f"every member's complete verbatim source frame, complete Rust "
            f"behavior, active macro behavior, source-defined integer, "
            f"allocation, ownership, output-mutation, error, and cleanup "
            f"semantics, plus focused SCC test `{test_name}` covering every "
            f"active branch; retain `INCHI❗` first-axis markers until the "
            "official C oracle passes exactly, and do not use a heuristic, "
            "fallback, placeholder, partial member, companion implementation, "
            "production FFI, or external command."
        )
    number = add_step(lines, number, False, port_action)
    number = add_read(lines, number)
    number = add_step(
        lines,
        number,
        False,
        f"Run `cargo test -p cosmolkit-inchi {test_name}` and require the main "
        "test harness to execute at least one matching focused test with zero "
        "failures; a harness reporting only filtered-out tests is not a pass.",
    )
    number = add_read(lines, number)
    number = add_step(
        lines,
        number,
        False,
        f"Add non-production official C oracle test `{oracle_name}` for {label} "
        "that exercises every active branch and compares every return value, "
        "output field, mutated byte, pointer/offset relation, error text, and "
        "cleanup-visible state that the source exposes; production Rust must "
        "not call the oracle.",
    )
    number = add_read(lines, number)
    number = add_step(
        lines,
        number,
        False,
        f"Run `cargo test -p cosmolkit-inchi {oracle_name}` and require the "
        "main test harness to execute at least one matching exact-parity test "
        "with zero failures and no ignored fields, tolerance, compatibility "
        "assertion, skipped malformed case, or input-specific exception.",
    )
    number = add_read(lines, number)
    return add_step(
        lines,
        number,
        False,
        f"Update the {label} source marker and its row in "
        "`dev/gap_reports/rdkit_inchi_core_scope_audit.md` only after the "
        "focused Rust test and official C oracle both pass exactly; preserve "
        "the independently reviewed second-axis performance marker.",
    )


def add_adapter_unit(
    lines: list[str], number: int, function: AdapterFunction
) -> int:
    label = f"`{function.name}`"
    source = f"`{function.path}:{function.line}`"
    oracle_name = f"rdkit_cpp_oracle__{normalize(function.name)}__exact"
    number = add_read(lines, number)
    number = add_step(
        lines,
        number,
        False,
        f"Port {label} from {source} into the toolkit-neutral RDKit InChI "
        f"adapter with its complete verbatim source frame, complete control "
        f"flow, graph mutation order, integer, ownership, exception, logging, "
        f"stereo, isotope, hydrogen, charge, and cleanup semantics as "
        f"applicable, plus focused test `{function.test_name}` covering every "
        f"active GCC/Linux branch; retain a `RDKit❗` first-axis marker until "
        "the pinned RDKit C++ oracle passes exactly, and do not use a "
        "heuristic, fallback, placeholder, partial branch, MolBlock/SMILES "
        "transit, production FFI, or external command.",
    )
    number = add_read(lines, number)
    number = add_step(
        lines,
        number,
        False,
        f"Run `cargo test -p cosmolkit-inchi {function.test_name}` and require "
        "the main test harness to execute at least one matching focused test "
        "with zero failures; a harness reporting only filtered-out tests is "
        "not a pass.",
    )
    number = add_read(lines, number)
    number = add_step(
        lines,
        number,
        False,
        f"Add non-production pinned RDKit `2026.03.1` C++ oracle test "
        f"`{oracle_name}` for {label} that compares every observable return, "
        "diagnostic, graph mutation, atom/bond/stereo field, property, "
        "iteration-order effect, and exception path exposed by the source; "
        "production Rust must not call the oracle.",
    )
    number = add_read(lines, number)
    number = add_step(
        lines,
        number,
        False,
        f"Run `cargo test -p cosmolkit-inchi {oracle_name}` and require the "
        "main test harness to execute at least one matching exact-parity test "
        "with zero failures and no ignored field, tolerance, compatibility "
        "assertion, skipped malformed case, or molecule-specific exception.",
    )
    number = add_read(lines, number)
    return add_step(
        lines,
        number,
        False,
        f"Update the {label} source marker and its row in "
        "`dev/gap_reports/rdkit_inchi_core_scope_audit.md` only after the "
        "focused Rust test and pinned RDKit C++ oracle both pass exactly; "
        "preserve the independently reviewed second-axis performance marker.",
    )


def add_task(
    lines: list[str],
    number: int,
    action: str,
    *,
    core: bool = False,
) -> int:
    number = add_read(lines, number, core=core)
    return add_step(lines, number, False, action)


def markdown_table(
    rows: list[tuple[str, ...]], headers: tuple[str, ...]
) -> list[str]:
    result = [
        "| " + " | ".join(headers) + " |",
        "|" + "---|" * len(headers),
    ]
    result.extend("| " + " | ".join(row) + " |" for row in rows)
    return result


def write_report(
    selected: set[Function],
    graph: dict[Function, set[Function]],
    callback_edges: dict[Function, set[Function]],
    done: set[tuple[str, int]],
    missing_components: list[list[Function]],
    all_functions: set[Function],
) -> None:
    reverse: dict[Function, set[Function]] = defaultdict(set)
    for caller in selected:
        for callee in graph[caller] & selected:
            reverse[callee].add(caller)
    missing = {
        function for component in missing_components for function in component
    }
    completed_selected = {
        function for function in selected if function.location in done
    }
    excluded = all_functions - selected
    completed_excluded = {
        function for function in excluded if function.location in done
    }
    unresolved_excluded = excluded - completed_excluded
    callback_targets = {
        target
        for caller in selected
        for target in callback_edges.get(caller, set())
        if target in selected
    }
    lines = [
        "# RDKit InChI Core Scope Audit",
        "",
        "## Scope Decision",
        "",
        "This audit replaces the all-export `libinchi` port objective with the "
        "exact engine and adapter closure required by four scalar RDKit APIs:",
        "",
        "- `Chem.MolToInchi(mol)`;",
        "- `Chem.MolToInchiKey(mol)`;",
        "- `InchiToInchiKey(inchi)`;",
        "- `Chem.MolFromInchi(inchi)`.",
        "",
        "MolBlock, SDF, V3000 Molfile, IXA, AuxInfo conversion, incremental "
        "INCHIGEN, version-query, CLI/demo, and unrelated export parity are "
        "outside this scope. Already ported outside-scope code remains private "
        "and frozen; this plan neither deletes it nor spends further parity "
        "work on it.",
        "",
        "## Compiler Evidence",
        "",
        "- Official target: GCC/Linux `libinchi` with `COMPILE_ANSI_ONLY`, "
        "`TARGET_API_LIB`, and `fPIC` from the official CMake target.",
        "- Definition and direct-edge evidence: the configured 60-entry "
        "`compile_commands.json`, GCC `-aux-info`, and GCC "
        "`-fdump-ipa-cgraph` artifacts.",
        "- Callback evidence: active GCC IPA `References: ... (addr)` records "
        "from selected callers; no function-name-only inference is used.",
        f"- Active configured definitions: `{len(all_functions)}`.",
        "- Selected official roots: `5`.",
        "- Compiler direct-edge selected closure: `1014` definitions.",
        f"- Callback-complete selected closure: `{len(selected)}` definitions.",
        f"- Compiler address-taken callback edges in the selected closure: "
        f"`{sum(len(value) for value in callback_edges.values())}`.",
        f"- Address-taken callback target definitions added to the direct "
        f"closure: `{len(selected) - 1014}`.",
        f"- Completed selected source locations: `{len(completed_selected)}`.",
        f"- Remaining official source functions: `{len(missing)}` in "
        f"`{len(missing_components)}` callee-first Port units.",
        f"- Active completed functions outside scope and frozen: "
        f"`{len(completed_excluded)}`.",
        f"- Active functions outside scope and not scheduled: "
        f"`{len(unresolved_excluded)}`.",
        "",
        "## Root Mapping",
        "",
    ]
    root_rows = [
        (
            "`Chem.MolToInchi`",
            "`GetINCHI`, `FreeINCHI`",
            "`inchi.cpp:2082`, `inchi.cpp:2100`",
        ),
        (
            "`Chem.MolToInchiKey`",
            "`GetINCHI`, `FreeINCHI`, `GetINCHIKeyFromINCHI`",
            "`inchi.h:107` composition",
        ),
        (
            "`InchiToInchiKey`",
            "`GetINCHIKeyFromINCHI`",
            "`inchi.cpp:2149`",
        ),
        (
            "`Chem.MolFromInchi`",
            "`GetStructFromINCHI`, `FreeStructFromINCHI`",
            "`inchi.cpp:1273`, `inchi.cpp:1371`, `inchi.cpp:1648`",
        ),
    ]
    lines.extend(
        markdown_table(
            root_rows, ("RDKit API", "Official roots", "RDKit evidence")
        )
    )
    lines.extend(
        [
            "",
            "## State Constraints",
            "",
            "- `GetINCHI` wraps plain `inchi_Input`, sets "
            "`extended_input.polymer = NULL` and "
            "`extended_input.v3000 = NULL`, then calls `GetINCHI1`.",
            "- `GetStructFromINCHI` uses the InChI string path and forces "
            "`INPUT_INCHI`; it does not consume a Molfile.",
            "- User option strings remain part of supported RDKit behavior, so "
            "their reachable normalization, canonicalization, reconstruction, "
            "stereo, FixedH, reconnected, warning, and error paths are retained.",
            "- Global `UserAction` and `ConsoleQuit` callbacks are initialized "
            "to null and have no configured RDKit API setter; no invented "
            "callback target is added for them.",
            "- Macro-only allocation behavior remains attached to callers and "
            "does not create fake functions for inactive `util.c` allocator "
            "bodies.",
            "- The abandoned import-only staging for the unimplemented "
            "`ReadMolfileToInpAtoms` step is not a completed port and is "
            "scheduled for removal before the first remaining core function; "
            "no completed Molfile/V3000 implementation is removed.",
            "",
            "## Remaining Official Functions",
            "",
        ]
    )
    remaining_rows: list[tuple[str, ...]] = []
    for index, component in enumerate(missing_components, 1):
        for function in component:
            callers = ", ".join(
                f"`{item.name}`" for item in sorted(reverse[function])
            )
            kind = (
                "mutual-recursion SCC"
                if len(component) > 1
                else "callback target"
                if function in callback_targets
                else "direct callee"
            )
            remaining_rows.append(
                (
                    str(index),
                    f"`{function.name}`",
                    f"`{function.path}:{function.line}`",
                    kind,
                    callers or "selected root",
                )
            )
    lines.extend(
        markdown_table(
            remaining_rows,
            ("Port unit", "Function", "Source", "Reason", "Selected callers"),
        )
    )
    if any(len(component) > 1 for component in missing_components):
        lines.extend(
            [
                "",
                "Every listed mutual-recursion SCC is one atomic Port unit; "
                "every other row is an individual function Port unit.",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "No remaining function belongs to an unfinished SCC. The "
                "completed `CheckINCHI`/`GetINCHIfromINCHI` SCC remains in the "
                "checked prefix and is not rescheduled.",
            ]
        )
    lines.extend(
        [
            "",
            "## Remaining RDKit Adapter Functions",
            "",
        ]
    )
    adapter_rows = [
        (
            str(index),
            f"`{item.name}`",
            f"`{item.path}:{item.line}`",
            item.responsibility,
        )
        for index, item in enumerate(ADAPTER_FUNCTIONS, 1)
    ]
    lines.extend(
        markdown_table(
            adapter_rows, ("Order", "Function", "Source", "Core responsibility")
        )
    )
    lines.extend(
        [
            "",
            "`MolBlockToInchi` and `getInchiVersion` are intentionally absent. "
            "The other 30 adapter functions form the source-level helper and "
            "public-entry closure of the four selected APIs.",
            "",
            "## Frozen And Excluded Active Definitions",
            "",
            "Every configured active definition outside the five-root "
            "callback-complete closure is frozen or unscheduled below. A "
            "`completed-frozen` row keeps its existing Rust code but receives "
            "no further Port or parity step. An `excluded-unported` row is not "
            "part of the new plan.",
            "",
        ]
    )
    excluded_rows = [
        (
            f"`{function.name}`",
            f"`{function.path}:{function.line}`",
            "completed-frozen"
            if function in completed_excluded
            else "excluded-unported",
        )
        for function in sorted(excluded)
    ]
    lines.extend(
        markdown_table(excluded_rows, ("Function", "Source", "Disposition"))
    )
    active_locations = {function.location for function in all_functions}
    inactive_historical = sorted(done - active_locations)
    lines.extend(
        [
            "",
            "## Historical Non-Function Anchors",
            "",
            "The following checked historical source locations do not "
            "correspond to active GCC function bodies under the selected "
            "target. They remain historical records and are not rescheduled:",
            "",
        ]
    )
    lines.extend(
        f"- `{path}:{line}`" for path, line in inactive_historical
    )
    lines.extend(
        [
            "",
            "## Completion Rule",
            "",
            f"No selected public API may be exposed until all {len(missing)} official "
            "functions and all 30 RDKit adapter functions have complete source "
            "frames, focused tests that actually execute, exact non-production "
            "oracles, and closed first-axis markers. Performance markers remain "
            "independent. Exact parity does not permit ignored fields, "
            "tolerances, compatibility assertions, skipped malformed inputs, "
            "or molecule-specific exceptions.",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_plan(
    original: str,
    missing_components: list[list[Function]],
) -> str:
    boundary = re.search(
        rf"^Step {PREFIX_END_STEP + 1} \[[ x]\]:.*$",
        original,
        re.MULTILINE,
    )
    if boundary is None:
        raise RuntimeError(f"Step {PREFIX_END_STEP + 1} not found")
    prefix = original[: boundary.start()].rstrip()
    lines = [prefix, ""]
    number = PREFIX_END_STEP + 1
    number = add_step(
        lines,
        number,
        True,
        "Read `dev/agent_plan_standard.md` and follow its one-action, "
        "mandatory-reading, immediate-focused-test, sequential-execution, "
        "artifact, and no-git requirements for the scope repair.",
    )
    number = add_read(lines, number, checked=True)
    number = add_step(
        lines,
        number,
        True,
        "Audit the configured GCC/Linux five-root official InChI closure, "
        "compiler address-taken callback closure, completed source locations, "
        "and four-API RDKit adapter closure and write "
        "`dev/gap_reports/rdkit_inchi_core_scope_audit.md` with exact remaining "
        "and frozen classifications.",
    )
    number = add_read(lines, number, checked=True)
    number = add_step(
        lines,
        number,
        True,
        "Rewrite every unfinished step from Step 4721 onward as the "
        "callee-first five-root core plan while preserving Steps 1-4720 "
        "verbatim, freezing completed outside-scope code, excluding unported "
        "outside-scope definitions, retaining mutual-recursion SCC atomicity, "
        "and scheduling each remaining official and RDKit adapter function "
        "with focused and exact-oracle validation.",
    )
    lines.extend(
        [
            "",
            "## RDKit InChI Core Execution Contract",
            "",
            "- This suffix replaces full `libinchi` export parity with the exact closure of `GetINCHI`, `FreeINCHI`, `GetStructFromINCHI`, `FreeStructFromINCHI`, and `GetINCHIKeyFromINCHI`, plus the RDKit adapter closure for the four selected scalar APIs.",
            "- Execute unchecked steps one at a time in order and continue until all steps are completed, blocked by an exact-parity failure, or the user interrupts.",
            "- Change only the corresponding `[ ]` to `[x]` after a step produces its required artifact and satisfies its stated result.",
            "- Every real task is immediately preceded by a fresh, actual reading of `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`; core-operation tasks also read `dev/README.md` in that same reading step.",
            "- A focused test command passes only when the main test harness actually executes at least one matching test and reports zero failures; output containing only filtered-out tests is not evidence.",
            "- Official C and pinned RDKit C++ code are non-production test oracles only; production Rust has no FFI, subprocess, or external-command dependency.",
            "- No heuristic, placeholder, partial port, compatibility threshold, ignored field, skipped malformed input, molecule-specific exception, SMILES/MolBlock transit, operation-local companion, duplicated callee, or silent fallback is permitted.",
            "- Already ported definitions outside the selected closure remain private and frozen. Do not delete, expose, extend, test for new parity claims, or use them to justify completion of a selected function.",
            "- Do not expose any selected public API until its complete official and RDKit source closure and exact aggregate oracles pass.",
            "",
            "## Remaining Official Engine Plan",
            "",
        ]
    )
    number = add_task(
        lines,
        number,
        "Remove only the abandoned `ReadMolfileToInpAtoms` import staging from "
        "`crates/cosmolkit-inchi/src/source/base/mol2atom.rs`, leaving every "
        "completed out-of-scope Molfile/V3000 function and source frame "
        "unchanged and frozen.",
    )
    number = add_task(
        lines,
        number,
        "Run `cargo check -p cosmolkit-inchi --tests` and require zero errors "
        "after removing the abandoned import-only staging.",
    )
    for component in missing_components:
        number = add_official_unit(lines, number, component)
    number = add_task(
        lines,
        number,
        "Add the non-production aggregate official C root oracle "
        "`official_c_oracle__rdkit_core_roots__exact` for `GetINCHI`, "
        "`FreeINCHI`, `GetStructFromINCHI`, `FreeStructFromINCHI`, and "
        "`GetINCHIKeyFromINCHI`, covering success, warning, error, malformed "
        "input, allocation failure, option branches, ordinary/FixedH/"
        "reconnected layers, isotope, charge, radicals, 0D/2D/3D stereo, "
        "coordinates, ownership, output-reset, and cleanup behavior with "
        "field-for-field and byte-for-byte Rust comparison.",
    )
    number = add_task(
        lines,
        number,
        "Run `cargo test -p cosmolkit-inchi "
        "official_c_oracle__rdkit_core_roots__exact` and require the main test "
        "harness to execute at least one matching test with exact parity and "
        "zero failures.",
    )
    lines.extend(["", "## Remaining RDKit Adapter Plan", ""])
    for function in ADAPTER_FUNCTIONS:
        number = add_adapter_unit(lines, number, function)
    number = add_task(
        lines,
        number,
        "Add the pinned RDKit `2026.03.1` aggregate scalar oracle "
        "`rdkit_inchi_core_scalar_exact` for `Chem.MolToInchi`, "
        "`Chem.MolToInchiKey`, `InchiToInchiKey`, and `Chem.MolFromInchi`, "
        "covering graph ordering, aromatic/kekule handling, all supported "
        "atom/bond/stereo/isotope/H/charge/radical states, coordinates, "
        "cleanup branches, options, diagnostics, malformed input, sanitize "
        "failure, unchanged source molecules, and concurrent calls.",
    )
    number = add_task(
        lines,
        number,
        "Run `cargo test -p cosmolkit-inchi rdkit_inchi_core_scalar_exact` in "
        "release mode and require exact output, graph-state, diagnostic, and "
        "error parity for every case with zero ignored fields or failures.",
    )
    number = add_task(
        lines,
        number,
        "Audit the completed official and RDKit closures and update "
        "`dev/gap_reports/rdkit_inchi_core_scope_audit.md` with every focused "
        "command, oracle command, exact result, final first-axis marker, "
        "preserved second-axis marker, and any remaining branch; stop if any "
        "selected branch remains open.",
    )
    lines.extend(["", "## Core Integration And Validation Plan", ""])
    number = add_task(
        lines,
        number,
        "Expose only the four completed toolkit-neutral scalar operations "
        "`mol_to_inchi`, `mol_to_inchi_key`, `inchi_to_inchi_key`, and "
        "`mol_from_inchi` with structured return diagnostics and explicit "
        "unsupported errors for inputs outside the audited state space.",
    )
    number = add_task(
        lines,
        number,
        "Add focused toolkit-neutral scalar API tests named "
        "`inchi_core_scalar_api` for success, warnings, errors, ownership, "
        "unchanged input, stereo, isotope, coordinates, and option behavior.",
    )
    number = add_task(
        lines,
        number,
        "Run `cargo test -p cosmolkit-inchi inchi_core_scalar_api` and require "
        "the main test harness to execute matching tests with zero failures.",
    )
    number = add_task(
        lines,
        number,
        "Implement the COSMolKit `Molecule` bridge for the four scalar InChI "
        "operations through narrowed read parts and registered operation "
        "machinery where topology mutation is required, preserving graph, "
        "coordinate, stereo, property, cache, copy-on-write, and structured "
        "error invariants.",
        core=True,
    )
    number = add_task(
        lines,
        number,
        "Add strict focused core bridge tests named `inchi_core_bridge` for "
        "graph index validity, coordinate alignment, atom/bond/stereo mapping, "
        "property behavior, cache invalidation, copy-on-write source "
        "immutability, unsupported-state errors, and exact RDKit parity.",
        core=True,
    )
    number = add_task(
        lines,
        number,
        "Run `cargo test -p cosmolkit-core --release --features "
        "op-contracts-strict inchi_core_bridge` and require matching tests to "
        "execute with zero failures.",
        core=True,
    )
    number = add_task(
        lines,
        number,
        "Expose the four completed value-style Rust facade APIs and add exact "
        "facade tests named `inchi_rust_facade` without exposing frozen "
        "Molfile, IXA, AuxInfo, INCHIGEN, or version-query APIs.",
    )
    number = add_task(
        lines,
        number,
        "Run `cargo test --workspace --release --features "
        "cosmolkit-core/op-contracts-strict inchi_rust_facade` and require "
        "matching tests to execute with zero failures.",
        core=True,
    )
    number = add_task(
        lines,
        number,
        "Expose the four completed scalar APIs through Python with generated "
        "stubs and add exact tests for `Chem.MolToInchi`, "
        "`Chem.MolToInchiKey`, `InchiToInchiKey`, and `Chem.MolFromInchi` "
        "semantics, diagnostics, options, and errors.",
    )
    number = add_task(
        lines,
        number,
        "Run `.venv/bin/maturin develop --manifest-path python/Cargo.toml` and "
        "then `.venv/bin/pytest python/tests/test_inchi.py`, requiring all "
        "selected InChI tests to execute with zero failures.",
    )
    number = add_task(
        lines,
        number,
        "Update the structured support registry, parity matrix, provenance, "
        "Rust documentation, Python documentation, and test-scope documents "
        "to claim support only for the four audited scalar APIs and to record "
        "MolBlock, SDF/V3000, IXA, AuxInfo, INCHIGEN, version query, and "
        "extended polymer input as frozen or unsupported as applicable.",
        core=True,
    )
    number = add_task(
        lines,
        number,
        "Run `cargo fmt --all -- --check`.",
    )
    number = add_task(
        lines,
        number,
        "Run `cargo check -p cosmolkit-core --features op-contracts-strict`.",
        core=True,
    )
    number = add_task(
        lines,
        number,
        "Run `cargo test -p cosmolkit-core --release --features "
        "op-contracts-strict`.",
        core=True,
    )
    number = add_task(
        lines,
        number,
        "Run `cargo test --workspace --release --features "
        "cosmolkit-core/op-contracts-strict`.",
        core=True,
    )
    number = add_task(
        lines,
        number,
        "Run `.venv/bin/python -m sphinx -b html python/docs/source "
        "python/docs/build/html` and `.venv/bin/basedpyright python/tests "
        "python/examples`, requiring both commands to complete without errors.",
    )
    number = add_task(
        lines,
        number,
        "Write `dev/gap_reports/rdkit_inchi_core_final_validation.md` with "
        "separate source-frame, focused-Rust-test, official-C-oracle, pinned-"
        "RDKit-oracle, remaining-selected-call-graph, frozen-outside-scope, "
        "unsupported-branch, performance-marker, and full-validation results; "
        "report exact parity only where every compared field and branch passes "
        "exactly.",
    )
    content = "\n".join(lines) + "\n"
    validate_plan(content, prefix, missing_components)
    return content


def validate_plan(
    content: str,
    prefix: str,
    missing_components: list[list[Function]],
) -> None:
    if not content.startswith(prefix):
        raise RuntimeError(f"Step 1-{PREFIX_END_STEP} prefix changed")
    step_pattern = re.compile(
        r"^Step (\d+) \[([ x])\]: (.+)$", re.MULTILINE
    )
    steps = [
        (int(match.group(1)), match.group(2), match.group(3))
        for match in step_pattern.finditer(content)
    ]
    if [item[0] for item in steps] != list(range(1, len(steps) + 1)):
        raise RuntimeError("step numbers are not contiguous")
    for number, checked, action in steps:
        if "\n" in action:
            raise RuntimeError(f"Step {number} is not one line")
        if number <= PREFIX_END_STEP + 5 and checked != "x":
            raise RuntimeError(f"Step {number} must be checked")
        if number > PREFIX_END_STEP + 5 and checked != " ":
            raise RuntimeError(f"Step {number} must remain unchecked")
    for index, (number, _, action) in enumerate(steps):
        if number <= PREFIX_END_STEP + 1:
            continue
        is_read = action.startswith("Read `dev/policy_invariants.md`")
        if is_read:
            continue
        previous_action = steps[index - 1][2]
        if not previous_action.startswith("Read `dev/policy_invariants.md`"):
            raise RuntimeError(
                f"Step {number} is not immediately preceded by mandatory reading"
            )
    suffix = content[len(prefix) :]
    for component in missing_components:
        for function in component:
            marker = f"`{function.path}:{function.line}`"
            if suffix.count(marker) != 1:
                raise RuntimeError(
                    f"remaining official function not represented once: {marker}"
                )
    for function in ADAPTER_FUNCTIONS:
        marker = f"`{function.path}:{function.line}`"
        if suffix.count(marker) != 1:
            raise RuntimeError(
                f"remaining adapter function not represented once: {marker}"
            )
    forbidden = (
        "Port all ",
        "Port subsystem",
        "Port file",
        "Use complete source-framed operation-local implementations",
        "Promote the operation-local",
    )
    for phrase in forbidden:
        if phrase in suffix:
            raise RuntimeError(f"forbidden plan phrase present: {phrase}")


def regenerate(dry_run: bool) -> None:
    original = PLAN_PATH.read_text(encoding="utf-8")
    done = completed_locations(original)
    functions_by_key, graph, callback_edges = load_graph()
    all_functions = set(functions_by_key.values())
    roots = {
        function
        for function in all_functions
        if function.name in set(ROOT_NAMES)
    }
    selected = closure(roots, graph)
    missing = {
        function for function in selected if function.location not in done
    }
    ordered = callee_first_components(selected, graph)
    missing_components = [
        [function for function in component if function in missing]
        for component in ordered
        if any(function in missing for function in component)
    ]
    for component, missing_component in zip(
        [item for item in ordered if any(f in missing for f in item)],
        missing_components,
    ):
        if len(component) > 1 and len(component) != len(missing_component):
            raise RuntimeError(
                "partially completed selected SCC: "
                + ", ".join(item.key for item in component)
            )
    if len(selected) != 1050:
        raise RuntimeError(
            f"expected 1050 callback-complete selected definitions, found {len(selected)}"
        )
    if len(missing) != 23:
        raise RuntimeError(
            f"expected 23 remaining official functions, found {len(missing)}"
        )
    if len(ADAPTER_FUNCTIONS) != 30:
        raise RuntimeError(
            f"expected 30 selected adapter functions, found {len(ADAPTER_FUNCTIONS)}"
        )
    plan = generate_plan(original, missing_components)
    print(f"selected official closure: {len(selected)}")
    print(f"completed selected functions: {len(selected) - len(missing)}")
    print(f"remaining official functions: {len(missing)}")
    print(f"remaining official Port units: {len(missing_components)}")
    print(f"remaining RDKit adapter functions: {len(ADAPTER_FUNCTIONS)}")
    print(
        "active completed outside scope: "
        f"{sum(item.location in done for item in all_functions - selected)}"
    )
    generated_steps = re.findall(r"^Step (\d+)", plan, re.MULTILINE)
    print(f"generated final step: {max(int(value) for value in generated_steps)}")
    if not dry_run:
        write_report(
            selected,
            graph,
            callback_edges,
            done,
            missing_components,
            all_functions,
        )
        PLAN_PATH.write_text(plan, encoding="utf-8")
        print(f"wrote {REPORT_PATH.relative_to(REPO_ROOT)}")
        print(f"wrote {PLAN_PATH.relative_to(REPO_ROOT)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    regenerate(args.dry_run)


if __name__ == "__main__":
    main()
