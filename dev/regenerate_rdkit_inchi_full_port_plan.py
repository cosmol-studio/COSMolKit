#!/usr/bin/env python3
"""Regenerate the unchecked InChI port-plan suffix from the GCC active audit."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
PLAN_PATH = REPO_ROOT / "dev" / "rdkit_inchi_full_port_plan.md"
AUDIT_PATH = (
    REPO_ROOT
    / "target"
    / "inchi-active-call-graph-audit"
    / "official_inchi_active_call_graph.json"
)
RDKIT_INVENTORY_PATH = (
    REPO_ROOT / "dev" / "gap_reports" / "rdkit_inchi_source_inventory.md"
)
READ_ACTION = (
    "Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` "
    "to reload and follow the required execution standard, source reproduction "
    "rules, artifact requirements, no-git rule, and completion criteria for the "
    "next task."
)


@dataclass(frozen=True, order=True)
class SourceFunction:
    key: str
    name: str
    path: str
    line: int
    reachable: bool
    test_name: str

    @property
    def location(self) -> tuple[str, int]:
        return self.path, self.line


def parse_active_audit() -> tuple[list[SourceFunction], dict[SourceFunction, set[SourceFunction]]]:
    data = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    functions = [
        SourceFunction(
            key=item["key"],
            name=item["name"],
            path=item["path"],
            line=item["line"],
            reachable=item["reachable"],
            test_name=item["test_name"],
        )
        for item in data["functions"]
    ]
    by_key = {function.key: function for function in functions}
    graph = {
        by_key[item["key"]]: {by_key[key] for key in item["callees"]}
        for item in data["functions"]
    }
    if len(functions) != 1364:
        raise RuntimeError(f"expected 1364 configured active definitions, found {len(functions)}")
    return functions, graph


def parse_rdkit_inventory() -> list[SourceFunction]:
    text = RDKIT_INVENTORY_PATH.read_text(encoding="utf-8")
    start = text.index("## Adapter Function Inventory")
    end = text.index("## Official Engine Dependencies")
    row = re.compile(
        r"^\| `(?P<name>[^`]+)` \| `(?P<path>inchi\.(?:cpp|h)):(?P<line>\d+)` \|",
        re.MULTILINE,
    )
    functions: list[SourceFunction] = []
    for match in row.finditer(text[start:end]):
        name = match.group("name")
        path = "third_party/rdkit/External/INCHI-API/" + match.group("path")
        line = int(match.group("line"))
        stem = Path(path).stem.lower()
        normalized = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
        functions.append(
            SourceFunction(
                key=f"{path}:{line}:{name}",
                name=name,
                path=path,
                line=line,
                reachable=True,
                test_name=f"source_port__{stem}__{normalized}__line_{line}",
            )
        )
    if len(functions) != 32:
        raise RuntimeError(f"expected 32 RDKit adapter definitions, found {len(functions)}")
    return functions


def components(graph: dict[SourceFunction, set[SourceFunction]]) -> list[list[SourceFunction]]:
    next_index = 0
    indices: dict[SourceFunction, int] = {}
    lowlinks: dict[SourceFunction, int] = {}
    stack: list[SourceFunction] = []
    on_stack: set[SourceFunction] = set()
    result: list[list[SourceFunction]] = []

    def visit(node: SourceFunction) -> None:
        nonlocal next_index
        indices[node] = next_index
        lowlinks[node] = next_index
        next_index += 1
        stack.append(node)
        on_stack.add(node)
        for callee in sorted(graph[node]):
            if callee not in indices:
                visit(callee)
                lowlinks[node] = min(lowlinks[node], lowlinks[callee])
            elif callee in on_stack:
                lowlinks[node] = min(lowlinks[node], indices[callee])
        if lowlinks[node] != indices[node]:
            return
        component: list[SourceFunction] = []
        while True:
            member = stack.pop()
            on_stack.remove(member)
            component.append(member)
            if member == node:
                break
        result.append(sorted(component))

    for function in sorted(graph):
        if function not in indices:
            visit(function)
    return result


def callee_first_components(
    graph: dict[SourceFunction, set[SourceFunction]],
) -> list[list[SourceFunction]]:
    grouped = components(graph)
    component_index = {
        function: index for index, component in enumerate(grouped) for function in component
    }
    dependencies: dict[int, set[int]] = defaultdict(set)
    for caller, callees in graph.items():
        caller_index = component_index[caller]
        for callee in callees:
            callee_index = component_index[callee]
            if caller_index != callee_index:
                dependencies[caller_index].add(callee_index)
    ordered: list[int] = []
    visited: set[int] = set()

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


def add_step(lines: list[str], number: int, checked: bool, action: str) -> int:
    lines.append(f"Step {number} [{'x' if checked else ' '}]: {action}")
    return number + 1


def add_read(lines: list[str], number: int, checked: bool = False) -> int:
    return add_step(lines, number, checked, READ_ACTION)


def add_port_and_test(
    lines: list[str], number: int, functions: list[SourceFunction], scope: str
) -> int:
    if len(functions) == 1:
        function = functions[0]
        reachability = (
            "configured export-root-reachable production behavior"
            if function.reachable
            else "configured active but export-unreachable target behavior; keep it private"
        )
        action = (
            f"Port `{function.name}` from `{function.path}:{function.line}` as {reachability} "
            f"into {scope} with its complete verbatim source frame, complete Rust behavior, "
            f"source-defined integer and ownership semantics, and focused test named "
            f"`{function.test_name}`; no companion copy, heuristic, fallback, placeholder, or "
            "partial implementation is permitted."
        )
        test_name = function.test_name
    else:
        members = "; ".join(
            f"`{function.name}` from `{function.path}:{function.line}`" for function in functions
        )
        first = functions[0]
        normalized = re.sub(r"[^a-z0-9]+", "_", first.name.lower()).strip("_")
        test_name = f"source_port__scc__{normalized}__line_{first.line}"
        action = (
            f"Port the configured mutual-recursion SCC as one atomic source-backed unit: {members}; "
            f"include every member's complete verbatim source frame and complete Rust behavior in "
            f"{scope}, plus focused SCC test `{test_name}` covering every member; no companion copy, "
            "promotion, heuristic, fallback, placeholder, or partial member is permitted."
        )
    number = add_step(lines, number, False, action)
    number = add_read(lines, number)
    return add_step(lines, number, False, f"Run `cargo test -p cosmolkit-inchi {test_name}`.")


def completed_locations(prefix: str) -> set[tuple[str, int]]:
    row = re.compile(
        r"^Step (\d+) \[x\]: Port `[^`]+` from `([^`]+):(\d+)`",
        re.MULTILINE,
    )
    return {
        (match.group(2), int(match.group(3)))
        for match in row.finditer(prefix)
        if int(match.group(1)) <= 104
    }


def validate(content: str, prefix: str, active: list[SourceFunction]) -> None:
    if not content.startswith(prefix):
        raise RuntimeError("Step 1-104 prefix changed")
    step = re.compile(r"^Step (\d+) \[([ x])\]: (.+)$", re.MULTILINE)
    parsed = [(int(item.group(1)), item.group(2), item.group(3)) for item in step.finditer(content)]
    if [item[0] for item in parsed] != list(range(1, len(parsed) + 1)):
        raise RuntimeError("step numbers are not contiguous")
    for number, checked, _ in parsed:
        if number <= 107 and checked != "x":
            raise RuntimeError(f"Step {number} must be checked")
        if number > 107 and checked != " ":
            raise RuntimeError(f"Step {number} must remain unchecked")
    done = completed_locations(prefix)
    remaining_locations = {item.location for item in active if item.location not in done}
    port_locations: set[tuple[str, int]] = set()
    suffix = content[len(prefix) :]
    for function in active:
        marker = f"`{function.name}` from `{function.path}:{function.line}`"
        if function.location in remaining_locations and suffix.count(marker) != 1:
            raise RuntimeError(f"active function is not represented exactly once: {marker}")
        if function.location in done and marker in suffix:
            raise RuntimeError(f"completed function was rescheduled: {marker}")
        if function.location in remaining_locations:
            port_locations.add(function.location)
    if port_locations != remaining_locations:
        raise RuntimeError("active port coverage mismatch")
    if (
        "Use complete source-framed operation-local implementations" in suffix
        or "Promote the operation-local" in suffix
    ):
        raise RuntimeError("forbidden companion/promotion mechanism in regenerated suffix")


def regenerate(dry_run: bool) -> None:
    active, graph = parse_active_audit()
    rdkit = parse_rdkit_inventory()
    original = PLAN_PATH.read_text(encoding="utf-8")
    step_105 = re.search(r"^Step 105 \[[ x]\]:.*$", original, re.MULTILINE)
    if step_105 is None:
        raise RuntimeError("Step 105 not found")
    prefix = original[: step_105.start()]
    done = completed_locations(prefix)

    by_name = {function.name: function for function in active}
    forced = [
        by_name[name]
        for name in (
            "already_have_this_message",
            "AddErrorMessage",
            "is_in_the_slist",
            "is_element_a_metal",
            "InchiToInchiAtom",
        )
    ]
    if any(function.location in done for function in forced):
        raise RuntimeError("forced repaired-prefix function is already marked complete")

    ordered = callee_first_components(graph)
    remaining_components: list[list[SourceFunction]] = []
    forced_set = set(forced)
    for component in ordered:
        unfinished = [
            function
            for function in component
            if function.location not in done and function not in forced_set
        ]
        if unfinished and len(unfinished) != len([item for item in component if item not in forced_set]):
            raise RuntimeError(
                "partially completed SCC: " + ", ".join(item.key for item in component)
            )
        if unfinished:
            remaining_components.append(unfinished)

    lines = [prefix.rstrip(), ""]
    number = 105
    number = add_step(
        lines,
        number,
        True,
        "Write the configured GCC/Linux `COMPILE_ANSI_ONLY` and `TARGET_API_LIB` active call graph audit to `dev/gap_reports/official_inchi_active_call_graph_audit.md`, including inactive, macro-only, header-defined, external, libc/intrinsic, CLI/demo, and export-unreachable classifications plus the complete compiler-derived `InchiToInchiAtom` callee closure.",
    )
    number = add_read(lines, number, checked=True)
    number = add_step(
        lines,
        number,
        True,
        "Rewrite every unfinished step from Step 105 onward from the configured compiler audit while preserving Step 1-104 text and completion markers exactly, excluding inactive-only function bodies, representing macro/header behavior without fake C functions, ordering active definitions callee-first, and making every mutual-recursion SCC an atomic source-backed unit.",
    )
    lines.extend(
        [
            "",
            "## Repaired Active-Graph Execution Contract",
            "",
            "- The GCC audit report and its machine-readable call graph are the source of truth for configured activity and direct call edges.",
            "- Steps 1-104 remain historical completed work. Their text and markers are unchanged; inactive `util.c` allocator anchors are corrected only through the explicit active `mode.h` macro behavior in the `InchiToInchiAtom` closure step.",
            "- Every active function has one explicit source-location Port step, except that every mutual-recursion SCC is one indivisible Port step listing every member and source location.",
            "- Every Port step includes complete source frames, complete Rust behavior, and a focused test, followed immediately by a mandatory reading step and the most specific test command.",
            "- Active but export-unreachable target functions stay private. Inactive-only raw bodies, declarations, libc/compiler intrinsics, CLI/demo-only behavior, and macro names do not become fake production function ports.",
            "- No operation-local companion, copied callee, later promotion, heuristic, SMILES/MolBlock transit, subprocess, production FFI, fallback, placeholder, or partial port is permitted.",
            "- No production API is exposed before its entire active reachable closure and exact official-C behavior oracle pass.",
            "",
            "## Repaired Plan",
            "",
        ]
    )

    for prerequisite in forced[:-1]:
        number = add_read(lines, number)
        number = add_port_and_test(lines, number, [prerequisite], "`cosmolkit-inchi`")
    number = add_read(lines, number)

    inchi = forced[-1]
    inchi_action = (
        f"Port `{inchi.name}` from `{inchi.path}:{inchi.line}` into `cosmolkit-inchi` with its "
        "complete verbatim source frame and complete Rust behavior; correct the closure's allocator "
        "anchors to the active `mode.h:1172-1181` malloc/calloc/realloc/free macros rather than the "
        "inactive `util.c` bodies; include focused test "
        f"`{inchi.test_name}` covering every active source branch: null/invalid input, "
        "`INCHI_IOSTREAM` string and file modes, single-line and cross-buffer reads, "
        "`FindToken`/`LoadLine`, SDF label/value, atom count and maximum capacity, bond count/type, "
        "2D/3D/undefined coordinates, stereo0D allocation/fill/free, ordinary and isotope H, formal "
        "charge, metal classification, parity/stereo, malformed input, fatal error, warning/error "
        "text, partial-allocation cleanup, and source-defined integer boundaries; no heuristic, "
        "fallback, placeholder, partial branch, or copied callee is permitted."
    )
    number = add_step(lines, number, False, inchi_action)
    number = add_read(lines, number)
    number = add_step(
        lines,
        number,
        False,
        f"Run `cargo test -p cosmolkit-inchi {inchi.test_name}`.",
    )
    number = add_read(lines, number)
    number = add_step(
        lines,
        number,
        False,
        "Add a non-production official C behavior oracle fixture for fixed `InchiToInchiAtom` inputs that emits return status, error code/text, atom and bond counts, every atom and bond field, coordinate dimensions, every stereo0D field, isotope-H, charge, and SDF labels/values; provenance-only checks are insufficient and production Rust must not call the oracle.",
    )
    number = add_read(lines, number)
    number = add_step(
        lines,
        number,
        False,
        "Run the official C `InchiToInchiAtom` behavior oracle and `cargo test -p cosmolkit-inchi official_c_oracle_inchi_to_inchi_atom_exact`, requiring Rust and official C to match every emitted field exactly with no ignored field, compatibility assertion, lowered threshold, exception skip, or molecule/line-number special case.",
    )

    for component in remaining_components:
        number = add_read(lines, number)
        number = add_port_and_test(lines, number, component, "`cosmolkit-inchi`")

    for function in rdkit:
        number = add_read(lines, number)
        number = add_port_and_test(lines, number, [function], "the toolkit-neutral adapter")

    integration_actions = [
        "Expose the completed toolkit-neutral scalar APIs only after every active reachable official and RDKit adapter closure has complete source frames, focused tests, and required exact official-C fixtures.",
        "Add complete scalar public-API tests for success, warning, error, ownership, and unchanged-input behavior.",
        "Run the complete `cosmolkit-inchi` scalar public-API tests.",
        "Add `INCHI_FEATURE` to the structured core support registry as unsupported while the COSMolKit bridge remains unexposed.",
        "Run the strict support-registry tests for the unsupported InChI state.",
        "Implement the complete COSMolKit `Molecule` conversion bridge using narrowed read parts and registered operation machinery where topology mutation is required.",
        "Add strict bridge tests for graph state, coordinates, stereo, properties, copy-on-write, and unchanged source molecules.",
        "Run the strict focused InChI bridge tests.",
        "Generate provenance-pinned official InChI and RDKit `2026.03.1` daily-small golden records.",
        "Add and run the exact daily-small official-engine and RDKit adapter parity test named `rdkit_inchi_parity_daily_small`.",
        "Generate provenance-pinned strict-5000 InChI parity records.",
        "Add and run the exact strict-5000 official-engine and RDKit adapter parity test named `rdkit_inchi_parity_strict_5000` in release mode.",
        "Migrate every selected RDKit `External/INCHI-API/test.cpp` case and official upstream fixture into the Rust regression suite and run it.",
        "Add and run concurrent scalar generation, parsing, diagnostics, and InChIKey tests.",
        "Change `INCHI_FEATURE` to `supported_with_rdkit_parity` only after every strict exact parity suite passes.",
        "Expose the completed value-style Rust `Molecule` APIs and add exact Rust facade tests named `inchi_rust_facade` for every API.",
        "Run `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict inchi_rust_facade`.",
        "Expose the completed scalar InChI APIs through Python with generated stubs, add exact Python parity tests, and run `.venv/bin/pytest python/tests/test_inchi.py`.",
        "Implement batch InChI APIs with stable order and indexed structured errors, then add and run `inchi_batch_contract`.",
        "Update Rust, Python, test-scope, support-matrix, provenance, and public InChI documentation.",
        "Run `cargo fmt --all`.",
        "Run `cargo check -p cosmolkit-core --features op-contracts-strict`.",
        "Run `cargo test -p cosmolkit-core --release --features op-contracts-strict`.",
        "Run `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict`.",
        "Run `.venv/bin/maturin develop --manifest-path python/Cargo.toml`, `.venv/bin/pytest`, the Sphinx HTML build, and `.venv/bin/basedpyright python/tests python/examples`.",
        "Write `dev/gap_reports/rdkit_inchi_final_parity_validation.md` with separate source-frame, focused-Rust-test, official-C-oracle, remaining-active-call-graph, unsupported-branch, and validation results; only exact field parity may be reported as exact parity.",
    ]
    for action in integration_actions:
        number = add_read(lines, number)
        number = add_step(lines, number, False, action)

    content = "\n".join(lines) + "\n"
    validate(content, prefix, active)
    print(f"configured active functions: {len(active)}")
    print(f"completed active locations retained: {sum(item.location in done for item in active)}")
    print(f"remaining active components: {len(remaining_components) + 3}")
    print(f"generated steps: {number - 1}")
    if not dry_run:
        PLAN_PATH.write_text(content, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    regenerate(args.dry_run)


if __name__ == "__main__":
    main()
