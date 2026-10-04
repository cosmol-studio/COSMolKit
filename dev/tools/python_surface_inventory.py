#!/usr/bin/env python3
"""Read-only source inventory of the registered historical PyO3 surface.

Print JSON; do not build/import either package or write generated references.
Names are declaration evidence, never evidence of behavior equivalence.
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess


def mask_literals_and_comments(source):
    out = list(source)
    raw_pattern = re.compile(r'(?:br|r)(#*)"')
    char_pattern = re.compile(r"'(?:\\(?:u\{[^}]+\}|x[0-9a-fA-F]{2}|.)|[^'\\\n])'")
    i = 0
    while i < len(source):
        start = i
        if source.startswith("//", i):
            i = source.find("\n", i)
            if i < 0:
                i = len(source)
        elif source.startswith("/*", i):
            i += 2
            level = 1
            while i < len(source) and level:
                if source.startswith("/*", i):
                    level += 1
                    i += 2
                elif source.startswith("*/", i):
                    level -= 1
                    i += 2
                else:
                    i += 1
        elif source[i] in "br" and (raw := raw_pattern.match(source, i)):
            ending = '"' + raw.group(1)
            end = source.find(ending, raw.end())
            if end < 0:
                raise ValueError("unterminated raw literal")
            i = end + len(ending)
        elif source[i] == '"':
            i += 1
            while i < len(source):
                if source[i] == "\\":
                    i += 2
                elif source[i] == '"':
                    i += 1
                    break
                else:
                    i += 1
        elif source[i] == "'" and (char := char_pattern.match(source, i)):
            i = char.end()
        else:
            i += 1
            continue
        for j in range(start, i):
            if source[j] != "\n":
                out[j] = " "
    return "".join(out)


def brace_end(masked, start):
    level = 0
    for i in range(start, len(masked)):
        if masked[i] == "{":
            level += 1
        elif masked[i] == "}":
            level -= 1
            if not level:
                return i
    raise ValueError("unbalanced body")


def line(source, offset):
    return source.count("\n", 0, offset) + 1


def surface(source, path, prefix=""):
    masked = mask_literals_and_comments(source)
    registered = set(re.findall(r"add_class\s*::\s*<\s*(\w+)\s*>", masked))
    exports = set(re.findall(r"wrap_pyfunction!\s*\(\s*(\w+)", masked))
    classes = {}
    entries = {}

    def add(owner, name, kind, offset, body=""):
        key = (owner + "." if owner else prefix) + name
        record = entries.setdefault(key, {
            "api": key, "kinds": [], "source": path,
            "lines": [], "legacy_core_reference": False,
        })
        if kind not in record["kinds"]:
            record["kinds"].append(kind)
        source_line = line(source, offset)
        if source_line not in record["lines"]:
            record["lines"].append(source_line)
        record["legacy_core_reference"] |= "cosmolkit_core" in body

    for attr in re.finditer(r"#\[pyclass\b", masked):
        match = re.search(r"\bstruct\s+(\w+)", masked[attr.end():])
        if not match:
            raise ValueError("pyclass without struct")
        offset = attr.end() + match.start()
        rust_name = match.group(1)
        if rust_name not in registered:
            continue
        attributes = source[attr.start():offset]
        named = re.search(r'\bname\s*=\s*"([^"]+)"', attributes)
        py_name = named.group(1) if named else rust_name
        classes[rust_name] = py_name
        add("", py_name, "class", offset)
        opening = masked.index("{", offset)
        closing = brace_end(masked, opening)
        fields = source[opening + 1:closing]
        if "get_all" in attributes or "set_all" in attributes:
            for field in re.finditer(r"(?m)^\s*(?:pub(?:\([^)]*\))?\s+)?(\w+)\s*:", mask_literals_and_comments(fields)):
                add(py_name, field.group(1), "property", opening + 1 + field.start())

    for attr in re.finditer(r"#\[pymethods\]", masked):
        impl = re.search(r"\bimpl(?:<[^\n]*?>)?\s+(\w+)", masked[attr.end():])
        if not impl:
            raise ValueError("pymethods without impl")
        offset = attr.end() + impl.start()
        owner = classes.get(impl.group(1))
        if owner is None:
            continue
        opening = masked.index("{", offset)
        closing = brace_end(masked, opening)
        cursor = opening + 1
        while cursor < closing:
            fn = re.search(r"\bfn\s+(\w+)\b", masked[cursor:closing])
            if not fn:
                break
            start = cursor + fn.start()
            header = source[cursor:start]
            body_open = masked.index("{", start)
            body_close = brace_end(masked, body_open)
            name = fn.group(1)
            renamed = re.search(r'\bname\s*=\s*"([^"]+)"', header)
            if renamed:
                name = renamed.group(1)
            if "#[new]" in header:
                name, kind = "__new__", "constructor"
            elif (prop := re.search(r"#\[(getter|setter)(?:\((\w+)\))?\]", header)):
                name = prop.group(2) or name
                if prop.group(1) == "setter" and not prop.group(2) and name.startswith("set_"):
                    name = name[4:]
                kind = "property"
            elif "#[classattr]" in header:
                kind = "class_attribute"
            elif name.startswith("__"):
                kind = "protocol"
            else:
                kind = "method"
            add(owner, name, kind, start, source[body_open:body_close + 1])
            cursor = body_close + 1
        block = source[opening + 1:closing]
        for const in re.finditer(r"#\[classattr\]\s*(?:pub\s+)?const\s+(\w+)", block):
            add(owner, const.group(1), "class_attribute", opening + 1 + const.start())

    for attr in re.finditer(r"#\[pyfunction\b", masked):
        fn = re.search(r"\bfn\s+(\$?\w+)\b", masked[attr.end():])
        start = attr.end() + fn.start()
        name = fn.group(1)
        if name.startswith("$"):
            continue
        if name not in exports:
            continue
        body_open = masked.index("{", start)
        body_close = brace_end(masked, body_open)
        renamed = re.search(r'\bname\s*=\s*"([^"]+)"', source[attr.start():start])
        add("", renamed.group(1) if renamed else name, "module_function", start, source[body_open:body_close + 1])
    # These four complete, inspected macros emit one PyO3 function per call.
    macro_names = (
        "python_infallible_float_descriptor", "python_float_descriptor",
        "python_count_descriptor", "python_fixed_chi_descriptor",
    )
    for macro in macro_names:
        for call in re.finditer(r"\b" + macro + r"!\s*\(\s*(\w+)", masked):
            name = call.group(1)
            if name in exports:
                add("", name, "module_function", call.start(), "cosmolkit_core")
    found_exports = {key.removeprefix(prefix) for key, entry in entries.items()
                     if "module_function" in entry["kinds"]}
    if found_exports != exports:
        raise ValueError("unresolved registered functions: " + repr(sorted(exports - found_exports)))
    return entries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--root", default=".")
    parser.add_argument("--summary", action="store_true")
    parser.add_argument("--tsv", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    items = {}
    sources = {}
    for path, prefix in [("python/src/lib.rs", ""), ("python/src/confseq_py.rs", "confseq.")]:
        source = subprocess.check_output(["git", "show", f"{args.baseline}:{path}"], cwd=root).decode()
        sources[path] = hashlib.sha256(source.encode()).hexdigest()
        current_path = root / path
        sources["current:" + path] = hashlib.sha256(current_path.read_bytes()).hexdigest()
        parsed = surface(source, path, prefix)
        for key in parsed:
            if key in items:
                raise ValueError("duplicate API " + key)
        items.update(parsed)
    if args.tsv:
        print("baseline_python_api\tkinds\tbaseline_source\tbaseline_lines")
        for entry in sorted(items.values(), key=lambda x: x["api"]):
            print("\t".join((entry["api"], ",".join(entry["kinds"]),
                             entry["source"], ",".join(map(str, entry["lines"])))))
        return
    counts = Counter(kind for item in items.values() for kind in item["kinds"])
    families = {}
    for item in items.values():
        owner = item["api"].split(".")[0] if "." in item["api"] else "module"
        families.setdefault(owner, Counter()).update(item["kinds"])
    result = {
        "baseline_commit": args.baseline, "current_root": str(root),
        "source_sha256": sources, "counts": dict(sorted(counts.items())),
        "total_unique_declaration_slots": len(items),
        "families": {k: dict(v) for k, v in sorted(families.items())},
    }
    if not args.summary:
        result["entries"] = sorted(items.values(), key=lambda x: x["api"])
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
