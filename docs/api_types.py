"""Supply native Python API types from the generated, checked stub."""

import ast
import json
from pathlib import Path

from sphinx import addnodes
from sphinx.domains.python import _parse_annotation


def read_stub(path):
    """Index getters/functions and resolve thin module-level factory aliases."""
    entries = {}
    aliases = {}

    def visit(body, prefix=""):
        for node in body:
            if isinstance(node, ast.ClassDef):
                visit(node.body, prefix + node.name + ".")
            elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                # A property setter returns None, not the property's value type.
                if any(isinstance(d, ast.Attribute) and d.attr == "setter"
                       for d in node.decorator_list):
                    continue
                if node.returns is not None:
                    entries.setdefault(prefix + node.name, []).append(node)
            elif isinstance(node, ast.Assign) and not prefix:
                if isinstance(node.value, (ast.Name, ast.Attribute)):
                    for target in node.targets:
                        if isinstance(target, ast.Name):
                            aliases[target.id] = ast.unparse(node.value)

    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    imported = {
        alias.asname or alias.name: f"{node.module}.{alias.name}"
        for node in tree.body
        if isinstance(node, ast.ImportFrom) and node.module and not node.level
        for alias in node.names
    }

    class QualifyTypes(ast.NodeTransformer):
        def visit_Name(self, node):
            if node.id in imported:
                return ast.copy_location(ast.parse(imported[node.id], mode="eval").body, node)
            return node

    # Resolve annotation imports, not enum members with the same short name.
    qualify = QualifyTypes()
    for node in ast.walk(tree):
        if isinstance(node, ast.arg) and node.annotation is not None:
            node.annotation = qualify.visit(node.annotation)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.returns is not None:
            node.returns = qualify.visit(node.returns)
    visit(tree.body)
    for alias, target in aliases.items():
        if target in entries:
            entries[alias] = entries[target]
    return entries


def annotations(app, name):
    if not name.startswith("cosmolkit."):
        return []
    return app._python_api_types.get(name.removeprefix("cosmolkit."), [])


def type_text(node):
    # Builtin names need no qualification; keep library/type namespaces intact.
    return ast.unparse(node).replace("builtins.", "")


def return_signature(app, what, name, obj, options, signature, return_annotation):
    if what not in {"function", "method"}:
        return None
    overloads = annotations(app, name)
    if not overloads:
        return None
    returns = {type_text(node.returns) for node in overloads}
    if len(returns) == 1:
        return signature if signature is not None else "()", returns.pop()

    # Different overload results (e.g. indexing versus slicing) must remain
    # paired with their own arguments, rather than becoming an ambiguous union.
    signatures = []
    for node in overloads:
        args = ast.unparse(node.args)
        if args.startswith("self, ") or args.startswith("cls, "):
            args = args.split(", ", 1)[1]
        elif args in {"self", "cls"}:
            args = ""
        signatures.append(f"({args}) -> {type_text(node.returns)}")
    return "\n".join(signatures), None


def property_type(app, domain, objtype, content):
    if domain != "py" or objtype not in {"attribute", "property"}:
        return
    for signature in content.parent.findall(addnodes.desc_signature):
        name = ".".join(filter(None, (signature.get("module"), signature.get("fullname"))))
        entries = annotations(app, name)
        already_typed = any(
            isinstance(child, addnodes.desc_annotation) and child.astext().lstrip().startswith(":")
            for child in signature
        )
        if entries and not already_typed:
            annotation = type_text(entries[0].returns)
            signature += addnodes.desc_annotation(
                annotation, "", addnodes.desc_sig_punctuation("", ":"),
                addnodes.desc_sig_space(), *_parse_annotation(annotation, app.env),
            )


def load_types(app):
    app._python_api_types = read_stub(app.config.python_api_stub)


def note_stub(app, docname, source):
    app.env.note_dependency(app.config.python_api_stub)


def public_description(app, what, name, obj, options, lines):
    """Omit explicit binding projection metadata, retaining the native prose."""
    if not name.startswith("cosmolkit."):
        return
    text = "\n".join(lines).lstrip()
    if not text.startswith("{"):
        return
    try:
        metadata, end = json.JSONDecoder().raw_decode(text)
    except ValueError:
        return
    if isinstance(metadata, dict) and "role" in metadata:
        lines[:] = text[end:].strip().splitlines()


def setup(app):
    app.add_config_value("python_api_stub", "", "env")
    app.connect("builder-inited", load_types)
    app.connect("source-read", note_stub)
    app.connect("autodoc-process-signature", return_signature)
    app.connect("autodoc-process-docstring", public_description)
    app.connect("object-description-transform", property_type)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
