"""Generate a public package from official PyO3 output without changing calls.

Inputs: an already compiled Rust extension, its official extracted stub, and an
output directory. This script does not compile Rust, run an extractor, build a
wheel, consult the registry, or execute runtime/type tests.

Changes to the public __init__.pyi:
1. Parse the official stub into a Python syntax tree (AST) and remove the
   fixed_projection._native / builtins qualifiers from type names.
2. Read JSON {"role": "overload", "target": "Class.method"} documentation on
   registered native branch functions (module targets omit Class). Group their
   official signatures by target, rename them and add @typing.overload, replacing
   the broad entrypoint declaration. Hide branch names from the public stub.
   These branches are real functions called by the native dispatcher, not dummy
   signature carriers. Do not infer options from Params fields or constructors;
   the binding owns initialization, call forms, defaults and exclusivity.
3. Apply explicit {"role": "return_type", "imports": [...], "type": "..."}
   binding declarations: import the named modules, parse the Python type
   expression and replace only that callable's return annotation. Remove its
   metadata docstring from the public stub, not from native runtime __doc__.
   This is a handwritten annotation, not Rust type inference; there are no
   NumPy-specific type names or conversion rules in this generator.
4. Split __getitem__ input/output unions into overloads: integer-like indexing
   returns the item type (including None if declared); slicing returns the same
   collection type. Preserve the official integer input protocol.
5. Remove positional-only separators from property setters as a BasedPyright
   workaround. Add typing imports and derive __all__ from the extension's public
   exports by loading the extension.

Executable __init__.py:
6. Re-export the existing native classes/functions without wrappers, replacement
   methods or new classes. All argument processing remains in the Rust binding.
   Indexing and property-setter adjustments above affect the stub only.

Files: write __init__.pyi and __init__.py, copy the compiled extension under its
importable _native name, and create py.typed. Do not modify the official stub or
copy a second _native.pyi into the package. AST serialization rewrites the whole
public stub's formatting rather than preserving its original source text.

Indexing convention: integer-like | slice input and item | collection output.
This is a contributor convention, not a claim to handle arbitrary indexing.
"""
from __future__ import annotations

import ast
import copy
from importlib.machinery import EXTENSION_SUFFIXES, ExtensionFileLoader
import importlib.util
import json
from pathlib import Path
import shutil
import sys


class PublicNames(ast.NodeTransformer):
    def visit_Attribute(self, node):
        name = ast.unparse(node)
        if name.startswith("fixed_projection._native.") or name.startswith("builtins."):
            return ast.copy_location(ast.Name(id=node.attr, ctx=ast.Load()), node)
        return self.generic_visit(node)


def positional(function):
    return function.args.posonlyargs + function.args.args


def overload(function):
    result = copy.deepcopy(function)
    result.decorator_list.insert(0, ast.parse("typing.overload", mode="eval").body)
    return result


def projection(function):
    doc = (ast.get_docstring(function) or "").strip()
    return json.loads(doc) if doc.startswith("{") else {}


def declared_return_types(tree):
    """Apply explicit type expressions without knowing their library or API."""
    imports = set()
    for function in ast.walk(tree):
        if not isinstance(function, ast.FunctionDef):
            continue
        metadata = projection(function)
        if metadata.get("role") != "return_type":
            continue
        function.returns = ast.parse(metadata["type"], mode="eval").body
        imports.update(metadata.get("imports", []))
        function.body = function.body[1:] or [ast.Expr(value=ast.Constant(value=Ellipsis))]
    tree.body[:0] = [ast.Import(names=[ast.alias(name=module)]) for module in sorted(imports)]


def union_parts(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr):
        return union_parts(node.left) + union_parts(node.right)
    return [node]


def indexing_forms(owner, method):
    inputs = union_parts(positional(method)[-1].annotation)
    integer = next(node for node in inputs if ast.unparse(node) != "slice")
    items = [node for node in union_parts(method.returns) if ast.unparse(node) != owner]
    item_type = items[0]
    for node in items[1:]:
        item_type = ast.BinOp(left=item_type, op=ast.BitOr(), right=node)
    forms = []
    for input_type, output_type in ((integer, item_type), (ast.Name(id="slice", ctx=ast.Load()), ast.Name(id=owner, ctx=ast.Load()))):
        form = overload(method)
        positional(form)[-1].annotation = copy.deepcopy(input_type)
        form.returns = copy.deepcopy(output_type)
        forms.append(form)
    return forms


def declared_overloads(tree):
    """Group officially extracted signatures using explicit binding metadata."""
    bodies = {"": tree.body}
    bodies.update({node.name: node.body for node in tree.body if isinstance(node, ast.ClassDef)})
    groups = {}
    for owner, body in bodies.items():
        for function in list(body):
            if not isinstance(function, ast.FunctionDef):
                continue
            metadata = projection(function)
            if metadata.get("role") != "overload":
                continue
            target = metadata["target"]
            form = overload(function)
            form.name = target.rsplit(".", 1)[-1]
            form.body = [ast.Expr(value=ast.Constant(value=Ellipsis))]
            groups.setdefault(target, []).append(form)
            body.remove(function)
    for target, forms in groups.items():
        owner, _, name = target.rpartition(".")
        body = bodies[owner]
        at = next(index for index, node in enumerate(body)
                  if isinstance(node, ast.FunctionDef) and node.name == name)
        body[at:at + 1] = forms


def generate(library: Path, official_stub: Path, output: Path):
    loader = ExtensionFileLoader("_native", str(library))
    spec = importlib.util.spec_from_file_location("_native", library, loader=loader)
    native = importlib.util.module_from_spec(spec)
    loader.exec_module(native)
    tree = PublicNames().visit(ast.parse(official_stub.read_text()))
    declared_overloads(tree)
    declared_return_types(tree)
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    for owner, cls in classes.items():
        body = cls.body
        for method in list(body):
            if isinstance(method, ast.FunctionDef) and method.name == "__getitem__":
                at = body.index(method)
                body[at:at + 1] = indexing_forms(owner, method)
    for cls in classes.values():
        for method in cls.body:
            if isinstance(method, ast.FunctionDef) and any(isinstance(tag, ast.Attribute) and tag.attr == "setter" for tag in method.decorator_list):
                # Work around BasedPyright's positional-only property-setter handling.
                method.args.args = method.args.posonlyargs + method.args.args
                method.args.posonlyargs = []
    exports = sorted(name for name in native.__all__ if not name.startswith("_"))
    tree.body = [node for node in tree.body if not isinstance(node, ast.Assign) or not any(isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets)]
    tree.body.insert(0, ast.Import(names=[ast.alias(name="typing")]))
    tree.body.append(ast.parse(f"__all__ = {exports!r}").body[0])

    package = output / "fixed_projection"
    package.mkdir(parents=True, exist_ok=True)
    shutil.copy2(library, package / ("_native" + EXTENSION_SUFFIXES[-1]))
    (package / "__init__.pyi").write_text(ast.unparse(ast.fix_missing_locations(tree)) + "\n")
    runtime = f"from . import _native\nfrom ._native import {', '.join(exports)}\n__all__ = {exports!r}\n"
    (package / "__init__.py").write_text(runtime)
    (package / "py.typed").touch()
    return package


if __name__ == "__main__":
    generate(*(Path(argument).resolve() for argument in sys.argv[1:]))
