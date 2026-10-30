"""Generate a public package from official PyO3 output without changing calls.

Inputs: an already compiled Rust extension, its official extracted stub, and an
output directory. This script does not compile Rust, run an extractor, build a
wheel, consult the registry, or execute runtime/type tests.

Changes to the public __init__.pyi:
1. Parse the official stub into a Python syntax tree (AST) and remove the
   cosmolkit / builtins qualifiers from type names.
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
   metadata prefix from the public stub, retaining the human-readable description.
   This is a handwritten annotation, not Rust type inference; there are no
   NumPy-specific type names or conversion rules in this generator.
4. Expand explicit integer_selection indexing declarations with a named selection
   result (or Self), validating integer/slice/indices/mask inputs and concrete
   item/selection outputs. Unmarked __getitem__ declarations remain unchanged;
   incompatible marked signatures fail generation. Preserve official types.
5. Remove positional-only separators from property setters as a BasedPyright
   workaround. Qualify property decorators as builtins.property so an actual
   member named property cannot shadow the decorator. Add imports and derive
   __all__ from the extension's public exports by loading the extension.
   Re-export explicit {"role": "private_export"} functions without adding them
   to public __all__; attribute them to the public package for stable pickle
   lookup. Their declarations remain officially extracted.
   For explicit error_projection declarations, use the official payload
   signatures as the public ValueError subclass. Remove the backend factory;
   turn explicitly marked exception-receiver getters into readonly properties.

Executable __init__.py:
6. Re-export the existing native classes/functions without wrappers, replacement
   methods or new classes. All argument processing remains in the Rust binding.
   Indexing and property-setter adjustments above affect the stub only.
   Emit the historical cosmolkit.cosmolkit reconstruction module for 0.5.0
   pickles; it re-exports only the existing binary reconstruction function.
   Emit __version__ from pyproject.toml in both package and stub.

Files: write __init__.pyi and __init__.py, copy the compiled extension under its
importable _native name, and create py.typed. Do not modify the official stub or
copy a second _native.pyi into the package. AST serialization rewrites the whole
public stub's formatting rather than preserving its original source text.

Indexing semantics belong to the binding's explicit protocol declaration, not
the method name or an inferred relationship between arbitrary union members.
"""
from __future__ import annotations

import ast
import tomllib
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
        if name.startswith("cosmolkit.") or name.startswith("builtins."):
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
    return json.JSONDecoder().raw_decode(doc)[0] if doc.startswith("{") else {}


def description(function):
    """Projection metadata is a JSON prefix; the remaining text is user help."""
    doc = (ast.get_docstring(function) or "").strip()
    if doc.startswith("{"):
        _, end = json.JSONDecoder().raw_decode(doc)
        return doc[end:].strip()
    return doc


def documented_body(text):
    return ([ast.Expr(value=ast.Constant(value=text))] if text else []) + [
        ast.Expr(value=ast.Constant(value=Ellipsis))
    ]


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
        function.body = documented_body(description(function))
    tree.body[:0] = [ast.Import(names=[ast.alias(name=module)]) for module in sorted(imports)]


def declared_errors(tree):
    """Project one explicitly declared backend; never infer errors by names."""
    for cls in tree.body:
        if not isinstance(cls, ast.ClassDef):
            continue
        markers = [node for node in cls.body if isinstance(node, ast.FunctionDef)
                   and projection(node).get('role') == 'error_projection']
        if not markers:
            continue
        if len(markers) != 1 or markers[0].name != '_exception_type':
            raise ValueError(f'{cls.name}: invalid error backend declaration')
        cls.body.remove(markers[0])
        cls.bases = [ast.Name(id='ValueError', ctx=ast.Load())]
        cls.decorator_list = [d for d in cls.decorator_list if ast.unparse(d) not in {'final', 'typing.final'}]
        for method in cls.body:
            if not isinstance(method, ast.FunctionDef):
                continue
            if any(ast.unparse(d) == 'staticmethod' for d in method.decorator_list):
                if projection(method).get('role') not in {'error_source', 'error_payload'}:
                    raise ValueError(f'{cls.name}.{method.name}: undeclared exception receiver')
                args = positional(method)
                if len(args) != 1 or method.args.kwonlyargs or method.args.vararg or method.args.kwarg:
                    raise ValueError(f'{cls.name}.{method.name}: expected one exception receiver')
                method.args = ast.arguments(posonlyargs=[ast.arg(arg='self')], args=[], kwonlyargs=[], kw_defaults=[], defaults=[])
                method.decorator_list = [ast.Name(id='property', ctx=ast.Load())]


def union_parts(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr):
        return union_parts(node.left) + union_parts(node.right)
    return [node]


def indexing_forms(owner, method):
    """Validate explicit integer/item and multi-selection/result correspondence."""
    metadata = projection(method)
    context = f"{owner}.{method.name}: integer_selection"
    if (set(metadata) != {"role", "protocol", "selection"}
            or metadata["role"] != "indexing" or metadata["protocol"] != "integer_selection"
            or not isinstance(metadata["selection"], str)
            or not metadata["selection"].isidentifier()
            or metadata["selection"] in {"Any", "PyAny", "Incomplete", "object", "None"}):
        raise ValueError(f"{context} requires the explicit supported indexing declaration")
    selected = owner if metadata["selection"] == "Self" else metadata["selection"]
    args = method.args
    parameters = positional(method)
    if (method.name != "__getitem__" or len(parameters) != 2
            or parameters[0].arg != "self" or args.vararg or args.kwarg
            or args.kwonlyargs or args.defaults):
        raise ValueError(f"{context} requires exactly (self, index), without defaults or extra parameters")
    if parameters[1].annotation is None or method.returns is None:
        raise ValueError(f"{context} requires explicit input and return annotations")
    inputs = union_parts(parameters[1].annotation)
    names = {ast.unparse(node) for node in inputs}
    integer_names = {"int", "SupportsIndex", "typing.SupportsIndex"}
    sequences = []
    sequence_items = set()
    for node in inputs:
        if not isinstance(node, ast.Subscript):
            continue
        container = ast.unparse(node.value)
        if container in {"Sequence", "typing.Sequence", "collections.abc.Sequence"}:
            item = node.slice
        elif (container in {"SupportsGetItem", "_typeshed.SupportsGetItem"}
              and isinstance(node.slice, ast.Tuple) and len(node.slice.elts) == 2
              and ast.unparse(node.slice.elts[0]) == "int"):
            # Official PyO3 Vec<T> input: an integer-indexed Python sequence.
            item = node.slice.elts[1]
        else:
            continue
        sequences.append(node)
        sequence_items.add(ast.unparse(item))
    if (len(inputs) != 4 or "slice" not in names or len(names & integer_names) != 1
            or len(sequences) != 2 or "bool" not in sequence_items
            or len(sequence_items & integer_names) != 1):
        raise ValueError(f"{context} requires integer, slice, integer-sequence and boolean-mask inputs")
    integer = next(node for node in inputs if ast.unparse(node) in integer_names)
    outputs = union_parts(method.returns)
    items = [node for node in outputs if ast.unparse(node) != selected]
    if len(items) == len(outputs) or not items or len(outputs) - len(items) != 1:
        raise ValueError(f"{context} return must contain an item type and exactly one {selected} branch")
    if any(isinstance(node, ast.Name) and node.id in {"Any", "PyAny", "Incomplete", "object"}
           or isinstance(node, ast.Attribute) and node.attr in {"Any", "PyAny", "Incomplete"}
           for item in items for node in ast.walk(item)):
        raise ValueError(f"{context} requires concrete item types")
    item_type = items[0]
    for node in items[1:]:
        item_type = ast.BinOp(left=item_type, op=ast.BitOr(), right=node)
    forms = []
    selection_input = ast.Name(id="slice", ctx=ast.Load())
    for sequence in sequences:
        selection_input = ast.BinOp(left=selection_input, op=ast.BitOr(), right=copy.deepcopy(sequence))
    for input_type, output_type in ((integer, item_type), (selection_input, ast.Name(id=selected, ctx=ast.Load()))):
        form = overload(method)
        positional(form)[-1].annotation = copy.deepcopy(input_type)
        form.returns = copy.deepcopy(output_type)
        form.body = documented_body(description(method))
        forms.append(form)
    return forms


def declared_indexing(tree):
    """No declaration means no indexing rewrite, including nested classes."""
    owners = {method: cls for cls in ast.walk(tree) if isinstance(cls, ast.ClassDef)
              for method in cls.body if isinstance(method, ast.FunctionDef)}
    for method in list(ast.walk(tree)):
        if not isinstance(method, ast.FunctionDef) or projection(method).get("role") != "indexing":
            continue
        cls = owners.get(method)
        if cls is None:
            raise ValueError("indexing declarations require a class __getitem__ method")
        at = cls.body.index(method)
        cls.body[at:at + 1] = indexing_forms(cls.name, method)


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
            form.body = documented_body(description(function))
            groups.setdefault(target, []).append(form)
            body.remove(function)
    for target, forms in groups.items():
        owner, _, name = target.rpartition(".")
        body = bodies[owner]
        at = next(index for index, node in enumerate(body)
                  if isinstance(node, ast.FunctionDef) and node.name == name)
        for form in forms:
            if not description(form):
                form.body = documented_body(description(body[at]))
        body[at:at + 1] = forms


def generate(library: Path, official_stub: Path, output: Path):
    loader = ExtensionFileLoader("_native", str(library))
    spec = importlib.util.spec_from_file_location("_native", library, loader=loader)
    native = importlib.util.module_from_spec(spec)
    loader.exec_module(native)
    tree = PublicNames().visit(ast.parse(official_stub.read_text()))
    declared_errors(tree)
    declared_overloads(tree)
    declared_return_types(tree)
    declared_indexing(tree)
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    for cls in classes.values():
        for method in cls.body:
            if isinstance(method, ast.FunctionDef):
                method.decorator_list = [
                    ast.Attribute(value=ast.Name(id="builtins", ctx=ast.Load()),
                                  attr="property", ctx=ast.Load())
                    if isinstance(tag, ast.Name) and tag.id == "property" else tag
                    for tag in method.decorator_list
                ]
            if isinstance(method, ast.FunctionDef) and any(isinstance(tag, ast.Attribute) and tag.attr == "setter" for tag in method.decorator_list):
                # Work around BasedPyright's positional-only property-setter handling.
                method.args.args = method.args.posonlyargs + method.args.args
                method.args.posonlyargs = []
    exports = sorted(name for name in native.__all__ if not name.startswith("_"))
    # Module initialization installs the explicitly declared exception backend.
    # The extractor conservatively adds an open-ended __getattr__; this package
    # instead has a closed, fully declared native export surface.
    declared = {node.name for node in tree.body if isinstance(node, (ast.ClassDef, ast.FunctionDef))}
    if "__getattr__" in declared and "__getattr__" not in vars(native):
        if set(exports) - declared:
            raise ValueError(f"Undeclared native exports: {sorted(set(exports) - declared)}")
        tree.body = [node for node in tree.body
                     if not isinstance(node, ast.FunctionDef) or node.name != "__getattr__"]
    private_exports = sorted(node.name for node in tree.body
                             if isinstance(node, ast.FunctionDef)
                             and projection(node).get("role") == "private_export")
    if any(not name.startswith("_") or name not in native.__all__ for name in private_exports):
        raise ValueError("private_export requires an explicitly exported native private function")
    tree.body = [node for node in tree.body if not isinstance(node, ast.Assign) or not any(isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets)]
    tree.body.insert(0, ast.Import(names=[ast.alias(name="builtins"), ast.alias(name="typing")]))
    version = tomllib.loads((Path(__file__).parent / "pyproject.toml").read_text())["project"]["version"]
    tree.body.append(ast.parse(f"__version__: str = {version!r}").body[0])
    tree.body.append(ast.parse(f"__all__ = {exports!r}").body[0])

    package = output / "cosmolkit"
    package.mkdir(parents=True, exist_ok=True)
    shutil.copy2(library, package / ("_native" + EXTENSION_SUFFIXES[-1]))
    (package / "__init__.pyi").write_text(ast.unparse(ast.fix_missing_locations(tree)) + "\n")
    runtime = f"from . import _native\nfrom ._native import {', '.join(exports + private_exports)}\n"
    runtime += "".join(f"{name}.__module__ = __name__\n" for name in private_exports)
    runtime += f"__version__ = {version!r}\n__all__ = {exports!r}\n"
    (package / "__init__.py").write_text(runtime)
    (package / "cosmolkit.py").write_text(
        '"""Private reconstruction path retained for 0.5.0 pickle files."""\n'
        'from . import _molecule_from_binary\n'
    )
    (package / "py.typed").touch()
    return package


if __name__ == "__main__":
    generate(*(Path(argument).resolve() for argument in sys.argv[1:]))
