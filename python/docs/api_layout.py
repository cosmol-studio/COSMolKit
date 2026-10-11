"""Present generated API objects with inline types and preserved anchors."""

from docutils import nodes
from sphinx import addnodes


API_PAGES = {"api", "javascript-api"}
CLASS_KINDS = {"class", "exception", "interface"}
MEMBER_ANCHOR_SCRIPT = """<script>
(() => {
    const reveal = () => {
        const target = document.querySelector(':target');
        const members = target?.closest('details.api-members');
        if (members && !members.open) {
            members.open = true;
            target.scrollIntoView();
        }
    };
    addEventListener('pageshow', reveal);
    addEventListener('hashchange', reveal);
})();
</script>"""


def _name(entry):
    signature = next(child for child in entry if isinstance(child, addnodes.desc_signature))
    return (signature.get("fullname") or signature["ids"][0]).rsplit(".", 1)[-1]


def _entries(doctree):
    return [entry for section in doctree.findall(nodes.section)
            if isinstance(section.parent, nodes.document)
            for entry in section if isinstance(entry, addnodes.desc)]


def format_javascript_signature(app, domain, objtype, content):
    """Move TypeDoc's types into signatures, retaining authored field prose."""
    if domain != "js" or objtype not in {"function", "class"}:
        return
    signature = next(child for child in content.parent
                     if isinstance(child, addnodes.desc_signature))
    parameters = {
        node.astext().split("=", 1)[0].removeprefix("...").removesuffix("?"): node
        for node in signature.findall(addnodes.desc_parameter)
    }
    # This event precedes Sphinx's field-list conversion, so field names still
    # contain their semantic roles ("type query", "param query", "returns").
    for fields in list(content):
        if not isinstance(fields, nodes.field_list):
            continue
        for field in list(fields):
            name, body = field
            role, _, argument = name.astext().partition(" ")
            if role in {"type", "paramtype"} and argument in parameters:
                parameter = parameters[argument]
                parameter += addnodes.desc_sig_punctuation("", ":")
                parameter += addnodes.desc_sig_space()
                for paragraph in body:
                    parameter.extend(child.deepcopy() for child in paragraph)
                fields.remove(field)
            elif role in {"param", "arg", "argument", "parameter"}:
                if not body.astext().strip() and argument in parameters:
                    fields.remove(field)
            elif role in {"returns", "return"} and objtype == "function":
                if len(body) == 0 or not isinstance(body[0], nodes.paragraph) or len(body[0]) == 0:
                    continue
                paragraph = body[0]
                result_type = paragraph[0]
                if not (isinstance(result_type, nodes.strong) or
                        isinstance(result_type, nodes.inline) and
                        "sphinx_js-type" in result_type.get("classes", [])):
                    continue
                signature += addnodes.desc_returns(
                    "", "", *(child.deepcopy() for child in result_type),
                )
                paragraph.remove(result_type)
                # sphinx-js joins a return type and its description with " -- ".
                if len(paragraph) and isinstance(paragraph[0], nodes.Text):
                    paragraph[0] = nodes.Text(str(paragraph[0]).removeprefix(" -- "))
                if not body.astext().strip():
                    fields.remove(field)
        if len(fields) == 0:
            content.remove(fields)


def format_api(app, doctree, docname):
    if docname not in API_PAGES:
        return

    priority = {name: index for index, name in enumerate(app.config.api_class_order)}

    def order(entry):
        kind, name = entry.get("objtype"), _name(entry)
        if kind == "function":
            return (0, 0, name.casefold())
        if kind in CLASS_KINDS:
            return (1, priority.get(name, len(priority)), name.casefold())
        return (2, 0, name.casefold())

    for section in doctree.findall(nodes.section):
        if not isinstance(section.parent, nodes.document):
            continue
        positions = [i for i, child in enumerate(section) if isinstance(child, addnodes.desc)]
        ordered = sorted((section[i] for i in positions), key=order)
        for index, entry in zip(positions, ordered):
            section[index] = entry

    for entry in list(doctree.findall(addnodes.desc)):
        if entry.get("objtype") not in CLASS_KINDS:
            continue
        content = next(child for child in entry if isinstance(child, addnodes.desc_content))
        members = [i for i, child in enumerate(content) if isinstance(child, addnodes.desc)]
        if not members:
            continue
        start, end = members[0], members[-1] + 1
        member_list = nodes.container(classes=["api-member-list"])
        member_list.extend(content.children[start:end])
        content[start:end] = [
            nodes.raw("", f'<details class="api-members"><summary>Members ({len(members)})</summary>', format="html"),
            member_list,
            nodes.raw("", "</details>", format="html"),
        ]
    # Chromium restores closed details on reload, even for a member fragment.
    # Keep normal pages collapsed and reveal only the explicitly linked member.
    if not any(node.get("api_member_anchor") for node in doctree.findall(nodes.raw)):
        doctree += nodes.raw("", MEMBER_ANCHOR_SCRIPT, format="html", api_member_anchor=True)


def api_toc(app, pagename, templatename, context, doctree):
    if pagename not in API_PAGES:
        return
    toc = nodes.bullet_list()
    for entry in _entries(doctree):
        kind = entry.get("objtype")
        if kind not in CLASS_KINDS | {"function"}:
            continue
        signature = next(child for child in entry if isinstance(child, addnodes.desc_signature))
        if not signature["ids"]:
            continue
        label = _name(entry) + ("()" if kind == "function" else "")
        link = nodes.reference("", "", internal=True, refuri="#" + signature["ids"][0])
        link += nodes.literal("", label)
        toc += nodes.list_item("", addnodes.compact_paragraph("", "", link))
    context["toc"] = app.builder.render_partial(toc)["fragment"]


def setup(app):
    app.add_config_value("api_class_order", [], "html")
    app.connect("object-description-transform", format_javascript_signature)
    app.connect("doctree-resolved", format_api)
    app.connect("html-page-context", api_toc)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
