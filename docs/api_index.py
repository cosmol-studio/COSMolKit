"""Generate independent language indexes from Sphinx's collected entries."""

from types import SimpleNamespace

from sphinx.environment.adapters.indexentries import IndexEntries


def language_indexes(app):
    entries = app.env.domains.index_domain.entries
    for language in ("python", "javascript"):
        page = f"{language}/genindex"
        # Give the standard index formatter a read-only subset, before it
        # groups names and subentries. Never mutate the shared environment.
        selected = {name: values for name, values in entries.items()
                    if name.startswith(language + "/")}
        environment = SimpleNamespace(domains=SimpleNamespace(
            index_domain=SimpleNamespace(entries=selected)))
        builder = SimpleNamespace(get_relative_uri=lambda _index, target:
                                  app.builder.get_relative_uri(page, target))
        index = IndexEntries(environment).create_index(builder)
        yield page, {
            "genindexentries": index,
            "genindexcounts": [sum(1 + len(subitems) for _, (_, subitems, _) in group)
                               for _, group in index],
            "split_index": False,
        }, "genindex.html"


def setup(app):
    app.connect("html-collect-pages", language_indexes)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
