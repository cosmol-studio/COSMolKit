// Independent of Dioxus hydration. Only the catalog is shared with latest;
// archive navigation and assets stay on the archive's own origin.
const catalogUrl = "https://kit.cosmol.org/versions.json";

function validateCatalog(catalog) {
    if (catalog?.schema_version !== 1 || !Array.isArray(catalog.versions)
        || catalog.versions.length === 0) throw new Error("Invalid docs version catalog");
    const names = new Set();
    const urls = new Set();
    for (const entry of catalog.versions) {
        if (typeof entry.version !== "string" || !entry.version.trim()
            || typeof entry.url !== "string") throw new Error("Invalid docs version entry");
        const url = new URL(entry.url);
        if (url.protocol !== "https:" || url.username || url.password
            || url.pathname !== "/" || url.search || url.hash
            || names.has(entry.version) || urls.has(url.href)) {
            throw new Error("Invalid or duplicate docs version URL");
        }
        names.add(entry.version);
        urls.add(url.href);
    }
    if (catalog.versions[0].version !== "latest"
        || catalog.versions[0].url !== "https://kit.cosmol.org/") {
        throw new Error("Docs catalog must start with the canonical latest site");
    }
    return catalog.versions;
}

function renderVersions(versions) {
    const options = document.getElementById("docs-version-options");
    const current = document.getElementById("docs-version-current");
    if (!options || !current) return;
    const selected = versions.find(entry => new URL(entry.url).origin === location.origin);
    // Development and unnamed preview deployments show latest, never redirect.
    // Keep a recognized built-in archive label if its snapshot URL was later
    // replaced in the shared catalog. An old deployment is still that version.
    current.textContent = selected?.version ?? current.textContent ?? "latest";
    options.replaceChildren(...versions.map(entry => {
        const link = document.createElement("a");
        link.textContent = entry.version;
        link.dataset.docsVersion = entry.version;
        // Version overviews are safe even when a topic does not exist there.
        link.setAttribute("href", entry === selected ? "/" : entry.url);
        if (entry === selected) link.setAttribute("aria-current", "true");
        return link;
    }));
}

async function refreshVersions() {
    const options = document.getElementById("docs-version-options");
    if (!options) return;
    renderVersions([...options.querySelectorAll("a[data-docs-version]")].map(link => ({
        version: link.dataset.docsVersion, url: link.href,
    })));
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 5000);
    try {
        const response = await fetch(catalogUrl, {
            mode: "cors", credentials: "omit", cache: "no-store", signal: controller.signal,
        });
        if (!response.ok) throw new Error("Docs catalog unavailable");
        renderVersions(validateCatalog(await response.json()));
    } catch (_) {
        // Retain the built-in links; do not redirect or load latest page assets.
    } finally {
        clearTimeout(timer);
    }
}

refreshVersions();
