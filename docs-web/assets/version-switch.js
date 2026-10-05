// HTML/CSS owns the menu; only the version catalog comes from latest.
const catalogUrl = "https://kit.cosmol.org/versions.json";
const options = document.getElementById("docs-version-options");
const current = document.getElementById("docs-version-current");

function renderVersions(versions) {
    const selected = versions.find(entry => new URL(entry.url).origin === location.origin);
    current.textContent = selected?.version ?? current.textContent;
    options.replaceChildren(...versions.map(entry => {
        const link = document.createElement("a");
        link.textContent = entry.version;
        link.dataset.docsVersion = entry.version;
        link.setAttribute("href", entry === selected ? "/" : entry.url);
        if (entry === selected) link.setAttribute("aria-current", "true");
        return link;
    }));
}

async function refreshVersions() {
    if (!options || !current) return;
    renderVersions([...options.querySelectorAll("a[data-docs-version]")].map(link => ({
        version: link.dataset.docsVersion, url: link.href,
    })));
    try {
        const response = await fetch(catalogUrl, {mode: "cors", credentials: "omit", cache: "no-store"});
        if (!response.ok) return;
        const catalog = await response.json(), versions = catalog?.versions;
        if (catalog?.schema_version !== 1 || !Array.isArray(versions)
            || versions[0]?.version !== "latest" || versions[0]?.url !== "https://kit.cosmol.org/") return;
        const urls = versions.map(entry => new URL(entry.url));
        if (!versions.every((entry, i) => typeof entry.version === "string" && entry.version.trim()
            && typeof entry.url === "string" && urls[i].protocol === "https:"
            && !urls[i].username && !urls[i].password && urls[i].pathname === "/" && !urls[i].search && !urls[i].hash)
            || new Set(versions.map(entry => entry.version)).size !== versions.length
            || new Set(urls.map(url => url.href)).size !== versions.length) return;
        renderVersions(versions);
    } catch (_) { /* Offline/invalid catalog: leave the static menu usable. */ }
}

refreshVersions();
