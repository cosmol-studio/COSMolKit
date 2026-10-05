// Execute the actual browser module with a small DOM/fetch boundary. No copies
// of production version-selection logic, network or third-party dependencies.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const root = path.join(__dirname, "..");
const source = fs.readFileSync(path.join(root, "assets/version-switch.js"), "utf8");
const catalog = JSON.parse(fs.readFileSync(path.join(root, "versions.json"), "utf8"));

async function run(origin, response = catalog, fail = false, ok = true, invalidJson = false) {
    class Link {
        constructor(entry = {}) {
            this.dataset = {docsVersion: entry.version};
            this.attrs = {href: entry.url};
            this.textContent = entry.version;
        }
        get href() { return new URL(this.attrs.href, origin).href; }
        setAttribute(key, value) { this.attrs[key] = value; }
    }
    const options = {
        links: catalog.versions.map(entry => new Link(entry)),
        querySelectorAll() { return this.links; },
        replaceChildren(...links) { this.links = links; },
    };
    const current = {textContent: "latest"};
    const requests = [];
    const context = vm.createContext({
        URL, location: {origin},
        document: {
            getElementById: id => ({"docs-version-options": options, "docs-version-current": current})[id],
            createElement: tag => { assert.equal(tag, "a"); return new Link(); },
        },
        fetch: async (url, init) => {
            requests.push({url, init});
            if (fail) throw new Error("offline");
            return {ok, json: async () => {
                if (invalidJson) throw new SyntaxError("invalid JSON");
                return response;
            }};
        },
    });
    // Wait for the module's real startup call, not a duplicate invocation.
    await vm.runInContext(source, context);
    assert.equal(requests.length, 1);
    assert.equal(requests[0].url, "https://kit.cosmol.org/versions.json");
    assert.equal(requests[0].init.credentials, "omit");
    assert.equal(requests[0].init.cache, "no-store");
    assert.equal(requests[0].init.mode, "cors");
    return {options, current};
}

(async () => {
    const latest = await run("https://kit.cosmol.org");
    assert.equal(latest.current.textContent, "latest");
    assert.equal(latest.options.links[0].attrs.href, "/");
    assert.equal(latest.options.links[1].attrs.href, catalog.versions[1].url);
    const archive = await run("https://c6862989.cosmolkit-docs-web.pages.dev");
    assert.equal(archive.current.textContent, "0.3.0");
    assert.equal(archive.options.links[1].attrs.href, "/");
    assert.equal(archive.options.links[1].attrs["aria-current"], "true");
    assert.equal(archive.options.links[0].attrs.href, catalog.versions[0].url);
    const expanded = structuredClone(catalog);
    expanded.versions.push({version: "0.4.0", url: "https://another-snapshot.pages.dev/"});
    const refreshed = await run("https://c6862989.cosmolkit-docs-web.pages.dev", expanded);
    assert.equal(refreshed.options.links.length, 3);
    assert.equal(refreshed.current.textContent, "0.3.0");
    assert.equal(refreshed.options.links[2].textContent, "0.4.0");
    const replaced = structuredClone(catalog);
    replaced.versions[1].url = "https://replacement-snapshot.pages.dev/";
    const olderSnapshot = await run("https://c6862989.cosmolkit-docs-web.pages.dev", replaced);
    assert.equal(olderSnapshot.current.textContent, "0.3.0");
    const offline = await run("https://c6862989.cosmolkit-docs-web.pages.dev", null, true);
    assert.equal(offline.current.textContent, "0.3.0");
    assert.equal(offline.options.links.length, 2);
    for (const result of [
        await run("https://c6862989.cosmolkit-docs-web.pages.dev", catalog, false, false),
        await run("https://c6862989.cosmolkit-docs-web.pages.dev", null, false, true, true),
    ]) {
        assert.equal(result.current.textContent, "0.3.0");
        assert.equal(result.options.links.length, 2);
    }
    for (const bad of [null, {}, {schema_version: 1, versions: []},
        {...catalog, versions: [...catalog.versions, {version: "evil", url: "javascript:alert(1)"}]},
        {...catalog, versions: [...catalog.versions, catalog.versions[1]]},
        {...catalog, versions: [{version: "latest", url: "https://wrong.example/"}]},
    ]) {
        const result = await run("https://c6862989.cosmolkit-docs-web.pages.dev", bad);
        assert.equal(result.options.links.length, 2);
        assert.equal(result.current.textContent, "0.3.0");
    }
    const local = await run("http://localhost:8080");
    assert.equal(local.current.textContent, "latest");
    console.log("PASS: latest/archive selection, new shared versions, offline/invalid fallback, relative current link, one catalog-only request");
})().catch(error => { console.error(error); process.exitCode = 1; });
