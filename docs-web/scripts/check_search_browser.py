"""Exercise the final static site's search UI in a local headless browser."""

import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Thread
from time import sleep
from urllib.parse import quote, urlsplit

from playwright.sync_api import sync_playwright

from check_ssg_output import check_output
from route_contract import PAGES, counterpart


class StaticPagesHandler(SimpleHTTPRequestHandler):
    def do_GET(self):
        # Local test harness for the emitted Pages rules, not a production router.
        rules = {}
        for line in (Path(self.directory) / "_redirects").read_text().splitlines():
            if line and not line.startswith("#"):
                source, destination, status = line.split()
                rules[source] = (destination, int(status))
        parsed = urlsplit(self.path)
        if parsed.path in rules:
            destination, status = rules[parsed.path]
            self.send_response(status)
            self.send_header("Location", destination + ("?" + parsed.query if parsed.query else ""))
            self.end_headers()
            return
        super().do_GET()

    def translate_path(self, path):
        candidate = Path(super().translate_path(path))
        if not candidate.suffix and candidate.with_suffix(".html").is_file():
            candidate = candidate.with_suffix(".html")
        return str(candidate)

    def log_message(self, format, *args):
        pass


def exercise_search(page, base_url, client=False):
    if client:
        page.add_init_script("""(() => {
            const add = EventTarget.prototype.addEventListener;
            const remove = EventTarget.prototype.removeEventListener;
            window.searchListeners = [];
            EventTarget.prototype.addEventListener = function(type, fn, options) {
                if (['docs-query', 'docs-search-form', 'search-results'].includes(this.id)
                    || (this === window && type === 'popstate' && document.getElementById('docs-search-form')))
                    searchListeners.push({target:this, type, fn, active:true});
                return add.call(this, type, fn, options);
            };
            EventTarget.prototype.removeEventListener = function(type, fn, options) {
                for (const item of searchListeners)
                    if (item.target === this && item.type === type && item.fn === fn) item.active = false;
                return remove.call(this, type, fn, options);
            };
        })();""")
    errors = []
    search_requests = []
    page.on("request", lambda request: search_requests.append(request.url) if "docs_search" in request.url else None)
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(base_url + "/")
    page.wait_for_selector('header a[href^="/python/search"]')
    assert not search_requests, search_requests
    page.evaluate("window.searchNavigationProbe = 42")
    page.locator('header a[href^="/python/search"]').click()
    page.wait_for_selector('#docs-search-form[data-ready="true"]')
    assert len([url for url in search_requests if url.endswith(".wasm")]) == 1, search_requests
    if client:
        assert page.evaluate("window.searchNavigationProbe") == 42, "expected client navigation"
    checked_pages = set()
    for query, has_results in (("Molecule", True), ("fingerprint", True), ("from_smiles", True), ("cosmolkitmissingsearchtoken", False), ("<script>alert(1)</script>", False)):
        page.locator('#docs-query').fill(query)
        page.wait_for_function("q => document.querySelector('#search-results')?.dataset.query === q", arg=query)
        count = page.locator("#search-results li").count()
        assert bool(count) == has_results, (query, page.locator('#docs-search-status').inner_text())
        assert count <= 30
        while page.locator('.docs-search-more').is_visible():
            page.locator('.docs-search-more').click()
        links = page.locator("#search-results li > a").evaluate_all("nodes => nodes.map(node => node.href)")
        for link in set(links):
            parsed = urlsplit(link)
            assert (parsed.path in ("/python", "/javascript") or parsed.path.startswith(("/python/", "/javascript/"))) and link.startswith(base_url + "/") and not parsed.path.endswith(".html"), link
            destination = parsed._replace(fragment="").geturl()
            # Hundreds of symbol anchors share one API document. Validate each
            # link above, but do not download the same multi-MB page per anchor.
            if destination not in checked_pages:
                assert page.request.get(destination).status == 200, link
                checked_pages.add(destination)
        assert not errors, errors
        print(f"PASS: search {query!r}: {len(links)} results; pagination and clean links")
    # Rapid input must discard stale work, and clearing must cancel pending input.
    page.locator('#docs-query').fill('Molecule')
    page.locator('#docs-query').fill('fingerprint')
    page.wait_for_function("document.querySelector('#search-results')?.dataset.query === 'fingerprint'")
    page.locator('#docs-query').fill('Molecule')
    page.locator('#docs-search-form button[type="reset"]').click()
    assert page.locator('#docs-query').input_value() == ''
    assert page.locator('#search-results li').count() == 0
    initial_requests = len(search_requests)
    page.locator('header a[href="/python/api"]').click()
    if client:
        page.wait_for_function("searchListeners.length >= 5 && searchListeners.every(item => !item.active)")
        print('PASS: leaving search releases form, results, and history listeners')
    page.locator('header a[href^="/python/search"]').click()
    page.wait_for_selector('#docs-search-form[data-ready="true"]')
    if client:
        assert len(search_requests) == initial_requests, "search module should be reused on client navigation"
    print('PASS: search WASM is lazy, input updates results, and stale queries are cancelled')
    page.goto(base_url + '/python/search?q=Molecule')
    page.wait_for_selector('#docs-search-form[data-ready="true"]')
    assert page.locator('#search-results li').count() == 30
    page.locator('#docs-query').fill('fingerprint')
    page.locator('#docs-query').press('Enter')
    page.go_back()
    page.wait_for_function("document.querySelector('#docs-query')?.value === 'Molecule'")
    assert page.locator('#search-results ul').count() == 1
    page.locator('header a[href="/python/api"]').click()
    page.wait_for_url('**/python/api')
    page.locator('header a[href^="/python/search"]').click()
    page.wait_for_selector('#docs-search-form[data-ready="true"]')
    page.locator('#docs-query').fill('Molecule')
    page.wait_for_function("document.querySelector('#search-results')?.dataset.query === 'Molecule'")
    assert page.locator('#search-results ul').count() == 1
    assert page.locator('#search-results li').count() == 30
    assert not errors, errors
    first = page.locator('#search-results li > a').first
    destination = first.get_attribute('href')
    first.click()
    page.wait_for_url(base_url + destination)
    fragment = urlsplit(destination).fragment
    if fragment:
        page.wait_for_function("id => document.getElementById(id) !== null", arg=fragment)
    print('PASS: query URL, clear, back navigation, and re-entering search')
    for query in ('molecular fingerprint', 'Molecule & atom', '\u5206\u5b50', '<script>alert(1)</script>'):
        page.goto(base_url + '/python/search?q=' + quote(query))
        page.wait_for_selector('#docs-search-form[data-ready="true"]')
        assert page.locator('#docs-query').input_value() == query
        assert not errors, errors
    print('PASS: bookmarked queries preserve spaces, ampersands, Unicode, and markup as text')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("public", type=Path, nargs="?")
    parser.add_argument("--url", help="Test an existing dx serve URL")
    parser.add_argument("--channel", help="Use an installed browser, for example msedge")
    args = parser.parse_args()
    if args.url:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(channel=args.channel)
            try:
                exercise_api_layout(browser.new_page(), args.url.rstrip('/'))
                exercise_listener_cleanup(browser.new_page(), args.url.rstrip('/'))
                exercise_navigation(browser.new_page(), args.url.rstrip('/'))
                exercise_search(browser.new_page(), args.url.rstrip('/'), client=True)
                exercise_load_failure(browser, args.url.rstrip('/'))
            finally:
                browser.close()
        return
    if not args.public:
        parser.error("provide PUBLIC_DIR or --url")
    check_output(args.public)
    server = ThreadingHTTPServer(("127.0.0.1", 0), partial(StaticPagesHandler, directory=str(args.public.resolve())))
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base_url = f"http://127.0.0.1:{server.server_port}"
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(channel=args.channel)
            try:
                page = browser.new_page()
                errors = []
                page.on("pageerror", lambda error: errors.append(str(error)))
                for route in PAGES:
                    response = page.goto(base_url + route["path"])
                    assert response.status == 200, route["path"]
                    switch = page.locator(".docs-language-switch")
                    if route["binding"] == "common":
                        assert switch.locator(".is-active").count() == 0
                    elif route["binding"] == "python":
                        assert urlsplit(switch.locator('a[aria-current="page"]').get_attribute("href")).path == route["path"]
                        destination = counterpart(PAGES, route, "javascript")
                        expected = destination["path"] if destination else "/javascript"
                        assert switch.locator(f'a[href="{expected}"]').count() == 1
                    assert not errors, errors
                for legacy in ("/api", "/api.html", "/api/"):
                    response = page.request.get(base_url + legacy + "?from=legacy", max_redirects=0)
                    assert response.status == 301
                    assert response.headers["location"] == "/python/api?from=legacy"
                    page.goto(base_url + legacy + "?from=legacy#cosmolkit.Molecule")
                    assert page.url == base_url + "/python/api?from=legacy#cosmolkit.Molecule"
                    assert page.locator('[id="cosmolkit.Molecule"]').count() == 1
                print("PASS: all declared pages and binding navigation; legacy redirects preserve query and fragment")
                exercise_api_layout(browser.new_page(), base_url)
                # Navigation and anchor positioning must also work with scripts disabled.
                static_context = browser.new_context(java_script_enabled=False)
                exercise_navigation(static_context.new_page(), base_url)
                static_context.close()
                exercise_search(page, base_url)
                exercise_load_failure(browser, base_url)
                browser.close()
            finally:
                if browser.is_connected():
                    browser.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def exercise_api_layout(page, base_url):
    """Compare real computed API styles, not just shared CSS class names."""
    try:
        for width in (1440, 1000, 390):
            page.set_viewport_size({"width": width, "height": 1000})
            snapshots = {}
            for language in ("python", "javascript"):
                page.goto(base_url + f"/{language}/api", wait_until="networkidle")
                snapshots[language] = page.evaluate("""() => {
                    const article = document.querySelector('.docs-article');
                    const accordions = [...article.querySelectorAll('details.api-members')];
                    if (!accordions.length || accordions.some(el => el.open))
                        throw new Error('Class members must be collapsed by default');
                    const entries = [...article.querySelector('section').querySelectorAll(':scope > dl')];
                    const kinds = entries.map(el => el.classList.contains('function') ? 'function' :
                        ['class', 'interface', 'exception'].some(kind => el.classList.contains(kind)) ? 'class' : 'other');
                    const firstClass = kinds.indexOf('class');
                    if (firstClass < 1 || kinds.slice(0, firstClass).some(kind => kind !== 'function') ||
                            kinds.slice(firstClass).includes('function'))
                        throw new Error('Top-level functions must precede classes');
                    const classes = entries.filter((el, i) => kinds[i] === 'class');
                    const names = classes.map(el => el.querySelector(':scope > dt .sig-name')?.textContent);
                    if (JSON.stringify(names.slice(0, 3)) !== JSON.stringify(['Molecule', 'BioStructure', 'Protein']))
                        throw new Error(`Primary class order differs: ${names.slice(0, 3)}`);
                    // Open only for computed-style comparisons; the initial state was checked above.
                    for (const accordion of accordions) accordion.open = true;
                    const cls = [...article.querySelectorAll('dl.class')].find(el =>
                        el.querySelector(':scope > dt .sig-name')?.textContent === 'Molecule');
                    const apiBlocks = [...article.querySelectorAll(
                        'dl[class]:not(.field-list):not(.option-list):not(.simple)')];
                    for (const block of apiBlocks)
                        if (getComputedStyle(block).margin !== '0px')
                            throw new Error(`API block margin differs from class: ${block.querySelector('dt')?.id}`);
                    const method = cls?.querySelector(':scope > dd > .api-members > .api-member-list > dl.method, :scope > dd > .api-members > .api-member-list > dl.function');
                    const attributes = [...article.querySelectorAll(
                        '.api-member-list > dl.property, .api-member-list > dl.attribute')];
                    const members = [...article.querySelectorAll('.api-member-list > dl')]
                        .filter(member => member.querySelector(':scope > dt.sig'));
                    const spacing = node => [node, node.querySelector(':scope > dt'), node.querySelector(':scope > dd')]
                        .map(el => ({margin: getComputedStyle(el).margin, padding: getComputedStyle(el).padding}));
                    if (!attributes.length) throw new Error('Missing API property samples');
                    const lastMethod = article.querySelector('.api-member-list > dl:last-child');
                    const properties = ['display', 'position', 'top', 'left', 'max-height', 'overflow-x',
                        'overflow-y', 'margin', 'padding', 'border-left', 'border-radius',
                        'font-family', 'font-size', 'font-weight', 'font-style', 'line-height', 'color',
                        'background-color', 'overflow-wrap', 'scrollbar-width', 'scrollbar-color',
                        'overscroll-behavior'];
                    const style = (el, pseudo = null) => {
                        if (!el) throw new Error('Missing API layout sample');
                        const computed = getComputedStyle(el, pseudo);
                        return Object.fromEntries(properties.map(key => [key, computed.getPropertyValue(key)]));
                    };
                    const memberStyle = node => [node, node.querySelector(':scope > dt'), node.querySelector(':scope > dd')]
                        .map(el => {
                            const computed = getComputedStyle(el);
                            return Object.fromEntries(['margin', 'padding', 'border', 'border-radius',
                                'background-color', 'color', 'font-family', 'font-size', 'font-weight',
                                'line-height'].map(key => [key, computed.getPropertyValue(key)]));
                        });
                    const baseline = JSON.stringify(memberStyle(method));
                    for (const member of members)
                        if (JSON.stringify(memberStyle(member)) !== baseline)
                            throw new Error(`Peer member styling differs: ${member.querySelector('dt').id}`);
                    const toc = document.querySelector('.docs-toc .docs-toc-tree');
                    const classLink = [...toc.querySelectorAll('a')].find(a => a.getAttribute('href')?.endsWith('.Molecule'));
                    const depth = el => {let n = 0; for (; el && el !== toc; el = el.parentElement) if (el.tagName === 'UL') n++; return n;};
                    const expectedIds = entries.filter((el, i) => kinds[i] !== 'other')
                        .map(el => '#' + el.querySelector(':scope > dt').id);
                    for (const tree of document.querySelectorAll('.docs-toc-tree')) {
                        const links = [...tree.querySelectorAll('a')];
                        if (tree.querySelector('ul ul') || JSON.stringify(links.map(a => a.getAttribute('href'))) !== JSON.stringify(expectedIds))
                            throw new Error('API TOC must list only top-level functions and classes in body order');
                    }
                    if (depth(classLink) !== 1) throw new Error('Redundant API page-title level');
                    const geometry = el => {const r = el.getBoundingClientRect(); return {...style(el), x: r.x, width: r.width};};
                    const shell = ['.docs-layout', '.docs-main', '.docs-article', '.docs-toc',
                        '.docs-toc-tree', '.docs-mobile-toc'].map(selector => {
                            const el = document.querySelector(selector), rect = el.getBoundingClientRect();
                            return [selector, {...style(el), x: rect.x, width: rect.width}];
                        });
                    const continuation = style(lastMethod, '::after');
                    // top:100% resolves to each method's content-dependent height.
                    // Require that boundary in both languages, then compare the style.
                    if (continuation.top !== getComputedStyle(lastMethod).height)
                        throw new Error('Method guide line does not continue from the bottom edge');
                    continuation.top = '100%';
                    return {shell, class: style(cls), signature: style(cls.querySelector(':scope > dt')),
                        body: style(cls.querySelector(':scope > dd')), method: style(method),
                        methodSignature: style(method.querySelector(':scope > dt')),
                        methodBody: style(method.querySelector(':scope > dd')),
                        propertySpacing: spacing(attributes[0]),
                        memberStyle: memberStyle(method),
                        continuation: {...continuation,
                            content: getComputedStyle(lastMethod, '::after').content,
                            bottom: getComputedStyle(lastMethod, '::after').bottom},
                        tocClass: geometry(classLink), tocDepth: depth(classLink),
                        overflows: document.documentElement.scrollWidth > innerWidth};
                }""")
            assert not snapshots["python"]["overflows"], (width, "Python page overflow")
            assert not snapshots["javascript"]["overflows"], (width, "JavaScript page overflow")
            # Python's existing, deliberately compact API spacing is the
            # baseline. Equality alone could let both languages regress.
            for name in ("class", "signature", "method", "methodSignature"):
                assert snapshots["python"][name]["margin"] == "0px", (
                    width, "original Python API margin", name, snapshots["python"][name])
            for name in ("body", "methodBody"):
                assert snapshots["python"][name]["margin"] == "0px 0px 0px 16px"
                assert snapshots["python"][name]["padding"] == "4px 0px 16px 16px"
            for name in snapshots["python"]:
                assert snapshots["python"][name] == snapshots["javascript"][name], (
                    width, name, snapshots["python"][name], snapshots["javascript"][name])
        # Native details must keep existing member deep links usable.
        for language in ('python', 'javascript'):
            page.goto(base_url + f'/{language}/api', wait_until='networkidle')
            owner = page.locator('dl.class').filter(has=page.locator('dt > .sig-name', has_text='Molecule')).first
            accordion = owner.locator(':scope > dd > details.api-members')
            assert not accordion.evaluate('(el) => el.open')
            accordion.locator(':scope > summary').click()
            assert accordion.evaluate('(el) => el.open')
            anchor = accordion.locator('.api-member-list > dl > dt[id]').first.get_attribute('id')
            accordion.locator(':scope > summary').click()
            assert not accordion.evaluate('(el) => el.open')
            page.goto(base_url + f'/{language}/api#{anchor}', wait_until='networkidle')
            page.wait_for_function('(id) => document.getElementById(id)?.closest("details")?.open', arg=anchor)
            page.reload(wait_until='networkidle')
            page.wait_for_function('(id) => document.getElementById(id)?.closest("details")?.open', arg=anchor)
        print('PASS: functions first, prioritized classes, collapsed members and flat TOC; member links work and Python/JavaScript styles match at desktop, tablet and mobile widths')
    finally:
        page.close()


def exercise_load_failure(browser, base_url):
    for pattern in ('**/docs_search-*.js', '**/docs_search_bg-*.wasm'):
        page = browser.new_page()
        try:
            page.route(pattern, lambda route: route.abort())
            page.goto(base_url + '/python/search?q=Molecule')
            page.wait_for_function("document.querySelector('#docs-search-status')?.textContent.includes('could not load')")
            page.unroute(pattern)
            page.reload()
            page.wait_for_selector('#docs-search-form[data-ready="true"]')
            assert page.locator('#search-results li').count() > 0
        finally:
            page.close()
    print('PASS: binding/WASM load failures are visible and recover after reload')


def exercise_listener_cleanup(page, base_url):
    # Delay CSS so this actually exercises the deferred scroll and its listeners.
    def delayed_css(route):
        sleep(0.1)
        route.continue_()
    page.route('**/*.css', delayed_css)
    page.add_init_script("""(() => {
        const add = EventTarget.prototype.addEventListener;
        const remove = EventTarget.prototype.removeEventListener;
        const listeners = [];
        window.anchorListeners = listeners;
        EventTarget.prototype.addEventListener = function(type, fn, options) {
            if (this instanceof HTMLLinkElement && type === 'load')
                listeners.push({target:this, fn, active:true});
            return add.call(this, type, fn, options);
        };
        EventTarget.prototype.removeEventListener = function(type, fn, options) {
            if (type === 'load')
                for (const item of listeners)
                    if (item.target === this && item.fn === fn) item.active = false;
            return remove.call(this, type, fn, options);
        };
    })();""")
    page.goto(base_url + '/python/api#cosmolkit.Molecule.atom_metadata')
    page.wait_for_function("anchorListeners.length > 0 && anchorListeners.every(item => !item.active)")
    page.wait_for_function("""() => {
        const y = document.getElementById('cosmolkit.Molecule.atom_metadata').getBoundingClientRect().top;
        return y >= 74 && y < innerHeight;
    }""")
    page.evaluate("""() => {
        scrollTo(0, 0);
        for (const item of anchorListeners) item.target.dispatchEvent(new Event('load'));
    }""")
    assert page.evaluate('scrollY') == 0, 'completed listeners must not scroll again'
    print('PASS: delayed stylesheet anchor positioning releases all load listeners')
    page.close()


def exercise_navigation(page, base_url):
    anchor = 'cosmolkit.Molecule.atom_metadata'
    target = '/python/api#' + anchor
    def assert_position(expected_anchor=anchor):
        page.wait_for_function("""id => {
            const target = document.getElementById(id);
            const header = document.querySelector('header');
            if (!target || !header) return false;
            const y = target.getBoundingClientRect().top;
            return y >= header.getBoundingClientRect().bottom && y < innerHeight;
        }""", arg=expected_anchor)

    page.goto(base_url + target)
    assert_position()
    page.reload()
    assert_position()
    page.goto(base_url + '/python/api')
    page.locator('.docs-toc a[href="#cosmolkit.Molecule"]').click()
    assert_position('cosmolkit.Molecule')
    page.locator('dl.class > dt[id="cosmolkit.Molecule"] + dd > details.api-members > summary').click()
    page.locator('[id="' + anchor + '"] > a.headerlink').click()
    assert_position()
    switch = page.locator('.docs-language-switch')
    assert switch.locator('a').count() == 2
    switch.locator('a[href="/javascript/api"]').click()
    page.wait_for_url(base_url + '/javascript/api')
    assert page.locator('.docs-language-switch a[aria-current="page"]').inner_text() == 'JavaScript'
    page.locator('.docs-language-switch a[href="/python/api"]').click()
    page.wait_for_url(base_url + '/python/api')
    page.go_back()
    page.wait_for_url(base_url + '/javascript/api')
    print('PASS: native language links and direct, reload, and on-page API anchor positioning')
    page.close()


if __name__ == "__main__":
    main()
