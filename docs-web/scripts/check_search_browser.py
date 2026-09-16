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
from route_contract import PAGES


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
            assert link.startswith(base_url + "/python") and not parsed.path.endswith(".html"), link
            assert page.request.get(link).status == 200, link
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
                        assert switch.locator('a[href="/javascript"]').count() == 1
                    assert not errors, errors
                for legacy in ("/api", "/api.html", "/api/"):
                    response = page.request.get(base_url + legacy + "?from=legacy", max_redirects=0)
                    assert response.status == 301
                    assert response.headers["location"] == "/python/api?from=legacy"
                    page.goto(base_url + legacy + "?from=legacy#cosmolkit.Molecule")
                    assert page.url == base_url + "/python/api?from=legacy#cosmolkit.Molecule"
                    assert page.locator('[id="cosmolkit.Molecule"]').count() == 1
                print("PASS: all declared pages and binding navigation; legacy redirects preserve query and fragment")
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
    page.goto(base_url + '/python/api#cosmolkit.Molecule.fingerprint_atom_pair')
    page.wait_for_function("anchorListeners.length > 0 && anchorListeners.every(item => !item.active)")
    page.wait_for_function("""() => {
        const y = document.getElementById('cosmolkit.Molecule.fingerprint_atom_pair').getBoundingClientRect().top;
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
    anchor = 'cosmolkit.Molecule.fingerprint_atom_pair'
    target = '/python/api#' + anchor
    def assert_position():
        page.wait_for_function("""id => {
            const target = document.getElementById(id);
            const header = document.querySelector('header');
            if (!target || !header) return false;
            const y = target.getBoundingClientRect().top;
            return y >= header.getBoundingClientRect().bottom && y < innerHeight;
        }""", arg=anchor)

    page.goto(base_url + target)
    assert_position()
    page.reload()
    assert_position()
    page.goto(base_url + '/python/api')
    page.locator('.docs-toc a[href="#' + anchor + '"]').click()
    assert_position()
    switch = page.locator('.docs-language-switch')
    assert switch.locator('a').count() == 2
    switch.locator('a[href="/javascript"]').click()
    page.wait_for_url(base_url + '/javascript')
    assert page.locator('.docs-language-switch a[aria-current="page"]').inner_text() == 'JavaScript'
    page.locator('.docs-language-switch a[href="/python"]').click()
    page.wait_for_url(base_url + '/python')
    page.go_back()
    page.wait_for_url(base_url + '/javascript')
    print('PASS: native language links and direct, reload, and on-page API anchor positioning')
    page.close()


if __name__ == "__main__":
    main()
