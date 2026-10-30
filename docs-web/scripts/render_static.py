"""Render the route contract using the native Dioxus build; no browser bundle."""

import os
from pathlib import Path
import socket
import subprocess
import sys
import time
from urllib.error import URLError
from urllib.request import ProxyHandler, build_opener

from route_contract import PAGES


def render(binary):
    binary = Path(binary).resolve()
    public = binary.parent / "public"
    public.mkdir(parents=True, exist_ok=True)
    # Dioxus replaces this title with each route's document::Title. Its fallback
    # template has no title slot and places the rendered title before <html>.
    (public / "index.html").write_text(
        '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1.0">'
        '<title>COSMolKit Documentation</title></head>'
        '<body><div id="main"></div></body></html>', encoding="utf-8"
    )
    with socket.socket() as listener:
        listener.bind(("127.0.0.1", 0))
        port = listener.getsockname()[1]
    environment = dict(os.environ, IP="127.0.0.1", PORT=str(port))
    process = subprocess.Popen([str(binary)], cwd=binary.parent, env=environment)
    try:
        origin = f"http://127.0.0.1:{port}"
        http = build_opener(ProxyHandler({}))
        deadline = time.monotonic() + 60
        while True:
            if process.poll() is not None:
                raise RuntimeError(f"Static renderer exited: {process.returncode}")
            try:
                with http.open(origin, timeout=2) as response:
                    response.read()
                break
            except URLError:
                if time.monotonic() >= deadline:
                    raise TimeoutError("Static renderer did not start")
                time.sleep(0.1)
        for page in PAGES:
            route = page["path"]
            with http.open(origin + route, timeout=60) as response:
                html = response.read()
            destination = public / route.lstrip("/") / "index.html"
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(html)
        print(f"Rendered {len(PAGES)} static pages to {public}")
    finally:
        if process.poll() is None:
            process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: render_static.py NATIVE_BINARY")
    render(sys.argv[1])
