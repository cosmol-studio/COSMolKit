"""Flatten Dioxus SSG directory indexes for Cloudflare Pages clean URLs."""

from pathlib import Path
import shutil
import sys


ROUTES = (
    "installation", "quickstart", "confseq", "molecule", "batch",
    "fingerprints", "descriptors", "protein", "io", "api", "search",
    "genindex", "py-modindex", "javascript", "python", "benchmarks", "validation",
)


def flatten_html_routes(public_dir: Path) -> int:
    # Only documentation routes belong here; Sphinx _modules and assets retain
    # their layout. Preflight every route before moving any content.
    pending = []
    for route in ROUTES:
        route_dir = public_dir / route
        destination = public_dir / f"{route}.html"
        if not route_dir.exists():
            continue
        index = route_dir / "index.html"
        if not route_dir.is_dir() or not index.is_file():
            raise ValueError(f"route must contain an SSG directory index: {route_dir}")
        if list(route_dir.iterdir()) != [index]:
            raise ValueError(f"route directory contains unexpected files: {route_dir}")
        if destination.exists():
            raise FileExistsError(f"cannot flatten {route_dir}: {destination} exists")
        pending.append((route_dir, index, destination))

    for route_dir, index, destination in pending:
        shutil.move(str(index), str(destination))
        route_dir.rmdir()
    return len(pending)


def main() -> None:
    default_dir = "target/dx/cosmolkit-docs-web/release/web/public"
    public_dir = Path(sys.argv[1] if len(sys.argv) > 1 else default_dir)
    if not public_dir.is_dir():
        raise SystemExit(f"public directory does not exist: {public_dir}")
    print(f"Flattened {flatten_html_routes(public_dir)} SSG routes in {public_dir}")


if __name__ == "__main__":
    main()
