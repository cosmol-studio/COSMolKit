"""Flatten Dioxus SSG directory indexes for Cloudflare Pages clean URLs."""

from pathlib import Path
import shutil
import sys


from route_contract import ROUTES


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
        children = {public_dir / other for other in ROUTES if Path(other).parent == Path(route)}
        flat_children = {child.with_suffix(".html") for child in children}
        if not route_dir.is_dir():
            raise ValueError(f"route must contain an SSG directory index: {route_dir}")
        # A flattened parent still owns a directory containing its child pages
        # and language-specific Sphinx resources. Never remove that container.
        resources = {route_dir / name for name in ("_static", "_sources", "_modules", "_images", "_downloads", "objects.inv")} if children else set()
        if set(route_dir.iterdir()) - {index} - children - flat_children - resources:
            raise ValueError(f"route directory contains unexpected files: {route_dir}")
        if not index.is_file():
            if children:
                continue
            raise ValueError(f"route must contain an SSG directory index: {route_dir}")
        if destination.exists():
            raise FileExistsError(f"cannot flatten {route_dir}: {destination} exists")
        pending.append((route_dir, index, destination))

    for route_dir, index, destination in pending:
        shutil.move(str(index), str(destination))
    for route_dir, _, _ in sorted(pending, key=lambda item: len(item[0].parts), reverse=True):
        if not any(route_dir.iterdir()):
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
