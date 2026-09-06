"""Flatten Dioxus SSG directories for legacy .html documentation URLs."""

from pathlib import Path
import shutil
import sys


def flatten_html_routes(public_dir: Path) -> int:
    flattened = 0
    for route_dir in sorted(public_dir.iterdir()):
        if not route_dir.is_dir() or not route_dir.name.endswith(".html"):
            continue

        index = route_dir / "index.html"
        if not index.is_file():
            continue

        children = list(route_dir.iterdir())
        if children != [index]:
            names = ", ".join(child.name for child in children)
            raise ValueError(
                f"route directory contains unexpected files: {route_dir} ({names})"
            )

        destination = public_dir / route_dir.name
        temporary = public_dir / f".{route_dir.name}.tmp"
        if temporary.exists():
            raise FileExistsError(f"temporary route already exists: {temporary}")
        existing = destination.is_file()

        shutil.move(str(index), str(temporary))
        route_dir.rmdir()
        if existing:
            if destination.read_bytes() != temporary.read_bytes():
                raise FileExistsError(
                    f"generated route conflicts with existing file: {destination}"
                )
            temporary.unlink()
        else:
            shutil.move(str(temporary), str(destination))
        flattened += 1

    return flattened


def main() -> None:
    default_dir = "target/dx/cosmolkit-docs-web/release/web/public"
    public_dir = Path(sys.argv[1] if len(sys.argv) > 1 else default_dir)
    if not public_dir.is_dir():
        raise SystemExit(f"public directory does not exist: {public_dir}")
    print(f"Flattened {flatten_html_routes(public_dir)} HTML routes in {public_dir}")


if __name__ == "__main__":
    main()
