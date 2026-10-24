"""Build once, generate the public package, then check ordinary usage."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from generate import generate

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extractor", default="pyo3-introspection")
    options = parser.parse_args()
    output = HERE / "artifacts/smoke"
    output.mkdir(parents=True, exist_ok=True)

    print("Build Rust extension", flush=True)
    result = subprocess.run([
        "cargo", "build", "--locked", "--manifest-path", str(HERE / "Cargo.toml"),
        "--target-dir", str(ROOT / "target/typed_facade_projection"),
        "--features", "extension", "--message-format=json",
    ], cwd=ROOT, env={**os.environ, "PYO3_PYTHON": sys.executable},
        text=True, stdout=subprocess.PIPE, check=True)
    artifacts = [json.loads(line) for line in result.stdout.splitlines()]
    extension = next(item for item in artifacts if item.get("reason") == "compiler-artifact"
                     and item["target"]["name"] == "_native")
    library = next(Path(name) for name in extension["filenames"]
                   if Path(name).suffix in {".so", ".dll", ".dylib"})

    print("Extract official stub", flush=True)
    official = output / "official"
    subprocess.run([options.extractor, str(library), "_native", str(official)], check=True)
    print("Generate public package", flush=True)
    package = generate(library, official / "__init__.pyi", output / "package")

    print("Check actual calls", flush=True)
    subprocess.run([sys.executable, str(HERE / "tests/runtime.py"), str(package.parent)], check=True)
    config = output / "typing.json"
    config.write_text(json.dumps({
        "include": [str(HERE / "tests/positive.py"), str(package / "__init__.pyi")],
        "extraPaths": [str(package.parent)], "pythonVersion": "3.12",
        "typeCheckingMode": "basic", "reportMissingModuleSource": False,
    }, indent=2))
    print("Check public usage types", flush=True)
    subprocess.run([sys.executable, "-m", "basedpyright", "--pythonpath", sys.executable,
                    "--project", str(config)], check=True)
    print(f"Generated package: {package}")


if __name__ == "__main__":
    main()
