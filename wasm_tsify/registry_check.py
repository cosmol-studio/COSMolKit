"""Check the built JS package and declarations against its compiled Rust registry."""
import argparse
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", type=Path, default=ROOT / "pkg/cosmolkit_tsify_probe.js")
    parser.add_argument("--report", type=Path, default=ROOT / "target/registry.json")
    args = parser.parse_args()
    command = shlex.join(["node", str(ROOT / "registry_check.cjs"),
                          str(args.package.resolve()), str(args.report.resolve())])
    # npm supplies the same pinned TypeScript parser used by the type consumers.
    result = subprocess.run(["npx", "--yes", "--package", "typescript@5.8.3", "-c", command])
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
