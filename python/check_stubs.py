"""Check runtime samples and types against the generated package; never build it."""
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def main():
    output = HERE / "artifacts/smoke"
    package = output / "package/cosmolkit"
    if not (package / "__init__.pyi").is_file():
        raise SystemExit("Run python/build_stubs.py before checking stubs.")

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


if __name__ == "__main__":
    main()
