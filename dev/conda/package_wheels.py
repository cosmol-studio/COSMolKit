"""Package release wheels as platform-specific conda archives."""

import argparse
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from conda_package_handling.api import create
from installer import install
from installer.destinations import SchemeDictionaryDestination
from installer.sources import WheelFile
from packaging.requirements import Requirement
from packaging.utils import parse_wheel_filename


RUNTIMES = {
    "linux-64": ["__glibc >=2.17", "libgcc-ng"],
    "linux-aarch64": ["__glibc >=2.17", "libgcc-ng"],
    "osx-64": ["__osx >=10.12"],
    "osx-arm64": ["__osx >=11.0"],
    "win-64": ["vc14_runtime", "ucrt"],
}


def subdirs(platform):
    if platform.startswith("manylinux"):
        if platform.endswith("aarch64"):
            return ["linux-aarch64"]
        if platform.endswith("x86_64"):
            return ["linux-64"]
    if platform.startswith("macosx"):
        if platform.endswith("universal2"):
            return ["osx-64", "osx-arm64"]
        if platform.endswith("x86_64"):
            return ["osx-64"]
        if platform.endswith("arm64"):
            return ["osx-arm64"]
    if platform == "win_amd64":
        return ["win-64"]
    raise ValueError(f"Unsupported wheel platform: {platform}")


def package(wheel, subdir, minor, output):
    python = f"3.{minor}"
    with TemporaryDirectory(prefix="cosmolkit-conda-") as folder:
        root = Path(folder)
        site = "Lib/site-packages" if subdir == "win-64" else f"lib/python{python}/site-packages"
        destination = SchemeDictionaryDestination(
            scheme_dict={
                "purelib": str(root / site), "platlib": str(root / site),
                "scripts": str(root / ("Scripts" if subdir == "win-64" else "bin")),
                "headers": str(root / "include"), "data": str(root),
            },
            interpreter="python.exe" if subdir == "win-64" else "/usr/bin/env python",
            script_kind="win-amd64" if subdir == "win-64" else "posix",
        )
        with WheelFile.open(wheel) as source:
            metadata = BytesParser().parsebytes(source.read_dist_info("METADATA").encode())
            install(source, destination, {"INSTALLER": b"conda\n"})
        requirements = [Requirement(value) for value in metadata.get_all("Requires-Dist") or []]
        if any(req.marker or req.extras or req.url for req in requirements):
            raise ValueError(f"Wheel dependencies require an explicit conda mapping: {wheel.name}")
        files = sorted(path.relative_to(root).as_posix() for path in root.rglob("*") if path.is_file())
        paths = []
        for name in files:
            path = root / name
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            paths.append({"_path": name, "path_type": "hardlink", "sha256": digest, "size_in_bytes": path.stat().st_size})
        info = root / "info"
        info.mkdir()
        index = {
            "name": "cosmolkit", "version": metadata["Version"],
            "build": f"py3{minor}_0", "build_number": 0, "subdir": subdir,
            "platform": subdir.split("-")[0],
            "arch": "x86_64" if subdir.endswith("-64") else subdir.split("-")[1],
            "license": "MIT",
            "depends": [f"python >=3.{minor},<3.{minor + 1}.0a0", f"python_abi 3.{minor}.* *_cp3{minor}",
                        *(f"{req.name} {req.specifier}".strip() for req in requirements), *RUNTIMES[subdir]],
        }
        (info / "index.json").write_text(json.dumps(index), encoding="utf-8")
        (info / "paths.json").write_text(json.dumps({"paths_version": 1, "paths": paths}), encoding="utf-8")
        (info / "files").write_text("\n".join(files) + "\n", encoding="utf-8")
        (info / "about.json").write_text(json.dumps({"license": "MIT", "summary": metadata["Summary"],
                                                    "home": "https://tools.cosmol.org/"}), encoding="utf-8")
        (info / "licenses").mkdir()
        shutil.copyfile(Path(__file__).resolve().parents[2] / "LICENSE", info / "licenses/LICENSE")
        target = output / subdir
        target.mkdir(parents=True, exist_ok=True)
        members = sorted(path.relative_to(root).as_posix() for path in root.rglob("*") if path.is_file())
        name = f"cosmolkit-{index['version']}-{index['build']}.conda"
        create(str(root), members, name, out_folder=str(target), compression_level=6)
        print(f"{wheel.name} -> {subdir}/{name}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheels", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    selected = {}
    for wheel in sorted(args.wheels.glob("cosmolkit-*.whl")):
        _, _, _, tags = parse_wheel_filename(wheel.name)
        for tag in tags:
            for subdir in subdirs(tag.platform):
                for minor in range(9, 15):
                    exact = tag.interpreter == f"cp3{minor}" and tag.abi == f"cp3{minor}"
                    stable = tag.abi == "abi3" and int(tag.interpreter[2:]) <= int(f"3{minor}")
                    if exact or stable:
                        key = (subdir, minor)
                        if key not in selected or exact:
                            selected[key] = (wheel, exact)
    if not selected:
        raise ValueError(f"No release wheels found in {args.wheels}")
    for (subdir, minor), (wheel, _) in sorted(selected.items()):
        package(wheel, subdir, minor, args.output)


if __name__ == "__main__":
    main()
