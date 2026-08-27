#!/usr/bin/env python3
"""Prepare and validate generated RDKit expected data used by Rust tests."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from _expected_schema import SCHEMAS, validate_jsonl_output


REPO_ROOT = Path(__file__).resolve().parents[3]
TOOL_DIR = Path(__file__).resolve().parent
REFERENCE_IDENTITY_PATH = REPO_ROOT / "testdata/reference/rdkit.json"
SCHEMA_PATH = TOOL_DIR / "_expected_schema.py"
PROFILE_INPUTS = {
    "ciplabeler_focused": REPO_ROOT
    / "testdata/stereo/fixtures/ciplabeler_focused.json",
    "atom_pair_focused": REPO_ROOT
    / "testdata/fingerprint/fixtures/rdkit/atom_pair_fingerprint_focused.smi",
    "layered_focused": REPO_ROOT
    / "testdata/fingerprint/fixtures/rdkit/layered_fingerprint_focused.smi",
    "pattern_focused": REPO_ROOT
    / "testdata/fingerprint/fixtures/rdkit/pattern_fingerprint_focused.smi",
    "smiles_small": REPO_ROOT / "testdata/smiles/corpus/smiles_small.smi",
    "smiles_5000": REPO_ROOT / "testdata/smiles/corpus/smiles_5000.smi",
    "smarts_source": REPO_ROOT / "testdata/smarts/corpus/rdkit_source_cases.json",
    "tautomer_focused": REPO_ROOT
    / "testdata/tautomer/fixtures/rdkit/tautomer_focused_cases.json",
    "tautomer_pcs_1k": REPO_ROOT
    / "testdata/tautomer/corpus/rdkit/1kPCS_tautomer.csv.gz",
    "tautomer_pcs_100k": REPO_ROOT
    / "testdata/tautomer/corpus/rdkit/100kPCS_tautomer.csv.gz",
}


@dataclass(frozen=True)
class GeneratorSpec:
    script: str
    output: str
    domain: str
    suites: frozenset[str]
    extra_inputs: tuple[str, ...] = ()
    generator_dependencies: tuple[str, ...] = ()
    deterministic_shards: int | None = None
    profiles: frozenset[str] | None = None


def spec(
    script: str,
    output: str,
    domain: str,
    suites: set[str],
    *,
    extra_inputs: tuple[str, ...] = (),
    generator_dependencies: tuple[str, ...] = (),
    deterministic_shards: int | None = None,
    profiles: set[str] | None = None,
) -> GeneratorSpec:
    return GeneratorSpec(
        script=f"_generate_{script}.py",
        output=output,
        domain=domain,
        suites=frozenset({*suites, "all", domain}),
        extra_inputs=extra_inputs,
        generator_dependencies=generator_dependencies,
        deterministic_shards=deterministic_shards,
        profiles=None if profiles is None else frozenset(profiles),
    )


GENERATOR_SPECS = [
    spec("inchi_golden", "inchi.jsonl", "inchi", {"default", "strict-corpus"}),
    spec(
        "v2000_minimal_golden",
        "molblock_v2000_minimal.jsonl",
        "molblock",
        {"default", "strict-corpus"},
    ),
    spec(
        "kekulize_molblock_golden",
        "molblock_v2000_kekulized.jsonl",
        "molblock",
        {"default", "strict-corpus"},
    ),
    spec(
        "kekulize_flags_golden",
        "kekulize_clear_flags_false.jsonl",
        "kekulize",
        {"default", "strict-corpus"},
    ),
    spec(
        "graph_features",
        "graph_features.jsonl",
        "graph",
        {"default", "strict-corpus"},
    ),
    spec(
        "molecular_descriptors_golden",
        "molecular_descriptors.jsonl",
        "descriptors",
        {"default", "strict-corpus"},
    ),
    spec(
        "delete_substructs_golden",
        "delete_substructs_onlyfrags_chirality.jsonl",
        "substructure",
        {"default", "strict-corpus", "delete-substructs"},
    ),
    spec(
        "smarts_golden",
        "smarts.jsonl",
        "smarts",
        {"smarts"},
        profiles={"smarts_source"},
    ),
    spec(
        "tetrahedral_stereo_geometry",
        "tetrahedral_stereo_geometry.jsonl",
        "stereo",
        {"iterative"},
    ),
    spec(
        "tetrahedral_stereo_geometry",
        "assign_atom_chiral_tags_from_structure.jsonl",
        "stereo",
        {"iterative"},
        extra_inputs=(
            "testdata/stereo/fixtures/assign_atom_chiral_tags_from_structure_cases.json",
        ),
    ),
    spec(
        "ciplabeler_golden",
        "ciplabeler.jsonl",
        "stereo",
        {"ciplabeler", "default", "strict-corpus"},
        deterministic_shards=16,
    ),
    spec(
        "dg_bounds_golden",
        "dg_bounds_matrix.jsonl",
        "distgeom",
        {"default", "strict-corpus"},
    ),
    spec(
        "conformer_generation_golden",
        "conformer_generation.jsonl",
        "conformer",
        {"iterative"},
        extra_inputs=("testdata/conformer/fixtures",),
    ),
    spec(
        "conformer_generation_library_golden",
        "conformer_generation_library.jsonl",
        "conformer",
        {"iterative"},
    ),
    spec(
        "confseq_embed_golden",
        "confseq_embed_template.jsonl",
        "conformer",
        {"iterative"},
    ),
    spec(
        "forcefield_params_golden",
        "forcefield_params.jsonl",
        "forcefield",
        {"iterative"},
    ),
    spec(
        "forcefield_coverage_golden",
        "forcefield_coverage.jsonl",
        "forcefield_coverage",
        {"strict-corpus"},
        generator_dependencies=("_generate_forcefield_params_golden.py",),
    ),
    spec(
        "mmff_builtin_golden",
        "mmff_builtin.jsonl",
        "forcefield",
        {"fixtures"},
        extra_inputs=("testdata/forcefield/fixtures/mmff/rdkit",),
    ),
    spec(
        "smiles_writer_golden",
        "smiles_writer.jsonl",
        "smiles",
        {"default", "strict-corpus"},
    ),
    spec(
        "isomeric_smiles_golden",
        "isomeric_smiles.jsonl",
        "smiles",
        {"default", "strict-corpus"},
    ),
    spec(
        "svg_golden",
        "svg_drawer.jsonl",
        "depiction",
        {"default", "strict-corpus"},
    ),
    spec(
        "prepared_draw_golden",
        "prepared_draw_molecule.jsonl",
        "depiction",
        {"default", "strict-corpus"},
    ),
    spec(
        "morgan_fingerprint_golden",
        "morgan_fingerprint.jsonl",
        "fingerprint",
        {"default", "strict-corpus", "fingerprint"},
    ),
    spec(
        "atom_pair_fingerprint_golden",
        "atom_pair_fingerprint.jsonl",
        "fingerprint",
        {"default", "strict-corpus", "fingerprint", "atom-pair"},
        generator_dependencies=("atom_pair_fingerprint_profile.json",),
    ),
    spec(
        "pattern_fingerprint_golden",
        "pattern_fingerprint.jsonl",
        "fingerprint",
        {"default", "strict-corpus", "fingerprint", "pattern"},
        generator_dependencies=("pattern_fingerprint_profile.json",),
    ),
    spec(
        "tautomer_golden",
        "tautomer.jsonl",
        "tautomer",
        {"tautomer"},
        generator_dependencies=("_tautomer_oracle.py", "tautomer_profile.json"),
        deterministic_shards=16,
        profiles={
            "tautomer_focused",
            "tautomer_pcs_1k",
            "tautomer_pcs_100k",
            "smiles_5000",
        },
    ),
    spec(
        "tautomer_catalog_golden",
        "tautomer_catalog.jsonl",
        "tautomer",
        {"tautomer"},
        extra_inputs=(
            "third_party/rdkit/Code/GraphMol/MolStandardize/TautomerCatalog/tautomerTransforms.in",
            "third_party/rdkit/Code/GraphMol/MolStandardize/TautomerCatalog/tautomerTransforms.v1.in",
        ),
        generator_dependencies=("_tautomer_oracle.py", "tautomer_profile.json"),
        profiles={"tautomer_focused"},
    ),
    spec(
        "maccs_fingerprint_golden",
        "maccs_fingerprint.jsonl",
        "fingerprint",
        {"default", "strict-corpus", "fingerprint"},
    ),
    spec(
        "rdkit_topological_fingerprint_golden",
        "rdkit_topological_fingerprint.jsonl",
        "fingerprint",
        {"default", "strict-corpus", "fingerprint"},
        generator_dependencies=("rdkit_topological_fingerprint_profile.json",),
    ),
    spec(
        "layered_fingerprint_golden",
        "layered_fingerprint.jsonl",
        "fingerprint",
        {"default", "strict-corpus", "fingerprint", "layered"},
        generator_dependencies=("layered_fingerprint_profile.json",),
    ),
    spec(
        "topological_torsion_fingerprint_golden",
        "topological_torsion_fingerprint.jsonl",
        "fingerprint",
        {"default", "strict-corpus", "fingerprint", "topological-torsion"},
        generator_dependencies=("topological_torsion_fingerprint_profile.json",),
    ),
    spec(
        "avalon_fingerprint_golden",
        "avalon_fingerprint.jsonl",
        "fingerprint",
        {"default", "strict-corpus", "fingerprint"},
        generator_dependencies=("avalon_fingerprint_profile.json",),
    ),
    spec(
        "sdf_write_golden",
        "sdf_write.jsonl",
        "sdf",
        {"default", "strict-corpus"},
    ),
    spec(
        "sdf_read_golden",
        "sdf_read.jsonl",
        "sdf",
        {"default", "strict-corpus"},
    ),
    spec(
        "mol2_read_golden",
        "mol2_read.jsonl",
        "mol2",
        {"fixtures"},
        extra_inputs=("testdata/mol2/fixtures/rdkit",),
    ),
    spec(
        "molfile_read_golden",
        "molfile_read.jsonl",
        "molblock",
        {"default", "strict-corpus"},
        extra_inputs=("testdata/rdkit_builtin/fixtures/Code/GraphMol/FileParsers",),
        generator_dependencies=("_generate_sdf_read_golden.py",),
    ),
    spec(
        "xyz_read_golden",
        "xyz_read.jsonl",
        "xyz",
        {"default", "strict-corpus"},
    ),
    spec(
        "builtin_fixture_migration_golden",
        "rdkit_builtin_fixture_migration.jsonl",
        "rdkit_builtin",
        {"fixtures"},
        extra_inputs=("testdata/rdkit_builtin/fixtures",),
    ),
]

SUITE_ALIASES = {
    "small": "default",
    "daily": "default",
    "strict": "strict-corpus",
    "5000": "strict-corpus",
    "inchi": "inchi",
}


def parse_args() -> argparse.Namespace:
    suites = sorted({suite for item in GENERATOR_SPECS for suite in item.suites})
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter used to run generator scripts",
    )
    parser.add_argument(
        "--profile",
        choices=sorted(PROFILE_INPUTS),
        default="smiles_small",
    )
    parser.add_argument(
        "--suite",
        choices=sorted({*suites, *SUITE_ALIASES}),
        default="default",
    )
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--clean", action="store_true")
    parser.add_argument(
        "--only",
        default=None,
        help="run one output stem, generator stem, or domain (for example ciplabeler)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="validate cached outputs and fail instead of regenerating stale data",
    )
    parser.add_argument("--jobs", type=int, default=None)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repository_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def expand_input(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    if path.is_dir():
        return sorted(item for item in path.rglob("*") if item.is_file())
    raise FileNotFoundError(f"required generation input does not exist: {path}")


def file_identities(paths: list[Path], checksums: dict[Path, str]) -> list[dict[str, str]]:
    identities = []
    for path in paths:
        resolved = path.resolve()
        if resolved not in checksums:
            checksums[resolved] = sha256_file(resolved)
        checksum = checksums[resolved]
        identities.append({"path": repository_path(resolved), "sha256": checksum})
    return identities


def reference_identity() -> dict[str, Any]:
    identity = json.loads(REFERENCE_IDENTITY_PATH.read_text(encoding="utf-8"))
    runtime = identity.get("runtime")
    if (
        identity.get("name") != "rdkit"
        or not identity.get("version")
        or not isinstance(runtime, dict)
        or runtime.get("implementation") != "cpython"
        or not runtime.get("version")
    ):
        raise SystemExit(f"invalid RDKit reference identity: {REFERENCE_IDENTITY_PATH}")
    return identity


def rdkit_version(python: str, expected: str) -> str:
    result = subprocess.run(
        [python, "-c", "import rdkit; print(rdkit.__version__)"],
        cwd=REPO_ROOT,
        text=True,
        check=True,
        stdout=subprocess.PIPE,
    )
    actual = result.stdout.strip()
    if actual != expected:
        raise SystemExit(
            f"RDKit version is {actual!r}, expected {expected!r}; "
            "run `uv sync --group dev`"
        )
    return actual


def python_runtime_identity(python: str, expected: dict[str, str]) -> dict[str, str]:
    script = (
        "import json, platform, sys; "
        "print(json.dumps({'implementation': sys.implementation.name, "
        "'version': platform.python_version()}, sort_keys=True))"
    )
    result = subprocess.run(
        [python, "-c", script],
        cwd=REPO_ROOT,
        text=True,
        check=True,
        stdout=subprocess.PIPE,
    )
    actual = json.loads(result.stdout)
    if actual != expected:
        raise SystemExit(
            f"Python reference runtime is {actual!r}, expected {expected!r}; "
            "RDKit Python descriptors can change bits across interpreter versions"
        )
    return actual


def platform_identity() -> dict[str, str]:
    return {
        "system": platform.system(),
        "machine": platform.machine(),
    }


def output_identity(
    item: GeneratorSpec,
    input_path: Path,
    profile: str,
    platform_id: dict[str, str],
    checksums: dict[Path, str],
) -> dict[str, Any]:
    generator_paths = [
        TOOL_DIR / "generate_all.py",
        SCHEMA_PATH,
        REFERENCE_IDENTITY_PATH,
        TOOL_DIR / item.script,
    ]
    generator_paths.extend(TOOL_DIR / name for name in item.generator_dependencies)
    input_paths = expand_input(input_path)
    if item.output == "ciplabeler.jsonl" and input_path.suffix.lower() == ".json":
        document = json.loads(input_path.read_text(encoding="utf-8"))
        for case in document.get("cases", []):
            if isinstance(case, dict) and isinstance(case.get("input_file"), str):
                input_paths.extend(expand_input(REPO_ROOT / case["input_file"]))
    for relative in item.extra_inputs:
        input_paths.extend(expand_input(REPO_ROOT / relative))
    arguments = ["--input", repository_path(input_path), "--output", item.output]
    if item.deterministic_shards is not None:
        arguments.extend(["--shards", str(item.deterministic_shards)])
    return {
        "path": item.output,
        "output_schema_version": 1,
        "generator": file_identities(generator_paths, checksums),
        "inputs": file_identities(input_paths, checksums),
        "options": {
            "profile": profile,
            "arguments": arguments,
        },
        "platform": platform_id,
    }


def expected_directory(domain: str, profile: str) -> Path:
    return REPO_ROOT / "testdata" / domain / "expected" / "rdkit" / profile


def manifest_identity(
    domain: str,
    profile: str,
    version: str,
    source_revision: str,
    runtime: dict[str, str],
    input_path: Path,
    checksums: dict[Path, str],
) -> dict[str, Any]:
    input_resolved = input_path.resolve()
    return {
        "schema_version": 1,
        "family": "rdkit",
        "domain": domain,
        "profile": profile,
        "reference_implementation": {
            "name": "rdkit",
            "version": version,
            "source_revision": source_revision,
        },
        "reference_runtime": runtime,
        "input": file_identities([input_resolved], checksums)[0],
    }


def validate_cached_domain(
    target: Path,
    top_identity: dict[str, Any],
    identities: list[dict[str, Any]],
) -> tuple[bool, str]:
    manifest_path = target / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return False, f"cannot read manifest: {error}"
    for key, expected in top_identity.items():
        if manifest.get(key) != expected:
            return False, f"manifest field {key!r} does not match"
    by_path = {entry.get("path"): entry for entry in manifest.get("outputs", [])}
    for identity in identities:
        entry = by_path.get(identity["path"])
        if entry is None:
            return False, f"manifest does not declare {identity['path']}"
        for key, expected in identity.items():
            if entry.get(key) != expected:
                return False, f"{identity['path']} identity field {key!r} does not match"
        output = target / identity["path"]
        try:
            actual_checksum, actual_records = validate_jsonl_output(
                output, identity["path"]
            )
        except (OSError, ValueError) as error:
            return False, f"cannot validate {output}: {error}"
        if entry.get("sha256") != actual_checksum:
            return False, f"{identity['path']} checksum does not match"
        if entry.get("records") != actual_records:
            return False, f"{identity['path']} record count does not match"
    return True, "exact identity and output checksums match"


def run_generator(
    python: str,
    item: GeneratorSpec,
    input_path: Path,
    output_path: Path,
    jobs: int,
) -> subprocess.CompletedProcess[str]:
    command = [
        python,
        str(TOOL_DIR / item.script),
        "--input",
        str(input_path),
        "--output",
        str(output_path),
    ]
    if item.deterministic_shards is not None:
        command.extend(
            ["--shards", str(item.deterministic_shards), "--jobs", str(jobs)]
        )
    return subprocess.run(
        command,
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def print_completed(item: GeneratorSpec, result: subprocess.CompletedProcess[str]) -> None:
    print(f"[expected] {item.script} -> {item.domain}/{item.output}")
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")
    if result.stderr:
        print(result.stderr, end="" if result.stderr.endswith("\n") else "\n", file=sys.stderr)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode,
            result.args,
            output=result.stdout,
            stderr=result.stderr,
        )


def publish_directory(staging: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    backup = target.with_name(f".{target.name}.backup-{uuid.uuid4().hex}")
    if target.exists():
        target.rename(backup)
    try:
        staging.rename(target)
    except BaseException:
        if backup.exists() and not target.exists():
            backup.rename(target)
        raise
    if backup.exists():
        shutil.rmtree(backup)


def main() -> None:
    args = parse_args()
    selected_suite = SUITE_ALIASES.get(args.suite, args.suite)
    if args.only is None:
        selected = [
            item
            for item in GENERATOR_SPECS
            if selected_suite in item.suites
            and (item.profiles is None or args.profile in item.profiles)
        ]
    else:
        selected = [
            item
            for item in GENERATOR_SPECS
            if (item.profiles is None or args.profile in item.profiles)
            and args.only
            in {
                item.domain,
                Path(item.output).stem,
                Path(item.script).stem.removeprefix("_generate_").removesuffix(
                    "_golden"
                ),
            }
        ]
    if not selected:
        selector = f"--only {args.only!r}" if args.only else f"--suite {args.suite!r}"
        raise SystemExit(f"{selector} did not select any generators")
    jobs = args.jobs or min(4, os.cpu_count() or 1, len(selected))
    if jobs < 1:
        raise SystemExit("--jobs must be >= 1")

    input_path = (args.input or PROFILE_INPUTS[args.profile]).resolve()
    if not input_path.is_file():
        raise SystemExit(f"input corpus does not exist: {input_path}")
    missing_schemas = sorted({item.output for item in GENERATOR_SPECS} - set(SCHEMAS))
    extra_schemas = sorted(set(SCHEMAS) - {item.output for item in GENERATOR_SPECS})
    if missing_schemas or extra_schemas:
        raise SystemExit(
            f"RDKit expected-data schema registry mismatch: "
            f"missing={missing_schemas}, extra={extra_schemas}"
        )
    reference = reference_identity()
    runtime = python_runtime_identity(args.python, reference["runtime"])
    version = rdkit_version(args.python, reference["version"])
    platform_id = platform_identity()
    checksums: dict[Path, str] = {}

    by_domain: dict[str, list[GeneratorSpec]] = {}
    for item in selected:
        by_domain.setdefault(item.domain, []).append(item)

    stale: dict[str, tuple[list[GeneratorSpec], dict[str, Any], list[dict[str, Any]]]] = {}
    for domain, items in by_domain.items():
        top_identity = manifest_identity(
            domain,
            args.profile,
            version,
            reference["source_revision"],
            runtime,
            input_path,
            checksums,
        )
        identities = [
            output_identity(item, input_path, args.profile, platform_id, checksums)
            for item in items
        ]
        target = expected_directory(domain, args.profile)
        valid, reason = validate_cached_domain(target, top_identity, identities)
        if valid and not args.clean:
            print(f"[expected] reuse {domain}/{args.profile}: {reason}")
        else:
            cause = "--clean requested" if args.clean else reason
            print(f"[expected] regenerate {domain}/{args.profile}: {cause}")
            stale[domain] = (items, top_identity, identities)

    if not stale:
        print(
            f"Validated {len(selected)} RDKit expected outputs for "
            f"profile {args.profile} suite {selected_suite}"
        )
        return

    if args.check:
        stale_domains = ", ".join(sorted(stale))
        raise SystemExit(f"generated RDKit expected data is stale: {stale_domains}")

    staging_dirs: dict[str, Path] = {}
    preserved_outputs: dict[str, list[dict[str, Any]]] = {}
    try:
        for domain, (items, top_identity, _) in stale.items():
            target = expected_directory(domain, args.profile)
            parent = target.parent
            parent.mkdir(parents=True, exist_ok=True)
            staging = Path(
                tempfile.mkdtemp(prefix=f".{args.profile}.prepare-", dir=parent)
            )
            staging_dirs[domain] = staging
            preserved_outputs[domain] = []
            if target.is_dir():
                try:
                    previous_manifest = json.loads(
                        (target / "manifest.json").read_text(encoding="utf-8")
                    )
                except (OSError, json.JSONDecodeError):
                    previous_manifest = {}
                if all(
                    previous_manifest.get(key) == value
                    for key, value in top_identity.items()
                ):
                    selected_paths = {item.output for item in items}
                    for entry in previous_manifest.get("outputs", []):
                        path = entry.get("path")
                        source = target / path if isinstance(path, str) else None
                        if (
                            isinstance(path, str)
                            and path not in selected_paths
                            and source is not None
                            and source.is_file()
                        ):
                            shutil.copy2(source, staging / path)
                            preserved_outputs[domain].append(entry)

        work = []
        for domain, (items, _, _) in stale.items():
            for item in items:
                work.append((item, staging_dirs[domain] / item.output))

        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
            futures = {
                executor.submit(
                    run_generator, args.python, item, input_path, output, jobs
                ): item
                for item, output in work
            }
            for future in concurrent.futures.as_completed(futures):
                item = futures[future]
                print_completed(item, future.result())

        for domain, (_, top_identity, identities) in stale.items():
            staging = staging_dirs[domain]
            outputs = list(preserved_outputs[domain])
            for identity in identities:
                output = staging / identity["path"]
                entry = dict(identity)
                entry["sha256"], entry["records"] = validate_jsonl_output(
                    output, identity["path"]
                )
                outputs.append(entry)
            manifest = {**top_identity, "outputs": outputs}
            (staging / "manifest.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            valid, reason = validate_cached_domain(staging, top_identity, identities)
            if not valid:
                raise RuntimeError(f"generated {domain} family failed validation: {reason}")

        for domain in sorted(stale):
            publish_directory(staging_dirs.pop(domain), expected_directory(domain, args.profile))
            print(f"[expected] published {domain}/{args.profile}")
    finally:
        for staging in staging_dirs.values():
            if staging.exists():
                shutil.rmtree(staging)

    print(
        f"Prepared {len(selected)} RDKit expected outputs for "
        f"profile {args.profile} suite {selected_suite}"
    )


if __name__ == "__main__":
    main()
