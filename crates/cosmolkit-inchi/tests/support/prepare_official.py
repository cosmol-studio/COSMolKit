"""Prepare local output oracles using the pinned, unmodified official library."""
import argparse
import ctypes
import hashlib
import json
import random
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
UPSTREAM = ROOT / "third_party/InChI"
TESTS = UPSTREAM / "INCHI-1-TEST"
sys.path.insert(0, str(TESTS / "tests/test_library/inchi_tests"))
sys.path.insert(0, str(TESTS / "src"))
from inchi_api import make_inchi_from_molfile_text, get_inchi_key_from_inchi
from sdf_pipeline.utils import read_records_from_gzipped_sdf, permute_molblock


def add_rdkit_reference(reference):
    """Add the independent toolkit pipeline without rerunning official InChI."""
    from rdkit import Chem, RDLogger, rdBase
    from rdkit.Chem import inchi
    if rdBase.rdkitVersion != "2026.09.1":
        raise RuntimeError(f"expected pinned RDKit 2026.09.1, got {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")
    digest = hashlib.sha256(Path(__file__).read_bytes())
    with reference.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    digest.update(rdBase.rdkitVersion.encode())
    output = reference.with_name(f"reference.rdkit-{digest.hexdigest()[:20]}.jsonl")
    if output.exists():
        return output
    pending = output.with_suffix(".partial.jsonl")
    with reference.open() as source, pending.open("w") as destination:
        for line in source:
            row = json.loads(line)
            if "molblock" in row:
                mol = Chem.MolFromMolBlock(row["molblock"], sanitize=False,
                                          removeHs=False, strictParsing=False)
                if mol is None:
                    result = {"stage": "sdf_rejected"}
                else:
                    try:
                        text = inchi.MolToInchi(mol)
                        result = {"stage": "output", "inchi": text,
                                  "key": inchi.InchiToInchiKey(text) if text else ""}
                    except Chem.AtomValenceException as error:
                        result = {"stage": "inchi_error", "kind": "valence",
                                  "detail": str(error)}
                    except Chem.KekulizeException as error:
                        # Preserve the implicated atoms, not incidental prose.
                        atoms = str(error).split("Unkekulized atoms:", 1)[1]
                        result = {"stage": "inchi_error", "kind": "kekulize",
                                  "atoms": [int(n) for n in re.findall(r"\d+", atoms)]}
                    # Unknown exceptions/native crashes abort preparation; they
                    # are not an oracle permitting arbitrary CK errors.
                row["rdkit"] = result
                row["rdkit_version"] = rdBase.rdkitVersion
            destination.write(json.dumps(row) + "\n")
    pending.rename(output)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=0, help="records per dataset; 0 means all")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--reuse-official", type=Path,
                        help="add RDKit outcomes to an existing official reference; no C rebuild")
    args = parser.parse_args()
    if args.limit < 0 or args.jobs < 1:
        parser.error("limit must be nonnegative and jobs positive")
    if args.reuse_official:
        print(add_rdkit_reference(args.reuse_official.resolve()))
        return
    import rdkit
    digest = hashlib.sha256()
    for path in sorted(UPSTREAM.rglob("*")):
        if path.is_file() and ".git" not in path.parts:
            digest.update(str(path.relative_to(UPSTREAM)).encode())
            digest.update(path.read_bytes())
    digest.update(Path(__file__).read_bytes())
    digest.update(f"{rdkit.__version__}:{args.limit}".encode())
    folder = ROOT / "target/inchi-official-special" / digest.hexdigest()[:20]
    output = folder / "reference.jsonl"
    if output.exists():
        print(add_rdkit_reference(output))
        return
    folder.mkdir(parents=True, exist_ok=True)
    build = folder / "build"
    with (folder / "build.log").open("w") as log:
        subprocess.run(["cmake", "-S", str(UPSTREAM / "INCHI-1-SRC/INCHI_API/libinchi/src"),
                        "-B", str(build)], stdout=log, stderr=log, check=True)
        subprocess.run(["cmake", "--build", str(build), "--parallel", str(args.jobs)],
                       stdout=log, stderr=log, check=True)
    libraries = [p for p in build.rglob("*") if p.is_file()
                 and p.suffix in (".so", ".dylib", ".dll") and "inchi" in p.name.lower()]
    if len(libraries) != 1:
        raise RuntimeError(f"expected one official library: {libraries}")
    library = ctypes.CDLL(str(libraries[0]))
    pending = folder / "reference.partial.jsonl"
    with pending.open("w") as stream:
        for dataset in ("inchi", "mcule"):
            source = TESTS / f"tests/test_library/data/ci/{dataset}.sdf.gz"
            for index, block in enumerate(read_records_from_gzipped_sdf(source)):
                if args.limit and index >= args.limit:
                    break
                # Official invariance_consumer resets the seed for each molecule.
                random.seed(42)
                for variant in range(11):
                    text = block if variant == 0 else permute_molblock(block)
                    if text is None:
                        row = dict(dataset=dataset, index=index, variant=variant, preparation_error="official permutation returned None")
                    else:
                        status, inchi, _, _, _ = make_inchi_from_molfile_text(library, text)
                        key_status, key = get_inchi_key_from_inchi(library, inchi) if inchi else (-1, "")
                        row = dict(dataset=dataset, index=index, variant=variant,
                                   molblock=text, status=status, inchi=inchi,
                                   key_status=key_status, key=key)
                    stream.write(json.dumps(row) + "\n")
                    stream.flush()
    pending.rename(output)
    print(add_rdkit_reference(output))


if __name__ == "__main__":
    main()
