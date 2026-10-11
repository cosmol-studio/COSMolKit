"""Bounded native preparation. Rust owns task/profile selection and cache keys."""
import concurrent.futures as futures
from collections import deque
import hashlib
import itertools
import json
import multiprocessing
import os
from pathlib import Path
import struct
import sys
import tempfile
import time
import zlib

# Prevent native numerical libraries from multiplying the chosen process count.
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdCIPLabeler
from fingerprints import UInt, morgan, pair_fingerprint, path_fingerprint
from isolate import isolated
from descriptors import descriptor
from chemistry import TASKS as CHEMISTRY_TASKS, parse as parse_chemistry, chemistry
from search_tautomer import search,tautomer
from geometry import TASKS as GEOMETRY_TASKS, geometry
import bio
import special

RDLogger.DisableLog("rdApp.*")
MAX_FRAME = 64 * 1024 * 1024
MAX_RECORD = MAX_FRAME - 4


def encode(value):
    if value is None:
        return b"\0"
    if isinstance(value, bool):
        return bytes((1, int(value)))
    if isinstance(value, UInt):
        return b"\3" + struct.pack("<Q", value)
    if isinstance(value, int):
        return b"\2" + struct.pack("<q", value)
    if isinstance(value, float):
        return b"\4" + struct.pack("<d", value)
    if isinstance(value, str):
        raw = value.encode("utf-8")
        return b"\5" + struct.pack("<I", len(raw)) + raw
    if isinstance(value, bytes):
        return b"\6" + struct.pack("<I", len(value)) + value
    if isinstance(value, (list, tuple)):
        return b"\7" + struct.pack("<I", len(value)) + b"".join(encode(item) for item in value)
    raise TypeError(f"unregistered reference output type {type(value)}")


def native(mol, task, profile):
    if task in GEOMETRY_TASKS: return geometry(mol,task,profile)
    if task=="substructure": return search(mol,profile)
    if task=="tautomer": return tautomer(mol,profile)
    if task in CHEMISTRY_TASKS:
        return chemistry(mol,task,profile)
    if task == "fingerprint_morgan":
        return morgan(mol, profile)
    if task in ("fingerprint_atom_pair", "fingerprint_torsion"):
        return pair_fingerprint(mol, task, profile)
    if task.startswith("fingerprint_"):
        return path_fingerprint(mol, task, profile)
    if task == "smiles_write":
        params = Chem.SmilesWriteParams()
        for public, native in (
            ("isomeric_smiles", "doIsomericSmiles"), ("kekule", "doKekule"),
            ("canonical", "canonical"), ("clean_stereo", "cleanStereo"),
            ("all_bonds_explicit", "allBondsExplicit"), ("all_hydrogens_explicit", "allHsExplicit"),
            ("include_dative_bonds", "includeDativeBonds"), ("ignore_atom_map_numbers", "ignoreAtomMapNumbers"),
        ):
            setattr(params, native, profile[public])
        params.rootedAtAtom = {"default": -1, "first": 0, "last": max(0, mol.GetNumAtoms() - 1)}[profile["root"]]
        return Chem.MolToSmiles(mol, params).encode("utf-8")
    if task == "valence":
        if profile["recalculate"]:
            mol.UpdatePropertyCache(strict=True)
        return [[atom.GetDegree(), atom.GetValence(Chem.ValenceType.EXPLICIT),
                 atom.GetNumImplicitHs(), atom.GetTotalNumHs(), atom.GetTotalValence()]
                for atom in mol.GetAtoms()]
    if task == "cip":
        labeled = Chem.Mol(mol)
        selection = profile["selection"]
        atoms = None if selection == "all" else []
        bonds = None if selection == "all" else []
        if selection == "first_tagged_atom":
            atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED][:1]
        if selection == "first_stereo_bond":
            bonds = [b.GetIdx() for b in mol.GetBonds() if b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)][:1]
        rdCIPLabeler.AssignCIPLabels(labeled, atoms, bonds, profile["max_recursive_iterations"])
        return [[item.GetProp("_CIPCode") if item.HasProp("_CIPCode") else None
                 for item in items] for items in (labeled.GetAtoms(), labeled.GetBonds())]
    return descriptor(mol,task,profile)


def iter_outcomes(task, profiles, smiles):
    for text in smiles:
        if task in special.TASKS:
            try:
                result = special.observe(task, text)
                # A complete MCS timeout is observable but not a successful comparison.
                status = 3 if task.startswith("mcs_") and result == [False] else 0
                value = b"\1" + encode("parse_rejected") if result is None else bytes((status,)) + encode(result)
            except Exception as error:
                value = b"\5" + encode(f"{type(error).__name__}: {error}")
            yield [value for _ in profiles]
            continue
        if task == "bio_pdb_output":
            structure = bio.parse(text)
            yield [b"\0" + encode(bio.pdb(structure, profile)) for profile in profiles]
            continue
        mol = None if task in CHEMISTRY_TASKS else Chem.MolFromSmiles(text)
        retained = None
        if task not in CHEMISTRY_TASKS and any(p.get("remove_hs") is False for p in profiles):
            params=Chem.SmilesParserParams()
            params.removeHs=False
            retained=Chem.MolFromSmiles(text,params)
        row = []
        for profile in profiles:
            if task in CHEMISTRY_TASKS:
                mol = parse_chemistry(text,task,profile)
            selected=retained if task not in CHEMISTRY_TASKS and profile.get("remove_hs") is False else mol
            if selected is None:
                output = b"\1" + encode("parse_rejected")
            else:
                try:
                    output = b"\0" + encode(native(selected, task, profile))
                except Exception as error:
                    # Never classify arbitrary exceptions as matching chemistry errors.
                    output = b"\5" + encode(f"{type(error).__name__}: {error}")
            row.append(output)
        yield row


def evaluate_rows(task, profiles, smiles):
    return list(iter_outcomes(task, profiles, smiles))


def compress_frame(raw, count):
    compressed = zlib.compress(raw, 1)
    if len(compressed) > MAX_FRAME:
        raise RuntimeError("compressed reference exceeds frame budget")
    return (struct.pack("<III", len(compressed), len(raw), count)
            + hashlib.sha256(compressed).digest() + compressed)


def calculate(job):
    task, profiles, smiles, folder = job
    outcomes = (isolated(evaluate_rows, task, profiles, smiles) if task == "fingerprint_layered"
                else iter_outcomes(task, profiles, smiles))
    raw = bytearray()
    frames, count = 0, 0
    # Keep pending results on disk, not in IPC queues or parent RAM. Each
    # worker retains only its current record/frame regardless of job length.
    with tempfile.NamedTemporaryFile(dir=folder, prefix="worker-", delete=False) as output:
        for index, values in enumerate(outcomes):
            row = bytearray(struct.pack("<I", len(profiles)))
            for value in values:
                row.extend(struct.pack("<I", len(value)))
                row.extend(value)
            if len(row) > MAX_RECORD:
                raise RuntimeError(f"{task} input {index}: reference row {len(row)} bytes exceeds {MAX_RECORD}; input={smiles[index]!r}")
            if len(raw) + len(row) + 4 > MAX_FRAME:
                output.write(compress_frame(raw,count))
                frames += 1
                raw.clear()
                count = 0
            raw.extend(struct.pack("<I", len(row)))
            raw.extend(row)
            count += 1
        if count:
            output.write(compress_frame(raw,count))
            frames += 1
    return frames, output.name


def inputs(request):
    with Path(request["corpus"]).open(encoding="utf-8") as stream:
        rows = (line.strip() for line in stream)
        structured = Path(request["corpus"]).suffix == ".jsonl"
        rows = ((line if structured else line.split()[0]) for line in rows if line and not line.startswith("#"))
        end = None if request["limit"] is None else request["offset"] + request["limit"]
        yield from itertools.islice(rows, request["offset"], end)


def generate(pool, request, task):
    request = dict(request, **{key:task[key] for key in ("corpus", "offset", "limit")})
    identity, folder = task["identity"], Path(task["folder"])
    name, profiles = identity["task"], identity["profiles"]
    folder.mkdir(parents=True, exist_ok=True)
    rows = iter(inputs(request))
    pending = deque()
    # Small selections still feed every worker; large selections amortize IPC.
    chunk = min(64, max(1, (identity["records"] + request["threads"] * 4 - 1)
                           // (request["threads"] * 4)))
    input_digest = hashlib.sha256(struct.pack("<Q", request["offset"]))
    # Results are bounded by jobs in flight, not by corpus length.
    def submit():
        batch = list(itertools.islice(rows, chunk))
        if batch:
            for text in batch:
                encoded = text.encode("utf-8")
                input_digest.update(struct.pack("<Q", len(encoded)))
                input_digest.update(encoded)
            pending.append((len(batch), pool.submit(calculate, (name, profiles, batch, str(folder)))))
    for _ in range(request["threads"] * 2):
        submit()
    started, updated = time.monotonic(), 0.0
    frames, count, byte_count = 0, 0, 8
    digest = hashlib.sha256(b"CKREF001")
    temporary = tempfile.NamedTemporaryFile(dir=folder, prefix="incomplete-", delete=False)
    try:
        with temporary as output:
            output.write(b"CKREF001")
            while pending:
                size, future = pending.popleft()
                frame_count, worker_file = future.result()
                with open(worker_file, "rb") as stream:
                    while frame := stream.read(1024 * 1024):
                        output.write(frame)
                        digest.update(frame)
                        byte_count += len(frame)
                os.unlink(worker_file)
                frames += frame_count
                count += size
                submit()
                now = time.monotonic()
                if now - updated >= 0.2 or count == identity["records"]:
                    width = 24
                    bars = width * count // identity["records"]
                    print(f"\r{name} [{'=' * bars}{' ' * (width - bars)}] {count}/{identity['records']} ({now - started:.1f}s)", end="", file=sys.stderr, flush=True)
                    updated = now
            output.flush()
            os.fsync(output.fileno())
        print(file=sys.stderr)
        if count != identity["records"] or input_digest.hexdigest() != identity["input_sha256"]:
            raise RuntimeError("corpus changed during preparation")
        # Publish data first; a manifest only exists for a complete generation.
        os.replace(temporary.name, folder / "reference.bin")
        manifest = {"format": 1, "identity": identity, "frames": frames,
                    "bytes": byte_count, "sha256": digest.hexdigest()}
        with tempfile.NamedTemporaryFile(mode="w", dir=folder, prefix="manifest-", delete=False) as output:
            json.dump(manifest, output, separators=(",", ":"), ensure_ascii=False)
            output.flush()
            os.fsync(output.fileno())
        os.replace(output.name, folder / "manifest.json")
    except BaseException:
        # Leave the partial file as evidence; it cannot be mistaken for a cache hit.
        print(f"\nIncomplete reference retained: {temporary.name}", file=sys.stderr)
        raise


def main():
    request = json.load(sys.stdin)
    if rdBase.rdkitVersion != request["rdkit"]:
        raise RuntimeError(f"need RDKit {request['rdkit']}, found {rdBase.rdkitVersion}")
    context = multiprocessing.get_context("fork" if sys.platform.startswith("linux") else "spawn")
    with futures.ProcessPoolExecutor(max_workers=request["threads"], mp_context=context) as pool:
        for index, task in enumerate(request["tasks"], 1):
            print(f"Task {index}/{len(request['tasks'])}: {task['identity']['task']} ({request['threads']} processes)", file=sys.stderr, flush=True)
            generate(pool, request, task)


if __name__ == "__main__":
    main()
