"""Preparation only: pinned reference calls, ordered output, bounded workers."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "tools/oracles/rdkit"), str(ROOT / "tools/testdata/rdkit")]


class Progress:
    def __init__(self, task):
        self.task = task
        self.started = time.monotonic()
        self.updated = 0

    def __call__(self, done, total):
        now = time.monotonic()
        if done not in (0, total) and now - self.updated < 0.2:
            return
        self.updated = now
        filled = 24 * done // max(1, total)
        bar = "=" * filled + " " * (24 - filled)
        prefix = "\r" if sys.stderr.isatty() else ""
        end = "" if sys.stderr.isatty() and done != total else "\n"
        print(f"{prefix}{self.task} [{bar}] {done}/{total} cases "
              f"({100 * done / max(1, total):5.1f}%) {now - self.started:.1f}s",
              file=sys.stderr, end=end, flush=True)


def collect(pool, function, rows, progress):
    progress(0, len(rows))
    results = [None] * len(rows)
    pending = {pool.submit(function, row): index for index, row in enumerate(rows)}
    for done, future in enumerate(as_completed(pending), 1):
        results[pending[future]] = future.result()
        progress(done, len(rows))
    return results


def parallel(function, rows, threads, progress):
    with ProcessPoolExecutor(max_workers=min(threads, max(1, len(rows)))) as pool:
        return collect(pool, function, rows, progress)


def descriptor_case(case):
    from _generate_molecular_descriptors_golden import build_record
    return build_record(case["smiles"])


def tautomer_case(payload):
    from _tautomer_oracle import parse_molecule, enumerate_branch, canonicalize_branch
    original = payload
    profile = original["Molecular"]["profile"]
    operation, parameters = next(iter(profile.items()))
    options = parameters["parameters"]
    molecule = parse_molecule({"smiles":original["Molecular"]["case"]["smiles"]})
    if molecule is None:
        observation = {"Error":{"stage":"Parse","detail":"MolFromSmiles returned None"}}
    else:
        branch = {**options, "catalog":"v1" if options["catalog"] == "v1" else "default"}
        result = (enumerate_branch if operation == "TautomerEnumeration" else canonicalize_branch)(molecule, branch)
        if result["ok"]:
            observation = {operation:{key:value for key,value in result.items() if key not in ("ok", "error")}}
        else:
            observation = {"Error":{"stage":"Operation","detail":json.dumps(result["error"],sort_keys=True)}}
    return {"input":original,"output":{"Molecular":observation}}


def batch_case(case):
    from rdkit import Chem
    molecule = Chem.MolFromSmiles(case["smiles"])
    return None if molecule is None else Chem.MolToSmiles(molecule)


def mmff_case(payload):
    import struct
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from fingerprint_values_pilot import prepare_forcefield_geometry
    case, profiles = payload
    bits = lambda value: struct.unpack(">Q", struct.pack(">d", value))[0]
    rows = []
    for profile in profiles:
        name, options = next(iter(profile.items()))
        row = {"case": case, "profile": profile, "preparation": None}
        stage = "Parse"
        try:
            if name == "Coverage":
                mol = Chem.MolFromSmiles(case["smiles"])
                if mol is None:
                    raise ValueError("RDKit MolFromSmiles returned None")
                stage = "Preparation"
                if options["add_hydrogens"]:
                    mol = Chem.AddHs(mol)
                stage = "Operation"
                output = {"Coverage": AllChem.MMFFHasAllMoleculeParams(mol)}
            else:
                # Reuse input preparation only; no UFF optimization is performed.
                common = {"case": case, "profile": {name: {**options, "add_hydrogens": True}}, "preparation": None}
                row["preparation"] = prepare_forcefield_geometry(common, label="MMFF")
                if "Rejected" in row["preparation"]:
                    rows.append({"input": {"Mmff": row}, "output": {"Mmff": {"Error": row["preparation"]["Rejected"]}}})
                    continue
                stage = "Preparation"
                geometry = row["preparation"]["Ready"]
                mol = Chem.MolFromMolBlock(geometry["molblock"], sanitize=True, removeHs=False)
                if mol is None or mol.GetNumAtoms() != geometry["atom_count"]:
                    raise ValueError("invalid prepared MMFF MolBlock")
                if name == "ConformerOptimization":
                    mol.RemoveAllConformers()
                    for prepared in geometry["coordinate_rows"]:
                        conformer = Chem.Conformer(mol.GetNumAtoms())
                        conformer.SetId(prepared["conformer_id"])
                        for atom, xyz in enumerate(prepared["xyz_bits"]):
                            conformer.SetAtomPosition(atom, [struct.unpack(">d", struct.pack(">Q", value))[0] for value in xyz])
                        mol.AddConformer(conformer, assignId=False)
                    stage = "Operation"
                    results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1,
                        maxIters=options["max_iterations"], mmffVariant=options["variant"],
                        nonBondedThresh=100., ignoreInterfragInteractions=True)
                    conformers = list(mol.GetConformers())
                    if len(results) != len(conformers):
                        raise RuntimeError("RDKit MMFF conformer result count changed")
                    output = {"OptimizedConformers": {"conformers": [
                        {"conformer_id": conf.GetId(), "status": status,
                         "energy_bits": bits(energy), "xyz_bits": [[bits(value) for value in conf.GetAtomPosition(i)] for i in range(mol.GetNumAtoms())]}
                        for conf, (status, energy) in zip(conformers, results, strict=True)]}}
                else:
                    stage = "Operation"
                    status = AllChem.MMFFOptimizeMolecule(mol, maxIters=options["max_iterations"],
                        mmffVariant=options["variant"], nonBondedThresh=100.,
                        confId=-1, ignoreInterfragInteractions=True)
                    energy = None
                    if status != -1:
                        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=options["variant"])
                        field = AllChem.MMFFGetMoleculeForceField(mol, props,
                            nonBondedThresh=100., confId=-1, ignoreInterfragInteractions=True)
                        if field is None:
                            raise RuntimeError("MMFF evaluation returned None after optimization")
                        energy = bits(field.CalcEnergy())
                    conf = mol.GetConformer()
                    output = {"Optimized": {"status": status, "energy_bits": energy,
                        "xyz_bits": [[bits(value) for value in conf.GetAtomPosition(i)] for i in range(mol.GetNumAtoms())]}}
        except (ValueError, RuntimeError) as error:
            output = {"Error": {"stage": stage, "detail": f"{type(error).__name__}: {error}"}}
        rows.append({"input": {"Mmff": row}, "output": {"Mmff": output}})
    return rows


def writer_case(payload):
    from rdkit import Chem
    from _generate_smiles_writer_golden import branch_result
    case, profiles = payload
    molecule = Chem.MolFromSmiles(case["smiles"])
    if molecule is None:
        return [{"Error": "Parse"} for _ in profiles]
    outcomes = []
    for profile in profiles:
        params = {**profile, "rooted_at_atom": None if profile["rooted_at_atom"] == "none" else profile["rooted_at_atom"]}
        result = branch_result(molecule, {"params": params})
        outcomes.append({"Smiles": result["smiles"]} if result["ok"] else {"Error": "Write"})
    return outcomes


def writer_rows(cases, profiles, threads, progress):
    from _generate_smiles_writer_golden import iter_branches
    existing = [branch["params"] for branch in iter_branches()]
    existing = [{**p, "rooted_at_atom": p["rooted_at_atom"] or "none"} for p in existing]
    if profiles != existing:
        raise ValueError("writer profiles must equal the existing ordered 768-branch matrix")
    return parallel(writer_case, [(case, profiles) for case in cases], threads, progress)


def structure_case(payload):
    from _generate_tetrahedral_stereo_geometry import deep_merge
    merged = deep_merge(payload["defaults"], payload["case"])
    environment = os.environ.copy()
    config = merged["environment"]
    name = "RDK_ENABLE_NONTETRAHEDRAL_STEREO"
    if config["mode"] == "unset":
        environment.pop(name, None)
    elif config["mode"] == "set":
        environment[name] = config["value"]
    else:
        raise ValueError("unknown structure-tag environment")
    result = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--structure-case"],
                            input=json.dumps(payload), text=True, capture_output=True, env=environment, check=True)
    return json.loads(result.stdout)


def bio_case(original):
    import gemmi
    case = original["BioPdbOutput"]["case"]
    profile = original["BioPdbOutput"]["profile"]
    with tempfile.TemporaryDirectory(prefix="fixed-gemmi-") as folder:
        source = Path(folder) / ("input." + case["format"])
        source.write_text(case["text"], encoding="utf-8")
        # The previous C++ reader did not merge chain parts. Python defaults to merging.
        structure = gemmi.read_structure(str(source), merge_chain_parts=False)
    options = gemmi.PdbWriteOptions()
    for key in ("ter_records", "numbered_ter", "ter_ignores_type", "preserve_serial", "end_record"):
        setattr(options, key, profile[key])
    options.minimal_file = True
    options.seqres_records = False
    options.ssbond_records = False
    options.link_records = False
    options.cispep_records = False
    records = ("ATOM  ", "HETATM", "ANISOU", "TER   ", "MODEL ", "ENDMDL", "END   ")
    text = "".join(line + "\n" for line in structure.make_pdb_string(options).split("\n")
                   if line[:6] in records)
    return {"input":original,"output":{"BioPdbOutput":{"text":text,"error":None}}}


def generate(request):
    threads = request["threads"]
    if not isinstance(threads,int) or threads < 1:
        raise ValueError("threads must be positive")
    kind = request["kind"]
    progress = Progress(request["task"])
    if kind == "corpus" and request["generator"] in ("generate_bio_pdb_output_pdb", "generate_bio_pdb_output_cif"):
        import gemmi
        pin = json.loads((PACKAGE / "testdata/reference/gemmi.json").read_text())
        if gemmi.__version__ != pin["version"]:
            raise RuntimeError(f"Gemmi version {gemmi.__version__} != {pin['version']}")
        return parallel(bio_case, request["input"], threads, progress)
    from rdkit import rdBase
    pin = json.loads((PACKAGE / "testdata/reference/rdkit.json").read_text())
    if rdBase.rdkitVersion != pin["version"]:
        raise RuntimeError(f"RDKit version {rdBase.rdkitVersion} != {pin['version']}")
    if kind == "descriptors":
        return parallel(descriptor_case, request["corpus"], threads, progress)
    if kind == "batch_smiles":
        recipes = request["input"]
        values = parallel(batch_case, recipes[0]["cases"], threads, progress)
        output = {"valid_mask":[value is not None for value in values],"smiles":values,"error_indices":[i for i,value in enumerate(values) if value is None]}
        return [{"input":recipe,"output":output} for recipe in recipes]
    if kind == "structure_tags":
        from _generate_tetrahedral_stereo_geometry import octahedral_case
        fixture = request["input"]
        cases = fixture["cases"] + [octahedral_case(case) for case in fixture["octahedral_switch_cases"]]
        with ThreadPoolExecutor(max_workers=threads) as pool:
            return collect(pool, structure_case, [{"defaults":fixture["defaults"],"case":case} for case in cases], progress)
    if kind == "tautomer_long_conjugated":
        from _tautomer_oracle import build_record
        fixture = request["input"]
        profile = {"branches":fixture["branches"]}
        rows = []
        progress(0, len(fixture["cases"]))
        for case in fixture["cases"]:
            rows.append(build_record(case, profile, ["default","v1"]))
            progress(len(rows), len(fixture["cases"]))
        return rows
    if kind != "corpus":
        raise ValueError(f"unknown recipe {kind}")
    generator = request["generator"]
    if generator == "generate_fingerprint":
        from fingerprints import fingerprint_case
        return parallel(fingerprint_case, request["input"], threads, progress)
    if generator in ("generate_mmff_has_all_molecule_params", "generate_mmff_optimize", "generate_mmff_optimize_conformers"):
        results = parallel(mmff_case, [(case, request["parameters"]) for case in request["corpus"]], threads, progress)
        return [row for batch in results for row in batch]
    if generator == "generate_smiles_write":
        cases, profiles = request["corpus"], request["parameters"]
        rows = writer_rows(cases, profiles, threads, progress)
        records = [{"input": {"SmilesWrite": {"case": case, "profile": profile}}, "output": {"SmilesWrite": outcome}}
                   for case, values in zip(cases, rows, strict=True)
                   for profile, outcome in zip(profiles, values, strict=True)]
        if [record["input"] for record in records] != request["input"]:
            raise ValueError("writer input order/parameters changed")
        return records
    if generator in ("generate_tautomer_enumeration", "generate_tautomer_canonicalization"):
        return parallel(tautomer_case, request["input"], threads, progress)
    from fingerprint_values_pilot import GENERATORS
    return GENERATORS[generator](request["corpus"], request["parameters"], threads, progress=progress)


if __name__ == "__main__":
    if sys.argv[1:] == ["--structure-case"]:
        from _generate_tetrahedral_stereo_geometry import run_assign_case
        json.dump(run_assign_case(json.load(sys.stdin)),sys.stdout)
    else:
        json.dump(generate(json.load(sys.stdin)), sys.stdout, allow_nan=False)
