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


def molalign_corpus_case(wrapped):
    from _generate_molalign_golden import corpus_case, record_for_call
    recipe = wrapped["MolAlign"]
    case = corpus_case(recipe["row"], recipe["case"]["smiles"])
    record = record_for_call(case, 0, case["calls"][0])
    prepared = dict(recipe, preparation={key: record[key] for key in
                    ("case_id", "call_index", "operation", "source", "parameters")})
    # RDKit's mutating result is compared with CK's returned value. CK's
    # non-mutating public API must additionally preserve its source molecule.
    record["source_after"] = record["before"]
    return {"input": {"MolAlign": prepared}, "output": {"MolAlign": record}}


def molalign_focused_call(payload):
    from _generate_molalign_golden import record_for_call
    return record_for_call(*payload)


def persistent_forcefield_case(wrapped):
    """One owned evaluator, common exact arbitrary positions, at most two steps."""
    import hashlib
    import random
    import struct
    from rdkit import Chem
    from rdkit.Chem import AllChem

    def bits(value):
        return struct.unpack("<Q", struct.pack("<d", value))[0]

    recipe = wrapped["PersistentForceField"]
    case = recipe["case"]
    seed = hashlib.sha256(
        f"{recipe['seed']}:{case['id']}:{case['smiles']}".encode()).digest()
    rng = random.Random(int.from_bytes(seed, "little"))
    molecule = Chem.MolFromSmiles(case["smiles"])
    if molecule is None:
        raise ValueError(f"{case['id']}: cannot prepare force-field molecule")
    molecule = Chem.AddHs(molecule)
    count = molecule.GetNumAtoms()
    if not count:
        raise ValueError(f"{case['id']}: empty force-field molecule")
    xyz = [[rng.uniform(-3.0, 3.0) for _ in range(3)] for _ in range(count)]
    conformer = Chem.Conformer(count)
    conformer.Set3D(True)
    for index, position in enumerate(xyz):
        conformer.SetAtomPosition(index, position)
    molecule.AddConformer(conformer, assignId=True)
    molblock = Chem.MolToMolBlock(molecule)
    # Use the same transported chemistry on both sides, but never its rounded
    # MolBlock coordinates: reinstall the original binary64 positions.
    molecule = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
    if molecule is None or molecule.GetNumAtoms() != count:
        raise ValueError(f"{case['id']}: force-field MolBlock transport failed")
    for index, position in enumerate(xyz):
        molecule.GetConformer(0).SetAtomPosition(index, position)
    prepared = dict(recipe, preparation={
        "molblock": molblock, "atom_count": count,
        "coordinate_rows": [{"conformer_id": 0, "xyz_bits": [[bits(v) for v in p] for p in xyz]}],
    })
    if recipe["kind"] == "Mmff":
        properties = AllChem.MMFFGetMoleculeProperties(molecule, mmffVariant="MMFF94")
        field = None if properties is None else AllChem.MMFFGetMoleculeForceField(
            molecule, properties, nonBondedThresh=100.0, confId=0,
            ignoreInterfragInteractions=True)
    elif recipe["kind"] == "Uff":
        try:
            field = AllChem.UFFGetMoleculeForceField(
                molecule, vdwThresh=10.0, confId=0, ignoreInterfragInteractions=True)
        except RuntimeError as error:
            # Pinned Builder.cpp:422 selects SP3D/degree-five centers; :324-327
            # passes their unchecked params to AngleBend.cpp:79. The native
            # parameter query identifies the missing center without a case-ID
            # whitelist or treating every RuntimeError as an expected rejection.
            lines = [line.strip() for line in str(error).splitlines() if line.strip()]
            if lines != ["Pre-condition Violation", "bad params pointer",
                         "Violation occurred on line 78 in file Code/ForceField/UFF/AngleBend.cpp",
                         "Failed Expression: at2Params", "RDKIT: 2026.03.6", "BOOST: 1_85"]:
                raise
            centers = [atom.GetIdx() for atom in molecule.GetAtoms()
                       if atom.GetHybridization() == Chem.HybridizationType.SP3D
                       and atom.GetDegree() == 5
                       and AllChem.GetUFFVdWParams(molecule, atom.GetIdx(), atom.GetIdx()) is None]
            if len(centers) != 1:
                raise
            return {"input": {"PersistentForceField": prepared}, "output": {
                "PersistentForceField": {"SourceTbpCenterParamsMissing": {"center_atom_index": centers[0]}}}}
    else:
        raise ValueError(f"unknown owned force-field kind {recipe['kind']}")
    if field is None:
        if recipe["kind"] != "Mmff":
            raise ValueError(f"{case['id']}: UFF factory unexpectedly returned None")
        output = "Unavailable"
    else:
        field.Initialize()

        def snapshot():
            positions = [list(molecule.GetConformer(0).GetAtomPosition(i)) for i in range(count)]
            flat = [v for p in positions for v in p]
            # Explicit current coordinates reset RDKit's distance matrix,
            # matching CK's freshly evaluated post-minimize energy/gradient.
            energy = field.CalcEnergy(flat)
            gradient = field.CalcGrad(flat)
            return {"energy_bits": bits(energy),
                    "gradient_bits": [[bits(v) for v in gradient[i:i + 3]] for i in range(0, len(gradient), 3)],
                    "positions_bits": [[bits(v) for v in p] for p in positions]}

        initial = snapshot()
        tolerance = struct.unpack("<d", struct.pack("<Q", recipe["force_tolerance_bits"]))[0]
        energy_tolerance = struct.unpack("<d", struct.pack("<Q", recipe["energy_tolerance_bits"]))[0]
        status = field.Minimize(maxIts=recipe["max_iterations"], forceTol=tolerance, energyTol=energy_tolerance)
        if status not in (0, 1):
            raise ValueError(f"{case['id']}: unexpected minimize status {status}")
        output = {"Evaluated": {"initial": initial, "final_state": snapshot(), "converged": status == 0}}
    return {"input": {"PersistentForceField": prepared}, "output": {"PersistentForceField": output}}


def tautomer_special_case(payload):
    from _tautomer_oracle import build_record
    case, branches = payload
    return build_record(case, {"branches": branches}, [branch["name"] for branch in branches])


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
    from forcefield_preparation import prepare_geometry
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
                row["preparation"] = prepare_geometry(common, label="MMFF")
                if "TimedOut" in row["preparation"]:
                    rows.append({"input": {"Mmff": row}, "output": {"Mmff": row["preparation"]}})
                    continue
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


def writer_native_profile(profile):
    """Project CK option names to the pinned oracle's internal parameter keys."""
    params = dict(profile)
    params["do_isomeric_smiles"] = params.pop("isomeric_smiles")
    params["do_kekule"] = params.pop("kekule")
    if params["rooted_at_atom"] == "none":
        params["rooted_at_atom"] = None
    return params


def writer_case(payload):
    from rdkit import Chem
    from _generate_smiles_writer_golden import branch_result
    case, profiles = payload
    molecule = Chem.MolFromSmiles(case["smiles"])
    if molecule is None:
        return [{"Error": "Parse"} for _ in profiles]
    outcomes = []
    for profile in profiles:
        params = writer_native_profile(profile)
        result = branch_result(molecule, {"params": params})
        outcomes.append({"Smiles": result["smiles"]} if result["ok"] else {"Error": "Write"})
    return outcomes


def writer_rows(cases, profiles, threads, progress):
    from _generate_smiles_writer_golden import iter_branches
    existing = [branch["params"] for branch in iter_branches()]
    if [writer_native_profile(profile) for profile in profiles] != existing:
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


def bio_mmcif_switch_case(payload):
    import gemmi
    case, all_groups, flag = payload
    with tempfile.TemporaryDirectory(prefix="fixed-gemmi-mmcif-") as folder:
        # A PDB structure's name comes from its filename and becomes the CIF block name.
        source = Path(folder) / Path(case["input"]).name
        source.write_text(case["text"], encoding="utf-8")
        # Match the original C++ read_structure_file (no merge_chain_parts).
        structure = gemmi.read_structure(str(source), merge_chain_parts=False)
        groups = gemmi.MmcifOutputGroups(all_groups)
        setattr(groups, flag, True if flag == "auth_all" else not all_groups)
        text = structure.make_mmcif_document(groups).as_string()
    value = True if flag == "auth_all" else not all_groups
    return {"case_id": case["case_id"], "all_groups": all_groups,
            "flag": flag, "value": value, "text": text}


def uff_case(payload):
    from fingerprint_values_pilot import uff
    from forcefield_preparation import prepare_geometry
    case, profiles, first_case_id = payload
    rows = []
    for profile in profiles:
        row = {"case": case, "profile": profile, "preparation": None}
        preparation = prepare_geometry(row, first_case_id=first_case_id)
        output = preparation if preparation is not None and "TimedOut" in preparation else uff(row, first_case_id=first_case_id)
        rows.append({"input": {"Uff": row}, "output": {"Uff": output}})
    return rows


def generate(request):
    threads = request["threads"]
    if not isinstance(threads,int) or threads < 1:
        raise ValueError("threads must be positive")
    kind = request["kind"]
    progress = Progress(request["task"])
    if kind == "bio_mmcif_switches" or (kind == "corpus" and request["generator"] in ("generate_bio_pdb_output_pdb", "generate_bio_pdb_output_cif")):
        import gemmi
        pin = json.loads((PACKAGE / "testdata/reference/gemmi.json").read_text())
        if gemmi.__version__ != pin["version"]:
            raise RuntimeError(f"Gemmi version {gemmi.__version__} != {pin['version']}")
        if kind == "bio_mmcif_switches":
            fixture = request["input"]
            return parallel(bio_mmcif_switch_case,
                            [(case, all_groups, flag) for case in fixture["cases"]
                             for all_groups in (False, True) for flag in fixture["flags"]],
                            threads, progress)
        return parallel(bio_case, request["input"], threads, progress)
    from rdkit import rdBase
    pin = json.loads((PACKAGE / "testdata/reference/rdkit.json").read_text())
    if rdBase.rdkitVersion != pin["version"]:
        raise RuntimeError(f"RDKit version {rdBase.rdkitVersion} != {pin['version']}")
    if kind == "descriptors":
        return parallel(descriptor_case, request["corpus"], threads, progress)
    if kind in ("forcefield_optimizers", "mmff_builtin"):
        from forcefield_regression import optimizer_case, builtin_case
        fixture = request["input"]
        function = optimizer_case if kind == "forcefield_optimizers" else builtin_case
        return parallel(function, fixture["cases"], threads, progress)
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
    if kind in ("tautomer_long_conjugated", "tautomer_focused"):
        fixture = request["input"]
        return parallel(tautomer_special_case,
                        [(case, fixture["branches"]) for case in fixture["cases"]],
                        threads, progress)
    if kind == "molalign_focused":
        return parallel(molalign_focused_call,
                        [(case, index, call) for case in request["input"]["cases"]
                         for index, call in enumerate(case["calls"])], threads, progress)
    if kind != "corpus":
        raise ValueError(f"unknown recipe {kind}")
    generator = request["generator"]
    if generator in ("generate_uff_has_all_molecule_params", "generate_uff_optimize", "generate_uff_optimize_conformers"):
        first_case_ids = {}
        for case in request["corpus"]:
            first_case_ids.setdefault(case["smiles"], case["id"])
        results = parallel(uff_case, [(case, request["parameters"], first_case_ids[case["smiles"]]) for case in request["corpus"]], threads, progress)
        return [row for batch in results for row in batch]
    if generator == "generate_persistent_force_field":
        return parallel(persistent_forcefield_case, request["input"], threads, progress)
    if generator == "generate_molalign":
        return parallel(molalign_corpus_case, request["input"], threads, progress)
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
