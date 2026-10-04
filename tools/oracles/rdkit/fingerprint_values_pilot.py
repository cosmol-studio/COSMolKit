"""Reference adapter only; Rust owns task selection, inputs and comparison."""
import json
import sys
import struct
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors, rdDepictor, rdFingerprintGenerator, AllChem
from rdkit.Chem.Draw import rdMolDraw2D


def topology(mol):
    return {"atoms": [{
        "id": a.GetIdx(), "atomic_number": a.GetAtomicNum(),
        "isotope": a.GetIsotope(), "charge": a.GetFormalCharge(),
        "explicit_h": a.GetNumExplicitHs(), "no_implicit": a.GetNoImplicit(),
        "radicals": a.GetNumRadicalElectrons(), "aromatic": a.GetIsAromatic(),
        "hybridization": int(a.GetHybridization()), "chiral": int(a.GetChiralTag()),
        "permutation": a.GetUnsignedProp("_chiralPermutation") if a.HasProp("_chiralPermutation") else 0,
        "atom_map": a.GetAtomMapNum(),
    } for a in mol.GetAtoms()], "bonds": [{
        "id": b.GetIdx(), "begin": b.GetBeginAtomIdx(), "end": b.GetEndAtomIdx(),
        "order": int(b.GetBondType()), "aromatic": b.GetIsAromatic(),
        "conjugated": b.GetIsConjugated(), "direction": int(b.GetBondDir()),
        "stereo": int(b.GetStereo()), "stereo_atoms": list(b.GetStereoAtoms()),
    } for b in mol.GetBonds()]}


def morgan_additional_output(output):
    atom_counts = output.GetAtomCounts()
    atom_to_bits = output.GetAtomToBits()
    bit_info_map = output.GetBitInfoMap()
    bit_paths = output.GetBitPaths()
    atoms_per_bit = output.GetAtomsPerBit()

    return {
        "atom_counts": None if atom_counts is None else [int(value) for value in atom_counts],
        "atom_to_bits": None if atom_to_bits is None else [
            [int(bit) for bit in atom_bits] for atom_bits in atom_to_bits
        ],
        "bit_info_map": None if bit_info_map is None else [
            [int(bit), [[int(atom), int(radius)] for atom, radius in environments]]
            for bit, environments in sorted(bit_info_map.items())
        ],
        "bit_paths": None if bit_paths is None else [
            [int(bit), [[int(index) for index in path] for path in paths]]
            for bit, paths in sorted(bit_paths.items())
        ],
        "atoms_per_bit": None if atoms_per_bit is None else [
            [int(bit), [[int(atom) for atom in atoms] for atoms in groups]]
            for bit, groups in sorted(atoms_per_bit.items())
        ],
    }


def morgan_fingerprint(mol, options):
    invariant_kind = options["invariants"]
    if invariant_kind == "Connectivity":
        atom_invariants = rdFingerprintGenerator.GetMorganAtomInvGen(True)
    elif invariant_kind == "Features":
        atom_invariants = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
    else:
        raise ValueError(f"unregistered Morgan invariant kind: {invariant_kind}")

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=options["radius"],
        countSimulation=options["count_simulation"],
        includeChirality=options["include_chirality"],
        useBondTypes=True,
        onlyNonzeroInvariants=False,
        includeRingMembership=True,
        countBounds=[1, 2, 4, 8],
        fpSize=2048,
        atomInvariantsGenerator=atom_invariants,
        includeRedundantEnvironments=False,
    )

    additional = rdFingerprintGenerator.AdditionalOutput()
    additional.AllocateAtomCounts()
    additional.AllocateAtomToBits()
    additional.AllocateBitInfoMap()
    additional.AllocateBitPaths()
    additional.AllocateAtomsPerBit()

    output = options["output"]
    if output == "DenseBits":
        fingerprint = generator.GetFingerprint(mol, additionalOutput=additional)
        return {"MorganDenseBits": {
            "length": int(fingerprint.GetNumBits()),
            "on_bits": [int(bit) for bit in fingerprint.GetOnBits()],
            "additional_output": morgan_additional_output(additional),
        }}
    if output == "SparseBits":
        fingerprint = generator.GetSparseFingerprint(mol, additionalOutput=additional)
        return {"MorganSparseBits": {
            "length": int(fingerprint.GetNumBits()),
            "on_bits": [int(bit) for bit in fingerprint.GetOnBits()],
            "additional_output": morgan_additional_output(additional),
        }}
    if output == "HashedCounts":
        fingerprint = generator.GetCountFingerprint(mol, additionalOutput=additional)
        return {"MorganHashedCounts": {
            "length": int(fingerprint.GetLength()),
            "entries": [[int(bit), int(count)] for bit, count in sorted(
                fingerprint.GetNonzeroElements().items()
            )],
            "additional_output": morgan_additional_output(additional),
        }}
    if output == "SparseCounts":
        fingerprint = generator.GetSparseCountFingerprint(mol, additionalOutput=additional)
        return {"MorganSparseCounts": {
            "length": int(fingerprint.GetLength()),
            "entries": [[int(bit), int(count)] for bit, count in sorted(
                fingerprint.GetNonzeroElements().items()
            )],
            "additional_output": morgan_additional_output(additional),
        }}
    raise ValueError(f"unregistered Morgan output kind: {output}")


def molecular(row):
    profile = row["profile"]
    name, options = (profile, {}) if isinstance(profile, str) else next(iter(profile.items()))
    stage = "Parse"
    try:
        params = Chem.SmilesParserParams()
        params.sanitize = options["sanitize"] if name == "SmilesRead" else name != "SanitizeAll"
        params.removeHs = options["remove_hydrogens"] if name in ("SmilesRead", "NumHeavyAtoms", "TotalAtomCount", "LipinskiHBA", "LipinskiHBD", "FractionCSP3", "NumHeteroatoms", "NumHba", "NumHbd", "NumRings", "NumHeterocycles", "NumAromaticRings", "NumSaturatedRings", "NumAliphaticRings", "NumAromaticHeterocycles", "NumAromaticCarbocycles", "NumAliphaticHeterocycles", "NumAliphaticCarbocycles", "NumSaturatedHeterocycles", "NumSaturatedCarbocycles") else name != "SanitizeAll"
        params.allowCXSMILES = True
        params.strictCXSMILES = True
        params.parseName = True
        # The pinned Python v1 parser does not expose v2 skipCleanup;
        # its C++ conversion leaves that field at the default false.
        mol = Chem.MolFromSmiles(row["case"]["smiles"], params)
        if mol is None:
            raise ValueError("RDKit MolFromSmiles returned None")
        stage = "Operation"
        if name == "SvgDefault":
            drawer = rdMolDraw2D.MolDraw2DSVG(300, 300, -1, -1, True)
            rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
            drawer.FinishDrawing()
            return {"Text": drawer.GetDrawingText()}
        if name == "DistanceMatrix":
            matrix = Chem.GetDistanceMatrix(mol, useBO=options["use_bond_order"],
                useAtomWts=options["use_atom_weights"], force=True)
            return {"Matrix": {"dimension": mol.GetNumAtoms(), "values_bits": [
                struct.unpack(">Q", struct.pack(">d", float(value)))[0]
                for value in matrix.flat]}}
        if name == "MolecularWeight":
            value = Descriptors.MolWt(mol, onlyHeavy=options["only_heavy"])
            return {"Float64Bits": struct.unpack(">Q", struct.pack(">d", value))[0]}
        if name == "ExactMolecularWeight":
            value = Descriptors.ExactMolWt(mol, onlyHeavy=options["only_heavy"])
            return {"Float64Bits": struct.unpack(">Q", struct.pack(">d", value))[0]}
        if name == "MolecularFormula":
            return {"Text": rdMolDescriptors.CalcMolFormula(mol,
                separateIsotopes=options["separate_isotopes"],
                abbreviateHIsotopes=options["abbreviate_h_isotopes"])}
        if name == "NumHeavyAtoms":
            return {"Unsigned": rdMolDescriptors.CalcNumHeavyAtoms(mol)}
        if name == "TotalAtomCount":
            return {"Unsigned": rdMolDescriptors.CalcNumAtoms(mol)}
        if name == "LipinskiHBA":
            return {"Unsigned": rdMolDescriptors.CalcNumLipinskiHBA(mol)}
        if name == "LipinskiHBD":
            return {"Unsigned": rdMolDescriptors.CalcNumLipinskiHBD(mol)}
        if name == "FractionCSP3":
            return {"Float64Bits": struct.unpack(">Q", struct.pack(">d", rdMolDescriptors.CalcFractionCSP3(mol)))[0]}
        if name == "NumHeteroatoms":
            return {"Unsigned": rdMolDescriptors.CalcNumHeteroatoms(mol)}
        if name == "NumHba":
            return {"Unsigned": rdMolDescriptors.CalcNumHBA(mol)}
        if name == "NumHbd":
            return {"Unsigned": rdMolDescriptors.CalcNumHBD(mol)}
        if name == "NumRings":
            return {"Unsigned": rdMolDescriptors.CalcNumRings(mol)}
        if name == "NumHeterocycles":
            return {"Unsigned": rdMolDescriptors.CalcNumHeterocycles(mol)}
        if name == "NumAromaticRings":
            return {"Unsigned": rdMolDescriptors.CalcNumAromaticRings(mol)}
        if name == "NumSaturatedRings":
            return {"Unsigned": rdMolDescriptors.CalcNumSaturatedRings(mol)}
        if name == "NumAliphaticRings":
            return {"Unsigned": rdMolDescriptors.CalcNumAliphaticRings(mol)}
        if name == "NumAromaticHeterocycles":
            return {"Unsigned": rdMolDescriptors.CalcNumAromaticHeterocycles(mol)}
        if name == "NumAromaticCarbocycles":
            return {"Unsigned": rdMolDescriptors.CalcNumAromaticCarbocycles(mol)}
        if name == "NumAliphaticHeterocycles":
            return {"Unsigned": rdMolDescriptors.CalcNumAliphaticHeterocycles(mol)}
        if name == "NumAliphaticCarbocycles":
            return {"Unsigned": rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)}
        if name == "NumSaturatedHeterocycles":
            return {"Unsigned": rdMolDescriptors.CalcNumSaturatedHeterocycles(mol)}
        if name == "NumSaturatedCarbocycles":
            return {"Unsigned": rdMolDescriptors.CalcNumSaturatedCarbocycles(mol)}
        if name == "Morgan":
            return morgan_fingerprint(mol, options)
        if name == "SanitizeAll":
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
        elif name == "Kekulize":
            Chem.Kekulize(mol, clearAromaticFlags=options["clear_aromatic_flags"])
        elif name == "AddHydrogens":
            params = Chem.AddHsParameters()
            params.explicitOnly = options["explicit_only"]
            params.addCoords = False
            params.addResidueInfo = False
            params.skipQueries = False
            mol = Chem.AddHs(mol, params)
        elif name == "RemoveHydrogens":
            stage = "Preparation"
            params = Chem.AddHsParameters()
            params.explicitOnly = False
            params.addCoords = False
            params.addResidueInfo = False
            params.skipQueries = False
            mol = Chem.AddHs(mol, params)
            stage = "Operation"
            params = Chem.RemoveHsParameters()
            for key in ("removeDegreeZero", "removeHigherDegrees", "removeOnlyHNeighbors",
                        "removeIsotopes", "removeAndTrackIsotopes", "removeDummyNeighbors",
                        "removeDefiningBondStereo", "removeWithQuery", "updateExplicitCount",
                        "removeHydrides", "removeNontetrahedralNeighbors"):
                setattr(params, key, False)
            for key in ("removeWithWedgedBond", "removeMapped", "removeInSGroups",
                        "showWarnings", "removeNonimplicit"):
                setattr(params, key, True)
            mol = Chem.RemoveHs(mol, params, sanitize=options["sanitize"])
        elif name == "Coordinates2dDefault":
            rdDepictor.SetPreferCoordGen(False)
            rdDepictor.Compute2DCoords(mol, canonOrient=False, clearConfs=True,
                coordMap={}, nFlipsPerSample=0, nSample=0, sampleSeed=0,
                permuteDeg4Nodes=False, forceRDKit=False, useRingTemplates=False)
            conf = mol.GetConformer()
            return {"Coordinates2d": {"topology": topology(mol), "xy_bits": [
                [struct.unpack(">Q", struct.pack(">d", coord))[0] for coord in
                 (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)]
                for i in range(mol.GetNumAtoms())]}}
        elif name != "SmilesRead":
            raise RuntimeError(f"unregistered molecular profile: {name}")
        return {"Topology": topology(mol)}
    except (ValueError, RuntimeError) as error:
        return {"Error": {"stage": stage, "detail": f"{type(error).__name__}: {error}"}}


common_geometry = {}
common_conformer_geometry = {}

def uff(row, first_case_id=None):
    name, options = next(iter(row["profile"].items()))
    stage = "Parse"
    if name in ("Optimization", "ConformerOptimization") and row["preparation"] is None:
        # Input preparation is separate from optimization. Source rejections
        # are explicit records, not absent artifacts or invented coordinates.
        key = (row["case"]["smiles"], options["add_hydrogens"])
        preparation_case_id = row["case"]["id"] if first_case_id is None else first_case_id
        if key not in common_geometry:
            try:
                seed_mol = Chem.MolFromSmiles(key[0])
                if seed_mol is None:
                    raise ValueError(f"UFF common geometry parse failed: {preparation_case_id}")
                stage = "Preparation"
                if key[1]:
                    seed_mol = Chem.AddHs(seed_mol)
                params = AllChem.ETKDGv3()
                params.randomSeed = 61453
                params.useRandomCoords = True
                params.numThreads = 1
                if AllChem.EmbedMolecule(seed_mol, params) != 0:
                    raise ValueError(f"UFF common geometry embedding failed: {preparation_case_id}")
                # Use the source writer's default Kekule representation.
                # The non-Kekule MolBlock loses pyrrolic [nH] on source reread;
                # this prepares a valid common input, never fits an FF output.
                block = Chem.MolToMolBlock(seed_mol, kekulize=True, forceV3000=True)
                common_geometry[key] = {"Ready": {
                    "molblock": block,
                    "atom_count": seed_mol.GetNumAtoms(),
                    "coordinate_rows": [],
                }}
            except (ValueError, RuntimeError) as error:
                common_geometry[key] = {"Rejected": {"stage": stage, "detail": f"{type(error).__name__}: {error}"}}
        if name == "Optimization":
            row["preparation"] = common_geometry[key]
        else:
            count = options["conformer_count"]
            conformer_key = (*key, count)
            if conformer_key not in common_conformer_geometry:
                base = common_geometry[key]
                if "Rejected" in base:
                    common_conformer_geometry[conformer_key] = base
                else:
                    try:
                        stage = "Preparation"
                        base_geometry = base["Ready"]
                        quantized = Chem.MolFromMolBlock(
                            base_geometry["molblock"], sanitize=True, removeHs=False)
                        if (quantized is None
                                or quantized.GetNumAtoms() != base_geometry["atom_count"]
                                or quantized.GetNumConformers() != 1):
                            raise ValueError("invalid quantized UFF base geometry")
                        base_conformer = quantized.GetConformer()
                        coordinate_rows = []
                        for conformer_id, factor in zip(
                                [7, 3, 11][:count], [1.0, 1.1, 1.2][:count]):
                            xyz_bits = []
                            for atom_index in range(quantized.GetNumAtoms()):
                                point = base_conformer.GetAtomPosition(atom_index)
                                xyz_bits.append([
                                    struct.unpack(">Q", struct.pack(">d", float(component) * factor))[0]
                                    for component in (point.x, point.y, point.z)
                                ])
                            coordinate_rows.append({
                                "conformer_id": conformer_id,
                                "xyz_bits": xyz_bits,
                            })
                        common_conformer_geometry[conformer_key] = {"Ready": {
                            "molblock": base_geometry["molblock"],
                            "atom_count": base_geometry["atom_count"],
                            "coordinate_rows": coordinate_rows,
                        }}
                    except (ValueError, RuntimeError) as error:
                        common_conformer_geometry[conformer_key] = {"Rejected": {
                            "stage": stage,
                            "detail": f"{type(error).__name__}: {error}",
                        }}
            row["preparation"] = common_conformer_geometry[conformer_key]
    if name in ("Optimization", "ConformerOptimization") and "Rejected" in row["preparation"]:
        return {"Error": row["preparation"]["Rejected"]}
    try:
        if name == "Coverage":
            mol = Chem.MolFromSmiles(row["case"]["smiles"])
            if mol is None:
                raise ValueError("RDKit MolFromSmiles returned None")
            stage = "Preparation"
            if options["add_hydrogens"]:
                mol = Chem.AddHs(mol)
            stage = "Operation"
            return {"Coverage": AllChem.UFFHasAllMoleculeParams(mol)}
        stage = "Preparation"
        # Reparse the serialized COMMON input, not the pre-serialization double
        # coordinates. CK reads precisely the same MolBlock and keeps explicit H.
        geometry = row["preparation"]["Ready"]
        mol = Chem.MolFromMolBlock(geometry["molblock"], sanitize=True, removeHs=False)
        if mol is None or mol.GetNumAtoms() != geometry["atom_count"]:
            raise ValueError("invalid prepared UFF MolBlock")
        if name == "ConformerOptimization":
            coordinate_rows = geometry["coordinate_rows"]
            if len(coordinate_rows) != options["conformer_count"]:
                raise ValueError("invalid prepared UFF conformer rows")
            mol.RemoveAllConformers()
            for row_coordinates in coordinate_rows:
                conformer = Chem.Conformer(mol.GetNumAtoms())
                conformer_id = row_coordinates["conformer_id"]
                conformer.SetId(conformer_id)
                if len(row_coordinates["xyz_bits"]) != mol.GetNumAtoms():
                    raise ValueError("invalid prepared UFF conformer atom rows")
                for atom_index, xyz_bits in enumerate(row_coordinates["xyz_bits"]):
                    xyz = tuple(
                        struct.unpack(">d", struct.pack(">Q", value))[0]
                        for value in xyz_bits
                    )
                    conformer.SetAtomPosition(atom_index, xyz)
                mol.AddConformer(conformer, assignId=False)
            if mol.GetNumConformers() != options["conformer_count"]:
                raise ValueError("prepared UFF conformer count changed")
            stage = "Operation"
            results = AllChem.UFFOptimizeMoleculeConfs(
                mol,
                numThreads=1,
                maxIters=options["max_iterations"],
                vdwThresh=options["vdw_threshold"],
                ignoreInterfragInteractions=options["ignore_interfragment_interactions"],
            )
            conformers = list(mol.GetConformers())
            if len(results) != len(conformers):
                raise RuntimeError("RDKit UFF conformer result count changed")
            bits = lambda value: struct.unpack(">Q", struct.pack(">d", value))[0]
            return {"OptimizedConformers": {"conformers": [
                {
                    "conformer_id": conformer.GetId(),
                    "status": status,
                    "energy_bits": bits(energy),
                    "xyz_bits": [
                        [bits(value) for value in conformer.GetAtomPosition(atom_index)]
                        for atom_index in range(mol.GetNumAtoms())
                    ],
                }
                for conformer, (status, energy) in zip(conformers, results)
            ]}}
        stage = "Operation"
        conf_id = -1 if options["conformer_id"] is None else options["conformer_id"]
        ff = AllChem.UFFGetMoleculeForceField(mol, vdwThresh=options["vdw_threshold"],
            confId=conf_id, ignoreInterfragInteractions=options["ignore_interfragment_interactions"])
        if ff is None:
            raise ValueError("UFF construction returned None")
        ff.Initialize()
        status = ff.Minimize(maxIts=options["max_iterations"], forceTol=1e-4, energyTol=1e-6)
        energy = ff.CalcEnergy()
        conf = mol.GetConformer(conf_id)
        bits = lambda value: struct.unpack(">Q", struct.pack(">d", value))[0]
        return {"Optimized": {"status": status, "energy_bits": bits(energy), "xyz_bits": [
            [bits(value) for value in conf.GetAtomPosition(i)] for i in range(mol.GetNumAtoms())]}}
    except (ValueError, RuntimeError) as error:
        return {"Error": {"stage": stage, "detail": f"{type(error).__name__}: {error}"}}


def fingerprint(envelope):
    row = envelope["Fingerprint"]
    ctor = {"U32": DataStructs.UIntSparseIntVect,
            "U64": DataStructs.ULongSparseIntVect}[row["width"]]
    case = row["case"]
    def build(entries):
        value = ctor(case["length"])
        for key, count in entries:
            if count == 0:
                value[key] = 1
        value -= 1
        for key, count in entries:
            if count != 0:
                value[key] = count
        return value
    left, right = build(case["left"]), build(case["right"])
    before = (dict(left.GetNonzeroElements()), dict(right.GetNonzeroElements()))
    op = row["operation"]
    if op == "FuzzyAnd":
        result = left & right
    elif op == "FuzzyOr":
        result = left | right
    else:
        raise ValueError(op)
    assert before == (dict(left.GetNonzeroElements()), dict(right.GetNonzeroElements()))
    return {"input": envelope, "output": {"Fingerprint": {
        "length": result.GetLength(),
        "entries": sorted(result.GetNonzeroElements().items()),
    }}}


def _molecular_case(case, parameters):
    rows = []
    for profile in parameters:
        envelope = {"Molecular": {"case": case, "profile": profile}}
        rows.append({"input": envelope, "output": {"Molecular": molecular(envelope["Molecular"])}})
    return rows


def _fingerprint_case(case, parameters):
    return [fingerprint({"Fingerprint": {"case": case, **parameter}}) for parameter in parameters]


def _uff_case(case, parameters, first_case_ids):
    rows = []
    first_case_id = first_case_ids[case["smiles"]]
    for profile in parameters:
        input_row = {"case": case, "profile": profile, "preparation": None}
        output = uff(input_row, first_case_id=first_case_id)
        rows.append({
            "input": {"Uff": input_row},
            "output": {"Uff": output},
        })
    return rows


def _generate(corpus, parameters, threads, worker):
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer")
    if not corpus or not parameters:
        raise ValueError("empty corpus or parameter matrix")
    work = partial(worker, parameters=parameters)
    if threads == 1:
        batches = map(work, corpus)
        return [row for batch in batches for row in batch]
    # Parallelize independent cases, keeping every parameter for one case
    # together. map preserves source order; completion order never labels rows.
    with ProcessPoolExecutor(max_workers=threads) as pool:
        batches = pool.map(work, corpus, chunksize=max(1, len(corpus) // (threads * 4)))
        return [row for batch in batches for row in batch]


def generate_fuzzy_and(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _fingerprint_case)


def generate_fuzzy_or(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _fingerprint_case)


def generate_smiles_read(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_sanitize(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_kekulize(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_molecular_weight(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_exact_molecular_weight(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_molecular_formula(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_heavy_atoms(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_total_atom_count(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_lipinski_hba(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_lipinski_hbd(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_fraction_csp3(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_heteroatoms(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_hba(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_hbd(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_rings(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_heterocycles(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_aromatic_rings(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_saturated_rings(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_aliphatic_rings(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_aromatic_heterocycles(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_aromatic_carbocycles(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_aliphatic_heterocycles(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_aliphatic_carbocycles(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_saturated_heterocycles(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_num_saturated_carbocycles(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_add_hydrogens(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_remove_hydrogens(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_coordinates_2d(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_svg(corpus, parameters, threads):
    if parameters != ["SvgDefault"]:
        raise ValueError("SVG requires exactly the frozen SvgDefault 300x300 profile")
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_distance_matrix(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def _generate_uff(corpus, parameters, threads):
    common_geometry.clear()
    common_conformer_geometry.clear()
    first_case_ids = {}
    for case in corpus:
        first_case_ids.setdefault(case["smiles"], case["id"])
    worker = partial(_uff_case, first_case_ids=first_case_ids)
    return _generate(corpus, parameters, threads, worker)


def generate_uff_has_all_molecule_params(corpus, parameters, threads):
    return _generate_uff(corpus, parameters, threads)


def generate_uff_optimize(corpus, parameters, threads):
    return _generate_uff(corpus, parameters, threads)


def generate_uff_optimize_conformers(corpus, parameters, threads):
    return _generate_uff(corpus, parameters, threads)
def generate_morgan_fingerprint(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_morgan_sparse_fingerprint(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_morgan_count_fingerprint(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


def generate_morgan_sparse_count_fingerprint(corpus, parameters, threads):
    return _generate(corpus, parameters, threads, _molecular_case)


GENERATORS = {
    "generate_fuzzy_and": generate_fuzzy_and,
    "generate_fuzzy_or": generate_fuzzy_or,
    "generate_smiles_read": generate_smiles_read,
    "generate_sanitize": generate_sanitize,
    "generate_kekulize": generate_kekulize,
    "generate_molecular_weight": generate_molecular_weight,
    "generate_exact_molecular_weight": generate_exact_molecular_weight,
    "generate_molecular_formula": generate_molecular_formula,
    "generate_num_heavy_atoms": generate_num_heavy_atoms,
    "generate_total_atom_count": generate_total_atom_count,
    "generate_lipinski_hba": generate_lipinski_hba,
    "generate_lipinski_hbd": generate_lipinski_hbd,
    "generate_fraction_csp3": generate_fraction_csp3,
    "generate_num_heteroatoms": generate_num_heteroatoms,
    "generate_num_hba": generate_num_hba,
    "generate_num_hbd": generate_num_hbd,
    "generate_num_rings": generate_num_rings,
    "generate_num_heterocycles": generate_num_heterocycles,
    "generate_num_aromatic_rings": generate_num_aromatic_rings,
    "generate_num_saturated_rings": generate_num_saturated_rings,
    "generate_num_aliphatic_rings": generate_num_aliphatic_rings,
    "generate_num_aromatic_heterocycles": generate_num_aromatic_heterocycles,
    "generate_num_aromatic_carbocycles": generate_num_aromatic_carbocycles,
    "generate_num_aliphatic_heterocycles": generate_num_aliphatic_heterocycles,
    "generate_num_aliphatic_carbocycles": generate_num_aliphatic_carbocycles,
    "generate_num_saturated_heterocycles": generate_num_saturated_heterocycles,
    "generate_num_saturated_carbocycles": generate_num_saturated_carbocycles,
    "generate_add_hydrogens": generate_add_hydrogens,
    "generate_remove_hydrogens": generate_remove_hydrogens,
    "generate_coordinates_2d": generate_coordinates_2d,
    "generate_svg": generate_svg,
    "generate_distance_matrix": generate_distance_matrix,
    "generate_uff_has_all_molecule_params": generate_uff_has_all_molecule_params,
    "generate_uff_optimize": generate_uff_optimize,
    "generate_uff_optimize_conformers": generate_uff_optimize_conformers,
    "generate_morgan_fingerprint": generate_morgan_fingerprint,
    "generate_morgan_sparse_fingerprint": generate_morgan_sparse_fingerprint,
    "generate_morgan_count_fingerprint": generate_morgan_count_fingerprint,
    "generate_morgan_sparse_count_fingerprint": generate_morgan_sparse_count_fingerprint,
}


def main():
    expected_version = sys.argv[1]
    if rdBase.rdkitVersion != expected_version:
        raise RuntimeError(f"RDKit {rdBase.rdkitVersion} != {expected_version}")
    request = json.load(sys.stdin)
    generator = GENERATORS[request["generator"]]
    rows = generator(request["corpus"], request["parameters"], request["threads"])
    json.dump(rows, sys.stdout, sort_keys=True)


if __name__ == "__main__":
    main()
