"""Direct pinned RDKit fingerprint calls; no dependency on the former runner."""
from functools import lru_cache
import json
from rdkit.Chem import rdFingerprintGenerator as F
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys


class UInt(int):
    """Unsigned wire integer (including the 64-bit sparse-vector length)."""


@lru_cache(maxsize=256)
def morgan_generator(serialized):
    p = json.loads(serialized)
    atom = (F.GetMorganFeatureAtomInvGen() if p["feature_invariants"] else
            F.GetMorganAtomInvGen(includeRingMembership=p["include_ring_membership"]))
    bond = None if p["bond_chirality"] is None else F.GetMorganBondInvGen(
        useBondTypes=p["use_bond_types"], useChirality=p["bond_chirality"])
    generator = F.GetMorganGenerator(
        radius=p["radius"], countSimulation=p["count_simulation"],
        includeChirality=p["include_chirality"], useBondTypes=p["use_bond_types"],
        onlyNonzeroInvariants=p["only_nonzero_invariants"],
        includeRingMembership=p["include_ring_membership"], countBounds=p["count_bounds"],
        fpSize=p["fp_size"], atomInvariantsGenerator=atom, bondInvariantsGenerator=bond,
        includeRedundantEnvironments=p["include_redundant_environments"])
    generator.GetOptions().numBitsPerFeature = p["bits_per_feature"]
    return generator


def fingerprint_value(fp, form):
    if form in ("bits", "sparse_bits"):
        return [UInt(fp.GetNumBits()), list(fp.GetOnBits())]
    return [UInt(fp.GetLength()), [[UInt(bit), count] for bit, count in sorted(fp.GetNonzeroElements().items())]]


def provenance(output):
    return [list(output.GetAtomCounts()),
            [[UInt(bit) for bit in row] for row in output.GetAtomToBits()],
            [[UInt(bit), list(rows)] for bit, rows in sorted(output.GetBitInfoMap().items())],
            [[UInt(bit), list(rows)] for bit, rows in sorted(output.GetAtomsPerBit().items())]]


def morgan(mol, p):
    generator = morgan_generator(json.dumps(p, sort_keys=True, separators=(",", ":")))
    kwargs = {}
    if p["from_first"]:
        kwargs["fromAtoms"] = [0] if mol.GetNumAtoms() else []
    if p["ignore_first"]:
        kwargs["ignoreAtoms"] = [0] if mol.GetNumAtoms() else []
    if p["custom_atoms"]:
        kwargs["customAtomInvariants"] = [0 if p["only_nonzero_invariants"] and i % p["zero_every"] == 0 else i + p["custom_atom_offset"]
                                          for i in range(mol.GetNumAtoms())]
    if p["custom_bonds"]:
        kwargs["customBondInvariants"] = [i + p["custom_bond_offset"] for i in range(mol.GetNumBonds())]
    output = None
    if p["additional_output"]:
        output = F.AdditionalOutput()
        for allocate in ("AllocateAtomCounts", "AllocateAtomToBits", "AllocateBitInfoMap", "AllocateAtomsPerBit"):
            getattr(output, allocate)()
        kwargs["additionalOutput"] = output
    method = {"bits": "GetFingerprint", "sparse_bits": "GetSparseFingerprint",
              "count": "GetCountFingerprint", "sparse_count": "GetSparseCountFingerprint"}[p["form"]]
    fp = getattr(generator, method)(mol, **kwargs)
    return [fingerprint_value(fp, p["form"]), None if output is None else provenance(output)]


@lru_cache(maxsize=256)
def pair_generator(task, serialized):
    p = json.loads(serialized)
    common = dict(includeChirality=p["include_chirality"], countSimulation=p["count_simulation"],
                  countBounds=p["count_bounds"], fpSize=p["fp_size"])
    if task == "fingerprint_atom_pair":
        generator = F.GetAtomPairGenerator(minDistance=p["min_distance"], maxDistance=p["max_distance"],
                                          use2D=p["use_2d"], **common)
    else:
        generator = F.GetTopologicalTorsionGenerator(torsionAtomCount=p["torsion_atom_count"], **common)
        generator.GetOptions().onlyShortestPaths = p["only_shortest_paths"]
    generator.GetOptions().numBitsPerFeature = p["bits_per_feature"]
    return generator


def pair_fingerprint(mol, task, p):
    generator = pair_generator(task, json.dumps(p, sort_keys=True, separators=(",", ":")))
    kwargs = {}
    if p["from_first"]:
        kwargs["fromAtoms"] = [0] if mol.GetNumAtoms() else []
    if p["ignore_first"]:
        kwargs["ignoreAtoms"] = [0] if mol.GetNumAtoms() else []
    if p["custom_atoms"]:
        offset = 11 if task == "fingerprint_atom_pair" else 17
        kwargs["customAtomInvariants"] = [i + offset for i in range(mol.GetNumAtoms())]
    output = None
    if p["additional_output"]:
        output = F.AdditionalOutput()
        for allocate in ("AllocateAtomCounts", "AllocateAtomToBits", "AllocateBitInfoMap", "AllocateAtomsPerBit", "AllocateBitPaths"):
            getattr(output, allocate)()
        kwargs["additionalOutput"] = output
    method = {"bits": "GetFingerprint", "sparse_bits": "GetSparseFingerprint",
              "count": "GetCountFingerprint", "sparse_count": "GetSparseCountFingerprint"}[p["form"]]
    fp = getattr(generator, method)(mol, **kwargs)
    extra = None if output is None else provenance(output) + [
        [[UInt(bit), list(rows)] for bit, rows in sorted(output.GetBitPaths().items())]]
    return [fingerprint_value(fp, p["form"]), extra]


def roots(mode, count):
    if mode == "all":
        return None
    if not count:
        return []
    return [0] if mode == "first" or count == 1 else [0, count - 1]


def path_fingerprint(mol, task, p):
    if task == "fingerprint_fuzzy":
        other = Chem.MolFromSmiles(p["other_smiles"])
        if other is None:
            raise ValueError("invalid explicit pair input")
        generator = F.GetMorganGenerator(radius=p["radius"],fpSize=p["fp_size"])
        def make(molecule):
            source = generator.GetCountFingerprint(molecule)
            if not p["wide"]:
                return source
            value = DataStructs.ULongSparseIntVect(1 << 40)
            for bit, count in source.GetNonzeroElements().items():
                value[bit + (1 << 32)] = count
            return value
        left, right = make(mol), make(other)
        return [fingerprint_value(left & right,"count"),fingerprint_value(left | right,"count")]
    if task == "fingerprint_topological":
        atom_bits, info = [], {}
        kwargs = dict(minPath=p["min_path"], maxPath=p["max_path"], fpSize=p["fp_size"],
                      nBitsPerHash=p["bits_per_feature"], useHs=p["use_hs"],
                      tgtDensity=p["target_density"], minSize=p["min_size"],
                      branchedPaths=p["branched_paths"], useBondOrder=p["use_bond_order"],
                      atomBits=atom_bits, bitInfo=info)
        if p["custom_atoms"]:
            kwargs["atomInvariants"] = list(range(1, mol.GetNumAtoms() + 1))
        selected = roots(p["roots"], mol.GetNumAtoms())
        if selected is not None:
            kwargs["fromAtoms"] = selected
        fp = Chem.RDKFingerprint(mol, **kwargs)
        return [fingerprint_value(fp, "bits"), [[UInt(bit) for bit in row] for row in atom_bits],
                [[UInt(bit), list(rows)] for bit, rows in sorted(info.items())]]
    if task == "fingerprint_layered":
        kwargs = dict(layerFlags=p["layers"], minPath=p["min_path"], maxPath=p["max_path"],
                      fpSize=p["fp_size"], branchedPaths=p["branched_paths"])
        selected = roots(p["roots"], mol.GetNumAtoms())
        if selected is not None:
            kwargs["fromAtoms"] = selected
        counts = None
        if p["counts"] != "absent":
            counts = [i + 10 if p["counts"] == "seeded" else 0 for i in range(mol.GetNumAtoms())]
            kwargs["atomCounts"] = counts
        mask = None
        if p["mask"] != "absent":
            mask = DataStructs.ExplicitBitVect(p["fp_size"])
            stride = {"even": 2, "mod_three": 3}.get(p["mask"])
            if stride:
                for i in range(0, p["fp_size"], stride):
                    mask.SetBit(i)
            kwargs["setOnlyBits"] = mask
        before = None if mask is None else tuple(mask.GetOnBits())
        fp = Chem.LayeredFingerprint(mol, **kwargs)
        if before != (None if mask is None else tuple(mask.GetOnBits())):
            raise RuntimeError("native Layered changed input mask")
        return [fingerprint_value(fp, "bits"), counts]
    if task == "fingerprint_pattern":
        kwargs = dict(fpSize=p["fp_size"], tautomerFingerprints=p["tautomeric"])
        # The source checks these optional vectors but does not use their contents.
        # Keep native-only option coverage; CK deliberately exposes no inert options.
        if p["native_counts"] != "absent":
            kwargs["atomCounts"] = [i + 11 if p["native_counts"] == "index_plus_11" else 0 for i in range(mol.GetNumAtoms())]
        if p["native_mask"] != "absent":
            mask = DataStructs.ExplicitBitVect(p["fp_size"])
            selected = range(p["fp_size"]) if p["native_mask"] == "all" else (
                range(0, p["fp_size"], 3) if p["native_mask"] == "sparse" else [])
            for i in selected:
                mask.SetBit(i)
            kwargs["setOnlyBits"] = mask
        return fingerprint_value(Chem.PatternFingerprint(mol, **kwargs), "bits")
    if task == "fingerprint_maccs":
        fp = MACCSkeys.GenMACCSKeys(mol)
        if p["raw"]:
            return fingerprint_value(fp, "bits")
        # maccs.rs::maccs_fingerprint: CK's explicit public projection removes
        # RDKit's unused bit 0. This is a contract projection, not output fitting.
        if fp.GetBit(0):
            raise RuntimeError("RDKit unexpectedly set the unused MACCS bit 0")
        return [UInt(166), [bit - 1 for bit in fp.GetOnBits()]]
    if task == "fingerprint_avalon":
        from rdkit.Avalon import pyAvalonTools
        fp = pyAvalonTools.GetAvalonFP(mol, nBits=p["fp_size"], isQuery=p["is_query"], bitFlags=p["bit_flags"])
        return fingerprint_value(fp, "bits")
    raise ValueError(f"unknown path fingerprint {task}")
