"""RDKit calls using the exact seeded parameters recorded in each input row."""
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator


def bits(fp):
    return {"length": fp.GetNumBits(), "on_bits": list(fp.GetOnBits())}


def roots(mode, count):
    if mode == "All":
        return None
    if mode == "First":
        return [0] if count else []
    if mode == "Terminals":
        return [0, count - 1] if count > 1 else list(range(count))
    raise ValueError(f"unknown root mode: {mode}")


def count_vector(fp):
    return {"length": fp.GetLength(),
            "entries": sorted([i, v] for i, v in fp.GetNonzeroElements().items())}


def fingerprint_case(input_row):
    row = input_row["Fingerprint"]
    mol = Chem.MolFromSmiles(row["case"]["smiles"])
    if mol is None:
        raise ValueError(f"invalid SMILES: {row['case']['id']}")
    profile = row["params"]
    if profile == "Maccs":
        raw = bits(MACCSkeys.GenMACCSKeys(mol))
        output = {"Maccs": {"raw": raw, "public": {
            "length": 166, "on_bits": [i - 1 for i in raw["on_bits"] if i]}}}
    else:
        kind, p = next(iter(profile.items()))
        if kind == "Topological":
            kwargs = dict(minPath=p["min_path"], maxPath=p["max_path"],
                          fpSize=p["fp_size"], nBitsPerHash=p["bits_per_feature"],
                          useHs=p["use_hs"], tgtDensity=p["density_milli"] / 1000,
                          minSize=p["min_size"], branchedPaths=p["branched"],
                          useBondOrder=p["bond_order"])
            selected = roots(p["roots"], mol.GetNumAtoms())
            if selected is not None:
                kwargs["fromAtoms"] = selected
            if p["custom_invariants"]:
                kwargs["atomInvariants"] = list(range(1, mol.GetNumAtoms() + 1))
            output = {"Bits": {"fingerprint": bits(Chem.RDKFingerprint(mol, **kwargs)),
                               "atom_counts": None}}
        elif kind == "Layered":
            counts = None if p["counts"] == "Absent" else [
                i + 10 if p["counts"] == "Seeded" else 0 for i in range(mol.GetNumAtoms())]
            mask = None
            if p["mask"] != "Absent":
                mask = DataStructs.ExplicitBitVect(p["fp_size"])
                if p["mask"] != "Empty":
                    for i in range(0, p["fp_size"], 2 if p["mask"] == "Even" else 3):
                        mask.SetBit(i)
            kwargs = dict(layerFlags=p["layers"], minPath=p["min_path"],
                          maxPath=p["max_path"], fpSize=p["fp_size"],
                          branchedPaths=p["branched"], setOnlyBits=mask)
            selected = roots(p["roots"], mol.GetNumAtoms())
            if selected is not None:
                kwargs["fromAtoms"] = selected
            if counts is not None:
                kwargs["atomCounts"] = counts
            mask_before = None if mask is None else bits(mask)
            fp = Chem.LayeredFingerprint(mol, **kwargs)
            if mask_before != (None if mask is None else bits(mask)):
                raise ValueError("RDKit Layered mutated mask")
            output = {"Bits": {"fingerprint": bits(fp), "atom_counts": counts}}
        elif kind == "Pattern":
            fp = Chem.PatternFingerprint(mol, fpSize=p["fp_size"],
                                         tautomerFingerprints=p["tautomeric"])
            output = {"Bits": {"fingerprint": bits(fp), "atom_counts": None}}
        elif kind == "Fuzzy":
            other = Chem.MolFromSmiles(row["right"]["smiles"])
            if other is None:
                raise ValueError(f"invalid right SMILES: {row['right']['id']}")
            generator = rdFingerprintGenerator.GetMorganGenerator(
                radius=p["radius"], fpSize=p["fp_size"])
            offset = (1 << 32) if p["wide"] else 0
            cls = DataStructs.ULongSparseIntVect if p["wide"] else DataStructs.UIntSparseIntVect

            def make(molecule):
                fp = cls(p["fp_size"] + offset)
                for i, v in generator.GetCountFingerprint(molecule).GetNonzeroElements().items():
                    fp[i + offset] = -v if p["signed"] and i % 2 else v
                return fp

            left, right = make(mol), make(other)
            before = (count_vector(left), count_vector(right))
            result = left | right if p["union"] else left & right
            if before != (count_vector(left), count_vector(right)):
                raise ValueError("RDKit fuzzy operation mutated an operand")
            output = {"Counts": {"left": before[0], "right": before[1],
                                 "result": count_vector(result)}}
        else:
            raise ValueError(f"unknown fingerprint profile: {kind}")
    return {"input": input_row, "output": {"Fingerprint": output}}
