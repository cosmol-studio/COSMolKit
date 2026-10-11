"""Historical Morgan native calls; ordinary reference dependencies stay intact."""
import json
import sys
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator as F
from reference import PACKAGE, Progress, parallel
from fingerprint_values_pilot import molecular, morgan_additional_output


def historical(mol, recipe):
    options = dict(radius=2, fpSize=2048)
    kwargs = {}
    if recipe == "chiral":
        options["includeChirality"] = True
    elif recipe == "radius3":
        options["radius"] = 3
    elif recipe == "no_bond_types":
        options.update(fpSize=1024, useBondTypes=False)
    elif recipe == "count_simulation":
        options.update(countSimulation=True, countBounds=[1, 2, 4, 8])
    elif recipe == "redundant_environments":
        options["includeRedundantEnvironments"] = True
    elif recipe == "feature_invariants":
        options["atomInvariantsGenerator"] = F.GetMorganFeatureAtomInvGen()
    elif recipe == "connectivity_without_ring_membership":
        options["atomInvariantsGenerator"] = F.GetMorganAtomInvGen(False)
    elif recipe == "chiral_bond_invariants":
        options["bondInvariantsGenerator"] = F.GetMorganBondInvGen(True, True)
    elif recipe == "from_first":
        kwargs["fromAtoms"] = [0]
    elif recipe == "ignore_first":
        kwargs["ignoreAtoms"] = [0]
    elif recipe == "custom_atom_invariants":
        kwargs["customAtomInvariants"] = list(range(1, mol.GetNumAtoms() + 1))
    elif recipe == "only_nonzero_custom_invariants":
        options["onlyNonzeroInvariants"] = True
        kwargs["customAtomInvariants"] = [0 if i % 2 == 0 else i + 1 for i in range(mol.GetNumAtoms())]
    elif recipe == "custom_bond_invariants":
        kwargs["customBondInvariants"] = list(range(7, mol.GetNumBonds() + 7))
    elif recipe not in ("default", "two_bits_per_feature"):
        raise ValueError(f"unknown historical Morgan recipe:{recipe}")
    if recipe in ("from_first", "ignore_first", "custom_atom_invariants", "only_nonzero_custom_invariants") and not mol.GetNumAtoms():
        raise ValueError("historical Morgan atom-conditional recipe without atoms")
    if recipe == "custom_bond_invariants" and not mol.GetNumBonds():
        raise ValueError("historical Morgan bond-conditional recipe without bonds")
    generator = F.GetMorganGenerator(**options)
    if recipe == "two_bits_per_feature":
        generator.GetOptions().numBitsPerFeature = 2
    extra = F.AdditionalOutput()
    extra.AllocateAtomCounts()
    extra.AllocateAtomToBits()
    extra.AllocateBitInfoMap()
    extra.AllocateBitPaths()
    extra.AllocateAtomsPerBit()
    fp = generator.GetFingerprint(mol, additionalOutput=extra, **kwargs)
    return {"atom_count": mol.GetNumAtoms(), "bond_count": mol.GetNumBonds(),
            "fingerprint": {"length": fp.GetNumBits(), "on_bits": list(fp.GetOnBits()),
                            "additional_output": morgan_additional_output(extra)}}


def morgan_case(input_row):
    row = input_row["Molecular"]
    p = row["profile"]
    if not isinstance(p, dict) or "MorganHistorical" not in p:
        return {"input": input_row, "output": {"Molecular": molecular(row)}}
    mol = Chem.MolFromSmiles(row["case"]["smiles"])
    if mol is None:
        return {"input": input_row, "output": {"Molecular": {
            "ParseRejected": {"detail": "RDKit MolFromSmiles returned None"}}}}
    return {"input": input_row, "output": {"Molecular": {
        "MorganHistorical": historical(mol, p["MorganHistorical"]["recipe"])}}}


def generate(request):
    pin = json.loads((PACKAGE / "testdata/reference/rdkit.json").read_text())
    if rdBase.rdkitVersion != pin["version"]:
        raise RuntimeError(f"RDKit version {rdBase.rdkitVersion} != {pin['version']}")
    threads = request["threads"]
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be positive")
    return parallel(morgan_case, request["input"], threads, Progress(request["task"]))


if __name__ == "__main__":
    json.dump(generate(json.load(sys.stdin)), sys.stdout, allow_nan=False)
