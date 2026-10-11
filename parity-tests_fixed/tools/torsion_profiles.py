"""TopologicalTorsion's four public forms against current pinned RDKit."""
import json
import sys
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator
from reference import PACKAGE, Progress, parallel
from fingerprints import count_vector


def torsion_case(input_row):
    row = input_row["Fingerprint"]
    p = row["params"]["TopologicalTorsion"]
    mol = Chem.MolFromSmiles(row["case"]["smiles"])
    if mol is None:
        return {"input": input_row, "output": {"Fingerprint": {"ParseRejected": {"right": False}}}}
    generator = rdFingerprintGenerator.GetTopologicalTorsionGenerator(
        includeChirality=p["include_chirality"], torsionAtomCount=p["torsion_atom_count"],
        countSimulation=p["count_simulation"], countBounds=p["count_bounds"], fpSize=p["fp_size"])
    generator.GetOptions().onlyShortestPaths = p["only_shortest_paths"]
    generator.GetOptions().numBitsPerFeature = p["bits_per_feature"]
    kwargs = {}
    first = [0] if mol.GetNumAtoms() else []
    if p["from_first"]:
        kwargs["fromAtoms"] = first
    if p["ignore_first"]:
        kwargs["ignoreAtoms"] = first
    if p["index_plus_17"]:
        kwargs["customAtomInvariants"] = list(range(17, mol.GetNumAtoms() + 17))
    additional = None
    if p["additional_output"]:
        additional = rdFingerprintGenerator.AdditionalOutput()
        additional.AllocateAtomCounts()
        additional.AllocateAtomToBits()
        additional.AllocateBitInfoMap()
        additional.AllocateBitPaths()
        additional.AllocateAtomsPerBit()
        kwargs["additionalOutput"] = additional
    method = {"SparseCount": "GetSparseCountFingerprint", "Count": "GetCountFingerprint",
              "SparseBits": "GetSparseFingerprint", "Bits": "GetFingerprint"}[p["form"]]
    # Unexpected native errors propagate; there is no success/error coalescing.
    fp = getattr(generator, method)(mol, **kwargs)
    extra = None
    if additional is not None:
        extra = {
            "atom_counts": list(additional.GetAtomCounts()),
            "atom_to_bits": [list(v) for v in additional.GetAtomToBits()],
            "bit_info_map": {str(k): [list(v) for v in values]
                             for k, values in additional.GetBitInfoMap().items()},
            "bit_paths": {str(k): [list(v) for v in values]
                          for k, values in additional.GetBitPaths().items()},
            "atoms_per_bit": {str(k): [list(v) for v in values]
                              for k, values in additional.GetAtomsPerBit().items()},
        }
    if p["form"] in ("Count", "SparseCount"):
        observed = {"Counts": {"fingerprint": count_vector(fp), "additional_output": extra}}
    else:
        observed = {"Bits": {"fingerprint": {
            "length": fp.GetNumBits(), "on_bits": sorted(i & 0xffffffff for i in fp.GetOnBits())},
            "additional_output": extra}}
    return {"input": input_row, "output": {"Fingerprint": {"TopologicalTorsion": observed}}}


def generate(request):
    pin = json.loads((PACKAGE / "testdata/reference/rdkit.json").read_text())
    if rdBase.rdkitVersion != pin["version"]:
        raise RuntimeError(f"RDKit version {rdBase.rdkitVersion} != {pin['version']}")
    threads = request["threads"]
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be positive")
    return parallel(torsion_case, request["input"], threads, Progress(request["task"]))


if __name__ == "__main__":
    json.dump(generate(json.load(sys.stdin)), sys.stdout, allow_nan=False)
