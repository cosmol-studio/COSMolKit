"""AtomPair's four public output forms against the current pinned RDKit."""
import json
import sys
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator
from reference import PACKAGE, Progress, parallel
from fingerprints import count_vector


def atom_pair_case(input_row):
    row = input_row["Fingerprint"]
    p = row["params"]["AtomPair"]
    mol = Chem.MolFromSmiles(row["case"]["smiles"])
    if mol is None:
        return {"input": input_row, "output": {"Fingerprint": {"ParseRejected": {"right": False}}}}
    generator = rdFingerprintGenerator.GetAtomPairGenerator(
        minDistance=p["min_distance"], maxDistance=p["max_distance"],
        includeChirality=p["include_chirality"], use2D=p["use_2d"],
        countSimulation=p["count_simulation"], countBounds=p["count_bounds"],
        fpSize=p["fp_size"])
    generator.GetOptions().numBitsPerFeature = p["bits_per_feature"]
    kwargs = {}
    first = [0] if mol.GetNumAtoms() else []
    if p["from_first"]:
        kwargs["fromAtoms"] = first
    if p["ignore_first"]:
        kwargs["ignoreAtoms"] = first
    if p["index_plus_11"]:
        kwargs["customAtomInvariants"] = list(range(11, mol.GetNumAtoms() + 11))
    additional = None
    if p["additional_output"]:
        additional = rdFingerprintGenerator.AdditionalOutput()
        additional.AllocateAtomCounts()
        additional.AllocateAtomToBits()
        additional.AllocateBitInfoMap()
        additional.AllocateAtomsPerBit()
        kwargs["additionalOutput"] = additional
    method = {"SparseCount": "GetSparseCountFingerprint", "Count": "GetCountFingerprint",
              "SparseBits": "GetSparseFingerprint", "Bits": "GetFingerprint"}[p["form"]]
    try:
        fp = getattr(generator, method)(mol, **kwargs)
    except IndexError as error:
        # FingerprintGenerator.cpp uses unmodded random extra bit IDs when
        # getSparseCountFingerprint() passes fpSize=0; SparseIntVect rejects
        # IDs outside AtomPairGenerator::getResultSize(). Keep the exact index.
        text = str(error)
        size = 1 << (23 + (4 if p["include_chirality"] else 0))
        if (p["form"] != "SparseCount" or p["bits_per_feature"] <= 1
                or not text.isdecimal() or int(text) < size):
            raise
        observed = {"IndexError": {"index": int(text), "size": size}}
    else:
        extra = None
        if additional is not None:
            extra = {
                "atom_counts": list(additional.GetAtomCounts()),
                "atom_to_bits": [list(v) for v in additional.GetAtomToBits()],
                "bit_info_map": {str(k): [list(v) for v in values]
                                 for k, values in additional.GetBitInfoMap().items()},
                "atoms_per_bit": {str(k): [list(v) for v in values]
                                  for k, values in additional.GetAtomsPerBit().items()},
            }
        if p["form"] in ("Count", "SparseCount"):
            observed = {"Counts": {"fingerprint": count_vector(fp), "additional_output": extra}}
        else:
            observed = {"Bits": {"fingerprint": {
                "length": fp.GetNumBits(), "on_bits": sorted(i & 0xffffffff for i in fp.GetOnBits())},
                "additional_output": extra}}
    return {"input": input_row, "output": {"Fingerprint": {"AtomPair": observed}}}


def generate(request):
    pin = json.loads((PACKAGE / "testdata/reference/rdkit.json").read_text())
    if rdBase.rdkitVersion != pin["version"]:
        raise RuntimeError(f"RDKit version {rdBase.rdkitVersion} != {pin['version']}")
    threads = request["threads"]
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be positive")
    return parallel(atom_pair_case, request["input"], threads, Progress(request["task"]))


if __name__ == "__main__":
    json.dump(generate(json.load(sys.stdin)), sys.stdout, allow_nan=False)
