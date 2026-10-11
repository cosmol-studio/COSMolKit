"""Historical Pattern calls against current pin; no historical outputs."""
import json
import sys
from rdkit import Chem, rdBase
from reference import PACKAGE, Progress, parallel
from fingerprints import bits, fingerprint_case
from _generate_pattern_fingerprint_golden import make_atom_counts, make_set_only_bits


def pattern_case(input_row):
    params = input_row["Fingerprint"]["params"]["Pattern"]
    options = params.get("source_options")
    if options is None:
        return fingerprint_case(input_row)
    mol = Chem.MolFromSmiles(input_row["Fingerprint"]["case"]["smiles"])
    if mol is None:
        return {"input": input_row, "output": {"Fingerprint": {"ParseRejected": {"right": False}}}}
    counts = make_atom_counts(mol, options["atom_counts"])
    mask = make_set_only_bits(params["fp_size"], options["set_only_bits"])
    before_counts = None if counts is None else list(counts)
    before_mask = None if mask is None else bits(mask)
    kwargs = {"fpSize": params["fp_size"], "tautomerFingerprints": params["tautomeric"]}
    if counts is not None:
        kwargs["atomCounts"] = counts
    if mask is not None:
        kwargs["setOnlyBits"] = mask
    # PatternFingerprints.cpp:167-170:
    # PRECONDITION(!atomCounts || atomCounts->size() >= mol.getNumAtoms(),
    #              "bad atomCounts size");
    # PRECONDITION(!setOnlyBits || setOnlyBits->getNumBits() == fpSize,
    #              "bad setOnlyBits size");
    # Neither argument is used after these preconditions in updatePatternFingerprint.
    # Wrap/MolOps.cpp copies atomCounts back unchanged. Preserve the actual
    # native calls, and fail if this audited source behavior ever changes.
    fp = Chem.PatternFingerprint(mol, **kwargs)
    if counts != before_counts or before_mask != (None if mask is None else bits(mask)):
        raise ValueError("RDKit Pattern unexpectedly mutated counts/mask")
    return {"input": input_row, "output": {"Fingerprint": {
        "Bits": {"fingerprint": bits(fp), "atom_counts": None}}}}


def generate(request):
    pin = json.loads((PACKAGE / "testdata/reference/rdkit.json").read_text())
    if rdBase.rdkitVersion != pin["version"]:
        raise RuntimeError(f"RDKit version {rdBase.rdkitVersion} != {pin['version']}")
    threads = request["threads"]
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be positive")
    return parallel(pattern_case, request["input"], threads, Progress(request["task"]))


if __name__ == "__main__":
    json.dump(generate(json.load(sys.stdin)), sys.stdout, allow_nan=False)
