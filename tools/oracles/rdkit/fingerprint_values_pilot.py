"""Reference adapter only; Rust owns task selection, inputs and comparison."""
import json
import sys
from rdkit import DataStructs, rdBase

expected_version = sys.argv[1]
if rdBase.rdkitVersion != expected_version:
    raise RuntimeError(f"RDKit {rdBase.rdkitVersion} != {expected_version}")

results = []
for row in json.load(sys.stdin):
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
    results.append({"input": row, "output": {
        "length": result.GetLength(),
        "entries": sorted(result.GetNonzeroElements().items()),
    }})
json.dump(results, sys.stdout, sort_keys=True)
