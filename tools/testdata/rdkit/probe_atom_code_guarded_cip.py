"""Pinned source-native guarded code and owning-molecule effects; no test imports."""
import itertools
import json
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def molecule(tagged, computed, prepare):
    mol = Chem.RWMol()
    for symbol in ["C", "F", "Cl", "Br", "I"]:
        mol.AddAtom(Chem.Atom(symbol))
    for endpoint in range(1, 5):
        mol.AddBond(0, endpoint, Chem.BondType.SINGLE)
    if tagged:
        mol.GetAtomWithIdx(0).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    mol.SetProp("user", "preserved")
    if computed is not None:
        mol.SetBoolProp("_CIPComputed", computed)
    if prepare:
        mol.UpdatePropertyCache(strict=False)
    return mol

def effects(mol):
    return {"marker": mol.GetProp("_CIPComputed") if mol.HasProp("_CIPComputed") else None,
            "user": mol.GetProp("user"),
            "atom_labels": [a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else None for a in mol.GetAtoms()],
            "marker_in_computed_list": "_CIPComputed" in tuple(mol.GetPropNames(includePrivate=True, includeComputed=True)),
            "marker_in_ordinary_list": "_CIPComputed" in tuple(mol.GetPropNames(includePrivate=True, includeComputed=False))}

for legacy, include, tagged, computed, subtract in itertools.product([False, True], [False, True], [False, True], [None, False, True], [0, 4, (1 << 32) - 1]):
    Chem.SetUseLegacyStereoPerception(legacy)
    mol = molecule(tagged, computed, True)
    record = {"legacy": legacy, "include": include, "tagged": tagged, "computed": computed, "subtract": subtract, "before": effects(mol)}
    try:
        record["code"] = rdMolDescriptors.GetAtomPairAtomCode(mol.GetAtomWithIdx(0), subtract, include)
    except Exception as error:
        record["error"] = {"type": type(error).__name__, "message": str(error)}
    record["after"] = effects(mol)
    print(json.dumps(record, sort_keys=True))

Chem.SetUseLegacyStereoPerception(False)
mol = molecule(True, None, False)
record = {"case": "missing-cache-before-cip", "before": effects(mol)}
try:
    record["code"] = rdMolDescriptors.GetAtomPairAtomCode(mol.GetAtomWithIdx(0), (1 << 32) - 1, True)
except Exception as error:
    record["error"] = {"type": type(error).__name__, "message": str(error)}
record["after"] = effects(mol)
print(json.dumps(record, sort_keys=True))
