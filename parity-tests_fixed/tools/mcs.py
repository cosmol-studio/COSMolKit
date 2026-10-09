"""Pinned RDKit calls only; no CK output or SMARTS normalization."""
import ctypes
import os

from rdkit import Chem
from rdkit.Chem import rdFMCS

FIELDS = {
    "store_all": "StoreAll", "maximize_bonds": "MaximizeBonds",
    "threshold": "Threshold", "timeout": "Timeout", "verbose": "Verbose",
    "initial_seed": "InitialSeed",
}
ATOM_FIELDS = {
    "match_valences": "MatchValences", "match_chiral_tag": "MatchChiralTag",
    "match_formal_charge": "MatchFormalCharge",
    "ring_matches_ring_only": "RingMatchesRingOnly",
    "complete_rings_only": "CompleteRingsOnly", "match_isotope": "MatchIsotope",
    "max_distance": "MaxDistance",
}
BOND_FIELDS = {
    "ring_matches_ring_only": "RingMatchesRingOnly",
    "complete_rings_only": "CompleteRingsOnly",
    "match_fused_rings": "MatchFusedRings",
    "match_fused_rings_strict": "MatchFusedRingsStrict", "match_stereo": "MatchStereo",
}
ATOMS = {
    "any": rdFMCS.AtomCompare.CompareAny,
    "elements": rdFMCS.AtomCompare.CompareElements,
    "isotopes": rdFMCS.AtomCompare.CompareIsotopes,
    "any_heavy_atom": rdFMCS.AtomCompare.CompareAnyHeavyAtom,
}
BONDS = {
    "any": rdFMCS.BondCompare.CompareAny,
    "order": rdFMCS.BondCompare.CompareOrder,
    "order_exact": rdFMCS.BondCompare.CompareOrderExact,
}


def parameters(values):
    params = rdFMCS.MCSParameters()
    for key, value in values.items():
        if key in FIELDS:
            setattr(params, FIELDS[key], value)
        elif key == "atom_comparator":
            params.AtomTyper = ATOMS[value]
        elif key == "bond_comparator":
            params.BondTyper = BONDS[value]
        elif key in ("atom_compare_parameters", "bond_compare_parameters"):
            atom = key == "atom_compare_parameters"
            target = params.AtomCompareParameters if atom else params.BondCompareParameters
            fields = ATOM_FIELDS if atom else BOND_FIELDS
            for name, item in value.items():
                setattr(target, fields[name], item)
        else:
            raise ValueError(f"unknown MCS parameter: {key}")
    return params


def molecule(recipe):
    if recipe["format"] == "smiles":
        params = Chem.SmilesParserParams()
        params.sanitize = recipe["sanitize"]
        params.removeHs = recipe["remove_hs"]
        mol = Chem.MolFromSmiles(recipe["text"], params)
    elif recipe["format"] == "mol":
        mol = Chem.MolFromMolBlock(recipe["text"], sanitize=recipe["sanitize"],
                                   removeHs=recipe["remove_hs"], strictParsing=True)
    else:
        raise ValueError(f"unknown MCS molecule format: {recipe['format']}")
    if mol is None:
        raise ValueError("RDKit rejected frozen MCS input")
    return mol


def mcs_case(wrapped):
    case, recipes = wrapped
    row = {key: case[key] for key in ("case_id", "inputs", "parameters", "source")}
    try:
        mols = [molecule(recipe) for recipe in recipes]
        params = parameters(case["parameters"])
        # Verbose C++ diagnostics must not corrupt the parent's JSON stdout.
        # Each worker owns its process descriptors; chemistry/options are unchanged.
        saved = os.dup(1)
        try:
            os.dup2(2, 1)
            result = rdFMCS.FindMCS(mols, params)
            ctypes.CDLL(None).fflush(None)
        finally:
            os.dup2(saved, 1)
            os.close(saved)
        query = result.queryMol
        row.update(status="ok", result={
            "atom_count": result.numAtoms, "bond_count": result.numBonds,
            "completed": not result.canceled, "smarts": result.smartsString,
            "degenerate": sorted(result.degenerateSmartsQueryMolDict),
            "query": None if query is None else {
                "atom_count": query.GetNumAtoms(), "bond_count": query.GetNumBonds(),
                "smarts": Chem.MolToSmarts(query),
            },
            "query_matches": None if query is None else [mol.HasSubstructMatch(query) for mol in mols],
        })
    except Exception as error:
        row.update(status="error", error={"type": type(error).__name__, "message": str(error)})
    return row
