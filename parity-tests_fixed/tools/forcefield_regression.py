"""Fixed force-field counterexamples: native preparation, never CK expectations."""
from rdkit import Chem
from rdkit.Chem import AllChem
from _generate_forcefield_params_golden import (
    forcefield_initial_energy_result,
    forcefield_optimized_result,
    forcefield_multi_optimized_result,
)
from _generate_mmff_builtin_golden import mmff_result


def optimizer_case(case):
    # Preserve the original fixed cases' seed, implicit-H graph and CXSMILES
    # coordinate roundtrip. Only the user-requested step budget changes to 2.
    mol = Chem.MolFromSmiles(case["smiles"])
    if mol is None:
        raise ValueError(f"{case['case_id']}: MolFromSmiles failed")
    status = AllChem.EmbedMolecule(mol, randomSeed=61453, useRandomCoords=True)
    if status != 0:
        raise ValueError(f"{case['case_id']}: EmbedMolecule returned {status}")
    cxsmiles = Chem.MolToCXSmiles(mol, isomericSmiles=True)
    mol = Chem.MolFromSmiles(cxsmiles)
    if mol is None:
        raise ValueError(f"{case['case_id']}: CXSMILES roundtrip failed")
    if case["forcefield"] == "mmff":
        def factory(m):
            return AllChem.MMFFGetMoleculeForceField(
                m, AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94"),
                nonBondedThresh=100.0, confId=0, ignoreInterfragInteractions=True)

        def optimize_conformers(m):
            return AllChem.MMFFOptimizeMoleculeConfs(
                m, numThreads=1, maxIters=2, mmffVariant="MMFF94",
                nonBondedThresh=100.0, ignoreInterfragInteractions=True)
        has_all = bool(AllChem.MMFFHasAllMoleculeParams(Chem.Mol(mol)))
    elif case["forcefield"] == "uff":
        def factory(m):
            return AllChem.UFFGetMoleculeForceField(
                m, vdwThresh=100.0, confId=0, ignoreInterfragInteractions=True)

        def optimize_conformers(m):
            return AllChem.UFFOptimizeMoleculeConfs(
                m, numThreads=1, maxIters=2, vdwThresh=100.0,
                ignoreInterfragInteractions=True)
        has_all = bool(AllChem.UFFHasAllMoleculeParams(Chem.Mol(mol)))
    else:
        raise ValueError(f"unknown force field: {case['forcefield']}")
    output = {
        "cxsmiles": cxsmiles,
        "coords": [[float(p[axis]) for axis in range(3)]
                   for p in mol.GetConformer().GetPositions()],
        "has_all": has_all,
        "initial": forcefield_initial_energy_result(factory, mol),
        "single": forcefield_optimized_result(factory, mol, 2),
        "multi": forcefield_multi_optimized_result(optimize_conformers, mol, 2),
    }
    return {"case_id": case["case_id"], "input": case, "output": output}


def builtin_case(case):
    return {"case_id": case["case_id"], "input": case,
            "output": mmff_result(case["smiles"], case["variant"])}
