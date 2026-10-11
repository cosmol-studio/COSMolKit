"""Core public results in original atom/bond order, without private serializers."""
from rdkit import Chem
from rdkit.Chem import rdMolHash

TASKS = {"smiles_read","sanitize","kekulize","add_hs","remove_hs","potential_stereo",
         "distance_matrix","fragments","largest_fragment","murcko_scaffold","net_scaffold","murcko_decompose"}

def state(mol):
    return [[[a.GetAtomicNum(),a.GetFormalCharge(),a.GetIsotope(),a.GetNumExplicitHs(),
              a.GetNoImplicit(),a.GetIsAromatic(),int(a.GetChiralTag()),int(a.GetHybridization()),
              a.GetNumRadicalElectrons(),a.GetAtomMapNum(),
              a.GetIntProp("_chiralPermutation") if a.HasProp("_chiralPermutation") else 0]
             for a in mol.GetAtoms()],
            [[b.GetBeginAtomIdx(),b.GetEndAtomIdx(),int(b.GetBondType()),b.GetIsAromatic(),
              b.GetIsConjugated(),int(b.GetBondDir()),int(b.GetStereo()),list(b.GetStereoAtoms())]
             for b in mol.GetBonds()]]

def parse(text,task,p):
    options=Chem.SmilesParserParams()
    options.sanitize=p["sanitize"] if task=="smiles_read" else task!="sanitize"
    options.removeHs=p["remove_hs"] if task=="smiles_read" else task not in ("sanitize","remove_hs")
    return Chem.MolFromSmiles(text,options)

def chemistry(mol,task,p):
    m=Chem.Mol(mol)  # Every profile receives an independent mutable native value.
    if task=="smiles_read": return state(m)
    if task=="sanitize":
        Chem.SanitizeMol(m)
        return state(m)
    if task=="kekulize":
        Chem.Kekulize(m,clearAromaticFlags=p["clear_aromatic_flags"])
        return state(m)
    if task=="add_hs":
        params=Chem.AddHsParameters()
        params.explicitOnly=p["explicit_only"]
        return state(Chem.AddHs(m,params,[0] if p["first_atom_only"] else None))
    if task=="remove_hs":
        params=Chem.RemoveHsParameters()
        params.removeIsotopes=p["remove_isotopes"]
        params.removeMapped=p["remove_mapped"]
        params.removeHydrides=p["remove_hydrides"]
        return state(Chem.RemoveHs(m,params,p["sanitize"]))
    if task=="distance_matrix":
        matrix=Chem.GetDistanceMatrix(m,useBO=p["use_bond_order"],useAtomWts=p["use_atom_weights"],force=True)
        return [len(matrix),[float(x) for x in matrix.flat]]
    if task=="potential_stereo":
        before=Chem.GetAllowNontetrahedralChirality()
        Chem.SetAllowNontetrahedralChirality(p["allow_nontetrahedral"])
        try:
            rows=Chem.FindPotentialStereo(m,cleanIt=p["clean"],flagPossible=p["flag_possible"])
            result=[[int(s.type),int(s.specified),s.centeredOn,int(s.descriptor),s.permutation,list(s.controllingAtoms)] for s in rows]
            return [result,state(m) if p["clean"] else None]
        finally: Chem.SetAllowNontetrahedralChirality(before)
    if task=="fragments": return [state(x) for x in Chem.GetMolFrags(m,asMols=True,sanitizeFrags=True)]
    if task=="largest_fragment":
        rows=Chem.GetMolFrags(m,asMols=True,sanitizeFrags=True)
        if not rows: raise ValueError("empty fragments")
        return state(max(enumerate(rows),key=lambda pair:(pair[1].GetNumAtoms(),pair[0]))[1])
    if task=="murcko_decompose": return state(Chem.MurckoDecompose(m))
    if task in ("murcko_scaffold","net_scaffold"):
        fn=rdMolHash.HashFunction.MurckoScaffold if task=="murcko_scaffold" else rdMolHash.HashFunction.ExtendedMurcko
        return rdMolHash.MolHash(m,fn).encode("utf-8")
    raise ValueError(f"unknown core task {task}")
