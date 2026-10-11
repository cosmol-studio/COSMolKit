"""Pinned native public search and tautomer calls, not copied runner machinery."""
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize as T
from chemistry import state
import json
import re
import struct
from pathlib import Path

@lru_cache(maxsize=256)
def query(text):
    q=Chem.MolFromSmarts(text)
    if q is None: raise ValueError(f"invalid query {text!r}")
    return q

def search(mol,p):
    q=query(p["smarts"])
    params=Chem.SubstructMatchParameters()
    for public,native in (
        ("max_matches","maxMatches"),("uniquify","uniquify"),("use_chirality","useChirality"),
        ("use_enhanced_stereo","useEnhancedStereo"),("specified_stereo_query_matches_unspecified","specifiedStereoQueryMatchesUnspecified"),
        ("use_query_query_matches","useQueryQueryMatches"),("recursion_possible","recursionPossible"),
        ("max_recursive_matches","maxRecursiveMatches"),("num_threads","numThreads"),
        ("aromatic_matches_conjugated","aromaticMatchesConjugated"),("aromatic_matches_single_or_double","aromaticMatchesSingleOrDouble"),
        ("atom_properties","atomProperties"),("bond_properties","bondProperties"),("use_generic_matchers","useGenericMatchers"),
    ): setattr(params,native,p[public])
    rows=mol.GetSubstructMatches(q,params)
    bonds=[[mol.GetBondBetweenAtoms(row[b.GetBeginAtomIdx()],row[b.GetEndAtomIdx()]).GetIdx() for b in q.GetBonds()] for row in rows]
    first=mol.GetSubstructMatch(q,params)
    return [[list(r) for r in rows],bonds,list(first) if first else None,
            mol.HasSubstructMatch(q,params),[list(r) for r in mol.GetSubstructMatches(q)] if p["compiled_default"] else None]

def score(m,e):
    return [T.ScoreRings(m),T.ScoreSubstructs(m),T.ScoreHeteroHs(m),e.ScoreTautomer(m)]

def tautomer_state(m):
    # Same explicit fixed-width graph fields as Rust, including CIP and the
    # complete stereo-reference vector. This is lossless integer packing.
    raw=bytearray(struct.pack("<I",m.GetNumAtoms()))
    for a in m.GetAtoms():
        raw.extend(struct.pack("<BhHBBBBBBII",a.GetAtomicNum(),a.GetFormalCharge(),a.GetIsotope(),
                               a.GetNumExplicitHs(),a.GetNoImplicit(),a.GetIsAromatic(),int(a.GetChiralTag()),
                               int(a.GetHybridization()),a.GetNumRadicalElectrons(),a.GetAtomMapNum(),
                               a.GetIntProp("_chiralPermutation") if a.HasProp("_chiralPermutation") else 0))
        cip=a.GetProp("_CIPCode").encode() if a.HasProp("_CIPCode") else None
        raw.extend(struct.pack("<I",len(cip) if cip is not None else 0xffffffff))
        if cip is not None:raw.extend(cip)
    raw.extend(struct.pack("<I",m.GetNumBonds()))
    for b in m.GetBonds():
        refs=b.GetStereoAtoms()
        raw.extend(struct.pack("<IIBBBBBI",b.GetBeginAtomIdx(),b.GetEndAtomIdx(),int(b.GetBondType()),
                               b.GetIsAromatic(),b.GetIsConjugated(),int(b.GetBondDir()),int(b.GetStereo()),len(refs)))
        for a in refs:raw.extend(struct.pack("<I",a))
    return bytes(raw)

@lru_cache(maxsize=1)
def v1_transforms():
    # Tautomer.h::getV1TautomerEnumerator uses defaultTautomerTransformsv1.
    # Read that pinned table, not a separately maintained copy of its chemistry.
    path=Path(__file__).resolve().parents[2]/"third_party/rdkit/Code/GraphMol/MolStandardize/TautomerCatalog/tautomerTransforms.v1.in"
    source=path.read_text()
    expressions=re.findall(r'std::string\(\s*((?:"(?:[^"\\]|\\.)*"\s*)+)\)',source)
    if len(expressions)!=source.count("std::make_tuple")*4:
        raise RuntimeError("unrecognized pinned V1 transform table layout")
    strings=["".join(json.loads(token) for token in re.findall(r'"(?:[^"\\]|\\.)*"',expression)) for expression in expressions]
    return [dict(zip(("name","smarts","bonds","charges"),strings[i:i+4])) for i in range(0,len(strings),4)]

def tautomer(mol,p):
    params=T.CleanupParameters()
    for field,native in (("max_tautomers","maxTautomers"),("max_transforms","maxTransforms"),
                         ("remove_sp3_stereo","tautomerRemoveSp3Stereo"),("remove_bond_stereo","tautomerRemoveBondStereo"),
                         ("remove_isotopic_hydrogens","tautomerRemoveIsotopicHs"),("reassign_stereo","tautomerReassignStereo")):
        setattr(params,native,p[field])
    if p["catalog"]=="v1":T.UpdateParamsFromJSON(params,json.dumps({"tautomerTransformData":v1_transforms()}))
    e=T.TautomerEnumerator(params)
    result=e.Enumerate(Chem.Mol(mol))
    entries=[[key.encode("utf-8"),Chem.MolToSmiles(m).encode("utf-8"),tautomer_state(m),score(m,e)] for key,m in zip(result.smiles,result.tautomers)]
    chosen=e.PickCanonical(result)
    direct=e.Canonicalize(Chem.Mol(mol))
    def selected(m):
        text=Chem.MolToSmiles(m).encode("utf-8")
        return [text,tautomer_state(m),score(m,e)]
    return [int(result.status),list(result.modifiedAtoms),list(result.modifiedBonds),entries,selected(chosen),selected(direct)]
