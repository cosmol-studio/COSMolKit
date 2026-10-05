"""Pinned source-native Morgan wrapper/provider and JSON boundary observations.

No ordinary tests imported; raw source values retained beside width transport.
"""
import argparse
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as fg

SMILES=("C1CCCCC1","c1ccccc1","C1CCOC1","CCO","")
PROFILES=("bare_default","bare_false","bare_true","explicit_false","explicit_true")
OUTPUTS=(("sparse_count","GetSparseCountFingerprint"),("sparse_bits","GetSparseFingerprint"),("count","GetCountFingerprint"),("bits","GetFingerprint"))

def generator(profile):
    if profile=="bare_default": return fg.GetMorganGenerator()
    if profile=="bare_false": return fg.GetMorganGenerator(includeRingMembership=False)
    if profile=="bare_true": return fg.GetMorganGenerator(includeRingMembership=True)
    if profile=="explicit_false": return fg.GetMorganGenerator(atomInvariantsGenerator=fg.GetMorganAtomInvGen(False))
    if profile=="explicit_true": return fg.GetMorganGenerator(atomInvariantsGenerator=fg.GetMorganAtomInvGen(True))
    raise ValueError(profile)

def output(owner,mol,method):
    ao=fg.AdditionalOutput()
    ao.AllocateAtomCounts();ao.AllocateAtomToBits();ao.AllocateBitInfoMap();ao.AllocateBitPaths()
    try:
        value=getattr(owner,method)(mol,additionalOutput=ao)
        result={"size":int(value.GetLength()) if hasattr(value,"GetLength") else int(value.GetNumBits())}
        if hasattr(value,"GetNonzeroElements"):
            result["nonzero"]=sorted([int(key),int(count)] for key,count in value.GetNonzeroElements().items())
        else:
            result["native_bits"]=[int(bit) for bit in value.GetOnBits()]
            result["bit_pattern_u32"]=sorted(bit & ((1<<32)-1) for bit in result["native_bits"])
        result["additional_output"]={"atom_counts":list(ao.GetAtomCounts()),"atom_to_bits":[list(bits) for bits in ao.GetAtomToBits()],"bit_info_map":sorted([int(key),[list(pair) for pair in pairs]] for key,pairs in ao.GetBitInfoMap().items()),"bit_paths":sorted([int(key),[list(path) for path in paths]] for key,paths in ao.GetBitPaths().items())}
        return result
    except Exception as error:
        return {"error":{"type":type(error).__name__,"message":str(error)}}

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,required=True);args=parser.parse_args();records=[]
    Chem.SetUseLegacyStereoPerception(True)
    for smiles in SMILES:
        mol=Chem.MolFromSmiles(smiles)
        for profile in PROFILES:
            owner=generator(profile)
            record={"case":"ring-wrapper-provider","smiles":smiles,"profile":profile,"json":json.loads(owner.ToJSON()),"info":owner.GetInfoString()}
            record["outputs"]={label:output(owner,mol,method) for label,method in OUTPUTS}
            records.append(record)
    mol=Chem.MolFromSmiles("CCO")
    feature=fg.GetMorganGenerator(radius=0,atomInvariantsGenerator=fg.GetMorganFeatureAtomInvGen())
    value=json.loads(feature.ToJSON());value["atomInvariantsGenerator"]["patternSMARTS"]=["[C]","[","[O]"]
    restored=fg.FingerprintGeneratorFromJSON(json.dumps(value))
    records.append({"case":"feature-invalid-source-null-pattern-skip","smiles":"CCO","input_json":value,"json":json.loads(restored.ToJSON()),"info":restored.GetInfoString(),"outputs":{label:output(restored,mol,method) for label,method in OUTPUTS}})
    with args.output.open("x") as stream:
        for record in records:stream.write(json.dumps(record,sort_keys=True)+"\n")
    print(json.dumps({"records":len(records),"output_calls":len(records)*len(OUTPUTS),"outputs_file":str(args.output),"test_count":0}))

if __name__=="__main__": main()
