"""The original generator's complete result/state schema, through public Python APIs."""
import base64
import json
import os
from pathlib import Path
import cosmolkit as ck
import pytest

ROOT = Path(os.environ["ALIGNMENT_ORACLE_ROOT"])
INPUTS = {(r["profile"],r["case_id"],r["call_index"]):r["inputs"] for r in map(json.loads,(ROOT.parent/"canonical-input-pickles.jsonl").read_text().splitlines())}

def state(m):
    return [dict(id=c.id(),is_3d=c.is_3d(),coordinates=c.coordinates()) for c in m.conformers_3d()]

def compare(actual, expected):
    if isinstance(expected,dict):
        assert actual.keys()==expected.keys()
        for k in expected: compare(actual[k],expected[k])
    elif isinstance(expected,list):
        assert len(actual)==len(expected)
        for a,b in zip(actual,expected):compare(a,b)
    elif isinstance(expected,float): assert actual==pytest.approx(expected,abs=1e-8,rel=0)
    else: assert actual==expected

def amap(values):return [ck.AlignmentAtomMap(*v) for v in values]
def parameters(typ,values):
    values=dict(values)
    values.pop("operation",None)  # generator record routing metadata, not a chemistry parameter
    if "atom_map" in values:values["atom_map"]=amap(values["atom_map"])
    if "atom_maps" in values:values["atom_maps"]=[amap(x) for x in values["atom_maps"]]
    return typ(**values)
def result(r):return dict(rmsd=r.rmsd(),transform=r.transform().matrix(),atom_map=[[m.probe_atom,m.reference_atom] for m in r.atom_map()])

def profile(name,count):
    rows=[json.loads(x) for x in (ROOT/name/"molalign.jsonl").read_text().splitlines()]
    assert len(rows)==count
    for r in rows:
        assert r["schema_version"]==1
        assert (r.get("error_type") is None)==(r["status"]=="ok")
        assert (r.get("error_message") is None)==(r["status"]=="ok")
        if r["operation"]=="input_parse":
            assert r["status"]=="error" and r["error_kind"]=="input_parse_error"
            with pytest.raises(ValueError):ck.Molecule.from_smiles(r["source"]["input_smiles"])
            continue
        pickles=INPUTS[(name,r["case_id"],r["call_index"])]
        molecules={k:ck.Molecule.from_binary(bytes(v)) for k,v in pickles.items()}
        for k,m in molecules.items():compare(state(m),r["before"][k])
        op=r["operation"];params=r["parameters"];actual=None
        try:
            if op=="alignment_transform":actual=result(molecules["probe"].alignment_transform_to(molecules["reference"],parameters(ck.AlignmentParameters,params)))
            elif op=="best_alignment":
                p=parameters(ck.BestAlignmentParameters,params); actual=result(molecules["probe"].best_alignment_to(molecules["reference"],p));compare(molecules["probe"].best_rmsd_to(molecules["reference"],p),actual["rmsd"])
            elif op=="coordinate_rmsd":actual=dict(rmsd=molecules["probe"].coordinate_rmsd_to(molecules["reference"],parameters(ck.CoordinateRmsdParameters,params)))
            elif op=="align_to":
                p=parameters(ck.AlignmentParameters,params); aligned,report=molecules["probe"].with_alignment_to(molecules["reference"],p);actual=result(report);compare(state(aligned),r["after"]["probe"]);compare(state(molecules["probe"]),r["before"]["probe"])
                report=molecules["probe"].align_to_(molecules["reference"],p);compare(result(report),actual)
            elif op=="all_conformer_best_rms":
                values=molecules["molecule"].all_conformer_best_rmsds(parameters(ck.AllConformerRmsdParameters,params));actual=dict(rmsds=[v.rmsd() for v in values],conformer_pairs=[[v.probe_conformer_id(),v.reference_conformer_id()] for v in values])
            elif op=="align_conformers":
                p=parameters(ck.ConformerAlignmentParameters,params);aligned,report=molecules["molecule"].with_aligned_conformers(p);actual=dict(rmsds=report.rmsds());compare(state(aligned),r["after"]["molecule"]);compare(state(molecules["molecule"]),r["before"]["molecule"]);report=molecules["molecule"].align_conformers_(p);compare(report.rmsds(),actual["rmsds"])
            else:raise AssertionError(op)
        except ValueError as e:
            assert r["status"]=="error",(r["case_id"],op,str(e))
            messages={"conformer_not_found":"conformer id", "weight_count_mismatch":"weights", "no_substructure_match":"no substructure match"}
            assert messages[r["error_kind"]] in str(e)
        else:
            assert r["status"]=="ok",(r["case_id"],op)
            compare(actual,r["result"])
        for k,m in molecules.items():compare(state(m),r["after"][k])

@pytest.mark.parametrize("name,count",[("molalign_focused",14),("smiles_small",152),("smiles_5000",5000)])
def test_original_molalign_complete_observables(name,count):profile(name,count)
