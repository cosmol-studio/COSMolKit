"""One-shot development diagnostic; never an oracle or fixture generator.

Run only with the existing main interpreter and -B. Exactly eight DrawMolecule
invocations; no parsing, sanitization, coordinate generation or second drawing.
PrepareMolForDrawing returns a separate copy, observed but never drawn here.
"""

import hashlib
import json
from pathlib import Path
import struct
import sys
import tempfile
import traceback
import xml.etree.ElementTree as ET

import rdkit
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D as draw


XY = ((1.5, 0.0), (0.75, 1.299038105676658),
      (-0.75, 1.299038105676658), (-1.5, 0.0),
      (-0.75, -1.299038105676658), (0.75, -1.299038105676658))
CASES = (("benzene", False, False), ("benzene", False, True),
         ("benzene", True, False), ("benzene", True, True),
         ("pyridine", False, False), ("pyridine", False, True),
         ("pyridine", True, False), ("pyridine", True, True))
MAIN = "/home/datahouse-raid5/wjt/COSMolKit/.venv/bin/python"
ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "testdata/depiction/fixtures/legacy_0_3_0"


def bits(value):
    return struct.unpack("<Q", struct.pack("<d", value))[0]


def snapshot(mol):
    # AllProps serialization covers properties and complete modeled RDKit state;
    # typed rows make chemistry and coordinate differences independently visible.
    return {
        "binary_all_props": mol.ToBinary(Chem.PropertyPickleOptions.AllProps).hex(),
        "atoms": [{"id": a.GetIdx(), "atomic_number": a.GetAtomicNum(),
                   "aromatic": a.GetIsAromatic(), "charge": a.GetFormalCharge(),
                   "isotope": a.GetIsotope(), "map": a.GetAtomMapNum(),
                   "explicit_h": a.GetNumExplicitHs(), "implicit_h": a.GetNumImplicitHs(),
                   "no_implicit": a.GetNoImplicit(), "radicals": a.GetNumRadicalElectrons(),
                   "chiral": str(a.GetChiralTag()), "hybridization": str(a.GetHybridization()),
                   "query": a.HasQuery()} for a in mol.GetAtoms()],
        "bonds": [{"id": b.GetIdx(), "begin": b.GetBeginAtomIdx(), "end": b.GetEndAtomIdx(),
                   "type": str(b.GetBondType()), "aromatic": b.GetIsAromatic(),
                   "conjugated": b.GetIsConjugated(), "direction": str(b.GetBondDir()),
                   "stereo": str(b.GetStereo()), "stereo_atoms": list(b.GetStereoAtoms()),
                   "query": b.HasQuery()} for b in mol.GetBonds()],
        "conformers": [{"id": c.GetId(), "is_3d": c.Is3D(),
                        "xyz_bits": [[bits(v) for v in c.GetAtomPosition(i)]
                                     for i in range(mol.GetNumAtoms())]}
                       for c in mol.GetConformers()],
        "atom_rings": [list(r) for r in mol.GetRingInfo().AtomRings()],
        "bond_rings": [list(r) for r in mol.GetRingInfo().BondRings()],
    }


def require_raw(mol, label, flag):
    assert mol.GetNumAtoms() == 6 and mol.GetNumBonds() == 6
    for i, atom in enumerate(mol.GetAtoms()):
        assert atom.GetIdx() == i and atom.GetIsAromatic() is True
        assert atom.GetAtomicNum() == (7 if label == "pyridine" and i == 0 else 6)
    for i, bond in enumerate(mol.GetBonds()):
        assert (bond.GetIdx(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) == (i, i, (i + 1) % 6)
        assert bond.GetBondType() == Chem.BondType.AROMATIC
        assert bond.GetIsAromatic() == flag
    assert mol.GetNumConformers() == 1
    conf = mol.GetConformer(0)
    assert conf.GetId() == 0 and not conf.Is3D()
    assert [[bits(v) for v in conf.GetAtomPosition(i)] for i in range(6)] == [
        [bits(x), bits(y), bits(0.0)] for x, y in XY]


def build_raw(label, flag):
    mol = Chem.RWMol()
    for i in range(6):
        mol.AddAtom(Chem.Atom(7 if label == "pyridine" and i == 0 else 6))
    for i in range(6):
        mol.AddBond(i, (i + 1) % 6, Chem.BondType.AROMATIC)
    after_add = {"atom_flags": [a.GetIsAromatic() for a in mol.GetAtoms()],
                 "bond_flags": [b.GetIsAromatic() for b in mol.GetBonds()]}
    # Source RWMol convenience AddBond sets flags: reset AFTER the entire loop.
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(True)
    for bond in mol.GetBonds():
        bond.SetIsAromatic(flag)
    conf = Chem.Conformer(6)
    conf.SetId(0)
    conf.Set3D(False)
    for i, (x, y) in enumerate(XY):
        conf.SetAtomPosition(i, (x, y, 0.0))
    mol.AddConformer(conf, assignId=False)
    # Derived cache prerequisites for raw no-auto rendering; no sanitize/normalization.
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSSSR(mol)
    require_raw(mol, label, flag)
    return mol, after_add


def svg_counts(svg):
    rows = []
    for i in range(6):
        paths = [p for p in ET.fromstring(svg).iter()
                 if p.tag.endswith("}path") and f"bond-{i}" in p.attrib.get("class", "").split()]
        rows.append({"bond": i, "paths": len(paths),
                     "dashed_paths": sum("stroke-dasharray" in p.attrib.get("style", "") for p in paths)})
    return rows


def write_new(path, data):
    with path.open("x", encoding="utf-8") as handle:
        handle.write(data)


def main():
    assert sys.executable == MAIN, (sys.executable, MAIN)
    # Reading the frozen JSON verifies supplied XY provenance without rewriting it.
    for label in ("benzene", "pyridine"):
        frozen = json.loads((FIXTURES / f"{label}.input.json").read_text())
        assert frozen["coordinates"] == [list(row) for row in XY]
        assert frozen["coordinate_bits"] == [[bits(x), bits(y)] for x, y in XY]
        assert (frozen["width"], frozen["height"]) == (300, 300)
    out = Path(tempfile.mkdtemp(prefix="drawing-aromatic-boundary-", dir=ROOT / "target"))
    identity = {"interpreter": sys.executable, "rdkit_version": rdkit.__version__,
                "rdkit_module": rdkit.__file__, "drawing_module": draw.__file__,
                "drawing_module_sha256": hashlib.sha256(Path(draw.__file__).read_bytes()).hexdigest(),
                "probe_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "prepare_signature": draw.PrepareMolForDrawing.__doc__,
                "draw_signature": draw.MolDraw2DSVG.DrawMolecule.__doc__,
                "scope": "supplemental wheel, not pinned native build",
                "raw_cache_prerequisites": "UpdatePropertyCache(strict=False), GetSSSR; no sanitize",
                "b1_limit": "independent PrepareMolForDrawing counterpart, private drawer state not exposed"}
    write_new(out / "identity.json", json.dumps(identity, indent=2) + "\n")
    print("OUTPUT", out, flush=True)
    print("IDENTITY", json.dumps(identity), flush=True)
    calls = 0
    outcomes = []
    for index, (label, flag, auto) in enumerate(CASES, 1):
        raw, after_add = build_raw(label, flag)
        before = snapshot(raw)
        prepare_input = Chem.Mol(raw)
        drawing_input = Chem.Mol(raw)
        record = {"call": index, "label": label, "bond_flag": flag, "automatic": auto,
                  "after_AddBond": after_add, "B0": before, "errors": []}
        require_raw(prepare_input, label, flag)
        prep_before = snapshot(prepare_input)
        if auto:
            try:
                prepared = draw.PrepareMolForDrawing(prepare_input, kekulize=True,
                           addChiralHs=True, wedgeBonds=True, forceCoords=False, wavyBonds=False)
                record["B1"] = snapshot(prepared)
                record["B1_kind"] = "explicit independent preparation counterpart"
            except Exception:
                record["B1"] = None
                record["errors"].append({"stage": "prepare", "traceback": traceback.format_exc()})
        else:
            record["B1"] = snapshot(drawing_input)
            record["B1_kind"] = "raw no-auto input"
        record["prepare_input_preserved"] = snapshot(prepare_input) == prep_before
        drawer = draw.MolDraw2DSVG(300, 300)
        options = drawer.drawOptions()
        options.prepareMolsBeforeDrawing = auto
        record["drawing_options"] = {name: getattr(options, name) for name in (
            "prepareMolsBeforeDrawing", "centreMoleculesBeforeDrawing", "useMolBlockWedging",
            "unspecifiedStereoIsUnknown", "simplifiedStereoGroupLabel", "addAtomIndices", "addBondIndices")}
        require_raw(raw, label, flag)
        require_raw(drawing_input, label, flag)
        draw_before = snapshot(drawing_input)
        assert draw_before == before
        try:
            drawer.DrawMolecule(drawing_input, confId=0)
        except Exception:
            record["errors"].append({"stage": "draw", "traceback": traceback.format_exc()})
        finally:
            calls += 1  # Count only the invocation just returned or raised.
            record["drawing_input_preserved"] = snapshot(drawing_input) == draw_before
            record["raw_input_preserved"] = snapshot(raw) == before
            record["B0_after_draw"] = snapshot(drawing_input)
        if not any(e["stage"] == "draw" for e in record["errors"]):
            try:
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                svg_path = out / f"{index}-{label}-flag{int(flag)}-auto{int(auto)}.svg"
                write_new(svg_path, svg)
                record["B2_svg"] = str(svg_path)
                record["B2_sha256"] = hashlib.sha256(svg.encode()).hexdigest()
                record["B2_counts"] = svg_counts(svg)
            except Exception:
                record["errors"].append({"stage": "finish_or_capture", "traceback": traceback.format_exc()})
        record["input_preserved_after_capture"] = snapshot(drawing_input) == draw_before
        record["outcome"] = "error" if record["errors"] else "drawn"
        outcomes.append(record)
        write_new(out / f"{index}-outcome.json", json.dumps(record, indent=2) + "\n")
        print("OUTCOME", json.dumps(record), flush=True)
    assert calls == 8 and len(outcomes) == 8
    preserved = all(all(r[k] for k in ("prepare_input_preserved", "drawing_input_preserved",
                     "raw_input_preserved", "input_preserved_after_capture")) for r in outcomes)
    census = {"actual_drawing_calls": calls, "outcomes": len(outcomes),
              "drawn": sum(r["outcome"] == "drawn" for r in outcomes),
              "errors": sum(r["outcome"] == "error" for r in outcomes),
              "all_per_call_inputs_preserved": preserved}
    write_new(out / "census.json", json.dumps(census, indent=2) + "\n")
    print("CENSUS", json.dumps(census), flush=True)
    return 0 if preserved and census["errors"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
