"""Capture current-pin XY bits for the existing 48-cell ordinary fixture.

Run manually from the repository root with .venv/bin/python and apply the
printed patch. This is not called by Cargo tests. Historical input/working-state
captures remain unchanged; coordinate_reference records the numerical pin.
Only RDKit computes expectations; COSMolKit is never imported.
"""
from pathlib import Path
import difflib
import hashlib
import json
import re
import struct
import sys

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor


ROOT = Path(__file__).resolve().parents[3]
FIXTURE = "testdata/depict_2d/expected/rdkit/prepared_property_presence.tsv"
PROVENANCE = FIXTURE.removesuffix(".tsv") + ".source.json"
CASES = "crates/cosmolkit/tests/fixtures/d2_prepared_property.rs"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def capture():
    pin = json.loads((ROOT / "testdata/reference/rdkit.json").read_text())
    if rdBase.rdkitVersion != pin["version"]:
        raise RuntimeError(f"Expected RDKit {pin['version']}, found {rdBase.rdkitVersion}")
    original = (ROOT / FIXTURE).read_text()
    provenance_text = (ROOT / PROVENANCE).read_text()
    provenance = json.loads(provenance_text)
    if digest(original.encode()) != provenance["fixture_sha256"]:
        raise RuntimeError("Fixture identity differs from its provenance")
    case_source = (ROOT / CASES).read_text()
    # Read only the three shared input literals, never CK-generated results.
    cases = dict(re.findall(r'\(\s*(\d+),\s*"([^"]+)"', case_source.split("fn unhex", 1)[0]))
    if set(map(int, cases)) != set(provenance["cases"]):
        raise RuntimeError("Fixed input census changed")
    header, *sections = original.split("BLOCK\t")
    blocks = []
    for section in sections:
        key, body = section.split("\n", 1)
        content, rest = body.split("END\n", 1)
        if rest.strip():
            raise RuntimeError("Unexpected fixture block trailer")
        if not content.startswith("XY\t"):
            blocks.append((key, content))
    cells = [row.split("\t") for row in header.splitlines()]
    if len(cells) != 48 or any(len(row) != 7 or row[0] != "CELL" for row in cells):
        raise RuntimeError("Expected the complete 48-cell matrix")
    Chem.SetUseLegacyStereoPerception(True)
    rdDepictor.SetPreferCoordGen(False)
    numeric = {}
    coordinate_rows = 0
    for cell in cells:
        match = re.fullmatch(r"line:(\d+):P([0-7]):O([01])", cell[1])
        if not match:
            raise RuntimeError(f"Unexpected cell {cell[1]}")
        line, state, orientation = match.groups()
        mol = Chem.MolFromSmiles(cases[line])
        if mol is None:
            raise RuntimeError(f"Native parse failed: {line}")
        mol.ClearProp("_StereochemDone")
        state = int(state)
        if state in (1, 2, 3):
            mol.SetIntProp("_StereochemDone", int(state == 2), state != 3)
        elif state in (4, 5, 6):
            mol.SetProp("_StereochemDone", {4: "", 5: "0", 6: "false"}[state], True)
        elif state == 7:
            mol.SetIntProp("_stereochemDone", 1)
        conf = rdDepictor.Compute2DCoords(
            mol, canonOrient=bool(int(orientation)), clearConfs=True,
            nFlipsPerSample=0, nSample=0, sampleSeed=0, permuteDeg4Nodes=False,
            forceRDKit=False, useRingTemplates=False,
        )
        rows = []
        for index in range(mol.GetNumAtoms()):
            point = mol.GetConformer(conf).GetAtomPosition(index)
            bits = [struct.unpack("<Q", struct.pack("<d", point[axis]))[0] for axis in range(2)]
            rows.append(f"XY\t{index}\t{bits[0]}\t{bits[1]}\n")
        content = "".join(rows)
        key = digest(content.encode())
        numeric[key] = content
        cell[6] = key
        coordinate_rows += len(rows)
    updated = "".join("\t".join(cell) + "\n" for cell in cells)
    updated += "".join(f"BLOCK\t{key}\n{content}END\n" for key, content in blocks + list(numeric.items()))
    provenance["fixture_sha256"] = digest(updated.encode())
    provenance["coordinate_reference"] = {
        "version": pin["version"], "commit": pin["source_revision"],
        "generator": str(Path(__file__).resolve().relative_to(ROOT)),
        "generator_sha256": digest(Path(__file__).read_bytes()),
        "input_cases": CASES, "input_selection": "CASES input literal table",
        "input_cases_sha256": digest(json.dumps(cases, sort_keys=True).encode()),
        "fixed_behavior_cells": len(cells), "coordinate_rows": coordinate_rows,
        "unique_coordinate_blocks": len(numeric),
        "source_parameters": "legacy=1; preferCoordGen=0; canonOrient=O; clearConfs=1; coordMap=empty; nFlipsPerSample=0; nSamples=0; sampleSeed=0; permuteDeg4Nodes=0; forceRDKit=0; useRingTemplates=0; usePathAngleExpansion=1 (native default)",
    }
    return [(FIXTURE, original, updated), (PROVENANCE, provenance_text, json.dumps(provenance, indent=2) + "\n")]


if __name__ == "__main__":
    print("*** Begin Patch")
    for name, before, after in capture():
        diff = list(difflib.unified_diff(before.splitlines(True), after.splitlines(True)))
        if not diff:
            continue
        print(f"*** Update File: {name}")
        for line in diff[2:]:
            sys.stdout.write("@@\n" if line.startswith("@@") else line)
    print("*** End Patch")
