"""Exact coordinate-record PDB projection through the pinned Gemmi API."""
import json
import gemmi

GEMMI = "0.7.5"


def parse(text):
    if gemmi.__version__ != GEMMI:
        raise RuntimeError(f"Gemmi {GEMMI} required, found {gemmi.__version__}")
    row = json.loads(text)
    formats = {"pdb": gemmi.CoorFormat.Pdb, "cif": gemmi.CoorFormat.Mmcif}
    return gemmi.read_structure_string(row["text"], merge_chain_parts=False,
                                       format=formats[row["format"]])


def pdb(structure, profile):
    if profile["gemmi_version"] != GEMMI:
        raise RuntimeError("BIO profile uses a different Gemmi version")
    options = gemmi.PdbWriteOptions(minimal=True)
    options.cryst1_record = False
    for name in ("ter_records", "numbered_ter", "ter_ignores_type", "preserve_serial", "end_record"):
        setattr(options, name, profile[name])
    # CK's public coordinate writer emits precisely these record categories.
    records = {"ATOM", "HETATM", "ANISOU", "TER", "MODEL", "ENDMDL", "END"}
    return "".join(line for line in structure.make_pdb_string(options).splitlines(keepends=True)
                   if line[:6].strip() in records).encode("utf-8")
