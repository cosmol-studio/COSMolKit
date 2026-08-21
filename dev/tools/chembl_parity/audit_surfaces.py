#!/usr/bin/env python3
# pyright: reportAttributeAccessIssue=false
"""Focused differential audits for expensive claimed-parity surfaces."""

from __future__ import annotations

import argparse
import itertools
import json
import time
import warnings
from pathlib import Path
from typing import Any, cast

import cosmolkit
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDistGeom
from rdkit.Chem.Draw import rdMolDraw2D

from .audit_core import (
    Audit,
    ck_atom_state,
    ck_bond_state,
    rd_atom_state,
    rd_bond_state,
    validate_runtime_versions,
)


RDLogger.DisableLog("rdApp.*")

BOOL_PARAMS = [
    "do_isomeric_smiles",
    "do_kekule",
    "canonical",
    "clean_stereo",
    "all_bonds_explicit",
    "all_hs_explicit",
    "include_dative_bonds",
    "ignore_atom_map_numbers",
]


def iter_records(path: Path, shard: int, shards: int):
    with path.open(encoding="utf-8") as source:
        for index, line in enumerate(source):
            if index % shards == shard:
                yield json.loads(line)


def parse_pair(record: dict[str, Any], *, sanitize: bool = True):
    try:
        rd_mol = Chem.MolFromSmiles(record["smiles"], sanitize=sanitize)
    except Exception:
        rd_mol = None
    try:
        ck_mol = cosmolkit.Molecule.from_smiles(
            record["smiles"], sanitize=sanitize
        )
    except Exception:
        ck_mol = None
    return rd_mol, ck_mol


def rd_graph_state(mol: Chem.Mol) -> dict[str, Any]:
    return {
        "atoms": [rd_atom_state(atom) for atom in mol.GetAtoms()],
        "bonds": [rd_bond_state(bond) for bond in mol.GetBonds()],
    }


def ck_graph_state(mol: Any) -> dict[str, Any]:
    return {
        "atoms": [ck_atom_state(atom) for atom in mol.atoms()],
        "bonds": [ck_bond_state(bond) for bond in mol.bonds()],
    }


def rd_raw_graph_state(mol: Chem.Mol) -> dict[str, Any]:
    return {
        "atoms": [
            (
                atom.GetAtomicNum(),
                atom.GetFormalCharge(),
                str(atom.GetChiralTag()),
                atom.GetIsotope() or None,
                atom.GetIsAromatic(),
                atom.GetNumExplicitHs(),
                atom.GetDegree(),
                atom.GetNumRadicalElectrons(),
            )
            for atom in mol.GetAtoms()
        ],
        "bonds": [rd_bond_state(bond) for bond in mol.GetBonds()],
    }


def ck_raw_graph_state(mol: Any) -> dict[str, Any]:
    return {
        "atoms": [
            (
                atom.atomic_num(),
                atom.formal_charge(),
                atom.chiral_tag_name(),
                atom.isotope(),
                atom.is_aromatic(),
                atom.explicit_hydrogens(),
                atom.degree(),
                atom.num_radical_electrons(),
            )
            for atom in mol.atoms()
        ],
        "bonds": [ck_bond_state(bond) for bond in mol.bonds()],
    }


def compare_graph_state(
    audit: Audit,
    feature: str,
    record: dict[str, Any],
    rd_mol: Chem.Mol,
    ck_mol: Any,
) -> None:
    rd_value = rd_graph_state(rd_mol)
    ck_value = ck_graph_state(ck_mol)
    if rd_value == ck_value:
        audit.count(f"match.{feature}.graph_state")
    else:
        audit.mismatch(f"{feature}.graph_state", record, rd_value, ck_value)


def check_value_source_unchanged(
    audit: Audit,
    feature: str,
    record: dict[str, Any],
    molecule: Any,
    before: bytes,
) -> None:
    after = molecule.mol_to_binary()
    if after == before:
        audit.count(f"match.{feature}.source_unchanged")
    else:
        audit.mismatch(
            f"{feature}.source_mutated",
            record,
            "byte-identical source molecule",
            "source molecule changed",
        )


def fresh_raw_sanitized_molecule(record: dict[str, Any]) -> Any:
    molecule = cosmolkit.Molecule.from_smiles(record["smiles"], sanitize=False)
    molecule.sanitize_()
    return molecule


def rd_smiles_params(mol: Chem.Mol, values: tuple[bool, ...], root: str | None):
    params = Chem.SmilesWriteParams()
    params.doIsomericSmiles = values[0]
    params.doKekule = values[1]
    params.canonical = values[2]
    params.cleanStereo = values[3]
    params.allBondsExplicit = values[4]
    params.allHsExplicit = values[5]
    params.includeDativeBonds = values[6]
    params.ignoreAtomMapNumbers = values[7]
    params.rootedAtAtom = (
        -1 if root is None else 0 if root == "first" else mol.GetNumAtoms() - 1
    )
    params.doRandom = False
    return params


def branch_name(values: tuple[bool, ...], root: str | None) -> str:
    bits = "".join("1" if value else "0" for value in values)
    return f"{bits}.root_{root or 'none'}"


def audit_smiles(audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any):
    for values in itertools.product([False, True], repeat=8):
        for root in (None, "first", "last"):
            name = branch_name(values, root)
            try:
                rd_value = Chem.MolToSmiles(rd_mol, rd_smiles_params(rd_mol, values, root))
            except Exception as error:
                rd_value = f"ERROR: {type(error).__name__}: {error}"
            rooted_at_atom = (
                None if root is None else 0 if root == "first" else rd_mol.GetNumAtoms() - 1
            )
            try:
                ck_value = ck_mol.to_smiles(
                    isomeric_smiles=values[0],
                    kekule=values[1],
                    canonical=values[2],
                    clean_stereo=values[3],
                    all_bonds_explicit=values[4],
                    all_hs_explicit=values[5],
                    include_dative_bonds=values[6],
                    ignore_atom_map_numbers=values[7],
                    rooted_at_atom=rooted_at_atom,
                )
            except Exception as error:
                ck_value = f"ERROR: {type(error).__name__}: {error}"
            if rd_value == ck_value:
                audit.count("match.smiles.matrix")
            else:
                audit.mismatch("smiles.matrix", record, rd_value, ck_value, name)


def audit_bounds(audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any):
    try:
        rd_value = np.asarray(rdDistGeom.GetMoleculeBoundsMatrix(rd_mol), dtype=np.float64)
    except Exception as error:
        audit.mismatch("bounds.rdkit_error", record, repr(error), None)
        return
    try:
        ck_value = np.asarray(ck_mol.dg_bounds_matrix(), dtype=np.float64)
    except Exception as error:
        audit.mismatch("bounds.cosmolkit_error", record, "success", repr(error))
        return
    if rd_value.shape != ck_value.shape:
        audit.mismatch("bounds.shape", record, rd_value.shape, ck_value.shape)
        return
    delta = np.abs(rd_value - ck_value)
    max_error = float(delta.max(initial=0.0))
    audit.count("bounds.entries", rd_value.size)
    if max_error <= 1.0e-8:
        audit.count("match.bounds.matrix")
    else:
        flat = int(np.argmax(delta))
        pos = np.unravel_index(flat, delta.shape)
        audit.mismatch(
            "bounds.matrix",
            record,
            float(rd_value[pos]),
            float(ck_value[pos]),
            f"max_error={max_error} position={tuple(int(x) for x in pos)}",
        )


def normalize_svg(svg: str) -> str:
    return (
        svg.replace(
            "xmlns:rdkit='http://www.rdkit.org/xml'",
            "xmlns:tool='__tool_namespace__'",
        )
        .replace(
            "xmlns:cosmolkit='https://www.cosmol.org'",
            "xmlns:tool='__tool_namespace__'",
        )
        .replace("rdkit:", "tool:")
        .replace("cosmolkit:", "tool:")
    )


def rd_svg(mol: Chem.Mol, width: int, height: int) -> str:
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height, -1, -1, True)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def audit_svg(audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any):
    try:
        rd_value = normalize_svg(rd_svg(rd_mol, 300, 300))
    except Exception as error:
        audit.mismatch("svg.rdkit_error", record, repr(error), None)
        return
    try:
        ck_value = normalize_svg(ck_mol.to_svg(300, 300))
    except Exception as error:
        audit.mismatch("svg.cosmolkit_error", record, "success", repr(error))
        return
    if rd_value == ck_value:
        audit.count("match.svg.final")
    else:
        first = next(
            (i for i, pair in enumerate(zip(rd_value, ck_value)) if pair[0] != pair[1]),
            min(len(rd_value), len(ck_value)),
        )
        audit.mismatch(
            "svg.final",
            record,
            rd_value[first : first + 240],
            ck_value[first : first + 240],
            f"first_offset={first} rd_len={len(rd_value)} ck_len={len(ck_value)}",
        )


def audit_conformer(
    audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any
):
    try:
        rd_with_h = Chem.AddHs(rd_mol)
        rd_params = rdDistGeom.ETKDGv3()
        rd_params.maxIterations = 3
        rd_params.randomSeed = 61453
        rd_params.numThreads = 1
        rd_params.timeout = 0
        rd_status = int(rdDistGeom.EmbedMolecule(rd_with_h, rd_params))
        rd_coords = (
            np.asarray(rd_with_h.GetConformer(rd_status).GetPositions(), dtype=np.float64)
            if rd_status >= 0
            else None
        )
    except Exception as error:
        audit.mismatch("conformer.rdkit_error", record, repr(error), None)
        return
    try:
        ck_with_h = ck_mol.with_hydrogens()
        ck_params = cosmolkit.EmbedParameters.etkdg_v3()
        ck_params.max_iterations = 3
        ck_params.random_seed = 61453
        ck_params.num_threads = 1
        ck_params.timeout = 0
        ck_result = ck_with_h.with_3d_conformer_result(ck_params)
        ck_status = int(ck_result.conf_id())
        ck_coords = (
            np.asarray(ck_result.molecule().coordinates_3d(), dtype=np.float64)
            if ck_status >= 0
            else None
        )
    except Exception as error:
        audit.mismatch("conformer.cosmolkit_error", record, rd_status, repr(error))
        return
    if rd_status != ck_status:
        audit.mismatch("conformer.status", record, rd_status, ck_status)
        return
    if rd_coords is None or ck_coords is None:
        audit.count("match.conformer.failed_status")
        return
    if rd_coords.shape != ck_coords.shape:
        audit.mismatch("conformer.shape", record, rd_coords.shape, ck_coords.shape)
        return
    max_error = float(np.abs(rd_coords - ck_coords).max(initial=0.0))
    if max_error <= 1.0e-6:
        audit.count("match.conformer.coordinates")
    else:
        audit.mismatch(
            "conformer.coordinates",
            record,
            rd_coords.tolist(),
            ck_coords.tolist(),
            f"max_error={max_error}",
        )


def molblock_body(block: str) -> str:
    return "\n".join(line for line in block.splitlines()[3:] if line != "$$$$")


def stereo_direction_states_equivalent(
    expected: list[tuple[Any, ...]], actual: list[tuple[Any, ...]]
) -> bool:
    """Compare slash directions modulo one inversion per stereo component."""
    if len(expected) != len(actual):
        return False

    directional = {"ENDUPRIGHT", "ENDDOWNRIGHT"}
    for expected_bond, actual_bond in zip(expected, actual):
        if expected_bond[:3] != actual_bond[:3] or expected_bond[4:] != actual_bond[4:]:
            return False
        expected_direction = expected_bond[3]
        actual_direction = actual_bond[3]
        if expected_direction == actual_direction:
            continue
        if {expected_direction, actual_direction} != directional:
            return False

    direction_bonds = {
        index
        for index, bond in enumerate(expected)
        if bond[3] in directional
    }
    parent = {index: index for index in direction_bonds}

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    constrained: set[int] = set()
    for bond in expected:
        if bond[2] != "DOUBLE" or bond[4] not in {"E", "Z", "CIS", "TRANS"}:
            continue
        endpoints = {bond[0], bond[1]}
        adjacent = [
            index
            for index in direction_bonds
            if expected[index][0] in endpoints or expected[index][1] in endpoints
        ]
        constrained.update(adjacent)
        for other in adjacent[1:]:
            union(adjacent[0], other)

    component_inversion: dict[int, bool] = {}
    for index in direction_bonds:
        inverted = expected[index][3] != actual[index][3]
        if inverted and index not in constrained:
            return False
        root = find(index)
        previous = component_inversion.setdefault(root, inverted)
        if previous != inverted:
            return False
    return True


def compare_parsed_source(
    audit: Audit,
    record: dict[str, Any],
    dimension: str,
    rd_source: Chem.Mol,
    ck_source: Any,
) -> None:
    rd_smiles = Chem.MolToSmiles(rd_source)
    ck_smiles = ck_source.to_smiles()
    if rd_smiles == ck_smiles:
        audit.count(f"match.io.{dimension}.parsed_smiles")
    else:
        audit.mismatch(f"io.{dimension}.parsed_smiles", record, rd_smiles, ck_smiles)
    rd_atoms = [rd_atom_state(atom) for atom in rd_source.GetAtoms()]
    ck_atoms = [ck_atom_state(atom) for atom in ck_source.atoms()]
    if rd_atoms == ck_atoms:
        audit.count(f"match.io.{dimension}.atom_state")
    else:
        audit.mismatch(f"io.{dimension}.atom_state", record, rd_atoms, ck_atoms)
    rd_bonds = [rd_bond_state(bond) for bond in rd_source.GetBonds()]
    ck_bonds = [ck_bond_state(bond) for bond in ck_source.bonds()]
    if stereo_direction_states_equivalent(rd_bonds, ck_bonds):
        audit.count(f"match.io.{dimension}.bond_state")
    else:
        audit.mismatch(f"io.{dimension}.bond_state", record, rd_bonds, ck_bonds)
    rd_conf = rd_source.GetConformer()
    rd_coords = np.asarray(rd_conf.GetPositions(), dtype=np.float64)
    ck_coords = np.asarray(
        ck_source.coordinates_2d()
        if dimension == "2d"
        else ck_source.coordinates_3d(),
        dtype=np.float64,
    )
    if dimension == "2d":
        rd_coords = rd_coords[:, :2]
        ck_coords = ck_coords[:, :2]
    max_error = float(np.abs(rd_coords - ck_coords).max(initial=0.0))
    if max_error <= 1.0e-8:
        audit.count(f"match.io.{dimension}.coordinates")
    else:
        audit.mismatch(
            f"io.{dimension}.coordinates",
            record,
            rd_coords.tolist(),
            ck_coords.tolist(),
            f"max_error={max_error}",
        )


def audit_io(audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any):
    del ck_mol
    sources: list[tuple[str, Chem.Mol]] = []
    rd_2d = Chem.Mol(rd_mol)
    try:
        AllChem.Compute2DCoords(rd_2d)
        Chem.WedgeMolBonds(rd_2d, rd_2d.GetConformer())
        sources.append(("2d", rd_2d))
    except Exception as error:
        audit.mismatch("io.2d.source_error", record, repr(error), None)

    rd_3d = Chem.Mol(rd_mol)
    try:
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 1
        if rdDistGeom.EmbedMolecule(rd_3d, params) == 0:
            sources.append(("3d", rd_3d))
        else:
            audit.count("skip.io.3d_embed")
    except Exception:
        audit.count("skip.io.3d_embed")

    for dimension, rd_source_before in sources:
        source_block = Chem.MolToMolBlock(
            rd_source_before, includeStereo=True, kekulize=True, forceV3000=True
        )
        rd_source = Chem.MolFromMolBlock(source_block, sanitize=True, removeHs=False)
        if rd_source is None:
            audit.mismatch(f"io.{dimension}.rdkit_reparse", record, "success", None)
            continue
        try:
            ck_source = cosmolkit.Molecule.read_mol_from_str(
                source_block,
                sanitize=True,
                coordinate_dim=dimension,
                remove_hs=False,
            )
        except Exception as error:
            audit.mismatch(
                f"io.{dimension}.cosmolkit_reparse", record, "success", repr(error)
            )
            continue
        if ck_source is None:
            audit.mismatch(f"io.{dimension}.cosmolkit_reparse", record, "success", None)
            continue
        compare_parsed_source(audit, record, dimension, rd_source, ck_source)

        for include_stereo, kekulize, force_v3000 in itertools.product(
            [False, True], repeat=3
        ):
            name = (
                f"{dimension}.stereo{int(include_stereo)}.kek{int(kekulize)}."
                f"v3k{int(force_v3000)}"
            )
            try:
                rd_block = Chem.MolToMolBlock(
                    rd_source,
                    includeStereo=include_stereo,
                    confId=-1,
                    kekulize=kekulize,
                    forceV3000=force_v3000,
                )
                rd_value = molblock_body(rd_block)
            except Exception as error:
                rd_value = f"ERROR: {type(error).__name__}: {error}"
            try:
                ck_block = (
                    ck_source.to_2d_sdf_string(
                        format="v3000" if force_v3000 else "v2000",
                        include_stereo=include_stereo,
                        kekulize=kekulize,
                    )
                    if dimension == "2d"
                    else ck_source.to_3d_sdf_string(
                        format="v3000" if force_v3000 else "v2000",
                        include_stereo=include_stereo,
                        kekulize=kekulize,
                    )
                )
                ck_value = molblock_body(ck_block)
            except Exception as error:
                ck_value = f"ERROR: {type(error).__name__}: {error}"
            if rd_value == ck_value:
                audit.count("match.io.writer_branch")
            else:
                audit.mismatch("io.writer_branch", record, rd_value, ck_value, name)


def audit_explicit_h(
    audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any
):
    try:
        rd_value = Chem.AddHs(rd_mol)
    except Exception as error:
        audit.mismatch("explicit_h.rdkit_error", record, repr(error), None)
        return
    try:
        ck_value = ck_mol.with_hydrogens()
    except Exception as error:
        audit.mismatch("explicit_h.cosmolkit_error", record, "success", repr(error))
        return
    pairs = {
        "num_atoms": (rd_value.GetNumAtoms(), ck_value.num_atoms()),
        "num_bonds": (rd_value.GetNumBonds(), ck_value.num_bonds()),
        "atom_state": (
            [rd_atom_state(atom) for atom in rd_value.GetAtoms()],
            [ck_atom_state(atom) for atom in ck_value.atoms()],
        ),
        "bond_state": (
            [rd_bond_state(bond) for bond in rd_value.GetBonds()],
            [ck_bond_state(bond) for bond in ck_value.bonds()],
        ),
    }
    for name, (rd_item, ck_item) in pairs.items():
        if rd_item == ck_item:
            audit.count(f"match.explicit_h.{name}")
        else:
            audit.mismatch(f"explicit_h.{name}", record, rd_item, ck_item)


def audit_topology_operations(
    audit: Audit, record: dict[str, Any], rd_raw: Chem.Mol, ck_raw: Any
) -> None:
    rd_raw_value = rd_raw_graph_state(rd_raw)
    ck_raw_value = ck_raw_graph_state(ck_raw)
    if rd_raw_value == ck_raw_value:
        audit.count("match.topology.raw_parse.graph_state")
    else:
        audit.mismatch(
            "topology.raw_parse.graph_state", record, rd_raw_value, ck_raw_value
        )

    rd_sanitized = Chem.Mol(rd_raw)
    try:
        Chem.SanitizeMol(rd_sanitized)
        rd_sanitize_error = None
    except Exception as error:
        rd_sanitize_error = f"{type(error).__name__}: {error}"
    ck_raw_before = ck_raw.mol_to_binary()
    try:
        ck_sanitized = ck_raw.sanitize()
        ck_sanitize_error = None
    except Exception as error:
        ck_sanitized = None
        ck_sanitize_error = f"{type(error).__name__}: {error}"
    check_value_source_unchanged(
        audit, "topology.sanitize.value", record, ck_raw, ck_raw_before
    )
    try:
        ck_sanitized_in_place = cosmolkit.Molecule.from_smiles(
            record["smiles"], sanitize=False
        )
        ck_sanitized_in_place.sanitize_()
        ck_in_place_error = None
    except Exception as error:
        ck_in_place_error = f"{type(error).__name__}: {error}"
    if (rd_sanitize_error is None) != (ck_sanitize_error is None):
        audit.mismatch(
            "topology.sanitize.status",
            record,
            rd_sanitize_error or "success",
            ck_sanitize_error or "success",
        )
        return
    if (rd_sanitize_error is None) != (ck_in_place_error is None):
        audit.mismatch(
            "topology.sanitize.in_place_status",
            record,
            rd_sanitize_error or "success",
            ck_in_place_error or "success",
        )
        return
    if rd_sanitize_error is not None:
        audit.count("match.topology.sanitize.rejected")
        audit.count("match.topology.sanitize.in_place_rejected")
        return
    assert ck_sanitized is not None
    compare_graph_state(
        audit,
        "topology.sanitize.value",
        record,
        rd_sanitized,
        ck_sanitized,
    )
    compare_graph_state(
        audit,
        "topology.sanitize.in_place",
        record,
        rd_sanitized,
        ck_sanitized_in_place,
    )

    rd_mol = rd_sanitized
    ck_mol = ck_sanitized

    for clear_aromatic_flags in (False, True):
        branch = f"clear_aromatic_flags_{str(clear_aromatic_flags).lower()}"
        rd_kekulized = Chem.Mol(rd_mol)
        try:
            Chem.Kekulize(
                rd_kekulized, clearAromaticFlags=clear_aromatic_flags
            )
            rd_error = None
        except Exception as error:
            rd_error = f"{type(error).__name__}: {error}"
        ck_source_before = ck_mol.mol_to_binary()
        try:
            ck_kekulized = ck_mol.with_kekulized_bonds(clear_aromatic_flags)
            ck_error = None
        except Exception as error:
            ck_kekulized = None
            ck_error = f"{type(error).__name__}: {error}"
        check_value_source_unchanged(
            audit, f"topology.kekulize.{branch}.value", record, ck_mol, ck_source_before
        )
        try:
            ck_kekulized_in_place = fresh_raw_sanitized_molecule(record)
            ck_kekulized_in_place.kekulize_(clear_aromatic_flags)
            ck_in_place_error = None
        except Exception as error:
            ck_kekulized_in_place = None
            ck_in_place_error = f"{type(error).__name__}: {error}"
        if (rd_error is None) != (ck_error is None):
            audit.mismatch(
                f"topology.kekulize.{branch}.status",
                record,
                rd_error or "success",
                ck_error or "success",
            )
            continue
        if (rd_error is None) != (ck_in_place_error is None):
            audit.mismatch(
                f"topology.kekulize.{branch}.in_place_status",
                record,
                rd_error or "success",
                ck_in_place_error or "success",
            )
            continue
        if rd_error is not None:
            audit.count(f"match.topology.kekulize.{branch}.rejected")
            audit.count(f"match.topology.kekulize.{branch}.in_place_rejected")
            continue
        assert ck_kekulized is not None and ck_kekulized_in_place is not None
        compare_graph_state(
            audit,
            f"topology.kekulize.{branch}.value",
            record,
            rd_kekulized,
            ck_kekulized,
        )
        compare_graph_state(
            audit,
            f"topology.kekulize.{branch}.in_place",
            record,
            rd_kekulized,
            ck_kekulized_in_place,
        )

    try:
        rd_with_h = Chem.AddHs(rd_mol)
        ck_with_h = ck_mol.with_hydrogens()
    except Exception as error:
        audit.mismatch(
            "topology.remove_h.initialization", record, "success", repr(error)
        )
        return
    for sanitize in (False, True):
        branch = f"sanitize_{str(sanitize).lower()}"
        try:
            rd_without_h = Chem.RemoveHs(Chem.Mol(rd_with_h), sanitize=sanitize)
            rd_error = None
        except Exception as error:
            rd_without_h = None
            rd_error = f"{type(error).__name__}: {error}"
        ck_source_before = ck_with_h.mol_to_binary()
        try:
            ck_without_h = ck_with_h.without_hydrogens(sanitize=sanitize)
            ck_error = None
        except Exception as error:
            ck_without_h = None
            ck_error = f"{type(error).__name__}: {error}"
        check_value_source_unchanged(
            audit,
            f"topology.remove_h.{branch}.value",
            record,
            ck_with_h,
            ck_source_before,
        )
        try:
            ck_without_h_in_place = fresh_raw_sanitized_molecule(
                record
            ).with_hydrogens()
            ck_without_h_in_place.remove_hydrogens_(sanitize=sanitize)
            ck_in_place_error = None
        except Exception as error:
            ck_without_h_in_place = None
            ck_in_place_error = f"{type(error).__name__}: {error}"
        if (rd_error is None) != (ck_error is None):
            audit.mismatch(
                f"topology.remove_h.{branch}.status",
                record,
                rd_error or "success",
                ck_error or "success",
            )
            continue
        if (rd_error is None) != (ck_in_place_error is None):
            audit.mismatch(
                f"topology.remove_h.{branch}.in_place_status",
                record,
                rd_error or "success",
                ck_in_place_error or "success",
            )
            continue
        if rd_error is not None:
            audit.count(f"match.topology.remove_h.{branch}.rejected")
            audit.count(f"match.topology.remove_h.{branch}.in_place_rejected")
            continue
        assert (
            rd_without_h is not None
            and ck_without_h is not None
            and ck_without_h_in_place is not None
        )
        compare_graph_state(
            audit,
            f"topology.remove_h.{branch}.value",
            record,
            rd_without_h,
            ck_without_h,
        )
        compare_graph_state(
            audit,
            f"topology.remove_h.{branch}.in_place",
            record,
            rd_without_h,
            ck_without_h_in_place,
        )


def audit_fragments(
    audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any
) -> None:
    try:
        rd_fragments = list(
            Chem.GetMolFrags(rd_mol, asMols=True, sanitizeFrags=True)
        )
    except Exception as error:
        audit.mismatch("fragments.rdkit_error", record, repr(error), None)
        return
    ck_source_before = ck_mol.mol_to_binary()
    try:
        ck_fragments = ck_mol.fragments()
    except Exception as error:
        audit.mismatch("fragments.cosmolkit_error", record, "success", repr(error))
        return
    check_value_source_unchanged(
        audit, "fragments.value", record, ck_mol, ck_source_before
    )
    if len(rd_fragments) != len(ck_fragments):
        audit.mismatch(
            "fragments.count", record, len(rd_fragments), len(ck_fragments)
        )
        return
    for index, (rd_fragment, ck_fragment) in enumerate(
        zip(rd_fragments, ck_fragments)
    ):
        rd_value = rd_graph_state(rd_fragment)
        ck_value = ck_graph_state(ck_fragment)
        if rd_value == ck_value:
            audit.count("match.fragments.graph_state")
        else:
            audit.mismatch(
                "fragments.graph_state",
                record,
                rd_value,
                ck_value,
                f"fragment_index={index}",
            )


def compare_2d_coordinates(
    audit: Audit,
    feature: str,
    record: dict[str, Any],
    rd_mol: Chem.Mol,
    ck_mol: Any,
) -> None:
    compare_graph_state(audit, feature, record, rd_mol, ck_mol)
    rd_coords = np.asarray(
        rd_mol.GetConformer().GetPositions(), dtype=np.float64
    )[:, :2]
    ck_coords = np.asarray(ck_mol.coordinates_2d(), dtype=np.float64)[:, :2]
    if rd_coords.shape != ck_coords.shape:
        audit.mismatch(
            f"{feature}.coordinate_shape", record, rd_coords.shape, ck_coords.shape
        )
        return
    max_error = float(np.abs(rd_coords - ck_coords).max(initial=0.0))
    if max_error <= 1.0e-8:
        audit.count(f"match.{feature}.coordinates")
    else:
        audit.mismatch(
            f"{feature}.coordinates",
            record,
            rd_coords.tolist(),
            ck_coords.tolist(),
            f"max_error={max_error}",
        )


def audit_coordinates_2d(
    audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any
) -> None:
    rd_with_2d = Chem.Mol(rd_mol)
    try:
        AllChem.Compute2DCoords(rd_with_2d)
    except Exception as error:
        audit.mismatch("coordinates_2d.rdkit_error", record, repr(error), None)
        return
    ck_source_before = ck_mol.mol_to_binary()
    try:
        ck_with_2d = ck_mol.with_2d_coordinates()
    except Exception as error:
        audit.mismatch("coordinates_2d.cosmolkit_error", record, "success", repr(error))
        return
    check_value_source_unchanged(
        audit, "coordinates_2d.value", record, ck_mol, ck_source_before
    )
    compare_2d_coordinates(
        audit, "coordinates_2d.value", record, rd_with_2d, ck_with_2d
    )
    try:
        ck_with_2d_in_place = cosmolkit.Molecule.from_smiles(record["smiles"])
        ck_with_2d_in_place.compute_2d_coordinates_()
    except Exception as error:
        audit.mismatch(
            "coordinates_2d.in_place_status",
            record,
            "success",
            f"{type(error).__name__}: {error}",
        )
        return
    compare_2d_coordinates(
        audit,
        "coordinates_2d.in_place",
        record,
        rd_with_2d,
        ck_with_2d_in_place,
    )


def capture_observation(function: Any) -> dict[str, Any]:
    try:
        return {"status": "ok", "value": function()}
    except Exception as error:
        return {
            "status": "error",
            "type": type(error).__name__,
            "message": str(error),
        }


def binary_observations(mol: Any) -> dict[str, Any]:
    return {
        "graph": capture_observation(lambda: ck_graph_state(mol)),
        "smiles": capture_observation(mol.to_smiles),
        "kekule_smiles": capture_observation(lambda: mol.to_smiles(kekule=True)),
        "hash": capture_observation(mol.hash),
        "morgan": capture_observation(lambda: mol.fingerprint_morgan().on_bits()),
        "mol_wt": capture_observation(lambda: float(cosmolkit.calc_mol_wt(mol)).hex()),
        "tpsa": capture_observation(lambda: float(cosmolkit.calc_tpsa(mol)).hex()),
        "has_2d": capture_observation(mol.has_2d_coordinates),
        "coordinates_2d": capture_observation(
            lambda: np.asarray(mol.coordinates_2d(), dtype=np.float64).tolist()
            if mol.has_2d_coordinates()
            else None
        ),
        "num_conformers": capture_observation(mol.num_conformers),
    }


def audit_binary_roundtrip(
    audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any
) -> None:
    del rd_mol
    try:
        source = ck_mol.with_2d_coordinates()
        before = binary_observations(source)
        payload = source.mol_to_binary()
        restored_method = cosmolkit.Molecule.mol_from_binary(payload)
        restored_function = cosmolkit.mol_from_binary(payload)
    except Exception as error:
        audit.mismatch("binary_roundtrip.error", record, "success", repr(error))
        return
    for branch, restored in (
        ("class_method", restored_method),
        ("module_function", restored_function),
    ):
        after = binary_observations(restored)
        for observation, expected in before.items():
            actual = after[observation]
            if actual == expected:
                audit.count(
                    f"match.binary_roundtrip.{branch}.{observation}"
                )
            else:
                audit.mismatch(
                    f"binary_roundtrip.{branch}.{observation}",
                    record,
                    expected,
                    actual,
                )
        restored_payload = restored.mol_to_binary()
        if restored_payload == payload:
            audit.count(f"match.binary_roundtrip.{branch}.deterministic_bytes")
        else:
            audit.mismatch(
                f"binary_roundtrip.{branch}.bytes",
                record,
                "byte-identical reserialization",
                "binary payload changed",
            )


def audit_forcefield(
    audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any
):
    rd_uff = bool(AllChem.UFFHasAllMoleculeParams(rd_mol))
    rd_mmff = bool(AllChem.MMFFHasAllMoleculeParams(rd_mol))
    try:
        ck_uff = bool(ck_mol.has_uff_params())
    except Exception as error:
        ck_uff = f"ERROR: {type(error).__name__}: {error}"
    try:
        ck_mmff = bool(ck_mol.has_mmff_params())
    except Exception as error:
        ck_mmff = f"ERROR: {type(error).__name__}: {error}"
    for name, rd_value, ck_value in (
        ("uff_params", rd_uff, ck_uff),
        ("mmff_params", rd_mmff, ck_mmff),
    ):
        if rd_value == ck_value:
            audit.count(f"match.forcefield.{name}")
        else:
            audit.mismatch(f"forcefield.{name}", record, rd_value, ck_value)

    try:
        rd_with_h = Chem.AddHs(rd_mol)
        ck_with_h = ck_mol.with_hydrogens()
        explicit_h_availability = (
            (
                "uff_params_explicit_h",
                bool(AllChem.UFFHasAllMoleculeParams(rd_with_h)),
                bool(ck_with_h.has_uff_params()),
            ),
            (
                "mmff_params_explicit_h",
                bool(AllChem.MMFFHasAllMoleculeParams(rd_with_h)),
                bool(ck_with_h.has_mmff_params()),
            ),
        )
        for name, rd_value, ck_value in explicit_h_availability:
            if rd_value == ck_value:
                audit.count(f"match.forcefield.{name}")
            else:
                audit.mismatch(f"forcefield.{name}", record, rd_value, ck_value)

        rd_initial = Chem.Mol(rd_mol)
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 61453
        params.numThreads = 1
        if rdDistGeom.EmbedMolecule(rd_initial, params) != 0:
            audit.count("skip.forcefield.embed")
            return
        initial_coords = np.asarray(rd_initial.GetConformer().GetPositions(), dtype=np.float64)
        ck_initial = ck_mol.with_only_3d_conformer(initial_coords)
    except Exception as error:
        audit.mismatch("forcefield.initialization", record, "success", repr(error))
        return

    if bool(AllChem.UFFHasAllMoleculeParams(rd_initial)) and ck_initial.has_uff_params():
        rd_opt = Chem.Mol(rd_initial)
        try:
            # UFFOptimizeMolecule() exposes only the status through Python, while
            # COSMolKit's source-level result also carries the energy returned by
            # ForceFieldsHelper::OptimizeMolecule(). The one-conformer batch API
            # exposes that exact source pair; rebuilding a force field after the
            # optimization can legitimately produce a different energy because
            # construction selects coordinate-dependent nonbonded terms.
            rd_status, rd_energy = AllChem.UFFOptimizeMoleculeConfs(
                rd_opt,
                numThreads=1,
                maxIters=200,
                vdwThresh=10.0,
                ignoreInterfragInteractions=True,
            )[0]
            rd_status = int(rd_status)
            rd_energy = float(rd_energy)
            rd_coords = np.asarray(rd_opt.GetConformer().GetPositions(), dtype=np.float64)
            ck_result = ck_initial.with_uff_optimized(max_iters=200, vdw_thresh=10.0)
            ck_status = int(ck_result.status_code())
            ck_energy = float(ck_result.energy())
            ck_coords = np.asarray(ck_result.molecule().coordinates_3d(), dtype=np.float64)
            max_coord_error = float(np.abs(rd_coords - ck_coords).max(initial=0.0))
            if (
                rd_status == ck_status
                and abs(rd_energy - ck_energy) <= 1.0e-6
                and max_coord_error <= 1.0e-6
            ):
                audit.count("match.forcefield.uff_optimized")
            else:
                audit.mismatch(
                    "forcefield.uff_optimized",
                    record,
                    {"status": rd_status, "energy": rd_energy},
                    {"status": ck_status, "energy": ck_energy},
                    f"max_coord_error={max_coord_error}",
                )
        except Exception as error:
            audit.mismatch("forcefield.uff_error", record, "success", repr(error))

    if bool(AllChem.MMFFHasAllMoleculeParams(rd_initial)) and ck_initial.has_mmff_params():
        rd_opt = Chem.Mol(rd_initial)
        rd_result: tuple[int, np.ndarray[Any, np.dtype[np.float64]]] | str
        try:
            rd_status = int(
                AllChem.MMFFOptimizeMolecule(
                    rd_opt,
                    mmffVariant="MMFF94",
                    maxIters=200,
                    nonBondedThresh=100.0,
                    confId=0,
                    ignoreInterfragInteractions=True,
                )
            )
            rd_coords = np.asarray(
                rd_opt.GetConformer().GetPositions(), dtype=np.float64
            )
            rd_result = (rd_status, rd_coords)
        except Exception as error:
            rd_result = f"{type(error).__name__}: {error}"

        ck_result_value: tuple[int, np.ndarray[Any, np.dtype[np.float64]]] | str
        try:
            ck_result = ck_initial.with_mmff_optimized(max_iters=200)
            ck_status = int(ck_result.status_code())
            ck_coords = np.asarray(
                ck_result.molecule().coordinates_3d(), dtype=np.float64
            )
            ck_result_value = (ck_status, ck_coords)
        except Exception as error:
            ck_result_value = f"{type(error).__name__}: {error}"

        if isinstance(rd_result, str) or isinstance(ck_result_value, str):
            if isinstance(rd_result, str) and isinstance(ck_result_value, str):
                audit.count("match.forcefield.mmff_rejected")
            else:
                audit.mismatch(
                    "forcefield.mmff_error",
                    record,
                    rd_result if isinstance(rd_result, str) else "success",
                    ck_result_value
                    if isinstance(ck_result_value, str)
                    else "success",
                )
        else:
            rd_status, rd_coords = rd_result
            ck_status, ck_coords = ck_result_value
            max_coord_error = float(np.abs(rd_coords - ck_coords).max(initial=0.0))
            if rd_status == ck_status and max_coord_error <= 1.0e-6:
                audit.count("match.forcefield.mmff_optimized")
            else:
                audit.mismatch(
                    "forcefield.mmff_optimized",
                    record,
                    {"status": rd_status},
                    {"status": ck_status},
                    f"max_coord_error={max_coord_error}",
                )


INCHI_OPTIONS = (
    "",
    "-AuxNone",
    "-FixedH",
    "-RecMet",
    "-SNon",
    "-SRel",
    "-SRac",
    "-SUU",
    "-SLUUD",
    "-FixedH -SUU -SLUUD",
)


def audit_inchi_options(
    audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any
):
    for options in INCHI_OPTIONS:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rd_value = Chem.MolToInchi(rd_mol, options=options)
        except Exception as error:
            rd_value = f"ERROR: {type(error).__name__}: {error}"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ck_value = ck_mol.to_inchi(options)
        except Exception as error:
            ck_value = f"ERROR: {type(error).__name__}: {error}"
        if rd_value == ck_value:
            audit.count("match.inchi.options.mol_to_inchi")
        else:
            audit.mismatch(
                "inchi.options.mol_to_inchi", record, rd_value, ck_value, options
            )

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rd_key = Chem.MolToInchiKey(rd_mol, options=options)
        except Exception as error:
            rd_key = f"ERROR: {type(error).__name__}: {error}"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ck_key = ck_mol.to_inchi_key(options)
        except Exception as error:
            ck_key = f"ERROR: {type(error).__name__}: {error}"
        if rd_key == ck_key:
            audit.count("match.inchi.options.mol_to_key")
        else:
            audit.mismatch("inchi.options.mol_to_key", record, rd_key, ck_key, options)


def inchi_parse_state(mol: Any, rdkit: bool) -> dict[str, Any] | None:
    if mol is None:
        return None
    if rdkit:
        return {
            "atoms": [rd_atom_state(atom) for atom in mol.GetAtoms()],
            "bonds": [rd_bond_state(bond) for bond in mol.GetBonds()],
        }
    return {
        "atoms": [ck_atom_state(atom) for atom in mol.atoms()],
        "bonds": [ck_bond_state(bond) for bond in mol.bonds()],
    }


def audit_inchi_parse(
    audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any
):
    del ck_mol
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            inchi = cast(str, cast(object, Chem.MolToInchi(rd_mol)))
    except Exception as error:
        audit.mismatch("inchi.parse.source_error", record, repr(error), None)
        return
    for sanitize, remove_hs in itertools.product([False, True], repeat=2):
        name = f"sanitize{int(sanitize)}.remove_hs{int(remove_hs)}"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rd_value = inchi_parse_state(
                    Chem.MolFromInchi(
                        inchi, sanitize=sanitize, removeHs=remove_hs
                    ),
                    True,
                )
        except Exception as error:
            rd_value = f"ERROR: {type(error).__name__}: {error}"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ck_value = inchi_parse_state(
                    cosmolkit.Molecule.from_inchi(
                        inchi, sanitize=sanitize, remove_hs=remove_hs
                    ),
                    False,
                )
        except Exception as error:
            ck_value = f"ERROR: {type(error).__name__}: {error}"
        if rd_value == ck_value:
            audit.count("match.inchi.parse_branch")
        else:
            audit.mismatch("inchi.parse_branch", record, rd_value, ck_value, name)


MODE_FUNCTIONS = {
    "binary_roundtrip": audit_binary_roundtrip,
    "smiles": audit_smiles,
    "bounds": audit_bounds,
    "coordinates_2d": audit_coordinates_2d,
    "svg": audit_svg,
    "conformer": audit_conformer,
    "fragments": audit_fragments,
    "io": audit_io,
    "explicit_h": audit_explicit_h,
    "forcefield": audit_forcefield,
    "inchi_options": audit_inchi_options,
    "inchi_parse": audit_inchi_parse,
    "topology_operations": audit_topology_operations,
}


def main() -> None:
    validate_runtime_versions()

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=sorted(MODE_FUNCTIONS), required=True)
    parser.add_argument("--limit", type=int, required=True)
    parser.add_argument("--max-atoms", type=int, default=80)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--max-examples", type=int, default=300)
    args = parser.parse_args()

    audit = Audit(args.output, args.max_examples)
    selected = 0
    scanned = 0
    started = time.monotonic()
    try:
        for record in iter_records(args.input, args.shard, args.shards):
            if selected >= args.limit:
                break
            scanned += 1
            rd_mol, ck_mol = parse_pair(
                record, sanitize=args.mode != "topology_operations"
            )
            if rd_mol is None and ck_mol is None:
                audit.count("skip.both_parse_rejected")
                continue
            if rd_mol is None:
                audit.mismatch("parse", record, "rejected", "COSMolKit accepted")
                continue
            if ck_mol is None:
                audit.mismatch("parse", record, "success", "COSMolKit parse failure")
                continue
            if rd_mol.GetNumAtoms() > args.max_atoms:
                audit.count("skip.max_atoms")
                continue
            selected += 1
            MODE_FUNCTIONS[args.mode](audit, record, rd_mol, ck_mol)
            if selected % 100 == 0:
                print(
                    json.dumps(
                        {
                            "mode": args.mode,
                            "selected": selected,
                            "scanned": scanned,
                            "elapsed": time.monotonic() - started,
                        }
                    ),
                    flush=True,
                )
    finally:
        audit.count("records.scanned", scanned)
        audit.finish(selected, args.mode, args.shard, args.shards)


if __name__ == "__main__":
    main()
