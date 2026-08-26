"""Explicit top-level schemas for generated RDKit JSONL expected data."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class OutputSchema:
    required: dict[str, frozenset[str]]
    optional: dict[str, frozenset[str]]


def _fields(specification: str) -> dict[str, frozenset[str]]:
    if not specification:
        return {}
    fields = {}
    for entry in specification.split(","):
        name, type_spec = entry.split(":", 1)
        fields[name] = frozenset(type_spec.split("|"))
    return fields


def _schema(required: str, optional: str = "") -> OutputSchema:
    return OutputSchema(required=_fields(required), optional=_fields(optional))


SCHEMAS = {
    "ciplabeler.jsonl": _schema(
        "calls:array,case_id:string,initial_state:null|object,parse_error:null|string,parse_status:string,schema_version:integer,source:object"
    ),
    "graph_features.jsonl": _schema(
        "addhs_removehs:null|object,chiral_centers:null|object,direct:null|object,error:null|string,possible_stereo:null|object,rdkit_ok:boolean,smiles:string,with_hs:null|object"
    ),
    "sdf_write.jsonl": _schema(
        "branches:object,error:null|string,rdkit_ok:boolean,smiles:string",
        "source_2d_error:null|string,source_2d_molblock:null|string,source_3d_error:null|string,source_3d_molblock:null|string",
    ),
    "sdf_read.jsonl": _schema(
        "case_id:string,dimension:string,error:null|string,format:string,rdkit_ok:boolean,sdf:null|string,smiles:null|string,stereo_markers:string",
        "api:string,atom_properties:array,atoms:array,bond_properties:array,bonds:array,chiral_tags:array,operation:string,positions:array,process_property_lists:boolean,properties:object,remove_hs:boolean,sanitize:boolean,smiles_out:object,strict_parsing:boolean",
    ),
    "kekulize_clear_flags_false.jsonl": _schema(
        "atom_is_aromatic:null|array,bond_is_aromatic:null|array,bond_types:null|array,kekulize_error:null|string,kekulize_ok:boolean,parse_error:null|string,parse_ok:boolean,smiles:string"
    ),
    "delete_substructs_onlyfrags_chirality.jsonl": _schema(
        "case:string,error:null,num_atoms:integer,num_bonds:integer,only_frags:boolean,rdkit_ok:boolean,result_smiles:string,smarts:string,smiles:string,use_chirality:boolean"
    ),
    "smarts.jsonl": _schema(
        "atom_mappings:null|array,case_id:string,kind:string,num_atoms:null|integer,num_bonds:null|integer,parse_ok:boolean,smarts:string,source:string,target_smiles:null|string,written_smarts:null|string"
    ),
    "smiles_writer.jsonl": _schema(
        "branches:object,error:null|string,rdkit_ok:boolean,smiles:string"
    ),
    "isomeric_smiles.jsonl": _schema(
        "error:null|string,isomeric_smiles:null|string,rdkit_ok:boolean,smiles:string"
    ),
    "dg_bounds_matrix.jsonl": _schema(
        "bounds:null|array,error:null|string,rdkit_ok:boolean,smiles:string"
    ),
    "conformer_generation.jsonl": _schema(
        "attrs:object,case_id:string,conformers:array,error:null,failure_counts:array,ids:null|array,mode:string,preset:string,rdkit_ok:boolean,source:string,source_kind:string,status:null|integer"
    ),
    "conformer_generation_library.jsonl": _schema(
        "coords:null|array,error:null|string,error_stage:null|string,max_iterations:integer,preset:string,rdkit_add_hs_ok:boolean,rdkit_embed_ok:boolean,rdkit_parse_ok:boolean,seed:integer,smiles:string,status:null|integer,timeout:integer"
    ),
    "confseq_embed_template.jsonl": _schema(
        "case_id:string,coords:array,error:null,in_smiles:string,optimize_with_uff:boolean,pre_add_hs_summary:object,preset:string,random_seed:integer,rdkit_ok:boolean,sdf:string,status:integer,td_smiles:string,with_h_summary:object"
    ),
    "forcefield_params.jsonl": _schema(
        "embedded:object,error:null|string,mmff:object,mmff_explicit_h:object,rdkit_ok:boolean,smiles:string,uff:object,uff_explicit_h:object"
    ),
    "forcefield_coverage.jsonl": _schema(
        "error:null|string,mmff:object,mmff_explicit_h:object,rdkit_ok:boolean,smiles:string,uff:object,uff_explicit_h:object"
    ),
    "mmff_builtin.jsonl": _schema(
        "atom_types:array,error:null,fixture:string,has_all:boolean,line_number:integer,num_atoms:integer,props_ok:boolean,rdkit_ok:boolean,row_name:string,smiles:string,variant:string"
    ),
    "inchi.jsonl": _schema(
        "error:null|string,inchi:null|string,mol_from_inchi_branches:object,rdkit_ok:boolean,row:integer,smiles:string"
    ),
    "rdkit_builtin_fixture_migration.jsonl": _schema(
        "fixture:string,kind:string",
        "atom_count:integer,atomic_numbers:array,bond_count:integer,bond_types:array,byte_len:integer,canonical_smiles:string,error:null|string,formal_charges:array,nonempty:boolean,rdkit_ok:boolean,record_index:integer,smiles:string",
    ),
    "mol2_read.jsonl": _schema(
        "atoms:null|array,bonds:null|array,case_id:string,chiral_tags:null|array,cleanup_substructures:boolean,error:null|string,fixture:string,mol2:string,positions:null|array,rdkit_null:boolean,rdkit_ok:boolean,remove_hs:boolean,sanitize:boolean,smiles_out:null|object,variant:string"
    ),
    "tetrahedral_stereo_geometry.jsonl": _schema(
        "centers:array,error:null|string,positions:null|array,rdkit_ok:boolean,smiles:string"
    ),
    "assign_atom_chiral_tags_from_structure.jsonl": _schema(
        "after:object,before:object,case_id:string,conf_id:integer,environment:object,error_text:null|string,error_type:null|string,replace_existing_tags:boolean,selected_conformer_id:null|integer,selection_reason:string,status:string"
    ),
    "morgan_fingerprint.jsonl": _schema(
        "branches:object,error:null|string,rdkit_ok:boolean,smiles:string"
    ),
    "atom_pair_fingerprint.jsonl": _schema(
        "branches:object,error:null|string,rdkit_ok:boolean,row:integer,smiles:string"
    ),
    "maccs_fingerprint.jsonl": _schema(
        "error:null|string,label:null|string,public_n_bits:integer,public_on_bits:null|array,raw_n_bits:integer,raw_on_bits:null|array,rdkit_ok:boolean,record_type:string,smiles:string"
    ),
    "rdkit_topological_fingerprint.jsonl": _schema(
        "branches:object,error:null|string,rdkit_ok:boolean,row:integer,smiles:string"
    ),
    "layered_fingerprint.jsonl": _schema(
        "branches:object,error:null|string,rdkit_ok:boolean,row:integer,smiles:string"
    ),
    "topological_torsion_fingerprint.jsonl": _schema(
        "error:null|string,errors:null|object,helpers:null|object,legacy:null|object,profiles:object,rdkit_ok:boolean,row:integer,smiles:string"
    ),
    "avalon_fingerprint.jsonl": _schema(
        "branches:object,row:integer,smiles:string"
    ),
    "molecular_descriptors.jsonl": _schema(
        "descriptor_bits:null|object,descriptor_option_bits:null|object,descriptor_options:null|object,descriptors:null|object,error:null|string,rdkit_ok:boolean,smiles:string",
    ),
    "svg_drawer.jsonl": _schema(
        "error:null|string,height:integer,rdkit_ok:boolean,smiles:string,svg:null|string,width:integer"
    ),
    "prepared_draw_molecule.jsonl": _schema(
        "atoms:null|array,bonds:null|array,error:null|string,rdkit_ok:boolean,smiles:string"
    ),
    "molblock_v2000_minimal.jsonl": _schema(
        "parse_error:null|string,parse_ok:boolean,smiles:string,v2000_body:null|string,v2000_error:null|string,v2000_ok:boolean,v3000_body:null|string,v3000_error:null|string,v3000_ok:boolean"
    ),
    "molblock_v2000_kekulized.jsonl": _schema(
        "kekulize_error:null|string,kekulize_ok:boolean,parse_error:null|string,parse_ok:boolean,smiles:string,v2000_body:null|string,v2000_error:null|string,v2000_ok:boolean,v3000_body:null|string,v3000_error:null|string,v3000_ok:boolean"
    ),
    "molfile_read.jsonl": _schema(
        "case_id:string,dimension:string,error:null|string,format:string,molblock:null|string,rdkit_ok:boolean,smiles:null|string,stereo_markers:string",
        "api:string,atoms:array,bonds:array,chiral_tags:array,operation:string,positions:array,remove_hs:boolean,sanitize:boolean,smiles_out:object,strict_parsing:boolean",
    ),
    "xyz_read.jsonl": _schema(
        "atomic_numbers:null|array,comment:null,coordinates:null|array,error:null|string,num_bonds:null|integer,rdkit_ok:boolean,smiles:string,xyz_block:null|string"
    ),
}


def _json_type(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "boolean"
    if isinstance(value, int):
        return "integer"
    if isinstance(value, float):
        return "number"
    if isinstance(value, str):
        return "string"
    if isinstance(value, list):
        return "array"
    if isinstance(value, dict):
        return "object"
    raise TypeError(f"unsupported JSON value type: {type(value).__name__}")


def _validate_record(output_name: str, line_no: int, record: Any) -> None:
    schema = SCHEMAS.get(output_name)
    if schema is None:
        raise ValueError(f"no expected-data schema is registered for {output_name}")
    if not isinstance(record, dict):
        raise ValueError(f"{output_name}:{line_no}: record must be a JSON object")
    actual = set(record)
    required = set(schema.required)
    allowed = required | set(schema.optional)
    missing = sorted(required - actual)
    unexpected = sorted(actual - allowed)
    if missing:
        raise ValueError(f"{output_name}:{line_no}: missing fields: {', '.join(missing)}")
    if unexpected:
        raise ValueError(
            f"{output_name}:{line_no}: unexpected fields: {', '.join(unexpected)}"
        )
    for name, value in record.items():
        allowed_types = schema.required.get(name) or schema.optional[name]
        actual_type = _json_type(value)
        if actual_type not in allowed_types:
            expected = "|".join(sorted(allowed_types))
            raise ValueError(
                f"{output_name}:{line_no}: field {name!r} has type {actual_type}, "
                f"expected {expected}"
            )


def validate_jsonl_output(path: Path, output_name: str) -> tuple[str, int]:
    """Validate one output and return its exact SHA-256 and record count."""
    if output_name not in SCHEMAS:
        raise ValueError(f"no expected-data schema is registered for {output_name}")
    digest = hashlib.sha256()
    records = 0
    with path.open("rb") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            digest.update(raw_line)
            stripped = raw_line.strip()
            if not stripped or raw_line.lstrip().startswith(b"#"):
                continue
            try:
                text = stripped.decode("utf-8")
            except UnicodeDecodeError as error:
                raise ValueError(f"{output_name}:{line_no}: invalid UTF-8: {error}") from error
            try:
                record = json.loads(text)
            except json.JSONDecodeError as error:
                raise ValueError(f"{output_name}:{line_no}: invalid JSON: {error}") from error
            _validate_record(output_name, line_no, record)
            records += 1
    return digest.hexdigest(), records
