#!/usr/bin/env python3
"""Build a source/evidence audit table from the already recorded 0.3.0 inventory.

This tool never invokes Git, builds, imports a chemistry extension, runs an
oracle, refreshes references, or changes the migration execution ledger.
Public Rust visibility comes from a separately, freshly generated rustdoc.
Manual correspondence rules below are audit proposals, not implemented APIs.
"""

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import importlib.util
import json
from pathlib import Path
import re


BASELINE = "d892ec3507c5b568c5ed5d86ae44e466f7d03855"
LABELS = {
    "alignment": "RMSD／对齐", "batch": "批处理", "bio": "BIO／词汇",
    "confseq": "ConfSeq", "core": "基础分子／编辑／立体化学",
    "depict": "2D／绘图", "descriptors": "描述符", "fingerprints": "指纹",
    "forcefields": "UFF／MMFF", "graph_scaffold_hash": "分片／scaffold／hash",
    "inchi": "InChI", "input_output": "分子 IO／SDF", "metadata": "版本信息",
    "search": "SMARTS／子结构", "serialization_interop": "binary／RDKit 交换",
    "three_d": "3D／距离几何", "support": "参数／结果／协议",
}
STATUS = {
    "public": "[x] 已有公开可编译接口（声明范围）",
    "partial": "[ ] 部分可用／旧合同未齐",
    "domain": "[ ] 仅领域实现／缺公共接入",
    "absent": "[ ] 未找到新版对应实现",
    "unsupported": "[ ] 对应领域边界 Unsupported／公共入口未齐",
    "projection": "[ ] Python 专属投影待迁移",
}
TYPE_NAMES = {
    "MoleculeEdit": "MoleculeBuilder", "PotentialStereoAnalysis": "PotentialStereoResult",
    "SmartsParserParams": "SmartsParseParams", "EmbedParameters": "EmbedParams",
    "ProteinAtom": "ProteinAtomRef", "ProteinChain": "ProteinChainRef",
    "ProteinResidue": "ProteinResidueRef", "StructureAtom": "BioAtomRow",
    "StructureChain": "BioChainRow", "StructureEntity": "BioEntityRow",
    "StructureModel": "BioModelRow", "StructureResidue": "BioResidueRow",
    "MmcifWriteOptions": "BioMmcifWriteParams",
    "TopologicalTorsionFingerprintOptions": "TopologicalTorsionFingerprintParams",
    "SubstructMatchResult": "MatchResult",
}
MOLECULE_NAMES = {
    "analyze_potential_stereo": "potential_stereo",
    "kekulize_": "kekulize_bonds_", "compute_2d_coordinates_": "compute_2d_coordinates_",
    "edit": "to_builder", "mol_from_binary": "from_binary", "mol_to_binary": "to_binary",
    "num_conformers": "num_3d_conformers", "read_mol_from_str": "from_mol",
    "read_mol2_from_str": "from_mol2", "read_sdf_from_str": "from_sdf",
    "from_pdb_block": "from_pdb", "from_mmcif_block": "from_mmcif",
    "from_xyz_block": "from_xyz", "to_pdb_block": "to_pdb",
    "to_2d_sdf_string": "to_sdf_2d", "to_3d_sdf_string": "to_sdf_3d",
    "has_uff_params": "uff_has_all_molecule_params",
    "has_mmff_params": "mmff_has_all_molecule_params",
    "read_sdf": "read_sdf", "write_sdf_to_directory": "write_sdf_files",
    "coordinates_3d": "conformers_3d", "stereoisomers": "enumerate_stereoisomers",
    "largest_fragment": "with_largest_fragment", "murcko_scaffold": "with_murcko_scaffold", "net_scaffold": "with_net_scaffold",
}
GLOBAL_NAMES = {
    "get_element_info": "element_info", "get_residue_info": "residue_info",
    "find_tabulated_residue": "find_residue_info", "find_tabulated_residue_idx": "find_residue_info_index",
    "residue_code_from_name": "residue_code", "element_from_symbol": "Element.from_symbol",
    "expand_protein_one_letter": "expand_one_letter", "expand_protein_one_letter_string": "expand_one_letter_sequence",
    "mol_from_binary": "Molecule.from_binary", "mol_to_binary": "Molecule.to_binary",
    "get_substruct_match": "Molecule.substruct_match", "get_substruct_matches": "Molecule.substruct_matches",
    "get_substruct_matches_with_params": "Molecule.substruct_matches_with_params",
    "has_substruct_match": "Molecule.has_substruct_match",
    "parse_smarts": "search.parse_smarts", "parse_smarts_with_params": "search.parse_smarts_with_params",
    "get_topological_torsion_generator": "TopologicalTorsionFingerprintGenerator.new",
    "topological_torsion_generator_from_json": "TopologicalTorsionFingerprintGenerator.from_json",
    "get_topological_torsion_fingerprint": "Molecule.fingerprint_topological_torsion_sparse_count",
    "get_hashed_topological_torsion_fingerprint": "Molecule.fingerprint_topological_torsion_count",
    "get_hashed_topological_torsion_fingerprint_as_bit_vect": "Molecule.fingerprint_topological_torsion",
    "get_topological_torsion_fingerprint_as_ids": "Molecule.topological_torsion_ids",
    "get_atom_code": "Molecule.atom_pair_atom_code", "py_score_path": "Molecule.topological_torsion_path_score",
    "explain_atom_code": "fingerprints.explain_atom_code", "explain_path_score": "fingerprints.explain_path_score",
    "inchi_to_key": "inchi_to_key",
    "uff_has_all_molecule_params": "Molecule.uff_has_all_molecule_params",
    "mmff_has_all_molecule_params": "Molecule.mmff_has_all_molecule_params",
    "uff_optimize_molecule": "Molecule.with_uff_optimized",
    "uff_optimize_molecule_confs": "Molecule.with_uff_optimized_confs",
    "mmff_optimize_molecule": "Molecule.with_mmff_optimized",
    "mmff_optimize_molecule_confs": "Molecule.with_mmff_optimized_confs",
}
ATOM_NAMES = {
    "idx": "id", "atomic_num": "atomic_number", "atom_map_num": "atom_map",
    "num_radical_electrons": "radical_electrons", "total_num_hs": "total_hydrogens",
}
BOND_NAMES = {
    "idx": "id", "begin_atom_idx": "begin", "end_atom_idx": "end",
    "bond_dir": "direction", "bond_type": "order", "bond_dir_code": "direction_code",
    "bond_dir_name": "direction_name", "bond_type_code": "order_code", "bond_type_name": "order_name",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class RootLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.items = {}

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "a" and attrs.get("class") in ("struct", "enum", "fn", "type", "constant", "mod"):
            href = attrs.get("href", "")
            if re.fullmatch(r"(?:struct|enum|fn|type|constant)\.\w+\.html", href):
                self.items[href.split(".")[1]] = href


def doc_inventory(root, doc_root):
    parser = RootLinks()
    parser.feed((doc_root / "index.html").read_text())
    inventory = {}
    for name, href in parser.items.items():
        path = doc_root / href
        text = path.read_text()
        inventory[name] = {
            "href": str(path.relative_to(root)), "sha256": sha(path),
            "methods": re.findall(r'id="method\.([\w]+)"', text),
            "inherent_methods": re.findall(r'id="method\.([\w]+)"', text.split('id="trait-implementations"')[0]),
            "fields": re.findall(r'id="structfield\.([\w]+)"', text),
            "constants": re.findall(r'id="associatedconstant\.([\w]+)"', text),
            "variants": re.findall(r'id="variant\.([\w]+)"', text),
            "kind": href.split(".")[0],
        }
    return inventory


def evidence_for(docs, path):
    """Return real public documentation anchors, never a guessed method."""
    parts = path.split(".")
    if parts[0] not in docs:
        return None
    entry = docs[parts[0]]
    if len(parts) == 1:
        return entry["href"]
    if len(parts) != 2:
        return None
    if parts[1] in entry["methods"]:
        return entry["href"] + "#method." + parts[1]
    if parts[1] in entry["fields"]:
        return entry["href"] + "#structfield." + parts[1]
    if parts[1] in entry["constants"]:
        return entry["href"] + "#associatedconstant." + parts[1]
    if parts[1] in entry["variants"]:
        return entry["href"] + "#variant." + parts[1]
    return None


def canonical(api):
    if api.startswith("calc_"):
        name = api[5:]
        name = {"mol_wt": "molecular_weight", "exact_mol_wt": "exact_molecular_weight",
                "mol_formula": "molecular_formula", "num_atoms": "total_atom_count"}.get(name, name)
        return "Molecule." + name
    if api in GLOBAL_NAMES:
        return GLOBAL_NAMES[api]
    if "." not in api:
        return TYPE_NAMES.get(api, api)
    owner, name = api.split(".", 1)
    if owner == "Molecule":
        name = MOLECULE_NAMES.get(name, name)
    elif owner == "Atom":
        name = ATOM_NAMES.get(name, name)
    elif owner == "Bond":
        name = BOND_NAMES.get(name, name)
    elif owner == "MoleculeEdit":
        name = {"commit": "build", "set_atom_charge": "set_atom_formal_charge"}.get(name, name)
    elif owner == "PotentialStereoAnalysis":
        name = {"molecule": "cleaned_molecule", "stereo_info": "stereo"}.get(name, name)
    elif owner in ("ProteinAtom",):
        name = {"atomic_num": "atomic_number", "index": "id"}.get(name, name)
    elif owner in ("ProteinChain", "ProteinResidue", "StructureAtom", "StructureChain", "StructureEntity", "StructureModel", "StructureResidue"):
        name = {"index": "id"}.get(name, name)
    elif owner == "TopologicalTorsionFingerprintGenerator" and name.startswith("get_"):
        name = name[4:]
        modes = {"fingerprint": "fingerprint_topological_torsion", "count_fingerprint": "fingerprint_topological_torsion_count",
                 "sparse_fingerprint": "fingerprint_topological_torsion_sparse_bits", "sparse_count_fingerprint": "fingerprint_topological_torsion_sparse_count"}
        if name in modes:
            return "Molecule." + modes[name] + "_with_params"
        if name.endswith("s") and name[:-1] in modes:
            return "MoleculeBatch." + modes[name[:-1]] + "_list_with_params"
        name = {"options": "params"}.get(name, name)
    elif owner == "TopologicalTorsionFingerprintOptions" and name == "set_count_bounds":
        name = "with_count_bounds"
    elif owner == "EmbedParameters" and name == "update_from_json":
        name = "with_json"
    elif owner == "MoleculeBatch":
        name = {"filter_valid": "with_valid_records", "read_sdf_records_from_str": "from_sdf_records"}.get(name, name)
    elif owner == "SparseCountFingerprint":
        name = {"get_value": "value", "size": "length"}.get(name, name)
    elif owner in ("Fingerprint", "SparseBitFingerprint"):
        name = {"n_bits": "length", "size": "length"}.get(name, name)
    if name == "__new__":
        name = "new"
    if owner in ("BioStructure", "Protein"):
        name = {"from_pdb_str": "from_pdb", "from_mmcif_str": "from_mmcif",
                "from_structure_str": "from_text", "from_pdb": "read_with_format", "from_mmcif": "read_with_format"}.get(name, name)
    return TYPE_NAMES.get(owner, owner) + "." + name


def old_declaration(root, current, source_cache):
    """Existing legacy binding header, explicitly not a newly fetched baseline."""
    if not current:
        return "当前未找到同名旧声明；未重新读取历史签名"
    path = root / current["source"]
    if path not in source_cache:
        source_cache[path] = path.read_text().splitlines()
    source = source_cache[path]
    n = current["lines"][0]
    snippet = []
    for text in source[n - 1:n + 22]:
        snippet.append(text.strip())
        if "{" in text or text.rstrip().endswith(";"):
            break
    return " ".join(snippet)[:1600]


def current_match(root, docs, api, kind, new, family, counts):
    """Source-reviewed correspondences; 'public' means callable, not parity."""
    note, paths = [], [new]
    state = None
    if api.startswith("calc_"):
        ordinary = counts["descriptor_breakdown"]
        if api in ordinary["public_rust_counterparts"]:
            name = ordinary["public_rust_counterparts"][api]
            paths = ["Molecule." + name]
            if name.endswith("_with_options"):
                paths.insert(0, new)
                note.append("默认方法已公开；现存 _with_options 需收敛参数对象／_with_params；选项不能丢失")
            if api == "calc_num_atoms":
                note.append("含隐式 H 的 descriptor；绝不映射到仅图行数的 Molecule.num_atoms")
            if api not in ("calc_mol_wt", "calc_exact_mol_wt", "calc_mol_formula", "calc_num_heavy_atoms", "calc_lipinski_hba", "calc_num_heteroatoms", "calc_num_atoms"):
                note.append("按公开合同需要已有 valence 或 ring 状态；裸构造／准备分支仍需与基线核对")
        elif api in ordinary["private_or_partial_domain_implementation_without_public_facade"]:
            state = "domain"
            name = new.split(".")[-1]
            if name.startswith(("slogp_vsa_", "smr_vsa_")):
                name = name.rsplit("_", 1)[0]
                note.append("逐箱标量必须从同一 VSA 数组结果读取，不能复制算法")
            name = {"crippen_descriptors": "crippen_totals", "labute_asa_contributions": "labute_contributions",
                    "tpsa": "tpsa_contributions"}.get(name, name)
            paths = ["cosmolkit_descriptors::" + name]
            note.append("领域层存在全部或部分 primitive；公共 Molecule 组装／缓存合同尚未接通")
        else:
            state, paths = "absent", []
            note.append("Chi／Hall–Kier／Kappa／Phi／QED／MQN 的新版实现未找到；旧绑定声明不计实现")
    elif family in ("alignment", "batch", "confseq", "serialization_interop", "graph_scaffold_hash", "inchi"):
        state = "unsupported" if family in ("alignment", "batch") and api.startswith(("Molecule.", "MoleculeBatch.")) else "absent"
        paths = {"alignment": ["cosmolkit_alignment::align"], "batch": ["cosmolkit_batch::process"],
                 "graph_scaffold_hash": ["cosmolkit_core::connected_components（仅 fragments primitive）"],
                 "confseq": [], "serialization_interop": [], "inchi": []}[family]
        if family in ("alignment", "batch") and state == "absent":
            paths = []
        if api in ("AlignmentResult.atom_map", "AlignmentResult.rmsd"):
            state, paths = "domain", ["cosmolkit_alignment::AlignmentResult::" + ("atom_pairs" if api.endswith("atom_map") else "rmsd")]
        if family == "graph_scaffold_hash":
            paths = ["cosmolkit_core::connected_components"] if api in ("Molecule.fragments", "Molecule.largest_fragment") else []
        if family == "batch" and api.split(".")[-1] in ("errors", "invalid_count", "invalid_mask", "parallel_jobs", "progress_bar", "to_list", "valid_count", "valid_mask", "with_parallel_jobs", "with_progress_bar", "filter_valid"):
            state, paths = "absent", []
        note.append({"alignment": "MAIN align 固定 Unsupported；结果只建模 atom_pairs/rmsd，缺 transform／conformer report",
                     "batch": "MAIN process 固定 Unsupported；结果读取／配置也缺唯一 MoleculeBatch 生命周期",
                     "confseq": "四个 decode 入口无新版 canonical owner",
                     "serialization_interop": "binary 旧格式／pickle／RDKit 交换公共兼容链未接通",
                     "graph_scaffold_hash": "已有图 primitive 不等于 fragments／scaffold／分子 hash 完整公开行为",
                     "inchi": "隔离审计／局部查找修复不等于可用构造、序列化或 InChIKey"}[family])
    elif family == "forcefields":
        state, paths = "absent", []
        if api.startswith("Molecule.") or "." not in api:
            state = "unsupported"
            name = "uff_has_all_molecule_params" if "uff" in api.lower() else ("mmff_optimize" if "optimize" in api.lower() else "mmff_has_all_molecule_params")
            paths = ["cosmolkit_forcefields::" + name]
            if "uff" in api.lower() and "optimize" in api.lower():
                state, paths = "absent", []
        note.append("MAIN 参数／优化入口未闭合；_1 UFF 交付另列，不等于 MAIN 或完整 MMFF 可用")
    elif family == "three_d":
        if api == "Molecule.coordinates_3d":
            paths, state = ["Molecule.conformers_3d", "Conformer3D.positions"], "partial"
            note.append("全部 3D 值可读；原选择 conformer_index／NumPy shape 的桥未实现，不能只取第一组代替")
        elif api == "Molecule.num_conformers":
            paths, state = ["Molecule.conformers_3d"], "partial"
            note.append("可借用 3D 列表并 len；目前没有对应 num_3d_conformers 方法")
        elif api in ("Molecule.with_3d_coordinates", "Molecule.with_added_3d_conformer", "Molecule.add_3d_conformer_"):
            paths, state = ["Molecule.to_builder", "MoleculeBuilder.add_3d_conformer", "MoleculeBuilder.build"], "partial"
            note.append("builder 的 detached 坐标构造已可用；旧值／原地操作名称与 contract 尚未接通")
        elif api.startswith("Molecule.") and any(term in api for term in ("embed_", "with_3d_conformer")):
            state, paths = "unsupported", ["cosmolkit_conformer::embed", "cosmolkit_conformer::embed_multiple"]
            note.append("embedding 固定 Unsupported；手动坐标、参数／结果值不能替代生成／bounds／选择更新链")
        else:
            state, paths = "absent", []
            note.append("未找到对应 bounds／坐标更新／EmbedParams预设或完整结果 schema；ConformerOptions只有max_conformers")
    elif api.startswith("Atom."):
        member = api.split(".")[1]
        if member in ("chiral_tag_code", "chiral_tag_name"):
            paths = ["Atom.chiral_tag", "ChiralTag.rdkit_" + member.rsplit("_", 1)[1]]
            note.append("代码／名称来自同一枚举，不要求建立平行映射表")
        elif member == "degree":
            state, paths = "partial", ["TopologyBlock.adjacency", "AdjacencyList.degree"]
            note.append("degree 需 topology 上下文；Atom 值本身没有该方法；待确定唯一公开读取投影")
        elif member in ("explicit_valence", "implicit_hydrogens", "total_num_hs", "total_valence"):
            state, paths = "domain", ["cosmolkit_core::assign_valence", "cosmolkit_core::total_hydrogen_count_from_validated"]
            note.append("旧 Atom 是含 prepared 字段的快照；新版 detached Atom 不持有有效 valence，需公开读取投影")
        elif member in ("cip_neighbor_order", "cip_rank"):
            state, paths = "partial", ["Atom.prop"]
            note.append("通用 typed property 存储可读；专用返回类型／computed／缺失语义的投影未接通")
        elif member == "idx":
            paths = ["Atom.id", "AtomId.index"]
    elif api.startswith("Bond."):
        member = api.split(".")[1]
        if member in ("bond_dir_code", "bond_dir_name", "bond_type_code", "bond_type_name", "stereo_code", "stereo_name"):
            getter, enum = ("direction", "BondDirection") if "dir" in member else (("order", "BondOrder") if "type" in member else ("stereo", "BondStereo"))
            paths = ["Bond." + getter, enum + ".rdkit_" + member.rsplit("_", 1)[1]]
        elif member == "cip_neighbor_order":
            state, paths = "partial", ["Bond.prop"]
            note.append("专用 CIP 顺序读取投影尚未接通；存储字段不能替代该返回合同")
        elif member in ("idx", "begin_atom_idx", "end_atom_idx"):
            paths = [new, ("BondId.index" if member == "idx" else "AtomId.index")]
    elif api.startswith("MoleculeEdit."):
        paths = [new]
        note.append("唯一 builder 路径；commit 映射 owned build，不能新增兼容 MoleculeEdit owner")
    elif family == "core" and api.startswith("Molecule."):
        if api in ("Molecule.stereoisomers", "Molecule.stereoisomer_count"):
            state, paths = "unsupported", ["cosmolkit_stereo::enumerate"]
            note.append("当前 stereo 枚举 owner 未 detached；基本 stereo／CIP 不替代枚举")
        elif api in ("Molecule.find_chiral_centers", "Molecule.tetrahedral_stereo", "Molecule.perceive_stereochemistry"):
            state, paths = "domain", ["cosmolkit_core::assign_legacy_stereochemistry", "cosmolkit_core::potential_stereo"]
            note.append("已有基本 stereo 算法；旧返回结构／参数／公共读取或变换调用链未接通")
        elif api == "Molecule.cip_computed":
            state, paths = "partial", ["Molecule.atoms", "Atom.cip_descriptor"]
            note.append("标签存在不等于原 cip_computed 的有效性合同；没有同名公开 query")
        elif api == "Molecule.set_2d_coordinates_":
            state, paths = "partial", ["MoleculeBuilder.set_2d_coordinates"]
            note.append("builder 可构造；尚无 live Molecule 原地 coordinate contract")
        elif api == "Molecule.sanitize_":
            state, paths = "partial", ["Molecule.sanitize", "Molecule.sanitize_with_params"]
            note.append("已有 value form；当前声明没有生成 in-place sanitize_，不能将返回新值算原地完成")
        elif api == "Molecule.analyze_potential_stereo":
            note.append("新版结果为 PotentialStereoResult；字段映射与 clean 参数／返回状态需保留")
    elif family == "depict":
        if api == "Molecule.with_2d_coordinates":
            state = "partial"
            note.append("MAIN 公开 value layout 可调用；现存 MAIN 5000报告仍208失败；_2的坐标修复未合入且完整状态红")
        elif api == "Molecule.has_2d_coordinates":
            state, paths = "partial", ["Molecule.coordinates_2d"]
            note.append("is_some 可作存在性基础；旧 empty-conformer 存在判定须核对，不能把空数组当不存在")
        elif api != "Molecule.coordinates_2d":
            state, paths = "absent", []
            note.append("MAIN 尚无该公共入口；_4绘图实现／Python选定profile另列")
    elif api.startswith(("BioStructure.", "Protein.")):
        owner, member = api.split(".")
        if member in ("from_pdb", "from_mmcif"):
            state, paths = "partial", [owner + ".read_with_format"]
            note.append("旧短名接收文件路径；新版短名接收文本；文件行为对应 read_with_format，不能静默改变输入解释")
        elif member in ("to_mmcif", "write_mmcif"):
            state = "partial"
            note.append("只支持坐标 profile；旧33组开关中28组缺公共组合输出／默认完整 categories")
        elif member == "to_molecule":
            state, paths = "absent", []
            note.append("BIO→化学 Molecule 桥未接通；不得把结构行直接当分子算法")
    elif api.startswith("ProteinAtom."):
        member = api.split(".")[1]
        paths = {"index": ["ProteinAtomRef.id", "BioAtomId.index"],
                 "atomic_num": ["ProteinAtomRef.element", "Element.atomic_number"],
                 "element_symbol": ["ProteinAtomRef.element", "Element.symbol"]}.get(member, [new])
        if member == "name":
            paths += ["AtomName.as_str"]
            note.append("AtomName 的借用字符串需保留旧名称去空格策略")
    elif api.startswith("ProteinChain."):
        member = api.split(".")[1]
        paths = [new, "BioChainId.index"] if member == "index" else [new]
    elif api.startswith("ProteinResidue."):
        member = api.split(".")[1]
        if member == "index":
            paths = [new, "BioResidueId.index"]
        elif member in ("canonical_one_letter_code", "parent_standard_code", "is_modified_amino_acid"):
            paths = ["ProteinResidueRef.info", "ResidueInfo." + member]
    elif api.startswith("Structure") and "." in api:
        owner, member = api.split(".")
        # Rows intentionally do not manufacture live hierarchy view objects.
        direct = {
            "StructureAtom": {"altloc": "altloc", "b_factor": "b_iso", "element": "element", "formal_charge": "formal_charge", "name": "name", "occupancy": "occupancy", "residue_index": "residue_id"},
            "StructureChain": {"entity_index": "entity_id", "kind": "kind", "model_index": "model_id"},
            "StructureEntity": {"kind": "kind", "polymer_kind": "polymer_kind", "sequence": "full_sequence", "subchains": "subchains"},
            "StructureModel": {"source_model_number": "source_model_number"},
            "StructureResidue": {"chain_index": "chain_id", "entity_kind": "entity_kind", "kind": "kind", "name": "name"},
        }
        if member in direct.get(owner, {}):
            paths = [TYPE_NAMES[owner] + "." + direct[owner][member]]
            if member.endswith("_index"):
                note.append("以具名 ID 读取；Python需显式投影 index／None")
        elif api == "StructureAtom.element_symbol":
            paths = ["BioAtomRow.element", "Element.symbol"]
        elif api == "StructureAtom.position":
            paths = ["BioStructure.atom_position"]
            note.append("坐标在结构独立 block；需保留所属 structure 与 BioAtomId 上下文")
        elif api == "StructureAtom.source_serial":
            paths = ["BioAtomRow.source", "AtomSourceIds.serial"]
        elif owner == "StructureChain" and member in ("auth_chain_id", "label_asym_id"):
            paths = ["BioChainRow.source", "ChainSourceIds." + member]
        elif api == "StructureEntity.source_id":
            paths = ["BioEntityRow.source", "EntitySourceIds.source_entity_id"]
        elif owner == "StructureResidue" and member in ("source_sequence_number", "insertion_code"):
            paths = ["BioResidueRow.source", "ResidueSourceIds.seq_id", "PdbSeqId." + ("seq_num" if member == "source_sequence_number" else "ins_code")]
        elif owner == "StructureResidue" and member in ("code", "info"):
            state, paths = "partial", ["BioResidueRow.name", "find_residue_info"]
            note.append("可组合现有 residue owner；尚无旧 view 的专用方法，须保持 unknown/name/kind 语义")
        else:
            state, paths = "partial", [TYPE_NAMES[owner]]
            note.append("行与 span 已有；旧 view 的 index／父子遍历需共享上下文投影，不是同名 row 方法")
    elif api.startswith("ResidueInfo."):
        if api.endswith(".kind_name"):
            state, paths = "projection", ["ResidueInfo.kind"]
            note.append("Python稳定显示名／枚举投影尚未迁移；Debug 文本不能替代合同")
    elif api in ("get_element_info",):
        state, paths = "domain", ["cosmolkit_core::element_info"]
        note.append("ElementInfo schema 已 re-export；产生元数据的 algorithm 函数未由顶层 re-export")
    elif api in ("expand_protein_one_letter", "expand_protein_one_letter_string"):
        state, paths = "partial", ["expand_one_letter" if api.endswith("one_letter") else "expand_one_letter_sequence"]
        note.append("旧蛋白专用入口隐含 ResidueInfoKind；规范入口需显式 kind，默认与错误投影待保留；不新增兼容别名")
    elif family == "fingerprints":
        if api.startswith("SparseCountFingerprint."):
            paths = [new]
            if api.endswith(".get_value"):
                note.append("与 value 合并为一个 canonical 入口；保留越界 Result 错误，不加兼容别名")
        elif api.startswith(("Fingerprint.", "SparseBitFingerprint.")):
            state = "domain"
            actual = {"Fingerprint.from_on_bits": "Fingerprint::from_on_bits", "Fingerprint.n_bits": "Fingerprint::n_bits",
                      "Fingerprint.on_bits": "Fingerprint::on_bits", "Fingerprint.tanimoto": "similarity::tanimoto",
                      "SparseBitFingerprint.on_bits": "SparseBitFingerprint::on_bits", "SparseBitFingerprint.size": "SparseBitFingerprint::n_bits"}
            paths = ["cosmolkit_fingerprints::" + actual[api]]
            note.append("值层存在；MAIN 未重导出 dense／sparse-bit 类型和相似度公共方法；现存n_bits与目标length需明确适配")
        else:
            state, paths = "absent", []
            note.append("MAIN 未找到生成器／additional output／完整结果链；值层不等于生成器已实现；Morgan _2另列")
            if api.startswith("TopologicalTorsionFingerprintGenerator.get_"):
                note.append("单分子生成改为Molecule receiver＋不可变参数对象；多分子改为MoleculeBatch；旧共享可变generator不是规范公共owner")
            if api == "TopologicalTorsionFingerprintOptions.set_count_bounds":
                note.append("参数必须不可变；with_count_bounds返回新配置，不保留旧共享可变generator引用")
    elif family == "search":
        state = "domain"
        name = api.split(".")[-1]
        if api.startswith("SubstructMatchResult."):
            paths = ["cosmolkit_search::MatchResult::" + ("atom_mapping" if name == "atom_pairs" else name)]
            if name == "atom_pairs":
                note.append("原 pairs 可由 query位置 enumerate(atom_mapping) 投影，领域值没有 atom_pairs 方法")
        else:
            name = {"parse_smarts_with_params": "parse_smarts", "to_smarts": "query_graph_to_smarts", "to_cx_smarts": "query_graph_to_cx_smarts",
                    "get_substruct_matches_with_params": "try_get_substruct_matches_with_params"}.get(name, name)
            paths = ["cosmolkit_search::" + name]
        note.append("detached QueryGraph parser／writer／matcher已存在，顶层公共组装／Python未接通")
        if api.startswith("parse_smarts"):
            note.append("旧返回 query-bearing Molecule；目标返回唯一 QueryGraph，此为架构所需类型迁移，不能偷偷 concrete 化")
        if api.startswith(("get_substruct", "has_substruct")):
            note.append("domain普通包装存在吞错/空结果路径；接入必须使用结构化Result并保留错误，不能直接照搬该fallback")
    elif family == "input_output":
        if api == "Molecule.read_sdf_from_str":
            paths = ["Molecule.from_sdf"]
            note.append("旧接口返回单个Molecule；新版同为单记录构造。旧sanitize/remove_hs/strict参数与缺失错误仍须核对")
            note.append("新版Molecule构造限具体图；query记录须用SdfRecord.graph/SdfGraph::Query，不能把旧query-bearing Molecule原样保留")
        elif api in ("Molecule.read_mol_from_str",):
            state, paths = "domain", ["cosmolkit_io::read_mol_block_detached_with_params"]
            note.append("molfile-only读取会忽略末尾SDF字段；不能以from_sdf替代该合同")
        elif api.startswith("SdfRecord."):
            member = api.split(".")[1]
            if member == "title":
                paths = ["SdfRecord.properties", "MoleculeProperties.name"]
            elif member == "molecule":
                note.append("新版返回Result<&Molecule>；query记录为WrongGraphKind，须走SdfGraph::Query；Python拥有值/错误投影待接通")
            elif member == "data_field":
                state, paths = "partial", ["SdfRecord.data_fields"]
                note.append("有顺序字段列表，缺data_field wrapper；现存旧绑定取首个同名键、缺失为None，不能改为最后键覆盖")
            elif member == "index":
                state, paths = "domain", ["cosmolkit_io::SdfRecordMetadata"]
                note.append("新公共 SdfRecord 没有流记录 index；需 reader provenance")
        else:
            state = "domain"
            actual = {
                "Molecule.from_pdb_block": "read_pdb_detached_with_params", "Molecule.from_xyz_block": "read_xyz_detached",
                "Molecule.read_mol": "read_mol_block_detached_with_params", "Molecule.read_mol2": "read_mol2_detached_with_params",
                "Molecule.read_mol2_from_str": "read_mol2_detached_with_params", "Molecule.read_sdf": "read_sdf_record_detached_with_params",
                "Molecule.to_2d_sdf_string": "write_sdf_record_detached", "Molecule.to_3d_sdf_string": "write_sdf_record_detached",
                "Molecule.to_pdb_block": "write_pdb_detached_with_params", "Molecule.write_sdf": "write_sdf_record_detached",
                "Molecule.write_sdf_to_directory": "write_sdf_record_detached",
                "SdfDataset.open": "SdfGraphDataset::open", "SdfDataset.metadata": "SdfGraphDataset::metadata",
                "SdfDataset.path": "SdfGraphDataset::path", "SdfReader.open": "SdfGraphReader::new",
            }
            if api.startswith("SdfRecordMetadata."):
                member = api.split(".")[1]
                fields = {"byte_range": ["byte_offset", "byte_len"], "line_range": ["line_offset", "line_len"]}.get(member, [member])
                paths = ["cosmolkit_io::SdfRecordMetadata::" + field for field in fields]
                note.append("范围须组合offset/len；metadata未由facade公开")
            elif api in actual:
                paths = ["cosmolkit_io::" + actual[api]]
            else:
                state, paths = "absent", []
            if api == "Molecule.from_mmcif_block":
                note.append("旧结构→化学Molecule转换依赖sanitize/remove_hs/flavor/proximity_bonding；BIO mmCIF结构解析不能替代该桥")
            if api == "Molecule.read_sdf":
                note.append("旧文件构造只返回首个Molecule，不能改名为read_sdf_records或改为列表")
            if api == "SdfReader.open":
                note.append("domain new需要BufRead；缺文件open与live分子包装")
            note.append("已有 detached IO／reader／dataset；live Molecule、文件／流／批次 wrapper未完整公共接入")
    if state is None:
        anchors = [evidence_for(docs, path) for path in paths]
        state = "public" if paths and all(anchors) else "partial"
        if state == "partial":
            note.append("存在基础值或候选组合，但对应公开成员未齐；不以名字匹配确认可用")
    if "__" in api and kind == "protocol":
        state = "projection"
        note.append("Python协议需要专门行为投影；Rust Clone／Debug／Iterator不自动证明旧协议")
    return state, paths, note


def support_match(docs, api, kind, new):
    paths, state, notes = [new], None, []
    if kind == "class":
        state = "public" if evidence_for(docs, new) else "absent"
        notes.append("类型 schema 存在不意味着其生产算法、字段或 Python 类已交付")
    elif kind == "constructor":
        state = "public" if evidence_for(docs, new) else "absent"
        notes.append("构造 defaults／校验／错误与旧版须逐项保留；new是Rust关联构造的命名")
    elif api.startswith("PotentialStereoAnalysis."):
        member = api.split(".")[1]
        paths = ["PotentialStereoResult." + {"molecule": "cleaned_molecule", "stereo_info": "stereo"}.get(member, member)]
        state = "partial"
        notes.append("新版 cleaned_molecule 是 Option；旧 analysis molecule／clean语义不直接等价")
    elif api.startswith("PotentialStereoInfo."):
        member = api.split(".")[1]
        paths = ["PotentialStereoInfo." + {"center_index": "centered_on", "center_kind": "centered_on"}.get(member, member)]
        notes.append("centered_on 是 AtomId／BondId 变体；旧 scalar 字段需显式投影")
    elif api.startswith("MmcifOutputGroups."):
        state, paths = "partial", ["BioMmcifWriteParams"]
        member = api.split(".")[1]
        notes.append("旧完整33组开关；" + ("坐标profile有部分对应，非完整可组合选项" if member in ("atoms", "auth_all", "block_name", "entry", "group_pdb") else "该组缺公共组合输出；私有 emitter 不能打勾"))
    elif api.startswith("MmcifWriteOptions."):
        state, paths = "partial", ["BioMmcifWriteParams"]
        notes.append("旧排版与组配置不是当前 coordinate-only 参数 schema；默认和组合行为未齐")
    elif api.startswith("SmartsParserParams."):
        state, paths = "domain", ["cosmolkit_search::SmartsParseParams"]
        notes.append("参数 owner 已有；顶层尚未 re-export／注册，Python类仍旧桥")
    elif api.startswith(("EmbedParameters.", "TopologicalTorsionFingerprintOptions.")):
        state, paths = "absent", []
        notes.append("规范参数对象不可变；旧setter共享generator/原地JSON更新需改成新配置值，不直接继承旧写入协议")
    elif kind == "class_attribute":
        state, paths = "absent", []
        notes.append("旧 atom-pair 编码常量公开合同需随唯一指纹 owner 投影，不能以表名存在算完成")
    elif kind == "protocol":
        state, paths = "projection", []
        notes.append("Python专属 repr／len／iter／getitem／pickle行为；需固定回归，Rust trait不是旧协议证明")
    else:
        state = "public" if evidence_for(docs, new) else "absent"
        notes.append("参数／结果字段缺口单列，不把字段数量当算法数量")
    if state is None:
        state = "public" if all(evidence_for(docs, path) for path in paths) else "partial"
    return state, paths, notes


def branch_note(api):
    if api in ("Molecule.fingerprint_morgan", "Molecule.fingerprint_morgan_with_output") or api.startswith(("MorganAdditionalOutput", "MorganFingerprintResult")):
        return "_2：morgan.rs 公共四路生成／additional output交付；未整合 MAIN", "../COSMolKit_2/crates/cosmolkit/src/morgan.rs"
    if "uff" in api.lower():
        return "_1：UFF公共链已有；一步5000仍38错误行，不是整包通过；未整合 MAIN", "../COSMolKit_1/dev/gap_reports/three_d_migration/parity_framework_migration.md"
    if api in ("Molecule.from_smiles", "Molecule.to_smiles", "Molecule.coordinates_2d", "Molecule.num_atoms", "Molecule.num_bonds", "Molecule.with_2d_coordinates", "Molecule.to_svg", "Molecule.to_png", "Molecule.write_svg", "Molecule.write_png"):
        return "_4：选定 drawing Python profile已交付；仅num_atoms/num_bonds无已知旧接口差异，其余参数／返回值／范围仍需核对；未整合 MAIN", "../COSMolKit_4/python/src/drawing_binding.rs"
    if api in ("Molecule.from_inchi", "Molecule.to_inchi", "Molecule.to_inchi_key", "inchi_to_key"):
        return "_3：隔离审计与局部修复；405 finding未闭合，非公共可用", "../COSMolKit_3/dev/gap_reports/inchi/element_lookup_repair.md"
    return "", ""


def corpus_tasks(api):
    name = api.split(".")[-1]
    special = {
        "calc_mol_wt": "molecular_weight", "calc_exact_mol_wt": "exact_molecular_weight",
        "calc_mol_formula": "molecular_formula", "calc_num_atoms": "total_atom_count",
        "Molecule.from_smiles": "smiles_read", "Molecule.with_hydrogens": "add_hydrogens",
        "Molecule.add_hydrogens_": "add_hydrogens", "Molecule.without_hydrogens": "remove_hydrogens",
        "Molecule.remove_hydrogens_": "remove_hydrogens", "Molecule.with_kekulized_bonds": "kekulize",
        "Molecule.kekulize_": "kekulize", "Molecule.with_2d_coordinates": "coordinates_2d",
        "Molecule.uff_has_all_molecule_params": "uff_has_all_molecule_params",
        "Molecule.has_uff_params": "uff_has_all_molecule_params", "uff_has_all_molecule_params": "uff_has_all_molecule_params",
        "Molecule.with_uff_optimized": "uff_optimize", "uff_optimize_molecule": "uff_optimize",
        "Molecule.with_uff_optimized_confs": "uff_optimize_conformers", "uff_optimize_molecule_confs": "uff_optimize_conformers",
        "Molecule.fingerprint_morgan": "morgan", "Molecule.fingerprint_morgan_with_output": "morgan",
    }
    name = special.get(api, name.removeprefix("calc_"))
    return [name + "_smiles"]


def historical_reports(root, audit_dir):
    summaries = []
    paths = [root / "target/parity-stage-final-5000/rust-report.json"]
    paths += list((root / "target/parity-tests").glob("*/rust-report.json"))
    # Only named supplemental folders; no broad target-tree or stale-binary inference.
    for sibling, folder in (("COSMolKit_1", "uff-one-step-5000"), ("COSMolKit_2", "morgan-separated-5000"),
                            ("COSMolKit_2", "d2-pair-current-5000"), ("COSMolKit_2", "d2-bcs-current-5000")):
        paths += [root.parent / sibling / "target/parity-tests" / folder / "rust-report.json"]
    for path in paths:
        if not path.exists():
            continue
        report = json.loads(path.read_text())
        for test in report.get("tests", []):
            details = Path(test["details"])
            record = {"aggregate": str(path), "aggregate_sha256": sha(path), **test,
                      "implementation_identity": "历史报告未绑定当前完整生产源树；不可自动继承当前验收"}
            if details.exists():
                rows = json.loads(details.read_text())
                if not isinstance(rows, list):
                    raise ValueError("unexpected comparison details: " + str(details))
                ids = {row["label"]["case_id"] for row in rows}
                failed = sum(row.get("matches") is not True for row in rows)
                record.update({"observed_rows": len(rows), "unique_case_ids": len(ids),
                               "observed_failures": failed, "details_sha256": sha(details),
                               "counts_verified": len(rows) == test["compared"] and failed == test["failed"],
                               "first_label": rows[0]["label"] if rows else None})
            else:
                record.update({"counts_verified": False, "missing_details": True})
            summaries.append(record)
    (audit_dir / "historical_corpus_evidence.json").write_text(json.dumps(summaries, ensure_ascii=False, indent=2) + "\n")
    return summaries


def corpus_status(api, kind, family, records):
    if kind not in ("method", "module_function") or ("." in api and not api.startswith(("Molecule.", "MoleculeBatch.", "confseq."))):
        return "测试未定义（值／参数／结果／协议；应做固定回归）", []
    if family == "metadata":
        return "测试未定义（版本元数据；无需5000化学语料）", []
    matching = [r for r in records if r["test"] in corpus_tasks(api) and r.get("unique_case_ids") == 5000 and r.get("counts_verified")]
    if matching:
        main = [r for r in matching if str(r["aggregate"]).startswith(str(Path.cwd() / "target"))]
        selected = main or matching
        failed = [r for r in selected if r["failed"]]
        desc = "; ".join(f"{r['test']}: {r['observed_rows']}比较／{r['unique_case_ids']}输入／{r['failed']}失败" for r in selected)
        state = "[ ] 历史未通过；当前待验收" if failed else "[ ] 历史Rust通过；当前／Python待验收"
        return state + "（" + desc + "）", [r["details"] for r in selected]
    if api in ("Molecule.to_svg", "Molecule.to_png", "Molecule.write_svg", "Molecule.write_png"):
        return "[ ] 未验收（_4只有固定binding回归，非5000图像语料）", []
    return "[ ] 未验收（未找到可核对的对应5000完整结果；待定义／运行）", []


def markdown_table(rows, columns):
    def cell(value):
        return str(value).replace("|", "\\|").replace("\n", " ")
    return "\n".join(["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"] +
                      ["| " + " | ".join(cell(row.get(c, "")) for c in columns) + " |" for row in rows])


def implementation_sources(root, paths, family):
    """Locate production definitions, rather than print invented domain symbols."""
    evidence = []
    for item in paths:
        if "::" not in item:
            continue
        package, *members = item.split("::")
        owner = root / "crates" / package.replace("_", "-") / "src"
        if not owner.exists():
            continue
        member = members[-1]
        preferred = None
        if len(members) >= 2:
            preferred = {("cosmolkit_io", "SdfGraphReader"): "sdf.rs", ("cosmolkit_io", "SdfGraphDataset"): "sdf.rs",
                         ("cosmolkit_io", "SdfRecordMetadata"): "sdf.rs", ("cosmolkit_fingerprints", "Fingerprint"): "values.rs",
                         ("cosmolkit_fingerprints", "SparseBitFingerprint"): "sparse_bits.rs"}.get((package, members[0]))
        pattern = re.compile(r"\bpub(?:\([^)]*\))?\s+(?:(?:const|async|unsafe)\s+)?(?:fn|struct|enum|type)\s+" + re.escape(member) + r"\b|\bpub\s+" + re.escape(member) + r"\s*:")
        for file in sorted(owner.rglob("*.rs")):
            if preferred and file.name != preferred:
                continue
            for i, line in enumerate(file.read_text().splitlines(), 1):
                if pattern.search(line):
                    evidence.append(str(file.relative_to(root)) + ":" + str(i))
    if not evidence:
        owner = {"alignment": "alignment", "batch": "batch", "confseq": "confseq", "forcefields": "forcefields",
                 "three_d": "conformer", "inchi": "inchi", "search": "search", "fingerprints": "fingerprints",
                 "input_output": "io", "descriptors": "descriptors", "graph_scaffold_hash": "core",
                 "bio": "bio", "depict": "depict", "core": "core", "serialization_interop": "io"}.get(family)
        if owner:
            path = root / "crates" / ("cosmolkit-" + owner) / "src/lib.rs"
            if path.exists():
                evidence = [str(path.relative_to(root)) + "（owner源码检索范围；不是对应函数存在证明）"]
    return list(dict.fromkeys(evidence))


def dynamic_surface(root, docs):
    """Supplement items missed by the inherited static PyO3 slot parser.

    Detail provenance is the retained legacy code, not a newly retrieved commit.
    The prior audit establishes eight enums/four InChI classes; additional
    members and BatchValidationError remain baseline-detail candidates.
    """
    rows = []
    source = (root / "python/src/lib.rs").read_text().splitlines()
    specs = [(name, name, [name], "动态词汇；旧IntEnum数值、名称、别名须独立投影核对") for name in
             ("Element", "BondOrder", "BondDirection", "BondStereo", "ChiralTag", "BatchErrorMode", "ResidueCode", "ResidueInfoKind")]
    maps = {"ELEMENT_MAP": "Element", "BOND_ORDER_MAP": "BondOrder", "BOND_DIRECTION_MAP": "BondDirection",
            "BOND_STEREO_MAP": "BondStereo", "CHIRAL_TAG_MAP": "ChiralTag", "BATCH_ERROR_MODE_MAP": "BatchErrorMode",
            "RESIDUE_CODE_MAP": "ResidueCode", "RESIDUE_INFO_KIND_MAP": "ResidueInfoKind"}
    specs += [(name, name, [owner], "Python只读MappingProxy；需保留MAP键与别名；不新增Rust平行映射算法") for name, owner in maps.items()]
    specs += [(name, name, [], "动态异常；继承与结构化错误字段投影尚未恢复") for name in
              ("BatchValidationError", "InchiError", "InchiAllocationError", "InchiUnsupportedStateError", "InchiDiagnosticWarning")]
    fields = {"BatchValidationError": ("error_count", "reason", "errors"), "InchiError": ("operation", "kind", "detail"),
              "InchiDiagnosticWarning": ("level", "message")}
    for owner, members in fields.items():
        specs += [(owner + "." + member, owner + "." + member, [], "动态结果/错误投影；errors必须返回副本；固定错误回归，不做5000独立算法测试") for member in members]
        specs += [(owner + ".__init__", owner + ".new", [], "Python异常构造；输入、defaults、继承与字段固定回归")]
    specs += [("__version__", "__version__（Python）／version（Rust）", ["version"], "Python模块元数据；不作为化学算法计数"),
              ("confseq（模块）", "confseq（模块）", [], "子模块；4个decode普通入口已包含在551项，不重复加算法")]
    for old, new, paths, note in specs:
        matching = [str(i) for i, line in enumerate(source, 1) if ('"' + old + '"') in line or (old.split(".")[0] in line and line.lstrip().startswith("class "))]
        enum_available = len(paths) == 1 and evidence_for(docs, paths[0])
        state = "public" if old == "__version__" and enum_available else ("partial" if enum_available else "projection")
        rows.append({"原始命名": old, "新版命名": new, "Rust接口是否实际可用": STATUS[state],
                     "Python接口是否实际可用": "[ ] MAIN默认包664编译错误，动态投影不可验收",
                     "是否通过5000语料测试": "测试未定义（模块／词汇／错误固定回归）",
                     "当前Rust对应": "; ".join(paths) or "无已确认对应公开投影",
                     "缺口与限制": note, "现存声明证据": "python/src/lib.rs:" + (matching[0] if matching else "932–987"),
                     "基线证据界限": "前ROOT确认8枚举/4个InChI类；细项按现存旧绑定补列，指定提交细项待重新核对",
                     "Rust状态码": state})
    return rows


def rust_surface_rows(docs, baseline_rows, dynamic):
    correspondence = defaultdict(list)
    for row in baseline_rows + dynamic:
        for path in row["当前Rust对应"].split("; "):
            if evidence_for(docs, path):
                correspondence[path].append(row["原始命名"])
    rows = []
    for name, item in sorted(docs.items()):
        paths = [(name, item["kind"])]
        paths += [(name + "." + member, "inherent_method") for member in dict.fromkeys(item["inherent_methods"])]
        paths += [(name + "." + field, "public_field") for field in dict.fromkeys(item["fields"])]
        paths += [(name + "." + field, "associated_constant") for field in dict.fromkeys(item["constants"])]
        paths += [(name + "." + field, "enum_variant") for field in dict.fromkeys(item["variants"])]
        for path, kind in paths:
            target = path.replace("_with_options", "_with_params")
            target = target.replace("CipLabelOptions", "CipLabelParams")
            if path == "BioAtomRow.calc_flag":
                target = "BioAtomRow.calculation_flag"
            type_or_metadata = kind in ("struct", "enum", "type", "constant", "public_field", "associated_constant", "enum_variant")
            limitation = "公开可访问值／schema，非算法交付证明" if type_or_metadata else "已有公开实现；未逐项运行旧参数与行为验收"
            if "Binding" in name or name in ("FunctionStatus", "ReceiverCategory", "CallableCategory") or name.isupper():
                limitation += "；框架声明/metadata，不计独立化学算法"
            rows.append({"原始命名": "; ".join(correspondence[path]) or "—（旧清单没有直接对应，可能为新配置/重构组合）",
                         "新版命名": target, "当前Rust对应": path, "条目种类": kind,
                         "Rust接口是否实际可用": "[x] 公开可访问（类型/字段非算法）" if type_or_metadata else STATUS["public"],
                         "Python接口是否实际可用": "[ ] MAIN默认包不可构建",
                         "是否通过5000语料测试": "测试未定义（配置/结果/metadata，固定回归）" if type_or_metadata else "[ ] 本轮未进行对应5000验收",
                         "缺口与限制": limitation, "命名调整": "规范命名建议；当前现名保留在对应列，本次未改代码" if target != path else "保持现名",
                         "公开接口证据": evidence_for(docs, path), "对应旧清单": bool(correspondence[path])})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--evidence", default="target/python-baseline-audit-20261004")
    parser.add_argument("--output-prefix", default="dev/gap_reports/python_0_3_0_current_api_matrix")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    audit_dir = root / args.evidence
    report_dir = root / "dev/gap_reports"
    counts_path = report_dir / "python_baseline_d892ec3_gap_counts.json"
    baseline_path = report_dir / "python_baseline_d892ec3_inventory.tsv"
    counts = json.loads(counts_path.read_text())
    assert counts["baseline"]["commit"] == BASELINE
    with baseline_path.open() as stream:
        baseline = list(csv.DictReader(stream, delimiter="\t"))
    public_baseline = [r for r in baseline if r["baseline_python_api"] != "_rebuild_molecule_from_pickle"]
    assert len(baseline) == 848 and len(public_baseline) == 847
    families = {api: family for family, item in counts["families"].items() for api in item["baseline_apis"]}
    assert len(families) == 551
    owner_families = defaultdict(Counter)
    for api, family in families.items():
        if "." in api:
            owner_families[api.split(".")[0]][family] += 1
    builds = json.loads((audit_dir / "build_results.json").read_text())
    doc_build = json.loads((audit_dir / "rust_doc_result.json").read_text())
    library_tests = json.loads((audit_dir / "rust_library_tests_result.json").read_text())
    assert builds[0]["exit_code"] == 0 and builds[1]["exit_code"] != 0 and doc_build["exit_code"] == 0
    diagnostic_lines = re.findall(r"^error(?:\[E\d+\])?:.*", (root / builds[1]["log"]).read_text(), re.M)
    diagnostic_count = sum(not line.startswith("error: could not compile") for line in diagnostic_lines)
    assert diagnostic_count == 664, "snapshot compiler count changed; review statuses before regeneration"
    assert library_tests["exit_code"] == 0 and any("122 passed; 0 failed; 0 ignored" in s for s in library_tests["summaries"])
    docs = doc_inventory(root, root / "target/doc/cosmolkit")
    (audit_dir / "current_rust_public_surface.json").write_text(json.dumps(docs, ensure_ascii=False, indent=2) + "\n")
    historical = historical_reports(root, audit_dir)
    current_python = json.loads((audit_dir / "current_python_surface.json").read_text())
    source_cache = {}
    rows = []
    for old in public_baseline:
        api, kind = old["baseline_python_api"], old["kinds"]
        family = families.get(api)
        if family is None:
            owner = api.split(".")[0]
            family = owner_families[owner].most_common(1)[0][0] if owner_families[owner] else "support"
        new = canonical(api)
        if kind in ("method", "module_function"):
            state, paths, notes = current_match(root, docs, api, kind, new, family, counts)
        else:
            state, paths, notes = support_match(docs, api, kind, new)
        anchors = [anchor for path in paths if (anchor := evidence_for(docs, path))]
        if state == "public":
            assert paths and len(anchors) == len(paths), (api, paths, anchors)
        if any(name in new for name in ("calc_", "mol_to_", "mol_from_", "get_", "_with_options")):
            raise ValueError("noncanonical proposal: " + new)
        corpus, corpus_evidence = corpus_status(api, kind, family, historical)
        branch, branch_path = branch_note(api)
        if branch_path and not (root / branch_path).resolve().exists():
            branch += "；所引文件当前缺失，仅保留历史通知，不视为本轮源码确认"
        naming = "现存公开名" if evidence_for(docs, new) else "规范映射建议（尚未实现或仅domain）"
        if api != new:
            notes.insert(0, "命名／类型映射；旧别名不是必须新增接口")
        if api in ("get_substruct_matches_with_params", "parse_smarts_with_params"):
            notes.append("short 与配置形态必须保留全部参数；不能只用默认行为替代")
        rows.append({
            "功能组": LABELS[family], "原始命名": api, "新版命名": new,
            "Rust接口是否实际可用": STATUS[state],
            "Python接口是否实际可用": "[ ] 默认包不可构建（本轮664错误）；旧同名声明不计支持",
            "是否通过5000语料测试": corpus,
            "当前Rust对应": "; ".join(paths) or "无已确认对应入口",
            "命名性质": naming, "缺口与限制": "；".join(notes),
            "其他目录补充": branch,
            "公开接口证据": "; ".join(anchors), "5000结果证据": "; ".join(corpus_evidence),
            "其他目录证据": branch_path,
            "基线来源": BASELINE + ":" + old["baseline_source"] + ":" + old["baseline_lines"] + "（历史行号，非当前源码行号）",
            "现存旧绑定声明": old_declaration(root, current_python.get(api), source_cache),
            "现存旧绑定位置": current_python[api]["source"] + ":" + ",".join(map(str, current_python[api]["lines"])) if api in current_python else "未找到",
            "实现/检索证据": "; ".join(implementation_sources(root, paths, family)),
            "条目种类": kind, "Rust状态码": state, "功能组代码": family,
            "当前同名Python声明": "存在" if api in current_python else "不存在",
        })
    ordinary = [r for r in rows if r["条目种类"] in ("method", "module_function")]
    support = [r for r in rows if r not in ordinary]
    assert len(ordinary) == 551 and len(support) == 296
    dynamic = dynamic_surface(root, docs)
    rust_surface = rust_surface_rows(docs, rows, dynamic)
    rust_unmapped = [r for r in rust_surface if not r["对应旧清单"]]
    prefix = root / args.output_prefix
    columns = list(rows[0])
    for suffix, selection in ((".tsv", rows), ("_callables.tsv", ordinary), ("_support.tsv", support)):
        with Path(str(prefix) + suffix).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(selection)
    for suffix, selection in (("_dynamic.tsv", dynamic), ("_rust_surface.tsv", rust_surface), ("_rust_unmapped.tsv", rust_unmapped)):
        with Path(str(prefix) + suffix).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, list(selection[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(selection)
    fingerprints = {}
    for base in (root / "crates", root / "python/src", root / "wasm/src"):
        for path in base.rglob("*.rs"):
            if "tests" not in path.parts:
                fingerprints[str(path.relative_to(root))] = sha(path)
    for path in [root / "Cargo.toml", root / "Cargo.lock", root / "python/Cargo.toml", counts_path, baseline_path,
                 root / "dev/public_api_design.md", root / "dev/crate_architecture.md", root / "dev/tools/python_surface_inventory.py", Path(__file__).resolve()]:
        fingerprints[str(path.relative_to(root))] = sha(path)
    family_counts = {}
    for family in LABELS:
        selected = [r for r in ordinary if r["功能组代码"] == family]
        if selected:
            family_counts[family] = {"label": LABELS[family], "ordinary_callables": len(selected),
                                     "rust_states": dict(Counter(r["Rust状态码"] for r in selected)),
                                     "python_usable": 0,
                                     "historical_5000_rows": sum(bool(r["5000结果证据"]) for r in selected),
                                     "current_5000_accepted": 0}
    audit = {"baseline_commit": BASELINE, "baseline_release": "0.3.0", "audited_utc": datetime.now(timezone.utc).isoformat(),
             "root": str(root), "scope": "MAIN defaults/full；siblings仅补充，不视为合并或当前验收",
             "baseline_provenance": "复用上一ROOT指定提交的inventory/counts；本轮未执行Git，不声称重新取出该提交",
             "baseline_inventory_sha256": sha(baseline_path), "source_sha256": fingerprints,
             "current_python_inventory_sha256": sha(audit_dir / "current_python_surface.json"),
             "rustdoc_index_sha256": sha(root / "target/doc/cosmolkit/index.html"),
             "python_compiler_diagnostics": diagnostic_count,
             "ordinary_callables": 551, "public_static_inventory_rows": 847, "support_rows": 296,
             "distinct_canonical_ordinary_names": len({r["新版命名"] for r in ordinary}),
             "excluded_private_helpers": ["_rebuild_molecule_from_pickle"],
             "builds": builds, "rustdoc": doc_build, "library_tests": library_tests, "families": family_counts,
             "dynamic_supplement_rows": len(dynamic), "public_root_rust_items": len(docs),
             "public_rust_surface_rows": len(rust_surface), "rust_surface_kind_counts": dict(Counter(r["条目种类"] for r in rust_surface)),
             "rust_unmapped_surface_rows": len(rust_unmapped),
             "rust_state_counts": dict(Counter(r["Rust状态码"] for r in ordinary)),
             "python_usable_main": 0, "current_5000_accepted": 0,
             "historical_corpus_records": len(historical),
             "limits": ["公开可编译的对应项不等于完整旧参数/结果/准备状态 parity", "当前没有新的5000验收；旧结果分 scope 保留，不自动打勾",
                        "没有迁移业务代码、改测试、改sole plan、启动agents或执行Git", "已安装旧extension／其他目录选定profile不等于当前MAIN Python可用"]}
    Path(str(prefix) + "_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
    md = ["# Python 0.3.0 → 新版公共功能对照表", "",
          f"基线：`{BASELINE}`（Python 0.3.0）。统计对象为 MAIN `{root}` 的实际工作树；其他工作目录单列。", "",
          "本文件是本轮步骤1–3的盘点证据，不是第二迁移执行计划或验收进度总账。未执行功能迁移、Git或agent恢复。", "",
          "## 统计口径与三个待办列", "",
          "- 主表551项＝433个类普通方法＋118个公开模块函数，包含结果读取方法，不是551个算法。",
          f"- 合并旧模块包装/别名后对应{audit['distinct_canonical_ordinary_names']}个规范目标名称；名称数仍不是算法数。",
          "- 附表296项＝66个类＋13个构造器＋136个属性＋71个协议＋10个类常量；静态公开总表847项。",
          f"- 另外补列{len(dynamic)}项动态表面；不把上一ROOT静态清单未展开的动态槽混入551/847。枚举/MAP按组计，成员与别名写入范围，不当成独立化学算法。",
          "- 原清单848项中的 `_rebuild_molecule_from_pickle` 是私有重建helper，不计普通公开功能；pickle协议本身仍保留。",
          "- Rust `[x]` 仅表示规范对应的公开成员／组合可由当前 `cosmolkit` 访问并编译，行为限于其声明范围；不承诺完整旧参数、边界、错误或parity。部分可用、仅domain、Unsupported与无实现保持 `[ ]`。",
          "- Python按当前默认MAIN包统计，本轮构建664错误，全部 `[ ]`；同名旧声明和已安装旧extension不算新版支持。",
          "- 5000列本轮没有新增验收勾选。已有Rust历史通过／失败写入每行，并核对真实comparison行数、case ID和失败数；缺少当前完整实现身份及Python投影验收，不自动继承。",
          "- 值／参数／结果／协议标“测试未定义”，意思是应制定固定回归或跟随生产者测试，不是已通过或删除覆盖。未实现算法的测试未定义仍属于待办。",
          "- “新版命名”是canonical目标；尚不存在的名称标为映射建议。“当前Rust对应”保留真实现名与组合，二者不能混淆。所有单分子算法移到Molecule语义receiver，移除calc_/mol_to_/mol_from_/冗余get_；原地Molecule尾部_保留。", "",
          "- 历史基线行号附带commit前缀；现存旧绑定签名/行号另列，只证明当前残留声明，不冒充0.3.0完整历史签名。Rust各列是本报告证据分类，不是FunctionStatus，也不修改正式执行ledger。", "",
          "## 分组统计", "",
          "| 功能组 | 旧普通入口 | MAIN公开可编译对应 | 部分可用 | 仅domain | Unsupported | 未找到实现 | Python投影专属 | MAIN Python可用 | 本轮5000验收 |",
          "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for family, item in family_counts.items():
        states = item["rust_states"]
        md.append("| " + " | ".join(map(str, [item["label"], item["ordinary_callables"],
                  *[states.get(s, 0) for s in ("public", "partial", "domain", "unsupported", "absent", "projection")], 0, 0])) + " |")
    md += ["", "## 实际检查", ""]
    for result in builds:
        md.append(f"- `{' '.join(result['command'])}`：exit {result['exit_code']}；日志 `{result['log']}`。")
    md += [f"- `{' '.join(doc_build['command'])}`：exit {doc_build['exit_code']}；用新生成的公开rustdoc检查成员与重导出，不把stale页面或私有方法算公开。",
           f"- `{' '.join(library_tests['command'])}`：exit {library_tests['exit_code']}；122 passed、0 failed、0 ignored（仅cosmolkit lib）。日志 `{library_tests['log']}`。",
           "- 编译/公开面检查与122个library测试均不等于完整workspace/集成/parity通过；本轮未运行5000语料、Python运行测试或完整workspace验收。", "",
           "## 必须保留的具体差异", "",
           "- 描述符78项：22个已有公开Rust对应、34个domain全部或部分primitive、22个未找到新版实现。`calc_num_atoms` 映射 `Molecule.total_atom_count`，不能与图原子数 `Molecule.num_atoms` 混同。",
           "- 当前存在 `molecular_weight_with_options`、`exact_molecular_weight_with_options`、`molecular_formula_with_options`，以及CIP `_with_options`。规范要求参数对象／`_with_params`；本次仅记录，未新增别名或改代码。",
           "- 旧BIO的 `from_pdb/from_mmcif` 接收路径，新版关联短构造接收文本。旧文件语义须映射到 `read_with_format`，不能悄悄改变解释。",
           "- 旧SMARTS解析返回query-bearing Molecule，新架构返回唯一QueryGraph；旧类方法／参数／结果桥必须明确迁移，不能用concrete转换代替。",
           "- 旧read_sdf/read_sdf_from_str均返回单个Molecule；不会改成多记录列表。TopologicalTorsion单分子generator方法与py_score_path迁至Molecule；参数对象改为不可变值，旧共享可变配置不直接继承。",
           "- MAIN mmCIF writer仅coordinate profile，旧33个组开关中的28个缺公共组合输出；私有emitter存在不等于旧组合功能完成。",
           "- MAIN batch/alignment/embedding及forcefield边界仍Unsupported；_1 UFF、_2 Morgan／2D、_4绘图binding不能直接计为MAIN完成。",
           "- _4选定绘图profile：旧300×300默认改为必填，from_smiles缺旧sanitize参数，to_smiles缺旧kwargs，coordinates_2d为list/None而非旧NumPy数组。这些都须保留在对照。",
           "- _2坐标5000/5000是1e-8容差，完整属性状态测试仍红；_1 UFF一步每任务4962成功exact、38错误行仍失败，不能记为5000全通过。", "",
           "## 主表：551个普通可调用入口", ""]
    display = ["原始命名", "新版命名", "Rust接口是否实际可用", "Python接口是否实际可用", "是否通过5000语料测试", "当前Rust对应", "缺口与限制", "其他目录补充"]
    for family in LABELS:
        selected = [r for r in ordinary if r["功能组代码"] == family]
        if selected:
            md += ["### " + LABELS[family], "", markdown_table(selected, display), ""]
    md += ["## 附表：296个静态类型／构造／属性／协议项", "",
           "附表与普通方法分开统计；条目详细源码、公开rustdoc、历史报告路径均保存在同名TSV／audit JSON。", ""]
    for kind in ("class", "constructor", "property", "protocol", "class_attribute"):
        selected = [r for r in support if r["条目种类"] == kind]
        md += ["### " + kind, "", markdown_table(selected, display[:-1]), ""]
    md += ["## 动态公开表面补表", "",
           "前ROOT确认8枚举和4个InChI动态类。这里进一步从现存旧绑定补出MAP、BatchValidationError、字段/构造及模块元数据；这些细项未重新读取指定提交，作为基线待核对补充，不虚增确定的847项。BatchValidationError.errors是额外动态方法候选，不重复计化学算法。", "",
           "Element须保留DUMMY和118种元素（119个值）及Uut/Uup别名；ResidueCode须保留表来源与TRY/WAT/H2O/+A等别名；BondOrder22、BondDirection7、BondStereo8及别名、ChiralTag9、ResidueInfoKind12、BatchErrorMode RAISE=1/KEEP=2的名字/数值/MAP只读语义均需固定回归。", "",
           markdown_table(dynamic, display[:-1]), "",
           "## 新版Rust公共面与旧清单未直接对应项", "",
           f"本轮fresh rustdoc根目录包含{len(docs)}个重导出/直接声明条目，展开inherent方法、公开字段、枚举成员、关联常量共{len(rust_surface)}槽；其中{len(rust_unmapped)}槽未被旧清单映射规则直接引用。不是新增{len(rust_unmapped)}个化学算法：包含模型/参数/错误/schema/框架metadata，也可能是旧行为拆分的新接口。", "",
           "完整可检索清单在 `_rust_surface.tsv`，未直接对应项在 `_rust_unmapped.tsv`；不重复展开Clone/Debug等自动trait和model模块别名。下表列Molecule、BIO及module函数的未直接对应项，其他类型详见TSV。", "",
           markdown_table([r for r in rust_unmapped if r['当前Rust对应'].split('.')[0] in ('Molecule', 'BioStructure', 'Protein') or r['条目种类'] == 'fn'], display[:-1]), "",
           "## 额外需求的边界", "",
           "SMIRKS是高速plan新增要求，不属于0.3.0旧公开功能差额；本次未实现。BIO新增PDB输出、CID selection、邻居／链序列与tautomer也不计旧基线缺口。", "",
           "## 证据与复现", "",
           f"- 基线清单：`{baseline_path.relative_to(root)}`；本轮复用指定提交统计，未用Git重新读取提交。",
           f"- 全部847行：`{prefix.name}.tsv`；普通551行：`{prefix.name}_callables.tsv`；支持296行：`{prefix.name}_support.tsv`。",
           f"- 动态补表：`{prefix.name}_dynamic.tsv`；Rust公共面：`{prefix.name}_rust_surface.tsv`；未直接对应项：`{prefix.name}_rust_unmapped.tsv`。",
           f"- 元数据、源码SHA256和命令退出码：`{prefix.name}_audit.json`。",
           f"- 本轮公开Rust目录证据：`{args.evidence}/current_rust_public_surface.json`。",
           f"- 核对后的历史5000结果：`{args.evidence}/historical_corpus_evidence.json`，生成数据仍留target，不复制大语料或reference。",
           f"- 重建表格：`python3 dev/tools/python_baseline_gap_table.py --root . --evidence {args.evidence}`；必须先提供本轮build／doc／current Python inventory材料。"]
    Path(str(prefix) + ".md").write_text("\n".join(md) + "\n")
    print(json.dumps({"ordinary": len(ordinary), "support": len(support), "rows": len(rows),
                      "states": audit["rust_state_counts"], "historical_records": len(historical),
                      "outputs": str(prefix)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
