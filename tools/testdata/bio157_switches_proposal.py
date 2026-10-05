"""Generate BIO157 native switch proposals; independent review remains required.

Compile bio157_switch_oracle_proposal.cpp against the pinned Gemmi source and
libgemmi_cpp, then pass that executable explicitly. This uses the existing
profile unchanged and writes fresh external evidence rather than repo data.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

FLAGS = ['atoms', 'block_name', 'entry', 'database_status', 'author', 'cell', 'symmetry', 'entity', 'entity_poly', 'struct_ref', 'chem_comp', 'exptl', 'diffrn', 'reflns', 'refine', 'title_keywords', 'ncs', 'struct_asym', 'origx', 'struct_conf', 'struct_sheet', 'struct_biol', 'assembly', 'conn', 'cis', 'modres', 'scale', 'atom_type', 'entity_poly_seq', 'tls', 'software', 'group_pdb', 'auth_all']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--oracle", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    out = args.output_directory.resolve()
    out.mkdir(parents=True, exist_ok=False)
    oracle = args.oracle.resolve(strict=True)
    profile = json.loads((root / "testdata/bio/gemmi_mmcif_writer_profile.json").read_text())
    rows = []
    for case in profile["cases"]:
        for all_groups in (False, True):
            for flag in FLAGS:
                value = not all_groups if flag != "auth_all" else True
                filename = f"{case['case_id']}-{int(all_groups)}-{flag}.cif"
                destination = out / filename
                argv = [str(oracle), str(root / case["input"]), str(destination),
                        str(int(all_groups)), flag, str(int(value))]
                subprocess.run(argv, check=True)
                data = destination.read_bytes()
                destination.chmod(0o444)
                rows.append(dict(case_id=case["case_id"], input=case["input"],
                                 all_groups=all_groups, flag=flag, value=value,
                                 argv=argv, exit=0, output=filename,
                                 bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
    manifest = out / "cases-proposal.json"
    manifest.write_text(json.dumps(rows, indent=2))
    manifest.chmod(0o444)
    print(json.dumps({"native_cases":len(rows), "manifest":str(manifest)}))

if __name__ == "__main__":
    main()
