"""Source-only pinned Torsions.py supplement; ordinary tests never run RDKit."""
import argparse
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.AtomPairs import Torsions

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert hashlib.sha256(args.corpus.read_bytes()).hexdigest() == "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
    assert hashlib.sha256(Path(Torsions.__file__).read_bytes()).hexdigest() == "27aa284d8ba4e51cb0059694bbe655e8309e860b4cb843a203a9bd9769f6b2d2"
    rows = [s.strip() for s in args.corpus.read_text().splitlines() if s.strip() and not s.strip().startswith("#")]
    assert len(rows) == 5000
    totals = {"inputs": 0, "profiles": 0, "score_calls": 0, "explanation_calls": 0, "errors": 0, "test_count": 0}
    with args.output.open("x") as output:
        for index, text in enumerate(rows):
            molecule = Chem.MolFromSmiles(text)
            assert molecule is not None, (index, text, "source parse returned null")
            atoms = molecule.GetNumAtoms()
            assert atoms > 0, (index, text, "source corpus empty molecule")
            path = list(range(min(4, atoms)))
            profiles = [
                ("first_default", path, len(path), None),
                ("reverse_empty_default", path[::-1], len(path), []),
                ("repeated_default", [i % atoms for i in range(4)], 4, None),
                ("first_custom_index_plus11", path, len(path), [i + 11 for i in range(atoms)]),
            ]
            record = {"index": index, "smiles": text, "atoms": atoms, "profiles": {}}
            for name, selected, size, codes in profiles:
                value = {"path": selected, "size": size, "atom_codes": codes}
                totals["score_calls"] += 1
                try:
                    score = Torsions.pyScorePath(molecule, selected, size, codes)
                    assert 0 <= score < (1 << 64)
                    value["score"] = score
                    totals["explanation_calls"] += 1
                    value["explanation"] = Torsions.ExplainPathScore(score, size)
                    value["error"] = None
                except Exception as error:
                    value["error"] = {"type": type(error).__name__, "message": str(error)}
                    totals["errors"] += 1
                record["profiles"][name] = value
                totals["profiles"] += 1
            output.write(json.dumps(record, sort_keys=True) + "\n")
            totals["inputs"] += 1
    print(json.dumps(totals), flush=True)

if __name__ == "__main__":
    main()
