"""Focused preparation-tool regressions; no corpus files are generated."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
import unittest

import reference

sys.path.insert(0, str(reference.ROOT / "tools/oracles/rdkit"))
import fingerprint_values_pilot as oracle


def twice(value):
    return 2 * value


class ProgressTests(unittest.TestCase):
    def test_mmcif_reference_retains_shared_bio_input_filename(self):
        import gemmi
        for name, format in (("sample.pdb", "pdb"), ("sample.cif", "mmcif")):
            path = reference.PACKAGE / "testdata/bio" / name
            case = {"case_id": name, "input": f"bio/{name}", "format": format,
                    "text": path.read_text()}
            groups = gemmi.MmcifOutputGroups(False)
            groups.block_name = True
            expected = gemmi.read_structure(str(path), merge_chain_parts=False)
            row = reference.bio_mmcif_switch_case((case, False, "block_name"))
            self.assertEqual(row["text"], expected.make_mmcif_document(groups).as_string())
            self.assertEqual(row["flag"], "block_name")
            self.assertTrue(row["value"])

    def test_progress_is_flushed_to_stderr_not_json_stdout(self):
        stderr, stdout = io.StringIO(), io.StringIO()
        with contextlib.redirect_stderr(stderr), contextlib.redirect_stdout(stdout):
            progress = reference.Progress("smiles_read_smiles")
            progress(0, 3)
            progress(3, 3)
        self.assertEqual(stdout.getvalue(), "")
        self.assertIn("smiles_read_smiles [", stderr.getvalue())
        self.assertIn("0/3 cases", stderr.getvalue())
        self.assertIn("3/3 cases (100.0%)", stderr.getvalue())

    def test_parallel_completion_preserves_input_order(self):
        from threading import Event
        later_finished = Event()

        def work(value):
            if value == 0:
                self.assertTrue(later_finished.wait(5))
            else:
                later_finished.set()
            return twice(value)

        updates = []
        with ThreadPoolExecutor(max_workers=3) as pool:
            values = reference.collect(pool, work, [0, 1, 2], lambda *u: updates.append(u))
        self.assertEqual(values, [0, 2, 4])
        self.assertEqual(updates, [(0, 3), (1, 3), (2, 3), (3, 3)])

    def test_failed_worker_never_reports_full_completion(self):
        updates = []

        def fail(_):
            raise ValueError("worker failed")

        with ThreadPoolExecutor(max_workers=1) as pool:
            with self.assertRaisesRegex(ValueError, "worker failed"):
                reference.collect(pool, fail, [0], lambda *u: updates.append(u))
        self.assertEqual(updates, [(0, 1)])

    def test_rdkit_progress_does_not_change_rows_or_order(self):
        cases = [{"id": str(i), "smiles": s} for i, s in enumerate(("CCO", "CC", "C"))]
        parameters = [{"SmilesRead": {"sanitize": True, "remove_hydrogens": True}}]
        expected = oracle.generate_smiles_read(cases, parameters, 1)
        for threads in (1, 2):
            updates = []
            actual = oracle.generate_smiles_read(
                cases, parameters, threads, progress=lambda *u: updates.append(u))
            self.assertEqual(actual, expected)
            self.assertEqual(updates, [(0, 3), (1, 3), (2, 3), (3, 3)])

    def test_real_generator_keeps_json_stdout_and_visible_task_progress(self):
        cases = [{"id": "ethanol", "smiles": "CCO"}]
        parameters = [{"SmilesRead": {"sanitize": True, "remove_hydrogens": True}}]
        result = subprocess.run(
            [sys.executable, str(Path(reference.__file__))],
            input=json.dumps({"kind": "corpus", "task": "smiles_read_smiles",
                             "generator": "generate_smiles_read", "corpus": cases,
                             "parameters": parameters, "threads": 2}),
            text=True, capture_output=True, check=True)
        self.assertEqual(json.loads(result.stdout), oracle.generate_smiles_read(cases, parameters, 1))
        self.assertIn("smiles_read_smiles [", result.stderr)
        self.assertIn("1/1 cases (100.0%)", result.stderr)


if __name__ == "__main__":
    unittest.main()
