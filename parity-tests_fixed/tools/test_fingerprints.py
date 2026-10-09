"""Reference-process boundaries; no corpus or expected writes."""
import copy
import os
import pickle
import unittest
from unittest.mock import MagicMock, patch

import fingerprints


def recipe():
    return {"Fingerprint": {"case": {"id": "fixed-CCO", "smiles": "CCO"},
                            "right": None, "seed": 0, "params": {"Layered": {
                                "layers": 63, "min_path": 1, "max_path": 7,
                                "fp_size": 2048, "branched": False, "roots": "All",
                                "counts": "Seeded", "mask": "Absent"}}}}


def process(exit_code, output=b"", diagnostics=b""):
    child = MagicMock()
    child.__enter__.return_value = child
    child.returncode = exit_code
    child.pid = 123
    child.communicate.return_value = (output, diagnostics)
    return child


class IsolationTests(unittest.TestCase):
    def test_success_preserves_complete_record_and_original_input(self):
        row = recipe()
        before = copy.deepcopy(row)
        result = {"input": row, "output": {"Fingerprint": {"Bits": {
            "fingerprint": {"length": 2048, "on_bits": [3, 8]}, "atom_counts": [12, 14, 12]}}}}
        child = process(0, pickle.dumps((True, result)))
        with patch.object(fingerprints.subprocess, "Popen", return_value=child) as launch:
            self.assertEqual(fingerprints.fingerprint_case(row), result)
        self.assertEqual(pickle.loads(child.communicate.call_args.args[0]), before)
        self.assertEqual(row, before)
        self.assertEqual(launch.call_args.kwargs["env"]["OPENBLAS_NUM_THREADS"], "1")

    def test_native_process_exits_keep_input_and_diagnostics(self):
        for code in [-11, -9, 0xC0000005]:
            with self.subTest(code=code):
                row = recipe()
                child = process(code, diagnostics=b"native diagnostic")
                with patch.object(fingerprints.subprocess, "Popen", return_value=child):
                    observed = fingerprints.fingerprint_case(row)
                self.assertEqual(observed["input"], row)
                self.assertEqual(observed["output"]["Fingerprint"]["ReferenceProcessFailure"],
                                 {"exit_code": code, "process_id": 123,
                                  "stderr": "native diagnostic"})

    def test_ordinary_source_exception_keeps_its_category(self):
        child = process(0, pickle.dumps((False, ValueError("source argument error"))))
        with patch.object(fingerprints.subprocess, "Popen", return_value=child):
            with self.assertRaisesRegex(ValueError, "source argument error"):
                fingerprints.fingerprint_case(recipe())

    def test_transport_exit_is_not_a_native_failure_observation(self):
        child = process(1, diagnostics=b"transport failure")
        with patch.object(fingerprints.subprocess, "Popen", return_value=child):
            with self.assertRaisesRegex(RuntimeError, "transport exited 1"):
                fingerprints.fingerprint_case(recipe())

    def test_spawn_error_is_not_a_native_failure_observation(self):
        with patch.object(fingerprints.subprocess, "Popen", side_effect=OSError("launch failed")):
            with self.assertRaisesRegex(OSError, "launch failed"):
                fingerprints.fingerprint_case(recipe())

    def test_safe_profiles_keep_direct_execution(self):
        for roots, branched in [("First", False), ("All", True)]:
            row = recipe()
            row["Fingerprint"]["params"]["Layered"].update(roots=roots, branched=branched)
            result = {"unchanged": row}
            with patch.object(fingerprints, "_fingerprint_case", return_value=result) as native:
                with patch.object(fingerprints.subprocess, "Popen") as launch:
                    self.assertIs(fingerprints.fingerprint_case(row), result)
            native.assert_called_once_with(row)
            launch.assert_not_called()

    @unittest.skipUnless(os.name == "posix", "fixed native signal assertion is POSIX-specific")
    def test_fixed_target_native_crash_is_isolated_and_retained(self):
        from rdkit import rdBase
        self.assertEqual(rdBase.rdkitVersion, "2026.03.6")
        row = recipe()
        observed = fingerprints.fingerprint_case(row)
        self.assertEqual(observed["input"], row)
        failure = observed["output"]["Fingerprint"]["ReferenceProcessFailure"]
        self.assertEqual(failure["exit_code"], -11)
        self.assertGreater(failure["process_id"], 0)
        self.assertIn("Segmentation fault", failure["stderr"])


if __name__ == "__main__":
    unittest.main()
