import hashlib
import struct
import unittest
import zlib
import tempfile
from pathlib import Path
from unittest.mock import patch
import prepare


class FrameTests(unittest.TestCase):
    def test_inchi_rejection_requires_confirmed_engine_status(self):
        for code in (-1, 0, 1, 2, 3, 4, 5):
            with self.subTest(code=code), patch.object(prepare.inchi, "MolToInchi", return_value=""), patch.object(prepare.rdinchi, "MolToInchi", return_value=("", code, "diagnostic", "", "")):
                rows = prepare.evaluate_rows("inchi_write", [{"options": ""}], ["C"])
                expected = b"\1" + prepare.encode(["InchiGenerationRejected", code])
                if code in (2, 3):
                    self.assertEqual(rows, [[expected]])
                else:
                    self.assertEqual(rows[0][0][0], 5)
        with patch.object(prepare.inchi, "MolToInchi", return_value="InChI=1S/CH4/h1H4"), patch.object(prepare.rdinchi, "MolToInchi") as raw:
            self.assertEqual(prepare.inchi_identifier(None), "InChI=1S/CH4/h1H4")
            raw.assert_not_called()

    def test_record_above_old_ceiling_roundtrips_without_truncation(self):
        value = b"x" * (64 * 1024 * 1024)
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(prepare, "iter_outcomes", return_value=iter([[value]])):
                count, path = prepare.calculate(("test", [{}], ["C"], folder))
            data = Path(path).read_bytes()
        self.assertEqual(count, 1)
        stored, expanded, records = struct.unpack_from("<III", data)
        self.assertEqual(records, 1)
        self.assertEqual(len(data), 44 + stored)
        self.assertEqual(hashlib.sha256(data[44:]).digest(), data[12:44])
        raw = zlib.decompress(data[44:])
        self.assertEqual(len(raw), expanded)
        self.assertEqual(struct.unpack_from("<III", raw), (8 + len(value), 1, len(value)))
        self.assertEqual(raw[12:], value)

    def test_record_above_ceiling_is_still_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(prepare, "iter_outcomes", return_value=iter([[b"x" * 60]])), patch.object(prepare, "MAX_RECORD", 64):
                with self.assertRaisesRegex(RuntimeError, "reference row.*exceeds"):
                    prepare.calculate(("test", [{}], ["C"], folder))

    def test_large_job_splits_frames_without_losing_input_or_profile(self):
        outputs = [[b"\0\0"] * 3 for _ in range(5)]
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(prepare, "iter_outcomes", return_value=iter(outputs)), patch.object(prepare, "MAX_FRAME", 64):
                count, path = prepare.calculate(("test", [{}, {}, {}], ["C"] * 5, folder))
            data = Path(path).read_bytes()
        self.assertEqual(count, 3)
        cursor, records = 0, 0
        while cursor < len(data):
            stored, expanded, size = struct.unpack_from("<III", data, cursor)
            checksum = data[cursor + 12:cursor + 44]
            payload = data[cursor + 44:cursor + 44 + stored]
            self.assertEqual(hashlib.sha256(payload).digest(), checksum)
            self.assertEqual(len(zlib.decompress(payload)), expanded)
            records += size
            cursor += 44 + stored
        self.assertEqual(records, 5)


if __name__ == "__main__":
    unittest.main()
