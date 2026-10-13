import hashlib
import struct
import unittest
import zlib
import tempfile
from pathlib import Path
from unittest.mock import patch
import prepare


class FrameTests(unittest.TestCase):
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
