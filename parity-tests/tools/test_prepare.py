import hashlib
import struct
import unittest
import zlib
import tempfile
from pathlib import Path
from unittest.mock import patch
import prepare


class FrameTests(unittest.TestCase):
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
