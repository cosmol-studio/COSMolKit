import os
import signal
import time
import unittest
from isolate import isolated


def cells(task, profiles, rows):
    output = []
    for row in rows:
        values = []
        for profile in profiles:
            if row == "bad" and profile:
                if task == "crash":
                    os.kill(os.getpid(), signal.SIGSEGV)
                time.sleep(10)
            values.append(b"\0\0")
        output.append(values)
    return output


class IsolationTests(unittest.TestCase):
    def test_crash_is_not_a_pass_and_only_failed_cell_is_classified(self):
        values = isolated(cells, "crash", [False, True], ["good", "bad"])
        self.assertEqual(values[0], [b"\0\0", b"\0\0"])
        self.assertEqual(values[1][0], b"\0\0")
        self.assertEqual(values[1][1][0], 4)

    def test_timeout_remains_distinct_from_crash(self):
        values = isolated(cells, "timeout", [True], ["bad"], timeout=0.1)
        self.assertEqual(values[0][0][0], 3)


if __name__ == "__main__":
    unittest.main()
