"""Deadline boundary tests; no RDKit calls or corpus generation."""
import multiprocessing
import os
import time
import unittest

from forcefield_preparation import GeometrySupervisor, GeometryWorkerError


def operation(action):
    if action == "block":
        time.sleep(30)
    if action == "native_timeout":
        raise TimeoutError("native exceeded timeout")
    if action == "ordinary_error":
        raise ValueError("ordinary source rejection")
    if action == "crash":
        os._exit(23)
    return {"Ready": {"action": action, "pid": os.getpid()}}


class GeometryDeadlineTests(unittest.TestCase):
    def setUp(self):
        self.supervisor = GeometrySupervisor(operation, timeout_seconds=5)
        self.addCleanup(self.supervisor.close)

    def test_success_reuses_child_and_preserves_result(self):
        first = self.supervisor.prepare("one")
        second = self.supervisor.prepare("two")
        self.assertEqual(first["Ready"]["action"], "one")
        self.assertEqual(second["Ready"]["action"], "two")
        self.assertEqual(first["Ready"]["pid"], second["Ready"]["pid"])

    def test_hard_deadline_reaps_worker_then_next_case_runs(self):
        pid = self.supervisor.prepare("warmup")["Ready"]["pid"]
        self.supervisor.timeout_seconds = 0.1
        started = time.monotonic()
        timed = self.supervisor.prepare("block")
        self.assertEqual(timed, {"TimedOut": {"limit_seconds": 0.1, "mechanism": "ProcessDeadline"}})
        self.assertLess(time.monotonic() - started, 3)
        self.assertIsNone(self.supervisor.process)
        self.assertNotIn(pid, [child.pid for child in multiprocessing.active_children()])
        self.supervisor.timeout_seconds = 5
        self.assertNotEqual(self.supervisor.prepare("next")["Ready"]["pid"], pid)

    def test_native_timeout_has_its_own_mechanism(self):
        self.assertEqual(self.supervisor.prepare("native_timeout"),
                         {"TimedOut": {"limit_seconds": 5, "mechanism": "Native"}})
        self.assertEqual(self.supervisor.prepare("after")["Ready"]["action"], "after")

    def test_ordinary_exception_is_not_timeout(self):
        with self.assertRaisesRegex(ValueError, "ordinary source rejection"):
            self.supervisor.prepare("ordinary_error")
        self.assertEqual(self.supervisor.prepare("after")["Ready"]["action"], "after")

    def test_worker_crash_is_not_timeout(self):
        with self.assertRaisesRegex(GeometryWorkerError, "worker exited without a result"):
            self.supervisor.prepare("crash")
        self.assertIsNone(self.supervisor.process)

    def test_invalid_deadline_is_rejected(self):
        with self.assertRaises(ValueError):
            GeometrySupervisor(operation, timeout_seconds=0)

    def test_worker_start_failure_is_not_a_source_error_or_timeout(self):
        broken = GeometrySupervisor(lambda: None)
        self.addCleanup(broken.close)
        with self.assertRaisesRegex(GeometryWorkerError, "could not start geometry worker"):
            broken.prepare()
        self.assertIsNone(broken.process)
        self.assertIsNone(broken.connection)


if __name__ == "__main__":
    unittest.main()
