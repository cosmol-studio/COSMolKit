"""User-selected 60-second deadline for common UFF/MMFF reference geometry.

The original geometry recipe stays in fingerprint_values_pilot. A supervised
process makes the deadline effective even inside native loops that do not poll
RDKit's cooperative timeout. Each existing reference worker reuses one child.
"""
from __future__ import annotations

import atexit
import multiprocessing
import time

TIMEOUT_SECONDS = 60


class GeometryWorkerError(Exception):
    """Reference infrastructure failure; never a chemistry error or timeout."""


def _geometry(row, first_case_id, label):
    from fingerprint_values_pilot import prepare_forcefield_geometry
    return prepare_forcefield_geometry(
        row, first_case_id=first_case_id, label=label,
        timeout_seconds=TIMEOUT_SECONDS)


def _serve(connection, operation):
    try:
        while True:
            try:
                arguments = connection.recv()
            except EOFError:
                return
            try:
                result = operation(*arguments)
            except Exception as error:
                connection.send((False, error))
            else:
                connection.send((True, result))
    finally:
        connection.close()


class GeometrySupervisor:
    def __init__(self, operation=_geometry, timeout_seconds=TIMEOUT_SECONDS):
        if timeout_seconds <= 0:
            raise ValueError("deadline must be positive")
        self.operation = operation
        self.timeout_seconds = timeout_seconds
        self.process = None
        self.connection = None

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None
        if self.process is not None:
            if self.process.is_alive():
                self.process.terminate()
            self.process.join(1)
            if self.process.is_alive():
                self.process.kill()
                self.process.join()
            self.process.close()
            self.process = None

    def prepare(self, *arguments):
        deadline = time.monotonic() + self.timeout_seconds
        if self.process is None:
            context = multiprocessing.get_context("spawn")
            self.connection, child = context.Pipe()
            self.process = context.Process(
                target=_serve, args=(child, self.operation), daemon=True)
            try:
                self.process.start()
            except BaseException as error:
                self.connection.close()
                self.connection = None
                self.process.close()
                self.process = None
                if isinstance(error, Exception):
                    raise GeometryWorkerError("could not start geometry worker") from error
                raise
            finally:
                child.close()
        try:
            self.connection.send(arguments)
            if not self.connection.poll(max(0, deadline - time.monotonic())):
                self.close()
                return {"TimedOut": {"limit_seconds": self.timeout_seconds,
                                     "mechanism": "ProcessDeadline"}}
            ok, value = self.connection.recv()
            if not ok:
                if isinstance(value, TimeoutError):
                    return {"TimedOut": {"limit_seconds": self.timeout_seconds,
                                         "mechanism": "Native"}}
                raise value
            return value
        except (EOFError, BrokenPipeError, ConnectionResetError) as error:
            self.close()
            raise GeometryWorkerError("geometry worker exited without a result") from error


_supervisor = None
_timeouts = {}


def prepare_geometry(row, first_case_id=None, label="UFF"):
    global _supervisor
    name, options = next(iter(row["profile"].items()))
    if name not in ("Optimization", "ConformerOptimization"):
        return row["preparation"]
    if row["preparation"] is not None:
        return row["preparation"]
    key = (row["case"]["smiles"], options["add_hydrogens"])
    if key in _timeouts:
        row["preparation"] = _timeouts[key]
        return row["preparation"]
    if _supervisor is None:
        _supervisor = GeometrySupervisor()
        atexit.register(_supervisor.close)
    row["preparation"] = _supervisor.prepare(row, first_case_id, label)
    if "TimedOut" in row["preparation"]:
        _timeouts[key] = row["preparation"]
    return row["preparation"]
