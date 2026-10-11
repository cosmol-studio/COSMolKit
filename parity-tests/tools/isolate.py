"""Isolate crash-prone native calls, bisecting only a genuinely failed job."""
import multiprocessing
import struct


def _child(connection, evaluate, task, profiles, rows):
    # A native crash must not create a multi-gigabyte core dump.
    try:
        import resource
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    except ImportError:
        pass
    try:
        connection.send((True, evaluate(task, profiles, rows)))
    except BaseException as error:
        connection.send((False, f"{type(error).__name__}: {error}"))
    finally:
        connection.close()


def isolated(evaluate, task, profiles, rows, timeout=120):
    """Return one encoded outcome per cell; infrastructure exceptions still fail."""
    context = multiprocessing.get_context("fork" if "fork" in multiprocessing.get_all_start_methods() else "spawn")
    receiver, sender = context.Pipe(duplex=False)
    child = context.Process(target=_child, args=(sender, evaluate, task, profiles, rows))
    child.start()
    sender.close()
    timed_out = False
    result = None
    try:
        if receiver.poll(timeout):
            try:
                result = receiver.recv()
            except EOFError:
                pass
        else:
            timed_out = True
            child.terminate()
        child.join(timeout=5)
        if child.is_alive():
            child.kill()
            child.join()
            timed_out = True
    finally:
        receiver.close()
    if result is not None:
        ok, value = result
        if not ok:
            raise RuntimeError(value)
        if child.exitcode != 0:
            raise RuntimeError(f"native worker exited {child.exitcode} after returning results")
        return value
    if not timed_out and (child.exitcode is None or child.exitcode >= 0):
        raise RuntimeError(f"native worker exited {child.exitcode} without returning results")
    if len(rows) > 1:
        middle = len(rows) // 2
        return isolated(evaluate, task, profiles, rows[:middle], timeout) + isolated(evaluate, task, profiles, rows[middle:], timeout)
    if len(profiles) > 1:
        middle = len(profiles) // 2
        left = isolated(evaluate, task, profiles[:middle], rows, timeout)
        right = isolated(evaluate, task, profiles[middle:], rows, timeout)
        return [left[0] + right[0]]
    status = 3 if timed_out else 4
    message = (f"native timeout after {timeout}s" if timed_out else f"native crash signal {-child.exitcode}").encode()
    return [[bytes((status, 5)) + struct.pack("<I", len(message)) + message]]
