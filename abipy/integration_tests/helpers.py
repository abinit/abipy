"""Helper functions shared by the integration tests."""

from __future__ import annotations

import io
import os
import warnings

import pytest

# Files in the task workdir that are useful to understand why a task failed.
_ERROR_FILES = ("run.err", "queue.qerr", "__ABI_MPIABORTFILE__", "run.log", "run.abo")


def assert_flow_ok(flow, scheduler=None, check_finalized=True) -> None:
    """
    Check that all the tasks in the flow completed successfully.

    If this is not the case, call pytest.fail with a report containing the status of the flow,
    the output of flow.debug and the exceptions raised by the scheduler (if any).
    Exceptions recorded by the scheduler are only reported as warnings if the flow completed successfully.

    Args:
        flow: Flow to check.
        scheduler: Scheduler used to run the flow (optional).
        check_finalized: True if all the works in the flow should be finalized.
    """
    excs = list(scheduler.exceptions) if scheduler is not None else []
    flow.check_status()

    if not flow.all_ok:
        buf = io.StringIO()
        flow.show_status(stream=buf)
        flow.debug(stream=buf)
        if excs:
            buf.write("\nExceptions raised by the scheduler:\n" + "\n".join(str(e) for e in excs))
        pytest.fail(f"Flow in {flow.workdir} did not complete successfully:\n{buf.getvalue()}", pytrace=False)

    if excs:
        warnings.warn(
            f"Flow in {flow.workdir} completed but the scheduler recorded {len(excs)} exception(s):\n"
            + "\n".join(str(e) for e in excs),
            stacklevel=2,
        )

    if check_finalized:
        not_finalized = [str(work) for work in flow if not work.finalized]
        if not_finalized:
            pytest.fail(f"Flow is ok but these works are not finalized: {not_finalized}", pytrace=False)


def run_flow(flow, scheduler=None, check_finalized=True):
    """
    Run the flow with the scheduler and check that all the tasks completed successfully.
    Use flow.make_scheduler() if scheduler is None. Return the scheduler.
    """
    if scheduler is None:
        scheduler = flow.make_scheduler()
    scheduler.start()
    assert_flow_ok(flow, scheduler=scheduler, check_finalized=check_finalized)
    return scheduler


def tail(filepath: str, nlines: int = 40) -> str:
    """Return the last `nlines` lines of a text file."""
    with open(filepath, errors="replace") as fh:
        return "".join(fh.readlines()[-nlines:])


def collect_task_errors(workdir: str, max_tasks: int = 5, nlines: int = 40) -> str:
    """
    Scan the task directories inside `workdir` and return a string with the tail
    of the files that signal an error (non-empty stderr, MPI abort file, ERROR/BUG in the log).
    """
    sections, ntasks = [], 0
    for dirpath, _, filenames in sorted(os.walk(workdir)):
        # Each task directory contains the input file run.abi.
        if "run.abi" not in filenames:
            continue

        paths = {f: os.path.join(dirpath, f) for f in _ERROR_FILES if f in filenames}
        has_error = any(
            os.path.getsize(paths[f]) > 0 for f in ("run.err", "queue.qerr", "__ABI_MPIABORTFILE__") if f in paths
        )
        if not has_error and "run.log" in paths:
            with open(paths["run.log"], errors="replace") as fh:
                log = fh.read()
            has_error = "--- !ERROR" in log or "--- !BUG" in log

        if not has_error:
            continue

        for path in paths.values():
            if os.path.getsize(path) > 0:
                sections.append(f"==> {path} <==\n{tail(path, nlines=nlines)}")

        ntasks += 1
        if ntasks >= max_tasks:
            sections.append(f"Output truncated after {max_tasks} tasks.")
            break

    return "\n".join(sections)
