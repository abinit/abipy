"""
Tools used to parallelize sections of python code.
"""
from __future__ import annotations

import os

from monty.collections import dict2namedtuple


_MAX_NPROCS = os.cpu_count()


def get_max_nprocs() -> int:
    """
    Return the maximum number of procs that can be used by AbiPy.
    """
    return _MAX_NPROCS


def set_max_nprocs(max_nprocs: int | None) -> int:
    """
    Set the maximum number of procs that can be used.
    If max_nprocs is None, os.cpu_count() is used.
    """
    global _MAX_NPROCS
    if max_nprocs is None:
        _MAX_NPROCS = os.cpu_count()
    else:
        _MAX_NPROCS = max_nprocs

    return _MAX_NPROCS


def pool_nprocs_pmode(nprocs: int | None, pmode: str):
    """
    Helper function that allows one to switch from ThreadPool to multiprocessing Pool.
    Returns named tuple with (nprocs, pool_class, using_msg)

    Args:
        nprocs: Number of procs to use. If None, use cpu_count.
        pmode: "threads": for ThreadPool.
               "processes" for multiprocessing Pool.
               "seq" for sequential execution (debugging)
    """
    from multiprocessing.pool import ThreadPool
    from multiprocessing import Pool

    max_nprocs = get_max_nprocs()

    if pmode == "seq":
        nprocs, pool_cls = 1, ThreadPool

    elif pmode == "threads":
        pool_cls = ThreadPool
        if nprocs is not None:
            pool_cls = ThreadPool if nprocs > 0 else Pool
            nprocs = min(abs(nprocs), max_nprocs)

    elif pmode == "processes":
        pool_cls = Pool
        if nprocs is not None:
            pool_cls = Pool if nprocs > 0 else ThreadPool
            nprocs = min(abs(nprocs), max_nprocs)

    else:
        raise ValueError(f"Invalid value of {pmode=}, it should be in ['seq', 'threads', 'processes']")

    return dict2namedtuple(nprocs=nprocs or max_nprocs,
                           pool_cls=pool_cls,
                           using_msg=f"using {nprocs=} with {pmode=} and Pool class: {pool_cls.__name__} ...",
                           )
