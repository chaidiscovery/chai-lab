# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""
Timeout utility for a function, creates a new process

Implementation modified from:
https://www.reddit.com/r/Python/comments/8t9bk4/the_absolutely_easiest_way_to_time_out_a_function/
"""

import multiprocessing
import queue as _queue
from enum import Enum
from functools import wraps
from multiprocessing import Process, Queue
from typing import Any

from typing_extensions import assert_never


# TODO: This is dangerous: revert once the underlying problem in rdkit is fixed
# RDKit Issue(https://github.com/rdkit/rdkit/discussions/7289)
class Undaemonize(object):
    """Context Manager to resolve AssertionError: daemonic processes are not allowed to have children
    See https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic"""

    def __init__(self):
        self.conf: dict = multiprocessing.process.current_process()._config  # type: ignore
        if "daemon" in self.conf:
            self.daemon_status_set = True
        else:
            self.daemon_status_set = False
        self.daemon_status_value = self.conf.get("daemon")

    def __enter__(self):
        if self.daemon_status_set:
            del self.conf["daemon"]

    def __exit__(self, type, value, traceback):
        if self.daemon_status_set:
            self.conf["daemon"] = self.daemon_status_value


class HandlerStatus(Enum):
    SUCCESS = 0
    EXCEPTION = 1


class ChildProcessException(Exception):
    pass


def timeout(timeout: float | int) -> Any:
    """Force function to timeout after 'seconds'.

    Returns:
        The return value of the wrapped function.
    Raises:
        TimeoutError if the function does not return before the timeout.
    """

    def handler(queue, func, args, kwargs) -> None:
        try:
            queue.put((HandlerStatus.SUCCESS, func(*args, **kwargs)))
        except Exception as e:
            queue.put((HandlerStatus.EXCEPTION, e))

    def decorator(func):
        @wraps(func)
        def new_fn(*args, **kwargs):
            queue: Queue = Queue()
            proc = Process(
                target=handler, args=(queue, func, args, kwargs), daemon=True
            )
            with Undaemonize():
                proc.start()
            proc.join(timeout=float(timeout))
            if proc.is_alive():
                proc.terminate()
                proc.join()
                raise TimeoutError(f"Function {func} timed out after {timeout} seconds")
            else:
                # When child process dies unexpectedly Queue.get waits indefinitely.
                # See Issue(https://bugs.python.org/issue43805)
                # prevent queue from hanging with another very short timeout
                try:
                    status, value = queue.get(timeout=0.1)
                except _queue.Empty:
                    # in this case, child process has died unexpectedly
                    raise ChildProcessException("Child process died unexpectedly")

                match status:
                    case HandlerStatus.SUCCESS:
                        return value
                    case HandlerStatus.EXCEPTION:
                        # Re-raise the exception we caught in the child process
                        raise value

                assert_never(status)

        return new_fn

    return decorator
