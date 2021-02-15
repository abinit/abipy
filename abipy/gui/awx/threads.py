"""Inspired to http://wiki.wxpython.org/LongRunningTasks."""
from __future__ import print_function, division

import time
import wx

from threading import Thread

__all__ = [
"WorkerThread",
]

# Define notification event for thread completion
EVT_THREAD_RESULT_ID = wx.NewId()

def EVT_THREAD_RESULT(win, func):
    """Define Thread_Result Event."""
    win.Connect(-1, -1, EVT_THREAD_RESULT_ID, func)


class ThreadResultEvent(wx.PyEvent):
    """Simple event to carry arbitrary result data."""

    def __init__(self, **data):
        """Init Result Event."""
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_THREAD_RESULT_ID)
        self.data = data


class WorkerThread(Thread):
    """Worker Thread Class that executes processing."""

    def __init__(self, notify_window, group=None, target=None, name=None, *args, **kwargs):
        """Init Worker Thread Class."""
        super(WorkerThread,self).__init__(group=group, target=target, name=name, args=args, kwargs=kwargs)

        self._notify_window = notify_window
        self._want_abort = 0

    def start(self):
        """Start Worker Thread."""
        super(WorkerThread, self).start()

    def run(self):
        """Run Worker Thread."""
        # This is the code executing in the new thread. Simulation of
        # a long process (well, 10s here) as a simple loop - you will
        # need to structure your processing so that you periodically
        # peek at the abort variable
        retcode = super(WorkerThread, self).run()
        
        #while True:
        #    time.sleep(1)
        #    if self._want_abort:
        #        # Use a result of None to acknowledge the abort (of
        #        # course you can use whatever you'd like or even a separate event type)
        #        wx.PostEvent(self._notify_window, ThreadResultEvent(retcode="abort"))
        #        return

        # Here's where the result would be returned (this is an
        # example fixed result of the number 10, but it could be any Python object)
        wx.PostEvent(self._notify_window, ThreadResultEvent(retcode=retcode))
        return retcode 

    def abort(self):
        """Abort worker thread."""
        # Method for use by main thread to signal an abort
        self._want_abort = 1

