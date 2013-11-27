from __future__ import print_function, division

import os
import wx
import time
import fnmatch

import abipy.gui.awx as awx
import wx.lib.agw.flatnotebook as fnb

from collections import OrderedDict
from pymatgen.io.abinitio.launcher import PyLauncher 
from abipy.gui.events import AbinitEventsFrame, AbinitEventsNotebookFrame
from abipy.gui.timer import MultiTimerFrame
from abipy.gui.browser import FileListFrame, viewerframe_from_filepath
from abipy.gui.editor import TextNotebookFrame, SimpleTextViewer


ID_SHOW_INPUTS = wx.NewId()
ID_SHOW_OUTPUTS = wx.NewId()
ID_SHOW_LOGS = wx.NewId()
ID_SHOW_JOB_SCRIPTS = wx.NewId()
ID_BROWSE = wx.NewId()
ID_SHOW_MAIN_EVENTS = wx.NewId()
ID_SHOW_LOG_EVENTS = wx.NewId()
ID_SHOW_TIMERS = wx.NewId()
ID_CHECK_STATUS = wx.NewId()


class FlowViewerFrame(awx.Frame):
    VERSION = "0.1"

    # Time in second after which we check the status of the tasks.
    REFRESH_INTERVAL = 15

    def __init__(self, parent, flow, **kwargs):
        """
        Args:
            parent:
                Parent window.
            flow:
                `AbinitFlow` with the list of `Workflow` objects.
        """
        if "title" not in kwargs:
            kwargs["title"] = self.codename
            
        super(FlowViewerFrame, self).__init__(parent, -1, **kwargs)

        self.statusbar = self.CreateStatusBar()

        menuBar = wx.MenuBar()

        file_menu = wx.Menu()
        #file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing WFK file")
        #file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the Viewer")
        #file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")
        menuBar.Append(file_menu, "File")

        #file_history = self.file_history = wx.FileHistory(8)
        #self.config = wx.Config(self.codename, style=wx.CONFIG_USE_LOCAL_FILE)
        #file_history.Load(self.config)
        #recent = wx.Menu()
        #file_history.UseMenu(recent)
        #file_history.AddFilesToMenu()
        #file_menu.AppendMenu(wx.ID_ANY, "&Recent Files", recent)
        #self.Bind(wx.EVT_MENU_RANGE, self.OnFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)

        help_menu = wx.Menu()
        help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menuBar.Append(help_menu, "Help")

        self.SetMenuBar(menuBar)

        # Create toolbar.
        self.toolbar = toolbar = self.CreateToolBar()

        #tsize = (48,48)
        #artBmp = wx.ArtProvider.GetBitmap
        toolbar.AddSimpleTool(ID_SHOW_INPUTS, wx.Bitmap(awx.path_img("in.png")), "Visualize the input files of the workflow.")
        toolbar.AddSimpleTool(ID_SHOW_OUTPUTS, wx.Bitmap(awx.path_img("out.png")), "Visualize the output files of the workflow..")
        toolbar.AddSimpleTool(ID_SHOW_LOGS, wx.Bitmap(awx.path_img("log.png")), "Visualize the log files of the workflow.")
        toolbar.AddSimpleTool(ID_SHOW_JOB_SCRIPTS, wx.Bitmap(awx.path_img("script.png")), "Visualize the scripts.")
        toolbar.AddSimpleTool(ID_BROWSE, wx.Bitmap(awx.path_img("browse.png")), "Browse all the files of the workflow.")
        toolbar.AddSimpleTool(ID_SHOW_MAIN_EVENTS, wx.Bitmap(awx.path_img("out_evt.png")), "Show the ABINIT events reported in the main output files.")
        toolbar.AddSimpleTool(ID_SHOW_LOG_EVENTS, wx.Bitmap(awx.path_img("log_evt.png")), "Show the ABINIT events reported in the log files.")
        toolbar.AddSimpleTool(ID_SHOW_TIMERS, wx.Bitmap(awx.path_img("timer.png")), "Show the ABINIT timers in the abo files.")

        toolbar.AddSeparator()
        toolbar.AddSimpleTool(ID_CHECK_STATUS, wx.Bitmap(awx.path_img("refresh.png")), "Check the status of the workflow(s).")

        #toolbar.AddSeparator()
        self.toolbar.Realize()
        self.Centre()

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            #(wx.ID_OPEN, self.OnOpen),
            #(wx.ID_CLOSE, self.OnClose),
            #(wx.ID_EXIT, self.OnExit),
            (wx.ID_ABOUT, self.OnAboutBox),
            #
            (ID_SHOW_INPUTS, self.OnShowInputs),
            (ID_SHOW_OUTPUTS, self.OnShowOutputs),
            (ID_SHOW_LOGS, self.OnShowLogs),
            (ID_SHOW_JOB_SCRIPTS, self.OnShowJobScripts),
            (ID_BROWSE, self.OnBrowse),
            (ID_SHOW_MAIN_EVENTS, self.OnShowMainEvents),
            (ID_SHOW_LOG_EVENTS, self.OnShowLogEvents),
            (ID_SHOW_TIMERS, self.OnShowTimers),
            (ID_CHECK_STATUS, self.OnCheckStatusButton),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        self.flow = flow

        self.check_launcher_file()

        #if filename is not None:
        #    self.ReadWorkflow(filename)
                                         
        self.BuildUi()

    def check_launcher_file(self, with_dialog=True):
        """
        Disable the launch button if we have a sheduler running, 
        since we don't want to have to processes modifying the flow.
        """
        self.disabled_launcher = False
        pid_file = fnmatch.filter(os.listdir(self.flow.workdir), "*.pid")

        if pid_file:
            self.disabled_launcher = True

            pid_file = os.path.join(self.flow.workdir, pid_file[0])

            with open(pid_file, "r") as fh:
                pid = int(fh.readline())

            message = ("Found pid file %s associated to an already running scheduler with pid %d. "
                       "Job submission has been disabled." % (pid_file, pid))

            if with_dialog:
                dialog = wx.MessageDialog(None, message=message, caption='Flow is being executed by a scheduler',  
                                          style=wx.OK | wx.ICON_EXCLAMATION)
                dialog.ShowModal()
                dialog.Destroy()

            else:
                self.statusbar.PushStatusText(message)

    @property
    def codename(self):
        """String with the code name."""
        return self.__class__.__name__ 

    def BuildUi(self):
        self.panel = panel = wx.Panel(self, -1)
        self.main_sizer = main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Here we create a panel and a notebook on the panel
        self.notebook = FlowNotebook(panel, self.flow)
        main_sizer.Add(self.notebook, 1, wx.EXPAND , 5)

        submit_button = wx.Button(panel, -1, label='Submit')
        submit_button.Bind(wx.EVT_BUTTON, self.OnSubmitButton)

        text = wx.StaticText(panel, -1, "Max nlaunch:")
        text.Wrap(-1)
        text.SetToolTipString("Maximum number of tasks that can be submitted. Use -1 for unlimited launches.")
        self.max_nlaunch = wx.SpinCtrl(panel, -1, value="1", min=-1)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(submit_button,  0,  wx.ALIGN_CENTER_HORIZONTAL, 5)
        hsizer.Add(text,  0,  wx.ALIGN_CENTER_HORIZONTAL, 5)
        hsizer.Add(self.max_nlaunch,  0,  wx.ALIGN_CENTER_HORIZONTAL, 5)
        main_sizer.Add(hsizer, 0,  wx.ALIGN_CENTER_HORIZONTAL, 5)                                                                                     

        panel.SetSizerAndFit(main_sizer)

        # Register this event when the GUI is IDLE
        self.last_refresh = time.time()
        self.Bind(wx.EVT_IDLE, self.OnIdle)

    def OnSubmitButton(self, event):
        """Submit up to max_nlauch tasks."""
        if self.disabled_launcher: return
        max_nlaunch = int(self.max_nlaunch.GetValue())
        nlaunch = PyLauncher(self.flow).rapidfire(max_nlaunch=max_nlaunch)

        self.statusbar.PushStatusText("Submitted %d tasks" % nlaunch)

    def GetSelectedWork(self):
        """
        Return the selected workflow namely that the 
        workflow associated to the active tab. None if list is empty.
        """
        return self.notebook.GetSelectedWork()

    def OnIdle(self, event):
        """Functoin executed when the GUI is idle."""
        self.check_launcher_file(with_dialog=False)

        now = time.time()
        if (now - self.last_refresh) > self.REFRESH_INTERVAL:
            self.CheckStatusAndRedraw()
            self.last_refresh = time.time()

    def OnCheckStatusButton(self, event):
        """Callback triggereed by the checkstatus button."""
        self.CheckStatusAndRedraw()

    def CheckStatusAndRedraw(self):
        """Check the status of all the workflows and redraw the panel."""
        self.statusbar.PushStatusText("Checking status...")

        self.flow.check_status()

        # Count the number of tasks with given status.
        counter = self.flow.status_counter

        # Save the active tab so that we can set it afterwards.
        old_selection = self.notebook.GetSelection()

        # Build new notebook and redraw the panel
        main_sizer = self.main_sizer
        main_sizer.Hide(0)
        main_sizer.Remove(0)
        new_notebook = FlowNotebook(self, self.flow)
        main_sizer.Insert(0, new_notebook, 1, wx.EXPAND, 5)
        self.notebook = new_notebook

        self.panel.Layout()

        # Reinstate the old selection
        self.notebook.SetSelection(old_selection)

        # Write number of jobs with given status.
        message = ", ".join("%s: %s" % (k, v) for (k, v) in counter.items())
        self.statusbar.PushStatusText(message)

    #def OnOpen(self, event):
    #    dlg = wx.FileDialog(self, message="Choose a __workflow__.pickle file", defaultDir=os.getcwd(),
    #                        wildcard="pickle files (*.pickle)|*.pickle",
    #                        style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
    #    )

    #    # Show the dialog and retrieve the user response. 
    #    # If it is the OK response, process the data.
    #    if dlg.ShowModal() == wx.ID_OK:
    #        filepath = dlg.GetPath()
    #        self.file_history.AddFileToHistory(filepath)
    #        self.file_history.Save(self.config)
    #        self.config.Flush()
    #        self.ReadWfkFile(filepath)

    #    dlg.Destroy()

    #def OnFileHistory(self, event):
    #    fileNum = event.GetId() - wx.ID_FILE1
    #    filepath = self.file_history.GetHistoryFile(fileNum)
    #    # move up the list
    #    self.file_history.AddFileToHistory(filepath)
    #    self.ReadWfkFile(filepath)

    #def ReadFlowFile(self, filepath):
    #    """Read the WFK file and build the UI."""
    #    self.statusbar.PushStatusText("Reading %s" % filepath)

    #    try:
    #        work = abiopen(filepath)
    #        self.work = wfkfile
    #        self.BuildUi()
    #        self.statusbar.PushStatusText("Workflow file %s loaded" % filepath)

    #    except Exception:
    #        awx.showErrorMessage(self)

    #def OnClose(self, event):
    #    print("onclose")
    #    self.flow.pickle_dump()
    #    self.Close()

    #def OnExit(self, event):
    #    """Save the status of the flow in the database and exit the GUI."""
    #    print("onexit")
    #    self.flow.pickle_dump()
    #    self.Close()

    def OnAboutBox(self, event):
        """"Info on the application."""
        awx.makeAboutBox(codename=self.codename, version=self.VERSION,
                         description="", developers="M. Giantomassi")

    def OnShowInputs(self, event):
        """Show all the input files of the selected `Workflow`."""
        work = self.GetSelectedWork() 
        if work is None: return
        frame = TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.abi")
        frame.Show()

    def OnShowOutputs(self, event):
        """Show all the output files of the selected `Workflow`."""
        work = self.GetSelectedWork() 
        if work is None: return
        frame = TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.abo")
        frame.Show()

    def OnShowJobScripts(self, event):
        """Show all the job script files of the selected `Workflow`."""
        work = self.GetSelectedWork() 
        if work is None: return
        frame = TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.sh")
        frame.Show()

    def OnShowLogs(self, event):
        """Show all the log files of the selected `Workflow`."""
        work = self.GetSelectedWork() 
        if work is None: return
        frame = TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.log")
        frame.Show()

    def OnBrowse(self, event):
        """Browse all the output files produced by the selected `Workflow`."""
        work = self.GetSelectedWork() 
        if work is None: return
        FileListFrame(self, dirpaths=work.workdir).Show()

    def OnShowMainEvents(self, event):
        """Browse all the main events of the tasks in the selected `Workflow`."""
        work = self.GetSelectedWork() 
        if work is None: return
        frame = AbinitEventsNotebookFrame(self, filenames=[task.output_file.path for task in work])
        frame.Show()

    def OnShowLogEvents(self, event):
        """Browse all the log events of the tasks in the selected `Workflow`."""
        work = self.GetSelectedWork() 
        if work is None: return
        frame = AbinitEventsNotebookFrame(self, [task.log_file.path for task in work])
        frame.Show()

    def OnShowTimers(self, event):
        """Analyze the timing data of all the output files of the selected `Workflow`."""
        work = self.GetSelectedWork() 
        if work is None: return
        timers = work.parse_timers()
        # Build the frame for analyzing multiple timers.
        MultiTimerFrame(self, timers).Show()


class FlowNotebook(fnb.FlatNotebook):
    """
    Notebook class
    """
    def __init__(self, parent, flow):
        super(FlowNotebook, self).__init__(parent, id=-1, style=fnb.FNB_NO_X_BUTTON | fnb.FNB_NAV_BUTTONS_WHEN_NEEDED)

        self.flow = flow
        for work in flow:
            tab = TabPanel(self, work)
            self.AddPage(tab, text=os.path.basename(work.workdir))

    def GetSelectedWork(self):
        """
        Return the selected workflow namely that the workflow associated to the 
        active tab. None is list is empy.
        """
        # FIXME: If we want to allow the user to remove some tab we have
        # to remove the corresponding work from workflows.
        # Easy if workflow_viewer can only read data but it might be source 
        # of bugs if we want to modify the object e.g. by submitting the calculation.
        # For the time-being, use this assertion to prevent users from removing pages.
        if self.GetPageCount() != len(self.flow):
            return awx.showErrorMessage(self, message="Bad user has removed pages from the notebook!")

        idx = self.GetSelection()
        if idx == -1:
            return None
                                                                                                       
        try:
            return self.flow[idx]
        except IndexError:
            return None


class TabPanel(wx.Panel):
    """
    Notebook tab.
    """
    def __init__(self, parent, work, **kwargs):
        wx.Panel.__init__(self, parent=parent, id=-1, **kwargs)

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # List Control with the individual tasks of the workflow.
        task_listctrl = TaskListCtrl(self, work)
        main_sizer.Add(task_listctrl, 1, wx.EXPAND, 5)

        label = wx.StaticText(self, -1, "Workflow class %s, status: %s, finalized: %s" % (
            work.__class__.__name__, work.status, work.finalized))
        label.Wrap(-1)
        main_sizer.Add(label, 0, wx.ALIGN_LEFT, 5)

        self.SetSizerAndFit(main_sizer)


class TaskListCtrl(wx.ListCtrl):
    """
    ListCtrl with the list of task in the `Workflow`.
    """
    def __init__(self, parent, work, **kwargs):
        """
        Args:
            parent:
                Parent window.
            work:
                `Workflow` containig the List of `Task` instances.
        """
        super(TaskListCtrl, self).__init__(parent, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN, **kwargs)

        self.work = work

        columns = ["Task", "Status", "Queue_id", 
                   "Errors", "Warnings", "Comments", 
                   "MPI", "OMP", 
                   "num_restarts", "Task Class",
                   ]

        for (index, col) in enumerate(columns):
            self.InsertColumn(index, col)

        # Used to store the Max width in pixels for the data in the column.
        column_widths = [awx.get_width_height(self, s)[0] for s in columns]

        for task in work:
            events = map(str, 3*["N/A"])

            try:
                report = task.get_event_report()
                if report is not None: 
                    events = map(str, [report.num_errors, report.num_warnings, report.num_comments])
            except:
                pass

            cpu_info = [task.mpi_ncpus, task.omp_ncpus]
            entry = map(str, [task.name, str(task.status), task.queue_id] + 
                              events + 
                              cpu_info + 
                              [task.num_restarts, task.__class__.__name__]
                        )
            w = [awx.get_width_height(self, s)[0] for s in entry]
            column_widths = map(max, zip(w, column_widths))

            self.Append(entry)

        # Set the width in pixel for each column.
        for (index, col) in enumerate(columns):
            self.SetColumnWidth(index, column_widths[index])
            #self.SetColumnWidth(index, wx.LIST_AUTOSIZE)

        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        #self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def OnRightClick(self, event):
        currentItem = event.m_itemIndex
        if currentItem == -1:
            return

        # Open the popup menu then destroy it to avoid mem leak.
        task = self.work[currentItem]
        menu = TaskPopupMenu(parent=self, task=task)
        self.PopupMenu(menu, event.GetPoint())
        menu.Destroy()

    #def OnItemActivated(self, event):
    #    currentItem = event.m_itemIndex
    #    task = self.work[currentItem]
    #    if task.can_run:
    #        task.start()
    #    # This is to update the database.
    #    self.flow.pickle_dump()


# Callbacks 

def show_task_main_output(parent, task):
    file = task.output_file

    if file.exists:
        frame = viewerframe_from_filepath(parent, file.path)
        frame.Show()
    else:
        awx.showErrorMessage(parent=parent, message="Output file %s does not exist" % file.path)


def show_task_log(parent, task):
    file = task.log_file

    if file.exists:
        frame = viewerframe_from_filepath(parent, file.path)
        frame.Show()
    else:
        awx.showErrorMessage(parent=parent, message="Output file %s does not exist" % file.path)


def show_task_main_events(parent, task):
    file = task.output_file

    if file.exists:
        AbinitEventsFrame(parent, file.path).Show()
    else:
        awx.showErrorMessage(parent=parent, message="Output file %s does not exist" % file.path)


def show_task_log_events(parent, task):
    file = task.log_file

    if file.exists:
        AbinitEventsFrame(parent, file.path).Show()
    else:
        awx.showErrorMessage(parent=parent, message="Log file %s does not exist" % file.path)


def browse_indir(parent, task):
    """Open a window that allows the user to browse the input files in indir."""
    FileListFrame(parent, dirpaths=task.indir.path).Show()


def browse_outdir(parent, task):
    """Open a window that allows the user to browse the output files in outdir."""
    FileListFrame(parent, dirpaths=task.outdir.path).Show()

def browse_tmpdir(parent, task):
    """Open a window that allows the user to browse the output files in outdir."""
    FileListFrame(parent, dirpaths=task.tmpdir.path).Show()

def show_history(parent, task):
    """Show the history of the task."""
    text = "\n".join(task.history)
    SimpleTextViewer(parent, text).Show()


def check_status_and_pickle(task):
    task.flow.check_status()
    task.flow.pickle_dump()


def task_restart(parent, task):
    """Restart the task."""
    task.restart()
    check_status_and_pickle(task)


def task_reset(parent, task):
    """Reset the status of the task."""
    task.reset()
    check_status_and_pickle(task)


def task_show_deps(parent, task):
    """Show the dependencies of the task."""
    text = task.str_deps()
    SimpleTextViewer(parent, text).Show()


def task_inspect(parent, task):
    """Inspect the results at runtime."""
    if hasattr(task, "inspect"):
        task.inspect()


class TaskPopupMenu(wx.Menu):
    """
    A `TaskPopupMenu` has a list of callback functions indexed by the menu title. 
    The signature of the callback function is func(parent, task) where parent is 
    the wx Window that will become the parent of the new frame created by the callback.
    and task is a `Task` instance.
    """
    MENU_TITLES = OrderedDict([
        ("output", show_task_main_output),
        ("log", show_task_log),
        ("main events", show_task_main_events),
        ("log events",  show_task_log_events),
        ("browse indir", browse_indir),
        ("browse outdir", browse_outdir),
        ("browse tmpdir", browse_tmpdir),
        ("history", show_history),
        ("restart", task_restart),
        ("reset", task_reset),
        ("dependencies", task_show_deps),
        ("inspect", task_inspect),
    ])

    def __init__(self, parent, task):
        super(TaskPopupMenu, self).__init__()
        self.parent, self.task = parent, task

        self._make_menu()

    def _make_menu(self):
        """Build the menu"""
        self.menu_title_by_id = OrderedDict()

        for title in self.MENU_TITLES:
            self.menu_title_by_id[wx.NewId()] = title

        for (id, title) in self.menu_title_by_id.items():
            # Register menu handlers with EVT_MENU, on the menu.
            self.Append(id, title)
            wx.EVT_MENU(self, id, self.OnMenuSelection)

    def _get_callback(self, title):
        return self.MENU_TITLES[title]

    def OnMenuSelection(self, event):
        title = self.menu_title_by_id[event.GetId()]
        callback = self._get_callback(title)
        #print("Calling callback %s on task %s" % (callback, self.task))
        try:
            callback(self.parent, self.task)
        except:
            awx.showErrorMessage(parent=self.parent)


def wxapp_flow_viewer(works):
    """Standalone application for `FlowViewerFrame"""
    app = awx.App()
    FlowViewerFrame(None, works).Show()
    return app

