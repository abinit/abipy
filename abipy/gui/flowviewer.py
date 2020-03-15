import os
import wx
import time
import fnmatch

import abipy.gui.awx as awx
import wx.lib.agw.flatnotebook as fnb

from collections import OrderedDict
from monty.dev import get_ncpus
from abipy.gui.events import AbinitEventsFrame, AbinitEventsNotebookFrame
from abipy.gui.timer import MultiTimerFrame, AbinitTimerFrame
from abipy.gui.browser import FileListFrame, DirBrowserFrame, frame_from_filepath, frameclass_from_filepath
from abipy.gui.editor import TextNotebookFrame, SimpleTextViewer
import abipy.flowtk as flowtk


def yaml_manager_dialog(parent):
    """
    Open a dialog that allows the user to select a YAML file with the taskmanager.
    Returns the new manager or None if the user clicked CANCEL or the specifed file is not valid.
    """
    dialog = wx.FileDialog(parent, message="Choose a taskmanager.yml file", defaultDir=os.getcwd(),
                           wildcard="YAML files (*.yml)|*.yml",
                           style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR)

    # Show the dialog and retrieve the user response.
    # If it is the OK response, process the data.
    if dialog.ShowModal() == wx.ID_CANCEL: return None
    filepath = dialog.GetPath()
    dialog.Destroy()

    try:
        return flowtk.TaskManager.from_file(filepath)
    except:
        awx.showErrorMessage(parent)
        return None


class FlowViewerFrame(awx.Frame):
    VERSION = "0.1"

    # Time in second after which we check the status of the tasks.
    REFRESH_INTERVAL = 120

    HELP_MSG = """\
Quick help:

  Task list:

     Left-Click:   Open directory with output files.
     Right-Click:  display popup menu with choices.

Also, these key bindings can be used
(For Mac OSX, replace 'Ctrl' with 'Apple'):

  Ctrl-Q:     quit"""

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

        # This combination of options for config seems to work on my Mac.
        self.config = wx.FileConfig(appName=self.codename, localFilename=self.codename + ".ini",
                                    style=wx.CONFIG_USE_LOCAL_FILE)

        # Build menu, toolbar and status bar.
        self.SetMenuBar(self.makeMenu())
        self.makeToolBar()
        self.statusbar = self.CreateStatusBar()
        self.Centre()

        self.flow = flow

        # Disable launch mode if we already executing the flow with the scheduler.
        self.check_launcher_file()

        # Build UI
        self.panel = panel = wx.Panel(self, -1)
        self.main_sizer = main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Here we create a panel and a notebook on the panel
        self.notebook = FlowNotebook(panel, self.flow)
        main_sizer.Add(self.notebook, 1, wx.EXPAND, 5)

        submit_button = wx.Button(panel, -1, label='Submit')
        submit_button.Bind(wx.EVT_BUTTON, self.OnSubmitButton)

        text = wx.StaticText(panel, -1, "Max nlaunch:")
        text.Wrap(-1)
        text.SetToolTipString("Maximum number of tasks that can be submitted. Use -1 for unlimited launches.")
        self.max_nlaunch = wx.SpinCtrl(panel, -1, value=str(get_ncpus()), min=-1)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(submit_button, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        hsizer.Add(text, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        hsizer.Add(self.max_nlaunch, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)
        main_sizer.Add(hsizer, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        panel.SetSizerAndFit(main_sizer)

        # Register this event when the GUI is IDLE
        # Not used anymore since the Yaml parser
        # is very slow when we use the pure python version.
        self.last_refresh = time.time()
        #self.Bind(wx.EVT_IDLE, self.OnIdle)

    @property
    def codename(self):
        """String with the code name."""
        return "Flow Viewer"

    def makeToolBar(self):
        """Create toolbar."""
        self.toolbar = toolbar = self.CreateToolBar()
        toolbar.SetToolBitmapSize(wx.Size(48, 48))

        def bitmap(path):
            return wx.Bitmap(awx.path_img(path))

        self.ID_SHOW_INPUTS = wx.NewId()
        self.ID_SHOW_OUTPUTS = wx.NewId()
        self.ID_SHOW_LOGS = wx.NewId()
        self.ID_SHOW_JOB_SCRIPTS = wx.NewId()
        self.ID_SHOW_JOB_OUTERRS = wx.NewId()
        self.ID_BROWSE = wx.NewId()
        self.ID_SHOW_MAIN_EVENTS = wx.NewId()
        self.ID_SHOW_LOG_EVENTS = wx.NewId()
        self.ID_SHOW_TIMERS = wx.NewId()
        self.ID_CHECK_STATUS = wx.NewId()

        toolbar.AddSimpleTool(self.ID_SHOW_INPUTS, bitmap("in.png"), "Visualize the input file(s) of the workflow.")
        toolbar.AddSimpleTool(self.ID_SHOW_OUTPUTS, bitmap("out.png"), "Visualize the output file(s) of the workflow..")
        toolbar.AddSimpleTool(self.ID_SHOW_LOGS, bitmap("log.png"), "Visualize the log file(s) of the workflow.")
        toolbar.AddSimpleTool(self.ID_SHOW_JOB_SCRIPTS, bitmap("script.png"), "Visualize the script(s).")
        toolbar.AddSimpleTool(self.ID_SHOW_JOB_OUTERRS, bitmap("script.png"), "Visualize the script(s).")
        toolbar.AddSimpleTool(self.ID_BROWSE, bitmap("browse.png"), "Browse all the files of the workflow.")
        toolbar.AddSimpleTool(self.ID_SHOW_MAIN_EVENTS, bitmap("out_evt.png"), "Show the ABINIT events reported in the main output files.")
        toolbar.AddSimpleTool(self.ID_SHOW_LOG_EVENTS, bitmap("log_evt.png"), "Show the ABINIT events reported in the log files.")
        toolbar.AddSimpleTool(self.ID_SHOW_TIMERS, bitmap("timer.png"), "Show the ABINIT timers in the abo files.")

        toolbar.AddSeparator()

        toolbar.AddSimpleTool(self.ID_CHECK_STATUS, bitmap("refresh.png"), "Check the status of the workflow(s).")

        #toolbar.AddSeparator()
        #self.file_selector = FileSelector(toolbar, self)
        #toolbar.AddControl(control=self.file_selector)

        toolbar.Realize()

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_SHOW_INPUTS, self.OnShowInputs),
            (self.ID_SHOW_OUTPUTS, self.OnShowOutputs),
            (self.ID_SHOW_LOGS, self.OnShowLogs),
            (self.ID_SHOW_JOB_SCRIPTS, self.OnShowJobScripts),
            (self.ID_SHOW_JOB_OUTERRS, self.OnShowJobOutErrs),
            (self.ID_BROWSE, self.OnBrowse),
            (self.ID_SHOW_MAIN_EVENTS, self.OnShowMainEvents),
            (self.ID_SHOW_LOG_EVENTS, self.OnShowLogEvents),
            (self.ID_SHOW_TIMERS, self.OnShowTimers),
            (self.ID_CHECK_STATUS, self.OnCheckStatusButton),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

    def makeMenu(self):
        """Make the menu bar."""
        menu_bar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open", help="Open an existing flow in a new frame")
        file_menu.Append(wx.ID_CLOSE, "&Close", help="Close the file associated to the active tab")
        file_menu.Append(wx.ID_EXIT, "&Quit", help="Exit the application")

        file_history = self.file_history = wx.FileHistory(8)
        file_history.Load(self.config)
        recent = wx.Menu()
        file_history.UseMenu(recent)
        file_history.AddFilesToMenu()
        file_menu.AppendMenu(-1, "&Recent Files", recent)
        self.Bind(wx.EVT_MENU_RANGE, self.OnFileHistory, id=wx.ID_FILE1, id2=wx.ID_FILE9)
        menu_bar.Append(file_menu, "File")

        flow_menu = wx.Menu()

        #self.ID_FLOW_CHANGE_MANAGER = wx.NewId()
        #flow_menu.Append(self.ID_FLOW_CHANGE_MANAGER, "Change TaskManager", help="")

        self.ID_FLOW_OPEN_OUTFILES = wx.NewId()
        flow_menu.Append(self.ID_FLOW_OPEN_OUTFILES, "Open output files", help="")

        self.ID_FLOW_TREE_VIEW = wx.NewId()
        flow_menu.Append(self.ID_FLOW_TREE_VIEW, "Tree view", help="Tree view of the tasks")

        menu_bar.Append(flow_menu, "Flow")

        help_menu = wx.Menu()
        help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menu_bar.Append(help_menu, "Help")

        self.ID_HELP_QUICKREF = wx.NewId()
        help_menu.Append(self.ID_HELP_QUICKREF, "Quick Reference ", help="Quick reference for " + self.codename)
        help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (wx.ID_OPEN, self.OnOpen),
            (wx.ID_CLOSE, self.OnClose),
            (wx.ID_EXIT, self.OnExit),
            (wx.ID_ABOUT, self.OnAboutBox),
            #
            #(self.ID_FLOW_CHANGE_MANAGER, self.onChangeManager),
            (self.ID_FLOW_OPEN_OUTFILES, self.onOpenOutFiles),
            (self.ID_FLOW_TREE_VIEW, self.onTaskTreeView),
            #
            (self.ID_HELP_QUICKREF, self.onQuickRef),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        return menu_bar

    def check_launcher_file(self, with_dialog=True):
        """
        Disable the launch button if we have a sheduler running,
        since we don't want to have processes modifying the flow.
        """
        self.disabled_launcher = False
        pid_files = fnmatch.filter(os.listdir(self.flow.workdir), "*.pid")

        if pid_files:
            pid_file = os.path.join(self.flow.workdir, pid_files[0])
            if not os.path.exists(pid_file): return
            self.disabled_launcher = True

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

    def OnSubmitButton(self, event):
        """Submit up to max_nlauch tasks."""
        if self.disabled_launcher: return
        max_nlaunch = int(self.max_nlaunch.GetValue())
        nlaunch = flowtk.PyLauncher(self.flow).rapidfire(max_nlaunch=max_nlaunch)

        self.statusbar.PushStatusText("Submitted %d tasks" % nlaunch)

    def GetSelectedWork(self):
        """
        Return the selected workflow namely that the
        workflow associated to the active tab. None if list is empty.
        """
        return self.notebook.GetSelectedWork()

    def OnIdle(self, event):
        """Function executed when the GUI is idle."""
        now = time.time()
        if (now - self.last_refresh) > self.REFRESH_INTERVAL:
            self.check_launcher_file(with_dialog=False)
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

    def read_file(self, filepath):
        self.statusbar.PushStatusText("Reading %s" % filepath)

        try:
            flow = flowtk.AbinitFlow.pickle_load(filepath)
            self.AddFileToHistory(filepath)
            return flow
        except:
            awx.showErrorMessage(self)
            return None

    def OnOpen(self, event):
        dialog = wx.FileDialog(
            self, message="Choose a __workflow__.pickle file", defaultDir=os.getcwd(),
            wildcard="pickle files (*.pickle)|*.pickle",
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR)

        # Show the dialog and retrieve the user response.
        # If it is the OK response, process the data.
        if dialog.ShowModal() == wx.ID_CANCEL: return

        filepath = dialog.GetPath()
        dialog.Destroy()

        # Add to the history.
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()
        flow = self.read_file(filepath)

        if flow is not None:
            # Open new frame.
            FlowViewerFrame(self, flow).Show()

    def AddFileToHistory(self, filepath):
        """Add the absolute filepath to the file history."""
        self.file_history.AddFileToHistory(filepath)
        self.file_history.Save(self.config)
        self.config.Flush()

    def OnFileHistory(self, event):
        fileNum = event.GetId() - wx.ID_FILE1
        filepath = self.file_history.GetHistoryFile(fileNum)
        self.file_history.AddFileToHistory(filepath)
        self.read_file(filepath)

    def OnClose(self, event):
        print("onclose")
        #self.flow.pickle_dump()
        #self.Close()

    def OnExit(self, event):
        """Save the status of the flow in the database and exit the GUI."""
        print("onexit")
        #self.flow.pickle_dump()
        self.Close()

    def OnAboutBox(self, event):
        """"Info on the application."""
        awx.makeAboutBox(codename=self.codename, version=self.VERSION,
                         description="", developers="M. Giantomassi")

    def onQuickRef(self, event=None):
        dialog = wx.MessageDialog(self, self.HELP_MSG, self.codename + " Quick Reference",
                               wx.OK | wx.ICON_INFORMATION)
        dialog.ShowModal()
        dialog.Destroy()

    def OnShowInputs(self, event):
        """Show all the input files of the selected `Workflow`."""
        work = self.GetSelectedWork()
        if work is None: return
        TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.abi").Show()

    def OnShowOutputs(self, event):
        """Show all the output files of the selected `Workflow`."""
        work = self.GetSelectedWork()
        if work is None: return
        TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.abo").Show()

    def OnShowJobScripts(self, event):
        """Show all the job script files of the selected `Workflow`."""
        work = self.GetSelectedWork()
        if work is None: return
        TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.sh").Show()

    def OnShowJobOutErrs(self, event):
        """Show all the queue output/error files files of the selected `Workflow`."""
        work = self.GetSelectedWork()
        if work is None: return
        TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.qout|*.qerr").Show()

    def OnShowLogs(self, event):
        """Show all the log files of the selected `Workflow`."""
        work = self.GetSelectedWork()
        if work is None: return
        TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.log").Show()

    def OnBrowse(self, event):
        """Browse all the output files produced by the selected `Workflow`."""
        work = self.GetSelectedWork()
        if work is None: return
        FileListFrame(self, dirpaths=work.workdir).Show()
        #DirBrowserFrame(self, dir=work.workdir).Show()

    def OnShowMainEvents(self, event):
        """Browse all the main events of the tasks in the selected `Workflow`."""
        work = self.GetSelectedWork()
        if work is None: return
        AbinitEventsNotebookFrame(self, filenames=[task.output_file.path for task in work]).Show()

    def OnShowLogEvents(self, event):
        """Browse all the log events of the tasks in the selected `Workflow`."""
        work = self.GetSelectedWork()
        if work is None: return
        AbinitEventsNotebookFrame(self, [task.log_file.path for task in work]).Show()

    def OnShowTimers(self, event):
        """Analyze the timing data of all the output files of the selected `Workflow`."""
        work = self.GetSelectedWork()
        if work is None: return
        timers = work.parse_timers()
        # Build the frame for analyzing multiple timers.
        MultiTimerFrame(self, timers).Show()

    def onTaskTreeView(self, event):
        TaskTreeView(self, self.flow).Show()

    #def onChangeManager(self, event):
    #    #ChangeTaskManager(self, self.flow).Show()
    #    new_manager = yaml_manager_dialog(self)
    #    if new_manager is None: return
    #    print(new_manager)
    #    #status_selected =  upper()
    #    #status = None if status_selected == "ALL" else status_selected
    #    # Change the manager of the errored tasks.
    #    #for task in flow.iflat_tasks(status="S_ERROR"):
    #    #    task.reset()
    #    #    task.set_manager(new_manager)

    def onOpenOutFiles(self, event):
        FileSelectorFrame(parent=self, viewer=self).Show()


class FlowNotebook(fnb.FlatNotebook):
    """
    Notebook class
    """
    def __init__(self, parent, flow):
        try:
            style = fnb.FNB_NO_X_BUTTON | fnb.FNB_NAV_BUTTONS_WHEN_NEEDED
        except AttributeError:
            style = fnb.FNB_NO_X_BUTTON

        super(FlowNotebook, self).__init__(parent, id=-1, style=style)

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
        if idx == -1: return None
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

        # List Control with the individual tasks of the workflow.
        task_listctrl = TaskListCtrl(self, work)

        label = wx.StaticText(self, -1, "Workflow class %s, status: %s, finalized: %s" % (
            work.__class__.__name__, work.status, work.finalized))
        label.Wrap(-1)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(label, 0, wx.ALIGN_LEFT, 5)
        main_sizer.Add(task_listctrl, 1, wx.EXPAND, 5)

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

        for index, col in enumerate(columns):
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

            cpu_info = [task.mpi_procs, task.omp_threads]
            entry = map(str, [task.name, str(task.status), task.queue_id] +
                              events +
                              cpu_info +
                              [task.num_restarts, task.__class__.__name__]
                        )
            w = [awx.get_width_height(self, s)[0] for s in entry]
            column_widths = map(max, zip(w, column_widths))

            self.Append(entry)

        # Set the width in pixel for each column.
        for index, col in enumerate(columns):
            self.SetColumnWidth(index, column_widths[index])

        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def OnRightClick(self, event):
        currentItem = event.m_itemIndex
        if currentItem == -1: return

        # Open the popup menu then destroy it to avoid mem leak.
        task = self.work[currentItem]
        menu = TaskPopupMenu(parent=self, task=task)
        self.PopupMenu(menu, event.GetPoint())
        menu.Destroy()

    def OnItemActivated(self, event):
        """
        Browse the outdir of the task
        """
        currentItem = event.m_itemIndex
        task = self.work[currentItem]
        FileListFrame(self, dirpaths=task.outdir.path, title=task.outdir.relpath).Show()


# Callbacks

def show_task_main_output(parent, task):
    """Show the main output file of the task."""
    file = task.output_file

    if file.exists:
        frame_from_filepath(parent, file.path).Show()
    else:
        awx.showErrorMessage(parent=parent, message="Output file %s does not exist" % file.path)


def show_task_log(parent, task):
    """Show the log file of the task."""
    file = task.log_file

    if file.exists:
        frame_from_filepath(parent, file.path).Show()
    else:
        awx.showErrorMessage(parent=parent, message="Output file %s does not exist" % file.path)


def show_task_main_events(parent, task):
    """Show the events reported in the main output file."""
    file = task.output_file

    if file.exists:
        AbinitEventsFrame(parent, file.path).Show()
    else:
        awx.showErrorMessage(parent=parent, message="Output file %s does not exist" % file.path)


def show_task_log_events(parent, task):
    """Show the events reported in the log file."""
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


def show_timer(parent, task):
    """Show timing data of the k."""
    try:
        frame = AbinitTimerFrame(parent, task.output_file.path)
        frame.Show()
    except awx.Error as exc:
        awx.showErrorMessage(parent, str(exc))


def check_status_and_pickle(task):
    """Check the status of the task and update the pickle database."""
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


def task_set_status(parent, task):
    """Reset the status of the task."""
    choices = [str(s) for s in flowtk.Node.ALL_STATUS]
    dialog = wx.SingleChoiceDialog(parent, message="Select new status", caption="", choices=choices)
    if dialog.ShowModal() == wx.ID_CANCEL: return None
    status = choices[dialog.GetSelection()]
    dialog.Destroy()

    task.set_status(status, info_msg="Status changed by user on %s" % time.asctime())
    #task.reset()
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
        ("set status", task_set_status),
        ("timer", show_timer),
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

#class ChangeTaskManager(awx.Frame):
#    def __init__(self, parent, flow, **kwargs)
#        super(ChangeTaskManager, self).__init__(self, parent, **kwargs)
#
#
#    def onOkButton(self, event):
#        new_manager = yaml_manager_dialog(self)
#        if new_manager is None: return
#        print(new_manager)
#        #status_selected =  upper()
#        #status = None if status_selected == "ALL" else status_selected
#        # Change the manager of the errored tasks.
#        #for task in flow.iflat_tasks(status="S_ERROR"):
#        #    task.reset()
#        #    task.set_manager(new_manager)


class TaskStatusTreePanel(awx.Panel):
    """
    Panel with a TreeCtrl that allows the user to navigate the tasks grouped by status.
    """
    def __init__(self, parent, flow, **kwargs):
        super(TaskStatusTreePanel, self).__init__(parent, -1, **kwargs)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        vbox = wx.BoxSizer(wx.VERTICAL)
        panel1 = awx.Panel(self, -1)
        panel2 = awx.Panel(self, -1)

        self.tree = tree = wx.TreeCtrl(panel1, 1, wx.DefaultPosition, (-1, -1), wx.TR_HIDE_ROOT | wx.TR_HAS_BUTTONS)

        root = self.tree.AddRoot('Task Status')

        # entry = collections.namedtuple("Entry", "task wi ti")
        #print(status2entries)
        self.status2entries = flow.groupby_status()

        self.status_nodes = []
        for status, entries in self.status2entries.items():
            node = tree.AppendItem(root, "%d %s tasks" % (len(entries), str(status)), data=wx.TreeItemData(status))
            self.status_nodes.append(node)
            for entry in entries:
                tree.AppendItem(node, "Task: " + str(entry.task), data=wx.TreeItemData(entry))

        tree.Bind(wx.EVT_TREE_SEL_CHANGED, self.onSelChanged)
        tree.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.onItemRightClick)

        self.display = wx.StaticText(panel2, -1, '', (10, 10), style=wx.ALIGN_LEFT)

        vbox.Add(self.tree, 1, wx.EXPAND)
        main_sizer.Add(panel1, 1, wx.EXPAND)
        main_sizer.Add(panel2, 1, wx.EXPAND)
        panel1.SetSizerAndFit(vbox)

        self.SetSizerAndFit(main_sizer)
        self.Centre()

    def onSelChanged(self, event):
        node = event.GetItem()
        if node in self.status_nodes: return

        proxy = self.tree.GetItemData(node)
        if proxy is None: return
        entry = proxy.GetData()

        task = entry.task
        s = str(task)
        s += "\nHistory:\n" + task.str_history

        self.display.SetLabel(s)

    def onItemRightClick(self, event):
        node = event.GetItem()
        proxy = self.tree.GetItemData(node)
        if proxy is None: return

        if node in self.status_nodes:
            status = proxy.GetData()
            print("received set of tasks with status %s" % status)
            popup_menu = self.makePopupMenuStatus()
            self.PopupMenu(popup_menu, event.GetPoint())
            popup_menu.Destroy()

        #print("click")
        #print("event", dir(event))

    def makePopupMenuStatus(self):
        self.ID_POPUP_CHANGE_MANAGER = wx.NewId()
        menu = wx.Menu()
        menu.Append(self.ID_POPUP_CHANGE_MANAGER, "Change manager")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_POPUP_CHANGE_MANAGER, self.onChangeManager),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        return menu

    def onChangeManager(self, event):
        print("changer manager")
        node = self.tree.GetSelection()
        status = self.tree.GetItemData(node).GetData()
        print("status", status)

        entries = self.status2entries[status]

        new_manager = yaml_manager_dialog(self)
        for e in entries:
            e.task.reset()
            e.task.set_manager(new_manager)

        #self.flow.build_and_pickle_dump()


class TaskTreeView(awx.Frame):
    """
    A frame with a TaskStatusTreePanel.
    """
    def __init__(self, parent, flow, **kwargs):
        if "title" not in kwargs:
            kwargs["title"] = "Task tree view: %s" % flow.workdir

        super(TaskTreeView, self).__init__(parent, **kwargs)

        panel = TaskStatusTreePanel(self, flow)


class FileSelector(awx.Panel):
#class FileSelector(wx.Control):
    """
    Control that allows the user to select multiple output files of the same type (either inside
    a `Workflow` on in the entire `Flow`. The open button will open a viewer to analyze
    the multiple files selected.
    """
    def __init__(self, parent, viewer, **kwargs):
        super(FileSelector, self).__init__(parent, -1, **kwargs)
        self.viewer = viewer
        panel = self #wx.Panel(self, -1)

        wcards = ["*GSR.nc", "*WFK-etsf.nc", "*SIGRES.nc", "*MDF.nc"]

        self.wcards_cbox = wx.ComboBox(panel, id=-1, name='File type', choices=wcards, value=wcards[0], style=wx.CB_READONLY)

        smodes = ["Selected Workflow", "Entire Flow"]
        self.select_rbox = wx.RadioBox(panel, id=1, name="Selection Mode", choices=smodes, style=wx.RA_SPECIFY_ROWS)

        open_button = wx.Button(panel, -1, label='Open files')
        open_button.Bind(wx.EVT_BUTTON, self.onOpenButton)

        #main_sizer = wx.BoxSizer(wx.HORIZONTAL)
        main_sizer = wx.GridBagSizer(10, 10)

        #vsizer = wx.BoxSizer(wx.VERTICAL)
        #vsizer.Add(self.wcards_cbox, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        #vsizer.Add(open_button, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        #main_sizer.Add(vsizer)
        #main_sizer.Add(self.select_rbox, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        main_sizer.Add(self.wcards_cbox, (0, 0), (1,1), wx.ALIGN_CENTER)
        main_sizer.Add(open_button, (1, 0), (1,1), wx.ALIGN_CENTER)
        main_sizer.Add(self.select_rbox, (0, 1), (3,3), wx.EXPAND)

        panel.SetSizerAndFit(main_sizer)

    def getWildCard(self):
        """Returns a string with the Abinit extension selected by the user."""
        return self.wcards_cbox.GetValue()

    def getSelectionMode(self):
        """Returns a string with the selection mode selected by the user."""
        index = self.select_rbox.GetSelection()
        return self.select_rbox.GetString(index)

    def onOpenButton(self, event):
        wcard = self.getWildCard()
        smode = self.getSelectionMode()
        print("wcard", wcard, "smode", smode)

        # Find the filepaths according to smode.
        filepaths = []

        if smode == "Selected Workflow":
            work = self.viewer.GetSelectedWork()
            for task in work:
                filepaths.extend(task.outdir.list_filepaths(wildcard=wcard))

        elif smode == "Entire Flow":
            for work in self.viewer.flow:
                for task in work:
                    filepaths.extend(task.outdir.list_filepaths(wildcard=wcard))

        else:
            return awx.showErrorMessage(self, "Wrong value of smode: %s" % smode)

        if not filepaths:
            return awx.showErrorMessage(self, "Cannot find any file matching the specified shell pattern")

        print("filepaths", filepaths)

        # Get the viewer class associated to these files, build the frame and show it.
        frame_class = frameclass_from_filepath(filepaths[0])
        if frame_class is None: return
        frame_class(self, filepaths).Show()


class FileSelectorFrame(wx.Frame):
    def __init__(self, parent, viewer, **kwargs):
        super(FileSelectorFrame, self).__init__(parent, -1, **kwargs)

        panel = wx.Panel(self, -1)
        file_selector = FileSelector(panel, viewer)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(file_selector, 1, wx.EXPAND , 5)
        panel.SetSizerAndFit(main_sizer)


#class StatusSelectorFrame(wx.Frame):
#    def __init__(self, parent, viewer, **kwargs):
#        super(StatusSelectorFrame, self).__init__(parent, -1, **kwargs)
#
#        panel = wx.Panel(self, -1)
#
#        #choices = ["S_OK", "*WFK-etsf.nc", "*SIGRES.nc", "*MDF.nc"]
#
#        #self.status_cbox = wx.ComboBox(panel, id=-1, name='File type', choices=choices, value=choices[0], style=wx.CB_READONLY)
#
#        #smodes = ["Selected Workflow", "Entire Flow"]
#        #self.select_rbox = wx.RadioBox(panel, id=1, name="Selection Mode", choices=smodes, style=wx.RA_SPECIFY_ROWS)
#
#        #open_button = wx.Button(panel, -1, label='Open files')
#        #open_button.Bind(wx.EVT_BUTTON, self.onOpenButton)
#
#        main_sizer = wx.BoxSizer(wx.VERTICAL)
#        main_sizer.Add(file_selector, 1, wx.EXPAND , 5)
#        panel.SetSizerAndFit(main_sizer)
#
#    def getSelectedStatus(self):
#        return self.choices_cbox.GetValue()


def wxapp_flow_viewer(works):
    """Standalone application for `FlowViewerFrame`."""
    app = awx.App()
    FlowViewerFrame(None, works).Show()
    return app
