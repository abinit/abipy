from __future__ import print_function, division

import os
import wx

import abipy.gui.awx as awx

from collections import OrderedDict
from abipy.gui.events import AbinitEventsFrame, AbinitEventsNotebookFrame
from abipy.gui.browser import FileListFrame, viewerframe_from_filepath
from abipy.gui.editor import TextNotebookFrame
import wx.lib.agw.flatnotebook as fnb

ID_SHOW_INPUTS = wx.NewId()
ID_SHOW_OUTPUTS = wx.NewId()
ID_SHOW_LOGS = wx.NewId()
ID_BROWSE = wx.NewId()
ID_SHOW_MAIN_EVENTS = wx.NewId()
ID_SHOW_LOG_EVENTS = wx.NewId()


_FRAME_SIZE = (720, 720)

class WorkflowViewerFrame(awx.Frame):
    VERSION = "0.1"

    def __init__(self, parent, workflows, **kwargs):
        """
        Args:
            workflows:
                List of `Workflow` objects or single `Workflow`.
        """
        if "size" not in kwargs:
            kwargs["size"] = _FRAME_SIZE

        super(WorkflowViewerFrame, self).__init__(parent, -1, self.codename, **kwargs)

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

        self.help_menu = wx.Menu()
        self.help_menu.Append(wx.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menuBar.Append(self.help_menu, "Help")

        self.SetMenuBar(menuBar)

        # Create toolbar.
        self.toolbar = toolbar = self.CreateToolBar()

        #tsize = (48,48)
        #artBmp = wx.ArtProvider.GetBitmap
        #toolbar.AddSimpleTool(wx.ID_OPEN, artBmp(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, tsize), "Open")
        toolbar.AddSimpleTool(ID_SHOW_INPUTS, wx.Bitmap(awx.path_img("in.png")), "Visualize the input files of the workflow.")
        toolbar.AddSimpleTool(ID_SHOW_OUTPUTS, wx.Bitmap(awx.path_img("out.png")), "Visualize the output files of the workflow..")
        toolbar.AddSimpleTool(ID_SHOW_LOGS, wx.Bitmap(awx.path_img("log.png")), "Visualize the log files of the workflow.")
        toolbar.AddSimpleTool(ID_BROWSE, wx.Bitmap(awx.path_img("browse.png")), "Browse all the files of the workflow.")
        toolbar.AddSimpleTool(ID_SHOW_MAIN_EVENTS, wx.Bitmap(awx.path_img("out_evt.png")), "Show the ABINIT events reported in the main output files.")
        toolbar.AddSimpleTool(ID_SHOW_LOG_EVENTS, wx.Bitmap(awx.path_img("log_evt.png")), "Show the ABINIT events reported in the log files.")

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
            (ID_BROWSE, self.OnBrowse),
            (ID_SHOW_MAIN_EVENTS, self.OnShowMainEvents),
            (ID_SHOW_LOG_EVENTS, self.OnShowLogEvents),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        self.workflows = workflows
        if not isinstance(workflows, (list, tuple)):
            self.workflows = [workflows]

        #if filename is not None:
        #    self.ReadWorkflow(filename)

        self.BuildUi()

    @property
    def codename(self):
        return self.__class__.__name__

    def BuildUi(self):
        self.panel = panel = wx.Panel(self, -1)
        self.main_sizer = main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Here we create a panel and a notebook on the panel
        self.notebook = Notebook(panel, self.workflows)

        main_sizer.Add(self.notebook, 1, wx.EXPAND, 5)

        check_button = wx.Button(panel, -1, label='Check Status')
        check_button.Bind(wx.EVT_BUTTON, self.OnCheckStatusButton)
        main_sizer.Add(check_button, 0,  wx.ALIGN_CENTER_HORIZONTAL, 5)                                                                                     

        panel.SetSizerAndFit(main_sizer)

        # Register this event when the GUI is IDLE
        # TODO: Use a timer to avoid too much overload.
        #self.Bind(wx.EVT_IDLE, self.OnIdle)

    def GetSelectedWork(self):
        """
        Return the selected workflow namely that the workflow associated to the 
        active tab. None is list is empy.
        """
        return self.notebook.GetSelectedWork()

    #def OnIdle(self, event):
    #    print("On Idle")
    #    self.CheckStatusAndRedraw()

    def OnCheckStatusButton(self, event):
        self.CheckStatusAndRedraw()

    def CheckStatusAndRedraw(self):
        for work in self.workflows:
            work.recheck_status()

        # Save the active tab so that we can set it afterwards.
        old_selection = self.notebook.GetSelection()

        # Redraw the panel
        main_sizer = self.main_sizer
        main_sizer.Hide(0)
        main_sizer.Remove(0)
        new_notebook = Notebook(self, self.workflows)
        main_sizer.Insert(0, new_notebook, 1, wx.EXPAND, 5)
        self.notebook = new_notebook

        self.panel.Layout()

        # Reinstate the old selection
        self.notebook.SetSelection(old_selection)

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

    #def ReadWorkflowFile(self, filepath):
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
    #    self.work = None
    #    self.DestroyPanel()

    #def OnExit(self, event):
    #    self.Close()

    def OnAboutBox(self, event):
        """"Info on the application."""
        awx.makeAboutBox(codename=self.codename, version=self.VERSION,
                         description="", developers="M. Giantomassi")

    def OnShowInputs(self, event):
        work = self.GetSelectedWork() 
        if work is None: return
        frame = TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.abi")
        frame.Show()

    def OnShowOutputs(self, event):
        work = self.GetSelectedWork() 
        if work is None: return
        frame = TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.abo")
        frame.Show()

    def OnShowLogs(self, event):
        work = self.GetSelectedWork() 
        if work is None: return
        frame = TextNotebookFrame.from_files_and_dir(self, dirpath=work.workdir, walk=True, wildcard="*.log")
        frame.Show()

    def OnBrowse(self, event):
        work = self.GetSelectedWork() 
        if work is None: return
        FileListFrame(self, dirpaths=work.workdir).Show()

    def OnShowMainEvents(self, event):
        work = self.GetSelectedWork() 
        if work is None: return
        frame = AbinitEventsNotebookFrame(self, filenames=[task.output_file.path for task in work])
        frame.Show()

    def OnShowLogEvents(self, event):
        work = self.GetSelectedWork() 
        if work is None: return
        frame = AbinitEventsNotebookFrame(self, [task.log_file.path for task in work])
        frame.Show()



class Notebook(fnb.FlatNotebook):
    """
    Notebook class
    """
    def __init__(self, parent, workflows):
        super(Notebook, self).__init__(parent, id=-1, style=fnb.FNB_NO_X_BUTTON | fnb.FNB_NAV_BUTTONS_WHEN_NEEDED)

        self.workflows = workflows
        for work in workflows:
            tab = TabPanel(self, work)
            self.AddPage(tab, text=os.path.basename(work.workdir))

        #self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
        #self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChanging)

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
        if self.GetPageCount() != len(self.workflows):
            return awx.showErrorMessage(self, message="Bad user has removed pages from the notebook!")
                                                                                                       
        idx = self.GetSelection()
        #print("work index is ", idx)
        if idx == -1:
            return None
                                                                                                       
        try:
            return self.workflows[idx]
        except IndexError:
            return None


class TabPanel(wx.Panel):
    """
    Notebook tab.
    """
    def __init__(self, parent, work, **kwargs):
        wx.Panel.__init__(self, parent=parent, id=-1, **kwargs)

        sizer = wx.BoxSizer(wx.VERTICAL)

        # List Control with the individual tasks of the workflow.
        task_listctrl = TaskListCtrl(self, work)
        sizer.Add(task_listctrl, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetSizer(sizer)


class TaskListCtrl(wx.ListCtrl):
    """
    """
    def __init__(self, parent, work, **kwargs):
        """
        Args:
            Task:
                List of `Task` instances.
        """
        super(TaskListCtrl, self).__init__(parent, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN, **kwargs)

        self.work = work

        columns = ["Task", "Status", "Queue_id", "Can run", "Errors", "Warnings", "Comments"]

        for (index, col) in enumerate(columns):
            self.InsertColumn(index, col)

        for task in work:

            events = map(str, 3*["N/A"])
            try:
                report = task.get_event_report()
                if report is not None: 
                    events = map(str, [report.num_errors, report.num_warnings, report.num_comments])
            except:
                pass

            entry = map(str, [task.short_name, task.str_status, task.queue_id, task.can_run] + events)
            self.Append(entry)

        for (index, col) in enumerate(columns):
            self.SetColumnWidth(index, wx.LIST_AUTOSIZE)

        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)
        #self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnItemSelected)

    def OnRightClick(self, event):
        currentItem = event.m_itemIndex
        #print("OnRightClick with currentItem %s" % str(currentItem))
        if currentItem == -1:
            return

        # Open the popup menu then destroy it to avoid mem leak.
        task = self.work[currentItem]
        menu = TaskPopupMenu(parent=self, task=task)
        self.PopupMenu(menu, event.GetPoint())
        menu.Destroy()

    def OnItemActivated(self, event):
        currentItem = event.m_itemIndex
        print("In OnItemActivated with currentItem %s" % str(currentItem))
        task = self.work[currentItem]

        if task.can_run:
            task.start()

    #def OnItemSelected(self, event):
    #    indices = [self.file_list.GetFirstSelected()]
    #    while indices[-1] != -1:
    #        indices.append(self.file_list.GetNextSelected(indices[-1]))
    #    indices = indices[:-1]
    #    print("in OnItemSelected with indices:", indices)

    #    # Plot multiple bands
    #    if len(indices) > 1:
    #        for index in indices:

    #def GetTask(self):
    #    """Returns the Task selected by the user."""
    #    task_idx = int(self.GetFirstSelected())
    #    return self.work[task_idx]


# Callbacks 
_FRAME_SIZE = (720, 720)

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
        message = "Output file %s does not exist" % file.path
        awx.showErrorMessage(parent=parent, message=message)


def show_task_log_events(parent, task):
    file = task.log_file

    if file.exists:
        AbinitEventsFrame(parent, file.path).Show()
    else:
        message = "Log file %s does not exist" % file.path
        awx.showErrorMessage(parent=parent, message=message)


def browse_outdir(parent, task):
    FileListFrame(parent, dirpaths=task.outdata_dir).Show()


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
        ("browse outdir", browse_outdir),
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


#class WorkflowViewerApp(awx.App):
#    def OnInit(self):
#        return True
#
#    def MacOpenFile(self, filename):
#        """Called for files droped on dock icon, or opened via finders context menu"""
#        if filename.endswith(".py"):
#            return
#        # Open filename in a new frame.
#        self.log("%s dropped on app %s" % (filename, self.appname))
#        WorkflowViewerFrame(parent=None, work=filename).Show()


def wxapp_workflow_viewer(workflows):
    app = wx.App()
    frame = WorkflowViewerFrame(None, workflows)
    frame.Show()
    return app
