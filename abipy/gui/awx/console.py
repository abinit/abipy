from __future__ import print_function, division

import wx

__all__ = [
    "ConsoleEvent",
    "ConsoleWindow",
]

#######################################################################
# Events for the console widget
#######################################################################

EVT_CONSOLE_TYPE = wx.NewEventType()
EVT_CONSOLE = wx.PyEventBinder(EVT_CONSOLE_TYPE, 1)


class ConsoleEvent(wx.PyEvent):
    """
    This event is triggering for messages that
    should be displayed within the console.
    """
    def __init__(self, kind, msg):
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_CONSOLE_TYPE)
        self.kind, self.msg = kind, msg


class ConsoleWindow(wx.TextCtrl):
    """
    This serves as the primary console for displaying commands and initialization information.
    The console window is color coded. To write to the console window, widgets should post console events.
    """
    # Message kinds
    INFO = "info"
    DEBUG = "debug"
    WARNING = "warning"
    CRITICAL = "critical"
    ERROR = 'error'

    ALL_KINDS =[INFO, DEBUG, WARNING, CRITICAL, ERROR]

    def __init__(self, parent):
        """Creates the console window."""
        wx.TextCtrl.__init__(self, parent, -1, style=wx.TE_MULTILINE | wx.TE_READONLY | wx.NO_BORDER | wx.TE_RICH2)

        self.idPOPUP_CLEAR = wx.NewId()

        # "kind > %s\n",
        self.format = {}
        for kind in self.ALL_KINDS:
            self.format[kind] = kind + "> %s\n"

        self.colors = {k: None for k in self.ALL_KINDS}
        self.wxColors = {k: None for k in self.ALL_KINDS}

        # Set ourselves to receive all console events
        self.Bind(EVT_CONSOLE, self.onConsoleEvent)

        # Set up the console menu handler
        self.Bind(wx.EVT_RIGHT_DOWN, self.onRightClick, self)

    def displayBanner(self, **kwargs):
        """Displays the application startup banner."""
        self.insertText('===============================================')
        self.insertText("Welcome to program %s" % version)
        self.insertText('===============================================')
        self.insertText('')

    def insertText(self, text, kind=None):
        """
        Writes text to the console window in the appropriate color for the kind of message.
        """
        self.SetDefaultStyle(wx.TextAttr(self.getColor(kind)))
        self.AppendText(self.format[kind] % text)

    def getColor(self, kind):
        """Gets the wx color object for the specified text kind."""
        #config = utils.getApplicationConfiguration()
        #current_color = config[kind + '_color']
        if kind is None:
            return None

        wx_color = None

        #current_color = "error"
        #if self.colors[kind] == current_color:
        #    wx_color = self.wxColors[kind]
        #else:
        #    wx_color = apply(wx.Colour, current_color)
        #    self.colors[kind] = current_color
        #    self.wxColors[kind] = wx_color

        return wx_color

    def onConsoleEvent(self, event):
        """
        Processes a console event by displaying the contents of the event to the console window.
        """
        self.insertText(event.msg, event.kind)

    def onRightClick(self, event):
        """Handles a right click by displaying the popupmenu at the click point."""
        menu = self.makePopupMenu()
        self.PopupMenu(menu, event.GetPosition())

    def makePopupMenu(self):
        """Creates the popup menu for operations on the console window."""
        menu = wx.Menu()
        menu.Append(self.idPOPUP_CLEAR, "Clear Console")

        self.Bind(wx.EVT_MENU, self.onClear, id=self.idPOPUP_CLEAR)
        return menu

    def onClear(self, event):
        """Processes a clearing of the console window event."""
        self.Clear()


class FrameWithConsole(wx.Frame):
    """
    Simple frame with a Choice control and two buttons: OK, Cancel
    Subclasses will usually redefine onOkButton.
    """
    def __init__(self, parent, choices, **kwargs):
        super(FrameWithConsole, self).__init__(parent, id=-1, **kwargs)
        self.parent = parent

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        panel = wx.Panel(self, -1)
        main_sizer.Add(panel, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL | wx.EXPAND, border=5)

        ok_button = wx.Button(self, wx.ID_OK, label='Ok')
        ok_button.Bind(wx.EVT_BUTTON, self.onOkButton)

        close_button = wx.Button(self, wx.ID_CANCEL, label='Cancel')
        close_button.Bind(wx.EVT_BUTTON, self.onCloseButton)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(ok_button)
        hbox.Add(close_button, flag=wx.LEFT, border=5)

        main_sizer.Add(hbox, flag=wx.ALIGN_CENTER | wx.ALL | wx.EXPAND)

        splitter = wx.SplitterWindow(self, -1, style=wx.SP_3DSASH)
        self.console = ConsoleWindow(splitter)
        csplitter = wx.SplitterWindow(splitter, -1, style=wx.SP_3DSASH)
        # Set the default control window size
        splitter.SplitHorizontally(self.console, csplitter, 150)

        main_sizer.Add(splitter, 1, wx.EXPAND)
        #main_sizer.Add(self.console, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL | wx.EXPAND, border=5)

        self.SetAutoLayout(True)
        self.SetSizerAndFit(main_sizer)

    def onOkButton(self, event):
        event = ConsoleEvent(kind="error", msg="Hello")
        #print("Firing %s" % event)
        wx.PostEvent(self.console, event)

    def onCloseButton(self, event):
        self.Destroy()


if __name__ == "__main__":
    app = wx.App()
    frame = FrameWithConsole(None, -1)
    frame.Show(True)
    app.SetTopWindow(frame)
    app.MainLoop()
