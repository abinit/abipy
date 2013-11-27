from __future__ import print_function, division

import abc
import sys
import wx


__all__ = [
    "App",
]


class App(wx.App):
    __metaclass__ = abc.ABCMeta

    def __init__(self, *args, **kwargs):
        wx.App.__init__(self, *args, **kwargs)

        # Initialize the logger 
        import logging
        loglevel = "WARNING"
        loglevel = "DEBUG"
        loglevel = "CRITICAL"

        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % loglevel)
        logging.basicConfig(level=numeric_level)
        
        # This catches events when the app is asked to activate by some other process
        self.Bind(wx.EVT_ACTIVATE_APP, self.OnActivate)

    def OnInit(self): 
        return True

    @property
    def appname(self):
        """Name of the application."""
        return self.__class__.__name__

    def __repr__(self):
        return "<%s at %s>" % (self.appname, id(self))

    def __str__(self):
        return self.__repr__()

    def BringWindowToFront(self):
        try: 
            # it's possible for this event to come when the frame is closed
            self.GetTopWindow().Raise()
        except:
            pass
        
    def OnActivate(self, event):
        # if this is an activate event, rather than something else, like iconize.
        if event.GetActive():
            self.BringWindowToFront()
        event.Skip()
    
    #@abc.abstractmethod
    #def MacOpenFile(self, filename):
    #    """Called for files droped on dock icon, or opened via finders context menu"""
    #    #if filename.endswith(".py"):
    #    #    return
    #    # Code to load filename.
    #    #self.log("%s dropped on app %s" % (filename, self.appname)) 
    #    
    #def MacReopenApp(self):
    #    """Called when the dock icon is clicked."""
    #    self.BringWindowToFront()

    #def MacNewFile(self):
    #    pass

    #def MacPrintFile(self, filepath):
    #    pass
