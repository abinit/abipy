#!/usr/bin/env python
"""MANDELBROT SET"""

import numpy as np

def mandelbrot(extent=None, ndivs=(1000,1000), iterations=100):
    """
    Based on http://cyrille.rossant.net/mandelbrot-set/

    Args:
        extent:
            (x_min, x_max, y_min, y_max))
        (nx,ny):
            Number of points sampled along x and y.
        iterations:
            number of iterations
    """
    if extent is None:
        extent = (-2, 1, -1.5, 1.5)

    x_min, x_max, y_min, y_max = extent
    nx, ny = ndivs

    # x, y are matrices containing the real and imaginary parts of all z values in the grid
    xvals = np.linspace(x_min, x_max, nx)
    yvals = np.linspace(y_min, y_max, ny)

    x, y = np.meshgrid(xvals, yvals)

    # we define c as c=x+iy, c is a nx x ny matrix.
    c = x + 1j*y

    # initially, z=c, we copy so that z and c are different objects in memory
    z = c.copy()

    # m is used to plot the fractal
    m = np.zeros((nx, ny))

    # iterations
    for n in range(iterations):
        print("Completed %d %%" % (100 * n/iterations))

        # indices of the numbers c such that |z(c)|<=10, with z = z_n
        indices = (np.abs(z) <= 10)

        # update z
        z[indices] = z[indices]**2 + c[indices]

        # update the values in m
        m[indices] = n

    return xvals, yvals, m


import wx
from wxmplot import ImageFrame, ImagePanel


class MyPanel(ImagePanel):
    """Pass"""
    def __init__(self, parent, **kwargs):
        super(MyPanel, self).__init__(parent, **kwargs)

    #def onLeftUp(self, event=None):
    #    super(MyPanel, self).onLeftUp(event=event)
    #    print("onLeftUp")
    #    print(self.zoom_ini)

    def zoom_leftup(self, event=None):
        super(MyPanel, self).zoom_leftup(event=event)
        print("zoom_leftup", self.zoom_lims[-1])

    #def redraw(self, col='int'):
    #    print("in redraw")
    #    super(MyPanel, self).redraw(col=col)

    def set_viewlimits(self, axes=None, autoscale=False):
        super(MyPanel, self).set_viewlimits(axes=None, autoscale=False)

        xmin, xmax = self.axes.get_xlim()
        ymin, ymax = self.axes.get_ylim()
        print(xmin, xmax)
        print(ymin, ymax)

        #x, y, m = mandelbrot(density=100)
        #self.clear()
        #self.display(m)

    #def lasso_leftup(self, event=None):
    #    """leftup event handler for lasso mode"""
    #    print("in lasso")

    #def zoom_motion(self, event=None):
    #    print("zoom_motion")
    #    super(MyPanel, self).zoom_motion(event=event)
    #    print(self.zoom_ini)


class FractalsFrame(ImageFrame):
    def __init__(self, parent, **kwargs):
        super(FractalsFrame, self).__init__(parent, **kwargs)

        ndivs = (1000,1000)
        x, y, m = mandelbrot(ndivs=ndvis)
        #print("x",x[0])
        #print("y",x[0])

        self.display(m) # x=x[0], y=x[0])

        self.panel.__class__ = MyPanel

"""
An example of how to use wx or wxagg in an application with the new
toolbar - comment out the setA_toolbar line for no toolbar
"""

# Used to guarantee to use at least Wx2.8
import wxversion
wxversion.ensureMinimal('2.8')



import matplotlib

from mandelbrot import mandelbrot
import matplotlib.pyplot as plt

# uncomment the following to use wx rather than wxagg
#matplotlib.use('WX')
#from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas

# comment out the following to use wx rather than wxagg
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure


class MyNavigationToolbar(NavigationToolbar2Wx):

    def __init__(self, frame):
        super(MyNavigationToolbar, self).__init__(frame.canvas)
        self.frame = frame

    def release_zoom(self, event):
        super(MyNavigationToolbar, self).release_zoom(event)
        print("in release_zoom",event)

        ax = self.canvas.figure.get_axes()[0]
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        print("xlim",xmin, xmax)
        print("ylim",ymin, ymax)

        extent = (xmin, xmax, ymin, ymax)
        self.frame.imshow_fractal(extent=extent)

        self.canvas.draw()
        self.canvas.Refresh(eraseBackground=False)


class ConfPanel(wx.Panel):

    def __init__(self, parent, **kwargs):
        super(ConfPanel, self).__init__(parent, -1, **kwargs)
        self.BuildUi()

    def BuildUi(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Get a list of the colormaps in matplotlib.  Ignore the ones that end with
        # '_r' because these are simply reversed versions of ones that don't end with '_r'
        colormap_list = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
        self.cmap_choice = wx.Choice(self, choices=colormap_list)

        main_sizer.Add(self.cmap_choice, 0, wx.LEFT, 5)

        self.SetSizerAndFit(main_sizer)


class CanvasFrame(wx.Frame):

    def __init__(self, parent):
        wx.Frame.__init__(self, parent, -1, 'CanvasFrame', size=(550,350))

        #main_sizer = wx.BoxSizer(wx.VERTICAL)
        #self.conf_panel = ConfPanel(self)
        #main_sizer.Add(self.conf_panel, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.SetBackgroundColour(wx.NamedColour("WHITE"))
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)

        self.cmap = plt.cm.hot
        self.ndivs = (1000, 1000)

        self.canvas = FigureCanvas(self, -1, self.figure)

        self.sizer = sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)

        self.toolbar = MyNavigationToolbar(self)
        self.toolbar.Realize()

        # On Windows platform, default window size is incorrect, so set toolbar width to figure width.
        tw, th = self.toolbar.GetSizeTuple()
        fw, fh = self.canvas.GetSizeTuple()
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version.
        self.toolbar.SetSize(wx.Size(fw, th))
        sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)

        # Update the axes menu on the toolbar
        self.toolbar.update()

        self.SetSizerAndFit(sizer)

        # Compute fractal and show picture.
        self.imshow_fractal()

    def imshow_fractal(self, extent=None):
        xvals, yvals, data = mandelbrot(extent=extent, ndivs=self.ndivs)

        if extent is None:
            extent = (xvals[0], xvals[-1], yvals[0], yvals[-1])

        self.axes.imshow(np.log(data), cmap=plt.cm.hot, origin="lower", extent=extent)

        # Save values
        self.xvals, self.yvals, self.data = xvals, yvals, data

    #def onclick(self, event):
    #    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
    #          event.button, event.x, event.y, event.xdata, event.ydata))


def wxapp_fractals():
    app = wx.App()
    #frame = FractalsFrame(None)
    frame = CanvasFrame(None)
    app.SetTopWindow(frame)
    frame.Show()
    return app


if __name__ == "__main__":
    import sys
    wxapp_fractals().MainLoop()
    sys.exit(0)

    x, y, m = mandelbrot()

    # we plot log(m)
    import matplotlib.pyplot as plt
    fig = plt.imshow(np.log(m), cmap=plt.cm.hot, extent=extent)

    plt.title('Mandelbrot Set')
    plt.xlabel('Re(z)')
    plt.ylabel('Im(z)')

    plt.show()
