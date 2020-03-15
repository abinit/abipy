from __future__ import division, print_function

import io
#import wxversion
#wxversion.ensureMinimal('2.8')
import wx
from wx.lib import buttons


__all__ = [
    "DisclosureCtrl"
]


class DisclosureCtrl(buttons.GenBitmapTextToggleButton):
    """Disclosure triangle button."""

    bmp0 = (b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\r\x00\x00'
            b'\x00\r\x08\x06\x00\x00\x00r\xeb\xe4|\x00\x00\x00\x04sBIT'
            b'\x08\x08\x08\x08|\x08d\x88\x00\x00\x00\xd5IDAT(\x91\x9d'
            b'\x921n\x83@\x10E\xdf.!\x96\x9c\x90P"\xd9\x86\x92\x82\x03p'
            b'\x02h\xa9i9A\xe2\xb38)9\x02w\xc0\t\x94I\xa0H\xd2\xc7'
            b'\x96(\x918\xc0\xa6Z\x17\x91\x05&_\x1a\xe9k4O\xfa3\x1a!'
            b'\xa4\xc1\\Im>\xde\xdf\xd4l\xa8,K\x9e\x9fv\xeax\xf8\x99\x84'
            b'\x85\x8e\xb7}|P\x00\xb6m\x13\x04\x01q\x1cc\xdd\xdd\x8bs\xd0'
            b'\xd5\xdfF\xdf\xf7TUE\xdb\xb6\xbc\xbe\xecU\x18\x86\x98\xd7'
            b'\x0b1\ni\r\xc3@Q\x14\xd4u\xcd\xf7\xd7\xa7r]\x97\xe5\xcd'
            b'\xad\x18\x85\xb4\xba\xae#\xcfs|\xdf?\xf5\xe4\xc8<\x00\x8e'
            b'\xe3\x90e\x19i\x9aN\xc7\xb3,\x8b(\x8a\xb8h\xa7Y\xd7'
            b'\xf3<\x0f\xd34I\x92\x84\xd5zsv\xf8$!\r\x844h\x9aFi?Y\xff'
            b'\xf9\xbd_\xd7\x8c7Z\xc0k\x8d8\x00\x00\x00\x00IEND\xaeB`\x82')

    bmp1 = (b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\r\x00\x00'
            b'\x00\r\x08\x06\x00\x00\x00r\xeb\xe4|\x00\x00\x00\x04sBIT'
            b'\x08\x08\x08\x08|\x08d\x88\x00\x00\x00\xc1IDAT(\x91\x9d\xd2;'
            b'\x8e\xc20\x10\x80\xe1\x7fB\x0c\x95\xc1\xa5%\xc0)S\xb8I\xd2'
            b'\xb8Le\x97\xf6\x99s\x81\xdd\xbdBN\xb2\xdb!Y\t\xac`\xbay|'
            b'\xa3)F\xa49\xf0n4o\x8bOQ\x0b\xf0\xf3\xfd\xf5\xbb,\x0b\xeb'
            b'\xba>\x1d\xec\xba\x8ey\x9e\x19\xc6I\x1a\x80a\x9cD)\xf5r'
            b'\xbbR\x8aa\x9c\xa4:\xaf\x94\x821f\x17\x18c(\xa5<\xf2\x07'
            b'\xba\xde\xee\xe2\xbd\xdfE\xde{\xae\xb7\xbbl\x10@J\t\xadu\x05'
            b'\xb4\xd6\xa4\x94\xaaZ\x85\xf4\xf9"1\xc6j \xc6\x88>_\xe4)\x02'
            b'\x08!`\xad\x05\xc0ZK\x08as\xee\x06\xa9\xe3Ir\xce\xb4mK\xce'
            b'\x19u<\xc9\xbf\x08\xc09G\xdf\xf78\xe7\xf6\xda\xc8\'\xbf\xf7'
            b'\x07\x13\x12\x18B\x17\x9fx\xa0\x00\x00\x00\x00IEND\xaeB`\x82')

    def __init__(self, parent, winid, label, *args, **kwds):
        kwds["style"] = wx.BORDER_NONE|wx.BU_EXACTFIT
        buttons.GenBitmapTextToggleButton.__init__(self, parent, winid, None,
                                                   label, *args, **kwds)
        if isinstance(self.bmp0, type(b'')):
            self.__class__.bmp0 = wx.BitmapFromImage(wx.ImageFromStream(
                io.BytesIO(self.bmp0)))
            self.__class__.bmp1 = wx.BitmapFromImage(wx.ImageFromStream(
                io.BytesIO(self.bmp1)))

        self.SetBitmapLabel(self.bmp0)
        self.SetBitmapSelected(self.bmp1)
        if not label:
            self.SetSize(self.bmp0.GetSize())
        else:
            self.SetBestSize()
        self.labelDelta = 0
        self.useFocusInd = False
        self.SetToolTipString('Show')
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)

    def OnEraseBackground(self, event):
        pass

    def Notify(self):
        wx.lib.buttons.GenBitmapTextToggleButton.Notify(self)
        self.SetToolTipString("%s" % ('Show' if self.up else 'Hide'))

    def DoGetBestSize(self):
        width, height, usemin = self._GetLabelSize()
        return width + 5, height + 4

    def OnPaint(self, event):
        width, height = self.GetClientSizeTuple()
        dc = wx.BufferedPaintDC(self)
        bgcol = self.GetBackgroundColour()
        brush = wx.Brush(bgcol, wx.SOLID)
        defattr = self.GetDefaultAttributes()
        if self.style & wx.BORDER_NONE and bgcol == defattr.colBg:
            defattr = self.GetParent().GetDefaultAttributes()
            if self.GetParent().GetBackgroundColour() == defattr.colBg:
                if wx.Platform == "__WXMSW__":
                    if self.DoEraseBackground(dc):
                        brush = None
                elif wx.Platform == "__WXMAC__":
                    brush.MacSetTheme(1)
            else:
                bgcol = self.GetParent().GetBackgroundColour()
                brush = wx.Brush(bgcol, wx.SOLID)
        if brush is not None:
            dc.SetBackground(brush)
            dc.Clear()
        self.DrawLabel(dc, width, height)

    def DrawLabel(self, dc, width, height, center=False):
        bmp = self.bmpLabel
        if bmp is not None:
            if self.bmpDisabled and not self.IsEnabled():
                bmp = self.bmpDisabled
            if self.bmpFocus and self.hasFocus:
                bmp = self.bmpFocus
            if self.bmpSelected and not self.up:
                bmp = self.bmpSelected
            bmpwidth, bmpheight = bmp.GetWidth(), bmp.GetHeight()
            hasmask = bmp.GetMask() is not None
        else:
            bmpwidth = bmpheight = 0

        dc.SetFont(self.GetFont())
        color = (self.GetForegroundColour() if self.IsEnabled() else
                 wx.SystemSettings.GetColour(wx.SYS_COLOUR_GRAYTEXT))
        dc.SetTextForeground(color)

        label = self.GetLabel()
        txtwidth, txtheight = dc.GetTextExtent(label)
        # center bitmap and text
        xpos = (width - bmpwidth - txtwidth) // 2 if center else 0
        if bmp is not None:
            dc.DrawBitmap(bmp, xpos, (height - bmpheight) // 2, hasmask)
            xpos += 5
        dc.DrawText(label, xpos + bmpwidth, (height - txtheight) // 2)
