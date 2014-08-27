#!/usr/bin/env python
# -*- coding: utf-8 -*-
# elements_gui.py

# Copyright (c) 2005-2014, Christoph Gohlke
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the copyright holders nor the names of any
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Periodic Table of Elements - A user interface for elements.py.

:Author: `Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>`_

:Version: 2013.03.18

Requirements
------------

* `CPython 2.7 <http://www.python.org>`_
* `wxPython 2.8 <http://www.wxpython.org>`_
* `Elements.py 2013.03.18 <http://www.lfd.uci.edu/~gohlke/>`_

"""
from __future__ import division, print_function

import sys
import math
import io
import webbrowser

import wxversion
#wxversion.ensureMinimal('2.8')
import wx
from wx.lib import fancytext, buttons, rcsizer

from elements import ELEMENTS, SERIES

from abipy.gui.awx.panels import ListCtrlFromTable
try:
    from pseudo_dojo.refdata.nist import database as nist
except ImportError:
    pass


class MainApp(wx.App):
    """Main application."""

    name = "Periodic Table of Elements"
    version = "2012.04.05"
    website = "http://www.lfd.uci.edu/~gohlke/"
    copyright = ("Christoph Gohlke\n"
                 "Laboratory for Fluorescence Dynamics\n"
                 "University of California, Irvine")
    icon = "elements"

    def OnInit(self):
        _ = wx.LogNull()
        wx.InitAllImageHandlers()
        mainframe = WxPeriodicTable(None, -1)
        self.SetTopWindow(mainframe)
        if wx.Platform != "__WXMAC__":
            mainframe.Centre()
        mainframe.ApplyLayout(False)
        mainframe.Show()
        return 1


class ElementButton(buttons.GenToggleButton):
    """Button representing chemical element."""

    def __init__(self, *args, **kwds):
        buttons.GenToggleButton.__init__(self, *args, **kwds)
        self.color = wx.Colour(255, 255, 255)
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        self.Bind(wx.EVT_RIGHT_DOWN, self.onRightDown)

    @property
    def Z(self):
        """Atomic number corresponding to this button."""
        return self.GetId() - 100

    def OnEraseBackground(self, event):
        pass

    def SetButtonColour(self, color):
        self.color = color

    def OnPaint(self, event):
        width, height = self.GetClientSizeTuple()
        dc = wx.BufferedPaintDC(self)
        brush = wx.Brush(self.GetBackgroundColour(), wx.SOLID)
        if wx.Platform == "__WXMAC__":
            brush.MacSetTheme(1)  # kThemeBrushDialogBackgroundActive
        dc.SetBackground(brush)
        dc.SetPen(wx.TRANSPARENT_PEN)
        dc.SetBrush(wx.Brush(self.color if self.up else 'WHITE', wx.SOLID))
        dc.Clear()
        dc.DrawRectanglePointSize((1, 1), (width - 2, height - 2))
        if self.up:
            self.DrawBezel(dc, 0, 0, width - 1, height - 1)
        self.DrawLabel(dc, width, height)
        if (self.hasFocus and self.useFocusInd):
            self.DrawFocusIndicator(dc, width, height)

    def DrawLabel(self, dc, width, height):
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        dc.SetFont(font)

        if self.IsEnabled():
            dc.SetTextForeground(self.GetForegroundColour())
        else:
            dc.SetTextForeground(wx.SystemSettings.GetColour(
                wx.SYS_COLOUR_GRAYTEXT))

        label = self.GetLabel()
        txtwidth, txtheight = dc.GetTextExtent(label)
        xpos = (width - txtwidth) // 2
        ypos = (height*0.75 - txtheight) // 2 - 1
        dc.DrawText(label, xpos, ypos)

        font.SetWeight(wx.FONTWEIGHT_LIGHT)
        font.SetPointSize((font.GetPointSize()*6) // 8)
        dc.SetFont(font)
        label = "%i" % (self.GetId() - 100)
        txtwidth, txtheight = dc.GetTextExtent(label)
        dc.DrawText(label, (width-txtwidth)//2, 4+ypos+(height-txtheight)//2)

    def makePopupMenu(self):
        """Build and return the popup menu."""
        menu = wx.Menu()

        self.ID_POPUP_NIST_LDA = wx.NewId()
        menu.Append(self.ID_POPUP_NIST_LDA, "NIST SCF data (LDA)")

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (self.ID_POPUP_NIST_LDA, self.onNistLda),
        ]
                                                            
        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)
                                                     
        return menu

    def onRightDown(self, event):
        """Called when right button is pressed."""
        popup_menu = self.makePopupMenu()
        self.PopupMenu(popup_menu, event.GetPosition())

    def onNistLda(self, event):
        """
        Show the LDA levels of the neutral atom.
        (useful to decide if semicore states should be included in the valence).
        """
        try:
            entry = nist.get_neutral_entry(self.Z)
        except KeyError:
            return

        table = entry.to_table()
        frame = wx.Frame(self, title="LDA levels for atom %s (NIST database) " % entry.symbol)
        ListCtrlFromTable(frame, table)
        frame.Show()


class PeriodicPanel(wx.Panel):
    """Periodic table of elements panel."""

    layout = """
        .  1  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  18 .
        1  H  2  .  .  .  .  .  .  .  .  .  .  13 14 15 16 17 He .
        2  Li Be .  .  .  .  .  .  .  .  .  .  B  C  N  O  F  Ne .
        3  Na Mg 3  4  5  6  7  8  9  10 11 12 Al Si P  S  Cl Ar .
        4  K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr .
        5  Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe .
        6  Cs Ba *  Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn .
        7  Fr Ra ** Rf Db Sg Bh Hs Mt .  .  .  .  .  .  .  .  .  .
        .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
        .  .  .  *  La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu .
        .  .  .  ** Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr .
    """
    # Class used to instanciates the buttons in the panel.
    element_button_class = ElementButton

    def __init__(self, *args, **kwds):
        kwds["style"] = wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, *args, **kwds)
        rows = len(self.layout.splitlines()) - 2
        cols = len(self.layout.splitlines()[1].split())
        self.sizer = wx.FlexGridSizer(rows, cols, 0, 0)
        self.buttons = list(range(0, len(ELEMENTS)))
        self.selected = -1

        self.info = ElementPanel(self, -1, pos=(0, 0))

        # create element buttons
        buttonsize = math.ceil(float(self.info.GetSize()[0] + 4) / 9.0)
        if buttonsize < 30:
            buttonsize = 30
        for row in self.layout.splitlines()[1:-1]:
            for col in row.split():
                if col == '.':
                    self.sizer.Add((SPACER, SPACER))
                elif col[0] in '123456789*':
                    static = wx.StaticText(self, -1, col,
                        style=wx.ALIGN_CENTER|wx.ALIGN_BOTTOM)
                    self.sizer.Add(static, 0,
                        (wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|
                         wx.ALIGN_CENTER_VERTICAL|wx.FIXED_MINSIZE), SPACER//2)
                else:
                    ele = ELEMENTS[col]
                    button = self.element_button_class(self, ele.number+100, ele.symbol,
                                           size=(buttonsize, buttonsize))
                    self.buttons[ele.number - 1] = button
                    button.SetBezelWidth(1)
                    button.SetToolTipString(ele.name)
                    col = COLORS[ele.series]
                    button.SetButtonColour(wx.Colour(col[0], col[1], col[2]))
                    self.sizer.Add(button, 0, (wx.LEFT|wx.BOTTOM|
                        wx.FIXED_MINSIZE|wx.ALIGN_CENTER_HORIZONTAL|
                        wx.ALIGN_CENTER_VERTICAL), 0)
                    self.Bind(wx.EVT_BUTTON, self.OnSelect, button)

        self.SetAutoLayout(True)
        self.SetSizer(self.sizer)
        self.sizer.SetSizeHints(self)
        self.sizer.Fit(self)

        # position element info panel
        cw = self.sizer.GetColWidths()
        rh = self.sizer.GetRowHeights()
        self.info.Move((sum(cw[:3])+cw[3]//2, (rh[0]-SPACER)//2-1))

        # legend of chemical series
        self.legendpos = (sum(cw[:13]), (rh[0]-SPACER)//2-1)
        self.legendsize = (sum(cw[13:17]) + cw[17]//2, -1)
        self.legend = wx.StaticText(self, -1, " Alkaline earth metals ",
            style=wx.ALIGN_CENTER, pos=self.legendpos, size=self.legendsize)
        self.legend.SetToolTipString("Chemical series")

        # blinking element button
        self.highlight = False
        self.timer = wx.Timer(self)
        self.timer.Start(750)
        self.Bind(wx.EVT_TIMER, self.OnTimer)

    def AddCtrl(self, ctrl, pos=200):
        self.sizer.Remove(pos)
        self.sizer.Insert(pos, ctrl, 0, (wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|
                                         wx.ALIGN_BOTTOM|wx.FIXED_MINSIZE), 0)
        self.Layout()

    def OnTimer(self, evt):
        button = self.buttons[self.selected]
        button.SetToggle(button.up)

    def OnSelect(self, evt):
        self.GetParent().SetSelection(evt.GetId() - 101)

    def SetSelection(self, select):
        """Set active element."""
        if self.selected == select:
            return
        # reset old selection
        self.buttons[self.selected].SetToggle(False)
        # highlight new selection
        self.selected = select
        self.buttons[select].SetToggle(True)
        ele = ELEMENTS[select + 1]
        col = COLORS[ele.series]
        self.legend.SetBackgroundColour(wx.Colour(col[0], col[1], col[2]))
        self.legend.SetLabel(SERIES[ele.series])
        self.legend.Move(self.legendpos)
        self.legend.SetSize(self.legendsize)
        self.info.SetSelection(select)


class ElementPanel(wx.Panel):
    """Element information panel."""

    def __init__(self, *args, **kwds):
        kwds["style"] = wx.NO_BORDER | wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, *args, **kwds)
        self.selected = -1

        # create controls
        self.number = wx.StaticText(self, -1, "808",
            style=wx.ALIGN_RIGHT)
        self.position = wx.StaticText(self, -1, "6, 88, 9",
            style=wx.ALIGN_LEFT)
        self.symbol = wx.StaticText(self, -1, "Mm",
            style=wx.ALIGN_CENTER_HORIZONTAL)
        self.name = wx.StaticText(self, -1, "Praseodymium ",
            style=wx.ALIGN_CENTER_HORIZONTAL)
        self.mass = wx.StaticText(self, -1, "123.4567890 ",
            style=wx.ALIGN_CENTER_HORIZONTAL)
        self.massnumber = wx.StaticText(self, -1, "123 A ",
            style=wx.ALIGN_RIGHT)
        self.protons = wx.StaticText(self, -1, "123 P ",
            style=wx.ALIGN_RIGHT)
        self.neutrons = wx.StaticText(self, -1, "123 N ",
            style=wx.ALIGN_RIGHT)
        self.electrons = wx.StaticText(self, -1, "123 e ",
            style=wx.ALIGN_RIGHT)
        self.eleshell = wx.StaticText(self, -1, "2, 8, 18, 32, 32, 15, 2",
            style=wx.ALIGN_LEFT)
        self.eleconfig = StaticFancyText(self, -1,
            "[Xe] 4f<sup>14</sup> 5d<sup>10</sup>"
            " 6s<sup>2</sup> 6p<sup>6</sup> ",
            style=wx.ALIGN_LEFT)
        self.oxistates = wx.StaticText(self, -1, "1*, 2, 3, 4, 5, 6, -7 ")
        self.atmrad = wx.StaticText(self, -1, "1.234 A ",
            style=wx.ALIGN_RIGHT)
        self.ionpot = wx.StaticText(self, -1, "123.4567890 eV ")
        self.eleneg = wx.StaticText(self, -1, "123.45678 ")

        # set control properties
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.number.SetFont(font)
        self.name.SetFont(font)
        font.SetPointSize(font.GetPointSize() * 1.9)
        self.symbol.SetFont(font)

        self.number.SetToolTipString("Atomic number")
        self.position.SetToolTipString("Group, Period, Block")
        self.symbol.SetToolTipString("Symbol")
        self.name.SetToolTipString("Name")
        self.mass.SetToolTipString("Relative atomic mass")
        self.eleshell.SetToolTipString("Electrons per shell")
        self.massnumber.SetToolTipString("Mass number (most abundant isotope)")
        self.protons.SetToolTipString("Protons")
        self.neutrons.SetToolTipString("Neutrons (most abundant isotope)")
        self.electrons.SetToolTipString("Electrons")
        self.eleconfig.SetToolTipString("Electron configuration")
        self.oxistates.SetToolTipString("Oxidation states")
        self.atmrad.SetToolTipString("Atomic radius")
        self.ionpot.SetToolTipString("Ionization potentials")
        self.eleneg.SetToolTipString("Electronegativity")

        # layout
        sizer = rcsizer.RowColSizer()
        sizer.col_w = SPACER
        sizer.row_h = SPACER
        sizer.Add(self.number, row=0, col=0,
            flag=wx.ALIGN_CENTER_HORIZONTAL|wx.FIXED_MINSIZE)
        sizer.Add(self.position, row=0, col=1,
            flag=wx.ALIGN_LEFT|wx.FIXED_MINSIZE)
        sizer.Add(self.symbol, row=1, col=0, rowspan=2, colspan=2,
            flag=(wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|
                  wx.FIXED_MINSIZE))
        sizer.Add(self.name, row=3, col=0, colspan=2,
            flag=wx.ALIGN_CENTER_HORIZONTAL|wx.FIXED_MINSIZE)
        sizer.Add(self.mass, row=4, col=0, colspan=2,
            flag=wx.ALIGN_CENTER_HORIZONTAL|wx.FIXED_MINSIZE)
        sizer.Add(self.massnumber, row=0, col=2,
            flag=wx.ALIGN_RIGHT|wx.FIXED_MINSIZE)
        sizer.Add(self.protons, row=1, col=2,
            flag=wx.ALIGN_RIGHT|wx.FIXED_MINSIZE)
        sizer.Add(self.neutrons, row=2, col=2,
            flag=wx.ALIGN_RIGHT|wx.FIXED_MINSIZE)
        sizer.Add(self.electrons, row=3, col=2,
            flag=wx.ALIGN_RIGHT|wx.FIXED_MINSIZE)
        sizer.Add(self.atmrad, row=4, col=2,
            flag=wx.ALIGN_RIGHT|wx.FIXED_MINSIZE)
        sizer.Add(self.eleconfig, row=0, col=4, flag=wx.ADJUST_MINSIZE)
        sizer.Add(self.eleshell, row=1, col=4, flag=wx.ADJUST_MINSIZE)
        sizer.Add(self.oxistates, row=2, col=4, flag=wx.ADJUST_MINSIZE)
        sizer.Add(self.ionpot, row=3, col=4, flag=wx.ADJUST_MINSIZE)
        sizer.Add(self.eleneg, row=4, col=4, flag=wx.ADJUST_MINSIZE)

        self.SetAutoLayout(True)
        self.SetSizer(sizer)
        sizer.Fit(self)
        sizer.SetSizeHints(self)
        self.Layout()

    def SetSelection(self, select):
        """Set active element."""
        if self.selected == select:
            return
        self.Freeze()
        self.selected = select
        ele = ELEMENTS[select + 1]

        self.number.SetLabel("%i" % ele.number)
        self.position.SetLabel("%i, %i, %s" % (
            ele.group, ele.period, ele.block))
        self.mass.SetLabel("%.10g" % ele.mass)
        self.eleshell.SetLabel(', '.join("%i" % i for i in ele.eleshells))
        self.massnumber.SetLabel('%i A ' % ele.nominalmass)
        self.protons.SetLabel('%i P ' % ele.protons)
        self.neutrons.SetLabel('%i N ' % ele.neutrons)
        self.electrons.SetLabel('%i e ' % ele.electrons)
        self.oxistates.SetLabel(ele.oxistates)
        self.atmrad.SetLabel(_u("%.10g \xc5 ") % ele.atmrad if ele.atmrad else
                             "")
        self.eleneg.SetLabel("%.10g" % ele.eleneg if ele.eleneg else "")
        self.ionpot.SetLabel(
            "%.10g eV" % ele.ionenergy[0] if ele.ionenergy else "")
        self.symbol.SetLabel(ele.symbol)
        self.name.SetLabel(ele.name)

        label = []
        for orb in ele.eleconfig.split():
            if not orb.startswith('[') and len(orb) > 2:
                orb = orb[:2] + '<sup>' + orb[2:] + '</sup>'
            label.append(orb)
        label.append("<sup> </sup>")  # fix ADJUST_MINSIZE
        self.eleconfig.SetLabel(' '.join(label))

        self.Thaw()
        self.Layout()


class DetailsPanel(wx.Panel):
    """Element details panel."""

    def __init__(self, *args, **kwds):
        kwds["style"] = wx.NO_BORDER|wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, *args, **kwds)
        self.selected = -1

        # create controls
        cb_style = wx.CB_READONLY|wx.CB_SORT
        if wx.Platform == "__WXMAC__":
            cb_style = wx.CB_READONLY


        self.names = LabeledCtrl(self, wx.ComboBox, "Element Name",
            choices=[p.name for p in ELEMENTS],
            style=cb_style, size=(1, -1))
        self.symbols = LabeledCtrl(self, wx.ComboBox, "Symbol", '',
            choices=[p.symbol for p in ELEMENTS],
            style=cb_style, size=(1, -1))
        self.numbers = LabeledCtrl(self, wx.ComboBox, "Number",
            choices=["%s" % p.number for p in ELEMENTS],
            style=wx.CB_READONLY, size=(1, -1))
        self.mass = LabeledCtrl(self, wx.ComboBox, "Relative Atomic Mass",
            choices=["%-.10g" % p.mass for p in ELEMENTS],
            style=wx.CB_READONLY, size=(1, -1))
        self.atmrad = LabeledCtrl(self, wx.ComboBox,
                                  _u("Atomic Radius (\xc5)"),
            choices=["%-.10g" % p.atmrad for p in ELEMENTS],
            style=wx.CB_READONLY, size=(1, -1))
        self.covrad = LabeledCtrl(self, wx.ComboBox,
                                  _u("Covalent Radius (\xc5)"),
            choices=["%-.10g" % p.covrad for p in ELEMENTS],
            style=wx.CB_READONLY, size=(1, -1))
        self.vdwrad = LabeledCtrl(self, wx.ComboBox,
                                  _u("V.d.Waals Radius (\xc5)"),
            choices=["%-.10g" % p.vdwrad for p in ELEMENTS],
            style=wx.CB_READONLY, size=(1, -1))
        self.eleneg = LabeledCtrl(self, wx.ComboBox, "Electronegativity",
            choices=["%-.10g" % p.eleneg for p in ELEMENTS],
            style=wx.CB_READONLY, size=(1, -1))
        self.eleconfig = LabeledCtrl(self, wx.ComboBox, "e- Config",
            choices=[p.eleconfig for p in ELEMENTS],
            style=wx.CB_READONLY, size=(1, -1))
        self.eleshells = LabeledCtrl(self, wx.ComboBox, "Electrons per Shell",
            choices=[', '.join("%i" % i for i in p.eleshells)
                     for p in ELEMENTS],
            style=wx.CB_READONLY, size=(1, -1))
        self.oxistates = LabeledCtrl(self, wx.ComboBox, "Oxidation States",
            choices=[p.oxistates for p in ELEMENTS],
            style=wx.CB_READONLY, size=(1, -1))
        self.ionpot = LabeledCtrl(self, wx.Choice,
            "Ionization Potentials (eV)", choices=[], size=(1, -1))
        self.isotopes = LabeledCtrl(self, wx.Choice, "Isotopes",
            choices=[], size=(1, -1))

        # layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer_top = wx.BoxSizer(wx.HORIZONTAL)
        sizer_left = wx.BoxSizer(wx.VERTICAL)
        sizer_right = wx.BoxSizer(wx.VERTICAL)
        sizer_num = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.RIGHT|wx.BOTTOM|wx.EXPAND|wx.ADJUST_MINSIZE
        sizer_left.Add(self.names, 0, style, SPACER)
        sizer_left.Add(self.mass, 0, style, SPACER)
        sizer_left.Add(self.atmrad, 0, style, SPACER)
        sizer_left.Add(self.covrad, 0, style, SPACER)
        sizer_left.Add(self.vdwrad, 0, style, SPACER)
        sizer_left.Add(self.eleneg, 0, style, SPACER)
        sizer_top.Add(sizer_left, 1, wx.LEFT|wx.RIGHT, 0)
        style = wx.BOTTOM|wx.EXPAND|wx.ADJUST_MINSIZE
        sizer_num.Add(self.symbols, 1, style, 0)
        sizer_num.Add((SPACER, 5), 0, 0, 0)
        sizer_num.Add(self.numbers, 1, style, 0)
        sizer_right.Add(sizer_num, 0, style, SPACER)
        sizer_right.Add(self.eleconfig, 0, style, SPACER)
        sizer_right.Add(self.eleshells, 0, style, SPACER)
        sizer_right.Add(self.oxistates, 0, style, SPACER)
        sizer_right.Add(self.ionpot, 0, style, SPACER)
        sizer_right.Add(self.isotopes, 0, style, SPACER)
        sizer_top.Add(sizer_right, 1, wx.TOP|wx.RIGHT, 0)
        sizer.Add(sizer_top, 1,
            wx.LEFT|wx.RIGHT|wx.TOP|wx.EXPAND|wx.ADJUST_MINSIZE, SPACER)
        self.SetAutoLayout(True)
        self.SetSizerAndFit(sizer, True)
        sizer.SetSizeHints(self)
        self.Layout()

        # bind events
        self.Bind(wx.EVT_COMBOBOX, self.OnSelectName, self.names.ctrl)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelectSymbol, self.symbols.ctrl)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelect, self.numbers.ctrl)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelect, self.mass.ctrl)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelect, self.atmrad.ctrl)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelect, self.covrad.ctrl)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelect, self.vdwrad.ctrl)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelect, self.eleshells.ctrl)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelect, self.eleneg.ctrl)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelect, self.eleconfig.ctrl)
        self.Bind(wx.EVT_COMBOBOX, self.OnSelect, self.oxistates.ctrl)

    def SetSelection(self, select):
        """Set active element."""
        if self.selected == select:
            return
        self.selected = select
        ele = ELEMENTS[select+1]

        self.names.ctrl.SetStringSelection(ele.name)
        self.symbols.ctrl.SetStringSelection(ele.symbol)
        self.numbers.ctrl.SetSelection(select)
        self.mass.ctrl.SetSelection(select)
        self.eleconfig.ctrl.SetSelection(select)
        self.atmrad.ctrl.SetSelection(select)
        self.covrad.ctrl.SetSelection(select)
        self.vdwrad.ctrl.SetSelection(select)
        self.eleneg.ctrl.SetSelection(select)
        self.eleshells.ctrl.SetSelection(select)
        self.oxistates.ctrl.SetSelection(select)

        self.isotopes.ctrl.Clear()
        for index, massnum in enumerate(sorted(ele.isotopes)):
            iso = ele.isotopes[massnum]
            self.isotopes.ctrl.Append("%3i:  %8.4f , %8.4f%%" % (
                massnum, iso.mass, iso.abundance*100.0))
            if massnum == ele.nominalmass:
                self.isotopes.ctrl.SetSelection(index)

        self.ionpot.ctrl.Clear()
        for ion in ele.ionenergy:
            self.ionpot.ctrl.Append("%8.4f" % ion)
        self.ionpot.ctrl.SetSelection(0)

    def OnSelect(self, evt):
        self.SetSelection(evt.GetSelection())
        event = SelectionEvent(pteEVT_ELE_CHANGED, self.GetId(), self.selected)
        self.GetEventHandler().ProcessEvent(event)
        evt.Skip()

    def OnSelectName(self, evt):
        name = self.names.ctrl.GetValue()
        self.SetSelection(ELEMENTS[name].number - 1)
        event = SelectionEvent(pteEVT_ELE_CHANGED, self.GetId(), self.selected)
        self.GetEventHandler().ProcessEvent(event)
        evt.Skip()

    def OnSelectSymbol(self, evt):
        name = self.symbols.ctrl.GetValue()
        self.SetSelection(ELEMENTS[name].number - 1)
        event = SelectionEvent(pteEVT_ELE_CHANGED, self.GetId(), self.selected)
        self.GetEventHandler().ProcessEvent(event)
        evt.Skip()


class DecriptionPanel(wx.Panel):
    """Element description panel."""

    def __init__(self, *args, **kwds):
        kwds["style"] = wx.NO_BORDER|wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, *args, **kwds)
        self.selected = -1

        self.description = wx.TextCtrl(self, -1, " \n \n",
            style=wx.TE_MULTILINE|wx.TE_READONLY)
        font = self.description.GetFont()
        font.SetPointSize((font.GetPointSize() + 1))
        self.description.SetFont(font)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.description, 1,
            wx.TOP|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND|wx.FIXED_MINSIZE,
            SPACER)

        self.SetAutoLayout(True)
        self.SetSizerAndFit(sizer, True)
        sizer.SetSizeHints(self)
        self.Layout()

    def SetSelection(self, select):
        """Set active element."""
        if self.selected == select:
            return
        self.selected = select

        ele = ELEMENTS[select + 1]
        self.description.SetValue(ele.description)


class LabeledCtrl(wx.BoxSizer):
    """BoxSizer containing label, control, and unit."""

    def __init__(self, parent, control, label, unit=None, space="  ",
                 *args, **kwds):
        wx.BoxSizer.__init__(self, wx.HORIZONTAL)
        self.label = wx.StaticText(parent, -1, label + space)
        self.ctrl = control(parent, -1, *args, **kwds)
        self.Add(self.label, 0, wx.ALIGN_CENTER_VERTICAL|wx.FIXED_MINSIZE)
        self.Add(self.ctrl, 1, (wx.LEFT|wx.EXPAND|wx.ALIGN_CENTER_VERTICAL|
                                wx.ALIGN_RIGHT|wx.ADJUST_MINSIZE), 0)
        if unit:
            self.unit = wx.StaticText(parent, -1, unit)
            self.Add(self.unit, 0,
                wx.RIGHT|wx.ALIGN_CENTER_VERTICAL|wx.FIXED_MINSIZE, 0)
        else:
            self.unit = None


class StaticFancyText(fancytext.StaticFancyText):
    """StaticFancyText with SetLabel function."""

    def SetLabel(self, label):
        bmp = fancytext.RenderToBitmap(
            label, wx.Brush(self.GetBackgroundColour(), wx.SOLID))
        self.SetBitmap(bmp)


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
        brush = None
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


class SelectionEvent(wx.PyCommandEvent):
    """Notification of changed element."""

    def __init__(self, evtType, winid, sel):
        wx.PyCommandEvent.__init__(self, evtType, winid)
        self.selection = sel

    def SetSelection(self, select):
        self.selection = select

    def GetSelection(self):
        return self.selection


class WxPeriodicTable(wx.Frame):
    """Main application window."""

    # Class used to instanciate the panel with the elements.
    periodic_panel_class = PeriodicPanel

    def __init__(self, *args, **kwds):
        kwds["style"] = (wx.DEFAULT_DIALOG_STYLE | wx.MINIMIZE_BOX | wx.TAB_TRAVERSAL)
        wx.Frame.__init__(self, *args, **kwds)
        self.selected = -1

        self.SetTitle(MainApp.name)
        icon = wx.EmptyIcon()
        #icon.CopyFromBitmap(wx.Bitmap(MainApp.icon + '.png',
        #                              wx.BITMAP_TYPE_ANY))
        self.SetIcon(icon)
        self.SetBackgroundColour(wx.SystemSettings_GetColour(
            wx.SYS_COLOUR_3DFACE))

        # create menu
        self.menu = wx.MenuBar()
        self.SetMenuBar(self.menu)
        menu = wx.Menu()
        menu.Append(wx.ID_EXIT, "Exit", "Exit the application", wx.ITEM_NORMAL)
        self.menu.Append(menu, "File")
        menu = wx.Menu()
        menu.Append(wx.ID_COPY,
            "Copy\tCtrl+C", "Copy selected element to the clipboard",
            wx.ITEM_NORMAL)
        self.menu.Append(menu, "Edit")
        menu = wx.Menu()
        menu.Append(wx.ID_VIEW_DETAILS, "Details", "Show or hide details",
            wx.ITEM_CHECK)
        self.menu.Append(menu, "View")
        menu = wx.Menu()
        menu.Append(wx.ID_ABOUT, "About...",
            "Display information about the program", wx.ITEM_NORMAL)
        self.menu.Append(menu, "Help")

        # create panels and controls
        self.notebook = wx.Notebook(self, -1, style=0)
        self.description = DecriptionPanel(self.notebook, -1)
        self.details = DetailsPanel(self.notebook, -1)

        self.table = self.periodic_panel_class(self, id=-1)

        self.disclose = DisclosureCtrl(self.table, -1, '')
        self.table.AddCtrl(self.disclose)
        self.notebook.AddPage(self.description, "Description")
        self.notebook.AddPage(self.details, "Properties")

        # event bindings
        self.Bind(wx.EVT_MENU, self.OnExit, id=wx.ID_EXIT)
        self.Bind(wx.EVT_MENU, self.OnCopy, id=wx.ID_COPY)
        self.Bind(wx.EVT_MENU, self.OnDetails, id=wx.ID_VIEW_DETAILS)
        self.Bind(wx.EVT_MENU, self.OnAbout, id=wx.ID_ABOUT)
        self.Bind(EVT_ELE_CHANGED, self.OnSelect, self.table)
        self.Bind(EVT_ELE_CHANGED, self.OnSelect, self.details)
        self.Bind(wx.EVT_BUTTON, self.OnDetails, self.disclose)

        # create sizers
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.table, 1, (wx.LEFT|wx.TOP|wx.RIGHT|wx.EXPAND|
            wx.ALIGN_CENTER_HORIZONTAL| wx.EXPAND|wx.ADJUST_MINSIZE), BORDER-5)
        self.sizer.Add((BORDER, BORDER))
        self.sizer.Add(self.notebook, 0, (wx.LEFT|wx.RIGHT|wx.BOTTOM|
            wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|
            wx.ADJUST_MINSIZE), BORDER)

        self.notebook.SetSelection(1)
        self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer, True)
        self.sizer.SetSizeHints(self)
        self.ApplyLayout(True)
        self.SetSelection(0)

    def ApplyLayout(self, show=False):
        self.show_info = show
        self.sizer.Show(2, show, True)
        self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer, True)
        self.sizer.SetSizeHints(self)
        self.Layout()
        self.menu.Check(wx.ID_VIEW_DETAILS, show)
        self.disclose.SetToggle(self.show_info)

    def OnDetails(self, evt):
        self.ApplyLayout(not self.show_info)

    def OnEraseBackground(self, event):
        pass

    def OnCopy(self, evt):
        dobj = wx.TextDataObject()
        dobj.SetText(repr(ELEMENTS[self.selected + 1]))
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(dobj)
            wx.TheClipboard.Close()

    def OnAbout(self, evt):
        info = wx.AboutDialogInfo()
        info.Name = MainApp.name
        info.Version = MainApp.version
        info.Copyright = MainApp.copyright
        info.WebSite = MainApp.website
        wx.AboutBox(info)

    def OnWebsite(self, evt):
        webbrowser.open(MainApp.website)

    def OnWikipedia(self, evt):
        webbrowser.open("http://en.wikipedia.org/wiki/%s" % (
            ELEMENTS[self.selected].name), 1)

    def OnWebElements(self, evt):
        webbrowser.open("http://www.webelements.com/%s/" % (
            ELEMENTS[self.selected].name.lower()))

    def OnSelect(self, evt):
        self.SetSelection(evt.GetSelection())

    def SetSelection(self, select):
        """Set active element."""
        if self.selected != select:
            self.selected = select
            self.description.SetSelection(select)
            self.table.SetSelection(select)
            self.details.SetSelection(select)

    def OnExit(self, evt):
        self.Close()
        if __name__ == "__main__":
            sys.exit(0)
        else:
            return self.selected


pteEVT_ELE_CHANGED = wx.NewEventType()
EVT_ELE_CHANGED = wx.PyEventBinder(pteEVT_ELE_CHANGED, 1)

SPACER = 12 if wx.Platform == "__WXMAC__" else 10
BORDER = 22 if wx.Platform == "__WXMAC__" else 10

COLORS = {
    1: (0x99, 0xff, 0x99),  # Nonmetals
    2: (0xc0, 0xff, 0xff),  # Noble gases
    3: (0xff, 0x99, 0x99),  # Alkali metals
    4: (0xff, 0xde, 0xad),  # Alkaline earth metals
    5: (0xcc, 0xcc, 0x99),  # Metalloids
    6: (0xff, 0xff, 0x99),  # Halogens
    7: (0xcc, 0xcc, 0xcc),  # Poor metals
    8: (0xff, 0xc0, 0xc0),  # Transition metals
    9: (0xff, 0xbf, 0xff),  # Lanthanides
    10: (0xff, 0x99, 0xcc),  # Actinides
}

_u = (lambda x: x.decode('latin-1')) if sys.version[0] == '2' else str

if __name__ == "__main__":
    MainApp(0).MainLoop()
