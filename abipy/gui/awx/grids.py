from __future__ import print_function, division

import sys
import wx

import wx.grid as gridlib

__all__ = [
    "SimpleGrid",
    "SimpleGridFrame",
]

class CopyPasteGrid(gridlib.Grid):
    """ 
    A Copy&Paste enabled grid class
    taken from http://wxpython-users.1045709.n5.nabble.com/copy-and-pasting-selections-in-wx-grid-cells-td2353289.html
    """
    def __init__(self, parent, **kwargs):
        super(CopyPasteGrid, self).__init__(parent, **kwargs)
        wx.EVT_KEY_DOWN(self, self.OnKey)

    def selection(self):
        # Show cell selection
        # If selection is cell...
        if self.GetSelectedCells():
            print("Selected cells " + str(self.GetSelectedCells()))
        # If selection is block...
        if self.GetSelectionBlockTopLeft():
            print("Selection block top left " + str(self.GetSelectionBlockTopLeft()))
        if self.GetSelectionBlockBottomRight():
            print("Selection block bottom right " + str(self.GetSelectionBlockBottomRight()))
        
        # If selection is col...
        if self.GetSelectedCols():
            print("Selected cols " + str(self.GetSelectedCols()))
        
        # If selection is row...
        if self.GetSelectedRows():
            print("Selected rows " + str(self.GetSelectedRows()))
    
    def currentcell(self):
        # Show cursor position
        row = self.GetGridCursorRow()
        col = self.GetGridCursorCol()
        cell = (row, col)
        print("Current cell " + str(cell))
        
    def OnKey(self, event):
        # If Ctrl+C is pressed...
        if event.ControlDown() and event.GetKeyCode() == 67:
            print("Ctrl+C")
            self.selection()
            # Call copy method
            self.copy()
            
        # If Ctrl+V is pressed...
        if event.ControlDown() and event.GetKeyCode() == 86:
            print("Ctrl+V")
            self.currentcell()
            # Call paste method
            self.paste()
            
        # If Supr is pressed
        if event.GetKeyCode() == 127:
            print("Supr")
            # Call delete method
            self.delete()
            
        # Skip other Key events
        if event.GetKeyCode():
            event.Skip()
            return

    def copy(self):
        print("Copy method")
        # Number of rows and cols
        rows = self.GetSelectionBlockBottomRight()[0][0] - self.GetSelectionBlockTopLeft()[0][0] + 1
        cols = self.GetSelectionBlockBottomRight()[0][1] - self.GetSelectionBlockTopLeft()[0][1] + 1
        
        # data variable contain text that must be set in the clipboard
        data = ''
        
        # For each cell in selected range append the cell value in the data variable
        # Tabs '\t' for cols and '\r' for rows
        for r in range(rows):
            for c in range(cols):
                data = data + str(self.GetCellValue(self.GetSelectionBlockTopLeft()[0][0] + r, self.GetSelectionBlockTopLeft()[0][1] + c))
                if c < cols - 1:
                    data = data + '\t'
            data = data + '\n'
        # Create text data object
        clipboard = wx.TextDataObject()
        # Set data object value
        clipboard.SetText(data)
        # Put the data in the clipboard
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()
        else:
            wx.MessageBox("Can't open the clipboard", "Error")
            
    def paste(self):
        print("Paste method")
        clipboard = wx.TextDataObject()
        if wx.TheClipboard.Open():
            wx.TheClipboard.GetData(clipboard)
            wx.TheClipboard.Close()
        else:
            wx.MessageBox("Can't open the clipboard", "Error")
        data = clipboard.GetText()
        table = []
        y = -1
        # Convert text in a array of lines
        for r in data.splitlines():
            y = y +1
            x = -1
            # Convert c in a array of text separated by tab
            for c in r.split('\t'):
                x = x +1
                self.SetCellValue(self.GetGridCursorRow() + y, self.GetGridCursorCol() + x, c)
                
    def delete(self):
        print("Delete method")
        # Number of rows and cols
        rows = self.GetSelectionBlockBottomRight()[0][0] - self.GetSelectionBlockTopLeft()[0][0] + 1
        cols = self.GetSelectionBlockBottomRight()[0][1] - self.GetSelectionBlockTopLeft()[0][1] + 1
        # Clear cells contents
        for r in range(rows):
            for c in range(cols):
                self.SetCellValue(self.GetSelectionBlockTopLeft()[0][0] + r, self.GetSelectionBlockTopLeft()[0][1] + c, '')
   


#class SimpleGrid(gridlib.Grid): 
class SimpleGrid(CopyPasteGrid):

    def __init__(self, parent, table, row_labels=None, col_labels=None, **kwargs):
        """
        Args:
            parent:
                parent window.
            table:
                List of string lists.
            row_labels:
                List of strings used to name the rows.
            col_labels:
                List of strings used to name the col.
        """
        super(SimpleGrid, self).__init__(parent, id=-1, **kwargs)
        self.log = sys.stdout

        self.moveTo = None
        self.Bind(wx.EVT_IDLE, self.OnIdle)
        
        self.nrows = nrows = len(table)
        dims = {len(row) for row in table}
        if len(dims) == 1:
            self.ncols = ncols = list(dims)[0]
        else:
            raise ValueError("Each row must have the same number of columns but dims %s" % str(dims))

        self.CreateGrid(nrows, ncols)
        
        attr = gridlib.GridCellAttr()
        attr.SetTextColour(wx.BLACK)
        attr.SetBackgroundColour(wx.WHITE)
        attr.SetFont(wx.Font(14, wx.SWISS, wx.NORMAL, wx.NORMAL))

        self.SetGridCellAttr(attr)

        if row_labels is not None:
            assert len(row_labels) == nrows
            for i, label in enumerate(row_labels):
                self.SetRowLabelValue(i, label)

        if col_labels is not None:
            assert len(col_labels) == ncols
            for i, label in enumerate(col_labels):
                self.SetColLabelValue(i, label)
            self.SetColLabelAlignment(wx.ALIGN_LEFT, wx.ALIGN_BOTTOM)

        # Cell formatting
        for r, row in enumerate(table):
            for c, col in enumerate(row):
                self.SetCellValue(r, c, table[r][c])
                self.SetReadOnly(r, c, True)

        self.AutoSize()
        self.ForceRefresh()

        # test all the events
        self.Bind(gridlib.EVT_GRID_CELL_LEFT_CLICK, self.OnCellLeftClick)
        self.Bind(gridlib.EVT_GRID_CELL_RIGHT_CLICK, self.OnCellRightClick)
        self.Bind(gridlib.EVT_GRID_CELL_LEFT_DCLICK, self.OnCellLeftDClick)
        self.Bind(gridlib.EVT_GRID_CELL_RIGHT_DCLICK, self.OnCellRightDClick)

        self.Bind(gridlib.EVT_GRID_LABEL_LEFT_CLICK, self.OnLabelLeftClick)
        self.Bind(gridlib.EVT_GRID_LABEL_RIGHT_CLICK, self.OnLabelRightClick)
        self.Bind(gridlib.EVT_GRID_LABEL_LEFT_DCLICK, self.OnLabelLeftDClick)
        self.Bind(gridlib.EVT_GRID_LABEL_RIGHT_DCLICK, self.OnLabelRightDClick)

        self.Bind(gridlib.EVT_GRID_ROW_SIZE, self.OnRowSize)
        self.Bind(gridlib.EVT_GRID_COL_SIZE, self.OnColSize)

        self.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.OnRangeSelect)
        self.Bind(gridlib.EVT_GRID_CELL_CHANGE, self.OnCellChange)
        self.Bind(gridlib.EVT_GRID_SELECT_CELL, self.OnSelectCell)

        self.Bind(gridlib.EVT_GRID_EDITOR_SHOWN, self.OnEditorShown)
        self.Bind(gridlib.EVT_GRID_EDITOR_HIDDEN, self.OnEditorHidden)
        self.Bind(gridlib.EVT_GRID_EDITOR_CREATED, self.OnEditorCreated)

    def SetGridCellAttr(self, attr):
        """
        Set cell attributes for the whole row.

        Args:
            attr:
                `gridlib.GridCellAttr`
        """
        # Note that GridCellAttr objects are reference counted, so attr.IncRef 
        # should be called every time Grid.Set*Attr(attr) is called. This is 
        # required to keep the Grid.Delete* methods from unexpectedly deleting the 
        # GridCellAttr object. 
        for row in range(self.nrows):
            attr.IncRef()
            self.SetRowAttr(row, attr)
        self.AutoSize()
        self.ForceRefresh()

    def OnCellLeftClick(self, evt):
        self.log.write("OnCellLeftClick: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnCellRightClick(self, evt):
        self.log.write("OnCellRightClick: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnCellLeftDClick(self, evt):
        self.log.write("OnCellLeftDClick: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnCellRightDClick(self, evt):
        self.log.write("OnCellRightDClick: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelLeftClick(self, evt):
        self.log.write("OnLabelLeftClick: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelRightClick(self, evt):
        self.log.write("OnLabelRightClick: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelLeftDClick(self, evt):
        self.log.write("OnLabelLeftDClick: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelRightDClick(self, evt):
        self.log.write("OnLabelRightDClick: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnRowSize(self, evt):
        self.log.write("OnRowSize: row %d, %s\n" % (evt.GetRowOrCol(), evt.GetPosition()))
        evt.Skip()

    def OnColSize(self, evt):
        self.log.write("OnColSize: col %d, %s\n" % (evt.GetRowOrCol(), evt.GetPosition()))
        evt.Skip()

    def OnRangeSelect(self, evt):
        if evt.Selecting():
            msg = 'Selected'
        else:
            msg = 'Deselected'
        self.log.write("OnRangeSelect: %s  top-left %s, bottom-right %s\n" %
                           (msg, evt.GetTopLeftCoords(), evt.GetBottomRightCoords()))
        evt.Skip()


    def OnCellChange(self, evt):
        self.log.write("OnCellChange: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))

        # Show how to stay in a cell that has bad data.  We can't just
        # call SetGridCursor here since we are nested inside one so it
        # won't have any effect.  Instead, set coordinates to move to in
        # idle time.
        value = self.GetCellValue(evt.GetRow(), evt.GetCol())

        if value == 'no good':
            self.moveTo = evt.GetRow(), evt.GetCol()

    def OnIdle(self, evt):
        if self.moveTo is not None:
            self.SetGridCursor(self.moveTo[0], self.moveTo[1])
            self.moveTo = None

        evt.Skip()

    def OnSelectCell(self, evt):
        if evt.Selecting():
            msg = 'Selected'
        else:
            msg = 'Deselected'
        self.log.write("OnSelectCell: %s (%d,%d) %s\n" % (msg, evt.GetRow(), evt.GetCol(), evt.GetPosition()))

        # Another way to stay in a cell that has a bad value...
        row = self.GetGridCursorRow()
        col = self.GetGridCursorCol()

        if self.IsCellEditControlEnabled():
            self.HideCellEditControl()
            self.DisableCellEditControl()

        value = self.GetCellValue(row, col)

        if value == 'no good 2':
            return  # cancels the cell selection

        evt.Skip()

    def OnEditorShown(self, evt):
        if evt.GetRow() == 6 and evt.GetCol() == 3 and \
           wx.MessageBox("Are you sure you wish to edit this cell?",
                        "Checking", wx.YES_NO) == wx.NO:
            evt.Veto()
            return

        self.log.write("OnEditorShown: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnEditorHidden(self, evt):
        if evt.GetRow() == 6 and evt.GetCol() == 3 and \
           wx.MessageBox("Are you sure you wish to  finish editing this cell?",
                        "Checking", wx.YES_NO) == wx.NO:
            evt.Veto()
            return

        self.log.write("OnEditorHidden: (%d,%d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnEditorCreated(self, evt):
        self.log.write("OnEditorCreated: (%d, %d) %s\n" % (evt.GetRow(), evt.GetCol(), evt.GetControl()))

    #def InstallGridHint(grid, rowcolhintcallback):
    #    prev_rowcol = [None,None]
    #    def OnMouseMotion(evt):
    #        # evt.GetRow() and evt.GetCol() would be nice to have here,
    #        # but as this is a mouse event, not a grid event, they are not
    #        # available and we need to compute them by hand.
    #        x, y = grid.CalcUnscrolledPosition(evt.GetPosition())
    #        row = grid.YToRow(y)
    #        col = grid.XToCol(x)
    #        if (row,col) != prev_rowcol and row >= 0 and col >= 0:
    #            prev_rowcol[:] = [row,col]
    #            hinttext = rowcolhintcallback(row, col)
    #            if hinttext is None:
    #                hinttext = ''
    #            grid.GetGridWindow().SetToolTipString(hinttext)
    #        evt.Skip()
    #    wx.EVT_MOTION(grid.GetGridWindow(), OnMouseMotion)


class SimpleGridFrame(wx.Frame):
    def __init__(self, parent, table, row_labels=None, col_labels=None, labels_from_table=False, **kwargs):
        """
        Args:
            labels_from_table:
                If True row_labes and col_labels are taken from the table. Not compatible with `row_labels` and `col_labels`.
        """
        super(SimpleGridFrame, self).__init__(parent, -1, **kwargs)

        if labels_from_table:
            # Extract labes from table and extract the subtable.
            assert not row_labels and not col_labels
            col_labels = table[0][1:]
            row_labels, new_table = [], []
            for row in table[1:]:
                row_labels.append(row[0])
                new_table.append(row[1:])
            table = new_table

        self.makeToolBar()

        self.panel = panel = wx.Panel(self, -1)
        self.grid = SimpleGrid(panel, table, row_labels=row_labels, col_labels=col_labels)

        self.main_sizer = main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self.grid, 1, wx.EXPAND | wx.ALL, 5)

        panel.SetSizerAndFit(main_sizer)

    def makeToolBar(self):
        """Creates the tool bar."""
        self.toolbar = toolbar = self.CreateToolBar()

        self.font_picker = wx.FontPickerCtrl(toolbar, -1)
        toolbar.AddControl(control=self.font_picker) 
        toolbar.Realize()

        #menu_bar = wx.MenuBar()
        # Associate menu/toolbar items with their handlers.
        #self.ID_ = wx.NewId()
        #menu_handlers = [
        #    (self.ID_FONT_CHANGED, self.onFontChanged),
        #]
        #for combo in menu_handlers:
        #    mid, handler = combo[:2]
        #    self.Bind(wx.EVT_MENU, handler, id=mid)

        self.Bind(wx.EVT_FONTPICKER_CHANGED, self.OnFontPickerChanged, self.font_picker)

    def OnFontPickerChanged(self, event):
        """Change the Font."""
        font = self.font_picker.GetSelectedFont()
        attr = gridlib.GridCellAttr()
        attr.SetFont(font)
        self.grid.SetGridCellAttr(attr)
        self.main_sizer.Fit(self.panel)


#class DemoSimpleGrid(wx.Frame):
#    def __init__(self, parent):
#        wx.Frame.__init__(self, parent, -1, "Simple Grid Demo")
#
#        row_labels = ["row1", "row2"]
#        col_labels = ["col1", "col2", "col3"]
#
#        table = [ 
#            ["1", "2", "3"],
#            ["4", "5", "6"]
#        ]
#        self.grid = SimpleGrid(self, table, row_labels, col_labels)
#
#
#if __name__ == '__main__':
#    app = wx.App()
#    frame = DemoSimpleGrid(None)
#    frame.Show(True)
#    app.MainLoop()
