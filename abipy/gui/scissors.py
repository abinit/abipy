#!/usr/bin/env python
from __future__ import print_function, division

import wx

from abipy.tools import AttrDict
from abipy.gui import awx
from abipy.electrons import ScissorsBuilder

class ScissorsParamsDialog(wx.Dialog):
#class ScissorsParamsDialog(wx.Frame):

    def __init__(self, parent, nsppol, num_domains, e0min, e0max, **kwargs):
        assert nsppol == 1

        super(ScissorsParamsDialog, self).__init__(parent, -1, **kwargs)
        self.SetTitle("Select scissors parameters")

        main_sizer = wx.BoxSizer( wx.VERTICAL )

        self.panel_spin = []
        for spin in range(nsppol):
            panel = ScissorsParamsPanel(self, num_domains, e0min, e0max)
            self.panel_spin.append(panel)

            main_sizer.Add(panel, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        hbox = wx.BoxSizer(wx.HORIZONTAL)
                                                                              
        ok_button = wx.Button(self, wx.ID_OK, label='Ok')
        close_button = wx.Button(self, wx.ID_CANCEL, label='Cancel')
                                                                              
        hbox.Add(ok_button)
        hbox.Add(close_button, flag=wx.LEFT, border=5)
                                                                              
        main_sizer.Add(hbox, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(main_sizer)
        self.Fit()

    @awx.verbose
    def GetScissorBuilderParams(self):
        domains_spin, bounds_spin = [], []
        for panel in self.panel_spin:
            p = panel.GetParams()
            domains_spin.append(p.domains)
            bounds_spin.append(p.bounds)

        return AttrDict(
            domains_spin=domains_spin,
            bounds_spin=bounds_spin,
        )


class ScissorsParamsPanel( wx.Panel ):

    def __init__(self, parent, num_domains, e0min, e0max, **kwargs):
        wx.Panel.__init__ ( self, parent, -1, **kwargs)

        self.bound_choices = ["c", "h"]

        self.emin = e0min - 0.1 * abs(e0min)
        self.emax = e0max + 0.1 * abs(e0max)
        self.num_domains = num_domains

        self.BuildUi()

    def BuildUi(self):

        main_sizer = wx.BoxSizer( wx.VERTICAL )

        self._control_domains = []

        for numdom in range(self.num_domains):
            static_boxsizer = self._AddControlBox(numdom)
            main_sizer.Add(static_boxsizer, 1, wx.EXPAND, 5 )

        self.SetSizer( main_sizer )
        self.Fit()

    def _AddControlBox(self, numdom):

        static_boxsizer = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, "Domain #%d" % numdom), wx.VERTICAL )

        vsizer1 = wx.BoxSizer( wx.VERTICAL )

        flexgrid_sizer1 = wx.FlexGridSizer( 0, 4, 0, 0 )
        flexgrid_sizer1.AddGrowableCol( 1 )
        flexgrid_sizer1.SetFlexibleDirection( wx.HORIZONTAL )
        flexgrid_sizer1.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )

        stext1 = wx.StaticText( self, wx.ID_ANY, "Emin:", wx.DefaultPosition, wx.DefaultSize, 0 )
        stext1.Wrap( -1 )
        flexgrid_sizer1.Add( stext1, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.TOP|wx.BOTTOM|wx.LEFT, 5 )

        min_spinctrl = wx.SpinCtrlDouble(self, id=-1, value=str(self.emin), min=self.emin, max=self.emax, inc=0.1)
        flexgrid_sizer1.Add(min_spinctrl, 1, wx.ALL, 5 )

        stext2 = wx.StaticText( self, wx.ID_ANY, "Bound_type:", wx.DefaultPosition, wx.DefaultSize, 0 )
        stext2.Wrap( -1 )
        flexgrid_sizer1.Add(stext2, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.TOP|wx.BOTTOM|wx.LEFT, 5 )

        minbound_choice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, self.bound_choices, 0 )
        minbound_choice.SetSelection( 0 )
        flexgrid_sizer1.Add(minbound_choice, 0, wx.ALL, 5 )

        stext3 = wx.StaticText( self, wx.ID_ANY, "Emax:", wx.DefaultPosition, wx.DefaultSize, 0 )
        stext3.Wrap( -1 )
        flexgrid_sizer1.Add(stext3, 0, wx.TOP|wx.BOTTOM|wx.LEFT|wx.ALIGN_CENTER_HORIZONTAL, 5 )

        max_spinctrl = wx.SpinCtrlDouble(self, id=-1, value=str(self.emax), min=self.emin, max=self.emax, inc=0.1)
        flexgrid_sizer1.Add(max_spinctrl, 1, wx.ALL, 5 )

        stext4 = wx.StaticText( self, wx.ID_ANY, "Bound_type:", wx.DefaultPosition, wx.DefaultSize, 0 )
        stext4.Wrap( -1 )
        flexgrid_sizer1.Add(stext4, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.TOP|wx.BOTTOM|wx.LEFT, 5 )

        maxbound_choice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, self.bound_choices, 0 )
        maxbound_choice.SetSelection( 0 )
        flexgrid_sizer1.Add(maxbound_choice, 0, wx.ALL, 5 )

        vsizer1.Add( flexgrid_sizer1, 1, wx.EXPAND, 5 )

        static_boxsizer.Add( vsizer1, 1, wx.EXPAND, 5 )

        class Control(AttrDict):
            def get_eminmax(self):
                return [float(self.min_spinctrl.GetValue()),
                        float(self.max_spinctrl.GetValue())]

            def get_bound(self):
                return [self.minbound_choice.GetStringSelection(),
                        self.maxbound_choice.GetStringSelection()]

        control = Control(
             min_spinctrl   =min_spinctrl,
             minbound_choice=minbound_choice,
             max_spinctrl   =max_spinctrl,
             maxbound_choice=maxbound_choice,
        )

        self._control_domains.append(control)

        return static_boxsizer

    def GetParams(self):
        domains = [c.get_eminmax() for c in self._control_domains]
        bounds = [c.get_bound() for c in self._control_domains]

        return AttrDict(
            domains=domains,
            bounds=bounds,
        )


class ScissorsBuilderFrame(wx.Frame):

    def __init__(self, parent, filepath=None, **kwargs):
        super(ScissorsBuilderFrame, self).__init__(parent, -1, **kwargs)

        self.scissors_builder = ScissorsBuilder.from_file(filepath)



        #self.OnPlotQps(None)

        self.BuildUi()

    @property
    def nsppol(self):
        return self.scissors_builder.nsppol

    @property
    def e0min(self):
        return self.scissors_builder.e0min
                            
    @property
    def e0max(self):
        return self.scissors_builder.e0max

    def BuildUi(self):
        hsizer = wx.BoxSizer( wx.HORIZONTAL )

        stext = wx.StaticText( self, wx.ID_ANY, "Numer of domains:", wx.DefaultPosition, wx.DefaultSize, 0 )
        stext.Wrap( -1 )
        hsizer.Add(stext, 0, wx.TOP|wx.BOTTOM|wx.LEFT, 5 )

        self.numdomains_ctrl = wx.SpinCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.SP_ARROW_KEYS, 1, 100, 1)
        hsizer.Add( self.numdomains_ctrl, 0, wx.ALL, 5 )

        startfit_button = wx.Button( self, wx.ID_ANY, "Start Fit", wx.DefaultPosition, wx.DefaultSize, 0 )
        startfit_button.Bind(wx.EVT_BUTTON, self.OnStartFitButton)
        hsizer.Add( startfit_button, 0, wx.ALL, 5 )

        self.SetSizer( hsizer )
        self.Fit()

    def OnPlotQps(self, event):
        self.scissors_builder.plot_qpe_vs_e0()

    def OnStartFitButton(self, event):
        num_domains = int(self.numdomains_ctrl.GetValue())

        build_dialog = ScissorsParamsDialog(self, self.nsppol, num_domains, self.e0min, self.e0max)

        if build_dialog.ShowModal() == wx.ID_OK:
            # Get the parameters of the scissors.
            p = build_dialog.GetScissorBuilderParams()
            awx.PRINT("scissor params",p)

            # Try the fit. 
            fitok = True
            try:
                self.scissors_builder.build(p.domains_spin, p.bounds_spin)
            except:
                awx.showErrorMessage(self)
                fitok = False

            # Plot the results and save the current contents in the file.
            if fitok:
                self.scissors_builder.plotfit()

                if awx.askUser(self, "Save data?"):
                    saveFileDialog = wx.FileDialog(self, "Save pickle file", "", "",
                        "Pickle files (*.pickle)|*.pickle", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

                    if saveFileDialog.ShowModal() == wx.ID_OK:
                        filepath = saveFileDialog.GetPath()
                        print("About to save in %s" % filepath)
                        self.scissors_builder.save_data(filepath)

                    saveFileDialog.Destroy()

        build_dialog.Destroy()


def wxapp_scissors(filepath):
    app = wx.App()
    frame = ScissorsBuilderFrame(None, filepath=filepath)
    frame.Show()
    app.MainLoop()
