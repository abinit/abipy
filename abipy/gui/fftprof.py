#!/usr/bin/env python
from __future__ import print_function, division

import wx

from collections import OrderedDict
from abipy.gui import awx
from abipy.abitools.fftprof import FFTProf

FFTALGS = OrderedDict([
    (112, "Goedecker"),
    (312, "FFTW3"),
    (412, "Goedecker-2002"),
    (512, "MKL-DFTI"),
])


class FftalgsPanel(wx.Panel):

    def __init__( self, parent ):
        wx.Panel.__init__(self, parent, id=-1)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        static_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "FFT Libraries"), wx.VERTICAL)
        grid_sizer = wx.GridSizer( 0, 2, 0, 0 )

        self.check_boxes = OrderedDict()

        for fftalg, fft_name in FFTALGS.items():
            cbox = wx.CheckBox(self, -1, fft_name)
            if fftalg == 112:
                cbox.SetValue(True)
            else:
                cbox.SetValue(False)

            self.check_boxes[fftalg] = cbox
            grid_sizer.Add(cbox, 0, wx.ALL, 5 )

        static_sizer.Add(grid_sizer, 0, wx.ALL | wx.EXPAND, 5)
        main_sizer.Add(static_sizer, 1, wx.EXPAND, 5)
        self.SetSizerAndFit(main_sizer)

    def GetValues(self):
        """Returns a list with the fftalgs selected by the user."""
        fftalgs = []
        for (fftalg, cbox) in self.check_boxes.items():
            if cbox.GetValue():
                fftalgs.append(fftalg)

        if not fftalgs: 
            fftalgs = [112]

        return fftalgs


class ControlPanel(wx.Panel):
    """
    !!   &CONTROL
    !!     tasks = string specifying the tasks to perform i.e. the routines that should be tested or profiled. 
    !!             allowed values: 
    !!                 fourdp --> Test FFT transforms of density and potentials on the full box.
    !!                 fourwf --> Test FFT transforms of wavefunctions using the zero-padded algorithm.
    !!                 gw_fft --> Test the FFT transforms used in the GW code.
    !!                 all    --> Test all FFT routines (DEFAULT)
    !!     fftalgs = list of fftalg values (used to select the FFT libraries to use, see abinit doc for more info)
    !!     ncalls = integer defining the number of calls for each tests. The final Wall time and CPU time
    !!              are computed by averaging the final results over ncalls executions.
    !!     max_nthreads = Maximum number of threads (DEFAULT 1, meaningful only if the code
    !!                    uses threaded external libraries or OpenMP parallelization)
    !!     ndat   = integer specifying how many FFT transforms should be executed for each call to the FFT routine
    !!              (same meaning as the ndat input variable passed to fourwf)
    !!     necut  = Used if tasks = "bench". Specifies the number of cutoff energies to profile (see also ecut_arth) 
    !!     ecut_arth = Used if tasks = "bench". Used to generate an arithmetic progression of cutoff energies 
    !!                 whose starting value is ecut_arth(1) and whose step is ecut_arth(2)
    """
    BENCHMARKS = ["bench_fourwf", "bench_fourdp", "bench_rhotwg"]

    def __init__( self, parent ):
        wx.Panel.__init__ ( self, parent, id =-1)

        main_sizer = wx.BoxSizer( wx.VERTICAL )

        bSizer5 = wx.BoxSizer( wx.HORIZONTAL )
                                                                                                        
        text = wx.StaticText( self, -1, "Tasks:", wx.DefaultPosition, wx.DefaultSize, 0 )
        text.Wrap( -1 )
        text.SetToolTipString("Select the routines for the benchmark.")
        self.task_choice = wx.Choice( self, -1, wx.DefaultPosition, wx.DefaultSize, self.BENCHMARKS, 0 )
        self.task_choice.SetSelection( 0 )

        bSizer5.Add(text, 0, wx.TOP|wx.BOTTOM|wx.LEFT, 5 )
        bSizer5.Add( self.task_choice, 0, wx.ALL, 5 )
                                                                                                        
        main_sizer.Add( bSizer5, 1, wx.EXPAND, 5 )

        hsz1 = wx.BoxSizer( wx.HORIZONTAL )

        text = wx.StaticText( self, -1, "ndat:", wx.DefaultPosition, wx.DefaultSize, 0 )
        text.Wrap( -1 )
        text.SetToolTipString("Number of FFT transform per call.")
        self.ndat_ctrl = wx.SpinCtrl( self, -1, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.SP_ARROW_KEYS, 0, 10, 1 )

        hsz1.Add(text, 0, wx.TOP|wx.BOTTOM|wx.LEFT, 5 )
        hsz1.Add( self.ndat_ctrl, 0, wx.ALL, 5 )

        text = wx.StaticText( self, -1, "ncalls:", wx.DefaultPosition, wx.DefaultSize, 0 )
        text.Wrap( -1 )
        text.SetToolTipString("Number of calls (used to get mean results)")
        self.ncalls_ctrl = wx.SpinCtrl( self, -1, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.SP_ARROW_KEYS, 0, 10, 5 )

        hsz1.Add(text, 0, wx.ALL, 5 )
        hsz1.Add( self.ncalls_ctrl, 0, wx.ALL, 5 )

        text = wx.StaticText( self, -1, "max_nthreads:", wx.DefaultPosition, wx.DefaultSize, 0 )
        text.Wrap( -1 )
        text.SetToolTipString("Maximum number of OpenMP threads")
        self.max_nthreads_ctrl = wx.SpinCtrl( self, -1, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.SP_ARROW_KEYS, 0, 10, 1 )

        hsz1.Add(text, 0, wx.ALL, 5 )
        hsz1.Add(self.max_nthreads_ctrl, 0, wx.ALL, 5 )

        #text = wx.StaticText( self, -1, "ecut_range:", wx.DefaultPosition, wx.DefaultSize, 0 )
        #text.Wrap( -1 )
        #hsz1.Add(text, 0, wx.ALL, 5 )
        #                                                                                                                           
        #self.ecutrange_ctrl = wx.SpinCtrl( self, -1, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.SP_ARROW_KEYS, 0, 10, 0 )
        #hsz1.Add(self.ecutrange_ctrl, 0, wx.ALL, 5 )

        main_sizer.Add(hsz1, 1, wx.EXPAND, 5 )

        self.fftalgs = FftalgsPanel(self)
        main_sizer.Add(self.fftalgs, 1, wx.EXPAND, 5 )

        self.SetSizerAndFit(main_sizer)

    def GetNamelist(self):
        """
        Returns the Fortran NAMELIST.
        """
        d = OrderedDict(
            tasks="'" + self.task_choice.GetStringSelection() + "'",
            fftalgs=",".join(str(v) for v in self.fftalgs.GetValues()),
            ncalls=self.ncalls_ctrl.GetValue()),
            max_nthreads=self.max_nthreads_ctrl.GetValue(),
            ndat=self.ndat_ctrl.GetValue(),
            necut=2,
            ecut_arth = "10, 5",
        )

        namelist = ["&CONTROL"]
        for k, v in d.items():
            namelist.append(k + " = " + str(v) + ",")
        namelist.append("/")

        return "\n".join(namelist)


class KpointChoice(wx.Choice):
    def __init__(self, parent):
        choices = map(str, [
            (0.1, 0.2, 0.3),
            (0.0, 0.0, 0.0),
            (1/2, 0 , 0 ),
            ( 0 , 0 ,1/2),
            (1/2, 0 ,1/2),
            ( 0 ,1/2, 0 ),
            (1/2,1/2, 0 ),
            ( 0 ,1/2,1/2),
            (1/2,1/2,1/2),
        ])

        wx.Choice.__init__(self, parent, -1, choices=choices)
        self.SetSelection(0)

    def GetStringSelection(self):
        s = super(KpointChoice, self).GetStringSelection()
        return s.replace("(", "").replace(")", "")


class SystemPanel(wx.Panel):
    """
    !!     ecut = cutoff energy for wavefunctions (real, Hartree units)
    !!     rprimd = Direct lattice vectors in Bohr. (3,3) matrix in Fortran column-major order
    !!     kpoint = real(3) vector specifying the reduced coordinates of the k-point of the wavefunction (used if tasks = "fourwf"). 
    !!               The value of the k-point defines the storage scheme (istwfk) of the u(G) coefficients and therefore
    !!               the FFT algorithm used to perform the transform u(G) <--> u(R) in fourwf.
    !!     osc_ecut = cutoff energy (Hartree) for the oscillator matrix elements computed in the GW code
    !!                 Corresponds to the input variables (ecuteps, ecutsigx) used in the main code.
    !!     nsym     =Number of symmetry operations (DEFAULT 1)
    !!     symrel(3,3,nsym) = Symmetry operation in real space used to select the FFT mesh in the routine getng (default: Identity matrix)
    """
    def __init__(self, parent):
        wx.Panel.__init__ ( self, parent, id =-1)

        main_sizer = wx.BoxSizer( wx.VERTICAL )

        hsz1 = wx.BoxSizer( wx.HORIZONTAL )

        text = wx.StaticText( self, -1, "kpoint:", wx.DefaultPosition, wx.DefaultSize, 0 )
        text.Wrap( -1 )
        text.SetToolTipString("K-point used in the FFT of the wavefunctions (istwfk).")
        hsz1.Add(text, 0, wx.TOP|wx.BOTTOM|wx.LEFT, 5 )
        self.kpoint_choice = KpointChoice(self)

        hsz1.Add( self.kpoint_choice, 0, wx.ALL, 5 )

        #text = wx.StaticText( self, -1, "ndat:", wx.DefaultPosition, wx.DefaultSize, 0 )
        #text.Wrap( -1 )
        #hsz1.Add(text, 0, wx.TOP|wx.BOTTOM|wx.LEFT, 5 )
        #                                                                                                                         
        #self.ndat_ctrl = wx.SpinCtrl( self, -1, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.SP_ARROW_KEYS, 0, 10, 4 )
        #hsz1.Add( self.ndat_ctrl, 0, wx.ALL, 5 )

        main_sizer.Add(hsz1, 1, wx.EXPAND, 5 )

        self.SetSizerAndFit(main_sizer)

    def GetNamelist(self):
        """
        Returns the Fortran NAMELIST.
        """
        d = OrderedDict(
            ecut=10,
            rprimd= "20, 0, 0, 0, 20, 0, 0, 0, 20",
            kpoint=self.kpoint_choice.GetStringSelection()
            #osc_ecut= TODO why here and not in control?
            #nsym=
            #symrel=
        )

        namelist = ["&SYSTEM"]
        for k, v in d.items():
            namelist.append(k + " = " + str(v) + ",")
        namelist.append("/")

        return "\n".join(namelist)


class FFTProfFrame(awx.Frame):
    def __init__(self, parent, **kwargs):
        super(FFTProfFrame, self).__init__(parent, -1, **kwargs)

        panel = wx.Panel(self, -1)

        self.control = ControlPanel(panel)
        self.system = SystemPanel(panel)

        start_button = wx.Button(panel, -1, label='Start Benchmark')
        start_button.Bind(wx.EVT_BUTTON, self.OnStartButton)

        main_sizer = wx.BoxSizer( wx.VERTICAL )
        main_sizer.Add(self.control, 1, wx.EXPAND, 5 )
        main_sizer.Add(self.system, 1, wx.EXPAND, 5 )
        main_sizer.Add(start_button,  0,  wx.ALIGN_CENTER_HORIZONTAL, 5)

        panel.SetSizerAndFit(main_sizer)

    def OnStartButton(self, event):
        """Run the benchmark."""
        # Build the input file for FFTPROF.
        fft_input  = self.control.GetNamelist() + "\n"
        fft_input += self.system.GetNamelist() + "\n"
        print(fft_input)

        # Run FFTPROF.
        try:
            fftprof = FFTProf(fft_input)
            returncode = fftprof.run()
            #if returncode != 0:

            # Analyze the results.
            fftprof.plot()

        except:
            return awx.showErrorMessage(self)


def wxapp_fftprof():
    """Stand alone application."""
    app = wx.App()
    frame = FFTProfFrame(None)
    frame.Show()
    return app


if __name__ == "__main__":
    wxapp_fftprof().MainLoop()
