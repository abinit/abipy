import wx

from collections import OrderedDict
from abipy.gui import awx

from abipy.tools.fftprof import FFTProf
from abipy.gui.editor import SimpleTextViewer
from .awx.panels import LinspaceControl


FFTALGS = OrderedDict([
    (112, "Goedecker"),
    (312, "FFTW3"),
    (412, "Goedecker-2002"),
    (512, "MKL-DFTI"),
])


class FftalgsPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, id=-1)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        static_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "FFT Libraries"), wx.VERTICAL)
        grid_sizer = wx.GridSizer(0, 2, 0, 0)

        self.check_boxes = OrderedDict()

        for fftalg, fft_name in FFTALGS.items():
            cbox = wx.CheckBox(self, -1, fft_name)
            if fftalg == 112:
                cbox.SetValue(True)
            else:
                cbox.SetValue(False)

            self.check_boxes[fftalg] = cbox
            grid_sizer.Add(cbox, 0, wx.ALL, 5)

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

    def __init__(self, parent):
        wx.Panel.__init__(self, parent, id=-1)

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsz1 = wx.BoxSizer(wx.HORIZONTAL)

        text = wx.StaticText(self, -1, "Tasks:")
        text.Wrap(-1)
        text.SetToolTipString("Select the routines for the benchmark.")
        self.task_choice = wx.Choice(self, -1, wx.DefaultPosition, wx.DefaultSize, self.BENCHMARKS, 0)
        self.task_choice.SetSelection(0)

        hsz1.Add(text, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        hsz1.Add(self.task_choice, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        text = wx.StaticText(self, -1, "ndat:")
        text.Wrap(-1)
        text.SetToolTipString("Number of FFT transform per call.")
        self.ndat_ctrl = wx.SpinCtrl(self, -1, value="1", min=1)

        hsz1.Add(text, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        hsz1.Add(self.ndat_ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        text = wx.StaticText(self, -1, "ncalls:")
        text.Wrap(-1)
        text.SetToolTipString("Number of calls (used to get mean results)")
        self.ncalls_ctrl = wx.SpinCtrl(self, -1, value="5", min=1)

        hsz1.Add(text, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        hsz1.Add(self.ncalls_ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        text = wx.StaticText(self, -1, "max_nthreads:")
        text.Wrap(-1)
        text.SetToolTipString("Maximum number of OpenMP threads")
        self.max_nthreads_ctrl = wx.SpinCtrl(self, -1, value="1", min=1)

        hsz1.Add(text, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        hsz1.Add(self.max_nthreads_ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        main_sizer.Add(hsz1, 1, wx.EXPAND, 5)

        # Control to read the range of cutoff energies.
        self.ecut_linspace = LinspaceControl(self, start=dict(value="20"), stop=dict(value="120"), num=dict(value="10"))

        static_sizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Ecut list:"), wx.VERTICAL)
        static_sizer.Add(self.ecut_linspace, 0, wx.ALL | wx.EXPAND, 5)

        # Control to specify fftalgs.
        self.fftalgs = FftalgsPanel(self)

        #vsz1 = wx.BoxSizer(wx.VERTICAL)
        #vsz1.Add(static_sizer, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        #vsz1.Add(self.fftalgs, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        #main_sizer.Add(vsz1, 1, wx.EXPAND, 5)

        main_sizer.Add(static_sizer, 1, wx.EXPAND, 5)
        main_sizer.Add(self.fftalgs, 1, wx.EXPAND, 5)

        self.SetSizerAndFit(main_sizer)

    def GetNamelist(self):
        """
        Build and returns a string with whe &CONTROL Fortran NAMELIST.
        """
        ecut_list = self.ecut_linspace.getValues()
        if len(ecut_list) <= 1:
            return awx.showErrorMessage(self, message="Ecut list is empty or contains only one item!")
        step = ecut_list[1] - ecut_list[0]

        d = OrderedDict(
            tasks="'" + self.task_choice.GetStringSelection() + "'",
            fftalgs=",".join(str(v) for v in self.fftalgs.GetValues()),
            ncalls=self.ncalls_ctrl.GetValue(),
            max_nthreads=self.max_nthreads_ctrl.GetValue(),
            ndat=self.ndat_ctrl.GetValue(),
            necut=len(ecut_list),
            ecut_arth="%s, %s" % (ecut_list[0], step)
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
            (1 / 2, 0, 0),
            (0, 0, 1 / 2),
            (1 / 2, 0, 1 / 2),
            (0, 1 / 2, 0),
            (1 / 2, 1 / 2, 0),
            (0, 1 / 2, 1 / 2),
            (1 / 2, 1 / 2, 1 / 2),
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
        wx.Panel.__init__(self, parent, id=-1)

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        hsz1 = wx.BoxSizer(wx.HORIZONTAL)

        text = wx.StaticText(self, -1, "kpoint:")
        text.Wrap(-1)
        text.SetToolTipString("K-point used in the FFT of the wavefunctions (istwfk).")
        self.kpoint_choice = KpointChoice(self)

        hsz1.Add(text, 0, wx.ALIGN_CENTER_VERTICAL | wx.TOP | wx.BOTTOM | wx.LEFT, 5)
        hsz1.Add(self.kpoint_choice, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        main_sizer.Add(hsz1, 1, wx.EXPAND, 5)

        self.SetSizerAndFit(main_sizer)

    def GetNamelist(self):
        """
        Build and returns a string with whe &SYSTEM Fortran NAMELIST.
        """
        d = OrderedDict(
            ecut=10,
            rprimd="20, 0, 0, 0, 20, 0, 0, 0, 20",
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

        show_button = wx.Button(panel, -1, label='Show Input')
        show_button.Bind(wx.EVT_BUTTON, self.OnShowInput)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self.control, 1, wx.EXPAND, 5)
        main_sizer.Add(self.system, 1, wx.EXPAND, 5)

        hsz = wx.BoxSizer(wx.HORIZONTAL)
        hsz.Add(start_button, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)
        hsz.Add(show_button, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        main_sizer.Add(hsz, 0, wx.ALIGN_CENTER_HORIZONTAL, 5)

        panel.SetSizerAndFit(main_sizer)

    def MakeInput(self):
        """Build the input file for FFTPROF."""
        fft_input = self.control.GetNamelist() + "\n"
        fft_input += self.system.GetNamelist() + "\n"
        return fft_input

    def OnShowInput(self, event):
        """Show the fftprof input in an external window."""
        fft_input = self.MakeInput()
        SimpleTextViewer(self, fft_input).Show()

    def OnStartButton(self, event):
        """Run the benchmark and plot the results with `matplotlib`."""
        fft_input = self.MakeInput()
        #return

        # Run FFTPROF.
        try:
            fftprof = FFTProf(fft_input)
            returncode = fftprof.run()

            if returncode != 0:
                raise RunTimeError("fftprof exited with returncode %s" % returncode)

            # Analyze the results.
            fftprof.plot()

        except:
            return awx.showErrorMessage(self)


def wxapp_fftprof():
    """Stand alone application."""
    app = awx.App()
    frame = FFTProfFrame(None)
    frame.Show()
    return app


if __name__ == "__main__":
    wxapp_fftprof().MainLoop()
