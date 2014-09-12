#!/usr/bin/env python
from __future__ import division, print_function

import os
import sys
import abc
import six
import numpy as np

from abipy import abilab
from abipy.electrons.gsr import GSR_Reader

from pymatgen.io.abinitio.flows import AbinitFlow

__all__ = [
    "EbandsRobot",
    "EosRobot",
    "FlowInspector",
    "SigresRobot",
]

@six.add_metaclass(abc.ABCMeta)
class FlowRobot(object):
    """
    The main function of a `FlowRobot` is facilitating the extraction of the 
    output data produced by an `AbinitFlow` or the runtime inspections of the 
    calculations. This is the base class from which all Robot classes should 
    derive. It provides helper functions to select the output files of the flow.

    Subclasses should provide the implementation of the two methods:
        
        - collect_data (optional) 
        - analyze_data (mandatory)

    Client code will use robot.run() to produce the results.
    The run method will call collect_data and analyze_data.
    """
    def __init__(self, path):
        """
        Base contructor.

        Args:
            path:
                Filename of the pickle file or directory name.
                If path is a directory, we walk through the directory
                tree starting from path and we read the `AbinitFlow`
                from the first pickle file found.
        """
        self.flow = AbinitFlow.pickle_load(path, disable_signals=True)

    #def walkdirs(self, dir="outdir"):
    #    for task, wi, ti in self.iflat_tasks():
    #        print(wi, ti, task)
    #        #if wi in self.exclude_works or (wi, ti) in self.exclude_tasks: 
    #        #directory = getattr(task, dir)
    #        #if type(task) == abilab.NscfTask:
    #        #    gsr_filepath = directory.has_abiext("GSR")

    def show_flow_status(self, stream=sys.stdout):
        """Print the status of the `AbinitFlow` on stream."""
        return self.flow.show_status(stream=stream)

    def all_outfiles(self, status=None, op="="):
        """
        List of abinit output files produced by the tasks in the flow.
        If status is not None, only the tasks whose status satisfies
        the condition: 

            `task.status op status`

        are selected.
        """
        all_outs = [task.output_file.path for task in 
            self.flow.iflat_tasks(status, op=op)]

        return all_outs

    def all_files_with_abiext(self, abiext, dirname="outdir", 
                              exclude_works=(), exclude_tasks=(),
                              include_works=(), include_tasks=()):
        """
        Returns a list of files with extenxions `abiext` located 
        in the drectory dirname
        """
        files = []
        for task, wi, ti in self.flow.iflat_tasks_wti():

            if wi in exclude_works or (wi, ti) in exclude_tasks: 
                continue

            if ( (include_works and wi not in include_works) or 
                 (include_tasks and (wi, ti) not in include_tasks) ): 
                continue

            directory = getattr(task, dirname)
            path = directory.has_abiext(abiext)
            if path:
                files.append(path)

        return files

    def collect_data(self, *args, **kwargs):
        """
        Collect the data. Subclasses should provide their own implementation (if needed).
        """
                                                     
    @abc.abstractmethod
    def analyze_data(self, *args, **kwargs):
        """
        Analyze the data collected in collect_data.
        """
                                                     
    def run(self, *args, **kwargs):
        """Entry point for client code."""
        self.collect_data(*args, **kwargs)
        return self.analyze_data(*args, **kwargs)


class EbandsRobot(FlowRobot):
    """
    This robot collects the band energies from the GSR files
    and returns an instance of `ElectronBandsPlotter`
    """
    def __init__(self, flow):
        super(EbandsRobot, self).__init__(flow)

        self.plotter = abilab.ElectronBandsPlotter()

    def analyze_data(self, *args, **kwargs):
        for task, in self.flow.iflat_tasks():
            directory = getattr(task, "outdir")

            if type(task) != abilab.NscfTask:
                continue

            gsr_filepath = directory.has_abiext("GSR")
            if gsr_filepath:
                self.plotter.add_ebands_from_file(gsr_filepath, label=gsr_filepath)

        return self.plotter


class EosRobot(FlowRobot):
    """
    This robot computes the equation of state by 
    fitting the total energy as function of the unit cell volume.
    """
    def collect_data(self, *args, **kwargs):
        #self.volumes = [13.72, 14.83, 16.0, 17.23, 18.52]
        #self.energies = [-56.29, -56.41, -56.46, -56.46, -56.42]

        volumes, etotals = [], []
        for task in self.flow.iflat_tasks():
            directory = getattr(task, "outdir")
            if type(task) != abilab.ScfTask:
                continue
                                                                                    
            gsr_filepath = directory.has_abiext("GSR")
            if gsr_filepath:
                with GSR_Reader(gsr_filepath) as r:
                   structure = r.read_structure()
                   # Volumes are in Ang^3
                   volumes.append(structure.volume)

                   # Read energy in Hartree
                   etotals.append(r.read_value("etotal"))

        self.etotals = abilab.ArrayWithUnit(etotals, "Ha").to("eV")
        self.volumes = np.array(volumes)

    def analyze_data(self, *args, **kwargs):
        # Extract volumes and energies from the output files of the calculation.
        # Here we use hardcoded values.

        fits = []
        for eos_name in abilab.EOS.MODELS:
           eos = abilab.EOS(eos_name=eos_name)
           # Note that eos.fit expects lengths in Angstrom, energies are in eV.
           # To specify different units use len_units and ene_units 
           fits.append(eos.fit(self.volumes, self.etotals, vol_unit="ang^3", ene_unit="eV"))

        return fits


class FlowInspector(FlowRobot):
    """
    This robot extract data from the main output files produced by Abinit
    and returns a plottable object.
    """
    def analyze_data(self, *args, **kwargs):
        """Analyze the collected data."""
        # List of SIGRES files computed with different values of nband.
        show = kwargs.get("show", True)

        from pymatgen.io.abinitio.tasks import Task
        from pymatgen.io.abinitio.abiinspect import plottable_from_outfile
        out_files = self.all_outfiles(status=Task.S_RUN, op=">=")

        figs = []
        for out in out_files:
            obj = plottable_from_outfile(out)
            if obj is not None:
                figs.append(obj.plot(title=os.path.relpath(out), show=show))

        return figs


class SigresRobot(FlowRobot):
    """
    This robot collect all the SIGRES files produced by the tasks
    and returns an instance of `SIGRES_Plotter`.
    """
    def collect_data(self, *args, **kwargs):
        """Collect the data"""
        # List of SIGRES files computed with different values of nband.
        self.sigres_files = self.all_files_with_abiext("SIGRES")

    def analyze_data(self, *args, **kwargs):
        """Analyze the collected data."""
        # Instantiate the plotter 
        from abipy.electrons.gw import SIGRES_Plotter
        plotter = SIGRES_Plotter()

        # Add the filepaths to the plotter.
        print(self.sigres_files)
        plotter.add_files(self.sigres_files)

        return plotter
        # Plot the convergence of the QP gaps.
        #plotter.plot_qpgaps(title="QP gaps vs sigma_nband", hspan=0.05)
        # Plot the convergence of the QP energies.
        #plotter.plot_qpenes(title="QP energies vs sigma_nband", hspan=0.05)

