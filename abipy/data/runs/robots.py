#!/usr/bin/env python
from __future__ import division, print_function

import os
import abc
import numpy as np

from abipy import abilab

from pymatgen.io.abinitio.flows import AbinitFlow


class FlowRobot(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, path):
        self.flow = abilab.AbinitFlow.pickle_load(path)

    #def walkdirs(self, dir="outdir"):
    #    for task, wi, ti in self.iflat_tasks():
    #        print(wi, ti, task)
    #        #if wi in self.exclude_works or (wi, ti) in self.exclude_tasks: 
    #        #directory = getattr(task, dir)
    #        #if type(task) == abilab.NscfTask:
    #        #    gsr_filepath = directory.has_abiext("GSR")

    def all_files_with_abiext(self, abiext, dirname="outdir", 
                              exclude_works=(), exclude_tasks=(),
                              include_works=(), include_tasks=()):

        files = []
        for task, wi, ti in self.flow.iflat_tasks():

            if wi in exclude_works or (wi, ti) in exclude_tasks: 
                continue

            if ( (include_works and wi not in include_works) or 
                 (include_tasks and (wi, ti) not in include_tasks) ): 
                continue

            directory = getattr(task, dirname)
            path = directory.has_abiext("SIGRES")
            if path:
                files.append(path)

        return files

    @abc.abstractmethod
    def collect_data(self, *args, **kwargs):
        """Collect the data"""
                                                     
    @abc.abstractmethod
    def analyze_data(self, *args, **kwargs):
        """Analyze the data collected previously."""
                                                     
    def run(self, *args, **kwargs):
        self.collect_data(*args, **kwargs)
        return self.analyze_data(*args, **kwargs)


class EbandsRobot(FlowRobot):
    def __init__(self, flow):
        super(EbandsRobot, self).__init__(flow)

        self.plotter = abilab.ElectronBandsPlotter()

    def collect_data(self, *args, **kwargs):

        for task, wi, ti in self.flow.iflat_tasks():
            directory = getattr(task, "outdir")

            if type(task) != abilab.NscfTask:
                continue

            gsr_filepath = directory.has_abiext("GSR")
            if gsr_filepath:
                self.plotter.add_ebands_from_file(gsr_filepath, label=gsr_filepath)

    def analyze_data(self, *args, **kwargs):
        return self.plotter


class EosRobot(FlowRobot):
    # This example shows how to compute the equation of state by 
    # fitting the total energy as function of the unit cell volume.
    #def __init__(self, flow, task_classes=None, work_classes=None):

    def collect_data(self, *args, **kwargs):
        self.volumes = [13.72, 14.83, 16.0, 17.23, 18.52]
        self.energies = [-56.29, -56.41, -56.46, -56.46, -56.42]

        volumes, energies = [], []
        #for task, wi, ti in self.flow.iflat_tasks():
        #    directory = getattr(task, "outdir")
        #    if type(task) != abilab.ScfTask:
        #        continue
        #                                                                            
        #    gsr_filepath = directory.has_abiext("GSR")
        #    if gsr_filepath:
        #        with GSR_Reader(gsr_filepath) as r:
        #           structure = r.read_structure()
        #           etotal = r.read_etotal()
        #           etotals.append(etotal)   
        #           volumes.append(structure.volume)
        #
        #self.etotals = np.array(etotals)
        #self.volumes = np.array(volumes)

    def analyze_data(self, *args, **kwargs):
        # Extract volumes and energies from the output files of the calculation.
        # Here we use hardcoded values.

        fits = []
        for eos_name in abilab.EOS.MODELS:
           eos = abilab.EOS(eos_name=eos_name)
           # Note that eos.fit expects lengths in Angstrom, energies are in eV.
           # To specify different units use len_units and ene_units 
           fits.append(eos.fit(self.volumes, self.energies, vol_unit="bohr^3", ene_unit="Ha"))

        return fits


class SigresRobot(FlowRobot):

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


if __name__ == "__main__":
    import sys
    path = sys.argv[1]

    plotter = SigresRobot(path).run()
    plotter.plot_qpgaps(title="QP gaps vs sigma_nband", hspan=0.05)

    #robot = EosRobot(path)
    #fits = robot.run()
    #for fit in fits: print(fit)

    #robot = EbandsRobot(path)
    #plotter = robot.run()
    #plotter.animate_files()
    #plotter.animate()

    sys.exit(0)
