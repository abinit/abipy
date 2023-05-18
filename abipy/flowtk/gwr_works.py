# coding: utf-8
"""
Works and Flows for GWR calculations (GW with supercells)
"""
from __future__ import annotations

import numpy as np

from abipy.abio.inputs import AbinitInput, RUNL, GWR_TASK
from .nodes import Node
from .tasks import ScfTask, NscfTask, TaskManager
from .works import Work


class DirectDiagoWork(Work):
    """
    This work performs the direct diagonalization of the KS Hamiltonian using the 
    density produced by the scf_task and produces a WFK file with empty states.
    in the outdir of the second task

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: DirectDiagoWork
    """

    @classmethod 
    def from_scf_input(cls, scf_input: AbinitInput, green_nband: int, 
                      manager: TaskManager=None) -> DirectDiagoWork:
        """
        Args:
            scf_input: Input for GS-SCF calculation.
            green_nband: Number of bands to compute in the direct diagonalizatio.
                A negative value activate full diagonalization with nband equal to
                the number of PWs.
        """
        work = cls(manager=manager)
        work.green_nband = green_nband
        gwr_task = GWR_TASK.HDIAGO_FULL if green_nband < 0 else GWR_TASK.HDIAGO
        work.scf_task = work.register_scf_task(scf_input)
        diago_input = scf_input.new_with_vars(optdriver=RUNL.GWR, 
                                              gwr_task=gwr_task,
                                              )
        work.diago_task = work.register_gwr_task(diago_input, deps={work.scf_task: "DEN"})
        return work


class _BaseGwrWork(Work):

    @classmethod                                                                     
    def from_varname_values(cls, varname_values: tuple, gwr_template: AbinitInput,            
                            den_node: Node, wfk_node: Node,                               
                            manager: TaskManager = None):
        """                                                                          
        Args:                                                                        
            varname_values:                                                                
            gwr_template:                                                        
            den_node:                                                                
            wfk_node:                                                                
            manager:                                                                 
    
        Example: 

            varname_values = ("nband", [8, 12, 14])
            Work.from_varname_values(varname_values, gwr_template, den_node, wfk_node)

            or:

            varname_values = [
                ("nband", [8, 12]),             
                ("gwr_ntau", [6, 8]),               
            ]             
            Work.from_varname_values(varname_values, gwr_template, den_node, wfk_node)
        """                                                                          
        work = cls(manager=manager)                                                  
        for new_inp in gwr_template.news_varname_values(varname_values):
            work.register_gwr_task(new_inp, deps={den_node: "DEN", wfk_node: "WFK"})
                                                                                     
        return work                                                                  
                     
        
class GwrSigmaConv(_BaseGwrWork):
    """

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GwrSigmaConv
    """


class GwrRpaEneConv(_BaseGwrWork):
    """
    This work computes the RPA energy for different number of points in 
    the minimax mesh
    """

    def on_all_ok(self):
        """
        """
        return super().on_all_ok()
