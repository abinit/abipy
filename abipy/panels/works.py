"""Panels to interact with the AbiPy tasks."""

from __future__ import annotations

from abipy.panels.nodes import NodeParameterized


class WorkPanel(NodeParameterized):
    """Panel to interact with an AbiPy Work."""

    def __init__(self, work, **params):
        """
        Args:
            work: |Work| object.
            params: Parameters passed to the parent class.
        """
        NodeParameterized.__init__(self, node=work, **params)
        self.work = work

    # def get_panel(self, as_dict=False, **kwargs):
    #    """Return tabs with widgets to interact with the flow."""

    #    return super().get_panel(as_dict=as_dict, **kwargs)
