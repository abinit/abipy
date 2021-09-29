""""Panels to interact with the AbiPy tasks."""
import param
import panel as pn
import panel.widgets as pnw

from abipy.panels.core import mpl, ply, dfc, depends_on_btn_click
from abipy.panels.nodes import NodeParameterized


class WorkPanel(NodeParameterized):
    """
    Panel to interact with an AbiPy Work
    """

    def __init__(self, work, **params):
        NodeParameterized.__init__(self, node=work, **params)
        self.work = work

    #def get_panel(self, as_dict=False, **kwargs):
    #    """Return tabs with widgets to interact with the flow."""

    #    return super().get_panel(as_dict=as_dict, **kwargs)
