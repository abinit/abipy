# coding: utf-8
""""""
from __future__ import division, print_function, absolute_import

from collections import namedtuple

import numpy as np
import pandas as pd


class MyDataFrame(pd.DataFrame):
    """
    Extends :class:`DataFrame` with methods for numerical analysis.
    """
    def quadfit(self, xname=None, yname=None):
        """
        Quadratic fit. Return namedtuple with the parameters of the fix a*x**2 + b*x + c,
        the position of the minimum in x0 and the value of the minimum in y0
        """
        xvals, yvals = self[xname], self[yname]

        a, b, c = np.polyfit(xvals, yvals, 2)
        x0 = -b/(2*a)
        y0 = a*x0**2 + b*x0 + c

        return namedtuple("quadfit_results", "a, b, c, x0, y0")(a=a, b=b, c=c, x0=x0, y0=y0)

    #def open_in_excel(self, index=True, excel_path="excel.exe", tmp_path='.'):
    #    """
    #    Open the dataframe in excel.

    #    Args:
    #       excel_path:path to your copy of excel
    #       index: True - export the index of the dataframe as the first columns
    #       tmp_path:  - directory to save the file in

    #    This creates a temporary file name, exports the dataframe to a csv of that file name,
    #    and then tells excel to open the file (in read only mode). (It uses df.to_csv instead
    #    of to_excel because if you don't have excel, you still get the csv.)

    #    Note - this does NOT delete the file when you exit. 
    #    """
    #    import tempfile
    #    f = tempfile.NamedTemporaryFile(delete=False, dir=tmp_path, suffix='.csv', prefix='tmp_')
    #    tmp_name = f.name
    #    f.close()

    #    self.to_csv(tmp_name, index=index)
    #    cmd = [excel_path, '/r', '/e', tmp_name]
    #    try:
    #        ret_val=subprocess.Popen(cmd).pid
    #    except:
    #        print("open_in_excel(): failed to open excel")
    #        print("filename = ", tmp_name)
    #        print("command line = ", cmd)
    #        print("Unexpected error:", sys.exc_info()[0])
