"""Tests for pandasw module."""
import abipy.display.pandasw as pdw

from abipy.core.testing import AbipyTest


class PandasWidgetTest(AbipyTest):

    def test_pandasw(self):
        """Test ipython widget for pandas."""
        if not self.has_ipywidgets():
            raise self.SkipTest("This test requires ipywidgets")

        import seaborn as sns
        data = sns.load_dataset("iris")
        assert pdw.plot(data)
