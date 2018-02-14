"""Tests for seabornw module."""
from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np
import pandas as pd
import seaborn as sns
import abipy.display.seabornw as snw

from abipy.core.testing import AbipyTest


class SeabornWidgetTest(AbipyTest):

    def test_seabornw(self):
        """Test ipython widget for seaborn."""
        if not self.has_ipywidgets():
            raise self.SkipTest("This test requires ipywidgets")

        tips = sns.load_dataset("tips")
        snw.jointplot(tips)

        titanic = sns.load_dataset("titanic")
        snw.countplot(titanic)

        tips = sns.load_dataset("tips")
        snw.jointplot(tips)
        snw.swarmplot(tips)
        snw.lmplot(tips)

        exercise = sns.load_dataset("exercise")
        snw.factorplot(exercise)
        snw.violinplot(tips)
        snw.stripplot(tips)
        snw.swarmplot(tips)
        snw.pointplot(tips)
        snw.barplot(tips)

        np.random.seed(0)
        uniform_data = np.random.rand(10, 12)
        snw.heatmap(uniform_data)

        flights = sns.load_dataset("flights")
        flights = flights.pivot("month", "year", "passengers")
        snw.clustermap(flights)

        # Generate a random dataset with strong simple effects and an interaction
        #n = 80
        #rs = np.random.RandomState(11)
        #x1 = rs.randn(n)
        #x2 = x1 / 5 + rs.randn(n)
        #b0, b1, b2, b3 = .5, .25, -1, 2
        #y = b0  + b1 * x1 + b2 * x2 + b3 * x1 * x2 + rs.randn(n)
        #df = pd.DataFrame(np.c_[x1, x2, y], columns=["x1", "x2", "y"])

        # Show a scatterplot of the predictors with the estimated model surface
        #snw.interactplot(df)
