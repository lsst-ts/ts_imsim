# This file is part of ts_imsim.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import unittest

from lsst.ts.imsim.utils.plotUtil import plotFwhmOfIters
from lsst.ts.imsim.utils.utility import getModulePath


class TestPlotUtil(unittest.TestCase):
    """Test the PlotUtil functions."""

    def setUp(self):
        modulePath = getModulePath()
        self.testData = os.path.join(modulePath, "tests", "testData")
        self.outFigFilePath = os.path.join(modulePath, "output", "img", "testFig.png")

    def tearDown(self):
        if os.path.exists(self.outFigFilePath):
            os.remove(self.outFigFilePath)

    def testPlotFwhmOfIters(self):
        iterDataDir = os.path.join(self.testData, "iterData")
        pssnFiles = [
            os.path.join(iterDataDir, "iter%d" % num, "img", "PSSN.txt")
            for num in range(5)
        ]

        plotFwhmOfIters(pssnFiles, saveToFilePath=self.outFigFilePath)
        self.assertTrue(os.path.exists(self.outFigFilePath))


if __name__ == "__main__":
    # Run the unit test
    unittest.main()
