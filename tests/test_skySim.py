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

from astroplan import FixedTarget, Observer
from astropy.time import Time
from lsst.ts.imsim.obsMetadata import ObsMetadata
from lsst.ts.imsim.skySim import SkySim
from lsst.ts.imsim.utils.utility import getModulePath


class TestSkySim(unittest.TestCase):
    def setUp(self):
        self.skySim = SkySim()

    def testSetCamera(self):
        self.skySim.setCamera("lsstfam")
        self.assertEqual(self.skySim._camera.getName(), "LSSTCam")

    def testCalcParallacticAngle(self):
        sirius = FixedTarget.from_name("sirius")
        t = Time(60000, format="mjd")
        obsMetadata = ObsMetadata(
            ra=sirius.ra.deg, dec=sirius.dec.deg, band="r", mjd=60000
        )
        rubin = Observer.at_site("cerro pachon")
        self.assertEqual(
            self.skySim.calcParallacticAngle(obsMetadata),
            rubin.parallactic_angle(t, sirius).deg,
        )

    def testAddStarByRaDecInDeg(self):
        self.skySim.addStarByRaDecInDeg(1, 2, 3, 4)
        self.assertEqual(len(self.skySim.starId), 1)

        self.skySim.addStarByRaDecInDeg(2, 2.1, 3, 4)
        self.assertEqual(len(self.skySim.starId), 2)
        self.assertEqual(self.skySim.starId[0], 1)
        self.assertEqual(self.skySim.starId[1], 2)

        # Try to add the same star Id again
        self.skySim.addStarByRaDecInDeg(2, 2.1, 3, 4)
        self.assertEqual(len(self.skySim.starId), 2)

    def testAddStarByFile(self):
        self._addStarByFile("wfsStar.txt")

        self.assertEqual(len(self.skySim.starId), 8)

        ra = self.skySim.ra
        dec = self.skySim.dec
        self.assertEqual(ra[2], -1.176)
        self.assertEqual(dec[2], 1.216)

        self.assertEqual(self.skySim.mag[2], 15.0)

    def _addStarByFile(self, skyFileName):
        skyFile = os.path.join(getModulePath(), "tests", "testData", "sky", skyFileName)
        self.skySim.addStarByFile(skyFile)

    def testAddStarByFileWithSglStar(self):
        self._addStarByFile("wfsSglStar.txt")

        self.assertEqual(len(self.skySim.starId), 1)

        ra = self.skySim.ra
        dec = self.skySim.dec
        self.assertEqual(ra[0], 1.196)
        self.assertEqual(dec[0], 1.176)

        self.assertEqual(self.skySim.mag[0], 17.0)

    def testExportSkyToFile(self):
        self._addStarByFile("wfsStar.txt")
        outputFilePath = os.path.join(getModulePath(), "output", "testSkyOutput.txt")

        self.skySim.exportSkyToFile(outputFilePath)
        self.assertTrue(os.path.isfile(outputFilePath))
        os.remove(outputFilePath)
