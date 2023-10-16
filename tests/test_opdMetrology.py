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

import lsst.obs.lsst as obs_lsst
import numpy as np
from astropy.io import fits
from lsst.afw.cameraGeom import FIELD_ANGLE
from lsst.ts.imsim.opdMetrology import OpdMetrology
from lsst.ts.imsim.utils.utility import getModulePath, getZkFromFile


class TestOpdMetrology(unittest.TestCase):
    """Test the OpdMetrology class."""

    def setUp(self):
        self.metr = OpdMetrology()
        self.testDataDir = os.path.join(getModulePath(), "tests", "testData")

    def testSetWeightingRatio(self):
        wt = [1, 2]
        self.metr.wt = wt

        wtInMetr = self.metr.wt
        self.assertEqual(len(wtInMetr), len(wt))
        self.assertEqual(np.sum(wtInMetr), 1)
        self.assertAlmostEqual(wtInMetr[1] / wtInMetr[0], 2)
        with self.assertRaises(ValueError) as context:
            self.metr.wt = [-1, 1]
        self.assertEqual(str(context.exception), "All weighting ratios should be >= 0.")

    def testSetWgtAndFieldXyOfGQLsstFam(self):
        self.metr.setWgtAndFieldXyOfGQ("lsstfam")

        fieldX = self.metr.fieldX
        self.assertEqual(len(fieldX), 189)

        wgt = self.metr.wt
        self.assertEqual(len(wgt), 189)

    def testSetWgtAndFieldXyOfGQErr(self):
        self.assertRaises(
            RuntimeError, self.metr.setWgtAndFieldXyOfGQ, "NoThisInstName"
        )

    def testSetDefaultLsstWfsGQ(self):
        self.metr.setDefaultLsstWfsGQ()

        fieldX = self.metr.fieldX
        fieldY = self.metr.fieldY
        # Get true values from obs_lsst
        camera = obs_lsst.LsstCam.getCamera()
        detIdMap = camera.getIdMap()
        trueX = []
        trueY = []
        for sens_id in [192, 196, 200, 204]:
            centerIntra = detIdMap[sens_id].getCenter(FIELD_ANGLE)
            centerExtra = detIdMap[sens_id - 1].getCenter(FIELD_ANGLE)
            centerDeg = np.degrees(np.array([centerIntra, centerExtra]))
            centerX = np.mean(centerDeg[:, 1])
            centerY = np.mean(centerDeg[:, 0])
            trueX.append(centerX)
            trueY.append(centerY)

        self.assertCountEqual(fieldX, trueX)
        self.assertCountEqual(
            fieldY,
            trueY,
        )

        wgt = self.metr.wt
        self.assertCountEqual(wgt, [0.25, 0.25, 0.25, 0.25])

    def testGetDefaultLsstWfsGQ(self):
        fieldWFSx, fieldWFSy, detIds = self.metr.getDefaultLsstWfsGQ()
        self.assertEqual(len(fieldWFSx), 4)
        self.assertEqual(len(detIds), 4)

    def _getOpdDir(self):
        opdFileDir = os.path.join(self.testDataDir, "opd")
        return opdFileDir

    def testGetZkFromOpd(self):
        opdDir = self._getOpdDir()
        zk = self.metr.getZkFromOpd(opdFitsFile=os.path.join(opdDir, "opd.fits"))[0]

        ansOpdFileName = "opd.zer"
        ansOpdFilePath = os.path.join(opdDir, ansOpdFileName)
        allOpdAns = getZkFromFile(ansOpdFilePath)
        self.assertLess(np.sum(np.abs(zk[3:] - allOpdAns[191])), 1e-5)

    def testRPTTfromOPD(self):
        """Test removal of piston (z1), x-tilt (z2), and y-tilt (z3)
        from the OPD map."""
        opdDir = self._getOpdDir()
        opdFilePath = os.path.join(opdDir, "opd.fits")
        opdMap = fits.getdata(opdFilePath, 1)
        opdRmPTT, opdx, opdy = self.metr.rmPTTfromOPD(opdMap=opdMap)

        # Flip OPD because it will be flipped inside getZkFromOpd
        zkRmPTT = self.metr.getZkFromOpd(opdMap=opdRmPTT)[0]
        zkRmPTTInUm = np.sum(np.abs(zkRmPTT[0:3])) / 1e3
        self.assertLess(zkRmPTTInUm, 9e-2)

    def testCalcPSSN(self):
        pssn = self._calcPssn()
        allData = self._getMetroAllAnsData()
        self.assertAlmostEqual(pssn, allData[0, 0])

    def _calcPssn(self):
        wavelengthInUm = 0.48
        opdFilePath = os.path.join(self._getOpdDir(), "opd.fits")
        opdMap = fits.getdata(opdFilePath, 0) * 1e-3
        pssn = self.metr.calcPSSN(wavelengthInUm, opdMap=opdMap)

        return pssn

    def _getMetroAllAnsData(self):
        ansAllDataFileName = "PSSN.txt"
        ansAllDataFilePath = os.path.join(self._getOpdDir(), ansAllDataFileName)
        allData = np.loadtxt(ansAllDataFilePath)

        return allData

    def testCalcFWHMeff(self):
        pssn = self._calcPssn()
        fwhm = self.metr.calcFWHMeff(pssn)

        allData = self._getMetroAllAnsData()
        self.assertAlmostEqual(fwhm, allData[1, 0])

    def testCalcGQvalue(self):
        self.metr.setDefaultLsstWfsGQ()
        allData = self._getMetroAllAnsData()
        valueList = allData[0, 0:4]

        GQvalue = self.metr.calcGQvalue(valueList)
        self.assertAlmostEqual(GQvalue, allData[0, -1])
