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

import numpy as np
from astropy.io import fits
from lsst.ts.imsim.opdMetrology import OpdMetrology
from lsst.ts.imsim.utils.utility import getModulePath, getZkFromFile


class TestOpdMetrology(unittest.TestCase):
    """Test the OpdMetrology class."""

    def setUp(self):
        self.metr = OpdMetrology()
        self.testDataDir = os.path.join(getModulePath(), "tests", "testData")

    def testSetCamera(self):
        self.metr.setCamera("lsstfam")
        self.assertEqual(self.metr._camera.getName(), "LSSTCam")
        self.metr.setCamera("comcam")
        self.assertEqual(self.metr._camera.getName(), "LSSTComCam")

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

    def testSetWgtAndFieldXyOfGQLsst(self):
        self.metr.setWgtAndFieldXyOfGQ("lsst")

        fieldXAns, fieldYAns, wgtAns = self._calcFieldXyAndWgtLsst()

        fieldX = self.metr.fieldX
        fieldY = self.metr.fieldY
        self.assertEqual(len(fieldX), 31)
        self.assertLess(np.sum(np.abs(fieldX - fieldXAns)), 1e-10)
        self.assertLess(np.sum(np.abs(fieldY - fieldYAns)), 1e-10)

        wgt = self.metr.wt
        self.assertEqual(len(wgt), 31)
        self.assertLess(np.sum(np.abs(wgt - wgtAns)), 1e-10)

    def _calcFieldXyAndWgtLsst(self):
        # The distance of point xi (used in Gaussian quadrature plane) to the
        # origin
        # This value is in [-1.75, 1.75]
        armLen = [0.379, 0.841, 1.237, 1.535, 1.708]

        # Weighting of point xi (used in Gaussian quadrature plane) for each
        # ring
        armW = [0.2369, 0.4786, 0.5689, 0.4786, 0.2369]

        # Number of points on each ring
        nArm = 6

        # Get the weighting for all field points (31 for lsst camera)
        # Consider the first element is center (0)
        wgt = np.concatenate([np.zeros(1), np.kron(armW, np.ones(nArm))])

        # Generate the fields point x, y coordinates
        pointAngle = np.arange(nArm) * (2 * np.pi) / nArm
        fieldX = np.concatenate([np.zeros(1), np.kron(armLen, np.cos(pointAngle))])
        fieldY = np.concatenate([np.zeros(1), np.kron(armLen, np.sin(pointAngle))])

        return fieldX, fieldY, wgt / np.sum(wgt)

    def testSetWgtAndFieldXyOfGQComCam(self):
        self.metr.setWgtAndFieldXyOfGQ("comcam")

        fieldXAns, fieldYAns, wgtAns = self._calcFieldXyAndWgtComCam()

        fieldX = self.metr.fieldX
        fieldY = self.metr.fieldY
        self.assertEqual(len(fieldX), 9)
        self.assertLess(np.sum(np.abs(fieldX - fieldXAns)), 1e-10)
        self.assertLess(np.sum(np.abs(fieldY - fieldYAns)), 1e-10)

        wgt = self.metr.wt
        self.assertEqual(len(wgt), 9)
        self.assertLess(np.sum(np.abs(wgt - wgtAns)), 1e-10)

    def _calcFieldXyAndWgtComCam(self):
        # ComCam is the cetral raft of LSST cam, which is composed of 3 x 3
        # CCDs.
        nRow = 3
        nCol = 3

        # Number of field points
        nField = nRow * nCol

        # Get the weighting for all field points (9 for comcam)
        wgt = np.ones(nField)

        # Distance to raft center in degree along x/y direction and the
        # related relative position
        sensorD = 0.2347
        coorComcam = sensorD * np.array([-1, 0, 1])

        # Generate the fields point x, y coordinates
        fieldX = np.kron(coorComcam, np.ones(nRow))
        fieldY = np.kron(np.ones(nCol), coorComcam)

        return fieldX, fieldY, wgt / np.sum(wgt)

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
        self.assertCountEqual(fieldX, [1.176, -1.176, -1.176, 1.176])
        self.assertCountEqual(fieldY, [1.176, -1.176, -1.176, 1.176])

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
        print(zk, allOpdAns)
        self.assertLess(np.sum(np.abs(zk[3:] - allOpdAns[203])), 1e-5)

    def testRmPTTfromOPD(self):
        opdDir = self._getOpdDir()
        opdFilePath = os.path.join(opdDir, "opd.fits")
        opdRmPTT, opdx, opdy = self.metr.rmPTTfromOPD(opdFitsFile=opdFilePath)

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
