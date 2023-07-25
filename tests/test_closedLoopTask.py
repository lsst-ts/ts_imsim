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
import argparse
import tempfile
import unittest

from lsst.ts.wep.utility import CamType, FilterType

from lsst.ts.imsim.obsMetadata import ObsMetadata
from lsst.ts.imsim.utils.utility import getModulePath
from lsst.ts.imsim.closedLoopTask import ClosedLoopTask


class TestclosedLoopTask(unittest.TestCase):
    """Test the closedLoopTask class."""

    def setUp(self):
        self.closedLoopTask = ClosedLoopTask()
        self.obsMetadata = ObsMetadata(0.0, 0.0, "r")

        rootTestDir = os.path.join(getModulePath(), "tests")
        self.testDir = tempfile.TemporaryDirectory(dir=rootTestDir)

    def tearDown(self):
        self.testDir.cleanup()

    def testConfigSkySimWithError(self):
        self.assertRaises(
            ValueError,
            self.closedLoopTask.configSkySim,
            "NoThisInstName",
            self.obsMetadata,
        )

    def testConfigSkySimNoSkyFileLSST(self):
        self.closedLoopTask.configSkySim("lsst", self.obsMetadata)

        skySim = self.closedLoopTask.skySim
        self.assertEqual(len(skySim.starId), 8)

    def testConfigSkySimWithSkyFile(self):
        testSkyFile = os.path.join(
            getModulePath(), "tests", "testData", "sky", "wfsStar.txt"
        )
        self.closedLoopTask.configSkySim(
            "lsst", self.obsMetadata, pathSkyFile=testSkyFile
        )

        skySim = self.closedLoopTask.skySim
        self.assertEqual(len(skySim.starId), 8)
        self.assertEqual(skySim.mag[0], 15)

    def testConfigOfcCalc(self):
        instName = "lsst"
        self.closedLoopTask.configOfcCalc(instName)

        ofcCalc = self.closedLoopTask.ofcCalc
        self.assertEqual(ofcCalc.ofc_data.name, instName)

    def testMapFilterRefToG(self):
        # test that the reference filter
        # gets mapped to g
        for filterTypeName in ["ref", ""]:
            mappedFilterName = self.closedLoopTask.mapFilterRefToG(filterTypeName)
            self.assertEqual(mappedFilterName, "g")

        # test that all other filters are
        # mapped to themselves
        for filterTypeName in "ugrizy":
            mappedFilterName = self.closedLoopTask.mapFilterRefToG(filterTypeName)
            self.assertEqual(mappedFilterName, filterTypeName)

    def testGetFilterTypeRef(self):
        filterType = self.closedLoopTask.getFilterType("ref")
        self.assertEqual(filterType, FilterType.REF)

    def testGetFilterTypeU(self):
        filterType = self.closedLoopTask.getFilterType("u")
        self.assertEqual(filterType, FilterType.LSST_U)

    def testGetFilterTypeG(self):
        filterType = self.closedLoopTask.getFilterType("g")
        self.assertEqual(filterType, FilterType.LSST_G)

    def testGetFilterTypeR(self):
        filterType = self.closedLoopTask.getFilterType("r")
        self.assertEqual(filterType, FilterType.LSST_R)

    def testGetFilterTypeI(self):
        filterType = self.closedLoopTask.getFilterType("i")
        self.assertEqual(filterType, FilterType.LSST_I)

    def testGetFilterTypeZ(self):
        filterType = self.closedLoopTask.getFilterType("z")
        self.assertEqual(filterType, FilterType.LSST_Z)

    def testGetFilterTypeY(self):
        filterType = self.closedLoopTask.getFilterType("y")
        self.assertEqual(filterType, FilterType.LSST_Y)

    def testGetFilterTypeErr(self):
        self.assertRaises(
            ValueError, self.closedLoopTask.getFilterType, "noThisFilterType"
        )

    def testGetCamTypeAndInstNameLsst(self):
        camType, instName = self.closedLoopTask.getCamTypeAndInstName("lsst")
        self.assertEqual(camType, CamType.LsstCam)
        self.assertEqual(instName, "lsst")

    def testGetCamTypeAndInstNameErr(self):
        self.assertRaises(
            ValueError, self.closedLoopTask.getCamTypeAndInstName, "noThisInst"
        )

    def testEraseDirectoryContent(self):
        # Make the temporary directory
        tempDir = os.path.join(self.testDir.name, "tempDir")
        os.mkdir(tempDir)
        files = os.listdir(self.testDir.name)
        self.assertEqual(len(files), 1)

        # Try to erase the content
        self.closedLoopTask.eraseDirectoryContent(self.testDir.name)

        files = os.listdir(self.testDir.name)
        self.assertEqual(len(files), 0)

    def testCheckBoresight(self):
        self.assertRaises(ValueError, self.closedLoopTask.checkBoresight, [-1, 0])
        self.assertRaises(ValueError, self.closedLoopTask.checkBoresight, [361, 0])
        self.assertRaises(ValueError, self.closedLoopTask.checkBoresight, [0, -91])
        self.assertRaises(ValueError, self.closedLoopTask.checkBoresight, [0, 91])

    def testCheckAndCreateBaseOutputDir(self):
        # check that the output dir is created
        # where the name is given
        baseOutputDir = os.path.join(self.testDir.name, "testBaseOutputDir")
        self.assertFalse(os.path.exists(baseOutputDir))
        self.closedLoopTask.checkAndCreateBaseOutputDir(baseOutputDir)
        self.assertTrue(os.path.exists(baseOutputDir))

    def testSetDefaultParser(self):
        parser = argparse.ArgumentParser()
        parser = ClosedLoopTask.setDefaultParser(parser)

        args = parser.parse_known_args()[0]
        self.assertEqual(args.inst, "lsst")
        self.assertEqual(args.filterType, "")
        self.assertEqual(args.rotCam, 0.0)
        self.assertEqual(args.iterNum, 5)
        self.assertEqual(args.output, "")
        self.assertFalse(args.clobber)
        self.assertEqual(args.configPointerFile, "")
        self.assertEqual(args.pipelineFile, "")
        self.assertEqual(args.skySeed, 42)
        self.assertEqual(args.pertSeed, 11)

    def testSetImgParser(self):
        parser = argparse.ArgumentParser()
        parser = ClosedLoopTask.setImgParser(parser)

        argsToTest = ["--boresightDeg", "1.2", "2.3"]

        args = parser.parse_known_args(args=argsToTest)[0]
        self.assertEqual(args.boresightDeg, [1.2, 2.3])
        self.assertEqual(args.skyFile, "")
        self.assertEqual(args.mjd, 59580)
        self.assertFalse(args.turnOffSkyBackground)
        self.assertFalse(args.turnOffAtmosphere)

    def testGetSensorNameListOfFieldsLsstWfs(self):
        sensorNameList = self.closedLoopTask.getSensorNameListOfFields("lsst")
        self.assertEqual(len(sensorNameList), 8)

        sensorNameListAns = [
            "R00_SW0",
            "R00_SW1",
            "R40_SW0",
            "R40_SW1",
            "R04_SW0",
            "R04_SW1",
            "R44_SW0",
            "R44_SW1",
        ]
        self.assertCountEqual(sensorNameList, sensorNameListAns)

    def testGetSensorIdListOfFields(self):
        sensorIdList = self.closedLoopTask.getSensorIdListOfFields("lsst")
        self.assertEqual(len(sensorIdList), 8)

        sensorIdListAns = [191, 192, 195, 196, 199, 200, 203, 204]
        self.assertCountEqual(sensorIdList, sensorIdListAns)
