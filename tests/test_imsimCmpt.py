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
import shutil
import unittest

import numpy as np
import yaml
from lsst.ts.imsim.imsimCmpt import ImsimCmpt
from lsst.ts.imsim.obsMetadata import ObsMetadata
from lsst.ts.imsim.skySim import SkySim
from lsst.ts.imsim.utils.sensorWavefrontError import SensorWavefrontError
from lsst.ts.imsim.utils.utility import (
    getCamera,
    getConfigDir,
    getModulePath,
    getZkFromFile,
)


class TestImsimCmpt(unittest.TestCase):
    """Test the imsim configuration module."""

    maxDiff = None

    def setUp(self):
        self.obsMetadataTest = ObsMetadata(ra=0.0, dec=0.0, band="r", mjd=59580.0)
        self.configPointerDefaultLsstCam = os.path.join(
            getConfigDir(), "lsstCamDefaultPointer.yaml"
        )
        with open(
            os.path.join(
                getModulePath(),
                "tests",
                "testData",
                "imsimConfig",
                "imsimConfigLsstCam.yaml",
            ),
            "r",
        ) as testFile:
            self.fullTestYaml = yaml.safe_load(testFile)
        self.imsimCmpt = ImsimCmpt()

        # Set the output directories
        self.outputDir = os.path.join(getModulePath(), "tests", "tmp")
        self.outputImgDir = os.path.join(self.outputDir, "img")

        self.imsimCmpt.outputDir = self.outputDir
        self.imsimCmpt.outputImgDir = self.outputImgDir

        # Set the file name of analyzed OPD data
        self.zkFileName = "opd.zer"
        self.pssnFileName = "PSSN.txt"

        # Create SkySim object
        self.skySim = SkySim()

    def tearDown(self):
        shutil.rmtree(self.outputDir)

    def _mapSensorNameAndId(self, sensorNameList):
        camera = getCamera("lsst")
        return dict(
            [
                (camera[detector].getName(), camera[detector].getId())
                for detector in sensorNameList
            ]
        )

    def testSetOutputDir(self):
        self.assertTrue(os.path.exists(self.imsimCmpt.outputDir))

    def testSetOutputImgDir(self):
        self.assertTrue(os.path.exists(self.imsimCmpt.outputImgDir))

    def testVerifyPointerFileRaisesError(self):
        requiredKeys = ["input", "gal", "image", "psf", "stamp", "output"]
        pointerWithoutInput = {key: "" for key in requiredKeys[1:]}
        with self.assertRaises(ValueError) as context:
            self.imsimCmpt._verifyPointerFile(pointerWithoutInput, requiredKeys)
        self.assertEqual(
            str(context.exception), "Config pointer file missing filepath for input."
        )

    def testAssembleConfig(self):
        fullConfigYaml = self.imsimCmpt.assembleConfigYaml(
            self.obsMetadataTest, self.configPointerDefaultLsstCam, "lsst"
        )
        self.fullTestYaml["output"]["dir"] = self.imsimCmpt.outputImgDir
        self.assertDictEqual(fullConfigYaml, self.fullTestYaml)

    def testConvertObsMetadataToText(self):
        obsVariablesText = self.imsimCmpt.convertObsMetadataToText(self.obsMetadataTest)
        with open(
            os.path.join(getConfigDir(), "obsVariablesDefault.yaml"), "r"
        ) as file:
            defaultVars = yaml.safe_load(file)
        self.assertDictEqual(yaml.safe_load(obsVariablesText), defaultVars)

    def testFormatOpdTextLsstCam(self):
        opdText = self.imsimCmpt.formatOpdText(self.obsMetadataTest, "lsst")
        self.assertTrue(opdText.startswith("  opd:"))
        self.assertTrue(opdText.endswith("- {thx: 1.176 deg, thy: -1.176 deg}\n"))

    def testAddConfigHeader(self):
        headerText = self.imsimCmpt.addConfigHeader(self.obsMetadataTest)
        headerYaml = yaml.safe_load(headerText)
        for key in list(headerYaml["header"].keys()):
            self.assertEqual(
                headerYaml["header"][key], self.fullTestYaml["output"]["header"][key]
            )

    def testGenInstanceCatalog(self):
        self.skySim.addStarByFile(
            os.path.join(getModulePath(), "tests", "testData", "sky", "wfsSglStar.txt")
        )
        instCat = self.imsimCmpt.genInstanceCatalog(self.skySim)
        expectedCat = "object  0\t 1.196000\t 1.176000 17.000000 "
        expectedCat += (
            "flatSED/sed_flat.txt.gz 0.0 0.0 0.0 0.0 0.0 0.0 point none none \n"
        )
        self.assertEqual(instCat, expectedCat)

    def testGenInstCatStars(self):
        self.skySim.addStarByFile(
            os.path.join(getModulePath(), "tests", "testData", "sky", "wfsSglStar.txt")
        )
        instCat = self.imsimCmpt.genInstCatStars(self.skySim)
        expectedCat = "object  0\t 1.196000\t 1.176000 17.000000 "
        expectedCat += (
            "flatSED/sed_flat.txt.gz 0.0 0.0 0.0 0.0 0.0 0.0 point none none \n"
        )
        self.assertEqual(instCat, expectedCat)

    def testGenerateStar(self):
        starId = 0
        ra = 1.0
        dec = 1.0
        magNorm = 2.0
        sedName = "flat.txt"
        content = self.imsimCmpt.generateStar(starId, ra, dec, magNorm, sedName)
        ansContent = "object  0\t 1.000000\t 1.000000  2.000000 "
        ansContent += "flat.txt 0.0 0.0 0.0 0.0 0.0 0.0 point "
        ansContent += "none none \n"
        self.assertEqual(content, ansContent)

    def testAnalyzeOpdData(self):
        self._analyzeLsstCamOpdData()

        zkFilePath = os.path.join(self.outputImgDir, self.zkFileName)
        pssnFilePath = os.path.join(self.outputImgDir, self.pssnFileName)
        self.assertTrue(os.path.exists(zkFilePath))
        self.assertTrue(os.path.exists(pssnFilePath))

        zk = getZkFromFile(zkFilePath)
        ansZkFilePath = os.path.join(self._getOpdFileDirOfLsstCam(), "opd.zer")
        ansZk = getZkFromFile(ansZkFilePath)

        for det in [191, 195, 199, 203]:
            delta = np.sum(np.abs(zk[det] - ansZk[det]))
            self.assertLess(delta, 1e-10)

        pssnData = np.loadtxt(pssnFilePath)
        pssn = pssnData[0, :]
        ansPssnFilePath = os.path.join(self._getOpdFileDirOfLsstCam(), "PSSN.txt")
        ansPssnData = np.loadtxt(ansPssnFilePath)
        ansPssn = ansPssnData[0, :]

        delta = np.sum(np.abs(pssn - ansPssn))
        self.assertLess(delta, 1e-10)

    def _analyzeLsstCamOpdData(self, rotOpdInDeg=0.0):
        shutil.copy(
            os.path.join(self._getOpdFileDirOfLsstCam(), "opd.fits"),
            self.imsimCmpt.outputImgDir,
        )
        self.imsimCmpt.analyzeOpdData(
            "lsst",
            zkFileName=self.zkFileName,
            rotOpdInDeg=rotOpdInDeg,
            pssnFileName=self.pssnFileName,
        )

    def _getOpdFileDirOfLsstCam(self):
        opdFileDir = os.path.join(getModulePath(), "tests", "testData", "opd")

        return opdFileDir

    def testGetOpdGqEffFwhmFromFile(self):
        self._analyzeLsstCamOpdData()

        gqEffFwhm = self.imsimCmpt.getOpdGqEffFwhmFromFile(self.pssnFileName)
        ansPssnFilePath = os.path.join(self._getOpdFileDirOfLsstCam(), "PSSN.txt")
        ansPssnData = np.loadtxt(ansPssnFilePath)
        ansGqEffFwhm = ansPssnData[1, -1]
        self.assertAlmostEqual(gqEffFwhm, ansGqEffFwhm, places=3)

    def _getRefSensorNameList(self):
        refSensorNameList = [
            "R00_SW0",
            "R44_SW0",
            "R04_SW0",
            "R40_SW0",
        ]

        return refSensorNameList

    def testMapOpdDataToListOfWfErr(self):
        self._analyzeLsstCamOpdData()

        refSensorNameList = self._getRefSensorNameList()
        mapSensorNameAndId = self._mapSensorNameAndId(refSensorNameList)
        ansSensorIdList = list(mapSensorNameAndId.values())

        listOfWfErr = self.imsimCmpt.mapOpdDataToListOfWfErr(
            self.zkFileName, ansSensorIdList, refSensorNameList
        )

        self.assertEqual(len(listOfWfErr), len(refSensorNameList))

        opdZk = getZkFromFile(
            os.path.join(self.imsimCmpt.outputImgDir, self.zkFileName)
        )
        mapSensorNameAndId = self._mapSensorNameAndId(refSensorNameList)
        for wfErr, refSensorName in zip(listOfWfErr, refSensorNameList):
            sensorId = wfErr.sensorId
            sensorName = wfErr.sensorName
            self.assertEqual(sensorName, refSensorName)
            self.assertEqual(sensorId, mapSensorNameAndId[refSensorName])

            zkInWfErr = wfErr.annularZernikePoly
            delta = np.sum(np.abs(zkInWfErr - opdZk[sensorId] / 1e3))
            self.assertEqual(delta, 0)

    def testGetListOfFwhmSensorData(self):
        self._analyzeLsstCamOpdData()
        refSensorNameList = self._getRefSensorNameList()
        mapSensorNameAndId = self._mapSensorNameAndId(refSensorNameList)
        ansSensorIdList = list(mapSensorNameAndId.values())

        (
            sensor_data_fwhm,
            sensor_data_sensor_id,
        ) = self.imsimCmpt.getListOfFwhmSensorData(self.pssnFileName, ansSensorIdList)
        self.assertEqual(len(sensor_data_fwhm), len(ansSensorIdList))

        ansData = self.imsimCmpt._getDataOfPssnFile(self.pssnFileName)
        ansFwhmData = ansData[1, :-1]

        for data_fwhm, sensor_id, ansSensorId, ansFwhm in zip(
            sensor_data_fwhm, sensor_data_sensor_id, ansSensorIdList, ansFwhmData
        ):
            self.assertEqual(sensor_id, ansSensorId)
            self.assertEqual(data_fwhm, ansFwhm)

    def testGetOpdPssnFromFile(self):
        self._analyzeLsstCamOpdData()

        # The correctness of values have been tested at the test case of
        # testAnalyzeOpdData
        pssn = self.imsimCmpt.getOpdPssnFromFile(self.pssnFileName)
        self.assertEqual(len(pssn), 4)

    def testReorderAndSaveWfErrFile(self):
        listOfWfErr = self._prepareListOfWfErr()

        refSensorNameList = self._getRefSensorNameList()
        mapSensorNameAndId = self._mapSensorNameAndId(refSensorNameList)
        sensorIdList = list(mapSensorNameAndId.values())

        zkFileName = "testZk.zer"
        camera = getCamera("lsst")
        self.imsimCmpt.reorderAndSaveWfErrFile(
            listOfWfErr, refSensorNameList, camera, zkFileName=zkFileName
        )

        zkFilePath = os.path.join(self.imsimCmpt.outputImgDir, zkFileName)
        zkInFile = getZkFromFile(zkFilePath)

        numOfZk = self.imsimCmpt.numOfZk
        self.assertEqual(len(zkInFile), len(refSensorNameList))
        self.assertEqual(
            len(zkInFile[sensorIdList[0]]),
            numOfZk,
        )

        self.assertEqual(np.sum(zkInFile[191]), 0)
        self.assertEqual(np.sum(zkInFile[203]), 0)

        delta = np.sum(np.abs(zkInFile[195] - listOfWfErr[0].annularZernikePoly))
        self.assertLess(delta, 1e-7)

        delta = np.sum(np.abs(zkInFile[199] - listOfWfErr[1].annularZernikePoly))
        self.assertLess(delta, 1e-7)

    def _prepareListOfWfErr(self):
        numOfZk = self.imsimCmpt.numOfZk

        sensorIdList = [195, 199]
        listOfWfErr = []
        for sensorId in sensorIdList:
            sensorWavefrontData = SensorWavefrontError()
            sensorWavefrontData.sensorId = sensorId

            wfErr = np.random.rand(numOfZk)
            sensorWavefrontData.annularZernikePoly = wfErr

            listOfWfErr.append(sensorWavefrontData)

        return listOfWfErr

    def testSaveDofInUmFileForNextIter(self):
        self.imsimCmpt.dofInUm = np.arange(50)
        dofInUmFileName = "dofPertInNextIter.mat"
        self.imsimCmpt.saveDofInUmFileForNextIter(dofInUmFileName=dofInUmFileName)

        filePath = os.path.join(self.outputDir, dofInUmFileName)
        self.assertTrue(os.path.exists(filePath))

        data = np.loadtxt(filePath)
        delta = np.sum(np.abs(self.imsimCmpt.dofInUm - data))
        self.assertEqual(delta, 0)
