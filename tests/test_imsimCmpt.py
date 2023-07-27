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
import yaml
import unittest
import shutil
import numpy as np

from lsst.ts.imsim.obsMetadata import ObsMetadata
from lsst.ts.imsim.utils.utility import getConfigDir, getModulePath, getCamera
from lsst.ts.imsim.imsimCmpt import ImsimCmpt
from lsst.ts.imsim.utils.sensorWavefrontError import SensorWavefrontError


class TestImsimCmpt(unittest.TestCase):
    """Test the imsim configuration module."""

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

    def tearDown(self):
        shutil.rmtree(self.outputDir)

    def _mapSensorNameAndId(self):
        return dict(
            [(detector.getId(), detector.getName()) for detector in getCamera("lsst")]
        )

    def testSetOutputDir(self):
        self.assertTrue(os.path.exists(self.imsimCmpt.outputDir))

    def testSetOutputImgDir(self):
        self.assertTrue(os.path.exists(self.imsimCmpt.outputImgDir))

    def _getRefSensorNameList(self):
        refSensorNameList = [
            "R22_S00",
            "R22_S01",
            "R22_S02",
            "R22_S10",
            "R22_S11",
            "R22_S12",
            "R22_S20",
            "R22_S21",
            "R22_S22",
        ]

        return refSensorNameList

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

    def testRunImsim(self):
        pass

    def testGenInstanceCatalog(self):
        pass

    def testGenInstCatStars(self):
        pass

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

        zk = np.loadtxt(zkFilePath)
        ansZkFilePath = os.path.join(
            self._getOpdFileDirOfComCam(), "sim7_iter0_opd.zer"
        )
        ansZk = np.loadtxt(ansZkFilePath)

        delta = np.sum(np.abs(zk - ansZk[:, 3:]))
        self.assertLess(delta, 1e-10)

        pssnData = np.loadtxt(pssnFilePath)
        pssn = pssnData[0, :]
        ansPssnFilePath = os.path.join(
            self._getOpdFileDirOfLsstCam(), "sim7_iter0_PSSN.txt"
        )
        ansPssnData = np.loadtxt(ansPssnFilePath)
        ansPssn = ansPssnData[0, :]

        delta = np.sum(np.abs(pssn - ansPssn))
        self.assertLess(delta, 1e-10)

    def _analyzeLsstCamOpdData(self, rotOpdInDeg=0.0):
        self.imsimCmpt.analyzeOpdData(
            "lsst",
            zkFileName=self.zkFileName,
            rotOpdInDeg=rotOpdInDeg,
            pssnFileName=self.pssnFileName,
        )

    def _getOpdFileDirOfLsstCam(self):
        opdFileDir = os.path.join(
            getModulePath(), "tests", "testData", "lsstOpdFile", "iter0"
        )

        return opdFileDir

    def testGetOpdGqEffFwhmFromFile(self):
        self._analyzeLsstCamOpdData()

        gqEffFwhm = self.imsimCmpt.getOpdGqEffFwhmFromFile(self.pssnFileName)
        self.assertAlmostEqual(gqEffFwhm, 0.5534, places=3)

    def testMapOpdDataToListOfWfErr(self):
        self._analyzeLsstCamOpdData()

        refSensorNameList = self._getRefSensorNameList()
        mapSensorNameAndId = self._mapSensorNameAndId()
        ansSensorIdList = mapSensorNameAndId.mapSensorNameToId(refSensorNameList)

        listOfWfErr = self.imsimCmpt.mapOpdDataToListOfWfErr(
            self.zkFileName, ansSensorIdList, refSensorNameList
        )

        self.assertEqual(len(listOfWfErr), len(refSensorNameList))

        opdZk = self.imsimCmpt._getZkFromFile(self.zkFileName)
        mapSensorNameAndId = self._mapSensorNameAndId()
        for wfErr, refSensorName, zk in zip(listOfWfErr, refSensorNameList, opdZk):
            sensorId = wfErr.sensorId
            sensorNameList = mapSensorNameAndId.mapSensorIdToName(sensorId)[0]
            self.assertEqual(sensorNameList[0], refSensorName)

            zkInWfErr = wfErr.getAnnularZernikePoly()
            delta = np.sum(np.abs(zkInWfErr - zk))
            self.assertEqual(delta, 0)

    def testGetListOfFwhmSensorData(self):
        self._analyzeLsstCamOpdData()
        refSensorNameList = self._getRefSensorNameList()
        mapSensorNameAndId = self._mapSensorNameAndId()
        ansSensorIdList = mapSensorNameAndId.mapSensorNameToId(refSensorNameList)

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
        # testAnalyzeComCamOpdData
        pssn = self.imsimCmpt.getOpdPssnFromFile(self.pssnFileName)
        self.assertEqual(len(pssn), 9)

    def testReorderAndSaveWfErrFile(self):
        listOfWfErr = self._prepareListOfWfErr()

        refSensorNameList = self._getRefSensorNameList()
        zkFileName = "testZk.zer"
        camera = getCamera("lsst")
        self.imsimCmpt.reorderAndSaveWfErrFile(
            listOfWfErr, refSensorNameList, camera, zkFileName=zkFileName
        )

        zkFilePath = os.path.join(self.imsimCmpt.outputImgDir, zkFileName)
        with open(zkFilePath, "r") as zkFile:
            zkInFile = yaml.safe_load(zkFile)
        print(zkInFile)

        numOfZk = self.imsimCmpt.numOfZk
        self.assertEqual(len(zkInFile), len(refSensorNameList))
        self.assertEqual(
            len(
                np.fromstring(
                    zkInFile[camera[refSensorNameList[0]].getId()][0], sep=" "
                )
            ),
            numOfZk,
        )

        self.assertEqual(np.sum(zkInFile[0, :]), 0)
        self.assertEqual(np.sum(zkInFile[3, :]), 0)

        delta = np.sum(np.abs(zkInFile[1, :] - listOfWfErr[2].annularZernikePoly))
        self.assertLess(delta, 1e-10)

        delta = np.sum(np.abs(zkInFile[2, :] - listOfWfErr[1].annularZernikePoly))
        self.assertLess(delta, 1e-10)

    def _prepareListOfWfErr(self):
        numOfZk = self.imsimCmpt.numOfZk

        sensorIdList = [2, 3, 1]
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
