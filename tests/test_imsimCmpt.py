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

from lsst.ts.imsim.obsMetadata import ObsMetadata
from lsst.ts.imsim.utils.utility import getConfigDir, getModulePath
from lsst.ts.imsim.imsimCmpt import ImsimCmpt


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
