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

from lsst.afw import cameraGeom
from lsst.obs.lsst import LsstCam
from lsst.ts.imsim.utils.utility import (
    getCamera,
    getConfigDir,
    getModulePath,
    getPolicyPath,
    getZkFromFile,
)


class TestUtility(unittest.TestCase):
    def testGetModulePath(self):
        self.assertEqual(os.environ["TS_IMSIM_DIR"], getModulePath())

    def testGetPolicyPath(self):
        self.assertEqual(
            os.path.join(os.environ["TS_IMSIM_DIR"], "policy"), getPolicyPath()
        )

    def testGetConfigDir(self):
        self.assertEqual(
            os.path.join(os.environ["TS_IMSIM_DIR"], "policy", "config"), getConfigDir()
        )

    def testGetCamera(self):
        lsstFamCam = getCamera("lsstfam")
        self.assertIsInstance(lsstFamCam, cameraGeom.Camera)
        self.assertEqual(lsstFamCam.getName(), LsstCam.getCamera().getName())

        lsstCam = getCamera("lsst")
        self.assertIsInstance(lsstCam, cameraGeom.Camera)
        self.assertEqual(lsstCam.getName(), LsstCam.getCamera().getName())

        with self.assertRaises(ValueError):
            getCamera("invalid")

    def testGetZkFromFile(self):
        zkFromFile = getZkFromFile(
            os.path.join(getModulePath(), "tests", "testData", "opd", "opd.zer")
        )
        self.assertCountEqual(zkFromFile.keys(), [191, 195, 199, 203])
