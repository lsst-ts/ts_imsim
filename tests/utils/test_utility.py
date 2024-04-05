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
from lsst.obs.lsst import LsstCam, LsstComCam
from lsst.ts.imsim.utils import (
    CamType,
    ModifiedEnvironment,
    get_camera,
    get_config_dir,
    get_module_path,
    get_policy_path,
    get_zk_from_file,
)


class TestUtility(unittest.TestCase):
    """Test the utility functions."""

    def setUp(self) -> None:
        self.test_data_dir = os.path.join(
            get_module_path(), "tests", "testData", "utils"
        )

    def test_get_module_path(self):
        self.assertEqual(os.environ["TS_IMSIM_DIR"], get_module_path())

    def test_get_policy_path(self):
        self.assertEqual(
            os.path.join(os.environ["TS_IMSIM_DIR"], "policy"), get_policy_path()
        )

    def test_get_config_dir(self):
        self.assertEqual(
            os.path.join(os.environ["TS_IMSIM_DIR"], "policy", "config"),
            get_config_dir(),
        )

    def test_get_camera(self):
        lsst_fam_cam = get_camera(CamType.LsstFamCam)
        self.assertIsInstance(lsst_fam_cam, cameraGeom.Camera)
        self.assertEqual(lsst_fam_cam.getName(), LsstCam.getCamera().getName())

        lsst_cam = get_camera(CamType.LsstCam)
        self.assertIsInstance(lsst_cam, cameraGeom.Camera)
        self.assertEqual(lsst_cam.getName(), LsstCam.getCamera().getName())

        com_cam = get_camera(CamType.ComCam)
        self.assertIsInstance(com_cam, cameraGeom.Camera)
        self.assertEqual(com_cam.getName(), LsstComCam.getCamera().getName())

        with self.assertRaises(ValueError):
            get_camera("invalid")

    def test_get_zk_from_file(self):
        zk_from_file = get_zk_from_file(
            os.path.join(get_module_path(), "tests", "testData", "opd", "opd.zer")
        )
        self.assertCountEqual(zk_from_file.keys(), [191, 195, 199, 203])

    def test_modified_environment(self):
        # Test adding a new environment variable
        self.assertNotIn("TEST_AOS_ROCKS_5772", os.environ)
        with ModifiedEnvironment(TEST_AOS_ROCKS_5772="parrot"):
            self.assertEqual(os.environ["TEST_AOS_ROCKS_5772"], "parrot")
        self.assertNotIn("TEST_AOS_ROCKS_5772", os.environ)

        # Test modifying an existing environment variable
        self.assertNotIn("TEST_AOS_ROCKS_57721", os.environ)
        os.environ["TEST_AOS_ROCKS_57721"] = "spam"
        self.assertEqual(os.environ["TEST_AOS_ROCKS_57721"], "spam")
        with ModifiedEnvironment(TEST_AOS_ROCKS_57721="eggs"):
            self.assertEqual(os.environ["TEST_AOS_ROCKS_57721"], "eggs")
        self.assertEqual(os.environ["TEST_AOS_ROCKS_57721"], "spam")

        # Test unsetting an existing environment variable
        with ModifiedEnvironment(TEST_AOS_ROCKS_57721=None):
            self.assertNotIn("TEST_AOS_ROCKS_57721", os.environ)
        self.assertEqual(os.environ["TEST_AOS_ROCKS_57721"], "spam")

        # Clean up
        del os.environ["TEST_AOS_ROCKS_57721"]
