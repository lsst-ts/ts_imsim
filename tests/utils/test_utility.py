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
from galsim.zernike import zernikeBasis
from lsst.afw import cameraGeom
from lsst.obs.lsst import LsstCam, LsstComCam
from lsst.ts.imsim.utils import (
    CamType,
    ModifiedEnvironment,
    get_cam_type,
    get_camera,
    get_config_dir,
    get_module_path,
    get_policy_path,
    get_zk_from_file,
    zernike_annular_fit,
    zernike_eval,
)


class TestUtility(unittest.TestCase):
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

    def test_get_cam_type(self):
        self.assertEqual(get_cam_type("lsst"), CamType.LsstCam)
        self.assertEqual(get_cam_type("lsstfam"), CamType.LsstFamCam)
        self.assertEqual(get_cam_type("comcam"), CamType.ComCam)
        instName = "telescope"
        assertMsg = f"Instrument name ({instName}) is not supported."
        with self.assertRaises(ValueError) as context:
            get_cam_type(instName)
        self.assertTrue(assertMsg in str(context.exception))

    def test_zernike_annular_fit(self):
        opdFitsFile = os.path.join(self.test_data_dir, "sim6_iter0_opd0.fits.gz")
        opd = fits.getdata(opdFitsFile)

        # x-, y-coordinate in the OPD image
        opdSize = opd.shape[0]
        opdGrid1d = np.linspace(-1, 1, opdSize)
        opdx, opdy = np.meshgrid(opdGrid1d, opdGrid1d)

        idx = opd != 0
        coef = zernike_annular_fit(opd[idx], opdx[idx], opdy[idx], 22, 0.61)

        ansOpdFileName = "sim6_iter0_opd.zer"
        ansOpdFilePath = os.path.join(self.test_data_dir, ansOpdFileName)
        allOpdAns = np.loadtxt(ansOpdFilePath)
        self.assertLess(np.sum(np.abs(coef - allOpdAns[0, :])), 1e-10)

    def test_zernike_eval(self):
        # Create a Zernike basis
        grid = np.linspace(-1, 1, 200)
        uGrid, vGrid = np.meshgrid(grid, grid)
        zkBasis = zernikeBasis(22, uGrid, vGrid, R_inner=0.61)[4:]

        # Evaluate each Zernike polynomial and compare to basis
        for i in range(len(zkBasis)):
            coeff = np.zeros(i + 1)
            coeff[-1] = 1
            z_eval = zernike_eval(coeff, uGrid, vGrid)
            np.allclose(z_eval, zkBasis[i])
