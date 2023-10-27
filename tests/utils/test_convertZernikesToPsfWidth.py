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
from lsst.ts.imsim.utils.convertZernikesToPsfWidth import (
    convertZernikesToPsfWidth,
    getPsfGradPerZernike,
)
from lsst.ts.imsim.utils.utility import getModulePath
from lsst.ts.wep.cwfs.instrument import Instrument
from lsst.ts.wep.utils import CamType


class TestConvertZernikesToPsfWidth(unittest.TestCase):
    """Test the convertZernikesToPsfWidth function."""

    def setUp(self):
        """Load saved conversion factors."""
        self.testDataDir = os.path.join(
            getModulePath(), "tests", "testData", "psfGradientsPerZernike"
        )

    def testLsstCam(self):
        """Test that the LsstCam values match the expected values."""
        # Calculate and compare conversion factors
        conversion_factors = getPsfGradPerZernike(jmax=37)
        expected_factors = np.genfromtxt(os.path.join(self.testDataDir, "lsstcam.txt"))
        self.assertTrue(np.allclose(conversion_factors, expected_factors, atol=1e-3))

        # Perform the comparison using convertZernikesToPsfWidth
        converted_amplitudes = convertZernikesToPsfWidth(np.ones(34))
        self.assertTrue(np.allclose(converted_amplitudes, expected_factors, atol=1e-3))

    def testAuxTel(self):
        """Test that the AuxTel values match the expected values.

        For AuxTel, we will use jmin=1.
        """
        # Setup the AuxTel instrument
        # Note the donut dimension shouldn't matter for this computation
        inst = Instrument()
        inst.configFromFile(160, CamType.AuxTel)
        R_outer = inst.apertureDiameter / 2
        R_inner = R_outer * inst.obscuration

        # Calculate and compare conversion factors
        conversion_factors = getPsfGradPerZernike(
            jmin=1,
            jmax=37,
            R_outer=R_outer,
            R_inner=R_inner,
        )
        expected_factors = np.genfromtxt(os.path.join(self.testDataDir, "auxtel.txt"))
        self.assertTrue(np.allclose(conversion_factors, expected_factors, atol=1e-3))

        # Perform the comparison using convertZernikesToPsfWidth
        converted_amplitudes = convertZernikesToPsfWidth(
            np.ones(37),
            jmin=1,
            R_outer=R_outer,
            R_inner=R_inner,
        )
        self.assertTrue(np.allclose(converted_amplitudes, expected_factors, atol=1e-3))


if __name__ == "__main__":
    # Run the unit test
    unittest.main()