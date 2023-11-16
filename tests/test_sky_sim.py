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

from lsst.ts.imsim import SkySim
from lsst.ts.imsim.utils.utility import get_module_path
from lsst.ts.wep.utils import CamType


class TestSkySim(unittest.TestCase):
    def setUp(self):
        self.sky_sim = SkySim()

    def test_set_camera(self):
        self.sky_sim.set_camera(CamType.LsstFamCam)
        self.assertEqual(self.sky_sim._camera.getName(), "LSSTCam")

    def test_add_star_by_ra_dec_in_deg(self):
        self.sky_sim.add_star_by_ra_dec_in_deg(1, 2, 3, 4)
        self.assertEqual(len(self.sky_sim.star_id), 1)

        self.sky_sim.add_star_by_ra_dec_in_deg(2, 2.1, 3, 4)
        self.assertEqual(len(self.sky_sim.star_id), 2)
        self.assertEqual(self.sky_sim.star_id[0], 1)
        self.assertEqual(self.sky_sim.star_id[1], 2)

        # Try to add the same star Id again
        self.sky_sim.add_star_by_ra_dec_in_deg(2, 2.1, 3, 4)
        self.assertEqual(len(self.sky_sim.star_id), 2)

    def test_add_star_by_file(self):
        self._add_star_by_file("wfsStar.txt")

        self.assertEqual(len(self.sky_sim.star_id), 8)

        ra = self.sky_sim.ra
        dec = self.sky_sim.dec
        self.assertEqual(ra[2], -1.176)
        self.assertEqual(dec[2], 1.216)

        self.assertEqual(self.sky_sim.mag[2], 15.0)

    def _add_star_by_file(self, sky_file_name):
        skyFile = os.path.join(
            get_module_path(), "tests", "testData", "sky", sky_file_name
        )
        self.sky_sim.add_star_by_file(skyFile)

    def test_add_star_by_file_with_sgl_star(self):
        self._add_star_by_file("wfsSglStar.txt")

        self.assertEqual(len(self.sky_sim.star_id), 1)

        ra = self.sky_sim.ra
        dec = self.sky_sim.dec
        self.assertEqual(ra[0], 1.196)
        self.assertEqual(dec[0], 1.176)

        self.assertEqual(self.sky_sim.mag[0], 17.0)

    def testExportSkyToFile(self):
        self._add_star_by_file("wfsStar.txt")
        output_file_path = os.path.join(
            get_module_path(), "output", "testSkyOutput.txt"
        )

        self.sky_sim.export_sky_to_file(output_file_path)
        self.assertTrue(os.path.isfile(output_file_path))
        os.remove(output_file_path)
