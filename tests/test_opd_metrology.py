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

import lsst.obs.lsst as obs_lsst
import numpy as np
from astropy.io import fits
from lsst.afw.cameraGeom import FIELD_ANGLE
from lsst.ts.imsim import OpdMetrology
from lsst.ts.imsim.utils.utility import get_module_path, get_zk_from_file
from lsst.ts.wep.utils import CamType


class TestOpdMetrology(unittest.TestCase):
    """Test the OpdMetrology class."""

    def setUp(self):
        self.metr = OpdMetrology()
        self.test_data_dir = os.path.join(get_module_path(), "tests", "testData")

    def test_set_weighting_ratio(self):
        wt = [1, 2]
        self.metr.wt = wt

        wt_in_metr = self.metr.wt
        self.assertEqual(len(wt_in_metr), len(wt))
        self.assertEqual(np.sum(wt_in_metr), 1)
        self.assertAlmostEqual(wt_in_metr[1] / wt_in_metr[0], 2)
        with self.assertRaises(ValueError) as context:
            self.metr.wt = [-1, 1]
        self.assertEqual(str(context.exception), "All weighting ratios should be >= 0.")

    def test_set_wgt_and_field_xy_of_gq_lsst_fam(self):
        self.metr.set_wgt_and_field_xy_of_gq(CamType.LsstFamCam)

        field_x = self.metr.field_x
        self.assertEqual(len(field_x), 189)

        wgt = self.metr.wt
        self.assertEqual(len(wgt), 189)

    def test_set_wgt_and_field_xy_of_gq_com_cam(self):
        self.metr.set_wgt_and_field_xy_of_gq(CamType.ComCam)

        field_x = self.metr.field_x
        self.assertEqual(len(field_x), 9)

        wgt = self.metr.wt
        self.assertEqual(len(wgt), 9)

    def test_set_wgt_and_field_xy_of_gq_err(self):
        self.assertRaises(
            ValueError, self.metr.set_wgt_and_field_xy_of_gq, "NoThisInstName"
        )

    def test_set_default_lsst_wfs_gq(self):
        self.metr.set_default_lsst_wfs_gq()

        field_x = self.metr.field_x
        field_y = self.metr.field_y
        # Get true values from obs_lsst
        camera = obs_lsst.LsstCam.getCamera()
        det_id_map = camera.getIdMap()
        true_x = []
        true_y = []
        for sens_id in [192, 196, 200, 204]:
            center_intra = det_id_map[sens_id].getCenter(FIELD_ANGLE)
            center_extra = det_id_map[sens_id - 1].getCenter(FIELD_ANGLE)
            center_det = np.degrees(np.array([center_intra, center_extra]))
            center_x = np.mean(center_det[:, 1])
            center_y = np.mean(center_det[:, 0])
            true_x.append(center_x)
            true_y.append(center_y)

        self.assertCountEqual(field_x, true_x)
        self.assertCountEqual(
            field_y,
            true_y,
        )

        wgt = self.metr.wt
        self.assertCountEqual(wgt, [0.25, 0.25, 0.25, 0.25])

    def test_get_default_lsst_wfs_gq(self):
        field_wfs_x, field_wfs_y, det_ids = self.metr.get_default_lsst_wfs_gq()
        self.assertEqual(len(field_wfs_x), 4)
        self.assertEqual(len(det_ids), 4)

    def _get_opd_dir(self):
        opd_file_dir = os.path.join(self.test_data_dir, "opd")
        return opd_file_dir

    def test_get_zk_from_opd(self):
        opd_dir = self._get_opd_dir()
        zk = self.metr.get_zk_from_opd(opd_fits_file=os.path.join(opd_dir, "opd.fits"))[
            0
        ]

        ans_opd_file_name = "opd.zer"
        ans_opd_file_path = os.path.join(opd_dir, ans_opd_file_name)
        all_opd_ans = get_zk_from_file(ans_opd_file_path)
        self.assertLess(np.sum(np.abs(zk[3:] - all_opd_ans[191])), 1e-5)

    def test_rm_piston_tip_tilt_from_opd(self):
        """Test removal of piston (z1), x-tilt (z2), and y-tilt (z3)
        from the OPD map."""
        opd_dir = self._get_opd_dir()
        opd_file_path = os.path.join(opd_dir, "opd.fits")
        opd_map = fits.getdata(opd_file_path, 1)
        opd_rm_ptt, opd_x, opd_y = self.metr.rm_piston_tip_tilt_from_opd(
            opd_map=opd_map
        )

        # Flip OPD because it will be flipped inside get_zk_from_opd
        zk_rm_ptt = self.metr.get_zk_from_opd(opd_map=opd_rm_ptt)[0]
        zk_rm_ptt_in_um = np.sum(np.abs(zk_rm_ptt[0:3])) / 1e3
        self.assertLess(zk_rm_ptt_in_um, 9e-2)

    def test_calc_pssn(self):
        pssn = self._calc_pssn()
        all_data = self._get_metro_all_ans_data()
        self.assertAlmostEqual(pssn, all_data[0, 0])

    def _calc_pssn(self):
        wavelength_in_um = 0.48
        opd_file_path = os.path.join(self._get_opd_dir(), "opd.fits")
        opd_map = fits.getdata(opd_file_path, 0) * 1e-3
        pssn = self.metr.calc_pssn(wavelength_in_um, opd_map=opd_map)

        return pssn

    def _get_metro_all_ans_data(self):
        ans_all_data_file_name = "PSSN.txt"
        ans_all_data_file_path = os.path.join(
            self._get_opd_dir(), ans_all_data_file_name
        )
        all_data = np.loadtxt(ans_all_data_file_path)

        return all_data

    def test_calc_fwhm_eff(self):
        pssn = self._calc_pssn()
        fwhm = self.metr.calc_fwhm_eff(pssn)

        all_data = self._get_metro_all_ans_data()
        self.assertAlmostEqual(fwhm, all_data[1, 0])

    def test_calc_gq_value(self):
        self.metr.set_default_lsst_wfs_gq()
        all_data = self._get_metro_all_ans_data()
        value_list = all_data[0, 0:4]

        gq_value = self.metr.calc_gq_value(value_list)
        self.assertAlmostEqual(gq_value, all_data[0, -1])
