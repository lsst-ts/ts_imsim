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

import galsim
import numpy as np
import yaml
from lsst.afw.cameraGeom import FIELD_ANGLE
from lsst.ts.imsim import ImsimCmpt, ObsMetadata, SkySim
from lsst.ts.imsim.utils import (
    CamType,
    SensorWavefrontError,
    get_camera,
    get_config_dir,
    get_module_path,
    get_zk_from_file,
)


class TestImsimCmpt(unittest.TestCase):
    """Test the imsim configuration module."""

    max_diff = None

    def setUp(self):
        self.obs_metadata_test = ObsMetadata(
            ra=0.0, dec=0.0, band="r", mjd=59580.0, raw_seeing=0.5
        )
        self.config_pointer_default_lsst_cam = os.path.join(
            get_config_dir(), "lsstCamDefaultPointer.yaml"
        )
        with open(
            os.path.join(
                get_module_path(),
                "tests",
                "testData",
                "imsimConfig",
                "imsimConfigLsstCam.yaml",
            )
        ) as test_file:
            self.full_test_yaml = yaml.safe_load(test_file)
        self.imsim_cmpt = ImsimCmpt()

        # Set the output directories
        self.output_dir = os.path.join(get_module_path(), "tests", "tmp")
        self.output_img_dir = os.path.join(self.output_dir, "img")
        self.opd_file_path = os.path.join(
            self._get_opd_file_dir_of_lsst_cam(), "opd.fits"
        )

        self.imsim_cmpt.output_dir = self.output_dir
        self.imsim_cmpt.output_img_dir = self.output_img_dir
        self.imsim_cmpt.opd_file_path = self.opd_file_path

        # Set the file name of analyzed OPD data
        self.zk_file_name = "opd.zer"
        self.pssn_file_name = "PSSN.txt"

        # Create SkySim object
        self.sky_sim = SkySim()

        # Store camera for tests
        self.camera = get_camera(CamType.LsstCam)

    def tearDown(self):
        shutil.rmtree(self.output_dir)

    def _map_sensor_name_and_id(self, sensor_name_list):
        return dict(
            [
                (self.camera[detector].getName(), self.camera[detector].getId())
                for detector in sensor_name_list
            ]
        )

    def test_set_output_dir(self):
        self.assertTrue(os.path.exists(self.imsim_cmpt.output_dir))

    def test_set_output_img_dir(self):
        self.assertTrue(os.path.exists(self.imsim_cmpt.output_img_dir))

    def test_verify_pointer_file_raises_error(self):
        required_keys = ["input", "gal", "image", "psf", "stamp", "output"]
        pointer_without_input = {key: "" for key in required_keys[1:]}
        with self.assertRaises(ValueError) as context:
            self.imsim_cmpt._verify_pointer_file(pointer_without_input, required_keys)
        self.assertEqual(
            str(context.exception), "Config pointer file missing filepath for input."
        )

    def test_assemble_config(self):
        full_config_yaml = self.imsim_cmpt.assemble_config_yaml(
            self.obs_metadata_test,
            self.config_pointer_default_lsst_cam,
            CamType.LsstCam,
        )
        self.full_test_yaml["output"]["dir"] = self.imsim_cmpt.output_img_dir
        self.assertDictEqual(full_config_yaml, self.full_test_yaml)

    def test_convert_obs_metadata_to_text(self):
        obs_variables_text = self.imsim_cmpt.convert_obs_metadata_to_text(
            self.obs_metadata_test
        )
        with open(os.path.join(get_config_dir(), "obsVariablesDefault.yaml")) as file:
            default_vars = yaml.safe_load(file)
        self.assertDictEqual(yaml.safe_load(obs_variables_text), default_vars)

    def test_format_opd_text_lsst_cam(self):
        opd_text = self.imsim_cmpt.format_opd_text(
            self.obs_metadata_test, CamType.LsstCam
        )

        # Should match values from obs_lsst
        test_intra_location = self.camera[204].getCenter(FIELD_ANGLE)
        test_extra_location = self.camera[203].getCenter(FIELD_ANGLE)
        test_y = np.degrees(test_intra_location[0] + test_extra_location[0]) / 2
        test_x = np.degrees(test_intra_location[1] + test_extra_location[1]) / 2

        self.assertTrue(opd_text.startswith("  opd:"))
        self.assertTrue(
            opd_text.endswith(f"- {{thx: {test_x} deg, thy: {test_y} deg}}\n")
        )

    def test_add_config_header(self):
        obs_info_text = self.imsim_cmpt.convert_obs_metadata_to_text(
            self.obs_metadata_test
        )
        header_text = self.imsim_cmpt.add_config_header(self.obs_metadata_test)
        header_yaml = yaml.safe_load(
            str(obs_info_text + "\n" + "output:\n" + header_text)
        )
        for key in list(header_yaml["output"]["header"].keys()):
            self.assertEqual(
                header_yaml["output"]["header"][key],
                self.full_test_yaml["output"]["header"][key],
            )

    def test_gen_instance_catalog(self):
        self.sky_sim.add_star_by_file(
            os.path.join(
                get_module_path(), "tests", "testData", "sky", "wfsSglStar.txt"
            )
        )
        inst_cat = self.imsim_cmpt.gen_instance_catalog(self.sky_sim)
        expected_cat = "object  0\t 1.196000\t 1.176000 17.000000 "
        expected_cat += (
            "flatSED/sed_flat.txt.gz 0.0 0.0 0.0 0.0 0.0 0.0 point none none \n"
        )
        self.assertEqual(inst_cat, expected_cat)

    def test_gen_inst_cat_stars(self):
        self.sky_sim.add_star_by_file(
            os.path.join(
                get_module_path(), "tests", "testData", "sky", "wfsSglStar.txt"
            )
        )
        inst_cat = self.imsim_cmpt.gen_inst_cat_stars(self.sky_sim)
        expected_cat = "object  0\t 1.196000\t 1.176000 17.000000 "
        expected_cat += (
            "flatSED/sed_flat.txt.gz 0.0 0.0 0.0 0.0 0.0 0.0 point none none \n"
        )
        self.assertEqual(inst_cat, expected_cat)

    def test_generate_star(self):
        star_id = 0
        ra = 1.0
        dec = 1.0
        mag_norm = 2.0
        sed_name = "flat.txt"
        content = self.imsim_cmpt.generate_star(star_id, ra, dec, mag_norm, sed_name)
        ans_content = "object  0\t 1.000000\t 1.000000  2.000000 "
        ans_content += "flat.txt 0.0 0.0 0.0 0.0 0.0 0.0 point "
        ans_content += "none none \n"
        self.assertEqual(content, ans_content)

    def test_analyze_opd_data(self):
        self._analyze_lsst_cam_opd_data()

        zk_file_path = os.path.join(self.output_img_dir, self.zk_file_name)
        pssn_file_path = os.path.join(self.output_img_dir, self.pssn_file_name)
        self.assertTrue(os.path.exists(zk_file_path))
        self.assertTrue(os.path.exists(pssn_file_path))

        zk = get_zk_from_file(zk_file_path)
        ans_zk_file_path = os.path.join(self._get_opd_file_dir_of_lsst_cam(), "opd.zer")
        ans_zk = get_zk_from_file(ans_zk_file_path)

        for det in [191, 195, 199, 203]:
            delta = np.sum(np.abs(zk[det] - ans_zk[det]))
            self.assertLess(delta, 1e-10)

        pssn_data = np.loadtxt(pssn_file_path)
        pssn = pssn_data[0, :]
        ans_pssn_file_path = os.path.join(
            self._get_opd_file_dir_of_lsst_cam(), "PSSN.txt"
        )
        ans_pssn_data = np.loadtxt(ans_pssn_file_path)
        ans_pssn = ans_pssn_data[0, :]

        delta = np.sum(np.abs(pssn - ans_pssn))
        self.assertLess(delta, 1e-10)

    def test_rotate_opd_data(self):
        """Test the rotation of OPD data."""
        rotation_angle = 20.0
        num_opd = 4

        zerinkes_rotation_0 = (
            self.imsim_cmpt._map_opd_to_zk(
                rot_opd_in_deg=0.0,
                num_opd=num_opd,
            )
            * 1e-3
        )

        zerinkes_rotation_20 = (
            self.imsim_cmpt._map_opd_to_zk(
                rot_opd_in_deg=rotation_angle,
                num_opd=num_opd,
            )
            * 1e-3
        )

        for idx in range(num_opd):
            padded_zernikes = np.pad(zerinkes_rotation_20[idx, :], (4, 0))

            galsim_zernikes = galsim.zernike.Zernike(
                padded_zernikes,
                R_outer=4.18,  # Outer radius obscuration
                R_inner=2.558,  # Inner radius obscuration
            )

            np.testing.assert_almost_equal(
                galsim_zernikes.rotate(np.deg2rad(rotation_angle)).coef[4:],
                zerinkes_rotation_0[idx, :],
                decimal=1,
            )

    def _analyze_lsst_cam_opd_data(self, rot_opd_in_deg=0.0):
        shutil.copy(
            os.path.join(self._get_opd_file_dir_of_lsst_cam(), "opd.fits"),
            self.imsim_cmpt.output_img_dir,
        )
        self.imsim_cmpt.analyze_opd_data(
            CamType.LsstCam,
            zk_file_name=self.zk_file_name,
            rot_opd_in_deg=rot_opd_in_deg,
            pssn_file_name=self.pssn_file_name,
        )

    def _get_opd_file_dir_of_lsst_cam(self):
        opd_file_dir = os.path.join(get_module_path(), "tests", "testData", "opd")

        return opd_file_dir

    def test_get_opd_gq_eff_fwhm_from_file(self):
        self._analyze_lsst_cam_opd_data()

        gq_eff_fwhm = self.imsim_cmpt.get_opd_gq_eff_fwhm_from_file(self.pssn_file_name)
        ans_pssn_file_path = os.path.join(
            self._get_opd_file_dir_of_lsst_cam(), "PSSN.txt"
        )
        ans_pssn_data = np.loadtxt(ans_pssn_file_path)
        ans_gq_eff_fwhm = ans_pssn_data[1, -1]
        self.assertAlmostEqual(gq_eff_fwhm, ans_gq_eff_fwhm, places=3)

    def _get_ref_sensor_name_list(self):
        ref_sensor_name_list = [
            "R00_SW0",
            "R44_SW0",
            "R04_SW0",
            "R40_SW0",
        ]

        return ref_sensor_name_list

    def test_map_opd_data_to_list_of_wf_err(self):
        self._analyze_lsst_cam_opd_data()

        ref_sensor_name_list = self._get_ref_sensor_name_list()
        map_sensor_name_and_id = self._map_sensor_name_and_id(ref_sensor_name_list)
        ans_sensor_id_list = list(map_sensor_name_and_id.values())

        list_of_wf_err = self.imsim_cmpt.map_opd_data_to_list_of_wf_err(
            self.zk_file_name, ans_sensor_id_list, ref_sensor_name_list
        )

        self.assertEqual(len(list_of_wf_err), len(ref_sensor_name_list))

        opd_zk = get_zk_from_file(
            os.path.join(self.imsim_cmpt.output_img_dir, self.zk_file_name)
        )
        map_sensor_name_and_id = self._map_sensor_name_and_id(ref_sensor_name_list)
        for wf_err, ref_sensor_name in zip(list_of_wf_err, ref_sensor_name_list):
            sensor_id = wf_err.sensor_id
            sensor_name = wf_err.sensor_name
            self.assertEqual(sensor_name, ref_sensor_name)
            self.assertEqual(sensor_id, map_sensor_name_and_id[ref_sensor_name])

            zk_in_wf_err = wf_err.annular_zernike_poly
            # imSim outputs OPD in nanometers so need to change to microns
            # to be consistent with what WEP expects.
            delta = np.sum(np.abs(zk_in_wf_err - opd_zk[sensor_id] / 1e3))
            self.assertEqual(delta, 0)

    def test_get_list_of_fwhm_sensor_data(self):
        self._analyze_lsst_cam_opd_data()

        (sensor_data_fwhm) = self.imsim_cmpt.get_list_of_fwhm_sensor_data(
            self.pssn_file_name
        )

        ans_data = self.imsim_cmpt._get_data_of_pssn_file(self.pssn_file_name)
        ans_fwhm_data = ans_data[1, :-1]

        for data_fwhm, ans_fwhm in zip(sensor_data_fwhm, ans_fwhm_data):
            self.assertEqual(data_fwhm, ans_fwhm)

    def test_get_opd_pssn_from_file(self):
        self._analyze_lsst_cam_opd_data()

        # The correctness of values have been tested at the test case of
        # testAnalyzeOpdData
        pssn = self.imsim_cmpt.get_opd_pssn_from_file(self.pssn_file_name)
        self.assertEqual(len(pssn), 4)

    def test_reorder_and_save_wf_err_file(self):
        list_of_wf_err = self._prepare_list_of_wf_err()

        ref_sensor_name_list = self._get_ref_sensor_name_list()
        map_sensor_name_and_id = self._map_sensor_name_and_id(ref_sensor_name_list)
        sensor_id_list = list(map_sensor_name_and_id.values())

        zk_file_name = "testZk.zer"
        self.imsim_cmpt.reorder_and_save_wf_err_file(
            list_of_wf_err, ref_sensor_name_list, self.camera, zk_file_name=zk_file_name
        )

        zk_file_path = os.path.join(self.imsim_cmpt.output_img_dir, zk_file_name)
        zk_in_file = get_zk_from_file(zk_file_path)

        num_of_zk = self.imsim_cmpt.num_of_zk
        self.assertEqual(len(zk_in_file), len(ref_sensor_name_list))
        self.assertEqual(
            len(zk_in_file[sensor_id_list[0]]),
            num_of_zk,
        )

        self.assertEqual(np.sum(zk_in_file[191]), 0)
        self.assertEqual(np.sum(zk_in_file[203]), 0)

        delta = np.sum(np.abs(zk_in_file[195] - list_of_wf_err[0].annular_zernike_poly))
        self.assertLess(delta, 1e-7)

        delta = np.sum(np.abs(zk_in_file[199] - list_of_wf_err[1].annular_zernike_poly))
        self.assertLess(delta, 1e-7)

    def _prepare_list_of_wf_err(self):
        num_of_zk = self.imsim_cmpt.num_of_zk

        sensor_id_list = [195, 199]
        list_of_wf_err = []
        for sensor_id in sensor_id_list:
            sensor_wavefront_data = SensorWavefrontError()
            sensor_wavefront_data.sensor_id = sensor_id

            wf_err = np.random.rand(num_of_zk)
            sensor_wavefront_data.annular_zernike_poly = wf_err

            list_of_wf_err.append(sensor_wavefront_data)

        return list_of_wf_err

    def test_save_dof_in_um_file_for_next_iter(self):
        self.imsim_cmpt.dof_in_um = np.arange(50)
        dof_in_um_file_name = "dofPertInNextIter.mat"
        self.imsim_cmpt.save_dof_in_um_file_for_next_iter(
            dof_in_um_file_name=dof_in_um_file_name
        )

        filePath = os.path.join(self.output_dir, dof_in_um_file_name)
        self.assertTrue(os.path.exists(filePath))

        data = np.loadtxt(filePath)
        delta = np.sum(np.abs(self.imsim_cmpt.dof_in_um - data))
        self.assertEqual(delta, 0)
