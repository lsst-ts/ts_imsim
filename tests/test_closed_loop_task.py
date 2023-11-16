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

import argparse
import os
import tempfile
import unittest

from lsst.ts.imsim import ClosedLoopTask, ObsMetadata
from lsst.ts.imsim.utils.utility import get_module_path
from lsst.ts.wep.utils import CamType, FilterType


class TestclosedLoopTask(unittest.TestCase):
    """Test the closedLoopTask class."""

    def setUp(self):
        self.closed_loop_task = ClosedLoopTask()
        self.obs_metadata = ObsMetadata(0.0, 0.0, "r")

        root_test_dir = os.path.join(get_module_path(), "tests")
        self.test_dir = tempfile.TemporaryDirectory(dir=root_test_dir)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_config_sky_sim_with_error(self):
        self.assertRaises(
            ValueError,
            self.closed_loop_task.config_sky_sim,
            "NoThisInstName",
            self.obs_metadata,
        )

    def test_config_sky_sim_no_sky_file_lsst(self):
        self.closed_loop_task.config_sky_sim(CamType.LsstCam, self.obs_metadata)

        sky_sim = self.closed_loop_task.sky_sim
        self.assertEqual(len(sky_sim.star_id), 8)

    def test_config_sky_sim_with_sky_file(self):
        test_sky_file = os.path.join(
            get_module_path(), "tests", "testData", "sky", "wfsStar.txt"
        )
        self.closed_loop_task.config_sky_sim(
            CamType.LsstCam, self.obs_metadata, path_sky_file=test_sky_file
        )

        sky_sim = self.closed_loop_task.sky_sim
        self.assertEqual(len(sky_sim.star_id), 8)
        self.assertEqual(sky_sim.mag[0], 15)

    def test_config_ofc_calc(self):
        cam_type = CamType.LsstCam
        self.closed_loop_task.config_ofc_calc(cam_type)

        ofc_calc = self.closed_loop_task.ofc_calc
        self.assertEqual(ofc_calc.ofc_data.name, "lsst")

    def test_map_filter_ref_to_g(self):
        # test that the reference filter
        # gets mapped to g
        for filter_type_name in ["ref", ""]:
            mapped_filter_name = self.closed_loop_task.map_filter_ref_to_g(
                filter_type_name
            )
            self.assertEqual(mapped_filter_name, "g")

        # test that all other filters are
        # mapped to themselves
        for filter_type_name in "ugrizy":
            mapped_filter_name = self.closed_loop_task.map_filter_ref_to_g(
                filter_type_name
            )
            self.assertEqual(mapped_filter_name, filter_type_name)

    def test_get_filter_type_ref(self):
        filter_type = self.closed_loop_task.get_filter_type("ref")
        self.assertEqual(filter_type, FilterType.REF)

    def test_get_filter_type_u(self):
        filter_type = self.closed_loop_task.get_filter_type("u")
        self.assertEqual(filter_type, FilterType.LSST_U)

    def test_get_filter_type_g(self):
        filter_type = self.closed_loop_task.get_filter_type("g")
        self.assertEqual(filter_type, FilterType.LSST_G)

    def test_get_filter_type_r(self):
        filter_type = self.closed_loop_task.get_filter_type("r")
        self.assertEqual(filter_type, FilterType.LSST_R)

    def test_get_filter_type_i(self):
        filter_type = self.closed_loop_task.get_filter_type("i")
        self.assertEqual(filter_type, FilterType.LSST_I)

    def test_get_filter_type_z(self):
        filter_type = self.closed_loop_task.get_filter_type("z")
        self.assertEqual(filter_type, FilterType.LSST_Z)

    def test_get_filter_type_y(self):
        filter_type = self.closed_loop_task.get_filter_type("y")
        self.assertEqual(filter_type, FilterType.LSST_Y)

    def test_get_filter_type_err(self):
        self.assertRaises(
            ValueError, self.closed_loop_task.get_filter_type, "noThisFilterType"
        )

    def test_erase_directory_content(self):
        # Make the temporary directory
        temp_dir = os.path.join(self.test_dir.name, "tempDir")
        os.mkdir(temp_dir)
        files = os.listdir(self.test_dir.name)
        self.assertEqual(len(files), 1)

        # Try to erase the content
        self.closed_loop_task.erase_directory_content(self.test_dir.name)

        files = os.listdir(self.test_dir.name)
        self.assertEqual(len(files), 0)

    def test_check_boresight(self):
        self.assertRaises(ValueError, self.closed_loop_task.check_boresight, [-1, 0])
        self.assertRaises(ValueError, self.closed_loop_task.check_boresight, [361, 0])
        self.assertRaises(ValueError, self.closed_loop_task.check_boresight, [0, -91])
        self.assertRaises(ValueError, self.closed_loop_task.check_boresight, [0, 91])

    def test_check_and_create_base_output_dir(self):
        # check that the output dir is created
        # where the name is given
        base_output_dir = os.path.join(self.test_dir.name, "testBaseOutputDir")
        self.assertFalse(os.path.exists(base_output_dir))
        self.closed_loop_task.check_and_create_base_output_dir(base_output_dir)
        self.assertTrue(os.path.exists(base_output_dir))

    def test_set_default_parser(self):
        parser = argparse.ArgumentParser()
        parser = ClosedLoopTask.set_default_parser(parser)

        args = parser.parse_known_args()[0]
        self.assertEqual(args.inst, "lsst")
        self.assertEqual(args.filter_type, "")
        self.assertEqual(args.rot_cam, 0.0)
        self.assertEqual(args.iter_num, 5)
        self.assertEqual(args.output, "")
        self.assertFalse(args.clobber)
        self.assertEqual(args.config_pointer_file, "")
        self.assertEqual(args.pipeline_file, "")
        self.assertEqual(args.sky_seed, 42)
        self.assertEqual(args.pert_seed, 11)

    def test_set_img_parser(self):
        parser = argparse.ArgumentParser()
        parser = ClosedLoopTask.set_img_parser(parser)

        args_to_test = ["--boresight_deg", "1.2", "2.3"]

        args = parser.parse_known_args(args=args_to_test)[0]
        self.assertEqual(args.boresight_deg, [1.2, 2.3])
        self.assertEqual(args.sky_file, "")
        self.assertEqual(args.mjd, 60115.33)
        self.assertEqual(args.raw_seeing, 0.5)
        self.assertFalse(args.turn_off_sky_background)
        self.assertFalse(args.turn_off_atmosphere)
        self.assertEqual(args.star_mag, 15.0)

    def test_get_sensor_name_list_of_field_lsst_wfs(self):
        sensor_name_list = self.closed_loop_task.get_sensor_name_list_of_fields(
            CamType.LsstCam
        )
        self.assertEqual(len(sensor_name_list), 8)

        sensor_name_list_ans = [
            "R00_SW0",
            "R00_SW1",
            "R40_SW0",
            "R40_SW1",
            "R04_SW0",
            "R04_SW1",
            "R44_SW0",
            "R44_SW1",
        ]
        self.assertCountEqual(sensor_name_list, sensor_name_list_ans)

    def test_get_sensor_id_list_of_fields(self):
        sensor_id_list = self.closed_loop_task.get_sensor_id_list_of_fields(
            CamType.LsstCam
        )
        self.assertEqual(len(sensor_id_list), 8)

        sensor_id_list_ans = [191, 192, 195, 196, 199, 200, 203, 204]
        self.assertCountEqual(sensor_id_list, sensor_id_list_ans)
