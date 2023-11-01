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

import unittest

import numpy as np
from lsst.ts.imsim.utils.sensorWavefrontError import SensorWavefrontError


class TestSensorWavefrontError(unittest.TestCase):
    """Test the SensorWavefrontError class."""

    def setUp(self):
        self.num_of_zk = 19
        self.sensor_wavefront_error = SensorWavefrontError(num_of_zk=self.num_of_zk)

        self.sensor_id = 999
        self.sensor_name = "R99_S99"

    def test_get_num_of_zk(self):
        self.assertEqual(self.sensor_wavefront_error.num_of_zk, self.num_of_zk)

    def test_get_sensor_id(self):
        sensor_id = self.sensor_wavefront_error.sensor_id
        self.assertEqual(sensor_id, self.sensor_id)
        self.assertTrue(isinstance(sensor_id, int))

    def testSetSensorId(self):
        sensor_id = 2
        self.sensor_wavefront_error.sensor_id = sensor_id

        sensor_id_in_object = self.sensor_wavefront_error.sensor_id
        self.assertEqual(sensor_id_in_object, sensor_id)

    def test_set_sensor_id_with_float_input(self):
        sensor_id = 2.0
        self.sensor_wavefront_error.sensor_id = sensor_id

        sensor_id_in_object = self.sensor_wavefront_error.sensor_id
        self.assertTrue(isinstance(sensor_id_in_object, int))
        self.assertEqual(sensor_id_in_object, sensor_id)

    def test_set_sensor_id_with_value_less_than_zero(self):
        with self.assertRaises(ValueError) as context:
            self.sensor_wavefront_error.sensor_id = -1
        self.assertEqual(str(context.exception), "sensor_id must be >= 0.")

    def test_get_sensor_name(self):
        sensor_name = self.sensor_wavefront_error.sensor_name
        self.assertEqual(sensor_name, self.sensor_name)
        self.assertTrue(isinstance(sensor_name, str))

    def test_set_sensor_name(self):
        sensor_name = "R42_S24"
        self.sensor_wavefront_error.sensor_name = sensor_name

        sensor_name_in_object = self.sensor_wavefront_error.sensor_name
        self.assertEqual(sensor_name_in_object, sensor_name)

    def test_get_annulary_zernike_poly(self):
        annular_zernike_poly = self.sensor_wavefront_error.annular_zernike_poly

        self.assertEqual(len(annular_zernike_poly), self.num_of_zk)
        self.assertTrue(isinstance(annular_zernike_poly, np.ndarray))

        delta = np.sum(np.abs(annular_zernike_poly))
        self.assertEqual(delta, 0)

    def test_set_annular_zernike_poly(self):
        rand_value = np.random.rand(self.num_of_zk)
        self.sensor_wavefront_error.annular_zernike_poly = rand_value

        value_in_obj = self.sensor_wavefront_error.annular_zernike_poly

        delta = np.sum(np.abs(rand_value - value_in_obj))
        self.assertEqual(delta, 0)

    def test_set_annular_zernike_poly_with_list_input(self):
        list_value = [1] * self.num_of_zk
        self.sensor_wavefront_error.annular_zernike_poly = list_value

        value_in_obj = self.sensor_wavefront_error.annular_zernike_poly
        self.assertEqual(np.sum(value_in_obj), self.num_of_zk)

    def test_set_annular_zernike_poly_with_wrong_length(self):
        wrong_value = np.ones(self.num_of_zk + 1)
        with self.assertRaises(ValueError) as context:
            self.sensor_wavefront_error.annular_zernike_poly = wrong_value
        self.assertEqual(
            str(context.exception),
            f"annular_zernike_poly must be an array of {self.num_of_zk} floats.",
        )


if __name__ == "__main__":
    # Run the unit test
    unittest.main()
