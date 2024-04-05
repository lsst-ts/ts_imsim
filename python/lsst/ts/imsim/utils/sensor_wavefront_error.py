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

__all__ = ["SensorWavefrontError"]

import numpy as np


class SensorWavefrontError:
    """Contains the wavefront errors for a single sensor."""

    def __init__(self, num_of_zk: int = 19) -> None:
        """Construct a sensor wavefront error.

        Parameters
        ----------
        num_of_zk : int, optional
            Number of annular Zernike polynomials. (the default is 19.)
        """
        # Sensor Id
        self._sensor_id = 999

        # Sensor Name
        self.sensor_name = "R99_S99"

        # Number of zk
        self.num_of_zk = int(num_of_zk)

        # Annular Zernike polynomials (zk)
        self._annular_zernike_poly = np.zeros(self.num_of_zk)

    @property
    def sensor_id(self) -> int:
        return self._sensor_id

    @sensor_id.setter
    def sensor_id(self, new_sensor_id: int) -> None:
        """Set the sensor Id.

        Parameters
        ----------
        new_sensor_id : int
            The Id of the sensor this wavefront error is for.

        Raises
        ------
        ValueError
            sensor_id must be >= 0.
        """
        if new_sensor_id < 0:
            raise ValueError("sensor_id must be >= 0.")
        self._sensor_id = int(new_sensor_id)

    @property
    def annular_zernike_poly(self) -> np.ndarray:
        return self._annular_zernike_poly

    @annular_zernike_poly.setter
    def annular_zernike_poly(self, new_annular_zernike_poly: np.ndarray) -> None:
        """Set the effective annular zernike poly.

        Parameters
        ----------
        new_annular_zernike_poly : numpy.ndarray[self.num_of_zk] (float)
            The poly describing the wavefront error in um.

        Raises
        ------
        ValueError
            annular_zernike_poly must be an array of self.num_of_zk floats.
        """
        if len(new_annular_zernike_poly) != self.num_of_zk:
            raise ValueError(
                "annular_zernike_poly must be an array of %d floats." % self.num_of_zk
            )
        self._annular_zernike_poly = np.array(new_annular_zernike_poly)


if __name__ == "__main__":
    pass
