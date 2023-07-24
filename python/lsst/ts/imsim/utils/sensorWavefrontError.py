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

import numpy as np


class SensorWavefrontError(object):
    """Contains the wavefront errors for a single sensor."""

    def __init__(self, numOfZk=19):
        """Constructs a sensor wavefront error.

        Parameters
        ----------
        numOfZk : int, optional
            Number of annular Zernike polynomials. (the default is 19.)
        """

        # Sensor Id
        self._sensorId = 999

        # Sensor Name
        self.sensorName = "R99_S99"

        # Number of zk
        self.numOfZk = int(numOfZk)

        # Annular Zernike polynomials (zk)
        self._annularZernikePoly = np.zeros(self.numOfZk)

    @property
    def sensorId(self):
        return self._sensorId

    @sensorId.setter
    def sensorId(self, newSensorId):
        """Set the sensor Id.

        Parameters
        ----------
        newSensorId : int
            The Id of the sensor this wavefront error is for.

        Raises
        ------
        ValueError
            sensorId must be >= 0.
        """

        if newSensorId < 0:
            raise ValueError("sensorId must be >= 0.")
        self._sensorId = int(newSensorId)

    @property
    def annularZernikePoly(self):
        return self._annularZernikePoly

    @annularZernikePoly.setter
    def annularZernikePoly(self, newAnnularZernikePoly):
        """Set the effective annular zernike poly.

        Parameters
        ----------
        newAnnularZernikePoly : numpy.ndarray[self.numOfZk] (float)
            The poly describing the wavefront error in um.

        Raises
        ------
        ValueError
            annularZernikePoly must be an array of self.numOfZk floats.
        """

        if len(newAnnularZernikePoly) != self.numOfZk:
            raise ValueError(
                "annularZernikePoly must be an array of %d floats." % self.numOfZk
            )
        self._annularZernikePoly = np.array(newAnnularZernikePoly)


if __name__ == "__main__":
    pass
