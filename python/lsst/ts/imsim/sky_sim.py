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

__all__ = ["SkySim"]

import numpy as np
from lsst.ts.imsim.utils import CamType, get_camera


class SkySim:
    """Initialization of sky simulator class."""

    def __init__(self) -> None:
        # Star ID
        self.star_id = np.array([], dtype=int)

        # Star RA
        self.ra = np.array([])

        # Star dec
        self.dec = np.array([])

        # Star magnitude
        self.mag = np.array([])

        # Camera
        self._camera = None

    def set_camera(self, cam_type: CamType) -> None:
        """Set the camera.

        Parameters
        ----------
        cam_type : lsst.ts.imsim.utils.CamType
            Camera type.
        """

        self._camera = get_camera(cam_type)

    def add_star_by_ra_dec_in_deg(
        self,
        star_id: int | list[int] | np.ndarray,
        ra_in_deg: float | list[float] | np.ndarray,
        dec_in_deg: float | list[float] | np.ndarray,
        mag: float | list[float] | np.ndarray,
    ) -> None:
        """Add the star information by (ra, dec) in degrees.

        Parameters
        ----------
        star_id : int, list[int], or numpy.ndarray[int]
            Star Id.
        ra_in_deg : float, list, or numpy.ndarray
            Star ra in degree.
        dec_in_deg : float, list, or numpy.ndarray
            Star dec in degree.
        mag : float, list, or numpy.ndarray
            Star magnitude.
        """

        # Check the inputs are list or not, and change the type if necessary
        star_id_list = self._change_to_list_if_necessary(star_id)
        ra_in_deg_list = self._change_to_list_if_necessary(ra_in_deg)
        dec_in_deg_list = self._change_to_list_if_necessary(dec_in_deg)
        mag_list = self._change_to_list_if_necessary(mag)

        # Add the stars
        for ii in range(len(star_id_list)):
            int_star_id = int(star_id_list[ii])
            if self._is_uniq_star_id(int_star_id):
                self.star_id = np.append(self.star_id, int_star_id)
                self.ra = np.append(self.ra, ra_in_deg_list[ii])
                self.dec = np.append(self.dec, dec_in_deg_list[ii])
                self.mag = np.append(self.mag, mag_list[ii])

    def _change_to_list_if_necessary(
        self, variable: int | float | list | np.ndarray
    ) -> list[int | float]:
        """Change the data type to list.

        Parameters
        ----------
        variable : int, float, list, or numpy.ndarray
            Variable.

        Returns
        -------
        list
            Variable as the list.
        """

        if isinstance(variable, (int, float)):
            return [variable]
        else:
            return variable

    def _is_uniq_star_id(self, star_id: int) -> bool:
        """Check the star ID is unique or not.

        Parameters
        ----------
        star_id : int
            Star Id.

        Returns
        -------
        bool
            True if the unique Id.
        """

        if star_id in self.star_id:
            is_unique = False
            print("StarId=%d is not unique." % star_id)
        else:
            is_unique = True

        return is_unique

    def add_star_by_file(self, read_file_path: str, skip_rows: int = 0) -> None:
        """Add the star data by reading the file.

        Parameters
        ----------
        read_file_path : str
            Star data file path.
        skip_rows : int, optional
            Skip the first "skiprows" lines. (the default is 0.)
        """

        data = np.loadtxt(read_file_path, skiprows=skip_rows)

        # Only consider the non-empty data
        if len(data) != 0:
            # Change to 2D array if the input is 1D array
            if data.ndim == 1:
                data = np.expand_dims(data, axis=0)

            for star in data:
                self.add_star_by_ra_dec_in_deg(star[0], star[1], star[2], star[3])

    def export_sky_to_file(self, output_file_path: str) -> None:
        """Export the star information into the file.

        Parameters
        ----------
        output_file_path : str
            Output file path.
        """

        # Add the header (star ID, ra, dec, magnitude)
        content = "# Id\t Ra\t\t Dec\t\t Mag\n"

        # Add the star information
        for ii in range(len(self.star_id)):
            content += "%d\t %3.6f\t %3.6f\t %3.6f\n" % (
                self.star_id[ii],
                self.ra[ii],
                self.dec[ii],
                self.mag[ii],
            )

        # Write into file
        file_out = open(output_file_path, "w")
        file_out.write(content)
        file_out.close()
