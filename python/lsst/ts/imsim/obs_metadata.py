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

__all__ = ["ObsMetadata"]

from dataclasses import dataclass, field

import astropy
from astroplan import Observer


@dataclass
class ObsMetadata:
    """Class for storing observation metadata."""

    ra: float  # ra in degrees
    dec: float  # dec in degrees
    band: str
    rotator_angle: float = 0.0
    exp_time: float = 30.0
    mjd: float = 59580.0
    seq_num: int = 1
    raw_seeing: int = 0.5
    obs_id: str = """$f"IM_P_{astropy.time.Time(mjd, format='mjd').strftime('%Y%m%d')}_{seqnum:06d}" """
    focus_z: float = 0.0  # Defocal distance in mm
    zenith: float = field(init=False)
    parallactic_angle: float = field(init=False)

    def __post_init__(self) -> None:
        """Populate the zenith and parallactic angles
        with values from init.
        """
        self.zenith = self.calc_zenith_angle()
        self.parallactic_angle = self.calc_parallactic_angle()

    def format_observation_info(self) -> (Observer, astropy.time, astropy.coordinates):
        """
        Get the observation info in the data structures needed
        to calculate observer information such as parallactic angle
        and zenith angle.

        Parameters
        ----------
        obs_metadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.

        Returns
        -------
        Observer
            astroplan.Observer for observer located at Rubin Observatory
        astropy.time
            Observation time
        astropy.coordinates
            Observation boresight
        """

        # Observer located at Rubin
        rubin_observer = Observer.at_site("cerro pachon")
        time = astropy.time.Time(self.mjd, format="mjd")
        boresight = astropy.coordinates.SkyCoord(
            f"{self.ra}d", f"{self.dec}d", frame="icrs"
        )

        return rubin_observer, time, boresight

    def calc_parallactic_angle(self) -> float:
        """Calculate the parallactic angle so we know the
        sky rotation angle on alt-az mount for the observation.

        Parameters
        ----------
        obs_metadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.

        Returns
        -------
        float
            Parallactic Angle in degrees.
        """

        observer, time, boresight = self.format_observation_info()

        return observer.parallactic_angle(time, boresight).deg

    def calc_zenith_angle(self) -> float:
        """Calculate the zenith angle so we can accurately
        populate the imsim configuration.

        Parameters
        ----------
        obs_metadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.

        Returns
        -------
        float
            Zenith Angle in degrees.
        """

        observer, time, boresight = self.format_observation_info()

        # Altitude refers to elevation angle up from horizon.
        # To get zenith we need to convert to the angle down
        # from zenith so we subtract altitude from 90 degrees.
        return 90.0 - observer.altaz(time, boresight).alt.deg
