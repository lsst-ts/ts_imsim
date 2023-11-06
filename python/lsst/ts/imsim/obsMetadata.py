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

from dataclasses import dataclass


@dataclass
class ObsMetadata:
    """Class for storing observation metadata."""

    ra: float  # ra in degrees
    dec: float  # dec in degrees
    band: str
    zenith: float = 0.0
    rotator_angle: float = 0.0
    exp_time: float = 30.0
    mjd: float = 59580.0
    seq_num: int = 1
    raw_seeing: int = 0.5
    obs_id: str = """$f"IM_P_{astropy.time.Time(mjd, format='mjd').strftime('%Y%m%d')}_{seqnum:06d}" """
    focus_z: float = 0.0  # Defocal distance in mm
