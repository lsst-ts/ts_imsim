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

from lsst.ts.imsim.obsMetadata import ObsMetadata


class TestObsMetadata(unittest.TestCase):
    """Test the ObsMetadata dataclass."""

    def testObsMetadata(self):
        obsMetaTest = ObsMetadata(ra=0.0, dec=0.0, band="r")
        self.assertEqual(obsMetaTest.ra, 0.0)
        self.assertEqual(obsMetaTest.dec, 0.0)
        self.assertEqual(obsMetaTest.band, "r")
        self.assertEqual(obsMetaTest.zenith, 0.0)
        self.assertEqual(obsMetaTest.rotatorAngle, 0.0)
        self.assertEqual(obsMetaTest.expTime, 30.0)
        self.assertEqual(obsMetaTest.mjd, 59580.0)
        self.assertEqual(obsMetaTest.seqNum, 1)
        self.assertEqual(
            obsMetaTest.obsId,
            "$f\"IM_P_{astropy.time.Time(mjd, format='mjd').strftime('%Y%m%d')}_{seqnum:06d}\" ",
        )
