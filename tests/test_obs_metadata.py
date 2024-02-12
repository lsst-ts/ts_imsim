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

from astroplan import FixedTarget, Observer
from astropy.time import Time
from lsst.ts.imsim import ObsMetadata


class TestObsMetadata(unittest.TestCase):
    """Test the ObsMetadata dataclass."""

    def setUp(self) -> None:
        self.obs_meta_test = ObsMetadata(ra=0.0, dec=0.0, band="r")

    def test_obs_metadata(self):
        self.assertEqual(self.obs_meta_test.ra, 0.0)
        self.assertEqual(self.obs_meta_test.dec, 0.0)
        self.assertEqual(self.obs_meta_test.band, "r")
        self.assertEqual(self.obs_meta_test.rotator_angle, 0.0)
        self.assertEqual(self.obs_meta_test.exp_time, 30.0)
        self.assertEqual(self.obs_meta_test.raw_seeing, 0.5)
        self.assertEqual(self.obs_meta_test.mjd, 60115.33)
        self.assertEqual(self.obs_meta_test.seq_num, 1)
        self.assertEqual(
            self.obs_meta_test.obs_id,
            "$f\"IM_P_{astropy.time.Time(mjd, format='mjd').strftime('%Y%m%d')}_{seqnum:06d}\" ",
        )
        self.assertEqual(self.obs_meta_test.focus_z, 0.0)
        self.assertAlmostEqual(self.obs_meta_test.zenith, 51.63121964610822)
        self.assertAlmostEqual(self.obs_meta_test.parallactic_angle, -130.1781932025212)

    def test_calc_parallactic_angle(self):
        sirius = FixedTarget.from_name("sirius")
        t = Time(60000, format="mjd")
        obs_metadata = ObsMetadata(
            ra=sirius.ra.deg, dec=sirius.dec.deg, band="r", mjd=60000
        )
        rubin = Observer.at_site("cerro pachon")
        self.assertAlmostEqual(
            obs_metadata.calc_parallactic_angle(),
            rubin.parallactic_angle(t, sirius).deg,
        )

    def test_calc_zenith_angle(self):
        sirius = FixedTarget.from_name("sirius")
        t = Time(59580.0, format="mjd")
        obs_metadata = ObsMetadata(
            ra=sirius.ra.deg, dec=sirius.dec.deg, band="r", mjd=59580.0
        )
        rubin = Observer.at_site("cerro pachon")
        self.assertAlmostEqual(
            obs_metadata.calc_zenith_angle(),
            90.0 - rubin.altaz(t, sirius).alt.deg,
        )

    def test_calc_alt_az(self):
        sirius = FixedTarget.from_name("sirius")
        t = Time(59580.0, format="mjd")
        obs_metadata = ObsMetadata(
            ra=sirius.ra.deg, dec=sirius.dec.deg, band="r", mjd=59580.0
        )
        rubin = Observer.at_site("cerro pachon")
        obs_alt, obs_az = obs_metadata.calc_alt_az()
        self.assertAlmostEqual(obs_alt, rubin.altaz(t, sirius).alt.deg)
        self.assertAlmostEqual(obs_az, rubin.altaz(t, sirius).az.deg)
