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
import unittest

from lsst.ts.imsim.utils.plotUtil import plot_fwhm_of_iters
from lsst.ts.imsim.utils.utility import get_module_path


class TestPlotUtil(unittest.TestCase):
    """Test the PlotUtil functions."""

    def setUp(self):
        module_path = get_module_path()
        self.test_data = os.path.join(module_path, "tests", "testData")
        self.out_fig_file_pat = os.path.join(
            module_path, "output", "img", "testFig.png"
        )

    def tearDown(self):
        if os.path.exists(self.out_fig_file_pat):
            os.remove(self.out_fig_file_pat)

    def test_plot_fwhm_of_iters(self):
        iter_data_dir = os.path.join(self.test_data, "iterData")
        pssn_files = [
            os.path.join(iter_data_dir, "iter%d" % num, "img", "PSSN.txt")
            for num in range(5)
        ]

        plot_fwhm_of_iters(pssn_files, save_to_file_path=self.out_fig_file_pat)
        self.assertTrue(os.path.exists(self.out_fig_file_pat))


if __name__ == "__main__":
    # Run the unit test
    unittest.main()
