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

from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np


def _save_fig(
    plt: plt, save_to_file_path: Optional[str] = None, dpi: Optional[int] = None
) -> None:
    """Save the figure.

    Parameters
    ----------
    plt : module
        Module of matplotlib.pyplot.
    save_to_file_path : str, optional
        File path to save the figure. If None, the figure will be showed. (the
        default is None.)
    dpi : int, optional
        The resolution in dots per inch. (the default is None.)
    """

    if save_to_file_path is not None:
        plt.savefig(save_to_file_path, dpi=dpi)
        plt.close()
    else:
        plt.show()


def plot_fwhm_of_iters(
    pssn_files: List[str], save_to_file_path: Optional[str] = None, dpi: None = None
) -> None:
    """Plot the FWHM of iteration.

    FWHM: Full width at half maximum.
    PSSN: Normalized point source sensitivity.

    Parameters
    ----------
    pssn_files : list
        List of PSSN files.
    save_to_file_path : str, optional
        File path to save the figure. If None, the figure will be showed. (the
        default is None.)
    dpi : int, optional
        The resolution in dots per inch. (the default is None.)
    """

    # Collect the FWHM data. The row is the FWHM for each field. The final row
    # is the GQ FWHM. The column is the iterations.
    num_of_iter = len(pssn_files)

    num_of_fwhm_data = 0
    fwhm_data_all = np.array([])
    for pssn_file in pssn_files:
        fwhm_data = np.loadtxt(pssn_file)[1, :]
        if num_of_fwhm_data == 0:
            num_of_fwhm_data = len(fwhm_data)

        fwhm_data_all = np.append(fwhm_data_all, fwhm_data)

    reshaped_fwhm_data = fwhm_data_all.reshape((num_of_iter, num_of_fwhm_data)).T

    # Plot the figure
    plt.figure()
    plt.plot(reshaped_fwhm_data[:-1, :].T, "bx-")
    plt.plot(reshaped_fwhm_data[-1, :], "ro-", label="GQ FWHM_eff")
    plt.xlabel("Iteration")
    plt.ylabel("Arcsec")
    plt.legend()

    _save_fig(plt, save_to_file_path=save_to_file_path, dpi=dpi)


if __name__ == "__main__":
    pass
