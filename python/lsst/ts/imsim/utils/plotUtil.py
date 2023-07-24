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
import matplotlib.pyplot as plt


def _saveFig(plt, saveToFilePath=None, dpi=None):
    """Save the figure.

    Parameters
    ----------
    plt : module
        Module of matplotlib.pyplot.
    saveToFilePath : str, optional
        File path to save the figure. If None, the figure will be showed. (the
        default is None.)
    dpi : int, optional
        The resolution in dots per inch. (the default is None.)
    """

    if saveToFilePath is not None:
        plt.savefig(saveToFilePath, dpi=dpi)
        plt.close()
    else:
        plt.show()


def plotFwhmOfIters(pssnFiles, saveToFilePath=None, dpi=None):
    """Plot the FWHM of iteration.

    FWHM: Full width at half maximum.
    PSSN: Normalized point source sensitivity.

    Parameters
    ----------
    pssnFiles : list
        List of PSSN files.
    saveToFilePath : str, optional
        File path to save the figure. If None, the figure will be showed. (the
        default is None.)
    dpi : int, optional
        The resolution in dots per inch. (the default is None.)
    """

    # Collect the FWHM data. The row is the FWHM for each field. The final row
    # is the GQ FWHM. The column is the iterations.
    numOfIter = len(pssnFiles)

    numOfFwhmData = 0
    fwhmDataAll = np.array([])
    for pssnFile in pssnFiles:
        fwhmData = np.loadtxt(pssnFile)[1, :]
        if numOfFwhmData == 0:
            numOfFwhmData = len(fwhmData)

        fwhmDataAll = np.append(fwhmDataAll, fwhmData)

    reshapedFwhmData = fwhmDataAll.reshape((numOfIter, numOfFwhmData)).T

    # Plot the figure
    plt.figure()
    plt.plot(reshapedFwhmData[:-1, :].T, "bx-")
    plt.plot(reshapedFwhmData[-1, :], "ro-", label="GQ FWHM_eff")
    plt.xlabel("Iteration")
    plt.ylabel("Arcsec")
    plt.legend()

    _saveFig(plt, saveToFilePath=saveToFilePath, dpi=dpi)


if __name__ == "__main__":
    pass
