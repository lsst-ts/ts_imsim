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

from lsst.utils import getPackageDir
from lsst.obs.lsst import LsstComCam, LsstCam


def getModulePath():
    """Get the path of module.

    Returns
    -------
    str
        Directory path of module.
    """

    return getPackageDir("ts_imsim")


def getPolicyPath():
    """Get the path of the policy directory.

    Returns
    -------
    str
        Directory of policy files.
    """

    return os.path.join(getModulePath(), "policy")


def getConfigDir():
    """Get the directory of configuration files.

    Returns
    -------
    str
        Directory of configuration files.
    """

    return os.path.join(getPolicyPath(), "config")


def getCamera(instName):
    """Returns a lsst instrument for a given instrument name.

    Parameters
    ----------
    instName : `str`
        Instrument name. Valid options are 'comcam or 'lsstfam'.

    Returns
    -------
    camera : `lsst.afw.cameraGeom.Camera`

    Raises
    ------
    ValueError
        If input `instName` is not valid.
    """
    # Check the input
    if (instName == "lsstfam") or (instName == "lsst"):
        return LsstCam().getCamera()
    elif instName == "comcam":
        return LsstComCam().getCamera()
    else:
        raise ValueError(
            f"This instrument name ({instName}) is not supported. Must be 'comcam' or 'lsstfam' or 'lsst'."
        )
