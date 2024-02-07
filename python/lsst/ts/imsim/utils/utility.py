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

__all__ = [
    "get_module_path",
    "get_policy_path",
    "get_config_dir",
    "get_camera",
    "make_dir",
    "get_zk_from_file",
    "ModifiedEnvironment",
    "CamType",
    "get_cam_type",
    "zernike_annular_fit",
    "zernike_eval",
]

import os
from enum import StrEnum

import numpy as np
import yaml
from galsim.zernike import Zernike as GSZernike
from lsst.afw import cameraGeom
from lsst.obs.lsst import LsstCam, LsstComCam
from lsst.utils import getPackageDir
from numpy import ndarray


class CamType(StrEnum):
    """
    Define allowed camera types.
    """

    LsstCam = "lsst"
    LsstFamCam = "lsstfam"
    ComCam = "comcam"


def get_module_path() -> str:
    """Get the path of module.

    Returns
    -------
    str
        Directory path of module.
    """

    return getPackageDir("ts_imsim")


def get_policy_path() -> str:
    """Get the path of the policy directory.

    Returns
    -------
    str
        Directory of policy files.
    """

    return os.path.join(get_module_path(), "policy")


def get_config_dir() -> str:
    """Get the directory of configuration files.

    Returns
    -------
    str
        Directory of configuration files.
    """

    return os.path.join(get_policy_path(), "config")


def get_camera(cam_type: CamType) -> cameraGeom.Camera:
    """Returns a lsst instrument for a given instrument name.

    Parameters
    ----------
    cam_type : lsst.ts.imsim.utils.CamType
        Camera type.

    Returns
    -------
    camera : `lsst.afw.cameraGeom.Camera`

    Raises
    ------
    ValueError
        If input `cam_type` is not valid.
    """
    # Check the input
    if cam_type in [CamType.LsstCam, CamType.LsstFamCam]:
        return LsstCam().getCamera()
    elif cam_type == CamType.ComCam:
        return LsstComCam().getCamera()
    else:
        raise ValueError(
            f"This CamType ({cam_type}) is not supported. Must be LsstCam, LsstFamCam or ComCam."
        )


def make_dir(new_dir: str, exist_ok: bool = True) -> None:
    """Make the new directory.

    Super-mkdir; create a leaf directory and all intermediate ones. Works
    like mkdir, except that any intermediate path segment (not just the
    rightmost) will be created if it does not exist.

    Parameters
    ----------
    new_dir : str
        New directory.
    exist_ok : bool, optional
        If the target directory already exists, raise an OSError if
        exist_ok is False. Otherwise no exception is raised. (the default
        is True.)
    """

    os.makedirs(new_dir, exist_ok=exist_ok)


def get_zk_from_file(zk_file_path: str) -> dict[int, ndarray]:
    """Get the zk (z4-z22) from file.

    Parameters
    ----------
    zk_file_path : str
        Zk file path.

    Returns
    -------
    numpy.ndarray
        zk matrix. The colunm is z4-z22. The raw is each data point.
    """

    with open(zk_file_path, "r") as file:
        zk = yaml.safe_load(file)
    for key, val in zk.items():
        zk[key] = np.fromstring(val[0], sep=" ")

    return zk


class ModifiedEnvironment:
    """Context manager to temporarily modify shell environment with specified
    overrides.

    Parameters
    ----------
    **kwargs : dict
        Dictionary of environment variables to override. Keys and values should
        be type str, unless the value is None, in which case the environment
        variable will be unset.
    """

    def __init__(self, **kwargs) -> None:
        self._overrides = kwargs
        self._originals = {}

    def __enter__(self) -> None:
        for key, value in self._overrides.items():
            # Save original value if exists
            if key in os.environ:
                self._originals[key] = os.environ[key]

            # Set new values or unset if new value is None
            if value is None:
                del os.environ[key]
            else:
                os.environ[key] = value

    def __exit__(self, exc_type: None, exc_val: None, exc_tb: None) -> None:
        for key in self._overrides:
            # Restore original values, delete or keep as they are
            # based on the original state.
            if key in self._originals:
                os.environ[key] = self._originals[key]
            else:
                del os.environ[key]


def get_cam_type(inst_name: str) -> CamType:
    """Get the camera type from instrument name.

    Parameters
    ----------
    inst_name : str
         Instrument name.

    Returns
    -------
    enum 'CamType'
        Camera type.

    Raises
    ------
    ValueError
        Instrument name is not supported.
    """
    if inst_name == "lsst":
        return CamType.LsstCam
    elif inst_name == "lsstfam":
        return CamType.LsstFamCam
    elif inst_name == "comcam":
        return CamType.ComCam
    else:
        raise ValueError(f"Instrument name ({inst_name}) is not supported.")


def zernike_annular_fit(
    s: np.ndarray, x: np.ndarray, y: np.ndarray, numTerms: int, e: float
) -> np.ndarray:
    """Get the coefficients of annular Zernike polynomials by fitting the
    wavefront surface.

    Parameters
    ----------
    s : numpy.ndarray
        Wavefront surface to be fitted.
    x : numpy.ndarray
        Normalized x coordinate between -1 and 1 (pupil coordinate).
    y : numpy.ndarray
        Normalized y coordinate between -1 and 1 (pupil coordinate).
    numTerms : int
        Number of annular Zernike terms used in the fit.
    e : float
        Obscuration ratio of annular Zernikes.

    Returns
    -------
    numpy.ndarray
        Coefficients of annular Zernike polynomials by the fitting.
    """
    # Check the dimensions of x and y are the same or not
    if x.shape != y.shape:
        print("x & y are not the same size")

    # Get the value that is finite
    sFinite = s[:].copy()
    xFinite = x[:].copy()
    yFinite = y[:].copy()

    finiteIndex = np.isfinite(sFinite + xFinite + yFinite)

    sFinite = sFinite[finiteIndex]
    xFinite = xFinite[finiteIndex]
    yFinite = yFinite[finiteIndex]

    # Do the fitting
    h = np.zeros([len(sFinite), numTerms])

    for ii in range(numTerms):
        z = np.zeros(numTerms)
        z[ii] = 1
        h[:, ii] = GSZernike(np.concatenate([[0], z]), R_inner=e)(xFinite, yFinite)

    # Solve the equation: H*Z = S => Z = H^(-1)S
    z = np.linalg.lstsq(h, sFinite, rcond=None)[0]

    return z


def zernike_eval(z: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Calculate the wavefront surface in the basis of Zernike polynomial.

    Parameters
    ----------
    z : numpy.ndarray
        Coefficient of Zernike polynomials.
    x : numpy.ndarray
        X coordinate on pupil plane.
    y : numpy.ndarray
        Y coordinate on pupil plane.

    Returns
    -------
    numpy.ndarray
        Wavefront surface.
    """
    # Calculate the wavefront surface
    # Use obscuration (e) = 0 for standard Zernike polynomials
    return GSZernike(np.concatenate([[0], z]), R_inner=0)(x, y)
