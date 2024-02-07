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

__all__ = ["OpdMetrology"]

import numpy as np
import yaml
from astropy.io import fits
from lsst.afw.cameraGeom import FIELD_ANGLE
from lsst.ts.imsim.utils import (
    CamType,
    calc_pssn,
    get_camera,
    zernike_annular_fit,
    zernike_eval,
)
from lsst.ts.ofc.utils import get_config_dir as getConfigDirOfc


class OpdMetrology:
    """Initialization of OPD metrology class.

    OPD: Optical path difference.
    """

    def __init__(self) -> None:
        self._wt = np.array([])
        self.field_x = np.array([])
        self.field_y = np.array([])
        self.sensor_ids = []

    @property
    def wt(self) -> list | np.ndarray:
        return self._wt

    @wt.setter
    def wt(self, new_wt: list | np.ndarray) -> None:
        """Set the weighting ratio used in Gaussian quadrature.

        Parameters
        ----------
        new_wt : list or numpy.ndarray
            Weighting ratio.

        Raises
        ------
        ValueError
            All weighting ratios should be >=0.
        """

        wt_arr = np.array(new_wt, dtype=float)
        if np.all(wt_arr >= 0):
            self._wt = wt_arr / np.sum(wt_arr)
        else:
            raise ValueError("All weighting ratios should be >= 0.")

    def set_default_lsst_wfs_gq(self) -> None:
        """Set default values for LSST WFS field X, Y
        and weighting ratio.
        """

        # Set equal full weights for each of the
        # four corner wavefront sensor pairs.
        self.wt = np.array([1.0, 1.0, 1.0, 1.0])
        wfs_field_x, wfs_field_y, sensor_ids = self.get_default_lsst_wfs_gq()
        self.field_x, self.field_y = (wfs_field_x, wfs_field_y)
        self.sensor_ids = sensor_ids

    def get_default_lsst_wfs_gq(self) -> tuple[list[float], list[float], list[int]]:
        """Get the default field X, Y of LSST WFS on GQ.

        WFS: Wavefront sensor.
        GQ: Gaussian quadrature

        Returns
        -------
        list
            Field X in degree.
        list
            Field Y in degree.
        list
            Detector IDs of extra-focal sensor in raft.
        """

        # Field x, y for 4 WFS in the Camera Coordinate System (CCS)
        # These will be chosen at the center of the extra-intra
        # focal pairs of wavefront sensors.
        det_ids = [191, 195, 199, 203]
        camera = get_camera(CamType.LsstCam)
        field_x = []
        field_y = []
        det_map = camera.getIdMap()
        for det_id in det_ids:
            det_extra_center = np.degrees(det_map[det_id].getCenter(FIELD_ANGLE))
            det_intra_center = np.degrees(det_map[det_id + 1].getCenter(FIELD_ANGLE))
            # Switch X,Y coordinates to convert from DVCS to CCS coords
            det_center = np.mean([det_extra_center, det_intra_center], axis=0)
            field_y.append(det_center[0])
            field_x.append(det_center[1])

        return field_x, field_y, det_ids

    def set_wgt_and_field_xy_of_gq(self, cam_type: CamType) -> None:
        """Set the GQ weighting ratio and field X, Y.

        GQ: Gaussian quadrature.

        Parameters
        ----------
        cam_type : lsst.ts.imsim.utils.CamType
            Camera type.
            Valid CamTypes are LsstCam, LsstFamCam, ComCam.
        """

        # Set camera and field ids for given instrument
        if cam_type == CamType.LsstCam:
            self.set_default_lsst_wfs_gq()
            return

        weight_dir_path = getConfigDirOfc() / "image_quality_weights"
        path_wgt_file = weight_dir_path / f"{cam_type.value}_weights.yaml"

        # Set the weighting ratio
        with open(path_wgt_file, "r") as file:
            wgt = yaml.safe_load(file)
        wgt_values = np.array(list(wgt.values()), dtype=float)
        # Normalize weights
        self.wt = wgt_values / np.sum(wgt_values)

        camera = get_camera(cam_type)
        if cam_type == CamType.LsstFamCam:
            self.sensor_ids = np.arange(189)
        elif cam_type == CamType.ComCam:
            self.sensor_ids = np.arange(9)

        field_x = []
        field_y = []
        det_map = camera.getIdMap()
        for det_id in self.sensor_ids:
            det_center = det_map[det_id].getCenter(FIELD_ANGLE)
            # Switch X,Y coordinates to convert from DVCS to CCS coords
            field_y.append(np.degrees(det_center[0]))
            field_x.append(np.degrees(det_center[1]))
        self.field_x = np.array(field_x)
        self.field_y = np.array(field_y)

    def get_zk_from_opd(
        self,
        opd_fits_file: str | None = None,
        opd_map: np.ndarray | None = None,
        zk_terms: int = 22,
        obscuration: float = 0.61,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get the wavefront error of OPD in the basis of annular Zernike
        polynomials.

        OPD: Optical path difference.

        Parameters
        ----------
        opd_fits_file : str, optional
            OPD FITS file. (the default is None.)
        opd_map : numpy.ndarray, optional
            OPD map data. (the default is None.)
        zk_terms : int, optional
            Number of terms of annular Zk (z1-z22 by default). (the default
            is 22.)
        obscuration : float, optional
            Obscuration of annular Zernike polynomial. (the default is 0.61.)

        Returns
        -------
        numpy.ndarray
            Annular Zernike polynomials. For imSim OPD, the unit is nm.
        numpy.ndarray
            OPD map.
        numpy.ndarray
            Meshgrid x in OPD map.
        numpy.ndarray
            Meshgrid y in OPD map.

        Raises
        ------
        ValueError
            The x, y dimensions of OPD are different.
        """

        # Get the OPD data (imSim OPD unit: nm)
        if opd_fits_file is not None:
            opd = fits.getdata(opd_fits_file)
        elif opd_map is not None:
            opd = opd_map.copy()

        # Check the x, y dimensions of OPD are the same
        if np.unique(opd.shape).size != 1:
            raise ValueError("The x, y dimensions of OPD are different.")

        # x-, y-coordinate in the OPD image
        opd_size = opd.shape[0]
        opd_grid_1d = np.linspace(-1, 1, opd_size)
        opd_x, opd_y = np.meshgrid(opd_grid_1d, opd_grid_1d)

        # Fit the OPD map with Zk and write into the file
        idx = ~np.isnan(opd)
        zk = zernike_annular_fit(
            opd[idx], opd_x[idx], opd_y[idx], zk_terms, obscuration
        )

        return zk, opd, opd_x, opd_y

    def rm_piston_tip_tilt_from_opd(
        self, opd_fits_file: str | None = None, opd_map: np.ndarray | None = None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Remove the piston (z1), x-tilt (z2), and y-tilt (z3)
        from the OPD map.

        OPD: Optical path difference.

        Parameters
        ----------
        opd_fits_file : str, optional
            OPD FITS file. (the default is None.)
        opd_map : numpy.ndarray, optional
            OPD map data. (the default is None.)

        Returns
        -------
        numpy.ndarray
            OPD map after removing the affection of z1-z3.
        numpy.ndarray
            Meshgrid x in OPD map.
        numpy.ndarray
            Meshgrid y in OPD map.
        """

        # Do the spherical Zernike fitting for the OPD map
        # Only fit the first three terms (z1-z3): piston, x-tilt, y-tilt
        zk, opd, opd_x, opd_y = self.get_zk_from_opd(
            opd_fits_file=opd_fits_file, opd_map=opd_map, zk_terms=3, obscuration=0
        )

        # Find the index that the value of OPD is not 0
        idx = ~np.isnan(opd)

        # Remove the PTT
        opd[idx] -= zernike_eval(zk, opd_x[idx], opd_y[idx])

        return opd, opd_x, opd_y

    def calc_pssn(
        self,
        wavelength_in_um: float,
        opd_fits_file: str | None = None,
        opd_map: np.ndarray | None = None,
        zen: int = 0,
        debug_level: int = 0,
    ) -> float:
        """Calculate the PSSN based on OPD map.

        PSSN: Normalized point source sensitivity.
        OPD: Optical path difference.

        Parameters
        ----------
        wavelength_in_um : float
            Wavelength in microns.
        opd_fits_file : str, optional
            OPD FITS file. Units need to be microns. (the default is None.)
        opd_map : numpy.ndarray, optional
            OPD map data. Units need to be microns. (the default is None.)
        zen : float, optional
            Telescope zenith angle in degree. (the default is 0.)
        debug_level : int, optional
            Debug level. The higher value gives more information. (the default
            is 0.)

        Returns
        -------
        float
            Calculated PSSN.
        """

        # Before calc_pssn,
        # (1) Remove PTT (piston, x-tilt, y-tilt),
        # (2) Make sure outside of pupil are all zeros
        opd_rm_ptt = self.rm_piston_tip_tilt_from_opd(
            opd_fits_file=opd_fits_file, opd_map=opd_map
        )[0]

        # Calculate the normalized point source sensitivity (PSSN)
        pssn = calc_pssn(opd_rm_ptt, wavelength_in_um, zen=zen, debug_level=debug_level)

        return pssn

    def calc_gq_value(self, value_list: list[float] | np.ndarray) -> float:
        """Calculate the GQ value.

        GQ: Gaussian quadrature

        Parameters
        ----------
        value_list : list or numpy.ndarray
            List of value (PSSN, effective FWHM, dm5, ellipticity).

        Returns
        -------
        float
            GQ value.

        Raises
        ------
        ValueError
            Length of wt ratio != length of value list.
        """
        # Check the lengths of weighting ratio and value list are the same
        if len(self.wt) != len(value_list):
            raise ValueError("Length of wt ratio != length of value list.")

        # Calculate the effective value on Gaussain quardure plane
        value_array = np.array(value_list, dtype=float)
        gq_value = np.sum(self.wt * value_array)

        return gq_value

    def calc_fwhm_eff(self, pssn: float) -> float:
        """Calculate the effective FWHM.

        FWHM: Full width at half maximum.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssn : float
            PSSN value.

        Returns
        -------
        float
            Effective FWHM.
        """

        # FWHMeff_sys = FWHMeff_atm * sqrt(1/PSSN - 1).
        # FWHMeff_atm = 0.6 arcsec.
        # Another correction factor (eta = 1.086) is used to account for the
        # difference between the
        # simple RSS and the more proper convolution.
        # Follow page 7 (section 7.2 FWHM) in document-17242 for more
        # information.
        eta = 1.086
        fwhm_atm = 0.6
        fwhm_eff = eta * fwhm_atm * np.sqrt(1 / pssn - 1)

        return fwhm_eff
