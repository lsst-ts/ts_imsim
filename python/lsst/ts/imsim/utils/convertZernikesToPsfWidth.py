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

import galsim
import numpy as np


def getPsfGradPerZernike(
    *,
    jmin: int = 4,
    jmax: int = 22,
    R_outer: float = 4.18,
    R_inner: float = 2.5498,
) -> np.ndarray:
    """Get the gradient of the PSF FWHM with respect to each Zernike.

    This function takes no positional arguments. All parameters must be passed
    by name (see the list of parameters below).

    Parameters
    ----------
    jmin : int
        The minimum Zernike Noll index, inclusive.
        (the default, 4, ignores piston, x & y offsets, and tilt.)
    jmax : int
        The max Zernike Noll index, inclusive. (the default is 22.)
    R_outer : float
        The outer radius of the telescope aperture, in meters.
        (the default, 4.18, corresponds to the LSST primary mirror.)
    R_inner : float
        The inner radius of the telescope aperture, in meters.
        (the default, 2.5498, corresponds to the LSST primary mirror.)

    Returns
    -------
    np.ndarray
        Gradient of the PSF FWHM with respect to the corresponding Zernike.
        Units are arcsec / micron.
    """
    # Calculate the conversion factors
    conversion_factors = np.zeros(jmax + 1)
    for i in range(jmin, jmax + 1):
        # Set coefficients for this Noll index: coefs = [0, 0, ..., 1]
        # Note the first coefficient is Noll index 0, which does not exist and
        # is therefore always ignored by galsim
        coefs = [0] * i + [1]

        # Create the Zernike polynomial with these coefficients
        Z = galsim.zernike.Zernike(coefs, R_outer=R_outer, R_inner=R_inner)

        # We can calculate the size of the PSF from the RMS of the gradient of
        # the wavefront. The gradient of the wavefront perturbs photon paths.
        # The RMS quantifies the size of the collective perturbation.
        # If we expand the wavefront gradient in another series of Zernike
        # polynomials, we can exploit the orthonormality of the Zernikes to
        # calculate the RMS from the Zernike coefficients.
        rms_tilt = np.sqrt(np.sum(Z.gradX.coef**2 + Z.gradY.coef**2) / 2)

        # Convert to arcsec per micron
        rms_tilt = np.rad2deg(rms_tilt * 1e-6) * 3600

        # Convert rms -> fwhm
        fwhm_tilt = 2 * np.sqrt(2 * np.log(2)) * rms_tilt

        # Save this conversion factor
        conversion_factors[i] = fwhm_tilt

    return conversion_factors[jmin:]


def convertZernikesToPsfWidth(
    zernikes: np.ndarray,
    jmin: int = 4,
    R_outer: float = 4.18,
    R_inner: float = 2.5498,
) -> np.ndarray:
    """Convert Zernike amplitudes to quadrature contribution to the PSF FWHM.

    Parameters
    ----------
    zernikes : np.ndarray
        Zernike amplitudes (in microns), starting with Noll index `jmin`.
        Either a 1D array of zernike amplitudes, or a 2D array, where each row
        corresponds to a different set of amplitudes.
    jmin : int
        The minimum Zernike Noll index, inclusive. The maximum Noll index is
        inferred from `jmin` and the length of `zernikes`.
        (the default is 4, which ignores piston, x & y offsets, and tilt.)
    R_outer : float
        The outer radius of the telescope aperture, in meters.
        (the default, 4.18, corresponds to the LSST primary mirror.)
    R_inner : float
        The inner radius of the telescope aperture, in meters.
        (the default, 2.5498, corresponds to the LSST primary mirror.)

    Returns
    -------
    dFWHM: np.ndarray
        Quadrature contribution of each Zernike vector to the PSF FWHM
        (in arcseconds).

    Notes
    -----
    Converting Zernike amplitudes to their quadrature contributions to the PSF
    FWHM allows for easier physical interpretation of Zernike amplitudes and
    the performance of the AOS system.

    For example, image we have a true set of zernikes, [Z4, Z5, Z6], such that
    ConvertZernikesToPsfWidth([Z4, Z5, Z6]) = [0.1, -0.2, 0.3] arcsecs.
    These Zernike perturbations increase the PSF FWHM by
    sqrt[(0.1)^2 + (-0.2)^2 + (0.3)^2] ~ 0.37 arcsecs.

    If the AOS perfectly corrects for these perturbations, the PSF FWHM will
    not increase in size. However, imagine the AOS estimates zernikes, such
    that ConvertZernikesToPsfWidth([Z4, Z5, Z6]) = [0.1, -0.3, 0.4] arcsecs.
    These estimated Zernikes, do not exactly match the true Zernikes above.
    Therefore, the post-correction PSF will still be degraded with respect to
    the optimal PSF. In particular, the PSF FWHM will be increased by
    sqrt[(0.1 - 0.1)^2 + (-0.2 - (-0.3))^2 + (0.3 - 0.4)^2] ~ 0.14 arcsecs.

    This conversion depends on a linear approximation that begins to break down
    for RSS(dFWHM) > 0.20 arcsecs. Beyond this point, the approximation tends
    to overestimate the PSF degradation. In other words, if
    sqrt(sum( dFWHM^2 )) > 0.20 arcsec, it is likely that dFWHM is
    over-estimated. However, the point beyond which this breakdown begins
    (and whether the approximation over- or under-estimates dFWHM) can change,
    depending on which Zernikes have large amplitudes. In general, if you have
    large Zernike amplitudes, proceed with caution!
    Note that if the amplitudes Z_est and Z_true are large, this is okay, as
    long as |Z_est - Z_true| is small.

    For a notebook demonstrating where the approximation breaks down:
    https://gist.github.com/jfcrenshaw/24056516cfa3ce0237e39507674a43e1
    """
    # Calculate jmax from jmin and the length of the zernike array
    jmax = jmin + zernikes.shape[-1] - 1

    # Calculate the conversion factors for each zernike
    conversion_factors = getPsfGradPerZernike(
        jmin=jmin,
        jmax=jmax,
        R_outer=R_outer,
        R_inner=R_inner,
    )

    # Convert the Zernike amplitudes from microns to their quadrature
    # contribution to the PSF FWHM
    dFWHM = conversion_factors * zernikes

    return dFWHM