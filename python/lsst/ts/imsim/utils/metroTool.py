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

import warnings

import numpy as np
import scipy.special as sp
from lsst.ts.wep.cwfs.tool import extractArray, padArray


def calc_pssn(
    array,
    wlum,
    aType="opd",
    D=8.36,
    r0inmRef=0.1382,
    zen=0,
    pmask=0,
    imagedelta=0,
    fno=1.2335,
    debugLevel=0,
):
    """Calculate the normalized point source sensitivity (PSSN).

    Parameters
    ----------
    array : numpy.ndarry
        Array that contains either opd or pdf. opd need to be in microns.
    wlum : float
        Wavelength in microns.
    aType : str, optional
        What is used to calculate pssn - either opd or psf. (the default is
        "opd".)
    D : float, optional
        Side length of OPD image in meter. (the default is 8.36.)
    r0inmRef : float, optional
        Fidicial atmosphere r0 @ 500nm in meter, Konstantinos uses 0.20. (the
        default is 0.1382.)
    zen : float, optional
        Telescope zenith angle in degree. (the default is 0.)
    pmask : int or numpy.ndarray[int], optional
        Pupil mask. when opd is used, it can be generated using opd image, we
        can put 0 or -1 or whatever here. When psf is used, this needs to be
        provided separately with same size as array. (the default is 0.)
    imagedelta : float, optional
        Only needed when psf is used. use 0 for opd. (the default is 0.)
    fno : float, optional
        Only needed when psf is used. use 0 for opd. (the default is 1.2335.)
    debugLevel : int, optional
        Debug level. The higher value gives more information. (the default
        is 0.)

    Returns
    -------
    float
        PSSN value.
    """

    # Only needed for psf: pmask, imagedelta, fno

    # THE INTERNAL RESOLUTION THAT FFTS OPERATE ON IS VERY IMPORTANT
    # TO THE ACCUARCY OF PSSN.
    # WHEN TYPE='OPD', NRESO=SIZE(ARRAY,1)
    # WHEN TYPE='PSF', NRESO=SIZE(PMASK,1)
    #    for the psf option, we can not first convert psf back to opd then
    #    start over,
    #    because psf=|exp(-2*OPD)|^2. information has been lost in the | |^2.
    #    we need to go forward with psf->mtf,
    #    and take care of the coordinates properly.

    # PSSN = (n_eff)_atm / (n_eff)_atm+sys
    # (n_eff))_atm = 1 / (int (PSF^2)_atm dOmega)
    # (n_eff))_atm+sys = 1 / (int (PSF^2)_atm+sys dOmega)

    # Check the type is "OPD" or "PSF"
    if aType not in ("opd", "psf"):
        raise ValueError("The type of %s is not allowed." % aType)

    # Squeeze the array if necessary
    if array.ndim == 3:
        array2D = array[0, :, :].squeeze()

    # Get the k value (magnification ratio used in creating MTF)
    if aType == "opd":
        try:
            m = max(array2D.shape)
        except NameError:
            m = max(array.shape)
        k = 1
    elif aType == "psf":
        m = max(pmask.shape)
        # Pupil needs to be padded k times larger to get imagedelta
        # Do not know where to find this formular. Check with Bo.
        k = fno * wlum / imagedelta

    # Get the modulation transfer function with the van Karman power spectrum
    mtfa = createMTFatm(D, m, k, wlum, zen, r0inmRef, model="vonK")

    # Get the pupil function
    if aType == "opd":
        try:
            iad = array2D != 0
        except NameError:
            iad = array != 0
    elif aType == "psf":
        # Add even number
        mk = int(m + np.rint((m * (k - 1) + 1e-5) / 2) * 2)
        # padArray(pmask, m)
        iad = pmask

    # OPD --> PSF --> OTF --> OTF' (OTF + atmosphere) --> PSF'
    # Check with Bo that we could get OTF' or PSF' from PhoSim or not directly.
    # The above question might not be a concern in the simulation.
    # However, for the real image, it loooks like this is hard to do
    # What should be the standard way to judge the PSSN in the real telescope?

    # OPD is zero for perfect telescope
    opdt = np.zeros((m, m))

    # OPD to PSF
    psft = opd2psf(
        opdt,
        iad,
        wlum,
        imagedelta=imagedelta,
        sensorFactor=1,
        fno=fno,
        debugLevel=debugLevel,
    )

    # PSF to optical transfer function (OTF)
    otft = psf2otf(psft)

    # Add atmosphere to perfect telescope
    otfa = otft * mtfa

    # OTF to PSF
    psfa = otf2psf(otfa)

    # Atmospheric PSS (point spread sensitivity) = 1/neff_atm
    pssa = np.sum(psfa**2)

    # Calculate PSF with error (atmosphere + system)
    if aType == "opd":
        if array.ndim == 2:
            ninst = 1
        else:
            ninst = array.shape[0]

        for ii in range(ninst):
            if array.ndim == 2:
                array2D = array
            else:
                array2D = array[ii, :, :].squeeze()

            psfei = opd2psf(array2D, iad, wlum, debugLevel=debugLevel)

            if ii == 0:
                psfe = psfei
            else:
                psfe += psfei

        # Do the normalization based on the number of instrument
        psfe = psfe / ninst

    elif aType == "psf":
        if array.shape[0] == mk:
            psfe = array

        elif array.shape[0] > mk:
            psfe = extractArray(array, mk)

        else:
            print(
                "calc_pssn: image provided too small, %d < %d x %6.4f."
                % (array.shape[0], m, k)
            )
            print("IQ is over-estimated !!!")
            psfe = padArray(array, mk)

        # Do the normalization of PSF
        psfe = psfe / np.sum(psfe) * np.sum(psft)

    # OTF with system error
    otfe = psf2otf(psfe)

    # Add the atmosphere error
    # OTF with system and atmosphere errors
    otftot = otfe * mtfa

    # PSF with system and atmosphere errors
    psftot = otf2psf(otftot)

    # atmospheric + error PSS
    pss = np.sum(psftot**2)

    # normalized PSS
    pssn = pss / pssa

    if debugLevel >= 3:
        print("pssn = %10.8e/%10.8e = %6.4f." % (pss, pssa, pssn))

    return pssn


def createMTFatm(D, m, k, wlum, zen, r0inmRef, model="vonK"):
    """Generate the modulation transfer function (MTF) for atmosphere.

    Parameters
    ----------
    D : float
        Side length of optical path difference (OPD) image in m.
    m : int
        Dimension of OPD image in pixel. The the number of pixel we want to
        have to cover the length of D.
    k : int
        Use a k-times bigger array to pad the MTF. Use k=1 for the same size.
    wlum : float
        Wavelength in um.
    zen : float
        Telescope zenith angle in degree.
    r0inmRef : float
        Reference r0 in meter at the wavelength of 0.5 um.
    model : str, optional
        Kolmogorov power spectrum ("Kolm") or van Karman power spectrum
        ("vonK"). (the default is "vonK".)

    Returns
    -------
    numpy.ndarray
        MTF at specific atmosphere model.
    """

    # Get the atmosphere phase structure function
    sfa = atmSF(D, m, wlum, zen, r0inmRef, model)

    # Get the modular transfer function for atmosphere
    mtfa = np.exp(-0.5 * sfa)

    # Add even number
    N = int(m + np.rint((m * (k - 1) + 1e-5) / 2) * 2)

    # Pad the matrix if necessary
    mtfa = padArray(mtfa, N)

    return mtfa


def atmSF(D, m, wlum, zen, r0inmRef, model):
    """Get the atmosphere phase structure function.

    Parameters
    ----------
    D : float
        Side length of optical path difference (OPD) image in m.
    m : int
        Dimension of OPD image in pixel.
    wlum : float
        Wavelength in um.
    zen : float
        Telescope zenith angle in degree.
    r0inmRef : float
        Reference r0 in meter at the wavelength of 0.5 um.
    model : str
        Kolmogorov power spectrum ("Kolm") or van Karman power spectrum
        ("vonK").

    Returns
    -------
    numpy.ndarray
        Atmosphere phase structure function.

    Raises
    ------
    ValueError
        The model type is not supported.
    """

    # Check the model
    if model not in ("Kolm", "vonK"):
        raise ValueError("Does not support %s atmosphere model." % model)

    # Get the atomosphere reference r0 in meter.
    r0a = r0Wz(r0inmRef, zen, wlum)

    # Round elements of the array to the nearest integer.
    m0 = np.rint(0.5 * (m + 1) + 1e-5)

    # Get the x, y coordinates index
    aa = np.arange(1, m + 1)
    x, y = np.meshgrid(aa, aa)

    # Frequency resolution in 1/rad
    dr = D / (m - 1)

    # Atmosphere r
    r = dr * np.sqrt((x - m0) ** 2 + (y - m0) ** 2)

    # Calculate the structure function

    # Kolmogorov power spectrum
    if model == "Kolm":
        # D(r) = 6.88 * (r/r0)^(5/3) in p.117, Chap. 11 of PhoSim referece
        sfa = 6.88 * (r / r0a) ** (5 / 3)

    # van Karman power spectrum
    elif model == "vonK":
        # Outer scale in meter
        L0 = 30

        # Gamma function is used
        sfa_c = (
            2
            * sp.gamma(11 / 6)
            / 2 ** (5 / 6)
            / np.pi ** (8 / 3)
            * (24 / 5 * sp.gamma(6 / 5)) ** (5 / 6)
            * (r0a / L0) ** (-5 / 3)
        )

        # Modified bessel of 2nd/3rd kind
        sfa_k = sp.kv(5 / 6, (2 * np.pi / L0 * r))

        # There is the undefined value (nan = 0 * inf, 0 from 'r' and inf from
        # 'sfa_k')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            sfa = sfa_c * (
                2 ** (-1 / 6) * sp.gamma(5 / 6)
                - (2 * np.pi / L0 * r) ** (5 / 6) * sfa_k
            )
            np.nan_to_num(sfa, copy=False)

    return sfa


def r0Wz(r0inmRef, zen, wlum):
    """Get the atomosphere reference r0, which is a function of zenith angle
    and wavelength.

    Parameters
    ----------
    r0inmRef : float
        Reference r0 in meter at the wavelength of 0.5 um.
    zen : float
        Telescope zenith angle in degree.
    wlum : float
        Wavelength in um.

    Returns
    -------
    float
        Atomosphere reference r0 in meter.
    """

    # Telescope zenith angle, change the unit from degree to radian
    zen = zen * np.pi / 180

    # Get the atmosphere reference r0
    r0aref = r0inmRef * np.cos(zen) ** 0.6

    # Atmosphere reference r0 at the specific wavelength in um
    # 0.5 um is the reference wavelength
    r0a = r0aref * (wlum / 0.5) ** 1.2

    return r0a


def psf2eAtmW(
    array,
    wlum,
    aType="opd",
    D=8.36,
    pmask=0,
    r0inmRef=0.1382,
    sensorFactor=1,
    zen=0,
    imagedelta=0.2,
    fno=1.2335,
    debugLevel=0,
):
    """Calculate the ellipticity with the error of atmosphere and weighting
    function.

    Parameters
    ----------
    array : numpy.ndarray
        Wavefront OPD in micron, or psf image.
    wlum : float
        Wavelength in microns.
    aType : str, optional
        Type of image ("opd" or "psf"). (the default is "opd".)
    D : float, optional
        Side length of optical path difference (OPD) image in m. (the default
        is 8.36.)
    pmask : int or numpy.ndarray[int], optional
        Pupil mask. (the default is 0.)
    r0inmRef : float, optional
        Fidicial atmosphere r0 @ 500nm in meter. (the default is 0.1382.)
    sensorFactor : float, optional
        Factor of sensor. (the default is 1.)
    zen : float, optional
        Telescope zenith angle in degree. (the default is 0.)
    imagedelta : float, optional
        Only needed when psf is used. 1 pixel = 0.2 arcsec. (the default is
        0.2.)
    fno : float, optional
        Only needed when psf is used. use 0 for opd. (the default is 1.2335.)
    debugLevel : int, optional
        The higher value gives more information. (the default is 0.)

    Returns
    -------
    float
        Ellipticity.
    numpy.ndarray
        Correlation function (XX).
    numpy.ndarray
        Correlation function (YY).
    numpy.ndarray
        Correlation function (XY).
    """

    # Unlike calc_pssn(), here imagedelta needs to be provided for type='opd'
    # because the ellipticity calculation operates on psf.

    # Get the k value
    k = fno * wlum / imagedelta

    # Get the PSF with the system error
    if aType == "opd":
        m = array.shape[0] / sensorFactor
        psfe = opd2psf(
            array,
            0,
            wlum,
            imagedelta=imagedelta,
            sensorFactor=sensorFactor,
            fno=fno,
            debugLevel=debugLevel,
        )
    else:
        m = max(pmask.shape)
        psfe = array

    # Opitcal transfer function (OTF) of system error
    otfe = psf2otf(psfe)

    # Modulation transfer function (MTF) with atmosphere
    mtfa = createMTFatm(D, m, k, wlum, zen, r0inmRef)

    # OTF with system and atmosphere errors
    otf = otfe * mtfa

    # PSF with system and atmosphere errors
    psf = otf2psf(otf)

    if debugLevel >= 3:
        print("Below from the Gaussian weigting function on ellipticity.")

    # Get the ellipticity and correlation function
    # The second input of psfeW should be pixeinum (1 pixel = 10 um).
    # Check this part with Bo.
    e, q11, q22, q12 = psf2eW(
        psf, imagedelta, wlum, atmModel="Gau", debugLevel=debugLevel
    )

    return e, q11, q22, q12


def psf2eW(psf, pixinum, wlum, atmModel="Gau", debugLevel=0):
    """Calculate the ellipticity with the weighting function.

    Parameters
    ----------
    psf : numpy.ndarray
        Point spread function (PSF).
    pixinum : float
        Pixel in um.
    wlum : float
        Wavelength in microns.
    atmModel : str, optional
        Atmosphere model ("Gau" or "2Gau"). (the default is "Gau".)
    debugLevel : int, optional
        The higher value gives more information. (the default is 0.)

    Returns
    -------
    float
        Ellipticity.
    numpy.ndarray
        Correlation function (XX).
    numpy.ndarray
        Correlation function (YY).
    numpy.ndarray
        Correlation function (XY).
    """

    # x, y positions
    x, y = np.meshgrid(np.arange(1, psf.shape[0] + 1), np.arange(1, psf.shape[1] + 1))

    # Average x and y
    xbar = np.sum(x * psf) / np.sum(psf)
    ybar = np.sum(y * psf) / np.sum(psf)

    # Show the averaged x and y
    if debugLevel >= 3:
        print("xbar=%6.3f, ybar=%6.3f" % (xbar, ybar))

    # Distance^2 to center
    r2 = (x - xbar) ** 2 + (y - ybar) ** 2

    # Weighting function based on the atmospheric model
    # FWHM is assigned to be 0.6 arcsec. Need to check with Bo for this.
    fwhminarcsec = 0.6
    oversample = 1
    W = createAtm(
        wlum,
        fwhminarcsec,
        r2,
        pixinum,
        oversample,
        model=atmModel,
        debugLevel=debugLevel,
    )

    # Apply the weighting function to PSF
    psf = psf * W

    # Correlation function
    Q11 = np.sum(((x - xbar) ** 2) * psf) / np.sum(psf)
    Q22 = np.sum(((y - ybar) ** 2) * psf) / np.sum(psf)
    Q12 = np.sum(((x - xbar) * (y - ybar)) * psf) / np.sum(psf)

    # Calculate the ellipticity
    T = Q11 + Q22
    if T > 1e-20:
        e1 = (Q11 - Q22) / T
        e2 = 2 * Q12 / T

        e = np.sqrt(e1**2 + e2**2)

    # No correlation
    else:
        e = 0

    return e, Q11, Q22, Q12


def createAtm(
    wlum, fwhminarcsec, gridsize, pixinum, oversample, model="Gau", debugLevel=0
):
    """Calculate the weighting function for a certain atmosphere model.

    Parameters
    ----------
    wlum : float
        Wavelength in microns.
    fwhminarcsec : float
        Full width in half maximum (FWHM) in arcsec.
    gridsize : int or numpy.ndarray[int]
        Size of grid. If it is the array, it should be (distance to center)^2.
        That means r2.
    pixinum : int
        Pixel in um.
    oversample : int
        k times of image resolution compared with the original one.
    model : str, optional
        Atmosphere model ("Gau" or "2Gau"). (the default is "Gau".)
    debugLevel : int, optional
        The higher value gives more information. (the default is 0.)

    Returns
    -------
    numpy.ndarray
        Weighting function of atmosphere.
    """

    # Get the weighting function

    # Distance^2 to center
    if isinstance(gridsize, (int)):
        nreso = gridsize * oversample

        # n for radius length
        nr = nreso / 2
        aa = np.linspace(-nr + 0.5, nr - 0.5, nreso)
        x, y = np.meshgrid(aa)

        r2 = x * x + y * y

    else:
        r2 = gridsize

    # FWHM in arcsec --> FWHM in um
    fwhminum = fwhminarcsec / 0.2 * 10

    # Calculate the weighting function
    if model == "Gau":
        # Sigma in um
        sig = fwhminum / 2 / np.sqrt(2 * np.log(2))
        sig = sig / (pixinum / oversample)

        z = np.exp(-r2 / 2 / sig**2)

    elif model == "2Gau":
        # Below is used to manually solve for sigma
        # let x = exp(-r^2/(2*alpha^2)), which results in 1/2*max
        # we want to get (1+.1)/2=0.55 from below
        # x=0.4673194304; printf('%20.10f\n'%x**.25*.1+x);
        sig = fwhminum / (2 * np.sqrt(-2 * np.log(0.4673194304)))

        # In (oversampled) pixel
        sig = sig / (pixinum / oversample)

        z = np.exp(-r2 / 2 / sig**2) + 0.4 / 4 * np.exp(-r2 / 8 / sig**2)

    if debugLevel >= 3:
        print("sigma1=%6.4f arcsec" % (sig * (pixinum / oversample) / 10 * 0.2))

    return z


def opd2psf(
    opd, pupil, wavelength, imagedelta=0, sensorFactor=1, fno=1.2335, debugLevel=0
):
    """Optical path difference (OPD) to point spread function (PSF).

    Parameters
    ----------
    opd : numpy.ndarray
        Optical path difference.
    pupil : float or numpy.ndarray
        Pupil function. If pupil is a number, not an array, we will get pupil
        geometry from OPD.
    wavelength : float
        Wavelength in um.
    imagedelta : float, optional
        Pixel size in um. Use 0 if pixel size is not specified. (the default
        is 0.)
    sensorFactor : float, optional
        Factor of sensor. Only need this if imagedelta != 0. (the default
        is 1.)
    fno : float, optional
         Only need this if imagedelta=0. (the default is 1.2335.)
    debugLevel : int, optional
        The higher value gives more information. (the default is 0.)

    Returns
    -------
    numpy.ndarray
        Normalized PSF.

    Raises
    ------
    ValueError
        Shapes of OPD and pupil are different.
    ValueError
        OPD shape is not square.
    ValueError
        Padding value is less than 1.
    """

    # Make sure all NaN in OPD to be 0
    opd[np.isnan(opd)] = 0

    # Get the pupil function from OPD if necessary
    if not isinstance(pupil, np.ndarray):
        pupil = opd != 0

    # Check the dimension of pupil and OPD should be the same
    if opd.shape != pupil.shape:
        raise ValueError("Shapes of OPD and pupil are different.")

    # For the PSF
    if imagedelta != 0:
        # Check the dimension of OPD
        if opd.shape[0] != opd.shape[1]:
            raise ValueError(
                "Error (opd2psf): OPD image size = (%d, %d)."
                % (opd.shape[0], opd.shape[1])
            )

        # Get the k value and the padding
        k = fno * wavelength / imagedelta
        padding = k / sensorFactor

        # Check the padding
        if padding < 1:
            errorMes = "opd2psf: Sampling too low, data inaccurate.\n"
            errorMes += "Imagedelta needs to be smaller than"
            errorMes += " fno * wlum = %4.2f um.\n" % (fno * wavelength)
            errorMes += "So that the padding factor > 1.\n"
            errorMes += "Otherwise we have to cut pupil to be < D."

            raise ValueError(errorMes)

        # Size of sensor
        sensorSamples = opd.shape[0]

        # Add even number for padding
        N = int(sensorSamples + np.rint(((padding - 1) * sensorSamples + 1e-5) / 2) * 2)
        pupil = padArray(pupil, N)
        opd = padArray(opd, N)

        # Show the padding information or not
        if debugLevel >= 3:
            print("padding = %8.6f." % padding)

    # If imagedelta = 0, we don't do any padding, and go with below
    z = pupil * np.exp(-2j * np.pi * opd / wavelength)
    z = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(z), s=z.shape))
    z = np.absolute(z**2)

    # Normalize the PSF
    z = z / np.sum(z)

    # Show the information of PSF from OPD
    if debugLevel >= 3:
        print("opd2psf(): imagedelta = %8.6f." % imagedelta, end="")

        if imagedelta == 0:
            print("0 means using OPD with padding as provided.")

        print("Verify psf has been normalized: %4.1f." % np.sum(z))

    return z


def psf2otf(psf):
    """Point spread function (PSF) to optical transfer function (OTF).

    Parameters
    ----------
    psf : numpy.ndarray
        Point spread function.

    Returns
    -------
    numpy.ndarray
        Optacal transfer function.
    """

    otf = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(psf), s=psf.shape))

    return otf


def otf2psf(otf):
    """Optical transfer function (OTF) to point spread function (PSF).

    Parameters
    ----------
    otf : numpy.ndarray
        Optical transfer function.

    Returns
    -------
    numpy.ndarray
        Point spread function.
    """

    psf = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(otf), s=otf.shape)))

    return psf


if __name__ == "__main__":
    pass
