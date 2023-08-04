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

import astropy
import numpy as np
from astroplan import Observer
from lsst.ts.imsim.utils.utility import getCamera


class SkySim:
    def __init__(self):
        """Initialization of sky simulator class."""

        # Star ID
        self.starId = np.array([], dtype=int)

        # Star RA
        self.ra = np.array([])

        # Star dec
        self.dec = np.array([])

        # Star magnitude
        self.mag = np.array([])

        # Camera
        self._camera = None

    def setCamera(self, instName):
        """Set the camera.

        Parameters
        ----------
        instName : `str`
            Instrument name. Valid options are 'comcam or 'lsstfam'.
        """

        self._camera = getCamera(instName)

    def calcParallacticAngle(self, obsMetadata):
        """Calculate the parallactic angle so we know the
        sky rotation angle on alt-az mount for the observation.

        Parameters
        ----------
        obsMetadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.

        Returns
        -------
        float
            Parallactic Angle in degrees.
        """
        time = astropy.time.Time(obsMetadata.mjd, format="mjd")
        rubin = Observer.at_site("cerro pachon")
        boresight = astropy.coordinates.SkyCoord(
            f"{obsMetadata.ra}d", f"{obsMetadata.dec}d", frame="icrs"
        )

        return rubin.parallactic_angle(time, boresight).deg

    def addStarByRaDecInDeg(self, starId, raInDeg, decInDeg, mag):
        """Add the star information by (ra, dec) in degrees.

        Parameters
        ----------
        starId : int, list[int], or numpy.ndarray[int]
            Star Id.
        raInDeg : float, list, or numpy.ndarray
            Star ra in degree.
        decInDeg : float, list, or numpy.ndarray
            Star dec in degree.
        mag : float, list, or numpy.ndarray
            Star magnitude.
        """

        # Check the inputs are list or not, and change the type if necessary
        starIdList = self._changeToListIfNecessary(starId)
        raInDegList = self._changeToListIfNecessary(raInDeg)
        decInDegList = self._changeToListIfNecessary(decInDeg)
        magList = self._changeToListIfNecessary(mag)

        # Add the stars
        for ii in range(len(starIdList)):
            intStarId = int(starIdList[ii])
            if self._isUniqStarId(intStarId):
                self.starId = np.append(self.starId, intStarId)
                self.ra = np.append(self.ra, raInDegList[ii])
                self.dec = np.append(self.dec, decInDegList[ii])
                self.mag = np.append(self.mag, magList[ii])

    def _changeToListIfNecessary(self, variable):
        """Change the data type to list.

        Parameters
        ----------
        variable : int, float, list, or numpy.ndarray
            Variable.

        Returns
        -------
        list
            Variable as the list.
        """

        if isinstance(variable, (int, float)):
            return [variable]
        else:
            return variable

    def _isUniqStarId(self, starId):
        """Check the star ID is unique or not.

        Parameters
        ----------
        starId : int
            Star Id.

        Returns
        -------
        bool
            True if the unique Id.
        """

        if starId in self.starId:
            isUnique = False
            print("StarId=%d is not unique." % starId)
        else:
            isUnique = True

        return isUnique

    def addStarByFile(self, readFilePath, skiprows=0):
        """Add the star data by reading the file.

        Parameters
        ----------
        readFilePath : str
            Star data file path.
        skiprows : int, optional
            Skip the first "skiprows" lines. (the default is 0.)
        """

        data = np.loadtxt(readFilePath, skiprows=skiprows)

        # Only consider the non-empty data
        if len(data) != 0:
            # Change to 2D array if the input is 1D array
            if data.ndim == 1:
                data = np.expand_dims(data, axis=0)

            for star in data:
                self.addStarByRaDecInDeg(star[0], star[1], star[2], star[3])

    def exportSkyToFile(self, outputFilePath):
        """Export the star information into the file.

        Parameters
        ----------
        outputFilePath : str
            Output file path.
        """

        # Add the header (star ID, ra, dec, magnitude)
        content = "# Id\t Ra\t\t Dec\t\t Mag\n"

        # Add the star information
        for ii in range(len(self.starId)):
            content += "%d\t %3.6f\t %3.6f\t %3.6f\n" % (
                self.starId[ii],
                self.ra[ii],
                self.dec[ii],
                self.mag[ii],
            )

        # Write into file
        fid = open(outputFilePath, "w")
        fid.write(content)
        fid.close()
