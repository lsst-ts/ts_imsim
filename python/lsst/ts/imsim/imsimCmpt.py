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
import yaml
import numpy as np
from scipy.ndimage import rotate
from astropy.io import fits

from lsst.ts.wep.utility import runProgram
from lsst.ts.ofc.ofc_data.base_ofc_data import BaseOFCData
from lsst.ts.imsim.opdMetrology import OpdMetrology
from lsst.ts.imsim.utils.utility import getConfigDir, makeDir


class ImsimCmpt:
    def __init__(self):
        """Class to take configurations for each imsim component and
        generate full imsim configuration files.
        """
        # Output directories
        self._outputDir = None
        self._outputImgDir = None

        # OPD information
        self.opdFilePath = None
        self.opdMetr = OpdMetrology()

        # Specify number of Zernikes
        self.numOfZk = 19

        # AOS Degrees of Freedom
        self.numOfDof = 50
        self.dofInUm = np.zeros(self.numOfDof, dtype=float)

    @property
    def outputDir(self):
        return self._outputDir

    @outputDir.setter
    def outputDir(self, newOutputDir):
        """
        Set the closed loop output directory and make it if it
        does not yet exits.

        Parameters
        ----------
        newOutputDir : str
            Path for output image directory.
        """
        makeDir(newOutputDir)
        self._outputDir = newOutputDir

    @property
    def outputImgDir(self):
        return self._outputImgDir

    @outputImgDir.setter
    def outputImgDir(self, newOutputImgDir):
        """
        Set the output image directory and make it if it
        does not yet exits.

        Parameters
        ----------
        newOutputImgDir : str
            Path for output image directory.
        """
        makeDir(newOutputImgDir)
        self._outputImgDir = newOutputImgDir

    def _verifyPointerFile(self, filePointerInfo, configSections):
        """Verify that pointer file has filepaths for all needed
        sections of a complete imsim configuration file.

        Parameters
        ----------
        filePointerInfo : dict
            Dictionary holding pointer types as keys and file path
            information as values.
        configSections : list
            List of required submodules for ImSim configuration.

        Raises
        ------
        ValueError
            Raised when filePointerInfo is missing a required submodule.
        """
        filePointerKeys = list(filePointerInfo.keys())
        for sectionName in configSections:
            if sectionName not in filePointerKeys:
                raise ValueError(
                    f"Config pointer file missing filepath for {sectionName}."
                )

    def assembleConfigYaml(
        self,
        obsMetadata,
        configPointerFile,
        instName,
        requiredModulesFile=None,
    ):
        """
        Parameters
        ----------
        obsMetadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.
        configPointerFile : str
            The location of the config pointer file that specifies
            the location of configurations files for each imsim
            component.
        instName : str
            Instrument name.
        requiredModulesFile : str or None
            Path to yaml file with required imsim modules. If None,
            then will use policy/requiredModulesDefault.yaml.
            (The default is None.)
        """

        if requiredModulesFile is None:
            requiredModulesFile = os.path.join(
                getConfigDir(), "requiredModulesDefault.yaml"
            )

        obsInfoText = self.convertObsMetadataToText(obsMetadata)

        configSections = ["input", "gal", "image", "psf", "stamp", "output"]

        with open(configPointerFile, "r") as f:
            filePointerInfo = yaml.safe_load(f)
        self._verifyPointerFile(filePointerInfo, configSections)

        with open(requiredModulesFile, "r") as requiredModules:
            requiredModulesText = requiredModules.read()

        # Assemble full configuration as text because of the
        # aliases in the individual yaml components that
        # fill in values from "evalVariables".
        fullConfigText = ""
        fullConfigText += requiredModulesText + "\n"
        fullConfigText += obsInfoText + "\n"

        # Treat input section first differently since it has multiple parts
        inputSectionText = "input:\n"
        for subsection in ["atm_psf", "sky_model", "telescope", "vignetting"]:
            with open(
                filePointerInfo["input"][subsection].format(**os.environ), "r"
            ) as subFile:
                for line in subFile.readlines():
                    inputSectionText += "  " + line
        fullConfigText += inputSectionText + "\n"
        # Now move on to other sections skipping input section
        for sectionName in configSections:
            if sectionName == "input":
                continue
            else:
                with open(
                    filePointerInfo[sectionName].format(**os.environ), "r"
                ) as sectionFile:
                    fullConfigText += sectionFile.read() + "\n"
            # Handle output specially since we need to add header info and OPD
            if sectionName == "output":
                fullConfigText += f"  dir: {self.outputImgDir}\n\n"
                # Add additional header text
                fullConfigText += self.addConfigHeader(obsMetadata) + "\n"
                # Add OPD
                fullConfigText += self.formatOpdText(obsMetadata, instName)

        # Assemble as yaml
        fullConfigYaml = yaml.safe_load(fullConfigText)

        # Add in OFC corrections in mm
        dofInMm = self.dofInUm * 1e-3
        dofInMm[[3, 4, 8, 9]] /= 1e-3  # Undo unit change for items in arcseconds
        fullConfigYaml["input"]["telescope"]["fea"]["aos_dof"] = {
            "dof": self.dofInUm.tolist(),
            "type": "List",
        }
        return fullConfigYaml

    def convertObsMetadataToText(self, obsMetadata):
        """
        Write out the evalVariables section of the config file.

        Parameters
        ----------
        obsMetadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.

        Returns
        -------
        str
            evalVariables for ImSim config
        """
        obsVariablesText = "eval_variables:\n"
        obsVariablesText += "  cboresight:\n"
        obsVariablesText += "    type: RADec\n"
        obsVariablesText += f"    ra: &ra {obsMetadata.ra} deg\n"
        obsVariablesText += f"    dec: &dec {obsMetadata.dec} deg\n"
        obsVariablesText += f"  sband: &band {obsMetadata.band}\n"
        obsVariablesText += f"  azenith: &zenith {obsMetadata.zenith} deg\n"
        obsVariablesText += f"  artp: &rtp {obsMetadata.rotatorAngle} deg\n"
        obsVariablesText += f"  fexptime: &exptime {obsMetadata.expTime}\n"
        obsVariablesText += f"  fmjd: &mjd {obsMetadata.mjd}\n"
        obsVariablesText += f"  iseqnum: &seqnum {obsMetadata.seqNum}\n"
        obsVariablesText += f"  sobsid: &obsid {obsMetadata.obsId}\n"

        return obsVariablesText

    def formatOpdText(self, obsMetadata, instName):
        """
        Write out the OPD section of the config file.

        Parameters
        ----------
        obsMetadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.
        instName: str
            Name of the instrument

        Returns
        -------
        str
            OPD information for ImSim config
        """
        opdSectionText = "  opd:\n"
        opdSectionText += "    file_name: opd.fits\n"
        opdSectionText += "    nx: 255\n"
        opdSectionText += "    rotTelPos: *rtp\n"
        opdSectionText += "    jmax: 28\n"
        opdSectionText += "    eps: 0.61\n"
        opdSectionText += "    projection: gnomonic\n"
        opdSectionText += f"    wavelength: {BaseOFCData().eff_wavelength[obsMetadata.band.upper()]*1e3}\n"
        opdSectionText += "    fields:\n"

        # Get the locations for the OPD from OPD Metrology
        if instName == "lsst":
            self.opdMetr.setDefaultLsstWfsGQ()
        else:
            self.opdMetr.setWgtAndFieldXyOfGQ(instName)
        for thx, thy in zip(self.opdMetr.fieldX, self.opdMetr.fieldY):
            opdSectionText += f"      - {{thx: {thx} deg, thy: {thy} deg}}\n"

        return opdSectionText

    def addConfigHeader(self, obsMetadata):
        """
        Write out the header section of the config file.

        Parameters
        ----------
        obsMetadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.

        Returns
        -------
        str
            Header information for ImSim config
        """
        headerText = "  header:\n"
        headerText += f"    mjd: {obsMetadata.mjd}\n"
        headerText += f"    observationStartMJD: {obsMetadata.mjd - (15/(60*60*24))}\n"
        headerText += f"    seqnum: {obsMetadata.seqNum}\n"
        headerText += f"    band: {obsMetadata.band}\n"
        headerText += f"    fieldRA: {obsMetadata.ra}\n"
        headerText += f"    fieldDec: {obsMetadata.dec}\n"
        headerText += f"    rotTelPos: {obsMetadata.rotatorAngle}\n"
        headerText += f"    airmass: {1.0/np.cos(obsMetadata.zenith)}\n"

        return headerText

    def runImsim(self, configFilePath):
        runProgram(f"galsim {configFilePath}")

    def genInstanceCatalog(self, skySim):
        """Generate the instance catalog."""
        # bandDict = {x[1]: x[0] for x in list(enumerate(['u', 'g', 'r', 'i', 'z', 'y']))}

        content = ""
        # content += f"mjd {obsMetadata.mjd}\n"
        # content += f"filter {bandDict[obsMetadata.band]}\n"
        # content += f"rightascension {obsMetadata.ra}\n"
        # content += f"declination {obsMetadata.dec}\n"
        # content += f"rotTelPos {obsMetadata.rotatorAngle}\n"
        # content += f"vistime 30.0\n"
        content += self.genInstCatStars(skySim)
        return content

    def genInstCatStars(self, skySim):
        content = ""
        for id, ra, dec, mag in zip(skySim.starId, skySim.ra, skySim.dec, skySim.mag):
            content += self.generateStar(id, ra, dec, mag)

        return content

    def generateStar(
        self,
        starId,
        ra,
        dec,
        magNorm,
        sedName="flatSED/sed_flat.txt.gz",
        redshift=0,
        gamma1=0,
        gamma2=0,
        kappa=0,
        deltaRa=0,
        deltaDec=0,
        sourceType="point",
    ):
        """Generate the star source.

        Parameters
        ----------
        starId : int
            Star Id.
        ra : float
            The right ascension of the center of the object or image in
            decimal degrees.
        dec : float
            The declination of the center of the object in decimal degrees.
        magNorm : float
            The normalization of the flux of the object in AB magnitudes
            at (500 nm)/(1+z) (which is roughly equivalent to V (AB) or
            g (AB)).
        sedName : str, optional
            The name of the SED file with a file path that is relative to the
            data directory in PhoSim.
            (The default is "flatSED/sed_flat.txt.gz")
        redshift : float, optional
            The redshift (or blueshift) of the object. Note that the SED does
            not need to be redshifted if using this. (the default is 0.)
        gamma1 : float, optional
            The value of the shear parameter gamma1 used in weak lensing.
            (the default is 0.)
        gamma2 : float, optional
            The value of the shear parameter gamma2 used in weak lensing.
            (the default is 0.)
        kappa : float, optional
            The value of the magnification parameter in weak lensing. (the
            default is 0.)
        deltaRa : float, optional
            The value of the declination offset in radians. This can be used
            either for weak lensing or objects that moved from another
            exposure if you do not want to change the source position in the
            first two columns. (the default is 0.)
        deltaDec : float, optional
            The value of the declination offset in radians. This can be used
            either for weak lensing or objects that moved from another
            exposure if you do not want to change the source position in the
            first two columns. (the default is 0.)
        sourceType : str, optional
            The name of the spatial model to be used as defined below. (the
            default is "point".)

        Returns
        -------
        str
            Perturbation command used in PhoSim.
        """

        content = "object %2d\t%9.6f\t%9.6f %9.6f %s " % (
            starId,
            ra,
            dec,
            magNorm,
            sedName,
        )
        content += "%.1f %.1f %.1f %.1f %.1f %.1f %s none none \n" % (
            redshift,
            gamma1,
            gamma2,
            kappa,
            deltaRa,
            deltaDec,
            sourceType,
        )

        return content

    def analyzeOpdData(
        self, instName, zkFileName="opd.zer", rotOpdInDeg=0.0, pssnFileName="PSSN.txt"
    ):
        """Analyze the OPD data.

        Rotate OPD to simulate the output by rotated camera. When anaylzing the
        PSSN, the unrotated OPD is used.

        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        instName : `str`
            Instrument name.
        zkFileName : str, optional
            OPD in zk file name. (the default is "opd.zer".)
        rotOpdInDeg : float, optional
            Rotate OPD in degree in the counter-clockwise direction. (the
            default is 0.0.)
        pssnFileName : str, optional
            PSSN file name. (the default is "PSSN.txt".)
        """

        # Set the weighting ratio and field positions of OPD
        if instName == "lsst":
            self.opdMetr.setDefaultLsstWfsGQ()
        else:
            self.opdMetr.setWgtAndFieldXyOfGQ(instName)
        numOpd = len(self.opdMetr.fieldX)

        self.opdFilePath = os.path.join(self.outputImgDir, "opd.fits")
        self._writeOpdZkFile(zkFileName, rotOpdInDeg, numOpd)
        self._writeOpdPssnFile(instName, pssnFileName, numOpd)

    def _writeOpdZkFile(self, zkFileName, rotOpdInDeg, numOpd):
        """Write the OPD in zk file.

        OPD: optical path difference.

        Parameters
        ----------
        zkFileName : str
            OPD in zk file name.
        rotOpdInDeg : float
            Rotate OPD in degree in the counter-clockwise direction.
        numOpd : int
            Number of OPD positions calculated.
        """

        filePath = os.path.join(self.outputImgDir, zkFileName)
        opdData = self._mapOpdToZk(rotOpdInDeg, numOpd)
        fileTxt = (
            "# The followings are OPD in rotation angle of %.2f degree in nm from z4 to z22:\n"
            % (rotOpdInDeg)
        )
        for sensorId, opdZk in zip(self.opdMetr.sensorIds, opdData):
            zkStr = f"{sensorId}: {opdZk}\n"
            fileTxt += zkStr
        with open(filePath, "w") as file:
            file.write(fileTxt)

    def _mapOpdToZk(self, rotOpdInDeg, numOpd):
        """Map the OPD to the basis of annular Zernike polynomial (Zk).

        OPD: optical path difference.

        Parameters
        ----------
        rotOpdInDeg : float
            Rotate OPD in degree in the counter-clockwise direction.
        numOpd : int
            Number of OPD positions calculated.

        Returns
        -------
        numpy.ndarray
            Zk data from OPD. This is a 2D array. The row is the OPD index and
            the column is z4 to z22 in um. The order of OPD index is based on
            the file name.
        """

        # Map the OPD to the Zk basis and do the collection
        # Get the number of OPD locations by looking at length of fieldX
        opdData = np.zeros((numOpd, self.numOfZk))
        for idx in range(numOpd):
            opd = fits.getdata(self.opdFilePath, idx)

            # Rotate OPD if needed
            if rotOpdInDeg != 0:
                opdRot = rotate(opd, rotOpdInDeg, reshape=False)
                opdRot[opd == 0] = 0
            else:
                opdRot = opd

            # z1 to z22 (22 terms)
            zk = self.opdMetr.getZkFromOpd(opdMap=opdRot)[0]

            # Only need to collect z4 to z22
            initIdx = 3
            opdData[idx, :] = zk[initIdx : initIdx + self.numOfZk]

        return opdData

    def _writeOpdPssnFile(self, instName, pssnFileName, numOpd):
        """Write the OPD PSSN in file.

        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        instName : `str`
            Instrument name.
        pssnFileName : str
            PSSN file name.
        numOpd : int
            Number of OPD positions calculated.
        """

        # Calculate the PSSN
        pssnList, gqEffPssn = self._calcPssnOpd(numOpd)

        # Calculate the FWHM
        effFwhmList, gqEffFwhm = self._calcEffFwhmOpd(pssnList)

        # Append the list to write the data into file
        pssnList.append(gqEffPssn)
        effFwhmList.append(gqEffFwhm)

        # Stack the data
        data = np.vstack((pssnList, effFwhmList))

        # Write to file
        filePath = os.path.join(self.outputImgDir, pssnFileName)
        header = "The followings are PSSN and FWHM (in arcsec) data. The final number is the GQ value."
        np.savetxt(filePath, data, header=header)

    def _calcPssnOpd(self, numOpd):
        """Calculate the PSSN of OPD.

        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.
        GQ: Gaussian quadrature.

        Parameters
        ----------
        numOpd : int
            Number of OPD positions calculated.

        Returns
        -------
        list
            PSSN list.
        float
            GQ effective PSSN.
        """

        pssnList = []
        for idx in range(numOpd):
            wavelengthInUm = fits.getheader(self.opdFilePath, idx)["WAVELEN"] * 1e-3
            # OPD data needs to be in microns. Imsim output is nm.
            pssn = self.opdMetr.calcPSSN(
                wavelengthInUm, opdMap=fits.getdata(self.opdFilePath, idx) * 1e-3
            )
            pssnList.append(pssn)

        # Calculate the GQ effectice PSSN
        gqEffPssn = self.opdMetr.calcGQvalue(pssnList)

        return pssnList, gqEffPssn

    def _calcEffFwhmOpd(self, pssnList):
        """Calculate the effective FWHM of OPD.

        FWHM: Full width and half maximum.
        PSSN: Normalized point source sensitivity.
        GQ: Gaussian quadrature.

        Parameters
        ----------
        pssnList : list
            List of PSSN.

        Returns
        -------
        list
            Effective FWHM list.
        float
            GQ effective FWHM.
        """

        # Calculate the list of effective FWHM
        effFwhmList = []
        for pssn in pssnList:
            effFwhm = self.opdMetr.calcFWHMeff(pssn)
            effFwhmList.append(effFwhm)

        # Calculate the GQ effectice FWHM
        gqEffFwhm = self.opdMetr.calcGQvalue(effFwhmList)

        return effFwhmList, gqEffFwhm

    def mapOpdDataToListOfWfErr(self, opdZkFileName, sensorIdList, sensorNameList):
        """Map the OPD data to the list of wavefront error.

        OPD: Optical path difference.

        Parameters
        ----------
        opdZkFileName : str
            OPD zk file name.
        sensorIdList : list
            Reference sensor ID list.
        sensorNameList : list
            Reference sensor name list.

        Returns
        -------
        list [lsst.ts.wep.ctrlIntf.SensorWavefrontError]
            List of SensorWavefrontError object.
        """

        opdZk = self._getZkFromFile(opdZkFileName)

        listOfWfErr = []
        for sensorId, sensorName, zk in zip(sensorIdList, sensorNameList, opdZk):
            sensorWavefrontData = SensorWavefrontError(numOfZk=self.getNumOfZk())
            sensorWavefrontData.setSensorId(sensorId)
            sensorWavefrontData.setSensorName(sensorName)
            sensorWavefrontData.setAnnularZernikePoly(zk)

            listOfWfErr.append(sensorWavefrontData)

        return listOfWfErr

    def _getZkFromFile(self, zkFileName):
        """Get the zk (z4-z22) from file.

        Parameters
        ----------
        zkFileName : str
            Zk file name.

        Returns
        -------
        numpy.ndarray
            zk matrix. The colunm is z4-z22. The raw is each data point.
        """

        filePath = os.path.join(self.outputImgDir, zkFileName)
        zk = np.loadtxt(filePath)

        return zk

    def getOpdGqEffFwhmFromFile(self, pssnFileName):
        """Get the OPD GQ effective FWHM from file.

        OPD: Optical path difference.
        GQ: Gaussian quadrature.
        FWHM: Full width at half maximum.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssnFileName : str
            PSSN file name.

        Returns
        -------
        float
            OPD GQ effective FWHM.
        """

        data = self._getDataOfPssnFile(pssnFileName)
        gqEffFwhm = data[1, -1]

        return gqEffFwhm

    def getListOfFwhmSensorData(self, pssnFileName, sensorIdList):
        """Get the list of FWHM sensor data based on the OPD PSSN file.

        FWHM: Full width at half maximum.
        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssnFileName : str
            PSSN file name.
        sensorIdList : list
            Reference sensor id list.

        Returns
        -------
        fwhmCollection : `np.ndarray [object]`
            Numpy array with fwhm data. This is a numpy array of arrays. The
            data type is `object` because each element may have different
            number of elements.
        sensor_id: `np.ndarray`
            Numpy array with sensor ids.
        """

        # Get the FWHM data from the PSSN file
        # The first row is the PSSN and the second one is the FWHM
        # The final element in each row is the GQ value
        data = self._getDataOfPssnFile(pssnFileName)
        fwhmData = data[1, :-1]

        sensor_id = np.array(sensorIdList, dtype=int)

        fwhmCollection = np.array([], dtype=object)
        for fwhm in fwhmData:
            fwhmCollection = np.append(fwhmCollection, fwhm)

        return fwhmCollection, sensor_id

    def getOpdPssnFromFile(self, pssnFileName):
        """Get the OPD PSSN from file.

        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssnFileName : str
            PSSN file name.

        Returns
        -------
        numpy.ndarray
            PSSN.
        """

        data = self._getDataOfPssnFile(pssnFileName)
        pssn = data[0, :-1]

        return pssn

    def _getDataOfPssnFile(self, pssnFileName):
        """Get the data of the PSSN file.

        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssnFileName : str
            PSSN file name.

        Returns
        -------
        numpy.ndarray
            Data of the PSSN file.
        """

        filePath = os.path.join(self.outputImgDir, pssnFileName)
        data = np.loadtxt(filePath)

        return data

    def reorderAndSaveWfErrFile(
        self, listOfWfErr, refSensorNameList, lsstCamera, zkFileName="wfs.zer"
    ):
        """Reorder the wavefront error in the wavefront error list according to
        the reference sensor name list and save to a file.

        The unexisted wavefront error will be a numpy zero array. The unit is
        um.

        Parameters
        ----------
        listOfWfErr : list [lsst.ts.wep.ctrlIntf.SensorWavefrontData]
            List of SensorWavefrontData object.
        refSensorNameList : list
            Reference sensor name list.
        lsstCamera : lsst.afw.cameraGeom.Camera
            Lsst instrument.
        zkFileName : str, optional
            Wavefront error file name. (the default is "wfs.zer".)
        """

        # Get the sensor name that in the wavefront error map
        wfErrMap = self._transListOfWfErrToMap(listOfWfErr, lsstCamera)
        nameListInWfErrMap = list(wfErrMap.keys())

        # Reorder the wavefront error map based on the reference sensor name
        # list.
        reorderedWfErrMap = dict()
        for sensorName in refSensorNameList:
            if sensorName in nameListInWfErrMap:
                wfErr = wfErrMap[sensorName]
            else:
                wfErr = np.zeros(self.numOfZk)
            reorderedWfErrMap[sensorName] = wfErr

        # Save the file
        filePath = os.path.join(self.outputImgDir, zkFileName)
        fileTxt = "# The followings are ZK in um from z4 to z22:\n"
        for key, val in reorderedWfErrMap.items():
            zkStr = f"{lsstCamera[key].getId()}: {val}\n"
            fileTxt += zkStr
        with open(filePath, "w") as file:
            file.write(fileTxt)

    def _transListOfWfErrToMap(self, listOfWfErr, lsstCamera):
        """Transform the list of wavefront error to map.

        Parameters
        ----------
        listOfWfErr : list [lsst.ts.wep.ctrlIntf.SensorWavefrontData]
            List of SensorWavefrontData object.
        lsstCamera : lsst.afw.cameraGeom.Camera
            Lsst instrument.

        Returns
        -------
        dict
            Calculated wavefront error. The dictionary key [str] is the
            abbreviated sensor name (e.g. R22_S11). The dictionary item
            [numpy.ndarray] is the averaged wavefront error (z4-z22) in um.
        """

        mapSensorNameAndId = dict(
            [(detector.getId(), detector.getName()) for detector in lsstCamera]
        )

        wfErrMap = dict()
        for sensorWavefrontData in listOfWfErr:
            sensorId = sensorWavefrontData.sensorId
            sensorName = mapSensorNameAndId[sensorId]

            avgErrInUm = sensorWavefrontData.annularZernikePoly

            wfErrMap[sensorName] = avgErrInUm

        return wfErrMap

    def _getWfErrValuesAndStackToMatrix(self, wfErrMap):
        """Get the wavefront errors and stack them to be a matrix.

        Parameters
        ----------
        wfErrMap : dict
            Calculated wavefront error. The dictionary key [str] is the
            abbreviated sensor name (e.g. R22_S11). The dictionary item
            [numpy.ndarray] is the averaged wavefront error (z4-z22) in um.

        Returns
        -------
        numpy.ndarray
            Wavefront errors as a matrix. The column is z4-z22 in um. The row
            is the individual sensor. The order is the same as the input of
            wfErrMap.
        """

        valueMatrix = np.empty((0, self.numOfZk))
        for wfErr in wfErrMap.values():
            valueMatrix = np.vstack((valueMatrix, wfErr))

        return valueMatrix

    def saveDofInUmFileForNextIter(self, dofInUmFileName="dofPertInNextIter.mat"):
        """Save the DOF in um data to file for the next iteration.

        DOF: degree of freedom.

        Parameters
        ----------
        dofInUmFileName : str, optional
            File name to save the DOF in um. (the default is
            "dofPertInNextIter.mat".)
        """

        filePath = os.path.join(self.outputDir, dofInUmFileName)
        header = "The followings are the DOF in um:"
        np.savetxt(filePath, np.transpose(self.dofInUm), header=header)
