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

import logging
import os
import shutil

import astropy
import numpy as np
from copy import deepcopy
from lsst.afw.cameraGeom import FIELD_ANGLE, DetectorType
from lsst.daf import butler as dafButler
from lsst.ts.imsim.imsimCmpt import ImsimCmpt
from lsst.ts.imsim.obsMetadata import ObsMetadata
from lsst.ts.imsim.opdMetrology import OpdMetrology
from lsst.ts.imsim.skySim import SkySim
from lsst.ts.imsim.utils.plotUtil import plotFwhmOfIters
from lsst.ts.imsim.utils.sensorWavefrontError import SensorWavefrontError
from lsst.ts.imsim.utils.utility import getCamera, getConfigDir, makeDir
from lsst.ts.ofc import OFC, OFCData
from lsst.ts.wep.utils import CamType, FilterType, rotMatrix, runProgram
from lsst.ts.wep.utils import getConfigDir as getWepConfigDir


class ClosedLoopTask:
    def __init__(self):
        """Initilization of the closed loop task class to
        run the simulation with imSim."""

        self.log = logging.getLogger(type(self).__name__)

        # Sky simulator
        self.skySim = None

        # OFC calculator
        self.ofcCalc = None

        # imSim Component
        self.imsimCmpt = None

        # OPD Metrology
        self.opdMetr = None

        # Ra/Dec/RotAng coordinates used in the simulation.
        self.boresightRa = None
        self.boresightDec = None
        self.boresightRotAng = None

        # Use CCD image
        self.useCcdImg = True

    def configSkySim(self, instName, obsMetadata, pathSkyFile="", starMag=15):
        """Configure the sky simulator.

        If the path of sky file is not provided, The defult OPD field positions
        will be used.

        OPD: Optical path difference.

        Parameters
        ----------
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.
        obsMetadata : lsst.ts.imsim.ObsMetadata object
            Observation metadata.
        pathSkyFile : str, optional
            Path to the sky file. (the default is "".)
        starMag : float, optional
            Default star magnitude if there is no sky file used. This is to
            pretend there are the stars at OPD field positions. (the default is
            15.)

        Raises
        ------
        ValueError
            This instrument name is not supported.
        """

        self.skySim = SkySim()
        self.skySim.setCamera(instName)
        if pathSkyFile == "":
            parAngle = self.skySim.calcParallacticAngle(obsMetadata)
            self._setSkySimBasedOnOpdFieldPos(instName, obsMetadata, parAngle, starMag)
        else:
            absSkyFilePath = os.path.abspath(pathSkyFile)
            self.skySim.addStarByFile(absSkyFilePath)

    def _setSkySimBasedOnOpdFieldPos(self, instName, obsMetadata, parAngle, starMag):
        """Set the sky simulator based on the OPD field positions.

        OPD: Optical path difference.

        Parameters
        ----------
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.
        obsMetadata : lsst.ts.imsim.ObsMetadata object
            Observation metadata.
        parAngle : float
            Parallactic angle.
        starMag : float
            Star magnitude. This is to pretend there are the stars at OPD field
            positions.

        Raises
        ------
        ValueError
            This instrument name is not supported.
        """

        self.log.info(
            "Use the default OPD field positions to be star positions. "
            f"The star magnitude is chosen to be {starMag}."
        )

        opdMetr = OpdMetrology()
        if instName in ["lsst", "lsstfam"]:
            fieldX, fieldY = list(), list()
            camera = getCamera(instName)
            for name in self.getSensorNameListOfFields(instName):
                detector = camera.get(name)
                xRad, yRad = detector.getCenter(FIELD_ANGLE)
                xDeg, yDeg = np.rad2deg(xRad), np.rad2deg(yRad)
                fieldY.append(xDeg)  # transpose for imSim
                fieldX.append(yDeg)
            opdMetr.fieldX = fieldX
            opdMetr.fieldY = fieldY
        else:
            raise ValueError(f"This instrument name ({instName}) is not supported.")

        starId = 0
        raInDegArr = np.array(opdMetr.fieldX)
        decInDegArr = np.array(opdMetr.fieldY)
        rotation = rotMatrix(obsMetadata.rotatorAngle - parAngle)
        for raInDeg, decInDeg in zip(raInDegArr, decInDegArr):
            # It is noted that the field position might be < 0. But it is
            # not the same case for ra (0 <= ra <= 360).
            raInDeg, decInDeg = np.dot(rotation, np.array([raInDeg, decInDeg]))
            raInDeg += obsMetadata.ra
            decInDeg += obsMetadata.dec
            if raInDeg < 0:
                raInDeg += 360.0
            self.skySim.addStarByRaDecInDeg(starId, raInDeg, decInDeg, starMag)
            starId += 1

    def configOfcCalc(self, instName):
        """Configure the OFC calculator.

        OFC: Optical feedback calculator.

        Parameters
        ----------
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.
        """

        ofc_data = OFCData(instName)
        ofc_data.xref = "x00"
        self.ofcCalc = OFC(ofc_data)

    def mapFilterRefToG(self, filterTypeName):
        """Map the reference filter to the G filter.

        Parameters
        ----------
        filterTypeName : str
            Filter type name: ref (or ''), u, g, r, i, z, or y.

        Returns
        -------
        filterTypeName : str
            Mapped filter type.
        """
        return "g" if filterTypeName in ("ref", "") else filterTypeName

    def checkBoresight(self, boresight):
        """Check the boresight.

        Parameters
        ----------
        boresight : list[float]
            Boresight [ra, dec] in degree.

        Raises
        ------
        ValueError
            The right ascension (RA) should be in [0, 360].
        ValueError
            The declination (Dec) should be in [-90, 90].
        """

        ra, dec = boresight
        if ra < 0 or ra > 360:
            raise ValueError("The right ascension (RA) should be in [0, 360].")

        if dec < -90 or dec > 90:
            raise ValueError("The declination (Dec) should be in [-90, 90].")

    def getSensorNameListOfFields(self, instName):
        """Get the list of sensor name of fields.

        The list will be sorted based on the field index.

        Parameters
        ----------
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.

        Returns
        -------
        list[str]
            List of sensor name.

        Raises
        ------
        ValueError
            This instrument name is not supported.
        """

        camera = getCamera(instName)
        detectorType = (
            DetectorType.WAVEFRONT if instName == "lsst" else DetectorType.SCIENCE
        )
        return [
            detector.getName()
            for detector in camera
            if detector.getType() == detectorType
        ]

    def getSensorIdListOfFields(self, instName):
        """Get the list of sensor ids of fields.

        The list will be sorted based on the field index.

        Parameters
        ----------
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.

        Returns
        -------
        list[int]
            List of sensor ids.

        Raises
        ------
        ValueError
            This instrument name is not supported.
        """

        camera = getCamera(instName)

        detectorType = (
            DetectorType.WAVEFRONT if instName == "lsst" else DetectorType.SCIENCE
        )
        return [
            detector.getId()
            for detector in camera
            if detector.getType() == detectorType
        ]

    def checkAndCreateBaseOutputDir(self, baseOutputDir):
        """Check and create the base output directory.

        This function will create the directory if it does not exist.

        Parameters
        ----------
        baseOutputDir : str
            Base output directory.

        Returns
        -------
        str
            Base output directory.
        """
        outputDir = baseOutputDir
        makeDir(outputDir, exist_ok=True)

        return outputDir

    def getCamTypeAndInstName(self, inst):
        """Get the camera type and instrument name.

        Parameters
        ----------
        inst : str
            Instrument to use: currently only lsst.

        Returns
        -------
        camType : enum 'CamType' in lsst.ts.wep.utility
            Camera type.
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.

        Raises
        ------
        ValueError
            This instrument is not supported.
        """

        if inst == "lsst":
            return CamType.LsstCam, "lsst"
        elif inst == "lsstfam":
            return CamType.LsstFamCam, "lsstfam"
        else:
            raise ValueError(f"This instrument ({inst}) is not supported.")

    def getFilterType(self, filterTypeName):
        """Get the filter type.

        Parameters
        ----------
        filterTypeName : str
            Filter type name: ref, u, g, r, i, z, or y.

        Returns
        -------
        filterType : enum 'FilterType' in lsst.ts.wep.utility
            Filter type.

        Raises
        ------
        ValueError
            This filter type is not supported.
        """

        if filterTypeName in {"", "ref"}:
            return FilterType.REF
        elif filterTypeName == "u":
            return FilterType.LSST_U
        elif filterTypeName == "g":
            return FilterType.LSST_G
        elif filterTypeName == "r":
            return FilterType.LSST_R
        elif filterTypeName == "i":
            return FilterType.LSST_I
        elif filterTypeName == "z":
            return FilterType.LSST_Z
        elif filterTypeName == "y":
            return FilterType.LSST_Y
        else:
            raise ValueError(f"This filter type ({filterTypeName}) is not supported.")

    def _runSim(
        self,
        camType,
        instName,
        obsMetadata,
        baseOutputDir,
        butlerRootPath,
        skySeed,
        pertSeed,
        iterNum,
        numPro=1,
        pipelineFile="",
        imsimConfigPointerFile="",
        turnOffSkyBackground=False,
        turnOffAtmosphere=False,
    ):
        """Run the simulation.

        Parameters
        ----------
        camType : enum 'CamType' in lsst.ts.wep.utility
            Camera type.
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.
        obsMetadata : lsst.ts.imsim.ObsMetadata object
            Observation metadata.
        baseOutputDir : str
            Base output directory.
        butlerRootPath : str
            Path to the butler gen 3 repository.
        skySeed : int
            Random seed for the sky background.
        pertSeed : int
            Random seed for the perturbations.
        iterNum : int
            Number of closed-loop iteration.
        numPro : int, optional
            Number of processors to use. (The default is 1.)
        pipelineFile : str, optional
            Path to existing pipeline yaml file to use. If empty string
            then the code will write its own default pipeline yaml.
            (The default is "".)
        imsimConfigPointerFile : str, optional
            Path to pointer file with locations of yaml configuration
            files for imsim submodules. If empty string then the code
            will use the default in policy/config for the given inst.
            (The default is "".)
        turnOffSkyBackground : bool, optional
            If set to True then the closed loop will simulate images
            without sky background. (The default is False.)
        turnOffAtmosphere : bool, optional
            If set to True then will turn off the imsim atmosphere.
            (The default is False.)
        """
        state0 = self.ofcCalc.ofc_controller.aggregated_state
        self.imsimCmpt.dofInUm = state0

        # If using wavefront sensors we measure one per pair
        # and the field
        if camType == CamType.LsstCam:
            cornerSensorNameList = self.getSensorNameListOfFields(instName)
            cornerSensorIdList = self.getSensorIdListOfFields(instName)
            refSensorNameList = []
            refSensorIdList = []
            for name, id in zip(cornerSensorNameList, cornerSensorIdList):
                if name.endswith("SW0"):
                    refSensorNameList.append(name)
                    refSensorIdList.append(id)
        else:
            refSensorNameList = self.getSensorNameListOfFields(instName)
            refSensorIdList = self.getSensorIdListOfFields(instName)

        # Common file and directory names
        opdZkFileName = "opd.zer"
        opdPssnFileName = "PSSN.txt"
        outputDirName = "pert"
        outputImgDirName = "img"
        iterDefaultDirName = "iter"
        dofInUmFileName = "dofPertInNextIter.mat"
        fwhmItersFileName = "fwhmIters.png"
        if pipelineFile == "":
            pipelineFile = None
        if imsimConfigPointerFile == "":
            imsimConfigPointerFile = None

        # Specific file names to the amplifier/eimage
        wfsZkFileName = "wfs.zer"

        # Do the iteration
        seqNum = 1000

        for iterCount in range(iterNum):
            # Set the observation sequence number
            obsMetadata.seqNum = seqNum + iterCount * 10

            # The iteration directory
            iterDirName = "%s%d" % (iterDefaultDirName, iterCount)

            # Set the output directory
            outputDir = os.path.join(baseOutputDir, iterDirName, outputDirName)
            makeDir(outputDir)
            self.imsimCmpt.outputDir = outputDir

            # Set the output image directory
            outputImgDir = os.path.join(baseOutputDir, iterDirName, outputImgDirName)
            makeDir(outputImgDir)
            self.imsimCmpt.outputImgDir = outputImgDir

            # Generate the sky images and calculate the wavefront error
            if camType == CamType.LsstCam:
                self._generateImages(
                    obsMetadata,
                    instName=instName,
                    skySeed=skySeed,
                    pertSeed=pertSeed,
                    numPro=numPro,
                    imsimConfigPointerFile=imsimConfigPointerFile,
                    turnOffSkyBackground=turnOffSkyBackground,
                    turnOffAtmosphere=turnOffAtmosphere,
                )
            elif camType == CamType.LsstFamCam:
                for focusZ in [-1.5, 1.5]:
                    obsMetadata.seqNum += 1
                    obsMetadata.focusZ = focusZ
                    self._generateImages(
                        obsMetadata,
                        instName=instName,
                        skySeed=skySeed,
                        pertSeed=pertSeed,
                        numPro=numPro,
                        imsimConfigPointerFile=imsimConfigPointerFile,
                        turnOffSkyBackground=turnOffSkyBackground,
                        turnOffAtmosphere=turnOffAtmosphere,
                    )

            # Analyze the OPD data
            # Rotate OPD in the reversed direction of camera
            # TODO: is the above still true in imsim?
            self.imsimCmpt.analyzeOpdData(
                instName,
                zkFileName=opdZkFileName,
                rotOpdInDeg=-obsMetadata.rotatorAngle,
                pssnFileName=opdPssnFileName,
            )

            if self.useCcdImg:
                if camType in [CamType.LsstCam, CamType.LsstFamCam]:
                    listOfWfErr = self._calcWfErrFromImg(
                        obsMetadata,
                        butlerRootPath=butlerRootPath,
                        instName=instName,
                        numPro=numPro,
                        pipelineFile=pipelineFile,
                    )
            else:
                listOfWfErr = self.imsimCmpt.mapOpdDataToListOfWfErr(
                    opdZkFileName, refSensorIdList, refSensorNameList
                )

            # Get the PSSN from file
            pssn = self.imsimCmpt.getOpdPssnFromFile(opdPssnFileName)
            self.log.info("Calculated PSSN is %s." % pssn)

            # Get the GQ effective FWHM from file
            gqEffFwhm = self.imsimCmpt.getOpdGqEffFwhmFromFile(opdPssnFileName)
            self.log.info("GQ effective FWHM is %.4f." % gqEffFwhm)

            # Set the FWHM data
            fwhm, sensor_id = self.imsimCmpt.getListOfFwhmSensorData(
                opdPssnFileName, refSensorIdList
            )

            self.imsimCmpt.reorderAndSaveWfErrFile(
                listOfWfErr,
                refSensorNameList,
                getCamera(instName),
                zkFileName=wfsZkFileName,
            )

            # Calculate the DOF
            wfe = np.array(
                [sensor_wfe.annularZernikePoly for sensor_wfe in listOfWfErr]
            )

            sensor_names = np.array(
                [sensor_wfe.sensorName for sensor_wfe in listOfWfErr]
            )
            field_idx = np.array(
                [
                    self.ofcCalc.ofc_data.field_idx[sensor_name]
                    for sensor_name in sensor_names
                ]
            )

            fwhmDict = {x: y for x, y in zip(sensor_id, fwhm)}
            ofc_fwhm = np.array(
                [fwhmDict[sensor_wfe.sensorId] for sensor_wfe in listOfWfErr]
            )

            if camType == CamType.LsstCam:
                # For the wavefront sensors the sensor ids
                # are different than the corresponding field row
                # index in the sensitivity matrix.
                self.ofcCalc.set_fwhm_data(ofc_fwhm, field_idx)
            else:
                self.ofcCalc.set_fwhm_data(fwhm, sensor_id)

            # Flip zernikes that are not symmetric across the y-axis
            # based upon flip in batoid coordinate system.
            # wfe[:, [1, 4, 6, 9, 11, 12, 14, 16]] *= -1.0

            self.ofcCalc.calculate_corrections(
                wfe=wfe,
                field_idx=field_idx,
                filter_name=obsMetadata.band.upper(),
                gain=-1,
                rot=obsMetadata.rotatorAngle,
            )

            # Set the new aggregated DOF to phosimCmpt
            dofInUm = self.ofcCalc.ofc_controller.aggregated_state
            self.imsimCmpt.dofInUm = dofInUm

            # Save the DOF file
            self.imsimCmpt.saveDofInUmFileForNextIter(dofInUmFileName=dofInUmFileName)

        # Summarize the FWHM
        pssnFiles = [
            os.path.join(
                baseOutputDir,
                "%s%d" % (iterDefaultDirName, num),
                outputImgDirName,
                opdPssnFileName,
            )
            for num in range(iterNum)
        ]
        saveToFilePath = os.path.join(baseOutputDir, fwhmItersFileName)
        plotFwhmOfIters(pssnFiles, saveToFilePath=saveToFilePath)

    def _generateImages(
        self,
        obsMetadata,
        instName,
        skySeed=42,
        pertSeed=11,
        numPro=1,
        imsimConfigPointerFile=None,
        turnOffSkyBackground=False,
        turnOffAtmosphere=False,
    ):
        """Calculate the wavefront error from the images generated by PhoSim.

        Parameters
        ----------
        obsMetadata : lsst.ts.imsim.ObsMetadata object
            Observation metadata.
        instName : str
            Instrument name.
        skySeed : int, optional
            Random seed for the sky background.
            (The default is 42.)
        pertSeed : int, optional
            Random seed for the perturbations.
            (The default is 11.)
        numPro : int, optional
            Number of processor to run imSim. (the default is 1.)
        imsimConfigPointerFile : str or None, optional
            Path to imsim config pointer file.
            If None then the code will use the default in policy directory.
            (The default is None.)
        turnOffSkyBackground : bool, optional
            If set to True then the closed loop will simulate images
            without sky background. (The default is False.)
        turnOffAtmosphere : bool, optional
            If set to True then will turn off the imsim atmosphere.
            (The default is False.)

        Returns
        -------
        list[lsst.ts.wep.ctrlIntf.SensorWavefrontError]
            List of SensorWavefrontError object.
        """

        # Generate the images
        if imsimConfigPointerFile is None:
            if instName == "lsst":
                imsimConfigPointerFile = os.path.join(
                    getConfigDir(), "lsstCamDefaultPointer.yaml"
                )
        baseConfigYaml = self.imsimCmpt.assembleConfigYaml(
            obsMetadata, imsimConfigPointerFile, instName
        )

        instCat = self.imsimCmpt.genInstanceCatalog(self.skySim)
        instCatPath = os.path.join(self.imsimCmpt.outputDir, "instCat.txt")
        with open(instCatPath, "w") as file:
            file.write(instCat)

        # Override imsim config defaults with instance catalog info
        baseConfigYaml["image"].pop("image_pos")
        baseConfigYaml["output"]["nproc"] = numPro
        baseConfigYaml["image"]["random_seed"] = skySeed
        baseConfigYaml["input"]["telescope"]["fea"]["m1m3_lut"]["seed"] = pertSeed
        if turnOffSkyBackground:
            baseConfigYaml["image"]["sky_level"] = 0
        if turnOffAtmosphere:
            baseConfigYaml["input"].pop("atm_psf")
            if {"type": "AtmosphericPSF"} in baseConfigYaml["psf"]["items"]:
                baseConfigYaml["psf"]["items"].remove({"type": "AtmosphericPSF"})
                baseConfigYaml["psf"]["items"].append(
                    {"type": "Kolmogorov", "fwhm": 0.7}
                )

        if instName == "lsst":
            imsimConfigYaml = self.imsimCmpt.addSourcesToConfig(
                baseConfigYaml, instCatPath, useCcdImg=self.useCcdImg
            )
            imsimConfigPath = os.path.join(
                self.imsimCmpt.outputDir, f"imsimConfig_{obsMetadata.seqNum}.yaml"
            )
            self.log.info(f"Writing Imsim Configuration file to {imsimConfigPath}")
            self.imsimCmpt.writeYamlAndRunImsim(imsimConfigPath, imsimConfigYaml)
        elif instName == "lsstfam":
            if self.useCcdImg:
                # Run once for OPD
                imsimOpdConfigPath = os.path.join(
                    self.imsimCmpt.outputDir, "imsimConfig_opd.yaml"
                )
                if not os.path.exists(imsimOpdConfigPath):
                    imsimConfigYaml = deepcopy(baseConfigYaml)
                    imsimConfigYaml = self.imsimCmpt.addSourcesToConfig(
                        imsimConfigYaml, instCatPath, useCcdImg=False
                    )
                    self.log.info(
                        f"Writing Imsim Configuration file to {imsimOpdConfigPath}"
                    )
                    self.imsimCmpt.writeYamlAndRunImsim(
                        imsimOpdConfigPath, imsimConfigYaml
                    )

                # Run CCD images
                imsimConfigYaml = self.imsimCmpt.addSourcesToConfig(
                    baseConfigYaml, instCatPath, useCcdImg=self.useCcdImg
                )

                # Add defocus
                imsimConfigYaml["input"]["telescope"]["focusZ"] = obsMetadata.focusZ * 1e-3
                # Remove OPD since we already created it
                imsimConfigYaml["output"].pop("opd")
                imsimConfigPath = os.path.join(
                    self.imsimCmpt.outputDir, f"imsimConfig_{obsMetadata.seqNum}.yaml"
                )
                self.log.info(f"Writing Imsim Configuration file to {imsimConfigPath}")
                self.imsimCmpt.writeYamlAndRunImsim(imsimConfigPath, imsimConfigYaml)
            else:
                # Run OPD only mode
                imsimConfigYaml = self.imsimCmpt.addSourcesToConfig(
                    baseConfigYaml, instCatPath, useCcdImg=False
                )
                imsimConfigPath = os.path.join(
                    self.imsimCmpt.outputDir, "imsimConfig.yaml"
                )
                imsimOpdPath = os.path.join(
                    self.imsimCmpt.outputImgDir,
                    imsimConfigYaml["output"]["opd"]["file_name"],
                )
                if os.path.exists(imsimOpdPath):
                    self.log.info(f"OPD already created, moving to analysis.")
                else:
                    self.log.info(
                        f"Writing Imsim Configuration file to {imsimConfigPath}"
                    )
                    self.imsimCmpt.writeYamlAndRunImsim(
                        imsimConfigPath, imsimConfigYaml
                    )

    def _calcWfErrFromImg(
        self,
        obsMetadata,
        butlerRootPath,
        instName,
        numPro=1,
        pipelineFile=None,
        filterTypeName="",
    ):
        """Calculate the wavefront error from the images generated by PhoSim.

        Parameters
        ----------
        obsMetadata : lsst.ts.imsim.ObsMetadata object
            Observation metadata.
        butlerRootPath : str
            Path to the butler repository.
        instName : str
            Instrument name.
        numPro : int, optional
            Number of processor to run DM pipeline. (the default is 1.)
        pipelineFile : str or None, optional
            Path to existing pipeline yaml file to use.
            If None then the code will write its own default pipeline yaml.
            (The default is None.)
        filterTypeName : str, optional
            Filter type name: ref (or ''), u, g, r, i, z, or y.

        Returns
        -------
        list[lsst.ts.wep.ctrlIntf.SensorWavefrontError]
            List of SensorWavefrontError object.
        """

        # Ingest images into butler gen3
        self.ingestData(butlerRootPath=butlerRootPath, instName=instName)

        listOfWfErr = self.runWep(
            obsMetadata.seqNum,
            butlerRootPath,
            instName,
            numPro=numPro,
            pipelineFile=pipelineFile,
            filterTypeName=filterTypeName,
        )

        return listOfWfErr

    def runWep(
        self,
        seqNum,
        butlerRootPath,
        instName,
        numPro=1,
        pipelineFile=None,
        filterTypeName="",
    ):
        """Run wavefront estimation pipeline task for wavefront sensors.

        Parameters
        ----------
        seqNum : int
            Observation id.
        butlerRootPath : str
            Path to the butler gen3 repos.
        instName : str
            Instrument name.
        numPro : int, optional
            Number of processor to run DM pipeline. (the default is 1.)
        pipelineFile : str or None, optional
            Path to existing pipeline yaml file to use.
            If None then the code will write its own default pipeline yaml.
            (The default is None.)
        filterTypeName : str, optional
            Filter type name: ref (or ''), u, g, r, i, z, or y.

        Returns
        -------
        list[lsst.ts.wep.ctrlIntf.SensorWavefrontError]
            List of SensorWavefrontError with the results of the wavefront
            estimation pipeline for each sensor.
        """

        butlerInstName = "Cam"
        if pipelineFile is None:
            pipelineYaml = f"{instName}Pipeline.yaml"
            pipelineYamlPath = os.path.join(butlerRootPath, pipelineYaml)
            self.writeWepConfiguration(instName, pipelineYamlPath, filterTypeName)
        else:
            pipelineYamlPath = pipelineFile

        butler = dafButler.Butler(butlerRootPath)

        if f"LSST{butlerInstName}/calib" not in butler.registry.queryCollections():
            self.log.info("Ingesting curated calibrations.")

            runProgram(
                f"butler write-curated-calibrations {butlerRootPath} lsst.obs.lsst.Lsst{butlerInstName}"
            )
        # Sequence number or seqNum is an integer number
        # associated with each image taken in a single day.
        # The limit for seqNum is 5 digits,
        # set by the expectation that no more than 100K images
        # could be taken in a single day (i.e. no more than 1/sec).
        if instName == "lsst":
            runProgram(
                f"pipetask run -b {butlerRootPath} "
                f"-i refcats,LSST{butlerInstName}/raw/all,LSST{butlerInstName}/calib/unbounded "
                f"--instrument lsst.obs.lsst.Lsst{butlerInstName} "
                f"--register-dataset-types --output-run ts_imsim_{seqNum} -p {pipelineYamlPath} -d "
                f'"visit.seq_num IN ({seqNum})" -j {numPro}'
            )
        elif instName == "lsstfam":
            runProgram(
                f"pipetask run -b {butlerRootPath} "
                f"-i refcats,LSST{butlerInstName}/raw/all,LSST{butlerInstName}/calib/unbounded "
                f"--instrument lsst.obs.lsst.Lsst{butlerInstName} "
                f"--register-dataset-types --output-run ts_imsim_{seqNum} -p {pipelineYamlPath} -d "
                f'"visit.seq_num IN ({seqNum-1}, {seqNum})" -j {numPro}'
            )

        # Need to redefine butler because the database changed.
        butler = dafButler.Butler(butlerRootPath)

        datasetRefs = butler.registry.queryDatasets(
            datasetType="zernikeEstimateAvg", collections=[f"ts_imsim_{seqNum}"]
        )

        # Get the map for detector Id to detector name
        camera = butler.get(
            "camera",
            {"instrument": f"LSST{butlerInstName}"},
            collections=[f"LSST{butlerInstName}/calib/unbounded"],
        )
        detIdMap = camera.getIdMap()
        detNameMap = camera.getNameMap()

        listOfWfErr = []

        for dataset in datasetRefs:
            dataId = {
                "instrument": dataset.dataId["instrument"],
                "detector": dataset.dataId["detector"],
                "visit": dataset.dataId["visit"],
            }

            zerCoeff = butler.get(
                "zernikeEstimateAvg",
                dataId=dataId,
                collections=[f"ts_imsim_{seqNum}"],
            )

            sensorWavefrontData = SensorWavefrontError()
            # Rotate from CCS to ZCS/PCS (Phosim Coordinate System)
            rotatedName = self.rotateDetNameCcsToZcs(
                detIdMap[dataset.dataId["detector"]].getName()
            )
            sensorWavefrontData.sensorName = rotatedName
            sensorWavefrontData.sensorId = detNameMap[rotatedName].getId()
            sensorWavefrontData.annularZernikePoly = zerCoeff

            listOfWfErr.append(sensorWavefrontData)

        return listOfWfErr

    def rotateDetNameCcsToZcs(self, detName):
        """Rotate the sensor name from CCS (Camera Coordinate System)
        to the ZCS (Zemax Coordinate System) that OFC expects.

        Parameters
        ----------
        detName : str
            Detector name in CCS.

        Returns
        -------
        str
            Detector expected at that location in ZCS.
        """
        raft, sensor = detName.split('_')
        camRX = int(raft[1]) - 2
        camRY = int(raft[2]) - 2

        zcsRX = -camRX
        zcsRY = camRY

        zcsRaft = f'R{zcsRX + 2}{zcsRY + 2}'

        if sensor.startswith('SW'):
            zcsSensor = sensor
        else:
            camSX = int(sensor[1]) - 1
            camSY = int(sensor[2]) - 1
            zcsSX = -camSX
            zcsSY = camSY
            zcsSensor = f'S{zcsSX + 1}{zcsSY + 1}'

        return f'{zcsRaft}_{zcsSensor}'

    def writeWepConfiguration(self, instName, pipelineYamlPath, filterTypeName):
        """Write wavefront estimation pipeline task configuration.

        Parameters
        ----------
        instName : str
            Name of the instrument this configuration is intended for.
        pipelineYamlPath : str
            Path where the pipeline task configuration yaml file
            should be saved.
        filterTypeName : str
            Filter type name: ref (or ''), u, g, r, i, z, or y.
        """

        butlerInstName = "ComCam" if instName == "comcam" else "Cam"

        # Remap reference filter
        filterTypeName = self.mapFilterRefToG(filterTypeName)

        with open(pipelineYamlPath, "w") as fp:
            fp.write(
                f"""# This yaml file is used to define the tasks and configuration of
# a Gen 3 pipeline used for testing in ts_wep.
description: wep basic processing test pipeline
# Here we specify the corresponding instrument for the data we
# will be using.
instrument: lsst.obs.lsst.Lsst{butlerInstName}
# Use imported instrument configuration
imports:
  - location: {getWepConfigDir()}/cwfs/instData/{instName}/instParamPipeConfig.yaml
# Then we can specify each task in our pipeline by a name
# and then specify the class name corresponding to that task
tasks:
  isr:
    class: lsst.ip.isr.isrTask.IsrTask
    # Below we specify the configuration settings we want to use
    # when running the task in this pipeline. Since our data doesn't
    # include bias or flats we only want to use doApplyGains and
    # doOverscan in our isr task.
    config:
      connections.outputExposure: 'postISRCCD'
      doBias: False
      doVariance: False
      doLinearize: False
      doCrosstalk: False
      doDefect: False
      doNanMasking: False
      doInterpolate: False
      doBrighterFatter: False
      doDark: False
      doFlat: False
      doApplyGains: True
      doFringe: False
      doOverscan: True
      python: OverscanCorrectionTask.ConfigClass.fitType = 'MEDIAN'
  generateDonutCatalogWcsTask:
    class: lsst.ts.wep.task.generateDonutCatalogWcsTask.GenerateDonutCatalogWcsTask
"""
            )

    def runImg(
        self,
        inst,
        filterTypeName,
        rotCamInDeg,
        boresight,
        mjd,
        baseOutputDir,
        pathSkyFile,
        doErsDirCont,
        skySeed,
        pertSeed,
        iterNum,
        pipelineFile,
        imsimConfigPointerFile,
        turnOffSkyBackground,
        turnOffAtmosphere,
        opdOnly,
        numPro,
    ):
        """Run the simulation of images.

        Parameters
        ----------
        inst : str
            Instrument to use: currently only lsst.
        filterTypeName : str
            Filter type name: ref, u, g, r, i, z, or y.
        rotCamInDeg : float
            The camera rotation angle in degree (-90 to 90).
        boresight : list[float]
            Boresight [ra, dec] in degree.
        mjd : float
            MJD of the observation.
        baseOutputDir : str
            Base output directory.
        pathSkyFile : str
            Path to the sky file.
        doErsDirCont : bool
            Do the erase of the content of base output directory or not.
        skySeed : int
            Random seed for the sky background.
        pertSeed : int
            Random seed for the perturbations.
        iterNum : int
            Number of closed-loop iteration.
        pipelineFile : str
            Path to existing pipeline yaml file to use. If empty string
            then the code will write its own default pipeline yaml.
        imsimConfigPointerFile : str
            Path to pointer file with locations of yaml configuration
            files for imsim submodules. If empty string then the code
            will use the default in policy/config for the given inst.
        turnOffSkyBackground : bool
            If set to True then the closed loop will simulate images
            without sky background.
        turnOffAtmosphere : bool
            If set to True then will turn off the imsim atmosphere.
        opdOnly : bool
            If set to True then will run the closed loop only with
            the OPD and not the actual simulated images.
        numPro : int
            Number of processors to use.
        """
        camType, instName = self.getCamTypeAndInstName(inst)
        baseOutputDir = self.checkAndCreateBaseOutputDir(baseOutputDir)
        if doErsDirCont:
            self.eraseDirectoryContent(baseOutputDir)
        self.checkBoresight(boresight)
        self.boresightRa = boresight[0]
        self.boresightDec = boresight[1]
        self.boresightRotAng = rotCamInDeg
        # Remap the reference filter to g
        filterTypeName = self.mapFilterRefToG(filterTypeName)

        if opdOnly is True:
            self.useCcdImg = False

        obsMetadata = ObsMetadata(
            ra=self.boresightRa,
            dec=self.boresightDec,
            band=filterTypeName.lower(),
            rotatorAngle=self.boresightRotAng,
            mjd=mjd,
        )

        # Configure the components
        self.opdMetr = OpdMetrology()
        self.opdMetr.setCamera(instName)
        self.configSkySim(instName, obsMetadata, pathSkyFile=pathSkyFile, starMag=15)
        self.configOfcCalc(instName)
        self.imsimCmpt = ImsimCmpt()

        # If pathSkyFile using default OPD positions write this to disk
        # so that the Butler can load it later
        if pathSkyFile == "":
            pathSkyFile = os.path.join(baseOutputDir, "sky_info.txt")
            self.skySim.exportSkyToFile(pathSkyFile)
            self.log.info(f"Wrote new sky file to {pathSkyFile}.")

        # generate butler gen3 repo if needed
        butlerRootPath = os.path.join(baseOutputDir, "imsimData")

        if self.useCcdImg:
            self.generateButler(butlerRootPath, instName)
            self.generateRefCatalog(
                instName=instName,
                butlerRootPath=butlerRootPath,
                pathSkyFile=pathSkyFile,
                filterTypeName=filterTypeName,
            )

        if instName == "lsst":
            # Append equal weights for CWFS fields to OFC data
            # Assign equal normalized weights to each of the
            # four corner wavefront sensor pairs.
            self.ofcCalc.ofc_data.normalized_image_quality_weight = np.append(
                self.ofcCalc.ofc_data.normalized_image_quality_weight,
                [0.25, 0.25, 0.25, 0.25],
            )

        self._runSim(
            camType,
            instName,
            obsMetadata,
            baseOutputDir,
            butlerRootPath,
            skySeed,
            pertSeed,
            iterNum,
            numPro=numPro,
            pipelineFile=pipelineFile,
            imsimConfigPointerFile=imsimConfigPointerFile,
            turnOffSkyBackground=turnOffSkyBackground,
            turnOffAtmosphere=turnOffAtmosphere,
        )

    def generateButler(self, butlerRootPath, instName):
        """Generate butler gen3.

        Parameters
        ----------
        butlerRootPath: `str`
            Path to where the butler repository should be created.
        instName: `str`
            Name of the instrument.
        """

        self.log.info(f"Generating butler gen3 in {butlerRootPath} for {instName}")

        runProgram(f"butler create {butlerRootPath}")

        if instName == "comcam":
            self.log.debug("Registering LsstComCam")
            runProgram(
                f"butler register-instrument {butlerRootPath} lsst.obs.lsst.LsstComCam"
            )
        else:
            self.log.debug("Registering LsstCam")
            runProgram(
                f"butler register-instrument {butlerRootPath} lsst.obs.lsst.LsstCam"
            )

    def generateRefCatalog(self, instName, butlerRootPath, pathSkyFile, filterTypeName):
        """Generate reference star catalog.

        Parameters
        ----------
        instName: `str`
            Name of the instrument.
        butlerRootPath: `str`
            Path to the butler gen3 repository.
        pathSkyFile: `str`
            Path to the catalog star file.
        filterTypeName : str
            Filter type name: ref (or ''), u, g, r, i, z, or y.
        """
        self.log.debug("Creating reference catalog.")

        catDir = os.path.join(butlerRootPath, "skydata")
        skyFilename = os.path.join(catDir, "sky_data.csv")
        catConfigFilename = os.path.join(catDir, "cat.cfg")
        skyEcsvFilename = os.path.join(catDir, "filename_to_htm.ecsv")
        catLogFilename = os.path.join(catDir, "convert.log")
        os.mkdir(catDir)

        # Read sky file and convert it to csv
        skyData = astropy.io.ascii.read(pathSkyFile)

        # Constructing the catalog of stars to use in the wavefront estimation
        # pipeline. It is used for target
        # selection, and affects magnitude limits
        # as set in generateDonutCatalogWcsTask pipeline yaml file
        skyData.rename_column("Mag", filterTypeName)

        skyData.write(skyFilename, format="csv", overwrite=True)

        with open(catConfigFilename, "w") as fp:
            fp.write(
                f"""config.ra_name='Ra'
config.dec_name='Dec'
config.id_name='Id'
config.mag_column_list=['{filterTypeName}']
config.dataset_config.ref_dataset_name='ref_cat'
"""
            )

        runProgram(
            f"convertReferenceCatalog {catDir} {catConfigFilename} {skyFilename} &> {catLogFilename}"
        )

        runProgram(
            f"butler register-dataset-type {butlerRootPath} cal_ref_cat SimpleCatalog htm7"
        )

        runProgram(
            f"butler ingest-files -t direct {butlerRootPath} cal_ref_cat refcats {skyEcsvFilename}"
        )

    def ingestData(self, butlerRootPath, instName):
        """Ingest data into a gen3 data Butler.

        Parameters
        ----------
        butlerRootPath : str
            Path to the butler repository.
        instName : str
            Instrument name.
        """
        outputImgDir = self.imsimCmpt.outputImgDir

        if instName in ["lsst", "lsstfam"]:
            runProgram(f"butler ingest-raws {butlerRootPath} {outputImgDir}/amp*")

        runProgram(f"butler define-visits {butlerRootPath} lsst.obs.lsst.LsstCam")

    def eraseDirectoryContent(self, targetDir):
        """Erase the directory content.

        Parameters
        ----------
        targetDir : str
            Target directory.
        """

        for theFile in os.listdir(targetDir):
            filePath = os.path.join(targetDir, theFile)
            if os.path.isfile(filePath):
                os.unlink(filePath)
            elif os.path.isdir(filePath):
                shutil.rmtree(filePath)

    @staticmethod
    def setDefaultParser(parser):
        """Set the default parser.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            Input parser.

        Returns
        -------
        argparse.ArgumentParser
            Updated parser.
        """

        parser.add_argument(
            "--inst",
            type=str,
            default="lsst",
            help="Instrument to use: currently only lsst. (default: lsst)",
        )

        parser.add_argument(
            "--filterType",
            type=str,
            default="",
            help="Filter type to use: U, G, R, I, Z, Y or empty string for "
            "reference wavelength. (default: '')",
        )

        parser.add_argument(
            "--rotCam",
            type=float,
            default=0.0,
            help="Rotate camera (degree) in counter-clockwise direction. (default: 0.0)",
        )

        parser.add_argument("--output", type=str, default="", help="Output directory.")

        parser.add_argument(
            "--log-level", type=int, default=logging.INFO, help="Log level."
        )

        parser.add_argument(
            "--clobber",
            default=False,
            action="store_true",
            help="Delete existing output directory.",
        )

        parser.add_argument(
            "--configPointerFile",
            type=str,
            default="",
            help="Imsim Configuration Pointer File.",
        )

        parser.add_argument(
            "--skySeed",
            type=int,
            default=42,
            help="Random seed for imsim sky (default: 42).",
        )

        parser.add_argument(
            "--pertSeed",
            type=int,
            default=11,
            help="Random seed for m1m3_lut fractional actuator random error. "
            "(Default: 11)",
        )

        parser.add_argument(
            "--iterNum",
            type=int,
            default=5,
            help="Number of closed-loop iterations. (default: 5)",
        )

        parser.add_argument(
            "--pipelineFile",
            type=str,
            default="",
            help="""
            Location of user-specified pipeline configuration file.
            If left as empty string the code will create a default file.
            (default: '')
            """,
        )

        parser.add_argument(
            "--numProc",
            type=int,
            default=1,
            help="Number of processor to run imSim and DM pipeline. (default: 1)",
        )

        return parser

    @staticmethod
    def setImgParser(parser):
        """Set the image-specific parser.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            Input parser.

        Returns
        -------
        argparse.ArgumentParser
            Updated parser.
        """

        parser.add_argument(
            "--boresightDeg",
            type=float,
            nargs=2,
            default=[0, 0],
            help="Boresight [ra, dec] in degree. The default is [0, 0].",
        )

        parser.add_argument(
            "--skyFile",
            type=str,
            default="",
            help="""
                 Text file contains the star Id, ra, dec, and magnitude.
                 The default is to use the OPD field positions with boresight
                 [ra, dec] = [0, 0].
                 """,
        )

        parser.add_argument(
            "--mjd", type=float, default=59580, help="Starting MJD of observation."
        )

        parser.add_argument(
            "--turnOffSkyBackground",
            action="store_true",
            help="Turn sky brightness model off.",
        )

        parser.add_argument(
            "--turnOffAtmosphere", action="store_true", help="Turn atmosphere off."
        )

        parser.add_argument(
            "--opdOnly", action="store_true", help="Turn atmosphere off."
        )

        return parser
