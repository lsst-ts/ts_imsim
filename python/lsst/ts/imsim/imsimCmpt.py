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

from lsst.ts.ofc.ofc_data.base_ofc_data import BaseOFCData
from lsst.ts.imsim.opdMetrology import OpdMetrology
from lsst.ts.imsim.utils.utility import getConfigDir


class ImsimCmpt:
    def __init__(self):
        """Class to take configurations for each imsim component and
        generate full imsim configuration files.
        """
        # Output directories
        self.outputDir = None
        self.outputImgDir = None

        # OPD information
        self.opdFilePath = None
        self.opdMetr = OpdMetrology()

        # Specify number of Zernikes
        self.numOfZk = 19

        # AOS Degrees of Freedom
        self.numOfDof = 50
        self.dofInUm = np.zeros(self.numOfDof, dtype=float)

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
        for subsection in ["atm_psf", "sky_model", "telescope"]:
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
