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
from lsst.ts.ofc.utils import get_config_dir as getConfigDirOfc
from lsst.ts.imsim.utils.utility import getCamera, getPolicyPath


class OpdMetrology:
    def __init__(self):
        """Initialization of OPD metrology class.

        OPD: Optical path difference.
        """

        self.wt = np.array([])
        self.fieldX = np.array([])
        self.fieldY = np.array([])
        self.sensorIds = []
        self._camera = None

    def setCamera(self, instName):
        """Set the camera.

        Parameters
        ----------
        instName : `str`
            Instrument name. Valid options are 'comcam or 'lsstfam'.
        """

        self._camera = getCamera(instName)

    def setDefaultLsstWfsGQ(self):
        """Set default values for LSST WFS field X, Y
        and weighting ratio.
        """

        # Set equal full weights for each of the
        # four corner wavefront sensor pairs.
        self.wt = np.array([0.25, 0.25, 0.25, 0.25])
        wfsFieldX, wfsFieldY, sensorIds = self.getDefaultLsstWfsGQ()
        self.fieldX, self.fieldY = (wfsFieldX, wfsFieldY)
        self.sensorIds = sensorIds

    def getDefaultLsstWfsGQ(self):
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

        # Field x, y for 4 WFS
        fieldWFSx = [1.176, -1.176, -1.176, 1.176]
        fieldWFSy = [1.176, 1.176, -1.176, -1.176]
        detIds = [203, 199, 191, 195]

        return fieldWFSx, fieldWFSy, detIds

    def setWgtAndFieldXyOfGQ(self, instName):
        """Set the GQ weighting ratio and field X, Y.

        GQ: Gaussian quadrature.

        Parameters
        ----------
        instName : `str`
            Instrument name.

        Raises
        ------
        RuntimeError
            If the instrument path does not exists.
            If fieldXy.yaml file does not exists in the instrument
            configuration directory.
        """

        instrumentPath = getConfigDirOfc() / instName

        if not instrumentPath.exists():
            raise RuntimeError(f"OFC instrument path does not exist: {instrumentPath}")

        # Set the weighting ratio
        pathWgtFile = instrumentPath / "imgQualWgt.yaml"
        with open(pathWgtFile, 'r') as file:
            wgt = yaml.safe_load(file)
        wgtValues = np.array(list(wgt.values()), dtype=float)
        # Normalize weights
        self.wt = wgtValues / np.sum(wgtValues)

        # Set the field (x, y)
        pathFieldXyFile = os.path.join(
            getPolicyPath(), "instrument", instName, "fieldXy.yaml"
        )

        if not os.path.exists(pathFieldXyFile):
            raise RuntimeError(f"Field xy file does not exists: {pathFieldXyFile}.")

        with open(pathFieldXyFile, 'r') as file:
            fieldXY = yaml.safe_load(file)
        fieldXY = np.array(fieldXY, dtype=float)
        self.fieldX, self.fieldY = (fieldXY[:, 0], fieldXY[:, 1])
