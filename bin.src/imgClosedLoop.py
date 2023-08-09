#!/usr/bin/env python

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

import argparse
import logging

from lsst.ts.imsim.closedLoopTask import ClosedLoopTask

if __name__ == "__main__":
    # Set the parser
    parser = argparse.ArgumentParser(
        description="Run AOS closed-loop simulation (default is amplifier files)."
    )
    parser = ClosedLoopTask.setDefaultParser(parser)
    parser = ClosedLoopTask.setImgParser(parser)

    # Get the arguments
    args = parser.parse_args()

    logging.basicConfig(format="%(levelname)s:%(message)s", level=args.log_level)

    # Run the simulation
    closeLoopTask = ClosedLoopTask()
    closeLoopTask.runImg(
        args.inst,
        args.filterType,
        args.rotCam,
        args.boresightDeg,
        args.mjd,
        args.output,
        args.skyFile,
        args.clobber,
        args.skySeed,
        args.pertSeed,
        args.iterNum,
        args.pipelineFile,
        args.configPointerFile,
        args.turnOffSkyBackground,
        args.turnOffAtmosphere,
        args.opdOnly,
    )
