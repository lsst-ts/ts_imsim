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

from lsst.ts.imsim.closed_loop_task import ClosedLoopTask

if __name__ == "__main__":
    # Set the parser
    parser = argparse.ArgumentParser(
        description="Run AOS closed-loop simulation (default is amplifier files)."
    )
    parser = ClosedLoopTask.set_default_parser(parser)
    parser = ClosedLoopTask.set_img_parser(parser)

    # Get the arguments
    args = parser.parse_args()

    logging.basicConfig(format="%(levelname)s:%(message)s", level=args.log_level)

    # Run the simulation
    closeLoopTask = ClosedLoopTask()
    closeLoopTask.run_img(
        args.inst,
        args.filter_type,
        args.rot_cam,
        args.boresight_deg,
        args.mjd,
        args.star_mag,
        args.output,
        args.sky_file,
        args.clobber,
        args.sky_seed,
        args.pert_seed,
        args.iter_num,
        args.pipeline_file,
        args.config_pointer_file,
        args.turn_off_sky_background,
        args.turn_off_atmosphere,
        args.turn_off_wavefront_estimates,
        args.num_proc,
        args.raw_seeing,
        args.imsim_log_file,
        args.wep_estimator,
    )
