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
from argparse import ArgumentParser
from copy import deepcopy
from glob import glob

import astropy
import numpy as np
from lsst.afw.cameraGeom import FIELD_ANGLE, DetectorType
from lsst.daf import butler as dafButler
from lsst.ts.imsim.imsimCmpt import ImsimCmpt
from lsst.ts.imsim.obsMetadata import ObsMetadata
from lsst.ts.imsim.opdMetrology import OpdMetrology
from lsst.ts.imsim.skySim import SkySim
from lsst.ts.imsim.utils import (
    SensorWavefrontError,
    get_camera,
    get_config_dir,
    make_dir,
    plot_fwhm_of_iters,
)
from lsst.ts.ofc import OFC, OFCData
from lsst.ts.wep.utils import CamType, FilterType
from lsst.ts.wep.utils import getConfigDir as getWepConfigDir
from lsst.ts.wep.utils import rotMatrix, runProgram


class ClosedLoopTask:
    def __init__(self) -> None:
        """Initilization of the closed loop task class to
        run the simulation with imSim."""

        self.log = logging.getLogger(type(self).__name__)

        # Sky simulator
        self.sky_sim = None

        # OFC calculator
        self.ofc_calc = None

        # imSim Component
        self.imsim_cmpt = None

        # Ra/Dec/RotAng coordinates used in the simulation.
        self.boresight_ra = None
        self.boresight_dec = None
        self.boresight_rot_ang = None

        # Use CCD image
        self.use_ccd_img = True

    def config_sky_sim(
        self,
        inst_name: str,
        obs_metadata: ObsMetadata,
        path_sky_file: str = "",
        star_mag: float = 15.0,
    ) -> None:
        """Configure the sky simulator.

        If the path of sky file is not provided, The defult OPD field positions
        will be used.

        OPD: Optical path difference.

        Parameters
        ----------
        inst_name : str
            Instrument name.
        obs_metadata : lsst.ts.imsim.ObsMetadata
            Observation metadata.
        path_sky_file : str, optional
            Path to the sky file. (the default is "".)
        star_mag : float, optional
            Default star magnitude if there is no sky file used. This is to
            pretend there are the stars at OPD field positions. (the default is
            15.)

        Raises
        ------
        ValueError
            This instrument name is not supported.
        """

        self.sky_sim = SkySim()
        self.sky_sim.set_camera(inst_name)
        if path_sky_file == "":
            par_angle = self.sky_sim.calc_parallactic_angle(obs_metadata)
            self._set_sky_sim_based_on_opd_field_pos(
                inst_name, obs_metadata, par_angle, star_mag
            )
        else:
            abs_sky_file_path = os.path.abspath(path_sky_file)
            self.sky_sim.add_star_by_file(abs_sky_file_path)

    def _set_sky_sim_based_on_opd_field_pos(
        self,
        inst_name: str,
        obs_metadata: ObsMetadata,
        par_angle: float,
        star_mag: float,
    ) -> None:
        """Set the sky simulator based on the OPD field positions.

        OPD: Optical path difference.

        Parameters
        ----------
        inst_name : str
            Instrument name.
        obs_metadata : lsst.ts.imsim.ObsMetadata object
            Observation metadata.
        par_angle : float
            Parallactic angle.
        star_mag : float
            Star magnitude. This is to pretend there are the stars at OPD field
            positions.

        Raises
        ------
        ValueError
            This instrument name is not supported.
        """

        self.log.info(
            "Use the default OPD field positions to be star positions. "
            f"The star magnitude is chosen to be {star_mag}."
        )

        opd_metr = OpdMetrology()
        if inst_name in ["lsst", "lsstfam"]:
            field_x, field_y = list(), list()
            camera = get_camera(inst_name)
            for name in self.get_sensor_name_list_of_fields(inst_name):
                detector = camera.get(name)
                x_rad, y_rad = detector.getCenter(FIELD_ANGLE)
                x_deg, y_deg = np.rad2deg(x_rad), np.rad2deg(y_rad)
                field_y.append(x_deg)  # transpose to convert from DVCS to CCS
                field_x.append(y_deg)
            opd_metr.field_x = field_x
            opd_metr.field_y = field_y
        else:
            raise ValueError(f"This instrument name ({inst_name}) is not supported.")

        star_id = 0
        ra_in_deg_arr = np.array(opd_metr.field_x)
        dec_in_deg_arr = np.array(opd_metr.field_y)
        # Extra 180 degree rotation based upon this note:
        # https://lsstc.slack.com/archives/CHXKSF3HC/p1651863987821319?thread_ts=1651863934.274719&cid=CHXKSF3HC
        # that shows photons farthest from Zenith on sky appear on "top"
        # of focal plane.
        rotation = rotMatrix(obs_metadata.rotator_angle - par_angle + 180)
        for ra_in_deg, dec_in_deg in zip(ra_in_deg_arr, dec_in_deg_arr):
            # It is noted that the field position might be < 0. But it is
            # not the same case for ra (0 <= ra <= 360).
            ra_in_deg, dec_in_deg = np.dot(rotation, np.array([ra_in_deg, dec_in_deg]))
            ra_in_deg += obs_metadata.ra
            dec_in_deg += obs_metadata.dec
            if ra_in_deg < 0:
                ra_in_deg += 360.0
            self.sky_sim.add_star_by_ra_dec_in_deg(
                star_id, ra_in_deg, dec_in_deg, star_mag
            )
            star_id += 1

    def config_ofc_calc(self, inst_name: str) -> None:
        """Configure the OFC calculator.

        OFC: Optical feedback calculator.

        Parameters
        ----------
        inst_name : str
            Instrument name.
        """

        self.ofc_calc = OFC(OFCData(inst_name))

    def map_filter_ref_to_g(self, filter_type_name: str) -> str:
        """Map the reference filter to the G filter.

        Parameters
        ----------
        filter_type_name : str
            Filter type name: ref (or ''), u, g, r, i, z, or y.

        Returns
        -------
        filter_type_name : str
            Mapped filter type.
        """
        return "g" if filter_type_name in ("ref", "") else filter_type_name

    def check_boresight(self, boresight: list[float]) -> None:
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

    def get_sensor_name_list_of_fields(self, inst_name: str) -> list[str]:
        """Get the list of sensor name of fields.

        The list will be sorted based on the field index.

        Parameters
        ----------
        inst_name : enum 'InstName' in lsst.ts.ofc.Utility
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

        camera = get_camera(inst_name)
        detector_type = (
            DetectorType.WAVEFRONT if inst_name == "lsst" else DetectorType.SCIENCE
        )
        return [
            detector.getName()
            for detector in camera
            if detector.getType() == detector_type
        ]

    def get_sensor_id_list_of_fields(self, inst_name: str) -> list[int]:
        """Get the list of sensor ids of fields.

        The list will be sorted based on the field index.

        Parameters
        ----------
        inst_name : enum 'InstName' in lsst.ts.ofc.Utility
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

        camera = get_camera(inst_name)

        detector_type = (
            DetectorType.WAVEFRONT if inst_name == "lsst" else DetectorType.SCIENCE
        )
        return [
            detector.getId()
            for detector in camera
            if detector.getType() == detector_type
        ]

    def check_and_create_base_output_dir(self, base_output_dir: str) -> str:
        """Check and create the base output directory.

        This function will create the directory if it does not exist.

        Parameters
        ----------
        base_output_dir : str
            Base output directory.

        Returns
        -------
        str
            Base output directory.
        """
        output_dir = base_output_dir
        make_dir(output_dir, exist_ok=True)

        return output_dir

    def get_cam_type_and_inst_name(self, inst: str) -> tuple[CamType, str]:
        """Get the camera type and instrument name.

        Parameters
        ----------
        inst : str
            Instrument to use: currently only lsst.

        Returns
        -------
        lsst.ts.wep.utils.enumUtils.CamType
            Camera type.
        str
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

    def get_filter_type(self, filter_type_name: str) -> FilterType:
        """Get the filter type.

        Parameters
        ----------
        filter_type_name : str
            Filter type name: ref, u, g, r, i, z, or y.

        Returns
        -------
        lsst.ts.wep.utils.FilterType
            Filter type.

        Raises
        ------
        ValueError
            This filter type is not supported.
        """

        if filter_type_name in {"", "ref"}:
            return FilterType.REF
        elif filter_type_name == "u":
            return FilterType.LSST_U
        elif filter_type_name == "g":
            return FilterType.LSST_G
        elif filter_type_name == "r":
            return FilterType.LSST_R
        elif filter_type_name == "i":
            return FilterType.LSST_I
        elif filter_type_name == "z":
            return FilterType.LSST_Z
        elif filter_type_name == "y":
            return FilterType.LSST_Y
        else:
            raise ValueError(f"This filter type ({filter_type_name}) is not supported.")

    def _run_sim(
        self,
        cam_type: CamType,
        inst_name: str,
        obs_metadata: ObsMetadata,
        base_output_dir: str,
        butler_root_path: str,
        sky_seed: int,
        pert_seed: int,
        iter_num: int,
        num_pro: int = 1,
        pipeline_file: str = "",
        imsim_config_pointer_file: str = "",
        turn_off_sky_background: bool = False,
        turn_off_atmosphere: bool = False,
    ) -> None:
        """Run the simulation.

        Parameters
        ----------
        cam_type : enum 'CamType' in lsst.ts.wep.utility
            Camera type.
        inst_name : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.
        obs_metadata : lsst.ts.imsim.ObsMetadata object
            Observation metadata.
        base_output_dir : str
            Base output directory.
        butler_root_path : str
            Path to the butler gen 3 repository.
        sky_seed : int
            Random seed for the sky background.
        pert_seed : int
            Random seed for the perturbations.
        iter_num : int
            Number of closed-loop iteration.
        num_pro : int, optional
            Number of processors to use. (The default is 1.)
        pipeline_file : str, optional
            Path to existing pipeline yaml file to use. If empty string
            then the code will write its own default pipeline yaml.
            (The default is "".)
        imsim_config_pointer_file : str, optional
            Path to pointer file with locations of yaml configuration
            files for imsim submodules. If empty string then the code
            will use the default in policy/config for the given inst.
            (The default is "".)
        turn_off_sky_background : bool, optional
            If set to True then the closed loop will simulate images
            without sky background. (The default is False.)
        turn_off_atmosphere : bool, optional
            If set to True then will turn off the imsim atmosphere.
            (The default is False.)
        """
        state_0 = self.ofc_calc.ofc_controller.aggregated_state
        self.imsim_cmpt.dof_in_um = state_0

        # If using wavefront sensors we measure one per pair
        # and the field
        if cam_type == CamType.LsstCam:
            corner_sensor_name_list = self.get_sensor_name_list_of_fields(inst_name)
            corner_sensor_id_list = self.get_sensor_id_list_of_fields(inst_name)
            ref_sensor_name_list = []
            ref_sensor_id_list = []
            for sens_name, sens_id in zip(
                corner_sensor_name_list, corner_sensor_id_list
            ):
                if sens_name.endswith("SW0"):
                    ref_sensor_name_list.append(sens_name)
                    ref_sensor_id_list.append(sens_id)
        else:
            ref_sensor_name_list = self.get_sensor_name_list_of_fields(inst_name)
            ref_sensor_id_list = self.get_sensor_id_list_of_fields(inst_name)

        # Common file and directory names
        opd_zk_file_name = "opd.zer"
        opd_pssn_file_name = "PSSN.txt"
        output_dir_name = "pert"
        output_img_dir_name = "img"
        iter_default_dir_name = "iter"
        dof_in_um_file_name = "dofPertInNextIter.mat"
        fwhm_iters_file_name = "fwhmIters.png"
        if pipeline_file == "":
            pipeline_file = None
        if imsim_config_pointer_file == "":
            imsim_config_pointer_file = None

        # Specific file names to the amplifier/eimage
        wfs_zk_file_name = "wfs.zer"

        # Do the iteration
        seq_num = 1000

        for iter_count in range(iter_num):
            # Set the observation sequence number
            obs_metadata.seq_num = seq_num + iter_count * 10

            # The iteration directory
            iter_dir_name = "%s%d" % (iter_default_dir_name, iter_count)

            # Set the output directory
            output_dir = os.path.join(base_output_dir, iter_dir_name, output_dir_name)
            make_dir(output_dir)
            self.imsim_cmpt.output_dir = output_dir

            # Set the output image directory
            output_img_dir = os.path.join(
                base_output_dir, iter_dir_name, output_img_dir_name
            )
            make_dir(output_img_dir)
            self.imsim_cmpt.output_img_dir = output_img_dir

            # Generate the sky images and calculate the wavefront error
            if cam_type == CamType.LsstCam:
                self._generate_images(
                    obs_metadata,
                    inst_name=inst_name,
                    sky_seed=sky_seed,
                    pert_seed=pert_seed,
                    num_pro=num_pro,
                    imsim_config_pointer_file=imsim_config_pointer_file,
                    turn_off_sky_background=turn_off_sky_background,
                    turn_off_atmosphere=turn_off_atmosphere,
                )
            elif cam_type == CamType.LsstFamCam:
                for focus_z in [-1.5, 1.5]:
                    obs_metadata.seq_num += 1
                    obs_metadata.focus_z = focus_z
                    self._generate_images(
                        obs_metadata,
                        inst_name=inst_name,
                        sky_seed=sky_seed,
                        pert_seed=pert_seed,
                        num_pro=num_pro,
                        imsim_config_pointer_file=imsim_config_pointer_file,
                        turn_off_sky_background=turn_off_sky_background,
                        turn_off_atmosphere=turn_off_atmosphere,
                    )

            # Analyze the OPD data
            self.imsim_cmpt.analyze_opd_data(
                inst_name,
                zk_file_name=opd_zk_file_name,
                rot_opd_in_deg=obs_metadata.rotator_angle,
                pssn_file_name=opd_pssn_file_name,
            )

            if self.use_ccd_img:
                if cam_type in [CamType.LsstCam, CamType.LsstFamCam]:
                    list_of_wf_err = self._calc_wf_err_from_img(
                        obs_metadata,
                        butler_root_path=butler_root_path,
                        inst_name=inst_name,
                        num_pro=num_pro,
                        pipeline_file=pipeline_file,
                    )
            else:
                list_of_wf_err = self.imsim_cmpt.map_opd_data_to_list_of_wf_err(
                    opd_zk_file_name, ref_sensor_id_list, ref_sensor_name_list
                )

            # Get the PSSN from file
            pssn = self.imsim_cmpt.get_opd_pssn_from_file(opd_pssn_file_name)
            self.log.info("Calculated PSSN is %s." % pssn)

            # Get the GQ effective FWHM from file
            gq_eff_fwhm = self.imsim_cmpt.get_opd_gq_eff_fwhm_from_file(
                opd_pssn_file_name
            )
            self.log.info("GQ effective FWHM is %.4f." % gq_eff_fwhm)

            # Set the FWHM data
            fwhm, sensor_id = self.imsim_cmpt.get_list_of_fwhm_sensor_data(
                opd_pssn_file_name, ref_sensor_id_list
            )

            self.imsim_cmpt.reorder_and_save_wf_err_file(
                list_of_wf_err,
                ref_sensor_name_list,
                get_camera(inst_name),
                zk_file_name=wfs_zk_file_name,
            )

            # Calculate the DOF
            wfe = np.array(
                [sensor_wfe.annular_zernike_poly for sensor_wfe in list_of_wf_err]
            )

            sensor_names = np.array(
                [sensor_wfe.sensor_name for sensor_wfe in list_of_wf_err]
            )
            field_idx = np.array(
                [
                    self.ofc_calc.ofc_data.field_idx[sensor_name]
                    for sensor_name in sensor_names
                ]
            )

            fwhm_dict = {x: y for x, y in zip(sensor_id, fwhm)}
            ofc_fwhm = np.array(
                [fwhm_dict[sensor_wfe.sensor_id] for sensor_wfe in list_of_wf_err]
            )

            if cam_type == CamType.LsstCam:
                # For the wavefront sensors the sensor ids
                # are different than the corresponding field row
                # index in the sensitivity matrix.
                self.ofc_calc.set_fwhm_data(ofc_fwhm, field_idx)
            else:
                self.ofc_calc.set_fwhm_data(fwhm, sensor_id)

            self.ofc_calc.calculate_corrections(
                wfe=wfe,
                field_idx=field_idx,
                filter_name=obs_metadata.band.upper(),
                gain=-1,
                rotation_angle=-obs_metadata.rotator_angle,
            )

            # Set the new aggregated DOF to phosimCmpt
            dof_in_um = self.ofc_calc.ofc_controller.aggregated_state
            self.imsim_cmpt.dof_in_um = dof_in_um

            # Save the DOF file
            self.imsim_cmpt.save_dof_in_um_file_for_next_iter(
                dof_in_um_file_name=dof_in_um_file_name
            )

        # Summarize the FWHM
        pssn_files = [
            os.path.join(
                base_output_dir,
                "%s%d" % (iter_default_dir_name, num),
                output_img_dir_name,
                opd_pssn_file_name,
            )
            for num in range(iter_num)
        ]
        save_to_file_path = os.path.join(base_output_dir, fwhm_iters_file_name)
        plot_fwhm_of_iters(pssn_files, save_to_file_path=save_to_file_path)

    def _generate_images(
        self,
        obs_metadata: ObsMetadata,
        inst_name: str,
        sky_seed: int = 42,
        pert_seed: int = 11,
        num_pro: int = 1,
        imsim_config_pointer_file: str | None = None,
        turn_off_sky_background: bool = False,
        turn_off_atmosphere: bool = False,
    ) -> None:
        """Calculate the wavefront error from the images generated by PhoSim.

        Parameters
        ----------
        obs_metadata : lsst.ts.imsim.ObsMetadata object
            Observation metadata.
        inst_name : str
            Instrument name.
        sky_seed : int, optional
            Random seed for the sky background.
            (The default is 42.)
        pert_seed : int, optional
            Random seed for the perturbations.
            (The default is 11.)
        num_pro : int, optional
            Number of processor to run imSim. (the default is 1.)
        imsim_config_pointer_file : str or None, optional
            Path to imsim config pointer file.
            If None then the code will use the default in policy directory.
            (The default is None.)
        turn_off_sky_background : bool, optional
            If set to True then the closed loop will simulate images
            without sky background. (The default is False.)
        turn_off_atmosphere : bool, optional
            If set to True then will turn off the imsim atmosphere.
            (The default is False.)
        """

        # Generate the images
        if imsim_config_pointer_file is None:
            if inst_name == "lsst":
                imsim_config_pointer_file = os.path.join(
                    get_config_dir(), "lsstCamDefaultPointer.yaml"
                )
        base_config_yaml = self.imsim_cmpt.assemble_config_yaml(
            obs_metadata, imsim_config_pointer_file, inst_name
        )

        inst_cat = self.imsim_cmpt.gen_instance_catalog(self.sky_sim)
        inst_cat_path = os.path.join(self.imsim_cmpt.output_dir, "instCat.txt")
        with open(inst_cat_path, "w") as file:
            file.write(inst_cat)

        # Override imsim config defaults with instance catalog info
        base_config_yaml["image"].pop("image_pos")
        base_config_yaml["output"]["nproc"] = num_pro
        base_config_yaml["image"]["random_seed"] = sky_seed
        base_config_yaml["input"]["telescope"]["fea"]["m1m3_lut"]["seed"] = pert_seed
        if turn_off_sky_background:
            base_config_yaml["image"]["sky_level"] = 0
        if turn_off_atmosphere:
            base_config_yaml["input"].pop("atm_psf")
            if {"type": "AtmosphericPSF"} in base_config_yaml["psf"]["items"]:
                base_config_yaml["psf"]["items"].remove({"type": "AtmosphericPSF"})
                base_config_yaml["psf"]["items"].append(
                    {"type": "Kolmogorov", "fwhm": 0.7}
                )

        if inst_name == "lsst":
            imsim_config_yaml = self.imsim_cmpt.add_sources_to_config(
                base_config_yaml, inst_cat_path, use_ccd_img=self.use_ccd_img
            )
            imsim_config_path = os.path.join(
                self.imsim_cmpt.output_dir, f"imsimConfig_{obs_metadata.seq_num}.yaml"
            )
            self.log.info(f"Writing Imsim Configuration file to {imsim_config_path}")
            self.imsim_cmpt.write_yaml_and_run_imsim(
                imsim_config_path, imsim_config_yaml
            )
        elif inst_name == "lsstfam":
            if self.use_ccd_img:
                # Run once for OPD
                imsim_opd_config_path = os.path.join(
                    self.imsim_cmpt.output_dir, "imsimConfig_opd.yaml"
                )
                if not os.path.exists(imsim_opd_config_path):
                    imsim_config_yaml = deepcopy(base_config_yaml)
                    imsim_config_yaml = self.imsim_cmpt.add_sources_to_config(
                        imsim_config_yaml, inst_cat_path, use_ccd_img=False
                    )
                    self.log.info(
                        f"Writing Imsim Configuration file to {imsim_opd_config_path}"
                    )
                    self.imsim_cmpt.write_yaml_and_run_imsim(
                        imsim_opd_config_path, imsim_config_yaml
                    )

                # Run CCD images
                imsim_config_yaml = self.imsim_cmpt.add_sources_to_config(
                    base_config_yaml, inst_cat_path, use_ccd_img=self.use_ccd_img
                )

                # Add defocus
                imsim_config_yaml["input"]["telescope"]["focusZ"] = (
                    obs_metadata.focus_z * 1e-3
                )
                # Remove OPD from config since we already created it
                imsim_config_yaml["output"].pop("opd")
                imsim_config_path = os.path.join(
                    self.imsim_cmpt.output_dir,
                    f"imsimConfig_{obs_metadata.seq_num}.yaml",
                )
                self.log.info(
                    f"Writing Imsim Configuration file to {imsim_config_path}"
                )
                self.imsim_cmpt.write_yaml_and_run_imsim(
                    imsim_config_path, imsim_config_yaml
                )
            else:
                # Run OPD only mode
                imsim_config_yaml = self.imsim_cmpt.add_sources_to_config(
                    base_config_yaml, inst_cat_path, use_ccd_img=False
                )
                imsim_config_path = os.path.join(
                    self.imsim_cmpt.output_dir, "imsimConfig.yaml"
                )
                imsimOpdPath = os.path.join(
                    self.imsim_cmpt.output_img_dir,
                    imsim_config_yaml["output"]["opd"]["file_name"],
                )
                if os.path.exists(imsimOpdPath):
                    self.log.info(f"OPD already created, moving to analysis.")
                else:
                    self.log.info(
                        f"Writing Imsim Configuration file to {imsim_config_path}"
                    )
                    self.imsim_cmpt.write_yaml_and_run_imsim(
                        imsim_config_path, imsim_config_yaml
                    )

    def _calc_wf_err_from_img(
        self,
        obs_metadata: ObsMetadata,
        butler_root_path: str,
        inst_name: str,
        num_pro: int = 1,
        pipeline_file: str | None = None,
        filter_type_name: str = "",
    ) -> list[SensorWavefrontError]:
        """Calculate the wavefront error from the images generated by PhoSim.

        Parameters
        ----------
        obs_metadata : lsst.ts.imsim.ObsMetadata object
            Observation metadata.
        butler_root_path : str
            Path to the butler repository.
        inst_name : str
            Instrument name.
        num_pro : int, optional
            Number of processor to run DM pipeline. (the default is 1.)
        pipeline_file : str or None, optional
            Path to existing pipeline yaml file to use.
            If None then the code will write its own default pipeline yaml.
            (The default is None.)
        filter_type_name : str, optional
            Filter type name: ref (or ''), u, g, r, i, z, or y.

        Returns
        -------
        list[lsst.ts.imsim.utils.SensorWavefrontError]
            List of SensorWavefrontError object.
        """

        # Ingest images into butler gen3
        self.ingest_data(butler_root_path=butler_root_path, inst_name=inst_name)

        list_of_wf_err = self.run_wep(
            obs_metadata.seq_num,
            butler_root_path,
            inst_name,
            num_pro=num_pro,
            pipeline_file=pipeline_file,
            filter_type_name=filter_type_name,
        )

        return list_of_wf_err

    def run_wep(
        self,
        seq_num: int,
        butler_root_path: str,
        inst_name: str,
        num_pro: int = 1,
        pipeline_file: str | None = None,
        filter_type_name: str = "",
    ) -> list[SensorWavefrontError]:
        """Run wavefront estimation pipeline task for wavefront sensors.

        Parameters
        ----------
        seq_num : int
            Observation id.
        butler_root_path : str
            Path to the butler gen3 repos.
        inst_name : str
            Instrument name.
        num_pro : int, optional
            Number of processor to run DM pipeline. (the default is 1.)
        pipeline_file : str or None, optional
            Path to existing pipeline yaml file to use.
            If None then the code will write its own default pipeline yaml.
            (The default is None.)
        filter_type_name : str, optional
            Filter type name: ref (or ''), u, g, r, i, z, or y.

        Returns
        -------
        list[lsst.ts.imsim.utils.SensorWavefrontError]
            List of SensorWavefrontError with the results of the wavefront
            estimation pipeline for each sensor.
        """

        butler_inst_name = "Cam"
        if pipeline_file is None:
            pipeline_yaml = f"{inst_name}Pipeline.yaml"
            pipeline_yaml_path = os.path.join(butler_root_path, pipeline_yaml)
            self.write_wep_configuration(
                inst_name, pipeline_yaml_path, filter_type_name
            )
        else:
            pipeline_yaml_path = pipeline_file

        butler = dafButler.Butler(butler_root_path)

        if f"LSST{butler_inst_name}/calib" not in butler.registry.queryCollections():
            self.log.info("Ingesting curated calibrations.")

            runProgram(
                f"butler write-curated-calibrations {butler_root_path} lsst.obs.lsst.Lsst{butler_inst_name}"
            )
        # Sequence number or seq_num is an integer number
        # associated with each image taken in a single day.
        # The limit for seq_num is 5 digits,
        # set by the expectation that no more than 100K images
        # could be taken in a single day (i.e. no more than 1/sec).
        if inst_name == "lsst":
            runProgram(
                f"pipetask run -b {butler_root_path} "
                f"-i refcats,LSST{butler_inst_name}/raw/all,LSST{butler_inst_name}/calib/unbounded "
                f"--instrument lsst.obs.lsst.Lsst{butler_inst_name} "
                f"--register-dataset-types --output-run ts_imsim_{seq_num} -p {pipeline_yaml_path} -d "
                f'"visit.seq_num IN ({seq_num})" -j {num_pro}'
            )
        elif inst_name == "lsstfam":
            runProgram(
                f"pipetask run -b {butler_root_path} "
                f"-i refcats,LSST{butler_inst_name}/raw/all,LSST{butler_inst_name}/calib/unbounded "
                f"--instrument lsst.obs.lsst.Lsst{butler_inst_name} "
                f"--register-dataset-types --output-run ts_imsim_{seq_num} -p {pipeline_yaml_path} -d "
                f'"visit.seq_num IN ({seq_num-1}, {seq_num})" -j {num_pro}'
            )

        # Need to redefine butler because the database changed.
        butler = dafButler.Butler(butler_root_path)

        dataset_refs = butler.registry.queryDatasets(
            datasetType="zernikeEstimateAvg", collections=[f"ts_imsim_{seq_num}"]
        )

        # Get the map for detector Id to detector name
        camera = butler.get(
            "camera",
            {"instrument": f"LSST{butler_inst_name}"},
            collections=[f"LSST{butler_inst_name}/calib/unbounded"],
        )
        det_id_map = camera.getIdMap()
        det_name_map = camera.getNameMap()

        list_of_wf_err = []

        for dataset in dataset_refs:
            data_id = {
                "instrument": dataset.dataId["instrument"],
                "detector": dataset.dataId["detector"],
                "visit": dataset.dataId["visit"],
            }

            zer_coeff = butler.get(
                "zernikeEstimateAvg",
                dataId=data_id,
                collections=[f"ts_imsim_{seq_num}"],
            )

            sensor_wavefront_data = SensorWavefrontError()
            # Rotate from CCS to ZCS/PCS (Phosim Coordinate System)
            rotated_name = self.rotate_det_name_ccs_to_zcs(
                det_id_map[dataset.dataId["detector"]].getName()
            )
            sensor_wavefront_data.sensor_name = rotated_name
            sensor_wavefront_data.sensor_id = det_name_map[rotated_name].getId()
            sensor_wavefront_data.annular_zernike_poly = zer_coeff

            list_of_wf_err.append(sensor_wavefront_data)

        return list_of_wf_err

    def rotate_det_name_ccs_to_zcs(self, det_name: str) -> str:
        """Rotate the sensor name from CCS (Camera Coordinate System)
        to the ZCS (Zemax Coordinate System) that OFC expects.

        Parameters
        ----------
        det_name : str
            Detector name in CCS.

        Returns
        -------
        str
            Detector expected at that location in ZCS.
        """
        raft, sensor = det_name.split("_")
        cam_rx = int(raft[1]) - 2
        cam_ry = int(raft[2]) - 2

        zcs_rx = -cam_rx
        zcs_ry = cam_ry

        zcs_raft = f"R{zcs_rx + 2}{zcs_ry + 2}"

        if sensor.startswith("SW"):
            zcs_sensor = sensor
        else:
            cam_sx = int(sensor[1]) - 1
            cam_sy = int(sensor[2]) - 1
            zcs_sx = -cam_sx
            zcs_sy = cam_sy
            zcs_sensor = f"S{zcs_sx + 1}{zcs_sy + 1}"

        return f"{zcs_raft}_{zcs_sensor}"

    def write_wep_configuration(
        self, inst_name: str, pipeline_yaml_path: str, filter_type_name: str
    ) -> None:
        """Write wavefront estimation pipeline task configuration.

        Parameters
        ----------
        inst_name : str
            Name of the instrument this configuration is intended for.
        pipeline_yaml_path : str
            Path where the pipeline task configuration yaml file
            should be saved.
        filter_type_name : str
            Filter type name: ref (or ''), u, g, r, i, z, or y.
        """

        butler_inst_name = "Cam"

        # Remap reference filter
        filter_type_name = self.map_filter_ref_to_g(filter_type_name)

        with open(pipeline_yaml_path, "w") as fp:
            fp.write(
                f"""# This yaml file is used to define the tasks and configuration of
# a Gen 3 pipeline used for testing in ts_wep.
description: wep basic processing test pipeline
# Here we specify the corresponding instrument for the data we
# will be using.
instrument: lsst.obs.lsst.Lsst{butler_inst_name}
# Use imported instrument configuration
imports:
  - location: {getWepConfigDir()}/cwfs/instData/{inst_name}/instParamPipeConfig.yaml
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

    def run_img(
        self,
        inst: str,
        filter_type_name: str,
        rot_cam_in_deg: float,
        boresight: list[float],
        mjd: float,
        base_output_dir: str,
        path_sky_file: str,
        do_erase_dir_content: bool,
        sky_seed: int,
        pert_seed: int,
        iter_num: int,
        pipeline_file: str,
        imsim_config_pointer_file: str,
        turn_off_sky_background: bool,
        turn_off_atmosphere: bool,
        turn_off_wavefront_estimates: bool,
        num_pro: int,
        raw_seeing: float,
    ) -> None:
        """Run the simulation of images.

        Parameters
        ----------
        inst : str
            Instrument to use: currently only lsst.
        filter_type_name : str
            Filter type name: ref, u, g, r, i, z, or y.
        rot_cam_in_deg : float
            The camera rotation angle in degree (-90 to 90).
        boresight : list[float]
            Boresight [ra, dec] in degree.
        mjd : float
            MJD of the observation.
        base_output_dir : str
            Base output directory.
        path_sky_file : str
            Path to the sky file.
        do_erase_dir_content : bool
            Do the erase of the content of base output directory or not.
        sky_seed : int
            Random seed for the sky background.
        pert_seed : int
            Random seed for the perturbations.
        iter_num : int
            Number of closed-loop iteration.
        pipeline_file : str
            Path to existing pipeline yaml file to use. If empty string
            then the code will write its own default pipeline yaml.
        imsim_config_pointer_file : str
            Path to pointer file with locations of yaml configuration
            files for imsim submodules. If empty string then the code
            will use the default in policy/config for the given inst.
        turn_off_sky_background : bool
            If set to True then the closed loop will simulate images
            without sky background.
        turn_off_atmosphere : bool
            If set to True then will turn off the imsim atmosphere.
        turn_off_wavefront_estimates : bool
            If set to True then will run the closed loop only with
            the OPD.fits files and not simulated images.
        num_pro : int
            Number of processors to use.
        raw_seeing : float
            Raw seeing in arcsec.
        """
        cam_type, inst_name = self.get_cam_type_and_inst_name(inst)
        base_output_dir = self.check_and_create_base_output_dir(base_output_dir)
        if do_erase_dir_content:
            self.erase_directory_content(base_output_dir)
        self.check_boresight(boresight)
        self.boresight_ra = boresight[0]
        self.boresight_dec = boresight[1]
        self.boresight_rot_ang = rot_cam_in_deg
        # Remap the reference filter to g
        filter_type_name = self.map_filter_ref_to_g(filter_type_name)

        if turn_off_wavefront_estimates is True:
            self.use_ccd_img = False

        obs_metadata = ObsMetadata(
            ra=self.boresight_ra,
            dec=self.boresight_dec,
            band=filter_type_name.lower(),
            rotator_angle=self.boresight_rot_ang,
            mjd=mjd,
            raw_seeing=raw_seeing,
        )

        # Configure the components
        self.config_sky_sim(
            inst_name, obs_metadata, path_sky_file=path_sky_file, star_mag=15
        )
        self.config_ofc_calc(inst_name)
        self.imsim_cmpt = ImsimCmpt()

        # If path_sky_file using default OPD positions write this to disk
        # so that the Butler can load it later
        if path_sky_file == "":
            path_sky_file = os.path.join(base_output_dir, "sky_info.txt")
            self.sky_sim.export_sky_to_file(path_sky_file)
            self.log.info(f"Wrote new sky file to {path_sky_file}.")

        # generate butler gen3 repo if needed
        butler_root_path = os.path.join(base_output_dir, "imsimData")

        if self.use_ccd_img:
            self.generate_butler(butler_root_path, inst_name)
            self.generate_ref_catalog(
                inst_name=inst_name,
                butler_root_path=butler_root_path,
                path_sky_file=path_sky_file,
                filter_type_name=filter_type_name,
            )

        if inst_name == "lsst":
            # Append equal weights for CWFS fields to OFC data
            # Assign equal normalized weights to each of the
            # four corner wavefront sensor pairs.
            self.ofc_calc.ofc_data.normalized_image_quality_weight = np.append(
                self.ofc_calc.ofc_data.normalized_image_quality_weight,
                [0.25, 0.25, 0.25, 0.25],
            )

        self._run_sim(
            cam_type,
            inst_name,
            obs_metadata,
            base_output_dir,
            butler_root_path,
            sky_seed,
            pert_seed,
            iter_num,
            num_pro=num_pro,
            pipeline_file=pipeline_file,
            imsim_config_pointer_file=imsim_config_pointer_file,
            turn_off_sky_background=turn_off_sky_background,
            turn_off_atmosphere=turn_off_atmosphere,
        )

    def generate_butler(self, butler_root_path: str, inst_name: str) -> None:
        """Generate butler gen3.

        Parameters
        ----------
        butler_root_path: `str`
            Path to where the butler repository should be created.
        inst_name: `str`
            Name of the instrument.
        """

        self.log.info(f"Generating butler gen3 in {butler_root_path} for {inst_name}")

        runProgram(f"butler create {butler_root_path}")

        self.log.debug("Registering LsstCam")
        runProgram(
            f"butler register-instrument {butler_root_path} lsst.obs.lsst.LsstCam"
        )

    def generate_ref_catalog(
        self, butler_root_path: str, path_sky_file: str, filter_type_name: str
    ) -> None:
        """Generate reference star catalog.

        Parameters
        ----------
        butler_root_path: `str`
            Path to the butler gen3 repository.
        path_sky_file: `str`
            Path to the catalog star file.
        filter_type_name : str
            Filter type name: ref (or ''), u, g, r, i, z, or y.
        """
        self.log.debug("Creating reference catalog.")

        cat_dir = os.path.join(butler_root_path, "skydata")
        sky_file_name = os.path.join(cat_dir, "sky_data.csv")
        cat_config_file_name = os.path.join(cat_dir, "cat.cfg")
        sky_ecsv_file_name = os.path.join(cat_dir, "filename_to_htm.ecsv")
        cat_log_file_name = os.path.join(cat_dir, "convert.log")
        os.mkdir(cat_dir)

        # Read sky file and convert it to csv
        sky_data = astropy.io.ascii.read(path_sky_file)

        # Constructing the catalog of stars to use in the wavefront estimation
        # pipeline. It is used for target
        # selection, and affects magnitude limits
        # as set in generateDonutCatalogWcsTask pipeline yaml file
        sky_data.rename_column("Mag", filter_type_name)

        sky_data.write(sky_file_name, format="csv", overwrite=True)

        with open(cat_config_file_name, "w") as fp:
            fp.write(
                f"""config.ra_name='Ra'
config.dec_name='Dec'
config.id_name='Id'
config.mag_column_list=['{filter_type_name}']
config.dataset_config.ref_dataset_name='ref_cat'
"""
            )

        runProgram(
            f"convertReferenceCatalog {cat_dir} {cat_config_file_name} {sky_file_name}",
            stdout=cat_log_file_name,
            stderr=cat_log_file_name,
        )

        runProgram(
            f"butler register-dataset-type {butler_root_path} cal_ref_cat SimpleCatalog htm7"
        )

        runProgram(
            f"butler ingest-files -t direct {butler_root_path} cal_ref_cat refcats {sky_ecsv_file_name}"
        )

    def ingest_data(self, butler_root_path: str, inst_name: str) -> None:
        """Ingest data into a gen3 data Butler.

        Parameters
        ----------
        butler_root_path : str
            Path to the butler repository.
        inst_name : str
            Instrument name.
        """
        output_img_dir = self.imsimCmpt.outputImgDir
        files = " ".join(glob(os.path.join(output_img_dir, "amp*")))

        if inst_name in ["lsst", "lsstfam"]:
            runProgram(f"butler ingest-raws {butler_root_path} {files}")

        runProgram(f"butler define-visits {butler_root_path} lsst.obs.lsst.LsstCam")

    def erase_directory_content(self, target_dir: str) -> None:
        """Erase the directory content.

        Parameters
        ----------
        target_dir : str
            Target directory.
        """

        for file_on in os.listdir(target_dir):
            file_path = os.path.join(target_dir, file_on)
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    @staticmethod
    def set_default_parser(parser: ArgumentParser) -> ArgumentParser:
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
            help="Instrument to use: currently lsst or lsstfam. (default: lsst)",
        )

        parser.add_argument(
            "--filter_type",
            type=str,
            default="",
            help="Filter type to use: u, g, r, i, z, y or empty string for "
            "reference wavelength. (default: '')",
        )

        parser.add_argument(
            "--rot_cam",
            type=float,
            default=0.0,
            help="Rotate camera (degree) in counter-clockwise direction. (default: 0.0)",
        )

        parser.add_argument("--output", type=str, default="", help="Output directory.")

        parser.add_argument(
            "--log_level", type=int, default=logging.INFO, help="Log level."
        )

        parser.add_argument(
            "--clobber",
            default=False,
            action="store_true",
            help="Delete existing output directory.",
        )

        parser.add_argument(
            "--config_pointer_file",
            type=str,
            default="",
            help="Imsim Configuration Pointer File.",
        )

        parser.add_argument(
            "--sky_seed",
            type=int,
            default=42,
            help="Random seed for imsim sky (default: 42).",
        )

        parser.add_argument(
            "--pert_seed",
            type=int,
            default=11,
            help="Random seed for m1m3_lut fractional actuator random error. "
            "(Default: 11)",
        )

        parser.add_argument(
            "--iter_num",
            type=int,
            default=5,
            help="Number of closed-loop iterations. (default: 5)",
        )

        parser.add_argument(
            "--pipeline_file",
            type=str,
            default="",
            help="""
            Location of user-specified pipeline configuration file.
            If left as empty string the code will create a default file.
            (default: '')
            """,
        )

        parser.add_argument(
            "--num_proc",
            type=int,
            default=1,
            help="Number of processor to run imSim and DM pipeline. (default: 1)",
        )

        return parser

    @staticmethod
    def set_img_parser(parser: ArgumentParser) -> ArgumentParser:
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
            "--boresight_deg",
            type=float,
            nargs=2,
            default=[0, 0],
            help="Boresight [ra, dec] in degree. The default is [0, 0].",
        )

        parser.add_argument(
            "--sky_file",
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
            "--turn_off_sky_background",
            action="store_true",
            help="Turn sky brightness model off.",
        )

        parser.add_argument(
            "--turn_off_atmosphere", action="store_true", help="Turn atmosphere off."
        )

        parser.add_argument(
            "--turn_off_wavefront_estimates",
            action="store_true",
            help="""
                 Run with true wavefront values only.
                 Turns off images generation and running CCDs through WEP.
                 """,
        )

        parser.add_argument(
            "--raw_seeing",
            type=float,
            default=0.5,
            help="Raw seeing in arcsec (default: 0.5).",
        )

        return parser
