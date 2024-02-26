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

__all__ = ["ImsimCmpt"]

import os
from typing import Any

import numpy as np
import yaml
from astropy.io import fits
from lsst.afw import cameraGeom
from lsst.ts.imsim.obs_metadata import ObsMetadata
from lsst.ts.imsim.opd_metrology import OpdMetrology
from lsst.ts.imsim.sky_sim import SkySim
from lsst.ts.imsim.utils import (
    ModifiedEnvironment,
    SensorWavefrontError,
    get_config_dir,
    get_zk_from_file,
    make_dir,
)
from lsst.ts.ofc.ofc_data.base_ofc_data import BaseOFCData
from lsst.ts.wep.utils import CamType, runProgram
from scipy.ndimage import rotate


class ImsimCmpt:
    """Class to take configurations for each imsim component and
    generate full imsim configuration files.
    """

    def __init__(self) -> None:
        # Output directories
        self._output_dir = None
        self._output_img_dir = None

        # OPD information
        self.opd_file_path = None
        self.opd_metr = OpdMetrology()

        # Specify number of Zernikes
        self.num_of_zk = 19

        # AOS Degrees of Freedom
        self.num_of_dof = 50
        self.dof_in_um = np.zeros(self.num_of_dof, dtype=float)

    @property
    def output_dir(self) -> str:
        return self._output_dir

    @output_dir.setter
    def output_dir(self, new_output_dir: str) -> None:
        """
        Set the closed loop output directory and make it if it
        does not yet exits.

        Parameters
        ----------
        new_output_dir : str
            Path for output image directory.
        """
        make_dir(new_output_dir)
        self._output_dir = new_output_dir

    @property
    def output_img_dir(self) -> str:
        return self._output_img_dir

    @output_img_dir.setter
    def output_img_dir(self, new_output_img_dir: str) -> None:
        """
        Set the output image directory and make it if it
        does not yet exits.

        Parameters
        ----------
        new_output_img_dir : str
            Path for output image directory.
        """
        make_dir(new_output_img_dir)
        self._output_img_dir = new_output_img_dir

    def _verify_pointer_file(
        self, file_pointer_info: dict, config_sections: list[str]
    ) -> None:
        """Verify that pointer file has filepaths for all needed
        sections of a complete imsim configuration file.

        Parameters
        ----------
        file_pointer_info : dict
            Dictionary holding pointer types as keys and file path
            information as values.
        config_sections : list
            List of required submodules for ImSim configuration.

        Raises
        ------
        ValueError
            Raised when file_pointer_info is missing a required submodule.
        """
        file_pointer_keys = list(file_pointer_info.keys())
        for section_name in config_sections:
            if section_name not in file_pointer_keys:
                raise ValueError(
                    f"Config pointer file missing filepath for {section_name}."
                )

    def assemble_config_yaml(
        self,
        obs_metadata: ObsMetadata,
        config_pointer_file: str,
        cam_type: CamType,
        required_modules_file: str | None = None,
    ) -> dict[str, Any]:
        """
        Parameters
        ----------
        obs_metadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.
        config_pointer_file : str
            The location of the config pointer file that specifies
            the location of configurations files for each imsim
            component.
        cam_type : lsst.ts.wep.utils.CamType
            Camera type.
        required_modules_file : str or None
            Path to yaml file with required imsim modules. If None,
            then will use policy/requiredModulesDefault.yaml.
            (The default is None.)
        """

        if required_modules_file is None:
            required_modules_file = os.path.join(
                get_config_dir(), "requiredModulesDefault.yaml"
            )

        obs_info_text = self.convert_obs_metadata_to_text(obs_metadata)

        config_sections = ["input", "gal", "image", "psf", "stamp", "output"]

        with open(config_pointer_file, "r") as f:
            file_pointer_info = yaml.safe_load(f)
        self._verify_pointer_file(file_pointer_info, config_sections)

        with open(required_modules_file, "r") as required_modules:
            required_modules_text = required_modules.read()

        # Assemble full configuration as text because of the
        # aliases in the individual yaml components that
        # fill in values from "evalVariables".
        full_config_text = ""
        full_config_text += required_modules_text + "\n"
        full_config_text += obs_info_text + "\n"

        # Treat input section first differently since it has multiple parts
        input_section_text = "input:\n"
        for sub_section in ["atm_psf", "sky_model", "telescope", "vignetting"]:
            with open(
                file_pointer_info["input"][sub_section].format(**os.environ), "r"
            ) as sub_file:
                for line in sub_file.readlines():
                    input_section_text += "  " + line
        full_config_text += input_section_text + "\n"
        # Now move on to other sections skipping input section
        for section_name in config_sections:
            if section_name == "input":
                continue
            else:
                with open(
                    file_pointer_info[section_name].format(**os.environ), "r"
                ) as section_file:
                    full_config_text += section_file.read() + "\n"
            # Handle output specially since we need to add header info and OPD
            if section_name == "output":
                full_config_text += f"  dir: {self.output_img_dir}\n\n"
                # Add additional header text
                full_config_text += self.add_config_header(obs_metadata) + "\n"
                # Add OPD
                full_config_text += self.format_opd_text(obs_metadata, cam_type)

        # Assemble as yaml
        full_config_yaml = yaml.safe_load(full_config_text)

        # Add in OFC corrections in um and arcsec
        dof_in_um = self.dof_in_um
        full_config_yaml["input"]["telescope"]["fea"]["aos_dof"] = {
            "dof": dof_in_um.tolist(),
            "type": "List",
        }
        return full_config_yaml

    def convert_obs_metadata_to_text(self, obs_metadata: ObsMetadata) -> str:
        """
        Write out the evalVariables section of the config file.

        Parameters
        ----------
        obs_metadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.

        Returns
        -------
        str
            eval_variables for ImSim config
        """
        obs_variables_text = "eval_variables:\n"
        obs_variables_text += "  cboresight:\n"
        obs_variables_text += "    type: RADec\n"
        obs_variables_text += f"    ra: &ra {obs_metadata.ra} deg\n"
        obs_variables_text += f"    dec: &dec {obs_metadata.dec} deg\n"
        obs_variables_text += f"  sband: &band {obs_metadata.band}\n"
        obs_variables_text += f"  azenith: &zenith {obs_metadata.zenith:.6f} deg\n"
        obs_variables_text += f"  artp: &rtp {obs_metadata.rotator_angle} deg\n"
        obs_variables_text += f"  fexptime: &exptime {obs_metadata.exp_time}\n"
        obs_variables_text += f"  fmjd: &mjd {obs_metadata.mjd}\n"
        obs_variables_text += f"  frawSeeing: &rawSeeing {obs_metadata.raw_seeing}\n"
        obs_variables_text += f"  iseqnum: &seqnum {obs_metadata.seq_num}\n"
        obs_variables_text += f"  sobsid: &obsid {obs_metadata.obs_id}\n"
        obs_variables_text += f"  aalt: &alt {obs_metadata.alt:.6f} deg\n"
        obs_variables_text += f"  aaz: &az {obs_metadata.az:.6f} deg\n"

        return obs_variables_text

    def format_opd_text(self, obs_metadata: ObsMetadata, cam_type: CamType) -> str:
        """
        Write out the OPD section of the config file.

        Parameters
        ----------
        obs_metadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.
        cam_type : lsst.ts.wep.utils.CamType
            Camera type.

        Returns
        -------
        str
            OPD information for ImSim config
        """
        opd_section_text = "  opd:\n"
        opd_section_text += "    file_name: opd.fits\n"
        opd_section_text += "    nx: 255\n"
        opd_section_text += "    rotTelPos: *rtp\n"
        opd_section_text += "    jmax: 28\n"
        opd_section_text += "    eps: 0.61\n"
        opd_section_text += "    projection: gnomonic\n"
        opd_section_text += f"    wavelength: {BaseOFCData().eff_wavelength[obs_metadata.band.upper()]*1e3}\n"
        opd_section_text += "    fields:\n"

        # Get the locations for the OPD from OPD Metrology
        self.opd_metr.set_wgt_and_field_xy_of_gq(cam_type)
        for thx, thy in zip(self.opd_metr.field_x, self.opd_metr.field_y):
            opd_section_text += f"      - {{thx: {thx} deg, thy: {thy} deg}}\n"

        return opd_section_text

    def add_config_header(self, obs_metadata: ObsMetadata) -> str:
        """
        Write out the header section of the config file.

        Parameters
        ----------
        obs_metadata : lsst.ts.imsim.ObsMetadata
            ObsMetadata dataclass object with observation information.

        Returns
        -------
        str
            Header information for ImSim config
        """
        header_text = "  header:\n"
        header_text += "    mjd: *mjd\n"
        header_text += (
            f"    observationStartMJD: {obs_metadata.mjd - (15/(60*60*24))}\n"
        )
        header_text += "    seqnum: *seqnum\n"
        header_text += "    band: *band\n"
        header_text += f"    fieldRA: {obs_metadata.ra}\n"
        header_text += f"    fieldDec: {obs_metadata.dec}\n"
        header_text += f"    rotTelPos: {obs_metadata.rotator_angle}\n"
        header_text += (
            f"    airmass: {1.0/np.cos(np.radians(obs_metadata.zenith)):.6f}\n"
        )
        header_text += f"    focusZ: {obs_metadata.focus_z}\n"
        header_text += f"    seeing: {obs_metadata.raw_seeing}\n"

        return header_text

    def run_imsim(self, config_file_path: str, imsim_log_file: str = "") -> None:
        """
        Run imSim with the given configuration file.

        Parameters
        ----------
        config_file_path : str
            Path to the imSim configuration yaml file.
        imsim_log_file : str, optional
            Path to save imSim log output. If empty string then the
            code will just default to sending it to stdout.
            (The default is "".)
        """
        # For mysterious reasons, having the KMP_INIT_AT_FORK environment
        # variable set causes problems in the imSim subprocess.  It may be
        # related to https://github.com/numpy/numpy/issues/11734.  Note that
        # this variable gets set when sklearn is imported, which happens inside
        # of ts_wep before ever getting to this method, so we can't just unset
        # it in the environment ahead of time.
        with ModifiedEnvironment(KMP_INIT_AT_FORK=None):
            if imsim_log_file == "":
                runProgram(f"galsim {config_file_path}")
            else:
                runProgram(
                    f"galsim {config_file_path} -l {os.path.join(self.output_dir, imsim_log_file)} -v 2"
                )

    def write_yaml_and_run_imsim(
        self,
        config_path: str,
        config_yaml: dict[str, Any],
        imsim_log_file: str = "",
    ) -> None:
        """Write yaml config file and run Imsim.

        Parameters
        ----------
        config_path : str
            Path to write config yaml file.
        config_yaml : dict
            Dictionary that contains imsim config details to write to yaml.
        imsim_log_file : str, optional
            Path to save imSim log output. If empty string then the
            code will just default to sending it to stdout.
            (The default is "".)
        """

        with open(config_path, "w") as file:
            yaml.safe_dump(config_yaml, file)
        self.run_imsim(config_path, imsim_log_file=imsim_log_file)

    def add_sources_to_config(
        self, config_yaml: dict[str, Any], inst_cat_path: str, use_ccd_img: bool = True
    ) -> dict[str, Any]:
        """Add source information to config. If using CCD it will add
        the instance catalog details. If only using OPD it will remove
        the instance catalog info so that we are not generating
        unnecessary images.

        Parameters
        ----------
        config_yaml : dict
            Dictionary that contains imsim config details to write to yaml.
        inst_cat_path : str
            Path to instance catalog file.
        use_ccd_img : bool, optional
            When outputting CCD images this is set to True so that the instance
            catalog information is included in the configuration. For OPD only
            this should be set to False. (The default is True.)

        Returns
        -------
        dict
            Updated yaml configuration.
        """

        if use_ccd_img:
            config_yaml["image"].pop("nobjects")
            config_yaml["image"]["world_pos"] = {"type": "InstCatWorldPos"}
            config_yaml["gal"] = {"type": "InstCatObj"}
            config_yaml["input"]["instance_catalog"] = {
                "file_name": inst_cat_path,
                "sed_dir": "$os.environ.get('SIMS_SED_LIBRARY_DIR')",
            }
        else:
            config_yaml["image"]["nobjects"] = 0
            # Only create one empty amplifier file and one OPD.fits file
            # when running OPD only
            config_yaml["output"]["nfiles"] = 1

        return config_yaml

    def gen_instance_catalog(self, sky_sim: SkySim) -> str:
        """
        Generate the instance catalog.

        Parameters
        ----------
        sky_sim : lsst.ts.imsim.skySim.SkySim
            The SkySim object with the stars to simulate.

        Returns
        -------
        str
            Text to write to instance catalog file.
        """

        content = ""
        content += self.gen_inst_cat_stars(sky_sim)
        return content

    def gen_inst_cat_stars(self, sky_sim: SkySim) -> str:
        """
        Add stars to the instance catalog.

        Parameters
        ----------
        sky_sim : lsst.ts.imsim.skySim.SkySim
            The SkySim object with the stars to simulate.

        Returns
        -------
        str
            Instance catalog text for the stars in the SkySim catalog.
        """
        content = ""
        for id, ra, dec, mag in zip(
            sky_sim.star_id, sky_sim.ra, sky_sim.dec, sky_sim.mag
        ):
            content += self.generate_star(id, ra, dec, mag)

        return content

    def generate_star(
        self,
        star_id: int,
        ra: float,
        dec: float,
        mag_norm: float,
        sed_name: str = "flatSED/sed_flat.txt.gz",
        redshift: float = 0,
        gamma_1: float = 0,
        gamma_2: float = 0,
        kappa: float = 0,
        delta_ra: float = 0,
        delta_dec: float = 0,
        source_type: str = "point",
    ) -> str:
        """Generate the star source.

        Parameters
        ----------
        star_id : int
            Star Id.
        ra : float
            The right ascension of the center of the object or image in
            decimal degrees.
        dec : float
            The declination of the center of the object in decimal degrees.
        mag_norm : float
            The normalization of the flux of the object in AB magnitudes
            at (500 nm)/(1+z) (which is roughly equivalent to V (AB) or
            g (AB)).
        sed_name : str, optional
            The name of the SED file with a file path that is relative to the
            $SIMS_SED_LIBRARY_DIR directory environment variable
            required for imSim..
            (The default is "flatSED/sed_flat.txt.gz")
        redshift : float, optional
            The redshift (or blueshift) of the object. Note that the SED does
            not need to be redshifted if using this. (the default is 0.)
        gamma_1 : float, optional
            The value of the shear parameter gamma1 used in weak lensing.
            (the default is 0.)
        gamma_2 : float, optional
            The value of the shear parameter gamma2 used in weak lensing.
            (the default is 0.)
        kappa : float, optional
            The value of the magnification parameter in weak lensing. (the
            default is 0.)
        delta_ra : float, optional
            The value of the declination offset in radians. This can be used
            either for weak lensing or objects that moved from another
            exposure if you do not want to change the source position in the
            first two columns. (the default is 0.)
        delta_dec : float, optional
            The value of the declination offset in radians. This can be used
            either for weak lensing or objects that moved from another
            exposure if you do not want to change the source position in the
            first two columns. (the default is 0.)
        source_type : str, optional
            The name of the spatial model to be used as defined below. (the
            default is "point".)

        Returns
        -------
        str
            Instance catalog entry for star.
        """

        content = "object %2d\t%9.6f\t%9.6f %9.6f %s " % (
            star_id,
            ra,
            dec,
            mag_norm,
            sed_name,
        )
        content += "%.1f %.1f %.1f %.1f %.1f %.1f %s none none \n" % (
            redshift,
            gamma_1,
            gamma_2,
            kappa,
            delta_ra,
            delta_dec,
            source_type,
        )

        return content

    def analyze_opd_data(
        self,
        cam_type: CamType,
        zk_file_name: str = "opd.zer",
        rot_opd_in_deg: float = 0.0,
        pssn_file_name: str = "PSSN.txt",
    ) -> None:
        """Analyze the OPD data.

        Rotate OPD to simulate the output by rotated camera. When anaylzing the
        PSSN, the unrotated OPD is used.

        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        cam_type : lsst.ts.wep.utils.CamType
            Camera type.
        zk_file_name : str, optional
            OPD in zk file name. (the default is "opd.zer".)
        rot_opd_in_deg : float, optional
            Rotate OPD in degree in the counter-clockwise direction. (the
            default is 0.0.)
        pssn_file_name : str, optional
            PSSN file name. (the default is "PSSN.txt".)
        """

        # Set the weighting ratio and field positions of OPD
        if cam_type == CamType.LsstCam:
            self.opd_metr.set_default_lsst_wfs_gq()
        else:
            self.opd_metr.set_wgt_and_field_xy_of_gq(cam_type)
        num_opd = len(self.opd_metr.field_x)

        self.opd_file_path = os.path.join(self.output_img_dir, "opd.fits")
        self._write_opd_zk_file(zk_file_name, rot_opd_in_deg, num_opd)
        self._write_opd_pssn_file(pssn_file_name, num_opd)

    def _write_opd_zk_file(
        self, zk_file_name: str, rot_opd_in_deg: float, num_opd: int
    ) -> None:
        """Write the OPD in zk file.

        OPD: optical path difference.

        Parameters
        ----------
        zk_file_name : str
            OPD in zk file name.
        rot_opd_in_deg : float
            Rotate OPD in degree in the counter-clockwise direction.
        num_opd : int
            Number of OPD positions calculated.
        """

        file_path = os.path.join(self.output_img_dir, zk_file_name)
        opd_data = self._map_opd_to_zk(rot_opd_in_deg, num_opd)
        file_txt = (
            "# The followings are OPD in rotation angle of %.2f degree in nm from z4 to z22:\n"
            % rot_opd_in_deg
        )
        for sensor_id, opd_zk in zip(self.opd_metr.sensor_ids, opd_data):
            zk_str = f"{sensor_id}: {opd_zk}\n"
            file_txt += zk_str
        with open(file_path, "w") as file:
            file.write(file_txt)

    def _map_opd_to_zk(self, rot_opd_in_deg: float, num_opd: int) -> np.ndarray:
        """Map the OPD to the basis of annular Zernike polynomial (Zk).

        OPD: optical path difference.

        Parameters
        ----------
        rot_opd_in_deg : float
            Rotate OPD in degree in the counter-clockwise direction.
        num_opd : int
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
        opd_data = np.zeros((num_opd, self.num_of_zk))
        for idx in range(num_opd):
            opd = fits.getdata(self.opd_file_path, idx)

            # Rotate OPD if needed
            if rot_opd_in_deg != 0:
                opd_rot = opd.copy()
                # Since to rotate the opd we need to substitue the nan values
                # for zeros, we need to find the minimum value of the opd
                # excluding the nan values. Then after the rotation we will
                # discard the values that are smaller than the minimum value.
                # Note that we use order = 0 to avoid interpolation errors.
                min_value = np.nanmin(np.abs(opd_rot))
                opd_rot[np.isnan(opd_rot)] = 0.0
                opd_rot = rotate(opd_rot, rot_opd_in_deg, reshape=False, order=0)
                opd_rot[np.abs(opd_rot) <= min_value] = np.nan
            else:
                opd_rot = opd

            # z1 to z22 (22 terms)
            zk = self.opd_metr.get_zk_from_opd(opd_map=opd_rot)[0]

            # Only need to collect z4 to z22
            init_idx = 3
            opd_data[idx, :] = zk[init_idx : init_idx + self.num_of_zk]

        return opd_data

    def _write_opd_pssn_file(self, pssn_file_name: str, num_opd: int) -> None:
        """Write the OPD PSSN in file.

        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssn_file_name : str
            PSSN file name.
        num_opd : int
            Number of OPD positions calculated.
        """

        # Calculate the PSSN
        pssn_list, gq_eff_pssn = self._calc_pssn_opd(num_opd)

        # Calculate the FWHM
        eff_fwhm_list, gq_eff_fwhm = self._calc_eff_fwhm_opd(pssn_list)

        # Append the list to write the data into file
        pssn_list.append(gq_eff_pssn)
        eff_fwhm_list.append(gq_eff_fwhm)

        # Stack the data
        data = np.vstack((pssn_list, eff_fwhm_list))

        # Write to file
        file_path = os.path.join(self.output_img_dir, pssn_file_name)
        header = "The followings are PSSN and FWHM (in arcsec) data. The final number is the GQ value."
        np.savetxt(file_path, data, header=header)

    def _calc_pssn_opd(self, num_opd: int) -> tuple[list[float], float]:
        """Calculate the PSSN of OPD.

        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.
        GQ: Gaussian quadrature.

        Parameters
        ----------
        num_opd : int
            Number of OPD positions calculated.

        Returns
        -------
        list
            PSSN list.
        float
            GQ effective PSSN.
        """

        pssn_list = []
        for idx in range(num_opd):
            wavelength_in_um = fits.getheader(self.opd_file_path, idx)["WAVELEN"] * 1e-3
            # OPD data needs to be in microns. Imsim output is nm.
            pssn = self.opd_metr.calc_pssn(
                wavelength_in_um, opd_map=fits.getdata(self.opd_file_path, idx) * 1e-3
            )
            pssn_list.append(pssn)

        # Calculate the GQ effectice PSSN
        gq_eff_pssn = self.opd_metr.calc_gq_value(pssn_list)

        return pssn_list, gq_eff_pssn

    def _calc_eff_fwhm_opd(self, pssn_list: list[float]) -> tuple[list[float], float]:
        """Calculate the effective FWHM of OPD.

        FWHM: Full width and half maximum.
        PSSN: Normalized point source sensitivity.
        GQ: Gaussian quadrature.

        Parameters
        ----------
        pssn_list : list
            List of PSSN.

        Returns
        -------
        list
            Effective FWHM list.
        float
            GQ effective FWHM.
        """

        # Calculate the list of effective FWHM
        eff_fwhm_list = []
        for pssn in pssn_list:
            eff_fwhm = self.opd_metr.calc_fwhm_eff(pssn)
            eff_fwhm_list.append(eff_fwhm)

        # Calculate the GQ effectice FWHM
        gq_eff_fwhm = self.opd_metr.calc_gq_value(eff_fwhm_list)

        return eff_fwhm_list, gq_eff_fwhm

    def map_opd_data_to_list_of_wf_err(
        self,
        opd_zk_file_name: str,
        sensor_id_list: list[int],
        sensor_name_list: list[str],
    ) -> list[SensorWavefrontError]:
        """Map the OPD data to the list of wavefront error.

        OPD: Optical path difference.

        Parameters
        ----------
        opd_zk_file_name : str
            OPD zk file name.
        sensor_id_list : list
            Reference sensor ID list.
        sensor_name_list : list
            Reference sensor name list.

        Returns
        -------
        list [lsst.ts.wep.ctrlIntf.SensorWavefrontError]
            List of SensorWavefrontError object.
        """

        opd_zk = get_zk_from_file(os.path.join(self.output_img_dir, opd_zk_file_name))

        list_of_wf_err = []
        for sensor_id, sensor_name in zip(sensor_id_list, sensor_name_list):
            sensor_wavefront_data = SensorWavefrontError(num_of_zk=self.num_of_zk)
            sensor_wavefront_data.sensor_id = sensor_id
            sensor_wavefront_data.sensor_name = sensor_name
            # imSim outputs OPD in nanometers so need to change to microns
            # to be consistent with what WEP expects.
            sensor_wavefront_data.annular_zernike_poly = opd_zk[sensor_id] / 1e3
            list_of_wf_err.append(sensor_wavefront_data)

        return list_of_wf_err

    def get_opd_gq_eff_fwhm_from_file(self, pssn_file_name: str) -> float:
        """Get the OPD GQ effective FWHM from file.

        OPD: Optical path difference.
        GQ: Gaussian quadrature.
        FWHM: Full width at half maximum.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssn_file_name : str
            PSSN file name.

        Returns
        -------
        float
            OPD GQ effective FWHM.
        """

        data = self._get_data_of_pssn_file(pssn_file_name)
        gq_eff_fwhm = data[1, -1]

        return gq_eff_fwhm

    def get_list_of_fwhm_sensor_data(
        self, pssn_file_name: str
    ) -> tuple[np.ndarray, np.ndarray]:
        """Get the list of FWHM sensor data based on the OPD PSSN file.

        FWHM: Full width at half maximum.
        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssn_file_name : str
            PSSN file name.
        sensor_id_list : list
            Reference sensor id list.

        Returns
        -------
        np.ndarray [object]
            Numpy array with fwhm data. This is a numpy array of arrays. The
            data type is `object` because each element may have different
            number of elements.
        """

        # Get the FWHM data from the PSSN file
        # The first row is the PSSN and the second one is the FWHM
        # The final element in each row is the GQ value
        data = self._get_data_of_pssn_file(pssn_file_name)
        fwhm_data = data[1, :-1]

        fwhm_collection = np.array([], dtype=object)
        for fwhm in fwhm_data:
            fwhm_collection = np.append(fwhm_collection, fwhm)

        return fwhm_collection

    def get_opd_pssn_from_file(self, pssn_file_name: str) -> np.ndarray:
        """Get the OPD PSSN from file.

        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssn_file_name : str
            PSSN file name.

        Returns
        -------
        numpy.ndarray
            PSSN.
        """

        data = self._get_data_of_pssn_file(pssn_file_name)
        pssn = data[0, :-1]

        return pssn

    def _get_data_of_pssn_file(self, pssn_file_name: str) -> np.ndarray:
        """Get the data of the PSSN file.

        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssn_file_name : str
            PSSN file name.

        Returns
        -------
        numpy.ndarray
            Data of the PSSN file.
        """

        file_path = os.path.join(self.output_img_dir, pssn_file_name)
        data = np.loadtxt(file_path)

        return data

    def reorder_and_save_wf_err_file(
        self,
        list_of_wf_err: list[SensorWavefrontError],
        ref_sensor_name_list: list[str],
        camera: cameraGeom.Camera,
        zk_file_name: str = "wfs.zer",
    ) -> None:
        """Reorder the wavefront error in the wavefront error list according to
        the reference sensor name list and save to a file.

        The unexisted wavefront error will be a numpy zero array. The unit is
        um.

        Parameters
        ----------
        list_of_wf_err : list [lsst.ts.wep.ctrlIntf.SensorWavefrontData]
            List of SensorWavefrontData object.
        ref_sensor_name_list : list
            Reference sensor name list.
        camera : lsst.afw.cameraGeom.Camera
            Instrument.
        zk_file_name : str, optional
            Wavefront error file name. (the default is "wfs.zer".)
        """

        # Get the sensor name that in the wavefront error map
        wf_err_map = self._trans_list_of_wf_err_to_map(list_of_wf_err, camera)
        name_list_in_wf_err_map = list(wf_err_map.keys())

        # Reorder the wavefront error map based on the reference sensor name
        # list.
        reordered_wf_err_map = dict()
        for sensor_name in ref_sensor_name_list:
            if sensor_name in name_list_in_wf_err_map:
                wf_err = wf_err_map[sensor_name]
            else:
                wf_err = np.zeros(self.num_of_zk)
            reordered_wf_err_map[sensor_name] = wf_err

        # Save the file
        file_path = os.path.join(self.output_img_dir, zk_file_name)
        file_txt = "# The followings are ZK in um from z4 to z22:\n"
        for key, val in reordered_wf_err_map.items():
            zk_str = f"{camera[key].getId()}: {val}\n"
            file_txt += zk_str
        with open(file_path, "w") as file:
            file.write(file_txt)

    def _trans_list_of_wf_err_to_map(
        self, list_of_wf_err: list[SensorWavefrontError], camera: cameraGeom.Camera
    ) -> dict[str, np.ndarray]:
        """Transform the list of wavefront error to map.

        Parameters
        ----------
        list_of_wf_err : list [lsst.ts.wep.ctrlIntf.SensorWavefrontData]
            List of SensorWavefrontData object.
        camera : lsst.afw.cameraGeom.Camera
            Instrument.

        Returns
        -------
        dict
            Calculated wavefront error. The dictionary key [str] is the
            abbreviated sensor name (e.g. R22_S11). The dictionary item
            [numpy.ndarray] is the averaged wavefront error (z4-z22) in um.
        """

        map_sensor_name_and_id = dict(
            [(detector.getId(), detector.getName()) for detector in camera]
        )

        wf_err_map = dict()
        for sensor_wf_data in list_of_wf_err:
            sensor_id = sensor_wf_data.sensor_id
            sensor_name = map_sensor_name_and_id[sensor_id]

            avg_err_in_um = sensor_wf_data.annular_zernike_poly

            wf_err_map[sensor_name] = avg_err_in_um

        return wf_err_map

    def _get_wf_err_values_and_stack_to_matrix(
        self, wf_err_map: dict[str, np.ndarray]
    ) -> np.ndarray:
        """Get the wavefront errors and stack them to be a matrix.

        Parameters
        ----------
        wf_err_map : dict
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

        value_matrix = np.empty((0, self.num_of_zk))
        for wf_err in wf_err_map.values():
            value_matrix = np.vstack((value_matrix, wf_err))

        return value_matrix

    def save_dof_in_um_file_for_next_iter(
        self, dof_in_um_file_name: str = "dofPertInNextIter.mat"
    ) -> None:
        """Save the DOF in um data to file for the next iteration.

        DOF: degree of freedom.

        Parameters
        ----------
        dof_in_um_file_name : str, optional
            File name to save the DOF in um. (the default is
            "dofPertInNextIter.mat".)
        """

        file_path = os.path.join(self.output_dir, dof_in_um_file_name)
        header = "The followings are the DOF in um:"
        np.savetxt(file_path, np.transpose(self.dof_in_um), header=header)
