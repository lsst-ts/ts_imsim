.. py:currentmodule:: lsst.ts.imsim

.. _lsst.ts.imsim-version_history:

##################
Version History
##################

-------------
1.2.0
-------------

* Update ts_wep API calls to run with ts_wep v9.0.
* Add some methods from old version of ts_wep to preserve functionality.

-------------
1.1.1
-------------

* Update wavefront sensor locations in tests with updated info from obs_lsst.

-------------
1.1.0
-------------

* Add alt, az calculation to obs_metadata.
* Change default imsim configuration to use LSST_Silicon stamp type and include RubinDiffractionOptics.
* Add imsim_log_file option to capture imsim output in a log file instead of printing to STDOUT.
* Add link to docs on README.md.

-------------
1.0.1
-------------

* Add PROGRAM to default closed loop headers so that science_program is populated when ingesting raws into butler repository.

-------------
1.0.0
-------------

* First Release Version

-------------
0.10.2
-------------

* Replace RAWSEEING in headers with SEEING to avoid FITS format error messages.
* Fix versionHistory header.

-------------
0.10.1
-------------

* Replace all inst_name calls with CamType.

-------------
0.10.0
-------------

* Add comcam support.
* Add new comcam configuration files.

-------------
0.9.0
-------------

* Remove ZCS coordinate system and transition to CCS.
* Update bending modes to latest versions in batoid_rubin.

-------------
0.8.4
-------------

* Update default MJD in obs_metadata to be consistent with default in closed_loop_task.

-------------
0.8.3
-------------

* Update default MJD to avoid twilight.

-------------
0.8.2
-------------

* Add .ruff.toml to .gitignore.

-------------
0.8.1
-------------

* Add correct zenith angle into imsim configuration files.
* Add parallactic angle and zenith angle calculation into obsMetadata.

-------------
0.8.0
-------------

* Add star_mag to img_closed_loop to specify default star magnitude.

-------------
0.7.0
-------------

* Add T&S pre-commit settings to ts_imsim.
* Change file names to snake_case.
* Move imgClosedLoop to img_closed_loop.

-------------
0.6.3
-------------

* Remove rawSeeing from stamp configuration yaml.
* Patching to fix compatibility with ts_ofc double Zernikes update.
* Rotation in closed loop now available and converging.

-------------
0.6.2
-------------

* Change bending modes to legacy bending modes until ts_ofc is updated.

-------------
0.6.1
-------------

* Add documentation with user and developer guides.

-------------
0.6.0
-------------

* Change from camelCase to snake_case.
* Add typing.
* Rename opdOnly to turn_off_wavefront_estimates.

-------------
0.5.4
-------------

* Adding 180 degree rotation in rotationMatrix to account for photons farthest from Zenith on sky appear on "top".

.. _lsst.ts.imsim-0.5.3:

-------------
0.5.3
-------------

* Fix rotation sign and interpolation approach when rotating opd.

.. _lsst.ts.imsim-0.5.2:

-------------
0.5.2
-------------

* Adding seeing as parameter for simulations.

.. _lsst.ts.imsim-0.5.1:

-------------
0.5.1
-------------

* Add MacOS support.

.. _lsst.ts.imsim-0.5.0:

-------------
0.5.0
-------------

* Add FAM support.
* Debug rotation problems.

.. _lsst.ts.imsim-0.4.2:

-------------
0.4.2
-------------

* Add config files for testing convergence with and without perturbations and fam testing files.

.. _lsst.ts.imsim-0.4.1:

-------------
0.4.1
-------------

* Update to use ts_wep v7.0.

.. _lsst.ts.imsim-0.4.0:

-------------
0.4.0
-------------

* Add closed loop OPD only mode.

.. _lsst.ts.imsim-0.3.0:

-------------
0.3.0
-------------

* Add closed loop infrastructure.
* Update README.
* Update Jenkinsfile to work with latest Jenkins environment changes.

.. _lsst.ts.imsim-0.2.0:

-------------
0.2.0
-------------

* Add configuration file creation for ImSim image generation.
* Update Jenkinsfile to run correctly.
* Add documentation stub to get Jenkins status checks to pass in github.

.. _lsst.ts.imsim-0.1.0:

-------------
0.1.0
-------------

* Initial stub of imsim repository.
