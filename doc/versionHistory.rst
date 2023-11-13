.. py:currentmodule:: lsst.ts.wep

.. _lsst.ts.wep-version_history:

##################
Version History
##################

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
