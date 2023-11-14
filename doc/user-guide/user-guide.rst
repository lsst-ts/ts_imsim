.. _User_Guide:

#####################
User Guide
#####################

The `ts_imsim` closed loop is designed to run as an executable through the command line as `img_closed_loop.py`.

.. _Installing_ts_imsim:

Installing `ts_imsim`
=====================

Needed Packages
---------------

- `ts_wep <https://github.com/lsst-ts/ts_wep>`_
- `ts_ofc <https://github.com/lsst-ts/ts_ofc>`_
- `imSim <https://github.com/LSSTDESC/imSim>`_
- `astroplan <https://github.com/astropy/astroplan>`_
- `black <https://github.com/psf/black>`_ (optional)
- `documenteer <https://github.com/lsst-sqre/documenteer>`_ (optional)
- `plantuml <http://plantuml.com>`_ (optional)
- `sphinxcontrib-plantuml <https://pypi.org/project/sphinxcontrib-plantuml/>`_ (optional)

Install ts_imsim
-----------------------------------------

1. Install and set up required packages from list above using `eups`.

2. Clone the `ts_imsim` package.

3. Setup the imSim environment by `eups`:

.. code-block:: bash

    cd $ts_imsim_directory
    setup -k -r .
    scons

.. _Using ts_imsim:

Using `ts_imsim`
================

Running the AOS Closed Loop
---------------------------

The closed loop is available from the command line after running `scons` with the command `img_closed_loop.py`.

A sample command to run 5 iterations of the closed loop with the lsst wavefront sensors is:

.. code-block:: python

    img_closed_loop.py --output cwfs_imsim_output/ --iter_num 5 --config_pointer_file ts_imsim/policy/config/lsstCamNoPertPointer.yaml --turn_off_atmosphere --turn_off_sky_background

where the following minimal set of command line options are used:

* **output**: Directory path for output.
  If it does not exist it will be made.
* **iter_num**: Number of iterations for the closed loop.
* **config_pointer_file**: The pointer file that contains the paths to all the individual configuration files for the simulation subsystems (input.telescope, output, stamp, etc.).
  See next section for more information.
* **turn_off_atmosphere**: Turn off the atmosphere in the images.
* **turn_off_sky_background**: Turn off the sky background. Helps to speed up image generation at expense of a realistic background.

Changing the Initial Telescope Configuration
--------------------------------------------

Module Configuration Files
**************************

`ts_imsim` is designed to take advantage of the modular nature of imSim by using configuration files linked together by a pointer file.
In the `policy/config` directory of `ts_imsim` there are folders for each of the following submodules with a default example configuration in each folder:

* **gal**: Profiles of point source and galaxy objects.
* **image**: Handles building an image with information such as bandpass, noise properties, sensor type.
* **input**: Models describing the telescope and atmosphere.
* **opd**: Optional. Provides information needed to create an output file showing the optical path difference at chosen locations in the focal plane.
* **output**: Define output files and image readout properties.
* **psf**: Describe the point spread function (PSF) model for the image.
* **stamp**: Describe the construction of a stamp for a single object.

Pointer Files
*************

To specify the version of each module configuration you wish to use one must point to each configuration file in a general pointer file written in yaml.
An example of what this looks like is provided at `policy/config/lsstCamDefaultPointer.yaml` and reproduced below:

.. code-block:: yaml

    # Default LSSTCam configuration files
    input:
      atm_psf: '{TS_IMSIM_DIR}/policy/config/input/atm_psf/atmPsfDefault.yaml'
      sky_model: '{TS_IMSIM_DIR}/policy/config/input/sky_model/skyModelDefault.yaml'
      telescope: '{TS_IMSIM_DIR}/policy/config/input/telescope/lsstCamDefault.yaml'
      vignetting: '{TS_IMSIM_DIR}/policy/config/input/vignetting/lsstCamDefault.yaml'
    gal: '{TS_IMSIM_DIR}/policy/config/gal/starDefault.yaml'
    image: '{TS_IMSIM_DIR}/policy/config/image/lsstCamDefault.yaml'
    psf: '{TS_IMSIM_DIR}/policy/config/psf/psfDefault.yaml'
    stamp: '{TS_IMSIM_DIR}/policy/config/stamp/lsstCamDefault.yaml'
    output: '{TS_IMSIM_DIR}/policy/config/output/lsstCamDefault.yaml'
    opd: '{TS_IMSIM_DIR}/policy/config/opd/lsstCamDefault.yaml'

The `ts_imsim` code understands the environmental variable `TS_IMSIM_DIR` after running scons and that can be used in these pointer files.
