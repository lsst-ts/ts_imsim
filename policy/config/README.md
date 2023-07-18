# Config directory

This directory contains default settings to run the closed loop with ImSim.
The way the `ts_imsim` closed loop runs is to create an ImSim configuration file based upon the individual yaml files for each submodule of ImSim as they are organized in this directory.
Then a "pointer file" is used to tell `ts_imsim` which versions of each of these input files to use when running the closed loop (see [lsstCamDefaultPointer.yaml](lsstCamDefaultPointer.yaml) for an example pointer file).
For certain configuration settings there are also command line arguments to override the settings in the files (run `imgClosedLoop.py -h` on command line for all options).

## General Configuration Files
* **lsstCamDefaultPointer.yaml**: This file is a default pointer file to run the closed loop on the LSST Camera Corner Wavefront Sensors.
* **obsVariablesDefault.yaml**: Default observation information settings for ImSim configuration.
* **requiredModulesDefault.yaml**: List of modules ImSim requires to run the closed loop.

## ImSim Submodule Configuration
Inside the **config** folder are subfolders for each piece of the `ts_imsim` configuration file.
The pointer file required by `ts_imsim` indicates which of these files to use for a given closed loop run.
