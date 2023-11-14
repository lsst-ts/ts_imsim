"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documentation builds.
"""

import os

import lsst.ts.imsim  # noqa
from documenteer.conf.pipelinespkg import *  # noqa

project = "ts_imsim"
html_theme_options["logotext"] = project  # noqa
html_title = project
html_short_title = project
doxylink = {}


# Support the sphinx extension of plantuml
extensions.append("sphinxcontrib.plantuml")  # noqa

# Put the path to plantuml.jar
plantuml_path = (
    os.path.expanduser("~/plantuml.jar")
    if os.getenv("PATH_PLANTUML") is None
    else os.getenv("PATH_PLANTUML")
)
plantuml = f"java -jar {plantuml_path}"
