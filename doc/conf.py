"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documentation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.cp.verify


_g = globals()
_g.update(build_package_configs(
    project_name='lsst.cp.verify',
    version=lsst.cp.verify.version.__version__))
