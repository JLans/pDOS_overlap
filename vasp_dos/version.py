from __future__ import absolute_import, division, print_function
from os.path import join as pjoin

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 1
_version_micro = ''  # use '' for first of series, number for 1 and above
_version_extra = 'dev'
_version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# Description should be a one-liner:
description = "vasp_dos: A module for extracting and plotting vasp density of states information"
# Long description will go up on the pypi page
long_description = """

vasp_dos
========
This module extracts density of states information.


License
=======
``vasp_dos`` is licensed under the terms of the MIT license. See the file
"LICENSE" for information on the history of this software, terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.

All trademarks referenced herein are property of their respective holders.

Copyright (c) 2020--, Joshua Lansford, The University of Delaware.
"""

NAME = "vasp_dos"
MAINTAINER = "Joshua Lansford"
MAINTAINER_EMAIL = "lansford.jl@gmail.com"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "http://github.com/JLANS/vasp_dos"
DOWNLOAD_URL = ""
LICENSE = "MIT"
AUTHOR = "Joshua Lansford"
AUTHOR_EMAIL = "lansford.jl@gmail.com"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGE_DATA = {'vasp_dos': [pjoin('data', '*'),pjoin('data', '*','*'),pjoin('data', '*','*','*')]}
REQUIRES = ["numpy", "matplotlib", "ase"]
PYTHON_REQUIRES = ">= 3.5"
