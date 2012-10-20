#!/usr/bin/python
"""Set up the FCRes module.

Instructions for installation on Linux:

1. Build the modules.
   $ ./setup.py build

2. Install the modules.
   $ sudo ./setup.py install
"""
from distutils.core import setup

# Install the core Python modules.
setup(name='FCRes',
      version='1.0',
      author='Kevin Davies',
      author_email='kdavies4@gmail.com',
      packages=['fcres'],
      scripts=['bin/loadres.py'],
      url='http://kdavies4.github.com/FCSys/',
      license = "Modelica License Version 2",
      description='Python utilities for the FCSys package',
      requires=['modelicares'],
      keywords = "modelica plot dymola mat fuel cell FC PEMFC FCSys",
      classifiers = ["Development Status :: 4 - Beta",
                     "Environment :: Console",
                     "License :: OSI Approved :: BSD License",
                     "Programming Language :: Python :: 2.7",
                     "Intended Audience :: Science/Research",
                     "Topic :: Scientific/Engineering",
                     "Topic :: Utilities",
                     ],
     )
# Note: modelicares is available at http://kdavies4.github.com/ModelicaRes/.
