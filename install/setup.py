#!/usr/bin/env python
"""Set up the FCSys package.

Instructions for installation on Linux:

1. Change to the directory with this file.
   $ cd FCSys/install

2. Build the modules.
   $ python ./setup.py build

3. Install the modules.
   $ sudo python ./setup.py install
"""

from distutils.core import setup

# Install the Python module to analyze and plot results from FCSys.
setup(name='fcsys',
      version='1.0',
      description='Python utilities for FCSys package',
      author='Kevin Davies',
      author_email='kld@alumni.carnegiemellon.edu',
      license = "Modelica License Version 2",
      url='http://modelica.org/libraries', # **Give full URL once available.
      requires=['scipy', 'matplotlib', 'modelicares'],
      package_dir = {'': '../resources/source/Python'},
      py_modules=['fcres'],
     )
# Note: modelicares, res, easywx, texunit, and arrow_line are all a part of the
# modelicares utilities, available at http://kdavies4.github.com/modelicares/.

# In order to run dymosim from the command line in Linux, add this line to
# /etc/environment or ~/.pam_environment:
# LD_LIBRARY_PATH=/opt/dymola/bin/lib
