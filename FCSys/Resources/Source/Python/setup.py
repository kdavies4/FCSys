#!/usr/bin/python
"""Set up the FCRes module.

See README.md for instructions.
"""
from distutils.core import setup
from glob import glob

import fcres # Only to read the version number

setup(name='FCRes',
      version=fcres.__version__,
      author='Kevin Davies',
      author_email='kdavies4@gmail.com',
      #credits=['Kevin Bandy'],
      packages=['fcres'],
      scripts=glob('bin/*'),
      url='http://kdavies4.github.com/FCSys/',
      license = "Modelica License Version 2",
      description='Python utilities for the FCSys package',
      long_description=open('README.md').read(),
      package_data={'fcres': ['fcres.ini']},
      provides=['fcres'],
      requires=['modelicares  (>=0.5)'],
      keywords=['Modelica', 'Dymola', 'matplotlib', 'proton exchange membrane',
                'polymer exchange membrane', 'fuel cell', 'PEMFC', 'FCSys',
                'electrochemistry'],
      classifiers=['Development Status :: 3 - Alpha',
                   'Operating System :: POSIX :: Linux',
                   'Operating System :: Microsoft :: Windows',
                   'Environment :: Console',
                   'License :: OSI Approved :: BSD License',
                   'Programming Language :: Python :: 2.7',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Utilities',
                     ],
     )
