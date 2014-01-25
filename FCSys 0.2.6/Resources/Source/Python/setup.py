#!/usr/bin/python
"""Set up the FCRes module.

See README.md for instructions.
"""
from distutils.core import setup
from glob import glob

setup(name='FCRes',
      version="0.2.6",
      author='Kevin Davies',
      author_email='kdavies4@gmail.com',
      #credits=['Kevin Bandy'],
      packages=['fcres'],
      scripts=glob('bin/*'),
      url='http://kdavies4.github.io/FCSys/',
      license = "Modelica License Version 2",
      description='Python utilities for the FCSys package',
      long_description=open('README.md').read(),
      package_data={'fcres': ['fcres.ini']},
      provides=['fcres'],
      requires=['modelicares  (>=0.7)'],
     )
