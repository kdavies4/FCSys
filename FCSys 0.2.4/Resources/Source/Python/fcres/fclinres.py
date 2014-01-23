#!/usr/bin/env python
"""Load and analyze results from FCSys_.

**update (moved to other script): The "fcres.py" file can be executed at the command line.  It will accept as an
argument the name of a results file or a directory with multiple files.  If no
arguments are provided, then it gives a dialog to choose the file or folder.
Finally, it provides a working session of `Python <http://www.python.org/>`_
with those results preloaded.

**Example:**

   .. code-block:: sh

      $ ./fcres.py examples/Polarization.mat
      Cell simulation results have been loaded from "examples/Polarization.mat".
      The CellSimRes instance is sim.
      In [1]:

.. _FCSys: http://kdavies4.github.com/FCSys/
.. _Modelica: http://www.modelica.org/
"""
__author__ = "Kevin Davies"
__email__ = "kdavies4@gmail.com"
__copyright__ = "Copyright 2012-2013, Georgia Tech Research Corporation"
__license__ = "BSD-compatible (see LICENSE.txt)"

from modelicares import LinRes

class FCLinRes(LinRes):
    """Fuel cell linearization results from FCSys_ and methods to analyze those
    results
    """
    pass # TODO: Add customizations here.

if __name__ == '__main__':
    """Test the contents of this file."""
    import doctest
    doctest.testmod()
    exit()
