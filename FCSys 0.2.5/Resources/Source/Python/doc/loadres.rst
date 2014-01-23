:mod:`loadres`
==============

Load results from FCSys simulation(s) and provide a Python_ interpreter to
analyze the results.

This file can be executed at the command line.  It will accept as an argument
the name of a results file or a directory with multiple files.  If no arguments
are provided, it gives a dialog to choose the file or folder.  Finally, it
provides a working session of `IPython <http://www.ipython.org/>`_ with those
results preloaded.  pylab is directly imported (``from pylab import *``) to
provide many functions of numpy_ and matplotlib_ (e.g., :meth:`sin` and
:meth:`plot`)

**Example:**

   .. code-block:: sh

      $ cd examples
      $ loadres.py Polarization.mat
      Simulation results have been loaded from "Polarization.mat".
      The SimRes instance is sim.
      In [1]:

.. _Modelica: http://www.modelica.org/
.. _Python: http://www.python.org/
