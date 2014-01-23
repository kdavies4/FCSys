#########
  FCRes
#########

**Python utilities to analyze and plot results from FCSys**

FCRes is an open-source tool to interpret the results of simulations from the
FCSys_ package and create publishable figures.  It is currently possible to

- Auto-generate simulation scripts,
- Run model executables with varying parameters,
- Browse data,
- Perform custom calculations, and
- Produce various plots and diagrams.

The figures are generated via matplotlib_, which offers a rich set of plotting
routines.  FCRes includes convenient functions to automatically pre-format
and label some figures, like xy plots, Bode and Nyquist plots, and Sankey
diagrams.  FCRes can be scripted or run from a Python_ interpreter with math and
matrix functions from NumPy_.

.. only:: html

  The links below describe the components of FCRes.  The top-level module,
  :mod:`fcres`, provides direct access to the most important classes and
  functions.  Others must be accessed through their submodules.  The
  :mod:`loadres` script loads data files and provides a Python_ interpreter to
  help analyze them.  The :mod:`fcres.fcsimres` submodule has classes to load,
  analyze, and plot results of fuel cell simulations.  The :mod:`fcres.fclinres`
  submodule has a class to load, analyze, and plot results from linearizing a
  fuel cell model.  The :mod:`fcres.simres` and :mod:`fcres.linres` submodules
  are similar to :mod:`fcres.simres` and :mod:`fcres.linres`, but are not
  specific to fuel cell simulations.  The :mod:`fcres.multi` submodule has
  functions to load and plot results from multiple data files at once.  The
  :mod:`fcres.exps` submodule has tools to set up and manage simulation
  experiments.  The :mod:`fcres.texunit` submodule has functions to translate
  Modelica_ *unit* and *displayUnit* strings into LaTeX_-formatted strings. The
  last submodule, :mod:`fcres.base`, has general supporting functions.

.. only:: latex

  The following chapters describe the components of FCRes.  The top-level
  module, :mod:`fcres`, provides direct access to the most important classes and
  functions.  Others must be accessed through their submodules.  The
  :mod:`loadres` script loads data files and provides a Python_ interpreter to
  help analyze them.  The :mod:`fcres.fcsimres` submodule has classes to load,
  analyze, and plot results of fuel cell simulations.  The :mod:`fcres.fclinres`
  submodule has a class to load, analyze, and plot results from linearizing a
  fuel cell model.  The :mod:`fcres.simres` and :mod:`fcres.linres` submodules
  are similar to :mod:`fcres.simres` and :mod:`fcres.linres`, but are not
  specific to fuel cell simulations.  The :mod:`fcres.multi` submodule has
  functions to load and plot results from multiple data files at once.  The
  :mod:`fcres.exps` submodule has tools to set up and manage simulation
  experiments.  The :mod:`fcres.texunit` submodule has functions to translate
  Modelica_ *unit* and *displayUnit* strings into LaTeX_-formatted strings. The
  last submodule, :mod:`fcres.base`, has general supporting functions.

  **Installation**

  First, download and install ModelicaRes_.

  Then install this package.  It can be downloaded as a part of FCSys_ or
  individually from the `PyPI page`_.  To install the package, download and
  extract it.  Then run the set up script (setup.py) from the base folder.  On
  Windows, use the following command::

     python setup.py install

  On Linux, use::

     sudo python setup.py install

  The matplotlibrc file in the base folder has some recommended revisions to
  matplotlib_'s defaults.  To use it, copy or move the file to the working
  directory or matplotlib_'s configuration directory.  See
  http://matplotlib.org/users/customizing.html for details.

  **Credits**

  Kevin Bandy helped in the development of this package.

  **License terms**

  FCRes_ is published under the terms of the BSD license (see LICENSE.txt).
  Please share any modifications you make (preferably on a Github fork from
  https://github.com/kdavies4/FCRes) so that others may benefit from your work.

  **See also**

  This module is based on `ModelicaRes
  <http://kdavies4.github.io/ModelicaRes/>`_.  The `pysimulator
  <https://code.google.com/p/pysimulator/>`_, `BuildingsPy
  <http://simulationresearch.lbl.gov/modelica/buildingspy/>`_, `DyMat
  <http://www.j-raedler.de/projects/dymat/>`_, and
  `awesim <https://github.com/saroele/awesim>`_ projects provide related Python_
  modules.  pysimulator_ includes an elaborate GUI and supports the Functional
  Model Interface (FMI).  BuildingsPy_ has a :class:`Tester` class that can be
  used for unit testing.  DyMat_ has functions to export Modelica_ simulation
  data to comma separated values (CSV), `Gnuplot <http://www.gnuplot.info/>`_,
  MATLAB:sup:`®`, and `Network Common Data Form (netCDF)
  <http://www.unidata.ucar.edu/software/netcdf/>`_.  awesim_ provides tools to
  help run simulation experiments and organize the results.

.. toctree::

  fcres
..  fcsimres
  fclinres
  simres
  linres
  multi
  exps
  texunit
  base

.. only:: html

  A PDF version of this documentation is available `here <FCRes.pdf>`_.

  **Installation**

  First, download and install ModelicaRes_.

  Then install this package.  It can be downloaded as a part of FCSys_ or
  individually from the `PyPI page`_.  To install the package, download and
  extract it.  Then run the set up script (setup.py) from the base folder.  On
  Windows, use the following command::

     python setup.py install

  On Linux, use::

     sudo python setup.py install

  The matplotlibrc file in the base folder has some recommended revisions to
  matplotlib_'s defaults.  To use it, copy or move the file to the working
  directory or matplotlib_'s configuration directory.  See
  http://matplotlib.org/users/customizing.html for details.

  **Credits**

  The main author is Kevin Davies.  Kevin Bandy also helped in the development.

  **License terms**

  FCRes_ is published under the terms of the BSD license (see LICENSE.txt).
  Please share any modifications you make (preferably on a Github fork from
  https://github.com/kdavies4/FCRes) so that others may benefit from your work.

  **See also**

  The `pysimulator`_, `BuildingsPy`_, DyMat_, and `awesim`_ projects provide
  related Python_ modules.  pysimulator_ includes an elaborate GUI and supports
  the Functional Model Interface (FMI).  BuildingsPy_ has a :class:`Tester`
  class that can be used for unit testing.  DyMat_ has functions to export
  Modelica_ simulation data to comma separated values (CSV), `Gnuplot`_,
  MATLAB:sup:`®`, and `Network Common Data Form (netCDF)`_.  awesim_ provides
  tools to help run simulation experiments and organize the results.


.. _Modelica: http://www.modelica.org/
.. _Python: http://www.python.org/
.. _matplotlib: http://www.matplotlib.org/
.. _NumPy: http://numpy.scipy.org/
.. _FCSys: http://kdavies4.github.io/FCSys/
.. _PyPI page: http://pypi.python.org/pypi/FCRes
.. _ModelicaRes: http://kdavies4.github.io/ModelicaRes/
.. _LaTeX: http://www.latex-project.org/
