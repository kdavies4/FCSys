#########
  FCRes
#########

**Warning!** This package is still in development. There are likely to
be many bugs. Please check back soon.

FCRes is a `Python <http://www.python.org>`_ package to load, analyze, and plot
simulation results from `FCSys <http://kdavies4.github.io/FCSys/>`_ - a
`Modelica <http://www.modelica.org>`_ library for modeling fuel cells and
electrochemistry. It is based on and requires ModelicaRes_.  For more
information about the functions of FCRes, see the HTML or PDF documentation in
the doc folder.

Credits
~~~~~~~

Kevin Bandy helped to develop this package.

Installation
~~~~~~~~~~~~

First, download and set up ModelicaRes_.  It is hosted in the `Python Package
Index (PyPI) <http://pypi.python.org/pypi/ModelicaRes>`_.  Follow the
installation instructions there.  Next, install FCRes by running the set-up
script (setup.py) from the base folder.  On Windows, use the following command::

    python setup.py install

On Linux, use::

    sudo python setup.py install

The matplotlibrc file in the base folder has some recommended revisions to
matplotlib_'s defaults. To use it, copy or move the file to the working
directory or matplotlib_'s configuration directory. See
http://matplotlib.org/users/customizing.html for details.


.. _ModelicaRes: http://kdavies4.github.io/ModelicaRes/
.. _matplotlib: http://www.matplotlib.org
