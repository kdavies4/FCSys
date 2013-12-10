FCRes
-----
**Warning!**  This package is still in development.  There are likely to be many
bugs.  Please check back soon.

FCRes is a [Python] package to load, analyze, and plot simulation results from
[FCSys] - a [Modelica] library for modeling fuel cells and electrochemistry.  It
is based on and requires [ModelicaRes].  For more information about the
functions of FCRes, see the HTML or PDF documentation in the [doc folder](doc).

### Credits

Kevin Bandy helped to develop this package.

## Installation

First, download and set up [ModelicaRes].  It is hosted in the
[Python Package Index (PyPI)](http://pypi.python.org/pypi/ModelicaRes).  Follow
the installation instructions there. Next, install FCRes by running the set up
script ([setup.py](setup.py)) from the base folder.  On Windows, use the
following command:

    python setup.py install

On Linux, use:

    sudo python setup.py install

The [matplotlibrc](matplotlibrc) file in the base folder has some recommended
revisions to [matplotlib]'s defaults.  To use it, copy or move the file to the
working directory or [matplotlib]'s configuration directory.  See
http://matplotlib.org/users/customizing.html for details.

[FCSys]: http://kdavies4.github.io/FCSys/
[Python]: http://www.python.org
[Modelica]: http://www.modelica.org
[ModelicaRes]: http://kdavies4.github.io/ModelicaRes/
[matplotlib]: http://www.matplotlib.org
