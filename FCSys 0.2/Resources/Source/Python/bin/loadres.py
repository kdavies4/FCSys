#!/usr/bin/env python
"""Load results from Modelica_ simulation(s) and provide a Python_ interpreter
to analyze the results.

This script can be executed at the command line.  It will accept as arguments
the names of result files or names of directories with result files.  The
filenames may contain wildcards.  If no arguments are given, the script
provides a dialog to choose a file or folder.  Finally, it provides working
session of `IPython <http://www.ipython.org/>`_ with the results preloaded.
PyLab_ is directly imported (``from pylab import *``) to provide many functions
of NumPy_ and matplotlib_ (e.g., :meth:`sin` and :meth:`plot`).  The essential
classes and functions of ModelicaRes are directly available as well.

**Example:**

.. code-block:: sh

   $ loadres examples
   Valid: SimRes('.../examples/ChuaCircuit.mat')
   Valid: SimRes('.../examples/ThreeTanks.mat')
   Valid: LinRes('.../examples/PID.mat')
   Simulation results have been loaded into sims[0] through sims[1].
   A linearization result has been loaded into lin.

where '...' depends on the local system.

You can now explore the simulation results or create plots using the methods in
:class:`modelicares.simres.SimRes`.  For example,

.. code-block:: python

   >>> sims[0].get_FV('L.v')
   -0.25352862
   >>> sims[0].get_unit('L.v')
   'V'

If a variable cannot be found, then suggestions are given:

.. code-block:: python

   >>> sims[0].get_description('L.vv')
   L.vv is not a valid variable name.
   <BLANKLINE>
   Did you mean one of the these?
          L.v
          L.p.v
          L.n.v
   >>> sims[0].get_description('L.v')
   'Voltage drop between the two pins (= p.v - n.v)'

To return all values of a variable, use its string as an index:

.. code-block:: python

   >>> sim['L.v']
   array([  0.00000000e+00, ... -2.53528625e-01], dtype=float32)

or an argument:

.. code-block:: python

   >>> sim('L.v')
   array([  0.00000000e+00, ... -2.53528625e-01], dtype=float32)

To see all the methods, use

   >>> help(sims[0])

To search for variables, use :meth:`modelicares.simres.SimRes.glob` with
wildcards:

   >>> sims[0].glob('L.p*')
   [u'L.p.i', u'L.p.v']

Likewise, you can explore the linearization result or create diagrams using the
methods in :class:`modelicares.linres.LinRes`:

.. code-block:: python

   >>> print lin
   Modelica linearization results from ".../examples/PID.mat"
   >>> lin.sys.A
   matrix([[   0.,    0.],
           [   0., -100.]])

.. _Modelica: http://www.modelica.org/
.. _Python: http://www.python.org/
.. _PyLab: http://www.scipy.org/PyLab
.. _NumPy: http://numpy.scipy.org/
.. _matplotlib: http://www.matplotlib.org/
"""
from fcres import CellSimRes, CellLinRes, multiload
from sys import argv
from time import sleep
from modelicares.gui import boolbox
#from wx import App, DirSelector, FileSelector
import os
from easygui import fileopenbox, diropenbox

if __name__ == '__main__':
    default_path = '../examples'

    def _local_exit(t=0.5):
        """Exit with a message and a delay."""
        print("Exiting...")
        sleep(t)
        exit()

    #choice = MessageDialog()

    # Determine if one file or multiple files should be loaded.
    if len(argv) == 1:
        choice = boolbox("Open one file or all files from a folder?",
                         title="Choice", choices=('File', 'Folder'), default=0)
        if choice is None:
            _local_exit()
        multi = choice == 1
    elif len(argv) == 2:
        if argv[1] in ['-test', '--test']:
            import doctest
            doctest.testmod()
            exit()
        else:
            multi = False
    else:
        multi = True

    # Load the file(s).
    if multi:
        if len(argv) == 1: # Prompt for a directory.
            loc = diropenbox("Choose the folder with the data file(s).",
                             default=default_path)
            # easygui is ugly, but wx seems to make the working session slow:
            #app = App()
            #loc = DirSelector("Choose the folder with the data file(s).",
            #    defaultPath=default_path)
            if loc == '': _local_exit()
        else:
            loc = argv[1:]
        sims, lins = multiload(loc)
        n_sim = len(sims)
        n_lin = len(lins)
        if n_sim == 0:
            if n_lin == 0:
                print("No files were loaded.")
                _local_exit()
            else:
                print("No cell simulation results have been loaded.")#, end="")
                      # TODO: Re-introduce the end="" for this and below when
                      # running Python 3 (to suppress the line ending).
                del sims
        elif n_sim == 1:
            print("A cell simulation result has been loaded into sims[0].")
        else:
            print("Cell simulation results have been loaded into sims[0] "
                  "through sims[%i]." % (n_sim-1))#, end="")
        if n_lin == 0:
            print("No cell linearization results have been loaded.")#, end="")
            del lins
        elif n_lin == 1:
            print("A cell linearization result has been loaded into lins[0].")
                  #, end="")
        else:
            print("Cell linearization results have been loaded into lins[0] "
                  "through lins[%i]." % (i-1))#, end="")
    else:
        if len(argv) == 1: # Prompt for file.
            fname = fileopenbox(msg="Select a data file.",
                default=os.path.join(default_path, "*.mat"),
                filetypes=['*.mat'])
            # easygui is ugly, but wx seems to make the working session slow:
            #app = App()
            #fname = FileSelector("Select a data file.",
            #    default_path=default_path, wildcard='*.mat')
            if fname == '': _local_exit()
        else:
            fname = argv[1]
        try:
            sim = CellSimRes(fname)
            print('Cell simulation results have been loaded from "%s".\n'
                  'The CellSimRes instance is sim.' % fname)
        except IOError:
            _local_exit()
        except:
            try:
                lin = LinRes(fname)
                print('Cell linearization results have been loaded from
                      '"%s".\nThe CellLinRes instance is lin.' % fname)
            except:
                print('The file, "%s", could not be opened as a cell '
                      'simulation or linearization result.' % fname)
                _local_exit()

    # Open the IPython or standard Python interpreter.
    #    http://writeonly.wordpress.com/2008/09/08/embedding-a-python-shell-in-a-python-script/,
    #    accessed 11/2/2010
    from pylab import *
    try:
        from IPython.Shell import IPShellEmbed
        IPShellEmbed(argv=['-noconfirm_exit'])()
        # Note: The -pylab option cannot be embedded (see
        # http://article.gmane.org/gmane.comp.python.ipython.user/1190/match=pylab)
    except ImportError:
        from code import InteractiveConsole
        # Calling this with globals ensures that we can see the environment.
        InteractiveConsole(globals()).interact()
