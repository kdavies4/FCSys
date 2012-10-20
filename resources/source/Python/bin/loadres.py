#!/usr/bin/env python
"""Load results from FCSys simulation(s) and provide a Python_ interpreter to
analyze the results.

This file can be executed at the command line.  It will accept as an argument
the name of a results file or a directory with multiple files.  If no arguments
are provided, it gives a dialog to choose the file or folder. Finally, it
provides a working session of `IPython <http://www.ipython.org/>`_ with those
results preloaded.

**Example:**

   .. code-block:: sh

      $ cd examples
      $ loadres.py Polarization.mat
      Simulation results have been loaded from "Polarization.mat".
      The SimRes instance is sim.
      In [1]:

.. _Modelica: http://www.modelica.org/
.. _Python: http://www.python.org/
"""
from fcres import CellSimRes, CellLinRes, multi_load
from pylab import *
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
            # easygui is ugly, but wx seems make the working session slow:
            #app = App()
            #loc = DirSelector("Choose the folder with the data file(s).",
            #    defaultPath=default_path)
            if loc == '': _local_exit()
        else:
            loc = argv[1:]
        sims, lins = multi_load(loc)
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
            # easygui is ugly, but wx seems to slow down the working session:
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
    try:
        from IPython.Shell import IPShellEmbed
        IPShellEmbed(argv=['-noconfirm_exit'])()
        # Note: The -pylab option cannot be embedded (see
        # http://article.gmane.org/gmane.comp.python.ipython.user/1190/match=pylab)
    except ImportError:
        from code import InteractiveConsole
        # Calling this with globals ensures that we can see the environment.
        InteractiveConsole(globals()).interact()
