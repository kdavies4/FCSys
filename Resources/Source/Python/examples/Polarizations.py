#!/usr/bin/env python
"""This a demo script to set up and run cell or general simulations in Dymola,
and then create plots of the results.
"""
__author__ = "Kevin Davies"
__version__ = "2012-10-11"

import os
import matplotlib.pyplot as plt

from res import save_figs
from modelicares import SimRes, SimSpec, gen_sim_script
from fcres import CellSimRes

# Settings
# Run the simulations? (otherwise, just plot)
run = False
# Name of the Modelica script (including the path)
fname = os.path.abspath('run-sims.mos')
# Directory where the working files (e.g., dslog.txt) should be placed
working_dir = os.path.abspath('../'*5) # The root of the FCSys folder is 4 levels above this one.
# List of Modelica packages that should be loaded (besides the Modelica Standard
# Library)
# Each may be a *.mo file or a path where a package.mo file resides, e.g.,
# "/opt/dymola/Modelica/Library/VehicleInterfaces 1.1.1".
packages = [os.path.abspath('../'*4), # FCSys
           ]
# List of simulations to run
simspecs = [SimSpec('FCSys.Subassemblies.Cells.Examples.Polarization', # Name of model
                startTime=None, # Start of simulation (default: 0)
                stopTime=1e-22, # End of simulation (default: 1)
                numberOfIntervals=None, # Number of output points (default: 0)
                outputInterval=None, # Distance between output points (default: 0)
                method=None, # Integration method (default: "Dassl")
                tolerance=1e-6, # Tolerance of integration (default: 0.0001)
                fixedstepsize=None, # Fixed step size for Euler (default: 0)
                resultFile='Polarization') # Where to store result
        ]
# Formats in which to save the figures, e.g., ['pdf', 'eps', 'svg', 'png']
# If the figures shouldn't be saved, specify an empty list.
formats = ['pdf', 'png']

if run:
    # Create the script to load the packages, simulate, and save the results.
    models, results_dir = gen_sim_script(simspecs, packages=packages,
                                         fname=fname, working_dir=working_dir)

    # Ask Dymola to run the script.
    os.system("bash /opt/dymola/bin/dymola.sh " + fname) # This is for Linux.
    # TODO: Add support for Windows.
else:
    models = [simspec.problem[simspec.problem.rfind('.')+1:] for simspec in simspecs]
    results_dir = os.path.split(fname)[0]

# Create plots.
# Note: The code between the '---' lines must be customized for each simulation
# and plot.
# ------------------------------

## Model 1
label = simspecs[0].resultFile
sim = CellSimRes(os.path.join(results_dir, label))
sim.plotfig(xname='cell.I',
            ynames1='cell.v', ylabel1="Potential", legends1="Average voltage",
            ynames2='cell.Wdot', ylabel2="Power", legends2="Power output",
            title="Cell Polarization", label=label)

# ------------------------------

# Save the plots.
save_figs(formats)
plt.show()
