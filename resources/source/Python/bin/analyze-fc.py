#!/usr/bin/env python
"""Produce basic plots of a fuel cell simulation.

Create a shortcut or link to this file from a folder that contains fuel cell
simulation results in a file named "Cell.mat".
"""
__author__ = "Kevin Davies"
__email__ = "kld@alumni.carnegiemellon.edu"
__version__ = "2011/10/29"

import matplotlib.pyplot as plt
import numpy as np

from fcres import CellSimRes

# Chemical components
SPECIES = ['eminus', 'Hplus', 'H2', 'H2O', 'N2', 'O2']
#SPECIES = ['eminus']
TEX = dict(eminus='$e^-$', Hplus='$H^+$', H2='$H_2$', H2O='$H_2O$', N2='$N_2$',
           O2='$O_2$')

# Load the results.
sim = CellSimRes("Cell.mat")

# Create plots specific to a single species
for species in SPECIES:
    #sim.plotfig_subregions(species + '.p',
    #                       title="Pressure of %s within Subregions over Time"
    #                             %TEX[species])
    #sim.plotfig_subregions(species + '.T',
    #                       title="Temperature of %s within Subregions over Time"
    #                             %TEX[species])
    sim.quiverfig_subregions(vect=species + '.center.Phi[%i].lower', n_rows=2,
                             times=np.linspace(0, sim.get_FV('Time'), 4),
                             title="Velocity Field of %s"%TEX[species],
                             xlabel="")

# Plot current distribution
sim.currdenfig()

plt.show()
