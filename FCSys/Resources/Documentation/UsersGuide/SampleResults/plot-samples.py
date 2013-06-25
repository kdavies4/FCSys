#!/usr/bin/python
"""Generate plots of FCSys examples (for FCSys.UsersGuide.SampleResults)
"""
__author__ = "Kevin Davies"
__version__ = "2013-6-29"

import os
import matplotlib.pyplot as plt

from fcres import SimRes, saveall

# Begin customize--------------------------------------------------------------

# Formats in which to save the figures (e.g., ['pdf', 'eps', 'svg', 'png'])
# If the figures shouldn't be saved, specify an empty list.
formats = ['svg', 'png']

# Input and output directories
in_dir = os.path.normpath('../'*6)
out_dir = os.path.normpath('.')

# Plot 1
name='SubregionsSound'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="Gas Velocity in Two Regions with an Initial Pressure Difference\n",
         label=os.path.join(out_dir, name),
         ynames1=['subregion1.gas.H2.faces[1, 1].phi[1]',
                  'subregion1.gas.H2.phi[1]',
                  'subregion2.gas.H2.faces[1, 1].phi[1]',
                  'subregion2.gas.H2.phi[1]'],
         legends1=['Outside walls', 'Region 1', 'Common boundary', 'Region 2'],
         ylabel1='Velocity', xunit='ms')
plt.xlim(xmax=0.3)

# Plot 2
name='ThermalConduction'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="One-Dimensional Thermal Conduction through Graphite\n",
         label=os.path.join(out_dir, name),
         ynames1=["subregion1.graphite.'C+'.T"] +
                 ["subregions[%i].graphite.'C+'.T" % i for i in range(1,9)] +
                 ["subregion2.graphite.'C+'.T"],
         legends1=['Region %i' % i for i in range(1, 11)],
         ylabel1='Temperature')

# Plot 3
name='ThermalConductionConvection'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="Gas Velocity Induced by Thermal Conduction\n",
         label=os.path.join(out_dir, name),
         ynames1=["subregion1.gas.N2.phi[1]"] +
                 ["subregions[%i].gas.N2.phi[1]" % i for i in range(1,9)] +
                 ["subregion2.gas.N2.phi[1]"],
         legends1=['Region %i' % i for i in range(1, 11)],
         leg1_kwargs=dict(loc='upper right', ncol=2),
         ylabel1='Velocity')
plt.xlim(xmax=15)

# Plot 4
name='SaturationPressure'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="Water Saturation Pressure\n",
         label=os.path.join(out_dir, name),
         xname="subregion.gas.H2O.T",
         ynames1=["subregion.gas.H2O.p", "p_sat"],
         legends1=["FCSys (from Gibbs equilibrium)",
                   "Modelica.Media (correlated function)"],
         ylabel1='Saturation pressure')

# Plot 5
name='SubregionEvaporation'
sim = SimRes(os.path.join(in_dir, name))
ax1, ax2 = sim.plot(title="Dynamic Evaporation and Condensation of Water at 25$^\circ$C\n",
                    label=os.path.join(out_dir, name),
                    ynames1=["subregion.gas.H2O.p", "p_sat"],
                    legends1=["Pressure",
                              "Saturation pressure\n(from Modelica.Media)"],
                    ylabel1='Pressure',
                    ynames2=["subregion.gas.H2O.physical.Ndot"],
                    legends2=["Rate of evaporation"],
                    ylabel2='Rate of evaporation', yunit2='umol/s')
ax1.set_ylim(ymin=0.5, ymax=6.0)
ax1_sat_frac=(3.1653 - 0.5)/5.5
ax2.set_ylim(ymin=-25*ax1_sat_frac, ymax=25*(1 - ax1_sat_frac))

# End customize----------------------------------------------------------------

# Save the plots.
saveall(formats)
plt.show()
