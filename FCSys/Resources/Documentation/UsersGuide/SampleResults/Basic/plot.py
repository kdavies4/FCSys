#!/usr/bin/python
"""Generate plots of FCSys examples (for FCSys.UsersGuide.SampleResults).
"""
__author__ = "Kevin Davies"
__version__ = "2013-6-29"

import os
import matplotlib.pyplot as plt

from bisect import bisect
from numpy import exp
from matplotlib import rcParams, rc_file_defaults
from matplotlib.transforms import blended_transform_factory
from fcres import SimRes, saveall

# Begin customize--------------------------------------------------------------

# Formats in which to save the figures (e.g., ['pdf', 'eps', 'svg', 'png'])
# If the figures shouldn't be saved, specify an empty list.
formats = ['svg', 'png']

# Input directory
in_dir = '../'*7

# Plot 1
name = 'AirColumn'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="Dynamic Pressure Distribution in a Column of Air",
         label=name,
         ynames1=['subregion2.gas.N2.p_faces[1, 2]',
                  'subregion2.gas.N2.p',
                  'subregions[3].gas.N2.p_faces[1, 2]',
                  'subregions[3].gas.N2.p',
                  'subregions[2].gas.N2.p_faces[1, 2]',
                  'subregions[2].gas.N2.p',
                  'subregions[1].gas.N2.p_faces[1, 2]',
                  'subregions[1].gas.N2.p',
                  'subregion1.gas.N2.p_faces[1, 2]',
                  'subregion1.gas.N2.p',
                  'subregion1.gas.N2.p_faces[1, 1]','environment.p'],
         legends1=['Upper boundary', 'Region 1', 'Interface 1', 'Region 2',
                   'Interface 2', 'Region 3', 'Interface 3', 'Region 4',
                   'Interface 3', 'Region 4', 'Lower boundary','atmosphere'],
         leg1_kwargs=dict(loc='right'),
         ylabel1='Pressure', yunit1='kPag')
plt.xlim(0, 1.4)

# Plot 2
name = 'Echo'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="""Gas Velocity in Two Regions with an Initial Pressure Difference
with Upstream Discretization""",
         label=name,
         ynames1=['subregion1.gas.H2.faces[1, 1].phi[1]',
                  'subregion1.gas.H2.phi[1]',
                  'subregion2.gas.H2.faces[1, 1].phi[1]',
                  'subregion2.gas.H2.phi[1]'],
         legends1=['Outside walls', 'Region 1', 'Common boundary', 'Region 2'],
         leg1_kwargs=dict(loc='upper right'),
         ylabel1='Velocity', xunit='ms')
plt.xlim(xmax=0.3)

# Plot 3
name = 'EchoCentral'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="""Gas Velocity in Two Regions with an Initial Pressure Difference
with Central Difference""",
         label=name,
         ynames1=['subregion1.gas.H2.faces[1, 1].phi[1]',
                  'subregion1.gas.H2.phi[1]',
                  'subregion2.gas.H2.faces[1, 1].phi[1]',
                  'subregion2.gas.H2.phi[1]'],
         legends1=['Outside walls', 'Region 1', 'Common boundary', 'Region 2'],
         leg1_kwargs=dict(loc='upper right'),
         ylabel1='Velocity', xunit='ms')
plt.xlim(xmax=0.3)

# Plot 4
name = 'ElectricalConduction'
sim = SimRes(os.path.join(in_dir, name))
rcParams['figure.subplot.bottom'] = 0.12
ax1, ax2 = sim.plot(title="Heating of an Electrical Resistor",
                    label=name,
                    ynames1=["subregion.graphite.'C+'.faces[1, 1].T",
                             "subregion.graphite.'C+'.T", "T_ex"],
                    legends1=["Boundary", "Bulk", "Expected"],
                    leg1_kwargs=dict(loc='right'),
                    ylabel1="Temperature")

# Expected time constant
tau = sim.get_IV("subregion.graphite.'C+'.tau_QT[1]", sim.to_unit('s'))
plt.axvline(tau, color='k', linestyle=':')
trans = blended_transform_factory(ax1.transData, ax1.transAxes)
plt.text(tau, -0.01, r"$\tau_{QT}$", ha='center', va='top', rotation=90,
         transform=trans)
# Actual time constant
T = sim["subregion.graphite.'C+'.T"]
T_target = T[0] + (1 - exp(-1))*(T[-1] - T[0])
i = bisect(T, T_target)
t = sim.get_times("subregion.graphite.'C+'.T")[i]
plt.axvline(t, color='k', linestyle=':')
plt.text(t, -0.01, "@63%", ha='center', va='top', rotation=90, transform=trans)

# Plot 5
name = 'Evaporation'
sim = SimRes(os.path.join(in_dir, name))
ax1, ax2 = sim.plot(title="Dynamic Evaporation and Condensation of Water at 25$^\circ\!$C",
                    label=name,
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

# Plot 6
name = 'HOR'
sim = SimRes(os.path.join(in_dir, name))
#**

# Plot 7
name = 'InternalFlow'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="Laminar Internal Flow with Additional Small Signal Velocity",
         label=name,
         ynames1=['Deltap', 'Deltap_Poiseuille'],
         legends1=['Model (with inertance)',
                   "Poiseuille's law (without inertance)"],
         leg1_kwargs=dict(loc='upper right'), ylabel1='Pressure difference')
plt.xlim(xmax=400)
plt.ylim(ymax=-2.5)


# Plot 7
name = 'ORR'
sim = SimRes(os.path.join(in_dir, name))
#**

# Plot 8
name = 'SaturationPressure'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="Water Saturation Pressure",
         label=name,
         xname="subregion.gas.H2O.T",
         ynames1=["subregion.gas.H2O.p", "p_sat"],
         legends1=["FCSys (via Gibbs equilibrium)",
                   "Modelica.Media (correlated function)"],
         ylabel1='Saturation pressure')

# Plot 9
name = 'ThermalConduction'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="One-Dimensional Thermal Conduction through Graphite",
         label=name,
         ynames1=["subregion1.graphite.'C+'.T"] +
                 ["subregions[%i].graphite.'C+'.T" % i for i in range(1,9)] +
                 ["subregion2.graphite.'C+'.T"],
         legends1=['Region %i' % i for i in range(1, 11)],
         leg1_kwargs=dict(loc='upper right'),         ylabel1='Temperature')

# Plot 10
name = 'ThermalConductionConvection'
sim = SimRes(os.path.join(in_dir, name))
sim.plot(title="Gas Velocity Induced by Thermal Conduction",
         label=name,
         ynames1=["subregion1.gas.N2.phi[1]"] +
                 ["subregions[%i].gas.N2.phi[1]" % i for i in range(1,9)] +
                 ["subregion2.gas.N2.phi[1]"],
         legends1=['Region %i' % i for i in range(1, 11)],
         leg1_kwargs=dict(loc='upper right', ncol=2),
         ylabel1='Velocity')
plt.xlim(xmax=8)

# End customize----------------------------------------------------------------

# Save the plots.
saveall(formats)
#plt.show()
