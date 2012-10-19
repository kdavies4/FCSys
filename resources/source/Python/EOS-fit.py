#!/usr/bin/env python
"""Fit the first virial coefficients of common gases as a function of pressure
rather than temperature.
"""

from pylab import *
from scipy import *

# If you experience problem "optimize not found", try to uncomment the following line. The problem is present at least at Ubuntu Lucid python scipy package
from scipy import optimize

# Global constants
kF = 96485.3399 # Faraday constant [C/mol] (NIST2009)
R = 8.314472/kF # gas constant [V/K] (NIST2009)
# Note that here the unit for the gas constant is V/K.  The unit of temperature
# is V (= K*R).  In other words, temperature is thermal voltage.


# H2O
# The data is from Dymond2002, p. 49-51.  Where there are multiple measurements
# at a temperature, the measurement with the smallest uncertainty is entered.
T_H2O = np.array([348.15, 360.65, 373.15, 380.56, 385.65, 392.34, 398.15, 404.47, 410.65, 423.15, 425.15, 435.65, 447.24, 448.15,])*R # Temperature [V]
B_H2O = np.array([-590.3, -515.0, -452.7, -374.0, -400.9, -350.3, -360.0, -357.6, -321.2, -320.9, -289.9, -267.0, -263.1, -227.9, -240.0])*1e6/kF # 1st virial coefficient [m3/C]
P = [(1.0 + B*Nhat)*Nhat*T for T, B in zip(T_H2O, B_H2O)]

# Make 2-D map of P to T & Nhat.
# Fit T = f(P, Nhat) to it.




exit()

P/Nhat = T + lim(P/(Nhat*T) - 1)*P/Nhat as Nhat->0
T = -lim(P/(Nhat*T))*P/Nhat as Nhat->0


T = -del(P/T)/del(Nhat)*P/Nhat
T = (P*del(T)/T - del(P))/T/del(Nhat)*P/Nhat
T^2 = (P*del(T)/T - del(P))/del(Nhat)*P/Nhat
T^2 = (P*del(ln(T)) - del(P))/del(Nhat)*P/Nhat
Let tau = ln(T)
exp(2*tau) = (P*del(tau) - del(P))/del(Nhat)*P/Nhat

del(ln(T)) = del(T)/T

dx/dy = lim((x1 - x0)/(y1 - y0)) as (y1 - y0)->0

dx/dy = lim((x1 - x0)/Nhat) as Nhat->0


# Generate data points with noise
num_points = 150
Tx = linspace(5., 8., num_points)
Ty = Tx

tX = 11.86*cos(2*pi/0.81*Tx-1.32) + 0.64*Tx+4*((0.5-rand(num_points))*exp(2*rand(num_points)**2))
tY = -32.14*cos(2*pi/0.8*Ty-1.94) + 0.15*Ty+7*((0.5-rand(num_points))*exp(2*rand(num_points)**2))




# Fit the first set
fitfunc = lambda p, x: p[0]*cos(2*pi/p[1]*x+p[2]) + p[3]*x # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
p0 = [-15., 0.8, 0., -1.] # Initial guess for the parameters
p1, success = optimize.leastsq(errfunc, p0[:], args=(Tx, tX))

time = linspace(Tx.min(), Tx.max(), 100)
plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-") # Plot of the data and the fit

# Fit the second set
p0 = [-15., 0.8, 0., -1.]
p2,success = optimize.leastsq(errfunc, p0[:], args=(Ty, tY))

time = linspace(Ty.min(), Ty.max(), 100)
plot(Ty, tY, "b^", time, fitfunc(p2, time), "b-")

# Legend the plot
title("Oscillations in the compressed trap")
xlabel("time [ms]")
ylabel("displacement [um]")
legend(('x position', 'x fit', 'y position', 'y fit'))

ax = axes()

text(0.8, 0.07,
     'x freq :  %.3f kHz \n y freq :  %.3f kHz' % (1/p1[1],1/p2[1]),
     fontsize=16,
     horizontalalignment='center',
     verticalalignment='center',
     transform=ax.transAxes)

show()
