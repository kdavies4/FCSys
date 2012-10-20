#!/usr/bin/env python
"""Load and analyze results from FCSys_.

The "fcres.py" file can be executed at the command line.  It will accept as an
argument the name of a results file or a directory with multiple files.  If no
arguments are provided, then it gives a dialog to choose the file or folder.
Finally, it provides a working session of `Python <http://www.python.org/>`_
with those results preloaded.

**Example:**

   .. code-block:: sh

      $ ./fcres.py examples/Polarization.mat
      Cell simulation results have been loaded from "examples/Polarization.mat".
      The CellSimRes instance is sim.
      In [1]:

.. _FCSys: http://www.modelica.org/libraries
.. _Modelica: http://www.modelica.org/
"""
# This file is part of FCSys (a Modelica package).
#
# Licensed by Kevin Davies under the Modelica License 2
# Copyright 2012--2012, Kevin Davies.
#
# This Modelica package is free software and the use is completely at your own
# risk; it can be redistributed and/or modified under the terms of the Modelica
# License 2.  For license conditions (including the disclaimer of warranty) see
# Modelica.UsersGuide.ModelicaLicense2 or visit
# http://www.modelica.org/licenses/ModelicaLicense2.
__author__ = "Kevin Davies"
__credits__ = "Kevin Bandy"
__version__ = "0.1"
__email__ = "kld@alumni.carnegiemellon.edu"
__status__ = "Development"

import os
import numpy as np
import matplotlib.pyplot as plt
import modelicares

from matplotlib import rcParams
from matplotlib.cbook import iterable
from collections import namedtuple
from texunit import unit2tex, label_number, label_quantity
from modelicares.helpers import (flatten, figure_w_label, setup_subplots,
                                 get_pow10, plot, convert)

# Create a container class for conditions to display in the subtitles.
Conditions = namedtuple('Conditions',
    ['cell_temp', 'composition', 'pressure', 'humidity', 'inlet_temp', 'flow'])

# Default units
# These units are used to display a variable's quantity if its displayUnit
# attribute is empty ('').  Each key is a dimension string and each entry is a
# unit string.  Both are formatted in Modelica unit notation.
# Generated from FCSys/resources/quantities.xls, 10/8/2012
default_units = {'l/(T.s)': 'cm/s2',
                 'l/T2': 'cm/s2',
                 'A': 'rad',
                 'A2': 'sr',
                 'l2': 'cm2',
                 'N2.T2/(l2.m)': 'uF',
                 'N2.T/(l2.m)': 'S',
                 'N/s': 'A',
                 'N/T': 'A',
                 'N/(l2.T)': 'A/cm2',
                 'N/(T.s)': 'A/s',
                 'N/T2': 'A/s',
                 'l2.m/T2': 'J',
                 'l2/T2': 'Sv',
                 'l.m/(T.s)': 'N',
                 'l.m/T2': 'N',
                 'l2.m/N2': 'uH',
                 'l': 'cm',
                 'l/N': 'm/mol',
                 'l2.m/(A.N.T)': 'Wb',
                 'm/(A.N.T)': 'T',
                 'A.N.T/(l2.m)': '1/Wb',
                 'm': 'g',
                 'm/N': 'g/mol',
                 'l2.m/(A.T)': 'J.s/rad',
                 's-1': '1/s',
                 '1/T': '1/s',
                 'N': 'C',
                 '1/N': '1/mol',
                 'N/l3': 'mol/cm3',
                 'l.m/N2': 'H/m',
                 'N2.T2/(l3.m)': 'F/m',
                 'l3.m/(N2.T2)': 'm/H',
                 'l2.m/(N.T2)': 'V',
                 'l3.m/(A.N.T2)': 'V.m/rad',
                 'l2.m/(N.T2.s)': 'V/s',
                 'l2.m/(N.T3)': 'V/s',
                 'l2.m/T3': 'W',
                 'l4.m/T3': 'W.m2',
                 'm/T3': 'W/m2',
                 'm.T5/l8': 'W/(m2.K4)',
                 'l2.m/(A2.T3)': 'cd',
                 'm/(l.T2)': 'kPa',
                 'm/(l.T2.s)': 'Pa/s',
                 'm/(l.T3)': 'Pa/s',
                 'T/N': '1/A',
                 'l2.m/(N2.T)': 'ohm',
                 'l.T/N': 'cm/A',
                 'T': 's',
                 'l/T': 'cm/s',
                 'l.m/(N.T)': 'kg.m/(C.s)',
                 'l3': 'cm3',
                 'l3/T': 'L/min',
                 'l3/N': 'cm3/C',
                 'l3/(N.s)': 'cm3/(mol.s)',
                 'l3/(N.T)': 'cm3/(mol.s)',
                 'A/l': 'rad/m'}

class SimRes(modelicares.SimRes):
    """Base class for Modelica_-based simulation results and methods to analyze
    those results

    .. :Note: This version of :class:`SimRes` is different than that in
       modelicares
       `http://kdavies4.github.com/modelicares/
       <http://kdavies4.github.com/modelicares/>`_.
       The Modelica_ *unit* attribute is used to indicate the dimension of a
       variable according to the method by Davies and Paredis 2012 (`ref
       <https://modelica.org/events/modelica2012/Proceedings>`_).
    """

    def get_unit(self, names):
        """Return the unit(s) for trajectory variable(s).

        The unit is the variable's *displayUnit* attribute unless it is empty
        ('').  In that case, the unit is taken from the global *default_unit*
        dictionary based on the variable's dimension.

        **Arguments:**

        - *names*: Name(s) of the variable(s) from which to get the unit(s)

             This may be a single string or (possibly nested) list of strings
             representing the names of the variables.

        If *names* is a string, then the output will be a single display unit.
        If *names* is a (optionally nested) list of strings, then the output
        will be a (nested) list of units.

        .. Note:: This method (:meth:`fcres.CellSimRes.get_unit`) is different
           than :meth:`modelicares.SimRes.get_unit`.  In FCSys_, the *unit*
           attribute is used to indicate the **dimension** of a quantity.  The
           **unit** is given by the variable's *displayUnit* attribute or the
           default unit based on the variable's dimension, as described above.

        **Example:**

           >>> from fcres import CellSimRes

           >>> sim = CellSimRes('examples/Polarization.mat')
           >>> sim.get_unit('defaults.p')
           'kPa'
           >>> sim.get_unit([[['defaults.p', 'defaults.T']]])
           [[['kPa', 'K']]]
        """
        def _get_unit(name):
            """Return the unit for a single variable.
            """
            global default_units

            if name == 'Time':
                # Time is special since it is predefined in Modelica.
                return 's'

            unit = self._traj[name].displayUnit
            if unit == '':
                try:
                    unit = default_units[self._traj[name].unit]
                except KeyError:
                    # Key error must be due to default_units; let unit = ''.
                    pass
            return unit

        return self._get(names, _get_unit)

    def _unitval(self, unitstr):
        """Convert a Modelica_ unit string to a numeric value using the base
        constants and units.

        For example, "m/s2" becomes m/s^2, where m and s are the values of keys
        'm' and 's' in the internal *U* dictionary.

        **Arguments:**

        - *unitstr*: Unit string in Modelica_ notation

             .. Seealso:: Modelica Specification, version 3.3, p. 235--236
                (https://www.modelica.org/documents)

                In summary, '.' indicates multiplication.  The denominator is
                enclosed in parentheses and begins with '/'.  Exponents
                directly follow the significand (e.g., no carat ('^')).
        """
        def _simple_unitval(unitstr):
            """Convert the numerator or denominator of a unit string to a
            number.
            """
            if unitstr == '1':
                return 1
            unitval = 1
            pointer = 0
            while pointer < len(unitstr):
                # Read the factor.
                start = pointer
                while pointer < len(unitstr) and unitstr[pointer].isalpha():
                    pointer += 1
                if start == pointer:
                    print('Unexpected character "%s" in unit substring "%s".'
                        '  It will be skipped.' % (unitstr[pointer], unitstr))
                    pointer += 1
                    continue
                try:
                    factor = self.U[unitstr[start:pointer]]
                except KeyError:
                    print('Factor "%s" in unit substring "%s" is not recognized.  '
                        'It will be skipped.'
                        % (unitstr[start:pointer], unitstr))
                    factor = 1

                # Handle the exponent, if any.
                start = pointer
                if (pointer < len(unitstr) and (unitstr[pointer] == '+' or
                    unitstr[pointer] == '-')):
                    pointer += 1
                    assert (pointer < len(unitstr) and
                            unitstr[pointer].isdigit()), ('"+" or "-" must be '
                           'followed by a digit in unit substring "%s".' %
                        unitstr)
                while pointer < len(unitstr) and unitstr[pointer].isdigit():
                    pointer += 1
                if pointer > start:
                    factor **= int(unitstr[start:pointer])

                # Combine the factor with the rest of the unit.
                unitval *= factor

                # Advance past the next dot ('.').
                if pointer < len(unitstr):
                    if unitstr[pointer] == '.':
                        pointer += 1
                    else:
                        print('Factors should be separated by dots (".").  '
                              'Unit substring "%s" is malformed.' % unitstr)

            return unitval

        # Split the numerator and the denominator (if present).
        parts = unitstr.split('/')
        if len(parts) > 2:
            print('Unit string "%s" has more than one division sign.' %
                  unitstr)

        # Remove parentheses.
        for i, part in enumerate(parts):
            if part.startswith('('):
                if part.endswith(')'):
                    parts[i] = part[1:-1]
                else:
                    print('Group "%s" in unit string "%s" begins with "(" but '
                          'does not end with ")".' % (part, unitstr))
                    parts[i] = part[1::]
        unitval = _simple_unitval(parts[0])
        for part in parts[1::]:
            unitval /= _simple_unitval(part)
        return unitval

    def unitmap(self, unitstr):
        """Return a function to map a quantity to a number in terms of a unit.

        **Arguments:**

        - *unitstr*: Unit string in Modelica_ notation

             .. Seealso:: Modelica Specification, version 3.3, p. 235--236
                (https://www.modelica.org/documents)

                In summary, '.' indicates multiplication.  The denominator is
                enclosed in parentheses and begins with '/'.  Exponents
                directly follow the significand (e.g., no carat ('^')).

        **Example:**

           >>> from fcres import CellSimRes
           >>> sim = CellSimRes('examples/Polarization.mat')
           >>> kPa = sim.unitmap('kPa')

           >>> # Method 1:
           >>> p = sim.get_IV('defaults.p')
           >>> print("The pressure is %.1f kPa." % kPa(p))
           The pressure is 149.6 kPa.

           >>> # Method 2:
           >>> p_kPa = sim.get_IV('defaults.p', kPa)
           >>> print("The pressure is %.1f kPa." % p_kPa)
           The pressure is 149.6 kPa.

           >>> # In other units:
           >>> print("The pressure is %.1f kPag." % sim.unitmap('kPag')(p))
           The pressure is 48.3 kPag.
           >>> print("The pressure is %.0f Pa." % sim.unitmap('Pa')(p))
           The pressure is 149600 Pa.
           >>> print("The pressure is %.3f atm." % sim.unitmap('atm')(p))
           The pressure is 1.476 atm.
           >>> print("The pressure is %.3f bar." % sim.unitmap('bar')(p))
           The pressure is 1.496 bar.
        """
        try:
            # Assume the unit involves an offset or a more complex function.
            return self.unitmaps[unitstr]
        except KeyError:
            # The unit is just a factor.
            return lambda x: x/self._unitval(unitstr)

    def plotfig(self, title="", label="xy",
                xname='Time', xlabel=None, xunit=None,
                ax1=None, ynames1=[], ylabel1=None, legends1=[], yunit1=None,
                leg1_kwargs=dict(loc='best'),
                ax2=None, ynames2=[], ylabel2=None, legends2=[], yunit2=None,
                leg2_kwargs=dict(loc='best'),
                prefix=False, **kwargs):
        """Create a figure (or subfigure) of 1D data as points and/or curves in
        2D Cartesian coordinates.

        **Arguments:**

        - *title*: Title of the figure

        - *label*: Label for the figure (ignored if ax is provided)

             This will be used as a base filename if the figure is saved.

        - *xname*: Name of the x-axis data

        - *xlabel*: Label for the x-axis

             If not provided, the variable's Modelica_ description string will
             be used.

        - *xunit*: String indicating the unit for the x axis

             If *xunit* is *None*, the Modelica_ variable's *displayUnit* or
             the default unit based on the variable's dimension will be used
             (in decreasing priority).

             .. Note:: Dimension checking is not currently performed, so it is
                important to ensure that a proper unit is chosen.

        - *ax1*: Primary y axes

             If *ax1* is not provided, then axes will be created in a new
             figure.

        - *ynames1*: Names of variables for the primary y axis

             If any names are invalid, then they will be skipped.

        - *ylabel1*: Label for the primary y axis

             If *ylabel1* is *None* (default) and all of the variables have
             the same Modelica_ description string, then the common description
             be used.

        - *yunit1*: String indicating the unit for the primary y-axis (see note
          for *xunit*)

             If *yunit1* is *None*, the Modelica_ *displayUnit* of the first
             entry of *ynames1* or the default unit based on that variable's
             dimension will be used (in decreasing priority).

        - *legends1*: List of legend entries for variables assigned to the
          primary y axis

             If *legends1* is an empty list ([]), ynames1 will be used.  If
             *legends1* is *None*, then no legend will be shown.

        - *leg1_kwargs*: Dictionary of keyword arguments for the primary legend

        - *ax2*, *ynames2*, *ylabel2*, *legends2*, *leg2_kwargs*: Similar to
          *ax1*, *ynames1*, *ylabel1*, *legends1*, and *leg1_kwargs*, but for
          the secondary y axis

        - *prefix*: If *True*, prefix the legend strings with the base filename
          of the class.

        - *\*\*kwargs*: Propagated to :meth:`modelicares.helpers.plot` (and
          thus to :meth:`matplotlib.pyplot.plot`)

             If both y axes are used (primary and secondary), then the *dashes*
             argument is ignored.  The curves on the primary axis will be solid
             and the curves on the secondary axis will be dotted.

        **Returns:**

        1. *ax1*: Primary y axes

        2. *ax2*: Secondary y axes

        **Example:**

           >>> from modelicares import helpers
           >>> from fcres import CellSimRes

           >>> sim = CellSimRes('examples/Polarization.mat')
           >>> sim.plotfig(xname='cell.I',
           ...             ynames1='cell.v', ylabel1="Potential",
           ...             legends1="Average voltage",
           ...             ynames2='cell.Wdot', ylabel2="Power",
           ...             legends2="Power output",
           ...             title="Cell Polarization",
           ...             label='examples/Polarization') # doctest: +ELLIPSIS
           (<matplotlib.axes.AxesSubplot object at 0x...>, <matplotlib.axes.AxesSubplot object at 0x...>)
           >>> helpers.save_figs()
           Saved examples/Polarization.pdf
           Saved examples/Polarization.png

        .. only:: html

           .. image:: examples/Polarization.png
              :scale: 50 %
              :alt: plot of cell profile

        .. only:: latex

           .. figure:: examples/Polarization.pdf
              :scale: 70 %

              Plot of cell profile
        """
        def _ystrings(ynames, ylabel, yunit, legends):
            """Generate a y-axis label and set of legend entries.
            """
            if ynames:
                if ylabel is None: # Try to create a reasonable axis label.
                    descriptions = self.get_description(ynames)
                    # If the descriptions are the same, label the y axis with
                    # the 1st one.
                    if len(set(descriptions)) == 1:
                        ylabel = descriptions[0]
                if legends == []:
                    legends = ynames
                if prefix:
                    legends = [self.fbase + ': ' + legend
                               for legend in legends]
                ylabel = label_number(ylabel, yunit)

            return ylabel, legends

        # Check the inputs (cast string to list, remove missing variables,
        # flatten)
        ynames1 = flatten([self.clean_names(ynames1)])
        ynames2 = flatten([self.clean_names(ynames2)])
        if ynames1 == [None]:
            ynames1 = None
        if ynames2 == [None]:
            ynames2 = None
        assert ynames1 is not None or ynames2 is not None, \
            "No signals were provided or they were invalid."
        if isinstance(legends1, basestring):
            legends1 = [legends1]
        if isinstance(legends2, basestring):
            legends2 = [legends2]

        # Create primary and secondary axes if necessary.
        if not ax1:
            fig = figure_w_label(label)
            ax1 = fig.add_subplot(111)
        if ynames2 and not ax2:
            ax2 = ax1.twinx()

        # Determine the units.
        if xunit is None:
            xunit = self.get_unit(xname)
        if yunit1 is None and ynames1:
            # Use unit of 1st variable, assume dimensionally consistent
            yunit1 = self.get_unit(ynames1[0])
        if yunit2 is None and ynames2:
            # Use unit of 1st variable, assume dimensionally consistent
            yunit2 = self.get_unit(ynames2[0])

        # Read the data.
        x = self.get_values(xname, f=self.unitmap(xunit))
        if xname == 'Time':
            y1 = self.get_values_at_times(ynames1, x, f=self.unitmap(yunit1))
            y2 = self.get_values_at_times(ynames2, x, f=self.unitmap(yunit2))
        else:
            times = self.get_times(xname)
            y1 = self.get_values_at_times(ynames1, times,
                                          f=self.unitmap(yunit1))
            y2 = self.get_values_at_times(ynames2, times,
                                          f=self.unitmap(yunit2))

        # Generate the x-axis label.
        if xlabel is None:
            xlabel = 'Time' if xname == 'Time' else self.get_description(xname)
        xlabel = label_number(xlabel, self.get_unit(xname))

        # Generate the y-axis labels and sets of legend entries.
        ylabel1, legends1 = _ystrings(ynames1, ylabel1, yunit1, legends1)
        ylabel2, legends2 = _ystrings(ynames2, ylabel2, yunit2, legends2)

        # Plot the data.
        if ynames1:
            if ynames2:
                # Use solid lines for primary axis and dotted lines for
                # secondary.
                try:
                    del kwargs['dashes']
                except:
                    pass
                p1 = plot(ax1, y1, x, label=legends1, dashes=[(1,0)], **kwargs)
                p2 = plot(ax2, y2, x, label=legends2, dashes=[(3,3)], **kwargs)
            else:
                p1 = plot(ax1, y1, x, label=legends1, **kwargs)
        if ynames2 and not ynames1:
            p2 = plot(ax2, y2, x, label=legends2, **kwargs)

        # Decorate the figure.
        ax1.set_title(title)
        ax1.set_xlabel(xlabel)
        #shift_scale_x(ax1)
        if ynames1:
            ax1.set_ylabel(ylabel1)
            #shift_scale_y(ax1)
        if ynames2:
            ax2.set_ylabel(ylabel2)
            #shift_scale_y(ax2)
        if legends1:
            if legends2:
                # Put the primary legend in the upper left and secondardy in
                # upper right.
                try:
                    del leg1_kwargs['loc']
                except:
                    pass
                try:
                    del leg2_kwargs['loc']
                except:
                    pass
                ax1.legend(loc=2, **leg1_kwargs)
                ax2.legend(loc=1, **leg2_kwargs)
            else:
                ax1.legend(**leg1_kwargs)
        if legends2 and not legends1:
                ax2.legend(**leg2_kwargs)

        return ax1, ax2

    def get_dimension(self, names):
        """Return the dimension(s) of trajectory variable(s).

        **Arguments:**

        - *names*: Name(s) of the variable(s) from which to get the
          dimension(s)

             This may be a single string or (possibly nested) list of strings
             representing the names of the variables.

        If *names* is a string, then the output will be a single unit.  If
        *names* is a (optionally nested) list of strings, then the output will
        be a (nested) list of units.

        **Example:**

           >>> from fcres import CellSimRes

           >>> sim = CellSimRes('examples/Polarization.mat')
           >>> sim.get_dimension('defaults.p')
           'm/(l.T2)'
           >>> sim.get_dimension([[['defaults.p', 'defaults.T']]])
           [[['m/(l.T2)', 'l2.m/(N.T2)']]]
        """
        return self._get(names, lambda name: self._traj[name].unit)

    def _set_constants(self):
        """Establish the values of constants and units.

        There are no arguments or return values.  All the constants and units
        are dependent on the base constants and units which are included in
        the defaults.base record.  See the Units package in FCSys/Units.mo and
        the Defaults model in FCSys/BCs.mo.
        """
        # Base constants and units
        rad = self.get_IV('defaults.base.rad')
        R_inf = self.get_IV('defaults.base.R_inf')
        c = self.get_IV('defaults.base.c')
        k_J = self.get_IV('defaults.base.k_J')
        R_K = self.get_IV('defaults.base.R_K')
        cd = self.get_IV("defaults.base.'cd'")
        k_F = self.get_IV('defaults.base.k_F')
        R = self.get_IV('defaults.base.R')
        self.U = dict(rad=rad, R_inf=R_inf, c=c, k_J=k_J, R_K=R_K, cd=cd,
                      k_F=k_F, R=R)

        # Mathematical constants
        import math
        pi = 2*math.acos(0)
        e = math.exp(1)
        self.U.update(pi=pi, e=e)

        # Empirical units
        m = 10973731.568539*rad/R_inf
        s = 299792458*m/c
        Wb = 483597.870e9/k_J
        S = 25812.8074434/R_K
        mol = 96485.3365*Wb*S/k_F
        K = 8.3144621*(Wb*rad)**2*S/(s*mol*R)
        self.U.update(m=m, s=s, Wb=Wb, S=S, mol=mol, K=K)

        # SI base units [BIPM2006, Table 1] and intermediate units
        V = Wb*rad/s
        A = V*S
        C = A*s
        J = V*C
        Sv = (m/s)**2
        kg = J/Sv
        self.U.update(V=V, A=A, C=C, J=J, Sv=Sv, kg=kg)

        # SI prefixes [BIPM2006, Table 5]
        yotta = 1e24,
        zetta = 1e21
        exa = 1e18
        peta = 1e15
        tera = 1e12
        giga = 1e9
        mega = 1e6
        kilo = 1e3
        hecto = 1e2
        deca = 1e1
        deci = 1e-1
        centi = 1e-2
        milli = 1e-3
        micro = 1e-6
        nano = 1e-9
        pico = 1e-12
        femto = 1e-15
        atto = 1e-18
        zepto = 1e-21
        yocto = 1e-24
        self.U.update(yotta=yotta, zetta=zetta, exa=exa, peta=peta, tera=tera,
                      giga=giga, mega=mega, kilo=kilo, hecto=hecto ,deca=deca,
                      deci=deci, centi=centi, milli=milli, micro=micro,
                      nano=nano, pico=pico, femto=femto, atto=atto,
                      zepto=zepto, yocto=yocto)

        # Coherent derived units in the SI with special names and symbols
        cyc = 2*pi*rad
        Hz = cyc/s
        sr = rad**2
        N = J/m
        Pa = N/m**2
        W = J/s
        F = C/V
        ohm = 1/S
        H = V*s/A
        T = Wb/m**2
        lm = cd*sr
        lx = lm/m**2
        Bq = Hz
        Gy = Sv
        kat = mol/s
        g = kg/kilo
        self.U.update(cyc=cyc, Hz=Hz, sr=sr, N=N, Pa=Pa, W=W, F=F, ohm=ohm,
                      H=H,  T=T, lm=lm, lx=lx, Bq=Bq, Gy=Gy, kat=kat, g=g)

        # Non-SI units accepted for use with SI units [BIPM2006, Table 6]
        min = 60*s
        hr = 60*min
        day = 24*hr
        degree = 2*pi*rad/360
        L = (0.1*m)**3
        self.U.update(min=min, hr=hr, day=day, degree=degree, L=L)

        # Derived physical constants
        # Electromagnetism
        G_0 = 2/R_K
        Phi_0 = 1/k_J
        q = G_0*Phi_0
        h = 2*q*Phi_0
        alpha = pi*1e-7*c*s*G_0/(m*S)
        Z_0 = 2*R_K*alpha
        mu_0 = Z_0/c
        epsilon_0 = 1/(Z_0*c)
        k_A = mu_0/(4*pi)
        k_e = k_A*c**2
        E_h = 2*R_inf*h*c
        eV = q*V
        self.U.update(G_0=G_0, Phi_0=Phi_0, q=q, h=h, alpha=alpha, Z_0=Z_0,
                      mu_0=mu_0, epsilon_0=epsilon_0, k_A=k_A, k_e=k_e,
                      E_h=E_h, eV=eV)
        # Electrochemistry
        N_A = k_F/q
        self.U.update(N_A=N_A)
        # Thermal physics
        k_B = R/N_A
        c_1 = cyc*h*c**2
        c_2 = h*c/k_B
        c_3_lambda = c_2/4.965114231744276
        c_3_f = 2.821439372122079*k_B/h
        sigma = 2*pi*(k_B*pi)**4/(15*(h*rad)**3*c**2)
        self.U.update(k_B=k_B, c_1=c_1, c_2=c_2, c_3_lambda=c_3_lambda,
                      c_3_f=c_3_f, sigma=sigma)

        # Selected other non-SI units from [BIPM2006, Table 8]
        bar = 1e5*Pa
        Aring = 0.1*nano*m
        self.U.update(bar=bar, Aring=Aring)

        # Additional units that are useful for fuel cells
        atm=101325*Pa
        kPa = kilo*Pa
        cm = centi*m
        mm = milli*m
        percent = centi
        self.U.update(atm=atm, kPa=kPa, cm=cm, mm=mm)
        self.U.update({'%': percent})

    def __init__(self, fname="Cell.mat", cell='cell'):
        """On initialization, load and preprocess the data.

        **Arguments:**

        - *fname*: Name of the Dymosim results trajectory file (\*.mat)

        - *cell*: Name of the cell model in the trajectory file.

             This should be relative to the top level of the simulated model
             and expressed in Modelica_ dot notation.  If the cell **is** the
             simulated model, then this should be an empty string.
        """
        self._load(fname)
        self._set_constants()

        # Save the base filename and the directory.
        self.dir, self.fbase = os.path.split(fname)
        self.fbase = os.path.splitext(self.fbase)[0]

def gen_subtitle_conditions(params, details):
    """Create a description of the operating conditions (to be used as a
    subtitle).

    params: TODO
    details: TODO
    """
    # TODO: Test and get this working.
    desc = ''
    needs_sep = False
    if details.cell_temp:
        desc += ('Cell: %s' % label_quantity(convert(params['T_cell']),
                                             params['T_cell'].unit,
                                             format='%.0f'))
        needs_sep = True
    if (details.composition or details.pressure or details.inlet_temp
        or details.humidity or idetails.humidity):
        desc += ('; ' if needs_sep else '') + 'An.|Ca.: '
        needs_sep = False
    if details.composition:
        desc += ', ' if needs_sep else ''
        if (params['anSource.Y_H2Dry'].number == 1.0
            and (params['caSource.Y_O2Dry'].number == 1.0 or
                 params['caSource.Y_O2Dry'].number == 0.21)):
            if params['caSource.Y_O2Dry'].number == 1.0:
                desc += '$H_2$|$O_2$'
            else:
                desc += '$H_2$|Air'
        else:
            desc += (label_quantity(convert(params['anSource.Y_H2Dry']),
                     params['anSource.Y_H2Dry'].unit, format='%.0f') +
                     ' $H_2$|' +
                     label_quantity(convert(params['caSource.Y_O2Dry']),
                                    params['caSource.Y_O2Dry'].unit,
                                    format='%.0f') +
                     ' $O_2$')
        needs_sep = True
    if details.pressure:
        desc += ((', ' if needs_sep else '') +
                 label_quantity(convert(params['anSink.P']), format='%.1f') +
                 '|' + label_quantity(convert(params['caSink.P']),
                 params['anSink.P'].unit, format='%.1f'))
        needs_sep = True
    if details.inlet_temp:
        desc += ((', ' if needs_sep else '') +
                 label_quantity(convert(params['anSource.T']), format='%.0f') +
                 '|' + label_quantity(convert(params['caSource.T']),
                 params['caSource.T'].unit, format='%.0f'))
        needs_sep = True
    if details.humidity:
        desc += ((', ' if needs_sep else '') +
                 label_quantity(convert(params['anSource.RH']), format='%.0f')
                 + '|' + label_quantity(convert(params['caSource.RH']),
                 params['caSource.RH'].unit, format='%.0f'))
        needs_sep = True
    if details.flow:
        desc += ((', ' if needs_sep else '') +
                 label_quantity(convert(params['anSource.SR']), format='%.1f')
                 + '|' + label_quantity(convert(params['caSource.SR']),
                 params['caSource.SR'].unit, format='%.1f'))
    return desc

def presuffix(items, prefix='', suffix=''):
    """Add a prefix and/or suffix to every entry in a list.

    The list may be multi-dimensional.  This is useful to join names of
    subregions with the name of a property.

    **Arguments:**

    - *items*: The list of items

    - *prefix*: The prefix to add

    - *suffix*: The suffix to add

    **Example:**

       >>> presuffix(['fix', 'historic', 'date', 'view'], 'pre')
       ['prefix', 'prehistoric', 'predate', 'preview']
    """
    if not prefix and not suffix:
        return items
    if iterable(items):
        if isinstance(items[0], basestring):
            return [prefix + item + suffix for item in items]
        else:
            return [presuffix(item, prefix, suffix) for item in items]

def multi_load(location):
    """Load multiple FCSys_ cell simulation and/or linearization results.

    **Arguments:**

    - *location*: Input directory, filename, or list of filenames

    **Returns:**

    1. List of cell simulations (:class:`CellSimRes` instances)

    2. List of cell linearizations (:class:`CellSimRes` instances)

    Either list may be empty.
    """
    # If updates are made here, consider also making them in modelicares.

    from glob import glob

    # Interpret the arguments.
    sims = [] # Simulation results
    lins = [] # Linearization results
    if type(location) is list:
        fnames = location
    else:
        try: # Assume the location is a directory.
            # Generate a list of files in the directory.
            fnames = glob(os.path.join(location, '*.mat'))
            fnames = [os.path.join(location, fname) for fname in fnames]
        except:
            print("Unrecognized location")
            raise

    # Load the files.
    for fname in fnames:
        try:
            sims.append(CellSimRes(fname))
            print('Cell simulation results have been loaded from "%s".' %
                  fname)
        except:
            try:
                lins.append(CellLinRes(fname))
                print('Linearization results have been loaded from "%s".' %
                      fname)
            except:
                print('Could not load cell simulation or linearization data '
                      'from "%s".  It will be skipped.' % fname)
    return sims, lins

def multi_plotfig(sims):
    """Plot a property across a list of simulations.

    **Arguments:**

    - *sims*: List of simulations
    """
    # TODO: Get this working.
    fig = plt.figure()
    #plt.setp(fig, 'label', fname)
    ax = fig.add_subplot(111)
    for sim in sims:
        sim.plotfig_subregions(ax=ax, prefix=True)
        #xy(title=title, fname=fname,
        #   ynames1=presuffix(self.storage_names,
        #                      suffix='.H2OCond.H2OlProps.T'),
        #   ylabel1="Temperature", legends1=self.storage_names, ax="1")
        #show()

class CellLinRes(modelicares.LinRes):
    """Fuel cell linearization results from FCSys_ and methods to analyze those
    results
    """
    pass # TODO: Add customizations here.

class CellSimRes(SimRes):
    """Fuel cell simulation results from FCSys_ and methods to analyze those
    results
    """

    # Global constants
    LAYERS = ['anFP', 'anGDL', 'anCL', 'pEM', 'caCL', 'caGDL', 'caFP']
    LAYER_INFO = {'anFP':"Anode flow plate",
                  'anGDL':"Anode GDL",
                  'anCL':"Anode catalyst layer",
                  'pEM':"PEM",
                  'caCL':"Cathode catalyst layer",
                  'caGDL':"Cathode GDL",
                  'caFP':"Cathode flow plate"}

    def animate_quiverfig_current(self, times, fname='current_eminus_xy'):
        """Create an animation of the x-direction current in the x-y plane.

        **Arguments:**

        - *times*: Times at which the frames should be generated

        - *fname*: Filename for the movie

            ".mpg" will be appended if necessary.
        """
        from res import animate, save_figs
        # TODO 10/28/11: Use the new movie utilities in matplotlib.

        # Set up.
        #[start_time, stop_time] = self.get_times('Time', [0,-1])
        #if not times:
        #    times = self.get_values
        #movie_dir = os.path.join(self.dir, movie_dir)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        # Plot and animate.
        for i, time in enumerate(times):
            figure_w_label(os.path.join(output_dir, '_tmp%03d' % i))
            #self.currentx_at_times(time=time)
            plt.plot(time, np.sin(time), '*') # Temporary; for debug
        save_figs('png')
        animate(fname)

    def colorfig(self, suffix, times=[0], n_rows=1, title="", subtitles=[],
                 label="color", slice_axis='z', slice_index=0,
                 xlabel="", xticklabels=[], xticks=[],
                 ylabel="", yticklabels=[], yticks=[],
                 clabel="", cbar_orientation='vertical',
                 margin_left=rcParams['figure.subplot.left'],
                 margin_right=1-rcParams['figure.subplot.right'],
                 margin_bottom=rcParams['figure.subplot.bottom'],
                 margin_top=1-rcParams['figure.subplot.top'],
                 margin_cbar=0.2,
                 wspace=0.1, hspace=0.25,
                 cbar_space=0.1, cbar_width=0.05,
                 **kwargs):
        """Create a figure with 2D scalar data at given time(s) on a color axis
        in 2D Cartesian coordinates.

        **Arguments:**

        - *suffix*: Name of the variable to be plotted (relative to the names
          of the subregions)

        - *times*: List of times at which the data should be sampled

             If multiple times are given, then subfigures will be generated.

        - *n_rows*: Number of rows of (sub)plots

        - *title*: Title for the figure

        - *subtitles*: List of subtitles (i.e., titles for each subplot)

             If not provided, "t = xx s" will be used, where xx is the time of
             each entry.  "(initial)" or "(final)" are appended if appropriate.

        - *label*: Label for the figure

             This will be used as a base filename if the figure is saved.

        - *slice_axis*: Axis normal to the screen or page ('x', 'y', or 'z')

        - *slice_index*: Position along slice_axis

        - *xlabel*: Label for the x-axes (only shown for the subplots on the
          bottom row)

        - *xticklabels*: Labels for the x-axis ticks (only shown for the
          subplots on the bottom row)

        - *xticks*: Positions of the x-axis ticks

        - *ylabel*: Label for the y axis (only shown for the subplots on the
          left column)

        - *yticklabels*: Labels for the y-axis ticks (only shown for the
          subplots on the left column)

        - *yticks*: Positions of the y-axis ticks

        - *clabel*: Label for the color- or c-bar axis

        - *cbar_orientation*: Orientation of the colorbar ("vertical" or
          "horizontal")

        - *margin_left*: Left margin

        - *margin_right*: Right margin (ignored if ``cbar_orientation ==
          'vertical'``)

        - *margin_bottom*: Bottom margin (ignored if ``cbar_orientation ==
          'horizontal'``)

        - *margin_top*: Top margin

        - *margin_cbar*: Margin reserved for the colorbar (right margin if
          ``cbar_orientation == 'vertical'`` and bottom margin if
          ``cbar_orientation == 'horizontal'``)

        - *wspace*: The amount of width reserved for blank space between
          subplots

        - *hspace*: The amount of height reserved for white space between
          subplots

        - *cbar_space*: Space between the subplot rectangles and the colorbar

        - *cbar_width*: Width of the colorbar if vertical (or height if
          horizontal)

        - *\*\*kwargs*: Propagated to :meth:`modelicares.helpers.color`
        """
        # The procedure for coordinating the images with a single colorbar was
        # copied and modified from
        # http://matplotlib.sourceforge.net/examples/pylab_examples/multi_image.html,
        # accessed 11/8/10.

        # 5/26/11: TODO: This method needs to be updated. Use quiverfig() as a
        # reference.

        from res import color

        # Get the data.
        n_plots = len(times)
        if slice_axis == 'x':
            names = presuffix(self.subregions[slice_index], suffix=suffix)
        elif slice_axis == 'y':
            names = presuffix(self.subregions[:][slice_index], suffix=suffix)
        else:
            names = presuffix(self.subregions[:][:][slice_index],
                              suffix=suffix)
        c = self.get_values_at_times(names, times)

        # Cast the data into a list of matrices (one at each sample time).
        if slice_axis == 'x':
            c = [np.array([[c[i_y + self.n_y*i_z][i]
                 for i_y in range(self.n_y-1,-1,-1)]
                 for i_z in range(self.n_z)]) for i in range(n_plots)]
        elif slice_axis == 'y':
            c = [np.array([[c[i_z + self.n_z*i_x][i]
                 for i_z in range(self.n_z-1,-1,-1)]
                 for i_x in range(self.n_x)]) for i in range(n_plots)]
        else:
            c = [np.array([[c[i_y + self.n_y*i_x][i]
                 for i_z in range(self.n_y-1,-1,-1)]
                 for i_x in range(self.n_x)]) for i in range(n_plots)]
        [start_time, stop_time] = self.get_times('Time', [0,-1])

        # Generate xlabel, xticks, and xticklabels.
        if not xlabel:
            if slice_axis == 'x':
                xlabel = 'z-axis index'
            elif slice_axis == 'y':
                xlabel = 'x-axis index'
            elif slice_axis == 'z':
                xlabel = 'y-axis index'
        if not xticklabels:
            if slice_axis == 'x':
                xticklabels = [str(i+1) for i in range(self.n_z)]
            elif slice_axis == 'y':
                xticklabels = [str(i+1) for i in range(self.n_x)]
            else:
                xticklabels = [str(i+1) for i in range(self.n_x)]
        if not xticks:
            if slice_axis == 'x':
                xticks = range(self.n_z)
            elif slice_axis == 'y':
                xticks = range(self.n_x)
            else:
                xticks = range(self.n_x)

        # Generate ylabel, yticks, and yticklabels.
        if not ylabel:
            if slice_axis == 'x':
                ylabel = 'y-axis index'
            elif slice_axis == 'y':
                ylabel = 'z-axis index'
            else:
                ylabel = 'x-axis index'
        if not yticklabels:
            if slice_axis == 'x':
                yticklabels = [str(i+1) for i in range(self.n_y)]
            elif slice_axis == 'y':
                yticklabels = [str(i+1) for i in range(self.n_z)]
            else:
                yticklabels = [str(i+1) for i in range(self.n_y)]
        if not yticks:
            if slice_axis == 'x':
                yticks = range(self.n_y)
            elif slice_axis == 'y':
                yticks = range(self.n_z)
            else:
                yticks = range(self.n_y)

        # Generate clabel.
        if not clabel:
            clabels = self.get_description(names)
            # If all of the descriptions are the same, use the first one.
            if len(set(clabels)) == 1:
                clabel = clabels[0]
            else:
                clabel = "Value"
        #units = self.get_unit(names)
        units = self.get_unit(names)
        if len(set(units)) == 1:
            clabel += label_number("", units[0])
        else:
            raise UserWarning("The variables have inconsistent units.  The "
                "colorbar unit will not match the units of all of the "
                "variables.")

        # Set up the subplots.
        if not subtitles:
            #unit = unit2tex(self.get_unit('Time'))
            unit = unit2tex(self.get_unit('Time'))
            subtitles = ["t = " + label_quantity(times[i], unit)
                for i in n_plots]
            for i, time in enumerate(times):
                if time == start_time:
                    subtitle += " (initial)"
                elif time == stop_time:
                    subtitle += " (final)"
        ax, cax = setup_subplots(n_plots=n_plots, n_rows=n_rows,
            title=title, subtitles=subtitles, label=label,
            xlabel=xlabel, xticklabels=xticklabels, xticks=xticks,
            ylabel=ylabel, yticklabels=yticklabels, yticks=yticks,
            ctype=cbar_orientation, clabel=clabel,
            margin_left=margin_left, margin_right=margin_right,
            margin_bottom=margin_bottom, margin_top=margin_top,
            margin_cbar=margin_cbar, wspace=wspace, hspace=hspace,
            cbar_space=cbar_space, cbar_width=cbar_width)

        # Create the plots.
        profiles = []
        c_min = np.inf
        c_max = -np.inf
        for i, time in enumerate(times):
            profiles.append(color(ax[i], c[i], **kwargs))

            # Find the minimum and maximum of the color-axis data.
            c_min = min(c_min, np.amin(c[i]))
            c_max = max(c_max, np.amax(c[i]))

        # Set the first image as the master, with all the others observing it
        # for changes in norm.
        class ImageFollower:
            """Update image in response to changes in clim on another image.
            """
            def __init__(self, follower):
                self.follower = follower
            def __call__(self, leader):
                self.follower.set_clim(leader.get_clim())
        norm = Normalize(c_min=c_min, c_max=c_max)
        for i, profile in enumerate(profiles):
            profile.set_norm(norm)
            if i > 0:
                profiles[0].callbacksSM.connect('changed',
                                                ImageFollower(profile))

        # Add the colorbar. (It is also based on the master image.)
        cbar = fig.colorbar(images[0], cax, orientation=cbar_orientation)

        # Scale the colorbar.
        if c_min == c_max:
            yticks = [c_min]
            for i in range(len(cbar.ax.get_yticks()) - 1):
                yticks.append(None)
            cbar.set_ticks(yticks)

    def colorfig_pressure(self, title="Pressure Profile",
                          label='pressureprofile', times=None):
        """Plot the pressure profiles of the subregions.

        **Arguments:**

        - *title*: Title for the figure

        - *label*: Label for the figure

             This will be used as a base filename if the figure is saved.

        - *times*: List of times at which the data should be sampled

             If multiple times are given, then subfigures will be generated.
        """
        if times is None:
            [start_time, stop_time] = self.get_values('Time', [0,-1])
            Delta_time = (stop_time - start_time)/5
            times = np.arange(start_time, stop_time + Delta_time/2, Delta_time)
        self.gen_profiles_fig(fname=fname, title=title,
            names=presuffix(self.storage_names, suffix='.p'),
            n_y=self.N_storagey, times=times,
            n_rows=int(np.ceil(np.sqrt(len(times)))),
            xlabel='', xticklabels=['AnFP','AnCL','CaCL','CaFP'],
            xticks=range(self.N_storagex),
            n_x_layers=self.N_storagex_layers,
            ylabel="y index", yticklabels=range(self.N_storagey-1,-1,-1),
            yticks=range(self.N_storagey), clabel=None)

        # TODO: Fix these later and make the code in gen_profiles_fig() more
        # general.
        #xticks = range(self.N_storagex)
        #xticklabels = ['AnFP','AnCL','CaCL','CaFP']

        # TODO: Create similar methods for temperature, density, etc.

    def currdenfig(self, title=None, times=[0], z_index=0, leg_kwargs={},
                   **kwargs):
        """Plot current densities of the cell segments at times.

        **Arguments:**

        - *title*: Title for the figure

        - *times*: List of times at which the data should be sampled

             If multiple times are given, then subfigures will be generated.

        - *z_index*: z-index at which the current densities are taken

        - *leg_kwargs*: Dictionary of keyword arguments for
          :meth:`matplotlib.pyplot.legend`

             If *leg_kwargs* is *None*, then no legend will be shown.

        - *\*\*kwargs*: Propagated to :meth:`barfig`
        """
        # Process the arguments.
        xlabel = kwargs.pop('xlabel', 'y-axis index (inlet to outlet)')
        name_template = self.cell + ('iprimeprime_seg[%i, ' + '%i]'
            % (z_index+1))
        #description = self.get_description(name_template % 1)
        #unit = self.get_unit(name_template % 1)
        unit = self.get_unit(name_template % 1)
        #ylabel = kwargs.pop('ylabel', label_number(description, unit))
        ylabel = kwargs.pop('ylabel', label_number("Current density", unit))

        # Create the plot.
        ax = self.barfig(names=[name_template % (i_y+1) for i_y in range(1)],
                         times=times, xlabel=xlabel, ylabel=ylabel,
                         leg_kwargs=None, **kwargs)
        for a, time in zip(ax, times):
            a.axhline(self.get_values('cell.iprimeprime', time),
                      linestyle='--', color='k', label='Entire cell')

        # Decorate.
        if title is None:
            if self.n_z == 1:
                plt.title("Current Distribution of Cell Segments")
            else:
                plt.title("Current Distribution of Cell Segments\n"
                          "z-axis index %i (of %i)" % (z_index+1, self.n_z))
        if leg_kwargs is not None:
            loc = leg_kwargs.pop('loc', 'best')
            if len(ax) == 1:
                ax[0].legend(loc=loc, **leg_kwargs)
            else:
                plt.figlegend(ax[0].lines, **leg_kwargs)

    def gen_subtitles_time(self, times):
        """Generate titles for subplots at a list of times.
        """
        time_unit = unit2tex(self.get_unit('Time'))
        subtitles = ["t = %s" % (time if time == 0 else
            label_quantity(time, time_unit)) for time in times]
        start_time, stop_time = self.get_times('Time', [0, -1])
        for i, time in enumerate(times):
            if time == start_time:
                subtitles[i] += " (initial)"
            elif time == stop_time:
                subtitles[i] += " (final)"
        return subtitles

    def get_vector_comps2D(self, axis, index, vect, times):
        """Retrieve 2-dimensional vector components at specified times.

        **Arguments:**

        - *axis*: Slice axis ('x', 'y', or 'z')

        - *index*: Position along *axis*

        - *vect*: Name of the vector to be extracted (relative to the
          subregions)

             *vect* should contain '%i' where the vector indices should be
             inserted; the vector will be indexed according to the choice of
             *slice_axis* (see below).  *vect* will be appended to the names of
             the subregions (separated by '.') in order to determine the
             full names of the vectors.

        - *times*: Times at which to sample the vectors

        **Returns:**

        1. Tuple of the x-axis vector component (u) and the y-axis vector
           component (v)

              Here, the x axis is in screen coordinates (left/right)---not
              model coordinates. The y axis is in also screen coordinates
              (bottom/top).

        2. Tuple of the number of subregions along the x- and y- axes

              This will not be the same as (*self.n_x*, *self.n_y*) if a
              chemical component is not modeled within one or more of the
              subregions.

        3. Unit of the vector components

        4. 2D boolean list indicating if a subregion actually contains the
           species)
        """
        # Create a 2D slice of the 3D space of subregions.
        names = self.subregions_slice2D(axis, index)

        # Numbers of subregions along the cross axes
        if axis == 'x':
            # Cross axes
            axes = ('y', 'z')
            n_x, n_y = self.n_y, self.n_z
        elif axis =='y':
            # Cross axes (wrap around through Cartesion space)
            axes = ('z', 'x')
            n_x, n_y = self.n_z, self.n_x
        else:
            axes = ('x', 'y')
            n_x, n_y = self.n_x, self.n_y

        # Generate the names of the variables for the u and v components of
        # the vector.  It is necessary to look up the model's indices for
        # those components.
        names_u = [[names[i_x][i_y] + '.' + vect % self.get_IV(names[i_x][i_y]
            + '.i_' + axes[0]) for i_y in range(n_y)] for i_x in range(n_x)]
        names_v = [[names[i_x][i_y] + '.' + vect % self.get_IV(names[i_x][i_y]
            + '.i_' + axes[1]) for i_y in range(n_y)] for i_x in range(n_x)]
        isvalid = [[(self.get_IV(names[i_x][i_y] + '.%sOpen' % axes[0]) == 1,
            self.get_IV(names[i_x][i_y] + '.%sOpen' % axes[1]) == 1)
            for i_y in range(n_y)]  for i_x in range(n_x)]

        # Extract the time sequences of the components at each position.
        u_traj = self.get_values_at_times(names_u, times)
        v_traj = self.get_values_at_times(names_v, times)

        # If both of the vector components are missing, it's a good indication
        # that the species isn't included in the subregion.  Mark it as blank.
        blanks = [[u_traj[i_x][i_y] is None and v_traj[i_x][i_y] is None
                  for i_y in range(n_y)] for i_x in range(n_x)]

        # Retrieve the vector components.
        u = []
        v = []
        zeros = np.zeros((n_y, n_x)) # y is 1st index (row)
        for i in range(len(times)):
            u.append(zeros)
            v.append(zeros)
            for i_x in range(n_x):
                for i_y in range(n_y):
                    if u_traj[i_x][i_y] is not None:
                        u[-1][i_y, i_x] = u_traj[i_x][i_y][i]
                    if v_traj[i_x][i_y] is not None:
                        v[-1][i_y, i_x] = v_traj[i_x][i_y][i]

        # Determine the unit of the data.  Assume that the unit of all the
        # variables are the same; there don't seem to be any cases where this
        # wouldn't be true.
        #unit = set(self.get_unit(names_u[0][0]))
        unit = self.get_unit(names_u[0][0])

        return (u, v), (n_x, n_y), unit, blanks

    def label_layers(self, axs, y=-0.11, ygap=0.03, shrink=0.2):
        """Label the layers along the x axis of plot(s).

        **Arguments:**

        - *axs*: List of axes which should be labeled

        - *y*: Vertical position of the horizontal grouping bar

        - *ygap*: Vertical gap between the grouping bar and the top of the
          labels

        - *shrink*: Fraction that each grouping bar should shrink away from
          touching its neighbor

             A value of 0.5 causes the bar to disappear.
        """
        # TODO: An attempt to make y automatic...
        #bboxes = []
        #for label in axs[-1].get_xticklabels():
        #    bbox = label.get_window_extent()
        #    # The figure transform goes from relative coords->pixels and we
        #    # want the inverse of that.
        #    bboxi = bbox.inverse_transformed(fig.transFigure)
        #    bboxes.append(bboxi)
        ## This is the bbox that bounds all the bboxes, again in relative
        ## figure coords.
        #bbox = mtransforms.Bbox.union(bboxes)
        #y = bbox.y0
        #

        # Generate a list of the layer names.
        n_x = self.n_x
        layers = [self.locate_layer(i)[0] for i in range(n_x)]
        layers.append('') # Add one more entry; be sure it is different.

        # Label the layers.
        for ax in axs:
            i_start = 0
            for i in range(self.n_x):
                if layers[i+1] != layers[i]:
                    ax.text(0.5*(i_start+i)/(n_x-1), y-ygap,
                            layers[i], ha='center', va='top',
                            transform=ax.transAxes)
                    ax.annotate('', xy=((i_start-0.5+shrink)/(n_x-1), y),
                                xytext=((i+0.5-shrink)/(n_x-1), y),
                                xycoords='axes fraction',
                                arrowprops=dict(arrowstyle='-', shrinkA=0,
                                shrinkB=0))
                    i_start = i + 1

    def layer(self, i):
        """Return the (abbreviated) name of the layer at a given x index
        (within the entire cell) and the x index within that layer.

        **Arguments:**

        - *i*: Index of the x axis (zero-based) within the entire cell---anode
          to cathode
        """
        n = 0
        for layer in self.LAYERS:
            if i < n + self.n_x_layers[layer]:
                return layer, i-n
            n += self.n_x_layers[layer]
        return None

    def plotfig_subregions(self, prop, **kwargs):
        """Plot a property within all subregions (by default, vs. time).

        **Arguments:**

        - *prop*: Name of the property

             This will be appended to the names of the subregions (separated by
             '.') in order to determine the names of the variables to plot on
             the primary y axis.

        - *\*\*kwargs*: Propagated to :meth:`modelicares.plotfig` (and thus to
          :meth:`modelicares.helpers.plot` and finally to
          :meth:`matplotlib.pyplot.plot`)
        """
        # Generate reasonable defaults.
        title = kwargs.pop('title', "%s within Subregions" % prop)
        label = kwargs.pop('label', "subregion_prop_%s" % prop)
        ynames1 = kwargs.pop('ynames1', flatten(self.subregions_w_prop(prop)))
        default_legends1 = [subregion.replace('.subregions', "")
                            for subregion in flatten(self.rel_subregions)]
        legends1 = kwargs.pop('legends1', default_legends1)

        # Create the plot.
        self.plotfig(title=title, label=label, ynames1=ynames1,
                     legends1=legends1, **kwargs)

    def quiverfig_subregions(self, vect, times=[0], n_rows=1,
                             title="", subtitles=None, label="quiver",
                             slice_axis='z', slice_index=0,
                             xlabel=None, xticklabels=None,
                             ylabel=None, yticklabels=None,
                             margin_left=rcParams['figure.subplot.left'],
                             margin_right=1-rcParams['figure.subplot.right'],
                             margin_bottom=rcParams['figure.subplot.bottom'],
                             margin_top=1-rcParams['figure.subplot.top'],
                             wspace=0.1, hspace=0.25,
                             **kwargs):
        """Create a figure with 2D vector data at given time(s) as arrows in 2D
        Cartesian coordinates.

        **Arguments:**

        - *vect*: Name of the vector to be plotted

             *vect* should contain "%i" where the vector indices should be
             inserted; the vector will be indexed according to the choice of
             *slice_axis* (see below).  *vect* will be appended to the names of
             the subregions (separated by '.') in order to determine the full
             names of the vectors.

        - *times*: List of times at which the data should be sampled

             If multiple times are given, then subfigures will be generated.

        - *n_rows*: Number of rows of (sub)plots

        - *title*: Title for the figure

        - *subtitles*: List of subtitles (i.e., titles for each subplot)

             If not provided, "t = *xx* s" will be used, where *xx* is  the
             time of each entry.  "(initial)" or "(final)" is appended if
             appropriate.

        - *label*: Label for the figure

             This will be used as a base filename if the figure is saved.

        - *slice_axis*: Axis normal to the screen or page ('x', 'y', or 'z')

        - *slice_index*: Position along *slice_axis*

        - *xlabel*: Label for the x-axes (only shown for the subplots on the
          bottom row)

             If *xlabel* is *None*, then the axis is labeled with "X Axis",
             "Y Axis", or "Z Axis" (as appropriate).

        - *xticklabels*: Labels for the x-axis ticks

             X-axis ticks are only shown for the subplots on the bottom row.
             If *xticklabels* is *None*, then the tick labels are "1", "2", ...

        - *ylabel*: Label for the y axis (only shown for the subplots on the
          left column)

             If *ylabel* is *None*, then the axis is labeled with "X Axis",
             "Y Axis", or "Z Axis" (as appropriate).

        - *yticklabels*: Labels for the y-axis ticks

             Y-axis ticks are only shown for the subplots on the left column.
             If *None*, then the tick labels are "1", "2", ...

        - *margin_left*: Left margin

        - *margin_right*: Right margin

        - *margin_bottom*: Bottom margin

        - *margin_top*: Top margin

        - *wspace*: The amount of width reserved for blank space between
          subplots

        - *hspace*: The amount of height reserved for white space between
          subplots

        - *\*\*kwargs*: Propagated to :meth:`modelicares.helpers.quiver` (and
          thus to :meth:`matplotlib.pyplot.quiver`)
        """
        # The quiver() scaling is static and manual; therefore, it shouldn't be
        # necessary to coordinate the scaling among subplots like it is for
        # colorfig().

        # TODO 10/31/11:
        # 1) Move the dots and the pivots to the centers of masses of the
        #    subregions.
        # 2) Use and show a grid that is in scale with the geometry of the
        #    subregions (create a separate method for this it so that it can be
        #    used for colorfig too; x and y scales would need to be different).
        # 3) Overlay the velocities at the faces with the velocities of the
        #    subregions (by calling this method again with hold on?).

        from res import quiver

        # Change the default quiver color.
        color = kwargs.pop('color', 'b')

        # Slice the subregions and retrieve the data.
        (us, vs), (n_x, n_y), unit, blanks = self.get_vector_comps2D(
            slice_axis, slice_index, vect, times)

        # Generate xlabel, xticks, and xticklabels.
        if xlabel is None:
            xlabel = ('y-axis index' if slice_axis == 'x' else 'z-axis index'
                      if slice_axis == 'y' else 'x-axis index')
        if xticklabels is None:
            xticklabels = [str(i+1) for i in range(n_x)]
        xticks = range(n_x)
        #xticks = [0,1,2,10] # TODO temp

        # Generate ylabel, yticks, and yticklabels.
        if ylabel is None:
            ylabel = ('z-axis index' if slice_axis == 'x' else 'x-axis index'
                      if slice_axis == 'y' else 'y-axis index')
        if yticklabels is None:
            yticklabels = [str(i+1) for i in range(n_y)]
        yticks = range(n_y)

        # Set up the subplots.
        n_plots = len(times) # Number of plots
        if not subtitles:
            subtitles = self.gen_subtitles_time(times)
        axs, n_cols = setup_subplots(n_plots=n_plots, n_rows=n_rows,
            title=title, subtitles=subtitles, label=label,
            xlabel=xlabel, xticklabels=xticklabels, xticks=xticks,
            ylabel=ylabel, yticklabels=yticklabels, yticks=yticks,
            margin_left=margin_left, margin_right=margin_right,
            margin_bottom=margin_bottom, margin_top=margin_top,
            wspace=wspace, hspace=hspace)

        # Create the plots.
        # scale = 1e-6 # 1 um/s; Zero may cause trouble wth quiverkey() TODO ?
        quivers = []
        scale = 0
        for ax, u, v in zip(axs, us, vs):
            quivers.append(quiver(ax, u, v, pad=2*0.5, color=color, **kwargs))
            # TODO test moving the locations
            #quivers.append(quiver(ax, x = [0,1,2,10], y = [0], u=u, v=v,
            #    **kwargs))

            # Add dots at the pivots.
            for i_x in range(n_x):
                for i_y in range(n_y):
                    ax.plot(i_x, i_y, marker='o', color='k',
                            markerfacecolor='w' if blanks[i_x][i_y] else 'k')

            # Find the maximum magnitude of the data.
            scale = max(scale, np.amax(np.abs(u)), np.amax(np.abs(v)))
        # TODO: Hack: For some reason, the last plot is badly scaled.  Update
        # it to match the one before it.
        if n_plots > 1:
            axs[-1].axis(axs[-2].axis())

        # Add the key.
        # Round the scale to one significant digit.
        pow10 = get_pow10(scale)
        scale = np.round(scale/10**pow10)*10**pow10
        axs[n_cols-1].quiverkey(Q=quivers[n_cols-1], X=1, Y=1.15, U=scale,
                                label=label_quantity(scale, unit),
                                labelpos='W')

        # Label the layers.
        self.label_layers(axs[len(axs)-n_cols:])

    def sankeyfig_energy(self):
        """TODO
        """
        pass

    def subregions_slice2D(self, axis, index):
        """Return a 2-dimensional slice of the names of the subregions.

        **Arguments:**

        - *axis*: The slice axis ('x', 'y', or 'z')

        - *index*: The index of the slice axis
        """
        # It would be possible to create a class for subregions that has a
        # __getiter__ method to allow direct slicing (e.g.,
        # http://stackoverflow.com/questions/2936863/python-implementing-slicing-in-getitem).
        # However, it is simpler for now to use the brute-force method.
        if axis == 'x':
            return [[self.subregions[index][i_y][i_z]
                    for i_z in range(self.n_z)] for i_y in range(self.n_y)]
        elif axis == 'y':
            return [[self.subregions[i_x][index][i_z]
                    for i_x in range(self.n_x)] for i_z in range(self.n_z)]
        else:
            return [[self.subregions[i_x][i_y][index]
                    for i_y in range(self.n_y)] for i_x in range(self.n_x)]

    def subregions_w_prop(self, prop):
        """Return a list of the names of the subregions appended with the name
        of a property.

        **Arguments:**

        - *prop*: Name of the property to be appended

             The '.' separator is automatically inserted (between the names of
             the subregion and property).
        """
        return presuffix(self.subregions, suffix='.' + prop)

# Define the reactions and species.
#REACTIONS = {'AqHO':'Aqeous H Oxidation', 'AqOR':'Aqueous O Reduction',
#             'HO':'H Oxidation', 'OR':'O Reduction',
#             'H2OCond':'H_2O Condensation'}
#SPECIES = {'eminus':'e^-', 'H2':'H_2', 'H2Og':'H_2O_{(g)}',
#           'H2Ol':'H_2O_{(l)}', 'H3Oplus':'H_3O^+', 'Hplus':'H^+', 'N2':'N_2',
#           'O2':'O_2'}

## Create the temporal plots.

# Plot the species concentration(s).
#if createAll or ishandle(curFig):
#dataStr = '.face.%s.N'
#for i in range(1):
#for i in range(len(SPECIES)):
#    N = getElementData(s, n, t, storageElement, sprintf(dataStr,
#                                                        SPECIES[i].name))
#    L_x = getElementData(s, n, t, storageElement, '.face.x.L')
#    L_y = getElementData(s, n, t, storageElement, '.face.y.L')
#    L_z = getElementData(s, n, t, storageElement, '.face.z.L')
#    conc = N/(L_x*L_y*L_z)
#    plotTrends(curFig, reshape(conc,size(conc,1),[]), t, storageElement,
#               label_number("Substance Concentration", 'mol/m3'),
#               ['Concentration of %s'SPECIES[i].text],
#               [fullfile(cd,'timeConc_'), SPECIES[i].name])

# Plot the species storage/consumption rate(s).
#if createAll or ishandle(curFig):
#    dataStr = '.face.%s.Ndot'
#for i in range(len(SPECIES)):
#    data = getElementData(s, n, t, storageElement,
#                          sprintf(dataStr, SPECIES[i].name))
#    plotTrends(curFig, reshape(data,size(data,1),[]), t, storageElement,
#               label_number("Substance Rate", 'mol/s'),
#               ['Storage/consumption Rate of %s' % SPECIES[i].text],
#               [fullfile(cd,'timeStorage_'), SPECIES[i].name])

# Plot the reaction rate(s).
#if createAll or ishandle(curFig):
#for i in range(len(reaction)):
#    dataStr = '.%s.Ndot_react'
#    data = getElementData(s, n, t, storageElement, sprintf(dataStr,
#                                                           reaction[i].name))
#    plotTrends(curFig, reshape(data,size(data,1),[]), t, storageElement,
#               label_number("Substance Rate", 'mol/s'),
#               [reaction[i].text, ' Rate'], [fullfile(cd,'timeRate_'),
#               reaction[i].name])


# Plot the species transfer rate(s).
#if createAll or ishandle(curFig):
#dataStr1 = '.face_n.%s.Ndot'
#dataStr2 = '.face_p.%s.Ndot'
#for i in range(1):
#for i in range(len(SPECIES)):
#    Ndot_n = getElementData(s, n, t, transferElement_x, sprintf(dataStr1,
#                            SPECIES[i].name))
#    Ndot_p = getElementData(s, n, t, transferElement_x, sprintf(dataStr2,
#                            SPECIES[i].name))
#    Ndot = Ndot_n - Ndot_p
#    plotTrends(curFig, reshape(Ndot,size(Ndot,1),[]), t, transferElement_x,
#               label_number("Substance Rate", 'mol/s'),
#               ['Transfer Rate of %s'%SPECIES[i].text, ' in the x direction'],
#               [fullfile(cd,'timeTransferx_'), SPECIES[i].name])

## Create the spatial plots.

# Subsample the time vector.
#t = [0:5]*t(end)/5
# The dimensions of the time vector determine the number of horizontal and
# vertical subplots.
#t = [t(1:3); t(4:6)];


# Plot the species concentration(s).
#if createAll or ishandle(curFig):
#dataStr = '.face.%s.N'
#for i in range(1):
#for i in range(len(SPECIES)):
#    N = getElementData(s, n, t, storageElement, sprintf(dataStr,
#                                                        SPECIES[i].name))
#    L_x = getElementData(s, n, t, storageElement, '.face.x.L')
#    L_y = getElementData(s, n, t, storageElement, '.face.y.L')
#    L_z = getElementData(s, n, t, storageElement, '.face.z.L')
#    conc = N./(L_x.*L_y.*L_z)
#    plotProfiles(curFig, conc, t, n_xStorage, 'Substance Concentration /mol '
#        'm^{-3}',['Concentration of ', SPECIES[i].text],
#        [fullfile(cd,'spaceConc_'),SPECIES[i].name])

# Plot the species storage/consumption rate(s).
#if createAll or ishandle(curFig):
#dataStr = '.face.%s.Ndot'
#for i in range(1):
#for i in range(len(SPECIES)):
#    data = getElementData(s, n, t, storageElement,
#        sprintf(dataStr,SPECIES[i].name))
#    plotProfiles(curFig, data, t, n_xStorage, 'Substance Rate /mol s^{-1}',
#        ['Storage/consumption Rate of ', SPECIES[i].text],
#        [fullfile(cd,'spaceStorage_'),SPECIES[i].name])

# Plot the reaction rate(s).
#if createAll or ishandle(curFig):
#for i in 3:3:
#for i = 1:length(reaction)
#    dataStr = '.%s.Ndot_react'
#    data = getElementData(s, n, t, storageElement,
#        sprintf(dataStr,reaction[i].name))
#    plotProfiles(curFig, data, t, n_xStorage, 'Substance Rate /mol s^{-1}',
#        [reaction[i].text, ' Rate'],[fullfile(cd,'spaceRate_'),
#        reaction[i].name])

# Plot the species transfer rate(s).
#if createAll or ishandle(curFig):
#dataStr1 = '.face_n.%s.Ndot'
#dataStr2 = '.face_p.%s.Ndot'
#for i in range(1):
#for i in range(len(SPECIES)):
#    Ndot_n = getElementData(s, n, t, transferElement_x,
#        sprintf(dataStr1,SPECIES[i].name))
#    Ndot_p = getElementData(s, n, t, transferElement_x,
#        sprintf(dataStr2,SPECIES[i].name))
#    Ndot = Ndot_n - Ndot_p
#    plotProfiles(curFig, Ndot, t, n_xTransfer_x, 'Substance Rate /mol s^{-1}',
#        ['Transfer Rate of ', SPECIES[i].text, ' in the x direction'],
#        [fullfile(cd,'spaceTransferx_'),SPECIES[i].name])


    def __init__(self, fname="Cell.mat", cell='cell'):
        """On initialization, load and preprocess the data.

        **Arguments:**

        - *fname*: Name of the Dymosim results trajectory file (\*.mat)

        - *cell*: Name of the cell model in the trajectory file.

             This should be relative to the top level of the simulated model
             and expressed in Modelica_ dot notation.  If the cell **is** the
             simulated model, then this should be an empty string.
        """
        self._load(fname)
        self._set_constants()

        # Save the base filename and the directory.
        self.dir, self.fbase = os.path.split(fname)
        self.fbase = os.path.splitext(self.fbase)[0]

        # Create unit maps that involve offsets or more complex mappings
        self.unitmaps = {
            # From thermodynamic temperature to degC
            'degC': lambda x: x/self._unitval('K') - 273.15,
            # From thermodynamic temperature to degF
            'degF': lambda x: (9/5.0)*(x/self._unitval('K') - 273.15) + 32,
            # From absolute pressure to kPag
            'kPag': lambda x: x/self._unitval('kPa') - 101.325,
            }
        # See "Conversions with offsets" in FCSys/configuration/units.mos.

        # Read the number of subregions in the y and z directions.
        if cell != "": cell += '.'
        self.cell = cell
        self.n_y = int(self.get_IV(cell + 'n_y'))
        self.n_z = int(self.get_IV(cell + 'n_z'))

        # Read the number of subregions in the x direction in each layer.
        self.n_x_layers = {}

        # Total number of subregions along the x axis (through-the-cell)
        self.n_x = 0
        for layer in self.LAYERS:
            n_x = 1
            while cell + layer + '.L_x[%i]' % (n_x+1) in self._traj:
                n_x += 1
            self.n_x_layers.update({layer:n_x})
            self.n_x += self.n_x_layers[layer]

        # Names of the subregions as a 3D list
        # Relative to the cell model
        self.rel_subregions = []
        for layer in self.LAYERS:
            if self.n_x_layers[layer] > 0:
                self.rel_subregions.extend([[[layer + ".subregions[%i, %i, %i]"
                    % (i_x+1, i_y+1, i_z+1) for i_z in range(self.n_z)]
                    for i_y in range(self.n_y)]
                    for i_x in range(self.n_x_layers[layer])])
        # Relative to the root of the simulated model
        self.subregions = presuffix(self.rel_subregions, prefix=cell)


if __name__ == "__main__":
    """Test the contents of this file."""
    import doctest
    doctest.testmod()
